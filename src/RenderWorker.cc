#include "RenderWorker.hh"
#include "BooleanTree.hh"

#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VSolid.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VisExtent.hh>

#include <QCollator>
#include <QColor>
#include <QImage>
#include <QPointF>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>

long recursivelySumNCalls(const std::vector<Element> &elts,
                          const ElemMutables e[], int idx) {
    long s = e[elts[idx].ecode].ngeocalls;
    for (size_t k = 0; k < elts[idx].children.size(); k++) {
        s += recursivelySumNCalls(elts, e, elts[idx].children[k]);
    }
    return s;
}

void recursivelyPrintNCalls(const std::vector<Element> &elts,
                            const ElemMutables e[], int depth, long net,
                            int idx) {
    if (net == 0) {
        net = recursivelySumNCalls(elts, e, idx);
    }
    char ws[256];
    int i = 0;
    for (; i < depth && i < 255; i++) {
        ws[i] = ' ';
    }
    ws[i] = '\0';
    qDebug("%s%s %.5f %ld", ws, elts[idx].name.data(),
           e[elts[idx].ecode].ngeocalls / double(net),
           e[elts[idx].ecode].ngeocalls);
    for (size_t k = 0; k < elts[idx].children.size(); k++) {
        recursivelyPrintNCalls(elts, e, depth + 1, net, elts[idx].children[k]);
    }
}

static G4ThreeVector condrot(const Element &e, const G4ThreeVector &vec) {
    if (e.rotated) {
        return e.rot * vec;
    }
    return vec;
}

static G4ThreeVector condirot(const Element &e, const G4ThreeVector &vec) {
    if (e.rotated) {
        return e.rot.inverse() * vec;
    }
    return vec;
}

static bool clipRay(const std::vector<Plane> &clipping_planes,
                    const G4ThreeVector &init, const G4ThreeVector &forward,
                    G4double &sdist, G4double &edist,
                    G4ThreeVector &entrynormal, G4ThreeVector &exitnormal) {
    sdist = 0.0;
    edist = kInfinity;
    // Computes esdist, entryexitnormal. Returns 1 iff a ray exists
    for (uint i = 0; i < clipping_planes.size(); i++) {
        const Plane &p = clipping_planes[i];
        G4double crs = p.normal * forward;
        if (crs > 0.0) {
            G4double ndist = (p.offset - p.normal * init) / crs;
            if (ndist > sdist) {
                entrynormal = -p.normal / p.normal.mag();
            }
            sdist = std::max(sdist, ndist);
        } else if (crs < 0.0) {
            G4double fdist = (p.offset - p.normal * init) / crs;
            if (fdist < edist) {
                exitnormal = p.normal / p.normal.mag();
            }
            edist = std::min(edist, fdist);
        } else {
            if (init * p.normal < p.offset) {
                // Either entire ray is in or out of plane
                return false;
            }
        }
    }
    if (edist < sdist) {
        // Planes clip everything out
        return false;
    }
    return true;
}

G4ThreeVector forwardDirection(const G4RotationMatrix &orientation) {
    return orientation.rowX();
}
G4ThreeVector initPoint(const QPointF &scpt, const ViewData &d) {
    const G4ThreeVector forward = d.orientation.rowX();
    const G4ThreeVector updown = d.orientation.rowY();
    const G4ThreeVector lright = d.orientation.rowZ();
    const G4ThreeVector init =
        d.camera + updown * scpt.y() * d.scale + lright * scpt.x() * d.scale;
    return init;
}

RayPoint traceRay(const G4ThreeVector &init, const G4ThreeVector &forward,
                  const std::vector<Element> &els,
                  const std::vector<Plane> &clipping_planes, Intersection *ints,
                  int maxhits, long iteration, ElemMutables mutables[],
                  bool first_visible_hit) {
    RayPoint ret;
    ret.N = 0;
    ret.intersections = ints;
    ret.front_clipped = false;
    ret.back_clipped = false;

    // Pseudorandom number to prevent consistent child ordering
    // so as to make overlapping children obvious
    const size_t rotfact = size_t(random());
    // Minimum feature size, yet still above double discrimination threshold
    //    const G4double epsilon = 1e-3 * CLHEP::nanometer;

    G4double sdist, edist;
    G4ThreeVector entrynormal;
    G4ThreeVector exitnormal;
    bool exists = clipRay(clipping_planes, init, forward, sdist, edist,
                          entrynormal, exitnormal);
    if (!exists) {
        return ret;
    }

    G4ThreeVector start = init + forward * sdist;
    const Element &root = els[0];
    ++mutables[root.ecode].ngeocalls;
    // Ensure ray starts in the root.
    bool clippable = false;
    if (!root.solid->Inside(start)) {
        G4double jdist = root.solid->DistanceToIn(start, forward);
        if (jdist + sdist >= edist) {
            // Root solid not reachable with available distance
            return ret;
        }
        sdist += jdist;
        start += forward * jdist;
        clippable = false;
    } else {
        clippable = true;
    }

    const int maxdepth = 10;
    // Assume no depth greater than depth
    // TODO: make all such buffers created but once and
    // passed in, so that stack elements stay close

    // Assume we are already in the root element.
    // We statically allocate, because many new/frees are too expensive
    const Element *stack[maxdepth];
    stack[0] = &root;
    size_t n = 1;
    G4ThreeVector local = start;
    for (int iter = 0; iter < 1000; iter++) {
        const Element &last = *stack[n - 1];

        // Strictly inclusion tests....
        // check inside on all children, append; stop when no longer dropping
        // levels
        bool found = false;
        for (size_t walk = 0; walk < last.children.size(); walk++) {
            // rotate/offset start point...
            const Element &elem =
                els[last.children[(walk + rotfact) % last.children.size()]];
            ++mutables[elem.ecode].ngeocalls;
            if (elem.solid->Inside(condrot(elem, (local + elem.offset)))) {
                stack[n] = &elem;
                ++n;
                found = true;
                local += elem.offset;
                // Picking the first intersection can lead to visual glitches
                // when two or more children overlap at the same point.
                // One solution: pick a random child order..
                break;
            }
        }
        if (!found) {
            break;
        }
    }
    G4ThreeVector offsets[maxdepth];
    offsets[0] = start;
    for (size_t i = 1; i < n; i++) {
        const Element &elem = *stack[i];
        offsets[i] = offsets[i - 1] + elem.offset;
    }
    const G4double isdist = sdist;

    if (clippable) {
        // Start point is in the world volume, and is an intersection
        bool store = !first_visible_hit || stack[n - 1]->visible;
        if (store) {
            ints[0].dist = sdist;
            ints[0].normal = entrynormal;
            ints[0].ecode = stack[n - 1]->ecode;
            ret.front_clipped = true;
            ret.N++;
            if (first_visible_hit) {
                return ret;
            }
        }
    } else {
        // Intersection on entry
        bool store = !first_visible_hit || root.visible;
        if (store) {
            ++mutables[root.ecode].ngeocalls;
            G4ThreeVector lnormal = root.solid->SurfaceNormal(local);
            ints[0].dist = sdist;
            ints[0].normal = lnormal;
            ints[0].ecode = root.ecode;
            ret.front_clipped = false;
            ret.N++;
            if (first_visible_hit) {
                return ret;
            }
        }
    }

    for (int iter = 0; iter < 1000; iter++) {
        const Element &curr = *stack[n - 1];
        const G4ThreeVector &pos = offsets[n - 1] + forward * (sdist - isdist);

        ElemMutables &cmu = mutables[curr.ecode];
        if (cmu.niter != iteration || cmu.abs_dist <= sdist /* epsilon ? */) {
            ++cmu.ngeocalls;
            cmu.abs_dist =
                sdist + curr.solid->DistanceToOut(condrot(curr, pos),
                                                  condrot(curr, forward));
            cmu.niter = iteration;
        }
        G4double exitdist = cmu.abs_dist - sdist;

        G4double closestDist = 2 * kInfinity;
        const Element *closest = NULL;
        G4ThreeVector closestPos;
        for (size_t walk = 0; walk < curr.children.size(); walk++) {
            const Element &elem =
                els[curr.children[(walk + rotfact) % curr.children.size()]];
            G4ThreeVector sub = pos + elem.offset;

            ElemMutables &emu = mutables[elem.ecode];
            if (emu.niter != iteration || emu.abs_dist <= sdist) {
                ++emu.ngeocalls;
                // Inside case happens when shapes are way out of bounds.
                // (Note: ought to make this optional, since I doubt
                //  Geant's regular walking handles this.)
                if (sdist > emu.abs_dist &&
                    elem.solid->Inside(condrot(elem, sub))) {
                    emu.abs_dist = sdist;
                } else {
                    if (sdist < emu.abs_dist) {
                        ++emu.ngeocalls;
                    }
                    emu.abs_dist =
                        sdist + elem.solid->DistanceToIn(
                                    condrot(elem, sub), condrot(elem, forward));
                }
                emu.niter = iteration;
            }
            G4double altdist = emu.abs_dist - sdist;
            if (altdist >= kInfinity) {
                // No intersection
                continue;
            } else if (altdist < closestDist) {
                closestPos = sub;
                closest = &elem;
                closestDist = altdist;
            }
        }
        G4double fdist = std::min(exitdist, closestDist);
        if (fdist > edist - sdist) {
            // End clip (counts as a hit!; always relevant)
            ints[ret.N].dist = edist;
            ints[ret.N].normal = exitnormal;
            ints[ret.N].ecode = -2;
            ret.back_clipped = true;
            ret.N++;
            return ret;
        } else if (exitdist < closestDist) {
            // Transition on leaving vol w/o intersections
            sdist += exitdist;
            G4ThreeVector lpos = offsets[n - 1] + forward * (sdist - isdist);
            // Drop from stack
            --n;

            // Record hit with normal
            bool store = !first_visible_hit || (n > 0 && stack[n - 1]->visible);
            if (store) {
                ++cmu.ngeocalls;
                G4ThreeVector lnormal =
                    curr.solid->SurfaceNormal(condrot(curr, lpos));
                ints[ret.N].dist = sdist;
                ints[ret.N].normal = condirot(curr, lnormal);
                ints[ret.N].ecode = n > 0 ? stack[n - 1]->ecode : -2;
                ret.N++;
                if (n == 0 || ret.N >= maxhits || first_visible_hit) {
                    if (first_visible_hit) {
                        qFatal("Unexpected discovery on exit");
                    }
                    return ret;
                }
            } else if (n == 0) {
                return ret;
            }
        } else {
            // Transition on visiting a child
            stack[n] = closest;
            offsets[n] = closest->offset + offsets[n - 1];
            sdist += closestDist;
            n++;
            G4ThreeVector lpos = offsets[n - 1] + forward * (sdist - isdist);

            // Record hit with normal
            bool store = !first_visible_hit || closest->visible;
            if (store) {
                ++mutables[closest->ecode].ngeocalls;
                G4ThreeVector lnormal =
                    closest->solid->SurfaceNormal(condrot(*closest, lpos));
                ints[ret.N].dist = sdist;
                ints[ret.N].normal = condirot(*closest, lnormal);
                ints[ret.N].ecode = (ret.N >= maxhits) ? -2 : closest->ecode;
                ret.N++;
                if (ret.N >= maxhits || first_visible_hit) {
                    return ret;
                }
            }
        }
    }
    return ret;
}

void compressTraces(RayPoint *ray, const std::vector<Element> &elts) {
    if (ray->N <= 1) {
        return;
    }
    Intersection *ints = ray->intersections;
    const G4double epsilon = 0.1 * CLHEP::nm;
    // First, nuke the empty slices
    int n = 0;
    for (int i = 0; i < ray->N - 1; i++) {
        G4double sep = std::abs(ints[i + 1].dist - ints[i].dist);
        if (sep > epsilon) {
            // Do both: evtly it overrides
            ints[n] = ints[i];
            ints[n + 1] = ints[i + 1];
            n++;
        }
    }

    // Finally, nuke intersections between similar material'd objects
    // (If one is visible/one isn't, don't contract
    int p = 0;
    // 001122334  =>  0022334
    // |x|x|y|x|  =>  |x|y|x|
    for (int i = 0; i < n; i++) {
        if (i == 0 ||
            elts[ints[i].ecode].ccode != elts[ints[i - 1].ecode].ccode ||
            elts[ints[i].ecode].visible != elts[ints[i - 1].ecode].visible) {
            ints[p] = ints[i];
            p++;
        }
    }
    ints[p] = ints[n];

    ray->N = p;
}

static FColor colorMap(const Intersection &intersection,
                       const G4ThreeVector &forward, const VColor &base,
                       const G4ThreeVector &position, double shade_scale,
                       bool is_clipping_plane) {
    const G4ThreeVector &normal = intersection.normal;

    // Opposed normals (i.e, for transp backsides) are mirrored
    G4double cx = std::abs(std::acos(-normal * forward) / CLHEP::pi);
    cx = 1.0 - std::max(0.0, std::min(1.0, 0.7 * cx));

    if (is_clipping_plane) {
        const G4ThreeVector &orthA = normal.orthogonal().unit();
        const G4ThreeVector &orthB = normal.cross(orthA).unit();

        double aslp = orthA * position;
        double bslp = orthB * position;

        double shade_factor = (aslp + bslp) * shade_scale;
        shade_factor *= 40.;
        double fm = std::fmod(shade_factor, 1.);
        if (fm < 0.)
            fm += 1.;
        if (fm < 0.5)
            cx *= 0.75;
    }

    return FColor(cx * base.redF(), cx * base.greenF(), cx * base.blueF());
}

void countTree(const std::vector<Element> &els, int index, int &treedepth,
               int &nelements) {
    // Count nelements & depth
    int t = 0;
    int n = 0;
    for (int s : els[index].children) {
        int at, an;
        countTree(els, s, at, an);
        t = std::max(at, t);
        n += an;
    }
    treedepth = t + 1;
    nelements = n + 1;
}

RayPoint rayAtPoint(const QPointF &pt, qreal radius,
                    const G4ThreeVector &forward, const ViewData &d, int &iter,
                    ElemMutables *mutables, Intersection *ints,
                    Intersection *altints, int M, int *ndevs) {
    RayPoint ray =
        traceRay(initPoint(pt, d), forward, d.elements, d.clipping_planes, ints,
                 M, iter, mutables, d.force_opaque);

    iter++;
    compressTraces(&ray, d.elements);

    if (d.level_of_detail <= -1) {
        for (int k = 0; k < ray.N; k++) {
            ndevs[k] = 0;
        }
        G4double seed = CLHEP::pi / 5 * (qrand() >> 16) / 16384.0;
        for (int k = 0; k < 10; k++) {
            QPointF off(radius * std::cos(seed + k * CLHEP::pi / 5),
                        radius * std::sin(seed + k * CLHEP::pi / 5));
            RayPoint aray = traceRay(initPoint(pt + off, d), forward,
                                     d.elements, d.clipping_planes, altints, M,
                                     iter, mutables, d.force_opaque);
            iter++;
            compressTraces(&aray, d.elements);
            // At which intersection have we disagreements?
            int minN = std::min(aray.N, ray.N);
            if (minN > 0) {
                for (int l = 0; l < minN; l++) {
                    Intersection intr = ray.intersections[l],
                                 aintr = aray.intersections[l];

                    const G4double jump =
                        1.0 * radius * d.scale /
                        std::abs(-G4ThreeVector(aintr.normal) *
                                 d.orientation.rowX());
                    bool diffmatbehind = false;
                    if (l < minN - 1) {
                        diffmatbehind = d.split_by_material
                                            ? d.elements[aintr.ecode].ccode !=
                                                  d.elements[intr.ecode].ccode
                                            : aintr.ecode != intr.ecode;
                    }
                    bool edgediff = (std::abs(G4ThreeVector(aintr.normal) *
                                              G4ThreeVector(intr.normal)) <
                                         0.3 || // Q: abs?
                                     std::abs(aintr.dist - intr.dist) > jump);
                    if (diffmatbehind || edgediff) {
                        ++ndevs[l];
                    }
                }
            }
            // Count where ray has intersections but aray doesn't
            for (int l = aray.N; l < ray.N; l++) {
                ++ndevs[l];
            }
        }
        const int devthresh = 2;
        for (int k = 0; k < ray.N; ++k) {
            if (ndevs[k] >= devthresh) {
                ray.N = k;
                if (k > 0) {
                    ray.intersections[k - 1].ecode = (-1);
                }
                break;
            }
        }
    }

    return ray;
}

QRgb colorForRay(const RayPoint &ray, QRgb trackcol, G4double trackdist,
                 const ViewData &d, const QPointF &pt,
                 const G4ThreeVector &forward) {
    bool rayoverride = false;
    int p = ray.N - 1;
    for (int k = 0; k < ray.N; k++) {
        if (ray.intersections[k].dist > trackdist) {
            p = k - 1;
            // WARNING: there exist exceptions!
            // NOT QUITE ACCURATE WITH LAYERING! TODO FIXME
            rayoverride = true;
            break;
        }
    }
    bool line = ray.N > 0 && (ray.intersections[ray.N - 1].ecode == (-1));
    if (line && p == ray.N - 1) {
        p = ray.N - 2;
    }
    FColor col =
        (line && !rayoverride) ? FColor(0., 0., 0., 1.f) : FColor(trackcol);
    // p indicates the first volume to use the color rule on
    // (p<0 indicates the line dominates)
    for (int k = p; k >= 0; --k) {
        // We use the intersection before the volume
        const Element &eback = d.elements[ray.intersections[k].ecode];
        const VColor &base_color = d.color_table[eback.ccode];
        const G4ThreeVector &intpos =
            initPoint(pt, d) + ray.intersections[k].dist * forward;
        const FColor altcol =
            colorMap(ray.intersections[k], forward, base_color, intpos,
                     1. / d.scale, (k == 0 && ray.front_clipped));
        double e = (d.force_opaque ? 1. : eback.alpha);
        if (!eback.visible) {
            continue;
        }
        col = FColor::blend(altcol, col, 1 - e);
    }
    return col.rgba();
}

static TrackMetaData *setupBallRadii(const TrackBlock *blocks, size_t ntracks) {
    TrackMetaData *meta = new TrackMetaData[ntracks];
    size_t i = 0;

    for (size_t z = 0; z < ntracks; z++) {
        const TrackHeader &h = blocks[i].h;
        const TrackBlock *pts = &blocks[i + 1];
        i += h.npts + 1;

        G4ThreeVector p0(pts[0].p.x, pts[0].p.y, pts[0].p.z);
        G4double radius = 0.0;
        for (int32_t j = 0; j < h.npts; j++) {
            const TrackPoint &tj = pts[j].p;
            G4ThreeVector pj(tj.x, tj.y, tj.z);
            radius = std::max(radius, (pj - p0).mag());
        }
        meta[z].ballRadius = radius;
    }
    return meta;
}

TrackData::TrackData() {}

TrackData::TrackData(const char *filename) {
    data = QSharedDataPointer<TrackPrivateData>(new TrackPrivateData(filename));
}

TrackData::TrackData(const TrackData &other) : data(other.data) {}

static void lineCuts(const std::vector<Plane> &clips, G4ThreeVector from,
                     G4ThreeVector to, G4double &lowE, G4double &highE) {
    G4double low = -kInfinity;
    G4double high = kInfinity;
    for (const Plane &p : clips) {
        G4double fpos = p.offset - p.normal * from;
        G4double tpos = p.offset - p.normal * to;
        if (fpos - tpos == 0.) {
            // TODO sign may be off
            if (fpos > 0) {
                // On cut side of parallel plane
                low = kInfinity;
                high = -kInfinity;
            } else {
                // No change, inside simplex
            }
            continue;
        }
        // [0,1] becomes [fpos,tpos]: inv img of 0.
        // x' = (x - 0) * (fpos - tpos) + fpos
        // x = (x' - fpos) /  (fpos - tpos)
        G4double pos = fpos / (fpos - tpos);
        if (tpos < fpos) {
            low = std::max(pos, low);
        } else {
            high = std::min(pos, high);
        }
    }
    lowE = low;
    highE = high;
}

static bool insideConvex(const std::vector<Plane> &clips, G4ThreeVector point) {
    // Optimizer will resolve to simple form anyway; & this is error free
    double lowE, highE;
    lineCuts(clips, point, point, lowE, highE);
    return lowE < highE;
}

static TrackPoint linmix(const TrackPoint &a, const TrackPoint &b, double mx) {
    TrackPoint c;
    double dx = 1.0 - mx;
    c.x = dx * a.x + mx * b.x;
    c.y = dx * a.y + mx * b.y;
    c.z = dx * a.z + mx * b.z;
    c.time = float(dx) * a.time + float(mx) * b.time;
    c.energy = float(dx) * a.energy + float(mx) * b.energy;
    return c;
}

TrackData::TrackData(const TrackData &other, const ViewData &vd,
                     const Range &trange, const Range &erange,
                     const IRange &selidxs, const IRange &genrange,
                     const QMap<int, bool> &type_active) {
    size_t otracks = other.getNTracks();
    const TrackBlock *oblocks = other.getBlocks();
    const TrackMetaData *ometa = other.getMeta();

    // Worst case allocation is one track (header+2 pts) per line segment
    TrackBlock *buf =
        (TrackBlock *)malloc(3 * other.getNBlocks() * sizeof(TrackBlock));
    size_t qtracks = 0;
    size_t qblocks = 0;

    float tlow = float(trange.low), thigh = float(trange.high);
    float elow = float(erange.low), ehigh = float(erange.high);
    size_t nlow = std::max(size_t(0), selidxs.low - 1);
    size_t nhigh = std::min(otracks, selidxs.high);
    // Scan up to track (selidxs.low-1)

    size_t i = 0;
    for (size_t z = 0; z < nlow; z++) {
        i += oblocks[i].h.npts + 1;
    }

    for (size_t z = nlow; z < nhigh; z++) {
        const TrackHeader &oheader = oblocks[i].h;
        const TrackPoint *seq = (TrackPoint *)&oblocks[i + 1];
        i += oheader.npts + 1;
        size_t qheader = -1;

        bool typekeep = type_active.value(oheader.ptype, type_active[0]);
        if (!typekeep || (genrange.low > ometa[z].generation ||
                          genrange.high < ometa[z].generation)) {
            continue;
        }

        G4ThreeVector fts(seq[0].x, seq[0].y, seq[0].z);
        bool intime = seq[0].time >= tlow && seq[0].time <= thigh;
        bool inenergy = seq[0].energy >= elow && seq[0].energy <= ehigh;
        bool inconvex = insideConvex(vd.clipping_planes, fts);
        bool started = false;
        if (inconvex && intime && inenergy) {
            TrackHeader h = oheader;
            h.npts = 1;
            qtracks++;
            qheader = qblocks;
            buf[qblocks].h = h;
            buf[qblocks + 1].p = seq[0];
            qblocks += 2;
            started = true;
        }

        for (int j = 0; j < oheader.npts - 1; j++) {
            double low, high;
            G4ThreeVector pl(seq[j].x, seq[j].y, seq[j].z);
            G4ThreeVector ph(seq[j + 1].x, seq[j + 1].y, seq[j + 1].z);
            lineCuts(vd.clipping_planes, pl, ph, low, high);
            float itrange = (seq[j + 1].time - seq[j].time == 0.f)
                                ? float(kInfinity)
                                : 1 / (seq[j + 1].time - seq[j].time);
            float tcutlow = (tlow - seq[j].time) * itrange;
            float tcuthigh = (thigh - seq[j].time) * itrange;
            float ierange = (seq[j + 1].energy - seq[j].energy == 0.f)
                                ? float(kInfinity)
                                : 1 / (seq[j + 1].energy - seq[j].energy);
            float ecutlow = (elow - seq[j].energy) * ierange;
            float ecuthigh = (ehigh - seq[j].energy) * ierange;
            low = std::max(double(std::max(tcutlow, ecutlow)), low);
            high = std::min(double(std::min(tcuthigh, ecuthigh)), high);

            // TODO: modify casework to localize results, not inputs
            if (low <= 0. && high >= 1.) {
                // Keep point
                if (!started) {
                    TrackHeader h = oheader;
                    h.npts = 1;
                    qheader = qblocks;
                    buf[qblocks].h = h;
                    buf[qblocks + 1].p = seq[j];
                    qblocks += 2;
                }
                buf[qheader].h.npts++;
                buf[qblocks].p = seq[j + 1];
                qblocks++;
            } else if (low <= 0. && high < 1.) {
                if (high > 0.) {
                    if (!started) {
                        TrackHeader h = oheader;
                        h.npts = 1;
                        qtracks++;
                        qheader = qblocks;
                        buf[qblocks].h = h;
                        buf[qblocks + 1].p = seq[j];
                        qblocks += 2;
                    }
                    buf[qheader].h.npts++;
                    buf[qblocks].p = linmix(seq[j], seq[j + 1], high);
                    qblocks++;
                } else {
                    // Not there entirely
                }
            } else if (low > 0. && high >= 1.) {
                if (low < 1.) {
                    // Starting new sequence
                    TrackHeader h = oheader;
                    h.npts = 2;
                    qtracks++;
                    qheader = qblocks;
                    buf[qblocks].h = h;
                    buf[qblocks + 1].p = linmix(seq[j], seq[j + 1], low);
                    buf[qblocks + 2].p = seq[j + 1];
                    qblocks += 3;
                } else {
                    // Not there entirely
                }
            } else {
                if (high < 1. && low > 0. && low < high) {
                    // low > 0, high < 1
                    // A whole new short segment...
                    TrackHeader h = oheader;
                    h.npts = 2;
                    qtracks++;
                    qheader = qblocks;
                    buf[qblocks].h = h;
                    buf[qblocks + 1].p = linmix(seq[j], seq[j + 1], low);
                    buf[qblocks + 2].p = linmix(seq[j], seq[j + 1], high);
                    qblocks += 3;
                }
            }
        }
    }

    // Shrink buffer as necessary
    buf = (TrackBlock *)realloc(buf, qblocks * sizeof(TrackBlock));
    TrackPrivateData *pd = new TrackPrivateData(qtracks, qblocks, buf);
    data = QSharedDataPointer<TrackPrivateData>(pd);
}

TrackData::~TrackData() {}

size_t TrackData::getNBlocks() const {
    if (!data) {
        return 0;
    }
    return data.constData()->nblocks;
}
size_t TrackData::getNTracks() const {
    if (!data) {
        return 0;
    }
    return data.constData()->ntracks;
}
const TrackBlock *TrackData::getBlocks() const {
    if (!data) {
        return NULL;
    }
    return data.constData()->data;
}

const TrackMetaData *TrackData::getMeta() const {
    if (!data) {
        return NULL;
    }
    return data.constData()->meta;
}

void TrackData::calcTimeBounds(double &lower, double &upper) const {
    upper = -kInfinity;
    lower = kInfinity;
    if (!data) {
        return;
    }
    TrackBlock *blocks = data.constData()->data;
    size_t ntracks = data.constData()->ntracks;
    size_t i = 0;
    for (size_t k = 0; k < ntracks; k++) {
        size_t npts = blocks[i].h.npts;
        i++;
        for (size_t j = 0; j < npts; j++) {
            upper = std::fmax(blocks[i].p.time, upper);
            lower = std::fmin(blocks[i].p.time, lower);
            i++;
        }
    }
}
void TrackData::calcEnergyBounds(double &lower, double &upper) const {
    upper = -kInfinity;
    lower = kInfinity;
    if (!data) {
        return;
    }
    TrackBlock *blocks = data.constData()->data;
    size_t ntracks = data.constData()->ntracks;
    size_t i = 0;
    for (size_t k = 0; k < ntracks; k++) {
        size_t npts = blocks[i].h.npts;
        i++;
        for (size_t j = 0; j < npts; j++) {
            float e = blocks[i].p.energy;
            upper = std::fmax(e, upper);
            if (e > 0.) {
                lower = std::fmin(e, lower);
            }
            i++;
        }
    }
}
int TrackData::calcMaxGenerations() const {
    const TrackPrivateData *tpd = data.constData();
    size_t ntracks = tpd->ntracks;
    TrackMetaData *meta = tpd->meta;
    TrackBlock *blocks = tpd->data;
    int maxgen = 0;
    size_t *index = new size_t[ntracks];
    size_t i = 0;
    for (size_t z = 0; z < ntracks; z++) {
        index[z] = i;
        i += blocks[i].h.npts + 1;
    }

    for (size_t z = 0; z < ntracks; z++) {
        size_t current = z;
        int ng = 1;
        while (current > 0 && blocks[index[current]].h.parent_id > 0) {
            size_t ncurrent = (size_t)-1;
            for (size_t y = current; y > 0; y--) {
                int hi = index[y - 1];
                if (blocks[hi].h.track_id ==
                    blocks[index[current]].h.parent_id) {
                    ncurrent = y - 1;
                    break;
                }
            }

            if (ncurrent == (size_t)-1) {
                break;
            }
            current = ncurrent;

            ng++;
        }
        meta[z].generation = (ushort)std::min(ng, 65535);
        maxgen = std::max(ng, maxgen);
    }
    delete[] index;
    return maxgen;
}

static std::map<float, int> mapdedup(const std::vector<float> &a) {
    std::map<float, int> c;
    for (float b : a) {
        if (c.count(b)) {
            c[b] += 1;
        } else {
            c[b] = 1;
        }
    }
    return c;
}

static void constructRangeHistogram(const std::vector<float> &starts,
                                    const std::vector<float> &ends,
                                    const std::vector<float> &spikes,
                                    QVector<QPointF> &pts) {
    std::map<float, int> mup = mapdedup(starts);
    std::map<float, int> mdown = mapdedup(ends);
    std::map<float, int> mdelt = mapdedup(spikes);
    mup[kInfinity] = 0;
    mdown[kInfinity] = 0;
    mdelt[kInfinity] = 0;
    std::map<float, int>::iterator upiter = mup.begin();
    std::map<float, int>::iterator downiter = mdown.begin();
    std::map<float, int>::iterator deltaiter = mdelt.begin();
    int height = 0;
    while (upiter != mup.end() || downiter != mdown.end() ||
           deltaiter != mdelt.end()) {
        // Pick next from all three; add least point
        std::pair<float, int> up = *upiter;
        std::pair<float, int> down = *downiter;
        std::pair<float, int> delta = *deltaiter;
        // Casework by number of ties for first
        float pos = std::min(up.first, std::min(down.first, delta.first));
        if (pos >= kInfinity) {
            break;
        }

        pts.push_back(QPointF(pos, height));
        if (up.first == pos) {
            height += up.second;
            ++upiter;
        }
        if (delta.first == pos) {
            pts.push_back(QPointF(pos, height + delta.second));
            ++deltaiter;
        }
        if (down.first == pos) {
            height -= down.second;
            ++downiter;
        }
        pts.push_back(QPointF(pos, height));
    }
}

void TrackData::constructRangeHistograms(QVector<QPointF> &tp,
                                         QVector<QPointF> &ep, const Range &tr,
                                         const Range &er) const {
    TrackBlock *blocks = data.constData()->data;
    size_t ntracks = data.constData()->ntracks;

    std::vector<float> tstarts, tends, tspikes;
    std::vector<float> estarts, eends, espikes;
    size_t i = 0;
    for (size_t m = 0; m < ntracks; m++) {
        const TrackHeader &h = blocks[i].h;
        i++;
        for (int32_t j = 0; j < h.npts - 1; j++) {
            const TrackPoint &ptA = blocks[i].p;
            const TrackPoint &ptB = blocks[i + 1].p;
            i++;

            float ta = ptA.time, tb = ptB.time;
            float ea = ptA.energy, eb = ptB.energy;

            float stl = std::max(0.f, (float(tr.low) - ta) /
                                          (ta == tb ? 1e38f : tb - ta));
            float sth = std::min(1.f, (float(tr.high) - ta) /
                                          (ta == tb ? 1e38f : tb - ta));
            float sel = std::max(0.f, (float(er.low) - ea) /
                                          (ea == eb ? 1e38f : eb - ea));
            float seh = std::min(1.f, (float(er.high) - ea) /
                                          (ea == eb ? 1e38f : eb - ea));
            float cta, ctb, cea, ceb;
            // restrict one by the other and vica versa
            cta = sel * ta + (1 - sel) * tb;
            ctb = (1 - seh) * ta + seh * tb;
            cea = stl * ea + (1 - stl) * eb;
            ceb = (1 - sth) * ea + sth * eb;
            //            cta = ta;
            //            ctb = tb;
            //            cea = ea;
            //            ceb = eb;
            if (cta == ctb) {
                //                tspikes.push_back(cta/CLHEP::ns);
            } else {
                tstarts.push_back(std::min(cta, ctb) / CLHEP::ns);
                tends.push_back(std::max(cta, ctb) / CLHEP::ns);
            }
            if (cea == ceb) {
                //                espikes.push_back(ea/CLHEP::eV);
            } else {
                estarts.push_back(std::min(cea, ceb) / CLHEP::eV);
                eends.push_back(std::max(cea, ceb) / CLHEP::eV);
            }
        }
        i++;
    }
    constructRangeHistogram(estarts, eends, espikes, ep);
    constructRangeHistogram(tstarts, tends, tspikes, tp);
}

TrackPrivateData::TrackPrivateData(const char *filename) {
    struct stat sb;
    int fd = open(filename, O_RDONLY);
    if (fd == -1) {
        qFatal("Invalid track file, '%s'", filename);
    }
    fstat(fd, &sb);

    mmapbytes = sb.st_size;
    char *buf = (char *)mmap(NULL, mmapbytes, PROT_READ, MAP_SHARED, fd, 0);
    close(fd);

    static_assert(sizeof(TrackPoint) == sizeof(TrackHeader) &&
                      sizeof(TrackPoint) == 32,
                  "Need uniform chunk size");
    nblocks = mmapbytes / sizeof(TrackPoint);

    data = reinterpret_cast<TrackBlock *>(buf);
    ntracks = 0;
    for (size_t i = 0; i < nblocks;) {
        i += data[i].h.npts + 1;
        ntracks++;
    }

    meta = setupBallRadii(data, ntracks);
}

TrackPrivateData::TrackPrivateData(size_t itracks, size_t iblocks,
                                   TrackBlock *idata) {
    ntracks = itracks;
    nblocks = iblocks;
    mmapbytes = 0;
    data = idata;
    meta = setupBallRadii(data, ntracks);
}
TrackPrivateData::TrackPrivateData(const TrackPrivateData &other)
    : QSharedData(other), ntracks(other.ntracks), nblocks(other.nblocks),
      data(other.data), meta(other.meta), mmapbytes(other.mmapbytes) {
    ref.ref();
}
TrackPrivateData::~TrackPrivateData() {
    if (mmapbytes > 0) {
        munmap(data, mmapbytes);
    } else {
        free(data);
    }
    if (meta) {
        delete[] meta;
    }
}

class ElemSort {
public:
    ElemSort(const std::vector<Element> &e) : elts(e) {}
    bool operator()(int i, int j) {
        const Element &l = elts[i];
        const Element &r = elts[j];
        static QCollator coll;
        coll.setNumericMode(true);

        const char *lname = l.name ? l.name.data() : "";
        const char *rname = r.name ? r.name.data() : "";
        return coll.compare(lname, rname) < 0;
    }

    const std::vector<Element> &elts;
};

int convertCreation(std::vector<Element> &elts, const G4VPhysicalVolume *phys,
                    G4RotationMatrix rot, int *counter) {
    int cc = 0;
    if (!counter) {
        counter = &cc;
    }

    G4ThreeVector offset = phys->GetFrameTranslation();
    offset = rot.inverse() * offset;
    const G4RotationMatrix &r = (phys->GetFrameRotation() != NULL)
                                    ? *phys->GetFrameRotation()
                                    : G4RotationMatrix();
    rot = r * rot;
    const G4LogicalVolume *log = phys->GetLogicalVolume();
    const G4Material *mat = log->GetMaterial();

    int k = elts.size();
    elts.push_back(Element());
    {
        // Reference only valid where elts is unmodified
        Element &m = elts[k];
        m.name = phys->GetName();

        // Only identity has a trace of +3 => norm2 of 0
        m.rotated = rot.norm2() > 1e-10;
        m.offset = offset;
        m.rot = rot;

        m.ccode = 0;
        m.material = mat;
        m.solid = log->GetSolid();
        if (0) {
            m.solid = BooleanTree::compile(m.solid);
        }
        m.visible = mat->GetDensity() > 0.01 * CLHEP::g / CLHEP::cm3;
        m.alpha = 0.8; // 1.0;// todo make basic alpha controllable
        m.ecode = *counter;
        (*counter)++;
    }

    std::vector<int> svi;
    for (int i = 0; i < log->GetNoDaughters(); i++) {
        svi.push_back(convertCreation(elts, log->GetDaughter(i), rot, counter));
    }
    ElemSort esort(elts);
    std::sort(svi.begin(), svi.end(), esort);

    elts[k].children = svi;

    return k;
}
