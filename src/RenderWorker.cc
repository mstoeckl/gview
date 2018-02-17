#include "RenderWorker.hh"
#include "LineCollection.hh"

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

static long recursivelySumNCalls(const Element &r, const ElemMutables e[]) {
    long s = e[r.ecode].ngeocalls;
    for (size_t k = 0; k < r.children.size(); k++) {
        s += recursivelySumNCalls(r.children[k], e);
    }
    return s;
}

static void recursivelyPrintNCalls(const Element &r, const ElemMutables e[],
                                   int depth = 0, long net = 0) {
    if (net == 0) {
        net = recursivelySumNCalls(r, e);
    }
    char ws[256];
    int i = 0;
    for (; i < depth && i < 255; i++) {
        ws[i] = ' ';
    }
    ws[i] = '\0';
    qDebug("%s%s %.5f %ld", ws, r.name.data(),
           e[r.ecode].ngeocalls / double(net), e[r.ecode].ngeocalls);
    for (size_t k = 0; k < r.children.size(); k++) {
        recursivelyPrintNCalls(r.children[k], e, depth + 1, net);
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

int traceRay(const QPointF &scpt, const ViewData &d, const Element *hits[],
             Intersection ints[], int maxhits, long iteration,
             ElemMutables mutables[]) {
    // Given a square camera, p in [0,1]X[0,1]. Rectangles increase one dim.
    // Return number of hits recorded.
    const G4ThreeVector forward = d.orientation.rowX();
    const G4ThreeVector updown = d.orientation.rowY();
    const G4ThreeVector lright = d.orientation.rowZ();

    const G4ThreeVector init =
        d.camera + updown * scpt.y() * d.scale + lright * scpt.x() * d.scale;

    // Pseudorandom number to prevent consistent child ordering
    // so as to make overlapping children obvious
    const size_t rotfact = size_t(random());
    // Minimum feature size, yet still above double discrimination threshold
    //    const G4double epsilon = 1e-3 * CLHEP::nanometer;

    G4double sdist, edist;
    G4ThreeVector entrynormal;
    G4ThreeVector exitnormal;
    bool exists = clipRay(d.clipping_planes, init, forward, sdist, edist,
                          entrynormal, exitnormal);
    if (!exists) {
        return 0;
    }

    G4ThreeVector start = init + forward * sdist;
    ++mutables[d.elements.ecode].ngeocalls;
    // Ensure ray starts in the root.
    bool clippable = false;
    if (!d.elements.solid->Inside(start)) {
        G4double jdist = d.elements.solid->DistanceToIn(start, forward);
        if (jdist + sdist >= edist) {
            // Root solid not reachable with available distance
            return 0;
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
    stack[0] = &d.elements;
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
                last.children[(walk + rotfact) % last.children.size()];
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

    int nhits = 1;
    if (clippable) {
        // Start point is in the world volume, and is an intersection
        hits[0] = stack[n - 1];
        ints[0].dist = sdist;
        ints[0].normal = entrynormal;
        ints[0].is_clipping_plane = true;
    } else {
        // Starting point
        ++mutables[d.elements.ecode].ngeocalls;
        G4ThreeVector lnormal = d.elements.solid->SurfaceNormal(local);
        hits[0] = &d.elements;
        ints[0].dist = sdist;
        ints[0].normal = lnormal;
        ints[0].is_clipping_plane = false;
        // Intersection on entry
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
                curr.children[(walk + rotfact) % curr.children.size()];
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
            // End clip (counts as a hit!)
            ints[nhits].dist = edist;
            ints[nhits].normal = exitnormal;
            ints[nhits].is_clipping_plane = true;
            return nhits;
        } else if (exitdist < closestDist) {
            // Transition on leaving vol w/o intersections
            sdist += exitdist;
            G4ThreeVector lpos = offsets[n - 1] + forward * (sdist - isdist);
            // Drop from stack
            --n;

            // Record hit with normal
            ++cmu.ngeocalls;
            G4ThreeVector lnormal =
                curr.solid->SurfaceNormal(condrot(curr, lpos));
            ints[nhits].dist = sdist;
            ints[nhits].normal = condirot(curr, lnormal);
            ints[nhits].is_clipping_plane = false;
            if (n == 0 || nhits >= maxhits) {
                return nhits;
            }
            hits[nhits] = stack[n - 1];
            nhits++;
        } else {
            // Transition on visiting a child
            stack[n] = closest;
            offsets[n] = closest->offset + offsets[n - 1];
            sdist += closestDist;
            n++;
            G4ThreeVector lpos = offsets[n - 1] + forward * (sdist - isdist);

            // Record hit with normal
            ++mutables[closest->ecode].ngeocalls;
            G4ThreeVector lnormal =
                closest->solid->SurfaceNormal(condrot(*closest, lpos));
            ints[nhits].dist = sdist;
            ints[nhits].normal = condirot(*closest, lnormal);
            ints[nhits].is_clipping_plane = false;
            if (nhits >= maxhits) {
                // Don't increment hit counter; ignore last materialc
                return nhits;
            }
            hits[nhits] = closest;
            nhits++;
        }
    }
    return nhits;
}

int compressTraces(const Element *hits[], Intersection ints[], int m) {
    if (m <= 1) {
        return m;
    }
    const G4double epsilon = 0.1 * CLHEP::nm;
    // First, nuke the empty slices
    int n = 0;
    for (int i = 0; i < m; i++) {
        G4double sep = std::abs(ints[i + 1].dist - ints[i].dist);
        if (sep > epsilon) {
            hits[n] = hits[i];
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
        if (i == 0 || hits[i]->ccode != hits[i - 1]->ccode ||
            hits[i]->visible != hits[i - 1]->visible) {
            hits[p] = hits[i];
            ints[p] = ints[i];
            p++;
        }
    }
    ints[p] = ints[n];

    //    qDebug("");
    //    for (int k = 0; k < p; k++) {
    //        G4double sep = std::abs(ints[k + 1].dist - ints[k].dist);
    //        qDebug("%d %d %s %g mm ", k, k + 1, hits[k]->name.data(),
    //               sep / CLHEP::mm);
    //    }

    return p;
}

static QRgb colorMap(const Intersection &intersection,
                     const G4RotationMatrix &frame, const VColor &base,
                     double shade_x, double shade_y) {
    const G4ThreeVector &forward = frame.rowX();
    const G4ThreeVector &normal = intersection.normal;

    // Opposed normals (i.e, for transp backsides) are mirrored
    G4double cx = std::abs(std::acos(-normal * forward) / CLHEP::pi);
    cx = 1.0 - std::max(0.0, std::min(1.0, 0.7 * cx));

    if (intersection.is_clipping_plane) {
        double aslp = (frame.rowY() * forward) + 1.0;
        double bslp = (frame.rowZ() * forward) + 1.0;

        double shade_factor = aslp * shade_x + bslp * shade_y;
        shade_factor *= 20.;
        if (std::fmod(shade_factor, 1.) < 0.5)
            cx *= 0.75;
    }

    return VColor::fromRgbF(cx * base.redF(), cx * base.greenF(),
                            cx * base.blueF())
        .rgb();
}

void countTree(const Element &e, int &treedepth, int &nelements) {
    // Count nelements & depth
    int t = 0;
    int n = 0;
    for (const Element &s : e.children) {
        int at, an;
        countTree(s, at, an);
        t = std::max(at, t);
        n += an;
    }
    treedepth = t + 1;
    nelements = n + 1;
}

static inline double iproject(const G4ThreeVector &point,
                              const G4ThreeVector &dcamera,
                              const G4double &dscale,
                              const G4RotationMatrix &ori, int w, int h,
                              int *dx, int *dy) {
    G4ThreeVector local = ori * (point - dcamera);
    G4double idscale = 1 / dscale;
    double fy = local.y() * idscale;
    double fx = local.z() * idscale;
    double off = local.x();
    int dmind = 2 * std::min(w, h);
    *dx = int(0.5 * w + fx * dmind);
    *dy = int(0.5 * h + fy * dmind);
    return off;
}

static void trackColors(const TrackHeader &h, const TrackPoint &pa,
                        const TrackPoint &pb, const G4ThreeVector &a,
                        const G4ThreeVector &b, const G4ThreeVector &forward,
                        QRgb &ca, QRgb &cb, float &wa, float &wb) {
    G4ThreeVector normal = b - a;
    double coa = (normal * forward) / normal.mag();

    //    double mtime = double(pa.time) * (1. - interp) * double(pb.time) *
    //    double ptime = (int(mtime / CLHEP::ns * 1024.) % 1024) / 1024.;

    double aloge = std::log10(double(pa.energy)) / 2.0;
    double patime = aloge - std::floor(aloge);
    ca = QColor::fromHsvF(patime, 0.8, 1.0 - 0.5 * std::abs(coa)).rgb();

    double bloge = std::log10(double(pb.energy)) / 2.0;
    double pbtime = bloge - std::floor(bloge);
    cb = QColor::fromHsvF(pbtime, 0.8, 1.0 - 0.5 * std::abs(coa)).rgb();

    if (h.ptype == 11) {
        // e-
        wa = 1.0f;
        wb = 1.0f;
    } else if (h.ptype == 22) {
        // gamma
        wa = 0.7f;
        wb = 0.7f;
    } else {
        // others (optiphoton)
        wa = 0.5f;
        wb = 0.5f;
    }
}

RenderRayTask::RenderRayTask(QRect p, RenderGraph &h, QSharedPointer<Context> c,
                             int i)
    : RenderGraphTask(p, h, c, i) {}

void RenderRayTask::run() {
    const ViewData &d = *ctx->viewdata;

    if (!d.elements.solid) {
        return;
    }
    int treedepth;
    int nelements;
    countTree(d.elements, treedepth, nelements);
    // ^ TODO: allocate traceray buffers to match!
    if (treedepth > 10) {
        qFatal("Excessive tree depth, fatal!");
    }
    QImage *next = &(*ctx->image);
    int w = next->width();
    int h = next->height();
#if 0
        QTime t = QTime::currentTime();
        qDebug("render lod %d w %d h %d | %d of %d", d.level_of_detail, w, h, slice,
               nslices);
#endif

    int mind = w > h ? h : w;

    int xl = pixels.left();
    int xh = pixels.right();
    int yl = pixels.top();
    int yh = pixels.bottom();

    const G4double radius = 0.8 / mind;
    const int M = 30;
    const Element *hits[M];
    Intersection ints[M + 1];
    const Element *althits[M];
    Intersection altints[M + 1];
    // Zero construct by default
    ElemMutables *mutables = new ElemMutables[nelements]();
    int iter = 1;

    QRgb *trackcolors = ctx->colors;
    G4double *trackdists = ctx->distances;

    for (int i = yl; i < yh; i++) {
        QRgb *pts = reinterpret_cast<QRgb *>(next->scanLine(i));
        for (int j = xl; j < xh; j++) {
            if (ctx->abort_flag) {
                delete[] mutables;
                return;
            }
            QPointF pt((j - w / 2.) / (2. * mind), (i - h / 2.) / (2. * mind));

            int m = traceRay(pt, d, hits, ints, M, iter, mutables);
            iter++;
            m = compressTraces(hits, ints, m);

            int p = m - 1;
            bool line = false;
            if (d.level_of_detail <= -1) {
                int ndevs[M + 1];
                for (int k = 0; k < M + 1; k++) {
                    ndevs[k] = 0;
                }
                G4double seed = CLHEP::pi / 5 * (qrand() >> 16) / 16384.0;
                for (int k = 0; k < 10; k++) {
                    QPointF off(radius * std::cos(seed + k * CLHEP::pi / 5),
                                radius * std::sin(seed + k * CLHEP::pi / 5));
                    int am = traceRay(pt + off, d, althits, altints, M, iter,
                                      mutables);
                    iter++;
                    am = compressTraces(althits, altints, am);
                    // At which intersection have we disagreements?
                    int cm = std::min(am, m);
                    for (int l = 0; l < cm + 1; l++) {
                        const G4double jump =
                            1.0 * radius * d.scale /
                            std::abs(-altints[l].normal * d.orientation.rowX());
                        bool diffmatbehind =
                            ((l < cm) && althits[l]->ccode != hits[l]->ccode);
                        bool edgediff =
                            (std::abs(altints[l].normal * ints[l].normal) <
                                 0.3 || // Q: abs?
                             std::abs(altints[l].dist - ints[l].dist) > jump);
                        if (diffmatbehind || edgediff) {
                            ++ndevs[l];
                        }
                    }
                    // Differing intersection count
                    for (int l = am; l < m; l++) {
                        ++ndevs[l];
                    }
                    for (int l = m; l < am; l++) {
                        ++ndevs[l];
                    }
                }
                const int devthresh = 2;
                for (int k = 0; k < m + 1; ++k) {
                    if (ndevs[k] >= devthresh) {
                        p = k - 1;
                        line = true;
                        break;
                    }
                }
            }

            int sidx = i * w + j;
            QRgb trackcol = trackcolors[sidx];
            // ^ trackcol includes background
            G4double trackdist = trackdists[sidx];

            bool rayoverride = false;
            for (int k = 0; k < p; k++) {
                if (ints[k].dist > trackdist) {
                    p = k - 1;
                    // WARNING: there exist exceptions!
                    // NOT QUITE ACCURATE WITH LAYERING! TODO FIXME
                    rayoverride = true;
                    break;
                }
            }

            QRgb col = (line && !rayoverride) ? qRgb(0, 0, 0) : trackcol;
            // p indicates the first volume to use the color rule on
            // (p<0 indicates the line dominates)
            for (int k = p; k >= 0; --k) {
                // We use the intersection before the volume
                const VColor &base_color = d.color_table[hits[k]->ccode];
                QRgb altcol = colorMap(ints[k], d.orientation, base_color,
                                       i / (double)mind, j / (double)mind);
                double e = hits[k]->alpha, f = (1. - hits[k]->alpha);
                if (!hits[k]->visible) {
                    continue;
                }
                col = qRgba(int(e * qRed(altcol) + f * qRed(col)),
                            int(e * qGreen(altcol) + f * qGreen(col)),
                            int(e * qBlue(altcol) + f * qBlue(col)),
                            int(e * qAlpha(altcol) + f * qAlpha(col)));
            }
            pts[j] = col;
        }
        QMetaObject::invokeMethod(&home, "progUpdate", Qt::QueuedConnection,
                                  Q_ARG(int, yh - yl + 1),
                                  Q_ARG(int, ctx->renderno));
    }
    if (d.level_of_detail <= -1) {
        recursivelyPrintNCalls(d.elements, mutables);
    }

    delete[] mutables;
    QMetaObject::invokeMethod(&home, "queueNext", Qt::QueuedConnection,
                              Q_ARG(int, this->id), Q_ARG(int, ctx->renderno));
}

RenderTrackTask::RenderTrackTask(QRect p, RenderGraph &h,
                                 QSharedPointer<Context> c, int i)
    : RenderGraphTask(p, h, c, i) {}

void RenderTrackTask::run() {
    int xl = pixels.left();
    int xh = pixels.right();
    int yl = pixels.top();
    int yh = pixels.bottom();

    double *dists = ctx->distances;
    QRgb *colors = ctx->colors;
    int w = ctx->image->width();
    int h = ctx->image->height();
    const ViewData &d = *ctx->viewdata;
    const TrackData &t = d.tracks;

    for (int y = yl; y < yh; y++) {
        for (int x = xl; x < xh; x++) {
            // Scale up to be beyond anything traceRay produces
            dists[y * w + x] = 4 * kInfinity;
            colors[y * w + x] = qRgba(255, 255, 255, 0.);
        }
    }

    size_t ntracks = t.getNTracks();
    const TrackHeader *headers = t.getHeaders();
    const TrackPoint *points = t.getPoints();

    int mind = w > h ? h : w;
    float rad = std::max(mind / 800.0f, 0.01f);

    for (size_t i = 0; i < ntracks; i++) {
        if (ctx->abort_flag) {
            return;
        }
        if (i % ((ntracks / 20) + 1) == 0) {
            QMetaObject::invokeMethod(&home, "progUpdate", Qt::QueuedConnection,
                                      Q_ARG(int, 20 > ntracks ? ntracks : 20),
                                      Q_ARG(int, ctx->renderno));
        }

        const TrackHeader &header = headers[i];
        const TrackPoint *tp = &points[header.offset];

        // Project calculations 1 point ahead
        TrackPoint sp = tp[0];
        G4ThreeVector spos(sp.x, sp.y, sp.z);
        int sdx, sdy;
        double soff;
        soff =
            iproject(spos, d.camera, d.scale, d.orientation, w, h, &sdx, &sdy);

        // Fast clipping for when the track is way out of view
        int rr = int(2 * mind * header.bballradius / d.scale);
        if (sdx + rr <= xl || sdx - rr > xh || sdy + rr <= yl ||
            sdy - rr > yh) {
            continue;
        }

        for (int j = 1; j < header.npts; j++) {
            // Adopt previous late point as early point
            TrackPoint ep = sp;
            int edx = sdx, edy = sdy;
            double eoff = soff;
            G4ThreeVector epos = spos;
            // Calculate new late point
            sp = tp[j];
            spos = G4ThreeVector(sp.x, sp.y, sp.z);
            soff = iproject(spos, d.camera, d.scale, d.orientation, w, h, &sdx,
                            &sdy);

            if ((edx < xl && sdx < xl) || (edy < yl && sdy < yl) ||
                (edx >= xh && sdx >= xh) || (edy >= yh && sdy >= yh)) {
                continue;
            }

            // NOTE: single point lines are acceptable
            float dy, dx;
            int n;
            if (std::abs(sdx - edx) > std::abs(sdy - edy)) {
                n = std::abs(sdx - edx) + 1;
                dx = sdx - edx > 0 ? 1.0 : -1.0;
                dy = float(sdy - edy) / float(n);
            } else {
                n = std::abs(sdy - edy) + 1;
                dy = sdy - edy > 0 ? 1.0 : -1.0;
                dx = float(sdx - edx) / float(n);
            }
            float ix = edx;
            float iy = edy;

            int uxl = dx == 0.0f ? 0 : int((xl - edx) / dx);
            int uyl = dy == 0.0f ? 0 : int((yl - edy) / dy);
            int uxh = dx == 0.0f ? INT_MAX : int((xh - edx) / dx);
            int uyh = dy == 0.0f ? INT_MAX : int((yh - edy) / dy);

            // Extra buffers are just in case
            int near = std::max(std::min(uxl, uxh), std::min(uyl, uyh)) - 1;
            int far = std::min(std::max(uxl, uxh), std::max(uyl, uyh)) + 1;

            if (std::max(near, 0) > std::min(far, n)) {
                continue;
            }

            float w1, w2;
            QRgb c1, c2;
            trackColors(header, ep, sp, epos, spos, d.orientation.rowX(), c1,
                        c2, w1, w2);

            for (int s = std::max(near, 0); s < std::min(far, n); s++) {
                float x = ix + s * dx;
                float y = iy + s * dy;

                float b = n >= 2 ? s / float(n - 1) : 0.5f;
                float bc = 1.f - b;
                double off = eoff * double(bc) + soff * double(b);
                float rws = w1 * (1.0f - b) + w2 * b;
                QRgb col = qRgb(int(bc * qRed(c1) + b * qRed(c2)),
                                int(bc * qGreen(c1) + b * qGreen(c2)),
                                int(bc * qBlue(c1) + b * qBlue(c2)));

                /* Draw circle */
                float spot = rad * rws;
                int ly = int(std::round(y - spot));
                int hy = int(std::round(y + spot));
                for (int ddy = std::max(yl, ly); ddy <= std::min(yh - 1, hy);
                     ddy++) {
                    float lr = std::sqrt(std::max(
                        0.f, spot * spot - float((ddy - y) * (ddy - y))));
                    int lx = int(std::round(x - lr));
                    int hx = int(std::round(x + lr));
                    for (int ddx = std::max(xl, lx);
                         ddx <= std::min(xh - 1, hx); ddx++) {
                        int sidx = ddy * w + ddx;
                        if (dists[sidx] > off) {
                            dists[sidx] = off;
                            colors[sidx] = col;
                        }
                    }
                }
            }
        }
    }

    QMetaObject::invokeMethod(&home, "queueNext", Qt::QueuedConnection,
                              Q_ARG(int, this->id), Q_ARG(int, ctx->renderno));
}

static void setupBallRadii(TrackHeader *headers, const TrackPoint *points,
                           size_t ntracks) {
    for (size_t i = 0; i < ntracks; i++) {
        TrackHeader &h = headers[i];
        const TrackPoint &t0 = points[h.offset];
        G4ThreeVector p0(t0.x, t0.y, t0.z);
        G4double radius = 0.0;
        for (size_t j = h.offset + 1; j < h.offset + size_t(h.npts); j++) {
            const TrackPoint &tj = points[j];
            G4ThreeVector pj(tj.x, tj.y, tj.z);
            radius = std::max(radius, (pj - p0).mag());
        }
        h.bballradius = radius;
    }
}

TrackData::TrackData() {}

TrackData::TrackData(const char *filename) {
    FILE *f = fopen(filename, "rb");
    fseek(f, 0, SEEK_END);
    size_t fsize = size_t(ftell(f));
    fseek(f, 0, SEEK_SET);
    char *buf = reinterpret_cast<char *>(malloc(fsize));
    fread(buf, fsize, 1, f);
    fclose(f);

    size_t blocksize = sizeof(TrackPoint);
    size_t n = fsize / blocksize;
    size_t ntracks = 0;
    for (size_t i = 0; i < n;) {
        const TrackHeader &h =
            *reinterpret_cast<TrackHeader *>(&buf[blocksize * i]);
        i += 1 + size_t(h.npts);
        ntracks++;
    }

    TrackPrivateData *pd = new TrackPrivateData(ntracks, n - ntracks);
    TrackHeader *headers = pd->headers;
    TrackPoint *points = pd->points;
    size_t j = 0;
    for (size_t i = 0; i < ntracks; i++) {
        const TrackHeader &h =
            *reinterpret_cast<TrackHeader *>(&buf[blocksize * j]);
        headers[i].npts = h.npts;
        headers[i].ptype = h.ptype;
        headers[i].offset = j - i;
        size_t start = j + size_t(h.npts);
        for (; j < start; j++) {
            if (j >= n - 1) {
                /* Truncated ! */
                headers[i].npts = j - i - headers[i].offset;
                qWarning("Track file truncated at %lu leaving run length %d",
                         j + 1, headers[i].npts);
                break;
            }
            points[j - i] =
                *reinterpret_cast<TrackPoint *>(&buf[blocksize * (j + 1)]);
        }
        j++;
    }
    free(buf);
    //    pd->tree = new LineCollection(headers, points, ntracks);
    pd->tree = NULL;
    setupBallRadii(headers, points, ntracks);

    data = QSharedDataPointer<TrackPrivateData>(pd);
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

TrackData::TrackData(const TrackData &other, ViewData &vd, Range trange,
                     Range erange, IRange selidxs) {
    size_t otracks = other.getNTracks();
    const TrackHeader *oheaders = other.getHeaders();
    const TrackPoint *opoints = other.getPoints();

    std::vector<TrackHeader> ch;
    std::vector<TrackPoint> cp;
    // Should yield extra space except for pathological cases
    ch.reserve(other.getNTracks());
    cp.reserve(other.getNPoints());
    float tlow = float(trange.low), thigh = float(trange.high);
    float elow = float(erange.low), ehigh = float(erange.high);
    size_t nlow = std::max(size_t(0), selidxs.low - 1);
    size_t nhigh = std::min(otracks, selidxs.high);
    for (size_t i = nlow; i < nhigh; i++) {
        const TrackPoint *seq = &opoints[oheaders[i].offset];
        int npts = oheaders[i].npts;
        G4ThreeVector fts(seq[0].x, seq[0].y, seq[0].z);
        bool intime = seq[0].time >= tlow && seq[0].time <= thigh;
        bool inenergy = seq[0].energy >= elow && seq[0].energy <= ehigh;
        bool inconvex = insideConvex(vd.clipping_planes, fts);
        bool started = false;
        if (inconvex && intime && inenergy) {
            TrackHeader h;
            h.offset = cp.size();
            h.ptype = oheaders[i].ptype;
            h.npts = 1;
            ch.push_back(h);
            cp.push_back(seq[0]);
            started = true;
        }

        for (int j = 0; j < npts - 1; j++) {
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
                    TrackHeader h;
                    h.offset = cp.size();
                    h.ptype = oheaders[i].ptype;
                    h.npts = 1;
                    ch.push_back(h);
                    cp.push_back(seq[j]);
                }
                ch[ch.size() - 1].npts++;
                cp.push_back(seq[j + 1]);
            } else if (low <= 0. && high < 1.) {
                if (high > 0.) {
                    if (!started) {
                        TrackHeader h;
                        h.offset = cp.size();
                        h.ptype = oheaders[i].ptype;
                        h.npts = 1;
                        ch.push_back(h);
                        cp.push_back(seq[j]);
                    }
                    ch[ch.size() - 1].npts++;
                    cp.push_back(linmix(seq[j], seq[j + 1], high));
                } else {
                    // Not there entirely
                }
            } else if (low > 0. && high >= 1.) {
                if (low < 1.) {
                    // Starting new sequence
                    TrackHeader h;
                    h.offset = cp.size();
                    h.ptype = oheaders[i].ptype;
                    h.npts = 2;
                    ch.push_back(h);
                    cp.push_back(linmix(seq[j], seq[j + 1], low));
                    cp.push_back(seq[j + 1]);
                } else {
                    // Not there entirely
                }
            } else {
                if (high < 1. && low > 0. && low < high) {
                    // low > 0, high < 1
                    // A whole new short segment...
                    TrackHeader h;
                    h.offset = cp.size();
                    h.ptype = oheaders[i].ptype;
                    h.npts = 2;
                    ch.push_back(h);
                    cp.push_back(linmix(seq[j], seq[j + 1], low));
                    cp.push_back(linmix(seq[j], seq[j + 1], high));
                }
            }
        }
    }

    size_t qtracks = ch.size();
    size_t qpoints = cp.size();
    TrackPrivateData *pd = new TrackPrivateData(qtracks, qpoints);
    TrackHeader *headers = pd->headers;
    TrackPoint *points = pd->points;
    memcpy(headers, ch.data(), sizeof(TrackHeader) * qtracks);
    memcpy(points, cp.data(), sizeof(TrackPoint) * qpoints);
    //    pd->tree = new LineCollection(headers, points, qtracks);
    pd->tree = NULL;
    setupBallRadii(headers, points, qtracks);
    data = QSharedDataPointer<TrackPrivateData>(pd);
}

TrackData::~TrackData() {}

size_t TrackData::getNPoints() const {
    if (!data) {
        return 0;
    }
    return data.constData()->npoints;
}
size_t TrackData::getNTracks() const {
    if (!data) {
        return 0;
    }
    return data.constData()->ntracks;
}
const TrackHeader *TrackData::getHeaders() const {
    if (!data) {
        return NULL;
    }
    return data.constData()->headers;
}
const TrackPoint *TrackData::getPoints() const {
    if (!data) {
        return NULL;
    }
    return data.constData()->points;
}
const LineCollection *TrackData::getTree() const {
    if (!data) {
        return NULL;
    }
    return data.constData()->tree;
}

void TrackData::calcTimeBounds(double &lower, double &upper) const {
    upper = -kInfinity;
    lower = kInfinity;
    if (!data) {
        return;
    }
    TrackPoint *pts = data.constData()->points;
    size_t npts = data.constData()->npoints;
    for (size_t i = 0; i < npts; i++) {
        upper = std::fmax(pts[i].time, upper);
        lower = std::fmin(pts[i].time, lower);
    }
}
void TrackData::calcEnergyBounds(double &lower, double &upper) const {
    upper = -kInfinity;
    lower = kInfinity;
    if (!data) {
        return;
    }
    TrackPoint *pts = data.constData()->points;
    size_t npts = data.constData()->npoints;
    for (size_t i = 0; i < npts; i++) {
        float e = pts[i].energy;
        upper = std::fmax(e, upper);
        if (e > 0.) {
            lower = std::fmin(pts[i].energy, lower);
        }
    }
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
    TrackPoint *pts = data.constData()->points;
    TrackHeader *headers = data.constData()->headers;
    size_t nheaders = data.constData()->ntracks;
    std::vector<float> tstarts, tends, tspikes;
    std::vector<float> estarts, eends, espikes;
    for (size_t i = 0; i < nheaders; i++) {
        const TrackHeader &h = headers[i];
        for (size_t j = i; j < i + h.npts - 1; j++) {
            float ta = pts[j].time, tb = pts[j + 1].time;
            float ea = pts[j].energy, eb = pts[j + 1].energy;
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
    }
    constructRangeHistogram(estarts, eends, espikes, ep);
    constructRangeHistogram(tstarts, tends, tspikes, tp);
}

TrackPrivateData::TrackPrivateData(size_t itracks, size_t ipoints) {
    ntracks = itracks;
    npoints = ipoints;
    headers = new TrackHeader[ntracks];
    points = new TrackPoint[npoints];
    tree = NULL;
}
TrackPrivateData::TrackPrivateData(const TrackPrivateData &other)
    : QSharedData(other), ntracks(other.ntracks), npoints(other.npoints),
      headers(other.headers), points(other.points) {
    ref.ref();
}
TrackPrivateData::~TrackPrivateData() {
    delete[] headers;
    delete[] points;
    if (tree) {
        delete tree;
    }
}

static bool sortbyname(const Element &l, const Element &r) {
    static QCollator coll;
    coll.setNumericMode(true);

    const char *lname = l.name ? l.name.data() : "";
    const char *rname = r.name ? r.name.data() : "";
    return coll.compare(lname, rname) < 0;
}

Element convertCreation(const G4VPhysicalVolume *phys, G4RotationMatrix rot,
                        int *counter) {
    int cc = 0;
    if (!counter) {
        counter = &cc;
    }
    Element m;
    m.name = phys->GetName();

    G4ThreeVector offset = phys->GetFrameTranslation();
    offset = rot.inverse() * offset;
    const G4RotationMatrix &r = (phys->GetFrameRotation() != NULL)
                                    ? *phys->GetFrameRotation()
                                    : G4RotationMatrix();
    rot = r * rot;

    // Only identity has a trace of +3 => norm2 of 0
    m.rotated = rot.norm2() > 1e-10;
    m.offset = offset;
    m.rot = rot;

    const G4LogicalVolume *log = phys->GetLogicalVolume();
    const G4Material *mat = log->GetMaterial();
    m.ccode = 0;
    m.material = mat;
    m.solid = log->GetSolid();
    m.visible = mat->GetDensity() > 0.01 * CLHEP::g / CLHEP::cm3;
    m.alpha = 0.8; // 1.0;// todo make basic alpha controllable
    m.ecode = *counter;
    (*counter)++;

    m.children = std::vector<Element>();
    for (int i = 0; i < log->GetNoDaughters(); i++) {
        m.children.push_back(
            convertCreation(log->GetDaughter(i), m.rot, counter));
    }
    std::sort(m.children.begin(), m.children.end(), sortbyname);
    return m;
}
