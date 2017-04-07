#include "RenderWorker.hh"
#include "LineCollection.hh"

#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VSolid.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VisExtent.hh>

#include <QCoreApplication>
#include <QProgressDialog>

static long recursivelySumNCalls(const Element &r) {
    long s = r.ngeocalls;
    for (size_t k = 0; k < r.children.size(); k++) {
        s += recursivelySumNCalls(r.children[k]);
    }
    return s;
}

static void recursivelyPrintNCalls(const Element &r, int depth = 0,
                                   long net = 0) {
    if (net == 0) {
        net = recursivelySumNCalls(r);
    }
    char ws[256];
    int i = 0;
    for (; i < depth && i < 255; i++) {
        ws[i] = ' ';
    }
    ws[i] = '\0';
    qDebug("%s%s %.5f %ld", ws, r.name.data(), r.ngeocalls / double(net),
           r.ngeocalls);
    for (size_t k = 0; k < r.children.size(); k++) {
        recursivelyPrintNCalls(r.children[k], depth + 1, net);
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

int traceRay(const QPointF &scpt, const ViewData &d, const Element *hits[],
             Intersection ints[], int maxhits, int iteration) {
    // Given a square camera, p in [0,1]X[0,1]. Rectangles increase one dim.
    // Return number of hits recorded.
    G4ThreeVector forward = d.orientation.rowX();
    G4ThreeVector updown = d.orientation.rowY();
    G4ThreeVector lright = d.orientation.rowZ();

    G4ThreeVector init =
        d.camera + updown * scpt.y() * d.scale + lright * scpt.x() * d.scale;

    // Pseudorandom number to prevent consistent child ordering
    // so as to make overlapping children obvious
    const size_t rotfact = size_t(random());
    // Minimum feature size, yet still above double discrimination threshold
    const G4double epsilon = 1e-3 * CLHEP::nanometer;

    G4double sdist = 0.0;
    G4double edist = kInfinity;
    G4ThreeVector entrynormal;
    G4ThreeVector exitnormal;
    for (uint i = 0; i < d.clipping_planes.size(); i++) {
        const Plane &p = d.clipping_planes[i];
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
                return 0;
            }
        }
    }
    if (edist < sdist) {
        // Planes clip everything out
        return 0;
    }

    G4ThreeVector start = init + forward * sdist;
    ++d.root.ngeocalls;
    // Ensure ray starts in the root.
    bool clippable = false;
    if (!d.root.solid->Inside(start)) {
        G4double jdist = d.root.solid->DistanceToIn(start, forward);
        if (jdist >= kInfinity) {
            // Root solid not reachable
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
    stack[0] = &d.root;
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
            ++elem.ngeocalls;
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
    G4ThreeVector positions[maxdepth];
    positions[0] = start;
    for (size_t i = 1; i < n; i++) {
        const Element &elem = *stack[i];
        positions[i] = positions[i - 1] + elem.offset;
    }
    int nhits = 1;
    if (clippable) {
        // Start point is in the world volume, and is an intersection
        hits[0] = stack[n - 1];
        ints[0].dist = sdist;
        ints[0].normal = entrynormal;
    } else {
        // Starting point
        ++d.root.ngeocalls;
        G4ThreeVector lnormal = d.root.solid->SurfaceNormal(local);
        hits[0] = &d.root;
        ints[0].dist = sdist;
        ints[0].normal = lnormal;
        // Intersection on entry
    }

    for (int iter = 0; iter < 1000; iter++) {
        const Element &curr = *stack[n - 1];
        const G4ThreeVector &pos = positions[n - 1];

        ++curr.ngeocalls;
        if (curr.niter != iteration || curr.abs_dist <= sdist /* epsilon ? */) {

            curr.abs_dist = sdist +
                            curr.solid->DistanceToOut(condrot(curr, pos),
                                                      condrot(curr, forward));
            curr.niter = iteration;
        }
        G4double exitdist = curr.abs_dist - sdist;

        G4double closestDist = 2 * kInfinity;
        const Element *closest = NULL;
        G4ThreeVector closestPos;
        for (size_t walk = 0; walk < curr.children.size(); walk++) {
            const Element &elem =
                curr.children[(walk + rotfact) % curr.children.size()];
            G4ThreeVector sub = pos + elem.offset;

            ++elem.ngeocalls;
            if (elem.niter != iteration || elem.abs_dist <= sdist) {
                elem.abs_dist = sdist +
                                elem.solid->DistanceToIn(
                                    condrot(elem, sub), condrot(elem, forward));
                elem.niter = iteration;
            }
            G4double altdist = elem.abs_dist - sdist;
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
            return nhits;
        } else if (exitdist < closestDist) {
            // Transition on leaving vol w/o intersections
            // Advance all levels, why not
            for (size_t j = 0; j < n; j++) {
                positions[j] += exitdist * forward;
            }
            sdist += exitdist;
            G4ThreeVector lpos = positions[n - 1];
            // Drop from stack
            --n;

            // Record hit with normal
            ++curr.ngeocalls;
            G4ThreeVector lnormal =
                curr.solid->SurfaceNormal(condrot(curr, lpos));
            ints[nhits].dist = sdist;
            ints[nhits].normal = condirot(curr, lnormal);
            if (n == 0 || nhits >= maxhits) {
                return nhits;
            }
            hits[nhits] = stack[n - 1];
            nhits++;
        } else {
            // Transition on visiting a child
            stack[n] = closest;
            positions[n] = closestPos;
            for (size_t j = 0; j < n + 1; j++) {
                positions[j] += closestDist * forward;
            }
            sdist += closestDist;
            n++;

            // Record hit with normal
            ++closest->ngeocalls;
            G4ThreeVector lnormal = closest->solid->SurfaceNormal(
                condrot(*closest, positions[n - 1]));
            ints[nhits].dist = sdist;
            ints[nhits].normal = condirot(*closest, lnormal);
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
        if (i == 0 || hits[i]->mat != hits[i - 1]->mat ||
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

static QRgb colorMap(const G4ThreeVector &normal, const G4ThreeVector &forward,
                     G4double hue, G4double dist) {
    // Opposed normals (i.e, for transp backsides) are mirrored
    G4double cx = 0.7 * std::abs(std::acos(-normal * forward) / CLHEP::pi);
    cx = std::max(0.0, std::min(1.0, cx));
    // note: need faster, less crashy hsv->rgb
    Q_UNUSED(dist);
    return QColor::fromHsvF(hue, 1.0, 1.0 - cx).rgb();
}

static int prepareTree(const Element &e) {
    int d = 0;
    e.niter = 0;
    e.abs_dist = 0.0;
    for (const Element &s : e.children) {
        d = std::max(prepareTree(s), d);
    }
    return d + 1;
}

static void sliceBounds(int slice, int nslices, int w, int h, int &xl, int &xh,
                        int &yl, int &yh) {
    // NOTE: tile system; step 1: prevent partial images at
    // all costs! step 2: use different tile systems for different
    // render layers. Requires a request graph/node system with
    // one entry point; is 1+N control threads :-(
    xl = 0;
    xh = w;
    yl = (slice * h / nslices);
    yh = slice == nslices - 1 ? h : ((slice + 1) * h / nslices);
}

RenderWorker::RenderWorker() { abort_task = false; }

RenderWorker::~RenderWorker() {}

void RenderWorker::coAbort() { abort_task = true; }
void RenderWorker::flushAbort() { abort_task = false; }

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

static QRgb trackColor(const TrackHeader &h, const TrackPoint &pa,
                       const TrackPoint &pb, double interp,
                       const G4ThreeVector &a, const G4ThreeVector &b,
                       const G4ThreeVector &forward) {
    Q_UNUSED(h);
    G4ThreeVector normal = b - a;
    double coa = (normal * forward) / normal.mag();
    // Q fails in fwd dir?? # of tracks changes? raw data modified?
    //    double mtime = double(pa.time) * (1. - interp) * double(pb.time) *
    //    interp;
    double menergy =
        double(pa.energy) * (1. - interp) * double(pb.energy) * interp;
    //    double ptime = (int(mtime / CLHEP::ns * 1024.) % 1024) / 1024.;
    //    double ptime = std::
    double loge = std::log10(menergy) / 2.0;
    double ptime = loge - std::floor(loge);
    return QColor::fromHsvF(ptime, 0.8, 1.0 - 0.5 * std::abs(coa)).rgb();
}

static bool trackTrace(const QPointF &scpt, const ViewData &d,
                       const TrackData &t, G4double &hitdist, QRgb &color) {
    // Ray tracing over all the lines. Can get nasty where
    // we have the source point, but that only happens once.
    // a median-based tree partitioner ought to resolve
    // the source point issue, although one lucky seed
    // cell gets 10K offsets given 200K lines...
    // -> 10K iterations, although all of them are
    // much cheaper...

    // TO test _cylinder_ intersection & perp dist,
    // (a -> b); rad; init; fwd

    // dy = ((b - a) * init) / ((b - a)*fwd)
    // if (dy notin [0,1]) skip
    // p = init - (b-a)*dy
    // closest point approach: min: ((a-p)+bt)^2, has t: (2*(a-p)/b^2) t = 0
    // (then t: is dist; min radius is ((a-p)+bt)^2, and we filter
    // have: p (along b-a),t (along fwd),r (along cross) (basically, a
    // coordinate transform)
    // and thus can decide yes/no.

    // Proper order: Cross-test first, to filter by R. Then filter by P;
    // Finally, calculate T.

    G4ThreeVector forward = d.orientation.rowX();
    G4ThreeVector updown = d.orientation.rowY();
    G4ThreeVector lright = d.orientation.rowZ();
    G4ThreeVector init =
        d.camera + updown * scpt.y() * d.scale + lright * scpt.x() * d.scale;

    const LineCollection *lc = t.getTree();
    const TrackHeader *headers = t.getHeaders();
    const TrackPoint *points = t.getPoints();
    SimplexIterator iter = lc->followRay(init, forward);
    do {
        const std::vector<std::pair<size_t, int>> &pairs =
            iter.getContainedLines();
        G4double nearest = kInfinity;
        std::pair<size_t, int> best;
        //        qDebug("pairs %lu",pairs.size());
        for (const std::pair<size_t, int> &p : pairs) {
            // Pt A, Pt b:
            const TrackHeader &h = headers[p.first];
            const TrackPoint &a = points[h.offset + size_t(p.first)];
            const TrackPoint &b = points[h.offset + size_t(p.first) + 1];
            G4ThreeVector pa(a.x, a.y, a.z);
            G4ThreeVector pb(b.x, b.y, b.z);
            G4ThreeVector dir = pb - pa;
            G4ThreeVector perp = dir.cross(forward);
            G4double rad = perp * (init - pa);
            if (std::abs(rad) < 1.0 * CLHEP::cm) {
                // TODO calculate distance to center plane!
                nearest = 10.0 * CLHEP::m;
                best = p;
            }

            // TODO: detect intersection, requires either 3 projections
            // or the use of a basis
        }
        if (nearest < kInfinity) {
            hitdist = nearest;
            color = qRgb(200, 100, 150);
            return true;
        }
    } while (iter.advance());

    // No intersection
    return false;
}

bool RenderWorker::renderTracks(const ViewData &d, const TrackData &t,
                                G4double *dists, QRgb *colors, int slice,
                                int nslices, int w, int h) {
    // Fixme: the lower edge sometimes displays blank
    // Exact cause unknown

    // TODO: apply line collection and resolve tracks
    // via raytracing...

    int xl, xh, yl, yh;
    sliceBounds(slice, nslices, w, h, xl, xh, yl, yh);
    for (int s = 0; s < (yh - yl) * (xh - xl); s++) {
        dists[s] = kInfinity;
        double r = 255. * s / ((yh - yl) * (xh - xl) - 1);
        Q_UNUSED(r);
        colors[s] = qRgb(255, 255, 255 /*r*/);
    }
    // temp fast abort...
    //    return true;

    size_t ntracks = t.getNTracks();
    const TrackHeader *headers = t.getHeaders();
    const TrackPoint *points = t.getPoints();

    int mind = w > h ? h : w;
    float rad = std::max(mind / 800.0f, 0.5f);
    // ^ ensure minimum 1 pixel brush width
    float r2 = rad * rad;

    for (size_t i = 0; i < ntracks; i++) {
        //        if (this->abort_task) {
        //            this->abort_task = false;
        //            emit aborted();
        //            return false;
        //        }

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

            for (int s = std::max(near, 0); s < std::min(far, n); s++) {
                float x = ix + s * dx;
                float y = iy + s * dy;
                double b = n >= 2 ? s / double(n - 1) : 0.5;
                // Technically, this is a fast circle rendering technique!
                QRgb col = qRgb(0, 0, 0);
                bool hascol = false;
                double off = eoff * (1 - b) + soff * b;
                int ur = int(std::ceil(rad));
                int cy = int(std::round(y));
                for (int ddy = std::max(yl, cy - ur);
                     ddy <= std::min(yh - 1, cy + ur); ddy++) {
                    float lr = std::sqrt(r2 - float((ddy - y) * (ddy - y)));
                    int lx = int(std::round(x - lr));
                    int hx = int(std::round(x + lr));
                    for (int ddx = std::max(xl, lx);
                         ddx <= std::min(xh - 1, hx); ddx++) {
                        int sidx = (ddy - yl) * (xh - xl) + (ddx - xl);
                        if (dists[sidx] > off) {
                            dists[sidx] = off;
                            if (!hascol) {
                                col = trackColor(header, ep, sp, b, epos, spos,
                                                 d.orientation.rowX());
                                hascol = true;
                            }
                            colors[sidx] = col;
                        }
                    }
                }
            }
        }
    }
    return true;
}

bool RenderWorker::render(ViewData d, TrackData t, QImage *next, int slice,
                          int nslices, QProgressDialog *progdiag) {
    if (!d.root.solid) {
        return false;
    }
    int treedepth = prepareTree(d.root);
    // ^ TODO: allocate traceray buffers to match!
    if (treedepth > 10) {
        qFatal("Excessive tree depth, fatal!");
    }
    int w = next->width();
    int h = next->height();
#if 0
    QTime t = QTime::currentTime();
    qDebug("render lod %d w %d h %d | %d of %d", d.level_of_detail, w, h, slice,
           nslices);
#endif

    int mind = w > h ? h : w;

    int xl, xh, yl, yh;
    sliceBounds(slice, nslices, w, h, xl, xh, yl, yh);

    const G4double radius = 0.8 / mind;
    const int M = 30;
    const Element *hits[M];
    Intersection ints[M + 1];
    const Element *althits[M];
    Intersection altints[M + 1];
    int iter = 0;
    G4double *trackdists = new G4double[(xh - xl) * (yh - yl)];
    QRgb *trackcolors = new QRgb[(xh - xl) * (yh - yl)];
    bool s = renderTracks(d, t, trackdists, trackcolors, slice, nslices, w, h);
    if (!s) {
        // aborted
        delete[] trackdists;
        delete[] trackcolors;
        return false;
    }

    for (int i = yl; i < yh; i++) {
        QRgb *pts = reinterpret_cast<QRgb *>(next->scanLine(i));
        for (int j = xl; j < xh; j++) {
            //            if (this->abort_task) {
            //                this->abort_task = false;
            //                emit aborted();
            //                delete[] trackdists;
            //                delete[] trackcolors;
            //                return false;
            //            }
            QPointF pt((j - w / 2.) / (2. * mind), (i - h / 2.) / (2. * mind));

            int m = traceRay(pt, d, hits, ints, M, iter);
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
                int am;
                for (int k = 0; k < 10; k++) {

                    QPointF off(radius * std::cos(seed + k * CLHEP::pi / 5),
                                radius * std::sin(seed + k * CLHEP::pi / 5));
                    am = traceRay(pt + off, d, althits, altints, M, iter);
                    iter++;
                    am = compressTraces(althits, altints, am);
                    // At which intersection have we disagreements?
                    int cm = std::min(am, m);
                    for (int l = 0; l < cm + 1; l++) {
                        const G4double jump =
                            1.0 * radius * d.scale /
                            std::abs(-altints[l].normal * d.orientation.rowX());
                        bool diffmatbehind =
                            ((l < cm) && althits[l]->mat != hits[l]->mat);
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

            int sidx = (i - yl) * (xh - xl) + (j - xl);
            QRgb trackcol = trackcolors[sidx];
            // ^ trackcol includes background
            G4double trackdist = trackdists[sidx];
            //            QRgb trackcol;
            //            G4double trackdist;
            //            bool hit = trackTrace(pt, d, t, trackdist, trackcol);
            //            if (!hit) {
            //                trackcol = qRgb(255,255,255);
            //                trackdist = kInfinity;
            //            }

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
                QRgb altcol = colorMap(ints[k].normal, d.orientation.rowX(),
                                       hits[k]->hue, ints[k].dist);
                double e = hits[k]->alpha, f = (1. - hits[k]->alpha);
                if (!hits[k]->visible) {
                    continue;
                }
                col = qRgb(int(e * qRed(altcol) + f * qRed(col)),
                           int(e * qGreen(altcol) + f * qGreen(col)),
                           int(e * qBlue(altcol) + f * qBlue(col)));
            }
            pts[j] = col;
        }
        if (progdiag) {
            progdiag->setValue(i + 1);
            QCoreApplication::processEvents();
        }
    }
#if 1
    if (d.level_of_detail <= -1) {
        recursivelyPrintNCalls(d.root);
    }
#endif
#if 0
    qDebug("portion done after %d ms", t.msecsTo(QTime::currentTime()));
#endif
    emit completed();
    delete[] trackdists;
    delete[] trackcolors;
    return true;
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

static bool insideConvex(std::vector<Plane> clips, G4ThreeVector point) {
    // The empty convex hull by default contains everything
    for (const Plane &p : clips) {
        if (p.normal * point - p.offset < 0) {
            return false;
        }
    }
    return true;
}

static void lineCuts(std::vector<Plane> clips, G4ThreeVector from,
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
    // TODO: first, clipping plane filtering. A line can be cut at most
    // twice.
    // by the convex hull.
    // So, iterate over lines, trace entry and exit points.
    // Create a std::vector<> of points, and one of headers.
    // Finally, reconstruct the whole....

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
        if (insideConvex(vd.clipping_planes, fts) && intime && inenergy) {
            TrackHeader h;
            h.offset = cp.size();
            h.ptype = oheaders[i].ptype;
            h.npts = 1;
            ch.push_back(h);
            cp.push_back(seq[0]);
        }

        for (int j = 0; j < npts - 1; j++) {
            double low, high;
            G4ThreeVector pl(seq[j].x, seq[j].y, seq[j].z);
            G4ThreeVector ph(seq[j + 1].x, seq[j + 1].y, seq[j + 1].z);
            lineCuts(vd.clipping_planes, pl, ph, low, high);
            float tcutlow =
                (tlow - seq[j].time) / (seq[j + 1].time - seq[j].time);
            float tcuthigh =
                (thigh - seq[j].time) / (seq[j + 1].time - seq[j].time);
            low = std::max(double(tcutlow), low);
            high = std::min(double(tcuthigh), high);
            float ecutlow =
                (elow - seq[j].energy) / (seq[j + 1].energy - seq[j].energy);
            float ecuthigh =
                (ehigh - seq[j].energy) / (seq[j + 1].energy - seq[j].energy);
            low = std::max(double(ecutlow), low);
            high = std::min(double(ecuthigh), high);

            if (low <= 0. && high >= 1.) {
                // Keep point
                ch[ch.size() - 1].npts++;
                cp.push_back(seq[j + 1]);
            } else if (low <= 0. && high < 1.) {
                if (high > 0.) {
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
                    h.npts = 1;
                    ch.push_back(h);
                    cp.push_back(linmix(seq[j], seq[j + 1], low));
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
