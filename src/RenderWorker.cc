#include "Viewer.hh"

#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VSolid.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VisExtent.hh>

#include <QCoreApplication>
#include <QProgressDialog>

typedef struct {
    G4double dist;
    G4ThreeVector normal;
} Intersection;

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

static int traceRay(const QPointF &scpt, const ViewData &d,
                    const Element *hits[], Intersection ints[], int maxhits) {
    // Given a square camera, p in [0,1]X[0,1]. Rectangles increase one dim.
    // Return number of hits recorded.

    G4ThreeVector forward = d.orientation.rowX();
    G4ThreeVector updown = d.orientation.rowY();
    G4ThreeVector lright = d.orientation.rowZ();

    // Q: proper camera dynamics...
    G4ThreeVector init =
        d.camera + updown * scpt.y() * d.scale + lright * scpt.x() * d.scale;

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

    // Assume we are already in the root element
    std::vector<const Element *> stack;
    stack.push_back(&d.root);
    G4ThreeVector local = start;
    for (int iter = 0; iter < 1000; iter++) {
        const Element &last = *stack[stack.size() - 1];

        // Strictly inclusion tests....
        // check inside on all children, append; stop when no longer dropping
        // levels
        bool found = false;
        for (const Element &elem : last.children) {
            G4ThreeVector sub;
            if (elem.rotated) {
                sub = elem.rot * (local + elem.offset);
            } else {
                sub = local + elem.offset;
            }
            // rotate/offset start point...
            ++elem.ngeocalls;
            if (elem.solid->Inside(sub)) {
                stack.push_back(&elem);
                found = true;
                local = sub;
                break;
            }
        }
        if (!found) {
            break;
        }
    }
    // Create direction vector sequence....
    // (we can tack on normals later; for now, defnormal works...)
    std::vector<G4ThreeVector> directions;
    directions.push_back(forward);
    std::vector<G4ThreeVector> positions;
    positions.push_back(start);
    for (size_t i = 1; i < stack.size(); i++) {
        const Element &elem = *stack[i];
        if (elem.rotated) {
            directions.push_back(elem.rot * directions[i - 1]);
            positions.push_back(elem.rot * (positions[i - 1] + elem.offset));
        } else {
            directions.push_back(directions[i - 1]);
            positions.push_back(positions[i - 1] + elem.offset);
        }
    }

    int nhits = 1;
    if (clippable) {
        // Start point is in the world volume, and is an intersection
        hits[0] = stack[stack.size() - 1];
        ints[0] = {sdist, entrynormal};
    } else {
        // Starting point
        ++d.root.ngeocalls;
        G4ThreeVector lnormal = d.root.solid->SurfaceNormal(local);
        hits[0] = &d.root;
        ints[0] = {sdist, lnormal};
        // Intersection on entry
    }

    for (int iter = 0; iter < 1000; iter++) {
        const Element &curr = *stack[stack.size() - 1];
        const G4ThreeVector &dir = directions[stack.size() - 1];
        const G4ThreeVector &pos = positions[stack.size() - 1];

        ++curr.ngeocalls;
        G4double exitdist = curr.solid->DistanceToOut(pos, dir);

        G4double closestDist = 2 * kInfinity;
        const Element *closest = NULL;
        G4ThreeVector closestDir;
        G4ThreeVector closestPos;
        for (const Element &elem : curr.children) {
            G4ThreeVector sub, subdir;
            if (elem.rotated) {
                sub = elem.rot * (pos + elem.offset);
                subdir = elem.rot * dir;
            } else {
                sub = pos + elem.offset;
                subdir = dir;
            }
            ++elem.ngeocalls;
            G4double altdist = elem.solid->DistanceToIn(sub, subdir);
            if (altdist >= kInfinity) {
                // No intersection
                continue;
            } else if (altdist < closestDist) {
                closestDir = subdir;
                closestPos = sub;
                closest = &elem;
                closestDist = altdist;
            }
        }
        G4double fdist = std::min(exitdist, closestDist);
        if (fdist > edist - sdist) {
            // End clip (counts as a hit!)
            ints[nhits] = {edist, exitnormal};
            return nhits;
        } else if (exitdist < closestDist) {
            // Transition on leaving vol w/o intersections
            // Advance all levels, why not
            for (size_t j = 0; j < positions.size(); j++) {
                positions[j] += exitdist * directions[j];
            }
            sdist += exitdist;
            G4ThreeVector lpos = positions[positions.size() - 1];
            stack.pop_back();
            directions.pop_back();
            positions.pop_back();

            // Only calc normal if need be
            ++curr.ngeocalls;
            G4ThreeVector lnormal = curr.solid->SurfaceNormal(lpos);
            ints[nhits] = {sdist, lnormal};
            if (stack.size() == 0 || nhits >= maxhits) {
                return nhits;
            }
            hits[nhits] = stack[stack.size() - 1];
            nhits++;
        } else {
            // Transition on visiting a child
            positions.push_back(closestPos);
            directions.push_back(closestDir);
            for (size_t j = 0; j < positions.size(); j++) {
                positions[j] += closestDist * directions[j];
            }
            stack.push_back(closest);
            sdist += closestDist;

            ++closest->ngeocalls;
            G4ThreeVector lnormal =
                closest->solid->SurfaceNormal(positions[positions.size() - 1]);
            ints[nhits] = {sdist, lnormal};
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

static int compressTraces(const Element *hits[], Intersection ints[], int m) {
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
    // Next, nuke invisible slices
    // TODO ALERT FIXME: otherwise we get weird compression bugs
    // (The alpha=0 fails because same-material invisibility creates holes
    // (really ought to do depth-absorption routines instead)

    // Finally, nuke intersections between similar material'd objects
    int p = 0;
    // 001122334  =>  0022334
    // |x|x|y|x|  =>  |x|y|x|
    for (int i = 0; i < n; i++) {
        if (i == 0 || hits[i]->mat != hits[i - 1]->mat) {
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

RenderWorker::RenderWorker() { abort_task = false; }

RenderWorker::~RenderWorker() {}

void RenderWorker::coAbort() { abort_task = true; }
void RenderWorker::flushAbort() { abort_task = false; }

bool RenderWorker::render(ViewData d, QImage *next, int slice, int nslices,
                          QProgressDialog *progdiag) {
    if (!d.root.solid) {
        return false;
    }
    int w = next->width();
    int h = next->height();
#if 0
    QTime t = QTime::currentTime();
    qDebug("render lod %d w %d h %d | %d of %d", d.level_of_detail, w, h, slice,
           nslices);
#endif

    int mind = w > h ? h : w;
    int hl = (slice * h / nslices);
    int hh = slice == nslices - 1 ? h : ((slice + 1) * h / nslices);
    const G4double radius = 0.8 / mind;
    for (int i = hl; i < hh; i++) {
        QRgb *pts = reinterpret_cast<QRgb *>(next->scanLine(i));
        for (int j = 0; j < w; j++) {
            if (this->abort_task) {
                this->abort_task = false;
                emit aborted();
                return false;
            }
            QPointF pt((j - w / 2.) / (2. * mind), (i - h / 2.) / (2. * mind));

            const int M = 30;
            const Element *hits[M];
            Intersection ints[M + 1];
            int m = traceRay(pt, d, hits, ints, M);
            m = compressTraces(hits, ints, m);

            // TODO: List editing on hits/dists/normals
            // to eliminate dummy transitions between objects of the
            // same material. (if three things in a row have the same mtl,
            // cut the middle one; do the same for the level_of_detail.
            // (Note: is a useful precursor step for the depth-based approach)
            // --> eliminates spurious lines of all sorts..
            int p = m - 1;
            bool line = false;
            if (d.level_of_detail <= -1) {
                int ndevs[M + 1];
                for (int k = 0; k < M + 1; k++) {
                    ndevs[k] = 0;
                }
                G4double seed = CLHEP::pi / 5 * (qrand() >> 16) / 16384.0;
                const Element *althits[M];
                Intersection altints[M + 1];
                int am;
                for (int k = 0; k < 10; k++) {

                    QPointF off(radius * std::cos(seed + k * CLHEP::pi / 5),
                                radius * std::sin(seed + k * CLHEP::pi / 5));
                    am = traceRay(pt + off, d, althits, altints, M);
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

            QRgb col = line ? qRgb(0, 0, 0) : qRgb(255, 255, 255);
            // p indicates the first volume to use the color rule on
            // (p<0 indicates the line dominates)
            for (int k = p; k >= 0; --k) {
                // We use the intersection before the volume
                QRgb altcol = colorMap(ints[k].normal, d.orientation.rowX(),
                                       hits[k]->hue, ints[k].dist);
                double e = hits[k]->alpha, f = (1. - hits[k]->alpha);
                if (!hits[k]->visible) {
                    e = 0.;
                    f = 1.;
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
    return true;
}
