/* SPDX-License-Identifier: GPL-3.0-only */
#include "RenderWorker.hh"
#include "BooleanTree.hh"
#include "Navigator.hh"

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
#include <QSet>

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

bool clipRay(const std::vector<Plane> &clipping_planes,
             const G4ThreeVector &init, const G4ThreeVector &forward,
             G4double &sdist, G4double &edist, G4ThreeVector &entrynormal,
             G4ThreeVector &exitnormal) {
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

RayPoint rayAtPoint(Navigator &nav, const QPointF &pt, qreal radius,
                    const G4ThreeVector &forward, const ViewData &d,
                    Intersection *ints, Intersection *altints, int M,
                    int *ndevs) {
    RayPoint ray =
        nav.traceRay(initPoint(pt, d), forward, ints, M, d.force_opaque);

    if (d.level_of_detail <= -1) {
        for (int k = 0; k < ray.N; k++) {
            ndevs[k] = 0;
        }
        G4double seed = CLHEP::pi / 5 * (qrand() >> 16) / 16384.0;
        for (int k = 0; k < 10; k++) {
            QPointF off(radius * std::cos(seed + k * CLHEP::pi / 5),
                        radius * std::sin(seed + k * CLHEP::pi / 5));
            RayPoint aray = nav.traceRay(initPoint(pt + off, d), forward,
                                         altints, M, d.force_opaque);
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
                    ray.intersections[k - 1].ecode = CODE_LINE;
                }
                break;
            }
        }
    }

    return ray;
}

void debugRayPoint(const RayPoint &ray) {
    qDebug("> Raypoint: N=%d front_clipped=%c back_clipped=%c", ray.N,
           ray.front_clipped ? 'Y' : 'n', ray.back_clipped ? 'Y' : 'n');
    for (int i = 0; i < ray.N; i++) {
        qDebug("  %d: ecode %d at %f cm", i, ray.intersections[i].ecode,
               ray.intersections[i].dist / CLHEP::cm);
    }
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
    bool line = ray.N > 0 && (ray.intersections[ray.N - 1].ecode == CODE_LINE);
    if (line && p == ray.N - 1) {
        p = ray.N - 2;
    }
    FColor col =
        (line && !rayoverride) ? FColor(0., 0., 0., 1.f) : FColor(trackcol);
    // p indicates the first volume to use the color rule on
    // (p<0 indicates the line dominates)
    for (int k = p; k >= 0; --k) {
        if (ray.intersections[k].ecode < 0) {
            // Either line or END_CODE
            col = FColor(1.0, 0., 0., 0);
            continue;
        }

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
        m.orig_vol = phys;

        // Only identity has a trace of +3 => norm2 of 0
        m.rotated = rot.norm2() > 1e-10;
        m.offset = offset;
        m.rot = rot;

        m.ccode = 0;
        m.material = mat;
        m.solid = log->GetSolid();
        m.cubicVolume = -std::numeric_limits<G4double>::infinity();
        m.surfaceArea = -std::numeric_limits<G4double>::infinity();
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
