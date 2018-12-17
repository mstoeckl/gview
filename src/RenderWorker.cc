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
                ray.intersections[k].ecode = CODE_LINE;
                ray.N = k + 1;
                break;
            }
        }
    }

    return ray;
}

void debugRayPoint(const RayPoint &ray, const std::vector<Element> &els) {
    qDebug("> Raypoint: N=%d front_clipped=%c back_clipped=%c", ray.N,
           ray.front_clipped ? 'Y' : 'n', ray.back_clipped ? 'Y' : 'n');
    for (int i = 0; i < ray.N; i++) {
        int ecode = ray.intersections[i].ecode;
        qDebug("  %d: at %f mm, ecode %d (%s)", i,
               ray.intersections[i].dist / CLHEP::mm, ecode,
               ecode < 0 ? "END"
                         : els[ray.intersections[i].ecode].name.c_str());
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
                    const G4String &suffix, const G4RotationMatrix &parent_rot,
                    const G4ThreeVector &parent_offset, int *counter) {
    int cc = 0;
    if (!counter) {
        counter = &cc;
    }

    G4ThreeVector offset = phys->GetFrameTranslation();
    offset = parent_rot.inverse() * offset;
    const G4RotationMatrix &r = (phys->GetFrameRotation() != NULL)
                                    ? *phys->GetFrameRotation()
                                    : G4RotationMatrix();
    G4RotationMatrix rot = r * parent_rot;
    const G4LogicalVolume *log = phys->GetLogicalVolume();
    const G4Material *mat = log->GetMaterial();

    int k = elts.size();
    elts.push_back(Element());
    {
        // Reference only valid where elts is unmodified
        Element &m = elts[k];
        m.name =
            phys->GetName().substr(0, phys->GetName().size() - suffix.size());
        m.orig_vol = phys;

        // Only identity has a trace of +3 => norm2 of 0
        m.rotated = rot.norm2() > 1e-10;
        m.offset = offset;
        m.global_offset = offset + parent_offset, m.rot = rot;

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
        svi.push_back(convertCreation(elts, log->GetDaughter(i), suffix, rot,
                                      offset + parent_offset, counter));
    }
    ElemSort esort(elts);
    std::sort(svi.begin(), svi.end(), esort);

    elts[k].children = svi;

    return k;
}

/**
 * Given a view-port distance, and a corresponding real-space distance, return
 * a view-port distance less than or equal to the provided distance, along with
 * a compact label for the associated real-space length.
 */
QPair<double, QString> ruler_distance(double real_distance,
                                      double vp_distance) {
    const char *labels[7] = {"fm", "pm", "nm", "μm", "mm", "m", "km"};

    int s = std::floor(std::log10(real_distance / CLHEP::fermi));
    double dec_floor = std::pow(10., s) * CLHEP::fermi;
    int uclass = std::max(0, std::min(6, s / 3));
    double unit_distance = std::pow(1000., uclass) * CLHEP::fermi;

    double frac_part = real_distance / dec_floor;
    double rule_real_length;
    if (frac_part > 5) {
        rule_real_length = 5 * dec_floor;
    } else {
        if (frac_part < 1) {
            qWarning("Ruler fraction below unit: %f", frac_part);
        }
        rule_real_length = dec_floor;
    }
    double rule_vp_length = (rule_real_length / real_distance) * vp_distance;
    double rule_unit_length = rule_real_length / unit_distance;

    const QString &desc = QStringLiteral("%1 %2")
                              .arg(rule_unit_length, 0, 'g')
                              .arg(labels[uclass]);
    return QPair<double, QString>(rule_vp_length, desc);
}
