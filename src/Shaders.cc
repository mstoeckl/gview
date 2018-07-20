/* SPDX-License-Identifier: GPL-3.0-only */
#include "Shaders.hh"

#include "Navigator.hh"

#include <QColor>

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

QRgb defaultColorForRay(const RayPoint &ray, QRgb trackcol, G4double trackdist,
                        const ViewData &d, const QPointF &pt,
                        const G4ThreeVector &forward) {
    /* Scan from front to back, as this gives the option of quitting early */
    /* Each intersection is a transition between two domains */
    FColor color(0., 0., 0., 0.);
    float weight = 1.0; // of the current component

    for (int k = 0; k < ray.N; k++) {
        const Intersection &ft = ray.intersections[k];
        int ecode = ft.ecode;
        if (trackdist < ft.dist || ecode == CODE_END) {
            // Early termination
            color = FColor::add(color, FColor(trackcol), weight);
            return color.rgba();
        } else if (ecode == CODE_LINE) {
            // Lines are black
            color = FColor::add(color, FColor(0., 0., 0., 1.), weight);
            return color.rgba();
        } else {
            const Element &eback = d.elements[ecode];
            if (!eback.visible) {
                continue;
            }
            const VColor &base_color = d.color_table[eback.ccode];
            const G4ThreeVector &intpos = initPoint(pt, d) + ft.dist * forward;
            const FColor altcol =
                colorMap(ft, forward, base_color, intpos, 1. / d.scale,
                         (k == 0 && ray.front_clipped) ||
                             (k == ray.N - 1 && ray.back_clipped));

            float cur_fraction = (d.force_opaque ? 1. : eback.alpha);
            color = FColor::add(color, altcol, weight * cur_fraction);
            weight = weight * (1. - cur_fraction);
            if (weight <= 0.) {
                // Early termination
                return color.rgba();
            }
        }
    }
    // Finally, add background color
    color = FColor::add(color, FColor(trackcol), weight);
    return color.rgba();
}

QRgb normalColorForRay(const RayPoint &ray, QRgb trackcol, G4double trackdist,
                       const ViewData &, const QPointF &,
                       const G4ThreeVector &) {
    if (ray.N <= 0) {
        return trackcol;
    }

    const Intersection &i = ray.intersections[0];
    if (i.dist > trackdist) {
        return trackcol;
    }

    FColor col((i.normal.x + 1.0) / 2, (i.normal.y + 1.0) / 2,
               (i.normal.z + 1.0) / 2);
    return col.rgba();
}

void rainbowColorForSegment(const TrackHeader &h, const TrackPoint &pa,
                            const TrackPoint &pb, const G4ThreeVector &a,
                            const G4ThreeVector &b,
                            const G4ThreeVector &forward, FColor &ca,
                            FColor &cb, float &wa, float &wb) {
    G4ThreeVector normal = b - a;
    double coa = (normal * forward) / normal.mag();

    //    double mtime = double(pa.time) * (1. - interp) * double(pb.time) *
    //    double ptime = (int(mtime / CLHEP::ns * 1024.) % 1024) / 1024.;

    double aloge = std::log10(double(pa.energy)) / 2.0;
    double patime = aloge - std::floor(aloge);
    ca = FColor(QColor::fromHsvF(patime, 0.8, 1.0 - 0.5 * std::abs(coa)).rgb());

    double bloge = std::log10(double(pb.energy)) / 2.0;
    double pbtime = bloge - std::floor(bloge);
    cb = FColor(QColor::fromHsvF(pbtime, 0.8, 1.0 - 0.5 * std::abs(coa)).rgb());

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

void typeColorForSegment(const TrackHeader &h, const TrackPoint &,
                         const TrackPoint &, const G4ThreeVector &a,
                         const G4ThreeVector &b, const G4ThreeVector &forward,
                         FColor &ca, FColor &cb, float &wa, float &wb) {
    G4ThreeVector normal = b - a;
    double coa = std::abs((normal * forward) / normal.mag());
    coa = 1 - coa / 2;

    FColor thue(0.2, 0.2, 0.2);
    if (h.ptype == 11) {
        thue = FColor(0.3, 1.0, 0.0);
    } else if (h.ptype == 22) {
        thue = FColor(1.0, 0.0, 1.0);
    } else if (h.ptype == 0) {
        thue = FColor(0., 0., 1.);
    }

    cb = ca =
        FColor(thue.redF() * coa, thue.greenF() * coa, thue.blueF() * coa);

    wa = 0.7f;
    wb = 0.7f;
}
