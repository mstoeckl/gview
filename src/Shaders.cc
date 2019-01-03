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
    cx = 1.0 - std::max(0.0, std::min(1.0, 0.9 * cx));

    if (is_clipping_plane) {
        const G4ThreeVector &orthA = normal.orthogonal().unit();
        const G4ThreeVector &orthB = normal.cross(orthA).unit();

        double aslp = orthA * position;
        double bslp = orthB * position;

        const double angle = -M_PI / 4 + -M_PI / 30;
        const double sin_r = std::sin(angle), cos_r = std::cos(angle);

        // project onto line
        double v = -aslp * sin_r + bslp * cos_r;

        double shade_factor = v * shade_scale;
        shade_factor *= 14;
        double fm = std::fmod(shade_factor, 1.);
        if (fm < 0.)
            fm += 1.;
        if (fm < 0.5)
            cx *= 0.75;
    }

    return FColor(cx * base.redF(), cx * base.greenF(), cx * base.blueF());
}

FColor defaultColorForRay(const RayPoint &ray, const FColor &trackcol,
                          G4double trackdist, double *voxel_cumulants,
                          const ViewData &d, const QPointF &pt,
                          const G4ThreeVector &forward, bool show_cut_marks) {
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
            return color;
        } else if (ecode == CODE_LINE) {
            // Lines are black
            color = FColor::add(color, FColor(0., 0., 0., 1.), weight);
            return color;
        } else {
            const Element &eback = d.elements[ecode];
            if (!eback.visible) {
                continue;
            }
            const VColor &base_color = d.color_table[eback.ccode];
            const G4ThreeVector &intpos = initPoint(pt, d) + ft.dist * forward;
            const FColor altcol = colorMap(
                ft, forward, base_color, intpos, 1. / d.scale,
                show_cut_marks && ((k == 0 && ray.front_clipped) ||
                                   (k == ray.N - 1 && ray.back_clipped)));

            float cur_fraction = (d.force_opaque ? 1. : eback.alpha);
            color = FColor::add(color, altcol, weight * cur_fraction);
            weight = weight * (1. - cur_fraction);
            if (weight <= 0.) {
                // Early termination
                return color;
            }
        }
    }
    // Finally, add background color
    return FColor::add(color, FColor(trackcol), weight);
}

FColor normalColorForRay(const RayPoint &ray, const FColor &trackcol,
                         G4double trackdist, double *voxel_cumulants,
                         const ViewData &d, const QPointF &,
                         const G4ThreeVector &, bool) {
    if (voxel_cumulants) {
        FColor voxc(0.5, 0.5, 0.5, 1.);
        FColor bgc(1., 1., 1., 0.);
        FColor line(0., 0., 0., 1.);

        // TODO: make last parameter & all color choices tunable
        double s = std::exp(-voxel_cumulants[0] / d.voxel_base_density);
        if (ray.N <= 0) {
            return FColor::blend(voxc, bgc, s);
        }
        const Intersection &i = ray.intersections[0];
        if (i.ecode == CODE_END) {
            return FColor::blend(voxc, bgc, s);
        }
        if (i.ecode == CODE_LINE) {
            return FColor::blend(voxc, line, s);
        }
        FColor col((i.normal.x + 1.0) / 2, (i.normal.y + 1.0) / 2,
                   (i.normal.z + 1.0) / 2);
        return FColor::blend(voxc, col, s);
    } else {
        if (ray.N <= 0) {
            return trackcol;
        }
        const Intersection &i = ray.intersections[0];
        if (i.dist > trackdist || i.ecode == CODE_END) {
            return trackcol;
        }
        if (i.ecode == CODE_LINE) {
            return FColor(0., 0., 0., 1.);
        }

        FColor col((i.normal.x + 1.0) / 2, (i.normal.y + 1.0) / 2,
                   (i.normal.z + 1.0) / 2);
        return col;
    }
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
