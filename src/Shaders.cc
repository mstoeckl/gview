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
