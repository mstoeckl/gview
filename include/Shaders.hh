/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#include "RenderWorker.hh"

// Q: just use a function pointer -- could save time & avoid allocproblem

typedef FColor GeoShader(const RayPoint &ray, const FColor &trackcol,
                         G4double trackdist, double *voxel_density,
                         const ViewData &d, const QPointF &pt,
                         const G4ThreeVector &forward, bool show_cut_marks);
typedef void TrackShader(const TrackHeader &h, const TrackPoint &pa,
                         const TrackPoint &pb, const G4ThreeVector &a,
                         const G4ThreeVector &b, const G4ThreeVector &forward,
                         FColor &ca, FColor &cb, float &wa, float &wb);

enum {
    gshaderDefault,
    gshaderNormal,
};
enum {
    tshaderRainbow,
    tshaderType,
};

GeoShader defaultColorForRay;
GeoShader normalColorForRay;

TrackShader rainbowColorForSegment;
TrackShader typeColorForSegment;

/* Placed in header for faster selection */
inline TrackShader *getTrackShader(int type) {
    switch (type) {
    case tshaderRainbow:
        return rainbowColorForSegment;
    case tshaderType:
        return typeColorForSegment;
    default:
        return NULL;
    }
}
inline GeoShader *getGeoShader(int type) {
    switch (type) {
    case gshaderDefault:
        return defaultColorForRay;
    case gshaderNormal:
        return normalColorForRay;
    default:
        return NULL;
    }
}
