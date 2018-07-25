/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#include "TrackData.hh"

#include <G4ThreeVector.hh>
#include <stddef.h>

typedef struct OctreeNode_s OctreeNode;
struct OctreeNode_s {
    // data?
    double density;
    OctreeNode *children; // pointer to an array of 8
};
typedef struct {
    G4ThreeVector min;
    G4ThreeVector max;
} Bounds;

typedef struct {
    size_t header_pos;
    int32_t i;
    float min, max;
} SegAddr;

typedef struct OctreeRoot_s {
    Bounds bounds;
    OctreeNode tree;
} OctreeRoot;

OctreeRoot *buildDensityOctree(const TrackBlock *blocks,
                               const QVector<SegAddr> &indices,
                               const Bounds &bounds);

typedef struct RayPoint_s RayPoint;
double traceDensityRay(const OctreeRoot &root, const G4ThreeVector &start,
                       const G4ThreeVector &direction, double distance);

void deleteOctree(const OctreeRoot *);
