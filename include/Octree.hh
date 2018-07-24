/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#include "TrackData.hh"

#include <G4ThreeVector.hh>
#include <stddef.h>

typedef struct OctreeNode_s OctreeNode;
struct OctreeNode_s {
    // data?
    double total;
    OctreeNode *children; // pointer to an array of 8
};
typedef struct {
    G4ThreeVector min;
    G4ThreeVector max;
} Bounds;

typedef struct {
    size_t track;
    int32_t i;
} SegAddr;

typedef struct OctreeRoot_s {
    Bounds bounds;
    OctreeNode tree;
} OctreeRoot;

OctreeRoot *buildDensityOctree(const TrackBlock *blocks,
                               const QVector<SegAddr> &indices,
                               const Bounds &bounds);

typedef struct RayPoint_s RayPoint;
FColor traceDensityRay(const OctreeRoot &root, const G4ThreeVector &start,
                       const G4ThreeVector &direction, double distance);

void deleteOctree(const OctreeRoot *);
