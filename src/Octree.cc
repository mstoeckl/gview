#/* SPDX-License-Identifier: GPL-3.0-only */
#include "Octree.hh"

#include <CLHEP/Units/SystemOfUnits.h>
#include <G4ThreeVector.hh>
#include <geomdefs.hh>

#define DEPTH_LIMIT 10

static G4ThreeVector centerPoint(const Bounds &b) {
    return G4ThreeVector(0.5 * (b.min[0] + b.max[0]),
                         0.5 * (b.min[1] + b.max[1]),
                         0.5 * (b.min[2] + b.max[2]));
}

static bool lineCrossesBlock(const G4ThreeVector &a, const G4ThreeVector &b,
                             const Bounds &z) {
    const G4ThreeVector d = b - a;
    // At least one axis is nonzero ; NaNs should cancel away
    const G4ThreeVector id(1. / d.x(), 1. / d.y(), 1. / d.z());
    double tmin = -kInfinity, tmax = kInfinity;

    for (int i = 0; i < 3; ++i) {
        double t1 = (z.min[i] - a[i]) * id[i];
        double t2 = (z.max[i] - b[i]) * id[i];

        tmin = std::max(tmin, std::min(t1, t2));
        tmax = std::min(tmax, std::max(t1, t2));
    }
    return tmax > std::max(tmin, 0.0);
}
static bool pointInBlock(const G4ThreeVector &a, const Bounds &bounds) {
    bool above =
        a[0] >= bounds.min[0] && a[1] >= bounds.min[1] && a[2] >= bounds.min[2];
    bool below =
        a[1] <= bounds.max[0] && a[1] <= bounds.max[1] && a[2] <= bounds.max[2];
    return above && below;
}

// static bool trackCrossesBlock(const TrackPoint *pts, int N,
//                              const Bounds &bounds) {
//    if (pointInBlock(G4ThreeVector(pts[0].x, pts[0].y, pts[0].z), bounds)) {
//        return true;
//    }
//    if (pointInBlock(G4ThreeVector(pts[N - 1].x, pts[N - 1].y, pts[N - 1].z),
//                     bounds)) {
//        return true;
//    }
//    for (int i = 1; i < N; i++) {
//        G4ThreeVector a(pts[i - 1].x, pts[i - 1].y, pts[i - 1].z);
//        G4ThreeVector b(pts[i].x, pts[i].y, pts[i].z);
//        if (lineCrossesBlock(a, b, bounds)) {
//            return true;
//        }
//    }
//    return false;
//}

static QVector<SegAddr> selectTracksInBlock(const TrackBlock *blocks,
                                            const QVector<SegAddr> &indices,
                                            const Bounds &subbounds) {
    QVector<SegAddr> ret;
    for (SegAddr addr : indices) {
        const TrackPoint *seq = &(blocks[addr.header_pos + 1].p);
        const G4ThreeVector va(seq[addr.i].x, seq[addr.i].y, seq[addr.i].z);
        const G4ThreeVector vb(seq[addr.i + 1].x, seq[addr.i + 1].y,
                               seq[addr.i + 1].z);

        // TODO: replace with a 'segment crosses block' primitive
        // Take note of repeated work
        if (pointInBlock(va, subbounds) || pointInBlock(vb, subbounds) ||
            lineCrossesBlock(va, vb, subbounds)) {
            ret.push_back(addr);
        }
    }
    return ret;
}

static double blockVol(const Bounds &bounds) {
    return (bounds.max[2] - bounds.min[2]) * (bounds.max[1] - bounds.min[1]) *
           (bounds.max[0] - bounds.min[0]);
}

typedef struct {
    bool below_is_A, below_is_B;
} PtSides;

static void initOctree(OctreeNode &node, const TrackBlock *blocks,
                       SegAddr **seg_buffer, int N_at_depth,
                       const Bounds &bounds, int split_thresh, int depth,
                       long &node_count, int nevt) {
    if (node_count % 1000 == 0) {
        qDebug("node %ld: depth=%d n=%d", node_count, depth, N_at_depth);
    }
    node_count++;
    // Calculate node density per unit. As lines typically span length, multiply
    // by estimated net line amount, around (V^(1/3)).
    node.density =
        N_at_depth / (double)nevt / std::pow(blockVol(bounds), 1. / 1.5);

    if (N_at_depth < split_thresh || depth >= DEPTH_LIMIT - 1) {
        node.children = NULL;
        return;
    } else {
        node.children = new OctreeNode[8];
        G4ThreeVector mid = centerPoint(bounds);
        bool by_planes = true;
        QVector<PtSides> sides[3];
        QVector<float> splits[3];
        SegAddr *seg_parent = seg_buffer[depth];
        SegAddr *seg_child = seg_buffer[depth + 1];
        if (by_planes) {
            for (int axis = 0; axis < 3; axis++) {
                sides[axis].resize(N_at_depth);
                splits[axis].resize(N_at_depth);

                // Along each axis, identify the line splitting point (if there
                // is one) and determine where the parts belong
                for (int k = 0; k < N_at_depth; k++) {
                    const SegAddr &addr = seg_parent[k];
                    const TrackBlock *seq = &blocks[addr.header_pos + 1];
                    const G4ThreeVector va(seq[addr.i].p.x, seq[addr.i].p.y,
                                           seq[addr.i].p.z);
                    const G4ThreeVector vb(seq[addr.i + 1].p.x,
                                           seq[addr.i + 1].p.y,
                                           seq[addr.i + 1].p.z);

                    G4ThreeVector pmin = va * (1. - addr.min) + addr.min * vb;
                    G4ThreeVector pmax = va * (1. - addr.max) + addr.max * vb;
                    PtSides sch = {pmin[axis] < mid[axis],
                                   pmax[axis] < mid[axis]};
                    if (sch.below_is_A != sch.below_is_B) {
                        // The segment spans the middle; we calculate the split
                        // point
                        double t =
                            (mid[axis] - va[axis]) / (vb[axis] - va[axis]);
                        splits[axis][k] = t;
                    }
                    sides[axis][k] = sch;
                }
            }
        }

        for (uint32_t i = 0; i < 8; i++) {
            bool upper[3] = {(bool)(i & 0x01), (bool)(i & 0x02),
                             (bool)(i & 0x04)};
            Bounds subbounds;
            for (int k = 0; k < 3; k++) {
                subbounds.min[k] = upper[k] ? mid[k] : bounds.min[k];
                subbounds.max[k] = upper[k] ? bounds.max[k] : mid[k];
            }

            // Based on the side classification, select all matching

            int N_sub = 0;
            if (by_planes) {
                for (int k = 0; k < N_at_depth; k++) {
                    SegAddr seg = seg_parent[k];
                    bool keep = true;
                    for (int axis = 0; axis < 3; axis++) {
                        PtSides sch = sides[axis][k];
                        if (upper[axis]) {
                            if (sch.below_is_A && sch.below_is_B) {
                                keep = false;
                                break;
                            } else if (!sch.below_is_A && !sch.below_is_B) {
                                // segment above, no action
                            } else if (sch.below_is_A && !sch.below_is_B) {
                                // segment moving up and in
                                seg.min = std::max(seg.min, splits[axis][k]);
                            } else if (!sch.below_is_A && sch.below_is_B) {
                                // segment moving down and out
                                seg.max = std::min(seg.max, splits[axis][k]);
                            }
                        } else {
                            if (sch.below_is_A && sch.below_is_B) {
                                // segment below, no action
                            } else if (!sch.below_is_A && !sch.below_is_B) {
                                keep = false;
                                break;
                            } else if (sch.below_is_A && !sch.below_is_B) {
                                // segment moving up and out
                                seg.max = std::min(seg.max, splits[axis][k]);
                            } else if (!sch.below_is_A && sch.below_is_B) {
                                // segment moving down and in
                                seg.min = std::max(seg.min, splits[axis][k]);
                            }
                        }
                    }
                    if (keep && seg.max > seg.min) {
                        seg_child[N_sub++] = seg;
                    }
                }
            } else {
                // subindices = selectTracksInBlock(blocks, indices, subbounds);
            }
            initOctree(node.children[i], blocks, seg_buffer, N_sub, subbounds,
                       split_thresh, depth + 1, node_count, nevt);
        }
    }
}

OctreeRoot *buildDensityOctree(const TrackBlock *blocks,
                               const QVector<SegAddr> &indices,
                               const Bounds &bounds) {
    OctreeRoot *root = new OctreeRoot();
    root->bounds = bounds;

    // Buffer setup (TODO: limit overallocation)
    SegAddr **temp_stacks = new SegAddr *[DEPTH_LIMIT + 1];
    temp_stacks[0] = (SegAddr *)indices.data();
    for (int i = 1; i < DEPTH_LIMIT + 1; i++) {
        temp_stacks[i] = new SegAddr[indices.size()];
    }

    long node_count = 0;
    int split_thresh = indices.size() / 10000;
    initOctree(root->tree, blocks, temp_stacks, indices.size(), root->bounds,
               split_thresh, 0, node_count, indices.size());
    long max_node_count =
        ((1L << (DEPTH_LIMIT * 3 + 3)) - 1) / ((1L << DEPTH_LIMIT) - 1);

    // Buffer cleanup
    for (int i = 1; i < DEPTH_LIMIT + 1; i++) {
        delete[] temp_stacks[i];
    }
    delete[] temp_stacks;

    qDebug("%ld/%ld total octree nodes", node_count, max_node_count);
    return root;
}
static void delOct(const OctreeNode &n) {
    if (n.children) {
        for (int i = 0; i < 8; i++) {
            delOct(n.children[i]);
        }
        delete[] n.children;
    }
}
void deleteOctree(const OctreeRoot *n) {
    if (n) {
        delOct(n->tree);
        delete n;
    }
}

static double crossDistance(const G4ThreeVector &p, const G4ThreeVector &dir,
                            const Bounds &b, double t_start, double t_stop) {
    for (int i = 0; i < 3; i++) {
        if (dir[i] == 0.) {
            if (p[i] < b.min[i] || p[i] > b.max[i]) {
                return 0.;
            }
        } else {
            double id = 1. / dir[i];
            double t_min = (b.min[i] - p[i]) * id;
            double t_max = (b.max[i] - p[i]) * id;
            double t_low = std::min(t_min, t_max);
            double t_high = std::max(t_min, t_max);
            t_start = std::max(t_low, t_start);
            t_stop = std::min(t_high, t_stop);
        }
    }
    if (t_stop < t_start) {
        return 0.;
    }
    return t_stop - t_start;
}

static void clipBox(const G4ThreeVector &p, const G4ThreeVector &dir,
                    const Bounds &b, double &t_start, double &t_stop,
                    int &fclipaxis, int &bclipaxis) {
    double t_low[3], t_high[3];
    for (int i = 0; i < 3; i++) {
        if (dir[i] == 0.) {
            // parallel to face so no clip
            // -- TODO: clip by source point
            if (p[i] < b.min[i] || p[i] > b.max[i]) {
                t_low[i] = INFINITY;
                t_high[i] = -INFINITY;
            } else {
                t_low[i] = -INFINITY;
                t_high[i] = INFINITY;
            }
        } else {
            double id = 1. / dir[i];
            double t_min = (b.min[i] - p[i]) * id;
            double t_max = (b.max[i] - p[i]) * id;
            t_low[i] = std::min(t_min, t_max),
            t_high[i] = std::max(t_min, t_max);
        }
    }

    // Intersect [t_start, t_stop] with [t_low[i], t_high[i]]
    for (int i = 0; i < 3; i++) {
        if (t_low[i] > t_start) {
            t_start = t_low[i];
            fclipaxis = i;
        }
        if (t_high[i] < t_stop) {
            t_stop = t_high[i];
            bclipaxis = i;
        }
    }
}

static double recursiveCubeIntegral(const OctreeNode &node,
                                    const Bounds &bounds,
                                    const G4ThreeVector &start,
                                    const G4ThreeVector &direction,
                                    double t_max) {
    double cxdist = crossDistance(start, direction, bounds, 0, t_max);
    if (cxdist <= 0.) {
        return 0.;
    }

    if (node.children) {
        double total = 0.;
        G4ThreeVector mid = centerPoint(bounds);
        for (uint32_t i = 0; i < 8; i++) {
            bool l[3] = {(bool)(i & 0x01), (bool)(i & 0x02), (bool)(i & 0x04)};
            Bounds subbounds;
            for (int k = 0; k < 3; k++) {
                subbounds.min[k] = l[k] ? mid[k] : bounds.min[k];
                subbounds.max[k] = l[k] ? bounds.max[k] : mid[k];
            }
            total += recursiveCubeIntegral(node.children[i], subbounds, start,
                                           direction, t_max);
        }
        return total;

    } else {
        return node.density * cxdist;
    }
}

double traceDensityRay(const OctreeRoot &root, const G4ThreeVector &start,
                       const G4ThreeVector &direction, double t_max) {
    // TODO: clip start/end offsets, and enter true recursive case

    // TODO: locate start location, & use recursive navigation method
    // (but with an *explicit* stack)
    //    int coords[DEPTH_LIMIT + 2];
    //    int depth = 0;

    // TODO: clip to bounds!
    double t_start = -kInfinity, t_end = kInfinity;
    int fclipaxis = -1, bclipaxis = -1;
    clipBox(start, direction, root.bounds, t_start, t_end, fclipaxis,
            bclipaxis);
    t_end = std::min(t_max, t_end);

    //    qDebug("%f %f | %d %d", t_start, t_end, fclipaxis, bclipaxis);

    if (t_start < 0) {
        return 0.; // FColor(0., 0, 0, 1);
    }
    if (t_end < t_start) {
        return 0.; // FColor(0.5, 0.5, 0.5, 1);
    }
    G4ThreeVector contact = start + direction * t_start;

    // recursive cube sum.

    double v = recursiveCubeIntegral(root.tree, root.bounds,
                                     start + t_start * direction, direction,
                                     t_end - t_start);

    // Normalize by the total tree density
    double xd = (root.bounds.max - root.bounds.min).mag();
    double scale_val = root.tree.density * xd;
    return v / scale_val;
}
