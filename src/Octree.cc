#/* SPDX-License-Identifier: GPL-3.0-only */
#include "Octree.hh"

#include <CLHEP/Units/SystemOfUnits.h>
#include <G4ThreeVector.hh>
#include <geomdefs.hh>

#define DEPTH_LIMIT 7

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
        const TrackPoint *seq = &(blocks[addr.track + 1].p);
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

static void initOctree(OctreeNode &node, const TrackBlock *blocks,
                       const QVector<SegAddr> &indices, const Bounds &bounds,
                       int split_thresh, int depth, long &node_count) {
    qDebug("node %ld: depth=%d n=%d", node_count, depth, indices.size());
    node_count++;
    // relative density (?)
    node.total = indices.size();
    if (indices.size() < split_thresh || depth >= DEPTH_LIMIT - 1) {
        node.children = NULL;
        return;
    } else {
        node.children = new OctreeNode[8];
        G4ThreeVector mid = centerPoint(bounds);
        for (uint32_t i = 0; i < 8; i++) {
            bool l[3] = {(bool)(i & 0x01), (bool)(i & 0x02), (bool)(i & 0x04)};
            Bounds subbounds;
            for (int k = 0; k < 3; k++) {
                subbounds.min[k] = l[k] ? mid[k] : bounds.min[k];
                subbounds.max[k] = l[k] ? bounds.max[k] : mid[k];
            }
            // TODO: classify the paths along each axis -- which halfbox?
            // after all, 6 halfboxes is easier than 8 corners.
            // also, perform class using the much shorter segments.

            // furthermore -- we *already* know that the line segment intersects
            // the box. So really, *if* we store tstart/tend, then
            // plane comparisons (against the composite vectors) are all that
            // are necessary; (these in turn, generate tmid-planeX, tmid-planeY,
            // etc.) then moving down, a simple sorting network resolves the
            // final T range.

            // TODO: how do we handle tight bundles of paths?
            const QVector<SegAddr> &subindices =
                selectTracksInBlock(blocks, indices, subbounds);
            initOctree(node.children[i], blocks, subindices, subbounds,
                       split_thresh, depth + 1, node_count);
        }
    }
}

OctreeRoot *buildDensityOctree(const TrackBlock *blocks,
                               const QVector<SegAddr> &indices,
                               const Bounds &bounds) {
    OctreeRoot *root = new OctreeRoot();
    root->bounds = bounds;

    long node_count = 0;
    int split_thresh = indices.size() / 1000;
    initOctree(root->tree, blocks, indices, root->bounds, split_thresh, 0,
               node_count);
    long max_node_count =
        ((1L << (DEPTH_LIMIT * 3 + 3)) - 1) / ((1L << DEPTH_LIMIT) - 1);

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
    delOct(n->tree);
    delete n;
}

static float hyptan(float f) {
    return (std::exp(f) - std::exp(-f)) / (std::exp(f) + std::exp(-f));
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
        //        double t_start = 0., t_end = t_max;
        //        int fclipaxis = -1, bclipaxis = -1;
        //        clipBox(start, direction, bounds, t_start, t_end, fclipaxis,
        //        bclipaxis); if (t_end <= t_start) {
        //            return 0;
        //        }
        //        return (t_end - t_start) * node.value;

        // sums can be greater than individual?
        return node.total;
        //        return cxdist * node.total;
    } /*


     // For now, trivial iteration.
     double total = 0.;
     G4ThreeVector mid = centerPoint(bounds);
     for (uint32_t i = 0; i < 8; i++) {
         bool l[3] = {(bool)(i & 0x01), (bool)(i & 0x02), (bool)(i & 0x04)};
         Bounds subbounds;
         for (int k = 0; k < 3; k++) {
             subbounds.min[k] = l[k] ? mid[k] : bounds.min[k];
             subbounds.max[k] = l[k] ? bounds.max[k] : mid[k];
         }

         double t_start = 0., t_end = t_max;
         int fclipaxis = -1, bclipaxis = -1;
         clipBox(start, direction, subbounds, t_start, t_end, fclipaxis,
                 bclipaxis);
         if (t_end <= t_start) {
             continue;
         }

         double amt = t_end - t_start;
         total += amt;
     }
     qDebug("%f %f", t_max, total);
     return total;*/
}

FColor traceDensityRay(const OctreeRoot &root, const G4ThreeVector &start,
                       const G4ThreeVector &direction, double t_max) {
    // TODO: clip start/end offsets, and enter true recursive case

    // TODO: locate start location, & use recursive navigation method
    // (but with an *explicit* stack)
    int coords[DEPTH_LIMIT + 2];
    int depth = 0;

    // TODO: clip to bounds!
    double t_start = -kInfinity, t_end = kInfinity;
    int fclipaxis = -1, bclipaxis = -1;
    clipBox(start, direction, root.bounds, t_start, t_end, fclipaxis,
            bclipaxis);
    t_end = std::min(t_max, t_end);

    //    qDebug("%f %f | %d %d", t_start, t_end, fclipaxis, bclipaxis);

    if (t_start < 0) {
        return FColor(0., 0, 0, 1);
    }
    if (t_end < t_start) {
        return FColor(0.5, 0.5, 0.5, 1);
    }
    G4ThreeVector contact = start + direction * t_start;

    // recursive cube sum.

    double v = recursiveCubeIntegral(root.tree, root.bounds,
                                     start + t_start * direction, direction,
                                     t_end - t_start);

    //    qDebug("%f %f %f mm", contact.x(), contact.y(), contact.z());
    //    if (t_start > t_end) {
    //        return FColor(1., 0, 0, 1);
    //    }

    //    double mtl = t_end - t_start;
    //    qDebug("%f | %f", v, mtl);
    //    qDebug("%f %f", t_start, t_end);

    //    double xd = (root.bounds.max - root.bounds.min).mag();
    //    double scale_val = 10000.0; // xd * root.tree.value ;
    double scale_val = root.tree.total;
    //    qDebug("%f %f %f", v, scale_val, xd);

    //    qDebug("%f %f", v, scale_val);
    float s = hyptan(v / scale_val);
    //    float k = hyptan(t_end / nat_dist);
    //    qDebug("%f %f -> %f %f", t_start / nat_dist, t_end / nat_dist, s, k);
    return FColor(1. - s, 1. - s, 1., 0.5);
}
