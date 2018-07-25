/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#include <QRect>
#include <QRgb>
#include <QSharedPointer>
#include <QVector>

class Context;
typedef struct RayPoint_s RayPoint;
typedef struct Intersection_s Intersection;

class FlatData {
public:
    QRgb *colors;
    double *distances;
    bool blank;
    int w, h;
    FlatData(int iw, int ih) {
        w = iw;
        h = ih;
        colors = (QRgb *)malloc(w * h * sizeof(QRgb));
        distances = (double *)malloc(w * h * sizeof(double));
    }
    ~FlatData() {
        free(colors);
        free(distances);
    }
};

class VoxData {
public:
    const int w, h;
    bool blank;
    double **voxtrails;
    VoxData(int iw, int ih) : w(iw), h(ih) {
        voxtrails = (double **)malloc(w * h * sizeof(double *));
        blank = true;
    }
    ~VoxData() { free(voxtrails); }
};

/**
 * RenderGraphNode is a single render graph node to compute
 * some quantity once. It aborts the calculation when no nodes
 * depend on it. Hence, when switching targets, add a new target before
 * removing the old once.
 */
class RenderGraphNode {
public:
    RenderGraphNode(const char *type);
    virtual ~RenderGraphNode();

    typedef enum { kWaiting, kActive, kComplete } WorkStatus;

    virtual void run(Context *c) = 0;

    // Shared pointers used to permit any-order unlinking
    void addDependency(QSharedPointer<RenderGraphNode>);
    void request() { nconsumers++; }
    void unrequest() { nconsumers--; }
    // True if waiting to run and all dependencies complete
    bool isReady() const;

public:
    volatile WorkStatus status;
    const char *const name;

    // Links to nodes the depends on
    QVector<QSharedPointer<RenderGraphNode>> reqs;
    // Number of nodes for which the result matters
protected:
    friend class RenderGraphHelper;
    volatile int nconsumers;
};

class RenderDummyTask : public RenderGraphNode {
public:
    RenderDummyTask();
    virtual void run(Context *) {}
};

class RenderRayTask : public RenderGraphNode {
public:
    RenderRayTask(QRect p, QSharedPointer<QImage> i,
                  QSharedPointer<FlatData> f);
    virtual void run(Context *);

private:
    QSharedPointer<QImage> image;
    QSharedPointer<FlatData> flat_data;
    const QRect domain;
};

class RenderTrackTask : public RenderGraphNode {
public:
    RenderTrackTask(QRect p, int s, int ns, QSharedPointer<FlatData> f);
    virtual void run(Context *);

    QSharedPointer<FlatData> part_data;

private:
    const QRect domain;
    const int shard, nshards;
};

class RenderMergeTask : public RenderGraphNode {
public:
    RenderMergeTask(QRect p, QVector<QSharedPointer<FlatData>> i,
                    QSharedPointer<FlatData> o);
    virtual void run(Context *);

    QSharedPointer<FlatData> flat_data;

private:
    QVector<QSharedPointer<FlatData>> part_data;
    const QRect domain;
};

class RenderRayBufferTask : public RenderGraphNode {
public:
    RenderRayBufferTask(QRect p, QSharedPointer<QVector<RayPoint>> r);
    virtual ~RenderRayBufferTask();
    virtual void run(Context *);

    QSharedPointer<QVector<RayPoint>> ray_data;

private:
    Intersection *intersection_store;
    const QRect domain;
};

class RenderColorTask : public RenderGraphNode {
public:
    RenderColorTask(QRect p, QSharedPointer<QVector<RayPoint>> r,
                    QSharedPointer<FlatData> f, QSharedPointer<VoxData> v,
                    QSharedPointer<QImage> i);
    virtual void run(Context *);

private:
    QSharedPointer<QVector<RayPoint>> ray_data;
    QSharedPointer<FlatData> flat_data;
    QSharedPointer<VoxData> vox_data;
    QSharedPointer<QImage> image;
    const QRect domain;
};

class RenderVoxelBufferTask : public RenderGraphNode {
public:
    RenderVoxelBufferTask(QRect p, QSharedPointer<VoxData> voxray_data,
                          QSharedPointer<QVector<RayPoint>> r);
    virtual ~RenderVoxelBufferTask();
    virtual void run(Context *);

    QSharedPointer<VoxData> voxray_data;

private:
    QSharedPointer<QVector<RayPoint>> ray_data;
    double *voxel_store;
    const QRect domain;
};
