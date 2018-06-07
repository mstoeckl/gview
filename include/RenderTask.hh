#pragma once

#include <QRect>
#include <QVector>

class Context;

class RenderGraphNode {
public:
    RenderGraphNode(const char *type);
    virtual ~RenderGraphNode();

    typedef enum { kWaiting, kActive, kComplete } WorkStatus;

    virtual void run(Context *c) const = 0;

public:
    QVector<RenderGraphNode *> reqs;
    volatile WorkStatus status;
    const char *name;
};

class RenderRayTask : public RenderGraphNode {
public:
    RenderRayTask(QRect p);
    virtual void run(Context *) const;

private:
    QRect domain;
};

class RenderTrackTask : public RenderGraphNode {
public:
    RenderTrackTask(QRect p, int shard);
    virtual void run(Context *) const;

private:
    QRect domain;
    int shard;
};

class RenderMergeTask : public RenderGraphNode {
public:
    RenderMergeTask(QRect p);
    virtual void run(Context *) const;

private:
    QRect domain;
};

class RenderRayBufferTask : public RenderGraphNode {
public:
    RenderRayBufferTask(QRect p, int shard);
    virtual void run(Context *) const;

private:
    QRect domain;
    int shard;
};

class RenderColorTask : public RenderGraphNode {
public:
    RenderColorTask(QRect p);
    virtual void run(Context *) const;

private:
    QRect domain;
};
