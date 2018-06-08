#pragma once

#include <QMap>
#include <QObject>
#include <QRect>
#include <QRgb>
#include <QRunnable>
#include <QSharedPointer>
#include <QVector>

class QElapsedTimer;
class QThreadPool;
class Context;
typedef struct ViewData_s ViewData;

class RenderGraphNode;
class RenderMergeTask;
class RenderRayTask;
class RenderTrackTask;
class RenderRayBufferTask;
class RenderColorTask;
class RenderDummyTask;

class Context {
public:
    Context(const ViewData &d, int iw, int ih);
    ~Context();

    const int w, h;
    const ViewData *viewdata;

private:
    Context(const Context &) = delete;
};

class RenderGraph : public QObject {
    Q_OBJECT
public:
    RenderGraph(int nthreads);
    ~RenderGraph();

public slots:

    void start(QSharedPointer<QImage> i, const ViewData &vd, int changed);
    void abort();

signals:
    void done(qreal elapsed_secs);
    void aborted();

    void progressed(int);

protected slots:
    void queueNext(RenderGraphNode *);

    friend class RenderGraphTask;

private:
    void doQueue(QSharedPointer<RenderGraphNode> node);

    QSharedPointer<Context> context;
    QElapsedTimer *timer;
    QThreadPool *pool;
    float progress;

    QVector<QSharedPointer<RenderGraphNode>> task_track;
    QVector<QSharedPointer<RenderGraphNode>> task_ray;
    QVector<QSharedPointer<RenderGraphNode>> task_merge;
    QVector<QSharedPointer<RenderGraphNode>> task_color;
    QVector<QSharedPointer<RenderGraphNode>> task_raybuf;
    QSharedPointer<RenderGraphNode> task_final;

    QVector<QSharedPointer<RenderGraphNode>> tasks;
};
