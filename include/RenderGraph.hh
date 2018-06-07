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
struct ViewData_s;
typedef struct ViewData_s ViewData;
typedef struct RayPoint_s RayPoint;
typedef struct Intersection_s Intersection;

typedef struct {
    QRgb *colors;
    double *distances;
    bool blank;
} FlatData;

class RenderGraphNode;

class Context {
public:
    Context(const ViewData &d, QSharedPointer<QImage> i, int nt, int seqno);
    ~Context();

    const int renderno;
    const int nthreads;
    const int w, h;

    const ViewData *viewdata;
    FlatData *partData;
    FlatData flatData;
    QSharedPointer<QImage> image;
    RayPoint *raydata;
    Intersection **intersection_store;

    volatile int abort_flag;

private:
    Context(const Context &) = delete;
};

class RenderGraph : public QObject {
    Q_OBJECT
public:
    RenderGraph(int nthreads);
    ~RenderGraph();

public slots:

    void start(QSharedPointer<QImage> i, const ViewData &vd);
    void abort();

signals:
    void done(qreal elapsed_secs);
    void aborted();

    void progressed(int);

protected slots:
    void queueNext(int);

    friend class RenderGraphTask;

private:
    void doQueue(QSharedPointer<RenderGraphNode> node);

    QSharedPointer<Context> context;
    QElapsedTimer *timer;
    QThreadPool *pool;
    float progress;
    int seqno;

    QVector<QSharedPointer<RenderGraphNode>> tasks;
};
