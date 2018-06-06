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

typedef struct {
    QRgb *colors;
    double *distances;
    bool blank;
} FlatData;

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
