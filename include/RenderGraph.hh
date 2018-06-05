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
class RenderGraphTask;
struct ViewData_s;
typedef struct ViewData_s ViewData;

typedef struct {
    QRgb *colors;
    double *distances;
} FlatData;

class Context {
public:
    Context(const ViewData &d, QSharedPointer<QImage> i, int nt);
    ~Context();
    const ViewData *viewdata;
    FlatData *partData;
    FlatData flatData;
    QSharedPointer<QImage> image;
    int renderno;
    int nthreads;
    volatile int abort_flag;
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
    void queueNext(int, int);
    void progUpdate(int, int);

    friend class RenderGraphTask;

private:
    void doQueue(int);

    QSharedPointer<Context> context;
    QElapsedTimer *timer;
    QThreadPool *pool;
    float progress;
    int seqno;

    typedef struct {
        QVector<int> reqs;
        QRect pxdm;
        int shard;
        int layer;
        bool inprogress;
    } Task;
    QMap<int, Task> tasks;
};

class RenderGraphTask : public QRunnable {
public:
    RenderGraphTask(QRect p, RenderGraph &h, QSharedPointer<Context> c, int id);
    virtual ~RenderGraphTask();
    int id;

protected:
    QRect pixels;
    RenderGraph &home;
    QSharedPointer<Context> ctx;
};
