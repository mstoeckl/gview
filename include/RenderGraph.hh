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

class Context {
public:
    Context(const ViewData &d, QSharedPointer<QImage> i, int w, int h);
    ~Context();
    const ViewData *viewdata;
    QRgb *colors;
    double *distances;
    QSharedPointer<QImage> image;
    int renderno;
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
    QSharedPointer<Context> context;
    QElapsedTimer *timer;
    QThreadPool *pool;
    float progress;
    int seqno;

    typedef struct {
        QVector<int> reqs;
        QRect pxdm;
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
