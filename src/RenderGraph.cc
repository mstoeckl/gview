#include "RenderGraph.hh"
#include "RenderTask.hh"
#include "RenderWorker.hh"

#include <QElapsedTimer>
#include <QImage>
#include <QSignalMapper>
#include <QThread>
#include <QThreadPool>

class RenderGraphHelper : public QRunnable {
public:
    RenderGraphHelper(QSharedPointer<RenderGraphNode> n,
                      QSharedPointer<Context> c, RenderGraph *h)
        : home(h), node(n), context(c) {
        setAutoDelete(true);
    }
    virtual ~RenderGraphHelper() {}

    virtual void run() {
        if (context->abort_flag) {
            return;
        }
        node->run(&*context);
        if (context->abort_flag) {
            return;
        }
        node->status = RenderGraphNode::kComplete;
        QMetaObject::invokeMethod(home, "queueNext", Qt::QueuedConnection,
                                  Q_ARG(int, context->renderno));
    }

private:
    RenderGraph *home;
    QSharedPointer<RenderGraphNode> node;
    QSharedPointer<Context> context;
};

Context::Context(const ViewData &d, QSharedPointer<QImage> im, int nt,
                 int seqno)
    : renderno(seqno), nthreads(nt), w(im->width()), h(im->height()) {
    viewdata = new ViewData(d);
    image = im;
    abort_flag = false;

    partData = new FlatData[nthreads];
    for (int k = 0; k < nthreads; k++) {
        partData[k].colors = new QRgb[w * h];
        partData[k].distances = new double[w * h];
        partData[k].blank = false;
    }
    flatData.colors = new QRgb[w * h];
    flatData.distances = new double[w * h];
    flatData.blank = false;
    raydata = new RayPoint[w * h];
    intersection_store = new Intersection *[nthreads * 4]();
}

Context::~Context() {
    delete viewdata;
    for (int k = 0; k < nthreads; k++) {
        delete[] partData[k].colors;
        delete[] partData[k].distances;
    }
    delete[] partData;
    delete[] flatData.colors;
    delete[] flatData.distances;
    delete[] raydata;
    for (int i = 0; i < nthreads * 4; i++) {
        if (intersection_store[i]) {
            free(intersection_store[i]);
        }
    }
    delete[] intersection_store;
}

RenderGraph::RenderGraph(int nthreads) {
    pool = new QThreadPool();
    pool->setMaxThreadCount(nthreads);

    timer = new QElapsedTimer();

    seqno = 0;
}

RenderGraph::~RenderGraph() {
    // Destruction requires cleanup hints
    if (!context.isNull()) {
        context->abort_flag = true;
    }
    // Reset graph; nodes are reference counted
    tasks.clear();
    context = QSharedPointer<Context>();
    pool->deleteLater();

    delete timer;
}

void RenderGraph::start(QSharedPointer<QImage> im, const ViewData &vd) {
    timer->start();

    int w = im->width();
    int h = im->height();

    int nthreads = pool->maxThreadCount();

    context = QSharedPointer<Context>(new Context(vd, im, nthreads, seqno));
    seqno++;

    progress = 0.f;

    tasks.clear();

    QRect fullWindow = QRect(0, 0, w + 1, h + 1);
    QVector<QRect> panes;

    const int Ny = 2 * nthreads;
    const int Nx = 2;
    for (int i = 0; i < Ny; i++) {
        int yl = (i * h / Ny);
        int yh = i == Ny - 1 ? h : ((i + 1) * h / Ny);
        for (int k = 0; k < Nx; k++) {
            int xl = (k * w / Nx);
            int xh = k == Nx - 1 ? w : ((k + 1) * w / Nx);
            panes.append(QRect(xl, yl, xh - xl + 1, yh - yl + 1));
        }
    }

    if (0) {
        // Unbuffered mode!
        RenderMergeTask *merge = new RenderMergeTask(fullWindow);
        for (int shard = 0; shard < nthreads; shard++) {
            RenderTrackTask *track = new RenderTrackTask(fullWindow, shard);
            merge->reqs.push_back(track);
            tasks.push_back(QSharedPointer<RenderGraphNode>(track));
        }
        tasks.push_back(QSharedPointer<RenderGraphNode>(merge));
        for (int i = 0; i < panes.size(); i++) {
            RenderRayTask *ray = new RenderRayTask(panes[i]);
            ray->reqs.push_back(merge);
            tasks.push_back(QSharedPointer<RenderGraphNode>(ray));
        }
    } else {
        // Buffered mode!
        RenderMergeTask *merge = new RenderMergeTask(fullWindow);
        RenderColorTask *color = new RenderColorTask(fullWindow);
        for (int shard = 0; shard < nthreads; shard++) {
            RenderTrackTask *track = new RenderTrackTask(fullWindow, shard);
            merge->reqs.push_back(track);
            tasks.push_back(QSharedPointer<RenderGraphNode>(track));
        }
        color->reqs.push_back(merge);
        for (int i = 0; i < panes.size(); i++) {
            RenderRayBufferTask *raybuf = new RenderRayBufferTask(panes[i], i);
            color->reqs.push_back(raybuf);
            tasks.push_back(QSharedPointer<RenderGraphNode>(raybuf));
        }
        tasks.push_back(QSharedPointer<RenderGraphNode>(merge));
        tasks.push_back(QSharedPointer<RenderGraphNode>(color));
    }

    for (QSharedPointer<RenderGraphNode> j : tasks) {
        j->status = RenderGraphNode::kWaiting;

        if (j->reqs.size() == 0) {
            doQueue(j);
        }
    }
}

void RenderGraph::doQueue(QSharedPointer<RenderGraphNode> node) {
    node->status = RenderGraphNode::kActive;
    if (1) {
        pool->start(new RenderGraphHelper(node, context, this));
    } else {
        RenderGraphHelper(node, context, this).run();
    }
}

void RenderGraph::abort() {
    if (context.isNull()) {
        return;
    }
    // Fast terminate currently running actions
    context->abort_flag = true;

    // Reset graph; nodes are reference counted
    tasks.clear();
    context = QSharedPointer<Context>();

    emit aborted();
}

void RenderGraph::queueNext(int renderno) {
    if (context.isNull() || renderno != context->renderno) {
        return;
    }

    QVector<QSharedPointer<RenderGraphNode>> nseeds;
    int n_incomplete = 0;
    for (QSharedPointer<RenderGraphNode> m : tasks) {
        if (m->status != RenderGraphNode::kComplete) {
            n_incomplete++;
        }

        if (m->status != RenderGraphNode::kWaiting) {
            continue;
        }

        bool reqsleft = false;
        for (RenderGraphNode *p : m->reqs) {
            if (p->status != RenderGraphNode::kComplete) {
                reqsleft = true;
                break;
            }
        }
        if (!reqsleft) {
            nseeds.append(m);
        }
    }

    if (n_incomplete <= 0) {
        // Reset graph
        context = QSharedPointer<Context>();

        qreal elapsed = 1e-9 * timer->nsecsElapsed();
        emit done(elapsed);
    } else {
        for (QSharedPointer<RenderGraphNode> t : nseeds) {
            doQueue(t);
        }
    }

    // Update progress
    float op = progress;
    progress += 1. / tasks.size();
    if (int(100 * op) != int(100 * progress)) {
        emit progressed(int(100 * progress));
    }
}

RenderGraphNode::RenderGraphNode(const char *type) {
    status = kWaiting;
    name = type;
}

RenderGraphNode::~RenderGraphNode() {}
