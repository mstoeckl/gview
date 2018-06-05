#include "RenderGraph.hh"
#include "RenderWorker.hh"

#include <QElapsedTimer>
#include <QImage>
#include <QThreadPool>

Context::Context(const ViewData &d, QSharedPointer<QImage> im, int nt) {
    viewdata = new ViewData(d);
    int w = im->width(), h = im->height();
    image = im;
    abort_flag = false;
    nthreads = nt;

    partData = new FlatData[nthreads];
    for (int k = 0; k < nthreads; k++) {
        partData[k].colors = new QRgb[w * h];
        partData[k].distances = new double[w * h];
    }
    flatData.colors = new QRgb[w * h];
    flatData.distances = new double[w * h];
}

Context::~Context() {
    delete viewdata;
    for (int k = 0; k < nthreads; k++) {
        delete[] partData[k].colors;
        delete[] partData[k].distances;
    }
    delete[] partData;
    delete flatData.colors;
    delete flatData.distances;
}

RenderGraph::RenderGraph(int nthreads) {
    pool = new QThreadPool();
    pool->setMaxThreadCount(nthreads);

    timer = new QElapsedTimer();

    qRegisterMetaType<RenderGraphTask *>("RenderGraphTask*");

    seqno = 0;
}

RenderGraph::~RenderGraph() {
    // Destruction requires cleanup hints
    if (!context.isNull()) {
        context->abort_flag = true;
    }
    // Reset graph
    tasks = QMap<int, Task>();
    context = QSharedPointer<Context>();
    seqno++;
    pool->deleteLater();

    delete timer;
}

void RenderGraph::start(QSharedPointer<QImage> im, const ViewData &vd) {
    timer->start();

    int w = im->width();
    int h = im->height();

    int nthreads = pool->maxThreadCount();

    context = QSharedPointer<Context>(new Context(vd, im, nthreads));
    context->renderno = seqno;
    progress = 0.f;

    tasks = QMap<int, Task>();

    int Ny = pool->maxThreadCount();
    int Nx = pool->maxThreadCount();
    int z = 0;
    for (int i = 0; i < Ny; i++) {
        // Vertical slice (tracks)
        Task g;
        g.layer = 0;
        g.pxdm = QRect(0, 0, w + 1, h + 1);
        g.reqs = QVector<int>();
        g.shard = i;
        g.inprogress = false;
        tasks[z] = g;
        z++;
    }

    Task tm;
    tm.layer = 1;
    tm.shard = 0;
    tm.pxdm = QRect(0, 0, w + 1, h + 1);
    tm.reqs = QVector<int>();
    for (int i = 0; i < Ny; i++) {
        tm.reqs.push_back(i);
    }
    tm.inprogress = false;
    int zstar = z;
    tasks[zstar] = tm;
    z++;

    for (int i = 0; i < Ny; i++) {
        int yl = (i * h / Ny);
        int yh = i == Ny - 1 ? h : ((i + 1) * h / Ny);
        for (int k = 0; k < Nx; k++) {
            // Horizontal slice (objects)
            int xl = (k * w / Nx);
            int xh = k == Nx - 1 ? w : ((k + 1) * w / Nx);
            Task j;
            j.layer = 2;
            j.pxdm = QRect(xl, yl, xh - xl + 1, yh - yl + 1);
            j.reqs = QVector<int>();
            j.reqs.append(zstar);
            j.inprogress = false;
            tasks[z] = j;
            z++;
        }
    }
    QVector<int> seeds;
    for (int i = 0; i < tasks.size(); i++) {
        const Task &j = tasks[i];
        if (j.reqs.size() == 0) {
            seeds.append(i);
        }
    }
    // Tasks may complete instantaneously, so we select first, start second
    for (int i : seeds) {
        doQueue(i);
    }
}

void RenderGraph::doQueue(int i) {
    Task &j = tasks[i];
    RenderGraphTask *r;
    if (j.layer == 0) {
        r = new RenderTrackTask(j.pxdm, j.shard, *this, context, i);
    } else if (j.layer == 1) {
        r = new RenderMergeTask(j.pxdm, *this, context, i);
    } else {
        r = new RenderRayTask(j.pxdm, *this, context, i);
    }
    j.inprogress = true;
    // Blocking start
    pool->start(r);
}

void RenderGraph::abort() {
    if (context.isNull()) {
        return;
    }
    // Fast terminate currently running actions
    context->abort_flag = true;

    // Reset graph
    tasks = QMap<int, Task>();
    context = QSharedPointer<Context>();
    seqno++;

    emit aborted();
}

void RenderGraph::queueNext(int taskid, int rno) {
    if (context.isNull() || rno != context->renderno || context->abort_flag) {
        return;
    }
    // Removing task marks as complete
    tasks.remove(taskid);
    QVector<int> nseeds;
    for (int i : tasks.keys()) {
        bool reqsleft = false;
        for (int j : tasks[i].reqs) {
            if (tasks.contains(j)) {
                reqsleft = true;
                break;
            }
        }
        if (!reqsleft && !tasks[i].inprogress) {
            nseeds.append(i);
            break;
        }
    }

    if (nseeds.size() <= 0 && tasks.size() <= 0) {
        // Reset graph
        tasks = QMap<int, Task>();
        context = QSharedPointer<Context>();
        seqno++;

        qreal elapsed = 1e-9 * timer->nsecsElapsed();
        emit done(elapsed);
    } else {
        for (int z : nseeds) {
            doQueue(z);
        }
    }
}

void RenderGraph::progUpdate(int steps, int rno) {
    if (context.isNull() || rno != context->renderno || context->abort_flag) {
        return;
    }

    // TODO: figure out ratcheting when
    // tasks complete in an efficient way...
    float op = progress;
    progress +=
        1. / (steps * pool->maxThreadCount() * (1 + pool->maxThreadCount()));
    if (int(100 * op) != int(100 * progress)) {
        emit progressed(int(100 * progress));
    }
}

RenderGraphTask::RenderGraphTask(QRect p, RenderGraph &h,
                                 QSharedPointer<Context> c, int i)
    : home(h) {
    pixels = p;
    ctx = c;
    id = i;
    this->setAutoDelete(true);
}

RenderGraphTask::~RenderGraphTask() {}
