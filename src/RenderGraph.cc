#include "RenderGraph.hh"
#include "RenderWorker.hh"

#include <QImage>
#include <QThreadPool>

Context::Context(const ViewData &d, QSharedPointer<QImage> i, int w, int h) {
    viewdata = new ViewData(d);
    colors = new QRgb[w * h];
    distances = new double[w * h];
    image = i;
    abort_flag = false;
}

Context::~Context() {
    delete viewdata;
    delete[] colors;
    delete[] distances;
}

RenderGraph::RenderGraph(int nthreads) {
    pool = new QThreadPool();
    pool->setMaxThreadCount(nthreads);

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
}

void RenderGraph::start(QSharedPointer<QImage> im, const ViewData &vd) {
    int w = im->width();
    int h = im->height();

    context = QSharedPointer<Context>(new Context(vd, im, w, h));
    context->renderno = seqno;
    progress = 0.f;

    tasks = QMap<int, Task>();

    int Ny = pool->maxThreadCount();
    int Nx = pool->maxThreadCount();
    int z = 0;
    for (int i = 0; i < Ny; i++) {
        // Vertical slice (tracks)
        int yl = (i * h / Ny);
        int yh = i == Ny - 1 ? h : ((i + 1) * h / Ny);
        Task g;
        g.layer = 0;
        g.pxdm = QRect(0, yl, w + 1, yh - yl + 1);
        g.reqs = QVector<int>();
        g.inprogress = false;
        tasks[z] = g;
        for (int k = 0; k < Nx; k++) {
            // Horizontal slice (objects)
            int xl = (k * w / Nx);
            int xh = k == Nx - 1 ? w : ((k + 1) * w / Nx);
            Task j;
            j.layer = 1;
            j.pxdm = QRect(xl, yl, xh - xl + 1, yh - yl + 1);
            j.reqs = QVector<int>();
            j.reqs.append(z);
            j.inprogress = false;
            tasks[z + k + 1] = j;
        }
        z += Nx + 1;
    }
    QVector<int> seeds;
    for (int i = 0; i < tasks.size(); i++) {
        const Task &j = tasks[i];
        if (j.reqs.size() == 0) {
            seeds.append(i);
            if (seeds.size() >= pool->maxThreadCount()) {
                break;
            }
        }
    }
    // Tasks may complete instantaneously, so we select first, start second
    for (int i : seeds) {
        Task &j = tasks[i];
        RenderGraphTask *r;
        if (j.layer == 0) {
            r = new RenderTrackTask(j.pxdm, *this, context, i);
        } else {
            r = new RenderRayTask(j.pxdm, *this, context, i);
        }
        j.inprogress = true;
        // Blocking start
        pool->start(r);
    }
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
    int sel = -1;
    for (int i : tasks.keys()) {
        bool reqsleft = false;
        for (int j : tasks[i].reqs) {
            if (tasks.contains(j)) {
                reqsleft = true;
                break;
            }
        }
        if (!reqsleft && !tasks[i].inprogress) {
            sel = i;
            break;
        }
    }

    if (sel < 0 && tasks.size() <= 0) {
        // Reset graph
        tasks = QMap<int, Task>();
        context = QSharedPointer<Context>();
        seqno++;

        emit done();
    } else if (sel >= 0) {
        Task &j = tasks[sel];
        RenderGraphTask *r;
        if (j.layer == 0) {
            r = new RenderTrackTask(j.pxdm, *this, context, sel);
        } else {
            r = new RenderRayTask(j.pxdm, *this, context, sel);
        }
        j.inprogress = true;
        // Blocking start
        pool->start(r);
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
