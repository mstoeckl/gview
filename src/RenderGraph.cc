/* SPDX-License-Identifier: GPL-3.0-only */
#include "RenderGraph.hh"
#include "Nursery.hh"
#include "RenderTask.hh"
#include "RenderWorker.hh"

#include <QElapsedTimer>
#include <QImage>
#include <QSignalMapper>

class RenderGraphHelper : public QRunnable {
public:
    RenderGraphHelper(QSharedPointer<RenderGraphNode> n,
                      QSharedPointer<Context> c, RenderGraph *h)
        : home(h), node(n), context(c) {
        setAutoDelete(true);
    }
    virtual ~RenderGraphHelper() {}

    virtual void run() {
        if (node->nconsumers <= 0) {
            return;
        }
        node->run(&*context);
        if (node->nconsumers <= 0) {
            return;
        }
        node->status = RenderGraphNode::kComplete;
        QMetaObject::invokeMethod(home, "queueNext", Qt::QueuedConnection,
                                  Q_ARG(RenderGraphNode *, node.data()));
    }

private:
    RenderGraph *home;
    QSharedPointer<RenderGraphNode> node;
    QSharedPointer<Context> context;
};

Context::Context(const ViewData &d, const GridSpec &g)
    : grid(g), viewdata(new ViewData(d)) {}
Context::~Context() { delete viewdata; }

RenderGraph::RenderGraph(int nthreads) : QObject(), nursery(nthreads) {
    timer = new QElapsedTimer();

    qRegisterMetaType<RenderGraphNode *>("RenderGraphNode*");
}

RenderGraph::~RenderGraph() {
    // Reset graph; nodes are reference counted and will self-terminate
    tasks.clear();
    task_color.clear();
    task_merge.clear();
    task_ray.clear();
    task_raybuf.clear();
    task_track.clear();

    context = QSharedPointer<Context>();

    delete timer;
}

void RenderGraph::start(QSharedPointer<QImage> im, const GridSpec &grid,
                        const ViewData &vd, int changed) {
    timer->start();
    bool color_changed = changed & CHANGE_COLOR;
    bool geo_changed = changed & CHANGE_GEO;
    bool tracks_changed = changed & CHANGE_TRACK;
    bool oneshot = changed & CHANGE_ONESHOT;

    int nthreads = nursery.idealThreadCount();

    context = QSharedPointer<Context>(new Context(vd, grid));

    progress = 0.f;

    const int w = grid.sampleWidth(), h = grid.sampleHeight();
    QRect fullWindow = QRect(0, 0, w + 1, h + 1);
    QVector<QRect> panes;
    QVector<QRect> rows;

    const int Ny = nthreads;
    const int Nx = 5;
    for (int i = 0; i < Ny; i++) {
        int yl = (i * h / Ny);
        int yh = i == Ny - 1 ? h : ((i + 1) * h / Ny);
        rows.append(QRect(0, yl, w + 1, yh - yl + 1));
        for (int k = 0; k < Nx; k++) {
            int xl = (k * w / Nx);
            int xh = k == Nx - 1 ? w : ((k + 1) * w / Nx);
            panes.append(QRect(xl, yl, xh - xl + 1, yh - yl + 1));
        }
    }

    // Replace graph parts (and dependents) which need changing
    if (oneshot) {
        task_voxel.clear();
        task_color.clear();
        task_merge.clear();
        task_ray.clear();
        task_raybuf.clear();
        task_track.clear();

        // Create large storage buffers
        QVector<QSharedPointer<FlatData>> part_data_list;
        for (int i = 0; i < nthreads; i++) {
            QSharedPointer<FlatData> part_data(new FlatData(w, h));
            part_data_list.push_back(part_data);
        }
        QSharedPointer<FlatData> flat_data(new FlatData(w, h));

        task_merge.push_back(QSharedPointer<RenderGraphNode>(
            new RenderMergeTask(fullWindow, part_data_list, flat_data)));
        for (int i = 0; i < nthreads; i++) {
            QSharedPointer<RenderGraphNode> track(new RenderTrackTask(
                fullWindow, i, nthreads, part_data_list[i]));
            for (const QSharedPointer<RenderGraphNode> &p : task_merge) {
                p->addDependency(track);
            }
            task_track.push_back(track);
        }
        for (int i = 0; i < panes.size(); i++) {
            QSharedPointer<RenderGraphNode> ray(
                new RenderRayTask(panes[i], im, flat_data));
            for (const QSharedPointer<RenderGraphNode> &p : task_merge) {
                ray->addDependency(p);
            }
            task_ray.push_back(ray);
        }
        task_final = QSharedPointer<RenderGraphNode>(new RenderDummyTask());
        for (const QSharedPointer<RenderGraphNode> &p : task_ray) {
            task_final->addDependency(p);
        }
        task_final->request();

        tasks = task_track + task_merge + task_ray;
        tasks.push_back(task_final);
    } else {
        // Only replace a segment if predecessors were replaced
        task_ray.clear();
        if (tracks_changed || task_track.size() <= 0) {
            tracks_changed = true;
            task_track.clear();
            task_merge.clear();

            QVector<QSharedPointer<FlatData>> part_data_list;
            for (int i = 0; i < nthreads; i++) {
                QSharedPointer<FlatData> part_data(new FlatData(w, h));
                part_data_list.push_back(part_data);
            }
            QSharedPointer<FlatData> flat_data(new FlatData(w, h));
            task_merge.push_back(QSharedPointer<RenderGraphNode>(
                new RenderMergeTask(fullWindow, part_data_list, flat_data)));
            for (int i = 0; i < nthreads; i++) {
                QSharedPointer<RenderGraphNode> track(new RenderTrackTask(
                    fullWindow, i, nthreads, part_data_list[i]));
                task_track.push_back(track);
                for (const QSharedPointer<RenderGraphNode> &p : task_merge) {
                    p->addDependency(track);
                }
            }
        }

        if (geo_changed || task_raybuf.size() <= 0) {
            geo_changed = true;
            task_raybuf.clear();

            QSharedPointer<QVector<RayPoint>> raypts(
                new QVector<RayPoint>(w * h));
            for (int i = 0; i < panes.size(); i++) {
                QSharedPointer<RenderGraphNode> ray(
                    new RenderRayBufferTask(panes[i], raypts));
                task_raybuf.push_back(ray);
            }
        }

        if (geo_changed || tracks_changed || task_voxel.size() <= 0) {
            task_voxel.clear();

            QSharedPointer<QVector<RayPoint>> ray_data =
                reinterpret_cast<RenderRayBufferTask *>(task_raybuf[0].data())
                    ->ray_data;
            QSharedPointer<VoxData> voxdata(new VoxData(w, h));
            for (int i = 0; i < panes.size(); i++) {
                QSharedPointer<RenderGraphNode> vox(
                    new RenderVoxelBufferTask(panes[i], voxdata, ray_data));
                // TODO: use rectangle intersection to find overlaps
                for (const QSharedPointer<RenderGraphNode> &p : task_raybuf) {
                    vox->addDependency(p);
                }
                task_voxel.push_back(vox);
            }
        }

        if (1) {
            (void)color_changed;
            task_color.clear();

            QSharedPointer<QVector<RayPoint>> ray_data =
                reinterpret_cast<RenderRayBufferTask *>(task_raybuf[0].data())
                    ->ray_data;
            QSharedPointer<VoxData> vox_data =
                reinterpret_cast<RenderVoxelBufferTask *>(task_voxel[0].data())
                    ->voxray_data;
            QSharedPointer<FlatData> flat_data =
                reinterpret_cast<RenderMergeTask *>(task_merge[0].data())
                    ->flat_data;
            if (ray_data.isNull() || flat_data.isNull() || vox_data.isNull()) {
                qFatal("Null color dependents");
            }

            for (int i = 0; i < rows.size(); i++) {
                QSharedPointer<RenderGraphNode> color(new RenderColorTask(
                    rows[i], ray_data, flat_data, vox_data, im));
                // TODO: use rectangle intersection test to find overlaps
                for (const QSharedPointer<RenderGraphNode> &p : task_merge) {
                    color->addDependency(p);
                }
                for (const QSharedPointer<RenderGraphNode> &p : task_raybuf) {
                    color->addDependency(p);
                }
                for (const QSharedPointer<RenderGraphNode> &p : task_voxel) {
                    color->addDependency(p);
                }
                task_color.push_back(color);
            }
        }

        task_final = QSharedPointer<RenderGraphNode>(new RenderDummyTask());
        for (const QSharedPointer<RenderGraphNode> &p : task_color) {
            task_final->addDependency(p);
        }
        task_final->request();

        tasks = task_track + task_merge + task_raybuf + task_color + task_voxel;
        tasks.push_back(task_final);
    }

    for (QSharedPointer<RenderGraphNode> j : tasks) {
        if (j.isNull()) {
            qFatal("Task was null, quitting");
        }
        if (j->isReady()) {
            doQueue(j);
        }
    }
}

void RenderGraph::doQueue(QSharedPointer<RenderGraphNode> node) {
    node->status = RenderGraphNode::kActive;
    if (1) {
        nursery.startSoon(new RenderGraphHelper(node, context, this),
                          kDeleteAfter);
    } else {
        RenderGraphHelper(node, context, this).run();
    }
}

void RenderGraph::abort() {
    // Remove the previous graph head
    task_final->unrequest();

    emit aborted();
}

void RenderGraph::queueNext(RenderGraphNode *node) {
    if (node == task_final.data()) {
        qreal elapsed = 1e-9 * timer->nsecsElapsed();
        emit done(elapsed);
    }

    QVector<QSharedPointer<RenderGraphNode>> nseeds;
    for (QSharedPointer<RenderGraphNode> m : tasks) {
        if (m->isReady()) {
            doQueue(m);
        }
    }

    // Update progress
    float op = progress;
    progress += 1. / tasks.size();
    if (int(100 * op) != int(100 * progress)) {
        emit progressed(int(100 * progress));
    }
}
