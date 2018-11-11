/* SPDX-License-Identifier: GPL-3.0-only */
#include "RenderTask.hh"

#include "General.hh"
#include "Navigator.hh"
#include "Octree.hh"
#include "RenderWorker.hh"
#include "Shaders.hh"

#include <QImage>
#include <geomdefs.hh>

#define MAX_BUFFER_DEPTH 30

RenderGraphNode::RenderGraphNode(const char *type)
    : status(kWaiting), name(type), reqs(), nconsumers(0) {}

RenderGraphNode::~RenderGraphNode() {
    for (QSharedPointer<RenderGraphNode> &p : reqs) {
        p->nconsumers--;
    }
}

void RenderGraphNode::addDependency(QSharedPointer<RenderGraphNode> p) {
    if (p.isNull()) {
        qFatal("Dependencies must not be null for task %s", name);
    }
    p->nconsumers++;
    reqs.push_back(p);
}
bool RenderGraphNode::isReady() const {
    if (status != kWaiting)
        return false;
    for (const QSharedPointer<RenderGraphNode> &p : reqs) {
        if (p->status != kComplete)
            return false;
    }
    return true;
}
RenderDummyTask::RenderDummyTask() : RenderGraphNode("dummy") {}

RenderRayTask::RenderRayTask(QRect p, QSharedPointer<QImage> i,
                             QSharedPointer<FlatData> f)
    : RenderGraphNode("ray"), image(i), flat_data(f), domain(p) {}

void RenderRayTask::run(Context *ctx) {
    const ViewData *d = ctx->viewdata;

    if (!d->elements[0].solid) {
        return;
    }
    int treedepth, nelements;
    countTree(d->elements, 0, treedepth, nelements);
    // ^ TODO: allocate traceray buffers to match!
    if (treedepth > 10) {
        qFatal("Excessive tree depth, fatal!");
    }
    const int w = ctx->grid.sampleWidth(), h = ctx->grid.sampleHeight();
    if (image->width() != w || image->height() != h) {
        qFatal("Image size mismatch, %d %d vs %d %d", image->width(),
               image->height(), w, h);
    }
#if 0
        QTime t = QTime::currentTime();
        qDebug("render lod %d w %d h %d | %d of %d", d.level_of_detail, w, h, slice,
               nslices);
#endif

    int mind = w > h ? h : w;

    int xl = domain.left();
    int xh = domain.right();
    int yl = domain.top();
    int yh = domain.bottom();

    const G4double radius = 3.0 / mind;
    // Zero construct by default
    Navigator *nav = Navigator::create(*d, d->navigator);
    GeoShader &shader = *getGeoShader(d->gshader);

    // reqs: of type?
    const FlatData *flatData = &(*flat_data);

    const G4ThreeVector &forward = forwardDirection(d->orientation);

    Intersection ints[MAX_BUFFER_DEPTH];
    Intersection aints[MAX_BUFFER_DEPTH];
    int ndevs[MAX_BUFFER_DEPTH + 1];

    for (int i = yl; i < yh; i++) {
        QRgb *pts = reinterpret_cast<QRgb *>(image->scanLine(i));
        for (int j = xl; j < xh; j++) {
            if (nconsumers <= 0) {
                delete nav;
                return;
            }

            // For merging at correct depth
            int sidx = i * w + j;

            QRgb trackcol;
            G4double trackdist;
            if (flatData->blank) {
                trackdist = 4 * kInfinity;
                trackcol = qRgba(255, 255, 255, 0);
            } else {
                trackcol = flatData->colors[sidx];
                trackdist = flatData->distances[sidx];
            }

            QPointF pt(ctx->grid.toViewCoord(j, i));
            RayPoint r = rayAtPoint(*nav, pt, radius, forward, *d, ints, aints,
                                    MAX_BUFFER_DEPTH, ndevs);
            FColor color = shader(r, FColor(trackcol), trackdist, NULL, *d, pt,
                                  forward, true);
            pts[j] = color.rgba();
        }
    }

    delete nav;
}

RenderRayBufferTask::RenderRayBufferTask(QRect p,
                                         QSharedPointer<QVector<RayPoint>> r)
    : RenderGraphNode("buffer"), ray_data(r), intersection_store(NULL),
      domain(p) {}
RenderRayBufferTask::~RenderRayBufferTask() {
    if (intersection_store)
        free(intersection_store);
}
void RenderRayBufferTask::run(Context *ctx) {
    const ViewData *d = ctx->viewdata;
    if (!d->elements[0].solid) {
        return;
    }
    int treedepth;
    int nelements;
    countTree(d->elements, 0, treedepth, nelements);
    // ^ TODO: allocate traceray buffers to match!
    if (treedepth > 10) {
        qFatal("Excessive tree depth, fatal!");
    }
    const int w = ctx->grid.sampleWidth(), h = ctx->grid.sampleHeight();

    int mind = std::min(w, h);

    int xl = domain.left();
    int xh = domain.right();
    int yl = domain.top();
    int yh = domain.bottom();

    const G4double radius = 3.0 / mind;
    // Zero construct by default
    Navigator *nav = Navigator::create(*d, d->navigator);

    const G4ThreeVector &forward = forwardDirection(d->orientation);

    Intersection ints[MAX_BUFFER_DEPTH];
    Intersection aints[MAX_BUFFER_DEPTH];
    int ndevs[MAX_BUFFER_DEPTH + 1];

    QVector<RayPoint> &rayData = *ray_data;

    size_t max_usage =
        MAX_BUFFER_DEPTH * sizeof(Intersection) * (yh - yl + 1) * (xh - xl + 1);
    Intersection *const storage = (Intersection *)malloc(max_usage);
    size_t istored = 0;

    for (int i = yl; i < yh; i++) {
        for (int j = xl; j < xh; j++) {
            if (nconsumers <= 0) {
                free(storage);
                delete nav;
                return;
            }

            QPointF pt(ctx->grid.toViewCoord(j, i));
            RayPoint r = rayAtPoint(*nav, pt, radius, forward, *d, ints, aints,
                                    MAX_BUFFER_DEPTH, ndevs);

            // Store ray, and copy its intersections to the buffer
            rayData[i * w + j] = r;
            if (r.N > 0) {
                for (int k = 0; k < r.N; k++) {
                    storage[istored + k] = r.intersections[k];
                }
                rayData[i * w + j].intersections = &storage[istored];
                istored += r.N;
            } else {
                rayData[i * w + j].intersections = NULL;
            }
        }
    }

    if (istored > 0) {
        // If any point had an intersection
        Intersection *verif = (Intersection *)realloc(
            (void *)storage, istored * sizeof(Intersection));
        if (verif != storage) {
            // Fix up pointers if the memory is moved
            for (int i = yl; i < yh; i++) {
                for (int j = xl; j < xh; j++) {
                    rayData[i * w + j].intersections += verif - storage;
                }
            }
        }
        intersection_store = verif;
    } else {
        free(storage);
        intersection_store = NULL;
    }

    delete nav;
}

static inline QPoint iproject(const G4ThreeVector &point,
                              const G4ThreeVector &dcamera,
                              const G4double &dscale,
                              const G4RotationMatrix &ori, const GridSpec &grid,
                              double &off) {
    // Map to F in smallest viewport box containing [-1,1]x[-1,1]
    G4ThreeVector local = ori * (point - dcamera);
    G4double idscale = 1 / dscale;
    QPointF F(local.z() * idscale, local.y() * idscale);
    off = local.x();
    QPointF sample_coord = grid.toSampleCoord(F);
    return QPoint(std::lrint(sample_coord.x()), std::lrint(sample_coord.y()));
}

RenderTrackTask::RenderTrackTask(QRect p, int s, int ns,
                                 QSharedPointer<FlatData> f)
    : RenderGraphNode("track"), part_data(f), domain(p), shard(s), nshards(ns) {
}

void RenderTrackTask::run(Context *ctx) {
    if (nconsumers <= 0) {
        return;
    }

    // Blank unless we make a change
    part_data->blank = true;

    int xl = domain.left();
    int xh = domain.right();
    int yl = domain.top();
    int yh = domain.bottom();

    const int w = ctx->grid.sampleWidth(), h = ctx->grid.sampleHeight();
    const ViewData &d = *ctx->viewdata;
    const TrackData &t = d.tracks;
    TrackShader &shader = *getTrackShader(d.tshader);

    size_t ntracks = t.getNTracks();
    size_t tfrom = (ntracks * shard) / nshards;
    size_t tupto = (ntracks * (shard + 1)) / nshards;
    if (tfrom == tupto) {
        return;
    }

    // Might have tracks, so we fill background
    double *dists = part_data->distances;
    QRgb *colors = part_data->colors;
    if (part_data->w != w || part_data->h != h) {
        qFatal("Bad size in RenderTrackTask::run");
    }

    for (int y = yl; y < yh; y++) {
        for (int x = xl; x < xh; x++) {
            // Scale up to be beyond anything traceRay produces
            dists[y * w + x] = 4 * kInfinity;
            colors[y * w + x] = qRgba(255, 255, 255, 0.);
        }
    }

    int mind = w > h ? h : w;
    float rad = std::max(mind / 800.0f, 0.01f);

    const TrackBlock *blocks = t.getBlocks();
    const TrackMetaData *metadata = t.getMeta();
    size_t i = 0;
    for (size_t z = 0; z < tfrom; z++) {
        i += blocks[i].h.npts + 1;
    }

    for (size_t z = tfrom; z < tupto; z++) {
        if (nconsumers <= 0) {
            return;
        }

        // Jump to next track
        const TrackHeader &header = blocks[i].h;
        const TrackPoint *tp = (TrackPoint *)&blocks[i + 1];
        i += header.npts + 1;

        // Project calculations 1 point ahead
        TrackPoint sp = tp[0];
        G4ThreeVector spos(sp.x, sp.y, sp.z);
        double soff;
        QPoint sd(
            iproject(spos, d.camera, d.scale, d.orientation, ctx->grid, soff));

        // Fast entire track clipping heuristic
        int rr = int(2 * mind * metadata[z].ballRadius / d.scale);
        if (sd.x() + rr <= xl || sd.x() - rr > xh || sd.y() + rr <= yl ||
            sd.y() - rr > yh) {
            continue;
        }

        part_data->blank = false;

        for (int j = 1; j < header.npts; j++) {
            // Adopt previous late point as early point
            TrackPoint ep = sp;
            QPoint ed(sd);
            double eoff = soff;
            G4ThreeVector epos = spos;
            // Calculate new late point
            sp = tp[j];
            spos = G4ThreeVector(sp.x, sp.y, sp.z);
            sd = iproject(spos, d.camera, d.scale, d.orientation, ctx->grid,
                          soff);

            if ((ed.x() < xl && sd.x() < xl) || (ed.y() < yl && sd.y() < yl) ||
                (ed.x() >= xh && sd.x() >= xh) ||
                (ed.y() >= yh && sd.y() >= yh)) {
                continue;
            }

            // NOTE: single point lines are acceptable
            float dy, dx;
            int n;
            QPoint smd = sd - ed;
            if (std::abs(smd.x()) > std::abs(smd.y())) {
                n = std::abs(smd.x()) + 1;
                dx = smd.x() > 0 ? 1.0 : -1.0;
                dy = float(smd.y()) / float(n);
            } else {
                n = std::abs(smd.y()) + 1;
                dy = smd.y() > 0 ? 1.0 : -1.0;
                dx = float(smd.x()) / float(n);
            }
            float ix = ed.x();
            float iy = ed.y();

            int uxl = dx == 0.0f ? 0 : int((xl - ed.x()) / dx);
            int uyl = dy == 0.0f ? 0 : int((yl - ed.y()) / dy);
            int uxh = dx == 0.0f ? INT_MAX : int((xh - ed.x()) / dx);
            int uyh = dy == 0.0f ? INT_MAX : int((yh - ed.y()) / dy);

            // Extra buffers are just in case
            int near = std::max(std::min(uxl, uxh), std::min(uyl, uyh)) - 1;
            int far = std::min(std::max(uxl, uxh), std::max(uyl, uyh)) + 1;

            if (std::max(near, 0) > std::min(far, n)) {
                continue;
            }

            float w1, w2;
            FColor c1, c2;
            shader(header, ep, sp, epos, spos, d.orientation.rowX(), c1, c2, w1,
                   w2);

            for (int s = std::max(near, 0); s < std::min(far, n); s++) {
                float x = ix + s * dx;
                float y = iy + s * dy;

                float b = n >= 2 ? s / float(n - 1) : 0.5f;
                double off = eoff * double(1.f - b) + soff * double(b);
                float rws = w1 * (1.0f - b) + w2 * b;
                FColor col = FColor::blend(c1, c2, b);

                // TODO: move to the merger phase, & multithread that!
                // (so, instead, just draw the center point?)

                /* Draw circle */
                float spot = rad * rws;
                int ly = int(std::lrint(y - spot));
                int hy = int(std::lrint(y + spot));
                for (int ddy = std::max(yl, ly); ddy <= std::min(yh - 1, hy);
                     ddy++) {
                    float lr = std::sqrt(std::max(
                        0.f, spot * spot - float((ddy - y) * (ddy - y))));
                    int lx = int(std::lrint(x - lr));
                    int hx = int(std::lrint(x + lr));
                    for (int ddx = std::max(xl, lx);
                         ddx <= std::min(xh - 1, hx); ddx++) {
                        int sidx = ddy * w + ddx;
                        if (dists[sidx] > off) {
                            dists[sidx] = off;
                            colors[sidx] = col.rgba();
                        }
                    }
                }
            }
        }
    }
}

RenderMergeTask::RenderMergeTask(QRect p, QVector<QSharedPointer<FlatData>> i,
                                 QSharedPointer<FlatData> o)
    : RenderGraphNode("merge"), flat_data(o), part_data(i), domain(p) {}

void RenderMergeTask::run(Context *ctx) {
    // Minimum merge over all partData
    if (nconsumers <= 0) {
        return;
    }

    bool all_blank = true;
    for (int k = 0; k < part_data.size(); k++) {
        all_blank &= part_data[k]->blank;
    }
    if (all_blank) {
        flat_data->blank = true;
        return;
    }
    flat_data->blank = false;

    int xl = domain.left();
    int xh = domain.right();
    int yl = domain.top();
    int yh = domain.bottom();
    const int w = ctx->grid.sampleWidth();
    FlatData &fd = *flat_data;
    for (int y = yl; y < yh; y++) {
        for (int x = xl; x < xh; x++) {
            // Scale up to be beyond anything traceRay produces
            fd.distances[y * w + x] = 4 * kInfinity;
            fd.colors[y * w + x] = qRgba(255, 255, 255, 0.);
        }
    }

    for (int k = 0; k < part_data.size(); k++) {
        FlatData &pd = *part_data[k];

        if (pd.blank) {
            continue;
        }
        for (int y = yl; y < yh; y++) {
            for (int x = xl; x < xh; x++) {
                int i = y * w + x;
                if (pd.distances[i] < fd.distances[i]) {
                    fd.distances[i] = pd.distances[i];
                    fd.colors[i] = pd.colors[i];
                }
            }
        }
    }
}

RenderColorTask::RenderColorTask(QRect p, QSharedPointer<QVector<RayPoint>> r,
                                 QSharedPointer<FlatData> f,
                                 QSharedPointer<VoxData> v,
                                 QSharedPointer<QImage> i)
    : RenderGraphNode("color"), ray_data(r), flat_data(f), vox_data(v),
      image(i), domain(p) {}
void RenderColorTask::run(Context *ctx) {
    int xl = domain.left();
    int xh = domain.right();
    int yl = domain.top();
    int yh = domain.bottom();
    const int w = ctx->grid.sampleWidth();

    const FlatData &flatData = *flat_data;
    const ViewData *viewData = ctx->viewdata;
    const QVector<RayPoint> &rayData = *ray_data;
    const VoxData &voxData = *vox_data;
    const G4ThreeVector &forward = forwardDirection(viewData->orientation);
    GeoShader &shader = *getGeoShader(viewData->gshader);

    for (int i = yl; i < yh; i++) {
        QRgb *pts = reinterpret_cast<QRgb *>(image->scanLine(i));
        for (int j = xl; j < xh; j++) {
            if (nconsumers <= 0) {
                return;
            }

            // For merging at correct depth
            int sidx = i * w + j;
            QRgb trackcol;
            G4double trackdist;
            if (flatData.blank) {
                trackdist = 4 * kInfinity;
                trackcol = qRgba(255, 255, 255, 0);
            } else {
                trackcol = flatData.colors[sidx];
                trackdist = flatData.distances[sidx];
            }

            const RayPoint &ray = rayData[sidx];
            QPointF pt(ctx->grid.toViewCoord(j, i));
            FColor color =
                shader(ray, FColor(trackcol), trackdist,
                       voxData.blank ? NULL : voxData.voxtrails[sidx],
                       *viewData, pt, forward, true);
            pts[j] = color.rgba();
        }
    }
}

RenderVoxelBufferTask::RenderVoxelBufferTask(
    QRect p, QSharedPointer<VoxData> iv, QSharedPointer<QVector<RayPoint>> r)
    : RenderGraphNode("voxel"), voxray_data(iv), ray_data(r),
      voxel_store(nullptr), domain(p) {}
RenderVoxelBufferTask::~RenderVoxelBufferTask() { free(voxel_store); }
void RenderVoxelBufferTask::run(Context *ctx) {
    int xl = domain.left();
    int xh = domain.right();
    int yl = domain.top();
    int yh = domain.bottom();

    const int w = ctx->grid.sampleWidth();
    const ViewData &vd = *ctx->viewdata;
    VoxData &voxData = *voxray_data;
    if (!ctx->viewdata->tracks.getOctree()) {
        // If no octree, do no work
        voxData.blank = true;
        return;
    }
    const OctreeRoot &octree = *ctx->viewdata->tracks.getOctree();
    voxData.blank = false;

    const QVector<RayPoint> &rayData = *ray_data;

    size_t max_usage =
        (MAX_BUFFER_DEPTH + 1) * sizeof(double) * (yh - yl + 1) * (xh - xl + 1);
    double *const storage = (double *)malloc(max_usage);
    size_t istored = 0;

    for (int i = yl; i < yh; i++) {
        for (int j = xl; j < xh; j++) {
            if (nconsumers <= 0) {
                return;
            }

            int sidx = i * w + j;

            // For merging at correct depth
            QPointF pt(ctx->grid.toViewCoord(j, i));
            const G4ThreeVector &fr = initPoint(pt, vd);
            const G4ThreeVector &dir = forwardDirection(vd.orientation);

            const RayPoint &r = rayData[sidx];
            double *start = &storage[istored];
            voxData.voxtrails[sidx] = start;
            double da = r.N ? r.intersections[0].dist : kInfinity;
            start[0] = traceDensityRay(octree, fr, dir, da);

            for (int z = 0; z < r.N - 1; z++) {
                double d0 = r.intersections[z].dist,
                       d1 = r.intersections[z + 1].dist;
                start[z + 1] =
                    traceDensityRay(octree, fr + dir * d1, dir, d1 - d0);
            }
            if (r.N) {
                start[r.N] = traceDensityRay(
                    octree, fr + dir * r.intersections[r.N - 1].dist, dir,
                    kInfinity);
            }
            istored += r.N + 1;
        }
    }

    //  Realloc and shrink everything. TODO: make a wrapper class!
    if (istored > 0) {
        // If any point had an intersection
        double *verif =
            (double *)realloc((void *)storage, istored * sizeof(double));
        if (verif != storage) {
            // Fix up pointers if the memory is moved
            for (int i = yl; i < yh; i++) {
                for (int j = xl; j < xh; j++) {
                    voxData.voxtrails[i * w + j] += verif - storage;
                }
            }
        }
        voxel_store = verif;
    } else {
        free(storage);
        voxel_store = NULL;
    }
}
