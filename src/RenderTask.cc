#include "RenderTask.hh"

#include "General.hh"
#include "RenderWorker.hh"

#include <QImage>
#include <geomdefs.hh>

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
    int w = ctx->w;
    int h = ctx->h;
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

    const G4double radius = 0.8 / mind;
    // Zero construct by default
    ElemMutables *mutables = new ElemMutables[nelements]();
    int iter = 1;

    // reqs: of type?
    const FlatData *flatData = &(*flat_data);

    const G4ThreeVector &forward = forwardDirection(d->orientation);

    const int M = 30;
    Intersection ints[M];
    Intersection aints[M];
    int ndevs[M + 1];

    for (int i = yl; i < yh; i++) {
        QRgb *pts = reinterpret_cast<QRgb *>(image->scanLine(i));
        for (int j = xl; j < xh; j++) {
            if (nconsumers <= 0) {
                delete[] mutables;
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

            QPointF pt((j - w / 2.) / (2. * mind), (i - h / 2.) / (2. * mind));
            RayPoint r = rayAtPoint(pt, radius, forward, *d, iter, mutables,
                                    ints, aints, M, ndevs);
            QRgb color = colorForRay(r, trackcol, trackdist, *d, pt, forward);
            pts[j] = color;
        }
    }
    if (d->level_of_detail <= -1) {
        recursivelyPrintNCalls(d->elements, mutables);
    }

    delete[] mutables;
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
    int w = ctx->w;
    int h = ctx->h;

    int mind = std::min(w, h);

    int xl = domain.left();
    int xh = domain.right();
    int yl = domain.top();
    int yh = domain.bottom();

    const G4double radius = 0.8 / mind;
    // Zero construct by default
    ElemMutables *mutables = new ElemMutables[nelements]();
    int iter = 1;

    const G4ThreeVector &forward = forwardDirection(d->orientation);

    const int M = 30;
    Intersection ints[M];
    Intersection aints[M];
    int ndevs[M + 1];

    QVector<RayPoint> &rayData = *ray_data;

    size_t max_usage = M * sizeof(Intersection) * (yh - yl + 1) * (xh - xl + 1);
    Intersection *const storage = (Intersection *)malloc(max_usage);
    size_t istored = 0;

    for (int i = yl; i < yh; i++) {
        for (int j = xl; j < xh; j++) {
            if (nconsumers <= 0) {
                delete[] mutables;
                return;
            }

            QPointF pt((j - w / 2.) / (2. * mind), (i - h / 2.) / (2. * mind));
            RayPoint r = rayAtPoint(pt, radius, forward, *d, iter, mutables,
                                    ints, aints, M, ndevs);

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
    if (d->level_of_detail <= -1) {
        recursivelyPrintNCalls(d->elements, mutables);
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

    delete[] mutables;
}

static inline double iproject(const G4ThreeVector &point,
                              const G4ThreeVector &dcamera,
                              const G4double &dscale,
                              const G4RotationMatrix &ori, int w, int h,
                              int *dx, int *dy) {
    G4ThreeVector local = ori * (point - dcamera);
    G4double idscale = 1 / dscale;
    double fy = local.y() * idscale;
    double fx = local.z() * idscale;
    double off = local.x();
    int dmind = 2 * std::min(w, h);
    *dx = int(0.5 * w + fx * dmind);
    *dy = int(0.5 * h + fy * dmind);
    return off;
}

static void trackColors(const TrackHeader &h, const TrackPoint &pa,
                        const TrackPoint &pb, const G4ThreeVector &a,
                        const G4ThreeVector &b, const G4ThreeVector &forward,
                        FColor &ca, FColor &cb, float &wa, float &wb) {
    G4ThreeVector normal = b - a;
    double coa = (normal * forward) / normal.mag();

    //    double mtime = double(pa.time) * (1. - interp) * double(pb.time) *
    //    double ptime = (int(mtime / CLHEP::ns * 1024.) % 1024) / 1024.;

    double aloge = std::log10(double(pa.energy)) / 2.0;
    double patime = aloge - std::floor(aloge);
    ca = FColor(QColor::fromHsvF(patime, 0.8, 1.0 - 0.5 * std::abs(coa)).rgb());

    double bloge = std::log10(double(pb.energy)) / 2.0;
    double pbtime = bloge - std::floor(bloge);
    cb = FColor(QColor::fromHsvF(pbtime, 0.8, 1.0 - 0.5 * std::abs(coa)).rgb());

    if (h.ptype == 11) {
        // e-
        wa = 1.0f;
        wb = 1.0f;
    } else if (h.ptype == 22) {
        // gamma
        wa = 0.7f;
        wb = 0.7f;
    } else {
        // others (optiphoton)
        wa = 0.5f;
        wb = 0.5f;
    }
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

    int w = ctx->w;
    int h = ctx->h;
    const ViewData &d = *ctx->viewdata;
    const TrackData &t = d.tracks;

    size_t ntracks = t.getNTracks();
    size_t tfrom = (ntracks * shard) / nshards;
    size_t tupto = (ntracks * (shard + 1)) / nshards;
    if (tfrom == tupto) {
        return;
    }

    // Might have tracks, so we fill background
    double *dists = part_data->distances;
    QRgb *colors = part_data->colors;
    if (part_data->w != ctx->w || part_data->h != ctx->h) {
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
        int sdx, sdy;
        double soff;
        soff =
            iproject(spos, d.camera, d.scale, d.orientation, w, h, &sdx, &sdy);

        // Fast entire track clipping heuristic
        int rr = int(2 * mind * metadata[z].ballRadius / d.scale);
        if (sdx + rr <= xl || sdx - rr > xh || sdy + rr <= yl ||
            sdy - rr > yh) {
            continue;
        }

        part_data->blank = false;

        for (int j = 1; j < header.npts; j++) {
            // Adopt previous late point as early point
            TrackPoint ep = sp;
            int edx = sdx, edy = sdy;
            double eoff = soff;
            G4ThreeVector epos = spos;
            // Calculate new late point
            sp = tp[j];
            spos = G4ThreeVector(sp.x, sp.y, sp.z);
            soff = iproject(spos, d.camera, d.scale, d.orientation, w, h, &sdx,
                            &sdy);

            if ((edx < xl && sdx < xl) || (edy < yl && sdy < yl) ||
                (edx >= xh && sdx >= xh) || (edy >= yh && sdy >= yh)) {
                continue;
            }

            // NOTE: single point lines are acceptable
            float dy, dx;
            int n;
            if (std::abs(sdx - edx) > std::abs(sdy - edy)) {
                n = std::abs(sdx - edx) + 1;
                dx = sdx - edx > 0 ? 1.0 : -1.0;
                dy = float(sdy - edy) / float(n);
            } else {
                n = std::abs(sdy - edy) + 1;
                dy = sdy - edy > 0 ? 1.0 : -1.0;
                dx = float(sdx - edx) / float(n);
            }
            float ix = edx;
            float iy = edy;

            int uxl = dx == 0.0f ? 0 : int((xl - edx) / dx);
            int uyl = dy == 0.0f ? 0 : int((yl - edy) / dy);
            int uxh = dx == 0.0f ? INT_MAX : int((xh - edx) / dx);
            int uyh = dy == 0.0f ? INT_MAX : int((yh - edy) / dy);

            // Extra buffers are just in case
            int near = std::max(std::min(uxl, uxh), std::min(uyl, uyh)) - 1;
            int far = std::min(std::max(uxl, uxh), std::max(uyl, uyh)) + 1;

            if (std::max(near, 0) > std::min(far, n)) {
                continue;
            }

            float w1, w2;
            FColor c1, c2;
            trackColors(header, ep, sp, epos, spos, d.orientation.rowX(), c1,
                        c2, w1, w2);

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
    int w = ctx->w;
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
                                 QSharedPointer<QImage> i)
    : RenderGraphNode("color"), ray_data(r), flat_data(f), image(i), domain(p) {
}
void RenderColorTask::run(Context *ctx) {
    int xl = domain.left();
    int xh = domain.right();
    int yl = domain.top();
    int yh = domain.bottom();
    int w = ctx->w;
    int h = ctx->h;
    int mind = std::min(w, h);

    const FlatData &flatData = *flat_data;
    const ViewData *viewData = ctx->viewdata;
    const QVector<RayPoint> &rayData = *ray_data;
    const G4ThreeVector &forward = forwardDirection(viewData->orientation);

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
            QPointF pt((j - w / 2.) / (2. * mind), (i - h / 2.) / (2. * mind));
            QRgb color =
                colorForRay(ray, trackcol, trackdist, *viewData, pt, forward);
            pts[j] = color;
        }
    }
}
