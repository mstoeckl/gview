#include "RenderTask.hh"

#include "General.hh"
#include "RenderWorker.hh"

#include <QImage>
#include <geomdefs.hh>

RenderRayTask::RenderRayTask(QRect p) : RenderGraphNode("ray"), domain(p) {}

void RenderRayTask::run(Context *ctx) const {
    const ViewData *d = ctx->viewdata;

    if (!d->elements.solid) {
        return;
    }
    int treedepth;
    int nelements;
    countTree(d->elements, treedepth, nelements);
    // ^ TODO: allocate traceray buffers to match!
    if (treedepth > 10) {
        qFatal("Excessive tree depth, fatal!");
    }
    int w = ctx->w;
    int h = ctx->h;
    if (ctx->image->width() != w || ctx->image->height() != h) {
        qFatal("Image size mismatch, %d %d vs %d %d", ctx->image->width(),
               ctx->image->height(), w, h);
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

    const FlatData *flatData = &ctx->flatData;

    const G4ThreeVector &forward = forwardDirection(d->orientation);

    const int M = 30;
    Intersection ints[M];
    Intersection aints[M];
    int ndevs[M + 1];

    for (int i = yl; i < yh; i++) {
        QRgb *pts = reinterpret_cast<QRgb *>(ctx->image->scanLine(i));
        for (int j = xl; j < xh; j++) {
            if (ctx->abort_flag) {
                delete[] mutables;
                return;
            }

            // For merging at correct depth
            int sidx = i * w + j;
            QRgb trackcol = flatData->colors[sidx];
            G4double trackdist = flatData->distances[sidx];

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

RenderRayBufferTask::RenderRayBufferTask(QRect p, int s)
    : RenderGraphNode("buffer"), domain(p), shard(s) {}
void RenderRayBufferTask::run(Context *ctx) const {
    const ViewData *d = ctx->viewdata;
    if (!d->elements.solid) {
        return;
    }
    int treedepth;
    int nelements;
    countTree(d->elements, treedepth, nelements);
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

    RayPoint *image = ctx->raydata;
    size_t max_usage = M * sizeof(Intersection) * (yh - yl + 1) * (xh - xl + 1);
    Intersection *const storage_start = (Intersection *)malloc(max_usage);
    Intersection *storage = storage_start;

    for (int i = yl; i < yh; i++) {
        for (int j = xl; j < xh; j++) {
            if (ctx->abort_flag) {
                delete[] mutables;
                return;
            }

            QPointF pt((j - w / 2.) / (2. * mind), (i - h / 2.) / (2. * mind));
            RayPoint r = rayAtPoint(pt, radius, forward, *d, iter, mutables,
                                    ints, aints, M, ndevs);

            // Store ray, and fill its intersections in the buffer
            image[i * w + j] = r;
            for (int k = 0; k < r.N; k++) {
                storage[k] = r.intersections[k];
            }
            image[i * w + j].intersections = storage;
            storage = &storage[r.N];
        }
    }
    if (d->level_of_detail <= -1) {
        recursivelyPrintNCalls(d->elements, mutables);
    }

    if (storage != storage_start) {
        // If any point had an intersection
        Intersection *verif = (Intersection *)realloc(
            (void *)storage_start, (char *)storage - (char *)storage_start);
        if (verif != storage_start) {
            qFatal("realloc should only shrink %lu %lu",
                   (char *)storage - (char *)storage_start, max_usage);
        }
        ctx->intersection_store[shard] = storage_start;
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

RenderTrackTask::RenderTrackTask(QRect p, int s)
    : RenderGraphNode("track"), domain(p), shard(s) {}

void RenderTrackTask::run(Context *ctx) const {
    if (ctx->abort_flag) {
        return;
    }

    // Blank unless we make a change
    ctx->partData[shard].blank = true;

    int xl = domain.left();
    int xh = domain.right();
    int yl = domain.top();
    int yh = domain.bottom();

    int w = ctx->w;
    int h = ctx->h;
    const ViewData &d = *ctx->viewdata;
    const TrackData &t = d.tracks;

    size_t ntracks = t.getNTracks();
    size_t tfrom = (ntracks * shard) / ctx->nthreads;
    size_t tupto = (ntracks * (shard + 1)) / ctx->nthreads;
    if (tfrom == tupto) {
        return;
    }

    // Might have tracks, so we fill background
    double *dists = ctx->partData[shard].distances;
    QRgb *colors = ctx->partData[shard].colors;
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
        if (ctx->abort_flag) {
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

        ctx->partData[shard].blank = false;

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

RenderMergeTask::RenderMergeTask(QRect p)
    : RenderGraphNode("merge"), domain(p) {}

void RenderMergeTask::run(Context *ctx) const {
    // Minimum merge over all partData
    Context &c = *ctx;
    if (c.abort_flag) {
        return;
    }

    int xl = domain.left();
    int xh = domain.right();
    int yl = domain.top();
    int yh = domain.bottom();
    int w = ctx->w;
    for (int y = yl; y < yh; y++) {
        for (int x = xl; x < xh; x++) {
            // Scale up to be beyond anything traceRay produces
            c.flatData.distances[y * w + x] = 4 * kInfinity;
            c.flatData.colors[y * w + x] = qRgba(255, 255, 255, 0.);
        }
    }

    for (int k = 0; k < c.nthreads; k++) {
        if (c.partData[k].blank) {
            continue;
        }
        for (int y = yl; y < yh; y++) {
            for (int x = xl; x < xh; x++) {
                int i = y * w + x;
                if (c.partData[k].distances[i] < c.flatData.distances[i]) {
                    c.flatData.distances[i] = c.partData[k].distances[i];
                    c.flatData.colors[i] = c.partData[k].colors[i];
                }
            }
        }
    }
}

RenderColorTask::RenderColorTask(QRect p)
    : RenderGraphNode("color"), domain(p) {}
void RenderColorTask::run(Context *ctx) const {
    int xl = domain.left();
    int xh = domain.right();
    int yl = domain.top();
    int yh = domain.bottom();
    int w = ctx->w;
    int h = ctx->h;
    int mind = std::min(w, h);

    const FlatData *flatData = &ctx->flatData;
    const ViewData *viewData = ctx->viewdata;
    const RayPoint *rayData = ctx->raydata;
    const G4ThreeVector &forward = forwardDirection(viewData->orientation);

    for (int i = yl; i < yh; i++) {
        QRgb *pts = reinterpret_cast<QRgb *>(ctx->image->scanLine(i));
        for (int j = xl; j < xh; j++) {
            if (ctx->abort_flag) {
                return;
            }

            // For merging at correct depth
            int sidx = i * w + j;
            QRgb trackcol = flatData->colors[sidx];
            G4double trackdist = flatData->distances[sidx];
            const RayPoint &ray = rayData[sidx];
            QPointF pt((j - w / 2.) / (2. * mind), (i - h / 2.) / (2. * mind));
            QRgb color =
                colorForRay(ray, trackcol, trackdist, *viewData, pt, forward);
            pts[j] = color;
        }
    }
}
