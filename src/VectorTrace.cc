/* SPDX-License-Identifier: GPL-3.0-only */
#include "VectorTrace.hh"

#include "Navigator.hh"
#include "Shaders.hh"

#include <QFile>
#include <QFileInfo>
#include <QPainter>
#include <QPen>
#include <QProcess>
#include <QSet>
#include <QTemporaryFile>
#include <QTextStream>
#include <QTime>

#include <G4Material.hh>

static FColor randColor() {
    return FColor(randint(65536) / 65535., randint(65536) / 65535.,
                  randint(65536) / 65535., 1.0);
}

static RayPoint dummy_point() {
    RayPoint r;
    r.N = 0;
    r.front_clipped = false;
    r.back_clipped = false;
    r.intersections = NULL;
    return r;
}
static RayPoint copy_alloc_point(RayPoint ray) {
    if (ray.intersections) {
        Intersection *oi = ray.intersections;
        ray.intersections = new Intersection[ray.N];
        for (int i = 0; i < ray.N; i++) {
            ray.intersections[i] = oi[i];
        }
    }
    return ray;
}

RenderPoint::RenderPoint()
    : coords(std::numeric_limits<double>::quiet_NaN(),
             std::numeric_limits<double>::quiet_NaN()),
      ray(dummy_point()) {
    ideal_color = FColor();
    region_class = -1;
    show_point = true;
}
RenderPoint::RenderPoint(QPointF spot, const RayPoint &rpt)
    : coords(spot), ray(copy_alloc_point(rpt)) {
    ideal_color = FColor();
    region_class = -1;
    subregion_class = -1;
    show_point = true;
}
RenderPoint::~RenderPoint() {
    if (ray.intersections)
        delete[] ray.intersections;
}
RenderPoint::RenderPoint(const RenderPoint &other)
    : coords(other.coords), ray(copy_alloc_point(other.ray)) {
    ideal_color = other.ideal_color;
    region_class = other.region_class;
    subregion_class = other.subregion_class;
    show_point = other.show_point;
}
RenderPoint &RenderPoint::operator=(RenderPoint copy_of_other) {
    swap(copy_of_other);
    return *this;
}
void RenderPoint::swap(RenderPoint &other) {
    std::swap(ray, other.ray);
    std::swap(region_class, other.region_class);
    std::swap(subregion_class, other.subregion_class);
    std::swap(coords, other.coords);
    std::swap(ideal_color, other.ideal_color);
    std::swap(show_point, other.show_point);
}

static QPointF grid_coord_to_point(const QPoint &pt, const QSize &grid_size) {
    GridSpec gs(grid_size.width(), grid_size.height(), 1);
    QPointF vpt = gs.toViewCoord(pt.x(), pt.y());
    return vpt;
}
static QPointF point_to_grid_coord(const QPointF &pt, const QSize &grid_size) {
    GridSpec gs(grid_size.width(), grid_size.height(), 1);
    QPointF gpt = gs.toSampleCoord(pt);
    return gpt;
}
static void recsetColorsByMaterial(std::vector<Element> &elts,
                                   std::vector<VColor> &color_table,
                                   QMap<QString, VColor> &color_map,
                                   QMap<QString, int> &idx_map, int idx) {
    // We hard-code color associations to be more consistent
    Element &elem = elts[idx];
    QString key(elem.material->GetName().c_str());
    if (!color_map.count(key)) {
        VColor c(randColor().rgba());
        color_map[key] = c;
        qWarning("Material `%s` assigned random color: rgb=(%f %f %f)",
                 key.toUtf8().constData(), c.redF(), c.blueF(), c.greenF());
    }

    if (!idx_map.count(key)) {
        int n = idx_map.size();
        idx_map[key] = n;
        color_table.push_back(color_map[key]);
    }
    elem.ccode = idx_map[key];
    for (int jdx : elem.children) {
        recsetColorsByMaterial(elts, color_table, color_map, idx_map, jdx);
    }
}

VectorTracer::VectorTracer(ViewData vd, TrackData td,
                           const QString &target_file, QObject *parent)
    : QObject(parent), view_data(vd), track_data(td) {
    step_next = Steps::sGrid;
    file_name = target_file;
    nqueries = 0;

    recolor();

    grid_size = QSize(100, 100);

    grid_points = NULL;
    grid_nclasses = 0;
}

VectorTracer::~VectorTracer() {
    if (grid_points)
        delete[] grid_points;
}

void VectorTracer::recolor() {
    element_colors.clear();
    for (VColor v : view_data.color_table) {
        element_colors.append(FColor(v.rgb()));
    }
}

QImage VectorTracer::preview(const QSize &sz) {
    QRgb *data = new QRgb[sz.width() * sz.height()];

    Navigator *nav = Navigator::create(view_data, view_data.navigator);
    const QPointF &corner_a = grid_coord_to_point(QPoint(0, 0), sz);
    const QPointF &corner_b = grid_coord_to_point(
        QPoint(grid_size.width() - 1, grid_size.height() - 1), sz);
    bound_low = QPointF(std::min(corner_a.x(), corner_b.x()),
                        std::min(corner_a.y(), corner_b.y()));
    bound_high = QPointF(std::max(corner_a.x(), corner_b.x()),
                         std::max(corner_a.y(), corner_b.y()));

    for (int x = 0; x < sz.width(); x++) {
        for (int y = 0; y < sz.height(); y++) {
            int idx = y * sz.width() + x;
            QPointF qt = grid_coord_to_point(QPoint(x, y), sz);
            RenderPoint pt = queryPoint(qt, nav);
            data[idx] = calculateInteriorColor(pt).rgba();
        }
    }
    delete nav;

    QImage img((uchar *)data, sz.width(), sz.height(), QImage::Format_ARGB32);
    img = img.copy();
    delete[] data;

    return img;
}

void VectorTracer::renderFull() {
    while (step_next != Steps::sDone)
        renderStep();
}
void VectorTracer::renderStep() {
    switch (step_next) {
    case Steps::sDone:
        break;
    case Steps::sGrid:
        computeGrid();
        step_next = Steps::sEdges;
        break;
    case Steps::sEdges:
        computeEdges();
        step_next = Steps::sCreases;
        break;
    case Steps::sCreases:
        computeCreases();
        step_next = Steps::sGradients;
        break;
    case Steps::sGradients:
        computeGradients();
        step_next = Steps::sDone;
        break;
    }
}
void VectorTracer::reset(const QSize &size, const QString &target_name) {
    step_next = Steps::sGrid;
    grid_size = size;
    file_name = target_name;
    qDebug("Reset tracer: size (%d,%d), file %s", grid_size.width(),
           grid_size.height(), file_name.toUtf8().constData());
}
int VectorTracer::faildepth(const RenderPoint &a, const RenderPoint &b) {
    // Dummy points always match each other
    const RayPoint &ra = a.ray;
    const RayPoint &rb = b.ray;

    if (!ra.intersections && !rb.intersections)
        return -1;

    if (ra.front_clipped != rb.front_clipped) {
        return 0;
    }

    for (int i = 0; i < std::min(ra.N, rb.N); i++) {
        if (ra.intersections[i].ecode != rb.intersections[i].ecode) {
            return i;
        }
    }
    if (ra.N != rb.N) {
        return std::min(ra.N, rb.N);
    }
    if (!ra.intersections || !rb.intersections)
        return 0;

    if (ra.back_clipped != rb.back_clipped) {
        // note ra.N == rb.N
        return ra.N - 1;
    }

    // no failure
    return -1;
}

RenderPoint VectorTracer::queryPoint(QPointF spot,
                                     Navigator *thread_specific_nav) {
    nqueries++;

    const double eps = 1e-10;
    if (spot.x() < bound_low.x() - eps || spot.y() < bound_low.y() - eps ||
        spot.x() > bound_high.x() + eps || spot.y() > bound_high.y() + eps) {
        // Spot out of bounds; return dummy point
        RenderPoint r(spot, dummy_point());
        return r;
    }

    const int dlimit = 100;
    Intersection ints[dlimit + 1];

    RayPoint rpt = thread_specific_nav->traceRay(
        initPoint(spot, view_data), forwardDirection(view_data.orientation),
        ints, dlimit, view_data.force_opaque);

    return RenderPoint(spot, rpt);
}
RenderPoint VectorTracer::getPoint(QPoint p) {
    if (p.x() < 0 || p.y() < 0 || p.x() >= grid_size.width() ||
        p.y() >= grid_size.height()) {
        RenderPoint q(grid_coord_to_point(p, grid_size), dummy_point());
        return q;
    } else {
        return grid_points[p.x() * grid_size.height() + p.y()];
    }
}

/* Note that this is fully symmetric */
void VectorTracer::bracketEdge(const RenderPoint &initial_inside,
                               const RenderPoint &initial_outside,
                               RenderPoint *result_inside,
                               RenderPoint *result_outside,
                               Navigator *thread_navigator) {
    if (typematch(initial_inside, initial_outside)) {
        if (!initial_inside.ray.N && !initial_outside.ray.N) {
            // If both regions are empty, fall back to midpoint.
            // (This might happen when the scene is empty)
            QPointF p_mid =
                0.5 * (initial_inside.coords + initial_outside.coords);
            *result_inside = queryPoint(p_mid, thread_navigator);
            *result_outside = *result_inside;
            return;
        } else {
            qFatal("Bracket must cross class boundary");
        }
    }
    const int nsubdivisions = 20;
    *result_inside = initial_inside;
    *result_outside = initial_outside;
    for (int k = 0; k < nsubdivisions; k++) {
        QPointF p_mid = 0.5 * (result_inside->coords + result_outside->coords);
        RenderPoint test = queryPoint(p_mid, thread_navigator);
        if (typematch(test, *result_inside)) {
            *result_inside = test;
        } else {
            *result_outside = test;
        }
    }
}

void VectorTracer::bracketCrease(const RenderPoint &initial_inside,
                                 const RenderPoint &initial_outside,
                                 RenderPoint *result_inside,
                                 RenderPoint *result_outside,
                                 Navigator *thread_nav) {
    const double cos_alpha = std::cos(M_PI / 12);
    const double min_jump = 1e-3;
    // Step 1: determine cause of difference;
    bool is_jump = false;
    int cd = crease_depth(initial_inside, initial_outside, cos_alpha, min_jump,
                          &is_jump);
    if (cd < 0) {
        qFatal("Crease invariant fail");
        // if cd >= 0, then fni(*rin),fni(*rout) nontrivial
    }
    *result_inside = initial_inside;
    *result_outside = initial_outside;

    // Step 2: bisect on cause of difference
    Intersection cx_inside = result_inside->ray.intersections[cd];
    Intersection cx_outside = result_outside->ray.intersections[cd];
    const int nsubdivisions = 20;
    for (int i = 0; i < nsubdivisions; i++) {
        RenderPoint mid = queryPoint(
            0.5 * (result_inside->coords + result_outside->coords), thread_nav);
        Intersection cx_mid;
        if (mid.ray.N <= cd) {
            // Size failure for midpoint; abort
            return;
        }
        cx_mid = mid.ray.intersections[cd];
        if (is_jump) {
            if (std::abs(cx_mid.dist - cx_inside.dist) <
                std::abs(cx_mid.dist - cx_outside.dist)) {
                *result_inside = mid;
            } else {
                *result_outside = mid;
            }
        } else {
            if (G4ThreeVector(cx_mid.normal).dot(cx_inside.normal) >
                G4ThreeVector(cx_mid.normal).dot(cx_outside.normal)) {
                *result_inside = mid;
            } else {
                *result_outside = mid;
            }
        }
    }
}

FColor VectorTracer::calculateInteriorColor(const RenderPoint &pt) {
    const G4ThreeVector &forward = forwardDirection(view_data.orientation);
    GeoShader &shader = *getGeoShader(view_data.gshader);
    FColor color = shader(pt.ray, FColor(1., 1., 1., 0), 9e99, NULL, view_data,
                          pt.coords, forward, false);
    return color;
}
FColor VectorTracer::calculateBoundaryColor(const RenderPoint &inside,
                                            const RenderPoint &outside) {

    int lim = faildepth(inside, outside);
    if (lim < 0) {
        qFatal("Faildepth should be nonnegative over a boundary");
    }
    if (lim > std::min(inside.ray.N, outside.ray.N)) {
        qFatal("Faildepth should not be greater by more than one than either "
               "option");
    }

    const G4ThreeVector &forward = forwardDirection(view_data.orientation);
    GeoShader &shader = *getGeoShader(view_data.gshader);

    const RayPoint &ri = inside.ray;
    const RayPoint &ro = outside.ray;

    RayPoint combo;
    combo.N = lim + 1;
    combo.intersections = new Intersection[combo.N];
    for (int i = 0; i < combo.N - 1; i++) {
        combo.intersections[i].dist =
            (ri.intersections[i].dist + ro.intersections[i].dist) / 2;
        combo.intersections[i].ecode = ri.intersections[i].ecode;
        combo.intersections[i].normal =
            CompactNormal((G4ThreeVector(ri.intersections[i].normal) +
                           G4ThreeVector(ro.intersections[i].normal)) /
                          2);
    }
    combo.intersections[combo.N - 1].ecode = CODE_LINE;
    if (combo.N > ri.N) {
        combo.intersections[combo.N - 1].dist =
            ro.N ? ro.intersections[combo.N - 1].dist : 9e99;
        combo.intersections[combo.N - 1].normal =
            ro.N ? ro.intersections[combo.N - 1].normal : CompactNormal();
    } else if (combo.N > ro.N) {
        combo.intersections[combo.N - 1].dist =
            ri.N ? ri.intersections[combo.N - 1].dist : 9e99;
        combo.intersections[combo.N - 1].normal =
            ro.N ? ri.intersections[combo.N - 1].normal : CompactNormal();
    } else {
        combo.intersections[combo.N - 1].dist =
            (ri.intersections[combo.N - 1].dist +
             ro.intersections[combo.N - 1].dist) /
            2;
        combo.intersections[combo.N - 1].normal = CompactNormal(
            (G4ThreeVector(ri.intersections[combo.N - 1].normal) +
             G4ThreeVector(ro.intersections[combo.N - 1].normal)) /
            2);
    }

    // TODO: record plane index instead, and zero on disagreement
    combo.front_clipped = inside.ray.front_clipped && outside.ray.front_clipped;
    combo.front_clipped = inside.ray.back_clipped && outside.ray.back_clipped;

    QPointF mco = (inside.coords + outside.coords) / 2;
    FColor color = shader(combo, FColor(1., 1., 1., 0), 9e99, NULL, view_data,
                          mco, forward, false);

    delete[] combo.intersections;

    return color;
}

void VectorTracer::computeGrid() {
    const int W = grid_size.width(), H = grid_size.height();

    qDebug("Starting grid computation");
    if (grid_points)
        delete[] grid_points;

    // Calculate bounds. Note there may be inversions
    const QPointF &corner_a = grid_coord_to_point(QPoint(0, 0), grid_size);
    const QPointF &corner_b = grid_coord_to_point(
        QPoint(grid_size.width() - 1, grid_size.height() - 1), grid_size);
    bound_low = QPointF(std::min(corner_a.x(), corner_b.x()),
                        std::min(corner_a.y(), corner_b.y()));
    bound_high = QPointF(std::max(corner_a.x(), corner_b.x()),
                         std::max(corner_a.y(), corner_b.y()));

    grid_points = new RenderPoint[W * H];
    Navigator *nav = Navigator::create(view_data, view_data.navigator);
    for (int i = 0; i < W; i++) {
        for (int j = 0; j < H; j++) {
            int idx = i * H + j;
            QPointF spot = grid_coord_to_point(QPoint(i, j), grid_size);
            grid_points[idx] = queryPoint(spot, nav);
            grid_points[idx].ideal_color =
                calculateInteriorColor(grid_points[idx]);
        }
    }
    delete nav;
    qDebug("Grid computed");

    qDebug("Identifying equivalence classes");
    typedef struct {
        // true = checked already
        bool neigh[4];
        int x, y;
    } nset;
    QVector<nset> stack;
    grid_nclasses = -1;
    const int xdeltas[4] = {1, 0, -1, 0};
    const int ydeltas[4] = {0, 1, 0, -1};
    for (int x = 0; x < W; x++) {
        for (int y = 0; y < H; y++) {
            int idx = x * H + y;
            if (grid_points[idx].region_class < 0) {
                grid_nclasses++;
                grid_points[idx].region_class = grid_nclasses;
                nset n;
                n.neigh[0] = false;
                n.neigh[1] = false;
                n.neigh[2] = false;
                n.neigh[3] = false;
                n.x = x;
                n.y = y;
                stack.push_back(n);
                while (stack.size()) {
                    nset s = stack.last();
                    bool found = false;
                    for (int i = 0; i < 4; i++) {
                        if (s.neigh[i])
                            continue;

                        s.neigh[i] = true;

                        int dx = s.x + xdeltas[i], dy = s.y + ydeltas[i];
                        if (dx < 0 || dx >= W || dy < 0 || dy >= H) {
                            continue;
                        }
                        int sdx = s.x * H + s.y;
                        int ddx = dx * H + dy;
                        if (grid_points[ddx].region_class >= 0)
                            continue;

                        if (typematch(grid_points[sdx], grid_points[ddx])) {
                            grid_points[ddx].region_class = grid_nclasses;
                            s.neigh[i] = true;

                            nset q;
                            q.neigh[0] = false;
                            q.neigh[1] = false;
                            q.neigh[2] = false;
                            q.neigh[3] = false;
                            q.x = dx;
                            q.y = dy;
                            stack.push_back(q);
                            // After this point any access to s may crash

                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        stack.pop_back();
                    }
                }
            }
        }
    }
    grid_nclasses++;
    qDebug("Found %d equivalence classes", grid_nclasses);
    {
        qDebug("Rendering classes to image");
        FColor *colors = new FColor[grid_nclasses];
        for (int i = 0; i < grid_nclasses; i++) {
            colors[i] = randColor();
        }
        QRgb *dat = new QRgb[W * H];
        for (int x = 0; x < W; x++) {
            for (int y = 0; y < H; y++) {
                dat[y * W + x] =
                    colors[grid_points[x * H + y].region_class].rgba();
            }
        }

        QImage img((uchar *)dat, W, H, QImage::Format_ARGB32);
        emit produceImagePhase(img.copy(), QString("Grid completed"), nqueries,
                               false);

        delete[] dat;
        delete[] colors;
    }
}
static inline uint qHash(QPointF t) {
    return qHash(t.x()) + 176 * qHash(t.y());
}
static inline uint qHash(QLineF t) {
    return qHash(t.p1()) + 133 * qHash(t.p2());
}
static inline bool operator<(const QPoint &a, const QPoint &b) {
    if (a.x() == b.x())
        return a.y() < b.y();
    return a.x() < b.x();
}
static QPolygon reverse_poly(const QPolygon &p) {
    QPolygon r;
    for (int k = 0; k < p.size(); k++) {
        r.push_back(p[p.size() - 1 - k]);
    }
    return r;
}
static QVector<RenderPoint> reverse_rloop(const QVector<RenderPoint> &p) {
    QVector<RenderPoint> r;
    for (int k = 0; k < p.size(); k++) {
        r.push_back(p[p.size() - 1 - k]);
    }
    return r;
}

QPolygonF scale_shift_rp_loop(const QVector<RenderPoint> &p, double scale,
                              const QPointF &then_shift,
                              const QSize &grid_size) {
    QPolygonF q;
    for (const RenderPoint &r : p) {
        QPointF gc = point_to_grid_coord(r.coords, grid_size);
        q.push_back(gc * scale + then_shift);
    }
    return q;
}

static void draw_boundaries_to(QPaintDevice *target,
                               const QVector<Region> &region_list, double S,
                               const QSize &grid_size) {

    QPainter p(target);
    p.fillRect(0, 0, target->width(), target->height(), Qt::white);
    for (const Region &b : region_list) {
        QPainterPath path;
        // Note: QPainter+QtSVG doesn't support true SVG closed paths
        // so eventually an SVG generator will need to be made
        QPolygonF poly =
            scale_shift_rp_loop(b.exterior, S, QPointF(S, S), grid_size);
        path.moveTo(poly[0]);
        for (int i = 1; i < poly.size(); i++) {
            path.lineTo(poly[i]);
        }
        path.closeSubpath();
        for (const QVector<RenderPoint> &loop : b.interior) {
            poly = scale_shift_rp_loop(loop, S, QPointF(S, S), grid_size);
            path.moveTo(poly[0]);
            for (int i = 1; i < poly.size(); i++) {
                path.lineTo(poly[i]);
            }
            path.closeSubpath();
        }
        path.setFillRule(Qt::OddEvenFill);

        p.setPen(QPen(Qt::black, 1, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin));
        p.setBrush(QColor(randColor().rgba()));
        p.drawPath(path);
    }
}

static int signed_area(const QPolygon &loop) {
    int area = 0;
    for (int k0 = 0; k0 < loop.size(); k0++) {
        int k1 = (k0 + 1) % loop.size();
        area += (loop[k1].x() - loop[k0].x()) * (loop[k1].y() + loop[k1].y());
    }
    return area;
}

static bool signed_area_cmp(const QPolygon &left, const QPolygon &right) {
    int larea = signed_area(left);
    int rarea = signed_area(right);
    return larea < rarea;
}

/**
 *  Given a list of edge segments, return region polygons sorted in size
 *  from largest to smallest area
 */
QVector<QPolygon> size_sorted_loops_from_seg(const QVector<QLine> &segs) {

    /* Stitch together line segments into loops, via... eqvlclass
     * routine!
     */
    typedef struct {
        int ids[2];
    } Lns;
    typedef struct {
        bool visited[2];
        int id;
    } Ens;
    QMap<QPoint, Lns> links;

    for (int i = 0; i < segs.size(); i++) {
        for (int k = 0; k < 2; k++) {
            QPoint p = k ? segs[i].p1() : segs[i].p2();
            if (links.count(p)) {
                links[p].ids[1] = i;
            } else {
                // the -1s will be replaced anyway
                Lns e;
                e.ids[0] = i;
                e.ids[1] = -1;
                links[p] = e;
            }
        }
    }

    QVector<bool> reached(segs.size(), false);
    QVector<QPolygon> loops;
    for (int i = 0; i < segs.size(); i++) {
        if (reached[i])
            continue;
        reached[i] = true;

        QVector<Ens> stack;
        QVector<int> iloop;
        {
            Ens e;
            e.visited[0] = false;
            e.visited[1] = false;
            e.id = i;
            stack.push_back(e);
            iloop.push_back(i);
        }
        while (stack.size()) {
            Ens &test = stack.last();
            bool next_found = false;
            for (int k = 0; k < 2; k++) {
                if (test.visited[k])
                    continue;
                test.visited[k] = true;

                const QLine &l = segs[test.id];
                QPoint ep = k ? l.p1() : l.p2();
                Lns alt = links[ep];
                int aid = (alt.ids[0] == test.id) ? alt.ids[1] : alt.ids[0];

                if (aid < 0) {
                    qFatal("encountered single point in loop formation (%d %d)",
                           ep.x(), ep.y());
                    continue;
                }
                if (reached[aid])
                    continue;
                reached[aid] = true;

                Ens f;
                f.visited[0] = false;
                f.visited[1] = false;
                f.id = aid;
                stack.push_back(f);
                // after stack.push_back, test is invalid
                iloop.push_back(aid);
                next_found = true;
                break;
            }
            if (!next_found) {
                stack.pop_back();
            }
        }

        // Note len(iloop) >= 4
        QPolygon loop;
        if (segs[iloop[0]].p1() == segs[iloop[1]].p1() ||
            segs[iloop[0]].p1() == segs[iloop[1]].p2()) {
            // p1 shared, start with p2
            loop.push_back(segs[iloop[0]].p2());
        } else {
            // p2 shared, start with p1
            loop.push_back(segs[iloop[0]].p1());
        }

        for (int k : iloop) {
            QLine l = segs[k];
            if (l.p1() == loop.last()) {
                loop.push_back(l.p2());
            } else if (l.p2() == loop.last()) {
                loop.push_back(l.p1());
            } else {
                qFatal("Invariant failure");
            }
        }
        if (loop.first() != loop.last()) {
            qFatal("Invariant failure");
        }
        loop.pop_back();

        loops.push_back(loop);
    }

    /* Normalize loop orientations to be positive. Self crossing will
     * never happen for marching squares results */
    for (int i = 0; i < loops.size(); i++) {
        int area = signed_area(loops[i]);

        if (area > 0) {
            loops[i] = reverse_poly(loops[i]);
        }
    }

    /* Loop with greatest extent is the containing loop */
    qSort(loops.begin(), loops.end(), signed_area_cmp);

    return loops;
}

static QVector<RenderPoint>
remove_rloop_duplicates(const QVector<RenderPoint> &bound, double epsilon) {
    QVector<RenderPoint> altbound;
    altbound.push_back(bound[0]);
    for (int i = 0; i < bound.size() - 1; i++) {
        if ((bound[i].coords - altbound.last().coords).manhattanLength() <
            epsilon) {
            continue;
        }
        altbound.push_back(bound[i]);
    }
    if ((bound.last().coords - altbound.last().coords).manhattanLength() <
        epsilon) {

    } else if ((bound.last().coords - altbound[0].coords).manhattanLength() <
               epsilon) {

    } else {
        altbound.push_back(bound.last());
    }
    return altbound;
}

typedef struct {
    QPointF inside;
    QPointF outside;
    int index;
} CornerCandidate;

static QVector<CornerCandidate>
find_rloop_corners(const QVector<RenderPoint> &bound, bool is_exterior_loop,
                   double epsilon) {

    QVector<CornerCandidate> corner_cands;

    const int npts = bound.size();
    for (int i1 = 0; i1 < npts; i1++) {
        int i0 = (i1 + npts - 1) % npts;
        int i2 = (i1 + 1) % npts;
        int i3 = (i1 + 2) % npts;

        QPointF dir_pre = bound[i1].coords - bound[i0].coords;
        QPointF dir_post = bound[i2].coords - bound[i3].coords;
        double len_pre =
            std::sqrt(dir_pre.x() * dir_pre.x() + dir_pre.y() * dir_pre.y());
        double len_post = std::sqrt(dir_post.x() * dir_post.x() +
                                    dir_post.y() * dir_post.y());
        if (len_pre <= 0. || len_post <= 0.) {
            continue;
        }
        dir_pre /= len_pre;
        dir_post /= len_post;
        double dotprod = QPointF::dotProduct(dir_pre, dir_post);
        // -1 is aligned; +1 is ultra spiky
        if (dotprod < std::cos(5 * M_PI / 6) ||
            dotprod > std::cos(1 * M_PI / 6)) {
            continue;
        }

        double det = dir_pre.y() * dir_post.x() - dir_pre.x() * dir_post.y();

        QPointF src_pre = bound[i1].coords, src_post = bound[i2].coords;
        double t = (dir_post.x() * (src_post.y() - src_pre.y()) +
                    dir_post.y() * (src_pre.x() - src_post.x())) /
                   det;
        double s = (dir_pre.x() * (src_post.y() - src_pre.y()) +
                    dir_pre.y() * (src_pre.x() - src_post.x())) /
                   det;

        // Real space min offset limit for corner
        double cepsilon = epsilon;

        if (s <= cepsilon || t <= cepsilon) {
            // If not exactly on the angle, can find intersection
            // point behind k0
            continue;
        }
        QPointF qavg = 0.5 * (src_pre + src_post);
        QPointF qcor_pre = src_pre + t * dir_pre;
        QPointF qcor_post = src_post + s * dir_post;

        QPointF qcor = 0.5 * (qcor_pre + qcor_post);
        QPointF qjump = 2 * qcor - qavg;

        CornerCandidate cand;
        cand.index = i1;
        // TODO: instead, decide using curve orientation as an additional guide
        // to convexity/concavity. Would require early loop reversal
        if (is_exterior_loop) {
            cand.inside = qavg;
            cand.outside = qjump;
        } else {
            cand.inside = qjump;
            cand.outside = qavg;
        }
        corner_cands.push_back(cand);
    }
    return corner_cands;
}

static QVector<RenderPoint> rloop_merge(const QVector<RenderPoint> &base,
                                        const QVector<RenderPoint> &mod,
                                        const QVector<int> &indexes) {
    QVector<RenderPoint> bmod = base;
    for (int i = 0; i < mod.size(); i++) {
        bmod.insert(indexes[i] + i + 1, mod[i]);
    }
    return bmod;
}

class Links {
public:
    Links() {
        p[0] = dummy;
        p[1] = dummy;
        p[2] = dummy;
        p[3] = dummy;
    }
    void add(QPoint r) {
        for (int i = 0; i < 4; i++) {
            if (p[i] == r) {
                return;
            }
            if (p[i] == dummy) {
                p[i] = r;
                return;
            }
        }
        qFatal("Overflow on link");
    }
    int count() const {
        for (int i = 0; i < 4; i++) {
            if (p[i] == dummy) {
                return i;
            }
        }
        return 4;
    }
    void replace_first(QPoint o, QPoint n) {
        for (int i = 0; i < 4; i++) {
            if (p[i] == o) {
                p[i] = n;
                return;
            }
        }
        qFatal("Failed to find original");
    }

    QPoint p[4];
    static constexpr QPoint dummy = QPoint(1 << 30, 1 << 30);
};

static void link_add(QMap<QPoint, Links> &links, const QPoint &a,
                     const QPoint &b) {
    links[a].add(b);
    links[b].add(a);
}
/* The obvious way to do this doesn't handle negative numbers correctly */
static void partition_sumpoint(const QPoint &p, QPoint *n, QPoint *f) {
    int fpos = 2 * std::max(std::abs(p.x()), std::abs(p.y()));
    int x = p.x() + 2 * fpos;
    int y = p.y() + 2 * fpos;
    *n = QPoint(x / 2 - fpos, y / 2 - fpos);
    *f = p - *n;
}
/* again, C standard isn't math standard */
static bool point_is_oo(const QPoint &p) {
    return (p.x() & 0x1) && (p.y() & 0x1);
}

void VectorTracer::computeEdges() {
    region_list.clear();
    int W = grid_size.width(), H = grid_size.height();
    Navigator *nav = Navigator::create(view_data, view_data.navigator);

    if (0) {
        for (int cls = 0; cls < grid_nclasses; cls++) {
            qDebug("Computing region boundaries %d", cls);
            // Compute a hull path/curve for the points, sampling new points as
            // needed. Right now, just a hull, with convex openings.
            int xmin = W, xmax = -1, ymin = H, ymax = -1;
            int lx = -1, ly = -1;
            for (int x = 0; x < W; x++) {
                for (int y = 0; y < H; y++) {
                    int idx = x * H + y;
                    if (grid_points[idx].region_class == cls) {
                        xmin = std::min(x, xmin);
                        xmax = std::max(x, xmax);
                        ymin = std::min(y, ymin);
                        ymax = std::max(y, ymax);
                        lx = x;
                        ly = y;
                    }
                }
            }

            QVector<QLine> segs;
            for (int x = xmin - 1; x < xmax + 1; x++) {
                for (int y = ymin - 1; y < ymax + 1; y++) {
                    int code = 0;
                    const int dxl[4] = {0, 1, 1, 0};
                    const int dyl[4] = {0, 0, 1, 1};
                    for (int i = 0; i < 4; i++) {
                        int nx = x + dxl[i], ny = y + dyl[i];
                        if (nx < 0 || nx >= W || ny < 0 || ny >= H)
                            continue;
                        if (grid_points[nx * H + ny].region_class == cls)
                            code |= 1 << i;
                    }

                    switch (code) {
                    case 0xF:
                    case 0x0:
                        // 0000=1111 All or no corners are in class
                        break;
                    case 0xE:
                    case 0x1:
                        // 0001=1110 Corner (x,y)
                        segs.push_back(
                            QLine(2 * x + 1, 2 * y, 2 * x, 2 * y + 1));
                        break;
                    case 0xD:
                    case 0x2:
                        // 0010=1101 Corner (x+1,y)
                        segs.push_back(
                            QLine(2 * x + 1, 2 * y, 2 * x + 2, 2 * y + 1));
                        break;
                    case 0xB:
                    case 0x4:
                        // 0100=1011 Corner (x+1,y+1)
                        segs.push_back(
                            QLine(2 * x + 1, 2 * y + 2, 2 * x + 2, 2 * y + 1));
                        break;
                    case 0x8:
                    case 0x7:
                        // 0111=1000 Corner (x,y+1)
                        segs.push_back(
                            QLine(2 * x + 1, 2 * y + 2, 2 * x, 2 * y + 1));
                        break;

                    case 0xC:
                    case 0x3:
                        // 0011=1100 Line (along x)
                        segs.push_back(
                            QLine(2 * x, 2 * y + 1, 2 * x + 2, 2 * y + 1));
                        break;
                    case 0x9:
                    case 0x6:
                        // 0110=1001 Line (along y)
                        segs.push_back(
                            QLine(2 * x + 1, 2 * y, 2 * x + 1, 2 * y + 2));
                        break;

                    case 0xA:
                    case 0x5:
                        // 0101=1010 Saddle; must test center point to decide
                        {
                            QPointF pt =
                                0.5 * (grid_points[x * H + y].coords +
                                       grid_points[(x + 1) * H + y + 1].coords);
                            RenderPoint rp = queryPoint(pt, nav);
                            const RenderPoint &class_rep =
                                (code == 0x5) ? grid_points[x * H + y]
                                              : grid_points[x * H + y + 1];

                            bool center_inclass = typematch(class_rep, rp);
                            if ((code == 0x5) ^ center_inclass) {
                                // Diagonals along (x,y+1)<->(x+1,y)
                                segs.push_back(
                                    QLine(2 * x + 1, 2 * y, 2 * x, 2 * y + 1));
                                segs.push_back(QLine(2 * x + 1, 2 * y + 2,
                                                     2 * x + 2, 2 * y + 1));
                            } else {
                                // Diagonals along (x,y)<->(x+1,y+1)
                                segs.push_back(QLine(2 * x + 1, 2 * y,
                                                     2 * x + 2, 2 * y + 1));
                                segs.push_back(QLine(2 * x + 1, 2 * y + 2,
                                                     2 * x, 2 * y + 1));
                            }
                        }

                        break;

                    default:
                        qFatal("Shouldn't happen");
                    }
                }
            }

            qDebug("Constructing region boundaries %d of %d: %d segments", cls,
                   grid_nclasses, segs.length());

            QVector<QPolygon> loops = size_sorted_loops_from_seg(segs);

            qDebug("Refining region boundaries %d: %d loops", cls,
                   loops.size());

            /* Refine all loops to stay barely within class. This _cannot_ be
             * shared between classes, as boundaries aren't perfect cliffs and
             * there may be third types in between. */
            const RenderPoint &class_representative = grid_points[lx * H + ly];

            Region region;
            region.class_no = cls;
            region.xmin = xmin;
            region.xmax = xmax;
            region.ymin = ymin;
            region.ymax = ymax;

            for (int i = 0; i < loops.size(); i++) {
                QVector<RenderPoint> qlp;
                // Loop refinement
                for (QPoint loc : loops[i]) {
                    // Binary search for the class point nearest the outside
                    QPoint ilow = loc / 2;
                    QPoint ihigh = loc - loc / 2;
                    bool low_inclass = true;
                    if (ilow.x() < 0 || ilow.x() >= W || ilow.y() < 0 ||
                        ilow.y() >= H) {
                        low_inclass = false;
                    } else if (grid_points[ilow.x() * H + ilow.y()]
                                   .region_class != cls) {
                        low_inclass = false;
                    }
                    QPoint ip_in = low_inclass ? ilow : ihigh;
                    QPoint ip_out = low_inclass ? ihigh : ilow;
                    RenderPoint in_point = getPoint(ip_in);
                    RenderPoint out_point = getPoint(ip_out);

                    RenderPoint in_limit, out_limit;
                    bracketEdge(in_point, out_point, &in_limit, &out_limit,
                                nav);
                    in_limit.ideal_color =
                        calculateBoundaryColor(in_limit, out_limit);

                    if (qlp.size() && (in_limit.coords - qlp.last().coords)
                                              .manhattanLength() < 1e-15) {
                        // In case bisection search, both times, reaches the
                        // grid point.
                        continue;
                    }
                    qlp.push_back(in_limit);
                }

                /* Loop corner insertion.
                 */
                double epsilon = (1. / std::max(W, H)) * 1e-5;
                qlp = remove_rloop_duplicates(qlp, epsilon);

                QVector<CornerCandidate> corner_cands =
                    find_rloop_corners(qlp, i == 0, epsilon);

                QVector<RenderPoint> pts;
                QVector<int> idxs;
                for (const CornerCandidate &corner : corner_cands) {
                    RenderPoint cin = queryPoint(corner.inside, nav);
                    RenderPoint cout = queryPoint(corner.outside, nav);
                    bool in_icls = typematch(cin, class_representative);
                    bool out_icls = typematch(cout, class_representative);

                    if (!in_icls || out_icls) {
                        continue;
                    }

                    RenderPoint plim_in, plim_out;
                    bracketEdge(cin, cout, &plim_in, &plim_out, nav);
                    plim_in.ideal_color =
                        calculateBoundaryColor(plim_in, plim_out);

                    pts.push_back(plim_in);
                    idxs.push_back(corner.index);
                }

                qlp = rloop_merge(qlp, pts, idxs);

                /* Place on boundary */
                if (i == 0) {
                    // Exterior and interior regions must have opposite
                    // orientation for most odd-even fill implementations to
                    // leave holes
                    region.exterior = reverse_rloop(qlp);
                } else {
                    region.interior.push_back(qlp);
                }
            }
            region_list.push_back(region);
        }
    }
    if (1) {
        /* Graph-style partitioning */
        QVector<QPoint> eooe_pts;
        QVector<QPoint> oo_multi_pts;
        QVector<QPoint> oo_pair_pts;
        QMap<QPoint, Links> link_map;
        for (int x = -1; x < W; x++) {
            for (int y = -1; y < H; y++) {
                // Iterate over loop
                const QPoint corners[4] = {
                    QPoint(x, y),
                    QPoint(x + 1, y),
                    QPoint(x + 1, y + 1),
                    QPoint(x, y + 1),
                };
                int ncls[4];
                for (int i = 0; i < 4; i++) {
                    if (corners[i].x() < 0 || corners[i].x() >= W ||
                        corners[i].y() < 0 || corners[i].y() >= H) {
                        ncls[i] = -1;
                    } else {
                        ncls[i] =
                            grid_points[corners[i].x() * H + corners[i].y()]
                                .region_class;
                    }
                }
                int cmin = std::min(std::min(ncls[0], ncls[1]),
                                    std::min(ncls[2], ncls[3]));
                int cmax = std::max(std::max(ncls[0], ncls[1]),
                                    std::max(ncls[2], ncls[3]));
                if (cmin == cmax) {
                    // All the same, do nothing
                    continue;
                }
                int dmin = 1 << 30, dmax = -(1 << 30);
                for (int i = 0; i < 4; i++) {
                    if (ncls[i] != cmin && ncls[i] != cmax) {
                        dmin = std::min(ncls[i], dmin);
                        dmax = std::max(ncls[i], dmax);
                    }
                }
                // ncls[0], splits[0], ncls[1], splits[1], ncls[2], ...
                QPoint splits[4] = {
                    corners[0] + corners[1], corners[1] + corners[2],
                    corners[2] + corners[3], corners[3] + corners[0]};
                QPoint center(2 * x + 1, 2 * y + 1);

                if (dmin > dmax) {
                    // Two different types, marchcube
                    bool islow[4] = {(ncls[0] == cmin), (ncls[1] == cmin),
                                     (ncls[2] == cmin), (ncls[3] == cmin)};
                    int nlow = islow[0] + islow[1] + islow[2] + islow[3];
                    if (nlow == 1) {
                        // Corner configuration, 1 low
                        for (int i = 0; i < 4; i++) {
                            if (islow[i]) {
                                link_add(link_map, splits[(i + 3) % 4],
                                         splits[i]);
                            }
                        }
                    } else if (nlow == 3) {
                        // Corner configuration, 1 high
                        for (int i = 0; i < 4; i++) {
                            if (!islow[i]) {
                                link_add(link_map, splits[(i + 3) % 4],
                                         splits[i]);
                            }
                        }
                    } else {
                        // Split configuration, nlow == 2
                        if (islow[0] == islow[2] || islow[1] == islow[3]) {
                            // Saddle
                            QPointF pt =
                                0.5 * (grid_points[x * H + y].coords +
                                       grid_points[(x + 1) * H + y + 1].coords);
                            RenderPoint rp = queryPoint(pt, nav);
                            const RenderPoint &rep0 =
                                grid_points[corners[0].x() * H +
                                            corners[0].y()];

                            bool center_inclass = typematch(rep0, rp);
                            // 0XaO1     0XaO1
                            //  dXb   or  dOb
                            // 3OcX2     3OcX2
                            if (center_inclass) {
                                link_add(link_map, splits[0], splits[1]);
                                link_add(link_map, splits[2], splits[3]);
                            } else {
                                link_add(link_map, splits[3], splits[0]);
                                link_add(link_map, splits[1], splits[2]);
                            }
                        } else if (islow[0] == islow[1]) {
                            // Axis aligned
                            link_add(link_map, splits[1], splits[3]);
                        } else if (islow[0] == islow[3]) {
                            // Axis aligned
                            link_add(link_map, splits[0], splits[2]);
                        }
                    }
                } else if (dmin == dmax) {
                    // three distinct types, T shape.
                    for (int i = 0; i < 4; i++) {
                        int j = (i + 1) % 4;
                        if (ncls[i] != ncls[j]) {
                            link_add(link_map, center, splits[i]);
                        }
                    }
                } else {
                    // four distinct types, cross shape
                    for (int i = 0; i < 4; i++) {
                        link_add(link_map, center, splits[i]);
                    }
                }
            }
        }
        for (const QPoint &p : link_map.keys()) {
            if (point_is_oo(p)) {
                oo_multi_pts.append(p);
            } else {
                eooe_pts.append(p);
            }
        }

        // Progressive boundary refinement: E/O points
        QMap<QPoint, QPair<RenderPoint, RenderPoint>> eo_refinement;
        for (QPoint &p : eooe_pts) {
            QPoint i0, i1;
            partition_sumpoint(p, &i0, &i1);

            RenderPoint p0 = getPoint(i0);
            RenderPoint p1 = getPoint(i1);

            QPair<RenderPoint, RenderPoint> &l = eo_refinement[p];
            RenderPoint &l0 = l.first, &l1 = l.second;
            bracketEdge(p0, p1, &l0, &l1, nav);

            l0.region_class = p0.region_class;
            l0.ideal_color = calculateBoundaryColor(l0, l1);

            l1.region_class = p1.region_class;
            l1.ideal_color = calculateBoundaryColor(l1, l0);
        }

        QMap<QPoint, QPair<RenderPoint, RenderPoint>> oo_corners;
        QMap<QPoint, QPair<QPoint, QPoint>> oo_insertions;
        // Iterate over EO/OE links
        for (const QPoint &qb : eooe_pts) {
            const Links &lb = link_map[qb];
            if (lb.count() != 2) {
                continue;
            }
            for (int j = 0; j < 2; j++) {
                const QPoint &qc = lb.p[j];
                const Links &lc = link_map[qc];
                if (lc.count() != 2) {
                    continue;
                }
                QPoint qa = lb.p[1 - j];
                QPoint qd = (lc.p[0] == qb) ? lc.p[1] : lc.p[0];
                if (point_is_oo(qa) || point_is_oo(qb) || point_is_oo(qc) ||
                    point_is_oo(qd)) {
                    continue;
                }

                // either along diagonal, or along line. beware saddle points.
                QPoint target;
                if (point_is_oo(qb + qc)) {
                    if (qb.x() & 0x1) {
                        target = QPoint(qb.x(), qc.y());
                    } else {
                        target = QPoint(qc.x(), qb.y());
                    }
                } else {
                    target = (qb + qc) / 2;
                }

                if (oo_corners.contains(target)) {
                    // already used
                    continue;
                }

                QPointF pa = (eo_refinement[qa].first.coords +
                              eo_refinement[qa].second.coords) /
                             2;
                QPointF pb = (eo_refinement[qb].first.coords +
                              eo_refinement[qb].second.coords) /
                             2;
                QPointF pc = (eo_refinement[qc].first.coords +
                              eo_refinement[qc].second.coords) /
                             2;
                QPointF pd = (eo_refinement[qd].first.coords +
                              eo_refinement[qd].second.coords) /
                             2;

                // A1,A2,B1,B2 ...
                QPointF dir_pre = pb - pa;
                QPointF dir_post = pc - pd;
                double len_pre = std::sqrt(dir_pre.x() * dir_pre.x() +
                                           dir_pre.y() * dir_pre.y());
                double len_post = std::sqrt(dir_post.x() * dir_post.x() +
                                            dir_post.y() * dir_post.y());
                if (len_pre <= 0. || len_post <= 0.) {
                    continue;
                }
                dir_pre /= len_pre;
                dir_post /= len_post;
                double dotprod = QPointF::dotProduct(dir_pre, dir_post);
                // -1 is aligned; +1 is ultra spiky
                if (dotprod < std::cos(5 * M_PI / 6) ||
                    dotprod > std::cos(1 * M_PI / 6)) {
                    continue;
                }

                double det =
                    dir_pre.y() * dir_post.x() - dir_pre.x() * dir_post.y();

                QPointF src_pre = pb, src_post = pc;
                double t = (dir_post.x() * (src_post.y() - src_pre.y()) +
                            dir_post.y() * (src_pre.x() - src_post.x())) /
                           det;
                double s = (dir_pre.x() * (src_post.y() - src_pre.y()) +
                            dir_pre.y() * (src_pre.x() - src_post.x())) /
                           det;

                // Real space min offset limit for corner
                double epsilon = (1. / std::max(W, H)) * 1e-5;
                if (s <= epsilon || t <= epsilon) {
                    // If not exactly on the angle, can find intersection
                    // point behind
                    continue;
                }

                // A series of 4 EO/OE points

                QPointF qavg = 0.5 * (src_pre + src_post);
                QPointF qcor_pre = src_pre + t * dir_pre;
                QPointF qcor_post = src_post + s * dir_post;

                QPointF qcor = 0.5 * (qcor_pre + qcor_post);
                QPointF qjump = 2 * qcor - qavg;

                // Do binary search between qcor and qjump
                // but only if qcor/qjump match classes as necessary

                RenderPoint pavg = queryPoint(qavg, nav);
                RenderPoint pjump = queryPoint(qjump, nav);

                if (typematch(pavg, eo_refinement[qb].first)) {
                    pavg.region_class = eo_refinement[qb].first.region_class;
                } else if (typematch(pavg, eo_refinement[qb].second)) {
                    pavg.region_class = eo_refinement[qb].second.region_class;
                } else {
                    // We jumped out of class, expect fine detail or noise,
                    // so refinement is futile
                    continue;
                }
                if (typematch(pjump, eo_refinement[qc].first)) {
                    pjump.region_class = eo_refinement[qc].first.region_class;
                } else if (typematch(pjump, eo_refinement[qc].second)) {
                    pjump.region_class = eo_refinement[qc].second.region_class;
                } else {
                    continue;
                }

                if (pavg.region_class == pjump.region_class) {
                    // Both points in the same region, no obvious corner to
                    // refine
                    continue;
                }

                RenderPoint l0, l1;
                bracketEdge(pavg, pjump, &l0, &l1, nav);

                l0.region_class = pavg.region_class;
                l0.ideal_color = calculateBoundaryColor(l0, l1);

                l1.region_class = pjump.region_class;
                l1.ideal_color = calculateBoundaryColor(l1, l0);

                // Note: there is, unlike with the EO/OE crossings, no
                // specific ordering to this sort of intersection.
                oo_corners[target] = QPair<RenderPoint, RenderPoint>(l0, l1);
                oo_insertions[target] = QPair<QPoint, QPoint>(qb, qc);
            }
        }
        for (const QPoint &inserted : oo_insertions.keys()) {
            const QPair<QPoint, QPoint> &vls = oo_insertions[inserted];
            qDebug("interpolating (%d,%d) to (%d,%d) with (%d,%d)",
                   vls.first.x(), vls.first.y(), vls.second.x(), vls.second.y(),
                   inserted.x(), inserted.y());
            Links &l0 = link_map[vls.first];
            Links &l1 = link_map[vls.second];
            l0.replace_first(vls.second, inserted);
            l1.replace_first(vls.first, inserted);
            link_map[inserted].add(vls.first);
            link_map[inserted].add(vls.second);
        }

        // TODO: FIXME, create a suitable subdivision algorithm for 3- and 4-
        // joints! (for instance: based on the average entry point, produce a
        // triangle whose three corners are very likely inside the constituent
        // regions; if so, perform contraction (try to move each corner inward
        // by 50% to the average point of the neighboring diagonal, as long as
        // the color remains the same.) as long as possible or until the
        // polytope has shrunk enough.

        // Every loop contains at least one E/O or O/E point, and those
        // have exactly two parents.
        QVector<QSet<QPoint>> class_points(grid_nclasses);
        for (QPoint p : eooe_pts) {
            QPoint n, f;
            partition_sumpoint(p, &n, &f);

            if (n.x() >= 0 && n.x() < W && n.y() >= 0 && n.y() < H) {
                class_points[grid_points[n.x() * H + n.y()].region_class]
                    .insert(p);
            }
            if (f.x() >= 0 && f.x() < W && f.y() >= 0 && f.y() < H) {
                class_points[grid_points[f.x() * H + f.y()].region_class]
                    .insert(p);
            }
        }

        QVector<Region> alt_reglist;
        for (int cls = 0; cls < grid_nclasses; cls++) {
            // Construct loops from class point seeds
            QSet<QPoint> &pts = class_points[cls];
            QVector<QVector<QPoint>> loops;
            while (pts.size()) {
                QPoint start = *pts.begin();
                QVector<QPoint> loop;
                // TODO: a single loop may hit a 4-junction with 3 classes
                // twice; must fix!

                // The solution would be an explicit winding direction heuristic
                // (i.e, we move along EO/OE nodes so that our class is
                // always to the left of our motion, picking up O/E nodes as
                // necessary.)

                QSet<QPoint> oo_seen;
                loop.append(start);
                pts.remove(start);
                while (true) {
                    Links l = link_map[start];
                    bool found = false;
                    for (int i = 0; i < l.count(); i++) {
                        QPoint a = l.p[i];
                        if (a == start) {
                            // no backtrack
                            continue;
                        }
                        // even => odd or even; odd => even
                        bool is_odd = point_is_oo(a);
                        if (is_odd) {
                            if (!oo_seen.contains(a)) {
                                start = a;
                                loop.append(start);
                                oo_seen.insert(start);
                                found = true;
                                break;
                            }
                        } else {
                            if (pts.contains(a)) {
                                // continue to an even point oif it borders
                                // class
                                start = a;
                                loop.append(start);
                                pts.remove(start);
                                found = true;
                                break;
                            }
                        }
                    }
                    if (!found) {
                        break;
                    }
                }
                loops.append(loop);
            }
            for (int i = 0; i < loops.size(); i++) {
                int area = signed_area(loops[i]);

                if (area > 0) {
                    loops[i] = reverse_poly(loops[i]);
                }
            }

            /* Loop with greatest extent is the containing loop */
            qSort(loops.begin(), loops.end(), signed_area_cmp);

            Region reg;
            reg.class_no = cls;
            reg.xmin = W;
            reg.xmax = -1;
            reg.ymin = H;
            reg.ymax = -1;

            // Interior/exterior loops -- make pairs. Then *skip* O/O points,
            // for now
            for (int i = 0; i < loops.size(); i++) {
                QVector<RenderPoint> rp;
                for (QPoint p : loops[i]) {
                    if (point_is_oo(p)) {
                        if (oo_corners.contains(p)) {
                            // Corner
                            const QPair<RenderPoint, RenderPoint> &cp =
                                oo_corners[p];
                            if (cp.first.region_class == cls) {
                                rp.append(cp.first);
                            } else {
                                rp.append(cp.second);
                            }
                        } else {
                            // This is a joint!
                            qDebug("Skipped joint at (%d,%d)", p.x(), p.y());

                            // TODO: get these working!
                            continue;
                        }
                    } else {
                        // Line portion
                        const QPair<RenderPoint, RenderPoint> &ep =
                            eo_refinement[p];
                        if (ep.first.region_class == cls) {
                            rp.append(ep.first);
                        } else {
                            rp.append(ep.second);
                        }
                    }
                }
                if (i == 0) {
                    reg.exterior = reverse_rloop(rp);
                } else {
                    reg.interior.append(rp);
                }
            }

            alt_reglist.append(reg);
        }

        // Set bounds
        for (int x = 0; x < W; x++) {
            for (int y = 0; y < H; y++) {
                int idx = x * H + y;
                int cls = grid_points[idx].region_class;
                alt_reglist[cls].xmin = std::min(x, alt_reglist[cls].xmin);
                alt_reglist[cls].xmax = std::max(x, alt_reglist[cls].xmax);
                alt_reglist[cls].ymin = std::min(y, alt_reglist[cls].ymin);
                alt_reglist[cls].ymax = std::max(y, alt_reglist[cls].ymax);
            }
        }

        if (0) {
            double S = 6;
            QImage disp(2 * S * W + 2 * S, 2 * S * H + 2 * S,
                        QImage::Format_ARGB32);
            disp.fill(0xFFFFFFFF);
            QPainter ptr(&disp);
            ptr.setPen(Qt::black);
            for (const QPoint &p : link_map.keys()) {
                const Links &l = link_map[p];
                for (int i = 0; i < 4; i++) {
                    if (l.p[i] != Links::dummy) {
                        ptr.drawLine(S * p + QPoint(3 * S / 2, 3 * S / 2),
                                     S * l.p[i] + QPoint(3 * S / 2, 3 * S / 2));
                    }
                }
            }
            ptr.setPen(Qt::red);
            for (const QPoint &p : link_map.keys()) {
                ptr.drawEllipse(QPointF(S * p + QPoint(3 * S / 2, 3 * S / 2)),
                                S / 3., S / 3.);
            }
            for (int i = 0; i < grid_nclasses; i++) {
                QColor cc = QColor(randColor().rgba());
                ptr.setPen(
                    QPen(cc, 1, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin));
                ptr.setBrush(cc);
                for (int x = 0; x < W; x++) {
                    for (int y = 0; y < H; y++) {
                        if (grid_points[x * H + y].region_class == i) {
                            QPoint p(x, y);
                            ptr.drawEllipse(
                                QPointF(2 * S * p +
                                        QPoint(3 * S / 2, 3 * S / 2)),
                                S / 4., S / 4.);
                        }
                    }
                }
            }
            emit produceImagePhase(disp, QString("Boundaries completed"),
                                   nqueries, false);
        }

        if (0) {
            // Draw them lines!
            int S = std::max(1, 1024 / std::max(W + 1, H + 1));
            QImage lineImage((W + 1) * S, (H + 1) * S, QImage::Format_ARGB32);
            draw_boundaries_to(&lineImage, alt_reglist, S, grid_size);
            emit produceImagePhase(lineImage, QString("Boundaries completed"),
                                   nqueries, false);
        }
        region_list = alt_reglist;
    }

    delete nav;
    qDebug("Number of regions: %d", region_list.size());

    // Draw them lines!
    int S = std::max(1, 1024 / std::max(W + 1, H + 1));
    QImage lineImage((W + 1) * S, (H + 1) * S, QImage::Format_ARGB32);
    draw_boundaries_to(&lineImage, region_list, S, grid_size);
    emit produceImagePhase(lineImage, QString("Boundaries completed"), nqueries,
                           false);
}

static bool crease_check(const Intersection &a0, const Intersection &a1,
                         const QPointF &c0, const QPointF &c1,
                         const ViewData &view_data, double cos_alpha,
                         double min_jump, bool *is_jump) {
    G4ThreeVector n0(a0.normal), n1(a1.normal);

    if (n0.mag() < 1e-10 || n1.mag() < 1e-10) {
        // empty normals
        return false;
    }
    if (std::abs(n0.mag() - 1.0) > 0.01 || std::abs(n1.mag() - 1.0) > 0.01) {
        qFatal("Non-normal nontrivial normals %.18f %.18f", n0.mag(), n1.mag());
    }

    if (n0.dot(n1) < cos_alpha) {
        if (is_jump)
            *is_jump = false;
        return true;
    }

    // Real space motion
    const G4ThreeVector &dy =
        view_data.scale * view_data.orientation.rowY() * (c1.y() - c0.y());
    const G4ThreeVector &dx =
        view_data.scale * view_data.orientation.rowZ() * (c1.x() - c0.x());
    const G4ThreeVector &F = view_data.orientation.rowX();
    G4ThreeVector disp = dx + dy;
    G4ThreeVector N = 0.5 * (n0 + n1);
    if (std::abs(N.dot(F)) < 1e-10)
        return false;

    double dexp = -N.dot(disp) / N.dot(F);
    double dact = a1.dist - a0.dist;
    double acceptable_error = min_jump * view_data.scale;

    // Jump discontinuity
    if (std::abs(dexp - dact) > acceptable_error) {
        if (is_jump)
            *is_jump = true;
        return true;
    }
    return false;
}
int VectorTracer::crease_depth(const RenderPoint &q0, const RenderPoint &q1,
                               double cos_alpha, double min_jump,
                               bool *is_jump) {
    // We assume no significant mismatches
    if (!q0.ray.intersections || !q1.ray.intersections)
        return -1;

    // Two clipping planes will always produce a crease
    if (q0.ray.N >= 1 && q1.ray.N >= 1 && q0.ray.front_clipped &&
        q1.ray.front_clipped &&
        q0.ray.intersections[0].normal != q1.ray.intersections[0].normal) {
        return 0;
    }

    for (int i = 0; i < std::min(q0.ray.N, q1.ray.N); i++) {
        const Intersection &a0 = q0.ray.intersections[i];
        const Intersection &a1 = q1.ray.intersections[i];
        if (crease_check(a0, a1, q0.coords, q1.coords, view_data, cos_alpha,
                         min_jump, is_jump)) {
            return i;
        }
    }

    return -1;
}

void VectorTracer::computeCreases() {
    qDebug("Computing creases and boundary color details");

    Navigator *nav = Navigator::create(view_data, view_data.navigator);
    const int W = grid_size.width(), H = grid_size.height();
    // Basically, there are creases inside each boundary

    // We also should determine the intensity of separating lines
    // as the point of separation may be set back a bit; in that case,
    // lines may vary in color as they progress

    // For now, we mark as creases any sharp (>30 degree jump) discontinuities
    // in angle, or >view_radius*1e-6 jumps in position.
    // Will need to binary-subdivide on each candidate line to verify this,
    // where candidates use looser (15 degree, 1e-3) conditions
    crease_edge_map.clear();
    int *flood_fill_data = new int[W * H];
    int *crease_edge_count = new int[W * H];
    for (int x = 0; x < W; x++) {
        for (int y = 0; y < H; y++) {
            crease_edge_count[x * H + y] = 0;
            flood_fill_data[x * H + y] = 0;
        }
    }
    const double cos_alpha = std::cos(M_PI / 12);
    const double min_jump = 1e-3;

    const int nsubclasses = W * H;
    for (Region &region : region_list) {
        QVector<QPoint> pt_list;
        for (int x = region.xmin; x <= region.xmax; x++) {
            for (int y = region.ymin; y <= region.ymax; y++) {
                for (int k = 0; k < 2; k++) {
                    // k=0, horiz; k=1, vert
                    QPoint line_marker(2 * x + (1 - k), 2 * y + k);
                    // A sharp yes/no on is_transition.
                    QPoint p0 = line_marker / 2;
                    if (p0.x() < 0 || p0.y() < 0 || p0.x() >= W ||
                        p0.y() >= H) {
                        continue;
                    }

                    const RenderPoint &q0 = grid_points[p0.x() * H + p0.y()];
                    if (q0.region_class != region.class_no) {
                        continue;
                    }
                    pt_list.push_back(p0);

                    QPoint p1 = line_marker - line_marker / 2;
                    if (p1.x() < 0 || p1.y() < 0 || p1.x() >= W ||
                        p1.y() >= H) {
                        continue;
                    }

                    const RenderPoint &q1 = grid_points[p1.x() * H + p1.y()];
                    if (q1.region_class != region.class_no) {
                        continue;
                    }

                    bool has_crease =
                        crease_depth(q0, q1, cos_alpha, min_jump) >= 0;

                    crease_edge_map[line_marker] = has_crease;
                    crease_edge_count[p0.x() * H + p0.y()]++;
                    crease_edge_count[p1.x() * H + p1.y()]++;
                }
            }
        }

        // Iterative region floodfill, works even with one subregion
        const QPoint offsets[4] = {QPoint(0, 1), QPoint(1, 0), QPoint(0, -1),
                                   QPoint(-1, 0)};
        QSet<QPoint> pt_set = QSet<QPoint>::fromList(pt_list.toList());
        int subreg_no = -1;
        QVector<QPoint> altered_list = pt_list;
        QVector<QPoint> subreg_representative_coords;
        while (pt_set.size()) {
            subreg_no += 1;
            pt_list = pt_set.toList().toVector();
            for (QPoint p : altered_list) {
                flood_fill_data[p.x() * H + p.y()] = 0;
            }
            altered_list.clear();

            QPoint p0 = pt_list[randint(pt_list.size())];
            grid_points[p0.x() * H + p0.y()].subregion_class = subreg_no;
            subreg_representative_coords.append(p0);
            flood_fill_data[p0.x() * H + p0.y()] = 1;
            pt_set.remove(p0);
            altered_list.push_back(p0);

            QVector<QPoint> queue;
            for (int i = 0; i < 4; i++) {
                queue.push_back(p0 + offsets[i]);
            }

            // Number of new positive candidates
            int nposadd = 1;
            for (int step = 2; nposadd && step < W * H; step++) {
                nposadd = 0;
                QSet<QPoint> nqueue;
                while (queue.size()) {
                    QPoint p = queue.last();
                    queue.pop_back();

                    if (p.x() < 0 || p.y() < 0 || p.x() >= W || p.y() >= H)
                        continue;
                    if (grid_points[p.x() * H + p.y()].region_class !=
                        region.class_no)
                        continue;
                    // No overwriting the original pair
                    if (flood_fill_data[p.x() * H + p.y()] != 0)
                        continue;
                    if (!pt_set.contains(p))
                        continue;

                    int nneg = 0;
                    int npos = 0;
                    bool active = false;
                    for (int k = 0; k < 4; k++) {
                        QPoint q = p + offsets[k];
                        QPoint line = p + q;
                        if (!crease_edge_map.count(line)) {
                            // Not a valid line
                            continue;
                        }
                        int qv = flood_fill_data[q.x() * H + q.y()];
                        // Not yet visited, so include it
                        bool crease_sep = crease_edge_map[line];
                        if (qv == 0) {
                            continue;
                        }
                        if (qv < 0 && crease_edge_count[q.x() * H + q.y()] &&
                            crease_edge_count[p.x() * H + p.y()]) {
                            // Don't accept negatives from crease edge point
                            // if self has one
                            continue;
                        }
                        // Only positives can cross a crease, becoming negative
                        // in the process.
                        if (crease_sep) {
                            if (qv < 0)
                                continue;
                            nneg += 1;
                        } else {
                            if (qv < 0) {
                                nneg += 1;
                            } else {
                                npos += 1;
                            }
                        }
                        active = true;
                    }
                    if (!active)
                        continue;

                    bool negative = nneg > npos;

                    flood_fill_data[p.x() * H + p.y()] =
                        negative ? -step : step;
                    altered_list.push_back(p);

                    // Only add neighbors if this cell changes
                    for (int k = 0; k < 4; k++) {
                        nqueue.insert(p + offsets[k]);
                        if (!negative)
                            nposadd++;
                    }

                    if (!negative) {
                        grid_points[p.x() * H + p.y()].subregion_class =
                            subreg_no;
                        pt_set.remove(p);
                    }
                }
                queue = nqueue.toList().toVector();
            }
        }
        subreg_no++;

        for (int subcls = 0; subcls < subreg_no; subcls++) {
            // Marching squares
            QVector<QLine> segs;
            for (int x = region.xmin - 1; x < region.xmax + 1; x++) {
                for (int y = region.ymin - 1; y < region.ymax + 1; y++) {
                    int code = 0;
                    const int dxl[4] = {0, 1, 1, 0};
                    const int dyl[4] = {0, 0, 1, 1};
                    for (int i = 0; i < 4; i++) {
                        int nx = x + dxl[i], ny = y + dyl[i];
                        if (nx < 0 || nx >= W || ny < 0 || ny >= H)
                            continue;
                        if (grid_points[nx * H + ny].region_class !=
                            region.class_no) {
                            continue;
                        }
                        if (grid_points[nx * H + ny].subregion_class == subcls)
                            code |= 1 << i;
                    }

                    switch (code) {
                    case 0xF:
                    case 0x0:
                        // 0000=1111 All or no corners are in class
                        break;
                    case 0xE:
                    case 0x1:
                        // 0001=1100 Corner (x,y)
                        segs.push_back(
                            QLine(2 * x + 1, 2 * y, 2 * x, 2 * y + 1));
                        break;
                    case 0xD:
                    case 0x2:
                        // 0010=1101 Corner (x+1,y)
                        segs.push_back(
                            QLine(2 * x + 1, 2 * y, 2 * x + 2, 2 * y + 1));
                        break;
                    case 0xB:
                    case 0x4:
                        // 0100=1011 Corner (x+1,y+1)
                        segs.push_back(
                            QLine(2 * x + 1, 2 * y + 2, 2 * x + 2, 2 * y + 1));
                        break;
                    case 0x8:
                    case 0x7:
                        // 0111=1000 Corner (x,y+1)
                        segs.push_back(
                            QLine(2 * x + 1, 2 * y + 2, 2 * x, 2 * y + 1));
                        break;

                    case 0xC:
                    case 0x3:
                        // 0011=1100 Line (along x)
                        segs.push_back(
                            QLine(2 * x, 2 * y + 1, 2 * x + 2, 2 * y + 1));
                        break;
                    case 0x9:
                    case 0x6:
                        // 0110=1001 Line (along y)
                        segs.push_back(
                            QLine(2 * x + 1, 2 * y, 2 * x + 1, 2 * y + 2));
                        break;

                    case 0xA:
                    case 0x5:
                        // 0101=1010 Saddle; must test center point to decide
                        {
                            QPointF pt =
                                0.5 * (grid_points[x * H + y].coords +
                                       grid_points[(x + 1) * H + y + 1].coords);
                            RenderPoint rp = queryPoint(pt, nav);
                            const RenderPoint &class_rep =
                                (code == 0x5) ? grid_points[x * H + y]
                                              : grid_points[x * H + y + 1];

                            bool center_inclass =
                                crease_depth(class_rep, rp, cos_alpha,
                                             min_jump) < 0;
                            if ((code == 0x5) ^ center_inclass) {
                                // Diagonals along (x,y+1)<->(x+1,y)
                                segs.push_back(
                                    QLine(2 * x + 1, 2 * y, 2 * x, 2 * y + 1));
                                segs.push_back(QLine(2 * x + 1, 2 * y + 2,
                                                     2 * x + 2, 2 * y + 1));
                            } else {
                                // Diagonals along (x,y)<->(x+1,y+1)
                                segs.push_back(QLine(2 * x + 1, 2 * y,
                                                     2 * x + 2, 2 * y + 1));
                                segs.push_back(QLine(2 * x + 1, 2 * y + 2,
                                                     2 * x, 2 * y + 1));
                            }
                        }

                        break;

                    default:
                        qFatal("Shouldn't happen");
                    }
                }
            }

            const QVector<QPolygon> &loops = size_sorted_loops_from_seg(segs);

            qDebug("class %d subclass %d nsegs %d nloops %d", region.class_no,
                   subcls, segs.size(), loops.size());

            Subregion subreg;
            subreg.subclass_no = subcls;
            // Default set up color info
            subreg.gradient_type = GradientType::gSolid;
            subreg.solid_color = qRgb(255, 255, 255);
            subreg.linear_angle = 0.;
            subreg.linear_start = 0.;
            subreg.linear_stop = 0.;
            subreg.linear_colors.clear();

            // Determine if subregion is clipped region
            subreg.is_clipped_patch = false;
            subreg.representative_coord = subreg_representative_coords[subcls];
            QPoint src = subreg_representative_coords[subcls];
            const RenderPoint &subreg_class_rep =
                grid_points[src.x() * H + src.y()];
            if (subreg_class_rep.ray.intersections) {
                subreg.is_clipped_patch = subreg_class_rep.ray.front_clipped;
            }
            if (subreg_class_rep.ray.N > 0) {
                int ec = subreg_class_rep.ray.intersections[0].ecode;
                if (ec < 0 || !view_data.elements[ec].visible)
                    subreg.is_clipped_patch = false;
            }

            subreg.boundaries.clear();
            for (int zzz = 0; zzz < loops.size(); zzz++) {
                const QPolygon &poly = loops[zzz];
                bool is_exterior_loop = (zzz == 0);
                QVector<RenderPoint> bound;

                // Only crease edges visible...
                for (QPoint line_marker : poly) {
                    QPoint p0 = line_marker / 2;
                    QPoint p1 = line_marker - line_marker / 2;

                    RenderPoint q0 = getPoint(p0);
                    RenderPoint q1 = getPoint(p1);

                    RenderPoint adj_point;
                    if (!crease_edge_map.contains(line_marker)) {
                        // Edge to out of region
                        // Pull typematch subdiv
                        RenderPoint lim_out;
                        bracketEdge(
                            q0.region_class == region.class_no ? q0 : q1,
                            q0.region_class == region.class_no ? q1 : q0,
                            &adj_point, &lim_out, nav);
                        adj_point.ideal_color =
                            calculateBoundaryColor(adj_point, lim_out);
                        adj_point.show_point = false;
                    } else if (!crease_edge_map[line_marker]) {
                        // Interior edge, not conflict
                        // Select midpoint of q0 and q1
                        adj_point =
                            queryPoint(0.5 * (q0.coords + q1.coords), nav);
                        adj_point.ideal_color =
                            calculateInteriorColor(adj_point);
                        adj_point.show_point = false;
                    } else {
                        // Interior edge, crease.
                        // Pick error type (normal or displacement)
                        // and bisect on that parameter
                        // Color is inside point color
                        RenderPoint lim_out;
                        bracketCrease(
                            q0.subregion_class == subreg.subclass_no ? q0 : q1,
                            q0.subregion_class == subreg.subclass_no ? q1 : q0,
                            &adj_point, &lim_out, nav);
                        adj_point.ideal_color =
                            calculateInteriorColor(adj_point);
                        adj_point.show_point = true;
                    }

                    FColor aic = adj_point.ideal_color;
                    adj_point.ideal_color =
                        FColor(aic.redF(), aic.greenF(), aic.blueF());
                    // Oh: rpoint, alpha=0 signifies invisible.
                    bound.push_back(adj_point);
                }

                // Filter out very close points
                double epsilon = (1. / std::max(W, H)) * 1e-5;
                bound = remove_rloop_duplicates(bound, epsilon);

                // Corner handling. First, we locate places for corners
                QVector<CornerCandidate> corner_cands =
                    find_rloop_corners(bound, is_exterior_loop, epsilon);

                QVector<RenderPoint> corner_points;
                QVector<int> corner_sidxs;
                for (int i = 0; i < corner_cands.size(); i++) {
                    int k = corner_cands[i].index;
                    const RenderPoint &p0 = bound[k];
                    const RenderPoint &p1 = bound[(k + 1) % bound.size()];
                    RenderPoint qA = queryPoint(corner_cands[i].inside, nav);

                    // Both types and subregion must match
                    if (!typematch(p0, qA) || !typematch(p1, qA)) {
                        continue;
                    }
                    bool with0 = crease_depth(qA, p0, cos_alpha, min_jump) < 0;
                    bool with1 = crease_depth(qA, p1, cos_alpha, min_jump) < 0;
                    if (!with0 || !with1) {
                        continue;
                    }

                    // Outside state checked only in relation to inside
                    RenderPoint qB = queryPoint(corner_cands[i].outside, nav);
                    RenderPoint adj_point;
                    if (!typematch(qB, qA)) {
                        // Region boundary
                        RenderPoint lim_out;
                        bracketEdge(qA, qB, &adj_point, &lim_out, nav);
                        adj_point.ideal_color =
                            calculateBoundaryColor(adj_point, lim_out);
                        adj_point.show_point = false;
                    } else if (crease_depth(qA, qB, cos_alpha, min_jump) < 0) {
                        // No crease found, average it
                        adj_point =
                            queryPoint(0.5 * (qA.coords + qB.coords), nav);
                        adj_point.ideal_color =
                            calculateInteriorColor(adj_point);
                        adj_point.show_point = false;
                    } else {
                        // Crease boundary
                        RenderPoint lim_out;
                        bracketCrease(qA, qB, &adj_point, &lim_out, nav);
                        adj_point.ideal_color =
                            calculateInteriorColor(adj_point);
                        adj_point.show_point = true;
                    }

                    corner_points.push_back(adj_point);
                    corner_sidxs.push_back(k);
                }

                // Merge the two

                bound = rloop_merge(bound, corner_points, corner_sidxs);

                if (!subreg.boundaries.size()) {
                    // Outside (largest) loop must be reversed relative to all
                    // else
                    bound = reverse_rloop(bound);
                }
                subreg.boundaries.push_back(bound);
            }
            region.subregions.push_back(subreg);
        }
    }
    delete nav;
    delete[] flood_fill_data;
    delete[] crease_edge_count;

    int S = std::max(10, 1024 / std::max(W - 1, H - 1));
    QImage image((W - 1) * S, (H - 1) * S, QImage::Format_ARGB32);
    {
        QRgb *colors = new QRgb[grid_nclasses];
        for (int i = 0; i < grid_nclasses; i++) {
            colors[i] = randColor().rgba();
        }
        QRgb *scolors = new QRgb[nsubclasses];
        for (int i = 0; i < nsubclasses; i++) {
            scolors[i] = randColor().rgba();
        }
        // fill colors?

        QPainter p(&image);
        p.fillRect(image.rect(), Qt::gray);
        for (int x = 0; x < W; x++) {
            for (int y = 0; y < H; y++) {
                QPointF center(S * x, S * y);
                double radius = S * 0.3;
                QRgb col = colors[grid_points[x * H + y].region_class];
                int sclass = grid_points[x * H + y].subregion_class;
                if (sclass < 0)
                    continue;
                QRgb scol = scolors[sclass];
                p.setPen(col);
                p.setBrush(QBrush(scol));
                p.drawEllipse(center, radius, radius);
            }
        }
        for (int x = 0; x < W; x++) {
            for (int y = 0; y < H; y++) {
                for (int k = 0; k < 2; k++) {
                    QPoint line_marker(2 * x + (1 - k), 2 * y + k);
                    if (crease_edge_map.count(line_marker)) {
                        QPoint f0 = line_marker / 2;
                        QPoint f1 = line_marker - f0;
                        bool crease = crease_edge_map[line_marker];
                        QPen pen;
                        if (crease) {
                            pen = QPen(Qt::white);
                            pen.setWidth(2.0);
                        } else {
                            pen = QPen(Qt::black);
                            pen.setWidth(1.0);
                        }

                        p.setPen(pen);
                        p.drawLine(f0 * S, f1 * S);
                    }
                }
            }
        }

        delete[] colors;
        delete[] scolors;
    }

    emit produceImagePhase(image, QString("Creases completed"), nqueries,
                           false);

    // Each line is identified via (2 * x + 1, 2 * y) or (2 * x, 2 * y + 1)
}

static QPolygonF shift_and_scale_loop(const QVector<RenderPoint> &loop,
                                      const QTransform &transform) {
    QPolygonF a;
    for (const RenderPoint &p : loop) {
        a.push_back(transform.map(p.coords));
    }
    return a;
}
static QVector<QPolygonF>
boundary_loops_for_region(const Region &region, const QTransform &transform) {
    QVector<QPolygonF> n;
    n.push_back(shift_and_scale_loop(region.exterior, transform));
    for (const QVector<RenderPoint> &loop : region.interior) {
        n.push_back(shift_and_scale_loop(loop, transform));
    }
    return n;
}
static QVector<QPolygonF>
boundary_loops_for_subregion(const Subregion &subreg,
                             const QTransform &transform) {
    QVector<QPolygonF> n;
    for (const QVector<RenderPoint> &loop : subreg.boundaries) {
        n.push_back(shift_and_scale_loop(loop, transform));
    }
    return n;
}

static QString color_hex_name_rgb(QRgb color) {
    int r = qRed(color);
    int g = qGreen(color);
    int b = qBlue(color);
    const char *numbers = "0123456789abcdef";
    return QString("#%1%2%3%4%5%6")
        .arg(numbers[r / 16])
        .arg(numbers[r % 16])
        .arg(numbers[g / 16])
        .arg(numbers[g % 16])
        .arg(numbers[b / 16])
        .arg(numbers[b % 16]);
}

static QPolygonF simplify_poly(const QPolygonF &orig, double max_error,
                               bool loop = true) {
    // TODO: bezier support; use rolling approach, error relative to
    // 50% subsample with 50% subsample using curvature tricks as well?
    if (0)
        return orig;

    // First, locate sharpest angle; will be our starting point
    const int n = orig.size();
    int sharpest = 0;
    if (loop) {
        qreal recsh = 1.;
        for (int i1 = 0; i1 < n; i1++) {
            int i0 = (i1 + n - 1) % n, i2 = (i1 + 1) % n;
            QPointF d01 = orig[i1] - orig[i0];
            QPointF d12 = orig[i2] - orig[i1];
            qreal len01 = d01.x() * d01.x() + d01.y() * d01.y();
            qreal len12 = d12.x() * d12.x() + d12.y() * d12.y();
            qreal lenpro = std::sqrt(len01 * len12);
            if (lenpro <= 0.)
                continue;

            double qr = QPointF::dotProduct(d01, d12) / lenpro;
            if (qr < recsh) {
                sharpest = i1;
                recsh = qr;
            }
        }
    } else {
        sharpest = 0;
    }

    // Then iterative lookahead, while calculating maximum offset
    QPolygonF poly;
    poly.push_back(orig[sharpest]);
    for (int is = 0; is < n; is++) {
        QPointF ps = orig[(is + sharpest) % n];
        int js = is + 1;
        for (; js <= n; js++) {
            // If line fails to match, js-=1, break
            QPointF pe = orig[(js + sharpest) % n];
            QPointF delta = pe - ps;
            qreal dlen =
                std::sqrt(delta.x() * delta.x() + delta.y() * delta.y());
            if (dlen <= 0.) {
                // Bad case; but we progress anyway
                break;
            }
            bool broken = false;
            // Check points in [is+1,js-1]
            for (int zs = is + 1; zs < js; zs++) {
                QPointF pc = orig[(zs + sharpest) % n];
                // Compute distance between [ps,pe] and pc
                qreal soff = (delta.y() * pc.x() - delta.x() * pc.y() +
                              pe.x() * ps.y() - pe.y() * ps.x()) /
                             dlen;
                if (std::abs(soff) > max_error) {
                    broken = true;
                    continue;
                }
            }
            if (broken) {
                js--;
                break;
            }
        }
        // Jump to ext pt
        if (js > n)
            js = n;
        poly.push_back(orig[(js + sharpest) % n]);
        is = js - 1;
    }

    if (!loop) {
        poly.pop_back();
    }

    return poly;
}

static QString svg_path_text_from_polygons(const QVector<QPolygonF> &loops,
                                           bool close_loops = true) {
    const int fprec = 7;
    QStringList path_string;
    for (const QPolygonF &poly : loops) {
        QPolygonF simpath = simplify_poly(poly, 1e-5, close_loops);

        QPointF s = simpath[0];
        path_string.append(QString("M%1,%2")
                               .arg(s.x(), 0, 'f', fprec)
                               .arg(s.y(), 0, 'f', fprec));
        for (int i = 1; i < simpath.size(); i++) {
            QPointF q = simpath[i];
            path_string.append(QString("L%1,%2")
                                   .arg(q.x(), 0, 'f', fprec)
                                   .arg(q.y(), 0, 'f', fprec));
        }
        // Note: A rx ry x-axis-rotation large-arc-flag sweep-flag x y
        // gives elliptical arc, very suitable for path compression
        if (close_loops) {
            path_string.append("Z");
        }
    }
    return path_string.join(" ");
}

static float compute_histogram_angle(const Region &region,
                                     const Subregion &subreg,
                                     RenderPoint *grid_points,
                                     const QSize &grid_size) {
    // TODO: ensure subregion bounding & point lookup is fast & local

    const int W = grid_size.width(), H = grid_size.height();

    // Here we compute the circular mean of the gradient line direction.
    // Thankfully RP1=S1.
    double sum_sin = 0., sum_cos = 0.;
    for (int x = region.xmin; x <= region.xmax; x++) {
        for (int y = region.ymin; y <= region.ymax; y++) {
            // There are four dx/dy kernels feasible;
            // we add them all together
            RenderPoint &base = grid_points[x * H + y];
            if (base.region_class != region.class_no)
                continue;
            if (base.subregion_class != subreg.subclass_no)
                continue;

            int dxs[4] = {1, 1, -1, -1};
            int dys[4] = {1, -1, -1, 1};
            for (int k = 0; k < 4; k++) {
                int xnx = x + dxs[k], yny = y + dys[k];
                if (xnx >= W || yny >= H || xnx < 0 || yny < 0)
                    continue;

                RenderPoint &pdx = grid_points[xnx * H + y];
                if (pdx.region_class != region.class_no)
                    continue;
                if (pdx.subregion_class != subreg.subclass_no)
                    continue;

                RenderPoint &pdy = grid_points[x * H + yny];
                if (pdy.region_class != region.class_no)
                    continue;
                if (pdy.subregion_class != subreg.subclass_no)
                    continue;

                // We consider derivatives in magnitude, as intra-subregion
                // color shifts tend to be scaling-only
                float mag_base = base.ideal_color.magnitude();
                float mag_dx = pdx.ideal_color.magnitude();
                float mag_dy = pdy.ideal_color.magnitude();

                double dx = (mag_dx - mag_base) / dxs[k];
                double dy = (mag_dy - mag_base) / dys[k];
                double mag2 = (dx * dx + dy * dy);
                if (mag2 <= 0.)
                    continue;
                // S1 -> RP1=S1 doubles the angle
                // Use identity (sin(2 * atan2(dy, dx))) = 2*dx*dy/mag2
                sum_sin += 2 * dx * dy / mag2;
                sum_cos += (dx * dx - dy * dy) / mag2;
            }
        }
    }
    if (sum_sin == 0. && sum_cos == 0.) {
        // Uniformly distributed angles, or just none at all
        return 0.;
    }
    double angle = std::atan2(sum_sin, sum_cos);
    // The half-pi shift gives an angle perpendicular to dominant
    // gradient direction
    return std::fmod(-angle / 2 + M_PI / 2, M_PI);
}

static float proj_angle(const QPointF &p, float angle) {
    float proj_x = -std::sin(angle);
    float proj_y = std::cos(angle);
    return p.x() * proj_x + p.y() * proj_y;
}

static void fill_linear_histogram_for_angle(
    Subregion &region, const QVector<QPair<FColor, QPointF>> &interior_colors,
    const QSize &gridsize, double angle) {
    int S = std::max(gridsize.width(), gridsize.height()) - 1;
    double grid_spacing = 1.0 / S;

    // Determine histogram parameters
    double tmin = std::numeric_limits<float>::infinity();
    double tmax = -std::numeric_limits<float>::infinity();
    for (const QPair<FColor, QPointF> &p : interior_colors) {
        double t = proj_angle(p.second, angle);
        tmin = std::min(tmin, t);
        tmax = std::max(tmax, t);
    }
    double min_spacing = 1.0 * grid_spacing;
    int nsteps = std::max(2, 1 + (int)std::floor((tmax - tmin) / min_spacing));

    region.linear_angle = angle;
    region.linear_start = tmin;
    region.linear_stop = tmax;
    region.linear_colors.clear();
    // Compute average color
    double *histogram_r = new double[nsteps];
    double *histogram_b = new double[nsteps];
    double *histogram_g = new double[nsteps];
    double *histogram_a = new double[nsteps];
    double *histogram_w = new double[nsteps];
    for (int i = 0; i < nsteps; i++) {
        histogram_r[i] = 0.;
        histogram_g[i] = 0.;
        histogram_b[i] = 0.;
        histogram_a[i] = 0.;
        histogram_w[i] = 0.;
    }
    for (const QPair<FColor, QPointF> &p : interior_colors) {
        QColor color(p.first.rgba());
        double t = proj_angle(p.second, angle);
        double cell = (nsteps - 1.) * (t - region.linear_start) /
                      (region.linear_stop - region.linear_start);
        int lidx = (int)std::floor(cell);
        lidx = std::max(0, std::min(nsteps - 2, lidx));
        double frac = std::max(0.0, std::min(1.0, cell - lidx));
        histogram_r[lidx] += color.redF() * (1. - frac);
        histogram_g[lidx] += color.greenF() * (1. - frac);
        histogram_b[lidx] += color.blueF() * (1. - frac);
        histogram_a[lidx] += color.alphaF() * (1. - frac);
        histogram_w[lidx] += 1. - frac;
        histogram_r[lidx + 1] += color.redF() * frac;
        histogram_g[lidx + 1] += color.greenF() * frac;
        histogram_b[lidx + 1] += color.blueF() * frac;
        histogram_a[lidx + 1] += color.alphaF() * frac;
        histogram_w[lidx + 1] += frac;
    }
    for (int i = 0; i < nsteps; i++) {
        if (histogram_w[i] > 0.) {
            double iw = 1. / histogram_w[i];
            double r = std::max(0.0, std::min(1.0, iw * histogram_r[i]));
            double g = std::max(0.0, std::min(1.0, iw * histogram_g[i]));
            double b = std::max(0.0, std::min(1.0, iw * histogram_b[i]));
            double a = std::max(0.0, std::min(1.0, iw * histogram_a[i]));
            QColor color = QColor::fromRgbF(r, g, b, a);
            region.linear_colors.push_back(color.rgba());
        } else {
            // TODO: average from neighbors
            region.linear_colors.push_back(qRgba(0, 0, 0, 0));
        }
    }
    delete[] histogram_r;
    delete[] histogram_b;
    delete[] histogram_g;
    delete[] histogram_a;
    delete[] histogram_w;
}

static int square(int s) { return s * s; }
static int rgba_distance2(QRgb a, QRgb b) {
    return square(qRed(a) - qRed(b)) + square(qGreen(a) - qGreen(b)) +
           square(qBlue(a) - qBlue(b)) + square(qAlpha(a) - qAlpha(b));
}

static long
compute_gradient_error(const Subregion &region,
                       const QVector<QPair<FColor, QPointF>> &interior_colors) {
    int nsteps = region.linear_colors.size();
    long sqe = 0.;
    for (const QPair<FColor, QPointF> &p : interior_colors) {
        QRgb color = p.first.rgba();

        QRgb pred;
        if (region.gradient_type == GradientType::gSolid) {
            pred = region.solid_color;
        } else if (region.gradient_type == GradientType::gLinear) {
            double t = proj_angle(p.second, region.linear_angle);

            double cell = (nsteps - 1.) * (t - region.linear_start) /
                          (region.linear_stop - region.linear_start);
            int lidx = (int)std::floor(cell);
            lidx = std::max(0, std::min(nsteps - 2, lidx));
            double frac = std::max(0.0, std::min(1.0, cell - lidx));

            QRgb a = region.linear_colors[lidx];
            QRgb b = region.linear_colors[lidx + 1];

            pred = qRgba(qRed(a) * (1. - frac) + qRed(b) * frac,
                         qGreen(a) * (1. - frac) + qGreen(b) * frac,
                         qBlue(a) * (1. - frac) + qBlue(b) * frac,
                         qAlpha(a) * (1. - frac) + qAlpha(b) * frac);
        } else {
            qFatal("Bad case");
        }

        sqe += rgba_distance2(color, pred);
    }
    return sqe;
}

static FColor compute_mean_color(const QVector<RenderPoint> &pts) {
    double bnet_r = 0, bnet_g = 0, bnet_b = 0, bnet_a = 0;
    for (const RenderPoint &rp : pts) {
        bnet_r += rp.ideal_color.redF();
        bnet_g += rp.ideal_color.greenF();
        bnet_b += rp.ideal_color.blueF();
        bnet_a += rp.ideal_color.alphaF();
    }
    int n = pts.size();
    return FColor(bnet_r / n, bnet_g / n, bnet_b / n, bnet_a / n);
}

static QVector<QVector<RenderPoint>>
extract_visible_boundary(const QVector<RenderPoint> &pts) {
    // We already assume there is at least one invisible point

    const int n = pts.size();
    int start = -1;
    for (int i = 0; i < n; i++) {
        if (!pts[i].show_point) {
            start = i;
            break;
        }
    }
    if (start < 0) {
        qWarning("Extract visible boundary, invariant fail");
        start = 0;
    }

    QVector<QVector<RenderPoint>> res;
    int i = 0;
    while (i < n) {
        for (; i < n; i++) {
            if (pts[(i + start) % n].show_point) {
                break;
            }
        }
        if (i >= n) {
            break;
        }

        QVector<RenderPoint> subpath;
        for (; i < n; i++) {
            if (pts[(i + start) % n].show_point) {
                subpath.push_back(pts[(i + start) % n]);
            } else {
                break;
            }
        }
        res.push_back(subpath);
    }
    return res;
}

static void run_command(const QStringList &command) {
    QProcess proc;
    proc.start(command.first(), command.mid(1));
    proc.waitForFinished();
}

void VectorTracer::computeGradients() {
    qDebug("Computing gradients");

    int W = grid_size.width(), H = grid_size.height();
    for (Region &region : region_list) {
        // First, compute best estimate region interior average color
        for (Subregion &subreg : region.subregions) {

            QVector<QPair<FColor, QPointF>> interior_colors;
            float net_r = 0, net_g = 0, net_b = 0, net_a = 0;
            int net_count = 0;
            for (int x = region.xmin; x <= region.xmax; x++) {
                for (int y = region.ymin; y <= region.ymax; y++) {
                    if (grid_points[x * H + y].region_class != region.class_no)
                        continue;
                    if (grid_points[x * H + y].subregion_class !=
                        subreg.subclass_no)
                        continue;
                    RenderPoint &pt = grid_points[x * H + y];

                    interior_colors.push_back(
                        QPair<FColor, QPointF>(pt.ideal_color, pt.coords));

                    net_r += pt.ideal_color.redF();
                    net_g += pt.ideal_color.greenF();
                    net_b += pt.ideal_color.blueF();
                    net_a += pt.ideal_color.alphaF();
                    net_count += 1;
                }
            }
            subreg.solid_color = FColor(net_r / net_count, net_g / net_count,
                                        net_b / net_count, net_a / net_count)
                                     .rgba();
            subreg.gradient_type = GradientType::gSolid;

            long cost = compute_gradient_error(subreg, interior_colors);
            if (cost > 0) {
                subreg.gradient_type = GradientType::gLinear;

                double angle = compute_histogram_angle(region, subreg,
                                                       grid_points, grid_size);
                fill_linear_histogram_for_angle(subreg, interior_colors,
                                                grid_size, angle);
                qDebug("class %d subclass %d angle %f nsteps %d",
                       region.class_no, subreg.subclass_no, angle,
                       subreg.linear_colors.size());
            }

            // Fall back to solid if uniform
            if (subreg.gradient_type == GradientType::gLinear) {
                QRgb icolor = subreg.linear_colors[0];
                bool different = false;
                for (int i = 1; i < subreg.linear_colors.size(); i++) {
                    if (icolor != subreg.linear_colors[i]) {
                        different = true;
                        break;
                    }
                }
                if (!different) {
                    subreg.gradient_type = GradientType::gSolid;
                }
            }
        }

        // Finally, compute boundary average colors
        region.meanExteriorColor = compute_mean_color(region.exterior).rgba();
        for (const QVector<RenderPoint> &loop : region.interior) {
            region.meanInteriorColors.push_back(
                compute_mean_color(loop).rgba());
        }
    }

    //
    // Since QtSvg is incapable of handling clip paths/pattern fills like every
    // other SVG handling code, we roll our own SVG generator.
    //
    // Note: As QPainter supports this, fixing it in QtSvg on
    // render/generator sides isn't very hard.
    //

    const int S = std::max(W, H);
    const double px_per_mm = 3.78;
    double T = 45 * px_per_mm;
    const QPointF offset =
        2 * QPointF(0.5 * W / S + 1. / S, (0.5 * H / S + 1. / S));
    QTransform transf(T / 2, 0, 0, -T / 2, offset.x() * T / 2,
                      offset.y() * T / 2);
    QTransform dirtransf(transf.m11(), transf.m12(), transf.m21(), transf.m22(),
                         0., 0.);

    // TODO: use `point_to_grid_coord` and then rescale the grid
    QRectF viewbox(0., 0., T * (W + 2) / (double)S, T * (H + 2) / (double)S);

    // Gradient -- a piecewise composite of linear and radial gradients.
    // It's easy to fit a single gradient on a local patch;
    // then add a dividing crease along the mismatch between regions,
    // where opposing candidates are closest. (Again, march-cube style?,
    // then refined ? or A/B decision style.

    QTemporaryFile tfile("/tmp/render.svg.XXXXXX");

    tfile.open();
    qDebug("Rendering final image, %d regions, to %s", region_list.size(),
           tfile.fileName().toUtf8().constData());
    {
        QTextStream s(&tfile);

        s << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n "
             "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";

        s << QString("<svg version=\"1.1\" width=\"%1\" height=\"%2\" >\n")
                 .arg(viewbox.width())
                 .arg(viewbox.height());

        s << "  <defs>\n";

        // Hatching fill overlays (by type)
        QMap<CompactNormal, int> normal_directions;
        for (const Region &region : region_list) {
            for (const Subregion &subreg : region.subregions) {
                // consider average normal over subregion
                if (subreg.is_clipped_patch) {
                    QPoint src = subreg.representative_coord;
                    const CompactNormal &normal =
                        grid_points[src.x() * H + src.y()]
                            .ray.intersections[0]
                            .normal;
                    if (!normal_directions.contains(normal)) {
                        normal_directions[normal] = normal_directions.size();
                    }
                }
            }
        }
        const G4ThreeVector &updown = view_data.orientation.rowY();
        const G4ThreeVector &lright = view_data.orientation.rowZ();
        for (CompactNormal nd : normal_directions.keys()) {
            int idx = normal_directions[nd];
            const G4ThreeVector &normal = nd;

            const G4ThreeVector &orthA = normal.orthogonal().unit();
            const G4ThreeVector &orthB = normal.cross(orthA).unit();
            const G4ThreeVector &pattern_dir = (orthA + orthB);
            // Convert to unit viewport space
            double udf = pattern_dir.dot(updown) / view_data.scale,
                   lrf = pattern_dir.dot(lright) / view_data.scale;
            double nudf, nlrf;
            // Convert to SVG space
            dirtransf.map(udf, lrf, &nudf, &nlrf);
            // Compute spacing/angle
            double spacing = 10 * std::sqrt(nudf * nudf + nlrf * nlrf);
            double angle = std::atan2(nudf, nlrf);

            //            double hatch_width = spacing / 2;

            double hatch_width = T / 40.;
            qDebug("spacing %f hatch_width %f", spacing, hatch_width);

            s << QString("    <pattern id=\"hatching%1\" width=\"%2\" "
                         "height=\"10\" patternTransform=\"rotate(%3 0 "
                         "0)\" patternUnits=\"userSpaceOnUse\">\n")
                     .arg(idx)
                     .arg(hatch_width)
                     .arg(angle * 180 / CLHEP::pi);
            s << QString("       <line x1=\"0\" y1=\"0\" x2=\"0\" "
                         "y2=\"10\" stroke=\"%1\" stroke-width=\"%2\"/>\n")
                     .arg(color_hex_name_rgb(qRgb(0, 0, 0)))
                     .arg(hatch_width);
            s << QString("    </pattern>\n");
        }

        // Gradients
        for (const Region &region : region_list) {
            for (const Subregion &subreg : region.subregions) {
                if (subreg.gradient_type == GradientType::gLinear) {
                    // Locate center point
                    QPointF sscenter;
                    for (const RenderPoint &p : subreg.boundaries[0]) {
                        sscenter += p.coords;
                    }
                    sscenter /= subreg.boundaries[0].size();

                    float ssf = proj_angle(sscenter, subreg.linear_angle);

                    float lsin = std::sin(subreg.linear_angle),
                          lcos = std::cos(subreg.linear_angle);
                    QPointF ssfstart = sscenter + (subreg.linear_start - ssf) *
                                                      QPointF(-lsin, lcos);
                    QPointF ssfstop = sscenter + (subreg.linear_stop - ssf) *
                                                     QPointF(-lsin, lcos);

                    ssfstart = transf.map(ssfstart);
                    ssfstop = transf.map(ssfstop);

                    s << QString("   <linearGradient id=\"gradient%1_%2\" "
                                 "x1=\"%3\" "
                                 "y1=\"%4\" x2=\"%5\" y2=\"%6\" "
                                 "spreadMethod=\"%7\" "
                                 "gradientUnits=\"userSpaceOnUse\">\n")
                             .arg(region.class_no)
                             .arg(subreg.subclass_no)
                             .arg(ssfstart.x())
                             .arg(ssfstart.y())
                             .arg(ssfstop.x())
                             .arg(ssfstop.y())
                             .arg(0 ? "repeat" : "pad");
                    int nsteps = subreg.linear_colors.size();
                    for (int i = 0; i < nsteps; i++) {
                        double stop_pos = i / (nsteps - 1.);

                        if (i > 0 && i < nsteps - 1 &&
                            subreg.linear_colors[i] ==
                                subreg.linear_colors[i - 1]) {
                            // Skip redundant points
                            continue;
                        }

                        double alpha = QColor(subreg.linear_colors[i]).alphaF();
                        s << QString(
                                 "    <stop offset=\"%1%\" stop-color=\"%2\" "
                                 "stop-opacity=\"%3\"/>\n")
                                 .arg(100 * stop_pos)
                                 .arg(color_hex_name_rgb(
                                     subreg.linear_colors[i]))
                                 .arg(alpha);
                    }
                    s << QString("  </linearGradient>\n");
                }
            }
        }
        s << "  </defs>\n";

        // Interior regions, clipped by boundaries
        s << QString("<g id=\"interiors\" fill-opacity=\"1\" stroke=\"none\"  "
                     "fill-rule=\"evenodd\">\n");
        for (const Region &region : region_list) {
            // Compute region limit
            for (const Subregion &subreg : region.subregions) {
                const QVector<QPolygonF> &subloops =
                    boundary_loops_for_subregion(subreg, transf);

                QString fill_desc;
                if (subreg.gradient_type == GradientType::gSolid) {
                    fill_desc = color_hex_name_rgb(subreg.solid_color);
                } else if (subreg.gradient_type == GradientType::gLinear) {
                    fill_desc = QString("url(#gradient%1_%2)")
                                    .arg(region.class_no)
                                    .arg(subreg.subclass_no);
                } else {
                    qFatal("Unsupported gradient type");
                }

                s << QString(
                         "  <path id=\"region%1_%2\" fill=\"%3\" d=\"%4\"/>\n")
                         .arg(region.class_no)
                         .arg(subreg.subclass_no)
                         .arg(fill_desc)
                         .arg(svg_path_text_from_polygons(subloops));

                // Overlay mix with hatching
                if (subreg.is_clipped_patch) {
                    QPoint src = subreg.representative_coord;
                    const CompactNormal &normal =
                        grid_points[src.x() * H + src.y()]
                            .ray.intersections[0]
                            .normal;
                    int idx = normal_directions[normal];
                    s << QString("  <path id=\"region_hatch%1_%2\" "
                                 "fill=\"url(#hatching%3)\" opacity=\"0.2\" "
                                 "d=\"%4\"/>\n")
                             .arg(region.class_no)
                             .arg(subreg.subclass_no)
                             .arg(idx)
                             .arg(svg_path_text_from_polygons(subloops));
                }
                if (subreg.is_clipped_patch) {
                }
            }
        }
        s << QString("</g>\n");

        // Boundaries
        s << QString("<g id=\"boundaries\" fill=\"none\" stroke-width=\"%1\" "
                     "fill-rule=\"evenodd\" stroke-linecap=\"square\" "
                     "stroke-linejoin=\"miter\">\n")
                 .arg(T * 0.003);
        for (const Region &region : region_list) {
            const QVector<QPolygonF> &loops =
                boundary_loops_for_region(region, transf);
            for (int i = 0; i < loops.size(); i++) {
                QVector<QPolygonF> solo;
                solo.push_back(loops[i]);
                QRgb c = i ? region.meanInteriorColors[i - 1]
                           : region.meanExteriorColor;
                s << QString("    <path id=\"edge%1_%2\" stroke=\"%3\" "
                             "d=\"%4\"/>\n")
                         .arg(region.class_no)
                         .arg(i)
                         .arg(color_hex_name_rgb(c))
                         .arg(svg_path_text_from_polygons(solo));
            }
            // Creases
            for (const Subregion &subreg : region.subregions) {
                for (int l = 0; l < subreg.boundaries.size(); l++) {
                    const QVector<RenderPoint> &bound = subreg.boundaries[l];
                    bool all_visible = true;
                    for (const RenderPoint &pt : bound) {
                        if (!pt.show_point) {
                            all_visible = false;
                            break;
                        }
                    }
                    if (all_visible) {
                        QPolygonF loop = shift_and_scale_loop(bound, transf);
                        QVector<QPolygonF> solo;
                        solo.push_back(loop);

                        FColor crease_col = compute_mean_color(bound);
                        crease_col = FColor(
                            crease_col.redF() / 2, crease_col.greenF() / 2,
                            crease_col.blueF() / 2, crease_col.alphaF());

                        s << QString(
                                 "    <path id=\"edge%1_%2_%3\" stroke=\"%4\" "
                                 "d=\"%5\"/>\n")
                                 .arg(region.class_no)
                                 .arg(subreg.subclass_no)
                                 .arg(l)
                                 .arg(color_hex_name_rgb(crease_col.rgba()))
                                 .arg(svg_path_text_from_polygons(solo));
                    } else {
                        const QVector<QVector<RenderPoint>> &visible_segments =
                            extract_visible_boundary(bound);
                        for (int m = 0; m < visible_segments.size(); m++) {
                            FColor crease_col =
                                compute_mean_color(visible_segments[m]);
                            crease_col = FColor(
                                crease_col.redF() / 2, crease_col.greenF() / 2,
                                crease_col.blueF() / 2, crease_col.alphaF());
                            QVector<QPolygonF> solo;
                            solo.push_back(shift_and_scale_loop(
                                visible_segments[m], transf));
                            s << QString("    <path id=\"edge%1_%2_%3_%4\" "
                                         "stroke=\"%5\" stroke-width=\"%6\" "
                                         "d=\"%7\"/>\n")
                                     .arg(region.class_no)
                                     .arg(subreg.subclass_no)
                                     .arg(l)
                                     .arg(m)
                                     .arg(color_hex_name_rgb(crease_col.rgba()))
                                     .arg(T * 0.002)
                                     .arg(svg_path_text_from_polygons(solo,
                                                                      false));
                        }
                    }
                }
            }
        }
        s << QString("</g>\n");

        s << "</svg>\n";
    }

    QTime time_start = QTime::currentTime();
    if (file_name.size()) {
        if (file_name.endsWith(".pdf")) {
            qDebug("Converting %s to %s", tfile.fileName().toUtf8().constData(),
                   file_name.toUtf8().constData());
            run_command(QStringList()
                        << "inkscape"
                        << "--export-pdf=" + file_name << tfile.fileName());
        } else {
            qDebug("Copying %s to %s", tfile.fileName().toUtf8().constData(),
                   file_name.toUtf8().constData());
            tfile.copy(file_name);
        }
    } else {
        // Automatically pick a name to copy the SVG to
        long i = 0;
        QString fmt("vector%1.svg");
        QString dst;
        do {
            dst = fmt.arg(i, 4, 10, QChar('0'));
            i++;
        } while (QFileInfo(dst).exists());
        qDebug("Copying %s to %s", tfile.fileName().toUtf8().constData(),
               dst.toUtf8().constData());
        tfile.copy(dst);
    }
    {
        // Use Inkscape to generate an image preview, because QtSVG is generally
        // incapable; this should take half a second or so
        QTemporaryFile png_file("/tmp/preview.png.XXXXXX");
        png_file.open();
        run_command(QStringList() << "inkscape"
                                  << "--export-png=" + png_file.fileName()
                                  << "--export-dpi=800" << tfile.fileName());

        QImage grad_image = QImage(png_file.fileName(), "png");

        emit produceImagePhase(grad_image, QString("Render completed"),
                               nqueries, true);
    }

    QTime time_end = QTime::currentTime();
    qDebug("Conversion took %f seconds", time_start.msecsTo(time_end) * 1e-3);
}
