/* SPDX-License-Identifier: GPL-3.0-only */
#include "VectorTrace.hh"

#include "Bezier.hh"
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
        emit produceImagePhase(img.copy(), QStringLiteral("Grid completed"),
                               nqueries, false);

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
        area += (loop[k1].x() - loop[k0].x()) * (loop[k1].y() + loop[k0].y());
    }
    return area;
}
static double signed_area(const QVector<RenderPoint> &loop) {
    double area = 0;
    for (int k0 = 0; k0 < loop.size(); k0++) {
        int k1 = (k0 + 1) % loop.size();
        area += (loop[k1].coords.x() - loop[k0].coords.x()) *
                (loop[k1].coords.y() + loop[k0].coords.y());
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

static constexpr QPoint dummypt = QPoint(1 << 30, 1 << 30);

class Links {
public:
    Links() {
        p[0] = dummypt;
        p[1] = dummypt;
        p[2] = dummypt;
        p[3] = dummypt;
    }
    void add(QPoint r) {
        for (int i = 0; i < 4; i++) {
            if (p[i] == r) {
                return;
            }
            if (p[i] == dummypt) {
                p[i] = r;
                return;
            }
        }
        qFatal("Overflow on link");
    }
    int count() const {
        for (int i = 0; i < 4; i++) {
            if (p[i] == dummypt) {
                return i;
            }
        }
        return 4;
    }
    bool contains(QPoint q) const {
        for (int i = 0; i < 4; i++) {
            if (p[i] == q) {
                return true;
            }
        }
        return false;
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
};
// constexpr QPoint Links::dummy = QPoint(1 << 30, 1 << 30);

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

    /* Graph-style partitioning */
    qDebug("Identifying local segments");
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
                    ncls[i] = grid_points[corners[i].x() * H + corners[i].y()]
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
                            link_add(link_map, splits[(i + 3) % 4], splits[i]);
                        }
                    }
                } else if (nlow == 3) {
                    // Corner configuration, 1 high
                    for (int i = 0; i < 4; i++) {
                        if (!islow[i]) {
                            link_add(link_map, splits[(i + 3) % 4], splits[i]);
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
                            grid_points[corners[0].x() * H + corners[0].y()];

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
    qDebug("Refining %d grid-aligned dividers", eooe_pts.size());
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

    qDebug("Introducing corners where appropriate");
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

            if (pavg.region_class == pjump.region_class ||
                typematch(pjump, pavg)) {
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
        Links &l0 = link_map[vls.first];
        Links &l1 = link_map[vls.second];
        l0.replace_first(vls.second, inserted);
        l1.replace_first(vls.first, inserted);
        link_map[inserted].add(vls.first);
        link_map[inserted].add(vls.second);
    }

    qDebug("Refining %d star joints", oo_multi_pts.size());
    QMap<QPoint, QVector<RenderPoint>> oo_poly;
    for (const QPoint &combo : oo_multi_pts) {
        const Links &l = link_map[combo];
        // iterate CCW and check ...

        QVector<QPoint> pts;
        const QPoint deltas[4] = {QPoint(0, 1), QPoint(1, 0), QPoint(0, -1),
                                  QPoint(-1, 0)};
        for (int i = 0; i < 4; i++) {
            if (l.contains(combo + deltas[i])) {
                pts.append(combo + deltas[i]);
            }
        }
        int m = pts.size();
        if (m <= 2) {
            qDebug("Need at least three points for poly");
            continue;
        }

        // order [A B|B C|C D|D A], with pairs close by
        QVector<const RenderPoint *> doubled_pts;
        const QPair<RenderPoint, RenderPoint> *last =
            &eo_refinement[pts.last()];
        for (int i = 0; i < m; i++) {
            const QPair<RenderPoint, RenderPoint> *pair =
                &eo_refinement[pts[i]];
            if (pair->first.region_class == last->first.region_class ||
                pair->first.region_class == last->second.region_class) {
                doubled_pts.append(&pair->first);
                doubled_pts.append(&pair->second);
            } else {
                doubled_pts.append(&pair->second);
                doubled_pts.append(&pair->first);
            }
            last = pair;
        }
        // take pairwise averages to form a polytope
        QVector<RenderPoint> corners(m);
        int invalid = 0;
        for (int i = 0; i < m; i += 1) {
            int i2 = 2 * i;
            int j2 = (2 * i + 2 * m - 1) % (2 * m);
            const RenderPoint *p0 = doubled_pts[i2];
            const RenderPoint *p1 = doubled_pts[j2];
            QPointF qavg = (p0->coords + p1->coords) / 2;
            corners[i] = queryPoint(qavg, nav);
            RenderPoint &ravg = corners[i];
            if (typematch(ravg, *p0) && typematch(ravg, *p1)) {
                ravg.region_class = p0->region_class;
            } else {
                invalid++;
            }
        }

        qDebug("(%d,%d) with %d links, %d merged %svalid", combo.x(), combo.y(),
               l.count(), m - invalid, invalid ? "in" : "");
        if (invalid) {
            continue;
        }

        // polytope minimization.
        const int nrounds = 200;
        int i = 0;
        for (int r = 0; r < nrounds; r++) {
            bool adv = false;
            for (int k = 0; k < m; k++) {
                int ip = i;
                i = (i + 1) % m;
                // Test a new point, and modify the corner list if
                // successful
                int in = (i + 1) % m;
                // QPointF mid = (corners[ip].coords + corners[in].coords) /
                // 2;
                QPointF mid =
                    (rand() & 1) ? corners[ip].coords : corners[in].coords;
                // TODO: approaching mean is *bad*, leads to degeneracy.
                // also, keep getting area sign flips this way, which is
                // *odd*

                // query slightly closer
                RenderPoint alt =
                    queryPoint((mid + corners[i].coords) / 2, nav);
                if (typematch(alt, corners[i])) {
                    alt.region_class = corners[i].region_class;
                    corners[i] = alt;
                    adv = true;
                    break;
                }
            }
            if (!adv) {
                continue;
            }
        }
        oo_poly[combo] = corners;
        qDebug("area of %g", signed_area(corners));
    }

    qDebug("Extracting loops");
    // Every loop contains at least one E/O or O/E point, and those
    // have exactly two parents.
    QVector<QSet<QPoint>> class_points(grid_nclasses);
    for (QPoint p : eooe_pts) {
        QPoint n, f;
        partition_sumpoint(p, &n, &f);

        if (n.x() >= 0 && n.x() < W && n.y() >= 0 && n.y() < H) {
            class_points[grid_points[n.x() * H + n.y()].region_class].insert(p);
        }
        if (f.x() >= 0 && f.x() < W && f.y() >= 0 && f.y() < H) {
            class_points[grid_points[f.x() * H + f.y()].region_class].insert(p);
        }
    }

    for (int cls = 0; cls < grid_nclasses; cls++) {
        // Construct loops from class point seeds
        QSet<QPoint> &pts = class_points[cls];
        QVector<QVector<QPoint>> loops;
        while (pts.size()) {
            QVector<QPoint> loop;
            const QPoint initial = *pts.begin();
            QPoint cur = initial;
            loop.append(cur);
            pts.remove(cur);
            QPoint prev = dummypt;
            while (true) {
                Links l = link_map[cur];
                bool at_even = !point_is_oo(cur);
                if (at_even) {
                    if (l.count() != 2) {
                        qDebug("Require that even points have exactly two "
                               "neighbors");
                    }
                    // Move to complementary node, possible even or odd
                    if (prev != dummypt && l.p[0] != prev && l.p[1] != prev) {
                        qFatal("invariant fail");
                    }
                    QPoint next = l.p[0] == prev ? l.p[1] : l.p[0];
                    if (next == loop[0]) {
                        // Loop has closed
                        break;
                    }
                    prev = cur;
                    cur = next;
                    loop.append(cur);
                    pts.remove(cur);
                } else {
                    // Move to a node which we haven't seen yet; it will be even
                    QPoint best = dummypt;
                    for (int i = 0; i < l.count(); i++) {
                        QPoint a = l.p[i];
                        if (a == prev) {
                            // backtracking forbidden
                            continue;
                        }
                        // todo: orientation condition
                        if (pts.contains(a)) {
                            best = a;
                            break;
                        }
                    }
                    if (best == dummypt) {
                        break;
                    }

                    prev = cur;
                    cur = best;
                    loop.append(cur);
                    pts.remove(cur);
                }
            }
            if (link_map[initial].p[0] != cur &&
                link_map[initial].p[1] != cur) {
                qWarning("Loop failed to close: (%d,%d) and (%d,%d) not linked",
                         initial.x(), initial.y(), cur.x(), cur.y());
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
            int lps = loops[i].size();
            for (int j = 0; j < lps; j++) {
                const QPoint &p = loops[i][j];
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
                    } else if (oo_poly.contains(p)) {
                        // This is a joint!
                        const QVector<RenderPoint> &rpcyc = oo_poly[p];
                        const QPoint &nxt = loops[i][(j + 1) % lps];
                        const QPoint &prv = loops[i][(j + lps - 1) % lps];
                        const QPair<RenderPoint, RenderPoint> &np =
                            eo_refinement[nxt];
                        const QPair<RenderPoint, RenderPoint> &pp =
                            eo_refinement[prv];
                        QPointF nc = np.first.region_class == cls
                                         ? np.first.coords
                                         : np.second.coords;
                        QPointF pc = pp.first.region_class == cls
                                         ? pp.first.coords
                                         : pp.second.coords;

                        QPointF avg = (nc + pc) / 2;

                        // Pick the closest point of our class to the
                        // surrounding segment average
                        int closest = -1;
                        qreal distance = std::numeric_limits<qreal>::infinity();
                        for (int k = 0; k < rpcyc.size(); k++) {
                            QPointF delta = (rpcyc[k].coords - avg);
                            double ad =
                                delta.x() * delta.x() + delta.y() * delta.y();
                            if (ad < distance && rpcyc[k].region_class == cls) {
                                closest = k;
                                distance = ad;
                            }
                        }
                        if (closest < 0) {
                            qDebug("Class not encountered in "
                                   "multilink");
                            continue;
                        }
                        rp.append(rpcyc[closest]);
                    } else {
                        // failed to create a joint here
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

        region_list.append(reg);
    }

    // Set bounds
    for (int x = 0; x < W; x++) {
        for (int y = 0; y < H; y++) {
            int idx = x * H + y;
            int cls = grid_points[idx].region_class;
            region_list[cls].xmin = std::min(x, region_list[cls].xmin);
            region_list[cls].xmax = std::max(x, region_list[cls].xmax);
            region_list[cls].ymin = std::min(y, region_list[cls].ymin);
            region_list[cls].ymax = std::max(y, region_list[cls].ymax);
        }
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
            subreg.linear_angle = std::numeric_limits<float>::quiet_NaN();
            subreg.radial_center =
                QPointF(std::numeric_limits<float>::quiet_NaN(),
                        std::numeric_limits<float>::quiet_NaN());
            subreg.gradient.min = std::numeric_limits<float>::infinity();
            subreg.gradient.max = -std::numeric_limits<float>::infinity();
            subreg.gradient.colors.clear();

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

    int S = std::max(4, 1024 / std::max(W - 1, H - 1));
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
    return QStringLiteral("#%1%2%3%4%5%6")
        .arg(numbers[r / 16])
        .arg(numbers[r % 16])
        .arg(numbers[g / 16])
        .arg(numbers[g % 16])
        .arg(numbers[b / 16])
        .arg(numbers[b % 16]);
}

static double compute_cosangle(const QPointF &a, const QPointF &corner,
                               const QPointF &b) {
    const QPointF &dc1 = corner - a;
    const QPointF &dc2 = corner - b;
    qreal len01 = dc1.x() * dc1.x() + dc1.y() * dc1.y();
    qreal len12 = dc2.x() * dc2.x() + dc2.y() * dc2.y();
    qreal lenpro = std::sqrt(len01 * len12);
    if (lenpro <= 0.) {
        // assume perfectly aligned
        return -1.;
    }

    return QPointF::dotProduct(dc1, dc2) / lenpro;
}

static void bezier_cubic_basis(double t, double *b0, double *b1, double *b2,
                               double *b3) {
    *b0 = (1. - t) * (1. - t) * (1. - t);
    *b1 = 3 * (1. - t) * (1. - t) * t;
    *b2 = 3 * (1. - t) * t * t;
    *b3 = t * t * t;
}
static double det_2d(double m11, double m12, //
                     double m21, double m22) {
    return m11 * m22 - m12 * m21;
}
static void solve_2d_lin(double m11, double m12, double v1, //
                         double m21, double m22, double v2, //
                         double *y1, double *y2) {
    // cramer's rule
    double base = det_2d(m11, m12, m21, m22);
    *y1 = det_2d(v1, m12, v2, m22) / base;
    *y2 = det_2d(m11, v1, m21, v2) / base;
}
static double dist_2d(const QPointF &p, const QPointF &q) {
    return std::sqrt(QPointF::dotProduct(p - q, p - q));
}

static QStringList svg_subpath_text(const QVector<QPointF> &seq, int estart,
                                    int iend, double acceptable_error) {
    if (iend <= estart) {
        qFatal("Invalid range: estart=%d !< iend=%d", estart, iend);
    }
    const int fprec = 7;
    const QPointF &pstart = seq[estart], &pend = seq[iend];
    if (iend == estart + 1) {
        QStringList lst;
        lst.append(QStringLiteral("L%1,%2")
                       .arg(pend.x(), 0, 'f', fprec)
                       .arg(pend.y(), 0, 'f', fprec));
        return lst;
    }

    // Attempt linear reduction
    {
        double max_deviation = 0.;
        const QPointF &delta = pend - pstart;
        for (int i = estart + 1; i < iend - 1; i++) {
            double num = (delta.y() * seq[i].x() - delta.x() * seq[i].y() +
                          pend.x() * pstart.y() - pend.y() * pstart.x());
            max_deviation = std::max(max_deviation, std::abs(num));
        }
        // normalize to line scale
        max_deviation /= std::sqrt(QPointF::dotProduct(delta, delta));

        if (max_deviation < acceptable_error) {
            // return linear motion if cost is low enough
            QStringList lst;
            lst.append(QStringLiteral("L%1,%2")
                           .arg(pend.x(), 0, 'f', fprec)
                           .arg(pend.y(), 0, 'f', fprec));
            return lst;
        }
    }

    // Attempt bezier reduction if there are enough intermediate points
    if (iend >= estart + 3) {
        double xctrl[4], yctrl[4];
        assert(sizeof(QPointF) == 2 * sizeof(double));
        int n = iend - estart + 1;
        double err = bezier_segment_fit((double *)(seq.data() + estart), n,
                                        xctrl, yctrl);
        if (err <= 2 * acceptable_error) {
            QStringList lst;
            lst.append(QStringLiteral("C%1,%2,%3,%4,%5,%6")
                           .arg(xctrl[1], 0, 'f', fprec)
                           .arg(yctrl[1], 0, 'f', fprec)
                           .arg(xctrl[2], 0, 'f', fprec)
                           .arg(yctrl[2], 0, 'f', fprec)
                           .arg(pend.x(), 0, 'f', fprec)
                           .arg(pend.y(), 0, 'f', fprec));
            return lst;
        }
    }

    // Split by sharpest point, unless they are all *very* flat, in which
    // case divide in two
    int isharpest = (iend + estart) / 2;
    double cosa = -1. + 1e-3;
    for (int i = estart + 1; i < iend - 1; i++) {
        double ca = compute_cosangle(seq[i - 1], seq[i], seq[i + 1]);
        if (ca > cosa) {
            isharpest = i;
            cosa = ca;
        }
    }

    return svg_subpath_text(seq, estart, isharpest, acceptable_error) +
           svg_subpath_text(seq, isharpest, iend, acceptable_error);
}

static QStringList svg_path_text_from_curve(const QPolygonF &pts, bool closed) {
    const int fprec = 7;
    QStringList strl;
    QPointF s = pts[0];
    strl.append(QStringLiteral("M%1,%2")
                    .arg(s.x(), 0, 'f', fprec)
                    .arg(s.y(), 0, 'f', fprec));
    if (pts.size() > 1) {
        strl << svg_subpath_text(pts, 0, pts.size() - 1, 5e-6);
    }
    if (closed) {
        strl.append("Z");
    }
    return strl;
}

static QString svg_path_text_from_polygons(const QVector<QPolygonF> &loops,
                                           bool close_loops = true) {
    QStringList path_string;
    for (const QPolygonF &poly : loops) {
        path_string << svg_path_text_from_curve(poly, close_loops);
    }
    return path_string.join(" ");
}

static QPointF l1_min(const QVector<QPointF> pts) {
    QPointF current;
    for (const QPointF &p : pts) {
        current += p;
    }
    current /= pts.size();

    for (int iter = 0; iter < 50; iter++) {
        int nskipped = 0;
        double net_weight = 0.;
        QPointF net_contrib;
        for (const QPointF &p : pts) {
            if (p == current) {
                nskipped++;
            } else {
                double weight = 1. / dist_2d(p, current);
                net_weight += weight;
                net_contrib += p * weight;
            }
        }
        if (net_weight == 0. || nskipped == pts.size()) {
            // eqvtly, all points are identical
            break;
        }
        QPointF ideal = net_contrib / net_weight;

        if (nskipped > 0) {
            QPointF remnant = net_contrib - current * net_weight;
            double weight =
                std::min(1., nskipped / dist_2d(QPointF(), remnant));
            ideal = (1 - weight) * ideal + weight * current;
        }

        current = ideal;

        double sl1 = 0.;
        for (const QPointF &p : pts) {
            sl1 += dist_2d(p, current);
        }
    }
    return current;
}

static bool sort_by_first(const QPair<double, double> &a,
                          const QPair<double, double> &b) {
    return a.first < b.first;
}

static QVector<QPair<double, double>>
medfilt(const QVector<QPair<double, double>> &v, int m) {
    int d = m / 2;
    QVector<QPair<double, double>> o;
    QVector<double> t(m);
    for (int i = 0; i < v.size() - m; i++) {
        for (int j = 0; j < m; j++) {
            t[j] = v[i + j].second;
        }
        qSort(t);
        o.append(QPair<double, double>(v[i + d].first, t[d]));
    }
    return o;
}

/* Compute candidates for the center of a radial distribution. This must
 * be done in some sort of global fashion, because aggregation of local
 * estimates (i.e, from circle centers predicted by four points) amplifies
 * gradient errors due to the grid spacing. */
static QVector<QPointF> compute_radial_centers(const Region &region,
                                               const Subregion &subreg,
                                               RenderPoint *grid_points,
                                               const QSize &grid_size) {
    const int W = grid_size.width(), H = grid_size.height();
    QVector<QPair<QPointF, QPointF>> gradients;
    // Step 1: compute all gradients
    for (int x = region.xmin; x <= region.xmax; x++) {
        for (int y = region.ymin; y <= region.ymax; y++) {
            // There are four dx/dy kernels feasible;
            // we add them all together
            RenderPoint &base = grid_points[x * H + y];
            if (base.region_class != region.class_no)
                continue;
            if (base.subregion_class != subreg.subclass_no)
                continue;

            double mag_base = base.ideal_color.magnitude();

            QPointF net_deriv;
            QPoint dxs[4] = {QPoint(1, 1), QPoint(1, -1), QPoint(-1, -1),
                             QPoint(-1, 1)};
            for (int k = 0; k < 4; k++) {
                int xnx = x + dxs[k].x(), yny = y + dxs[k].y();
                if (xnx >= W || yny >= H || xnx < 0 || yny < 0) {
                    continue;
                }

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

                float mag_dx = pdx.ideal_color.magnitude();
                float mag_dy = pdy.ideal_color.magnitude();

                double dx =
                    (mag_dx - mag_base) / (pdx.coords.x() - base.coords.x());
                double dy =
                    (mag_dy - mag_base) / (pdy.coords.y() - base.coords.y());
                net_deriv += QPointF(dx, dy);
            }
            if (net_deriv.x() != 0 && net_deriv.y() != 0) {
                net_deriv /= dist_2d(QPointF(), net_deriv);
                gradients.append(
                    QPair<QPointF, QPointF>(base.coords, net_deriv));
            }
        }
    }
    if (gradients.size() < 5) {
        return QVector<QPointF>();
    }

    // Step 2: sample over angles
    const int n_angles =
        12; /* effectively controls minimum detectable slice angle */
    double offsets[n_angles];
    for (int i = 0; i < n_angles; i++) {
        double angle = i * CLHEP::pi / n_angles;
        double sa = std::sin(angle), ca = std::cos(angle);
        QVector<QPair<double, double>> tcp;
        for (const QPair<QPointF, QPointF> &p : gradients) {
            double t = -sa * p.first.x() + ca * p.first.y();
            double contra = std::abs(ca * p.second.x() + sa * p.second.y());
            tcp.push_back(QPair<double, double>(t, contra));
        }
        qSort(tcp.begin(), tcp.end(), sort_by_first);

        // median filter over cb
        int nmedlength = 5;
        QVector<QPair<double, double>> tcpfilt = medfilt(tcp, nmedlength);

        // select average position of peak (medfilt -> exist multiple)
        int nbest = 0;
        double best = 0.;
        double mxs = -1.;
        for (int k = 0; k < tcpfilt.size(); k++) {
            if (tcpfilt[k].second > mxs) {
                mxs = tcpfilt[k].second;
                best = tcpfilt[k].first;
                nbest = 1;
            } else if (tcpfilt[k].second == mxs) {
                best += tcpfilt[k].first;
                nbest += 1;
            }
        }
        offsets[i] = best / nbest;
    }

    // Step 3: compute pairwise intersections
    QVector<QPointF> ints;
    for (int i = 0; i < n_angles; i++) {
        int k = (i + 1) % n_angles;
        double angle_i = i * CLHEP::pi / n_angles;
        double angle_k = k * CLHEP::pi / n_angles;
        double si = std::sin(angle_i), ci = std::cos(angle_i);
        double sk = std::sin(angle_k), ck = std::cos(angle_k);

        QPointF pi(-si * offsets[i], ci * offsets[i]);
        QPointF pk(-sk * offsets[k], ck * offsets[k]);
        QPointF di(ci, si);
        QPointF dk(ck, sk);

        double s, t;
        solve_2d_lin(di.x(), dk.x(), pk.x() - pi.x(), di.y(), dk.y(),
                     pk.y() - pi.y(), &s, &t);
        QPointF crossi = pi + s * di;
        QPointF crossk = pk - t * dk;
        QPointF cx = (crossi + crossk) / 2;
        ints.append(cx);
    }

    // Add the mean point, for better high-symmetry perf
    ints.push_back(l1_min(ints));
    return ints;
}

static float compute_histogram_angle(const Region &region,
                                     const Subregion &subreg,
                                     RenderPoint *grid_points,
                                     const QSize &grid_size) {
    // TODO: ensure subregion bounding & point lookup is fast & local
    // TODO: figure out how to do this given an unstructured set of points

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

static int square(int s) { return s * s; }
static int rgba_distance2(QRgb a, QRgb b) {
    return square(qRed(a) - qRed(b)) + square(qGreen(a) - qGreen(b)) +
           square(qBlue(a) - qBlue(b)) + square(qAlpha(a) - qAlpha(b));
}

static Gradient
compute_gradient(const QVector<QPair<FColor, float>> &colors_at_position,
                 float min_spacing, long *error_post_quant) {
    Gradient g;
    g.min = std::numeric_limits<float>::infinity();
    g.max = -std::numeric_limits<float>::infinity();
    for (const QPair<FColor, float> &p : colors_at_position) {
        g.min = std::min(g.min, p.second);
        g.max = std::max(g.max, p.second);
    }
    int nsteps =
        std::max(2, 1 + (int)std::floor((g.max - g.min) / min_spacing));

    FColor *histogram_c = new FColor[nsteps]();
    double *histogram_w = new double[nsteps]();
    for (const QPair<FColor, float> &p : colors_at_position) {
        double cell = (nsteps - 1.) * (p.second - g.min) / (g.max - g.min);
        int lidx = (int)std::floor(cell);
        lidx = std::max(0, std::min(nsteps - 2, lidx));
        double frac = std::max(0.0, std::min(1.0, cell - lidx));
        histogram_c[lidx] =
            FColor::add(histogram_c[lidx], p.first, (1. - frac));
        histogram_c[lidx + 1] =
            FColor::add(histogram_c[lidx + 1], p.first, frac);
        histogram_w[lidx] += 1. - frac;
        histogram_w[lidx + 1] += frac;
    }
    QVector<QPair<int, FColor>> known_stops;
    for (int i = 0; i < nsteps; i++) {
        if (histogram_w[i] > 0.) {
            double iw = 1. / histogram_w[i];
            double rv =
                std::max(0.0, std::min(1.0, iw * histogram_c[i].redF()));
            double gv =
                std::max(0.0, std::min(1.0, iw * histogram_c[i].greenF()));
            double bv =
                std::max(0.0, std::min(1.0, iw * histogram_c[i].blueF()));
            double av =
                std::max(0.0, std::min(1.0, iw * histogram_c[i].alphaF()));
            known_stops.push_back(
                QPair<int, FColor>(i, FColor(rv, gv, bv, av)));
        }
    }
    delete[] histogram_c;
    delete[] histogram_w;

    int j = 0;
    for (int i = 0; i < nsteps; i++) {
        if (j == 0 && known_stops[j].first > i) {
            g.colors.push_back(known_stops[j].second.rgba());
        } else if (j == known_stops.size() - 1) {
            g.colors.push_back(known_stops[j].second.rgba());
        } else {
            const QPair<int, FColor> &left = known_stops[j];
            const QPair<int, FColor> &right = known_stops[j + 1];
            FColor mix = FColor::blend(left.second, right.second,
                                       (i - left.first) /
                                           (double)(right.first - left.first));
            g.colors.push_back(mix.rgbaRound());
            if (right.first > i) {
                j++;
            }
        }
    }
    *error_post_quant = 0;
    for (const QPair<FColor, float> &p : colors_at_position) {
        double cell = (nsteps - 1.) * (p.second - g.min) / (g.max - g.min);
        int lidx = (int)std::floor(cell);
        lidx = std::max(0, std::min(nsteps - 2, lidx));
        double frac = std::max(0.0, std::min(1.0, cell - lidx));
        const FColor &prediction = FColor::blend(
            FColor(g.colors[lidx]), FColor(g.colors[lidx + 1]), frac);
        *error_post_quant +=
            rgba_distance2(p.first.rgbaRound(), prediction.rgbaRound());
    }
    return g;
}

static Gradient
compute_linear_gradient(const QVector<QPair<FColor, QPointF>> &interior_colors,
                        const QSize &gridsize, double angle, long *error) {
    if (std::isnan(angle)) {
        *error = std::numeric_limits<long>::max();
        Gradient g;
        return g;
    }
    int S = std::max(gridsize.width(), gridsize.height()) - 1;
    double grid_spacing = 1.0 / S;
    QVector<QPair<FColor, float>> colors_at_position;
    for (const QPair<FColor, QPointF> &p : interior_colors) {
        colors_at_position.push_back(
            QPair<FColor, float>(p.first, proj_angle(p.second, angle)));
    }
    return compute_gradient(colors_at_position, grid_spacing, error);
}

static Gradient
compute_radial_gradient(const QVector<QPair<FColor, QPointF>> &interior_colors,
                        const QSize &gridsize, const QPointF &center,
                        long *error) {
    if (std::isnan(center.x()) || std::isnan(center.y())) {
        *error = std::numeric_limits<long>::max();
        Gradient g;
        return g;
    }
    int S = std::max(gridsize.width(), gridsize.height()) - 1;
    double grid_spacing = 1.0 / S;
    QVector<QPair<FColor, float>> colors_at_position;
    for (const QPair<FColor, QPointF> &p : interior_colors) {
        colors_at_position.push_back(
            QPair<FColor, float>(p.first, dist_2d(p.second, center)));
    }
    return compute_gradient(colors_at_position, grid_spacing, error);
}

static long
compute_solid_error(FColor unif_pred,
                    const QVector<QPair<FColor, QPointF>> &interior_colors) {
    long sqe = 0.;
    for (const QPair<FColor, QPointF> &p : interior_colors) {
        sqe += rgba_distance2(unif_pred.rgbaRound(), p.first.rgbaRound());
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

static QVector<QRgb> hatch_dim(const QVector<QRgb> &a) {
    QVector<QRgb> r(a.size());
    for (int i = 0; i < a.size(); i++) {
        r[i] = qRgba(4 * qRed(a[i]) / 5, 4 * qGreen(a[i]) / 5,
                     4 * qBlue(a[i]) / 5, qAlpha(a[i]));
    }
    return r;
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

static void add_gradient_stops(QTextStream &s, const QVector<QRgb> &v,
                               float minv, float maxv) {
    if (v.size() < 2) {
        qFatal("Invalid color vector, length %d", v.size());
    }

    int nsteps = v.size();
    for (int i = 0; i < nsteps; i++) {
        double stop_pos = (maxv - minv) * i / (nsteps - 1.) + minv;

        if (i > 0 && i < nsteps - 1 && v[i] == v[i - 1]) {
            // Skip redundant points
            continue;
        }
        if (qAlpha(v[i]) < 255) {
            double alpha = qAlpha(v[i]) / 255.;
            s << QStringLiteral("  <stop offset=\"%1\" stop-color=\"%2\" "
                                "stop-opacity=\"%3\"/>\n")
                     .arg(stop_pos, 0, 'f', 5)
                     .arg(color_hex_name_rgb(v[i]))
                     .arg(alpha, 0, 'f', 3);
        } else {
            // SVG specifies a default stop-opacity of 1
            s << QStringLiteral("  <stop offset=\"%1\" stop-color=\"%2\"/>\n")
                     .arg(stop_pos, 0, 'f', 5)
                     .arg(color_hex_name_rgb(v[i]));
        }
    }
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
            FColor meanc = FColor(net_r / net_count, net_g / net_count,
                                  net_b / net_count, net_a / net_count);
            subreg.solid_color = meanc.rgbaRound();
            long best_err = compute_solid_error(meanc, interior_colors);
            subreg.gradient_type = GradientType::gSolid;
            if (best_err > 0) {
                long this_err = std::numeric_limits<long>::max();
                subreg.linear_angle = compute_histogram_angle(
                    region, subreg, grid_points, grid_size);
                Gradient g = compute_linear_gradient(
                    interior_colors, grid_size, subreg.linear_angle, &this_err);
                if (this_err < best_err) {
                    subreg.gradient_type = GradientType::gLinear;
                    subreg.gradient = g;
                    best_err = this_err;
                }
            }
            if (best_err > 0) {
                long this_err = std::numeric_limits<long>::max();
                for (const QPointF &c : compute_radial_centers(
                         region, subreg, grid_points, grid_size)) {
                    Gradient g = compute_radial_gradient(
                        interior_colors, grid_size, c, &this_err);
                    if (this_err < best_err) {
                        subreg.radial_center = c;
                        subreg.gradient_type = GradientType::gRadial;
                        subreg.gradient = g;
                        best_err = this_err;
                    }
                }
            }

            const char *type =
                subreg.gradient_type == GradientType::gSolid
                    ? "sol"
                    : (subreg.gradient_type == GradientType::gRadial ? "rad"
                                                                     : "lin");
            qDebug("class %d subclass %d type=%s nsteps %d angle %f center "
                   "(%f,%f)",
                   region.class_no, subreg.subclass_no, type,
                   subreg.gradient.colors.size(), subreg.linear_angle,
                   subreg.radial_center.x(), subreg.radial_center.y());
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

    QTemporaryFile tfile(QStringLiteral("/tmp/render.svg.XXXXXX"));

    tfile.open();
    qDebug("Rendering final image, %d regions, to %s", region_list.size(),
           tfile.fileName().toUtf8().constData());
    {
        QTextStream s(&tfile);

        s << QStringLiteral(
            "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n "
            "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");

        s << QStringLiteral("<svg width=\"%1\" height=\"%2\" ersion=\"1.1\" "
                            "xmlns=\"http://www.w3.org/2000/svg\" "
                            "xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n")
                 .arg(viewbox.width())
                 .arg(viewbox.height());

        s << QStringLiteral("<defs>\n");

        // Hatching fill textures
        QMap<CompactNormal, int> normal_directions;
        for (const Region &region : region_list) {
            for (const Subregion &subreg : region.subregions) {
                // consider average normal over subregion
                if (subreg.is_clipped_patch) {
                    QPoint src = subreg.representative_coord;
                    const RenderPoint &r = grid_points[src.x() * H + src.y()];
                    const CompactNormal &normal = r.ray.intersections[0].normal;
                    if (!normal_directions.contains(normal)) {
                        normal_directions[normal] = normal_directions.size();
                    }
                }
            }
        }
        const G4ThreeVector &updown = view_data.orientation.rowY();
        const G4ThreeVector &lright = view_data.orientation.rowZ();
        for (const CompactNormal &nd : normal_directions.keys()) {
            int idx = normal_directions[nd];
            const G4ThreeVector &normal = nd;

            const G4ThreeVector &orthA = normal.orthogonal().unit();
            const G4ThreeVector &orthB = normal.cross(orthA).unit();
            const G4ThreeVector &pattern_dir = (orthA + orthB);
            // Use unit viewport space
            double udf = pattern_dir.dot(updown), lrf = pattern_dir.dot(lright);
            // Compute spacing/angle
            double spacing = 0.0175 * std::sqrt(udf * udf + lrf * lrf);
            // Viewed perpendicularly, the stripe angle should not be pi/4
            double angle = std::atan2(udf, lrf) - M_PI / 30;

            // we do not want infinite density
            spacing = std::max(spacing, 0.001);

            // We produce only one clip path, for the darker strips. This
            // will overlay the lighter region
            s << QStringLiteral(" <clipPath id=\"hatching%1\">").arg(idx);

            double L = std::sqrt(2);
            int n = std::ceil(L / spacing / 2);
            QPointF prim(std::cos(angle), std::sin(angle));
            QPointF sec(-std::sin(angle), std::cos(angle));
            prim *= spacing;
            sec *= L;
            QStringList lpath;
            for (int x = -n; x <= n; x++) {
                QPointF corners[4] = {
                    prim * (2 * x - 0.5) + sec, prim * (2 * x - 0.5) - sec,
                    prim * (2 * x + 0.5) - sec, prim * (2 * x + 0.5) + sec};
                for (int j = 0; j < 4; j++) {
                    QPointF wco = transf.map(corners[j]);
                    lpath.push_back(QStringLiteral("%1%2,%3")
                                        .arg("MLLL"[j])
                                        .arg(wco.x(), 0, 'g', 6)
                                        .arg(wco.y(), 0, 'g', 6));
                }
                lpath.push_back("Z");
            }

            s << QStringLiteral("  <path d=\"%1\"/>\n").arg(lpath.join(" "));
            s << QStringLiteral(" </clipPath>\n");
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
                    QPointF ssfstart = sscenter + (subreg.gradient.min - ssf) *
                                                      QPointF(-lsin, lcos);
                    QPointF ssfstop = sscenter + (subreg.gradient.max - ssf) *
                                                     QPointF(-lsin, lcos);

                    ssfstart = transf.map(ssfstart);
                    ssfstop = transf.map(ssfstop);

                    for (int v = 0; v < (subreg.is_clipped_patch ? 2 : 1);
                         v++) {
                        s << QStringLiteral(
                                 " <linearGradient id=\"gradient%1_%2%3\" "
                                 "x1=\"%4\" "
                                 "y1=\"%5\" x2=\"%6\" y2=\"%7\" "
                                 "spreadMethod=\"%8\" "
                                 "gradientUnits=\"userSpaceOnUse\">\n")
                                 .arg(region.class_no)
                                 .arg(subreg.subclass_no)
                                 .arg(subreg.is_clipped_patch
                                          ? (v ? "_dim" : "")
                                          : "")
                                 .arg(ssfstart.x())
                                 .arg(ssfstart.y())
                                 .arg(ssfstop.x())
                                 .arg(ssfstop.y())
                                 .arg(0 ? "repeat" : "pad");
                        bool dim = subreg.is_clipped_patch && v == 1;
                        add_gradient_stops(
                            s,
                            dim ? hatch_dim(subreg.gradient.colors)
                                : subreg.gradient.colors,
                            0., 1.);
                        s << QStringLiteral(" </linearGradient>\n");
                    }
                } else if (subreg.gradient_type == GradientType::gRadial) {
                    QPointF tc = transf.map(subreg.radial_center);
                    double rs = transf.m11() * subreg.gradient.max;

                    for (int v = 0; v < (subreg.is_clipped_patch ? 2 : 1);
                         v++) {
                        s << QStringLiteral(
                                 " <radialGradient id=\"gradient%1_%2%3\" "
                                 "cx=\"%4\" "
                                 "cy=\"%5\" r=\"%6\" "
                                 "spreadMethod=\"pad\" "
                                 "gradientUnits=\"userSpaceOnUse\">\n")
                                 .arg(region.class_no)
                                 .arg(subreg.subclass_no)
                                 .arg(subreg.is_clipped_patch
                                          ? (v ? "_dim" : "")
                                          : "")
                                 .arg(tc.x())
                                 .arg(tc.y())
                                 .arg(rs);
                        bool dim = subreg.is_clipped_patch && v == 1;
                        add_gradient_stops(
                            s,
                            dim ? hatch_dim(subreg.gradient.colors)
                                : subreg.gradient.colors,
                            subreg.gradient.min / subreg.gradient.max, 1.);
                        s << QStringLiteral(" </radialGradient>\n");
                    }
                }
            }
        }
        s << QStringLiteral("</defs>\n");

        // Interior regions, clipped by boundaries
        s << QStringLiteral(
            "<g id=\"interiors\" fill-opacity=\"1\" stroke=\"none\" "
            "fill-rule=\"evenodd\">\n");
        for (const Region &region : region_list) {
            // Compute region limit
            for (const Subregion &subreg : region.subregions) {
                const QVector<QPolygonF> &subloops =
                    boundary_loops_for_subregion(subreg, transf);

                const QString &btext = svg_path_text_from_polygons(subloops);

                QString fill_desc, fill_desc_dim;
                if (subreg.gradient_type == GradientType::gSolid) {
                    fill_desc = color_hex_name_rgb(subreg.solid_color);
                    QRgb c = qRgba(4 * qRed(subreg.solid_color) / 5,
                                   4 * qGreen(subreg.solid_color) / 5,
                                   4 * qBlue(subreg.solid_color) / 5,
                                   qAlpha(subreg.solid_color));
                    fill_desc_dim = color_hex_name_rgb(c);
                } else if (subreg.gradient_type == GradientType::gLinear ||
                           subreg.gradient_type == GradientType::gRadial) {
                    fill_desc = QStringLiteral("url(#gradient%1_%2)")
                                    .arg(region.class_no)
                                    .arg(subreg.subclass_no);
                    fill_desc_dim = QStringLiteral("url(#gradient%1_%2_dim)")
                                        .arg(region.class_no)
                                        .arg(subreg.subclass_no);
                } else {
                    qFatal("Unsupported gradient type");
                }

                s << QStringLiteral(
                         " <path id=\"region%1_%2\" fill=\"%3\" d=\"%4\"/>\n")
                         .arg(region.class_no)
                         .arg(subreg.subclass_no)
                         .arg(fill_desc)
                         .arg(btext);

                if (subreg.is_clipped_patch) {
                    QPoint src = subreg.representative_coord;
                    const RenderPoint &r = grid_points[src.x() * H + src.y()];
                    const CompactNormal &normal = r.ray.intersections[0].normal;
                    int idx = normal_directions[normal];
                    s << QStringLiteral(
                             " <path id=\"region%1_%2_hatch\" fill=\"%3\" "
                             "clip-path=\"url(#hatching%4)\" d=\"%5\"/>\n")
                             .arg(region.class_no)
                             .arg(subreg.subclass_no)
                             .arg(fill_desc_dim)
                             .arg(idx)
                             .arg(btext);
                }
            }
        }
        s << QStringLiteral("</g>\n");

        // Boundaries
        s << QStringLiteral(
                 "<g id=\"boundaries\" fill=\"none\" stroke-width=\"%1\" "
                 "stroke-linecap=\"square\" stroke-linejoin=\"miter\">\n")
                 .arg(T * 0.003);
        for (const Region &region : region_list) {
            const QVector<QPolygonF> &loops =
                boundary_loops_for_region(region, transf);
            for (int i = 0; i < loops.size(); i++) {
                QVector<QPolygonF> solo;
                solo.push_back(loops[i]);
                QRgb c = i ? region.meanInteriorColors[i - 1]
                           : region.meanExteriorColor;
                s << QStringLiteral(" <path id=\"edge%1_%2\" stroke=\"%3\" "
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

                        s << QStringLiteral(
                                 " <path id=\"edge%1_%2_%3\" stroke=\"%4\" "
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
                            s << QStringLiteral(
                                     " <path id=\"edge%1_%2_%3_%4\" "
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
        s << QStringLiteral("</g>\n");

        s << QStringLiteral("<g text-anchor=\"middle\" font-size=\"3.0px\" "
                            "stroke-width=\"%1\">\n")
                 .arg(T * 0.0025, 0, 'g');
        {
            double viewport_width =
                2 * view_data.scale * (W * 1. / std::min(W, H));
            const QPair<double, QString> &rule =
                ruler_distance(viewport_width * 0.3, 0.3 * T);

            s << QStringLiteral(" <path fill=\"none\" stroke=\"#000000\" "
                                "d=\"M%1,%4 L%1,%3 %2,%3 %2,%4\" "
                                "id=\"ruler_bracket\"/>\n")
                     .arg(0.05 * T, 0, 'f', 5)
                     .arg(0.05 * T + rule.first, 0, 'f', 5)
                     .arg(0.95 * T, 0, 'f', 5)
                     .arg(0.93 * T, 0, 'f', 5);

            s << QStringLiteral(" <text fill=\"#000000\" "
                                "x=\"%1\" "
                                "y=\"%2\" id=\"zero\">0</text>\n")
                     .arg(0.05 * T, 0, 'f', 5)
                     .arg(0.975 * T, 0, 'f', 5);
            s << QStringLiteral(" <text fill=\"#000000\" "
                                "x=\"%1\" "
                                "y=\"%2\" id=\"zero\">%3</text>\n")
                     .arg(0.05 * T + rule.first, 0, 'f', 5)
                     .arg(0.975 * T, 0, 'f', 5)
                     .arg(rule.second);
            s << QStringLiteral("</g>\n");
        }
        s << QStringLiteral("</svg>\n");
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
        QString fmt(QStringLiteral("vector%1.svg"));
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

        emit produceImagePhase(grad_image, QStringLiteral("Render completed"),
                               nqueries, true);
    }

    QTime time_end = QTime::currentTime();
    qDebug("Conversion took %f seconds", time_start.msecsTo(time_end) * 1e-3);
}
