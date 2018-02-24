#include "VectorTrace.hh"

#include <QFile>
#include <QPainter>
#include <QPen>
#include <QProcess>
#include <QSet>
#include <QTextStream>
#include <QTime>

#include <G4Material.hh>

static uint32_t randint(uint32_t excl_upper) {
    return ((uint32_t)qrand()) % excl_upper;
}

static FColor randColor() {
    return FColor(randint(65536) / 65535., randint(65536) / 65535.,
                  randint(65536) / 65535., 1.0);
}

FColor::FColor() : r(0.), g(0.), b(0.), a(0.) {}
FColor::FColor(float ir, float ig, float ib, float ia)
    : r(ir), g(ig), b(ib), a(ia) {}
QRgb FColor::rgba() const { return QColor::fromRgbF(r, g, b, a).rgba(); }
float FColor::magnitude() const { return a * std::max(r, std::max(b, g)); }
QString FColor::hexName() const {
    QRgb color = rgba();
    int cr = qRed(color);
    int cg = qGreen(color);
    int cb = qBlue(color);
    const char *numbers = "0123456789abcdef";
    return QString("#%1%2%3%4%5%6")
        .arg(numbers[cr / 16])
        .arg(numbers[cr % 16])
        .arg(numbers[cg / 16])
        .arg(numbers[cg % 16])
        .arg(numbers[cb / 16])
        .arg(numbers[cb % 16]);
}
FColor FColor::blend(const FColor &a, const FColor &b, float s) {
    float t = 1. - s;
    return FColor(t * a.r + s * b.r, t * a.g + s * b.g, t * a.b + s * b.b,
                  t * a.a + s * b.a);
}

RenderPoint::RenderPoint() {
    // Initialize to nan to break use when uninitialized
    double nan = std::numeric_limits<double>::quiet_NaN();
    coords = QPointF(nan, nan);
    intersections = NULL;
    elements = NULL;
    nhits = 0;
    ideal_color = FColor();
    region_class = -1;
}
RenderPoint::RenderPoint(QPointF spot, int inhits, const Intersection *srcints,
                         const Element **srcelems) {
    coords = spot;
    nhits = inhits;
    elements = new Element *[nhits];
    intersections = new Intersection[nhits + 1];
    memcpy(elements, srcelems, nhits * sizeof(Element *));
    memcpy(intersections, srcints, (nhits + 1) * sizeof(Intersection));
    ideal_color = FColor();
    region_class = -1;
    subregion_class = -1;
}
RenderPoint::~RenderPoint() {
    if (elements)
        delete[] elements;
    if (intersections)
        delete[] intersections;
}
RenderPoint::RenderPoint(const RenderPoint &other) {
    nhits = other.nhits;
    coords = other.coords;
    ideal_color = other.ideal_color;
    region_class = other.region_class;
    subregion_class = other.subregion_class;
    if (other.intersections) {
        intersections = new Intersection[nhits + 1];
        memcpy(intersections, other.intersections,
               (nhits + 1) * sizeof(Intersection));
    } else {
        intersections = NULL;
    }
    if (other.elements) {
        elements = new Element *[nhits];
        memcpy(elements, other.elements, nhits * sizeof(Element *));
    } else {
        elements = NULL;
    }
}
RenderPoint &RenderPoint::operator=(RenderPoint copy_of_other) {
    swap(copy_of_other);
    return *this;
}
void RenderPoint::swap(RenderPoint &other) {
    std::swap(intersections, other.intersections);
    std::swap(elements, other.elements);
    std::swap(nhits, other.nhits);
    std::swap(region_class, other.region_class);
    std::swap(subregion_class, other.subregion_class);
    std::swap(coords, other.coords);
    std::swap(ideal_color, other.ideal_color);
}

static void recsetColorsByMaterial(Element &elem, QVector<FColor> &color_table,
                                   QMap<QString, FColor> &color_map,
                                   QMap<QString, int> &idx_map) {
    // We hard-code color associations to be more consistent
    QString key(elem.material->GetName().c_str());
    if (!color_map.count(key)) {
        FColor c = randColor();
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
    for (Element &e : elem.children) {
        recsetColorsByMaterial(e, color_table, color_map, idx_map);
    }
}

VectorTracer::VectorTracer(ViewData vd, TrackData td,
                           const QString &target_file, bool transparency,
                           QObject *parent)
    : QObject(parent), view_data(vd), track_data(td) {
    step_next = Steps::sGrid;
    file_name = target_file;
    nqueries = 0;
    transparent_volumes = transparency;

    QMap<QString, FColor> color_map;
    color_map["ArGas"] = FColor(0.8, 0.8, 0.8);
    color_map["Argas"] = color_map["ArGas"];
    color_map["ArgonGas"] = color_map["ArgonGas"];
    color_map["Air"] = FColor(0.8, 0.8, 0.8);

    color_map["Al6061"] = FColor(1.0, 0.2, 0.2);
    color_map["Al5083"] = FColor(0.8, 0.4, 0.4);
    // Often have nickel plated lead
    color_map["G4_Pb"] = FColor(0.4, 0.0, 0.7);
    color_map["G4_Ni"] = FColor(0.5, 0.1, 0.8);

    color_map["G4_Cu"] = FColor(0.0, 0.9, 0.7);
    color_map["BaF2"] = FColor(1.0, 1.0, 1.0);
    QMap<QString, int> idx_map;
    element_colors.clear();
    recsetColorsByMaterial(view_data.elements, element_colors, color_map,
                           idx_map);

    grid_size = QSize(1000, 1000);
    //    grid_size = QSize(300, 300);
    //    grid_size = QSize(100, 100);
    //    grid_size = QSize(30, 30);

    grid_points = NULL;
    grid_nclasses = 0;
    ray_mutables = NULL;
    ray_iteration = 0;
}

VectorTracer::~VectorTracer() {
    if (grid_points)
        delete[] grid_points;
    if (ray_mutables)
        delete[] ray_mutables;
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
void VectorTracer::reset(bool transp) {
    step_next = Steps::sGrid;
    transparent_volumes = transp;
    qDebug("Reset tracer: %s mode", transp ? "transparent" : "opaque");
}
int VectorTracer::faildepth(const RenderPoint &a, const RenderPoint &b) {
    // Dummy points always match each other
    if (!a.intersections && !b.intersections)
        return -1;

    if (!transparent_volumes) {
        // Consider only the first disagreement in visible volumes
        // Q: should compressTraces automatically vanish invisible volumes?
        // and use a null element gap instead?
        if (!a.intersections || !b.intersections)
            return 0;
        Element *first_a = NULL, *first_b = NULL;
        Intersection *first_ia = NULL, *first_ib = NULL;
        for (int i = 0; i < a.nhits; i++) {
            if (a.elements[i]->visible) {
                first_a = a.elements[i];
                first_ia = &a.intersections[i];
                break;
            }
        }
        for (int i = 0; i < b.nhits; i++) {
            if (b.elements[i]->visible) {
                first_b = b.elements[i];
                first_ib = &b.intersections[i];
                break;
            }
        }
        if (first_a != first_b) {
            return 0;
        }
        if (first_ia && first_ib &&
            first_ia->is_clipping_plane != first_ib->is_clipping_plane) {
            return 0;
        }
        return -1;
    } else {
        for (int i = 0; i < std::min(a.nhits, b.nhits); i++) {
            if (a.elements[i] != b.elements[i] ||
                a.intersections[i].is_clipping_plane !=
                    b.intersections[i].is_clipping_plane) {
                return i;
            }
        }
        if (a.nhits != b.nhits) {
            return std::min(a.nhits, b.nhits);
        }
        if (!a.intersections || !b.intersections)
            return 0;

        int n = a.nhits; // also = b.nhits
        if (a.intersections[n].is_clipping_plane !=
            b.intersections[n].is_clipping_plane) {
            return n;
        }

        // no failure
        return -1;
    }
}

static Intersection *first_nontrivial_intersection(const RenderPoint &p) {
    for (int i = 0; i <= p.nhits; i++) {
        if ((i > 0 && p.elements[i - 1]->visible) ||
            (i < p.nhits && p.elements[i]->visible)) {
            return &p.intersections[i];
            break;
        }
    }
    return NULL;
}
static QPointF grid_coord_to_point(const QPoint &pt, const QSize &grid_size) {
    const int W = grid_size.width() - 1, H = grid_size.height() - 1;
    const int S = std::max(W, H);

    QPointF spot((pt.x() - 0.5 * W) / S, (pt.y() - 0.5 * H) / S);
    return spot;
}
static QPointF point_to_grid_coord(const QPointF &pt, const QSize &grid_size) {
    const int W = grid_size.width() - 1, H = grid_size.height() - 1;
    const int S = std::max(W, H);

    QPointF spot(pt.x() * S + 0.5 * W, pt.y() * S + 0.5 * H);
    return spot;
}

RenderPoint VectorTracer::queryPoint(QPointF spot) {
    nqueries++;

    QPointF bound_low = grid_coord_to_point(QPoint(0, 0), grid_size);
    QPointF bound_high = grid_coord_to_point(
        QPoint(grid_size.width() - 1, grid_size.height() - 1), grid_size);
    const double eps = 1e-10;
    if (spot.x() < bound_low.x() - eps || spot.y() < bound_low.y() - eps ||
        spot.x() > bound_high.x() + eps || spot.y() > bound_high.y() + eps) {
        // Spot out of bounds; return dummy point
        RenderPoint r;
        r.coords = spot;
        return r;
    }

    if (!ray_mutables) {
        int treedepth;
        int nelements;
        countTree(view_data.elements, treedepth, nelements);
        ray_mutables = new ElemMutables[nelements]();
        if (treedepth > 10) {
            qFatal("Excessive tree depth, fatal!");
        }
    }
    const int dlimit = 100;
    const Element *hits[dlimit];
    Intersection ints[dlimit + 1];

    int m = traceRay(spot, view_data, hits, ints, dlimit, ray_iteration,
                     ray_mutables);
    ray_iteration++;
    m = compressTraces(hits, ints, m);

    return RenderPoint(spot, m, ints, hits);
}
RenderPoint VectorTracer::getPoint(QPoint p) {
    if (p.x() < 0 || p.y() < 0 || p.x() >= grid_size.width() ||
        p.y() >= grid_size.height()) {
        RenderPoint q;
        q.coords = grid_coord_to_point(p, grid_size);
        return q;
    } else {
        return grid_points[p.x() * grid_size.height() + p.y()];
    }
}
void VectorTracer::bracketEdge(const RenderPoint &initial_inside,
                               const RenderPoint &initial_outside,
                               RenderPoint *result_inside,
                               RenderPoint *result_outside) {
    if (typematch(initial_inside, initial_outside)) {
        qFatal("Bracket must cross class boundary");
    }
    const int nsubdivisions = 20;
    *result_inside = initial_inside;
    *result_outside = initial_outside;
    for (int k = 0; k < nsubdivisions; k++) {
        QPointF p_mid = 0.5 * (result_inside->coords + result_outside->coords);
        RenderPoint test = queryPoint(p_mid);
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
                                 RenderPoint *result_outside) {
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
    Intersection cx_inside =
        transparent_volumes ? result_inside->intersections[cd]
                            : *first_nontrivial_intersection(*result_inside);
    Intersection cx_outside =
        transparent_volumes ? result_outside->intersections[cd]
                            : *first_nontrivial_intersection(*result_outside);
    const int nsubdivisions = 20;
    for (int i = 0; i < nsubdivisions; i++) {
        RenderPoint mid =
            queryPoint(0.5 * (result_inside->coords + result_outside->coords));
        Intersection cx_mid;
        if (transparent_volumes) {
            if (mid.nhits < cd) {
                // Size failure for midpoint; abort
                return;
            }
            cx_mid = mid.intersections[cd];
        } else {
            Intersection *k = first_nontrivial_intersection(mid);
            if (!k) {
                // Lack of visible zone for midpoint; abort
                return;
            }
            cx_mid = *k;
        }
        if (is_jump) {
            if (std::abs(cx_mid.dist - cx_inside.dist) <
                std::abs(cx_mid.dist - cx_outside.dist)) {
                *result_inside = mid;
            } else {
                *result_outside = mid;
            }
        } else {
            if (cx_mid.normal.dot(cx_inside.normal) >
                cx_mid.normal.dot(cx_outside.normal)) {
                *result_inside = mid;
            } else {
                *result_outside = mid;
            }
        }
    }
}
static FColor intersectionColor(const G4ThreeVector &normal,
                                const G4ThreeVector &forward,
                                const FColor &base) {
    // Opposed normals (i.e, for transp backsides) are mirrored
    float cx = std::abs(std::acos(-normal * forward) / CLHEP::pi);
    cx = 1.0 - std::max(0.0, std::min(1.0, 0.7 * cx));

    return FColor(cx * base.redF(), cx * base.greenF(), cx * base.blueF(), 1.0);
}

static FColor retractionMerge(const RenderPoint &pt, FColor seed,
                              const G4ThreeVector &normal,
                              const QVector<FColor> &color_table,
                              int startpoint, bool transparent_volumes) {
    for (int k = startpoint; k >= 0; --k) {
        if (!pt.elements[k]->visible) {
            continue;
        }

        // We use the intersection before the volume
        const FColor &base_color = color_table[pt.elements[k]->ccode];
        FColor alt_color =
            intersectionColor(pt.intersections[k].normal, normal, base_color);
        double e = transparent_volumes ? pt.elements[k]->alpha : 1.0;

        seed = FColor::blend(seed, alt_color, e);
    }
    return seed;
}
FColor VectorTracer::calculateInteriorColor(const RenderPoint &pt) {

    FColor color(0.9, 0.9, 0.9, 0.1);
    return retractionMerge(pt, color, view_data.orientation.rowX(),
                           element_colors, pt.nhits - 1, transparent_volumes);
}
FColor VectorTracer::calculateBoundaryColor(const RenderPoint &inside,
                                            const RenderPoint &outside) {
    // lim necessary <= inside.nhits
    int lim = faildepth(inside, outside);
    FColor color(0., 0., 0., 1.);
    if (!transparent_volumes)
        return color;
    // It is possible that lim=inside.nhits if outside has an extra solid
    return retractionMerge(inside, color, view_data.orientation.rowX(),
                           element_colors, lim - 1, true);
}

void VectorTracer::computeGrid() {
    const int W = grid_size.width(), H = grid_size.height();

    qDebug("Starting grid computation");
    if (grid_points)
        delete[] grid_points;
    grid_points = new RenderPoint[W * H];
    for (int i = 0; i < W; i++) {
        for (int j = 0; j < H; j++) {
            int idx = i * H + j;
            QPointF spot = grid_coord_to_point(QPoint(i, j), grid_size);
            grid_points[idx] = queryPoint(spot);
            grid_points[idx].ideal_color =
                calculateInteriorColor(grid_points[idx]);
        }
    }
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

void VectorTracer::computeEdges() {
    region_list.clear();
    int W = grid_size.width(), H = grid_size.height();
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
                    // 0001=1100 Corner (x,y)
                    segs.push_back(QLine(2 * x + 1, 2 * y, 2 * x, 2 * y + 1));
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
                        RenderPoint rp = queryPoint(pt);
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
                            segs.push_back(
                                QLine(2 * x + 1, 2 * y, 2 * x + 2, 2 * y + 1));
                            segs.push_back(
                                QLine(2 * x + 1, 2 * y + 2, 2 * x, 2 * y + 1));
                        }
                    }

                    break;

                default:
                    qFatal("Shouldn't happen");
                }
            }
        }

        qDebug("Construction region boundaries %d: %d segments", cls,
               segs.length());

        QVector<QPolygon> loops = size_sorted_loops_from_seg(segs);

        qDebug("Refining region boundaries %d: %d loops", cls, loops.size());

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
        region.is_clipped_patch = false;
        if (class_representative.intersections) {
            region.is_clipped_patch =
                class_representative.intersections[0].is_clipping_plane;
        }
        if (class_representative.nhits > 0) {
            if (!class_representative.elements[0]->visible)
                region.is_clipped_patch = false;
        }

        for (int i = 0; i < loops.size(); i++) {
            QVector<RenderPoint> qlp;
            // Loop refinement
            for (QPoint loc : loops[i]) {
                // Binary search for the class point nearest the outside
                QPoint ilow((loc.x() - loc.x() % 2) / 2,
                            (loc.y() - loc.y() % 2) / 2);
                QPoint ihigh((loc.x() + loc.x() % 2) / 2,
                             (loc.y() + loc.y() % 2) / 2);
                bool low_inclass = true;
                if (ilow.x() < 0 || ilow.x() >= W || ilow.y() < 0 ||
                    ilow.y() >= H) {
                    low_inclass = false;
                } else if (grid_points[ilow.x() * H + ilow.y()].region_class !=
                           cls) {
                    low_inclass = false;
                }
                QPoint ip_in = low_inclass ? ilow : ihigh;
                QPoint ip_out = low_inclass ? ihigh : ilow;
                RenderPoint in_point = getPoint(ip_in);
                RenderPoint out_point = getPoint(ip_out);

                RenderPoint in_limit, out_limit;
                bracketEdge(in_point, out_point, &in_limit, &out_limit);
                in_limit.ideal_color =
                    calculateBoundaryColor(in_limit, out_limit);

                if (qlp.size() &&
                    (in_limit.coords - qlp.last().coords).manhattanLength() <
                        1e-15) {
                    // In case bisection search, both times, reaches the
                    // grid point.
                    continue;
                }
                qlp.push_back(in_limit);
            }

            /* Loop corner insertion -- when there's a sharp double angle
             * change, there most likely is a class point at its extension.
             */
            for (int k = 0; k < qlp.size(); k++) {
                int k0 = k, k1 = (k + 1) % qlp.size(),
                    k2 = (k + 2) % qlp.size(), k3 = (k + 3) % qlp.size();
                QPointF dir_pre = qlp[k1].coords - qlp[k0].coords;
                QPointF dir_post = qlp[k2].coords - qlp[k3].coords;
                double len_pre = std::sqrt(dir_pre.x() * dir_pre.x() +
                                           dir_pre.y() * dir_pre.y());
                double len_post = std::sqrt(dir_post.x() * dir_post.x() +
                                            dir_post.y() * dir_post.y());
                if (len_pre <= 0. || len_post <= 0.) {
                    qWarning(
                        "Pt. overlap %f %f | %f %f | %f %f | %f %f | %f %f",
                        len_pre, len_post, qlp[k0].coords.x(),
                        qlp[k0].coords.y(), qlp[k1].coords.x(),
                        qlp[k1].coords.y(), qlp[k2].coords.x(),
                        qlp[k2].coords.y(), qlp[k3].coords.x(),
                        qlp[k3].coords.y());
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

                QPointF src_pre = qlp[k1].coords, src_post = qlp[k2].coords;
                double t = (dir_post.x() * (src_post.y() - src_pre.y()) +
                            dir_post.y() * (src_pre.x() - src_post.x())) /
                           det;
                double s = (dir_pre.x() * (src_post.y() - src_pre.y()) +
                            dir_pre.y() * (src_pre.x() - src_post.x())) /
                           det;
                if (s < 0 || t < 0) {
                    // If not exactly on the angle, can find intersection
                    // point behind k0
                    continue;
                }
                QPointF qavg = 0.5 * (src_pre + src_post);
                QPointF qcor_pre = src_pre + t * dir_pre;
                QPointF qcor_post = src_post + s * dir_post;
                QPointF qcor = 0.5 * (qcor_pre + qcor_post);
                QPointF qjump = 2 * qcor - qavg;

                RenderPoint pavg = queryPoint(qavg);
                RenderPoint pjump = queryPoint(qjump);
                bool avg_inclass = typematch(pavg, class_representative);
                bool jump_inclass = typematch(pjump, class_representative);

                if (!avg_inclass || jump_inclass) {
                    continue;
                }

                // TODO: this doesn't get as close to the corner as we'd
                // like it might be off slightly
                RenderPoint plim_in, plim_out;
                bracketEdge(pavg, pjump, &plim_in, &plim_out);
                plim_in.ideal_color = calculateBoundaryColor(plim_in, plim_out);

                if ((plim_in.coords - src_pre).manhattanLength() < 1e-15 ||
                    (plim_in.coords - src_post).manhattanLength() < 1e-15) {
                    qDebug("too close to existing");
                    continue;
                }

                // Add point 2 ahead, and skip 3.
                // Note: this has a bug sometimes, possible re. wrapping
                qlp.insert(k2, plim_in);
                k += 2;
            }

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
    if (a0.normal.mag() < 1e-10 || a1.normal.mag() < 1e-10) {
        // empty normals
        return false;
    }
    if (std::abs(a0.normal.mag() - 1.0) > 1e-10 ||
        std::abs(a1.normal.mag() - 1.0) > 1e-10) {
        qFatal("Non-normal nontrivial normals %f %f", a0.normal.mag(),
               a1.normal.mag());
    }

    if (a0.normal.dot(a1.normal) < cos_alpha) {
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
    G4ThreeVector N = 0.5 * (a0.normal + a1.normal);
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
    if (!q0.intersections || !q1.intersections)
        return -1;

    if (transparent_volumes) {
        for (int i = 0; i <= std::min(q0.nhits, q1.nhits); i++) {
            const Intersection &a0 = q0.intersections[i];
            const Intersection &a1 = q1.intersections[i];
            if (crease_check(a0, a1, q0.coords, q1.coords, view_data, cos_alpha,
                             min_jump, is_jump)) {
                return i;
            }
        }
    } else {
        Intersection *first0 = first_nontrivial_intersection(q0);
        Intersection *first1 = first_nontrivial_intersection(q1);
        if (!first0 || !first1)
            return -1;

        if (crease_check(*first0, *first1, q0.coords, q1.coords, view_data,
                         cos_alpha, min_jump, is_jump)) {
            return 0;
        }
    }

    return -1;
}

void VectorTracer::computeCreases() {
    qDebug("Computing creases and boundary color details");

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
        while (pt_set.size()) {
            subreg_no += 1;
            pt_list = pt_set.toList().toVector();
            for (QPoint p : altered_list) {
                flood_fill_data[p.x() * H + p.y()] = 0;
            }
            altered_list.clear();

            QPoint p0 = pt_list[randint(pt_list.size())];
            grid_points[p0.x() * H + p0.y()].subregion_class = subreg_no;
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
                            RenderPoint rp = queryPoint(pt);
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
            subreg.linear_nsteps = 0;
            subreg.linear_colors.clear();

            subreg.boundaries.clear();
            for (const QPolygon &poly : loops) {
                QVector<RenderPoint> bound;

                // Only crease edges visible...
                for (QPoint line_marker : poly) {
                    QPoint p0 = line_marker / 2;
                    QPoint p1 = line_marker - line_marker / 2;

                    RenderPoint q0 = getPoint(p0);
                    RenderPoint q1 = getPoint(p1);

                    RenderPoint adj_point;
                    bool is_crease = false;
                    if (!crease_edge_map.contains(line_marker)) {
                        // Edge to out of region
                        // Pull typematch subdiv
                        RenderPoint lim_in, lim_out;
                        bracketEdge(
                            q0.region_class == region.class_no ? q0 : q1,
                            q0.region_class == region.class_no ? q1 : q0,
                            &lim_in, &lim_out);
                        lim_in.ideal_color =
                            calculateBoundaryColor(lim_in, lim_out);
                        adj_point = lim_in;
                    } else if (!crease_edge_map[line_marker]) {
                        // Interior edge, not conflict
                        // Select midpoint of q0 and q1
                        adj_point = queryPoint(0.5 * (q0.coords + q1.coords));
                        adj_point.ideal_color =
                            calculateInteriorColor(adj_point);
                    } else {
                        // Interior edge, crease.
                        // Pick error type (normal or displacement)
                        // and bisect on that parameter
                        // Color is inside point color
                        RenderPoint lim_in, lim_out;
                        bracketCrease(
                            q0.subregion_class == subreg.subclass_no ? q0 : q1,
                            q0.subregion_class == subreg.subclass_no ? q1 : q0,
                            &lim_in, &lim_out);
                        lim_in.ideal_color = calculateInteriorColor(lim_in);
                        adj_point = lim_in;
                        is_crease = true;
                    }

                    FColor aic = adj_point.ideal_color;
                    adj_point.ideal_color =
                        FColor(aic.redF(), aic.greenF(), aic.blueF(),
                               is_crease ? 1.0 : 0.0);
                    // Oh: rpoint, alpha=0 signifies invisible.
                    bound.push_back(adj_point);
                }

                // TODO: an edge expansion pass, identical to the one from
                // the parent boundary

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
                                      const QPointF &offset, double T) {
    QPolygonF a;
    for (const RenderPoint &p : loop) {
        a.push_back((p.coords + offset) * T);
    }
    return a;
}
static QVector<QPolygonF> boundary_loops_for_region(const Region &region,
                                                    const QPointF &offset,
                                                    double T) {
    QVector<QPolygonF> n;
    n.push_back(shift_and_scale_loop(region.exterior, offset, T));
    for (const QVector<RenderPoint> &loop : region.interior) {
        n.push_back(shift_and_scale_loop(loop, offset, T));
    }
    return n;
}
static QVector<QPolygonF> boundary_loops_for_subregion(const Subregion &subreg,
                                                       const QPointF &offset,
                                                       double T) {
    QVector<QPolygonF> n;
    for (const QVector<RenderPoint> &loop : subreg.boundaries) {
        n.push_back(shift_and_scale_loop(loop, offset, T));
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

static QPolygonF simplify_poly(const QPolygonF &orig, double max_error) {
    // First, locate sharpest angle; will be our starting point
    const int n = orig.size();
    int sharpest = 0;
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

    return poly;
}

static QString svg_path_from_polygons(const QVector<QPolygonF> &loops,
                                      QRgb color, const QString &id) {
    const int fprec = 10;
    QStringList path_string;
    for (const QPolygonF &poly : loops) {
        QPolygonF simpath = simplify_poly(poly, 1e-6);

        QPointF s = simpath[0];
        path_string.append(QString("M%1,%2")
                               .arg(s.x(), 0, 'g', fprec)
                               .arg(s.y(), 0, 'g', fprec));
        for (int i = 1; i < simpath.size(); i++) {
            QPointF q = simpath[i];
            path_string.append(QString("L%1,%2")
                                   .arg(q.x(), 0, 'g', fprec)
                                   .arg(q.y(), 0, 'g', fprec));
        }
        // Note: A rx ry x-axis-rotation large-arc-flag sweep-flag x y
        // gives elliptical arc, very suitable for path compression
        path_string.append("Z");
    }
    return QString("<path id=\"%3\" stroke=\"%1\" fill-rule=\"evenodd\" "
                   "d=\"%2\"/>\n")
        .arg(color_hex_name_rgb(color))
        .arg(path_string.join(" "))
        .arg(id);
}

static QRectF bounding_rect_for_boundary(const QVector<RenderPoint> &loop,
                                         QPointF offset, double T) {
    qreal xmax = -std::numeric_limits<qreal>::infinity();
    qreal xmin = +std::numeric_limits<qreal>::infinity();
    qreal ymax = -std::numeric_limits<qreal>::infinity();
    qreal ymin = +std::numeric_limits<qreal>::infinity();
    for (const RenderPoint &r : loop) {
        QPointF h = (r.coords + offset) * T;
        xmax = std::max(xmax, h.x());
        xmin = std::min(xmin, h.x());
        ymax = std::max(ymax, h.y());
        ymin = std::min(ymin, h.y());
    }
    return QRectF(xmin, ymin, xmax - xmin, ymax - ymin);
}

static double compute_histogram_angle(const Region &region,
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

                float dx = (mag_dx - mag_base) / dxs[k];
                float dy = (mag_dy - mag_base) / dys[k];
                if (dx == 0. && dy == 0.)
                    continue;
                // S1 -> RP1=S1. Via angle is clear, but could be optimized
                float angle = std::atan2(dy, dx);
                angle = angle * 2;
                sum_sin += std::sin(angle);
                sum_cos += std::cos(angle);
            }
        }
    }
    if (sum_sin == 0. && sum_cos == 0.) {
        // Uniformly distributed angles, or just none at all
        return 0.;
    }
    float angle = std::atan2(sum_sin, sum_cos);
    // The half-pi shift gives an angle perpendicular to dominant
    // gradient direction
    return std::fmod(-angle / 2 + M_PI / 2, M_PI);
}

static void fill_linear_histogram_for_angle(
    Subregion &region, const QVector<QPair<FColor, QPointF>> &interior_colors,
    const QSize &gridsize, double angle) {
    int S = std::max(gridsize.width(), gridsize.height()) - 1;
    double grid_spacing = 1.0 / S;

    // Determine histogram parameters
    double tmin = std::numeric_limits<float>::infinity();
    double tmax = -std::numeric_limits<float>::infinity();
    double proj_x = std::sin(angle);
    double proj_y = std::cos(angle);
    for (const QPair<FColor, QPointF> &p : interior_colors) {
        double t = proj_x * p.second.x() + proj_y * p.second.y();
        tmin = std::min(tmin, t);
        tmax = std::max(tmax, t);
    }
    double min_spacing = 1.0 * grid_spacing;
    int nsteps = std::max(2, 1 + (int)std::floor((tmax - tmin) / min_spacing));

    region.linear_angle = angle;
    region.linear_start = tmin;
    region.linear_stop = tmax;
    region.linear_nsteps = nsteps;
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
        double t = proj_x * p.second.x() + proj_y * p.second.y();
        double cell = (region.linear_nsteps - 1.) * (t - region.linear_start) /
                      (region.linear_stop - region.linear_start);
        int lidx = (int)std::floor(cell);
        lidx = std::max(0, std::min(region.linear_nsteps - 2, lidx));
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
    long sqe = 0.;
    for (const QPair<FColor, QPointF> &p : interior_colors) {
        QRgb color = p.first.rgba();

        QRgb pred;
        if (region.gradient_type == GradientType::gSolid) {
            pred = region.solid_color;
        } else if (region.gradient_type == GradientType::gLinear) {
            double proj_x = std::sin(region.linear_angle);
            double proj_y = std::cos(region.linear_angle);
            double t = proj_x * p.second.x() + proj_y * p.second.y();

            double cell = (region.linear_nsteps - 1.) *
                          (t - region.linear_start) /
                          (region.linear_stop - region.linear_start);
            int lidx = (int)std::floor(cell);
            lidx = std::max(0, std::min(region.linear_nsteps - 2, lidx));
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
                       subreg.linear_nsteps);
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
        float bnet_r = 0, bnet_g = 0, bnet_b = 0, bnet_a = 0;
        for (const RenderPoint &rp : region.exterior) {
            bnet_r += rp.ideal_color.redF();
            bnet_g += rp.ideal_color.greenF();
            bnet_b += rp.ideal_color.blueF();
            bnet_a += rp.ideal_color.alphaF();
        }
        int n = region.exterior.size();
        region.meanExteriorColor =
            qRgba(bnet_r / n, bnet_g / n, bnet_b / n, bnet_a / n);
        for (const QVector<RenderPoint> &loop : region.interior) {
            int inet_r = 0, inet_g = 0, inet_b = 0, inet_a = 0;
            int m = loop.size();
            for (const RenderPoint &rp : loop) {
                inet_r += rp.ideal_color.redF();
                inet_g += rp.ideal_color.greenF();
                inet_b += rp.ideal_color.blueF();
                inet_a += rp.ideal_color.alphaF();
            }
            region.meanInteriorColors.push_back(
                qRgba(inet_r / m, inet_g / m, inet_b / m, inet_a / m));
        }
    }

    //
    // Since QtSvg is incapable of handling clip paths like every
    // other SVG handling code, we roll our own SVG generator.
    //
    // Note: As QPainter supports this, fixing it in QtSvg on
    // render/generator sides isn't very hard.
    //

    const int S = std::max(W, H);
    const double px_per_mm = 3.78;
    double T = 45 * px_per_mm;
    QRectF viewbox(0., 0., T * (W + 2) / (double)S, T * (H + 2) / (double)S);

    // Gradient -- a piecewise composite of linear and radial gradients.
    // It's easy to fit a single gradient on a local patch;
    // then add a dividing crease along the mismatch between regions,
    // where opposing candidates are closest. (Again, march-cube style?,
    // then refined ? or A/B decision style.

    qDebug("Rendering final image, %d regions", region_list.size());
    {
        QFile ofile(file_name);
        ofile.open(QFile::WriteOnly | QFile::Truncate);
        QTextStream s(&ofile);

        s << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n "
             "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";

        s << QString("<svg version=\"1.1\" width=\"%1\" height=\"%2\" >\n")
                 .arg(viewbox.width())
                 .arg(viewbox.height());

        s << "  <defs>\n";

        const QPointF offset(0.5 * W / S + 1. / S, 0.5 * H / S + 1. / S);
        // Clipping paths
        for (const Region &region : region_list) {
            const QVector<QPolygonF> &loops =
                boundary_loops_for_region(region, offset, T);

            s << QString("    <clipPath id=\"boundary%1\">\n")
                     .arg(region.class_no);
            s << "        "
              << svg_path_from_polygons(
                     loops, qRgb(0, 0, 0),
                     QString("region_bound%1").arg(region.class_no));
            s << QString("    </clipPath>\n");
            for (const Subregion &subreg : region.subregions) {
                s << QString("    <clipPath id=\"boundary%1_%2\">\n")
                         .arg(region.class_no)
                         .arg(subreg.subclass_no);
                const QVector<QPolygonF> &subloops =
                    boundary_loops_for_subregion(subreg, offset, T);
                s << "        "
                  << svg_path_from_polygons(subloops, qRgb(0, 0, 0),
                                            QString("subregion_bound%1_%2")
                                                .arg(region.class_no)
                                                .arg(subreg.subclass_no));
                s << QString("    </clipPath>\n");
            }
        }
        // Hatching fill overlays
        for (const Region &region : region_list) {
            // TODO: per-material settings? normal angle dependence
            if (region.is_clipped_patch) {
                double angle = 60.0;
                double hatch_width = T / 40.;

                s << QString("    <pattern id=\"hatching%1\" width=\"%2\" "
                             "height=\"10\" patternTransform=\"rotate(%3 0 "
                             "0)\" patternUnits=\"userSpaceOnUse\">\n")
                         .arg(region.class_no)
                         .arg(hatch_width)
                         .arg(angle);
                s << QString("       <line x1=\"0\" y1=\"0\" x2=\"0\" "
                             "y2=\"10\" style=\"stroke:%1; stroke-width:%2\" "
                             "/>\n")
                         .arg(color_hex_name_rgb(qRgb(0, 0, 0)))
                         .arg(hatch_width);
                s << QString("    </pattern>\n");
            }
        }

        // Gradients
        for (const Region &region : region_list) {
            for (const Subregion &subreg : region.subregions) {
                if (subreg.gradient_type == GradientType::gLinear) {
                    s << QString("   <linearGradient id=\"gradient%1_%2\" "
                                 "x1=\"%3\" "
                                 "y1=\"%4\" x2=\"%5\" y2=\"%6\" "
                                 "spreadMethod=\"%7\" "
                                 "gradientUnits=\"userSpaceOnUse\">\n")
                             .arg(region.class_no)
                             .arg(subreg.subclass_no)
                             .arg(viewbox.center().x() -
                                  T * std::sin(subreg.linear_angle))
                             .arg(viewbox.center().y() -
                                  T * std::cos(subreg.linear_angle))
                             .arg(viewbox.center().x() +
                                  T * std::sin(subreg.linear_angle))
                             .arg(viewbox.center().y() +
                                  T * std::cos(subreg.linear_angle))
                             .arg(1 ? "repeat" : "pad");
                    for (int i = 0; i < subreg.linear_colors.size(); i++) {
                        double stop_pos =
                            (subreg.linear_start +
                             i * (subreg.linear_stop - subreg.linear_start) /
                                 (subreg.linear_nsteps - 1.));
                        if ((i <= 0 || subreg.linear_colors[i] ==
                                           subreg.linear_colors[i - 1]) &&
                            (i >= subreg.linear_colors.size() - 1 ||
                             subreg.linear_colors[i] ==
                                 subreg.linear_colors[i + 1])) {
                            // Skip redundant points
                            continue;
                        }

                        // in range [-0.5, 0.5] subset [-1.0,1.0]
                        stop_pos = (stop_pos + 1.0) / 2.0;
                        // in range [0.0, 1.0]

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
        s << QString(
            "<g fill-opacity=\"1\" stroke=\"none\" id=\"interiors\">\n");
        for (const Region &region : region_list) {
            // Compute region limit
            for (const Subregion &subreg : region.subregions) {

                QRectF subreg_limit =
                    bounding_rect_for_boundary(subreg.boundaries[0], offset, T);

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

                s << QString("  <rect id=\"region%5_%6\" x=\"%1\" y=\"%2\" "
                             "width=\"%3\" "
                             "height=\"%4\" "
                             "clip-path=\"url(#boundary%5_%6)\" "
                             "fill=\"%7\" />\n")
                         .arg(subreg_limit.x())
                         .arg(subreg_limit.y())
                         .arg(subreg_limit.width())
                         .arg(subreg_limit.height())
                         .arg(region.class_no)
                         .arg(subreg.subclass_no)
                         .arg(fill_desc);
            }
            QRectF region_limit =
                bounding_rect_for_boundary(region.exterior, offset, T);
            // Overlay mix with hatching
            if (region.is_clipped_patch) {
                s << QString("  <rect id=\"region_hatch%5\" x=\"%1\" y=\"%2\" "
                             "width=\"%3\" "
                             "height=\"%4\" style=\"fill: url(#hatching%5) "
                             "#fff;\" opacity=\"0.2\" "
                             "clip-path=\"url(#boundary%6)\"/>\n")
                         .arg(region_limit.x())
                         .arg(region_limit.y())
                         .arg(region_limit.width())
                         .arg(region_limit.height())
                         .arg(region.class_no)
                         .arg(region.class_no);
            }
        }
        s << QString("</g>\n");

        // Boundaries
        s << QString("<g fill=\"none\" stroke=\"black\" stroke-width=\"%1\" "
                     "fill-rule=\"evenodd\" stroke-linecap=\"square\" "
                     "stroke-linejoin=\"miter\" id=\"boundaries\">\n")
                 .arg(T * 0.003);
        for (const Region &region : region_list) {
            const QVector<QPolygonF> &loops =
                boundary_loops_for_region(region, offset, T);
            for (int i = 0; i < loops.size(); i++) {
                QVector<QPolygonF> solo;
                solo.push_back(loops[i]);
                QRgb c = i ? region.meanInteriorColors[i - 1]
                           : region.meanExteriorColor;
                s << "    "
                  << svg_path_from_polygons(
                         solo, c,
                         QString("edge%1_%2").arg(region.class_no).arg(i));
            }
        }
        s << QString("</g>\n");

        s << "</svg>\n";
    }

    // We proxy via inkscape because QtSvg barely supports SVG.
    // This takes half a second to run
    QTime time_start = QTime::currentTime();
    QStringList command;
    command << "inkscape";
    command << "--export-png=/tmp/write.png";
    command << "--export-dpi=800";
    command << file_name;
    QProcess proc;
    proc.start(command.first(), command.mid(1));
    proc.waitForFinished();
    QImage grad_image = QImage("/tmp/write.png");
    emit produceImagePhase(grad_image, QString("Render completed"), nqueries,
                           true);
    QTime time_end = QTime::currentTime();
    qDebug("Conversion took %f seconds", time_start.msecsTo(time_end) * 1e-3);
}
