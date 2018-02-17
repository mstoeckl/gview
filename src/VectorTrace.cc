#include "VectorTrace.hh"

#include <QApplication>
#include <QBitmap>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QImage>
#include <QLabel>
#include <QMenuBar>
#include <QPixmap>
#include <QPushButton>
#include <QSet>
#include <QStatusBar>
#include <QVBoxLayout>
#include <QtSvg/QSvgGenerator>

#include <G4VSolid.hh>
#include <G4VisExtent.hh>

ImageWidget::ImageWidget() : QWidget(), image() {
    QSizePolicy exp(QSizePolicy::Expanding, QSizePolicy::Expanding);
    this->setSizePolicy(exp);
}
ImageWidget::~ImageWidget() {}
void ImageWidget::setImage(QImage im) {
    // duplicate image to avoid case where data
    // is changed underneath it
    image = im.copy();
    this->update();
}
void ImageWidget::paintEvent(QPaintEvent *) {
    if (this->height() <= 0 || this->width() <= 0 || image.isNull()) {
        return;
    }
    int s = std::min(this->width() / image.width(),
                     this->height() / image.height());

    QPainter q(this);
    q.fillRect(rect(), Qt::gray);
    QImage mvd;
    if (s > 0) {
        mvd = image.scaled(image.width() * s, image.height() * s,
                           Qt::IgnoreAspectRatio, Qt::FastTransformation);

    } else {
        mvd = image.scaled(this->width(), this->height(), Qt::KeepAspectRatio,
                           Qt::SmoothTransformation);
    }
    q.drawImage(
        this->rect().center() - QPoint(mvd.width() / 2, mvd.height() / 2), mvd);
}

RenderPoint::RenderPoint() {
    coords = QPointF(0., 0.);
    intersections = NULL;
    elements = NULL;
    nhits = 0;
    ideal_color = qRgb(0., 0., 0.);
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
    ideal_color = qRgb(0., 0., 0.);
    region_class = -1;
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
    intersections = new Intersection[nhits + 1];
    elements = new Element *[nhits];
    memcpy(elements, other.elements, nhits * sizeof(Element *));
    memcpy(intersections, other.intersections,
           (nhits + 1) * sizeof(Intersection));
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
    std::swap(coords, other.coords);
    std::swap(ideal_color, other.ideal_color);
}

VectorTracer::VectorTracer(GeoOption option) : QMainWindow() {
    this->setWindowTitle(option.name.c_str());

    /* Initialize functional stuff */

    view_data.elements = convertCreation(option.vol);
    view_data.scene_radius =
        view_data.elements.solid->GetExtent().GetExtentRadius();
    view_data.scale = view_data.scene_radius / 20;
    view_data.camera = G4ThreeVector(-2 * view_data.scene_radius, 0, 0);
    view_data.orientation = CLHEP::HepRotationX(45 * CLHEP::deg);
    view_data.level_of_detail = 0;
    view_data.split_by_material = true;
    view_data.clipping_planes = std::vector<Plane>();
    Plane p;
    p.normal = G4ThreeVector(1, 0, 0);
    p.offset = 0.;
    view_data.clipping_planes.push_back(p);

    grid_size = QSize(999, 1001);
    grid_points = NULL;
    grid_nclasses = 0;
    ray_mutables = NULL;
    ray_iteration = 0;

    /* Initialize UI */
    this->menuBar()->addAction("Exit", this, SLOT(closeProgram()));

    step_next = Steps::sGrid;

    button_full = new QPushButton("Render Full");
    connect(button_full, SIGNAL(pressed()), SLOT(renderFull()));
    button_step = new QPushButton("Render Step");
    connect(button_step, SIGNAL(pressed()), SLOT(renderStep()));
    label_step = new QLabel("Step: ---");

    image_grid = new ImageWidget();
    image_edge = new ImageWidget();
    image_crease = new ImageWidget();
    image_gradient = new ImageWidget();
    image_final = new QGraphicsView();

    QHBoxLayout *layout_columns = new QHBoxLayout();

    QVBoxLayout *layout_control = new QVBoxLayout();
    layout_control->addWidget(button_full);
    layout_control->addWidget(button_step);
    layout_control->addWidget(label_step);
    layout_columns->addLayout(layout_control, 0);

    QGridLayout *layout_steps = new QGridLayout();
    layout_steps->addWidget(image_grid, 0, 0);
    layout_steps->addWidget(image_edge, 0, 1);
    layout_steps->addWidget(image_crease, 1, 0);
    layout_steps->addWidget(image_gradient, 1, 1);
    layout_columns->addLayout(layout_steps, 2);

    layout_columns->addWidget(image_final, 1);
    QWidget *central_widget = new QWidget();
    central_widget->setLayout(layout_columns);
    this->setCentralWidget(central_widget);

    this->statusBar()->showMessage("Status:");

    this->show();
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
        label_step->setText("Step: Grid");
        break;
    case Steps::sEdges:
        computeEdges();
        step_next = Steps::sCreases;
        label_step->setText("Step: Edges");
        break;
    case Steps::sCreases:
        computeCreases();
        step_next = Steps::sGradients;
        label_step->setText("Step: Creases");
        break;
    case Steps::sGradients:
        computeGradients();
        label_step->setText("Step: Gradients");
        step_next = Steps::sDone;
        break;
    }
}
static bool typematch(const RenderPoint &a, const RenderPoint &b) {
    if (a.nhits != b.nhits) {
        return false;
    }
    for (int i = 0; i < a.nhits; i++) {
        if (a.elements[i] != b.elements[i])
            return false;
    }
    return true;
}

RenderPoint VectorTracer::queryPoint(QPointF spot) {
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

QPointF grid_coord_to_point(const QPoint &pt, const QSize &grid_size) {
    const int W = grid_size.width() - 1, H = grid_size.height() - 1;
    QPointF spot((pt.x() / (double)W) - 0.5, (pt.y() / (double)H) - 0.5);
    return spot;
}
QPointF point_to_grid_coord(const QPointF &pt, const QSize &grid_size) {
    const int W = grid_size.width() - 1, H = grid_size.height() - 1;
    QPointF spot((pt.x() + 0.5) * W, (pt.y() + 0.5) * H);
    return spot;
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
        QRgb *colors = new QRgb[grid_nclasses];
        for (int i = 0; i < grid_nclasses; i++) {
            colors[i] = qRgb(qrand() % 256, qrand() % 256, qrand() % 256);
        }
        QRgb *dat = new QRgb[W * H];
        for (int x = 0; x < W; x++) {
            for (int y = 0; y < H; y++) {
                dat[y * W + x] = colors[grid_points[x * H + y].region_class];
            }
        }

        QImage img((uchar *)dat, W, H, QImage::Format_ARGB32);
        image_grid->setImage(img);
        delete[] dat;
        delete[] colors;
    }
}
inline uint qHash(QPointF t) { return qHash(t.x()) + 176 * qHash(t.y()); }
inline uint qHash(QLineF t) { return qHash(t.p1()) + 133 * qHash(t.p2()); }
inline bool operator<(const QPoint &a, const QPoint &b) {
    if (a.x() == b.x())
        return a.y() < b.y();
    return a.x() < b.x();
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

QColor randColor() {
    return QColor::fromRgb(qRgb(qrand() % 255, qrand() % 255, qrand() % 255));
}

void draw_boundaries_to(QPaintDevice *target,
                        const QVector<Boundary> &edge_boundaries, double S,
                        const QSize &grid_size) {

    QPainter p(target);
    p.fillRect(0, 0, target->width(), target->height(), Qt::white);
    for (const Boundary &b : edge_boundaries) {
        QPainterPath path;
        path.addPolygon(
            scale_shift_rp_loop(b.exterior, S, QPointF(S, S), grid_size));
        for (const QVector<RenderPoint> &loop : b.interior) {
            path.addPolygon(
                scale_shift_rp_loop(loop, S, QPointF(S, S), grid_size));
        }
        path.setFillRule(Qt::OddEvenFill);

        p.setPen(
            QPen(randColor(), 1, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin));
        p.setBrush(randColor());
        p.drawPath(path);
    }
}

void VectorTracer::computeEdges() {
    edge_boundaries.clear();
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

        /* Stitch together line segments into loops, via... eqvlclass routine!
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
                    Lns alt = links[k ? l.p1() : l.p2()];
                    int aid = (alt.ids[0] == test.id) ? alt.ids[1] : alt.ids[0];

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
                } else {
                    loop.push_back(l.p1());
                }
            }
            if (loop.first() != loop.last()) {
                qFatal("Invariant failure");
            }
            loop.pop_back();

            loops.push_back(loop);
        }

        /* Loop with greatest extent is the containing loop */
        int mxsz = 0;
        int bc = 0;
        for (int i = 0; i < loops.size(); i++) {
            // TODO: handle internal holes when it comes to it.
            const QRect &bound = loops[i].boundingRect();
            if (bound.width() * bound.height() > mxsz) {
                mxsz = bound.width() * bound.height();
                bc = i;
            }
        }

        qDebug("Refining region boundaries %d: %d loops", cls, loops.size());

        /* Refine all loops to stay barely within class. This _cannot_ be shared
         * between classes, as boundaries aren't perfect cliffs and there may be
         * third types in between. */
        const int nsubdivisions = 10;
        const RenderPoint &class_representative = grid_points[lx * H + ly];
        Boundary boundary;
        for (int i = 0; i < loops.size(); i++) {
            QVector<RenderPoint> qlp;
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

                QPointF p_in = grid_coord_to_point(ip_in, grid_size);
                QPointF p_out =
                    grid_coord_to_point(low_inclass ? ihigh : ilow, grid_size);
                qDebug("from (%f,%f) (%f,%f)", p_in.x(), p_in.y(), p_out.x(),
                       p_out.y());

                // need a class element
                RenderPoint in_point;
                in_point = grid_points[ip_in.x() * H + ip_in.y()];
                if (in_point.region_class != cls) {
                    qFatal("Invariant failure -- in_point not in class");
                }
                for (int k = 0; k < nsubdivisions; k++) {
                    QPointF p_mid = 0.5 * (p_in + p_out);
                    RenderPoint test = queryPoint(p_mid);
                    if (typematch(test, class_representative)) {
                        p_in = p_mid;
                        test.region_class = cls;
                        in_point = test;
                    } else {
                        p_out = p_mid;
                    }
                }
                qlp.push_back(in_point);
            }
            if (i == bc) {
                boundary.exterior = qlp;
            } else {
                boundary.interior.push_back(qlp);
            }
        }
        edge_boundaries.push_back(boundary);
    }

    qDebug("%d %d", edge_boundaries.size(), grid_nclasses);

    // Draw them lines!
    int S = 8;
    QImage lineImage((W + 1) * S, (H + 1) * S, QImage::Format_ARGB32);
    draw_boundaries_to(&lineImage, edge_boundaries, S, grid_size);
    image_edge->setImage(lineImage);

    QSvgGenerator target;
    target.setFileName("out.svg");
    draw_boundaries_to(&target, edge_boundaries, 1, grid_size);

    //    image_final->fitInView(scene->itemsBoundingRect());
}
void VectorTracer::computeCreases() {}
void VectorTracer::computeGradients() {}
void VectorTracer::closeProgram() { QApplication::quit(); }
