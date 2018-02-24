#include "Viewer.hh"

#include "ColorConfig.hh"
#include "CustomWidgets.hh"
#include "Overview.hh"
#include "RenderWidget.hh"
#include "VectorTrace.hh"

#include <G4GDMLParser.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VSolid.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VisExtent.hh>

#include <QAction>
#include <QActionGroup>
#include <QCheckBox>
#include <QCollator>
#include <QDockWidget>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QGridLayout>
#include <QHeaderView>
#include <QItemSelection>
#include <QLabel>
#include <QListWidget>
#include <QMenuBar>
#include <QPair>
#include <QPushButton>
#include <QResizeEvent>
#include <QScrollArea>
#include <QSignalMapper>
#include <QTableView>
#include <QTableWidget>
#include <QTreeView>
#include <QVBoxLayout>

static int exp10(int exp) {
    int r = 1;
    int base = 10;
    while (exp) {
        if (exp & 1)
            r *= base;
        exp >>= 1;
        base *= base;
    }
    return r;
}

static void includeMaterials(const G4LogicalVolume *p,
                             QSet<const G4Material *> &mtls) {
    const G4Material *mtl = p->GetMaterial();
    mtls.insert(mtl);
    for (int i = 0; i < p->GetNoDaughters(); i++) {
        includeMaterials(p->GetDaughter(i)->GetLogicalVolume(), mtls);
    }
}

std::vector<const G4Material *>
constructMaterialList(const std::vector<GeoOption> &geo_options) {
    QSet<const G4Material *> materials;
    for (const GeoOption &g : geo_options) {
        includeMaterials(g.vol->GetLogicalVolume(), materials);
    }
    QList<QPair<QString, const G4Material *>> mlist;
    for (const G4Material *g : materials) {
        mlist.append(QPair<QString, const G4Material *>(
            QString(g->GetName().data()), g));
    }
    qStableSort(mlist);
    std::vector<const G4Material *> mtl_list;
    for (const QPair<QString, const G4Material *> &p : mlist) {
        mtl_list.push_back(p.second);
    }
    return mtl_list;
}

Viewer::Viewer(const std::vector<GeoOption> &options,
               const std::vector<TrackData> &trackopts)
    : QMainWindow() {
    srand(1000 * QTime::currentTime().second() + QTime::currentTime().msec());

    which_geo = 0;
    geo_options = options;
    track_options = trackopts;
    vd.elements = convertCreation(geo_options[which_geo].vol);
    vd.scene_radius = vd.elements.solid->GetExtent().GetExtentRadius();
    vd.scale = 2 * vd.scene_radius;
    vd.camera = G4ThreeVector(-4 * vd.scene_radius, 0, 0);
    vd.orientation = G4RotationMatrix();
    vd.level_of_detail = 0; // 0 is full; 1 is 1/9, 2 is 1/81; depends on timing
    vd.split_by_material = true;
    vd.clipping_planes = std::vector<Plane>();
    which_tracks = track_options.size() > 0 ? 1 : 0;

    TrackData td =
        which_tracks == 0 ? TrackData() : track_options[which_tracks - 1];
    for (size_t i = 0; i < track_options.size(); i++) {
        TrackRestriction rests;
        track_options[i].calcTimeBounds(rests.time.low, rests.time.high);
        rests.time.high += 0.01 * CLHEP::ns;
        track_options[i].calcEnergyBounds(rests.energy.low, rests.energy.high);
        rests.energy.high *= 1.2;
        rests.energy.low /= 1.2;
        rests.energy.low = std::max(1.0 * CLHEP::eV, rests.energy.low);
        rests.seqno.low = 1;
        rests.seqno.high = track_options[i].getNTracks();
        track_res_actual.push_back(rests);
        track_res_bounds.push_back(rests);
    }
    if (which_tracks == 0) {
        trackdata = TrackData();
    } else {
        TrackRestriction current = track_res_actual[which_tracks - 1];
        trackdata =
            TrackData(td, vd, current.time, current.energy, current.seqno);
    }
    rayiter = 0;

    QAction *clipAction = new QAction("Clipping", this);
    clipAction->setToolTip("Edit geometry clipping planes and track limits");
    connect(clipAction, SIGNAL(triggered()), this, SLOT(restClip()));

    QAction *treeAction = new QAction("Tree", this);
    treeAction->setToolTip("View object hierarchy");
    connect(treeAction, SIGNAL(triggered()), this, SLOT(restTree()));

    QAction *infoAction = new QAction("Info", this);
    infoAction->setToolTip("Misc information on selected volume");
    connect(infoAction, SIGNAL(triggered()), this, SLOT(restInfo()));

    QAction *rayAction = new QAction("Ray", this);
    rayAction->setToolTip("List currently-hovered-over visible objects");
    connect(rayAction, SIGNAL(triggered()), this, SLOT(restRay()));

    QAction *mtlAction = new QAction("Color", this);
    mtlAction->setToolTip("Control volume colors");
    connect(mtlAction, SIGNAL(triggered()), this, SLOT(restColor()));

    QAction *oriAction = new QAction("Orientation", this);
    oriAction->setToolTip("Select a camera direction");
    connect(oriAction, SIGNAL(triggered()), this, SLOT(restOrient()));

    QAction *screenAction = new QAction("Screenshot", this);
    screenAction->setToolTip("Take a screenshot of active scene");
    connect(screenAction, SIGNAL(triggered()), this, SLOT(screenshot()));

    QAction *vectorTAction = new QAction("Vector (Transparent)", this);
    vectorTAction->setToolTip(
        "Take a transparent vector screenshot of active scene");
    connect(vectorTAction, SIGNAL(triggered()), this,
            SLOT(vectorTScreenshot()));

    QAction *vectorOAction = new QAction("Vector (Opaque)", this);
    vectorOAction->setToolTip(
        "Take an opaque vector screenshot of active scene");
    connect(vectorOAction, SIGNAL(triggered()), this,
            SLOT(vectorOScreenshot()));

    QAction *screen4Action = new QAction("ScrSht4x", this);
    screen4Action->setToolTip(
        "Take a screenshot of active scene; at 4x sampling");
    QSignalMapper *amp4 = new QSignalMapper(screen4Action);
    connect(screen4Action, SIGNAL(triggered()), amp4, SLOT(map()));
    amp4->setMapping(screen4Action, 4);
    connect(amp4, SIGNAL(mapped(int)), this, SLOT(screenshot(int)));

    clicked = false;
    shift = false;

    // TODO: icons for all actions, esp. screenshot actions
    // & so on. Bundle?
    gpicker_menu = this->menuBar()->addMenu("Choose Geometry");
    tpicker_menu = this->menuBar()->addMenu("Choose Tracks");
    reloadChoiceMenus();
    // TODO: checkbox by visibility
    QMenu *sub = this->menuBar()->addMenu("View");
    sub->addAction(clipAction);
    sub->addAction(treeAction);
    sub->addAction(infoAction);
    sub->addAction(rayAction);
    sub->addAction(mtlAction);
    sub->addAction(oriAction);
    this->menuBar()->addSeparator();
    QMenu *imgmenu = this->menuBar()->addMenu("Screenshot");
    imgmenu->addAction(screenAction);
    imgmenu->addAction(screen4Action);
    imgmenu->addAction(vectorTAction);
    imgmenu->addAction(vectorOAction);

    rwidget = new RenderWidget(vd, trackdata);
    this->setCentralWidget(rwidget);
    rwidget->setFocusPolicy(Qt::WheelFocus);

    // Clipping plane control
    dock_clip = new QDockWidget("Clipping", this);
    QWidget *cont = new QWidget();
    QVBoxLayout *vb = new QVBoxLayout();
    Plane iplanes[3];
    iplanes[0].offset = 0;
    iplanes[1].offset = 0;
    iplanes[2].offset = 0;
    iplanes[0].normal = G4ThreeVector(1, 0, 0);
    iplanes[1].normal = G4ThreeVector(0, 1, 0);
    iplanes[2].normal = G4ThreeVector(0, 0, 1);
    for (int i = 0; i < 3; i++) {
        plane_edit[i] = new PlaneEdit(iplanes[i]);
        vb->addWidget(plane_edit[i], 0);
        connect(plane_edit[i], SIGNAL(updated()), this, SLOT(updatePlanes()));
    }
    times_range = new HistogrammicRangeSlider(false);
    times_range->setSuffix("ns");
    energy_range = new HistogrammicRangeSlider(true);
    energy_range->setSuffix("eV");
    count_lower = new QSpinBox();
    count_upper = new QSpinBox();
    if (which_tracks > 0) {
        const TrackRestriction &res = track_res_actual[which_tracks - 1];
        times_range->setRange(res.time.low / CLHEP::ns,
                              res.time.high / CLHEP::ns);
        times_range->setValue(res.time.low / CLHEP::ns,
                              res.time.high / CLHEP::ns);
        energy_range->setRange(res.energy.low / CLHEP::eV,
                               res.energy.high / CLHEP::eV);
        energy_range->setValue(res.energy.low / CLHEP::eV,
                               res.energy.high / CLHEP::eV);
        count_lower->setRange(1, int(res.seqno.high));
        count_upper->setRange(1, int(res.seqno.high));
        count_lower->setValue(0);
        count_upper->setValue(int(res.seqno.high));
        int step = exp10(
            std::max(0, int(std::floor(std::log10(res.seqno.high / 15.)))));
        count_lower->setSingleStep(step);
        count_upper->setSingleStep(step);
    } else {
        times_range->setDisabled(true);
        energy_range->setDisabled(true);
        count_lower->setDisabled(true);
        count_upper->setDisabled(true);
    }
    connect(count_lower, SIGNAL(valueChanged(int)), this, SLOT(updatePlanes()));
    connect(count_upper, SIGNAL(valueChanged(int)), this, SLOT(updatePlanes()));
    connect(times_range, SIGNAL(valueChanged()), this, SLOT(updatePlanes()));
    connect(energy_range, SIGNAL(valueChanged()), this, SLOT(updatePlanes()));
    vb->addWidget(times_range);
    vb->addWidget(energy_range);
    times_range->setMaximumWidth(plane_edit[0]->sizeHint().width());
    energy_range->setMaximumWidth(plane_edit[0]->sizeHint().width());
    QHBoxLayout *crb = new QHBoxLayout();
    crb->addWidget(count_lower);
    crb->addWidget(count_upper);
    vb->addLayout(crb, 0);
    linecount_label = new QLabel("Lines: 0");
    vb->addWidget(linecount_label);
    vb->addStretch(1);
    cont->setLayout(vb);
    QScrollArea *asf = new QScrollArea();
    asf->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    asf->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    asf->setWidget(cont);
    asf->setWidgetResizable(true);
    cont->setSizePolicy(
        QSizePolicy(QSizePolicy::Maximum, QSizePolicy::Maximum));
    asf->setSizePolicy(
        QSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred));
    dock_clip->setWidget(asf);

    // Tree view (with vis/novis, hue control
    dock_tree = new QDockWidget("Element tree", this);
    tree_view = new QTreeView();
    tree_model = new OverView(vd, this);
    tree_view->setSortingEnabled(false);
    tree_view->setModel(tree_model);
    tree_view->setRootIndex(QModelIndex());
    tree_view->setSelectionMode(QAbstractItemView::SingleSelection);
    tree_view->setSelectionBehavior(QAbstractItemView::SelectRows);
    tree_view->setIndentation(20);
    tree_view->collapseAll();
    tree_view->expandToDepth(1);
    tree_view->header()->resizeSections(QHeaderView::ResizeToContents);
    tree_view->header()->setStretchLastSection(true);
    tree_view->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    tree_view->setItemDelegateForColumn(2,
                                        new AlphaBoxDelegate(tree_model, this));
    connect(tree_view, SIGNAL(activated(const QModelIndex &)), tree_model,
            SLOT(respToActive(const QModelIndex &)));
    connect(tree_model, SIGNAL(colorChange()), rwidget, SLOT(rerender()));

    connect(
        tree_view->selectionModel(),
        SIGNAL(
            selectionChanged(const QItemSelection &, const QItemSelection &)),
        tree_model,
        SLOT(respToSelection(const QItemSelection &, const QItemSelection &)));
    connect(tree_model, SIGNAL(selectedElement(const Element *)), this,
            SLOT(indicateElement(const Element *)));
    dock_tree->setWidget(tree_view);

    // Info for a particular item (mass, volume, mtl specs; table...)
    dock_info = new QDockWidget("Element info", this);
    info_table = new QTableView();
    info_model = new InfoModel(this);
    info_table->setModel(info_model);
    info_table->horizontalHeader()->setStretchLastSection(true);
    info_table->horizontalHeader()->hide();
    dock_info->setWidget(info_table);

    dock_ray = new QDockWidget("Ray table", this);
    ray_table = new QListWidget();
    dock_ray->setWidget(ray_table);
    connect(ray_table, SIGNAL(itemSelectionChanged()), this, SLOT(rayLookup()));

    dock_color = new QDockWidget("Color", this);
    QWidget *mtlc = new QWidget();
    QVBoxLayout *mla = new QVBoxLayout();
    mtl_showlines = new QCheckBox("Render lines");
    mtl_showlines->setCheckState(Qt::Unchecked);
    connect(mtl_showlines, SIGNAL(stateChanged(int)), this,
            SLOT(updateShowLines()));
    std::vector<const G4Material *> mtl_list =
        constructMaterialList(geo_options);
    color_config = new ColorConfig(vd, mtl_list);
    color_config->reassignColors();
    connect(color_config, SIGNAL(colorChange()), this, SLOT(updateColors()));
    mla->addWidget(mtl_showlines);
    mla->addWidget(color_config);
    mtlc->setLayout(mla);
    dock_color->setWidget(mtlc);

    dock_orient = new QDockWidget("Orientation", this);
    QPushButton *orp1 = new QPushButton("X+");
    QPushButton *orp2 = new QPushButton("Y+");
    QPushButton *orp3 = new QPushButton("Z+");
    QPushButton *orm1 = new QPushButton("X-");
    QPushButton *orm2 = new QPushButton("Y-");
    QPushButton *orm3 = new QPushButton("Z-");
    QPushButton *orr1 = new QPushButton("+45");
    QPushButton *orr2 = new QPushButton("-45");
    QHBoxLayout *prow = new QHBoxLayout();
    QSignalMapper *sqm = new QSignalMapper(this);
    connect(sqm, SIGNAL(mapped(int)), this, SLOT(setViewRotation(int)));
    sqm->setMapping(orp1, 0);
    sqm->setMapping(orp2, 1);
    sqm->setMapping(orp3, 2);
    sqm->setMapping(orm1, 3);
    sqm->setMapping(orm2, 4);
    sqm->setMapping(orm3, 5);
    sqm->setMapping(orr1, 6);
    sqm->setMapping(orr2, 7);
    connect(orp1, SIGNAL(clicked()), sqm, SLOT(map()));
    connect(orp2, SIGNAL(clicked()), sqm, SLOT(map()));
    connect(orp3, SIGNAL(clicked()), sqm, SLOT(map()));
    connect(orm1, SIGNAL(clicked()), sqm, SLOT(map()));
    connect(orm2, SIGNAL(clicked()), sqm, SLOT(map()));
    connect(orm3, SIGNAL(clicked()), sqm, SLOT(map()));
    connect(orr1, SIGNAL(clicked()), sqm, SLOT(map()));
    connect(orr2, SIGNAL(clicked()), sqm, SLOT(map()));
    prow->addWidget(orp1);
    prow->addWidget(orp2);
    prow->addWidget(orp3);
    QHBoxLayout *mrow = new QHBoxLayout();
    mrow->addWidget(orm1);
    mrow->addWidget(orm2);
    mrow->addWidget(orm3);
    QHBoxLayout *rrow = new QHBoxLayout();
    rrow->addWidget(orr1);
    rrow->addWidget(orr2);
    QVBoxLayout *vbr = new QVBoxLayout();
    vbr->addLayout(prow, 0);
    vbr->addLayout(mrow, 0);
    vbr->addLayout(rrow, 0);
    vbr->addStretch(5);
    QWidget *ornt = new QWidget();
    ornt->setLayout(vbr);
    dock_orient->setWidget(ornt);

    QDockWidget::DockWidgetFeatures feat = QDockWidget::DockWidgetClosable |
                                           QDockWidget::DockWidgetMovable |
                                           QDockWidget::DockWidgetFloatable;
    addDockWidget(Qt::LeftDockWidgetArea, dock_clip);
    dock_clip->setAllowedAreas(Qt::AllDockWidgetAreas);
    dock_clip->setFeatures(feat);
    dock_clip->setVisible(false);

    addDockWidget(Qt::LeftDockWidgetArea, dock_tree);
    dock_tree->setAllowedAreas(Qt::AllDockWidgetAreas);
    dock_tree->setFeatures(feat);
    dock_tree->setVisible(false);

    addDockWidget(Qt::LeftDockWidgetArea, dock_info);
    dock_info->setAllowedAreas(Qt::AllDockWidgetAreas);
    dock_info->setFeatures(feat);
    dock_info->setVisible(false);

    addDockWidget(Qt::LeftDockWidgetArea, dock_ray);
    dock_ray->setAllowedAreas(Qt::AllDockWidgetAreas);
    dock_ray->setFeatures(feat);
    dock_ray->setVisible(false);

    addDockWidget(Qt::LeftDockWidgetArea, dock_color);
    dock_color->setAllowedAreas(Qt::AllDockWidgetAreas);
    dock_color->setFeatures(feat);
    dock_color->setVisible(false);

    addDockWidget(Qt::LeftDockWidgetArea, dock_orient);
    dock_orient->setAllowedAreas(Qt::AllDockWidgetAreas);
    dock_orient->setFeatures(feat);
    dock_orient->setVisible(false);

    setMouseTracking(true);

    this->setMinimumSize(QSize(400, 400));
    this->show();
}

Viewer::~Viewer() {}

void Viewer::restClip() { dock_clip->setVisible(true); }
void Viewer::restTree() { dock_tree->setVisible(true); }
void Viewer::restInfo() { dock_info->setVisible(true); }
void Viewer::restRay() { dock_ray->setVisible(true); }
void Viewer::restColor() { dock_color->setVisible(true); }
void Viewer::restOrient() { dock_orient->setVisible(true); }

void Viewer::keyPressEvent(QKeyEvent *event) {
    G4double mvd = 0.1;

    G4ThreeVector trans;
    G4RotationMatrix rot;
    switch (event->key()) {
    case Qt::Key_Up:
        trans = vd.scale * vd.orientation.rowY() * mvd;
        break;
    case Qt::Key_Down:
        trans = -vd.scale * vd.orientation.rowY() * mvd;
        break;
    case Qt::Key_Left:
        trans = vd.scale * vd.orientation.rowZ() * mvd;
        break;
    case Qt::Key_Right:
        trans = -vd.scale * vd.orientation.rowZ() * mvd;
        break;
    // Note: corotate the camera vector
    case Qt::Key_W:
        rot = CLHEP::HepRotationZ(-atan2(vd.scale, vd.scene_radius) / 12);
        break;
    case Qt::Key_S:
        rot = CLHEP::HepRotationZ(atan2(vd.scale, vd.scene_radius) / 12);
        break;
    case Qt::Key_A:
        rot = CLHEP::HepRotationY(-atan2(vd.scale, vd.scene_radius) / 12);
        break;
    case Qt::Key_D:
        rot = CLHEP::HepRotationY(atan2(vd.scale, vd.scene_radius) / 12);
        break;
    case Qt::Key_PageUp:
        vd.scale *= std::pow(2, 1 / 12.);
        break;
    case Qt::Key_PageDown:
        vd.scale /= std::pow(2, 1 / 12.);
        break;
    default:
        return;
    }
    vd.camera -= trans;
    vd.camera =
        (vd.orientation.inverse() * rot.inverse() * vd.orientation) * vd.camera;
    vd.orientation = rot * vd.orientation;

#if 0
    G4cout << "Scale: " << vd.scale << G4endl;
    G4cout << "Camera: " << vd.camera << G4endl;
    G4cout << "Ori: " << vd.orientation << G4endl;
#endif
    rwidget->rerender();
}

void Viewer::mousePressEvent(QMouseEvent *event) {
    if (event->modifiers() & Qt::ShiftModifier) {
        shift = true;
    } else {
        shift = false;
    }
    clicked = true;
    clickpt = event->pos();
    lastpt = clickpt;
}

void Viewer::mouseReleaseEvent(QMouseEvent *) {
    // or mouse exit?
    clicked = false;
}

void Viewer::mouseMoveEvent(QMouseEvent *event) {
    // Ray tracking
    int h = rwidget->geometry().height();
    int w = rwidget->geometry().width();
    int mind = std::min(w, h);
    QPoint coord = rwidget->mapFromGlobal(event->globalPos());
    QPointF pt((coord.x() - w / 2.) / (2. * mind),
               (coord.y() - h / 2.) / (2. * mind));
    const int M = 50;
    Intersection ints[M + 1];
    const Element *elems[M];
    int td, nelem;
    countTree(vd.elements, td, nelem);
    Q_UNUSED(td);
    ElemMutables *mutables = new ElemMutables[nelem]();
    int m = traceRay(pt, vd, elems, ints, M, 1, mutables);
    delete[] mutables;
    rayiter++;
    m = compressTraces(elems, ints, m);
    ray_table->clear();
    ray_list.clear();
    for (int j = 0; j < m; j++) {
        QString name(elems[j]->name.data());
        ray_table->addItem(name);
        ray_list.push_back(elems[j]);
    }

    if (!clicked)
        return;
    int dmm = std::min(this->width(), this->height());
    if (shift) {
        // Translate with mouse
        QPoint delta = event->pos() - lastpt;
        QPointF dp = vd.scale * QPointF(delta) / dmm * 0.5;
        vd.camera -=
            vd.orientation.rowZ() * dp.x() + vd.orientation.rowY() * dp.y();
        lastpt = event->pos();
    } else {
        // All rotations in progress are relative to start point
        G4double step = atan2(vd.scale, vd.scene_radius) * 4.0;
        QPointF nalph = QPointF(event->pos() - clickpt) / dmm * step;
        QPointF palph = QPointF(lastpt - clickpt) / dmm * step;

        G4RotationMatrix next =
            CLHEP::HepRotationY(nalph.x()) * CLHEP::HepRotationZ(-nalph.y());
        G4RotationMatrix prev =
            CLHEP::HepRotationY(palph.x()) * CLHEP::HepRotationZ(-palph.y());

        G4RotationMatrix rot = prev.inverse() * next;
        // rotate
        vd.camera =
            (vd.orientation.inverse() * rot.inverse() * vd.orientation) *
            vd.camera;
        vd.orientation = rot * vd.orientation;
        lastpt = event->pos();
    }
    rwidget->rerender();
}

void Viewer::wheelEvent(QWheelEvent *event) {
    // Argh, still not good....
    vd.scale *= std::pow(2, event->delta() / 4800.);
    if (vd.scale > 8 * vd.scene_radius) {
        vd.scale = 8 * vd.scene_radius;
    }
    if (vd.scale < vd.scene_radius / 1048576.) {
        vd.scale = vd.scene_radius / 1048576.;
    }
    rwidget->rerender();
}

void Viewer::updatePlanes() {
    vd.clipping_planes = std::vector<Plane>();
    for (int i = 0; i < 3; i++) {
        Plane p = plane_edit[i]->getPlane();
        if (p.normal != G4ThreeVector()) {
            vd.clipping_planes.push_back(p);
        }
    }
    if (which_tracks > 0) {
        double tlow, thigh;
        times_range->value(tlow, thigh);
        tlow *= CLHEP::ns;
        thigh *= CLHEP::ns;
        double elow, ehigh;
        energy_range->value(elow, ehigh);
        elow *= CLHEP::eV;
        ehigh *= CLHEP::eV;
        size_t nlow = size_t(count_lower->value());
        size_t nhigh = size_t(count_upper->value());

        TrackRestriction &res = track_res_actual[which_tracks - 1];
        res.energy = {std::min(elow, ehigh), std::max(elow, ehigh)};
        res.time = {std::min(tlow, thigh), std::max(tlow, thigh)};
        res.seqno = {std::min(nlow, nhigh), std::max(nlow, nhigh)};
        trackdata = TrackData(track_options[which_tracks - 1], vd, res.time,
                              res.energy, res.seqno);
        if (0) {
            // Temp disabled on grounds of lag
            QVector<QPointF> ep, tp;
            track_options[which_tracks - 1].constructRangeHistograms(
                tp, ep, res.time, res.energy);
            energy_range->setHistogram(ep);
            times_range->setHistogram(tp);
        }
    } else {
        // Q: how to pull QSharedData on the Elements as well
        trackdata = TrackData();
    }
    linecount_label->setText(QString("Lines: %1").arg(trackdata.getNTracks()));
    rwidget->rerender();
}

void Viewer::updateColors() {
    // ColorConfig handles ViewData updates
    color_config->reassignColors();
    info_model->setElement(info_model->curE(), vd);
    rwidget->rerender();
}

void Viewer::changeGeometry(QAction *act) {
    QString s = act->toolTip();
    G4String s2(s.toUtf8().constData());
    for (size_t i = 0; i < geo_options.size(); i++) {
        if (s2 == geo_options[i].name) {
            if (which_geo != i) {
                // Change geometry
                which_geo = i;
                vd.elements = convertCreation(geo_options[which_geo].vol);
                vd.scene_radius =
                    vd.elements.solid->GetExtent().GetExtentRadius();
                if (4 * vd.scene_radius > vd.camera.mag()) {
                    vd.camera *= 4 * vd.scene_radius / vd.camera.mag();
                }
                std::vector<const G4Material *> mtl_list =
                    constructMaterialList(geo_options);
                color_config->mergeMaterials(mtl_list);
                color_config->reassignColors();
                tree_model->recalculate();
                tree_view->collapseAll();
                tree_view->expandToDepth(1);
                indicateElement(NULL);
                rwidget->rerender();
                return;
            }
        }
    }
}

void Viewer::changeTracks(QAction *act) {
    const QList<QAction *> &o = act->actionGroup()->actions();
    size_t c = 0;
    for (QAction *q : o) {
        if (q == act) {
            break;
        }
        c++;
    }
    which_tracks = c;
    bool active = which_tracks > 0;
    if (active) {
        times_range->blockSignals(true);
        energy_range->blockSignals(true);
        count_lower->blockSignals(true);
        count_upper->blockSignals(true);
        const TrackRestriction &res = track_res_bounds[which_tracks - 1];
        const TrackRestriction &cur = track_res_actual[which_tracks - 1];
        times_range->setRange(res.time.low / CLHEP::ns,
                              res.time.high / CLHEP::ns);
        times_range->setValue(cur.time.low / CLHEP::ns,
                              cur.time.high / CLHEP::ns);
        energy_range->setRange(res.energy.low / CLHEP::eV,
                               res.energy.high / CLHEP::eV);
        energy_range->setValue(cur.energy.low / CLHEP::eV,
                               cur.energy.high / CLHEP::eV);
        count_lower->setRange(1, int(res.seqno.high));
        count_upper->setRange(1, int(res.seqno.high));
        count_lower->setValue(int(cur.seqno.low));
        count_upper->setValue(int(cur.seqno.high));
        int step = exp10(
            std::max(0, int(std::floor(std::log10(res.seqno.high / 15.)))));
        count_lower->setSingleStep(step);
        count_upper->setSingleStep(step);
        times_range->blockSignals(false);
        energy_range->blockSignals(false);
        count_lower->blockSignals(false);
        count_upper->blockSignals(false);
    }
    times_range->setEnabled(active);
    energy_range->setEnabled(active);
    count_lower->setEnabled(active);
    count_upper->setEnabled(active);
    updatePlanes();
}

void Viewer::indicateElement(const Element *e) {
    info_model->setElement(e, vd);
}

void Viewer::rayLookup() {
    QList<QListWidgetItem *> li = ray_table->selectedItems();
    if (li.size() != 1) {
        return;
    }
    int r = ray_table->row(li[0]);
    const Element *e = ray_list[r];

    QModelIndex index = tree_model->indexFromElement(e);
    tree_view->scrollTo(index, QAbstractItemView::PositionAtCenter);
    tree_view->selectionModel()->clearSelection();
    tree_view->selectionModel()->select(index, QItemSelectionModel::Select);
    // ^ Triggers indicateElement
}

void Viewer::updateShowLines() {
    rwidget->setFullDetail(mtl_showlines->isChecked());
}

void Viewer::screenshot(int sx) {
    int w = rwidget->width() * sx, h = rwidget->height() * sx;
    RenderSaveObject *rso = new RenderSaveObject(vd, trackdata, w, h);
    rso->start();
}
void Viewer::vectorTScreenshot() {
    VectorTracer *vt =
        new VectorTracer(vd, trackdata, "vector_transparent.svg", true);
    vt->reset(true, QSize(1000, 1000));
    vt->renderFull();
    delete vt;
}
void Viewer::vectorOScreenshot() {
    VectorTracer *vt =
        new VectorTracer(vd, trackdata, "vector_opaque.svg", false);
    vt->reset(false, QSize(1000, 1000));
    vt->renderFull();
    delete vt;
}
void Viewer::reloadChoiceMenus() {
    gpicker_menu->clear();
    tpicker_menu->clear();

    QActionGroup *opts = new QActionGroup(this);
    for (size_t i = 0; i < geo_options.size(); i++) {
        QString s = QString(geo_options[i].name.c_str());
        QAction *ch = new QAction(s, this);
        ch->setToolTip(s);
        ch->setCheckable(true);
        opts->addAction(ch);
        if (i == which_geo) {
            ch->setChecked(true);
        }
        gpicker_menu->addAction(ch);
    }
    QAction *gadd = gpicker_menu->addAction("Open Geometry");
    connect(opts, SIGNAL(triggered(QAction *)), this,
            SLOT(changeGeometry(QAction *)));
    QActionGroup *topts = new QActionGroup(this);
    QAction *base = new QAction("None", this);
    base->setToolTip(base->text());
    base->setCheckable(true);
    if (track_options.size() == 0) {
        base->setChecked(true);
    }
    topts->addAction(base);
    tpicker_menu->addAction(base);
    if (which_tracks == 0) {
        base->setChecked(true);
    }
    for (size_t i = 0; i < track_options.size(); i++) {
        QString s = QString("%1").arg(i + 1);
        QAction *ch = new QAction(s, this);
        ch->setToolTip(s);
        ch->setCheckable(true);
        topts->addAction(ch);
        if (i == which_tracks - 1) {
            ch->setChecked(true);
        }
        tpicker_menu->addAction(ch);
    }
    QAction *tadd = tpicker_menu->addAction("Open Tracks");
    connect(topts, SIGNAL(triggered(QAction *)), this,
            SLOT(changeTracks(QAction *)));
    connect(gadd, SIGNAL(triggered()), this, SLOT(openGeometry()));
    connect(tadd, SIGNAL(triggered()), this, SLOT(openTracks()));
}

void Viewer::openGeometry() {
    QString selected = "GDML (*.gdml *.gdml.gz)";
    QString fn = QFileDialog::getOpenFileName(
        NULL, "Open Geometry", QDir::currentPath(),
        "All files (*.*);;GDML (*.gdml *.gdml.gz)", &selected);
    if (fn.size() >= 8 && fn.right(8) == ".gdml.gz") {
        system("rm -f /tmp/copy.gdml.gz");
        QString ar = "cp " + fn + " /tmp/copy.gdml.gz";
        system(ar.toUtf8().constData());
        QString gu = "gzip -df /tmp/copy.gdml.gz";
        system(gu.toUtf8().constData());
        fn = "/tmp/copy.gdml";
    }
    if (fn.size() >= 5 && fn.right(5) == ".gdml") {

        G4GDMLParser p;
        p.SetAddPointerToName(true);
        G4cout << "Started reading (may take a while)..." << G4endl;
        GeoOption g;
        g.name = G4String(fn.toUtf8().constData());
        p.Read(g.name, false);
        G4cout << "Done reading..." << G4endl;
        // Need to modify volume name to prevent collisions in lookup
        g.vol = p.GetWorldVolume();
        char buf[30];
        sprintf(buf, "-%d", int(geo_options.size()));
        G4String name = g.vol->GetName() + buf;
        g.vol->SetName(name);
        g.vol->GetLogicalVolume()->SetName(name);
        geo_options.push_back(g);
        G4cout << "Done converting..." << G4endl;
        p.Clear();
        reloadChoiceMenus();
    } else {
        qDebug("File not a GDML file!");
    }
}

void Viewer::openTracks() {
    QString selected = "Track data (*.dat *.dat.gz *.track *.track.gz)";
    QString fn =
        QFileDialog::getOpenFileName(NULL, "Open Tracks", QDir::currentPath(),
                                     "All files (*.*);;" + selected, &selected);
    if (fn.size() >= 3 && fn.right(3) == ".gz") {
        system("rm -f /tmp/copy.dat.gz");
        QString ar = "cp " + fn + " /tmp/copy.dat.gz";
        system(ar.toUtf8().constData());
        QString gu = "gzip -df /tmp/copy.dat.gz";
        system(gu.toUtf8().constData());
        fn = "/tmp/copy.dat";
    }
    if ((fn.size() >= 4 && fn.right(4) == ".dat") ||
        (fn.size() >= 6 && fn.right(4) == ".track")) {
        TrackData nxt = TrackData(fn.toUtf8().constData());
        TrackRestriction rests;
        nxt.calcTimeBounds(rests.time.low, rests.time.high);
        rests.time.high += 0.01 * CLHEP::ns;
        nxt.calcEnergyBounds(rests.energy.low, rests.energy.high);
        rests.energy.high *= 1.2;
        rests.energy.low /= 1.2;
        rests.energy.low = std::max(1.0 * CLHEP::eV, rests.energy.low);
        rests.seqno.low = 1;
        rests.seqno.high = nxt.getNTracks();
        track_res_actual.push_back(rests);
        track_res_bounds.push_back(rests);
        track_options.push_back(nxt);
        reloadChoiceMenus();
    } else {
        qDebug("File not a track data file!");
    }
}

void Viewer::setViewRotation(int sel) {
    if (sel < 6) {
        const CLHEP::Hep3Vector atv[6] = {
            G4ThreeVector(1, 0, 0),  G4ThreeVector(0, 1, 0),
            G4ThreeVector(0, 0, 1),  G4ThreeVector(-1, 0, 0),
            G4ThreeVector(0, -1, 0), G4ThreeVector(0, 0, -1)};
        const CLHEP::Hep3Vector upv[6] = {
            G4ThreeVector(0, 1, 0),  G4ThreeVector(0, 0, 1),
            G4ThreeVector(1, 0, 0),  G4ThreeVector(0, -1, 0),
            G4ThreeVector(0, 0, -1), G4ThreeVector(-1, 0, 0)};
        G4ThreeVector a = atv[sel];
        G4ThreeVector b = upv[sel];
        G4ThreeVector c = atv[sel].cross(upv[sel]);
        vd.orientation = CLHEP::HepRotation(a, b, c).inverse();
        vd.camera = -a * 2 * vd.scale;
    } else {
        G4RotationMatrix rot;
        if (sel == 6) {
            rot = CLHEP::HepRotationX(CLHEP::pi / 4);
        } else {
            rot = CLHEP::HepRotationX(-CLHEP::pi / 4);
        }
        vd.camera =
            (vd.orientation.inverse() * rot.inverse() * vd.orientation) *
            vd.camera;
        vd.orientation = rot * vd.orientation;
    }

    rwidget->rerender();
}
