#include "Viewer.hh"

#include "CustomWidgets.hh"
#include "Overview.hh"
#include "RenderWidget.hh"

#include <G4BooleanSolid.hh>
#include <G4DisplacedSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VSolid.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VisExtent.hh>

#include <QAction>
#include <QActionGroup>
#include <QCheckBox>
#include <QDockWidget>
#include <QDoubleSpinBox>
#include <QGridLayout>
#include <QHeaderView>
#include <QItemSelection>
#include <QListWidget>
#include <QMenuBar>
#include <QPair>
#include <QPushButton>
#include <QResizeEvent>
#include <QSignalMapper>
#include <QTableView>
#include <QTableWidget>
#include <QTreeView>
#include <QVBoxLayout>

static const G4RotationMatrix identityRotation = G4RotationMatrix();
Element convertCreation(const G4VPhysicalVolume *phys,
                        std::map<const G4Material *, int> &idxs,
                        G4RotationMatrix rot = identityRotation,
                        int *counter = NULL) {
    int cc = 0;
    if (!counter) {
        counter = &cc;
    }
    Element m;
    m.name = phys->GetName();

    G4ThreeVector offset = phys->GetFrameTranslation();
    offset = rot.inverse() * offset;
    const G4RotationMatrix &r = (phys->GetFrameRotation() != NULL)
                                    ? *phys->GetFrameRotation()
                                    : identityRotation;
    rot = r * rot;

    // Only identity has a trace of +3 => norm2 of 0
    m.rotated = rot.norm2() > 1e-10;
    m.offset = offset;
    m.rot = rot;

    const G4LogicalVolume *log = phys->GetLogicalVolume();
    const G4Material *mat = log->GetMaterial();
    m.matcode = idxs[mat];
    m.solid = log->GetSolid();
    m.visible = mat->GetDensity() > 0.1 * CLHEP::g / CLHEP::cm3;
    m.alpha = 0.8;
    m.ecode = *counter;
    (*counter)++;

    m.children = std::vector<Element>();
    for (int i = 0; i < log->GetNoDaughters(); i++) {
        m.children.push_back(
            convertCreation(log->GetDaughter(i), idxs, m.rot, counter));
    }
    return m;
}

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

Viewer::Viewer(const std::vector<GeoOption> &options,
               const std::vector<TrackData> &trackopts)
    : QMainWindow() {
    srand(1000 * QTime::currentTime().second() + QTime::currentTime().msec());

    which_geo = 0;
    geo_options = options;
    track_options = trackopts;
    // Construct list of valid materials & lookups
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
    int k = 0;
    vd.matinfo = std::vector<MaterialInfo>();
    for (const QPair<QString, const G4Material *> &p : mlist) {
        MaterialInfo mi;
        mi.mtl = p.second;
        mi.hue = rand() / float(RAND_MAX - 1);
        vd.matinfo.push_back(mi);
        vd.matcode_map[p.second] = k;
        k++;
    }
    vd.elements = convertCreation(geo_options[which_geo].vol, vd.matcode_map);
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

    QMenu *picker = new QMenu("Choose Geometry");
    QActionGroup *opts = new QActionGroup(this);
    for (size_t i = 0; i < geo_options.size(); i++) {
        QString s = QString(geo_options[i].name.c_str());
        QAction *ch = new QAction(s);
        ch->setToolTip(s);
        ch->setCheckable(true);
        opts->addAction(ch);
        if (i == 0) {
            ch->setChecked(true);
        }
        picker->addAction(ch);
    }
    connect(opts, SIGNAL(triggered(QAction *)), this,
            SLOT(changeGeometry(QAction *)));
    QMenu *tpicker = new QMenu("Choose Tracks");
    QActionGroup *topts = new QActionGroup(this);
    QAction *base = new QAction("None");
    base->setToolTip(base->text());
    base->setCheckable(true);
    if (track_options.size() == 0) {
        base->setChecked(true);
    }
    topts->addAction(base);
    tpicker->addAction(base);
    for (size_t i = 0; i < track_options.size(); i++) {
        QString s = QString("%1").arg(i + 1);
        QAction *ch = new QAction(s);
        ch->setToolTip(s);
        ch->setCheckable(true);
        topts->addAction(ch);
        if (i == 0) {
            ch->setChecked(true);
        }
        tpicker->addAction(ch);
    }
    connect(topts, SIGNAL(triggered(QAction *)), this,
            SLOT(changeTracks(QAction *)));

    QAction *clipAction = new QAction("Clipping");
    clipAction->setToolTip("Edit geometry clipping planes and track limits");
    connect(clipAction, SIGNAL(triggered()), this, SLOT(restClip()));

    QAction *treeAction = new QAction("Tree");
    treeAction->setToolTip("View object hierarchy");
    connect(treeAction, SIGNAL(triggered()), this, SLOT(restTree()));

    QAction *infoAction = new QAction("Info");
    infoAction->setToolTip("Misc information on selected volume");
    connect(infoAction, SIGNAL(triggered()), this, SLOT(restInfo()));

    QAction *rayAction = new QAction("Ray");
    rayAction->setToolTip("List currently-hovered-over visible objects");
    connect(rayAction, SIGNAL(triggered()), this, SLOT(restRay()));

    QAction *mtlAction = new QAction("Material");
    mtlAction->setToolTip("Control material view properties");
    connect(mtlAction, SIGNAL(triggered()), this, SLOT(restMtl()));

    QAction *screenAction = new QAction("Screenshot");
    screenAction->setToolTip("Take a screenshot of active scene");
    connect(screenAction, SIGNAL(triggered()), this, SLOT(screenshot()));

    QAction *screen4Action = new QAction("ScrSht4x");
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
    this->menuBar()->addMenu(picker);
    this->menuBar()->addMenu(tpicker);
    // TODO: checkbox by visibility
    QMenu *sub = this->menuBar()->addMenu("View");
    sub->addAction(clipAction);
    sub->addAction(treeAction);
    sub->addAction(infoAction);
    sub->addAction(rayAction);
    sub->addAction(mtlAction);
    this->menuBar()->addSeparator();
    this->menuBar()->addAction(screenAction);
    this->menuBar()->addAction(screen4Action);

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
    times_lower = new QDoubleSpinBox();
    times_upper = new QDoubleSpinBox();
    times_lower->setSuffix("ns");
    times_upper->setSuffix("ns");
    times_lower->setDecimals(3);
    times_upper->setDecimals(3);
    times_lower->setSingleStep(0.1);
    times_upper->setSingleStep(0.1);
    energy_lower = new ExpoSpinBox();
    energy_upper = new ExpoSpinBox();
    energy_lower->setSingleStep(10);
    energy_upper->setSingleStep(10);
    energy_lower->setSuffix("eV");
    energy_upper->setSuffix("eV");
    count_lower = new QSpinBox();
    count_upper = new QSpinBox();
    if (which_tracks > 0) {
        const TrackRestriction &res = track_res_actual[which_tracks - 1];
        times_lower->setRange(res.time.low / CLHEP::ns,
                              res.time.high / CLHEP::ns);
        times_upper->setRange(res.time.low / CLHEP::ns,
                              res.time.high / CLHEP::ns);
        times_lower->setValue(res.time.low / CLHEP::ns);
        times_upper->setValue(res.time.high / CLHEP::ns);
        int el = energy_lower->nearestIntFromExp(res.energy.low / CLHEP::eV);
        int eh = energy_lower->nearestIntFromExp(res.energy.high / CLHEP::eV);
        energy_lower->setRange(el, eh);
        energy_upper->setRange(el, eh);
        energy_lower->setValue(el);
        energy_upper->setValue(eh);
        count_lower->setRange(1, int(res.seqno.high));
        count_upper->setRange(1, int(res.seqno.high));
        count_lower->setValue(0);
        count_upper->setValue(int(res.seqno.high));
        int step = exp10(
            std::max(0, int(std::floor(std::log10(res.seqno.high / 15.)))));
        count_lower->setSingleStep(step);
        count_upper->setSingleStep(step);
    } else {
        times_lower->setDisabled(true);
        times_upper->setDisabled(true);
        energy_lower->setRange(0, 900 - 1);
        energy_upper->setRange(0, 900 - 1);
        energy_lower->setDisabled(true);
        energy_upper->setDisabled(true);
        count_lower->setDisabled(true);
        count_upper->setDisabled(true);
    }
    connect(count_lower, SIGNAL(valueChanged(int)), this, SLOT(updatePlanes()));
    connect(count_upper, SIGNAL(valueChanged(int)), this, SLOT(updatePlanes()));
    connect(times_lower, SIGNAL(valueChanged(double)), this,
            SLOT(updatePlanes()));
    connect(times_upper, SIGNAL(valueChanged(double)), this,
            SLOT(updatePlanes()));
    connect(energy_lower, SIGNAL(valueChanged(int)), this,
            SLOT(updatePlanes()));
    connect(energy_upper, SIGNAL(valueChanged(int)), this,
            SLOT(updatePlanes()));
    QHBoxLayout *trb = new QHBoxLayout();
    trb->addWidget(times_lower);
    trb->addWidget(times_upper);
    vb->addLayout(trb, 0);
    QHBoxLayout *erb = new QHBoxLayout();
    erb->addWidget(energy_lower);
    erb->addWidget(energy_upper);
    vb->addLayout(erb, 0);
    QHBoxLayout *crb = new QHBoxLayout();
    crb->addWidget(count_lower);
    crb->addWidget(count_upper);
    vb->addLayout(crb, 0);
    vb->addStretch(1);
    cont->setLayout(vb);
    dock_clip->setWidget(cont);

    // Tree view (with vis/novis, hue control
    dock_tree = new QDockWidget("Element tree", this);
    tree_view = new QTreeView();
    tree_model = new OverView(vd);
    tree_view->setSortingEnabled(false);
    tree_view->setModel(tree_model);
    tree_view->setRootIndex(QModelIndex());
    tree_view->setSelectionMode(QAbstractItemView::SingleSelection);
    tree_view->setSelectionBehavior(QAbstractItemView::SelectRows);
    tree_view->setIndentation(20);
    tree_view->collapseAll();
    tree_view->header()->resizeSections(QHeaderView::ResizeToContents);
    tree_view->header()->setStretchLastSection(true);
    tree_view->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    tree_view->setItemDelegateForColumn(2, new AlphaBoxDelegate(tree_model));
    connect(tree_view, SIGNAL(activated(const QModelIndex &)), tree_model,
            SLOT(respToActive(const QModelIndex &)));
    connect(tree_model, SIGNAL(colorChange()), rwidget, SLOT(rerender()));

    connect(tree_view->selectionModel(),
            SIGNAL(selectionChanged(const QItemSelection &,
                                    const QItemSelection &)),
            tree_model, SLOT(respToSelection(const QItemSelection &,
                                             const QItemSelection &)));
    connect(tree_model, SIGNAL(selectedElement(Element *)), this,
            SLOT(indicateElement(Element *)));
    dock_tree->setWidget(tree_view);

    // Info for a particular item (mass, volume, mtl specs; table...)
    dock_info = new QDockWidget("Element info", this);
    info_table = new QTableWidget();
    info_table->setColumnCount(1);
    // TODO: tooltips, require a model & headerData()
    QStringList keys;
    keys << "Name"
         << "Material"
         << "Density"
         << "Volume"
         << "Mass"
         << "Surface Area"
         << "Bool roots"
         << "Bool depth"
         << "Bool splits";
    info_table->setEditTriggers(QTableWidget::NoEditTriggers);
    info_table->setRowCount(keys.size());
    info_table->setVerticalHeaderLabels(keys);
    QStringList kv = QStringList() << "Value";
    info_table->setHorizontalHeaderLabels(kv);
    info_table->horizontalHeader()->setStretchLastSection(true);
    for (int i = 0; i < keys.size(); i++) {
        info_table->setItem(i, 0, new QTableWidgetItem(""));
    }
    dock_info->setWidget(info_table);

    dock_ray = new QDockWidget("Ray table", this);
    ray_table = new QListWidget();
    dock_ray->setWidget(ray_table);

    dock_mtl = new QDockWidget("Materials", this);
    QWidget *mtlc = new QWidget();
    QVBoxLayout *mla = new QVBoxLayout();
    mtl_divchk = new QCheckBox("Split by material");
    mtl_divchk->setCheckState(vd.split_by_material ? Qt::Checked
                                                   : Qt::Unchecked);
    connect(mtl_divchk, SIGNAL(stateChanged(int)), this,
            SLOT(updateMaterials()));
    mtl_showlines = new QCheckBox("Render lines");
    mtl_showlines->setCheckState(Qt::Unchecked);
    connect(mtl_showlines, SIGNAL(stateChanged(int)), this,
            SLOT(updateShowLines()));
    mtl_table = new QTableView();
    mtl_table->setSelectionBehavior(QAbstractItemView::SelectItems);
    mtl_table->setSelectionMode(QAbstractItemView::NoSelection);
    MaterialModel *mmod = new MaterialModel(vd);
    mtl_table->setModel(mmod);
    mtl_table->horizontalHeader()->setStretchLastSection(true);
    mtl_table->setItemDelegateForColumn(0, new HueSpinBoxDelegate(mmod));
    connect(mmod, SIGNAL(colorChange()), this, SLOT(updateMaterials()));
    mla->addWidget(mtl_divchk);
    mla->addWidget(mtl_showlines);
    mla->addWidget(mtl_table);
    mtlc->setLayout(mla);
    dock_mtl->setWidget(mtlc);

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

    addDockWidget(Qt::LeftDockWidgetArea, dock_mtl);
    dock_mtl->setAllowedAreas(Qt::AllDockWidgetAreas);
    dock_mtl->setFeatures(feat);
    dock_mtl->setVisible(false);

    setMouseTracking(true);

    this->show();
}

Viewer::~Viewer() {}

void Viewer::restClip() { dock_clip->setVisible(true); }
void Viewer::restTree() { dock_tree->setVisible(true); }
void Viewer::restInfo() { dock_info->setVisible(true); }
void Viewer::restRay() { dock_ray->setVisible(true); }
void Viewer::restMtl() { dock_mtl->setVisible(true); }

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
    for (int j = 0; j < m; j++) {
        QString name(elems[j]->name.data());
        ray_table->addItem(name);
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
        double tlow = times_lower->value() * CLHEP::ns;
        double thigh = times_upper->value() * CLHEP::ns;
        double elow =
            energy_lower->expFromInt(energy_lower->value()) * CLHEP::eV;
        double ehigh =
            energy_upper->expFromInt(energy_upper->value()) * CLHEP::eV;
        size_t nlow = size_t(count_lower->value());
        size_t nhigh = size_t(count_upper->value());

        TrackRestriction &res = track_res_actual[which_tracks - 1];
        res.energy = {std::min(elow, ehigh), std::max(elow, ehigh)};
        res.time = {std::min(tlow, thigh), std::max(tlow, thigh)};
        res.seqno = {std::min(nlow, nhigh), std::max(nlow, nhigh)};
        trackdata = TrackData(track_options[which_tracks - 1], vd, res.time,
                              res.energy, res.seqno);
    } else {
        // Q: how to pull QSharedData on the Elements as well
        trackdata = TrackData();
    }
    rwidget->rerender();
}

void Viewer::updateMaterials() {
    vd.split_by_material = mtl_divchk->isChecked();
    // Model autoupdates others
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
                vd.elements =
                    convertCreation(geo_options[which_geo].vol, vd.matcode_map);
                vd.scene_radius =
                    vd.elements.solid->GetExtent().GetExtentRadius();
                if (4 * vd.scene_radius > vd.camera.mag()) {
                    vd.camera *= 4 * vd.scene_radius / vd.camera.mag();
                }
                tree_model->recalculate();
                tree_view->collapseAll();
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
        times_lower->blockSignals(true);
        times_upper->blockSignals(true);
        energy_lower->blockSignals(true);
        energy_upper->blockSignals(true);
        count_lower->blockSignals(true);
        count_upper->blockSignals(true);
        const TrackRestriction &res = track_res_bounds[which_tracks - 1];
        const TrackRestriction &cur = track_res_actual[which_tracks - 1];
        times_lower->setRange(res.time.low / CLHEP::ns,
                              res.time.high / CLHEP::ns);
        times_upper->setRange(res.time.low / CLHEP::ns,
                              res.time.high / CLHEP::ns);
        times_lower->setValue(cur.time.low / CLHEP::ns);
        times_upper->setValue(cur.time.high / CLHEP::ns);
        int el = energy_lower->nearestIntFromExp(res.energy.low / CLHEP::eV);
        int eh = energy_lower->nearestIntFromExp(res.energy.high / CLHEP::eV);
        int cl = energy_lower->nearestIntFromExp(cur.energy.low / CLHEP::eV);
        int ch = energy_lower->nearestIntFromExp(cur.energy.high / CLHEP::eV);
        energy_lower->setRange(el, eh);
        energy_upper->setRange(el, eh);
        energy_lower->setValue(cl);
        energy_upper->setValue(ch);
        count_lower->setRange(1, int(res.seqno.high));
        count_upper->setRange(1, int(res.seqno.high));
        count_lower->setValue(int(cur.seqno.low));
        count_upper->setValue(int(cur.seqno.high));
        int step = exp10(
            std::max(0, int(std::floor(std::log10(res.seqno.high / 15.)))));
        count_lower->setSingleStep(step);
        count_upper->setSingleStep(step);
        times_lower->blockSignals(false);
        times_upper->blockSignals(false);
        energy_lower->blockSignals(false);
        energy_upper->blockSignals(false);
        count_lower->blockSignals(false);
        count_upper->blockSignals(false);
    }
    times_lower->setEnabled(active);
    times_upper->setEnabled(active);
    energy_lower->setEnabled(active);
    energy_upper->setEnabled(active);
    count_lower->setEnabled(active);
    count_upper->setEnabled(active);
    updatePlanes();
}

void calculateBooleanProperties(const G4VSolid *sol,
                                QSet<const G4VSolid *> &roots, int &treedepth,
                                int &nbooleans, int depth = 0) {
    const G4BooleanSolid *b = dynamic_cast<const G4BooleanSolid *>(sol);
    treedepth = std::max(treedepth, depth);
    if (b) {
        calculateBooleanProperties(b->GetConstituentSolid(0), roots, treedepth,
                                   nbooleans, depth + 1);
        calculateBooleanProperties(b->GetConstituentSolid(1), roots, treedepth,
                                   nbooleans, depth + 1);
        nbooleans++;
        return;
    }
    const G4DisplacedSolid *d = dynamic_cast<const G4DisplacedSolid *>(sol);
    if (d) {
        calculateBooleanProperties(d->GetConstituentMovedSolid(), roots,
                                   treedepth, nbooleans, depth + 1);
        return;
    }
    roots.insert(sol);
}

void Viewer::indicateElement(Element *e) {
    if (!e) {
        for (int i = 0; i < info_table->rowCount(); i++) {
            info_table->item(i, 0)->setText("");
        }
    } else {
        // fill info table
        info_table->item(0, 0)->setText(e->name.c_str());
        const G4Material *mat = vd.matinfo[e->matcode].mtl;
        info_table->item(1, 0)->setText(mat->GetName().c_str());
        G4double dens = mat->GetDensity();
        info_table->item(2, 0)->setText(
            QString::number(dens / (CLHEP::g / CLHEP::cm3), 'g', 4) + " g/cm3");
        G4double volume = e->solid->GetCubicVolume();
        info_table->item(3, 0)->setText(
            QString::number(volume / CLHEP::cm3, 'g', 4) + " cm3");
        info_table->item(4, 0)->setText(
            QString::number(volume * dens / CLHEP::kg, 'g', 4) + " kg");
        G4double surf = e->solid->GetSurfaceArea();
        info_table->item(5, 0)->setText(
            QString::number(surf / CLHEP::cm2, 'g', 4) + " cm2");

        QSet<const G4VSolid *> roots;
        int treedepth = 0;
        int nbooleans = 0;
        calculateBooleanProperties(e->solid, roots, treedepth, nbooleans);
        info_table->item(6, 0)->setText(QString::number(roots.size()));
        info_table->item(7, 0)->setText(QString::number(treedepth));
        info_table->item(8, 0)->setText(QString::number(nbooleans));
    }
}

void Viewer::updateShowLines() {
    rwidget->setFullDetail(mtl_showlines->isChecked());
}

void Viewer::screenshot(int sx) {
    int w = rwidget->width() * sx, h = rwidget->height() * sx;
    RenderSaveObject *rso = new RenderSaveObject(vd, trackdata, w, h);
    rso->start();
}
