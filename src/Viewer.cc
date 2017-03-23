#include "Viewer.hh"

#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VSolid.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VisExtent.hh>

#include <QAction>
#include <QActionGroup>
#include <QComboBox>
#include <QDir>
#include <QDockWidget>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QGridLayout>
#include <QHeaderView>
#include <QImage>
#include <QImageWriter>
#include <QItemSelection>
#include <QListWidget>
#include <QMenuBar>
#include <QMessageBox>
#include <QPainter>
#include <QProgressDialog>
#include <QPushButton>
#include <QResizeEvent>
#include <QSignalMapper>
#include <QTableWidget>
#include <QThread>
#include <QTime>
#include <QTreeView>
#include <QVBoxLayout>

RenderWidget::RenderWidget(ViewData &v) : QWidget(), currView(v) {
    setAttribute(Qt::WA_OpaquePaintEvent, true);

    back = QImage(50, 50, QImage::Format_RGB32);
    {
        QPainter q(&back);
        q.fillRect(back.rect(), Qt::blue);
        q.fillRect(QRect(5, 5, 25, 25), Qt::red);
    }

    t = QVector<QThread *>();
    w = QVector<RenderWorker *>();
    qRegisterMetaType<ViewData>("ViewData");
    qRegisterMetaType<QImage *>("QImage*");
    int itc = QThread::idealThreadCount();
    itc = itc < 0 ? 2 : itc;
#if 0
    itc = 1; // Debug override
#endif
    for (int i = 0; i < itc; i++) {
        QThread *tt = new QThread();
        RenderWorker *ww = new RenderWorker();
        ww->moveToThread(tt);
        connect(ww, SIGNAL(completed()), this, SLOT(completed()));
        connect(ww, SIGNAL(aborted()), this, SLOT(aborted()));
        tt->start();
        t.append(tt);
        w.append(ww);
    }
    state = NONE;
    last_level_of_detail = 10000;
}

RenderWidget::~RenderWidget() {}

void RenderWidget::rerender() {
    currView.level_of_detail = 2;
    rerender_priv();
}

void RenderWidget::rerender_priv() {
    if (currView.level_of_detail > last_level_of_detail && state != NONE) {
        // don't know which workers have already completed, so abort all,
        // and queue an abort-clearer on all that ensures the next run
        // isn't aborted
        for (RenderWorker *ww : w) {
            ww->abort_task = true;
            QMetaObject::invokeMethod(ww, "flushAbort", Qt::QueuedConnection);
        }
// current response count indicates completed workers
// the aborted ones should bring the total up
#if 0
        qDebug("abort: response count %d", response_count);
#endif
        state = ACTIVE_AND_QUEUED;
        last_level_of_detail = 10000;
        return;
    }
    if (state == ACTIVE || state == ACTIVE_AND_QUEUED) {
        state = ACTIVE_AND_QUEUED;
        return;
    }

    int scl = int(std::pow(3.0, std::max(currView.level_of_detail, 0)));
    if (currView.level_of_detail == 0) {
        request_time = QTime::currentTime();
        arrived = aReqd;
    }
    next =
        QImage(this->width() / scl, this->height() / scl, QImage::Format_RGB32);
    for (int i = 0; i < w.size(); i++) {
        QMetaObject::invokeMethod(
            w[i], "render", Qt::QueuedConnection, Q_ARG(ViewData, currView),
            Q_ARG(QImage *, &next), Q_ARG(int, i), Q_ARG(int, w.size()));
    }
    response_count = 0;

    last_level_of_detail = currView.level_of_detail;
    if (currView.level_of_detail <= 0) {
        state = ACTIVE;
    } else {
        state = ACTIVE_AND_QUEUED;
        currView.level_of_detail--;
    }
}

void RenderWidget::completed() {
    if (response_count < w.size()) {
        response_count++;
    }
    if (response_count < w.size()) {
        return;
    }
    back = next;
    if (state == ACTIVE) {
        state = NONE;
    } else if (state == ACTIVE_AND_QUEUED) {
        state = NONE;
        rerender_priv();
    }
    if (arrived == aReqd && back.size() == this->size()) {
        arrived = aCompl;
    }
    this->repaint();
}
void RenderWidget::aborted() {
    if (response_count < w.size()) {
        response_count++;
    }
    if (response_count < w.size()) {
        return;
    }
    state = NONE;
    rerender_priv();
}

void RenderWidget::resizeEvent(QResizeEvent *) {
    currView.level_of_detail = 2;
    rerender_priv();
}

void RenderWidget::paintEvent(QPaintEvent *) {
    QPainter q(this);
    QRect sz = back.rect();
    QImage mvd;
    if (sz.height() == this->height() && sz.width() == this->width()) {
        mvd = back;
    } else if (this->width() * sz.height() >= this->height() * sz.width()) {
        mvd = back.scaledToHeight(this->height(), Qt::FastTransformation);
    } else {
        mvd = back.scaledToWidth(this->width(), Qt::FastTransformation);
    }
    q.drawImage(
        this->rect().center() - QPoint(mvd.width() / 2, mvd.height() / 2), mvd);
    if (arrived == aCompl) {
        arrived = aThere;
#if 0
        qDebug("img completed after %d ms",
               request_time.msecsTo(QTime::currentTime()));
#endif
    }
}

PlaneEdit::PlaneEdit(Plane p) {
    nx = new QDoubleSpinBox();
    ny = new QDoubleSpinBox();
    nz = new QDoubleSpinBox();
    d = new QDoubleSpinBox();
    nx->setRange(-10, 10);
    ny->setRange(-10, 10);
    nz->setRange(-10, 10);
    d->setRange(-100, 100);
    d->setDecimals(1);
    nx->setValue(p.normal.x());
    ny->setValue(p.normal.y());
    nz->setValue(p.normal.z());
    d->setValue(p.offset);

    act = new QPushButton("ON");
    act->setCheckable(true);
    connect(act, SIGNAL(toggled(bool)), this, SLOT(setActive(bool)));

    unit = new QComboBox();
    unit->addItem("um");
    unit->addItem("mm");
    unit->addItem("cm");
    unit->addItem("m");
    unit->setCurrentIndex(2);

    connect(nx, SIGNAL(valueChanged(double)), this, SIGNAL(updated()));
    connect(ny, SIGNAL(valueChanged(double)), this, SIGNAL(updated()));
    connect(nz, SIGNAL(valueChanged(double)), this, SIGNAL(updated()));
    connect(d, SIGNAL(valueChanged(double)), this, SIGNAL(updated()));
    connect(act, SIGNAL(toggled(bool)), this, SIGNAL(updated()));
    connect(unit, SIGNAL(currentIndexChanged(int)), this, SIGNAL(updated()));

    QGridLayout *grid = new QGridLayout();
    grid->addWidget(nx, 0, 0);
    grid->addWidget(ny, 0, 1);
    grid->addWidget(nz, 0, 2);
    grid->addWidget(d, 1, 0);
    grid->addWidget(unit, 1, 1);
    grid->addWidget(act, 1, 2);

    this->setLayout(grid);
    setActive(false);
}

PlaneEdit::~PlaneEdit() {}

void PlaneEdit::setActive(bool active) {
    nx->setEnabled(active);
    ny->setEnabled(active);
    nz->setEnabled(active);
    d->setEnabled(active);
    unit->setEnabled(active);
    act->setText(active ? "OFF" : "ON");
}

Plane PlaneEdit::getPlane() {
    Plane p;
    p.normal = act->isChecked()
                   ? G4ThreeVector(nx->value(), ny->value(), nz->value())
                   : G4ThreeVector();
    p.offset = d->value();
    return p;
}

ViewData initVD(const Element &root) {
    ViewData d;
    d.root = root;
    d.scene_radius = d.root.solid->GetExtent().GetExtentRadius();
    d.scale = 2 * d.scene_radius;
    d.camera = G4ThreeVector(-4 * d.scene_radius, 0, 0);
    d.orientation = G4RotationMatrix();
    d.level_of_detail = 0; // 0 is full; 1 is 1/9, 2 is 1/81; depends on timing
    d.clipping_planes = std::vector<Plane>();
    return d;
}

Element convertCreation(G4VPhysicalVolume *phys) {
    Element m;
    m.name = phys->GetName();

    m.rotated = phys->GetFrameRotation() != NULL;
    m.offset = phys->GetFrameTranslation();
    m.rot = phys->GetFrameRotation() ? *phys->GetFrameRotation()
                                     : G4RotationMatrix();

    m.mat = phys->GetLogicalVolume()->GetMaterial();
    m.solid = phys->GetLogicalVolume()->GetSolid();
    m.visible = m.mat->GetDensity() > 0.1 * CLHEP::g / CLHEP::cm3;
    m.hue = reinterpret_cast<long>(m.mat) % 16384 / 16384.0;
    m.alpha = 0.8;
    m.ngeocalls = 0;

    m.children = std::vector<Element>();
    for (int i = 0; i < phys->GetLogicalVolume()->GetNoDaughters(); i++) {
        m.children.push_back(
            convertCreation(phys->GetLogicalVolume()->GetDaughter(i)));
    }
    return m;
}

HueSpinBoxDelegate::HueSpinBoxDelegate(OverView *model, QObject *parent)
    : QItemDelegate(parent) {
    oneTrueModel = model;
}

QWidget *HueSpinBoxDelegate::createEditor(QWidget *parent,
                                          const QStyleOptionViewItem &,
                                          const QModelIndex &idx) const {
    QDoubleSpinBox *editor = new QDoubleSpinBox(parent);
    // TODO: override virts. dbl valueFromText(QStr), QStr textFromValue(dbl)
    // to prevent 0 being added to the tail realtime. (Alternatively;
    // insert mode editing
    editor->setDecimals(2);
    editor->setRange(0., 1.);
    editor->setSingleStep(0.05);
    editor->setWrapping(true);
    editor->setStyleSheet("QDoubleSpinBox { color : blue; } ");
    // Insert extra parameters in unused fields
    editor->setToolTipDuration(
        int(reinterpret_cast<long>(idx.internalPointer())));
    editor->setMaximumHeight(1000 + idx.row());
    QSignalMapper *q = new QSignalMapper(editor);
    connect(editor, SIGNAL(valueChanged(double)), q, SLOT(map()));
    q->setMapping(editor, editor);
    connect(q, SIGNAL(mapped(QWidget *)), oneTrueModel,
            SLOT(hueUpdate(QWidget *)));
    return editor;
}

void HueSpinBoxDelegate::setEditorData(QWidget *editor,
                                       const QModelIndex &index) const {
    QDoubleSpinBox *spinBox = static_cast<QDoubleSpinBox *>(editor);
    spinBox->blockSignals(true);
    spinBox->setValue(index.model()->data(index, Qt::EditRole).toDouble());
    spinBox->blockSignals(false);
}
void HueSpinBoxDelegate::setModelData(QWidget *editor,
                                      QAbstractItemModel *model,
                                      const QModelIndex &index) const {
    QDoubleSpinBox *spinBox = static_cast<QDoubleSpinBox *>(editor);
    spinBox->interpretText();
    model->setData(index, spinBox->value(), Qt::EditRole);
}

void HueSpinBoxDelegate::updateEditorGeometry(
    QWidget *editor, const QStyleOptionViewItem &option,
    const QModelIndex &) const {
    editor->setGeometry(option.rect);
}

AlphaBoxDelegate::AlphaBoxDelegate(OverView *model, QObject *parent)
    : QItemDelegate(parent) {
    oneTrueModel = model;
}

QWidget *AlphaBoxDelegate::createEditor(QWidget *parent,
                                        const QStyleOptionViewItem &,
                                        const QModelIndex &idx) const {
    QDoubleSpinBox *editor = new QDoubleSpinBox(parent);
    editor->setDecimals(2);
    editor->setRange(0., 1.);
    editor->setSingleStep(0.2);
    // Insert extra parameters in unused fields
    editor->setToolTipDuration(
        int(reinterpret_cast<long>(idx.internalPointer())));
    editor->setMaximumHeight(1000 + idx.row());
    QSignalMapper *q = new QSignalMapper(editor);
    connect(editor, SIGNAL(valueChanged(double)), q, SLOT(map()));
    q->setMapping(editor, editor);
    connect(q, SIGNAL(mapped(QWidget *)), oneTrueModel,
            SLOT(alphaUpdate(QWidget *)));
    return editor;
}

void AlphaBoxDelegate::setEditorData(QWidget *editor,
                                     const QModelIndex &index) const {
    QDoubleSpinBox *spinBox = static_cast<QDoubleSpinBox *>(editor);
    spinBox->blockSignals(true);
    spinBox->setValue(index.model()->data(index, Qt::EditRole).toDouble());
    spinBox->blockSignals(false);
}
void AlphaBoxDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
                                    const QModelIndex &index) const {
    QDoubleSpinBox *spinBox = static_cast<QDoubleSpinBox *>(editor);
    spinBox->interpretText();
    model->setData(index, spinBox->value(), Qt::EditRole);
}

void AlphaBoxDelegate::updateEditorGeometry(QWidget *editor,
                                            const QStyleOptionViewItem &option,
                                            const QModelIndex &) const {
    editor->setGeometry(option.rect);
}

OverView::OverView(ViewData &c) : currView(c) { recalculate(); }

OverView::~OverView() {}

QModelIndex OverView::index(int r, int c, const QModelIndex &p) const {
    if (!hasIndex(r, c, p))
        return QModelIndex();

    // Note: root item is _above_ idx=0, has invalid index, etc.
    if (p.isValid()) {
        const Node &n = link[int(reinterpret_cast<long>(p.internalPointer()))];
        if (r < n.lexi.size()) {
            int cidx = n.sub[n.lexi[r]];
            return createIndex(r, c, reinterpret_cast<void *>(long(cidx)));
        } else {
            return QModelIndex();
        }
    } else {
        if (r == 0) {
            return createIndex(0, c, reinterpret_cast<void *>(long(0)));
        }
        return QModelIndex();
    }
}
QModelIndex OverView::parent(const QModelIndex &index) const {
    if (!index.isValid()) {
        return QModelIndex();
    }
    int idx = int(reinterpret_cast<long>(index.internalPointer()));
    if (idx == 0) {
        return QModelIndex();
    }
    int pidx = link[idx].parent;
    if (pidx == 0) {
        return createIndex(0, 0, reinterpret_cast<void *>(long(0)));
    }
    int gidx = link[pidx].parent;
    for (int j = 0; j < link[gidx].lexi.size(); j++) {
        int sidx = link[gidx].sub[int(link[gidx].lexi[j])];
        if (sidx == pidx) {
            return createIndex(j, 0, reinterpret_cast<void *>(long(pidx)));
        }
    }
    return QModelIndex();
}
int OverView::rowCount(const QModelIndex &index) const {
    if (index.isValid()) {
        if (index.column() > 0) {
            return 0;
        }

        const Node &n =
            link[int(reinterpret_cast<long>(index.internalPointer()))];
        return n.sub.size();
    } else {
        return 1;
    }
}
int OverView::columnCount(const QModelIndex &) const { return 4; }
QVariant OverView::headerData(int section, Qt::Orientation orientation,
                              int role) const {
    if (orientation == Qt::Vertical) {
        return QVariant();
    }
    if (role == Qt::DisplayRole) {
        if (section == 0) {
            return QVariant("Name");
        } else if (section == 1) {
            return QVariant("Vis");
        } else if (section == 2) {
            return QVariant("Hue");
        } else if (section == 3) {
            return QVariant("Alpha");
        }
    }
    return QVariant();
}
Qt::ItemFlags OverView::flags(const QModelIndex &index) const {
    if (!index.isValid()) {
        return 0;
    }

    Qt::ItemFlags base = Qt::ItemIsSelectable | Qt::ItemIsEnabled;
    if (index.column() == 1) {
        // ItemIsEnabled = 32,
        // ItemIsAutoTristate = 64, govern state
        return base | Qt::ItemIsUserCheckable;
    }
    if (index.column() == 2 || index.column() == 3) {
        return base | Qt::ItemIsEditable;
    }
    return base;
}
QVariant OverView::data(const QModelIndex &index, int role) const {
    if (!index.isValid()) {
        return QVariant();
    }
    const Node &n = link[int(reinterpret_cast<long>(index.internalPointer()))];

    // Assume only valid indices...
    if (role == Qt::DisplayRole) {
        if (index.column() == 0) {
            return QVariant(QString(n.elem->name.c_str()));
        }
        if (index.column() == 2) {
            return QVariant(QString::number(n.elem->hue, 'f', 2));
        }
        if (index.column() == 3) {
            return QVariant(QString::number(n.elem->alpha, 'f', 2));
        }
    }
    if (role == Qt::CheckStateRole && index.column() == 1) {
        return QVariant(n.elem->visible ? Qt::Checked : Qt::Unchecked);
    }
    if (role == Qt::EditRole) {
        if (index.column() == 2) {
            return QVariant(n.elem->hue);
        }
        if (index.column() == 3) {
            return QVariant(n.elem->alpha);
        }
    }

    return QVariant();
}
bool OverView::setData(const QModelIndex &index, const QVariant &value,
                       int role) {
    if (!index.isValid()) {
        return false;
    }
    const Node &n = link[int(reinterpret_cast<long>(index.internalPointer()))];

    if (role == Qt::EditRole) {
        if (index.column() == 2) {
            n.elem->hue = value.toDouble();
            QVector<int> roles;
            roles.push_back(Qt::DisplayRole);
            emit dataChanged(index, index, roles);
            emit colorChange();
            return true;
        }
        if (index.column() == 3) {
            n.elem->alpha = value.toDouble();
            QVector<int> roles;
            roles.push_back(Qt::DisplayRole);
            emit dataChanged(index, index, roles);
            emit colorChange();
            return true;
        }
    }
    return false;
}
void OverView::respToActive(const QModelIndex &index) {
    if (index.column() == 1) {
        int idx = int(reinterpret_cast<long>(index.internalPointer()));
        link[idx].elem->visible ^= true;
        QVector<int> roles;
        roles.push_back(Qt::CheckStateRole);
        emit dataChanged(index, index, roles);
        emit colorChange();
    }
}
void OverView::hueUpdate(QWidget *w) {
    QDoubleSpinBox *e = static_cast<QDoubleSpinBox *>(w);
    if (!e) {
        return;
    }
    int idx = e->toolTipDuration();
    int row = e->maximumHeight() - 1000;
    link[idx].elem->hue = e->value();
    QModelIndex index =
        createIndex(row, 2, reinterpret_cast<void *>(long(idx)));
    QVector<int> roles;
    roles.push_back(Qt::EditRole);
    emit dataChanged(index, index, roles);
    emit colorChange();
}
void OverView::alphaUpdate(QWidget *w) {
    QDoubleSpinBox *e = static_cast<QDoubleSpinBox *>(w);
    if (!e) {
        return;
    }
    int idx = e->toolTipDuration();
    int row = e->maximumHeight() - 1000;
    link[idx].elem->alpha = e->value();
    QModelIndex index =
        createIndex(row, 3, reinterpret_cast<void *>(long(idx)));
    QVector<int> roles;
    roles.push_back(Qt::EditRole);
    emit dataChanged(index, index, roles);
    emit colorChange();
}

void OverView::respToSelection(const QItemSelection &ns,
                               const QItemSelection &) {
    QModelIndexList idxs = ns.indexes();
    if (idxs.size() && idxs[0].isValid()) {
        // Target first item only
        const Node &n =
            link[int(reinterpret_cast<long>(idxs[0].internalPointer()))];
        emit selectedElement(n.elem);
    }
}

void OverView::recalculate() {
    beginResetModel();

    link = QVector<Node>();
    QVector<int> stack;
    stack.push_back(0);
    while (stack.size()) {
        Node n;
        n.eaddr = stack;
        QVector<int> chain;
        if (!link.size()) {
            n.parent = 0; // parent is self :-)
            n.elem = &currView.root;
            chain.push_back(0);
        } else {
            chain.push_back(0);
            for (int i = 1; i < stack.length() - 1; i++) {
                int o = stack[i];
                const Node &cur = link[chain[chain.size() - 1]];
                chain.push_back(cur.sub[cur.lexi[o]]);
            }
            n.parent = chain[chain.length() - 1];
            int realidx = stack[stack.length() - 1];
            n.elem = &link[n.parent].elem->children[size_t(realidx)];
            link[n.parent].sub[link[n.parent].lexi[realidx]] = link.size();
            chain.push_back(link.size());
        }

        n.lexi = QVector<int>();
        n.sub = QVector<int>();

        size_t nkids = n.elem->children.size();
        if (!nkids) {
            // Roll back if no kids
            while (chain.size()) {
                int npos = stack.last() + 1;
                stack.pop_back();
                chain.pop_back();
                if (!chain.size()) {
                    break;
                }
                if (npos < link[chain.last()].sub.size()) {
                    stack.push_back(npos);
                    break;
                }
            }
            link.push_back(n);
            continue;
        }

        // Determine lexical order for kids // todo: get working...
        QStringList names;
        for (size_t j = 0; j < nkids; j++) {
            names << QString(n.elem->children[j].name.c_str());
            n.sub.push_back(-1); // placeholder
        }
        names.sort(Qt::CaseInsensitive);
        // assuming no duplicate names... (todo: a true indirect sort..)
        for (int i = 0; i < int(nkids); i++) {
            for (size_t j = 0; j < nkids; j++) {
                if (names[i] == QString(n.elem->children[j].name.c_str())) {
                    n.lexi.push_back(i);
                    break;
                }
            }
        }
        // Visit first child
        stack.push_back(0);
        link.push_back(n);
    }
    // Reconstruct
    endResetModel();
}

Viewer::Viewer(std::vector<GeoOption> options, size_t idx) : QMainWindow() {
    which_geo = idx;
    geo_options = options;
    if (!geo_options[which_geo].cache) {
        geo_options[which_geo].cache = geo_options[which_geo].cons->Construct();
    }
    vd = initVD(convertCreation(geo_options[which_geo].cache));

    QMenu *picker = new QMenu("Choose Geometry");
    QActionGroup *opts = new QActionGroup(this);
    for (size_t i = 0; i < geo_options.size(); i++) {
        QString s = QString(geo_options[i].name.c_str());
        QAction *ch = new QAction(s);
        ch->setToolTip(s);
        ch->setCheckable(true);
        opts->addAction(ch);
        if (i == idx) {
            ch->setChecked(true);
        }
        picker->addAction(ch);
    }
    connect(opts, SIGNAL(triggered(QAction *)), this,
            SLOT(changeGeometry(QAction *)));

    QAction *clipAction = new QAction("Clipping");
    clipAction->setToolTip("Edit clipping planes");
    connect(clipAction, SIGNAL(triggered()), this, SLOT(restClip()));

    QAction *treeAction = new QAction("Tree");
    treeAction->setToolTip("View tree");
    connect(treeAction, SIGNAL(triggered()), this, SLOT(restTree()));

    QAction *infoAction = new QAction("Info");
    infoAction->setToolTip("Volume Info");
    connect(infoAction, SIGNAL(triggered()), this, SLOT(restInfo()));

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

    this->menuBar()->addMenu(picker);
    this->menuBar()->addAction(clipAction);
    this->menuBar()->addAction(treeAction);
    this->menuBar()->addAction(infoAction);
    this->menuBar()->addAction(screenAction);
    this->menuBar()->addAction(screen4Action);

    rwidget = new RenderWidget(vd);
    this->setCentralWidget(rwidget);
    rwidget->setFocusPolicy(Qt::WheelFocus);

    // Clipping plane control
    dock_clip = new QDockWidget(this);
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
    vb->addStretch(1);
    cont->setLayout(vb);
    dock_clip->setWidget(cont);

    // Tree view (with vis/novis, hue control
    dock_tree = new QDockWidget(this);
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
    tree_view->header()->setStretchLastSection(false);
    tree_view->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    tree_view->setItemDelegateForColumn(2, new HueSpinBoxDelegate(tree_model));
    tree_view->setItemDelegateForColumn(3, new AlphaBoxDelegate(tree_model));
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
    dock_info = new QDockWidget(this);
    info_table = new QTableWidget();
    info_table->setColumnCount(1);
    QStringList keys;
    keys << "Name"
         << "Material"
         << "Density"
         << "Volume"
         << "Mass"
         << "Geocalls";
    info_table->setRowCount(keys.size());
    info_table->setVerticalHeaderLabels(keys);
    QStringList kv = QStringList() << "Value";
    info_table->setHorizontalHeaderLabels(kv);
    info_table->horizontalHeader()->setStretchLastSection(true);
    for (int i = 0; i < keys.size(); i++) {
        info_table->setItem(i, 0, new QTableWidgetItem(""));
    }

    dock_info->setWidget(info_table);

    addDockWidget(Qt::LeftDockWidgetArea, dock_clip);
    dock_clip->setAllowedAreas(Qt::AllDockWidgetAreas);
    dock_clip->setFeatures(QDockWidget::DockWidgetClosable |
                           QDockWidget::DockWidgetMovable);
    dock_clip->setVisible(false);

    addDockWidget(Qt::LeftDockWidgetArea, dock_tree);
    dock_tree->setAllowedAreas(Qt::AllDockWidgetAreas);
    dock_tree->setFeatures(QDockWidget::DockWidgetClosable |
                           QDockWidget::DockWidgetMovable);
    dock_tree->setVisible(false);

    addDockWidget(Qt::LeftDockWidgetArea, dock_info);
    dock_info->setAllowedAreas(Qt::AllDockWidgetAreas);
    dock_info->setFeatures(QDockWidget::DockWidgetClosable |
                           QDockWidget::DockWidgetMovable |
                           QDockWidget::DockWidgetFloatable);
    dock_info->setVisible(false);

    // set Layout, etc; use an autoshrinking list

    setMouseTracking(true);

    this->show();
}

Viewer::~Viewer() {}

void Viewer::restClip() { dock_clip->setVisible(true); }
void Viewer::restTree() { dock_tree->setVisible(true); }
void Viewer::restInfo() { dock_info->setVisible(true); }

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
                if (!geo_options[which_geo].cache) {
                    geo_options[which_geo].cache =
                        geo_options[which_geo].cons->Construct();
                }
                vd.root = convertCreation(geo_options[which_geo].cache);
                vd.scene_radius = vd.root.solid->GetExtent().GetExtentRadius();
                if (4 * vd.scene_radius > vd.camera.mag()) {
                    vd.camera *= 4 * vd.scene_radius / vd.camera.mag();
                }
                tree_model->recalculate();
                tree_view->expandToDepth(2);
                indicateElement(NULL);
                rwidget->rerender();
                return;
            }
        }
    }
}

void Viewer::indicateElement(Element *e) {
    if (!e) {
        for (int i = 0; i < info_table->rowCount(); i++) {
            info_table->item(i, 0)->setText("");
        }
    } else {
        // fill info table
        info_table->item(0, 0)->setText(e->name.c_str());
        info_table->item(1, 0)->setText(e->mat->GetName().c_str());
        G4double dens = e->mat->GetDensity();
        info_table->item(2, 0)->setText(
            QString::number(dens / (CLHEP::g / CLHEP::cm3), 'g', 4) + " g/cm3");
        G4double volume = e->solid->GetCubicVolume();
        info_table->item(3, 0)->setText(
            QString::number(volume / CLHEP::cm3, 'g', 4) + " cm3");
        info_table->item(4, 0)->setText(
            QString::number(volume * dens / CLHEP::kg, 'g', 4) + " kg");
        info_table->item(5, 0)->setText(QString::number(e->ngeocalls));
    }
}

void Viewer::screenshot(int sx) {
    RenderWorker w;
    QProgressDialog d("Rendering image", "Cancel render", 0,
                      rwidget->height() * sx);
    d.setMinimumDuration(1000);
    ViewData sccopy = vd;
    sccopy.level_of_detail = -1;
    connect(&d, SIGNAL(canceled()), &w, SLOT(coAbort()));
    QImage im(rwidget->width() * sx, rwidget->height() * sx,
              QImage::Format_RGB32);
    bool worked = w.render(sccopy, &im, 0, 1, &d);
    d.close();
    if (!worked) {
        return;
    }
    QString n = QFileDialog::getSaveFileName(this, "Save screenshot",
                                             QDir::current().canonicalPath());
    if (n.isEmpty()) {
        // Didn't choose a file...
        return;
    }
    QImageWriter r(n);
    r.setDescription("Screenshot of GEANT4 scene");
    r.write(im);
}
