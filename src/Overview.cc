/* SPDX-License-Identifier: GPL-3.0-only */
#include "Overview.hh"

#include "BooleanTree.hh"
#include "RenderWorker.hh"

#include "G4Material.hh"
#include <G4BooleanSolid.hh>
#include <G4DisplacedSolid.hh>

#include <QDoubleSpinBox>
#include <QItemSelection>
#include <QSignalMapper>

InfoModel::InfoModel(QObject *parent) : QAbstractTableModel(parent) {
    QStringList keys;
    keys << "Name"
         << "Material"
         << "Density"
         << "Volume"
         << "Mass"
         << "Surface Area"
         << "Bool roots"
         << "Bool depth"
         << "Bool splits"
         << "Color"
         << "Center";
    opts = keys.toVector();
    QStringList tool;
    tool << "Element name (from physical volume)"
         << "Material associated with the logical volume"
         << "Density taken from material property"
         << "Estimated solid volume"
         << "Product of estimated volume and material density"
         << "Estimated solid surface area"
         << "Number of distinct leaves in the tree formed by boolean solid "
            "composition"
         << "Depth of the boolean solid tree (incl. displacements) for this "
            "element"
         << "Number of boolean operations (counting duplicates) in boolean "
            "solid tree"
         << "Current color with which the object is displayed"
         << "Global coordinates for the solid zero point";
    tooltips = tool.toVector();
    vals = QVector<QString>(opts.size());
    col = QColor(Qt::white);
    last = NULL;
}
InfoModel::~InfoModel() {}

void calculateBooleanProperties(const G4VSolid *sol,
                                QSet<const G4VSolid *> &roots, int &treedepth,
                                int &nbooleans, int depth) {
    const BooleanTree *bso = dynamic_cast<const BooleanTree *>(sol);
    if (bso) {
        sol = bso->GetOriginal();
    }

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

void InfoModel::setElement(const Element *e, const ViewData &vd) {
    if (!e) {
        vals = QVector<QString>(opts.size());
        col = QColor(Qt::white);
    } else {
        // fill info table
        vals[0] = e->name;
        const G4Material *mat = e->material;
        vals[1] = mat->GetName();
        G4double dens = mat->GetDensity();
        vals[2] =
            QString::number(dens / (CLHEP::g / CLHEP::cm3), 'g', 4) + " g/cm3";
        G4double volume = const_cast<G4VSolid *>(e->solid)->GetCubicVolume();
        vals[3] = QString::number(volume / CLHEP::cm3, 'g', 4) + " cm3";
        vals[4] = QString::number(volume * dens / CLHEP::kg, 'g', 4) + " kg";
        G4double surf = const_cast<G4VSolid *>(e->solid)->GetSurfaceArea();
        vals[5] = QString::number(surf / CLHEP::cm2, 'g', 4) + " cm2";

        QSet<const G4VSolid *> roots;
        int treedepth = 0;
        int nbooleans = 0;
        calculateBooleanProperties(e->solid, roots, treedepth, nbooleans);
        vals[6] = QString::number(roots.size());
        vals[7] = QString::number(treedepth);
        vals[8] = QString::number(nbooleans);
        vals[10] =
            QString::number(e->global_offset.x() / CLHEP::cm, 'g', 4) + ", " +
            QString::number(e->global_offset.y() / CLHEP::cm, 'g', 4) + ", " +
            QString::number(e->global_offset.z() / CLHEP::cm, 'g', 4) + " cm";
        col = QColor::fromRgb(vd.color_table[e->ccode].rgb());
    }
    last = e;
    QVector<int> roles;
    roles.push_back(Qt::DisplayRole);
    emit dataChanged(index(0, 0), index(0, opts.size() - 1), roles);
}
int InfoModel::rowCount(const QModelIndex &) const { return int(opts.size()); }
int InfoModel::columnCount(const QModelIndex &) const { return 1; }
QVariant InfoModel::data(const QModelIndex &index, int role) const {
    if (!index.isValid()) {
        return QVariant();
    }
    // Assume only valid indices...
    if (role == Qt::DisplayRole) {
        if (index.column() == 0) {
            if (index.row() != 9) {
                return QVariant(vals[index.row()]);
            } else {
                return QVariant("");
            }
        }
    }
    if (index.row() == 9 && index.column() == 0) {
        if (role == Qt::BackgroundColorRole) {
            return QVariant(col);
        }
    }

    return QVariant();
}
bool InfoModel::setData(const QModelIndex &, const QVariant &, int) {
    return false;
}

QVariant InfoModel::headerData(int section, Qt::Orientation orientation,
                               int role) const {
    if (orientation == Qt::Horizontal) {
        return QVariant();
    }
    if (role == Qt::DisplayRole) {
        return QVariant(opts[section]);
    }
    if (role == Qt::ToolTipRole) {
        return QVariant(tooltips[section]);
    }
    return QVariant();
}
Qt::ItemFlags InfoModel::flags(const QModelIndex &) const {
    return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}

HueSpinBoxDelegate::HueSpinBoxDelegate(MaterialModel *model, QObject *parent)
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
    // Insert extra parameters in unused fields
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
    double nhue = index.model()->data(index, Qt::EditRole).toDouble();
    f3 o_srgb = rainbow_nhue(nhue);
    const QColor &c = QColor::fromRgbF(o_srgb[0], o_srgb[1], o_srgb[2]);

    QDoubleSpinBox *spinBox = static_cast<QDoubleSpinBox *>(editor);
    spinBox->blockSignals(true);
    spinBox->setValue(std::round(nhue * 100.) / 100.);
    spinBox->setStyleSheet("background-color: " + c.name());
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

MaterialModel::MaterialModel(std::vector<VColor> &c,
                             std::vector<const G4Material *> &m,
                             QObject *parent)
    : QAbstractTableModel(parent), colors(c), materials(m) {}

MaterialModel::~MaterialModel() {}
void MaterialModel::recalculate() {
    beginResetModel();
    endResetModel();
}

int MaterialModel::rowCount(const QModelIndex &) const {
    return int(colors.size());
}
int MaterialModel::columnCount(const QModelIndex &) const { return 1; }
QVariant MaterialModel::data(const QModelIndex &index, int role) const {
    if (!index.isValid()) {
        return QVariant();
    }
    int r = index.row();

    f3 c_srgb(colors[size_t(r)].redF(), colors[size_t(r)].greenF(),
              colors[size_t(r)].blueF());
    qreal nhue = color_srgb_to_nhue(c_srgb);

    // Assume only valid indices...
    if (role == Qt::DisplayRole) {
        if (index.column() == 0) {
            return QVariant(QString::number(nhue, 'f', 2));
        }
    }
    if (role == Qt::BackgroundColorRole) {
        if (index.column() == 0) {
            f3 o_srgb = rainbow_nhue(nhue);
            return QVariant(QColor::fromRgbF(o_srgb[0], o_srgb[1], o_srgb[2]));
        }
    }
    if (role == Qt::EditRole) {
        if (index.column() == 0) {
            return QVariant(nhue);
        }
    }
    return QVariant();
}
bool MaterialModel::setData(const QModelIndex &index, const QVariant &value,
                            int role) {
    if (!index.isValid()) {
        return false;
    }
    int r = index.row();
    if (role == Qt::EditRole) {
        if (index.column() == 0) {
            f3 o_srgb = rainbow_nhue(value.toDouble());

            colors[size_t(r)] =
                VColor(QColor::fromRgbF(o_srgb[0], o_srgb[1], o_srgb[2]).rgb());
            QVector<int> roles;
            roles.push_back(Qt::DisplayRole);
            emit dataChanged(index, index, roles);
            emit colorChange();
            return true;
        }
    }
    return false;
}

QVariant MaterialModel::headerData(int section, Qt::Orientation orientation,
                                   int role) const {
    if (role != Qt::DisplayRole) {
        return QVariant();
    }
    if (orientation == Qt::Horizontal) {
        return QVariant("Hue");
    }
    const G4Material *m = materials[size_t(section)];
    if (m) {
        return QVariant(m->GetName().data());
    } else {
        return QVariant("(null)");
    }
}
Qt::ItemFlags MaterialModel::flags(const QModelIndex &) const {
    return Qt::ItemIsEditable | Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}
void MaterialModel::hueUpdate(QWidget *w) {
    QDoubleSpinBox *e = static_cast<QDoubleSpinBox *>(w);
    if (!e) {
        return;
    }
    int row = e->maximumHeight() - 1000;
    double nhue = e->value();
    f3 o_srgb = rainbow_nhue(nhue);
    colors[size_t(row)] =
        VColor(QColor::fromRgbF(o_srgb[0], o_srgb[1], o_srgb[2]).rgb());
    QModelIndex idx = index(row, 0);
    QVector<int> roles;
    roles.push_back(Qt::EditRole);
    emit dataChanged(idx, idx, roles);
    emit colorChange();
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
    double v = index.model()->data(index, Qt::EditRole).toDouble();
    spinBox->setValue(v);
    QColor c = QColor::fromRgbF((3. + v) / 4., (3. + v) / 4., (3. + v) / 4.);
    spinBox->setStyleSheet("background-color: " + c.name());
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

OverView::OverView(ViewData &c, QObject *parent)
    : QAbstractItemModel(parent), currView(c) {
    recalculate();
}

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
    if (index.column() == 2) {
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
            return QVariant(QString::number(n.elem->alpha, 'f', 2));
        }
    }
    if (role == Qt::CheckStateRole && index.column() == 1) {
        return QVariant(n.elem->visible ? Qt::Checked : Qt::Unchecked);
    }
    if (role == Qt::EditRole) {
        if (index.column() == 2) {
            return QVariant(n.elem->alpha);
        }
    }
    if (role == Qt::BackgroundColorRole) {
        if (index.column() == 2) {
            return QVariant(QColor::fromRgbF((3. + n.elem->alpha) / 4.,
                                             (3. + n.elem->alpha) / 4.,
                                             (3. + n.elem->alpha) / 4.));
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

QModelIndex OverView::indexFromElement(const Element *e) {
    for (int i = 0; i < link.size(); i++) {
        const Node &n = link[i];
        if (n.elem == e) {
            const Node &p = link[n.parent];
            // Locate position relative to parent
            for (int k = 0; k < p.lexi.size(); k++) {
                int sidx = p.sub[p.lexi[k]];
                if (sidx == i) {
                    return createIndex(k, 0, reinterpret_cast<void *>(long(i)));
                }
            }
        }
    }
    return QModelIndex();
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
            n.elem = &currView.elements[0];
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
            n.elem =
                &currView
                     .elements[link[n.parent].elem->children[size_t(realidx)]];
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

        // Determine lexical order for kids
        QStringList names;
        for (size_t j = 0; j < nkids; j++) {
            names << QString(
                currView.elements[n.elem->children[j]].name.c_str());
            n.sub.push_back(-1); // placeholder
        }
        names.sort(Qt::CaseInsensitive);
        // assuming no duplicate names... (todo: a true indirect sort..)
        for (int i = 0; i < int(nkids); i++) {
            for (size_t j = 0; j < nkids; j++) {
                if (names[i] ==
                    QString(
                        currView.elements[n.elem->children[j]].name.c_str())) {
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
