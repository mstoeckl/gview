#include "Overview.hh"

#include "G4Material.hh"
#include "RenderWorker.hh"

#include <QDoubleSpinBox>
#include <QItemSelection>
#include <QSignalMapper>

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
    QDoubleSpinBox *spinBox = static_cast<QDoubleSpinBox *>(editor);
    spinBox->blockSignals(true);
    double v = index.model()->data(index, Qt::EditRole).toDouble();
    spinBox->setValue(v);
    spinBox->setStyleSheet("background-color: " +
                           QColor::fromHslF(v, 0.5, 0.7).name());
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

MaterialModel::MaterialModel(ViewData &c) : currView(c) {
    // calculate on/off corresponding to current element tree ?
}

MaterialModel::~MaterialModel() {}

int MaterialModel::rowCount(const QModelIndex &) const {
    return int(currView.matinfo.size());
}
int MaterialModel::columnCount(const QModelIndex &) const { return 1; }
QVariant MaterialModel::data(const QModelIndex &index, int role) const {
    if (!index.isValid()) {
        return QVariant();
    }
    int r = index.row();
    const MaterialInfo &m = currView.matinfo[size_t(r)];

    // Assume only valid indices...
    if (role == Qt::DisplayRole) {
        if (index.column() == 0) {
            return QVariant(QString::number(m.hue, 'f', 2));
        }
    }
    if (role == Qt::BackgroundColorRole) {
        if (index.column() == 0) {
            return QVariant(QColor::fromHslF(m.hue, 0.5, 0.7));
        }
    }
    if (role == Qt::EditRole) {
        if (index.column() == 0) {
            return QVariant(m.hue);
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
            currView.matinfo[size_t(r)].hue = value.toDouble();
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
    return QVariant(currView.matinfo[size_t(section)].mtl->GetName().data());
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
    currView.matinfo[size_t(row)].hue = e->value();
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
            n.elem = &currView.elements;
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

        // Determine lexical order for kids
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
