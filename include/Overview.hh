/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#include <QAbstractItemModel>
#include <QItemDelegate>

struct Element_s;
typedef struct Element_s Element;
typedef struct ViewData_s ViewData;
class OverView;
class MaterialModel;
class QItemSelection;
class G4Material;
class G4VSolid;
class VColor;

void calculateBooleanProperties(const G4VSolid *sol,
                                QSet<const G4VSolid *> &roots, int &treedepth,
                                int &nbooleans, int depth = 0);

class InfoModel : public QAbstractTableModel {
    Q_OBJECT
public:
    InfoModel(QObject *parent = 0);
    virtual ~InfoModel();
    void setElement(const Element *e, const ViewData &vd);
    virtual int rowCount(const QModelIndex &p = QModelIndex()) const;
    virtual int columnCount(const QModelIndex &p = QModelIndex()) const;
    virtual QVariant headerData(int section, Qt::Orientation orientation,
                                int role = Qt::DisplayRole) const;
    virtual QVariant data(const QModelIndex &index,
                          int role = Qt::DisplayRole) const;
    virtual bool setData(const QModelIndex &index, const QVariant &value,
                         int role = Qt::EditRole);
    Qt::ItemFlags flags(const QModelIndex &index) const;
    const Element *curE() const { return last; }

private:
    QVector<QString> opts;
    QVector<QString> tooltips;
    QVector<QString> vals;
    QColor col;
    const Element *last;
};

class HueSpinBoxDelegate : public QItemDelegate {
    Q_OBJECT
public:
    HueSpinBoxDelegate(MaterialModel *model, QObject *parent = 0);
    virtual ~HueSpinBoxDelegate() {}

    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option,
                          const QModelIndex &index) const;

    void setEditorData(QWidget *editor, const QModelIndex &index) const;
    void setModelData(QWidget *editor, QAbstractItemModel *model,
                      const QModelIndex &index) const;

    void updateEditorGeometry(QWidget *editor,
                              const QStyleOptionViewItem &option,
                              const QModelIndex &index) const;

private:
    MaterialModel *oneTrueModel;
};

class MaterialModel : public QAbstractTableModel {
    Q_OBJECT
public:
    MaterialModel(std::vector<VColor> &colors,
                  std::vector<const G4Material *> &materials,
                  QObject *parent = 0);
    virtual ~MaterialModel();

    virtual int rowCount(const QModelIndex &p = QModelIndex()) const;
    virtual int columnCount(const QModelIndex &p = QModelIndex()) const;
    virtual QVariant headerData(int section, Qt::Orientation orientation,
                                int role = Qt::DisplayRole) const;
    virtual QVariant data(const QModelIndex &index,
                          int role = Qt::DisplayRole) const;
    virtual bool setData(const QModelIndex &index, const QVariant &value,
                         int role = Qt::EditRole);
    Qt::ItemFlags flags(const QModelIndex &index) const;
    void recalculate();
signals:
    void colorChange();
public slots:
    void hueUpdate(QWidget *);

private:
    std::vector<VColor> &colors;
    std::vector<const G4Material *> &materials;
};

class AlphaBoxDelegate : public QItemDelegate {
    Q_OBJECT
public:
    AlphaBoxDelegate(OverView *model, QObject *parent = 0);
    virtual ~AlphaBoxDelegate() {}

    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option,
                          const QModelIndex &index) const;

    void setEditorData(QWidget *editor, const QModelIndex &index) const;
    void setModelData(QWidget *editor, QAbstractItemModel *model,
                      const QModelIndex &index) const;

    void updateEditorGeometry(QWidget *editor,
                              const QStyleOptionViewItem &option,
                              const QModelIndex &index) const;

private:
    OverView *oneTrueModel;
};

class OverView : public QAbstractItemModel {
    Q_OBJECT
public:
    OverView(struct ViewData_s &c, QObject *parent = 0);
    virtual ~OverView();

    virtual QModelIndex index(int r, int c,
                              const QModelIndex &p = QModelIndex()) const;
    virtual QModelIndex parent(const QModelIndex &chld) const;
    virtual Qt::ItemFlags flags(const QModelIndex &index) const;
    virtual int rowCount(const QModelIndex &p = QModelIndex()) const;
    virtual int columnCount(const QModelIndex &p = QModelIndex()) const;
    virtual QVariant headerData(int section, Qt::Orientation orientation,
                                int role = Qt::DisplayRole) const;
    virtual QVariant data(const QModelIndex &index,
                          int role = Qt::DisplayRole) const;
    virtual bool setData(const QModelIndex &index, const QVariant &value,
                         int role = Qt::EditRole);
    QModelIndex indexFromElement(const Element *e);

    void recalculate();
signals:
    void colorChange();
    void selectedElement(const Element *e);
public slots:
    void respToActive(const QModelIndex &index);
    void respToSelection(const QItemSelection &, const QItemSelection &);
    void alphaUpdate(QWidget *);

private:
    struct ViewData_s &currView;
    typedef struct {
        struct Element_s *elem;
        QVector<int> eaddr;

        int parent;
        QVector<int> sub;
        QVector<int> lexi;
    } Node;
    QVector<Node> link;
};
