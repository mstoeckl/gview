#pragma once

#include <QAbstractItemModel>
#include <QItemDelegate>

struct ViewData_s;
struct Element_s;
typedef struct Element_s Element;
class OverView;
class QItemSelection;

class HueSpinBoxDelegate : public QItemDelegate {
    Q_OBJECT
public:
    HueSpinBoxDelegate(OverView *model, QObject *parent = 0);
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
    OverView *oneTrueModel;
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
    OverView(struct ViewData_s &c);
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

    void recalculate();
signals:
    void colorChange();
    void selectedElement(Element *e);
public slots:
    void respToActive(const QModelIndex &index);
    void respToSelection(const QItemSelection &, const QItemSelection &);
    void hueUpdate(QWidget *);
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
