#pragma once

#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>

#include <QAbstractItemModel>
#include <QImage>
#include <QItemDelegate>
#include <QMainWindow>
#include <QTime>

#include <vector>

class G4VPhysicalVolume;
class G4VUserDetectorConstruction;
class G4Material;
class G4VSolid;

class QKeyEvent;
class QDoubleSpinBox;
class QPushButton;
class QComboBox;
class QTreeView;
class QItemSelection;
class QTableWidget;
class QStyleOptionViewItem;
class QProgressDialog;

class PlaneEdit;
class OverView;

typedef struct {
    G4String name;
    G4VUserDetectorConstruction *cons;
    G4VPhysicalVolume *cache;
} GeoOption;

typedef struct {
    G4ThreeVector normal;
    G4double offset;
    // Keep n*x >= o; Drop n*x < o
} Plane;

typedef struct Element_s {
    G4String name;
    // Note: replicas not yet available
    G4ThreeVector offset;
    G4RotationMatrix rot;

    G4Material *mat;
    G4VSolid *solid;
    bool rotated;

    // only three changeable fields
    bool visible;
    double hue;
    double alpha;
    // Alpha = 1.0 : opaque; alpha < 1.0, we do linear sequential color merging
    // (yes, resulting colors may be weird. Not as good as exponential
    // color influence falloff, but you can't have everything.
    // (the background color is white!)

    // statistics, frequently updated, nolock
    mutable long ngeocalls;
    // Caching for acceleration
    mutable int niter;
    mutable double abs_dist;

    std::vector<struct Element_s> children;
} Element;

typedef struct {
    Element root;
    G4ThreeVector camera;
    G4RotationMatrix orientation;
    std::vector<Plane> clipping_planes;
    G4double scale;
    G4double scene_radius;
    int level_of_detail;
} ViewData;

class RenderWorker : public QObject {
    Q_OBJECT
public:
    RenderWorker();
    ~RenderWorker();
    bool abort_task;
public slots:
    bool render(ViewData p, QImage *i, int slice, int nslices,
                QProgressDialog *d = NULL);
    void coAbort();
    void flushAbort();
signals:
    void completed();
    void aborted();
};

class RenderWidget : public QWidget {
    Q_OBJECT
public:
    RenderWidget(ViewData &v);
    virtual ~RenderWidget();

    virtual void resizeEvent(QResizeEvent *evt);
    virtual void paintEvent(QPaintEvent *evt);
public slots:
    void completed();
    void aborted();

    void rerender();

private:
    void rerender_priv();
    ViewData &currView;
    enum { NONE, ACTIVE, ACTIVE_AND_QUEUED } state;
    int last_level_of_detail;

    QVector<QThread *> t;
    QVector<RenderWorker *> w;
    int response_count;

    QImage back;
    QImage next;
    enum { aReqd, aCompl, aThere } arrived;
    QTime request_time;
};

class PlaneEdit : public QWidget {
    Q_OBJECT
public:
    PlaneEdit(Plane p);
    virtual ~PlaneEdit();
    Plane getPlane();
public slots:
    void setActive(bool active);
signals:
    void updated();

private:
    QDoubleSpinBox *nx;
    QDoubleSpinBox *ny;
    QDoubleSpinBox *nz;
    QDoubleSpinBox *d;
    QComboBox *unit;
    QPushButton *act;
};

class HueSpinBoxDelegate : public QItemDelegate {
    Q_OBJECT
public:
    HueSpinBoxDelegate(OverView *model, QObject *parent = 0);

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
    OverView(ViewData &c);
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
    ViewData &currView;
    typedef struct {
        Element *elem;
        QVector<int> eaddr;

        int parent;
        QVector<int> sub;
        QVector<int> lexi;
    } Node;
    QVector<Node> link;
};

class Viewer : public QMainWindow {
    Q_OBJECT
    // viewer, etc; need a render method that fills a QImage...
public:
    Viewer(std::vector<GeoOption> options, size_t idx);
    virtual ~Viewer();

    virtual void keyPressEvent(QKeyEvent *event);
    virtual void mousePressEvent(QMouseEvent *event);
    virtual void mouseReleaseEvent(QMouseEvent *event);
    virtual void mouseMoveEvent(QMouseEvent *event);
    virtual void wheelEvent(QWheelEvent *event);
public slots:
    void restClip();
    void restTree();
    void restInfo();
    void updatePlanes();
    void screenshot(int sx = 1);
    void changeGeometry(QAction *);
    void indicateElement(Element *);

private:
    std::vector<GeoOption> geo_options;
    size_t which_geo;
    ViewData vd;
    RenderWidget *rwidget;
    QDockWidget *dock_clip;
    QDockWidget *dock_tree;
    QDockWidget *dock_info;
    PlaneEdit *plane_edit[3];
    QTreeView *tree_view;
    OverView *tree_model;
    QTableWidget *info_table;

    QPoint clickpt;
    QPoint lastpt;
    bool clicked;
    bool shift;
};
