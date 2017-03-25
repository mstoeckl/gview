#pragma once

#include "RenderWorker.hh"

#include <QImage>
#include <QMainWindow>
#include <QTime>

#include <vector>

class G4VPhysicalVolume;
class G4VUserDetectorConstruction;

class QTableWidget;
class QTreeView;
class QDoubleSpinBox;
class QComboBox;
class QPushButton;

class PlaneEdit;
class OverView;

typedef struct {
    G4String name;
    G4VUserDetectorConstruction *cons;
    G4VPhysicalVolume *cache;
} GeoOption;

class RenderWidget : public QWidget {
    Q_OBJECT
public:
    RenderWidget(ViewData &v, const TrackData &t);
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
    const TrackData &trackdata;
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
    TrackData trackdata;
    TrackData origtrackdata;
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
