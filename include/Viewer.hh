#pragma once

#include "RenderWidget.hh"
#include "RenderWorker.hh"

#include <QMainWindow>

#include <vector>

class G4VPhysicalVolume;
class G4VUserDetectorConstruction;

class QTableWidget;
class QTreeView;
class QDoubleSpinBox;
class QComboBox;
class QPushButton;
class QListWidget;
class QSpinBox;

class PlaneEdit;
class OverView;
class ExpoSpinBox;

typedef struct {
    G4String name;
    G4VUserDetectorConstruction *cons;
    G4VPhysicalVolume *cache;
} GeoOption;

typedef struct {
    Range energy;
    Range time;
    IRange seqno;
} TrackRestriction;

class Viewer : public QMainWindow {
    Q_OBJECT
public:
    Viewer(const std::vector<GeoOption> &options,
           const std::vector<TrackData> &trackopts);
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
    void restRay();
    void updatePlanes();
    void screenshot(int sx = 1);
    void changeGeometry(QAction *);
    void changeTracks(QAction *);
    void indicateElement(Element *);

private:
    std::vector<GeoOption> geo_options;
    std::vector<TrackData> track_options;
    std::vector<TrackRestriction> track_res_bounds;
    std::vector<TrackRestriction> track_res_actual;
    size_t which_geo;
    size_t which_tracks;
    ViewData vd;
    TrackData trackdata;
    RenderWidget *rwidget;
    QDockWidget *dock_clip;
    QDockWidget *dock_tree;
    QDockWidget *dock_info;
    QDockWidget *dock_ray;
    PlaneEdit *plane_edit[3];
    QDoubleSpinBox *times_lower;
    QDoubleSpinBox *times_upper;
    ExpoSpinBox *energy_lower;
    ExpoSpinBox *energy_upper;
    QSpinBox *count_lower;
    QSpinBox *count_upper;
    QTreeView *tree_view;
    OverView *tree_model;
    QTableWidget *info_table;
    QListWidget *ray_table;

    QPoint clickpt;
    QPoint lastpt;
    bool clicked;
    bool shift;
    int rayiter;
};
