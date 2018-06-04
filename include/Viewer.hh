#pragma once

#include "RenderWorker.hh"

#include <QMainWindow>
#include <QMap>

#include <vector>

class G4VPhysicalVolume;
class G4VUserDetectorConstruction;
class G4Material;

class QTableWidget;
class QTreeView;
class QDoubleSpinBox;
class QComboBox;
class QPushButton;
class QListWidget;
class QSpinBox;
class QCheckBox;
class QLabel;
class QTableView;
class ColorConfig;

class PlaneEdit;
class OverView;
class ExpoSpinBox;
class RenderWidget;
class InfoModel;
class HistogrammicRangeSlider;

typedef struct {
    G4String name;
    G4VPhysicalVolume *vol;
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
    void restColor();
    void restOrient();
    void updatePlanes();
    void updateColors();
    void updateShowLines();
    void screenshot(int sx = 1);
    void vectorTScreenshot();
    void vectorOScreenshot();
    void vectorPreview();
    void changeGeometry(QAction *);
    void changeTracks(QAction *);
    void indicateElement(const Element *);
    void rayLookup();
    void openGeometry();
    void openTracks();
    void setViewRotation(int);

private:
    void reloadChoiceMenus();

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
    QDockWidget *dock_color;
    QDockWidget *dock_orient;
    PlaneEdit *plane_edit[3];
    HistogrammicRangeSlider *times_range;
    HistogrammicRangeSlider *energy_range;
    QSpinBox *count_lower;
    QSpinBox *count_upper;
    QTreeView *tree_view;
    OverView *tree_model;
    QTableView *info_table;
    InfoModel *info_model;
    QListWidget *ray_table;
    QVector<const Element *> ray_list;
    QCheckBox *mtl_showlines;
    ColorConfig *color_config;

    QMenu *gpicker_menu;
    QMenu *tpicker_menu;
    QLabel *linecount_label;
    QListWidget *line_type_selection;

    QPoint clickpt;
    QPoint lastpt;
    bool clicked;
    bool shift;
    int rayiter;
};
