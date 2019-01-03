/* SPDX-License-Identifier: GPL-3.0-only */
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
class NameSelector;
class HistogrammicRangeSlider;

typedef struct {
    G4String name;
    G4VPhysicalVolume *vol;
    G4String suffix;
} GeoOption;

class Viewer : public QMainWindow {
    Q_OBJECT
public:
    Viewer(const std::vector<GeoOption> &options,
           const std::vector<TrackData> &trackopts);
    virtual ~Viewer();

public slots:
    void processKey(QKeyEvent *e);
    void processMouse(QMouseEvent *);
    void processWheel(QWheelEvent *);
    void processContextMenu(QContextMenuEvent *);
    void processResize(QResizeEvent *);

    void showFrameTime(qreal);
    void restClip();
    void restTree();
    void restInfo();
    void restRay();
    void restColor();
    void restOrient();
    void updatePlanes();
    void updateTracks(bool plane_change = false);
    void updateColors();
    void updateNavigator();
    void updateGShader();
    void updateTShader();
    void updateVoxDens();
    void updateShowLines();
    void screenshot(int sx = 1);
    void vectorScreenshot();
    void vectorPreview();
    void changeGeometry(QAction *);
    void changeTracks(QAction *);
    void indicateElement(const Element *);
    void rayLookup();
    void openGeometry();
    void openTracks();
    void setViewRotation(int);
    void updatePivot();

private:
    void reloadChoiceMenus();
    void reloadLineTypeSelection();

    std::vector<GeoOption> geo_options;
    std::vector<TrackData> track_options;
    std::vector<TrackRestriction> track_res_bounds;
    std::vector<TrackRestriction> track_res_actual;
    size_t which_geo;
    size_t which_tracks;
    ViewData vd;
    TrackData trackdata;
    RenderWidget *rwidget;
    QAction *frame_time_display;
    QAction *screen_action;
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
    QSpinBox *nanc_lower;
    QSpinBox *nanc_upper;
    QTreeView *tree_view;
    OverView *tree_model;
    QTableView *info_table;
    InfoModel *info_model;
    QListWidget *ray_table;
    QVector<const Element *> ray_list;
    QCheckBox *mtl_showlines;
    QComboBox *navig_sel, *gshader_sel, *tshader_sel;
    ColorConfig *color_config;
    QDoubleSpinBox *vox_density;
    QComboBox *pivot_volume;

    QMenu *gpicker_menu, *tpicker_menu, *view_menu, *screenshot_menu;
    QLabel *linecount_label;
    QListWidget *line_type_selection;

    QPointF clickpt;
    QPointF lastpt;
    bool clicked;
    bool shift;
    long rayiter;
};

void recursiveNameAppend(G4VPhysicalVolume *vp, const char *suffix);
