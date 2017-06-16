#pragma once

#include "RenderWorker.hh"

#include <QSet>
#include <QWidget>

class QCheckBox;
class QComboBox;
class QTableView;
class QVBoxLayout;
class QPushButton;
class QLabel;
class QListWidget;
class MaterialModel;

typedef struct {
    double inflow_val;
    double inflow_err;
    double deposit_val;
    double deposit_err;
    long nsamples;
} FlowData;

class NameSelector : public QWidget {
    Q_OBJECT
public:
    NameSelector(QString label, QWidget *parent = Q_NULLPTR);
    ~NameSelector();

    void setNames(const QSet<QString> &);
    QSet<QString> getSelected();
signals:
    void selectionChanged();
public slots:
    void clear();
private slots:
    void addElement(int);

private:
    QList<QString> names;
    QComboBox *search;
    QPushButton *wipe;
    QListWidget *collected;
};

class ColorConfig : public QWidget {
    Q_OBJECT
public:
    typedef enum {
        ColorByMaterial,
        ColorByProperty,
        ColorFromFlowmap
    } ColorMode;

    ColorConfig(ViewData &ivd, const std::vector<const G4Material *> &mtl_list);
    ~ColorConfig();
    void mergeMaterials(const std::vector<const G4Material *> &mtl_list);
public slots:
    void reassignColors();
private slots:
    void changeMode();
    void loadFlowMap();
    void updatePropBaseColor();
    void updatePropTargetColor();
signals:
    void colorChange();

private:
    ViewData &vd;

    QCheckBox *div_by_class;
    QComboBox *mode_chooser;
    ColorMode active_mode;
    QVBoxLayout *superlayout;
    QWidget *stretch_widget;

    QTableView *mtl_table;
    MaterialModel *mtl_model;
    std::vector<const G4Material *> material_list;
    std::vector<QColor> mtl_color_table;

    QComboBox *prop_select;
    QColor prop_base;
    QColor prop_target;
    QPushButton *prop_base_button;
    QPushButton *prop_target_button;

    QPushButton *flow_load;
    QLabel *flow_label;
    QVector<QPair<QStringList, FlowData>> flow_db;
    QSet<QString> flow_names;
    NameSelector *flow_target;
    NameSelector *flow_skip;
    NameSelector *flow_require;
    long flow_base_n;
};
