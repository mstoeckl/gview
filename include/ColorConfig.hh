/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#include "RenderWorker.hh"

#include <QLineEdit>
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
class QListView;

typedef struct {
    double inflow_val;
    double inflow_err;
    double deposit_val;
    double deposit_err;
    long nsamples;
} FlowData;

class NameComp;
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
    void addElements(const QStringList &);

    void showPopup();
    void filterPopup(const QString &);
    void applyChoice();

private:
    static NameComp *nc;
    QList<QString> names;
    QLineEdit *search;
    QLineEdit *search_prime;
    QFrame *container;
    QListWidget *list;

    QPushButton *wipe;
    QListWidget *collected;
};

class ColorConfig : public QWidget {
    Q_OBJECT
public:
    typedef enum {
        ColorByMaterial,
        ColorByMPreset,
        ColorByProperty,
        ColorFromFlowmap
    } ColorMode;

    ColorConfig(ViewData &ivd, const std::vector<const G4Material *> &mtl_list);
    ~ColorConfig();
    void mergeMaterials(const std::vector<const G4Material *> &mtl_list);
public slots:
    int reassignColors();
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
    QCheckBox *force_opaque;
    QComboBox *mode_chooser;
    ColorMode active_mode;
    QVBoxLayout *superlayout;
    QWidget *stretch_widget;

    QTableView *mtl_table;
    MaterialModel *mtl_model;
    std::vector<const G4Material *> material_list;
    std::vector<VColor> mtl_color_table;

    QComboBox *prop_select;
    QColor prop_base;
    QColor prop_target;
    QPushButton *prop_base_button;
    QPushButton *prop_target_button;

    QPushButton *flow_load;
    QLabel *flow_label;
    QVector<QPair<QVector<short>, FlowData>> flow_db;
    QMap<QString, short> flow_names;
    NameSelector *flow_target;
    NameSelector *flow_skip;
    NameSelector *flow_require;
    long flow_base_n;
};
