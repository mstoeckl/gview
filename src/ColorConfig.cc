#include "ColorConfig.hh"

#include "Overview.hh"

#include <G4Material.hh>
#include <G4VSolid.hh>

#include <QCheckBox>
#include <QColorDialog>
#include <QComboBox>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QListView>
#include <QPushButton>
#include <QTableView>
#include <QVBoxLayout>

QIcon iconFromColor(QColor c) {
    QPixmap p(QSize(256, 256));
    p.fill(c);
    return QIcon(p);
}

ColorConfig::ColorConfig(ViewData &ivd,
                         const std::vector<const G4Material *> &mtl_list)
    : vd(ivd) {
    material_list = mtl_list;
    active_mode = ColorByMaterial;
    mode_chooser = new QComboBox();
    mode_chooser->addItem("Material");
    mode_chooser->addItem("Property");
    mode_chooser->addItem("Flow map");
    mode_chooser->setCurrentIndex(0);
    connect(mode_chooser, SIGNAL(currentIndexChanged(int)), this,
            SLOT(changeMode()));

    div_by_class = new QCheckBox("Split by material");
    div_by_class->setCheckState(vd.split_by_material ? Qt::Checked
                                                     : Qt::Unchecked);
    connect(div_by_class, SIGNAL(stateChanged(int)), this, SLOT(changeMode()));

    mtl_table = new QTableView();
    mtl_table->setSelectionBehavior(QAbstractItemView::SelectItems);
    mtl_table->setSelectionMode(QAbstractItemView::NoSelection);

    for (int i = 0; i < int(material_list.size()); i++) {
        mtl_color_table.push_back(
            QColor::fromHslF(rand() / float(RAND_MAX - 1), 1., 0.5));
    }
    MaterialModel *mmod = new MaterialModel(mtl_color_table, material_list);
    mtl_table->setModel(mmod);
    mtl_table->horizontalHeader()->setStretchLastSection(true);
    mtl_table->setItemDelegateForColumn(0, new HueSpinBoxDelegate(mmod));
    connect(mmod, SIGNAL(colorChange()), this, SIGNAL(colorChange()));

    prop_select = new QComboBox();
    prop_select->addItems(QStringList() << "Density"
                                        << "Log Volume"
                                        << "Bool Tree");
    connect(prop_select, SIGNAL(currentIndexChanged(int)), this,
            SIGNAL(colorChange()));
    prop_base = QColor::fromRgbF(1, 1, 1);
    prop_target = QColor::fromRgbF(0, 0, 1);
    prop_base_button = new QPushButton(iconFromColor(prop_base), "Base color");
    prop_target_button =
        new QPushButton(iconFromColor(prop_target), "Target color");
    connect(prop_base_button, SIGNAL(pressed()), this,
            SLOT(updatePropBaseColor()));
    connect(prop_target_button, SIGNAL(pressed()), this,
            SLOT(updatePropTargetColor()));

    flow_load = new QPushButton("Load flow map...");
    connect(flow_load, SIGNAL(pressed()), this, SLOT(loadFlowMap()));
    flow_label = new QLabel("No paths loaded.");

    /* simulateneously layout everything, modifying visibility as needed */
    superlayout = new QVBoxLayout();

    superlayout->addWidget(mtl_table, 1);

    superlayout->addWidget(prop_select);
    superlayout->addWidget(prop_base_button);
    superlayout->addWidget(prop_target_button);

    superlayout->addWidget(flow_load);
    superlayout->addWidget(flow_label);

    stretch_widget = new QWidget();
    superlayout->addWidget(stretch_widget, 1);

    superlayout->addStretch(0);
    changeMode();

    QVBoxLayout *lyt = new QVBoxLayout();
    lyt->addWidget(div_by_class, 0);
    lyt->addWidget(mode_chooser, 0);
    lyt->addLayout(superlayout, 10);
    setLayout(lyt);
}
ColorConfig::~ColorConfig() {}

void ColorConfig::updatePropBaseColor() {
    QColor d = QColorDialog::getColor(prop_base, this, "Select base color");
    if (d.isValid()) {
        prop_base_button->setIcon(iconFromColor(d));
        prop_base = d;
        emit colorChange();
    }
}

void ColorConfig::updatePropTargetColor() {
    QColor d = QColorDialog::getColor(prop_target, this, "Select target color");
    if (d.isValid()) {
        prop_target_button->setIcon(iconFromColor(d));
        prop_target = d;
        emit colorChange();
    }
}

void recsetMtlColors(Element &e, const std::map<const G4Material *, int> &m) {
    e.ccode = m.at(e.material);
    for (Element &d : e.children) {
        recsetMtlColors(d, m);
    }
}
inline qreal mix(qreal a, qreal b, qreal i) {
    i = std::min(1., std::max(0., i));
    return (1. - i) * a + i * b;
}

void recsetPropColors(Element &e, int m, std::vector<QColor> &c, QColor f,
                      QColor t) {
    e.ccode = c.size();
    double amp = 0;
    switch (m) {
    default:
    case 0:
        /* todo: make a two-pass system that has sets equal-density classes....
         */
        amp = 1 / 20. * e.material->GetDensity() / CLHEP::g * CLHEP::cm3;
        break;
    case 1:
        amp = 1 / 5. *
              std::log10(const_cast<G4VSolid *>(e.solid)->GetCubicVolume() /
                         CLHEP::cm3);
        break;
    case 2: {
        QSet<const G4VSolid *> roots;
        int td = 0, nb = 0;
        calculateBooleanProperties(e.solid, roots, td, nb);
        amp = nb / 13.;
    } break;
    }
    /* linear interplation; looks bad for lots of combinations,
     * but works well for classes white-to-shade, black-to-shade */
    c.push_back(QColor::fromRgbF(mix(f.redF(), t.redF(), amp),
                                 mix(f.greenF(), t.greenF(), amp),
                                 mix(f.blueF(), t.blueF(), amp)));

    for (Element &d : e.children) {
        recsetPropColors(d, m, c, f, t);
    }
}

void recsetFlowColors(Element &e) {
    e.ccode = 0;
    for (Element &d : e.children) {
        recsetFlowColors(d);
    }
}

void ColorConfig::reassignColors() {
    switch (active_mode) {
    case ColorByMaterial: {
        vd.color_table = mtl_color_table;
        std::map<const G4Material *, int> idxs;
        for (int i = 0; i < int(material_list.size()); i++) {
            idxs[material_list[i]] = i;
        }
        recsetMtlColors(vd.elements, idxs);
    } break;
    case ColorByProperty: {
        vd.color_table.clear();
        recsetPropColors(vd.elements, prop_select->currentIndex(),
                         vd.color_table, prop_base, prop_target);
    } break;
    case ColorFromFlowmap: {
        vd.color_table.clear();
        vd.color_table.push_back(QColor(Qt::white));

        // uses: only the flowdb; for each elem,
        // calc inflow given restr. table
        recsetFlowColors(vd.elements);
    } break;
    }
}

void ColorConfig::changeMode() {
    for (int i = superlayout->count() - 1; i >= 0; i--) {
        if (superlayout->itemAt(i)->widget()) {
            superlayout->itemAt(i)->widget()->setVisible(false);
        }
    }
    switch (mode_chooser->currentIndex()) {
    case 0: {
        active_mode = ColorByMaterial;

        mtl_table->setVisible(true);
    } break;
    case 1: {
        active_mode = ColorByProperty;

        prop_select->setVisible(true);
        prop_base_button->setVisible(true);
        prop_target_button->setVisible(true);
        stretch_widget->setVisible(true);
    } break;
    case 2: {
        active_mode = ColorFromFlowmap;

        flow_load->setVisible(true);
        flow_label->setVisible(true);
        stretch_widget->setVisible(true);
    } break;
    }
    emit colorChange();
}
void ColorConfig::loadFlowMap() {
    /* The really fun part: we note deposition, for one */
    QString pth = QFileDialog::getOpenFileName(NULL, "Load Flow Map");
    if (pth.isEmpty()) {
        return;
    }

    QFile tf(pth);
    tf.open(QFile::ReadOnly);

    flow_base_n = 0;
    flow_db.clear();
    flow_names.clear();
    while (!tf.atEnd()) {
        QString line = tf.readLine();
        QStringList segments = line.split(" ");
        if (segments.size() != 10) {
            continue;
        }
        FlowData f;
        flow_base_n = std::max(
            flow_base_n, segments[4].mid(1, segments[4].size() - 2).toLong());
        f.inflow_val = segments[1].toDouble();
        f.inflow_err = segments[3].toDouble();
        f.deposit_val = segments[5].toDouble();
        f.deposit_err = segments[7].toDouble();
        f.nsamples = segments[9].mid(1, segments[9].size() - 2).toLong();

        QStringList key = segments[0].mid(2, segments[0].size() - 3).split('>');
        flow_db[key] = f;
        for (QString s : key) {
            flow_names.insert(s);
        }
    }
    flow_label->setText(
        QString("Paths: %1; N %2").arg(flow_db.size()).arg(flow_base_n));

    qDebug("Loaded flow map with %d paths and %d names", int(flow_base_n),
           flow_names.size());
}
