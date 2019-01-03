/* SPDX-License-Identifier: GPL-3.0-only */
#include "ColorConfig.hh"

#include "Overview.hh"

#include <G4Material.hh>
#include <G4VSolid.hh>

#include <QCheckBox>
#include <QCollator>
#include <QColorDialog>
#include <QComboBox>
#include <QCompleter>
#include <QFileDialog>
#include <QFocusEvent>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QListWidgetItem>
#include <QPushButton>
#include <QTableView>
#include <QVBoxLayout>

#include <zlib.h>

class NameComp {
private:
    QCollator q;

public:
    NameComp() { q.setNumericMode(true); }
    bool operator()(const QString &l, const QString &r) const {
        return q.compare(l, r) < 0;
    }
};
NameComp *NameSelector::nc = NULL;

NameSelector::NameSelector(QString label, QWidget *parent) : QWidget(parent) {
    if (!nc) {
        nc = new NameComp();
    }

    search = new QLineEdit();
    search_prime = new QLineEdit();
    container = new QFrame(this);
    list = new QListWidget();
    QVBoxLayout *hh = new QVBoxLayout(container);
    hh->addWidget(search_prime, 0);
    hh->addWidget(list, 10);
    hh->setContentsMargins(0, 0, 0, 0);
    container->setLayout(hh);
    connect(search, SIGNAL(textEdited(const QString &)), this,
            SLOT(showPopup()));
    connect(search_prime, SIGNAL(textEdited(const QString &)), this,
            SLOT(filterPopup(const QString &)));
    connect(search_prime, SIGNAL(returnPressed()), this, SLOT(applyChoice()));

    container->setWindowFlags(Qt::Popup);
    container->setFrameShadow(QFrame::Plain);
    container->setLineWidth(0);
    container->setMidLineWidth(0);
    list->setSelectionBehavior(QAbstractItemView::SelectRows);
    list->setSelectionMode(QAbstractItemView::ExtendedSelection);
    list->setFocusProxy(search_prime);
    list->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    list->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);

    wipe = new QPushButton("X");
    collected = new QListWidget();

    QHBoxLayout *hl = new QHBoxLayout();
    hl->addWidget(new QLabel(label), 0);
    hl->addWidget(search, 1);
    hl->addWidget(wipe, 0);
    connect(wipe, SIGNAL(pressed()), this, SLOT(clear()));

    QVBoxLayout *vl = new QVBoxLayout();
    vl->addLayout(hl);
    vl->addWidget(collected);
    setLayout(vl);
}
NameSelector::~NameSelector() {}

void NameSelector::showPopup() {
    filterPopup(search->text());
    QPoint corner = search->mapToGlobal(QPoint(0, 0));
    QSize sz(search->width(), search->height() * 8);
    QRect rct(corner, sz);
    container->setGeometry(rct);
    container->raise();
    container->show();
    search_prime->setFocus(Qt::OtherFocusReason);
    search_prime->setText(search->text());
}

void NameSelector::applyChoice() {
    QList<QListWidgetItem *> se = list->selectedItems();
    QStringList net;
    if (se.size()) {
        for (QListWidgetItem *s : se) {
            net.push_back(s->text());
        }
    } else {
        for (int i = 0; i < list->count(); i++) {
            net.push_back(list->item(i)->text());
        }
    }
    search_prime->setText("");
    search->setText("");
    container->hide();

    addElements(net);
}

void NameSelector::filterPopup(const QString &s) {
    search->setText(s);
    QRegExp r(s.split("").join(".*"));
    QStringList nn;
    for (QString q : names) {
        if (r.exactMatch(q)) {
            nn.append(q);
        }
    }
    list->clear();
    list->addItems(nn);
    QSet<QString> n;
    for (int i = 0; i < collected->count(); i++) {
        n.insert(collected->item(i)->text());
    }
    for (int i = 0; i < list->count(); i++) {
        QString w = list->item(i)->text();
        if (n.contains(w)) {
            list->item(i)->setTextColor(Qt::gray);
        }
    }
}

void NameSelector::setNames(const QSet<QString> &n) {
    QVector<QString> a = n.toList().toVector();
    qSort(a.begin(), a.end(), *nc);
    names = a.toList();
    filterPopup(search_prime->text());
    /* We don't reset the existing names, as unused
     * names don't harm anyone */
}
void NameSelector::clear() {
    if (collected->count()) {
        collected->clear();
        emit selectionChanged();
    }
}

void NameSelector::addElements(const QStringList &cands) {
    QSet<QString> sel = getSelected();
    int sz = sel.size();
    for (QString target : cands) {
        if (sel.contains(target)) {
            continue;
        }
        sel.insert(target);
    }
    if (sz == sel.size()) {
        return;
    }
    QVector<QString> s = sel.toList().toVector();
    qSort(s.begin(), s.end(), *nc);
    collected->clear();
    collected->addItems(s.toList());
    emit selectionChanged();
}

QSet<QString> NameSelector::getSelected() {
    QSet<QString> n;
    for (int i = 0; i < collected->count(); i++) {
        n.insert(collected->item(i)->text());
    }
    return n;
}

QIcon iconFromColor(QColor c) {
    QPixmap p(QSize(256, 256));
    p.fill(c);
    return QIcon(p);
}

ColorConfig::ColorConfig(ViewData &ivd,
                         const std::vector<const G4Material *> &mtl_list)
    : vd(ivd) {
    active_mode = ColorByMaterial;
    mode_chooser = new QComboBox();
    mode_chooser->addItem("Material", ColorByMaterial);
    mode_chooser->addItem("MPreset", ColorByMPreset);
    mode_chooser->addItem("Property", ColorByProperty);
    mode_chooser->addItem("Flow map", ColorByProperty);
    mode_chooser->setCurrentIndex(0);
    connect(mode_chooser, SIGNAL(currentIndexChanged(int)), this,
            SLOT(changeMode()));

    div_by_class = new QCheckBox("Split by type");
    div_by_class->setCheckState(vd.split_by_material ? Qt::Checked
                                                     : Qt::Unchecked);
    connect(div_by_class, SIGNAL(stateChanged(int)), this,
            SIGNAL(colorChange()));

    force_opaque = new QCheckBox("Force opaque");
    div_by_class->setCheckState(vd.force_opaque ? Qt::Checked : Qt::Unchecked);
    connect(force_opaque, SIGNAL(stateChanged(int)), this,
            SIGNAL(colorChange()));

    mtl_table = new QTableView();
    mtl_table->setSelectionBehavior(QAbstractItemView::SelectItems);
    mtl_table->setSelectionMode(QAbstractItemView::NoSelection);
    mtl_model = new MaterialModel(mtl_color_table, material_list, this);
    mergeMaterials(mtl_list);
    mtl_table->setModel(mtl_model);
    mtl_table->horizontalHeader()->setStretchLastSection(true);
    mtl_table->setItemDelegateForColumn(
        0, new HueSpinBoxDelegate(mtl_model, this));
    connect(mtl_model, SIGNAL(colorChange()), this, SIGNAL(colorChange()));

    prop_select = new QComboBox();
    prop_select->addItems(QStringList() << "Density"
                                        << "Log Volume"
                                        << "Bool Tree"
                                        << "Atomic Z"
                                        << "SA3/V2"
                                        << "Neighbors");
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
    flow_target = new NameSelector("Target");
    flow_skip = new NameSelector("Skip");
    flow_require = new NameSelector("Require");
    connect(flow_target, SIGNAL(selectionChanged()), this,
            SIGNAL(colorChange()));
    connect(flow_skip, SIGNAL(selectionChanged()), this, SIGNAL(colorChange()));
    connect(flow_require, SIGNAL(selectionChanged()), this,
            SIGNAL(colorChange()));
    flow_base_n = 0;

    /* simulateneously layout everything, modifying visibility as needed */
    superlayout = new QVBoxLayout();

    superlayout->addWidget(mtl_table, 1);

    superlayout->addWidget(prop_select);
    superlayout->addWidget(prop_base_button);
    superlayout->addWidget(prop_target_button);

    superlayout->addWidget(flow_load, 0);
    superlayout->addWidget(flow_label, 0);
    superlayout->addWidget(flow_target, 1);
    superlayout->addWidget(flow_skip, 1);
    superlayout->addWidget(flow_require, 1);

    stretch_widget = new QWidget();
    superlayout->addWidget(stretch_widget, 1);

    superlayout->addStretch(0);
    changeMode();

    QVBoxLayout *lyt = new QVBoxLayout();
    lyt->addWidget(div_by_class, 0);
    lyt->addWidget(force_opaque, 0);
    QHBoxLayout *mv = new QHBoxLayout();
    mv->addWidget(new QLabel("Mode:"), 0);
    mv->addWidget(mode_chooser, 1);
    lyt->addLayout(mv, 0);
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

static void recsetMtlColors(std::vector<Element> &elts,
                            const std::map<const G4Material *, int> &m, int i) {
    Element &e = elts[i];
    e.ccode = m.at(e.material);
    for (int j : e.children) {
        recsetMtlColors(elts, m, j);
    }
}
inline qreal mix(qreal a, qreal b, qreal i) {
    i = std::min(1., std::max(0., i));
    return (1. - i) * a + i * b;
}

class VCalc : public QRunnable {
public:
    VCalc(Element &e) : elem(e) {}
    void run() {
        elem.cubicVolume = const_cast<G4VSolid *>(elem.solid)->GetCubicVolume();
    }

private:
    Element &elem;
};
class SCalc : public QRunnable {
public:
    SCalc(Element &e) : elem(e) {}
    void run() {
        elem.surfaceArea = const_cast<G4VSolid *>(elem.solid)->GetSurfaceArea();
    }

private:
    Element &elem;
};

static void calcCachedProps(std::vector<Element> &elts, int m) {
    Nursery nursery;

    for (Element &e : elts) {
        if ((m == 1 || m == 4) && e.cubicVolume < 0) {
            nursery.startSoon(new VCalc(e), kDeleteAfter);
        }
        if (m == 4 && e.surfaceArea < 0) {
            nursery.startSoon(new SCalc(e), kDeleteAfter);
        }
    }
}

static void recgetProp(std::vector<Element> &elts, int e, int p, int m,
                       double *v) {
    double amp = 0;
    switch (m) {
    default:
    case 0:
        /* Tungsten has density 19.25 g/cm3, and Osmium wins at 22.59 g/cm3 */
        amp = elts[e].material->GetDensity() / CLHEP::g * CLHEP::cm3;
        break;
    case 1:
        amp = std::log10(elts[e].cubicVolume / CLHEP::cm3);
        break;
    case 2: {
        QSet<const G4VSolid *> roots;
        int td = 0, nb = 0;
        calculateBooleanProperties(elts[e].solid, roots, td, nb);
        amp = nb;
        break;
    }
    case 3: {
        /* Z is important w.r.t shielding */
        double net = 0;
        double nz = 0;
        for (size_t i = 0; i < elts[e].material->GetNumberOfElements(); i++) {
            G4Element *h = elts[e].material->GetElementVector()->at(i);
            double w = elts[e].material->GetFractionVector()[i] /
                       h->GetAtomicMassAmu();
            net += w;
            nz += h->GetZ() * w;
        }
        amp = (nz / net);
        break;
    }
    case 4: {
        /* SA3/V2 ratio as (horrible) roundness measure */
        double vol = elts[e].cubicVolume;
        double sa = elts[e].surfaceArea;
        double sav = sa * sa * sa / vol / vol;
        amp = std::log10(sav);
        break;
    }
    case 5: {
        /* Neighbors */
        amp = (elts[e].children.size() + elts[p].children.size());
        break;
    }
    }
    v[elts[e].ecode] = amp;
    for (int d : elts[e].children) {
        recgetProp(elts, d, e, m, v);
    }
}

static void recsetProp(std::vector<Element> &elts, int i, const int *m) {
    elts[i].ccode = m[elts[i].ecode];
    for (int j : elts[i].children) {
        recsetProp(elts, j, m);
    }
}

static void recsetMatchColors(std::vector<Element> &elts, int i,
                              const QMap<QString, short> &n) {
    elts[i].ccode = n.contains(elts[i].name.data()) ? 1 : 0;
    for (int j : elts[i].children) {
        recsetMatchColors(elts, j, n);
    }
}

class SortedStaticArraySet {
public:
    /* Very efficient for a small number of total elements
     * Alternatively, use sorted unique arrays */
    SortedStaticArraySet(const QVector<short> &in, short code_max) {
        sz = code_max;
        buf = new uint8_t[code_max];
        cnt = 0;
        look = 0;
        setSet(in);
    }
    SortedStaticArraySet(short code_max) {
        sz = code_max;
        buf = new uint8_t[code_max];
        cnt = 0;
        look = 0;
    }
    void setSet(const QVector<short> &in) {
        memset(buf, 0, sz * sizeof(uint8_t));
        for (short s : in) {
            buf[s] = 1;
        }
        cnt = 0;
        for (int i = 0; i < sz; i++) {
            cnt += buf[i];
        }
    }
    ~SortedStaticArraySet() {
        if (buf)
            delete[] buf;
        if (look)
            delete[] look;
    }
    bool contains(short e) const { return buf[e] > 0; }
    bool contains(const SortedStaticArraySet &e) const {
        if (look) {
            for (int i = 0; i < cnt; i++) {
                if (!e.buf[look[i]])
                    return false;
            }
        } else {
            for (int i = 0; i < sz; i++) {
                if (e.buf[i] && !buf[i]) {
                    return false;
                }
            }
        }
        return true;
    }
    bool intersects(const SortedStaticArraySet &e) const {
        if (look) {
            for (int i = 0; i < cnt; i++) {
                if (e.buf[look[i]])
                    return true;
            }
        } else {
            for (int i = 0; i < sz; i++) {
                if (e.buf[i] && buf[i]) {
                    return true;
                }
            }
        }
        return false;
    }
    int size() const { return cnt; }
    SortedStaticArraySet &fast() {
        if (look)
            return *this;
        look = new short[cnt];
        for (int i = 0, m = 0; i < sz; i++) {
            if (buf[i]) {
                look[m] = i;
                m++;
            }
        }
        return *this;
    }

private:
    uint8_t *buf;
    short *look;
    int sz;
    int cnt;
};

static void recsetFlowColors(std::vector<Element> &elts, int i,
                             const QMap<QString, short> &names,
                             const double *deps,
                             const SortedStaticArraySet &targets,
                             const SortedStaticArraySet &skip,
                             const SortedStaticArraySet &reqd, double total,
                             std::vector<VColor> &colors) {
    Element &e = elts[i];
    QString lstr(e.name.data());
    short label = -1;
    if (names.count(lstr)) {
        label = names[lstr];
    }
    if (label > 0 && skip.contains(label)) {
        e.ccode = 1;
    } else if (label > 0) {
        double tv = deps[label];
        /* span 5 orders of magnitude */
        double rv = 1.0 - (-std::log10(tv / total)) / 5.;
        QColor f = QColor(Qt::white);
        Qt::GlobalColor dt = Qt::green;
        if (targets.contains(label))
            dt = Qt::cyan;
        else if (reqd.contains(label))
            dt = Qt::yellow;
        QColor t = QColor(dt);
        e.ccode = colors.size();
        colors.push_back(VColor::fromRgbF(mix(f.redF(), t.redF(), rv),
                                          mix(f.greenF(), t.greenF(), rv),
                                          mix(f.blueF(), t.blueF(), rv)));
    } else {
        e.ccode = 0;
    }
    for (int j : e.children) {
        recsetFlowColors(elts, j, names, deps, targets, skip, reqd, total,
                         colors);
    }
}

static FColor randColor() {
    return FColor(randint(65536) / 65535., randint(65536) / 65535.,
                  randint(65536) / 65535., 1.0);
}

static std::vector<VColor>
preset_colors(const std::vector<const G4Material *> &material_list) {
    QMap<QString, VColor> color_map;
    color_map["ArGas"] = VColor::fromRgbF(0.8, 0.8, 0.8);
    color_map["Argas"] = color_map["ArGas"];
    color_map["ArgonGas"] = color_map["ArgonGas"];
    color_map["Air"] = VColor::fromRgbF(0.8, 0.8, 0.8);

    color_map["Al6061"] = VColor(255, 41, 58); // illu
    color_map["Al5083"] = VColor(159, 30, 45); // illu
    color_map["G4_Al"] = VColor(170, 0, 0);    // auto
    color_map["G4_C"] = color_map["G4_Al"];    // dummy

    // Often have nickel plated lead
    color_map["G4_Pb"] = VColor(146, 42, 153); // illu
    color_map["G4_Ni"] = VColor(180, 80, 187); // illu, twisted

    color_map["G4_Cu"] = VColor(98, 200, 156); // illu
    color_map["G4_Au"] = VColor::fromRgbF(1.0, 0.9, 0.0);
    color_map["Gold"] = color_map["G4_Au"];
    color_map["BaF2"] = VColor::fromRgbF(1.0, 1.0, 1.0);

    color_map["Teflon"] = VColor(94, 188, 236); // illu
    color_map["C2H2O"] = VColor::fromRgbF(0.9, 0.9, 0.5);
    color_map["Lexan"] = VColor(203, 250, 13); // illu, twisted
    color_map["Mylar"] = VColor::fromRgbF(0.7, 0.9, 0.7);

    color_map["PolyimideFoam"] = VColor(250, 178, 14); // illu
    color_map["Delrin"] = VColor::fromRgbF(0.6, 0.7, 0.9);
    color_map["UVSilica"] = VColor::fromRgbF(0.6, 0.8, 1.0);
    color_map["UHMW"] = VColor(250, 178, 14); // illu

    color_map["Polycarb"] = VColor::fromRgbF(0.7, 0.8, 0.3);
    color_map["Hevimet"] = VColor(38, 77, 164); // illu
    color_map["Kapton"] = VColor::fromRgbF(0.7, 0.7, 0.6);

    color_map["TungstenCarbide"] = VColor::fromRgbF(0.3, 0.3, 0.3);
    color_map["G4_W"] = VColor::fromRgbF(0.4, 0.3, 0.3);

    color_map["Stainless304"] = VColor::fromRgbF(0.6, 0.0, 0.0);
    color_map["Vespel"] = VColor::fromRgbF(0.7, 0.5, 0.3);

    std::vector<VColor> colors;
    for (const G4Material *m : material_list) {
        const char *cn = m->GetName().c_str();
        if (color_map.count(cn)) {
            colors.push_back(color_map[cn]);
        } else {
            qWarning("No preset color available for `%s`", cn);
            colors.push_back(VColor(randColor().rgba()));
        }
    }
    return colors;
}

int ColorConfig::reassignColors() {
    int render_change = CHANGE_COLOR;
    if (force_opaque->isChecked() != vd.force_opaque ||
        vd.split_by_material != div_by_class->isChecked()) {
        render_change |= CHANGE_GEO;
    }
    vd.force_opaque = force_opaque->isChecked();
    vd.split_by_material = div_by_class->isChecked();

    switch (active_mode) {
    case ColorByMPreset:
    case ColorByMaterial: {
        if (active_mode == ColorByMaterial) {
            vd.color_table = mtl_color_table;
        } else {
            vd.color_table.clear();
            vd.color_table = preset_colors(material_list);
        }
        std::map<const G4Material *, int> idxs;
        for (int i = 0; i < int(material_list.size()); i++) {
            idxs[material_list[i]] = i;
        }
        recsetMtlColors(vd.elements, idxs, 0);
    } break;
    case ColorByProperty: {
        int td = 0, ne = 0;
        countTree(vd.elements, 0, td, ne);
        double *vals = new double[ne]();
        int ci = prop_select->currentIndex();
        calcCachedProps(vd.elements, ci);
        recgetProp(vd.elements, 0, 0, ci, vals);
        double mn = kInfinity, mx = -kInfinity;
        for (int i = 0; i < ne; i++) {
            mn = std::min(vals[i], mn);
            mx = std::max(vals[i], mx);
        }
        if (mx - mn <= 0.) {
            mx = 1.;
            mn = 0.;
        }
        if (ci == 0 || ci == 2 || ci == 3 || ci == 5) {
            /* if not logarithmic, base at 0 */
            mn = 0.;
        }
        /* Q: establish classes for equal vals ? */
        vd.color_table.clear();

        int *refs = new int[ne]();
        QMap<double, int> vmatches;
        for (int i = 0; i < ne; i++) {
            if (vmatches.count(vals[i])) {
                refs[i] = vmatches[vals[i]];
                continue;
            }
            vmatches[vals[i]] = vd.color_table.size();
            refs[i] = vd.color_table.size();
            double v = (vals[i] - mn) / (mx - mn);
            /* linear interplation; looks bad for lots of combinations,
             * but works well for classes white-to-shade, black-to-shade */
            vd.color_table.push_back(VColor::fromRgbF(
                mix(prop_base.redF(), prop_target.redF(), v),
                mix(prop_base.greenF(), prop_target.greenF(), v),
                mix(prop_base.blueF(), prop_target.blueF(), v)));
        }
        recsetProp(vd.elements, 0, refs);
        delete[] vals;
        delete[] refs;
    } break;
    case ColorFromFlowmap: {
        /* TODO: single pass; construct cmap on load; then during pass,
         * estimate dose for everything in cmap; then, only in recursive
         * pass, div by dose if in table. Is O(nk) not O(nk^2). */
        QSet<short> ltargets, lskip, lreq;
        for (QString s : flow_target->getSelected())
            ltargets.insert(flow_names[s]);
        for (QString s : flow_skip->getSelected())
            lskip.insert(flow_names[s]);
        for (QString s : flow_require->getSelected())
            lreq.insert(flow_names[s]);
        SortedStaticArraySet targets(ltargets.toList().toVector(),
                                     flow_names.size() + 1);
        SortedStaticArraySet skip(lskip.toList().toVector(),
                                  flow_names.size() + 1);
        SortedStaticArraySet reqd(lreq.toList().toVector(),
                                  flow_names.size() + 1);
        targets.fast();
        skip.fast();
        reqd.fast();

        double total = 0.;
        double *deps = new double[flow_names.size() + 1]();
        SortedStaticArraySet as(flow_names.size() + 1);
        for (const QPair<QVector<short>, FlowData> &p : flow_db) {
            as.setSet(p.first);
            if ((!targets.size() || targets.contains(p.first.last())) &&
                !skip.intersects(as) && reqd.contains(as)) {
                for (short s : p.first) {
                    deps[s] += p.second.deposit_val;
                }
                total += p.second.deposit_val;
            }
        }
        if (total <= 0.) {
            vd.color_table.clear();
            vd.color_table.push_back(VColor::fromRgbF(0.4, 0.4, 0.4));
            vd.color_table.push_back(VColor::fromRgbF(0.5, 0.5, 1.));
            recsetMatchColors(vd.elements, 0, flow_names);
        } else {
            /* begin with unidentified color, then skip color */
            vd.color_table.clear();
            vd.color_table.push_back(VColor::fromRgbF(0.4, 0.4, 0.4));
            vd.color_table.push_back(VColor::fromRgbF(1.0, 0., 0.));
            recsetFlowColors(vd.elements, 0, flow_names, deps, targets, skip,
                             reqd, total, vd.color_table);
        }
        delete[] deps;
    } break;
    }

    return render_change;
}

void ColorConfig::changeMode() {
    for (int i = superlayout->count() - 1; i >= 0; i--) {
        if (superlayout->itemAt(i)->widget()) {
            superlayout->itemAt(i)->widget()->setVisible(false);
        }
    }
    active_mode = (ColorMode)mode_chooser->currentData().toInt();
    switch (active_mode) {
    case ColorByMaterial: {
        mtl_table->setVisible(true);
    } break;
    case ColorByMPreset:
        break;
    case ColorByProperty: {
        prop_select->setVisible(true);
        prop_base_button->setVisible(true);
        prop_target_button->setVisible(true);
        stretch_widget->setVisible(true);
    } break;
    case ColorFromFlowmap: {
        flow_load->setVisible(true);
        flow_label->setVisible(true);
        flow_target->setVisible(true);
        flow_skip->setVisible(true);
        flow_require->setVisible(true);
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

    /* Undo GZIP */
    system("rm -f /tmp/flow.gz");
    QString ar = "cp " + pth + " /tmp/flow.dat.gz";
    system(ar.toUtf8().constData());
    QString gu = "gzip -df /tmp/flow.dat.gz";
    system(gu.toUtf8().constData());
    pth = "/tmp/flow.dat";

    /* Actually load the file */
    QFile tf(pth);
    tf.open(QFile::ReadOnly);
    QByteArray header = tf.readLine();
    QByteArray names = tf.readLine();
    int j = 0;
    for (QByteArray b : names.mid(1, names.size() - 2).split(' ')) {
        /* Duplicated names grab the last code */
        QString s(b);
        flow_names[s] = j;
        j++;
    }
    flow_base_n = header.split('|')
                      .last()
                      .split(' ')
                      .last()
                      .replace(" ", "")
                      .replace("\n", "")
                      .toLong();

    QByteArray rest = tf.readAll();
    const char *rst = (const char *)rest.constData();
    const char *const orst = &rst[rest.size()];
    while (rst < orst) {
        QVector<short> tkey;
        for (const short *sp = (const short *)rst; *sp; ++sp) {
            tkey.push_back(*sp - 1);
        }
        rst += sizeof(short) * (tkey.size() + 1);

        FlowData f;
        const float *dp = (const float *)rst;
        // [inflow, deposit]x[electron, gamma, other]x[v,e]
        f.inflow_val = dp[0] + dp[2];
        f.deposit_val = dp[6] + dp[8];
        rst += sizeof(float) * 12;

        flow_db.append(QPair<QVector<short>, FlowData>(tkey, f));
    }
    flow_label->setText(
        QString("Paths: %1; N %2").arg(flow_db.size()).arg(flow_base_n));
    flow_require->setNames(flow_names.keys().toSet());
    flow_target->setNames(flow_names.keys().toSet());
    flow_skip->setNames(flow_names.keys().toSet());

    qDebug("Loaded flow map %s with %d samples and %d names",
           pth.toUtf8().constData(), int(flow_base_n), flow_names.size());
    emit colorChange();
}

void ColorConfig::mergeMaterials(
    const std::vector<const G4Material *> &mtl_list) {
    std::vector<const G4Material *> old_list = material_list;
    material_list = mtl_list;
    std::vector<VColor> old_colors = mtl_color_table;
    mtl_color_table.clear();

    for (size_t i = 0; i < material_list.size(); i++) {
        bool found = false;
        for (size_t j = 0; j < old_list.size(); j++) {
            if (old_list[j] == material_list[i]) {
                mtl_color_table.push_back(old_colors[j]);
                found = true;
                break;
            }
        }
        if (!found) {
            f3 color = rainbow_nhue(rand() / float(RAND_MAX - 1));
            mtl_color_table.push_back(
                VColor::fromRgbF(color[0], color[1], color[2]));
        }
    }
    mtl_model->recalculate();
}
