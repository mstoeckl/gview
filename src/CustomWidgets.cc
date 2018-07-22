/* SPDX-License-Identifier: GPL-3.0-only */
#include "CustomWidgets.hh"

#include <QApplication>
#include <QComboBox>
#include <QFocusEvent>
#include <QGraphicsLineItem>
#include <QGraphicsPathItem>
#include <QGraphicsRectItem>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGridLayout>
#include <QLineEdit>
#include <QMenu>
#include <QMenuBar>
#include <QPushButton>
#include <QSignalMapper>

DistanceSpinBox::DistanceSpinBox(QWidget *parent) : QDoubleSpinBox(parent) {
    setRange(-10 * CLHEP::m, 10 * CLHEP::m);
    // TODO: inherit QAbstractSpinBox instead, and establish proper stepping
    // (outward increases units; inward never decreases them; from 0 uses
    // default)
}
DistanceSpinBox::~DistanceSpinBox() {}

QValidator::State DistanceSpinBox::validate(QString &text, int &pos) const {
    if (text.size() < 3) {
        return QValidator::Intermediate;
    }
    if (!text.endsWith("cm")) {
        return QValidator::Invalid;
    }

    QString cpy = text.left(text.size() - 2);
    QValidator::State s = QDoubleSpinBox::validate(cpy, pos);
    text = cpy + "cm";
    return s;
}
double DistanceSpinBox::valueFromText(const QString &text) const {
    return QDoubleSpinBox::valueFromText(text.left(text.size() - 2)) *
           CLHEP::cm;
}
QString DistanceSpinBox::textFromValue(double val) const {
    return QDoubleSpinBox::textFromValue(val / CLHEP::cm) + "cm";
}
NormalAxisSpinBox::NormalAxisSpinBox(NormalSelector *parent, int i)
    : QAbstractSpinBox(parent), link(parent), index(i) {
    setSpecialValueText(i == 0 ? "x" : (i == 1 ? "y" : "z"));
    setButtonSymbols(QAbstractSpinBox::UpDownArrows);
    setCorrectionMode(QAbstractSpinBox::CorrectToNearestValue);
    setWrapping(false);
    setKeyboardTracking(true);
    setAlignment(Qt::AlignRight);
    setInputMethodHints(Qt::ImhFormattedNumbersOnly);

    QLineEdit *le = lineEdit();
    le->setClearButtonEnabled(false);
    le->setMaxLength(6);
    le->setInputMask("#9.999");
    le->setText("+0.000");

    connect(le, SIGNAL(textEdited(const QString &)), this,
            SLOT(handleUpdate()));
}
NormalAxisSpinBox::~NormalAxisSpinBox() {}

QValidator::State NormalAxisSpinBox::validate(QString &input, int &) const {
    if (input.size() != 6) {
        // as per input mask, length *should* always be 6
        return QValidator::Invalid;
    }
    if (input.contains(QChar(' '))) {
        return QValidator::Invalid;
    }
    if (!(input[0] == '+' || input[0] == '-')) {
        return QValidator::Invalid;
    }
    if (!(input[1] == '0' || input[1] == '1')) {
        return QValidator::Invalid;
    }

    bool lim = input[1] == '1';
    if (lim && input.right(3) != "000") {
        // abs larger than 1 -- requires fixup
        return QValidator::Intermediate;
    }
    return QValidator::Acceptable;
}
void NormalAxisSpinBox::fixup(QString &input) const {
    // Clean up from Intermediate state
    bool lim = input[1] == '1';
    if (lim) {
        input[3] = QChar('0');
        input[4] = QChar('0');
        input[5] = QChar('0');
    }
}
void NormalAxisSpinBox::stepBy(int steps) { link->stepBy(steps, index); }
QAbstractSpinBox::StepEnabled NormalAxisSpinBox::stepEnabled() const {
    double v = link->current[index];
    return (v >= 1. ? StepNone : StepUpEnabled) |
           (v <= -1. ? StepNone : StepDownEnabled);
}
void NormalAxisSpinBox::displayValue(double val) {
    QString r = QString("%1").arg(val, 6, 'f', 3, '+');
    lineEdit()->setText(r);
}
void NormalAxisSpinBox::handleUpdate() { link->handleUpdate(); }
double NormalAxisSpinBox::apparentValue(bool &ok) const {
    const QLineEdit *le = lineEdit();
    int pos = lineEdit()->cursorPosition();
    QString text = le->text();
    QValidator::State state = validate(text, pos);
    if (state == QValidator::Invalid) {
        ok = false;
        return 0.;
    }
    if (state == QValidator::Intermediate) {
        fixup(text);
    }
    ok = true;
    return text.toDouble();
}
NormalSelector::NormalSelector(QWidget *parent) : QWidget(parent) {
    qRegisterMetaType<G4ThreeVector>("G4ThreeVector");
    setFocusPolicy(Qt::StrongFocus);

    sx = new NormalAxisSpinBox(this, 0);
    sy = new NormalAxisSpinBox(this, 1);
    sz = new NormalAxisSpinBox(this, 2);
    current = G4ThreeVector();

    QHBoxLayout *layout = new QHBoxLayout();
    layout->addWidget(sx, 1);
    layout->addWidget(sy, 1);
    layout->addWidget(sz, 1);
    layout->setContentsMargins(0, 0, 0, 0);
    this->setLayout(layout);

    connect(qApp, SIGNAL(focusChanged(QWidget *, QWidget *)), this,
            SLOT(trackFocusChange(QWidget *, QWidget *)));
}
NormalSelector::~NormalSelector() {}

void NormalSelector::setValue(const G4ThreeVector &value) {
    if (value.mag2() <= 0.) {
        qWarning("NormalSelector expected a nonzero vector");
        return;
    }
    G4ThreeVector repl = value.unit();
    if (current != repl) {
        // Only emit if change identified
        current = repl;
        sx->displayValue(current.x());
        sy->displayValue(current.y());
        sz->displayValue(current.z());
        emit valueChanged(current);
    }
}
G4ThreeVector NormalSelector::value() const { return current; }
void NormalSelector::stepBy(int dir, int axis) {
    // We move by astep toward the given pole.
    // If we are within eps of the opposite pole, we flip
    const double astep = M_PI / 12;
    const double cos_astep = std::cos(astep);
    const double eps = 1e-4;

    bool x0 = std::abs(current.x()) < eps;
    bool y0 = std::abs(current.y()) < eps;
    bool z0 = std::abs(current.z()) < eps;

    const G4ThreeVector &uc(current.unit());
    if (axis == 0) {
        if (dir * uc.x() >= cos_astep || (y0 && z0)) {
            current = dir * CLHEP::HepXHat;
        } else {
            G4ThreeVector perp = current.cross(CLHEP::HepXHat);
            CLHEP::HepRotation rot(perp.unit(), dir * astep);
            current = (rot * current).unit();
        }
    } else if (axis == 1) {
        if (dir * uc.y() >= cos_astep || (x0 && z0)) {
            current = dir * CLHEP::HepYHat;
        } else {
            G4ThreeVector perp = current.cross(CLHEP::HepYHat);
            CLHEP::HepRotation rot(perp.unit(), dir * astep);
            current = (rot * current).unit();
        }
    } else if (axis == 2) {
        if (dir * uc.z() >= cos_astep || (x0 && y0)) {
            current = dir * CLHEP::HepZHat;
        } else {
            G4ThreeVector perp = current.cross(CLHEP::HepZHat);
            CLHEP::HepRotation rot(perp.unit(), dir * astep);
            current = (rot * current).unit();
        }
    }
    sx->displayValue(current.x());
    sy->displayValue(current.y());
    sz->displayValue(current.z());
    emit valueChanged(current);
}
void NormalSelector::handleUpdate() {
    bool xok = true, yok = true, zok = true;
    G4ThreeVector a(sx->apparentValue(xok), sy->apparentValue(yok),
                    sx->apparentValue(zok));
    if (xok && yok && zok && a.mag2() > 0.) {
        // We have a valid state, up to scale.
        current = a.unit();
        emit valueChanged(current);
    }
}
void NormalSelector::trackFocusChange(QWidget *ow, QWidget *nw) {
    const bool was_inside = ow == sx || ow == sy || ow == sz;
    const bool is_inside = nw == sx || nw == sy || nw == sz;
    if (!was_inside && is_inside) {
        // focus in
    } else if (was_inside && !is_inside) {
        // focus out: display normalized values
        sx->displayValue(current.x());
        sy->displayValue(current.y());
        sz->displayValue(current.z());
    }
}

static void register_helper(QAction *act, QSignalMapper *map, int index) {
    map->setMapping(act, index);
    QObject::connect(act, SIGNAL(triggered()), map, SLOT(map()));
}

PlaneEdit::PlaneEdit(const Plane &p) {
    n = new NormalSelector(this);
    n->setValue(p.normal);
    d = new DistanceSpinBox();
    d->setValue(p.offset);

    act = new QPushButton("ON");
    act->setCheckable(true);
    connect(act, SIGNAL(toggled(bool)), this, SLOT(setActive(bool)));

    QSignalMapper *sigmap_action = new QSignalMapper(this);
    connect(sigmap_action, SIGNAL(mapped(int)), this, SLOT(doAction(int)));
    QMenuBar *mbar = new QMenuBar();
    QMenu *menu = new QMenu(" ... ", this);
    register_helper(menu->addAction("Flip"), sigmap_action, 0);
    register_helper(menu->addAction("Local X<0"), sigmap_action, 1);
    register_helper(menu->addAction("Local X>0"), sigmap_action, 2);
    register_helper(menu->addAction("Local Y<0"), sigmap_action, 3);
    register_helper(menu->addAction("Local Y>0"), sigmap_action, 4);
    register_helper(menu->addAction("Local Z<0"), sigmap_action, 5);
    register_helper(menu->addAction("Local Z>0"), sigmap_action, 6);
    mbar->addMenu(menu);

    connect(n, SIGNAL(valueChanged(G4ThreeVector)), this, SIGNAL(updated()));
    connect(d, SIGNAL(valueChanged(double)), this, SIGNAL(updated()));
    connect(act, SIGNAL(toggled(bool)), this, SIGNAL(updated()));

    QVBoxLayout *sup = new QVBoxLayout();
    sup->addWidget(n);
    QHBoxLayout *row = new QHBoxLayout();
    row->addWidget(d, 0, Qt::AlignLeft);
    row->addWidget(mbar, 1, Qt::AlignHCenter);
    row->addWidget(act, 0, Qt::AlignRight);
    sup->addLayout(row);

    this->setLayout(sup);
    setActive(false);
}

PlaneEdit::~PlaneEdit() {}

void PlaneEdit::doAction(int choice) {
    if (choice == 0) {
        n->setValue(-n->value());
        emit updated();
    } else if (choice <= 6) {
        G4ThreeVector o[6] = {CLHEP::HepXHat, -CLHEP::HepXHat,
                              CLHEP::HepYHat, -CLHEP::HepYHat,
                              CLHEP::HepZHat, -CLHEP::HepZHat};
        bool changed = false;
        if (d->value() != 0) {
            d->setValue(0.);
            changed = true;
        }
        if (n->value() != o[choice - 1]) {
            n->setValue(o[choice - 1]);
            changed = true;
        }
        if (changed) {
            emit updated();
        }
    }
}

void PlaneEdit::setActive(bool active) {
    n->setEnabled(active);
    d->setEnabled(active);
    act->setText(active ? "OFF" : "ON");
}

Plane PlaneEdit::getPlane() {
    Plane p;
    p.normal = act->isChecked() ? n->value() : G4ThreeVector();
    p.offset = d->value();
    return p;
}

ExpoSpinBox::ExpoSpinBox() {}
ExpoSpinBox::~ExpoSpinBox() {}

QValidator::State ExpoSpinBox::validate(QString &text, int &) const {
    static const QRegularExpression regExp(tr("(\\d).(\\d)e(\\d)"));
    Q_ASSERT(regExp.isValid());
    const QRegularExpressionMatch match = regExp.match(text);
    if (match.isValid()) {
        return QValidator::Acceptable;
    }
    // note: invalid are wrong character classes; may want to cull those
    // TODO: evtly, optimize to work precisely & error correct
    return QValidator::Intermediate;
}

int ExpoSpinBox::valueFromText(const QString &text) const {
    static const QRegularExpression regExp(tr("(\\d).(\\d)e(\\d)"));
    Q_ASSERT(regExp.isValid());

    const QRegularExpressionMatch match = regExp.match(text);
    if (match.isValid()) {
        int a = match.captured(1).toInt();
        int b = match.captured(2).toInt();
        int c = match.captured(3).toInt();
        return c * 90 + (a - 1) * 10 + b;
    }
    return 0;
}

QString ExpoSpinBox::textFromValue(int i) const {
    return QString("%1.%2e%3").arg((i % 90) / 10 + 1).arg(i % 10).arg(i / 90);
}

double ExpoSpinBox::expFromInt(int r) {
    int e = r / 90;
    double s = (r % 90) / 10. + 1.;
    return s * std::pow(10., e);
}
int ExpoSpinBox::nearestIntFromExp(double d) {
    int base = int(std::floor(std::log10(d)));
    double mant = d / std::pow(10., base);
    int imt = int(std::round(mant * 10.)) - 10;
    return 90 * base + imt;
}

HistogrammicRangeSlider::HistogrammicRangeSlider(bool exponential) {
    isexp = exponential;
    if (exponential) {
        espinLow = new ExpoSpinBox();
        espinHigh = new ExpoSpinBox();
        espinLow->setSingleStep(10);
        espinHigh->setSingleStep(10);
    } else {
        lspinLow = new QDoubleSpinBox();
        lspinHigh = new QDoubleSpinBox();
        lspinLow->setDecimals(3);
        lspinHigh->setDecimals(3);
        lspinLow->setSingleStep(0.1);
        lspinHigh->setSingleStep(0.1);
    }
    gv = new QGraphicsView();
    gv->setBaseSize(QSize(3, 1));
    scene = new QGraphicsScene();
    scene->setSceneRect(QRectF(0., 0., 200., 100.0));
    gv->setScene(scene);
    gv->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    gv->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    gv->setFixedHeight(2 * espinLow->minimumSizeHint().height());
    gv->fitInView(scene->sceneRect(), Qt::IgnoreAspectRatio);
    gv->setMinimumSize(QSize(0, 0));
    gv->setBaseSize(QSize(0, 0));
    gv->setSizePolicy(
        QSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding));

    // TODO: construct and add a histogram...
    scene->addRect(QRectF(1., 1., 20., 50.), QPen(Qt::black));
    // TODO: create a custom range - cover element
    // which slides if you click/drag center, and resizes
    // if you get act close to end points, and limits motion
    // at X and Y coordinates. Must be _semitransparent_
    limitLeft = scene->addLine(QLineF(QPointF(30., -300.), QPointF(30., 300.)),
                               QPen(Qt::red));
    limitRight = scene->addLine(QLineF(QPointF(70., -300.), QPointF(70., 300.)),
                                QPen(Qt::blue));
    limitLeft->setFlag(QGraphicsItem::ItemIsMovable, true);
    limitLeft->setFlag(QGraphicsItem::ItemIsSelectable, true);
    limitRight->setFlag(QGraphicsItem::ItemIsMovable, true);
    limitRight->setFlag(QGraphicsItem::ItemIsSelectable, true);

    QHBoxLayout *hlayout = new QHBoxLayout();
    hlayout->addWidget(espinLow);
    hlayout->addWidget(espinHigh);

    QVBoxLayout *vlayout = new QVBoxLayout();
    vlayout->addWidget(gv);
    vlayout->addLayout(hlayout);

    this->setLayout(vlayout);

    if (isexp) {
        connect(espinLow, SIGNAL(valueChanged(int)), this,
                SLOT(handleUpdate()));
        connect(espinHigh, SIGNAL(valueChanged(int)), this,
                SLOT(handleUpdate()));
    } else {
        connect(lspinLow, SIGNAL(valueChanged(double)), this,
                SLOT(handleUpdate()));
        connect(lspinHigh, SIGNAL(valueChanged(double)), this,
                SLOT(handleUpdate()));
    }
}
HistogrammicRangeSlider::~HistogrammicRangeSlider() { delete scene; }
void HistogrammicRangeSlider::setRange(double left, double right) {
    double low = std::min(left, right);
    double high = std::max(left, right);
    if (isexp) {
        espinLow->setRange(ExpoSpinBox::nearestIntFromExp(low),
                           ExpoSpinBox::nearestIntFromExp(high));
        espinHigh->setRange(ExpoSpinBox::nearestIntFromExp(low),
                            ExpoSpinBox::nearestIntFromExp(high));
    } else {
        lspinLow->setRange(low, high);
        lspinHigh->setRange(low, high);
    }
}
void HistogrammicRangeSlider::setHistogram(const QVector<QPointF> &pts) {
    // Linear transform the points & display
    scene->clear();

    double low, high;
    range(low, high);
    double mh = 0;
    double my = 0;
    for (QPointF p : pts) {
        mh = std::max(mh, p.y());
        my = std::max(my, p.x());
    }
    if (isexp) {
        low = std::log(low);
        high = std::log(high);
    }

    QPainterPath pp;
    pp.moveTo(QPointF(0, 0));
    for (QPointF p : pts) {
        if (isexp && p.x() <= 0.) {
            continue;
        }
        qreal z = isexp ? std::log(p.x()) : p.x();
        qreal x = 0.05 + 0.9 * scene->width() * (z - low) / (high - low);
        qreal y = scene->height() - 0.9 * scene->height() * p.y() / mh;
        pp.lineTo(QPointF(x, y));
    }

    scene->addPath(pp, QPen(Qt::black));

    // Q: need to resize viewport to actually show?
    gv->setUpdatesEnabled(true);
    gv->update();
}
void HistogrammicRangeSlider::value(double &left, double &right) const {
    double a, b;
    if (isexp) {
        a = ExpoSpinBox::expFromInt(espinLow->value());
        b = ExpoSpinBox::expFromInt(espinHigh->value());
    } else {
        a = lspinLow->value();
        b = lspinHigh->value();
    }
    left = std::min(a, b);
    right = std::max(a, b);
}
void HistogrammicRangeSlider::range(double &left, double &right) const {
    double a, b;
    if (isexp) {
        a = ExpoSpinBox::expFromInt(espinLow->minimum());
        b = ExpoSpinBox::expFromInt(espinHigh->maximum());
    } else {
        a = lspinLow->minimum();
        b = lspinHigh->maximum();
    }
    left = std::min(a, b);
    right = std::max(a, b);
}
void HistogrammicRangeSlider::setValue(double left, double right) {
    if (isexp) {
        int a = ExpoSpinBox::nearestIntFromExp(left);
        int b = ExpoSpinBox::nearestIntFromExp(right);
        espinLow->setValue(std::min(a, b));
        espinHigh->setValue(std::max(a, b));
    } else {
        lspinLow->setValue(std::min(left, right));
        lspinHigh->setValue(std::max(left, right));
    }
}
void HistogrammicRangeSlider::setSuffix(const QString &suff) {
    if (isexp) {
        espinLow->setSuffix(suff);
        espinHigh->setSuffix(suff);
    } else {
        lspinLow->setSuffix(suff);
        lspinHigh->setSuffix(suff);
    }
}
void HistogrammicRangeSlider::handleUpdate() {
    if (!this->signalsBlocked()) {
        emit valueChanged();
    }
}
void HistogrammicRangeSlider::resizeEvent(QResizeEvent *) {
    gv->fitInView(scene->sceneRect(), Qt::IgnoreAspectRatio);
}
