#include "CustomWidgets.hh"

#include <QComboBox>
#include <QGraphicsLineItem>
#include <QGraphicsPathItem>
#include <QGraphicsRectItem>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QGridLayout>
#include <QPushButton>

PlaneEdit::PlaneEdit(Plane p) {
    nx = new QDoubleSpinBox();
    ny = new QDoubleSpinBox();
    nz = new QDoubleSpinBox();
    d = new QDoubleSpinBox();
    nx->setRange(-10, 10);
    ny->setRange(-10, 10);
    nz->setRange(-10, 10);
    d->setRange(-100, 100);
    d->setDecimals(1);
    nx->setValue(p.normal.x());
    ny->setValue(p.normal.y());
    nz->setValue(p.normal.z());
    d->setValue(p.offset);

    act = new QPushButton("ON");
    act->setCheckable(true);
    connect(act, SIGNAL(toggled(bool)), this, SLOT(setActive(bool)));

    unit = new QComboBox();
    unit->addItem("um");
    unit->addItem("mm");
    unit->addItem("cm");
    unit->addItem("m");
    unit->setCurrentIndex(2);

    connect(nx, SIGNAL(valueChanged(double)), this, SIGNAL(updated()));
    connect(ny, SIGNAL(valueChanged(double)), this, SIGNAL(updated()));
    connect(nz, SIGNAL(valueChanged(double)), this, SIGNAL(updated()));
    connect(d, SIGNAL(valueChanged(double)), this, SIGNAL(updated()));
    connect(act, SIGNAL(toggled(bool)), this, SIGNAL(updated()));
    connect(unit, SIGNAL(currentIndexChanged(int)), this, SIGNAL(updated()));

    QGridLayout *grid = new QGridLayout();
    grid->addWidget(nx, 0, 0);
    grid->addWidget(ny, 0, 1);
    grid->addWidget(nz, 0, 2);
    grid->addWidget(d, 1, 0);
    grid->addWidget(unit, 1, 1);
    grid->addWidget(act, 1, 2);

    this->setLayout(grid);
    setActive(false);
}

PlaneEdit::~PlaneEdit() {}

void PlaneEdit::setActive(bool active) {
    nx->setEnabled(active);
    ny->setEnabled(active);
    nz->setEnabled(active);
    d->setEnabled(active);
    unit->setEnabled(active);
    act->setText(active ? "OFF" : "ON");
}

Plane PlaneEdit::getPlane() {
    Plane p;
    p.normal = act->isChecked()
                   ? G4ThreeVector(nx->value(), ny->value(), nz->value())
                   : G4ThreeVector();
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
