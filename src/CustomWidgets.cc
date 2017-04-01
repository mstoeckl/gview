#include "CustomWidgets.hh"

#include <QComboBox>
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

double ExpoSpinBox::expFromInt(int r) const {
    int e = r / 90;
    double s = (r % 90) / 10. + 1.;
    return s * std::pow(10., e);
}
int ExpoSpinBox::nearestIntFromExp(double d) const {
    int base = int(std::floor(std::log10(d)));
    double mant = d / std::pow(10., base);
    int imt = int(std::round(mant * 10.)) - 10;
    return 90 * base + imt;
}
