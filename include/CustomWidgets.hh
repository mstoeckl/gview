#pragma once

#include "RenderWorker.hh"

#include <QSpinBox>

class QComboBox;
class QDoubleSpinBox;
class QPushButton;

class PlaneEdit : public QWidget {
    Q_OBJECT
public:
    PlaneEdit(Plane p);
    virtual ~PlaneEdit();
    Plane getPlane();
public slots:
    void setActive(bool active);
signals:
    void updated();

private:
    QDoubleSpinBox *nx;
    QDoubleSpinBox *ny;
    QDoubleSpinBox *nz;
    QDoubleSpinBox *d;
    QComboBox *unit;
    QPushButton *act;
};

class ExpoSpinBox : public QSpinBox {
    Q_OBJECT
public:
    ExpoSpinBox();
    virtual ~ExpoSpinBox();
    QValidator::State validate(QString &text, int &pos) const;
    virtual int valueFromText(const QString &) const;
    virtual QString textFromValue(int) const;
    double expFromInt(int) const;
    int nearestIntFromExp(double) const;
};
