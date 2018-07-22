/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#include "RenderWorker.hh"

#include <QSpinBox>

class QComboBox;
class QPushButton;
class QGraphicsScene;
class QGraphicsView;
class QGraphicsLineItem;
class QGraphicsPathItem;

/**
 * Distance spin box, from 0.001nm to 999.9m. Displays 4 significant digits.
 *
 * Valid: 1.000nm, 32.3nm, 443.2nm, 560.3nm, 50.00m
 * Intermediate: 1e-7m, 1.234e5mm
 * Invalid: 16, seven, 1meter
 */
class DistanceSpinBox : public QAbstractSpinBox {
    Q_OBJECT
public:
    DistanceSpinBox(QWidget *parent = nullptr);
    virtual ~DistanceSpinBox();

    QValidator::State validate(QString &text, int &) const override;
    virtual void fixup(QString &text) const override;
    virtual void stepBy(int steps) override;
    virtual StepEnabled stepEnabled() const override;

    void setValue(double val);
    double value() const;

signals:
    void valueChanged(double);

private slots:
    void handleUpdate(const QString &);

private:
    static long parseValue(const QString &text, bool &ok, int &which_unit);
    static QString formatValue(long val, int unit);
    double current;
    long internal;
    int last_unit;
};

class NormalSelector;
/**
 * Helper class for QWidget. Format is precisely [+-]d.ddd
 */
class NormalAxisSpinBox : public QAbstractSpinBox {
    Q_OBJECT
public:
    NormalAxisSpinBox(NormalSelector *parent, int i);
    virtual ~NormalAxisSpinBox();

    virtual QValidator::State validate(QString &input, int &pos) const override;
    virtual void fixup(QString &input) const override;
    virtual void stepBy(int steps) override;
    virtual StepEnabled stepEnabled() const override;

    void displayValue(double val);
    double apparentValue(bool &ok) const;

private slots:
    void handleUpdate();

private:
    NormalSelector *link;
    int index;
};

/**
 * An efficient normalized direction selection widget
 */
class NormalSelector : public QWidget {
    Q_OBJECT
public:
    NormalSelector(QWidget *parent = nullptr);
    virtual ~NormalSelector();

    void setValue(const G4ThreeVector &);
    G4ThreeVector value() const;
signals:
    void valueChanged(G4ThreeVector);

protected:
    friend class NormalAxisSpinBox;
    void stepBy(int step, int axis);
    void handleUpdate();

    G4ThreeVector current, trial;
private slots:
    void trackFocusChange(QWidget *, QWidget *);

private:
    NormalAxisSpinBox *sx, *sy, *sz;
};

/*
 * Clipping plane editor. Designed for compact distance-angle representation.
 */
class PlaneEdit : public QWidget {
    Q_OBJECT
public:
    PlaneEdit(const Plane &p);
    virtual ~PlaneEdit();
    Plane getPlane();
    void setPlane(const Plane &);
    bool isActive();
public slots:
    void setActive(bool active);
signals:
    void updated();
private slots:
    void doAction(int choice);

private:
    NormalSelector *n;
    DistanceSpinBox *d;
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
    static double expFromInt(int);
    static int nearestIntFromExp(double);
};

class HistogrammicRangeSlider : public QWidget {
    Q_OBJECT
public:
    HistogrammicRangeSlider(bool exponential);
    virtual ~HistogrammicRangeSlider();
    void setRange(double left, double right);
    void range(double &left, double &right) const;
    void setSuffix(const QString &suff);
    void setHistogram(const QVector<QPointF> &p);
    void value(double &left, double &right) const;
    void setValue(double left, double right);
    virtual void resizeEvent(QResizeEvent *event);
signals:
    void valueChanged();
private slots:
    void handleUpdate();

private:
    QGraphicsView *gv;
    QGraphicsScene *scene;
    QGraphicsLineItem *limitRight;
    QGraphicsLineItem *limitLeft;
    union {
        QDoubleSpinBox *lspinLow;
        ExpoSpinBox *espinLow;
    };
    union {
        QDoubleSpinBox *lspinHigh;
        ExpoSpinBox *espinHigh;
    };
    bool isexp;
};
