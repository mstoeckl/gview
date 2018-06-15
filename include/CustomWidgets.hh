/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#include "RenderWorker.hh"

#include <QSpinBox>

class QComboBox;
class QDoubleSpinBox;
class QPushButton;
class QGraphicsScene;
class QGraphicsView;
class QGraphicsLineItem;
class QGraphicsPathItem;

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
