/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#include "RenderWorker.hh"
#include <QMainWindow>

class VectorTracer;
class QPushButton;
class G4VPhysicalVolume;
class QRadioButton;
class QLabel;
class QComboBox;
class QLineEdit;

class ImageWidget : public QWidget {
    Q_OBJECT
public:
    ImageWidget();
    virtual ~ImageWidget();
    void setImage(QImage im);
    bool hasImage() { return !image.isNull(); }

private:
    virtual void paintEvent(QPaintEvent *evt);
    QImage image;
};

class VectorPreview : public QMainWindow {
    Q_OBJECT
public:
    VectorPreview(ViewData vd, TrackData td);
    virtual ~VectorPreview();
private slots:
    void updateSettings();
    void updateColors();
    void queueRender();
    void acceptImage(QImage, QString, int, bool);

private:
    VectorTracer *tracer;
    QThread *thread;

    QPushButton *button_render;
    QPushButton *button_reroll;
    QLineEdit *line_name;
    QComboBox *combo_resolution;

    ImageWidget *display;
};

class VectorWindow : public QMainWindow {
    Q_OBJECT
public:
    VectorWindow(const char *name, G4VPhysicalVolume *vol);
    virtual ~VectorWindow();
public slots:
    void handleImageUpdate(QImage, QString, int, bool);
private slots:
    void closeProgram();
    void reset();

private:
    /* UI */
    QPushButton *button_full;
    QPushButton *button_step;
    QPushButton *button_reset;
    QPushButton *button_reroll;
    QLabel *label_step;
    QComboBox *list_resolution;

    ImageWidget *image_grid;
    ImageWidget *image_edge;
    ImageWidget *image_crease;
    ImageWidget *image_gradient;

    VectorTracer *vtracer;
    QThread *thread;
};
