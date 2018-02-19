#pragma once

#include <QMainWindow>

#include "Viewer.hh"

class VectorTracer;

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

class VectorWindow : public QMainWindow {
    Q_OBJECT
public:
    VectorWindow(GeoOption option);
    virtual ~VectorWindow();
public slots:
    void handleImageUpdate(QImage, QString, int, bool);
private slots:
    void closeProgram();

private:
    /* UI */ QPushButton *button_full;
    QPushButton *button_step;
    QLabel *label_step;

    ImageWidget *image_grid;
    ImageWidget *image_edge;
    ImageWidget *image_crease;
    ImageWidget *image_gradient;

    VectorTracer *vtracer;
    QThread *thread;
};
