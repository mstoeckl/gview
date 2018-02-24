#pragma once

#include <QMainWindow>

#include <G4String.hh>

class VectorTracer;
class QPushButton;
class G4VPhysicalVolume;
class QRadioButton;
class QLabel;
class QComboBox;

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
    VectorWindow(G4String name, G4VPhysicalVolume *vol);
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
    QRadioButton *choice_transparent;
    QRadioButton *choice_opaque;
    QComboBox *list_resolution;

    ImageWidget *image_grid;
    ImageWidget *image_edge;
    ImageWidget *image_crease;
    ImageWidget *image_gradient;

    VectorTracer *vtracer;
    QThread *thread;
};
