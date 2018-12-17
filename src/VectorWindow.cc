/* SPDX-License-Identifier: GPL-3.0-only */
#include "VectorWindow.hh"

#include "VectorTrace.hh"

#include <QApplication>
#include <QBitmap>
#include <QButtonGroup>
#include <QComboBox>
#include <QFile>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QImage>
#include <QLabel>
#include <QLineEdit>
#include <QMenuBar>
#include <QPainter>
#include <QPixmap>
#include <QProcess>
#include <QPushButton>
#include <QRadioButton>
#include <QSet>
#include <QSettings>
#include <QStatusBar>
#include <QTextStream>
#include <QThread>
#include <QTime>
#include <QVBoxLayout>

#include <G4VSolid.hh>
#include <G4VisExtent.hh>

ImageWidget::ImageWidget() : QWidget(), image() {
    QSizePolicy exp(QSizePolicy::Expanding, QSizePolicy::Expanding);
    this->setSizePolicy(exp);
}
ImageWidget::~ImageWidget() {}
void ImageWidget::setImage(QImage im) {
    // duplicate image to avoid case where data
    // is changed underneath it
    image = im.copy();
    this->update();
}
void ImageWidget::paintEvent(QPaintEvent *) {
    if (this->height() <= 0 || this->width() <= 0 || image.isNull()) {
        QPainter q(this);
        q.fillRect(rect(), Qt::white);
        return;
    }
    int s = std::min(this->width() / image.width(),
                     this->height() / image.height());

    QPainter q(this);
    q.fillRect(rect(), Qt::darkGray);
    QImage mvd;
    if (s > 0) {
        mvd = image.scaled(image.width() * s, image.height() * s,
                           Qt::IgnoreAspectRatio, Qt::FastTransformation);

    } else {
        mvd = image.scaled(this->width(), this->height(), Qt::KeepAspectRatio,
                           Qt::SmoothTransformation);
    }
    QPoint sz(mvd.width(), mvd.height());
    QRect target(this->rect().center() - sz / 2,
                 this->rect().center() + (sz - sz / 2));
    q.fillRect(target, Qt::lightGray);
    q.drawImage(target, mvd);
}

VectorPreview::VectorPreview(ViewData vd, TrackData td) {
    display = new ImageWidget();
    display->setMinimumSize(QSize(300, 300));

    button_render = new QPushButton("Render");
    button_reroll = new QPushButton("Recolor");
    line_name = new QLineEdit();
    line_name->setPlaceholderText("<automatic.svg>");
    combo_resolution = new QComboBox();
    combo_resolution->addItem("25x25");
    combo_resolution->addItem("100x100");
    combo_resolution->addItem("400x400");
    combo_resolution->addItem("1000x1000");
    combo_resolution->addItem("2500x2500");
    combo_resolution->setCurrentIndex(1);

    connect(button_reroll, SIGNAL(pressed()), this, SLOT(updateColors()));
    connect(combo_resolution, SIGNAL(currentIndexChanged(int)), this,
            SLOT(updateSettings()));
    connect(line_name, SIGNAL(textEdited(const QString &)), this,
            SLOT(updateSettings()));

    tracer = new VectorTracer(vd, td, QString());
    thread = new QThread();

    tracer->moveToThread(thread);
    connect(tracer, SIGNAL(produceImagePhase(QImage, QString, int, bool)), this,
            SLOT(acceptImage(QImage, QString, int, bool)));
    connect(button_render, SIGNAL(pressed()), this, SLOT(queueRender()));
    thread->start();
    // todo: termination handling -- only when window closes?

    QHBoxLayout *layout_columns = new QHBoxLayout();

    QVBoxLayout *layout_control = new QVBoxLayout();
    layout_control->addWidget(button_render);
    layout_control->addStrut(10);
    layout_control->addWidget(button_reroll);
    layout_control->addWidget(line_name);
    layout_control->addWidget(combo_resolution);
    layout_control->addStretch(1);

    layout_columns->addLayout(layout_control, 0);
    layout_columns->addWidget(display, 1);

    QWidget *central_widget = new QWidget();
    central_widget->setLayout(layout_columns);
    this->setCentralWidget(central_widget);

    updateSettings();

    this->setWindowFlag(Qt::Tool);
}
VectorPreview::~VectorPreview() { delete tracer; }
void VectorPreview::queueRender() {
    updateSettings();
    button_render->setEnabled(false);
    QMetaObject::invokeMethod(tracer, "renderFull", Qt::QueuedConnection);
}
void VectorPreview::updateSettings() {
    // We use an empty string to denote an automatically named file
    QString name = line_name->text();
    if (name.size()) {
        if (!name.endsWith(".pdf") && !name.endsWith(".svg")) {
            name = name + ".svg";
        }
    }

    if (button_render->isEnabled()) {
        QSize szs[5] = {QSize(25, 25), QSize(100, 100), QSize(400, 400),
                        QSize(1000, 1000), QSize(2500, 2500)};
        int i = combo_resolution->currentIndex();
        tracer->reset(szs[i], name);
        display->setImage(tracer->preview(QSize(100, 100)));
    }
}
void VectorPreview::updateColors() {
    tracer->recolor();
    updateSettings();
}
void VectorPreview::acceptImage(QImage im, QString, int, bool done) {
    display->setImage(im);
    if (done) {
        button_render->setEnabled(true);
    }
}

VectorWindow::VectorWindow(const char *name, G4VPhysicalVolume *vol)
    : QMainWindow() {
    this->setWindowTitle(name);

    /* Initialize functional stuff */

    ViewData view_data;
    view_data.orig_vol = vol;
    view_data.elements.clear();
    convertCreation(view_data.elements, vol, "");
    view_data.scene_radius =
        view_data.elements[0].solid->GetExtent().GetExtentRadius();
    view_data.scale = view_data.scene_radius; // *0.05 for hxrd3
    view_data.camera = G4ThreeVector(-2 * view_data.scene_radius, 0, 0);
    //    view_data.orientation = CLHEP::HepRotationX(37.44 * CLHEP::deg);
    view_data.orientation = G4RotationMatrix();
    view_data.level_of_detail = 0;
    view_data.split_by_material = true;
    view_data.clipping_planes = std::vector<Plane>();
    Plane p;
    p.normal = G4ThreeVector(1, 0, 0);
    p.offset = 0.;
    //    view_data.clipping_planes.push_back(p);
    view_data.color_table.clear();
    TrackData td;

    qsrand(QTime::currentTime().msec());
    vtracer = new VectorTracer(view_data, td, QString(), this);
    connect(vtracer, SIGNAL(produceImagePhase(QImage, QString, int, bool)),
            SLOT(handleImageUpdate(QImage, QString, int, bool)));

    //    thread = new QThread(this);
    //    vtracer->moveToThread(thread);
    //    thread->start();

    /* Initialize UI */
    this->menuBar()->addAction("Exit", this, SLOT(closeProgram()));

    button_full = new QPushButton("Render Full");
    connect(button_full, SIGNAL(pressed()), vtracer, SLOT(renderFull()));
    button_step = new QPushButton("Render Step");
    connect(button_step, SIGNAL(pressed()), vtracer, SLOT(renderStep()));
    label_step = new QLabel("Step: ---");
    button_reset = new QPushButton("Reset");
    connect(button_reset, SIGNAL(pressed()), this, SLOT(reset()));
    button_reroll = new QPushButton("Reroll colors");
    connect(button_reroll, SIGNAL(pressed()), vtracer, SLOT(recolor()));

    list_resolution = new QComboBox();
    list_resolution->addItem("1000", QVariant(QSize(1000, 1000)));
    list_resolution->addItem("300", QVariant(QSize(300, 300)));
    list_resolution->addItem("100", QVariant(QSize(100, 100)));
    list_resolution->addItem("30", QVariant(QSize(30, 30)));
    list_resolution->addItem("10", QVariant(QSize(10, 10)));

    QSettings s("gview", "gview");
    int idx = s.value("resol_index", 2).toInt();
    list_resolution->setCurrentIndex(idx);

    image_grid = new ImageWidget();
    image_edge = new ImageWidget();
    image_crease = new ImageWidget();
    image_gradient = new ImageWidget();

    QHBoxLayout *layout_columns = new QHBoxLayout();

    QVBoxLayout *layout_control = new QVBoxLayout();
    layout_control->addWidget(button_full);
    layout_control->addWidget(button_step);
    layout_control->addWidget(label_step);
    layout_control->addWidget(button_reset);
    layout_control->addWidget(button_reroll);
    layout_control->addWidget(list_resolution);
    layout_control->addStretch(1);
    layout_columns->addLayout(layout_control, 0);

    QGridLayout *layout_steps = new QGridLayout();
    layout_steps->addWidget(image_grid, 0, 0);
    layout_steps->addWidget(image_edge, 0, 1);
    layout_steps->addWidget(image_crease, 1, 0);
    layout_steps->addWidget(image_gradient, 1, 1);
    layout_columns->addLayout(layout_steps, 1);

    QWidget *central_widget = new QWidget();
    central_widget->setLayout(layout_columns);
    this->setCentralWidget(central_widget);

    this->statusBar()->showMessage("Status:");

    this->show();

    reset();
}
VectorWindow::~VectorWindow() {}

void VectorWindow::closeProgram() {
    //    QMetaObject::invokeMethod(thread, SLOT(quit()), Qt::QueuedConnection);
    QApplication::quit();
}
void VectorWindow::handleImageUpdate(QImage img, QString s, int nqueries,
                                     bool done) {
    if (!image_grid->hasImage()) {
        image_grid->setImage(img);
    } else if (!image_edge->hasImage()) {
        image_edge->setImage(img);
    } else if (!image_crease->hasImage()) {
        image_crease->setImage(img);
    } else if (!image_gradient->hasImage()) {
        image_gradient->setImage(img);
    }
    label_step->setText(s);
    statusBar()->showMessage(QStringLiteral("Queries: %1").arg(nqueries));
    if (done) {
        button_full->setEnabled(false);
        button_step->setEnabled(false);
    }
}
void VectorWindow::reset() {
    QSize resol = list_resolution->currentData().toSize();
    vtracer->reset(resol, QString());
    image_grid->setImage(QImage());
    image_edge->setImage(QImage());
    image_crease->setImage(QImage());
    image_gradient->setImage(QImage());
    button_full->setEnabled(true);
    button_step->setEnabled(true);

    QSettings s("gview", "gview");
    s.setValue("resol_index", list_resolution->currentIndex());
}
