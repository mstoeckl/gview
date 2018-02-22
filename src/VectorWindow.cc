#include "VectorWindow.hh"

#include "VectorTrace.hh"

#include <QApplication>
#include <QBitmap>
#include <QButtonGroup>
#include <QFile>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QImage>
#include <QLabel>
#include <QMenuBar>
#include <QPainter>
#include <QPixmap>
#include <QProcess>
#include <QPushButton>
#include <QRadioButton>
#include <QSet>
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

VectorWindow::VectorWindow(G4String name, G4VPhysicalVolume *vol)
    : QMainWindow() {
    this->setWindowTitle(name.c_str());

    /* Initialize functional stuff */

    ViewData view_data;
    view_data.elements = convertCreation(vol);
    view_data.scene_radius =
        view_data.elements.solid->GetExtent().GetExtentRadius();
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
    QString target = "vector.svg";

    qsrand(QTime::currentTime().msec());
    bool transparent = false;
    vtracer = new VectorTracer(view_data, td, target, transparent, this);
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

    choice_opaque = new QRadioButton("Opaque");
    choice_transparent = new QRadioButton("Transparent");
    QButtonGroup *choice_group = new QButtonGroup(this);
    choice_group->addButton(choice_opaque);
    choice_group->addButton(choice_transparent);
    (transparent ? choice_transparent : choice_opaque)->setChecked(true);
    connect(choice_opaque, SIGNAL(toggled(bool)), this, SLOT(reset()));
    connect(choice_transparent, SIGNAL(toggled(bool)), this, SLOT(reset()));

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
    layout_control->addWidget(choice_opaque);
    layout_control->addWidget(choice_transparent);
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
    statusBar()->showMessage(QString("Queries: %1").arg(nqueries));
    if (done) {
        button_full->setEnabled(false);
        button_step->setEnabled(false);
    }
}
void VectorWindow::reset() {
    bool is_transp = choice_transparent->isChecked();
    vtracer->reset(is_transp);
    image_grid->setImage(QImage());
    image_edge->setImage(QImage());
    image_crease->setImage(QImage());
    image_gradient->setImage(QImage());
    button_full->setEnabled(true);
    button_step->setEnabled(true);
}
