#include "RenderWidget.hh"

#include <QDir>
#include <QFileDialog>
#include <QImageWriter>
#include <QPainter>
#include <QProgressDialog>
#include <QResizeEvent>
#include <QRunnable>
#include <QThread>
#include <QThreadPool>

const G4double GOAL_FRAME_TIME = 0.030;
const int DOWNSCALE_BASE = 2;

static int threadCount() {
    int itc = QThread::idealThreadCount();
    itc = itc < 0 ? 2 : itc;
#if 0
    itc = 1; // Debug override
#endif
    return itc;
}

static int ipow(int b, int e) {
    int x = 1;
    while (e > 0) {
        x *= b;
        e--;
    }
    return x;
}
static int ilog(int b, int x) {
    int e = 0;
    while (x > 1) {
        x /= b;
        e++;
    }
    return e;
}

static int isqrt(int x) { return std::sqrt(x); }

RenderWidget::RenderWidget(ViewData &v, const TrackData &tdr)
    : QWidget(), currView(v), trackdata(tdr), graph(threadCount()) {
    setAttribute(Qt::WA_OpaquePaintEvent, true);

    back = QSharedPointer<QImage>(new QImage(50, 50, QImage::Format_RGB32));
    back->fill(QColor::fromHslF(0.3, 0.5, 0.7));

    state = NONE;
    to_full_detail = false;
    last_level_of_detail = 10000;
    immediate_lod = 16;

    connect(&graph, SIGNAL(done(qreal)), this, SLOT(completed(qreal)));
    connect(&graph, SIGNAL(aborted()), this, SLOT(aborted()));
}

RenderWidget::~RenderWidget() {}

void RenderWidget::setFullDetail(bool b) {
    if (b == to_full_detail) {
        return;
    }
    to_full_detail = b;
    if (state == NONE) {
        rerender_priv();
    }
}

void RenderWidget::rerender() {
    currView.level_of_detail = ilog(DOWNSCALE_BASE, immediate_lod);
    rerender_priv();
}

void RenderWidget::rerender_priv() {
    if (currView.level_of_detail > last_level_of_detail && state != NONE) {
        graph.abort();
        state = ACTIVE_AND_QUEUED;
        last_level_of_detail = 10000;
        return;
    }
    if (state == ACTIVE || state == ACTIVE_AND_QUEUED) {
        state = ACTIVE_AND_QUEUED;
        return;
    }

    int ideal = ipow(DOWNSCALE_BASE, std::max(currView.level_of_detail, 0));
    int super = ipow(DOWNSCALE_BASE, ilog(DOWNSCALE_BASE, immediate_lod));
    int scl = std::max(1, (immediate_lod * ideal) / super);

    if (currView.level_of_detail == (to_full_detail ? -1 : 0)) {
        request_time = QTime::currentTime();
        arrived = aReqd;
    }
    next = QSharedPointer<QImage>(new QImage(
        this->width() / scl, this->height() / scl, QImage::Format_RGB32));

    TrackData d = currView.tracks;
    currView.tracks = trackdata;
    graph.start(next, currView);
    currView.tracks = d;

    last_level_of_detail = currView.level_of_detail;
    if (currView.level_of_detail <= (to_full_detail ? -1 : 0)) {
        state = ACTIVE;
    } else {
        state = ACTIVE_AND_QUEUED;
        currView.level_of_detail--;
    }
}

void RenderWidget::completed(qreal time_secs) {
    back = next;
    if (state == ACTIVE) {
        state = NONE;
    } else if (state == ACTIVE_AND_QUEUED) {
        state = NONE;
        rerender_priv();
    }
    if (arrived == aReqd && back->size() == this->size()) {
        arrived = aCompl;
    }
    this->repaint();

    int scale = std::max(width() / back->width(), height() / back->height());
    if (scale == immediate_lod) {
        /* Only change immediate_lod via its result */
        qreal time_per_pixel = time_secs / (back->width() * back->height());
        qreal target_pixels = GOAL_FRAME_TIME / time_per_pixel;
        qreal target_scale = isqrt(width() * height() / target_pixels);
        /* smooth changes */
        immediate_lod = (std::max((int)target_scale, 1) + immediate_lod) / 2;
    }
}
void RenderWidget::aborted() {
    state = NONE;
    rerender_priv();
}

void RenderWidget::resizeEvent(QResizeEvent *evt) {
    if (evt->size().width() <= 0 || evt->size().height() <= 0) {
        // Don't bother rendering empty images
        return;
    }
    currView.level_of_detail = ilog(DOWNSCALE_BASE, immediate_lod);
    rerender_priv();
}

void RenderWidget::paintEvent(QPaintEvent *) {
    if (this->height() <= 0 || this->width() <= 0) {
        return;
    }
    QPainter q(this);
    QRect sz = back->rect();
    QImage mvd;
    if (sz.height() == this->height() && sz.width() == this->width()) {
        mvd = *back;
    } else if (this->width() * sz.height() >= this->height() * sz.width()) {
        mvd = back->scaledToHeight(this->height(), Qt::FastTransformation);
    } else {
        mvd = back->scaledToWidth(this->width(), Qt::FastTransformation);
    }
    q.drawImage(
        this->rect().center() - QPoint(mvd.width() / 2, mvd.height() / 2), mvd);
    if (arrived == aCompl) {
        arrived = aThere;
#if 0
        qDebug("img completed after %d ms",
               request_time.msecsTo(QTime::currentTime()));
#endif
    }
}

RenderSaveObject::RenderSaveObject(ViewData &v, const TrackData &t, int w,
                                   int h)
    : viewdata(v), trackdata(t), graph(QThread::idealThreadCount() / 4 + 1) {
    target = QSharedPointer<QImage>(new QImage(w, h, QImage::Format_ARGB32));
    progress = new QProgressDialog("Rendering image", "Cancel render", 0, h);
    progress->setMinimumDuration(1000);
    progress->setRange(0, 100);
    connect(&graph, SIGNAL(aborted()), this, SLOT(aborted()));
    connect(&graph, SIGNAL(done(qreal)), this, SLOT(completed()));
    connect(&graph, SIGNAL(progressed(int)), progress, SLOT(setValue(int)));

    connect(progress, SIGNAL(canceled()), this, SLOT(abort()));
}

RenderSaveObject::~RenderSaveObject() { progress->deleteLater(); }

void RenderSaveObject::start() {
    // Change happens in main thread, so nothing bad inbetween
    int lod = viewdata.level_of_detail;
    viewdata.level_of_detail = -1;
    // Q: are copy costs too high? make LOD render side effect?
    TrackData d = viewdata.tracks;
    viewdata.tracks = trackdata;
    graph.start(target, viewdata);
    viewdata.tracks = d;

    viewdata.level_of_detail = lod;
}

void RenderSaveObject::abort() {
    graph.abort();
    progress->setVisible(false);
}

void RenderSaveObject::aborted() { this->deleteLater(); }

void RenderSaveObject::completed() {
    progress->setVisible(false);
    QString n = QFileDialog::getSaveFileName(NULL, "Save screenshot",
                                             QDir::current().canonicalPath());
    if (!n.isEmpty()) {
        // Chose a file
        QImageWriter r(n);
        r.setDescription("Screenshot of GEANT4 scene");
        r.write(*target);
    }
    this->deleteLater();
}
