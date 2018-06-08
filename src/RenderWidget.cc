#include "RenderWidget.hh"

#include <QDir>
#include <QFileDialog>
#include <QImageWriter>
#include <QPainter>
#include <QProgressDialog>
#include <QResizeEvent>
#include <QThread>

const G4double GOAL_FRAME_TIME = 0.020;
const int DOWNSCALE_BASE = 2;
const int MAX_LOD = 20; /* maximum render pixel size in pixels */

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

static QImage fastIntegerScale(const QImage &base, size_t S) {
    if (S <= 1) {
        return base;
    }
    if (base.pixelFormat().bitsPerPixel() != 32) {
        qFatal("Improper input pixel size, %d bits != 32 bits",
               base.pixelFormat().bitsPerPixel());
    }

    QElapsedTimer timer;
    timer.start();

    const size_t W = base.width(), H = base.height();
    QImage scaled(W * S, H * S, base.format());
    size_t slinelength = scaled.bytesPerLine();
    size_t blinelength = base.bytesPerLine();

    // Resampling + row copy is fastest & branch free
    const uchar *__restrict__ baseData = base.constBits();
    uchar *__restrict__ scaledData = scaled.bits();
    static size_t map[16384];
    if (W * S > 16384) {
        qFatal("Upsampling into huge image.");
    }
    for (size_t j = 0; j < W; j++) {
        for (size_t k = 0; k < S; k++) {
            map[j * S + k] = j;
        }
    }
    const size_t *__restrict__ remap = map;
    for (size_t k = 0; k < H; k++) {
        const uint32_t *__restrict__ source_line =
            (const uint32_t *)&baseData[k * blinelength];
        uint32_t *__restrict__ fill_line =
            (uint32_t * __restrict__) & scaledData[S * k * slinelength];
        for (size_t j = 0; j < W * S; j++) {
            fill_line[j] = source_line[remap[j]];
        }
        for (size_t i = 1; i < S; i++) {
            uint32_t *__restrict__ copy_line =
                (uint32_t * __restrict__) &
                scaledData[(S * k + i) * slinelength];
            memcpy(copy_line, fill_line, slinelength);
        }
    }

    return scaled;
}

RenderWidget::RenderWidget(ViewData &v, const TrackData &tdr)
    : QWidget(), currView(v), trackdata(tdr), graph(threadCount()) {
    setAttribute(Qt::WA_OpaquePaintEvent, true);

    cached = QImage(50, 50, QImage::Format_RGB32);
    cached.fill(QColor::fromHslF(0.3, 0.5, 0.7));
    back_scale_factor = 1;
    back_request_timer = QElapsedTimer();

    state = NONE;
    to_full_detail = false;
    last_level_of_detail = 10000;
    immediate_lod = 16;
    changed_inputs = CHANGE_COLOR | CHANGE_GEO | CHANGE_TRACK;

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
        changed_inputs |= CHANGE_GEO;
        rerender_priv();
    }
}

void RenderWidget::rerender(int what_changed) {
    if (what_changed != CHANGE_COLOR) {
        // Pure color changes are fast; all other types are slow
        currView.level_of_detail = ilog(DOWNSCALE_BASE, immediate_lod);
        changed_inputs |= CHANGE_VIEWPORT;
    }
    changed_inputs |= what_changed;
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
        arrived = aReqd;
    }
    next = QSharedPointer<QImage>(new QImage(
        this->width() / scl, this->height() / scl, QImage::Format_RGB32));
    next_scale_factor = scl;
    next_request_timer = QElapsedTimer();
    next_request_timer.start();

    TrackData d = currView.tracks;
    currView.tracks = trackdata;
    graph.start(next, currView, changed_inputs);
    changed_inputs = 0;
    currView.tracks = d;

    last_level_of_detail = currView.level_of_detail;
    if (currView.level_of_detail <= (to_full_detail ? -1 : 0)) {
        state = ACTIVE;
    } else {
        state = ACTIVE_AND_QUEUED;
        changed_inputs |= CHANGE_VIEWPORT;
        currView.level_of_detail--;
    }
}

void RenderWidget::completed(qreal time_secs) {
    back = next;
    back_scale_factor = next_scale_factor;
    back_request_timer = next_request_timer;

    if (state == ACTIVE) {
        state = NONE;
    } else if (state == ACTIVE_AND_QUEUED) {
        state = NONE;
        rerender_priv();
    }
    if (arrived == aReqd && back->size() == this->size()) {
        arrived = aCompl;
    }

    cached = fastIntegerScale(*back, back_scale_factor);
    this->repaint();

    if (back_scale_factor == immediate_lod) {
        /* Only change immediate_lod via its result */
        qreal time_per_pixel = time_secs / (back->width() * back->height());
        qreal target_pixels = GOAL_FRAME_TIME / time_per_pixel;
        qreal target_scale = isqrt(width() * height() / target_pixels);
        /* smooth changes */
        immediate_lod = (std::max((int)target_scale, 1) + immediate_lod) / 2;
        immediate_lod = std::min(immediate_lod, MAX_LOD);
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
    changed_inputs |= CHANGE_VIEWPORT;
    rerender_priv();
}

void RenderWidget::paintEvent(QPaintEvent *) {
    if (this->height() <= 0 || this->width() <= 0) {
        return;
    }

    QElapsedTimer paintTimer;
    paintTimer.start();

    QPainter q(this);
    QPoint corner =
        this->rect().center() - QPoint(cached.width() / 2, cached.height() / 2);
    q.drawImage(corner, cached);
    if (arrived == aCompl) {
        arrived = aThere;
    }
    qint64 im = paintTimer.nsecsElapsed();
    drawRuler(q);

#if 0
    qDebug("img completed after %f ms; %f ms paint (%f img)",
           back_request_timer.nsecsElapsed() * 1e-6,
           paintTimer.nsecsElapsed() * 1e-6, im * 1e-6);
#else
    (void)im;
#endif
}

void RenderWidget::drawRuler(QPainter &q) {
    // draw ruler in bottom left corner
    int max_ruler_length = std::max(this->width() / 3, 50);
    // Max ruler length in real space nm
    double ds =
        (0.5 * currView.scale) * max_ruler_length / this->width() / CLHEP::nm;

    int s = std::floor(std::log10(ds));
    int unit = std::min(6, std::max(0, s / 3 + 2));
    const char *labels[7] = {"fm", "pm", "nm", "μm", "mm", "m", "km"};

    double fracpart = ds / std::pow(10., s);

    int ruler_length = (max_ruler_length / fracpart);

    int dh = std::max(10, this->height() / 40);

    q.setPen(QColor::fromRgbF(0., 0., 0., 0.7));
    q.drawLine(dh, this->height() - dh, dh + ruler_length, this->height() - dh);
    q.drawLine(dh, this->height() - dh, dh, this->height() - 2 * dh);
    q.drawLine(dh + ruler_length, this->height() - dh, dh + ruler_length,
               this->height() - 2 * dh);

    QTextOption topt;
    topt.setAlignment(Qt::AlignTop | Qt::AlignRight);
    QPoint anchor = QPoint(dh, this->height() - dh);
    q.drawText(QRect(anchor - QPoint(dh, 0), anchor + QPoint(0, dh)),
               QString("0"), topt);

    QTextOption ropt;
    ropt.setAlignment(Qt::AlignTop | Qt::AlignLeft);
    QPoint ranchor = QPoint(dh + ruler_length, this->height() - dh);
    q.drawText(
        QRect(ranchor - QPoint(dh, 0), ranchor + QPoint(max_ruler_length, dh)),
        QString("%1 ").arg(ipow(10, s - 3 * (s / 3))) + labels[unit], ropt);
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
    graph.start(target, viewdata, CHANGE_ONESHOT);
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
