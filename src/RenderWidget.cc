#include "RenderWidget.hh"

#include <QPainter>
#include <QThread>

RenderWidget::RenderWidget(ViewData &v, const TrackData &tdr)
    : QWidget(), currView(v), trackdata(tdr) {
    setAttribute(Qt::WA_OpaquePaintEvent, true);

    back = QImage(50, 50, QImage::Format_RGB32);
    back.fill(QColor::fromHslF(0.3, 0.5, 0.7));

    t = QVector<QThread *>();
    w = QVector<RenderWorker *>();
    qRegisterMetaType<ViewData>("ViewData");
    qRegisterMetaType<TrackData *>("TrackData*");
    qRegisterMetaType<TrackData>("TrackData");
    qRegisterMetaType<QImage *>("QImage*");
    int itc = QThread::idealThreadCount();
    itc = itc < 0 ? 2 : itc;
#if 0
    itc = 1; // Debug override
#endif
    for (int i = 0; i < itc; i++) {
        QThread *tt = new QThread();
        RenderWorker *ww = new RenderWorker();
        ww->moveToThread(tt);
        connect(ww, SIGNAL(completed()), this, SLOT(completed()));
        connect(ww, SIGNAL(aborted()), this, SLOT(aborted()));
        tt->start();
        t.append(tt);
        w.append(ww);
    }
    state = NONE;
    last_level_of_detail = 10000;
}

RenderWidget::~RenderWidget() {}

void RenderWidget::rerender() {
    currView.level_of_detail = 2;
    rerender_priv();
}

void RenderWidget::rerender_priv() {
    if (currView.level_of_detail > last_level_of_detail && state != NONE) {
        // don't know which workers have already completed, so abort all,
        // and queue an abort-clearer on all that ensures the next run
        // isn't aborted
        for (RenderWorker *ww : w) {
            ww->abort_task = true;
            QMetaObject::invokeMethod(ww, "flushAbort", Qt::QueuedConnection);
        }
// current response count indicates completed workers
// the aborted ones should bring the total up
#if 0
        qDebug("abort: response count %d", response_count);
#endif
        state = ACTIVE_AND_QUEUED;
        last_level_of_detail = 10000;
        return;
    }
    if (state == ACTIVE || state == ACTIVE_AND_QUEUED) {
        state = ACTIVE_AND_QUEUED;
        return;
    }

    int scl = int(std::pow(3.0, std::max(currView.level_of_detail, 0)));
    if (currView.level_of_detail == 0) {
        request_time = QTime::currentTime();
        arrived = aReqd;
    }
    // Q: SharedData on image so delete'ing Renderwidget doesn't crash Worker?
    next =
        QImage(this->width() / scl, this->height() / scl, QImage::Format_RGB32);
    for (int i = 0; i < w.size(); i++) {
        QMetaObject::invokeMethod(
            w[i], "render", Qt::QueuedConnection, Q_ARG(ViewData, currView),
            Q_ARG(TrackData, trackdata), Q_ARG(QImage *, &next), Q_ARG(int, i),
            Q_ARG(int, w.size()));
    }
    response_count = 0;

    last_level_of_detail = currView.level_of_detail;
    if (currView.level_of_detail <= 0) {
        state = ACTIVE;
    } else {
        state = ACTIVE_AND_QUEUED;
        currView.level_of_detail--;
    }
}

void RenderWidget::completed() {
    if (response_count < w.size()) {
        response_count++;
    }
    if (response_count < w.size()) {
        return;
    }
    back = next;
    if (state == ACTIVE) {
        state = NONE;
    } else if (state == ACTIVE_AND_QUEUED) {
        state = NONE;
        rerender_priv();
    }
    if (arrived == aReqd && back.size() == this->size()) {
        arrived = aCompl;
    }
    this->repaint();
}
void RenderWidget::aborted() {
    if (response_count < w.size()) {
        response_count++;
    }
    if (response_count < w.size()) {
        return;
    }
    state = NONE;
    rerender_priv();
}

void RenderWidget::resizeEvent(QResizeEvent *) {
    currView.level_of_detail = 2;
    rerender_priv();
}

void RenderWidget::paintEvent(QPaintEvent *) {
    QPainter q(this);
    QRect sz = back.rect();
    QImage mvd;
    if (sz.height() == this->height() && sz.width() == this->width()) {
        mvd = back;
    } else if (this->width() * sz.height() >= this->height() * sz.width()) {
        mvd = back.scaledToHeight(this->height(), Qt::FastTransformation);
    } else {
        mvd = back.scaledToWidth(this->width(), Qt::FastTransformation);
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
