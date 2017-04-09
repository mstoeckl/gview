#pragma once

#include "RenderWorker.hh"

#include <QTime>
#include <QWidget>

class RenderWidget : public QWidget {
    Q_OBJECT
public:
    RenderWidget(ViewData &v, const TrackData &t);
    virtual ~RenderWidget();

    virtual void resizeEvent(QResizeEvent *evt);
    virtual void paintEvent(QPaintEvent *evt);
public slots:
    void completed();
    void aborted();

    void rerender();

private:
    void rerender_priv();
    ViewData &currView;
    const TrackData &trackdata;
    enum { NONE, ACTIVE, ACTIVE_AND_QUEUED } state;
    int last_level_of_detail;

    QVector<QThread *> t;
    QVector<RenderWorker *> w;
    int response_count;

    QSharedPointer<QImage> back;
    QSharedPointer<QImage> next;
    enum { aReqd, aCompl, aThere } arrived;
    QTime request_time;
};

class RenderSaveObject : public QObject {
    Q_OBJECT
public:
    RenderSaveObject(ViewData &v, const TrackData &t, int w, int h);
    virtual ~RenderSaveObject();
    void start();
public slots:
    void abort();
    void aborted();
    void completed();

private:
    ViewData &viewdata;
    const TrackData &trackdata;
    QSharedPointer<QImage> target;
    QThread *thread;
    RenderWorker *worker;
    QProgressDialog *progress;
};
