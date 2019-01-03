/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#include "RenderWorker.hh"

#include <QElapsedTimer>
#include <QWidget>

class RenderWidget : public QWidget {
    Q_OBJECT
public:
    RenderWidget(ViewData &v, const TrackData &t);
    virtual ~RenderWidget();

signals:
    void frameTime(qreal);
    void forwardKey(QKeyEvent *e);
    void forwardMouse(QMouseEvent *);
    void forwardWheel(QWheelEvent *);
    void forwardContextMenu(QContextMenuEvent *);
    void forwardResize(QResizeEvent *evt);

public slots:
    void setFullDetail(bool);

    void completed(qreal);
    void aborted();

    void rerender(int changes);

protected:
    virtual void resizeEvent(QResizeEvent *evt) override;
    virtual void paintEvent(QPaintEvent *evt) override;
    virtual void keyPressEvent(QKeyEvent *) override;
    virtual void keyReleaseEvent(QKeyEvent *) override;
    virtual void mousePressEvent(QMouseEvent *) override;
    virtual void mouseReleaseEvent(QMouseEvent *) override;
    virtual void mouseMoveEvent(QMouseEvent *) override;
    virtual void wheelEvent(QWheelEvent *) override;
    virtual void contextMenuEvent(QContextMenuEvent *) override;

private:
    void rerender_priv();
    void drawRuler(QPainter &);
    ViewData &currView;
    const TrackData &trackdata;
    enum { NONE, ACTIVE, ACTIVE_AND_QUEUED } state;
    int last_level_of_detail;
    bool to_full_detail;
    int immediate_lod;

    RenderGraph graph;
    int changed_inputs;

    QImage cached;
    QSharedPointer<QImage> back;
    GridSpec back_grid;
    qint64 back_request_time;
    QElapsedTimer back_request_timer;
    QSharedPointer<QImage> next;
    GridSpec next_grid;
    QElapsedTimer next_request_timer;

    enum { aReqd, aCompl, aThere } arrived;
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
    RenderGraph graph;
    QProgressDialog *progress;
};
