#pragma once

#include "Viewer.hh"

class RenderPoint {
public:
    RenderPoint(QPointF spot, int nhits, const Intersection *srcints,
                const Element **srcelems);
    RenderPoint();
    ~RenderPoint();

    RenderPoint(const RenderPoint &);
    RenderPoint &operator=(RenderPoint);

    QPointF coords;
    int nhits;
    Intersection *intersections;
    Element **elements;
    QRgb ideal_color;
    int region_class;

private:
    void swap(RenderPoint &);
};

typedef struct {
    int class_no;
    // Bounds within the grid
    int xmin, xmax, ymin, ymax;
    // Boundary
    QVector<RenderPoint> exterior;
    QVector<QVector<RenderPoint>> interior;
    // Render cached properties
    bool is_clipped_patch;
    QRgb meanColor;
    QRgb meanExteriorColor;
    QVector<QRgb> meanInteriorColors;
} Region;

enum class Steps { sGrid, sEdges, sCreases, sGradients, sDone };

class VectorTracer : public QObject {
    Q_OBJECT
public:
    VectorTracer(ViewData view_data, TrackData track_data,
                 const QString &target_file, bool transparency = true,
                 QObject *parent = NULL);
    virtual ~VectorTracer();
public slots:
    void renderFull();
    void renderStep();

    void computeGrid();
    void computeEdges();
    void computeCreases();
    void computeGradients();
signals:
    void produceImagePhase(QImage, QString message, int nqueries, bool done);

private:
    /* Function */
    RenderPoint queryPoint(QPointF);
    void bracketEdge(const RenderPoint &initial_inside,
                     const RenderPoint &initial_outside,
                     RenderPoint *result_inside, RenderPoint *result_outside);
    QColor calculateInteriorColor(const RenderPoint &pt);
    QColor calculateBoundaryColor(const RenderPoint &inside,
                                  const RenderPoint &outside);
    int faildepth(const RenderPoint &a, const RenderPoint &b);
    inline bool typematch(const RenderPoint &a, const RenderPoint &b) {
        return faildepth(a, b) < 0;
    }
    Steps step_next;
    QString file_name;
    long nqueries;
    bool transparent_volumes;

    ViewData view_data;
    TrackData track_data;
    QSize grid_size;
    RenderPoint *grid_points;
    int grid_nclasses;
    QVector<Region> region_list;
    QMap<QPoint, RenderPoint> edge_refinements;
    ElemMutables *ray_mutables;
    long ray_iteration;
};
