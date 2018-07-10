/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#include "General.hh"
#include "Viewer.hh"

class RenderPoint {
public:
    RenderPoint(QPointF spot, const RayPoint &rpt);
    RenderPoint();
    ~RenderPoint();

    RenderPoint(const RenderPoint &);
    RenderPoint &operator=(RenderPoint);

    QPointF coords;
    RayPoint ray;

    FColor ideal_color;
    int region_class;
    int subregion_class;
    bool show_point;

private:
    void swap(RenderPoint &);
};

enum class GradientType { gSolid, gLinear, gRadial };

typedef struct {
    int subclass_no;
    // zero alpha signifies invisible
    QVector<QVector<RenderPoint>> boundaries;

    // Color parameters hold for subregions
    GradientType gradient_type;
    // No gradient parameter
    QRgb solid_color;
    // Linear gradient parameters
    float linear_angle;
    float linear_start;
    float linear_stop;
    int linear_nsteps;
    QVector<QRgb> linear_colors;
} Subregion;

typedef struct {
    int class_no;
    // Bounds within the grid
    int xmin, xmax, ymin, ymax;
    // Boundary
    QVector<RenderPoint> exterior;
    QVector<QVector<RenderPoint>> interior;
    // Add clipping plane marks
    bool is_clipped_patch;
    // Exterior lines
    QRgb meanExteriorColor;
    QVector<QRgb> meanInteriorColors;
    // Subregion color details
    QVector<Subregion> subregions;
} Region;

enum class Steps { sGrid, sEdges, sCreases, sGradients, sDone };

class VectorTracer : public QObject {
    Q_OBJECT
public:
    VectorTracer(ViewData view_data, TrackData track_data,
                 const QString &target_file, bool transparency = true,
                 QObject *parent = NULL);
    virtual ~VectorTracer();
    QImage preview(const QSize &sz);

public slots:
    void renderFull();
    void renderStep();
    void reset(bool transparent, const QSize &grid_size,
               const QString &target_name);
    void recolor();

    void computeGrid();
    void computeEdges();
    void computeCreases();
    void computeGradients();
signals:
    void produceImagePhase(QImage, QString message, int nqueries, bool done);

private:
    /* Function */
    RenderPoint queryPoint(QPointF);
    RenderPoint getPoint(QPoint);
    void bracketEdge(const RenderPoint &initial_inside,
                     const RenderPoint &initial_outside,
                     RenderPoint *result_inside, RenderPoint *result_outside);
    void bracketCrease(const RenderPoint &initial_inside,
                       const RenderPoint &initial_outside,
                       RenderPoint *result_inside, RenderPoint *result_outside);
    FColor calculateInteriorColor(const RenderPoint &pt);
    FColor calculateBoundaryColor(const RenderPoint &inside,
                                  const RenderPoint &outside);
    int faildepth(const RenderPoint &a, const RenderPoint &b);
    inline bool typematch(const RenderPoint &a, const RenderPoint &b) {
        return faildepth(a, b) < 0;
    }
    int crease_depth(const RenderPoint &q0, const RenderPoint &q1,
                     double cos_alpha, double min_jump, bool *is_jump = NULL);
    Steps step_next;
    QString file_name;
    long nqueries;
    bool transparent_volumes;
    QVector<FColor> element_colors;

    ViewData view_data;
    TrackData track_data;
    QSize grid_size;
    RenderPoint *grid_points;
    int grid_nclasses;
    QVector<Region> region_list;
    QMap<QPoint, RenderPoint> edge_refinements;
    QMap<QPoint, bool> crease_edge_map;
};
