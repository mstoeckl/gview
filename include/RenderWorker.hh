#pragma once

#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>

#include <QObject>
#include <QRgb>
#include <QSharedData>
#include <QSharedPointer>

#include <vector>

class LineCollection;
class G4Material;
class G4VSolid;
class QProgressDialog;

typedef struct ViewData_s ViewData;

typedef struct {
    G4ThreeVector normal;
    G4double offset;
    // Keep n*x >= o; Drop n*x < o
} Plane;

typedef struct {
    G4ThreeVector normal;
    G4double dist;
} Intersection;

typedef struct {
    int npts;
    int ptype;
    size_t offset;
    G4double bballradius;
} TrackHeader;

typedef struct {
    float energy;
    float time;
    double x;
    double y;
    double z;
} TrackPoint;

typedef struct { double low, high; } Range;
typedef struct { size_t low, high; } IRange;

class TrackPrivateData : public QSharedData {
public:
    size_t ntracks;
    size_t npoints;
    TrackHeader *headers;
    TrackPoint *points;
    LineCollection *tree;
    TrackPrivateData(size_t itracks, size_t ipoints);
    TrackPrivateData(const TrackPrivateData &other);
    ~TrackPrivateData();
};

class TrackData {
public:
    TrackData();
    TrackData(const char *filename);
    TrackData(const TrackData &other);
    TrackData(const TrackData &other, ViewData &viewrestr, Range seltimes,
              Range selenergies, IRange selidxs);
    ~TrackData();
    size_t getNPoints() const;
    size_t getNTracks() const;
    const TrackHeader *getHeaders() const;
    const TrackPoint *getPoints() const;
    const LineCollection *getTree() const;
    void calcTimeBounds(double &lower, double &upper) const;
    void calcEnergyBounds(double &lower, double &upper) const;

private:
    QSharedDataPointer<TrackPrivateData> data;
};

typedef struct {
    double abs_dist;
    long ngeocalls;
    int niter;
} ElemMutables;

typedef struct Element_s {
    G4String name;
    // Note: replicas not yet available
    G4ThreeVector offset;
    G4RotationMatrix rot;

    G4VSolid *solid;

    // To index element mutables and material properties
    int ecode;
    int matcode;

    // Is rotation matrix nontrivial
    bool rotated;

    // only three changeable fields
    bool visible;
    double alpha;

    // Alpha = 1.0 : opaque; alpha < 1.0, we do linear sequential color merging
    // (yes, resulting colors may be weird. Not as good as exponential
    // color influence falloff, but you can't have everything.
    // (the background color is white!)

    // statistics, frequently updated, nolock
    //    mutable long ngeocalls;
    //    // Caching for acceleration
    //    mutable int niter;
    //    mutable double abs_dist;

    std::vector<struct Element_s> children;
} Element;

typedef struct {
    const G4Material *mtl;
    double hue;
} MaterialInfo;

typedef struct ViewData_s {
    // What is being viewed
    Element elements;
    TrackData tracks;
    // How it is viewed
    G4ThreeVector camera;
    G4RotationMatrix orientation;
    std::vector<Plane> clipping_planes;
    G4double scale;
    G4double scene_radius;
    std::vector<MaterialInfo> matinfo;
    std::map<const G4Material *, int> matcode_map;
    bool split_by_material;
    // Simplification level
    int level_of_detail;
} ViewData;

class RenderWorker : public QObject {
    Q_OBJECT
public:
    RenderWorker();
    ~RenderWorker();
    bool abort_task;
public slots:
    bool render(ViewData p, TrackData t, QSharedPointer<QImage> i, int slice,
                int nslices);
    void coAbort();
    void flushAbort();
    void selfDestruct();
signals:
    void progressed(int);
    void completed();
    void aborted();

private:
    bool renderTracks(const ViewData &d, const TrackData &t, G4double *dists,
                      QRgb *colors, int slice, int nslices, int w, int h);
};

void countTree(const Element &e, int &treedepth, int &nelements);
int traceRay(const QPointF &scpt, const ViewData &d, const Element *hits[],
             Intersection ints[], int maxhits, int iteration,
             ElemMutables mutables[]);
int compressTraces(const Element *hits[], Intersection ints[], int m);
