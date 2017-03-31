#pragma once

#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>

#include <QObject>
#include <QRgb>
#include <QSharedData>

#include <vector>

class G4Material;
class G4VSolid;
class QProgressDialog;

typedef struct {
    G4ThreeVector normal;
    G4double offset;
    // Keep n*x >= o; Drop n*x < o
} Plane;

typedef struct Element_s {
    G4String name;
    // Note: replicas not yet available
    G4ThreeVector offset;
    G4RotationMatrix rot;

    G4Material *mat;
    G4VSolid *solid;
    bool rotated;

    // only three changeable fields
    bool visible;
    double hue;
    double alpha;
    // Alpha = 1.0 : opaque; alpha < 1.0, we do linear sequential color merging
    // (yes, resulting colors may be weird. Not as good as exponential
    // color influence falloff, but you can't have everything.
    // (the background color is white!)

    // statistics, frequently updated, nolock
    mutable long ngeocalls;
    // Caching for acceleration
    mutable int niter;
    mutable double abs_dist;

    std::vector<struct Element_s> children;
} Element;

typedef struct ViewData_s {
    Element root;
    G4ThreeVector camera;
    G4RotationMatrix orientation;
    std::vector<Plane> clipping_planes;
    G4double scale;
    G4double scene_radius;
    int level_of_detail;
} ViewData;

typedef struct {
    int npts;
    int ptype;
    size_t offset;
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
    void calcTimeBounds(double &lower, double &upper) const;
    void calcEnergyBounds(double &lower, double &upper) const;

private:
    QSharedDataPointer<TrackPrivateData> data;
};

class RenderWorker : public QObject {
    Q_OBJECT
public:
    RenderWorker();
    ~RenderWorker();
    bool abort_task;
public slots:
    bool render(ViewData p, TrackData t, QImage *i, int slice, int nslices,
                QProgressDialog *d = NULL);
    void coAbort();
    void flushAbort();
signals:
    void completed();
    void aborted();

private:
    bool renderTracks(const ViewData &d, const TrackData &t, G4double *dists,
                      QRgb *colors, int slice, int nslices, int w, int h);
};
