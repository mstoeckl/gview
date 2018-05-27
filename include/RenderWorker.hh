#pragma once

#include "RenderGraph.hh"

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
    bool is_clipping_plane;
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

typedef struct {
    double low, high;
} Range;
typedef struct {
    size_t low, high;
} IRange;

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
    void constructRangeHistograms(QVector<QPointF> &tp, QVector<QPointF> &ep,
                                  const Range &tr, const Range &er) const;

private:
    QSharedDataPointer<TrackPrivateData> data;
};

typedef struct {
    double abs_dist;
    long ngeocalls;
    long niter;
} ElemMutables;

typedef struct Element_s {
    G4String name;
    // Note: replicas not yet available
    G4ThreeVector offset;
    G4RotationMatrix rot;

    const G4VSolid *solid;
    const G4Material *material;

    // To index element mutables and color properties
    int ecode;
    int ccode;

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

// Very simple and fast color handling
class VColor {
public:
    VColor(uint8_t r, uint8_t g, uint8_t b) : pr(r), pg(g), pb(b) {}
    VColor(QRgb rgb) : pr(qRed(rgb)), pg(qGreen(rgb)), pb(qBlue(rgb)) {}
    static VColor fromRgbF(float r, float g, float b) {
        return VColor(255.f * r, 255.f * g, 255.f * b);
    }
    QRgb rgb() const { return qRgb(pr, pg, pb); }
    float redF() const { return pr / 255.f; }
    float greenF() const { return pg / 255.f; }
    float blueF() const { return pb / 255.f; }

private:
    uint8_t pr, pg, pb;
};

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
    std::vector<VColor> color_table;
    bool split_by_material;
    // Simplification level
    int level_of_detail;
} ViewData;

class RenderRayTask : public RenderGraphTask {
public:
    RenderRayTask(QRect p, RenderGraph &h, QSharedPointer<Context> c, int id);
    virtual void run();
};

class RenderTrackTask : public RenderGraphTask {
public:
    RenderTrackTask(QRect p, RenderGraph &h, QSharedPointer<Context> c, int id);
    virtual void run();
};

void countTree(const Element &e, int &treedepth, int &nelements);
G4ThreeVector forwardDirection(const G4RotationMatrix &);
G4ThreeVector initPoint(const QPointF &, const ViewData &);
int traceRay(const G4ThreeVector &init, const G4ThreeVector &forward,
             const ViewData &d, const Element *hits[], Intersection ints[],
             int maxhits, long iteration, ElemMutables mutables[]);
int compressTraces(const Element *hits[], Intersection ints[], int m);

class G4VPhysicalVolume;
Element convertCreation(const G4VPhysicalVolume *phys,
                        G4RotationMatrix rot = G4RotationMatrix(),
                        int *counter = NULL);
