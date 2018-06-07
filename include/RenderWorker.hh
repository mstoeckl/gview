#pragma once

#include "General.hh"
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

typedef struct Element_s Element;
typedef struct ViewData_s ViewData;

typedef struct {
    G4ThreeVector normal;
    G4double offset;
    // Keep n*x >= o; Drop n*x < o
} Plane;

class CompactNormal {
public:
    CompactNormal() : x(0.), y(0.), z(0.) {}
    CompactNormal(G4ThreeVector n) : x(n.x()), y(n.y()), z(n.z()) {}
    operator G4ThreeVector() const { return G4ThreeVector(x, y, z); }

    float x;
    float y;
    float z;
};

typedef struct {
    // The normal at the intersection
    CompactNormal normal;
    // Distance from ray start to intersection
    G4double dist;
    // Volume behind the intersection
    const Element *elem;
} Intersection;

typedef struct {
    Intersection *intersections;
    int N;
    // TODO: record clipping plane index
    bool front_clipped;
    bool back_clipped;
} RayPoint;

typedef struct {
    int32_t npts;
    int32_t ptype;
    int64_t event_id;
    int32_t track_id;
    int32_t parent_id;
    uint64_t unused;
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

typedef union {
    TrackHeader h;
    TrackPoint p;
} TrackBlock;

typedef struct {
    double ballRadius;
} TrackMetaData;

class TrackPrivateData : public QSharedData {
public:
    size_t ntracks;
    size_t nblocks;
    TrackBlock *data;
    TrackMetaData *meta;

    TrackPrivateData(size_t itracks, size_t iblocks, TrackBlock *idata);
    explicit TrackPrivateData(const char *filename);
    TrackPrivateData(const TrackPrivateData &other);
    ~TrackPrivateData();

private:
    size_t mmapbytes;
};

class TrackData {
public:
    TrackData();
    TrackData(const char *filename);
    TrackData(const TrackData &other);
    TrackData(const TrackData &other, ViewData &viewrestr, Range seltimes,
              Range selenergies, IRange selidxs,
              const QMap<int, bool> &type_active);
    ~TrackData();
    size_t getNBlocks() const;
    size_t getNTracks() const;
    const TrackBlock *getBlocks() const;
    const TrackMetaData *getMeta() const;
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

    std::vector<Element> children;
} Element;

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
    bool force_opaque;
    // Simplification level
    int level_of_detail;
} ViewData;

class RenderRayTask : public RenderGraphNode {
public:
    RenderRayTask(QRect p);
    virtual void run(Context *) const;

private:
    QRect domain;
};

class RenderTrackTask : public RenderGraphNode {
public:
    RenderTrackTask(QRect p, int shard);
    virtual void run(Context *) const;

private:
    QRect domain;
    int shard;
};

class RenderMergeTask : public RenderGraphNode {
public:
    RenderMergeTask(QRect p);
    virtual void run(Context *) const;

private:
    QRect domain;
};

void countTree(const Element &e, int &treedepth, int &nelements);
G4ThreeVector forwardDirection(const G4RotationMatrix &);
G4ThreeVector initPoint(const QPointF &, const ViewData &);
RayPoint traceRay(const G4ThreeVector &init, const G4ThreeVector &forward,
                  const Element &root,
                  const std::vector<Plane> &clipping_planes,
                  Intersection *intersections, int maxhits, long iteration,
                  ElemMutables mutables[], bool first_visible_hit = false);
void compressTraces(RayPoint *pt);

class G4VPhysicalVolume;
Element convertCreation(const G4VPhysicalVolume *phys,
                        G4RotationMatrix rot = G4RotationMatrix(),
                        int *counter = NULL);
