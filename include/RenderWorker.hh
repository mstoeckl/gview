/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#include "General.hh"
#include "RenderGraph.hh"
#include "TrackData.hh"

#include <G4RotationMatrix.hh>
#include <G4ThreeVector.hh>

#include <QObject>
#include <QRgb>
#include <QSharedData>
#include <QSharedPointer>

#include <vector>

enum {
    CHANGE_ONESHOT = 1 << 1,
    CHANGE_TRACK = 1 << 2,
    CHANGE_GEO = 1 << 3,
    CHANGE_COLOR = 1 << 4,
    CHANGE_VIEWPORT = CHANGE_GEO | CHANGE_TRACK | CHANGE_COLOR
};

class LineCollection;
class G4Material;
class G4VSolid;
class QProgressDialog;
class Navigator;
class G4VPhysicalVolume;

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
    bool operator==(const CompactNormal &o) const {
        return x == o.x && y == o.y && z == o.z;
    }
    bool operator!=(const CompactNormal &o) const { return !(o == *this); }
    bool operator<(const CompactNormal &o) const {
        return x < o.x || (x == o.x && (y < o.y || (y == o.y && z < o.z)));
    }

    float x;
    float y;
    float z;
};

typedef struct Intersection_s {
    // The normal at the intersection
    CompactNormal normal;
    // Volume behind the intersection
    int ecode;
    // Distance from ray start to intersection
    G4double dist;
} Intersection;

typedef struct RayPoint_s {
    Intersection *intersections;
    int N;
    // TODO: record clipping plane index
    bool front_clipped;
    bool back_clipped;
} RayPoint;

typedef struct {
    double abs_dist;
    long ngeocalls;
    long niter;
} ElemMutables;

typedef struct Element_s {
    G4String name;
    // Note: replicas not yet available
    G4ThreeVector offset; // local
    G4RotationMatrix rot; // global

    const G4VPhysicalVolume *orig_vol;
    const G4VSolid *solid;
    const G4Material *material;
    double cubicVolume;
    double surfaceArea;

    // To index element mutables and color properties
    int ecode;

    // Is rotation matrix nontrivial
    bool rotated;

    // Display control fields
    int ccode;
    bool visible;
    double alpha;

    // Relative to the root node
    G4ThreeVector global_offset;

    // Alpha = 1.0 : opaque; alpha < 1.0, we do linear sequential color
    // merging (yes, resulting colors may be weird. Not as good as
    // exponential color influence falloff, but you can't have everything.
    // (the background color is white!)

    std::vector<int> children;
} Element;

typedef struct ViewData_s {
    G4VPhysicalVolume *orig_vol;
    void *octree;
    // What is being viewed. elements[0] is root
    std::vector<Element> elements;
    TrackData tracks;
    // How it is viewed
    G4ThreeVector camera;
    G4RotationMatrix orientation;
    // Frame for rotations
    G4ThreeVector base_offset;
    G4RotationMatrix base_rotation;
    G4double scale;
    G4double scene_radius;
    std::vector<Plane> clipping_planes;
    std::vector<VColor> color_table;
    double voxel_base_density;
    bool split_by_material;
    int navigator, gshader, tshader;
    bool force_opaque;
    // Simplification level
    int level_of_detail;
} ViewData;

void countTree(const std::vector<Element> &els, int index, int &treedepth,
               int &nelements);
G4ThreeVector forwardDirection(const G4RotationMatrix &);
G4ThreeVector initPoint(const QPointF &, const ViewData &);
bool clipRay(const std::vector<Plane> &clipping_planes,
             const G4ThreeVector &init, const G4ThreeVector &forward,
             G4double &sdist, G4double &edist, G4ThreeVector &entrynormal,
             G4ThreeVector &exitnormal);
RayPoint rayAtPoint(Navigator &nav, const QPointF &pt, qreal radius,
                    const G4ThreeVector &forward, const ViewData &d,
                    Intersection *ints, Intersection *altints, int M,
                    int *ndevs);
void debugRayPoint(const RayPoint &ray, const std::vector<Element> &els);

class G4VPhysicalVolume;
int convertCreation(std::vector<Element> &elts, const G4VPhysicalVolume *phys,
                    const G4String &suffix,
                    const G4RotationMatrix &parent_rot = G4RotationMatrix(),
                    const G4ThreeVector &parent_offset = G4ThreeVector(),
                    int *counter = NULL);
long recursivelySumNCalls(const std::vector<Element> &elts,
                          const ElemMutables e[]);
void recursivelyPrintNCalls(const std::vector<Element> &elts,
                            const ElemMutables e[], int depth = 0, long net = 0,
                            int idx = 0);

QPair<double, QString> ruler_distance(double scale, double max_width_frac);
