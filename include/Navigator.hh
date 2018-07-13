/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#include "RenderWorker.hh"

#include <QMap>

class G4Navigator;

#define CODE_LINE -1
#define CODE_END -2

enum {
    nFastVolNav = 0,
    nGeantNav = 1,
    nVoxelNav = 2,
};

/* Get constructed per render task, and permit raytrace exchanges */
class Navigator {
public:
    virtual ~Navigator();
    /* Returns a RayPoint which uses `intersections` for storage. */
    virtual RayPoint traceRay(const G4ThreeVector &init,
                              const G4ThreeVector &forward,
                              Intersection *intersections, int maxhits,
                              bool first_visible_hit = false) = 0;
    static Navigator *create(const ViewData &vd, int navtype);
};

class FastVolNavigator : public Navigator {
public:
    FastVolNavigator(const std::vector<Element> &els,
                     const std::vector<Plane> &clipping_planes);
    virtual ~FastVolNavigator();
    RayPoint traceRay(const G4ThreeVector &init, const G4ThreeVector &forward,
                      Intersection *intersections, int maxhits,
                      bool first_visible_hit = false) override;

private:
    const std::vector<Element> &els;
    const std::vector<Plane> &clipping_planes;
    ElemMutables *mutables;
    long iteration;
};

class GeantNavigator : public Navigator {
public:
    GeantNavigator(G4VPhysicalVolume *world, const std::vector<Element> &els,
                   const std::vector<Plane> &clipping_planes);
    virtual ~GeantNavigator();
    RayPoint traceRay(const G4ThreeVector &init, const G4ThreeVector &forward,
                      Intersection *intersections, int maxhits,
                      bool first_visible_hit = false) override;

private:
    G4VPhysicalVolume *world;
    G4Navigator *nav;
    const std::vector<Element> &els;
    const std::vector<Plane> &clipping_planes;
    QMap<const G4VPhysicalVolume *, int> ecode_map;
};
