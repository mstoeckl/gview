/* SPDX-License-Identifier: GPL-3.0-only */
#include "Navigator.hh"

#include <G4Box.hh>
#include <G4Navigator.hh>
#include <G4PVPlacement.hh>
#include <G4StateManager.hh>
#include <G4VExceptionHandler.hh>
#include <G4VSolid.hh>

#include <QThread>

static G4ThreeVector condrot(const Element &e, const G4ThreeVector &vec) {
    if (e.rotated) {
        return e.rot * vec;
    }
    return vec;
}

static G4ThreeVector condirot(const Element &e, const G4ThreeVector &vec) {
    if (e.rotated) {
        return e.rot.inverse() * vec;
    }
    return vec;
}

static RayPoint ftraceRay(const G4ThreeVector &init,
                          const G4ThreeVector &forward,
                          const std::vector<Element> &els,
                          const std::vector<Plane> &clipping_planes,
                          Intersection *ints, int maxhits, long iteration,
                          ElemMutables mutables[], bool first_visible_hit) {
    RayPoint ret;
    ret.N = 0;
    ret.intersections = ints;
    ret.front_clipped = false;
    ret.back_clipped = false;

    // Pseudorandom number to prevent consistent child ordering
    // so as to make overlapping children obvious
    const size_t rotfact = size_t(random());
    // Minimum feature size, yet still above double discrimination threshold
    //    const G4double epsilon = 1e-3 * CLHEP::nanometer;

    // We set edist as kInfinity/4 so that later increments of kInfinity/2
    // (which is less than the distanceToIn failure result) surpass it.
    G4double sdist = 0., edist = kInfinity / 4;
    G4ThreeVector entrynormal;
    G4ThreeVector exitnormal;
    bool exists = clipRay(clipping_planes, init, forward, sdist, edist,
                          entrynormal, exitnormal);
    if (!exists) {
        return ret;
    }

    G4ThreeVector start = init + forward * sdist;
    const Element &root = els[0];
    ++mutables[root.ecode].ngeocalls;
    // Ensure ray starts in the root.
    bool clippable = false;
    if (!root.solid->Inside(start)) {
        G4double jdist = root.solid->DistanceToIn(start, forward);
        if (jdist + sdist >= edist) {
            // Root solid not reachable with available distance
            return ret;
        }
        sdist += jdist;
        start += forward * jdist;
        clippable = false;
    } else {
        clippable = true;
    }

    const int maxdepth = 10;
    // Assume no depth greater than depth
    // TODO: make all such buffers created but once and
    // passed in, so that stack elements stay close

    const double epsilon = 1e-6 * CLHEP::nm;

    // Assume we are already in the root element.
    // We statically allocate, because many new/frees are too expensive
    const Element *stack[maxdepth];
    stack[0] = &root;
    size_t n = 1;
    G4ThreeVector local = start;
    for (int iter = 0; iter < 1000; iter++) {
        const Element &last = *stack[n - 1];

        // Strictly inclusion tests....
        // check inside on all children, append; stop when no longer dropping
        // levels
        bool found = false;
        for (size_t walk = 0; walk < last.children.size(); walk++) {
            // rotate/offset start point...
            const Element &elem =
                els[last.children[(walk + rotfact) % last.children.size()]];
            ++mutables[elem.ecode].ngeocalls;
            if (elem.solid->Inside(condrot(elem, (local + elem.offset)))) {
                stack[n] = &elem;
                ++n;
                found = true;
                local += elem.offset;
                // Picking the first intersection can lead to visual glitches
                // when two or more children overlap at the same point.
                // One solution: pick a random child order..
                break;
            }
        }
        if (!found) {
            break;
        }
    }
    G4ThreeVector offsets[maxdepth];
    offsets[0] = start;
    for (size_t i = 1; i < n; i++) {
        const Element &elem = *stack[i];
        offsets[i] = offsets[i - 1] + elem.offset;
    }
    const G4double isdist = sdist;

    if (clippable) {
        // Start point is in the world volume, and is an intersection
        if (!first_visible_hit || stack[n - 1]->visible) {
            ints[0].dist = sdist;
            ints[0].normal = entrynormal;
            ints[0].ecode = stack[n - 1]->ecode;
            ret.front_clipped = true;
            ret.N++;
            if (first_visible_hit) {
                return ret;
            }
        }
    } else {
        // Intersection on entry
        if (!first_visible_hit || root.visible) {
            ++mutables[root.ecode].ngeocalls;
            G4ThreeVector lnormal = root.solid->SurfaceNormal(local);
            ints[0].dist = sdist;
            ints[0].normal = lnormal;
            ints[0].ecode = root.ecode;
            ret.front_clipped = false;
            ret.N++;
            if (first_visible_hit) {
                return ret;
            }
        }
    }

    for (int iter = 0; iter < 1000; iter++) {
        const Element &curr = *stack[n - 1];
        const G4ThreeVector &pos = offsets[n - 1] + forward * (sdist - isdist);

        ElemMutables &cmu = mutables[curr.ecode];
        if (cmu.niter != iteration || cmu.abs_dist <= sdist /* epsilon ? */) {
            ++cmu.ngeocalls;
            cmu.abs_dist =
                sdist + curr.solid->DistanceToOut(condrot(curr, pos),
                                                  condrot(curr, forward));
            cmu.niter = iteration;
        }
        G4double exitdist = std::min(kInfinity / 4, cmu.abs_dist - sdist);

        // Typically, no-collision implies a distance of kInfinity,
        // but floating point error may reduce this slightly. Thus if
        // the distance-to-in is kInfinity, we automatically get no-intersection
        G4double closestDist = kInfinity / 2;
        const Element *closest = NULL;
        G4ThreeVector closestPos;
        for (size_t walk = 0; walk < curr.children.size(); walk++) {
            const Element &elem =
                els[curr.children[(walk + rotfact) % curr.children.size()]];
            G4ThreeVector sub = pos + elem.offset;

            ElemMutables &emu = mutables[elem.ecode];
            if (emu.niter != iteration || emu.abs_dist <= sdist) {
                ++emu.ngeocalls;
                // Inside case happens when shapes are way out of bounds.
                // (Note: ought to make this optional, since I doubt
                //  Geant's regular walking handles this.)
                if (sdist > emu.abs_dist &&
                    elem.solid->Inside(condrot(elem, sub))) {
                    emu.abs_dist = sdist;
                } else {
                    if (sdist < emu.abs_dist) {
                        ++emu.ngeocalls;
                    }
                    emu.abs_dist =
                        sdist + elem.solid->DistanceToIn(
                                    condrot(elem, sub), condrot(elem, forward));
                }
                emu.niter = iteration;
            }
            G4double altdist = emu.abs_dist - sdist;
            if (altdist < closestDist) {
                closestPos = sub;
                closest = &elem;
                closestDist = altdist;
            }
        }
        G4double fdist = std::min(exitdist, closestDist);
        if (fdist > edist - sdist) {
            // End clip (counts as a hit!; always relevant)
            ints[ret.N].dist = edist;
            ints[ret.N].normal = exitnormal;
            ints[ret.N].ecode = CODE_END;
            ret.back_clipped = true;
            ret.N++;
            return ret;
        } else if (exitdist < closestDist) {
            // Transition on leaving vol w/o intersections
            sdist += exitdist;
            G4ThreeVector lpos = offsets[n - 1] + forward * (sdist - isdist);
            // Drop from stack
            --n;

            // Record hit with normal
            bool store = !first_visible_hit || (n > 0 && stack[n - 1]->visible);
            if (store) {
                ++cmu.ngeocalls;
                G4ThreeVector lnormal =
                    curr.solid->SurfaceNormal(condrot(curr, lpos));
                if (ret.N > 0 && sdist <= ints[ret.N - 1].dist + epsilon) {
                    // On zero-width gap, replace with latest volume
                    ints[ret.N - 1].dist = (sdist + ints[ret.N - 1].dist) / 2;
                    ints[ret.N - 1].normal = condirot(curr, lnormal);
                    ints[ret.N - 1].ecode =
                        n > 0 ? stack[n - 1]->ecode : CODE_END;
                } else {
                    // New volume
                    ints[ret.N].dist = sdist;
                    ints[ret.N].normal = condirot(curr, lnormal);
                    ints[ret.N].ecode = n > 0 ? stack[n - 1]->ecode : CODE_END;
                    ret.N++;
                }
                if (n == 0 || ret.N >= maxhits || first_visible_hit) {
                    return ret;
                }
            } else if (n == 0) {
                return ret;
            }
        } else {
            // Transition on visiting a child
            stack[n] = closest;
            offsets[n] = closest->offset + offsets[n - 1];
            sdist += closestDist;
            n++;
            G4ThreeVector lpos = offsets[n - 1] + forward * (sdist - isdist);

            // Record hit with normal
            bool store = !first_visible_hit || closest->visible;
            if (store) {
                ++mutables[closest->ecode].ngeocalls;
                G4ThreeVector lnormal =
                    closest->solid->SurfaceNormal(condrot(*closest, lpos));
                if (ret.N > 0 && sdist <= ints[ret.N - 1].dist + epsilon) {
                    // On zero-width gap, replace with latest volume
                    ints[ret.N - 1].dist = (sdist + ints[ret.N - 1].dist) / 2;
                    ints[ret.N - 1].normal = condirot(*closest, lnormal);
                    ints[ret.N - 1].ecode = closest->ecode;
                } else {
                    // New volume
                    ints[ret.N].dist = sdist;
                    ints[ret.N].normal = condirot(*closest, lnormal);
                    ints[ret.N].ecode =
                        (ret.N >= maxhits) ? CODE_END : closest->ecode;
                    ret.N++;
                }
                if (ret.N >= maxhits || first_visible_hit) {
                    return ret;
                }
            }
        }
    }
    return ret;
}

Navigator::~Navigator() {}
Navigator *Navigator::create(const ViewData &vd, int navtype) {
    switch (navtype) {
    case nFastVolNav:
        return new FastVolNavigator(vd.elements, vd.clipping_planes);
    case nGeantNav:
        return new GeantNavigator(vd.orig_vol, vd.elements, vd.clipping_planes);
    default:
    case nVoxelNav:
        return NULL;
    }
}

FastVolNavigator::FastVolNavigator(const std::vector<Element> &iels,
                                   const std::vector<Plane> &iclipping_planes)
    : els(iels), clipping_planes(iclipping_planes), mutables(nullptr),
      iteration(1L) {
    mutables = new ElemMutables[els.size()]();
}
FastVolNavigator::~FastVolNavigator() { delete[] mutables; }
RayPoint FastVolNavigator::traceRay(const G4ThreeVector &init,
                                    const G4ThreeVector &forward,
                                    Intersection *intersections, int maxhits,
                                    bool first_visible_hit) {
    RayPoint rpt = ftraceRay(init, forward, els, clipping_planes, intersections,
                             maxhits, iteration, mutables, first_visible_hit);
    iteration++;
    return rpt;
}

class GExh : public G4VExceptionHandler {
public:
    GExh() {}

    virtual G4bool Notify(const char *originOfException,
                          const char *exceptionCode,
                          G4ExceptionSeverity severity,
                          const char *description) {
        (void)originOfException;
        (void)exceptionCode;
        (void)severity;
        (void)description;
        // TODO: logarithmic duplicate counting
        // print on 1st, 10th, 100th, etc.

        //        qDebug("EXC: %s | %s | %d | %s", originOfException,
        //        exceptionCode,
        //               severity, description);
        return false;
    }
};

GeantNavigator::GeantNavigator(G4VPhysicalVolume *iworld,
                               const std::vector<Element> &iels,
                               const std::vector<Plane> &iclipping_planes)
    : world(iworld), els(iels), clipping_planes(iclipping_planes) {
    if (!iworld) {
        qFatal("Passed in null world");
    }

    // Set the threadlocal exceptionhandler
    G4VExceptionHandler *ohandle =
        G4StateManager::GetStateManager()->GetExceptionHandler();
    if (ohandle) {
        delete ohandle;
    }
    G4StateManager::GetStateManager()->SetExceptionHandler(new GExh());

    // Ensure that threadlocals have copies in this thread?
    G4PVManager *pvm =
        const_cast<G4PVManager *>(&G4VPhysicalVolume::GetSubInstanceManager());
    pvm->SlaveCopySubInstanceArray();
    G4LVManager *lvm =
        const_cast<G4LVManager *>(&G4LogicalVolume::GetSubInstanceManager());
    lvm->SlaveCopySubInstanceArray();

    for (const Element &n : els) {
        ecode_map[n.orig_vol] = n.ecode;
    }

    nav = new G4Navigator();
    nav->SetWorldVolume(world);
}
GeantNavigator::~GeantNavigator() { delete nav; }

static const char *volname(G4VPhysicalVolume *v) {
    return v ? v->GetName().c_str() : "(null)";
}

RayPoint GeantNavigator::traceRay(const G4ThreeVector &init,
                                  const G4ThreeVector &forward,
                                  Intersection *intersections, int maxhits,
                                  bool first_visible_hit) {
    // The key bit
    nav->ResetStackAndState();

    RayPoint r;
    r.N = 0;
    r.front_clipped = false;
    r.back_clipped = false;
    r.intersections = intersections;
    G4double sdist = 0., edist = kInfinity / 2;
    G4ThreeVector entrynormal;
    G4ThreeVector exitnormal;
    bool exists = clipRay(clipping_planes, init, forward, sdist, edist,
                          entrynormal, exitnormal);
    if (!exists) {
        return r;
    }

    {
        const G4ThreeVector first_pos = init + sdist * forward;
        const G4VSolid *root_solid = world->GetLogicalVolume()->GetSolid();
        EInside init_pos = root_solid->Inside(first_pos);
        if (init_pos == kInside) {
            // We are already inside the world. find out which subvolume
            G4VPhysicalVolume *last_vol = nav->LocateGlobalPointAndSetup(
                init + sdist * forward, &forward);
            int ec = ecode_map.value(last_vol, CODE_END);
            if (!first_visible_hit || (ec >= 0 && els[ec].visible)) {
                Intersection &ins = r.intersections[r.N++];
                ins.dist = sdist;
                ins.normal = entrynormal;
                ins.ecode = ec;
                r.front_clipped = true;
                if (first_visible_hit) {
                    return r;
                }
            }
        } else {
            // Are outside the world. Must perform advance-in by hand.
            sdist += root_solid->DistanceToIn(first_pos, forward);
            if (sdist >= edist) {
                // nope. end clip hit first
                return r;
            }
            // determine the region we just entered
            G4VPhysicalVolume *last_vol = nav->LocateGlobalPointAndSetup(
                init + sdist * forward, &forward);
            int ec = ecode_map.value(last_vol, CODE_END);
            if (!first_visible_hit || (ec >= 0 && els[ec].visible)) {
                Intersection &ins = r.intersections[r.N++];
                ins.dist = sdist;
                ins.normal = CompactNormal(
                    root_solid->SurfaceNormal(init + sdist * forward));
                ins.ecode = ec;
                if (first_visible_hit) {
                    return r;
                }
            }
        }
    }

    // Now we begin, definitely inside the world
    for (int i = 0; i < 100; i++) {
        /* After any position change, we locate ourselves, and then
         * determine forward step & local normal */

        G4double safety;
        bool valid;
        const G4ThreeVector here = init + sdist * forward;
        double length = nav->ComputeStep(here, forward, kInfinity, safety);
        // step was always geometrically limited
        nav->SetGeometricallyLimitedStep();
        const G4ThreeVector &normal = nav->GetGlobalExitNormal(here, &valid);
        G4VPhysicalVolume *v = nav->LocateGlobalPointAndSetup(here, &forward);

        sdist += length;

        if (sdist > edist) {
            int nc = ecode_map.value(v, CODE_END);
            if (!first_visible_hit || (nc >= 0 && els[nc].visible)) {
                Intersection &ins = r.intersections[r.N++];
                ins.dist = sdist;
                ins.normal = exitnormal;
                ins.ecode = nc;
                r.back_clipped = true;
                if (r.N >= maxhits || first_visible_hit) {
                    return r;
                }
            }
            break;
        } else if (length > 0) {
            int nc = ecode_map.value(v, CODE_END);
            // points out of previous volume, into current volume
            if (!first_visible_hit || (nc >= 0 && els[nc].visible)) {
                Intersection &ins = r.intersections[r.N++];
                ins.dist = sdist;
                ins.normal = CompactNormal(-normal);
                ins.ecode = nc;
                if (r.N >= maxhits || first_visible_hit) {
                    return r;
                }
            }
        }

        if (v == NULL) {
            return r;
        }
    }

    return r;
}
