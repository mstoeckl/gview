#pragma once

#include "RenderWorker.hh"

class LCData;
class LCNode;
class LineCollection;
class SimplexIterator {
public:
    SimplexIterator(const LineCollection &in, const G4ThreeVector &,
                    const G4ThreeVector &);
    // alternately, size_t*, int*, and N are a faster return method;
    // Requires a Freezing step for the Line Collection
    std::vector<std::pair<size_t, int>> getContainedLines() const;
    bool advance();

private:
    const LineCollection &lc;
    const LCData *current;
    const G4ThreeVector &src;
    const G4ThreeVector &dir;
    G4double t;
};

// Divide lines per BSP tree to minimize the maximum # of lines/node
class LineCollection {
public:
    LineCollection(const TrackHeader *, const TrackPoint *, size_t);
    ~LineCollection();
    SimplexIterator followRay(const G4ThreeVector &,
                              const G4ThreeVector &) const;

private:
    LCData *d;
    const TrackHeader *mheaders;
    const TrackPoint *mpoints;
    friend class SimplexIterator;
};
