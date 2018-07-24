/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#include "General.hh"

#include <QMap>
#include <QPointF>
#include <QSharedData>
#include <QVector>

class G4ParticleDefinition;
typedef struct ViewData_s ViewData;
typedef struct OctreeRoot_s OctreeRoot;

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

typedef union {
    TrackHeader h;
    TrackPoint p;
} TrackBlock;

typedef struct {
    double ballRadius;
    ushort generation;
} TrackMetaData;

typedef struct {
    Range energy;
    Range time;
    IRange seqno;
    IRange ngen;
    QMap<int32_t, bool> type_visible;
    QVector<QPair<int32_t, const G4ParticleDefinition *>> type_ids;
    QVector<const G4ParticleDefinition *> types;
} TrackRestriction;

class TrackPrivateData : public QSharedData {
public:
    size_t ntracks;
    size_t nblocks;
    TrackBlock *data;
    TrackMetaData *meta;
    OctreeRoot *octree;

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
    TrackData(const TrackData &other, const ViewData &viewrestr,
              const TrackRestriction &restr);
    ~TrackData();
    size_t getNBlocks() const;
    size_t getNTracks() const;
    const TrackBlock *getBlocks() const;
    const TrackMetaData *getMeta() const;
    const OctreeRoot *getOctree() const;
    const QMap<int32_t, const G4ParticleDefinition *> calcTypes() const;
    void calcTimeBounds(double &lower, double &upper) const;
    void calcEnergyBounds(double &lower, double &upper) const;
    void constructRangeHistograms(QVector<QPointF> &tp, QVector<QPointF> &ep,
                                  const Range &tr, const Range &er) const;
    int calcMaxGenerations() const;

private:
    QSharedDataPointer<TrackPrivateData> data;
};
