/* SPDX-License-Identifier: GPL-3.0-only */
#include "TrackData.hh"

#include "RenderWorker.hh"
#include <geomdefs.hh>

#include <G4Deuteron.hh>
#include <G4Electron.hh>
#include <G4Gamma.hh>
#include <G4GenericIon.hh>
#include <G4Neutron.hh>
#include <G4OpticalPhoton.hh>
#include <G4Positron.hh>
#include <G4Proton.hh>
#include <G4Triton.hh>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

TrackData::TrackData() {}

TrackData::TrackData(const char *filename) {
    data = QSharedDataPointer<TrackPrivateData>(new TrackPrivateData(filename));
}

TrackData::TrackData(const TrackData &other) : data(other.data) {}

static void lineCuts(const std::vector<Plane> &clips, G4ThreeVector from,
                     G4ThreeVector to, G4double &lowE, G4double &highE) {
    G4double low = -kInfinity;
    G4double high = kInfinity;
    for (const Plane &p : clips) {
        G4double fpos = p.offset - p.normal * from;
        G4double tpos = p.offset - p.normal * to;
        if (fpos - tpos == 0.) {
            // TODO sign may be off
            if (fpos > 0) {
                // On cut side of parallel plane
                low = kInfinity;
                high = -kInfinity;
            } else {
                // No change, inside simplex
            }
            continue;
        }
        // [0,1] becomes [fpos,tpos]: inv img of 0.
        // x' = (x - 0) * (fpos - tpos) + fpos
        // x = (x' - fpos) /  (fpos - tpos)
        G4double pos = fpos / (fpos - tpos);
        if (tpos < fpos) {
            low = std::max(pos, low);
        } else {
            high = std::min(pos, high);
        }
    }
    lowE = low;
    highE = high;
}

static bool insideConvex(const std::vector<Plane> &clips, G4ThreeVector point) {
    // Optimizer will resolve to simple form anyway; & this is error free
    double lowE, highE;
    lineCuts(clips, point, point, lowE, highE);
    return lowE < highE;
}

static TrackPoint linmix(const TrackPoint &a, const TrackPoint &b, double mx) {
    TrackPoint c;
    double dx = 1.0 - mx;
    c.x = dx * a.x + mx * b.x;
    c.y = dx * a.y + mx * b.y;
    c.z = dx * a.z + mx * b.z;
    c.time = float(dx) * a.time + float(mx) * b.time;
    c.energy = float(dx) * a.energy + float(mx) * b.energy;
    return c;
}

TrackData::TrackData(const TrackData &other, const ViewData &vd,
                     const TrackRestriction &rstr) {
    size_t otracks = other.getNTracks();
    const TrackBlock *oblocks = other.getBlocks();
    const TrackMetaData *ometa = other.getMeta();
    const Range &trange = rstr.time;
    const Range &erange = rstr.energy;
    const IRange &selidxs = rstr.seqno;
    const IRange &genrange = rstr.ngen;
    const QMap<int32_t, bool> type_active = rstr.type_visible;

    // Worst case allocation is one track (header+2 pts) per line segment
    TrackBlock *buf =
        (TrackBlock *)malloc(3 * other.getNBlocks() * sizeof(TrackBlock));
    size_t qtracks = 0;
    size_t qblocks = 0;

    float tlow = float(trange.low), thigh = float(trange.high);
    float elow = float(erange.low), ehigh = float(erange.high);
    size_t nlow = std::max(size_t(0), selidxs.low - 1);
    size_t nhigh = std::min(otracks, selidxs.high);
    // Scan up to track (selidxs.low-1)

    size_t i = 0;
    for (size_t z = 0; z < nlow; z++) {
        i += oblocks[i].h.npts + 1;
    }

    for (size_t z = nlow; z < nhigh; z++) {
        const TrackHeader &oheader = oblocks[i].h;
        const TrackPoint *seq = (TrackPoint *)&oblocks[i + 1];
        i += oheader.npts + 1;
        size_t qheader = -1;

        bool typekeep = type_active.value(oheader.ptype, false);
        if (!typekeep || (genrange.low > ometa[z].generation ||
                          genrange.high < ometa[z].generation)) {
            continue;
        }

        G4ThreeVector fts(seq[0].x, seq[0].y, seq[0].z);
        bool intime = seq[0].time >= tlow && seq[0].time <= thigh;
        bool inenergy = seq[0].energy >= elow && seq[0].energy <= ehigh;
        bool inconvex = insideConvex(vd.clipping_planes, fts);
        bool started = false;
        if (inconvex && intime && inenergy) {
            TrackHeader h = oheader;
            h.npts = 1;
            qtracks++;
            qheader = qblocks;
            buf[qblocks].h = h;
            buf[qblocks + 1].p = seq[0];
            qblocks += 2;
            started = true;
        }

        for (int j = 0; j < oheader.npts - 1; j++) {
            double low, high;
            G4ThreeVector pl(seq[j].x, seq[j].y, seq[j].z);
            G4ThreeVector ph(seq[j + 1].x, seq[j + 1].y, seq[j + 1].z);
            lineCuts(vd.clipping_planes, pl, ph, low, high);
            float itrange = (seq[j + 1].time - seq[j].time == 0.f)
                                ? float(kInfinity)
                                : 1 / (seq[j + 1].time - seq[j].time);
            float tcutlow = (tlow - seq[j].time) * itrange;
            float tcuthigh = (thigh - seq[j].time) * itrange;
            float ierange = (seq[j + 1].energy - seq[j].energy == 0.f)
                                ? float(kInfinity)
                                : 1 / (seq[j + 1].energy - seq[j].energy);
            float ecutlow = (elow - seq[j].energy) * ierange;
            float ecuthigh = (ehigh - seq[j].energy) * ierange;
            low = std::max(double(std::max(tcutlow, ecutlow)), low);
            high = std::min(double(std::min(tcuthigh, ecuthigh)), high);

            // TODO: modify casework to localize results, not inputs
            if (low <= 0. && high >= 1.) {
                // Keep point
                if (!started) {
                    TrackHeader h = oheader;
                    h.npts = 1;
                    qheader = qblocks;
                    buf[qblocks].h = h;
                    buf[qblocks + 1].p = seq[j];
                    qblocks += 2;
                }
                buf[qheader].h.npts++;
                buf[qblocks].p = seq[j + 1];
                qblocks++;
            } else if (low <= 0. && high < 1.) {
                if (high > 0.) {
                    if (!started) {
                        TrackHeader h = oheader;
                        h.npts = 1;
                        qtracks++;
                        qheader = qblocks;
                        buf[qblocks].h = h;
                        buf[qblocks + 1].p = seq[j];
                        qblocks += 2;
                    }
                    buf[qheader].h.npts++;
                    buf[qblocks].p = linmix(seq[j], seq[j + 1], high);
                    qblocks++;
                } else {
                    // Not there entirely
                }
            } else if (low > 0. && high >= 1.) {
                if (low < 1.) {
                    // Starting new sequence
                    TrackHeader h = oheader;
                    h.npts = 2;
                    qtracks++;
                    qheader = qblocks;
                    buf[qblocks].h = h;
                    buf[qblocks + 1].p = linmix(seq[j], seq[j + 1], low);
                    buf[qblocks + 2].p = seq[j + 1];
                    qblocks += 3;
                } else {
                    // Not there entirely
                }
            } else {
                if (high < 1. && low > 0. && low < high) {
                    // low > 0, high < 1
                    // A whole new short segment...
                    TrackHeader h = oheader;
                    h.npts = 2;
                    qtracks++;
                    qheader = qblocks;
                    buf[qblocks].h = h;
                    buf[qblocks + 1].p = linmix(seq[j], seq[j + 1], low);
                    buf[qblocks + 2].p = linmix(seq[j], seq[j + 1], high);
                    qblocks += 3;
                }
            }
        }
    }

    // Shrink buffer as necessary
    buf = (TrackBlock *)realloc(buf, qblocks * sizeof(TrackBlock));
    TrackPrivateData *pd = new TrackPrivateData(qtracks, qblocks, buf);
    data = QSharedDataPointer<TrackPrivateData>(pd);
}

TrackData::~TrackData() {}

size_t TrackData::getNBlocks() const {
    if (!data) {
        return 0;
    }
    return data.constData()->nblocks;
}
size_t TrackData::getNTracks() const {
    if (!data) {
        return 0;
    }
    return data.constData()->ntracks;
}
const TrackBlock *TrackData::getBlocks() const {
    if (!data) {
        return NULL;
    }
    return data.constData()->data;
}

const TrackMetaData *TrackData::getMeta() const {
    if (!data) {
        return NULL;
    }
    return data.constData()->meta;
}

void TrackData::calcTimeBounds(double &lower, double &upper) const {
    upper = -kInfinity;
    lower = kInfinity;
    if (!data) {
        return;
    }
    TrackBlock *blocks = data.constData()->data;
    size_t ntracks = data.constData()->ntracks;
    size_t i = 0;
    for (size_t k = 0; k < ntracks; k++) {
        size_t npts = blocks[i].h.npts;
        i++;
        for (size_t j = 0; j < npts; j++) {
            upper = std::fmax(blocks[i].p.time, upper);
            lower = std::fmin(blocks[i].p.time, lower);
            i++;
        }
    }
}
void TrackData::calcEnergyBounds(double &lower, double &upper) const {
    upper = -kInfinity;
    lower = kInfinity;
    if (!data) {
        return;
    }
    TrackBlock *blocks = data.constData()->data;
    size_t ntracks = data.constData()->ntracks;
    size_t i = 0;
    for (size_t k = 0; k < ntracks; k++) {
        size_t npts = blocks[i].h.npts;
        i++;
        for (size_t j = 0; j < npts; j++) {
            float e = blocks[i].p.energy;
            upper = std::fmax(e, upper);
            if (e > 0.) {
                lower = std::fmin(e, lower);
            }
            i++;
        }
    }
}
int TrackData::calcMaxGenerations() const {
    const TrackPrivateData *tpd = data.constData();
    size_t ntracks = tpd->ntracks;
    TrackMetaData *meta = tpd->meta;
    TrackBlock *blocks = tpd->data;
    int maxgen = 0;
    size_t *index = new size_t[ntracks];
    size_t i = 0;
    for (size_t z = 0; z < ntracks; z++) {
        index[z] = i;
        i += blocks[i].h.npts + 1;
    }

    for (size_t z = 0; z < ntracks; z++) {
        size_t current = z;
        int ng = 1;
        while (current > 0 && blocks[index[current]].h.parent_id > 0) {
            size_t ncurrent = (size_t)-1;
            for (size_t y = current; y > 0; y--) {
                int hi = index[y - 1];
                if (blocks[hi].h.track_id ==
                    blocks[index[current]].h.parent_id) {
                    ncurrent = y - 1;
                    break;
                }
            }

            if (ncurrent == (size_t)-1) {
                break;
            }
            current = ncurrent;

            ng++;
        }
        meta[z].generation = (ushort)std::min(ng, 65535);
        maxgen = std::max(ng, maxgen);
    }
    delete[] index;
    return maxgen;
}
const QMap<int32_t, const G4ParticleDefinition *> TrackData::calcTypes() const {
    QMap<int32_t, const G4ParticleDefinition *> s;

    const int ntypes = 8;
    const G4ParticleDefinition *defs[ntypes] = {
        G4Gamma::Definition(),    G4Neutron::Definition(),
        G4Proton::Definition(),   G4Electron::Definition(),
        G4Positron::Definition(), G4OpticalPhoton::Definition(),
        G4Deuteron::Definition(), G4Triton::Definition()};

    const TrackPrivateData *tpd = data.constData();
    size_t ntracks = tpd->ntracks;
    TrackBlock *blocks = tpd->data;
    size_t i = 0;
    for (size_t z = 0; z < ntracks; z++) {
        int32_t type = blocks[i].h.ptype;
        i += blocks[i].h.npts + 1;

        if (!s.contains(type)) {
            const G4ParticleDefinition *found = nullptr;
            for (int j = 0; j < ntypes; j++) {
                if (type == defs[j]->GetPDGEncoding()) {
                    found = defs[j];
                    break;
                }
            }
            if (!found && type > 1000020000) {
                found = G4GenericIon::Definition();
            }
            if (found) {
                s[type] = found;
            } else {
                qWarning("Unidentified particle encoding %d", type);
            }
        }
    }

    return s;
}

static std::map<float, int> mapdedup(const std::vector<float> &a) {
    std::map<float, int> c;
    for (float b : a) {
        if (c.count(b)) {
            c[b] += 1;
        } else {
            c[b] = 1;
        }
    }
    return c;
}

static TrackMetaData *setupBallRadii(const TrackBlock *blocks, size_t ntracks) {
    TrackMetaData *meta = new TrackMetaData[ntracks];
    size_t i = 0;

    for (size_t z = 0; z < ntracks; z++) {
        const TrackHeader &h = blocks[i].h;
        const TrackBlock *pts = &blocks[i + 1];
        i += h.npts + 1;

        G4ThreeVector p0(pts[0].p.x, pts[0].p.y, pts[0].p.z);
        G4double radius = 0.0;
        for (int32_t j = 0; j < h.npts; j++) {
            const TrackPoint &tj = pts[j].p;
            G4ThreeVector pj(tj.x, tj.y, tj.z);
            radius = std::max(radius, (pj - p0).mag());
        }
        meta[z].ballRadius = radius;
    }
    return meta;
}

static void constructRangeHistogram(const std::vector<float> &starts,
                                    const std::vector<float> &ends,
                                    const std::vector<float> &spikes,
                                    QVector<QPointF> &pts) {
    std::map<float, int> mup = mapdedup(starts);
    std::map<float, int> mdown = mapdedup(ends);
    std::map<float, int> mdelt = mapdedup(spikes);
    mup[kInfinity] = 0;
    mdown[kInfinity] = 0;
    mdelt[kInfinity] = 0;
    std::map<float, int>::iterator upiter = mup.begin();
    std::map<float, int>::iterator downiter = mdown.begin();
    std::map<float, int>::iterator deltaiter = mdelt.begin();
    int height = 0;
    while (upiter != mup.end() || downiter != mdown.end() ||
           deltaiter != mdelt.end()) {
        // Pick next from all three; add least point
        std::pair<float, int> up = *upiter;
        std::pair<float, int> down = *downiter;
        std::pair<float, int> delta = *deltaiter;
        // Casework by number of ties for first
        float pos = std::min(up.first, std::min(down.first, delta.first));
        if (pos >= kInfinity) {
            break;
        }

        pts.push_back(QPointF(pos, height));
        if (up.first == pos) {
            height += up.second;
            ++upiter;
        }
        if (delta.first == pos) {
            pts.push_back(QPointF(pos, height + delta.second));
            ++deltaiter;
        }
        if (down.first == pos) {
            height -= down.second;
            ++downiter;
        }
        pts.push_back(QPointF(pos, height));
    }
}

void TrackData::constructRangeHistograms(QVector<QPointF> &tp,
                                         QVector<QPointF> &ep, const Range &tr,
                                         const Range &er) const {
    TrackBlock *blocks = data.constData()->data;
    size_t ntracks = data.constData()->ntracks;

    std::vector<float> tstarts, tends, tspikes;
    std::vector<float> estarts, eends, espikes;
    size_t i = 0;
    for (size_t m = 0; m < ntracks; m++) {
        const TrackHeader &h = blocks[i].h;
        i++;
        for (int32_t j = 0; j < h.npts - 1; j++) {
            const TrackPoint &ptA = blocks[i].p;
            const TrackPoint &ptB = blocks[i + 1].p;
            i++;

            float ta = ptA.time, tb = ptB.time;
            float ea = ptA.energy, eb = ptB.energy;

            float stl = std::max(0.f, (float(tr.low) - ta) /
                                          (ta == tb ? 1e38f : tb - ta));
            float sth = std::min(1.f, (float(tr.high) - ta) /
                                          (ta == tb ? 1e38f : tb - ta));
            float sel = std::max(0.f, (float(er.low) - ea) /
                                          (ea == eb ? 1e38f : eb - ea));
            float seh = std::min(1.f, (float(er.high) - ea) /
                                          (ea == eb ? 1e38f : eb - ea));
            float cta, ctb, cea, ceb;
            // restrict one by the other and vica versa
            cta = sel * ta + (1 - sel) * tb;
            ctb = (1 - seh) * ta + seh * tb;
            cea = stl * ea + (1 - stl) * eb;
            ceb = (1 - sth) * ea + sth * eb;
            //            cta = ta;
            //            ctb = tb;
            //            cea = ea;
            //            ceb = eb;
            if (cta == ctb) {
                //                tspikes.push_back(cta/CLHEP::ns);
            } else {
                tstarts.push_back(std::min(cta, ctb) / CLHEP::ns);
                tends.push_back(std::max(cta, ctb) / CLHEP::ns);
            }
            if (cea == ceb) {
                //                espikes.push_back(ea/CLHEP::eV);
            } else {
                estarts.push_back(std::min(cea, ceb) / CLHEP::eV);
                eends.push_back(std::max(cea, ceb) / CLHEP::eV);
            }
        }
        i++;
    }
    constructRangeHistogram(estarts, eends, espikes, ep);
    constructRangeHistogram(tstarts, tends, tspikes, tp);
}

TrackPrivateData::TrackPrivateData(const char *filename) {
    struct stat sb;
    int fd = open(filename, O_RDONLY);
    if (fd == -1) {
        qFatal("Invalid track file, '%s'", filename);
    }
    fstat(fd, &sb);

    mmapbytes = sb.st_size;
    char *buf = (char *)mmap(NULL, mmapbytes, PROT_READ, MAP_SHARED, fd, 0);
    close(fd);

    static_assert(sizeof(TrackPoint) == sizeof(TrackHeader) &&
                      sizeof(TrackPoint) == 32,
                  "Need uniform chunk size");
    nblocks = mmapbytes / sizeof(TrackPoint);

    data = reinterpret_cast<TrackBlock *>(buf);
    ntracks = 0;
    for (size_t i = 0; i < nblocks;) {
        i += data[i].h.npts + 1;
        ntracks++;
    }

    meta = setupBallRadii(data, ntracks);
}

TrackPrivateData::TrackPrivateData(size_t itracks, size_t iblocks,
                                   TrackBlock *idata) {
    ntracks = itracks;
    nblocks = iblocks;
    mmapbytes = 0;
    data = idata;
    meta = setupBallRadii(data, ntracks);
}
TrackPrivateData::TrackPrivateData(const TrackPrivateData &other)
    : QSharedData(other), ntracks(other.ntracks), nblocks(other.nblocks),
      data(other.data), meta(other.meta), mmapbytes(other.mmapbytes) {
    ref.ref();
}
TrackPrivateData::~TrackPrivateData() {
    if (mmapbytes > 0) {
        munmap(data, mmapbytes);
    } else {
        free(data);
    }
    if (meta) {
        delete[] meta;
    }
}
