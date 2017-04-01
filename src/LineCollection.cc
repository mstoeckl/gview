#include "LineCollection.hh"

typedef struct {
    size_t center;
    size_t xpt;
    size_t ypt;
} Trindx;

class LCData {
public:
    // Same recursive construction method
    // as before; only counting rules are nasty..
    LCData() {
        left = NULL;
        right = NULL;
        up = NULL;
    }
    ~LCData() {
        if (left)
            delete left;
        if (right)
            delete right;
    }
    LCData *up;
    LCData *left;
    LCData *right;
    // Note: 1 size_t, 1 length
    // are eqvt.
    std::vector<std::pair<size_t, int>> contents;
    Trindx plane;
};

static void getCenterAndNormal(const LCData* node, const TrackPoint* points, G4ThreeVector& center, G4ThreeVector& normal) {
    const TrackPoint& cpt = points[node->plane.center];
    const TrackPoint& xpt = points[node->plane.xpt];
    const TrackPoint& ypt = points[node->plane.ypt];
    G4ThreeVector pc(cpt.x,cpt.y,cpt.z);
    G4ThreeVector px(xpt.x,xpt.y,xpt.z);
    G4ThreeVector py(ypt.x,ypt.y,ypt.z);
    normal = (px - pc).cross(py - pc);
    center = pc;
}

SimplexIterator::SimplexIterator(const LineCollection& iin, const G4ThreeVector &isrc,
                                 const G4ThreeVector &idir)
    : lc(iin), src(isrc), dir(idir) {
    current = lc.d;
    t = 0;
    // Descend tree to locate current element
    while (current->left) {
        G4ThreeVector center, normal;
        getCenterAndNormal(current, lc.mpoints, center, normal);
        if ((src - center)*normal > 0) {
            current = current->left;
        } else {
            current = current->right;
        }
    }
}

std::vector<std::pair<size_t, int>> SimplexIterator::getContainedLines() const {
    if (!current) {
        G4Exception("SimplexIterator::getContainedLines",
                    "Iterator left tree hierarchy", FatalException,
                    "Description");
    }
    return current->contents;
}
bool SimplexIterator::advance() {
    if (!current) {
        G4Exception("SimplexIterator::getContainedLines",
                    "Iterator left tree hierarchy", FatalException,
                    "Description");
    }
    // TODO: casework!!



    // Move to next node in tree.... | go up, fwd, etc...
    return false;
}

static size_t getMaxCount(LCData *d) {
    if (d->left) {
        return getMaxCount(d->left);
    }
    return d->contents.size();
}

static size_t randint(size_t exclmax) {
    // If we were to have 2**30+1 lines, needed
    size_t mx = (size_t(RAND_MAX) + 1);
    if (Q_UNLIKELY(exclmax >= RAND_MAX)) {
        mx = mx * mx;
        size_t cap = mx - mx % exclmax;
        size_t cc;
        do {
            cc = size_t(RAND_MAX) * size_t(random()) + size_t(random());
        } while (cc >= cap);
        return cc % exclmax;
    } else {
        size_t cap = mx - mx % exclmax;
        size_t cc;
        do {
            cc = size_t(random());
        } while (cc >= cap);
        return cc % exclmax;
    }
}

static bool splitNode(LCData *root, const TrackHeader *headers,
                      const TrackPoint *points) {
    // Construct divplane so as to divide lines among children
    // So, pick N sets of random triplets of points
    // from line collection. (Is: 3x choice.. N=1..2K)
    size_t np = root->contents.size();
    double x = 0., y = 0., z = 0.;
    for (const std::pair<size_t, int> &line : root->contents) {
        const TrackHeader &h = headers[line.first];
        const TrackPoint &a = points[h.offset + size_t(line.second)];
        const TrackPoint &b = points[h.offset + size_t(line.second) + 1];
        x += a.x + b.x;
        y += a.y + b.y;
        z += a.z + b.z;
    }
    G4ThreeVector centroid(0.5 * x / np, 0.5 * y / np, 0.5 * z / np);

    size_t besti = 0, bestj = 0, bestk = 0;
    size_t score = np;
    for (int tr = 0; tr < 20; tr++) {
        size_t i = randint(2 * np), j = randint(2 * np), k = randint(2 * np);
        if (i == k || i == j || j == k) {
            continue;
        }
        // Let I be center and be closest to average location
        size_t ni = i / 2, di = i % 2;
        size_t nj = j / 2, dj = j % 2;
        size_t nk = k / 2, dk = k % 2;
        size_t ci = headers[root->contents[ni].first].offset +
                    size_t(root->contents[ni].second) + di;
        size_t cj = headers[root->contents[nj].first].offset +
                    size_t(root->contents[nj].second) + dj;
        size_t ck = headers[root->contents[nk].first].offset +
                    size_t(root->contents[nk].second) + dk;
        const TrackPoint &ti = points[ci];
        const TrackPoint &tj = points[cj];
        const TrackPoint &tk = points[ck];
        G4ThreeVector pi(ti.x, ti.y, ti.z);
        G4ThreeVector pj(tj.x, tj.y, tj.z);
        G4ThreeVector pk(tk.x, tk.y, tk.z);
        if ((pi - centroid).mag2() > (pj - centroid).mag2()) {
            std::swap(pi, pj);
            std::swap(ci, cj);
        }
        if ((pi - centroid).mag2() > (pk - centroid).mag2()) {
            std::swap(pi, pk);
            std::swap(ci, ck);
        }
        G4ThreeVector normal = (pj - pi).cross(pk - pi);
        if (normal.mag2() <= 0.) {
            // Plane needs noncollinear points
            continue;
        }
        // Now PI is center point!
        size_t counts[3] = {0, 0, 0}; // R,C,L
        for (const std::pair<size_t, int> &line : root->contents) {
            const TrackHeader &h = headers[line.first];
            const TrackPoint &a = points[h.offset + size_t(line.second)];
            const TrackPoint &b = points[h.offset + size_t(line.second) + 1];
            G4ThreeVector pa(a.x, a.y, a.z);
            G4ThreeVector pb(b.x, b.y, b.z);
            G4double sgna = (pa - pi) * normal;
            G4double sgnb = (pb - pi) * normal;
            if (sgna == 0. && sgnb == 0.) {
                // Line over plane: Use a heuristic splitting function
                if (line.first % 2 == 0) {
                    ++counts[0];
                } else {
                    ++counts[2];
                }
            }
            if (sgna >= 0. && sgnb >= 0.) {
                // Entirely on left side
                ++counts[0];
            } else if (sgna <= 0. && sgnb <= 0.) {
                // Entirely on right side
                ++counts[2];
            } else {
                // Crossing
                ++counts[1];
            }
        }
        // Goal: minimize the maximum count
        size_t maxcap = counts[1] + std::max(counts[0], counts[2]);
        if (maxcap <= score) {
            score = maxcap;
            besti = ci;
            bestj = cj;
            bestk = ck;
        }
    }
    if (score >= np) {
        // Failed to find a dividing plane that minimizes
        // the number of lines crossing the plane.
        return false;
    }
    root->plane.center = besti;
    root->plane.xpt = bestj;
    root->plane.ypt = bestk;
    root->left = new LCData();
    root->right = new LCData();
    root->left->up = root;
    root->right->up = root;
    const TrackPoint &tc = points[root->plane.center];
    const TrackPoint &tx = points[root->plane.xpt];
    const TrackPoint &ty = points[root->plane.ypt];
    G4ThreeVector pc(tc.x, tc.y, tc.z);
    G4ThreeVector px(tx.x, tx.y, tx.z);
    G4ThreeVector py(ty.x, ty.y, ty.z);
    G4ThreeVector normal = (px - pc).cross(py - pc);

    // TODO: abuse the fact that we know exactly
    // how many children there will be, so as to
    // minimize allocation costs. Perhaps malloc(),malloc(),N ?
    for (const std::pair<size_t, int> &line : root->contents) {
        const TrackHeader &h = headers[line.first];
        const TrackPoint &a = points[h.offset + size_t(line.second)];
        const TrackPoint &b = points[h.offset + size_t(line.second) + 1];
        G4ThreeVector pa(a.x, a.y, a.z);
        G4ThreeVector pb(b.x, b.y, b.z);
        G4double sgna = (pa - pc) * normal;
        G4double sgnb = (pb - pc) * normal;
        if (sgna == 0. && sgnb == 0.) {
            // Line over plane: Use a heuristic splitting function
            if (line.first % 2 == 0) {
                root->left->contents.push_back(line);
            } else {
                root->right->contents.push_back(line);
            }
        } else if (sgna >= 0. && sgnb >= 0.) {
            // Entirely on left side
            root->left->contents.push_back(line);
        } else if (sgna <= 0. && sgnb <= 0.) {
            // Entirely on right side
            root->right->contents.push_back(line);
        } else {
            // Crossing
            root->left->contents.push_back(line);
            root->right->contents.push_back(line);
        }
    }
    // Wipe original contents
    root->contents = std::vector<std::pair<size_t, int>>();
    return true;
}

static size_t recursiveCount(LCData *r) {
    if (!r->left)
        return r->contents.size();
    return recursiveCount(r->left) + recursiveCount(r->right);
}
static size_t recursiveLeafMax(LCData *r) {
    if (!r->left)
        return r->contents.size();
    return std::max(recursiveLeafMax(r->left), recursiveLeafMax(r->right));
}

static void recursivePrintTree(LCData *root, int depth = 0) {
    char ws[256];
    int i = 0;
    for (; i < depth && i < 255; i++) {
        ws[i] = ' ';
    }
    ws[i] = '\0';
    if (!root->left) {
        qDebug("%s - %lu", ws, root->contents.size());
    } else {
        size_t kids = recursiveCount(root);
        size_t maxleaf = recursiveLeafMax(root);
        qDebug("%s===> %lu %lu", ws, kids, maxleaf);
        recursivePrintTree(root->left, depth + 1);
        recursivePrintTree(root->right, depth + 1);
    }
}

LineCollection::LineCollection(const TrackHeader *headers,
                               const TrackPoint *points, size_t nheaders) {
    mheaders = headers;
    mpoints = points;
    d = new LCData();
    for (size_t i = 0; i < nheaders; i++) {
        const TrackHeader &h = headers[i];
        for (int j = 0; j < h.npts - 1; j++) {
            // Skip stationary lines, although
            // there probably are none
            const TrackPoint &a = points[h.offset + size_t(j)];
            const TrackPoint &b = points[h.offset + size_t(j) + 1];
            if (a.x - b.x == 0. && a.y - b.y == 0. && a.z - b.z == 0.) {
                continue;
            }
            d->contents.push_back(std::pair<size_t, int>(i, j));
        }
    }
    const size_t lthresh = 50;
    size_t failcount = 0;
    for (int niter = 0; niter < 100000; niter++) {
        if (niter % 1000 == 0) {
            qDebug("ct: %lu (%g)", getMaxCount(d), failcount / double(niter));
        }
        // Locate leftmost leaf
        LCData *lcr = d;
        while (lcr->left)
            lcr = lcr->left;
        if (lcr->contents.size() <= lthresh) {
            // Goal achieved
            break;
        }
        // Divide node
        bool wassplit = splitNode(lcr, headers, points);
        if (!wassplit) {
            failcount++;
            continue;
        }
        // Reorient tree, swapping L/R as need be
        while (lcr) {
            if (getMaxCount(lcr->left) < getMaxCount(lcr->right)) {
                std::swap(lcr->left, lcr->right);
                std::swap(lcr->plane.xpt, lcr->plane.ypt);
            }
            lcr = lcr->up;
        }
    }

    if (0) {
        recursivePrintTree(d);
    }

    // Issue: above algorithm is kinda slow :-(
    // (basically, N^2log(N) if we were to reach 1-leaf leaves
}

LineCollection::~LineCollection() { delete d; }

SimplexIterator LineCollection::followRay(const G4ThreeVector &src,
                                          const G4ThreeVector &dir) const {
    return SimplexIterator(*this, src, dir);
}
