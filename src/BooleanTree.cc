/* SPDX-License-Identifier: GPL-3.0-only */
#include "BooleanTree.hh"

#include <G4Box.hh>
#include <G4DisplacedSolid.hh>
#include <G4IntersectionSolid.hh>
#include <G4ScaleTransform.hh>
#include <G4ScaledSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>

#include <QDebug>

const G4VSolid *BooleanTree::compile(const G4VSolid *input) {
    if (dynamic_cast<const G4UnionSolid *>(input) ||
        dynamic_cast<const G4SubtractionSolid *>(input) ||
        dynamic_cast<const G4IntersectionSolid *>(input) ||
        dynamic_cast<const G4DisplacedSolid *>(input) ||
        dynamic_cast<const G4ScaledSolid *>(input)) {
        return new BooleanTree(input);
    } else {
        // We can't optimize this
        return input;
    }
}

static int rootCount(const G4VSolid *s) {
    // Multi
    const G4BooleanSolid *as_boolean = dynamic_cast<const G4BooleanSolid *>(s);
    const G4DisplacedSolid *as_displaced =
        dynamic_cast<const G4DisplacedSolid *>(s);
    const G4ScaledSolid *as_scaled = dynamic_cast<const G4ScaledSolid *>(s);
    if (as_boolean) {
        return rootCount(as_boolean->GetConstituentSolid(0)) +
               rootCount(as_boolean->GetConstituentSolid(1));
    } else if (as_displaced) {
        return rootCount(as_displaced->GetConstituentMovedSolid());
    } else if (as_scaled) {
        return rootCount(as_scaled->GetUnscaledSolid());
    } else {
        return 1;
    }
}

static int rectreefill(const G4VSolid *cur, BooleanTreeRoot *roots,
                       BooleanTreeNode *nodes, int *rootlen, int *nodelen,
                       const G4Transform3D &progtrans) {
    const G4UnionSolid *as_union = dynamic_cast<const G4UnionSolid *>(cur);
    const G4IntersectionSolid *as_inter =
        dynamic_cast<const G4IntersectionSolid *>(cur);
    const G4SubtractionSolid *as_sub =
        dynamic_cast<const G4SubtractionSolid *>(cur);
    const G4DisplacedSolid *as_displaced =
        dynamic_cast<const G4DisplacedSolid *>(cur);
    const G4ScaledSolid *as_scaled = dynamic_cast<const G4ScaledSolid *>(cur);
    if (as_union || as_inter || as_sub) {
        int i = *nodelen;
        (*nodelen)++;

        const G4BooleanSolid *as_boolean =
            dynamic_cast<const G4BooleanSolid *>(cur);
        const G4VSolid *va = as_boolean->GetConstituentSolid(0);
        const G4VSolid *vb = as_boolean->GetConstituentSolid(1);
        int ka = rectreefill(va, roots, nodes, rootlen, nodelen, progtrans);
        int kb = rectreefill(vb, roots, nodes, rootlen, nodelen, progtrans);

        nodes[i].kidA = ka;
        nodes[i].kidB = kb;
        nodes[i].debug = as_boolean;

        if (as_union) {
            nodes[i].op = OpUnion;
        } else if (as_inter) {
            nodes[i].op = OpIntersection;
        } else {
            nodes[i].op = OpDifference;
        }
        return -i - 1;
    } else if (as_scaled || as_displaced) {
        G4Transform3D ntrans;
        const G4VSolid *child;
        if (as_scaled) {
            ntrans = progtrans * as_scaled->GetScaleTransform();
            child = as_scaled->GetUnscaledSolid();
        } else {
            ntrans = progtrans * (G4Transform3D)as_displaced->GetTransform();
            child = as_displaced->GetConstituentMovedSolid();
        }
        return rectreefill(child, roots, nodes, rootlen, nodelen, ntrans);
    } else {
        int i = *rootlen;
        (*rootlen)++;
        roots[i].root = cur;
        roots[i].trans_data = progtrans;

        const G4Box *as_box = dynamic_cast<const G4Box *>(cur);
        const G4Tubs *as_tubs = dynamic_cast<const G4Tubs *>(cur);
        if (as_box) {
            roots[i].type = RootBox;
        } else if (as_tubs) {
            roots[i].type = RootTubs;
        } else {
            roots[i].type = RootOther;
        }
        HepGeom::Scale3D scale;
        HepGeom::Rotate3D rotation;
        HepGeom::Translate3D translation;
        roots[i].trans_data.getDecomposition(scale, rotation, translation);

        G4ThreeVector sc(scale.xx(), scale.yy(), scale.zz());
        G4double sce2 = (sc - G4ThreeVector(1., 1., 1.)).mag();
        G4double re =
            rotation.getRotation().distance2(CLHEP::HepRotation::IDENTITY);
        G4double te2 = translation.getTranslation().mag();
        bool ntrivrot = re > 0. || sce2 > 0.;
        bool ntrivtra = te2 > 0.;
        if (ntrivrot) {
            roots[i].trans_type = TransformFull;
        } else if (ntrivtra) {
            roots[i].trans_type = TransformOffset;
        } else {
            roots[i].trans_type = TransformNone;
        }
        return i + 1;
    }
}

BooleanTree::BooleanTree(const G4VSolid *s)
    : G4VSolid(s->GetName()), original(s), N(rootCount(s)) {
    roots = new BooleanTreeRoot[N];
    nodes = new BooleanTreeNode[N - 1];
    int nroots = 0;
    int nnodes = 0;
    G4Translate3D trans;
    start = rectreefill(s, roots, nodes, &nroots, &nnodes, trans);
}
BooleanTree::~BooleanTree() {
    delete[] roots;
    delete[] nodes;
}

G4bool BooleanTree::CalculateExtent(const EAxis pAxis,
                                    const G4VoxelLimits &pVoxelLimit,
                                    const G4AffineTransform &pTransform,
                                    G4double &pMin, G4double &pMax) const {
    return original->CalculateExtent(pAxis, pVoxelLimit, pTransform, pMin,
                                     pMax);
}
EInside BooleanTree::Inside(const G4ThreeVector &p) const {
    return original->Inside(p);
}
G4ThreeVector BooleanTree::SurfaceNormal(const G4ThreeVector &p) const {
    return original->SurfaceNormal(p);
}

typedef enum { rAbove, rBelow, rInside } RangeState;

typedef struct {
    double v;
    RangeState s;
} Query;

/*
 * Compute the query distance D.
 */
static Query recDTI(int id, BooleanTreeRoot *roots, BooleanTreeNode *nodes,
                    const G4ThreeVector &p, const G4ThreeVector &v,
                    double lowThresh = -kInfinity,
                    double highThresh = kInfinity) {

    // We return true iff the final distance D is lowThresh<= d <=highThresh
    // If we return false, Q: tooHigh, tooLow -- distinct cases matter!

    // Note: lowThresh can be explicitly tested for. Also, we can skip forward
    // via an Inside() then DTO/DTI call, which *may* be cheaper than a DTI
    // +retry, and definitely is cheaper than DTO followed by DTI, if we have a
    // prior leading edge distance. Hence may as well always skip forward (?)

    if (id < 0) {
        id = -(id + 1);

        const BooleanTreeNode &t = nodes[id];
        switch (t.op) {
        case OpUnion: {
            Query qA =
                recDTI(t.kidA, roots, nodes, p, v, lowThresh, highThresh);
            if (qA.s == rBelow) {
                return qA;
            } else if (qA.s == rAbove) {
                return recDTI(t.kidB, roots, nodes, p, v, lowThresh,
                              highThresh);
            } else {
                Query qB = recDTI(t.kidB, roots, nodes, p, v, lowThresh, qA.v);
                if (qB.s == rInside) {
                    Query r;
                    r.s = rInside;
                    r.v = std::min(qB.v, qA.v);
                    return r;
                } else if (qB.s == rAbove) {
                    return qA;
                } else {
                    return qB;
                }
            }
        } break;
        case OpIntersection: {
            // TODO: expand alt search
            Query q;
            const G4IntersectionSolid *s = (const G4IntersectionSolid *)t.debug;
            q.v = s->G4IntersectionSolid::DistanceToIn(p, v);
            q.s = q.v < lowThresh ? rBelow
                                  : (q.v > highThresh ? rAbove : rInside);
            return q;
        } break;
        case OpDifference: {
            // TODO: expand
            const G4SubtractionSolid *s = (const G4SubtractionSolid *)t.debug;
            Query q;
            q.v = s->G4SubtractionSolid::DistanceToIn(p, v);
            q.s = q.v < lowThresh ? rBelow
                                  : (q.v > highThresh ? rAbove : rInside);
            return q;
        } break;
        default:
            qFatal("Bad case");
            return Query();
            break;
        }
    } else {
        id = id - 1;
        const BooleanTreeRoot &t = roots[id];
        double dist = 0.;
        if (t.trans_type != TransformNone) {
            // TransformPoint; TransformAxis
            const G4ThreeVector &newPoint =
                t.trans_data.getRotation() * p + t.trans_data.getTranslation();
            const G4ThreeVector &newDirection = t.trans_data.getRotation() * v;
            dist = t.root->DistanceToIn(newPoint, newDirection);
        } else {
            dist = t.root->DistanceToIn(p, v);
        }
        Query q;
        q.v = dist;
        q.s =
            dist < lowThresh ? rBelow : (dist > highThresh ? rAbove : rInside);
        return q;
    }
}

G4double BooleanTree::DistanceToIn(const G4ThreeVector &p,
                                   const G4ThreeVector &v) const {
    // The rare simple translation/scale case
    if (N <= 1)
        return original->DistanceToIn(p, v);

    // We use region searching to prune unnecessary calls
    Query q = recDTI(start, roots, nodes, p, v);
    G4double distance = q.v;

    G4double comp = original->DistanceToIn(p, v);
    if (((comp - distance) / (comp + distance)) > 1e-10) {
        qDebug("EEK %g %g", comp, distance);
    }
    return comp;
}
G4double BooleanTree::DistanceToIn(const G4ThreeVector &p) const {
    return original->DistanceToIn(p);
}
G4double BooleanTree::DistanceToOut(const G4ThreeVector &p,
                                    const G4ThreeVector &v,
                                    const G4bool calcNorm, G4bool *validNorm,
                                    G4ThreeVector *n) const {
    // TODO: heavily used/replace

    return original->DistanceToOut(p, v, calcNorm, validNorm, n);
}
G4double BooleanTree::DistanceToOut(const G4ThreeVector &p) const {
    return original->DistanceToOut(p);
}
G4GeometryType BooleanTree::GetEntityType() const {
    return original->GetEntityType();
}
std::ostream &BooleanTree::StreamInfo(std::ostream &os) const {
    return original->StreamInfo(os);
}
void BooleanTree::DescribeYourselfTo(G4VGraphicsScene &scene) const {
    original->DescribeYourselfTo(scene);
}
