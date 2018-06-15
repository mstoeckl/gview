/* SPDX-License-Identifier: GPL-3.0-only */
#pragma once

#include <G4Transform3D.hh>
#include <G4VSolid.hh>

class G4BooleanSolid;

typedef enum { OpUnion, OpIntersection, OpDifference } BooleanOperation;

typedef struct {
    BooleanOperation op;
    int kidA;
    int kidB;
    const G4BooleanSolid *debug;
} BooleanTreeNode;

typedef enum {
    RootBox,
    RootTubs,
    RootOther,
} BooleanTreeRootTypes;

typedef enum {
    TransformNone,
    TransformOffset,
    TransformFull,
} TransformType;

typedef struct {
    const G4VSolid *root;
    G4Transform3D trans_data;
    TransformType trans_type;
    BooleanTreeRootTypes type;
} BooleanTreeRoot;

class BooleanTree : public G4VSolid {
public:
    static const G4VSolid *compile(const G4VSolid *input);

    virtual G4bool CalculateExtent(const EAxis pAxis,
                                   const G4VoxelLimits &pVoxelLimit,
                                   const G4AffineTransform &pTransform,
                                   G4double &pMin,
                                   G4double &pMax) const override;
    virtual EInside Inside(const G4ThreeVector &p) const override;
    virtual G4ThreeVector SurfaceNormal(const G4ThreeVector &p) const override;
    virtual G4double DistanceToIn(const G4ThreeVector &p,
                                  const G4ThreeVector &v) const override;
    virtual G4double DistanceToIn(const G4ThreeVector &p) const override;
    virtual G4double DistanceToOut(const G4ThreeVector &p,
                                   const G4ThreeVector &v,
                                   const G4bool calcNorm = false,
                                   G4bool *validNorm = 0,
                                   G4ThreeVector *n = 0) const override;
    virtual G4double DistanceToOut(const G4ThreeVector &p) const override;
    virtual G4GeometryType GetEntityType() const override;
    virtual std::ostream &StreamInfo(std::ostream &os) const override;
    virtual void DescribeYourselfTo(G4VGraphicsScene &scene) const override;

    const G4VSolid *GetOriginal() const { return original; }

private:
    BooleanTree(const G4VSolid *);
    ~BooleanTree();

    const G4VSolid *original;
    const int N;

    BooleanTreeRoot *roots;
    BooleanTreeNode *nodes;
    int start;
};
