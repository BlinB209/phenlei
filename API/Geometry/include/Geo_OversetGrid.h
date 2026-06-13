//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      Geo_OversetGrid.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include "Geo_Grid.h"
#include "Geo_Interface.h"
#include "Geo_Interpoint.h"
#pragma warning (disable:4100)
// Composite Overlapping Meshes for the Solution of Partial Differential Equations
// G.CHESSGIRE AND W. D. HENSHAW
// JOURNAL OF COMPUTATIONAL PHYSICS 90, 1-64(1990)

// The algorithm for the generation of a composite grid must detect regions of overlap
// and determine the points to be used for interpolation between component grids.
// In addition, the algorithm should recognize regions of grid points which are not 
// needed in the computation.

// A valid composite grid G consists of a set of component grids Gk, labelled by k, 
// k = 1, ..., ng. The particular composite grid which is generated will depend upon a 
// number of parameters including
// (1) the width of the interpolation formula,
// (2) the width of the discretization formula,
// (3) the minimum allowable overlap,
// (4) the ordering of the component grids.

// Each point(i,j,k) on a valid composite grid must be one of the following:
// (1) Discretization point.
//     A discretization point is either an interior point or a boundary point.
//     Point(i,j,k) is called an interior point if it can be discretized, to
//     the required order, in terms of points on component grid Gk which are
//     interior points, boundary points, or interpolation points. The discretization
//     is assumed to be centered so that a rectangle of points is required;
//     the width and height of the rectangle being specified by the user.
//     Point(i,j,k) is called a boundary point if it lies on the true boundary
//     and can be discretized to the required order in terms of points on grid Gk
//     which are interior points, boundary points, or interpolation points.
//     (Boundary condition discretizations may be different from interior point
//     discretizations.)
// (2) Interpolation point.
//     Point(i,j,k) is an interpolation point if it can be interpolated from 
//     discretization or interpolation points on another component grid Gk'
//     with k'<> k to the required order.
// (3) Exterior or unused point.
//     Point(i,j,k) is an exterior or unused point if it is not a discretization
//     or interpolation point.

// It is possible that a given point could belong to more than one of the above
// categories. For example, some discretization points might just as well be 
// interpolation points. However, for efficiency we try to create a composite grid 
// with a minimum number of interpolation points. It is possible that the only valid
// composits grid is one consisting entirely of unused points, a null grid. This 
// situtation can be avoided by having enough overlap between component grids.

// The basic idea in our composite grid construction algorithm is to think of 
// ordering the component grids, k = 1, 2, ..., ng so that higher-numbered grids 
// cover over parts of lower-numbered component grids. Points which are removed will
// be in general lie underneath a higher numbered component grid.

// We now outline the composite grid algorithm for constructing a composite grid from
// a set of component grids. We combine an explanation in words with brief sections of
// pseudo-code. The composite grid is described by a flag array, kr(i,j,k)(where kr 
// stands for coordinates). This array contains a code for each point on the composite grid, 
// indicating whether the point(i,j,k) is a discretization, exterior of interpolation point.
// By the end of the algorithm the values in this array will have the following meaning

//            = k     if (i,j,k) is a discretization point
//  kr(i,j,k) = - k'  if (i,j,k) is interpolated from grid k' <> k
//            = 0     if (i,j,k) is an exterior point

namespace PHSPACE
{
class Grid;
class StructGrid;
class UnstructGrid;

enum OversetPointType
{ 
    ExteriorPoint = 0
};

class Geo_OversetGrid
{
private:
    Grid *grid;
    RDouble *xb, *yb, *zb;
    uint_t numberOfBoundaryPoints;
    int rank;
    int nBlocks;

    OversetInfoProxy *oversetInfoProxy;
public:
    Geo_OversetGrid(Grid *grid_in, int rank, int nBlocks);
    ~Geo_OversetGrid();
public:
    bool IsExteriorPoint(int flagOfPoint);

    int GetZoneID() const { return GetGrid()->GetZoneID(); }

    OversetInfoProxy * GetOversetInfoProxy() { return oversetInfoProxy; }

    Grid * GetGrid() const { return grid; }
    void SetX(RDouble *xb) { this->xb = xb; }
    void SetY(RDouble *yb) { this->yb = yb; }
    void SetZ(RDouble *zb) { this->zb = zb; }
    RDouble * GetX() const { return xb; }
    RDouble * GetY() const { return yb; }
    RDouble * GetZ() const { return zb; }

    void SetNumberOfBoundaryPoints(uint_t numberOfBoundaryPoints) { this->numberOfBoundaryPoints = numberOfBoundaryPoints; }
    uint_t GetNumberOfBoundaryPoints() const { return numberOfBoundaryPoints; }

    void SetRank(int rank) { this->rank = rank; }
    int  GetRank() const { return rank; }
    int  GetNBlocks() const { return nBlocks; }
    void InitReceivingInformation(vector<int> &cellIndexRef1, vector<int> &zoneIndexRef2, vector<int> &cellIndexRef2);
    void PrepareSendingInformation(Geo_OversetGrid **oversetGrids);
    void DumpOversetGrid(fstream &file);
public:
    virtual void InitialFlagOfGridPoints(int numberOfOversetGrid) {};
    virtual void MarkNonBoundaryPoints(int numberOfBoundaryPoints, RDouble *xBoundary, RDouble *yBoundary, RDouble *zBoundary) {};
    virtual bool FindGridOfEachPointCanBeInterpolatedFrom(Geo_OversetGrid **oversetGrids) { return true; }
    virtual void MarkAllPointsOnLowerGridNeededByHigherGrid(Geo_OversetGrid **oversetGrids) {};
    virtual void MarkAllPointsOnHigherGridNeededByLowerGrid(Geo_OversetGrid **oversetGrids) {};
    virtual void NormalizeAllValidPoints(Geo_OversetGrid **oversetGrids) {};
public:
    virtual void PrepareReceivingInformation(Geo_OversetGrid **oversetGrids) {};
    virtual void MarkNegativeSignOfIndexFlag(int index) {};
    virtual bool FindInterplatePoints(RDouble xm, RDouble ym, RDouble zm, int &index, RDouble &xit, RDouble &yit, RDouble &zit) = 0;
    virtual void SetFlag(vector<int> &idlist, int value) {};
    virtual void DumpFlags(fstream &file) {};
    virtual void PreProcess();
    virtual void PostProcess(Geo_OversetGrid **oversetGrids);
    virtual void GetCellCenterCoor(vector<int> &index, vector<RDouble> &xcc, vector<RDouble> &ycc, vector<RDouble> &zcc) {};
};

bool p_On_RHSide(RDouble &x, RDouble &y, RDouble &z, RDouble *p1, RDouble *p2);
bool InBox(RDouble &xb, RDouble &yb, RDouble &zb, RDouble *pmin, RDouble *pmax);
bool RightSide(RDouble &xm, RDouble &ym, RDouble &zm, RDouble &xfc, RDouble &yfc, RDouble &zfc, RDouble &xfn, RDouble &yfn, RDouble &zfn);

vector<string> GetTecplotTitle(const string &field_name);
void DumpField(UnstructGrid *grid, fstream &file, RFloat *field, const string &name);
void Visualization(UnstructGrid *grid, RFloat *field, const string &filename, const string &fieldname);

}
