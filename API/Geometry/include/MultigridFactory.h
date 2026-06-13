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
//! @file      MultigridFactory.h
//! @brief     Analyse relationship between multigrid.
//! @author    Guo Yongheng.

#pragma once
#include "PHMpi.h"
#include "Geo_Grid.h"
#include "Geo_StructGrid.h"
#include "OversetGridFactory.h"

namespace PHSPACE
{
class PHSurface
{
private:
    int ni, nj, nk;
    int ist, jst, kst;
    int np;
    int *data;
public:
     PHSurface(Range I, Range J, Range K);
    ~PHSurface();
public:
    int & operator()(int i, int j, int k);
    int * GetData();
    int   GetNumberOfFaces();
};

class StructFace
{
private:
    int ni, nj, nk;
    vector< int > numberOfFaceElements;
    Range IST, ISE, IED;
    Range JST, JSE, JED;
    Range KST, KSE, KED;
    StructGrid *structGrid;

    vector< PHSurface * > *neighborZoneIndexContainer;
    vector< PHSurface * > *neighborZoneLocalInterfaceIndexContainer;
    vector< PHSurface * > *localInterfaceIndexContainer;
    vector< PHSurface * > *localBoundaryFaceIndexContainer;

    vector< StructFace * > *structFaceGroup;
public:
    StructFace(StructGrid *structGrid);
    ~StructFace();
public:
    void PaintSelf();
    void PaintNeighbor();
    void SetStructFaceGroup(vector< StructFace * > *structFaceGroup) { this->structFaceGroup = structFaceGroup; }
    void CollectInterfaceInformation(InterfaceInfo *interfaceInformation);
    int  ComputeNumberOfInterfaces();
private:
    void Init();
    void GenerateSubspace(vector< PHSurface * > *micro);

    PHSurface & GetNeighborZoneIndexContainer(int m) { return * (*neighborZoneIndexContainer)[m]; }
    PHSurface & GetNeighborZoneLocalInterfaceIndexContainer(int m) { return * (*neighborZoneLocalInterfaceIndexContainer)[m]; }
    PHSurface & GetLocalInterfaceIndexContainer(int m) { return * (*localInterfaceIndexContainer)[m]; }
    PHSurface & GetLocalBoundaryFaceIndexContainer(int m) { return * (*localBoundaryFaceIndexContainer)[m]; }

    int * GetNeighborZoneIndex(int m) { return (*neighborZoneIndexContainer)[m]->GetData(); }
    int * GetNeighborZoneLocalInterfaceIndex(int m) { return (*neighborZoneLocalInterfaceIndexContainer)[m]->GetData(); }
    int * GetLocalInterfaceIndex(int m) { return (*localInterfaceIndexContainer)[m]->GetData(); }
    int * GetLocalBoundaryFaceIndex(int m) { return (*localBoundaryFaceIndexContainer)[m]->GetData(); }

    StructFace * GetStructFaceGroup(int m) { return (*structFaceGroup)[m]; }
};

class MultigridManager
{
public:
    MultigridManager(Grid **OrdinaryGrid, int isOverset, int globalBlockNumber);
    ~MultigridManager();
public:
    void Run();
protected:
    Grid **OrdinaryGrid;
    int isOverset;
public:
    void GenerateGlobalNodeLinkStructure();
    void GenerateGlobalCellLinkStructure();
    void ReadGridCoordinate();
    void GridExtend();
    void ComputeCellCenter();
    void FaceExtend();
    void PrismLineExtend();
    void CornerExtend();
    void GetCrossPoint(int mPatch, int nPatch, RDouble x, RDouble y, RDouble z, RDouble &xx, RDouble &yy, RDouble &zz, int &nFlap, int &nPointBy, int &nBlock, int &iPoint, int &jPoint, int &kPoint);
    void GenerateGlobalBlockSpace();
    void GenerateGlobalPointSpace();
    void GenerateGlobalCenterSpace();
    void GenerateLocalPointSpace();
    void ReadBoundaryCondition();
    void InitialNodeLinksData();
    void SetNewPatchModel();
    void ComputeContravariantComponent();
    void OutputLinkData();
    void GenerateGlobalFaceSpace();
    void GenerateLocalFaceSpace();
    void InitializeLocalFaceSpace(PHInt2D &myInt2D, int value);
    void GetOtherBlockDirections(int id1, int &id2, int &id3);
    string GetLinkFileName();
public:
    void GenerateLinkStructure();
    void GenerateHoleStructure();
    void ReadOriginalHoleFile();
    
    void ComputeLocalCoordinate();
    void OutputNewHoleFile();
protected:
    int numberOfHoles;
    PHInt1D holeZoneStart, holeZoneEnd, holeZoneAttached, numberOfHoleSurfaces, numberOfHolePoints, zoneInverseIndex;
    PHInt2D nsi, nsj, spst;
    PHDouble2D xsurf, ysurf, zsurf;
protected:
    int globalBlockNumber, globalPointNumber, globalCenterNumber, maxLocalFaceNumber, specialBlockLabel;
    PHInt1D overallFaceNumber, externalFaceNumber, internalFaceNumber, processLabel;
    PHInt1D startPointLabel, iDimension, jDimension, kDimension;
    PHInt1D patchIPoint, patchJPoint, patchKPoint;
    PHInt2D iminMaster, imaxMaster, jminMaster, jmaxMaster, kminMaster, kmaxMaster;
    PHInt2D iminSlave, imaxSlave, jminSlave, jmaxSlave, kminSlave, kmaxSlave;
    PHInt2D iminInner, imaxInner, jminInner, jmaxInner, kminInner, kmaxInner;
    PHInt2D tarPatch, faceType;
    PHInt1D startCenterLabel, iCenterDimension, jCenterDimension, kCenterDimension;
    PHInt1D blockLocation, iLocation, jLocation, kLocation;
    PHInt1D characterPoint;
    PHDouble1D coordinateX, coordinateY, coordinateZ;
    PHDouble1D cellCenterX, cellCenterY, cellCenterZ;

public:
    void GenerateHoleGridStructure();
protected:
    void ReadZoneInverseMapping();
    void ReadOriginalHoleGrid();
    void ComputeLocalHoleGrid();
    void OutputNewHoleGrid();
protected:
    vector<HoleGrid *> holeGridList;
};

void Standardize(vector< int > &st, vector< int > &ed, int side);
string GetHoleBasicFileName();
void GaussElimination(PHDouble2D &a, PHDouble1D &b, int n, PHDouble1D &x, int &flag);
}

