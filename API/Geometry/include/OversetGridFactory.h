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
//! @file      OversetGridFactory.h
//! @brief     Probe overlap cells.
//! @author    Guo Yongheng.

#pragma once
#include "PHMpi.h"
#include "Geo_StructGrid.h"

namespace PHSPACE
{
const int IN_HOLE_COLOR = 0;
const int GENERIC_COLOR = 1;
const int OVERSET_COLOR = 2;
const int SECOND_COLOR  = 3;

class HoleGrid
{
protected:
    int zone_st, zone_ed, zone_at;
    int ns, np;
    vector< int > *nsi, *nsj, *nsp;
    vector< RDouble > *x, *y, *z;
    vector< RDouble > *box;
    vector< vector< RDouble > * > *capsule;
public:
    LIB_EXPORT  HoleGrid();
    LIB_EXPORT ~HoleGrid();

    void SetZoneStart      (int iZone);
    void SetZoneEnd        (int iZone);
    void SetZoneAttached   (int iZone);
    void SetNumberOfSurface(int nSurf);
    void SetNumberOfPoints (int np);
    
    int GetZoneStart();
    int GetZoneEnd();
    int GetZoneAttached();
    int GetNumberOfSurface();
    int GetNumberOfPoints();

    int * GetNSI();
    int * GetNSJ();
    int * GetNSP();

    RDouble * GetX();
    RDouble * GetY();
    RDouble * GetZ();

    RDouble * GetBox();

    LIB_EXPORT void DigHoles();
    LIB_EXPORT void GenerateSurfaceSpace();
    LIB_EXPORT void GeneratePointSpace();
    LIB_EXPORT int GetCellColor(RDouble x, RDouble y, RDouble z, RDouble far);
private:
    void GenerateCapsualSpace();
    void ProbeBorderCell(int iType, int jType);
    void MeasureBoxAndCapsual();
    void ProbeCellsInHole();
    void CopyNodeCoordinate(vector< vector< RDouble > > &r, vector< vector< RDouble > > &s, vector< int > &node);
    void ProduceCasing(vector< vector< RDouble > > &nodeCoordinate, int i);
    void Radial_Triangle_Intersect(vector< RDouble > &pst, vector< RDouble > &ped, vector< RDouble > &xqd, vector< RDouble > &yqd, vector< RDouble > &zqd, vector< RDouble > &pp_ist, int &pCounter);
    void GaussElimination(vector< vector< RDouble > > &a, vector< RDouble > &b, int n, vector< RDouble > &x, int &flag);
private:
    vector< vector< RDouble > > phy;
};

class Stainer
{
private:
    int numberOfHoles;
    vector< HoleGrid* > *holeGridContainer;
public:
    LIB_EXPORT  Stainer();
    LIB_EXPORT ~Stainer();

    LIB_EXPORT void ReadHoleData();
    LIB_EXPORT void ProbeOversetCells();
private:
    void PaintGenericColor();
    void ProbeInternalOverlapCells();
    void ProbeExternalOverlapCells();
    void EraseSecondColor();
};

class BackgroundTree 
{
private:
    bool isEmptyBinaryTree;
    int zone_st, zone_ed, numberOfUnitCells;
    vector< RDouble > *unitCellCapsule, *inf_plane, *sup_plane;
    int *zone_location, *i_location, *j_location, *k_location, *fruitLevel;
    vector< int > *fruitInheritMapping;
    vector< int > fruitIndexStorage;
public:
    LIB_EXPORT  BackgroundTree();
    LIB_EXPORT ~BackgroundTree();

    LIB_EXPORT void Generate(int zone_st, int zone_ed);
    LIB_EXPORT void SetUnitCellParameter(int iUnitCell, vector< RDouble > &xx, vector< RDouble > &yy, vector< RDouble > &zz);
    LIB_EXPORT void SearchNodes(OversetCellCollector *oversetCellCollector, RDouble xx, RDouble yy, RDouble zz, int ip, int iElement);
    LIB_EXPORT void ProbeNearestPoint(RDouble xx, RDouble yy, RDouble zz, RDouble &distance, int &index);

    void SetZoneViewRange(int st, int ed);
    void SetUnitCellLocation(int iUnitCell, int iZone, int iCore, int jCore, int kCore);
    bool IsEmpty();
private:
    void ComputeNumberOfUnitCells();
    void GenerateSourceSpace();
    void ComputeUnitCellCapsule();
    void InitFruitLevelAndInheritMapping();
    void ComputeRootSpaceParameter();
    void AppendNewFruit();
    void AppendNewFruit(int iUnitCell);
        
    void SearchNodes(RDouble xx, RDouble yy, RDouble zz, int n);
    void ComputeLocalCoordinateByNewtonMethod(RDouble xx, RDouble yy, RDouble zz, vector< RDouble > &ss, vector< int > &node);

    bool UnitCellCapsuleBelongsToNewNode(vector< RDouble > &capsule, vector< RDouble > &inf_plane, vector< RDouble > &sup_plane);
    bool PointBelongsToUnitCellCapsule(RDouble xx, RDouble yy, RDouble zz, vector< RDouble > &capsule);
    bool PointBelongsToNewNode(RDouble xx, RDouble yy, RDouble zz, int kNode);

    void ProduceNewLeftNode (int kNode, vector< RDouble > &inf_left,  vector< RDouble > &sup_left);
    void ProduceNewRightNode(int kNode, vector< RDouble > &inf_right, vector< RDouble > &sup_right);

    void AppendNewLeftNode (int iUnitCell, int kNode, vector< RDouble > &inf_new_node, vector< RDouble > &sup_new_node);
    void AppendNewRightNode(int iUnitCell, int kNode, vector< RDouble > &inf_new_node, vector< RDouble > &sup_new_node);
};

class OversetCellCollector 
{
private:
    int zone_st, zone_ed;
    vector< int > elementNumberRecorder;
    vector< int > sourceProcessorIndex;

    vector< int > *localCellIndex, *zoneIndex;
    vector< RDouble > *xCenter, *yCenter, *zCenter;

    vector< int > *spaceLattice;
    vector< RDouble > *xLocal, *yLocal, *zLocal, *radius;
public:
    LIB_EXPORT  OversetCellCollector();
    LIB_EXPORT ~OversetCellCollector();

    LIB_EXPORT void RunGeometricAnalysis(int zone_st, int zone_ed);
    LIB_EXPORT void SearchOptimalUnitCell(BackgroundTree *backgroundTree);
    LIB_EXPORT void ResizeSourceIndexContainer();
    LIB_EXPORT void InsertOversetParameter(vector< int > *zoneContainer, vector< int > *cellContainer);
    LIB_EXPORT void InsertServiceParameter(vector< RDouble > *ksai, vector< int > *waiter);
    LIB_EXPORT void CollectOversetCells(int iZone, vector< int > &cell, vector< RDouble > &xx, vector< RDouble > &yy, vector< RDouble > &zz);
    LIB_EXPORT void SetNearestPointParameter(int ip, int iElement, RDouble distance, int index);
    LIB_EXPORT void SetOptimalPointParameter(int ip, int iElement, vector< RDouble > &ss, vector< int > &node);
    
    void SetZoneViewRange(int st, int ed);
    
private:
    void GenerateTargetSpace();
    void CollectOversetCells();
    void BcastOversetCells();
    void GenerateServiceSpace();
    void SetVirtualPointParameter();
    void ProbeNearestPoint(BackgroundTree *backgroundTree, RDouble xx, RDouble yy, RDouble zz, int ip, int iElement);
};

class OversetStructGrid
{
private:
    int searchTimes;
    vector< int > sourceZoneStart, sourceZoneEnd, targetZoneStart, targetZoneEnd;
    vector< int > zoneInverseIndex;

    vector< RDouble > *ksai;
    vector< int > *waiter;
    vector< int > indexGroup, processGroup;

    vector< int > *zoneContainer, *cellContainer;
    vector< int > *nodeSend, *nodeRecv;
    
    vector< string > solverNameContainer;
    vector< int > equationNumberContainer;
    vector< RDouble4D * > *field;
    vector< RDouble > **dataSend, **dataRecv;
public:
    LIB_EXPORT  OversetStructGrid();
    LIB_EXPORT ~OversetStructGrid();

    LIB_EXPORT void RunGeometricAnalysis();

    LIB_EXPORT void FillSolverContainer();
    LIB_EXPORT void InitFieldNew();
    LIB_EXPORT void InitField();

    LIB_EXPORT void CommunicateOversetField(int iSolver);
    LIB_EXPORT void CorrectOversetField(int iSolver);

    int GetZoneInverseIndex(int i_view);
private:
    void ReadBasicParameter();
    void ProbeOversetCells();
    void GenerateSourceAndTargetSpace();
    void ProbeRelativeOversetLocation();
    void FillNodeIndexContainer();
    void InsertOversetAndServiceParameter(OversetCellCollector *oversetCellCollector);
    void LinearInsertValue(int numberOfEquations, vector< RDouble > &backdrop, RDouble xx, RDouble yy, RDouble zz, vector< RDouble > &newCellValue);
};

#include "OversetGridFactory.hxx"

int GetCodeOfDigHoles();
int GetCodeOfTurbulenceModel();

LIB_EXPORT void GenerateStructGridTable();
LIB_EXPORT void DeleteStructGridTable();
LIB_EXPORT StructGrid * GetStructGrid(int i);

LIB_EXPORT OversetStructGrid * GetOversetStructGrid();
LIB_EXPORT void CreateOversetStructGrid();
LIB_EXPORT void DeleteOversetStructGrid();

LIB_EXPORT int GetZoneInverseIndex(int i_view);

LIB_EXPORT void ComputePhysicalCoordinate(vector< vector< RDouble > > &frame, RDouble x, RDouble y, RDouble z, vector< RDouble > &r);

}
