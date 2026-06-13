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
//! @file      UnstructuredOversetConfig.h
//! @brief     Unstructured grid overset configuration.
//! @author    Chang Xinhua, He Kun, Zhang Yong, Zhang Laiping.

#pragma once
#include "PHMpi.h"
#include "OversetKDTree.h"
using namespace std;
#define perp(u,v) ((u).x * (v).y - (u).y * (v).x) 
namespace PHSPACE
{
class OversetConfigFactoryMulti;
class OversetInformationProxy;
class UnstructGrid;
class OversetConfig;
class InterCells;
class ZoneSearchNodes;
class ZoneSearchLines;
class Pre_WalldistCompute;
class TimeSpan;

void RunOversetGridConfig();

void RunUnstructuredOversetGridConfig();    //! Perform an oversetConfig individually.

void InitializeOversetGridConfig();

void InitializeExplicitMethod();

void ReConfigUnstructOversetGrid();

void ConfigOversetSlip();

void ExplicitHoleCut();

void ImplicitHoleCut();

void OutUnstructOversetGridConfig();

void FreeUnstructOversetGridConfig();

OversetConfigFactoryMulti *GetOversetConfigFactoryMulti();

class OversetConfigFactoryMulti    //! Call OversetConfig as a whole, and complete the information communication between each zone.
{
public:
    OversetConfigFactoryMulti();
    ~OversetConfigFactoryMulti();

private:
    int numberOfBlocks;    //! Number of overlapping grids.
    int numberOfZones;     //! Total number of partitions in the grid.
    OversetConfig **oversetConfigMultiZone;    //! Corresponding to each grid partition.
    int geometricDimension;
    bool hasCalMinMax;
private:
    ZoneSearchNodes **zonesSearchNodes;    //! Set of all calculate points.
    ZoneSearchLines **zonesSearchLines;    //! Set of all calculate lines.
    InterCells **zonesInterpolateCells;    //! Set of all interpolation cells.
    TimeSpan *testTime;

private:
    vector< UnstructGrid * > innerGrid;

private:    //! This part of data is used to complete the attribute communication of the points on the interface.
    vector< int > numberOfGlobalNodes;    //! Number of calculate points in the whole grid.
    vector< vector< int > > keyActiveOfGlobalInterNodes;    //! Property of points on the interface.

    int *hasSymmetryBCInBlocks;
    int *numberOfInterfaceNodesInBlocks;
    int **keyOfInterfaceNodesInBlocks;
private:    //! Smallest box of each Block is used to improve the query efficiency.
    RDouble **pMinOfBlocks;
    RDouble **pMaxOfBlocks;

    RDouble **pMinOfBlocksByLayer;
    RDouble **pMaxOfBlocksByLayer;
public:
    void Initialize();
    void OversetExConfigMulti();
    void OversetImConfigMulti();
    void OversetSlipConfigMulti();
    void FreeKeyBasicActiveNodes();
    void BuildNodesMinDistanceInOtherBlocks();
    void BuildNodesMinDistanceInOtherBlocks2();
    void BuildLocalToGlobalInterfaceNodesMaps();

    void NodesOperator();
    void CellsOperator();
    void OutPutOversetFiles();
    void OutPutOversetInformation();
    void OutPutOversetVisualizations(const string &name);
public:
    RDouble **GetMinBoxOfBlocks()
    {
        return this->pMinOfBlocks;
    }
    RDouble **GetMaxBoxOfBlocks()
    {
        return this->pMaxOfBlocks;
    }
    RDouble **GetMinBoxOfBlocksByLayer()
    {
        return this->pMinOfBlocksByLayer;
    }
    RDouble **GetMaxBoxOfBlocksByLayer()
    {
        return this->pMaxOfBlocksByLayer;
    }
    InterCells **GetZonesInterpolateCells()
    {
        return this->zonesInterpolateCells;
    }
    ZoneSearchNodes **GetZonesSearchNodes()
    {
        return this->zonesSearchNodes;
    }
    ZoneSearchLines **GetZonesSearchLines()
    {
        return this->zonesSearchLines;
    }
    OversetConfig **GetOversetConfig()
    {
        return this->oversetConfigMultiZone;
    }
    int **GetKeyOfInterfaceNodesInBlocks()
    {
        return this->keyOfInterfaceNodesInBlocks;
    }
    int *GetHasSymmetryBCInBlocks()
    {
        return this->hasSymmetryBCInBlocks;
    }
private:    //! Initialize.
    void BuildNodesMinDistanceInThisBlocks(vector< RDouble > *&wallNodeCoordinateXIn,
        vector< RDouble > *&wallNodeCoordinateYIn,
        vector< RDouble > *&wallNodeCoordinateZIn);
    void ReadInInnerGrid();
    void SetNodesMinDistanceInThisBlocksByWalldist();
    void SetNodesMinDistanceInOtherBlocks();

    bool ReadInInnerGridH5(const string &fileName, UnstructGrid *gridIn);
    void BuildLocalToGlobalInterfaceNodesMap(int iBlockIn, int *numberOfNodesInBlocksIn);

    void GetGlobalInterfaceNodes();
    void GetGlobalInterfaceNodes(int iBlock);
    void ComputeNodeWalldist();
    void ComputeNodeWalldistInOtherBlock();
    void ReComputeWallDistance();
    void ComputeNodeAndCellTopology();

private:    //! Core process: points.
    void SetInitialValue();
    void CalMinBoxForSearch();
    void SetBufferRegionOnOversetBoundary();
    void BuildKDTreesBySearchRegions();
    void BuildNodesMinDistanceInOtherBlocksByKDTrees();
    void SetActivePropertyForSearchRegions();
    void SetSearchRegionByActivePropertys();
    void BuildMinMaxBoxOfZoneBySearchRegions();
    void BuildKDTreeByMinMaxBoxOfZones();
    void BuildInBlockCellsAndNodesByInnerAuxiliaryGridKDTrees(vector< GridKDTree * > &innerAuxiliaryGridKDTreesIn);

    void BuildKDTreeByActiveCells();
    void CollectEffectiveDonorCells();
    void BuildDefectCellRelationships();

    void BuildMinMaxBoxOfBlocks();

    void BuildZonesSearchNodes();
    void SearchZonesSearchNodesInLineBoxes();
    void SearchZonesSearchNodesInLineBoxesByLayer();
    void SearchZonesSearchNodesInKDTrees();
    void AllReduceNodesMinDistanceInOtherBlock();
    void SetNodesMinDistanceInOtherBlockByZonesNodes();

    void SpecifyInitialInterNodesByDistance();
    void SpecifySpecicalCellsInteractWithWallBoundary();
    void SpecifyCellsInBodyByAuxiliaryInnerGrid();
    void FillTheActiveRegion();
    void SetNodeCharacterOnPhysicalBoundary();
    void SetNegativeForBufferRegion();
    void InterFaceBoundaryOptimize();
    void SetNegativeForInnerNodes();
    void SetNegativeForInnerAndInitialInterNodes();

private:    //! Core process: cells.
    void FindInterCellsByBoundary();
    void FindActiveCellsByBoundary();
    //! Assign the centroid coordinates of the contributing element to the interpolation element.
    void ResetCellCentersByBoundary();
    void FindActiveCells();
    void FindInterCells();
    void CreatGlobalInterCells();
    void AllReduceZonesInterpolateCells();
    void AllReduceZonesDefectCells();

private:    //! Core process: cells.
    void BuildInterpolateRelationships();
    void FindInterpolateRelationships();
    void FindDefectRelationships();

private:    //! Check.
    void CheckInterpolateRelationships();
    bool CheckDonorZoneAndDonorCellServer();
    void PrintCellInformation();

private:    //! Communication of interface.
    void SetCellLayersByActiveNodes();
    void CharacterNodeToCell();
    void CharacterCellToNode();
    int  CharacterNodeToValidNodeByCell();

    void CommunicateNodePropertyByKeyOfInterfaceNodesInBlocks();
    void CommunicateNodeCharacterOnInterface();
    void CommunicateCellIBlank();

private:    //! Reconstruction of interpolation information.
    void ResetOversetInformation();
    void GetInterAndDonorCellManagers();
    void GetInterAndDonorCellManagersNew();
    void GetNeighborZoneIndexContainerForTargetZone(int iZone, PHVectorInt1D &neighborZoneList);

private:    //! Output file.
    void OutKeyActive();
    void OutPutTecplot();
    void OutPutBoundaryCells();
    void OutInterCellInformation();
    void OutGlobalInterCells();
    void OutGlobalNodes();

public:
    int GetNumberOfBlocks()
    {
        return numberOfBlocks;
    }

    vector< UnstructGrid * > &GetInnerGrid()
    {
        return innerGrid;
    }
};

class ZoneSearchLines
{
public:
    ZoneSearchLines(int numberOfSearchNodesIn, int numberOfSearchLinesIn, int blockIndexIn);
    ~ZoneSearchLines();

private:
    int numberOfNodes;
    int numberOfLines;
    int blockIndex;

    RDouble *x;
    RDouble *y;
    RDouble *z;

    int *line2Node;
    int *lineIndex;

    int *nodesActive;
    int *linesActive;

    RDouble *pMin;
    RDouble *pMax;

public:
    void SetNumberOfNodes(int numberOfNodesIn)
    {
        this->numberOfNodes = numberOfNodesIn;
    }
    void SetBlockIndex(int blockIndex)
    {
        this->blockIndex = blockIndex;
    }
    int GetNumberOfNodes()
    {
        return numberOfNodes;
    }
    int GetNumberOfLines()
    {
        return numberOfLines;
    }
    int GetBlockIndex()
    {
        return blockIndex;
    }
    RDouble *GetX()
    {
        return x;
    }
    RDouble *GetY()
    {
        return y;
    }
    RDouble *GetZ()
    {
        return z;
    }
    int *GetNodesActive()
    {
        return nodesActive;
    }
    int *GetLinesActive()
    {
        return linesActive;
    }
    int *GetLine2Node()
    {
        return line2Node;
    }
    int *GetLineIndex()
    {
        return lineIndex;
    }
    RDouble *GetPMin()
    {
        return pMin;
    }
    RDouble *GetPMax()
    {
        return pMax;
    }
};

class ZoneSearchNodes
{
public:
    ZoneSearchNodes(int numberOfSearchNodesIn, int blockIndexIn);
    ~ZoneSearchNodes();

private:
    int numberOfNodes;
    int blockIndex;

    RDouble *x;
    RDouble *y;
    RDouble *z;

    RDouble *nodesMinDistanceToOtherBlock;

    RDouble *pMin;
    RDouble *pMax;

public:

    void SetNumberOfNodes(int numberOfNodesIn)
    {
        this->numberOfNodes = numberOfNodesIn;
    }
    void SetBlockIndex(int blockIndex)
    {
        this->blockIndex = blockIndex;
    }

    int GetNumberOfNodes()
    {
        return numberOfNodes;
    }
    int GetBlockIndex()
    {
        return blockIndex;
    }

    RDouble *GetX()
    {
        return x;
    }
    RDouble *GetY()
    {
        return y;
    }
    RDouble *GetZ()
    {
        return z;
    }
    RDouble *GetNodesMinDistanceToOtherBlock()
    {
        return nodesMinDistanceToOtherBlock;
    }

    RDouble *GetPMin()
    {
        return pMin;
    }
    RDouble *GetPMax()
    {
        return pMax;
    }
};

class InterCells
{
public:
    InterCells(int numberOfInterpolateCellsIn, int blockIndexIn);
    ~InterCells();

private:
    int numberOfInterCells;
    int blockIndex;

    int *donorZoneIndex;
    int *donorCellIndex;
    int *isDefectCell;

    int *donorLevel;
    RDouble *cellMinDistance;

    RDouble *x;
    RDouble *y;
    RDouble *z;

public:
    void FreeUselessMemery();

    int GetBlockIndex()
    {
        return blockIndex;
    }

    int GetNumberOfInterCells()
    {
        return numberOfInterCells;
    }

    int *GetDonorZoneIndex()
    {
        return donorZoneIndex;
    }

    int *GetDonorCellIndex()
    {
        return donorCellIndex;
    }

    int *IsDefectCell()
    {
        return isDefectCell;
    }

    RDouble *GetX()
    {
        return x;
    }
    RDouble *GetY()
    {
        return y;
    }
    RDouble *GetZ()
    {
        return z;
    }

    RDouble *GetMinDist()
    {
        return cellMinDistance;
    }
    int *GetKey()
    {
        return donorLevel;
    }

    void SetCellCenterX(RDouble *xIn)
    {
        this->x = xIn;
    }
    void SetCellCenterY(RDouble *yIn)
    {
        this->y = yIn;
    }
    void SetCellCenterZ(RDouble *zIn)
    {
        this->z = zIn;
    }
};

class OversetConfig
{
    typedef void ( OversetConfig:: *BuildInBlockCellsAndNodesMethod )( vector< GridKDTree * > & );
public:
    OversetConfig();
    ~OversetConfig();

public:    //! Necessary global information.
    int numberOfZones;
    int blockIndex;
    int zoneIndex;
    int geometricDimension;
    int donorCellMethod;
    int nLayerOfCellsOfCutting;

    UnstructGrid *grid;

    OversetInformationProxy *oversetInformationProxy;    //! Information of interpolation calculated.
    BuildInBlockCellsAndNodesMethod buildInBlockCellsAndNodesMethod;

public:    //! Improve efficiency of query.
    GridKDTree *gridKDTree;
    RDouble *minBoxOfZone;
    RDouble *maxBoxOfZone;
    int *minBoxNodeID;
    int *maxBoxNodeID;
    bool hasCalMinMaxForZone;

    DataStruct_KDTree *cellKDTree;

    RDouble *minBoxOfZoneByLayer;
    RDouble *maxBoxOfZoneByLayer;

    RDouble *minBoxOfLine;
    RDouble *maxBoxOfLine;

public:
    int numberOfInterCells;

    int *localToGlobalInterfaceNodesMap;
    int *keyActiveOfNodes;
    int *keyActiveOfLines;
    int *keyActiveOfFaces;

    vector< int > interpolateCellList;      //!cells to be interpolated
    vector< int > interpolateFaceList;      //!cells to be interpolated
    vector< int > defectCellList;           //!cells without effective donorCell
    set   < int > interpolateExpendList;    //!cells which are donorCells
    set   < int > donorFaceList;            //!cells which are donorCells

    int *keyBufferOfNodes;

    int *keyBufferOfCells;

    int *keySearchOfNodes;    //! keySearch = 1£ºcells and points used for query.

    int *keySearchOfCells;    //! keySearch = 0£ºpoints in the non overlapping areas (For example, point in the far field).

    int *keySearchOfLines;    //! keySearch = 0£ºpoints in the non overlapping areas (For example, point in the far field).

    int *keySearchOfNodesByLayer;    //! keySearch = 1£ºcells and points used for query.

    int *keySearchOfCellsByLayer;

    int *keyBasicActiveNodes;

    int *nLayerOfCellsFromWall;

    //! the mark of the node which is inside of the body of other block
    //! 0 - outside
    //! 1 - inside
    int *keyInBlockNodes;

    RDouble *nodesMinDistanceInThisBlock;

    RDouble *nodesMinDistanceInOtherBlock;

public:
    vector< RDouble > nodeWalldistDistance;

public:    //! Boundary condition.
    vector< int > interfaceBoundaryNodeList;
    vector< int > physicalBoundaryNodeList;
    vector< int > innerBoundaryNodeList;

    set< int > wallBoundaryNodeSet;
    set< int > intersectLinesWithWall;
    set< int > cellIdBeenSearched;

    vector< int > localToGlobalInterNodesIndex;    //! From local interNode to global interNode.
    vector< int > localToGlobalNodesIndex;         //! From local interNode to global Node.
public:
    void BuildBoundaryNodeList(vector< RDouble > &wallNodeCoordinateXIn,
        vector< RDouble > &wallNodeCoordinateYIn,
        vector< RDouble > &wallNodeCoordinateZIn);
    void BuildNodesMinDistanceInThisBlock(vector< RDouble > *wallNodeCoordinateXInBlocksIn,
        vector< RDouble > *wallNodeCoordinateYInBlocksIn,
        vector< RDouble > *wallNodeCoordinateZInBlocksIn);

    void BuildBoundaryNodeList();
    void SetNodesMinDistanceInThisBlocksByWalldist();
    void SetNodesMinDistanceInOtherBlocks();

    bool hasBeenSearched(int cellIndex)
    {
        set< int >::iterator iter = cellIdBeenSearched.find(cellIndex);
        if (iter != cellIdBeenSearched.end())
        {
            return true;
        }
        return false;
    }

public:    //! Initialize.
    void Initialize(UnstructGrid *gridIn, int nZone, int iZone, int iBlock);
    void SetInitialValue();
    void FreeKeyBasicActiveNode();
    void FreeKeyInBlockNodes();
    void FreeKDTree();
    void SetInitialValueForBufferRigion();
    void SetKeyBufferByKeyActive();

public:    //! Core operations: points.
    void ComputeNodeAndCellTopology();
    void CopyMinDist(RDouble *nodesMinDistanceInOtherBlockIn);
    void CopyNodesActive(int* nodesActiveIn, int* nodes2NodeIn, int numberOfLinesIn);
    void BuildKDTreeByMinMaxBoxOfZone();
    void BuildKDTreeByActiveCell();
    void SearchZonesSearchNodesInLineBox();
    void SearchZonesSearchNodesInLineBoxByLayer();
    void CollectEffectiveDonorCell();
    void SearchZonesSearchNodesInKDTree();

    RDouble ComputeDistanceInDonorZone(RDouble *nodeCoordinateIn, int donorCell, bool wallDist = true);
    RDouble ComputeDistanceInDonorZoneForFace(RDouble *nodeCoordinateIn, int donorFace);
    void SpecifyInitialInterNodesByDistance();
    void SpecifySpecicalCellsInteractWithWallBoundary();
    void SetNegativeForBufferRegion();
    void SetNegativeForInnerNodes();
    void SetNegativeForInnerAndInitialInterNodes();
    void SetNodeCharacterOnPhysicalBoundary();
    void FillTheRemainingNodes();

public:    //! Core operations: globalNodes.
    void BuildMinMaxBoxOfZone();

public:    //! Core operations: points.
    void SetActivePropertyForSearchRegion();
    void SetSearchRegionByActiveProperty();
    void BuildMinMaxBoxOfZoneBySearchRegion();

public:    //! Core operations: cells.

    void DefineInterCellsByBoundary();
    void DefineActCellsByBoundary();
    //! Assign the centroid coordinates of the contributing element to the interpolation element.
    void ResetCellCentersByInterCell();
    void FindInterpolateExpend();
    void FindActiveCells();
    void FindInterCells();
    int  FindDonorCell(RDouble *nodeCoordinateIn);
    void FindInterpolateRelationship();
    void FindDefectRelationship();
    void FindDefectRelationshipByFaces();
    void FindDefectRelationshipByExpend();

    void CheckInterpolateRelationship(int iBlockIn,
        InterCells *zoneInterpolateCellsIn,
        int &numberOfInterpolateCellsInThisBlockIn,
        int &numberOfNegativeCellsInThisBlock_1In,
        int &numberOfNegativeCellsInThisBlock_2In,
        int &numberOfNegativeCellsInThisBlock_3In);

public:    //! Reconstruction of interpolation information.
    void InitializeOversetInformation();
    void PostProcessOversetInformation();

public:    //! Communication.
    void MarkGlobalInterfaceNodes(int iBlock, vector< int > &keyNodesGlobal);
    void CalLocalToGlobalInterfaceNodes(int iBlock, vector< int > &keyNodesGlobal);

    void SetCellLayersByActiveNodes();

    void CharacterNodeToCell();
    void CharacterCellToNode();
    void CharacterInterfaceToNode();
    int  CharacterNodeToValidNodeByCell();
    void UpdateInterfaceNodeCharacter();
    void DownInterfaceNodeCharacter();

public:
    void BuildInBlockCellsAndNodesByInnerAuxiliaryGridKDTree(vector< GridKDTree * > &innerAuxiliaryGridKDTreesIn);

public:    //! Check&Output.
    void OutPutOversetVisualization(const string &name);
    void OutPutOversetVisualizationHeader(fstream &file);
    void OutPutOversetVisualizationCell(fstream &file);
    void OutPutOversetVisualizationPoint(fstream &file);
    void OutPutOversetVisualizationTail(fstream &file);
    void ShowGrid();
    void ShowBoundaryCells();
    void ShowDifferentLayersCell(int iLayer);
    void ShowIntersectCells();
    void ShowIntersectLines(fstream &flowFile);
    void PrintGridByOverSetConditionPoint(fstream &file, int key);
    void PrintGridByOverSetConditionCell(fstream &file, int key);
    void OutKeyBuffer();
    void OutKeyActive();
    void OutputLocalNodeWalldistance();

    inline bool IfZoneSearchNodesIntersectionWithThisZone(RDouble *minBoxOfZoneSearchNodesIn, RDouble *maxBoxOfZoneSearchNodesIn)
    {
        for (int iDimension = 0; iDimension < geometricDimension; ++ iDimension)
        {
            if (maxBoxOfZoneSearchNodesIn[iDimension] < minBoxOfZone[iDimension]) return false;
            if (minBoxOfZoneSearchNodesIn[iDimension] > maxBoxOfZone[iDimension]) return false;
        }
        return true;
    }

    inline bool IfWallFacesIntersectionWithThisZone(RDouble *minBoxOfZoneSearchNodesIn, RDouble *maxBoxOfZoneSearchNodesIn)
    {
        for (int iDimension = 0; iDimension < geometricDimension; ++ iDimension)
        {
            if (maxBoxOfZoneSearchNodesIn[iDimension] < minBoxOfLine[iDimension]) return false;
            if (minBoxOfZoneSearchNodesIn[iDimension] > maxBoxOfLine[iDimension]) return false;
        }
        return true;
    }

    inline bool IfWallFaceIntersectionWithThisZone(Vector3D s1P0, Vector3D s1P1, Vector3D s2P0, Vector3D s2P1)
    {
        Vector3D u = s1P1 - s1P0;
        Vector3D v = s2P1 - s2P0;
        Vector3D w = s1P0 - s2P0;
        Vector3D I0, I1;
        RDouble D = perp(u,v);
        //! test if they are parallel (includes either being a point)
        if (fabs(D) < 0.00000001 )
        {   //! S1 and S2 are parallel
            if (perp(u,w) != 0 || perp(v,w) != 0)
            {
                return 0;    //! they are NOT collinear
            }
            //! they are collinear or degenerate
            //! check if they are degenerate points
            RDouble du = Dot(u,u);
            RDouble dv = Dot(v,v);
            if (du==0 && dv==0)
            {   //! both segments are points
                //! they are distinct points
                if (s1P0 != s2P0)
                {
                    return 0;
                }
                else
                {
                    I0 = s1P0;
                    //! they are the same point
                    return 1;
                }
            }
            if (du==0)
            {   //! S1 is a single point
                if (inSegment(s2P0, s2P1, s1P0) == 0)    //! but is not in S2
                {
                    return 0;
                }
                else
                {
                    I0 = s1P0;
                    return 1;
                }
            }
            if (dv==0)
            {   //! S2 a single point
                if (inSegment(s1P0, s1P1, s2P0) == 0)    //! but is not in S1
                {
                    return 0;
                }
                else
                {
                    I0 = s2P0;
                    return 1;
                }
            }
            //! they are collinear segments - get overlap (or not)
            RDouble t0, t1;    //! endpoints of S1 in eqn for S2
            Vector3D w2 = s1P1 - s2P0;
            if (v.x != 0)
            {
                t0 = w.x / v.x;
                t1 = w2.x / v.x;
            }
            else
            {
                t0 = w.y / v.y;
                t1 = w2.y / v.y;
            }
            if (t0 > t1)
            {   //! must have t0 smaller than t1
                RDouble t=t0; t0=t1; t1=t;    //! swap if not
            }
            if (t0 > 1 || t1 < 0)
            {
                return 0;    //! NO overlap
            }
            t0 = t0<0? 0 : t0;    //! clip to min 0
            t1 = t1>1? 1 : t1;    //! clip to max 1
            if (t0 == t1)
            {   //! intersect is a point
                I0 = s2P0 + v * t0;
                return 1;
            }
            //! they overlap in a valid subsegment
            I0 = s2P0 + v * t0;
            I1 = s2P0 + v * t1;
            return 2;
        }
        //! the segments are skew and may intersect in a point
        //! get the intersect parameter for S1
        RDouble sI = perp(v,w) / D;
        if (sI < 0 - EPSILON || sI > 1 + EPSILON)    //! no intersect with S1
        {
            return 0;
        }
        //! get the intersect parameter for S2
        RDouble tI = perp(u,w) / D;
        if (tI < 0 - EPSILON || tI > 1+ EPSILON)    //! no intersect with S2
        {
            return 0;
        }

        I0 = s1P0 + u * sI;
        return 1;

    }
    int inSegment(Vector3D s1P0, Vector3D s1P1, Vector3D P)
    {
        if (s1P0.x != s1P1.x)
        {   //! S is not vertical
            if (s1P0.x <= P.x && P.x <= s1P1.x)
                return 1;
            if (s1P0.x >= P.x && P.x >= s1P1.x)
                return 1;
        }
        else
        {   //! S is vertical, so test y coordinate
            if (s1P0.y <= P.y && P.y <= s1P1.y)
                return 1;
            if (s1P0.y >= P.y && P.y >= s1P1.y)
                return 1;
        }
        return 0;
    }

    inline int LineIntersectionWithFace(Vector3D line2Node1, Vector3D line2Node2, Vector3D face2Node1, Vector3D face2Node2, Vector3D face2Node3)
    {
        vector< RDouble> nodeCoordinate1;
        nodeCoordinate1.push_back(face2Node1.x);
        nodeCoordinate1.push_back(face2Node1.y);
        nodeCoordinate1.push_back(face2Node1.z);

        vector< RDouble> nodeCoordinate2;
        nodeCoordinate2.push_back(face2Node2.x);
        nodeCoordinate2.push_back(face2Node2.y);
        nodeCoordinate2.push_back(face2Node2.z);

        vector< RDouble> nodeCoordinate3;
        nodeCoordinate3.push_back(face2Node3.x);
        nodeCoordinate3.push_back(face2Node3.y);
        nodeCoordinate3.push_back(face2Node3.z);

        RDouble *minBoxOfTemp = new RDouble[geometricDimension];
        RDouble *maxBoxOfTemp = new RDouble[geometricDimension];

        SetField(minBoxOfTemp, LARGE, geometricDimension);
        SetField(maxBoxOfTemp, -LARGE, geometricDimension);
        for (int iDimension = 0; iDimension < geometricDimension; ++iDimension)
        {
            RDouble& minInThisDimension = minBoxOfTemp[iDimension];
            RDouble& maxInThisDimension = maxBoxOfTemp[iDimension];

            minInThisDimension = MIN(minInThisDimension, nodeCoordinate1[iDimension]);
            maxInThisDimension = MAX(maxInThisDimension, nodeCoordinate1[iDimension]);

            minInThisDimension = MIN(minInThisDimension, nodeCoordinate2[iDimension]);
            maxInThisDimension = MAX(maxInThisDimension, nodeCoordinate2[iDimension]);

            minInThisDimension = MIN(minInThisDimension, nodeCoordinate3[iDimension]);
            maxInThisDimension = MAX(maxInThisDimension, nodeCoordinate3[iDimension]);

            /*minInThisDimension -= eps;
            maxInThisDimension += eps;*/
        }

        if (IfWallFacesIntersectionWithThisZone(minBoxOfTemp, maxBoxOfTemp) == false)
        {
            return 0;
        }

        RDouble tmp1=0.0, tmp2=0.0, tmp3=0.0, area = 0.0;
        Vector3D fn, lineDistance, distanceNode1, distanceNode2, nodeCoordinate;
        RDouble *fnTemp = new RDouble[3]();

        area = ComputeTriangleArea(face2Node1, face2Node2, face2Node3, fn, fnTemp);

        fn.x = fnTemp[0];
        fn.y = fnTemp[1];
        fn.z = fnTemp[2];

        lineDistance = line2Node1 - line2Node2;

        tmp3 = lineDistance * fn;

        RDouble lineDistanceSquare = lineDistance * lineDistance;

        if (tmp3 * tmp3 / lineDistanceSquare < EPSILON)
        {
            return 0;
        }

        distanceNode1 = line2Node1 - face2Node1;
        distanceNode2 = line2Node2 - face2Node1;

        tmp1 = distanceNode1 * fn;
        tmp2 = distanceNode2 * fn;
        if (tmp1 * tmp2 > 0)
        {
            return 0;
        }

        tmp3 = tmp1/tmp3;
        nodeCoordinate = line2Node1 - lineDistance * tmp3;

        RDouble areaSum = 0.0;
        Vector3D tn;
        RDouble *tnTemp = new RDouble[3]();
        areaSum += ComputeTriangleArea(face2Node1, face2Node2, nodeCoordinate, tn, tnTemp);
        areaSum += ComputeTriangleArea(face2Node2, face2Node3, nodeCoordinate, tn, tnTemp);
        areaSum += ComputeTriangleArea(face2Node3, face2Node1, nodeCoordinate, tn, tnTemp);

        if (areaSum > area * (1.0 + EPSILON) )
        {
            return 0;
        }
        else
        {
            return 1;
        }

        delete [] minBoxOfTemp;    minBoxOfTemp = nullptr;
        delete [] maxBoxOfTemp;    maxBoxOfTemp = nullptr;

        delete [] fnTemp;    fnTemp = nullptr;
        delete [] fnTemp;    fnTemp = nullptr;

        delete [] tnTemp;    tnTemp = nullptr;
        delete [] tnTemp;    tnTemp = nullptr;

    }
    inline RDouble ComputeTriangleArea(Vector3D face2Node1, Vector3D face2Node2, Vector3D face2Node3, Vector3D fn,  RDouble* fnTemp)
    {
        Vector3D r1, r2;
        RDouble area = 0.0;

        r1 = face2Node2 - face2Node1;
        r2 = face2Node3 - face2Node1;

        fn = Cross(r1, r2);
        area = sqrt(fn * fn);
        fn /= area;

        area *= 0.5;

        fnTemp[0] = fn.x;
        fnTemp[1] = fn.y;
        fnTemp[2] = fn.z;

        return area;
    }

public:
    int *GetKeyActiveOfNodes()
    {
        return keyActiveOfNodes;
    }
    int *GetKeyActiveOfCells();
    int *GetKeyBufferOfCells()
    {
        return keyBufferOfCells;
    }
    int *GetKeyBufferOfNodes()
    {
        return keyBufferOfNodes;
    }
    int *GetKeySearchOfCells()
    {
        return keySearchOfCells;
    }
    int *GetKeySearchOfNodes()
    {
        return keySearchOfNodes;
    }

    int *GetKeySearchOfLines()
    {
        return keySearchOfLines;
    }

    int *GetLayerOfCells()
    {
        return nLayerOfCellsFromWall;
    }

    vector< int > &GetInterfaceBoundaryNodeList()
    {
        return interfaceBoundaryNodeList;
    }

    int *GetLocalToGlobalInterfaceNodesMap()
    {
        return localToGlobalInterfaceNodesMap;
    }

    RDouble *GetPMin()
    {
        return minBoxOfZone;
    }
    RDouble *GetPMax()
    {
        return maxBoxOfZone;
    }
    RDouble *GetNodeWalldistDistance()
    {
        return nodesMinDistanceInOtherBlock;
    }

    int GetBlockIndex()
    {
        return blockIndex;
    }

    int GetNumberOfInterCells()
    {
        return numberOfInterCells;
    }

    set< int > &GetInterCellExpendList()
    { 
        return this->interpolateExpendList;
    }

    set< int > &GetDonorFaceListList()
    { 
        return this->donorFaceList;
    }

    vector< int > &GetInterpolateCellList()
    {
        return this->interpolateCellList;
    }

    vector< int > &GetInterpolateFaceList()
    {
        return this->interpolateFaceList;
    }

    vector< int > &GetDefectCellList()
    {
        return this->defectCellList;
    }

    OversetInformationProxy *GetOversetInformationProxy()
    {
        return oversetInformationProxy;
    }
    UnstructGrid *GetGrid()
    {
        return grid;
    }
private:
    inline void OutPutValueName(fstream &file, const string &name)
    {
        file << "# " << name << endl;
    }

    inline void OutPutNullValue(fstream &file, int length)
    {
        int numberOfWordsInEachLine = 5;
        int iPrint = 0;
        for (int iValue = 0; iValue < length; ++ iValue)
        {
            file << 0 << " ";
            ++ iPrint;
            if (iPrint % numberOfWordsInEachLine == 0) file << endl;
        }
        if (iPrint % numberOfWordsInEachLine != 0) file << endl;
        file << endl;
    }

    inline void OutPutVectorValue(fstream &file, int length, vector< int > &value)
    {
        bool *valueFull = new bool[length];
        SetField(valueFull, false, length);
        int numberOfValue = static_cast<int>(value.size());
        for (int iValue = 0; iValue < numberOfValue; ++ iValue)
        {
            valueFull[value[iValue]] = true;
        }

        int numberOfWordsInEachLine = 5;
        int iPrint = 0;
        for (int iValue = 0; iValue < length; ++ iValue)
        {
            file << valueFull[iValue] << " ";
            ++ iPrint;
            if (iPrint % numberOfWordsInEachLine == 0) file << endl;
        }
        if (iPrint % numberOfWordsInEachLine != 0) file << endl;
        file << endl;

        delete [] valueFull;    valueFull = nullptr;
    }

    template < typename T >
    inline void OutPutPointerValue(fstream &file, int length, T *value)
    {
        if (NULL == value)
        {
            OutPutNullValue(file, length);
        }
        else
        {
            int numberOfWordsInEachLine = 5;
            int iPrint = 0;
            for (int iValue = 0; iValue < length; ++ iValue)
            {
                file << value[iValue] << " ";
                ++ iPrint;
                if (iPrint % numberOfWordsInEachLine == 0) file << endl;
            }
            if (iPrint % numberOfWordsInEachLine != 0) file << endl;
            file << endl;
        }
    }
};

RDouble GetToleranceForOversetBox();
int GetKeyEnlargeOfInterBoundary();
int GetNCommunicateCellIBlank();
int GetKeyEnlargeOfNodes();
int GetKeyEnlargeOfSearchNodes();
int GetReadInAuxiliaryInnerGrid();
int GetReadInAuxiliaryOuterGrid();
int GetReadInSklFileOrNot();
int GetSymetryOrNot();
int GetIfOutPutOversetVisualization();

template< typename T >
void PH_AllreduceSepModeForVector(vector< T > &sendingBuffer, vector< T > &receivingBuffer, int numberOfElements, MPI_Datatype mpiDataType, PH_Op op);

}