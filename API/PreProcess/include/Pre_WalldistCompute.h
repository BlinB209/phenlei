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
//! @file      Pre_WalldistCompute.h
//! @brief     Wall distance compute for both structured and unstructured grid.
//! @author    Bell.

#pragma once
#include "GridType.h"
#include "TK_Time.h"
#include "Math_BasisFunction.h"
#include "Geo_StructGrid.h"
#include "Geo_UnstructGrid.h"
#include "PHHeader.h"

using namespace std;

namespace PHSPACE
{
class WallStructure;

//! @brief Pre_WalldistCompute class defines the method of 
//! wall distance compute for both structured and unstructured grid.
//! Wall Distance: the most nearest distance of each cell center to the solid wall.
class Pre_WalldistCompute
{
    typedef pair<std::size_t, int> FaceID;
    typedef DataStruct_AdtNode<FaceID, RDouble> WallFaceNode;
    typedef DataStruct_AdtTree<FaceID, RDouble> WallFaceTree;
    typedef DataStruct_AdtNode<int, RDouble>    WallStructNode;
    typedef DataStruct_AdtTree<int, RDouble>    WallStructTree;
    typedef DataStruct_KDTree                   WallStructKDTree;
    typedef vector<RDouble> value_type;
    static const int ACCURATE_METHOD   = 0;
    static const int FAST_METHOD       = 1;
    static const int SUPER_FAST_METHOD = 2;
    static const int KDTREE_METHOD     = 3;

    static const int CELL_METHOD = 0;
    static const int NODE_METHOD = 1;

private:
    //! Wall distance computing method.
    //! -# 0: fast but not accurate enough.
    //! -# 1: accurate but not fast enough.
    //! -# 3: k-dimension method.
    int walldistComputeMethod;

    int cellMethodOrNodeMethod;

    //! Grid level in multi-grid method.
    //! 0 represents the default initial finest computational mesh.
    //! 1, 2, ..., N is the coarse level mesh.
    int level;

    int presentBlock;

public:
    //! @param[in] walldistComputeMethod    Wall distance computing method.
    //!                                     -# 0: accurate but not fast enough.
    //!                                     -# 1: fast but not accurate enough.
    //! @param[in] level    Grid level in multi-grid method.
    //! 0 represents the default initial finest computational mesh.
    //! 1, 2, ..., N is the coarse level mesh.
    LIB_EXPORT Pre_WalldistCompute(int walldistComputeMethod, int level);

    LIB_EXPORT ~Pre_WalldistCompute();

public:
    //! Run to compute the wall distance of all grids.
    LIB_EXPORT void Run();

    //! Compute the wall distance of nodes for overlap config.
    LIB_EXPORT void ComputeNodeDistance();

    LIB_EXPORT void ComputeNodeDistanceInOtherBlock();

    LIB_EXPORT void ComputeWallDistForNodesList(int blockIndex, int nTotalNode, RDouble *x, RDouble *y, RDouble *z, RDouble *dist);

    vector<WallStructure *> * GetWallstructureList();
private:
    //! Fill wall structure by server collection-bcast method.
    void FillWallStructureNew();
    void ServerCollection();
    void ServerBcast();

    void InitWallDistance();
    void InitWallDistance(Grid *grid);
    void InitWallDistance(StructGrid *grid);
    void InitWallDistance(UnstructGrid *grid);

    void CheckWallDistance();
    void CheckWallDistance(UnstructGrid *grid);
    void CheckWallDistanceByNodeMethod(UnstructGrid *grid);

    void PostWallDistance();
    void PostWallDistance(Grid *grid);
    void PostWallDistance(StructGrid *grid);
    void PostWallDistance(UnstructGrid *grid);
    void PostWallDistanceByNodeMethod(UnstructGrid *grid);

    bool IsWalldistComputOver();
    bool IsWalldistComputOver(Grid *grid);
    bool IsWalldistComputOver(StructGrid *grid);
    bool IsWalldistComputOver(UnstructGrid *grid);
    bool IsNodeWalldistComputOver();

    //! Fill wall structure by bcast wall struct to each other.
    void FillWallStructure();
    void ComputeWallDistance();
    void ComputeWallDistance(Grid *grid);
    void ComputeWallDistance(StructGrid *grid);
    void ComputeWallDistanceKDTree(StructGrid *grid);
    void ComputeWallDistance(UnstructGrid *grid);
    void ComputeWallDist2DByBinarySearching(UnstructGrid *grid);
    void ComputeWallDist3DByProjectingAccurate(UnstructGrid *grid);
    void ComputeWallDist3DByProjectingAccurateAndNodeMethod(UnstructGrid *grid);
    void ComputeWallDist3DByProjectingFast(UnstructGrid *grid);
    void ComputeWallDist3DByProjectingFastAndNodeMethod(UnstructGrid *grid);
    void ComputeWallDist3DByProjectingFast(StructGrid *grid);
    void ComputeNodeWallDist3DByProjectingFast(int nTotalNode, RDouble *x, RDouble *y, RDouble *z, RDouble *nodeWalldist);
    void ComputeNodeWallDist3DByKDTreeMethod(int nTotalNode, RDouble *x, RDouble *y, RDouble *z, RDouble *nodeWalldist);
    void CompressWallStructure(DataContainer *&cdata, Grid *grid);
    void CompressWallStructure(DataContainer *&cdata, StructGrid *grid);
    void CompressWallStructure(DataContainer *&cdata, UnstructGrid *grid);

    void FilterWallStructureList();
    void FindNearstWallStructureOld(RDouble *cellCenter, vector<WallStructure *> &nearestWallStruct);

    void DeCompressWallStructure(DataContainer *&cdata, const bool &ifMoveToBegin);

    void BuildWallStructTree();
    void BuildWallStructKDTree();
    void BuildWallStructKDTreeNode();
    void BuildWallStructKDTreeNodeInOtherBlock();

    RDouble ComputeDistanceOfANodeToSegment(const RDouble xcc, const RDouble ycc, const RDouble zcc,
                                            const RDouble xStart, const RDouble yStart, const RDouble zStart,
                                            const RDouble xMid, const RDouble yMid, const RDouble zMid,
                                            const RDouble xEnd, const RDouble yEnd, const RDouble zEnd);
    RDouble ComputeDistanceOfANodeToAFace(RDouble *cellCenter, RDouble *faceCenter, RDouble *faceNormal, 
                                          const vector<int> &face2node, value_type &x, value_type &y, value_type &z);

    void RotateAxis(const RDouble *oldR, RDouble *newR);

    RDouble GetGlobalGridSize(int type = 0);

    void ComputeNodeDistanceByDirectMethod();
    void ComputeNodeDistanceByADTMethod();
    void ComputeNodeDistanceByKDTreeMethod();
    void ComputeNodeDistanceByKDTreeMethodInOtherBlock();
    void ComputeNodeDistanceByDirectMethod(UnstructGrid *grid);
    void ComputeNodeDistanceForSphere(UnstructGrid *grid);
    void ComputeNodeDistanceByADTMethodForBlock(int iBlock);
    void ComputeNodeDistanceByKDTreeMethodForBlock();
    void ComputeNodeDistanceByKDTreeMethodForBlockInOtherBlock();
    void ComputeWallDistForNodesListByADT(int blockIndex, int nTotalNode, RDouble *x, RDouble *y, RDouble *z, RDouble *dist);
    void ComputeWallDistForNodesListByKDT(int blockIndex, int nTotalNode, RDouble *x, RDouble *y, RDouble *z, RDouble *dist);
    void ComputeWallDistForNodesListByDirectMethod(int blockIndex, int nTotalNode, RDouble *x, RDouble *y, RDouble *z, RDouble *dist);
    void ComputeWallDistForNodesListForSphere(int blockIndex, int nTotalNode, RDouble *x, RDouble *y, RDouble *z, RDouble *dist);

    int  JudgeInitialRotateAxis();
private:
    //! The data structure to store the wall faces and nodes, and the face-nodes connection.
    vector<WallStructure *> *wallstructureList;

    //! The current processor index.
    int myid;

    //! Number of total wall faces.
    int nTWFace;

    //! Boundary box of wall structure.
    RDouble wallMin[6], wallMax[6];

    //! ADT Tree of wall structures.
    WallFaceTree *wallFaceTree;

    //! K-Dimendion Tree of wall structures.
    WallStructKDTree *wallFaceKDTree;

    //! K-Dimendion Tree of wall structures.
    WallStructKDTree *wallNodeKDTree;

    //! Tolerance of wall.
    RDouble wallTolerance;

    //! Number of total faces on wall.
    size_t nTWallFace;

    //! Rotate Axis:
    //!   0 -- withot rotate.
    //!   1 -- axis X.
    //!   2 -- axis Y.
    //!   3 -- axis Z.
    int rotateAxis;

    //! Rotate degree.
    RDouble rotateDegree;

    //! Time statistical.
    TimeSpan timeSpan;
};

//! @brief CWallDist class defines the basic wall face data structure:
//! 1: wall face center, for both structured and unstructured grid.
//! 2: wall points coordinate and the face-to-points connection, for unstructured grid only.
class WallStructure
{
public:
    typedef vector<RDouble> value_type;
    WallStructure();
    ~WallStructure();

private:
    //! Coordinates of the points on the solid wall.
    value_type x, y, z;

    //! Center of the faces on the solid wall.
    value_type xfc, yfc, zfc;

    //! Normal of the faces on the solid wall.
    value_type xfn, yfn, zfn;

    //! Number of the points per faces on the solid wall.
    int *nPointPerFace;

    //! The points list of faces on the solid wall.
    int *wallFace2Node;

    //! Node position of each face.
    //! wallFace2Node[nodePosition[iFace]]: start point ID.
    //! wallFace2Node[nodePosition[iFace + 1]]: end point ID.
    int *nodePosition;

    //! Boundary box.
    RDouble box[6];

    //! Distance between box boundary and a filed point.
    RDouble distance;

    //! The block that this wallstructure belongs to.
    int iBlock;

public:
    std::size_t GetNumberOfWallFaces();
    std::size_t GetNumberOfWallPoints();

    value_type & GetX();
    value_type & GetY();
    value_type & GetZ();

    value_type & GetXFaceCenter();
    value_type & GetYFaceCenter();
    value_type & GetZFaceCenter();

    value_type & GetXFaceNormal();
    value_type & GetYFaceNormal();
    value_type & GetZFaceNormal();

    RDouble * GetBox();
    RDouble GetDistance();

    int * GetWallFace2Node();
    int * GetnPointPerFace();
    int * GetNodePosition();

    void SetWallFace2Node(int *wallFace2Node);
    void SetnPointPerFace(int *nPointPerFace);
    void SetNodePosition(int *nodePosition);
    void SetDistance(RDouble dist);

    void Box();

    void Encode(DataContainer *&cdata, const bool &ifMoveToBegin);
    void Decode(DataContainer *&cdata, const bool &ifMoveToBegin);
    
    vector<vector<int> > & GetFace2NodeVector();
    void GetFace2NodeVector(vector<int> &face2node, int faceID);

    int  GetIBlock();
    void SetIBlock(int iBlockIn);
private:
    //! Vector type of face to node data structure.
    vector<vector<int> > face2node;
};


class AleWalldistManager
{
private:
    int storeAleWalldis;
    int iterNumOneCycle;

public:
    AleWalldistManager();
    ~AleWalldistManager();

public:
    void Initialize();
    void Run();

private:
    bool JudgeWalldistExist();
    void ReadWalldist();
    void GenerateAleWalldist();

private:
    int iterFileIndex;
    string walldistFile;
};

const int DUMPAllWALLDISTINFO = 0;
const int DUMPWALLDISTONLY    = 1;

//! Process of generate wall distance.
void GenerateWalldist();

void AleWalldistGenerate();

#ifdef AI_RandomForestRegressor
void DumpCellToWallNormal(vector<vector<RDouble>> allcelltowallnormal);

void DumpMaxBoundaryLayerWalldist();
#endif

//!
void AllocateWalldist();

//! Compute wall distance.
void ComputeWalldist();

//! Write wall distance.
void WriteWalldist(int dumpData = DUMPAllWALLDISTINFO);

void CompMaxBoundaryLayerWalldist();

bool IsBoundaryLayerExist();

//!
void WalldistStatics();

//! Dump out the wall distance to file.
void DumpWalldistToFile(ActionKey *actkey, int iZone);

#ifdef USE_TecplotLib
//! Visualization of wall distance.
void VisualWalldist();
#endif

//! If read wall distance file or not.
bool IsReadWalldist();

#include "Pre_WalldistCompute.hxx"

}