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
//! @file      Pre_Block_Struct.h
//! @brief     Structured grid partition: block definition.
//! @author    Guo Yongheng.

#pragma once
#include "Geo_StructBC.h"

using namespace std;
namespace PHSPACE
{
class Pre_Block_Struct;
class Pre_Patch_Struct;

//! @brief Pre_BcFace_Struct defines boundary face for struct grid.
class Pre_BcFace_Struct
{
private:
    //! Boundary type of current face.
    int boundaryType;

    string bcName;

    int s_dir3d[3];

    //! Range of node index and local axis label.
    vector<int> nodeStart, nodeEnd, axisLabel;

    //! Boundary Direction.
    int bcDirection;

    //! Unit matrix which is used to describe axis mapping between source block to target one.
    vector<vector<int> > unitMatrix;

    //! Pointer of attached struct grid.
    Pre_Block_Struct *simpleBlock;

    //! Pointer of correlative patch face.
    Pre_Patch_Struct *next;

    //! Children of current boundary face.
    vector<Pre_BcFace_Struct *> *child;
public:
    LIB_EXPORT  Pre_BcFace_Struct();
    LIB_EXPORT ~Pre_BcFace_Struct();
public:
    //! Recursive function for freeing space of Bc faces and patch faces.
    LIB_EXPORT void Free();

    //! Erect vertex index mapping between patch faces.
    LIB_EXPORT void ErectVertexMapping(int geometricDimension);

    //! Get absolute value of node index and sort them.
    LIB_EXPORT void Normalize(int geometricDimension);

    //! Compute local node index according to block parameter.
    LIB_EXPORT void ComputeLocalCoordinate();

    //! Output current boundary face information to file.
    LIB_EXPORT void OutputBoundaryInformation(fstream &file, int geometricDimension);
    LIB_EXPORT void OutputBoundaryInformation(StructBC *bcregion, int geometricDimension);

    //! Extract leaves of Bc faces which can not be split any more.
    LIB_EXPORT vector<Pre_BcFace_Struct *> GetLeaves();

    void SetRegion(vector<int> &nodeStart, vector<int> &nodeEnd);
    void SetDirection(int bcDirection) { this->bcDirection =bcDirection; }
    void SetAxisLabel(vector<int> &axisLabel);
    void SetSimpleBlock(Pre_Block_Struct *simpleBlock);
    void SetNext(Pre_Patch_Struct *next);
    void SetChild(vector<Pre_BcFace_Struct *> *child);
    void SetBoundaryType(int boundaryType);
    void SetBoundaryName(string bcName);
    void SetFaceMatchingTargetDirIndex(int *dir3dIn);
    vector<int> & GetNodeStart();
    vector<int> & GetNodeEnd();
    vector<int> & GetAxisLabel();
    Pre_Block_Struct * GetSimpleBlock();
    Pre_Patch_Struct * GetNext();
    int GetBoundaryType();
    int GetDirection() { return bcDirection; }
    string GetBoundaryName();
    int * GetFaceMatchingTargetDirIndex();
    vector<Pre_BcFace_Struct *> * GetChild();
private:
    void GetSurfaceDirection(vector<int> &st, vector<int> &ed, vector<int> &nd, int bcDirection, int geometricDimension);
    void InitUnitMatrix();
};

//! @brief Pre_Patch_Struct defines patch face for current boundary face.
class Pre_Patch_Struct
{
private:
    //! Correlative block index.
    int targetBlockLabel;

    //! Range of node index and local axis label.
    vector<int> nodeStart, nodeEnd, axisLabel;

    //! Boundary Direction.
    int bcDirection;

    //! Node index mapping operator.
    vector<vector<int> > frameOfAxes;

    //! Pointer of target struct block;
    Pre_Block_Struct *simpleBlock;
public:
    LIB_EXPORT  Pre_Patch_Struct();
    LIB_EXPORT ~Pre_Patch_Struct();
public:
    void SetRegion(vector<int> &nodeStart, vector<int> &nodeEnd);
    void SetNodeMapping(vector<vector<int> > &frameOfAxes);
    void SetAxisLabel(vector<int> &axisLabel);
    void SetSimpleBlock(Pre_Block_Struct *simpleBlock);
    void SetTargetBlockLabel(int i);
    void SetDirection(int bcDirection) { this->bcDirection =bcDirection; }
    int GetDirection() { return bcDirection; }

    vector<int> & GetNodeStart();
    vector<int> & GetNodeEnd();
    vector<int> & GetAxisLabel();

    Pre_Block_Struct * GetSimpleBlock();
    vector<vector<int> > & GetFrameOfAxes();
};

//! @brief Pre_Block_Struct defines basal structured grid.
class Pre_Block_Struct
{
private:
    //! For multigrid case, numberOfUnitCells equals 2**(numberOfMultigrid - 1).
    static int numberOfUnitCells;

    //! Current zone index.
    int zoneIndex;
    
    //! Original zone index.
    int originalZoneIndex;
    
    //! Processor index of current zone.
    int processorIndex;

    //! Parent pointer of current zone.
    Pre_Block_Struct *parent;

    //! Child pointer of current zone.
    vector<Pre_Block_Struct *> *child;

    //! Bc face pointer group of current zone.
    vector<Pre_BcFace_Struct *> boundaryFaceList;

    //! Relative start and end node index in parent zone.
    vector<int> originalIndex;

    //! For current zone, node number along i,j,k directions.
    vector<int> nodeDimension;

    //! For original zone, ni, nj, nk are node number along i, j, k directions.
    int ni, nj, nk;

    //! For original zone, x, y, z are node coordinates with circulation k, j, i.
    vector<RDouble> x, y, z;
public:
    LIB_EXPORT  Pre_Block_Struct();
    LIB_EXPORT ~Pre_Block_Struct();
public:
    //! Calculate number of unit cells for multi-grid case.
    LIB_EXPORT static void CalculateNumberOfUnitCells(int numberOfMultigrid);

    //! Generate space for node coordinates, i.e. x, y and z.
    LIB_EXPORT void GenerateCoordinateSpace();

    //! Split current block by the parameter "numberOfCells" and compute relative locations for the left and right block.
    LIB_EXPORT void Split(RDouble numberOfCells, vector<int> &leftDimension, vector<int> &rightDimension, vector<int> &leftOriginalIndex, vector<int> &rightOriginalIndex);
    
    //! Create boundary face for new zones.
    LIB_EXPORT void CreateBoundaryFace(int geometricDimension);

    //! Normalize every boundary face for current zone.
    LIB_EXPORT void NormalizeBoundaryFace(int geometricDimension);

    void PatchChild(vector<int> &t_st, vector<int> &t_ed, Pre_Block_Struct *sourceBlock, int geometricDimension);
    void CreateBoundaryFace(Pre_BcFace_Struct *boundaryFace, int geometricDimension);
    void GetAbsoluteCoordinate(vector<int> &st, vector<int> &ed);
    void GetLocalCoordinate(vector<int> &st, vector<int> &ed, vector<int> &local_st, vector<int> &local_ed);
    void OutputBoundaryInformation(fstream &file, int geometricDimension);
    void GetRootInformation(Pre_Block_Struct *&root, vector<int> &ref);
    void GenerateBoundaryFaceList(int numberOfBoundaryFaces);
    void WriteGridCoordinate(fstream &file, vector<int> &nst, vector<int> &ned);
    void CorrectSplitDimension(int axisSelect, RDouble numberOfCells, int geometricDimension, int &isplit);

    static int GetNumberOfUnitCells();
    int GetNumberOfNodes();
    int GetNumberOfCells();
    void AddBoundaryFace(Pre_BcFace_Struct *boundaryFace);

    void SetZoneIndex(int iZone);
    void SetProcessorIndex(int iProcessor);
    void SetOriginalZoneIndex(int iZone);
    void SetParentBlock(Pre_Block_Struct *parent);
    void SetBoundaryFace(int i, Pre_BcFace_Struct *boundaryFace);

    int GetZoneIndex();
    int GetOriginalZoneIndex();
    Pre_BcFace_Struct * GetBoundaryFace(int i);
    vector<Pre_BcFace_Struct *> GetBoundaryFaceList() { return this->boundaryFaceList; }

    void SetOriginalIndex(vector<int> &od);
    void SetNodeDimension(vector<int> &nd);

    vector<int> & GetOriginalIndex();
    vector<int> & GetNodeDimension();

    int GetNI() const;
    int GetNJ() const;
    int GetNK() const;

    void SetNI(int ni);
    void SetNJ(int nj);
    void SetNK(int nk);

    RDouble * GetX();
    RDouble * GetY();
    RDouble * GetZ();
private:
    bool IsOverlapBoundaryFace(vector<int> &st, vector<int> &ed, int geometricDimension);
    bool IsPoleBoundaryFace(int iDimension);
    bool PatchChild(Pre_BcFace_Struct *boundaryFace, vector<int> &t_st, vector<int> &t_ed, Pre_Block_Struct *sourceBlock, int geometricDimension);
    void GerneralSplit(Pre_BcFace_Struct *boundaryFace, Pre_BcFace_Struct *overlapBoundaryFace, int geometricDimension);
    void PatchSplit(Pre_BcFace_Struct *boundaryFace, Pre_BcFace_Struct *overlapBoundaryFace, int geometricDimension);
    void Free();
    Pre_BcFace_Struct * GetOverlapBoundaryFace(vector<int> &st, vector<int> &ed, vector<int> &st_ref, vector<int> &ed_ref, int geometricDimension);
};

#include "Pre_Block_Struct.hxx"

int GetSurfaceLabel(vector<int> &st, vector<int> &ed);
int SignFunction(int value);
void CalculateVertexIndexOnTargetBlock(vector<vector<int> > &frameOfAxes, vector<int> &p1, vector<int> &p2);
void GetIJKRegion(vector<int> &blockNodeStart, vector<int> &blockNodeEnd, vector<int> &faceNodeStart, vector<int> &faceNodeEnd, int localFaceLabel);
}