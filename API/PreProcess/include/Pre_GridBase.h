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
//! @file      Pre_GridBase.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include <vector>
#include <set>
#include <map>
#include <string>
#include "Precision.h"
#include "cgnslib.h"
#include "Geo_UnstructGrid.h"

using namespace std;

namespace PHSPACE
{

class CGNSRawCoor
{
public:
    void * GetCoor(int icoor) { return coor[icoor]; }
private:
    cgsize_t nTotalNode;
    DataType_t *type;
    void *coor[3];
public:
    CGNSRawCoor();
    ~CGNSRawCoor();
public:
    void SetAllData(RDouble *x, RDouble *y, RDouble *z);
    void AllocateData(int icoor, cgsize_t nTotalNode, DataType_t type);
private:
    void SetData(int icoor, DataType_t type, RDouble *x);
    void DeAllocateData();
};

class RawGrid
{
public:
    RawGrid();
    ~RawGrid();
private:
    cgsize_t nTotalNode;
    RDouble *x, *y, *z;
public:
    cgsize_t GetNTotalNode() const { return nTotalNode; }
    void SetNTotalNode(cgsize_t nTotalNode) { this->nTotalNode = nTotalNode; }
    RDouble * GetX() const { return x; }
    RDouble * GetY() const { return y; }
    RDouble * GetZ() const { return z; }
    void SetX(RDouble *x) { this->x = x; }
    void SetY(RDouble *y) { this->y = y; }
    void SetZ(RDouble *z) { this->z = z; }
};

class MultiIndex
{
public:
    int size, id;
    int *data;
public:
    MultiIndex();
    MultiIndex(int size, int id = 0);
    MultiIndex(const MultiIndex &rhs);
    MultiIndex & operator = (const MultiIndex &rhs);
    ~MultiIndex();
public:
    bool operator < (const MultiIndex &rhs) const;
    void SetData(int *data);
};

class CElemFace
{
public:
    CElemFace();
    ~CElemFace();
private:
    int face_num;
    int *face_clkwise;
    int *facept_num;
    int **face_id;
public:
    int ** GetFaceID() const { return face_id; }
    int * GetFacePointNumber() const { return facept_num; }
    int GetFaceNumber() const { return face_num; }
public:
    void Init(int elem_type);
};

class BaseElement
{
public:
    BaseElement();
    ~BaseElement();
private:
    int size;
    CElemFace *elementface;
    cgsize_t *elementpoint;
public:
    cgsize_t * GetElementPoint() const { return elementpoint; }
    CElemFace * GetElementFace() const { return elementface; }
private:
    void InitElementFace();
    void InitElementPoint();
};

class BaseBCType
{
public:
    BaseBCType();
    ~BaseBCType();
    void Init();
private:
    int ibcmax;
    map<int, int> bctypemap;
    map<int, int> fts2cgns;
public:
    int GetBCType(int bctype) { return bctypemap[bctype]; }
    int Fantasy2CGNS(int bctype) { return fts2cgns[bctype]; }
};

class Base_Grid_Conn;

class CGNSBase
{
public:
    CGNSBase();
    ~CGNSBase();
public:
    void SetNTotalCell(cgsize_t nTotalCell) { this->nTotalCell = nTotalCell; }
    cgsize_t GetNTotalCell() const { return nTotalCell; }

    void SetNTotalNode(cgsize_t nTotalNode) { this->nTotalNode = nTotalNode; }
    cgsize_t GetNTotalNode() const { return nTotalNode; }

    void SetNBoundFace(cgsize_t nBoundFace) { this->nBoundFace = nBoundFace; }
    cgsize_t GetNBoundFace() const { return nBoundFace; }
private:
    cgsize_t nTotalNode, nTotalCell, nBoundFace;
    int      nSections, nBCRegions;
    string   *boundaryName;
    string   volumeName;
    cgsize_t *iStart, *iEnd;
    cgsize_t *elementDataSize;
    int      *elementType;
    cgsize_t **connList;
    int      *nBCElem;
    cgsize_t **mixedConnID;

    int *bcType;
    int *boundaryElementType, * bcGridLocation;
    cgsize_t **bcConnList;
    int vcType;
    cgsize_t **connectOffSet;
    int isCartesian3D, isCartesian2D;
    int cellStartIndex;

    int nIFaceRegion;
    int *nIFaceElem;
    int *interFaceType;
    int *interFaceElementType, *interFaceGridLocation;
    cgsize_t **interFaceConnList;
    string   *interFaceName;

public:
    string *GetBCName() const { return boundaryName; }
    void   SetBCName(string *nameIn) { boundaryName = nameIn; }

    string GetVCName() { return volumeName; }
    void   SetVCName(string nameIn) { volumeName = nameIn; }

    int  GetVCType() { return vcType; }
    void SetVCType(int vcTypeIn) { vcType = vcTypeIn; }

    cgsize_t **GetBoundaryElementConnectionList() const { return bcConnList; }
    void     SetBoundaryElementConnectionList(cgsize_t **bcConnList) { this->bcConnList = bcConnList; }

    int  *GetBCType() const { return bcType; }
    void SetBCType(int *bcType) { this->bcType = bcType; }

    int  *GetBoundaryGridLocation() const { return bcGridLocation; }
    void SetBoundaryGridLocation(int *bcGridLocation) { this->bcGridLocation = bcGridLocation; }

    int  *GetBoundaryElementType() const { return boundaryElementType; }
    void SetBoundaryElementType(int *boundaryElementType) { this->boundaryElementType = boundaryElementType; }

    int  *GetNumberOfBCElements() const { return nBCElem; }
    void SetNumberOfBCElements(int *nBCElem) { this->nBCElem = nBCElem; }

    cgsize_t **GetElementConnectionList() const { return connList; }
    cgsize_t **GetMixedConnectionIndex() const { return mixedConnID; }
    int      GetNumberOfSections() const { return nSections; }
    void     SetNumberOfSections(int nSections) { this->nSections = nSections; }

    int  GetNumberOfBCRegions() const { return nBCRegions; }
    void SetNumberOfBCRegions(int nBCRegions) { this->nBCRegions = nBCRegions; }

    cgsize_t *GetIndexOfStart() const { return iStart; }
    cgsize_t *GetIndexOfEnd() const { return iEnd; }
    cgsize_t *GetElementDataSize() const { return elementDataSize; }

    cgsize_t **GetConnectOffSet() const { return connectOffSet; }

    void SetIsCartesian2D(int isCartesian2D) { this->isCartesian2D = isCartesian2D; }
    void SetIsCartesian3D(int isCartesian3D) { this->isCartesian3D = isCartesian3D; }
    void SetCellStartIndex(int cellStartIndex) { this->cellStartIndex = cellStartIndex; }
    int GetIsCartesian2D() { return isCartesian2D; }
    int GetIsCartesian3D() { return isCartesian3D; }
    int GetCellStartIndex() { return cellStartIndex; }

    int  *GetNumberOfInterFaceElements() { return nIFaceElem; }
    void SetNumberOfInterFaceElements(int *nIFaceElem) { this->nIFaceElem = nIFaceElem; }

    int  *GetInterFaceElementType() { return interFaceElementType; }
    void SetInterFaceElementType(int *interFaceElementType) { this->interFaceElementType = interFaceElementType; }

    int  *GetInterFaceGridLocation() { return interFaceGridLocation; }
    void SetInterFaceGridLocation(int *interFaceGridLocation) { this->interFaceGridLocation = interFaceGridLocation; }

    int  *GetInterFaceType() { return interFaceType; }
    void SetInterFaceType(int *interFaceType) { this->interFaceType = interFaceType; }

    cgsize_t **GetInterFaceElementConnectionList() { return interFaceConnList; }
    void     SetInterFaceElementConnectionList(cgsize_t **interFaceConnList) { this->interFaceConnList = interFaceConnList; }

    string *GetIFaceName() const { return interFaceName; }
    void   SetIFaceName(string *nameIn) { interFaceName = nameIn; }

    int *GetElementType() const { return elementType; }
public:
    void CreateElement(int nSections);
    void CreateBC(int nBCRegions);
    void CreateInterFace(int nIFaceRegions);
    void PreProcess(Base_Grid_Conn *gConn, BaseElement *baseElem);
    void PostProcess(Base_Grid_Conn *gConn);
    void GetISection(cgsize_t iCell, int &isection);

    bool CheckBCFace(int *bc_face, int nfacept, set<cgsize_t> &bc_vertex);

    void ComputeGridConn(Base_Grid_Conn *gConn, BaseElement *baseElem);
    void GetElem(cgsize_t *elem_source, cgsize_t iCell, int ptnum, cgsize_t *&elem);
    void GetBCFaceMixed(cgsize_t *elem_source, cgsize_t iCell, int isection, cgsize_t *&elem);
    void GetElem(BaseElement *base_elem, cgsize_t iCell, int &etype, cgsize_t *&elem);
    void GetFace(BaseElement *base_elem, cgsize_t iCell, int &facept_num, cgsize_t *&bcface);
    void GetFace(Base_Grid_Conn *gconn, cgsize_t *elem, cgsize_t iCell, CElemFace *elemface);
    void GetFace(cgsize_t *elem, int loc, int *bc_face, int &nfacept, CElemFace *elemface);

    //! Using for NGON_n and NFACE_n type of sections to get 3D Cartesian grid face information of face2node and face2cell.
    void GetFaceInformationofFace2Node(Base_Grid_Conn *gConn, int sectionID);
    void GetFaceInformationofFace2Cell(Base_Grid_Conn *gConn, int sectionID);

    //! Construct element face for 2D Cartesian grid.
    void GetCartesian2DFace(Base_Grid_Conn *gConn, cgsize_t *elem, cgsize_t iCell, int nPoint);
    void GetCartesian2DBCFace(cgsize_t iCell, int sectionID, int &facePtNum, cgsize_t *&bcFace);
    void ProcessBCFaceCartesian2D(Base_Grid_Conn *gConn, cgsize_t iElement, int sectionID, BCType_t bcType, const string &BCName);

    void ComputeElemFace(Base_Grid_Conn *gConn, BaseElement *baseElem);
    void ComputeBCFace(Base_Grid_Conn *gConn, BaseElement *baseElem);
    void ComputeCellToNode(Base_Grid_Conn *gConn, BaseElement *baseElem);
    bool MarkBCFace(Base_Grid_Conn *gconn, cgsize_t *bcface, int facept_num, BCType_t bctype, const string &BCName);
    void MarkBCInterface(Base_Grid_Conn *gconn);

    void GetBCFace(BaseElement *base_elem, cgsize_t iCell, int &facept_num, cgsize_t *&bcface);
    void ProcessBCFace(Base_Grid_Conn *gconn, BaseElement *base_elem, cgsize_t iElement, BCType_t bctype, const string &BCName);
    void ProcessBCVertex(Base_Grid_Conn *gconn, BaseElement *base_elem, int bc_elem_set_type, int nBCElem, cgsize_t *bc_conn_list, int bctype, const string &BCName);
    void ProcessBCVertex(Base_Grid_Conn *gconn, BaseElement *base_elem, set<cgsize_t> &bc_vertex, int bctype, const string &BCName);
};

class CGNS_BCFace
{
public:
    CGNS_BCFace(int faceID, int bcType, const string &bcName)
    {
        this->faceID = faceID;
        this->bcType = bcType;
        this->bcName = bcName;
    }

private:
    int faceID;
    int bcType;
    string bcName;
    //int bcID;

public:
    int GetFaceID() const { return this->faceID; }
    int GetBCType() const { return this->bcType; }
    const string & GetBCName() const { return this->bcName; }

    void SetFaceID(int faceID) { this->faceID = faceID; }
    void SetBCType(int bcType) { this->bcType = bcType; }
    void SetBCName(const string &bcName) { this->bcName = bcName; }

    CGNS_BCFace & operator = (const CGNS_BCFace &rhs)
    {
        if (this == &rhs) return *this;

        this->faceID = rhs.GetFaceID();
        this->bcType = rhs.GetBCType();
        this->bcName = rhs.GetBCName();

        return *this;
    }
};

class Base_Grid_Conn
{
public:
    Base_Grid_Conn();
    ~Base_Grid_Conn();
public:
    void SetNTotalFace(cgsize_t nTotalFace) { this->nTotalFace = nTotalFace; }
    cgsize_t GetNTotalFace() const { return nTotalFace; }

    void SetNBoundFace(cgsize_t nBoundFace) { this->nBoundFace = nBoundFace; }
    cgsize_t GetNBoundFace() const { return nBoundFace; }

    void SetNTotalCell(cgsize_t nTotalCell) { this->nTotalCell = nTotalCell; }
    cgsize_t GetNTotalCell() const { return nTotalCell; }

    void SetNTotalNode(cgsize_t nTotalNode) { this->nTotalNode = nTotalNode; }
    cgsize_t GetNTotalNode() const { return nTotalNode; }

    vector< vector<cgsize_t> > * GetFace2Node() const { return face2node; }
    vector<CGNS_BCFace*> * GetBCFace() const { return this->bcFace; }
    vector<cgsize_t> * GetLeftCellOfFace() const { return left_cell_of_face; }
    vector<cgsize_t> * GetRightCellOfFace() const { return right_cell_of_face; }
    vector<int> * GetFaceLocation() const { return face_location; }
    vector<int> * GetNodeNumberOfEachFace() const { return node_number_of_each_face; }
    set<MultiIndex> * GetFaceSet() const { return faceset; }
    vector<int> * GetNodeNumberOfEachCell() const { return node_number_of_each_cell; }
    vector< vector<cgsize_t> > * GetCell2Node() const { return cell2node; }

    void FreeFaceSet();
    void FreeFace2Node();
    void FreeCell2Node();
    void FreeCellOfFace();
    void FreeNodeNumberOfEachFace();
    void FreeNodeNumberOfEachCell();

private:
    cgsize_t nTotalNode, nTotalCell, nBoundFace, nTotalFace;
    set<MultiIndex> *faceset;
    vector<int> *face_location;
    vector<CGNS_BCFace *> *bcFace;
    vector<cgsize_t> *left_cell_of_face, *right_cell_of_face;
    vector<int> *node_number_of_each_face;
    vector< vector<cgsize_t> > *face2node;
    vector<int> *node_number_of_each_cell;
    vector< vector<cgsize_t> > *cell2node;
};

class CGNSFactory
{
private:
    int nBlocks;
    CGNSBase       **baseCGNS;
    Base_Grid_Conn **baseConn;
    BaseElement    *baseElem;
    BaseBCType     *fantasyBC;
    RawGrid        **rawGrid;

public:
    CGNSFactory(int nBlocks);
    ~CGNSFactory();
public:
    Base_Grid_Conn **GetBaseGridConnection() const { return baseConn; }
    Base_Grid_Conn *GetBaseGridConnection(int iZone) const { return baseConn[iZone]; }
    BaseElement *GetBaseElement() const { return baseElem; }
    BaseBCType *GetBaseBCType() const { return fantasyBC; }
    RawGrid  *GetRawGrid(int iZone) { return rawGrid[iZone]; }
    CGNSBase *GetCGNSBase(int iZone = 0) { return baseCGNS[iZone]; }
    int GetNumberofBlocks() const { return nBlocks; }
    void ConvertGrid2Fantasy();
    void DumpSelfDefinePartInformation();
    void SetNumberofBlocks(int nBlocksIn) { this->nBlocks = nBlocksIn; }
private:
    void ComputeGridConn(CGNSBase *cgnsGrid, Base_Grid_Conn *gConn);
    void ReLableSelfDefineBC();
};

void CGNS2UnsGrid(CGNSFactory *factoryCGNS, int iZone, UnstructGrid *grid);

}