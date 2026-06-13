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
//! @file      Pre_FluentToPHengLEI_Unstrut.h
//! @brief     Grid conversion from fluent format(.cas) to PHengLEI(.fts).
//! @author    He Xin.

#pragma once
#include "Pre_GridConversion.h"
#include "Pre_GridBase.h"

using namespace std;

namespace PHSPACE
{
class GridFormat_Fluent;
//! @brief Pre_FluentToPHengLEI_Unstrut class defines the method of 
//! Grid conversion from fluent format(.cas) to PHengLEI(.fts).
class Pre_FluentToPHengLEI_Unstruct : public Pre_GridConversion
{
public:
    //! @param[in] Grid file name of fluent format grid.
    LIB_EXPORT Pre_FluentToPHengLEI_Unstruct(const string &gridFileName);

    LIB_EXPORT ~Pre_FluentToPHengLEI_Unstruct();

    //! Write additional information besides grid file, e.g face BC file and BC name.
    //! @param[in] Output grid file name of PHengLEI.
    LIB_EXPORT void WriteAdditionalInformation(const string &out_grid_file);

private:
    void ReadGrid();
    void Conversion();
    void ReadFluent(GridFormat_Fluent *&fluent_factory, int &nBlocks, const string &casefile);
    void AnalysisFluentKeyWord(GridFormat_Fluent *fluent_factory, fstream &file, string &line);
    void ReadFluentNodes(RawGrid *rawgrid, int first_index, int last_index, int dimension, fstream &file);
    void ReadFluentFace(GridFormat_Fluent *fluent_factory, int zone_id, int first_index, int last_index, int bc_type, int face_type, fstream &file);
    void ReadFluentCell(int first_index, int last_index, fstream &file);
    void ReadFluentNodesBinary(RawGrid *rawgrid, int first_index, int last_index, int dimension, fstream &file);
    void ReadFluentFaceBinary(GridFormat_Fluent *fluent_factory, int zone_id, int first_index, int last_index, int bc_type, int face_type, fstream &file);
    void ReadFluentCellBinary(int first_index, int last_index, fstream &file);
    void ReadFluentEdge(int first_index, int last_index, fstream &file);
    void ReadFluentPeriodicface(int first_index, int last_index, fstream &file);
    void ReadFluentCellTree(int first_index, int last_index, fstream &file);
    void ReadFluentFaceTree(int first_index, int last_index, fstream &file);

private:
    GridFormat_Fluent *fluent_factory;
};

class FluentFacePatch
{
private:
    int zone_id, bc_type, face_type;
    int nTotalFace;
    vector<int> left_cell_of_face, right_cell_of_face;
    vector<vector<int> > face2node;
    string patchName;
public:
    int  GetZoneID() const { return zone_id; }
    void SetZoneID(int zone_id) { this->zone_id = zone_id; }
    int  GetFaceType() const { return face_type; }
    void SetFaceType(int face_type) { this->face_type = face_type; }
    int  GetBCType() const { return bc_type; }
    void SetBCType(int bc_type) { this->bc_type = bc_type; }
    int  GetNTotalFace() const { return nTotalFace; }
    void SetNTotalFace(int nTotalFace) { this->nTotalFace = nTotalFace; }
    void SetPatchName(const string &patchNameIn) { this->patchName = patchNameIn; }
    string GetPatchName() const { return this->patchName; }
public:
    vector<vector<int> > & GetFace2Node() { return face2node; }
    vector<int> & GetLeftCellOfFace()  { return left_cell_of_face; }
    vector<int> & GetRightCellOfFace() { return right_cell_of_face; }
};

class GridFormat_Fluent
{
public:
    GridFormat_Fluent();
    ~GridFormat_Fluent();
private:
    RawGrid *rawgrid;
    int dimension;
    int nTotalNode, nTotalCell, nTotalFace, nBoundFace;
    int nPatchs;
    vector<FluentFacePatch *> face_patchs;
    map<int, int> bcmap;
public:
    RawGrid * GetRawGrid() const { return rawgrid; }
    int  GetNTotalNode() const { return nTotalNode; }
    void SetNTotalNode(int nTotalNode) { this->nTotalNode = nTotalNode; }
    int  GetNTotalCell() const { return nTotalCell; }
    void SetNTotalCell(int nTotalCell) { this->nTotalCell = nTotalCell; }
    int  GetNTotalFace() const { return nTotalFace; }
    void SetNTotalFace(int nTotalFace) { this->nTotalFace = nTotalFace; }
    int  GetNBoundFace() { return nBoundFace; }
    int  GetDimension() const { return dimension; }
    void SetDimension(int dimension) { this->dimension = dimension; }
    vector<FluentFacePatch *> & GetFluentFacePatchS() { return face_patchs; }
public:
    void Process();
    void ReorderFacePatchs();
    void BuildBCMap();
    int  Fluent2FantasyBC(int bctype);
    void SetBoundaryName(int zoneIndex, const string &BCName);
    void ResetBCTypeOfSelfDefinePart(int zoneIndex);
    void DumpSelfDefinePartInformation();
    void Convert2UnsGrid(UnstructGrid *grid);
};

}