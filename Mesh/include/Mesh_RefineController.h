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
//! @file      Mesh_RefineController.h
//! @brief     Single unstructured grid refine.
//! @author    Baka, Bell.

#pragma once
#include "Mesh_CellTopo.h"
#include "Geo_PointFactory.h"
#include "Mesh_FaceTopo.h"

namespace PHSPACE
{

//! @brief Mesh_RefineController achieve single unstructured grid refine.
class Mesh_RefineController
{
private:
    //! Current grid index.
    int zoneIndex;

    //! Node number.
    int nTotalNode;

    //! Boundary face number
    int nBoundaryFace;

    //! Cell number.
    int oldNumberOfCell;

    //! Cell inoformation of grid.
    Mesh_CellTopo *cellTopo;

    //! Face inoformation of grid.
    Mesh_FaceTopo *faceTopo;

    //! Node inoformation of grid.
    Geo_PointFactory *nodeTopo;

    //! Interface inoformation of grid.
    InterfaceInfo *interfaceInformation;

    //! Face inoformation of refined grid.
    Mesh_FaceTopo *computationalFaceTopo;

    //! If current grid has been refined.
    bool hasBeenRefined;

public:
    Mesh_RefineController();
    ~Mesh_RefineController();

public:
    //! Return cell, face and node information of grid.
    Mesh_CellTopo    * GetCellFactory() { return this->cellTopo; }
    Mesh_FaceTopo    * GetFaceFactory() { return this->faceTopo; }
    Geo_PointFactory * GetNodeFactory() { return this->nodeTopo; }

    //! Set interface information of grid.
    void SetInterfaceInformation(InterfaceInfo *interfaceInformationIn) { this->interfaceInformation = interfaceInformationIn; }

    //! Construct grid topo information, node to cell, cell to cell.
    void Initialize();
    void ConstructNode2Cell();
    void ConstructCell2Cell();

    //! Set refine type of each cell.
    void InitAnisoRefineType();
    void FixAnisoRefineType();

    //! Grid refine.
    void RefineGrid();

    //! Split all cell.
    void SplitAllElements();
    void SplitElement(int cellIndex);
    void ReactivateChildElement(int cellIndex);
    void AddChildElement(int parentElementIndex);

    //! Split all face.
    void SplitElementFace(int cellIndex);
    void SplitFace(PHVectorInt1D &iCompositeFace2Node, int compositeFaceType);

    //! Construct grid information after refined.
    void ConstructBasicData();
    void ScanAllElementsToReconstructFaceDataStructure();

    //! Reorder cell index after refined.
    void GenerateCellIndex();

    //! Construct cell and face index map between refined cells and global cells.
    void GenerateGlobalToComputationalCellFaceIndex();
    void GenerateGlobalToComputationalCellFaceIndex(int iFace, int &iComputationalFaceCount, int &iBoundaryFace);
    void GenerateGlobalToComputationalCellFaceIndexNoParentCase(int iFace, int &iComputationalFaceCount, int &iBoundaryFace);
    void GenerateGlobalToComputationalCellFaceIndexWithParentCase(int iFace, int &iComputationalFaceCount, int &iBoundaryFace);
    void GenerateComputationalCellFaceInformation();
    void GenerateComputationalElementFaceInformationFirstPart();
    void GenerateComputationalElementFaceInformationSecondPart();

    //! Construct face topo information after refined, include face to cell and face to node.
    void ProcessBoundaryFace(int iFaceForComputation);
    void ProcessGereralFace(int iFaceForComputation);
    void ProcessHangingFace(int iFaceForComputation);
    void GenerateComputationalFaceNodeIndex();

    //! Construct middle points of each edge.
    void GenerateEdgesWithMiddlePoints(int iCompositeFace);

    //! Set grid information.
    void SetZoneIndex(int zoneIndexIn) { this->zoneIndex = zoneIndexIn; }
    void SetNTotalNode(int nTotalNodeIn) { this->nTotalNode = nTotalNodeIn; }
    void SetNBoundaryFace(int nBoundaryFacesIn) { this->nBoundaryFace = nBoundaryFacesIn; }
    void SetOldNumberOfCell(int oldNumberOfCellIn) { this->oldNumberOfCell = oldNumberOfCellIn; }
    void SetRefineStatus(bool statusIn) { this->hasBeenRefined = statusIn; }

    //! Return grid information.
    int GetNBoundaryFace() { return this->nBoundaryFace; }
    int GetOldNumberOfCell() { return this->oldNumberOfCell; }

    //! Generate refined grid.
    void GenerateComputationalGrid(Grid *gridIn);
    void GenerateComputationalUnstructuredGrid(UnstructGrid *grid);
    void ReorderFaceToNodeIndexes(UnstructGrid *grid);
    void ReorderFaceToCellIndexes(UnstructGrid *grid);
    void GenerateCellToNode(UnstructGrid *grid);

    //! Communicate interface information, when parallel grid refining.
    void CompressRefineType(DataContainer *dataContainer, int neighborZoneIndex);
    void DecompressRefineType(DataContainer *dataContainer, int neighborZoneIndex);
    bool IsNeedInfectInterfaceIsotropicRefineType();
    void GetSourceCellIndex(int iFace, int &sourceIndex);
    void GetTargetIndex(int iFace, int &targetIndex);
};

}