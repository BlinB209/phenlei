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
//! @file      ComputationalGrid.h
//! @brief     Construct grid interface information, used in Parallel unstructured grid refining. 
//! @author    Baka, Bell.

#pragma once
#include "TypeDefine.h"
#include "DataContainer.h"
#include "Geo_Interface.h"
#include "Region.h"

namespace PHSPACE
{
class InterfaceLinkStructure;

//! @brief ZoneInterface store base information of grid.
class ZoneInterface
{
public:
    ZoneInterface();
    ~ZoneInterface();

private:
    //! Grid index.
    int zoneIndex;

    //! Nodes number of grid.
    int numberOfNodes;

    //! Min-max box of the nodes coordinates.
    RDouble *pmin, *pmax;
    RDouble disMin, disMax;
    
    //! Coordinates.
    RDouble *x, *y, *z;

    //! Face topo information.
    PHVectorInt2D faceNodeIndexContainer;

public:
    void SetNumberOfNodes(int numIn) { this->numberOfNodes = numIn; }
    int  GetNumberOfNodes() const { return this->numberOfNodes; }

    uint_t GetNumberOfFaces() const { return this->faceNodeIndexContainer.size(); }

    void SetZoneIndex(int indexIn) { this->zoneIndex = indexIn; }
    int  GetZoneIndex() const { return this->zoneIndex; }

    void SetMinBox(RDouble *boxIn) { this->pmin = boxIn; }
    void SetMaxBox(RDouble *boxIn) { this->pmax = boxIn; }
    RDouble * GetMinBox() const { return this->pmin; }
    RDouble * GetMaxBox() const { return this->pmax; }

    void SetMinMaxDistance(RDouble &dismin, RDouble &dismax) { this->disMin = dismin; this->disMax = dismax; }
    void GetMinMaxDistance(RDouble &dismin, RDouble &dismax) { dismin = this->disMin; dismax = this->disMax; }

    void SetX(RDouble *x) { this->x = x; }
    void SetY(RDouble *y) { this->y = y; }
    void SetZ(RDouble *z) { this->z = z; }

    RDouble * GetX() const { return x; }
    RDouble * GetY() const { return y; }
    RDouble * GetZ() const { return z; }

    PHVectorInt2D & GetFaceNodeIndexContainer() { return this->faceNodeIndexContainer; }

    //! Compress and decompress data, used to communicate grid information.
    void Encode(DataContainer *&dataContainer);
    void Decode(DataContainer *&dataContainer);
};

//! @brief ComputationalGridProcess construct grid interface information, used in Parallel unstructured grid refining. 
class ComputationalGridProcess
{
public:
    ComputationalGridProcess(Region *region);
    ~ComputationalGridProcess();
protected:
    //! Number of grid.
    int numberOfZones;

    //! Interface information.
    InterfaceLinkStructure *interfaceLinkStructure;

    //! Interface information of neighbor grid.
    vector< ZoneInterface * > neighborZoneInterfaces;
    vector< ZoneInterface * > zoneInterfacesTopologyList;

    //! current grid.
    int currentZoneID;
    UnstructGrid *currentGrid;

    //! boundary information of face.
    vector< PHVectorInt1D >    boundaryTypeVectorList;
    vector< vector< string > > boundaryNameVectorList;

    //! Index map between loacl interface and global interface.
    vector< PHVectorInt1D > newLocalInterfaceIndexToGlobalInterfaceIndexMappingList;
    
    //! Left cell index of boundary face.
    vector< PHVectorInt1D > newBoundaryLeftCellIndexList;

    //! Face topo information.
    vector< PHVectorInt1D > newLeftCellIndexVectorList, newRightCellIndexVectorList;
    vector< PHVectorInt2D > faceNodeIndexVectorList;

 public:
    //! Initialize: allocate memory.
    void Initialize(Grid **gridContainer, int nZones);

    //! Construct grid interface information. 
    void PostProcessMultiZoneComputationalGrids();

protected:
    //! Construct neighbor grid interface topo information
    void CreateNeighborZoneInterfacesTopology();

    //! Construct grid interface topo information
    void BuildZoneInterfacesOfThisProcessor();
    void BuildZoneInterfacesTopology(UnstructGrid *unstructuredGrid, int iZone);

    //! Init neighbor grid interfaces.
    void InitZoneInterfacesOfThisProcessor();
    void AddZoneInterfacesTopology(ZoneInterface *zoneInterfaceIn);

    //! Communicate interface information.
    void SwapZoneInterfaceTopologyToNeighborZones();
    void SwapZoneInterfaceTopologyToNeighborZonesNonBlocking();
    void SwapZoneInterfaceTopologyToNeighborZonesBlocking();
    void SwapInterfaceTopologyBetweenTwoNeighboringZone(int zoneIndex, int sendProcessorIndex, int receiveProcessorIndex);
    void CompressZoneInterfaceIntoDataContainer(DataContainer *&dataContainer, int globalZoneIndex);
    void AnasysZoneInterfaceFromDataContainer(DataContainer *&dataContainer);
    void DecodeZoneInterfacesTopology(DataContainer *&dataContainer);

    //! Get the length of neighbors.
    void SwapZoneInterfaceLengthToNeighborZones(vector<CharVecSizeType> &neighborLength);

    //! Create interfaceLink structure.
    void CreateInterfaceLinkStructure();

    //! Get neighbor grid coordinates information.
    void GetNeighborZoneBoundingBox   (RDouble *pmin, RDouble *pmax);
    void GetNeighborZoneMinMaxDistance(RDouble &mindis, RDouble &maxdis);

    //! Covert NO_BOUNDARY_CONDITION to INTERFACE.
    void PreprocessBoundaryCondition();
    void ChangeBoundaryConditionTypeIfNeeded();

    //! Construct grid interface information map between local grid and neighbor frid.
    void GenerateLocalGlobalInterfaceMappingByNeighborZones();
    void GenerateLocalGlobalInterfaceMappingByNeighborZones(int iZone);

    //! Build interface information for each grid.
    void CreateInterfaceInformation();
    void CreateLocalInterfaceToBoundaryFaceIndexesMapping();

    //! Reconstruct interface topology, only used to QUAD_4 face.
    void ReconstructInterfaceTopology();

    //! Reconstruct face topology, after node index changed.
    void ReGenerateLocalGlobalInterfaceMapping();
    void InitializeNewLocalGlobalInterfaceMapping();
    void ModifyFaceNodesIndexes               (InterfaceLinkStructure *interfaceLinkStructure);
    void ModifyBoundaryInformation            (InterfaceLinkStructure *interfaceLinkStructure);
    void SetBoundaryFaceNodeIndexVector       (InterfaceLinkStructure *interfaceLinkStructure);
    void SetNewFaceToCellIndexesContainer     (InterfaceLinkStructure *interfaceLinkStructure);
    void ResetNumberOfBoundaryCondition       (InterfaceLinkStructure *interfaceLinkStructure, PHVectorInt1D &numberOfChildFacesOfInterface);
    void ProcessGlobalToLocalInterfaceTerm    (InterfaceLinkStructure *interfaceLinkStructure);
    void DoFurtherWorkUsingBoundaryTypeVector (InterfaceLinkStructure *interfaceLinkStructure);
    void ReGenerateLocalGlobalInterfaceMapping(InterfaceLinkStructure *interfaceLinkStructure);

    //! Update interface information after reconstruct face topology.
    void UpdateLocalGlobalInterfaceMapping();
    void UpdateOtherTopologyTerm();
    void UpdateOtherTopologyTerm(InterfaceLinkStructure *interfaceLinkStructure);
    void UpdateOtherBoundaryTopologyTerm();
    void UpdateOtherGridFaceTopologyTerm();

    //! Search neighbor grid index and interface index of in neighbor grid of each interface in current grid.
    void MatchInterfaceTopology();

private:
    Grid ** GetGridContainer() { return gridContainer; }
    Grid  * GetGrid(int iZone) { return gridContainer[iZone]; }
    InterfaceLinkStructure * GetInterfaceLinkStructure() { return interfaceLinkStructure; }
    vector< ZoneInterface* > & GetNeighborZoneInterfaces() { return neighborZoneInterfaces; }

private:
    Region *region;
    Grid **gridContainer;

};

//! Return if interface is exist in grid.
bool ExistInterfaces(UnstructGrid *unstructuredGrid);

//! Return node coordinate in face.
void GetFaceCoordinateList(PHVectorInt1D &faceNodeIndexes, int numberOfPoints, RDouble *xList, RDouble *yList, RDouble *zList, RDouble *x, RDouble *y, RDouble *z);

//! Construct ADTNode by node coordinate in face, then insert in ADTree.
int  GetCoordinateIndex(DataStruct_AdtTree< int, RDouble > *adtTree, RDouble *coordinate, RDouble *minWindow, RDouble *maxWindow, int &pointCount, int &pointIndex);
void GetCoordinateIndexList(InterfaceLinkStructure *interfaceLinkStructure, int &pointCount, RDouble *xList, RDouble *yList, RDouble *zList, int numberOfPoints, PHVectorInt1D &pointIndexes);

}