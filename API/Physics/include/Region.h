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
//! @file      Region.h
//! @brief     It defines the region, that is defined as the computational domain.
//!            It is used to manage the geometry, solvers, fields and so on.
//!            It is one of the most important classes in PHengLEI.
//! @author    Bell, He Xin.

#pragma once
#include "Pointer.h"
#include "TK_Exit.h"
#include "Geo_StructBC.h"
#include "MultigridFactory.h"
#include "MixGrid.h"
#include "IO_GridReader.h"
namespace PHSPACE
{
//! @brief Region class is used to manage the geometry, solvers, fields and so on.
//!        It is one of the most important classes in PHengLEI.
//!        In PHengLEI, the whole computational domain is divided in to many sub-domains(zones) using
//!        different partition methods on structured/unstructured grid. Each sub-domain is defined as 'zone'.
//!        'Grid', which are computational mesh include point/face/cell, are built on each 'zone'.
//!        The zones belong to the current processor are defined as 'Region'.
//!        The relationship of these object:\n
//!        structured/unstructured mesh --> grid --> zone --> region.
class Zone;
class Region
{
public:
    LIB_EXPORT Region();
    LIB_EXPORT ~Region();

private:
    //! The main data structure, which include the 'Zones' (Referenced as Zone class).
    //! The size of zoneContainer is equal to the number of global partition zones.
    //! IMPORTANT: only the zones that belong to the current processor are NON-ZERO Entity,
    //! the zones that belong to other processors are set to be ZERO in this container.
    vector <Zone *> *zoneContainer;
public:
    //! Return the iZone-th zone.
    Zone * GetZone(int iZone);

    //! If the iZone-th zone belong to current processor.\n
    //! @return
    //! -# <em>false</em>:  the zone with globalZoneIndex dose not belong to current processor.\n
    //! -# <em>true</em>:   the zone with globalZoneIndex belongs to current processor.
    bool ZoneBelongToCurrentProcessor(int globalZoneIndex);

    //! Return the global zone index of 'localZoneIndex'.
    int  GetProcessLocalZoneIndexToGlobalZoneIndex(int localZoneIndex);

    void SetIBlock(int iBlock) { this->iBlock = iBlock; }
    //! Read grid from file.
    LIB_EXPORT void ReadGrid();

    LIB_EXPORT void InitGridForOversetConfig();

    //! Return the local zone index of 'globalZoneIndex'.
    LIB_EXPORT int  GetProcessGlobalZoneIndexToLocalZoneIndex(int globalZoneIndex);
    
    //! Initialize the geometry (grid actually), the initialization tasks include:\n
    //!   1: Set up Multi-Block relationship, that is the MPI communication relationship.\n
    //!   2: Calculate the grid metrics, such as volume, area, center and so on.\n
    //!   3: Agglomerate coarse grids, if multi-grid method is used.\n
    //!   4: Calculate the wall distance, that is the nearest distance between each cell center and solid wall.\n
    //!   5: Change boundary condition if necessary.\n
    //!   6: Initialize dynamic(moving) grid.
    LIB_EXPORT void InitGeometry();


    //! If need, change the default boundary condition to other type according to the 
    //! input parameter file named 'BCMap.hypara'. Details reference to the file.
    LIB_EXPORT void ChangeBoundaryCondition();

    //! Communicate the grid metrics on the interfaces.\n
    //! And calculate the coefficient matrix if LSQ method is used.
    LIB_EXPORT void SwapGeometry();
    
    //! Dump out some informations after whole of the simulation.
    LIB_EXPORT void PostSimulation();

    //**************************************************************
    //! The following functions should not be called unless special
    //! reasons.
    //**************************************************************

    //! Read ordinary grid from file. this function is not allowed to call. 
    //! If the grid wanted to be read, call the 'ReadGrid' function above, Please!
    LIB_EXPORT void ReadOrdinaryGridP2PMode(PHVectorString1D & gridGroupNameList);
    //! Read ordinary grid from file. this function is not allowed to call. 
    LIB_EXPORT void ReadOrdinaryGridCSMode(PHVectorString1D & gridGroupNameList);
    //! Read ordinary grid from file. this function is not allowed to call. 
    //! If the grid wanted to be read, call the 'ReadGrid' function above, Please!
    LIB_EXPORT void ReadOrdinaryGrid(string fileName, int fileMode);

    //! Read ordinary grid from file. this function is not allowed to call. 
    LIB_EXPORT void ReadOrdinaryGrid(string fileName);
    LIB_EXPORT void ReadOrdinaryGrid(PHVectorString1D & gridGroupNameList);

    //! Read overset / overlapping grid. this function is not allowed to call. 
    //! If the grid wanted to be read, call the 'ReadGrid' function above, Please!
    LIB_EXPORT void ReadOversetGrid();

    LIB_EXPORT void SetOversetGrid();

    LIB_EXPORT void ReviseInterZoneIndexForOverSetGrids(int iZone);

    //! Build structured/unstructured grid if HYBRID Solver is used.
    LIB_EXPORT void BuildAndDumpMixedGrid();

    LIB_EXPORT void WriteFaceBoundaryName(Grid **grids, int nZones, const string &gridFileName);

    //! Set up Multi-Block relationship, that is the MPI communication relationship.\n
    LIB_EXPORT void SetupMultiBlockRelationship();

    //! Agglomerate coarse grids and set up multi-grid relationship, if multi-grid method is used.
    LIB_EXPORT void MultiGridSetup();

    //! Initialize the grid,include:\n
    //!   1: Set the number of level of Multi-Grid.\n
    //!   2: Calculate the grid metrics.\n
    //!   3: Initialize moving grid.
    LIB_EXPORT void InitGrids();

#ifdef AI_RandomForestRegressor
    void DumpCellCenterMetrics();
#endif

    void DumpWallFaceCenter();

    //! Update all of parameters in the running data base of each zones according to the global_para.
    LIB_EXPORT void UpdateAllZonesPara();

    //! Read single grid from file stream, create zone, and then insert in to zoneContainer.
    LIB_EXPORT void ReadSingleGrid(fstream &file, int iZone);

    //! Read layout information from file.\n
    //! 'Lay out' is defined as the multi-block information on both structured/unstructured grid, include:\n
    //!   1: The zone index of each global zone.\n
    //!   2: The assigned processor index of each global zone.\n
    //!   3: The grid type of each global zone. 'grid type': STRUCTURED grid or UNSTRUCTURED grid.
    LIB_EXPORT void PreprocessLayoutInformation(vector< string > &filenamelist);

    //! Fix the grid type according to the lay out information read from file.
    LIB_EXPORT void ClassifyGridSystem();

    //! Monitor the data according to location. 
    LIB_EXPORT void DataMonitor();

    void InitVariableWallTemperature();

    void ReSetBoundaryConditionByGlobalBC();

    void CommunicateCellIBlank();

    void AddZone(Zone *zone);

    void CommGridInfoOnlyforStructHighOrder(int mission);

    void InitOriginalGrid();
private:
    void DumpProcessorGridNumber();
    void CreateGridFileNameList(vector< string > &filenamelist);
    void CreateGridFileNameListCSMode(vector< string > &filenamelist);
    void CreateGridFileNameListP2PMode(vector< string > &filenamelist);
    void GetZoneLayoutInformation(vector< string > &filenamelist);
    void GetZoneLayoutInformation(const string &filename);
    void ReadNumberOfZones(fstream &file);
    void InitZoneLayoutInformation(vector< string > &filenamelist);
    void InitZoneLayoutInformation(const string &filename);
    void AllocateZoneLayoutInformation();

    void MultiBlocksInterface();

    //! build the connection in multiblocks.
    void MultiBlocksInterpoint();

    void MultiBlocksOverInterface();

    void InitZoneNeighborsInformation();

    //! build the zone Neighbor information for point.
    void InitZoneNeighborsForPointInformation();

    void SwapZoneNeighborsInformation();

    //! swap the zone Neighbor information for point.
    void SwapZoneNeighborsForPointInformation();

    void SwapNeighborsSendContent();

    //! swap the zone Neighbor information for send.
    void SwapNeighborsSendContentForPoint();

    void SwapNeighborsSendContentNonblocking();

    void BroadcastZoneNeighborInformation();

    //! swap the zone Neighbor information for point.
    void BroadcastZoneNeighborForPointInformation();

    void BroadcastZoneNeighborInformationToEachProcessor();
    void BroadcastZoneNeighborInformationToNeighborZones();
    void SwapZoneConnectivityBetweenNeighborZones();

    //! swap the zone Connectivity information for point in neighbor zones.
    void SwapZoneConnectivityForPointBetweenNeighborZones();

    //! New method to bcast zone neighbor information.\n
    //! step1: collecting from each other processor by server.
    //! step2: bcasting to each other processor by server.
    void BroadcastZoneNeighborInformationCollectBcast();
    void ServerCollectionNeighborInfor();
    void ServerBcastNeighborInfor();
    void CleanNeighborInfor();

    void GetProcessorZones(int iZone, ZoneNeighbor * zoneNeighbor, vector< int > * neighborProcessors, vector< set< int > > * processorZones);

    void CheckMeshMultigridByFTS();
    void CoarseGrids();
    void ShowCoarseGrids();

    void SetNGrids();
    void ComputeMetrics(int level = 0);
    void InitMovingGrids(int level = 0);
    void InitCellBoundaryType(int level = 0);

    //! Global grid information static, for structured and unstructured.
    void GlobalGridStatic(int level = 0);
    //! Compute global min&&max cell number in all zones.
    void ComputeGlobalCellSize(int level = 0);
    //! Compute global total number of cells and total volume for vencat limiter of type 2 && 3.
    void ComputeGlobalGridSize(int level = 0);
    void buildGlobalCellsArrary(int level = 0);

#ifdef USE_GMRESSOLVER
    //! GMRESParallel
    //! get shift of cells in the local zone with respect to the global mesh for GMRES
    void GetCurrentShiftCellSizeToGlobalSize(int level = 0);
 #endif   
 
    //! MixGrid
    void ComputeGlobalGridCellSize(int level = 0);

    void CommunicateCellCenterData();
    void ComputeMetrics_BC();
    void ComputeMetricsandJacobinWithBCforStructHighOrder();
    void UpdateDiffOrderforStructHighOrder();
    void UpdateMetricsWithBCforStructHighOrder();
    void UpdateJacobianWithBCforStructHighOrder();

    void ComputeWeight();
    void TestReconstruction();

    int *GetNode2All(int iZone);
    int *GetFace2All(int iZone);
    int *GetCell2All(int iZone);

    void ReadStructuredOversetGrid();
    void ReadGridTopologyData(const string & fileName);
    void DeleteLabelContainer();

    void ReadUnstructuredOversetGrid();
    void ReadGrid(fstream &file);
    void ReadOversetGrid(const string &gridfile);
    void ReadSingleOversetGrid(fstream &file, int iZone);
    void AllocateOversetGrid();
    void AddUnStructInfoToStructGrid();

    void CreateZonesFromGridReader(IO_GridReader * gridReader);
    void CreateZonesFromGridGroups(IO_GridReader * gridReader);
    void ExcangeProcessorIDForGridProcessor();
    void CreateZoneAndGrid(int iZone);
    void RecoverGridFromFile(fstream &file, int iZone);
    void InterpretDataToGrid(DataContainer *cdata, int iZone);
    //! Check grid validity.
    void GridValidityCheck();

private:
    int numberOfProcessorStructgridUsed;
    int numberOfZoneInLocalProcess;
    int nGrids;
    int iBlock;
    int numberOfTotalZones;
    int zoneStart;
};

void ResetAllZonesOversetGrid();
#include "Region.hxx"

}
