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
//! @file      Pre_HDF5File.h
//! @brief     It defines HDF5 file read and write method.
//! @author    Baka.

#pragma once
#include "IO_HDF5File.h"
#include "Region.h"

namespace PHSPACE
{

class IO_HDF5Write
{
private:
    Grid **grids;
    int nBlocks;
    string gridFileName;

public:
    //! Dump nBolcks grids to file.
    //! @param[in] gridFileName    Grid file name.
    //! @param[in] grids           The grids that need to be dumped.
    //! @param[in] nBlocks         The number of the grid zones.
    IO_HDF5Write(const string &gridFileName, Grid **grids, int nBlocks);
    ~IO_HDF5Write();

public:
    void Run();

private:
    hid_t CreateFile(const string &filename);
    void WriteVersionInfo(hid_t loc_id);
    void CreateTotalInfo();
    void WriteTotalInfo(hid_t loc_id);
    void WriteGridCoordinates(hid_t loc_id, Grid *gridIn);
    void WriteInterfaceInfo(hid_t loc_id, Grid *gridIn);
    void WriteLnkInfo(hid_t loc_id, Grid *gridIn);
    void WriteVolumeCondition(hid_t loc_id, Grid *gridIn);

    void WriteStructGrid(hid_t loc_id, Grid *gridIn);
    void WriteBCRegion(hid_t loc_id, Grid *gridIn);
    void WriteStrBCName(hid_t loc_id, Grid *gridIn);
    void WriteStrPartitionInfo(hid_t loc_id, Grid *gridIn);

    void WriteUnstructGrid(hid_t loc_id, Grid *gridIn, int iZone);
    void WriteFaceTopology(hid_t loc_id, Grid *gridIn);
    void WriteCellTopology(hid_t loc_id, Grid *gridIn);
    void WriteUnstrBCName(hid_t loc_id, Grid *gridIn);
    void WriteInterpointInfo(hid_t loc_id, Grid *gridIn);
    void WriteUnstrPartitionInfo(hid_t loc_id, Grid *gridIn);

private:
    //! The grids int current processor
    Grid **gridsLocal;
    int *nPartOfEachProcessor;

    int *block_proc;
    int *block_idx;
    int *block_type;
    int *file_index;
};

class IO_HDF5Read
{
private:
    PHVectorString1D gridNameList;
    double version;

public:
    //! Read grid from file.
    //! @param[in] gridFileName    Grid file name.
    //! @param[in] region_in       The class to manage the geometry.
    IO_HDF5Read(PHVectorString1D &gridNameListIn, Region *region_in);
    ~IO_HDF5Read();

public:
    bool Run();

private:
    void CheckNumberofFiles(int iGrid);
    void InitGridFileNameListAndFileProcessor();
    hid_t OpenH5File(const string &filename);

    void ReadVersionInfo(hid_t loc_id);
    void ReadTotalInfoOld(hid_t loc_id);
    void ReadNumberOfZonesOld(hid_t loc_id);    //! Useless.
    void ReadEachGridOld(hid_t loc_id);

    void ReadTotalInfo();
    void ReadNumberOfZones();
    void ReadEachGrid();

    void ReadGridCoordinates(hid_t loc_id, Grid *gridIn);
    void ReadInterfaceInfo(hid_t loc_id, Grid *gridIn);
    void ReadVolumeCondition(hid_t loc_id, Grid *gridIn);

    void ReadStructGrid(hid_t loc_id, int iZone);
    void ReadBCRegion(hid_t loc_id, Grid *gridIn);
    void ReadLnkInfo(hid_t loc_id, Grid *gridIn);
    void ReadLnk(hid_t loc_id, int iZone);
    void ReadStrBCName(hid_t loc_id, Grid *gridIn);
    void ReadStrPartitionInfo(hid_t loc_id, Grid *gridIn);

    void ReadUnstructGrid(hid_t loc_id, int iZone);
    void ReadFaceTopology(hid_t loc_id, Grid *gridIn);
    void ReadCellTopology(hid_t loc_id, Grid *gridIn);
    void ReadUnstrBCName(hid_t loc_id, Grid *gridIn);
    void ReadInterpointInfo(hid_t loc_id, Grid *gridIn);
    void ReadUnstrPartitionInfo(hid_t loc_id, Grid *gridIn);

    void SetNZones(int nZonesIn) { this->numberOfZones = nZonesIn; }
    int GetNZones() const { return numberOfZones; }

    void RedistributeZonesIfNeed();
    bool IsZoneLayoutInitialized();
    void SetGlobalZoneLayout();
    void ClassifyGridSystem();
    void AddGrid(Grid *grid);
    void CreateZones();
    void ReadWallDist();
    void SetMixGridInfo();

    void SetBlockProcDump(int *blockProcDumpIn) { this->blockProcDump = blockProcDumpIn; }
    void SetBlockProcGrid(int *blockProcGridIn) { this->blockProcGrid = blockProcGridIn; }
    void SetBlockProc(int *blockProcIn) { this->blockProc = blockProcIn; }
    void SetBlockType(int *blockTypeIn) { this->blockType = blockTypeIn; }
    void SetBlockIndex(int *blockIdxIn) { this->blockIdx = blockIdxIn; }
    int * GetBlockProcDump() { return this->blockProcDump; }
    int * GetBlockProc()     { return this->blockProc; }
    int * GetBlockProcGrid() { return this->blockProcGrid; }
    int * GetBlockType()     { return this->blockType; }
    int * GetBlockIndex()    { return this->blockIdx; }

private:
    //! The file index which the current processor charged.
    int fileID;
    bool isFileProcessor;
    bool isP2PMode;

    //! The processor index of the current one.
    int myProcID;

    //! The processor index of the serve.
    int serverProcID;

    //! Number of global grid files.
    int globalFileNumber;
    vector < vector<string> > localGridNameList;

    //! The total number of global zones.
    int numberOfZones;

    //! Number of zones in each .fts file.
    VVInt nZonesList;

    uint_t numberOfGridGroups;

    int zoneStart;
    int presentGridGroup;
    int *blockProc;
    int *blockType;
    int *blockIdx ;
    int *blockProcDump;
    int *blockProcGrid;
    int *blockFileProc;
    int taskType;

    vector<Grid *> grids;
    Region *region;
};

//! Different from class IO_HDF5Read, class IO_ReadHDF5Grid only provide function interface for grid information reading.
class IO_ReadHDF5Grid
{
private:
    //! Grid file name.
    string gridFileName;

    //! File identifier.
    hid_t fileID;

    //! Number of zones in this grid file
    int numberOfZones;

    //! Zone index in this grid file.
    int *zoneIndexs;

    //! Zone type in this grid file.
    int *zoneTypes;

    //! Default processor index of zone
    int *zoneProcs;

public:
    //! @param[in] gridFileNameIn    Grid file name.
    IO_ReadHDF5Grid(const string &gridFileNameIn);
    ~IO_ReadHDF5Grid();

public:
    //! Open grid file.
    void OpenFile();

    //! Close grid file.
    void CloseFile();

    //! Get number of zones in this grid file.
    //! @param[out] nZones    Number of zones.
    void GetNumberOfZones(int &nZones);

    //! Get zone index in this grid file.
    //! @param[out] zoneIndexOut    Index of each zone, size = nZones.
    void GetZoneIndex(int *zoneIndexOut);

    //! Get zone type in this grid file.
    //! @param[out] zoneTypeOut     Type of each zone, size = nZones.
    void GetZoneType(int *zoneTypeOut);

    //! Get default processor index of zone in this grid file.
    //! @param[out] zoneProcOut     Default processor of each zone, size = nZones.
    void GetZoneProc(int *zoneProcOut);

    //! Get single grid information.
    void GetGridInfo(Grid *gridOut);

    //! Get single grid size.
    //! @param[in]  zoneIndex          Index of zone.
    //! @param[out] iDimensionsIn      Zone size.
    //!             iDimensionsIn[0]   I(Structed grid) / Total node(Unstructed grid).
    //!             iDimensionsIn[1]   J(Structed grid) / Total face(Unstructed grid).
    //!             iDimensionsIn[2]   K(Structed grid) / Total cell(Unstructed grid).
    void GetZoneSize(int zoneIndex, int *iDimensionsOut);

    //! Get single grid coordinate.
    //! @param[in]  zoneIndex      Index of zone.
    //! @param[in]  iCoord         Coordinate index.
    //!             iCoord = 0,    X coordinate.
    //!             iCoord = 1,    Y coordinate.
    //!             iCoord = 2,    Z coordinate.
    //! @param[in]  dataSize       Total node(Unstructed grid) / I * J * K(Structed grid).
    //! @param[out] coordOut       Coordinate of zone, size = total node(Unstructed grid) / I * J * K(Structed grid).
    void GetZoneCoord(int zoneIndex, int iCoord, int dataSize, RDouble *coordOut);

    //! Get number of interface.
    //! @param[in]  zoneIndex      Index of zone.
    //! @param[out] nIFaceOut      Number of interface.
    void GetNIFace(int zoneIndex, int &nIFaceOut);

    //! Get interface information.
    //! @param[in]  zoneIndex        Index of zone.
    //! @param[out] interfaceInfo    Interface information.
    void GetInterfaceInfo(int zoneIndex, InterfaceInfo *interfaceInfo);

    //! Get number of boundary(Structed grid).
    //! @param[in]  zoneIndex      Index of zone.
    //! @param[out] nBCRegionOut   Number of boundary.
    void GetNBCRegion(int zoneIndex, int &nBCRegionOut);

    //! Get boundary information(Structed grid).
    //! @param[in]  zoneIndex           Index of zone.
    //! @param[out] compositeBCRegion   Boundary information.
    void GetBCRegion(int zoneIndex, StructBCSet *compositeBCRegion);
    void GetBCName(int zoneIndex, StructBCSet *compositeBCRegion);

    //! Get face topology(Unstructed grid)
    //! @param[in]  zoneIndex                 Index of zone.
    //! @param[in]  dataSize                  Number of total face.
    //! @param[out] nodeNumberOfEachFaceOut   Number of nodes in each face.
    void GetUnsNodeNumOfEachFace(int zoneIndex, int dataSize, int *nodeNumberOfEachFaceOut);

    //! Get face topology(Unstructed grid)
    //! @param[in]  zoneIndex         Index of zone.
    //! @param[in]  dataSize          Sum of array nodeNumberOfEachFace.
    //! @param[out] face2NodeOut      Index of nodes in each face.
    void GetUnsFace2Node(int zoneIndex, int dataSize, int *face2NodeOut);

    //! Get face topology(Unstructed grid)
    //! @param[in]  zoneIndex            Index of zone.
    //! @param[in]  dataSize             Number of total face.
    //! @param[out] leftCellOfFaceOut    Index of cell on the left of each face.
    //! @param[out] rightCellOfFaceOut   Index of cell on the right of each face.
    void GetUnsLeftCellOfFace(int zoneIndex, int dataSize, int *leftCellOfFaceOut);
    void GetUnsRightCellOfFace(int zoneIndex, int dataSize, int *rightCellOfFaceOut);

    //! Get cell topology(Unstructed grid)
    //! @param[in]  zoneIndex                   Index of zone.
    //! @param[in]  dataSize                    Number of total cell.
    //! @param[out] nodeNumberOfEachCellOut     Number of nodes in each cell.
    void GetUnsNodeNumOfEachCell(int zoneIndex, int dataSize, int *nodeNumberOfEachCellOut);

    //! Get cell topology(Unstructed grid)
    //! @param[in]  zoneIndex           Index of zone.
    //! @param[in]  dataSize            Sum of array nodeNumberOfEachCell.
    //! @param[out] cell2NodeOut        Index of nodes in each cell.
    void GetUnsCell2Node(int zoneIndex, int dataSize, int *cell2NodeOut);

    //! Get number of boundary face(Unstructed grid).
    //! @param[in]  zoneIndex           Index of zone.
    //! @param[out] nBoundFaceOut       Number of boundary face.
    void GetUnsNBoundFace(int zoneIndex, int &nBoundFaceOut);

    //! Get boundary information(Unstructed grid).
    //! @param[in]  zoneIndex           Index of zone.
    //! @param[in]  dataSize            Number of boundary face.
    //! @param[out] bcInfoOut           Boundary information.
    void GetUnsBCInfo(int zoneIndex, int dataSize, UnstructBCSet **bcInfoOut);
    void GetUnsBCName(int zoneIndex, int dataSize, UnstructBCSet **bcInfoOut);

    //! Get ordinary grid information(Structed grid)
    void ReadStrOrdinaryGridIndex(int zoneIndex, int &ordinaryGridIndexOut);
    void ReadStrOrdinaryDimStart(int zoneIndex, int *ordinaryDimStartOut);
    void ReadStrOrdinaryDimEnd(int zoneIndex, int *ordinaryDimEndOut);

    //! Get ordinary grid information(Unstructed grid)
    void ReadUnsOrdinaryGridIndex(int zoneIndex, int &ordinaryGridIndexOut);
    void ReadUnsOrdinaryNodeIndex(int zoneIndex, int *ordinaryNodeIndexOut, int nTotalNode);
    void ReadUnsOrdinaryFaceIndex(int zoneIndex, int *ordinaryFaceIndexOut, int nTotalFace);
    void ReadUnsOrdinaryCellIndex(int zoneIndex, int *ordinaryCellIndexOut, int nTotalCell);

private:
    void GetStructGridInfo(Grid *gridOut);
    void GetUnstructGridInfo(Grid *gridOut);

private:
    bool fileOpened;
};

LIB_EXPORT bool ReadHDF5File(Region *region_in, PHVectorString1D &gridGroupNameList);
LIB_EXPORT void DumpHDF5Grid(const string &gridFileName, Grid **grids, int nBlocks);

}