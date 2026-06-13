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
//! @file      IO_GridReader.h
//! @brief     It defines the grid reading method.
//! @author    Bell.

#pragma once
#include "GridGroup.h"
using namespace std;

namespace PHSPACE
{

//! @brief IO_GridReader class defines the method of grid reader.\n
class IO_GridReader
{
private:
    //! The grid group name list.
    vector< string > gridGroupNameList;

    //! Grid dimension.
    //! -# 2 : two dimensional grid.
    //! -# 3 : three dimensional grid.
    int dimension;

    //! Grid group to store the grids.
    GridGroup * gridGroup;

    //! Number of global grid files.
    int globalFileNumber;

    //! Grid groups for multi-files
    vector< GridGroup * > gridGroups;
public:
    //! Construct the GridReader.
    //! @param[in] gridFileName    the grid file name.
    //! @param[in] dimension       Grid dimension.
    //!                            -# 2 : two dimensional grid.
    //!                            -# 3 : three dimensional grid.
    //! @param[in] gridType    Grid type.
    //!            -# PHSPACE::UNSTRUCTGRID (0): unstructured grid.
    //!            -# PHSPACE::STRUCTGRID   (1): structured grid.
    //!            -# PHSPACE::MIXGRID      (2): mixed grid, this type is not support for other developer.
    LIB_EXPORT IO_GridReader(PHVectorString1D & gridGroupNameListIn, int dimension);

    LIB_EXPORT ~IO_GridReader();

public:
    //! Read ordinary grid.
    LIB_EXPORT void ReadGridsByP2PMode();

    //! Reading in Multi Files (need further revise)
    LIB_EXPORT void ReadGridsByCSMode();

    //! Read the cell to node topology by cs mode.
    LIB_EXPORT void ReadCellToNodeByCSMode();

    //! Read the interpoint to interpoint info by cs mode.
    LIB_EXPORT void ReadInterpointInfoByCSMode();

    //! Read the cell to node topology by p2p mode.
    LIB_EXPORT void ReadCellToNodeByP2PMode();

    //! Read the interpoint to interpoint info by P2P mode.
    LIB_EXPORT void ReadInterpointInfoByP2PMode();

    //! Read the walldist information by cs mode.
    LIB_EXPORT void ReadWallDistByCSMode();

    //! Read the walldist information by p2p mode.
    LIB_EXPORT void ReadWallDistByP2PMode();

    //! Read the Face bc topology by cs mode.
    LIB_EXPORT void ReadFaceBCByCSMode();

    //! Read the face direction with matching target face by cs mode.
    LIB_EXPORT void ReadFaceBcDirByCSMode();

    //! Read the Face bc topology by p2p mode.
    LIB_EXPORT void ReadFaceBCByP2PMode();

    //! Read the face direction with matching target face by p2p mode.
    LIB_EXPORT void ReadFaceBcDirByP2PMode();

    //! Return number of grid zones.
    LIB_EXPORT int GetNumberOfGrids();

    //! Return number of total Global simulation grid zones.
    LIB_EXPORT int GetNumberofGlobalZones();

    //! Return the i-th grid.
    LIB_EXPORT Grid * GetGrid(int iGrid);

    LIB_EXPORT GridGroup * GetGridGroup();

    LIB_EXPORT vector< GridGroup * > & GetGridGroups();

    //! Return the number of grid files.
    LIB_EXPORT std::size_t GetNumberOfGridFile() const;


    void InitGridFileNameListAndFileProcessor();

private:
    void ReadNumberOfTotalZonesAndInitialGlobalValue();
    int  CheckNumberofFiles();
    bool CheckGridConvertionMixGrid();
    void ReadGrid(const string & gridFileName, int startZoneIndex);
    void ReadGridGroup(const string & gridFileName, int startZoneIndex);
    void ReadNumberOfZones(fstream & file, int & nZones);
    void ReadOneZoneGrid(fstream & file, int & iZone);
    void ReadZoneLayout(fstream & file, int nZones);
    void ReadZoneLayoutCSMode(fstream & file, int nZones, int startZoneIndex);
    void RedistributeZonesIfNeed(int startZoneIndex);
    void ClassifyGridSystem();
    void DecompressDataToGrid(DataContainer *cdata, int iZone);
    void GetSendRecvProcID(int iZone, int &send_proc, int &recv_proc);
    void GetSendRecvProcID(int localZoneID, int &send_proc, int &recv_proc, int &recv_proc_grid);
    void Read_Bcast(fstream &file, void *data, int size);
    void Read(fstream &file, void *data, int size);
    void Bcast(void *data, int size, int proc, int tag);
    void ReadCompressedData(fstream &file, DataContainer *cdata, int send_proc, int recv_proc, int tag = 0);
    void Open(fstream &file, const string &filename, const ios_base::openmode &openmode);
    void OpenSerialFile(fstream &file, const string &filename, const ios_base::openmode &openmode);
    void Close(fstream &file);

    void DecompressDataToCellNode(DataContainer *cdata, int iZone);

    //! Decompress the interpoint information from the DataContainer.
    //! @param[in] cdata               the DataContainer which contains the interpoint information.
    //! @param[in] iZone               the zone number.\n
    void DecompressInterpointInfo(DataContainer *cdata, int iZone);
    void DecompressDataToFaceBC(DataContainer *cdata, int iZone);
    void DecompressDataToFaceBcDir(DataContainer *cdata, int iZone);
    void DecompressDataToWallDist(DataContainer *cdata, int iZone);

private:
    //! The file index which the current processor charged.
    int fileID;

    //! The processor index of the current one.
    int myProcID;

    //! The processor index of the serve.
    int serverProcID;

    //! Number of zones in each .fts file.
    int * nZonesList;
    vector < int > numberOfZones;

    //! If the cell to node file exist.
    bool cellnodeFileExist;

    //! If the cell to node file exist.
    bool interpointInfoFileExist;

    //! If the face BC file exist.
    bool faceBCFileExist;

    //! If the walldist file exist.
    bool walldistFileExist;

    int numberOfTotalZones;

    int presentGridGroup;
};

#include "IO_GridReader.hxx"





}