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
//! @file      Geo_UnstructGridPartition.h
//! @brief     Unstructured grid partition: serial and parallel.
//! @author    Bell.

#pragma once
#include "Geo_UnstructGrid.h"
#include "TK_Exit.h"
#ifdef USE_WINDOWS_X64
#include "parmetis64.h"
#else
#include "parmetis.h"
#endif
using namespace std;

const int SERIAL_PARTITION    = 0;
const int PARMETIS_METHOD     = 1;
const int METIS_METHOD        = 2;
const int METIS_OPENMP_METHOD = 3;

namespace PHSPACE
{
//! @brief Pre_UnstructGridPartition class defines the method of unstructured grid partition.
//! Serial grid partition using single processor.
//! Parallel grid partition using many processors.
class Pre_UnstructGridPartition
{
    static const int NOT_IN_THIS_ZONE = -2;
    static const int IN_THIS_ZONE     = -1;
private:
    //! The original source grid that want to be decomposed.
    UnstructGrid *gridSource;

    //! The partition grids, which is divided from original grid.
    Grid **grids;

    //! The number of new partition wanted.
    //! For serial partition, it is equal to the max processor.
    //! For parallel partition, it is equal to the GLOBAL max processor.
    int maxProcessor;

    //! Original grid file name, if the original grid does not converted from grid file,
    //! it could be NON, that is mean orgGridFileName = "".
    string orgGridFileName;

    //! Parallel partition method:
    //!   1 -- using ParMETIS for homogeneous MPI.
    //!   2 -- using METIS for homogeneous MPI.
    //!   3 -- using METIS partition for homogeneous OpenMP.
    int parallelPartitionMethod;

    //! OpenMP partition method:
    //!   0 -- NON used.
    //!   1 -- Equally partition for each zone.
    int ompPartMethod;
public:
    //! @param[in] gridSource          Original source grid need to be decomposed.
    //! @param[in] nPartitionBlocks    The number of GLOBAL blocks wanted. It is equal to the number of decomposed 
    //!                                zones for serial unstructured grid partition.
    //! @param[in] buildInterfaceMethod  Interfaces building method, used default value of '1' PLEASE!
    //! @param[in] parallelPartitionMethod Parallel partition method:
    //!  - # 1 using ParMETIS for homogeneous MPI.
    //!  - # 2 using METIS for homogeneous MPI.
    //!  - # 3 using METIS partition for homogeneous OpenMP.
    LIB_EXPORT Pre_UnstructGridPartition(UnstructGrid *gridSource, int nPartitionBlocks, int buildInterfaceMethod, 
                                         string orgGridFileName = "", int parallelPartitionMethod = SERIAL_PARTITION);
    LIB_EXPORT ~Pre_UnstructGridPartition();

public:
    //! Unstructured grid partition serially.
    LIB_EXPORT void PartitionGridSerial(const string &out_grid_file);

    //! Unstructured grid partition parallel.
    LIB_EXPORT void PartitionGridParallel(const string &out_grid_file);

    //! Dump out the partition grid to file of .fts format.
    LIB_EXPORT void DumpPartitionGrids(const string &out_grid_file);

    //! Return the partition grids.
    Grid ** GetPartitionGrids();

    //! Set parallel partition method:
    //!   1 -- using ParMETIS.
    //!   2 -- using METIS.
    void SetParallelPartitionMethod(int parallelPartitionMethod);

    //void SetOMPPartMethod(int ompPartMethod);

private:
    void Init();
    Grid ** CreateAllPartGrid();
    void CreatePartition();
    void ComputeGridForPartition(UnstructGrid *gridpart, int iZone);
    void GetVertexDist();
    bool IsParallelPartition();
    int GetMaxProcessorIndexInThisZone();
    idx_t * CreatePartitionSerial(UnstructGrid *grid, int maxProcessor);
    idx_t * CreatePartitionParallel(UnstructGrid *grid, int maxProcessor, idx_t *vtxdist);
    void MetisGrid(idx_t nTotalCell, idx_t *xadj, idx_t *adjncy, idx_t nparts, idx_t *part);
    void ParMetisGrid(idx_t *xadj, idx_t *adjncy, idx_t nparts, idx_t *part, idx_t *vtxdist);
    void ConvertPartNumberToGlobal(idx_t *part, int nTotalCell);
    void CommunicateNeiborCell(int *neighborCell);
    void SwapInterface(int *neighborCell, int iZone, int jZone);
    bool FindMatchPart(UnstructGrid *grid, int izone_t, int iface_t, int ie_t, int *interface_interfaceid_t);
    void ReconInterfaces(Grid **grids, int iZone);
    void BuildInterfaceLink(Grid **grids, int maxProcessor);

    //! Get neighborzone and ghost cell for periodic boundary, used for interface reconstruct.
    void SetPeriodInfo(Grid **grids, int iZone);
    void MatchPeriodicFace(Grid **grids, int iZone, int iBFace, RDouble *xlist, RDouble *ylist, RDouble *zlist, int &FaceIndex);
    void ResetFace2Face(Grid **grids);

    void FillnTCellofEachZone(int *&nTCellOfZones);

    void DumpParallelPartitionGrid(const string &gridFileName, Grid **grids, int nBlocks);
    void DumpSerialPartitionGrid(const string &gridFileName, Grid **grids);
    void DumpParallelPartitionGrid(const string &gridFileName, Grid **grids, int nBlocks, int *block_index_in);

    //! Compute the processor index that each cell belongs to.
    void ComputeCellMark();

    void SetProcessorList();
    bool DoseProcExist(int iproc);

    void FreeLocal2All();
    void CreateLocal2All(int iZone, UnstructGrid *gridpart);

    void ComputeNodeMark(UnstructGrid *gridpart);
    void ComputeFaceMark(int iZone, UnstructGrid *gridpart);

    int * GetNode2All(UnstructGrid *gridpart);
    int * GetFace2All(UnstructGrid *gridpart);
    int * GetCell2All(int iZone, UnstructGrid *gridpart);

    //! Compute the point to zone ID.
    void ComputePoint2ZoneID();

    //! Compute the total cell number of the node in computing the node value.
    void ComputeCellNumberForNodeValue();

    void SetCoorPart(UnstructGrid *gridpart);
    void SetLinkPart(int iZone, UnstructGrid *gridpart);
    void SetWallDist(UnstructGrid *gridpart);
    void SetPartitionInfo(UnstructGrid *gridpart);
    void SetVolumeCondition(UnstructGrid *gridpart);

    void SetNodeNumberOfEachFaceAndFace2Node(UnstructGrid *gridpart);
    void SetFace2CellAndBC(int iZone, UnstructGrid *gridpart);
    void SetInterfaceInfoPart(int iZone, UnstructGrid *gridpart);

    //! Set the interpoint information in the partition for iZone.
    //! @param[in] iZone       the number of zone.
    //! @param[in] gridpart    the grid after partition.
    void SetInterpointInfoPart(int iZone, UnstructGrid *gridpart);

    //! Modify the interpoint information in the partition for iZone.
    //! @param[in] iZone    the number of zone.
    //! @param[in] grids    the grid after partition.
    void ModifyInterpointInfoPart(Grid **grids);
    void SetInterfaceInfoPartForParallelPartition(int iZone, UnstructGrid *gridpart);
    void SetNTotalCell(int iZone, UnstructGrid *gridpart);
    int  GetNTotalCell(int iZone, idx_t *part);

    void CommunicateInterfaceTopo();
    void SwapInterfaceTopo(int iZone, int jZone);
    void SetInterfaceInfoPartInOrgBC(Grid **grids);
    void SetInterfaceInfoPartInOrgBC(int iZone, UnstructGrid *gridpart);
    void BuildInterfaceLinkForParallelPartition(Grid **grids);

    void CheckDoseNeedPartitionAdditionalData();
    //void PartitionCellToNode();
    void WriteAdditionalData(const string &out_grid_file);
    void WriteCellToNode(const string &out_grid_file);
    void WriteCellToNode(fstream &file, UnstructGrid *gridPart);

    //! Write the interpoint information into the file.
    //! @param[in] out_grid_file    the file name for output.
    void WriteInterpointInfo(const string &out_grid_file);

    //! Write the interpoint information into the file.
    //! @param[in] file        the file  for output.
    //! @param[in] gridPart    the grid after partition.
    void WriteInterpointInfo(fstream &file, UnstructGrid *gridPart);

    void WriteFaceBC(const string &out_grid_file);
    void WriteFaceBC(fstream &file, UnstructGrid *gridPart);
    void WriteWalldist(const string &out_grid_file);
    void WriteWalldist(ActionKey *actkey, UnstructGrid *gridPart);

    void ComputeGlobalToLocalCellToNodeMapping(UnstructGrid *gridpart, int iZone);

    int GetNumberofLocalZones();

    bool IsBelongtoThisZone(int faceMark);

private:
    //! The node/face/cell index mapping from new partition to old zone.
    int *node2all, *face2all, *cell2all;

    //! pointNumber2zoneID: the sequence of grid points that can be repeated, corresponding to zoneID.
    //! point2ZoneID: According to the sequence of grid points pointNumber2zoneID, the corresponding number of zoneID.
    //! pointPosition2zoneID: Position of point number in pointNumber2zoneID, pointNumber2zoneID is repeatable.
    int *pointNumber2zoneID, *point2ZoneID, *pointPosition2zoneID;

    //! cellNumberForNodeValue is the number of cells used to compute the value of each node.
    int *cellNumberForNodeValue;

    //! The part(processor) index that each cell belongs to.
    idx_t *part;

    //! The local cell/node/face index in the new partition.
    //! In other words, the index mapping from old zone to new partition.
    int *cell_mark, *node_mark, *face_mark;

    //! Interface re-construct method.
    //! -# 1   : Construct by local search in the original grid.
    //! -# else: Construct by global search among all of the decomposed grids.
    int npartmethod;

    //! Vertex distribution.
    //! For grid partition, each cell is a vertex in the partition graph.
    idx_t *vtxdist;

    //! The processor list that plan to be partition.
    //! For serial partition, the list up through 0 to maxProcessor (but not included).
    //! For parallel partition, the list include the processors in this zone.
    set<int> processorList;

    //! The following several pointers are temporary used.
    //! Local interface index of the faces that belong to the original BC interfaces.
    //! Local BC face index of the faces that belong to the original BC interfaces.
    int *interfaceIndexInBC_send;
    int *interfaceIndexInBC_recv, *partIndexInNewPartion_recv;
    int *part2GlobalZoneIndex, *part2LocalZoneIndex;
    int *nIFaceInOrgBC;

    //! Does need partition cell to node.
    bool partitionCellToNode;

    //! Does need partition face BC information.
    bool partitionFaceBC;

    //! Does need partition WallDist information.
    bool partitionWallDist;

    //! Number of parts in each processor, if metis is used to parallel partition.
    int nPartInEachProcIfMetis;
};

#include "Pre_UnstructGridPartition.hxx"

}