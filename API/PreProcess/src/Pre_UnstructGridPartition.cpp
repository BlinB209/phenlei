#include "Pre_UnstructGridPartition.h"
#include "IO_FileName.h"
#include "Glb_Dimension.h"
#ifdef USE_WINDOWS_X64
#else
#include "metis.h"
#endif
#include "Geo_UnstructBC.h"
#include "PHIO.h"
#include "Math_BasisFunction.h"
#include "Pre_GridConversion.h"
#include "Constants.h"
#include "IO_HDF5File.h"
#include "TK_Log.h"

using namespace std;

namespace PHSPACE
{
LIB_EXPORT Pre_UnstructGridPartition::Pre_UnstructGridPartition(UnstructGrid *gridSource, int nPartitionBlocks,
                                 int buildInterfaceMethod, string orgGridFileName, int parallelPartitionMethod)
{
    // if using ParMetis, the number of partitioned blocks is equal to the input.
    // else, it total nBlock = nProcessors * maxProcessor, where, maxProcessor in
    // this processor is the number of partitioned blocks.
    SetParallelPartitionMethod(parallelPartitionMethod);
    if (parallelPartitionMethod == METIS_METHOD || parallelPartitionMethod == METIS_OPENMP_METHOD)
    {
        this->maxProcessor = nPartitionBlocks;
        int nProcessors = PHMPI::GetNumberOfProcessor();
        if (nPartitionBlocks % nProcessors != 0)
        {
            ostringstream oss;
            oss << "Error: if using Metis to parallel partition, the number of zones wanted must INTER times "
                << "of the number of processors used !" << endl;
            PrintToWindow(oss);
        }
        this->nPartInEachProcIfMetis = nPartitionBlocks / nProcessors;
    }
    else
    {
        this->maxProcessor = nPartitionBlocks;
    }
    
    grids = 0;

    this->gridSource  = gridSource;
    this->npartmethod = buildInterfaceMethod;

    node2all = 0;
    face2all = 0;
    cell2all = 0;
    part = 0;
    cell_mark = 0;
    node_mark = 0;
    face_mark = 0;

    vtxdist = 0;
    interfaceIndexInBC_send    = 0;
    interfaceIndexInBC_recv    = 0;
    partIndexInNewPartion_recv = 0;
    part2GlobalZoneIndex       = 0;
    part2LocalZoneIndex        = 0;
    nIFaceInOrgBC              = new int[maxProcessor];

    this->orgGridFileName = orgGridFileName;

    parallelPartitionMethod = 1;

    partitionCellToNode = false;
    partitionFaceBC     = false;
    partitionWallDist   = false;
}

LIB_EXPORT Pre_UnstructGridPartition::~Pre_UnstructGridPartition()
{
    delete [] cell_mark;    cell_mark = NULL;
    delete [] face_mark;    face_mark = NULL;
    delete [] node_mark;    node_mark = NULL;
    delete [] part;    part = NULL;

    delete [] vtxdist;    vtxdist = NULL;

    delete [] interfaceIndexInBC_send;    interfaceIndexInBC_send = NULL;
    delete [] interfaceIndexInBC_recv;    interfaceIndexInBC_recv = NULL;
    delete [] partIndexInNewPartion_recv;    partIndexInNewPartion_recv = NULL;
    delete [] part2GlobalZoneIndex;    part2GlobalZoneIndex = NULL;
    delete [] part2LocalZoneIndex;    part2LocalZoneIndex = NULL;
    delete [] nIFaceInOrgBC;    nIFaceInOrgBC = NULL;

    delete [] pointNumber2zoneID;    pointNumber2zoneID = NULL;
    delete [] point2ZoneID;    point2ZoneID = NULL;
    delete [] pointPosition2zoneID;    pointPosition2zoneID = NULL;
    delete [] cellNumberForNodeValue;    cellNumberForNodeValue = NULL;

    FreeLocal2All();
}

void Pre_UnstructGridPartition::Init()
{
    CheckDoseNeedPartitionAdditionalData();
}

void Pre_UnstructGridPartition::CheckDoseNeedPartitionAdditionalData()
{
    if (this->orgGridFileName == "")
    {
        partitionCellToNode = false;
        partitionFaceBC     = false;
        return;
    }

    string originalCellToNodeFileName, originalFaceBCFileName, originalWallDistFileName;
    if (PHMPI::GetNumberOfProcessor() == 1)
    {
        originalCellToNodeFileName = ChangeExtensionOfFileName(orgGridFileName, "c2n");
        originalFaceBCFileName     = ChangeExtensionOfFileName(orgGridFileName, "bc");
        originalWallDistFileName   = ChangeExtensionOfFileName(orgGridFileName, "wdt");
    }
    else
    {
        if (parallelPartitionMethod == PARMETIS_METHOD)
        {
            partitionCellToNode = false;
            partitionFaceBC     = false;
            partitionWallDist   = false;
            return;
        }
        else
        {
            // Check if the first .c2n file exist or not.
            originalCellToNodeFileName = ChangeExtensionOfFileName(orgGridFileName, "c2n");
            originalFaceBCFileName     = ChangeExtensionOfFileName(orgGridFileName, "bc");
            originalWallDistFileName   = ChangeExtensionOfFileName(orgGridFileName, "wdt");
        }
    }
    
    //! Check if c2n file exist.
    partitionCellToNode = this->gridSource->ExistCell2Node();

    //! Check if face BC file exist.
    partitionFaceBC = SerialFileExist(originalFaceBCFileName);

    //! Check if walldist file exist.
    if (!this->gridSource->GetWallDist())
    {
        partitionWallDist = false;
    }
    else
    {
        partitionWallDist = true;
    }
}

Grid ** Pre_UnstructGridPartition::CreateAllPartGrid()
{
    this->CreatePartition();
    
    this->ComputeCellMark();
    
    this->ComputePoint2ZoneID();
    
    gridSource->computeBoundaryPointLabel();
    
    Grid **grids = new Grid * [maxProcessor];

    int nLocalZones = this->GetNumberofLocalZones();
    int interval = GetProgressInterval(nLocalZones, 10);
    int myid = PHMPI::GetCurrentProcessorID();
    int count = 0;
    
    ComputeCellNumberForNodeValue();

    for (int iproc = 0; iproc < maxProcessor; ++ iproc)
    {
        if (!DoseProcExist(iproc))
        {
            grids[iproc] = 0;
            continue;
        }

        ++ count;
        if (myid == PHMPI::GetServerProcessorID() && count > 0 && count % interval == 0)
        {
            ProgressMonitorToLogFile(count, nLocalZones, "grid partition reconstructing");
            ProgressMonitorToWindows(count, nLocalZones, "grid partition reconstructing");
        }

        GridID *index = new GridID(iproc);
        grids[iproc] = CreateGridGeneral(UNSTRUCTGRID, index, 0, GetDim());
        
        ComputeGridForPartition(static_cast<UnstructGrid *>(grids[iproc]), iproc);
        
        ComputeGlobalToLocalCellToNodeMapping(static_cast<UnstructGrid *>(grids[iproc]), iproc);
    }

    if (!IsParallelPartition())
    {
        ModifyInterpointInfoPart(grids);
    }
    else
    {
        PrintToWindow("Warning: Haven't consider inter-point for parallel partition!\n");
    }
    
    SetInterfaceInfoPartInOrgBC(grids);
    
    BuildInterfaceLinkForParallelPartition(grids);
    
    return grids;
}

void Pre_UnstructGridPartition::FreeLocal2All()
{
    delete [] node2all;
    delete [] face2all;
    delete [] cell2all;

    node2all = 0;
    face2all = 0;
    cell2all = 0;
}
#pragma warning(disable:4100)
LIB_EXPORT void Pre_UnstructGridPartition::PartitionGridSerial(const string &out_grid_file)
{
    Init();

    grids = this->CreateAllPartGrid();
    
    BuildInterfaceLink(grids, maxProcessor);
}

LIB_EXPORT void Pre_UnstructGridPartition::PartitionGridParallel(const string &out_grid_file)
{
    Init();

    GetVertexDist();

    grids = CreateAllPartGrid();

    // The following function is need only for serial partition.\n
    // This has been don in BuildInterfaceLinkForParallelPartition.
    //BuildInterfaceLink(grids, maxProcessor); 
}
#pragma warning(default:4100)
void Pre_UnstructGridPartition::GetVertexDist()
{
    //if (parallelPartitionMethod == METIS_METHOD) return;

#if IDXTYPEWIDTH == 32
    PrintToWindow("  Using 32-bit METIS/PARMETIS lib ...\n");
#else
    PrintToWindow("  Using 64-bit METIS/PARMETIS lib ...\n");
#endif

    int nZones = PHMPI::GetNumberofGlobalZones();
    int nProc  = PHMPI::GetNumberOfProcessor();
    ASSERT(nProc == nZones);
    int *nTCellofEachZone = new int[nProc];

    FillnTCellofEachZone(nTCellofEachZone);

    vtxdist = new idx_t[nZones + 1];
    vtxdist[0] = 0;

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        vtxdist[iZone + 1] = vtxdist[iZone] + nTCellofEachZone[iZone];
    }

    if (PHMPI::GetCurrentProcessorID() == PHMPI::GetServerProcessorID())
    {
        int minNTCell = 10000000, maxNTCell = 0;
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            minNTCell = MIN(minNTCell, nTCellofEachZone[iZone]);
            maxNTCell = MAX(maxNTCell, nTCellofEachZone[iZone]);
        }

        ostringstream oss;
        oss << "  Before partition, Min && Max number of cell : " << minNTCell << ", " << maxNTCell << endl;
        oss << "                    Min to Max nTCell ratio : " << (minNTCell * 1.0) / (maxNTCell * 1.0) << endl;
        PrintToWindow(oss);
        WriteLogFile(oss);
    }

    delete [] nTCellofEachZone;
}

void Pre_UnstructGridPartition::CreatePartition()
{
    int nTotalCell = gridSource->GetNTotalCell();
    int nTotalFace = gridSource->GetNTotalFace();
    int nTotalNode = gridSource->GetNTotalNode();

    cell_mark = new int[nTotalCell];    //all 2 local
    face_mark = new int[nTotalFace];    //all 2 local
    node_mark = new int[nTotalNode];    //all 2 local

    if (!IsParallelPartition())
    {
        part = CreatePartitionSerial(gridSource, maxProcessor);
    }
    else
    {
        part = CreatePartitionParallel(gridSource, maxProcessor, vtxdist);
    }

    // The following section is used in parallel partition.
    // This pointers are used in communication.
    int nBoundFace = gridSource->GetNBoundFace();
    interfaceIndexInBC_send    = new int[nBoundFace];
    interfaceIndexInBC_recv    = new int[nBoundFace];
    partIndexInNewPartion_recv = new int[nBoundFace];
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        interfaceIndexInBC_send[iFace]    = -1;
        interfaceIndexInBC_recv[iFace]    = -1;
        partIndexInNewPartion_recv[iFace] = -1;
    }

    SetProcessorList();

}

idx_t * Pre_UnstructGridPartition::CreatePartitionSerial(UnstructGrid *grid, int maxProcessor)
{
    if (maxProcessor < 2)
    {
        TK_Exit::ExceptionExit("The number of partitions should be greater than 1!\n");
    }

    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();

    //! The data of large scale grid may exceed the value range of int type, so we change the
    //! data type from int to long-long-int. 
    long long int nTotalFaceLongType = nTotalFace;
    long long int nBoundFaceLongType = nBoundFace;
    long long int size = 2 * (nTotalFaceLongType - nBoundFaceLongType);

    idx_t *xadj   = new idx_t[nTotalCell + 1];
    idx_t *adjncy = new idx_t[size];
    idx_t *part   = new idx_t[nTotalCell];

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        part[iCell] = -1;
    }

    Get_Xadj_Adjncy(grid, xadj, adjncy);
    MetisGrid(nTotalCell, xadj, adjncy, maxProcessor, part);

    delete [] xadj;
    delete [] adjncy;

    return part;
}

idx_t * Pre_UnstructGridPartition::CreatePartitionParallel(UnstructGrid *grid, int maxProcessor, idx_t *vtxdist)
{
    if (maxProcessor < 2)
    {
        TK_Exit::ExceptionExit("The number of partitions should be greater than 1!\n");
    }

    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nIFace = grid->GetNIFace();

    idx_t *xadj   = new idx_t[nTotalCell + 1];
    idx_t *adjncy = new idx_t[2 * (nTotalFace - nBoundFace + nIFace)];
    idx_t *part   = new idx_t[nTotalCell];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        part[iCell] = -1;
    }

    int *neighborCell = new int[nIFace];
    CommunicateNeiborCell(neighborCell);

    if (parallelPartitionMethod == PARMETIS_METHOD)
    {
        Get_Xadj_Adjncy(grid, xadj, adjncy, neighborCell, vtxdist);
    }
    else
    {
        Get_Xadj_Adjncy(grid, xadj, adjncy, 0, 0);
    }
    
    if (parallelPartitionMethod == PARMETIS_METHOD)
    {
        PrintToWindow("Parallel grid partition by ParMETIS ...\n");

        ParMetisGrid(xadj, adjncy, maxProcessor, part, vtxdist);
    }
    else if (parallelPartitionMethod == METIS_METHOD || parallelPartitionMethod == METIS_OPENMP_METHOD)
    {
        PrintToWindow("Parallel grid partition by METIS ...\n");

        MetisGrid(nTotalCell, xadj, adjncy, nPartInEachProcIfMetis, part);
        ConvertPartNumberToGlobal(part, nTotalCell);
    }

    delete [] xadj;
    delete [] adjncy;
    delete [] neighborCell;

    return part;
}

void Pre_UnstructGridPartition::MetisGrid(idx_t nTotalCell, idx_t *xadj, idx_t *adjncy, idx_t nparts, idx_t *part)
{
    WriteLogFile("Now beginning partition graph by METIS!");

    idx_t ncon     = 1;
    idx_t *vwgt    = 0;
    idx_t *vsize   = 0;
    idx_t *adjwgt  = 0;
    real_t *tpwgts = 0;
    real_t *ubvec  = 0;
    idx_t options[METIS_NOPTIONS];
    idx_t objval = 0;

    METIS_SetDefaultOptions(options);
  
    PrintToWindow("Now begining partition graph!\n");
  
    if (nparts > 8)
    {
        PrintToWindow("Using K-way Partitioning!\n");
        //METIS4.0
        //METIS_PartGraphKway(&nTotalCell,xadj,adjncy,vwgt,adjwgt,&wgtflag, &numflag, 
        //                         &nparts,options,&edgecut,part);

        //METIS5.0
        METIS_PartGraphKway(&nTotalCell, &ncon, xadj, adjncy, vwgt, vsize, adjwgt, 
                           &nparts, tpwgts, ubvec, options, &objval, part);

        //METIS PartGraphVKway (int *n, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, int *wgtflag,
        //int *numflag, int *nparts, int *options, int *volume, idxtype *part);

        //int METIS_PartGraphKway(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy, 
        //                        idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *nparts, 
        //                        real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *objval, 
        //                        idx_t *part)
    }
    else
    {
        PrintToWindow("Using Recursive Partitioning!\n");
        //METIS4.0
        //METIS_PartGraphRecursive(&nTotalCell,xadj,adjncy,vwgt,adjwgt,&wgtflag, &numflag, 
        //                         &nparts,options,&edgecut,part);

        METIS_PartGraphRecursive(&nTotalCell, &ncon, xadj, adjncy, vwgt, vsize, adjwgt, 
                                 &nparts, tpwgts, ubvec, options, &objval, part);

        //METIS5.0
        //int METIS_PartGraphRecursive(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, 
        //                             idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, 
        //                             idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, 
        //                             idx_t *objval, idx_t *part)
    }

    ostringstream oss;
    oss << "The interface number: " << objval << endl; 
    oss << "Partition is finished!\n";
    PrintToWindow(oss);
}

void Pre_UnstructGridPartition::ParMetisGrid(idx_t *xadj, idx_t *adjncy, idx_t nparts, idx_t *part, idx_t *vtxdist)
{
    WriteLogFile("Now beginning parallel partition graph by ParMETIS!");
    PrintToWindow("Now beginning parallel partition graph by ParMETIS!\n");

    idx_t ncon     = 1;
    idx_t *vwgt    = 0;
    idx_t *adjwgt  = 0;
    real_t *tpwgts = new real_t[nparts];
    for (int iPart = 0; iPart < nparts; ++ iPart)
    {
        tpwgts[iPart] = static_cast<real_t> (1.0 / nparts);
    }

    //! parmetisBalance: is used to specify the imbalance tolerance.
    //! 1     : perfect balance;
    //! nparts: perfect imbalance;
    //! 1.05  :recommended.
    RDouble parmetisBalance = GlobalDataBase::GetDoubleParaFromDB("parmetisBalance");

    real_t *ubvec = new real_t[ncon];
    for (int icon = 0; icon < ncon; ++ icon)
    {
        ubvec[icon] = static_cast<real_t> (parmetisBalance);    //modified by LIJIN,2019-2-19
    }
    idx_t options[METIS_NOPTIONS];
    options[0] = 0;
    idx_t wgtflag = 0;
    idx_t numflag = 0;
    idx_t edgecut = 0;

    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    
    if (PHMPI::GetNumberOfProcessor() == 1)
    {
        TK_Exit::ExceptionExit("Error: number of processors must great than one for parallel partition!\n");
    }
    else
    {
        ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, &ncon,
            &nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
    }

    delete [] ubvec;
    delete [] tpwgts;

    WriteLogFile("Partition is finished by Parallel Metis!");
    PrintToWindow("Partition is finished by Parallel Metis!\n");
}

void Pre_UnstructGridPartition::ConvertPartNumberToGlobal(idx_t *part, int nTotalCell)
{
    // Each Processor owns one zone, and each one is partition to be same number of zones.
    // so, the start zone index in this processor(zone) is start from:
    // number of pre-zones * nPartInEachProcIfMetis in each original zones.
    int myid = PHMPI::GetCurrentProcessorID();
    int startZoneID = myid * this->nPartInEachProcIfMetis; 

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        part[iCell] += startZoneID;
    }
}

void Pre_UnstructGridPartition::CommunicateNeiborCell(int *neighborCell)
{
    WriteLogFile("Start to Communicate Neighbor Cells");

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    // Each processor runs along the same order, which is from 0 to nZones.
    // Each simulation zone is looped over to judge if the communication is necessary.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();
        for (std::size_t ineighbor = 0; ineighbor < numberOfNeighbor; ++ ineighbor)
        {
            int neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor(ineighbor);
            SwapInterface(neighborCell, iZone, neighborZoneIndex);
        }
    }

    WriteLogFile("End to Communicate Neighbor Cells");
}

void Pre_UnstructGridPartition::SwapInterface(int *neighborCell, int iZone, int jZone)
{
    using namespace PHMPI;

    //! iZone is the current zone, and jZone is the target zone.
    int send_proc = GetZoneProcessorID(iZone);
    int recv_proc = GetZoneProcessorID(jZone);

    int tag = iZone;
    int myid = GetCurrentProcessorID();

    DataContainer *cdata = new DataContainer();
    cdata->MoveToBegin();

    if (myid == send_proc)
    {
        UnstructGrid  *grid  = UnstructGridCast(GetGrid(iZone, 0));
        InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
        if (interfaceInfo)
        {
            int ineighbor = interfaceInfo->FindIthNeighbor(jZone);
            int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
            int *faceIndexForSend = interfaceInfo->GetFaceIndexForSend(ineighbor);
            int *qtmp = new int[nIFaceOfNeighbor];
            int iFace, sourceCell;
            for (int iLocalFace = 0; iLocalFace < nIFaceOfNeighbor; ++ iLocalFace)
            {
                iFace = faceIndexForSend[iLocalFace];
                grid->GetSourceIndex(iFace, 1, sourceCell);
                qtmp[iLocalFace] = sourceCell;
            }
            cdata->Write(qtmp, nIFaceOfNeighbor * sizeof(int));
            delete [] qtmp;    qtmp = nullptr;
        }
    }

    PH_Trade(cdata, send_proc, recv_proc, tag); 

    if (myid == recv_proc)
    {
        UnstructGrid *grid = UnstructGridCast(GetGrid(jZone, 0));
        InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
        if (interfaceInfo)
        {
            int ineighbor = interfaceInfo->FindIthNeighbor(iZone);
            int iFace;
            int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
            int *faceIndexForRecv = interfaceInfo->GetFaceIndexForRecv(ineighbor);
            cdata->MoveToBegin();
            for (int iLocalFace = 0; iLocalFace < nIFaceOfNeighbor; ++ iLocalFace)
            {
                iFace = faceIndexForRecv[iLocalFace];
                int cellIndex;
                cdata->Read(&cellIndex, sizeof(int));
                neighborCell[iFace] = cellIndex;
            }
        }
    }

    delete cdata;    cdata = nullptr;
}

void Pre_UnstructGridPartition::SetProcessorList()
{
    if (!IsParallelPartition()) return;

    int nTotalCell = this->gridSource->GetNTotalCell();
    
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int partIndex = static_cast<int>(part[iCell]);

        processorList.insert(partIndex);
    }

    ostringstream oss;
    oss << "  Part index in processor " << PHMPI::GetCurrentProcessorID() << ": ";
    for (set<int>::iterator iter = processorList.begin(); iter != processorList.end(); ++ iter)
    {
        int partID = *iter;
        if (partID < 0)
        {
            TK_Exit::ExceptionExit("Error: partID < 0!");
        }
        oss << partID << " ";
    }
    oss << endl;
    WriteLogFile(oss.str());
    
    //! Compute the global zone index of each cell.
    int nLocalZones = this->GetNumberofLocalZones();

    int nZones = PHMPI::GetNumberofGlobalZones();
    int nProc  = PHMPI::GetNumberOfProcessor();
    ASSERT(nProc == nZones);
    int myid = PHMPI::GetCurrentProcessorID();
    int *nPartOfEachProcessor = new int[nZones];

    if (nProc == nZones)
    {
        //! Each processors has ONE zone.
        int *localSizeTemp = new int [nProc];
        SetField(localSizeTemp, 1, nProc);
        int **nPartListTemp = new int * [nProc];

        PHMPI::FillGobalDataCollectBcast(&nLocalZones, localSizeTemp, nPartListTemp);

        for (int iProc = 0; iProc < nProc; ++ iProc)
        {
            nPartOfEachProcessor[iProc] = nPartListTemp[iProc][0];
        }

        DelPointer2(nPartListTemp);
        DelPointer(localSizeTemp);
    }
    else
    {
        TK_Exit::ExceptionExit("number of processors != number of zones.");
        //PHMPI::FillGobalData(nLocalZones, nPartOfEachProcessor);
    }

    part2LocalZoneIndex  = new int[maxProcessor];
    part2GlobalZoneIndex = new int [maxProcessor];
    for (int iproc = 0; iproc < maxProcessor; ++ iproc)
    {
        part2LocalZoneIndex[iproc]  = -1;
        part2GlobalZoneIndex[iproc] = -1;
    }

    int localIndex = 0;
    for (int iproc = 0; iproc < maxProcessor; ++ iproc)
    {
        if (processorList.find(iproc) != processorList.end())
        {
            part2LocalZoneIndex[iproc] = localIndex ++;
        }
    }
    
    int count = 0;
    for (int iProc = 0; iProc < myid; ++ iProc)
    {
        count += nPartOfEachProcessor[iProc];
    }
    for (set<int>::iterator iter = processorList.begin(); iter != processorList.end(); ++ iter)
    {
        int partIndex = *iter;
        localIndex = part2LocalZoneIndex[partIndex];
        part2GlobalZoneIndex[partIndex] = localIndex + count;
    }

    delete [] nPartOfEachProcessor;
}

int Pre_UnstructGridPartition::GetNumberofLocalZones()
{
    int nLocalZones = maxProcessor;
    if (parallelPartitionMethod == PARMETIS_METHOD)
    {
        nLocalZones = static_cast<int>(processorList.size());
    }
    else if (parallelPartitionMethod == METIS_METHOD || parallelPartitionMethod == METIS_OPENMP_METHOD)
    {
        nLocalZones = static_cast<int>(processorList.size());
    }

    return nLocalZones;
}

bool Pre_UnstructGridPartition::DoseProcExist(int iproc)
{
    if (!IsParallelPartition()) return true;

    if (processorList.find(iproc) == processorList.end())
    {
        return false;
    }
    else
    {
        return true;
    }
}

int Pre_UnstructGridPartition::GetMaxProcessorIndexInThisZone()
{
    int nTotalCell = this->gridSource->GetNTotalCell();
    int maxNumber = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int partIndex = static_cast<int>(part[iCell]);
        maxNumber = MAX(maxNumber, partIndex);
    }
    return (maxNumber + 1);
}

void Pre_UnstructGridPartition::ComputeCellMark()
{
    int nTotalCell = gridSource->GetNTotalCell();

    int  maxProcessor = GetMaxProcessorIndexInThisZone();

    int *count_array = new int[maxProcessor];

    for (int i = 0; i < maxProcessor; ++ i)
    {
        count_array[i] = 0;
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        cell_mark[iCell] = count_array[ part[iCell] ] ++;
    }

    delete [] count_array;
}


//! Set the interface information of the interfaces that belong to the original BC faces.\n
//! The interfaces that belong to original interior faces have been set by SetInterfaceInfoPart().
void Pre_UnstructGridPartition::SetInterfaceInfoPartInOrgBC(Grid **grids)
{
    if (!IsParallelPartition()) return;

    if (npartmethod != 1) return;
    
    CommunicateInterfaceTopo();
    
    for (int iproc = 0; iproc < maxProcessor; ++ iproc)
    {
        if (!DoseProcExist(iproc))
        {
            grids[iproc] = 0;
            continue;
        }

        SetInterfaceInfoPartInOrgBC(iproc, UnstructGridCast(grids[iproc]));
    }
}

void Pre_UnstructGridPartition::SetInterfaceInfoPartInOrgBC(int iZone, UnstructGrid *gridpart)
{
    if (npartmethod != 1) return;
    
    int nBoundFace = gridSource->GetNBoundFace();
    int nTotalCell = gridSource->GetNTotalCell();
    int *left_cell_of_face  = gridSource->GetLeftCellOfFace();
    int *right_cell_of_face = gridSource->GetRightCellOfFace();
    int le, re;

    UnstructBCSet **bcr = gridSource->GetBCRecord();

    //! Re-compute bc_face_mark for each sub-partition.
    //! The original face_mark have been overwritten.
    int *bc_face_mark = new int [nBoundFace];
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        bc_face_mark[iFace] = -2;
    }
    
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = left_cell_of_face [iFace];
        re = right_cell_of_face[iFace];
        if (part[le] == iZone)
        {
            bc_face_mark[iFace] = -1;
        }
        else if (re < nTotalCell && part[re] == iZone)
        {
            bc_face_mark[iFace] = -1;
        }
    }
    
    //! Only deal with interfaces that belong to original BC faces.
    InterfaceInfo *local_iinfo = gridpart->GetInterfaceInfo();
    int *local_interFace2ZoneID       = local_iinfo->GetInterFace2ZoneID();
    int *local_interFace2InterFaceID  = local_iinfo->GetInterFace2InterFaceID();
    int *local_interFace2BoundaryFace = local_iinfo->GetInterFace2BoundaryFace();

    int i_nTFace = 0;
    int i_nIFace = 0;
    
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        if (bc_face_mark[iFace] == -1)
        {
            bc_face_mark[iFace] = i_nTFace ++;

            int bcType = bcr[iFace]->GetKey();
            if (bcType == PHENGLEI::INTERFACE)
            {
                local_interFace2ZoneID[i_nIFace]       = partIndexInNewPartion_recv[iFace];
                local_interFace2InterFaceID[i_nIFace]  = interfaceIndexInBC_recv[iFace];
                local_interFace2BoundaryFace[i_nIFace] = bc_face_mark[iFace];

                i_nIFace ++;
            }
        }
    }

    delete [] bc_face_mark;
}

void Pre_UnstructGridPartition::BuildInterfaceLinkForParallelPartition(Grid **grids)
{
    if (!IsParallelPartition()) return;

    for (int iZone = 0; iZone < maxProcessor; ++ iZone)
    {
        if (!grids[iZone]) continue;

        InterfaceInfo *interfaceInfo = grids[iZone]->GetInterfaceInfo();
        if (!interfaceInfo) continue;

        int *interFace2ZoneID      = interfaceInfo->GetInterFace2ZoneID();
        int *interFace2InterFaceID = interfaceInfo->GetInterFace2InterFaceID();

        UnstructGrid *grid = UnstructGridCast(grids[iZone]);

        int nBoundFace = grid->GetNBoundFace();
        int nIFace = interfaceInfo->GetNIFace();
        int nphysi = nBoundFace - nIFace;    //! number of physical boundary = total number boundary face - number of interface.

        int *left_cell_of_face  = grid->GetLeftCellOfFace();
        int *right_cell_of_face = grid->GetRightCellOfFace();

        bool *isNeighborInterfaceFound = interfaceInfo->GetIsNeighborInterfaceFound();

        int nbz, le, re, ie;

        //! Interfaces count from 0.
        for (int iFace = 0; iFace < nIFace; ++ iFace)
        {
            if (iFace < nIFaceInOrgBC[iZone])
            {
                //! The interfaces in the original BC \n
                //! have been built in function SetInterfaceInfoPartInOrgBC.
                continue;
            }

            if (isNeighborInterfaceFound[iFace] == true)
            {
                //! The neighbor of this face has been found and paired, so skip it.
                continue;
            }
            
            //! Neighbor zone index.
            nbz = interFace2ZoneID[iFace];
            le  = left_cell_of_face [iFace + nphysi];
            re  = right_cell_of_face[iFace + nphysi];
            ie  = MAX(le, re);
            if (nbz >= iZone)
            {
                UnstructGrid *grid_t = UnstructGridCast(grids[nbz]);
                if (!FindMatchPart(grid_t, iZone, iFace, ie, interFace2InterFaceID))
                {
                    TK_Exit::ExceptionExit("!FindMatchPart\n");
                }
                else
                {
                    isNeighborInterfaceFound[iFace] = true;
                }
            }
        }
    }

    //! Map the local part index to global zone index.
    for (int iZone = 0; iZone < maxProcessor; ++ iZone)
    {
        if (!grids[iZone]) continue;

        InterfaceInfo *interfaceInfo = grids[iZone]->GetInterfaceInfo();
        if (!interfaceInfo) continue;

        int *interFace2ZoneID = interfaceInfo->GetInterFace2ZoneID();
        int nIFace = interfaceInfo->GetNIFace();
        for (int iFace = 0; iFace < nIFace; ++ iFace)
        {
            if (iFace < nIFaceInOrgBC[iZone])
            {
                //! The interfaces in the original BC \n
                //! have been built in function SetInterfaceInfoPartInOrgBC.
                continue;
            }
            
            int localPartIndex = interFace2ZoneID[iFace];
            interFace2ZoneID[iFace] = part2GlobalZoneIndex[localPartIndex];
        }
    }
}

void Pre_UnstructGridPartition::CommunicateInterfaceTopo()
{
    if (!IsParallelPartition()) return;

    int nZones = PHMPI::GetNumberofGlobalZones();

    //! Each processor runs along the same order, which is from 0 to nZones.
    //! Each simulation zone is looped over to judge if the communication is necessary.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();

        for (std::size_t ineighbor = 0; ineighbor < numberOfNeighbor; ++ ineighbor)
        {
            int neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor(ineighbor);

            SwapInterfaceTopo(iZone, neighborZoneIndex);
        }
    }
}

void Pre_UnstructGridPartition::SwapInterfaceTopo(int iZone, int jZone)
{
    using namespace PHMPI;

    //! iZone is the current zone, and jZone is the target zone.
    int send_proc = GetZoneProcessorID(iZone);
    int recv_proc = GetZoneProcessorID(jZone);

    int tag = iZone;
    int myid = GetCurrentProcessorID();

    DataContainer *cdata = new DataContainer();

    if (myid == send_proc)
    {
        UnstructGrid  *grid = UnstructGridCast(GetGrid(iZone, 0));
        InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
        if (interfaceInfo)
        {
            int ineighbor = interfaceInfo->FindIthNeighbor(jZone);
            int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
            int *faceIndexForSend = interfaceInfo->GetFaceIndexForSend(ineighbor);
            int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();
            int iFace, le, partIndex, localInterfaceIndex, bcFaceIndex;
            cdata->MoveToBegin();
            int nm = 2;
            int nlen = nm * nIFaceOfNeighbor;
            int ntmp = 0;
            int *qtmp = new int [nlen];
            for (int iLocalFace = 0; iLocalFace < nIFaceOfNeighbor; ++ iLocalFace)
            {
                iFace = faceIndexForSend[iLocalFace];
                bcFaceIndex = interFace2BoundaryFace[iFace];    // BC face index in old grid.
                grid->GetSourceIndex(iFace, 1, le);

                partIndex = part2GlobalZoneIndex[part[le]];    // global zone index in new partition.
                localInterfaceIndex = interfaceIndexInBC_send[bcFaceIndex]; // local interface index in new partition.

                qtmp[ntmp++] = partIndex;              // interFace2ZoneID
                qtmp[ntmp++] = localInterfaceIndex;    // interFace2InterFaceID
            }
            cdata->Write(qtmp, nlen * sizeof(int));
            delete [] qtmp;    qtmp = nullptr;
        }
    }

    PH_Trade(cdata, send_proc, recv_proc, tag);   

    if (myid == recv_proc)
    {
        UnstructGrid *grid = UnstructGridCast(GetGrid(jZone, 0));
        InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
        if (interfaceInfo)
        {
            int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

            int ineighbor = interfaceInfo->FindIthNeighbor(iZone);
            int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
            int *faceIndexForRecv = interfaceInfo->GetFaceIndexForRecv(ineighbor);
            int iFace, globalZoneIndexInNeighbor, localInterfaceIndexInNeighbor, bcFaceIndex;
            cdata->MoveToBegin();
            for (int iLocalFace = 0; iLocalFace < nIFaceOfNeighbor; ++ iLocalFace)
            {
                iFace = faceIndexForRecv[iLocalFace];
                bcFaceIndex = interFace2BoundaryFace[iFace];
                
                cdata->Read(&globalZoneIndexInNeighbor, sizeof(int));
                cdata->Read(&localInterfaceIndexInNeighbor, sizeof(int));

                partIndexInNewPartion_recv[bcFaceIndex] = globalZoneIndexInNeighbor;
                interfaceIndexInBC_recv[bcFaceIndex] = localInterfaceIndexInNeighbor;
            }
        }
    }

    delete cdata;    cdata = nullptr;
}

bool Pre_UnstructGridPartition::FindMatchPart(UnstructGrid *grid, int izone_t, int iface_t, int ie_t, int *interface_interfaceid_t)
{
    bool found = false;

    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return found;

    int nBoundFace = grid->GetNBoundFace();
    int nIFace = interfaceInfo->GetNIFace();
    int nphysi = nBoundFace - nIFace;    //! number of physical boundary = total number boundary face - number of interface.

    int *left_cell_of_face  = grid->GetLeftCellOfFace();
    int *right_cell_of_face = grid->GetRightCellOfFace();

    int *interFace2ZoneID      = interfaceInfo->GetInterFace2ZoneID();
    int *interFace2InterFaceID = interfaceInfo->GetInterFace2InterFaceID();
    
    bool *isNeighborInterfaceFound = interfaceInfo->GetIsNeighborInterfaceFound();

    int le, re, ie;

    //! Interfaces count from 0.
    for (int iFace = 0; iFace < nIFace; ++ iFace)
    {
        le = left_cell_of_face [iFace + nphysi];
        re = right_cell_of_face[iFace + nphysi];
        ie = MAX(le, re);

        if (isNeighborInterfaceFound[iFace] == true)
        {
            //! The neighbor of this face has been found and paired, so skip it.
            continue;
        }

        //! The interface iFace target zone == iZone.
        if ((interFace2ZoneID[iFace] == izone_t) &&
            (interFace2InterFaceID[iFace] == ie_t) && 
            (ie == interface_interfaceid_t[iface_t]))
        {
            //! The jFace of current zone and the iFace corresponding to target zone are all count from 0.
            interFace2InterFaceID  [iFace  ] = iface_t;
            interface_interfaceid_t[iface_t] = iFace;

            isNeighborInterfaceFound[iFace] = true;

            found = true;
            break;
        } 
    }

    return found;
}

//! For each interface in the current zone,\n
//! Find out the interface-index in the neighbor zone.
void Pre_UnstructGridPartition::BuildInterfaceLink(Grid **grids, int maxProcessor)
{
    if (npartmethod == 1)
    {
        for (int iProc = 0; iProc < maxProcessor; ++ iProc)
        {
            //! To build interface link for periodic boundary.
            SetPeriodInfo(grids, iProc);
            ReconInterfaces(grids, iProc);
        }
        ResetFace2Face(grids);
    }
    else
    {
        ConstructGlobalInterfaceLink(grids, maxProcessor);
    }
}

void Pre_UnstructGridPartition::SetPeriodInfo(Grid **grids, int iZone)
{   
    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");
    if (periodicType == NO_PERIODICITY)
    {
        return;
    }
    UnstructGrid *gridpart = UnstructGridCast(grids[iZone]);

    UnstructBCSet *unstructBCSet = gridSource->GetUnstructBCSet();
    InterfaceInfo *GlobalInterfaceInfo = gridSource->GetInterfaceInfo();
    if (!GlobalInterfaceInfo) return;

    int *GlobalinterFace2BoundaryFace = GlobalInterfaceInfo->GetInterFace2BoundaryFace();
    int *GlobalinterFace2ZoneID       = GlobalInterfaceInfo->GetInterFace2ZoneID();
    int *bcRegionIDofBCFace           = unstructBCSet->GetBCFaceInfo();
    int *rightCellOfFace              = gridSource->GetRightCellOfFace();
    int nBoundFace                    = gridSource->GetNBoundFace();

    RDouble *xfc = gridSource->GetFaceCenterX();
    RDouble *yfc = gridSource->GetFaceCenterY();
    RDouble *zfc = gridSource->GetFaceCenterZ();

    RDouble *localxfc = gridpart->GetFaceCenterX();
    RDouble *localyfc = gridpart->GetFaceCenterY();
    RDouble *localzfc = gridpart->GetFaceCenterZ();

    int *node_number_of_each_face = gridSource->GetNodeNumberOfEachFace();
    int *face2node = gridSource->GetFace2Node();

    UnstructBCSet *partBCSet = gridpart->GetUnstructBCSet();
    int i_nBFace = gridpart->GetNBoundFace();

    InterfaceInfo *interfaceInfo = gridpart->GetInterfaceInfo();
    int i_nIFace = interfaceInfo->GetNIFace();
    int *interFace2ZoneID = interfaceInfo->GetInterFace2ZoneID();
    int *interFace2CellID = interfaceInfo->GetInterFace2CellID();
    int *interFace2InterFaceID = interfaceInfo->GetInterFace2InterFaceID();
    int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();
    int *interFaceDirection = interfaceInfo->GetInterFaceDirection();

    bool *isNeighborInterfaceFound = interfaceInfo->GetIsNeighborInterfaceFound();

    RDouble translationLength[3] = { 0.0 };
    GlobalDataBase::GetData("translationLength", &translationLength, PHDOUBLE, 3);
    double rotationAngle = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("rotationAngle");
    rotationAngle = rotationAngle * PI / 180.0;

    int npartBound = gridpart->GetNBoundFace();
    int nIFace = interfaceInfo->GetNIFace();
    int nphysi = npartBound - nIFace;

    int count = 0;
    int Globalpos = 0;
    int maxlist = 0;

    //! The size of ii_local here is equal to periodic faces in current zone.
    int ii_local = 0;

    //! Find periodic boundary information in source grid.
    for (int iBFace = 0; iBFace < nBoundFace; iBFace++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iBFace]);
        string bcName = bcRegion->GetBCName();
        int bcType = bcRegion->GetBCType();

        if (bcName != "Periodic_up" && bcName != "Periodic_down")
        {
            Globalpos += node_number_of_each_face[iBFace];
            continue;
        }

        maxlist = MAX(maxlist, node_number_of_each_face[iBFace]);
        vector<int> index(maxlist), pindex(maxlist);

        RDouble *xlist = new RDouble[maxlist];
        RDouble *ylist = new RDouble[maxlist];
        RDouble *zlist = new RDouble[maxlist];

        RDouble *Globalxlist = new RDouble[maxlist];
        RDouble *Globalylist = new RDouble[maxlist];
        RDouble *Globalzlist = new RDouble[maxlist];

        RDouble *Globalx = gridSource->GetX();
        RDouble *Globaly = gridSource->GetY();
        RDouble *Globalz = gridSource->GetZ();

        for (int m = 0; m < node_number_of_each_face[iBFace]; ++m)
        {
            //! get node index.
            index[m] = face2node[Globalpos + m];
        }
        Globalpos += node_number_of_each_face[iBFace];

        int nlist = node_number_of_each_face[iBFace];
        pindex.resize(nlist);

        GetFaceCoorList(index, nlist, xlist, ylist, zlist, Globalx, Globaly, Globalz);

        //! Correct coordinates to match target coordinate.
        if (periodicType == TRANSLATIONAL_PERIODICITY)
        {
            if (bcRegion->GetBCName() == "Periodic_up")
            {
                for (int m = 0; m < node_number_of_each_face[iBFace]; ++m)
                {
                    Globalxlist[m] = xlist[m] + translationLength[0];
                    Globalylist[m] = ylist[m] + translationLength[1];
                    Globalzlist[m] = zlist[m] + translationLength[2];
                }
            }
            else
            {
                for (int m = 0; m < node_number_of_each_face[iBFace]; ++m)
                {
                    Globalxlist[m] = xlist[m] - translationLength[0];
                    Globalylist[m] = ylist[m] - translationLength[1];
                    Globalzlist[m] = zlist[m] - translationLength[2];
                }
            }
        }
        else if (periodicType == ROTATIONAL_PERIODICITY)
        {
            if (bcRegion->GetBCName() == "Periodic_down")
            {
                for (int m = 0; m < node_number_of_each_face[iBFace]; ++m)
                {
                    RDouble rotPoint[3] = { 0.0, 0.0, 0.0 };
                    rotPoint[0] = xlist[m];
                    rotPoint[1] = ylist[m] * cos(2.0 * PI - rotationAngle) - zlist[m] * sin(2.0 * PI - rotationAngle);
                    rotPoint[2] = ylist[m] * sin(2.0 * PI - rotationAngle) + zlist[m] * cos(2.0 * PI - rotationAngle);
                    Globalxlist[m] = rotPoint[0];
                    Globalylist[m] = rotPoint[1];
                    Globalzlist[m] = rotPoint[2];
                }
            }
            if (bcRegion->GetBCName() == "Periodic_up")
            {
                for (int m = 0; m < node_number_of_each_face[iBFace]; ++m)
                {
                    RDouble rotPoint[3] = { 0.0, 0.0, 0.0 };
                    rotPoint[0] = xlist[m];
                    rotPoint[1] = ylist[m] * cos(rotationAngle) - zlist[m] * sin(rotationAngle);
                    rotPoint[2] = ylist[m] * sin(rotationAngle) + zlist[m] * cos(rotationAngle);
                    Globalxlist[m] = rotPoint[0];
                    Globalylist[m] = rotPoint[1];
                    Globalzlist[m] = rotPoint[2];
                }
            }
        }

        if (IsInterface(bcType))
        {
            if (bcRegion->GetBCName() == "Periodic_up" || bcRegion->GetBCName() == "Periodic_down")
            {
                int localIndex = -1;
                //! Find local face index for iBFace.
                MatchPeriodicFace(grids, iZone, iBFace, xlist, ylist, zlist, localIndex);

                //! Build periodic boundary link for faces that locate in part grid.

                if (localIndex < i_nBFace && localIndex > -1)
                {
                    int Targetzone;
                    int TargetFace = -1;

                    for (Targetzone = 0; Targetzone < maxProcessor; Targetzone++)
                    {
                        //! Use corrected coordinate to find target face index.
                        MatchPeriodicFace(grids, Targetzone, iBFace, Globalxlist, Globalylist, Globalzlist, TargetFace);

                        //! If find target face index in target zone, stop.
                        if (TargetFace >= 0)
                        {
                            //! Due to the periodic boundary, left cell of TargetFace take the place of right cell of iBFace. 
                            UnstructGrid *Targetgrid = UnstructGridCast(grids[Targetzone]);
                            int *left_cell_of_face = Targetgrid->GetLeftCellOfFace();
                            int *right_cell_of_face = gridpart->GetRightCellOfFace();

                            int Targetcell = left_cell_of_face[TargetFace];

                            //! Manually set neighbors for periodic boundary.
                            interFaceDirection[ii_local] = 1;
                            interFace2ZoneID[ii_local] = static_cast<int>(Targetzone);
                            interFace2CellID[ii_local] = cell_mark[Targetcell];
                            interFace2InterFaceID[ii_local] = TargetFace;
                            interFace2BoundaryFace[ii_local] = localIndex;

                            isNeighborInterfaceFound[ii_local] = true;
                            ii_local++;
                            break;
                        }
                    }
                }
            }
        }
    }
}

//! Find FaceIndex in Targetzone.
//! iBFace: global face index.
//! xlist, ylist, zlist: node cooridinate of the face.
//! FaceIndex: local face index.
void Pre_UnstructGridPartition::MatchPeriodicFace(Grid **grids, int Targetzone, int iBFace, RDouble *xlist, RDouble *ylist, RDouble *zlist, int &FaceIndex)
{
    UnstructGrid *Targetgrid = UnstructGridCast(grids[Targetzone]);
    UnstructBCSet *TargetBCSet = Targetgrid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = TargetBCSet->GetBCFaceInfo();
    int nBoundFace = Targetgrid->GetNBoundFace();
    int *face2node = Targetgrid->GetFace2Node();
    int *node_number_of_each_face = gridSource->GetNodeNumberOfEachFace();

    int nodepos = 0;
    int count = 0;

    for (int iFace = 0; iFace < nBoundFace; iFace++)
    {
        UnstructBC *TargetbcRegion = TargetBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        string TargetbcName = TargetbcRegion->GetBCName();
        int bcType = TargetbcRegion->GetBCType();
        int *Target_node_number_of_each_face = Targetgrid->GetNodeNumberOfEachFace();

        if (!IsInterface(bcType))
        {
            nodepos += Target_node_number_of_each_face[iFace];
            continue;
        }

        if (TargetbcName != "Periodic_up" && TargetbcName != "Periodic_down" )
        {
            nodepos += Target_node_number_of_each_face[iFace];
            continue;
        }

        int Target_maxlist = 0;

        Target_maxlist = MAX(Target_maxlist, Target_node_number_of_each_face[iFace]);
        vector<int> Target_index(Target_maxlist), Target_pindex(Target_maxlist);

        RDouble *Targetxlist = new RDouble[Target_maxlist];
        RDouble *Targetylist = new RDouble[Target_maxlist];
        RDouble *Targetzlist = new RDouble[Target_maxlist];

        RDouble *Targetx = Targetgrid->GetX();
        RDouble *Targety = Targetgrid->GetY();
        RDouble *Targetz = Targetgrid->GetZ();

        for (int m = 0; m < Target_node_number_of_each_face[iFace]; ++m)
        {
            //! get node index
            Target_index[m] = face2node[nodepos + m];
        }
        nodepos += Target_node_number_of_each_face[iFace];

        int Target_nlist = Target_node_number_of_each_face[iFace];
        Target_pindex.resize(Target_nlist);

        GetFaceCoorList(Target_index, Target_nlist, Targetxlist, Targetylist, Targetzlist, Targetx, Targety, Targetz);
        //! every node coordinate of iBFace should equal to iFace.
        for (int m = 0; m < node_number_of_each_face[iBFace]; ++m)
        {
            for (int n = 0; n < Target_node_number_of_each_face[iFace]; ++n)
            {
                RDouble errx, erry, errz;

                //! now tolerance is set to 1e-6.
                errx = xlist[m] - Targetxlist[n];
                erry = ylist[m] - Targetylist[n];
                errz = zlist[m] - Targetzlist[n];

                errx = sqrt(errx * errx);
                erry = sqrt(erry * erry);
                errz = sqrt(errz * errz);

                if (errx<1e-6 && erry<1e-6 && errz<1e-6)
                {
                    count++;
                    break;
                }
            }

            if (m == node_number_of_each_face[iBFace] - 1)
            {
                if (count == node_number_of_each_face[iBFace])
                {
                    FaceIndex = iFace;
                }
                //! Initial count after finishing search for one face. 
                count = 0;
            }
        }
    }
}

//! interFace2interface data: turn faceindex into facecount.
void Pre_UnstructGridPartition::ResetFace2Face(Grid **grids)
{
    int maxFace = 0;
    for (int Zone = 0; Zone < maxProcessor; Zone++)
    {
        UnstructGrid *iGrid = UnstructGridCast(grids[Zone]);
        InterfaceInfo *iinterfaceInfo = iGrid->GetInterfaceInfo();
        int inIFace = iinterfaceInfo->GetNIFace();

        if (maxFace < inIFace)
        {
            maxFace = inIFace;
        }
    }

    //! Define two-dimensional array to store "interFace2InterFaceID[ii_local]".
    int **Face2Face = new int * [maxProcessor];
    for (int iProcessor = 0; iProcessor < maxProcessor; iProcessor++)
    {
        Face2Face[iProcessor] = new int[maxFace];
    }

    for (int ZoneNumber = 0; ZoneNumber < maxProcessor; ZoneNumber++)
    {
        UnstructGrid *Grid_ori = UnstructGridCast(grids[ZoneNumber]);
        UnstructBCSet *BCSet_ori = Grid_ori->GetUnstructBCSet();
        int *bcRegionIDofBCFace_ori = BCSet_ori->GetBCFaceInfo();

        InterfaceInfo *interfaceInfo_ori = Grid_ori->GetInterfaceInfo();
        int nIFace_ori = interfaceInfo_ori->GetNIFace();
        int * Face2NeighborZone_ori = interfaceInfo_ori->GetInterFace2ZoneID();
        int * Face2NeighborFace_ori = interfaceInfo_ori->GetInterFace2InterFaceID();
        int * Face2BoundaryFace_ori = interfaceInfo_ori->GetInterFace2BoundaryFace();

        for (int FaceNumber = 0; FaceNumber < nIFace_ori; FaceNumber++)
        {
            int iBoundaryFace = Face2BoundaryFace_ori[FaceNumber];
            int iNeighborFace = Face2NeighborFace_ori[FaceNumber];
            int NeighborZone = Face2NeighborZone_ori[FaceNumber];

            UnstructBC *bcRegion_ori = BCSet_ori->GetBCRegion(bcRegionIDofBCFace_ori[iBoundaryFace]);
            string bcName_ori = bcRegion_ori->GetBCName();

            if (bcName_ori == "Periodic_up" || bcName_ori == "Periodic_down")
            {
                UnstructGrid *Grid_target = UnstructGridCast(grids[NeighborZone]);
                InterfaceInfo *interfaceInfo_target = Grid_target->GetInterfaceInfo();

                int nIFace_target = interfaceInfo_target->GetNIFace();
                int *Face2BoundaryFace_target = interfaceInfo_target->GetInterFace2BoundaryFace();

                //! Store interface2interface data in Face2Face for periodic face.
                for (int iTargetFace = 0; iTargetFace < nIFace_target; iTargetFace++)
                {
                    if (Face2BoundaryFace_target[iTargetFace] == iNeighborFace)
                    {
                        Face2Face[ZoneNumber][FaceNumber] = iTargetFace;
                    }
                }
            }
        }
    }

    //! Store interface2interface data in interfaceInfo for periodic face.
    for (int ZoneNumber = 0; ZoneNumber < maxProcessor; ZoneNumber++)
    {
        UnstructGrid *i_Grid = UnstructGridCast(grids[ZoneNumber]);
        InterfaceInfo *i_interfaceInfo = i_Grid->GetInterfaceInfo();

        UnstructBCSet *i_BCSet = i_Grid->GetUnstructBCSet();
        int *i_bcRegionIDofBCFace = i_BCSet->GetBCFaceInfo();

        int i_nIFace = i_interfaceInfo->GetNIFace();
        int *i_interFace2InterFaceID = i_interfaceInfo->GetInterFace2InterFaceID();
        int *i_Face2BoundaryFace = i_interfaceInfo->GetInterFace2BoundaryFace();

        for (int FaceNumber = 0; FaceNumber < i_nIFace; FaceNumber++)
        {
            int iBoundaryFace = i_Face2BoundaryFace[FaceNumber];
            UnstructBC *i_bcRegion = i_BCSet->GetBCRegion(i_bcRegionIDofBCFace[iBoundaryFace]);
            string i_bcName = i_bcRegion->GetBCName();

            if (i_bcName == "Periodic_up" || i_bcName == "Periodic_down")
            {
                i_interFace2InterFaceID[FaceNumber] = Face2Face[ZoneNumber][FaceNumber];
            }
        }
    }

    //! DeAllocate.
    for (int iProcessor = 0; iProcessor < maxProcessor; iProcessor++)
    {
        delete []Face2Face[iProcessor];
    }
    delete []Face2Face;
}

//! For each interface in the current zone,\n
//! Find out the interface-index in the neighbor zone.
void Pre_UnstructGridPartition::ReconInterfaces(Grid **grids, int iZone)
{
    if (!grids[iZone]) return;

    InterfaceInfo *interfaceInfo = grids[iZone]->GetInterfaceInfo();
    if (!interfaceInfo) return;

    int *interFace2ZoneID      = interfaceInfo->GetInterFace2ZoneID();
    int *interFace2InterFaceID = interfaceInfo->GetInterFace2InterFaceID();

    UnstructGrid *grid = UnstructGridCast(grids[iZone]);

    int nBoundFace = grid->GetNBoundFace();
    int nIFace = interfaceInfo->GetNIFace();
    int nphysi = nBoundFace - nIFace;    //! number of physical boundary = total number boundary face - number of interface.

    int *left_cell_of_face  = grid->GetLeftCellOfFace();
    int *right_cell_of_face = grid->GetRightCellOfFace();

    bool *isNeighborInterfaceFound = interfaceInfo->GetIsNeighborInterfaceFound();

    int nbz, le, re, ie;

    //! Interfaces count from 0.
    for (int iFace = 0; iFace < nIFace; ++ iFace)
    {
        if (isNeighborInterfaceFound[iFace] == true)
        {
            //! The neighbor of this face has been found and paired, so skip it.
            continue;
        }

        //! Neighbor zone index.
        nbz = interFace2ZoneID[iFace];
        le  = left_cell_of_face [iFace + nphysi];
        re  = right_cell_of_face[iFace + nphysi];
        ie  = MAX(le, re);

        if (nbz >= iZone)
        {
            UnstructGrid *grid_t = UnstructGridCast(grids[nbz]);
            if (!FindMatchPart(grid_t, iZone, iFace, ie, interFace2InterFaceID))
            {
                TK_Exit::ExceptionExit("!FindMatchPart\n");
            }
            else
            {
                isNeighborInterfaceFound[iFace] = true;
            }
        }
    }
}

void Pre_UnstructGridPartition::SetCoorPart(UnstructGrid *gridpart)
{
    int nTotalNode = gridpart->GetNTotalNode();
    RDouble *x = new RDouble[nTotalNode];
    RDouble *y = new RDouble[nTotalNode];
    RDouble *z = new RDouble[nTotalNode];

    gridpart->SetX(x);
    gridpart->SetY(y);
    gridpart->SetZ(z);

    int nTNodeSrc = gridSource->GetNTotalNode();
    RDouble *x_src = gridSource->GetX();
    RDouble *y_src = gridSource->GetY();
    RDouble *z_src = gridSource->GetZ();

    int count = 0;
    for (int iNode = 0; iNode < nTNodeSrc; ++ iNode)
    {
        if (node_mark[iNode] > -1)
        {
            x[count] = x_src[iNode];
            y[count] = y_src[iNode];
            z[count] = z_src[iNode];
            ++ count;
        }
    }

    if (count != nTotalNode)
    {
        cout << "error in SetCoorPart\n";
    }
}

void Pre_UnstructGridPartition::SetNodeNumberOfEachFaceAndFace2Node(UnstructGrid *gridpart)
{
    int *node_number_of_each_face_src = gridSource->GetNodeNumberOfEachFace();
    int *face2node_src = gridSource->GetFace2Node();

    int nTotalFace = gridpart->GetNTotalFace();
    int *node_number_of_each_face = new int[nTotalFace];
    gridpart->SetNodeNumberOfEachFace(node_number_of_each_face);

    //! Extract node_number_of_each_face from the original unstructured grid.
    long long int count = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int faceid = face2all[iFace];
        node_number_of_each_face[iFace] = node_number_of_each_face_src[faceid];
        count += node_number_of_each_face[iFace];
    }

    int *face2node = new int[count];
    gridpart->SetFace2Node(face2node);

    long long int *face2NodeSubscript = gridSource->GetFace2NodeSubscript();

    count = 0;
    long long int jnode;
    jnode = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        //! faceid is the index of global face.
        int faceid = face2all[iFace];
        for (jnode = face2NodeSubscript[faceid]; jnode < face2NodeSubscript[faceid + 1]; ++ jnode)
        {
            face2node[count ++] = node_mark[face2node_src[jnode]];
        }
    }
}

void Pre_UnstructGridPartition::SetFace2CellAndBC(int iZone, UnstructGrid *gridpart)
{
    int nBFace_src = gridSource->GetNBoundFace();
    int *left_cell_of_face_src  = gridSource->GetLeftCellOfFace();
    int *right_cell_of_face_src = gridSource->GetRightCellOfFace();

    int nTotalFace = gridpart->GetNTotalFace();
    int nBoundFace = gridpart->GetNBoundFace();

    //! Face cell.
    int *left_cell_of_face  = new int[nTotalFace]();
    int *right_cell_of_face = new int[nTotalFace]();
    gridpart->SetLeftCellOfFace (left_cell_of_face);
    gridpart->SetRightCellOfFace(right_cell_of_face);

    UnstructBCSet **bcr = new UnstructBCSet * [nBoundFace];
    gridpart->SetBCRecord(bcr);

    UnstructBCSet *unstructBCSetSrc = gridSource->GetUnstructBCSet();
    int *bcRegionIDofBCFaceSrc = unstructBCSetSrc->GetBCFaceInfo();

    vector <string> bcNameList;
    vector <int> bcTypeList;
    set<string> bcNameSet;

    int le_src = 0, re_src = 0, le = 0, re = 0, bctype = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int faceid = face2all[iFace];

        le_src = left_cell_of_face_src [faceid];
        re_src = right_cell_of_face_src[faceid];

        string bcName = " ";
        if (faceid < nBFace_src)
        {     
            UnstructBC *bcRegion = unstructBCSetSrc->GetBCRegion(bcRegionIDofBCFaceSrc[faceid]);

            re = - 1;
            le = cell_mark[le_src];
            bctype = bcRegion->GetBCType();
            bcName = bcRegion->GetBCName();
        }
        else
        {
            bctype = -1;
            bcName = "Interface";
            //! int face.
            if (part[le_src] == iZone)
            {
                le = cell_mark[le_src];
                re = - 1;
            }
            else if (part[re_src] == iZone)
            {
                re = cell_mark[re_src];
                le = - 1;
            }
            else
            {
                cout << "error in SetLinkPart\n";
            }
        }

        bcTypeList.push_back(bctype);
        bcNameList.push_back(bcName);
        bcNameSet.insert(bcNameList[iFace]);

        bcr[iFace] = new UnstructBCSet();
        bcr[iFace]->SetKey(bctype);
        bcr[iFace]->SetBCName(bcName);

        left_cell_of_face [iFace] = le;
        right_cell_of_face[iFace] = re;
    }
    int nBCRegionUnstruct = static_cast<int>(bcNameSet.size());
    gridpart->CreateUnstructBCSet(nBCRegionUnstruct);

    UnstructBCSet *unstructBCSet = gridpart->GetUnstructBCSet();
    int *bcRegionIDofBCFace = new int[nBoundFace];
    set<string>::iterator iter;

    int count = 0;
    for (iter = bcNameSet.begin(); iter != bcNameSet.end(); iter++)
    {
        UnstructBC *unstructBC = new UnstructBC(count);
        unstructBCSet->SetBCRegion(count, unstructBC);
        unstructBC->SetBCName(*iter);

        for (int iFace = 0; iFace < nBoundFace; iFace++)
        {
            if (bcNameList[iFace] == unstructBC->GetBCName())
            {
                unstructBC->SetBCType(bcTypeList[iFace]);
                bcRegionIDofBCFace[iFace] = count;
                vector<int> *faceIndex = unstructBC->GetFaceIndex();
                faceIndex->push_back(iFace);
            }

        }
        count++;
    }
    unstructBCSet->SetBCFaceInfo(bcRegionIDofBCFace);

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int faceid = face2all[iFace];
        le_src = left_cell_of_face_src[faceid];
        re_src = right_cell_of_face_src[faceid];

        left_cell_of_face [iFace] = cell_mark[le_src];
        right_cell_of_face[iFace] = cell_mark[re_src];
    }
}

void Pre_UnstructGridPartition::SetLinkPart(int iZone, UnstructGrid *gridpart)
{
    SetNodeNumberOfEachFaceAndFace2Node(gridpart);
    SetFace2CellAndBC(iZone, gridpart);
}

void Pre_UnstructGridPartition::SetInterfaceInfoPart(int iZone, UnstructGrid *gridpart)
{
    if (IsParallelPartition())
    {
        SetInterfaceInfoPartForParallelPartition(iZone, gridpart);
        return;
    }

    if (npartmethod != 1) return;

    InterfaceInfo *interfaceInfo = gridpart->GetInterfaceInfo();

    int i_nIFace = interfaceInfo->GetNIFace();
    int i_nBFace = gridpart->GetNBoundFace();

    int nTotalFace = gridSource->GetNTotalFace();
    int nBoundFace = gridSource->GetNBoundFace();
    int *left_cell_of_face  = gridSource->GetLeftCellOfFace();
    int *right_cell_of_face = gridSource->GetRightCellOfFace();

    int le, re;

    int *interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
    int *interFace2CellID       = interfaceInfo->GetInterFace2CellID();
    int *interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();
    int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();
    int *interFaceDirection     = interfaceInfo->GetInterFaceDirection();

    int nphysi = i_nBFace - i_nIFace;

    //! Here,the interFace2InterFaceID temporarily stores the cells index of neighbor zone.

    int count = 0;
    //! Interfaces in the original interior field.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        //! face_mark: the corresponding relationship from global face index to local face index.
        //! The following judgment condition means that this face is the interface of the zone.
        //! the index of a face is greater than 0 and less than i_nBFace in new partition zone, and it is the internal face of the original grid.
        //! So it must be interface.
        if (face_mark[iFace] < i_nBFace && face_mark[iFace] > -1)
        {
            ++ count;
            le = left_cell_of_face [iFace];
            re = right_cell_of_face[iFace];

            //! ii_local: the index of interface counts from 0.
            int ii_local = face_mark[iFace] - nphysi;

            //! interface_interfaceidtemporarily records local index of the neighbor interface cell.
            //! It will be modified to the interface index of neighbor zone.
            if (part[le] == iZone)
            {
                interFaceDirection    [ii_local] = 1;
                interFace2ZoneID      [ii_local] = static_cast<int>(part[re]);
                interFace2CellID      [ii_local] = cell_mark[re];
                interFace2InterFaceID [ii_local] = cell_mark[re];
                interFace2BoundaryFace[ii_local] = face_mark[iFace];
            }
            else if (part[re] == iZone)
            {
                interFaceDirection    [ii_local] = - 1;
                interFace2ZoneID      [ii_local] = static_cast<int>(part[le]);
                interFace2CellID      [ii_local] = cell_mark[le];
                interFace2InterFaceID [ii_local] = cell_mark[le];
                interFace2BoundaryFace[ii_local] = face_mark[iFace];
            }
            else
            {
                TK_Exit::ExceptionExit("Error: CAN NOT find neighbor interface cell!");
            }
        }
    }
}

void Pre_UnstructGridPartition::SetInterpointInfoPart(int iZone, UnstructGrid *gridpart)
{
    if (npartmethod != 1) return;

    gridpart->computeBoundaryPointLabel();
    map<int, int> boundaryPointLabel_gridpart = gridpart->GetBoundaryPointLabel();
    map<int, int> boundaryPointLabel_gridSource = gridSource->GetBoundaryPointLabel();

    int i_nTNode = gridpart->GetNTotalNode();

    int i_nIPoint = 0;

    //! Compute the number of interpoints in current zone.
    for (int i = 0; i < i_nTNode; i ++)
    {
        int numberofnode2all = node2all[i];
        i_nIPoint += pointPosition2zoneID[numberofnode2all+1] - pointPosition2zoneID[numberofnode2all] -1;
    }

    if (i_nIPoint == 0)
    {
        return;
    }

    InterpointInformation *interpointInformation = new InterpointInformation(i_nIPoint);
    gridpart->SetInterpointInfo(interpointInformation);

    //! Only partition procedure need the data of isNeighborInterpointFound, so new it here.
    bool *isNeighborInterpointFound = new bool [i_nIPoint];
    PHSPACE::SetField(isNeighborInterpointFound, false, i_nIPoint);
    interpointInformation->SetIsNeighborInterpointFound(isNeighborInterpointFound);

    int *interPoint2ZoneID = interpointInformation->GetInterPoint2ZoneID();
    int *interPoint2InterPointID = interpointInformation->GetInterPoint2InterPointID();
    int *interPoint2GlobalPoint = interpointInformation->GetInterPoint2GlobalPoint();
    int *cellNumberOfInterPoint = interpointInformation->GetCellNumberOfInterPoint();
    int *totalZonesOfInterPoint = interpointInformation->GetTotalZonesOfInterPoint();
    int *labelOfInterPoint = interpointInformation->GetLabelOfInterPoint();

    int count = 0;
    for (int iNode = 0; iNode < i_nTNode; iNode ++)
    {
        int numberofnode2all = node2all[iNode];
        if (pointPosition2zoneID[numberofnode2all+1] - pointPosition2zoneID[numberofnode2all] == 1) continue;
        for (int i = pointPosition2zoneID[numberofnode2all]; i < pointPosition2zoneID[numberofnode2all+1]; i ++)
        {
            if (point2ZoneID[i] != iZone)
            {
                interPoint2ZoneID[count] = point2ZoneID[i];
                interPoint2InterPointID[count] = node2all[iNode];    //! Temporarily to stores the index of global points,it is modified to the index of interpoints later.
                interPoint2GlobalPoint[count] = iNode;
                cellNumberOfInterPoint[count] = cellNumberForNodeValue[numberofnode2all];
                totalZonesOfInterPoint[count] = pointPosition2zoneID[numberofnode2all+1] - pointPosition2zoneID[numberofnode2all];
                if (boundaryPointLabel_gridpart[iNode] != boundaryPointLabel_gridSource[numberofnode2all])
                {
                    labelOfInterPoint[count] = 0;
                }
                else
                {
                    labelOfInterPoint[count] = 1;
                }
                count ++;
            }
        }
    }
}

void Pre_UnstructGridPartition::ModifyInterpointInfoPart(Grid **grids)
{
    for (int iproc = 0; iproc < maxProcessor; ++ iproc)
    {
        UnstructGrid *gridpart = static_cast<UnstructGrid *>(grids[iproc]);
        InterpointInformation *interpointInformation = gridpart->GetInterpointInfo();

        int i_nIPoint = interpointInformation->GetNumberOfInterpoints();
        int *interPoint2ZoneID = interpointInformation->GetInterPoint2ZoneID();
        int *interPoint2InterPointID = interpointInformation->GetInterPoint2InterPointID();

        for (int iPoint = 0; iPoint < i_nIPoint; iPoint ++)
        {
            int neighbor_zone = interPoint2ZoneID[iPoint];
            if (iproc > neighbor_zone) continue;
            UnstructGrid *gridpart_neighbor = static_cast<UnstructGrid *>(grids[neighbor_zone]);
            int i_nIPoint_neightbor = gridpart_neighbor->GetInterpointInfo()->GetNumberOfInterpoints();
            int *interPoint2InterPointID_neightbor = gridpart_neighbor->GetInterpointInfo()->GetInterPoint2InterPointID();
            int *interPoint2ZoneID_neightbor = gridpart_neighbor->GetInterpointInfo()->GetInterPoint2ZoneID();
            
            for (int jPoint = 0; jPoint <i_nIPoint_neightbor; jPoint ++)
            {
                if (interPoint2InterPointID[iPoint] == interPoint2InterPointID_neightbor[jPoint] &&
                    interPoint2ZoneID_neightbor[jPoint] == iproc)
                {
                    interPoint2InterPointID[iPoint] = jPoint;
                    interPoint2InterPointID_neightbor[jPoint] = iPoint;
                    break;
                }
            }
        }
    }
}

void Pre_UnstructGridPartition::SetInterfaceInfoPartForParallelPartition(int iZone, UnstructGrid *gridpart)
{
    if (npartmethod != 1) return;

    InterfaceInfo *interfaceInfo = gridpart->GetInterfaceInfo();

    int i_nIFace = interfaceInfo->GetNIFace();
    int i_nBFace = gridpart->GetNBoundFace();

    int nTotalFace = gridSource->GetNTotalFace();
    int nBoundFace = gridSource->GetNBoundFace();
    int *left_cell_of_face  = gridSource->GetLeftCellOfFace();
    int *right_cell_of_face = gridSource->GetRightCellOfFace();

    int le, re;

    int *interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
    int *interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();
    int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();
    int *interFaceDirection     = interfaceInfo->GetInterFaceDirection();

    int nphysi = i_nBFace - i_nIFace;

    //! Set default value for parallel partition.
    for (int iFace = 0; iFace < i_nIFace; ++ iFace)
    {
        interFace2BoundaryFace[iFace] = -1;
        interFace2InterFaceID[iFace]  = -1;
        interFace2ZoneID[iFace]       = -1;
    }

    //! Here,the interFace2InterFaceID temporarily stores the cells index of neighbor zone.

    int count = 0;
    //! Interfaces in the original interior field.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        //! face_mark: the corresponding relationship from global face index to local face index.
        //! The following judgment condition means that this face is the interface of the zone.
        //! the index of a face is greater than 0 and less than i_nBFace in new partition zone, and it is the internal face of the original grid.
        //! So it must be interface.
        if (face_mark[iFace] < i_nBFace && face_mark[iFace] > -1)
        {
            ++ count;
            le = left_cell_of_face [iFace];
            re = right_cell_of_face[iFace];

            //! ii_local: the index of interface counts from 0.
            int ii_local = face_mark[iFace] - nphysi;

            //! interface_interfaceidtemporarily records local index of the neighbor interface cell.
            //! It will be modified to the interface index of neighbor zone.
            if (part[le] == iZone)
            {
                interFaceDirection    [ii_local] = 1;
                interFace2ZoneID      [ii_local] = static_cast<int>(part[re]);
                interFace2InterFaceID [ii_local] = cell_mark[re];
                interFace2BoundaryFace[ii_local] = face_mark[iFace];
            }
            else if (part[re] == iZone)
            {
                interFaceDirection    [ii_local] = -1;
                interFace2ZoneID      [ii_local] = static_cast<int>(part[le]);
                interFace2InterFaceID [ii_local] = cell_mark[le];
                interFace2BoundaryFace[ii_local] = face_mark[iFace];
            }
            else
            {
                TK_Exit::ExceptionExit("Error: CAN NOT find neighbor interface cell!");
            }
        }
    }
}

void Pre_UnstructGridPartition::ComputeNodeMark(UnstructGrid *gridpart)
{
    int nTotalFace = gridSource->GetNTotalFace();
    int nTotalNode = gridSource->GetNTotalNode();
    int *node_number_of_each_face  = gridSource->GetNodeNumberOfEachFace();
    int *face2node = gridSource->GetFace2Node();

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        node_mark[iNode] = NOT_IN_THIS_ZONE;
    }

    //! Set all the node_mark belong to iZone to -1.
   long long int count = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        if (face_mark[iFace] > -1)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                node_mark[face2node[count ++]] = IN_THIS_ZONE;
            }
        }
        else
        {
            count += node_number_of_each_face[iFace];
        }
    }

    int i_nTNode = 0;
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (IsBelongtoThisZone(node_mark[iNode]))
        {
            node_mark[iNode] = i_nTNode ++;
        }
    }

    gridpart->SetNTotalNode(i_nTNode);
}

void Pre_UnstructGridPartition::ComputeFaceMark(int iZone, UnstructGrid *gridpart)
{
    int nTotalCell = gridSource->GetNTotalCell();
    int nTotalFace = gridSource->GetNTotalFace();
    int nBoundFace = gridSource->GetNBoundFace();
    int *left_cell_of_face  = gridSource->GetLeftCellOfFace();
    int *right_cell_of_face = gridSource->GetRightCellOfFace();
    
    UnstructBCSet *unstructBCSet = gridSource->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    int le = 0, re = 0;

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        face_mark[iFace] = NOT_IN_THIS_ZONE;
    }

    //! The face belong to this zone is set to be 'IN_THIS_ZONE'.
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        le = left_cell_of_face [iFace];
        re = right_cell_of_face[iFace];

        if (part[le] == iZone)
        {
            face_mark[iFace] = IN_THIS_ZONE;
        }
        else if (re < nTotalCell && part[re] == iZone)
        {
            face_mark[iFace] = IN_THIS_ZONE;
        }
    }

    int i_nTFace = 0;
    int i_nIFace = 0;

    //! BC faces.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        if (IsBelongtoThisZone(face_mark[iFace]))
        {
            face_mark[iFace] = i_nTFace ++;

            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);

            int bcType = bcRegion->GetBCType();
            if (bcType == PHENGLEI::INTERFACE)
            {
                //! Parallel partition.
                interfaceIndexInBC_send[iFace] = i_nIFace;
                i_nIFace ++;
            }
        }
    }

    nIFaceInOrgBC[iZone] = i_nIFace;

    //! Interior interfaces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        if (IsBelongtoThisZone(face_mark[iFace]))
        {
            le = left_cell_of_face [iFace];
            re = right_cell_of_face[iFace];
            if (part[le] != part[re])
            {
                face_mark[iFace] = i_nTFace ++;
                i_nIFace ++;
            }
        }
    }

    int i_nBFace = i_nTFace;

    //! Interior NON-interface.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        if (face_mark[iFace] < 0)
        {
            //! Has not been set.
            if (IsBelongtoThisZone(face_mark[iFace]))
            {
                le = left_cell_of_face [iFace];
                re = right_cell_of_face[iFace];

                if (part[le] == part[re])
                {
                    face_mark[iFace] = i_nTFace ++;
                }
            }
        }
    }

    gridpart->SetNTotalFace(i_nTFace);
    gridpart->SetNBoundFace(i_nBFace);
    if (i_nIFace > 0)
    {
        InterfaceInfo *interfaceInfo = new InterfaceInfo(i_nIFace);
        gridpart->SetInterfaceInfo(interfaceInfo);

        //! Only partition procedure need the data of isNeighborInterfaceFound,
        //! so new it here.
        bool *isNeighborInterfaceFound = new bool [i_nIFace];
        PHSPACE::SetField(isNeighborInterfaceFound, false, i_nIFace);
        interfaceInfo->SetIsNeighborInterfaceFound(isNeighborInterfaceFound);
    }
}

int Pre_UnstructGridPartition::GetNTotalCell(int iZone, idx_t *part)
{
    int nTotalCell = gridSource->GetNTotalCell();
    int count = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (part[iCell] == iZone)
        {
            count ++;
        }
    }
    return count;
}

int * Pre_UnstructGridPartition::GetNode2All(UnstructGrid *gridpart)
{
    int i_nTNode = gridpart->GetNTotalNode();

    int *node2all = new int[i_nTNode];    //! node2all, record the connection relationship from nodes in partition zone to all nodes.

    int nTotalNode = gridSource->GetNTotalNode();

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (node_mark[iNode] > -1)
        {
            node2all[node_mark[iNode]] = iNode;
        }
    }

    return node2all;
}

int * Pre_UnstructGridPartition::GetFace2All(UnstructGrid *gridpart)
{
    int i_nTFace = gridpart->GetNTotalFace();

    int *face2all = new int[i_nTFace];    //! face2all,record the connection relationship from faces in partition zone to all faces(after ordering).

    int nTotalFace = gridSource->GetNTotalFace();

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        if (face_mark[iFace] > -1)
        {
            face2all[face_mark[iFace]] = iFace;
        }
    }

    return face2all;
}

int * Pre_UnstructGridPartition::GetCell2All(int iZone, UnstructGrid *gridpart)
{
    int nTotalCell = gridpart->GetNTotalCell();

    int *cell2all = new int[nTotalCell];    //! cell2all,record the connection relationship from cells in partition zone to all cells.

    int nTotalCell_All = gridSource->GetNTotalCell();

    int count = 0;
    for (int iCell = 0; iCell < nTotalCell_All; ++ iCell)
    {
        if (part[iCell] == iZone)
        {
            cell2all[count] = iCell;
            ++ count;
        }
    }

    return cell2all;
}

void Pre_UnstructGridPartition::ComputePoint2ZoneID()
{
    int *n2c = gridSource->GetNode2Cell();
    int nTotalCell = gridSource->GetNTotalCell();
    int *nCPN = gridSource->GetCellNumberOfEachNode();
    int nTotalNode = gridSource->GetNTotalNode();

    set<int> *point2ZoneIDSet = new set<int>[nTotalNode];
    vector <int> pointsForPartition;
    vector <int> point2ZoneID_Vector;

    pointPosition2zoneID = new int [nTotalNode + 1];
    pointPosition2zoneID[0] = 0;

    long long int count = 0;
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        //! node2ZoneSet is used to store and filter the repeated corresponding relationship from nodes to zone.
        for (int jcell = 0; jcell < nCPN[iNode]; ++ jcell)
        {
            if (!(n2c[count + jcell] < nTotalCell))
            {
                continue;
            }
                point2ZoneIDSet[iNode].insert(static_cast<int>(part[n2c[count + jcell]]));
        }
        count += nCPN[iNode];

        //! The size of the corresponding relationship from nodes to zone can not be determined,we can not creat static array to store.
        //! So it can only be stored with vector.
        for (set<int>::iterator iter = point2ZoneIDSet[iNode].begin(); iter != point2ZoneIDSet[iNode].end(); ++ iter)
        {
            pointsForPartition.push_back(iNode);
            point2ZoneID_Vector.push_back((*iter));
        }
        pointPosition2zoneID[iNode+1] = pointPosition2zoneID[iNode] + static_cast<int>(point2ZoneIDSet[iNode].size());
    }

    int pointsSize = static_cast<int>(pointsForPartition.size());
    pointNumber2zoneID = new int [pointsSize];
    point2ZoneID = new int [pointsSize];
    
    for (int iPoint =0; iPoint < pointsSize; iPoint ++)
    {
        pointNumber2zoneID[iPoint] = pointsForPartition[iPoint];
        point2ZoneID[iPoint] = point2ZoneID_Vector[iPoint];
    }

    delete [] point2ZoneIDSet;
}

void Pre_UnstructGridPartition::ComputeCellNumberForNodeValue()
{
    int nTotalNode = gridSource->GetNTotalNode();
    cellNumberForNodeValue = new int[nTotalNode];

    int nBoundFace = gridSource->GetNBoundFace();
    int nTotalFace = gridSource->GetNTotalFace();

    int *face2node = gridSource->GetFace2Node();
    int *node_number_of_each_face = gridSource->GetNodeNumberOfEachFace();

    for (int iNode = 0; iNode < nTotalNode; iNode ++)
    {
        cellNumberForNodeValue[iNode] = 0;
    }

    int pt;
    long long int nodepos = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            pt = face2node[nodepos + j];
            //! From left.
            cellNumberForNodeValue[pt] += 1;
            //! From right.
            cellNumberForNodeValue[pt] += 1;
        }
        nodepos += node_number_of_each_face[iFace];
    }

    nodepos = 0;
    UnstructBCSet **bcr = gridSource->GetBCRecord();
    int bcType = 0;

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        bcType = bcr[iFace]->GetKey();
        if (bcType == PHENGLEI::SYMMETRY)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt = face2node[nodepos + j];
                cellNumberForNodeValue[pt] = 0;
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

    nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        bcType = bcr[iFace]->GetKey();
        if (bcType == PHENGLEI::SYMMETRY)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt = face2node[nodepos + j];
                //! From left.
                cellNumberForNodeValue[pt] += 1;
                //! From right.
                cellNumberForNodeValue[pt] += 1;
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

    nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        bcType = bcr[iFace]->GetKey();
        if (bcType == PHENGLEI::FARFIELD)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt = face2node[nodepos + j];

                cellNumberForNodeValue[pt] = 0;
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

    nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        bcType = bcr[iFace]->GetKey();
        if (bcType == PHENGLEI::FARFIELD)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt = face2node[nodepos + j];
                //! From left.
                cellNumberForNodeValue[pt] += 1;
                //! From right.
                cellNumberForNodeValue[pt] += 1;
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

    nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        bcType = bcr[iFace]->GetKey();
        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt = face2node[nodepos + j];

                cellNumberForNodeValue[pt] = 0;
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

    nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        bcType = bcr[iFace]->GetKey();
        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt = face2node[nodepos + j];
                //! From left.
                cellNumberForNodeValue[pt] += 1;
                //! From right.
                cellNumberForNodeValue[pt] += 1;
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

}

void Pre_UnstructGridPartition::SetNTotalCell(int iZone, UnstructGrid *gridpart)
{
    int i_nTotalCell = GetNTotalCell(iZone, part);
    gridpart->SetNTotalCell(i_nTotalCell);
}

void Pre_UnstructGridPartition::CreateLocal2All(int iZone, UnstructGrid *gridpart)
{
    node2all = GetNode2All(gridpart);
    face2all = GetFace2All(gridpart);
    cell2all = GetCell2All(iZone, gridpart);
}

void Pre_UnstructGridPartition::ComputeGridForPartition(UnstructGrid *gridpart, int iZone)
{
    SetNTotalCell(iZone, gridpart);
    ComputeFaceMark(iZone, gridpart);
    ComputeNodeMark(gridpart);

    CreateLocal2All(iZone, gridpart);

    SetCoorPart(gridpart);
    SetLinkPart(iZone, gridpart);
    SetInterfaceInfoPart(iZone, gridpart);
    SetInterpointInfoPart(iZone, gridpart);
    SetVolumeCondition(gridpart);
    SetWallDist(gridpart);
    SetPartitionInfo(gridpart);

    FreeLocal2All();
}

void Pre_UnstructGridPartition::SetVolumeCondition(UnstructGrid *gridpart)
{
    SimpleVC *volumeCondition = gridSource->GetVolumeConditionIn();
    if (!volumeCondition)
    {
        return;
    }

    int hasVolumeCondition = 1;
    GlobalDataBase::UpdateData("hasVolumeCondition", &hasVolumeCondition, PHINT, 1);

    string vcName = volumeCondition->GetVCName();
    int vcType = volumeCondition->GetVCType();

    SimpleVC *volumeConditionPart = new SimpleVC();
    volumeConditionPart->SetVCName(vcName);
    volumeConditionPart->SetVCType(vcType);
    gridpart->SetVolumeCondition(volumeConditionPart);
}

void Pre_UnstructGridPartition::SetWallDist(UnstructGrid *gridpart)
{
    if (!partitionWallDist) return;
    
    int nTotalCell = gridpart->GetNTotalCell();
    gridpart->InitWallDist();
    
    RDouble *walldist = gridpart->GetWallDist();
    RDouble *walldist_src = gridSource->GetWallDist();

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int iCell2ALL = cell2all[iCell];
        walldist[iCell] = walldist_src[iCell2ALL];
    }
}

void Pre_UnstructGridPartition::SetPartitionInfo(UnstructGrid *gridpart)
{
    int nTotalNode = gridpart->GetNTotalNode();
    int nTotalFace = gridpart->GetNTotalFace();
    int nTotalCell = gridpart->GetNTotalCell();

    int *ordinaryNodeIndex = new int[nTotalNode];
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        ordinaryNodeIndex[iNode] = node2all[iNode];
    }

    int *ordinaryFaceIndex = new int[nTotalFace];
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        ordinaryFaceIndex[iFace] = face2all[iFace];
    }

    int *ordinaryCellIndex = new int[nTotalCell];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        ordinaryCellIndex[iCell] = cell2all[iCell];
    }

    gridpart->SetOrdinaryGridIndex(gridSource->GetZoneID());
    gridpart->SetOrdinaryNodeIndex(ordinaryNodeIndex);
    gridpart->SetOrdinaryFaceIndex(ordinaryFaceIndex);
    gridpart->SetOrdinaryCellIndex(ordinaryCellIndex);
}

void Pre_UnstructGridPartition::WriteAdditionalData(const string &out_grid_file)
{
    WriteInterpointInfo(out_grid_file);

    //! Write cell to node connection.
    WriteCellToNode(out_grid_file);

    //! Write face boundary condition information.
    WriteFaceBC(out_grid_file);

    //! Write walldist information.
    WriteWalldist(out_grid_file);
}

void Pre_UnstructGridPartition::WriteWalldist(const string &out_grid_file)
{
    if (!partitionWallDist) return;
    string walldistfile = ChangeExtensionOfFileName(out_grid_file, "wdt");

    ActionKey *actkey = new ActionKey();
    actkey->filename = (walldistfile);

    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    string finalFileName = PHSPACE::AddSymbolToFileName(walldistfile, '_', currentProcessorIndex);
    actkey->filepos = H5Fcreate(finalFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    using namespace PHMPI;
    for (int iZone = 0; iZone < maxProcessor; ++ iZone)
    {
        if (!grids[iZone]) continue;

        UnstructGrid *grid = UnstructGridCast(grids[iZone]);
        WriteWalldist(actkey, grid);
    }
    
    H5Fclose(actkey->filepos);
    delete actkey;    actkey = nullptr;
}

void Pre_UnstructGridPartition::WriteWalldist(ActionKey *actkey, UnstructGrid *gridPart)
{
    ostringstream oss;
    oss << "Grid-" << gridPart->GetZoneID();
    string dataName = oss.str();

    int nTotalCell = gridPart->GetNTotalCell();
    RDouble *walldist = gridPart->GetWallDist();

    CreateAndWriteData(actkey->filepos, dataName, 1, nTotalCell, PHDOUBLE, walldist);
}

void Pre_UnstructGridPartition::WriteCellToNode(const string &out_grid_file)
{
    if (!partitionCellToNode) return;
    string partitionedCellToNodeFileName = ChangeExtensionOfFileName(out_grid_file, "c2n");
    fstream file;
    PHSPACE::OpenFile(file, partitionedCellToNodeFileName, ios_base::out|ios_base::binary);

    PHWrite(file, maxProcessor);
    for (int iproc = 0; iproc < maxProcessor; ++ iproc)
    {
        if (!grids[iproc]) continue;

        UnstructGrid *gridPart = UnstructGridCast(grids[iproc]);
        WriteCellToNode(file, gridPart);
    }

    PHSPACE::CloseFile(file);
}

void Pre_UnstructGridPartition::WriteCellToNode(fstream &file, UnstructGrid *gridPart)
{
    VirtualFile *virtualFile = new VirtualFile(&file);
    virtualFile->BeginWriteWork();

    int *cellNodeIndexContainer  = gridPart->GetCell2Node();
    int *cellNodeNumberContainer = gridPart->GetNodeNumberOfEachCell();
    int iZone = gridPart->GetGridID()->GetIndex();

    int numberOfCells = gridPart->GetNTotalCell();

    PHWrite(virtualFile, iZone);
    PHWrite(virtualFile, numberOfCells);

    int iPosition = 0;
    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        int nodeNumber = cellNodeNumberContainer[iCell];

        PHWrite(virtualFile, nodeNumber);

        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            PHWrite(virtualFile, cellNodeIndexContainer[iPosition ++]);
        }
    }

    virtualFile->EndWriteWork();

    delete virtualFile;
}

void Pre_UnstructGridPartition::WriteInterpointInfo(const string &out_grid_file)
{
    if (IsParallelPartition() && (parallelPartitionMethod == PARMETIS_METHOD))
    {
        return;
    }

    string interpointInfoFileName = ChangeExtensionOfFileName(out_grid_file, "interpoint");
    fstream file;
    PHSPACE::OpenFile(file, interpointInfoFileName, ios_base::out|ios_base::binary);

    PHWrite(file, maxProcessor);
    for (int iproc = 0; iproc < maxProcessor; ++ iproc)
    {
        if (!grids[iproc])
        {
            continue;
        }
        UnstructGrid *gridPart = UnstructGridCast(grids[iproc]);
        WriteInterpointInfo(file, gridPart);
    }

    PHSPACE::CloseFile(file);
}

void Pre_UnstructGridPartition::WriteInterpointInfo(fstream &file, UnstructGrid *gridPart)
{
    VirtualFile *virtualFile = new VirtualFile(&file);
    virtualFile->BeginWriteWork();
    
    int iZone = gridPart->GetGridID()->GetIndex();
    
    InterpointInformation *interpointInformation = gridPart->GetInterpointInfo();
    
    int i_nIPoint = interpointInformation->GetNumberOfInterpoints();
    int *interPoint2ZoneID = interpointInformation->GetInterPoint2ZoneID();
    int *interPoint2InterPointID = interpointInformation->GetInterPoint2InterPointID();
    int *interPoint2GlobalPoint = interpointInformation->GetInterPoint2GlobalPoint();
    int *cellNumberOfInterPoint = interpointInformation->GetCellNumberOfInterPoint();
    int *totalZonesOfInterPoint = interpointInformation->GetTotalZonesOfInterPoint();
    int *labelOfInterPoint = interpointInformation->GetLabelOfInterPoint();

    PHWrite(virtualFile, iZone);
    PHWrite(virtualFile, i_nIPoint);
    
    for (int iPoint = 0; iPoint < i_nIPoint; iPoint ++)
    {
        PHWrite(virtualFile, interPoint2ZoneID[iPoint]);
        PHWrite(virtualFile, interPoint2InterPointID[iPoint]);
        PHWrite(virtualFile, interPoint2GlobalPoint[iPoint]);
        PHWrite(virtualFile, cellNumberOfInterPoint[iPoint]);
        PHWrite(virtualFile, totalZonesOfInterPoint[iPoint]);
        PHWrite(virtualFile, labelOfInterPoint[iPoint]);
    }

    virtualFile->EndWriteWork();

    delete virtualFile;
}

void Pre_UnstructGridPartition::WriteFaceBC(const string &out_grid_file)
{
    if (!partitionFaceBC) return;
    string partitionedFaceBCFileName = ChangeExtensionOfFileName(out_grid_file, "bc");
    fstream file;
    PHSPACE::OpenFile(file, partitionedFaceBCFileName, ios_base::out|ios_base::binary);

    PHWrite(file, maxProcessor);
    for (int iproc = 0; iproc < maxProcessor; ++ iproc)
    {
        if (!grids[iproc]) continue;

        UnstructGrid *gridPart = UnstructGridCast(grids[iproc]);
        WriteFaceBC(file, gridPart);
    }

    PHSPACE::CloseFile(file);
}

void Pre_UnstructGridPartition::WriteFaceBC(fstream &file, UnstructGrid *gridPart)
{
    gridPart->WriteFaceBC(file);
}

void Pre_UnstructGridPartition::ComputeGlobalToLocalCellToNodeMapping(UnstructGrid *gridpart, int iZone)
{
    if (!partitionCellToNode) return;

    int *cellToNodeGlobal       = gridSource->GetCell2Node();
    int *nodeNumberOfCellGlobal = gridSource->GetNodeNumberOfEachCell();

    int numberOfCells = gridpart->GetNTotalCell();
    
    int *nodeNumberOfCellLocal = new int [numberOfCells];
    gridpart->SetNodeNumberOfEachCell(nodeNumberOfCellLocal);

    int numberOfGlobalCells = gridSource->GetNTotalCell();
    int *globalToLocalNodeIndexMapping = this->node_mark;
    idx_t *globalCellIndexToLocalZoneIndexMapping = this->part;

    int iCount = 0;
    long long int count = 0;
    for (int iCell = 0; iCell < numberOfGlobalCells; ++ iCell)
    {
        if (globalCellIndexToLocalZoneIndexMapping[iCell] == iZone)
        {
            int nodeNumber = nodeNumberOfCellGlobal[iCell];
            nodeNumberOfCellLocal[iCount ++] = nodeNumber;
            count += nodeNumber;
        }
    }

    int *cellToNodeLocal = new int [count];
    gridpart->SetCell2Node(cellToNodeLocal);

    long long int iPositionGlobal = 0;
    count = 0;
    for (int iCell = 0; iCell < numberOfGlobalCells; ++ iCell)
    {
        int nodeNumber = nodeNumberOfCellGlobal[iCell];

        if (globalCellIndexToLocalZoneIndexMapping[iCell] == iZone)
        {
            for (int iNode = 0; iNode < nodeNumber; ++ iNode)
            {
                int globalNodeIndex = cellToNodeGlobal[iPositionGlobal + iNode];
                int nodeIndex = globalToLocalNodeIndexMapping[globalNodeIndex];
                cellToNodeLocal[count ++] = nodeIndex;
            }
        }
        iPositionGlobal += nodeNumber;
    }
    gridpart->SetCell2Node(cellToNodeLocal);
    gridpart->SetNodeNumberOfEachCell(nodeNumberOfCellLocal);
}

bool Pre_UnstructGridPartition::IsParallelPartition()
{
    //if (PHMPI::GetNumberofGlobalZones() > 1)
    if (vtxdist)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void Pre_UnstructGridPartition::FillnTCellofEachZone(int *&nTCellOfZones)
{
    WriteLogFile("Start to Fill nTotalCell of Each Zone");

    int nZones = PHMPI::GetNumberofGlobalZones();
    int myid   = PHMPI::GetCurrentProcessorID();

    int nTProcessor = PHMPI::GetNumberOfProcessor();

    if (nTProcessor == nZones)
    {
        //! Each processors has ONE zone.
        int nTCellInThisZone = gridSource->GetNTotalCell();

        int *localSizeTemp = new int [nTProcessor];
        SetField(localSizeTemp, 1, nTProcessor);
        int **nTCellListTemp = new int * [nTProcessor];

        PHMPI::FillGobalDataCollectBcast(&nTCellInThisZone, localSizeTemp, nTCellListTemp);

        for (int iProc = 0; iProc < nTProcessor; ++ iProc)
        {
            nTCellOfZones[iProc] = nTCellListTemp[iProc][0];
        }

        DelPointer2(nTCellListTemp);
        DelPointer(localSizeTemp);
    }
    else
    {
        //! Each processors has multiple zone.
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int procID = PHMPI::GetZoneProcessorID(iZone);
            int nTCellInThisZone;

            if (myid == procID)
            {
                nTCellInThisZone = GetGrid(iZone, 0)->GetNTotalCell();
            }

            PHMPI::PH_Bcast(&nTCellInThisZone, sizeof(int), procID);

            nTCellOfZones[iZone] = nTCellInThisZone;
        }
    }

    WriteLogFile("End to Fill nTCell of Each Zone");
}

LIB_EXPORT void Pre_UnstructGridPartition::DumpPartitionGrids(const string &out_grid_file)
{
    if (IsParallelPartition())
    {
        DumpParallelPartitionGrid(out_grid_file, grids, maxProcessor);
    }
    else
    {
        DumpSerialPartitionGrid(out_grid_file, grids);
    }
}

void Pre_UnstructGridPartition::DumpSerialPartitionGrid(const string &gridFileName, Grid **grids)
{
    DumpGrid(gridFileName, grids, maxProcessor);

    string partitionedGridFileName = AddSymbolToFileName(gridFileName, '_', PHMPI::GetCurrentProcessorID());

    int dumpOldGrid = GlobalDataBase::GetIntParaFromDB("dumpOldGrid");
    if (dumpOldGrid)
    {
        WriteAdditionalData(partitionedGridFileName);
    }

    WriteWalldist(gridFileName);
}

void Pre_UnstructGridPartition::DumpParallelPartitionGrid(const string &gridFileName, Grid **grids, int nBlocks)
{
    int dumpOldGrid = GlobalDataBase::GetIntParaFromDB("dumpOldGrid");
    if (!dumpOldGrid)
    {
        DumpGrid(gridFileName, grids, nBlocks);

        WriteWalldist(gridFileName);
        return;
    }

    //! Compute the actual number of parts in this processor.
    int nBlocksActual = 0;
    for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
    {
        if (grids[iBlock])
        {
            ++ nBlocksActual;
        }
    }
    
    Grid **gridsActual = new Grid * [nBlocksActual];
    nBlocksActual = 0;
    for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
    {
        if (grids[iBlock])
        {
            gridsActual[nBlocksActual] = grids[iBlock];
            ++ nBlocksActual;
        }
    }
    
    int nZones = PHMPI::GetNumberofGlobalZones();
    int nProc  = PHMPI::GetNumberOfProcessor();
    ASSERT(nProc == nZones);
    int myid = PHMPI::GetCurrentProcessorID();
    int *nPartOfEachProcessor = new int[nZones];

    if (nProc == nZones)
    {
        //! Each processors has ONE zone.
        int *localSizeTemp = new int [nProc];
        SetField(localSizeTemp, 1, nProc);
        int **nPartListTemp = new int * [nProc];

        PHMPI::FillGobalDataCollectBcast(&nBlocksActual, localSizeTemp, nPartListTemp);

        for (int iProc = 0; iProc < nProc; ++ iProc)
        {
            nPartOfEachProcessor[iProc] = nPartListTemp[iProc][0];
        }

        DelPointer2(nPartListTemp);
        DelPointer(localSizeTemp);
    }
    else
    {
        PHMPI::FillGobalData(nBlocksActual, nPartOfEachProcessor);
    }

    int *block_index = new int[nBlocksActual];
    int count = 0;
    for (int iProc = 0; iProc < myid; ++ iProc)
    {
        count += nPartOfEachProcessor[iProc];
    }
    for (int iBlock = 0; iBlock < nBlocksActual; ++ iBlock)
    {
        block_index[iBlock] = iBlock + count;
    }

    string partitionedGridFileName = AddSymbolToFileName(gridFileName, '_', PHMPI::GetCurrentProcessorID());
    DumpParallelPartitionGrid(partitionedGridFileName, gridsActual, nBlocksActual, block_index);
    WriteAdditionalData(partitionedGridFileName);

    delete [] gridsActual;
    delete [] nPartOfEachProcessor;
    delete [] block_index;
}

void Pre_UnstructGridPartition::DumpParallelPartitionGrid(const string &gridFileName, Grid **grids, int nBlocks, int *block_index_in)
{
    if (gridFileName == "") return;
    fstream file;
    PHSPACE::OpenFile(file, gridFileName, ios_base::out|ios_base::binary|ios_base::trunc);

    int *block_proc     = new int [nBlocks];
    int *block_idx      = new int [nBlocks];
    int *block_type     = new int [nBlocks];
    int *block_proc_inp = new int [nBlocks];

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        block_proc[iZone] = grids[iZone]->GetGridID()->GetIndex();  // Processor ID actually.
        block_idx [iZone] = block_index_in[iZone];
        block_type[iZone] = grids[iZone]->Type();
    }

    if (parallelPartitionMethod == METIS_OPENMP_METHOD)
    {
        //! For OpenMP case, reset the processor ID of sub-zone to the source grid's.
        int zoneID = gridSource->GetGridID()->GetIndex();
        int processorID = PHMPI::GetZoneProcessorID(zoneID);
        for (int iZone = 0; iZone < nBlocks; ++ iZone)
        {
            block_proc[iZone] = processorID;
        }
    }

    file.write(reinterpret_cast<char *>(&nBlocks), sizeof(int));
    file.write(reinterpret_cast<char *>(block_proc), nBlocks * sizeof(int));
    file.write(reinterpret_cast<char *>(block_idx), nBlocks * sizeof(int));
    file.write(reinterpret_cast<char *>(block_type), nBlocks * sizeof(int));
    
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        grids[iZone]->WriteGrid(file);
    }

    file.close();
    file.clear();

    delete [] block_idx;
    delete [] block_proc;
    delete [] block_proc_inp;
    delete [] block_type;
}

}
