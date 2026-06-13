#include "Pre_HDF5File.h"
#include "Math_BasisFunction.h"
#include "Glb_Dimension.h"
#include "IO_FileName.h"
#include "Zone.h"
#include "GridType.h"
#include "Geo_UnstructBC.h"
#include "PHIO.h"
#include "GlobalDataBase.h"
#include "TK_Parse.h"

namespace PHSPACE
{
using namespace PHMPI;

IO_HDF5Write::IO_HDF5Write(const string &gridFileName, Grid **grids, int nBlocks)
{
    this->gridFileName = gridFileName;
    this->grids = grids;
    this->nBlocks = nBlocks;

    gridsLocal = 0;
    nPartOfEachProcessor = 0;

    block_proc = 0;
    block_idx  = 0;
    block_type = 0;
    file_index = 0;
}

IO_HDF5Write::~IO_HDF5Write()
{
    if (gridsLocal)
    {
        delete [] gridsLocal;    gridsLocal = nullptr;
        delete [] nPartOfEachProcessor;    nPartOfEachProcessor = nullptr;
    }

    if (block_proc)
    {
        delete [] block_proc;    block_proc = nullptr;
        delete [] block_idx;    block_idx = nullptr;
        delete [] block_type;    block_type = nullptr;
        delete [] file_index;    file_index = nullptr;
    }
}

hid_t IO_HDF5Write::CreateFile(const string &filename)
{
    hid_t file;
    file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    return file;
}

void IO_HDF5Write::Run()
{
    hid_t file;
    int tasktype = GetTaskCode();
    int numberOfMultifile = PHSPACE::GlobalDataBase::GetIntParaFromDB("numberOfMultifile");
    CreateTotalInfo();

    int myid = PHMPI::GetCurrentProcessorID();
    if (tasktype == PHSPACE::PARTITION_GRID && numberOfMultifile > 1)
    {
        file = CreateFile(gridFileName);
    }
    else
    {
        file = CreateFile(AddSymbolToFileName(gridFileName, "_", myid));
    }

    WriteVersionInfo(file);
    WriteTotalInfo(file);

    int nBlocksLocal = nPartOfEachProcessor[myid];
    for (int iZone = 0; iZone < nBlocksLocal; ++ iZone)
    {
        ostringstream oss;
        oss << "      Dumping zone " << iZone << " ....\n";
        PrintToWindow(oss);

        if (block_type[iZone] == STRUCTGRID)
        {
            WriteStructGrid(file, gridsLocal[iZone]);
        }
        else
        {
            WriteUnstructGrid(file, gridsLocal[iZone], iZone);
        }
    }

    H5Fclose(file);
}

void IO_HDF5Write::WriteVersionInfo(hid_t loc_id)
{
    hid_t grpData;

    grpData = GetGroup(loc_id, "Version");

    RDouble Version = 2.0;
    CreateAndWriteData(grpData, "VersionData", 1, 1, PHDOUBLE, &Version);

    H5Gclose(grpData);
}

void IO_HDF5Write::CreateTotalInfo()
{
    int nBlocksLocal = 0;

    for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
    {
        if (grids[iBlock])
        {
            ++ nBlocksLocal;
        }
    }

    gridsLocal = new Grid *[nBlocksLocal]();

    nBlocksLocal = 0;
    for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
    {
        if (grids[iBlock])
        {
            gridsLocal[nBlocksLocal] = grids[iBlock];
            ++ nBlocksLocal;
        }
    }

    int nProc = PHMPI::GetNumberOfProcessor();
    int myid  = PHMPI::GetCurrentProcessorID();

    nPartOfEachProcessor = new int[nProc]();
    PH_AllGather(&nBlocksLocal, 1, nPartOfEachProcessor, 1);

    block_proc = new int[nBlocksLocal]();
    block_idx  = new int[nBlocksLocal]();
    block_type = new int[nBlocksLocal]();
    file_index = new int[nBlocksLocal]();

    int count = 0;
    for (int iProc = 0; iProc < myid; ++ iProc)
    {
        count += nPartOfEachProcessor[iProc];
    }

    for (int iBlock = 0; iBlock < nBlocksLocal; ++ iBlock)
    {
        block_proc[iBlock] = gridsLocal[iBlock]->GetGridID()->GetIndex();

        block_idx [iBlock] = gridsLocal[iBlock]->GetGridID()->GetIndex();
        if (GlobalDataBase::GetIntParaFromDB("parallelPartitionMethod") == 1)
        {
            block_idx [iBlock] = iBlock + count;
        }

        block_type[iBlock] = gridsLocal[iBlock]->Type();
        file_index[iBlock] = gridsLocal[iBlock]->GetFileIndex();
    }

    int *block_proc_inp = GetZoneProcessorID_INP();
    if (block_proc_inp)
    {
        count = 0;
        for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
        {
            if (grids[iBlock])
            {
                block_proc[count] = block_proc_inp[iBlock];
                count++;
            }
        }
    }

    /*if (parallelPartitionMethod == METIS_OPENMP_METHOD)
    {
        // For OpenMP case, reset the processor ID of sub-zone to the source grid's.
        int zoneID = gridSource->GetGridID()->GetIndex();
        int processorID = PHMPI::GetZoneProcessorID(zoneID);
        for (int iZone = 0; iZone < nBlocks; ++ iZone)
        {
            block_proc[iZone] = processorID;     
        }
    }*/
}

void IO_HDF5Write::WriteTotalInfo(hid_t loc_id)
{
    hid_t grpData;
    grpData = GetGroup(loc_id, "Information");

    int myid = PHMPI::GetCurrentProcessorID();
    int nBlocksLocal = nPartOfEachProcessor[myid];

    CreateAndWriteData(grpData, "nBlocks"   , 1, 1           , PHINT, &nBlocksLocal);
    CreateAndWriteData(grpData, "block_proc", 1, nBlocksLocal, PHINT, block_proc);
    CreateAndWriteData(grpData, "block_idx" , 1, nBlocksLocal, PHINT, block_idx);
    CreateAndWriteData(grpData, "block_type", 1, nBlocksLocal, PHINT, block_type);
    CreateAndWriteData(grpData, "file_index", 1, nBlocksLocal, PHINT, file_index);

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteUnstructGrid(hid_t loc_id, Grid *gridIn, int iZone)
{
    hid_t grpData;
    string grpName;

    ostringstream oss;
    int gridID = block_idx[iZone];
    oss << "Grid-" << gridID;
    grpName = oss.str();
    grpData = GetGroup(loc_id, grpName);

    int iDimensions[3];
    iDimensions[0] = gridIn->GetNTotalNode();
    iDimensions[1] = gridIn->GetNTotalFace();
    iDimensions[2] = gridIn->GetNTotalCell();
    CreateAndWriteData(grpData, "iDimensions", 1, 3, PHINT, iDimensions);

    WriteGridCoordinates(grpData, gridIn);

    WriteFaceTopology(grpData, gridIn);

    WriteInterfaceInfo(grpData, gridIn);

    WriteCellTopology(grpData, gridIn);

    WriteUnstrBCName(grpData, gridIn);

    WriteInterpointInfo(grpData, gridIn);

    WriteUnstrPartitionInfo(grpData, gridIn);

    WriteVolumeCondition(grpData, gridIn);

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteInterpointInfo(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;
    grpData = GetGroup(loc_id, "InterpointInfo");

    int nIPoint = 0;
    InterpointInformation *interpointInfo = gridIn->GetInterpointInfo();
    if (!interpointInfo)
    {
        CreateAndWriteData(grpData, "nIPoint", 1, 1, PHINT, &nIPoint);
        H5Gclose(grpData);
        return;
    }

    nIPoint = interpointInfo->GetNumberOfInterpoints();
    CreateAndWriteData(grpData, "nIPoint", 1, 1, PHINT, &nIPoint);

    int *interPoint2ZoneID       = interpointInfo->GetInterPoint2ZoneID();
    int *interPoint2InterPointID = interpointInfo->GetInterPoint2InterPointID();
    int *interPoint2GlobalPoint  = interpointInfo->GetInterPoint2GlobalPoint();
    int *cellNumberOfInterPoint  = interpointInfo->GetCellNumberOfInterPoint();
    int *totalZonesOfInterPoint  = interpointInfo->GetTotalZonesOfInterPoint();
    int *labelOfInterPoint       = interpointInfo->GetLabelOfInterPoint();

    CreateAndWriteData(grpData, "interPoint2ZoneID"      , 1, nIPoint, PHINT, interPoint2ZoneID);
    CreateAndWriteData(grpData, "interPoint2InterPointID", 1, nIPoint, PHINT, interPoint2InterPointID);
    CreateAndWriteData(grpData, "interPoint2GlobalPoint" , 1, nIPoint, PHINT, interPoint2GlobalPoint);
    CreateAndWriteData(grpData, "cellNumberOfInterPoint" , 1, nIPoint, PHINT, cellNumberOfInterPoint);
    CreateAndWriteData(grpData, "totalZonesOfInterPoint" , 1, nIPoint, PHINT, totalZonesOfInterPoint);
    CreateAndWriteData(grpData, "labelOfInterPoint"      , 1, nIPoint, PHINT, labelOfInterPoint);

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteStructGrid(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;
    string grpName;

    ostringstream oss;
    int gridID = gridIn->GetZoneID();
    oss << "Grid-" << gridID;
    grpName = oss.str();
    grpData = GetGroup(loc_id, grpName);

    StructGrid *grid = StructGridCast(gridIn);

    int iDimensions[3];
    iDimensions[0] = grid->GetNI();
    iDimensions[1] = grid->GetNJ();
    iDimensions[2] = grid->GetNK();
    CreateAndWriteData(grpData, "iDimensions", 1, 3, PHINT, iDimensions);

    WriteGridCoordinates(grpData, gridIn);
    WriteBCRegion(grpData, gridIn);
    WriteInterfaceInfo(grpData, gridIn);
    WriteStrBCName(grpData, gridIn);
    WriteLnkInfo(grpData, gridIn);
    WriteStrPartitionInfo(grpData, gridIn);
    WriteVolumeCondition(grpData, gridIn);

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteGridCoordinates(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;

    grpData = GetGroup(loc_id, "GridCoordinates");

    RDouble *x = gridIn->GetX();
    RDouble *y = gridIn->GetY();
    RDouble *z = gridIn->GetZ();
    int nTotalNode = gridIn->GetNTotalNode();

    gridIn->RotateAboutAxis();

    CreateAndWriteData(grpData, "CoordinateX", 1, nTotalNode, PHDOUBLE, x);
    CreateAndWriteData(grpData, "CoordinateY", 1, nTotalNode, PHDOUBLE, y);
    CreateAndWriteData(grpData, "CoordinateZ", 1, nTotalNode, PHDOUBLE, z);

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteFaceTopology(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;

    grpData = GetGroup(loc_id, "FaceTopology");

    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalFace = grid->GetNTotalFace();
    int *nodeNumberOfEachFace = grid->GetNodeNumberOfEachFace();

    long long int faceCount = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        faceCount += nodeNumberOfEachFace[iFace];
    }

    int *face2Node = grid->GetFace2Node();
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    int *nodePosi = new int[nTotalFace + 1];
    int nodeSum = 0;
    nodePosi[0] = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        nodeSum += nodeNumberOfEachFace[iFace];
        nodePosi[iFace+1] = nodeSum;
    }

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        if (leftCellOfFace[iFace] < 0)
        {
            std::reverse(face2Node + nodePosi[iFace], face2Node + nodePosi[iFace+1]);
            SWAP(leftCellOfFace[iFace], rightCellOfFace[iFace]);
        }
    }
    delete [] nodePosi;    nodePosi = nullptr;

    CreateAndWriteData(grpData, "nodeNumberOfEachFace", 1, nTotalFace, PHINT, nodeNumberOfEachFace);
    CreateAndWriteData(grpData, "face2Node"           , 1, faceCount , PHINT, face2Node);
    CreateAndWriteData(grpData, "leftCellOfFace"      , 1, nTotalFace, PHINT, leftCellOfFace);
    CreateAndWriteData(grpData, "rightCellOfFace"     , 1, nTotalFace, PHINT, rightCellOfFace);

    int nBoundFace = grid->GetNBoundFace();
    CreateAndWriteData(grpData, "nBoundFace", 1, 1, PHINT, &nBoundFace);

    UnstructBCSet **bcr = grid->GetBCRecord();

    int *tmpBcType = new int[nBoundFace];
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        tmpBcType[iFace] = bcr[iFace]->GetKey();
    }
    CreateAndWriteData(grpData, "BcType", 1, nBoundFace, PHINT, tmpBcType);

    delete [] tmpBcType;    tmpBcType = nullptr;
    H5Gclose(grpData);
}

void IO_HDF5Write::WriteCellTopology(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;
    grpData = GetGroup(loc_id, "CellTopology");

    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int *nodeNumberOfEachCell = grid->GetNodeNumberOfEachCell();
    CreateAndWriteData(grpData, "nodeNumberOfEachCell", 1, nTotalCell, PHINT, nodeNumberOfEachCell);

    int *cell2Node = grid->GetCell2Node();
    long long int nodeCount = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        nodeCount += nodeNumberOfEachCell[iCell];
    }
    CreateAndWriteData(grpData, "cell2Node", 1, nodeCount, PHINT, cell2Node);

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteUnstrBCName(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;

    grpData = GetGroup(loc_id, "BCName");

    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nBoundFace = grid->GetNBoundFace();
    string *bcNameList = new string[nBoundFace];

    uint_long totalSize = 0;
   
    UnstructBCSet **bcr = grid->GetBCRecord();

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        const string &bcName = bcr[iFace]->GetBCName();
        bcNameList[iFace] = bcName;

        totalSize += static_cast< uint_long > (bcName.size());
        totalSize += 1;
    }

    char *bcNameChar = new char[totalSize];
    unsigned int count = 0;
    for (int iFace = 0; iFace < nBoundFace; iFace ++)
    {
        string &bcName = bcNameList[iFace];
        streamsize nameSize = bcName.size();
        for (unsigned int iChar = 0; iChar < nameSize; ++ iChar)
        {
            bcNameChar[count ++] = bcName[iChar];
        }
        bcNameChar[count ++] = '\0';
    }

    CreateAndWriteData(grpData, "BCName", 1, totalSize, PHSTRING, bcNameChar);

    delete [] bcNameList;    bcNameList = nullptr;
    delete [] bcNameChar;    bcNameChar = nullptr;
    H5Gclose(grpData);
}

void IO_HDF5Write::WriteBCRegion(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;
    int iMin, iMax, jMin, jMax, kMin, kMax, bcType;
    int sLr, sNd;

    grpData = GetGroup(loc_id, "BCRegion");

    StructGrid *grid = StructGridCast(gridIn);
    StructBCSet *structBCSet = grid->GetStructBCSet();

    int nBCRegion = structBCSet->GetnBCRegion();
    CreateAndWriteData(grpData, "nBCRegion", 1, 1, PHINT, &nBCRegion);

    int nBCTotal = nBCRegion;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
        bcType = bcregion->GetBCType();
        if (bcType == -1)
        {
            nBCTotal ++;
        }
    }

    int **BCInfor = NewPointer2<int>(nBCTotal, 9);
    int iBCinfo = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        int index = 0;
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
        bcregion->GetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);
        bcType = bcregion->GetBCType();
        sLr = bcregion->GetFaceLeftOrRightIndex();
        sNd = bcregion->GetFaceDirection();

        BCInfor[iBCinfo][index++] = iMin;
        BCInfor[iBCinfo][index++] = iMax;
        BCInfor[iBCinfo][index++] = jMin;
        BCInfor[iBCinfo][index++] = jMax;
        BCInfor[iBCinfo][index++] = kMin;
        BCInfor[iBCinfo][index++] = kMax;
        BCInfor[iBCinfo][index++] = bcType;
        BCInfor[iBCinfo][index++] = sLr;
        BCInfor[iBCinfo][index++] = sNd;
        iBCinfo ++;

        if (bcType == -1)
        {
            index = 0;
            int *t_st = bcregion->GetTargetStart();
            int *t_ed = bcregion->GetTargetEnd();
            int neighborZoneIndex = bcregion->GetTargetRegionBlock();

            BCInfor[iBCinfo][index++] = t_st[0];
            BCInfor[iBCinfo][index++] = t_ed[0];
            BCInfor[iBCinfo][index++] = t_st[1];
            BCInfor[iBCinfo][index++] = t_ed[1];
            BCInfor[iBCinfo][index++] = t_st[2];
            BCInfor[iBCinfo][index++] = t_ed[2];
            BCInfor[iBCinfo][index++] = neighborZoneIndex;
            BCInfor[iBCinfo][index++] = 0;
            BCInfor[iBCinfo][index++] = 0;
            iBCinfo ++;
        }
    }

    CreateAndWriteData(grpData, "BCInfor", nBCTotal, 9, PHINT, BCInfor[0]);
    DelPointer2(BCInfor);

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteStrBCName(hid_t loc_id, Grid *gridIn)
{
     hid_t grpData;

    grpData = GetGroup(loc_id, "BCName");

    string boundaryName[100] = {"NO_BOUNDARY_CONDITION", "EXTRAPOLATION", "SOLID_SURFACE", "SYMMETRY", 
                                "FARFIELD", "INFLOW", "OUTFLOW", "POLE", "Wall_8", "Wall_9", "Wall_10", "Wall_11", 
                                "Wall_12", "Wall_13", "Wall_14", "Wall_15", "Wall_16", "Wall_17", "Wall_18", "Wall_19",
                                "Wall_20", "Wall_21", "Wall_22", "Wall_23", "Wall_24", "Wall_25", "Wall_26", "Wall_27", 
                                "Wall_28", "Wall_29", "Wall_30", "Wall_31", "Wall_32", "Wall_33", "Wall_34", "Wall_35",
                                "Wall_36", "Wall_37", "Wall_38", "Wall_39", "Wall_40", "Wall_41", "Wall_42", "Wall_43",
                                "Wall_44", "Wall_45", "Wall_46", "Wall_47", "Wall_48", "Wall_49", "Wall_50"};
    boundaryName[52] = "PRESSURE_INLET";
    boundaryName[62] = "PRESSURE_OUTLET";
    boundaryName[71] = "POLE1";
    boundaryName[72] = "POLE2";
    boundaryName[73] = "POLE3";

    StructGrid *grid = StructGridCast(gridIn);
    StructBCSet *structBCSet = grid->GetStructBCSet();

    int nBCRegion = static_cast<int> (structBCSet->GetnBCRegion());
    string *bcNameList = new string[nBCRegion];
    uint_long totalSize = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion ++)
    {
        StructBC *bcRegion = structBCSet->GetBCRegion(iBCRegion);

        string bcName = bcRegion->GetBCName();
        int bctype = bcRegion->GetBCType();
        if (bcName == "")
        {
            if (bctype == -1)
            {
                bcName = "Connection";
            }
            else
            {
                bcName = boundaryName[bctype];
                if ((bctype > 7) && (bctype < 51))
                {
                    bctype = 2;
                    bcRegion->SetBCType(bctype);
                }
            }
            bcRegion->SetBCName(bcName);
        }
        bcNameList[iBCRegion] = bcName;

        totalSize += static_cast< uint_long > (bcNameList[iBCRegion].size());
        totalSize += 1;
    }

    char *bcNameChar = new char[totalSize];
    unsigned int count = 0;
    for (int boco = 0; boco < nBCRegion; ++ boco)
    {
        string &bcName = bcNameList[boco];
        streamsize nameSize = bcName.size();
        for (unsigned int iChar = 0; iChar < nameSize; ++ iChar)
        {
            bcNameChar[count ++] = bcName[iChar];
        }
        bcNameChar[count ++] = '\0';
    }

    CreateAndWriteData(grpData, "BCName", 1, totalSize, PHSTRING, bcNameChar);

    delete [] bcNameChar;    bcNameChar = nullptr;
    delete [] bcNameList;    bcNameList = nullptr;
    H5Gclose(grpData);
}

void IO_HDF5Write::WriteVolumeCondition(hid_t loc_id, Grid *gridIn)
{
    SimpleVC *volumeCondition = gridIn->GetVolumeConditionIn();
    if (!volumeCondition)
    {
        return;
    }

    hid_t grpData;
    grpData = GetGroup(loc_id, "VolumeCondition");

    string vcName = volumeCondition->GetVCName();
    uint_long nameSize = vcName.size() + 1;

    char *vcNameChar = new char[nameSize];
    for (unsigned int iChar = 0; iChar < nameSize - 1; ++ iChar)
    {
        vcNameChar[iChar] = vcName[iChar];
    }
    vcNameChar[nameSize - 1] = '\0';

    CreateAndWriteData(grpData, "VCName", 1, nameSize, PHSTRING, vcNameChar);
    delete [] vcNameChar;    vcNameChar = nullptr;

    int vcType = volumeCondition->GetVCType();
    CreateAndWriteData(grpData, "VCType", 1, 1, PHINT, &vcType);

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteInterfaceInfo(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;

    grpData = GetGroup(loc_id, "InterfaceInfo");

    InterfaceInfo *interfaceInfo = gridIn->GetInterfaceInfo();

    int nIFace = 0;
    if (interfaceInfo)
    {
        nIFace = interfaceInfo->GetNIFace();
    }
    CreateAndWriteData(grpData, "nIFace", 1, 1, PHINT, &nIFace);

    if (nIFace)
    {
        int *interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
        int *interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();
        int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

        CreateAndWriteData(grpData, "interFace2ZoneID"      , 1, nIFace, PHINT, interFace2ZoneID);
        CreateAndWriteData(grpData, "interFace2InterFaceID" , 1, nIFace, PHINT, interFace2InterFaceID);
        CreateAndWriteData(grpData, "interFace2BoundaryFace", 1, nIFace, PHINT, interFace2BoundaryFace);
    }

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteLnkInfo(hid_t loc_id, Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    StructBCSet *structBCSet = grid->GetStructBCSet();

    int nBCLnkTotal = 0;
    int nBCRegion = structBCSet->GetnBCRegion();
    for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion ++)
    {
        StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);

        int bctype = structBC->GetBCType();
        string bcName = structBC->GetBCName();
        if (bctype == PHENGLEI::INTERFACE && !(bcName == "Periodic_up" || bcName == "Periodic_down"))
        {
            nBCLnkTotal ++;
        }
    }

    if (!nBCLnkTotal)
    {
        return;
    }

    hid_t grpData;
    grpData = GetGroup(loc_id, "LnkInfo");

    int **LnkInfor = NewPointer2<int>(nBCLnkTotal * 2, 7);
    int iLnkinfo = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion ++)
    {
        StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);

        int bctype = structBC->GetBCType();
        string bcName = structBC->GetBCName();
        if (!IsInterface(bctype) || bcName == "Periodic_up" || bcName == "Periodic_down") continue;

        int index = 0;
        int icount = 0;
        int nbt = structBC->GetTargetRegionBlock() + 1;
        int *lnkInfo = structBC->GetInkInfo();

        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = bctype;

        iLnkinfo ++;
        index = 0;

        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = nbt;

        iLnkinfo ++;
    }

    CreateAndWriteData(grpData, "LnkInfor", nBCLnkTotal * 2, 7, PHINT, LnkInfor[0]);
    DelPointer2(LnkInfor);

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteStrPartitionInfo(hid_t loc_id, Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ordinaryGridIndex = grid->GetOrdinaryGridIndex();
    if (ordinaryGridIndex == -1)
    {
        return;
    }

    hid_t grpData;
    grpData = GetGroup(loc_id, "PartitionInfo");

    int *ordinaryDimStartIndex = grid->GetOrdinaryDimStartIndex();
    int *ordinaryDimEndIndex   = grid->GetOrdinaryDimEndIndex();

    CreateAndWriteData(grpData, "OrdinaryGrid"    , 1, 1, PHINT, &ordinaryGridIndex);
    CreateAndWriteData(grpData, "OrdinaryDimStart", 1, 3, PHINT, ordinaryDimStartIndex);
    CreateAndWriteData(grpData, "OrdinaryDimEnd"  , 1, 3, PHINT, ordinaryDimEndIndex  );

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteUnstrPartitionInfo(hid_t loc_id, Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int ordinaryGridIndex = grid->GetOrdinaryGridIndex();
    if (ordinaryGridIndex == -1)
    {
        return;
    }

    hid_t grpData;
    grpData = GetGroup(loc_id, "PartitionInfo");

    int nTotalNode = grid->GetNTotalNode();
    int nTotalFace = grid->GetNTotalFace();
    int nTotalCell = grid->GetNTotalCell();
    int *ordinaryNodeIndex = grid->GetOrdinaryNodeIndex();
    int *ordinaryFaceIndex = grid->GetOrdinaryFaceIndex();
    int *ordinaryCellIndex = grid->GetOrdinaryCellIndex();

    CreateAndWriteData(grpData, "OrdinaryGrid"     , 1, 1         , PHINT, &ordinaryGridIndex);
    CreateAndWriteData(grpData, "OrdinaryNodeIndex", 1, nTotalNode, PHINT, ordinaryNodeIndex);
    CreateAndWriteData(grpData, "OrdinaryFaceIndex", 1, nTotalFace, PHINT, ordinaryFaceIndex);
    CreateAndWriteData(grpData, "OrdinaryCellIndex", 1, nTotalCell, PHINT, ordinaryCellIndex);

    H5Gclose(grpData);
}

IO_HDF5Read::IO_HDF5Read(PHVectorString1D &gridNameListIn, Region *region_in)
{
    this->gridNameList = gridNameListIn;
    this->region = region_in;

    this->numberOfGridGroups = gridNameList.size();
    this->grids.resize(0);
    this->myProcID = PHMPI::GetCurrentProcessorID();
    this->serverProcID = 0;
    this->zoneStart = 0;
    this->presentGridGroup = 0;
    this->fileID = -1;
    this->isFileProcessor = false;
    this->taskType = GetTaskCode();
    this->isP2PMode = false;

    blockProcDump = 0;
    numberOfZones = 0;
    blockProc     = 0;
    blockType     = 0;
    blockIdx      = 0;
    blockProcGrid = 0;
    blockFileProc = 0;
    nZonesList.resize(numberOfGridGroups);
    localGridNameList.resize(numberOfGridGroups);
    version = 0.0;
    globalFileNumber = 0;
}

IO_HDF5Read::~IO_HDF5Read()
{
    delete [] blockProcDump;    blockProcDump = nullptr;
    delete [] blockProc;    blockProc = nullptr;
    delete [] blockType;    blockType = nullptr;
    delete [] blockIdx;    blockIdx = nullptr;
    delete [] blockProcGrid;    blockProcGrid = nullptr;
    delete [] blockFileProc;    blockFileProc = nullptr;
}

void IO_HDF5Read::CheckNumberofFiles(int iGrid)
{
    string gridFileName = gridNameList[iGrid];

    string gridfileProc = gridFileName;
    gridfileProc = AddSymbolToFileName(gridfileProc, "_", myProcID);
    bool fileExist = FileExist(gridfileProc);

    int nFileProcessors = 0;
    if (!fileExist)
    {
        if (myProcID == serverProcID)
        {
            //! The server processor must be file processor.
            TK_Exit::FileOpenErrorExit(gridfileProc);
        }
    }
    else
    {
        nFileProcessors = 1;
    }

    PH_AllReduce(&nFileProcessors, &globalFileNumber, 1, MPI_SUM);

    PHMPI::SetNumberofFiles(globalFileNumber);
    WriteLogFile("Global grid files number: ", globalFileNumber);

    nZonesList[iGrid].resize(globalFileNumber);

    localGridNameList[iGrid].resize(0);
    for (int iFile = 0; iFile < globalFileNumber; ++ iFile)
    {
        string gridfile = gridFileName;
        gridfile = AddSymbolToFileName(gridfile, "_", iFile);
        localGridNameList[iGrid].push_back(gridfile);
    }

    if (globalFileNumber > 1)
    {
        isP2PMode = true;
    }
}

void IO_HDF5Read::InitGridFileNameListAndFileProcessor()
{
    for (int iGrid = 0; iGrid < numberOfGridGroups; ++ iGrid)
    {
        CheckNumberofFiles(iGrid);
    }

    if (PHMPI::GetNumberOfProcessor() % PHMPI::GetNumberofFiles() != 0)
    {
        ostringstream oss;
        oss << "   ---- numberOfProcessor = " << PHMPI::GetNumberOfProcessor() << ", numberofFiles = " << PHMPI::GetNumberofFiles() << " ----\n" << "Abnormally Exit Program!!!\n";
        TK_Exit::ExceptionExit(oss.str());
    }

    int fileProcessorInterval = PHMPI::GetIntervalofProcessorFiles();
    isFileProcessor = (myProcID % fileProcessorInterval == 0);
    //isFileProcessor = (myProcID % fileProcessorInterval == 0 && myProcID < globalFileNumber * fileProcessorInterval);

    this->fileID = PHMPI::GetFileIndexofCurrentProcessor();
}

hid_t IO_HDF5Read::OpenH5File(const string &filename)
{
    hid_t file;

    file = OpenHDF5File(AddSymbolToFileName(filename, "_", this->fileID));
    return file;
}

bool IO_HDF5Read::Run()
{
    bool HDF5FileRead = false;

    InitGridFileNameListAndFileProcessor();
    ReadNumberOfZones();
    ReadTotalInfo();
    SetMixGridInfo();
    RedistributeZonesIfNeed();
    SetGlobalZoneLayout();
    ClassifyGridSystem();
    ReadEachGrid();
    CreateZones();
    ReadWallDist();

    HDF5FileRead = true;
    return HDF5FileRead;
}

void IO_HDF5Read::ReadWallDist()
{
    string WallDistFileName = ChangeExtensionOfFileName(gridNameList.front(), "wdt");

    bool fileExist = true;
    if (isFileProcessor)
    {
        fileExist = FileExist(AddSymbolToFileName(WallDistFileName, "_", this->fileID));
    }

    bool fileCheck = true;
    PH_AllReduce(&fileExist, &fileCheck, 1, MPI_MIN);
    if (!fileCheck)
    {
        return;
    }

    hid_t file;
    int zoneCount = 0;

    zoneStart = 0;
    for (int iGrid = 0; iGrid < numberOfGridGroups; ++ iGrid)
    {
        int nZonesOfGridGroup = 0;
        for (int iFile = 0; iFile < localGridNameList[iGrid].size(); ++ iFile)
        {
            file = OpenHDF5File(ChangeExtensionOfFileName(localGridNameList[0][iFile], "wdt"));
            if (file < 0)
            {
                return;
            }

            for (int iZone = 0; iZone < nZonesList[iGrid][iFile]; ++ iZone)
            {
                int gridIndex = blockIdx[zoneCount ++];
                if (myProcID != blockProc[gridIndex + zoneStart])
                {
                    continue;
                }

                ostringstream oss;
                oss << "Grid-" << gridIndex + zoneStart;
                string dataName = oss.str();

                Grid *grid = GetGrid(gridIndex + zoneStart, 0);
                grid->InitWallDist();

                int nTotalCell = grid->GetNTotalCell();
                double *walldist = new double[nTotalCell];
                ReadData(file, walldist, dataName);

                RDouble *wdt = new RDouble[nTotalCell];
                for (int iCell = 0; iCell < nTotalCell; ++ iCell)
                {
                    wdt[iCell] = static_cast<RDouble>(walldist[iCell]);
                }
                grid->ReadWalldist(wdt);

                string nearestwallfacenormalxname = dataName + "_nwfnx";
                string nearestwallfacenormalyname = dataName + "_nwfny";
                string nearestwallfacenormalzname = dataName + "_nwfnz";
                if (CheckDataExist(file, nearestwallfacenormalxname))
                {
                    double *nwfnx = new double[nTotalCell];
                    ReadData(file, nwfnx, nearestwallfacenormalxname);

                    RDouble *nearestwallfacenormalx = new RDouble[nTotalCell];
                    for (int iCell = 0; iCell < nTotalCell; ++iCell)
                    {
                        nearestwallfacenormalx[iCell] = static_cast<RDouble>(nwfnx[iCell]);
                        //cout << nearestwallfacenormalx[iCell] << endl;
                    }

                    double *nwfny = new double[nTotalCell];
                    ReadData(file, nwfny, nearestwallfacenormalyname);

                    RDouble *nearestwallfacenormaly = new RDouble[nTotalCell];
                    for (int iCell = 0; iCell < nTotalCell; ++iCell)
                    {
                        nearestwallfacenormaly[iCell] = static_cast<RDouble>(nwfny[iCell]);
                    }

                    double *nwfnz = new double[nTotalCell];
                    ReadData(file, nwfnz, nearestwallfacenormalzname);

                    RDouble *nearestwallfacenormalz = new RDouble[nTotalCell];
                    for (int iCell = 0; iCell < nTotalCell; ++iCell)
                    {
                        nearestwallfacenormalz[iCell] = static_cast<RDouble>(nwfnz[iCell]);
                    }

                    grid->ReadNearestwallfacenormalx(nearestwallfacenormalx);
                    grid->ReadNearestwallfacenormaly(nearestwallfacenormaly);
                    grid->ReadNearestwallfacenormalz(nearestwallfacenormalz);

                    delete [] nearestwallfacenormalx;    nearestwallfacenormalx = nullptr;
                    delete [] nearestwallfacenormaly;    nearestwallfacenormaly = nullptr;
                    delete [] nearestwallfacenormalz;    nearestwallfacenormalz = nullptr;

                    delete [] nwfnx;    nwfnx = nullptr;
                    delete [] nwfny;    nwfny = nullptr;
                    delete [] nwfnz;    nwfnz = nullptr;
                }

                delete [] walldist;    walldist = nullptr;
                delete [] wdt;    wdt = nullptr;
            }

            nZonesOfGridGroup += nZonesList[iGrid][iFile];
            H5Fclose(file);
        }

        zoneStart += nZonesOfGridGroup;
    }
}

void IO_HDF5Read::ReadVersionInfo(hid_t loc_id)
{
    hid_t grpData;

    grpData = OpenGroup(loc_id, "Version");

    ReadData(grpData, &this->version, "VersionData");

    H5Gclose(grpData);
}

void IO_HDF5Read::ReadNumberOfZones()
{
    hid_t file;

    for (int iGrid = 0; iGrid < numberOfGridGroups; ++ iGrid)
    {
        for (int iFile = 0; iFile < globalFileNumber; ++ iFile)
        {
            file = OpenHDF5File(localGridNameList[iGrid][iFile]);
            if (file < 0)
            {
                TK_Exit::ExceptionExit("Grid file can not open !");
            }

            hid_t grpData;
            grpData = OpenGroup(file, "Information");

            int nZones = 0;
            ReadData(grpData, &nZones, "nBlocks");
            numberOfZones += nZones;
            nZonesList[iGrid][iFile] = nZones;

            H5Gclose(grpData);
            H5Fclose(file);
        }
    }

    PHMPI::SetNumberOfGlobalZones   (numberOfZones);
    PHMPI::CreateZoneProcessorIDDump(numberOfZones);
    PHMPI::CreateZoneProcessorID    (numberOfZones);
    PHMPI::CreateZoneGridID         (numberOfZones);
    PHMPI::CreateZoneGridType       (numberOfZones);
    PHMPI::CreateZoneFileID         (numberOfZones);

    int numberOfProcessor = PHMPI::GetNumberOfProcessor();
    if (numberOfZones < numberOfProcessor)
    {
        ostringstream oss;
        oss << "   ---- numberOfProcessor = " << numberOfProcessor << ", numberOfTotalZones = " << numberOfZones << " ----\n" << "Abnormally Exit Program!!!\n";
        TK_Exit::ExceptionExit(oss.str());
    }
}

void IO_HDF5Read::ReadNumberOfZonesOld(hid_t loc_id)
{
    hid_t grpData;
    grpData = OpenGroup(loc_id, "Information");

    int *nZonesListTemp = new int[globalFileNumber];
    SetField(nZonesListTemp, 0, globalFileNumber);

    int nZones = 0;
    if (isFileProcessor)
    {
        ReadData(grpData, &nZones, "nBlocks");
        nZonesListTemp[fileID] = nZones;
        //MPI_Allgather(&nZones, 1, MPI_INT, nZonesList, 1, MPI_INT, MPI_COMM_WORLD);
    }

    PH_AllReduce(&nZones, &numberOfZones, 1, MPI_SUM);
    //PH_AllReduce(nZonesListTemp, nZonesList, globalFileNumber, MPI_SUM);
    delete [] nZonesListTemp;    nZonesListTemp = nullptr;

    PHMPI::SetNumberOfGlobalZones   (numberOfZones);
    PHMPI::CreateZoneProcessorIDDump(numberOfZones);
    PHMPI::CreateZoneProcessorID    (numberOfZones);
    PHMPI::CreateZoneGridID         (numberOfZones);
    PHMPI::CreateZoneGridType       (numberOfZones);
    PHMPI::CreateZoneFileID         (numberOfZones);

    int numberOfProcessor = PHMPI::GetNumberOfProcessor();
    if (numberOfZones < numberOfProcessor)
    {
        ostringstream oss;
        oss << "   ---- numberOfProcessor = " << numberOfProcessor << ", numberOfTotalZones = " << numberOfZones << " ----\n" << "Abnormally Exit Program!!!\n";
        TK_Exit::ExceptionExit(oss.str());
    }

    H5Gclose(grpData);
}

void IO_HDF5Read::SetMixGridInfo()
{
    int nsimutask = GlobalDataBase::GetIntParaFromDB("nsimutask");
    if (!(nsimutask == CREATE_GRID || nsimutask == GRIDCHECK_SOLVER))
    {
        return;
    }

    if (nsimutask == CREATE_GRID)
    {
        int gridObj = GlobalDataBase::GetIntParaFromDB("gridobj");
        if (gridObj != PHSPACE::GRID_CONVERSION && gridObj != PHSPACE::GRID_STRTOUNS)
        {
            return;
        }

        int gridType = GlobalDataBase::GetIntParaFromDB("gridtype");
        if (gridType != PHSPACE::MIXGRID && gridObj != PHSPACE::GRID_STRTOUNS)
        {
            return;
        }
    }

    CreatetZoneProcessorID_INP(numberOfZones);
    int *global_block_proc_inp = GetZoneProcessorID_INP();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        global_block_proc_inp[iZone] = blockProc[iZone];
    }
}

void IO_HDF5Read::ReadTotalInfo()
{
    blockProc     = new int[numberOfZones];
    blockType     = new int[numberOfZones];
    blockIdx      = new int[numberOfZones];
    blockFileProc = new int[numberOfZones];

    SetField(blockProc    , 0, numberOfZones);
    SetField(blockType    , 0, numberOfZones);
    SetField(blockIdx     , 0, numberOfZones);
    SetField(blockFileProc, 0, numberOfZones);

    int processorNumInEachFile = PHMPI::GetNumberOfProcessor() / PHMPI::GetNumberofFiles();

    hid_t file;
    int count = 0;

    for (int iGrid = 0; iGrid < numberOfGridGroups; ++ iGrid)
    {
        for (int iFile = 0; iFile < globalFileNumber; ++ iFile)
        {
            file = OpenHDF5File(localGridNameList[iGrid][iFile]);

            hid_t grpData;
            grpData = OpenGroup(file, "Information");

            int *blockProc_local = ReadIntOneRow(grpData, "block_proc");
            int *blockIdx_local  = ReadIntOneRow(grpData, "block_idx");
            int *blockType_local = ReadIntOneRow(grpData, "block_type");

            for (int iZone = 0; iZone < nZonesList[iGrid][iFile]; ++ iZone)
            {
                blockProc[count + iZone] = blockProc_local[iZone];
                blockType[count + iZone] = blockType_local[iZone];
                blockIdx [count + iZone] = blockIdx_local [iZone];

                blockFileProc[count + iZone] = processorNumInEachFile * iFile;
            }
            delete [] blockProc_local;    blockProc_local = nullptr;
            delete [] blockIdx_local ;    blockIdx_local = nullptr;
            delete [] blockType_local;    blockType_local = nullptr;

            count += nZonesList[iGrid][iFile];
            H5Gclose(grpData);
            H5Fclose(file);
        }
    }

    blockProcDump = new int[numberOfZones];
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        blockProcDump[iZone] = blockProc[iZone];
    }
}

void IO_HDF5Read::ReadTotalInfoOld(hid_t loc_id)
{
    hid_t grpData;
    grpData = OpenGroup(loc_id, "Information");

    int *blockProcTemp     = new int[numberOfZones];
    int *blockTypeTemp     = new int[numberOfZones];
    int *blockIdxTemp      = new int[numberOfZones];
    int *blockFileProcTemp = new int[numberOfZones];

    blockProc     = new int[numberOfZones];
    blockType     = new int[numberOfZones];
    blockIdx      = new int[numberOfZones];
    blockFileProc = new int[numberOfZones];

    SetField(blockProcTemp    , 0, numberOfZones);
    SetField(blockTypeTemp    , 0, numberOfZones);
    SetField(blockIdxTemp     , 0, numberOfZones);
    SetField(blockFileProcTemp, 0, numberOfZones);

    if (isFileProcessor)
    {
        int *blockProc_local = ReadIntOneRow(grpData, "block_proc");
        int *blockIdx_local  = ReadIntOneRow(grpData, "block_idx");
        int *blockType_local = ReadIntOneRow(grpData, "block_type");

        int count = 0;
        for (int ifile = 0; ifile < fileID; ++ ifile)
        {
            count += nZonesList[0][ifile];
        }

        int myid = PHMPI::GetCurrentProcessorID();
        for (int iZone = 0; iZone < nZonesList[0][fileID]; ++ iZone)
        {
            blockProcTemp[count + iZone] = blockProc_local[iZone];
            blockTypeTemp[count + iZone] = blockType_local[iZone];
            blockIdxTemp [count + iZone] = blockIdx_local [iZone];

            blockFileProcTemp[count + iZone] = myid;
        }

        delete [] blockProc_local;    blockProc_local = nullptr;
        delete [] blockIdx_local ;    blockIdx_local = nullptr;
        delete [] blockType_local;    blockType_local = nullptr;
    }

    PH_AllReduce(blockProcTemp    , blockProc    , numberOfZones, MPI_SUM);
    PH_AllReduce(blockIdxTemp     , blockIdx     , numberOfZones, MPI_SUM);
    PH_AllReduce(blockTypeTemp    , blockType    , numberOfZones, MPI_SUM);
    PH_AllReduce(blockFileProcTemp, blockFileProc, numberOfZones, MPI_SUM);

    delete [] blockProcTemp;    blockProcTemp = nullptr;
    delete [] blockIdxTemp ;    blockIdxTemp = nullptr;
    delete [] blockTypeTemp;    blockTypeTemp = nullptr;
    delete [] blockFileProcTemp;    blockFileProcTemp = nullptr;

    blockProcDump = new int[numberOfZones];
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        blockProcDump[iZone] = blockProc[iZone];
    }

    H5Gclose(grpData);
}

void IO_HDF5Read::ReadEachGrid()
{
    int *blockProc = GetBlockProc();
    int *blockType = GetBlockType();

    hid_t file;
    zoneStart = 0;
    for (int iGrid = 0; iGrid < numberOfGridGroups; ++ iGrid)
    {
        presentGridGroup = iGrid;
        int nZonesOfGridGroup = 0;
        for (int iFile = 0; iFile < globalFileNumber; ++ iFile)
        {
            file = OpenHDF5File(localGridNameList[iGrid][iFile]);

            hid_t grpData;
            grpData = OpenGroup(file, "Information");
            int *blockIdx_local = ReadIntOneRow(grpData, "block_idx");
            H5Gclose(grpData);

            for (int iZone = 0; iZone < nZonesList[iGrid][iFile]; ++ iZone)
            {
                int gridIndex = blockIdx_local[iZone];
                if (myProcID != blockProc[gridIndex + zoneStart])
                {
                    AddGrid(0);
                    continue;
                }

                int gridType = blockType[gridIndex + zoneStart];
                if (gridType == PHSPACE::UNSTRUCTGRID)
                {
                    ReadUnstructGrid(file, gridIndex);
                }
                else
                {
                    ReadStructGrid(file, gridIndex);
                }
            }

            nZonesOfGridGroup += nZonesList[iGrid][iFile];
            delete [] blockIdx_local;    blockIdx_local = nullptr;
            H5Fclose(file);
        }
        zoneStart += nZonesOfGridGroup;
    }
}

void IO_HDF5Read::ReadEachGridOld(hid_t loc_id)
{
    int *blockProc = GetBlockProc();
    int *blockType = GetBlockType();

    int nZones = GetNZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        if (myProcID != blockProc[iZone])
        {
            AddGrid(0);
            continue;
        }

        int gridType = blockType[iZone];
        if (gridType == PHSPACE::UNSTRUCTGRID)
        {
            ReadUnstructGrid(loc_id, iZone);
        }
        else
        {
            ReadStructGrid(loc_id, iZone);
        }
    }
}

void IO_HDF5Read::ReadGridCoordinates(hid_t loc_id, Grid *gridIn)
{
    PrintToWindow("    Reading coordinates ...\n");
    hid_t grpData;

    grpData = OpenGroup(loc_id, "GridCoordinates");

    RDouble *CoordinateX = ReadDoubleOneRow(grpData, "CoordinateX");
    RDouble *CoordinateY = ReadDoubleOneRow(grpData, "CoordinateY");
    RDouble *CoordinateZ = ReadDoubleOneRow(grpData, "CoordinateZ");

    gridIn->SetX(CoordinateX);
    gridIn->SetY(CoordinateY);
    gridIn->SetZ(CoordinateZ);

    gridIn->RotateAboutAxis();

    gridIn->ComputeMinMaxBox();

    H5Gclose(grpData);
}

void IO_HDF5Read::ReadInterfaceInfo(hid_t loc_id, Grid *gridIn)
{
    PrintToWindow("    Reading interface ...\n");
    hid_t grpData;

    grpData = OpenGroup(loc_id, "InterfaceInfo");

    InterfaceInfo *interfaceInfo = 0;
    gridIn->SetInterfaceInfo(interfaceInfo);

    int nIFace = 0;
    ReadData(grpData, &nIFace, "nIFace");

    if (nIFace > 0)
    {
        interfaceInfo = new InterfaceInfo(nIFace, gridIn);
        gridIn->SetInterfaceInfo(interfaceInfo);

        int *interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
        int *interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();
        int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

        int *interFace2ZoneIDTemp       = ReadIntOneRow(grpData, "interFace2ZoneID");
        int *interFace2InterFaceIDTemp  = ReadIntOneRow(grpData, "interFace2InterFaceID");
        int *interFace2BoundaryFaceTemp = ReadIntOneRow(grpData, "interFace2BoundaryFace");

        for (int iFace = 0; iFace < nIFace; ++ iFace)
        {
            interFace2ZoneID[iFace]       = interFace2ZoneIDTemp[iFace];
            interFace2InterFaceID[iFace]  = interFace2InterFaceIDTemp[iFace];
            interFace2BoundaryFace[iFace] = interFace2BoundaryFaceTemp[iFace];
        }

        delete [] interFace2ZoneIDTemp;    interFace2ZoneIDTemp = nullptr;
        delete [] interFace2InterFaceIDTemp;    interFace2InterFaceIDTemp = nullptr;
        delete [] interFace2BoundaryFaceTemp;    interFace2BoundaryFaceTemp = nullptr;
    }

    H5Gclose(grpData);
}

void IO_HDF5Read::ReadVolumeCondition(hid_t loc_id, Grid *gridIn)
{
    bool dataExist = CheckDataExist(loc_id, "VolumeCondition");
    if (!dataExist)
    {
        return;
    }

    hid_t grpData;
    grpData = OpenGroup(loc_id, "VolumeCondition");
    char *VCName = ReadCharData(grpData, "VCName");

    SimpleVC *volumeCondition = new SimpleVC();

    int vcType = -1;
    ReadData(grpData, &vcType, "VCType");
    volumeCondition->SetVCType(vcType);

    volumeCondition->SetVCName(VCName);
    gridIn->SetVolumeCondition(volumeCondition);

    delete [] VCName;    VCName = nullptr;
    H5Gclose(grpData);
}

void IO_HDF5Read::ReadStructGrid(hid_t loc_id, int iZone)
{
    PrintToWindow("Reading Structured Grid Of Zone ", iZone, "...\n");

    hid_t  grpData;
    string grpName;

    GridID *index = new GridID(iZone + zoneStart);
    Grid *grid = new StructGrid();

    grid->InitGrid(index, 0, PHSPACE::GetDim(), STRUCTGRID);
    grid->SetIBlock(presentGridGroup);

    ostringstream oss;
    oss << "Grid-" << iZone;
    grpName = oss.str();
    grpData = OpenGroup(loc_id, grpName);

    StructGrid *strGrid = StructGridCast(grid);
    int *iDimensions = ReadIntOneRow(grpData, "iDimensions");
    strGrid->SetNI(iDimensions[0]);
    strGrid->SetNJ(iDimensions[1]);
    strGrid->SetNK(iDimensions[2]);
    strGrid->SetBasicDimension();

    ReadGridCoordinates(grpData, grid);
    strGrid->SetArrayLayout();

    ReadStrBCName(grpData, grid);
    ReadBCRegion(grpData, grid);
    ReadStrPartitionInfo(grpData, grid);
    ReadInterfaceInfo(grpData, grid);
    ReadVolumeCondition(grpData, grid);

    grid->ReviseNeighborZoneIndex(zoneStart);
    AddGrid(grid);
    delete [] iDimensions;    iDimensions = nullptr;
    H5Gclose(grpData);
}

void IO_HDF5Read::ReadLnk(hid_t loc_id, int iZone)
{
    int nZones = PHMPI::GetNumberofGlobalZones();
    Grid **OrdinaryGrid = new Grid *[nZones];

    hid_t  grpData;
    string grpName;

    ostringstream oss;
    oss << "Grid-" << iZone;
    grpName = oss.str();
    grpData = OpenGroup(loc_id, grpName);

    OrdinaryGrid[iZone] = GetGrid(iZone, 0);

    ReadLnkInfo(grpData, OrdinaryGrid[iZone]);
    H5Gclose(grpData);
}

void IO_HDF5Read::ReadBCRegion(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;

    grpData = OpenGroup(loc_id, "BCRegion");

    StructGrid *grid = StructGridCast(gridIn);
    StructBCSet *structBCSet = grid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();

    int **BCInfor = ReadIntTwoRow(grpData, "BCInfor");

    int iMin, iMax, jMin, jMax, kMin, kMax, bcType;
    int sLr, sNd;
    int iBCinfo = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        int index = 0;
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);

        iMin   = BCInfor[iBCinfo][index++];
        iMax   = BCInfor[iBCinfo][index++];
        jMin   = BCInfor[iBCinfo][index++];
        jMax   = BCInfor[iBCinfo][index++];
        kMin   = BCInfor[iBCinfo][index++];
        kMax   = BCInfor[iBCinfo][index++];
        bcType = BCInfor[iBCinfo][index++];
        sLr    = BCInfor[iBCinfo][index++];
        sNd    = BCInfor[iBCinfo][index++];
        iBCinfo ++;

        //! Restore periodic boundary bcType and correct read information only for PARTITION_GRID task.
        if (GetTaskCode() == PARTITION_GRID && PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType") != NO_PERIODICITY)
        {
            if (bcregion->GetBCName() == "Periodic_up"|| bcregion->GetBCName() == "Periodic_down")
            {
                bcType = PHENGLEI::USER_DEFINED;
                iBCinfo ++;
            }
        }

        bcregion->SetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);
        bcregion->SetBCType(bcType);
        bcregion->SetFaceLeftOrRightIndex(sLr);
        bcregion->SetFaceDirection(sNd);
        bcregion->InitFaceDirectionIndex();

        if (bcType == -1)
        {
            index = 0;

            iMin = BCInfor[iBCinfo][index++];
            iMax = BCInfor[iBCinfo][index++];
            jMin = BCInfor[iBCinfo][index++];
            jMax = BCInfor[iBCinfo][index++];
            kMin = BCInfor[iBCinfo][index++];
            kMax = BCInfor[iBCinfo][index++];
            bcregion->SetTargetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);

            int neighborZoneIndex = BCInfor[iBCinfo][index++];
            bcregion->SetTargetRegionBlock(neighborZoneIndex);

            iBCinfo ++;
        }
    }
    grid->SetBCFaceInfo();

    DelPointer2(BCInfor);
    H5Gclose(grpData);
}

void IO_HDF5Read::ReadLnkInfo(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;

    grpData = OpenGroup(loc_id, "LnkInfo");

    StructGrid *grid = StructGridCast(gridIn);
    StructBCSet *structBCSet = grid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();

    int **LnkInfor = ReadIntTwoRow(grpData, "LnkInfor");

    int iMin, iMax, jMin, jMax, kMin, kMax, nbt;
    int iLnkinfo = 0;
    int ndim = PHSPACE::GlobalDataBase::GetIntParaFromDB("ndim");

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        int index = 0;
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);

        int bctype = bcregion->GetBCType();
        if ((!IsInterface(bctype))) continue;

        iMin = LnkInfor[iLnkinfo][index++];
        iMax = LnkInfor[iLnkinfo][index++];
        jMin = LnkInfor[iLnkinfo][index++];
        jMax = LnkInfor[iLnkinfo][index++];

        if (ndim == THREE_D)
        {
            kMin = LnkInfor[iLnkinfo][index++];
            kMax = LnkInfor[iLnkinfo][index++];
        }
        else
        {
            kMin = 1;
            kMax = 1;
        }
        iLnkinfo++;

        bcregion->SetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);

        index = 0;

        iMin = LnkInfor[iLnkinfo][index++];
        iMax = LnkInfor[iLnkinfo][index++];
        jMin = LnkInfor[iLnkinfo][index++];
        jMax = LnkInfor[iLnkinfo][index++];

        if (ndim == THREE_D)
        {
            kMin = LnkInfor[iLnkinfo][index++];
            kMax = LnkInfor[iLnkinfo][index++];
        }
        else
        {
            kMin = 1;
            kMax = 1;
        }

        nbt = LnkInfor[iLnkinfo][index++];
        bcregion->SetTargetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);
        bcregion->SetTargetRegionBlock(nbt - 1);
        if (GetDim() == THREE_D)
        {
            bcregion->ComputeRelativeParameters();
        }

        iLnkinfo++;
    }

    DelPointer2(LnkInfor);
    H5Gclose(grpData);
}

void IO_HDF5Read::ReadStrBCName(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;
    
    grpData = OpenGroup(loc_id, "BCRegion");
    StructGrid *grid = StructGridCast(gridIn);

    int nBCRegion = 0;
    ReadData(grpData, &nBCRegion, "nBCRegion");

    grid->CreateCompositeBCRegion(nBCRegion);
    H5Gclose(grpData);

    grpData = OpenGroup(loc_id, "BCName");
    int gridIndex = grid->GetZoneID();
    StructBCSet *structBCSet = grid->GetStructBCSet();

    char *BCName = ReadCharData(grpData, "BCName");
    string *BCNameList = new string[nBCRegion];
    uint_long count = 0;

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        while (BCName[count] != '\0')
        {
            BCNameList[iBCRegion].push_back(BCName[count]);
            ++ count;
        }
        ++ count;
    }

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcregion = new StructBC(gridIndex, iBCRegion);
        structBCSet->SetBCRegion(iBCRegion, bcregion);

        bcregion->SetBCName(BCNameList[iBCRegion]);
    }

    delete [] BCNameList;    BCNameList = nullptr;
    delete [] BCName;    BCName = nullptr;
    H5Gclose(grpData);
}

void IO_HDF5Read::ReadStrPartitionInfo(hid_t loc_id, Grid *gridIn)
{
    bool dataExist = CheckDataExist(loc_id, "PartitionInfo");
    if (!dataExist)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int ni    = grid->GetNI();
        int nj    = grid->GetNJ();
        int nk    = grid->GetNK();
        int iZone = grid->GetZoneID();

        int ordinaryDimStartIndex[3], ordinaryDimEndIndex[3];
        ordinaryDimStartIndex[0] = 1; ordinaryDimEndIndex[0] = ni;
        ordinaryDimStartIndex[1] = 1; ordinaryDimEndIndex[1] = nj;
        ordinaryDimStartIndex[2] = 1; ordinaryDimEndIndex[2] = nk;

        grid->SetOrdinaryGridIndex(iZone);
        grid->SetOrdinaryDimStartIndex(&ordinaryDimStartIndex[0]);
        grid->SetOrdinaryDimEndIndex(&ordinaryDimEndIndex[0]);

        return;
    }

    hid_t grpData;
    grpData = OpenGroup(loc_id, "PartitionInfo");
    StructGrid *grid = StructGridCast(gridIn);

    int ordinaryGridIndex = -1;
    int *ordinaryDimStartIndex = new int[3];
    int *ordinaryDimEndIndex   = new int[3];

    ReadData(grpData, &ordinaryGridIndex   , "OrdinaryGrid");
    ReadData(grpData, ordinaryDimStartIndex, "OrdinaryDimStart");
    ReadData(grpData, ordinaryDimEndIndex  , "OrdinaryDimEnd");

    grid->SetOrdinaryGridIndex(ordinaryGridIndex);
    grid->SetOrdinaryDimStartIndex(ordinaryDimStartIndex);
    grid->SetOrdinaryDimEndIndex(ordinaryDimEndIndex);

    delete [] ordinaryDimStartIndex;    ordinaryDimStartIndex = nullptr;
    delete [] ordinaryDimEndIndex;    ordinaryDimEndIndex = nullptr;

    H5Gclose(grpData);
}

void IO_HDF5Read::ReadUnstrPartitionInfo(hid_t loc_id, Grid *gridIn)
{
    bool dataExist = CheckDataExist(loc_id, "PartitionInfo");
    if (!dataExist)
    {
        return;
    }

    hid_t grpData;
    grpData = OpenGroup(loc_id, "PartitionInfo");
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int ordinaryGridIndex = -1;
    ReadData(grpData, &ordinaryGridIndex, "OrdinaryGrid");
    grid->SetOrdinaryGridIndex(ordinaryGridIndex);

    int *ordinaryNodeIndex  = ReadIntOneRow(grpData, "OrdinaryNodeIndex");
    int *ordinaryFaceIndex  = ReadIntOneRow(grpData, "OrdinaryFaceIndex");
    int *ordinaryCellIndex  = ReadIntOneRow(grpData, "OrdinaryCellIndex");
    grid->SetOrdinaryNodeIndex(ordinaryNodeIndex);
    grid->SetOrdinaryFaceIndex(ordinaryFaceIndex);
    grid->SetOrdinaryCellIndex(ordinaryCellIndex);

    H5Gclose(grpData);
}

void IO_HDF5Read::ReadUnstructGrid(hid_t loc_id, int iZone)
{
    PrintToWindow("Reading Unstructured Grid Of Zone ", iZone, "...\n");

    hid_t  grpData;
    string grpName;

    GridID *index = new GridID(iZone + zoneStart);
    Grid *grid = new UnstructGrid();

    grid->InitGrid(index, 0, PHSPACE::GetDim(), UNSTRUCTGRID);
    grid->SetIBlock(0);

    ostringstream oss;
    oss << "Grid-" << iZone;
    grpName = oss.str();
    grpData = OpenGroup(loc_id, grpName);

    int *iDimensions = ReadIntOneRow(grpData, "iDimensions");
    grid->SetNTotalNode(iDimensions[0]);
    grid->SetNTotalFace(iDimensions[1]);
    grid->SetNTotalCell(iDimensions[2]);

    PrintToWindow("  Grid Dimension  : ", grid->GetDim(), "\n");
    PrintToWindow("  Number of Points: ", iDimensions[0], "\n");
    PrintToWindow("  Number of Faces : ", iDimensions[1], "\n");
    PrintToWindow("  Number of Cells : ", iDimensions[2], "\n");

    ReadGridCoordinates(grpData, grid);
    ReadFaceTopology(grpData, grid);
    ReadInterfaceInfo(grpData, grid);
    ReadCellTopology(grpData, grid);
    ReadUnstrBCName(grpData, grid);
    ReadInterpointInfo(grpData, grid);
    ReadUnstrPartitionInfo(grpData, grid);
    ReadVolumeCondition(grpData, grid);

    grid->SetIBlock(presentGridGroup);
    grid->ReviseNeighborZoneIndex(zoneStart);

    AddGrid(grid);

    delete [] iDimensions;    iDimensions = nullptr;
    H5Gclose(grpData);
}

void IO_HDF5Read::ReadFaceTopology(hid_t loc_id, Grid *gridIn)
{
    PrintToWindow("    Reading face topology ...\n");
    hid_t grpData;

    grpData = OpenGroup(loc_id, "FaceTopology");

    UnstructGrid *grid = UnstructGridCast(gridIn);

    int *nodeNumberOfEachFace = ReadIntOneRow(grpData, "nodeNumberOfEachFace");
    grid->SetNodeNumberOfEachFace(nodeNumberOfEachFace);

    int *face2Node = ReadIntOneRow(grpData, "face2Node");
    grid->SetFace2Node(face2Node);

    int *leftCellOfFace  = ReadIntOneRow(grpData, "leftCellOfFace");
    int *rightCellOfFace = ReadIntOneRow(grpData, "rightCellOfFace");
    grid->SetLeftCellOfFace (leftCellOfFace);
    grid->SetRightCellOfFace(rightCellOfFace);

    int nTotalFace =grid->GetNTotalFace();
    int *nodePosi = new int[nTotalFace + 1];
    int nodeSum = 0;
    nodePosi[0] = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        nodeSum += nodeNumberOfEachFace[iFace];
        nodePosi[iFace+1] = nodeSum;
    }

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        if (leftCellOfFace[iFace] < 0)
        {
            std::reverse(face2Node + nodePosi[iFace], face2Node + nodePosi[iFace+1]);
            SWAP(leftCellOfFace[iFace], rightCellOfFace[iFace]);
        }
    }
    delete [] nodePosi;    nodePosi = nullptr;

    int nBoundFace = 0;
    ReadData(grpData, &nBoundFace, "nBoundFace");
    grid->SetNBoundFace(nBoundFace);
    PrintToWindow("  Number of BC Faces : ", nBoundFace, "\n" );

    grid->SpecifyRightCellofBC();

    H5Gclose(grpData);
}

void IO_HDF5Read::ReadCellTopology(hid_t loc_id, Grid *gridIn)
{
    PrintToWindow("    Reading cell-to-node ...\n");
    hid_t grpData;

    grpData = OpenGroup(loc_id, "CellTopology");

    UnstructGrid *grid = UnstructGridCast(gridIn);

    int *cell2Node = ReadIntOneRow(grpData, "cell2Node");
    int *nodeNumberOfEachCell = ReadIntOneRow(grpData, "nodeNumberOfEachCell");

    grid->SetNodeNumberOfEachCell(nodeNumberOfEachCell);
    grid->SetCell2Node(cell2Node);

    H5Gclose(grpData);
}

void IO_HDF5Read::ReadUnstrBCName(hid_t loc_id, Grid *gridIn)
{
    PrintToWindow("    Reading boundary name ...\n");
    hid_t grpData;

    grpData = OpenGroup(loc_id, "BCName");

    UnstructGrid *grid = UnstructGridCast(gridIn);
    char *BCName = ReadCharData(grpData, "BCName");

    grpData = OpenGroup(loc_id, "FaceTopology");
    int *BcType = ReadIntOneRow(grpData, "BcType");

    int nBoundFace = grid->GetNBoundFace();
    string *BCNameList = new string[nBoundFace];
   
    int count = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        while (BCName[count] != '\0')
        {
            BCNameList[iFace].push_back(BCName[count]);
            ++ count;
        }
        ++ count;
    }

    UnstructBCSet **bcr = new UnstructBCSet *[nBoundFace];

    //! Create UnstructbcRegion.
    set<string> bcNameSet;
    for (int iFace = 0; iFace < nBoundFace; iFace++)
    {   
        bcr[iFace] = new UnstructBCSet();
        bcr[iFace]->SetKey(BcType[iFace]);
        bcr[iFace]->SetBCName(BCNameList[iFace]);

        bcNameSet.insert(BCNameList[iFace]);
    }
    grid->SetBCRecord(bcr);

    int nBCRegionUnstruct = static_cast<int>(bcNameSet.size());
    
    grid->CreateUnstructBCSet(nBCRegionUnstruct);

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = new int[nBoundFace];
    set<string>::iterator iter;
    count = 0;
    for (iter = bcNameSet.begin(); iter != bcNameSet.end(); iter++)
    {
        UnstructBC *unstructBC = new UnstructBC(count);
        unstructBCSet->SetBCRegion(count, unstructBC);
        unstructBC->SetBCName(*iter);

        for (int iFace = 0; iFace < nBoundFace; iFace++)
        {
            if (BCNameList[iFace] == unstructBC->GetBCName())
            {
                unstructBC->SetBCType(BcType[iFace]);
                bcRegionIDofBCFace[iFace] = count;
                vector<int> *faceIndex = unstructBC->GetFaceIndex();
                faceIndex->push_back(iFace);
            }
           
        }
        count++;
    }
    unstructBCSet->SetBCFaceInfo(bcRegionIDofBCFace);
    delete [] BcType;    BcType = nullptr;
    delete [] BCNameList;    BCNameList = nullptr;
    delete [] BCName;    BCName = nullptr;
    H5Gclose(grpData);
}

void IO_HDF5Read::ReadInterpointInfo(hid_t loc_id, Grid *gridIn)
{
    if (isP2PMode)
    {
        return;
    }

    if (numberOfGridGroups != 1)
    {
        return;
    }

    InterpointInformation *ipointinfo = 0;
    gridIn->SetInterpointInfo(ipointinfo);

    hid_t grpData;
    grpData = OpenGroup(loc_id, "InterpointInfo");

    int nIPoint;
    ReadData(grpData, &nIPoint, "nIPoint");
    if (!nIPoint)
    {
        H5Gclose(grpData);
        return;
    }

    PrintToWindow("    Reading interpoint ...\n");

    ipointinfo = new InterpointInformation(nIPoint,gridIn);
    gridIn->SetInterpointInfo(ipointinfo);

    int *interPoint2ZoneID = ipointinfo->GetInterPoint2ZoneID();
    int *interPoint2InterPointID = ipointinfo->GetInterPoint2InterPointID();
    int *interPoint2GlobalPoint = ipointinfo->GetInterPoint2GlobalPoint();
    int *cellNumberOfInterPoint = ipointinfo->GetCellNumberOfInterPoint();
    int *totalZonesOfInterPoint = ipointinfo->GetTotalZonesOfInterPoint();
    int *labelOfInterPoint = ipointinfo->GetLabelOfInterPoint();

    ReadData(grpData, interPoint2ZoneID      , "interPoint2ZoneID");
    ReadData(grpData, interPoint2InterPointID, "interPoint2InterPointID");
    ReadData(grpData, interPoint2GlobalPoint , "interPoint2GlobalPoint");
    ReadData(grpData, cellNumberOfInterPoint , "cellNumberOfInterPoint");
    ReadData(grpData, totalZonesOfInterPoint , "totalZonesOfInterPoint");
    ReadData(grpData, labelOfInterPoint      , "labelOfInterPoint");

    H5Gclose(grpData);

    zoneConnectivityForPoint = new ZoneConnectivityForPoint();
}

bool IO_HDF5Read::IsZoneLayoutInitialized()
{
    if (this->blockIdx || this->blockProc || this->blockType)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void IO_HDF5Read::RedistributeZonesIfNeed()
{
    int nZones = GetNZones();
    int numberOfProcessor = PHMPI::GetNumberOfProcessor();

    int *blockType     = GetBlockType();
    int *blockProc     = GetBlockProc();
    int *blockProcGrid = GetBlockProcGrid();
    int *blockProcDump = GetBlockProcDump();

    int zoneDistributionMethod = 0;
    int maxProc = blockProc[0];

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        maxProc = MAX(maxProc, blockProc[iZone]);
    }
    maxProc += 1;

    if (numberOfProcessor == maxProc)
    {
        zoneDistributionMethod = DETERMIN_BY_PARTITION;
    }
    else
    {
        zoneDistributionMethod = REDISTRIBUTION;
    }
    PHMPI::SetZoneDistributionMethod(zoneDistributionMethod);

    ostringstream oss;
    oss << "Max processor ID : " << maxProc << " of total zone " << nZones;
    if (zoneDistributionMethod == DETERMIN_BY_PARTITION)
    {
        oss << ", Zones' distribution is specified by partition ...";
    }
    else
    {
        oss << ", Zones' distribution is specified randomly ...";
    }
    oss << endl;
    WriteLogFile(oss);

    if (zoneDistributionMethod == REDISTRIBUTION)
    {
        if (nZones % numberOfProcessor == 0)
        {
            int nZonesPerProcessor = nZones / numberOfProcessor;
            int processorID = 0;
            for (int iZone = 0; iZone < nZones; ++ iZone)
            {
                blockProc[iZone] = processorID;
                if ((iZone + 1) % nZonesPerProcessor == 0)
                {
                    processorID ++;
                }
            }
        }
        else
        {
            PrintToWindow("Waining: nZones % numberOfProcessor != 0 !!");
            for (int iZone = 0; iZone < nZones; ++ iZone)
            {
                blockProc[iZone] = iZone % numberOfProcessor;
            }
        }
    }
    else
    {
        if (numberOfProcessor == 1)
        {
            if (IsConvertGridToMixGrid())
            {
                if (blockType[0] == PHSPACE::STRUCTGRID)
                {
                    int maxproc = blockProc[0];
                    for (int iZone = 0; iZone < nZones; ++ iZone)
                    {
                        if (blockProc[iZone] > maxproc)
                        {
                            maxproc = blockProc[iZone];
                        }
                    }
                    maxproc += 1;
                    PHMPI::SetNumberOfProcStructUsed(maxproc);
                }
                else
                {
                    for (int iZone = 0; iZone < nZones; ++ iZone)
                    {
                        blockProcDump[iZone] += PHMPI::GetNumberOfProcStructUsed();
                    }
                }
            }

            for (int iZone = 0; iZone < nZones; ++ iZone)
            {
                blockProc[iZone] = 0;
            }
        }
    }

    if (!blockProcGrid)
    {
        blockProcGrid = new int[nZones];
    }
    SetBlockProcGrid(blockProcGrid);

    int numberOfGridProcessor = PHMPI::GetNumberOfGridProcessor();
    if (numberOfGridProcessor != 0)
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            blockProcGrid[iZone] = iZone % numberOfGridProcessor + numberOfProcessor;
        }
    }
    else
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            blockProcGrid[iZone] = -1;
        }
    }
}

void IO_HDF5Read::SetGlobalZoneLayout()
{
    if (!PHMPI::GetZoneGridType())
    {
        PHMPI::SetNumberOfGlobalZones(numberOfZones);
        PHMPI::CreateZoneProcessorIDDump(numberOfZones);
        PHMPI::CreateZoneProcessorID(numberOfZones);
        PHMPI::CreateZoneGridID(numberOfZones);
        PHMPI::CreateZoneGridType(numberOfZones);
        PHMPI::CreateZoneFileID(numberOfZones);
    }

    int *globalBlockProc     = PHMPI::GetZoneProcessorID();
    int *globalBlockProcGrid = PHMPI::GetZoneProcessorIDForGrid();
    int *globalBlockType     = PHMPI::GetZoneGridType();
    int *globalBlockIdx      = PHMPI::GetZoneGridID();
    int *globalBlockFile     = PHMPI::GetZoneFileID();
    int *globalBlockProcDump = PHMPI::GetZoneProcessorIDDump();

    set<int> zoneIDinCurrentProc;
    int myid = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        globalBlockProc    [iZone] = blockProc    [iZone];
        globalBlockProcGrid[iZone] = blockProcGrid[iZone];
        globalBlockType    [iZone] = blockType    [iZone];
        globalBlockIdx     [iZone] = blockIdx     [iZone];
        globalBlockFile    [iZone] = blockFileProc[iZone];

        if (blockProc[iZone] == myid || blockProcGrid[iZone] == myid)
        {
            zoneIDinCurrentProc.insert(iZone);
        }
    }

    ostringstream oss;
    oss << "Zone index in current processor : ";
    for (set<int>::iterator iter = zoneIDinCurrentProc.begin(); iter != zoneIDinCurrentProc.end(); ++ iter)
    {
        int partID = *iter;
        oss << partID << " ";
    }
    oss << endl;
    WriteLogFile(oss.str());

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        globalBlockProcDump[iZone] = blockProcDump[iZone];
    }
}

void IO_HDF5Read::ClassifyGridSystem()
{
    int nZones = PHMPI::GetNumberofGlobalZones();

    set<int> gridType;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        gridType.insert(PHMPI::GetZoneGridType(iZone));
    }

    int sysGridType = GetSystemGridType(gridType);
    GlobalDataBase::UpdateData("sys_gridtype", &sysGridType, PHINT, 1);

    vector<string> sysGridTypeName;
    sysGridTypeName.push_back("UNSTRUCTGRID");
    sysGridTypeName.push_back("STRUCTGRID");
    sysGridTypeName.push_back("MIXGRID");

    PrintToWindow("Grid Type: ", sysGridTypeName[sysGridType], "\n");
}

void IO_HDF5Read::AddGrid(Grid *grid)
{
    grids.push_back(grid);
}

void IO_HDF5Read::CreateZones()
{
    zoneStart = 0;
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int globalZoneIndex = zoneStart + iZone;

        int recvProc     = PHMPI::GetZoneProcessorID       (globalZoneIndex);
        int recvProcGrid = PHMPI::GetZoneProcessorIDForGrid(globalZoneIndex);

        if (myProcID == recvProc || myProcID == recvProcGrid)
        {
            int nLocalZones = PHMPI::GetNumberofLocalZones();

            Grid *grid = grids[iZone];
            grid->GetGridID()->SetLocalIndex(nLocalZones);

            Zone *zone = new Zone(globalZoneIndex);
            zone->AddGridAndCopyZonePara(grid);
            zone->AddGridToGlobal(globalZoneIndex);

            region->AddZone(zone);
        }
        else
        {
            region->AddZone(0);
            continue;
        }
    }
}

IO_ReadHDF5Grid::IO_ReadHDF5Grid(const string &gridFileNameIn)
{
    this->gridFileName  = gridFileNameIn;
    this->fileID        = -1;
    this->numberOfZones = 0;
    this->fileOpened    = false;
    zoneIndexs          = 0;
    zoneTypes           = 0;
    zoneProcs           = 0;
}

IO_ReadHDF5Grid::~IO_ReadHDF5Grid()
{
    if (fileOpened)
    {
        CloseFile();
    }

    if (zoneIndexs)
    {
        delete [] zoneIndexs;    zoneIndexs = nullptr;
    }

    if (zoneTypes)
    {
        delete [] zoneTypes;     zoneTypes = nullptr;
    }

    if (zoneProcs)
    {
        delete [] zoneProcs;     zoneProcs = nullptr;
    }

}

void IO_ReadHDF5Grid::OpenFile()
{
    if (fileOpened)
    {
        ostringstream errorMess;
        errorMess << "  Error: This grid file have been opened, file name = " << gridFileName << endl;
        TK_Exit::ExceptionExit(errorMess.str());
    }

    fileID = OpenHDF5File(this->gridFileName);
    if (fileID < 0)
    {
        ostringstream errorMess;
        errorMess << "  Error: This grid file can not open, file name = " << gridFileName << endl;
        TK_Exit::ExceptionExit(errorMess.str());
    }
    fileOpened = true;
}

void IO_ReadHDF5Grid::CloseFile()
{
    if (!fileOpened)
    {
        ostringstream errorMess;
        errorMess << "  Error: This grid file is not opened, file name = " << gridFileName << endl;
        TK_Exit::ExceptionExit(errorMess.str());
    }
    H5Fclose(fileID);
    fileOpened = false;
}

void IO_ReadHDF5Grid::GetNumberOfZones(int &nZones)
{
    if (!fileOpened)
    {
        ostringstream errorMess;
        errorMess << "  Error: This grid file is not opened, file name = " << gridFileName << endl;
        TK_Exit::ExceptionExit(errorMess.str());
    }

    hid_t grpData;
    grpData = OpenGroup(fileID, "Information");
    ReadData(grpData, &numberOfZones, "nBlocks");
    nZones = numberOfZones;
    H5Gclose(grpData);
}

void IO_ReadHDF5Grid::GetZoneIndex(int *zoneIndexOut)
{
    if (!fileOpened)
    {
        ostringstream errorMess;
        errorMess << "  Error: This grid file is not opened, file name = " << gridFileName << endl;
        TK_Exit::ExceptionExit(errorMess.str());
    }

    hid_t grpData;
    grpData = OpenGroup(fileID, "Information");
    if (zoneIndexs)
    {
        delete [] zoneIndexs;    zoneIndexs = nullptr;
    }
    zoneIndexs = ReadIntOneRow(grpData, "block_idx");
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        zoneIndexOut[iZone] = zoneIndexs[iZone];
    }
    H5Gclose(grpData);
}

void IO_ReadHDF5Grid::GetZoneType(int *zoneTypeOut)
{
    if (!fileOpened)
    {
        ostringstream errorMess;
        errorMess << "  Error: This grid file is not opened, file name = " << gridFileName << endl;
        TK_Exit::ExceptionExit(errorMess.str());
    }

    hid_t grpData;
    grpData = OpenGroup(fileID, "Information");
    if (zoneTypes)
    {
        delete [] zoneTypes;    zoneTypes = nullptr;
    }
    zoneTypes = ReadIntOneRow(grpData, "block_type");
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        zoneTypeOut[iZone] = zoneTypes[iZone];
    }
    H5Gclose(grpData);
}

void IO_ReadHDF5Grid::GetZoneProc(int *zoneProcOut)
{
    if (!fileOpened)
    {
        ostringstream errorMess;
        errorMess << "  Error: This grid file is not opened, file name = " << gridFileName << endl;
        TK_Exit::ExceptionExit(errorMess.str());
    }

    hid_t grpData;
    grpData = OpenGroup(fileID, "Information");
    if (zoneProcs)
    {
        delete [] zoneProcs;    zoneProcs = nullptr;
    }
    zoneProcs = ReadIntOneRow(grpData, "block_proc");
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        zoneProcOut[iZone] = zoneProcs[iZone];
    }
    H5Gclose(grpData);
}

void IO_ReadHDF5Grid::GetZoneSize(int zoneIndex, int *iDimensionsOut)
{
    hid_t  gridGroup;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);
    int *iDimensions = ReadIntOneRow(gridGroup, "iDimensions");
    for (int iSize = 0; iSize < THREE_D; ++ iSize)
    {
        iDimensionsOut[iSize] = iDimensions[iSize];
    }
    delete [] iDimensions;    iDimensions = nullptr;
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetZoneCoord(int zoneIndex, int iCoord, int dataSize, RDouble *coordOut)
{
    hid_t  gridGroup, grpData;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);
    grpData   = OpenGroup(gridGroup, "GridCoordinates");

    RDouble *coordinate = 0;
    if (X_DIR == iCoord + 1)
    {
        coordinate = ReadDoubleOneRow(grpData, "CoordinateX");
    }
    else if (Y_DIR == iCoord + 1)
    {
        coordinate = ReadDoubleOneRow(grpData, "CoordinateY");
    }
    else
    {
        coordinate = ReadDoubleOneRow(grpData, "CoordinateZ");
    }
    for (int iData = 0; iData < dataSize; ++ iData)
    {
        coordOut[iData] = coordinate[iData];
    }
    delete [] coordinate;    coordinate = nullptr;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetNIFace(int zoneIndex, int &nIFaceOut)
{
    hid_t  gridGroup, grpData;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);
    grpData   = OpenGroup(gridGroup, "InterfaceInfo");
    int nIFace = 0;
    ReadData(grpData, &nIFace, "nIFace");
    nIFaceOut = nIFace;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetInterfaceInfo(int zoneIndex, InterfaceInfo *interfaceInfo)
{
    hid_t  gridGroup, grpData;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);
    grpData   = OpenGroup(gridGroup, "InterfaceInfo");

    int *interFace2ZoneID           = interfaceInfo->GetInterFace2ZoneID();
    int *interFace2InterFaceID      = interfaceInfo->GetInterFace2InterFaceID();
    int *interFace2BoundaryFace     = interfaceInfo->GetInterFace2BoundaryFace();
    int *interFace2ZoneIDTemp       = ReadIntOneRow(grpData, "interFace2ZoneID");
    int *interFace2InterFaceIDTemp  = ReadIntOneRow(grpData, "interFace2InterFaceID");
    int *interFace2BoundaryFaceTemp = ReadIntOneRow(grpData, "interFace2BoundaryFace");
    int nIFace                      = interfaceInfo->GetNIFace();
    for (int iFace = 0; iFace < nIFace; ++ iFace)
    {
        interFace2ZoneID[iFace]       = interFace2ZoneIDTemp[iFace];
        interFace2InterFaceID[iFace]  = interFace2InterFaceIDTemp[iFace];
        interFace2BoundaryFace[iFace] = interFace2BoundaryFaceTemp[iFace];
    }

    delete [] interFace2ZoneIDTemp;          interFace2ZoneIDTemp       = nullptr;
    delete [] interFace2InterFaceIDTemp;     interFace2InterFaceIDTemp  = nullptr;
    delete [] interFace2BoundaryFaceTemp;    interFace2BoundaryFaceTemp = nullptr;

    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetNBCRegion(int zoneIndex, int &nBCRegionOut)
{
    hid_t  gridGroup, grpData;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);
    grpData   = OpenGroup(gridGroup, "BCRegion");
    int nBCRegion = 0;
    ReadData(grpData, &nBCRegion, "nBCRegion");
    nBCRegionOut = nBCRegion;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetBCName(int zoneIndex, StructBCSet *compositeBCRegion)
{
    hid_t  gridGroup, grpData;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);
    grpData   = OpenGroup(gridGroup, "BCName");
    char *BCName   = ReadCharData(grpData, "BCName");
    int  nBCRegion = compositeBCRegion->GetnBCRegion();
    string *BCNameList = new string [nBCRegion];
    uint_long count = 0;

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        while (BCName[count] != '\0')
        {
            BCNameList[iBCRegion].push_back(BCName[count]);
            ++ count;
        }
        ++ count;
    }

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcRegion = new StructBC(zoneIndex, iBCRegion);
        compositeBCRegion->SetBCRegion(iBCRegion, bcRegion);

        bcRegion->SetBCName(BCNameList[iBCRegion]);
    }

    delete [] BCNameList;    BCNameList = nullptr;
    delete [] BCName;        BCName     = nullptr;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetBCRegion(int zoneIndex, StructBCSet *compositeBCRegion)
{
    hid_t  gridGroup, grpData;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);
    grpData   = OpenGroup(gridGroup, "BCRegion");
    int nBCRegion = compositeBCRegion->GetnBCRegion();
    int **BCInfor = ReadIntTwoRow(grpData, "BCInfor");
    int iMin, iMax, jMin, jMax, kMin, kMax, bcType;
    int sLr, sNd;
    int iBCinfo = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        int index = 0;
        StructBC *bcRegion = compositeBCRegion->GetBCRegion(iBCRegion);
        iMin   = BCInfor[iBCinfo][index ++];
        iMax   = BCInfor[iBCinfo][index ++];
        jMin   = BCInfor[iBCinfo][index ++];
        jMax   = BCInfor[iBCinfo][index ++];
        kMin   = BCInfor[iBCinfo][index ++];
        kMax   = BCInfor[iBCinfo][index ++];
        bcType = BCInfor[iBCinfo][index ++];
        sLr    = BCInfor[iBCinfo][index ++];
        sNd    = BCInfor[iBCinfo][index ++];
        iBCinfo ++;

        //! Restore periodic boundary bcType and correct read information only for PARTITION_GRID task.
        if (GetTaskCode() == PARTITION_GRID && PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType") != 0)
        {
            if (bcRegion->GetBCName() == "Periodic_up" || bcRegion->GetBCName() == "Periodic_down")
            {
                bcType = PHENGLEI::PERIODIC;
                iBCinfo ++;
            }
        }

        bcRegion->SetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);
        bcRegion->SetBCType(bcType);
        bcRegion->SetFaceLeftOrRightIndex(sLr);
        bcRegion->SetFaceDirection(sNd);
        bcRegion->InitFaceDirectionIndex();

        if (-1 == bcType)
        {
            index = 0;
            iMin  = BCInfor[iBCinfo][index ++];
            iMax  = BCInfor[iBCinfo][index ++];
            jMin  = BCInfor[iBCinfo][index ++];
            jMax  = BCInfor[iBCinfo][index ++];
            kMin  = BCInfor[iBCinfo][index ++];
            kMax  = BCInfor[iBCinfo][index ++];
            bcRegion->SetTargetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);
            int neighborZoneIndex = BCInfor[iBCinfo][index ++];
            bcRegion->SetTargetRegionBlock(neighborZoneIndex);
            iBCinfo ++;
        }
    }

    DelPointer2(BCInfor);
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetUnsNodeNumOfEachFace(int zoneIndex, int dataSize, int *nodeNumberOfEachFaceOut)
{
    hid_t  gridGroup, grpData;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);
    grpData   = OpenGroup(gridGroup, "FaceTopology");
    int *nodeNumberOfEachFace = ReadIntOneRow(grpData, "nodeNumberOfEachFace");
    for (int iData = 0; iData < dataSize; ++ iData)
    {
        nodeNumberOfEachFaceOut[iData] = nodeNumberOfEachFace[iData];
    }
    delete [] nodeNumberOfEachFace;    nodeNumberOfEachFace = nullptr;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetUnsFace2Node(int zoneIndex, int dataSize, int *face2NodeOut)
{
    hid_t  gridGroup, grpData;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);
    grpData   = OpenGroup(gridGroup, "FaceTopology");
    int *face2Node = ReadIntOneRow(grpData, "face2Node");
    for (int iData = 0; iData < dataSize; ++ iData)
    {
        face2NodeOut[iData] = face2Node[iData];
    }
    delete [] face2Node;    face2Node = nullptr;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetUnsLeftCellOfFace(int zoneIndex, int dataSize, int *leftCellOfFaceOut)
{
    hid_t  gridGroup, grpData;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);
    grpData   = OpenGroup(gridGroup, "FaceTopology");
    int *leftCellOfFace = ReadIntOneRow(grpData, "leftCellOfFace");
    for (int iData = 0; iData < dataSize; ++ iData)
    {
        leftCellOfFaceOut[iData] = leftCellOfFace[iData];
    }
    delete [] leftCellOfFace;    leftCellOfFace = nullptr;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetUnsRightCellOfFace(int zoneIndex, int dataSize, int *rightCellOfFaceOut)
{
    hid_t  gridGroup, grpData;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);
    grpData   = OpenGroup(gridGroup, "FaceTopology");
    int *rightCellOfFace = ReadIntOneRow(grpData, "rightCellOfFace");
    for (int iData = 0; iData < dataSize; ++ iData)
    {
        rightCellOfFaceOut[iData] = rightCellOfFace[iData];
    }
    delete [] rightCellOfFace;    rightCellOfFace = nullptr;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetUnsNBoundFace(int zoneIndex, int &nBoundFaceOut)
{
    hid_t  gridGroup, grpData;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);
    grpData   = OpenGroup(gridGroup, "FaceTopology");
    int nBoundFace = 0;
    ReadData(grpData, &nBoundFace, "nBoundFace");
    nBoundFaceOut = nBoundFace;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetUnsBCInfo(int zoneIndex, int dataSize, UnstructBCSet **bcInfoOut)
{
    hid_t  gridGroup, grpData;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);
    grpData   = OpenGroup(gridGroup, "FaceTopology");
    int *bcType = ReadIntOneRow(grpData, "BcType");
    for (int iData = 0; iData < dataSize; ++ iData)
    {
        bcInfoOut[iData]->SetKey(bcType[iData]);
    }
    delete [] bcType;    bcType = nullptr;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetUnsBCName(int zoneIndex, int dataSize, UnstructBCSet **bcInfoOut)
{
    hid_t  gridGroup, grpData;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);
    grpData   = OpenGroup(gridGroup, "BCName");
    char *bcName = ReadCharData(grpData, "BCName");
    string *bcNameList = new string [dataSize];

    uint_long count = 0;
    for (int iFace = 0; iFace < dataSize; ++ iFace)
    {
        while (bcName[count] != '\0')
        {
            bcNameList[iFace].push_back(bcName[count]);
            ++ count;
        }
        ++ count;
    }
    for (int iData = 0; iData < dataSize; ++ iData)
    {
        bcInfoOut[iData]->SetBCName(bcNameList[iData]);
    }

    delete [] bcName;        bcName     = nullptr;
    delete [] bcNameList;    bcNameList = nullptr;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetUnsNodeNumOfEachCell(int zoneIndex, int dataSize, int *nodeNumberOfEachCellOut)
{
    hid_t  gridGroup, grpData;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);
    grpData   = OpenGroup(gridGroup, "CellTopology");
    int *nodeNumberOfEachCell = ReadIntOneRow(grpData, "nodeNumberOfEachCell");
    for (int iData = 0; iData < dataSize; ++ iData)
    {
        nodeNumberOfEachCellOut[iData] = nodeNumberOfEachCell[iData];
    }
    delete [] nodeNumberOfEachCell;    nodeNumberOfEachCell = nullptr;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetUnsCell2Node(int zoneIndex, int dataSize, int *cell2NodeOut)
{
    hid_t  gridGroup, grpData;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);
    grpData   = OpenGroup(gridGroup, "CellTopology");
    int *cell2Node = ReadIntOneRow(grpData, "cell2Node");
    for (int iData = 0; iData < dataSize; ++ iData)
    {
        cell2NodeOut[iData] = cell2Node[iData];
    }
    delete [] cell2Node;    cell2Node = nullptr;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetUnstructGridInfo(Grid *gridOut)
{
    int          gridIndex = gridOut->GetZoneID();
    UnstructGrid *unsGrid  = UnstructGridCast(gridOut);

    //! Read grid size.
    int *iDimensions = new int [THREE_D];
    GetZoneSize(gridIndex, iDimensions);
    int nTotalNode = iDimensions [0];
    int nTotalFace = iDimensions [1];
    int nTotalCell = iDimensions [2];
    unsGrid->SetNTotalNode(nTotalNode);
    unsGrid->SetNTotalFace(nTotalFace);
    unsGrid->SetNTotalCell(nTotalCell);
    delete [] iDimensions;    iDimensions = nullptr;

    //! Read grid coordiate.
    RDouble **coordinates = new RDouble * [THREE_D];
    for (int iDim = 0; iDim < THREE_D; ++ iDim)
    {
        coordinates[iDim] = new RDouble [nTotalNode];
        GetZoneCoord(gridIndex, iDim, nTotalNode, coordinates[iDim]);
    }
    unsGrid->SetX(coordinates[0]);
    unsGrid->SetY(coordinates[1]);
    unsGrid->SetZ(coordinates[2]);
    unsGrid->RotateAboutAxis();
    unsGrid->ComputeMinMaxBox();
    delete [] coordinates;    coordinates = nullptr;

    //! Read face topo.
    int *nodeNumberOfEachFace = new int [nTotalFace];
    GetUnsNodeNumOfEachFace(gridIndex, nTotalFace, nodeNumberOfEachFace);
    unsGrid->SetNodeNumberOfEachFace(nodeNumberOfEachFace);

    int face2NodeLength = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        face2NodeLength += nodeNumberOfEachFace[iFace];
    }
    int *face2Node = new int [face2NodeLength];
    GetUnsFace2Node(gridIndex, face2NodeLength, face2Node);
    unsGrid->SetFace2Node(face2Node);

    int *leftCellOfFace = new int [nTotalFace];
    GetUnsLeftCellOfFace(gridIndex, nTotalFace, leftCellOfFace);
    unsGrid->SetLeftCellOfFace(leftCellOfFace);

    int *rightCellOfFace = new int [nTotalFace];
    GetUnsRightCellOfFace(gridIndex, nTotalFace, rightCellOfFace);
    unsGrid->SetRightCellOfFace(rightCellOfFace);

    //! Read cell topo.
    int *nodeNumberOfEachCell = new int [nTotalCell];
    GetUnsNodeNumOfEachCell(gridIndex, nTotalCell, nodeNumberOfEachCell);
    unsGrid->SetNodeNumberOfEachCell(nodeNumberOfEachCell);
    int cell2NodeLength = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        cell2NodeLength += nodeNumberOfEachCell[iCell];
    }
    int *cell2Node = new int [cell2NodeLength];
    GetUnsCell2Node(gridIndex, cell2NodeLength, cell2Node);
    unsGrid->SetCell2Node(cell2Node);

    //! Read interface information.
    int nIFace = 0;
    GetNIFace(gridIndex, nIFace);
    unsGrid->SetNIFace(nIFace);
    InterfaceInfo *interfaceInfo = 0;
    if (nIFace > 0)
    {
        interfaceInfo = new InterfaceInfo(nIFace, unsGrid);
        GetInterfaceInfo(gridIndex, interfaceInfo);
    }
    unsGrid->SetInterfaceInfo(interfaceInfo);

    //! Read boundary information.
    int nBoundFace = 0;
    GetUnsNBoundFace(gridIndex, nBoundFace);
    unsGrid->SetNBoundFace(nBoundFace);
    UnstructBCSet **bcInfo = new UnstructBCSet * [nBoundFace];
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        bcInfo[iFace] = new UnstructBCSet();
    }
    GetUnsBCInfo(gridIndex, nBoundFace, bcInfo);
    GetUnsBCName(gridIndex, nBoundFace, bcInfo);
    unsGrid->SetBCRecord(bcInfo);
    unsGrid->ConstructBCRegion();
    unsGrid->SpecifyRightCellofBC();
    int *keyActiveOfCells = new int [nTotalCell];
    SetField(keyActiveOfCells, 1, nTotalCell);
    unsGrid->SetBlankIndex(keyActiveOfCells);

    int ordinaryGridIndex = -1;
    int *ordinaryNodeIndex = new int [nTotalNode];
    int *ordinaryFaceIndex = new int [nTotalFace];
    int *ordinaryCellIndex = new int [nTotalCell];

    ReadUnsOrdinaryGridIndex(gridIndex, ordinaryGridIndex);
    ReadUnsOrdinaryNodeIndex(gridIndex, ordinaryNodeIndex, nTotalNode);
    ReadUnsOrdinaryFaceIndex(gridIndex, ordinaryFaceIndex, nTotalFace);
    ReadUnsOrdinaryCellIndex(gridIndex, ordinaryCellIndex, nTotalCell);

    unsGrid->SetOrdinaryGridIndex(ordinaryGridIndex);
    unsGrid->SetOrdinaryNodeIndex(ordinaryNodeIndex);
    unsGrid->SetOrdinaryFaceIndex(ordinaryFaceIndex);
    unsGrid->SetOrdinaryCellIndex(ordinaryCellIndex);
}

void IO_ReadHDF5Grid::ReadStrOrdinaryGridIndex(int zoneIndex, int &ordinaryGridIndexOut)
{
    hid_t  gridGroup;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);

    bool dataExist = CheckDataExist(gridGroup, "PartitionInfo");
    if (!dataExist)
    {
        ordinaryGridIndexOut = -1;

        H5Gclose(gridGroup);
        return;
    }

    hid_t grpData;
    grpData = OpenGroup(gridGroup, "PartitionInfo");
    int ordinaryGridIndex = -1;
    ReadData(grpData, &ordinaryGridIndex, "OrdinaryGrid");
    ordinaryGridIndexOut = ordinaryGridIndex;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::ReadStrOrdinaryDimStart(int zoneIndex, int *ordinaryDimStartOut)
{
    hid_t  gridGroup;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);

    bool dataExist = CheckDataExist(gridGroup, "PartitionInfo");
    if (!dataExist)
    {
        ordinaryDimStartOut[0] = 0;
        ordinaryDimStartOut[1] = 0;
        ordinaryDimStartOut[2] = 0;
        H5Gclose(gridGroup);
        return;
    }

    hid_t grpData;
    grpData = OpenGroup(gridGroup, "PartitionInfo");
    int *ordinaryDimStartIndex = new int [3];
    ReadData(grpData, ordinaryDimStartIndex, "OrdinaryDimStart");
    ordinaryDimStartOut[0] = ordinaryDimStartIndex[0];
    ordinaryDimStartOut[1] = ordinaryDimStartIndex[1];
    ordinaryDimStartOut[2] = ordinaryDimStartIndex[2];
    delete [] ordinaryDimStartIndex;    ordinaryDimStartIndex = nullptr;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::ReadStrOrdinaryDimEnd(int zoneIndex, int *ordinaryDimEndOut)
{
    hid_t  gridGroup;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName   = oss.str();
    gridGroup = OpenGroup(fileID, grpName);

    bool dataExist = CheckDataExist(gridGroup, "PartitionInfo");
    if (!dataExist)
    {
        ordinaryDimEndOut[0] = 0;
        ordinaryDimEndOut[1] = 0;
        ordinaryDimEndOut[2] = 0;

        H5Gclose(gridGroup);
        return;
    }

    hid_t grpData;
    grpData = OpenGroup(gridGroup, "PartitionInfo");
    int *ordinaryDimEndIndex = new int[3];
    ReadData(grpData, ordinaryDimEndIndex, "OrdinaryDimEnd");
    ordinaryDimEndOut[0] = ordinaryDimEndIndex[0];
    ordinaryDimEndOut[1] = ordinaryDimEndIndex[1];
    ordinaryDimEndOut[2] = ordinaryDimEndIndex[2];
    delete [] ordinaryDimEndIndex;    ordinaryDimEndIndex = nullptr;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::ReadUnsOrdinaryGridIndex(int zoneIndex, int &ordinaryGridIndexOut)
{
    hid_t  gridGroup;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName = oss.str();
    gridGroup = OpenGroup(fileID, grpName);

    bool dataExist = CheckDataExist(gridGroup, "PartitionInfo");
    if (!dataExist)
    {
        ordinaryGridIndexOut = -1;
        H5Gclose(gridGroup);
        return;
    }

    hid_t grpData;
    grpData = OpenGroup(gridGroup, "PartitionInfo");
    int ordinaryGridIndex = -1;
    ReadData(grpData, &ordinaryGridIndex, "OrdinaryGrid");
    ordinaryGridIndexOut = ordinaryGridIndex;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::ReadUnsOrdinaryNodeIndex(int zoneIndex, int *ordinaryNodeIndexOut, int nTotalNode)
{
    hid_t  gridGroup;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName = oss.str();
    gridGroup = OpenGroup(fileID, grpName);

    bool dataExist = CheckDataExist(gridGroup, "PartitionInfo");
    if (!dataExist)
    {
        delete [] ordinaryNodeIndexOut;    ordinaryNodeIndexOut = nullptr;
        H5Gclose(gridGroup);
        return;
    }

    hid_t grpData;
    grpData = OpenGroup(gridGroup, "PartitionInfo");
    int *ordinaryNodeIndex = new int [nTotalNode];
    ReadData(grpData, ordinaryNodeIndex, "OrdinaryNodeIndex");
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        ordinaryNodeIndexOut[iNode] = ordinaryNodeIndex[iNode];
    }
    delete [] ordinaryNodeIndex;    ordinaryNodeIndex = nullptr;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::ReadUnsOrdinaryFaceIndex(int zoneIndex, int *ordinaryFaceIndexOut, int nTotalFace)
{
    hid_t  gridGroup;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName = oss.str();
    gridGroup = OpenGroup(fileID, grpName);

    bool dataExist = CheckDataExist(gridGroup, "PartitionInfo");
    if (!dataExist)
    {
        delete [] ordinaryFaceIndexOut;    ordinaryFaceIndexOut = nullptr;
        H5Gclose(gridGroup);
        return;
    }

    hid_t grpData;
    grpData = OpenGroup(gridGroup, "PartitionInfo");
    int *ordinaryFaceIndex = new int [nTotalFace];
    ReadData(grpData, ordinaryFaceIndex, "OrdinaryFaceIndex");
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        ordinaryFaceIndexOut[iFace] = ordinaryFaceIndex[iFace];
    }
    delete [] ordinaryFaceIndex;    ordinaryFaceIndex = nullptr;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::ReadUnsOrdinaryCellIndex(int zoneIndex, int *ordinaryCellIndexOut, int nTotalCell)
{
    hid_t  gridGroup;
    string grpName;
    ostringstream oss;
    oss << "Grid-" << zoneIndex;
    grpName = oss.str();
    gridGroup = OpenGroup(fileID, grpName);

    bool dataExist = CheckDataExist(gridGroup, "PartitionInfo");
    if (!dataExist)
    {
        delete [] ordinaryCellIndexOut;    ordinaryCellIndexOut = nullptr;
        H5Gclose(gridGroup);
        return;
    }

    hid_t grpData;
    grpData = OpenGroup(gridGroup, "PartitionInfo");
    int *ordinaryCellIndex = new int [nTotalCell];
    ReadData(grpData, ordinaryCellIndex, "OrdinaryCellIndex");
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        ordinaryCellIndexOut[iCell] = ordinaryCellIndex[iCell];
    }
    delete [] ordinaryCellIndex;    ordinaryCellIndex = nullptr;
    H5Gclose(grpData);
    H5Gclose(gridGroup);
}

void IO_ReadHDF5Grid::GetStructGridInfo(Grid *gridOut)
{
    int        gridIndex = gridOut->GetZoneID();
    StructGrid *strGrid  = StructGridCast(gridOut);

    //! Read grid size.
    int *iDimensions = new int[THREE_D];
    GetZoneSize(gridIndex, iDimensions);
    strGrid->SetNI(iDimensions[0]);
    strGrid->SetNJ(iDimensions[1]);
    strGrid->SetNK(iDimensions[2]);
    strGrid->SetBasicDimension();

    //! Read grid coordiate.
    int nTotalNode = iDimensions[0] * iDimensions[1] * iDimensions[2];
    RDouble **coordinates = new RDouble * [THREE_D];
    for (int iDim = 0; iDim < THREE_D; ++ iDim)
    {
        coordinates[iDim] = new RDouble[nTotalNode];
        GetZoneCoord(gridIndex, iDim, nTotalNode, coordinates[iDim]);
    }
    strGrid->SetX(coordinates[0]);
    strGrid->SetY(coordinates[1]);
    strGrid->SetZ(coordinates[2]);
    strGrid->RotateAboutAxis();
    strGrid->ComputeMinMaxBox();
    strGrid->SetArrayLayout();

    //! Read boundary information.
    int nBCRegion = 0;
    GetNBCRegion(gridIndex, nBCRegion);
    strGrid->CreateCompositeBCRegion(nBCRegion);
    StructBCSet *compositeBCRegion = strGrid->GetStructBCSet();
    GetBCName(gridIndex, compositeBCRegion);
    GetBCRegion(gridIndex, compositeBCRegion);
    strGrid->SetBCFaceInfo();

    //! Read interface information.
    int nIFace = 0;
    GetNIFace(gridIndex, nIFace);
    strGrid->SetNIFace(nIFace);
    InterfaceInfo *interfaceInfo = 0;
    if (nIFace > 0)
    {
        interfaceInfo = new InterfaceInfo(nIFace, strGrid);
        GetInterfaceInfo(gridIndex, interfaceInfo);
    }
    strGrid->SetInterfaceInfo(interfaceInfo);

    int ordinaryGridIndex = -1;
    int *ordinaryDimStartIndex = new int[3];
    int *ordinaryDimEndIndex   = new int[3];

    ReadStrOrdinaryGridIndex(gridIndex, ordinaryGridIndex);
    ReadStrOrdinaryDimStart(gridIndex, ordinaryDimStartIndex);
    ReadStrOrdinaryDimEnd(gridIndex, ordinaryDimEndIndex);

    strGrid->SetOrdinaryGridIndex(ordinaryGridIndex);
    strGrid->SetOrdinaryDimStartIndex(ordinaryDimStartIndex);
    strGrid->SetOrdinaryDimEndIndex(ordinaryDimEndIndex);

    delete [] ordinaryDimStartIndex;    ordinaryDimStartIndex = nullptr;
    delete [] ordinaryDimEndIndex;      ordinaryDimEndIndex   = nullptr;
    delete [] coordinates;              coordinates           = nullptr;
    delete [] iDimensions;              iDimensions           = nullptr;
}

void IO_ReadHDF5Grid::GetGridInfo(Grid *gridOut)
{
    int gridType = gridOut->Type();
    if (gridType == UNSTRUCTGRID)
    {
        GetUnstructGridInfo(gridOut);
    }
    else
    {
        GetStructGridInfo(gridOut);
    }
}

LIB_EXPORT void DumpHDF5Grid(const string &gridFileName, Grid **grids, int nBlocks)
{
    PrintToWindow("    Writing new HDF5 version grid to file ....\n");
    IO_HDF5Write writeGrid(gridFileName, grids, nBlocks);
    writeGrid.Run();
}

LIB_EXPORT bool ReadHDF5File(Region *region_in, PHVectorString1D &gridGroupNameList)
{
    bool HDF5FileRead = false;
    fstream file;
    string gridFileName = gridGroupNameList[0];
    OpenSerialFile(file, gridFileName, ios_base::in|ios_base::binary);
    int nZone;
    int sizeInt = sizeof(int);
    file.read((char *)&nZone, sizeInt);
    CloseFile(file);

    if (nZone <= static_cast <int> (1e8))
    {
        PrintToWindow("This grid file is old fts version!\n");
        WriteLogFile("This grid file is old fts version!\n");
        return HDF5FileRead;
    }
    else
    {
        PrintToWindow("This grid file is new HDF5 version!\n");
        WriteLogFile("This grid file is new HDF5 version!\n");
    }

    IO_HDF5Read readGrid(gridGroupNameList, region_in);
    HDF5FileRead = readGrid.Run();
    int isOldGrid = 0;
    GlobalDataBase::UpdateData("isOldGrid", &isOldGrid, PHINT, 1);
    return HDF5FileRead;
}

}