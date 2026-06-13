#include <cmath>
#include "Pre_Block_Struct.h"
#include "Pre_StructGridPartition.h"
#include "Glb_Dimension.h"
#include "IO_HDF5File.h"
#include "Pre_GridTranslator.h"
#include "IO_FileName.h"
#include "Pre_GridConversion.h"
#include "PHIO.h"
#include "TK_Log.h"
#include "Pre_HDF5File.h"

namespace PHSPACE
{
LIB_EXPORT Pre_StrGridPartition::Pre_StrGridPartition(int numberOfProcessors_in)
{
    this->numberOfProcessors = numberOfProcessors_in;

    geometricDimension = PHSPACE::GlobalDataBase::GetIntParaFromDB("ndim");

    boundaryContainer = new DataContainer();

    numberOfBasicBlocks = 0;

    LnkFileExist = false;
    partitionWallDist = false;

    cellIndexOfMark = new int[3];
    traceMark = PHSPACE::GlobalDataBase::GetIntParaFromDB("traceMark");
    if (traceMark)
    {
        blockIndexOfMark = PHSPACE::GlobalDataBase::GetIntParaFromDB("blockIndexOfMark");
        GlobalDataBase::GetData("cellIndexOfMark", cellIndexOfMark, PHINT, 3);
        if (geometricDimension == 2)
        {
            cellIndexOfMark[2] = 1;
        }
    }
}

LIB_EXPORT Pre_StrGridPartition::~Pre_StrGridPartition()
{
    delete boundaryContainer;
    delete [] cellIndexOfMark;    cellIndexOfMark = NULL;
    for (set<Pre_Block_Struct *>::iterator iter = generalBlockGroup.begin(); iter != generalBlockGroup.end(); ++ iter)
    {
        delete *iter;
    }
    return;
}

LIB_EXPORT void Pre_StrGridPartition::Run()
{
    string original_grid_file = "grid.fts";
    GlobalDataBase::GetData("original_grid_file", &original_grid_file, PHSTRING, 1);
    GlobalDataBase::UpdateData("gridfile", &original_grid_file, PHSTRING, 1);

    string partition_grid_file = "metis.fts";
    GlobalDataBase::GetData("partition_grid_file", &partition_grid_file, PHSTRING, 1);

    int isOverset = 0;
    GlobalDataBase::GetData("codeOfOversetGrid", &isOverset, PHINT, 1);
    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");

    Region *region = new Region();
    region->ReadGrid();

    numberOfOriginalBlocks = PHMPI::GetNumberofGlobalZones();
    OrdinaryGrid = new Grid *[numberOfOriginalBlocks];

    for (int iBlock = 0; iBlock < numberOfOriginalBlocks; iBlock ++)
    {
        OrdinaryGrid[iBlock] = GetGrid(iBlock, 0);
    }

    ProcessBCInfo();

    int isOldGrid = GlobalDataBase::GetIntParaFromDB("isOldGrid");
    if (isOldGrid)
    {
        ReadLnkFile(original_grid_file);
    }
    else
    {
        ReadLnkInfor(original_grid_file);
    }

    ReadBoundaryData();

    InitSplittingParameter();

    BinarySplitBlocks();

    SetBlockIndex();

    InitParititionGrids();

    OutputProcessorInformation();

    OutputPartitionedGrid();
    
    for (int iZone = 0; iZone < numberOfBasicBlocks; ++ iZone)
    {
        parititionGrids[iZone]->SetLnkInfo();
    }

    int dumpOldGrid = GlobalDataBase::GetIntParaFromDB("dumpOldGrid");
    int ndim = PHSPACE::GetDim();
    if (ndim == TWO_D || isOverset || (ndim == THREE_D && periodicType != NO_PERIODICITY))
    {
        CheckMeshMultigrid(parititionGrids, numberOfBasicBlocks);

        for (int iZone = 0; iZone < numberOfBasicBlocks; ++ iZone)
        {
            parititionGrids[iZone]->ProcessBCInfo();
        }

        ConstructGlobalInterfaceLink(parititionGrids, numberOfBasicBlocks);

        if (dumpOldGrid)
        {
            WriteLnkFile(partition_grid_file, parititionGrids, numberOfBasicBlocks);
        }

        int tasktype = GetTaskCode();
        int numberOfMultifile = PHSPACE::GlobalDataBase::GetIntParaFromDB("numberOfMultifile");

        if (tasktype == PHSPACE::PARTITION_GRID && numberOfMultifile > 1)
        {
            int maximumProcess = PHSPACE::GlobalDataBase::GetIntParaFromDB("maxproc");
            

            if (maximumProcess < numberOfMultifile)
            {
                ostringstream oss;
                oss << " Error: the numberOfMultifile's value is larger than maximumProcess's value " << endl;
                TK_Exit::ExceptionExit(oss.str());
                
            }
            else if (maximumProcess % numberOfMultifile != 0)
            {
                ostringstream oss;
                oss << " Error: the maximumProcess's value is not an integral multiple of numberOfMultifile's value " << endl;
                TK_Exit::ExceptionExit(oss.str());
            }

            int procNumberOfEachfile = maximumProcess / numberOfMultifile;
            for (int iFile = 0; iFile < numberOfMultifile; ++iFile)
            {
                vector<Grid *> tempGrids;
                tempGrids.resize(0);

                for (int iBlock = 0; iBlock < numberOfBasicBlocks; ++iBlock)
                {
                    int procIndex = parititionGrids[iBlock]->GetGridID()->GetIndex();
                    int *block_proc_inp = PHMPI::GetZoneProcessorID_INP();
                    if (block_proc_inp)
                    {
                        procIndex = block_proc_inp[iBlock];
                    }

                    if ((procIndex / procNumberOfEachfile) == iFile)
                    {
                        tempGrids.push_back(parititionGrids[iBlock]);
                    }
                    else
                    {
                        Grid *tempZero = nullptr;
                        tempGrids.push_back(tempZero);
                    }
                }

                string out_gfile_temp = AddSymbolToFileName(partition_grid_file, "_", iFile);
                DumpGrid(out_gfile_temp, &tempGrids[0], numberOfBasicBlocks);
            }
        }
        else
        {
            DumpGrid(partition_grid_file, parititionGrids, numberOfBasicBlocks);
        }
        
        string faceBCFileName = ChangeExtensionOfFileName(original_grid_file, "bc");
        if (SerialFileExist(faceBCFileName))
        {
            WriteBcFileFromFTS(partition_grid_file, parititionGrids, numberOfBasicBlocks);
        }
    }
    else
    {
        //! Output grid link relation before translate the grid.
        //! It's necessary to check if it work properly.
        if (dumpOldGrid)
        {
            WriteLnkFile(partition_grid_file, parititionGrids, numberOfBasicBlocks);
        }

        CheckMeshMultigrid(parititionGrids, numberOfBasicBlocks);
        GridTranslator *gridTranslator = new GridTranslator(parititionGrids, numberOfBasicBlocks);
        gridTranslator->Run();
 
        string faceBCFileName = ChangeExtensionOfFileName(original_grid_file, "bc");
        if (SerialFileExist(faceBCFileName)) WriteBcFileFromFTS(partition_grid_file, parititionGrids, numberOfBasicBlocks);
 
        delete gridTranslator;
    }

    WriteWalldist(partition_grid_file, parititionGrids);
    DumpMarkInfo();

    for (int i = 0; i < numberOfBasicBlocks; ++ i)
    {
        delete parititionGrids[i];
    }
    delete [] parititionGrids;    parititionGrids = NULL;

    delete region;
}

void Pre_StrGridPartition::WriteWalldist(const string &out_grid_file, Grid **grids_in)
{
    if (!partitionWallDist) return;
    string walldistfile = ChangeExtensionOfFileName(out_grid_file, "wdt");

    ActionKey *actkey = new ActionKey();
    actkey->filename = (walldistfile);

    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    string finalFileName = PHSPACE::AddSymbolToFileName(walldistfile, '_', currentProcessorIndex);
    actkey->filepos = H5Fcreate(finalFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    using namespace PHMPI;
    for (int iZone = 0; iZone < numberOfBasicBlocks; ++ iZone)
    {
        if (!grids_in[iZone]) continue;

        StructGrid *grid = StructGridCast(grids_in[iZone]);
        WriteWalldist(actkey, grid);
    }
    
    H5Fclose(actkey->filepos);
    delete actkey;    actkey = nullptr;
}

void Pre_StrGridPartition::WriteWalldist(ActionKey *actkey, StructGrid *grid)
{
    ostringstream oss;
    oss << "Grid-" << grid->GetZoneID();
    string dataName = oss.str();

    int nTotalCell = grid->GetNTotalCell();
    RDouble3D &walldist = * grid->GetWallDist();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    CreateAndWriteData(actkey->filepos, dataName, 1, nTotalCell, PHDOUBLE, &walldist(ist, jst, kst));
}

void Pre_StrGridPartition::ProcessBCInfo()
{
    for (int iBlock = 0; iBlock < numberOfOriginalBlocks; iBlock ++)
    {
        StructGrid *grid = StructGridCast(OrdinaryGrid[iBlock]);
        StructBCSet *structBCSet = grid->GetStructBCSet();

        int numberOfBoundaryFaces = structBCSet->GetnBCRegion();
        for (int iBCRegion = 0; iBCRegion < numberOfBoundaryFaces; iBCRegion ++)
        {
            StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);
           

            int *s_st = structBC->GetStartPoint();
            int *s_ed = structBC->GetEndPoint();
            int *s_lr3d = structBC->GetFaceDirectionIndex();

            for (int m = 0; m < geometricDimension; ++ m)
            {
                //! It may be true in two (or more) different directions after splited, now using the value already existed.
                if (structBC->GetFaceDirection() == m)
                {
                    if (s_st[m] > 1 || s_lr3d[m] == 1)
                    {
                        s_st[m] += 1;
                        s_ed[m] += 1;
                    }
                }
                else
                {
                    s_ed[m] += 1;
                }
            }
        }
    }
}

void Pre_StrGridPartition::SetBlockIndex()
{
    int zoneIndex = 0;
    bool blockIndexOfMarkChanged = false;

    for (int iProcessor = 0; iProcessor < numberOfProcessors; ++ iProcessor)
    {
        uint_t nZones = specialBlockTable[iProcessor].size();
        for (int j = 0; j < nZones; ++ j)
        {
            int oldZoneIndex = specialBlockTable[iProcessor][j]->GetZoneIndex();
            if (!blockIndexOfMarkChanged && blockIndexOfMark == oldZoneIndex)
            {
                blockIndexOfMark = zoneIndex;
                blockIndexOfMarkChanged = true;
            }

            specialBlockTable[iProcessor][j]->SetZoneIndex(zoneIndex);
            specialBlockTable[iProcessor][j]->SetProcessorIndex(iProcessor);
            ++ zoneIndex;
        }
    }

    return;
}

void Pre_StrGridPartition::OutputPartitionedGrid()
{
    ExtractBoundaryInformation();
    ExtractOriginalGridInformation();
    OutputPartitionedGridFile();
    SetWallDist();
}

void Pre_StrGridPartition::SetWallDist()
{
    boundaryContainer->MoveToBegin();

    int numberOfBasicBlocks = 0;
    PHSPACE::PHRead(boundaryContainer, numberOfBasicBlocks);

    for (int i = 0; i < numberOfBasicBlocks; ++ i)
    {
        vector< int > nodeDimension(3);
        PHSPACE::PHRead(boundaryContainer, nodeDimension, 3);
    }

    int iBasic;
    vector< int > nst(3), ned(3);

    for (int i = 0; i < numberOfBasicBlocks; ++ i)
    {
        StructGrid *grid = StructGridCast(parititionGrids[i]);

        PHSPACE::PHRead(boundaryContainer, iBasic);

        PHSPACE::PHRead(boundaryContainer, nst[0]);
        PHSPACE::PHRead(boundaryContainer, ned[0]);
        PHSPACE::PHRead(boundaryContainer, nst[1]);
        PHSPACE::PHRead(boundaryContainer, ned[1]);
        PHSPACE::PHRead(boundaryContainer, nst[2]);
        PHSPACE::PHRead(boundaryContainer, ned[2]);

        StructGrid *Ordinary_grid = StructGridCast(OrdinaryGrid[iBasic]);

        WriteGridWallDist(grid, Ordinary_grid, nst);

        grid->SetOrdinaryGridIndex(iBasic);
        grid->SetOrdinaryDimStartIndex(&nst[0]);
        grid->SetOrdinaryDimEndIndex(&ned[0]);
    }

//    if (partitionWallDist)
//    {
        cout << "\nWrite ParititionGrids' WallDist completed ... \n" << endl;
//    }
}

void Pre_StrGridPartition::WriteGridWallDist(StructGrid *grid, StructGrid *Ordinary_grid, vector< int > &nst)
{
    if (!Ordinary_grid->GetWallDist()) return;

    grid->InitWallDist();

    RDouble3D &walldist = * grid->GetWallDist();
    RDouble3D &Ordinary_walldist = * Ordinary_grid->GetWallDist();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int kk = nst[2];
    for (int k = kst; k <= ked; ++ k)
    {
        int jj = nst[1];
        for (int j = jst; j <= jed; ++ j)
        {
            int ii = nst[0];
            for (int i = ist; i <= ied; ++ i)
            {
                walldist(i, j, k) = Ordinary_walldist(ii, jj, kk);
                ii ++;
            }
            jj ++;
        }
        kk ++;
    }

    partitionWallDist = true;
}

void Pre_StrGridPartition::OutputPartitionedGridFile()
{
    boundaryContainer->MoveToBegin();

    int numberOfBasicBlocks = 0;
    PHSPACE::PHRead(boundaryContainer, numberOfBasicBlocks);

    for (int i = 0; i < numberOfBasicBlocks; ++ i)
    {
        vector< int > nodeDimension(3);
        PHSPACE::PHRead(boundaryContainer, nodeDimension, 3);
    }

    cout << "Writing ParititionGrids' Coordinate" << endl;

    int iProcessors = 0;
    uint_t BlocksOfiProcessors = specialBlockTable[iProcessors].size();
    int count = 0;
    int myid = PHMPI::GetCurrentProcessorID();

    int iBasic;
    vector< int > nst(3), ned(3);
    for (int i = 0; i < numberOfBasicBlocks; ++ i)
    {
        StructGrid *grid = StructGridCast(parititionGrids[i]);

        PHSPACE::PHRead(boundaryContainer, iBasic);

        PHSPACE::PHRead(boundaryContainer, nst[0]);
        PHSPACE::PHRead(boundaryContainer, ned[0]);
        PHSPACE::PHRead(boundaryContainer, nst[1]);
        PHSPACE::PHRead(boundaryContainer, ned[1]);
        PHSPACE::PHRead(boundaryContainer, nst[2]);
        PHSPACE::PHRead(boundaryContainer, ned[2]);
        
        Pre_Block_Struct *simpleBlock = referencedBlockContainer[iBasic];

        WriteGridCoordinate(grid, simpleBlock, nst, ned);

        ++ count;
        if (myid == PHMPI::GetServerProcessorID() && count > 0 && count == BlocksOfiProcessors)
        {
            ++ iProcessors;
            ProgressMonitorToLogFile(iProcessors, numberOfProcessors, "grid partition reconstructing");
            ProgressMonitorToWindows(iProcessors, numberOfProcessors, "grid partition reconstructing");

            if (iProcessors == numberOfProcessors) continue;

            BlocksOfiProcessors = specialBlockTable[iProcessors].size();
            count = 0;
        }
    }
}

void Pre_StrGridPartition::WriteGridCoordinate(StructGrid *grid, Pre_Block_Struct *simpleBlock, vector< int > &nst, vector< int > &ned)
{
    vector< RDouble * > s(3, NULL);

    s[0] = simpleBlock->GetX();
    s[1] = simpleBlock->GetY();
    s[2] = simpleBlock->GetZ();

    int nx = ned[0] - nst[0] + 1;
    int ny = ned[1] - nst[1] + 1;
    int nz = ned[2] - nst[2] + 1;
    
    grid->SetNI(nx);
    grid->SetNJ(ny);
    grid->SetNK(nz);

    Range IC(1, nx-1);
    Range JC(1, ny-1);
    Range KC(1, nz-1);

    if (grid->GetDim() == TWO_D) KC.setRange(1, 1);

    int nTotalNode, nTotalFace, nTotalCell;

    nTotalNode = nx * ny * nz;
    nTotalCell = IC.length() * JC.length() * KC.length();

    nTotalFace = (ny - 1) * (nz - 1) * nx + (nz - 1) * (nx - 1) * ny + (nx - 1) * (ny - 1) * nz;

    if (grid->GetDim() == TWO_D)
    {
        nTotalFace = (ny - 1) * nx + (nx - 1) * ny;
    }

    grid->SetNTotalNode(nTotalNode);
    grid->SetNTotalFace(nTotalFace);
    grid->SetNTotalCell(nTotalCell);

    int ni = simpleBlock->GetNI();
    int nj = simpleBlock->GetNJ();

    RDouble *x = new RDouble [nTotalNode];
    RDouble *y = new RDouble [nTotalNode];
    RDouble *z = new RDouble [nTotalNode];
    grid->SetX(x);
    grid->SetY(y);
    grid->SetZ(z);
    grid->SetArrayLayout();

    RDouble3D &xx = * grid->GetStructX();
    RDouble3D &yy = * grid->GetStructY();
    RDouble3D &zz = * grid->GetStructZ();

    Range I(1, nx);
    Range J(1, ny);
    Range K(1, nz);
    if (nz == 1) K.setRange(1, 1);
    
    int ii, jj, kk;
    for (int k = nst[2]-1; k < ned[2]; ++ k)
    {
        kk = k - (nst[2] - 1) + 1;
        for (int j = nst[1]-1; j < ned[1]; ++ j)
        {
            jj = j - (nst[1] - 1) + 1;
            for (int i = nst[0]-1; i < ned[0]; ++ i)
            {
                ii = i - (nst[0] - 1) + 1;
                int id = ni * nj * k + ni * j + i;
                xx(ii, jj, kk) = s[0][id];
                yy(ii, jj, kk) = s[1][id];
                zz(ii, jj, kk) = s[2][id];
            }
        }
    }

    grid->ComputeMinMaxBox();
}

void Pre_StrGridPartition::ExtractBoundaryInformation()
{
    boundaryContainer->MoveToBegin();
    PHSPACE::PHWrite(boundaryContainer, numberOfBasicBlocks);

    for (int iProcessor = 0; iProcessor < numberOfProcessors; ++ iProcessor)
    {
        uint_t gp_size = specialBlockTable[iProcessor].size();
        for (int j = 0; j < gp_size; ++ j)
        {
            vector< int > &nodeDimension = specialBlockTable[iProcessor][j]->GetNodeDimension();
            PHSPACE::PHWrite(boundaryContainer, nodeDimension, 3);
        }
    }

    int originalZoneIndex, ist, ied, jst, jed, kst, ked;
    vector< int > ref;
    Pre_Block_Struct *simpleBlock = 0;

    for (int iProcessor = 0; iProcessor < numberOfProcessors; ++ iProcessor)
    {
        uint_t gp_size = specialBlockTable[iProcessor].size();
        for (int j = 0; j < gp_size; ++ j)
        {
            vector< int > &nodeDimension = specialBlockTable[iProcessor][j]->GetNodeDimension();
            vector< int >(3, 0).swap(ref);

            specialBlockTable[iProcessor][j]->GetRootInformation(simpleBlock, ref);

            originalZoneIndex = simpleBlock->GetOriginalZoneIndex();
            PHSPACE::PHWrite(boundaryContainer, originalZoneIndex);

            ist = ref[0] + 1;
            ied = ref[0] + nodeDimension[0];
            jst = ref[1] + 1;
            jed = ref[1] + nodeDimension[1];
            kst = ref[2] + 1;
            ked = ref[2] + nodeDimension[2];

            PHSPACE::PHWrite(boundaryContainer, ist);
            PHSPACE::PHWrite(boundaryContainer, ied);
            PHSPACE::PHWrite(boundaryContainer, jst);
            PHSPACE::PHWrite(boundaryContainer, jed);
            PHSPACE::PHWrite(boundaryContainer, kst);
            PHSPACE::PHWrite(boundaryContainer, ked);
        }
    }
    return;
}

void Pre_StrGridPartition::ExtractOriginalGridInformation()
{
    for (int iBlock = 0; iBlock < numberOfOriginalBlocks; ++ iBlock)
    {
        Pre_Block_Struct *originalBlock = referencedBlockContainer[iBlock];
        int ni, nj, nk;

        StructGrid *grid = StructGridCast(OrdinaryGrid[iBlock]);

        ni = grid->GetNI();
        nj = grid->GetNJ();
        nk = grid->GetNK();
        originalBlock->SetNI(ni);
        originalBlock->SetNJ(nj);
        originalBlock->SetNK(nk);

        originalBlock->GenerateCoordinateSpace();
    }

    for (int iBlock = 0; iBlock < numberOfOriginalBlocks; ++ iBlock)
    {
        Pre_Block_Struct *originalBlock = referencedBlockContainer[iBlock];
        
        StructGrid *grid = StructGridCast(OrdinaryGrid[iBlock]);

        int nTotalNode = grid->GetNTotalNode();

        RDouble *x = grid->GetX();
        RDouble *y = grid->GetY();
        RDouble *z = grid->GetZ();
        
        RDouble *xx = originalBlock->GetX();
        RDouble *yy = originalBlock->GetY();
        RDouble *zz = originalBlock->GetZ();

        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            xx[iNode] = x[iNode];
            yy[iNode] = y[iNode];
            zz[iNode] = z[iNode];
        }
    }

    return;
}

void Pre_StrGridPartition::OutputProcessorInformation()
{
    string cs;
    std::ostringstream oss;

    int *block_proc_inp = PHMPI::GetZoneProcessorID_INP();

    oss << "\n";
    oss << "      ----- Grid Partition -----\n"
        << " Number of final partitions: " << numberOfProcessors << "\n";
        //<< " Number of processors = " << numberOfProcessors << "\n";

    int procwidth = GetIntegerDigitWidth(numberOfProcessors);

    int mincellsize = 0;
    for (std::size_t j = 0; j < specialBlockTable[0].size(); ++ j)
    {
        mincellsize += specialBlockTable[0][j]->GetNumberOfCells();
    }
    int maxcellsize = mincellsize;

    int totalcellsize = 0;

    int numberOfBlocks = 0;
    for (int i = 0; i < numberOfProcessors; ++ i)
    {
        int count = 0;
        int localcellsize = 0;

        for (std::size_t j = 0; j < specialBlockTable[i].size(); ++ j)
        {
            cout << " proc = " << i << " maxproc = " << numberOfProcessors;
            cout << " block " << numberOfBlocks << " size = " << specialBlockTable[i][j]->GetNumberOfNodes() << "\n";
            count += specialBlockTable[i][j]->GetNumberOfNodes();
            block_proc_inp[numberOfBlocks] = i;
            ++ numberOfBlocks;
            localcellsize += specialBlockTable[i][j]->GetNumberOfCells();
        }
        cout << " total size = " << count << "\n";

        oss << " Processor " << setw(procwidth) << i << ": block = " << specialBlockTable[i].size() 
            << ", size = " << count << ", cell = " << localcellsize << "\n";

        totalcellsize += localcellsize;

        if (mincellsize > localcellsize)
        {
            mincellsize = localcellsize;
        }

        if (maxcellsize < localcellsize)
        {
            maxcellsize = localcellsize;
        }
    }

    oss.setf(ios::fixed);
    oss.precision(2);
    oss << " Total   Cell Size = " << totalcellsize << "\n" 
        << " Average Cell Size = " << static_cast<RDouble>(totalcellsize) / numberOfProcessors << "\n"
        << " Minimum Cell Size = " << mincellsize << "(" 
        << 100. * (static_cast<RDouble>(mincellsize * numberOfProcessors) / totalcellsize - 1.) << "%)\n"
        << " Maximum Cell Size = " << maxcellsize << "("
        << 100. * (static_cast<RDouble>(maxcellsize * numberOfProcessors) / totalcellsize - 1.) << "%)\n";
    (void)oss.flags();

    cs = oss.str();
    cout << cs << endl;
    WriteLogFile(cs);

    int count = 0;
    for (int i = 0; i < numberOfProcessors; ++ i)
    {
        for (std::size_t j = 0; j < specialBlockTable[i].size(); ++ j)
        {
            StructGrid *grid = StructGridCast(parititionGrids[count]);
            count ++;

            OutputBoundaryInformation(specialBlockTable[i][j], grid, i);
        }
    }
    return;
}

void Pre_StrGridPartition::OutputBoundaryInformation(Pre_Block_Struct *Pre_Block_Struct_in, StructGrid *grid, int iZone)
{
    vector< Pre_BcFace_Struct * > leafBoundaryFace;
    vector< Pre_BcFace_Struct * > boundaryFaceList = Pre_Block_Struct_in->GetBoundaryFaceList();

    uint_t numberOfBoundaryFaces = boundaryFaceList.size();
    for (int i = 0; i < numberOfBoundaryFaces; ++ i)
    {
        vector< Pre_BcFace_Struct * > cl = boundaryFaceList[i]->GetLeaves();
        uint_t numberOfLeaves = cl.size();
        for (int j = 0; j < numberOfLeaves; ++ j)
        {
            leafBoundaryFace.push_back(cl[j]);
        }
    }

    int numberOfLeafBoundaryFaces = static_cast<int>(leafBoundaryFace.size());
    grid->CreateCompositeBCRegion(numberOfLeafBoundaryFaces);
    StructBCSet *structBCSet = grid->GetStructBCSet();

    for (int i = 0; i < numberOfLeafBoundaryFaces; ++ i)
    {
        StructBC *bcregion = new StructBC(iZone, i);
        structBCSet->SetBCRegion(i, bcregion);
        leafBoundaryFace[i]->OutputBoundaryInformation(bcregion, geometricDimension);
    }

    return;
}

void Pre_StrGridPartition::InitParititionGrids()
{
    for (int iProcessor = 0; iProcessor < numberOfProcessors; ++ iProcessor)
    {
        numberOfBasicBlocks += static_cast<int>(specialBlockTable[iProcessor].size());
    }

    PHMPI::CreatetZoneProcessorID_INP(numberOfBasicBlocks);

    parititionGrids = new Grid * [numberOfBasicBlocks];
    for (int iZone = 0; iZone < numberOfBasicBlocks; ++ iZone)
    {
        GridID *index = new GridID(iZone);
        Grid *grid;
        grid = CreateGridGeneral(STRUCTGRID, index, 0, geometricDimension);
        parititionGrids[iZone] = grid;
    }
}

void Pre_StrGridPartition::BinarySplitBlocks()
{
    const RDouble nmin = 5 * 5 * 2;

    while (! unallottedBlockGroup.empty())
    {
        cout << "unallottedBlockGroup.size = " << unallottedBlockGroup.size() << "\n";
        int first = GetMinimumElementIndex(reserves);

        Pre_Block_Struct *largestBlock = 0;
        RDouble largestNumberOfCells = 0.0;
        ProbeLargestBlock(largestBlock, largestNumberOfCells);

        RDouble tau     = reserves[first] + largestNumberOfCells;
        RDouble nrest   = tau - idealTime;
        RDouble ntarget = largestNumberOfCells - nrest;

        if (tau > idealTime * (1.0 + epsilon) && nrest > nmin && ntarget > nmin)
        {
            RDouble ratio = ntarget / largestNumberOfCells;

            if (abs(ratio - 0.5) < 0.25)
            {
                SplitSingleBlock(largestBlock, ntarget, first, true);
            }
            else
            {
                RDouble nhalf = largestNumberOfCells / 2.0;
                SplitSingleBlock(largestBlock, nhalf, first, false);
            }
        }
        else
        {
            AssignSingleBlock(largestBlock, first, largestNumberOfCells);
        }
    }

    if (unallottedBlockGroup.empty())
    {
        cout << "All BlockGroup have been allotted!\n";
    }
    return;
}

void Pre_StrGridPartition::SplitSingleBlock(Pre_Block_Struct *simpleBlock, RDouble ntarget, int ip, bool allotted)
{
    vector< int > leftDimension(3), rightDimension(3), leftOriginalIndex(3), rightOriginalIndex(3);
    simpleBlock->Split(ntarget, leftDimension, rightDimension, leftOriginalIndex, rightOriginalIndex);

    Pre_Block_Struct *leftBlock = new Pre_Block_Struct();
    leftBlock->SetNodeDimension(leftDimension);
    leftBlock->SetOriginalIndex(leftOriginalIndex);

    Pre_Block_Struct *rightBlock = new Pre_Block_Struct();
    rightBlock->SetNodeDimension(rightDimension);
    rightBlock->SetOriginalIndex(rightOriginalIndex);

    leftBlock->SetZoneIndex(numberOfGeneralBlocks);
    numberOfGeneralBlocks += 1;

    rightBlock->SetZoneIndex(numberOfGeneralBlocks);
    numberOfGeneralBlocks += 1;

    leftBlock->SetParentBlock(simpleBlock);
    rightBlock->SetParentBlock(simpleBlock);

    leftBlock->CreateBoundaryFace(geometricDimension);
    rightBlock->CreateBoundaryFace(geometricDimension);

    TraceMark(leftBlock, rightBlock, simpleBlock);

    CreateInterfaceBoundaryFace(leftBlock, rightBlock);

    unallottedBlockGroup.erase(simpleBlock);
    unallottedBlockGroup.insert(leftBlock);
    unallottedBlockGroup.insert(rightBlock);

    generalBlockGroup.insert(leftBlock);
    generalBlockGroup.insert(rightBlock);

    if (allotted)
    {
        AssignSingleBlock(leftBlock, ip, ntarget);
    }

    return;
}

void Pre_StrGridPartition::TraceMark(Pre_Block_Struct *leftBlock, Pre_Block_Struct *rightBlock, Pre_Block_Struct *parentBlock)
{
    if (!traceMark)
    {
        return;
    }
    else if (parentBlock->GetZoneIndex() != blockIndexOfMark)
    {
        return;
    }

    vector<int> leftDimension = leftBlock->GetNodeDimension();
    vector<int> rightOriginalIndex = rightBlock->GetOriginalIndex();

    int axisSelect = 0;
    for (int iDimension = 0; iDimension < geometricDimension; ++ iDimension)
    {
        if (rightOriginalIndex[iDimension])
        {
            axisSelect = iDimension;
        }
    }

    if (cellIndexOfMark[axisSelect] < leftDimension[axisSelect])
    {
        blockIndexOfMark = leftBlock->GetZoneIndex();
    }
    else
    {
        blockIndexOfMark = rightBlock->GetZoneIndex();
        cellIndexOfMark[axisSelect] = cellIndexOfMark[axisSelect] - leftDimension[axisSelect] + 1;
    }
}

void Pre_StrGridPartition::DumpMarkInfo()
{
    if (!traceMark)
    {
        return;
    }

    std::ostringstream oss;
    oss << endl;
    oss << "blockIndexOfMark = " << blockIndexOfMark << endl;
    oss << "cellIndexOfMark =";

    for (int iDimension = 0; iDimension < geometricDimension; ++ iDimension)
    {
        oss << "  " << cellIndexOfMark[iDimension];
    }
    oss << endl;

    WriteLogFile(oss);
}

void Pre_StrGridPartition::AssignSingleBlock(Pre_Block_Struct *simpleBlock, int ip, RDouble numberOfCells)
{
    reserves[ip] += numberOfCells;
    unallottedBlockGroup.erase(simpleBlock);
    specialBlockTable[ip].push_back(simpleBlock);

    return;
}

void Pre_StrGridPartition::CreateInterfaceBoundaryFace(Pre_Block_Struct *leftBlock, Pre_Block_Struct *rightBlock)
{
    vector< int > l_st(3), l_ed(3), r_st(3), r_ed(3);
    leftBlock->GetAbsoluteCoordinate(l_st, l_ed);
    rightBlock->GetAbsoluteCoordinate(r_st, r_ed);

    vector< int > st(3), ed(3);

    if (GetInterfaceCoordinate(l_st, l_ed, r_st, r_ed, st, ed))
    {
        CreateInterfaceBoundaryFace(leftBlock, rightBlock, st, ed);
    }

    return;
}

void Pre_StrGridPartition::CreateInterfaceBoundaryFace(Pre_Block_Struct *leftBlock, Pre_Block_Struct *rightBlock, vector< int > &st, vector< int > &ed)
{
    //! Note that Normalize is not independent ,it is ordering!
    vector< int > l_st(3), l_ed(3), r_st(3), r_ed(3);
    leftBlock->GetLocalCoordinate(st, ed, l_st, l_ed);
    rightBlock->GetLocalCoordinate(st, ed, r_st, r_ed);

    Pre_BcFace_Struct *lbc = new Pre_BcFace_Struct();
    Pre_BcFace_Struct *rbc = new Pre_BcFace_Struct();

    lbc->SetSimpleBlock(leftBlock);
    rbc->SetSimpleBlock(rightBlock);

    MarkPatchFace(l_st, l_ed, geometricDimension);
    MarkPatchFace(r_st, r_ed, geometricDimension);

    lbc->SetRegion(l_st, l_ed);
    lbc->SetBoundaryType(-1);
    lbc->SetNext(new Pre_Patch_Struct());
    lbc->GetNext()->SetSimpleBlock(rightBlock);
    lbc->GetNext()->SetRegion(r_st, r_ed);

    rbc->SetRegion(r_st, r_ed);
    rbc->SetBoundaryType(-1);
    rbc->SetNext(new Pre_Patch_Struct());
    rbc->GetNext()->SetSimpleBlock(leftBlock);
    rbc->GetNext()->SetRegion(l_st, l_ed);

    lbc->Normalize(geometricDimension);
    rbc->Normalize(geometricDimension);

    leftBlock->AddBoundaryFace(lbc);
    rightBlock->AddBoundaryFace(rbc);

    return;
}

void Pre_StrGridPartition::MarkPatchFace(vector< int > &st, vector< int > &ed, int geometricDimension)
{
    int nsw  [] = { 0, 1, 2, 0, 1 };
    int nsw2d[] = { 0, 1, 0, 1 };
    int nsurf = GetSurfaceLabel(st, ed);

    int m = nsw[nsurf+1];

    if (geometricDimension == 2) m = nsw2d[nsurf+1];

    st[m] = - abs(st[m]);
    ed[m] = - abs(ed[m]);

    return;
}

bool Pre_StrGridPartition::GetInterfaceCoordinate(vector< int > &l_st, vector< int > &l_ed, vector< int > &r_st, vector< int > &r_ed, vector< int > &st, vector< int > &ed)
{
    vector< int > ll_st(3), ll_ed(3), rr_st(3), rr_ed(3);
    for (int m = 0; m < 6; ++ m)
    {
        GetIJKRegion(l_st, l_ed, ll_st, ll_ed, m);
        for (int n = 0; n < 6; ++ n)
        {
            GetIJKRegion(r_st, r_ed, rr_st, rr_ed, n);
            if (IsSameDomain(ll_st, ll_ed, rr_st, rr_ed))
            {
                st = ll_st;
                ed = ll_ed;
                return true;
            }
        }
    }

    return false;
}

bool Pre_StrGridPartition::IsSameDomain(vector< int > &l_st, vector< int > &l_ed, vector< int > &r_st, vector< int > &r_ed)
{
    for (int m = 0; m < 3; ++ m)
    {
        if (l_st[m] != r_st[m] || l_ed[m] != r_ed[m])
        {
            return false;
        }
    }
    return true;
}

void Pre_StrGridPartition::ProbeLargestBlock(Pre_Block_Struct *&simpleBlock, RDouble &numberOfCells)
{
    numberOfCells = 0.0;

    for (set< Pre_Block_Struct * >::iterator iter = unallottedBlockGroup.begin(); iter != unallottedBlockGroup.end(); ++ iter)
    {
        RDouble nsize = (*iter)->GetNumberOfCells();
        if (numberOfCells < nsize)
        {
            numberOfCells = nsize;
            simpleBlock   = (*iter);
        }
    }
    return;
}

int Pre_StrGridPartition::GetMinimumElementIndex(vector<RDouble> &x)
{
    RDouble value = x[0];
    int index = 0;
    uint_t numberOfElements = x.size();

    for (int iElement = 1; iElement < numberOfElements; ++ iElement)
    {
        if (x[iElement] < value)
        {
            index = iElement;
            value = x[iElement];
        }
    }

    return index;
}

void Pre_StrGridPartition::InitSplittingParameter()
{
    numberOfGeneralBlocks = numberOfOriginalBlocks;
    delta = 1.e-5;
    reserves.resize(numberOfProcessors, delta);

    RDouble sumOfOriginalBlockCells = 0.0;
    for (int iBlock = 0; iBlock < numberOfOriginalBlocks; ++ iBlock)
    {
        sumOfOriginalBlockCells += referencedBlockContainer[iBlock]->GetNumberOfCells();

        unallottedBlockGroup.insert(referencedBlockContainer[iBlock]);

        generalBlockGroup.insert(referencedBlockContainer[iBlock]);
    }

    idealTime = sumOfOriginalBlockCells / numberOfProcessors;
    epsilon = 1.e-2;

    specialBlockTable.resize(numberOfProcessors);
}

void Pre_StrGridPartition::ReadLnkFile(const string &fileName_in)
{
    string LnkFileName = ChangeExtensionOfFileName(fileName_in, "link");
    if (!SerialFileExist(LnkFileName))
    {
        TK_Exit::ExceptionExit("There is no link file!");
        return;
    }
    LnkFileExist = true;

    fstream file;
    OpenSerialFile(file, LnkFileName, ios_base::in|ios_base::binary);

    DataContainer *LnkData = new DataContainer();
    ReadFile(file, LnkData);

    VirtualFile *virtualFile = new VirtualFile(LnkData);
    virtualFile->BeginReadWork();

    Grid **OrdinaryGrid = GetOrdinaryGrid();
    for (int iBlock = 0; iBlock < numberOfOriginalBlocks; iBlock ++)
    {
        StructGrid *grid = StructGridCast(OrdinaryGrid[iBlock]);

        StructBCSet *structBCSet = grid->GetStructBCSet();

        int nBCRegion = structBCSet->GetnBCRegion();
        for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion ++)
        {
            StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);

            int bctype = structBC->GetBCType();
            if ((!IsInterface(bctype))) continue;

            int imin, imax, jmin, jmax, kmin, kmax, nbt;

            PHRead(virtualFile, imin);
            PHRead(virtualFile, imax);
            PHRead(virtualFile, jmin);
            PHRead(virtualFile, jmax);

            if (geometricDimension == THREE_D)
            {
                PHRead(virtualFile, kmin);
                PHRead(virtualFile, kmax);
            }
            else
            {
                kmin = 1;
                kmax = 1;
            }

            structBC->SetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);

            PHRead(virtualFile, imin);
            PHRead(virtualFile, imax);
            PHRead(virtualFile, jmin);
            PHRead(virtualFile, jmax);

            if (geometricDimension == THREE_D)
            {
                PHRead(virtualFile, kmin);
                PHRead(virtualFile, kmax);
            }
            else
            {
                kmin = 1;
                kmax = 1;
            }

            PHRead(virtualFile, nbt);

            structBC->SetTargetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);
            structBC->SetTargetRegionBlock(nbt - 1);
            if (GetDim() == THREE_D)
            {
                structBC->ComputeRelativeParameters();
            }
        }
    }

    virtualFile->EndReadWork();

    delete virtualFile;
    delete LnkData;

    file.close();
    file.clear();
}

void Pre_StrGridPartition::ReadLnkInfor(const string &fileName_in)
{
    hid_t file;

    int fileID = PHMPI::GetFileIndexofCurrentProcessor();
    string actualFileName = AddSymbolToFileName(fileName_in, "_", fileID);
    file = OpenHDF5File(actualFileName);

    int nZones = PHMPI::GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        OrdinaryGrid[iZone] = GetGrid(iZone, 0);

        StructGrid *grid = StructGridCast(OrdinaryGrid[iZone]);
        StructBCSet *structBCSet = grid->GetStructBCSet();

        int nBCLnkTotal = 0;
        int nBCRegion = structBCSet->GetnBCRegion();
        for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion ++)
        {
            StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);

            int bctype = structBC->GetBCType();
            if (bctype == PHENGLEI::INTERFACE)
            {
                nBCLnkTotal++;
            }
        }

        if (!nBCLnkTotal)
        {
            continue;
        }

        hid_t grpGrid, grpData;
        string grpName;

        ostringstream oss;
        oss << "Grid-" << iZone;
        grpName = oss.str();
        grpGrid = OpenGroup(file, grpName);

        grpData = OpenGroup(grpGrid, "LnkInfo");
        int **LnkInfor = ReadIntTwoRow(grpData, "LnkInfor");

        int iMin, iMax, jMin, jMax, kMin, kMax, nbt;
        int iLnkinfo = 0;
        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            int index = 0;
            StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);

            int bctype = bcregion->GetBCType();
            if (!IsInterface(bctype)) continue;

            iMin = LnkInfor[iLnkinfo][index++];
            iMax = LnkInfor[iLnkinfo][index++];
            jMin = LnkInfor[iLnkinfo][index++];
            jMax = LnkInfor[iLnkinfo][index++];
            kMin = LnkInfor[iLnkinfo][index++];
            kMax = LnkInfor[iLnkinfo][index++];
            bcregion->SetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);

            iLnkinfo ++;
            index = 0;

            iMin = LnkInfor[iLnkinfo][index++];
            iMax = LnkInfor[iLnkinfo][index++];
            jMin = LnkInfor[iLnkinfo][index++];
            jMax = LnkInfor[iLnkinfo][index++];
            kMin = LnkInfor[iLnkinfo][index++];
            kMax = LnkInfor[iLnkinfo][index++];

            nbt = LnkInfor[iLnkinfo][index++];
            bcregion->SetTargetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);
            bcregion->SetTargetRegionBlock(nbt - 1);
            if (GetDim() == THREE_D)
            {
                bcregion->ComputeRelativeParameters();
            }

            iLnkinfo ++;
        }

        DelPointer2(LnkInfor);

        H5Gclose(grpData);
        H5Gclose(grpGrid);
    }

    H5Fclose(file);
}

void Pre_StrGridPartition::ReadBoundaryData()
{
    int numberOfMultigrid = PHSPACE::GlobalDataBase::GetIntParaFromDB("numberOfMultigrid");

    cout << "numberOfMultigrid = " << numberOfMultigrid << endl;

    Pre_Block_Struct::CalculateNumberOfUnitCells(numberOfMultigrid);

    cout << "numberOfOriginalBlocks = " << numberOfOriginalBlocks << "\n";
    referencedBlockContainer.resize(numberOfOriginalBlocks);

    int numberOfUnitCells = Pre_Block_Struct::GetNumberOfUnitCells();

    Grid **OrdinaryGrid = GetOrdinaryGrid();

    for (int iBlock = 0; iBlock < numberOfOriginalBlocks; ++ iBlock)
    {
        referencedBlockContainer[iBlock] = new Pre_Block_Struct();
        referencedBlockContainer[iBlock]->SetOriginalZoneIndex(iBlock);
        referencedBlockContainer[iBlock]->SetZoneIndex(iBlock);
        referencedBlockContainer[iBlock]->SetProcessorIndex(0);
    }

    for (int iBlock = 0; iBlock < numberOfOriginalBlocks; ++ iBlock)
    {
        StructGrid *grid = StructGridCast(OrdinaryGrid[iBlock]);

        vector< int > nodeDimension(3);
        nodeDimension[0] = grid->GetNI();
        nodeDimension[1] = grid->GetNJ();
        nodeDimension[2] = grid->GetNK();

        Pre_Block_Struct *simpleBlock = referencedBlockContainer[iBlock];
        simpleBlock->SetNodeDimension(nodeDimension);
        simpleBlock->SetNI(nodeDimension[0]);
        simpleBlock->SetNJ(nodeDimension[1]);
        simpleBlock->SetNK(nodeDimension[2]);

        StructBCSet *structBCSet = grid->GetStructBCSet();
        int numberOfBoundaryFaces = structBCSet->GetnBCRegion();
        simpleBlock->GenerateBoundaryFaceList(numberOfBoundaryFaces);

        for (int iBCRegion = 0; iBCRegion < numberOfBoundaryFaces; iBCRegion ++)
        {
            StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);

            vector< int > nst(3), ned(3);
            int bctype;
            string bcname;

            for (int i = 0; i < 3; i ++)
            {
                nst[i] = structBC->GetStartPoint(i);
                ned[i] = structBC->GetEndPoint(i);
            }

            CheckMultigrid(nst, ned, numberOfUnitCells);

            bctype = structBC->GetBCType();
            bcname = structBC->GetBCName();
            int bcdirection = structBC->GetFaceDirection();

            Pre_BcFace_Struct *boundaryFace = new Pre_BcFace_Struct();
            simpleBlock->SetBoundaryFace(iBCRegion, boundaryFace);
            boundaryFace->SetSimpleBlock(simpleBlock);
            boundaryFace->SetBoundaryType(bctype);
            boundaryFace->SetBoundaryName(bcname);
            boundaryFace->SetRegion(nst, ned);
            boundaryFace->SetDirection(bcdirection);

            if (IsInterface(bctype))
            {
                int *t_sd = structBC->GetTargetStart();
                int *t_ed = structBC->GetTargetEnd();

                for (int i = 0; i < 3; i ++)
                {
                    nst[i] = t_sd[i];
                    ned[i] = t_ed[i];
                }

                CheckMultigrid(nst, ned, numberOfUnitCells);

                int nbt = structBC->GetTargetRegionBlock();
                int targetbcdirection = structBC->GetFaceDirectionOfTargetBlock();

                boundaryFace->SetNext(new Pre_Patch_Struct());
                Pre_Patch_Struct *simplePatch = boundaryFace->GetNext();

                simplePatch->SetSimpleBlock(referencedBlockContainer[nbt]);

                simplePatch->SetRegion(nst, ned);
                simplePatch->SetTargetBlockLabel(simpleBlock->GetZoneIndex());
                simplePatch->SetDirection(targetbcdirection);
            }
        }

        simpleBlock->NormalizeBoundaryFace(geometricDimension);
    }
}

void Pre_StrGridPartition::CheckMultigrid(vector< int > nst, vector< int > ned, int numberOfUnitCells)
{
    int resi = (nst[0] - ned[0]) % numberOfUnitCells;
    int resj = (nst[1] - ned[1]) % numberOfUnitCells;
    int resk = (nst[2] - ned[2]) % numberOfUnitCells;

    if (abs(resi) + abs(resj) + abs(resk) != 0)
    {
        TK_Exit::PrintDebugInfoExit("Multigrid condition is not satisfied !");
    }

    return;
}

int GetIntegerDigitWidth(int integerData)
{
    return 1 + static_cast< int >(floor(log10(static_cast< RDouble > (max(1, integerData)))));
}

}

