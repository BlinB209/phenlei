#include "Glb_Dimension.h"
#include "GridPartition.h"
#include "Pre_UnstructGridPartition.h"
#include "TK_Parse.h"
#include "GridgenIn.h"
#include "Pre_StructGridPartition.h"
#include "GridType.h"
#include "TK_Log.h"
#include "IO_FileName.h"

using namespace std;

#pragma warning (disable:913)
namespace PHSPACE
{
void PartitionGrid()
{
    int axisup = 1;
    GlobalDataBase::UpdateData("axisup", &axisup, PHINT, 1);

    int numberOfGridFile = GlobalDataBase::GetIntParaFromDB("numberOfGridFile");
    if (numberOfGridFile > 1)
    {
        MultiGridPartition();
        return;
    }

    int partitionedGridType = GlobalDataBase::GetIntParaFromDB("pgridtype");
    if (partitionedGridType == UNSTRUCTGRID)
    {
        PartitionUnstructGrid();
    }
    else if (partitionedGridType == STRUCTGRID)
    {
        PartitionStructGrid();
    }
}

void MultiGridPartition()
{
    PHVectorString1D gridFileNameList;
    PHVectorString1D outFileNameList;
    PHVectorInt1D gridTypeList;
    PHVectorInt1D maxProcList;

    GetPartitionPara(gridFileNameList, outFileNameList, gridTypeList, maxProcList);
    int axisup = 1;
    GlobalDataBase::UpdateData("axisup", &axisup, PHINT, 1);

    int numberOfGridFile = GlobalDataBase::GetIntParaFromDB("numberOfGridFile");
    for (int iGridFile = 0; iGridFile < numberOfGridFile; ++ iGridFile)
    {
        int pgridtype = gridTypeList[iGridFile];
        int maxproc = maxProcList[iGridFile];
        string original_grid_file  = gridFileNameList[iGridFile];
        string partition_grid_file = outFileNameList[iGridFile];

        GlobalDataBase::UpdateData("pgridtype", &pgridtype, PHINT, 1);
        GlobalDataBase::UpdateData("maxproc", &maxproc, PHINT, 1);
        GlobalDataBase::UpdateData("original_grid_file", &original_grid_file, PHSTRING, 1);
        GlobalDataBase::UpdateData("partition_grid_file", &partition_grid_file, PHSTRING, 1);

        if (pgridtype == UNSTRUCTGRID)
        {
            PartitionUnstructGrid();
        }
        else if (pgridtype == STRUCTGRID)
        {
            PartitionStructGrid();
        }

        DeleteGGrids();
        PHMPI::SetNumberofLocalZones(0);
    }
}

void GetPartitionPara(PHVectorString1D &gridFileNameList, PHVectorString1D &outFileNameList, PHVectorInt1D &gridTypeList, PHVectorInt1D &maxProcList)
{
    gridFileNameList.push_back(PHSPACE::GlobalDataBase::GetStrParaFromDB("original_grid_file"));
    gridTypeList.push_back(PHSPACE::GlobalDataBase::GetIntParaFromDB("pgridtype"));
    maxProcList.push_back(PHSPACE::GlobalDataBase::GetIntParaFromDB("maxproc"));
    outFileNameList.push_back(AddSymbolToFileName(gridFileNameList[0], "__", maxProcList[0]));

    int numberOfGridFile = GlobalDataBase::GetIntParaFromDB("numberOfGridFile");
    for (int iGridFile = 1; iGridFile < numberOfGridFile; ++ iGridFile)
    {
        ostringstream ossFileName;
        ossFileName << "original_grid_file" << iGridFile;
        gridFileNameList.push_back(PHSPACE::GlobalDataBase::GetStrParaFromDB(ossFileName.str()));

        ostringstream ossGridType;
        ossGridType << "pgridtype" << iGridFile;
        gridTypeList.push_back(PHSPACE::GlobalDataBase::GetIntParaFromDB(ossGridType.str()));

        ostringstream ossMaxProc;
        ossMaxProc << "maxproc" << iGridFile;
        maxProcList.push_back(PHSPACE::GlobalDataBase::GetIntParaFromDB(ossMaxProc.str()));

        outFileNameList.push_back(AddSymbolToFileName(gridFileNameList[iGridFile], "__", maxProcList[iGridFile]));
    }
}

void PartitionStructGrid()
{
    int maximumProcess = PHSPACE::GlobalDataBase::GetIntParaFromDB("maxproc");

    if (maximumProcess < 2)
    {
        return;
    }

    Pre_StrGridPartition blockSplittingManager(maximumProcess);

    blockSplittingManager.Run();
}

void PartitionUnstructGrid()
{
    string originalGridFile = GlobalDataBase::GetStrParaFromDB("original_grid_file");
    GlobalDataBase::UpdateData("gridfile", &originalGridFile, PHSTRING, 1);

    string partitionGridFile = GlobalDataBase::GetStrParaFromDB("partition_grid_file");

    int nProcessor = GlobalDataBase::GetIntParaFromDB("maxproc");
    int partMethod = GlobalDataBase::GetIntParaFromDB("npartmethod");

    Grid **parititionGrids = 0;

    Region *region = new Region();

    int numberOfProcessor = PHMPI::GetNumberOfProcessor();
    if (numberOfProcessor == 1)
    {
        //! Serial partition.
        //! Step1: Read the original grid.
        region->ReadGrid();

        int nZones = PHMPI::GetNumberofGlobalZones();
        PrintToWindow("Number of original && final partitions: ", nZones, nProcessor, "\n");
        if (nZones > 1)
        {
            TK_Exit::ExceptionExit("Error: could not partition Multi Zones serially, using parallel please!");
        }

        //! Step2: partition grid serially.
        Grid *grid = GetGrid(0, 0);

        RDouble nCellsPerProc = 100000.0;
        if (GlobalDataBase::IsExist("nCellsPerProc", PHDOUBLE, 1))
        {
            GlobalDataBase::GetData("nCellsPerProc", &nCellsPerProc, PHDOUBLE, 1);
            int nTotalCell = grid->GetNTotalCell();
            RDouble nprocessor = static_cast<double>(nTotalCell) / nCellsPerProc;
            nProcessor = static_cast<int>(round(nprocessor));
            if (nProcessor < 2)
            {
                PrintToWindow("Number of original && final partitions due to nCellsPerProc: ", nZones, nProcessor, "\n");
                return;
            }
            partitionGridFile = originalGridFile;
            partitionGridFile = AddSymbolToFileName(partitionGridFile, "__", nProcessor);
            PrintToWindow("Number of original && final partitions due to nCellsPerProc: ", nZones, nProcessor, "\n");
        }

        Pre_UnstructGridPartition unstructGridParition(UnstructGridCast(grid), nProcessor, partMethod, originalGridFile);

        unstructGridParition.PartitionGridSerial(partitionGridFile);

        unstructGridParition.DumpPartitionGrids(partitionGridFile);

        parititionGrids = unstructGridParition.GetPartitionGrids();
    }
    else
    {
        int parallelPartitionMethod = METIS_METHOD;
        GlobalDataBase::GetData("parallelPartitionMethod", &parallelPartitionMethod, PHINT, 1);

        //! Parallel partition.
        //! Step1: Read the original grid.
        region->ReadGrid();

        int nZones = PHMPI::GetNumberofGlobalZones();
        PrintToWindow("Number of original && final partitions: ", nZones, nProcessor, "\n");
        if (nZones != numberOfProcessor)
        {
            TK_Exit::ExceptionExit("Error: nZones must be equal to nProcessor for parallel partition!");
        }

        region->SetupMultiBlockRelationship();

        //! Find out the grid in this processor.
        //! Warning: only ONE grid in each processor in case of parallel partition.
        Grid *gridOfThisProcessor = 0;
        int myID = PHMPI::GetCurrentProcessorID();
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int procID = PHMPI::GetZoneProcessorID(iZone);

            if (myID == procID)
            {
                gridOfThisProcessor = GetGrid(iZone, 0);
                break;
            }
        }

        //! Step2: partition grid parallel.
        Pre_UnstructGridPartition parallelUnstructGridParition(UnstructGridCast(gridOfThisProcessor), nProcessor, partMethod, originalGridFile, parallelPartitionMethod);

        parallelUnstructGridParition.PartitionGridParallel(partitionGridFile);

        parallelUnstructGridParition.DumpPartitionGrids(partitionGridFile);

        parititionGrids = parallelUnstructGridParition.GetPartitionGrids();
    }

    //! Free the partition grids.
    for (int iProc = 0; iProc < nProcessor; ++ iProc)
    {
        delete parititionGrids[iProc];
    }
    delete [] parititionGrids;

    delete region;
}

void RefineStructGrid()
{
    string coordinateFile = GlobalDataBase::GetStrParaFromDB("coor_file");
    string layoutFile = GlobalDataBase::GetStrParaFromDB("layout_file");

    PHSPACE::RefineStructGrid(coordinateFile, layoutFile);
}

void ReadGridgenBC(const string &BCFile, DataContainer *cData)
{
    fstream file;
    file.open(BCFile.c_str(),ios_base::in);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(BCFile);
    }

    ReadGridgenBC(file, cData);

    file.close();
    file.clear();
}

void ReadGridgenBC(fstream &file, DataContainer *cDdata)
{
    string line, word;
    string separator = " =\t\r\n#$,;";

    CharVecSizeType sizeOfIntType = sizeof(int);

    ReadNewLine(file, line);
    line = FindNextWord(line, word, separator);
    int flowSolverId;
    from_string <int> (flowSolverId, word, std::dec);
    cDdata->Write(&flowSolverId, sizeOfIntType);
    cout << "flow_solver_id = " << flowSolverId << "\n";

    ReadNewLine(file, line);
    line = FindNextWord(line, word, separator);
    int numberOfBlocks;
    from_string <int> (numberOfBlocks, word, std::dec);

    cDdata->Write(&numberOfBlocks, sizeOfIntType);

    cout << "number_of_blocks = " << numberOfBlocks << "\n";

    int nDim = GetDim();

    for (int i = 0; i < numberOfBlocks; ++ i)
    {
        int ni,nj,nk;

        ReadNewLine(file, line);
        line = FindNextWord(line, word, separator);
        from_string <int> (ni, word, std::dec);
        line = FindNextWord(line, word, separator);
        from_string <int> (nj, word, std::dec);
        nk = 1;
        if (nDim == THREE_D)
        {
            line = FindNextWord(line, word, separator);
            from_string <int> (nk, word, std::dec);
        }

        cDdata->Write(&ni, sizeOfIntType);
        cDdata->Write(&nj, sizeOfIntType);
        cDdata->Write(&nk, sizeOfIntType);

        int wordwidth = 8;

        cout << setiosflags(ios::right);
        cout << " i = " << setw(wordwidth) << i;
        cout << " nBlocks = " << setw(wordwidth) << numberOfBlocks;
        cout << " ni, nj, nk = ";
        cout << setw(wordwidth) << ni;
        cout << setw(wordwidth) << nj;
        cout << setw(wordwidth) << nk;
        cout << "\n";

        string blockName;
        ReadNewLine(file, line);
        line = FindNextWord(line, blockName, separator);
        cDdata->WriteString(blockName);

        int nBCRegion;
        ReadNewLine(file, line);
        line = FindNextWord(line, word, separator);
        from_string <int> (nBCRegion, word, std::dec);
        cDdata->Write(&nBCRegion, sizeOfIntType);

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            int imin, imax, jmin, jmax, kmin, kmax, BCType;
            ReadNewLine(file, line);
            line = FindNextWord(line, word, separator);
            from_string <int> (imin, word, std::dec);
            line = FindNextWord(line, word, separator);
            from_string <int> (imax, word, std::dec);

            line = FindNextWord(line, word, separator);
            from_string <int> (jmin, word, std::dec);
            line = FindNextWord(line, word, separator);
            from_string <int> (jmax, word, std::dec);

            if (nDim == THREE_D)
            {
                line = FindNextWord(line, word, separator);
                from_string <int> (kmin, word, std::dec);
                line = FindNextWord(line, word, separator);
                from_string <int> (kmax, word, std::dec);
            }
            else
            {
                kmin = 1;
                kmax = 1;
            }

            cDdata->Write(&imin, sizeOfIntType);
            cDdata->Write(&imax, sizeOfIntType);
            cDdata->Write(&jmin, sizeOfIntType);
            cDdata->Write(&jmax, sizeOfIntType);
            cDdata->Write(&kmin, sizeOfIntType);
            cDdata->Write(&kmax, sizeOfIntType);

            line = FindNextWord(line, word, separator);
            from_string <int> (BCType, word, std::dec);

            cDdata->Write(&BCType, sizeOfIntType);

            if (IsInterface(BCType))
            {
                int nbt;
                ReadNewLine(file, line);
                line = FindNextWord(line, word, separator);
                from_string <int> (imin, word, std::dec);
                line = FindNextWord(line, word, separator);
                from_string <int> (imax, word, std::dec);

                line = FindNextWord(line, word, separator);
                from_string <int> (jmin, word, std::dec);
                line = FindNextWord(line, word, separator);
                from_string <int> (jmax, word, std::dec);

                if (nDim == THREE_D)
                {
                    line = FindNextWord(line, word, separator);
                    from_string <int> (kmin, word, std::dec);
                    line = FindNextWord(line, word, separator);
                    from_string <int> (kmax, word, std::dec);
                }
                else
                {
                    kmin = 1;
                    kmax = 1;
                }

                line = FindNextWord(line, word, separator);
                from_string <int> (nbt, word, std::dec);

                cDdata->Write(&imin, sizeOfIntType);
                cDdata->Write(&imax, sizeOfIntType);
                cDdata->Write(&jmin, sizeOfIntType);
                cDdata->Write(&jmax, sizeOfIntType);
                cDdata->Write(&kmin, sizeOfIntType);
                cDdata->Write(&kmax, sizeOfIntType);

                cDdata->Write(&nbt, sizeOfIntType);
            }
        }
    }
}

void BlockRefine::SetNewDimension()
{
    int ni = dim[0];
    int nj = dim[1];
    int nk = dim[2];

    oldDim[0] = dim[0];
    oldDim[1] = dim[1];
    oldDim[2] = dim[2];

    dim[0] = 2 * (ni - 1) + 1;
    dim[1] = 2 * (nj - 1) + 1;
    dim[2] = 2 * (nk - 1) + 1;
}

void BlockRefine::Refine()
{
    uint_t nBCRegion = BCList.size();

    for (uint_t iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        BCList[iBCRegion]->Refine();
    }
}

void BlockRefine::RefineCoordinate(RDouble3D &x, RDouble3D &xold, const int &ni, const int &nj, const int &nk)
{
    int istp = 2;
    int jstp = 2;
    int kstp = 2;

    int i, j, k, ii, jj, kk;
    for (k = 1; k <= nk; k += kstp)
    {
        for (j = 1; j <= nj; j += jstp)
        {
            for (i = 1; i <= ni; i += istp)
            {
                ii = (i + 1) / 2;
                jj = (j + 1) / 2;
                kk = (k + 1) / 2;
                x(i, j, k) = xold(ii, jj, kk);
            }
        }
    }

    for (k = 1; k <= nk; k += kstp)
    {
        for (j = 1; j <= nj; j += jstp)
        {
            for (i = 2; i <= ni; i += istp)
            {
                x(i, j, k) = half * (x(i-1, j, k) + x(i+1, j, k));
            }
        }
    }

    for (k = 1; k <= nk; k += kstp)
    {
        for (j = 2; j <= nj; j += jstp)
        {
            for (i = 1; i <= ni; i += istp)
            {
                x(i, j, k) = half * (x(i, j-1, k) + x(i, j+1, k));
            }
        }
    }

    for (k = 2; k <= nk; k += kstp)
    {
        for (j = 1; j <= nj; j += jstp)
        {
            for (i = 1; i <= ni; i += istp)
            {
                x(i, j, k) = half * (x(i, j, k-1) + x(i, j, k+1));
            }
        }
    }

    for (k = 1; k <= nk; k += kstp)
    {
        for (j = 2; j <= nj; j += jstp)
        {
            for (i = 2; i <= ni; i += istp)
            {
                x(i, j, k) = fourth * (x(i-1, j-1, k) + x(i+1, j+1, k) 
                                     + x(i+1, j-1, k) + x(i-1, j+1, k));
            }
        }
    }

    for (k = 2; k <= nk; k += kstp)
    {
        for (j = 1; j <= nj; j += jstp)
        {
            for (i = 2; i <= ni; i += istp)
            {
                x(i, j, k) = fourth * (x(i-1, j, k-1) + x(i+1, j, k+1) 
                                     + x(i+1, j, k-1) + x(i-1, j, k+1));
            }
        }
    }

    for (k = 2; k <= nk; k += kstp)
    {
        for (j = 2; j <= nj; j += jstp)
        {
            for (i = 1; i <= ni; i += istp)
            {
                x(i, j, k) = fourth * (x(i, j-1 ,k-1) + x(i, j+1, k+1) 
                                     + x(i, j+1, k-1) + x(i, j-1, k+1));
            }
        }
    }

    for (k = 2; k <= nk; k += kstp)
    {
        for (j = 2; j <= nj; j += jstp)
        {
            for (i = 2; i <= ni; i += istp)
            {
                x(i, j, k) = eighth * (x(i-1, j-1, k-1) + x(i+1, j+1, k-1) 
                                     + x(i-1, j+1, k-1) + x(i+1, j-1, k-1)
                                     + x(i-1, j-1, k+1) + x(i+1, j+1, k+1) 
                                     + x(i-1, j+1, k+1) + x(i+1, j-1, k+1));
            }
        }
    }
}

StructGrid *BlockRefine::CreateStructGrid(const int &ni, const int &nj, const int &nk)
{
    GridID *index = new GridID(0);
    StructGrid *grid = StructGridCast(CreateGridGeneral(STRUCTGRID, index, 0, GetDim()));

    grid->SetNI(ni);
    grid->SetNJ(nj);
    grid->SetNK(nk);

    grid->SetBasicDimension();

    int nTotalNode = grid->GetNTotalNode();
    RDouble *x = new RDouble[nTotalNode];
    RDouble *y = new RDouble[nTotalNode];
    RDouble *z = new RDouble[nTotalNode];
    grid->SetX(x);
    grid->SetY(y);
    grid->SetZ(z);
    grid->SetArrayLayout();

    return grid;
}

void BlockRefine::DumpGrid(fstream &file, Grid *gridIn)
{
    StructGrid *gridOld = StructGridCast(gridIn);

    //! ni,nj,nk are fine grid dimension
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    StructGrid *grid = CreateStructGrid(ni,nj,nk);

    RDouble3D &xx = *(grid->GetStructX());
    RDouble3D &yy = *(grid->GetStructY());
    RDouble3D &zz = *(grid->GetStructZ());

    RDouble3D &xxOld = *(gridOld->GetStructX());
    RDouble3D &yyOld = *(gridOld->GetStructY());
    RDouble3D &zzOld = *(gridOld->GetStructZ());

    RefineCoordinate(xx, xxOld, ni, nj, nk);
    RefineCoordinate(yy, yyOld, ni, nj, nk);
    RefineCoordinate(zz, zzOld, ni, nj, nk);

    grid->WriteXYZ(file, 1, ni, 1, nj, 1, nk);

    delete grid;
}

void BlockBCRefine::RefineBC()
{
    for (int m = 0; m < 3; ++ m)
    {
        st[m] = SIGN(1, st[m]) * (2 * ABS(st[m]) - 1);
        ed[m] = SIGN(1, ed[m]) * (2 * ABS(ed[m]) - 1);
    }
}

void BlockBCRefine::Refine()
{
    RefineBC();
    if (next)
    {
        next->RefineBC();
    }
}

void BlockBCRefine::SetContent(const int &BCType, const int &imin, const int &imax, 
                                 const int &jmin, const int &jmax, const int &kmin, const int &kmax)
{
    this->BCType = BCType;

    this->st[0]  = imin;
    this->ed[0]  = imax;

    this->st[1]  = jmin;
    this->ed[1]  = jmax;

    this->st[2]  = kmin;
    this->ed[2]  = kmax;
}

void BlockBCRefine::SetPatchContent(const int &nb, const int &imin, const int &imax, 
                                      const int &jmin, const int &jmax, const int &kmin, const int &kmax)
{
    if (!this->next)
    {
        this->next = new BlockBCPatchRefine();
    }

    this->next->SetContent(nb, imin, imax, jmin, jmax, kmin, kmax);
}

void BlockBCPatchRefine::RefineBC()
{
    for (int m = 0; m < 3; ++ m)
    {
        st[m] = SIGN(1, st[m]) * (2 * ABS(st[m]) - 1);
        ed[m] = SIGN(1, ed[m]) * (2 * ABS(ed[m]) - 1);
    }
}

void BlockBCPatchRefine::SetContent(const int &BCType, const int &imin, const int &imax, 
                                      const int &jmin, const int &jmax, const int &kmin, const int &kmax)
{
    this->nb    = BCType;

    this->st[0] = imin;
    this->ed[0] = imax;

    this->st[1] = jmin;
    this->ed[1] = jmax;

    this->st[2] = kmin;
    this->ed[2] = kmax;
}

void GridRefine::Init(DataContainer *cDdata)
{
    CharVecSizeType sizeOfIntType = sizeof(int);

    cDdata->MoveToBegin();
    int flowSolverId;
    cDdata->Read(&flowSolverId, sizeOfIntType);

    cDdata->Read(&numberOfBlocks, sizeOfIntType);

    cout << "number_of_blocks = " << numberOfBlocks << "\n";

    refblock.resize(numberOfBlocks);

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        refblock[iZone] = new BlockRefine();
        refblock[iZone]->SetID(iZone);
    }

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        int ni, nj, nk;

        cDdata->Read(&ni, sizeOfIntType);
        cDdata->Read(&nj, sizeOfIntType);
        cDdata->Read(&nk, sizeOfIntType);

        string blockname;
        cDdata->ReadString(blockname);

        BlockRefine *block = refblock[iZone];

        block->SetNI(ni);
        block->SetNJ(nj);
        block->SetNK(nk);

        int nBCRegion;
        cDdata->Read(&nBCRegion, sizeOfIntType);

        vector <BlockBCRefine *> &BCList = block->GetBCList();

        BCList.resize(nBCRegion);

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            int imin, imax, jmin, jmax, kmin, kmax, BCType;

            cDdata->Read(&imin, sizeOfIntType);
            cDdata->Read(&imax, sizeOfIntType);
            cDdata->Read(&jmin, sizeOfIntType);
            cDdata->Read(&jmax, sizeOfIntType);
            cDdata->Read(&kmin, sizeOfIntType);
            cDdata->Read(&kmax, sizeOfIntType);
            cDdata->Read(&BCType, sizeOfIntType);

            BlockBCRefine *BCRefined = new BlockBCRefine();
            BCList[iBCRegion] = BCRefined;

            BCRefined->SetContent(BCType, imin, imax, jmin, jmax, kmin, kmax);

            if (IsInterface(BCType))
            {
                int indexOfTargetBoundary;
                cDdata->Read(&imin, sizeOfIntType);
                cDdata->Read(&imax, sizeOfIntType);
                cDdata->Read(&jmin, sizeOfIntType);
                cDdata->Read(&jmax, sizeOfIntType);
                cDdata->Read(&kmin, sizeOfIntType);
                cDdata->Read(&kmax, sizeOfIntType);
                cDdata->Read(&indexOfTargetBoundary , sizeOfIntType);

                BCRefined->SetPatchContent(BCType, imin, imax, jmin, jmax, kmin, kmax);
            }
        }
    }
}

void GridRefine::SetNewDimension()
{
    uint_t numberOfBlocks = refblock.size();
    for (uint_t iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        BlockRefine *blockToRefin = refblock[iZone];
        blockToRefin->SetNewDimension();
    }
}

void GridRefine::Refine()
{
    SetNewDimension();
    uint_t numberOfBlocks = refblock.size();
    for (uint_t iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        BlockRefine *blockToRefin = refblock[iZone];
        blockToRefin->Refine();
    }
}

int GridRefine::MapBC(int i)
{
    return SIGN(1,i) * (2 * ABS(i) - 1);
}

void GridRefine::DumpLayout(DataContainer *cDdata)
{
    fstream file;
    string partLayoutFile = GlobalDataBase::GetStrParaFromDB("part_layout_file");

    file.open(partLayoutFile.c_str(), ios_base::out);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(partLayoutFile);
    }

    CharVecSizeType sizeOfIntType = sizeof(int);

    cDdata->MoveToBegin();

    file << setiosflags(ios::right);
    int wordwidth = 8;

    int flowSolverID;
    cDdata->Read(&flowSolverID, sizeOfIntType);

    file << setw(wordwidth) << flowSolverID << "\n";

    int numberOfBlocks;
    cDdata->Read(&numberOfBlocks, sizeOfIntType);

    file << setw(wordwidth) << numberOfBlocks << "\n";

    cout << "number_of_blocks = " << numberOfBlocks << "\n";

    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        int ni, nj, nk;

        cDdata->Read(&ni, sizeOfIntType);
        cDdata->Read(&nj, sizeOfIntType);
        cDdata->Read(&nk, sizeOfIntType);

        file << setw(wordwidth) << 2 * ni - 1;
        file << setw(wordwidth) << 2 * nj - 1;

        if (GetDim() == THREE_D)
        {
            file << setw(wordwidth) << 2 * nk - 1;
        }

        file << "\n";

        string blockName;
        cDdata->ReadString(blockName);

        file << blockName << "\n";

        int nBCRegion;
        cDdata->Read(&nBCRegion, sizeOfIntType);

        file << setw(wordwidth) << nBCRegion << "\n";

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            int imin, imax, jmin, jmax, kmin, kmax, BCType;

            cDdata->Read(&imin, sizeOfIntType);
            cDdata->Read(&imax, sizeOfIntType);
            cDdata->Read(&jmin, sizeOfIntType);
            cDdata->Read(&jmax, sizeOfIntType);
            cDdata->Read(&kmin, sizeOfIntType);
            cDdata->Read(&kmax, sizeOfIntType);
            cDdata->Read(&BCType, sizeOfIntType);

            file << setw(wordwidth) << MapBC(imin);
            file << setw(wordwidth) << MapBC(imax);
            file << setw(wordwidth) << MapBC(jmin);
            file << setw(wordwidth) << MapBC(jmax);

            if (GetDim() == THREE_D)
            {
                file << setw(wordwidth) << MapBC(kmin);
                file << setw(wordwidth) << MapBC(kmax);
            }

            file << setw(wordwidth) << BCType << "\n";

            if (IsInterface(BCType))
            {
                int nbt;
                cDdata->Read(&imin, sizeOfIntType);
                cDdata->Read(&imax, sizeOfIntType);
                cDdata->Read(&jmin, sizeOfIntType);
                cDdata->Read(&jmax, sizeOfIntType);
                cDdata->Read(&kmin, sizeOfIntType);
                cDdata->Read(&kmax, sizeOfIntType);
                cDdata->Read(&nbt , sizeOfIntType);

                file << setw(wordwidth) << MapBC(imin);
                file << setw(wordwidth) << MapBC(imax);
                file << setw(wordwidth) << MapBC(jmin);
                file << setw(wordwidth) << MapBC(jmax);

                if (GetDim() == THREE_D)
                {
                    file << setw(wordwidth) << MapBC(kmin);
                    file << setw(wordwidth) << MapBC(kmax);
                }
                file << setw(wordwidth) << nbt << "\n";
            }
        }
    }

    file.close();
    file.clear();
}

void GridRefine::DumpGrid()
{
    string coordinateFile            = GlobalDataBase::GetStrParaFromDB("coor_file");
    string partitionedCoordinateFile = GlobalDataBase::GetStrParaFromDB("part_coor_file");

    CharVecSizeType sizeOfIntType = sizeof(int);

    uint_t numberOfBlocks = refblock.size();
    int nBlocks;
    Grid **structGrids;
    ReadPlot3DCoorBinary(coordinateFile, nBlocks, structGrids);

    fstream file;
    file.open(partitionedCoordinateFile.c_str(),ios_base::out|ios_base::binary|ios_base::trunc);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(partitionedCoordinateFile);
    }

    file.write(reinterpret_cast <char *> (&nBlocks), sizeOfIntType);

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        int ni = refblock[iZone]->GetNI();
        int nj = refblock[iZone]->GetNJ();
        int nk = refblock[iZone]->GetNK();
        file.write(reinterpret_cast <char *> (&ni), sizeOfIntType);
        file.write(reinterpret_cast <char *> (&nj), sizeOfIntType);
        file.write(reinterpret_cast <char *> (&nk), sizeOfIntType);
    }

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        BlockRefine *blockRefined = refblock[iZone];
        blockRefined->DumpGrid(file,structGrids[iZone]);
    }

    file.close();

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        delete structGrids[iZone];
    }
    delete [] structGrids;
}

void GridRefine::RefineStructGrid()
{
    string layoutFile = GlobalDataBase::GetStrParaFromDB("layout_file");

    DataContainer *cDdata = new DataContainer();

    ReadGridgenBC(layoutFile, cDdata);

    this->Init(cDdata);
    this->Refine();
    this->DumpLayout(cDdata);
    this->DumpGrid();

    delete cDdata;
}
#pragma warning(disable:4100)
void RefineStructGrid(const string &coordinateFile, const string &layoutFile)
{
    GridRefine *gridToRefine = new GridRefine;

    gridToRefine->RefineStructGrid();

    delete gridToRefine;
}
#pragma warning(default:4100)
int in_set(int *a, int n, int value)
{
    for (int i = 0; i < n; ++ i)
    {
        if (a[i] == value)
        {
            return 1;
        }
    }
    return 0;
}

void GetIJKRegion(int *dim, int &ist, int &ied, int &jst, int &jed, int &kst, int &ked, int m)
{
    ist = 1;
    ied = dim[0];

    jst = 1;
    jed = dim[1];

    kst = 1;
    ked = dim[2];

    if (m == 0)
    {
        ied = 1;
    }
    else if (m == 1)
    {
        ist = dim[0];
    }
    else if (m == 2)
    {
        jed = 1;
    }  
    else if (m == 3)
    {
        jst = dim[1];
    }
    else if (m == 4)
    {
        ked = 1;
    }
    else if (m == 5)
    {
        kst = dim[2];
    }
}

void GetIJKRegion(int *st, int *ed, int &ist, int &ied, int &jst, int &jed, int &kst, int &ked, int m)
{
    ist = st[0];
    ied = ed[0];

    jst = st[1];
    jed = ed[1];

    kst = st[2];
    ked = ed[2];

    if (m == 0)
    {
        ied = st[0];
    }
    else if (m == 1)
    {
        ist = ed[0];
    }
    else if (m == 2)
    {
        jed = st[1];
    }
    else if (m == 3)
    {
        jst = ed[1];
    }
    else if (m == 4)
    {
        ked = st[2];
    }
    else if (m == 5)
    {
        kst = ed[2];
    }
}

void GetIJKRegion(int *box_st, int *box_ed, int *st, int *ed, int m)
{
    GetIJKRegion(box_st, box_ed, st[0], ed[0], st[1], ed[1], st[2], ed[2], m);
}

}
