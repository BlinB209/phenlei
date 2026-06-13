#include "Pre_Plot3DToPHengLEI_Struct.h"
#include "TK_Exit.h"
#include "Glb_Dimension.h"
#include "Pre_GridTranslator.h"
#include "TK_Parse.h"
#include "TK_Log.h"

namespace PHSPACE
{
LIB_EXPORT Pre_Plot3DToPHengLEI_Struct::Pre_Plot3DToPHengLEI_Struct(const string &gridFileName) : 
        Pre_GridConversion(gridFileName)
{

}

LIB_EXPORT Pre_Plot3DToPHengLEI_Struct::~Pre_Plot3DToPHengLEI_Struct()
{

}

void Pre_Plot3DToPHengLEI_Struct::ReadGrid()
{
    int fileFormat = GlobalDataBase::GetIntParaFromDB("fileformat");
    //! fileformat : Grid file format.
    //!              1 --- ASCII
    //!              0 --- BINARY
    fileFormat = CheckPlot3DFileIfASCII();
    GlobalDataBase::UpdateData("fileformat", &fileFormat, PHINT, 1);

    if (fileFormat == 1)
    {
        int isPlot3D = 1;
        isPlot3D = CheckGridgenCoorASCIIFileType();
        if (isPlot3D)
        {
            ReadGridgenCoorASCIIPlot3D();
        }
        else
        {
            ReadGridgenCoorASCIIGridgen();
        }
    }
    else if (fileFormat == 0)
    {
        ReadGridgenCoor();
    }
    ReadGridgenBC();
}

int Pre_Plot3DToPHengLEI_Struct::CheckPlot3DFileIfASCII()
{
    fstream file;
    file.open(gridFileName.c_str(), ios_base::in|ios_base::binary);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(gridFileName);
    }

    CharVecSizeType sizeOfIntType = sizeof(int);

    int numberOfBlocks = 0;
    file.read(reinterpret_cast <char *> (&numberOfBlocks), sizeOfIntType);

    file.close();
    file.clear();

    if (numberOfBlocks <= static_cast <int> (1e8))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

int Pre_Plot3DToPHengLEI_Struct::CheckGridgenCoorASCIIFileType()
{
    fstream file;
    file.open(gridFileName.c_str(), ios_base::in);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(gridFileName);
    }

    int numberOfBlocks = 0;
    file >> numberOfBlocks;

    int isPlot3D = 1;
    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        int ni, nj, nk;
        file >> ni;
        file >> nj;
        file >> nk;

        if (ni < 0 || nj < 0 || nk < 0)
        {
            isPlot3D = 0;
            break;
        }
    }

    file.close();
    file.clear();

    return isPlot3D;
}

void Pre_Plot3DToPHengLEI_Struct::ReadGridgenCoorASCIIPlot3D()
{
    fstream file;
    file.open(gridFileName.c_str(), ios_base::in);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(gridFileName);
    }
    file >> nBlocks;

    grids = new Grid *[nBlocks];
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        int ni, nj, nk;
        file >> ni;
        file >> nj;
        file >> nk;

        GridID *index = new GridID(iZone);
        StructGrid *grid = StructGridCast(CreateGridGeneral(STRUCTGRID, index, 0, GetDim()));
        grids[iZone] = grid;
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

        int wordWidth = 8;

        cout << setiosflags(ios::right);
        cout << " iZone = " << setw(wordWidth) << iZone;
        cout << " numberOfBlocks = " << setw(wordWidth) << nBlocks;
        cout << " ni, nj, nk = ";
        cout << setw(wordWidth) << ni;
        cout << setw(wordWidth) << nj;
        cout << setw(wordWidth) << nk;
        cout << "\n";
    }

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        StructGrid *grid = StructGridCast(grids[iZone]);
        grid->ReadXYZ(file);
    }
    file.close();
    file.clear();
}

void Pre_Plot3DToPHengLEI_Struct::ReadGridgenCoorASCIIGridgen()
{
    fstream file;
    file.open(gridFileName.c_str(), ios_base::in);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(gridFileName);
    }
    file >> nBlocks;

    grids = new Grid *[nBlocks];
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        int ni, nj, nk;
        file >> ni;
        file >> nj;
        file >> nk;

        GridID *index = new GridID(iZone);
        StructGrid *grid = StructGridCast(CreateGridGeneral(STRUCTGRID, index, 0, GetDim()));
        grids[iZone] = grid;
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

        int wordWidth = 8;

        cout << setiosflags(ios::right);
        cout << " iZone = " << setw(wordWidth) << iZone;
        cout << " numberOfBlocks = " << setw(wordWidth) << nBlocks;
        cout << " ni, nj, nk = ";
        cout << setw(wordWidth) << ni;
        cout << setw(wordWidth) << nj;
        cout << setw(wordWidth) << nk;
        cout << "\n";

        grid->ReadXYZ(file);
    }

    file.close();
    file.clear();
}

void Pre_Plot3DToPHengLEI_Struct::ReadGridgenCoor()
{
    fstream file;
    file.open(gridFileName.c_str(), ios_base::in|ios_base::binary);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(gridFileName);
    }

    CharVecSizeType sizeOfIntType = sizeof(int);

    file.read(reinterpret_cast <char *> (&nBlocks), sizeOfIntType);

    grids = new Grid *[nBlocks];
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        int ni, nj, nk;
        file.read(reinterpret_cast <char *> (&ni), sizeOfIntType);
        file.read(reinterpret_cast <char *> (&nj), sizeOfIntType);
        file.read(reinterpret_cast <char *> (&nk), sizeOfIntType);

        GridID *index = new GridID(iZone);
        StructGrid *grid = StructGridCast(CreateGridGeneral(STRUCTGRID, index, 0, GetDim()));
        grids[iZone] = grid;
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

        int wordWidth = 8;

        cout << setiosflags(ios::right);
        cout << " iZone = " << setw(wordWidth) << iZone;
        cout << " numberOfBlocks = " << setw(wordWidth) << nBlocks;
        cout << " ni, nj, nk = ";
        cout << setw(wordWidth) << ni;
        cout << setw(wordWidth) << nj;
        cout << setw(wordWidth) << nk;
        cout << "\n";
    }

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        StructGrid *grid = StructGridCast(grids[iZone]);
        grid->ReadXYZ(file);
    }
    file.close();
    file.clear();
}

void Pre_Plot3DToPHengLEI_Struct::ReadGridgenBC()
{
    string boundaryFile = GlobalDataBase::GetStrParaFromDB("bnd_file");
    fstream file;
    file.open(boundaryFile.c_str(), ios_base::in);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(boundaryFile);
    }

    string line, word;
    string separator = " =\t\r\n#$,;";

    ReadNewLine(file, line);
    line = FindNextWord(line, word, separator);
    int flowSolverID;
    from_string <int> (flowSolverID, word, std::dec);
    cout << "flow_solver_id = " << flowSolverID << "\n";

    ReadNewLine(file, line);
    line = FindNextWord(line, word, separator);
    int numberOfBlocks;
    from_string <int> (numberOfBlocks, word, std::dec);
    cout << "number_of_blocks = " << numberOfBlocks << "\n";

    if (numberOfBlocks != nBlocks)
    {
        TK_Exit::ExceptionExit("Block in the boundary file do not match that in the grid file ! \n");
    }

    int nDim = GetDim();

    bool blockProcessFlag = false;

    using namespace PHMPI;

    CreatetZoneProcessorID_INP(numberOfBlocks);
    int *blockProcessorINP = GetZoneProcessorID_INP();

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

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        StructGrid *grid = StructGridCast(grids[iZone]);

        int ni, nj, nk;
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

        cout << "ni,nj,nk = " << ni << " " << nj << " " << nk << "\n";
        int ndif = ABS(ni - grid->GetNI()) + ABS(nj - grid->GetNJ()) + ABS(nk - grid->GetNK());
        if (ndif != 0)
        {
            cout << "Dimensions of zone " << iZone + 1 << " in the boundary file do not match those in the grid file ! \n";
        }

        line = FindNextWord(line, word, separator);

        if (word == "proc")
        {
            blockProcessFlag = true;
            line = FindNextWord(line, word, separator);
            int proc = 0;
            from_string <int> (proc, word, std::dec);
            blockProcessorINP[iZone] = proc;
        }

        string blockname;
        ReadNewLine(file, line);
        line = FindNextWord(line, blockname, separator);
        int numberOfBCRegion;
        ReadNewLine(file, line);
        line = FindNextWord(line, word, separator);
        from_string <int> (numberOfBCRegion, word, std::dec);

        grid->CreateCompositeBCRegion(numberOfBCRegion);
        StructBCSet *structBCSet = grid->GetStructBCSet();
        for (int iBCRegion = 0; iBCRegion < numberOfBCRegion; ++ iBCRegion)
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

            line = FindNextWord(line, word, separator);
            from_string <int> (BCType, word, std::dec);

            StructBC *BCRegions = new StructBC(iZone, iBCRegion);
            structBCSet->SetBCRegion(iBCRegion, BCRegions);
            BCRegions->SetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);
            BCRegions->SetBCType(BCType);

            string bcName;
            line = FindNextWord(line, bcName, separator);
            if (bcName == "")
            {
                if (BCType == -1)
                {
                    bcName = "Connection";
                }
                else
                {
                    bcName = boundaryName[BCType];
                    if ((BCType >7) && (BCType <51))
                    {
                        BCType = 2;
                        BCRegions->SetBCType(BCType);
                    }
                }
                BCRegions->SetBCName(bcName);
            }
            else
            {
                BCRegions->SetBCName(bcName);
            }

            if (IsInterface(BCType))
            {
                int nbt, originalBlockIndex;
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
                from_string <int> (originalBlockIndex, word, std::dec);
                nbt = originalBlockIndex - 1;
                BCRegions->SetTargetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);
                BCRegions->SetTargetRegionBlock(nbt);
                if (nDim == THREE_D)
                {
                    BCRegions->ComputeRelativeParameters();
                }
            }
        }
    }

    if (!blockProcessFlag)
    {
        delete [] blockProcessorINP;    blockProcessorINP = nullptr;
        SetZoneProcessID_INP(0);
    }

    file.close();
    file.clear();
}

void Pre_Plot3DToPHengLEI_Struct::Conversion()
{
    CheckMeshMultigrid();

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        Grid *grid = grids[iZone];
        grid->SetLnkInfo();
        grid->ProcessBCInfo();
    }
}

void Pre_Plot3DToPHengLEI_Struct::CheckMeshMultigrid()
{
    int km, NN;
    int km_grid = 1024;

    for (int iZone = 0; iZone < nBlocks; iZone ++)
    {
        StructGrid *grid = StructGridCast(grids[iZone]);

        StructBCSet *structBCSet = grid->GetStructBCSet();

        int nBCRegion = structBCSet->GetnBCRegion();
        for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion ++)
        {
            km = 1;
            NN = 2;

            StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);
            int *s_st = structBC->GetStartPoint();
            int *s_ed = structBC->GetEndPoint();

            while ((((abs(s_st[0])-1)%NN)==0) && (((abs(s_ed[0])-1)%NN)==0) && (((abs(s_st[1])-1)%NN)==0) && (((abs(s_ed[1])-1)%NN)==0) && (((abs(s_st[2])-1)%NN)==0) && (((abs(s_ed[2])-1)%NN)==0))
            {
                NN = NN * 2;
                km = km + 1;
            }
            km_grid = min(km_grid, km);
        }
    }

    PrintToWindow("The most valid multi-grid level is: ", km_grid, "\n");
}

}