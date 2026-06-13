#include <string>
#include "GridgenIn.h"
#include "Geo_StructBC.h"
#include "TK_Parse.h"
#include "Glb_Dimension.h"
#include "Math_BasisFunction.h"
#include "TK_Exit.h"
#include "Geo_StructGrid.h"
#include "GridType.h"

#pragma warning (disable:913)
namespace PHSPACE
{
void ReadPlot3dGrid(Grid **&structGrids, int &numberOfBlocks, const string &gridgenFile, const string &BCFile)
{
    int fileFormat = GlobalDataBase::GetIntParaFromDB("fileformat");
    //! Fileformat: Grid file format.
    //!             1 --- ASCII
    //!             0 --- BINARY
    fileFormat = CheckCoorFileIfASCII(gridgenFile);
    GlobalDataBase::UpdateData("fileformat", &fileFormat, PHINT, 1);

    if (fileFormat == 1)
    {
        ReadPlot3DCoorASCII(gridgenFile, numberOfBlocks, structGrids);
    }
    else if (fileFormat == 0)
    {
        ReadPlot3DCoorBinary(gridgenFile, numberOfBlocks, structGrids);
    }
    ReadBoundaryFile(structGrids, numberOfBlocks, BCFile);
}

void ReadGridgenGrid(Grid **&structGrids, int &numberOfBlocks, const string &gridgenFile, const string &BCFile)
{
    int fileFormat = GlobalDataBase::GetIntParaFromDB("fileformat");
    //! Fileformat: Grid file format.
    //!             1 --- ASCII
    //!             0 --- BINARY
    fileFormat = CheckCoorFileIfASCII(gridgenFile);
    GlobalDataBase::UpdateData("fileformat", &fileFormat, PHINT, 1);

    if (fileFormat == 1)
    {
        ReadGridgenCoorASCII(gridgenFile, numberOfBlocks, structGrids);
    }
    else
    {
        ostringstream oss;
        oss << "  Error: this situation has not been considered, for file type = BINARY, when from_gtype = Gridgen." << endl;
        TK_Exit::ExceptionExit(oss.str());
    }

    ReadBoundaryFile(structGrids, numberOfBlocks, BCFile);
}

void ReadPlot3DCoorBinary(const string &coordinateFile, int &numberOfBlocks, Grid **&structGrids)
{
    fstream file;
    file.open(coordinateFile.c_str(), ios_base::in|ios_base::binary);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(coordinateFile);
    }

    CharVecSizeType sizeOfIntType = sizeof(int);

    file.read(reinterpret_cast <char *> (&numberOfBlocks), sizeOfIntType);

    structGrids = new Grid *[numberOfBlocks];
    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        int ni, nj, nk;
        file.read(reinterpret_cast <char *> (&ni), sizeOfIntType);
        file.read(reinterpret_cast <char *> (&nj), sizeOfIntType);
        file.read(reinterpret_cast <char *> (&nk), sizeOfIntType);

        GridID *index = new GridID(iZone);
        StructGrid *grid = StructGridCast(CreateGridGeneral(STRUCTGRID, index, 0, GetDim()));
        structGrids[iZone] = grid;
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
        cout << " numberOfBlocks = " << setw(wordWidth) << numberOfBlocks;
        cout << " ni, nj, nk = ";
        cout << setw(wordWidth) << ni;
        cout << setw(wordWidth) << nj;
        cout << setw(wordWidth) << nk;
        cout << "\n";
        //cin.get();

        //delete index;
    }

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        StructGrid *grid = StructGridCast(structGrids[iZone]);
        grid->ReadXYZ(file);
    }
    file.close();
    file.clear();
}

void ReadGridgenCoorASCII(const string &coordinateFile, int &numberOfBlocks, Grid **&structGrids)
{
    fstream file;
    file.open(coordinateFile.c_str(),ios_base::in);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(coordinateFile);
    }
    file >> numberOfBlocks;

    structGrids = new Grid *[numberOfBlocks];
    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        int ni, nj, nk;
        file >> ni;
        file >> nj;
        file >> nk;

        GridID *index = new GridID(iZone);
        StructGrid *grid = StructGridCast(CreateGridGeneral(STRUCTGRID, index, 0, GetDim()));
        structGrids[iZone] = grid;
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
        cout << " numberOfBlocks = " << setw(wordWidth) << numberOfBlocks;
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

int CheckCoorFileIfASCII(const string &coordinateFile)
{
    fstream file;
    file.open(coordinateFile.c_str(), ios_base::in|ios_base::binary);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(coordinateFile);
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

void ReadPlot3DCoorASCII(const string &coordinateFile, int &numberOfBlocks, Grid **&structGrids)
{
    fstream file;
    file.open(coordinateFile.c_str(), ios_base::in);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(coordinateFile);
    }
    file >> numberOfBlocks;

    structGrids = new Grid *[numberOfBlocks];
    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        int ni, nj, nk;
        file >> ni;
        file >> nj;
        file >> nk;

        GridID *index = new GridID(iZone);
        StructGrid *grid = StructGridCast(CreateGridGeneral(STRUCTGRID, index, 0, GetDim()));
        structGrids[iZone] = grid;
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
        cout << " numberOfBlocks = " << setw(wordWidth) << numberOfBlocks;
        cout << " ni, nj, nk = ";
        cout << setw(wordWidth) << ni;
        cout << setw(wordWidth) << nj;
        cout << setw(wordWidth) << nk;
        cout << "\n";
    }

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        StructGrid *grid = StructGridCast(structGrids[iZone]);
        grid->ReadXYZ(file);
    }
    file.close();
    file.clear();
}

void ReadBoundaryFile(Grid **&structGrids, int &inputNumberOfBlocks, const string &BCFile)
{
    fstream file;
    file.open(BCFile.c_str(), ios_base::in);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(BCFile);
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

    if (numberOfBlocks != inputNumberOfBlocks)
    {
        TK_Exit::ExceptionExit("The number of blocks in boundary file is inconsistent with the number of blocks in grid file!\n");
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
        StructGrid *grid = StructGridCast(structGrids[iZone]);

        int ni, nj, nk;
        ReadNewLine(file, line);
        line = FindNextWord(line, word, separator);
        from_string <int> (ni,word, std::dec);
        line = FindNextWord(line, word, separator);
        from_string <int> (nj, word, std::dec);
        nk = 1;
        if (nDim == THREE_D)
        {
            line = FindNextWord(line, word, separator);
            from_string <int> (nk, word, std::dec);
        }

        cout << "ni,nj,nk = " << ni << " " << nj << " " << nk << "\n";
        int ndif = ABS(ni - grid->GetNI()) + ABS(nj - grid->GetNJ()) + ABS(nk - grid->GetNK()) ;
        if (ndif != 0)
        {
            cout << "No." << iZone + 1 << "The number of IJK-Dimension in boundary file is inconsistent with the number of IJK-Dimension in grid file!\n";
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
                    if ((BCType > 7) && (BCType < 51))
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
                from_string <int> (originalBlockIndex , word, std::dec);
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
        delete [] blockProcessorINP;
        SetZoneProcessID_INP(0);
    }

    file.close();
    file.clear();
}

void ReadPartitionFileForOverset(int &inputNumberOfBlocks, const string &PartitionFile)
{
    fstream file;
    file.open(PartitionFile.c_str(), ios_base::in);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(PartitionFile);
    }

    string line, word;
    string separator = " =\t\r\n#$,;";

    string str1;
    ReadNewLine(file, line);
    line = FindNextWord(line, str1, separator);

    ReadNewLine(file, line);
    line = FindNextWord(line, word, separator);
    int numberOfBlocks;
    from_string <int> (numberOfBlocks, word, std::dec);
    cout << "number_of_blocks = " << numberOfBlocks << "\n";

    if (numberOfBlocks != inputNumberOfBlocks)
    {
        TK_Exit::ExceptionExit(" the numberOfBlocks of partition files does not agree with in mesh file\n");
    }

    ReadNewLine(file, line);    //! Read the third line.

    bool blockProcessFlag = true;

    using namespace PHMPI;

    CreatetZoneProcessorID_INP(numberOfBlocks);
    CreateZoneGroupID(numberOfBlocks);
    int *blockProcessorINP = GetZoneProcessorID_INP();
    int *blockGroupID = GetZoneGroupID();

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        int zone, rank, group;
        ReadNewLine(file, line);
        line = FindNextWord(line, word, separator);
        from_string <int> (zone, word, std::dec);
        line = FindNextWord(line, word, separator);
        from_string <int> (rank, word, std::dec);
        line = FindNextWord(line, word, separator);
        from_string <int> (group, word, std::dec);

        cout << "zone, rank, group = " << zone << " " << rank << " " << group << "\n";

        if (iZone != zone)
        {
            cout << "partition.ppsh NO." << iZone + 1 << "block data error!\n";
        }

        blockProcessorINP[iZone] = rank;
        blockGroupID[iZone] = group;
    }

    if (!blockProcessFlag)
    {
        delete [] blockProcessorINP;
        delete [] blockGroupID;
        SetZoneProcessID_INP(0);
    }

    file.close();
    file.clear();
}

}
