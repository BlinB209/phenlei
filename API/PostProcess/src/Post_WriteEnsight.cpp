#include "Post_WriteEnsight.h"
#include "PHIO.h"
#include "TK_Exit.h"
#include "TK_Parse.h"
#include "TK_Log.h"
#include "IO_FileName.h"
#include "Glb_Dimension.h"
#pragma warning(disable : 4996)

namespace PHSPACE
{

Post_WriteEnsight::Post_WriteEnsight()
{
    this->dataOnElement      = false;
    this->writeFieldEachIter = false;

    int nProcessor = PHMPI::GetNumberOfProcessor();

    if (nProcessor > 1)
    {
        fileMethod = MPI_PARALLEL_IO;
    }
    else
    {
        fileMethod = STDIO_SERIAL;
    }

    int visualfileType = GlobalDataBase::GetIntParaFromDB("visualfileType");
    if (visualfileType == Ensight_Binary || fileMethod)
    {
        fileType = 0;
    }
    else
    {
        fileType = 1;
    }
}

Post_WriteEnsight::~Post_WriteEnsight()
{

}

void Post_WriteEnsight::Initialize()
{
    string visualfile = "tecflow.dat";
    GlobalDataBase::GetData("visualfile", &visualfile, PHSTRING, 1);

    visualFileName = ChangeExtensionOfFileName(visualfile, "case");;
    if (flowType == AverageFlow)
    {
        visualFileName = AddSymbolToFileName(visualfile,"_Average");
    }

    if (flowType == AverageReynoldsStress)
    {
        visualFileName = AddSymbolToFileName(visualfile,"_ReynoldsStress");
    }

    string sentinelfilename = "results/sentinel";
    fstream sentinelfile;
    sentinelfile.open(sentinelfilename.c_str(), ios::in);
    VTKvisual = false;
    if (!sentinelfile && !visualfile.empty() && (flowType == 0))
    {
        VTKFileName = visualfile + ".bak";
        VTKvisual = true;
    }

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady)
    {
        int outnstep         = GlobalDataBase::GetIntParaFromDB("outnstep");
        int intervalStepPlot = GlobalDataBase::GetIntParaFromDB("intervalStepPlot");

        this->currentIter = outnstep / intervalStepPlot;
        stringstream ss;
        ss << setw(6) << setfill('0') << currentIter;
        string str;
        ss >> str;

        writeFieldEachIter = true;
        this->currentIterStr = str;
    }

    meshFileName = ChangeExtensionOfFileName(visualfile, "geo");
}

void Post_WriteEnsight::WriteFile()
{
    int intervalStepPlot = GlobalDataBase::GetIntParaFromDB("intervalStepPlot");
    int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (((outIterStep + 1) % intervalStepPlot == 0) && (outIterStep % intervalStepPlot != 0))
    {
        return;
    }

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    if (VTKvisual == true && myid == server)
    {
        bool CharacteristicBoundary = false;
        DumpToVTK DumpToVTK(VTKDataList, VTKFileName, CharacteristicBoundary, this->flowType);
        DumpToVTK.Run();

        WriteSentinelFile();
    }

    WriteGeometryFile();

    WriteFieldData();

    UpdateCaseFile();
}

void Post_WriteEnsight::WriteGeometryFile()
{
    if (!gridChanged)
    {
        return;
    }

    GetFileHandle(meshFileName.c_str());

    if (!fileType)
    {
        WriteString("C Binary");
    }

    WriteString("Geometry File");
    WriteString("Output by Ensight");
    WriteString("node id assign");
    WriteString("element id assign");

    if (fileMethod > STDIO_SERIAL)
    {
        WriteGeometryFileParallel();
    }
    else
    {
        WriteGeometryFileSerial();
    }

    ClearFileHandle();
    this->gridChanged = false;
}

void Post_WriteEnsight::WriteGeometryFileSerial()
{
    int ipart  = 1;

    //! Write boundary grid.
    for (int iBoundary = 0; iBoundary < boundaryGrid.size(); ++ iBoundary)
    {
        Grid *grid        = boundaryGrid[iBoundary];
        int  oriGridIndex = grid->GetZoneID();

        ostringstream ossBlockName;
        ossBlockName << "Zone" << oriGridIndex << "_" << boundaryName[iBoundary];
        string blockName = ossBlockName.str();

        if (grid->Type() == UNSTRUCTGRID)
        {
            UnstructGrid *grid_uns      = UnstructGridCast(grid);
            int           numberOfNodes = grid_uns->GetNTotalNode();
            int           numberOfCells = grid_uns->GetNTotalCell();
            RDouble      *x             = grid_uns->GetX();
            RDouble      *y             = grid_uns->GetY();
            RDouble      *z             = grid_uns->GetZ();

            int *cell2node              = grid_uns->GetCell2Node();
            int *nodeNumberOfEachCell   = grid_uns->GetNodeNumberOfEachCell();

            if (fileType)
            {
                fprintf(serialFilehandle, "%s\n", "part");
                fprintf(serialFilehandle, "%10d\n", ipart);
                fprintf(serialFilehandle, "%s\n", blockName.c_str());
                ipart += 1;

                fprintf(serialFilehandle, "%s\n", "coordinates");
                fprintf(serialFilehandle, "%10d\n", numberOfNodes);

                for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                {
                    float vertex_buf = static_cast<float>(x[iNode]);
                    fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                }
                for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                {
                    float vertex_buf = static_cast<float>(y[iNode]);
                    fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                }
                for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                {
                    float vertex_buf = static_cast<float>(z[iNode]);
                    fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                }

                fprintf(serialFilehandle, "%s\n", "quad4");
                fprintf(serialFilehandle, "%10d\n", numberOfCells);

                int subscriptIndex = 0;
                for (int iCell = 0; iCell < numberOfCells; ++ iCell)
                {
                    for (int iNode = 0; iNode < nodeNumberOfEachCell[iCell]; ++ iNode)
                    {
                        fprintf(serialFilehandle, "%10d", cell2node[subscriptIndex ++] + 1);
                    }
                    fprintf(serialFilehandle, "\n");
                }
            }
            else
            {
                WriteString("part");
                WriteInt(ipart);
                ipart += 1;

                WriteString(blockName.c_str());
                WriteString("coordinates");

                WriteInt(numberOfNodes);
                for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                {
                    float vertex_buf = static_cast<float>(x[iNode]);
                    fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                }
                for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                {
                    float vertex_buf = static_cast<float>(y[iNode]);
                    fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                }
                for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                {
                    float vertex_buf = static_cast<float>(z[iNode]);
                    fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                }

                WriteString("quad4");
                WriteInt(numberOfCells);

                vector<int> buf(4);
                int subscriptIndex = 0;
                for (int iCell = 0; iCell < numberOfCells; ++ iCell)
                {
                    for (int iNode = 0; iNode < nodeNumberOfEachCell[iCell]; ++ iNode)
                    {
                        buf[iNode] = cell2node[subscriptIndex ++] + 1;
                    }

                    fwrite(buf.data(), sizeof(int), nodeNumberOfEachCell[iCell], serialFilehandle);
                }
            }
        }
        else
        {
            StructGrid *grid_str   = StructGridCast(grid);
            int         nTotalNode = grid_str->GetNTotalNode();

            int     ni = grid_str->GetNI();
            int     nj = grid_str->GetNJ();
            int     nk = grid_str->GetNK();
            RDouble *x = grid_str->GetX();
            RDouble *y = grid_str->GetY();
            RDouble *z = grid_str->GetZ();

            if (fileType)
            {
                fprintf(serialFilehandle, "%s\n", "part");
                fprintf(serialFilehandle, "%10d\n", ipart);
                fprintf(serialFilehandle, "%s\n", blockName.c_str());
                ipart += 1;

                fprintf(serialFilehandle, "%s\n", "block");
                fprintf(serialFilehandle, "%10d%10d%10d\n", ni, nj, nk);
                for (int iNode = 0; iNode < nTotalNode; ++iNode)
                {
                    float vertex_buf = static_cast<float>(x[iNode]);
                    fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                }
                for (int iNode = 0; iNode < nTotalNode; ++iNode)
                {
                    float vertex_buf = static_cast<float>(y[iNode]);
                    fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                }
                for (int iNode = 0; iNode < nTotalNode; ++iNode)
                {
                    float vertex_buf = static_cast<float>(z[iNode]);
                    fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                }
            }
            else
            {
                WriteString("part");
                WriteInt(ipart);
                ipart += 1;

                WriteString(blockName.c_str());

                WriteString("block");
                WriteInt(ni);
                WriteInt(nj);
                WriteInt(nk);

                for (int iNode = 0; iNode < nTotalNode; ++iNode)
                {
                    float vertex_buf = static_cast<float>(x[iNode]);
                    fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                }
                for (int iNode = 0; iNode < nTotalNode; ++iNode)
                {
                    float vertex_buf = static_cast<float>(y[iNode]);
                    fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                }
                for (int iNode = 0; iNode < nTotalNode; ++iNode)
                {
                    float vertex_buf = static_cast<float>(z[iNode]);
                    fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                }
            }
        }
    }


    //! Write block grid.
    int nZones = PHMPI::GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Grid *grid = flowFieldDataOnNode[iZone]->GetGrid();
        if (!WantVisualField(grid))
        {
            continue;
        }

        int gridIndex = grid->GetZoneID();
        ostringstream ossBlockName;
        ossBlockName << "Zone" << gridIndex << "_" << "BLK";
        string blockName = ossBlockName.str();

        if (grid->Type() == UNSTRUCTGRID)
        {
            UnstructGrid *grid_uns = UnstructGridCast(grid);

            int   numberOfNodes           = grid_uns->GetNTotalNode();
            int   numberOfCells           = grid_uns->GetNTotalCell();
            int  *keyActiveOfCells        = grid_uns->GetBlankIndex();
            int **cellNodeIndexContainer  = grid_uns->GetCell2NodeArray();
            int  *cellNodeNumberContainer = grid_uns->GetNodeNumberOfEachCell();
            RDouble *x                    = grid_uns->GetX();
            RDouble *y                    = grid_uns->GetY();
            RDouble *z                    = grid_uns->GetZ();

            if (grid_uns->GetDim() != THREE_D)
            {
                vector <int> cell_type_count(2, 0);
                for (int iCell = 0; iCell < numberOfCells; ++iCell)
                {
                    if (keyActiveOfCells[iCell] == ACTIVE)
                    {
                        if (cellNodeNumberContainer[iCell] == 3)
                            cell_type_count[0] += 1;
                        else
                            cell_type_count[1] += 1;
                    }
                }

                vector < vector<int> > cell_type_list(2);
                for (int i = 0; i < 2; ++i)
                {
                    cell_type_list[i].resize(cell_type_count[i]);
                    cell_type_count[i] = 0;
                }

                for (int iCell = 0; iCell < numberOfCells; ++iCell)
                {
                    if (keyActiveOfCells[iCell] == ACTIVE)
                    {
                        if (cellNodeNumberContainer[iCell] == 3)
                        {
                            cell_type_list[0][cell_type_count[0]] = iCell;
                            cell_type_count[0] += 1;
                        }
                        else
                        {
                            cell_type_list[1][cell_type_count[1]] = iCell;
                            cell_type_count[1] += 1;
                        }
                    }
                }

                if (fileType)
                {
                    fprintf(serialFilehandle, "%s\n", "part");
                    fprintf(serialFilehandle, "%10d\n", ipart);
                    fprintf(serialFilehandle, "%s\n", blockName.c_str());
                    ipart += 1;

                    fprintf(serialFilehandle, "%s\n", "coordinates");
                    fprintf(serialFilehandle, "%10d\n", numberOfNodes);
                    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                    {
                        float vertex_buf = static_cast<float>(x[iNode]);
                        fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                    }
                    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                    {
                        float vertex_buf = static_cast<float>(y[iNode]);
                        fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                    }
                    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                    {
                        float vertex_buf = static_cast<float>(z[iNode]);
                        fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                    }

                    vector<string> cell_type_name;
                    cell_type_name.push_back("tria3");
                    cell_type_name.push_back("quad4");
                    for (int i = 0; i < 2; ++i)
                    {
                        if (cell_type_count[i] > 0)
                        {
                            fprintf(serialFilehandle, "%s\n", cell_type_name[i].c_str());
                            fprintf(serialFilehandle, "%10d\n", cell_type_count[i]);

                            for (int iCell = 0; iCell < cell_type_count[i]; ++ iCell)
                            {
                                for (int iNode = 0; iNode < cellNodeNumberContainer[cell_type_list[i][iCell]]; ++ iNode)
                                {
                                    fprintf(serialFilehandle, "%10d", cellNodeIndexContainer[cell_type_list[i][iCell]][iNode] + 1);
                                }
                                fprintf(serialFilehandle, "\n");
                            }
                        }
                    }
                }
                else
                {
                    WriteString("part");
                    WriteInt(ipart);
                    ipart += 1;

                    WriteString(blockName.c_str());
                    WriteString("coordinates");
                    WriteInt(numberOfNodes);
                    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                    {
                        float vertex_buf = static_cast<float>(x[iNode]);
                        fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                    }
                    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                    {
                        float vertex_buf = static_cast<float>(y[iNode]);
                        fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                    }
                    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                    {
                        float vertex_buf = static_cast<float>(z[iNode]);
                        fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                    }

                    vector<string> cell_type_name;
                    cell_type_name.push_back("tria3");
                    cell_type_name.push_back("quad4");
                    vector<int> buf(4);
                    for (int i = 0; i < 2; ++i)
                    {
                        if (cell_type_count[i] > 0)
                        {
                            WriteString(cell_type_name[i].c_str());
                            WriteInt(cell_type_count[i]);
                            for (int iCell = 0; iCell < cell_type_count[i]; ++ iCell)
                            {
                                for (int iNode = 0; iNode < cellNodeNumberContainer[cell_type_list[i][iCell]]; ++ iNode)
                                {
                                    buf[iNode] = cellNodeIndexContainer[cell_type_list[i][iCell]][iNode] + 1;
                                }

                                fwrite(buf.data(), sizeof(int), cellNodeNumberContainer[cell_type_list[i][iCell]], serialFilehandle);
                            }
                        }
                    }
                }
            }
            else
            {
                vector<int> cell_type_count(4, 0);
                for (int iCell = 0; iCell < numberOfCells; ++ iCell)
                {
                    if (keyActiveOfCells[iCell] == ACTIVE)
                    {
                        if (cellNodeNumberContainer[iCell] == 4)
                            cell_type_count[0] += 1;
                        else if (cellNodeNumberContainer[iCell] == 5)
                            cell_type_count[1] += 1;
                        else if (cellNodeNumberContainer[iCell] == 6)
                            cell_type_count[2] += 1;
                        else if (cellNodeNumberContainer[iCell] == 8)
                            cell_type_count[3] += 1;
                        else
                            continue;
                    }
                }

                vector< vector<int> > cell_type_list(4);
                for (int i = 0; i < 4; ++ i)
                {
                    cell_type_list[i].resize(cell_type_count[i]);
                    cell_type_count[i] = 0;
                }

                for (int iCell = 0; iCell < numberOfCells; ++ iCell)
                {
                    if (keyActiveOfCells[iCell] == ACTIVE)
                    {
                        if (cellNodeNumberContainer[iCell] == 4)
                        {
                            cell_type_list[0][cell_type_count[0]] = iCell;
                            cell_type_count[0] += 1;
                        }
                        else if (cellNodeNumberContainer[iCell] == 5)
                        {
                            cell_type_list[1][cell_type_count[1]] = iCell;
                            cell_type_count[1] += 1;
                        }
                        else if (cellNodeNumberContainer[iCell] == 6)
                        {
                            cell_type_list[2][cell_type_count[2]] = iCell;
                            cell_type_count[2] += 1;
                        }
                        else if (cellNodeNumberContainer[iCell] == 8)
                        {
                            cell_type_list[3][cell_type_count[3]] = iCell;
                            cell_type_count[3] += 1;
                        }
                        else
                            continue;
                    }
                }

                string element[4] = {"tetra4", "pyramid5", "penta6", "hexa8"};
                if (fileType)
                {
                    fprintf(serialFilehandle, "%s\n", "part");
                    fprintf(serialFilehandle, "%10d\n", ipart);
                    fprintf(serialFilehandle, "%s\n", blockName.c_str());
                    ipart += 1;

                    fprintf(serialFilehandle, "%s\n", "coordinates");
                    fprintf(serialFilehandle, "%10d\n", numberOfNodes);
                    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
                    {
                        float vertex_buf = static_cast<float>(x[iNode]);
                        fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                    }
                    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
                    {
                        float vertex_buf = static_cast<float>(y[iNode]);
                        fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                    }
                    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
                    {
                        float vertex_buf = static_cast<float>(z[iNode]);
                        fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                    }

                    for (int i = 0; i < 4; ++ i)
                    {
                        if (cell_type_count[i] > 0)
                        {
                            fprintf(serialFilehandle, "%s\n", element[i].c_str());
                            fprintf(serialFilehandle, "%10d\n", cell_type_count[i]);

                            for (int iCell = 0; iCell < cell_type_count[i]; ++ iCell)
                            {
                                for (int iNode = 0; iNode < cellNodeNumberContainer[cell_type_list[i][iCell]]; ++ iNode)
                                {
                                    fprintf(serialFilehandle, "%10d", cellNodeIndexContainer[cell_type_list[i][iCell]][iNode] + 1);
                                }
                                fprintf(serialFilehandle, "\n");
                            }
                        }
                    }
                }
                else
                {
                    WriteString("part");
                    WriteInt(ipart);
                    ipart += 1;

                    WriteString(blockName.c_str());

                    WriteString("coordinates");
                    WriteInt(numberOfNodes);

                    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
                    {
                        float vertex_buf = static_cast<float>(x[iNode]);
                        fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                    }
                    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
                    {
                        float vertex_buf = static_cast<float>(y[iNode]);
                        fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                    }
                    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
                    {
                        float vertex_buf = static_cast<float>(z[iNode]);
                        fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                    }

                    vector<int> buf(8);
                    for (int i = 0; i < 4; ++ i)
                    {
                        if (cell_type_count[i] > 0)
                        {
                            WriteString(element[i].c_str());
                            WriteInt(cell_type_count[i]);
                            for (int iCell = 0; iCell < cell_type_count[i]; ++ iCell)
                            {
                                for (int iNode = 0; iNode < cellNodeNumberContainer[cell_type_list[i][iCell]]; ++ iNode)
                                {
                                    buf[iNode] = cellNodeIndexContainer[cell_type_list[i][iCell]][iNode] + 1;
                                }

                                fwrite(buf.data(), sizeof(int), cellNodeNumberContainer[cell_type_list[i][iCell]], serialFilehandle);
                            }
                        }
                    }
                }
            }
        }
        else
        {
            StructGrid *grid_str   = StructGridCast(grid);
            int         nTotalNode = grid_str->GetNTotalNode();

            int     ni = grid_str->GetNI();
            int     nj = grid_str->GetNJ();
            int     nk = grid_str->GetNK();
            RDouble *x = grid_str->GetX();
            RDouble *y = grid_str->GetY();
            RDouble *z = grid_str->GetZ();

            if (fileType)
            {
                fprintf(serialFilehandle, "%s\n", "part");
                fprintf(serialFilehandle, "%10d\n", ipart);
                fprintf(serialFilehandle, "%s\n", blockName.c_str());
                ipart += 1;

                fprintf(serialFilehandle, "%s\n", "block");
                fprintf(serialFilehandle, "%10d%10d%10d\n", ni, nj, nk);
                for (int iNode = 0; iNode < nTotalNode; ++iNode)
                {
                    float vertex_buf = static_cast<float>(x[iNode]);
                    fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                }
                for (int iNode = 0; iNode < nTotalNode; ++iNode)
                {
                    float vertex_buf = static_cast<float>(y[iNode]);
                    fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                }
                for (int iNode = 0; iNode < nTotalNode; ++iNode)
                {
                    float vertex_buf = static_cast<float>(z[iNode]);
                    fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                }
            }
            else
            {
                WriteString("part");
                WriteInt(ipart);
                ipart += 1;

                WriteString(blockName.c_str());

                WriteString("block");
                WriteInt(ni);
                WriteInt(nj);
                WriteInt(nk);

                for (int iNode = 0; iNode < nTotalNode; ++iNode)
                {
                    float vertex_buf = static_cast<float>(x[iNode]);
                    fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                }
                for (int iNode = 0; iNode < nTotalNode; ++iNode)
                {
                    float vertex_buf = static_cast<float>(y[iNode]);
                    fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                }
                for (int iNode = 0; iNode < nTotalNode; ++iNode)
                {
                    float vertex_buf = static_cast<float>(z[iNode]);
                    fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                }
            }
        }
    }
}

void Post_WriteEnsight::WriteGeometryFileParallel()
{
    using namespace PHMPI;
    MPI_Offset offset_my_proc = 0;
    int        my_zone_num    = 0;

    //! Colloect data length.
    for (int iBoundary = 0; iBoundary < boundaryGrid.size(); ++ iBoundary)
    {
        my_zone_num += 1;
        Grid *grid = boundaryGrid[iBoundary];
        if (grid->Type() == UNSTRUCTGRID)
        {
            UnstructGrid *grid_uns          = UnstructGridCast(grid);
            int           numberOfNodes           = grid_uns->GetNTotalNode();
            int           numberOfCells           = grid_uns->GetNTotalCell();

            offset_my_proc += 80 + sizeof(int) + 80 
                           +  80 + sizeof(int) + numberOfNodes * sizeof(float) * 3
                           +  80 + sizeof(int) + numberOfCells * 4 * sizeof(int);
        }
        else
        {
            StructGrid *grid_str   = StructGridCast(grid);
            int         nTotalNode = grid_str->GetNTotalNode();
            offset_my_proc += 80 + sizeof(int) + 80 + 80 + sizeof(int) * 3 +
                                      nTotalNode * sizeof(float) * 3;
        }
    }

    int myid   = GetCurrentProcessorID();
    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int proc_id = GetZoneProcessorID(iZone);
        if (myid != proc_id)
        {
            continue;
        }

        Grid *grid = flowFieldDataOnNode[iZone]->GetGrid();
        if (!WantVisualField(grid))
        {
            continue;
        }

        my_zone_num += 1;
        if (grid->Type() == UNSTRUCTGRID)
        {
            UnstructGrid *grid_uns                = UnstructGridCast(grid);
            int           numberOfNodes           = grid_uns->GetNTotalNode();
            int           numberOfCells           = grid_uns->GetNTotalCell();
            int          *keyActiveOfCells        = grid_uns->GetBlankIndex();
            int          *cellNodeNumberContainer = grid_uns->GetNodeNumberOfEachCell();
            if (grid_uns->GetDim() != THREE_D)
            {
                std::vector<int> cell_type_count(2, 0);
                for (int iCell = 0; iCell < numberOfCells; ++iCell)
                {
                    if (keyActiveOfCells[iCell] == ACTIVE)
                    {
                        if (cellNodeNumberContainer[iCell] == 3)
                            cell_type_count[0] += 1;
                        else
                            cell_type_count[1] += 1;
                    }
                }

                offset_my_proc += 80 + sizeof(int) + 80 + 80 + sizeof(int) +
                    numberOfNodes * sizeof(float) * 3 +
                    (cell_type_count[0] > 0 ? 1 : 0) * (80 + sizeof(int)) +
                    cell_type_count[0] * 3 * sizeof(int) +
                    (cell_type_count[1] > 0 ? 1 : 0) * (80 + sizeof(int)) +
                    cell_type_count[1] * 4 * sizeof(int);
            }
            else if (grid_uns->GetDim() == THREE_D)
            {
                std::vector<int> cell_type_count(4, 0);
                for (int iCell = 0; iCell < numberOfCells; ++iCell)
                {
                    if (keyActiveOfCells[iCell] == ACTIVE)
                    {
                        if (cellNodeNumberContainer[iCell] == 4)
                            cell_type_count[0] += 1;
                        else if (cellNodeNumberContainer[iCell] == 5)
                            cell_type_count[1] += 1;
                        else if (cellNodeNumberContainer[iCell] == 6)
                            cell_type_count[2] += 1;
                        else if (cellNodeNumberContainer[iCell] == 8)
                            cell_type_count[3] += 1;
                        else
                            continue;
                    }
                }
                offset_my_proc += 80 + sizeof(int) + 80 + 80 + sizeof(int) +
                    numberOfNodes * sizeof(float) * 3 +
                    (cell_type_count[0] > 0 ? 1 : 0) * (80 + sizeof(int)) +
                    cell_type_count[0] * 4 * sizeof(int) +
                    (cell_type_count[1] > 0 ? 1 : 0) * (80 + sizeof(int)) +
                    cell_type_count[1] * 5 * sizeof(int) +
                    (cell_type_count[2] > 0 ? 1 : 0) * (80 + sizeof(int)) +
                    cell_type_count[2] * 6 * sizeof(int) +
                    (cell_type_count[3] > 0 ? 1 : 0) * (80 + sizeof(int)) +
                    cell_type_count[3] * 8 * sizeof(int);
            }
        }
        else
        {
            StructGrid *grid_str   = StructGridCast(grid);
            int         nTotalNode = grid_str->GetNTotalNode();
            offset_my_proc += 80 + sizeof(int) + 80 + 80 + sizeof(int) * 3 +
                nTotalNode * sizeof(float) * 3;
        }
    }


    //! mpi scan total offset
    MPI_Offset offset_write;
    MPI_Scan(&offset_my_proc, &offset_write, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    offset_write = offset_write - offset_my_proc + fileOffset;

    int currentZone;
    MPI_Scan(&my_zone_num, &currentZone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    currentZone = currentZone - my_zone_num + 1;


    //! Write boundary grid.
    for (int iBoundary = 0; iBoundary < boundaryGrid.size(); ++ iBoundary)
    {
        Grid *grid = boundaryGrid[iBoundary];
        int  oriGridIndex = grid->GetZoneID();

        ostringstream ossBlockName;
        ossBlockName << "Zone" << oriGridIndex << "_" << boundaryName[iBoundary];
        string blockName = ossBlockName.str();

        if (grid->Type() == UNSTRUCTGRID)
        {
            UnstructGrid *grid_uns      = UnstructGridCast(grid);
            int           numberOfNodes = grid_uns->GetNTotalNode();
            int           numberOfCells = grid_uns->GetNTotalCell();
            RDouble      *x             = grid_uns->GetX();
            RDouble      *y             = grid_uns->GetY();
            RDouble      *z             = grid_uns->GetZ();

            int *cell2node              = grid_uns->GetCell2Node();
            int *nodeNumberOfEachCell   = grid_uns->GetNodeNumberOfEachCell();

            MpiWriteString("part", offset_write);
            offset_write += 80;

            MpiWriteInt(currentZone, offset_write);
            offset_write += sizeof(int);
            currentZone  += 1;

            MpiWriteString(blockName.c_str(), offset_write);
            offset_write += 80;

            MpiWriteString("coordinates", offset_write);
            offset_write += 80;

            MpiWriteInt(numberOfNodes, offset_write);
            offset_write += sizeof(int);

            MPI_Status status;
            for (int iNode = 0; iNode < numberOfNodes; ++iNode)
            {
                float      vertex_buf  = static_cast<float>(x[iNode]);
                MPI_Offset offset_node = offset_write + iNode * sizeof(float);
                MPI_File_write_at(mpiFilehandle, offset_node, &vertex_buf, sizeof(float), MPI_BYTE, &status);
            }
            for (int iNode = 0; iNode < numberOfNodes; ++iNode)
            {
                float      vertex_buf  = static_cast<float>(y[iNode]);
                MPI_Offset offset_node = offset_write + (numberOfNodes + iNode) * sizeof(float);
                MPI_File_write_at(mpiFilehandle, offset_node, &vertex_buf, sizeof(float), MPI_BYTE, &status);
            }
            for (int iNode = 0; iNode < numberOfNodes; ++iNode)
            {
                float      vertex_buf  = static_cast<float>(z[iNode]);
                MPI_Offset offset_node = offset_write + (numberOfNodes * 2 + iNode) * sizeof(float);
                MPI_File_write_at(mpiFilehandle, offset_node, &vertex_buf, sizeof(float), MPI_BYTE, &status);
            }
            offset_write += numberOfNodes * 3 * sizeof(float);

            MpiWriteString("quad4", offset_write);
            offset_write += 80;
            MpiWriteInt(numberOfCells, offset_write);
            offset_write += sizeof(int);

            vector<int> buf(4);
            int subscriptIndex = 0;
            for (int iCell = 0; iCell < numberOfCells; ++iCell)
            {
                for (int iNode = 0; iNode < nodeNumberOfEachCell[iCell]; ++ iNode)
                {
                    buf[iNode] = cell2node[subscriptIndex ++] + 1;
                }

                MPI_Offset offset_cell = offset_write + iCell * sizeof(int) * nodeNumberOfEachCell[iCell];
                MPI_File_write_at(mpiFilehandle, offset_cell, buf.data(), sizeof(int) * nodeNumberOfEachCell[iCell], MPI_BYTE, &status);
            }
            offset_write += numberOfCells * 4 * sizeof(int);
        }
        else
        {
            StructGrid *grid_str   = StructGridCast(grid);
            int         nTotalNode = grid_str->GetNTotalNode();

            int      ni = grid_str->GetNI();
            int      nj = grid_str->GetNJ();
            int      nk = grid_str->GetNK();
            RDouble *x  = grid_str->GetX();
            RDouble *y  = grid_str->GetY();
            RDouble *z  = grid_str->GetZ();

            MpiWriteString("part", offset_write);
            offset_write += 80;

            MpiWriteInt(currentZone, offset_write);
            offset_write += sizeof(int);
            currentZone  += 1;

            MpiWriteString(blockName.c_str(), offset_write);
            offset_write += 80;

            MpiWriteString("block", offset_write);
            offset_write += 80;

            MpiWriteInt(ni, offset_write);
            offset_write += sizeof(int);
            MpiWriteInt(nj, offset_write);
            offset_write += sizeof(int);
            MpiWriteInt(nk, offset_write);
            offset_write += sizeof(int);

            MPI_Status status;
            for (int iNode = 0; iNode < nTotalNode; ++iNode)
            {
                float      vertex_buf  = static_cast<float>(x[iNode]);
                MPI_Offset offset_node = offset_write + iNode * sizeof(float);
                MPI_File_write_at(mpiFilehandle, offset_node, &vertex_buf, sizeof(float), MPI_BYTE, &status);
            }
            for (int iNode = 0; iNode < nTotalNode; ++iNode)
            {
                float      vertex_buf  = static_cast<float>(y[iNode]);
                MPI_Offset offset_node = offset_write + (nTotalNode + iNode) * sizeof(float);
                MPI_File_write_at(mpiFilehandle, offset_node, &vertex_buf, sizeof(float), MPI_BYTE, &status);
            }
            for (int iNode = 0; iNode < nTotalNode; ++iNode)
            {
                float      vertex_buf  = static_cast<float>(z[iNode]);
                MPI_Offset offset_node = offset_write + (nTotalNode * 2 + iNode) * sizeof(float);
                MPI_File_write_at(mpiFilehandle, offset_node, &vertex_buf, sizeof(float), MPI_BYTE, &status);
            }
            offset_write += nTotalNode * 3 * sizeof(float);
        }
    }


    //! Write block grid.
    for (int iZone = 0; iZone < flowFieldDataOnNode.size(); ++ iZone)
    {
        int proc_id = GetZoneProcessorID(iZone);
        if (myid != proc_id)
        {
            continue;
        }

        Grid *grid = flowFieldDataOnNode[iZone]->GetGrid();
        if (!WantVisualField(grid))
        {
            continue;
        }

        int gridIndex = grid->GetZoneID();
        ostringstream ossBlockName;
        ossBlockName << "Zone" << gridIndex << "_" << "BLK";
        string blockName = ossBlockName.str();

        if (grid->Type() == UNSTRUCTGRID)
        {
            UnstructGrid *grid_uns = UnstructGridCast(grid);

            int   numberOfNodes           = grid_uns->GetNTotalNode();
            int   numberOfCells           = grid_uns->GetNTotalCell();
            int  *keyActiveOfCells        = grid_uns->GetBlankIndex();
            int **cellNodeIndexContainer  = grid_uns->GetCell2NodeArray();
            int  *cellNodeNumberContainer = grid_uns->GetNodeNumberOfEachCell();
            RDouble *x                    = grid_uns->GetX();
            RDouble *y                    = grid_uns->GetY();
            RDouble *z                    = grid_uns->GetZ();

            if (grid_uns->GetDim() != THREE_D)
            {
                vector <int> cell_type_count(2, 0);
                for (int iCell = 0; iCell < numberOfCells; ++iCell)
                {
                    if (keyActiveOfCells[iCell] == ACTIVE)
                    {
                        if (cellNodeNumberContainer[iCell] == 3)
                            cell_type_count[0] += 1;
                        else
                            cell_type_count[1] += 1;
                    }
                }

                vector < vector<int> > cell_type_list(2);
                for (int i = 0; i < 2; ++i)
                {
                    cell_type_list[i].resize(cell_type_count[i]);
                    cell_type_count[i] = 0;
                }

                for (int iCell = 0; iCell < numberOfCells; ++iCell)
                {
                    if (keyActiveOfCells[iCell] == ACTIVE)
                    {
                        if (cellNodeNumberContainer[iCell] == 3)
                        {
                            cell_type_list[0][cell_type_count[0]] = iCell;
                            cell_type_count[0] += 1;
                        }
                        else
                        {
                            cell_type_list[1][cell_type_count[1]] = iCell;
                            cell_type_count[1] += 1;
                        }
                    }
                }

                MpiWriteString("part", offset_write);
                offset_write += 80;

                MpiWriteInt(currentZone, offset_write);
                offset_write += sizeof(int);
                currentZone  += 1;

                MpiWriteString(blockName.c_str(), offset_write);
                offset_write += 80;

                MpiWriteString("coordinates", offset_write);
                offset_write += 80;

                MpiWriteInt(numberOfNodes, offset_write);
                offset_write += sizeof(int);

                MPI_Status status;
                for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                {
                    float      vertex_buf  = static_cast<float>(x[iNode]);
                    MPI_Offset offset_node = offset_write + iNode * sizeof(float);
                    MPI_File_write_at(mpiFilehandle, offset_node, &vertex_buf, sizeof(float), MPI_BYTE, &status);
                }
                for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                {
                    float      vertex_buf  = static_cast<float>(y[iNode]);
                    MPI_Offset offset_node = offset_write + (numberOfNodes + iNode) * sizeof(float);
                    MPI_File_write_at(mpiFilehandle, offset_node, &vertex_buf, sizeof(float), MPI_BYTE, &status);
                }
                for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                {
                    float      vertex_buf  = static_cast<float>(z[iNode]);
                    MPI_Offset offset_node = offset_write + (numberOfNodes * 2 + iNode) * sizeof(float);
                    MPI_File_write_at(mpiFilehandle, offset_node, &vertex_buf, sizeof(float), MPI_BYTE, &status);
                }
                offset_write += numberOfNodes * 3 * sizeof(float);

                vector<string> cell_type_name;
                cell_type_name.push_back("tria3");
                cell_type_name.push_back("quad4");

                vector<int> buf(4);
                for (int i = 0; i < 2; ++i)
                {
                    if (cell_type_count[i] > 0)
                    {
                        MpiWriteString(cell_type_name[i].c_str(), offset_write);
                        offset_write += 80;
                        MpiWriteInt(cell_type_count[i], offset_write);
                        offset_write += sizeof(int);

                        for (int iCell = 0; iCell < cell_type_count[i]; ++iCell)
                        {
                            for (int iNode = 0; iNode < cellNodeNumberContainer[cell_type_list[i][iCell]]; ++ iNode)
                            {
                                buf[iNode] = cellNodeIndexContainer[cell_type_list[i][iCell]][iNode] + 1;
                            }

                            MPI_Offset offset_cell = offset_write + iCell * sizeof(int) * cellNodeNumberContainer[cell_type_list[i][iCell]];
                            MPI_File_write_at(mpiFilehandle, offset_cell, buf.data(), sizeof(int) * cellNodeNumberContainer[cell_type_list[i][iCell]], MPI_BYTE, &status);
                        }
                        offset_write += cell_type_count[i] * (i == 0 ? 3 : 4) * sizeof(int);
                    }
                }
            }
            else
            {
                vector<int> cell_type_count(4, 0);
                for (int iCell = 0; iCell < numberOfCells; ++ iCell)
                {
                    if (keyActiveOfCells[iCell] == ACTIVE)
                    {
                        if (cellNodeNumberContainer[iCell] == 4)
                            cell_type_count[0] += 1;
                        else if (cellNodeNumberContainer[iCell] == 5)
                            cell_type_count[1] += 1;
                        else if (cellNodeNumberContainer[iCell] == 6)
                            cell_type_count[2] += 1;
                        else if (cellNodeNumberContainer[iCell] == 8)
                            cell_type_count[3] += 1;
                        else
                            continue;
                    }
                }

                vector< vector<int> > cell_type_list(4);
                for (int i = 0; i < 4; ++ i)
                {
                    cell_type_list[i].resize(cell_type_count[i]);
                    cell_type_count[i] = 0;
                }

                for (int iCell = 0; iCell < numberOfCells; ++ iCell)
                {
                    if (keyActiveOfCells[iCell] == ACTIVE)
                    {
                        if (cellNodeNumberContainer[iCell] == 4)
                        {
                            cell_type_list[0][cell_type_count[0]] = iCell;
                            cell_type_count[0] += 1;
                        }
                        else if (cellNodeNumberContainer[iCell] == 5)
                        {
                            cell_type_list[1][cell_type_count[1]] = iCell;
                            cell_type_count[1] += 1;
                        }
                        else if (cellNodeNumberContainer[iCell] == 6)
                        {
                            cell_type_list[2][cell_type_count[2]] = iCell;
                            cell_type_count[2] += 1;
                        }
                        else if (cellNodeNumberContainer[iCell] == 8)
                        {
                            cell_type_list[3][cell_type_count[3]] = iCell;
                            cell_type_count[3] += 1;
                        }
                        else
                            continue;
                    }
                }

                MpiWriteString("part", offset_write);
                offset_write += 80;

                MpiWriteInt(currentZone, offset_write);
                offset_write += sizeof(int);
                currentZone  += 1;

                MpiWriteString(blockName.c_str(), offset_write);
                offset_write += 80;

                MpiWriteString("coordinates", offset_write);
                offset_write += 80;

                MpiWriteInt(numberOfNodes, offset_write);
                offset_write += sizeof(int);

                MPI_Status status;
                for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                {
                    float      vertex_buf  = static_cast<float>(x[iNode]);
                    MPI_Offset offset_node = offset_write + iNode * sizeof(float);
                    MPI_File_write_at(mpiFilehandle, offset_node, &vertex_buf, sizeof(float), MPI_BYTE, &status);
                }
                for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                {
                    float      vertex_buf  = static_cast<float>(y[iNode]);
                    MPI_Offset offset_node = offset_write + (numberOfNodes + iNode) * sizeof(float);
                    MPI_File_write_at(mpiFilehandle, offset_node, &vertex_buf, sizeof(float), MPI_BYTE, &status);
                }
                for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                {
                    float      vertex_buf  = static_cast<float>(z[iNode]);
                    MPI_Offset offset_node = offset_write + (numberOfNodes * 2 + iNode) * sizeof(float);
                    MPI_File_write_at(mpiFilehandle, offset_node, &vertex_buf, sizeof(float), MPI_BYTE, &status);
                }
                offset_write += numberOfNodes * 3 * sizeof(float);

                string element[4] = {"tetra4", "pyramid5", "penta6", "hexa8"};
                int     element_vertex_num[4] = {4, 5, 6, 8};
                vector<int> buf(8);
                for (int i = 0; i < 4; ++ i)
                {
                    if (cell_type_count[i] > 0)
                    {
                        MpiWriteString(element[i].c_str(), offset_write);
                        offset_write += 80;
                        MpiWriteInt(cell_type_count[i], offset_write);
                        offset_write += sizeof(int);

                        for (int iCell = 0; iCell < cell_type_count[i]; ++ iCell)
                        {
                            for (int iNode = 0; iNode < cellNodeNumberContainer[cell_type_list[i][iCell]]; ++ iNode)
                            {
                                buf[iNode] = cellNodeIndexContainer[cell_type_list[i][iCell]][iNode] + 1;
                            }

                            MPI_Offset offset_cell = offset_write + iCell * sizeof(int) * cellNodeNumberContainer[cell_type_list[i][iCell]];
                            MPI_File_write_at(mpiFilehandle, offset_cell, buf.data(), sizeof(int) * cellNodeNumberContainer[cell_type_list[i][iCell]], MPI_BYTE, &status);
                        }
                        offset_write += cell_type_count[i] * element_vertex_num[i] * sizeof(int);
                    }
                }
            }
        }
        else
        {
            StructGrid *grid_str   = StructGridCast(grid);
            int         nTotalNode = grid_str->GetNTotalNode();

            int     ni = grid_str->GetNI();
            int     nj = grid_str->GetNJ();
            int     nk = grid_str->GetNK();
            RDouble *x = grid_str->GetX();
            RDouble *y = grid_str->GetY();
            RDouble *z = grid_str->GetZ();

            MpiWriteString("part", offset_write);
            offset_write += 80;

            MpiWriteInt(currentZone, offset_write);
            offset_write += sizeof(int);
            currentZone  += 1;

            MpiWriteString(blockName.c_str(), offset_write);
            offset_write += 80;

            MpiWriteString("block", offset_write);
            offset_write += 80;

            MpiWriteInt(ni, offset_write);
            offset_write += sizeof(int);
            MpiWriteInt(nj, offset_write);
            offset_write += sizeof(int);
            MpiWriteInt(nk, offset_write);
            offset_write += sizeof(int);

            MPI_Status status;
            for (int iNode = 0; iNode < nTotalNode; ++iNode)
            {
                float      vertex_buf  = static_cast<float>(x[iNode]);
                MPI_Offset offset_node = offset_write + iNode * sizeof(float);
                MPI_File_write_at(mpiFilehandle, offset_node, &vertex_buf, sizeof(float), MPI_BYTE, &status);
            }
            for (int iNode = 0; iNode < nTotalNode; ++iNode)
            {
                float      vertex_buf  = static_cast<float>(y[iNode]);
                MPI_Offset offset_node = offset_write + (nTotalNode + iNode) * sizeof(float);
                MPI_File_write_at(mpiFilehandle, offset_node, &vertex_buf, sizeof(float), MPI_BYTE, &status);
            }
            for (int iNode = 0; iNode < nTotalNode; ++iNode)
            {
                float      vertex_buf  = static_cast<float>(z[iNode]);
                MPI_Offset offset_node = offset_write + (nTotalNode * 2 + iNode) * sizeof(float);
                MPI_File_write_at(mpiFilehandle, offset_node, &vertex_buf, sizeof(float), MPI_BYTE, &status);
            }
            offset_write += nTotalNode * 3 * sizeof(float);
        }
    }
}

void Post_WriteEnsight::WriteFieldData()
{
    if (fileMethod > STDIO_SERIAL)
    {
        WriteFieldDataParallel();
    }
    else
    {
        WriteFieldDataSerial();
    }
}

void Post_WriteEnsight::WriteFieldDataSerial()
{
    set<int> visualVariables = flowFieldDataOnNode[0]->GetVisualVariables();
    for (set<int>::iterator varIter = visualVariables.begin(); varIter != visualVariables.end(); ++ varIter)
    {
        int     variableType = *varIter;
        string       varName = flowFieldDataOnNode[0]->GetVariableName(variableType);
        string fieldFileName = ChangeExtensionOfFileName(visualFileName, varName);
        if (writeFieldEachIter)
        {
            fieldFileName = fieldFileName + currentIterStr;
        }

        GetFileHandle(fieldFileName.c_str());
        WriteString("field File");

        int ipart  = 1;

        //! Write boundary grid.
        for (int iBoundary = 0; iBoundary < boundaryGrid.size(); ++ iBoundary)
        {
            Grid *grid = boundaryGrid[iBoundary];
            if (grid->Type() == UNSTRUCTGRID)
            {
                UnstructGrid *grid_uns      = UnstructGridCast(grid);
                int           numberOfNodes = grid_uns->GetNTotalNode();
                int           oriGridIndex  = grid->GetZoneID();
                RDouble      *qn            = (RDouble *)flowFieldDataOnNode[oriGridIndex]->GetVisualNodeVarPtr(varName);

                if (!dataOnElement)
                {
                    if (fileType)
                    {
                        fprintf(serialFilehandle, "%s\n", "part");
                        fprintf(serialFilehandle, "%10d\n", ipart);
                        ipart += 1;

                        fprintf(serialFilehandle, "%s\n", "coordinates");
                        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
                        {
                            int originaNodeIndex = nodeMap[iBoundary][0][iNode];
                            fprintf(serialFilehandle, "%12.4E\n", qn[originaNodeIndex]);
                        }
                    }
                    else
                    {
                        WriteString("part");
                        WriteInt(ipart);
                        ipart += 1;

                        WriteString("coordinates");
                        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
                        {
                            int originaNodeIndex = nodeMap[iBoundary][0][iNode];
                            float buf = static_cast<float>(qn[originaNodeIndex]);
                            fwrite(&buf, sizeof(float), 1, serialFilehandle);
                        }
                    }
                }
                else
                {
                    TK_Exit::UnexpectedVarValue("dataOnElement", dataOnElement);
                }
            }
            else
            {
                StructGrid *grid_str     = StructGridCast(grid);
                int         nTotalNode   = grid_str->GetNTotalNode();
                int         oriGridIndex = grid->GetZoneID();
                RDouble4D  *qn           = (RDouble4D *)flowFieldDataOnNode[oriGridIndex]->GetVisualNodeVarPtr(varName);

                if (!dataOnElement)
                {
                    if (fileType)
                    {
                        fprintf(serialFilehandle, "%s\n", "part");
                        fprintf(serialFilehandle, "%10d\n", ipart);
                        ipart += 1;
                        fprintf(serialFilehandle, "%s\n", "block");

                        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
                        {
                            int originaI = nodeMap[iBoundary][0][iNode];
                            int originaJ = nodeMap[iBoundary][1][iNode];
                            int originaK = nodeMap[iBoundary][2][iNode];

                            RDouble vertex_buf = (*qn)(originaI, originaJ, originaK, 0);
                            fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                        }
                    }
                    else
                    {
                        WriteString("part");
                        WriteInt(ipart);
                        ipart += 1;
                        WriteString("block");

                        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
                        {
                            int originaI = nodeMap[iBoundary][0][iNode];
                            int originaJ = nodeMap[iBoundary][1][iNode];
                            int originaK = nodeMap[iBoundary][2][iNode];

                            float vertex_buf = static_cast<float>((*qn)(originaI, originaJ, originaK, 0));
                            fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                        }
                    }
                }
            }
        }


        //! Write block grid.
        for (int iZone = 0; iZone < flowFieldDataOnNode.size(); ++iZone)
        {
            Grid *grid = flowFieldDataOnNode[iZone]->GetGrid();
            if (!WantVisualField(grid))
            {
                continue;
            }

            if (grid->Type() == UNSTRUCTGRID)
            {
                UnstructGrid *grid_uns      = UnstructGridCast(grid);
                int           numberOfNodes = grid_uns->GetNTotalNode();
                RDouble      *qn            = (RDouble *)flowFieldDataOnNode[iZone]->GetVisualNodeVarPtr(varName);

                if (!dataOnElement)
                {
                    if (fileType)
                    {
                        fprintf(serialFilehandle, "%s\n", "part");
                        fprintf(serialFilehandle, "%10d\n", ipart);
                        ipart += 1;

                        fprintf(serialFilehandle, "%s\n", "coordinates");
                        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
                        {
                            fprintf(serialFilehandle, "%12.4E\n", qn[iNode]);
                        }
                    }
                    else
                    {
                        WriteString("part");
                        WriteInt(ipart);
                        ipart += 1;

                        WriteString("coordinates");
                        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
                        {
                            float buf = static_cast<float>(qn[iNode]);
                            fwrite(&buf, sizeof(float), 1, serialFilehandle);
                        }
                    }
                }
                else
                {
                    TK_Exit::UnexpectedVarValue("dataOnElement", dataOnElement);
                }
            }
            else
            {
                StructGrid *grid_str = StructGridCast(grid);
                RDouble4D  *qn       = (RDouble4D *)flowFieldDataOnNode[iZone]->GetVisualNodeVarPtr(varName);

                int ni = grid_str->GetNI();
                int nj = grid_str->GetNJ();
                int nk = grid_str->GetNK();

                if (!dataOnElement)
                {
                    if (fileType)
                    {
                        fprintf(serialFilehandle, "%s\n", "part");
                        fprintf(serialFilehandle, "%10d\n", ipart);
                        ipart += 1;
                        fprintf(serialFilehandle, "%s\n", "block");

                        for (int k = 1; k <= nk; ++ k)
                        {
                            for (int j = 1; j <= nj; ++ j)
                            {
                                for (int i = 1; i <= ni; ++ i)
                                {
                                    RDouble vertex_buf = (*qn)(i, j, k, 0);
                                    fprintf(serialFilehandle, "%12.4E\n", vertex_buf);
                                }
                            }
                        }
                    }
                    else
                    {
                        WriteString("part");
                        WriteInt(ipart);
                        ipart += 1;
                        WriteString("block");

                        for (int k = 1; k <= nk; ++ k)
                        {
                            for (int j = 1; j <= nj; ++ j)
                            {
                                for (int i = 1; i <= ni; ++ i)
                                {
                                    float vertex_buf = static_cast<float>((*qn)(i, j, k, 0));
                                    fwrite(&vertex_buf, sizeof(float), 1, serialFilehandle);
                                }
                            }
                        }
                    }
                }
                else
                {
                    TK_Exit::UnexpectedVarValue("dataOnElement", dataOnElement);
                }
            }
        }

        ClearFileHandle();
    }
}

void Post_WriteEnsight::WriteFieldDataParallel()
{
    //! Colloect data length.
    int nZones = static_cast<int>(flowFieldDataOnNode.size());

    MPI_Offset offset_my_proc = 0;
    int        my_zone_num    = 0;

    for (int iBoundary = 0; iBoundary < boundaryGrid.size(); ++ iBoundary)
    {
        my_zone_num += 1;
        Grid *grid = boundaryGrid[iBoundary];
        if (grid->Type() == UNSTRUCTGRID)
        {
            UnstructGrid *grid_uns      = UnstructGridCast(grid);
            int           numberOfNodes = grid_uns->GetNTotalNode();

            offset_my_proc += 80 + sizeof(int) + 80 + numberOfNodes * sizeof(float);
        }
        else
        {
            StructGrid *grid_str   = StructGridCast(grid);
            int         nTotalNode = grid_str->GetNTotalNode();
            int         nTotalCell = grid_str->GetNTotalCell();

            if (!dataOnElement)
            {
                offset_my_proc += 80 + sizeof(int) + 80 + nTotalNode * sizeof(float);
            }
            else
            {
                offset_my_proc += 80 + sizeof(int) + 80 + nTotalCell * sizeof(float);
            }
        }
    }

    for (int iZone = 0; iZone < nZones; ++iZone)
    {
        if (!flowFieldDataOnNode[iZone])
        {
            continue;
        }

        Grid *grid = flowFieldDataOnNode[iZone]->GetGrid();
        if (!WantVisualField(grid))
        {
            continue;
        }

        my_zone_num += 1;
        if (grid->Type() == UNSTRUCTGRID)
        {
            UnstructGrid *grid_uns                = UnstructGridCast(grid);
            int           numberOfNodes           = grid_uns->GetNTotalNode();

            offset_my_proc += 80 + sizeof(int) + 80 + numberOfNodes * sizeof(float);
        }
        else
        {
            StructGrid *grid_str   = StructGridCast(grid);
            int         nTotalNode = grid_str->GetNTotalNode();
            int         nTotalCell = grid_str->GetNTotalCell();

            if (!dataOnElement)
            {
                offset_my_proc += 80 + sizeof(int) + 80 + nTotalNode * sizeof(float);
            }
            else
            {
                offset_my_proc += 80 + sizeof(int) + 80 + nTotalCell * sizeof(float);
            }
        }
    }

    MPI_Offset offsetStart;
    MPI_Scan(&offset_my_proc, &offsetStart, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    offsetStart = offsetStart - offset_my_proc + 80;
    //offset_write = offset_write - offset_my_proc + fileOffset;

    int zoneStartIndex;
    MPI_Scan(&my_zone_num, &zoneStartIndex, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    zoneStartIndex = zoneStartIndex - my_zone_num + 1;

    Post_Visual *visualData = 0;
    for (int iZone = 0; iZone < flowFieldDataOnNode.size(); ++ iZone)
    {
        if (flowFieldDataOnNode[iZone])
        {
            visualData = flowFieldDataOnNode[iZone];
            break;
        }
    }

    set<int> visualVariables = visualData->GetVisualVariables();
    for (set<int>::iterator varIter = visualVariables.begin(); varIter != visualVariables.end(); ++ varIter)
    {
        int     variableType = *varIter;
        string       varName = visualData->GetVariableName(variableType);
        string fieldFileName = ChangeExtensionOfFileName(visualFileName, varName);
        if (writeFieldEachIter)
        {
            fieldFileName = fieldFileName + currentIterStr;
        }

        GetFileHandle(fieldFileName.c_str());
        WriteString("field File");

        int        currentZone  = zoneStartIndex;
        MPI_Offset offset_write = offsetStart;
        for (int iBoundary = 0; iBoundary < boundaryGrid.size(); ++ iBoundary)
        {
            Grid *grid = boundaryGrid[iBoundary];
            if (grid->Type() == UNSTRUCTGRID)
            {
                UnstructGrid *grid_uns      = UnstructGridCast(grid);
                int           numberOfNodes = grid_uns->GetNTotalNode();
                int           oriGridIndex  = grid->GetZoneID();
                RDouble      *qn            = (RDouble *)flowFieldDataOnNode[oriGridIndex]->GetVisualNodeVarPtr(varName);

                MpiWriteString("part", offset_write);
                offset_write += 80;

                MpiWriteInt(currentZone, offset_write);
                offset_write += sizeof(int);
                currentZone  += 1;

                MpiWriteString("coordinates", offset_write);
                offset_write += 80;

                MPI_Status status;
                for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
                {
                    int originaNodeIndex = nodeMap[iBoundary][0][iNode];
                    float buf = static_cast<float>(qn[originaNodeIndex]);

                    MPI_Offset offset_cell = offset_write + iNode * sizeof(float);
                    MPI_File_write_at(mpiFilehandle, offset_cell, &buf, sizeof(float), MPI_BYTE, &status);
                }
                offset_write += numberOfNodes * sizeof(float);
            }
            else
            {
                StructGrid *grid_str     = StructGridCast(grid);
                int         nTotalNode   = grid_str->GetNTotalNode();
                int         oriGridIndex = grid->GetZoneID();
                RDouble4D  *qn           = (RDouble4D *)flowFieldDataOnNode[oriGridIndex]->GetVisualNodeVarPtr(varName);

                MpiWriteString("part", offset_write);
                offset_write += 80;

                MpiWriteInt(currentZone, offset_write);
                offset_write += sizeof(int);
                currentZone  += 1;

                MpiWriteString("block", offset_write);
                offset_write += 80;

                MPI_Status status;
                for (int iNode = 0; iNode < nTotalNode; ++ iNode)
                {
                    int originaI = nodeMap[iBoundary][0][iNode];
                    int originaJ = nodeMap[iBoundary][1][iNode];
                    int originaK = nodeMap[iBoundary][2][iNode];

                    float vertex_buf = static_cast<float>((*qn)(originaI, originaJ, originaK, 0));

                    MPI_Offset offset_node = offset_write + iNode * sizeof(float);
                    MPI_File_write_at(mpiFilehandle, offset_node, &vertex_buf, sizeof(float), MPI_BYTE, &status);
                }
                offset_write += nTotalNode * sizeof(float);
            }
        }

        for (int iZone = 0; iZone < nZones; ++iZone)
        {
            if (!flowFieldDataOnNode[iZone])
            {
                continue;
            }

            Grid *grid = flowFieldDataOnNode[iZone]->GetGrid();
            if (!WantVisualField(grid))
            {
                continue;
            }

            if (grid->Type() == UNSTRUCTGRID)
            {
                UnstructGrid *grid_uns      = UnstructGridCast(grid);
                int           numberOfNodes = grid_uns->GetNTotalNode();
                RDouble      *qn            = (RDouble *)flowFieldDataOnNode[iZone]->GetVisualNodeVarPtr(varName);

                MpiWriteString("part", offset_write);
                offset_write += 80;

                MpiWriteInt(currentZone, offset_write);
                offset_write += sizeof(int);
                currentZone  += 1;

                MpiWriteString("coordinates", offset_write);
                offset_write += 80;

                MPI_Status status;
                for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
                {
                    float buf = static_cast<float>(qn[iNode]);

                    MPI_Offset offset_cell = offset_write + iNode * sizeof(float);
                    MPI_File_write_at(mpiFilehandle, offset_cell, &buf, sizeof(float), MPI_BYTE, &status);
                }
                offset_write += numberOfNodes * sizeof(float);
            }
            else
            {
                StructGrid *grid_str = StructGridCast(grid);
                RDouble4D  *qn       = (RDouble4D *)flowFieldDataOnNode[iZone]->GetVisualNodeVarPtr(varName);

                int ni = grid_str->GetNI();
                int nj = grid_str->GetNJ();
                int nk = grid_str->GetNK();

                MpiWriteString("part", offset_write);
                offset_write += 80;

                MpiWriteInt(currentZone, offset_write);
                offset_write += sizeof(int);
                currentZone  += 1;

                MpiWriteString("block", offset_write);
                offset_write += 80;

                MPI_Status status;
                for (int k = 1; k <= nk; ++ k)
                {
                    for (int j = 1; j <= nj; ++ j)
                    {
                        for (int i = 1; i <= ni; ++ i)
                        {
                            float vertex_buf = static_cast<float>((*qn)(i, j, k, 0));

                            MPI_File_write_at(mpiFilehandle, offset_write, &vertex_buf, sizeof(float), MPI_BYTE, &status);
                            offset_write += sizeof(float);
                        }
                    }
                }
            }
        }

        ClearFileHandle();
    }
}

void Post_WriteEnsight::UpdateCaseFile()
{
    int myid = PHMPI::GetCurrentProcessorID();

    if (myid == 0)
    {

        ofstream fout;
        fout.open(visualFileName.c_str(), ios::out | ios::trunc);

        string fname, fext;
        GetNameExt(visualFileName, fname, fext, "/");

        if (fout.is_open())
        {
            fout << "# CASE for domain visualization\n" << endl;
            fout << "FORMAT" << endl;
            fout << "type:			ensight gold\n" << endl;
            fout << "GEOMETRY" << std::endl;
            fout << "model:			" + ChangeExtensionOfFileName(fext, "geo") << endl
                << endl;

            fout << "VARIABLE" << endl;

            set<int> visualVariables = flowFieldDataOnNode[boundaryGrid[0]->GetZoneID()]->GetVisualVariables();
            for (set<int>::iterator varIter = visualVariables.begin(); varIter != visualVariables.end(); ++ varIter)
            {
                int     variableType = *varIter;
                string       varName = flowFieldDataOnNode[0]->GetVariableName(variableType);
                string fieldFileName = ChangeExtensionOfFileName(fext, varName);
                if (writeFieldEachIter)
                {
                    fieldFileName = fieldFileName + "******";
                }

                if (!dataOnElement)
                {
                    fout << "scalar per node:		1       " 
                         << varName
                         << "      "
                         << fieldFileName
                         << std::endl;
                }
                else
                {
                    fout << "scalar per element:		1       " 
                         << varName
                         << "      "
                         << fieldFileName
                         << std::endl;
                }
            }

            if (!writeFieldEachIter)
            {
                fout.close();
                return;
            }

            fout << "\nTIME" << endl;
            fout << "time set:       1" << endl;
            fout << "number of steps:        " << currentIter + 1 << endl;
            fout << "filename start number:      0" << endl;
            fout << "filename increment:         1" << endl;
            fout << "time values:" << endl;
            for (int iter = 0; iter <= currentIter; ++ iter)
            {
                fout << iter << endl;
            }
            fout.close();
        }
    }
}

void Post_WriteEnsight::ClearFieldData()
{
    for (int iData = 0; iData < flowFieldDataOnNode.size(); ++ iData)
    {
        if (flowFieldDataOnNode[iData])
        {
            FreePointer(flowFieldDataOnNode[iData]);
        }
    }

    for (unsigned int iZone = 0; iZone < VTKDataList.size(); ++ iZone)
    {
        delete VTKDataList[iZone];
    }
    VTKDataList.resize(0);
}


void Post_WriteEnsight::GetFileHandle(const char *fileName)
{
    if (fileMethod > STDIO_SERIAL)
    {

        int ret = MPI_File_open(MPI_COMM_WORLD, fileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiFilehandle);
        if (ret == MPI_SUCCESS)
        {
            ret = MPI_File_get_position(mpiFilehandle, &fileOffset);
        }
    }
    else
    {
        serialFilehandle= fopen(fileName, "wb");
    }
}

void Post_WriteEnsight::ClearFileHandle()
{
    if (fileMethod > STDIO_SERIAL)
    {
        MPI_File_close(&mpiFilehandle);
    }
    else
    {
        fclose(serialFilehandle);
    }
}

void Post_WriteEnsight::WriteString(const char *str)
{
    char buf[82];
    strncpy(buf, str, 80);
    for (int i = static_cast<int>(strlen(buf)); i < 81; ++i)
        buf[i] = '\0';

    if (fileMethod > STDIO_SERIAL)
    {

        MPI_Status status;
        MPI_File_write_at_all(mpiFilehandle, fileOffset, buf, 80, MPI_BYTE, &status);

        fileOffset += 80;

    }
    else
    {
        if (fileType)
        {
            fprintf(serialFilehandle, "%s\n", str);
        }
        else
        {
            fwrite(buf, sizeof(char), 80, serialFilehandle);
        }
    }
}

void Post_WriteEnsight::MpiWriteString(const char *str, MPI_Offset offset)
{
    char buf[82];
    strncpy(buf, str, 80);
    for (int i = static_cast<int>(strlen(buf)); i < 81; ++i)
        buf[i] = '\0';

    MPI_Status status;
    MPI_File_write_at(mpiFilehandle, offset, buf, 80, MPI_BYTE, &status);
}

void Post_WriteEnsight::MpiWriteInt(int num, MPI_Offset offset)
{
    unsigned char  _copybuf[1024];
    unsigned char *copybuf = _copybuf;
    memcpy(copybuf, &num, sizeof(int));

    MPI_Status status;
    MPI_File_write_at(mpiFilehandle, offset, copybuf, sizeof(int), MPI_BYTE, &status);
}

void Post_WriteEnsight::WriteInt(int num)
{
    unsigned char  _copybuf[1024];
    unsigned char *copybuf = _copybuf;
    memcpy(copybuf, &num, sizeof(int));

    if (fileMethod > STDIO_SERIAL)
    {
        MPI_Status status;
        MPI_File_write_at_all(mpiFilehandle, fileOffset, copybuf, sizeof(int), MPI_BYTE, &status);

        fileOffset += sizeof(int);
    }
    else
    {
        if (fileType)
        {
            fprintf(serialFilehandle, "%10d\n", num);
        }
        else
        {
            fwrite(copybuf, sizeof(int), 1, serialFilehandle);
        }
    }
}










}