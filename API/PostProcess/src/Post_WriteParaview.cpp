#include "Post_WriteParaview.h"
#include "TK_Parse.h"
#include "IO_FileName.h"
#include "PHIO.h"
#include "TK_Log.h"
#include "Glb_Dimension.h"

namespace PHSPACE
{
Post_WriteParaview::Post_WriteParaview()
{
    PHMakeDir("results/FlowField");
    this->pvdDataList.resize(0);
}

Post_WriteParaview::~Post_WriteParaview()
{

}

void Post_WriteParaview::Initialize()
{
    string visualfile = "tecflow.dat";
    GlobalDataBase::GetData("visualfile", &visualfile, PHSTRING, 1);

    string fname, fext;
    GetNameExt(visualfile, fname, fext, "/");

    visualFileName = fname + "/FlowField/" + ".vtu";
    pvdFileName    = ChangeExtensionOfFileName(visualfile, "pvd");
    VTKFileName    = visualfile;
    VTKvisual      = false;

    if (flowType == AverageFlow)
    {
        visualFileName = AddSymbolToFileName(visualFileName, "Average");
        pvdFileName    = AddSymbolToFileName(pvdFileName   , "_Average");
    }

    if (flowType == AverageReynoldsStress)
    {
        visualFileName = AddSymbolToFileName(visualFileName, "ReynoldsStress");
        pvdFileName    = AddSymbolToFileName(pvdFileName   , "_ReynoldsStress");
    }

    string sentinelfilename = "results/sentinel";
    fstream sentinelfile;
    sentinelfile.open(sentinelfilename.c_str(), ios::in);

    if (!sentinelfile && !visualfile.empty() && (flowType == 0))
    {
        VTKFileName = visualfile + ".bak";
        VTKvisual   = true;
    }

    int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
    this->iter = outnstep;

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady)
    {
        visualFileName = AddSymbolToFileName(visualFileName, outnstep, '_');
        pvdFileName    = AddSymbolToFileName(pvdFileName   , '_'     , outnstep);
    }
}

void Post_WriteParaview::StorePVDData()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        DataContainer *cdata = new DataContainer();
        cdata->MoveToBegin();

        int send_proc = GetZoneProcessorID(iZone);
        int recv_proc = GetServerProcessorID();

        if (myid == send_proc)
        {
            StorePVDData(iZone, cdata);

            int TecioMission = Writecomplete;
            PHWrite(cdata, TecioMission);
        }

        PH_Trade(cdata, send_proc, recv_proc, iZone);

        if (myid == recv_proc)
        {
            pvdDataList.push_back(cdata);
        }
        else
        {
            delete cdata;    cdata = nullptr;
        }
    }
}

void Post_WriteParaview::StorePVDData(int zoneIndex, DataContainer *cdata)
{
    PHWrite(cdata, zoneIndex);

    Grid *grid = flowFieldDataOnNode[zoneIndex]->GetGrid();
    int gridType = grid->Type();
    PHWrite(cdata, gridType);

    for (int iBoundary = 0; iBoundary < boundaryGrid.size(); ++ iBoundary)
    {
        if (boundaryGrid[iBoundary]->GetZoneID() != zoneIndex)
        {
            continue;
        }

        int TecioMission = WriteBoundary;
        PHWrite(cdata, TecioMission);

        string bcName = boundaryName[iBoundary];
        cdata->WriteString(bcName);
    }

    if (!WantVisualField(grid))
    {
        return;
    }

    int TecioMission = WriteBlock;
    PHWrite(cdata, TecioMission);

    string blockName = "BLK";
    cdata->WriteString(blockName);
}

void Post_WriteParaview::WriteFile()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    if (VTKvisual == true && myid == server)
    {
        bool CharacteristicBoundary = false;
        DumpToVTK DumpToVTK(VTKDataList, VTKFileName, CharacteristicBoundary, this->flowType);
        DumpToVTK.Run();

        WriteSentinelFile();
    }

    WriteFieldData();

    StorePVDData();

    WritePVDFile();
}

void Post_WriteParaview::WritePVDFile()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid != server)
    {
        return;
    }

    fstream pvdFile;
    OpenFile(pvdFile, pvdFileName, ios_base::out | ios_base::trunc);

    pvdFile.precision(16);
    pvdFile << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
    pvdFile << "    <Collection>" << endl;

    int bcCount = 0;

    uint_t nZones = pvdDataList.size();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        DataContainer *cdata = pvdDataList[iZone];
        cdata->MoveToBegin();

        int gridID, gridType, dumpMission;

        PHRead(cdata, gridID);
        PHRead(cdata, gridType);
        PHRead(cdata, dumpMission);

        int localPart = 0;
        while (true)
        {
            if (dumpMission == Writecomplete)
            {
                break;
            }

            string bcName;
            cdata->ReadString(bcName);

            ostringstream blockName;
            blockName << "Grid" << gridID << "_Part" << localPart<< "_" << bcName;
            string blkName = blockName.str();
            localPart ++;

            string fieldFileName = "FlowField/.vtu";
            fieldFileName = AddSymbolToFileName(fieldFileName, blkName);

            if (gridType == UNSTRUCTGRID)
            {
                fieldFileName = ChangeExtensionOfFileName(fieldFileName, "vtu");
            }
            else
            {
                fieldFileName = ChangeExtensionOfFileName(fieldFileName, "vts");
            }
            pvdFile << "		<DataSet " << "timestep=\"" << iter << "\" part=\"" << bcCount << "\" file=\"" << fieldFileName << "\" ";
            pvdFile << "name = \"" << blkName << "\"/>" << endl;

            PHRead(cdata, dumpMission);
            bcCount ++;
        }
    }

    pvdFile << "    </Collection>" << endl;
    pvdFile << "</VTKFile>\n";

    CloseFile(pvdFile);
}

void Post_WriteParaview::WriteFieldData()
{
    for (int iZone = 0; iZone < flowFieldDataOnNode.size(); ++ iZone)
    {
        if (!flowFieldDataOnNode[iZone])
        {
            continue;
        }

        int localPart = 0;
        //! Write boundary grid.
        for (int iBoundary = 0; iBoundary < boundaryGrid.size(); ++ iBoundary)
        {
            Grid *grid = boundaryGrid[iBoundary];
            int oriGridIndex = grid->GetZoneID();
            if (oriGridIndex != iZone)
            {
                continue;
            }

            string bcName = boundaryName[iBoundary];
            string fieldFileName = AddSymbolToFileName(visualFileName, "Grid", oriGridIndex);
            fieldFileName = AddSymbolToFileName(fieldFileName, "_Part", localPart);
            fieldFileName = AddSymbolToFileName(fieldFileName, "_", bcName);
            localPart ++;

            fstream flowFile;

            if (grid->Type() == UNSTRUCTGRID)
            {
                fieldFileName = ChangeExtensionOfFileName(fieldFileName, "vtu");
                OpenFile(flowFile, fieldFileName, ios_base::out | ios_base::trunc);

                WriteUnsBoundaryFieldData(flowFile, iBoundary);
            }
            else
            {
                fieldFileName = ChangeExtensionOfFileName(fieldFileName, "vts");
                OpenFile(flowFile, fieldFileName, ios_base::out | ios_base::trunc);

                WriteStrBoundaryFieldData(flowFile, iBoundary);
            }

            CloseFile(flowFile);
        }


        //! Write block grid.
        Grid *grid = flowFieldDataOnNode[iZone]->GetGrid();
        if (!WantVisualField(grid))
        {
            continue;
        }

        int oriGridIndex = grid->GetZoneID();
        string fieldFileName = AddSymbolToFileName(visualFileName, "Grid", oriGridIndex);
        fieldFileName = AddSymbolToFileName(fieldFileName, "_Part", localPart);
        fieldFileName = AddSymbolToFileName(fieldFileName, "_", "BLK");
        localPart ++;

        fstream flowFile;

        if (grid->Type() == UNSTRUCTGRID)
        {
            fieldFileName = ChangeExtensionOfFileName(fieldFileName, "vtu");
            OpenFile(flowFile, fieldFileName, ios_base::out | ios_base::trunc);

            WriteUnsBlockFieldData(flowFile, iZone);
        }
        else
        {
            fieldFileName = ChangeExtensionOfFileName(fieldFileName, "vts");
            OpenFile(flowFile, fieldFileName, ios_base::out | ios_base::trunc);

            WriteStrBlockFieldData(flowFile, iZone);
        }

        CloseFile(flowFile);
    }
}

void Post_WriteParaview::WriteUnsBoundaryFieldData(fstream &flowFile, int iBoundary)
{
    flowFile.precision(16);
    flowFile << "<?xml version=\"1.0\" ?>" << endl;
    flowFile << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
    flowFile << "<!-- ITER " << iter << " -->" << endl;
    flowFile << "	<UnstructuredGrid>" << endl;

    Grid *grid          = boundaryGrid[iBoundary];
    int   numberOfNodes = grid->GetNTotalNode();
    int   numberOfCells = grid->GetNTotalCell();

    flowFile << "		<Piece NumberOfPoints=\"" << numberOfNodes << "\" NumberOfCells=\"" << numberOfCells << "\">" << endl;
    flowFile << "			<PointData>" << endl;

    int oriGridIndex = grid->GetZoneID();
    set<int> visualVariables = flowFieldDataOnNode[oriGridIndex]->GetVisualVariables();


    //! Write velocity
    bool visualU = true;
    bool visualV = true;
    bool visualW = true;

    RDouble *qU = nullptr;
    if (flowFieldDataOnNode[oriGridIndex]->IsNeedVisualization(VISUAL_U))
    {
        string varName = flowFieldDataOnNode[oriGridIndex]->GetVariableName(VISUAL_U);
        qU = (RDouble *)flowFieldDataOnNode[oriGridIndex]->GetVisualNodeVarPtr(varName);
    }
    else
    {
        qU = new RDouble[numberOfNodes];
        SetField(qU, 0.0, numberOfNodes);

        visualU = false;
    }

    RDouble *qV = nullptr;
    if (flowFieldDataOnNode[oriGridIndex]->IsNeedVisualization(VISUAL_V))
    {
        string varName = flowFieldDataOnNode[oriGridIndex]->GetVariableName(VISUAL_V);
        qV = (RDouble *)flowFieldDataOnNode[oriGridIndex]->GetVisualNodeVarPtr(varName);
    }
    else
    {
        qV = new RDouble[numberOfNodes];
        SetField(qV, 0.0, numberOfNodes);

        visualV = false;
    }

    RDouble *qW = nullptr;
    if (flowFieldDataOnNode[oriGridIndex]->IsNeedVisualization(VISUAL_W))
    {
        string varName = flowFieldDataOnNode[oriGridIndex]->GetVariableName(VISUAL_W);
        qW = (RDouble *)flowFieldDataOnNode[oriGridIndex]->GetVisualNodeVarPtr(varName);
    }
    else
    {
        qW = new RDouble[numberOfNodes];
        SetField(qW, 0.0, numberOfNodes);

        visualW = false;
    }

    flowFile << "				<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">" << endl;
    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        int originaNodeIndex = nodeMap[iBoundary][0][iNode];
        flowFile << qU[originaNodeIndex] << " " << qV[originaNodeIndex] << " " << qW[originaNodeIndex] << " ";
    }
    flowFile << endl;
    flowFile << "				</DataArray>" << endl;

    if (!visualU)
    {
        delete [] qU;
        qU = nullptr;
    }

    if (!visualV)
    {
        delete [] qV;
        qV = nullptr;
    }

    if (!visualW)
    {
        delete [] qW;
        qW = nullptr;
    }


    //! Write other data.
    for (set<int>::iterator varIter = visualVariables.begin(); varIter != visualVariables.end(); ++ varIter)
    {
        int variableType = *varIter;
        if (variableType == VISUAL_U || variableType == VISUAL_V || variableType == VISUAL_W)
        {
            continue;
        }

        string  varName = flowFieldDataOnNode[oriGridIndex]->GetVariableName(variableType);
        RDouble *qn = (RDouble *)flowFieldDataOnNode[oriGridIndex]->GetVisualNodeVarPtr(varName);

        flowFile << "				<DataArray type=\"Float64\" Name=\"" << varName << "\" format=\"ascii\">" << endl;
        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            int originaNodeIndex = nodeMap[iBoundary][0][iNode];
            flowFile << qn[originaNodeIndex] << " ";
        }
        flowFile << endl;
        flowFile << "				</DataArray>" << endl;
    }
    flowFile << "			</PointData>" << endl;


    //! Points coordinates
    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    flowFile << "			<Points>" << endl;
    flowFile << "				<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        flowFile << x[iNode] << " " << y[iNode] << " " << z[iNode] << " " << endl;
    }
    flowFile << endl;
    flowFile << "				</DataArray>" << endl;
    flowFile << "			</Points>" << endl;


    //! Write out Cell data: connectivity, offsets, element types.
    flowFile << "			<Cells>" << endl;

    //! Write connectivity array
    int *cell2node            = UnstructGridCast(grid)->GetCell2Node();
    int *nodeNumberOfEachCell = UnstructGridCast(grid)->GetNodeNumberOfEachCell();

    flowFile << "				<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
    int count = 0;
    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        for (int iNode = 0; iNode < nodeNumberOfEachCell[iCell]; ++ iNode)
        {
            flowFile << cell2node[count] << "  ";
            count ++;
        }
        flowFile << endl;
    }
    flowFile << "				</DataArray>" << endl;

    //! Write cell-node offsets
    int nvPerCell = 4;
    flowFile << "				<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        flowFile << (iCell + 1) * nvPerCell << " ";
    }
    flowFile << endl;
    flowFile << "				</DataArray>" << endl;

    //! Write VTK element type
    int eType;
    if (PHSPACE::GetDim() == TWO_D)
    {
        eType = VTK_LINE;
    }
    else
    {
        eType = VTK_QUAD;
    }
    flowFile << "				<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        flowFile << eType << " ";
    }
    flowFile << endl;
    flowFile << "				</DataArray>" << endl;

    flowFile << "			</Cells>" << endl;
    flowFile << "		</Piece>" << endl;
    flowFile << "	</UnstructuredGrid>" << endl;
    flowFile << "</VTKFile>" << endl;
}

void Post_WriteParaview::WriteStrBoundaryFieldData(fstream &flowFile, int iBoundary)
{
    StructGrid *grid_str = StructGridCast(boundaryGrid[iBoundary]);

    flowFile.precision(16);
    flowFile << "<?xml version=\"1.0\" ?>" << endl;
    flowFile << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
    flowFile << "<!-- ITER " << iter << " -->" << endl;

    int ni         = grid_str->GetNI();
    int nj         = grid_str->GetNJ();
    int nk         = grid_str->GetNK();
    int nTotalNode = grid_str->GetNTotalNode();

    flowFile << "    <StructuredGrid WholeExtent=\""
             << 1 << " " << ni << " "
             << 1 << " " << nj << " "
             << 1 << " " << nk << " "
             << "\">" << endl;

    flowFile << "        <Piece Extent=\""
             << 1 << " " << ni << " "
             << 1 << " " << nj << " "
             << 1 << " " << nk << " "
             << "\">" << endl;

    flowFile << "			<PointData>" << endl;

    int      oriGridIndex    = grid_str->GetZoneID();
    set<int> visualVariables = flowFieldDataOnNode[oriGridIndex]->GetVisualVariables();


    //! Write velocity
    bool visualU = true;
    bool visualV = true;
    bool visualW = true;

    Range I(-1, ni + 1);
    Range J(-1, nj + 1);
    Range K(-1, nk + 1);
    if (nk == 1) K.setRange(1, 1);
    Range M(0, 0);

    RDouble4D *qU = nullptr;
    if (flowFieldDataOnNode[oriGridIndex]->IsNeedVisualization(VISUAL_U))
    {
        string varName = flowFieldDataOnNode[oriGridIndex]->GetVariableName(VISUAL_U);
        qU = (RDouble4D*)flowFieldDataOnNode[oriGridIndex]->GetVisualNodeVarPtr(varName);
    }
    else
    {
        qU = new RDouble4D(I, J, K, M, fortranArray);
        *qU = 0.0;

        visualU = false;
    }

    RDouble4D *qV = nullptr;
    if (flowFieldDataOnNode[oriGridIndex]->IsNeedVisualization(VISUAL_V))
    {
        string varName = flowFieldDataOnNode[oriGridIndex]->GetVariableName(VISUAL_V);
        qV = (RDouble4D*)flowFieldDataOnNode[oriGridIndex]->GetVisualNodeVarPtr(varName);
    }
    else
    {
        qV = new RDouble4D(I, J, K, M, fortranArray);
        *qV = 0.0;

        visualV = false;
    }

    RDouble4D *qW = nullptr;
    if (flowFieldDataOnNode[oriGridIndex]->IsNeedVisualization(VISUAL_W))
    {
        string varName = flowFieldDataOnNode[oriGridIndex]->GetVariableName(VISUAL_W);
        qW = (RDouble4D*)flowFieldDataOnNode[oriGridIndex]->GetVisualNodeVarPtr(varName);
    }
    else
    {
        qW = new RDouble4D(I, J, K, M, fortranArray);
        *qW = 0.0;

        visualW = false;
    }

    flowFile << "				<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">" << endl;
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        int originaI = nodeMap[iBoundary][0][iNode];
        int originaJ = nodeMap[iBoundary][1][iNode];
        int originaK = nodeMap[iBoundary][2][iNode];

        RDouble dataU = (*qU)(originaI, originaJ, originaK, 0);
        RDouble dataV = (*qV)(originaI, originaJ, originaK, 0);
        RDouble dataW = (*qW)(originaI, originaJ, originaK, 0);

        flowFile << dataU << " " << dataV << " " << dataW << " ";
    }
    flowFile << endl;
    flowFile << "				</DataArray>" << endl;

    if (!visualU)
    {
        delete qU;
        qU = nullptr;
    }

    if (!visualV)
    {
        delete qV;
        qV = nullptr;
    }

    if (!visualW)
    {
        delete qW;
        qW = nullptr;
    }


    //! Write other data.
    for (set<int>::iterator varIter = visualVariables.begin(); varIter != visualVariables.end(); ++ varIter)
    {
        int variableType = *varIter;
        if (variableType == VISUAL_U || variableType == VISUAL_V || variableType == VISUAL_W)
        {
            continue;
        }

        string     varName = flowFieldDataOnNode[oriGridIndex]->GetVariableName(variableType);
        RDouble4D *qn      = (RDouble4D*)flowFieldDataOnNode[oriGridIndex]->GetVisualNodeVarPtr(varName);

        flowFile << "				<DataArray type=\"Float64\" Name=\"" << varName << "\" format=\"ascii\">" << endl;
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            int originaI = nodeMap[iBoundary][0][iNode];
            int originaJ = nodeMap[iBoundary][1][iNode];
            int originaK = nodeMap[iBoundary][2][iNode];

            RDouble writeData = (*qn)(originaI, originaJ, originaK, 0);
            flowFile << writeData << " ";
        }

        flowFile << endl;
        flowFile << "				</DataArray>" << endl;
    }
    flowFile << "			</PointData>" << endl;


    //! Points coordinates
    RDouble *x = grid_str->GetX();
    RDouble *y = grid_str->GetY();
    RDouble *z = grid_str->GetZ();

    flowFile << "			<Points>" << endl;
    flowFile << "				<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        flowFile << x[iNode] << " " << y[iNode] << " " << z[iNode] << " ";
    }

    flowFile << endl;
    flowFile << "				</DataArray>" << endl;
    flowFile << "			</Points>" << endl;
    flowFile << "        </Piece>" << std::endl;
    flowFile << "    </StructuredGrid>" << endl;
    flowFile << "</VTKFile>" << endl;
}

void Post_WriteParaview::WriteUnsBlockFieldData(fstream &flowFile, int iZone)
{
    UnstructGrid *grid_uns = UnstructGridCast(flowFieldDataOnNode[iZone]->GetGrid());

    flowFile.precision(16);
    flowFile << "<?xml version=\"1.0\" ?>" << endl;
    flowFile << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
    flowFile << "	<UnstructuredGrid>" << endl;

    int numberOfNodes = grid_uns->GetNTotalNode();
    int numberOfCells = grid_uns->GetNTotalCell();

    flowFile << "		<Piece NumberOfPoints=\"" << numberOfNodes << "\" NumberOfCells=\"" << numberOfCells << "\">" << endl;
    flowFile << "			<PointData>" << endl;

    set<int> visualVariables = flowFieldDataOnNode[iZone]->GetVisualVariables();


    //! Write velocity
    bool visualU = true;
    bool visualV = true;
    bool visualW = true;

    RDouble *qU = nullptr;
    if (flowFieldDataOnNode[iZone]->IsNeedVisualization(VISUAL_U))
    {
        string varName = flowFieldDataOnNode[iZone]->GetVariableName(VISUAL_U);
        qU = (RDouble *)flowFieldDataOnNode[iZone]->GetVisualNodeVarPtr(varName);
    }
    else
    {
        qU = new RDouble[numberOfNodes];
        SetField(qU, 0.0, numberOfNodes);

        visualU = false;
    }

    RDouble *qV = nullptr;
    if (flowFieldDataOnNode[iZone]->IsNeedVisualization(VISUAL_V))
    {
        string varName = flowFieldDataOnNode[iZone]->GetVariableName(VISUAL_V);
        qV = (RDouble *)flowFieldDataOnNode[iZone]->GetVisualNodeVarPtr(varName);
    }
    else
    {
        qV = new RDouble[numberOfNodes];
        SetField(qV, 0.0, numberOfNodes);

        visualV = false;
    }

    RDouble *qW = nullptr;
    if (flowFieldDataOnNode[iZone]->IsNeedVisualization(VISUAL_W))
    {
        string varName = flowFieldDataOnNode[iZone]->GetVariableName(VISUAL_W);
        qW = (RDouble *)flowFieldDataOnNode[iZone]->GetVisualNodeVarPtr(varName);
    }
    else
    {
        qW = new RDouble[numberOfNodes];
        SetField(qW, 0.0, numberOfNodes);

        visualW = false;
    }

    flowFile << "				<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">" << endl;
    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        flowFile << qU[iNode] << " " << qV[iNode] << " " << qW[iNode] << " ";
    }
    flowFile << endl;
    flowFile << "				</DataArray>" << endl;

    if (!visualU)
    {
        delete [] qU;
        qU = nullptr;
    }

    if (!visualV)
    {
        delete [] qV;
        qV = nullptr;
    }

    if (!visualW)
    {
        delete [] qW;
        qW = nullptr;
    }


    //! Write other data.
    for (set<int>::iterator varIter = visualVariables.begin(); varIter != visualVariables.end(); ++ varIter)
    {
        int variableType = *varIter;
        if (variableType == VISUAL_U || variableType == VISUAL_V || variableType == VISUAL_W)
        {
            continue;
        }

        string varName = flowFieldDataOnNode[iZone]->GetVariableName(variableType);
        RDouble *qn = (RDouble *)flowFieldDataOnNode[iZone]->GetVisualNodeVarPtr(varName);

        flowFile << "				<DataArray type=\"Float64\" Name=\"" << varName << "\" format=\"ascii\">" << endl;
        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            flowFile << qn[iNode] << " ";
        }
        flowFile << endl;
        flowFile << "				</DataArray>" << endl;
    }
    flowFile << "			</PointData>" << endl;


    //! Points coordinates
    RDouble *x = grid_uns->GetX();
    RDouble *y = grid_uns->GetY();
    RDouble *z = grid_uns->GetZ();

    flowFile << "			<Points>" << endl;
    flowFile << "				<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        flowFile << x[iNode] << " " << y[iNode] << " " << z[iNode] << " " << endl;
    }
    flowFile << endl;
    flowFile << "				</DataArray>" << endl;
    flowFile << "			</Points>" << endl;


    //! Write out Cell data: connectivity, offsets, element types.
    flowFile << "			<Cells>" << endl;

    //! Write connectivity array
    int *cell2node            = grid_uns->GetCell2Node();
    int *nodeNumberOfEachCell = grid_uns->GetNodeNumberOfEachCell();

    flowFile << "				<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
    int count = 0;
    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        for (int iNode = 0; iNode < nodeNumberOfEachCell[iCell]; ++ iNode)
        {
            flowFile << cell2node[count] << "  ";
            count ++;
        }
        flowFile << endl;
    }
    flowFile << "				</DataArray>" << endl;

    //! Write cell-node offsets
    vector<int> offset(numberOfCells);
    int tmp = 0;
    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        tmp += nodeNumberOfEachCell[iCell];
        offset[iCell] = tmp;
    }

    flowFile << "				<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\""
             << " RangeMin=\"" << offset[0] << "\" RangeMax=\"" << offset[numberOfCells - 1] << "\">" << endl;
    flowFile << "					";

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        flowFile << offset[iCell] << " ";
    }
    flowFile << endl;
    flowFile << "				</DataArray>" << endl;

    //! Write VTK element type
    if (PHSPACE::GetDim() == TWO_D)
    {
        flowFile << "				<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" RangeMin=\"5\""
                 << " RangeMax=\"9\">" << endl;
    }
    else if (PHSPACE::GetDim() == THREE_D)
    {
        flowFile << "				<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\" RangeMin=\"12\""
                 << " RangeMax=\"13\">" << endl;
    }

    int eType = 0;
    flowFile << "					";
    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        if (PHSPACE::GetDim() == TWO_D)
        {
            if (nodeNumberOfEachCell[iCell] == 4)
            {
                eType = VTK_QUAD;
                flowFile << eType << " ";
            }
            else if (nodeNumberOfEachCell[iCell] == 3)
            {
                eType = VTK_TRIANGLE;
                flowFile << eType << " ";
            }
        }
        else if (PHSPACE::GetDim() == THREE_D)
        {
            if (nodeNumberOfEachCell[iCell] == 8)
            {
                eType = VTK_HEXAHEDRON;
                flowFile << eType << " ";
            }
            else if (nodeNumberOfEachCell[iCell] == 6)
            {
                eType = VTK_WEDGE;
                flowFile << eType << " ";
            }
        }
    }
    flowFile << endl;
    flowFile << "				</DataArray>" << endl;

    flowFile << "			</Cells>" << endl;
    flowFile << "		</Piece>" << endl;
    flowFile << "	</UnstructuredGrid>" << endl;
    flowFile << "</VTKFile>" << endl;
}

void Post_WriteParaview::WriteStrBlockFieldData(fstream &flowFile, int iZone)
{
    StructGrid *grid_str = StructGridCast(flowFieldDataOnNode[iZone]->GetGrid());

    flowFile.precision(16);
    flowFile << "<?xml version=\"1.0\" ?>" << endl;
    flowFile << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
    flowFile << "<!-- ITER " << iter << " -->" << endl;

    int ni         = grid_str->GetNI();
    int nj         = grid_str->GetNJ();
    int nk         = grid_str->GetNK();

    flowFile << "    <StructuredGrid WholeExtent=\""
        << 1 << " " << ni << " "
        << 1 << " " << nj << " "
        << 1 << " " << nk << " "
        << "\">" << endl;

    flowFile << "        <Piece Extent=\""
        << 1 << " " << ni << " "
        << 1 << " " << nj << " "
        << 1 << " " << nk << " "
        << "\">" << endl;

    flowFile << "			<PointData>" << endl;
    set<int> visualVariables = flowFieldDataOnNode[iZone]->GetVisualVariables();


    //! Write velocity
    bool visualU = true;
    bool visualV = true;
    bool visualW = true;

    Range I(-1, ni + 1);
    Range J(-1, nj + 1);
    Range K(-1, nk + 1);
    if (nk == 1) K.setRange(1, 1);
    Range M(0, 0);

    RDouble4D *qU = nullptr;
    if (flowFieldDataOnNode[iZone]->IsNeedVisualization(VISUAL_U))
    {
        string varName = flowFieldDataOnNode[iZone]->GetVariableName(VISUAL_U);
        qU = (RDouble4D*)flowFieldDataOnNode[iZone]->GetVisualNodeVarPtr(varName);
    }
    else
    {
        qU = new RDouble4D(I, J, K, M, fortranArray);
        *qU = 0.0;

        visualU = false;
    }

    RDouble4D *qV = nullptr;
    if (flowFieldDataOnNode[iZone]->IsNeedVisualization(VISUAL_V))
    {
        string varName = flowFieldDataOnNode[iZone]->GetVariableName(VISUAL_V);
        qV = (RDouble4D*)flowFieldDataOnNode[iZone]->GetVisualNodeVarPtr(varName);
    }
    else
    {
        qV = new RDouble4D(I, J, K, M, fortranArray);
        *qV = 0.0;

        visualV = false;
    }

    RDouble4D *qW = nullptr;
    if (flowFieldDataOnNode[iZone]->IsNeedVisualization(VISUAL_W))
    {
        string varName = flowFieldDataOnNode[iZone]->GetVariableName(VISUAL_W);
        qW = (RDouble4D*)flowFieldDataOnNode[iZone]->GetVisualNodeVarPtr(varName);
    }
    else
    {
        qW = new RDouble4D(I, J, K, M, fortranArray);
        *qW = 0.0;

        visualW = false;
    }

    flowFile << "				<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"Velocity\" format=\"ascii\">" << endl;
    for (int k = 1; k <= nk; ++ k)
    {
        for (int j = 1; j <= nj; ++ j)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                RDouble dataU = (*qU)(i, j, k, 0);
                RDouble dataV = (*qV)(i, j, k, 0);
                RDouble dataW = (*qW)(i, j, k, 0);

                flowFile << dataU << " " << dataV << " " << dataW << " ";
            }
        }
    }
    flowFile << endl;
    flowFile << "				</DataArray>" << endl;

    if (!visualU)
    {
        delete qU;
        qU = nullptr;
    }

    if (!visualV)
    {
        delete qV;
        qV = nullptr;
    }

    if (!visualW)
    {
        delete qW;
        qW = nullptr;
    }


    //! Write other data.
    for (set<int>::iterator varIter = visualVariables.begin(); varIter != visualVariables.end(); ++ varIter)
    {
        int variableType = *varIter;
        if (variableType == VISUAL_U || variableType == VISUAL_V || variableType == VISUAL_W)
        {
            continue;
        }

        string     varName = flowFieldDataOnNode[iZone]->GetVariableName(variableType);
        RDouble4D *qn = (RDouble4D*)flowFieldDataOnNode[iZone]->GetVisualNodeVarPtr(varName);

        flowFile << "				<DataArray type=\"Float64\" Name=\"" << varName << "\" format=\"ascii\">" << endl;
        for (int k = 1; k <= nk; ++ k)
        {
            for (int j = 1; j <= nj; ++ j)
            {
                for (int i = 1; i <= ni; ++ i)
                {
                    RDouble writeData = (*qn)(i, j, k, 0);
                    flowFile << writeData << " ";
                }
            }
        }
        flowFile << endl;
        flowFile << "				</DataArray>" << endl;
    }
    flowFile << "			</PointData>" << endl;


    //! Points coordinates
    RDouble3D &x = *grid_str->GetStructX();
    RDouble3D &y = *grid_str->GetStructY();
    RDouble3D &z = *grid_str->GetStructZ();

    flowFile << "			<Points>" << endl;
    flowFile << "				<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

    for (int k = 1; k <= nk; ++ k)
    {
        for (int j = 1; j <= nj; ++ j)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                flowFile << x(i, j, k) << " " << y(i, j, k) << " " << z(i, j, k) << " ";
            }
        }
    }

    flowFile << endl;
    flowFile << "				</DataArray>" << endl;
    flowFile << "			</Points>" << endl;
    flowFile << "        </Piece>" << std::endl;
    flowFile << "    </StructuredGrid>" << endl;
    flowFile << "</VTKFile>" << endl;
}

void Post_WriteParaview::ClearFieldData()
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

    for (unsigned int iZone = 0; iZone < pvdDataList.size(); ++ iZone)
    {
        delete pvdDataList[iZone];
    }
    pvdDataList.resize(0);
}

}