#include "Post_Probes.h"
#include "TK_Exit.h"
#include "Pointer.h"
#include "IO_FileName.h"
#include "PHIO.h"
#include "Constants.h"
#include "Glb_Dimension.h"
#include "TK_Log.h"
#include "TK_Time.h"
#include "Gas.h"

using namespace std;
namespace PHSPACE
{
LIB_EXPORT Post_Probes::Post_Probes(int nProbeVariables_in, int *probeVariables_in, Grid *grid_in)
{
    probesVariablesMap.insert(pair<int, string>(PROBE_DENSITY    , "density"));
    probesVariablesMap.insert(pair<int, string>(PROBE_U          , "u"));
    probesVariablesMap.insert(pair<int, string>(PROBE_V          , "v"));
    probesVariablesMap.insert(pair<int, string>(PROBE_W          , "w"));
    probesVariablesMap.insert(pair<int, string>(PROBE_PRESSURE   , "pressure"));
    probesVariablesMap.insert(pair<int, string>(PROBE_TEMPERATURE, "temperature"));
    probesVariablesMap.insert(pair<int, string>(PROBE_MACH, "mach"));
    probesVariablesMap.insert(pair<int, string>(PROBE_DIMENSIONAL_DENSITY,     "density_Dim"));
    probesVariablesMap.insert(pair<int, string>(PROBE_DIMENSIONAL_U      ,     "u_Dim"));
    probesVariablesMap.insert(pair<int, string>(PROBE_DIMENSIONAL_V      ,     "v_Dim"));
    probesVariablesMap.insert(pair<int, string>(PROBE_DIMENSIONAL_W      ,     "w_Dim"));
    probesVariablesMap.insert(pair<int, string>(PROBE_VELOCITY_MAGNITUDE ,     "Vmag_Dim"));
    probesVariablesMap.insert(pair<int, string>(PROBE_DIMENSIONAL_PRESSURE,    "pressure_Dim"));
    probesVariablesMap.insert(pair<int, string>(PROBE_DIMENSIONAL_TEMPERATURE, "temperature_Dim"));
    probesVariablesMap.insert(pair<int, string>(PROBE_VISCOSITY_LAMINAR  , "viscosityLaminar"));
    probesVariablesMap.insert(pair<int, string>(PROBE_VISCOSITY_TURBULENT, "viscosityTurbulence"));
    probesVariablesMap.insert(pair<int, string>(PROBE_MODELED_TKE        , "modeledTKE"));
    probesVariablesMap.insert(pair<int, string>(PROBE_MODELED_DISSIPATION, "modeleddissipationrate"));

    nProbeVariables = nProbeVariables_in;
    int nCount = 0;
    for (int iVar = 0; iVar < nProbeVariables; ++ iVar)
    {
        int varTemp = probeVariables_in[iVar];
        if (probesVariablesMap.find(varTemp) == probesVariablesMap.end())
        {
            continue;
        }
        probeVariables.insert(varTemp);
        ++ nCount;
    }
    nProbeVariables = nCount;

    probeVariablesPtr = new Data_Field();

    this->grid = grid_in;

    if (grid_in)
    {
        gridType = grid_in->Type();
    }
    else
    {
        gridType = -1;
    }
}

LIB_EXPORT Post_Probes::~Post_Probes()
{
    FreeAllProbeVarPtr();
    delete probeVariablesPtr;
}

LIB_EXPORT void Post_Probes::GetASCIIFileHeader(vector<string> &title_tecplot)
{
    if (nProbeVariables <= 0) return;
    int dataMonitorType = GlobalDataBase::GetIntParaFromDB("dataMonitorType");

    if (dataMonitorType == PROBESMONITOR)
    {
        GetProbesASCIIFileHeader(title_tecplot);
    }
    else if (dataMonitorType == LINESMONITOR || dataMonitorType == SURFACESMONITOR)
    {
        GetLinesOrSurfacesASCIIFileHeader(title_tecplot);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("dataMonitorType", dataMonitorType);
    }
}

LIB_EXPORT void Post_Probes::GetProbesASCIIFileHeader(vector<string> &title_tecplot)
{
    int samplefileMode = GlobalDataBase::GetIntParaFromDB("samplefileMode");
    ostringstream oss;
    oss << "title = \"Probes Flow Fields Data of PHengLEI\"" << endl;
    oss << "Variables =" << endl;
    oss << "\"iter\"" << endl;

    if (samplefileMode != 0)
    {
        oss << "\"x\"" << endl;
        oss << "\"y\"" << endl;
        oss << "\"z\"" << endl;
    }

    set<int>::iterator varIter;
    for (varIter = probeVariables.begin(); varIter != probeVariables.end(); ++ varIter)
    {
        int variableType = *varIter;
        string varName = probesVariablesMap[variableType];
        string varNameTemp = "\"" + varName + "\"";
        oss << varNameTemp << endl;
    }
    //! Chemical species must be the last!
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        //! Chemical turned on.
        using namespace GAS_SPACE;
        int nm    = GlobalDataBase::GetIntParaFromDB("nm");
        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int neqn  = nl + nchem;

        string *varname = gas->GetNameOfSpecies();
        for (int m = nm; m < neqn; ++ m)
        {
            string varName = "massfraction-" + varname[m-nm];
            string species_name = "\"" + varName + "\"";
            oss << species_name << endl;
        }
        for (int m = nm; m < neqn; ++ m)
        {
            string varName = "molefraction-" + varname[m-nm];
            string species_name = "\"" + varName + "\"";
            oss << species_name << endl;
        }
    }

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady)
    {
        oss << "\"physicalTime\"" << endl;
    }
    title_tecplot.push_back(oss.str());
    return;
}

LIB_EXPORT void Post_Probes::GetLinesOrSurfacesASCIIFileHeader(vector<string> &title_tecplot)
{
    ostringstream oss;
    oss << "title = \"Probes Flow Fields Data of PHengLEI\"" << endl;
    oss << "Variables =" << endl;
    oss << "\"iter\"" << endl;
    oss << "\"x\"" << endl;
    oss << "\"y\"" << endl;
    oss << "\"z\"" << endl;
    set<int>::iterator varIter;
    for (varIter = probeVariables.begin(); varIter != probeVariables.end(); ++ varIter)
    {
        int variableType = *varIter;
        string varName = probesVariablesMap[variableType];
        string varNameTemp = "\"" + varName + "\"";
        oss << varNameTemp << endl;
    }

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady)
    {
        oss << "\"physicalTime\"" << endl;
    }
    title_tecplot.push_back(oss.str());
    return;
}

LIB_EXPORT void Post_Probes::UpdateProbesVarPtr(const string &name, RDouble *cellCenterData)
{
    if (!grid)
    {
        TK_Exit::ExceptionExit("Post_Probe::UpdateProbesVarPtr needs Geometry infomation.");
    }

    int probeVariablesInterpolationMethod = GlobalDataBase::GetIntParaFromDB("probeVariablesInterpolationMethod");

    if (probeVariablesInterpolationMethod == CELLVALUE)
    {
        UpdateProbesVarPtrM1(name, cellCenterData);
    }
    else if (probeVariablesInterpolationMethod == CELLSTOPROBE)
    {
        UpdateProbesVarPtrM2(name, cellCenterData);
    }
    else if (probeVariablesInterpolationMethod == NODESTOPROBE)
    {
        UpdateProbesVarPtrM3(name, cellCenterData);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("probeVariablesInterpolationMethod", probeVariablesInterpolationMethod);
    }
}

LIB_EXPORT void Post_Probes::UpdateProbesVarPtrM1(const string &name, RDouble *cellCenterData)
{
    UnstructGrid *unstructGrid = UnstructGridCast(grid);
    int probesLocalNumber = unstructGrid->GetZoneProbesNumber();
    if (probesLocalNumber == 0) return;

    vector<int> probesCellID = unstructGrid->GetZoneProbesCellID();

    RDouble *q_probes = new RDouble [probesLocalNumber];
    SetField(q_probes, 0.0, probesLocalNumber);

    for (int iProbe = 0; iProbe < probesLocalNumber; ++ iProbe)
    {
        int probeCell = probesCellID[iProbe];
        q_probes[iProbe] = cellCenterData[probeCell];
    }
    probeVariablesPtr->UpdateDataPtr(name, q_probes);
}

LIB_EXPORT void Post_Probes::UpdateProbesVarPtrM2(const string &name, RDouble *cellCenterData)
{
    UnstructGrid *unstructGrid = UnstructGridCast(grid);
    int nBoundFace = unstructGrid->GetNBoundFace();
    int nTotalCell = unstructGrid->GetNTotalCell();
    int **cell2face = unstructGrid->GetCell2Face();
    int *faceNumberOfEachCell = unstructGrid->GetFaceNumberOfEachCell();

    int probesLocalNumber = unstructGrid->GetZoneProbesNumber();
    if (probesLocalNumber == 0) return;

    vector<int> *cell2Cell = unstructGrid->GetCell2Cell();
    vector<int> probesCellID = unstructGrid->GetZoneProbesCellID();

    RDouble *q_probes = new RDouble [probesLocalNumber];
    SetField(q_probes, 0.0, probesLocalNumber);

    for (int iProbe = 0; iProbe < probesLocalNumber; ++ iProbe)
    {
        vector<int> interpolationCellsID;
        vector<RDouble> numerator;

        int probeCell = probesCellID[iProbe];
        NecessaryVarForInterpolation(unstructGrid, iProbe, probeCell, numerator,interpolationCellsID, true);

        int nNeighbor = static_cast<int>(cell2Cell[probeCell].size());
        for (int iNeighbor = 0; iNeighbor < nNeighbor; ++ iNeighbor)
        {
            int neighborCell = cell2Cell[probeCell][iNeighbor];
            NecessaryVarForInterpolation(unstructGrid, iProbe, neighborCell, numerator,interpolationCellsID, true);
        }

        for (int iFace = 0; iFace < faceNumberOfEachCell[probeCell]; ++ iFace)
        {
            int faceID = cell2face[probeCell][iFace];

            if (faceID < nBoundFace)
            {
                int neighborGhostCell = faceID + nTotalCell;
                NecessaryVarForInterpolation(unstructGrid, iProbe, neighborGhostCell, numerator,interpolationCellsID, true);
            }
        }

        int nCells = static_cast<int>(interpolationCellsID.size());
        RDouble denominator = SUMMARY(numerator, nCells);

        for (int iCell = 0; iCell < nCells; ++ iCell)
        {
            RDouble weightCoef = numerator[iCell] / denominator;
            int cellID = interpolationCellsID[iCell];
            q_probes[iProbe] += cellCenterData[cellID] * weightCoef;
        }
    }
    probeVariablesPtr->UpdateDataPtr(name, q_probes);
}

LIB_EXPORT void Post_Probes::UpdateProbesVarPtrM3(const string &name, RDouble *cellCenterData)
{
    UnstructGrid *unstructGrid = UnstructGridCast(grid);
    int probesLocalNumber = unstructGrid->GetZoneProbesNumber();
    if (probesLocalNumber == 0) return;

    int nTotalNode = unstructGrid->GetNTotalNode();
    int **cell2Node = unstructGrid->GetCell2NodeArray();
    int *nodeNumberOfEachCell = unstructGrid->GetNodeNumberOfEachCell();

    RDouble *qn = new RDouble [nTotalNode];
    CompNodeVar(unstructGrid, qn, cellCenterData);

    vector<int> probesCellID = unstructGrid->GetZoneProbesCellID();

    RDouble *q_probes = new RDouble [probesLocalNumber];
    SetField(q_probes, 0.0, probesLocalNumber);

    for (int iProbe = 0; iProbe < probesLocalNumber; ++ iProbe)
    {
        vector<int> interpolationNodesID;
        vector<RDouble> numerator;

        int probeCell = probesCellID[iProbe];

        for (int iNode = 0; iNode < nodeNumberOfEachCell[probeCell]; ++ iNode)
        {
            int node = cell2Node[probeCell][iNode];
            NecessaryVarForInterpolation(unstructGrid, iProbe, node, numerator, interpolationNodesID, false);
        }

        int nNodes = static_cast<int>(interpolationNodesID.size());
        RDouble denominator = SUMMARY(numerator, nNodes);

        for (int iNode = 0; iNode < nNodes; ++ iNode)
        {
            RDouble weightCoef = numerator[iNode] / denominator;
            int nodeID = interpolationNodesID[iNode];
            q_probes[iProbe] += qn[nodeID] * weightCoef;
        }
    }

    probeVariablesPtr->UpdateDataPtr(name, q_probes);

    delete [] qn;
}

LIB_EXPORT void Post_Probes::UpdateProbesVarPtr(const string &name, RDouble4D &cellCenterData, int index)
{
    if (!grid)
    {
        TK_Exit::ExceptionExit("Post_Probe::UpdateProbesVarPtr needs Geometry infomation.");
    }

    int probeVariablesInterpolationMethod = GlobalDataBase::GetIntParaFromDB("probeVariablesInterpolationMethod");

    if (probeVariablesInterpolationMethod == CELLVALUE)
    {
        UpdateProbesVarPtrM1(name, cellCenterData, index);
    }
    else if (probeVariablesInterpolationMethod == CELLSTOPROBE)
    {
        UpdateProbesVarPtrM2(name, cellCenterData, index);
    }
    else if (probeVariablesInterpolationMethod == NODESTOPROBE)
    {
        UpdateProbesVarPtrM3(name, cellCenterData, index);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("probeVariablesInterpolationMethod", probeVariablesInterpolationMethod);
    }
}

LIB_EXPORT void Post_Probes::UpdateProbesVarPtrM1(const string &name, RDouble4D &cellCenterData, int index)
{
    StructGrid *structGrid = StructGridCast(grid);

    int probesLocalNumber = structGrid->GetZoneProbesNumber();
    if (probesLocalNumber == 0) return;

    vector<int> probesCellNI = structGrid->GetZoneProbesCellNI();
    vector<int> probesCellNJ = structGrid->GetZoneProbesCellNJ();
    vector<int> probesCellNK = structGrid->GetZoneProbesCellNK();

    RDouble *q_probes = new RDouble [probesLocalNumber];
    SetField(q_probes, 0.0, probesLocalNumber);

    for (int iProbe = 0; iProbe < probesLocalNumber; ++ iProbe)
    {
        int probeNI = probesCellNI[iProbe];
        int probeNJ = probesCellNJ[iProbe];
        int probeNK = probesCellNK[iProbe];
        q_probes[iProbe] = cellCenterData(probeNI, probeNJ, probeNK, index);
    }

    probeVariablesPtr->UpdateDataPtr(name, q_probes);
}

LIB_EXPORT void Post_Probes::UpdateProbesVarPtrM2(const string &name, RDouble4D &cellCenterData, int index)
{
    StructGrid *structGrid = StructGridCast(grid);

    int probesLocalNumber = structGrid->GetZoneProbesNumber();
    if (probesLocalNumber == 0) return;

    vector<int> probesCellNI = structGrid->GetZoneProbesCellNI();
    vector<int> probesCellNJ = structGrid->GetZoneProbesCellNJ();
    vector<int> probesCellNK = structGrid->GetZoneProbesCellNK();

    RDouble *q_probes = new RDouble [probesLocalNumber];
    SetField(q_probes, 0.0, probesLocalNumber);

    int nDim = GetDim();

    for (int iProbe = 0; iProbe < probesLocalNumber; ++ iProbe)
    {
        vector<vector<int> > interpolationCellsID;
        vector<RDouble> numerator;

        int probeNI = probesCellNI[iProbe];
        int probeNJ = probesCellNJ[iProbe];
        int probeNK = probesCellNK[iProbe];

        NecessaryVarForInterpolation(structGrid, iProbe, probeNI, probeNJ, probeNK, numerator, interpolationCellsID, true);

        for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
        {
            int i1, j1, k1;
            GetNsurfIndex(iSurface, i1, j1, k1);

            int il = probeNI - i1;
            int jl = probeNJ - j1;
            int kl = probeNK - k1;

            int ir = probeNI + i1;
            int jr = probeNJ + j1;
            int kr = probeNK + k1;

            NecessaryVarForInterpolation(structGrid, iProbe, il, jl, kl, numerator, interpolationCellsID, true);
            NecessaryVarForInterpolation(structGrid, iProbe, ir, jr, kr, numerator, interpolationCellsID, true);
        }

        int nCells = static_cast<int>(interpolationCellsID.size());
        RDouble denominator = SUMMARY(numerator, nCells);

        for (int iCell = 0; iCell < nCells; ++ iCell)
        {
            RDouble weightCoef = numerator[iCell] / denominator;
            int cellID_i = interpolationCellsID[iCell][0];
            int cellID_j = interpolationCellsID[iCell][1];
            int cellID_k = interpolationCellsID[iCell][2];
            q_probes[iProbe] += cellCenterData(cellID_i, cellID_j, cellID_k, index) * weightCoef;
        }
    }

    probeVariablesPtr->UpdateDataPtr(name, q_probes);
}

LIB_EXPORT void Post_Probes::UpdateProbesVarPtrM3(const string &name, RDouble4D &cellCenterData, int index)
{
    StructGrid *structGrid = StructGridCast(grid);

    int probesLocalNumber = structGrid->GetZoneProbesNumber();
    if (probesLocalNumber == 0) return;

    int ni = structGrid->GetNI();
    int nj = structGrid->GetNJ();
    int nk = structGrid->GetNK();

    Range I(-1, ni+1);
    Range J(-1, nj+1);
    Range K(-1, nk+1);
    if (nk == 1) K.setRange(1, 1);
    Range M(0, 0);

    RDouble4D *qn = new RDouble4D(I, J, K, M, fortranArray);
    CompNodeVar(structGrid, *qn, 0, cellCenterData, index);

    vector<int> probesCellNI = structGrid->GetZoneProbesCellNI();
    vector<int> probesCellNJ = structGrid->GetZoneProbesCellNJ();
    vector<int> probesCellNK = structGrid->GetZoneProbesCellNK();

    RDouble *q_probes = new RDouble [probesLocalNumber];
    SetField(q_probes, 0.0, probesLocalNumber);

    int i1 = 1;
    int j1 = 1;
    int k1 = 1;
    if (nk == 1) k1 = 0;

    for (int iProbe = 0; iProbe < probesLocalNumber; ++ iProbe)
    {
        vector<vector<int> > interpolationNodesID;
        vector<RDouble> numerator;

        int probeNI = probesCellNI[iProbe];
        int probeNJ = probesCellNJ[iProbe];
        int probeNK = probesCellNK[iProbe];

        for (int k = probeNK; k <= probeNK + k1; ++ k)
        {
            for (int j = probeNJ; j <= probeNJ + j1; ++ j)
            {
                for (int i = probeNI; i <= probeNI + i1; ++ i)
                {
                    NecessaryVarForInterpolation(structGrid, iProbe, i, j, k, numerator, interpolationNodesID, false);
                }
            }
        }

        int nNodes = static_cast<int>(interpolationNodesID.size());
        RDouble denominator = SUMMARY(numerator, nNodes);

        for (int iNode = 0; iNode < nNodes; ++ iNode)
        {
            RDouble weightCoef = numerator[iNode] / denominator;
            int nodeID_i = interpolationNodesID[iNode][0];
            int nodeID_j = interpolationNodesID[iNode][1];
            int nodeID_k = interpolationNodesID[iNode][2];
            q_probes[iProbe] += (*qn)(nodeID_i, nodeID_j, nodeID_k, 0) * weightCoef;
        }
    }

    probeVariablesPtr->UpdateDataPtr(name, q_probes);
    delete qn;
}

LIB_EXPORT void Post_Probes::NecessaryVarForInterpolation(UnstructGrid *grid, int iprobe, int index, vector<RDouble> &numerator, vector<int> &interpolationIndexs, bool flag)
{
    vector<vector<RDouble> > probesCoordinates = grid->GetZoneProbesCoordinates();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *xn = grid->GetX();
    RDouble *yn = grid->GetY();
    RDouble *zn = grid->GetZ();

    RDouble distance = 0;
    RDouble odistance = 0;

    //! The range of r is from 0.5 to 3.0.
    const int r = 1;

    if (flag)
    {
        distance = DISTANCE(xcc[index] - probesCoordinates[iprobe][0], ycc[index] - probesCoordinates[iprobe][1], zcc[index] - probesCoordinates[iprobe][2]) + TINY;
    }
    else
    {
        distance = DISTANCE(xn[index] - probesCoordinates[iprobe][0], yn[index] - probesCoordinates[iprobe][1], zn[index] - probesCoordinates[iprobe][2]) + TINY;
    }

    odistance = pow(1.0 / distance, r);
    numerator.push_back(odistance);
    interpolationIndexs.push_back(index);
}

LIB_EXPORT void Post_Probes::NecessaryVarForInterpolation(StructGrid *grid, int iprobe, int i, int j, int k, vector<RDouble> &numerator, vector<vector<int> > &interpolationIndexs, bool flag)
{
    vector<vector<RDouble> > probesCoordinates = grid->GetZoneProbesCoordinates();
    RDouble3D &xcc = *(grid->GetCellCenterX());
    RDouble3D &ycc = *(grid->GetCellCenterY());
    RDouble3D &zcc = *(grid->GetCellCenterZ());

    RDouble3D &xn = * grid->GetStructX();
    RDouble3D &yn = * grid->GetStructY();
    RDouble3D &zn = * grid->GetStructZ();

    RDouble distance = 0;
    RDouble odistance = 0;
    vector<int> index;

    //! The range of r is from 0.5 to 3.0.
    const int r = 1;

    if (flag)
    {
        distance = DISTANCE(xcc(i, j, k) - probesCoordinates[iprobe][0], ycc(i, j, k) - probesCoordinates[iprobe][1], zcc(i, j, k) - probesCoordinates[iprobe][2]) + TINY;
    }
    else
    {
        distance = DISTANCE(xn(i, j, k) - probesCoordinates[iprobe][0], yn(i, j, k) - probesCoordinates[iprobe][1], zn(i, j, k) - probesCoordinates[iprobe][2]) + TINY;
    }

    odistance = pow(1.0 / distance, r);
    numerator.push_back(odistance);
    index.push_back(i);
    index.push_back(j);
    index.push_back(k);
    interpolationIndexs.push_back(index);
    index.resize(0);
}

LIB_EXPORT void Post_Probes::GetAllProbesVarPtr(RDouble **q_probes)
{
    int varCount = 0;

    for (set<int>::iterator varIter = probeVariables.begin(); varIter != probeVariables.end(); ++ varIter)
    {
        int variableType = *varIter;
        string varName = GetProbesVariableName(variableType);
        q_probes[varCount ++] = (RDouble *)GetProbeVarPtr(varName);
    }
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        using namespace GAS_SPACE;

        int nm    = GlobalDataBase::GetIntParaFromDB("nm");
        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int neqn  = nl + nchem;

        string *varname = gas->GetNameOfSpecies();
        for (int m = nm; m < neqn; ++ m)
        {
            string varName = "massfraction-" + varname[m-nm];
            q_probes[varCount++] = (RDouble *)GetProbeVarPtr(varName);
        }
        for (int m = nm; m < neqn; ++ m)
        {
            string varName = "molefraction-" + varname[m-nm];
            q_probes[varCount++] = (RDouble *)GetProbeVarPtr(varName);
        }
    }
}

void Post_Probes::FreeAllProbeVarPtr()
{
    for (set<int>::iterator varIter = probeVariables.begin(); varIter != probeVariables.end(); ++ varIter)
    {
        int variableType = *varIter;
        string varName = GetProbesVariableName(variableType);
        RDouble *q_probes = (RDouble *)GetProbeVarPtr(varName);
        if (q_probes)
        {
            delete [] q_probes;
        }
    }
}

DumpToFileASCII::DumpToFileASCII(vector<DataContainer *> datalist_in, string filename_in)
{
    this->dataList = datalist_in;
    this->fileName = filename_in;

    for (unsigned int iZone = 0; iZone < dataList.size(); ++ iZone)
    {
        dataList[iZone]->MoveToBegin();
    }
    nVariables = 0;
}

DumpToFileASCII::~DumpToFileASCII()
{

}

void DumpToFileASCII::InitProbesVariables()
{
    int nProbeVariables = 6;
    int probeVariables[100] = {0, 1, 2, 3, 4, 5};
    nProbeVariables = GlobalDataBase::GetIntParaFromDB("nProbeVariables");
    GlobalDataBase::GetData("probeVariables", probeVariables, PHINT, nProbeVariables);

    Post_Probes *postProbes = new Post_Probes(nProbeVariables, probeVariables);

    int nVariablesToDump = postProbes->GetNumberofProbeVariables();

    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        using namespace GAS_SPACE;

        int nm    = GlobalDataBase::GetIntParaFromDB("nm");
        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int neqn  = nl + nchem;

        nVariablesToDump = nVariablesToDump + 2 * (neqn - nm);
    }
    SetnVariablesToDump(nVariablesToDump);

    vector<string> title_tecplot;
    postProbes->GetASCIIFileHeader(title_tecplot);
    SetVariables(title_tecplot.front());

    FreePointer(postProbes);
}

void DumpToFileASCII::Run()
{
    InitProbesVariables();
    WriteFlowFile();
}

void DumpToFileASCII::WriteFlowFile()
{
    int dataMonitorType = GlobalDataBase::GetIntParaFromDB("dataMonitorType");
    int samplefileMode = GlobalDataBase::GetIntParaFromDB("samplefileMode");
    int outerStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");

    if (dataMonitorType == PROBESMONITOR)
    {
        if (samplefileMode != 0)
        {
            WriteProbesFlowFile(outerStep);
        }
        else
        {
            WriteProbesFlowFile();
        }
    }
    else if (dataMonitorType == LINESMONITOR)
    {
        if (samplefileMode != 0)
        {
            WriteLinesFlowFile(outerStep);
        }
        else
        {
            WriteLinesFlowFile();
        }
    }
    else if (dataMonitorType == SURFACESMONITOR)
    {
        if (samplefileMode != 0)
        {
            WriteSurfacesFlowFile(outerStep);
        }
        else
        {
            WriteSurfacesFlowFile();
        }
    }
    else
    {
        TK_Exit::UnexpectedVarValue("dataMonitorType", dataMonitorType);
    }
}

void DumpToFileASCII::WriteProbesFlowFile()
{
    uint_t nzones = dataList.size();

    for (int iZone = 0; iZone < nzones; ++ iZone)
    {
        DataContainer *cdata = dataList[iZone];

        int probesLocalNumber;
        PHRead(cdata, &probesLocalNumber, 1);

        if (probesLocalNumber == 0) continue;

        vector<vector<RDouble> > probesCoordinates;
        vector<int> probesGlobalID, probesLineID, probesSurfaceID;
        probesCoordinates.resize(probesLocalNumber);
        probesGlobalID.resize(probesLocalNumber);
        probesLineID.resize(probesLocalNumber);
        probesSurfaceID.resize(probesLocalNumber);
        int nDim;
        RDouble **qProbes = NewPointer2<RDouble> (nVariables, probesLocalNumber);

        for (int iProbe = 0; iProbe < probesLocalNumber; ++ iProbe)
        {
            if (GetFilename() == "")
            {
                return;
            }
            ostringstream ossFileName;
            string finalFileName = GetFilename();

            PHRead(cdata, &nDim, 1);
            PHRead(cdata, &probesGlobalID[iProbe], 1);
            PHRead(cdata, &probesLineID[iProbe], 1);
            PHRead(cdata, &probesSurfaceID[iProbe], 1);
            probesCoordinates[iProbe].resize(nDim);

            for (int iDim = 0; iDim < nDim; ++ iDim)
            {
                PHRead(cdata, &probesCoordinates[iProbe][iDim], 1);
                finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, '_', probesCoordinates[iProbe][iDim]);
            }

            if (PHMPI::IsParallelRun())
            {
                int fileID = PHMPI::GetFileIndexofCurrentProcessor();
                finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, '_', fileID);
            }

            ossFileName << finalFileName;

            string probesFlowFileName = ossFileName.str();
            fstream probesFlowFile;
            OpenFile(probesFlowFile, probesFlowFileName, ios_base::out|ios_base::app);
            probesFlowFile << setiosflags(ios::left);
            probesFlowFile << setiosflags(ios::scientific);
            probesFlowFile << setprecision(10);

            if (IfFileEmpty(probesFlowFile))
            {
                probesFlowFile << variables << endl;
            }

            int wordWidth = 8;

            int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");

            probesFlowFile << setw(wordWidth) << outIterStep << "	";

            for (int iVariable = 0; iVariable < nVariables; ++ iVariable)
            {
                PHRead(cdata, qProbes[iVariable][iProbe]);
                probesFlowFile << setw(wordWidth) << qProbes[iVariable][iProbe] << "	";
            }

            int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
            if (isUnsteady)
            {
                int outerStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
                RDouble physicalTimeStepDimensional = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStepDimensional");
                RDouble real_time = physicalTimeStepDimensional * outerStep;
                probesFlowFile << setw(wordWidth) << real_time << endl;
            }
            else
            {
                probesFlowFile << endl;
            }
            CloseFile(probesFlowFile);
        }

        DelPointer2(qProbes);
    }
}

void DumpToFileASCII::WriteLinesFlowFile()
{
    int nLines = GlobalDataBase::GetIntParaFromDB("nLines");
    int *nProbesOfLine = reinterpret_cast< int * > (GlobalDataBase::GetDataPtr("nProbesOfLine"));
    vector<vector<RDouble> > *probesReorderCoordinates = new vector<vector<RDouble> > [nLines];
    vector<vector<RDouble> > *qGlobalProbes = new vector<vector<RDouble> > [nLines];

    for (int iLine = 0; iLine < nLines; ++ iLine)
    {
        int nProbes = nProbesOfLine[iLine];
        probesReorderCoordinates[iLine].resize(nProbes);
        qGlobalProbes[iLine].resize(nProbes);
    }

    uint_t nzones = dataList.size();
    for (int iZone = 0; iZone < nzones; ++ iZone)
    {
        DataContainer *cdata = dataList[iZone];

        int probesLocalNumber;
        PHRead(cdata, &probesLocalNumber, 1);

        if (probesLocalNumber == 0) continue;

        vector<vector<RDouble> > probesCoordinates;
        vector<int> probesGlobalID, probesLineID, probesSurfaceID;
        probesCoordinates.resize(probesLocalNumber);
        probesGlobalID.resize(probesLocalNumber);
        probesLineID.resize(probesLocalNumber);
        probesSurfaceID.resize(probesLocalNumber);
        int nDim;
        RDouble **qProbes = NewPointer2<RDouble> (nVariables, probesLocalNumber);

        for (int iProbe = 0; iProbe < probesLocalNumber; ++ iProbe)
        {
            PHRead(cdata, &nDim, 1);
            PHRead(cdata, &probesGlobalID[iProbe], 1);
            PHRead(cdata, &probesLineID[iProbe], 1);
            PHRead(cdata, &probesSurfaceID[iProbe], 1);
            probesCoordinates[iProbe].resize(nDim);

            int iProbeLineID = probesLineID[iProbe];
            int iProbeGlobalID = probesGlobalID[iProbe];
            probesReorderCoordinates[iProbeLineID][iProbeGlobalID].resize(nDim);
            qGlobalProbes[iProbeLineID][iProbeGlobalID].resize(nVariables);

            for (int iDim = 0; iDim < nDim; ++ iDim)
            {
                PHRead(cdata, &probesCoordinates[iProbe][iDim], 1);
                probesReorderCoordinates[iProbeLineID][iProbeGlobalID][iDim] = probesCoordinates[iProbe][iDim];
            }

            for (int iVariable = 0; iVariable < nVariables; ++ iVariable)
            {
                PHRead(cdata, qProbes[iVariable][iProbe]);
                qGlobalProbes[iProbeLineID][iProbeGlobalID][iVariable] = qProbes[iVariable][iProbe];
            }
        }

        DelPointer2(qProbes);
    }

    for (int iLine = 0; iLine < nLines; ++ iLine)
    {
        int nProbes = nProbesOfLine[iLine];

        if (GetFilename() == "")
        {
            return;
        }

        ostringstream ossFileName;
        string finalFileName = GetFilename();
        finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, "_", "line", iLine);

        if (PHMPI::IsParallelRun())
        {
            int fileID = PHMPI::GetFileIndexofCurrentProcessor();
            finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, '_', fileID);
        }

        ossFileName << finalFileName;

        string linesFlowFileName = ossFileName.str();
        fstream linesFlowFile;
        OpenFile(linesFlowFile, linesFlowFileName, ios_base::out|ios_base::app);
        linesFlowFile << setiosflags(ios::left);
        linesFlowFile << setiosflags(ios::scientific);
        linesFlowFile << setprecision(10);

        if (IfFileEmpty(linesFlowFile))
        {
            linesFlowFile << variables << endl;
        }

        int wordWidth = 8;

        for (int iProbe = 0; iProbe < nProbes; ++ iProbe)
        {
            int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
            linesFlowFile << setw(wordWidth) << outIterStep << "	";

            for (int iDim = 0; iDim < probesReorderCoordinates[iLine][iProbe].size(); ++ iDim)
            {
                linesFlowFile << setw(wordWidth) << probesReorderCoordinates[iLine][iProbe][iDim] << "	";
            }

            for (int iVariable = 0; iVariable < nVariables; ++ iVariable)
            {
                linesFlowFile << setw(wordWidth) << qGlobalProbes[iLine][iProbe][iVariable] << "	";
            }

            int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
            if (isUnsteady)
            {
                int outerStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
                RDouble physicalTimeStepDimensional = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStepDimensional");
                RDouble real_time = physicalTimeStepDimensional * outerStep;
                linesFlowFile << setw(wordWidth) << real_time << endl;
            }
            else
            {
                linesFlowFile << endl;
            }
        }

        CloseFile(linesFlowFile);
    }

    delete [] probesReorderCoordinates;
    delete [] qGlobalProbes;
    probesReorderCoordinates = NULL;
    qGlobalProbes = NULL;
}

void DumpToFileASCII::WriteSurfacesFlowFile()
{
    int nSurfaces= GlobalDataBase::GetIntParaFromDB("nSurface");
    int *nProbesOfSurface = reinterpret_cast< int * > (GlobalDataBase::GetDataPtr("nProbesOfSurface"));

    vector<vector<RDouble> > *probesReorderCoordinates = new vector<vector<RDouble> > [nSurfaces];
    vector<vector<RDouble> > *qGlobalProbes = new vector<vector<RDouble> > [nSurfaces];

    for (int iSurface = 0; iSurface < nSurfaces; ++ iSurface)
    {
        int nProbes = nProbesOfSurface[iSurface];
        probesReorderCoordinates[iSurface].resize(nProbes);
        qGlobalProbes[iSurface].resize(nProbes);
    }

    uint_t nzones = dataList.size();
    for (int iZone = 0; iZone < nzones; ++ iZone)
    {
        DataContainer *cdata = dataList[iZone];

        int probesLocalNumber;
        PHRead(cdata, &probesLocalNumber, 1);

        if (probesLocalNumber == 0) continue;

        vector<vector<RDouble> > probesCoordinates;
        vector<int> probesGlobalID, probesLineID, probesSurfaceID;
        probesCoordinates.resize(probesLocalNumber);
        probesGlobalID.resize(probesLocalNumber);
        probesLineID.resize(probesLocalNumber);
        probesSurfaceID.resize(probesLocalNumber);
        int nDim;
        RDouble **qProbes = NewPointer2<RDouble> (nVariables, probesLocalNumber);

        for (int iProbe = 0; iProbe < probesLocalNumber; ++ iProbe)
        {
            PHRead(cdata, &nDim, 1);
            PHRead(cdata, &probesGlobalID[iProbe], 1);
            PHRead(cdata, &probesLineID[iProbe], 1);
            PHRead(cdata, &probesSurfaceID[iProbe], 1);
            probesCoordinates[iProbe].resize(nDim);

            int iProbeSurfaceID = probesSurfaceID[iProbe];
            int iProbeGlobalID = probesGlobalID[iProbe];
            probesReorderCoordinates[iProbeSurfaceID][iProbeGlobalID].resize(nDim);
            qGlobalProbes[iProbeSurfaceID][iProbeGlobalID].resize(nVariables);

            for (int iDim = 0; iDim < nDim; ++ iDim)
            {
                PHRead(cdata, &probesCoordinates[iProbe][iDim], 1);
                probesReorderCoordinates[iProbeSurfaceID][iProbeGlobalID][iDim] = probesCoordinates[iProbe][iDim];
            }

            for (int iVariable = 0; iVariable < nVariables; ++ iVariable)
            {
                PHRead(cdata, qProbes[iVariable][iProbe]);
                qGlobalProbes[iProbeSurfaceID][iProbeGlobalID][iVariable] = qProbes[iVariable][iProbe];
            }
        }

        DelPointer2(qProbes);
    }

    for (int iSurface = 0; iSurface < nSurfaces; ++ iSurface)
    {
        int nProbes = nProbesOfSurface[iSurface];

        if (GetFilename() == "")
        {
            return;
        }

        ostringstream ossFileName;
        string finalFileName = GetFilename();
        finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, "_", "surface", iSurface);

        if (PHMPI::IsParallelRun())
        {
            int fileID = PHMPI::GetFileIndexofCurrentProcessor();
            finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, '_', fileID);
        }

        ossFileName << finalFileName;

        string surfacesFlowFileName = ossFileName.str();
        fstream surfacesFlowFile;
        OpenFile(surfacesFlowFile, surfacesFlowFileName, ios_base::out|ios_base::app);
        surfacesFlowFile << setiosflags(ios::left);
        surfacesFlowFile << setiosflags(ios::scientific);
        surfacesFlowFile << setprecision(10);

        if (IfFileEmpty(surfacesFlowFile))
        {
            surfacesFlowFile << variables << endl;
        }

        int wordWidth = 8;

        for (int iProbe = 0; iProbe < nProbes; ++ iProbe)
        {
            int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
            surfacesFlowFile << setw(wordWidth) << outIterStep << "	";

            for (int iDim = 0; iDim < probesReorderCoordinates[iSurface][iProbe].size(); ++ iDim)
            {
                surfacesFlowFile << setw(wordWidth) << probesReorderCoordinates[iSurface][iProbe][iDim] << "	";
            }

            for (int iVariable = 0; iVariable < nVariables; ++ iVariable)
            {
                surfacesFlowFile << setw(wordWidth) << qGlobalProbes[iSurface][iProbe][iVariable] << "	";
            }

            int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
            if (isUnsteady)
            {
                int outerStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
                RDouble physicalTimeStepDimensional = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStepDimensional");
                RDouble real_time = physicalTimeStepDimensional * outerStep;
                surfacesFlowFile << setw(wordWidth) << real_time << endl;
            }
            else
            {
                surfacesFlowFile << endl;
            }

        }

        CloseFile(surfacesFlowFile);
    }

    delete [] probesReorderCoordinates;
    delete [] qGlobalProbes;
    probesReorderCoordinates = NULL;
    qGlobalProbes = NULL;
}

void DumpToFileASCII::WriteProbesFlowFile(int outerStep)
{
    int nTotalProbes = GlobalDataBase::GetIntParaFromDB("nTotalProbes");
    vector<vector<RDouble> > probesReorderCoordinates;
    probesReorderCoordinates.resize(nTotalProbes);
    vector<vector<RDouble> > qGlobalProbes;
    qGlobalProbes.resize(nTotalProbes);

    uint_t nzones = dataList.size();
    for (int iZone = 0; iZone < nzones; ++ iZone)
    {
        DataContainer *cdata = dataList[iZone];

        int probesLocalNumber;
        PHRead(cdata, &probesLocalNumber, 1);

        if (probesLocalNumber == 0) continue;

        vector<vector<RDouble> > probesCoordinates;
        vector<int> probesGlobalID, probesLineID, probesSurfaceID;
        probesCoordinates.resize(probesLocalNumber);
        probesGlobalID.resize(probesLocalNumber);
        probesLineID.resize(probesLocalNumber);
        probesSurfaceID.resize(probesLocalNumber);
        int nDim;
        RDouble **qProbes = NewPointer2<RDouble> (nVariables, probesLocalNumber);

        for (int iProbe = 0; iProbe < probesLocalNumber; ++ iProbe)
        {
            PHRead(cdata, &nDim, 1);
            PHRead(cdata, &probesGlobalID[iProbe], 1);
            PHRead(cdata, &probesLineID[iProbe], 1);
            PHRead(cdata, &probesSurfaceID[iProbe], 1);
            probesCoordinates[iProbe].resize(nDim);

            int iProbeGlobalID =probesGlobalID[iProbe];
            probesReorderCoordinates[iProbeGlobalID].resize(nDim);
            qGlobalProbes[iProbeGlobalID].resize(nVariables);

            for (int iDim = 0; iDim < nDim; ++ iDim)
            {
                PHRead(cdata, &probesCoordinates[iProbe][iDim], 1);
                probesReorderCoordinates[iProbeGlobalID][iDim] = probesCoordinates[iProbe][iDim];
            }

            for (int iVariable = 0; iVariable < nVariables; ++ iVariable)
            {
                PHRead(cdata, qProbes[iVariable][iProbe]);
                qGlobalProbes[iProbeGlobalID][iVariable] = qProbes[iVariable][iProbe];
            }
        }

        DelPointer2(qProbes);
    }

    for (int iProbe = 0; iProbe < nTotalProbes; ++ iProbe)
    {
        if (GetFilename() == "")
        {
            return;
        }
        ostringstream ossFileName;
        string finalFileName = GetFilename();

        finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, "_", "probes");
        finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, '_', outerStep);
        if (PHMPI::IsParallelRun())
        {
            int fileID = PHMPI::GetFileIndexofCurrentProcessor();
            finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, '_', fileID);
        }

        ossFileName << finalFileName;

        string probesFlowFileName = ossFileName.str();
        fstream probesFlowFile;
        OpenFile(probesFlowFile, probesFlowFileName, ios_base::out|ios_base::app);
        probesFlowFile << setiosflags(ios::left);
        probesFlowFile << setiosflags(ios::scientific);
        probesFlowFile << setprecision(10);

        if (IfFileEmpty(probesFlowFile))
        {
            probesFlowFile << variables << endl;
        }

        int wordWidth = 8;

        probesFlowFile << setw(wordWidth) << outerStep << "	";

        for (int iDim = 0; iDim < probesReorderCoordinates[iProbe].size(); ++ iDim)
        {
            probesFlowFile << setw(wordWidth) << probesReorderCoordinates[iProbe][iDim] << "	";
        }

        for (int iVariable = 0; iVariable < nVariables; ++ iVariable)
        {
            probesFlowFile << setw(wordWidth) << qGlobalProbes[iProbe][iVariable] << "	";
        }

        int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
        if (isUnsteady)
        {
            int outerStep1 = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
            RDouble physicalTimeStepDimensional = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStepDimensional");
            RDouble real_time = physicalTimeStepDimensional * outerStep1;
            probesFlowFile << setw(wordWidth) << real_time << endl;
        }
        else
        {
            probesFlowFile << endl;
        }

        CloseFile(probesFlowFile);
    }
}

void DumpToFileASCII::WriteLinesFlowFile(int outerStep)
{
    int nLines = GlobalDataBase::GetIntParaFromDB("nLines");
    int *nProbesOfLine = reinterpret_cast< int * > (GlobalDataBase::GetDataPtr("nProbesOfLine"));
    vector<vector<RDouble> > *probesReorderCoordinates = new vector<vector<RDouble> > [nLines];
    vector<vector<RDouble> > *qGlobalProbes = new vector<vector<RDouble> > [nLines];

    for (int iLine = 0; iLine < nLines; ++ iLine)
    {
        int nProbes = nProbesOfLine[iLine];
        probesReorderCoordinates[iLine].resize(nProbes);
        qGlobalProbes[iLine].resize(nProbes);
    }

    uint_t nzones = dataList.size();
    for (int iZone = 0; iZone < nzones; ++ iZone)
    {
        DataContainer *cdata = dataList[iZone];

        int probesLocalNumber;
        PHRead(cdata, &probesLocalNumber, 1);

        if (probesLocalNumber == 0) continue;

        vector<vector<RDouble> > probesCoordinates;
        vector<int> probesGlobalID, probesLineID, probesSurfaceID;
        probesCoordinates.resize(probesLocalNumber);
        probesGlobalID.resize(probesLocalNumber);
        probesLineID.resize(probesLocalNumber);
        probesSurfaceID.resize(probesLocalNumber);
        int nDim;
        RDouble **qProbes = NewPointer2<RDouble> (nVariables, probesLocalNumber);

        for (int iProbe = 0; iProbe < probesLocalNumber; ++ iProbe)
        {
            PHRead(cdata, &nDim, 1);
            PHRead(cdata, &probesGlobalID[iProbe], 1);
            PHRead(cdata, &probesLineID[iProbe], 1);
            PHRead(cdata, &probesSurfaceID[iProbe], 1);
            probesCoordinates[iProbe].resize(nDim);

            int iProbeLineID = probesLineID[iProbe];
            int iProbeGlobalID = probesGlobalID[iProbe];
            probesReorderCoordinates[iProbeLineID][iProbeGlobalID].resize(nDim);
            qGlobalProbes[iProbeLineID][iProbeGlobalID].resize(nVariables);

            for (int iDim = 0; iDim < nDim; ++ iDim)
            {
                PHRead(cdata, &probesCoordinates[iProbe][iDim], 1);
                probesReorderCoordinates[iProbeLineID][iProbeGlobalID][iDim] = probesCoordinates[iProbe][iDim];
            }

            for (int iVariable = 0; iVariable < nVariables; ++ iVariable)
            {
                PHRead(cdata, qProbes[iVariable][iProbe]);
                qGlobalProbes[iProbeLineID][iProbeGlobalID][iVariable] = qProbes[iVariable][iProbe];
            }
        }

        DelPointer2(qProbes);
    }

    for (int iLine = 0; iLine < nLines; ++ iLine)
    {
        int nProbes = nProbesOfLine[iLine];

        if (GetFilename() == "")
        {
            return;
        }

        ostringstream ossFileName;
        string finalFileName = GetFilename();
        finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, "_", "line", iLine);
        finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, '_', outerStep);

        if (PHMPI::IsParallelRun())
        {
            int fileID = PHMPI::GetFileIndexofCurrentProcessor();
            finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, '_', fileID);
        }

        ossFileName << finalFileName;

        string linesFlowFileName = ossFileName.str();
        fstream linesFlowFile;
        OpenFile(linesFlowFile, linesFlowFileName, ios_base::out|ios_base::app);
        linesFlowFile << setiosflags(ios::left);
        linesFlowFile << setiosflags(ios::scientific);
        linesFlowFile << setprecision(10);

        if (IfFileEmpty(linesFlowFile))
        {
            linesFlowFile << variables << endl;
        }

        int wordWidth = 8;

        for (int iProbe = 0; iProbe < nProbes; ++ iProbe)
        {
            linesFlowFile << setw(wordWidth) << outerStep << "	";

            for (int iDim = 0; iDim < probesReorderCoordinates[iLine][iProbe].size(); ++ iDim)
            {
                linesFlowFile << setw(wordWidth) << probesReorderCoordinates[iLine][iProbe][iDim] << "	";
            }

            for (int iVariable = 0; iVariable < nVariables; ++ iVariable)
            {
                linesFlowFile << setw(wordWidth) << qGlobalProbes[iLine][iProbe][iVariable] << "	";
            }

            int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
            if (isUnsteady)
            {
                int outerStep1 = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
                RDouble physicalTimeStepDimensional = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStepDimensional");
                RDouble real_time = physicalTimeStepDimensional * outerStep1;
                linesFlowFile << setw(wordWidth) << real_time << endl;
            }
            else
            {
                linesFlowFile << endl;
            }
        }

        CloseFile(linesFlowFile);
    }

    delete [] probesReorderCoordinates;
    delete [] qGlobalProbes;
    probesReorderCoordinates = NULL;
    qGlobalProbes = NULL;
}

void DumpToFileASCII::WriteSurfacesFlowFile(int outerStep)
{
    int nSurfaces= GlobalDataBase::GetIntParaFromDB("nSurface");
    int *nProbesOfSurface = reinterpret_cast< int * > (GlobalDataBase::GetDataPtr("nProbesOfSurface"));

    vector<vector<RDouble> > *probesReorderCoordinates = new vector<vector<RDouble> > [nSurfaces];
    vector<vector<RDouble> > *qGlobalProbes = new vector<vector<RDouble> > [nSurfaces];

    for (int iSurface = 0; iSurface < nSurfaces; ++ iSurface)
    {
        int nProbes = nProbesOfSurface[iSurface];
        probesReorderCoordinates[iSurface].resize(nProbes);
        qGlobalProbes[iSurface].resize(nProbes);
    }

    uint_t nzones = dataList.size();
    for (int iZone = 0; iZone < nzones; ++ iZone)
    {
        DataContainer *cdata = dataList[iZone];

        int probesLocalNumber;
        PHRead(cdata, &probesLocalNumber, 1);

        if (probesLocalNumber == 0) continue;

        vector<vector<RDouble> > probesCoordinates;
        vector<int> probesGlobalID, probesLineID, probesSurfaceID;
        probesCoordinates.resize(probesLocalNumber);
        probesGlobalID.resize(probesLocalNumber);
        probesLineID.resize(probesLocalNumber);
        probesSurfaceID.resize(probesLocalNumber);
        int nDim;
        RDouble **qProbes = NewPointer2<RDouble> (nVariables, probesLocalNumber);

        for (int iProbe = 0; iProbe < probesLocalNumber; ++ iProbe)
        {
            PHRead(cdata, &nDim, 1);
            PHRead(cdata, &probesGlobalID[iProbe], 1);
            PHRead(cdata, &probesLineID[iProbe], 1);
            PHRead(cdata, &probesSurfaceID[iProbe], 1);
            probesCoordinates[iProbe].resize(nDim);

            int iProbeSurfaceID = probesSurfaceID[iProbe];
            int iProbeGlobalID = probesGlobalID[iProbe];
            probesReorderCoordinates[iProbeSurfaceID][iProbeGlobalID].resize(nDim);
            qGlobalProbes[iProbeSurfaceID][iProbeGlobalID].resize(nVariables);

            for (int iDim = 0; iDim < nDim; ++ iDim)
            {
                PHRead(cdata, &probesCoordinates[iProbe][iDim], 1);
                probesReorderCoordinates[iProbeSurfaceID][iProbeGlobalID][iDim] = probesCoordinates[iProbe][iDim];
            }

            for (int iVariable = 0; iVariable < nVariables; ++ iVariable)
            {
                PHRead(cdata, qProbes[iVariable][iProbe]);
                qGlobalProbes[iProbeSurfaceID][iProbeGlobalID][iVariable] = qProbes[iVariable][iProbe];
            }
        }

        DelPointer2(qProbes);
    }

    for (int iSurface = 0; iSurface < nSurfaces; ++ iSurface)
    {
        int nProbes = nProbesOfSurface[iSurface];

        if (GetFilename() == "")
        {
            return;
        }

        ostringstream ossFileName;
        string finalFileName = GetFilename();
        finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, "_", "surface", iSurface);
        finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, '_', outerStep);

        if (PHMPI::IsParallelRun())
        {
            int fileID = PHMPI::GetFileIndexofCurrentProcessor();
            finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, '_', fileID);
        }

        ossFileName << finalFileName;

        string surfacesFlowFileName = ossFileName.str();
        fstream surfacesFlowFile;
        OpenFile(surfacesFlowFile, surfacesFlowFileName, ios_base::out|ios_base::app);
        surfacesFlowFile << setiosflags(ios::left);
        surfacesFlowFile << setiosflags(ios::scientific);
        surfacesFlowFile << setprecision(10);

        if (IfFileEmpty(surfacesFlowFile))
        {
            surfacesFlowFile << variables << endl;
        }

        int wordWidth = 8;

        for (int iProbe = 0; iProbe < nProbes; ++ iProbe)
        {
            surfacesFlowFile << setw(wordWidth) << outerStep << "	";

            for (int iDim = 0; iDim < probesReorderCoordinates[iSurface][iProbe].size(); ++ iDim)
            {
                surfacesFlowFile << setw(wordWidth) << probesReorderCoordinates[iSurface][iProbe][iDim] << "	";
            }

            for (int iVariable = 0; iVariable < nVariables; ++ iVariable)
            {
                surfacesFlowFile << setw(wordWidth) << qGlobalProbes[iSurface][iProbe][iVariable] << "	";
            }

            int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
            if (isUnsteady)
            {
                int outerStep1 = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
                RDouble physicalTimeStepDimensional = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStepDimensional");
                RDouble real_time = physicalTimeStepDimensional * outerStep1;
                surfacesFlowFile << setw(wordWidth) << real_time << endl;
            }
            else
            {
                surfacesFlowFile << endl;
            }
        }

        CloseFile(surfacesFlowFile);
    }

    delete [] probesReorderCoordinates;
    delete [] qGlobalProbes;
    probesReorderCoordinates = NULL;
    qGlobalProbes = NULL;
}

}
