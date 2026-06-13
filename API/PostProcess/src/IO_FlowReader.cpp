#include "IO_FlowReader.h"
#include "Glb_Dimension.h"
#include "IO_FileName.h"
#include "GridType.h"
#include "PHIO.h"

using namespace std;
namespace PHSPACE
{

IO_FlowReader::IO_FlowReader(Region *region_in)
{
    this->region = region_in;

    this->nEquation           = 0;
    this->nTemperatureModel   = 0;
    this->nTurbulenceEquation = 0;
    this->nTransitionEquation = 0;
}

IO_FlowReader::~IO_FlowReader()
{

}

void IO_FlowReader::Run()
{
    Init();
    ReadFlowFile();
    ReadTurbFile();
    ReadTransitionFile();
}

void IO_FlowReader::Init()
{
    using namespace PHMPI;

    restartNSFile = GlobalDataBase::GetStrParaFromDB("restartNSFile");
    turbFile = GlobalDataBase::GetStrParaFromDB("turbFile");
    transitionFile = GlobalDataBase::GetStrParaFromDB("transitionFile");

    int nChemical           = GlobalDataBase::GetIntParaFromDB("nchem");
    int nNSEquation         = GlobalDataBase::GetIntParaFromDB("nm");
    this->nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    this->nEquation         = nNSEquation + nChemical + nTemperatureModel - 1;

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    if (viscousType == ONE_EQU)
    {
        this->nTurbulenceEquation = 1;
    }
    else if (viscousType == TWO_EQU)
    {
        this->nTurbulenceEquation = 2;
        int transitionType = GlobalDataBase::GetIntParaFromDB("transitionType");
        if (transitionType > 0)
        {
            this->nTransitionEquation = transitionType;
        }
    }

    int nLocalZones = GetNumberofLocalZones();
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Grid *grid = GetGrid(zoneID);
        int gridType = grid->Type();
        if (gridType == STRUCTGRID)
        {
            InitMemoryStruct(grid);
        }
        else
        {
            InitMemoryUnstruct(grid);
        }
    }
}

void IO_FlowReader::InitMemoryStruct(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Range M(0, nEquation - 1);
    Range T(0, nTemperatureModel - 1);

    RDouble4D *q = new RDouble4D(I, J, K, M, fortranArray);
    grid->UpdateDataPtr("q", q);

    RDouble4D *t = new RDouble4D(I, J, K, T, fortranArray);
    grid->UpdateDataPtr("t", t);

    RDouble3D *visl = new RDouble3D(I, J, K, fortranArray);
    grid->UpdateDataPtr("visl", visl);

    RDouble3D *viscousTurbulent = new RDouble3D(I, J, K, fortranArray);
    grid->UpdateDataPtr("vist", viscousTurbulent);
    *viscousTurbulent = 0.0;

    if (nTurbulenceEquation > 0)
    {
        Range Mturb(0, nTurbulenceEquation - 1);
        RDouble4D *q_turb = new RDouble4D(I, J, K, Mturb, fortranArray);
        grid->UpdateDataPtr("q_turb", q_turb);
    }

    if (nTransitionEquation > 0)
    {
        Range Mturb(0, nTransitionEquation - 1);
        RDouble4D *qTransition = new RDouble4D(I, J, K, Mturb, fortranArray);
        grid->UpdateDataPtr("q_transition", qTransition);
    }
}

void IO_FlowReader::InitMemoryUnstruct(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q   = NewPointer2<RDouble>(nEquation, nTotal);
    grid->UpdateDataPtr("q", q);

    RDouble **gradPrimtiveVarX = NewPointer2<RDouble>(nEquation, nTotal);
    RDouble **gradPrimtiveVarY = NewPointer2<RDouble>(nEquation, nTotal);
    RDouble **gradPrimtiveVarZ = NewPointer2<RDouble>(nEquation, nTotal);

    grid->UpdateDataPtr("gradPrimtiveVarX", gradPrimtiveVarX);
    grid->UpdateDataPtr("gradPrimtiveVarY", gradPrimtiveVarY);
    grid->UpdateDataPtr("gradPrimtiveVarZ", gradPrimtiveVarZ);

    if (nTurbulenceEquation > 0)
    {
        RDouble **qTurbulence = NewPointer2<RDouble>(nTurbulenceEquation, nTotal);
        grid->UpdateDataPtr("q_turb", qTurbulence);
    }

    if (nTransitionEquation > 0)
    {
        RDouble **qTransition = NewPointer2<RDouble>(nTransitionEquation, nTotal);
        grid->UpdateDataPtr("q_transition", qTransition);
    }
}

void IO_FlowReader::ReadFlowFile()
{
    hid_t file;
    file = OpenHDF5File(restartNSFile);

    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Grid *grid = GetGrid(zoneID);
        int gridType = grid->Type();
        if (gridType == STRUCTGRID)
        {
            ReadFlowFileStr(file, grid);
        }
        else
        {
            ReadFlowFileUnStr(file, grid);
        }
    }

    H5Fclose(file);
}

void IO_FlowReader::ReadFlowFileUnStr(hid_t loc_id, Grid *gridIn)
{
    RDouble **q = reinterpret_cast< RDouble ** > ( gridIn->GetDataPtr("q") );

    hid_t grploc;
    string grpName;

    int outIterStep = 0;
    ReadData(loc_id, &outIterStep, "outnstep");
    GlobalDataBase::UpdateData("outnstep", &outIterStep, PHINT, 1);

    ostringstream oss;
    oss << "Group" << gridIn->GetZoneID();
    grpName = oss.str();

    grploc = OpenGroup(loc_id, grpName);
    ReadData(grploc, q[0], "q");

    H5Gclose(grploc);
}

void IO_FlowReader::ReadFlowFileStr(hid_t loc_id, Grid *gridIn)
{
    RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (gridIn->GetDataPtr("q"));

    StructGrid *strGrid = StructGridCast(gridIn);
    int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
    strGrid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

    hid_t grploc;
    string grpName;

    int outIterStep = 0;
    ReadData(loc_id, &outIterStep, "outnstep");
    GlobalDataBase::UpdateData("outnstep", &outIterStep, PHINT, 1);

    ostringstream oss;
    oss << "Group" << gridIn->GetZoneID();
    grpName = oss.str();

    grploc = OpenGroup(loc_id, grpName);
    ReadData(grploc, &primitiveVariables(iCellStart, jCellStart, kCellStart, 0), "q");

    H5Gclose(grploc);
}

void IO_FlowReader::ReadTurbFile()
{
    if (!nTurbulenceEquation)
    {
        return;
    }

    hid_t file;
    file = OpenHDF5File(turbFile);

    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Grid *grid = GetGrid(zoneID);
        int gridType = grid->Type();
        if (gridType == STRUCTGRID)
        {
            ReadTurbFileStr(file, grid);
        }
        else
        {
            ReadTurbFileUnStr(file, grid);
        }
    }

    H5Fclose(file);
}

void IO_FlowReader::ReadTurbFileStr(hid_t loc_id, Grid *gridIn)
{
    StructGrid *strGrid = StructGridCast(gridIn);
    int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
    strGrid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (strGrid->GetDataPtr("q_turb"));
    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (strGrid->GetDataPtr("vist"));

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss << "Group" << gridIn->GetZoneID();
    grpName = oss.str();

    grploc = OpenGroup(loc_id, grpName);
    ReadData(grploc, &qTurbulent(iCellStart, jCellStart, kCellStart, 0), "q_turb");
    ReadData(grploc, &viscousTurbulent(iCellStart, jCellStart, kCellStart), "visturb");

    H5Gclose(grploc);
}

void IO_FlowReader::ReadTurbFileUnStr(hid_t loc_id, Grid *gridIn)
{
    RDouble **qTurbulence     = reinterpret_cast<RDouble **> (gridIn->GetDataPtr("q_turb"));

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss << "Group" << gridIn->GetZoneID();
    grpName = oss.str();

    grploc = OpenGroup(loc_id, grpName);
    ReadData(grploc, qTurbulence[0], "q_turb");

    H5Gclose(grploc);
}

void IO_FlowReader::ReadTransitionFile()
{
    if (!nTransitionEquation)
    {
        return;
    }

    hid_t file;
    file = OpenHDF5File(transitionFile);

    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Grid *grid = GetGrid(zoneID);
        int gridType = grid->Type();
        if (gridType == STRUCTGRID)
        {
            ReadTransitionFileStr(file, grid);
        }
        else
        {
            ReadTransitionFileUnStr(file, grid);
        }
    }

    H5Fclose(file);
}

void IO_FlowReader::ReadTransitionFileStr(hid_t loc_id, Grid *gridIn)
{
    StructGrid *strGrid = StructGridCast(gridIn);
    int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
    strGrid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

    RDouble4D &qTransition = *reinterpret_cast<RDouble4D *> (strGrid->GetDataPtr("q_transition"));

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss << "Group" << gridIn->GetZoneID();
    grpName = oss.str();

    grploc = OpenGroup(loc_id, grpName);
    ReadData(grploc, &qTransition(iCellStart, jCellStart, kCellStart, 0), "q_transition");

    H5Gclose(grploc);
}

void IO_FlowReader::ReadTransitionFileUnStr(hid_t loc_id, Grid *gridIn)
{
    RDouble **qTransition     = reinterpret_cast<RDouble **> (gridIn->GetDataPtr("q_transition"));

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss << "Group" << gridIn->GetZoneID();
    grpName = oss.str();

    grploc = OpenGroup(loc_id, grpName);
    ReadData(grploc, qTransition[0], "q_transition");

    H5Gclose(grploc);
}

}