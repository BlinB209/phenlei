#include "PostProcess.h"
#include "ForceProcessStruct.h"
#include "ForceProcessUnstruct.h"
#include "GridType.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "Glb_Dimension.h"
#include <cmath>

using namespace std;
namespace PHSPACE
{

PostProcess::PostProcess(Region *region_in)
{
    this->region = region_in;
}

PostProcess::~PostProcess()
{

}

void PostProcess::Run()
{
    ReadGridFile();
    ReadFlowFile();

    int postTask = GlobalDataBase::GetIntParaFromDB("postTask");
    switch (postTask)
    {
    case COMPAIRFORCE:
        CompAirForceCoef();
        break;
    default:
        TK_Exit::UnexpectedVarValue("postTask = ", postTask);
        break;
    }
}

void PostProcess::ReadGridFile()
{
    region->ReadGrid();

    cout << endl << "Read grid finish ..." << endl;
}

void PostProcess::ComputeMetrics()
{
    using namespace PHMPI;

    ActionKey *actkey = new ActionKey();
    actkey->action = COMPUTEMETRICS;

    int myid = PHMPI::GetCurrentProcessorID();
    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int send_proc = GetZoneProcessorIDSepMode(iZone);
        if (myid == send_proc)
        {
            Grid *grid = GetGrid(iZone, 0);
            grid->ComputeMetrics(actkey);
        }
    }
    delete actkey;    actkey = nullptr;

}

void PostProcess::ReadFlowFile()
{
    cout << endl<< "Begin to read flow data ..." << endl << endl;

    IO_FlowReader *flowReader = new IO_FlowReader(this->region);
    flowReader->Run();

    delete flowReader;

    cout << "Read flow data finish ..." << endl;
}

void PostProcess::CompAirForceCoef()
{
    cout << endl << "Computer air force coefficient ..." << endl;

    InitGlobalValues();
    ReadBCFile();
    ComputeMetrics();

    Grid *grid = GetGrid(0);
    int gridtype = grid->Type();
    if (gridtype == STRUCTGRID)
    {
        ForceProcessStruct *forceprocess = new ForceProcessStruct();
        forceprocess->Run();

        delete forceprocess;
    }
    else
    {
        ForceProcessUnstruct *forceprocess = new ForceProcessUnstruct();
        forceprocess->Run();

        delete forceprocess;
    }

    cout << endl << "Computer finish ..." << endl;
}

void PostProcess::ReadBCFile()
{
    GlobalBoundaryCondition::ReadGlobalBoundaryCondition();

    GlobalBoundaryCondition::SetGlobalBCByGlobalDataBase();
}

void PostProcess::InitGlobalValues()
{
    ComputeCoefficientOfStateEquation();
    ComputeReferenceTsuth();
    ComputeReferenceLengthDimensional();

    int nm    = GlobalDataBase::GetIntParaFromDB("nm");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int neqn  = nm + nchem;

    int nl    = nm;
    GlobalDataBase::UpdateData("nl", &nl, PHINT, 1);

    RDouble *prim_inf = new RDouble[neqn];
    GlobalDataBase::UpdateDataPtr("prim_inf", prim_inf);

    RDouble attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    RDouble angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");
    RDouble coefficientOfStateEquation = GlobalDataBase::GetDoubleParaFromDB("coefficientOfStateEquation");

    RDouble attack   = attackd * PI / 180.0;
    RDouble sideslip = angleSlide * PI / 180.0;

    RDouble reference_density     = 1.0;
    RDouble reference_temperature = 1.0;

    RDouble reference_pressure, massoo_ref;
    massoo_ref = 1.0;

    reference_pressure = coefficientOfStateEquation * reference_density * reference_temperature / massoo_ref;

    using namespace IDX;

    prim_inf[IR] = reference_density;
    prim_inf[IU] = 1.0 * cos(attack) * cos(sideslip);
    prim_inf[IV] = 1.0 * sin(attack) * cos(sideslip);
    prim_inf[IW] = 1.0 *               sin(sideslip);
    prim_inf[IP] = reference_pressure;
}

void PostProcess::ComputeCoefficientOfStateEquation()
{
    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");

    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");

    RDouble coefficientOfStateEquation = 1.0 / (refGama * refMachNumber * refMachNumber);
    GlobalDataBase::UpdateData("coefficientOfStateEquation", &coefficientOfStateEquation, PHDOUBLE, 1);
}

void PostProcess::ComputeReferenceTsuth()
{
    RDouble tsi  = 110.4;

    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

    RDouble tsuth = tsi / refDimensionalTemperature;
    GlobalDataBase::UpdateData("tsuth", &tsuth, PHDOUBLE, 1);
}

void PostProcess::ComputeReferenceLengthDimensional()
{
    RDouble reynoldsReferenceLengthDimensional = 1.0;

    RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

    reynoldsReferenceLengthDimensional = gridScaleFactor;

    GlobalDataBase::UpdateData("reynoldsReferenceLengthDimensional", &reynoldsReferenceLengthDimensional, PHDOUBLE, 1);
}





}