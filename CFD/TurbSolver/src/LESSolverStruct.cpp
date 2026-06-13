#include "PHMpi.h"
#include "LESSolverStruct.h"
#include "GradientOperation.h"
#include "Math_BasisFunction.h"
#include "Constants.h"
#include "Param_LESSolverStruct.h"
#include "TK_Exit.h"
#include "Glb_Dimension.h"
#include "AleModel.h"
#include "NSSolverStructFD.h"
#include "IO_FileName.h"
#pragma warning(disable:6385)
#pragma warning(disable:26451)
using namespace std;

namespace PHSPACE
{

LESSolverStruct::LESSolverStruct()
{

}

LESSolverStruct::~LESSolverStruct()
{
    DeAllocateGlobalVariables();
    FreeControlParameters();
}

void LESSolverStruct::AllocateGlobalVar(Grid *gridIn)
{
    Param_LESSolverStruct *parameters = GetControlParameters();

    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Range IFACE, JFACE, KFACE;
    GetRange(ni, nj, nk, 0, 0, IFACE, JFACE, KFACE);

    Range Direction(1, 3);
    Range N(1, 6);
    Range M(1, 3);

    RDouble3D *cellLengthScale = new RDouble3D(Range(1, ni-1), Range(1, nj-1), Range(1, nk-1), fortranArray);    //! Length scale of cell.
    *cellLengthScale = 0.0;
    grid->UpdateDataPtr("cellLengthScale", cellLengthScale);    //! Length scale of cell.

    RDouble3D *wallFunction = new RDouble3D(I, J, K, fortranArray); //! Wall damping function.
    *wallFunction = 1.0;
    grid->UpdateDataPtr("wallFunction", wallFunction);    //! Wall damping function.

    RDouble3D *anisotropicConstant = new RDouble3D(I, J, K, fortranArray);    //! Model coefficients of the anisotropic SGS stresses using in dynamic model.
    *anisotropicConstant = 0.0;
    grid->UpdateDataPtr("anisotropicConstant", anisotropicConstant);    //! Model coefficients of the anisotropic SGS stresses using in dynamic model.

    RDouble3D *isotropicConstant = new RDouble3D(I, J, K, fortranArray);    //! Model coefficients of the isotropic SGS stresses using in dynamic model.
    *isotropicConstant = 0.0;
    grid->UpdateDataPtr("isotropicConstant", isotropicConstant);    //! Model coefficients of the isotropic SGS stresses using in dynamic model.

    RDouble3D *turbulentPrandtlNumber = new RDouble3D(I, J, K, fortranArray);    //! Turbulence Prandtl number using in dynamic model.
    *turbulentPrandtlNumber = 0.9;
    grid->UpdateDataPtr("turbulentPrandtlNumber", turbulentPrandtlNumber);    //! Turbulence Prandtl number using in dynamic model.

    RDouble4D *strainRateTensor = new RDouble4D(I, J, K, Range(0, 9), fortranArray);    //! Strain Rate Tensor Sij and Omegaij
    *strainRateTensor = 0.0;
    grid->UpdateDataPtr("strainRateTensor", strainRateTensor);    //! Strain Rate Tensor Sij and Omegaij.

    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (isFVMOrFDM == FD_METHOD)
    {
        Range I4(- 3, ni + 3);
        Range J4(- 3, nj + 3);
        Range K4(- 3, nk + 3);
        UVWTproxy = new FieldProxy();    //! Flow field variables u/v/w/t used in calculating gradient of u/v/w/t;
        UVWTfaceIproxy = new FieldProxy();    //! Flow field variables u/v/w/t defined on I faces.
        UVWTfaceJproxy = new FieldProxy();    //! Flow field variables u/v/w/t defined on J faces.
        UVWTfaceKproxy = new FieldProxy();    //! Flow field variables u/v/w/t defined on K faces.
        UVWTproxy->SetField_STR(new RDouble4D(I4, J4, K4, Range(0, 3), fortranArray), true);
        UVWTfaceIproxy->SetField_STR(new RDouble4D(I, J, K, Range(0, 3), fortranArray), true); 
        UVWTfaceJproxy->SetField_STR(new RDouble4D(I, J, K, Range(0, 3), fortranArray), true);
        UVWTfaceKproxy->SetField_STR(new RDouble4D(I, J, K, Range(0, 3), fortranArray), true);
    }

    int subgridScaleModel = parameters->GetSubgridScaleModel();
    if (subgridScaleModel == 2 || subgridScaleModel == 4)
    {
        RDouble4D *leonardStress = new RDouble4D(I, J, K, N, fortranArray);    //! Leonard stresses using in dynamic model.
        RDouble4D *modeledStress = new RDouble4D(I, J, K, N, fortranArray);    //! Modeled stresses using in dynamic model.
        RDouble4D *leonardTemperatureFlux = new RDouble4D(I, J, K, Range(0, 2), fortranArray);    //! Leonard heat flux using in dynamic model.
        RDouble4D *modeledTemperatureFlux = new RDouble4D(I, J, K, Range(0, 2), fortranArray);    //! Modeled heat flux using in dynamic model.
        RDouble4D *testFilteredQ = new RDouble4D(I, J, K, Range(0, 3), fortranArray);    //! Flow field variables after test filtering.
        RDouble3D *testFilteredTemperature = new RDouble3D(I, J, K, fortranArray);    //! Temperature after test filtering.
        RDouble4D *testFilteredStrainRateTensor = new RDouble4D(I, J, K, Range(0, 6), fortranArray);    //! Strain rate tensor after test filtering.
        RDouble4D *testFilteredRhoUiUj = new RDouble4D(I, J, K, N, fortranArray);    //! rhoUiUj after test filtering.
        RDouble4D *testFilteredRhoUT = new RDouble4D(I, J, K, M, fortranArray);    //! rhoUT after test filtering.
        RDouble4D *testFilteredAlphaIJ = new RDouble4D(I, J, K, Range(0, 6), fortranArray);
        RDouble4D *betaIJ = new RDouble4D(I, J, K, N, fortranArray);
        RDouble4D *gradT = new RDouble4D(I, J, K, Range(0, 2), fortranArray);    //! Temperature gradient.
        RDouble4D *testFilteredGradT = new RDouble4D(I, J, K, Range(0, 2), fortranArray);    //! Temperature gradient after test filtering.
        RDouble4D *testFilteredRhoStrainRateMagnitudeGradT = new RDouble4D(I, J, K, Range(0, 2), fortranArray);    //! rho|S|gradT after test filtering.

        *leonardStress = 0.0;
        *modeledStress= 0.0;
        *leonardTemperatureFlux = 0.0;
        *modeledTemperatureFlux = 0.0;
        *testFilteredQ = 0.0;
        *testFilteredTemperature = 0.0;
        *testFilteredStrainRateTensor = 0.0;
        *testFilteredRhoUiUj = 0.0;
        *testFilteredRhoUT = 0.0;
        *testFilteredAlphaIJ = 0.0;
        *betaIJ = 0.0;
        *gradT = 0.0;
        *testFilteredGradT = 0.0;
        *testFilteredRhoStrainRateMagnitudeGradT = 0.0;
 
        grid->UpdateDataPtr("leonardStress", leonardStress);    //! Leonard stresses using in dynamic model.
        grid->UpdateDataPtr("modeledStress", modeledStress);    //! Modeled stresses using in dynamic model.
        grid->UpdateDataPtr("leonardTemperatureFlux", leonardTemperatureFlux);    //! Leonard heat flux using in dynamic model.
        grid->UpdateDataPtr("modeledTemperatureFlux", modeledTemperatureFlux);    //! Modeled heat flux using in dynamic model.
        grid->UpdateDataPtr("testFilteredQ", testFilteredQ);    //! Flow field variables after test filtering.
        grid->UpdateDataPtr("testFilteredTemperature", testFilteredTemperature);    //! Temperature after test filtering.
        grid->UpdateDataPtr("testFilteredStrainRateTensor", testFilteredStrainRateTensor);    //! Strain rate tensor after test filtering.
        grid->UpdateDataPtr("testFilteredRhoUiUj", testFilteredRhoUiUj);    //! rhoUiUj after test filtering.
        grid->UpdateDataPtr("testFilteredRhoUT", testFilteredRhoUT);    //! rhoUT after test filtering.
        grid->UpdateDataPtr("testFilteredAlphaIJ", testFilteredAlphaIJ);
        grid->UpdateDataPtr("betaIJ", betaIJ);
        grid->UpdateDataPtr("gradT", gradT);    //! Temperature gradient.
        grid->UpdateDataPtr("testFilteredGradT", testFilteredGradT);    //! Temperature gradient after test filtering.
        grid->UpdateDataPtr("testFilteredRhoStrainRateMagnitudeGradT", testFilteredRhoStrainRateMagnitudeGradT);    //! rho|S|gradT after test filtering.
    }

    Init(grid);
}

void LESSolverStruct::DeAllocateGlobalVar(Grid *gridIn)
{
    Param_LESSolverStruct *parameters = GetControlParameters();
    StructGrid *grid = StructGridCast(gridIn);

    RDouble3D *cellLengthScale = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("cellLengthScale"));
    RDouble3D *wallFunction = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));
    delete cellLengthScale;
    delete wallFunction;

    RDouble3D *anisotropicConstant = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));
    RDouble3D *isotropicConstant = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("isotropicConstant"));
    RDouble3D *turbulentPrandtlNumber = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));
    RDouble4D *strainRateTensor = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));

    delete isotropicConstant;
    delete anisotropicConstant;
    delete turbulentPrandtlNumber;
    delete strainRateTensor;

    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (isFVMOrFDM == FD_METHOD)
    {
        delete UVWTproxy;
        delete UVWTfaceIproxy;
        delete UVWTfaceJproxy;
        delete UVWTfaceKproxy;
    }

    int subgridScaleModel = parameters->GetSubgridScaleModel();
    if (subgridScaleModel == 2 || subgridScaleModel == 4)
    { 
        RDouble4D *leonardStress = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("leonardStress"));
        RDouble4D *modeledStress = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("modeledStress"));
        RDouble4D *leonardTemperatureFlux = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("leonardTemperatureFlux"));
        RDouble4D *modeledTemperatureFlux = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("modeledTemperatureFlux"));
        RDouble4D *testFilteredQ = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));
        RDouble3D *testFilteredTemperature = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("testFilteredTemperature"));
        RDouble4D *testFilteredStrainRateTensor = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredStrainRateTensor"));
        RDouble4D *testFilteredRhoUiUj = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoUiUj"));
        RDouble4D *testFilteredRhoUT = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoUT"));
        RDouble4D *testFilteredAlphaIJ = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredAlphaIJ"));
        RDouble4D *betaIJ = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("betaIJ"));
        RDouble4D *gradT = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradT"));
        RDouble4D *testFilteredGradT = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredGradT"));
        RDouble4D *testFilteredRhoStrainRateMagnitudeGradT = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoStrainRateMagnitudeGradT"));

        delete leonardStress;
        delete modeledStress;
        delete leonardTemperatureFlux;
        delete modeledTemperatureFlux;
        delete testFilteredQ;
        delete testFilteredTemperature;
        delete testFilteredStrainRateTensor;
        delete testFilteredRhoUiUj;
        delete testFilteredRhoUT;
        delete testFilteredAlphaIJ;
        delete betaIJ;
        delete gradT;
        delete testFilteredGradT;
        delete testFilteredRhoStrainRateMagnitudeGradT;
    }
}

LIB_EXPORT void LESSolverStruct::InitControlParameters()
{
    FreeControlParameters();
    controlParameters = new Param_LESSolverStruct();
    controlParameters->Init();
}

LIB_EXPORT Param_LESSolverStruct * LESSolverStruct::GetControlParameters()
{
    return static_cast<Param_LESSolverStruct *> (controlParameters);
}

bool LESSolverStruct::JudgeIfRestart()
{
    string turbFile = ".\results\turb.dat";
    GlobalDataBase::GetData("turbfile", &turbFile, PHSTRING, 1);

    if (PHMPI::IsParallelRun())
    {
        turbFile = PHSPACE::AddSymbolToFileName(turbFile, "_", 0);
    }

    bool restart_flag = false;

    ifstream infile(turbFile.c_str(), ios::in);
    if (infile)
    {
        restart_flag = true;
        infile.close();
        infile.clear();
    }

    return restart_flag;
}

bool LESSolverStruct::JudgeIfReadAverage()
{
    string turbFile = ".\results\turb.dat";
    GlobalDataBase::GetData("turbfile", &turbFile, PHSTRING, 1);
    string averageTurbFile = AddSymbolToFileName(turbFile, "_", "Average");

    if (PHMPI::IsParallelRun())
    {
        averageTurbFile = PHSPACE::AddSymbolToFileName(averageTurbFile, "_", 0);
    }

    bool restart_flag = false;

    ifstream infile(averageTurbFile.c_str(), ios::in);
    if (infile)
    {
        restart_flag = true;
        infile.close();
        infile.clear();
    }

    return restart_flag;
}

void LESSolverStruct::DumpRestartH5(ActionKey *actkey)
{
    using namespace PHMPI;
    int currentProcessorID = GetCurrentProcessorID();

    int version = 1;
    int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (currentProcessorID == GetServerProcessorID())
    {
        WriteData(actkey->filepos, &version, "Version");
        WriteData(actkey->filepos, &outIterStep, "outnstep");

        if (IsNeedStatistics())
        {
            int nStatisticalStep = 0;
            GlobalDataBase::GetData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
            WriteData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
        }
    }

    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    StructGrid *strGrid = StructGridCast(GetGrid(actkey->level));
    int nTotalCell = strGrid->GetNTotalCell();

    int ist, ied, jst, jed, kst, ked;
    strGrid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss << "Group" << gridID;
    grpName = oss.str();

    grploc = OpenGroup(actkey->filepos, grpName);

    WriteData(grploc, &nTotalCell, "nTotalCell");
    WriteData(grploc, &viscousTurbulent(ist, jst, kst), "visturb");

    H5Gclose(grploc);
}

void LESSolverStruct::ReadRestartH5(ActionKey *actkey)
{
    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    StructGrid *strGrid = StructGridCast(GetGrid(actkey->level));
    int nTotalCell = strGrid->GetNTotalCell();

    int ist, ied, jst, jed, kst, ked;
    strGrid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    int outIterstepofNS = GlobalDataBase::GetIntParaFromDB("outnstep");
    int outIterStepofTurb = 0;
    ReadData(actkey->filepos, &outIterStepofTurb, "outnstep");

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss << "Group" << gridID;
    grpName = oss.str();

    grploc = OpenGroup(actkey->filepos, grpName);

    int nTotalCellRestart = 0;
    ReadData(grploc, &nTotalCellRestart, "nTotalCell");
    if (nTotalCellRestart != nTotalCell)
    {
        ostringstream erroeInfo;
        erroeInfo << " Error: the cell number in turb.dat is not equal to the cell number in grid file !" << endl;
        TK_Exit::ExceptionExit(erroeInfo.str());
    }

    ReadData(grploc, &viscousTurbulent(ist, jst, kst), "visturb");

    H5Gclose(grploc);

    CompareOutStepOfFlowfieldFile(outIterstepofNS, outIterStepofTurb);
}

void LESSolverStruct::InitFlowAsRestart()
{
    Param_LESSolverStruct *parameters = GetControlParameters();

    StructGrid *grid = StructGridCast(GetGrid(0));

    RDouble3D &viscousTurbulence = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble3D &subgridScaleEnergy = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("subgridScaleEnergy"));
    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));

    RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                viscousTurbulence(i, j, k) = freeStreamViscosity;

                subgridScaleEnergy(i, j, k) = 0.0;

                turbulentPrandtlNumber(i, j, k) = 0.9;
            }
        }
    }
}

void LESSolverStruct::Init(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    ComputeCellLengthScale(grid);
}

void LESSolverStruct::UploadInterfaceData(ActionKey *actkey)
{
    UploadInterfaceValue(actkey);
}

void LESSolverStruct::UploadInterfaceValue(ActionKey *actkey)
{
    int level = actkey->level;
    StructGrid *grid = reinterpret_cast<StructGrid *> (GetGrid(level));
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo)
    {
        return;
    }

    RDouble3D *viscousTurbulence = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    PHSPACE::UploadInterfaceValue(grid, viscousTurbulence, "vist");

    RDouble3D *subgridScaleEnergy = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("subgridScaleEnergy"));
    PHSPACE::UploadInterfaceValue(grid, subgridScaleEnergy, "subgridScaleEnergy");

    RDouble3D *turbulentPrandtlNumber = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));
    PHSPACE::UploadInterfaceValue(grid, turbulentPrandtlNumber, "turbulentPrandtlNumber");
}

void LESSolverStruct::DownloadInterfaceData(ActionKey *actkey)
{
    DownloadInterfaceValue(actkey);
}

void LESSolverStruct::DownloadInterfaceValue(ActionKey *actkey)
{
    int level = actkey->level;
    StructGrid *grid = reinterpret_cast<StructGrid *> (GetGrid(level));
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo)
    {
        return;
    }

    RDouble3D *viscousTurbulence = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    PHSPACE::DownloadInterfaceValue(grid, viscousTurbulence, "vist");

    RDouble3D *subgridScaleEnergy = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("subgridScaleEnergy"));
    PHSPACE::DownloadInterfaceValue(grid, subgridScaleEnergy, "subgridScaleEnergy");

    RDouble3D *turbulentPrandtlNumber = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));
    PHSPACE::DownloadInterfaceValue(grid, turbulentPrandtlNumber, "turbulentPrandtlNumber");
}

void LESSolverStruct::RotateVectorFromInterface(Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
{
    StructGrid *grid = StructGridCast(gridIn);
    int level = grid->GetLevel();
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();
    StructBCSet * structBCSet = finestGrid->GetStructBCSet();
    RDouble4D &rotNSgradValueX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotNSgradValueX"));
    RDouble4D &rotNSgradValueY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotNSgradValueY"));
    RDouble4D &rotNSgradValueZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotNSgradValueZ"));

    rotNSgradValueX = 0.0;
    rotNSgradValueY = 0.0;
    rotNSgradValueZ = 0.0;

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);
    RDouble rotationAngle = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("rotationAngle");
    rotationAngle = rotationAngle*PI/180;
    if (nEquation > 0)
    {
        RDouble4D &fieldRecvNSY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
        RDouble4D &fieldRecvNSZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                int it, jt, kt;
                finestGrid->GetTargetIndexIJK(iFace, iGhostLayer + 1, it, jt, kt);
                finestGrid->RemapMultigridIJK(level, it, jt, kt);
                int *ibcregions = structBCSet->GetIFaceInfo();
                int iBCRegion = ibcregions[iFace];
                StructBC * bcregion = structBCSet->GetBCRegion(iBCRegion);

                string bcName = bcregion->GetBCName();
                if (bcName == "Periodic_up")
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        rotNSgradValueY(it, jt, kt, m) = fieldRecvNSY(it, jt, kt, m) * cos(2*PI-rotationAngle) - fieldRecvNSZ(it, jt, kt, m) * sin(2*PI-rotationAngle);
                        rotNSgradValueZ(it, jt, kt, m) = fieldRecvNSY(it, jt, kt, m) * sin(2*PI-rotationAngle) + fieldRecvNSZ(it, jt, kt, m) * cos(2*PI-rotationAngle);
                    }
                }
                else if (bcName == "Periodic_down")
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        rotNSgradValueY(it, jt, kt, m) = fieldRecvNSY(it, jt, kt, m) * cos(rotationAngle) - fieldRecvNSZ(it, jt, kt, m) * sin(rotationAngle);
                        rotNSgradValueZ(it, jt, kt, m) = fieldRecvNSY(it, jt, kt, m) * sin(rotationAngle) + fieldRecvNSZ(it, jt, kt, m) * cos(rotationAngle);
                    }
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    fieldRecvNSY(it, jt, kt, m) = rotNSgradValueY(it, jt, kt, m);
                    fieldRecvNSZ(it, jt, kt, m) = rotNSgradValueZ(it, jt, kt, m);
                }
            }
        }
    }
}

/*
void LESSolverStruct::ComputeCellLengthScale(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble3D &cellLengthScale = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("cellLengthScale"));

    RDouble4D &area = *(grid->GetFaceArea());
    RDouble3D &vol  = *(grid->GetCellVolume());

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int nDim = GetDim();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);
    //double areaL[3], areaR[3], delta[3];

    int i, j, k;

    //! inner cells
    for (k = kst; k <= ked; ++ k)
    {
        for (j = jst; j <= jed; ++ j)
        {
            for (i = ist; i <= ied; ++ i)
            {
                RDouble len3 = 1.0;

                for (int iSurface = 1; iSurface <= nDim; ++iSurface)
                {
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);

                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    //areaL[iSurface - 1] = area(i, j, k, iSurface);
                    //areaR[iSurface - 1] = area(il, jl, kl, iSurface);
                    RDouble areaL = area(i, j, k, iSurface);
                    RDouble areaR = area(il, jl, kl, iSurface);

                    RDouble volume = vol(i, j, k);

                    //delta[iSurface - 1] = 2.0 * volume / (areaL[iSurface - 1] + areaR[iSurface - 1]);
                    RDouble delta = 2.0 * volume / (areaL + areaR);

                    len3 = len3 * delta;
                }
                
                cellLengthScale(i, j, k) = pow(len3, third);
            }
        }
    }

    //! ghost cells, for i = 0 and ni;
    i = 0;

    for (k = kst; k <= ked; ++ k)
    {
        for (j = jst; j <= jed; ++ j)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }

    i = ni;

    for (k = kst; k <= ked; ++ k)
    {
        for (j = jst; j <= jed; ++ j)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }

    //! ghost cells, for j = 0 and nj;
    j = 0;

    for (k = kst; k <= ked; ++ k)
    {
        for (i = ist; i <= ied; ++ i)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }

    j = nj;

    for (k = kst; k <= ked; ++ k)
    {
        for (i = ist; i <= ied; ++ i)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }

    //! ghost cells, for k = 0 and nk;
    k = 0;

    for (j = jst; j <= jed; ++ j)
    {
        for (i = ist; i <= ied; ++ i)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }

    k = nk;

    for (j = jst; j <= jed; ++ j)
    {
        for (i = ist; i <= ied; ++ i)
        {
            RDouble len3 = 1.0;

            for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
            {
                int il1, jl1, kl1;
                GetNsurfIndex(iSurface, il1, jl1, kl1);

                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble areaL = area(i, j, k, iSurface);
                RDouble areaR = area(il, jl, kl, iSurface);

                RDouble volume = vol(i, j, k);

                RDouble delta = 2.0 * volume / (areaL + areaR);

                len3 = len3 * delta;
            }

            cellLengthScale(i, j, k) = pow(len3, third);
        }
    }
}
*/

void LESSolverStruct::ComputeCellLengthScale(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble3D &cellLengthScale = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("cellLengthScale"));
/*
    RDouble4D &area = *(grid->GetFaceArea());
    RDouble3D &vol  = *(grid->GetCellVolume());

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int nDim = GetDim();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);
    //double areaL[3], areaR[3], delta[3];

    int i, j, k;

    //! inner cells
    for (k = kst; k <= ked; ++ k)
    {
        for (j = jst; j <= jed; ++ j)
        {
            for (i = ist; i <= ied; ++ i)
            {
                RDouble len3 = 1.0;

                for (int iSurface = 1; iSurface <= nDim; ++iSurface)
                {
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);

                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    //areaL[iSurface - 1] = area(i, j, k, iSurface);
                    //areaR[iSurface - 1] = area(il, jl, kl, iSurface);
                    RDouble areaL = area(i, j, k, iSurface);
                    RDouble areaR = area(il, jl, kl, iSurface);

                    RDouble volume = vol(i, j, k);

                    //delta[iSurface - 1] = 2.0 * volume / (areaL[iSurface - 1] + areaR[iSurface - 1]);
                    RDouble delta = 2.0 * volume / (areaL + areaR);

                    len3 = len3 * delta;
                }
                
                cellLengthScale(i, j, k) = pow(len3, third);
            }
        }
    }

    GhostCell3D(cellLengthScale, ni, nj, nk);
*/
    RDouble3D *largestLocalGridLength = grid->GetLargestLocalGridLength();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                cellLengthScale(i, j, k) = (*largestLocalGridLength)(i, j, k);
            }
        }
    }
}

void LESSolverStruct::ComputeWallFunction(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    
    Param_LESSolverStruct *parameters = GetControlParameters();
    int wallDampingFunctionType = parameters->GetWallDampingFunctionType();

    RDouble3D &wallFunction = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));

    int iterationStep;
    int nMGLevel = parameters->GetNMGLevel();        
    if (nMGLevel <= 1)
    {
        //! Finest grid.
        GlobalDataBase::GetData("outnstep", &iterationStep, PHINT, 1);
    }
    else
    {
        //! Coarse grid.
        GlobalDataBase::GetData("newnstep", &iterationStep, PHINT, 1);
    }

    int cflNstep = parameters->GetCFLVaryStep();
    if (iterationStep < cflNstep)
    {
        WallFunctionofVanDriest(grid);

        return;
    }

    if (wallDampingFunctionType == 0)
    {
        wallFunction = 1.0;
    }
    else if (wallDampingFunctionType == 1)
    {
        WallFunctionofVanDriest(grid);
    }
    else if (wallDampingFunctionType == 2)
    {
        WallFunctionofDXB(grid);
    }
    else if (wallDampingFunctionType == 3)
    {
        WallFunctionofPiomelli(grid);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("wallDampingFunctionType", wallDampingFunctionType);
    }
}

void LESSolverStruct::WallFunctionofVanDriest(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    RDouble3D & wallDistance = * StructGridCast(grid)->GetWallDist();

    RDouble4D & q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble3D & visl = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble3D &wallFunction = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("wallFunction"));

    RDouble4D &strainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));

    Param_LESSolverStruct *parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();

    const double kappa = 0.41;
    const double Ap = 25.0;

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble wallDist = wallDistance(i, j, k);

                RDouble mul = visl(i, j, k);

                RDouble rho = q(i, j, k, IR);

                RDouble strain = strainRateTensor(i, j, k, 0);

                RDouble yplus = strain * (kappa * wallDist) * (kappa * wallDist) * rho / mul * refReNumber;

                RDouble wf = 1.0 - exp(-yplus / Ap);
                wallFunction(i, j, k) = wf;
            }
        }
    }
}

void LESSolverStruct::WallFunctionofDXB(Grid *gridIn)
{
    using namespace IDX;

    StructGrid * grid = StructGridCast(gridIn);

    RDouble3D & wallDistance = * StructGridCast(grid)->GetWallDist();

    RDouble4D & q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble3D & visl = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble3D &wallFunction = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));

    RDouble4D &strainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));

    Param_LESSolverStruct *parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble rgam = parameters->GetRgam();

    const double kappa = 0.41;
    const double Ap = 25.0;
    const double yploge = 1000.0;

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble wallDist = wallDistance(i, j, k);

                RDouble mul = visl(i, j, k);

                RDouble rho = q(i, j, k, IR);

                RDouble strain = strainRateTensor(i, j, k, 0);

                RDouble yplus = strain * (kappa * wallDist) * (kappa * wallDist) * rho / mul * refReNumber;

                RDouble ffd = tanh((yploge/yplus) * (yploge/yplus) * (yploge/yplus));
                RDouble vanDriestInner = rgam * MIN((1.0 - exp(-pow(yplus/Ap, 7.0/6.0))), (1.0 - exp(-yplus / Ap)));
                RDouble vanDriestOut = 1.0;

                wallFunction(i, j, k) = vanDriestInner * ffd + vanDriestOut * (1.0 - ffd);
            }
        }
    }
}

void LESSolverStruct::WallFunctionofPiomelli(Grid *gridIn)
{
    using namespace IDX;

    StructGrid * grid = StructGridCast(gridIn);

    RDouble3D & wallDistance = * StructGridCast(grid)->GetWallDist();

    RDouble4D & q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble3D & visl = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble3D &wallFunction = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));

    RDouble4D &strainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));

    Param_LESSolverStruct *parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();

    const double kappa = 0.41;
    const double Ap = 25.0;

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble wallDist = wallDistance(i, j, k);

                RDouble mul = visl(i, j, k);

                RDouble rho = q(i, j, k, IR);

                RDouble strain = strainRateTensor(i, j, k, 0);

                RDouble yplus = strain * (kappa * wallDist) * (kappa * wallDist) * rho / mul * refReNumber;

                wallFunction(i, j, k) = sqrt(1.0 - exp(-(yplus/Ap) * (yplus/Ap) * (yplus/Ap)));
            }
        }
    }
}

void LESSolverStruct::ComputeViscousCoeff(Grid *gridIn)
{
}

void LESSolverStruct::GetDependentVariablesforStructHighOrder(Grid *gridIn)
{
    ObtainGradientCellCenterHighOrder(gridIn);
}

void LESSolverStruct::ComputeGradientCellCenter(Grid *gridIn)
{
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (isFVMOrFDM == FD_METHOD) return;

    StructGrid *strgrid = StructGridCast(gridIn);
    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    RDouble4D &xfv = *(strgrid->GetFaceVectorX());
    RDouble4D &yfv = *(strgrid->GetFaceVectorY());
    RDouble4D &zfv = *(strgrid->GetFaceVectorZ());
    RDouble3D &volume      = *(strgrid->GetCellVolume());

    int ndim = GetDim();
    int nTemperatureModel = 1;

    RDouble4D &q = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("q"));
    RDouble4D &temperature = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("t"));

    RDouble4D &dqdx = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &dqdy = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &dqdz = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterZ"));

    dqdx = 0.0;
    dqdy = 0.0;
    dqdz = 0.0;

    for (int iSurface = 1; iSurface <= ndim; ++ iSurface)
    {
        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);

        int ist = 1;
        int ied = ni - 1 + il1;
        int jst = 1;
        int jed = nj - 1 + jl1;
        int kst = 1;
        int ked = nk - 1 + kl1;

        if (ndim == TWO_D)
        {
            ked = 1;
        }
        int il, jl, kl;

        //! The gradient of u,v,w.
        for (int m = 0; m < ndim; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        il = i - il1;
                        jl = j - jl1;
                        kl = k - kl1;

                        RDouble phis = q(i,j,k,m+1) + q(il,jl,kl,m+1);

                        RDouble ddx = phis * xfv(i, j, k, iSurface);
                        RDouble ddy = phis * yfv(i, j, k, iSurface);
                        RDouble ddz = phis * zfv(i, j, k, iSurface);

                        dqdx(i ,j ,k ,m) -= ddx;
                        dqdy(i ,j ,k ,m) -= ddy;
                        dqdz(i ,j ,k ,m) -= ddz;

                        dqdx(il, jl, kl, m) += ddx;
                        dqdy(il, jl, kl, m) += ddy;
                        dqdz(il, jl, kl, m) += ddz;
                    }
                }
            }
        }

        //! The gradient of temperature.
        for (int m = 0; m < nTemperatureModel; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        il = i - il1;
                        jl = j - jl1;
                        kl = k - kl1;

                        RDouble phis = temperature(i, j, k, m) + temperature(il, jl, kl, m);

                        RDouble ddx = phis * xfv(i, j, k, iSurface);
                        RDouble ddy = phis * yfv(i, j, k, iSurface);
                        RDouble ddz = phis * zfv(i, j, k, iSurface);

                        dqdx(i, j, k, m + 3) -= ddx;
                        dqdy(i, j, k, m + 3) -= ddy;
                        dqdz(i, j, k, m + 3) -= ddz;

                        dqdx(il, jl, kl, m + 3) += ddx;
                        dqdy(il, jl, kl, m + 3) += ddy;
                        dqdz(il, jl, kl, m + 3) += ddz;
                    }
                }
            }
        }
    }

    int ist = 1;
    int ied = ni-1;
    int jst = 1;
    int jed = nj-1;
    int kst = 1;
    int ked = nk-1;

    if (ndim == TWO_D) ked = 1;

    int nTotalVariable = nTemperatureModel + 3;
    for (int m = 0; m < nTotalVariable; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble oov = half / volume(i, j, k);
                    dqdx(i, j, k, m) *= oov;
                    dqdy(i, j, k, m) *= oov;
                    dqdz(i, j, k, m) *= oov;
                }
            }
        }
    }
}

void LESSolverStruct::ObtainGradientCellCenterHighOrder(Grid *gridIn)
{
    StructGrid *strgrid = StructGridCast(gridIn);

    GetUVWTproxy(strgrid);

    int nDim = GetDim();
    for (int nsurf = 1; nsurf <= nDim; ++nsurf)
    {
        GetUVWTfaceproxy(strgrid, nsurf);
        GetdqdkxietactaCellCenter(strgrid, nsurf);
    }

    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    RDouble4D &xfv  = *(strgrid->GetFaceVectorX());
    RDouble4D &yfv  = *(strgrid->GetFaceVectorY());
    RDouble4D &zfv  = *(strgrid->GetFaceVectorZ());
    RDouble3D &vol  = *(strgrid->GetCellVolume());  

    RDouble4D &dqdx = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &dqdy = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &dqdz = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterZ"));  

    RDouble4D &dqdkxi = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdkxi"));
    RDouble4D &dqdeta = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdeta"));
    RDouble4D &dqdcta = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdcta"));

    int ist = 1;
    int jst = 1;
    int kst = 1;
    int ied = ni - 1;
    int jed = nj - 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    } 

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble kxix = half * (xfv(i, j, k, 1) + xfv(i + 1, j, k, 1)); 
                RDouble kxiy = half * (yfv(i, j, k, 1) + yfv(i + 1, j, k, 1));
                RDouble kxiz = half * (zfv(i, j, k, 1) + zfv(i + 1, j, k, 1));

                RDouble etax = half * (xfv(i, j, k, 2) + xfv(i, j + 1, k, 2)); 
                RDouble etay = half * (yfv(i, j, k, 2) + yfv(i, j + 1, k, 2));
                RDouble etaz = half * (zfv(i, j, k, 2) + zfv(i, j + 1, k, 2));
                for (int m = 0; m <= 3; ++ m)
                {
                    dqdx(i, j, k, m) = kxix * dqdkxi(i, j, k, m) + etax * dqdeta(i, j, k, m);
                    dqdy(i, j, k, m) = kxiy * dqdkxi(i, j, k, m) + etay * dqdeta(i, j, k, m);
                    dqdz(i, j, k, m) = kxiz * dqdkxi(i, j, k, m) + etaz * dqdeta(i, j, k, m);
                }
            }
        }
    }

    if (nDim == 3)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble ctax = half * (xfv(i, j, k, 3) + xfv(i, j, k + 1, 3));
                    RDouble ctay = half * (yfv(i, j, k, 3) + yfv(i, j, k + 1, 3));
                    RDouble ctaz = half * (zfv(i, j, k, 3) + zfv(i, j, k + 1, 3));
                    for (int m = 0; m <= 3; ++ m)
                    {
                        dqdx(i, j, k, m) += ctax * dqdcta(i, j, k, m);
                        dqdy(i, j, k, m) += ctay * dqdcta(i, j, k, m);
                        dqdz(i, j, k, m) += ctaz * dqdcta(i, j, k, m);
                    }
                }
            }
        }
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m <= 3; ++ m)
                {
                    dqdx(i, j, k, m) /= vol(i, j, k);
                    dqdy(i, j, k, m) /= vol(i, j, k);
                    dqdz(i, j, k, m) /= vol(i, j, k);
                }
            }
        }
    }

    GhostCell3D(dqdx, ni, nj, nk, 4);
    GhostCell3D(dqdy, ni, nj, nk, 4);
    GhostCell3D(dqdz, ni, nj, nk, 4);
}

void LESSolverStruct::GetUVWTproxy(Grid *gridIn)
{
    StructGrid *strgrid = StructGridCast(gridIn);
    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    RDouble4D &q = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("q_FD"));
    RDouble4D &t = *reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("t_FD"));
    RDouble4D &UVWT = UVWTproxy->GetField_STR();

    int ist = - 3;
    int jst = - 3;
    int kst = - 3;
    int ied = ni + 3;
    int jed = nj + 3;
    int ked = nk + 3;
    
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                UVWT(i, j, k, 0) = q(i, j, k, 1);
                UVWT(i, j, k, 1) = q(i, j, k, 2);
                UVWT(i, j, k, 2) = q(i, j, k, 3);
                UVWT(i, j, k, 3) = t(i, j, k, 0);
            }
        }
    }
}

void LESSolverStruct::GetUVWTfaceproxy(Grid *gridIn, int nsurf)
{
    string highordersolvername = GlobalDataBase::GetStrParaFromDB("str_highorder_solver");
    if (highordersolvername.substr(0, 4) == "HDCS")
    {
        GetUVWTfaceproxy_HDCS(gridIn, nsurf);
        return;
    }

    StructGrid *strgrid = StructGridCast(gridIn);
    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    RDouble4D &UVWT = UVWTproxy->GetField_STR();

    RDouble4D *UVWTface;
    Int1D *DiffOrder;
    if (nsurf == 1)
    {
        UVWTface = &(UVWTfaceIproxy->GetField_STR());
        DiffOrder = strgrid->DiscretePrecisionKXI;
    }
    else if (nsurf == 2)
    {
        UVWTface = &(UVWTfaceJproxy->GetField_STR());
        DiffOrder = strgrid->DiscretePrecisionETA;
    }
    else
    {
        UVWTface = &(UVWTfaceKproxy->GetField_STR());
        DiffOrder = strgrid->DiscretePrecisionCTA;
    }
    
    int il1, jl1, kl1;
    GetNsurfIndex(nsurf, il1, jl1, kl1);

    int ist = 1 - il1;
    int jst = 1 - jl1;
    int kst = 1 - kl1;
    int ied = ni - 1 + 2*il1;
    int jed = nj - 1 + 2*jl1;
    int ked = nk - 1 + 2*kl1;    
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int localcellID = i * il1 + j * jl1 + k * kl1;
                
                if ((*DiffOrder)(localcellID) == 2 || (*DiffOrder)(localcellID - 1) == 2)
                {
                    for (int m = 0; m <= 3; ++m)
                    {                    
                        RDouble b = UVWT(i - il1,   j - jl1,   k - kl1,   m);
                        RDouble c = UVWT(i,         j,         k,         m);
                        (*UVWTface)(i, j, k, m) = half * (b + c);
                    }
                }
                else
                {
                    for (int m = 0; m <= 3; ++m)
                    {                    
                        RDouble a = UVWT(i - il1*2, j - jl1*2, k - kl1*2, m);
                        RDouble b = UVWT(i - il1,   j - jl1,   k - kl1,   m);
                        RDouble c = UVWT(i,         j,         k,         m);
                        RDouble d = UVWT(i + il1,   j + jl1,   k + kl1,   m);
                        (*UVWTface)(i, j, k, m) = (- a + 9.0 * b + 9.0 * c - d) / 16.0;
                    }
                    if ((*UVWTface)(i, j, k, 3) < SMALL)
                    {
                        (*UVWTface)(i, j, k, 3) = half * (UVWT(i - il1, j - jl1, k - kl1,  3) + UVWT(i, j, k, 3));
                    }
                }                
            }
        }
    }   
}

void LESSolverStruct::GetUVWTfaceproxy_HDCS(Grid *gridIn, int nsurf)
{
    StructGrid *strgrid = StructGridCast(gridIn);
    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    RDouble4D &UVWT = UVWTproxy->GetField_STR();

    RDouble4D *UVWTface;
    if (nsurf == 1)
    {
        UVWTface = &(UVWTfaceIproxy->GetField_STR());
    }
    else if (nsurf == 2)
    {
        UVWTface = &(UVWTfaceJproxy->GetField_STR());
    }
    else
    {
        UVWTface = &(UVWTfaceKproxy->GetField_STR());
    }
    
    int il1, jl1, kl1;
    GetNsurfIndex(nsurf, il1, jl1, kl1);

    int ist = 1;
    int jst = 1;
    int kst = 1;
    int ied = ni - 1 + il1;
    int jed = nj - 1 + jl1;
    int ked = nk - 1 + kl1;    
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m <= 3; ++m)
                {
                    RDouble a = UVWT(i - il1*3, j - jl1*3, k - kl1*3, m);
                    RDouble b = UVWT(i - il1*2, j - jl1*2, k - kl1*2, m);
                    RDouble c = UVWT(i - il1,   j - jl1,   k - kl1,   m);
                    RDouble d = UVWT(i,         j,         k,         m);
                    RDouble e = UVWT(i + il1,   j + jl1,   k + kl1,   m);
                    RDouble f = UVWT(i + il1*2, j + jl1*2, k + kl1*2, m);
                    (*UVWTface)(i, j, k, m) = (150.0 * (c + d) - 25.0 * (b + e) + 3.0 * (a + f)) / 256.0;
                }
                if ((*UVWTface)(i, j, k, 3) < SMALL)
                {
                    (*UVWTface)(i, j, k, 3) = half * (UVWT(i - il1, j - jl1, k - kl1,  3) + UVWT(i, j, k, 3));
                }
            }
        }
    }   
}

void LESSolverStruct::GetdqdkxietactaCellCenter(Grid *gridIn, int nsurf)
{
    string highordersolvername = GlobalDataBase::GetStrParaFromDB("str_highorder_solver");
    if (highordersolvername.substr(0, 4) == "HDCS")
    {
        GetdqdkxietactaCellCenter_HDCS(gridIn, nsurf);
        return;
    }

    StructGrid *strgrid = StructGridCast(gridIn);
    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    RDouble4D *UVWTface;
    Int1D *DiffOrder;
    RDouble4D *dqdkxietacta;
    if (nsurf == 1)
    {
        UVWTface = &(UVWTfaceIproxy->GetField_STR());
        DiffOrder = strgrid->DiscretePrecisionKXI;
        dqdkxietacta = &(*reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdkxi")));
    }
    else if (nsurf == 2)
    {
        UVWTface = &(UVWTfaceJproxy->GetField_STR());
        DiffOrder = strgrid->DiscretePrecisionETA;
        dqdkxietacta = &(*reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdeta")));
    }
    else
    {
        UVWTface = &(UVWTfaceKproxy->GetField_STR());
        DiffOrder = strgrid->DiscretePrecisionCTA;
        dqdkxietacta = &(*reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdcta")));
    }

    int ist = 1;
    int jst = 1;
    int kst = 1;
    int ied = ni - 1;
    int jed = nj - 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    int il1, jl1, kl1;
    GetNsurfIndex(nsurf, il1, jl1, kl1);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int localcellID = i * il1 + j * jl1 + k * kl1;
                
                if ((*DiffOrder)(localcellID) == 2)
                {
                    for (int m = 0; m <= 3; ++ m)
                    {
                        RDouble q1 = (*UVWTface)(i      , j      , k      , m);
                        RDouble q2 = (*UVWTface)(i + il1, j + jl1, k + kl1, m);

                        (*dqdkxietacta)(i, j, k, m) = q2 - q1;
                    }
                }
                else
                {
                    for (int m = 0; m <= 3; ++ m)
                    {
                        RDouble q0 = (*UVWTface)(i - il1    , j - jl1    , k - kl1    , m);
                        RDouble q1 = (*UVWTface)(i          , j          , k          , m);
                        RDouble q2 = (*UVWTface)(i + il1    , j + jl1    , k + kl1    , m);
                        RDouble q3 = (*UVWTface)(i + il1 * 2, j + jl1 * 2, k + kl1 * 2, m);

                        (*dqdkxietacta)(i, j, k, m) = (27.0 * (q2 - q1) - (q3 - q0)) / 24.0;
                    }
                }
            }
        }
    }   
}

void LESSolverStruct::GetdqdkxietactaCellCenter_HDCS(Grid *gridIn, int nsurf)
{
    StructGrid *strgrid = StructGridCast(gridIn);
    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    const RDouble alfa_sixth =  64.0 / 45.0;
    const RDouble    a_sixth = - 2.0 / 9.0;
    const RDouble    b_sixth =   1.0 / 180.0;

    RDouble4D &UVWT = UVWTproxy->GetField_STR();
    RDouble4D *UVWTface;
    RDouble4D *dqdkxietacta;
    if (nsurf == 1)
    {
        UVWTface = &(UVWTfaceIproxy->GetField_STR());
        dqdkxietacta = &(*reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdkxi")));
    }
    else if (nsurf == 2)
    {
        UVWTface = &(UVWTfaceJproxy->GetField_STR());
        dqdkxietacta = &(*reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdeta")));
    }
    else
    {
        UVWTface = &(UVWTfaceKproxy->GetField_STR());
        dqdkxietacta = &(*reinterpret_cast <RDouble4D *> (strgrid->GetDataPtr("dqdcta")));
    }

    int ist = 1;
    int jst = 1;
    int kst = 1;
    int ied = ni - 1;
    int jed = nj - 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    int il1, jl1, kl1;
    GetNsurfIndex(nsurf, il1, jl1, kl1);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m <= 3; ++ m)
                {
                    RDouble fL = (*UVWTface)(i          , j          , k          , m);
                    RDouble fR = (*UVWTface)(i + il1    , j + jl1    , k + kl1    , m);

                    RDouble CL2 = UVWT(i - il1*2, j - jl1*2, k - kl1*2, m);
                    RDouble CL1 = UVWT(i - il1,   j - jl1,   k - kl1,   m);                    
                    RDouble CR1 = UVWT(i + il1,   j + jl1,   k + kl1,   m);
                    RDouble CR2 = UVWT(i + il1*2, j + jl1*2, k + kl1*2, m);
                    
                    (*dqdkxietacta)(i, j, k, m) = alfa_sixth * (fR - fL) + a_sixth * (CR1 - CL1) + b_sixth * (CR2 - CL2);
                }
            }
        }
    }   
}

void LESSolverStruct::UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy)
{

    StructGrid *grid = StructGridCast(gridIn);

    ObtainViscosity(grid);

    ObtainBoundaryValue(grid);
}

void LESSolverStruct::ObtainViscosity(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    ComputeStrainRateTensor(grid);

    ComputeWallFunction(grid);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int subgridScaleModel = parameters->GetSubgridScaleModel();

    int iterationStep;
    int nMGLevel = parameters->GetNMGLevel();        
    if (nMGLevel <= 1)
    {
        //! Finest grid.
        GlobalDataBase::GetData("outnstep", &iterationStep, PHINT, 1);
    }
    else
    {
        //! Coarse grid.
        GlobalDataBase::GetData("newnstep", &iterationStep, PHINT, 1);
    }

    int cflNstep = parameters->GetCFLVaryStep();

    if (iterationStep < cflNstep)
    {
        Smagorinsky(grid);

        return;
    }

    if (subgridScaleModel == 1)
    {
        Smagorinsky(grid);
    }
    else if (subgridScaleModel == 2)
    {
        DynamicSmagViscosityCompressibleFD(grid);
    }
    else if (subgridScaleModel == 3)
    {
        WALE(grid);
    }
    else if (subgridScaleModel == 6)
    {
        Sigma(grid);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("subgridScaleModel", subgridScaleModel);
    }
}

void LESSolverStruct::Smagorinsky(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    RDouble3D &cellLengthScale = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("cellLengthScale"));
    RDouble3D &wallFunction = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));

    RDouble3D &viscousTurbulence = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble3D &subgridScaleEnergy = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("subgridScaleEnergy"));
    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));

    RDouble4D &q = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble4D &strainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));

    Param_LESSolverStruct *parameters = GetControlParameters();

    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    int monitorVistmax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");

    RDouble refReNumber  = parameters->GetRefReNumber();

    RDouble smagConstant = parameters->GetSmagConstant();

    RDouble iConstant = parameters->GetIsotropicConstant();

    RDouble viscousTurbulenceMaximum = 0.0;
    int imax = 1;
    int jmax = 1;
    int kmax = 1;

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble deltaBar = cellLengthScale(i, j, k);

                RDouble dampingFunction = wallFunction(i, j, k);

                RDouble rho = ABS(q(i, j, k, IR)) + SMALL;

                RDouble Strain = strainRateTensor(i, j, k, 0);

                RDouble mut = (smagConstant * deltaBar * dampingFunction) * (smagConstant * deltaBar * dampingFunction) * rho * Strain * refReNumber;
                RDouble kSGS = 2.0 * dampingFunction * dampingFunction * iConstant * rho * Strain * Strain * refReNumber;

                viscousTurbulence(i, j, k) = MIN(eddyViscosityLimit, mut);
                subgridScaleEnergy(i, j, k) = kSGS;

                if (viscousTurbulenceMaximum < mut)
                {
                    viscousTurbulenceMaximum = mut;
                    imax = i;
                    jmax = j;
                    kmax = k;
                }
            }
        }
    }

    RDouble prandtlTurbulence = parameters->GetPrandtlTurbulence();

    turbulentPrandtlNumber = prandtlTurbulence;

    RDouble vist_max = viscousTurbulenceMaximum;
    grid->UpdateData("vist_max", &vist_max, PHDOUBLE, 1);

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorVistmax)
        {
            cout << "vistmax = " << viscousTurbulenceMaximum << " " << imax << " " << jmax << " " << kmax << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif
}

void LESSolverStruct::DynamicSmagViscosity(Grid *grid_in)
{
    using namespace IDX;

    RDouble u1, u2, u3;
    RDouble s11, s22, s33, s12, s13, s23;
    RDouble sij2, smod;
    RDouble dudx, dudy, dudz;
    RDouble dvdx, dvdy, dvdz;
    RDouble dwdx, dwdy, dwdz;
    RDouble rho, mut, mul;

    RDouble rL11, rL22, rL33, rL12, rL13, rL23;
    RDouble rM11, rM22, rM33, rM12, rM13, rM23;
    RDouble tfss, tfsijtfsij;

    StructGrid *grid = StructGridCast(grid_in);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();
    int nis1 = ni - 1;
    int njs1 = nj - 1;
    int nks1 = nk - 1;

    RFloat4D & q    = * reinterpret_cast<RFloat4D *> (grid->GetDataPtr("q" ));
    RFloat3D & visl = * reinterpret_cast<RFloat3D *> (grid->GetDataPtr("visl"));
    RFloat3D & vist = * reinterpret_cast<RFloat3D *> (grid->GetDataPtr("vist"));

    Param_LESSolverStruct *parameters = GetControlParameters();

    RDouble refReNumber = parameters->GetRefReNumber();

    RDouble rfwr = parameters->GetTestFilterScale();

    int turbViscousCutType = parameters->GetTurbViscousCutType();

    int *filterDirection = parameters->GetFilterDirection();

    int *averageDirection = parameters->GetAverageDirection();

    Range I(1, ni-1);
    Range J(1, nj-1);
    Range K(1, nk-1);
    if (nk == 1) K.setRange(1, 1);

    int ist, ied, jst, jed, kst, ked;
    GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

    Range II(-1, ni+1);
    Range JJ(-1, nj+1);
    Range KK(-1, nk+1);

    RFloat3D gradux(II, JJ, KK, fortranArray);
    RFloat3D graduy(II, JJ, KK, fortranArray);
    RFloat3D graduz(II, JJ, KK, fortranArray);

    RFloat3D gradvx(II, JJ, KK, fortranArray);
    RFloat3D gradvy(II, JJ, KK, fortranArray);
    RFloat3D gradvz(II, JJ, KK, fortranArray);

    RFloat3D gradwx(II, JJ, KK, fortranArray);
    RFloat3D gradwy(II, JJ, KK, fortranArray);
    RFloat3D gradwz(II, JJ, KK, fortranArray);

    GradCenter(grid, q, gradux, graduy, graduz, IU);
    GradCenter(grid, q, gradvx, gradvy, gradvz, IV);
    GradCenter(grid, q, gradwx, gradwy, gradwz, IW);

    Range III(0, ni);
    Range JJJ(0, nj);
    Range KKK(0, nk);

    RFloat3D *tfu1, *tfu2, *tfu3;
    RFloat3D *tfuifuj11, *tfuifuj22, *tfuifuj33, *tfuifuj12, *tfuifuj13, *tfuifuj23;
    RFloat3D *tfs11, *tfs22, *tfs33, *tfs12, *tfs13, *tfs23;
    RFloat3D *tfsfs11, *tfsfs22, *tfsfs33, *tfsfs12, *tfsfs13, *tfsfs23;
    RFloat3D *rLijMij, *rMijMij;
    RFloat3D *tfsmod;
    RFloat3D *rLijMijavg, *rMijMijavg;
    RFloat3D *cdelta2;

    RFloat1D *fs11, *fs22, *fs33, *fs12, *fs13, *fs23;
    RFloat1D *fsfs11, *fsfs22, *fsfs33, *fsfs12, *fsfs13, *fsfs23;
    RFloat1D *fu1, *fu2, *fu3;
    RFloat1D *fuifuj11, *fuifuj22, *fuifuj33, *fuifuj12, *fuifuj13, *fuifuj23;

    tfu1 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfu2 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfu3 = new RFloat3D(III, JJJ, KKK, fortranArray);

    tfuifuj11 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfuifuj22 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfuifuj33 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfuifuj12 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfuifuj13 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfuifuj23 = new RFloat3D(III, JJJ, KKK, fortranArray);

    tfs11 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfs22 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfs33 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfs12 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfs13 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfs23 = new RFloat3D(III, JJJ, KKK, fortranArray);

    tfsfs11 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfsfs22 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfsfs33 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfsfs12 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfsfs13 = new RFloat3D(III, JJJ, KKK, fortranArray);
    tfsfs23 = new RFloat3D(III, JJJ, KKK, fortranArray);

    rLijMij = new RFloat3D(III, JJJ, KKK, fortranArray);
    rMijMij = new RFloat3D(III, JJJ, KKK, fortranArray);

    tfsmod = new RFloat3D(III, JJJ, KKK, fortranArray);

    rLijMijavg = new RFloat3D(III, JJJ, KKK, fortranArray);
    rMijMijavg = new RFloat3D(III, JJJ, KKK, fortranArray);

    cdelta2 = new RFloat3D(III, JJJ, KKK, fortranArray);

    Range IJKMAX (0, MAX(MAX(ni, nj), nk));

    fs11 = new RFloat1D(IJKMAX, fortranArray);
    fs22 = new RFloat1D(IJKMAX, fortranArray);
    fs33 = new RFloat1D(IJKMAX, fortranArray);
    fs12 = new RFloat1D(IJKMAX, fortranArray);
    fs13 = new RFloat1D(IJKMAX, fortranArray);
    fs23 = new RFloat1D(IJKMAX, fortranArray);
    
    fsfs11 = new RFloat1D(IJKMAX, fortranArray);
    fsfs22 = new RFloat1D(IJKMAX, fortranArray);
    fsfs33 = new RFloat1D(IJKMAX, fortranArray);
    fsfs12 = new RFloat1D(IJKMAX, fortranArray);
    fsfs13 = new RFloat1D(IJKMAX, fortranArray);
    fsfs23 = new RFloat1D(IJKMAX, fortranArray);

    fu1 = new RFloat1D(IJKMAX, fortranArray);
    fu2 = new RFloat1D(IJKMAX, fortranArray);
    fu3 = new RFloat1D(IJKMAX, fortranArray);

    fuifuj11 = new RFloat1D(IJKMAX, fortranArray);
    fuifuj22 = new RFloat1D(IJKMAX, fortranArray);
    fuifuj33 = new RFloat1D(IJKMAX, fortranArray);
    fuifuj12 = new RFloat1D(IJKMAX, fortranArray);
    fuifuj13 = new RFloat1D(IJKMAX, fortranArray);
    fuifuj23 = new RFloat1D(IJKMAX, fortranArray);

    //! Prepare variables to be filtered
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {

                u1 = q(i, j, k, IU);    //! q(j,k,i,2:4) = F(u(1:3)) --- IMPLICIT Filter!!!!
                u2 = q(i, j, k, IV);
                u3 = q(i, j, k, IW);

                dudx = gradux(i, j, k);
                dudy = graduy(i, j, k);
                dudz = graduz(i, j, k);

                dvdx = gradvx(i, j, k);
                dvdy = gradvy(i, j, k);
                dvdz = gradvz(i, j, k);

                dwdx = gradwx(i, j, k);
                dwdy = gradwy(i, j, k);
                dwdz = gradwz(i, j, k);

                s11 = dudx;
                s22 = dvdy;
                s33 = dwdz;
                s12 = half * (dudy + dvdx);
                s13 = half * (dudz + dwdx);
                s23 = half * (dvdz + dwdy);

                sij2 = s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23);
                smod = sqrt(2.0 * sij2);
                //! F(ui)
                (*tfu1)(i, j, k) = u1;
                (*tfu2)(i, j, k) = u2;
                (*tfu3)(i, j, k) = u3;
                //! F(ui)F(uj)
                (*tfuifuj11)(i, j, k) = u1 * u1;
                (*tfuifuj22)(i, j, k) = u2 * u2;
                (*tfuifuj33)(i, j, k) = u3 * u3;
                (*tfuifuj12)(i, j, k) = u1 * u2;
                (*tfuifuj13)(i, j, k) = u1 * u3;
                (*tfuifuj23)(i, j, k) = u2 * u3;
                //! F(Sij)
                (*tfs11)(i, j, k) = s11;
                (*tfs22)(i, j, k) = s22;
                (*tfs33)(i, j, k) = s33;
                (*tfs12)(i, j, k) = s12;
                (*tfs13)(i, j, k) = s13;
                (*tfs23)(i, j, k) = s23;
                //! F(|S|)
                (*tfsmod)(i, j, k) = smod;
                //! F(|S|)F(Sij)
                (*tfsfs11)(i, j, k) = smod * s11;
                (*tfsfs22)(i, j, k) = smod * s22;
                (*tfsfs33)(i, j, k) = smod * s33;
                (*tfsfs12)(i, j, k) = smod * s12;
                (*tfsfs13)(i, j, k) = smod * s13;
                (*tfsfs23)(i, j, k) = smod * s23;
            }
        }
    }

    RDouble a = (*tfs11)(1, 1, 1);
    RDouble b = (*tfs22)(1, 1, 1);
    RDouble c = (*tfs33)(1, 1, 1);
    RDouble d = (*tfs12)(1, 1, 1);
    RDouble e = (*tfs13)(1, 1, 1);
    RDouble f = (*tfs23)(1, 1, 1);
    RDouble g = (*tfsfs11)(1, 1, 1);
    RDouble h = (*tfsfs22)(1, 1, 1);
    RDouble p = (*tfsfs33)(1, 1, 1);
    RDouble r = (*tfsfs12)(1, 1, 1);
    RDouble s = (*tfsfs13)(1, 1, 1);
    
    //! Get filtered variables in I direction
    if (filterDirection[0] == 1)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    //! F(ui)
                    (*fu1)(i) = (*tfu1)(i, j, k);
                    (*fu2)(i) = (*tfu2)(i, j, k);
                    (*fu3)(i) = (*tfu3)(i, j, k);
                    //! F(ui)F(uj)
                    (*fuifuj11)(i) = (*tfuifuj11)(i, j, k);
                    (*fuifuj22)(i) = (*tfuifuj22)(i, j, k);
                    (*fuifuj33)(i) = (*tfuifuj33)(i, j, k);
                    (*fuifuj12)(i) = (*tfuifuj12)(i, j, k);
                    (*fuifuj13)(i) = (*tfuifuj13)(i, j, k);
                    (*fuifuj23)(i) = (*tfuifuj23)(i, j, k);
                    //! F(Sij)
                    (*fs11)(i) = (*tfs11)(i, j, k);
                    (*fs22)(i) = (*tfs22)(i, j, k);
                    (*fs33)(i) = (*tfs33)(i, j, k);
                    (*fs12)(i) = (*tfs12)(i, j, k);
                    (*fs13)(i) = (*tfs13)(i, j, k);
                    (*fs23)(i) = (*tfs23)(i, j, k);
                    //! F(|S|)F(Sij)
                    (*fsfs11)(i) = (*tfsfs11)(i, j, k);
                    (*fsfs22)(i) = (*tfsfs22)(i, j, k);
                    (*fsfs33)(i) = (*tfsfs33)(i, j, k);
                    (*fsfs12)(i) = (*tfsfs12)(i, j, k);
                    (*fsfs13)(i) = (*tfsfs13)(i, j, k);
                    (*fsfs23)(i) = (*tfsfs23)(i, j, k);
                }
                //! left boundary:
                int i0 = 0;

                u1 = q(i0, j, k, IU);
                u2 = q(i0, j, k, IV);
                u3 = q(i0, j, k, IW);
                //! F(ui)
                (*fu1)(i0) = u1;
                (*fu2)(i0) = u2;
                (*fu3)(i0) = u3;
                //! F(ui)F(uj)
                (*fuifuj11)(i0) = u1 * u1;
                (*fuifuj22)(i0) = u2 * u2;
                (*fuifuj33)(i0) = u3 * u3;
                (*fuifuj12)(i0) = u1 * u2;
                (*fuifuj13)(i0) = u1 * u3;
                (*fuifuj23)(i0) = u2 * u3;
                //! F(Sij) 
                (*fs11)(i0) = (*fs11)(1);
                (*fs22)(i0) = (*fs22)(1);
                (*fs33)(i0) = (*fs33)(1);
                (*fs12)(i0) = (*fs12)(1);
                (*fs13)(i0) = (*fs13)(1);
                (*fs23)(i0) = (*fs23)(1);
                //! F(|S|)F(Sij)
                (*fsfs11)(i0) = (*fsfs11)(1);
                (*fsfs22)(i0) = (*fsfs22)(1);
                (*fsfs33)(i0) = (*fsfs33)(1);
                (*fsfs12)(i0) = (*fsfs12)(1);
                (*fsfs13)(i0) = (*fsfs13)(1);
                (*fsfs23)(i0) = (*fsfs23)(1);
                //! Right boundary:
                u1 = q(ni, j, k, IU);
                u2 = q(ni, j, k, IV);
                u3 = q(ni, j, k, IW);
                //! F(ui)
                (*fu1)(ni) = u1;
                (*fu2)(ni) = u2;
                (*fu3)(ni) = u3;
                //! F(ui)F(uj)
                (*fuifuj11)(ni) = u1 * u1;
                (*fuifuj22)(ni) = u2 * u2;
                (*fuifuj33)(ni) = u3 * u3;
                (*fuifuj12)(ni) = u1 * u2;
                (*fuifuj13)(ni) = u1 * u3;
                (*fuifuj23)(ni) = u2 * u3;
                //! F(Sij) 
                (*fs11)(ni) = (*fs11)(nis1);
                (*fs22)(ni) = (*fs22)(nis1);
                (*fs33)(ni) = (*fs33)(nis1);
                (*fs12)(ni) = (*fs12)(nis1);
                (*fs13)(ni) = (*fs13)(nis1);
                (*fs23)(ni) = (*fs23)(nis1);
                //! F(|S|)F(Sij)
                (*fsfs11)(ni) = (*fsfs11)(nis1);
                (*fsfs22)(ni) = (*fsfs22)(nis1);
                (*fsfs33)(ni) = (*fsfs33)(nis1);
                (*fsfs12)(ni) = (*fsfs12)(nis1);
                (*fsfs13)(ni) = (*fsfs13)(nis1);
                (*fsfs23)(ni) = (*fsfs23)(nis1);
                //! Filter
                filter(ni, fu1, fu2, fu3, fuifuj11, fuifuj22, fuifuj33, fuifuj12, fuifuj13, fuifuj23, fs11, fs22, fs33, fs12, fs13, fs23, fsfs11, fsfs22, fsfs33, fsfs12, fsfs13, fsfs23);

                for (int i = ist; i <= ied; ++ i)
                {
                    //! Fill TF(ui)
                    (*tfu1)(i, j, k) = (*fu1)(i);
                    (*tfu2)(i, j, k) = (*fu2)(i);
                    (*tfu3)(i, j, k) = (*fu3)(i);
                    //! Fill T(F(ui)F(uj))
                    (*tfuifuj11)(i, j, k) = (*fuifuj11)(i);
                    (*tfuifuj22)(i, j, k) = (*fuifuj22)(i);
                    (*tfuifuj33)(i, j, k) = (*fuifuj33)(i);
                    (*tfuifuj12)(i, j, k) = (*fuifuj12)(i);
                    (*tfuifuj13)(i, j, k) = (*fuifuj13)(i);
                    (*tfuifuj23)(i, j, k) = (*fuifuj23)(i);
                    //! Get T(F(Sij))
                    (*tfs11)(i, j, k) = (*fs11)(i);
                    (*tfs22)(i, j, k) = (*fs22)(i);
                    (*tfs33)(i, j, k) = (*fs33)(i);
                    (*tfs12)(i, j, k) = (*fs12)(i);
                    (*tfs13)(i, j, k) = (*fs13)(i);
                    (*tfs23)(i, j, k) = (*fs23)(i);
                    //! Get T(F(|S|)F(Sij))
                    (*tfsfs11)(i, j, k) = (*fsfs11)(i);
                    (*tfsfs22)(i, j, k) = (*fsfs22)(i);
                    (*tfsfs33)(i, j, k) = (*fsfs33)(i);
                    (*tfsfs12)(i, j, k) = (*fsfs12)(i);
                    (*tfsfs13)(i, j, k) = (*fsfs13)(i);
                    (*tfsfs23)(i, j, k) = (*fsfs23)(i);
                }
            }
        }
    } //! filterDirection[0] == 1

    a = (*tfs11)(1, 1, 1);
    b = (*tfs22)(1, 1, 1);
    c = (*tfs33)(1, 1, 1);
    d = (*tfs12)(1, 1, 1);
    e = (*tfs13)(1, 1, 1);
    f = (*tfs23)(1, 1, 1);
    g = (*tfsfs11)(1, 1, 1);
    h = (*tfsfs22)(1, 1, 1);
    p = (*tfsfs33)(1, 1, 1);
    r = (*tfsfs12)(1, 1, 1);
    s = (*tfsfs13)(1, 1, 1);

// get filtered variables in J direction
    if (filterDirection[1] == 1) 
    {
        for (int i = ist; i <= ied; ++ i)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    //! F(ui)
                    (*fu1)(j) = (*tfu1)(i, j, k);
                    (*fu2)(j) = (*tfu2)(i, j, k);
                    (*fu3)(j) = (*tfu3)(i, j, k);
                    //! F(ui)F(uj)
                    (*fuifuj11)(j) = (*tfuifuj11)(i, j, k);
                    (*fuifuj22)(j) = (*tfuifuj22)(i, j, k);
                    (*fuifuj33)(j) = (*tfuifuj33)(i, j, k);
                    (*fuifuj12)(j) = (*tfuifuj12)(i, j, k);
                    (*fuifuj13)(j) = (*tfuifuj13)(i, j, k);
                    (*fuifuj23)(j) = (*tfuifuj23)(i, j, k);
                    //! F(Sij)
                    (*fs11)(j) = (*tfs11)(i, j, k);
                    (*fs22)(j) = (*tfs22)(i, j, k);
                    (*fs33)(j) = (*tfs33)(i, j, k);
                    (*fs12)(j) = (*tfs12)(i, j, k);
                    (*fs13)(j) = (*tfs13)(i, j, k);
                    (*fs23)(j) = (*tfs23)(i, j, k);
                    //! F(|S|)F(Sij)
                    (*fsfs11)(j) = (*tfsfs11)(i, j, k);
                    (*fsfs22)(j) = (*tfsfs22)(i, j, k);
                    (*fsfs33)(j) = (*tfsfs33)(i, j, k);
                    (*fsfs12)(j) = (*tfsfs12)(i, j, k);
                    (*fsfs13)(j) = (*tfsfs13)(i, j, k);
                    (*fsfs23)(j) = (*tfsfs23)(i, j, k);
                }
                //! Left boundary:
                int j0 = 0;

                u1 = q(i, j0, k, IU);
                u2 = q(i, j0, k, IV);
                u3 = q(i, j0, k, IW);
                //! F(ui)
                (*fu1)(j0) = u1;
                (*fu2)(j0) = u2;
                (*fu3)(j0) = u3;
                //! F(ui)F(uj)
                (*fuifuj11)(j0) = u1 * u1;
                (*fuifuj22)(j0) = u2 * u2;
                (*fuifuj33)(j0) = u3 * u3;
                (*fuifuj12)(j0) = u1 * u2;
                (*fuifuj13)(j0) = u1 * u3;
                (*fuifuj23)(j0) = u2 * u3;
                //! F(Sij) 
                (*fs11)(j0) = (*fs11)(1);
                (*fs22)(j0) = (*fs22)(1);
                (*fs33)(j0) = (*fs33)(1);
                (*fs12)(j0) = (*fs12)(1);
                (*fs13)(j0) = (*fs13)(1);
                (*fs23)(j0) = (*fs23)(1);
                //! F(|S|)F(Sij)
                (*fsfs11)(j0) = (*fsfs11)(1);
                (*fsfs22)(j0) = (*fsfs22)(1);
                (*fsfs33)(j0) = (*fsfs33)(1);
                (*fsfs12)(j0) = (*fsfs12)(1);
                (*fsfs13)(j0) = (*fsfs13)(1);
                (*fsfs23)(j0) = (*fsfs23)(1);
                //! Right boundary:
                u1 = q(i, nj, k, IU);
                u2 = q(i, nj, k, IV);
                u3 = q(i, nj, k, IW);
                //! F(ui)
                (*fu1)(nj) = u1;
                (*fu2)(nj) = u2;
                (*fu3)(nj) = u3;
                //! F(ui)F(uj)
                (*fuifuj11)(nj) = u1 * u1;
                (*fuifuj22)(nj) = u2 * u2;
                (*fuifuj33)(nj) = u3 * u3;
                (*fuifuj12)(nj) = u1 * u2;
                (*fuifuj13)(nj) = u1 * u3;
                (*fuifuj23)(nj) = u2 * u3;
                //! F(Sij)
                (*fs11)(nj) = (*fs11)(njs1);
                (*fs22)(nj) = (*fs22)(njs1);
                (*fs33)(nj) = (*fs33)(njs1);
                (*fs12)(nj) = (*fs12)(njs1);
                (*fs13)(nj) = (*fs13)(njs1);
                (*fs23)(nj) = (*fs23)(njs1);
                //! F(|S|)F(Sij)
                (*fsfs11)(nj) = (*fsfs11)(njs1);
                (*fsfs22)(nj) = (*fsfs22)(njs1);
                (*fsfs33)(nj) = (*fsfs33)(njs1);
                (*fsfs12)(nj) = (*fsfs12)(njs1);
                (*fsfs13)(nj) = (*fsfs13)(njs1);
                (*fsfs23)(nj) = (*fsfs23)(njs1);        
                //! Filter
                filter(nj, fu1, fu2, fu3, fuifuj11, fuifuj22, fuifuj33, fuifuj12, fuifuj13, fuifuj23, fs11, fs22, fs33, fs12, fs13, fs23, fsfs11, fsfs22, fsfs33, fsfs12, fsfs13, fsfs23);

                for (int j = jst; j <= jed; ++ j)
                {
                    //! Fill TF(ui)
                    (*tfu1)(i,j,k) = (*fu1)(j);
                    (*tfu2)(i,j,k) = (*fu2)(j);
                    (*tfu3)(i,j,k) = (*fu3)(j);
                    //! Fill T(F(ui)F(uj))
                    (*tfuifuj11)(i,j,k) = (*fuifuj11)(j);
                    (*tfuifuj22)(i,j,k) = (*fuifuj22)(j);
                    (*tfuifuj33)(i,j,k) = (*fuifuj33)(j);
                    (*tfuifuj12)(i,j,k) = (*fuifuj12)(j);
                    (*tfuifuj13)(i,j,k) = (*fuifuj13)(j);
                    (*tfuifuj23)(i,j,k) = (*fuifuj23)(j);
                    //! Get T(F(Sij))
                    (*tfs11)(i,j,k) = (*fs11)(j);
                    (*tfs22)(i,j,k) = (*fs22)(j);
                    (*tfs33)(i,j,k) = (*fs33)(j);
                    (*tfs12)(i,j,k) = (*fs12)(j);
                    (*tfs13)(i,j,k) = (*fs13)(j);
                    (*tfs23)(i,j,k) = (*fs23)(j);
                    //! Get T(F(|S|)F(Sij))
                    (*tfsfs11)(i,j,k) = (*fsfs11)(j);
                    (*tfsfs22)(i,j,k) = (*fsfs22)(j);
                    (*tfsfs33)(i,j,k) = (*fsfs33)(j);
                    (*tfsfs12)(i,j,k) = (*fsfs12)(j);
                    (*tfsfs13)(i,j,k) = (*fsfs13)(j);
                    (*tfsfs23)(i,j,k) = (*fsfs23)(j);
                }
            }
        }

    } //! filterDirection[1] == 0

    a = (*tfs11)(1, 1, 1);
    b = (*tfs22)(1, 1, 1);
    c = (*tfs33)(1, 1, 1);
    d = (*tfs12)(1, 1, 1);
    e = (*tfs13)(1, 1, 1);
    f = (*tfs23)(1, 1, 1);
    g = (*tfsfs11)(1, 1, 1);
    h = (*tfsfs22)(1, 1, 1);
    p = (*tfsfs33)(1, 1, 1);
    r = (*tfsfs12)(1, 1, 1);
    s = (*tfsfs13)(1, 1, 1);

    //! Get filtered variables in K direction
    if (filterDirection[2] == 1) 
    {
        for (int i = ist; i <= ied; ++ i)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    //! F(ui)
                    (*fu1)(k) = (*tfu1)(i,j,k);
                    (*fu2)(k) = (*tfu2)(i,j,k);
                    (*fu3)(k) = (*tfu3)(i,j,k);
                    //! F(ui)F(uj)
                    (*fuifuj11)(k) = (*tfuifuj11)(i,j,k);
                    (*fuifuj22)(k) = (*tfuifuj22)(i,j,k);
                    (*fuifuj33)(k) = (*tfuifuj33)(i,j,k);
                    (*fuifuj12)(k) = (*tfuifuj12)(i,j,k);
                    (*fuifuj13)(k) = (*tfuifuj13)(i,j,k);
                    (*fuifuj23)(k) = (*tfuifuj23)(i,j,k);
                    //! F(Sij)
                    (*fs11)(k) = (*tfs11)(i,j,k);
                    (*fs22)(k) = (*tfs22)(i,j,k);
                    (*fs33)(k) = (*tfs33)(i,j,k);
                    (*fs12)(k) = (*tfs12)(i,j,k);
                    (*fs13)(k) = (*tfs13)(i,j,k);
                    (*fs23)(k) = (*tfs23)(i,j,k);
                    //! F(|S|)F(Sij)
                    (*fsfs11)(k) = (*tfsfs11)(i,j,k);
                    (*fsfs22)(k) = (*tfsfs22)(i,j,k);
                    (*fsfs33)(k) = (*tfsfs33)(i,j,k);
                    (*fsfs12)(k) = (*tfsfs12)(i,j,k);
                    (*fsfs13)(k) = (*tfsfs13)(i,j,k);
                    (*fsfs23)(k) = (*tfsfs23)(i,j,k);
                }
                //! left boundary:
                int k0 = 0;

                u1 = q(i,j,k0,IU);
                u2 = q(i,j,k0,IV);
                u3 = q(i,j,k0,IW);
                //! F(ui)
                (*fu1)(k0) = u1;
                (*fu2)(k0) = u2;
                (*fu3)(k0) = u3;
                //! F(ui)F(uj)
                (*fuifuj11)(k0) = u1*u1;
                (*fuifuj22)(k0) = u2*u2;
                (*fuifuj33)(k0) = u3*u3;
                (*fuifuj12)(k0) = u1*u2;
                (*fuifuj13)(k0) = u1*u3;
                (*fuifuj23)(k0) = u2*u3;
                //! F(Sij) 
                (*fs11)(k0) = (*fs11)(1);
                (*fs22)(k0) = (*fs22)(1);
                (*fs33)(k0) = (*fs33)(1);
                (*fs12)(k0) = (*fs12)(1);
                (*fs13)(k0) = (*fs13)(1);
                (*fs23)(k0) = (*fs23)(1);
                //! F(|S|)F(Sij)
                (*fsfs11)(k0) = (*fsfs11)(1);
                (*fsfs22)(k0) = (*fsfs22)(1);
                (*fsfs33)(k0) = (*fsfs33)(1);
                (*fsfs12)(k0) = (*fsfs12)(1);
                (*fsfs13)(k0) = (*fsfs13)(1);
                (*fsfs23)(k0) = (*fsfs23)(1);
                //! Right boundary:
                u1 = q(i,j,nk,IU);
                u2 = q(i,j,nk,IV);
                u3 = q(i,j,nk,IW);
                //! F(ui)
                (*fu1)(nk) = u1;
                (*fu2)(nk) = u2;
                (*fu3)(nk) = u3;
                //! F(ui)F(uj)
                (*fuifuj11)(nk) = u1*u1;
                (*fuifuj22)(nk) = u2*u2;
                (*fuifuj33)(nk) = u3*u3;
                (*fuifuj12)(nk) = u1*u2;
                (*fuifuj13)(nk) = u1*u3;
                (*fuifuj23)(nk) = u2*u3;
                //! F(Sij)
                (*fs11)(nk) = (*fs11)(nks1);
                (*fs22)(nk) = (*fs22)(nks1);
                (*fs33)(nk) = (*fs33)(nks1);
                (*fs12)(nk) = (*fs12)(nks1);
                (*fs13)(nk) = (*fs13)(nks1);
                (*fs23)(nk) = (*fs23)(nks1);
                //! F(|S|)F(Sij)
                (*fsfs11)(nk) = (*fsfs11)(nks1);
                (*fsfs22)(nk) = (*fsfs22)(nks1);
                (*fsfs33)(nk) = (*fsfs33)(nks1);
                (*fsfs12)(nk) = (*fsfs12)(nks1);
                (*fsfs13)(nk) = (*fsfs13)(nks1);
                (*fsfs23)(nk) = (*fsfs23)(nks1);        
                //! Filter
                filter(nk,fu1,fu2,fu3,fuifuj11,fuifuj22,fuifuj33,fuifuj12,fuifuj13,fuifuj23,fs11,fs22,fs33,fs12,fs13,fs23,fsfs11,fsfs22,fsfs33,fsfs12,fsfs13,fsfs23);
        
                for (int k = kst; k <= ked; ++ k)
                {
                    //! Fill TF(ui)
                    (*tfu1)(i,j,k) = (*fu1)(k);
                    (*tfu2)(i,j,k) = (*fu2)(k);
                    (*tfu3)(i,j,k) = (*fu3)(k);
                    //! Fill T(F(ui)F(uj))
                    (*tfuifuj11)(i,j,k) = (*fuifuj11)(k);
                    (*tfuifuj22)(i,j,k) = (*fuifuj22)(k);
                    (*tfuifuj33)(i,j,k) = (*fuifuj33)(k);
                    (*tfuifuj12)(i,j,k) = (*fuifuj12)(k);
                    (*tfuifuj13)(i,j,k) = (*fuifuj13)(k);
                    (*tfuifuj23)(i,j,k) = (*fuifuj23)(k);
                    //! Get T(F(Sij))
                    (*tfs11)(i,j,k) = (*fs11)(k);
                    (*tfs22)(i,j,k) = (*fs22)(k);
                    (*tfs33)(i,j,k) = (*fs33)(k);
                    (*tfs12)(i,j,k) = (*fs12)(k);
                    (*tfs13)(i,j,k) = (*fs13)(k);
                    (*tfs23)(i,j,k) = (*fs23)(k);
                    //! Get T(F(|S|)F(Sij))
                    (*tfsfs11)(i,j,k) = (*fsfs11)(k);
                    (*tfsfs22)(i,j,k) = (*fsfs22)(k);
                    (*tfsfs33)(i,j,k) = (*fsfs33)(k);
                    (*tfsfs12)(i,j,k) = (*fsfs12)(k);
                    (*tfsfs13)(i,j,k) = (*fsfs13)(k);
                    (*tfsfs23)(i,j,k) = (*fsfs23)(k);
                }
            }
        }

    } //! filterDirection[2] == 1

    //! Get LijMij, MijMij
    for (int i = ist; i <= ied; ++ i) 
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                //! Lij = T(F(ui)F(uj)) - TF(ui)*TF(uj)
                rL11 = (*tfuifuj11)(i,j,k) - (*tfu1)(i,j,k) * (*tfu1)(i,j,k);
                rL22 = (*tfuifuj22)(i,j,k) - (*tfu2)(i,j,k) * (*tfu2)(i,j,k);
                rL33 = (*tfuifuj33)(i,j,k) - (*tfu3)(i,j,k) * (*tfu3)(i,j,k);
                rL12 = (*tfuifuj12)(i,j,k) - (*tfu1)(i,j,k) * (*tfu2)(i,j,k);
                rL13 = (*tfuifuj13)(i,j,k) - (*tfu1)(i,j,k) * (*tfu3)(i,j,k);
                rL23 = (*tfuifuj23)(i,j,k) - (*tfu2)(i,j,k) * (*tfu3)(i,j,k);
                //! Mij = (TF(Detla)/F(Delta))^2*|TF(S)|*T(F(Sij)) -  T(F(|S|)F(Sij))
                //! TF(|S|)
                tfsijtfsij = (*tfs11)(i,j,k) * (*tfs11)(i,j,k) 
                           + (*tfs22)(i,j,k) * (*tfs22)(i,j,k)
                           + (*tfs33)(i,j,k) * (*tfs33)(i,j,k)
                           + 2.0*( (*tfs12)(i,j,k) * (*tfs12)(i,j,k)
                                  + (*tfs13)(i,j,k) * (*tfs13)(i,j,k)
                                  + (*tfs23)(i,j,k) * (*tfs23)(i,j,k));
                tfss = sqrt(2.0*tfsijtfsij);

                rM11 = rfwr * rfwr * tfss * (*tfs11)(i,j,k) - (*tfsfs11)(i,j,k);
                rM22 = rfwr * rfwr * tfss * (*tfs22)(i,j,k) - (*tfsfs22)(i,j,k);
                rM33 = rfwr * rfwr * tfss * (*tfs33)(i,j,k) - (*tfsfs33)(i,j,k);
                rM12 = rfwr * rfwr * tfss * (*tfs12)(i,j,k) - (*tfsfs12)(i,j,k);
                rM13 = rfwr * rfwr * tfss * (*tfs13)(i,j,k) - (*tfsfs13)(i,j,k);
                rM23 = rfwr * rfwr * tfss * (*tfs23)(i,j,k) - (*tfsfs23)(i,j,k);

                (*rLijMij)(i,j,k) = rL11*rM11 +rL22*rM22 +rL33*rM33
                                  + 2.0*(rL12*rM12 +rL13*rM13 +rL23*rM23);

                (*rMijMij)(i,j,k) = rM11*rM11 +rM22*rM22 +rM33*rM33
                                  + 2.0*(rM12*rM12 +rM13*rM13 +rM23*rM23);
            }
        }
    }

    a = (*tfs11)(1, 1, 1);
    b = (*tfsfs11)(1, 1, 1);

    //! Get <LijMij>, <MijMij> in I direction
    
    if ((averageDirection[0] == 0) && (averageDirection[1] == 0) && (averageDirection[2] == 0))
    {
        int averageWidth = parameters->GetAverageWidth();

        int naw = averageWidth;
        //! Boundary points and corner points were averaged only with points inside this block!
        //! This should be OK, since the foregoing filter operations were executed only on points in this block also!
        for (int i = ist; i <= ied; ++ i)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    RDouble rLMavg = 0.0;
                    RDouble rMMavg = 0.0;
                    int nnavg = 0;
                    for (int idi = -naw; idi <= naw; ++ idi)
                    {
                        for (int idj = -naw; j <= naw; ++ idj)
                        {
                            for (int idk = -naw; k <= naw; ++ idk)
                            {
                                int ii = i + idi;
                                int jj = j + idj;
                                int kk = k + idk;
                                bool Lvalid = (ii>=1 && ii<=nis1) && (jj>=1 && jj<=njs1) && (kk>=1 && kk<=nks1);
                                if (Lvalid) 
                                {
                                    nnavg = nnavg + 1;
                                    rLMavg = rLMavg + (*rLijMij)(jj,kk,ii);
                                    rMMavg = rMMavg + (*rMijMij)(jj,kk,ii);
                                }
                            }
                        }
                    }
                    (*rLijMijavg)(i,j,k) = rLMavg/(static_cast<RDouble>(nnavg)+EPSILON);
                    (*rMijMijavg)(i,j,k) = rMMavg/(static_cast<RDouble>(nnavg) + EPSILON);
                }
            }
        }

        for (int i = ist; i <= ied; ++ i)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    (*rLijMij)(i,j,k) = (*rLijMijavg)(i,j,k);
                    (*rMijMij)(i,j,k) = (*rMijMijavg)(i,j,k);
                }
            }
        }
    }
    else
    {
        if (averageDirection[0] != 0)
        //! Get <LijMij>, <MijMij> in I direction
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    RDouble rLMavg = 0.0;
                    RDouble rMMavg = 0.0;
                    for (int i = ist; i <= ied; ++ i)
                    {
                        rLMavg = rLMavg + (*rLijMij)(i,j,k);
                        rMMavg = rMMavg + (*rMijMij)(i,j,k);
                    }
                    for (int i = ist; i <= ied; ++ i)
                    {
                        (*rLijMij)(i,j,k) = rLMavg/nis1;
                        (*rMijMij)(i,j,k) = rMMavg/nis1;
                    }
                }
            }
        }

        RDouble lm = (*rLijMij)(1, 1, 1);
        RDouble mm = (*rMijMij)(1, 1, 1);

        if (averageDirection[1] != 0)
        //! Get <LijMij>, <MijMij> in J direction
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble rLMavg = 0.0;
                    RDouble rMMavg = 0.0;
                    for (int j = ist; j <= jed; ++ j)
                    {
                        rLMavg = rLMavg + (*rLijMij)(i,j,k);
                        rMMavg = rMMavg + (*rMijMij)(i,j,k);
                    }
                    for (int j = jst; j <= jed; ++ j)
                    {
                        (*rLijMij)(i,j,k) = rLMavg/njs1;
                        (*rMijMij)(i,j,k) = rMMavg/njs1;
                    }
                }
            }
        }

        lm = (*rLijMij)(1, 1, 1);
        mm = (*rMijMij)(1, 1, 1);

        if (averageDirection[2] != 0)
        //! Get <LijMij>, <MijMij> in K direction
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    RDouble rLMavg = 0.0;
                    RDouble rMMavg = 0.0;
                    for (int k = kst; k <= ked; ++ k)
                    {
                        rLMavg = rLMavg + (*rLijMij)(i,j,k);
                        rMMavg = rMMavg + (*rMijMij)(i,j,k);
                    }
                    for (int k = kst; k <= ked; ++ k)
                    {
                        (*rLijMij)(i,j,k) = rLMavg/nks1;
                        (*rMijMij)(i,j,k) = rMMavg/nks1;
                    }
                }
            }
        }
    }


    //! Get cdelta2, get vist(i,j,k)

    int nNegativeCr = 0;
    int nNegativeMu = 0;
    int nLowerBoundCr = 0;
    int nUpperBoundCr = 0;

    for (int i = ist; i <= ied; ++ i)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                RDouble lijmij = (*rLijMij)(i,j,k);
                RDouble mijmij = (*rMijMij)(i,j,k);
                RDouble cr = -0.5 * lijmij / mijmij;

                (*cdelta2)(i,j,k) = cr;

                u1 = q(i,j,k,IU);    //! q(j,k,i,2:4) = F(u(1:3)) --- IMPLICIT Filter!!!!
                u2 = q(i,j,k,IV);
                u3 = q(i,j,k,IW);

                dudx = gradux(i,j,k);
                dudy = graduy(i,j,k);
                dudz = graduz(i,j,k);

                dvdx = gradvx(i,j,k);
                dvdy = gradvy(i,j,k);
                dvdz = gradvz(i,j,k);

                dwdx = gradwx(i,j,k);
                dwdy = gradwy(i,j,k);
                dwdz = gradwz(i,j,k);

                s11 = dudx;
                s22 = dvdy;
                s33 = dwdz;
                s12 = half * (dudy + dvdx);
                s13 = half * (dudz + dwdx);
                s23 = half * (dvdz + dwdy);

                sij2 = s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23);
                smod = sqrt(2.0*sij2);
                rho = q(i,j,k,IR);
                mul = visl(i,j,k);

                mut = rho * refReNumber * (*cdelta2)(i,j,k) * smod;
                //! Smagorinsky nut

                if (cr < 0.0)
                {
                    nNegativeCr = nNegativeCr + 1;
                }

                if (cr < -0.01)
                {
                    //cr = -0.01;
                    nLowerBoundCr = nLowerBoundCr + 1;
                }

                if (cr > 0.5)
                {
                    //cr = 0.5;
                    nUpperBoundCr = nUpperBoundCr + 1;
                }

                if (mut < -mul)
                {
                    nNegativeMu = nNegativeMu + 1;
                }

                if (turbViscousCutType == 1)
                {
                    mut = MAX(mut - mul, 0.0);
                }
                else if (turbViscousCutType == 2)
                {
                    mut = MAX(mut, 0.0);
                }

                vist(i,j,k) = mut;
            }
        }
    }

    grid->UpdateData("nNegativeCr", &nNegativeCr, PHINT, 1);
    grid->UpdateData("nNegativeMu", &nNegativeMu, PHINT, 1);

    int monitorNegativeConstant = GlobalDataBase::GetIntParaFromDB("monitorNegativeConstant");

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorNegativeConstant)
        {
            RDouble percentNegativeCr = nNegativeCr / static_cast<RDouble> ((ni-1) * (nj-1) * (nk-1));
            RDouble percentNegativeMu = nNegativeMu / static_cast<RDouble> ((ni-1) * (nj-1) * (nk-1));
            RDouble percentLowerBoundCr = nLowerBoundCr / static_cast<RDouble> ((ni-1) * (nj-1) * (nk-1));
            RDouble percentUpperBoundCr = nUpperBoundCr / static_cast<RDouble> ((ni-1) * (nj-1) * (nk-1));

            cout << "percentNegativeCr = " << percentNegativeCr << "\n";
            //cout << "nNegativeCi = " << nNegativeCi << "\n";
            //cout << "nNegativePrt = " << nNegativePrt << "\n";
            cout << "percentNegativeMu = " << percentNegativeMu << "\n";
            cout << "percentLowerBoundCr = " << percentLowerBoundCr << "\n";
            cout << "percentUpperBoundCr = " << percentUpperBoundCr << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif

    delete tfu1      ; delete tfu2      ; delete tfu3;
    delete tfuifuj11 ; delete tfuifuj22 ; delete tfuifuj33; delete tfuifuj12; delete tfuifuj13; delete tfuifuj23;
    delete tfs11     ; delete tfs22     ; delete tfs33    ; delete tfs12    ; delete tfs13    ; delete tfs23    ;
    delete tfsfs11   ; delete tfsfs22   ; delete tfsfs33  ; delete tfsfs12  ; delete tfsfs13  ; delete tfsfs23  ;
    delete rLijMij   ; delete rMijMij   ;
    delete tfsmod    ;
    delete rLijMijavg; delete rMijMijavg;
    delete cdelta2   ;

    delete fs11      ; delete fs22      ; delete fs33     ; delete fs12     ; delete fs13     ; delete fs23     ;
    delete fsfs11    ; delete fsfs22    ; delete fsfs33   ; delete fsfs12   ; delete fsfs13   ; delete fsfs23   ;
    delete fu1       ; delete fu2       ; delete fu3      ;
    delete fuifuj11  ; delete fuifuj22  ; delete fuifuj33 ; delete fuifuj12 ; delete fuifuj13 ; delete fuifuj23 ;
}


void LESSolverStruct::DynamicSmagViscosityCompressibleOld(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_LESSolverStruct *parameters = GetControlParameters();

    RDouble refReNumber = parameters->GetRefReNumber();

    int turbViscousCutType = parameters->GetTurbViscousCutType();

    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    int monitorVistmax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");

    int monitorNegativeConstant = GlobalDataBase::GetIntParaFromDB("monitorNegativeConstant");

    int prtCutType = GlobalDataBase::GetIntParaFromDB("prtCutType");

    RDouble oprt = parameters->GetoPrandtlTurbulence();

    RDouble refGama = parameters->GetRefGama();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble specificHeatAtConstantPressure = 1.0 / ((refGama - 1.0) * refMachNumber * refMachNumber);

    RDouble4D &wallFunction = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("wallFunction"));

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &visl = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble5D &gradPrimtiveVarFaceXofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceXofLES"));
    RDouble5D &gradPrimtiveVarFaceYofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceYofLES"));
    RDouble5D &gradPrimtiveVarFaceZofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceZofLES"));

    RDouble5D &gradtx = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceXofLES"));
    RDouble5D &gradty = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceYofLES"));
    RDouble5D &gradtz = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceZofLES"));

    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));
    RDouble3D &anisotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));
    RDouble3D &isotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("isotropicConstant"));

    RDouble5D &subgridStress = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("subgridStress"));
    //RDouble4D &subgridEnergy = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("kSGS"));
    RDouble5D &subgridHeatFlux = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("subgridHeatFlux"));

    ComputeAnisotropicConstant(grid);

    ComputeIsotropicConstant(grid);

    ComputePrandtlNumber(grid);

    double vistmax  = 0.0;
    double vistmin  = 1.0e30;

    int nNegativeCr = 0;
    int nLowerBoundCr = 0;
    int nUpperBoundCr = 0;
    int nNegativeCi = 0;
    int nNegativePrt = 0;
    int nNegativeMu = 0;

    int imax = 1;
    int jmax = 1;
    int kmax = 1;

    int imin = 1;
    int jmin = 1;
    int kmin = 1;

    int nDim = GetDim();

    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        int ist, ied, jst, jed, kst, ked;
        grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    RDouble dudx = gradPrimtiveVarFaceXofLES(il, jl, kl, IU, iSurface);
                    RDouble dudy = gradPrimtiveVarFaceYofLES(il, jl, kl, IU, iSurface);
                    RDouble dudz = gradPrimtiveVarFaceZofLES(il, jl, kl, IU, iSurface);

                    RDouble dvdx = gradPrimtiveVarFaceXofLES(il, jl, kl, IV, iSurface);
                    RDouble dvdy = gradPrimtiveVarFaceYofLES(il, jl, kl, IV, iSurface);
                    RDouble dvdz = gradPrimtiveVarFaceZofLES(il, jl, kl, IV, iSurface);

                    RDouble dwdx = gradPrimtiveVarFaceXofLES(il, jl, kl, IW, iSurface);
                    RDouble dwdy = gradPrimtiveVarFaceYofLES(il, jl, kl, IW, iSurface);
                    RDouble dwdz = gradPrimtiveVarFaceZofLES(il, jl, kl, IW, iSurface);

                    RDouble dtdx = gradtx(il, jl, kl, ITT, iSurface);
                    RDouble dtdy = gradty(il, jl, kl, ITT, iSurface);
                    RDouble dtdz = gradtz(il, jl, kl, ITT, iSurface);

                    RDouble s11 = dudx;
                    RDouble s22 = dvdy;
                    RDouble s33 = dwdz;
                    RDouble s12 = half * (dudy + dvdx);
                    RDouble s13 = half * (dudz + dwdx);
                    RDouble s23 = half * (dvdz + dwdy);

                    RDouble sijsij = s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23);
                    RDouble ss = sqrt(two * sijsij);
                    RDouble skk = s11 + s22 + s33;

                    RDouble rho = 0.5 * (q(i, j, k, IR) + q(il, jl, kl, IR));
                    RDouble mul = 0.5 * (visl(i, j, k) + visl(il, jl, kl));
                    RDouble dampingFunction = wallFunction(il, jl, kl, iSurface);

                    RDouble cr = 0.5 * (anisotropicConstant(il, jl, kl) + anisotropicConstant(i, j, k));

                    if (cr < -0.01)
                    {
                        //cr = -0.01;
                        nLowerBoundCr = nLowerBoundCr + 1;
                    }

                    /*
                    if (cr > 0.5)
                    {
                        cr = 0.5;
                        nUpperBoundCr = nUpperBoundCr + 1;
                    }
                    */

                    RDouble mut = cr * dampingFunction * dampingFunction * rho * ss * refReNumber;

                    RDouble ci = 0.5 * (isotropicConstant(i, j, k) + isotropicConstant(il, jl, kl));
                    //ci = MIN(MAX(ci, 0.0), 0.01);

                    RDouble kSGS = 2.0 * ci * dampingFunction * dampingFunction * rho * ss * ss * refReNumber;

                    oprt = 2.0 / (turbulentPrandtlNumber(il, jl, kl) + turbulentPrandtlNumber(i, j, k));
                    //kSGS = 0.0;

                    RDouble hfxSGS = -mut * specificHeatAtConstantPressure * oprt * dtdx;
                    RDouble hfySGS = -mut * specificHeatAtConstantPressure * oprt * dtdy;
                    RDouble hfzSGS = -mut * specificHeatAtConstantPressure * oprt * dtdz;

                    subgridHeatFlux(il, jl, kl, 1, iSurface) = hfxSGS;
                    subgridHeatFlux(il, jl, kl, 2, iSurface) = hfySGS;
                    subgridHeatFlux(il, jl, kl, 3, iSurface) = hfzSGS;

                    if (mut < -mul)
                    {
                        nNegativeMu = nNegativeMu + 1;
                    }

                    if (turbViscousCutType == 1)
                    {
                        mut = MAX(mut - mul, 0.0);
                    }
                    else if (turbViscousCutType == 2)
                    {
                        mut = MAX(mut, 0.0);
                    }
                    else if (turbViscousCutType == 3)
                    {
                        mut = MAX(mut, -mul);
                    }

                    if (prtCutType == 1)
                    {
                        oprt = MAX(oprt, 0.0);
                    }

                    mut = MIN(eddyViscosityLimit, mut);

                    RDouble txxSGS = -2.0 * mut * (s11 - 1.0/3.0 * skk) + 1.0/3.0 * kSGS;
                    RDouble tyySGS = -2.0 * mut * (s22 - 1.0/3.0 * skk) + 1.0/3.0 * kSGS;
                    RDouble tzzSGS = -2.0 * mut * (s33 - 1.0/3.0 * skk) + 1.0/3.0 * kSGS;
                    RDouble txySGS = -2.0 * mut * s12;
                    RDouble txzSGS = -2.0 * mut * s13;
                    RDouble tyzSGS = -2.0 * mut * s23;

                    subgridStress(il, jl, kl, 1, iSurface) = txxSGS;
                    subgridStress(il, jl, kl, 2, iSurface) = tyySGS;
                    subgridStress(il, jl, kl, 3, iSurface) = tzzSGS;
                    subgridStress(il, jl, kl, 4, iSurface) = txySGS;
                    subgridStress(il, jl, kl, 5, iSurface) = txzSGS;
                    subgridStress(il, jl, kl, 6, iSurface) = tyzSGS;     

                    if (vistmax < mut)
                    {
                        vistmax = mut;
                        imax = i;
                        jmax = j;
                        kmax = k;
                    }

                    if (vistmin > mut)
                    {
                        vistmin = mut;
                        imin = i;
                        jmin = j;
                        kmin = k;
                    }

                    if (cr < 0.0)
                    {
                        nNegativeCr = nNegativeCr + 1;
                    }

                    if (ci < 0.0)
                    {
                        nNegativeCi = nNegativeCi + 1;
                    }

                    if (oprt < 0.0)
                    {
                        nNegativePrt = nNegativePrt + 1;
                    }
                }
            }
        }

    }

    double vist_max = vistmax;
    double vist_min = vistmin;

    grid->UpdateData("vist_max", &vist_max, PHDOUBLE, 1);
    grid->UpdateData("vist_min", &vist_min, PHDOUBLE, 1);
    grid->UpdateData("nNegativeCr", &nNegativeCr, PHINT, 1);
    grid->UpdateData("nNegativeCi", &nNegativeCi, PHINT, 1);
    grid->UpdateData("nNegativePrt", &nNegativePrt, PHINT, 1);
    grid->UpdateData("nNegativeMu", &nNegativeMu, PHINT, 1);
    grid->UpdateData("nLowerBoundCr", &nLowerBoundCr, PHINT, 1);
    grid->UpdateData("nUpperBoundCr", &nUpperBoundCr, PHINT, 1);

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorVistmax)
        {
            cout << "vistmax = " << vistmax << " " << imax << " " << jmax << " " << kmax << "\n";
            cout << "vistmin = " << vistmin << " " << imin << " " << jmin << " " << kmin << "\n";
        }

        if (monitorNegativeConstant)
        {
            RDouble percentNegativeCr = nNegativeCr / static_cast<RDouble> (ni * nj * nk * nDim);
            RDouble percentNegativeMu = nNegativeMu / static_cast<RDouble> (ni * nj * nk * nDim);
            RDouble percentLowerBoundCr = nLowerBoundCr / static_cast<RDouble> (ni * nj * nk * nDim);
            RDouble percentUpperBoundCr = nUpperBoundCr / static_cast<RDouble> (ni * nj * nk * nDim);

            cout << "percentNegativeCr = " << percentNegativeCr << "\n";
            //cout << "nNegativeCi = " << nNegativeCi << "\n";
            //cout << "nNegativePrt = " << nNegativePrt << "\n";
            cout << "percentNegativeMu = " << percentNegativeMu << "\n";
            cout << "percentLowerBoundCr = " << percentLowerBoundCr << "\n";
            cout << "percentUpperBoundCr = " << percentUpperBoundCr << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif
}

void LESSolverStruct::DynamicSmagViscosityCompressible(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int monitorNegativeConstant = GlobalDataBase::GetIntParaFromDB("monitorNegativeConstant");

    RDouble3D &wallFunction = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));

    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));
    RDouble3D &anisotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));
    RDouble3D &isotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("isotropicConstant"));

    ComputeAnisotropicConstant(grid);

    ComputeIsotropicConstant(grid);

    ComputePrandtlNumber(grid);

    int nNegativeCr = 0;
    int nLowerBoundCr = 0;
    int nUpperBoundCr = 0;
    int nNegativeCi = 0;
    int nNegativePrt = 0;

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble dampingFunction = wallFunction(i, j, k);

                RDouble cr = anisotropicConstant(i, j, k);

                RDouble ci = isotropicConstant(i, j, k);

                RDouble prt = turbulentPrandtlNumber(i, j, k);

                if (cr < -0.01)
                {
                    //cr = -0.01;
                    nLowerBoundCr = nLowerBoundCr + 1;
                }

                if (cr > 0.5)
                {
                    //cr = 0.5;
                    nUpperBoundCr = nUpperBoundCr + 1;
                }

                if (cr < 0.0)
                {
                    nNegativeCr = nNegativeCr + 1;
                }

                if (ci < 0.0)
                {
                    nNegativeCi = nNegativeCi + 1;
                }

                if (prt < 0.0)
                {
                    nNegativePrt = nNegativePrt + 1;
                }

                anisotropicConstant(i, j, k) = cr * dampingFunction * dampingFunction;

                isotropicConstant(i, j, k) = ci * dampingFunction * dampingFunction;

                //turbulentPrandtlNumber(i, j, k) = prt;
            }
        }
    }

    grid->UpdateData("nNegativeCr", &nNegativeCr, PHINT, 1);
    grid->UpdateData("nNegativeCi", &nNegativeCi, PHINT, 1);
    grid->UpdateData("nNegativePrt", &nNegativePrt, PHINT, 1);
    grid->UpdateData("nLowerBoundCr", &nLowerBoundCr, PHINT, 1);
    grid->UpdateData("nUpperBoundCr", &nUpperBoundCr, PHINT, 1);

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorNegativeConstant)
        {
            RDouble percentNegativeCr = nNegativeCr / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));
            RDouble percentLowerBoundCr = nLowerBoundCr / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));
            RDouble percentUpperBoundCr = nUpperBoundCr / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));

            cout << "percentNegativeCr = " << percentNegativeCr << "\n";
            cout << "percentNegativeCi = " << nNegativeCi << "\n";
            cout << "percentNegativePrt = " << nNegativePrt << "\n";
            cout << "percentLowerBoundCr = " << percentLowerBoundCr << "\n";
            cout << "percentUpperBoundCr = " << percentUpperBoundCr << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif
}

void LESSolverStruct::DynamicSmagViscosityCompressibleFD(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &q = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble4D &gradUVWTCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &gradUVWTCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    Param_LESSolverStruct *parameters = GetControlParameters();

    RDouble refReNumber = parameters->GetRefReNumber();

    int monitorNegativeConstant = GlobalDataBase::GetIntParaFromDB("monitorNegativeConstant");

    RDouble3D &visl = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble3D &wallFunction = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));
    RDouble3D &viscousTurbulence = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble3D &subgridScaleEnergy = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("subgridScaleEnergy"));

    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));
    RDouble3D &anisotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));
    RDouble3D &isotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("isotropicConstant"));

    ComputeAnisotropicConstant(grid);

    ComputeIsotropicConstant(grid);

    ComputePrandtlNumber(grid);

    int nNegativeCr = 0;
    int nLowerBoundCr = 0;
    int nUpperBoundCr = 0;
    int nNegativeCi = 0;
    int nNegativePrt = 0;

    int turbViscousCutType = parameters->GetTurbViscousCutType();

    RDouble viscousTurbulenceMaximum = 0.0;
    int imax = 1;
    int jmax = 1;
    int kmax = 1;

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble cr = anisotropicConstant(i, j, k);
                RDouble ci = isotropicConstant(i, j, k);
                RDouble prt = turbulentPrandtlNumber(i, j, k);

                RDouble dampingFunction = wallFunction(i, j, k);
                //RDouble dampingFunction = 1.0;

                RDouble rho = ABS(q(i, j, k, IR)) + SMALL;

                RDouble dudx = gradUVWTCellCenterX(i, j, k, 0);
                RDouble dudy = gradUVWTCellCenterY(i, j, k, 0);
                RDouble dudz = gradUVWTCellCenterZ(i, j, k, 0);

                RDouble dvdx = gradUVWTCellCenterX(i, j, k, 1);
                RDouble dvdy = gradUVWTCellCenterY(i, j, k, 1);
                RDouble dvdz = gradUVWTCellCenterZ(i, j, k, 1);

                RDouble dwdx = gradUVWTCellCenterX(i, j, k, 2);
                RDouble dwdy = gradUVWTCellCenterY(i, j, k, 2);
                RDouble dwdz = gradUVWTCellCenterZ(i, j, k, 2);

                RDouble s11 = dudx; 
                RDouble s22 = dvdy; 
                RDouble s33 = dwdz;
                RDouble s12 = half * (dudy + dvdx); 
                RDouble s13 = half * (dudz + dwdx); 
                RDouble s23 = half * (dvdz + dwdy);

                RDouble sij2 = two * (s11*s11 + s22*s22 + s33*s33 + two * (s12*s12 + s13*s13 + s23*s23));  //! modulus of S.
                RDouble Strain = sqrt(sij2);

                RDouble mut = dampingFunction * dampingFunction * cr * rho * Strain * refReNumber;
                RDouble kSGS = 2.0 * dampingFunction * dampingFunction * ci * rho * Strain * Strain * refReNumber;

                RDouble mul = visl(i, j, k);
                if (turbViscousCutType == 1)
                {
                    mut = MAX(mut - mul, 0.0);
                    kSGS = 0.0;
                    prt = 0.9;
                }
                else if (turbViscousCutType == 2)
                {
                    mut = MAX(mut, 0.0);
                    kSGS = 0.0;
                    prt = 0.9;
                }
                else if (turbViscousCutType == 3)
                {
                    mut = MAX(mut, -mul);
                    kSGS = 0.0;
                    prt = 0.9;
                }

                viscousTurbulence(i, j, k) = mut;
                subgridScaleEnergy(i, j, k) = kSGS;
                turbulentPrandtlNumber(i, j, k) = prt;

                if (cr < -0.01)
                {
                    //cr = -0.01;
                    nLowerBoundCr = nLowerBoundCr + 1;
                }

                if (cr > 0.5)
                {
                    //cr = 0.5;
                    nUpperBoundCr = nUpperBoundCr + 1;
                }

                if (cr < 0.0)
                {
                    nNegativeCr = nNegativeCr + 1;
                }

                if (ci < 0.0)
                {
                    nNegativeCi = nNegativeCi + 1;
                }

                if (prt < 0.0)
                {
                    nNegativePrt = nNegativePrt + 1;
                }

                if (viscousTurbulenceMaximum < mut)
                {
                    viscousTurbulenceMaximum = mut;
                    imax = i;
                    jmax = j;
                    kmax = k;
                }
            }
        }
    }

    grid->UpdateData("nNegativeCr", &nNegativeCr, PHINT, 1);
    grid->UpdateData("nNegativeCi", &nNegativeCi, PHINT, 1);
    grid->UpdateData("nNegativePrt", &nNegativePrt, PHINT, 1);
    grid->UpdateData("nLowerBoundCr", &nLowerBoundCr, PHINT, 1);
    grid->UpdateData("nUpperBoundCr", &nUpperBoundCr, PHINT, 1);

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorNegativeConstant)
        {
            RDouble percentNegativeCr = nNegativeCr / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));
            RDouble percentLowerBoundCr = nLowerBoundCr / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));
            RDouble percentUpperBoundCr = nUpperBoundCr / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));

            cout << "percentNegativeCr = " << percentNegativeCr << "\n";
            cout << "percentNegativeCi = " << nNegativeCi << "\n";
            cout << "percentNegativePrt = " << nNegativePrt << "\n";
            cout << "percentLowerBoundCr = " << percentLowerBoundCr << "\n";
            cout << "percentUpperBoundCr = " << percentUpperBoundCr << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif
}

void LESSolverStruct::DynamicSmagViscosityIncompressibleOld(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_LESSolverStruct *parameters = GetControlParameters();

    RDouble refReNumber = parameters->GetRefReNumber();

    int turbViscousCutType = parameters->GetTurbViscousCutType();

    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    int monitorVistmax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");

    int monitorNegativeConstant = GlobalDataBase::GetIntParaFromDB("monitorNegativeConstant");

    int prtCutType = GlobalDataBase::GetIntParaFromDB("prtCutType");

    RDouble oprt = parameters->GetoPrandtlTurbulence();

    RDouble refGama = parameters->GetRefGama();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble specificHeatAtConstantPressure = 1.0 / ((refGama - 1.0) * refMachNumber * refMachNumber);

    RDouble4D &wallFunction = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("wallFunction"));

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &visl = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble5D &gradPrimtiveVarFaceXofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceXofLES"));
    RDouble5D &gradPrimtiveVarFaceYofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceYofLES"));
    RDouble5D &gradPrimtiveVarFaceZofLES = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradPrimtiveVarFaceZofLES"));

    RDouble5D &gradtx = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceXofLES"));
    RDouble5D &gradty = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceYofLES"));
    RDouble5D &gradtz = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("gradTemperatureFaceZofLES"));

    RDouble3D &anisotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));
    RDouble3D &isotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("isotropicConstant"));

    RDouble5D &subgridStress = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("subgridStress"));
    RDouble5D &subgridHeatFlux = * reinterpret_cast<RDouble5D *> (grid->GetDataPtr("subgridHeatFlux"));

    ComputeAnisotropicConstant(grid);

    //ComputePrandtlNumber(grid);

    double vistmax  = 0.0;
    double vistmin  = 1.0e30;

    int nNegativeCr = 0;
    int nLowerBoundCr = 0;
    int nUpperBoundCr = 0;
    int nNegativeCi = 0;
    int nNegativePrt = 0;
    int nNegativeMu = 0;

    int imax = 1;
    int jmax = 1;
    int kmax = 1;

    int imin = 1;
    int jmin = 1;
    int kmin = 1;

    int nDim = GetDim();

    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        int ist, ied, jst, jed, kst, ked;
        grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    RDouble dudx = gradPrimtiveVarFaceXofLES(il, jl, kl, IU, iSurface);
                    RDouble dudy = gradPrimtiveVarFaceYofLES(il, jl, kl, IU, iSurface);
                    RDouble dudz = gradPrimtiveVarFaceZofLES(il, jl, kl, IU, iSurface);

                    RDouble dvdx = gradPrimtiveVarFaceXofLES(il, jl, kl, IV, iSurface);
                    RDouble dvdy = gradPrimtiveVarFaceYofLES(il, jl, kl, IV, iSurface);
                    RDouble dvdz = gradPrimtiveVarFaceZofLES(il, jl, kl, IV, iSurface);

                    RDouble dwdx = gradPrimtiveVarFaceXofLES(il, jl, kl, IW, iSurface);
                    RDouble dwdy = gradPrimtiveVarFaceYofLES(il, jl, kl, IW, iSurface);
                    RDouble dwdz = gradPrimtiveVarFaceZofLES(il, jl, kl, IW, iSurface);

                    RDouble dtdx = gradtx(il, jl, kl, ITT, iSurface);
                    RDouble dtdy = gradty(il, jl, kl, ITT, iSurface);
                    RDouble dtdz = gradtz(il, jl, kl, ITT, iSurface);

                    RDouble s11 = dudx;
                    RDouble s22 = dvdy;
                    RDouble s33 = dwdz;
                    RDouble s12 = half * (dudy + dvdx);
                    RDouble s13 = half * (dudz + dwdx);
                    RDouble s23 = half * (dvdz + dwdy);

                    RDouble sijsij = s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23);
                    RDouble ss = sqrt(two * sijsij);

                    RDouble rho = 0.5 * (q(i, j, k, IR) + q(il, jl, kl, IR));
                    RDouble mul = 0.5 * (visl(i, j, k) + visl(il, jl, kl));
                    RDouble dampingFunction = wallFunction(il, jl, kl, iSurface);

                    RDouble cr = 0.5 * (anisotropicConstant(il, jl, kl) + anisotropicConstant(i, j, k));

                    if (cr < -0.01)
                    {
                        //cr = -0.01;
                        nLowerBoundCr = nLowerBoundCr + 1;
                    }

                    if (cr > 0.5)
                    {
                        //cr = 0.5;
                        nUpperBoundCr = nUpperBoundCr + 1;
                    }                   

                    RDouble mut = cr * dampingFunction * dampingFunction * rho * ss * refReNumber;

                    RDouble ci = 0.5 * (isotropicConstant(i, j, k) + isotropicConstant(il, jl, kl));
                    //ci = MIN(MAX(ci, 0.0), 0.01);

                    RDouble kSGS = 2.0 * ci * dampingFunction * dampingFunction * rho * ss * ss * refReNumber;

                    //oprt = 2.0 / (turbulentPrandtlNumber(il, jl, kl) + turbulentPrandtlNumber(i, j, k));
                    oprt = 1.111111;

                    RDouble hfxSGS = -mut * specificHeatAtConstantPressure * oprt * dtdx;
                    RDouble hfySGS = -mut * specificHeatAtConstantPressure * oprt * dtdy;
                    RDouble hfzSGS = -mut * specificHeatAtConstantPressure * oprt * dtdz;

                    subgridHeatFlux(il, jl, kl, 1, iSurface) = hfxSGS;
                    subgridHeatFlux(il, jl, kl, 2, iSurface) = hfySGS;
                    subgridHeatFlux(il, jl, kl, 3, iSurface) = hfzSGS;

                    if (mut < -mul)
                    {
                        nNegativeMu = nNegativeMu + 1;
                    }

                    if (turbViscousCutType == 1)
                    {
                        mut = MAX(mut - mul, 0.0);
                    }
                    else if (turbViscousCutType == 2)
                    {
                        mut = MAX(mut, 0.0);
                    }
                    else if (turbViscousCutType == 3)
                    {
                        mut = MAX(mut, -mul);
                    }

                    if (prtCutType == 1)
                    {
                        oprt = MAX(oprt, 0.0);
                    }

                    mut = MIN(eddyViscosityLimit, mut);

                    RDouble txxSGS = -2.0 * mut * s11 + 1.0/3.0 * kSGS;
                    RDouble tyySGS = -2.0 * mut * s22 + 1.0/3.0 * kSGS;
                    RDouble tzzSGS = -2.0 * mut * s33 + 1.0/3.0 * kSGS;
                    RDouble txySGS = -2.0 * mut * s12;
                    RDouble txzSGS = -2.0 * mut * s13;
                    RDouble tyzSGS = -2.0 * mut * s23;

                    subgridStress(il, jl, kl, 1, iSurface) = txxSGS;
                    subgridStress(il, jl, kl, 2, iSurface) = tyySGS;
                    subgridStress(il, jl, kl, 3, iSurface) = tzzSGS;
                    subgridStress(il, jl, kl, 4, iSurface) = txySGS;
                    subgridStress(il, jl, kl, 5, iSurface) = txzSGS;
                    subgridStress(il, jl, kl, 6, iSurface) = tyzSGS;     

                    if (vistmax < mut)
                    {
                        vistmax = mut;
                        imax = i;
                        jmax = j;
                        kmax = k;
                    }

                    if (vistmin > mut)
                    {
                        vistmin = mut;
                        imin = i;
                        jmin = j;
                        kmin = k;
                    }

                    if (cr < 0.0)
                    {
                        nNegativeCr = nNegativeCr + 1;
                    }

                    if (ci < 0.0)
                    {
                        nNegativeCi = nNegativeCi + 1;
                    }

                    if (oprt < 0.0)
                    {
                        nNegativePrt = nNegativePrt + 1;
                    }
                }
            }
        }

    }

    double vist_max = vistmax;
    double vist_min = vistmin;

    grid->UpdateData("vist_max", &vist_max, PHDOUBLE, 1);
    grid->UpdateData("vist_min", &vist_min, PHDOUBLE, 1);
    grid->UpdateData("nNegativeCr", &nNegativeCr, PHINT, 1);
    grid->UpdateData("nNegativeCi", &nNegativeCi, PHINT, 1);
    grid->UpdateData("nNegativePrt", &nNegativePrt, PHINT, 1);
    grid->UpdateData("nNegativeMu", &nNegativeMu, PHINT, 1);
    grid->UpdateData("nLowerBoundCr", &nLowerBoundCr, PHINT, 1);
    grid->UpdateData("nUpperBoundCr", &nUpperBoundCr, PHINT, 1);

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorVistmax)
        {
            cout << "vistmax = " << vistmax << " " << imax << " " << jmax << " " << kmax << "\n";
            cout << "vistmin = " << vistmin << " " << imin << " " << jmin << " " << kmin << "\n";
        }

        if (monitorNegativeConstant)
        {
            RDouble percentNegativeCr = nNegativeCr / static_cast<RDouble> (ni * nj * nk * nDim);
            RDouble percentNegativeMu = nNegativeMu / static_cast<RDouble> (ni * nj * nk * nDim);
            RDouble percentLowerBoundCr = nLowerBoundCr / static_cast<RDouble> (ni * nj * nk * nDim);
            RDouble percentUpperBoundCr = nUpperBoundCr / static_cast<RDouble> (ni * nj * nk * nDim);

            cout << "percentNegativeCr = " << percentNegativeCr << "\n";
            //cout << "nNegativeCi = " << nNegativeCi << "\n";
            //cout << "nNegativePrt = " << nNegativePrt << "\n";
            cout << "percentNegativeMu = " << percentNegativeMu << "\n";
            cout << "percentLowerBoundCr = " << percentLowerBoundCr << "\n";
            cout << "percentUpperBoundCr = " << percentUpperBoundCr << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif
}

void LESSolverStruct::DynamicSmagViscosityIncompressible(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_LESSolverStruct *parameters = GetControlParameters();

    //RDouble refReNumber  = parameters->GetRefReNumber();

    int monitorNegativeConstant = GlobalDataBase::GetIntParaFromDB("monitorNegativeConstant");

    RDouble3D &wallFunction = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));

    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));
    RDouble3D &anisotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));
    RDouble3D &isotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("isotropicConstant"));

    ComputeAnisotropicConstant(grid);

    //ComputePrandtlNumber(grid);

    int nNegativeCr = 0;
    int nLowerBoundCr = 0;
    int nUpperBoundCr = 0;

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble dampingFunction = wallFunction(i, j, k);

                RDouble cr = anisotropicConstant(i, j, k);

                if (cr < -0.01)
                {
                    //cr = -0.01;
                    nLowerBoundCr = nLowerBoundCr + 1;
                }

                if (cr > 0.5)
                {
                    //cr = 0.5;
                    nUpperBoundCr = nUpperBoundCr + 1;
                }                    

                if (cr < 0.0)
                {
                    nNegativeCr = nNegativeCr + 1;
                }

                anisotropicConstant(i, j, k) = cr * dampingFunction * dampingFunction;
            }
        }
    }

    RDouble prandtlTurbulence = parameters->GetPrandtlTurbulence();

    turbulentPrandtlNumber = prandtlTurbulence;

    isotropicConstant = 0.0;

    grid->UpdateData("nNegativeCr", &nNegativeCr, PHINT, 1);
    grid->UpdateData("nLowerBoundCr", &nLowerBoundCr, PHINT, 1);
    grid->UpdateData("nUpperBoundCr", &nUpperBoundCr, PHINT, 1);

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorNegativeConstant)
        {
            RDouble percentNegativeCr = nNegativeCr / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));
            RDouble percentLowerBoundCr = nLowerBoundCr / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));
            RDouble percentUpperBoundCr = nUpperBoundCr / static_cast<RDouble> ((ni + 1) * (nj + 1) * (nk + 1));

            cout << "percentNegativeCr = " << percentNegativeCr << "\n";
            cout << "percentLowerBoundCr = " << percentLowerBoundCr << "\n";
            cout << "percentUpperBoundCr = " << percentUpperBoundCr << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif
}

void LESSolverStruct::ComputeTestFilteredQ(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    testFilteredQ(I, J, K, Range(IR, IW)) = q(I, J, K, Range(IR, IW));

    //! Get conservative variables for test filtering
    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    for (int m = IU; m <= IW; ++m)
    {
        for (int k = kst; k <= ked; ++k)
        {
            for (int j = jst; j <= jed; ++j)
            {
                for (int i = ist; i <= ied; ++i)
                {
                    RDouble rho = q(i, j, k, IR);
                    RDouble primitiveVariable = q(i, j, k, m);
                    testFilteredQ(i, j, k, m) = rho * primitiveVariable;
                }
            }
        }
    }

    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);
    //! filter Q in I direction
    if (filterDirection[0] == 1) 
    {
        for (int n = 0; n <= 3; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredQ(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredQ(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredQ(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    //! filter Q in J direction
    if (filterDirection[1] == 1)
    {
        for (int n = 0; n <= 3; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredQ(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredQ(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredQ(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    //! filter Q in K direction
    if (filterDirection[2] == 1) 
    {
        for (int n = 0; n <= 3; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredQ(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredQ(i, j, K, n);
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredQ(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    //! Recover primitive variables after test filtering
    for (int m = IU; m <= IW; ++m)
    {
        for (int k = kst; k <= ked; ++k)
        {
            for (int j = jst; j <= jed; ++j)
            {
                for (int i = ist; i <= ied; ++i)
                {
                    RDouble rho = testFilteredQ(i, j, k, IR);
                    RDouble conservativeVariable = testFilteredQ(i, j, k, m);
                    testFilteredQ(i, j, k, m) = conservativeVariable / rho;
                }
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeTestFilteredQFD(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);    

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_FD"));

    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    testFilteredQ(I, J, K, Range(IR, IW)) = q(I, J, K, Range(IR, IW));

    // !get conservative variables for test filtering
    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    for (int m = IU; m <= IW; ++m)
    {
        for (int k = kst; k <= ked; ++k)
        {
            for (int j = jst; j <= jed; ++j)
            {
                for (int i = ist; i <= ied; ++i)
                {
                    RDouble rho = q(i, j, k, IR);
                    RDouble primitiveVariable = q(i, j, k, m);
                    testFilteredQ(i, j, k, m) = rho * primitiveVariable;
                }
            }
        }
    }

    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);
    //! filter Q in I direction
    if (filterDirection[0] == 1) 
    {
        for (int n = 0; n <= 3; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredQ(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredQ(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredQ(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    //! filter Q in J direction
    if (filterDirection[1] == 1)
    {
        for (int n = 0; n <= 3; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredQ(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredQ(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredQ(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    //! filter Q in K direction
    if (filterDirection[2] == 1) 
    {
        for (int n = 0; n <= 3; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredQ(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredQ(i, j, K, n);
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredQ(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    //! Recover primitive variables after test filtering
    for (int m = IU; m <= IW; ++m)
    {
        for (int k = kst; k <= ked; ++k)
        {
            for (int j = jst; j <= jed; ++j)
            {
                for (int i = ist; i <= ied; ++i)
                {
                    RDouble rho = testFilteredQ(i, j, k, IR);
                    RDouble conservativeVariable = testFilteredQ(i, j, k, m);
                    testFilteredQ(i, j, k, m) = conservativeVariable / rho;
                }
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

/*
void LESSolverStruct::ComputeTestFilteredDensity(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble3D &testFilteredDensity = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("testFilteredDensity"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    testFilteredDensity(I, J, K) = q(I, J, K, IR);

    //! filter density in I direction
    if (filterDirection[0] == 1) 
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                (*testFilteredVariableX)(I) = testFilteredDensity(I, j, k);
                (*filteredVariableX)(I) = testFilteredDensity(I, j, k);
                TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                testFilteredDensity(I, j, k) = (*testFilteredVariableX)(I);
            }
        }
    }

    //! filter density in J direction
    if (filterDirection[1] == 1)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                (*testFilteredVariableY)(J) = testFilteredDensity(i, J, k);
                (*filteredVariableY)(J) = testFilteredDensity(i, J, k);
                TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                testFilteredDensity(i, J, k) = (*testFilteredVariableY)(J);
            }
        }
    }

    //! filter density in K direction
    if (filterDirection[2] == 1)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                (*testFilteredVariableZ)(K) = testFilteredDensity(i, j, K);
                (*filteredVariableZ)(K) = testFilteredDensity(i, j, K);
                TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                testFilteredDensity(i, j, K) = (*testFilteredVariableZ)(K);
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}
*/

void LESSolverStruct::ComputeTestFilteredTemperature(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &temperatures = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("t"));
    RDouble4D &q = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble3D &testFilteredTemperature = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("testFilteredTemperature"));
    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);
  
    //! Get conservative variables for test filtering
    testFilteredTemperature(I, J, K) = temperatures(I, J, K, ITT) * q(I, J, K, IR);

    //! filter temperature in I direction
    if (filterDirection[0] == 1)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                (*testFilteredVariableX)(I) = testFilteredTemperature(I, j, k);
                (*filteredVariableX)(I) = testFilteredTemperature(I, j, k);
                TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                testFilteredTemperature(I, j, k) = (*testFilteredVariableX)(I);
            }
        }
    }

    //! filter temperature in J direction
    if (filterDirection[1] == 1)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                (*testFilteredVariableY)(J) = testFilteredTemperature(i, J, k);
                (*filteredVariableY)(J) = testFilteredTemperature(i, J, k);
                TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                testFilteredTemperature(i, J, k) = (*testFilteredVariableY)(J);
            }
        }
    }

    //! filter temperature in K direction
    if (filterDirection[2] == 1)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                (*testFilteredVariableZ)(K) = testFilteredTemperature(i, j, K);
                (*filteredVariableZ)(K) = testFilteredTemperature(i, j, K);
                TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                testFilteredTemperature(i, j, K) = (*testFilteredVariableZ)(K);
            }
        }
    }

    //! Recover primitive variables after test filtering
    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                RDouble rho = testFilteredQ(i, j, k, IR);
                RDouble conservativeVariable = testFilteredQ(i, j, k, IR);
                testFilteredTemperature(i, j, k) = conservativeVariable / rho;
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeTestFilteredTemperatureFD(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &temperatures = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("t_FD"));
    RDouble4D &q = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_FD"));

    RDouble3D &testFilteredTemperature = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("testFilteredTemperature"));
    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);
  
    //! Get conservative variables for test filtering
    testFilteredTemperature(I, J, K) = temperatures(I, J, K, ITT) * q(I, J, K, IR);

    //! filter temperature in I direction
    if (filterDirection[0] == 1)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                (*testFilteredVariableX)(I) = testFilteredTemperature(I, j, k);
                (*filteredVariableX)(I) = testFilteredTemperature(I, j, k);
                TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                testFilteredTemperature(I, j, k) = (*testFilteredVariableX)(I);
            }
        }
    }

    //! filter temperature in J direction
    if (filterDirection[1] == 1)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                (*testFilteredVariableY)(J) = testFilteredTemperature(i, J, k);
                (*filteredVariableY)(J) = testFilteredTemperature(i, J, k);
                TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                testFilteredTemperature(i, J, k) = (*testFilteredVariableY)(J);
            }
        }
    }

    //! filter temperature in K direction
    if (filterDirection[2] == 1)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                (*testFilteredVariableZ)(K) = testFilteredTemperature(i, j, K);
                (*filteredVariableZ)(K) = testFilteredTemperature(i, j, K);
                TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                testFilteredTemperature(i, j, K) = (*testFilteredVariableZ)(K);
            }
        }
    }

    //! Recover primitive variables after test filtering
    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                RDouble rho = testFilteredQ(i, j, k, IR);
                RDouble conservativeVariable = testFilteredQ(i, j, k, IR);
                testFilteredTemperature(i, j, k) = conservativeVariable / rho;
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeStrainRateTensor(Grid *gridIn)
{
    Param_LESSolverStruct *parameters = GetControlParameters();
    int subgridScaleModel = parameters->GetSubgridScaleModel();
    if (subgridScaleModel == 6)
    {
        return;
    }

    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);

    RDouble4D &strainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));

    RDouble4D &gradUVWTCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &gradUVWTCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble dudx = gradUVWTCellCenterX(i, j, k, 0);
                RDouble dudy = gradUVWTCellCenterY(i, j, k, 0);
                RDouble dudz = gradUVWTCellCenterZ(i, j, k, 0);

                RDouble dvdx = gradUVWTCellCenterX(i, j, k, 1);
                RDouble dvdy = gradUVWTCellCenterY(i, j, k, 1);
                RDouble dvdz = gradUVWTCellCenterZ(i, j, k, 1);

                RDouble dwdx = gradUVWTCellCenterX(i, j, k, 2);
                RDouble dwdy = gradUVWTCellCenterY(i, j, k, 2);
                RDouble dwdz = gradUVWTCellCenterZ(i, j, k, 2);

                RDouble s11 = dudx; 
                RDouble s22 = dvdy; 
                RDouble s33 = dwdz;
                RDouble s12 = half * (dudy + dvdx); 
                RDouble s13 = half * (dudz + dwdx); 
                RDouble s23 = half * (dvdz + dwdy);

                RDouble s = 2.0 * (s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23));

                RDouble omega12 = 0.5 * (dudy - dvdx);
                RDouble omega13 = 0.5 * (dudz - dwdx);
                RDouble omega23 = 0.5 * (dvdz - dwdy);

                strainRateTensor(i, j, k, 0) = sqrt(s);

                strainRateTensor(i, j, k, 1) = s11;
                strainRateTensor(i, j, k, 2) = s22;
                strainRateTensor(i, j, k, 3) = s33;
                strainRateTensor(i, j, k, 4) = s12;
                strainRateTensor(i, j, k, 5) = s13;
                strainRateTensor(i, j, k, 6) = s23;

                strainRateTensor(i, j, k, 7) = omega12;
                strainRateTensor(i, j, k, 8) = omega13;
                strainRateTensor(i, j, k, 9) = omega23;

            }
        }
    }
}

void LESSolverStruct::ComputeGradT(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &gradUVWTCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &gradUVWTCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    RDouble4D &gradT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradT"));

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    gradT(I, J, K, 0) = gradUVWTCellCenterX(I, J, K, 3);
    gradT(I, J, K, 1) = gradUVWTCellCenterY(I, J, K, 3);
    gradT(I, J, K, 2) = gradUVWTCellCenterZ(I, J, K, 3);
}

void LESSolverStruct::ComputeGradTFD(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &gradUVWTCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &gradUVWTCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    RDouble4D &gradT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradT"));

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    gradT(I, J, K, 0) = gradUVWTCellCenterX(I, J, K, 3);
    gradT(I, J, K, 1) = gradUVWTCellCenterY(I, J, K, 3);
    gradT(I, J, K, 2) = gradUVWTCellCenterZ(I, J, K, 3);
}

void LESSolverStruct::ComputeTestFilteredStrainRateTensor(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Range N(0, 6);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &strainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));
    RDouble4D &testFilteredStrainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredStrainRateTensor"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    testFilteredStrainRateTensor(I, J, K, N) = strainRateTensor(I, J, K, N);

    //! filter Sij in I direction
    if (filterDirection[0] == 1)
    {
        for (int n = 1; n <= 6; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredStrainRateTensor(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredStrainRateTensor(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredStrainRateTensor(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    //! filter Sij in J direction
    if (filterDirection[1] == 1)
    {
        for (int n = 1; n <= 6; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredStrainRateTensor(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredStrainRateTensor(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredStrainRateTensor(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    //! filter Sij in K direction
    if (filterDirection[2] == 1)
    {
        for (int n = 1; n <= 6; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredStrainRateTensor(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredStrainRateTensor(i, j, K, n);
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredStrainRateTensor(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    //! Compute |S|
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble tFS11 = testFilteredStrainRateTensor(i, j, k, 1);
                RDouble tFS22 = testFilteredStrainRateTensor(i, j, k, 2);
                RDouble tFS33 = testFilteredStrainRateTensor(i, j, k, 3);
                RDouble tFS12 = testFilteredStrainRateTensor(i, j, k, 4);
                RDouble tFS13 = testFilteredStrainRateTensor(i, j, k, 5);
                RDouble tFS23 = testFilteredStrainRateTensor(i, j, k, 6);

                RDouble tFS = 2.0 * (tFS11 * tFS11 + tFS22 * tFS22 + tFS33 * tFS33 + 2.0 * (tFS12 * tFS12 + tFS13 * tFS13 + tFS23 * tFS23));
                testFilteredStrainRateTensor(i, j, k, 0) = sqrt(tFS);
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeTestFilteredStrainRateTensorFD(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Range N(0, 6);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    RDouble4D &testFilteredStrainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredStrainRateTensor"));

    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));

    RDouble3D gradux(I, J, K, fortranArray);
    RDouble3D graduy(I, J, K, fortranArray);
    RDouble3D graduz(I, J, K, fortranArray);

    RDouble3D gradvx(I, J, K, fortranArray);
    RDouble3D gradvy(I, J, K, fortranArray);
    RDouble3D gradvz(I, J, K, fortranArray);

    RDouble3D gradwx(I, J, K, fortranArray);
    RDouble3D gradwy(I, J, K, fortranArray);
    RDouble3D gradwz(I, J, K, fortranArray);

    GradCenter(grid, testFilteredQ, gradux, graduy, graduz, IU);
    GradCenter(grid, testFilteredQ, gradvx, gradvy, gradvz, IV);
    GradCenter(grid, testFilteredQ, gradwx, gradwy, gradwz, IW);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble dudx = gradux(i, j, k);
                RDouble dudy = graduy(i, j, k);
                RDouble dudz = graduz(i, j, k);

                RDouble dvdx = gradvx(i, j, k);
                RDouble dvdy = gradvy(i, j, k);
                RDouble dvdz = gradvz(i, j, k);

                RDouble dwdx = gradwx(i, j, k);
                RDouble dwdy = gradwy(i, j, k);
                RDouble dwdz = gradwz(i, j, k);

                RDouble s11 = dudx; 
                RDouble s22 = dvdy; 
                RDouble s33 = dwdz;
                RDouble s12 = half * (dudy + dvdx); 
                RDouble s13 = half * (dudz + dwdx); 
                RDouble s23 = half * (dvdz + dwdy);

                RDouble sij2 = two * (s11*s11 + s22*s22 + s33*s33 + two * (s12*s12 + s13*s13 + s23*s23));  //! modulus of S.
                RDouble Strain = sqrt(sij2);

                testFilteredStrainRateTensor(i, j, k, 0) = Strain;
                testFilteredStrainRateTensor(i, j, k, 1) = s11;
                testFilteredStrainRateTensor(i, j, k, 2) = s22;
                testFilteredStrainRateTensor(i, j, k, 3) = s33;
                testFilteredStrainRateTensor(i, j, k, 4) = s12;
                testFilteredStrainRateTensor(i, j, k, 5) = s13;
                testFilteredStrainRateTensor(i, j, k, 6) = s23;
            }
        }
    }
}

void LESSolverStruct::ComputeTestFilteredRhoUiUj(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();    
    int *filterDirection = parameters->GetFilterDirection();
    int subgridScaleModel  = parameters->GetSubgridScaleModel();
    
    RDouble4D &q = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    
    RDouble4D &testFilteredRhoUiUj = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoUiUj"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    if (subgridScaleModel == 2)
    {
        testFilteredRhoUiUj(I, J, K, 1) = q(I, J, K, IR) * q(I, J, K, IU) * q(I, J, K, IU);
        testFilteredRhoUiUj(I, J, K, 2) = q(I, J, K, IR) * q(I, J, K, IV) * q(I, J, K, IV);
        testFilteredRhoUiUj(I, J, K, 3) = q(I, J, K, IR) * q(I, J, K, IW) * q(I, J, K, IW);
        testFilteredRhoUiUj(I, J, K, 4) = q(I, J, K, IR) * q(I, J, K, IU) * q(I, J, K, IV);
        testFilteredRhoUiUj(I, J, K, 5) = q(I, J, K, IR) * q(I, J, K, IU) * q(I, J, K, IW);
        testFilteredRhoUiUj(I, J, K, 6) = q(I, J, K, IR) * q(I, J, K, IV) * q(I, J, K, IW);
    }
    else if (subgridScaleModel == 4)
    {
        testFilteredRhoUiUj(I, J, K, 1) = q(I, J, K, IU) * q(I, J, K, IU);
        testFilteredRhoUiUj(I, J, K, 2) = q(I, J, K, IV) * q(I, J, K, IV);
        testFilteredRhoUiUj(I, J, K, 3) = q(I, J, K, IW) * q(I, J, K, IW);
        testFilteredRhoUiUj(I, J, K, 4) = q(I, J, K, IU) * q(I, J, K, IV);
        testFilteredRhoUiUj(I, J, K, 5) = q(I, J, K, IU) * q(I, J, K, IW);
        testFilteredRhoUiUj(I, J, K, 6) = q(I, J, K, IV) * q(I, J, K, IW);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("subgridScaleModel", subgridScaleModel);
    }

    //! filter RhoUiUj in I direction
    if (filterDirection[0] == 1)
    {
        for (int n = 1; n <= 6; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredRhoUiUj(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredRhoUiUj(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredRhoUiUj(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    //! filter RhoUiUj in J direction
    if (filterDirection[1] == 1)
    {
        for (int n = 1; n <= 6; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredRhoUiUj(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredRhoUiUj(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredRhoUiUj(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    //! filter RhoUiUj in K direction
    if (filterDirection[2] == 1)
    {
        for (int n = 1; n <= 6; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredRhoUiUj(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredRhoUiUj(i, j, K, n);                
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredRhoUiUj(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeTestFilteredRhoUiUjFD(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();    
    int *filterDirection = parameters->GetFilterDirection();
    int subgridScaleModel  = parameters->GetSubgridScaleModel();
    
    RDouble4D &q = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_FD"));
    
    RDouble4D &testFilteredRhoUiUj = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoUiUj"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    if (subgridScaleModel == 2)
    {
        testFilteredRhoUiUj(I, J, K, 1) = q(I, J, K, IR) * q(I, J, K, IU) * q(I, J, K, IU);
        testFilteredRhoUiUj(I, J, K, 2) = q(I, J, K, IR) * q(I, J, K, IV) * q(I, J, K, IV);
        testFilteredRhoUiUj(I, J, K, 3) = q(I, J, K, IR) * q(I, J, K, IW) * q(I, J, K, IW);
        testFilteredRhoUiUj(I, J, K, 4) = q(I, J, K, IR) * q(I, J, K, IU) * q(I, J, K, IV);
        testFilteredRhoUiUj(I, J, K, 5) = q(I, J, K, IR) * q(I, J, K, IU) * q(I, J, K, IW);
        testFilteredRhoUiUj(I, J, K, 6) = q(I, J, K, IR) * q(I, J, K, IV) * q(I, J, K, IW);
    }
    else if (subgridScaleModel == 4)
    {
        testFilteredRhoUiUj(I, J, K, 1) = q(I, J, K, IU) * q(I, J, K, IU);
        testFilteredRhoUiUj(I, J, K, 2) = q(I, J, K, IV) * q(I, J, K, IV);
        testFilteredRhoUiUj(I, J, K, 3) = q(I, J, K, IW) * q(I, J, K, IW);
        testFilteredRhoUiUj(I, J, K, 4) = q(I, J, K, IU) * q(I, J, K, IV);
        testFilteredRhoUiUj(I, J, K, 5) = q(I, J, K, IU) * q(I, J, K, IW);
        testFilteredRhoUiUj(I, J, K, 6) = q(I, J, K, IV) * q(I, J, K, IW);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("subgridScaleModel", subgridScaleModel);
    }

    //! filter RhoUiUj in I direction
    if (filterDirection[0] == 1)
    {
        for (int n = 1; n <= 6; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredRhoUiUj(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredRhoUiUj(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredRhoUiUj(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    //! filter RhoUiUj in J direction
    if (filterDirection[1] == 1)
    {
        for (int n = 1; n <= 6; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredRhoUiUj(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredRhoUiUj(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredRhoUiUj(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    //! filter RhoUiUj in K direction
    if (filterDirection[2] == 1)
    {
        for (int n = 1; n <= 6; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredRhoUiUj(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredRhoUiUj(i, j, K, n);                
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredRhoUiUj(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeLeonardStress(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    RDouble4D &leonardStress = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("leonardStress"));
    RDouble4D &testFilteredRhoUiUj = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoUiUj"));
    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));

    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (isFVMOrFDM == FV_METHOD)
    {
        ComputeTestFilteredQ(grid);
    }
    else
    {
        ComputeTestFilteredQFD(grid);
    }

    Param_LESSolverStruct *parameters = GetControlParameters();
    int subgridScaleModel = parameters->GetSubgridScaleModel();

    if (subgridScaleModel == 4)
    {
        testFilteredQ(I, J, K, 0) = 1.0;
    }

    if (isFVMOrFDM == FV_METHOD)
    {
        ComputeTestFilteredRhoUiUj(grid);
    }
    else
    {
        ComputeTestFilteredRhoUiUjFD(grid);
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble tFRho = testFilteredQ(i, j, k, 0);

                RDouble tFU = testFilteredQ(i, j, k, 1);
                RDouble tFV = testFilteredQ(i, j, k, 2);
                RDouble tFW = testFilteredQ(i, j, k, 3);

                RDouble tFRhoUU = testFilteredRhoUiUj(i, j, k, 1);
                RDouble tFRhoVV = testFilteredRhoUiUj(i, j, k, 2);
                RDouble tFRhoWW = testFilteredRhoUiUj(i, j, k, 3);
                RDouble tFRhoUV = testFilteredRhoUiUj(i, j, k, 4);
                RDouble tFRhoUW = testFilteredRhoUiUj(i, j, k, 5);
                RDouble tFRhoVW = testFilteredRhoUiUj(i, j, k, 6);

                leonardStress(i, j, k, 1) = tFRhoUU - tFRho * tFU * tFU;
                leonardStress(i, j, k, 2) = tFRhoVV - tFRho * tFV * tFV;
                leonardStress(i, j, k, 3) = tFRhoWW - tFRho * tFW * tFW;
                leonardStress(i, j, k, 4) = tFRhoUV - tFRho * tFU * tFV;
                leonardStress(i, j, k, 5) = tFRhoUW - tFRho * tFU * tFW;
                leonardStress(i, j, k, 6) = tFRhoVW - tFRho * tFV * tFW;
            }
        }
    }
}

void LESSolverStruct::ComputeModeledStress(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Range N(1, 6);

    RDouble4D &modeledStress = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("modeledStress"));
    RDouble4D &testFilteredAlphaIJ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredAlphaIJ"));
    RDouble4D &betaIJ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("betaIJ"));

    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (isFVMOrFDM == FV_METHOD)
    {
        ComputeTestFilteredAlphaIJ(grid);
    }
    else
    {
        ComputeTestFilteredAlphaIJFD(grid);
    }

    ComputeBetaIJ(grid);

    for (int m = 1; m <= 6; ++m)
    {
        for (int k = kst; k <= ked; ++k)
        {
            for (int j = jst; j <= jed; ++j)
            {
                for (int i = ist; i <= ied; ++i)
                {
                    RDouble beta = betaIJ(i, j, k, m);
                    RDouble alpha = testFilteredAlphaIJ(i, j, k, m); 
                    modeledStress(i, j, k, m) = beta - alpha; 
                }
            }
        }
    }
}

void LESSolverStruct::ComputeTestFilteredAlphaIJ(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();
    int subgridScaleModel = parameters->GetSubgridScaleModel();

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble4D &strainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));

    RDouble4D &testFilteredAlphaIJ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("testFilteredAlphaIJ"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);

    if (subgridScaleModel == 2)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble rho = q(i, j, k, IR);
                    RDouble s = strainRateTensor(i, j, k, 0);
                    RDouble s11 = strainRateTensor(i, j, k, 1);
                    RDouble s22 = strainRateTensor(i, j, k, 2);
                    RDouble s33 = strainRateTensor(i, j, k, 3);
                    RDouble s12 = strainRateTensor(i, j, k, 4);
                    RDouble s13 = strainRateTensor(i, j, k, 5);
                    RDouble s23 = strainRateTensor(i, j, k, 6);

                    testFilteredAlphaIJ(i, j, k, 1) = -2.0 * rho * s * (s11 - 1.0/3.0 * (s11 + s22 + s33));
                    testFilteredAlphaIJ(i, j, k, 2) = -2.0 * rho * s * (s22 - 1.0/3.0 * (s11 + s22 + s33));
                    testFilteredAlphaIJ(i, j, k, 3) = -2.0 * rho * s * (s33 - 1.0/3.0 * (s11 + s22 + s33));
                    testFilteredAlphaIJ(i, j, k, 4) = -2.0 * rho * s * s12;
                    testFilteredAlphaIJ(i, j, k, 5) = -2.0 * rho * s * s13;
                    testFilteredAlphaIJ(i, j, k, 6) = -2.0 * rho * s * s23;

                    testFilteredAlphaIJ(i, j, k, 0) = 2.0 * rho * s * s;
                }
            }
        }
    }
    else if (subgridScaleModel == 4)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble s = strainRateTensor(i, j, k, 0);
                    RDouble s11 = strainRateTensor(i, j, k, 1);
                    RDouble s22 = strainRateTensor(i, j, k, 2);
                    RDouble s33 = strainRateTensor(i, j, k, 3);
                    RDouble s12 = strainRateTensor(i, j, k, 4);
                    RDouble s13 = strainRateTensor(i, j, k, 5);
                    RDouble s23 = strainRateTensor(i, j, k, 6);

                    testFilteredAlphaIJ(i, j, k, 1) = -2.0 * s * s11;
                    testFilteredAlphaIJ(i, j, k, 2) = -2.0 * s * s22;
                    testFilteredAlphaIJ(i, j, k, 3) = -2.0 * s * s33;
                    testFilteredAlphaIJ(i, j, k, 4) = -2.0 * s * s12;
                    testFilteredAlphaIJ(i, j, k, 5) = -2.0 * s * s13;
                    testFilteredAlphaIJ(i, j, k, 6) = -2.0 * s * s23;

                    testFilteredAlphaIJ(i, j, k, 0) = 2.0 * s * s;
                }
            }
        }
    }
    else
    {
        TK_Exit::UnexpectedVarValue("subgridScaleModel", subgridScaleModel);
    }

    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    //! filter AlphaIJ in I direction
    if (filterDirection[0] == 1)
    {
        for (int n = 0; n <= 6; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredAlphaIJ(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredAlphaIJ(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredAlphaIJ(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    //! filter AlphaIJ in J direction
    if (filterDirection[1] == 1)
    {
        for (int n = 0; n <= 6; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredAlphaIJ(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredAlphaIJ(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredAlphaIJ(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    //! filter AlphaIJ in K direction
    if (filterDirection[2] == 1)
    {
        for (int n = 0; n <= 6; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredAlphaIJ(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredAlphaIJ(i, j, K, n);                
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredAlphaIJ(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeTestFilteredAlphaIJFD(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();
    int subgridScaleModel = parameters->GetSubgridScaleModel();

    RDouble4D &q = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_FD"));
    RDouble4D &strainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));

    RDouble4D &testFilteredAlphaIJ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("testFilteredAlphaIJ"));

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);

    if (subgridScaleModel == 2)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble rho = q(i, j, k, IR);
                    RDouble s = strainRateTensor(i, j, k, 0);
                    RDouble s11 = strainRateTensor(i, j, k, 1);
                    RDouble s22 = strainRateTensor(i, j, k, 2);
                    RDouble s33 = strainRateTensor(i, j, k, 3);
                    RDouble s12 = strainRateTensor(i, j, k, 4);
                    RDouble s13 = strainRateTensor(i, j, k, 5);
                    RDouble s23 = strainRateTensor(i, j, k, 6);

                    testFilteredAlphaIJ(i, j, k, 1) = -2.0 * rho * s * (s11 - 1.0/3.0 * (s11 + s22 + s33));
                    testFilteredAlphaIJ(i, j, k, 2) = -2.0 * rho * s * (s22 - 1.0/3.0 * (s11 + s22 + s33));
                    testFilteredAlphaIJ(i, j, k, 3) = -2.0 * rho * s * (s33 - 1.0/3.0 * (s11 + s22 + s33));
                    testFilteredAlphaIJ(i, j, k, 4) = -2.0 * rho * s * s12;
                    testFilteredAlphaIJ(i, j, k, 5) = -2.0 * rho * s * s13;
                    testFilteredAlphaIJ(i, j, k, 6) = -2.0 * rho * s * s23;

                    testFilteredAlphaIJ(i, j, k, 0) = 2.0 * rho * s * s;
                }
            }
        }
    }
    else
    {
        TK_Exit::UnexpectedVarValue("subgridScaleModel", subgridScaleModel);
    }

    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);
    //! filter AlphaIJ in I direction
    if (filterDirection[0] == 1)
    {
        for (int n = 0; n <= 6; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredAlphaIJ(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredAlphaIJ(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredAlphaIJ(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    //! filter AlphaIJ in J direction
    if (filterDirection[1] == 1)
    {
        for (int n = 0; n <= 6; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredAlphaIJ(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredAlphaIJ(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredAlphaIJ(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    //! filter AlphaIJ in K direction
    if (filterDirection[2] == 1)
    {
        for (int n = 0; n <= 6; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredAlphaIJ(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredAlphaIJ(i, j, K, n);                
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredAlphaIJ(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeBetaIJ(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Param_LESSolverStruct *parameters = GetControlParameters();
    RDouble testFilterScale = parameters->GetTestFilterScale();
    int subgridScaleModel = parameters->GetSubgridScaleModel();

    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));
    RDouble4D &testFilteredStrainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredStrainRateTensor"));

    RDouble4D &betaIJ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("betaIJ"));
    
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (isFVMOrFDM == FV_METHOD)
    {
        ComputeTestFilteredStrainRateTensor(grid);
    }
    else
    {
        ComputeTestFilteredStrainRateTensor(grid);
    }

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    if (subgridScaleModel == 2)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble rho = testFilteredQ(i, j, k, 0);
                    RDouble s = testFilteredStrainRateTensor(i, j, k, 0);
                    RDouble s11 = testFilteredStrainRateTensor(i, j, k, 1);
                    RDouble s22 = testFilteredStrainRateTensor(i, j, k, 2);
                    RDouble s33 = testFilteredStrainRateTensor(i, j, k, 3);
                    RDouble s12 = testFilteredStrainRateTensor(i, j, k, 4);
                    RDouble s13 = testFilteredStrainRateTensor(i, j, k, 5);
                    RDouble s23 = testFilteredStrainRateTensor(i, j, k, 6);

                    betaIJ(i, j, k, 1) = -2.0 * testFilterScale * testFilterScale * rho * s * (s11 - 1.0/3.0 * (s11 + s22 + s33));
                    betaIJ(i, j, k, 2) = -2.0 * testFilterScale * testFilterScale * rho * s * (s22 - 1.0/3.0 * (s11 + s22 + s33));
                    betaIJ(i, j, k, 3) = -2.0 * testFilterScale * testFilterScale * rho * s * (s33 - 1.0/3.0 * (s11 + s22 + s33));
                    betaIJ(i, j, k, 4) = -2.0 * testFilterScale * testFilterScale * rho * s * s12;
                    betaIJ(i, j, k, 5) = -2.0 * testFilterScale * testFilterScale * rho * s * s13;
                    betaIJ(i, j, k, 6) = -2.0 * testFilterScale * testFilterScale * rho * s * s23;
                }
            }
        }
    }
    else if (subgridScaleModel == 4)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble s = testFilteredStrainRateTensor(i, j, k, 0);
                    RDouble s11 = testFilteredStrainRateTensor(i, j, k, 1);
                    RDouble s22 = testFilteredStrainRateTensor(i, j, k, 2);
                    RDouble s33 = testFilteredStrainRateTensor(i, j, k, 3);
                    RDouble s12 = testFilteredStrainRateTensor(i, j, k, 4);
                    RDouble s13 = testFilteredStrainRateTensor(i, j, k, 5);
                    RDouble s23 = testFilteredStrainRateTensor(i, j, k, 6);

                    betaIJ(i, j, k, 1) = -2.0 * testFilterScale * testFilterScale * s * s11;
                    betaIJ(i, j, k, 2) = -2.0 * testFilterScale * testFilterScale * s * s22;
                    betaIJ(i, j, k, 3) = -2.0 * testFilterScale * testFilterScale * s * s33;
                    betaIJ(i, j, k, 4) = -2.0 * testFilterScale * testFilterScale * s * s12;
                    betaIJ(i, j, k, 5) = -2.0 * testFilterScale * testFilterScale * s * s13;
                    betaIJ(i, j, k, 6) = -2.0 * testFilterScale * testFilterScale * s * s23;
                }
            }
        }
    }
    else
    {
        TK_Exit::UnexpectedVarValue("subgridScaleModel", subgridScaleModel);
    }
}

void LESSolverStruct::ComputeAnisotropicConstant(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int subgridScaleModel = parameters->GetSubgridScaleModel();

    RDouble3D &anisotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));
    RDouble4D &leonardStress = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("leonardStress"));
    RDouble4D &modeledStress = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("modeledStress"));

    ComputeLeonardStress(grid);

    ComputeModeledStress(grid);

    RDouble3D *LijMij = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *MijMij = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *LmmMnn = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *LijMijAvg = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *MijMijAvg = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *LmmMnnAvg = new RDouble3D(I, J, K, fortranArray);

    if (subgridScaleModel == 2)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble l11 = leonardStress(i, j, k, 1);
                    RDouble l22 = leonardStress(i, j, k, 2);
                    RDouble l33 = leonardStress(i, j, k, 3);
                    RDouble l12 = leonardStress(i, j, k, 4);
                    RDouble l13 = leonardStress(i, j, k, 5);
                    RDouble l23 = leonardStress(i, j, k, 6);

                    RDouble m11 = modeledStress(i, j, k, 1);
                    RDouble m22 = modeledStress(i, j, k, 2);
                    RDouble m33 = modeledStress(i, j, k, 3);
                    RDouble m12 = modeledStress(i, j, k, 4);
                    RDouble m13 = modeledStress(i, j, k, 5);
                    RDouble m23 = modeledStress(i, j, k, 6);

                    RDouble lm = l11 * m11 + l22 * m22 + l33 * m33 + 2.0 * (l12 * m12 + l13 * m13 + l23 * m23);
                    RDouble mm = m11 * m11 + m22 * m22 + m33 * m33 + 2.0 * (m12 * m12 + m13 * m13 + m23 * m23);
                    RDouble lmmn = (l11 + l22 + l33) * (m11 + m22 + m33);

                    (*LijMij)(i, j, k) = lm;
                    (*MijMij)(i, j, k) = mm;
                    (*LmmMnn)(i, j, k) = lmmn;
                }
            }
        }
    }
    else if (subgridScaleModel == 4)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble l11 = leonardStress(i, j, k, 1);
                    RDouble l22 = leonardStress(i, j, k, 2);
                    RDouble l33 = leonardStress(i, j, k, 3);
                    RDouble l12 = leonardStress(i, j, k, 4);
                    RDouble l13 = leonardStress(i, j, k, 5);
                    RDouble l23 = leonardStress(i, j, k, 6);

                    RDouble m11 = modeledStress(i, j, k, 1);
                    RDouble m22 = modeledStress(i, j, k, 2);
                    RDouble m33 = modeledStress(i, j, k, 3);
                    RDouble m12 = modeledStress(i, j, k, 4);
                    RDouble m13 = modeledStress(i, j, k, 5);
                    RDouble m23 = modeledStress(i, j, k, 6);

                    RDouble lm = l11 * m11 + l22 * m22 + l33 * m33 + 2.0 * (l12 * m12 + l13 * m13 + l23 * m23);
                    RDouble mm = m11 * m11 + m22 * m22 + m33 * m33 + 2.0 * (m12 * m12 + m13 * m13 + m23 * m23);

                    (*LijMij)(i, j, k) = lm;
                    (*MijMij)(i, j, k) = mm;
                    (*LmmMnn)(i, j, k) = 0.0;
                }
            }
        }
    }
    else
    {
        TK_Exit::UnexpectedVarValue("subgridScaleModel", subgridScaleModel);
    }

    (*LijMijAvg)(I, J, K) = (*LijMij)(I, J, K);
    (*MijMijAvg)(I, J, K) = (*MijMij)(I, J, K);
    (*LmmMnnAvg)(I, J, K) = (*LmmMnn)(I, J, K);

    int *averageDirection = parameters->GetAverageDirection();
    
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    if (averageDirection[0] == 0 && averageDirection[1] == 0 && averageDirection[2] == 0)
    {
        int averageWidth = parameters->GetAverageWidth();

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble LMAvg = 0.0;
                    RDouble MMAvg = 0.0;
                    RDouble LmMnAvg = 0.0;

                    int nAvg = 0;

                    for (int ka = -averageWidth; ka <= averageWidth; ++ ka)
                    {
                        for (int ja = -averageWidth; ja <= averageWidth; ++ ja)
                        {
                            for (int ia = -averageWidth; ia <= averageWidth; ++ ia)
                            {
                                int kk = k + ka;
                                int jj = j + ja;
                                int ii = i + ia;

                                bool ifAverage = (ii >= ist && ii <= ied) && (jj >= jst && jj <= jed) && (kk >= kst && kk <= ked);

                                if (ifAverage)
                                {
                                    LMAvg = LMAvg + (*LijMij)(ii, jj, kk);
                                    MMAvg = MMAvg + (*MijMij)(ii, jj, kk);
                                    LmMnAvg = LmMnAvg + (*LmmMnn)(ii, jj, kk);
                                    nAvg = nAvg + 1;
                                }
                            }
                        }
                    }

                    (*LijMijAvg)(i, j, k) = LMAvg / static_cast<RDouble>(nAvg);
                    (*MijMijAvg)(i, j, k) = MMAvg / static_cast<RDouble>(nAvg);
                    (*LmmMnnAvg)(i, j, k) = LmMnAvg / static_cast<RDouble>(nAvg);
                }
            }
        }
    }
    else
    {
        //! average in I direction
        if (averageDirection[0] != 0)
        {
            int nAvg = ied - ist + 1;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    RDouble LMAvg = 0.0;
                    RDouble MMAvg = 0.0;
                    RDouble LmMnAvg = 0.0;

                    for (int i = ist; i <= ied; ++ i)
                    {
                        LMAvg = LMAvg + (*LijMij)(i, j, k);
                        MMAvg = MMAvg + (*MijMij)(i, j, k);
                        LmMnAvg = LmMnAvg + (*LmmMnn)(i, j, k);
                    }

                    for (int i = ist; i <= ied; ++ i)
                    {
                        (*LijMijAvg)(i, j, k) = LMAvg / static_cast <RDouble>(nAvg);
                        (*MijMijAvg)(i, j, k) = MMAvg / static_cast <RDouble>(nAvg);
                        (*LmmMnnAvg)(i, j, k) = LmMnAvg / static_cast <RDouble>(nAvg);
                    }

                    for (int i = ist; i <= ied; ++ i)
                    {
                        (*LijMij)(i, j, k) = (*LijMijAvg)(i, j, k);
                        (*MijMij)(i, j, k) = (*MijMijAvg)(i, j, k);
                        (*LmmMnn)(i, j, k) = (*LmmMnnAvg)(i, j, k);
                    }
                }
            }
        }

        RDouble lm = (*LijMijAvg)(1, 1, 1);
        RDouble mm = (*MijMijAvg)(1, 1, 1);

        //! average in J direction
        if (averageDirection[1] != 0)
        {
            int nAvg = jed - jst + 1;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble LMAvg = 0.0;
                    RDouble MMAvg = 0.0;
                    RDouble LmMnAvg = 0.0;

                    for (int j = jst; j <= jed; ++ j)
                    {
                        LMAvg = LMAvg + (*LijMij)(i, j, k);
                        MMAvg = MMAvg + (*MijMij)(i, j, k);
                        LmMnAvg = LmMnAvg + (*LmmMnn)(i, j, k);
                    }

                    for (int j = jst; j <= jed; ++ j)
                    {
                        (*LijMijAvg)(i, j, k) = LMAvg / static_cast <RDouble>(nAvg);
                        (*MijMijAvg)(i, j, k) = MMAvg / static_cast <RDouble>(nAvg);
                        (*LmmMnnAvg)(i, j, k) = LmMnAvg / static_cast <RDouble>(nAvg);
                    }

                    for (int j = jst; j <= jed; ++ j)
                    {
                        (*LijMij)(i, j, k) = (*LijMijAvg)(i, j, k);
                        (*MijMij)(i, j, k) = (*MijMijAvg)(i, j, k);
                        (*LmmMnn)(i, j, k) = (*LmmMnnAvg)(i, j, k);
                    }
                }
            }
        }

        lm = (*LijMijAvg)(1, 1, 1);
        mm = (*MijMijAvg)(1, 1, 1);

        //! average in K direction
        if (averageDirection[2] != 0)
        {
            int nAvg = ked - kst + 1;

            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble LMAvg = 0.0;
                    RDouble MMAvg = 0.0;
                    RDouble LmMnAvg = 0.0;

                    for (int k = kst; k <= ked; ++ k)
                    {
                        LMAvg = LMAvg + (*LijMij)(i, j, k);
                        MMAvg = MMAvg + (*MijMij)(i, j, k);
                        LmMnAvg = LmMnAvg + (*LmmMnn)(i, j, k);
                    }

                    for (int k = kst; k <= ked; ++ k)
                    {
                        (*LijMijAvg)(i, j, k) = LMAvg / static_cast <RDouble>(nAvg);
                        (*MijMijAvg)(i, j, k) = MMAvg / static_cast <RDouble>(nAvg);
                        (*LmmMnnAvg)(i, j, k) = LmMnAvg / static_cast <RDouble>(nAvg);
                    }

                    for (int k = kst; k <= ked; ++ k)
                    {
                        (*LijMij)(i, j, k) = (*LijMijAvg)(i, j, k);
                        (*MijMij)(i, j, k) = (*MijMijAvg)(i, j, k);
                        (*LmmMnn)(i, j, k) = (*LmmMnnAvg)(i, j, k);
                    }
                }
            }
        }
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble lm = (*LijMijAvg)(i, j, k);
                RDouble mm = (*MijMijAvg)(i, j, k);
                RDouble lmMn = (*LmmMnnAvg)(i, j, k);

                anisotropicConstant(i, j, k) = (lm - 1.0/3.0 * lmMn) / mm; 
            }
        }
    }

    //GhostCell3D(anisotropicConstant, ni, nj, nk);

    delete LijMij;
    delete MijMij;
    delete LmmMnn;
    delete LijMijAvg;
    delete MijMijAvg;
    delete LmmMnnAvg;
}

void LESSolverStruct::ComputeIsotropicConstant(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    RDouble testFilterScale = parameters->GetTestFilterScale();
    int *averageDirection = parameters->GetAverageDirection();

    RDouble3D &isotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("isotropicConstant"));
    RDouble4D &leonardStress = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("leonardStress"));
    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));
    RDouble4D &testFilteredAlphaIJ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredAlphaIJ"));
    RDouble4D &testFilteredStrainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredStrainRateTensor"));

    RDouble3D *numerator = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *denominator = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *numeratorAvg = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *denominatorAvg = new RDouble3D(I, J, K, fortranArray);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble l11 = leonardStress(i, j, k, 1);
                RDouble l22 = leonardStress(i, j, k, 2);
                RDouble l33 = leonardStress(i, j, k, 3);
                (*numerator)(i, j, k) = l11 + l22 + l33;

                RDouble tFRho = testFilteredQ(i, j, k, 0);
                RDouble testFilteredStranRateMagnitude = testFilteredStrainRateTensor(i, j, k, 0);
                RDouble testFilteredAlpha = testFilteredAlphaIJ(i, j, k, 0);
                (*denominator)(i, j, k) = 2.0 * testFilterScale * testFilterScale * tFRho * testFilteredStranRateMagnitude * testFilteredStranRateMagnitude - testFilteredAlpha;
            }
        }
    }

    (*numeratorAvg)(I, J, K) = (*numerator)(I, J, K);
    (*denominatorAvg)(I, J, K) = (*denominator)(I, J, K);

    if (averageDirection[0] == 0 && averageDirection[1] == 0 && averageDirection[2] == 0)
    {
        int averageWidth = parameters->GetAverageWidth();
        
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble num = 0.0;
                    RDouble den = 0.0;

                    int nAvg = 0;

                    for (int ka = -averageWidth; ka <= averageWidth; ++ ka)
                    {
                        for (int ja = -averageWidth; ja <= averageWidth; ++ ja)
                        {
                            for (int ia = -averageWidth; ia <= averageWidth; ++ ia)
                            {
                                int kk = k + ka;
                                int jj = j + ja;
                                int ii = i + ia;

                                bool ifAverage = (ii >= ist && ii <= ied) && (jj >= jst && jj <= jed) && (kk >= kst && kk <= ked);

                                if (ifAverage)
                                {
                                    num = num + (*numerator)(ii, jj, kk);
                                    den = den + (*denominator)(ii, jj, kk);
                                    nAvg = nAvg + 1;
                                }
                            }
                        }
                    }

                    (*numeratorAvg)(i, j, k) = num / nAvg;
                    (*denominatorAvg)(i, j, k) = den / nAvg;
                }
            }
        }
    }
    else
    {
        //! average in I direction
        if (averageDirection[0] != 0)
        {
            int nAvg = ied - ist + 1;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    RDouble num = 0.0;
                    RDouble den = 0.0;

                    for (int i = ist; i <= ied; ++ i)
                    {
                        num = num + (*numerator)(i, j, k);
                        den = den + (*denominator)(i, j, k);
                    }

                    for (int i = ist; i <= ied; ++ i)
                    {
                        (*numeratorAvg)(i, j, k) = num / static_cast<RDouble>(nAvg);
                        (*denominatorAvg)(i, j, k) = den / static_cast<RDouble>(nAvg);
                    }

                    for (int i = ist; i <= ied; ++ i)
                    {
                        (*numerator)(i, j, k) = (*numeratorAvg)(i, j, k);
                        (*denominator)(i, j, k) = (*denominatorAvg)(i, j, k);
                    }
                }
            }
        }

        //! average in J direction
        if (averageDirection[1] != 0)
        {
            int nAvg = jed - jst + 1;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble num = 0.0;
                    RDouble den = 0.0;

                    for (int j = jst; j <= jed; ++ j)
                    {
                        num = num + (*numerator)(i, j, k);
                        den = den + (*denominator)(i, j, k);
                    }

                    for (int j = jst; j <= jed; ++ j)
                    {
                        (*numeratorAvg)(i, j, k) = num / static_cast<RDouble>(nAvg);
                        (*denominatorAvg)(i, j, k) = den / static_cast<RDouble>(nAvg);
                    }

                    for (int j = jst; j <= jed; ++ j)
                    {
                        (*numerator)(i, j, k) = (*numeratorAvg)(i, j, k);
                        (*denominator)(i, j, k) = (*denominatorAvg)(i, j, k);
                    }
                }
            }
        }

        //! average in K direction
        if (averageDirection[2] != 0)
        {
            int nAvg = ked - kst + 1;

            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble num = 0.0;
                    RDouble den = 0.0;

                    for (int k = kst; k <= ked; ++ k)
                    {
                        num = num + (*numerator)(i, j, k);
                        den = den + (*denominator)(i, j, k);
                    }

                    for (int k = kst; k <= ked; ++ k)
                    {
                        (*numeratorAvg)(i, j, k) = num / static_cast<RDouble>(nAvg);
                        (*denominatorAvg)(i, j, k) = den / static_cast<RDouble>(nAvg);
                    }

                    for (int k = kst; k <= ked; ++ k)
                    {
                        (*numerator)(i, j, k) = (*numeratorAvg)(i, j, k);
                        (*denominator)(i, j, k) = (*denominatorAvg)(i, j, k);
                    }
                }
            }
        }
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble num = (*numeratorAvg)(i, j, k);
                RDouble den = (*denominatorAvg)(i, j, k);
                isotropicConstant(i, j, k) = num / den;
            }
        }
    }

    delete numerator;
    delete denominator;
    delete numeratorAvg;
    delete denominatorAvg;
}

void LESSolverStruct::ComputePrandtlNumber(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *averageDirection = parameters->GetAverageDirection();

    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));

    RDouble3D &anisotropicConstant = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("anisotropicConstant"));

    RDouble4D &leonardTemperatureFlux = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("leonardTemperatureFlux"));
    RDouble4D &modeledTemperatureFlux = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("modeledTemperatureFlux"));

    RDouble3D *numerator = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *denominator = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *numeratorAvg = new RDouble3D(I, J, K, fortranArray);
    RDouble3D *denominatorAvg = new RDouble3D(I, J, K, fortranArray);

    ComputeLeonardTemperatureTerm(grid);

    ComputeModeledTemperatureTerm(grid);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble mTF0 = modeledTemperatureFlux(i, j, k, 0);
                RDouble mTF1 = modeledTemperatureFlux(i, j, k, 1);
                RDouble mTF2 = modeledTemperatureFlux(i, j, k, 2);

                RDouble lTF0 = leonardTemperatureFlux(i, j, k, 0);
                RDouble lTF1 = leonardTemperatureFlux(i, j, k, 1);
                RDouble lTF2 = leonardTemperatureFlux(i, j, k, 2);

                (*numerator)(i, j, k) = mTF0 * mTF0 + mTF1 * mTF1 + mTF2 * mTF2;
                (*denominator)(i, j, k) = lTF0 * mTF0 + lTF1 * mTF1 + lTF2 * mTF2;
            }
        }
    }

    (*numeratorAvg)(I, J, K) = (*numerator)(I, J, K);
    (*denominatorAvg)(I, J, K) = (*denominator)(I, J, K);

    if (averageDirection[0] == 0 && averageDirection[1] == 0 && averageDirection[2] == 0)
    {
        int averageWidth = parameters->GetAverageWidth();

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble num = 0.0;
                    RDouble den = 0.0;

                    int nAvg = 0;

                    for (int ka = -averageWidth; ka <= averageWidth; ++ ka)
                    {
                        for (int ja = -averageWidth; ja <= averageWidth; ++ ja)
                        {
                            for (int ia = -averageWidth; ia <= averageWidth; ++ ia)
                            {
                                int kk = k + ka;
                                int jj = j + ja;
                                int ii = i + ia;

                                bool ifAverage = (ii >= ist && ii <= ied) && (jj >= jst && jj <= jed) && (kk >= kst && kk <= ked);

                                if (ifAverage)
                                {
                                    num = num + (*numerator)(ii, jj, kk);
                                    den = den + (*denominator)(ii, jj, kk);
                                    nAvg = nAvg + 1;
                                }
                            }
                        }
                    }

                    (*numeratorAvg)(i, j, k) = num / nAvg;
                    (*denominatorAvg)(i, j, k) = den / nAvg;
                }
            }
        }
    }
    else
    {
        //! average in I direction
        if (averageDirection[0] != 0)
        {
            int nAvg = ied - ist + 1;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    RDouble num = 0.0;
                    RDouble den = 0.0;

                    for (int i = ist; i <= ied; ++ i)
                    {
                        num = num + (*numerator)(i, j, k);
                        den = den + (*denominator)(i, j, k);
                    }

                    for (int i = ist; i <= ied; ++ i)
                    {
                        (*numeratorAvg)(i, j, k) = num / nAvg;
                        (*denominatorAvg)(i, j, k) = den / nAvg;
                    }

                    for (int i = ist; i <= ied; ++ i)
                    {
                        (*numerator)(i, j, k) = (*numeratorAvg)(i, j, k);
                        (*denominator)(i, j, k) = (*denominatorAvg)(i, j, k);
                    }
                }
            }
        }

        //! average in J direction
        if (averageDirection[1] != 0)
        {
            int nAvg = jed - jst + 1;

            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble num = 0.0;
                    RDouble den = 0.0;

                    for (int j = jst; j <= jed; ++ j)
                    {
                        num = num + (*numerator)(i, j, k);
                        den = den + (*denominator)(i, j, k);
                    }

                    for (int j = jst; j <= jed; ++ j)
                    {
                        (*numeratorAvg)(i, j, k) = num / nAvg;
                        (*denominatorAvg)(i, j, k) = den / nAvg;
                    }

                    for (int j = jst; j <= jed; ++ j)
                    {
                        (*numerator)(i, j, k) = (*numeratorAvg)(i, j, k);
                        (*denominator)(i, j, k) = (*denominatorAvg)(i, j, k);
                    }
                }
            }
        }

        //! average in K direction
        if (averageDirection[2] != 0)
        {
            int nAvg = ked - kst + 1;

            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble num = 0.0;
                    RDouble den = 0.0;

                    for (int k = kst; k <= ked; ++ k)
                    {
                        num = num + (*numerator)(i, j, k);
                        den = den + (*denominator)(i, j, k);
                    }

                    for (int k = kst; k <= ked; ++ k)
                    {
                        (*numeratorAvg)(i, j, k) = num / nAvg;
                        (*denominatorAvg)(i, j, k) = den / nAvg;
                    }

                    for (int k = kst; k <= ked; ++ k)
                    {
                        (*numerator)(i, j, k) = (*numeratorAvg)(i, j, k);
                        (*denominator)(i, j, k) = (*denominatorAvg)(i, j, k);
                    }
                }
            }
        }
    }

    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble cr = anisotropicConstant(i, j, k);

                RDouble KjKj = (*numeratorAvg)(i, j, k);
                RDouble NjNj = (*denominatorAvg)(i, j, k);

                RDouble prt = cr * KjKj / NjNj;
                turbulentPrandtlNumber(i, j, k) = prt;
            }
        }
    }

    delete numerator;
    delete denominator;
    delete numeratorAvg;
    delete denominatorAvg;
}

void LESSolverStruct::ComputeLeonardTemperatureTerm(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    RDouble4D &leonardTemperatureFlux = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("leonardTemperatureFlux"));

    RDouble4D &testFilteredRhoUT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoUT"));
    RDouble3D &testFilteredTemperature = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("testFilteredTemperature"));
    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));

    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (isFVMOrFDM == FV_METHOD)
    {
        ComputeTestFilteredRhoUT(grid);

        ComputeTestFilteredTemperature(grid);
    }
    else
    {
        ComputeTestFilteredRhoUTFD(grid);

        ComputeTestFilteredTemperatureFD(grid);
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble tFRhoUT = testFilteredRhoUT(i, j, k, 1);
                RDouble tFRhoVT = testFilteredRhoUT(i, j, k, 2);
                RDouble tFRhoWT = testFilteredRhoUT(i, j, k, 3);

                RDouble tFRho = testFilteredQ(i, j, k, 0);

                RDouble tFU = testFilteredQ(i, j, k, 1);
                RDouble tFV = testFilteredQ(i, j, k, 2);
                RDouble tFW = testFilteredQ(i, j, k, 3);

                RDouble tFT = testFilteredTemperature(i, j, k);

                leonardTemperatureFlux(i, j, k, 0) = tFRhoUT - tFRho * tFU * tFT;
                leonardTemperatureFlux(i, j, k, 1) = tFRhoVT - tFRho * tFV * tFT;
                leonardTemperatureFlux(i, j, k, 2) = tFRhoWT - tFRho * tFW * tFT;
            }
        }
    }
}

void LESSolverStruct::ComputeModeledTemperatureTerm(Grid *gridIn)
{
    using namespace IDX;
 
    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    RDouble testFilterScale = parameters->GetTestFilterScale();

    RDouble4D &modeledTemperatureFlux = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("modeledTemperatureFlux"));

    RDouble4D &testFilteredQ = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredQ"));
    RDouble4D &testFilteredRhoStrainRateMagnitudeGradT = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoStrainRateMagnitudeGradT"));
    RDouble4D &testFilteredGradT = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredGradT"));
    RDouble4D &testFilteredStrainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredStrainRateTensor"));

    ComputeTestFilteredGradT(grid);

    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (isFVMOrFDM == FV_METHOD)
    {
        ComputeTestFilteredRhoStrainRateMagnitudeGradT(grid);
    }
    else
    {
        ComputeTestFilteredRhoStrainRateMagnitudeGradTFD(grid);
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble tFRSGradTx = testFilteredRhoStrainRateMagnitudeGradT(i, j, k, 0);
                RDouble tFRSGradTy = testFilteredRhoStrainRateMagnitudeGradT(i, j, k, 1);
                RDouble tFRSGradTz = testFilteredRhoStrainRateMagnitudeGradT(i, j, k, 2);

                RDouble tFGradTx = testFilteredGradT(i, j, k, 0);
                RDouble tFGradTy = testFilteredGradT(i, j, k, 1);
                RDouble tFGradTz = testFilteredGradT(i, j, k, 2);

                RDouble tFRho = testFilteredQ(i, j, k, 0);
                RDouble tFS = testFilteredStrainRateTensor(i, j, k, 0);

                modeledTemperatureFlux(i, j, k, 0) = tFRSGradTx - testFilterScale * testFilterScale * tFRho * tFS * tFGradTx;
                modeledTemperatureFlux(i, j, k, 1) = tFRSGradTy - testFilterScale * testFilterScale * tFRho * tFS * tFGradTy;
                modeledTemperatureFlux(i, j, k, 2) = tFRSGradTz - testFilterScale * testFilterScale * tFRho * tFS * tFGradTz;
            }
        }
    }
}

void LESSolverStruct::ComputeTestFilteredRhoUT(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &q = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble4D &temperatures = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("t"));

    RDouble4D &testFilteredRhoUT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoUT"));

    testFilteredRhoUT(I, J, K, 1) = q(I, J, K, IR) * q(I, J, K, IU) * temperatures(I, J, K, ITT);
    testFilteredRhoUT(I, J, K, 2) = q(I, J, K, IR) * q(I, J, K, IV) * temperatures(I, J, K, ITT);
    testFilteredRhoUT(I, J, K, 3) = q(I, J, K, IR) * q(I, J, K, IW) * temperatures(I, J, K, ITT);

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    //! filter rhoUT in I direction
    if (filterDirection[0] == 1)
    {
        for (int n = 1; n <= 3; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredRhoUT(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredRhoUT(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredRhoUT(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    //! filter rhoUT in J direction
    if (filterDirection[1] == 1)
    {
        for (int n = 1; n <= 3; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredRhoUT(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredRhoUT(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredRhoUT(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    //! filter rhoUT in K direction
    if (filterDirection[2] == 1)
    {
        for (int n = 1; n <= 3; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredRhoUT(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredRhoUT(i, j, K, n);                
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredRhoUT(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeTestFilteredRhoUTFD(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &q = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_FD"));
    RDouble4D &temperatures = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("t_FD"));

    RDouble4D &testFilteredRhoUT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoUT"));

    testFilteredRhoUT(I, J, K, 1) = q(I, J, K, IR) * q(I, J, K, IU) * temperatures(I, J, K, ITT);
    testFilteredRhoUT(I, J, K, 2) = q(I, J, K, IR) * q(I, J, K, IV) * temperatures(I, J, K, ITT);
    testFilteredRhoUT(I, J, K, 3) = q(I, J, K, IR) * q(I, J, K, IW) * temperatures(I, J, K, ITT);

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    //! filter rhoUT in I direction
    if (filterDirection[0] == 1)
    {
        for (int n = 1; n <= 3; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredRhoUT(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredRhoUT(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredRhoUT(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    //! filter rhoUT in J direction
    if (filterDirection[1] == 1)
    {
        for (int n = 1; n <= 3; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredRhoUT(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredRhoUT(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredRhoUT(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    //! filter rhoUT in K direction
    if (filterDirection[2] == 1)
    {
        for (int n = 1; n <= 3; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredRhoUT(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredRhoUT(i, j, K, n);                
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredRhoUT(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeTestFilteredGradT(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    //grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &gradT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradT"));
    RDouble4D &testFilteredGradT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredGradT"));

    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (isFVMOrFDM == FV_METHOD)
    {
        ComputeGradT(grid);
    }
    else
    {
        ComputeGradTFD(grid);
    }

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    testFilteredGradT(I, J, K, Range(0, 2)) = gradT(I, J, K, Range(0, 2));

    //! filter gradT in I direction
    if (filterDirection[0] == 1)
    {
        for (int n = 0; n <= 2; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredGradT(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredGradT(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredGradT(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    //! filter gradT in J direction
    if (filterDirection[1] == 1)
    {
        for (int n = 1; n <= 2; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredGradT(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredGradT(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredGradT(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    //! filter gradT in K direction
    if (filterDirection[2] == 1)
    {
        for (int n = 1; n <= 2; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredGradT(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredGradT(i, j, K, n);
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredGradT(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeTestFilteredRhoStrainRateMagnitudeGradT(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &q = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble4D &testFilteredRhoStrainRateMagnitudeGradT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoStrainRateMagnitudeGradT"));    
    RDouble4D &strainRateTensor = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));
    RDouble4D &gradT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradT"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble rho = q(i, j, k, IR);
                RDouble strainRateMagnitude = strainRateTensor(i, j, k, 0);
                RDouble dtdx = gradT(i, j, k, 0);
                RDouble dtdy = gradT(i, j, k, 1);
                RDouble dtdz = gradT(i, j, k, 2);

                testFilteredRhoStrainRateMagnitudeGradT(i, j, k, 0) = rho * strainRateMagnitude * dtdx;
                testFilteredRhoStrainRateMagnitudeGradT(i, j, k, 1) = rho * strainRateMagnitude * dtdy;
                testFilteredRhoStrainRateMagnitudeGradT(i, j, k, 2) = rho * strainRateMagnitude * dtdz;
            }
        }
    }

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);
    //! filter strainRateMagnitudeGradT in I direction
    if (filterDirection[0] == 1)
    {
        for (int n = 0; n <= 2; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredRhoStrainRateMagnitudeGradT(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredRhoStrainRateMagnitudeGradT(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredRhoStrainRateMagnitudeGradT(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    //! filter strainRateMagnitudeGradT in J direction
    if (filterDirection[1] == 1) 
    {
        for (int n = 0; n <= 2; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredRhoStrainRateMagnitudeGradT(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredRhoStrainRateMagnitudeGradT(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredRhoStrainRateMagnitudeGradT(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    //! filter strainRateMagnitudeGradT in K direction
    if (filterDirection[2] == 1)
    {
        for (int n = 0; n <= 2; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredRhoStrainRateMagnitudeGradT(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredRhoStrainRateMagnitudeGradT(i, j, K, n);
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredRhoStrainRateMagnitudeGradT(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::ComputeTestFilteredRhoStrainRateMagnitudeGradTFD(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Param_LESSolverStruct *parameters = GetControlParameters();
    int *filterDirection = parameters->GetFilterDirection();

    RDouble4D &q = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_FD"));

    RDouble4D &testFilteredRhoStrainRateMagnitudeGradT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("testFilteredRhoStrainRateMagnitudeGradT"));    
    RDouble4D &strainRateTensor = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));
    RDouble4D &gradT = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradT"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble rho = q(i, j, k, IR);
                RDouble strainRateMagnitude = strainRateTensor(i, j, k, 0);
                RDouble dtdx = gradT(i, j, k, 0);
                RDouble dtdy = gradT(i, j, k, 1);
                RDouble dtdz = gradT(i, j, k, 2);

                testFilteredRhoStrainRateMagnitudeGradT(i, j, k, 0) = rho * strainRateMagnitude * dtdx;
                testFilteredRhoStrainRateMagnitudeGradT(i, j, k, 1) = rho * strainRateMagnitude * dtdy;
                testFilteredRhoStrainRateMagnitudeGradT(i, j, k, 2) = rho * strainRateMagnitude * dtdz;
            }
        }
    }

    RDouble1D *testFilteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *testFilteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *testFilteredVariableZ = new RDouble1D(K, fortranArray);
    RDouble1D *filteredVariableX = new RDouble1D(I, fortranArray);
    RDouble1D *filteredVariableY = new RDouble1D(J, fortranArray);
    RDouble1D *filteredVariableZ = new RDouble1D(K, fortranArray);

    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);
    //! filter strainRateMagnitudeGradT in I direction
    if (filterDirection[0] == 1)
    {
        for (int n = 0; n <= 2; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    (*testFilteredVariableX)(I) = testFilteredRhoStrainRateMagnitudeGradT(I, j, k, n);
                    (*filteredVariableX)(I) = testFilteredRhoStrainRateMagnitudeGradT(I, j, k, n);
                    TestFilter(ist, ied, testFilteredVariableX, filteredVariableX);
                    testFilteredRhoStrainRateMagnitudeGradT(I, j, k, n) = (*testFilteredVariableX)(I);
                }
            }
        }
    }

    //! filter strainRateMagnitudeGradT in J direction
    if (filterDirection[1] == 1) 
    {
        for (int n = 0; n <= 2; ++ n)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableY)(J) = testFilteredRhoStrainRateMagnitudeGradT(i, J, k, n);
                    (*filteredVariableY)(J) = testFilteredRhoStrainRateMagnitudeGradT(i, J, k, n);
                    TestFilter(jst, jed, testFilteredVariableY, filteredVariableY);
                    testFilteredRhoStrainRateMagnitudeGradT(i, J, k, n) = (*testFilteredVariableY)(J);
                }
            }
        }
    }

    //! filter strainRateMagnitudeGradT in K direction
    if (filterDirection[2] == 1)
    {
        for (int n = 0; n <= 2; ++ n)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*testFilteredVariableZ)(K) = testFilteredRhoStrainRateMagnitudeGradT(i, j, K, n);
                    (*filteredVariableZ)(K) = testFilteredRhoStrainRateMagnitudeGradT(i, j, K, n);
                    TestFilter(kst, ked, testFilteredVariableZ, filteredVariableZ);
                    testFilteredRhoStrainRateMagnitudeGradT(i, j, K, n) = (*testFilteredVariableZ)(K);
                }
            }
        }
    }

    delete testFilteredVariableX;
    delete testFilteredVariableY;
    delete testFilteredVariableZ;
    delete filteredVariableX;
    delete filteredVariableY;
    delete filteredVariableZ;
}

void LESSolverStruct::WALE(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    RDouble3D &cellLengthScale = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("cellLengthScale"));
    RDouble3D &wallFunction = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));

    Param_LESSolverStruct *parameters = GetControlParameters();

    RDouble waleConstant = parameters->GetWaleConstant();

    RDouble iConstant = parameters->GetIsotropicConstant();

    RDouble refReNumber  = parameters->GetRefReNumber();

    RDouble4D &q = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble3D &viscousTurbulence = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble3D &subgridScaleEnergy = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("subgridScaleEnergy"));
    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));

    RDouble4D &strainRateTensor = * reinterpret_cast<RDouble4D *> (grid->GetDataPtr("strainRateTensor"));

    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {                
                RDouble deltaBar = cellLengthScale(i, j, k);

                RDouble dampingFunction = wallFunction(i, j, k);

                RDouble rho = ABS(q(i, j, k, IR)) + SMALL;

                RDouble s11 = strainRateTensor(i, j, k, 1);
                RDouble s22 = strainRateTensor(i, j, k, 2);
                RDouble s33 = strainRateTensor(i, j, k, 3);
                RDouble s12 = strainRateTensor(i, j, k, 4);
                RDouble s13 = strainRateTensor(i, j, k, 5);
                RDouble s23 = strainRateTensor(i, j, k, 6);

                RDouble omega12 = strainRateTensor(i, j, k, 7);
                RDouble omega13 = strainRateTensor(i, j, k, 8);
                RDouble omega23 = strainRateTensor(i, j, k, 9);

                RDouble sij2  = s11 * s11 + s22 * s22 + s33 * s33 + 2.0 * (s12 * s12 + s13 * s13 + s23 * s23);
                RDouble omegaij2 = 2.0 * (omega12 * omega12 + omega13 * omega13 + omega23 * omega23);
                RDouble sij2SubOmegaij2 = 1.0/3.0 * (sij2 - omegaij2);

                RDouble sd11 = s11 * s11 + s12 * s12 + s13 * s13 - omega12 * omega12 - omega13 * omega13 - sij2SubOmegaij2;
                RDouble sd22 = s12 * s12 + s22 * s22 + s23 * s23 - omega12 * omega12 - omega23 * omega23 - sij2SubOmegaij2;
                RDouble sd33 = s13 * s13 + s23 * s23 + s33 * s33 - omega13 * omega13 - omega23 * omega23 - sij2SubOmegaij2;
                RDouble sd12 = s11 * s12 + s12 * s22 + s13 * s23 - omega13 * omega23;
                RDouble sd13 = s11 * s13 + s12 * s23 + s13 * s33 + omega12 * omega23;
                RDouble sd23 = s12 * s13 + s22 * s23 + s23 * s33 - omega12 * omega13;

                RDouble sdij2 = sd11 * sd11 + sd22 * sd22 + sd33 * sd33 + 2.0 * (sd12 * sd12 + sd13 * sd13 + sd23 * sd23);

                RDouble ss = pow(sdij2, 1.5) / (pow(sij2, 2.5) + pow(sdij2, 1.25) + SMALL);

                RDouble mut = (waleConstant * deltaBar * dampingFunction) * (waleConstant * deltaBar * dampingFunction) * rho * ss * refReNumber;
                mut = MIN(eddyViscosityLimit, mut);
                RDouble kSGS = 2.0 * iConstant / rho * (mut /deltaBar) * (mut /deltaBar) / refReNumber;  //! Divide by one refReNumber here, and divide by another refReNumber in the Viscousflux.

                viscousTurbulence(i, j, k) = mut;
                subgridScaleEnergy(i, j, k) = kSGS;
            }
        }
    }

    RDouble prandtlTurbulence = parameters->GetPrandtlTurbulence();

    turbulentPrandtlNumber = prandtlTurbulence;
}

void LESSolverStruct::Sigma(Grid *gridIn)
{
    using namespace IDX;

    StructGrid *grid = StructGridCast(gridIn);

    RDouble3D &cellLengthScale = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("cellLengthScale"));
    RDouble3D &wallFunction = * reinterpret_cast<RDouble3D *> (grid->GetDataPtr("wallFunction"));

    Param_LESSolverStruct *parameters = GetControlParameters();

    RDouble sigmaConstant = parameters->GetSigmaConstant();

    RDouble refReNumber  = parameters->GetRefReNumber();

    RDouble4D &q = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble3D &viscousTurbulence = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble3D &subgridScaleEnergy = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("subgridScaleEnergy"));
    RDouble3D &turbulentPrandtlNumber = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("turbulentPrandtlNumber"));

    RDouble4D &gradUVWTCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &gradUVWTCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    int imax = 1;
    int jmax = 1;
    int kmax = 1;

    int imin = 1;
    int jmin = 1;
    int kmin = 1;

    RDouble viscousTurbulenceMaximum = -1.0e30;
    RDouble viscousTurbulenceMinimum = 1.0e30;

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble ** velocityGradTensor = AleModel :: CreateMatrix();

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                velocityGradTensor[0][0] = gradUVWTCellCenterX(i, j, k, 0);
                velocityGradTensor[0][1] = gradUVWTCellCenterY(i, j, k, 0);
                velocityGradTensor[0][2] = gradUVWTCellCenterZ(i, j, k, 0);

                velocityGradTensor[1][0] = gradUVWTCellCenterX(i, j, k, 1);
                velocityGradTensor[1][1] = gradUVWTCellCenterY(i, j, k, 1);
                velocityGradTensor[1][2] = gradUVWTCellCenterZ(i, j, k, 1);

                velocityGradTensor[2][0] = gradUVWTCellCenterX(i, j, k, 2);
                velocityGradTensor[2][1] = gradUVWTCellCenterY(i, j, k, 2);
                velocityGradTensor[2][2] = gradUVWTCellCenterZ(i, j, k, 2);

                RDouble differentialSigma = 0.0;
                ComputeDifferentialSigma(velocityGradTensor, differentialSigma);

                RDouble deltaBar = cellLengthScale(i, j, k);
                RDouble rho = ABS(q(i, j, k, IR)) + SMALL;
                RDouble dampingFunction = wallFunction(i, j, k);

                RDouble mut = rho * (sigmaConstant * deltaBar * dampingFunction) * (sigmaConstant * deltaBar * dampingFunction) * differentialSigma * refReNumber;
                mut = MIN(eddyViscosityLimit, mut);

                viscousTurbulence(i, j, k) = mut;
                subgridScaleEnergy(i, j, k) = 0.0;

                if (viscousTurbulenceMaximum < mut)
                {
                    viscousTurbulenceMaximum = mut;
                    imax = i;
                    jmax = j;
                    kmax = k;
                }

                if (viscousTurbulenceMinimum > mut)
                {
                    viscousTurbulenceMinimum = mut;
                    imin = i;
                    jmin = j;
                    kmin = k;
                }
            }
        }
    }

    RDouble prandtlTurbulence = parameters->GetPrandtlTurbulence();

    turbulentPrandtlNumber = prandtlTurbulence;

    AleModel::DestroyMatrix(velocityGradTensor);
}

void LESSolverStruct::ComputeDifferentialSigma(RDouble ** velocityGradTensor, RDouble &differentialSigma)
{
    RDouble ** velocityGradTensorTranspose = AleModel :: CreateMatrix();
    RDouble ** gij = AleModel :: CreateMatrix();
    RDouble ** gij2 = AleModel :: CreateMatrix();

    AleModel::SetMatrix(velocityGradTensorTranspose, velocityGradTensor);

    AleModel::TransposeMatrix(velocityGradTensorTranspose);

    AleModel::MatrixMultiply(velocityGradTensorTranspose, velocityGradTensor, gij);

    AleModel::MatrixMultiply(gij, gij, gij2);

    RDouble tracegij = gij[0][0] + gij[1][1] + gij[2][2];
    RDouble tracegij2 = gij2[0][0] + gij2[1][1] + gij2[2][2];

    RDouble detgij = gij[0][0] * (gij[1][1] * gij[2][2] - gij[1][2] * gij[2][1])
                   - gij[0][1] * (gij[1][0] * gij[2][2] - gij[1][2] * gij[2][0])
                   + gij[0][2] * (gij[1][0] * gij[2][1] - gij[1][1] * gij[2][0]);

    RDouble t1 = tracegij;
    RDouble t2 = 0.5 * (tracegij * tracegij - tracegij2);
    RDouble t3 = detgij;

    RDouble alpha1 = t1 * t1 / 9.0 - t2 / 3.0;
    RDouble alpha2 = t1 * t1 * t1 / 27.0 - t1 * t2 / 6.0 + t3 / 2.0;
    RDouble alpha3 = 1.0 /3.0 * acos(MAX(MIN(alpha2 / sqrt(alpha1 * alpha1 * alpha1), 1.0), -1.0));

    RDouble sigma1 = sqrt(MAX(t1 / 3.0 + 2.0 * sqrt(alpha1) * cos(alpha3), 0.0));
    RDouble sigma2 = sqrt(MAX(t1 / 3.0 - 2.0 * sqrt(alpha1) * cos(acos(0.5) + alpha3), 0.0));
    RDouble sigma3 = sqrt(MAX(t1 / 3.0 - 2.0 * sqrt(alpha1) * cos(acos(0.5) - alpha3), 0.0));

    differentialSigma = sigma3 * (sigma1 - sigma2) * (sigma2 - sigma3) / sigma1 / sigma1;
    differentialSigma = MAX(differentialSigma, 0.0);

    AleModel::DestroyMatrix(velocityGradTensorTranspose);
    AleModel::DestroyMatrix(gij);
    AleModel::DestroyMatrix(gij2);
}

void LESSolverStruct::filter(const int & nn
                            , RFloat1D * fu1, RFloat1D * fu2, RFloat1D * fu3
                            , RFloat1D * fuifuj11, RFloat1D * fuifuj22, RFloat1D * fuifuj33, RFloat1D * fuifuj12, RFloat1D * fuifuj13, RFloat1D * fuifuj23
                            , RFloat1D * fs11    , RFloat1D * fs22    , RFloat1D * fs33    , RFloat1D * fs12    , RFloat1D * fs13    , RFloat1D * fs23    
                            , RFloat1D * fsfs11  , RFloat1D * fsfs22  , RFloat1D * fsfs33  , RFloat1D * fsfs12  , RFloat1D * fsfs13  , RFloat1D * fsfs23)
{
    //! Get TF(ui)
    TestFilter(nn, fu1);
    TestFilter(nn, fu2);
    TestFilter(nn, fu3);
    //! Get T(F(ui)F(uj))
    TestFilter(nn, fuifuj11);
    TestFilter(nn, fuifuj22);
    TestFilter(nn, fuifuj33);
    TestFilter(nn, fuifuj12);
    TestFilter(nn, fuifuj13);
    TestFilter(nn, fuifuj23);
    //! Get T(F(Sij))
    TestFilter(nn, fs11);
    TestFilter(nn, fs22);
    TestFilter(nn, fs33);
    TestFilter(nn, fs12);
    TestFilter(nn, fs13);
    TestFilter(nn, fs23);
    //! Get T(F(|S|)F(Sij))
    TestFilter(nn, fsfs11);
    TestFilter(nn, fsfs22);
    TestFilter(nn, fsfs33);
    TestFilter(nn, fsfs12);
    TestFilter(nn, fsfs13);
    TestFilter(nn, fsfs23);
}

void LESSolverStruct::TestFilter(const int &nn, RFloat1D *fu)
{
    RFloat1D *tfu;
    Range II(1,nn);
    tfu = new RFloat1D(II, fortranArray);
    (* tfu) = 0.0;

    Param_LESSolverStruct *parameters = GetControlParameters();

    int testFilterWidth = parameters->GetTestFilterWidth();

    RDouble * testFilterCoef = parameters->GetTestFilterCoef();

    for (int i = 1; i <= (nn-1); ++ i)
    {
        (*tfu)(i) = 0.0;
        for (int nw = -testFilterWidth; nw <= testFilterWidth; ++ nw)
        {
            RDouble coef = testFilterCoef[nw + testFilterWidth];
            RDouble f = (*fu)(i+nw);

            (*tfu)(i) = (*tfu)(i) + coef * f;
        }
    }

    for (int i = 1; i <= (nn-1); ++ i)
    {
        (*fu)(i) = (*tfu)(i);
    }

    delete tfu;
}

void LESSolverStruct::TestFilter(const int &n0, const int &nn, RDouble1D *&tf, RDouble1D *&f)
{
    Param_LESSolverStruct *parameters = GetControlParameters();

    int testFilterWidth = parameters->GetTestFilterWidth();

    RDouble * testFilterCoef = parameters->GetTestFilterCoef();

    for (int i = n0; i <= nn; ++ i)
    {
        RDouble sum = 0.0;
        for (int nw = -testFilterWidth; nw <= testFilterWidth; ++ nw)
        {
            RDouble coef = testFilterCoef[nw + testFilterWidth];
            RDouble x = (*f)(i + nw);
            sum = sum + coef * x;

            (*tf)(i) = sum;
        }
    }
}

void LESSolverStruct::TestFilter(const int &n0, const int &nn, RDouble *tf, RDouble *f)
{
    Param_LESSolverStruct *parameters = GetControlParameters();

    int testFilterWidth = parameters->GetTestFilterWidth();

    RDouble * testFilterCoef = parameters->GetTestFilterCoef();

    for (int i = n0; i <= nn; ++ i)
    {
        tf[i] = 0.0;
        for (int nw = -testFilterWidth; nw <= testFilterWidth; ++ nw)
        {
            tf[i] = tf[i] + testFilterCoef[nw + testFilterWidth] * f[i+nw];
        }
    }
}

void LESSolverStruct::GetDeltaBar(const int deltaFunctionType, const RDouble deltai, const RDouble deltaj, const RDouble deltak, RDouble &deltaBar)
{
    if (deltaFunctionType == 1)
    {
        deltaBar = MAX(MAX(deltai, deltaj), deltak);
    }
    else if (deltaFunctionType == 2)
    {
        deltaBar = pow((deltaj * deltak * deltai), third);
    }
    else if (deltaFunctionType == 3)
    {
        GetFilterWidthOfScotti(deltai, deltaj, deltak, deltaBar);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("deltaFunctionType", deltaFunctionType);
    }
}

void LESSolverStruct::GetFilterWidthOfScotti(const RDouble deltai, const RDouble deltaj, const RDouble deltak, RDouble &deltaBar)
{
    RDouble deltaMax = MAX(MAX(deltai, deltaj), deltak);
    RDouble a1 = deltai / deltaMax;
    RDouble a2 = deltaj / deltaMax;
    RDouble a3 = deltak / deltaMax;

    if (a2 < a1)
    {
        RDouble atmp = a1;
        a1 = a2;
        a2 = atmp;
    }
    if (a3 < a1)
    {
        RDouble atmp = a1;
        a1 = a3;
        a3 = atmp;
    }
    if (a3 < a2)
    {
        RDouble atmp = a2;
        a2 = a3;
        a3 = atmp;
    }

    deltaBar = cosh(sqrt(4.0 / 27.0) * (log(a1) * log(a1) - log(a1) * log(a2) + log(a2) * log(a2)));
}

void LESSolverStruct::Boundary(Grid *gridIn)
{
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (isFVMOrFDM == FD_METHOD) return;

    using namespace PHENGLEI;
    StructGrid *grid = StructGridCast(gridIn);

    StructBCSet *structBCSet = grid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);
        int BCType = structBC->GetBCType();

        if (IsInterface(BCType))
        {
            continue;
        }
        else if (BCType == EXTRAPOLATION)
        {
            OutFlowBC(grid, structBC);
        }
        else if (IsWall(BCType))
        {
            VisWall(grid, structBC);
        }
        else if (BCType == SYMMETRY)
        {
            SymmetryBC(grid, structBC);
        }
        else if (BCType == FARFIELD)
        {
            FarFieldBC(grid, structBC);
        }
        else if (BCType == INFLOW)
        {
            InFlowBC(grid, structBC);
        }
        else if (BCType == PRESSURE_INLET)
        {
            InFlowBC(grid, structBC);
        }
        else if (BCType == MASS_FLOW_INLET)
        {
            InFlowBC(grid, structBC);
        }
        else if (BCType == OUTFLOW)
        {
            OutFlowBC(grid, structBC);
        }
        else if (BCType == PRESSURE_OUTLET)
        {
            OutFlowBC(grid, structBC);
        }
        else if (BCType == MASS_FLOW_OUTLET)
        {
            OutFlowBC(grid, structBC);
        }
        else if (BCType == POLE || BCType / 10 == POLE)
        {
            OutFlowBC(grid, structBC);
        }
        else if (BCType == EXTERNAL_BC)
        {
            OutFlowBC(grid, structBC);
        }
        else
        {
            ostringstream oss;
            oss << "Error : Illegal BCtype ID " << BCType << endl;
            TK_Exit::ExceptionExit(oss);
        }
    }

    CornerPoint(grid);
}

void LESSolverStruct::ObtainBoundaryValue(Grid *gridIn)
{
    using namespace PHENGLEI;
    StructGrid *grid = StructGridCast(gridIn);

    StructBCSet *structBCSet = grid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);
        int BCType = structBC->GetBCType();

        if (IsInterface(BCType))
        {
            continue;
        }
        else if (BCType == EXTRAPOLATION)
        {
            OutFlowBC(grid, structBC);
        }
        else if (IsWall(BCType))
        {
            VisWall(grid, structBC);
        }
        else if (BCType == SYMMETRY)
        {
            SymmetryBC(grid, structBC);
        }
        else if (BCType == FARFIELD)
        {
            FarFieldBC(grid, structBC);
        }
        else if (BCType == INFLOW)
        {
            InFlowBC(grid, structBC);
        }
        else if (BCType == PRESSURE_INLET)
        {
            InFlowBC(grid, structBC);
        }
        else if (BCType == MASS_FLOW_INLET)
        {
            InFlowBC(grid, structBC);
        }
        else if (BCType == OUTFLOW)
        {
            OutFlowBC(grid, structBC);
        }
        else if (BCType == PRESSURE_OUTLET)
        {
            OutFlowBC(grid, structBC);
        }
        else if (BCType == MASS_FLOW_OUTLET)
        {
            OutFlowBC(grid, structBC);
        }
        else if (BCType == POLE || BCType / 10 == POLE)
        {
            OutFlowBC(grid, structBC);
        }

        else if (BCType == EXTERNAL_BC)
        {
            OutFlowBC(grid, structBC);
        }
        else
        {
            ostringstream oss;
            oss << "Error : Illegal BCtype ID " << BCType << endl;
            TK_Exit::ExceptionExit(oss);
        }
    }

    CornerPoint(grid);
}

void LESSolverStruct::CornerPoint(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble3D & viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));
    FillCornerPoint3D(viscousTurbulence, ni, nj, nk);
}

void LESSolverStruct::VisWall(Grid *gridIn, StructBC *structBC)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble3D &viscousTurbulence   = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble3D &subgridScaleEnergy = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("subgridScaleEnergy"));

    int ist, ied, jst, jed, kst, ked;
    structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    using namespace IDX;

    int is1, js1, ks1, is2, js2, ks2;
    int it1, jt1, kt1, it2, jt2, kt2;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                structBC->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                structBC->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                RestrictIndex(is1, js1, ks1, ni-1, nj-1, nk-1);
                RestrictIndex(is2, js2, ks2, ni-1, nj-1, nk-1);

                structBC->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                viscousTurbulence(it1, jt1, kt1) = -viscousTurbulence(is1, js1, ks1);
                viscousTurbulence(it2, jt2, kt2) = -viscousTurbulence(is2, js2, ks2);

                subgridScaleEnergy(it1, jt1, kt1) = -subgridScaleEnergy(is1, js1, ks1);
                subgridScaleEnergy(it2, jt2, kt2) = -subgridScaleEnergy(is2, js2, ks2);
            }
        }
    }
}

void LESSolverStruct::FarFieldBC(Grid *gridIn, StructBC *structBC)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &faceNormalX  = *(grid->GetFaceNormalX());
    RDouble4D &faceNormalY  = *(grid->GetFaceNormalY());
    RDouble4D &faceNormalZ  = *(grid->GetFaceNormalZ());

    Param_LESSolverStruct *parameters = GetControlParameters();

    RDouble4D &qLaminar         = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &gamma            = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("gama"));
    RDouble3D &viscousTurbulence = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble3D &subgridScaleEnergy = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("subgridScaleEnergy"));

    int nm = GlobalDataBase::GetIntParaFromDB("nm");

    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();

    using namespace IDX;

    RDouble *primitiveVariablesInfinity = reinterpret_cast<RDouble *> (GlobalDataBase::GetDataPtr("prim_inf"));
    RDouble rFarfield = primitiveVariablesInfinity[IR];
    RDouble uFarfield = primitiveVariablesInfinity[IU];
    RDouble vFarfield = primitiveVariablesInfinity[IV];
    RDouble wFarfield = primitiveVariablesInfinity[IW];
    RDouble pFarfield = primitiveVariablesInfinity[IP];

    int ist, ied, jst, jed, kst, ked;
    structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int iSurface = structBC->GetFaceDirection() + 1;
    int leftOrRightIndex    = structBC->GetFaceLeftOrRightIndex();

    RDouble *prims1 = new RDouble [nm];

    int is1, js1, ks1, is2, js2, ks2;
    int it1, jt1, kt1, it2, jt2, kt2;

    using namespace IDX;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int in, jn, kn;
                structBC->GetBoundaryFaceIndex(i, j, k, in, jn, kn);
                RDouble nxs = leftOrRightIndex * faceNormalX(in, jn, kn, iSurface);
                RDouble nys = leftOrRightIndex * faceNormalY(in, jn, kn, iSurface);
                RDouble nzs = leftOrRightIndex * faceNormalZ(in, jn, kn, iSurface);

                structBC->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                structBC->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                RestrictIndex(is1, js1, ks1, ni-1, nj-1, nk-1);
                RestrictIndex(is2, js2, ks2, ni-1, nj-1, nk-1);

                structBC->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                for (int m = 0; m < nm; ++ m)
                {
                    prims1[m] = qLaminar(is1, js1, ks1, m);
                }

                RDouble rin = prims1[IR];
                RDouble uin = prims1[IU];
                RDouble vin = prims1[IV];
                RDouble win = prims1[IW];
                RDouble pin = prims1[IP];

                RDouble vno = nxs * uFarfield + nys * vFarfield + nzs * wFarfield;
                RDouble vni = nxs * uin + nys * vin + nzs * win;
                RDouble vei = sqrt(uin * uin + vin * vin + win * win);

                RDouble gama = gamma(is1, js1, ks1);

                RDouble coo = sqrt(ABS(refGama * pFarfield / rFarfield));
                RDouble cIN = sqrt(ABS(gama  * pin / rin));

                if (vei > cIN)
                {
                    if (vni >= 0.0)
                    {
                        viscousTurbulence(it1, jt1, kt1) = viscousTurbulence(is1, js1, ks1);
                        viscousTurbulence(it2, jt2, kt2) = viscousTurbulence(is2, js2, ks2);

                        subgridScaleEnergy(it1, jt1, kt1) = subgridScaleEnergy(is1, js1, ks1);
                        subgridScaleEnergy(it2, jt2, kt2) = subgridScaleEnergy(is2, js2, ks2);
                    }
                    else
                    {
                        viscousTurbulence(it1, jt1, kt1) = freeStreamViscosity;
                        viscousTurbulence(it2, jt2, kt2) = freeStreamViscosity;

                        subgridScaleEnergy(it1, jt1, kt1) = 0.0;
                        subgridScaleEnergy(it2, jt2, kt2) = 0.0;
                    }
                    continue;
                }

                RDouble riemp = vni + 2.0 * cIN / (gama  - 1.0);
                RDouble riemm = vno - 2.0 * coo / (refGama - 1.0);
                RDouble vnb = half * (riemp + riemm);

                if (vnb  >= 0.0)
                {
                    //! Exit.
                    viscousTurbulence(it1, jt1, kt1) = viscousTurbulence(is1, js1, ks1);
                    viscousTurbulence(it2, jt2, kt2) = viscousTurbulence(it1, jt1, kt1);

                    subgridScaleEnergy(it1, jt1, kt1) = subgridScaleEnergy(is1, js1, ks1);
                    subgridScaleEnergy(it2, jt2, kt2) = subgridScaleEnergy(it1, jt1, kt1);
                }
                else
                {
                    //! Inlet.
                    viscousTurbulence(it1, jt1, kt1) = freeStreamViscosity;
                    viscousTurbulence(it2, jt2, kt2) = freeStreamViscosity;

                    subgridScaleEnergy(it1, jt1, kt1) = 0.0;
                    subgridScaleEnergy(it2, jt2, kt2) = 0.0;
                }
            }
        }
    }
    delete [] prims1;
}

void LESSolverStruct::SymmetryBC(Grid *gridIn, StructBC *structBC)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble3D &viscousTurbulence = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble3D &subgridScaleEnergy = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("subgridScaleEnergy"));

    int ist, ied, jst, jed, kst, ked;
    structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int is, js, ks, it, jt, kt;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int layer = 1; layer <= 2; ++ layer)
                {
                    structBC->GetInsideCellIndex(i, j, k, is, js, ks, layer);
                    RestrictIndex(is, js, ks, ni-1, nj-1, nk-1);
                    structBC->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                    viscousTurbulence(it, jt, kt) = viscousTurbulence(is, js, ks);
                    subgridScaleEnergy(it, jt, kt) = subgridScaleEnergy(is, js, ks);
                }
            }
        }
    }
}

void LESSolverStruct::InFlowBC(Grid *grid, StructBC *structBC)
{
    Param_LESSolverStruct *parameters = GetControlParameters();
    RDouble3D &viscousTurbulence = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble3D &subgridScaleEnergy = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("subgridScaleEnergy"));
    RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();

    int ist, ied, jst, jed, kst, ked;
    structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int it, jt, kt;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int layer = 1; layer <= 2; ++ layer)
                {
                    structBC->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                    viscousTurbulence(it, jt, kt) = freeStreamViscosity;
                    subgridScaleEnergy(it, jt, kt) = 0.0;
                }
            }
        }
    }
}

void LESSolverStruct::OutFlowBC(Grid *gridIn, StructBC *structBC)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble3D &viscousTurbulence = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble3D &subgridScaleEnergy = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("subgridScaleEnergy"));

    int ist, ied, jst, jed, kst, ked;
    structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int is1, js1, ks1, is2, js2, ks2;
    int it1, jt1, kt1, it2, jt2, kt2;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                structBC->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                structBC->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                RestrictIndex(is1, js1, ks1, ni-1, nj-1, nk-1);
                RestrictIndex(is2, js2, ks2, ni-1, nj-1, nk-1);

                structBC->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                structBC->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                viscousTurbulence(it1, jt1, kt1) = viscousTurbulence(is1, js1, ks1);
                viscousTurbulence(it2, jt2, kt2) = viscousTurbulence(it1, jt1, kt1);

                subgridScaleEnergy(it1, jt1, kt1) = subgridScaleEnergy(is1, js1, ks1);
                subgridScaleEnergy(it2, jt2, kt2) = subgridScaleEnergy(it1, jt1, kt1);
            }
        }
    }
}

}