#include "LESSolverUnstr.h"
#include "Geo_UnstructBC.h"
#include "IO_FileName.h"
#include "MultiGridOperation.h"
#include "PHMatrix.h"
#include "TK_Exit.h"
#include "Param_LESSolverUnstruct.h"
#include "Residual.h"
#include "FieldProxy.h"
#include "Math_BasisFunction.h"
#include "Constants.h"
#include "TK_Log.h"
#include "IO_FileName.h"
#include "Geo_Grid.h"

using namespace std;

namespace PHSPACE
{

LESSolverUnstr::LESSolverUnstr()
{

}

LESSolverUnstr::~LESSolverUnstr()
{

}

void LESSolverUnstr::AllocateGlobalVar(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfBoundaryFace + numberOfTotalCell;

    RDouble *cellLengthScale = new RDouble[numberOfTotalCell];    //! Length scale of cell.
    grid->UpdateDataPtr("cellLengthScale", cellLengthScale);

    RDouble **strainRateTensor = NewPointer2<RDouble>(10, numberOfTotal);    //! Strain Rate Tensor Sij and Omegaij
    grid->UpdateDataPtr("strainRateTensor", strainRateTensor);

    RDouble *wallFunction = new RDouble[numberOfTotalCell];    //! Wall damping function.
    grid->UpdateDataPtr("wallFunction", wallFunction);

    Init(grid);
}

void LESSolverUnstr::DeAllocateGlobalVar(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **strainRateTensor = reinterpret_cast<RDouble **> (grid->GetDataPtr("strainRateTensor"));
    DelPointer2(strainRateTensor);

    RDouble *cellLengthScale  = reinterpret_cast<RDouble *> (grid->GetDataPtr("cellLengthScale"));
    delete [] cellLengthScale;

    RDouble *wallFunction  = reinterpret_cast<RDouble *> (grid->GetDataPtr("wallFunction"));
    delete [] wallFunction;
}

void LESSolverUnstr::InitFlowAsRestart()
{
    UnstructGrid *grid = UnstructGridCast(GetGrid(0));

    int numberOfTotalCell    = grid->GetNTotalCell();
    Param_LESSolverUnstruct *parameters = GetControlParameters();
    RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();

    RDouble *viscousTurbulence = reinterpret_cast<RDouble  *> (grid->GetDataPtr("vist"));
    RDouble *subgridScaleEnergy = reinterpret_cast<RDouble  *> (grid->GetDataPtr("subgridScaleEnergy"));

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        viscousTurbulence[iCell] = freeStreamViscosity;
        subgridScaleEnergy[iCell] = 0.0;
    }
}

void LESSolverUnstr::Init(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    ComputeCellLengthScale(grid);
}

void LESSolverUnstr::ComputeCellLengthScale(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int numberOfTotalCell = grid->GetNTotalCell();

    RDouble *vol  = grid->GetCellVolume();
    RDouble *cellLengthScale = reinterpret_cast<RDouble *> (grid->GetDataPtr("cellLengthScale"));

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        RDouble volume = vol[iCell];

        cellLengthScale[iCell] = pow(volume, third);
    }
}

void LESSolverUnstr::Boundary(Grid *gridIn)
{
    /*
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *wallDistance = grid->GetWallDist();

    RDouble  *nxs  = grid->GetFaceNormalX();
    RDouble  *nys  = grid->GetFaceNormalY();
    RDouble  *nzs  = grid->GetFaceNormalZ();

    UnstructBCSet **bcRecord = grid->GetBCRecord();
    int numberOfBoundaryFace = grid->GetNBoundFace();

    Param_LESSolverUnstruct *parameters = GetControlParameters();
    RDouble **q      = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *gamma = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble  *> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble  *> (grid->GetDataPtr("vist"));

    RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();

    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble refGama = parameters->GetRefGama();

    using namespace IDX;
    using namespace PHENGLEI;

    RDouble *primitiveVariableInflow = reinterpret_cast<RDouble *> (GlobalDataBase::GetDataPtr("prim_inf"));
    RDouble roo = primitiveVariableInflow[IR];
    RDouble uoo = primitiveVariableInflow[IU];
    RDouble voo = primitiveVariableInflow[IV];
    RDouble woo = primitiveVariableInflow[IW];
    RDouble poo = primitiveVariableInflow[IP];

    RDouble gamm1 = refGama - 1.0;

    int count = 0;

    for (int iFace = 0; iFace < numberOfBoundaryFace; ++ iFace)
    {
        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];
        int bcType = bcRecord[iFace]->GetKey();

        if (IsInterface(bcType))
        {
            //! This statement is very important, and related to the way this part is written.
            continue;
        }
        else if (bcType == EXTRAPOLATION)
        {
            viscousTurbulence[re] = viscousTurbulence[le];
        }
        else if (IsWall(bcType))
        {
            viscousTurbulence[re] = -viscousTurbulence[le];
        }
        else if (bcType == SYMMETRY)
        {
            viscousTurbulence[re] = viscousTurbulence[le];
        }
        else if (bcType == FARFIELD)
        {
            RDouble rin,uin,vin,win,pin;
            RDouble vno,vni,vei,coo,cin;
            RDouble riemp,riemm,cb,vnb;

            rin = q[IR][le];
            uin = q[IU][le];
            vin = q[IV][le];
            win = q[IW][le];
            pin = q[IP][le];

            vno = nxs[iFace] * uoo + nys[iFace] * voo + nzs[iFace] * woo;
            vni = nxs[iFace] * uin + nys[iFace] * vin + nzs[iFace] * win;
            vei = sqrt(uin * uin + vin * vin + win * win);

            RDouble gama = gamma[le];

            coo = sqrt(ABS(refGama * poo / roo));
            cin = sqrt(ABS(gama  * pin / rin));

            //! supersonic
            if (vei > cin)
            {
                if (vni >= 0.0)
                {
                    //! exit
                    viscousTurbulence[re] = viscousTurbulence[le];
                }
                else
                {
                    //! inlet
                    viscousTurbulence[re] = freeStreamViscosity;
                }
                continue;
            }

            //! subsonic
            riemp = vni + 2.0 * cin / (gama  - 1.0);
            riemm = vno - 2.0 * coo / (refGama - 1.0);
            vnb   = half   * (riemp + riemm);
            cb    = fourth * (riemp - riemm) * gama;

            if (vnb >= 0.0)
            {
                //! exit
                viscousTurbulence[re] = viscousTurbulence[le];
            }
            else
            {
                //! inlet
                viscousTurbulence[re] = freeStreamViscosity;
            }

        }
        else if (bcType == INFLOW)
        {
            viscousTurbulence[re] = freeStreamViscosity;
        }
        else if (bcType == PRESSURE_INLET)
        {
            viscousTurbulence[re] = freeStreamViscosity;
        }
        else if (bcType == MASS_FLOW_INLET)
        {
            viscousTurbulence[re] = freeStreamViscosity;
        }
        else if (bcType == OUTFLOW)
        {
            viscousTurbulence[re] = viscousTurbulence[le];
        }
        else if (bcType == PRESSURE_OUTLET)
        {
            viscousTurbulence[re] = viscousTurbulence[le];
        }
        else if (bcType == MASS_FLOW_OUTLET)
        {
            viscousTurbulence[re] = viscousTurbulence[le];
        }
        else
        {
            viscousTurbulence[re] = viscousTurbulence[le];
        }
    }
    */
}

void LESSolverUnstr::ObtainBoundaryValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble  *nxs  = grid->GetFaceNormalX();
    RDouble  *nys  = grid->GetFaceNormalY();
    RDouble  *nzs  = grid->GetFaceNormalZ();

    UnstructBCSet **bcRecord = grid->GetBCRecord();
    int numberOfBoundaryFace = grid->GetNBoundFace();

    Param_LESSolverUnstruct *parameters = GetControlParameters();
    RDouble **q      = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *gamma = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble  *> (grid->GetDataPtr("vist"));
    RDouble *subgridScaleEnergy = reinterpret_cast<RDouble  *> (grid->GetDataPtr("subgridScaleEnergy"));

    RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();

    RDouble refGama = parameters->GetRefGama();

    using namespace IDX;
    using namespace PHENGLEI;

    RDouble *primitiveVariableInflow = reinterpret_cast<RDouble *> (GlobalDataBase::GetDataPtr("prim_inf"));
    RDouble roo = primitiveVariableInflow[IR];
    RDouble uoo = primitiveVariableInflow[IU];
    RDouble voo = primitiveVariableInflow[IV];
    RDouble woo = primitiveVariableInflow[IW];
    RDouble poo = primitiveVariableInflow[IP];

    for (int iFace = 0; iFace < numberOfBoundaryFace; ++ iFace)
    {
        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];
        int bcType = bcRecord[iFace]->GetKey();

        if (IsInterface(bcType))
        {
            //! This statement is very important, and related to the way this part is written.
            continue;
        }
        else if (bcType == EXTRAPOLATION)
        {
            viscousTurbulence[re] = viscousTurbulence[le];
            subgridScaleEnergy[re] = subgridScaleEnergy[le];
        }
        else if (IsWall(bcType))
        {
            viscousTurbulence[re] = -viscousTurbulence[le];
            subgridScaleEnergy[re] = -subgridScaleEnergy[le];
        }
        else if (bcType == SYMMETRY)
        {
            viscousTurbulence[re] = viscousTurbulence[le];
            subgridScaleEnergy[re] = subgridScaleEnergy[le];
        }
        else if (bcType == FARFIELD)
        {
            RDouble rin,uin,vin,win,pin;
            RDouble vno,vni,vei,coo,cIN;
            RDouble riemp,riemm,cb,vnb;

            rin = q[IR][le];
            uin = q[IU][le];
            vin = q[IV][le];
            win = q[IW][le];
            pin = q[IP][le];

            vno = nxs[iFace] * uoo + nys[iFace] * voo + nzs[iFace] * woo;
            vni = nxs[iFace] * uin + nys[iFace] * vin + nzs[iFace] * win;
            vei = sqrt(uin * uin + vin * vin + win * win);

            RDouble gama = gamma[le];

            coo = sqrt(ABS(refGama * poo / roo));
            cIN = sqrt(ABS(gama  * pin / rin));

            //! supersonic
            if (vei > cIN)
            {
                if (vni >= 0.0)
                {
                    //! exit
                    viscousTurbulence[re] = viscousTurbulence[le];
                    subgridScaleEnergy[re] = subgridScaleEnergy[le];
                }
                else
                {
                    //! inlet
                    viscousTurbulence[re] = freeStreamViscosity;
                    subgridScaleEnergy[re] = 0.0;
                }
                continue;
            }

            //! subsonic
            riemp = vni + 2.0 * cIN / (gama  - 1.0);
            riemm = vno - 2.0 * coo / (refGama - 1.0);
            vnb   = half   * (riemp + riemm);
            cb    = fourth * (riemp - riemm) * gama;

            if (vnb >= 0.0)
            {
                //! exit
                viscousTurbulence[re] = viscousTurbulence[le];
                subgridScaleEnergy[re] = subgridScaleEnergy[le];
            }
            else
            {
                //! inlet
                viscousTurbulence[re] = freeStreamViscosity;
                subgridScaleEnergy[re] = 0.0;
            }

        }
        else if (bcType == INFLOW)
        {
            viscousTurbulence[re] = freeStreamViscosity;
            subgridScaleEnergy[re] = 0.0;
        }
        else if (bcType == PRESSURE_INLET)
        {
            viscousTurbulence[re] = freeStreamViscosity;
            subgridScaleEnergy[re] = 0.0;
        }
        else if (bcType == MASS_FLOW_INLET)
        {
            viscousTurbulence[re] = freeStreamViscosity;
            subgridScaleEnergy[re] = 0.0;
        }
        else if (bcType == OUTFLOW)
        {
            viscousTurbulence[re] = viscousTurbulence[le];
            subgridScaleEnergy[re] = subgridScaleEnergy[le];
        }
        else if (bcType == PRESSURE_OUTLET)
        {
            viscousTurbulence[re] = viscousTurbulence[le];
            subgridScaleEnergy[re] = subgridScaleEnergy[le];
        }
        else if (bcType == MASS_FLOW_OUTLET)
        {
            viscousTurbulence[re] = viscousTurbulence[le];
            subgridScaleEnergy[re] = subgridScaleEnergy[le];
        }
        else
        {
            viscousTurbulence[re] = viscousTurbulence[le];
            subgridScaleEnergy[re] = subgridScaleEnergy[le];
        }
    }
}

void LESSolverUnstr::ComputeGradient(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    if (grid->IsFinestGrid())
    {
        GetGradientField(grid);
    }
    else
    {
        return;
    }
}

void LESSolverUnstr::GetGradientField(Grid *gridIn)
{ 
    UnstructGrid *gridUnstruct = UnstructGridCast(gridIn);
    int nEquation = GetNumberOfEquations();

    RDouble **gradPrimtiveVarX = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradPrimtiveVarY = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradPrimtiveVarZ = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarZ"));
    FieldProxy *qProxy = GetFieldProxy(gridIn, "q");
    RDouble **q     = qProxy->GetField_UNS();
    RDouble **qnode = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("qnode"));

    for (int m = 0; m < 5; ++ m)
    {
        string gradientName = GlobalDataBase::GetStrParaFromDB("gradientName");

        if (gradientName == "ggnode" || gradientName == "ggnodelaplacian")
        {
            gridUnstruct->CompGradientGGNode_NEW(q[m], qnode[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
            //gridUnstruct->CompGradientGGNode(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "ggcell")
        {
            gridUnstruct->CompGradientGGCell(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "lsq")
        {
            gridUnstruct->CompGradientLSQ(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }    
        else if (gradientName == "ggnode_weight")
        {
            gridUnstruct->CompGradientGGNodeWeight(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "gg_m2")
        {
            gridUnstruct->CompGradientGGModified2(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "ggcellnew")
        {
            gridUnstruct->CompGradientGGCellNew(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "ggcellw")
        {
            gridUnstruct->CompGradientGGCellW(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else
        {
            TK_Exit::ExceptionExit("No reconstruction method has been choosed ! /n");
        }
    }

    PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradPrimtiveVarX, "dqdx", nEquation);
    PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradPrimtiveVarY, "dqdy", nEquation);
    PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradPrimtiveVarZ, "dqdz", nEquation);

    int nTemperatureModel = 1;

    RDouble **gradTemperatureX = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTemperatureX"));
    RDouble **gradTemperatureY = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTemperatureY"));
    RDouble **gradTemperatureZ = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTemperatureZ"));
    RDouble **t     = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("t"));
    RDouble **tnode = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("tnode"));

    for (int m = 0; m < nTemperatureModel; ++ m)
    {
        string gradientName = GlobalDataBase::GetStrParaFromDB("gradientName");

        if (gradientName == "ggnode" || gradientName == "ggnodelaplacian")
        {
            gridUnstruct->CompGradientGGNode_NEW(t[m], tnode[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
            //gridUnstruct->CompGradientGGNode(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else if (gradientName == "ggcell")
        {
            gridUnstruct->CompGradientGGCell(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else if (gradientName == "lsq")
        {
            gridUnstruct->CompGradientLSQ(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }    
        else if (gradientName == "ggnode_weight")
        {
            gridUnstruct->CompGradientGGNodeWeight(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else if (gradientName == "gg_m2")
        {
            gridUnstruct->CompGradientGGModified2(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else if (gradientName == "ggcellnew")
        {
            gridUnstruct->CompGradientGGCellNew(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else if (gradientName == "ggcellw")
        {
            gridUnstruct->CompGradientGGCellW(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else
        {
            TK_Exit::ExceptionExit("No reconstruction method has been choosed ! /n");
        }
    }

    PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradTemperatureX, "dtdx", nTemperatureModel);
    PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradTemperatureY, "dtdy", nTemperatureModel);
    PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradTemperatureZ, "dtdz", nTemperatureModel);

    delete qProxy;
}

void LESSolverUnstr::ComputeViscousCoeff(Grid *gridIn)
{
    /*
    UnstructGrid *grid = UnstructGridCast(gridIn);

    ComputeStrainRateTensor(grid);

    ComputeWallFunction(grid);

    Param_LESSolverUnstruct *parameters = GetControlParameters();
    int subgridScaleModel = parameters->GetSubgridScaleModel();

    if (subgridScaleModel == 1)
    {
        Smagorinsky(grid);
    }
    else if (subgridScaleModel == 3)
    {
        WALE(grid);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("subgridScaleModel", subgridScaleModel);
    }
    */
}

void LESSolverUnstr::UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy)
{
    UnstructGrid *grid  = UnstructGridCast(gridIn);

    ObtainViscosity(grid);

    ObtainBoundaryValue(grid);
}

void LESSolverUnstr::ObtainViscosity(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    ComputeStrainRateTensor(grid);

    ComputeWallFunction(grid);

    Param_LESSolverUnstruct *parameters = GetControlParameters();
    int subgridScaleModel = parameters->GetSubgridScaleModel();

    if (subgridScaleModel == 1)
    {
        Smagorinsky(grid);
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

void LESSolverUnstr::Smagorinsky(Grid *gridIn)
{
    using namespace IDX;

    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble *cellLengthScale = reinterpret_cast<RDouble *> (grid->GetDataPtr("cellLengthScale"));
    RDouble *wallFunction = reinterpret_cast<RDouble * > (grid->GetDataPtr("wallFunction"));
    RDouble **strainRateTensor = reinterpret_cast<RDouble **> (grid->GetDataPtr("strainRateTensor"));

    RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));
    RDouble *subgridScaleEnergy = reinterpret_cast<RDouble * > (grid->GetDataPtr("subgridScaleEnergy"));

    Param_LESSolverUnstruct *parameters = GetControlParameters();

    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    RDouble refReNumber  = parameters->GetRefReNumber();

    RDouble smagConstant = parameters->GetSmagConstant();

    RDouble iConstant = parameters->GetIsotropicConstant();

    int numberOfTotalCell = grid->GetNTotalCell();

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        RDouble rho = ABS(q[IR][iCell]) + SMALL;

        RDouble deltaBar = cellLengthScale[iCell];

        RDouble dampingFunction = wallFunction[iCell];

        RDouble strain = strainRateTensor[0][iCell];

        RDouble mut = (smagConstant * deltaBar * dampingFunction) * (smagConstant * deltaBar * dampingFunction) * rho * strain * refReNumber;
        RDouble kSGS = 2.0 * dampingFunction * dampingFunction * iConstant * rho * strain * strain * refReNumber;

        viscousTurbulence[iCell] = MIN(eddyViscosityLimit, mut);
        subgridScaleEnergy[iCell] = kSGS;
    }
}

void LESSolverUnstr::WALE(Grid *gridIn)
{
    using namespace IDX;

    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble *cellLengthScale = reinterpret_cast<RDouble *> (grid->GetDataPtr("cellLengthScale"));
    RDouble *wallFunction = reinterpret_cast<RDouble * > (grid->GetDataPtr("wallFunction"));
    RDouble **strainRateTensor = reinterpret_cast<RDouble **> (grid->GetDataPtr("strainRateTensor"));

    RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));
    RDouble *subgridScaleEnergy = reinterpret_cast<RDouble * > (grid->GetDataPtr("subgridScaleEnergy"));

    Param_LESSolverUnstruct *parameters = GetControlParameters();

    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    RDouble refReNumber  = parameters->GetRefReNumber();

    RDouble waleConstant = parameters->GetWaleConstant();

    RDouble iConstant = parameters->GetIsotropicConstant();

    int numberOfTotalCell = grid->GetNTotalCell();

    RDouble viscousTurbulenceMax = 0.0;
    RDouble viscousTurbulenceMin = 1.0e30;

    int iMax = 0;
    int iMin = 0;

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        RDouble rho = ABS(q[IR][iCell]) + SMALL;

        RDouble deltaBar = cellLengthScale[iCell];

        RDouble dampingFunction = wallFunction[iCell];

        RDouble s11 = strainRateTensor[1][iCell];
        RDouble s22 = strainRateTensor[2][iCell];
        RDouble s33 = strainRateTensor[3][iCell];
        RDouble s12 = strainRateTensor[4][iCell];
        RDouble s13 = strainRateTensor[5][iCell];
        RDouble s23 = strainRateTensor[6][iCell];

        RDouble omega12 = strainRateTensor[7][iCell];
        RDouble omega13 = strainRateTensor[8][iCell];
        RDouble omega23 = strainRateTensor[9][iCell];

        RDouble sij2  = s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23);                   
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
        RDouble kSGS = 2.0 * iConstant /rho * (mut /deltaBar) * (mut /deltaBar) / refReNumber;  //! Divide by one refReNumber here, and divide by another refReNumber in the Viscousflux.

        viscousTurbulence[iCell] = mut;
        subgridScaleEnergy[iCell] = kSGS;

        if (viscousTurbulenceMax < viscousTurbulence[iCell]) 
        {
            viscousTurbulenceMax = viscousTurbulence[iCell];
            iMax = iCell;
        }
        if (viscousTurbulenceMin > viscousTurbulence[iCell])
        {
            viscousTurbulenceMin = viscousTurbulence[iCell];
            iMin = iCell;
        }
    }
}

void LESSolverUnstr::Sigma(Grid *gridIn)
{
    using namespace IDX;

    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble *cellLengthScale = reinterpret_cast<RDouble *> (grid->GetDataPtr("cellLengthScale"));
    RDouble *wallFunction = reinterpret_cast<RDouble * > (grid->GetDataPtr("wallFunction"));

    RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));
    RDouble *subgridScaleEnergy = reinterpret_cast<RDouble * > (grid->GetDataPtr("subgridScaleEnergy"));

    RDouble **gradientPrimtiveVariableX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradientPrimtiveVariableY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradientPrimtiveVariableZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble *dudx = gradientPrimtiveVariableX[1];
    RDouble *dudy = gradientPrimtiveVariableY[1];
    RDouble *dudz = gradientPrimtiveVariableZ[1];

    RDouble *dvdx = gradientPrimtiveVariableX[2];
    RDouble *dvdy = gradientPrimtiveVariableY[2];
    RDouble *dvdz = gradientPrimtiveVariableZ[2];

    RDouble *dwdx = gradientPrimtiveVariableX[3];
    RDouble *dwdy = gradientPrimtiveVariableY[3];
    RDouble *dwdz = gradientPrimtiveVariableZ[3];

    Param_LESSolverUnstruct *parameters = GetControlParameters();

    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    RDouble refReNumber  = parameters->GetRefReNumber();

    RDouble sigmaConstant = parameters->GetSigmaConstant();

    int numberOfTotalCell = grid->GetNTotalCell();

    RDouble viscousTurbulenceMax = 0.0;
    RDouble viscousTurbulenceMin = 1.0e30;

    int iMax = 0;
    int iMin = 0;

    RDouble ** velocityGradTensor = AleModel :: CreateMatrix();

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        velocityGradTensor[0][0] = dudx[iCell];
        velocityGradTensor[0][1] = dudy[iCell];
        velocityGradTensor[0][2] = dudz[iCell];

        velocityGradTensor[1][0] = dvdx[iCell];
        velocityGradTensor[1][1] = dvdy[iCell];
        velocityGradTensor[1][2] = dvdz[iCell];

        velocityGradTensor[2][0] = dwdx[iCell];
        velocityGradTensor[2][1] = dwdy[iCell];
        velocityGradTensor[2][2] = dwdz[iCell];

        RDouble differentialSigma = 0.0;
        ComputeDifferentialSigma(velocityGradTensor, differentialSigma);

        RDouble rho = ABS(q[IR][iCell]) + SMALL;

        RDouble deltaBar = cellLengthScale[iCell];

        RDouble dampingFunction = wallFunction[iCell];

        RDouble mut = rho * (sigmaConstant * deltaBar * dampingFunction) * (sigmaConstant * deltaBar * dampingFunction) * differentialSigma * refReNumber;
        mut = MIN(eddyViscosityLimit, mut);

        viscousTurbulence[iCell] = mut;
        subgridScaleEnergy[iCell] = 0.0;

        if (viscousTurbulenceMax < viscousTurbulence[iCell]) 
        {
            viscousTurbulenceMax = viscousTurbulence[iCell];
            iMax = iCell;
        }
        if (viscousTurbulenceMin > viscousTurbulence[iCell])
        {
            viscousTurbulenceMin = viscousTurbulence[iCell];
            iMin = iCell;
        }
    }
}

void LESSolverUnstr::ComputeDifferentialSigma(RDouble ** velocityGradTensor, RDouble &differentialSigma)
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

void LESSolverUnstr::ComputeStrainRateTensor(Grid *gridIn)
{
    Param_LESSolverUnstruct *parameters = GetControlParameters();
    int subgridScaleModel = parameters->GetSubgridScaleModel();
    if (subgridScaleModel == 6)
    {
        return;
    }

    using namespace IDX;

    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **strainRateTensor = reinterpret_cast<RDouble **> (grid->GetDataPtr("strainRateTensor"));

    RDouble **gradientPrimtiveVariableX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradientPrimtiveVariableY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradientPrimtiveVariableZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble *dudx = gradientPrimtiveVariableX[1];
    RDouble *dudy = gradientPrimtiveVariableY[1];
    RDouble *dudz = gradientPrimtiveVariableZ[1];

    RDouble *dvdx = gradientPrimtiveVariableX[2];
    RDouble *dvdy = gradientPrimtiveVariableY[2];
    RDouble *dvdz = gradientPrimtiveVariableZ[2];

    RDouble *dwdx = gradientPrimtiveVariableX[3];
    RDouble *dwdy = gradientPrimtiveVariableY[3];
    RDouble *dwdz = gradientPrimtiveVariableZ[3];

    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfBoundaryFace + numberOfTotalCell;

    for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
    {
        RDouble s11 = dudx[iCell];
        RDouble s22 = dvdy[iCell];
        RDouble s33 = dwdz[iCell];
        RDouble s12 = half * (dudy[iCell] + dvdx[iCell]);
        RDouble s13 = half * (dudz[iCell] + dwdx[iCell]);
        RDouble s23 = half * (dvdz[iCell] + dwdy[iCell]);

        RDouble omega12 = 0.5 * (dudy[iCell] - dvdx[iCell]);
        RDouble omega13 = 0.5 * (dudz[iCell] - dwdx[iCell]);
        RDouble omega23 = 0.5 * (dvdz[iCell] - dwdy[iCell]);

        RDouble sij2 = two * (PHSPACE::SQR(s11, s22, s33) + two * PHSPACE::SQR(s12, s13, s23));
        RDouble ss = sqrt(sij2);

        strainRateTensor[0][iCell] = ss;

        strainRateTensor[1][iCell] = s11;
        strainRateTensor[2][iCell] = s22;
        strainRateTensor[3][iCell] = s33;
        strainRateTensor[4][iCell] = s12;
        strainRateTensor[5][iCell] = s13;
        strainRateTensor[6][iCell] = s23;

        strainRateTensor[7][iCell] = omega12;
        strainRateTensor[8][iCell] = omega13;
        strainRateTensor[9][iCell] = omega23;
    }
}

void LESSolverUnstr::ComputeWallFunction(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    
    Param_LESSolverUnstruct *parameters = GetControlParameters();
    int wallDampingFunctionType = parameters->GetWallDampingFunctionType();

    RDouble *wallFunction = reinterpret_cast<RDouble * > (grid->GetDataPtr("wallFunction"));

    if (wallDampingFunctionType == 0)
    {
        int numberOfTotalCell = grid->GetNTotalCell();

        for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
        {
            wallFunction[iCell] = 1.0;
        }
    }
    else if (wallDampingFunctionType == 1)
    {
        WallFunctionofVanDriest(grid);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("wallDampingFunctionType", wallDampingFunctionType);
    }
}

void LESSolverUnstr::WallFunctionofVanDriest(Grid *gridIn)
{
    using namespace IDX;

    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble *wallDistance = grid->GetWallDist();

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble *visl = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));

    RDouble *wallFunction = reinterpret_cast <RDouble *> (grid->GetDataPtr("wallFunction"));

    RDouble **strainRateTensor = reinterpret_cast<RDouble **> (grid->GetDataPtr("strainRateTensor"));

    Param_LESSolverUnstruct *parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();

    const double kappa = 0.41;
    const double Ap = 25.0;

    int numberOfTotalCell = grid->GetNTotalCell();

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        RDouble wallDist = wallDistance[iCell];

        RDouble mul = visl[iCell];

        RDouble rho = q[IR][iCell];

        RDouble ss = strainRateTensor[0][iCell];

        RDouble yplus = ss * (kappa * wallDist) * (kappa * wallDist) * rho / mul * refReNumber;

        wallFunction[iCell] = 1.0 - exp(-yplus / Ap);
    }
}

void LESSolverUnstr::UploadInterfaceData(ActionKey *actkey)
{
    UploadInterfaceValue(actkey);
}

void LESSolverUnstr::UploadInterfaceValue(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (interfaceInformation == 0) return;

    RDouble *viscousTurbulence = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    PHSPACE::UploadInterfaceValue(grid, viscousTurbulence, "vist");

    RDouble *subgridScaleEnergy = reinterpret_cast<RDouble *>(grid->GetDataPtr("subgridScaleEnergy"));
    PHSPACE::UploadInterfaceValue(grid, subgridScaleEnergy, "subgridScaleEnergy");
}

void LESSolverUnstr::DownloadInterfaceData(ActionKey *actkey)
{
    DownloadInterfaceValue(actkey);
}

void LESSolverUnstr::DownloadInterfaceValue(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (interfaceInformation == 0) return;

    RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
    PHSPACE::DownloadInterfaceValue(grid, viscousTurbulence, "vist");

    RDouble *subgridScaleEnergy = reinterpret_cast<RDouble *> (grid->GetDataPtr("subgridScaleEnergy"));
    PHSPACE::DownloadInterfaceValue(grid, subgridScaleEnergy, "subgridScaleEnergy");
}

LIB_EXPORT void LESSolverUnstr::InitControlParameters()
{
    FreeControlParameters();
    controlParameters = new Param_LESSolverUnstruct();
    controlParameters->Init();
}

LIB_EXPORT Param_LESSolverUnstruct * LESSolverUnstr::GetControlParameters()
{
    return static_cast<Param_LESSolverUnstruct *>(controlParameters);
}

FieldProxy * LESSolverUnstr::GetFieldProxy(Grid *gridIn, const string &field_name)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble ** field = reinterpret_cast< RDouble ** > (grid->GetDataPtr(field_name));

    FieldProxy * field_proxy = new FieldProxy();

    field_proxy->SetField_UNS(field);

    return field_proxy;
}

}