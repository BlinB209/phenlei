#include "PerfectGas.h"
#include "Constants.h"
#include "Geo_SimpleBC.h"
#include "GlobalDataBase.h"
#include "TK_Exit.h"
#include "PHIO.h"
#include "Glb_Dimension.h"
#include "TK_Log.h"
#include "Math_BasisFunction.h"
#include "IO_FileName.h"
#include "TK_Parse.h"

namespace PHSPACE
{
namespace GAS_SPACE
{

PerfectGas::PerfectGas()
{
    nGasModel = 0;
    if (GlobalDataBase::IsExist("nGasModel", PHINT, 1))
    {
        nGasModel = GlobalDataBase::GetIntParaFromDB("nGasModel");
    }
    molecularDiameter = 3.6e-10; //m
}

PerfectGas::~PerfectGas()
{
    RDouble *prim_inf = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("prim_inf"));
    delete [] prim_inf;    prim_inf = NULL;
    GlobalDataBase::DeleteDataPtr("prim_inf");

    //Restore the following variables.
    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");
    RDouble forceReferenceLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
    RDouble forceReferenceLengthSpanWise = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");
    RDouble forceReferenceArea = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");

    RDouble TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    RDouble TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    RDouble TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

    if (GetDim() == THREE_D)
    {
        forceReferenceArea *= reynoldsReferenceLengthDimensional*reynoldsReferenceLengthDimensional;
    }
    else if (GetDim() == TWO_D)
    {
        forceReferenceArea *= reynoldsReferenceLengthDimensional;
    }
    forceReferenceLength *= reynoldsReferenceLengthDimensional;
    forceReferenceLengthSpanWise *= reynoldsReferenceLengthDimensional;
    TorqueRefX *= reynoldsReferenceLengthDimensional;
    TorqueRefY *= reynoldsReferenceLengthDimensional;
    TorqueRefZ *= reynoldsReferenceLengthDimensional;

    GlobalDataBase::UpdateData("forceReferenceLength", &forceReferenceLength, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("forceReferenceLengthSpanWise", &forceReferenceLengthSpanWise, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("forceReferenceArea", &forceReferenceArea, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("TorqueRefZ", &TorqueRefZ, PHDOUBLE, 1);

    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    refReNumber /= reynoldsReferenceLengthDimensional;
    GlobalDataBase::UpdateData("refReNumber", &refReNumber, PHDOUBLE, 1);
}
#ifdef USE_GMRESSOLVER
// GMRESPV
void PerfectGas::dPrimitive2dConservative(RDouble *prim, RDouble gama, RDouble** dqdcv)
{
    using namespace IDX;

    RDouble& r     = prim[IR];
    RDouble& u     = prim[IU];
    RDouble& v     = prim[IV];
    RDouble& w     = prim[IW];
    RDouble& p     = prim[IP];

    dqdcv[IR][IR] = 1.0;
    dqdcv[IU][IR] = -1.0*u/r;
    dqdcv[IV][IR] = -1.0*v/r;
    dqdcv[IW][IR] = -1.0*w/r;
    dqdcv[IP][IR] = 0.5*(gama-1.0)*(u*u + v*v + w*w);

    dqdcv[IU][IU] = 1.0/r;
    dqdcv[IP][IU] = -1.0*(gama-1)*u;

    dqdcv[IV][IV] = 1.0/r;
    dqdcv[IP][IV] = -1.0*(gama-1)*v;

    dqdcv[IW][IW] = 1.0/r;
    dqdcv[IP][IW] = -1.0*(gama-1)*w;

    dqdcv[IP][IP] = gama-1;


}

// GMRESVis
// the derivative of primitive variables w.r.t conservative variables
void PerfectGas::dGradient2dPrimitive(RDouble *prim, int sign, RDouble** dgraddq, char direction, RDouble nxs, RDouble nys, RDouble nzs, RDouble ns, RDouble vol, int nPara)
{
    using namespace IDX;
    RDouble& r     = prim[IR];
    RDouble& p     = prim[IP];

    if(direction == 'x')
    {
        dgraddq[IR][IR] = sign*0.5/vol*ns*nxs;
        dgraddq[IU][IR] = 0;
        dgraddq[IV][IR] = 0;
        dgraddq[IW][IR] = 0;
        dgraddq[IP][IR] = 0;

        dgraddq[IR][IU] = 0;
        dgraddq[IU][IU] = sign*0.5/vol*ns*nxs;
        dgraddq[IV][IU] = 0;
        dgraddq[IW][IU] = 0;
        dgraddq[IP][IU] = 0;

        dgraddq[IR][IV] = 0;
        dgraddq[IU][IV] = 0;
        dgraddq[IV][IV] = sign*0.5/vol*ns*nxs;
        dgraddq[IW][IV] = 0;
        dgraddq[IP][IV] = 0;

        dgraddq[IR][IW] = 0;
        dgraddq[IU][IW] = 0;
        dgraddq[IV][IW] = 0;
        dgraddq[IW][IW] = sign*0.5/vol*ns*nxs;
        dgraddq[IP][IW] = 0;

        dgraddq[IR][IP] = 0;
        dgraddq[IU][IP] = 0;
        dgraddq[IV][IP] = 0;
        dgraddq[IW][IP] = 0;
        dgraddq[IP][IP] = sign*0.5/vol*ns*nxs;

        if( nPara == 6)
        {
            dgraddq[IP+1][IR] = sign*0.5/vol*ns*nxs*(-1.0*p/(coefficientOfStateEquation*r*r));
            dgraddq[IP+1][IU] = 0;
            dgraddq[IP+1][IV] = 0;
            dgraddq[IP+1][IW] = 0;
            dgraddq[IP+1][IP] = sign*0.5/vol*ns*nxs*1.0/(coefficientOfStateEquation*r);
        }
    }
    else if(direction == 'y')
    {
        dgraddq[IR][IR] = sign*0.5/vol*ns*nys;
        dgraddq[IU][IR] = 0;
        dgraddq[IV][IR] = 0;
        dgraddq[IW][IR] = 0;
        dgraddq[IP][IR] = 0;

        dgraddq[IR][IU] = 0;
        dgraddq[IU][IU] = sign*0.5/vol*ns*nys;
        dgraddq[IV][IU] = 0;
        dgraddq[IW][IU] = 0;
        dgraddq[IP][IU] = 0;

        dgraddq[IR][IV] = 0;
        dgraddq[IU][IV] = 0;
        dgraddq[IV][IV] = sign*0.5/vol*ns*nys;
        dgraddq[IW][IV] = 0;
        dgraddq[IP][IV] = 0;

        dgraddq[IR][IW] = 0;
        dgraddq[IU][IW] = 0;
        dgraddq[IV][IW] = 0;
        dgraddq[IW][IW] = sign*0.5/vol*ns*nys;
        dgraddq[IP][IW] = 0;

        dgraddq[IR][IP] = 0;
        dgraddq[IU][IP] = 0;
        dgraddq[IV][IP] = 0;
        dgraddq[IW][IP] = 0;
        dgraddq[IP][IP] = sign*0.5/vol*ns*nys;


        if( nPara == 6)
        {
            dgraddq[IP+1][IR] = sign*0.5/vol*ns*nys*(-1.0*p/(coefficientOfStateEquation*r*r));
            dgraddq[IP+1][IU] = 0;
            dgraddq[IP+1][IV] = 0;
            dgraddq[IP+1][IW] = 0;
            dgraddq[IP+1][IP] = sign*0.5/vol*ns*nys*1.0/(coefficientOfStateEquation*r);
        }

    }
    else if(direction == 'z')
    {
        dgraddq[IR][IR] = sign*0.5/vol*ns*nzs;
        dgraddq[IU][IR] = 0;
        dgraddq[IV][IR] = 0;
        dgraddq[IW][IR] = 0;
        dgraddq[IP][IR] = 0;

        dgraddq[IR][IU] = 0;
        dgraddq[IU][IU] = sign*0.5/vol*ns*nzs;
        dgraddq[IV][IU] = 0;
        dgraddq[IW][IU] = 0;
        dgraddq[IP][IU] = 0;

        dgraddq[IR][IV] = 0;
        dgraddq[IU][IV] = 0;
        dgraddq[IV][IV] = sign*0.5/vol*ns*nzs;
        dgraddq[IW][IV] = 0;
        dgraddq[IP][IV] = 0;

        dgraddq[IR][IW] = 0;
        dgraddq[IU][IW] = 0;
        dgraddq[IV][IW] = 0;
        dgraddq[IW][IW] = sign*0.5/vol*ns*nzs;
        dgraddq[IP][IW] = 0;

        dgraddq[IR][IP] = 0;
        dgraddq[IU][IP] = 0;
        dgraddq[IV][IP] = 0;
        dgraddq[IW][IP] = 0;
        dgraddq[IP][IP] = sign*0.5/vol*ns*nzs;

        if( nPara == 6)
        {
            dgraddq[IP+1][IR] = sign*0.5/vol*ns*nzs*(-1.0*p/(coefficientOfStateEquation*r*r));
            dgraddq[IP+1][IU] = 0;
            dgraddq[IP+1][IV] = 0;
            dgraddq[IP+1][IW] = 0;
            dgraddq[IP+1][IP] = sign*0.5/vol*ns*nzs*1.0/(coefficientOfStateEquation*r);
            
        }
    }

}

// GMRESnolim
void PerfectGas::dGradient2dPrimitive4limiter(int sign, RDouble** dgraddq, char direction, RDouble nxs, RDouble nys, RDouble nzs, RDouble ns, RDouble vol)
{
    using namespace IDX;

    if(direction == 'x')
    {
        dgraddq[IR][IR] = sign*0.5/vol*ns*nxs;
       
        dgraddq[IU][IU] = sign*0.5/vol*ns*nxs;
        
        dgraddq[IV][IV] = sign*0.5/vol*ns*nxs;
     
        dgraddq[IW][IW] = sign*0.5/vol*ns*nxs;
    
        dgraddq[IP][IP] = sign*0.5/vol*ns*nxs;
    }
    else if(direction == 'y')
    {
        dgraddq[IR][IR] = sign*0.5/vol*ns*nys;
  
        dgraddq[IU][IU] = sign*0.5/vol*ns*nys;
     
        dgraddq[IV][IV] = sign*0.5/vol*ns*nys;

        dgraddq[IW][IW] = sign*0.5/vol*ns*nys;

        dgraddq[IP][IP] = sign*0.5/vol*ns*nys;

    }
    else if(direction == 'z')
    {
        dgraddq[IR][IR] = sign*0.5/vol*ns*nzs;
  
        dgraddq[IU][IU] = sign*0.5/vol*ns*nzs;
       
        dgraddq[IV][IV] = sign*0.5/vol*ns*nzs;
        
        dgraddq[IW][IW] = sign*0.5/vol*ns*nzs;
        
        dgraddq[IP][IP] = sign*0.5/vol*ns*nzs;
    }

}  
#endif

void PerfectGas::Primitive2Conservative(RDouble *prim, RDouble gama, RDouble Tv, RDouble Te, RDouble *q)
{
    using namespace IDX;
    RDouble em;

    //! obtain the total internal energy Em.
    ComputeInternalEnergy(prim, gama, Tv, Te, em);

    q[IR ] = prim[IR];
    q[IRU] = prim[IR] * prim[IU];
    q[IRV] = prim[IR] * prim[IV];
    q[IRW] = prim[IR] * prim[IW];
    q[IRE] = prim[IR] * em;

    int isPorousZone = GlobalDataBase::GetIntParaFromDB("isPorousZone");
    if (isPorousZone != NON_POROUS)
    {
        RDouble porosity = GlobalDataBase::GetDoubleParaFromDB("porosity");
        RDouble densitySolid = GlobalDataBase::GetDoubleParaFromDB("densitySolid");
        RDouble cpSolid = GlobalDataBase::GetDoubleParaFromDB("cpSolid");
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble refVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
        RDouble refTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble referenceEnergy = refVelocity * refVelocity;
        RDouble tm = gama * refMachNumber * refMachNumber * prim[IP] / prim[IR];
        q[IRE] = q[IRE] * porosity + densitySolid * cpSolid * tm * refTemperature / referenceEnergy * (1.0 - porosity);
    }
}

void PerfectGas::ComputeInternalEnergy(RDouble *primitiveVariables, RDouble gama, RDouble vibrationTemperature, RDouble electronTemperature, RDouble &totalEnergy)
{
    using namespace IDX;

    RDouble &um = primitiveVariables[IU];
    RDouble &vm = primitiveVariables[IV];
    RDouble &wm = primitiveVariables[IW];
    RDouble &density = primitiveVariables[IR];
    RDouble &pressure = primitiveVariables[IP];
    RDouble squareVelocity = um * um + vm * vm + wm * wm;

    totalEnergy = (pressure / density) / (gama - one) + half * squareVelocity;
}

void PerfectGas::Conservative2Primitive(RDouble *q, RDouble gama, RDouble *prim, RDouble *temperature)
{
    using namespace IDX;
    RDouble density, oDensity, um, vm, wm, internalEnergy;
    density = q[IR];
    if (density <= zero)
    {
        prim[IR] = density;
        return;
    }

    oDensity = 1.0 / density;
    um = q[IRU] * oDensity;
    vm = q[IRV] * oDensity;
    wm = q[IRW] * oDensity;
    internalEnergy = q[IRE] - half * density * (um * um + vm * vm + wm * wm);

    prim[IR] = density;
    prim[IU] = um;
    prim[IV] = vm;
    prim[IW] = wm;
    prim[IP] = (gama - one) * internalEnergy;

    //! Compute temperature.
    RDouble gasConstant = this->GetCoefficientOfStateEquation();
    temperature[ITT] = prim[IP] / (prim[IR] * gasConstant);
    temperature[ITV] = temperature[ITT];
    temperature[ITE] = temperature[ITT];
}

void PerfectGas::GetSpecificHeatRatio(RDouble *prim, RDouble &gama)
{
    gama = referenceParameterDimensional.GetAverageSpecificHeatRatio();
}

void PerfectGas::ComputeEnthalpyByPrimitive(RDouble *primitiveVariables, RDouble &gama, RDouble &enthalpy, RDouble* temperatures)
{
    using namespace IDX;
    enthalpy = (gama / (gama - one)) * (primitiveVariables[IP] / primitiveVariables[IR]);
}

void PerfectGas::ComputeTotalEnthalpyAndDH(RDouble *primitiveVariables, RDouble &gama, RDouble *deltaQ, RDouble &totalEnthalpy, RDouble &deltaEnthalpy)
{
    //! This function is used for perfect gas and single temperature model.
    //! To compute the total enthalpy and the variable dh in the function MXDQ_STD(). dh=b2, b2 denotes the coefficient \n
    //! of the vector M*dQ which can be referred to the forumla (A.7) and (A.8) in the appendix A of the PHengLEI Theory manual.
    using namespace IDX;

    //! Obtain the primitive variables.
    RDouble &um = primitiveVariables[IU];
    RDouble &vm = primitiveVariables[IV];
    RDouble &wm = primitiveVariables[IW];
    RDouble &density = primitiveVariables[IR];
    RDouble &pressure = primitiveVariables[IP];

    //! c1 is product coefficient of the variable dq[IR]£¬gama1 = gama-1¡£
    RDouble squareVelocity, c1, gama1;
    squareVelocity = um * um + vm * vm + wm * wm;
    gama1 = gama - 1.0;
    c1 = half * gama1 * squareVelocity;
    deltaEnthalpy = -gama1 * (um * deltaQ[IRU] + vm * deltaQ[IRV] + wm * deltaQ[IRW] - deltaQ[IRE]);

    //! enthalpy is static enthalpy.
    RDouble enthalpy = (gama / gama1)* (pressure / density);

    deltaEnthalpy += c1 * deltaQ[IR];
    //! h0 = h+v*v/2, total enthalpy equals to static enthalpy plus half of the squared velocity.
    totalEnthalpy = enthalpy + half * squareVelocity;
}

void PerfectGas::ComputeDHAndTotalEnthalpy(RDouble *primitiveVariables, RDouble &gama, RDouble *deltaQ, RDouble &deltaEnthalpy, RDouble &totalEnthalpy)
{
    using namespace IDX;

    RDouble squareVelocity, ae, af;
    RDouble enthalpy;
    //! Obtain the primitive variables.
    RDouble &um = primitiveVariables[IU];
    RDouble &vm = primitiveVariables[IV];
    RDouble &wm = primitiveVariables[IW];
    RDouble &density = primitiveVariables[IR];
    RDouble &pressure = primitiveVariables[IP];

    squareVelocity = um * um + vm * vm + wm * wm;
    ae = gama - one;
    af = half * ae * squareVelocity;
    deltaEnthalpy = -ae * (um * deltaQ[IRU] + vm * deltaQ[IRV] + wm * deltaQ[IRW] - deltaQ[IRE]);
    enthalpy = (gama / (gama - one)) * (pressure / density);
  
    deltaEnthalpy += af * deltaQ[IR];
    // ht = h+v*v/2
    totalEnthalpy = enthalpy + half * squareVelocity;
}

RDouble PerfectGas::GetMixedGasEnthalpy(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature)
{
    RDouble gamma = referenceParameterDimensional.GetAverageSpecificHeatRatio();
    RDouble  staticEnthaly = (gamma / (gamma - 1.0)) * primitiveVars[IDX::IP] / primitiveVars[IDX::IR];
    return staticEnthaly;
}

void PerfectGas::ComputeViscosityByPrimitive(RDouble *primitiveVariables, RDouble &transRotationTemperature, RDouble &electronTemperature, RDouble &suthTemperature, RDouble &viscosity)
{
    if (this->nGasModel > 1)    //! The one molecule of the gas.
    {
        //! Referred to the "2008-Lofthouse-Nonequilibrium Hypersonic Aerothermodynamics Using the Direct Simulation Monte Carlo And Navier-Stokes Models".
        RDouble T = transRotationTemperature * referenceParameterDimensional.GetTemperature();
        //! The default is the Argon.
        RDouble Omega = 0.734, Tref = 1000.0, Dref = 3.595e-10;
        if (this->nGasModel == 3)    //! The pure gas of Nitrogen.
        {
            Omega = 0.7;
            Tref = 290.0;
            Dref = 4.11e-10;
        }
        RDouble mass = referenceParameterDimensional.GetAverageMolecularWeight() / AVOGADRO_CONSTANT;
        RDouble refViscosity = 15 * sqrt(PI * mass * BOLTZMANN_CONSTANT * Tref);
        refViscosity = refViscosity / (2 * PI * Dref * Dref) / (5 - 2 * Omega) / (7 - 2 * Omega);

        viscosity = refViscosity * pow(T / Tref, Omega);
        viscosity = viscosity / referenceParameterDimensional.GetViscosity();
    }
    else
    {
        viscosity = transRotationTemperature * sqrt(transRotationTemperature) * (1.0 + suthTemperature) / (transRotationTemperature + suthTemperature);
    }
}

void PerfectGas::InitGlobalParameterOfNSEquation()
{
    //! Obtain the variable of "reynoldsReferenceLengthDimensional" according to "gridScaleFactor".
    InitCommonParameter();

    InitReferenceParameter();

    ComputeReferenceParameter();

    PrintInflowConditionToWindow();

    ComputeProperReynoldsNumberForGrid();

    //! The five reference variables for the whole aerodynamic coefficients integration, which are from the file of "cfd_para.hypara", should be updated with "reynoldsReferenceLengthDimensional".
    ComputeOtherProperParameterForGrid();

    InitOtherParameter();

    GlobalBoundaryCondition::SetGlobalBCByGlobalDataBase();

    GlobalBoundaryCondition::InitMassFlowBoundary();

    GlobalBoundaryCondition::SetBCDataBaseByGlobalBC();
}

void PerfectGas::InitNumberOfSpecies()
{
    int numberOfSpecies = 0;
    GlobalDataBase::UpdateData("numberOfSpecies", &numberOfSpecies, PHINT, 1);

    int nm = GlobalDataBase::GetIntParaFromDB("nm");
    int nl = nm;
    GlobalDataBase::UpdateData("nl", &nl, PHINT, 1);
}

void PerfectGas::ComputeProperReynoldsNumberForGrid()
{
    RDouble refReNumber = referenceParameterDimensional.GetReynoldsNumber();
    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

    refReNumber *= reynoldsReferenceLengthDimensional;
    RDouble oreynolds = 1.0 / refReNumber;

    referenceParameterDimensional.SetReynoldsNumber(refReNumber);
    GlobalDataBase::UpdateData("refReNumber", &refReNumber, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("oreynolds", &oreynolds, PHDOUBLE, 1);
}

void PerfectGas::InitCommonParameter()
{
    RDouble reynoldsReferenceLengthDimensional = 1.0;

    RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

    reynoldsReferenceLengthDimensional = gridScaleFactor;

    GlobalDataBase::UpdateData("reynoldsReferenceLengthDimensional", &reynoldsReferenceLengthDimensional, PHDOUBLE, 1);

    //!for AleManager classs variables
    GlobalDataBase::UpdateData("referenceLength", &reynoldsReferenceLengthDimensional, PHDOUBLE, 1);

    RDouble prl = GlobalDataBase::GetDoubleParaFromDB("prl");
    RDouble prt = GlobalDataBase::GetDoubleParaFromDB("prt");
    RDouble gamma = GlobalDataBase::GetDoubleParaFromDB("refGama");

    if (this->nGasModel == 2)    //! The Argon gas.
    {
        //! Referred to the "2008-Lofthouse-Nonequilibrium Hypersonic Aerothermodynamics Using the Direct Simulation Monte Carlo And Navier-Stokes Models".
        gamma = 5.0 / 3.0;
        GlobalDataBase::UpdateData("refGama", &gamma, PHDOUBLE, 1);
        prl = 4.0 * gamma / (9.0 * gamma - 5.0);
        GlobalDataBase::UpdateData("prl", &prl, PHDOUBLE, 1);
        GlobalDataBase::UpdateData("refGama", &gamma, PHDOUBLE, 1);
    }
    else if (this->nGasModel == 3)    //! The Nitrogen gas.
    {
        gamma = 1.4;
        prl = 4.0 * gamma / (9.0 * gamma - 5.0);
        GlobalDataBase::UpdateData("prl", &prl, PHDOUBLE, 1);
        GlobalDataBase::UpdateData("refGama", &gamma, PHDOUBLE, 1);
    }

    RDouble oprl = 1.0 / prl;
    GlobalDataBase::UpdateData("oprl", &oprl, PHDOUBLE, 1);
    RDouble oprt = 1.0 / prt;
    GlobalDataBase::UpdateData("oprt", &oprt, PHDOUBLE, 1);
}

void PerfectGas::InitOtherParameter()
{
    int nolstress, nrokplus;
    nolstress = -1;
    nrokplus  = -1;
    GlobalDataBase::UpdateData("nolstress", &nolstress, PHINT, 1);
    GlobalDataBase::UpdateData("nrokplus" , &nrokplus , PHINT, 1);

    int outnstep = 0;
    GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);

    int innstep = 0;
    GlobalDataBase::UpdateData("innstep", &innstep, PHINT, 1);

    int subTotalIterStep = 0;
    GlobalDataBase::UpdateData("subTotalIterStep", &subTotalIterStep, PHINT, 1);

    RDouble physicalTime = 0.0;
    GlobalDataBase::UpdateData("physicalTime", &physicalTime, PHDOUBLE, 1);

    int nIndexTT = 0;
    int nIndexTV = 0;
    int nIndexTE = 0;
    GlobalDataBase::UpdateData("indexOfTT", &nIndexTT, PHINT, 1);
    GlobalDataBase::UpdateData("indexOfTV", &nIndexTV, PHINT, 1);
    GlobalDataBase::UpdateData("indexOfTE", &nIndexTE, PHINT, 1);
}

void PerfectGas::InitReferenceParameter()
{
    InitNumberOfSpecies();

    int inflowParaType = GlobalDataBase::GetIntParaFromDB("inflowParaType");

    if (inflowParaType == NONDIMENSIONCONDITION)
    {
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble refReNumber   = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

        referenceParameterDimensional.SetMachNumber(refMachNumber);
        referenceParameterDimensional.SetReynoldsNumber(refReNumber);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
    }
    else if (inflowParaType == FLIGHTCONDITION)
    {
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");

        referenceParameterDimensional.SetMachNumber(refMachNumber);

        SetAirInformationByDataBase();
    }
    else if (inflowParaType == EXPERIMENTCONDITION)
    {
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble refDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");

        referenceParameterDimensional.SetMachNumber(refMachNumber);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
        referenceParameterDimensional.SetPressure(refDimensionalPressure);

        SetAirInformationByExpData();
    }
    else if (inflowParaType == TEMPERATURE_DENSITY)    //! where temperature and density are fixed.
    {
        RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");

        referenceParameterDimensional.SetVelocity(refDimensionalVelocity);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
        referenceParameterDimensional.SetDensity(refDimensionalDensity);
    }
    else if (inflowParaType == TEMPERATURE_PRESSURE || inflowParaType == WINDSPEEDPROFILE)    //! where temperature and pressure are fixed.
    {
        RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble refDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");

        referenceParameterDimensional.SetVelocity(refDimensionalVelocity);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
        referenceParameterDimensional.SetPressure(refDimensionalPressure);
    }
    else if (inflowParaType == MACH_TEMP_PRE)    //! where temperature and pressure are fixed.
    {
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble refDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");

        referenceParameterDimensional.SetMachNumber(refMachNumber);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
        referenceParameterDimensional.SetPressure(refDimensionalPressure);
    }
    else if (inflowParaType == WEATHERCONDITION)
    {
        SetAirInformationByWRFDataBase();
    }
    else
    {
        TK_Exit::UnexpectedVarValue("inflowParaType", inflowParaType);
    }
}

void PerfectGas::SetAirInformationByDataBase()
{
    RDouble height = GlobalDataBase::GetDoubleParaFromDB("height");

    RDouble refDimensionalTemperature = 1.0;
    RDouble refDimensionalPressure = 1.0;
    RDouble refDimensionalDensity = 1.0;
    RDouble refDimensionalSonicSpeed=1.0;

    GetAirInfo(height, refDimensionalTemperature, refDimensionalPressure, refDimensionalDensity, refDimensionalSonicSpeed);
    
    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    referenceParameterDimensional.SetPressure(refDimensionalPressure);
    referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
    referenceParameterDimensional.SetSoundVelocity(refDimensionalSonicSpeed);

    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalTemperature", &refDimensionalTemperature, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalSonicSpeed", &refDimensionalSonicSpeed, PHDOUBLE, 1);
}

void PerfectGas::SetAirInformationByExpData()
{
    RDouble refDimensionalTemperature;
    RDouble refDimensionalPressure;
    RDouble refDimensionalDensity;
    RDouble refDimensionalSonicSpeed;

    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    RDouble T0 = referenceParameterDimensional.GetTemperature();
    RDouble P0 = referenceParameterDimensional.GetPressure();

    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble refMolecularWeight = GlobalDataBase::GetDoubleParaFromDB("refMolecularWeight");
    int nGasModel = GlobalDataBase::GetIntParaFromDB("nGasModel");

    RDouble R = 287.053;
    if(nGasModel > 0)
    {
        R = 1000.0 * GeneralGasConstant / refMolecularWeight;
    }

    RDouble mi   = refGama / (refGama - 1.0);
    RDouble coef = 1.0 + (refGama - 1.0) * refMachNumber * refMachNumber / 2.0;

    RDouble T, P, rho, c;
    T = T0 / coef;
    P = P0 / pow(coef, mi);
    rho = P / R / T;
    c = sqrt(refGama * R * T);

    refDimensionalDensity = rho;
    refDimensionalPressure = P;
    refDimensionalTemperature = T;
    refDimensionalSonicSpeed = c;

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    referenceParameterDimensional.SetPressure(refDimensionalPressure);
    referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
    referenceParameterDimensional.SetSoundVelocity(refDimensionalSonicSpeed);

    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalTemperature", &refDimensionalTemperature, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalSonicSpeed", &refDimensionalSonicSpeed, PHDOUBLE, 1);
}

void PerfectGas::SetAirInformationByExpData(int useSetting)
{
    RDouble refDimensionalTemperature;
    RDouble refDimensionalPressure;
    RDouble refDimensionalDensity;
    RDouble refDimensionalSonicSpeed;

    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    RDouble T0 = referenceParameterDimensional.GetTemperature();
    RDouble P0 = referenceParameterDimensional.GetPressure();

    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble refMolecularWeight = GlobalDataBase::GetDoubleParaFromDB("refMolecularWeight");
    int nGasModel = GlobalDataBase::GetIntParaFromDB("nGasModel");

    RDouble R = 287.053;
    RDouble mi   = refGama / (refGama - 1.0);
    RDouble coef = 1.0 + (refGama - 1.0) * refMachNumber * refMachNumber / 2.0;

    RDouble T, P, rho, c;
    T = T0 / coef;
    P = P0 / pow(coef, mi);
    rho = P / R / T;
    c = sqrt(refGama * R * T);

    refDimensionalDensity = rho;
    refDimensionalPressure = P;
    refDimensionalTemperature = T;
    refDimensionalSonicSpeed = c;

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    referenceParameterDimensional.SetPressure(refDimensionalPressure);
    referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
    referenceParameterDimensional.SetSoundVelocity(refDimensionalSonicSpeed);

    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalTemperature", &refDimensionalTemperature, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalSonicSpeed", &refDimensionalSonicSpeed, PHDOUBLE, 1);
}

void PerfectGas::GetAirInfo(const RDouble &h1, RDouble &t, RDouble &p, RDouble &den, RDouble &a)
{
    //USSA 76 Model for atmosphere.
    RDouble h      = h1 * 1000.0;
    RDouble r      = 287.053;
    RDouble g0     = 9.80665;
    RDouble rp     = 6.37111e6;
    RDouble g      = (rp/(rp+h))*(rp/(rp+h))*g0;
    RDouble t0     = 288.15;
    RDouble p0     = 10.1325e2;
    RDouble rho0   = 1.225;
    RDouble t11    = 216.65;
    RDouble p11    = 2.2632e2;
    RDouble rho11  = 3.6392e-1;
    RDouble t20    = t11;
    RDouble p20    = 5.4747e1;
    RDouble rho20  = 8.8035e-2;
    RDouble t32    = 228.65;
    RDouble p32    = 8.6789;
    RDouble rho32  = 1.3225e-2;
    RDouble t47    = 270.65;
    RDouble p47    = 1.1090;
    RDouble rho47  = 1.4275e-3;
    RDouble t52    = t47;
    RDouble p52    = 5.8997e-1;
    RDouble rho52  = 7.5943e-4;
    RDouble t61    = 252.65;
    RDouble p61    = 1.8209e-1;
    RDouble rho61  = 2.5109e-4;
    RDouble t79    = 180.65;
    RDouble p79    = 1.0376e-2;
    RDouble rho79  = 2.0010e-5;
    RDouble t90    = t79;
    RDouble p90    = 1.6437e-3;
    RDouble rho90  = 3.4165e-6;
    RDouble t100   = 210.02;
    RDouble p100   = 3.0070e-4;
    RDouble rho100 = 5.6044e-7;
    RDouble t110   = 257.00;
    RDouble p110   = 7.3527e-5;
    RDouble rho110 = 9.7081e-8;
    RDouble t120   = 349.49;
    RDouble p120   = 2.5209e-5;
    RDouble rho120 = 2.2222e-8;
    RDouble t150   = 892.79;
    RDouble p150   = 5.0599e-6;
    RDouble rho150 = 2.0752e-9;
    RDouble t160   = 1022.20;
    RDouble p160   = 3.6929e-6;
    RDouble rho160 = 1.2336e-9;
    RDouble t170   = 1103.40;
    RDouble p170   = 2.7915e-6;
    RDouble rho170 = 7.8155e-10;
    RDouble t190   = 1205.40;
    RDouble p190   = 1.6845e-6;
    RDouble rho190 = 3.5807e-10;
    RDouble t230   = 1322.30;
    RDouble p230   = 6.7138e-7;
    RDouble rho230 = 1.5640e-10;
    RDouble t300   = 1432.10;
    RDouble p300   = 1.8828e-7;
    RDouble rho300 = 1.9159e-11;
    RDouble t400   = 1487.40;
    RDouble p400   = 4.0278e-8;
    RDouble rho400 = 2.8028e-12;
    RDouble t500   = 1499.20;
    RDouble p500   = 1.0949e-8;
    RDouble rho500 = 5.2148e-13;
    RDouble t600   = 1506.10;
    RDouble p600   = 3.4475e-9;
    RDouble rho600 = 1.1367e-13;
    RDouble t700   = 1507.60;

    RDouble rho = 0.0;

    if (h <= 11019.0) 
    {
        RDouble al1 = (t11 - t0)/11019.0;
        t   = t0 + al1 * h;
        p   = p0 * pow(t/t0, -g/(r * al1));
        rho = rho0 * pow(t/t0, -1.0 - g/(r*al1));
    }
    else if (h <= 20063.0)
    {
        t   = t11;
        p   = p11   * exp(-g*(h-11019.0)/(r*t11));
        rho = rho11 * exp(-g*(h-11019.0)/(r*t11));
    }
    else if (h <= 32162.0)
    {
        RDouble al2 = (t32-t20)/(32162.0-20063.0);
        t   = t11+al2*(h-20063.0);
        p   = p20 * pow(t/t11, -g/(r*al2));
        rho = rho20 * pow(t/t11, -1.0-g/(r*al2));
    }
    else if (h <= 47350.0)
    {
        RDouble al3 = (t47-t32)/(47350.0-32162.0);
        t   = t32+al3*(h-32162.0);
        p   = p32 * pow(t/t32, -g/(r*al3));
        rho = rho32 * pow(t/t32, -1.0-g/(r*al3));
    }

    else if (h <= 52429.0)
    {
        t   = t47;
        p   = p47*exp(-g*(h-47350.0)/(r*t47));
        rho = rho47*exp(-g*(h-47350.0)/(r*t47));
    }
    else if(h <= 61591.0)
    {
        RDouble al4 =(t61-t52)/(61591.0-52429.0);
        t = t47+al4*(h-52429.0);
        p = p52 * pow(t/t47, -g/(r*al4));
        rho = rho52* pow(t/t47,-1.0-g/(r*al4));
    }
    else if(h <= 79994.0)
    {
        RDouble al5=(t79-t61)/(79994.0-61591.0);
        t=t61+al5*(h-61591.0);
        p=p61*pow(t/t61, -g/(r*al5));
        rho=rho61*pow(t/t61,-1.0-g/(r*al5));
    }
    else if(h <= 90000.0)
    {
        t=t79;
        p=p79*exp(-g*(h-79994.0)/(r*t79));
        rho=rho79*exp(-g*(h-79994.0)/(r*t79));
    }
    else if(h <= 100000.0)
    {
        RDouble  al6=(t100-t90)/10000.0;
        t=t79+al6*(h-90000.0);
        p=p90* pow(t/t79,-g/(r*al6));
        rho=rho90* pow(t/t79,-1.0-g/(r*al6));
    }
    else if(h <= 110000.0)
    {
        RDouble al7=(t110-t100)/10000.0;
        t=t100+al7*(h-100000.0);
        p=p100*pow(t/t100, -g/(r*al7));
        rho=rho100*(t/t100, -1.0-g/(r*al7));
    }
    else if(h <= 120000.0)
    {
        RDouble al8=(t120-t110)/10000.0;
        t=t110+al8*(h-110000.0);
        p=p110*pow(t/t110, -g/(r*al8));
        rho=rho110*pow(t/t110, -1.0-g/(r*al8));
    }
    else if(h <= 150000.0)
    {
        RDouble al9=(t150-t120)/30000.0;
        t=t120+al9*(h-120000.0);
        p=p120*pow(t/t120,-g/(r*al9));
        rho=rho120*pow(t/t120,-1.0-g/(r*al9));
    }
    else if(h <= 160000.0)
    {
        RDouble al10=(t160-t150)/10000.0;
        t=t150+al10*(h-150000.0);
        p=p150*pow(t/t150, -g/(r*al10));
        rho=rho150*pow(t/t150, -1.0-g/(r*al10));
    }
    else if(h <= 170000.0)
    {
        RDouble al11=(t170-t160)/10000.0;
        t=t160+al11*(h-160000.0);
        p=p160*pow(t/t160, -g/(r*al11));
        rho=rho160*pow(t/t160,-1.0-g/(r*al11));
    }
    else if(h <= 190000.0)
    {
        RDouble al12=(t190-t170)/20000.0;
        t=t170+al12*(h-170000.0);
        p=p170*pow(t/t170, -g/(r*al12));
        rho=rho170*pow(t/t170,-1.0-g/(r*al12));
    }
    else if(h <= 230000.0)
    {
        RDouble al13=(t230-t190)/40000.0;
        t=t190+al13*(h-190000.0);
        p=p190*pow(t/t190,-g/(r*al13));
        rho=rho190*pow(t/t190,-1.0-g/(r*al13));
    }
    else if(h <= 300000.0)
    {
        RDouble al14=(t300-t230)/70000.0;
        t=t230+al14*(h-230000.0);
        p=p230*pow(t/t230, -g/(r*al14));
        rho=rho230*pow(t/t230, -1.0-g/(r*al14));
    }
    else if(h <= 400000.0)
    {
        RDouble al15=(t400-t300)/100000.0;
        t=t300+al15*(h-300000.0);
        p=p300*pow(t/t300, -g/(r*al15));
        rho=rho300*pow(t/t300, -1.0-g/(r*al15));
    }
    else if(h <= 500000.0)
    {
        RDouble al16=(t500-t400)/100000.0;
        t=t400+al16*(h-400000.0);
        p=p400*pow(t/t400, -g/(r*al16));
        rho=rho400*pow(t/t400, -1.0-g/(r*al16));
    }
    else if(h <= 600000.0)
    {
        RDouble al17=(t600-t500)/100000.0;
        t=t500+al17*(h-500000.0);
        p=p500*pow(t/t500, -g/(r*al17));
        rho=rho500*pow(t/t500, -1.0-g/(r*al17));
    }
    else if(h <= 700000.0)
    {
        RDouble al18=(t700-t600)/100000.0;
        t=t600+al18*(h-600000.0);
        p=p600*pow(t/t600, -g/(r*al18));
        rho=rho600*pow(t/t600, -1.0-g/(r*al18));
    }
    else
    {
        TK_Exit::ExceptionExit("Error: height must be less than 700km!!!\n", true);
    }

    a   = sqrt(1.4*r*t);
    p   = p * 100.0;
    den = rho;
}

//! Set air parameter information through weather database.
void PerfectGas::SetAirInformationByWRFDataBase()
{
    RDouble angleSlide =0.0;

    RDouble refDimensionalTemperature;
    RDouble refDimensionalPressure;
    RDouble refDimensionalDensity;
    RDouble refDimensionalVelocity;
    RDouble angle;
    RDouble refDimensionalSonicSpeed;

    string filePath = GlobalDataBase::GetStrParaFromDB("weatherDataFilePath");
    RDouble longitude = GlobalDataBase::GetDoubleParaFromDB("longitude");
    RDouble latitude =  GlobalDataBase::GetDoubleParaFromDB("latitude");

    //! Read air parameter information from weather data file.
    ReadAirInfo(filePath, longitude, latitude, refDimensionalTemperature, refDimensionalPressure, refDimensionalDensity, angle, refDimensionalVelocity, refDimensionalSonicSpeed);

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    referenceParameterDimensional.SetPressure(refDimensionalPressure);
    referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
    referenceParameterDimensional.SetVelocity(refDimensionalVelocity);
    referenceParameterDimensional.SetSoundVelocity(refDimensionalSonicSpeed);

    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalTemperature", &refDimensionalTemperature, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalSonicSpeed", &refDimensionalSonicSpeed, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("attackd", &angle, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("angleSlide", &angleSlide, PHDOUBLE, 1);
}

void PerfectGas::ReadAirInfo(string &filepath, RDouble &lon, RDouble &lat, RDouble &t, RDouble &p, RDouble &den, RDouble &g, RDouble &v, RDouble &a)
{
    RDouble r      = 287.053;
    RDouble t0     = 273.15;
    RDouble p0     = 1.01325e5;
    RDouble rho0   = 1.292;

    vector <RDouble> pressure;
    vector <RDouble> temperature;
    vector <RDouble> angle;
    vector <RDouble> velocity;

    string var1 = "P";
    string var2 = "T";
    string var3 = "WIND";

    int nfile = 10;
    int heightInternal = 10;
    int fileid = 0;

    for (int ifile = 0; ifile < nfile; ++ ifile)
    {
        fileid = (ifile + 1) * heightInternal;
        ReadWRFSingleData(filepath, var1, fileid, lon, lat, pressure);
        ReadWRFSingleData(filepath, var2, fileid, lon, lat, temperature);
        ReadWRFDoubleData(filepath, var3, fileid, lon, lat, angle, velocity);
    }

    RDouble pp = ComputeVariablesAverage(pressure);
    RDouble tt = ComputeVariablesAverage(temperature);
    RDouble gg = ComputeVariablesAverage(angle);
    RDouble vv = ComputeVariablesAverage(velocity);

    p = pp * 100.0;
    t = t0 + tt;
    //g = 90 + gg;
    g = 270 - gg;
    v = vv;
    a = sqrt(1.4*r*t);
    den = p * rho0 * t0 / ( p0 * t);
}

//! Reading a kind of variable from one file.
void PerfectGas::ReadWRFSingleData(string &filepath, string &varname,int &fileindex, RDouble &longitude, RDouble &latitude, vector< RDouble > & variable)
{
    fstream file;
    string relativeFilename = "15120108.dat";
    relativeFilename = PHSPACE::AddSymbolToFileName(relativeFilename, '_', varname);
    relativeFilename = PHSPACE::AddSymbolToFileName(relativeFilename, '_', fileindex);
    string filename = filepath + relativeFilename;
    file.open(filename.c_str(), ios_base::in);

    string line, word;
    string separator  = " =\t\r\n#$,;\"'";

    RDouble longitudeInterval, latitudeInterval, longitude1, longitude2, latitude1, latitude2;
    int longitudeVertices, latitudeVertices, timeSeries;
    int nlongitudeVertice,nlatitudeVertice;

    SkipLines(file, 1);

    getline(file, line);

    SkipWords(line,5);

    line = FindNextWord(line, word, separator);
    from_string <RDouble> (longitudeInterval, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string <RDouble> (latitudeInterval, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string <RDouble> (longitude1, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string <RDouble> (longitude2, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string <RDouble> (latitude1, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string <RDouble> (latitude2, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string <int> (longitudeVertices, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string <int> (latitudeVertices, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string <int> (timeSeries, word, std::dec);

    bool LocationInMinMaxBox = CheckIfLocationInMinMaxBox(longitude, latitude, longitude1, longitude2, latitude1, latitude2);

    if(!LocationInMinMaxBox)
    {
        PrintToWindow("Geographic scope read from meteorological data : ","\n");
        PrintToWindow("Longitude range : ", longitude1, "--", longitude2, "\n");
        PrintToWindow("Latitude range : ", latitude1, "--", latitude2, "\n");
        TK_Exit::ExceptionExit("Error: the longitude or latitude you set is not within the specified range, please check and reset!!!\n", true);
    }
    else
    {
        RDouble longitudeRatio = (longitude - longitude1)/longitudeInterval + one;
        nlongitudeVertice = static_cast<int>(ROUND(longitudeRatio));
        RDouble latitudeRatio = (latitude - latitude1)/latitudeInterval + one;
        nlatitudeVertice = static_cast<int>(ROUND(latitudeRatio));

        int time  = 12;
        int nline = nlatitudeVertice + (time - 1) * latitudeVertices - 1;
        SkipLines(file, nline);

        getline(file, line);

        int nword = nlongitudeVertice - 1;
        SkipWords(line, nword);

        line = FindNextWord(line, word, separator);
        RDouble dataIn;
        from_string<RDouble>(dataIn, word, std::dec);
        variable.push_back(dataIn);

        file.close();
        file.clear();
    }
}

//! Reading several kinds of variable from one file.
void PerfectGas::ReadWRFDoubleData(string &filepath, string &varname,int &fileindex, RDouble &longitude, RDouble &latitude, vector< RDouble > & variable1, vector< RDouble > & variable2)
{
    fstream file;
    string relativeFilename = "15120108.dat";
    relativeFilename = PHSPACE::AddSymbolToFileName(relativeFilename, '_', varname);
    relativeFilename = PHSPACE::AddSymbolToFileName(relativeFilename, '_', fileindex);
    string filename = filepath + relativeFilename;
    file.open(filename.c_str(), ios_base::in);

    string line, word;
    string separator  = " =\t\r\n#$,;\"'";

    RDouble longitudeInterval, latitudeInterval, longitude1, longitude2, latitude1, latitude2;
    int longitudeVertices, latitudeVertices, timeSeries;
    int nlongitudeVertice,nlatitudeVertice;

    SkipLines(file, 1);

    getline(file, line);

    SkipWords(line,5);

    line = FindNextWord(line, word, separator);
    from_string <RDouble> (longitudeInterval, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string <RDouble> (latitudeInterval, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string <RDouble> (longitude1, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string <RDouble> (longitude2, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string <RDouble> (latitude1, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string <RDouble> (latitude2, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string <int> (longitudeVertices, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string <int> (latitudeVertices, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string <int> (timeSeries, word, std::dec);

    bool LocationInMinMaxBox = CheckIfLocationInMinMaxBox(longitude, latitude, longitude1, longitude2, latitude1, latitude2);

    if(!LocationInMinMaxBox)
    {
        PrintToWindow("Geographic scope read from meteorological data : ","\n");
        PrintToWindow("Longitude range : ", longitude1, "--", longitude2, "\n");
        PrintToWindow("Latitude range : ", latitude1, "--", latitude2, "\n");
        TK_Exit::ExceptionExit("Error: the longitude or latitude is not within the specified range!!!\n", true);
    }
    else
    {
        RDouble longitudeRatio = (longitude - longitude1)/longitudeInterval + 1.5;
        nlongitudeVertice = static_cast<int>(longitudeRatio);
        RDouble latitudeRatio = (latitude - latitude1)/latitudeInterval + 1.5;
        nlatitudeVertice = static_cast<int>(latitudeRatio);

        int time  = 12;
        int nline = nlatitudeVertice + (time - 1) * latitudeVertices - 1;
        SkipLines(file, nline);

        getline(file, line);

        int nword = 2 * (nlongitudeVertice - 1);
        SkipWords(line, nword);

        RDouble dataIn;
        line = FindNextWord(line, word, separator);
        from_string<RDouble>(dataIn, word, std::dec);
        variable1.push_back(dataIn);

        line = FindNextWord(line, word, separator);
        from_string<RDouble>(dataIn, word, std::dec);
        variable2.push_back(dataIn);

        file.close();
        file.clear();
    }

}

RDouble PerfectGas::ComputeVariablesAverage(vector< RDouble > & variable)
{
    RDouble vv = 0.0;
    for ( int m = 0; m < variable.size(); ++ m )
    {
        vv += variable[m] / variable.size();
    }
    return vv;
}

void PerfectGas::ComputeReferenceParameter()
{
    //! The first two variables is temperature independent, 
    //! depend on molecular weight and mass fraction of species if chemical reactions exist.
    ComputeReferenceMolecularInformation();

    ComputeReferenceGeneralGasConstant();

    //! The following five variables depend on temperature, 
    //! and mass fraction of species if chemical reactions exist.
    ComputeReferenceViscosityDimensional();

    ComputeReferenceSpecificHeatRatio();

    ComputeReferenceSoundVelocity();

    ComputeReferenceVelocity();

    ComputeCoefficientOfStateEquation();
    //! The above seven variables should be computed firstly.

    //! The variables in this function are computed by different methods according to "inflowParaType", 
    //! including density, pressure or reynolds.
    ComputeReferenceGasInformation();

    ComputeReferencePrimitive();

    ComputePTotal();
}

void PerfectGas::ComputeReferencePrimitive()
{
    RDouble refMachNumber = referenceParameterDimensional.GetMachNumber();
    RDouble reference_sound_velocity, reference_pressure;

    RDouble massoo_ref = 1.0;

    RDouble reference_density     = 1.0;
    RDouble reference_temperature = 1.0;

    reference_sound_velocity = 1.0 / refMachNumber;

    reference_pressure = coefficientOfStateEquation * 
    reference_density * reference_temperature / massoo_ref;

    using namespace IDX;
    int nLaminar = GlobalDataBase::GetIntParaFromDB("nl");
    RDouble *prim_inf = new RDouble[nLaminar];

    prim_inf[IR] = reference_density;

    GlobalBoundaryCondition::SetSpeedDrictionForReferenceVar(one, prim_inf);

    prim_inf[IP] = reference_pressure;

    GlobalDataBase::UpdateDataPtr("prim_inf", prim_inf);
}

void PerfectGas::ComputeReferenceGasInformation()
{
    int inflowParaType = GlobalDataBase::GetIntParaFromDB("inflowParaType");

    if (inflowParaType == NONDIMENSIONCONDITION)    //! Obtain the referce density and referce pressure for nondimension condition.
    {
        ComputeDensityWithReynoldsNumber();
        ComputePressureInGasStateEquation();
    }
    else if (inflowParaType == FLIGHTCONDITION || inflowParaType == EXPERIMENTCONDITION)    //! obtain the referce density for nondimension, and the Re.
    {
        NormalizeAirInformation();
        ComputeReferenceReynoldsNumber();
    }
    else if (inflowParaType == TEMPERATURE_DENSITY)
    {
        RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();
        RDouble refDimensionalDensity = referenceParameterDimensional.GetDensity();
        RDouble gasConstant = referenceParameterDimensional.GetAverageGeneralGasConstant();

        RDouble refDimensionalPressure = refDimensionalDensity * gasConstant * refDimensionalTemperature;
        referenceParameterDimensional.SetPressure(refDimensionalPressure);
        GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);

        ComputeReferenceReynoldsNumber();
    }
    else if (inflowParaType == TEMPERATURE_PRESSURE || inflowParaType == MACH_TEMP_PRE || inflowParaType == WEATHERCONDITION || inflowParaType == WINDSPEEDPROFILE)
    {
        NormalizeAirInformation();
        ComputeReferenceReynoldsNumber();
    }
    else
    {
        TK_Exit::UnexpectedVarValue("inflowParaType", inflowParaType);
    }
}

void PerfectGas::ComputeReferenceReynoldsNumber()
{
    RDouble refDimensionalDensity = referenceParameterDimensional.GetDensity();
    RDouble reference_velocity_dimensional  = referenceParameterDimensional.GetVelocity();
    RDouble reference_viscosity_dimensional = referenceParameterDimensional.GetViscosity();

    RDouble refReNumber = refDimensionalDensity * reference_velocity_dimensional / reference_viscosity_dimensional;
    RDouble oreynolds = 1.0 / refReNumber;

    referenceParameterDimensional.SetReynoldsNumber(refReNumber);
    GlobalDataBase::UpdateData("refReNumber", &refReNumber, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("oreynolds", &oreynolds, PHDOUBLE, 1);
}

void PerfectGas::NormalizeAirInformation()
{
    RDouble refDimensionalPressure = referenceParameterDimensional.GetPressure();
    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();
    RDouble reference_average_generanl_gas_constant = referenceParameterDimensional.GetAverageGeneralGasConstant();

    RDouble refDimensionalDensity = refDimensionalPressure / reference_average_generanl_gas_constant / refDimensionalTemperature;

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
}

void PerfectGas::ComputePressureInGasStateEquation()
{
    RDouble refDimensionalDensity = referenceParameterDimensional.GetDensity();
    RDouble reference_average_generanl_gas_constant = referenceParameterDimensional.GetAverageGeneralGasConstant();
    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();

    RDouble refDimensionalPressure = refDimensionalDensity * reference_average_generanl_gas_constant * refDimensionalTemperature;

    referenceParameterDimensional.SetPressure(refDimensionalPressure);
    GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);
}

void PerfectGas::ComputeDensityWithReynoldsNumber()
{
    RDouble refReNumber = referenceParameterDimensional.GetReynoldsNumber();
    RDouble reference_viscosity_dimensional = referenceParameterDimensional.GetViscosity();
    RDouble reference_velocity_dimensional = referenceParameterDimensional.GetVelocity();

    RDouble refDimensionalDensity = refReNumber * reference_viscosity_dimensional / reference_velocity_dimensional;

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
    
    //!for AleManager classs variables
    GlobalDataBase::UpdateData("referenceDensity", &refDimensionalDensity, PHDOUBLE, 1);
}

void PerfectGas::ComputeCoefficientOfStateEquation()
{
    //! No matter whether chemical reactions exist or not.
    RDouble gama0 = referenceParameterDimensional.GetAverageSpecificHeatRatio();
    RDouble refMachNumber = referenceParameterDimensional.GetMachNumber();

    //! The variable is the non-dimensional value of 8.314/M??M denotes the non-dimensioanl molecular weight of the mixture gas.
    coefficientOfStateEquation = 1.0 / (gama0 * refMachNumber * refMachNumber);    //! universal gas constant.
    GlobalDataBase::UpdateData("coefficientOfStateEquation", &coefficientOfStateEquation, PHDOUBLE, 1);
}

void PerfectGas::ComputeReferenceVelocity()
{
    RDouble refDimensionalSonicSpeed = referenceParameterDimensional.GetSoundVelocity();

    int inflowParaType = GlobalDataBase::GetIntParaFromDB("inflowParaType");

    if (inflowParaType == TEMPERATURE_DENSITY || inflowParaType == TEMPERATURE_PRESSURE || inflowParaType == WINDSPEEDPROFILE)
    {
        RDouble refDimensionalVelocity = referenceParameterDimensional.GetVelocity();
        RDouble refMachNumber = refDimensionalVelocity / refDimensionalSonicSpeed;
        referenceParameterDimensional.SetMachNumber(refMachNumber);
        GlobalDataBase::UpdateData("refMachNumber", &refMachNumber, PHDOUBLE, 1);
    }
    else if (inflowParaType == WEATHERCONDITION)
    {
        RDouble refDimensionalVelocity = referenceParameterDimensional.GetVelocity();
        GlobalDataBase::UpdateData("refDimensionalVelocity", &refDimensionalVelocity, PHDOUBLE, 1);
        RDouble refMachNumber = refDimensionalVelocity / refDimensionalSonicSpeed;
        referenceParameterDimensional.SetMachNumber(refMachNumber);
        GlobalDataBase::UpdateData("refMachNumber", &refMachNumber, PHDOUBLE, 1);
    }
    else
    {
        RDouble refMachNumber = referenceParameterDimensional.GetMachNumber();
        RDouble refDimensionalVelocity = refDimensionalSonicSpeed * refMachNumber;
        referenceParameterDimensional.SetVelocity(refDimensionalVelocity);
        GlobalDataBase::UpdateData("refDimensionalVelocity", &refDimensionalVelocity, PHDOUBLE, 1);
        //!for AleManager classs variables
        GlobalDataBase::UpdateData("referenceVelocity", &refDimensionalVelocity, PHDOUBLE, 1);
    }
}

void PerfectGas::ComputeReferenceSoundVelocity()
{
    RDouble gama0 = referenceParameterDimensional.GetAverageSpecificHeatRatio();
    RDouble reference_general_gas_constant = referenceParameterDimensional.GetAverageGeneralGasConstant();
    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();

    //! unit: m/s.
    RDouble refDimensionalSonicSpeed = sqrt(gama0 * reference_general_gas_constant * refDimensionalTemperature);

    referenceParameterDimensional.SetSoundVelocity(refDimensionalSonicSpeed);
    GlobalDataBase::UpdateData("refDimensionalSonicSpeed", &refDimensionalSonicSpeed, PHDOUBLE, 1);
}

void PerfectGas::ComputeReferenceSpecificHeatRatio()
{
    RDouble refGama = 1.4;
    refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    referenceParameterDimensional.SetAverageSpecificHeatRatio(refGama);
}

void PerfectGas::ComputeReferenceViscosityDimensional()
{
    if (this->nGasModel > 1)    //! The one molecule of the gas.
    {
        //! Referred to the "2008-Lofthouse-Nonequilibrium Hypersonic Aerothermodynamics Using the Direct Simulation Monte Carlo And Navier-Stokes Models".
        RDouble T = referenceParameterDimensional.GetTemperature();
        //! The default is the Argon.
        RDouble Omega = 0.734, Tref = 1000.0, Dref = 3.595e-10;
        if (this->nGasModel == 3)    //! The pure gas of Nitrogen.
        {
            Omega = 0.70;
            Tref = 290.0;
            Dref = 4.11e-10;
        }
        RDouble mass = referenceParameterDimensional.GetAverageMolecularWeight() / AVOGADRO_CONSTANT;
        RDouble refViscosity = 15 * sqrt(PI * mass * BOLTZMANN_CONSTANT * Tref);
        refViscosity = refViscosity / (2 * PI * Dref * Dref) / (5 - 2 * Omega) / (7 - 2 * Omega);

        RDouble reference_viscosity_dimensional = refViscosity * pow(T / Tref, Omega);
        referenceParameterDimensional.SetViscosity(reference_viscosity_dimensional);

        RDouble tsi  = 110.4;
        RDouble tsuth = tsi / T;
        GlobalDataBase::UpdateData("tsuth", &tsuth, PHDOUBLE, 1);
    }
    else
    {
        RDouble t0i  = 288.15;
        RDouble tsi  = 110.4;
        RDouble mu0i = 1.7894E-05;

        RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();
        RDouble reference_viscosity_dimensional = (t0i + tsi) / (refDimensionalTemperature + tsi) * pow(refDimensionalTemperature / t0i, 1.5) * mu0i;
        referenceParameterDimensional.SetViscosity(reference_viscosity_dimensional);

        RDouble tsuth = tsi / refDimensionalTemperature;
        GlobalDataBase::UpdateData("tsuth", &tsuth, PHDOUBLE, 1);
    }
}

void PerfectGas::ComputeReferenceMolecularInformation()
{
    //! unit: kg/mol.
    RDouble refMolecularWeight = 28.9644;

    molecularDiameter = 3.6e-10;

    if (nGasModel == 1)    //! The Carbon dioxide.
    {
        refMolecularWeight = 44.01;
        molecularDiameter = 3.370e-10;    //! unit: m
    }
    else if (nGasModel == 2)    //! The Argon.
    {
        refMolecularWeight = 40.0;
        molecularDiameter = 3.595e-10;    //! unit: m
    }
    else if (nGasModel == 3)    //! The Nitrogen.
    {
        refMolecularWeight = 28.01293;
        molecularDiameter = 4.110e-10;    //! unit: m
    }

    GlobalDataBase::UpdateData("refMolecularWeight", &refMolecularWeight, PHDOUBLE, 1);
    RDouble reference_average_molecular_weight_dimensional = refMolecularWeight * 1.0e-3;
    referenceParameterDimensional.SetAverageMolecularWeight(reference_average_molecular_weight_dimensional);
}

void PerfectGas::ComputeReferenceGeneralGasConstant()
{
    RDouble reference_average_molecular_weight_dimensional = referenceParameterDimensional.GetAverageMolecularWeight();
    RDouble reference_average_general_gas_constant         = rjmk / reference_average_molecular_weight_dimensional;

    referenceParameterDimensional.SetAverageGeneralGasConstant(reference_average_general_gas_constant);
    GlobalDataBase::UpdateData("reference_average_general_gas_constant", &reference_average_general_gas_constant, PHDOUBLE, 1);
}

void PerfectGas::ComputeOtherProperParameterForGrid()
{
    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

    RDouble forceReferenceLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
    RDouble forceReferenceLengthSpanWise = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");
    RDouble forceReferenceArea = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
    RDouble TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    RDouble TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    RDouble TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

    if (GetDim() == THREE_D)
    {
        forceReferenceArea /= reynoldsReferenceLengthDimensional*reynoldsReferenceLengthDimensional;
    }
    else if (GetDim() == TWO_D)
    {
        forceReferenceArea /= reynoldsReferenceLengthDimensional;
    }
    forceReferenceLength /= reynoldsReferenceLengthDimensional;
    forceReferenceLengthSpanWise /= reynoldsReferenceLengthDimensional;
    TorqueRefX /= reynoldsReferenceLengthDimensional;
    TorqueRefY /= reynoldsReferenceLengthDimensional;
    TorqueRefZ /= reynoldsReferenceLengthDimensional;

    GlobalDataBase::UpdateData("forceReferenceLength", &forceReferenceLength, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("forceReferenceLengthSpanWise", &forceReferenceLengthSpanWise, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("forceReferenceArea", &forceReferenceArea, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("TorqueRefZ", &TorqueRefZ, PHDOUBLE, 1);

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady == 1)
    {
        PrintToWindow("Unsteady Informations" , "\n");
        RDouble reference_velocity_dimensional = referenceParameterDimensional.GetVelocity();
        RDouble unit = reynoldsReferenceLengthDimensional / reference_velocity_dimensional;
        RDouble physicalTimeStepDimensional = 0.0;
        RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");
        if (GlobalDataBase::IsExist("physicalTimeStepDimensional", PHDOUBLE, 1))
        {
            physicalTimeStepDimensional = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStepDimensional");
            if (physicalTimeStepDimensional >= 0.0)
            {
                physicalTimeStep = physicalTimeStepDimensional / unit;
            }
        }
        physicalTimeStepDimensional = physicalTimeStep * unit;
        GlobalDataBase::UpdateData("physicalTimeStepDimensional", &physicalTimeStepDimensional, PHDOUBLE, 1);
        //! The non-dimensional physical time step based on grid scale.
        physicalTimeStep /= reynoldsReferenceLengthDimensional;
        GlobalDataBase::UpdateData("physicalTimeStep", &physicalTimeStep, PHDOUBLE, 1);

        PrintToWindow("  Physical Time Step       :    ", physicalTimeStepDimensional , "s", "\n");
    }
}

RDouble PerfectGas::GetUniversalGasConstant()
{
    //! Return the non-dimensional universal gas constant.
    return coefficientOfStateEquation;
}

RDouble PerfectGas::ComputeReferenceTotalEnthalpy()
{
    RDouble totalEnthalpy = 0.0;
    RDouble gamma = referenceParameterDimensional.GetAverageSpecificHeatRatio();
    RDouble density = referenceParameterDimensional.GetDensity();
    RDouble pressure = referenceParameterDimensional.GetPressure();
    RDouble velocity = referenceParameterDimensional.GetVelocity();
    totalEnthalpy = gamma / (gamma - 1.0) * (pressure / density) + 0.5 * velocity * velocity;
    return totalEnthalpy;
}

void PerfectGas::ComputePTotal()
{
    RDouble gama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble refMachNumber = referenceParameterDimensional.GetMachNumber();
    RDouble referenceSoundVelocity, referencePressure;

    RDouble massooRef = 1.0;

    RDouble referenceDensity = 1.0;
    RDouble referenceTemperature = 1.0;

    RDouble pTotal;

    referenceSoundVelocity = 1.0 / refMachNumber;

    referencePressure = coefficientOfStateEquation * referenceDensity * referenceTemperature / massooRef;
    RDouble m12 = refMachNumber * refMachNumber;
    RDouble gamp1 = gama + one;
    RDouble gamm1 = gama - one;
    RDouble c1 = two / gamm1 + m12;
    RDouble c2 = two * gama / gamm1 * m12 - one;

    RDouble cc = gama / (gama - one);
    RDouble t0t = one + half * gamm1 * c1 / c2;

    if (refMachNumber > one)
    {
        referencePressure *= (two * gama / gamp1 * m12 - gamm1 / gamp1);
    }

    pTotal = referencePressure * pow(t0t, cc);

    GlobalDataBase::UpdateData("pTotal", &pTotal, PHDOUBLE, 1);
}

}

//! Check if the longitude and latitude are within the specified range in the weather data file
bool CheckIfLocationInMinMaxBox(RDouble longitude,RDouble latitude,RDouble a,RDouble b,RDouble c,RDouble d)
{
    if ( longitude < MIN(a, b) || longitude > MAX(a, b) ) return false;
    if ( latitude < MIN(c, d) || latitude > MAX(c, d) ) return false;
    return true;
}
}
