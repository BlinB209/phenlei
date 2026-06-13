//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      Transition.h
//! @brief     Record paramters of TransitionSolver.
//! @author    He Kun.

#pragma once
#include "LIB_Macro.h"
#include "Precision.h"
#include "Math_BasisFunction.h"

using namespace std;

namespace PHSPACE
{
    RDouble ReynoldsNumberBasedOnStrainRate(RDouble density, RDouble distanceToWall, RDouble viscosity,
        RDouble strainrate, RDouble reynoldsNumberInflow);

    RDouble ReynoldsNumberBasedOnDissipation(RDouble density, RDouble distanceToWall, RDouble viscosity,
        RDouble dissipationRate, RDouble reynoldsNumberInflow);

    RDouble ReynoldsDissipation(RDouble density, RDouble distanceToWall, RDouble viscosity,
        RDouble dissipationRate, RDouble reynoldsNumberInflow);

    RDouble ReynoldsNumberBasedOnScf(RDouble heightCrossFlow, RDouble velocity, RDouble density, 
        RDouble viscosityLaminar, RDouble viscosityTurbulence, RDouble reynoldsNumberInflow);

    RDouble TimeScaleInSourceTerm(RDouble density, RDouble velocity, RDouble viscosityLaminar, RDouble viscosityTurbulence, RDouble reynoldsNumberInflow);

    RDouble MomentumThickness(RDouble density, RDouble velocity, RDouble viscosity, RDouble onsetReynoldsOnMomentumThickness,
        RDouble reynoldsNumberInflow);

    RDouble FlengthGivenByLangtry(RDouble localTransitiononsetReynoldsOnMomentumThickness);
    RDouble FlengthCaliBrated(RDouble localTransitiononsetReynoldsOnMomentumThickness);
    RDouble FlengthUnstrCaliBrated(RDouble localTransitiononsetReynoldsOnMomentumThickness);

    RDouble HighReynoldsCorrectionOfFlength(RDouble ReynoldsNumberBasedOnDissipation, RDouble Flengthori);

    RDouble ViscosityRatio(RDouble density, RDouble viscosity, RDouble kineticEnergy, RDouble dissipationRate, RDouble reynoldsInflow);

    RDouble TransitionOnsetFunction(RDouble Rev, RDouble Rectac, RDouble RT);
    RDouble TransitionOnsetMomentumThicknessReynolds(RDouble Rectabar);
    RDouble TransitionOnsetMomentumThicknessReynoldsCaliBrated(RDouble Rectabar);
    RDouble TransitionOnsetMomentumThicknessReynoldsUnstrCaliBrated(RDouble Rectabar);    
    RDouble SeparationCorrectionOfIntermittency(RDouble gmori, RDouble Rev, RDouble Rectac, RDouble RT, RDouble Fctat, RDouble s1, RDouble refMaNumber);
    RDouble TurbulenceIntensity(RDouble velocity, RDouble kineticEnergy);
    RDouble ControlFunctionFturb(RDouble RT);
    void CorrectionOfBlendingFunctionInSST(RDouble density, RDouble distance, RDouble viscosity, RDouble kineticEnergy, RDouble reynoldsInflow, RDouble &F1);

    RDouble CorrectionOfDestructionInKEquation(RDouble gmeff, RDouble Dk);
    RDouble CorrectionOfProductionInKEquation(RDouble gmeff, RDouble Pk);
    RDouble EmpiricalCorrelationOfFlamdacta(RDouble Tu, RDouble lamdacta);
    RDouble EmpiricalCorrelationOfRectat(RDouble Tu, RDouble Flamdacta);
    RDouble AccelerationAlongStreamline(RDouble u, RDouble v, RDouble w, RDouble dudx, RDouble dudy, RDouble dudz,
        RDouble dvdx, RDouble dvdy, RDouble dvdz, RDouble dwdx, RDouble dwdy, RDouble dwdz);
    void BlendingFunctionOfFctat(RDouble gm, RDouble Rew, RDouble Rectabar, RDouble vorticity, RDouble viscosity, RDouble density,
        RDouble absu, RDouble distance, RDouble reynoldsInflow, RDouble ce2, RDouble &Fctat, RDouble &Fctatcf);

    RDouble PressureGradientFunction(RDouble density, RDouble momentumThickness, RDouble viscosity, RDouble dUds, RDouble reynoldsInflow);
}