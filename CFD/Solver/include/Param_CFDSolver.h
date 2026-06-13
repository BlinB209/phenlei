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
//! @file      Param_CFDSolver.h
//! @brief     Record parameters of CFD Solver.
//! @author    Bell, Zhang Jian, Wan Yunbo, Meng Liyuan.

#pragma once
#include "LIB_Macro.h"
#include "Precision.h"

namespace PHSPACE
{

class Param_CFDSolver
{
public:
    LIB_EXPORT Param_CFDSolver();

    LIB_EXPORT virtual ~Param_CFDSolver();

public:
    //! Get restart flow file name.
    virtual const string & GetRestartNSFile() const { return restartNSFile; }

    //! Get res file name.
    virtual const string & GetResSaveFile() const { return resSaveFile; }

    //! Get flow field name.
    virtual const string & GetFlowFieldTecFile() const { return flowFieldTecFile; };

    //! Init all parameters.
    LIB_EXPORT virtual void Init();

    //! Get reference mach number.
    RDouble GetRefMachNumber() const;

    //! Get reference reynolds number.
    RDouble GetRefReNumber() const;

    //! Get the count backwards of reference reynolds number.
    RDouble GetoRefReNumber() const;

    //! Get if need plot volume field.
    int GetPlotFieldType() const;

    //! Get lower range of plotting field.
    const RDouble * GetLowerPlotFieldBox() const;

    //! Get upper range of plotting field.
    const RDouble * GetUpperPlotFieldBox() const;

    //! Get viscous model type.
    int GetViscousType() const;

    //! Get viscous model name
    const string & GetViscousName() const;

    //! Return the flow is viscous or not.
    bool IsViscous() const;

    //! Get if unsteady of steady.
    int GetIsUnsteady() const;

    //! Get interval steps of dumping residuals.
    int GetIntervalStepRes() const;

    //! Get interval steps of dumping flow field.
    int GetIntervalStepPlot() const;

    //! Get if using ALE method.
    int GetIsCodeOfAleModel() const;

    //! Get CFL when computing start.
    RDouble GetCFLStart() const;

    //! Get final CFL.
    RDouble GetCFLEnd() const;

    RDouble GetFineCFL() const;

    RDouble GetPMax() const;
    RDouble GetPMin() const;
    RDouble GetDeltaMax() const;
    RDouble GetMagnifyFactor() const;
    RDouble GetReduceFactor() const;

    //! The number of step when cfl increase from CFLStart to CFLEnd.
    int GetCFLVaryStep() const;

    //! Get number of MultiGrid level.
    int GetNMGLevel() const;

    //! Get LUSGS sweeps.
    int GetNLUSGSSweeps() const;

    //! Get tolerance of LUSGS scheme.
    RDouble GetLUSGSTolerance() const;
    int GetIntervalStepFlow() const;

    int GetIfLocalTimeStep() const;

    RDouble GetTorqueRefX() const;
    RDouble GetTorqueRefY() const;
    RDouble GetTorqueRefZ() const;

    //! Get reference gama value.
    RDouble GetRefGama() const;

    //! Get reference dimensional density.
    RDouble GetRefDimensionalPressure() const;

    //! Get reference dimensional density.
    RDouble GetRefDimensionalDensity() const;

    //! Get reference dimensional velocity.
    RDouble GetRefDimensionalVelocity() const;

    //! Get reference dimensional sonic speed.
    RDouble GetRefDimensionalSonicSpeed() const;

    //! Get reference dimensional temperature.
    RDouble GetRefDimensionalTemperature() const;

     //! Get AusmpwPlusLimiter.
    RDouble GetAusmpwPlusLimiter() const;

    //! Get trTemperatureMin
    RDouble GetTrTemperatureMinNonDim() const;
    //! Get flag of starting from steady result.
    int GetIfStartFromSteadyResults() const;

    //! Get flag of statics flow field.
    int GetIfStaticsFlowField() const;

    //! Get flag of statics reynolds stress.
    int GetIfStaticsReynoldsStress() const;

    //! Get flag of starting with changing inflow.
    int GetIsRestartChangeInflow() const;

    //! Get method of computing statistical reynolds stress.
    int GetStatisticMethod() const;

    int GetIsOverLapping() const;

    RDouble GetSkewnessAngle() const;
    RDouble GetTurbSkewnessAngle() const;

    //! Get parameter judge whether modify Wenn Scheme.
    int GetWennSchemeFlag() const;

    int GetnDiagonalModifiedTurb() const;

    //! Get if carry out precondition or not.
    int GetIfLowSpeedPrecon() const;

    int GetnNumberOfSpeedStep() const;
    int *GetSpeedVaryStep() const;
    RDouble *GetSpeedVaryCoef() const;

    int GetnNumberOfCFLStep() const;
    int *GetCFLVaryMultiStep() const;
    RDouble *GetCFLVaryCoef() const;

private:
    string restartNSFile;
    string resSaveFile;
    string flowFieldTecFile;

    RDouble refMachNumber;
    RDouble refReNumber;
    RDouble oRefReNumber;

    int plotFieldType;

    RDouble lowerPlotFieldBox[3];
    RDouble upperPlotFieldBox[3];

    int viscousType;
    string viscousName;

    int isUnsteady;
    int intervalStepRes;
    int intervalStepPlot;

    int isAle;

    RDouble CFLStart;
    RDouble CFLEnd;
    RDouble fineCFL;
    int CFLVaryStep;
    RDouble pMax, pMin, deltaMax, magnifyFactor, reduceFactor;

    int nNumberOfSpeedStep;
    int *speedVaryStep;
    RDouble *speedVaryCoef;

    int nNumberOfCFLStep;
    int *CFLVaryMultiStep;
    RDouble *CFLVaryCoef;

    int nMGLevel;

    int nLUSGSSweeps;
    RDouble LUSGSTolerance; 

    RDouble refGama;
    RDouble refDimensionalPressure;
    RDouble refDimensionalDensity;      
    RDouble refDimensionalSonicSpeed;  
    RDouble refDimensionalTemperature;
    RDouble refDimensionalVelocity;

    RDouble AusmpwPlusLimiter;

    RDouble trTemperatureMinNonDim;

    int ifStartFromSteadyResults;
    int ifStaticsFlowField;
    int ifStaticsReynoldsStress;
    int isRestartChangeInflow;
    int statisticMethod;
    int intervalStepFlow;
    int ifLocalTimeStep;

    RDouble TorqueRefX;
    RDouble TorqueRefY;
    RDouble TorqueRefZ;

    int isOverset;
    int isWennScheme;

    RDouble skewnessAngle;
    RDouble turbSkewnessAngle;

    int nDiagonalModifiedTurb;

    //! For incompressible flow (mach < 0.3), using precondition method or not.
    int ifLowSpeedPrecon;
};

#include "Param_CFDSolver.hxx"
}