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
//! @file      PerfectGas.h
//! @brief     Explain this file briefly.
//! @author    Li Peng, He Xin.

#pragma once
#include "Precision.h"
#include "Gas.h"
#include <string>
using namespace std;

namespace PHSPACE
{
class DataContainer;
namespace GAS_SPACE
{

class PerfectGas :public Gas
{
private:
    FluidParameter referenceParameterDimensional;
private:
    RDouble coefficientOfStateEquation;    //! universal gas constant, R
    int nGasModel;    //! The gas model, 0 is the earth atmosphere, 1 the Mars atmosphere, 2 the Argon, 3 the Nitrogen.
    RDouble molecularDiameter;
public:
    PerfectGas();
    ~PerfectGas();

public:
    void GetReferenceParameters(FluidParameter &refParam){refParam = referenceParameterDimensional;};
    RDouble GetMolecularDiameter() {return this->molecularDiameter;};
    int GetGasModelType() {return this->nGasModel;};
public:
    //! The primitive variables are transfered to the conservative variables.
    //! @param[in ]: q prim indicates the primitive variables.
    //! @param[in ]: gama denotes the specific heat ratio, which is valid in the perfect gas.
    //! @param[in ]: Tv denotes the vibrational temperature(The non-dimensional value).
    //! @param[in ]: Te denotes the electron temperature(The non-dimensional value).
    //! @param[out]: indicates the conservative variables.
    void Primitive2Conservative(RDouble *prim, RDouble gama, RDouble Tv, RDouble Te, RDouble *q);

#ifdef USE_GMRESSOLVER
    //! GMRESPV 
    //! the derivative of primitive variables w.r.t conservative variables
    void dPrimitive2dConservative(RDouble *prim, RDouble gama, RDouble** dqdcv);

    //! GMRESVis
    //! the derivative of primitive variables w.r.t conservative variables
    void dGradient2dPrimitive(RDouble *prim, int sign, RDouble** dgraddq, char direction, RDouble nxs, RDouble nys, RDouble nzs, RDouble ns, RDouble vol,int nPara = 6);
    
    //! GMRESnolim
    void dGradient2dPrimitive4limiter(int sign, RDouble** dgraddq, char direction, RDouble nxs, RDouble nys, RDouble nzs, RDouble ns, RDouble vol);
#endif

    //! To obtain the total energy E = e + 1/2*V2.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables, gama is the specific heat ratio of the mixture gas.
    //! @param[in ]: vibrationTemperature is non-dimensional value of vibration temperature.
    //! @param[in ]: electronTemperature is non-dimensional value of electron  temperature.
    //! @param[OUT]: totalEnergy denotes the total internal energy of the mixed gas.
    void ComputeInternalEnergy(RDouble *primitiveVariables, RDouble gama, RDouble vibrationTemperature, RDouble electronTemperature, RDouble &totalEnergy);

    //! The conservative variables are transfered to the primitive variables.
    //! @param[in ]: q indicates the conservative variables.
    //! @param[in ]: gama denotes the specific heat ratio, which is valid in the perfect gas.
    //! @param[out]: prim indicates the primitive variables.
    //! @param[in/out]: temperature denotes the array of temperatures in the current time-advancing step and next time-advancing step.
    //! As input parameter, temperature indicates the temperature of the current time-advancing step. Meanwhile,
    //! it denotes the temperature of the next time-advancing step as output parameter.
    void Conservative2Primitive(RDouble *q, RDouble gama, RDouble *prim, RDouble *temperature);

    virtual void GetSpecificHeatRatio(RDouble *prim, RDouble &gama);

    //! Compute the static enthalpy with the primitive variables.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: gama is the specific heat ratio.
    //! @param[out]: enthalpy denotes the static enthalpy.
    //! @param[in ]: temperatures stores the translation-rotation temperature, vibration temperature and the electron temperature.
    void ComputeEnthalpyByPrimitive(RDouble *primitiveVariables, RDouble &gama, RDouble &enthalpy, RDouble* temperatures);

    //! This function is used for perfect gas and single temperature model.
    //! To compute the total enthalpy and the variable dh in the function MXDQ_STD(). dh=b2, b2 denotes the coefficient \n
    //! of the vector M*dQ which can be referred to the forumla (A.7) and (A.8) in the appendix A of the PHengLEI Theory manual.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: gama is the specific heat ratio.
    //! @param[in ]: deltaQ is an array to saving the differences of the primitive variables.
    //! @param[out]: totalEnthalpy denotes the total enthalpy.
    //! @param[out]: deltaEnthalpy denotes the difference of static enthalpy.
    void ComputeTotalEnthalpyAndDH(RDouble *primitiveVariables, RDouble &gama, RDouble *deltaQ, RDouble &totalEnthalpy, RDouble &deltaEnthalpy);


    //! To compute the total enthalpy and the difference of the static enthalpy.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: gama is the specific heat ratio.
    //! @param[in ]: deltaQ is an array to saving the differences of the primitive variables.
    //! @param[out]: deltaEnthalpy denotes the difference of static enthalpy.
    //! @param[out]: totalEnthalpy denotes the total enthalpy.
    void ComputeDHAndTotalEnthalpy(RDouble *primitiveVariables, RDouble &gama, RDouble *deltaQ, RDouble &deltaEnthalpy, RDouble &totalEnthalpy);

    //! To compute the static enthalpy of the mixed gas.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: transRotationTemperature denotes the temperature of translation-rotation mode.
    //! @param[in ]: vibrationTemperature denotes the temperature of vibration mode.
    //! @param[in ]: electronTemperature denotes the temperature of electron mode.
    RDouble GetMixedGasEnthalpy(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature);

    //! To obtain the non-dimensional viscosity of the mixture gas according to the primary variables.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: transRotationTemperature is the non-dimensional translation-rotation temperature of the mixture gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[in ]: suthTemperature denotes the temperature value in the Sutherland fomula.
    //! @param[out]: viscosity is the non-dimensional viscosity of the mixture gas [N*s/m2].
    void ComputeViscosityByPrimitive(RDouble *primitiveVariables, RDouble &transRotationTemperature, RDouble &electronTemperature, RDouble &suthTemperature, RDouble &viscosity);

    void InitGlobalParameterOfNSEquation();

    void InitCommonParameter();

    void InitReferenceParameter();

    void InitNumberOfSpecies();

    void InitOtherParameter();

    //! Obtain the T, p, rho and c by the height (Rg=287.053, gamma=1.4 are used). However, only the T and p is useful. 
    //! The other two variables will be recomputed to be compatible with the condition of chemical reactions existing.
    void SetAirInformationByDataBase();

    //! Obtain the T, p, rho and c by the T0, P0, Ma (Rg=287.053, gamma=1.4 are used). However, only the T and p is useful. 
    //! The other two variables will be recomputed to be compatible with the condition of chemical reactions existing.
    void SetAirInformationByExpData(int useSetting);
    void SetAirInformationByExpData();

    void GetAirInfo(const RDouble &height, RDouble &temperature, RDouble &pressure, RDouble &density, RDouble &soundspeed);

    void SetAirInformationByWRFDataBase();
    void ReadAirInfo(string &filepath, RDouble &longitude, RDouble &latitude, RDouble &temperature, RDouble &pressure, RDouble &density, RDouble &angle, RDouble &velocity, RDouble &soundspeed);
    void ReadWRFSingleData(string &filepath, string &varname,int &fileindex, RDouble &longitude, RDouble &latitude, vector< RDouble > & variable);
    void ReadWRFDoubleData(string &filepath, string &varname,int &fileindex, RDouble &longitude, RDouble &latitude, vector< RDouble > & variable1, vector< RDouble > & variable2);
    RDouble ComputeVariablesAverage(vector< RDouble > & variable);


    void ComputeReferenceParameter();

    void ComputeDensityWithReynoldsNumber();

    void ComputeReferencePrimitive();

    void ComputeReferenceGasInformation();

    void ComputeReferenceReynoldsNumber();

    //! This function is mainly for the average density of the mixture in chemical reactions, 
    //! the value of mixture density has little difference with that of no chemical reaction, which can be ignored.
    void NormalizeAirInformation();

    void ComputePressureInGasStateEquation();

    void ComputeCoefficientOfStateEquation();

    void ComputeReferenceVelocity();

    void ComputeReferenceSoundVelocity();

    void ComputeReferenceSpecificHeatRatio();

    void ComputeReferenceViscosityDimensional();

    void ComputeProperReynoldsNumberForGrid();

    void ComputeOtherProperParameterForGrid();

    void ComputeReferenceGeneralGasConstant();

    void ComputeReferenceMolecularInformation();

    RDouble GetCoefficientOfStateEquation() { return coefficientOfStateEquation; };

    RDouble GetUniversalGasConstant();

    RDouble ComputeReferenceTotalEnthalpy();

    void ComputePTotal();

};

}
bool CheckIfLocationInMinMaxBox(RDouble longitude,RDouble latitude,RDouble a,RDouble b,RDouble c,RDouble d);
}
