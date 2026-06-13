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
//! @file      MixedGas.h
//! @brief     Explain this file briefly.
//! @author    Meng Liyuan, Li Peng.

#pragma once
#include "Precision.h"
#include "Gas.h"
#include "Geo_StructGrid.h"
#include "Geo_UnstructGrid.h"
#include "ChemModel.h"

#include <string>
using namespace std;

namespace PHSPACE
{
class DataContainer;

namespace GAS_SPACE
{
const RDouble chem_temp_diff = 1.0e-8;
const int chem_pres_iter_max = 1000;

class ThermodynamicFunction
{
public:
    ThermodynamicFunction();
    ~ThermodynamicFunction();
private:
    int number_of_polynomial_coefficient;
    int number_of_temperature_interval;
    RDouble **polynomial_coefficient;
    RDouble *enthalpyFitCoefficient;
    RDouble **polynomial_coefficientT2;
public:
    int  GetNumberOfPolynomialCoefficient() const { return number_of_polynomial_coefficient; }
    void SetNumberOfPolynomialCoefficient(int number_of_polynomial_coefficient) { this->number_of_polynomial_coefficient = number_of_polynomial_coefficient; }
    void Init(int number_of_temperature_interval, int number_of_polynomial_coefficient);
    RDouble * GetPolynomialCoefficient(int i_temperature_interval) { return polynomial_coefficient[i_temperature_interval]; }
    RDouble * GetPolynomialCoefficientT2(int i_temperature_interval) { return polynomial_coefficientT2[i_temperature_interval]; }
    RDouble * GetEnthalpyFitCoefficient() { return enthalpyFitCoefficient; }
public:
    void ReadPolynomialCoefficient(fstream &file);
    void CreatePolynomialCoefficient(RDouble **curveFitData, RDouble **curveFitDataT2, RDouble *enthalpyFitData);
    void ReadPolynomialCoefficient(DataContainer *cdata);
    void WritePolynomialCoefficient(DataContainer *cdata);
};

class ThermodynamicManager
{
public:
    ThermodynamicManager();
    ~ThermodynamicManager();
private:
    int numberOfSpecies;
    int numberOfTemperatureInterval, numberofpolynomialcoefficient;
    ThermodynamicFunction **thermoDynamicFunction;
    RDouble *temperatureRange;
    RDouble *temperatureStart;
public:
    void Init(int numberOfSpecies);
    void Read(fstream &file);
    void GenerateData(vector<string> &namelist);
    void Read(DataContainer *cdata);
    void Write(DataContainer *cdata);
public:
    RDouble * GetEnthalpyFitCoefficient(int iSpecies);
    RDouble * GetPolynomialCoefficient(int iSpecies, int itemperature);
    RDouble * GetPolynomialCoefficientT2(int iSpecies, int itemperature);
    void GetTemperatureRangeIndex(const RDouble &temperature_dimensional, int &indexOfTemperatureRange);

    RDouble GetTemperatureMin();
    RDouble GetTemperatureMax();

    int GetNumberofPolynomialCoefficient();
    RDouble GetBenchmarkTemperature(int iSpecies);
};

class ReactionRate
{
public:
    ReactionRate();
    ~ReactionRate();
private:
    int numberOfReaction;
    RDouble *afr, *bfr, *cfr, *abr, *bbr, *cbr;
public:
    void Init(int numberOfReaction);
    void Read(fstream &file);
    //! To generate the reaction rate coefficient of curve fitting method.
    //! Param[in]: reationCoef denotes the array of saving the chemical reaction coefficient of each reaction equation.
    void Create(RDouble **reationCoef);
    void Read(DataContainer *cdata);
    void Write(DataContainer *cdata);
    RDouble * GetForwardReactionRateCoefficientA() { return afr; }
    RDouble * GetForwardReactionRateCoefficientB() { return bfr; }
    RDouble * GetForwardReactionRateCoefficientC() { return cfr; }

    RDouble * GetBackwardReactionRateCoefficientA() { return abr; }
    RDouble * GetBackwardReactionRateCoefficientB() { return bbr; }
    RDouble * GetBackwardReactionRateCoefficientC() { return cbr; }
};

class StoichiometricEquation
{
public:
    StoichiometricEquation();
    ~StoichiometricEquation();
private:
    int numberOfReaction, numberOfSpecies;
    int **cvp, **cvn;
    RDouble **cpz;
    int *reactionType;
    bool *isCollisionType;
public:
    void Init(int numberOfReaction, int numberOfSpecies);
    void Read(fstream &file);
    void Create(int **forwardCoef, int **backwardCoef, RDouble **collisionCoef, int **reactionFlag);
    void Read(DataContainer *cdata);
    void Write(DataContainer *cdata);

    int * GetChemicalReactionType() { return reactionType; }
    bool * GetCollisionReactionFlag() { return isCollisionType; }

    int ** GetForwardReactionStoichiometricCoefficient() const { return cvp; }
    int ** GetBackwardReactionStoichiometricCoefficient() const { return cvn; }
    RDouble ** GetReactionStoichiometricThirdBodyCoefficient() const { return cpz; }
};

class CurveFitsMethod
{
public:
    CurveFitsMethod();
    ~CurveFitsMethod();
private:
    int numberOfSpecies;
    //The curve fit coefficients of blottner method.
    RDouble *blotter_ai, *blotter_bi, *blotter_ci, *blotter_di, *blotter_ei;
    //The curve fit coefficients of average collision areas of species couples.
    RDouble *omega_ai, *omega_bi, *omega_ci, *omega_di;
    //The curve fit coefficients of vibrational excitation rate of electron collision-N2.
    RDouble *curvefit_av, *curvefit_bv, *curvefit_cv;
    //The curve fit coefficients of effective collision cross-sectional area between the electron and neutral heavy particle.
    //The sorting sequence of the neutral heavy particle is O, O2, NO, N and N2.
    RDouble *curvefit_aes, *curvefit_bes, *curvefit_ces;
public:
    void Init(int numberOfSpecies, int nGasModel = 0);
    void Read(fstream &file);
    void GenerateData(vector<string> &namelist);
    void Read(DataContainer *cdata);
    void Write(DataContainer *cdata);
public:
    //!Obtain the array of curve fit coefficients that is used to compute the viscosities of species.
    RDouble * GetBlotterCurveFittingCoefficientA() const { return blotter_ai; }
    RDouble * GetBlotterCurveFittingCoefficientB() const { return blotter_bi; }
    RDouble * GetBlotterCurveFittingCoefficientC() const { return blotter_ci; }
    RDouble * GetBlotterCurveFittingCoefficientD() const { return blotter_di; }
    RDouble * GetBlotterCurveFittingCoefficientE() const { return blotter_ei; }
    //! Set the curve fitting data of Blotter viscosity.
    //! @param[in ]: speciesName is an array of saving the names of species.
    //! @param[in ]: nSpecies is the number of the total species.
    //void SetBlotterCurveFittingCoefficients(string *speciesName, int nSpecies);
    //void SetGuptaCurveFittingCoefficients(string *speciesName, int nSpecies);

    //!Obtain the array of curve fit coefficients that is used to compute the average collision areas of species couples.
    RDouble * GetOmegaCurveFittingCoefficientA() const { return omega_ai; }
    RDouble * GetOmegaCurveFittingCoefficientB() const { return omega_bi; }
    RDouble * GetOmegaCurveFittingCoefficientC() const { return omega_ci; }
    RDouble * GetOmegaCurveFittingCoefficientD() const { return omega_di; }

    //!Set the curve fit coefficients of vibrational excitation rate of electron collision-N2.
    void SetNitrogenVibrationalExcitationRate();
    //!Obtain the array of curve fit coefficients that is used to compute the vibrational excitation rate of electron collision-N2.
    RDouble * GetNitrogenVibrationalExcitationRateA() const { return curvefit_av; }
    RDouble * GetNitrogenVibrationalExcitationRateB() const { return curvefit_bv; }
    RDouble * GetNitrogenVibrationalExcitationRateC() const { return curvefit_cv; }

    //!Set the curve fit coefficients of effective collision cross-sectional area between the electron and neutral heavy particle.
    void SetEffectiveCollisionCrossSectionalAreaCurveFits(int nGasModel = 0);
    void SetEffectiveCollisionCrossSectionalAreaCurveFits(vector<string> &namelist);

    //!Obtain the array of curve fit coefficients that is used to compute the effective collision cross-sectional area between the electron and neutral heavy particle.
    RDouble * GetEffectiveCollisionCrossSectionalAreaCurveFitsA() const { return curvefit_aes; }
    RDouble * GetEffectiveCollisionCrossSectionalAreaCurveFitsB() const { return curvefit_bes; }
    RDouble * GetEffectiveCollisionCrossSectionalAreaCurveFitsC() const { return curvefit_ces; }
};

class SchmidtNumber
{
public:
    SchmidtNumber();
    ~SchmidtNumber();
private:
    int numberOfSpecies;
    RDouble *scl, *sct, *oscl, *osct;
public:
    void Init(int numberOfSpecies);
    //void Read(fstream &file);
    void ComputeSpeciesSchmidtNumber(int *ionTypeOfSpecies);
    void Read(DataContainer *cdata);
    void Write(DataContainer *cdata);
public:
    RDouble * GetLaminarSchmidtNumber() const { return scl; }
    RDouble * GetTurbulentSchmidtNumber() const { return sct; }

    RDouble * GetLaminarSchmidtNumberReciprocal() const { return oscl; }
    RDouble * GetTurbulentSchmidtNumberReciprocal() const { return osct; }
};

//! Define the struct to store thermodynamic properties.
typedef struct{
    //! The type is used to identify electron, monatom and polyatom.
    int nType;
    //! The degeneracy.
    int gs0, gs1;
    //! The characteristic temperature of vibration mode.
    //RDouble Tvs;
    VibrationModeData Tvs;
    //! The characteristic temperature of electron mode.
    RDouble Tes;
    //! The enthalpy of formation.
    RDouble hs0;
}Thermo_Param;


class MolecularProperty
{
public:
    MolecularProperty();
    ~MolecularProperty();
private:
    int numberOfSpecies;
    RDouble *molecularWeightDimensional, *molecularWeight, *oMolecularWeightDimensional, *oMolecularWeight;
    RDouble *characteristicTemperatureOfSpecies;
    RDouble *mass_fraction_reference;
    string *name_of_species;
    int    *ionTypeOfSpecies;
    RDouble *collisionCrossSection;
    SchmidtNumber *schmidt_number;
private:
    //! The thermodynamic parameters of species in the mixed gas.
    Thermo_Param *thermodynamicProperties;
public:
    void Init(int numberOfSpecies);
    void Read(fstream &file);
    void GenerateData(vector<string> &namelist, vector<RDouble> &initMassFractions);
    void Read(DataContainer *cdata);
    void Write(DataContainer *cdata);
public:
    //! Obtain the nondimensional molecular weight of the mixture gas.
    RDouble * GetMolecularWeight() const { return molecularWeight; }
    RDouble * GetMolecularWeightReciprocal() const { return oMolecularWeight; }

    //! Obtain the dimensional molecular weight of the mixture gas whose unit is kg/mol.
    RDouble * GetMolecularWeightDimensional() const { return molecularWeightDimensional; }
    RDouble * GetMolecularWeightReciprocalDimensional() const { return oMolecularWeightDimensional; }

    RDouble * GetInitialMassFraction() const { return mass_fraction_reference; }

    string * GetNameOfSpecies() const { return name_of_species; }

    RDouble * GetLaminarSchmidtNumberReciprocal() const;
    RDouble * GetTurbulentSchmidtNumberReciprocal() const;

    //! To obtain the thermodynamic properties of species.
    Thermo_Param * GetThermodynamicProperties() const { return thermodynamicProperties; }

    int * GetIonTypeOfSpecies() const { return ionTypeOfSpecies; }
};

class MixedGas : public Gas
{
public:
    MixedGas();
    ~MixedGas();
private:
    FluidParameter referenceParameterDimensional;    //! Ref
    ThermodynamicManager  *thermodynamicManager;     //! Ppolynomial fitting for Species Energy
    ReactionRate *reactionRate;                      //! kb,kf
    StoichiometricEquation *stoichiometricEquation;  //! Equation coef for matrix
    CurveFitsMethod *CurveFits;                      //! Curve fitting for viscosity,collision section
    MolecularProperty *molecularProperty;            //! Species characteristic parameter
private:
    int nm, nl, numberOfSpecies, numberOfReaction, nchem, ntmodel, nEquation;
    int nTEnergyModel;  //! 0 indicates Energy level method ,  1 is for  curve fitting method.
    int nSpeciesEquation; //! The number of species equation.
    int nElectronIndex; //! The index of electron in the array.
    int nNitrogenIndex; //! The index of nitrogen in the array.
    int nLastSpeciesIndex; //! The index of the key species with the maximum mass fraction in the array.
    int nGasModel;    //! nGasModel = 0 indicates the gas of Earth, nGasModel = 1 is for gas of Mars.
    int nSpeciesLimit ;
    int mTT;
    int mTV;
    int mTE;
    int nEnergyRecycle;
    int nIsChemicalFreeze;
    int nIsSuperCatalytic;
    int nChemcalSourceModified;
    int nChemcalSourceEsMethod;
    int nDensityForWallMethod;

    RDouble electronMass;
    RDouble electronMassDimension;
    RDouble reynoldsReferenceLengthDimensional;
    RDouble coefficientOfStateEquation; //universal gas constant, R
    RDouble oprl, oprt;
    RDouble referenceEnergy;
    RDouble referenceSpecificHeat;
    //! The power of translational-rotational temperature in the Park's V-D(Vibration-Dissociation) coupling model.
    RDouble parkVDPower; //the value is in range of [0.0, 1.0].
    RDouble chemicalRelaxCorf;
    RDouble veTemperatureMin;
    RDouble veTemperatureMinNonDimensional;
    RDouble trTemperatureMin;
    RDouble trTemperatureMinNonDimensional;
    RDouble maxViscous;
    RDouble lnMaximum;
  //  int nMaxStepTemperature, precisionTemperature;
public:
    void GetReferenceParameters(FluidParameter &refParam) { refParam = referenceParameterDimensional; }
    int GetSpeciesEquationNumber() const { return nSpeciesEquation; }
    int GetNitrogenIndex() const { return nNitrogenIndex; }
    int GetElectronIndex() const { return nElectronIndex; }
    void SetSpeciesEquationNumber(int ns) { nSpeciesEquation = ns; }
    void SetNitrogenIndex(int nn) { nNitrogenIndex = nn; }
    void SetElectronIndex(int ne) { nElectronIndex = ne; }
    void SetElectronMolecularWeight(RDouble me) { electronMass = me; }
    RDouble GetElectronMolecularWeight() const { return electronMass; }
    RDouble GetDimensionalElectronMolecularWeight() const { return electronMassDimension; }
    int GetNSEquationNumber() const { return nm; }
    int GetLaminarNumber() const { return nl; }
    int GetTemperatureModel() const { return ntmodel; }
    int GetChemicalType() const { return nchem; }
    int GetnLastSpeciesIndex() const { return nLastSpeciesIndex; }
    RDouble * GetInitMassFraction() const { return molecularProperty->GetInitialMassFraction(); }
    virtual RDouble * GetMolecularWeightDimensional() const { return molecularProperty->GetMolecularWeightDimensional(); }
    RDouble ComputeReferenceTotalEnthalpy();
    int GetmTT() const { return mTT; }
    int GetmTV() const { return mTV; }
    int GetmTE() const { return mTE; }
    int GetnEnergyRecycle() const { return nEnergyRecycle; }
    int GetnChemcalSourceModified() const { return nChemcalSourceModified; }
    int GetGasModelType() { return this->nGasModel; }
    int GetnDensityForWallMethod() const { return nDensityForWallMethod; }
private:
    RDouble *moleFractions;
    RDouble *massFractions;
    RDouble *speciesViscosity;
    RDouble *speciesWeight;
    RDouble *speciesCp;
    RDouble *speciesCv;
    RDouble *speciesEnthalpy;
    RDouble *speciesEnergyTr;
    RDouble *speciesEnergyVb;
    RDouble *speciesEnergyEe;
    RDouble *speciesConductivity;
    RDouble *workSpecies;

    RDouble *dfsdx;
    RDouble *dfsdy;
    RDouble *dfsdz;
    RDouble *fvis;
    Thermo_Energy Thermo_Energy_temparay;    //! Energy parameter.
    
    RDouble *maximumSpecies;    //! The maximum of each Specie.
    RDouble **polynomialCoefT2; //! The interpolation coefficients for the polynomial in the two-temperature model.
private:
    RDouble GetTemperatureMin();
    RDouble GetTemperatureMax();
    Thermo_Energy *GetThermo_Energy_temparay() { return &Thermo_Energy_temparay; }
private:
    void AllocateWorkingSpace();
    void DeAllocateWorkingSpace();
    void InitCommonParameter();
    void InitOtherParameter();
    void ComputeReferenceParameter();
    void ComputeReferenceTsuth();

    void ComputeReferenceViscosityDimensionalChemical();
    void ComputeReferenceSpecificHeatRatio();
    void ComputeReferenceGeneralGasConstant();

    void ComputeReferenceSpecificHeatRatioWithChemical();
    void ComputeReferenceViscosityDimensionalNoChemical();
    void ComputeReferenceViscosityDimensional();
    void ComputeReferenceGasInformation();
    void ComputeReferenceSoundVelocity();
    void ComputeReferenceInitMassFractionInformation();
    void ComputeReferenceMolecularInformation();
    void ComputeReferenceMolecularInformationWithChemical();
    void ComputeReferenceMolecularInformationWithoutChemical();
    void ComputeReferenceVelocity();
    void ComputeReferenceReynoldsNumber();
    void ComputeProperReynoldsNumberForGrid();
    void ComputeOtherProperParameterForGrid();
    RDouble ComputeReferenceReynoldsNumberFromTotalPT(void);
    void ComputeDensityWithReynoldsNumber();
    void ComputePressureInGasStateEquation();
    void InitReferenceParameter();
    void InitNumberOfSpecies();
    void ComputeCoefficientOfStateEquation();
    void Print();

    //! Obtain the T, p, rho and c by the height (Rg=287.053, gamma=1.4 are used). However, only the T and p is useful. 
    //! The other two variables will be recomputed to be compatible with the condition of chemical reactions existing.
    void SetAirInformationByDataBase();

    //! Obtain the T, p, rho and c by the T0, P0, Ma (Rg=287.053, gamma=1.4 are used). However, only the T and p is useful. 
    //!The other two variables will be recomputed to be compatible with the condition of chemical reactions existing.
    void SetAirInformationByExpData();

    //! This function is mainly for the average density of the mixture in chemical reactions, 
    //! The value of mixture density has little difference with that of no chemical reaction, which can be ignored.
    void NormalizeAirInformation();

    void ComputeReferencePrimitive();

    void CreateAllDataClass(int numberOfSpecies, int numberOfReaction);
    void Read(fstream &file);
    void Read(DataContainer *cdata);
    void Write(DataContainer *cdata);

    void SetMaximumSpecies(RDouble *mass_fraction_ref);
public:
    void CompressData(DataContainer *&cdata);
    void DecompressData(DataContainer *cdata);
private:
    RDouble * GetLaminarSchmidtNumberReciprocal() const;
    RDouble * GetTurbulentSchmidtNumberReciprocal() const;
public:
    //! To obtain the laminar Schmidt number of each component.
    void GetLaminarSchmidtNumber(RDouble *Scl) const;
    //! To obtain the turbulent Schmidt number of each component.
    void GetTurbulentSchmidtNumber(RDouble *Sct) const;
    //! To obtain the molecular weight of each component.
    RDouble * GetMolecularWeight() const { return molecularProperty->GetMolecularWeight(); }

    int  GetNumberOfSpecies() const { return numberOfSpecies; }
    void SetNumberOfSpecies(int numberOfSpecies) { this->numberOfSpecies = numberOfSpecies; }
    void InitGlobalParameterOfNSEquation();
    void ReadGasModel(fstream &file);
    //! To Generate the gas model.
    void GenerateGasModel(string chemicalModel);
    void GetNameList(vector<string> &namelist, vector<RDouble> &initMassFraction);

    void InitGasModel();
    //! To import the data of gas model from file.
    void ReadGasModel();

    void GetTemperatureRangeIndex(const RDouble &temperature_dimensional, int &indexOfTemperatureRange);
    RDouble GetCoefficientOfStateEquation() { return coefficientOfStateEquation; }
    string * GetNameOfSpecies() const;

    void GetElectronMolecularWeightAndMassFraction(RDouble *massFraction, RDouble &ce, RDouble &Me);

    //! Compute pressure when the total interal energy e is known.
    //! @param[in ]: primitiveVars is the array that stores the primitive variables.
    //! @param[in ]: gama is the specific heat ratio.
    //! @param[in ]: internalEnergy denotes the total internal energy that equals to the density multiplies the specific internal energy.
    //! @param[in ]: temperature denotes the translation-rotation temperature in the one-temperature model.
    //! @param[out]: pressure is the returned value that denotes the pressure of the mixed gas.
    void GetPressure(RDouble *primitiveVars, RDouble gama, RDouble internalEnergy, RDouble temperature, RDouble &pressure);

    //! Compute temperatures and pressure in the next time-advancing step.
    //! @param[in ]: primitiveVars denotes the array that stores the primitive variables.
    //! @param[in ]: internalEnergy denotes the total internal energy that equals to the density multiplies the specific internal energy.
    //! @param[in/out]: temperature denotes the array of temperatures in the current time-advancing step and next time-advancing step.
    //! As input parameter, temperature indicates the temperature of the current time-advancing step. Meanwhile,
    //! it denotes the temperature of the next time-advancing step as output parameter.
    //! @param[out]: pressure denotes the pressure of the mixed gas.
    void GetTemperatureAndPressure(RDouble *primitiveVars, RDouble internalEnergy, RDouble *temperature, RDouble &pressure);

    //! To obtain the non-dimensional enthalpy of each species in accordance with the temperature of mixture gas.
    //! @param[in ]: temperature denotes the non-dimensional temperature of the mixture gas.
    //! @param[out]: speciesEnthaly saves the non-dimensional static enthalpy of each species.
    void ComputeSpeciesEnthalpy(RDouble temperature, RDouble *speciesEnthaly);

    //! To obtain the dimensional total enthalpy of each species in accordance with the temperature of mixture gas.
    //! @param[in ]: temperature denotes the dimensional temperature of the mixture gas.
    //! @param[out]: speciesEnthaly saves the dimensional total enthalpy of each species.
    void ComputeSpeciesEnthalpyDimensional(RDouble temperature, RDouble *speciesEnthaly);

    void ComputeSpeciesEnergy(RDouble temperature, RDouble *speciesEnergy);
    void ComputeSpeciesEnergyDimensional(RDouble temperature, RDouble *speciesEnergy);

    //! To obtain the non-dimensional molecular weight of the mixture gas in accordace with the mass fraction of each species.
    //! @param[in ]: massFraction is an array saving the mass fraction of each species.
    RDouble GetMixedGasMolecularWeight(RDouble *massFraction);

    // !To obtain the non-dimensional molecular weight of the mixture gas eliminating electron.
    // !@param[in ]: massFraction is an array saving the mass fraction of each species.
    // !@param[out]: ceDivideMe is the value of the electron mass fraction dividing the electron molecular weight.
    RDouble GetMixedGasMolecularWeightReciprocalWithoutElectron(RDouble *massFraction, RDouble &ceDivideMe);

    //! To obtain the non-dimensional universal gas constant called R = 8.314J/(mol*K).
    RDouble GetUniversalGasConstant();

    //! To obtain the static enthalpy of the mixed gas.
    //! @param[in ]: massFraction is an array saving the mass fraction of each species.
    //! @param[in ]: temperature denotes the non-dimensional temperature of the mixture gas.
    //! @param[out]: return the non-dimensional value of static enthalpy of mixed gas.
    RDouble GetMixtureGasEnthalpy(RDouble *massFraction, RDouble temperature);

    //! To obtain the non-dimensional constant pressure specific heat called cp according to the mass fraction of each species.
    //! @param[in ]: massFraction is an array saving the mass fraction of each species.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! @param[out]: mixedcp denotes the non-dimensional constant pressure specific heat of the mixture gas.
    void ComputeConstantPressureSpecificHeatByMassFraction(RDouble *massFraction, RDouble temperature, RDouble &mixedcp);

    //! To obtain the non-dimensional constant pressure specific heat called cp of each species by fitting formula.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! @param[out]: cpSpecies is an array saving the non-dimensional constant pressure specific heat of each species.
    void ComputeSpeciesConstantPressureSpecificHeat(RDouble temperature, RDouble *cpSpecies);

    //! To obtain the dimensional constant pressure specific heat called cp[J/(kg*K)] of each species by fitting formula.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! @param[out]: cpSpeciesDimensional is an array saving the dimensional constant pressure specific heat of each species.
    void ComputeSpeciesConstantPressureSpecificHeatDimensional(RDouble temperature, RDouble *cpSpeciesDimensional);

    //! To obtain the non-dimensional constant volume specific heat called cv according to the mass fraction of each species.
    //! @param[in ]: massFraction is an array saving the mass fraction of each species.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! @param[out]: mixedGasCv denotes the non-dimensional constant volume specific heat of the mixture gas.
    void ComputeConstantVolumeSpecificHeatByMassFraction(RDouble *massFraction, RDouble temperature, RDouble &mixedGasCv);

    //! To obtain the dimensional constant volume specific heat called cv[J/(kg*K)] of each species by fitting formula.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! @param[out]: cvSpeciesDimensional is an array saving the dimensional constant volume specific heat of each species.
    void ComputeSpeciesConstantVolumeSpecificHeatDimensional(RDouble temperature, RDouble *cvSpeciesDimensional);

    //! To obtain the non-dimensional constant volume specific heat called cv of each species by fitting formula.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! @param[out]: cvSpecies is an array saving the non-dimensional constant volume specific heat of each species.
    void ComputeSpeciesConstantVolumeSpecificHeat(RDouble temperature, RDouble *cvSpecies);

    //! To obtain the non-dimensional viscosity of each species according to the Blotter fitting formula.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! @param[out]: viscosityOfSpecies is an array saving the non-dimensional viscosity of each species.
    void ComputeOneTemperatureModelSpeciesViscosity(RDouble temperature, RDouble *viscosityOfSpecies);

    //! To obtain the dimensional viscosity of each species according to the Blotter fitting formula.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! @param[out]: viscosityOfSpeciesDimensional is an array saving the dimensional viscosity of each species [N*s/m2]
    void ComputeOneTemperatureModelSpeciesViscosityDimensional(RDouble temperature, RDouble *viscosityOfSpeciesDimensional);

    void ComputeSpeciesEnthalpy(RDouble temperature, RDouble *speciesCv, RDouble *speciesE, RDouble *speciesH);

    void ComputeMixtureByMassFraction(RDouble *massFraction, RDouble *mixtureOfSpecies, RDouble &mixture);
    void ComputeMixtureByPrimitive(RDouble *prim, RDouble *mixtureOfSpecies, RDouble &mixture);
    void ComputeMixtureCoefficientByWilkeFormula(RDouble *moleFractionOfSpecies, RDouble *chem_var, RDouble *chem_phi);

    void ComputeMixtureByWilkeFormula(RDouble *moleFractionOfSpecies, RDouble *mixtureOfSpecies, RDouble *phiOfSpecies, RDouble &mixture);

    //! To obtain the mass fraction of each species according to other fraction of each species.
    //! @param[in ]: fractionOfSpecies is an array of saving the fraction of each specie,such as mole fraction,volume fraction.
    void MassFractionConversion(RDouble *fractionOfSpecies);
    //! To obtain the mole fraction of each species according to the mass fraction of each species.
    //! @param[in ]: massFractionOfSpecies is an array of saving the mass fraction of each species.
    //! @param[out]: moleFractionOfSpecies is an array of saving the mole fraction of each species.
    void ComputeMoleFractionByMassFraction(RDouble *massFractionOfSpecies, RDouble *moleFractionOfSpecies);

    //! To compute the non-dimensional temperature of the mixture gas in accordance with primary variables saving in the array of prim.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[out]: temperature is the non-dimensional value of temperature.
    void ComputeTemperature(RDouble *primitiveVariables, RDouble &temperature);

    //! To compute the non-dimensional temperature of the mixture gas in accordance with primary variables saving in the array of prim.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[out]: transRotationTemperature is the non-dimensional value of translation-rotation temperature.
    //! @param[out]: vibrationTemperature is the non-dimensional value of vibration temperature.
    //! @param[out]: electronTemperature is the non-dimensional value of electron temperature.
    void GetTemperature(RDouble *primitiveVariables, RDouble &transRotationTemperature, RDouble &vibrationTemperature, RDouble &electronTemperature);
    void GetTemperatureR(RDouble *primitiveVariables, RDouble &transRotationTemperature, RDouble &vibrationTemperature, RDouble &electronTemperature);

    //! To obtain the mass fraction of the last species with the label ns-1.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    void NormalizePrimitive(RDouble *primitiveVariables);

    //! Compute the molecular weight reciprocal of the mixed gas.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[out]: massReciprocal denotes the reciprocal of the molecular weight.
    void ComputeMolecularWeightReciprocal(RDouble *primitiveVariables, RDouble &massReciprocal);

    void ComputeMolecularWeightReciprocalDimensional(RDouble *prim, RDouble &omavDimensional);

    //! Compute the specific heat at constant pressure using the Eucken method.
    //! @param[out]: cpSpecies is an array to saving the specific heat at constant pressure of each species.
    void ComputeSpeciesConstantPressureSpecificHeatByEuckenFormula(RDouble *cpSpecies);

    //! Compute the static enthalpy with the primitive variables.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: gama is the specific heat ratio.
    //! @param[out]: enthalpy denotes the static enthalpy.
    //! @param[in ]: temperatures stores the translation-rotation temperature, vibration temperature and the electron temperature.
    void ComputeEnthalpyByPrimitive(RDouble *primitiveVariables, RDouble gama, RDouble &enthalpy, RDouble *temperatures);

    //! To obtain the total energy E = e + 1/2*V2.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables, gama is the specific heat ratio of the mixture gas.
    //! @param[in ]: vibrationTemperature is non-dimensional value of vibration temperature.
    //! @param[in ]: electronTemperature is non-dimensional value of electron  temperature.
    //! @param[OUT]: totalEnergy denotes the total internal energy of the mixed gas.
    void ComputeInternalEnergy(RDouble *primitiveVariables, RDouble gama, RDouble vibrationTemperature, RDouble electronTemperature, RDouble &totalEnergy);

    //! To compute the total enthalpy and the difference of the static enthalpy.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: gama is the specific heat ratio.
    //! @param[in ]: deltaQ is an array to saving the differences of the primitive variables.
    //! @param[out]: deltaEnthalpy denotes the difference of static enthalpy.
    //! @param[out]: totalEnthalpy denotes the total enthalpy.
    void ComputeDHAndTotalEnthalpy(RDouble *primitiveVariables, RDouble &gama, RDouble *deltaQ, RDouble &deltaEnthalpy, RDouble &totalEnthalpy);

    //! This function is used for perfect gas and single temperature model.
    //! To compute the total enthalpy and the variable dh in the function MXDQ_STD(). dh=b2, b2 denotes the coefficient \n
    //! of the vector M*dQ which can be referred to the forumla (A.7) and (A.8) in the appendix A of the PHengLEI Theory manual.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: gama is the specific heat ratio.
    //! @param[in ]: deltaQ is an array to saving the differences of the primitive variables.
    //! @param[out]: totalEnthalpy denotes the total enthalpy.
    //! @param[out]: deltaEnthalpy denotes the difference of static enthalpy.
    void ComputeTotalEnthalpyAndDH(RDouble *primitiveVariables, RDouble &gama, RDouble *deltaQ, RDouble &totalEnthalpy, RDouble &deltaEnthalpy);

    //! To obtain the non-dimensional viscosity of the mixture gas according to the primary variables.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: transRotationTemperature is the non-dimensional translation-rotation temperature of the mixture gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[in ]: suthTemperature denotes the temperature value in the Sutherland fomula.
    //! @param[out]: viscosity is the non-dimensional viscosity of the mixture gas [N*s/m2].
    void ComputeViscosityByPrimitive(RDouble *primitiveVariables, RDouble &transRotationTemperature, RDouble &electronTemperature, RDouble &suthTemperature, RDouble &viscosity);

    //! To Compute the viscosity, heat conductivity and mass diffusion coefficient.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: temperature is the array of temperature for each temperature model.
    //! @param[out]: viscosity is the non-dimensional viscosity of the mixture gas.
    //! @param[out]: heatConductivity denotes the array of heat conductivity for each mode.
    //! @param[out]: rhoDs denotes the array of mass diffusion coefficients for each species.
    void ComputeTransportCoefficients(RDouble *primitiveVariables, RDouble *temperature, RDouble visTurbulence, RDouble &viscosity, RDouble *heatConductivity, RDouble *rhoDs);
    void ComputeTransportCoefficientsR(RDouble *primitiveVariables, RDouble *temperature, RDouble *speciesCps, RDouble *speciesCvvs, RDouble *speciesCves, RDouble visTurbulence, RDouble &viscosity, RDouble *heatConductivity, RDouble *rhoDs);

    //! To obtain the non-dimensional viscosity of the mixture gas according to the mass fraction of each species.
    //! @param[in ]: density denotes the non-dimensional value of density of the mixed gas.
    //! @param[in ]: massFractionOfSpecies is an array of saving the mass fraction of each species.
    //! @param[in ]: transRotationTemperature is the non-dimensional translation-rotation temperature of the mixture gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[in ]: suthTemperature denotes the temperature value in the Sutherland fomula.
    //! @param[out]: viscosity is the non-dimensional viscosity of the mixture gas [N*s/m2].  
    void ComputeViscosity(RDouble density, RDouble *massFractionOfSpecies, RDouble &transRotationTemperature, RDouble &electronTemperature, RDouble &suthTemperature, RDouble &viscosity);

    //! To compute the heat conductivity of each species using the Eucken formula.
    //! @param[in ]: viscosity is array of saving the species viscosity.
    //! @param[in ]: cps is the array of saving the species values of specific heat at constant pressure.
    //! @param[out]: conductivity is an array of saving the heat conductivity of each species.
    void ComputeSpeciesHeatConductivityByEuckenFormula(RDouble *viscosity, RDouble *cps, RDouble *conductivity);

    //! To compute heat conductivity of mixture gas using Eucken formula and Wassilewa relation.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! The function returns the heat conductivity of mixture gas.
    RDouble ComputeMixtureGasHeatConductivityByWassilewaFormula(RDouble *primitiveVariables, RDouble temperature);

    //! To obtain the partition function values based on the targets, e.g conductivity or viscosity.
    //! @param[in ]: target is array of viscosity or heat conductivity of each species.
    //! @param[in ]: speciesMass is array of molecular weight.
    //! @param[in ]: moleFraction is array of mole fraction.
    //! @param[out]: phi is an array of saving the partition function value of each species.
    void GetPartitionFunctionValues(RDouble *target, RDouble *speciesMass, RDouble *moleFraction, RDouble *phi);

    //! Compute the source terms of chemical reaction equations.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[out]: chemicalSourceTerms returns the source terms of chemical reaction equations in the specified control cells.
    //! @param[in ]: iCell denotes the number of the specified control cell.
    void ChemicalSource(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, int iCell);
    //! Compute the source terms of vibration-electron energy equations.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[in ]: chemicalSourceTerms stores the source terms of chemical reaction equations in the specified control cells.
    //! @param[out]: energySourceTerms returns the source terms of vibration-electron energy equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    void EnergySource(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble **energySourceTerms, int nCellNumber);

    //! Compute the source term of vibration and electron energy equation.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: omega is the source term of chemical reaction equation.
    //! @param[in ]: temperature is non-dimensional value of temperatures.
    //! @param[out]: return the dimensional value of source term of vibration and electron energy equation.
    RDouble ComputeVibrationAndElectronEnergySourceTerm(RDouble *primitiveVariables, RDouble *omega, RDouble *temperature);

    //! Compute the source term of vibration energy equation.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: omega is the source term of chemical reaction equation.
    //! @param[in ]: temperature is non-dimensional value of temperatures.
    //! @param[out]: return the dimensional value of source term of vibration energy equation.
    RDouble ComputeVibrationEnergySourceTerm(RDouble *primitiveVariables, RDouble *omega, RDouble *temperature);

    //! Compute relaxation time of translation and vibration.The parameters are dimensional values.
    //! @param[in ]: thermodynamicProperties saves the thermodynamic properties of species.
    //! @param[in ]: massFraction saves the mass fractions of species.
    //! @param[in ]: speciesMass saves the molecular weight of species.
    //! @param[in ]: temperature, pressure, and density denote the temperature of translation and rotation, pressure and density, respectively.
    //! @param[out]: tvs saves the relaxation time of translation and vibration of species.
    void GetTranslationAndVibrationRelaxationTime(Thermo_Param *thermodynamicProperties, RDouble *massFraction, RDouble *speciesMass, 
                                                  RDouble temperature, RDouble pressure, RDouble density, RDouble *tvs);
    void GetTranslationAndVibrationRelaxationTimeR(Thermo_Param *thermodynamicProperties, RDouble *massFraction, RDouble *speciesMass, 
                                                  RDouble temperature, RDouble pressure, RDouble density, RDouble *tvs);

    //! Compute relaxation time of electron and vibration.The parameters are dimensional values.
    //! @param[in ]: thermodynamicProperties saves the thermodynamic properties of species.
    //! @param[in ]: electronPressure is the pressure of electron.
    //! @param[in ]: electronTemperature is the temperature of electron.
    //! @param[out]: tes saves the relaxation time of electron and vibration of species.
    void GetElectronAndVibrationRelaxationTime(Thermo_Param *thermodynamicProperties, RDouble electronPressure, RDouble electronTemperature, RDouble *tes);

    //! Compute the source term of electron energy equation.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: omega is the source term of chemical reaction equation.
    //! @param[in ]: temperature is non-dimensional value of temperatures.
    //! @param[out]: return the dimensional value of source term of vibration energy equation.
    RDouble ComputeElectronEnergySourceTerm(RDouble *primitiveVariables, RDouble *omega, RDouble *temperature);

    //! Compute the effective collision cross-sectional area of the electron and the s-th heavy species(the element except for electron).
    //! @param[in ]: thermodynamicProperties saves the thermodynamic properties of species.
    //! @param[in ]: electronTemperature is the dimensional value of temperature of electron.
    //! @param[in ]: electronNumberDensity denotes the number density of electron.
    //! @param[out]: sigma saves the effective collision cross-sectional areas of each heavy species.
    void GetEffectiveCollisionCrossSectionalArea(Thermo_Param *thermodynamicProperties, RDouble electronTemperature, RDouble electronNumberDensity, RDouble *sigma);

    void ChemicalSpectrumRadius(RDouble **q, RDouble **t, RDouble *vol, RDouble **src, int iCell);
    //! To obtain the approximating polynomial coefficients from the polynomial curve fits in five temperature ranges. \n
    //! The approximating polynomial is input for evaluating the species specific heats and enthalpies, and the temperature\n
    //! boundary is smoothed by linearly averaging the polynomial cofficients.
    //! @param[in ]: nSpecies denotes the species index in the collection.
    //! @param[in ]: temperature is the dimensional temperature value.
    //! @param[out]: polynomialCoef is an array of saving the polynomial coefficients.
    void GetLinearlyAveragedPolynomialCoefficients(int nSpecies, RDouble temperature, RDouble *polynomialCoef);
    void GetLinearlyAveragedPolynomialCoefficientsT2(RDouble temperature, RDouble **polynomialCoef);
    void GetLinearlyAveragedPolynomialCoefficientsIndex(RDouble temperature, RDouble &coef, int &nIndex1, int &nIndex2);
    void GetLinearlyAveragedPolynomialCoefficientsByIndex(int nSpecies, RDouble *polynomialCoef, RDouble coef, int nIndex1, int nIndex2);

    //! Compute chemical source term and its Jacobian matrix. The first function is batch processing of data, the second is single processing of data.
    void ComputeChemicalSourceDerivatives(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber);
    void ComputeChemicalSourceJacobian(RDouble *primitiveVariables, RDouble *temperatures, RDouble cellVolume, RDouble *chemicalSourceTerms, RDouble **sourceDerivatives);
    void ComputeChemicalSourceAndJacobian(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber);
    void ComputeChemicalSourceAndJacobianD3(int *Inzone, RDouble4D &primitiveVariables, RDouble4D &temperatures, RDouble4D &speciesEs, RDouble4D &speciesCvs, RDouble3D &totalCv, RDouble3D &cellVolume, RDouble4D &chemicalSourceTerms, RDouble5D &sourceDerivatives);

    //! Compute the partial derivatives of source terms of chemical reaction equations without electron.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[out]: chemicalSourceTerms returns the source terms of chemical reaction equations in the specified control cells.
    //! @param[out]: sourceDerivatives returns the elements of the Jacobian matrix in the source terms of chemical reaction equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    void ComputeChemicalSourceDerivatives1(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber);

    //! Compute the partial derivatives of source terms of chemical reaction equations, it is used for the reactions including electron.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[out]: chemicalSourceTerms returns the source terms of chemical reaction equations in the specified control cells.
    //! @param[out]: sourceDerivatives returns the elements of the Jacobian matrix in the source terms of chemical reaction equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    void ComputeChemicalSourceDerivatives2(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber);

    //! Compute the partial derivatives of source terms of chemical reaction equations. The function is used for two-temperature model.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[in ]: chemicalSourceTerms stores the source terms of chemical reaction equations in the specified control cells.
    //! @param[in ]: chemicalSourceDerivatives stores the partial derivatives of the source terms of chemical reaction equations in the specified control cells.
    //! @param[out]: energySourceDerivatives returns the partial derivatives of the source terms of energy equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    void ComputeEnergySourceDerivatives1(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
                                        RDouble ***chemicalSourceDerivatives, RDouble ***energySourceDerivatives, int nCellNumber);

    //! Compute the partial derivatives of source terms of vibration-electron energy equations. The function is used for two-temperature model.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[in ]: chemicalSourceTerms stores the source terms of chemical reaction equations in the specified control cells.
    //! @param[in ]: chemicalSourceDerivatives stores the partial derivatives of the source terms of chemical reaction equations in the specified control cells.
    //! @param[out]: energySourceDerivatives returns the diagonal element of partial derivatives to the source terms of energy equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    void ComputeEnergySourceDiagonalDerivatives1(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
                                                 RDouble ***chemicalSourceDerivatives, RDouble **energySourceDerivatives, int nCellNumber);

    //! Compute the partial derivatives of source terms of vibration-electron energy equations. The function is used for three-temperature model.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[in ]: chemicalSourceTerms stores the source terms of chemical reaction equations in the specified control cells.
    //! @param[in ]: chemicalSourceDerivatives stores the partial derivatives of the source terms of chemical reaction equations in the specified control cells.
    //! @param[out]: energySourceDerivatives returns the diagonal element of partial derivatives to the source terms of energy equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    void ComputeEnergySourceDiagonalDerivatives2(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
                                                 RDouble ***chemicalSourceDerivatives, RDouble **energySourceDerivatives, int nCellNumber);

    //! Compute energy source term and its Jacobian. The function is for two-temperature model, it is equal to ComputeEnergySourceDiagonalDerivatives1() which is the batch processing.
    void ComputeEnergySourceJacobianT2(RDouble *primitiveVariables, RDouble *temperatures, RDouble cellVolume, RDouble *chemicalSourceTerms, RDouble **sourceDerivatives);

    //! Compute energy source term and its Jacobian. The function is for three-temperature model, it is equal to ComputeEnergySourceDiagonalDerivatives2() which is the batch processing.
    void ComputeEnergySourceJacobianT3(RDouble *primitiveVariables, RDouble *temperatures, RDouble cellVolume, RDouble *chemicalSourceTerms, RDouble **sourceDerivatives);

    void ComputeEnergySourceAndJacobianT2(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber);

    void ComputeEnergySourceAndJacobianT3(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber);

    void ComputeEnergySourceAndJacobianT2D3(int *nCell, RDouble4D &primitiveVariables, RDouble4D &temperatures, RDouble4D &speciesCvvs, RDouble4D &speciesCves, RDouble4D &speciesEvs, RDouble4D &speciesEes, RDouble3D &totalCvtr, RDouble3D &totalCvv, RDouble3D &totalCve, RDouble3D &cellVolume, RDouble4D &chemicalSourceTerms, RDouble5D &sourceDerivatives);

    void ComputeEnergySourceAndJacobianT3D3(int *nCell, RDouble4D &primitiveVariables, RDouble4D &temperatures, RDouble4D &speciesCvvs, RDouble4D &speciesCves, RDouble4D &speciesEvs, RDouble4D &speciesEes, RDouble3D &totalCvtr, RDouble3D &totalCvv, RDouble3D &totalCve, RDouble3D &cellVolume, RDouble4D &chemicalSourceTerms, RDouble5D &sourceDerivatives);

public:

    //! The primitive variables are transfered to the conservative variables.
    //! @param[in ]: q prim indicates the primitive variables.
    //! @param[in ]: gama denotes the specific heat ratio, which is valid in the perfect gas.
    //! @param[in ]: Tv denotes the vibrational temperature(The non-dimensional value).
    //! @param[in ]: Te denotes the electron temperature(The non-dimensional value).
    //! @param[out]: indicates the conservative variables.
    void Primitive2Conservative(RDouble *prim, RDouble gama, RDouble Tv, RDouble Te, RDouble *q);
    void Primitive2ConservativeR(RDouble *prim, RDouble gama, RDouble Tv, RDouble Te, RDouble staticE, RDouble *q);

    //! The conservative variables are transfered to the primitive variables.
    //! @param[in ]: q indicates the conservative variables.
    //! @param[in ]: gama denotes the specific heat ratio, which is valid in the perfect gas.
    //! @param[out]: prim indicates the primitive variables.
    //! @param[in/out]: temperature denotes the array of temperatures in the current time-advancing step and next time-advancing step.
    //! As input parameter, temperature indicates the temperature of the current time-advancing step. Meanwhile,
    //! It denotes the temperature of the next time-advancing step as output parameter.
    void Conservative2Primitive(RDouble *q, RDouble gama, RDouble *prim, RDouble *temperature);
    void Conservative2PrimitiveR(RDouble *q, RDouble gama, RDouble *prim, RDouble *temperature);

    void GetSpecificHeatRatio(RDouble *prim, RDouble &gama);

    //! Compute the specific heat ratio and the temperatures.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[out]: gama denotes the specific heat ratio of the mixed gas.
    //! @param[out]: transRotationTemperature denotes the translation-rotation temperature of the mixed gas.
    //! @param[out]: vibrationTemperature denotes the vibration temperature of the mixed gas.
    //! @param[out]: electronTemperature denotes the electron temperature of the mixed gas.
    void GetSpecificHeatRatioAndTemperatute(RDouble *primitiveVariables, RDouble &gama, RDouble &transRotationTemperature, RDouble &vibrationTemperature, RDouble &electronTemperature);

    void ComputeDensityDiffusionAndKCP(RDouble temperature, RDouble *primface, RDouble mul, RDouble mut, RDouble *rho_ds_face, RDouble *hintSpeciesOfFace, RDouble &kcp);
    void GetAirInfo(const RDouble &height, RDouble &temperature, RDouble &pressure, RDouble &density, RDouble &soundspeed);

    //! To obtain the coefficients of computing the M*dQ in the function MXDQ_STD_LP(), added by LiPeng on Jan. 10, 2019.
    //! @param[in ]: primitiveVars is the array that stores the primitive variables of the current cell.
    //! @param[out]: alpha, beta and speciesBeta are referred to the fomular (A.5) in the appendix A of the PHengLEI Theory manual.
    void ComputeCoefficientInMXDQ(RDouble *primitiveVars, RDouble &alpha, RDouble &beta, RDouble *speciesBeta);
    void ComputeCoefficientInMXDQR(RDouble *primitiveVars, RDouble trTemperature, RDouble squareVelocity, RDouble *speciesCvs, RDouble totalCv, RDouble &alpha, RDouble &beta, RDouble *speciesBeta);

    //!To obtain the coefficients of computing the M*dQ in the function MVXDQ_STD_LP(), added by LiPeng on Jan. 21, 2019.
    //!K=visl/Prl + vist/Prt denotes the thermal conductivity coefficient, rD=visl/Scl + vist/Sct is the species diffusion coefficient,\n
    //!Ds is an array of saving the diffusion coefficient of each species.
    void ComputeCoefficientInMVXDQ(RDouble *prim, RDouble K, RDouble rD, RDouble *theta_s, RDouble *phi_s, RDouble *Ds);

public:

    //! Compute the mass fractions of Earth mixture for fully catalytic wall condition.
    //! @param[in ]: fs is array of mass fractions of the components in the multi-species.
    //! @param[out]: fcw is array of mass fractions under fully catalytic condition.
    void GetEarthFullyCatalyticMassFraction(RDouble *fs, RDouble *fcw);

    //! Compute the mass fractions of Mars mixture for fully catalytic wall condition.
    //! @param[in ]: fs is array of mass fractions of the components in the multi-species.
    //! @param[out]: fcw is array of mass fractions under fully catalytic condition.
    void GetMarsFullyCatalyticMassFraction(RDouble *fs, RDouble *fcw);

    //! To find the chemical element in the species, if found, return the coefficient of the element.
    //! @param[in ]: species_name denotes the name of the species.
    //! @param[out]: element_name denotes the name of the specified element.
    //! @param[out]: coef is the stoichiometric coefficient of the element.
    bool FindChemicalElement(string species_name, char element_name, int &coef);

    //! To obtain the index of the specified species in the sequence in accordance with its name.
    //! @param[in ]: species_name denotes the name of the species.
    //! @param[out]: the function returns the index of the species if found, otherwise, returns value of -1.
    int GetSpeciesIndex(string species_name);

public:

    //! To compute the translation and rotation specific heat at constant volume.
    //! @param[out]: speciesTransRotationCv is an array saving the translation and rotation specific heat at constant volume of each species.
    void ComputeTranslationAndRotationSpecificHeatAtConstantVolume(RDouble *speciesTransRotationCv);

    //! To compute the translation and rotation specific heat at constant volume.
    //! @param[out]: speciesTransRotationCv is an array saving the dimensional value of translation and rotation specific heat at constant volume of each species.
    void ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(RDouble *speciesTransRotationCv);

    //! To compute the the translation and rotation energy.
    //! @param[in ]: transRotationTemperature is the non-dimensional translation-rotation temperature of the mixture gas.
    //! @param[out]: speciesTransRotationEnergy is an array saving the translation and rotation energy of each species.
    void ComputeTranslationAndRotationEnergy(RDouble transRotationTemperature, RDouble *speciesTransRotationEnergy);

    //! To compute the vibration specific heat at constant volume.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[out]: speciesVibrationCv is an array saving the vibration specific heat at constant volume of each species.
    void ComputeVibrationSpecificHeatAtConstantVolume(RDouble vibrationTemperature, RDouble *speciesVibrationCv);

    //! To compute the vibration specific heat at constant volume.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[out]: speciesVibrationCv is an array saving the dimensional value of vibration specific heat at constant volume of each species.
    void ComputeDimensionalVibrationSpecificHeatAtConstantVolume(RDouble vibrationTemperature, RDouble *speciesVibrationCv);

    //! To compute the dimensional value of the vibration energy (J/kg).
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[out]: speciesVibrationEnergy is an array saving the vibration energy of each species.
    void ComputeDimensionalVibrationEnergy(RDouble vibrationTemperature, RDouble *speciesVibrationEnergy);

    //! To compute the vibration energy.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[out]: speciesVibrationEnergy is an array saving the vibration energy of each species.
    void ComputeVibrationEnergy(RDouble vibrationTemperature, RDouble *speciesVibrationEnergy);

    //! To compute the electron specific heat at constant volume.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[out]: speciesElectronCv is an array saving the electron specific heat at constant volume of each species.
    void ComputeElectronSpecificHeatAtConstantVolume(RDouble electronTemperature, RDouble *speciesElectronCv);

    //! To compute the electron specific heat at constant volume.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[out]: speciesElectronCv is an array saving the dimensional value of electron specific heat at constant volume of each species.
    void ComputeDimensionalElectronSpecificHeatAtConstantVolume(RDouble electronTemperature, RDouble *speciesElectronCv);

    //! To compute the dimensional value of electron energy (J/kg).
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[out]: speciesElectronEnergy is an array saving the electron energy of each species.
    void ComputeDimensionalElectronEnergy(RDouble electronTemperature, RDouble *speciesElectronEnergy);

    //! To compute the electron energy.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[out]: speciesElectronEnergy is an array saving the electron energy of each species.
    void ComputeElectronEnergy(RDouble electronTemperature, RDouble *speciesElectronEnergy);

    //! To compute the enthalpies of species.
    //! @param[in ]: transRotationTemperature is the non-dimensional translation-rotation temperature of the mixture gas.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[out]: speciesEnthalpy is an array saving the enthalpies of each species.
    void ComputeSpeciesEnthalpy(RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature, RDouble *speciesEnthalpy);

    //! To compute the translation-rotation specific heat at constant volume of the mixed gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the total translation-rotation specific heat at constant volume.
    RDouble GetMixedGasTranslationAndRotationSpecificHeatAtConstantVolume(RDouble *massFraction);

    //! To compute the translation-rotation specific heat at constant volume of the mixed gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the dimensional value of total translation-rotation specific heat at constant volume.
    RDouble GetMixedGasDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(RDouble *massFraction);

    //! To compute the vibration specific heat at constant volume of the mixed gas.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the total vibration specific heat at constant volume.
    RDouble GetMixedGasVibrationSpecificHeatAtConstantVolume(RDouble vibrationTemperature, RDouble *massFraction);

    //! To compute the vibration specific heat at constant volume of the mixed gas.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the dimensional value of total vibration specific heat at constant volume.
    RDouble GetMixedGasDimensionalVibrationSpecificHeatAtConstantVolume(RDouble vibrationTemperature, RDouble *massFraction);

    //! To compute the electron specific heat at constant volume of the mixed gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the total electron specific heat at constant volume.
    RDouble GetMixedGasElectronSpecificHeatAtConstantVolume(RDouble electronTemperature, RDouble *massFraction);

    //! To compute the electron specific heat at constant volume of the mixed gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the dimensional value of total electron specific heat at constant volume.
    RDouble GetMixedGasDimensionalElectronSpecificHeatAtConstantVolume(RDouble electronTemperature, RDouble *massFraction);

    //! To compute the translation-rotation energy of the mixed gas.
    //! @param[in ]: transRotationTemperature is the non-dimensional translation-rotation temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the total translation-rotation energy.
    RDouble GetMixedGasTranslationAndRotationEnergy(RDouble transRotationTemperature, RDouble *massFraction);

    //! To compute the dimensional value of vibration energy of the mixed gas.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the total vibration energy.
    RDouble GetMixedGasDimensionalVibrationEnergy(RDouble vibrationTemperature, RDouble *massFraction);

    //! To compute the vibration energy of the mixed gas.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the total vibration energy.
    RDouble GetMixedGasVibrationEnergy(RDouble vibrationTemperature, RDouble *massFraction);

    //! To compute the dimensional value of electron energy of the mixed gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the total electron energy.
    RDouble GetMixedGasDimensionalElectronEnergy(RDouble electronTemperature, RDouble *massFraction);

    //! To compute the electron energy of the mixed gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the total electron energy.
    RDouble GetMixedGasElectronEnergy(RDouble electronTemperature, RDouble *massFraction);

    //! To compute the total enthalpy of the mixed gas.
    //! @param[in ]: transRotationTemperature is the non-dimensional translation-rotation temperature of the mixture gas.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[out]: return the total enthalpy of the mixed gas.
    RDouble ComputeMixedGasEnthalpy(RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature, RDouble *massFraction);

    //! To compute the total enthalpy of the mixed gas.
    //! @param[in ]: primitiveVars indicates the array of primitive variables.
    //! @param[in ]: temperatures indicates the array of temperatures.
    RDouble ComputeMixedGasEnthalpy(RDouble *primitiveVars, RDouble *temperatures);

    //! To compute the total enthalpy of formation.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species.
    //! @param[out]: return the non-dimensional value of the formation enthalpy.
    RDouble GetTotalFormationEnthalpy(RDouble *massFraction);

    //! To compute translation-rotation temperature where the pressure is known.
    //! @param[in ]: primitiveVariables is an array of saving the primary variables.
    //! @param[in ]: electronPressure is the pressure of electron.
    //! @param[out]: return the non-dimensional value of translation-rotation temperature.
    RDouble ComputeTranslationAndRotationTemperatureByPrimitiveVariables(RDouble *primitiveVariables, RDouble electronPressure = 0.0);

    //! To compute translation-rotation temperature.
    //! @param[in ]: transRotationEnergy is the total internal energy of translation and rotation mode.
    //! @param[in ]: massFraction is an array of saving the mass fractions of species.
    //! @param[out]: return the non-dimensional value of translation-rotation temperature.
    RDouble ComputeTranslationAndRotationTemperature(RDouble *massFraction, RDouble transRotationEnergy);

    //! To compute vibration temperature.
    //! @param[in ]: massFraction is an array of saving the mass fractions of species.
    //! @param[in ]: vibrationEnergy is the total internal energy of vibration mode.
    //! @param[in ]: initTemperature is the initial temperature of the newton iteration method.
    //! @param[out]: return the non-dimensional value of vibration temperature.
    RDouble ComputeVibrationTemperature(RDouble *massFraction, RDouble vibrationEnergy, RDouble initTemperature);

    //! To compute electron temperature.
    //! @param[in ]: massFraction is an array of saving the mass fractions of species.
    //! @param[in ]: electronEnergy is the total internal energy of electron mode.
    //! @param[in ]: initTemperature is the initial temperature of the newton iteration method.
    //! @param[out]: return the non-dimensional value of electron temperature.
    RDouble ComputeElectronTemperature(RDouble *massFraction, RDouble electronEnergy, RDouble initTemperature);

    //! To compute vibration-electron temperature.
    //! @param[in ]: massFraction is an array of saving the mass fractions of species.
    //! @param[in ]: vibrationElectronEnergy is the total internal energy of vibration and electron mode.
    //! @param[in ]: initTemperature is the initial temperature of the newton iteration method.
    //! @param[out]: return the non-dimensional value of vibration-electron temperature.
    RDouble ComputeVibrationAndElectronTemperature(RDouble *massFraction, RDouble vibrationElectronEnergy, RDouble initTemperature);

    //! To compute the temperature of one-temperature model.
    //! @param[in ]: massFraction is an array of mass fractions of species.
    //! @param[in ]: internalEnergy is total internal energy of the mixed gas.
    //! @param[out]: return the non-dimensional value of temperature.
    RDouble ComputeOneTemperatureModelTemperature(RDouble *massFraction, RDouble internalEnergy);

    //! To compute non-equilibrium temperature of chemical reaction flow for one-temperature model using the curve fitting method.
    //! @param[in ]: massFraction is an array of mass fractions of species.
    //! @param[in ]: internalEnergy is total internal energy of the mixed gas.
    //! @param[out]: return the non-dimensional value of temperature.
    RDouble ComputeNonequilibrumTemperatureViaBisectionMethod(RDouble *massFraction, RDouble internalEnergy);

    //! To compute non-equilibrium temperature of chemical reaction flow for one-temperature model using the newton method.
    //! @param[in ]: massFraction is an array of mass fractions of species.
    //! @param[in ]: internalEnergy is total internal energy of the mixed gas.
    //! @param[in ]: initTemperature is the initial temperature of the newton iteration method.
    //! @param[out]: return the non-dimensional value of temperature.
    RDouble ComputeNonequilibrumTemperature(RDouble *massFraction, RDouble internalEnergy, RDouble initTemperature);

    //! To obtain the electron pressure.
    //! @param[in ]: primitiveVariables is an array of saving the primary variables.
    //! @param[in ]: electronTemperature is non-dimensional value of electron temperature.
    //! @param[out]: return the non-dimensional value of electron pressure.
    RDouble GetElectronPressure(RDouble *primitiveVariables, RDouble electronTemperature);

    //! To compute the heat conductivity of translation and rotation using the Eucken formula.
    //! @param[in ]: viscosity is array of saving the species viscosity.
    //! @param[out]: conductivity is an array of saving the translation-rotation heat conductivity of each species.
    void ComputeSpeciesTranslationAndRotationHeatConductivity(RDouble *viscosity, RDouble *conductivity);

    //! To compute the heat conductivity of vibration using the Eucken formula.
    //! @param[in ]: viscosity is array of saving the species viscosity.
    //! @param[in ]: speciesVibrationCv is the array of saving the vibration specific heat at constant pressure.
    //! @param[out]: conductivity is an array of saving the vibration heat conductivity of each species.
    void ComputeSpeciesVibrationHeatConductivity(RDouble *viscosity, RDouble *speciesVibrationCv, RDouble *conductivity);

    //! To compute the heat conductivity of electron using the Eucken formula.
    //! @param[in ]: viscosity is array of saving the species viscosity.
    //! @param[in ]: speciesElectronCv is the array of saving the electron specific heat at constant pressure.
    //! @param[out]: conductivity is an array of saving the electron heat conductivity of each species.
    void ComputeSpeciesElectronHeatConductivity(RDouble *viscosity, RDouble *speciesElectronCv, RDouble *conductivity);

    //! To compute the translation-rotation heat conductivity of the mixed gas using Wilke method.
    //! @param[in ]: primitiveVariables is an array of saving the primary variables.
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronTemperature denotes the non-dimensional value of electron temperature.
    //! @param[out]: return the translation-rotation heat conductivity of the mixed gas.
    RDouble GetMixedGasTranslationAndRotationHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble electronTemperature);

    //! To compute the vibration heat conductivity of the mixed gas using Wilke method.
    //! @param[in ]: primitiveVariables is an array of saving the primary variables.
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: vibrationTemperature denotes the non-dimensional value of vibration temperature.
    //! @param[in ]: electronTemperature denotes the non-dimensional value of electron temperature.
    //! @param[out]: return the vibration heat conductivity of the mixed gas.
    RDouble GetMixedGasVibrationHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature);

    //! To compute the electron heat conductivity of the mixed gas using Wilke method.
    //! @param[in ]: primitiveVariables is an array of saving the primary variables.
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronTemperature denotes the non-dimensional value of electron temperature.
    //! @param[out]: return the electron heat conductivity of the mixed gas.
    RDouble GetMixedGasElectronHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble electronTemperature);

    void ComputeMixedGasTranslationVibrationElectronHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature, RDouble &transRotationConductivity, RDouble &vibrationConductivity, RDouble &electronConductivity);

    //! To obtain the heat conductivity of the mixed gas using Wilke method, the value is identified by the type of the heat conductivity.
    //! @param[in ]: moleFraction is an array of saving the mole fractions of species.
    //! @param[in ]: conductivity is an array of saving the heat conductivity of species. The heat conductivity could be translation-rotation\n
    //! heat conductivity, vibration heat conductivity or electron heat conductivity.
    //! @param[in ]: Phi denotes the array of the partion function.
    //! @param[out]: return the heat conductivity of the mixed gas.
    RDouble GetMixedGasHeatConductivityWithWilkeFormula(RDouble *moleFraction, RDouble *conductivity, RDouble *phi);

    //! To compute the non-dimensional values of viscosities of species using Gupta or Blotter method. 
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronTemperature denotes the non-dimensional value of electron temperature.
    //! @param[out]: speciesViscosity is an array to storing the non-dimensional values of viscosities of species.
    void ComputeSpeciesViscosityWithCurveFitMethod(RDouble transRotationTemperature, RDouble electronTemperature, RDouble *speciesViscosity);

    //! To compute the dimensional values of viscosities of species using Gupta or Blotter method. 
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronTemperature denotes the non-dimensional value of electron temperature.
    //! @param[out]: speciesViscosity is an array to storing the dimensional values of viscosities of species(Unit: kg/(m*s)).
    void ComputeSpeciesDimensionalViscosityWithCurveFitMethod(RDouble transRotationTemperature, RDouble electronTemperature, RDouble *speciesViscosity);

    //! To compute the non-dimensional values of viscosities of species using Lennard-Jones method. 
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronPressure denotes the non-dimensional value of electron pressure.
    //! @param[out]: speciesViscosity is an array to storing the non-dimensional values of viscosities of species.
    void ComputeSpeciesViscosityWithLennardJonesMethod(RDouble transRotationTemperature, RDouble electronPressure, RDouble *speciesViscosity);

    //! To compute the dimensional values of viscosities of species using using Lennard-Jones method. 
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronPressure denotes the non-dimensional value of electron pressure.
    //! @param[out]: speciesViscosity is an array to storing the dimensional values of viscosities of species(Unit: kg/(m*s)).
    void ComputeSpeciesDimensionalViscosityWithLennardJonesMethod(RDouble transRotationTemperature, RDouble electronPressure, RDouble *speciesViscosity);

    //! To compute the collision area between the same species using the curve fit method.
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronPressure denotes the non-dimensional value of electron pressure.
    //! @param[out]: piOmega is an array to storing the dimensional values of collision area of species(Unit: m2).
    void ComputeAverageCollisionAreaOmega22(RDouble transRotationTemperature, RDouble electronPressure, RDouble *piOmega);

    //! To obtain the viscosity of the mixed gas using Wilke method.
    //! @param[in ]: moleFraction is an array of saving the mole fractions of species.
    //! @param[in ]: viscosity is an array of saving the viscosity of species.
    //! @param[in ]: phi denotes the array of the partion function.
    //! @param[out]: return the viscosity of the mixed gas.
    RDouble GetMixedGasViscosityWithWilkeFormula(RDouble *moleFraction, RDouble *viscosity, RDouble *phi);

    //! To obtain the diffusion coefficient of each species. added by LiPeng in Mar. 11,2019.
    //! @param[in ]: massFraction is an array of saving the primary variables.
    //! @param[in ]: viscosityLaminar is the laminar viscosity, and viscosityTurbulence is the turbulent viscosity.
    //! @param[out]: speciesDiffusionCoef stores the species mass diffusion coefficient of each component.
    //! The function returns the diffusion coefficients of each species which saved in rDs = (1 - Cs)/(1- Xs) * (mul/Scl + mut/Sct).
    void ComputeSpeciesMassDiffusionCoefficient(RDouble *massFraction, RDouble viscosityLaminar, RDouble viscosityTurbulence, RDouble *speciesDiffusionCoef);

    //! To obtain the mass diffusion coefficient with CLN method, the obtained value is enjoyed by all of the species.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronTemperature denotes the non-dimensional value of electron temperature.
    //! @param[out]: return the mass diffusion coefficient.
    RDouble GetMassDiffusionCoefficientWithCLNModel(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble electronTemperature);

    //! Compute the coefficient called gama system in the multi-temperature model.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: electronTemperature is the non-dimensional value of vibration or electron temperature.
    //! @param[out]: gama are partial derivative of electron pressure to vibration or electron energy.
    //! @param[out]: speciesGama is an array of partial derivatives of electron pressure to species density.
    void GetElectronPressurePartialDerivatives(RDouble *primitiveVars, RDouble electronTemperature, RDouble &gama, RDouble *speciesGama);

    //! Compute the coefficient called alpha system in the multi-temperature model.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: transRotationTemperature is the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronTemperature is the non-dimensional value of vibration or electron temperature.
    //! @param[out]: alpha is partial derivative of pressure to internal energy.
    //! @param[out]: beta is partial derivative of pressure to total density.
    //! @param[out]: speciesBeta is an array of partial derivatives of pressure to species density..
    void GetPressurePartialDerivatives(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble electronTemperature,
                                       RDouble &alpha, RDouble &beta, RDouble *speciesBeta);

    //! Compute the main partial derivatives in the multi-temperature model.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: transRotationTemperature is the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronTemperature is the non-dimensional value of vibration or electron temperature.
    //! @param[out]: alpha is the partial derivative of pressure to total energy.
    //! @param[out]: beta is the partial derivative of pressure to total density.
    //! @param[out]: gama is the partial derivative of electron pressure to vibration-electron enegy.
    //! @param[out]: speciesBeta,speciesGama are the partial derivatives of pressure and electron pressure to species density.
    void GetMultiTemperatureModelPartialDerivatives(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble electronTemperature,
                                                    RDouble &alpha, RDouble &beta, RDouble &gama, RDouble *speciesBeta, RDouble *speciesGama);
    void GetMultiTemperatureModelPartialDerivativesR(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble electronTemperature,
                                                    RDouble &alpha, RDouble &beta, RDouble &gama, RDouble *speciesBeta, RDouble *speciesGama,
                                                    RDouble squareVelocity, RDouble *speciesEtrs, RDouble totalCvtr, RDouble totalCvv, RDouble totalCve);

    //! Compute partial derivative of pressure to total density.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[out]: return the partial derivative called alpha.
    RDouble ComputePressureDerivativeToTotalDensity(RDouble *primitiveVars);

    //! To compute the static enthalpy of the mixed gas.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: transRotationTemperature denotes the temperature of translation-rotation mode.
    //! @param[in ]: vibrationTemperature denotes the temperature of vibration mode.
    //! @param[in ]: electronTemperature denotes the temperature of electron mode.
    RDouble GetMixedGasEnthalpy(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature);

    RDouble ComputeMolecularWeightReciprocal(RDouble *prim);

    RDouble ComputeMolecularWeight(RDouble *prim);

    //! To compute heat conductivity of mixture gas using Eucken formula and Wassilewa relation.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: transRotationTemperature is the non-dimensional value of translation-rotation temperature of the mixture gas.
    //! @param[in ]: vibrationTemperature is the non-dimensional value of vibration temperature of the mixture gas.
    //! @param[in ]: electronTemperature is the non-dimensional value of electron temperature of the mixture gas.
    //! @param[out]: transRotationConductivity is heat conductivity of translation-rotation of mixied gas.
    //! @param[out]: vibrationConductivity is heat conductivity of vibration of mixied gas.
    //! @param[out]: electronConductivity is heat conductivity of electron of mixied gas.
    void ComputeMixtureGasHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature,
        RDouble &transRotationConductivity, RDouble &vibrationConductivity, RDouble &electronConductivity);

    void ComputeMixtureGasSpecificHeatAtConstantVolume(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature,
        RDouble &cvtr, RDouble &cvv, RDouble &cve);

    //! Compute the partial derivatives of source terms of chemical reaction equations.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[in ]: chemicalSourceTerms stores the source terms of chemical reaction equations in the specified control cells.
    //! @param[in ]: chemicalSourceDerivatives stores the partial derivatives of the source terms of chemical reaction equations in the specified control cells.
    //! @param[out]: energySourceDerivatives returns the partial derivatives of the source terms of energy equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    void ComputeEnergySourceDerivatives(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
        RDouble ***chemicalSourceDerivatives, RDouble ***energySourceDerivatives, int nCellNumber);

    //! Compute the partial derivatives of source terms of chemical reaction equations.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[in ]: chemicalSourceTerms stores the source terms of chemical reaction equations in the specified control cells.
    //! @param[in ]: chemicalSourceDerivatives stores the partial derivatives of the source terms of chemical reaction equations in the specified control cells.
    //! @param[out]: energySourceDerivatives returns the diagonal element of partial derivatives to the source terms of energy equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    void ComputeEnergySourceDiagonalDerivatives(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
        RDouble ***chemicalSourceDerivatives, RDouble **energySourceDerivatives, int nCellNumber);

    //! To obtain the maximum mass diffusion coefficient among the species. added by LiPeng in Jan. 2019.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: viscosityLaminar is the laminar viscosity, and viscosityTurbulence is the turbulent viscosity.
    //! The function returns the maximum mass diffusion coefficient called (1 - Cs)/(1- Xs) * (mul/Scl + mut/Sct).
    RDouble GetMaximumSpeciesMassDiffusionCoef(RDouble *primitiveVars, RDouble viscosityLaminar, RDouble viscosityTurbulence);

    RDouble ComputeMolecularWeightReciprocalWithoutElectron(RDouble *prim, RDouble &ce_div_me);

    //! To obtain the diffusion coefficient of each species. added by LiPeng in Mar. 11,2019.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: viscosityLaminar is the laminar viscosity, and viscosityTurbulence is the turbulent viscosity.
    //! @param[out]: speciesDiffusionCoef stores the species mass diffusion coefficient of each component.
    //! The function returns the diffusion coefficients of each species which saved in rDs = (1 - Cs)/(1- Xs) * (mul/Scl + mut/Sct).
    void GetSpeciesMassDiffusionCoef(RDouble *primitiveVars, RDouble viscosityLaminar, RDouble viscosityTurbulence, RDouble *speciesDiffusionCoef);

    //! To obtain the non-dimensional constant pressure specific heat in accordance with the primary variables.
    //! prim is an array of saving the primary variables, the non-dimensional temperature called Tm and the non-dimensional \n
    //! Constant pressure specific heat called cp are returned.
    void GetTemperatureAndCP(RDouble *prim, RDouble &Tm, RDouble &cp);
    
    //! To obtain the non-dimensional enthalpy of each species, the variable Tm is the non-dimensional temperature of the mixture gas,\n
    //! The enthalpy of each species saved in the array of hs is returned.
    void GetEverySpeciesEnthalpy(RDouble trTemperature, RDouble vTemperature, RDouble eTemperature, RDouble *hs);

    //! To compute the electron energy of the mixed gas.
    //! @param[in ]: prim is an array of saving the primary variables.
    //! @param[in ]: Te is the non-dimensional electron temperature of the mixture gas.
    //! @param[out]: return the total electron energy.
    RDouble GetMixedGasElectronEnergy(RDouble *prim, RDouble Te);

    //! To compute the vibration specific heat at constant volume of the mixed gas.
    //! @param[in ]: prim is an array of saving the primary variables.
    //! @param[in ]: Tv is the non-dimensional vibration temperature of the mixture gas.
    //! @param[out]: return the total vibration specific heat at constant volume.
    RDouble GetMixedGasVibrationSpecificHeatAtConstantVolume(RDouble *prim, RDouble Tv);

    //! To compute the vibration specific heat at constant volume of the mixed gas.
    //! @param[in ]: prim is an array of saving the primary variables.
    //! @param[in ]: Te is the non-dimensional electron temperature of the mixture gas.
    //! @param[out]: return the dimensional value of total vibration specific heat at constant volume.
    RDouble GetMixedGasElectronSpecificHeatAtConstantVolume(RDouble *prim, RDouble Te);

    //! To compute the translation-rotation specific heat at constant volume of the mixed gas.
    //! @param[in ]: prim is an array of saving the primary variables. 
    //! @param[out]: return the dimensional value of total translation-rotation specific heat at constant volume.
    RDouble GetMixedGasTransAndRotatSpecHeatAtConstVolume(RDouble *prim);

    //! To compute the translation-rotation energy of the mixed gas.
    //! @param[in ]: prim is an array of saving the primary variables.
    //! @param[in ]: Ttr is the non-dimensional translation-rotation temperature of the mixture gas. 
    //! @param[out]: return the total translation-rotation energy.
    RDouble GetMixedGasTranslationAndRotationEnergy(RDouble *prim, RDouble Ttr);

    //! To compute the vibration energy of the mixed gas.
    //! @param[in ]: prim is an array of saving the primary variables.
    //! @param[in ]: Tv is the non-dimensional vibration temperature of the mixture gas.
    //! @param[out]: return the total vibration energy.
    RDouble GetMixedGasVibrationEnergy(RDouble *prim, RDouble Tv);

    //! Compute the frozen sound speed in the multi-temperature model.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[out]: return the frozen sound speed.
    RDouble ComputeFrozenSoundSpeed(RDouble *primitiveVars);

    //! To obtain the non-dimensional temperature according to the primary variables,\n
    //! prim is an array of saving the primary variables, the non-dimensional temperature is returned.
    RDouble GetTemperature(RDouble *prim);

public:

    //! Compute the specific heat of species.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: temperatures indicates the array of temperatures including translational-rotational temperature, vibrational temperature and electron temperature.
    //! @param[out]: spciesCvtr denotes the array of translational-rotational specific heat for species.
    //! @param[out]: spciesCvv  denotes the array of vibrational specific heat for species.
    //! @param[out]: totalCvtr and totalCvv denotes the total specific heat of mixture gas.
    void ComputeSpecificHeat(RDouble *primitiveVars, RDouble *temperatures, RDouble *spciesCvtr, RDouble *spciesCvv, RDouble &totalCvtr, RDouble &totalCvv);

    //! Compute the energies of species.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: temperatures indicates the array of temperatures including translational-rotational temperature, vibrational temperature and electron temperature.
    //! @param[out]: spciesEtr denotes the array of translational-rotational energy for species.
    //! @param[out]: spciesEv  denotes the array of vibrational energy for species.
    //! @param[out]: spciesEnthalpy denotes the array of enthalpies for species.
    //! @param[out]: totalEtr and totalEv denotes the total energy of mixture gas.
    //! @param[out]: totalEnthalpy denotes the total enthalpy of mixture gas.
    void ComputeEnergy(RDouble *primitiveVars, RDouble *temperatures, RDouble *spciesEtr, RDouble *spciesEv, RDouble *spciesEnthalpy, RDouble &totalEtr, RDouble &totalEv, RDouble &totalEnthalpy);

    //! Compute the partial derivatives of pressure.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: temperatures indicates the array of temperatures including translational-rotational temperature, vibrational temperature and electron temperature.
    //! @param[out]: alpha indicates the partial derivative of pressure vs. total energy(dp/d(rho*E)).
    //! @param[out]: beta indicates the partial derivative of pressure vs. density(dp/d(rho)).
    //! @param[out]: eta indicates the partial derivative of pressure vs. vibration energy(dp/d(rho * Ev)).
    //! @param[out]: speciesBeta indicates the partial derivatives of pressure vs. species density(dp/d(rho * ci)).
    void ComputePressureDerivatives(RDouble *primitiveVars, RDouble *temperatures, RDouble &alpha, RDouble &beta, RDouble &eta, RDouble *speciesBeta);

    void ComputeSpecificHeatByFitMethod(RDouble temperature, RDouble *cvSpecies, bool isVibrationModel = true);
    void ComputeEnthalpyByFitMethod(RDouble temperature, RDouble *hSpecies, bool isVibrationModel = true);

    void ComputeSpecificHeatByPolynomialFitMethod(RDouble temperature, RDouble *cvSpecies, bool isVibrationModel = true);
    void ComputeEnthalpyByPolynomialFitMethod(RDouble temperature, RDouble *hSpecies, bool isVibrationModel = true);
    void ComputeTemperatureFromTotalTemperature(RDouble *massF, RDouble totalTemperature, RDouble mach, RDouble &Temperature, RDouble &gama, int Dimensional = 0);

};

void ChemCompressData(DataContainer *&cdata);

void ChemDecompressData(DataContainer *cdata);

//! The function computes the formula (M +- R) * dQ, where M and R denote the Jacobian matrix and spectrum radius of inviscid flux,respectively.
//! @param[in ]: primitiveVars is the array that stores the primitive variables of the current cell.
//! @param[in ]: temperature is the array that stores the temperature values of the current cell.
//! @param[in ]: nx, ny and nz are components of the average normal vector located in the center of the current cell.
//! @param[in ]: faceArea denotes the avarage area located in the center of the current cell.
//! @param[in ]: moveVelocity denotes the movement velocity located in the center of the current cell, and it is used for dynamic mesh.
//! @param[in ]: deltaQ is the increment of conservative vector.
//! @param[out]: flux is the inviscid flux that equals to the Jacobian matrix multiplies the difference values of the conservative variables.
//! @param[in ]: radius denotes the spectrum radius of invsicid Jacobian matrix.
//! @param[in ]: iSign denotes the computation symbol, -1 denotes M substracts R, 1 denotes M plus R.
void ChemicalMXDQ(RDouble *primitiveVars, RDouble *temperature, RDouble nx, RDouble ny, RDouble nz, RDouble faceArea, RDouble moveVelocity,
                  RDouble *deltaQ, RDouble *flux, int nNSEquation, int nLaminar, int nTemperatureModel, RDouble radius, int iSign,int nElectronIndex);
void ChemicalMXDQR(RDouble *primitiveVars, RDouble *temperature, RDouble nx, RDouble ny, RDouble nz, RDouble faceArea, RDouble moveVelocity,
                   RDouble *deltaQ, RDouble *flux, int nNSEquation, int nLaminar, int nTemperatureModel, RDouble radius, int iSign, int nElectronIndex,
                   RDouble *speciesCvs, RDouble *speciesEtrs, RDouble totalH,  RDouble totalCv, RDouble totalCvtr, RDouble totalCvv, RDouble totalCve);
}

}
