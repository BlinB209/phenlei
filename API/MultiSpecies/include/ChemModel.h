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
//! @file      ChemModel.h
//! @brief     Define chemical reactions.
//! @author    Meng Liyuan, Li Peng.

#pragma once
#include <string>
#include "Precision.h"
#include "TypeDefine.h"
using namespace std;

namespace PHSPACE
{
namespace GAS_SPACE
{
//! Define the class of ChemicalReactions that is used to construct the chemical reactions and their reaction rates.
class ChemicalReactions{
private:
    int nReactionType; 
    //nType = 0 indicates the Gupta Model, nType = 1 indicates the Park Model, and nType = 2 indicates the Dunn-Kang Model.
    //nType = 3 indicates the Park Model of Mars, and nType = 4 indicates the McKenzie Model of Mars.
    int nSpecies, nReaction;
    string *speciesName; //The names of the species in the mixture.
private:
    int **stoichForward, **stoichBackward;
    RDouble **collisionCoef;
    RDouble **reactionData; //curve fitting data of reaction rate.

    //The type of chemical reaction, reactionType[i][j], i indicates the number of the reaction, j=0 or j=1.
    //reactionType[i][0]=0 for reaction without the third collision bodies.
    //reactionType[i][0]=1 for collision reaction. e.g. AB + M = A + B + M.
    //reactionType[i][1]=0 for dissociation reaction(AB+M=A+B+M)
    //reactionType[i][1]=1 for dissociation reaction with electron collision(AB+e=A+B+e).
    //reactionType[i][1]=2 for neutral exchange reaction(AB+C=AC+B)
    //reactionType[i][1]=3 for charge exchange reaction(AB+ + C=A+ + BC).
    //reactionType[i][1]=4 for recombination ionization reaction(A+B=AB+ + e)
    //reactionType[i][1]=5 for ionization reaction with electron collision(A+e=A+ + e + e).
    int **reactionType;

public:
    ChemicalReactions();
    ChemicalReactions(int ns, int nr);
    ~ChemicalReactions();
    void AllocateMemories();
    void DeAllocateMemories();
    int GetChemicalReactionType() const { return nReactionType; }
    void SetChemicalReactionType(int ntype) { this->nReactionType = ntype; }
    void GetChemicalModelInfo(int &ns, int &nr) { ns = nSpecies; nr = nReaction; }
    void SetChemicalModelInfo(int ns, int nr);
    string * GetSpeciesNames() const { return this->speciesName; }
    int ** GetForwardStoichiometricCoefficient() const { return stoichForward; }
    int ** GetBackwardStoichiometricCoefficient() const { return stoichBackward; }
    RDouble ** GetThirdBodyCollisionCoefficient() const { return collisionCoef; }
    RDouble ** GetReactionRateCurveFittingData() const { return reactionData; }
    int ** GetReactionType() const { return reactionType; }
public:
    //! To obtain the data of the Gupta 11s20r reaction model.
    //! Param[in/out]: reactionCoef denotes the stoichiometric coefficients of forward and backward reaction equation.
    //! Param[in/out]: collisionCoef indicates the collision coefficients of third body in the reaction equation.
    //! Param[in/out]: rateCoef indicates the curve fitting data of reaction rates.
    void GenerateGupta11S20RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

    //! To obtain the data of the Park 11s48r reaction model.
    //! Param[in/out]: reactionCoef denotes the stoichiometric coefficients of forward and backward reaction equation.
    //! Param[in/out]: collisionCoef indicates the collision coefficients of third body in the reaction equation.
    //! Param[in/out]: rateCoef indicates the curve fitting data of reaction rates.
    void GeneratePark11S48RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

    //! To obtain the data of the Dunn-Kang 11s26r reaction model.
    //! Param[in/out]: reactionCoef denotes the stoichiometric coefficients of forward and backward reaction equation.
    //! Param[in/out]: collisionCoef indicates the collision coefficients of third body in the reaction equation.
    //! Param[in/out]: rateCoef indicates the curve fitting data of reaction rates.
    void GenerateDunnKang11S26RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

    //! To obtain the data of the Park 12s12r reaction model of Mars.
    //! Param[in/out]: reactionCoef denotes the stoichiometric coefficients of forward and backward reaction equation.
    //! Param[in/out]: collisionCoef indicates the collision coefficients of third body in the reaction equation.
    //! Param[in/out]: rateCoef indicates the curve fitting data of reaction rates.
    void GenerateMarsPark8S12RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

    //! To obtain the data of the McKenzie 12s12r reaction model of Mars.
    //! Param[in/out]: reactionCoef denotes the stoichiometric coefficients of forward and backward reaction equation.
    //! Param[in/out]: collisionCoef indicates the collision coefficients of third body in the reaction equation.
    //! Param[in/out]: rateCoef indicates the curve fitting data of reaction rates.
    void GenerateMarsMcKenzie8S12RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

    void GenerateGuptaAblation16S36RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

    void GenerateParkAblation16S38RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

    //void GenerateDunnKangAblation17S44RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

    void GenerateDunnKangAblation16S42RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

    void GenerateCombustionGas12S20RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);
    
    void GenerateGasMixture2S1RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

    //! Initialization for chemical reaction model.
    void Initialization();

    //! Import data of basic chemical reaction model.
    //! Param[in ]: nType indicates the reaction type.
    void CreateBasicChemicalReactionModel(int nType = 0);

    //! Create a new chemical reaction model according to the name list.
    //! Param[in ]: namelist indicates the vector of the species names in the new model.
    //! Param[out]: pNewModel indicates the new chemical reaction model.
    void CreateChemicalReactionModel(vector<string> &namelist, ChemicalReactions *pNewModel);
};

//Find the species in the name list, if found, return the index of the species in the name list.
int SearchSpeciesInNameList(string strName, vector<string> &namelist);

//! Define the class of ChemkinFitData that is used to compute the specific heat of species via curve fitting method.
class ChemkinFitData{
private:
    int nSpecies;//number of species;
    string *speciesName; //The names of the species in the mixture.
    RDouble ***curveFittingData;
    RDouble **enthalpyFitData; //curve fit coefficients for enthalpy.
    RDouble ***curveFittingDataT2; //Curve fit coefficients for two-temperature model.
    int nInterval, nCoefNum;
    RDouble *temperatureRange;
    RDouble *temperatureStart; //The start temperature when the curve fit method is valid.
public:
    ChemkinFitData();
    ~ChemkinFitData();
    void AllocateMemories();
    void DeAllocateMemories();
    void GetInterpolationInfo(int &ns, int &m, int &n) { ns = nSpecies; m = nInterval; n = nCoefNum; }
    void SetInterpolationInfo(int ns, int m = 5, int n = 6);
    string * GetSpeciesNames() const { return this->speciesName; }
    RDouble *** GetCurveFittingData() const { return curveFittingData; }
    RDouble * GetTemperatureRange() const { return temperatureRange; }
    RDouble ** GetEnthalpyFitData() const { return enthalpyFitData; }
    RDouble * GetStartTemperature() const { return temperatureStart; }
    RDouble *** GetCurveFittingDataT2() const { return curveFittingDataT2; }
public:
    void SetCurveFittingDataT2();

    //! To Generate the curve fitting data in the chemkin method.
    void GenerateFullChemkinData();

    //! Create a new ChemkinFitData according to the name list.
    //! Param[in ]: namelist indicates the vector of the species names in the new model.
    //! Param[out]: pNewData indicates the new ChemkinFitData.
    void CreateChemkinFitData(vector<string> &namelist, ChemkinFitData *pNewData);
};

//! Define the class of TransportFitData that is used to compute the transport coefficients such as viscosity.
class TransportFitData{
private:
    int nSpecies; //number of species;
    string *speciesName; //The names of the species in the mixture.
    RDouble **viscosityFitData; //the curve fitting data of viscosity.
    RDouble **collisionAreaFitData; // the curve fitting data of the average collision area.
public:
    TransportFitData();
    ~TransportFitData();
    void AllocateMemories();
    void DeAllocateMemories();
    void SetSpeciesNumber(int ns);
    int GetSpeciesNumber() const { return nSpecies; }
    string * GetSpeciesNames() const { return this->speciesName; }
    RDouble ** GetViscosityCurveFittingData() const { return viscosityFitData; }
    RDouble ** GetCollisionAreaCurveFittingData() const { return collisionAreaFitData; }
public:
    //! To Generate the full curve fitting data of 16 species.
    void GenerateFullCurveFittingData();

    //! Create a new TransportFitData according to the name list.
    //! Param[in ]: namelist indicates the vector of the species names in the new model.
    //! Param[out]: pNewData indicates the new TransportFitData.
    void CreateTransportFitData(vector<string> &namelist, TransportFitData *pNewData);
};

struct VibrationModeData{
    int nMode; //The number of vibrational mode.
    int g[10]; //The degeneracy of each vibrational mode.
    RDouble Tv[10]; //The characteristic vibrational temperature of each vibrational mode, (K).
};

//! Define the class of AtmosphereData that is used to save the physical and chemical data.
class AtmosphereData{
private:
    int nSpecies; //number of species;
    string *speciesName; // the name of species.
    int *nSpeciesType; // the type of species.
    int *nSpeciesCharge; // the charge of species.
    int *nSpeciesG0; // the degeneracy of species.
    int *nSpeciesG1; // the degeneracy of species.
    RDouble *speciesMass; // the molecular weight of species, (g/mol).
    RDouble *speciesTtr; // the characteristic temperature of species, (K).
    VibrationModeData *speciesVibrationalModeData; //The data of vibrational mode.
    RDouble *speciesTes; // the characteristic electronic temperature of species, (K).
    RDouble *speciesEnthalpy; // the formation enthalpy of species, (J/kg).
    RDouble *speciesSigma; // the collision area of species, (m2).
public:
    AtmosphereData();
    ~AtmosphereData();
    void AllocateMemories();
    void DeAllocateMemories();
    void SetSpeciesNumber(int ns);
    string * GetSpeciesNames() const { return speciesName; }
    int GetSpeciesNumber() const { return nSpecies; }
    int * GetSpeciesType() const { return nSpeciesType; }
    int * GetSpeciesCharge() const { return nSpeciesCharge; }
    RDouble * GetSpeciesMolecularWeight() const{ return speciesMass; }
    int * GetSpeciesDegeneracy0() const { return nSpeciesG0; }
    int * GetSpeciesDegeneracy1() const { return nSpeciesG1; }
    RDouble * GetSpeciesTemperature() const { return speciesTtr; }
    VibrationModeData * GetSpeciesVibrationalModeData() const { return speciesVibrationalModeData; }
    RDouble * GetSpeciesElectronTemperature() const { return speciesTes; }
    RDouble * GetSpeciesFormationEnthalpy() const { return speciesEnthalpy; }
    RDouble * GetSpeciesCollisionArea() const { return speciesSigma; }
public:
    //! To Generate the physical and chemicla data of 16 species.
    void GenerateFullPhysicalChemicalData();

    //! Create a new AtmosphereData according to the name list.
    //! Param[in ]: namelist indicates the vector of the species names in the new model.
    //! Param[out]: pNewData indicates the new AtmosphereData.
    void CreateAtmosphereData(vector<string> &namelist, AtmosphereData *pNewData);
};

}// end of GAS_SPACE
}// end of PHSPACE
