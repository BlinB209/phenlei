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
//! @file      Force.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include "DataContainer.h"
using namespace std;

namespace PHSPACE
{

class ActionKey;
class DataContainer;
class AerodynamicForceManager;
class AerodynamicForce;

class AerodynamicForceBody
{
public:
    AerodynamicForceBody();
    ~AerodynamicForceBody();
protected:
    static vector< std::set< string > > bodies;
    static bool initialize;
public:
    static void Initialize();
    static bool FaceInAerodynamicForceBodies(string bodyName, int aerodynamicForceBodyIndex);
};

bool FaceInAerodynamicForceBodies(string bodyName, int aerodynamicForceBodyIndex);

class BasicAerodynamicForce
{
protected:
    RDouble forceX, forceY, forceZ;
    RDouble forceXByPressure, forceYByPressure, forceZByPressure;
    RDouble forceXByViscosity, forceYByViscosity, forceZByViscosity;
protected:
    RDouble momentX, momentY, momentZ;
    RDouble momentXByPressure, momentYByPressure, momentZByPressure;
    RDouble momentXByViscosity, momentYByViscosity, momentZByViscosity;
protected:
    RDouble momentReferenceX, momentReferenceY, momentReferenceZ;
public:
    BasicAerodynamicForce();
    ~BasicAerodynamicForce();
protected:
    void ZeroBasicAerodynamicForce();
public:
    void SetMomentReferencePoint(RDouble momentReferenceX, RDouble momentReferenceY, RDouble momentReferenceZ);
    void SetPressureAerodynamicForce(RDouble forceXByPressure, RDouble forceYByPressure, RDouble forceZByPressure);
    void SetViscousAerodynamicForce(RDouble forceXByViscosity, RDouble forceYByViscosity, RDouble forceZByViscosity);
    void ComputeResultantAerodynamicForce();
    void SetPressureAerodynamicMoment(RDouble faceCenterX, RDouble faceCenterY, RDouble faceCenterZ);
    void SetViscousAerodynamicMoment(RDouble faceCenterX, RDouble faceCenterY, RDouble faceCenterZ);
    void ComputeAerodynamicMoment(RDouble faceCenterX, RDouble faceCenterY, RDouble faceCenterZ);
    void ComputeResultantAerodynamicMoment();
public:
    void WriteAerodynamicForceToDataContainer(DataContainer *dataContainer);
    void ReadAerodynamicForceFromDataContainer(DataContainer *dataContainer);
    void AddAerodynamicForce(BasicAerodynamicForce *basicAerodynamicForce);
    //void StatisticalForce(BasicAerodynamicForce *basicAerodynamicForce, int numberOfAverage);
};

class BasicAleForce
{
protected:
    RDouble forceX, forceY, forceZ;
    RDouble momentX, momentY, momentZ;
    RDouble powerConsumption;
public:
    BasicAleForce();
    ~BasicAleForce();
public:
    void ZeroBasicAleForce();
    void SetBasicAleForceByAerodynamicForce(AerodynamicForce *aerodynamicForce);
public:
    RDouble GetForceX() { return forceX; }
    RDouble GetForceY() { return forceY; }
    RDouble GetForceZ() { return forceZ; }
    RDouble GetMomentX() { return momentX; }
    RDouble GetMomentY() { return momentY; }
    RDouble GetMomentZ() { return momentZ; }
    RDouble GetPowerConsumption() { return powerConsumption; }
public:
    void ReadBasicAleForceFromFile(fstream &file);
    void ReadBasicAleForceFromBasicAleForce(BasicAleForce *basicAleForce);
    void WriteBasicAleForceToFile(fstream &file);
};

class AleForceOfBody
{
public:
    AleForceOfBody();
    ~AleForceOfBody();
protected:
    int bodyIndex;
    int numberOfTimeLayers;
    vector< BasicAleForce * > basicAleForceContainer;
public:
    int  GetBodyIndex() const { return bodyIndex; }
    void SetBodyIndex(int bodyIndex) { this->bodyIndex = bodyIndex; }

    int  GetNumberOfTimeLayers() const { return numberOfTimeLayers; }
    void SetNumberOfTimeLayers(int numberOfTimeLayers) { this->numberOfTimeLayers = numberOfTimeLayers; }
public:
    void AllocateMemory();
    void DeallocateMemory();
public:
    BasicAleForce *GetBasicAleForce(int timeLayerIndex);
    void DumpRestartAleForceOfBodyContainer(fstream &file);
    void InitializeAleForceOfBodyContainer();
    void ReadAleForceOfBodyContainer();
    void UpdateAleForceOfBodyOfNewTimeLayer();
    void UpdateUnsteadyAleForceTerm();
};

class AleForceManager
{
public:
    AleForceManager();
    ~AleForceManager();
protected:
    vector< AleForceOfBody * > aleForceOfBodyContainer;
    int numberOfMovingBodies;

    RDouble dimensionalForceCoefficient;
    RDouble dimensionalMomentCoefficient;
    RDouble dimensionalPowerCoefficient;
    int numberOfTimeLayers;
public:
    void Initialize();
    void Restart();
    void Run();
    void Post();

    void DumpRestartAleForceOfBodyContainer(fstream &file);

    BasicAleForce *GetBasicAleForce(int bodyIndex, int iTimeLayer);
public:
    int GetNumberOfTimeLayers() const { return numberOfTimeLayers; }
    void SetNumberOfTimeLayers(int numberOfTimeLayers) { this->numberOfTimeLayers = numberOfTimeLayers; }
    RDouble GetDimensionalForceCoefficient() { return dimensionalForceCoefficient; }
    RDouble GetDimensionalPowerCoefficient() { return dimensionalPowerCoefficient; }
    RDouble GetDimensionalMomentCoefficient() { return dimensionalMomentCoefficient; }
protected:
    AleForceOfBody *GetAleForceOfBody(int bodyIndex) { return aleForceOfBodyContainer[bodyIndex]; }
protected:
    void AllocateAleForceOfBodyContainer();
    void DeallocateAleForceOfBodyContainer();
    void InitializeDimensionalReferenceVariables();
    void InitializeAleForceOfBodyContainer();
    void ReadAleForceOfBodyContainer();
};

void CopyBasicAleForce(BasicAleForce *basicAleForce, BasicAleForce *basicAleForceIn);

class AerodynamicForce
{
public:
    AerodynamicForce();
    ~AerodynamicForce();
protected:
    RDouble forceX, forceY, forceZ;
    RDouble forceXByPressure, forceYByPressure, forceZByPressure;
    RDouble forceXByViscosity, forceYByViscosity, forceZByViscosity;

    RDouble momentX, momentY, momentZ;
    RDouble momentXByPressure, momentYByPressure, momentZByPressure;
    RDouble momentXByViscosity, momentYByViscosity, momentZByViscosity;
    RDouble hingeMoment;
protected:
    RDouble coefForceX, coefForceY, coefForceZ;
    RDouble coefForceXByPressure, coefForceYByPressure, coefForceZByPressure;
    RDouble coefForceXByViscosity, coefForceYByViscosity, coefForceZByViscosity;

    RDouble coefMomentX, coefMomentY, coefMomentZ;
    RDouble coefMomentXByPressure, coefMomentYByPressure, coefMomentZByPressure;
    RDouble coefMomentXByViscosity, coefMomentYByViscosity, coefMomentZByViscosity;
protected:
    RDouble totalLiftCoefficient, totalDragCoefficient, totalCrossCoefficient;
    RDouble totalPressureDragCoefficient, totalViscosityDragCoefficient;
    RDouble pressureCenter;
    RDouble DragCoefLiftCoef2piar;
    RDouble angleOfAttackDenotedByRadian, angleOfSideslipDenotedByRadian;
    RDouble coefHingeMoment;
protected:
    string bodyName;
    RDouble referenceArea, referenceLength, referenceSpainLength;
    RDouble torqueRefX, torqueRefY, torqueRefZ;
    RDouble localCoordAxis0[3];
    RDouble localCoordAxis1[3];
protected:
    void Zero();
public:
    void Initialize(RDouble angleOfAttackDenotedByRadian, RDouble angleOfSideslipDenotedByRadian, RDouble referenceAreaOfAerodynamicForce, RDouble referenceLengthOfAerodynamicMoment,
        RDouble referenceSpainLengthOfAerodynamicMoment);
    void Add(DataContainer *dataContainer);
    void DumpCoefficient(ActionKey *actionKey);
    void DumpComponentCoefficient(ActionKey *actionKey);
    void DumpCoefficient(fstream &file);
    void ComputeCoefficient();

    //! Dump global force data.
    void DumpAllGlobalForceData(const string &partName, int partID, std::ostringstream &oss);

    //! Dump local force data.
    void DumpAllLocalForceData(const string &partName, int partID, std::ostringstream &oss);

    void CollectOfAllProcessors();
public:
    RDouble GetForceCoefficientX() { return coefForceX; }
    RDouble GetForceCoefficientY() { return coefForceY; }
    RDouble GetForceCoefficientZ() { return coefForceZ; }

    RDouble GetMomentCoefficientX() { return coefMomentX; }
    RDouble GetMomentCoefficientY() { return coefMomentY; }
    RDouble GetMomentCoefficientZ() { return coefMomentZ; }
public:
    RDouble GetAerodynamicForceX() { return forceX; }
    RDouble GetAerodynamicForceY() { return forceY; }
    RDouble GetAerodynamicForceZ() { return forceZ; }

    RDouble GetAerodynamicMomentX() { return momentX; }
    RDouble GetAerodynamicMomentY() { return momentY; }
    RDouble GetAerodynamicMomentZ() { return momentZ; }

public:
    void SetReferenceArea(RDouble referenceArea) { this->referenceArea = referenceArea; }
    void SetReferenceLength(RDouble referenceLength) { this->referenceLength = referenceLength; }
    void SetReferenceSpainLength(RDouble referenceSpainLength) { this->referenceSpainLength = referenceSpainLength; }
    void SetTorqueRef(RDouble torqueRefX, RDouble torqueRefY, RDouble torqueRefZ);
    void SetHingeMomentRef(RDouble *localCoordAxis0, RDouble *localCoordAxis1);

    void SetBodyName(string bodyName) { this->bodyName = bodyName; }
    string GetBodyName() { return bodyName; }

    void GetMomentByOriginOfCoordinate(RDouble &aerodynamicMomentX0, RDouble &aerodynamicMomentY0, RDouble &aerodynamicMomentZ0);
};

void InitializeAerodynamicForce();
void AddAerodynamicForce(ActionKey *actionKey);
void CollectAerodynamicForceOfAllProcessors();
void ComputeAerodynamicForceForALE();
void DumpAerodynamicForceCoefficient(ActionKey *actionKey);

void WriteAerodynamicForceToActionKey(ActionKey *actionKey, BasicAerodynamicForce **basicAerodynamicForceContainer);
void WriteAerodynamicForceToActionKey(ActionKey *actionKey, vector< BasicAerodynamicForce * > *basicAerodynamicForceVector);

BasicAerodynamicForce **CreateBasicAerodynamicForceContainer();
void FreeBasicAerodynamicForceContainer(BasicAerodynamicForce **basicAerodynamicForceContainer);

vector< BasicAerodynamicForce * > *CreateBasicAerodynamicForceVector();
void FreeBasicAerodynamicForceVector(vector< BasicAerodynamicForce * > *basicAerodynamicForceVector);

class AerodynamicForceManager
{
public:
    AerodynamicForceManager();
    ~AerodynamicForceManager();
protected:
    size_t numberOfAerodynamicForceBodies;
    vector<AerodynamicForce *> aerodynamicForceContainer;
    bool isAllocateMemory;

    //! The corresponding Boundary Condition ID of each moving body.
    int *bodyID;
protected:
    void AllocateMemory();
    void DeallocateMemory();
public:
    void InitializeAerodynamicForceManager();
    void AddAerodynamicForce(ActionKey *actionKey);
    void CollectAerodynamicForceOfAllProcessors();

    void DumpAerodynamicForceCoefficient(ActionKey *actionKey);

    void SetNumberOfAerodynamicForces(int numberOfAerodynamicForceBodies) { this->numberOfAerodynamicForceBodies = numberOfAerodynamicForceBodies; }

    //! Dump all force data in forcefile.
    void DumpAllForceData(ActionKey *actkey);

public:
    AerodynamicForce *GetAerodynamicForce(int bodyIndex) { return aerodynamicForceContainer[bodyIndex]; }
    AerodynamicForce *GetAerodynamicForce() { return aerodynamicForceContainer[numberOfAerodynamicForceBodies]; }
protected:
    void ExtractForceDataFromActionKey(ActionKey *actionKey, vector< DataContainer * > *superContainerArray);
    vector< DataContainer * > *CreateDataContainerArray();
    void FreeDataContainerArray(vector< DataContainer * > *superContainerArray);
};

AerodynamicForceManager *GetAerodynamicForceManager();
AerodynamicForce *GetCurrentAerodynamicForce(int bodyIndex);

BasicAerodynamicForce *GetBasicAerodynamicForce(int iBoundaryConditionRegion, BasicAerodynamicForce **basicAerodynamicForceContainer);
BasicAerodynamicForce *GetBasicAerodynamicForce(int iBoundaryConditionRegion, vector< BasicAerodynamicForce * > *basicAerodynamicForceVector);

}
