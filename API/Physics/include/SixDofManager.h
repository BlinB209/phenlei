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
//! @file      SixDofManager.h
//! @brief     Explain this file briefly.
//! @author    xxx.

#pragma once

using namespace std;

namespace PHSPACE
{

namespace KINETIC_SCHEME
{
    const int ADMAS_BASHFORTH_FIRST_ORDER  = 0;
    const int ADMAS_BASHFORTH_SECOND_ORDER = 1;
    const int IMPLICIT_EULER_FIRST_ORDER   = 2;
    const int IMPLICIT_EULER_SECOND_ORDER  = 3;
    const int ADMAS_MOULTON_SECOND_ORDER   = 4;
    const int ADMAS_MOULTON_THIRD_ORDER    = 5;
}

class MassCharacter
{
protected:
    RDouble mass;
    RDouble momentOfInertiaTensorXX, momentOfInertiaTensorYY, momentOfInertiaTensorZZ, momentOfInertiaTensorXY, momentOfInertiaTensorXZ, momentOfInertiaTensorYZ;
public:
    MassCharacter() { mass = momentOfInertiaTensorXX = momentOfInertiaTensorYY = momentOfInertiaTensorZZ = momentOfInertiaTensorXY = momentOfInertiaTensorXZ = momentOfInertiaTensorYZ = 0.0; };
    MassCharacter(RDouble mass, RDouble momentOfInertiaTensorXX, RDouble momentOfInertiaTensorYY, RDouble momentOfInertiaTensorZZ, RDouble momentOfInertiaTensorXY, RDouble momentOfInertiaTensorXZ, RDouble momentOfInertiaTensorYZ);

    ~MassCharacter() {};
    
    void SetMassCharacter(RDouble mass, vector< RDouble > &massMatrix);
    void SetMassCharacter(RDouble mass, RDouble momentOfInertiaTensorXX, RDouble momentOfInertiaTensorYY, RDouble momentOfInertiaTensorZZ, RDouble momentOfInertiaTensorXY, RDouble momentOfInertiaTensorXZ, RDouble momentOfInertiaTensorYZ);
    void GetMassCharacter(RDouble &mass, RDouble &momentOfInertiaTensorXX, RDouble &momentOfInertiaTensorYY, RDouble &momentOfInertiaTensorZZ, RDouble &momentOfInertiaTensorXY, RDouble &momentOfInertiaTensorXZ, RDouble &momentOfInertiaTensorYZ);
    void GetMomentInertia(RDouble &momentOfInertiaTensorXX, RDouble &momentOfInertiaTensorYY, RDouble &momentOfInertiaTensorZZ, RDouble &momentOfInertiaTensorXY,RDouble &momentOfInertiaTensorXZ,RDouble &momentOfInertiaTensorYZ);
    void GetMomentInertia(vector< RDouble > &momentOfInertia);

    RDouble GetMass(){ return mass; };
    RDouble GetMomentOfInertiaTensorXX() { return momentOfInertiaTensorXX; };
    RDouble GetMomentOfInertiaTensorYY() { return momentOfInertiaTensorYY; };
    RDouble GetMomentOfInertiaTensorZZ() { return momentOfInertiaTensorZZ; };
    RDouble GetMomentOfInertiaTensorXY() { return momentOfInertiaTensorXY; };
    RDouble GetMomentOfInertiaTensorXZ() { return momentOfInertiaTensorXZ; };
    RDouble GetMomentOfInertiaTensorYZ() { return momentOfInertiaTensorYZ; };
public:
    RDouble ** GetMomentOfInertiaTensor();
    RDouble ** GetInverseMomentOfInertiaTensor();
};

class SixDofParameter
{
protected:
    RDouble referenceLength;
    RDouble massCenterX, massCenterY, massCenterZ;
    RDouble massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ;
    RDouble eulerAngleAlpha, eulerAngleBeta, eulerAngleGama;
    RDouble angularVelocityX, angularVelocityY, angularVelocityZ;
public:
    SixDofParameter();
    ~SixDofParameter() {};

    void SetMassCenterX(RDouble massCenterX){ this->massCenterX = massCenterX; };
    void SetMassCenterY(RDouble massCenterY){ this->massCenterY = massCenterY; };
    void SetMassCenterZ(RDouble massCenterZ){ this->massCenterZ = massCenterZ; };
    void SetMassCenter (RDouble massCenterX, RDouble massCenterY, RDouble massCenterZ)
    {
        this->massCenterX = massCenterX; 
        this->massCenterY = massCenterY; 
        this->massCenterZ = massCenterZ; 
    };

    void SetMassCenterVelocityX(RDouble massCenterVelocityX){ this->massCenterVelocityX = massCenterVelocityX; };
    void SetMassCenterVelocityY(RDouble massCenterVelocityY){ this->massCenterVelocityY = massCenterVelocityY; };
    void SetMassCenterVelocityZ(RDouble massCenterVelocityZ){ this->massCenterVelocityZ = massCenterVelocityZ; };
    void SetMassCenterVelocity (RDouble massCenterVelocityX, RDouble massCenterVelocityY, RDouble massCenterVelocityZ)
    {
        this->massCenterVelocityX = massCenterVelocityX;
        this->massCenterVelocityY = massCenterVelocityY;
        this->massCenterVelocityZ = massCenterVelocityZ;
    };

    void SetEulerAngleAlpha(RDouble eulerAngleAlpha) { this->eulerAngleAlpha = eulerAngleAlpha; };
    void SetEulerAngleBeta (RDouble eulerAngleBeta)  { this->eulerAngleBeta  = eulerAngleBeta; };
    void SetEulerAngleGama (RDouble eulerAngleGama)  { this->eulerAngleGama  = eulerAngleGama; };

    void SetEulerAngle     (RDouble eulerAngleGama, RDouble eulerAngleBeta, RDouble eulerAngleAlpha) 
    {
        this->eulerAngleGama  = eulerAngleGama;
        this->eulerAngleBeta  = eulerAngleBeta;
        this->eulerAngleAlpha = eulerAngleAlpha;
    };

    void SetAngularVelocityX(RDouble angularVelocityX) { this->angularVelocityX = angularVelocityX; };
    void SetAngularVelocityY(RDouble angularVelocityY) { this->angularVelocityY = angularVelocityY; };
    void SetAngularVelocityZ(RDouble angularVelocityZ) { this->angularVelocityZ = angularVelocityZ; };

    void SetAngularVelocity(RDouble angularVelocityX, RDouble angularVelocityY, RDouble angularVelocityZ)
    { 
        this->angularVelocityX = angularVelocityX;
        this->angularVelocityY = angularVelocityY;
        this->angularVelocityZ = angularVelocityZ;
    };

    RDouble GetMassCenterX() { return massCenterX; };
    RDouble GetMassCenterY() { return massCenterY; };
    RDouble GetMassCenterZ() { return massCenterZ; };
    void GetMassCenter(RDouble &massCenterX, RDouble &massCenterY, RDouble &massCenterZ)
    {
        massCenterX = this->massCenterX;
        massCenterY = this->massCenterY;
        massCenterZ = this->massCenterZ; 
    };

    RDouble GetNondimensionalMassCenterX();
    RDouble GetNondimensionalMassCenterY();
    RDouble GetNondimensionalMassCenterZ();
    void GetNondimensionalMassCenterVector(vector< RDouble > &nondimensionalMassCenterVector);
    void GetNondimensionalMassCenter(RDouble &dimensionlessMassCenterX, RDouble &dimensionlessMassCenterY, RDouble &dimensionlessMassCenterZ);

    RDouble GetMassCenterVelocityX() { return massCenterVelocityX; };
    RDouble GetMassCenterVelocityY() { return massCenterVelocityY; };
    RDouble GetMassCenterVelocityZ() { return massCenterVelocityZ; };

    void GetMassCenterVelocity(RDouble &massCenterVelocityX, RDouble &massCenterVelocityY, RDouble &massCenterVelocityZ)
    {
        massCenterVelocityX = this->massCenterVelocityX;
        massCenterVelocityY = this->massCenterVelocityY;
        massCenterVelocityZ = this->massCenterVelocityZ; 
    };

    RDouble GetEulerAngleAlpha() { return eulerAngleAlpha; };
    RDouble GetEulerAngleBeta()  { return eulerAngleBeta ; };
    RDouble GetEulerAngleGama()  { return eulerAngleGama ; };

    void GetEulerAngle(RDouble &eulerAngleGama, RDouble &eulerAngleBeta, RDouble &eulerAngleAlpha) 
    { 
        eulerAngleGama  = this->eulerAngleGama; 
        eulerAngleBeta  = this->eulerAngleBeta; 
        eulerAngleAlpha = this->eulerAngleAlpha;
    };

    void GetEulerAngleContainer(vector< RDouble > &eulerAngleContainer)
    {
        return this->GetEulerAngle(eulerAngleContainer[ 0 ], eulerAngleContainer[ 1 ], eulerAngleContainer[ 2 ]);
    };

    RDouble GetAngularVelocityX() { return angularVelocityX; };
    RDouble GetAngularVelocityY() { return angularVelocityY; };
    RDouble GetAngularVelocityZ() { return angularVelocityZ; };
    void GetAngularVelocity(RDouble &angularVelocityX, RDouble &angularVelocityY, RDouble &angularVelocityZ)
    {
        angularVelocityX = this->angularVelocityX;
        angularVelocityY = this->angularVelocityY;
        angularVelocityZ = this->angularVelocityZ;
    };

    void GetAngularVelocityContainer(vector< RDouble > & angularVelocityContainer)
    {
        return this->GetAngularVelocity(angularVelocityContainer[ 0 ], angularVelocityContainer[ 1 ], angularVelocityContainer[ 2 ]);
    }
public:
    void GetGeneralVariableContainer(vector< RDouble > &generalVariableContainer);
    void SetGeneralVariableContainer(RDouble *generalVariableContainer);
    void SetGeneralVariableContainer(vector< RDouble > &massCenter, vector< RDouble > &attitudeAngle, 
        vector< RDouble > &massCenterVelocity, vector< RDouble > &angularVelocity);
    void UpdateEulerAngle           (RDouble *eulerAngleVariationContainer);
    void UpdateAngularVelocity      (RDouble *angularVelocityVariationContainer);
};

class SixDofSolver;
class SixDofSolverManager
{
protected:
    int numberOfMovingBodies;
    SixDofSolver **sixDofSolverContainer;

public:
    SixDofSolverManager ();
    ~SixDofSolverManager();

public:
    void Initialize();
    void Restart();
    void Run();
    void Post();
    void Dump();
    void DumpRestart(fstream &file);
    void ResetMassCenter();
protected:
    void AllocateSixDofSolverContainer();

public:
    MassCharacter   * GetMassCharacter    (int bodyIndex);
    SixDofParameter * GetSixDofParameterN1(int bodyIndex);
    SixDofParameter * GetSixDofParameterN (int bodyIndex);
    SixDofParameter * GetSixDofParameter  (int bodyIndex);

    SixDofSolver * GetSixDofSolver(int bodyIndex) { return sixDofSolverContainer[ bodyIndex ]; }
};

class SixDofSolver
{
public:
    SixDofSolver(SixDofSolverManager *sixDofSolverManager);
    virtual ~SixDofSolver();

protected:
    int bodyIndex;
    int numberOfTimeLayers;
    int resetMassCenter;
    
    SixDofSolverManager *sixDofSolverManager;

    vector< SixDofParameter * > sixDofParameterContainer;
    vector< vector< RDouble > > sixDofForceContainer;
    MassCharacter * massCharacter;

    //Ω«∂»ª°∂»ªªÀ„ 
    RDouble DegToRad, RadToDeg;

public:
    virtual void Initialize(int bodyIndexIn);
    void Restart();
    void Run();
    void Post();
    void UpdateUnsteadyTerm();

protected:
    void AllocateAllVariables();
    void AllocateSixDofParameter();
    void AllocateSixDofForceContainer();
    void AllocateMassCharacter();

    void InitializeMassCharacter();
    void InitializeSixDofByParameter();

    void DeallocateAllVariables();
public:
    virtual void InitializeSixDofByFiles(fstream &file);
    virtual void DumpRestartSixDofSolver(fstream &file);

    virtual void ResetMassCenter();
public:
    SixDofParameter * GetSixDofParameter(int iTimeLayer) { return sixDofParameterContainer[ iTimeLayer ]; };
    vector< RDouble > & GetForceVector    (int iTimeLayer) { return sixDofForceContainer    [ iTimeLayer ]; }
    MassCharacter   * GetMassCharacter() { return massCharacter; };

    int  GetNumberOfTimeLayers() { return numberOfTimeLayers; }
    void SetNumberOfTimeLayers(int numberOfTimeLayers) { this->numberOfTimeLayers = numberOfTimeLayers; }

protected:
    void MovingByFile                                    ();
    void Stationary                                      ();
    void SolveSixDofEquation                             ();
    void ModifyForceAndSixDofInformationForXYSymmetryCase();

    void OscillatingInXDirection();
    void OscillatingInYDirection();
    void Pitching               ();
    void Yawing                 ();
    void Rolling                ();
    void XDirectionOnly         ();
    void YDirectionOnly         ();
    void RotationInXDirection   ();
    void RotationInZDirection   ();
    
    void ModifyForceAndSixDofInformationForPitchingOnly          ();
    void ModifyForceAndSixDofInformationForRollingOnly           ();
    void ModifyForceAndSixDofInformationForPitchingAndRollingOnly();
    
    void ReadSixDofParameterFromFile();
    void Rotation_Rudder_Test1      ();

protected:
    void GetSixDofFlux(SixDofParameter *sixDofParameter, MassCharacter *massCharacter, vector< RDouble > &force, vector< RDouble > &generalSixDofFlux);
    
    void ComputeRightHandSideOfRigidBodyTranslationalMotionEquation  (SixDofParameter *sixDofParameter, MassCharacter *massCharacter, vector< RDouble > &force, vector< RDouble > &generalSixDofFlux);
    void ComputeRightHandSideOfRigidBodyTranslationalDynamicsEquation(SixDofParameter *sixDofParameter, MassCharacter *massCharacter, vector< RDouble > &force, vector< RDouble > &generalSixDofFlux);
    void ComputeRightHandSideOfRigidBodyRotationalMotionEquation     (SixDofParameter *sixDofParameter, MassCharacter *massCharacter, vector< RDouble > &force, vector< RDouble > &generalSixDofFlux);
    void ComputeRightHandSideOfRigidBodyRotationalDynamicsEquation   (SixDofParameter *sixDofParameter, MassCharacter *massCharacter, vector< RDouble > &force, vector< RDouble > &generalSixDofFlux);

    void ComputeSixDofForceTimeSequenceContainer();
    void GetDimensionalTotalForceAndMomentOfBody(vector< RDouble > &dimensionalGeneralForceContainer, int iTimeLayer);
    void GetMassCenterAndEulerAngleContainer(vector< RDouble > &massCenterContainer, vector< RDouble > &eulerAngleContainer, int iTimeLayer);
    void GetMomentAboutMassCenterInInertialFrame(vector< RDouble > &massCenterContainer, vector< RDouble > &totalForceVector, vector< RDouble > &momentVector, vector< RDouble > &momentAboutMassCenterInInertialFrame);
    void GetDimensionalAdditionalForceAndMomentOfBody(vector< RDouble > &dimensionalGeneralForceContainer, int iTimeLayer);
    void GetDimensionalAdditionalForceAndMomentOfBodyByParameter(vector< RDouble > &dimensionalGeneralForceContainer, int iTimeLayer);
    void ComputeCurrentSixDofVariables(vector< vector< RDouble > > &qContainer, vector< vector< RDouble > > &fluxContainer, int numberOfEquations);
    void ComputeCurrentSixDofVariables(RDouble *q, RDouble *qN, RDouble *qN1, RDouble *fluxM, RDouble *fluxN, RDouble *fluxN1, int numberOfEquations);
    void GetDimensionalAerodynamicForceAndMoment(vector< RDouble > &dimensionalAerodynamicForce, vector< RDouble > &dimensionalAerodynamicMoment, int bodyIndex, int iTimeLayer);
public:
    void GetDimensionalAerodynamicForceAndMomentOfBody   (vector< RDouble > &dimensionalGeneralAerodynamicForceContainer, int iTimeLayer);
    void GetDimensionalAerodynamicForceAndMomentBodyFrame(vector< RDouble > &dimensionalGeneralAerodynamicForceContainer, int iTimeLayer);
    void GetDimensionalAerodynamicForceAndMomentInerFrame(vector< RDouble > &dimensionalGeneralAerodynamicForceContainer, int iTimeLayer);

public:
    void DumpSixDofSolverInformation();

protected:
    void DumpSixDofParameterInformation();
    void DumpSixDofParameterInformation(fstream &file);

    void DumpAleForceInformation();
    void DumpAleForceInformationBodyFrame(fstream &file);
    void DumpAleForceInformationInerFrame(fstream &file);
};

SixDofParameter * GetSixDofParameter(int bodyIndex);

void GetSixDofCoefficients(RDouble &a, RDouble &b, RDouble &c, RDouble &d, RDouble &e, RFloat &relaxParameter);

void ConvertVectorFromInertialCoordinateSystemToBodyCoordinateSystem(vector< RDouble > &vectorInInertialFrame, vector< RDouble > &angleVector, vector< RDouble > &massCenter, vector< RDouble > &vectorInBodyFrame);
void ConvertVectorFromBodyCoordinateSystemToInertialCoordinateSystem(vector< RDouble > &vectorInBodyFrame,     vector< RDouble > &angleVector, vector< RDouble > &massCenter, vector< RDouble > &vectorInInertialFrame);
void ConvertVectorFromInertialCoordinateSystemToBodyCoordinateSystem(vector< RDouble > &vectorInInertialFrame, vector< RDouble > &angleVector, vector< RDouble > &vectorInBodyFrame);
void ConvertVectorFromBodyCoordinateSystemToInertialCoordinateSystem(vector< RDouble > &vectorInBodyFrame,     vector< RDouble > &angleVector, vector< RDouble > &vectorInInertialFrame);

void RotateAboutAxisX(vector< RDouble > &vectorIn, RDouble rotateAngle, vector< RDouble > &vectorOut);
void RotateAboutAxisY(vector< RDouble > &vectorIn, RDouble rotateAngle, vector< RDouble > &vectorOut);
void RotateAboutAxisZ(vector< RDouble > &vectorIn, RDouble rotateAngle, vector< RDouble > &vectorOut);

}