#include "PHMpi.h"
#include "AleManager.h"
#include "Force.h"
#include "SixDofManager.h"
#include "Pointer.h"
#include "PHMatrix.h"
#include "GlobalDataBase.h"
#include "Constants.h"
#include "Force.h"
#include "AleForceManager.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "IO_FileReader.h"
#include "TK_Exit.h"
#include "TK_Log.h"
#include "Math_BasisFunction.h"

namespace PHSPACE
{

RDouble ** MassCharacter::GetMomentOfInertiaTensor()
{
    RDouble ** momentOfInertiaTensor = PHSPACE::NewPointer2< RDouble > (3, 3);

    momentOfInertiaTensor[ 0 ][ 0 ] =   momentOfInertiaTensorXX;
    momentOfInertiaTensor[ 0 ][ 1 ] = - momentOfInertiaTensorXY;
    momentOfInertiaTensor[ 0 ][ 2 ] = - momentOfInertiaTensorXZ;

    momentOfInertiaTensor[ 1 ][ 0 ] = - momentOfInertiaTensorXY;
    momentOfInertiaTensor[ 1 ][ 1 ] =   momentOfInertiaTensorYY;
    momentOfInertiaTensor[ 1 ][ 2 ] = - momentOfInertiaTensorYZ;

    momentOfInertiaTensor[ 2 ][ 0 ] = - momentOfInertiaTensorXZ;
    momentOfInertiaTensor[ 2 ][ 1 ] = - momentOfInertiaTensorYZ;
    momentOfInertiaTensor[ 2 ][ 2 ] =   momentOfInertiaTensorZZ;

    return momentOfInertiaTensor;
}

RDouble ** MassCharacter::GetInverseMomentOfInertiaTensor()
{
    RDouble ** inverseMomentOfInertiaTensor = PHSPACE::NewPointer2< RDouble > (3, 3);

    RDouble d   =       momentOfInertiaTensorXX * momentOfInertiaTensorYY * momentOfInertiaTensorZZ -
                        momentOfInertiaTensorXX * momentOfInertiaTensorYZ * momentOfInertiaTensorYZ -
                        momentOfInertiaTensorYY * momentOfInertiaTensorXZ * momentOfInertiaTensorXZ -
                        momentOfInertiaTensorZZ * momentOfInertiaTensorXY * momentOfInertiaTensorXY -
                  2.0 * momentOfInertiaTensorXY * momentOfInertiaTensorXZ * momentOfInertiaTensorYZ;

    RDouble a11 = momentOfInertiaTensorYY * momentOfInertiaTensorZZ - momentOfInertiaTensorYZ * momentOfInertiaTensorYZ;
    RDouble a12 = momentOfInertiaTensorXY * momentOfInertiaTensorZZ + momentOfInertiaTensorYZ * momentOfInertiaTensorXZ;
    RDouble a13 = momentOfInertiaTensorXY * momentOfInertiaTensorYZ + momentOfInertiaTensorXZ * momentOfInertiaTensorYY;

    RDouble a21 = momentOfInertiaTensorXY * momentOfInertiaTensorZZ + momentOfInertiaTensorYZ * momentOfInertiaTensorXZ;
    RDouble a22 = momentOfInertiaTensorXX * momentOfInertiaTensorZZ - momentOfInertiaTensorXZ * momentOfInertiaTensorXZ;
    RDouble a23 = momentOfInertiaTensorXX * momentOfInertiaTensorYZ + momentOfInertiaTensorXZ * momentOfInertiaTensorXY;

    RDouble a31 = momentOfInertiaTensorXY * momentOfInertiaTensorYZ + momentOfInertiaTensorXZ * momentOfInertiaTensorYY;
    RDouble a32 = momentOfInertiaTensorXX * momentOfInertiaTensorYZ + momentOfInertiaTensorXZ * momentOfInertiaTensorXY;
    RDouble a33 = momentOfInertiaTensorXX * momentOfInertiaTensorYY - momentOfInertiaTensorXY * momentOfInertiaTensorXY;

    inverseMomentOfInertiaTensor[ 0 ][ 0 ] = a11 / d;
    inverseMomentOfInertiaTensor[ 0 ][ 1 ] = a12 / d;
    inverseMomentOfInertiaTensor[ 0 ][ 2 ] = a13 / d;

    inverseMomentOfInertiaTensor[ 1 ][ 0 ] = a21 / d;
    inverseMomentOfInertiaTensor[ 1 ][ 1 ] = a22 / d;
    inverseMomentOfInertiaTensor[ 1 ][ 2 ] = a23 / d;

    inverseMomentOfInertiaTensor[ 2 ][ 0 ] = a31 / d;
    inverseMomentOfInertiaTensor[ 2 ][ 1 ] = a32 / d;
    inverseMomentOfInertiaTensor[ 2 ][ 2 ] = a33 / d;

    return inverseMomentOfInertiaTensor;
}

MassCharacter::MassCharacter(RDouble mass, RDouble momentOfInertiaTensorXX, RDouble momentOfInertiaTensorYY, RDouble momentOfInertiaTensorZZ, RDouble momentOfInertiaTensorXY, RDouble momentOfInertiaTensorXZ, RDouble momentOfInertiaTensorYZ)
{
    this->mass = mass;
    this->momentOfInertiaTensorXX  = momentOfInertiaTensorXX;
    this->momentOfInertiaTensorYY  = momentOfInertiaTensorYY;
    this->momentOfInertiaTensorZZ  = momentOfInertiaTensorZZ;
    this->momentOfInertiaTensorXY  = momentOfInertiaTensorXY;
    this->momentOfInertiaTensorXZ  = momentOfInertiaTensorXZ;
    this->momentOfInertiaTensorYZ  = momentOfInertiaTensorYZ;
}

void MassCharacter::SetMassCharacter(RDouble mass, vector< RDouble > & massMatrix)
{
    this->mass = mass;
    this->momentOfInertiaTensorXX  = massMatrix[ 0 ];
    this->momentOfInertiaTensorYY  = massMatrix[ 1 ];
    this->momentOfInertiaTensorZZ  = massMatrix[ 2 ];
    this->momentOfInertiaTensorXY  = massMatrix[ 3 ];
    this->momentOfInertiaTensorXZ  = massMatrix[ 4 ];
    this->momentOfInertiaTensorYZ  = massMatrix[ 5 ];
}

void MassCharacter::SetMassCharacter(RDouble mass, RDouble momentOfInertiaTensorXX, RDouble momentOfInertiaTensorYY, RDouble momentOfInertiaTensorZZ, RDouble momentOfInertiaTensorXY, RDouble momentOfInertiaTensorXZ, RDouble momentOfInertiaTensorYZ)
{
    this->mass = mass;
    this->momentOfInertiaTensorXX  = momentOfInertiaTensorXX;
    this->momentOfInertiaTensorYY  = momentOfInertiaTensorYY;
    this->momentOfInertiaTensorZZ  = momentOfInertiaTensorZZ;
    this->momentOfInertiaTensorXY  = momentOfInertiaTensorXY;
    this->momentOfInertiaTensorXZ  = momentOfInertiaTensorXZ;
    this->momentOfInertiaTensorYZ  = momentOfInertiaTensorYZ;
}

void MassCharacter::GetMassCharacter(RDouble & mass, RDouble & momentOfInertiaTensorXX, RDouble & momentOfInertiaTensorYY, RDouble & momentOfInertiaTensorZZ,RDouble & momentOfInertiaTensorXY,RDouble & momentOfInertiaTensorXZ,RDouble & momentOfInertiaTensorYZ)
{
    mass = this->mass;
    momentOfInertiaTensorXX  = this->momentOfInertiaTensorXX;
    momentOfInertiaTensorYY  = this->momentOfInertiaTensorYY;
    momentOfInertiaTensorZZ  = this->momentOfInertiaTensorZZ;
    momentOfInertiaTensorXY  = this->momentOfInertiaTensorXY;
    momentOfInertiaTensorXZ  = this->momentOfInertiaTensorXZ;
    momentOfInertiaTensorYZ  = this->momentOfInertiaTensorYZ;
}

void MassCharacter::GetMomentInertia(vector< RDouble > & momentOfInertia)
{
    momentOfInertia[ 0 ]  = this->momentOfInertiaTensorXX;
    momentOfInertia[ 1 ]  = this->momentOfInertiaTensorYY;
    momentOfInertia[ 2 ]  = this->momentOfInertiaTensorZZ;
    momentOfInertia[ 3 ]  = this->momentOfInertiaTensorXY;
    momentOfInertia[ 4 ]  = this->momentOfInertiaTensorXZ;
    momentOfInertia[ 5 ]  = this->momentOfInertiaTensorYZ;
}

void MassCharacter::GetMomentInertia(RDouble & momentOfInertiaTensorXX, RDouble & momentOfInertiaTensorYY, RDouble & momentOfInertiaTensorZZ, RDouble & momentOfInertiaTensorXY, RDouble & momentOfInertiaTensorXZ, RDouble & momentOfInertiaTensorYZ)
{
    momentOfInertiaTensorXX  = this->momentOfInertiaTensorXX;
    momentOfInertiaTensorYY  = this->momentOfInertiaTensorYY;
    momentOfInertiaTensorZZ  = this->momentOfInertiaTensorZZ;
    momentOfInertiaTensorXY  = this->momentOfInertiaTensorXY;
    momentOfInertiaTensorXZ  = this->momentOfInertiaTensorXZ;
    momentOfInertiaTensorYZ  = this->momentOfInertiaTensorYZ;
}

SixDofParameter::SixDofParameter()
{
    massCenterX = 0.0;
    massCenterY = 0.0;
    massCenterZ = 0.0;

    eulerAngleAlpha = 0.0;
    eulerAngleBeta  = 0.0;
    eulerAngleGama  = 0.0; 

    massCenterVelocityX    = 0.0;
    massCenterVelocityY    = 0.0;
    massCenterVelocityZ    = 0.0;

    angularVelocityX = 0.0;
    angularVelocityY = 0.0;
    angularVelocityZ = 0.0; 

    referenceLength   = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("referenceLength");
}

void SixDofParameter::GetNondimensionalMassCenter(RDouble & dimensionlessMassCenterX, RDouble & dimensionlessMassCenterY, RDouble & dimensionlessMassCenterZ)
{
    dimensionlessMassCenterX = this->massCenterX / referenceLength;
    dimensionlessMassCenterY = this->massCenterY / referenceLength;
    dimensionlessMassCenterZ = this->massCenterZ / referenceLength;
}

RDouble SixDofParameter::GetNondimensionalMassCenterX()
{
    return this->massCenterX / referenceLength;
}

RDouble SixDofParameter::GetNondimensionalMassCenterY()
{
    return this->massCenterY / referenceLength;
}

RDouble SixDofParameter::GetNondimensionalMassCenterZ()
{
    return this->massCenterZ / referenceLength;
}

void SixDofParameter::GetNondimensionalMassCenterVector(vector< RDouble > & nondimensionalMassCenterVector)
{
    nondimensionalMassCenterVector[ 0 ] = this->massCenterX / referenceLength;
    nondimensionalMassCenterVector[ 1 ] = this->massCenterY / referenceLength;
    nondimensionalMassCenterVector[ 2 ] = this->massCenterZ / referenceLength;
}

void SixDofParameter::GetGeneralVariableContainer(vector< RDouble > & generalVariableContainer)
{
    this->GetMassCenter        (generalVariableContainer[ 0 ], generalVariableContainer[ 1  ], generalVariableContainer[ 2  ]);
    this->GetMassCenterVelocity(generalVariableContainer[ 3 ], generalVariableContainer[ 4  ], generalVariableContainer[ 5  ]);
    this->GetEulerAngle        (generalVariableContainer[ 6 ], generalVariableContainer[ 7  ], generalVariableContainer[ 8  ]);
    this->GetAngularVelocity   (generalVariableContainer[ 9 ], generalVariableContainer[ 10 ], generalVariableContainer[ 11 ]);
}

void SixDofParameter::SetGeneralVariableContainer(RDouble * generalVariableContainer)
{
    this->SetMassCenter        (generalVariableContainer[ 0 ], generalVariableContainer[ 1  ], generalVariableContainer[ 2  ]);
    this->SetMassCenterVelocity(generalVariableContainer[ 3 ], generalVariableContainer[ 4  ], generalVariableContainer[ 5  ]);
    this->SetEulerAngle        (generalVariableContainer[ 6 ], generalVariableContainer[ 7  ], generalVariableContainer[ 8  ]);
    this->SetAngularVelocity   (generalVariableContainer[ 9 ], generalVariableContainer[ 10 ], generalVariableContainer[ 11 ]);
}

void SixDofParameter::SetGeneralVariableContainer(vector< RDouble > & massCenter, vector< RDouble >& attitudeAngle, 
    vector< RDouble >& massCenterVelocity, vector< RDouble >& angularVelocity)
{
    this->SetMassCenter        (massCenter        [0], massCenter        [1],         massCenter[2]);
    this->SetEulerAngle        (attitudeAngle     [0], attitudeAngle     [1],      attitudeAngle[2]);
    this->SetMassCenterVelocity(massCenterVelocity[0], massCenterVelocity[1], massCenterVelocity[2]);
    this->SetAngularVelocity   (angularVelocity   [0], angularVelocity   [1],    angularVelocity[2]);
}

void SixDofParameter::UpdateEulerAngle(RDouble * eulerAngleVariationContainer)
{
    eulerAngleGama  += eulerAngleVariationContainer[ 0 ];
    eulerAngleBeta  += eulerAngleVariationContainer[ 1 ];
    eulerAngleAlpha += eulerAngleVariationContainer[ 2 ];
}

void SixDofParameter::UpdateAngularVelocity(RDouble * angularVelocityVariationContainer)
{
    angularVelocityX += angularVelocityVariationContainer[ 0 ];
    angularVelocityY += angularVelocityVariationContainer[ 1 ];
    angularVelocityZ += angularVelocityVariationContainer[ 2 ];
}

SixDofSolver::SixDofSolver(SixDofSolverManager * sixDofSolverManager)
{
    this->sixDofSolverManager = sixDofSolverManager;

    bodyIndex = 0;

    resetMassCenter = 0;

    AllocateAllVariables();

    DegToRad = PI / 180.0;
    RadToDeg = 180.0 / PI;
}

SixDofSolver::~SixDofSolver()
{
    DeallocateAllVariables();
}

void SixDofSolver::AllocateAllVariables()
{
    AllocateSixDofParameter();

    AllocateSixDofForceContainer();

    AllocateMassCharacter();
}

void SixDofSolver::DeallocateAllVariables()
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();
    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        delete sixDofParameterContainer[ iTimeLayer ];
        sixDofParameterContainer[iTimeLayer] = nullptr;
    }

    delete massCharacter;
    massCharacter = nullptr;
}

void SixDofSolver::AllocateSixDofParameter()
{
    int numberOfTimeLayers = 3;
    SetNumberOfTimeLayers(numberOfTimeLayers);

    sixDofParameterContainer.resize(numberOfTimeLayers);

    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        sixDofParameterContainer[ iTimeLayer ] = new SixDofParameter();
    }
}

void SixDofSolver::AllocateMassCharacter()
{
    massCharacter = new MassCharacter();
}

void SixDofSolver::AllocateSixDofForceContainer()
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();

    sixDofForceContainer.resize(numberOfTimeLayers);

    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        sixDofForceContainer[ iTimeLayer ].resize(6);
    }
}

void SixDofSolver::InitializeMassCharacter()
{
    MassCharacter * massCharacter = this->GetMassCharacter();

    RDouble mass = GetRDoubleParameterFromDataBase(bodyIndex, "mass");

    vector< RDouble > massMatrix;
    GetRDoubleVectorFromDataBase(massMatrix, bodyIndex, "massMatrix", 6);

    massCharacter->SetMassCharacter(mass, massMatrix);
}

void SixDofSolver::Initialize(int bodyIndexIn)
{
    this->bodyIndex = bodyIndexIn;

    this->InitializeMassCharacter();

    this->InitializeSixDofByParameter();
}

void SixDofSolver::Restart()
{
    fstream & file = PHSPACE::GetCurrentAleFile();
    this->InitializeSixDofByFiles(file);
}

void SixDofSolver::InitializeSixDofByParameter()
{
    vector< RDouble > massCenter;
    GetRDoubleVectorFromDataBase(massCenter, bodyIndex, "massCenter"  , 3);

    vector< RDouble > attitudeAngle;
    GetRDoubleVectorFromDataBase(attitudeAngle, bodyIndex, "attitudeAngle", 3);

    vector< RDouble > massCenterVelocity;
    GetRDoubleVectorFromDataBase(massCenterVelocity, bodyIndex, "massCenterVelocity", 3);

    vector< RDouble > angularVelocity;
    GetRDoubleVectorFromDataBase(angularVelocity, bodyIndex, "angularVelocity", 3);

    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

    massCenter[ 0 ] /= reynoldsReferenceLengthDimensional;
    massCenter[ 1 ] /= reynoldsReferenceLengthDimensional;
    massCenter[ 2 ] /= reynoldsReferenceLengthDimensional;

    attitudeAngle[ 0 ] *= PI / 180;
    attitudeAngle[ 1 ] *= PI / 180;
    attitudeAngle[ 2 ] *= PI / 180;

    angularVelocity[ 0 ] *= PI / 180;
    angularVelocity[ 1 ] *= PI / 180;
    angularVelocity[ 2 ] *= PI / 180;

    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        SixDofParameter * sixDofParameter = this->GetSixDofParameter(iTimeLayer);
        sixDofParameter->SetGeneralVariableContainer(massCenter, attitudeAngle, massCenterVelocity, angularVelocity);
    }
}

void SixDofSolver::InitializeSixDofByFiles(fstream & file)
{
    const int numberOfEquations = 12;
    RDouble generalVariableContainer[ numberOfEquations ];

    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        //PHSPACE::PHRead(file, generalVariableContainer, numberOfEquations);
        for (int iVar = 0; iVar < numberOfEquations; ++ iVar)
        {
            file >> generalVariableContainer[iVar];
        }

        sixDofParameterContainer[iTimeLayer]->SetGeneralVariableContainer(generalVariableContainer);
    }
}

void SixDofSolver::ResetMassCenter()
{
    int resetMassCenter = GetIntegerParameterFromDataBaseIfExist(bodyIndex, "resetMassCenter", 0);
    this -> resetMassCenter = resetMassCenter;
    if (resetMassCenter == 0)
    {
        return;
    }

    vector< RDouble > massCenterDxyz;
    GetRDoubleVectorFromDataBase(massCenterDxyz, bodyIndex, "massCenterDxyz", 3);

    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

    massCenterDxyz[0] /= reynoldsReferenceLengthDimensional;
    massCenterDxyz[1] /= reynoldsReferenceLengthDimensional;
    massCenterDxyz[2] /= reynoldsReferenceLengthDimensional;

    vector< RDouble > massCenter;
    vector< RDouble > angleVector;
    vector< RDouble > massCenterVelocity;
    vector< RDouble > angularVelocity;

    massCenter.resize(3);
    angleVector.resize(3);
    massCenterVelocity.resize(3);
    angularVelocity.resize(3);

    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        SixDofParameter *sixDofParameter = this->GetSixDofParameter(iTimeLayer);

        sixDofParameter->GetEulerAngle(angleVector[0], angleVector[1], angleVector[2]);
        sixDofParameter->GetMassCenter(massCenter[0], massCenter[1], massCenter[2]);

        sixDofParameter->GetMassCenterVelocity(massCenterVelocity[0], massCenterVelocity[1], massCenterVelocity[2]);
        sixDofParameter->GetAngularVelocity(angularVelocity[0], angularVelocity[1], angularVelocity[2]);

        ConvertVectorFromBodyCoordinateSystemToInertialCoordinateSystem(massCenterDxyz, angleVector, massCenterDxyz);
        massCenter[0] += massCenterDxyz[0];
        massCenter[1] += massCenterDxyz[1];
        massCenter[2] += massCenterDxyz[2];
        sixDofParameter->SetMassCenter(massCenter[0], massCenter[1], massCenter[2]);

        massCenterVelocity[0] += angularVelocity[1] * massCenterDxyz[2] - angularVelocity[2] * massCenterDxyz[1] ;
        massCenterVelocity[1] += angularVelocity[2] * massCenterDxyz[0] - angularVelocity[0] * massCenterDxyz[2] ;
        massCenterVelocity[2] += angularVelocity[0] * massCenterDxyz[1] - angularVelocity[1] * massCenterDxyz[0];
        sixDofParameter->SetMassCenterVelocity(massCenterVelocity[0], massCenterVelocity[1], massCenterVelocity[2]);
    }
}

void SixDofSolver::Run()
{
    ComputeSixDofForceTimeSequenceContainer();

    int RBDMethod = PHSPACE::GetIntegerParameterFromDataBase(bodyIndex, "RBDMethod");

    if (RBDMethod == - 1)//! given motion partten.
    {
        MovingByFile();
    }
    else if (RBDMethod == 0  )//! still.
    {
        Stationary();
    }
    else if (RBDMethod == 1  )//! six DOF motion.
    {
        SolveSixDofEquation();
    }
    else if (RBDMethod == 2  )//! three DOF motion.
    {
        ModifyForceAndSixDofInformationForXYSymmetryCase();
        SolveSixDofEquation();
    }
    else if (RBDMethod == 11 )//! X-axis forced motion.
    {
        OscillatingInXDirection();
    }
    else if (RBDMethod == 12 )//! Y-axis forced motion.
    {
        OscillatingInYDirection();
    }
    else if (RBDMethod == 13 )//! Z-axis forced motion.
    {
        TK_Exit::ExceptionExit("There is no Z-axis forced motion in temporary.");
        //OscillatingInZDirection();
    }
    else if (RBDMethod == 14 )//! forced pitch motion.
    {
        Pitching();
        //SolveSixDofEquation();
    }
    else if (RBDMethod == 15 )//! forced yaw motion.
    {
        Yawing();
    }
    else if (RBDMethod == 16 )//! forced roll motion.
    {
        Rolling();
    }
    else if (RBDMethod == 21 )//! single DOF motion.
    {
        XDirectionOnly();
        //SolveSixDofEquation();
    }
    else if ( RBDMethod == 17  )//! Forced rotation along X periodically;
    {
        RotationInXDirection();
    }
    else if (RBDMethod == 19)//! Forced rotation along Z periodically;
    {
        RotationInZDirection();
    }
    else if (RBDMethod == 22 )//! single DOF Y-axis motion
    {
        YDirectionOnly();
        //ModifyForceAndSixDofInformationForYDirectionOnly();
        //SolveSixDofEquation();
    }
    else if (RBDMethod == 23 )//! single DOF Z-axis motion
    {
        TK_Exit::ExceptionExit("There is no single DOF Z-axis motion in temporary.");
        //ModifyForceAndSixDofInformationForZDirectionOnly();
        //SolveSixDofEquation();
    }
    else if (RBDMethod == 24 )//! single DOF pitch motion.
    {
        ModifyForceAndSixDofInformationForPitchingOnly();
        SolveSixDofEquation();
    }
    else if (RBDMethod == 25 )//! single DOF yaw motion.
    {
        TK_Exit::ExceptionExit("There is no single DOF yaw motion. in temporary.");
        //ModifyForceAndSixDofInformationForYawingOnly();
        //SolveSixDofEquation();
    }
    else if (RBDMethod == 26 )//! single DOF roll motion.
    {
        ModifyForceAndSixDofInformationForRollingOnly();
        SolveSixDofEquation();
    }
    else if (RBDMethod == 31 )//! single DOF pitch and roll motion.
    {
        ModifyForceAndSixDofInformationForPitchingAndRollingOnly();
        SolveSixDofEquation();
    }
    else
    {
        TK_Exit::ExceptionExit("There is no this motion");
    }
}

void SixDofSolver::Post()
{
    UpdateUnsteadyTerm();
    DumpSixDofSolverInformation();
}

void SixDofSolver::UpdateUnsteadyTerm()
{
    SixDofParameter * sixDofParameter   = this->GetSixDofParameter(0);
    SixDofParameter * sixDofParameterN  = this->GetSixDofParameter(1);
    SixDofParameter * sixDofParameterN1 = this->GetSixDofParameter(2);

    (* sixDofParameterN1) = (* sixDofParameterN);
    (* sixDofParameterN ) = (* sixDofParameter );
}

void SixDofSolver::ComputeSixDofForceTimeSequenceContainer()
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();

    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        vector< RDouble > & forceVector = this->GetForceVector(iTimeLayer);
        this->GetDimensionalTotalForceAndMomentOfBody(forceVector, iTimeLayer);
    }
}

void SixDofSolver::GetDimensionalTotalForceAndMomentOfBody(vector< RDouble > & dimensionalGeneralForceContainer, int iTimeLayer)
{
    this->GetDimensionalAerodynamicForceAndMomentOfBody(dimensionalGeneralForceContainer, iTimeLayer);
    this->GetDimensionalAdditionalForceAndMomentOfBody (dimensionalGeneralForceContainer, iTimeLayer);
}

void SixDofSolver::GetMassCenterAndEulerAngleContainer(vector< RDouble > & massCenterContainer, vector< RDouble > & eulerAngleContainer, int iTimeLayer)
{
    SixDofParameter * sixDofParameter = this->GetSixDofParameter(iTimeLayer);

    sixDofParameter->GetMassCenter(massCenterContainer[ 0 ], massCenterContainer[ 1 ], massCenterContainer[ 2 ]);
    sixDofParameter->GetEulerAngleContainer(eulerAngleContainer);
}

void SixDofSolver::GetMomentAboutMassCenterInInertialFrame(vector< RDouble > & massCenterContainer, vector< RDouble > & totalForceVector, vector< RDouble > & momentVector, vector< RDouble > & momentAboutMassCenterInInertialFrame)
{
    RDouble momentChange [ 3 ];
    PHSPACE::CrossProduct(& massCenterContainer[ 0 ], & totalForceVector[ 0 ], momentChange);

    momentAboutMassCenterInInertialFrame[ 0 ] = momentVector[ 0 ] - momentChange[ 0 ];
    momentAboutMassCenterInInertialFrame[ 1 ] = momentVector[ 1 ] - momentChange[ 1 ];
    momentAboutMassCenterInInertialFrame[ 2 ] = momentVector[ 2 ] - momentChange[ 2 ];
}

void SixDofSolver::GetDimensionalAerodynamicForceAndMomentOfBody(vector< RDouble > & dimensionalGeneralAerodynamicForceContainer, int iTimeLayer)
{
    vector< RDouble > massCenterContainer(3);
    vector< RDouble > eulerAngleContainer(3);

    this->GetMassCenterAndEulerAngleContainer(massCenterContainer, eulerAngleContainer, iTimeLayer);

    vector< RDouble > dimensionalAerodynamicForce (3);
    vector< RDouble > dimensionalAerodynamicMoment(3);
    
    GetDimensionalAerodynamicForceAndMoment(dimensionalAerodynamicForce, dimensionalAerodynamicMoment, bodyIndex, iTimeLayer);

    vector< RDouble > dimensionalAerodynamicMomentAboutMassCenter(3);
    this->GetMomentAboutMassCenterInInertialFrame(massCenterContainer, dimensionalAerodynamicForce, dimensionalAerodynamicMoment, dimensionalAerodynamicMomentAboutMassCenter);

    vector< RDouble > dimensionalAerodynamicForceInBodyFrame(3);
    vector< RDouble > dimensionalAerodynamicMomentAboutMassCenterInBodyFrame(3);
    PHSPACE::ConvertVectorFromInertialCoordinateSystemToBodyCoordinateSystem(dimensionalAerodynamicForce, eulerAngleContainer, dimensionalAerodynamicForceInBodyFrame);
    PHSPACE::ConvertVectorFromInertialCoordinateSystemToBodyCoordinateSystem(dimensionalAerodynamicMomentAboutMassCenter, eulerAngleContainer, dimensionalAerodynamicMomentAboutMassCenterInBodyFrame);

    dimensionalGeneralAerodynamicForceContainer[ 0 ] = dimensionalAerodynamicForceInBodyFrame[ 0 ];
    dimensionalGeneralAerodynamicForceContainer[ 1 ] = dimensionalAerodynamicForceInBodyFrame[ 1 ];
    dimensionalGeneralAerodynamicForceContainer[ 2 ] = dimensionalAerodynamicForceInBodyFrame[ 2 ];

    dimensionalGeneralAerodynamicForceContainer[ 3 ] = dimensionalAerodynamicMomentAboutMassCenterInBodyFrame[ 0 ];
    dimensionalGeneralAerodynamicForceContainer[ 4 ] = dimensionalAerodynamicMomentAboutMassCenterInBodyFrame[ 1 ];
    dimensionalGeneralAerodynamicForceContainer[ 5 ] = dimensionalAerodynamicMomentAboutMassCenterInBodyFrame[ 2 ];
}

void SixDofSolver::GetDimensionalAerodynamicForceAndMomentInerFrame(vector< RDouble > & dimensionalGeneralAerodynamicForceContainer, int iTimeLayer)
{
    vector< RDouble > massCenterContainer(3);
    vector< RDouble > eulerAngleContainer(3);

    this->GetMassCenterAndEulerAngleContainer(massCenterContainer, eulerAngleContainer, iTimeLayer);

    vector< RDouble > dimensionalAerodynamicForce (3);
    vector< RDouble > dimensionalAerodynamicMoment(3);
    
    GetDimensionalAerodynamicForceAndMoment(dimensionalAerodynamicForce, dimensionalAerodynamicMoment, bodyIndex, iTimeLayer);

    vector< RDouble > dimensionalAerodynamicMomentAboutMassCenter(3);
    this->GetMomentAboutMassCenterInInertialFrame(massCenterContainer, dimensionalAerodynamicForce, dimensionalAerodynamicMoment, dimensionalAerodynamicMomentAboutMassCenter);
    
    dimensionalGeneralAerodynamicForceContainer[ 0 ] = dimensionalAerodynamicForce[ 0 ];
    dimensionalGeneralAerodynamicForceContainer[ 1 ] = dimensionalAerodynamicForce[ 1 ];
    dimensionalGeneralAerodynamicForceContainer[ 2 ] = dimensionalAerodynamicForce[ 2 ];

    dimensionalGeneralAerodynamicForceContainer[ 3 ] = dimensionalAerodynamicMomentAboutMassCenter[ 0 ];
    dimensionalGeneralAerodynamicForceContainer[ 4 ] = dimensionalAerodynamicMomentAboutMassCenter[ 1 ];
    dimensionalGeneralAerodynamicForceContainer[ 5 ] = dimensionalAerodynamicMomentAboutMassCenter[ 2 ];
}

void SixDofSolver::GetDimensionalAerodynamicForceAndMomentBodyFrame(vector< RDouble > & dimensionalGeneralAerodynamicForceContainer, int iTimeLayer)
{
    vector< RDouble > massCenterContainer(3);
    vector< RDouble > eulerAngleContainer(3);

    this->GetMassCenterAndEulerAngleContainer(massCenterContainer, eulerAngleContainer, iTimeLayer);

    vector< RDouble > dimensionalAerodynamicForce (3);
    vector< RDouble > dimensionalAerodynamicMoment(3);
    
    GetDimensionalAerodynamicForceAndMoment(dimensionalAerodynamicForce, dimensionalAerodynamicMoment, bodyIndex, iTimeLayer);

    vector< RDouble > dimensionalAerodynamicMomentAboutMassCenter(3);
    this->GetMomentAboutMassCenterInInertialFrame(massCenterContainer, dimensionalAerodynamicForce, dimensionalAerodynamicMoment, dimensionalAerodynamicMomentAboutMassCenter);
    
    vector< RDouble > dimensionalAerodynamicForceInBodyFrame(3);
    vector< RDouble > dimensionalAerodynamicMomentAboutMassCenterInBodyFrame(3);

    PHSPACE::ConvertVectorFromInertialCoordinateSystemToBodyCoordinateSystem(dimensionalAerodynamicForce, eulerAngleContainer, dimensionalAerodynamicForceInBodyFrame);
    PHSPACE::ConvertVectorFromInertialCoordinateSystemToBodyCoordinateSystem(dimensionalAerodynamicMomentAboutMassCenter, eulerAngleContainer, dimensionalAerodynamicMomentAboutMassCenterInBodyFrame);

    dimensionalGeneralAerodynamicForceContainer[ 0 ] = dimensionalAerodynamicForceInBodyFrame[ 0 ];
    dimensionalGeneralAerodynamicForceContainer[ 1 ] = dimensionalAerodynamicForceInBodyFrame[ 1 ];
    dimensionalGeneralAerodynamicForceContainer[ 2 ] = dimensionalAerodynamicForceInBodyFrame[ 2 ];

    dimensionalGeneralAerodynamicForceContainer[ 3 ] = dimensionalAerodynamicMomentAboutMassCenterInBodyFrame[ 0 ];
    dimensionalGeneralAerodynamicForceContainer[ 4 ] = dimensionalAerodynamicMomentAboutMassCenterInBodyFrame[ 1 ];
    dimensionalGeneralAerodynamicForceContainer[ 5 ] = dimensionalAerodynamicMomentAboutMassCenterInBodyFrame[ 2 ];
}

void SixDofSolver::GetDimensionalAdditionalForceAndMomentOfBodyByParameter(vector< RDouble > & dimensionalGeneralForceContainer, int iTimeLayer)
{
    vector< RDouble > massCenterContainer(3);
    vector< RDouble > eulerAngleContainer(3);

    this->GetMassCenterAndEulerAngleContainer(massCenterContainer, eulerAngleContainer, iTimeLayer);

    vector< RDouble > dimensionalAdditionalForce(3);
    vector< RDouble > dimensionalAdditionalMoment(3);

    for (int iDim = 0; iDim < 3; ++ iDim)
    {
        dimensionalAdditionalForce[iDim] = zero;
        dimensionalAdditionalMoment[iDim] = zero;
    }

    RDouble aleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("aleSimulationTime");
    RDouble  timeOfNstep = aleSimulationTime - GetAleTimeStep();

    RDouble addedForceTime = GetRDoubleParameterFromDataBaseIfExist(bodyIndex, "addedForceTime", zero);
    vector< RDouble > addedForcePosition(3);
    GetRDoubleVectorFromDataBase(addedForcePosition, bodyIndex, "addedForcePosition", 3);

    SixDofParameter *sixDofParameter = this->GetSixDofParameter(0);
    RDouble MassCenterY = sixDofParameter->GetMassCenterY();

    if (timeOfNstep < addedForceTime || MassCenterY > addedForcePosition[1])
    {
        GetRDoubleVectorFromDataBase(dimensionalAdditionalForce, bodyIndex, "addedForce", 3);

        GetRDoubleVectorFromDataBase(dimensionalAdditionalMoment, bodyIndex, "addedMoment", 3);
    }

    RDouble mass = GetRDoubleParameterFromDataBase(bodyIndex, "mass");

    RDouble gravity = GetRDoubleParameterFromDataBaseIfExist(bodyIndex, "gravity", zero);

    dimensionalAdditionalForce[1] -= gravity * mass;

    vector< RDouble > dimensionalAdditionalForceInBodyFrame(3);
    vector< RDouble > dimensionalAdditionalMomentInBodyFrame(3);

    PHSPACE::ConvertVectorFromInertialCoordinateSystemToBodyCoordinateSystem(dimensionalAdditionalForce, eulerAngleContainer, dimensionalAdditionalForceInBodyFrame);
    PHSPACE::ConvertVectorFromInertialCoordinateSystemToBodyCoordinateSystem(dimensionalAdditionalMoment, eulerAngleContainer, dimensionalAdditionalMomentInBodyFrame);

    dimensionalGeneralForceContainer[ 0 ] += dimensionalAdditionalForceInBodyFrame[ 0 ];
    dimensionalGeneralForceContainer[ 1 ] += dimensionalAdditionalForceInBodyFrame[ 1 ];
    dimensionalGeneralForceContainer[ 2 ] += dimensionalAdditionalForceInBodyFrame[ 2 ];
    dimensionalGeneralForceContainer[ 3 ] += dimensionalAdditionalMomentInBodyFrame[ 0 ];
    dimensionalGeneralForceContainer[ 4 ] += dimensionalAdditionalMomentInBodyFrame[ 1 ];
    dimensionalGeneralForceContainer[ 5 ] += dimensionalAdditionalMomentInBodyFrame[ 2 ];
}

void SixDofSolver::GetDimensionalAdditionalForceAndMomentOfBody(vector< RDouble > & dimensionalGeneralForceContainer, int iTimeLayer)
{
    GetDimensionalAdditionalForceAndMomentOfBodyByParameter(dimensionalGeneralForceContainer, iTimeLayer);
}

void SixDofSolver::ReadSixDofParameterFromFile()
{
    string uDFSixDofFileName = GetStringParameterFromDataBase(bodyIndex, "uDFSixDofFileName");

    fstream infile(uDFSixDofFileName.c_str(), ios_base::in);
    
    string temp;
    int outerIterationSteps = GlobalDataBase::GetIntParaFromDB("outnstep");

    RDouble MassCenterX, MassCenterY, MassCenterZ;
    RDouble EulerAngleAlpha, EulerAngleBeta, EulerAngleGama;

    for (int iLine = 0; iLine < outerIterationSteps - 1; ++ iLine)
    {
        PHSPACE::ReadNextNonEmptyLine(infile, temp);
    }

    infile >> MassCenterX;
    infile >> MassCenterY;
    infile >> MassCenterZ;

    for (int i = 0; i < 3; ++ i)
    {
        infile >> temp;
    }

    infile >> EulerAngleAlpha;
    infile >> EulerAngleBeta;
    infile >> EulerAngleGama;

    for (int i = 0; i < 3; ++ i)
    {
        infile >> temp;
    }

    SixDofParameter * sixDofParameter = this->GetSixDofParameter(0);

    sixDofParameter->SetMassCenterX(MassCenterX);
    sixDofParameter->SetMassCenterY(MassCenterY);
    sixDofParameter->SetMassCenterZ(MassCenterZ);
    
    sixDofParameter->SetEulerAngleAlpha(EulerAngleAlpha);
    sixDofParameter->SetEulerAngleBeta (EulerAngleBeta);
    sixDofParameter->SetEulerAngleGama (EulerAngleGama);
}

void SixDofSolver::Stationary()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter(0);
    SixDofParameter * sixDofParameter1 = this->GetSixDofParameter(1);

    RDouble oldMassCenterX = sixDofParameter1->GetMassCenterX();
    RDouble oldMassCenterY = sixDofParameter1->GetMassCenterY();
    RDouble oldMassCenterZ = sixDofParameter1->GetMassCenterZ();

    RDouble oldEulerAngleAlpha = sixDofParameter1->GetEulerAngleAlpha();
    RDouble oldEulerAngleBeta  = sixDofParameter1->GetEulerAngleBeta();
    RDouble oldEulerAngleGama  = sixDofParameter1->GetEulerAngleGama();

    sixDofParameter0->SetMassCenterX(oldMassCenterX);
    sixDofParameter0->SetMassCenterY(oldMassCenterY);
    sixDofParameter0->SetMassCenterZ(oldMassCenterZ);

    sixDofParameter0->SetEulerAngleAlpha(oldEulerAngleAlpha);
    sixDofParameter0->SetEulerAngleBeta (oldEulerAngleBeta);
    sixDofParameter0->SetEulerAngleGama (oldEulerAngleGama);
}

void SixDofSolver::Rolling()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter(0);
    SixDofParameter * sixDofParameter1 = this->GetSixDofParameter(1);

    RDouble oldEulerAngleAlpha = sixDofParameter1->GetEulerAngleAlpha();
    RDouble oldEulerAngleBeta = sixDofParameter1->GetEulerAngleBeta();

    RDouble oldAngularVelocityX = sixDofParameter1->GetAngularVelocityX();
    RDouble oldAngularVelocityY = sixDofParameter1->GetAngularVelocityY();

    RDouble newEulerAngleAlpha = sixDofParameter0->GetEulerAngleAlpha();
    RDouble newEulerAngleBeta  = sixDofParameter0->GetEulerAngleBeta();
    RDouble newEulerAngleGama  = sixDofParameter0->GetEulerAngleGama();

    RDouble newAngularVelocityX = sixDofParameter0->GetAngularVelocityX();
    RDouble newAngularVelocityY = sixDofParameter0->GetAngularVelocityY();
    RDouble newAngularVelocityZ = sixDofParameter0->GetAngularVelocityZ();

    RDouble aleAmplitude = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "amplitude") * PI / 180.0;
    RDouble aleReduceFrequency = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "reduceFrequency");

    RDouble aleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("aleSimulationTime");
    RDouble dimensionalReferenceTime = GlobalDataBase::GetDoubleParaFromDB("dimensionalReferenceTime");
    aleSimulationTime /= dimensionalReferenceTime;

    RDouble gama = sin(2.0 * aleSimulationTime * aleReduceFrequency);
    gama = aleAmplitude * gama;
    RDouble angularVelocity = cos(2.0 * aleSimulationTime * aleReduceFrequency);
    angularVelocity = 2.0 * aleAmplitude * aleReduceFrequency * angularVelocity;

    newEulerAngleAlpha = oldEulerAngleAlpha;
    newEulerAngleBeta  = oldEulerAngleBeta;
    newEulerAngleGama  = gama;

    newAngularVelocityX = oldAngularVelocityX;
    newAngularVelocityY = oldAngularVelocityY;
    newAngularVelocityZ = angularVelocity;

    sixDofParameter0->SetEulerAngleAlpha(newEulerAngleAlpha);
    sixDofParameter0->SetEulerAngleBeta (newEulerAngleBeta);
    sixDofParameter0->SetEulerAngleGama (newEulerAngleGama);

    sixDofParameter0->SetAngularVelocityX(newAngularVelocityX);
    sixDofParameter0->SetAngularVelocityY(newAngularVelocityY);
    sixDofParameter0->SetAngularVelocityZ(newAngularVelocityZ);
}

void SixDofSolver::Yawing()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter(0);
    SixDofParameter * sixDofParameter1 = this->GetSixDofParameter(1);

    RDouble oldEulerAngleAlpha = sixDofParameter1->GetEulerAngleAlpha();
    RDouble oldEulerAngleGama  = sixDofParameter1->GetEulerAngleGama();

    RDouble oldAngularVelocityX = sixDofParameter1->GetAngularVelocityX();
    RDouble oldAngularVelocityZ = sixDofParameter1->GetAngularVelocityZ();

    RDouble newEulerAngleAlpha = sixDofParameter0->GetEulerAngleAlpha();
    RDouble newEulerAngleBeta  = sixDofParameter0->GetEulerAngleBeta();
    RDouble newEulerAngleGama  = sixDofParameter0->GetEulerAngleGama();

    RDouble newAngularVelocityX = sixDofParameter0->GetAngularVelocityX();
    RDouble newAngularVelocityY = sixDofParameter0->GetAngularVelocityY();
    RDouble newAngularVelocityZ = sixDofParameter0->GetAngularVelocityZ();

    RDouble aleAmplitude = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "amplitude") * PI / 180.0;
    RDouble aleReduceFrequency = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "reduceFrequency");

    RDouble aleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("aleSimulationTime");
    RDouble dimensionalReferenceTime = GlobalDataBase::GetDoubleParaFromDB("dimensionalReferenceTime");
    aleSimulationTime /= dimensionalReferenceTime;

    RDouble beta = sin(2.0 * aleSimulationTime * aleReduceFrequency);
    beta = aleAmplitude * beta;
    RDouble angularVelocity = cos(2.0 * aleSimulationTime * aleReduceFrequency);
    angularVelocity = 2.0 * aleAmplitude * aleReduceFrequency * angularVelocity;

    newEulerAngleAlpha = oldEulerAngleAlpha;
    newEulerAngleBeta  = beta;
    newEulerAngleGama  = oldEulerAngleGama;

    newAngularVelocityX = oldAngularVelocityX;
    newAngularVelocityY = angularVelocity;
    newAngularVelocityZ = oldAngularVelocityZ;

    sixDofParameter0->SetEulerAngleAlpha(newEulerAngleAlpha);
    sixDofParameter0->SetEulerAngleBeta (newEulerAngleBeta);
    sixDofParameter0->SetEulerAngleGama (newEulerAngleGama);

    sixDofParameter0->SetAngularVelocityX(newAngularVelocityX);
    sixDofParameter0->SetAngularVelocityY(newAngularVelocityY);
    sixDofParameter0->SetAngularVelocityZ(newAngularVelocityZ);
}

void SixDofSolver::Pitching()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter(0);
    SixDofParameter * sixDofParameter1 = this->GetSixDofParameter(1);

    RDouble oldEulerAngleBeta = sixDofParameter1->GetEulerAngleBeta();
    RDouble oldEulerAngleGama  = sixDofParameter1->GetEulerAngleGama();

    RDouble oldAngularVelocityY = sixDofParameter1->GetAngularVelocityY();
    RDouble oldAngularVelocityZ = sixDofParameter1->GetAngularVelocityZ();

    RDouble newEulerAngleAlpha = sixDofParameter0->GetEulerAngleAlpha();
    RDouble newEulerAngleBeta  = sixDofParameter0->GetEulerAngleBeta();
    RDouble newEulerAngleGama  = sixDofParameter0->GetEulerAngleGama();

    RDouble newAngularVelocityX = sixDofParameter0->GetAngularVelocityX();
    RDouble newAngularVelocityY = sixDofParameter0->GetAngularVelocityY();
    RDouble newAngularVelocityZ = sixDofParameter0->GetAngularVelocityZ();

    RDouble aleAmplitude = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "amplitude") * PI / 180.0;
    RDouble aleReduceFrequency = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "reduceFrequency");

    RDouble aleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("aleSimulationTime");
    RDouble dimensionalReferenceTime = GlobalDataBase::GetDoubleParaFromDB("dimensionalReferenceTime");
    aleSimulationTime /= dimensionalReferenceTime;

    RDouble alpha = sin(2.0 * aleSimulationTime * aleReduceFrequency);
    alpha = aleAmplitude * alpha;
    RDouble angularVelocity = cos(2.0 * aleSimulationTime * aleReduceFrequency);
    angularVelocity = 2.0 * aleAmplitude * aleReduceFrequency * angularVelocity;

    newEulerAngleAlpha = alpha;
    newEulerAngleBeta  = oldEulerAngleBeta;
    newEulerAngleGama  = oldEulerAngleGama;

    newAngularVelocityX = angularVelocity;
    newAngularVelocityY = oldAngularVelocityY;
    newAngularVelocityZ = oldAngularVelocityZ;

    sixDofParameter0->SetEulerAngleAlpha(newEulerAngleAlpha);
    sixDofParameter0->SetEulerAngleBeta (newEulerAngleBeta);
    sixDofParameter0->SetEulerAngleGama (newEulerAngleGama);

    sixDofParameter0->SetAngularVelocityX(newAngularVelocityX);
    sixDofParameter0->SetAngularVelocityY(newAngularVelocityY);
    sixDofParameter0->SetAngularVelocityZ(newAngularVelocityZ);
}

void SixDofSolver::RotationInXDirection()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter(0);
    SixDofParameter * sixDofParameter1 = this->GetSixDofParameter(1);

    RDouble oldEulerAngleAlpha = sixDofParameter1->GetEulerAngleAlpha();
    RDouble oldEulerAngleBeta  = sixDofParameter1->GetEulerAngleBeta();

    RDouble oldAngularVelocityX = sixDofParameter1->GetAngularVelocityX();
    RDouble oldAngularVelocityY = sixDofParameter1->GetAngularVelocityY();

    RDouble newEulerAngleAlpha = sixDofParameter0->GetEulerAngleAlpha();
    RDouble newEulerAngleBeta  = sixDofParameter0->GetEulerAngleBeta();
    RDouble newEulerAngleGama  = sixDofParameter0->GetEulerAngleGama();

    RDouble newAngularVelocityX = sixDofParameter0->GetAngularVelocityX();
    RDouble newAngularVelocityY = sixDofParameter0->GetAngularVelocityY();
    RDouble newAngularVelocityZ = sixDofParameter0->GetAngularVelocityZ();

    int aleDirection  = PHSPACE::GetIntegerParameterFromDataBase( bodyIndex, "direction" );
    RDouble aleRotateFrequency = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "rotateFrequency");
    RDouble referenceVelocity = GlobalDataBase::GetDoubleParaFromDB("referenceVelocity");
    aleRotateFrequency /=  referenceVelocity;

    RDouble aleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("aleSimulationTime");
    RDouble dimensionalReferenceTime = GlobalDataBase::GetDoubleParaFromDB("dimensionalReferenceTime");
    aleSimulationTime /= dimensionalReferenceTime;

    RDouble deltGama =  2.0 * PI * aleRotateFrequency * aleDirection;
    RDouble gama = aleSimulationTime * deltGama;

    RDouble angularVelocity = deltGama;

    newEulerAngleAlpha = oldEulerAngleAlpha;
    newEulerAngleBeta  = oldEulerAngleBeta;
    newEulerAngleGama  = gama;

    newAngularVelocityX = oldAngularVelocityX;
    newAngularVelocityY = oldAngularVelocityY;
    newAngularVelocityZ = angularVelocity;

    sixDofParameter0->SetEulerAngleAlpha(newEulerAngleAlpha);
    sixDofParameter0->SetEulerAngleBeta (newEulerAngleBeta);
    sixDofParameter0->SetEulerAngleGama (newEulerAngleGama);

    sixDofParameter0->SetAngularVelocityX(newAngularVelocityX);
    sixDofParameter0->SetAngularVelocityY(newAngularVelocityY);
    sixDofParameter0->SetAngularVelocityZ(newAngularVelocityZ);
}

void SixDofSolver::RotationInZDirection()
{
    SixDofParameter* sixDofParameter0 = this->GetSixDofParameter(0);
    SixDofParameter* sixDofParameter1 = this->GetSixDofParameter(1);

    RDouble oldEulerAngleBeta = sixDofParameter1->GetEulerAngleBeta();
    RDouble oldEulerAngleGama = sixDofParameter1->GetEulerAngleGama();

    RDouble oldAngularVelocityY = sixDofParameter1->GetAngularVelocityY();
    RDouble oldAngularVelocityZ = sixDofParameter1->GetAngularVelocityZ();

    RDouble newEulerAngleAlpha = sixDofParameter0->GetEulerAngleAlpha();
    RDouble newEulerAngleBeta = sixDofParameter0->GetEulerAngleBeta();
    RDouble newEulerAngleGama = sixDofParameter0->GetEulerAngleGama();

    RDouble newAngularVelocityX = sixDofParameter0->GetAngularVelocityX();
    RDouble newAngularVelocityY = sixDofParameter0->GetAngularVelocityY();
    RDouble newAngularVelocityZ = sixDofParameter0->GetAngularVelocityZ();

    RDouble aleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("aleSimulationTime");
    RDouble dimensionalReferenceTime = GlobalDataBase::GetDoubleParaFromDB("dimensionalReferenceTime");
    aleSimulationTime /= dimensionalReferenceTime;

    int aleDirection = PHSPACE::GetIntegerParameterFromDataBase(bodyIndex, "direction");
    RDouble aleRotateFrequency = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "rotateFrequency");
    RDouble referenceVelocity = GlobalDataBase::GetDoubleParaFromDB("referenceVelocity");
    aleRotateFrequency /=  referenceVelocity;

    RDouble deltAlpha =  2.0 * PI * aleRotateFrequency * aleDirection;
    RDouble alpha = deltAlpha * aleSimulationTime; 

    RDouble angularVelocity = deltAlpha;

    newEulerAngleAlpha = alpha;
    newEulerAngleBeta  = oldEulerAngleBeta;
    newEulerAngleGama  = oldEulerAngleGama;

    newAngularVelocityX = angularVelocity;
    newAngularVelocityY = oldAngularVelocityY;
    newAngularVelocityZ = oldAngularVelocityZ;

    sixDofParameter0->SetEulerAngleAlpha(newEulerAngleAlpha);
    sixDofParameter0->SetEulerAngleBeta(newEulerAngleBeta);
    sixDofParameter0->SetEulerAngleGama(newEulerAngleGama);

    sixDofParameter0->SetAngularVelocityX(newAngularVelocityX);
    sixDofParameter0->SetAngularVelocityY(newAngularVelocityY);
    sixDofParameter0->SetAngularVelocityZ(newAngularVelocityZ);
}

void SixDofSolver::Rotation_Rudder_Test1()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter(0);
    SixDofParameter * sixDofParameter1 = this->GetSixDofParameter(1);

    RDouble oldMassCenterX     = sixDofParameter1->GetMassCenterX();
    RDouble oldMassCenterY     = sixDofParameter1->GetMassCenterY();
    RDouble oldMassCenterZ     = sixDofParameter1->GetMassCenterZ();

    RDouble oldEulerAngleBeta  = sixDofParameter1->GetEulerAngleBeta ();
    RDouble oldEulerAngleGama  = sixDofParameter1->GetEulerAngleGama ();

    RDouble newMassCenterX     = sixDofParameter0->GetMassCenterX();
    RDouble newMassCenterY     = sixDofParameter0->GetMassCenterY();
    RDouble newMassCenterZ     = sixDofParameter0->GetMassCenterZ();

    RDouble newEulerAngleAlpha = sixDofParameter0->GetEulerAngleAlpha();
    RDouble newEulerAngleBeta  = sixDofParameter0->GetEulerAngleBeta ();
    RDouble newEulerAngleGama  = sixDofParameter0->GetEulerAngleGama ();

    RDouble aleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("aleSimulationTime");
    RDouble dimensionalReferenceTime = GlobalDataBase::GetDoubleParaFromDB("dimensionalReferenceTime");
    aleSimulationTime /= dimensionalReferenceTime;

    vector< RDouble >                   controlLawParameters;
    GetRDoubleVectorFromDataBase(       controlLawParameters, 1, "controlLawParameters", 2);

    RDouble maximumRudderAngle        = controlLawParameters[ 0 ] * PI / 180.0;
    RDouble imposedRudderRotationTime = controlLawParameters[ 1 ];

    RDouble rudderAngle = maximumRudderAngle * 0.5 * 
                        (1.0 + cos(PI * (aleSimulationTime / imposedRudderRotationTime  - 1.0)));
    
    if (aleSimulationTime > imposedRudderRotationTime)
    {
        rudderAngle = maximumRudderAngle;
    }

    cout << "rudderAngle = " << rudderAngle * RadToDeg << endl;

    newMassCenterX     = oldMassCenterX;
    newMassCenterY     = oldMassCenterY;
    newMassCenterZ     = oldMassCenterZ;

    newEulerAngleAlpha = rudderAngle;
    newEulerAngleBeta  = oldEulerAngleBeta;
    newEulerAngleGama  = oldEulerAngleGama;

    sixDofParameter0->SetMassCenterX(newMassCenterX);
    sixDofParameter0->SetMassCenterY(newMassCenterY);
    sixDofParameter0->SetMassCenterZ(newMassCenterZ);

    sixDofParameter0->SetEulerAngleAlpha(newEulerAngleAlpha);
    sixDofParameter0->SetEulerAngleBeta (newEulerAngleBeta );
    sixDofParameter0->SetEulerAngleGama (newEulerAngleGama );
}

void SixDofSolver::OscillatingInXDirection()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter(0);
    SixDofParameter * sixDofParameter1 = this->GetSixDofParameter(1);

    RDouble oldMassCenterX     = sixDofParameter1->GetMassCenterX();
    RDouble oldMassCenterY     = sixDofParameter1->GetMassCenterY();
    RDouble oldMassCenterZ     = sixDofParameter1->GetMassCenterZ();

    RDouble oldEulerAngleAlpha = sixDofParameter1->GetEulerAngleAlpha();
    RDouble oldEulerAngleBeta  = sixDofParameter1->GetEulerAngleBeta();
    RDouble oldEulerAngleGama  = sixDofParameter1->GetEulerAngleGama();

    RDouble newMassCenterX     = sixDofParameter0->GetMassCenterX();
    RDouble newMassCenterY     = sixDofParameter0->GetMassCenterY();
    RDouble newMassCenterZ     = sixDofParameter0->GetMassCenterZ();

    RDouble newEulerAngleAlpha = sixDofParameter0->GetEulerAngleAlpha();
    RDouble newEulerAngleBeta  = sixDofParameter0->GetEulerAngleBeta();
    RDouble newEulerAngleGama  = sixDofParameter0->GetEulerAngleGama();

    RDouble aleAmplitude = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "amplitude");
    RDouble alePeriod    = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "period"   );

    RDouble aleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("aleSimulationTime");
    RDouble aleTimeStep      = PHSPACE::GetAleTimeStep();
    RDouble dimensionalReferenceTime = GlobalDataBase::GetDoubleParaFromDB("dimensionalReferenceTime");
    aleSimulationTime /= dimensionalReferenceTime;
    aleTimeStep /= dimensionalReferenceTime;

    RDouble ds = aleAmplitude * (sin(2.0 * PI * aleSimulationTime / alePeriod) - sin(2.0 * PI * ((aleSimulationTime - aleTimeStep) / alePeriod)));

    //ds = -1.0 * dimensionalAleTimeStep;

    RDouble dx = ds;
    RDouble dy = 0.0;
    RDouble dz = 0.0;

    newMassCenterX     = oldMassCenterX + dx;
    newMassCenterY     = oldMassCenterY + dy;
    newMassCenterZ     = oldMassCenterZ + dz;

    newEulerAngleAlpha = oldEulerAngleAlpha;
    newEulerAngleBeta  = oldEulerAngleBeta;
    newEulerAngleGama  = oldEulerAngleGama;

    sixDofParameter0->SetMassCenterX(newMassCenterX);
    sixDofParameter0->SetMassCenterY(newMassCenterY);
    sixDofParameter0->SetMassCenterZ(newMassCenterZ);

    sixDofParameter0->SetEulerAngleAlpha(newEulerAngleAlpha);
    sixDofParameter0->SetEulerAngleBeta (newEulerAngleBeta );
    sixDofParameter0->SetEulerAngleGama (newEulerAngleGama );

    RDouble TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    TorqueRefX += dx;
    GlobalDataBase::UpdateData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
}

void SixDofSolver::OscillatingInYDirection()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter(0);
    SixDofParameter * sixDofParameter1 = this->GetSixDofParameter(1);

    RDouble oldMassCenterX = sixDofParameter1->GetMassCenterX();
    RDouble oldMassCenterY = sixDofParameter1->GetMassCenterY();
    RDouble oldMassCenterZ = sixDofParameter1->GetMassCenterZ();

    RDouble oldEulerAngleAlpha = sixDofParameter1->GetEulerAngleAlpha();
    RDouble oldEulerAngleBeta = sixDofParameter1->GetEulerAngleBeta();
    RDouble oldEulerAngleGama = sixDofParameter1->GetEulerAngleGama();

    RDouble newMassCenterX = sixDofParameter0->GetMassCenterX();
    RDouble newMassCenterY = sixDofParameter0->GetMassCenterY();
    RDouble newMassCenterZ = sixDofParameter0->GetMassCenterZ();

    RDouble newEulerAngleAlpha = sixDofParameter0->GetEulerAngleAlpha();
    RDouble newEulerAngleBeta = sixDofParameter0->GetEulerAngleBeta();
    RDouble newEulerAngleGama = sixDofParameter0->GetEulerAngleGama();

    RDouble aleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("aleSimulationTime");
    RDouble aleTimeStep = PHSPACE::GetAleTimeStep();
    RDouble dimensionalReferenceTime = GlobalDataBase::GetDoubleParaFromDB("dimensionalReferenceTime");
    aleSimulationTime /= dimensionalReferenceTime;
    aleTimeStep /= dimensionalReferenceTime;

    RDouble aleAmplitude = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "amplitude") * PI / 180.0;
    RDouble aleReduceFrequency = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "reduceFrequency");

    RDouble attackd = GlobalDataBase::GetDoubleParaFromDB("attackd")* PI / 180.0;

    RDouble alpha = sin(2.0 * aleSimulationTime * aleReduceFrequency);
    alpha = aleAmplitude * alpha;

    RDouble ds = (tan(attackd + alpha) * cos(attackd) - sin(attackd))* aleTimeStep;

    RDouble dx = 0.0;
    RDouble dy = ds;
    RDouble dz = 0.0;

    newMassCenterX = oldMassCenterX + dx;
    newMassCenterY = oldMassCenterY + dy;
    newMassCenterZ = oldMassCenterZ + dz;

    newEulerAngleAlpha = oldEulerAngleAlpha;
    newEulerAngleBeta = oldEulerAngleBeta;
    newEulerAngleGama = oldEulerAngleGama;

    sixDofParameter0->SetMassCenterX(newMassCenterX);
    sixDofParameter0->SetMassCenterY(newMassCenterY);
    sixDofParameter0->SetMassCenterZ(newMassCenterZ);

    sixDofParameter0->SetEulerAngleAlpha(newEulerAngleAlpha);
    sixDofParameter0->SetEulerAngleBeta(newEulerAngleBeta);
    sixDofParameter0->SetEulerAngleGama(newEulerAngleGama);

    RDouble TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    TorqueRefY += dy;
    GlobalDataBase::UpdateData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
}

void SixDofSolver::XDirectionOnly()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter(0);
    SixDofParameter * sixDofParameter1 = this->GetSixDofParameter(1);

    RDouble oldMassCenterX     = sixDofParameter1->GetMassCenterX();
    RDouble oldMassCenterY     = sixDofParameter1->GetMassCenterY();
    RDouble oldMassCenterZ     = sixDofParameter1->GetMassCenterZ();

    RDouble oldEulerAngleAlpha = sixDofParameter1->GetEulerAngleAlpha();
    RDouble oldEulerAngleBeta  = sixDofParameter1->GetEulerAngleBeta();
    RDouble oldEulerAngleGama  = sixDofParameter1->GetEulerAngleGama();

    RDouble newMassCenterX     = sixDofParameter0->GetMassCenterX();
    RDouble newMassCenterY     = sixDofParameter0->GetMassCenterY();
    RDouble newMassCenterZ     = sixDofParameter0->GetMassCenterZ();

    RDouble newEulerAngleAlpha = sixDofParameter0->GetEulerAngleAlpha();
    RDouble newEulerAngleBeta  = sixDofParameter0->GetEulerAngleBeta();
    RDouble newEulerAngleGama  = sixDofParameter0->GetEulerAngleGama();

    RDouble aleAmplitude = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "amplitude");

    RDouble aleTimeStep      = PHSPACE::GetAleTimeStep();
    RDouble dimensionalReferenceTime = GlobalDataBase::GetDoubleParaFromDB("dimensionalReferenceTime");
    aleTimeStep /= dimensionalReferenceTime;

    RDouble ds = aleAmplitude * aleTimeStep;

    newMassCenterX     = oldMassCenterX + ds;
    newMassCenterY     = oldMassCenterY;
    newMassCenterZ     = oldMassCenterZ;

    newEulerAngleAlpha = oldEulerAngleAlpha;
    newEulerAngleBeta  = oldEulerAngleBeta;
    newEulerAngleGama  = oldEulerAngleGama;

    sixDofParameter0->SetMassCenterX(newMassCenterX);
    sixDofParameter0->SetMassCenterY(newMassCenterY);
    sixDofParameter0->SetMassCenterZ(newMassCenterZ);

    sixDofParameter0->SetEulerAngleAlpha(newEulerAngleAlpha);
    sixDofParameter0->SetEulerAngleBeta (newEulerAngleBeta );
    sixDofParameter0->SetEulerAngleGama (newEulerAngleGama );

    RDouble TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    TorqueRefX += ds;
    GlobalDataBase::UpdateData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
}

void SixDofSolver::YDirectionOnly()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter(0);
    SixDofParameter * sixDofParameter1 = this->GetSixDofParameter(1);

    RDouble oldMassCenterX     = sixDofParameter1->GetMassCenterX();
    RDouble oldMassCenterY     = sixDofParameter1->GetMassCenterY();
    RDouble oldMassCenterZ     = sixDofParameter1->GetMassCenterZ();

    RDouble oldEulerAngleAlpha = sixDofParameter1->GetEulerAngleAlpha();
    RDouble oldEulerAngleBeta  = sixDofParameter1->GetEulerAngleBeta();
    RDouble oldEulerAngleGama  = sixDofParameter1->GetEulerAngleGama();

    RDouble newMassCenterX     = sixDofParameter0->GetMassCenterX();
    RDouble newMassCenterY     = sixDofParameter0->GetMassCenterY();
    RDouble newMassCenterZ     = sixDofParameter0->GetMassCenterZ();

    RDouble newEulerAngleAlpha = sixDofParameter0->GetEulerAngleAlpha();
    RDouble newEulerAngleBeta  = sixDofParameter0->GetEulerAngleBeta();
    RDouble newEulerAngleGama  = sixDofParameter0->GetEulerAngleGama();

    RDouble aleAmplitude = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "amplitude");

    RDouble aleTimeStep      = PHSPACE::GetAleTimeStep();
    RDouble dimensionalReferenceTime = GlobalDataBase::GetDoubleParaFromDB("dimensionalReferenceTime");
    aleTimeStep /= dimensionalReferenceTime;

    RDouble ds = aleAmplitude * aleTimeStep;

    newMassCenterX     = oldMassCenterX;
    newMassCenterY     = oldMassCenterY + ds;
    newMassCenterZ     = oldMassCenterZ;

    newEulerAngleAlpha = oldEulerAngleAlpha;
    newEulerAngleBeta  = oldEulerAngleBeta;
    newEulerAngleGama  = oldEulerAngleGama;

    sixDofParameter0->SetMassCenterX(newMassCenterX);
    sixDofParameter0->SetMassCenterY(newMassCenterY);
    sixDofParameter0->SetMassCenterZ(newMassCenterZ);

    sixDofParameter0->SetEulerAngleAlpha(newEulerAngleAlpha);
    sixDofParameter0->SetEulerAngleBeta (newEulerAngleBeta );
    sixDofParameter0->SetEulerAngleGama (newEulerAngleGama );

    RDouble TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    TorqueRefY += ds;
    GlobalDataBase::UpdateData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
}

void SixDofSolver::ModifyForceAndSixDofInformationForXYSymmetryCase()
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();

    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        vector< RDouble > & forceVector = this->GetForceVector(iTimeLayer);
        forceVector[ 2 ] = 0;
        forceVector[ 3 ] = 0;
        forceVector[ 4 ] = 0;
    }

    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        SixDofParameter * sixDofParameter = this->GetSixDofParameter(iTimeLayer);

        RDouble massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ;
        RDouble angularVelocityX, angularVelocityY, angularVelocityZ;

        sixDofParameter->GetMassCenterVelocity(massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ);
        sixDofParameter->GetAngularVelocity(angularVelocityX, angularVelocityY, angularVelocityZ);

        massCenterVelocityZ = 0.0;
        angularVelocityX    = 0.0;
        angularVelocityY    = 0.0;

        sixDofParameter->SetMassCenterVelocity(massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ);
        sixDofParameter->SetAngularVelocity(angularVelocityX, angularVelocityY, angularVelocityZ);
    }
}

void SixDofSolver::ModifyForceAndSixDofInformationForPitchingOnly()
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();

    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        vector< RDouble > & forceVector = this->GetForceVector(iTimeLayer);
        forceVector[ 0 ] = 0;
        forceVector[ 1 ] = 0;
        forceVector[ 2 ] = 0;
        forceVector[ 3 ] = 0;
        forceVector[ 4 ] = 0;
    }

    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        SixDofParameter * sixDofParameter = this->GetSixDofParameter(iTimeLayer);

        RDouble massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ;
        RDouble angularVelocityX, angularVelocityY, angularVelocityZ;

        sixDofParameter->GetMassCenterVelocity(massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ);
        sixDofParameter->GetAngularVelocity(angularVelocityX, angularVelocityY, angularVelocityZ);
        
        massCenterVelocityX = 0.0;
        massCenterVelocityY = 0.0;
        massCenterVelocityZ = 0.0;
        angularVelocityX    = 0.0;
        angularVelocityY    = 0.0;

        sixDofParameter->SetMassCenterVelocity(massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ);
        sixDofParameter->SetAngularVelocity(angularVelocityX, angularVelocityY, angularVelocityZ);
    }
}

void SixDofSolver::ModifyForceAndSixDofInformationForRollingOnly()
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();

    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        vector< RDouble > & forceVector = this->GetForceVector(iTimeLayer);
        forceVector[ 0 ] = 0;
        forceVector[ 1 ] = 0;
        forceVector[ 2 ] = 0;
        forceVector[ 4 ] = 0;
        forceVector[ 5 ] = 0;
    }

    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        SixDofParameter * sixDofParameter = this->GetSixDofParameter(iTimeLayer);

        RDouble massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ;
        RDouble angularVelocityX, angularVelocityY, angularVelocityZ;

        sixDofParameter->GetMassCenterVelocity(massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ);
        sixDofParameter->GetAngularVelocity(angularVelocityX, angularVelocityY, angularVelocityZ);
        
        massCenterVelocityX = 0.0;
        massCenterVelocityY = 0.0;
        massCenterVelocityZ = 0.0;
        angularVelocityY    = 0.0;
        angularVelocityZ    = 0.0;

        sixDofParameter->SetMassCenterVelocity(massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ);
        sixDofParameter->SetAngularVelocity(angularVelocityX, angularVelocityY, angularVelocityZ);
    }
}

void SixDofSolver::ModifyForceAndSixDofInformationForPitchingAndRollingOnly()
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();

    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        vector< RDouble > & forceVector = this->GetForceVector(iTimeLayer);
        forceVector[ 0 ] = 0;
        forceVector[ 1 ] = 0;
        forceVector[ 2 ] = 0;
        forceVector[ 4 ] = 0;
    }

    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        SixDofParameter * sixDofParameter = this->GetSixDofParameter(iTimeLayer);

        RDouble massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ;
        RDouble angularVelocityX, angularVelocityY, angularVelocityZ;

        sixDofParameter->GetMassCenterVelocity(massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ);
        sixDofParameter->GetAngularVelocity(angularVelocityX, angularVelocityY, angularVelocityZ);
        
        massCenterVelocityX = 0.0;
        massCenterVelocityY = 0.0;
        massCenterVelocityZ = 0.0;
        angularVelocityY    = 0.0;

        sixDofParameter->SetMassCenterVelocity(massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ);
        sixDofParameter->SetAngularVelocity(angularVelocityX, angularVelocityY, angularVelocityZ);
    }
}

void SixDofSolver::ComputeCurrentSixDofVariables(vector< vector< RDouble > > & qContainer, vector< vector< RDouble > > & fluxContainer, int numberOfEquations)
{
    this->ComputeCurrentSixDofVariables(PHSPACE::HeaderPointer(qContainer   [ 0 ]), 
                                         PHSPACE::HeaderPointer(qContainer   [ 1 ]),
                                         PHSPACE::HeaderPointer(qContainer   [ 2 ]),
                                         PHSPACE::HeaderPointer(fluxContainer[ 0 ]),
                                         PHSPACE::HeaderPointer(fluxContainer[ 1 ]), 
                                         PHSPACE::HeaderPointer(fluxContainer[ 2 ]),
                                         numberOfEquations);
}

void SixDofSolver::ComputeCurrentSixDofVariables(RDouble * q, RDouble * qN, RDouble * qN1, RDouble * fluxM, RDouble * fluxN, RDouble * fluxN1, int numberOfEquations)
{
    RDouble a, b, c, d, e;
    RDouble relaxParameter = GlobalDataBase::GetDoubleParaFromDB("relaxParameterOfKinetic");
    PHSPACE::GetSixDofCoefficients(a, b, c, d, e, relaxParameter);

    RDouble coef1 = (1 - relaxParameter);
    RDouble coef2 = relaxParameter * a;
    RDouble coef3 = relaxParameter * b;
    RDouble coef4 = relaxParameter * c;
    RDouble coef5 = relaxParameter * d;
    RDouble coef6 = relaxParameter * e;

    for (int iEquation = 0; iEquation < numberOfEquations; ++ iEquation)
    {
        q[ iEquation ] = coef1 * q     [ iEquation ] + 
                         coef2 * qN    [ iEquation ] + 
                         coef3 * qN1   [ iEquation ] + 
                         coef4 * fluxM [ iEquation ] + 
                         coef5 * fluxN [ iEquation ] + 
                         coef6 * fluxN1[ iEquation ];
    }
}

void SixDofSolver::SolveSixDofEquation()
{
    const int numberOfEquations = 12;
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();

    vector< vector< RDouble > > sixDofVariableContainer;
    vector< vector< RDouble > > sixDofFluxContainer;

    AllocateVector(sixDofVariableContainer, numberOfTimeLayers, numberOfEquations);
    AllocateVector(sixDofFluxContainer    , numberOfTimeLayers, numberOfEquations);

    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        vector< RDouble > & forceVector     = this->GetForceVector    (iTimeLayer);
        SixDofParameter * sixDofParameter = this->GetSixDofParameter(iTimeLayer);

        sixDofParameter->GetGeneralVariableContainer(sixDofVariableContainer[ iTimeLayer ]);

        GetSixDofFlux(sixDofParameter, massCharacter, forceVector, sixDofFluxContainer[ iTimeLayer ]);
    }

    ComputeCurrentSixDofVariables(sixDofVariableContainer, sixDofFluxContainer, numberOfEquations);

    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter(0);

    sixDofParameter0->SetGeneralVariableContainer(& sixDofVariableContainer[ 0 ][ 0 ]);
}

void SixDofSolver::MovingByFile()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter(0);

    int     outIterationStep = GlobalDataBase::GetIntParaFromDB("outnstep");

    int aleStartStrategy = GlobalDataBase::GetIntParaFromDB("aleStartStrategy");
    if (aleStartStrategy == -1)
    {
        outIterationStep = 1;
    }

    RDouble newMassCenterX = 0.0;
    RDouble newMassCenterY = 0.0;
    RDouble newMassCenterZ = 0.0;

    RDouble newEulerAngleGama  = 0.0;
    RDouble newEulerAngleBeta  = 0.0;
    RDouble newEulerAngleAlpha = 0.0;

    RDouble tmp;

    string uDFSixDofFileName = GetStringParameterFromDataBase(bodyIndex, "uDFSixDofFileName");

    ios_base::openmode openMode = ios_base::in;
    fstream file;

    PHSPACE::OpenFile(file, uDFSixDofFileName, openMode);

    for (int iStep = 0; iStep < outIterationStep; ++ iStep)
    {
        file >> tmp;
        file >> newMassCenterX;
        file >> newMassCenterY;
        file >> newMassCenterZ;
        file >> newEulerAngleGama ;
        file >> newEulerAngleBeta ;
        file >> newEulerAngleAlpha;
    }

    newEulerAngleGama  *= PI / 180.0;
    newEulerAngleBeta  *= PI / 180.0;
    newEulerAngleAlpha *= PI / 180.0;

    sixDofParameter0->SetMassCenterX(newMassCenterX);
    sixDofParameter0->SetMassCenterY(newMassCenterY);
    sixDofParameter0->SetMassCenterZ(newMassCenterZ);

    sixDofParameter0->SetEulerAngleAlpha(newEulerAngleAlpha);
    sixDofParameter0->SetEulerAngleBeta (newEulerAngleBeta );
    sixDofParameter0->SetEulerAngleGama (newEulerAngleGama );
}

void SixDofSolver::GetDimensionalAerodynamicForceAndMoment(vector< RDouble > & dimensionalAerodynamicForce, vector< RDouble > & dimensionalAerodynamicMoment, int bodyIndex, int iTimeLayer)
{
    AleForceManager * aleForceManager = GetAleManager()->GetAleForceManager();

    BasicAleForce * basicAleForce = aleForceManager->GetBasicAleForce(bodyIndex, iTimeLayer);

    RDouble dimensionalForceCoefficient  = aleForceManager->GetDimensionalForceCoefficient ();
    RDouble dimensionalMomentCoefficient = aleForceManager->GetDimensionalMomentCoefficient();
    dimensionalAerodynamicForce [0] = basicAleForce->GetForceX() * dimensionalForceCoefficient;
    dimensionalAerodynamicForce [1] = basicAleForce->GetForceY() * dimensionalForceCoefficient;
    dimensionalAerodynamicForce [2] = basicAleForce->GetForceZ() * dimensionalForceCoefficient;

    dimensionalAerodynamicMoment[0] = basicAleForce->GetMomentX() * dimensionalMomentCoefficient;
    dimensionalAerodynamicMoment[1] = basicAleForce->GetMomentY() * dimensionalMomentCoefficient;
    dimensionalAerodynamicMoment[2] = basicAleForce->GetMomentZ() * dimensionalMomentCoefficient;
}

void SixDofSolver::DumpSixDofSolverInformation()
{
    this->DumpSixDofParameterInformation();
    this->DumpAleForceInformation();
}

void SixDofSolver::DumpSixDofParameterInformation()
{
    string fileName = GlobalDataBase::GetStrParaFromDB("sixDofFileName");
    string sixDofFileName = PHSPACE::AddSymbolToFileName(fileName, bodyIndex + 1);

    ios_base::openmode openMode = ios_base::out|ios_base::app;

    fstream file;
    PHSPACE::OpenFile(file, sixDofFileName, openMode);

    this->DumpSixDofParameterInformation(file);

    PHSPACE::CloseFile(file);
}

void SixDofSolver::DumpAleForceInformation()
{
    string fileName = GlobalDataBase::GetStrParaFromDB("sixDofFileName");
    string aleForceFileName = PHSPACE::AddSymbolToFileName(fileName, bodyIndex + 1);

    ios_base::openmode openMode = ios_base::out|ios_base::app;
    fstream file;

    string fileName1 = PHSPACE::AddSymbolToFileName(aleForceFileName, "_body");
    PHSPACE::OpenFile(file, fileName1, openMode);
    this->DumpAleForceInformationBodyFrame(file);
    PHSPACE::CloseFile(file);

    string fileName2 = PHSPACE::AddSymbolToFileName(aleForceFileName, "_iner");
    PHSPACE::OpenFile(file, fileName2, openMode);
    this->DumpAleForceInformationInerFrame(file);
    PHSPACE::CloseFile(file);
}

void SixDofSolver::DumpSixDofParameterInformation(fstream & file)
{
    RDouble massCenterX, massCenterY, massCenterZ;
    RDouble massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ;
    RDouble eulerAngleAlpha, eulerAngleBeta, eulerAngleGama;
    RDouble angularVelocityX, angularVelocityY, angularVelocityZ;

    int iTimeLayer = 0;
    SixDofParameter * sixDofParameter =  this->GetSixDofParameter(iTimeLayer);

    sixDofParameter->GetMassCenter(massCenterX, massCenterY, massCenterZ);
    sixDofParameter->GetMassCenterVelocity(massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ);
    sixDofParameter->GetEulerAngle(eulerAngleGama, eulerAngleBeta,eulerAngleAlpha);
    sixDofParameter->GetAngularVelocity(angularVelocityX, angularVelocityY, angularVelocityZ);

    file << setiosflags(ios::left);
    file << setiosflags(ios::scientific);
    file << setprecision(6);

    RDouble aleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("aleSimulationTime");

    file << aleSimulationTime << "\t";

    file << massCenterX << "\t";
    file << massCenterY << "\t";
    file << massCenterZ << "\t";
    file << massCenterVelocityX << "\t";
    file << massCenterVelocityY << "\t";
    file << massCenterVelocityZ << "\t";
    file << eulerAngleGama << "\t";
    file << eulerAngleBeta << "\t";
    file << eulerAngleAlpha << "\t";
    file << angularVelocityX << "\t";
    file << angularVelocityY << "\t";
    file << angularVelocityZ << "\t";

    file << "\n";
}

void SixDofSolver::DumpAleForceInformationBodyFrame(fstream & file)
{
    int numberOfAleForceComponents = 6;
    vector< RDouble > aleForceContainer(numberOfAleForceComponents);

    this->GetDimensionalAerodynamicForceAndMomentBodyFrame(aleForceContainer, 0);

    file << setiosflags(ios::left);
    file << setiosflags(ios::scientific);
    file << setprecision(6);

    RDouble aleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("aleSimulationTime");

    file << aleSimulationTime << "\t";

    for (int iAleForceComponent = 0; iAleForceComponent < numberOfAleForceComponents; ++ iAleForceComponent)
    {
        file << aleForceContainer[ iAleForceComponent ] << "\t";
    }
    file << "\n";
}

void SixDofSolver::DumpAleForceInformationInerFrame(fstream & file)
{
    int numberOfAleForceComponents = 6;
    vector< RDouble > aleForceContainer(numberOfAleForceComponents);

    this->GetDimensionalAerodynamicForceAndMomentInerFrame(aleForceContainer, 0);

    file << setiosflags(ios::left);
    file << setiosflags(ios::scientific);
    file << setprecision(6);

    RDouble aleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("aleSimulationTime");

    file << aleSimulationTime << "\t";

    for (int iAleForceComponent = 0; iAleForceComponent < numberOfAleForceComponents; ++ iAleForceComponent)
    {
        file << aleForceContainer[ iAleForceComponent ] << "\t";
    }
    file << "\n";
}

void SixDofSolver::DumpRestartSixDofSolver(fstream & file)
{
    const int numberOfEquations = 12;
    vector< RDouble >  generalVariableContainer(numberOfEquations);

    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer)
    {
        SixDofParameter * sixDofParameter = this->GetSixDofParameter(iTimeLayer);
        sixDofParameter->GetGeneralVariableContainer(generalVariableContainer);
        if (resetMassCenter == 1)
        {
        }
        for (int iVar = 0; iVar < numberOfEquations; ++ iVar)
        {
            file << generalVariableContainer[iVar] << "   ";
        }
        file << "\n";
        //PHSPACE::PHWrite(file, generalVariableContainer, numberOfEquations);
    }
}

void GetSixDofCoefficients(RDouble & a, RDouble & b, RDouble & c, RDouble & d, RDouble & e, RDouble & relaxParameter)
{
    int methodForKineticEquation = PHSPACE::GlobalDataBase::GetIntParaFromDB("methodForKineticEquation");
    relaxParameter = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("relaxParameterOfKinetic");

    RDouble aleTimeStep  = PHSPACE::GetAleTimeStep();

    RDouble thet = 0.0, xi = 0.0, phi = 0.0;

    if (     methodForKineticEquation == KINETIC_SCHEME::ADMAS_BASHFORTH_FIRST_ORDER )
    {
        thet = 0.0;
        xi   = 0.0;
        phi  = 0.0;
        relaxParameter = 1.0;
    }
    else if (methodForKineticEquation == KINETIC_SCHEME::ADMAS_BASHFORTH_SECOND_ORDER)
    {
        thet = 0.0;
        xi   = 0.0;
        phi  = 1.0 / 2.0;
        relaxParameter = 1.0;
    }
    else if (methodForKineticEquation == KINETIC_SCHEME::IMPLICIT_EULER_FIRST_ORDER )
    {
        thet = 1.0;
        xi   = 0.0;
        phi  = 0.0;
    }
    else if (methodForKineticEquation == KINETIC_SCHEME::IMPLICIT_EULER_SECOND_ORDER )
    {
        thet = 1.0;
        xi   = 1.0 / 2.0;
        phi  = 0.0;
    }
    else if (methodForKineticEquation == KINETIC_SCHEME::ADMAS_MOULTON_SECOND_ORDER )
    {
        thet = 1.0 / 2.0;
        xi   = 0.0;
        phi  = 0.0;
    }
    else if (methodForKineticEquation == KINETIC_SCHEME::ADMAS_MOULTON_THIRD_ORDER)
    {
        thet = 5.0 / 12.0;
        xi   = 0.0;
        phi  = 1.0 / 12.0;
    }

    a = aleTimeStep / (1.0 + xi) * ((1.0 + xi) / aleTimeStep + xi / aleTimeStep);
    b = aleTimeStep / (1.0 + xi) * (- xi / aleTimeStep);
    c = aleTimeStep / (1.0 + xi) * (thet);
    d = aleTimeStep / (1.0 + xi) * (1.0 - thet + phi);
    e = aleTimeStep / (1.0 + xi) * (- phi);
}

void SixDofSolver::GetSixDofFlux(SixDofParameter * sixDofParameter, MassCharacter * massCharacter, vector< RDouble > & force, vector< RDouble > & generalSixDofFlux)
{
    this->ComputeRightHandSideOfRigidBodyTranslationalMotionEquation  (sixDofParameter, massCharacter, force, generalSixDofFlux);
    this->ComputeRightHandSideOfRigidBodyTranslationalDynamicsEquation(sixDofParameter, massCharacter, force, generalSixDofFlux);
    this->ComputeRightHandSideOfRigidBodyRotationalMotionEquation     (sixDofParameter, massCharacter, force, generalSixDofFlux);
    this->ComputeRightHandSideOfRigidBodyRotationalDynamicsEquation   (sixDofParameter, massCharacter, force, generalSixDofFlux);
}

void SixDofSolver::ComputeRightHandSideOfRigidBodyTranslationalMotionEquation(SixDofParameter * sixDofParameter, MassCharacter * massCharacter, vector< RDouble > & force, vector< RDouble > & generalSixDofFlux)
{
    RDouble massCenterVelocityX,  massCenterVelocityY, massCenterVelocityZ;
    sixDofParameter->GetMassCenterVelocity(massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ);

    //! location
    generalSixDofFlux[ 0 ] = massCenterVelocityX;
    generalSixDofFlux[ 1 ] = massCenterVelocityY;
    generalSixDofFlux[ 2 ] = massCenterVelocityZ;
}

void SixDofSolver::ComputeRightHandSideOfRigidBodyTranslationalDynamicsEquation(SixDofParameter * sixDofParameter, MassCharacter * massCharacter, vector< RDouble > & force, vector< RDouble > & generalSixDofFlux)
{
    RDouble forceX = force[ 0 ];
    RDouble forceY = force[ 1 ];
    RDouble forcez = force[ 2 ];

    RDouble mass = massCharacter->GetMass();

    //! velocity
    generalSixDofFlux[ 3 ] = forceX / mass;
    generalSixDofFlux[ 4 ] = forceY / mass;
    generalSixDofFlux[ 5 ] = forcez / mass;
}

void SixDofSolver::ComputeRightHandSideOfRigidBodyRotationalMotionEquation(SixDofParameter * sixDofParameter, MassCharacter * massCharacter, vector< RDouble > & force, vector< RDouble > & generalSixDofFlux)
{
    RDouble angularVelocityX, angularVelocityY, angularVelocityZ;
    RDouble eulerAngleAlpha, eulerAngleBeta, eulerAngleGama;

    sixDofParameter->GetEulerAngle(eulerAngleGama, eulerAngleBeta, eulerAngleAlpha);
    sixDofParameter->GetAngularVelocity(angularVelocityX, angularVelocityY, angularVelocityZ);

    //! attitude
    generalSixDofFlux[ 8 ] = angularVelocityY * sin(eulerAngleGama )                          + angularVelocityZ * cos(eulerAngleGama );
    generalSixDofFlux[ 7 ] = angularVelocityY * cos(eulerAngleGama ) / cos(eulerAngleAlpha) - angularVelocityZ * sin(eulerAngleGama ) / cos(eulerAngleAlpha);
    generalSixDofFlux[ 6 ] = angularVelocityX - generalSixDofFlux[ 7 ] * sin(eulerAngleAlpha);
}

void SixDofSolver::ComputeRightHandSideOfRigidBodyRotationalDynamicsEquation(SixDofParameter * sixDofParameter, MassCharacter * massCharacter, vector< RDouble > & force, vector< RDouble > & generalSixDofFlux)
{
    RDouble ** momentOfInertiaTensor        = massCharacter->GetMomentOfInertiaTensor();
    RDouble ** inverseMomentOfInertiaTensor = massCharacter->GetInverseMomentOfInertiaTensor();

    vector< RDouble > angularVelocityVector(3);
    sixDofParameter->GetAngularVelocityContainer(angularVelocityVector);

    vector< RDouble > momentumMoment(3);
    PHSPACE::MatrixMultiply(momentOfInertiaTensor, angularVelocityVector, momentumMoment, 3, 3);

    RDouble momentVector[ 3 ];
    momentVector[ 0 ] = force[ 3 ];
    momentVector[ 1 ] = force[ 4 ];
    momentVector[ 2 ] = force[ 5 ];

    vector< RDouble > temp1(3);
    PHSPACE::CrossProduct(angularVelocityVector, momentumMoment, temp1);

    vector< RDouble > temp2(3);
    temp2[ 0 ] = momentVector[ 0 ] - temp1[ 0 ];
    temp2[ 1 ] = momentVector[ 1 ] - temp1[ 1 ];
    temp2[ 2 ] = momentVector[ 2 ] - temp1[ 2 ];

    vector< RDouble > dRotationDtVector(3);

    PHSPACE::MatrixMultiply(inverseMomentOfInertiaTensor, temp2, dRotationDtVector, 3, 3);

    generalSixDofFlux[ 9  ] = dRotationDtVector[ 0 ];
    generalSixDofFlux[ 10 ] = dRotationDtVector[ 1 ];
    generalSixDofFlux[ 11 ] = dRotationDtVector[ 2 ];

    PHSPACE::DelPointer2(momentOfInertiaTensor       );
    PHSPACE::DelPointer2(inverseMomentOfInertiaTensor);
}

SixDofSolverManager::SixDofSolverManager()
{
    sixDofSolverContainer = 0;
    numberOfMovingBodies  = 0;
}

SixDofSolverManager::~SixDofSolverManager() 
{ 
    if (sixDofSolverContainer)
    {
        for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
        {
            delete sixDofSolverContainer[ iMovingBody ];
            sixDofSolverContainer[iMovingBody] = nullptr;
        }
    }
    delete [] sixDofSolverContainer;
    sixDofSolverContainer = nullptr;
}

void SixDofSolverManager::AllocateSixDofSolverContainer()
{
    numberOfMovingBodies = PHSPACE::GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    sixDofSolverContainer = new SixDofSolver * [ numberOfMovingBodies ];
    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        int RBDMethod = PHSPACE::GetIntegerParameterFromDataBase(iMovingBody, "RBDMethod");
        
        if (RBDMethod <= 100)
        {
            sixDofSolverContainer[ iMovingBody ] = new SixDofSolver(this);
        }
        else
        {
            sixDofSolverContainer[ iMovingBody ] = new SixDofSolver(this);
        }
    }
}

void SixDofSolverManager::Initialize()
{
    AllocateSixDofSolverContainer();

    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        sixDofSolverContainer[ iMovingBody ]->Initialize(iMovingBody);
    }
}

void SixDofSolverManager::Restart()
{
    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        sixDofSolverContainer[ iMovingBody ]->Restart();
    }
}

void SixDofSolverManager::DumpRestart(fstream & file)
{
    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        sixDofSolverContainer[ iMovingBody ]->DumpRestartSixDofSolver(file);
    }
}

void SixDofSolverManager::ResetMassCenter()
{
    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        sixDofSolverContainer[iMovingBody]->ResetMassCenter();
    }
}

void SixDofSolverManager::Run()
{
    using namespace PHMPI;
    int currentProcessor = GetCurrentProcessorID();
    int serverProcessor = GetServerProcessorID();

    vector< RDouble > sixDofVariableContainer;

    sixDofVariableContainer.resize(12);

    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        SixDofParameter * sixDofParameter   = sixDofSolverContainer[ iMovingBody ]->GetSixDofParameter(0);

        if (currentProcessor == serverProcessor)
        {
            sixDofSolverContainer[ iMovingBody ]->Run();
            sixDofParameter->GetGeneralVariableContainer(sixDofVariableContainer);
        }

        PHMPI::PH_Bcast(&sixDofVariableContainer[ 0 ], 12*sizeof(RDouble), serverProcessor);

        sixDofParameter->SetGeneralVariableContainer(&sixDofVariableContainer[ 0 ]);
    }
}

void SixDofSolverManager::Post()
{
    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        SixDofSolver * sixDofSolver = this->GetSixDofSolver(iMovingBody);
        sixDofSolver->UpdateUnsteadyTerm();
    }
    
    Dump();
}

void SixDofSolverManager::Dump()
{
    using namespace PHMPI;
    int currentProcessor = GetCurrentProcessorID();
    int serverProcessor = GetServerProcessorID();
    if (currentProcessor != serverProcessor)
    {
        return;
    }

    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        SixDofSolver * sixDofSolver = this->GetSixDofSolver(iMovingBody);
        sixDofSolver->DumpSixDofSolverInformation();
    }
}

MassCharacter * SixDofSolverManager::GetMassCharacter(int bodyIndex)
{
    SixDofSolver * sixDofSolver = this->GetSixDofSolver(bodyIndex);
    return sixDofSolver->GetMassCharacter();
}

SixDofParameter * SixDofSolverManager::GetSixDofParameterN1(int bodyIndex)
{
    SixDofSolver * sixDofSolver = this->GetSixDofSolver(bodyIndex);
    return sixDofSolver->GetSixDofParameter(2);
}

SixDofParameter * SixDofSolverManager::GetSixDofParameterN (int bodyIndex)
{
    SixDofSolver * sixDofSolver = this->GetSixDofSolver(bodyIndex);
    return sixDofSolver->GetSixDofParameter(1);
}

SixDofParameter * SixDofSolverManager::GetSixDofParameter  (int bodyIndex)
{
    SixDofSolver * sixDofSolver = this->GetSixDofSolver(bodyIndex);
    return sixDofSolver->GetSixDofParameter(0);
}

SixDofParameter * GetSixDofParameter(int bodyIndex)
{
    SixDofSolverManager * sixDofSolverManager = GetAleManager()->GetSixDofSolverManager();

    return sixDofSolverManager->GetSixDofParameter(bodyIndex);
}

void ConvertVectorFromInertialCoordinateSystemToBodyCoordinateSystem(vector< RDouble > & vectorInInertialFrame, vector< RDouble > & angleVector, vector< RDouble > & massCenter, vector< RDouble > & vectorInBodyFrame)
{
    vectorInBodyFrame[ 0 ] = vectorInInertialFrame[ 0 ] - massCenter[ 0 ];
    vectorInBodyFrame[ 1 ] = vectorInInertialFrame[ 1 ] - massCenter[ 1 ];
    vectorInBodyFrame[ 2 ] = vectorInInertialFrame[ 2 ] - massCenter[ 2 ];

    PHSPACE::RotateAboutAxisY(vectorInBodyFrame, angleVector[ 1 ], vectorInBodyFrame);
    PHSPACE::RotateAboutAxisZ(vectorInBodyFrame, angleVector[ 2 ], vectorInBodyFrame);
    PHSPACE::RotateAboutAxisX(vectorInBodyFrame, angleVector[ 0 ], vectorInBodyFrame);
}

void ConvertVectorFromInertialCoordinateSystemToBodyCoordinateSystem(vector< RDouble > & vectorInInertialFrame, vector< RDouble > & angleVector, vector< RDouble > & vectorInBodyFrame)
{
    PHSPACE::RotateAboutAxisY(vectorInInertialFrame, angleVector[ 1 ], vectorInBodyFrame);
    PHSPACE::RotateAboutAxisZ(vectorInBodyFrame    , angleVector[ 2 ], vectorInBodyFrame);
    PHSPACE::RotateAboutAxisX(vectorInBodyFrame    , angleVector[ 0 ], vectorInBodyFrame);
}

void ConvertVectorFromBodyCoordinateSystemToInertialCoordinateSystem(vector< RDouble > & vectorInBodyFrame, vector< RDouble > & angleVector, vector< RDouble > & massCenter,  vector< RDouble > & vectorInInertialFrame)
{
    PHSPACE::RotateAboutAxisX(vectorInBodyFrame    , - angleVector[ 0 ], vectorInInertialFrame);
    PHSPACE::RotateAboutAxisZ(vectorInInertialFrame, - angleVector[ 2 ], vectorInInertialFrame);
    PHSPACE::RotateAboutAxisY(vectorInInertialFrame, - angleVector[ 1 ], vectorInInertialFrame);

    vectorInInertialFrame[ 0 ] += massCenter[ 0 ];
    vectorInInertialFrame[ 1 ] += massCenter[ 1 ];
    vectorInInertialFrame[ 2 ] += massCenter[ 2 ];
}

void ConvertVectorFromBodyCoordinateSystemToInertialCoordinateSystem(vector< RDouble > & vectorInBodyFrame, vector< RDouble > & angleVector, vector< RDouble > & vectorInInertialFrame)
{
    PHSPACE::RotateAboutAxisX(vectorInBodyFrame    , - angleVector[ 0 ], vectorInInertialFrame);
    PHSPACE::RotateAboutAxisZ(vectorInInertialFrame, - angleVector[ 2 ], vectorInInertialFrame);
    PHSPACE::RotateAboutAxisY(vectorInInertialFrame, - angleVector[ 1 ], vectorInInertialFrame);
}

void RotateAboutAxisX(vector< RDouble > & vectorIn, RDouble rotateAngle, vector< RDouble > & vectorOut)
{
    RDouble oldX = vectorIn[ 0 ];
    RDouble oldY = vectorIn[ 1 ];
    RDouble oldZ = vectorIn[ 2 ];
    vectorOut[ 0 ] = 1.0      * oldX + 0.0                * oldY + 0.0                * oldZ;
    vectorOut[ 1 ] = 0.0      * oldX + cos(rotateAngle) * oldY + sin(rotateAngle) * oldZ;
    vectorOut[ 2 ] = 0.0      * oldX - sin(rotateAngle) * oldY + cos(rotateAngle) * oldZ;
}

void RotateAboutAxisY(vector< RDouble > & vectorIn, RDouble rotateAngle, vector< RDouble > & vectorOut)
{
    RDouble oldX = vectorIn[ 0 ];
    RDouble oldY = vectorIn[ 1 ];
    RDouble oldZ = vectorIn[ 2 ];
    vectorOut[ 0 ] = cos(rotateAngle) * oldX + 0.0      * oldY - sin(rotateAngle) * oldZ;
    vectorOut[ 1 ] = 0.0                * oldX + 1.0      * oldY + 0.0                * oldZ;
    vectorOut[ 2 ] = sin(rotateAngle) * oldX + 0.0      * oldY + cos(rotateAngle) * oldZ;
}

void RotateAboutAxisZ(vector< RDouble > & vectorIn, RDouble rotateAngle, vector< RDouble > & vectorOut)
{
    RDouble oldX = vectorIn[ 0 ];
    RDouble oldY = vectorIn[ 1 ];
    RDouble oldZ = vectorIn[ 2 ];
    vectorOut[ 0 ] =   cos(rotateAngle) * oldX + sin(rotateAngle) * oldY + 0.0      * oldZ;
    vectorOut[ 1 ] = - sin(rotateAngle) * oldX + cos(rotateAngle) * oldY + 0.0      * oldZ;
    vectorOut[ 2 ] =   0.0                * oldX + 0.0                * oldY + 1.0      * oldZ;
}

}