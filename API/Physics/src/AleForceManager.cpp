#include "AleForceManager.h"
#include "Math_BasisFunction.h"
#include "Constants.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "Geo_SimpleBC.h"
#include "TK_Time.h"
#include "TK_Parse.h"
#include "TK_Exit.h"
#include "Glb_Dimension.h"
#include "GlobalDataBase.h"
#include "PHHeader.h"
#include "PHMpi.h"
#include "AleManager.h"
#include "Solver.h"

using namespace std;

namespace PHSPACE
{
bool AerodynamicForceBody::initialize = false;
vector< std::set< string > >  AerodynamicForceBody::bodies;

AerodynamicForceBody::AerodynamicForceBody()
{
}

AerodynamicForceBody::~AerodynamicForceBody()
{
}

void AerodynamicForceBody::Initialize()
{
    if (!AerodynamicForceBody::initialize)
    {
        int numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");
        AerodynamicForceBody::bodies.resize(numberOfMovingBodies);

        set < string > *bodyList = GlobalBoundaryCondition::GetBodyNameList();
        set < string >::iterator bodyListIter = bodyList->begin();
        for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
        {
            string bodyName = *bodyListIter;
            AerodynamicForceBody::bodies[iMovingBody].insert(bodyName);
            bodyListIter ++;
        }
        AerodynamicForceBody::initialize = true;
    }
}

bool AerodynamicForceBody::FaceInAerodynamicForceBodies(string bodyName, int aerodynamicForceBodyIndex)
{
    AerodynamicForceBody::Initialize();
    set< string >::iterator iter = AerodynamicForceBody::bodies[aerodynamicForceBodyIndex].find(bodyName);
    return iter != AerodynamicForceBody::bodies[aerodynamicForceBodyIndex].end();
}

bool FaceInAerodynamicForceBodies(string bodyName, int aerodynamicForceBodyIndex)
{
    return AerodynamicForceBody::FaceInAerodynamicForceBodies(bodyName, aerodynamicForceBodyIndex);
}

AerodynamicForceManager aerodynamicForceManager;

void InitializeAerodynamicForce()
{
    AerodynamicForceManager *aerodynamicForceManager1 = PHSPACE::GetAerodynamicForceManager();
    aerodynamicForceManager1->InitializeAerodynamicForceManager();
}

void ComputeAerodynamicForceForALE()
{
    ActionKey *actkeyComputeAirForceCoef = new ActionKey();
    actkeyComputeAirForceCoef->action = COMPUTE_AERODYNAMIC_FORCE;

    PHSPACE::InitializeAerodynamicForce();

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    actkeyComputeAirForceCoef->GetData()->MoveToBegin();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int ZoneProcessorID = GetZoneProcessorID(iZone);
        int currentProcessorID = GetCurrentProcessorID();

        if (currentProcessorID == ZoneProcessorID)
        {
            //! If zone i belongs to the current process.
            PHSolver *solver = static_cast <PHSolver *> (GlobalSolvers::GetSolver(iZone, actkeyComputeAirForceCoef->solverID));
            solver->AirForceCoefBodies(actkeyComputeAirForceCoef);
            PHSPACE::AddAerodynamicForce(actkeyComputeAirForceCoef);
        }
    }

    //! Third: Collecting AirForceCoef on Each Zone to serve.
    PHSPACE::CollectAerodynamicForceOfAllProcessors();

    delete actkeyComputeAirForceCoef;
}

void AddAerodynamicForce(ActionKey *actionKey)
{
    AerodynamicForceManager *aerodynamicForceManager2 = PHSPACE::GetAerodynamicForceManager();
    aerodynamicForceManager2->AddAerodynamicForce(actionKey);
}

void CollectAerodynamicForceOfAllProcessors()
{
    AerodynamicForceManager *aerodynamicForceManager3 = PHSPACE::GetAerodynamicForceManager();
    aerodynamicForceManager3->CollectAerodynamicForceOfAllProcessors();
}

BasicAerodynamicForce **CreateBasicAerodynamicForceContainer()
{
    int numberOfAerodynamicForceBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    BasicAerodynamicForce **basicAerodynamicForceContainer = new BasicAerodynamicForce * [numberOfAerodynamicForceBodies];

    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodies; ++ iMovingBody)
    {
        basicAerodynamicForceContainer[iMovingBody] = new BasicAerodynamicForce();
    }

    return basicAerodynamicForceContainer;
}

vector< BasicAerodynamicForce * > *CreateBasicAerodynamicForceVector()
{
    int numberOfAerodynamicForceBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    vector< BasicAerodynamicForce * > *basicAerodynamicForceVector = new vector< BasicAerodynamicForce * >(numberOfAerodynamicForceBodies);

    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodies; ++ iMovingBody)
    {
        (*basicAerodynamicForceVector)[iMovingBody] = new BasicAerodynamicForce();
    }

    return basicAerodynamicForceVector;
}

void FreeBasicAerodynamicForceVector(vector<BasicAerodynamicForce *> *basicAerodynamicForceVector)
{
    int numberOfAerodynamicForceBodies = static_cast<int>(basicAerodynamicForceVector->size());

    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodies; ++ iMovingBody)
    {
        delete (*basicAerodynamicForceVector)[iMovingBody];
    }

    delete  basicAerodynamicForceVector;
}

void FreeBasicAerodynamicForceContainer(BasicAerodynamicForce **basicAerodynamicForceContainer)
{
    int numberOfAerodynamicForceBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodies; ++ iMovingBody)
    {
        delete basicAerodynamicForceContainer[iMovingBody];
    }

    delete [] basicAerodynamicForceContainer;
}

void WriteAerodynamicForceToActionKey(ActionKey *actionKey, BasicAerodynamicForce **basicAerodynamicForceContainer)
{
    int numberOfAerodynamicForceBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    DataContainer *dataContainer = actionKey->GetData();
    dataContainer->MoveToBegin();

    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodies; ++ iMovingBody)
    {
        basicAerodynamicForceContainer[iMovingBody]->WriteAerodynamicForceToDataContainer(dataContainer);
    }
}

void WriteAerodynamicForceToActionKey(ActionKey *actionKey, vector< BasicAerodynamicForce * > *basicAerodynamicForceVector)
{
    int numberOfAerodynamicForceBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    DataContainer *dataContainer = actionKey->GetData();
    dataContainer->MoveToBegin();

    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodies; ++ iMovingBody)
    {
        (*basicAerodynamicForceVector)[iMovingBody]->WriteAerodynamicForceToDataContainer(dataContainer);
    }
}

BasicAerodynamicForce::BasicAerodynamicForce()
{
    this->ZeroBasicAerodynamicForce();
}

BasicAerodynamicForce::~BasicAerodynamicForce()
{
}

void BasicAerodynamicForce::ZeroBasicAerodynamicForce()
{
    forceX = 0.0;
    forceY = 0.0;
    forceZ = 0.0;

    forceXByPressure = 0.0;
    forceYByPressure = 0.0;
    forceZByPressure = 0.0;

    forceXByViscosity = 0.0;
    forceYByViscosity = 0.0;
    forceZByViscosity = 0.0;

    momentX = 0.0;
    momentY = 0.0;
    momentZ = 0.0;

    momentXByPressure = 0.0;
    momentYByPressure = 0.0;
    momentZByPressure = 0.0;

    momentXByViscosity = 0.0;
    momentYByViscosity = 0.0;
    momentZByViscosity = 0.0;
}

void BasicAerodynamicForce::SetMomentReferencePoint(RDouble momentReferenceX, RDouble momentReferenceY, RDouble momentReferenceZ)
{
    this->momentReferenceX = momentReferenceX;
    this->momentReferenceY = momentReferenceY;
    this->momentReferenceZ = momentReferenceZ;
}

void BasicAerodynamicForce::SetPressureAerodynamicForce(RDouble forceXByPressure, RDouble forceYByPressure, RDouble forceZByPressure)
{
    this->forceXByPressure = forceXByPressure;
    this->forceYByPressure = forceYByPressure;
    this->forceZByPressure = forceZByPressure;
}

void BasicAerodynamicForce::SetViscousAerodynamicForce(RDouble forceXByViscosity, RDouble forceYByViscosity, RDouble forceZByViscosity)
{
    this->forceXByViscosity = forceXByViscosity;
    this->forceYByViscosity = forceYByViscosity;
    this->forceZByViscosity = forceZByViscosity;
}

void BasicAerodynamicForce::ComputeResultantAerodynamicForce()
{
    forceX = forceXByPressure + forceXByViscosity;
    forceY = forceYByPressure + forceYByViscosity;
    forceZ = forceZByPressure + forceZByViscosity;
}

void BasicAerodynamicForce::ComputeResultantAerodynamicMoment()
{
    momentX = momentXByPressure + momentXByViscosity;
    momentY = momentYByPressure + momentYByViscosity;
    momentZ = momentZByPressure + momentZByViscosity;
}

void BasicAerodynamicForce::SetPressureAerodynamicMoment(RDouble faceCenterX, RDouble faceCenterY, RDouble faceCenterZ)
{
    momentXByPressure = (faceCenterY - momentReferenceY) * forceZByPressure - (faceCenterZ - momentReferenceZ) * forceYByPressure;
    momentYByPressure = (faceCenterZ - momentReferenceZ) * forceXByPressure - (faceCenterX - momentReferenceX) * forceZByPressure;
    momentZByPressure = (faceCenterX - momentReferenceX) * forceYByPressure - (faceCenterY - momentReferenceY) * forceXByPressure;
}

void BasicAerodynamicForce::SetViscousAerodynamicMoment(RDouble faceCenterX, RDouble faceCenterY, RDouble faceCenterZ)
{
    momentXByViscosity = (faceCenterY - momentReferenceY) * forceZByViscosity - (faceCenterZ - momentReferenceZ) * forceYByViscosity;
    momentYByViscosity = (faceCenterZ - momentReferenceZ) * forceXByViscosity - (faceCenterX - momentReferenceX) * forceZByViscosity;
    momentZByViscosity = (faceCenterX - momentReferenceX) * forceYByViscosity - (faceCenterY - momentReferenceY) * forceXByViscosity;
}

void BasicAerodynamicForce::ComputeAerodynamicMoment(RDouble faceCenterX, RDouble faceCenterY, RDouble faceCenterZ)
{
    momentX = (faceCenterY - momentReferenceY) * forceZ - (faceCenterZ - momentReferenceZ) * forceY;
    momentY = (faceCenterZ - momentReferenceZ) * forceX - (faceCenterX - momentReferenceX) * forceZ;
    momentZ = (faceCenterX - momentReferenceX) * forceY - (faceCenterY - momentReferenceY) * forceX;
}

void BasicAerodynamicForce::AddAerodynamicForce(BasicAerodynamicForce *basicAerodynamicForce)
{
    this->forceX += basicAerodynamicForce->forceX;
    this->forceY += basicAerodynamicForce->forceY;
    this->forceZ += basicAerodynamicForce->forceZ;

    this->forceXByPressure += basicAerodynamicForce->forceXByPressure;
    this->forceYByPressure += basicAerodynamicForce->forceYByPressure;
    this->forceZByPressure += basicAerodynamicForce->forceZByPressure;

    this->forceXByViscosity += basicAerodynamicForce->forceXByViscosity;
    this->forceYByViscosity += basicAerodynamicForce->forceYByViscosity;
    this->forceZByViscosity += basicAerodynamicForce->forceZByViscosity;

    this->momentX += basicAerodynamicForce->momentX;
    this->momentY += basicAerodynamicForce->momentY;
    this->momentZ += basicAerodynamicForce->momentZ;

    this->momentXByPressure += basicAerodynamicForce->momentXByPressure;
    this->momentYByPressure += basicAerodynamicForce->momentYByPressure;
    this->momentZByPressure += basicAerodynamicForce->momentZByPressure;

    this->momentXByViscosity += basicAerodynamicForce->momentXByViscosity;
    this->momentYByViscosity += basicAerodynamicForce->momentYByViscosity;
    this->momentZByViscosity += basicAerodynamicForce->momentZByViscosity;
}

void BasicAerodynamicForce::WriteAerodynamicForceToDataContainer(DataContainer *dataContainer)
{
    PHWrite(dataContainer, forceX);
    PHWrite(dataContainer, forceY);
    PHWrite(dataContainer, forceZ);

    PHWrite(dataContainer, forceXByPressure);
    PHWrite(dataContainer, forceYByPressure);
    PHWrite(dataContainer, forceZByPressure);

    PHWrite(dataContainer, forceXByViscosity);
    PHWrite(dataContainer, forceYByViscosity);
    PHWrite(dataContainer, forceZByViscosity);

    PHWrite(dataContainer, momentX);
    PHWrite(dataContainer, momentY);
    PHWrite(dataContainer, momentZ);

    PHWrite(dataContainer, momentXByPressure);
    PHWrite(dataContainer, momentYByPressure);
    PHWrite(dataContainer, momentZByPressure);

    PHWrite(dataContainer, momentXByViscosity);
    PHWrite(dataContainer, momentYByViscosity);
    PHWrite(dataContainer, momentZByViscosity);
}

void BasicAerodynamicForce::ReadAerodynamicForceFromDataContainer(DataContainer *dataContainer)
{
    PHRead(dataContainer, forceX);
    PHRead(dataContainer, forceY);
    PHRead(dataContainer, forceZ);

    PHRead(dataContainer, forceXByPressure);
    PHRead(dataContainer, forceYByPressure);
    PHRead(dataContainer, forceZByPressure);

    PHRead(dataContainer, forceXByViscosity);
    PHRead(dataContainer, forceYByViscosity);
    PHRead(dataContainer, forceZByViscosity);

    PHRead(dataContainer, momentX);
    PHRead(dataContainer, momentY);
    PHRead(dataContainer, momentZ);

    PHRead(dataContainer, momentXByPressure);
    PHRead(dataContainer, momentYByPressure);
    PHRead(dataContainer, momentZByPressure);

    PHRead(dataContainer, momentXByViscosity);
    PHRead(dataContainer, momentYByViscosity);
    PHRead(dataContainer, momentZByViscosity);
}

BasicAleForce::BasicAleForce()
{
    this->ZeroBasicAleForce();
}

BasicAleForce::~BasicAleForce()
{
}

void BasicAleForce::ZeroBasicAleForce()
{
    forceX = 0.0;
    forceY = 0.0;
    forceZ = 0.0;

    momentX = 0.0;
    momentY = 0.0;
    momentZ = 0.0;

    powerConsumption = 0.0;
}

void BasicAleForce::SetBasicAleForceByAerodynamicForce(AerodynamicForce *aerodynamicForce)
{
    forceX = aerodynamicForce->GetAerodynamicForceX();
    forceY = aerodynamicForce->GetAerodynamicForceY();
    forceZ = aerodynamicForce->GetAerodynamicForceZ();
    aerodynamicForce->GetMomentByOriginOfCoordinate(momentX, momentY, momentZ);
}

void BasicAleForce::ReadBasicAleForceFromBasicAleForce(BasicAleForce *basicAleForce)
{
    (*this) = (*basicAleForce);
}

void BasicAleForce::ReadBasicAleForceFromFile(fstream &file)
{
    //PHSPACE::PHRead(file, forceX);
    //PHSPACE::PHRead(file, forceY);
    //PHSPACE::PHRead(file, forceZ);
    //
    //PHSPACE::PHRead(file, momentX);
    //PHSPACE::PHRead(file, momentY);
    //PHSPACE::PHRead(file, momentZ);
    //
    //PHSPACE::PHRead(file, powerConsumption);
    file >> forceX >> forceY >> forceZ;
    file >> momentX >> momentY >> momentZ;
    file >> powerConsumption;
}

void BasicAleForce::WriteBasicAleForceToFile(fstream &file)
{
    //PHSPACE::PHWrite(file, forceX);
    //PHSPACE::PHWrite(file, forceY);
    //PHSPACE::PHWrite(file, forceZ);
    //
    //PHSPACE::PHWrite(file, momentX);
    //PHSPACE::PHWrite(file, momentY);
    //PHSPACE::PHWrite(file, momentZ);
    //
    //PHSPACE::PHWrite(file, powerConsumption);
    file << forceX << "   " << forceY << "   " << forceZ << "\n";
    
    file << momentX << "   " << momentY << "   " << momentZ << "\n";

    file << powerConsumption << "\n";
}

AleForceOfBody::AleForceOfBody()
{
    bodyIndex = 0;
    numberOfTimeLayers = 0;
}

AleForceOfBody::~AleForceOfBody()
{
    this->DeallocateMemory();
}

void AleForceOfBody::AllocateMemory()
{
    basicAleForceContainer.resize(numberOfTimeLayers);
    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++iTimeLayer)
    {
        basicAleForceContainer[iTimeLayer] = new BasicAleForce();
    }
}

void AleForceOfBody::DeallocateMemory()
{
    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++iTimeLayer)
    {
        delete basicAleForceContainer[iTimeLayer];
    }
}

BasicAleForce *AleForceOfBody::GetBasicAleForce(int timeLayerIndex)
{
    return basicAleForceContainer[timeLayerIndex];
}

void AleForceOfBody::DumpRestartAleForceOfBodyContainer(fstream &file)
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();
    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++iTimeLayer)
    {
        BasicAleForce *basicAleForce = this->GetBasicAleForce(iTimeLayer);
        basicAleForce->WriteBasicAleForceToFile(file);
    }
}

void AleForceOfBody::InitializeAleForceOfBodyContainer()
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();
    int bodyIndex = this->GetBodyIndex();
    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++iTimeLayer)
    {
        BasicAleForce *basicAleForce = this->GetBasicAleForce(iTimeLayer);
        AerodynamicForce *aerodynamicForce = PHSPACE::GetCurrentAerodynamicForce(bodyIndex);
        basicAleForce->SetBasicAleForceByAerodynamicForce(aerodynamicForce);
    }
}

void AleForceOfBody::ReadAleForceOfBodyContainer()
{
    fstream &aleFile = PHSPACE::GetCurrentAleFile();

    int numberOfTimeLayers = this->GetNumberOfTimeLayers();
    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++iTimeLayer)
    {
        BasicAleForce *basicAleForce = this->GetBasicAleForce(iTimeLayer);
        fstream &file = aleFile;
        basicAleForce->ReadBasicAleForceFromFile(file);
    }
}

void AleForceOfBody::UpdateAleForceOfBodyOfNewTimeLayer()
{
    int iTimeLayer = 0;
    int bodyIndex = this->GetBodyIndex();

    AerodynamicForce *aerodynamicForce = PHSPACE::GetCurrentAerodynamicForce(bodyIndex);

    BasicAleForce *basicAleForce = this->GetBasicAleForce(iTimeLayer);
    basicAleForce->SetBasicAleForceByAerodynamicForce(aerodynamicForce);
}

void AleForceOfBody::UpdateUnsteadyAleForceTerm()
{
    int iTimeLevel0 = 0;
    int iTimeLevel1 = 1;
    int iTimeLevel2 = 2;

    BasicAleForce *basicAleForce = this->GetBasicAleForce(iTimeLevel0);
    BasicAleForce *basicAleForceNTimeLevel = this->GetBasicAleForce(iTimeLevel1);
    BasicAleForce *basicAleForceN1TimeLevel = this->GetBasicAleForce(iTimeLevel2);

    PHSPACE::CopyBasicAleForce(basicAleForceN1TimeLevel, basicAleForceNTimeLevel);
    PHSPACE::CopyBasicAleForce(basicAleForceNTimeLevel, basicAleForce);
}

AleForceManager::AleForceManager()
{
    numberOfMovingBodies = 0;

    dimensionalForceCoefficient = 1.0;
    dimensionalMomentCoefficient = 1.0;
    dimensionalPowerCoefficient = 1.0;
    numberOfTimeLayers = 0;
}

AleForceManager::~AleForceManager()
{
    this->DeallocateAleForceOfBodyContainer();
}

void AleForceManager::Initialize()
{
    this->AllocateAleForceOfBodyContainer();
    this->InitializeAleForceOfBodyContainer();
}

void AleForceManager::Restart()
{
    this->ReadAleForceOfBodyContainer();
}

void AleForceManager::InitializeDimensionalReferenceVariables()
{
    RDouble referenceDensity = GlobalDataBase::GetDoubleParaFromDB("referenceDensity");
    RDouble referenceVelocity = GlobalDataBase::GetDoubleParaFromDB("referenceVelocity");
    RDouble referenceLength = GlobalDataBase::GetDoubleParaFromDB("referenceLength");

    int geometricDimension = GetDim();

    RDouble dynamicPressure = referenceDensity * PHSPACE::SQR(referenceVelocity);

    dimensionalForceCoefficient = half * dynamicPressure * PHSPACE::SQR(referenceLength);
    dimensionalMomentCoefficient = half * dynamicPressure * PHSPACE::POWER3(referenceLength);
    dimensionalPowerCoefficient = dimensionalForceCoefficient * referenceVelocity;

    if (geometricDimension == TWO_D)
    {
        dimensionalForceCoefficient = half * dynamicPressure * referenceLength;
        dimensionalMomentCoefficient = half * dynamicPressure * PHSPACE::SQR(referenceLength);
    }
}

void AleForceManager::AllocateAleForceOfBodyContainer()
{
    numberOfTimeLayers = 3;
    numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");
    this->InitializeDimensionalReferenceVariables();

    aleForceOfBodyContainer.resize(numberOfMovingBodies);

    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        AleForceOfBody *aleForceOfBody = new AleForceOfBody();
        aleForceOfBodyContainer[iMovingBody] = aleForceOfBody;

        aleForceOfBody->SetBodyIndex(iMovingBody);
        aleForceOfBody->SetNumberOfTimeLayers(numberOfTimeLayers);
        aleForceOfBody->AllocateMemory();
    }
}

void AleForceManager::DeallocateAleForceOfBodyContainer()
{
    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        delete aleForceOfBodyContainer[iMovingBody];
    }
}

void AleForceManager::InitializeAleForceOfBodyContainer()
{
    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        AleForceOfBody *aleForceOfBody = this->GetAleForceOfBody(iMovingBody);

        aleForceOfBody->InitializeAleForceOfBodyContainer();
    }
}

void AleForceManager::ReadAleForceOfBodyContainer()
{

    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        AleForceOfBody *aleForceOfBody = this->GetAleForceOfBody(iMovingBody);
        aleForceOfBody->ReadAleForceOfBodyContainer();

        int resetMassCenter = GetIntegerParameterFromDataBaseIfExist(iMovingBody, "resetMassCenter", 0);
        if (resetMassCenter == 1)
        {
            aleForceOfBody->InitializeAleForceOfBodyContainer();
        }
    }
}

void AleForceManager::Run()
{
    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        aleForceOfBodyContainer[iMovingBody]->UpdateAleForceOfBodyOfNewTimeLayer();
    }
}

void AleForceManager::Post()
{
    int numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        AleForceOfBody *aleForceOfBody = this->GetAleForceOfBody(iMovingBody);
        aleForceOfBody->UpdateUnsteadyAleForceTerm();
    }
}

void AleForceManager::DumpRestartAleForceOfBodyContainer(fstream &file)
{
    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        AleForceOfBody *aleForceOfBody = this->GetAleForceOfBody(iMovingBody);
        aleForceOfBody->DumpRestartAleForceOfBodyContainer(file);
    }
}

BasicAleForce *AleForceManager::GetBasicAleForce(int bodyIndex, int iTimeLayer)
{
    AleForceOfBody *aleForceOfBody = this->GetAleForceOfBody(bodyIndex);
    return aleForceOfBody->GetBasicAleForce(iTimeLayer);
}

void CopyBasicAleForce(BasicAleForce *basicAleForce, BasicAleForce *basicAleForceIn)
{
    *basicAleForce = (*basicAleForceIn);
}

AerodynamicForce::AerodynamicForce()
{
    bodyName = "AllWall";
    this->Zero();
}

AerodynamicForce::~AerodynamicForce()
{
}

void AerodynamicForce::Zero()
{
    forceX = 0.0;
    forceY = 0.0;
    forceZ = 0.0;

    forceXByPressure = 0.0;
    forceYByPressure = 0.0;
    forceZByPressure = 0.0;

    forceXByViscosity = 0.0;
    forceYByViscosity = 0.0;
    forceZByViscosity = 0.0;

    momentX = 0.0;
    momentY = 0.0;
    momentZ = 0.0;

    momentXByPressure = 0.0;
    momentYByPressure = 0.0;
    momentZByPressure = 0.0;

    momentXByViscosity = 0.0;
    momentYByViscosity = 0.0;
    momentZByViscosity = 0.0;

    coefForceX = 0.0;
    coefForceY = 0.0;
    coefForceZ = 0.0;

    coefForceXByPressure = 0.0;
    coefForceYByPressure = 0.0;
    coefForceZByPressure = 0.0;

    coefForceXByViscosity = 0.0;
    coefForceYByViscosity = 0.0;
    coefForceZByViscosity = 0.0;

    coefMomentX = 0.0;
    coefMomentY = 0.0;
    coefMomentZ = 0.0;

    coefMomentXByPressure = 0.0;
    coefMomentYByPressure = 0.0;
    coefMomentZByPressure = 0.0;

    coefMomentXByViscosity = 0.0;
    coefMomentYByViscosity = 0.0;
    coefMomentZByViscosity = 0.0;
}

void AerodynamicForce::Initialize(RDouble angleOfAttackDenotedByRadian, RDouble angleOfSideslipDenotedByRadian, RDouble referenceAreaOfAerodynamicForce, RDouble referenceLengthOfAerodynamicMoment, RDouble referenceSpainLengthOfAerodynamicMoment)
{
    this->angleOfAttackDenotedByRadian = angleOfAttackDenotedByRadian;
    this->angleOfSideslipDenotedByRadian = angleOfSideslipDenotedByRadian;
    this->referenceArea = referenceAreaOfAerodynamicForce;
    this->referenceLength = referenceLengthOfAerodynamicMoment;
    this->referenceSpainLength = referenceSpainLengthOfAerodynamicMoment;

    this->Zero();
}

void AerodynamicForce::SetTorqueRef(RDouble torqueRefX, RDouble torqueRefY, RDouble torqueRefZ)
{
    this->torqueRefX = torqueRefX;
    this->torqueRefY = torqueRefY;
    this->torqueRefZ = torqueRefZ;
};

void AerodynamicForce::SetHingeMomentRef(RDouble *localCoordAxis0, RDouble *localCoordAxis1)
{
    for (int m = 0; m < 3; ++m)
    {
        this->localCoordAxis0[m] = localCoordAxis0[m];
        this->localCoordAxis1[m] = localCoordAxis1[m];
    }
};

void AerodynamicForce::Add(DataContainer *dataContainer)
{
    const int numberOfForceVariables = 18;
    RDouble forceMomentContainer[numberOfForceVariables];

    dataContainer->MoveToBegin();

    PHSPACE::PHRead(dataContainer, forceMomentContainer, numberOfForceVariables);

    forceX += forceMomentContainer[0];
    forceY += forceMomentContainer[1];
    forceZ += forceMomentContainer[2];

    forceXByPressure += forceMomentContainer[3];
    forceYByPressure += forceMomentContainer[4];
    forceZByPressure += forceMomentContainer[5];

    forceXByViscosity += forceMomentContainer[6];
    forceYByViscosity += forceMomentContainer[7];
    forceZByViscosity += forceMomentContainer[8];

    momentX += forceMomentContainer[9];
    momentY += forceMomentContainer[10];
    momentZ += forceMomentContainer[11];

    momentXByPressure += forceMomentContainer[12];
    momentYByPressure += forceMomentContainer[13];
    momentZByPressure += forceMomentContainer[14];

    momentXByViscosity += forceMomentContainer[15];
    momentYByViscosity += forceMomentContainer[16];
    momentZByViscosity += forceMomentContainer[17];
}

void AerodynamicForce::GetMomentByOriginOfCoordinate(RDouble &aerodynamicMomentX0, RDouble &aerodynamicMomentY0, RDouble &aerodynamicMomentZ0)
{
    RDouble momentReferencePoint[3];

    momentReferencePoint[0] = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    momentReferencePoint[1] = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    momentReferencePoint[2] = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

    RDouble aerodynamicForce[3];

    aerodynamicForce[0] = forceX;
    aerodynamicForce[1] = forceY;
    aerodynamicForce[2] = forceZ;

    RDouble changeOfAerodynamicMoment[3];
    PHSPACE::CrossProduct(momentReferencePoint, aerodynamicForce, changeOfAerodynamicMoment);

    aerodynamicMomentX0 = momentX + changeOfAerodynamicMoment[0];
    aerodynamicMomentY0 = momentY + changeOfAerodynamicMoment[1];
    aerodynamicMomentZ0 = momentZ + changeOfAerodynamicMoment[2];
}

void AerodynamicForce::DumpCoefficient(ActionKey *actionKey)
{
    using namespace PHMPI;
    int myID = GetCurrentProcessorID();
    int serverTmp = GetServerProcessorID();
    if (myID != serverTmp) return;

    PHSPACE::ParallelOpenFile(actionKey);

    fstream &file = *actionKey->GetFile();

    this->DumpCoefficient(file);

    PHSPACE::ParallelCloseFile(file);
}

void AerodynamicForce::DumpComponentCoefficient(ActionKey *actionKey)
{
    using namespace PHMPI;
    int myID = GetCurrentProcessorID();
    int serverTmp = GetServerProcessorID();
    if (myID != serverTmp) return;

    string totalForceFileName = actionKey->GetFileName();
    string partAircoefFileName = PHSPACE::AddSymbolToFileName(totalForceFileName, "_", bodyName);

    fstream file;
    file.open(partAircoefFileName.c_str(), ios_base::out | ios_base::app);

    this->DumpCoefficient(file);

    file.close();
    file.clear();
}

void AerodynamicForce::DumpCoefficient(fstream &file)
{
    this->ComputeCoefficient();

    ostringstream oss;

    int outerIterationSteps = GlobalDataBase::GetIntParaFromDB("outnstep");

    oss << setiosflags(ios::right);
    oss << setiosflags(ios::scientific);

    if (PHSPACE::IfFileEmpty(file))
    {
        vector<string> titleOfTecplot;
        titleOfTecplot.push_back("Title=\"THE AIRCOEF\"");
        titleOfTecplot.push_back("Variables=");
        titleOfTecplot.push_back("iter");
        titleOfTecplot.push_back("<i>C<sub>L</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>D</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>Z</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>D_p</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>D_f</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>D</sub></i>-<i>C<sub>L</sub></i><sup>2</sup>/(<greek>p</greek><math>4</math>A<sub>r</sub>)</i>");
        titleOfTecplot.push_back("<i>X<sub>cp</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>A</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>N</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>Z1</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>ml</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>mn</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>m</sub></i>");
        if (bodyName == "AllWall")
        {
            titleOfTecplot.push_back("\"WallTime\"");
        }

        for (std::size_t iTitle = 0; iTitle < titleOfTecplot.size(); ++iTitle)
        {
            oss << titleOfTecplot[iTitle] << "\n";
        }
    }

    oss << outerIterationSteps << "\t";

    oss << setprecision(10);

    oss << totalLiftCoefficient << "\t";
    oss << totalDragCoefficient << "\t";
    oss << totalCrossCoefficient << "\t";
    oss << totalPressureDragCoefficient << "\t";
    oss << totalViscosityDragCoefficient << "\t";
    oss << DragCoefLiftCoef2piar << "\t";
    oss << pressureCenter << "\t";

    oss << coefForceX << "\t";
    oss << coefForceY << "\t";
    oss << coefForceZ << "\t";
    oss << coefMomentX << "\t";
    oss << coefMomentY << "\t";
    oss << coefMomentZ << "\t";

    if (bodyName == "AllWall")
    {
        RDouble wallTime = PHSPACE::TIME_SPACE::GetWallTime();
        oss << wallTime;
    }

    oss << "\n";
    PHSPACE::WriteASCIIFile(file, oss.str());
}

void AerodynamicForce::ComputeCoefficient()
{
    RDouble sina = sin(angleOfAttackDenotedByRadian);
    RDouble cosa = cos(angleOfAttackDenotedByRadian);
    RDouble sinb = sin(angleOfSideslipDenotedByRadian);
    RDouble cosb = cos(angleOfSideslipDenotedByRadian);

    const RDouble AR = 9.5;

    RDouble referenceCoefficientOfForce = 1.0 / referenceArea;
    RDouble referenceCoefficientOfMoment = 1.0 / (referenceArea * referenceLength);
    RDouble referenceCoefficientOfMomentSpain = 1.0 / (referenceArea * referenceSpainLength);

    coefForceX = forceX * referenceCoefficientOfForce;
    coefForceY = forceY * referenceCoefficientOfForce;
    coefForceZ = forceZ * referenceCoefficientOfForce;

    coefForceXByPressure = forceXByPressure * referenceCoefficientOfForce;
    coefForceYByPressure = forceYByPressure * referenceCoefficientOfForce;
    coefForceZByPressure = forceZByPressure * referenceCoefficientOfForce;

    coefForceXByViscosity = forceXByViscosity * referenceCoefficientOfForce;
    coefForceYByViscosity = forceYByViscosity * referenceCoefficientOfForce;
    coefForceZByViscosity = forceZByViscosity * referenceCoefficientOfForce;

    coefMomentX = momentX * referenceCoefficientOfMomentSpain;
    coefMomentY = momentY * referenceCoefficientOfMomentSpain;
    coefMomentZ = momentZ * referenceCoefficientOfMoment;

    coefMomentXByPressure = momentXByPressure * referenceCoefficientOfMomentSpain;
    coefMomentYByPressure = momentYByPressure * referenceCoefficientOfMomentSpain;
    coefMomentZByPressure = momentZByPressure * referenceCoefficientOfMoment;

    coefMomentXByViscosity = momentXByViscosity * referenceCoefficientOfMomentSpain;
    coefMomentYByViscosity = momentYByViscosity * referenceCoefficientOfMomentSpain;
    coefMomentZByViscosity = momentZByViscosity * referenceCoefficientOfMoment;

    totalLiftCoefficient = -sina * coefForceX + cosa * coefForceY;
    totalDragCoefficient = cosa * cosb * coefForceX + sina * cosb * coefForceY + sinb * coefForceZ;
    totalCrossCoefficient = -cosa * sinb * coefForceX - sina * sinb * coefForceY + cosb * coefForceZ;
    totalPressureDragCoefficient = cosa * cosb * coefForceXByPressure + sina * cosb * coefForceYByPressure + sinb * coefForceZByPressure;
    totalViscosityDragCoefficient = cosa * cosb * coefForceXByViscosity + sina * cosb * coefForceYByViscosity + sinb * coefForceZByViscosity;
    DragCoefLiftCoef2piar = totalDragCoefficient - totalLiftCoefficient * totalLiftCoefficient / (PI * AR);
    pressureCenter = PHSPACE::SIGN(1.0, coefForceY) * coefMomentZ / (PHSPACE::ABS(coefForceY) + SMALL);
}

void AerodynamicForce::CollectOfAllProcessors()
{
    using namespace PHMPI;
    RDouble forceXSum = forceX;
    RDouble forceYSum = forceY;
    RDouble forceZSum = forceZ;

    RDouble forceXByPressureSum = forceXByPressure;
    RDouble forceYByPressureSum = forceYByPressure;
    RDouble forceZByPressureSum = forceZByPressure;

    RDouble forceXByViscositySum = forceXByViscosity;
    RDouble forceYByViscositySum = forceYByViscosity;
    RDouble forceZByViscositySum = forceZByViscosity;

    RDouble momentXSum = momentX;
    RDouble momentYSum = momentY;
    RDouble momentZSum = momentZ;

    RDouble momentXByPressureSum = momentXByPressure;
    RDouble momentYByPressureSum = momentYByPressure;
    RDouble momentZByPressureSum = momentZByPressure;

    RDouble momentXByViscositySum = momentXByViscosity;
    RDouble momentYByViscositySum = momentYByViscosity;
    RDouble momentZByViscositySum = momentZByViscosity;

    PH_Reduce(&forceX, &forceXSum, 1, PH_SUM);
    PH_Reduce(&forceY, &forceYSum, 1, PH_SUM);
    PH_Reduce(&forceZ, &forceZSum, 1, PH_SUM);

    PH_Reduce(&forceXByPressure, &forceXByPressureSum, 1, PH_SUM);
    PH_Reduce(&forceYByPressure, &forceYByPressureSum, 1, PH_SUM);
    PH_Reduce(&forceZByPressure, &forceZByPressureSum, 1, PH_SUM);

    PH_Reduce(&forceXByViscosity, &forceXByViscositySum, 1, PH_SUM);
    PH_Reduce(&forceYByViscosity, &forceYByViscositySum, 1, PH_SUM);
    PH_Reduce(&forceZByViscosity, &forceZByViscositySum, 1, PH_SUM);

    PH_Reduce(&momentX, &momentXSum, 1, PH_SUM);
    PH_Reduce(&momentY, &momentYSum, 1, PH_SUM);
    PH_Reduce(&momentZ, &momentZSum, 1, PH_SUM);

    PH_Reduce(&momentXByPressure, &momentXByPressureSum, 1, PH_SUM);
    PH_Reduce(&momentYByPressure, &momentYByPressureSum, 1, PH_SUM);
    PH_Reduce(&momentZByPressure, &momentZByPressureSum, 1, PH_SUM);

    PH_Reduce(&momentXByViscosity, &momentXByViscositySum, 1, PH_SUM);
    PH_Reduce(&momentYByViscosity, &momentYByViscositySum, 1, PH_SUM);
    PH_Reduce(&momentZByViscosity, &momentZByViscositySum, 1, PH_SUM);

    forceX = forceXSum;
    forceY = forceYSum;
    forceZ = forceZSum;

    forceXByPressure = forceXByPressureSum;
    forceYByPressure = forceYByPressureSum;
    forceZByPressure = forceZByPressureSum;

    forceXByViscosity = forceXByViscositySum;
    forceYByViscosity = forceYByViscositySum;
    forceZByViscosity = forceZByViscositySum;

    momentX = momentXSum;
    momentY = momentYSum;
    momentZ = momentZSum;

    momentXByPressure = momentXByPressureSum;
    momentYByPressure = momentYByPressureSum;
    momentZByPressure = momentZByPressureSum;

    momentXByViscosity = momentXByViscositySum;
    momentYByViscosity = momentYByViscositySum;
    momentZByViscosity = momentZByViscositySum;
}

AerodynamicForceManager::AerodynamicForceManager()
{
    isAllocateMemory = false;
    numberOfAerodynamicForceBodies = 0;
    bodyID = 0;
}

AerodynamicForceManager::~AerodynamicForceManager()
{
    this->DeallocateMemory();
}

void AerodynamicForceManager::AllocateMemory()
{
    if (isAllocateMemory) return;

    isAllocateMemory = true;

    aerodynamicForceContainer.resize(numberOfAerodynamicForceBodies);

    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodies; ++ iMovingBody)
    {
        aerodynamicForceContainer[iMovingBody] = new AerodynamicForce();
    }
}

void AerodynamicForceManager::DeallocateMemory()
{
    if (!isAllocateMemory) return;

    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodies; ++ iMovingBody)
    {
        delete aerodynamicForceContainer[iMovingBody];
    }
}

void AerodynamicForceManager::InitializeAerodynamicForceManager()
{
    numberOfAerodynamicForceBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");
    this->AllocateMemory();

    RDouble angleOfAttackDenotedByRadian = GlobalDataBase::GetDoubleParaFromDB("attack");
    RDouble angleOfSideslipDenotedByRadian = GlobalDataBase::GetDoubleParaFromDB("sideslip");

    RDouble referenceAreaOfAerodynamicForce = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
    RDouble referenceLengthOfAerodynamicMoment = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
    RDouble referenceSpainLengthOfAerodynamicMoment = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");

    set < string > *bodyList = GlobalBoundaryCondition::GetBodyNameList();
    set < string >::iterator bodyListIter = bodyList->begin();
    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodies; ++ iMovingBody)
    {
        AerodynamicForce *aerodynamicForce = this->GetAerodynamicForce(iMovingBody);
        aerodynamicForce->Initialize(angleOfAttackDenotedByRadian, angleOfSideslipDenotedByRadian, referenceAreaOfAerodynamicForce,
            referenceLengthOfAerodynamicMoment, referenceSpainLengthOfAerodynamicMoment);

        string bodyName = *bodyListIter;
        aerodynamicForce->SetBodyName(bodyName);
        bodyListIter ++;
    }
}

void AerodynamicForceManager::ExtractForceDataFromActionKey(ActionKey *actionKey, vector< DataContainer * > *dataContainerArray)
{
    DataContainer *dataContainer = actionKey->GetData();
    dataContainer->MoveToBegin();

    const int numberOfData = 18;
    RDouble bufferForSwap[numberOfData];

    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodies; ++ iMovingBody)
    {
        PHSPACE::PHRead(dataContainer, bufferForSwap, numberOfData);
        PHSPACE::PHWrite((*dataContainerArray)[iMovingBody], bufferForSwap, numberOfData);
    }
}

vector< DataContainer * > *AerodynamicForceManager::CreateDataContainerArray()
{
    vector< DataContainer * > *dataContainerArray = new vector< DataContainer * >(numberOfAerodynamicForceBodies + 1);
    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodies; ++ iMovingBody)
    {
        (*dataContainerArray)[iMovingBody] = new DataContainer();
    }
    return dataContainerArray;
}

void AerodynamicForceManager::FreeDataContainerArray(vector< DataContainer * > *dataContainerArray)
{
    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodies; ++ iMovingBody)
    {
        delete (*dataContainerArray)[iMovingBody];
    }
    delete dataContainerArray;
}

void AerodynamicForceManager::AddAerodynamicForce(ActionKey *actionKey)
{
    vector< DataContainer * > *dataContainerArray = this->CreateDataContainerArray();

    this->ExtractForceDataFromActionKey(actionKey, dataContainerArray);

    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodies; ++ iMovingBody)
    {
        DataContainer *dataContainer = (*dataContainerArray)[iMovingBody];
        aerodynamicForceContainer[iMovingBody]->Add(dataContainer);
    }
    this->FreeDataContainerArray(dataContainerArray);
}

void AerodynamicForceManager::CollectAerodynamicForceOfAllProcessors()
{
    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodies; ++ iMovingBody)
    {
        aerodynamicForceContainer[iMovingBody]->CollectOfAllProcessors();
    }
}

AerodynamicForceManager *GetAerodynamicForceManager()
{
    return &aerodynamicForceManager;
}

AerodynamicForce *GetCurrentAerodynamicForce(int bodyIndex)
{
    return GetAerodynamicForceManager()->GetAerodynamicForce(bodyIndex);
}

BasicAerodynamicForce *GetBasicAerodynamicForce(int iBoundaryConditionRegion, BasicAerodynamicForce **basicAerodynamicForceContainer)
{
    return basicAerodynamicForceContainer[iBoundaryConditionRegion];
}

BasicAerodynamicForce *GetBasicAerodynamicForce(int iBoundaryConditionRegion, vector<BasicAerodynamicForce *> *basicAerodynamicForceVector)
{
    return (*basicAerodynamicForceVector)[iBoundaryConditionRegion];
}

bool OutPutSubIterationForceHeader()
{
    int outerIterationSteps = GlobalDataBase::GetIntParaFromDB("outnstep");
    int innerIterationSteps = GlobalDataBase::GetIntParaFromDB("innstep");
    int numberOfForceSaveSteps = GlobalDataBase::GetIntParaFromDB("intervalStepForce");

    if (outerIterationSteps / numberOfForceSaveSteps == 1 && innerIterationSteps / numberOfForceSaveSteps == 1)
    {
        return true;
    }

    return false;
}

}
