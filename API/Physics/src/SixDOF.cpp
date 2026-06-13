#include "PHMpi.h"
#include "SixDOF.h"
#include "Force.h"
#include "AleForceManager.h"
#include "Zone.h"
#include "GlobalDataBase.h"
#include "GridType.h"
#include "PHHeader.h"
#include "Math_BasisFunction.h"
#include "PHSuperMatrix.h"
#include "TK_Log.h"


namespace PHSPACE
{
PHRigidBodyMovement RigidBodyMovement;

PHRigidBodyMovement::PHRigidBodyMovement()
{
    massCenterOfAircraft    = new PHSuperMatrix(3);
    inertiaTensorOfAircraft = new PHSuperMatrix(3, 3);   //  moment of inertia of aircraft
    reverseOfInertiaTensor  = new PHSuperMatrix(3, 3);   //  the inverse matrix of moment of inertia

    // the data changed with the flow solvers.
    // input£ºforce and moment coefficient.
    aerodynamicForce  = new PHSuperMatrix(3);  //! aerodynamic force.
    aerodynamicMoment = new PHSuperMatrix(3);  //! aerodynamic moment.
    otherForceSum     = new PHSuperMatrix(3);
    otherMomentSum    = new PHSuperMatrix(3);

    eulerAngle             = new PHSuperMatrix(3);
    eulerAngleRate         = new PHSuperMatrix(3, 3);
    eulerAngleAcceleration = new PHSuperMatrix(3);

    transformMatrixCurrent = new PHSuperMatrix(3, 3);
    transformMatrixLastTime    = new PHSuperMatrix(3, 3);

    momentOfInertia        = new PHSuperMatrix(3, 3);
    bodyAngleRate          = new PHSuperMatrix(3);
    bodyAngleAcceleration  = new PHSuperMatrix(3);
}

PHRigidBodyMovement::~PHRigidBodyMovement()
{
    delete massCenterOfAircraft;
    delete inertiaTensorOfAircraft;
    delete reverseOfInertiaTensor;

    delete aerodynamicForce;
    delete aerodynamicMoment;
    delete otherForceSum;
    delete otherMomentSum;
    delete eulerAngle;
    delete eulerAngleRate;
    delete eulerAngleAcceleration;

    delete transformMatrixCurrent;
    delete transformMatrixLastTime;
    delete momentOfInertia;
    delete bodyAngleRate;
    delete bodyAngleAcceleration;
}

void PHRigidBodyMovement::UpdateCoordinate(Grid * gridIn)
{
    StructGrid * grid = static_cast<StructGrid* >(gridIn);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble3D& x = *grid->GetStructX();
    RDouble3D& y = *grid->GetStructY();
    RDouble3D& z = *grid->GetStructZ();

    PHSuperMatrix xyzc = * massCenterOfAircraft;
    PHSuperMatrix dxyz(3);
    PHSuperMatrix xyz(3);

    PHSuperMatrix & transformMatrixCurrent  = * this->transformMatrixCurrent;
    PHSuperMatrix & transformMatrixLastTime = * this->transformMatrixLastTime;

    PHSuperMatrix matrix = PHSPACE::TransposeMatrix(transformMatrixCurrent) * transformMatrixLastTime;

    for (int k = 1; k <= nk; ++ k)
    {
        for (int j = 1; j <= nj; ++ j)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                dxyz(1) = x(i, j, k) - xyzc(1);
                dxyz(2) = y(i, j, k) - xyzc(2);
                dxyz(3) = z(i, j, k) - xyzc(3);

                xyz = matrix * dxyz;

                x(i, j, k) = xyz(1) + xyzc(1);
                y(i, j, k) = xyz(2) + xyzc(2);
                z(i, j, k) = xyz(3) + xyzc(3);
            }
        }
    }

    return;
}

void PHRigidBodyMovement::UpdateGridVelocity(Grid *gridIn)
{
    StructGrid *grid = static_cast<StructGrid* >(gridIn);

    PHSuperMatrix velocityCellface(3);

    RDouble4D& faceNormalX = *(grid->GetFaceNormalX());
    RDouble4D& faceNormalY = *(grid->GetFaceNormalY());
    RDouble4D& faceNormalZ = *(grid->GetFaceNormalZ());
    RDouble4D& faceNormalVelocity = *(grid->GetFaceNormalVelocity());

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble xt, yt, zt, kt;
    int nd;

    nd = 1;
    for (int k = 1; k < nk; ++ k)
    {
        for (int j = 1; j < nj; ++ j)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                velocityCellface = CalculateGridCellFaceVelocity(grid, i, j, k, nd);

                xt = velocityCellface(1);
                yt = velocityCellface(2);
                zt = velocityCellface(3);

                kt = - (xt * faceNormalX(i, j, k, nd) + 
                         yt * faceNormalY(i, j, k, nd) + 
                         zt * faceNormalZ(i, j, k, nd));
                    
                faceNormalVelocity(i, j, k, nd) = kt;
            }
        }
    }

    nd = 2;
    for (int k = 1; k < nk; ++ k)
    {
        for (int i = 1; i < ni; ++ i)
        {
            for (int j = 1; j <= nj; ++ j)
            {
                velocityCellface = CalculateGridCellFaceVelocity(grid, i, j, k, nd);

                xt = velocityCellface(1);
                yt = velocityCellface(2);
                zt = velocityCellface(3);

                kt = - (xt * faceNormalX(i, j, k, nd) + 
                         yt * faceNormalY(i, j, k, nd) + 
                         zt * faceNormalZ(i, j, k, nd));
                    
                faceNormalVelocity(i, j, k, nd) = kt;
            }
        }
    }

    nd = 3;
    for (int j = 1; j < nj; ++ j)
    {
        for (int i = 1; i < ni; ++ i)
        {
            for (int k = 1; k <= nk; ++ k)
            {
                velocityCellface = CalculateGridCellFaceVelocity(grid, i, j, k, nd);

                xt = velocityCellface(1);
                yt = velocityCellface(2);
                zt = velocityCellface(3);

                kt = - (xt * faceNormalX(i, j, k, nd) +
                         yt * faceNormalY(i, j, k, nd) + 
                         zt * faceNormalZ(i, j, k, nd));
                    
                faceNormalVelocity(i, j, k, nd) = kt;
            }
        }
    }
}

void PHRigidBodyMovement::InitializeAircraftProperty()
{
    //! init.
    InitializeAerodynamicForce();

    RDouble massCenter[3];
    GlobalDataBase::GetData("massCenter", &massCenter, PHDOUBLE, 3);
    PHSuperMatrix& massCenterOfAircraft = *this->massCenterOfAircraft;
    massCenterOfAircraft(1) = massCenter[0];
    massCenterOfAircraft(2) = massCenter[1];
    massCenterOfAircraft(3) = massCenter[2];

    RDouble massMatrix[6];
    GlobalDataBase::GetData("massMatrix", &massMatrix, PHDOUBLE, 6);

    PHSuperMatrix& inertiaTensorOfAircraft = *this->inertiaTensorOfAircraft;
    inertiaTensorOfAircraft(1, 1) = massMatrix[0];
    inertiaTensorOfAircraft(2, 2) = massMatrix[1];
    inertiaTensorOfAircraft(3, 3) = massMatrix[2];

    //Ixy = Iyx
    inertiaTensorOfAircraft(1, 2) = massMatrix[3];
    inertiaTensorOfAircraft(2, 1) = massMatrix[3];

    //Ixz = Izx
    inertiaTensorOfAircraft(1, 3) = massMatrix[4];
    inertiaTensorOfAircraft(3, 1) = massMatrix[4];

    //Iyz = Izy
    inertiaTensorOfAircraft(2, 3) = massMatrix[5];
    inertiaTensorOfAircraft(3, 2) = massMatrix[5];

    massOfAircraft = GlobalDataBase::GetDoubleParaFromDB("mass");

    RDouble attitudeAngle[3];
    GlobalDataBase::GetData("attitudeAngle", &attitudeAngle, PHDOUBLE, 3);
    massCenterOfAircraft(1) = attitudeAngle[0];
    massCenterOfAircraft(2) = attitudeAngle[1];
    massCenterOfAircraft(3) = attitudeAngle[2];

    RDouble referenceAreaOfAerodynamicForce    = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
    RDouble referenceLengthOfAerodynamicMoment = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");

    RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");

    densityRef = refDimensionalDensity;
    areaRef    = referenceAreaOfAerodynamicForce;
    lengthRef  = referenceLengthOfAerodynamicMoment;
}

void PHRigidBodyMovement::InitializeEulerAngleInformation()
{
    //! init the Euler angle information.
    PHSuperMatrix & eulerAngle             = * this->eulerAngle;
    PHSuperMatrix & eulerAngleRate         = * this->eulerAngleRate;
    PHSuperMatrix & eulerAngleAcceleration = * this->eulerAngleAcceleration;

    for (int i = 1; i <= 3 ; ++ i)
    {
        eulerAngle(i) = 0.0; 
        for (int j = 1; j <= 3; ++ j)
        {
            eulerAngleRate(i, j) = 0.0;
        }
        eulerAngleAcceleration(i) = 0.0;      
    }
}

void PHRigidBodyMovement::CalculationOfMatrix()
{
    PHSuperMatrix & eulerAngle = * this->eulerAngle;

    RDouble pitchAngle = eulerAngle(1);
    RDouble yawAngle   = eulerAngle(2);
    RDouble rollAngle  = eulerAngle(3);

    PHSuperMatrix & transformMatrixCurrent = * this->transformMatrixCurrent;
    PHSuperMatrix & transformMatrixLastTime    = * this->transformMatrixLastTime;

    transformMatrixCurrent = 
        BaseMatrixOfRotate(rollAngle ,"x") * 
        BaseMatrixOfRotate(pitchAngle,"z") * 
        BaseMatrixOfRotate(yawAngle  ,"y");

    transformMatrixLastTime = transformMatrixCurrent;

    AngleRateTransformMatrixBodyToGlobal(eulerAngle,"y-z-x");

}

void PHRigidBodyMovement::NondimensionalOf6DOF()
{
    //! non-dimensional six DOF equations.
    RDouble MassRef, InertiaRef;

    MassRef    = 0.5 * densityRef * areaRef * lengthRef;
    InertiaRef = MassRef;

    massOfAircraft /= MassRef;

    PHSuperMatrix & inertiaTensorOfAircraft = * this->inertiaTensorOfAircraft;

    for (int i = 1; i <= 3; ++ i)
    {
        for (int j = 1; j <= 3; ++ j)
        {
            inertiaTensorOfAircraft(i, j) /= InertiaRef;
        }
    }

    PHSuperMatrix & reverseOfInertiaTensor = * this->reverseOfInertiaTensor;
    reverseOfInertiaTensor = Matrix3X3Inverse(inertiaTensorOfAircraft);
}

void PHRigidBodyMovement::GetReferenceValue()
{
    //! Get non-dimensional reference values.
    RDouble referenceAreaOfAerodynamicForce = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
    RDouble referenceLengthOfAerodynamicMoment = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");

    RDouble referenceDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");

    densityRef = referenceDimensionalDensity;
    areaRef    = referenceAreaOfAerodynamicForce;
    lengthRef  = referenceLengthOfAerodynamicMoment;
}

void PHRigidBodyMovement::AerodynamicsInCFDSolver()
{
    int bodyIndex = 0;
    AerodynamicForce * aerodynamicForce = PHSPACE::GetCurrentAerodynamicForce(bodyIndex);

    //! Get aerodynamic coefficients.
    PHSuperMatrix & aerodynamicForceCoefficient = * this->aerodynamicForce;
    aerodynamicForceCoefficient(1) = aerodynamicForce->GetForceCoefficientX();
    aerodynamicForceCoefficient(2) = aerodynamicForce->GetForceCoefficientY();
    aerodynamicForceCoefficient(3) = aerodynamicForce->GetForceCoefficientZ();

    PHSuperMatrix & aerodynamicMomentCoefficient = * this->aerodynamicMoment;

    aerodynamicMomentCoefficient(1) = aerodynamicForce->GetMomentCoefficientX();
    aerodynamicMomentCoefficient(2) = aerodynamicForce->GetMomentCoefficientY();
    aerodynamicMomentCoefficient(3) = aerodynamicForce->GetMomentCoefficientZ();
}

void PHRigidBodyMovement::UpdateDynamicInfo6DOF()
{
    //! update the six DOF dynamic information, include Euler angle solution and coordinates.
    UpdateOfInertia();
    CalculateAngleRateOfBody();
    CalculateEulerAngleRate();
    CalculateEulerAngle();
    CalculateEulerAngleAcceleration();
    CalculateAngleOfAttackAndSlide();
    CalculationOfMatrix();
}

void PHRigidBodyMovement::InitializeInertia()          
{
    PHSuperMatrix newMomentOfInertia(3);

    PHSuperMatrix & inertiaTensorOfAircraft = * this->inertiaTensorOfAircraft;
    PHSuperMatrix & momentOfInertia = * this->momentOfInertia;
    PHSuperMatrix & bodyAngleRate   = * this->bodyAngleRate;

    newMomentOfInertia = inertiaTensorOfAircraft * bodyAngleRate;

    for (int i = 1; i <= 3; ++ i)
    {
        momentOfInertia(i, 1) = newMomentOfInertia(i);
        momentOfInertia(i, 2) = momentOfInertia(i, 1);
        momentOfInertia(i, 3) = momentOfInertia(i, 2);
    }
}

void PHRigidBodyMovement::UpdateOfInertia()          
{
    //! update the angular momentum information.
    RDouble c1, c2, c3;

    c1 =  1.5 / currentMarchTimeStep;
    c2 = -2.0 / currentMarchTimeStep;
    c3 =  0.5 / currentMarchTimeStep;

    PHSuperMatrix currentInertia(3);
    PHSuperMatrix newMomentOfInertia(3);

    PHSuperMatrix & momentOfInertia = * this->momentOfInertia;
    PHSuperMatrix & bodyAngleRate   = * this->bodyAngleRate;

    for (int i = 1; i <= 3; ++ i)
    {
        currentInertia(i) = momentOfInertia(i, 1);
    }

    PHSuperMatrix crossResult = CrossProductVector(bodyAngleRate, currentInertia);

    PHSuperMatrix & aerodynamicMoment = * this->aerodynamicMoment;
    PHSuperMatrix & otherMomentSum    = * this->otherMomentSum;

    for (int i = 1; i <= 3; ++ i)
    {
        newMomentOfInertia(i) = ( aerodynamicMoment(i) + otherMomentSum(i) 
                                    - crossResult(i)
                                    - c2 * currentInertia(i) - c3 * momentOfInertia(2, 1)
                                   ) / c1;
    }

    for (int i = 1; i <= 3; ++ i)
    {
        momentOfInertia(i, 3) = momentOfInertia(i, 2);
        momentOfInertia(i, 2) = momentOfInertia(i, 1);
        momentOfInertia(i, 1) = newMomentOfInertia(i);
    }
}

void PHRigidBodyMovement::CalculateAngleRateOfBody()          
{
    //! calculate the angular velocity of the projectile coordinate system.

    PHMatrix< RDouble > currentInertia(3);

    PHSuperMatrix & momentOfInertia = * this->momentOfInertia;

    for (int i = 1; i <= 3; ++ i)
    {
        currentInertia(i) = momentOfInertia(i, 1);
    }

    PHSuperMatrix & inertiaTensorOfAircraft = * this->inertiaTensorOfAircraft;
    PHSuperMatrix & bodyAngleRate   = * this->bodyAngleRate;

    bodyAngleRate = Matrix3X3Inverse(inertiaTensorOfAircraft) * currentInertia;

}

void PHRigidBodyMovement::CalculateEulerAngleRate()          
{
    //! calculate the angular velocity of the coordinate system.

    PHSuperMatrix newAnglerate(3);

    PHSuperMatrix & transformMatrixCurrent = * this->transformMatrixCurrent;
    PHSuperMatrix & bodyAngleRate   = * this->bodyAngleRate;

    newAnglerate = PHSPACE::TransposeMatrix(transformMatrixCurrent) * bodyAngleRate;

    PHSuperMatrix & eulerAngleRate = * this->eulerAngleRate;

    for (int i = 1; i <= 3; ++ i)
    {
        eulerAngleRate(i, 3) = eulerAngleRate(i, 2);
        eulerAngleRate(i, 2) = eulerAngleRate(i, 1);
        eulerAngleRate(i, 1) = newAnglerate(i);
    }
}

void PHRigidBodyMovement::CalculateEulerAngle()
{
    //! Calculate the euler angle.

    PHSuperMatrix & eulerAngle = * this->eulerAngle;
    PHSuperMatrix & eulerAngleRate = * this->eulerAngleRate;

    for (int i = 1; i <= 3; ++ i)
    {
        eulerAngle(i) += 0.5 * (eulerAngleRate(i, 1) + eulerAngleRate(i, 2)) * currentMarchTimeStep;
    }
}

void PHRigidBodyMovement::CalculateEulerAngleAcceleration()          
{
    //! Calculate the euler angle acceleration of the computing axis systems.
    RDouble c1, c2, c3;

    c1 =  1.5 / currentMarchTimeStep;
    c2 = -2.0 / currentMarchTimeStep;
    c3 =  0.5 / currentMarchTimeStep;

    PHSuperMatrix & eulerAngleRate         = * this->eulerAngleRate;
    PHSuperMatrix & eulerAngleAcceleration = * this->eulerAngleAcceleration;

    for (int i = 1; i <= 3; ++ i)
    {
        eulerAngleAcceleration(i) = c1 * eulerAngleRate(i, 1) + 
                                      c2 * eulerAngleRate(i, 2) + 
                                      c3 * eulerAngleRate(i, 3) ;
    }
}

void PHRigidBodyMovement::CalculateAngleOfAttackAndSlide()          
{
    //! Calculate the initiative angle of attack and angle of slid in projectile axis system.
    PHSuperMatrix velocityOfInflow(3);
    PHSuperMatrix velocityInBody(3);

    velocityOfInflow(1) = cos(initiativeAngleOfAttack) * cos(initiativeAngleOfSlide);
    velocityOfInflow(2) = cos(initiativeAngleOfAttack) * sin(initiativeAngleOfSlide);
    velocityOfInflow(3) = sin(initiativeAngleOfAttack) ;
    PHSuperMatrix & velocityOfInflow1 = velocityOfInflow;

    RDouble ARC2DEG = 180.0 / acos(- 1.0);

    temporaryAngleOfAttack = ARC2DEG * atan(velocityOfInflow(2) / PHSPACE::MAX(velocityOfInflow(1), 1.0e-20));
    temporaryAngleOfSlide  = ARC2DEG * asin(velocityOfInflow(3));

    PHSuperMatrix & transformMatrixCurrent = * this->transformMatrixCurrent;

    velocityInBody = transformMatrixCurrent * velocityOfInflow1;
}

void PHRigidBodyMovement::UpdateSingleZoneGridInformation(Grid * gridIn)
{
    if (PHSPACE::UNSTRUCTGRID == gridIn->Type())
    {
        return;
    }

    //! update the mesh information, include the mesh coordinate and derivatives.
    this->UpdateCoordinate(gridIn);
    this->UpdateCellArea(gridIn);
    this->UpdateGridVelocity(gridIn);
}

void PHRigidBodyMovement::UpdateCellArea(Grid *gridIn)
{
    //! update the mesh derivatives.
    StructGrid *grid = static_cast<StructGrid*>(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D& faceNormalX = *(grid->GetFaceNormalX());
    RDouble4D& faceNormalY = *(grid->GetFaceNormalY());
    RDouble4D& faceNormalZ = *(grid->GetFaceNormalZ());

    PHSuperMatrix kxyz(3);

    PHSuperMatrix & transformMatrixCurrent = *(this->transformMatrixCurrent);
    PHSuperMatrix & transformMatrixLastTime    = *(this->transformMatrixLastTime   );

    PHMatrix< RDouble > matrix = PHSPACE::TransposeMatrix(transformMatrixCurrent) * transformMatrixLastTime;

    for (int nd = 1; nd <= 3; ++ nd)
    {
        for (int k = 1; k <= nk; ++ k)
        {
            for (int j = 1; j <= nj; ++ j)
            {
                for (int i = 1; i <= ni; ++ i)
                {
                    kxyz(1) = faceNormalX(i, j, k, nd);
                    kxyz(2) = faceNormalY(i, j, k, nd);
                    kxyz(3) = faceNormalZ(i, j, k, nd);

                    kxyz = matrix * kxyz;

                    faceNormalX(i, j, k, nd) = kxyz(1);
                    faceNormalY(i, j, k, nd) = kxyz(2);
                    faceNormalZ(i, j, k, nd) = kxyz(3);
                }
            }
        }
    }
    return;
}

PHMatrix< RDouble > PHRigidBodyMovement::CalculateGridCellFaceVelocity(Grid *gridIn, int i, int j, int k, int nd)
{
    StructGrid *grid = static_cast<StructGrid*>(gridIn);
    //! Calculate grid cell face velocity.
    PHSuperMatrix xyzc = GetCellFaceCoordinate(grid, i, j, k, nd);
    PHSuperMatrix currentAngelRate(3);

    PHSuperMatrix & eulerAngleRate = *this->eulerAngleRate;

    for (int ii = 1; ii <= 3; ++ii)
    {
        xyzc(ii) = (xyzc(ii) - (* massCenterOfAircraft)(ii)) / lengthRef;
        currentAngelRate(ii) = eulerAngleRate(ii, 1);
    }
    return CrossProductVector(currentAngelRate, xyzc);
}

PHMatrix< RDouble > PHRigidBodyMovement::GetCellFaceCoordinate(Grid *gridIn, int i, int j, int k, int nd)
{
    StructGrid *grid = static_cast<StructGrid*>(gridIn);
    RDouble3D &structuredX = *(grid->GetStructX());
    RDouble3D &structuredY = *(grid->GetStructY());
    RDouble3D &structuredZ = *(grid->GetStructZ());

    PHMatrix< RDouble > xyzc(3);

    if (nd == 1) 
    {
        xyzc(1) = (  structuredX(i  , j    , k) + structuredX(i  , j    , k + 1)   
                      + structuredX(i  , j + 1, k) + structuredX(i  , j + 1, k + 1) ) * 0.25 ;
        xyzc(2) = (  structuredY(i  , j    , k) + structuredY(i  , j    , k + 1)   
                      + structuredY(i  , j + 1, k) + structuredY(i  , j + 1, k + 1) ) * 0.25 ;
        xyzc(3) = (  structuredZ(i  , j    , k) + structuredZ(i  , j    , k + 1)   
                      + structuredZ(i  , j + 1, k) + structuredZ(i  , j + 1, k + 1) ) * 0.25 ;
    }
    else if (nd == 2) 
    {
        xyzc(1) = (  structuredX(i    , j  , k) + structuredX(i    , j  , k + 1)   
                      + structuredX(i + 1, j  , k) + structuredX(i + 1, j  , k + 1) ) * 0.25 ;
        xyzc(2) = (  structuredY(i    , j  , k) + structuredY(i    , j  , k + 1)   
                      + structuredY(i + 1, j  , k) + structuredY(i + 1, j  , k + 1) ) * 0.25 ;
        xyzc(3) = (  structuredZ(i    , j  , k) + structuredZ(i    , j  , k + 1)   
                      + structuredZ(i + 1, j  , k) + structuredZ(i + 1, j  , k + 1) ) * 0.25 ;
    }
    else
    {
        xyzc(1) = (  structuredX(i    , j  , k) + structuredX(i    , j + 1, k )   
                      + structuredX(i + 1, j  , k) + structuredX(i + 1, j + 1, k ) ) * 0.25 ;
        xyzc(2) = (  structuredY(i    , j  , k) + structuredY(i    , j + 1, k )   
                      + structuredY(i + 1, j  , k) + structuredY(i + 1, j + 1, k ) ) * 0.25 ;
        xyzc(3) = (  structuredZ(i    , j  , k) + structuredZ(i    , j + 1, k )   
                      + structuredZ(i + 1, j  , k) + structuredZ(i + 1, j + 1, k ) ) * 0.25 ;
    }
    return xyzc;
}

void PHRigidBodyMovement::InitializeSixDof()
{
    this->InitializeAircraftProperty();
    this->NondimensionalOf6DOF();
    this->InitializeEulerAngleInformation();
    this->CalculationOfMatrix();
}

void PHRigidBodyMovement::MainCallOf6DOF()
{
    //! six DOF module main functions.
    RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

    currentMarchTimeStep = physicalTimeStep;

    this->AerodynamicsInCFDSolver();

    this->UpdateDynamicInfo6DOF();

    int nZones = PHMPI::GetNumberofGlobalZones();

    int gridLevelIndex = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Grid *grid = PHSPACE::GetGrid(iZone, gridLevelIndex);
        if (! grid) continue;
        this->UpdateSingleZoneGridInformation(grid);
    }
}

void InitializeSixDof()
{
    RigidBodyMovement.InitializeSixDof();
}

void MainCallOf6DOF()
{
    RigidBodyMovement.MainCallOf6DOF();
}

}


