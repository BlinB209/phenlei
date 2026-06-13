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
//! @file      SixDOF.h
//! @brief     Explain this file briefly.
//! @author    xxx.

#pragma once
#include "Constants.h"
#include "Geo_StructGrid.h"
#include "PHSuperMatrix.h"
#include "Force.h"

namespace PHSPACE
{

class PHRigidBodyMovement
{
public:
    typedef PHMatrix< RDouble > PHSuperMatrix;
protected:
    //! variables declaration.
    //! the attributes of aircraft.
    RDouble massOfAircraft;                           //! the mass of aircraft.
    PHSuperMatrix *massCenterOfAircraft;              //! the centroid position of aircraft.
    PHSuperMatrix *inertiaTensorOfAircraft;           //! the moment of inertia of aircraft.
    PHSuperMatrix *reverseOfInertiaTensor;            //! the inverse matrix of moment of inertia

    //! the data changed with the flow solvers.
    //! input£ºforce and moment coefficient.
    PHSuperMatrix *aerodynamicForce;                  //! aerodynamic force.
    PHSuperMatrix *aerodynamicMoment;                 //! aerodynamic moment.
    PHSuperMatrix *otherForceSum;                     //! sum of other forces.(gravity, thrust.)
    PHSuperMatrix *otherMomentSum;                    //! sum of other moment.
    //! input£ºthe physical nondimensional time step of flow solvers.
    RDouble currentMarchTimeStep;                     //! marching time step.
    //! input£ºfreestream AOA, slide angle.
    RDouble initiativeAngleOfAttack;                  //! init angle of attack.
    RDouble initiativeAngleOfSlide;                   //! init angle of slide.
    //! output£ºtemporary angle of attack, temporary angle of slide.
    RDouble temporaryAngleOfAttack;                   //! temporary angle of attack.
    RDouble temporaryAngleOfSlide;                    //! temporary angle of slide.
    //! input/output£ºmesh information.
    StructGrid *grid;
    //! int ni, nj, nk;
    //! changed coordinate data.
    PHSuperMatrix *eulerAngle;                        //! Euler angle,
    PHSuperMatrix *eulerAngleRate;                    //! euler angle rate in time n-2, n-1 and n.
    PHSuperMatrix *eulerAngleAcceleration;            //! euler angle acceleration.

    PHSuperMatrix *transformMatrixCurrent;            //! the current time system. -> transform matrix of body axis.
    PHSuperMatrix *transformMatrixLastTime;           //! the last time system. -> transform matrix of body axis.
    //! the amount of motion.
    PHSuperMatrix *momentOfInertia;                   //! moment Of Inertia in time n-2, n-1 and n.
    PHSuperMatrix *bodyAngleRate;                     //! body angle rate
    PHSuperMatrix *bodyAngleAcceleration;

    //! reference values.
    //! input£ºforce and moment coefficient.
    RDouble densityRef;                               //! reference desity.
    RDouble areaRef;                                  //! reference area.
    RDouble lengthRef;                                //! reference length.
public:
    PHRigidBodyMovement();
    ~PHRigidBodyMovement();
protected:
    void UpdateCoordinate               (Grid *gridIn);
    void UpdateGridVelocity             (Grid *gridIn);
    void InitializeAircraftProperty     ();
    void InitializeEulerAngleInformation();
    void InitializeInertia      ();
    void CalculationOfMatrix    ();
    void NondimensionalOf6DOF   ();
    void GetReferenceValue      ();
    void AerodynamicsInCFDSolver();
    void UpdateDynamicInfo6DOF();
    void UpdateMomentOfInertia();
    void UpdateOfInertia();
    void CalculateAngleRateOfBody();
    void CalculateEulerAngleRate();
    void CalculateEulerAngle();
    void CalculateEulerAngleAcceleration();
    void CalculateAngleOfAttackAndSlide();
    PHSuperMatrix CalculateGridCellFaceVelocity(Grid *gridIn, int i, int j, int k, int nd);
    PHSuperMatrix GetCellFaceCoordinate(Grid *gridIn, int i, int j, int k, int nd);
    void UpdateCellArea(Grid *gridIn);
    void UpdateSingleZoneGridInformation(Grid *gridIn);
public:
    void InitializeSixDof();
    void MainCallOf6DOF();
};

void InitializeSixDof();
void MainCallOf6DOF();

}

