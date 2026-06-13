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
//! @file      AleModel.h
//! @brief     Explain this file briefly.
//! @author    xxx.

#pragma once
#include "Precision.h"

namespace PHSPACE
{
const int LMSM_IMPEUL    = 1;
const int LMSM_IMPTRZ    = 2;
const int LMSM_IMPBD2    = 3;
const int LMSM_IMPTHR    = 4;
const int LMSM_IMPFOR    = 5;
const int LMSM_NODUAL    = 6;

void lmsmcoef(int lmsmflag, int outnstep, RDouble *outerdt, RDouble &xi, RDouble &theta, RDouble &phi);
void ordcoef(int order, RDouble &xi, RDouble &theta, RDouble &phi, RDouble *outerdt, int iact);

class Grid;
class AleParameter
{
public:
    int isUnsteady, isAle;
    int order;
    //! Linear Multi Step Methods flag.
    int lmsmflag; 
    //! Outer iteration step.
    int outnstep;

    RDouble theta;
    RDouble phi;
    RDouble xi;
    RDouble *outerdt;
    RDouble dtau;
    RDouble dtsave;
    RDouble simutime;
    RDouble timemax;
public:
    int & GetLMSMFlag() { return lmsmflag; };
    void SetLMSMFlag();
    void Init();
    AleParameter();
    ~AleParameter();
};

class AleModel
{
public:
    //! Parameters for the Arbitrary Lagrangian Eulerian Method (ALE).
    //! Euler's angles and reference point position at time n, n-1, n-2
    //! and their times derivatives at time n (0)

    //! x,y,z position.
    //! coord0(1,1),coord0(2,1),coord0(3,1)
    //! x,y,z velocity.
    //! coord0(1,2),coord0(2,2),coord0(3,2)
    //! x,y,z accelerations.
    //! coord0(1,3),coord0(2,3),coord0(3,3)
    
    int maxale;
    static RDouble * euler_angle;    //! Rotation euler angle.
    static RDouble **trans_ref;      //! Reference point, velocity, accelerations.
    static RDouble **xyz_ref;        //! Reference point n, n-1, n-2.
    static RDouble **rotate_ref;     //! Reference point, velocity, accelerations.
    static RDouble **angle_ref;    
    static RDouble **matrix;         //! Matrix.
    int * ialepar;
    static AleParameter *param;
public:
    static void Init();             //! Initialization.
    static void Allocate();         //! Allocate memory.
    static void Free();

    static RDouble ** CreateMatrix();
    static void DestroyMatrix(RDouble **matrix);

    static RDouble * CreateVector();
    static void DestroyVector(RDouble *v);

    static void CrossMultiply(RDouble *omega, RDouble *dr, RDouble *vel);

    static void SetMatrix(RDouble ** matrix_x, RDouble ** matrix_y);
    static void TransposeMatrix(RDouble ** matrix);
    static void ComputeRelativeRotateMatrix();
    static void ComputeRelativeRotateMatrix(RDouble ** matrix_now, RDouble ** matrix_old);

    static void MatrixMultiply(RDouble **x, RDouble **y, RDouble **z);
    static void MatrixRotate  (RDouble *angle, RDouble **matrix);
    static void MatrixRotateX (RDouble  angle, RDouble **matrix);
    static void MatrixRotateY (RDouble  angle, RDouble **matrix);
    static void MatrixRotateZ (RDouble  angle, RDouble **matrix);

    static void GetGridFileName(string &filename);

    static void EulerAngle();
    static void DynamicGrid();
    static void GenerateDynamicGrid();
    static void BackUpOldGrid();
    static void ReadDynamicGrid();
    static void ReadDynamicGrid(fstream &file);
    static void ReadGridCoordinate(fstream &file, int iZone);
    static void ComputeGridVelocity();
    static void UpdateGridInfo();
    static void RotateTranslation();
    static void RotateTranslation(Grid *grid);
    static void GridVerticeVelocity       (Grid *grid, RDouble *xt, RDouble *yt, RDouble *zt);
    static void GridVerticeVelocityRigid  (Grid *grid, RDouble *xt, RDouble *yt, RDouble *zt);
    static void GridVerticeVelocityGeneral(Grid *grid, RDouble *xt, RDouble *yt, RDouble *zt);
    static void GridSurfaceVelocity       (Grid *grid, RDouble *xt, RDouble *yt, RDouble *zt);
    static void CheckMovement();
    static void Deform();
    static void MoveGrid();
    static void LinearMultiStepMethodsCoeffcient(AleParameter *param);
    static void GetOuterdt(AleParameter *param);
    static void UnsteadyPreSolve();
    static void GetOuterSrc();
};

void chckmovmt(int *ialepar, RDouble **angle0, RDouble **coord0, RDouble **angle, RDouble **coord,RDouble *outerdt, 
               int lmsmflag, int outnstep);
//void deftrans(RDouble **coord, RDouble time, RDouble *outerdt, int *ialepar, RDouble *ralepar);
void chckdeform();

}

