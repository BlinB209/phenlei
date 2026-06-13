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
//! @file      MatrixLUSGS.h
//! @brief     Implicit iterative method for Linearized N-S equations based on exact Jacobian matrix.
//! @author    LiuJian, Yang Yufeng.

#include "Geo_StructGrid.h"
#include "FieldProxy.h"
#include "FaceProxy.h"
#include "PHMatrix.h"

namespace PHSPACE
{

//! To execute the forward sweep of the Matrix LU-SGS method.
//! @param[in ]: gridIn denotes the current computational regions of the grids.
//! @param[out]: dqProxy is an array that store the difference values between the conservative variables at n and n+1 time step.
void SolveMatrixLUSGSForwardSweep(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);

//! To execute the backward sweep of the Matrix LU-SGS method.
//! @param[in ]: gridIn denotes the current computational regions of the grids.
//! @param[out]: dqProxy and dqStarProxy are arrays that store the difference values between the conservative variables at n and n+1 time step.
void SolveMatrixLUSGSBackwardSweep(Grid *gridIn, FieldProxy *dq_proxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);

//! Compute and save the diagonal matrix on grids that are used for matrix LU-SGS method.
//! @param[in ]: gridIn denotes the current computational regions of the grids.
void SolveDiagonalMatrix(StructGrid *gridIn);

//! Construct the Jacobian matrix of the inviscid flux item.
//! @param[in ]: primitiveVars is the array that stores the primitive variables of the current cell.
//! @param[in ]: temperature is the array that stores the temperature values of the current cell.
//! @param[in ]: xfn, yfn and zfn are components of the average normal vector located in the center of the current cell.
//! @param[in ]: area denotes the avarage area located in the center of the current cell.
//! @param[in ]: vgn denotes the moment velocity located in the center of the current cell, and it is used for dynamic mesh.
//! @param[in ]: gama denotes the specific heat ratio.
//! @param[out]: tempMatrix is an matrix that stores the Jacobian matrix of the inviscid flux item.
void ComputeInviscidMatrix(RDouble *primitiveVars, RDouble xfn, RDouble yfn, RDouble zfn, RDouble area, RDouble vgn, RDouble gama, RDouble2D *tempMatrix);

//! Construct the roe matrix of the inviscid flux item.
//! @param[in ]: priml is the array that stores the primitive variables of the left cell.
//! @param[in ]: primr is the array that stores the primitive variables of the current cell.
//! @param[in ]: xfn, yfn and zfn are components of the average normal vector located in the center of the current cell.
//! @param[in ]: gama denotes the specific heat ratio.
//! @param[in ]: area denotes the avarage area located in the center of the current cell.
//! @param[in ]: rtemL and rtemR are the parameters of Entropy fix, add by liujian 2019 From JCP ,1995, Lin H-C
//! @param[out]: jacobianMatrix is an matrix that stores the Roe matrix of the inviscid flux item.
void SolveRoeMatrix(RDouble *priml, RDouble *primr, RDouble xfn, RDouble yfn, RDouble zfn, RDouble gama, RDouble area, RDouble rtemL, RDouble rtemR, RDouble2D *jacobianMatrix);

//! Construct the Jacobian matrix of the viscous flux item.
//! @param[in ]: primitiveVars is the array that stores the primitive variables of the current cell.
//! @param[in ]: temperature is the array that stores the temperature values of the current cell.
//! @param[in ]: xfn, yfn and zfn are components of the average normal vector located in the center of the current cell.
//! @param[in ]: area denotes the avarage area located in the center of the current cell.
//! @param[in ]: vol denotes the avarage volume located in the center of the current cell.
//! @param[in ]: vgn denotes the moment velocity located in the center of the current cell, and it is used for dynamic mesh.
//! @param[in ]: gama denotes the specific heat ratio.
//! @param[in ]: laminarViscosity and turbulentViscosity denote the laminar viscosity and turbulent viscosity, respectively.
//! @param[out]: tempMatrix is an matrix that stores the Jacobian matrix of the viscous flux item.
void ComputeViscousMatrix(RDouble *primitiveVars, RDouble *temperature, RDouble xfn, RDouble yfn, RDouble zfn, RDouble area, RDouble vol, 
                          RDouble vgn, RDouble gama, RDouble laminarViscosity, RDouble turbulentViscosity, RDouble2D *tempMatrix);

//! Compute and save the inverse matrix of the diagonal matrix by LU decomposition.
//! @param[in ]: diagonalMatrixD is the diagonal matrix used for the matrix LU-SGS method.
//! @param[in ]: nEquation is the 
//! @param[out]: inverseDiaMatrix is the inverse matrix of the diagonal matrix.
void ObtainInverseDiagonalMatrix(RDouble2D *diagonalMatrixD, int nEquation, RDouble2D *inverseDiaMatrix);

}