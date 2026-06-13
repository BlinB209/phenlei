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
//! @file      Jacobian.h
//! @brief     Explain this file briefly.
//! @author    He Xin, Bell, He Kun, He Xianyao, Liu Jian, Xu Qingxin, etc.

#pragma once
#include "Precision.h"
namespace PHSPACE
{
//! The fuction computes the formula (M +- R) * dQ, M and R denote the Jacobian matrix and spectrum radius of inviscid flux,respectively,\n
//! dQ is the increment of conservative vector.
void MXDQ_STD(RDouble *prim, RDouble nxs, RDouble nys, RDouble nzs, RDouble area, RDouble vgn, RDouble gama, RDouble *dq, RDouble *f, 
              int nm, int nl, RDouble radius, int ipn);

void MXDQ_STD(RDouble *prim, RDouble *gykb, RDouble *dq, RDouble *df, RDouble gama, int nm, int nl, RDouble radius, int ipn);

//! The fuction computes the formula (Mv +- Rv) * dQ, Mv and Rv denote the Jacobian matrix and spectrum radius of viscous flux,respectively,\n
//! dQ is the increment of conservative vector.
void MVXDQ_STD_LP(RDouble *prim, RDouble nxs, RDouble nys, RDouble nzs, RDouble area, RDouble vol, RDouble visl, RDouble vist, RDouble gama, 
                  RDouble *dq, RDouble *f, int nm, int nl, RDouble radius, int ipn);

//! The functions for low speed precondition.
//! Compute conservative precondition matrix.
//! The function computes flux,respectively,\n
void ComputeConservativePreconMatrix(int sign, RDouble **MG, RDouble *prim, RDouble gama, RDouble preconCoeff, int nm, int neqn);
void ConservativePreconMatrixVariables(RDouble *field, RDouble *preconMatrix, int neqn);
void ConservativePreconMatrixVariables(RDouble *field, RDouble *preconMatrix, RDouble timeCoeff, int neqn);

void MatrixFreeDF(RDouble *q, RDouble *dq, RDouble *nxyz, RDouble gama, RDouble *df, int nm, int neqn, int ipn);
void PreconMatrixFreeDF(RDouble *q, RDouble *dq, RDouble *nxyz, RDouble **MG, RDouble gama, RDouble *df, int nm);
void PreconMatrixFreeDF(RDouble *q, RDouble *dq, RDouble *nxyz, RDouble *preconMatrix, RDouble gama, RDouble preconCoeff, RDouble *df, int nm, int neqn, int ipn);
void PreconMatrixFreeDF(RDouble *q, RDouble *dq, RDouble *nxyz, RDouble *preconMatrix, RDouble gama, RDouble preconCoeff, RDouble timeCoeff, RDouble *df, int nm, int neqn, int ipn);
}
