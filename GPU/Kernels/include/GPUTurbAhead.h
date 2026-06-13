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
//! @file      GPUTurbAhead.h
//! @brief     Turbulent functions for turbulent solver.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
namespace GPUKernels
{
    void CallGPUTurbZeroResiduals(const int numberOfGeneralCells, const int numberOfEquations);
    void CallGPUTurbStoreRhsByResidual(const int nTotal, const int n_turb);
    void CallGPUTurbLoadQ(const int nTotalCell, const int nTotal, const int n_turb);
    void CallGPUTurbFillField(const int nTotal, const int neqn);
    void CallGPUTurbLoadResiduals(const int nTotal, const int n_turb);
    void CallGPUTurbLhs(const int nTotalCell, const int nTotal, const int n_turb, const double coef);
    void CallGPUTurbUpdateFlowField(const int nTotalCell, const int nTotal, const int n_turb, const double coef,
                                    const double turb_relax, const int nnegtive_max, const int ZoneID,
                                    const double SMALL, const double cv13, const double mytmax, const double my_nu_max);
    void CallGPUTurbUpdateLUSGSFlowField(const int nTotalCell, const int nTotal, const int n_turb, const double coef,
                                         const double turb_relax, const int nnegtive_max, const int ZoneID,
                                         const double SMALL, const double cv13, const double mytmax,
                                         const double my_nu_max);
    void CallGPURecoverResidual(const int nTotal, const int n_turb, const int a);
    void CallGPUTurbViscosity(const int nTotal, const int nTotalCell, const double mytmax, const double turb_cv13,
                              const double SMALL);
    __global__ void GPUTurbZeroResiduals(const int numberOfGeneralCells, const int numberOfEquations, RFloat *res);
    __global__ void GPUTurbStoreRhsByResidual(const int nTotal, const int n_turb, const RFloat *res, RFloat *rhs);
    __global__ void GPUTurbLoadQ(const int nTotalCell, const int nTotal, const int n_turb, const RFloat *q_turb,
                                 RFloat *qq_turb);
    __global__ void GPUTurbFillField(const int nTotal, const int neqn, RFloat *field1, const RFloat *field2);
    __global__ void GPUTurbLoadResiduals(const int nTotal, const int neqn, RFloat *res, const RFloat *rhs);
    __global__ void GPUTurbLhs(const int nTotalCell, const int nTotal, const int n_turb, const double coef,
                               const RFloat *dt, RFloat *dq_turb);
    __global__ void
    GPUTurbUpdateFlowField(const int nTotalCell, const int nTotal, const int n_turb, const int *cell2face,
                           const int *posiCell2Face, const int *face_number_of_each_cell, const int *left_cell_of_face,
                           const int *right_cell_of_face, const RFloat *q_turb, const double coef, const double *turboo,
                           const RFloat *visl, const double turb_relax, const int nnegtive_max, const int ZoneID,
                           const double SMALL, const double cv13, const double mytmax, const double my_nu_max,
                           int *n_neg, int *n_pos, const RFloat *t, const RFloat *qpmv, const RFloat *dq, RFloat *q);
    __global__ void GPUSetNegPosZero(int *n_neg, int *n_pos);
    __device__ void GPUSMoothTurbPoint(const int nTotalCell, const int nTotal, const int n_turb, const int *cell2face,
                                       const int *posiCell2Face, const int *left_cell_of_face,
                                       const int *right_cell_of_face, const int *face_number_of_each_cell,
                                       const RFloat *q_turb, const int i, RFloat *prim_var);
    __global__ void GPURecoverResidual(const int nTotal, const int n_turb, RFloat *res, const RFloat *rhs);
    __global__ void GPUTurbViscosity(const int nTotal, const int nTotalCell, const double mytmax,
                                     const double turb_cv13, const double SMALL, const RFloat *q, const RFloat *q_turb,
                                     const RFloat *visl, RFloat *vist, double *vistmax, double *vistmin);
    __global__ void GPUTurbViscosity_S1(const int nTotal, const int nTotalCell, const double mytmax,
                                        const double turb_cv13, const double SMALL, const RFloat *q,
                                        const RFloat *q_turb, const RFloat *visl, RFloat *vist);
    __global__ void GPUTurbUpdateFlowField_S1(const int nTotalCell, const int nTotal, const int n_turb,
                                              const int *cell2face, const int *posiCell2Face,
                                              const int *face_number_of_each_cell, const int *left_cell_of_face,
                                              const int *right_cell_of_face, const RFloat *q_turb, const double coef,
                                              const double *turboo, const RFloat *visl, const double turb_relax,
                                              const int nnegtive_max, const int ZoneID, const double SMALL,
                                              const double cv13, const double mytmax, const double my_nu_max,
                                              int *n_neg, int *n_pos, const RFloat *t, const RFloat *qpmv,
                                              const RFloat *dq, RFloat *q);
} //! namespace GPUKernels
