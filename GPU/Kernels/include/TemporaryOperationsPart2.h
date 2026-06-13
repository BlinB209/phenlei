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
//! @file      TemporaryOperationsPart2.h
//! @brief     Copy functions for device.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include <iostream>
#include "Precision.h"
#include "cudaErrorHandle.h"

using namespace PHSPACE;
namespace TemporaryOperations
{
    void GPUGradientVarsCopy(string q_name, int nTotal, int nTVar, RFloat **q2d, RFloat **dqdx2d, RFloat **dqdy2d,
                             RFloat **dqdz2d);
    void GPUViscousCoefCopy(const RFloat *h_visl, const RFloat *h_vist, const int nTotal);
    void GPUResCopy(RFloat **res, const int nl, const int nTotal);
    void GPUResCopyBack(RFloat **res, const int nl, const int nTotal);
    void GPUGradientCopyBack(const string q_name, const int nTotal, const int nTVar, RFloat **dqdx2d, RFloat **dqdy2d,
                             RFloat **dqdz2d);
    void GPUStoreBoundGradCopyBack(const int nBoundFace, const int nTVar, RFloat **bdqx, RFloat **bdqy, RFloat **bdqz);
    void GPUQ_turbCopy(RFloat **q_turb, const int nTotal, const int n_turb);
    void GPUFluxTurbCopy(RFloat **flux, RFloat **flux_sub, const int n_turb);
    void GPUResTurbCopy(RFloat **res_turb, const int nTotal, const int n_turb);
    void GPUResTurbCopyBack(RFloat **res_turb, const int n_turb, const int nTotal);
    void GPUSpecTurbCopy(const int nTotalCell, const int nBoundFace, const int n_turb, RFloat **spec_turb);
    void GPUDtCopy(const int nTotalCell, const int nBoundFace, const RFloat *dt);
    void GPUSpecTurbCopyBack(RFloat **spec_turb, const int n_turb, const int nTotal);
    void GPUMatTurbCopyBack(RFloat **mat_turbl, RFloat **mat_turbr, const int n_turb, const int nTotalFace);
    void GPUQ_turbCopyBack(RFloat **q_turb, const int nTotal, const int n_turb);
    void GPUTurbRhsCopy(RFloat **rhs, const int nTotal, const int n_turb);
    void GPUQProxyTempTurbCopyBack(RFloat **field1, const int nTotal, const int neqn);
    void GPUDQ_turbCopyBack(RFloat **dq_turb, const int nTotal, const int n_turb);
    void GPUSpecCopyBack(RFloat *spec, const int nTotal);
} //! namespace TemporaryOperations
