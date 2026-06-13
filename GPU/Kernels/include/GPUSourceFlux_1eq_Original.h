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
//! @file      GPUSourceFlux_1eq_Original.h
//! @brief     Source flux.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
namespace GPUKernels
{
    void CallGPUSpecTurbInit(const int nTotal, const int n_turb);

    __global__ void GPUSpecTurbInit(const int nTotal, const int n_turb, RFloat *spec_turb);
    void CallGPUSourceFlux_1eq_Original(const int nTotalCell, const int nBoundFace, const int isplt, const double SMALL,
                                        const double oreynolds, const RFloat cv13, const RFloat sac2, const RFloat sac3,
                                        const RFloat cw2, const RFloat cw36, const RFloat or6, const RFloat cb1,
                                        const RFloat cb2s, const RFloat cw1k, const RFloat sdilim,
                                        const bool IfModOrgSaMode, const RFloat rkap2, const RFloat cw1);
    __global__ void GPUSourceFlux_1eq_Original(const int nTotalCell, const int nBoundFace, const int isplt,
                                               const double SMALL, const double oreynolds, const RFloat *dVdx,
                                               const RFloat *dVdy, const RFloat *dVdz, const RFloat *q_turb,
                                               const RFloat *q, const RFloat *dqdx, const RFloat *dqdy,
                                               const RFloat *dqdz, const RFloat *visl, const RDouble *vol,
                                               const RDouble *lengthScale, const RFloat cv13, const RFloat sac2,
                                               const RFloat sac3, const RFloat cw2, const RFloat cw36, const RFloat or6,
                                               const RFloat cb1, const RFloat cb2s, const RFloat cw1k,
                                               const RFloat sdilim, const bool IfModOrgSaMode, const RFloat rkap2,
                                               const RFloat cw1, RFloat *spec_turb, RFloat *res_turb);
} //！ namespace GPUKernels
