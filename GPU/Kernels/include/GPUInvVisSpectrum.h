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
//! @file      GPUInvVisSpectrum.h
//! @brief     Inviscid spectral radius.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
#include "Param_NSSolverUnstruct.h"
namespace GPUKernels
{
    void CallGPUSetFieldInvSpectrum(int nTotalCell, const RFloat value);
    void CallGPUSetFieldVisSpectrum(int nTotalCell, const RFloat value);
    //void CallGPUInvSpectrumRadiusCalculate(const int nTotalCell, const int nTotal, const int nTotalFace);
    void CallGPUSpectrumRadiusInviscid(Grid *gridIn);
    //void CallGPUVisSpectrumRadiusCalculate(const int nTotalCell, const int nTotal, const int nTotalFace, const double refReNumber, const double SMALL);
    void            CallGPUSpectrumRadiusViscous(Grid *gridIn, Param_NSSolverUnstruct *parameters);
    __global__ void GPUSetFieldInvSpectrum(const int nTotalCell, RFloat *invSpectrum, const RFloat value);
    __global__ void GPUSetFieldVisSpectrum(const int nTotalCell, RFloat *visSpectrum, const RFloat value);
    __global__ void GPUInvSpectrumRadiusCalculate(const int nTotalCell, const int nTotal, const int nTotalFace,
                                                  const int *left_cell_of_face, int *right_cell_of_face,
                                                  const RDouble *xfn, const RDouble *yfn, const RDouble *zfn,
                                                  const RDouble *area, const RFloat *q, const RFloat *gama,
                                                  RFloat *invSpectrum);
    __global__ void GPUVisSpectrumRadiusCalculate(const int nTotalCell, const int nTotal, const int nTotalFace,
                                                  const double refReNumber, const double SMALL,
                                                  const int *left_cell_of_face, const int *right_cell_of_face,
                                                  const RDouble *xcc, const RDouble *ycc, const RDouble *zcc,
                                                  const RDouble *xfn, const RDouble *yfn, const RDouble *zfn,
                                                  const RDouble *area, const RDouble *q, const RFloat *visl,
                                                  const RFloat *vist, RFloat *visSpectrumRadius);
} //！ namespace GPUKernels
