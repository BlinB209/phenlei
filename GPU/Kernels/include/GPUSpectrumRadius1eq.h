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
//! @file      GPUSpectrumRadius1eq.h
//! @brief     Spectrum radius for turbulent solver.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
namespace GPUKernels
{
    void CallGPUSpecTurbCell(const int nTotalCell, const int nBoundFace, const RFloat dualTimeSpectrumC1,
                             const RFloat dualTimeSpectrumC2, const RFloat cflturb, const double SMALL);

    __global__ void GPUSpecTurbCell(const int nTotalCell, const int nTotal, const RFloat dualTimeSpectrumC1,
                                    const RFloat dualTimeSpectrumC2, const RFloat cflturb, const double SMALL,
                                    const RDouble *vol, const RFloat *dt, RFloat *spec_turb);

    void CallGPUSpecTurbMatTurbFaces(const int nTotalFace, const int nBoundFace, const int nTotalCell,
                                     const double SMALL, const RFloat oreynolds, const double osigma, const double cb2,
                                     const bool IfModOrgSaMode);

    __global__ void GPUSpecTurbMatTurbFaces(const int nTotalFace, const int nBoundFace, const int nTotalCell,
                                            const double SMALL, const RFloat oreynolds, const double osigma,
                                            const bool IfModOrgSaMode, const int *left_cell_of_face,
                                            const int *right_cell_of_face, const RDouble *xfn, const RDouble *yfn,
                                            const RDouble *zfn, const RDouble *area, const RFloat *q,
                                            const RFloat *gama, const RFloat *q_turb, const RFloat *visl,
                                            const RDouble *vol, const RFloat cb2, RFloat *spec_turb, RFloat *mat_turbl,
                                            RFloat *mat_turbr);
} //! namespace GPUKernels
