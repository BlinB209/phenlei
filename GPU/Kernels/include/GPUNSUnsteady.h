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
//! @file      GPUNSUnsteady.h
//! @brief     Unsteady functions.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "BasicDeviceVariables.h"
namespace GPUKernels
{
    void CallGPUNSUnsteadyInitSpec(const int nTotal, const RFloat dualTimeSpectrumC1, const RFloat dualTimeSpectrumC2);
    __global__ void GPUNSUnsteadyInitSpec(const int nTotal, const RFloat dualTimeSpectrumC1,
                                          const RFloat dualTimeSpectrumC2, const double *vol, const RFloat *dt,
                                          RFloat *spec);

    void CallGPUNSUpdateUnsteadyFlow(const int nTotal, const int nEquation, const int ifStaticsFlowField);

    __global__ void GPUNSUpdateUnsteadyFlow(const int nTotal, const int nEquation, const RFloat *q,
                                            const RFloat *resTmp, RFloat *qn1, RFloat *qn2, RFloat *resn1,
                                            RFloat *resn2);
    void CallGPUNSUnsteadyConvergence(const int nTotal, const int nTotalCell, const int nEquation, const int nchem,
                                      const int nm, const int nl, const int ntmodel, const double refGama);
    __global__ void setSumAsZero(RFloat *sum1, RFloat *sum2);
    __global__ void GPUNSUnsteadyConvergence(const int nTotal, const int nTotalCell, const int nEquation,
                                             const int nchem, const int nm, const int nl, const int ntmodel,
                                             const double refGama, const RFloat *q, const RFloat *qn1,
                                             const RFloat *qn2, const RFloat *res, RFloat *sum1, RFloat *sum2);
    void CallGPUNSDualTimeSourceRes(const int nTotal, const int nTotalCell, const int nEquation, const int nchem,
                                    const int nm, const int nl, const int ntmodel, const RDouble refGama,
                                    const RFloat dualTimeResC1, const RFloat dualTimeResC2, const RFloat dualTimeResC3,
                                    const RFloat dualTimeQC1, const RFloat dualTimeQC2, const RFloat dualTimeQC3);
    __global__ void GPUNSDualTimeSourceRes(const int nTotal, const int nTotalCell, const int nEquation, const int nchem,
                                           const int nm, const int nl, const int ntmodel, const RDouble refGama,
                                           const RFloat dualTimeResC1, const RFloat dualTimeResC2,
                                           const RFloat dualTimeResC3, const RFloat dualTimeQC1,
                                           const RFloat dualTimeQC2, const RFloat dualTimeQC3, const RDouble *vol,
                                           const RFloat *q, const RFloat *qn1, const RFloat *qn2, const RFloat *resn1,
                                           const RFloat *resn2, RFloat *res);

    void CallGPUSetNSResTmpbyRes(const int nTotal, const int nTotalCell, const int nEquation);

    __global__ void GPUSetNSResTmpbyRes(const int nTotal, const int nTotalCell, const int nEquation, const RFloat *res,
                                        RFloat *resTmp);
} //! namespace GPUKernels
