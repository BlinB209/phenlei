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
//! @file      TemporaryOperations.h
//! @brief     Copy functions for device.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "Geo_Grid.h"
#include "Limiter.h"
#include "Gradient.h"
#include "GPUTestFunctions.h"

using namespace PHSPACE;

namespace TemporaryOperations
{
    void GPUQCopy(Grid *grid_in, const int neqn, const int nTotal);
    void GPUTurbUnsteadyCopy(Grid *grid_in);
    void GPUNSUnsteadyCopy(Grid *grid_in);
    void GPUQCopy(Grid *grid_in);
    void GPUQ_NSCopy(RFloat **q, const int neqn_ns, const int nTotal);
    void GPUQl_NSCopy(RFloat **ql, const int neqn_ns, const int len_ql);
    void GPUQr_NSCopy(RFloat **qr, const int neqn_ns, const int len_qr);
    void GPUQ_NSCopyBack(RFloat **q, const int neqn_ns, const int nTotal);
    void CallGPUFlowVariablesCopyBack(Grid *grid_in);
    void GPUGamaCopy(Grid *grid_in, const int nTotal);

    RFloat **GPUQl_NSCopyBack(RFloat *d_ql_ns, const int neqn_ns, const int len_ql);
    RFloat **GPUQr_NSCopyBack(RFloat *d_qr_ns, const int neqn_ns, const int len_qr);
    RFloat  *GPUGamaL_NSCopyBack(RFloat *d_gamaL_ns, const int len_gamaL);
    RFloat  *GPUGamaR_NSCopyBack(RFloat *d_gamaR_ns, const int len_gamaR);

    //! copy limiter onto device
    void GPULimitCopy(Limiter *limiter, const int neqn, const int nTotal);
    //! copy LIMIT onto device
    void GPULIMITCopy(RFloat **LIMIT, const int neqn, const int nTotal);
    //! copy dqdx ... onto device
    void GPUDqdxDqDyDqDzCopy(Gradient *gradient, const int neqn, const int nTotal);
    //! copy xtn .. onto device
    void GPUXtnYtnZtnCopy(Grid *grid_in, const int nTotalFace);
    //! copy vgn ont to device
    void GPUVgnCopy(Grid *grid_in, const int nTotalFace);
    //! copy res onto device
    void GPUResCopy(Grid *grid_in);
    //! copy epsCell onto device
    void GPUepsCellCopy(RFloat *epsCell, const int nTotalCell);
    void GPUQlQrCopyBack(FaceProxy *face_proxy, const int SEG_LEN);
    void GPUGamaLRCopyBack(FaceProxy *face_proxy, const int SEG_LEN);
    void GPUFluxCopyBack(FaceProxy *face_proxy, const int SEG_LEN);
    //! copy res back to cpu
    void GPUResCopyBack(Grid *grid_in, const int nTotal);
    //! copy d_epsCell back to device
    void GPUepsCellCopyBack(RFloat *epsCell, RFloat *d_epsCell, const int nTotalCell);
    //! copy dmax dmin onto device
    void GPUDmaxDminCopy(RFloat *dmax, RFloat *dmin, const int nTotal);
    void GPUQCopy(RFloat **q, RFloat *d_q, const int neqn, const int nTotal);
    void GPUQlQrTurbCopyBack(FaceProxy *face_proxy, const int n_turb, const int SEG_LEN);
    void GPUFluxTurbCopyBack(RFloat **flux, RFloat *d_flux, const int neqn, const int SEG_LEN);
    void GPUFluxSubTurbCopy(RFloat **flux_sub, RFloat *d_flux_sub, const int neqn, const int SEG_LEN);
    void GPUVislCopy(RFloat *visl, RFloat *d_visl, const int nTotal);
    void GPUVislCopy(RFloat **res_turb, RFloat *d_res_turb, const int neqn, const int nTotal);
    void GPUResTurbCopyBack(RFloat **res_turb, RFloat *d_res_turb, const int neqn, const int nTotal);
    void GPUGamaCopyBack(Grid *grid_in);
    void GPUTCopyBack(Grid *grid_in);
    //!void GPUPrim_infCopy(RFloat *prim_inf, const int neqn);
    void GPUPrim_infCopy();
    void GPUQ_proxyCopyBack(RFloat **q, const int m, const int n);
    void GPUQ_proxy_tmpCopyBack(RFloat **q, const int m, const int n);
    void GPUrhsCopy(RFloat **rhs, const int neqn, const int nTotal);
    void GPUrhsCopyBack(RFloat **rhs, const int neqn, const int nTotal);
    void GPUDtCopy(RFloat *dt, const int nTotal);
    void GPUResCopyBack(Grid *grid_in);
    void GPUInvSpectrumRadiusCopyBack(RFloat *invSpectrumRadius, const int nTotalCell);
    void GPUInvSpectrumRadiusCopy(RFloat *invSpectrumRadius, const int nTotalCell);
    void GPUVisSpectrumRadiusCopyBack(RFloat *visSpectrumRadius, const int nTotalCell);
    void GPUVisSpectrumRadiusCopy(RFloat *visSpectrumRadius, const int nTotalCell);
    void GPUDtCopyBack(RFloat *dt, const int nTotal);
    void GPUDtvCopyBack(RFloat *dtv, const int nTotalCell);
    void GPUMDtCopyBack(RFloat &minDt, RFloat &maxDt);
    void GPUVislCopyBack(RFloat *visl, const int nTotal);
    void GPUQ_TurbCopy(RFloat **q, const int m, const int n);
    void GPUVislCopy(RFloat *visl, const int nTotal);
    void GPUVistCopy(RFloat *vist, const int nTotal);
    void GPUMVistCopyBack(RFloat &vistmin, RFloat &vistmax);
    void GPUVistCopyBack(RFloat *vist, const int nTotal);
    void GPUDQNSCopy(RFloat **dq, const int nl, const int nTotal);
    void GPUDQNSCopyBack(RFloat **dq, const int nl, const int nTotal);
    void GPUResTurbCopyBack(Grid *grid_in);
    void GPUDQ_turbCopy(RFloat **dq_turb, const int nTotal, const int n_turb);
}
