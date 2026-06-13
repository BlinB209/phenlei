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
//! @file      GPUInviscidFlux.h
//! @brief     Inviscid flux.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "Geo_Grid.h"
#include "NSSolver.h"
#include "FaceProxy.h"
#include "GeomProxy.h"
#include "Param_NSSolverUnstruct.h"
#include "Limiter.h"

using namespace PHSPACE;

namespace GPUNSSolverUnstruct
{
    //! call function for kernel GPUGetQlQr
    void CallGPUGetQlQr(Grid *grid_in, const int nst, const int ned, const int SEG_LEN);
    //! GetQlQr GPU kernel
    __global__ void GPUGetQlQr(int *d_left_cell_of_face, int *d_right_cell_of_face, RFloat *d_q_ns, RFloat *d_ql_ns,
                               RFloat *d_qr_ns, const int nst, const int ned, const int d_neqn, const int d_nTotal,
                               const int d_SEG_LEN);
    //! call function for kernel GPUGetGamaLR
    void CallGPUGetGamaLR(const int nst, const int ned);
    //! GetQlQr GPU kernel
    __global__ void GPUGetGamaLR(int *d_left_cell_of_face, int *d_right_cell_of_face, RFloat *d_gama_ns,
                                 RFloat *d_gamaL_ns, RFloat *d_gamaR_ns, const int nst, const int ned);
    //! call function for kernel GPUReConstructFaceValueLoop1
    void CallGPUReConstructFaceValue(Grid *gridIn, Limiter *limiter, int localStart, int localEnd, const int SEG_LEN,
                                     Param_NSSolverUnstruct *parameters);
    void CallGPUReConstructFaceValueLoop1(Grid *grid_in, const int limit_vec, int nst, int ned, const int SEG_LEN);
    //! GPUReConstructFaceValue kernel for loop1
    __global__ void GPUReConstructFaceValueLoop1(RDouble *d_xfc, RDouble *d_yfc, RDouble *d_zfc,
                                                 int *d_left_cell_of_face, int *d_right_cell_of_face, RDouble *d_xcc,
                                                 RDouble *d_ycc, RDouble *d_zcc, int limit_vec, RFloat *d_limit,
                                                 RFloat *d_LIMIT, RFloat *d_ql_ns, RFloat *d_qr_ns,
                                                 RFloat *d_dqdx_proxy, RFloat *d_dqdy_proxy, RFloat *d_dqdz_proxy,
                                                 int nst, int ned, int d_neqn_ns, int d_nTotal, int d_SEG_LEN);
    //! call function for kernel GPUReConstructFaceValueLoop2
    void CallGPUReConstructFaceValueLoop2(Grid *grid_in, const int limit_vec, int nst, int ned, const int SEG_LEN);
    //! GPUReConstructFaceValue kernel for loop2
    __global__ void GPUReConstructFaceValueLoop2(RDouble *d_xfc, RDouble *d_yfc, RDouble *d_zfc,
                                                 int *d_left_cell_of_face, int *d_right_cell_of_face, RDouble *d_xcc,
                                                 RDouble *d_ycc, RDouble *d_zcc, int limit_vec, RFloat *d_limit,
                                                 RFloat *d_LIMIT, RFloat *d_ql_ns, RFloat *d_qr_ns,
                                                 RFloat *d_dqdx_proxy, RFloat *d_dqdy_proxy, RFloat *d_dqdz_proxy,
                                                 int nst, int ned, int d_neqn_ns, int d_nTotal, int d_SEG_LEN);
    //! GPUReConstructFaceValue kernel for second loop
    __global__ void GPUReConstructFaceValueLoop2_S1(RDouble *d_xfc, RDouble *d_yfc, RDouble *d_zfc,
                                                    int *d_left_cell_of_face, int *d_right_cell_of_face, RDouble *d_xcc,
                                                    RDouble *d_ycc, RDouble *d_zcc, int limit_vec, RFloat *d_limit,
                                                    RFloat *d_LIMIT, RFloat *d_ql_ns, RFloat *d_qr_ns,
                                                    RFloat *d_dqdx_proxy, RFloat *d_dqdy_proxy, RFloat *d_dqdz_proxy,
                                                    int nst, int ned, int d_neqn_ns, int d_nTotal, int d_SEG_LEN);
    //! call function for BoundaryQlQrFix loop
    void CallGPUBoundaryQlQrFix(Grid *grid_in, FaceProxy *face_proxy, Param_NSSolverUnstruct *parameters, int nst,
                                int ned, const int SEG_LEN);
    //! kernel function for BoundaryQlQrFix loop
    __global__ void GPUBoundaryQlQrFixOpt(RFloat *d_ql_ns, RFloat *d_qr_ns, int *d_left_cell_of_face,
                                          int *d_right_cell_of_face, RFloat *d_q_ns, RDouble *d_xtn, RDouble *d_ytn,
                                          RDouble *d_ztn, RDouble *d_vgn, int *d_boundaryType, RDouble twall,
                                          int iviscous, double refDimensionalTemperature,
                                          double coefficientOfStateEquation, const int nst, const int nMid,
                                          const int neqn, const int nTotal, const int SEG_LEN, const int SYMMETRY,
                                          const int SOLID_SURFACE, RDouble *xfn, RDouble *yfn, RDouble *zfn);
    __global__ void GPUBoundaryQlQrFix(RFloat *d_ql_ns, RFloat *d_qr_ns, int *d_left_cell_of_face,
                                       int *d_right_cell_of_face, RFloat *d_q_ns, RDouble *d_xtn, RDouble *d_ytn,
                                       RDouble *d_ztn, int *d_boundaryType, RDouble twall, int iviscous,
                                       double refDimensionalTemperature, double coefficientOfStateEquation,
                                       const int nst, const int nMid, const int neqn, const int nTotal,
                                       const int SEG_LEN, const int SYMMETRY, const int SOLID_SURFACE, RDouble *xfn,
                                       RDouble *yfn, RDouble *zfn);
    //! call function for Roe_Scheme_Old
    void CallGPUinviscidScheme(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara, const int nst,
                               const int ned);
    //! kernel function for Roe_Scheme_Old
    __global__ void GPURoe_Scheme_Old(RFloat *d_ql_ns, RFloat *d_qr_ns, RDouble *xfn, RDouble *yfn, RDouble *zfn,
                                      RDouble *vgn, RFloat *gamal, RFloat *gamar, RFloat *flux, const int neqn,
                                      const int nl, const int nlen, const int nm, const int RoeEntropyFixMethod,
                                      const int nchem, const int SEG_LEN, const int nst, const int ned,
                                      const int precon, RFloat alf_l, RFloat alf_n, RFloat mach2);
    __global__ void GPURoe_Scheme_Old_S1(RFloat *d_ql_ns, RFloat *d_qr_ns, RDouble *xfn, RDouble *yfn, RDouble *zfn,
                                         RDouble *vgn, RFloat *gamal, RFloat *gamar, RFloat *flux, const int neqn,
                                         const int nl, const int nlen, const int nm, const int RoeEntropyFixMethod,
                                         const int nchem, const int SEG_LEN, const int nst, const int ned,
                                         const int precon, RFloat alf_l, RFloat alf_n, RFloat mach2);
    //! call kernel for GPUInviscidFluxWrap's first for loop
    void CallGPUInviscidFluxWrapLoop1(const int nl, const int nst, const int ned, const int SEG_LEN);
    //! kernel for GPUInviscidFluxWrap's first for loop
    __global__ void GPUInviscidFluxWrapLoop1(RFloat *d_flux, RDouble *d_area, const int nl, const int nst,
                                             const int ned, const int SEG_LEN);
    //! call kernel for LoadFlux
    void CallGPULoadFlux(Grid *grid_in, FaceProxy *faceProxy, int nst, int ned);
    //! call kernel for GPULoadFluxLoop1
    void CallGPULoadFluxLoop1(const int nl, const int nTotal, const int SEG_LEN, const int nst, const int nMid);
    __global__ void GPULoadFluxLoop1(RFloat *d_res_ns, RFloat *d_flux, int *d_left_cell_of_face, const int nl,
                                     const int nTotal, const int SEG_LEN, const int nst, const int nMid);
    //! call kernel for GPULoadFluxLoop2
    void CallGPULoadFluxLoop2(const int nl, const int nTotal, const int SEG_LEN, const int nst, const int nMid,
                              const int ned);
    __global__ void GPULoadFluxLoop2(RFloat *d_res_ns, RFloat *d_flux, int *d_left_cell_of_face,
                                     int *d_right_cell_of_face, const int nl, const int nTotal, const int SEG_LEN,
                                     const int nst, const int nMid, const int ned);

    //! void CallGPUCompGamaAndTField_S(double refGama,
    //!                                 RFloat coefficientOfStateEquation, const int nl,
    //!                                 const int nchem, const int nTotal);
    void CallGPUComputeGamaAndTemperature(Grid *grid_in, Param_NSSolverUnstruct *parameters);

    __global__ void GPUCompGamaAndTField_S(RFloat *q, RFloat *gama, RFloat *t, double refGama,
                                           RFloat coefficientOfStateEquation, RFloat omav, const int IR, const int IP,
                                           const int ITT, const int nTotal);

    void CallGPUSetFluxSubZero(const int SEG_LEN);

    __global__ void GPUSetFluxSubZero(const int SEG_LEN, RFloat *flux_sub);
} //! namespace GPUNSSolverUnstruct

//! check the rho and pressure
__device__ bool GPUPositiveCheck(RFloat *qtry);
