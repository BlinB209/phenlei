#include <iostream>
#include "Geo_SimpleBC.h"
#include "BasicDeviceVariables.h"
#include "Constants.h"
#include "Gas.h"
#include "GPUBasicFunctions.h"
#include "GPUDeviceControl.h"
#include "GPUInviscidFlux.h"
#include "cudaErrorHandle.h"
#if __CUDA_ARCH__ < 600
#include "GPUPortability.h"
#endif

namespace GPUNSSolverUnstruct
{
    //! set d_flux_sub_turb zero for LoadFlux in InviscidFlux
    void CallGPUSetFluxSubZero(const int SEG_LEN)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        //! GPUSetFluxSubZero<<<160, 1024>>>(SEG_LEN, d_flux_sub_turb);
        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = SEG_LEN;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUSetFluxSubZero, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUSetFluxSubZero, 0, gridSize, blockSize);
#endif

        GPUSetFluxSubZero<<<gridSize, blockSize>>>(SEG_LEN, d_flux_sub_turb);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    //! call function for kernel GetQlQrGPU
    void CallGPUGetQlQr(Grid *grid_in, const int nst, const int ned, const int SEG_LEN)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid       = UnstructGridCast(grid_in);
        const int     nTotalCell = grid->GetNTotalCell();
        const int     nBoundFace = grid->GetNBoundFace();
        int           nTotal     = nTotalCell + nBoundFace;
        int           nl;
        GlobalDataBase::GetData("nl", &nl, PHINT, 1);
        int nchem;
        GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);
        int neqn = nl + nchem;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = ned - nst;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUGetQlQr, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUGetQlQr, 0, gridSize, blockSize);
#endif

        //! GetQlQr kernel
        GPUGetQlQr<<<gridSize, blockSize>>>(d_left_cell_of_face, d_right_cell_of_face, d_q_ns, d_ql_ns, d_qr_ns, nst,
                                            ned, neqn, nTotal, SEG_LEN);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    //! GetQlQr GPU kernel
    __global__ void GPUGetQlQr(int *d_left_cell_of_face, int *d_right_cell_of_face, RFloat *d_q_ns, RFloat *d_ql_ns,
                               RFloat *d_qr_ns, const int nst, const int ned, const int d_neqn, const int d_nTotal,
                               const int d_SEG_LEN)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int iFace = nst + bidx * blockDim.x + tidx; iFace < ned; iFace += blockDim.x * gridDim.x)
        {
            int le, re, j;
            le = d_left_cell_of_face[iFace];
            re = d_right_cell_of_face[iFace];
            j  = iFace - nst;

            for (int m = 0; m < d_neqn; ++m)
            {
                d_ql_ns[m * d_SEG_LEN + j] = d_q_ns[m * d_nTotal + le];
                d_qr_ns[m * d_SEG_LEN + j] = d_q_ns[m * d_nTotal + re];
            }
        }
    }
    //! call function for kernel GPUGetGamaLR
    void CallGPUGetGamaLR(const int nst, const int ned)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = ned - nst;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUGetGamaLR, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUGetGamaLR, 0, gridSize, blockSize);
#endif

        //! GPUGetGamaLR kernel
        GPUGetGamaLR<<<gridSize, blockSize>>>(d_left_cell_of_face, d_right_cell_of_face, d_gama_ns, d_gamaL_ns,
                                              d_gamaR_ns, nst, ned);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    //! GetGamaLR GPU kernel
    __global__ void GPUGetGamaLR(int *d_left_cell_of_face, int *d_right_cell_of_face, RFloat *d_gama_ns,
                                 RFloat *d_gamaL_ns, RFloat *d_gamaR_ns, const int nst, const int ned)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int iFace = nst + bidx * blockDim.x + tidx; iFace < ned; iFace += blockDim.x * gridDim.x)
        {
            int le, re, j;
            le = d_left_cell_of_face[iFace];
            re = d_right_cell_of_face[iFace];
            j  = iFace - nst;

            d_gamaL_ns[j] = d_gama_ns[le];
            d_gamaR_ns[j] = d_gama_ns[re];
        }
    }

    void CallGPUReConstructFaceValue(Grid *gridIn, Limiter *limiter, int localStart, int localEnd, const int SEG_LEN,
                                     Param_NSSolverUnstruct *parameters)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid = UnstructGridCast(gridIn);

        if (!grid->IsFinestGrid())
        {
            //! Reconstruction is not necessary for coarse grid.
            return;
        }

        //!Param_NSSolverUnstruct *parameters = GetControlParameters();
        int  limiterType  = parameters->GetLimiterType();
        int  isInIniting  = GlobalDataBase::GetIntParaFromDB("isInIniting");
        bool isFirstOrder = (limiterType == ILMT_FIRST) || isInIniting;
        if (isFirstOrder)
        {
            return;
        }

        //!int nEquation = GetNumberOfEquations();

        bool isStructuredLimiter = (limiterType >= ILMT_STRUCT);
        if (isStructuredLimiter)
        {
            //! Limiter for structured limiter.//!gpu not use
            //! ReConstructQlQr_STR(gridIn, faceProxy, nEquation, localStart, localEnd);
            return;
        }

        int       reconstructMethod                    = GlobalDataBase::GetIntParaFromDB("reconmeth");
        const int RECONSTRUCT_USING_RESPECTIVE_LIMITER = 0;

        int limiterVector = parameters->GetLimitVector();
        if (reconstructMethod == RECONSTRUCT_USING_RESPECTIVE_LIMITER)
        {
            CallGPUReConstructFaceValueLoop1(gridIn, limiterVector, localStart, localEnd, SEGCTION_LENGTH);
        }
        else
        {
            CallGPUReConstructFaceValueLoop2(gridIn, limiterVector, localStart, localEnd, SEGCTION_LENGTH);
        }
    }
    //! call function for kernel GPUReConstructFaceValueLoop1
    void CallGPUReConstructFaceValueLoop1(Grid *grid_in, const int limit_vec, int nst, int ned, const int SEG_LEN)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid       = UnstructGridCast(grid_in);
        int           nBoundFace = grid->GetNBoundFace();
        int           nTotalCell = grid->GetNTotalCell();
        int           nTotal     = nBoundFace + nTotalCell;

        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        int neqn  = nl + nchem;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = ned - nst;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUReConstructFaceValueLoop1, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUReConstructFaceValueLoop1, 0, gridSize, blockSize);
#endif
        GPUReConstructFaceValueLoop1<<<gridSize, blockSize>>>(
            d_xfc, d_yfc, d_zfc, d_left_cell_of_face, d_right_cell_of_face, d_xcc, d_ycc, d_zcc, limit_vec, d_limit,
            d_LIMIT, d_ql_ns, d_qr_ns, d_dqdx_proxy, d_dqdy_proxy, d_dqdz_proxy, nst, ned, neqn, nTotal, SEG_LEN);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    //! GPUReConstructFaceValue kernel for loop1
    __global__ void GPUReConstructFaceValueLoop1(RDouble *d_xfc, RDouble *d_yfc, RDouble *d_zfc,
                                                 int *d_left_cell_of_face, int *d_right_cell_of_face, RDouble *d_xcc,
                                                 RDouble *d_ycc, RDouble *d_zcc, int limit_vec, RFloat *d_limit,
                                                 RFloat *d_LIMIT, RFloat *d_ql_ns, RFloat *d_qr_ns,
                                                 RFloat *d_dqdx_proxy, RFloat *d_dqdy_proxy, RFloat *d_dqdz_proxy,
                                                 int nst, int ned, int d_neqn_ns, int d_nTotal, int d_SEG_LEN)
    {
        RDouble dx, dy, dz;
        RFloat *qtry = new RFloat[d_neqn_ns];
        RFloat *d_limit_tmp;
        d_limit_tmp = d_limit;

        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int i = nst + bidx * blockDim.x + tidx; i < ned; i += blockDim.x * gridDim.x)
        {
            int le, re, j;
            j  = i - nst;
            le = d_left_cell_of_face[i];
            re = d_right_cell_of_face[i];

            dx = d_xfc[i] - d_xcc[le];
            dy = d_yfc[i] - d_ycc[le];
            dz = d_zfc[i] - d_zcc[le];

            for (int m = 0; m < d_neqn_ns; ++m)
            {
                qtry[m] = d_ql_ns[m * d_SEG_LEN + j];
            }

            for (int m = 0; m < d_neqn_ns; ++m)
            {
                if (limit_vec == 1)
                {
                    d_limit_tmp = &(d_LIMIT[m * d_nTotal]);
                }
                //! qtry[m] += limit[le] * ( dqdx[m][le] * dx + dqdy[m][le] * dy +
                //! dqdz[m][le] * dz );
                qtry[m] += d_limit_tmp[le]
                           * (d_dqdx_proxy[m * d_nTotal + le] * dx + d_dqdy_proxy[m * d_nTotal + le] * dy
                              + d_dqdz_proxy[m * d_nTotal + le] * dz);
            }

            if (GPUPositiveCheck(qtry))
            {
                for (int m = 0; m < d_neqn_ns; ++m)
                {
                    d_ql_ns[m * d_SEG_LEN + j] = qtry[m];
                }
            }

            dx = d_xfc[i] - d_xcc[re];
            dy = d_yfc[i] - d_ycc[re];
            dz = d_zfc[i] - d_zcc[re];

            for (int m = 0; m < d_neqn_ns; ++m)
            {
                qtry[m] = d_qr_ns[m * d_SEG_LEN + j];
            }

            for (int m = 0; m < d_neqn_ns; ++m)
            {
                if (limit_vec == 1)
                {
                    d_limit_tmp = &(d_LIMIT[m * d_nTotal]);
                }
                //! qtry[m] += limit[re] * ( dqdx[m][re] * dx + dqdy[m][re] * dy +
                //! dqdz[m][re] * dz );
                qtry[m] += d_limit_tmp[re]
                           * (d_dqdx_proxy[m * d_nTotal + re] * dx + d_dqdy_proxy[m * d_nTotal + re] * dy
                              + d_dqdz_proxy[m * d_nTotal + re] * dz);
            }

            if (GPUPositiveCheck(qtry))
            {
                for (int m = 0; m < d_neqn_ns; ++m)
                {
                    d_qr_ns[m * d_SEG_LEN + j] = qtry[m];
                }
            }
        }

        delete[] qtry;
    }
    //! call function for kernel GPUReConstructFaceValueLoop2
    void CallGPUReConstructFaceValueLoop2(Grid *grid_in, const int limit_vec, int nst, int ned, const int SEG_LEN)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid = UnstructGridCast(grid_in);

        int nBoundFace = grid->GetNBoundFace();
        int nTotalCell = grid->GetNTotalCell();
        int nTotal     = nBoundFace + nTotalCell;

        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        int neqn  = nl + nchem;

        int gridSize        = 1;
        int blockSize       = 1;
        int loopLen         = ned - nst;
        int regsPerThread   = 80;
        int residentBlockNo = 4;

        size_t dsMemPerThread = neqn * sizeof(RFloat); //! calculate by the intended shared memory use in the kernel
        size_t dsMemPerBlock;
        KernelLaunchPara((void *)GPUReConstructFaceValueLoop2_S1, loopLen, regsPerThread, residentBlockNo,
                         dsMemPerThread, dsMemPerBlock, gridSize, blockSize, GPUProp);

//! KernelLaunchPara((void*)GPUReConstructFaceValueLoop2, loopLen, 0, gridSize,
//! blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUReConstructFaceValueLoop2_S1, dsMemPerBlock, gridSize, blockSize);
#endif
        GPUReConstructFaceValueLoop2_S1<<<gridSize, blockSize, dsMemPerBlock>>>(
            d_xfc, d_yfc, d_zfc, d_left_cell_of_face, d_right_cell_of_face, d_xcc, d_ycc, d_zcc, limit_vec, d_limit,
            d_LIMIT, d_ql_ns, d_qr_ns, d_dqdx_proxy, d_dqdy_proxy, d_dqdz_proxy, nst, ned, neqn, nTotal, SEG_LEN);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    //! GPUReConstructFaceValue kernel for second loop
    __global__ void GPUReConstructFaceValueLoop2(RDouble *d_xfc, RDouble *d_yfc, RDouble *d_zfc,
                                                 int *d_left_cell_of_face, int *d_right_cell_of_face, RDouble *d_xcc,
                                                 RDouble *d_ycc, RDouble *d_zcc, int limit_vec, RFloat *d_limit,
                                                 RFloat *d_LIMIT, RFloat *d_ql_ns, RFloat *d_qr_ns,
                                                 RFloat *d_dqdx_proxy, RFloat *d_dqdy_proxy, RFloat *d_dqdz_proxy,
                                                 int nst, int ned, int d_neqn_ns, int d_nTotal, int d_SEG_LEN)
    {
        RDouble dx, dy, dz;
        RFloat *qtry = new RFloat[d_neqn_ns];
        RFloat *d_limit_tmp;
        d_limit_tmp = d_limit;

        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int i = nst + bidx * blockDim.x + tidx; i < ned; i += blockDim.x * gridDim.x)
        {
            int le, re, j;
            j  = i - nst;
            le = d_left_cell_of_face[i];
            re = d_right_cell_of_face[i];

            dx = d_xfc[i] - d_xcc[le];
            dy = d_yfc[i] - d_ycc[le];
            dz = d_zfc[i] - d_zcc[le];

            for (int m = 0; m < d_neqn_ns; ++m)
            {
                qtry[m] = d_ql_ns[m * d_SEG_LEN + j];
            }

            for (int m = 0; m < d_neqn_ns; ++m)
            {
                if (limit_vec == 1)
                {
                    d_limit_tmp = &(d_LIMIT[m * d_nTotal]);
                }

                qtry[m] += GPUMIN(d_limit_tmp[le], d_limit_tmp[re])
                           * (d_dqdx_proxy[m * d_nTotal + le] * dx + d_dqdy_proxy[m * d_nTotal + le] * dy
                              + d_dqdz_proxy[m * d_nTotal + le] * dz);
            }

            if (GPUPositiveCheck(qtry))
            {
                for (int m = 0; m < d_neqn_ns; ++m)
                {
                    d_ql_ns[m * d_SEG_LEN + j] = qtry[m];
                }
            }

            dx = d_xfc[i] - d_xcc[re];
            dy = d_yfc[i] - d_ycc[re];
            dz = d_zfc[i] - d_zcc[re];

            for (int m = 0; m < d_neqn_ns; ++m)
            {
                qtry[m] = d_qr_ns[m * d_SEG_LEN + j];
            }

            for (int m = 0; m < d_neqn_ns; ++m)
            {
                if (limit_vec == 1)
                {
                    d_limit_tmp = &(d_LIMIT[m * d_nTotal]);
                }

                qtry[m] += GPUMIN(d_limit_tmp[le], d_limit_tmp[re])
                           * (d_dqdx_proxy[m * d_nTotal + re] * dx + d_dqdy_proxy[m * d_nTotal + re] * dy
                              + d_dqdz_proxy[m * d_nTotal + re] * dz);
            }

            if (GPUPositiveCheck(qtry))
            {
                for (int m = 0; m < d_neqn_ns; ++m)
                {
                    d_qr_ns[m * d_SEG_LEN + j] = qtry[m];
                }
            }
        }

        delete[] qtry;
    }
    //! GPUReConstructFaceValue kernel for second loop
    __global__ void GPUReConstructFaceValueLoop2_S1(RDouble *d_xfc, RDouble *d_yfc, RDouble *d_zfc,
                                                    int *d_left_cell_of_face, int *d_right_cell_of_face, RDouble *d_xcc,
                                                    RDouble *d_ycc, RDouble *d_zcc, int limit_vec, RFloat *d_limit,
                                                    RFloat *d_LIMIT, RFloat *d_ql_ns, RFloat *d_qr_ns,
                                                    RFloat *d_dqdx_proxy, RFloat *d_dqdy_proxy, RFloat *d_dqdz_proxy,
                                                    int nst, int ned, int d_neqn_ns, int d_nTotal, int d_SEG_LEN)
    {
        extern __shared__ RFloat array[];
        RFloat                  *qtry = array;

        RDouble dx, dy, dz;
        RFloat *d_limit_tmp;
        d_limit_tmp = d_limit;

        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int i = nst + bidx * blockDim.x + tidx; i < ned; i += blockDim.x * gridDim.x)
        {
            int le, re, j;
            j  = i - nst;
            le = d_left_cell_of_face[i];
            re = d_right_cell_of_face[i];

            dx = d_xfc[i] - d_xcc[le];
            dy = d_yfc[i] - d_ycc[le];
            dz = d_zfc[i] - d_zcc[le];

            for (int m = 0; m < d_neqn_ns; ++m)
            {
                qtry[m * blockDim.x + tidx] = d_ql_ns[m * d_SEG_LEN + j];
            }

            for (int m = 0; m < d_neqn_ns; ++m)
            {
                if (limit_vec == 1)
                {
                    d_limit_tmp = &(d_LIMIT[m * d_nTotal]);
                }

                qtry[m * blockDim.x + tidx] +=
                    GPUMIN(d_limit_tmp[le], d_limit_tmp[re])
                    * (d_dqdx_proxy[m * d_nTotal + le] * dx + d_dqdy_proxy[m * d_nTotal + le] * dy
                       + d_dqdz_proxy[m * d_nTotal + le] * dz);
            }

            //! if ( GPUPositiveCheck(qtry) )
            if (qtry[0 * blockDim.x + tidx] > 0.0 && qtry[4 * blockDim.x + tidx] > 0.0)
            {
                for (int m = 0; m < d_neqn_ns; ++m)
                {
                    d_ql_ns[m * d_SEG_LEN + j] = qtry[m * blockDim.x + tidx];
                }
            }

            dx = d_xfc[i] - d_xcc[re];
            dy = d_yfc[i] - d_ycc[re];
            dz = d_zfc[i] - d_zcc[re];

            for (int m = 0; m < d_neqn_ns; ++m)
            {
                qtry[m * blockDim.x + tidx] = d_qr_ns[m * d_SEG_LEN + j];
            }

            for (int m = 0; m < d_neqn_ns; ++m)
            {
                if (limit_vec == 1)
                {
                    d_limit_tmp = &(d_LIMIT[m * d_nTotal]);
                }

                qtry[m * blockDim.x + tidx] +=
                    GPUMIN(d_limit_tmp[le], d_limit_tmp[re])
                    * (d_dqdx_proxy[m * d_nTotal + re] * dx + d_dqdy_proxy[m * d_nTotal + re] * dy
                       + d_dqdz_proxy[m * d_nTotal + re] * dz);
            }

            //! if ( GPUPositiveCheck(qtry) )
            if (qtry[0 * blockDim.x + tidx] > 0.0 && qtry[4 * blockDim.x + tidx] > 0.0)
            {
                for (int m = 0; m < d_neqn_ns; ++m)
                {
                    d_qr_ns[m * d_SEG_LEN + j] = qtry[m * blockDim.x + tidx];
                }
            }
        }

        //! delete []qtry;
    }
    //! call function for BoundaryQlQrFix loop
    void CallGPUBoundaryQlQrFix(Grid *grid_in, FaceProxy *face_proxy, Param_NSSolverUnstruct *parameters, int nst,
                                int ned, const int SEG_LEN)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid = UnstructGridCast(grid_in);

        const int nTotalCell = grid->GetNTotalCell();
        const int nBoundFace = grid->GetNBoundFace();
        int       nTotal     = nTotalCell + nBoundFace;
        int       nMid;

        //! Check if there are boundary faces. If no, return.
        if (nst >= nBoundFace) return;

        if (ned <= nBoundFace)
        {
            //! If they are all boundary faces
            nMid = ned;
        }
        else
        {
            //! Part of them are boundary faces
            nMid = nBoundFace;
        }

        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        int neqn  = nl + nchem;

        int     iviscous                   = parameters->GetViscousType();
        RDouble twall                      = parameters->GetWallTemperature();
        double  refDimensionalTemperature  = parameters->GetRefDimensionalTemperature();
        double  refGama                    = parameters->GetRefGama();
        double  coefficientOfStateEquation = GlobalDataBase::GetDoubleParaFromDB("coefficientOfStateEquation");

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nMid - nst;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUBoundaryQlQrFixOpt, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUBoundaryQlQrFixOpt, 0, gridSize, blockSize);
#endif

        GPUBoundaryQlQrFixOpt<<<gridSize, blockSize>>>(
            d_ql_ns, d_qr_ns, d_left_cell_of_face, d_right_cell_of_face, d_q_ns, d_xtn, d_ytn, d_ztn, d_vgn,
            d_boundaryType, twall, iviscous, refDimensionalTemperature, coefficientOfStateEquation, nst, nMid, neqn,
            nTotal, SEG_LEN, PHENGLEI::SYMMETRY, PHENGLEI::SOLID_SURFACE, d_xfn, d_yfn, d_zfn);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    //! kernel function for BoundaryQlQrFix loop
    __global__ void GPUBoundaryQlQrFixOpt(RFloat *d_ql_ns, RFloat *d_qr_ns, int *d_left_cell_of_face,
                                          int *d_right_cell_of_face, RFloat *d_q_ns, RDouble *d_xtn, RDouble *d_ytn,
                                          RDouble *d_ztn, RDouble *d_vgn, int *d_boundaryType, RDouble twall,
                                          int iviscous, double refDimensionalTemperature,
                                          double coefficientOfStateEquation, const int nst, const int nMid,
                                          const int neqn, const int nTotal, const int SEG_LEN, const int SYMMETRY,
                                          const int SOLID_SURFACE, RDouble *xfn, RDouble *yfn, RDouble *zfn)
    {
        RFloat prim[5];

        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int i = nst + bidx * blockDim.x + tidx; i < nMid; i += blockDim.x * gridDim.x)
        {
            int j = i - nst;
            //! PHENGLEI::INTERFACE
            if (d_boundaryType[i] != SOLID_SURFACE && d_boundaryType[i] != SYMMETRY) continue;

            int le = d_left_cell_of_face[i];
            int re = d_right_cell_of_face[i];

            for (int m = 0; m < neqn; ++m)
            {
                RFloat temp = 0.5 * (d_q_ns[m * nTotal + le] + d_q_ns[m * nTotal + re]);

                d_ql_ns[m * SEG_LEN + j] = temp;
                d_qr_ns[m * SEG_LEN + j] = temp;
            }

            //! SYMMETRY
            if (d_boundaryType[i] == SYMMETRY)
            {
                RDouble vn = 2.0
                             * (xfn[i] * d_ql_ns[1 * SEG_LEN + j] + yfn[i] * d_ql_ns[2 * SEG_LEN + j]
                                + zfn[i] * d_ql_ns[3 * SEG_LEN + j] - d_vgn[i]);

                d_qr_ns[0 * SEG_LEN + j] = d_ql_ns[0 * SEG_LEN + j];
                d_qr_ns[1 * SEG_LEN + j] = d_ql_ns[1 * SEG_LEN + j] - vn * xfn[i];
                d_qr_ns[2 * SEG_LEN + j] = d_ql_ns[2 * SEG_LEN + j] - vn * yfn[i];
                d_qr_ns[3 * SEG_LEN + j] = d_ql_ns[3 * SEG_LEN + j] - vn * zfn[i];
                d_qr_ns[4 * SEG_LEN + j] = d_ql_ns[4 * SEG_LEN + j];
            }
            //! PHENGLEI::SOLID_SURFACE
            if (d_boundaryType[i] == SOLID_SURFACE && iviscous)
            {
                //! IU 1, IV 2, IW 3
                            /*
                d_ql_ns[1*SEG_LEN+j] = d_xtn[i];
                d_ql_ns[2*SEG_LEN+j] = d_ytn[i];
                d_ql_ns[3*SEG_LEN+j] = d_ztn[i];

                d_qr_ns[1*SEG_LEN+j] = d_xtn[i];
                d_qr_ns[2*SEG_LEN+j] = d_ytn[i];
                d_qr_ns[3*SEG_LEN+j] = d_ztn[i];
                    */
                //!d_qr_ns[1 * SEG_LEN + j] = -d_ql_ns[1 * SEG_LEN + j] + 2.0 * d_xtn[i];
                //!d_qr_ns[2 * SEG_LEN + j] = -d_ql_ns[2 * SEG_LEN + j] + 2.0 * d_ytn[i];
                //!d_qr_ns[3 * SEG_LEN + j] = -d_ql_ns[3 * SEG_LEN + j] + 2.0 * d_ztn[i];

                d_ql_ns[1 * SEG_LEN + j] = 0.0;
                d_ql_ns[2 * SEG_LEN + j] = 0.0;
                d_ql_ns[3 * SEG_LEN + j] = 0.0;

                d_qr_ns[1 * SEG_LEN + j] = 0.0;
                d_qr_ns[2 * SEG_LEN + j] = 0.0;
                d_qr_ns[3 * SEG_LEN + j] = 0.0;

                //! d_qr_ns[1 * SEG_LEN + j] = 0.0;
                //! d_qr_ns[2 * SEG_LEN + j] = 0.0;
                //! d_qr_ns[3 * SEG_LEN + j] = 0.0;

                //! d_ql_ns[1 * SEG_LEN + j] = 0.0;
                //! d_ql_ns[2 * SEG_LEN + j] = 0.0;
                //! d_ql_ns[3 * SEG_LEN + j] = 0.0;

                if (twall <= 0.0)
                {
                    //! Viscous WALL, adiabatic
                }
                else
                {
                    //! iso-thermal wall.
                    RFloat tw = twall / refDimensionalTemperature;

                    for (int m = 0; m < neqn; ++m)
                    {
                        prim[m] = d_ql_ns[m * SEG_LEN + j];
                    }
                    //! RFloat omav = one;
                    RFloat omav = 1.0;
                    //! this chemical part will be neglected firstly
                    //! using namespace CHEMICAL_SPACE;
                    /*
                    if ( nchem == 1 )
                    {
                                                            if(i == nst)
                                                            {
                                                                    printf("chemisty part is
                    not supported in cuda mode\n"); exit(1);
                                                            }
                            //!ComputeMolecularWeightReciprocal(prim, omav);
                    }
                    */
                    //! RFloat presw = prim[IP];
                    RFloat presw = prim[4];
                    RFloat rhow  = presw / (coefficientOfStateEquation * tw * omav);
                    //! ql[IR][j] = rhow;
                    //! qr[IR][j] = rhow;
                    d_ql_ns[0 * SEG_LEN + j] = rhow;
                    d_qr_ns[0 * SEG_LEN + j] = rhow;
                }
            }
        }
    }

    __global__ void GPUBoundaryQlQrFix(RFloat *d_ql_ns, RFloat *d_qr_ns, int *d_left_cell_of_face,
                                       int *d_right_cell_of_face, RFloat *d_q_ns, RDouble *d_xtn, RDouble *d_ytn,
                                       RDouble *d_ztn, int *d_boundaryType, RDouble twall, int iviscous,
                                       double refDimensionalTemperature, double coefficientOfStateEquation,
                                       const int nst, const int nMid, const int neqn, const int nTotal,
                                       const int SEG_LEN, const int SYMMETRY, const int SOLID_SURFACE, RDouble *xfn,
                                       RDouble *yfn, RDouble *zfn)
    {
        RFloat *prim = new RFloat[neqn];

        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int i = nst + bidx * blockDim.x + tidx; i < nMid; i += blockDim.x * gridDim.x)
        {
            int j = i - nst;
            //! PHENGLEI::INTERFACE
            if (d_boundaryType[i] != SOLID_SURFACE && d_boundaryType[i] != SYMMETRY) continue;

            int le = d_left_cell_of_face[i];
            int re = d_right_cell_of_face[i];
            /*
            for ( int m = 0; m < neqn; ++ m )
            {
            RFloat temp =  0.5 * ( d_q_ns[m*nTotal+le] + d_q_ns[m*nTotal+re] );

            d_ql_ns[m*SEG_LEN+j] = temp;
            d_qr_ns[m*SEG_LEN+j] = temp;
            }
                */
            //! SYMMETRY
            if (d_boundaryType[i] == SYMMETRY)
            {
                RDouble vn = 2.0
                             * (xfn[i] * d_ql_ns[1 * SEG_LEN + j] + yfn[i] * d_ql_ns[2 * SEG_LEN + j]
                                + zfn[i] * d_ql_ns[3 * SEG_LEN + j]);

                d_qr_ns[0 * SEG_LEN + j] = d_ql_ns[0 * SEG_LEN + j];
                d_qr_ns[1 * SEG_LEN + j] = d_ql_ns[1 * SEG_LEN + j] - vn * xfn[i];
                d_qr_ns[2 * SEG_LEN + j] = d_ql_ns[2 * SEG_LEN + j] - vn * yfn[i];
                d_qr_ns[3 * SEG_LEN + j] = d_ql_ns[3 * SEG_LEN + j] - vn * zfn[i];
                d_qr_ns[4 * SEG_LEN + j] = d_ql_ns[4 * SEG_LEN + j];
            }
            //! PHENGLEI::SOLID_SURFACE
            if (d_boundaryType[i] == 2 && iviscous)
            {
                //! IU 1, IV 2, IW 3
                /*
                d_ql_ns[1*SEG_LEN+j] = d_xtn[i];
                d_ql_ns[2*SEG_LEN+j] = d_ytn[i];
                d_ql_ns[3*SEG_LEN+j] = d_ztn[i];

                d_qr_ns[1*SEG_LEN+j] = d_xtn[i];
                d_qr_ns[2*SEG_LEN+j] = d_ytn[i];
                d_qr_ns[3*SEG_LEN+j] = d_ztn[i];
                    */
                d_qr_ns[1 * SEG_LEN + j] = -d_ql_ns[1 * SEG_LEN + j] + 2.0 * d_xtn[i];
                d_qr_ns[2 * SEG_LEN + j] = -d_ql_ns[2 * SEG_LEN + j] + 2.0 * d_ytn[i];
                d_qr_ns[3 * SEG_LEN + j] = -d_ql_ns[3 * SEG_LEN + j] + 2.0 * d_ztn[i];

                if (twall <= 0.0)
                {
                    //! Viscous WALL, adiabatic
                }
                else
                {
                    //! iso-thermal wall.
                    RFloat tw = twall / refDimensionalTemperature;

                    for (int m = 0; m < neqn; ++m)
                    {
                        prim[m] = d_ql_ns[m * SEG_LEN + j];
                    }
                    //! RFloat omav = one;
                    RFloat omav = 1.0;
                    //! this chemical part will be neglected firstly
                    //! using namespace CHEMICAL_SPACE;
                    /*
                    if ( nchem == 1 )
                    {
                                                            if(i == nst)
                                                            {
                                                                    printf("chemisty part is
                    not supported in cuda mode\n"); exit(1);
                                                            }
                            //!ComputeMolecularWeightReciprocal(prim, omav);
                    }
                    */
                    //! RFloat presw = prim[IP];
                    RFloat presw = prim[4];
                    RFloat rhow  = presw / (coefficientOfStateEquation * tw * omav);
                    //! ql[IR][j] = rhow;
                    //! qr[IR][j] = rhow;
                    d_ql_ns[0 * SEG_LEN + j] = rhow;
                    d_qr_ns[0 * SEG_LEN + j] = rhow;
                }
            }
        }
        delete[] prim;
    }
    //! call function for Roe_Scheme_Old
    void CallGPUinviscidScheme(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara, const int nst,
                               const int ned)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int neqn                = invSchemePara->GetNumberOfTotalEquation();
        int nl                  = invSchemePara->GetNumberOfLaminar();
        int nlen                = invSchemePara->GetLength();
        int nm                  = invSchemePara->Getnm();
        int RoeEntropyFixMethod = invSchemePara->GetEntropyFixMethod();
        //!int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        //!int precon = GlobalDataBase::GetIntParaFromDB("precon");
        int    nchem         = invSchemePara->GetNumberOfChemical();
        int    precon        = invSchemePara->GetIfPrecondition();
        RFloat refMachNumber = invSchemePara->GetMachNumber();
        RFloat alf_l         = invSchemePara->GetLeftEntropyFixCoefficient();
        RFloat alf_n         = invSchemePara->GetRightEntropyFixCoefficient();
        RFloat mach2         = refMachNumber * refMachNumber;

        if (precon == 1 && refMachNumber < 0.6)
        {
            alf_l = 1.0e-8;
            alf_n = 1.0e-8;
        }
        /*
        int gridSize= 1;
        int blockSize =1;
        int loopLen= ned-nst;
        int blockSizeLimit= 128;
        */
        int    gridSize        = 1;
        int    blockSize       = 1;
        int    loopLen         = ned - nst;
        int    regsPerThread   = 224;
        int    residentBlockNo = 2;
        size_t dsMemPerThread = 8 * neqn * sizeof(RFloat); //! calculate by the intended shared memory use in the kernel
        size_t dsMemPerBlock;
        //! int blockSizeLimit;  //! calculated by the register use

        int uns_scheme = GlobalDataBase::GetIntParaFromDB("uns_scheme");
        if (uns_scheme == ISCHEME_ROE)
        {
            KernelLaunchPara((void *)GPURoe_Scheme_Old_S1, loopLen, regsPerThread, residentBlockNo, dsMemPerThread,
                             dsMemPerBlock, gridSize, blockSize, GPUProp);
            //! KernelLaunchPara((void*)GPURoe_Scheme_Old, loopLen, 0, gridSize, blockSize,
            //! blockSizeLimit);
            //!  adjust the parameters
            //! gridSize= 4;
#ifdef KERNELLAUNCHTEST
            //! ReportKernelPara((void*)GPURoe_Scheme_Old, 0, gridSize, blockSize);
            ReportKernelPara((void *)GPURoe_Scheme_Old_S1, dsMemPerBlock, gridSize, blockSize);
#endif
            GPURoe_Scheme_Old_S1<<<gridSize, blockSize, dsMemPerBlock>>>(
                d_ql_ns, d_qr_ns, d_xfn, d_yfn, d_zfn, d_vgn, d_gamaL_ns, d_gamaR_ns, d_flux, neqn, nl, nlen, nm,
                RoeEntropyFixMethod, nchem, d_SEG_LEN, nst, ned, precon, alf_l, alf_n, mach2);
#ifdef KERNELLAUNCHTEST
            HANDLE_KERNEL_ERR();
#endif
        }
        else
        {
            printf("the function for uns_scheme != ISCHEME_ROE is not optimized for "
                   "gpu\n");
            exit(0);
        }
    }
    //! kernel function for Roe_Scheme_Old
    __global__ void GPURoe_Scheme_Old(RFloat *d_ql_ns, RFloat *d_qr_ns, RDouble *xfn, RDouble *yfn, RDouble *zfn,
                                      RDouble *vgn, RFloat *gamal, RFloat *gamar, RFloat *flux, const int neqn,
                                      const int nl, const int nlen, const int nm, const int RoeEntropyFixMethod,
                                      const int nchem, const int SEG_LEN, const int nst, const int ned,
                                      const int precon, RFloat alf_l, RFloat alf_n, RFloat mach2)
    {
        RFloat gama, gamm1;
        RFloat vSquareL, vnL, hint_L, hL, rvnL, EL, gmL;
        RFloat vSquareR, vnR, hint_R, hR, rvnR, ER, gmR;
        RFloat rRoeAverage, uRoeAverage, vRoeAverage, wRoeAverage, pRoeAverage, vSquareRoeAverage, hRoeAverage, vn, c2,
            cm;
        RFloat Ur2, cnew2, cnew, vnew, beta, alpha;
        RFloat ratio, coef;
        RFloat tmp0, tmp1;
        RFloat eigv1, eigv2, eigv3, eigv_max;           //! cL, cR;
        RFloat theta, xi1, xi2, xi3, xi4, dh, dc, c2dc; //! ae,
        RFloat xsn, ysn, zsn, vb;

        RFloat *prim       = new RFloat[neqn];
        RFloat *qConserveL = new RFloat[nl];
        RFloat *qConserveR = new RFloat[nl];
        RFloat *dq         = new RFloat[nl];

        RFloat *primL = new RFloat[neqn];
        RFloat *primR = new RFloat[neqn];
        RFloat *fluxL = new RFloat[nl];
        RFloat *fluxR = new RFloat[nl];

        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int iFace = bidx * blockDim.x + tidx; iFace < (ned - nst); iFace += blockDim.x * gridDim.x)
        {
            for (int m = 0; m < neqn; ++m)
            {
                primL[m] = d_ql_ns[m * SEG_LEN + iFace];
                primR[m] = d_qr_ns[m * SEG_LEN + iFace];
            }

            RFloat &rL = primL[0]; //! ql[0][iFace];
            RFloat &uL = primL[1]; //! ql[1][iFace];
            RFloat &vL = primL[2]; //! ql[2][iFace];
            RFloat &wL = primL[3]; //! ql[3][iFace];
            RFloat &pL = primL[4]; //! ql[IP][iFace];

            RFloat &rR = primR[0]; //! qr[0][iFace];
            RFloat &uR = primR[1]; //! qr[1][iFace];
            RFloat &vR = primR[2]; //! qr[2][iFace];
            RFloat &wR = primR[3]; //! qr[3][iFace];
            RFloat &pR = primR[4]; //! qr[IP][iFace];

            xsn = xfn[nst + iFace];
            ysn = yfn[nst + iFace];
            zsn = zfn[nst + iFace];
            vb  = vgn[nst + iFace];
            gmL = gamal[iFace];
            gmR = gamar[iFace];

            //!###################################################################//!
            //! to get the flux ÓÉ×óÓÒ½çÃæµÄÖµÇóÍ¨Á¿

            vSquareL = uL * uL + vL * vL + wL * wL;
            vSquareR = uR * uR + vR * vR + wR * wR;

            hint_L = (gmL / (gmL - 1.0)) * (pL / rL);
            hint_R = (gmR / (gmR - 1.0)) * (pR / rR);
            //! chemical part is not supported in the kernel by sunxu 20191105
            /*
    if ( nchem != 0 )
    {
            ComputeEnthalpyByPrimitive(primL, gmL, hint_L);
            ComputeEnthalpyByPrimitive(primR, gmR, hint_R);
    }
    */
            hL = hint_L + 0.5 * vSquareL;
            hR = hint_R + 0.5 * vSquareR;

            vnL  = xsn * uL + ysn * vL + zsn * wL - vb;
            vnR  = xsn * uR + ysn * vR + zsn * wR - vb;
            rvnL = rL * vnL;
            rvnR = rR * vnR;

            flux[0 * SEG_LEN + iFace] = 0.5 * (rvnL + rvnR);
            flux[1 * SEG_LEN + iFace] = 0.5 * (rvnL * uL + xsn * pL + rvnR * uR + xsn * pR);
            flux[2 * SEG_LEN + iFace] = 0.5 * (rvnL * vL + ysn * pL + rvnR * vR + ysn * pR);
            flux[3 * SEG_LEN + iFace] = 0.5 * (rvnL * wL + zsn * pL + rvnR * wR + zsn * pR);
            flux[4 * SEG_LEN + iFace] = 0.5 * (rvnL * hL + vb * pL + rvnR * hR + vb * pR);

            for (int m = nm; m < nl; ++m)
            {
                //!            fluxL[m] = 0.5 * ( primL[m] * fluxL[0 ] );
                //!            fluxR[m] = 0.5 * ( primR[m] * fluxR[0 ] );
                flux[m * SEG_LEN + iFace] = 0.5 * (primL[m] * rvnL + primR[m] * rvnR);
            }

            /*  Ö±½Ó¼Óµ½Í¨Á¿
            fluxL[0 ] = 0.5 *   rvnL                         ;
            fluxL[1] = 0.5 * ( rvnL * uL + xfn[iFace] * pL );
            fluxL[2] = 0.5 * ( rvnL * vL + yfn[iFace] * pL );
            fluxL[3] = 0.5 * ( rvnL * wL + zfn[iFace] * pL );
            fluxL[4] = 0.5 * ( rvnL * hL + vgn[iFace] * pL );

            fluxR[0 ] = 0.5 *   rvnR                       ;
            fluxR[1] = 0.5 * ( rvnR * uR + xfn[iFace] * pR );
            fluxR[2] = 0.5 * ( rvnR * vR + yfn[iFace] * pR );
            fluxR[3] = 0.5 * ( rvnR * wR + zfn[iFace] * pR );
            fluxR[4] = 0.5 * ( rvnR * hR + vgn[iFace] * pR );

            //! 1/2 * ( F(QL) + F(QR) ).Ç°ÃæÒÑ¾­³ËÒÔ1/2
            for ( int m = 0; m < nl; ++ m )
            {
                    flux[m*SEG_LEN+iFace] = 0.5 * ( fluxL[m] + fluxR[m] );
            }
            */
            //!µ½´Ë£¬Í¨Á¿Çó½âÍê±Ï£¬ÏÂÃæ¿ªÊ¼¼ÆËãºÄÉ¢Ïî
            //!##########################################################################################//!

            //!ÏÂÃæÖ÷Òª¼ÆËãºÄÉ¢Ïî

            //!ÏÈ½øÐÐRoeÆ½¾ù Roe-average operation£¬using density, three velocity
            //! components and pressure
            ratio = sqrt(rR / rL);
            //! test
            //! printf("On device, ratio= %f\n", ratio);
            coef = 1.0 / (1.0 + ratio);

            rRoeAverage = sqrt(rL * rR);
            uRoeAverage = (uL + uR * ratio) * coef;
            vRoeAverage = (vL + vR * ratio) * coef;
            wRoeAverage = (wL + wR * ratio) * coef;
            pRoeAverage = (pL + pR * ratio) * coef;

            gama = (gmL + gmR * ratio) * coef; //! 0.5 * ( gamal[iFace] + gamar[iFace]
                                               //! );
            gamm1             = gama - 1.0;
            vSquareRoeAverage = uRoeAverage * uRoeAverage + vRoeAverage * vRoeAverage + wRoeAverage * wRoeAverage;
            hRoeAverage =
                gama / gamm1 * pRoeAverage / rRoeAverage + 0.5 * vSquareRoeAverage; //!( hL + hR * ratio ) * coef;

            vn    = xsn * uRoeAverage + ysn * vRoeAverage + zsn * wRoeAverage - vb;
            theta = vn + vb;

            //! eigenvalue
            c2 = gamm1 * (hRoeAverage - 0.5 * vSquareRoeAverage); //! sound speed gama * pRoeAverage /
                                                                  //! rRoeAverage   ;
            cm = sqrt(GPUABS(c2));

            eigv1 = GPUABS(vn);

            if (precon == 0)
            {
                eigv2 = GPUABS(vn + cm);
                eigv3 = GPUABS(vn - cm);
            }
            else
            {
                beta  = GPUMIN(GPUMAX(vSquareRoeAverage / c2, 3.0 * mach2), 1.0);
                vnew  = 0.5 * vn * (1 + beta);
                cnew  = 0.5 * sqrt(((beta - 1) * vn) * ((beta - 1) * vn) + 4 * beta * c2);
                cnew2 = cnew * cnew;
                eigv2 = GPUABS(vnew + cnew);
                eigv3 = GPUABS(vnew - cnew);
            }
            //!Îªns·½³ÌÌí¼ÓµÄdq
            EL = (1.0 / (gmL - 1.0)) * (pL / rL) + 0.5 * vSquareL;
            ER = (1.0 / (gmR - 1.0)) * (pR / rR) + 0.5 * vSquareR;

            //! qConserveL[0] = rL; qConserveL[1] = rL*uL; qConserveL[2] = rL*vL;
            //! qConserveL[3] = rL*wL; qConserveL[4] = rL*EL; qConserveR[0] = rR;
            //! qConserveR[1] = rR*uR; qConserveR[2] = rR*vR; qConserveR[3] = rR*wR;
            //! qConserveR[4] = rR*ER;
            dq[0] = rR - rL;
            dq[1] = rR * uR - rL * uL;
            dq[2] = rR * vR - rL * vL;
            dq[3] = rR * wR - rL * wL;
            dq[4] = rR * ER - rL * EL;
            //! chemeical part is not suppoerted in the kernel by sunxu 20191105
            /*
            if ( nchem != 0 )
            {
                    Primitive2Conservative(primL, gmL, qConserveL);
                    Primitive2Conservative(primR, gmR, qConserveR);
            }
            */
            //! Ö±½ÓÐ´³öÀ´
            for (int m = nm; m < nl; ++m)
            {
                dq[m] = qConserveR[m] - qConserveL[m];
            }

            //!        cL = sqrt( gamal[iFace] * pL / rL );
            //!       cR = sqrt( gamar[iFace] * pR / rR );

            //! Entropy fix
            if (RoeEntropyFixMethod == 1)
            {
                eigv_max = GPUMAX(eigv2, eigv3);
                tmp0     = eigv_max * alf_l;
                tmp1     = eigv_max * alf_n;
                eigv1    = GPUMAX(tmp0, eigv1);
                eigv2    = GPUMAX(tmp1, eigv2);
                eigv3    = GPUMAX(tmp1, eigv3);
            }
            //! bell 20130326 add
            else if (RoeEntropyFixMethod == 2)
            {
                //! 2: org Harten's method,
            }
            else if (RoeEntropyFixMethod == 3)
            {
                //! 3: Method in ZYB 's Doctor thesis.
            }
            else if (RoeEntropyFixMethod == 4)
            {
                //! 4: Method in Ustar
                eigv_max = GPUMAX(eigv2, eigv3);

                ratio = GPUMIN(1.0, fabs(vn) / cm);
                tmp0  = GPUMAX(eigv_max * alf_l * ratio, eigv_max * alf_l / 100.0);
                tmp1  = GPUMAX(eigv_max * alf_n * ratio, eigv_max * alf_n / 100.0);

                eigv1 = GPUMAX(tmp0, eigv1);
                eigv2 = GPUMAX(tmp1, eigv2);
                eigv3 = GPUMAX(tmp1, eigv3);
            }
            else
            {
                //! TK_Exit::ExceptionExit("Error, no entrfix method is adopt! \n");
            }

            if (precon == 0)
            {
                xi1 = (two * eigv1 - eigv2 - eigv3) / (two * c2);
                xi2 = (eigv2 - eigv3) / (two * cm);
            }
            else
            {
                xi1   = (two * eigv1 - eigv2 - eigv3) / (two * cnew2);
                xi2   = (eigv2 - eigv3) / (two * cnew);
                alpha = 0.5 * (1 - beta);
                xi3   = two * alpha * eigv1 / cnew2;
                xi4   = alpha * vn * xi2 / cnew2;
            }
            //! prim[0] = rRoeAverage; prim[1] = uRoeAverage; prim[2] = vRoeAverage;
            //! prim[3] = wRoeAverage; prim[IP] = pRoeAverage;

            /*
    for ( int m = nm; m < neqn; ++ m )
    {
                prim[m] = 0.5 * ( primL[m] + primR[m] );
    }
    */
            //! ae   = gama - 1.0;
            dc = theta * dq[0] - xsn * dq[1] - ysn * dq[2] - zsn * dq[3];
            if (precon == 0)
            {
                c2dc = c2 * dc;
            }
            else
            {
                Ur2  = beta * c2;
                c2dc = cnew2 * dc;
            }
            //!        vSquare = prim[1]*prim[1] + prim[2]*prim[2] + prim[3]*prim[3];
            //!        hRoeAverage = gama/gamm1 * pRoeAverage/rRoeAverage + 0.5 *
            //! vSquareRoeAverage;

            dh = -gamm1 * (uRoeAverage * dq[1] + vRoeAverage * dq[2] + wRoeAverage * dq[3] - dq[4])
                 + 0.5 * gamm1 * vSquareRoeAverage * dq[0];
            //! chemical part is not supported in the kernel by sunxu 20191105
            /*
            if ( nchem != 0 )
            {
                    ComputeDHAndTotalEnthalpy(prim, gama, dq, dh, hRoeAverage);
            }
            */
            if (precon == 0)
            {
                flux[0 * SEG_LEN + iFace] -= 0.5 * (eigv1 * dq[0] - dh * xi1 - dc * xi2);
                flux[1 * SEG_LEN + iFace] -=
                    0.5 * (eigv1 * dq[1] + (xsn * c2dc - uRoeAverage * dh) * xi1 + (xsn * dh - uRoeAverage * dc) * xi2);
                flux[2 * SEG_LEN + iFace] -=
                    0.5 * (eigv1 * dq[2] + (ysn * c2dc - vRoeAverage * dh) * xi1 + (ysn * dh - vRoeAverage * dc) * xi2);
                flux[3 * SEG_LEN + iFace] -=
                    0.5 * (eigv1 * dq[3] + (zsn * c2dc - wRoeAverage * dh) * xi1 + (zsn * dh - wRoeAverage * dc) * xi2);
                flux[4 * SEG_LEN + iFace] -=
                    0.5
                    * (eigv1 * dq[4] + (theta * c2dc - hRoeAverage * dh) * xi1 + (theta * dh - hRoeAverage * dc) * xi2);

                for (int m = nm; m < nl; ++m)
                {
                    prim[m] = 0.5 * (primL[m] + primR[m]);
                    flux[m * SEG_LEN + iFace] -= 0.5 * (eigv1 * dq[m] - prim[m] * (dh * xi1 + dc * xi2));
                }
            }
            else
            {
                flux[0 * SEG_LEN + iFace] -= 0.5
                                             * (eigv1 * dq[0] - dh * cnew2 / Ur2 * xi1 - dc * xi2
                                                + dh * cnew2 / Ur2 * xi3 - dh * cnew2 / Ur2 * xi4);
                flux[1 * SEG_LEN + iFace] -=
                    0.5
                    * (eigv1 * dq[1] + (xsn * c2dc - uRoeAverage * dh * cnew2 / Ur2) * xi1
                       + (xsn * dh - uRoeAverage * dc) * xi2 + uRoeAverage * dh * cnew2 / Ur2 * xi3
                       - (xsn * c2dc + uRoeAverage * dh * cnew2 / Ur2) * xi4);
                flux[2 * SEG_LEN + iFace] -=
                    0.5
                    * (eigv1 * dq[2] + (ysn * c2dc - vRoeAverage * dh * cnew2 / Ur2) * xi1
                       + (ysn * dh - vRoeAverage * dc) * xi2 + vRoeAverage * dh * cnew2 / Ur2 * xi3
                       - (ysn * c2dc + vRoeAverage * dh * cnew2 / Ur2) * xi4);
                flux[3 * SEG_LEN + iFace] -=
                    0.5
                    * (eigv1 * dq[3] + (zsn * c2dc - wRoeAverage * dh * cnew2 / Ur2) * xi1
                       + (zsn * dh - wRoeAverage * dc) * xi2 + wRoeAverage * dh * cnew2 / Ur2 * xi3
                       - (zsn * c2dc + wRoeAverage * dh * cnew2 / Ur2) * xi4);
                flux[4 * SEG_LEN + iFace] -=
                    0.5
                    * (eigv1 * dq[4] + (theta * c2dc - hRoeAverage * dh * cnew2 / Ur2) * xi1
                       + (theta * dh - hRoeAverage * dc) * xi2 + hRoeAverage * dh * cnew2 / Ur2 * xi3
                       - (theta * c2dc + hRoeAverage * dh * cnew2 / Ur2) * xi4);

                for (int m = nm; m < nl; ++m)
                {
                    prim[m] = 0.5 * (primL[m] + primR[m]);
                    flux[m * SEG_LEN + iFace] -=
                        0.5 * (eigv1 * dq[m] - prim[m] * (dh * cnew2 / Ur2 * (xi1 - xi3 + xi4) + dc * xi2));
                }
            }

            /*Ç°ÃæÒÑ¾­³ËÒÔ1/2´Ë´¦²»ÓÃÔÙÑ­»·
            for ( int m = 0; m < nl; ++ m )
            {
                flux[m][iFace] *= 0.5;
            }
    */
        }

        delete[] prim;
        delete[] primL;
        delete[] primR;
        delete[] qConserveL;
        delete[] qConserveR;
        delete[] dq;
        delete[] fluxL;
        delete[] fluxR;
    }

    //! kernel function for Roe_Scheme_Old
    __global__ void GPURoe_Scheme_Old_S1(RFloat *d_ql_ns, RFloat *d_qr_ns, RDouble *xfn, RDouble *yfn, RDouble *zfn,
                                         RDouble *vgn, RFloat *gamal, RFloat *gamar, RFloat *flux, const int neqn,
                                         const int nl, const int nlen, const int nm, const int RoeEntropyFixMethod,
                                         const int nchem, const int SEG_LEN, const int nst, const int ned,
                                         const int precon, RFloat alf_l, RFloat alf_n, RFloat mach2)
    {
        RFloat gama, gamm1;
        RFloat vSquareL, vnL, hint_L, hL, rvnL, EL, gmL;
        RFloat vSquareR, vnR, hint_R, hR, rvnR, ER, gmR;
        RFloat rRoeAverage, uRoeAverage, vRoeAverage, wRoeAverage, pRoeAverage, vSquareRoeAverage, hRoeAverage, vn, c2,
            cm;
        RFloat Ur2, cnew2, cnew, vnew, beta, alpha;
        RFloat ratio, coef;
        RFloat tmp0, tmp1;
        RFloat eigv1, eigv2, eigv3, eigv_max;           //! cL, cR;
        RFloat theta, xi1, xi2, xi3, xi4, dh, dc, c2dc; //! ae,
        RFloat xsn, ysn, zsn, vb;
        //! use shared memory for the temporary variables
        extern __shared__ RFloat array[];
        RFloat                  *prim       = array;
        RFloat                  *qConserveL = array + neqn * blockDim.x;
        RFloat                  *qConserveR = array + neqn * blockDim.x * 2;
        RFloat                  *dq         = array + neqn * blockDim.x * 3;

        RFloat *primL = array + neqn * blockDim.x * 4;
        RFloat *primR = array + neqn * blockDim.x * 5;
        RFloat *fluxL = array + neqn * blockDim.x * 6;
        RFloat *fluxR = array + neqn * blockDim.x * 7;

        const RDouble LITTLE = 1.0E-10;
        /*
  RFloat * prim       = new RFloat[neqn];
  RFloat * qConserveL = new RFloat[nl];
  RFloat * qConserveR = new RFloat[nl];
  RFloat * dq         = new RFloat[nl];

  RFloat * primL = new RFloat[neqn];
  RFloat * primR = new RFloat[neqn];
  RFloat * fluxL = new RFloat[nl];
  RFloat * fluxR = new RFloat[nl];
  */
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int iFace = bidx * blockDim.x + tidx; iFace < (ned - nst); iFace += blockDim.x * gridDim.x)
        {
            for (int m = 0; m < neqn; ++m)
            {
                primL[m * blockDim.x + tidx] = d_ql_ns[m * SEG_LEN + iFace];
                primR[m * blockDim.x + tidx] = d_qr_ns[m * SEG_LEN + iFace];
            }

            RFloat &rL = primL[0 * blockDim.x + tidx]; //! ql[0][iFace];
            RFloat &uL = primL[1 * blockDim.x + tidx]; //! ql[1][iFace];
            RFloat &vL = primL[2 * blockDim.x + tidx]; //! ql[2][iFace];
            RFloat &wL = primL[3 * blockDim.x + tidx]; //! ql[3][iFace];
            RFloat &pL = primL[4 * blockDim.x + tidx]; //! ql[IP][iFace];

            RFloat &rR = primR[0 * blockDim.x + tidx]; //! qr[0][iFace];
            RFloat &uR = primR[1 * blockDim.x + tidx]; //! qr[1][iFace];
            RFloat &vR = primR[2 * blockDim.x + tidx]; //! qr[2][iFace];
            RFloat &wR = primR[3 * blockDim.x + tidx]; //! qr[3][iFace];
            RFloat &pR = primR[4 * blockDim.x + tidx]; //! qr[IP][iFace];

            xsn = xfn[nst + iFace];
            ysn = yfn[nst + iFace];
            zsn = zfn[nst + iFace];
            vb  = vgn[nst + iFace];
            gmL = gamal[iFace];
            gmR = gamar[iFace];

            rL = GPUMAX(rL, LITTLE);
            rR = GPUMAX(rR, LITTLE);
            pL = GPUMAX(pL, LITTLE);
            pR = GPUMAX(pR, LITTLE);

            //!###################################################################//!
            //! to get the flux ÓÉ×óÓÒ½çÃæµÄÖµÇóÍ¨Á¿

            vSquareL = uL * uL + vL * vL + wL * wL;
            vSquareR = uR * uR + vR * vR + wR * wR;

            hint_L = (gmL / (gmL - 1.0)) * (pL / rL);
            hint_R = (gmR / (gmR - 1.0)) * (pR / rR);
            //! chemical part is not supported in the kernel by sunxu 20191105
            /*
    if ( nchem != 0 )
    {
            ComputeEnthalpyByPrimitive(primL, gmL, hint_L);
            ComputeEnthalpyByPrimitive(primR, gmR, hint_R);
    }
    */
            hL = hint_L + 0.5 * vSquareL;
            hR = hint_R + 0.5 * vSquareR;

            vnL  = xsn * uL + ysn * vL + zsn * wL - vb;
            vnR  = xsn * uR + ysn * vR + zsn * wR - vb;
            rvnL = rL * vnL;
            rvnR = rR * vnR;

            flux[0 * SEG_LEN + iFace] = 0.5 * (rvnL + rvnR);
            flux[1 * SEG_LEN + iFace] = 0.5 * (rvnL * uL + xsn * pL + rvnR * uR + xsn * pR);
            flux[2 * SEG_LEN + iFace] = 0.5 * (rvnL * vL + ysn * pL + rvnR * vR + ysn * pR);
            flux[3 * SEG_LEN + iFace] = 0.5 * (rvnL * wL + zsn * pL + rvnR * wR + zsn * pR);
            flux[4 * SEG_LEN + iFace] = 0.5 * (rvnL * hL + vb * pL + rvnR * hR + vb * pR);

            for (int m = nm; m < nl; ++m)
            {
                //!            fluxL[m*blockDim.x+tidx] = 0.5 * (
                //! primL[m*blockDim.x+tidx] * fluxL[0 ] );
                //! fluxR[m*blockDim.x+tidx] = 0.5 * ( primR[m*blockDim.x+tidx] * fluxR[0 ]
                //!);
                //! This loop is useless
                //!flux[m * SEG_LEN + iFace] = 0.5 * (primL[m * blockDim.x + tidx] * rvnL +
                //!                                  primR[m * blockDim.x + tidx] * rvnR);
            }

            /*  Ö±½Ó¼Óµ½Í¨Á¿
    fluxL[0 ] = 0.5 *   rvnL                         ;
    fluxL[1] = 0.5 * ( rvnL * uL + xfn[iFace] * pL );
    fluxL[2] = 0.5 * ( rvnL * vL + yfn[iFace] * pL );
    fluxL[3] = 0.5 * ( rvnL * wL + zfn[iFace] * pL );
    fluxL[4] = 0.5 * ( rvnL * hL + vgn[iFace] * pL );

    fluxR[0 ] = 0.5 *   rvnR                       ;
    fluxR[1] = 0.5 * ( rvnR * uR + xfn[iFace] * pR );
    fluxR[2] = 0.5 * ( rvnR * vR + yfn[iFace] * pR );
    fluxR[3] = 0.5 * ( rvnR * wR + zfn[iFace] * pR );
    fluxR[4] = 0.5 * ( rvnR * hR + vgn[iFace] * pR );

    //! 1/2 * ( F(QL) + F(QR) ).Ç°ÃæÒÑ¾­³ËÒÔ1/2
    for ( int m = 0; m < nl; ++ m )
    {
            flux[m*SEG_LEN+iFace] = 0.5 * ( fluxL[m*blockDim.x+tidx] +
    fluxR[m*blockDim.x+tidx] );
    }
    */
            //!µ½´Ë£¬Í¨Á¿Çó½âÍê±Ï£¬ÏÂÃæ¿ªÊ¼¼ÆËãºÄÉ¢Ïî
            //!##########################################################################################//!

            //!ÏÂÃæÖ÷Òª¼ÆËãºÄÉ¢Ïî

            //!ÏÈ½øÐÐRoeÆ½¾ù Roe-average operation£¬using density, three velocity
            //! components and pressure
            ratio = sqrt(rR / rL);
            //! test
            //! printf("On device, ratio= %f\n", ratio);
            coef = 1.0 / (1.0 + ratio);

            rRoeAverage = sqrt(rL * rR);
            uRoeAverage = (uL + uR * ratio) * coef;
            vRoeAverage = (vL + vR * ratio) * coef;
            wRoeAverage = (wL + wR * ratio) * coef;
            pRoeAverage = (pL + pR * ratio) * coef;

            gama = (gmL + gmR * ratio) * coef; //! 0.5 * ( gamal[iFace] + gamar[iFace]
                                               //! );
            gamm1             = gama - 1.0;
            vSquareRoeAverage = uRoeAverage * uRoeAverage + vRoeAverage * vRoeAverage + wRoeAverage * wRoeAverage;
            hRoeAverage =
                gama / gamm1 * pRoeAverage / rRoeAverage + 0.5 * vSquareRoeAverage; //!( hL + hR * ratio ) * coef;

            vn    = xsn * uRoeAverage + ysn * vRoeAverage + zsn * wRoeAverage - vb;
            theta = vn + vb;

            //! eigenvalue
            c2 = gamm1 * (hRoeAverage - 0.5 * vSquareRoeAverage); //! sound speed gama * pRoeAverage /
                                                                  //! rRoeAverage   ;
            cm = sqrt(GPUABS(c2));

            eigv1 = GPUABS(vn);

            if (precon == 0)
            {
                eigv2 = GPUABS(vn + cm);
                eigv3 = GPUABS(vn - cm);
            }
            else
            {
                beta  = GPUMIN(GPUMAX(vSquareRoeAverage / c2, 3.0 * mach2), 1.0);
                vnew  = 0.5 * vn * (1 + beta);
                cnew  = 0.5 * sqrt(((beta - 1) * vn) * ((beta - 1) * vn) + 4 * beta * c2);
                cnew2 = cnew * cnew;
                eigv2 = GPUABS(vnew + cnew);
                eigv3 = GPUABS(vnew - cnew);
            }
            //!Îªns·½³ÌÌí¼ÓµÄdq
            EL = (1.0 / (gmL - 1.0)) * (pL / rL) + 0.5 * vSquareL;
            ER = (1.0 / (gmR - 1.0)) * (pR / rR) + 0.5 * vSquareR;

            //! qConserveL[0] = rL; qConserveL[1] = rL*uL; qConserveL[2] = rL*vL;
            //! qConserveL[3] = rL*wL; qConserveL[4] = rL*EL; qConserveR[0] = rR;
            //! qConserveR[1] = rR*uR; qConserveR[2] = rR*vR; qConserveR[3] = rR*wR;
            //! qConserveR[4] = rR*ER;
            dq[0 * blockDim.x + tidx] = rR - rL;
            dq[1 * blockDim.x + tidx] = rR * uR - rL * uL;
            dq[2 * blockDim.x + tidx] = rR * vR - rL * vL;
            dq[3 * blockDim.x + tidx] = rR * wR - rL * wL;
            dq[4 * blockDim.x + tidx] = rR * ER - rL * EL;
            //! chemeical part is not suppoerted in the kernel by sunxu 20191105
            /*
    if ( nchem != 0 )
    {
            Primitive2Conservative(primL, gmL, qConserveL);
            Primitive2Conservative(primR, gmR, qConserveR);
    }
    */
            //! Ö±½ÓÐ´³öÀ´
            for (int m = nm; m < nl; ++m)
            {
                dq[m * blockDim.x + tidx] = qConserveR[m * blockDim.x + tidx] - qConserveL[m * blockDim.x + tidx];
            }

            //!        cL = sqrt( gamal[iFace] * pL / rL );
            //!       cR = sqrt( gamar[iFace] * pR / rR );

            //! Entropy fix
            if (RoeEntropyFixMethod == 1)
            {
                eigv_max = GPUMAX(eigv2, eigv3);
                tmp0     = eigv_max * alf_l;
                tmp1     = eigv_max * alf_n;
                eigv1    = GPUMAX(tmp0, eigv1);
                eigv2    = GPUMAX(tmp1, eigv2);
                eigv3    = GPUMAX(tmp1, eigv3);
            }
            //! bell 20130326 add
            else if (RoeEntropyFixMethod == 2)
            {
                //! 2: org Harten's method,
            }
            else if (RoeEntropyFixMethod == 3)
            {
                //! 3: Method in ZYB 's Doctor thesis.
            }
            else if (RoeEntropyFixMethod == 4)
            {
                //! 4: Method in Ustar
                eigv_max = GPUMAX(eigv2, eigv3);

                ratio = GPUMIN(1.0, fabs(vn) / cm);
                tmp0  = GPUMAX(eigv_max * alf_l * ratio, eigv_max * alf_l / 100.0);
                tmp1  = GPUMAX(eigv_max * alf_n * ratio, eigv_max * alf_n / 100.0);

                eigv1 = GPUMAX(tmp0, eigv1);
                eigv2 = GPUMAX(tmp1, eigv2);
                eigv3 = GPUMAX(tmp1, eigv3);
            }
            else
            {
                //! TK_Exit::ExceptionExit("Error, no entrfix method is adopt! \n");
            }

            if (precon == 0)
            {
                xi1 = (two * eigv1 - eigv2 - eigv3) / (two * c2);
                xi2 = (eigv2 - eigv3) / (two * cm);
            }
            else
            {
                xi1   = (two * eigv1 - eigv2 - eigv3) / (two * cnew2);
                xi2   = (eigv2 - eigv3) / (two * cnew);
                alpha = 0.5 * (1 - beta);
                xi3   = two * alpha * eigv1 / cnew2;
                xi4   = alpha * vn * xi2 / cnew2;
            }
            //! prim[0] = rRoeAverage; prim[1] = uRoeAverage; prim[2] = vRoeAverage;
            //! prim[3] = wRoeAverage; prim[IP] = pRoeAverage;

            /*
    for ( int m = nm; m < neqn; ++ m )
    {
                prim[m*blockDim.x+tidx] = 0.5 * ( primL[m*blockDim.x+tidx] +
    primR[m*blockDim.x+tidx] );
    }
    */
            //! ae   = gama - 1.0;
            dc = theta * dq[0 * blockDim.x + tidx] - xsn * dq[1 * blockDim.x + tidx] - ysn * dq[2 * blockDim.x + tidx]
                 - zsn * dq[3 * blockDim.x + tidx];
            c2dc = 0.0;
            if (precon == 0)
            {
                c2dc = c2 * dc;
            }
            else
            {
                Ur2  = beta * c2;
                c2dc = cnew2 * dc;
            }
            //!        vSquare = prim[1]*prim[1] + prim[2]*prim[2] + prim[3]*prim[3];
            //!        hRoeAverage = gama/gamm1 * pRoeAverage/rRoeAverage + 0.5 *
            //! vSquareRoeAverage;

            dh = -gamm1
                     * (uRoeAverage * dq[1 * blockDim.x + tidx] + vRoeAverage * dq[2 * blockDim.x + tidx]
                        + wRoeAverage * dq[3 * blockDim.x + tidx] - dq[4 * blockDim.x + tidx])
                 + 0.5 * gamm1 * vSquareRoeAverage * dq[0 * blockDim.x + tidx];
            //! chemical part is not supported in the kernel by sunxu 20191105
            /*
    if ( nchem != 0 )
    {
            ComputeDHAndTotalEnthalpy(prim, gama, dq, dh, hRoeAverage);
    }
    */
            if (precon == 0)
            {
                flux[0 * SEG_LEN + iFace] -= 0.5 * (eigv1 * dq[0 * blockDim.x + tidx] - dh * xi1 - dc * xi2);
                flux[1 * SEG_LEN + iFace] -=
                    0.5
                    * (eigv1 * dq[1 * blockDim.x + tidx] + (xsn * c2dc - uRoeAverage * dh) * xi1
                       + (xsn * dh - uRoeAverage * dc) * xi2);
                flux[2 * SEG_LEN + iFace] -=
                    0.5
                    * (eigv1 * dq[2 * blockDim.x + tidx] + (ysn * c2dc - vRoeAverage * dh) * xi1
                       + (ysn * dh - vRoeAverage * dc) * xi2);
                flux[3 * SEG_LEN + iFace] -=
                    0.5
                    * (eigv1 * dq[3 * blockDim.x + tidx] + (zsn * c2dc - wRoeAverage * dh) * xi1
                       + (zsn * dh - wRoeAverage * dc) * xi2);
                flux[4 * SEG_LEN + iFace] -=
                    0.5
                    * (eigv1 * dq[4 * blockDim.x + tidx] + (theta * c2dc - hRoeAverage * dh) * xi1
                       + (theta * dh - hRoeAverage * dc) * xi2);

                for (int m = nm; m < nl; ++m)
                {
                    prim[m * blockDim.x + tidx] = 0.5 * (primL[m * blockDim.x + tidx] + primR[m * blockDim.x + tidx]);
                    flux[m * SEG_LEN + iFace] -=
                        0.5 * (eigv1 * dq[m * blockDim.x + tidx] - prim[m * blockDim.x + tidx] * (dh * xi1 + dc * xi2));
                }
            }
            else
            {
                flux[0 * SEG_LEN + iFace] -= 0.5
                                             * (eigv1 * dq[0 * blockDim.x + tidx] - dh * cnew2 / Ur2 * xi1 - dc * xi2
                                                + dh * cnew2 / Ur2 * xi3 - dh * cnew2 / Ur2 * xi4);
                flux[1 * SEG_LEN + iFace] -=
                    0.5
                    * (eigv1 * dq[1 * blockDim.x + tidx] + (xsn * c2dc - uRoeAverage * dh * cnew2 / Ur2) * xi1
                       + (xsn * dh - uRoeAverage * dc) * xi2 + uRoeAverage * dh * cnew2 / Ur2 * xi3
                       - (xsn * c2dc + uRoeAverage * dh * cnew2 / Ur2) * xi4);
                flux[2 * SEG_LEN + iFace] -=
                    0.5
                    * (eigv1 * dq[2 * blockDim.x + tidx] + (ysn * c2dc - vRoeAverage * dh * cnew2 / Ur2) * xi1
                       + (ysn * dh - vRoeAverage * dc) * xi2 + vRoeAverage * dh * cnew2 / Ur2 * xi3
                       - (ysn * c2dc + vRoeAverage * dh * cnew2 / Ur2) * xi4);
                flux[3 * SEG_LEN + iFace] -=
                    0.5
                    * (eigv1 * dq[3 * blockDim.x + tidx] + (zsn * c2dc - wRoeAverage * dh * cnew2 / Ur2) * xi1
                       + (zsn * dh - wRoeAverage * dc) * xi2 + wRoeAverage * dh * cnew2 / Ur2 * xi3
                       - (zsn * c2dc + wRoeAverage * dh * cnew2 / Ur2) * xi4);
                flux[4 * SEG_LEN + iFace] -=
                    0.5
                    * (eigv1 * dq[4 * blockDim.x + tidx] + (theta * c2dc - hRoeAverage * dh * cnew2 / Ur2) * xi1
                       + (theta * dh - hRoeAverage * dc) * xi2 + hRoeAverage * dh * cnew2 / Ur2 * xi3
                       - (theta * c2dc + hRoeAverage * dh * cnew2 / Ur2) * xi4);

                for (int m = nm; m < nl; ++m)
                {
                    prim[m * blockDim.x + tidx] = 0.5 * (primL[m * blockDim.x + tidx] + primR[m * blockDim.x + tidx]);
                    flux[m * SEG_LEN + iFace] -=
                        0.5
                        * (eigv1 * dq[m * blockDim.x + tidx]
                           - prim[m * blockDim.x + tidx] * (dh * cnew2 / Ur2 * (xi1 - xi3 + xi4) + dc * xi2));
                }
            }

            /*Ç°ÃæÒÑ¾­³ËÒÔ1/2´Ë´¦²»ÓÃÔÙÑ­»·
            for ( int m = 0; m < nl; ++ m )
            {
                flux[m][iFace] *= 0.5;
            }
    */
        }
        /*
  delete [] prim;
  delete [] primL;
  delete [] primR;
  delete [] qConserveL;
  delete [] qConserveR;
  delete [] dq;
  delete [] fluxL;
  delete [] fluxR;
  */
    }
    //! call kernel for GPUInviscidFluxWrap's first for loop
    void CallGPUInviscidFluxWrapLoop1(const int nl, const int nst, const int ned, const int SEG_LEN)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = ned - nst;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUInviscidFluxWrapLoop1, loopLen, 0, gridSize, blockSize, blockSizeLimit);
        //! adjust gridSize
        gridSize = 50;
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUInviscidFluxWrapLoop1, 0, gridSize, blockSize);
#endif
        GPUInviscidFluxWrapLoop1<<<gridSize, blockSize>>>(d_flux, d_area, nl, nst, ned, SEG_LEN);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    //! kernel for GPUInviscidFluxWrap's first for loop
    __global__ void GPUInviscidFluxWrapLoop1(RFloat *d_flux, RDouble *d_area, const int nl, const int nst,
                                             const int ned, const int SEG_LEN)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        for (int iFace = nst + bidx * blockDim.x + tidx; iFace < ned; iFace += blockDim.x * gridDim.x)
        {
            int j = iFace - nst;

            RDouble &areas = d_area[iFace];

            for (int m = 0; m < nl; ++m)
            {
                d_flux[m * SEG_LEN + j] *= areas;
            }
        }
    }
    //! call kernel for LoadFlux
    void CallGPULoadFlux(Grid *grid_in, FaceProxy *faceProxy, int nst, int ned)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        UnstructGrid *grid       = UnstructGridCast(grid_in);
        int           nBoundFace = grid->GetNBoundFace();
        int           nTotalCell = grid->GetNTotalCell();
        int           nTotal     = nTotalCell + nBoundFace;

        int nl = GlobalDataBase::GetIntParaFromDB("nl");

        //! Determine if there are boundary faces.
        int nMid = nst;
        if (ned <= nBoundFace)
        {
            //! If all boundary faces
            nMid = ned;
        }
        else if (nst < nBoundFace)
        {
            //! Part of them are boundary faces
            nMid = nBoundFace;
        }

        //! call kernel for GPULoadFluxLoop1
        CallGPULoadFluxLoop1(nl, nTotal, d_SEG_LEN, nst, nMid);

        //! call kernel for GPULoadFluxLoop2
        CallGPULoadFluxLoop2(nl, nTotal, d_SEG_LEN, nst, nMid, ned);
    }
    //! call kernel for GPULoadFluxLoop1
    void CallGPULoadFluxLoop1(const int nl, const int nTotal, const int SEG_LEN, const int nst, const int nMid)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize  = 1;
        int blockSize = 1;
        int loopLen   = nMid - nst;
        if (0 == loopLen) return;
        KernelLaunchPara((void *)GPULoadFluxLoop1, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPULoadFluxLoop1, 0, gridSize, blockSize);
#endif

        //! printf("In CallGPULoadFluxLoop1, nst= %d, nMid= %d\n", nst, nMid);
#ifndef CUDADEBUGATOMIC
        GPULoadFluxLoop1<<<gridSize, blockSize>>>(d_res_ns, d_flux, d_left_cell_of_face, nl, nTotal, SEG_LEN, nst,
                                                  nMid);
#else
        GPULoadFluxLoop1<<<1, 1>>>(d_res_ns, d_flux, d_left_cell_of_face, nl, nTotal, SEG_LEN, nst, nMid);
#endif
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    __global__ void GPULoadFluxLoop1(RFloat *d_res_ns, RFloat *d_flux, int *d_left_cell_of_face, const int nl,
                                     const int nTotal, const int SEG_LEN, const int nst, const int nMid)
    {
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        //! For boundary faces, remeber re is ghost cell
        for (int i = nst + bidx * blockDim.x + tidx; i < nMid; i += blockDim.x * gridDim.x)
        {
            int le, j;
            le = d_left_cell_of_face[i];
            j  = i - nst;
            //! if(le== 44352) printf("le== 44352\n");
            for (int m = 0; m < nl; ++m)
            {
#if __CUDA_ARCH__ < 600
                using namespace GPUKernels;
                atomicAddTest(&(d_res_ns[m * nTotal + le]), -(d_flux[m * SEG_LEN + j]));
#else
                atomicAdd(&(d_res_ns[m * nTotal + le]), -(d_flux[m * SEG_LEN + j]));
#endif
                //! d_res_ns[m*nTotal+le] -= d_flux[m*SEG_LEN+j];
                //!  test
                //! if(m== 4 && le== 44352) printf("Loop1, d_flux[%d]= %f, d_res_ns[%d]=
                //! %f\n",m*SEG_LEN+j, d_flux[m*SEG_LEN+j], m*nTotal+le,
                //! d_res_ns[m*nTotal+le]);
            }
        }
    }
    //! call kernel for GPULoadFluxLoop2
    void CallGPULoadFluxLoop2(const int nl, const int nTotal, const int SEG_LEN, const int nst, const int nMid,
                              const int ned)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize  = 1;
        int blockSize = 1;
        int loopLen   = ned - nMid;
        if (0 == loopLen) return;
        KernelLaunchPara((void *)GPULoadFluxLoop2, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPULoadFluxLoop2, 0, gridSize, blockSize);
#endif
#ifndef CUDADEBUGATOMIC
        GPULoadFluxLoop2<<<gridSize, blockSize>>>(d_res_ns, d_flux, d_left_cell_of_face, d_right_cell_of_face, nl,
                                                  nTotal, SEG_LEN, nst, nMid, ned);
#else
        GPULoadFluxLoop2<<<1, 1>>>(d_res_ns, d_flux, d_left_cell_of_face, d_right_cell_of_face, nl, nTotal, SEG_LEN,
                                   nst, nMid, ned);
#endif
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    __global__ void GPULoadFluxLoop2(RFloat *d_res_ns, RFloat *d_flux, int *d_left_cell_of_face,
                                     int *d_right_cell_of_face, const int nl, const int nTotal, const int SEG_LEN,
                                     const int nst, const int nMid, const int ned)
    {
        //! volatile RFloat *d_res_ns_tmp=;
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;

        //! Interior faces
        for (int i = nMid + bidx * blockDim.x + tidx; i < ned; i += blockDim.x * gridDim.x)
        {
            int le, re, j;
            le = d_left_cell_of_face[i];
            re = d_right_cell_of_face[i];
            j  = i - nst;

            for (int m = 0; m < nl; ++m)
            {
#if __CUDA_ARCH__ < 600
                using namespace GPUKernels;
                atomicAddTest(&(d_res_ns[m * nTotal + le]), -(d_flux[m * SEG_LEN + j]));
                atomicAddTest(&(d_res_ns[m * nTotal + re]), d_flux[m * SEG_LEN + j]);
#else
                atomicAdd(&(d_res_ns[m * nTotal + le]), -(d_flux[m * SEG_LEN + j]));
                //! d_res_ns[m*nTotal+le] -= d_flux[m*SEG_LEN+j];
                atomicAdd(&(d_res_ns[m * nTotal + re]), d_flux[m * SEG_LEN + j]);
                //! d_res_ns[m*nTotal+re] += d_flux[m*SEG_LEN+j];
#endif
            }
        }
    }
    /*
        void CallGPUCompGamaAndTField(double refGama, RFloat
   coefficientOfStateEquation, RFloat omav, const int IR, const int IP, const
   int ITT, const int nTotal)
        {
                using namespace GPUMemory;
                using namespace GPUFlowVariables;
                using namespace GPUGeomVariables;

                int gridSize= 1;
                int blockSize =1;
                int loopLen= nTotal;
                int blockSizeLimit= 0;
                KernelLaunchPara((void*)GPUCompGamaAndTField, loopLen, 0,
   gridSize, blockSize, blockSizeLimit); #ifdef KERNELLAUNCHTEST
                        ReportKernelPara((void*)GPUCompGamaAndTField, 0,
   gridSize, blockSize); #endif GPUCompGamaAndTField<<<gridSize,
   blockSize>>>(d_q_ns, d_gama_ns, d_t_proxy, refGama,
   coefficientOfStateEquation, omav, IR, IP, ITT, nTotal); #ifdef
   KERNELLAUNCHTEST HANDLE_KERNEL_ERR(); #endif

        }
        __global__ void GPUCompGamaAndTField(RFloat *q, RFloat *gama,  RFloat
   *t, double refGama, RFloat coefficientOfStateEquation, RFloat omav, const int
   IR, const int IP, const int ITT, const int nTotal)
        {
                //! omit the initialize of t
                const int tidx= threadIdx.x;
                const int bidx= blockIdx.x;
                for (int icell = bidx * blockDim.x + tidx;
                         icell < nTotal;
                         icell += blockDim.x * gridDim.x )
                {
                        RFloat &rm = q[IR*nTotal+ icell];
                        RFloat &pm = q[IP*nTotal+ icell];
                        gama[icell] = refGama;
                        //!t[ITT*nTotal+ icell] = pm / (
   coefficientOfStateEquation * rm * omav );
               }
        }
*/
    //! void CallGPUCompGamaAndTField_S(double refGama,
    //!                                 RFloat coefficientOfStateEquation, const int nl,
    //!                                 const int nchem, const int nTotal)
    void CallGPUComputeGamaAndTemperature(Grid *grid_in, Param_NSSolverUnstruct *parameters)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        using namespace GAS_SPACE;
        using namespace IDX;

        UnstructGrid *grid = UnstructGridCast(grid_in);

        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal     = nTotalCell + nBoundFace;

        RDouble refGama                    = parameters->GetRefGama();
        RDouble coefficientOfStateEquation = gas->GetCoefficientOfStateEquation();

        d_refGama                    = refGama;
        d_coefficientOfStateEquation = gas->GetCoefficientOfStateEquation();
        d_ntmodel                    = parameters->GetTemperatureModel();

        using namespace PHSPACE;
        using namespace IDX;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotal;
        int blockSizeLimit = 0;

        //!int neqn = nl + nchem;
        RFloat omav = one;
        KernelLaunchPara((void *)GPUCompGamaAndTField_S, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompGamaAndTField_S, 0, gridSize, blockSize);
#endif
        GPUCompGamaAndTField_S<<<gridSize, blockSize>>>(d_q_ns, d_gama_ns, d_t_proxy, refGama,
                                                        coefficientOfStateEquation, omav, IR, IP, ITT, nTotal);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUCompGamaAndTField_S(RFloat *q, RFloat *gama, RFloat *t, double refGama,
                                           RFloat coefficientOfStateEquation, RFloat omav, const int IR, const int IP,
                                           const int ITT, const int nTotal)
    {
        //! omit the initialize of t
        const int tidx = threadIdx.x;
        const int bidx = blockIdx.x;
        for (int icell = bidx * blockDim.x + tidx; icell < nTotal; icell += blockDim.x * gridDim.x)
        {
            RFloat &rm              = q[IR * nTotal + icell];
            RFloat &pm              = q[IP * nTotal + icell];
            gama[icell]             = refGama;
            t[ITT * nTotal + icell] = pm / (coefficientOfStateEquation * rm * omav);
        }
    }

    __global__ void GPUSetFluxSubZero(const int SEG_LEN, RFloat *flux_sub)
    {
        int i     = blockIdx.x * blockDim.x + threadIdx.x;
        int iface = 0;
        for (iface = i; iface < SEG_LEN; iface += blockDim.x * gridDim.x)
        {
            flux_sub[iface] = 0.0;
        }
    }
} //! namespace GPUNSSolverUnstruct

//! check the rho and pressure
__device__ bool GPUPositiveCheck(RFloat *qtry)
{
    if (qtry[0] <= 0.0 || qtry[4] <= 0.0) return false;
    return true;
}
