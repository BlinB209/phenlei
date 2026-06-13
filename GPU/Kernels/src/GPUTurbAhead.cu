#include "GPUTurbAhead.h"
#include "GPUBasicFunctions.h"
#include "GPUDeviceControl.h"
#include "GPUNSSolver.h"

namespace GPUKernels
{
    void CallGPUTurbZeroResiduals(const int numberOfGeneralCells, const int numberOfEquations)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = numberOfGeneralCells;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUTurbZeroResiduals, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUTurbZeroResiduals, 0, gridSize, blockSize);
#endif

        GPUTurbZeroResiduals<<<gridSize, blockSize>>>(numberOfGeneralCells, numberOfEquations, d_res_turb);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUTurbStoreRhsByResidual(const int nTotal, const int n_turb)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotal;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUTurbStoreRhsByResidual, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUTurbStoreRhsByResidual, 0, gridSize, blockSize);
#endif
        GPUTurbStoreRhsByResidual<<<gridSize, blockSize>>>(nTotal, n_turb, d_rhs_turb, d_res_turb);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUTurbLoadQ(const int nTotalCell, const int nTotal, const int n_turb)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUTurbLoadQ, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUTurbLoadQ, 0, gridSize, blockSize);
#endif

        //!GPUTurbLoadQ<<<gridSize, blockSize>>>(nTotalCell, nTotal, n_turb, d_q_turb_proxy, d_q_proxy_turb);
        GPUTurbLoadQ<<<gridSize, blockSize>>>(nTotalCell, nTotal, n_turb, d_q_turb_proxy, d_q_turb_proxy);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUTurbFillField(const int nTotal, const int neqn)
    {
        using namespace GPUMemory;
        using namespace GPUControlVariables;
        using namespace GPUFlowVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotal;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUTurbFillField, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUTurbFillField, 0, gridSize, blockSize);
#endif

        if (d_FillField_turb == "q_proxyToq_proxy_temp")
        {
            //!FillField(grid, q_proxy_temp, q_proxy) in CFDSolver::RungeKutta
            //!GPUTurbFillField<<<gridSize, blockSize>>>(nTotal, neqn, d_q_proxy_temp_turb, d_q_proxy_turb);
            GPUTurbFillField<<<gridSize, blockSize>>>(nTotal, neqn, d_q_proxy_temp_turb, d_q_turb_proxy);
#ifdef KERNELLAUNCHTEST
            HANDLE_KERNEL_ERR();
#endif
        }
        else if (d_FillField_turb == "q_proxy_tempToq_proxy")
        {
            //!FillField(grid, q_proxy, q_proxy_temp) in CFDSolver::RungeKutta
            //!GPUTurbFillField<<<gridSize, blockSize>>>(nTotal, neqn, d_q_proxy_turb, d_q_proxy_temp_turb);
            GPUTurbFillField<<<gridSize, blockSize>>>(nTotal, neqn, d_q_turb_proxy, d_q_proxy_temp_turb);
#ifdef KERNELLAUNCHTEST
            HANDLE_KERNEL_ERR();
#endif
        }
        else if (d_FillField_turb == "resTodq")
        {
            size_t sizeDq = nTotal * neqn * sizeof(RFloat);
            HANDLE_API_ERR(cudaMemcpy(d_dq_turb, d_res_turb, sizeDq, cudaMemcpyDeviceToDevice));
        }
        else
        {
            cout << "Error : d_FillField_turb owns a wrong value " << d_FillField_turb
                 << ", which should be q_proxyToq_proxy_temp or q_proxy_tempToq_proxy" << endl;
            exit(1);
        }
    }

    void CallGPUTurbLoadResiduals(const int nTotal, const int n_turb)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotal;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUTurbLoadResiduals, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUTurbLoadResiduals, 0, gridSize, blockSize);
#endif
        GPUTurbLoadResiduals<<<gridSize, blockSize>>>(nTotal, n_turb, d_res_turb, d_rhs_turb);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUTurbLhs(const int nTotalCell, const int nTotal, const int n_turb, const double coef)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUTurbLhs, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUTurbLhs, 0, gridSize, blockSize);
#endif

        GPUTurbLhs<<<gridSize, blockSize>>>(nTotalCell, nTotal, n_turb, coef, d_dt, d_res_turb);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUTurbUpdateFlowField(const int nTotalCell, const int nTotal, const int n_turb, const double coef,
                                    const double turb_relax, const int nnegtive_max, const int ZoneID,
                                    const double SMALL, const double cv13, const double mytmax, const double my_nu_max)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        GPUSetNegPosZero<<<1, 1>>>(d_n_neg, d_n_pos);
        /*    
        GPUTurbUpdateFlowField<<<1024, 128>>>(nTotalCell, nTotal, n_turb, d_cell2Face, 
                    d_posiCell2Face, d_face_number_of_each_cell, d_left_cell_of_face, 
                    d_right_cell_of_face, d_q_turb_proxy, coef, d_turboo, d_visl, turb_relax,
                    nnegtive_max, ZoneID, SMALL, cv13, mytmax, my_nu_max, d_n_neg, d_n_pos, 
                    d_t_proxy, d_q_ns, d_res_turb, d_q_proxy_turb);
        
        
        */
        int    gridSize        = 1;
        int    blockSize       = 1;
        int    loopLen         = nTotalCell;
        int    regsPerThread   = 64;
        int    residentBlockNo = 4;
        size_t dsMemPerThread  = n_turb * sizeof(RFloat); //! calculate by the intended shared memory use in the kernel
        size_t dsMemPerBlock;
        KernelLaunchPara((void *)GPUTurbUpdateFlowField_S1, loopLen, regsPerThread, residentBlockNo, dsMemPerThread,
                         dsMemPerBlock, gridSize, blockSize, GPUProp);

//!KernelLaunchPara((void*)GPUReConstructFaceValueLoop2, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUTurbUpdateFlowField_S1, dsMemPerBlock, gridSize, blockSize);
#endif
#ifndef CUDADEBUGATOMIC
        GPUTurbUpdateFlowField_S1<<<gridSize, blockSize>>>(
            nTotalCell, nTotal, n_turb, d_cell2Face, d_posiCell2Face, d_face_number_of_each_cell, d_left_cell_of_face,
            d_right_cell_of_face, d_q_turb_proxy, coef, d_turboo, d_visl, turb_relax, nnegtive_max, ZoneID, SMALL, cv13,
            mytmax, my_nu_max, d_n_neg, d_n_pos,
            //!d_t_proxy, d_q_ns, d_res_turb, d_q_proxy_turb);
            d_t_proxy, d_q_ns, d_res_turb, d_q_turb_proxy);
#else
        GPUTurbUpdateFlowField_S1<<<1, 1>>>(nTotalCell, nTotal, n_turb, d_cell2Face, d_posiCell2Face,
                                            d_face_number_of_each_cell, d_left_cell_of_face, d_right_cell_of_face,
                                            d_q_turb_proxy, coef, d_turboo, d_visl, turb_relax, nnegtive_max, ZoneID,
                                            SMALL, cv13, mytmax, my_nu_max, d_n_neg, d_n_pos,
                                            //!d_t_proxy, d_q_ns, d_res_turb, d_q_proxy_turb);
                                            d_t_proxy, d_q_ns, d_res_turb, d_q_turb_proxy);
#endif

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUTurbUpdateLUSGSFlowField(const int nTotalCell, const int nTotal, const int n_turb, const double coef,
                                         const double turb_relax, const int nnegtive_max, const int ZoneID,
                                         const double SMALL, const double cv13, const double mytmax,
                                         const double my_nu_max)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        GPUSetNegPosZero<<<1, 1>>>(d_n_neg, d_n_pos);
        int    gridSize        = 1;
        int    blockSize       = 1;
        int    loopLen         = nTotalCell;
        int    regsPerThread   = 64;
        int    residentBlockNo = 4;
        size_t dsMemPerThread  = n_turb * sizeof(RFloat); //! calculate by the intended shared memory use in the kernel
        size_t dsMemPerBlock;
        KernelLaunchPara((void *)GPUTurbUpdateFlowField_S1, loopLen, regsPerThread, residentBlockNo, dsMemPerThread,
                         dsMemPerBlock, gridSize, blockSize, GPUProp);

#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUTurbUpdateFlowField_S1, dsMemPerBlock, gridSize, blockSize);
#endif
#ifndef CUDADEBUGATOMIC
        GPUTurbUpdateFlowField_S1<<<gridSize, blockSize>>>(
            nTotalCell, nTotal, n_turb, d_cell2Face, d_posiCell2Face, d_face_number_of_each_cell, d_left_cell_of_face,
            d_right_cell_of_face, d_q_turb_proxy, coef, d_turboo, d_visl, turb_relax, nnegtive_max, ZoneID, SMALL, cv13,
            mytmax, my_nu_max, d_n_neg, d_n_pos, d_t_proxy, d_q_ns, d_dq_turb, d_q_turb_proxy);
#else
        GPUTurbUpdateFlowField_S1<<<1, 1>>>(
            nTotalCell, nTotal, n_turb, d_cell2Face, d_posiCell2Face, d_face_number_of_each_cell, d_left_cell_of_face,
            d_right_cell_of_face, d_q_turb_proxy, coef, d_turboo, d_visl, turb_relax, nnegtive_max, ZoneID, SMALL, cv13,
            mytmax, my_nu_max, d_n_neg, d_n_pos, d_t_proxy, d_q_ns, d_dq_turb, d_q_turb_proxy);
#endif

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPURecoverResidual(const int nTotal, const int n_turb, const int a)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotal;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPURecoverResidual, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPURecoverResidual, 0, gridSize, blockSize);
#endif

        GPURecoverResidual<<<gridSize, blockSize>>>(nTotal, n_turb, d_res_turb, d_rhs_turb);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUTurbViscosity(const int nTotal, const int nTotalCell, const double mytmax, const double turb_cv13,
                              const double SMALL)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUNSSolverUnstruct;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotal;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUTurbViscosity, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUTurbViscosity, 0, gridSize, blockSize);
#endif
        /*
        gridSize= 1;
        blockSize =1;    //! it need to use reduction operations by sunxu, or use single thread
        GPUTurbViscosity<<<gridSize, blockSize>>>(nTotal, nTotalCell, mytmax, turb_cv13, SMALL, d_q_ns, d_q_turb_proxy, d_visl, d_vist, d_vistmax, d_vistmin);
        */
        GPUTurbViscosity_S1<<<gridSize, blockSize>>>(nTotal, nTotalCell, mytmax, turb_cv13, SMALL, d_q_ns,
                                                     d_q_turb_proxy, d_visl, d_vist);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        loopLen                = nTotalCell;
        int    regsPerThread   = 32;
        int    residentBlockNo = 4;
        size_t dsMemPerThread  = 2 * sizeof(RFloat); //! calculate by the intended shared memory use in the kernel
        size_t dsMemPerBlock;
        KernelLaunchPara((void *)GPUComputeMinMaxStep1, loopLen, regsPerThread, residentBlockNo, dsMemPerThread,
                         dsMemPerBlock, gridSize, blockSize, GPUProp);

        //! calculate the gridSize with register and shared memory limitation,
        //! because it will be used as blockSize in the second reduction.
        regsPerThread = 26; //! register for GPUComputeMinMaxStep2
        int gridSize_tmp;
        KernelLaunchPara((void *)GPUComputeMinMaxStep2, loopLen, regsPerThread, residentBlockNo, dsMemPerThread,
                         dsMemPerBlock, gridSize_tmp, gridSize, GPUProp);

        //!RFloat *bmin_tmp, *bmax_tmp;
        size_t bSize = gridSize * sizeof(RFloat);
        //!HANDLE_API_ERR(cudaMalloc((void**)&bmin_tmp, bSize));
        //!HANDLE_API_ERR(cudaMalloc((void**)&bmax_tmp, bSize));

//!KernelLaunchPara((void*)GPUReConstructFaceValueLoop2, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUComputeMinMaxStep1, dsMemPerBlock, gridSize, blockSize);
#endif

        //!GPUComputeMinMaxStep1<<<gridSize, blockSize, dsMemPerBlock>>>(d_vist, bmin_tmp, bmax_tmp, nTotalCell);
        GPUComputeMinMaxStep1<<<gridSize, blockSize, dsMemPerBlock>>>(d_vist, d_minLocalTurbVis, d_maxLocalTurbVis,
                                                                      nTotalCell);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        //!GPUComputeMinMaxStep2<<<1, gridSize, dsMemPerBlock>>>(bmin_tmp, bmax_tmp, d_vistmin, d_vistmax);
        GPUComputeMinMaxStep2<<<1, gridSize, dsMemPerBlock>>>(d_minLocalTurbVis, d_maxLocalTurbVis, d_vistmin,
                                                              d_vistmax);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        //!HANDLE_API_ERR(cudaFree(bmin_tmp));
        //!HANDLE_API_ERR(cudaFree(bmax_tmp));

        /*
        gridSize= 128;
                blockSize= 256;
                size_t dsMemPerBlock= 2*blockSize* sizeof(RFloat);

                RFloat *bmin_tmp, *bmax_tmp;
                size_t bSize= gridSize* sizeof(RFloat);
                HANDLE_API_ERR(cudaMalloc((void**)&bmin_tmp, bSize));
                HANDLE_API_ERR(cudaMalloc((void**)&bmax_tmp, bSize));

                GPUComputeMinMaxStep1<<<gridSize, blockSize, dsMemPerBlock>>>(d_vist, bmin_tmp, bmax_tmp, nTotalCell);
                #ifdef KERNELLAUNCHTEST
                HANDLE_KERNEL_ERR();
                #endif
                dsMemPerBlock= 2*gridSize* sizeof(RFloat);
                GPUComputeMinMaxStep2<<<1, gridSize, dsMemPerBlock>>>(bmin_tmp, bmax_tmp, d_vistmin, d_vistmax);
                #ifdef KERNELLAUNCHTEST
                HANDLE_KERNEL_ERR();
                #endif
                HANDLE_API_ERR(cudaFree(bmin_tmp));
                HANDLE_API_ERR(cudaFree(bmax_tmp));
        */
    }

    __global__ void GPUTurbZeroResiduals(const int numberOfGeneralCells, const int numberOfEquations, RFloat *res)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int iCell = 0;
        for (iCell = i; iCell < numberOfGeneralCells; iCell += blockDim.x * gridDim.x)
        {
            for (int iEquation = 0; iEquation < numberOfEquations; ++iEquation)
            {
                res[iEquation * numberOfGeneralCells + iCell] = 0.0;
            }
        }
    }

    __global__ void GPUTurbStoreRhsByResidual(const int nTotal, const int n_turb, const RFloat *res, RFloat *rhs)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int icell = 0;
        int m;
        for (icell = i; icell < nTotal; icell += blockDim.x * gridDim.x)
        {
            for (m = 0; m < n_turb; ++m)
            {
                rhs[m * nTotal + icell] = res[m * nTotal + icell];
            }
        }
    }

    __global__ void GPUTurbLoadQ(const int nTotalCell, const int nTotal, const int n_turb, const RFloat *q_turb,
                                 RFloat *qq_turb)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int icell = 0;
        for (icell = i; icell < nTotalCell; icell += blockDim.x * gridDim.x)
        {
            for (int m = 0; m < n_turb; ++m)
            {
                qq_turb[m * nTotal + icell] = q_turb[m * nTotal + icell];
            }
        }
    }

    __global__ void GPUTurbFillField(const int nTotal, const int neqn, RFloat *field1, const RFloat *field2)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int icell = 0;
        for (icell = i; icell < nTotal; icell += blockDim.x * gridDim.x)
        {
            for (int m = 0; m < neqn; ++m)
            {
                field1[m * nTotal + icell] = field2[m * nTotal + icell];
            }
        }
    }
    __global__ void GPUTurbLoadResiduals(const int nTotal, const int neqn, RFloat *res, const RFloat *rhs)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int icell = 0;
        for (icell = i; icell < nTotal; icell += blockDim.x * gridDim.x)
        {
            for (int m = 0; m < neqn; ++m)
            {
                res[m * nTotal + icell] = -rhs[m * nTotal + icell];
            }
        }
    }

    __global__ void GPUTurbLhs(const int nTotalCell, const int nTotal, const int n_turb, const double coef,
                               const RFloat *dt, RFloat *dq_turb)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int icell = 0;
        int m;
        for (icell = i; icell < nTotalCell; icell += blockDim.x * gridDim.x)
        {
            for (m = 0; m < n_turb; ++m)
            {
                dq_turb[m * nTotal + icell] = dt[icell] * coef * dq_turb[m * nTotal + icell];
            }
        }
    }

    __global__ void GPUTurbUpdateFlowField(const int nTotalCell, const int nTotal, const int n_turb,
                                           const int *cell2face, const int *posiCell2Face,
                                           const int *face_number_of_each_cell, const int *left_cell_of_face,
                                           const int *right_cell_of_face, const RFloat *q_turb, const double coef,
                                           const double *turboo, const RFloat *visl, const double turb_relax,
                                           const int nnegtive_max, const int ZoneID, const double SMALL,
                                           const double cv13, const double mytmax, const double my_nu_max, int *n_neg,
                                           int *n_pos, const RFloat *t, const RFloat *qpmv, const RFloat *dq, RFloat *q)
    {
        int     i     = blockDim.x * blockIdx.x + threadIdx.x;
        int     icell = 0;
        RFloat  qold, qnew;
        double  ld, fv1, ld3, nu_max;
        int     i_max;
        RFloat *prim_var = new RFloat[n_turb];
        for (icell = i; icell < nTotalCell; icell += blockDim.x * gridDim.x)
        {
            //!q[0][icell] = MAX( q[0][icell] + dq[0][icell], static_cast<RFloat>(zero) );
            qold = q[0 * nTotal + icell];
            qnew = q[0 * nTotal + icell] + coef * dq[0 * nTotal + icell];
            qnew = (1.0 - turb_relax) * qold + turb_relax * qnew;
            if (qnew < 0.0)
            {
                //!++ n_neg;
                atomicAdd(n_neg, 1);
                if ((*n_neg) <= nnegtive_max && true)
                {
                    printf("Zone = %d, cell = %d, nTotalCell = %d, q = %e, dq = %e\n", ZoneID, icell, nTotalCell,
                           q[0 * nTotal + icell], dq[0 * nTotal + icell]);
                }

                qnew = turboo[0];
            }
            q[0 * nTotal + icell] = qnew;

            ld     = GPUABS(q[0 * nTotal + icell]) / (visl[icell] / qpmv[0 * nTotal + icell] + SMALL);
            ld3    = ld * ld * ld;
            fv1    = ld3 / (ld3 + cv13);
            nu_max = visl[icell] * mytmax / (qpmv[0 * nTotal + icell] * fv1);
            if (q[0 * nTotal + icell] > nu_max)
            {
                //!SMoothTurbPoint(grid, prim_var, icell);
                GPUSMoothTurbPoint(nTotalCell, nTotal, n_turb, cell2face, posiCell2Face, left_cell_of_face,
                                   right_cell_of_face, face_number_of_each_cell, q_turb, icell, prim_var);
                q[0 * nTotal + icell] = prim_var[0];

                if (q[0 * nTotal + icell] > nu_max)
                {
                    q[0 * nTotal + icell] = nu_max;
                }
                if ((*n_pos) < nnegtive_max)
                {
                    printf("t[0][%d] = %e\n", icell, t[0 * nTotal + icell]);
                    printf("qpmv[0][%d] = %e\n", icell, qpmv[0 * nTotal + icell]);
                    printf("qpmv[1][%d] = %e\n", icell, qpmv[1 * nTotal + icell]);
                    printf("qpmv[2][%d] = %e\n", icell, qpmv[2 * nTotal + icell]);
                    printf("qpmv[3][%d] = %e\n", icell, qpmv[3 * nTotal + icell]);
                    printf("qpmv[4][%d] = %e\n", icell, qpmv[4 * nTotal + icell]);
                    printf("nu_max = %e\n", nu_max);
                    printf("visl[%d] = %e\n", icell, visl[icell]);
                    printf("fv1 = %e, ld = %e\n", fv1, ld);
                    printf("icell = %d, dq[0][%d] = ", icell, dq[0 * nTotal + icell]);
                }

                //!++ n_pos;
                atomicAdd(n_pos, 1);
            }
            //!some operation should be used like allreduce max
            /*
                    if ( my_nu_max < q[0 * nTotal + icell] ) {
                        my_nu_max = q[ 0 * nTotal + icell];
                        i_max = icell;
                    }
*/
        }

        delete[] prim_var;
    }
    __global__ void GPUTurbUpdateFlowField_S1(const int nTotalCell, const int nTotal, const int n_turb,
                                              const int *cell2face, const int *posiCell2Face,
                                              const int *face_number_of_each_cell, const int *left_cell_of_face,
                                              const int *right_cell_of_face, const RFloat *q_turb, const double coef,
                                              const double *turboo, const RFloat *visl, const double turb_relax,
                                              const int nnegtive_max, const int ZoneID, const double SMALL,
                                              const double cv13, const double mytmax, const double my_nu_max,
                                              int *n_neg, int *n_pos, const RFloat *t, const RFloat *qpmv,
                                              const RFloat *dq, RFloat *q)
    {
        int                      i     = blockDim.x * blockIdx.x + threadIdx.x;
        int                      icell = 0;
        RFloat                   qold, qnew;
        double                   ld, fv1, ld3, nu_max;
        int                      i_max;
        extern __shared__ RFloat array[];
        RFloat                  *prim_var = array;
        //!RFloat *prim_var = new RFloat [n_turb];

        for (icell = i; icell < nTotalCell; icell += blockDim.x * gridDim.x)
        {
            //!q[0][icell] = MAX( q[0][icell] + dq[0][icell], static_cast<RFloat>(zero) );
            qold = q[0 * nTotal + icell];
            qnew = q[0 * nTotal + icell] + coef * dq[0 * nTotal + icell];
            qnew = (1.0 - turb_relax) * qold + turb_relax * qnew;
            if (qnew < 0.0)
            {
                //!++ n_neg;
                atomicAdd(n_neg, 1);
                if ((*n_neg) <= nnegtive_max && true)
                {
                    printf("Zone = %d, cell = %d, nTotalCell = %d, q = %e, dq = %e\n", ZoneID, icell, nTotalCell,
                           q[0 * nTotal + icell], dq[0 * nTotal + icell]);
                }

                qnew = turboo[0];
            }
            q[0 * nTotal + icell] = qnew;

            ld     = GPUABS(q[0 * nTotal + icell]) / (visl[icell] / qpmv[0 * nTotal + icell] + SMALL);
            ld3    = ld * ld * ld;
            fv1    = ld3 / (ld3 + cv13);
            nu_max = visl[icell] * mytmax / (qpmv[0 * nTotal + icell] * fv1);
            if (q[0 * nTotal + icell] > nu_max)
            {
                //!SMoothTurbPoint(grid, prim_var, icell);
                GPUSMoothTurbPoint(nTotalCell, nTotal, n_turb, cell2face, posiCell2Face, left_cell_of_face,
                                   right_cell_of_face, face_number_of_each_cell, q_turb, icell,
                                   prim_var + n_turb * threadIdx.x);
                q[0 * nTotal + icell] = prim_var[0 + n_turb * threadIdx.x];

                if (q[0 * nTotal + icell] > nu_max)
                {
                    q[0 * nTotal + icell] = nu_max;
                }
                if ((*n_pos) < nnegtive_max)
                {
                    printf("t[0][%d] = %e\n", icell, t[0 * nTotal + icell]);
                    printf("qpmv[0][%d] = %e\n", icell, qpmv[0 * nTotal + icell]);
                    printf("qpmv[1][%d] = %e\n", icell, qpmv[1 * nTotal + icell]);
                    printf("qpmv[2][%d] = %e\n", icell, qpmv[2 * nTotal + icell]);
                    printf("qpmv[3][%d] = %e\n", icell, qpmv[3 * nTotal + icell]);
                    printf("qpmv[4][%d] = %e\n", icell, qpmv[4 * nTotal + icell]);
                    printf("nu_max = %e\n", nu_max);
                    printf("visl[%d] = %e\n", icell, visl[icell]);
                    printf("fv1 = %e, ld = %e\n", fv1, ld);
                    printf("icell = %d, dq[0][%d] = ", icell, dq[0 * nTotal + icell]);
                }

                //!++ n_pos;
                atomicAdd(n_pos, 1);
            }
            //!some operation should be used like allreduce max
            /*
                    if ( my_nu_max < q[0 * nTotal + icell] ) {
                        my_nu_max = q[ 0 * nTotal + icell];
                        i_max = icell;
                    }
*/
        }

        //!delete []prim_var;
    }

    __global__ void GPUSetNegPosZero(int *n_neg, int *n_pos)
    {
        (*n_neg) = 0;
        (*n_pos) = 0;
    }

    __device__ void GPUSMoothTurbPoint(const int nTotalCell, const int nTotal, const int n_turb, const int *cell2face,
                                       const int *posiCell2Face, const int *left_cell_of_face,
                                       const int *right_cell_of_face, const int *face_number_of_each_cell,
                                       const RFloat *q_turb, const int i, RFloat *prim_var)
    {
        int    m, j, face, le, re, ie;
        RFloat vol_sum = 0.0;
        for (m = 0; m < n_turb; ++m)
        {
            prim_var[m] = 0.0;
        }

        for (j = 0; j < face_number_of_each_cell[i]; ++j)
        {
            face = cell2face[posiCell2Face[i] + j];
            le   = left_cell_of_face[face];
            re   = right_cell_of_face[face];

            if (le != i)
            {
                ie = le;
            }
            else
            {
                ie = re;
            }

            if (ie >= nTotalCell) continue;

            for (m = 0; m < n_turb; ++m)
            {
                prim_var[m] += q_turb[m * nTotal + ie];
            }
            vol_sum += 1.0;
        }

        for (m = 0; m < n_turb; ++m)
        {
            prim_var[m] /= vol_sum;
        }
    }

    __global__ void GPURecoverResidual(const int nTotal, const int n_turb, RFloat *res, const RFloat *rhs)
    {
        int i = blockDim.x * blockIdx.x + threadIdx.x;
        int icell, m;
        for (icell = i; icell < nTotal; icell += blockDim.x * gridDim.x)
        {
            for (m = 0; m < n_turb; ++m)
            {
                res[m * nTotal + icell] = rhs[m * nTotal + icell];
            }
        }
    }

    __global__ void GPUTurbViscosity(const int nTotal, const int nTotalCell, const double mytmax,
                                     const double turb_cv13, const double SMALL, const RFloat *q, const RFloat *q_turb,
                                     const RFloat *visl, RFloat *vist, double *vistmax, double *vistmin)
    {
        int    i = blockDim.x * blockIdx.x + threadIdx.x;
        int    icell;
        double ld, ld3, fv1;
        vistmax[0] = 0.0;
        vistmin[0] = 1.0e30;

        for (icell = i; icell < nTotal; icell += blockDim.x * gridDim.x)
        {
            ld          = q_turb[0 * nTotal + icell] * q[0 * nTotal + icell] / (visl[icell] + SMALL);
            ld3         = ld * ld * ld;
            fv1         = ld3 / (ld3 + turb_cv13);
            vist[icell] = q[0 * nTotal + icell] * fv1 * q_turb[0 * nTotal + icell];
            vist[icell] = GPUMIN(mytmax, vist[icell]);          //!visl[icell] * coef_kvist *
            if (icell < nTotalCell && vistmax[0] < vist[icell]) //! bell 20130527 add icell < nTotalCell
            {
                vistmax[0] = vist[icell];
                //!imax = icell;
            }
            if (icell < nTotalCell && vistmin[0] > vist[icell]) //! bell 20130527 add icell < nTotalCell
            {
                vistmin[0] = vist[icell];
                //!imin = icell;
            }
        }
    }
    __global__ void GPUTurbViscosity_S1(const int nTotal, const int nTotalCell, const double mytmax,
                                        const double turb_cv13, const double SMALL, const RFloat *q,
                                        const RFloat *q_turb, const RFloat *visl, RFloat *vist)
    {
        int    i = blockDim.x * blockIdx.x + threadIdx.x;
        int    icell;
        double ld, ld3, fv1;

        for (icell = i; icell < nTotal; icell += blockDim.x * gridDim.x)
        {
            //!if (icell == 0) printf("In GPU, on term %d, q_turb = %.30e\n", icell, q_turb[0 * nTotal + icell]);
            ld          = q_turb[0 * nTotal + icell] * q[0 * nTotal + icell] / (visl[icell] + SMALL);
            ld3         = ld * ld * ld;
            fv1         = ld3 / (ld3 + turb_cv13);
            vist[icell] = q[0 * nTotal + icell] * fv1 * q_turb[0 * nTotal + icell];
            vist[icell] = GPUMIN(mytmax, vist[icell]); //!visl[icell] * coef_kvist *
                                                       /*
            if ( icell < nTotalCell && vistmax[0] < vist[icell] )          //! bell 20130527 add icell < nTotalCell
                        {
                            vistmax[0] =  vist[icell];
                            //!imax = icell;
                        }
                    if ( icell < nTotalCell && vistmin[0] > vist[icell] )          //! bell 20130527 add icell < nTotalCell
                        {
                            vistmin[0] =  vist[icell];
                            //!imin = icell;
                        }
            */
        }
    }
} //! namespace GPUKernels
