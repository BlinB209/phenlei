#include "GPUNSUnsteady.h"
#include "cudaErrorHandle.h"
#include "GPUDeviceControl.h"
#include "GPULUSGS.h" //!for Primitive2Conservative
#if __CUDA_ARCH__ < 600
#include "GPUPortability.h"
#endif

namespace GPUKernels
{
    void CallGPUNSUnsteadyInitSpec(const int nTotal, const RFloat dualTimeSpectrumC1, const RFloat dualTimeSpectrumC2)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotal;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUNSUnsteadyInitSpec, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUNSUnsteadyInitSpec, 0, gridSize, blockSize);
#endif
        GPUNSUnsteadyInitSpec<<<gridSize, blockSize>>>(nTotal, dualTimeSpectrumC1, dualTimeSpectrumC2, d_vol, d_dt,
                                                       d_spec_ns);

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUNSUnsteadyInitSpec(const int nTotal, const RFloat dualTimeSpectrumC1,
                                          const RFloat dualTimeSpectrumC2, const double *vol, const RFloat *dt,
                                          RFloat *spec)
    {
        int threadID = blockDim.x * blockIdx.x + threadIdx.x;

        for (int iCell = threadID; iCell < nTotal; iCell += blockDim.x * gridDim.x)
        {
            spec[iCell] = dualTimeSpectrumC1 * vol[iCell] + (dualTimeSpectrumC2 / dt[iCell]);
            //! notice: the dt[iCell] has divided vol, as same as above.
        }
    }

    void CallGPUNSUpdateUnsteadyFlow(const int nTotal, const int nEquation, const int ifStaticsFlowField)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        if (ifStaticsFlowField != 0)
        {
            printf("Error in CallNSUpdateUnsteadyFlow, ifStaticsFlowField != 0\n");
            exit(1);
        }

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotal;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUNSUpdateUnsteadyFlow, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUNSUpdateUnsteadyFlow, 0, gridSize, blockSize);
#endif
        GPUNSUpdateUnsteadyFlow<<<gridSize, blockSize>>>(nTotal, nEquation, d_q_ns, d_res_ns_unsteady_tmp,
                                                         d_q_ns_unsteady_n1, d_q_ns_unsteady_n2, d_res_ns_unsteady_n1,
                                                         d_res_ns_unsteady_n2);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUNSUpdateUnsteadyFlow(const int nTotal, const int nEquation, const RFloat *q,
                                            const RFloat *resTmp, RFloat *qn1, RFloat *qn2, RFloat *resn1,
                                            RFloat *resn2)
    {
        int threadID = blockDim.x * blockIdx.x + threadIdx.x;
        for (int iCell = threadID; iCell < nTotal; iCell += blockDim.x * gridDim.x)
        {
            for (int m = 0; m < nEquation; ++m)
            {
                qn2[m * nTotal + iCell]   = qn1[m * nTotal + iCell];
                qn1[m * nTotal + iCell]   = q[m * nTotal + iCell];
                resn2[m * nTotal + iCell] = resn1[m * nTotal + iCell];
                resn1[m * nTotal + iCell] = resTmp
                    [m * nTotal
                     + iCell]; //! Here the current outnstep is over, the value of the stored resTmp should be assigned to resn1 for the next outnstep. It should be noticed that resTmp only contain the inviscid and viscous flux.
            }
        }
    }

    void CallGPUNSUnsteadyConvergence(const int nTotal, const int nTotalCell, const int nEquation, const int nchem,
                                      const int nm, const int nl, const int ntmodel, const double refGama)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        setSumAsZero<<<1, 1>>>(d_sum1_ns_unsteady, d_sum2_ns_unsteady);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUNSUnsteadyConvergence, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUNSUnsteadyConvergence, 0, gridSize, blockSize);
#endif

#ifndef CUDADEBUGATOMIC
        GPUNSUnsteadyConvergence<<<gridSize, blockSize>>>(nTotal, nTotalCell, nEquation, nchem, nm, nl, ntmodel,
                                                          refGama, d_q_ns, d_q_ns_unsteady_n1, d_q_ns_unsteady_n2,
                                                          d_res_ns, d_sum1_ns_unsteady, d_sum2_ns_unsteady);
#else
        GPUNSUnsteadyConvergence<<<1, 1>>>(nTotal, nTotalCell, nEquation, nchem, nm, nl, ntmodel, refGama, d_q_ns,
                                           d_q_ns_unsteady_n1, d_q_ns_unsteady_n2, d_res_ns, d_sum1_ns_unsteady,
                                           d_sum2_ns_unsteady);
#endif
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void setSumAsZero(RFloat *sum1, RFloat *sum2)
    {
        *sum1 = 0.0;
        *sum2 = 0.0;
    }

    __global__ void GPUNSUnsteadyConvergence(const int nTotal, const int nTotalCell, const int nEquation,
                                             const int nchem, const int nm, const int nl, const int ntmodel,
                                             const double refGama, const RFloat *q, const RFloat *qn1,
                                             const RFloat *qn2, const RFloat *res, RFloat *sum1, RFloat *sum2)
    {
        int    threadID = blockDim.x * blockIdx.x + threadIdx.x;
        RFloat prim0[5], prim1[5], prim2[5];
        RFloat qcsv0[5], qcsv1[5], qcsv2[5];
        RFloat Tv[3], Te[3];
        Tv[0] = 0.0;
        Tv[1] = 0.0;
        Tv[2] = 0.0;
        Te[0] = 0.0;
        Te[1] = 0.0;
        Te[2] = 0.0;
        for (int iCell = threadID; iCell < nTotalCell; iCell += blockDim.x * gridDim.x)
        {
            for (int m = 0; m < nEquation; ++m)
            {
                prim0[m] = q[m * nTotal + iCell];
                prim1[m] = qn1[m * nTotal + iCell];
                prim2[m] = qn2[m * nTotal + iCell];
            }

            if (ntmodel > 1)
            {
                printf("Error: in GPUNSUnsteadyConvergence, ntmodel>1\n");
            }

            //!Primitive2ConservativeÊµ¼ÊÉÏ²»ÐèÒªÊäÈëgama£¬ºóÃæÐèÒª¸ÄÐ´Õâ¸öÎÊÌâ
            Primitive2Conservative(nchem, nm, nl, ntmodel, prim0, refGama, Tv[0], Te[0], qcsv0);
            Primitive2Conservative(nchem, nm, nl, ntmodel, prim1, refGama, Tv[1], Te[1], qcsv1);
            Primitive2Conservative(nchem, nm, nl, ntmodel, prim2, refGama, Tv[2], Te[2], qcsv2);

            //!ÕâÀïÐèÒªnl
            for (int m = 0; m < nl; ++m)
            {
                //!´ËÊ±res±ØÐëÒÑ¾­×ª»¯Îªdq£¬¼´´ËÊ±µÄ²Ð²î²»ÊÇÓÒ¶ËÏî£¬¶øÊÇdq
                RFloat dq_p = res[m * nTotal + iCell]; //! qn+1,p+1 - qn+1,p
                RFloat dq_n = qcsv0[m] - qcsv1[m];     //! qn+1,p+1 - qn
//!data race, many threads write into the same sum1 and sum2
#if __CUDA_ARCH__ < 600
                atomicAddTest(sum1, dq_p * dq_p);
                atomicAddTest(sum2, dq_n * dq_n);
#else
                atomicAdd(sum1, dq_p * dq_p);
                atomicAdd(sum2, dq_n * dq_n);
#endif
            }
        }
    }

    void CallGPUNSDualTimeSourceRes(const int nTotal, const int nTotalCell, const int nEquation, const int nchem,
                                    const int nm, const int nl, const int ntmodel, const RDouble refGama,
                                    const RFloat dualTimeResC1, const RFloat dualTimeResC2, const RFloat dualTimeResC3,
                                    const RFloat dualTimeQC1, const RFloat dualTimeQC2, const RFloat dualTimeQC3)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUNSDualTimeSourceRes, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUNSDualTimeSourceRes, 0, gridSize, blockSize);
#endif
        GPUNSDualTimeSourceRes<<<gridSize, blockSize>>>(
            nTotal, nTotalCell, nEquation, nchem, nm, nl, ntmodel, refGama, dualTimeResC1, dualTimeResC2, dualTimeResC3,
            dualTimeQC1, dualTimeQC2, dualTimeQC3, d_vol, d_q_ns, d_q_ns_unsteady_n1, d_q_ns_unsteady_n2,
            d_res_ns_unsteady_n1, d_res_ns_unsteady_n2, d_res_ns);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUNSDualTimeSourceRes(const int nTotal, const int nTotalCell, const int nEquation, const int nchem,
                                           const int nm, const int nl, const int ntmodel, const RDouble refGama,
                                           const RFloat dualTimeResC1, const RFloat dualTimeResC2,
                                           const RFloat dualTimeResC3, const RFloat dualTimeQC1,
                                           const RFloat dualTimeQC2, const RFloat dualTimeQC3, const RDouble *vol,
                                           const RFloat *q, const RFloat *qn1, const RFloat *qn2, const RFloat *resn1,
                                           const RFloat *resn2, RFloat *res)
    {
        int    threadID = blockIdx.x * blockDim.x + threadIdx.x;
        RFloat prim0[5], prim1[5], prim2[5];
        RFloat qcsv0[5], qcsv1[5], qcsv2[5];
        RFloat Tv[3], Te[3];
        Tv[0] = 0.0;
        Tv[1] = 0.0;
        Tv[2] = 0.0;
        Te[0] = 0.0;
        Te[1] = 0.0;
        Te[2] = 0.0;

        for (int iCell = threadID; iCell < nTotalCell; iCell += blockDim.x * gridDim.x)
        {
            for (int m = 0; m < nEquation; ++m)
            {
                prim0[m] = q[m * nTotal + iCell];
                prim1[m] = qn1[m * nTotal + iCell];
                prim2[m] = qn2[m * nTotal + iCell];
            }
            if (ntmodel > 1)
            {
                printf("Error in GPUNSDualTimeSourceRes, ntmodel>1\n");
            }

            //! In fact, there is no need to input the variable of refGama for the function of Primitive2Conservative, it is to be improved.
            Primitive2Conservative(nchem, nm, nl, ntmodel, prim0, refGama, Tv[0], Te[0], qcsv0);
            Primitive2Conservative(nchem, nm, nl, ntmodel, prim1, refGama, Tv[1], Te[1], qcsv1);
            Primitive2Conservative(nchem, nm, nl, ntmodel, prim2, refGama, Tv[2], Te[2], qcsv2);

            for (int m = 0; m < nEquation; ++m)
            {
                RFloat dualSrcRes = dualTimeResC1 * res[m * nTotal + iCell] + dualTimeResC2 * resn1[m * nTotal + iCell]
                                    + dualTimeResC3 * resn2[m * nTotal + iCell];

                RFloat dualSrcQ = dualTimeQC1 * qcsv0[m] * vol[iCell] + dualTimeQC2 * qcsv1[m] * vol[iCell]
                                  + dualTimeQC3 * qcsv2[m] * vol[iCell];

                res[m * nTotal + iCell] = dualSrcRes + dualSrcQ;
            }
        }
    }

    void CallGPUSetNSResTmpbyRes(const int nTotal, const int nTotalCell, const int nEquation)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUSetNSResTmpbyRes, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUSetNSResTmpbyRes, 0, gridSize, blockSize);
#endif
        GPUSetNSResTmpbyRes<<<gridSize, blockSize>>>(nTotal, nTotalCell, nEquation, d_res_ns, d_res_ns_unsteady_tmp);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUSetNSResTmpbyRes(const int nTotal, const int nTotalCell, const int nEquation, const RFloat *res,
                                        RFloat *resTmp)
    {
        int threadID = blockDim.x * blockIdx.x + threadIdx.x;
        for (int cellID = threadID; cellID < nTotalCell; cellID += blockDim.x * gridDim.x)
        {
            for (int m = 0; m < nEquation; m++)
            {
                resTmp[m * nTotal + cellID] = res[m * nTotal + cellID];
            }
        }
    }

} //! namespace GPUKernels
