#include "GPUCompGradientGGCell.h"
#include "cudaErrorHandle.h"
#include "GPUDeviceControl.h"
#if __CUDA_ARCH__ < 600
#include "GPUPortability.h"
#endif

using namespace std;
namespace GPUKernels
{
    void CallGPUCompGradientGGCell(const string d_gradient_field_proxy, const int d_index, const int nTotalCell,
                                   const int nBoundFace, const int nTotalFace)
    {
        RFloat *dqdx;
        RFloat *dqdy;
        RFloat *dqdz;
        RFloat *q;
        int     index;
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        if (d_gradient_field_proxy == "d_q_proxy")
        {
            dqdx  = d_dqdx_proxy;
            dqdy  = d_dqdy_proxy;
            dqdz  = d_dqdz_proxy;
            q     = d_q_ns;
            index = d_index;
        }
        else if (d_gradient_field_proxy == "d_t_proxy")
        {
            dqdx  = d_dtdx_proxy;
            dqdy  = d_dtdy_proxy;
            dqdz  = d_dtdz_proxy;
            q     = d_t_proxy;
            index = d_index;
        }
        else if (d_gradient_field_proxy == "d_q_turb_proxy")
        {
            dqdx  = d_dq_turbdx_proxy;
            dqdy  = d_dq_turbdy_proxy;
            dqdz  = d_dq_turbdz_proxy;
            q     = d_q_turb_proxy;
            index = d_index;
        }
        else if (d_gradient_field_proxy == "d_velocity_proxy")
        {
            dqdx  = d_dveldx_proxy;
            dqdy  = d_dveldy_proxy;
            dqdz  = d_dveldz_proxy;
            q     = d_vel_proxy;
            index = d_index;
        }
        else
        {
            cout << "Unknown gradient_field_proxy for kernel GPUCompGradientGGCell" << endl;
            exit(1);
            return;
        }

        int nTotal         = nTotalCell + nBoundFace;
        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotal;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUCompGradientGGCellInit, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompGradientGGCellInit, 0, gridSize, blockSize);
#endif

        //！set zero into dqdx, dqdy, dqdz
        GPUCompGradientGGCellInit<<<gridSize, blockSize>>>(index, nTotal, dqdx, dqdy, dqdz);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        loopLen = nBoundFace;
        KernelLaunchPara((void *)GPUCompGradientGGCellBoundaryFaceCal, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompGradientGGCellBoundaryFaceCal, 0, gridSize, blockSize);
#endif

        GPUCompGradientGGCellBoundaryFaceCal<<<gridSize, blockSize>>>(
            index, nTotalCell, nBoundFace, d_left_cell_of_face, d_xfn, d_yfn, d_zfn, d_area, q, dqdx, dqdy, dqdz);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        loopLen = nTotalFace - nBoundFace;
        KernelLaunchPara((void *)GPUCompGradientGGCellInteriorFaceCal, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompGradientGGCellInteriorFaceCal, 0, gridSize, blockSize);
#endif

        GPUCompGradientGGCellInteriorFaceCal<<<gridSize, blockSize>>>(index, nTotalCell, nTotalFace, nBoundFace,
                                                                      d_left_cell_of_face, d_right_cell_of_face, d_xfn,
                                                                      d_yfn, d_zfn, d_area, q, dqdx, dqdy, dqdz);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        loopLen = nTotalCell;
        KernelLaunchPara((void *)GPUCompGradientGGCellDiVol, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompGradientGGCellDiVol, 0, gridSize, blockSize);
#endif

        GPUCompGradientGGCellDiVol<<<gridSize, blockSize>>>(index, nTotalCell, nTotal, dqdx, dqdy, dqdz, d_vol);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
        loopLen = nBoundFace;
        KernelLaunchPara((void *)GPUCompGradientGGCellGhostSet, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompGradientGGCellGhostSet, 0, gridSize, blockSize);
#endif

        GPUCompGradientGGCellGhostSet<<<gridSize, blockSize>>>(index, nTotalCell, nTotal, nBoundFace,
                                                               d_left_cell_of_face, dqdx, dqdy, dqdz);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUCompGradientGGCellInit(const int index, const int nTotal, RFloat *dqdx, RFloat *dqdy,
                                              RFloat *dqdz)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int icell = 0;
        for (icell = i; icell < nTotal; icell += blockDim.x * gridDim.x)
        {
            dqdx[icell + index * nTotal] = 0.0;
            dqdy[icell + index * nTotal] = 0.0;
            dqdz[icell + index * nTotal] = 0.0;
        }
    }

    __global__ void GPUCompGradientGGCellBoundaryFaceCal(const int index, const int nTotalCell, const int nBoundFace,
                                                         const int *left_cell_of_face, const RDouble *nxs,
                                                         const RDouble *nys, const RDouble *nzs, const RDouble *ns,
                                                         RFloat *q, RFloat *dqdx, RFloat *dqdy, RFloat *dqdz)
    {
        const int tidx   = threadIdx.x;
        const int bidx   = blockIdx.x;
        int       nTotal = nTotalCell + nBoundFace;
        for (int iface = bidx * blockDim.x + tidx; iface < nBoundFace; iface += blockDim.x * gridDim.x)
        {
            int le = left_cell_of_face[iface];
            int re = iface + nTotalCell;

            RFloat  qfc = half * (q[le + index * nTotal] + q[re + index * nTotal]);
            RDouble nx  = nxs[iface] * ns[iface];
            RDouble ny  = nys[iface] * ns[iface];
            RDouble nz  = nzs[iface] * ns[iface];

#if __CUDA_ARCH__ < 600
            atomicAddTest(dqdx + (nTotalCell + nBoundFace) * index + le, qfc * nx);
            atomicAddTest(dqdy + (nTotalCell + nBoundFace) * index + le, qfc * ny);
            atomicAddTest(dqdz + (nTotalCell + nBoundFace) * index + le, qfc * nz);
#else
            atomicAdd(dqdx + (nTotalCell + nBoundFace) * index + le, qfc * nx);
            atomicAdd(dqdy + (nTotalCell + nBoundFace) * index + le, qfc * ny);
            atomicAdd(dqdz + (nTotalCell + nBoundFace) * index + le, qfc * nz);
#endif
        }
    }

    __global__ void GPUCompGradientGGCellInteriorFaceCal(const int index, const int nTotalCell, const int nTotalFace,
                                                         const int nBoundFace, const int *left_cell_of_face,
                                                         const int *right_cell_of_face, const RDouble *nxs,
                                                         const RDouble *nys, const RDouble *nzs, const RDouble *ns,
                                                         RFloat *q, RFloat *dqdx, RFloat *dqdy, RFloat *dqdz)
    {
        const int tidx   = threadIdx.x;
        const int bidx   = blockIdx.x;
        int       nTotal = nTotalCell + nBoundFace;
        for (int iface = nBoundFace + bidx * blockDim.x + tidx; iface < nTotalFace; iface += blockDim.x * gridDim.x)
        {
            int le = left_cell_of_face[iface];
            int re = right_cell_of_face[iface];

            RFloat  qfc = half * (q[le + index * nTotal] + q[re + index * nTotal]);
            RDouble nx  = nxs[iface] * ns[iface];
            RDouble ny  = nys[iface] * ns[iface];
            RDouble nz  = nzs[iface] * ns[iface];
/*
            dqdx[le] += qfc * nx;
            dqdy[le] += qfc * ny;
            dqdz[le] += qfc * nz;
            dqdx[re] -= qfc * nx;
            dqdy[re] -= qfc * ny;
            dqdz[re] -= qfc * nz;
            */
#if __CUDA_ARCH__ < 600
            atomicAddTest(dqdx + (nTotalCell + nBoundFace) * index + le, qfc * nx);
            atomicAddTest(dqdy + (nTotalCell + nBoundFace) * index + le, qfc * ny);
            atomicAddTest(dqdz + (nTotalCell + nBoundFace) * index + le, qfc * nz);

            atomicAddTest(dqdx + (nTotalCell + nBoundFace) * index + re, -qfc * nx);
            atomicAddTest(dqdy + (nTotalCell + nBoundFace) * index + re, -qfc * ny);
            atomicAddTest(dqdz + (nTotalCell + nBoundFace) * index + re, -qfc * nz);
#else
            atomicAdd(dqdx + (nTotalCell + nBoundFace) * index + le, qfc * nx);
            atomicAdd(dqdy + (nTotalCell + nBoundFace) * index + le, qfc * ny);
            atomicAdd(dqdz + (nTotalCell + nBoundFace) * index + le, qfc * nz);

            atomicAdd(dqdx + (nTotalCell + nBoundFace) * index + re, -qfc * nx);
            atomicAdd(dqdy + (nTotalCell + nBoundFace) * index + re, -qfc * ny);
            atomicAdd(dqdz + (nTotalCell + nBoundFace) * index + re, -qfc * nz);
#endif
        }
    }

    __global__ void GPUCompGradientGGCellDiVol(const int index, const int nTotalCell, const int nTotal, RFloat *dqdx,
                                               RFloat *dqdy, RFloat *dqdz, RDouble *vol)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int icell = 0;
        for (icell = i; icell < nTotalCell; icell += blockDim.x * gridDim.x)
        {
            //！just for test
            //！if ((icell == 0)||(icell == 1)) printf("vol[%d]=%e\n", icell, vol[icell]);
            RDouble ovol = 1.0 / vol[icell];
            dqdx[icell + index * nTotal] *= ovol;
            dqdy[icell + index * nTotal] *= ovol;
            dqdz[icell + index * nTotal] *= ovol;
        }
    }

    __global__ void GPUCompGradientGGCellGhostSet(const int index, const int nTotalCell, const int nTotal,
                                                  const int nBoundFace, const int *left_cell_of_face, RFloat *dqdx,
                                                  RFloat *dqdy, RFloat *dqdz)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int iface = 0;
        for (iface = i; iface < nBoundFace; iface += blockDim.x * gridDim.x)
        {
            int le = left_cell_of_face[iface];
            int re = iface + nTotalCell;

            dqdx[re + index * nTotal] = dqdx[le + index * nTotal];
            dqdy[re + index * nTotal] = dqdy[le + index * nTotal];
            dqdz[re + index * nTotal] = dqdz[le + index * nTotal];
        }
    }
} //！ namespace GPUKernels
