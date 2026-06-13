#include "GPULoadFlux.h"
#include "cudaErrorHandle.h"
#include "GPUDeviceControl.h"
#if __CUDA_ARCH__ < 600
#include "GPUPortability.h"
#endif

using namespace std;
namespace GPUKernels
{
    void CallGPULoadFlux(Grid *grid_in, int nst, int ned)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;
        using namespace GPUControlVariables;

        UnstructGrid *grid = UnstructGridCast(grid_in);

        int nBoundFace = grid->GetNBoundFace();
        int nTotalCell = grid->GetNTotalCell();
        int nTotal     = nTotalCell + nBoundFace;

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

        /* original kernel size
                gridSize= 2092;
                blockSize= 640;
                */
        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nMid - nst;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPULoadFluxOnBound, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPULoadFluxOnBound, 0, gridSize, blockSize);
#endif

        //!if (d_visflux_solver != "NSSolverUnstruct") return;
        GPULoadFluxOnBound<<<gridSize, blockSize>>>(nst, nMid, nTotalCell, nBoundFace, nl, d_SEG_LEN,
                                                    d_left_cell_of_face, d_flux, d_res_ns);
//! #ifndef CUDADEBUGATOMIC
//!     GPULoadFluxOnBound<<<gridSize, blockSize>>>(nst, nMid, nTotalCell, nBoundFace, nl, d_SEG_LEN, d_left_cell_of_face, d_flux, d_res_ns);
//! #else
//!     GPULoadFluxOnBound<<<1, 1>>>(nst, nMid, nTotalCell, nBoundFace, nl, d_SEG_LEN, d_left_cell_of_face, d_flux, d_res_ns);
//! #endif
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        loopLen = ned - nMid;
#ifndef NESTLOOPOPT
        KernelLaunchPara((void *)GPULoadFluxInterior, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPULoadFluxInterior, 0, gridSize, blockSize);
#endif
#else
        KernelLaunchPara((void *)GPULoadFluxInteriorOptSep, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPULoadFluxInteriorOptSep, 0, gridSize, blockSize);
#endif
#endif

#ifndef CUDADEBUGATOMIC
#ifndef NESTLOOPOPT
        GPULoadFluxInterior<<<gridSize, blockSize>>>(nst, nMid, ned, nTotalCell, nBoundFace, nl, d_SEG_LEN,
                                                     d_left_cell_of_face, d_right_cell_of_face, d_flux, d_res_ns);
#else
        for (int equationID = 0; equationID < nl; equationID++)
        {
            GPULoadFluxInteriorOptSep<<<gridSize, blockSize>>>(equationID, nst, nMid, ned, nTotalCell, nBoundFace, nl,
                                                               d_SEG_LEN, d_left_cell_of_face, d_right_cell_of_face,
                                                               d_flux, d_res_ns);
#ifdef KERNELLAUNCHTEST
            HANDLE_KERNEL_ERR();
#endif
        }
#endif
#else
#ifndef NESTLOOPOPT
        GPULoadFluxInterior<<<1, 1>>>(nst, nMid, ned, nTotalCell, nBoundFace, nl, d_SEG_LEN, d_left_cell_of_face,
                                      d_right_cell_of_face, d_flux, d_res_ns);
#else
        for (int equationID = 0; equationID < nl; equationID++)
        {
            GPULoadFluxInteriorOptSep<<<1, 1>>>(equationID, nst, nMid, ned, nTotalCell, nBoundFace, nl, d_SEG_LEN,
                                                d_left_cell_of_face, d_right_cell_of_face, d_flux, d_res_ns);
#ifdef KERNELLAUNCHTEST
            HANDLE_KERNEL_ERR();
#endif
        }

#endif
#endif
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPULoadFluxOnBound(const int nst, const int nMid, const int nTotalCell, const int nBoundFace,
                                       const int nl, const int len, const int *left_cell_of_face, const RFloat *flux,
                                       RFloat *res)
    {
        int k = blockIdx.x * blockDim.x + threadIdx.x;
        int i = 0;
        int le, j;
        int nTotal = nTotalCell + nBoundFace;
        for (i = k + nst; i < nMid; i += blockDim.x * gridDim.x)
        {
            le = left_cell_of_face[i];
            j  = i - nst;
            for (int m = 0; m < nl; ++m)
            {
//!res[m][le] -= flux[m][j]
#if __CUDA_ARCH__ < 600
                atomicAddTest(res + m * nTotal + le, -flux[m * len + j]);
#else
                atomicAdd(res + m * nTotal + le, -flux[m * len + j]);
#endif
            }
        }
    }

    __global__ void GPULoadFluxInterior(const int nst, const int nMid, const int ned, const int nTotalCell,
                                        const int nBoundFace, const int nl, const int len, const int *left_cell_of_face,
                                        const int *right_cell_of_face, const RFloat *flux, RFloat *res)
    {
        //! Interior faces
        int k      = blockDim.x * blockIdx.x + threadIdx.x;
        int nTotal = nTotalCell + nBoundFace;
        int le, re, j;
        for (int i = nMid + k; i < ned; i += blockDim.x * gridDim.x)
        {
            le = left_cell_of_face[i];
            re = right_cell_of_face[i];
            j  = i - nst;

            for (int m = 0; m < nl; ++m)
            {
                //!res[m][le] -= flux[m][j];
                //!res[m][re] += flux[m][j];
#if __CUDA_ARCH__ < 600
                atomicAddTest(res + m * nTotal + le, -flux[m * len + j]);
                atomicAddTest(res + m * nTotal + re, flux[m * len + j]);
#else
                atomicAdd(res + m * nTotal + le, -flux[m * len + j]);
                atomicAdd(res + m * nTotal + re, flux[m * len + j]);
#endif
            }
        }
    }

    __global__ void GPULoadFluxInteriorOptSep(const int equationID, const int nst, const int nMid, const int ned,
                                              const int nTotalCell, const int nBoundFace, const int nl, const int len,
                                              const int *left_cell_of_face, const int *right_cell_of_face,
                                              const RFloat *flux, RFloat *res)
    {
        //! Interior faces
        int k      = blockDim.x * blockIdx.x + threadIdx.x;
        int nTotal = nTotalCell + nBoundFace;
        int le, re, j;
        for (int i = nMid + k; i < ned; i += blockDim.x * gridDim.x)
        {
            le = left_cell_of_face[i];
            re = right_cell_of_face[i];
            j  = i - nst;

#if __CUDA_ARCH__ < 600
            atomicAddTest(res + equationID * nTotal + le, -flux[equationID * len + j]);
            atomicAddTest(res + equationID * nTotal + re, flux[equationID * len + j]);
#else
            atomicAdd(res + equationID * nTotal + le, -flux[equationID * len + j]);
            atomicAdd(res + equationID * nTotal + re, flux[equationID * len + j]);
#endif
        }
    }
} //! namespace GPUKernels
