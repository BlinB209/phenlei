#include "GPUStoreBoundGrad.h"
#include "cudaErrorHandle.h"
#include "GPUDeviceControl.h"

using namespace std;
namespace GPUKernels
{
    void CallGPUStoreBoundGrad(int nTotalCell, int nBoundFace, int nTVar)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nBoundFace;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUStoreBoundGrad, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUStoreBoundGrad, 0, gridSize, blockSize);
#endif

        GPUStoreBoundGrad<<<gridSize, blockSize>>>(nTotalCell, nBoundFace, nTVar, d_left_cell_of_face, d_dqdx_proxy,
                                                   d_dqdy_proxy, d_dqdz_proxy, d_bdqx, d_bdqy, d_bdqz);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUStoreBoundGrad(const int nTotalCell, const int nBoundFace, const int nTVar,
                                      const int *left_cell_of_face, RFloat *d_dqdx_proxy, RFloat *d_dqdy_proxy,
                                      RFloat *d_dqdz_proxy, RFloat *d_bdqx, RFloat *d_bdqy, RFloat *d_bdqz)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int iface = 0;
        int le, j;
        for (iface = i; iface < nBoundFace; iface += blockDim.x * gridDim.x)
        {
            le = left_cell_of_face[iface];
            for (j = 0; j < nTVar; j++)
            {
                d_bdqx[j * nBoundFace + iface] = d_dqdx_proxy[j * (nTotalCell + nBoundFace) + le];
                d_bdqy[j * nBoundFace + iface] = d_dqdy_proxy[j * (nTotalCell + nBoundFace) + le];
                d_bdqz[j * nBoundFace + iface] = d_dqdz_proxy[j * (nTotalCell + nBoundFace) + le];
            }
        }
    }
} //! namespace GPUKernels
