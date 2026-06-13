#include "GPUSpectrumRadius.h"
#include "cudaErrorHandle.h"
#include "GPUDeviceControl.h"
using namespace std;
namespace GPUKernels
{
    void CallGPUInitSpec(const int nTotal)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotal;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUInitSpec, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUInitSpec, 0, gridSize, blockSize);
#endif

        GPUInitSpec<<<gridSize, blockSize>>>(nTotal, d_spec_ns, d_dt);

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUAddInvSpec(const int nTotalCell)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUAddInvSpec, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUAddInvSpec, 0, gridSize, blockSize);
#endif

        GPUAddInvSpec<<<gridSize, blockSize>>>(nTotalCell, d_spec_ns, d_invSpectrumRadius);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    void CallGPUAddVisSpec(const int nTotalCell)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUGeomVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nTotalCell;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUAddVisSpec, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUAddVisSpec, 0, gridSize, blockSize);
#endif

        GPUAddVisSpec<<<gridSize, blockSize>>>(nTotalCell, d_spec_ns, d_visSpectrumRadius);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }
    __global__ void GPUInitSpec(const int nTotal, RFloat *spec, const RFloat *dt)
    {
        int threadID = blockIdx.x * blockDim.x + threadIdx.x;
        for (int cellID = threadID; cellID < nTotal; cellID += gridDim.x * blockDim.x)
        {
            spec[cellID] = 1.0 / dt[cellID];
        }
    }

    __global__ void GPUAddInvSpec(const int nTotalCell, RFloat *spec, const RFloat *invSpectrumRadius)
    {
        int threadID = blockIdx.x * blockDim.x + threadIdx.x;
        for (int cellID = threadID; cellID < nTotalCell; cellID += gridDim.x * blockDim.x)
        {
            spec[cellID] += invSpectrumRadius[cellID];
        }
    }

    __global__ void GPUAddVisSpec(const int nTotalCell, RFloat *spec, const RFloat *visSpectrumRadius)
    {
        int threadID = blockIdx.x * blockDim.x + threadIdx.x;
        for (int cellID = threadID; cellID < nTotalCell; cellID += gridDim.x * blockDim.x)
        {
            spec[cellID] += visSpectrumRadius[cellID];
        }
    }
} //! namespace GPUKernels
