#include "GPUCompViscousCoef.h"
#include "cudaErrorHandle.h"
#include "GPUDeviceControl.h"

using namespace std;
namespace GPUKernels
{
    //!I think it is no use!!
    void CallGPUCompViscousCoef(const int nTotal, const int ITT, const double tsuth, const double visl_min)
    {
        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        //!GPUCompViscousCoef<<<2092 ,640>>>(nTotal, ITT, tsuth, visl_min, d_visl, d_t_proxy);
        int gridSize  = 1;
        int blockSize = 1;
        int loopLen   = nTotal;
        KernelLaunchPara((void *)GPUCompViscousCoef, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUCompViscousCoef, 0, gridSize, blockSize);
#endif

        GPUCompViscousCoef<<<gridSize, blockSize>>>(nTotal, ITT, tsuth, visl_min, d_t_proxy);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUCompViscousCoef(const int nTotal, const int ITT, const double tsuth, const double visl_min,
                                       RFloat *t)
    {
        //!RFloat * visl, RFloat * t){
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int icell = 0;
        //!t[icell] = 0;
        for (icell = i; icell < nTotal; icell += blockDim.x * gridDim.x)
        {
            //!visl[icell] = t[ITT*nTotal+icell]*sqrt(t[ITT*nTotal+icell]);
            t[icell] = 0.0;
        }
    }
} //! namespace GPUKernels
