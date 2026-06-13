#include "GPUFaceWeight.h"
#include "cudaErrorHandle.h"
#include "GPUDeviceControl.h"
#include "Geo_SimpleBC.h"
using namespace std;
using namespace PHENGLEI;
namespace GPUKernels
{

    void CallGPUFaceWeight(const int nst, const int ned, const int nBoundFace)
    {
        using namespace GPUMemory;
        using namespace GPUControlVariables;
        using namespace GPUGeomVariables;
        using namespace GPUFlowVariables;

        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = ned - nst;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUFaceWeight, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUFaceWeight, 0, gridSize, blockSize);
#endif

        GPUFaceWeight<<<gridSize, blockSize>>>(nst, ned, d_left_cell_of_face, d_right_cell_of_face, d_xfc, d_yfc, d_zfc,
                                               d_xcc, d_ycc, d_zcc, d_deltl, d_deltr);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        if (nst < nBoundFace)
        {
            int nMid = nBoundFace;
            if (ned < nBoundFace) nMid = ned;
            KernelLaunchPara((void *)GPUFaceWeightOnBound, loopLen, 0, gridSize, blockSize, blockSizeLimit);
#ifdef KERNELLAUNCHTEST
            ReportKernelPara((void *)GPUFaceWeightOnBound, 0, gridSize, blockSize);
#endif

            GPUFaceWeightOnBound<<<gridSize, blockSize>>>(nst, nMid, d_deltl, d_deltr, d_boundaryType,
                                                          PHENGLEI::INTERFACE);
#ifdef KERNELLAUNCHTEST
            HANDLE_KERNEL_ERR();
#endif
        }
    }

    __global__ void GPUFaceWeight(const int nst, const int ned, const int *left_cell_of_face,
                                  const int *right_cell_of_face, const RDouble *xfc, const RDouble *yfc,
                                  const RDouble *zfc, const RDouble *xcc, const RDouble *ycc, const RDouble *zcc,
                                  RDouble *deltl, RDouble *deltr)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int iface = 0;
        int le, re, j;

        RDouble delt1, delt2, delta;
        RDouble dxl, dyl, dzl, dxr, dyr, dzr;

        for (iface = nst + i; iface < ned; iface += blockDim.x * gridDim.x)
        {
            le = left_cell_of_face[iface];
            re = right_cell_of_face[iface];
            j  = iface - nst;

            dxl = xfc[iface] - xcc[le];
            dyl = yfc[iface] - ycc[le];
            dzl = zfc[iface] - zcc[le];

            dxr = xfc[iface] - xcc[re];
            dyr = yfc[iface] - ycc[re];
            dzr = zfc[iface] - zcc[re];

            //! Left
            delt1    = sqrt(dxl * dxl + dyl * dyl + dzl * dzl);
            delt2    = sqrt(dxr * dxr + dyr * dyr + dzr * dzr);
            delta    = 1.0 / (delt1 + delt2 + SMALL);
            deltl[j] = delt2 * delta;
            deltr[j] = delt1 * delta;
        }
    }

    __global__ void GPUFaceWeightOnBound(const int nst, const int nMid, RDouble *deltl, RDouble *deltr,
                                         const int *boundaryType, int bctype_in)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int iface = 0;
        int j;
        for (iface = nst + i; iface < nMid; iface += blockDim.x * gridDim.x)
        {
            if (!(boundaryType[iface] == bctype_in || boundaryType[iface] < 0))
            {
                j        = iface - nst;
                deltl[j] = 1.0;
                deltr[j] = 0.0;
            }
        }
    }
} //! namespace GPUKernels
