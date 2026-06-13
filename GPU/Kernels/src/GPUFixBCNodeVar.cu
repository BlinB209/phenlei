#include "GPUFixBCNodeVar.h"
#include "GPUDeviceControl.h"
#include "Geo_SimpleBC.h"
#if __CUDA_ARCH__ < 600
#include "GPUPortability.h"
#endif

using namespace std;
namespace GPUKernels
{
    void CallGPUFixBCNodeVarByCompNodeVar(const int index, const int nTotalNode, const int nBoundFace,
                                          const int nTotalCell, const int *left_cell_of_face,
                                          const int *right_cell_of_face, const int *face2node,
                                          const int *node_number_of_each_face, const int *boundaryType,
                                          const long long int *nodepos, RFloat *q, RFloat *q_n, int *n_count,
                                          const int bctype_in, const bool twoside)
    { //!Try different type on q_n
#ifdef KERNELLAUNCHTEST
        cout << "Test CallGPUFixBCNodeVarByCompNodeVar with index=" << index << ", bctype_in=" << bctype_in
             << ", twoside=" << twoside << endl;
#endif

        //!bctype_in == FANTASY::INTERFACE is not considered, because only "Pressure Far Field", "Symmetry", and "Wall" are used, corresponding to bcType = 4, 3, 2
        /*
        int gridSize = 0;
        int blockSize = 0;
        int minGridSize = 0;
        HANDLE_API_ERR(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, GPUFixBCNodeVarInit, 0, nBoundFace));
        gridSize = (minGridSize + blockSize - 1) / blockSize;
        #ifdef KERNELLAUNCHTEST
        //!    KernelActiveOccupancy("GPUFixBCNodeVarInit", gridSize, blockSize);
            printf("For GPUFixBCNodeVarInit, gridSize = %d, blockSize = %d\n", gridSize, blockSize);
        #endif
        //!GPUFixBCNodeVarInit<<<2092, 640>>>(nBoundFace, nTotalCell, left_cell_of_face, 
        */
        int gridSize       = 1;
        int blockSize      = 1;
        int loopLen        = nBoundFace;
        int blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUFixBCNodeVarInit, loopLen, 0, gridSize, blockSize, blockSizeLimit);

#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUFixBCNodeVarInit, 0, gridSize, blockSize);
#endif

        GPUFixBCNodeVarInit<<<gridSize, blockSize>>>(nBoundFace, nTotalCell, left_cell_of_face, face2node,
                                                     node_number_of_each_face, boundaryType, nodepos, q_n, n_count,
                                                     bctype_in);

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
        /*
        HANDLE_API_ERR(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, GPUFixBCNodeVarCal, 0, nBoundFace));
        gridSize = (minGridSize + blockSize - 1) / blockSize;

        #ifdef KERNELLAUNCHTEST
            printf("For GPUFixBCNodeVarCal, gridSize = %d, blockSize = %d\n", gridSize, blockSize);
        #endif
        //!GPUFixBCNodeVarCal<<<2092, 640>>>(index, nBoundFace, nTotalCell, left_cell_of_face, 
        */
        loopLen        = nBoundFace;
        blockSizeLimit = 0;
        KernelLaunchPara((void *)GPUFixBCNodeVarCal, loopLen, 0, gridSize, blockSize, blockSizeLimit);

#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUFixBCNodeVarCal, 0, gridSize, blockSize);
#endif

        GPUFixBCNodeVarCal<<<gridSize, blockSize>>>(index, nBoundFace, nTotalCell, left_cell_of_face, face2node,
                                                    node_number_of_each_face, boundaryType, nodepos, q, q_n, n_count,
                                                    bctype_in, twoside);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUFixBCNodeVarInit(const int nBoundFace, const int nTotalCell, const int *left_cell_of_face,
                                        const int *face2node, const int *node_number_of_each_face,
                                        const int *boundaryType, const long long int *nodepos, RFloat *q_n,
                                        int *n_count, const int bctype_in)
    { //!Try different type on q_n

        int le, re, pt, bctype;
        int i = blockDim.x * blockIdx.x + threadIdx.x;
        for (int iface = i; iface < nBoundFace; iface += blockDim.x * gridDim.x)
        {
            le     = left_cell_of_face[iface];
            re     = iface + nTotalCell;
            bctype = boundaryType[iface];
            if (bctype == bctype_in)
            {
                for (int j = 0; j < node_number_of_each_face[iface]; ++j)
                {
                    pt = face2node[nodepos[iface] + j];

                    //!The judge is deleted because INTERFACE does not exist in the case
                    //!if ( ( bctype_in == FANTASY::INTERFACE ) &&
                    //!( corner_point[ pt ] == 2 ) ) continue;

                    q_n[pt]     = 0.0;
                    n_count[pt] = 0;
                }
            }
        }
    }

    __global__ void GPUFixBCNodeVarCal(const int index, const int nBoundFace, const int nTotalCell,
                                       const int *left_cell_of_face, const int *face2node,
                                       const int *node_number_of_each_face, const int *boundaryType,
                                       const long long int *nodepos, RFloat *q, RFloat *q_n, int *n_count,
                                       const int bctype_in, const bool twoside)
    {
        //!Try different type on q_n
        int le, re, pt, bctype;
        //!double qle, qre;
        RFloat qle, qre; //!Try different type
        int    i = blockDim.x * blockIdx.x + threadIdx.x;
        for (int iface = i; iface < nBoundFace; iface += blockDim.x * gridDim.x)
        {
            le     = left_cell_of_face[iface];
            re     = iface + nTotalCell;
            bctype = boundaryType[iface];
            if (bctype == bctype_in)
            {
                for (int j = 0; j < node_number_of_each_face[iface]; ++j)
                {
                    pt = face2node[nodepos[iface] + j];

                    //!The judge is deleted because INTERFACE does not exist in the case
                    //!if ( ( bctype_in == FANTASY::INTERFACE ) &&
                    //!( corner_point[ pt ] == 2 ) ) continue;

                    qle = q[le + index * (nTotalCell + nBoundFace)];
#if __CUDA_ARCH__ < 600
                    atomicAddTest(q_n + pt, qle);
#else
                    atomicAdd(q_n + pt, qle);
#endif
                    atomicAdd(n_count + pt, 1);
                    if (twoside)
                    {
                        //!q_n[ pt ] += q[ re ];
                        //!++ n_count[ pt ];
                        qre = q[re + index * (nTotalCell + nBoundFace)];
#if __CUDA_ARCH__ < 600
                        atomicAddTest(q_n + pt, qre);
#else
                        atomicAdd(q_n + pt, qre);
#endif
                        atomicAdd(n_count + pt, 1);
                    }
                }
            }
        }
    }
} //! namespace GPUKernels
