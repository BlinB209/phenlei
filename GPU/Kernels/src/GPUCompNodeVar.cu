#include "GPUCompNodeVar.h"
#include "GPUFixBCNodeVar.h"
#include "cudaErrorHandle.h"
#include "GPUDeviceControl.h"
#include "Geo_SimpleBC.h"
#if __CUDA_ARCH__ < 600
#include "GPUPortability.h"
#endif
//!using namespace GPUMemory;
//!using namespace GPUFlowVariables;
//!using namespace GPUControlVariables;
using namespace std;
using namespace PHENGLEI;
namespace GPUKernels
{
    void CallGPUCompNodeVarByGradient(const string d_gradient_field_proxy, const int d_gradient_var_index,
                                      const int nTotalNode, const int nBoundFace, const int nTotalFace,
                                      const int nTotalCell)
    {
#ifdef KERNELLAUNCHTEST
        cout << "Test in CallGPUCompNodeVarByGradient with gradient_field_proxy=" << d_gradient_field_proxy
             << "and d_gradient_var_index=" << d_gradient_var_index << endl;
#endif
        RFloat *q;
        int     index;
        //!RFloat *q_n;
        //!int *n_count;
        //!q_n = d_q_n_tmp;
        //!n_count = d_n_count_tmp;

        using namespace GPUMemory;
        using namespace GPUFlowVariables;
        using namespace GPUControlVariables;
        using namespace GPUGeomVariables;
        //!GPUNcountQnInit<<<2092, 640>>>(nTotalNode, q_n, n_count);
        /*
        int gridSize = 0;
        int blockSize = 0;
        int minGridSize = 0;
        HANDLE_API_ERR(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, GPUNcountQnInit, 0, nTotalNode));
        gridSize = (minGridSize + blockSize - 1) / blockSize;
        #ifdef KERNELLAUNCHTEST
            printf("For GPUNcountQnInit, gridSize = %d, blockSize = %d\n", gridSize, blockSize);
        #endif
        //!GPUNcountQnInit<<<2092, 640>>>(nTotalNode, d_q_n_tmp, d_n_count_tmp);
        */
        int gridSize  = 1;
        int blockSize = 1;
        int loopLen   = nTotalNode;
        KernelLaunchPara((void *)GPUNcountQnInit, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUNcountQnInit, 0, gridSize, blockSize);
#endif

        GPUNcountQnInit<<<gridSize, blockSize>>>(nTotalNode, d_q_n_tmp, d_n_count_tmp);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        //!set zero on d_q_n_double and d_n_count_tmp
        //!GPUNcountQnInit<<<2092, 640>>>(nTotalNode, d_q_n_double, d_n_count_tmp);
        //!
        if (d_gradient_field_proxy == "d_q_proxy")
        {
            q     = d_q_ns;
            index = d_gradient_var_index;
        }
        else if (d_gradient_field_proxy == "d_t_proxy")
        {
            q     = d_t_proxy;
            index = d_gradient_var_index;
        }
        else if (d_gradient_field_proxy == "d_q_turb_proxy")
        {
            q     = d_q_turb_proxy;
            index = d_gradient_var_index;
        }
        else if (d_gradient_field_proxy == "d_velocity_proxy")
        {
            q     = d_vel_proxy;
            index = d_gradient_var_index;
        }
        else
        {
            return;
        }
        //!calculate node values on boundary faces
        /*
        HANDLE_API_ERR(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, GPUBoundaryFaceNodeCal, 0, nBoundFace));
        gridSize = (minGridSize + blockSize - 1) / blockSize;
        #ifdef KERNELLAUNCHTEST
            printf("For GPUBoundaryFaceNodeCal, gridSize = %d, blockSize = %d\n", gridSize, blockSize);
        #endif
        //!GPUBoundaryFaceNodeCal<<<2092, 640>>>(index, nBoundFace, 
        */
        loopLen = nBoundFace;
        KernelLaunchPara((void *)GPUBoundaryFaceNodeCal, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUBoundaryFaceNodeCal, 0, gridSize, blockSize);
#endif

        GPUBoundaryFaceNodeCal<<<gridSize, blockSize>>>(index, nBoundFace, nTotalCell, d_left_cell_of_face,
                                                        d_node_number_of_each_face, d_face2node, d_nodePosiFace, d_x,
                                                        d_y, d_z, d_xcc, d_ycc, d_zcc, d_q_n_tmp,
                                                        //!d_xcc, d_ycc, d_zcc, d_q_n_double,
                                                        d_n_count_tmp, q);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        //!calculate node values on interior faces
        /*
        HANDLE_API_ERR(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, GPUInteriorFaceNodeCal, 0, nTotalFace-nBoundFace));
        gridSize = (minGridSize + blockSize - 1) / blockSize;
        #ifdef KERNELLAUNCHTEST
            printf("For GPUInteriorFaceNodeCal, gridSize = %d, blockSize = %d\n", gridSize, blockSize);
        #endif
        //!GPUInteriorFaceNodeCal<<<2092, 640>>>(index, nBoundFace, 
        */
        loopLen = nTotalFace - nBoundFace;
        KernelLaunchPara((void *)GPUInteriorFaceNodeCal, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUInteriorFaceNodeCal, 0, gridSize, blockSize);
#endif

        GPUInteriorFaceNodeCal<<<gridSize, blockSize>>>(index, nBoundFace, nTotalCell, nTotalFace, d_left_cell_of_face,
                                                        d_right_cell_of_face, d_node_number_of_each_face, d_face2node,
                                                        //!d_nodePosiFace, d_q_n_double,
                                                        d_nodePosiFace, d_q_n_tmp, d_n_count_tmp, q);

#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif

        CallGPUFixBCNodeVarByCompNodeVar(index, nTotalNode, nBoundFace, nTotalCell, d_left_cell_of_face,
                                         d_right_cell_of_face, d_face2node, d_node_number_of_each_face,
                                         //!d_boundaryType, d_nodePosiFace, q, d_q_n_double, d_n_count_tmp,
                                         d_boundaryType, d_nodePosiFace, q, d_q_n_tmp, d_n_count_tmp,
                                         PHENGLEI::SYMMETRY, true);

        CallGPUFixBCNodeVarByCompNodeVar(index, nTotalNode, nBoundFace, nTotalCell, d_left_cell_of_face,
                                         d_right_cell_of_face, d_face2node, d_node_number_of_each_face,
                                         //!d_boundaryType, d_nodePosiFace, q, d_q_n_double, d_n_count_tmp,
                                         d_boundaryType, d_nodePosiFace, q, d_q_n_tmp, d_n_count_tmp,
                                         PHENGLEI::SOLID_SURFACE, true);

        CallGPUFixBCNodeVarByCompNodeVar(index, nTotalNode, nBoundFace, nTotalCell, d_left_cell_of_face,
                                         d_right_cell_of_face, d_face2node, d_node_number_of_each_face,
                                         //!d_boundaryType, d_nodePosiFace, q, d_q_n_double, d_n_count_tmp,
                                         d_boundaryType, d_nodePosiFace, q, d_q_n_tmp, d_n_count_tmp,
                                         PHENGLEI::INTERFACE, true);

        CallGPUFixBCNodeVarByCompNodeVar(index, nTotalNode, nBoundFace, nTotalCell, d_left_cell_of_face,
                                         d_right_cell_of_face, d_face2node, d_node_number_of_each_face,
                                         //!d_boundaryType, d_nodePosiFace, q, d_q_n_double, d_n_count_tmp,
                                         d_boundaryType, d_nodePosiFace, q, d_q_n_tmp, d_n_count_tmp,
                                         PHENGLEI::FARFIELD, true);

        //!GPUNodeVarAve<<<2092, 640>>>(nTotalNode, d_q_n_double, d_n_count_tmp, SMALL);
        loopLen = nTotalNode;
        KernelLaunchPara((void *)GPUNodeVarAve, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUNodeVarAve, 0, gridSize, blockSize);
#endif
        GPUNodeVarAve<<<gridSize, blockSize>>>(nTotalNode, d_q_n_tmp, d_n_count_tmp, SMALL);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
        loopLen = nBoundFace;
        KernelLaunchPara((void *)GPUNodeBoundaryCalc, loopLen, 0, gridSize, blockSize);
#ifdef KERNELLAUNCHTEST
        ReportKernelPara((void *)GPUNodeBoundaryCalc, 0, gridSize, blockSize);
#endif
        gridSize  = 1;
        blockSize = 1;
        GPUNodeBoundaryCalc<<<gridSize, blockSize>>>(
            nBoundFace, d_left_cell_of_face, d_right_cell_of_face, q, d_boundaryType, PHENGLEI::SOLID_SURFACE,
            PHENGLEI::SYMMETRY, d_face2node, d_q_n_tmp, d_node_number_of_each_face, d_nodePosiFace, nTotalCell, index);
#ifdef KERNELLAUNCHTEST
        HANDLE_KERNEL_ERR();
#endif
    }

    __global__ void GPUNcountQnInit(const int nTotalNode, RFloat *q_n, int *n_count)
    {
        //!    __global__ void GPUNcountQnInit(const int nTotalNode, double *q_n, int *n_count){
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        int j = 0;
        for (j = i; j < nTotalNode; j += blockDim.x * gridDim.x)
        {
            q_n[j]     = 0.0;
            n_count[j] = 0;
        }
    }

    __global__ void GPUBoundaryFaceNodeCal(const int index, const int nBoundFace, const int nTotalCell,
                                           const int *left_cell_of_face, const int *node_number_of_each_face,
                                           const int *face2node, const long long int *nodepos, RDouble *x, RDouble *y,
                                           RDouble *z, RDouble *xcc, RDouble *ycc, RDouble *zcc, RFloat *q_n,
                                           //!RDouble * xcc, RDouble * ycc, RDouble * zcc, double * q_n,
                                           int *n_count, RFloat *q)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int iface = 0;
        int le, re, pt, j;
        //!double qle, qre;
        RFloat qle, qre;
        for (iface = i; iface < nBoundFace; iface += gridDim.x * blockDim.x)
        {
            le = left_cell_of_face[iface];
            re = iface + nTotalCell;
            //!qle = q[le];
            //!qre = q[re];
            qle = q[le + index * (nTotalCell + nBoundFace)];
            qre = q[re + index * (nTotalCell + nBoundFace)];
            for (j = 0; j < node_number_of_each_face[iface]; j++)
            {
                pt = face2node[nodepos[iface] + j];
#if __CUDA_ARCH__ < 600
                atomicAddTest(q_n + pt, qle);
#else
                atomicAdd(q_n + pt, qle);
#endif
                atomicAdd(n_count + pt, 1);
#if __CUDA_ARCH__ < 600
                atomicAddTest(q_n + pt, qre);
#else
                atomicAdd(q_n + pt, qre);
#endif
                atomicAdd(n_count + pt, 1);
            }
        }
    }

    __global__ void GPUInteriorFaceNodeCal(const int index, const int nBoundFace, const int nTotalCell,
                                           const int nTotalFace, const int *left_cell_of_face,
                                           const int *right_cell_of_face, const int *node_number_of_each_face,
                                           const int *face2node, const long long int *nodepos,
                                           //!double * q_n, int * n_count, RFloat * q){
                                           RFloat *q_n, int *n_count, RFloat *q)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int iface = 0;
        int le, re, pt, j;
        //!double qle, qre;
        RFloat qle, qre;
        for (iface = i + nBoundFace; iface < nTotalFace; iface += gridDim.x * blockDim.x)
        {
            le  = left_cell_of_face[iface];
            re  = right_cell_of_face[iface];
            qle = q[le + index * (nTotalCell + nBoundFace)];
            qre = q[re + index * (nTotalCell + nBoundFace)];
            for (j = 0; j < node_number_of_each_face[iface]; j++)
            {
                pt = face2node[nodepos[iface] + j];

#if __CUDA_ARCH__ < 600
                atomicAddTest(q_n + pt, qle);
#else
                atomicAdd(q_n + pt, qle);
#endif
                atomicAdd(n_count + pt, 1);
#if __CUDA_ARCH__ < 600
                atomicAddTest(q_n + pt, qre);
#else
                atomicAdd(q_n + pt, qre);
#endif
                atomicAdd(n_count + pt, 1);
            }
        }
    }

    //!__global__ void GPUNodeVarAve(const int nTotalNode, double * q_n, int * n_count, const double SMALL){
    __global__ void GPUNodeVarAve(const int nTotalNode, RFloat *q_n, int *n_count, const double SMALL)
    {
        int i     = blockDim.x * blockIdx.x + threadIdx.x;
        int inode = 0;
        for (inode = i; inode < nTotalNode; inode += blockDim.x * gridDim.x)
        {
            q_n[inode] /= (n_count[inode] + SMALL);
        }
    }
    __global__ void GPUNodeBoundaryCalc(const int nBoundFace, const int *leftCellOfFace, const int *rightCellOfFace,
                                        const RFloat *q, const int *d_boundaryType, const int SOLID_SURFACE,
                                        const int SYMMETRY, const int *face2node, RFloat *qNode,
                                        const int *nodeNumberOfEachFace, const long long int *nodepos,
                                        const int nTotalCell, int index)
    {
        //! Modify the velocity node value on the wall.
        const RDouble VELOCITY_ZERO = 1.0e-16;
        //!BCRecord **bcRecord = grid->GetBCRecord();
        //!int nodePosition = 0;
        int i = blockDim.x * blockIdx.x + threadIdx.x;
        for (int iFace = i; iFace < nBoundFace; iFace += blockDim.x * gridDim.x)
        {
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];

            RDouble leftValue  = q[le + index * (nTotalCell + nBoundFace)];
            RDouble rightValue = q[re + index * (nTotalCell + nBoundFace)];

            int bcType = d_boundaryType[iFace];
            if (bcType == SOLID_SURFACE)
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++jNode)
                {
                    int point = face2node[nodepos[iFace] + jNode];
                    if (leftValue * rightValue <= VELOCITY_ZERO)
                    {
                        //! it is velocity.
                        qNode[point] = 0;
                    }
                    else
                    {
                        qNode[point] = leftValue;
                    }
                }
            }
            else if (bcType == SYMMETRY)
            {
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++jNode)
                {
                    int point = face2node[nodepos[iFace] + jNode];

                    qNode[point] = 0.5 * (leftValue + rightValue);
                }
            }
            //!atomicAdd(&nodePosition,nodeNumberOfEachFace[iFace]);
        }
    }

} //! namespace GPUKernels
