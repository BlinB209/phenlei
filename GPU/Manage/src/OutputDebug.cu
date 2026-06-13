#include <unistd.h>

#include <iostream>
#include "OutputDebug.h"
#include "TK_Exit.h"
#include "TK_Log.h"
#include "mpi.h"
#include "stdio.h"
using namespace PHSPACE;
using namespace PHMPI;
using namespace std;

int  processNum   = 0;
int  d_iterStep   = 0;
int  d_rkStep     = 0;
char hostName[60] = {'\0'};
int  globalRank;
int  localRank;
int  globalZoneID;
int  localZoneID;

__constant__ int d_globalLocalRank[4];

FILE *debugOutput;
void  startDebugOutput()
{
    if ((debugOutput = fopen("outputDebugResult", "w")) == NULL)
    {
        TK_Exit::ExceptionExit("Error: fail to open debugOutput!\n");
    }
}
//! get host name for multi process, the function should be placed after MPI_init
void getCodeID(const int numberOfProcessors)
{
    processNum = numberOfProcessors;
    int h_globalLocalRank[4];
    //! get host name firstly
    gethostname(hostName, sizeof(hostName));
    //! get global rank
    MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);
    //! get local rank
    localRank = -1;
    {
        MPI_Comm local_comm;
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, globalRank, MPI_INFO_NULL, &local_comm);
        MPI_Comm_rank(local_comm, &localRank);
        MPI_Comm_free(&local_comm);
    }
    h_globalLocalRank[0] = globalRank;
    h_globalLocalRank[1] = localRank;
    int nLocalZones      = GetNumberofLocalZones();

    for (int izone = 0; izone < nLocalZones; izone++)
    {
        int zoneID   = GetLocalZoneIDToGlobalZoneID(izone);
        globalZoneID = zoneID;
        localZoneID  = izone;
        ostringstream oss;
        oss << "Program is running in getZoneID, " << "globalRank = " << globalRank << ", localRank = " << localRank
            << ", globalZoneID = " << globalZoneID << ", localZoneID = " << localZoneID
            << ", nLocalZones = " << nLocalZones << "\n";
        PrintToWindow(oss);
    }
    ostringstream oss;
    oss << "Program is running in getCodeID, "
        << "globalRank = " << globalRank << ", localRank = " << localRank << ", globalZoneID = " << globalZoneID
        << ", localZoneID = " << localZoneID << ", nLocalZones = " << nLocalZones << "\n";
    PrintToWindow(oss);

    h_globalLocalRank[2] = globalZoneID;
    h_globalLocalRank[3] = localZoneID;
    cudaMemcpyToSymbol(d_globalLocalRank, h_globalLocalRank, 4 * sizeof(int));
    kernelOutputRank<<<1, 1>>>();
}

// void getZoneID() {
//   using namespace PHSPACE;  // namespace is used in PHengLEI. Otherwise,
//                             // function cannot be found!!
//   using namespace PHMPI;
//   int nLocalZones = GetNumberofLocalZones();
//   for (int izone = 0; izone < nLocalZones; izone++) {
//     int zoneID = GetLocalZoneIDToGlobalZoneID(izone);
//     printf(
//         "Program is running in getZoneID, zoneID = %d, izone = %d, localRank
//         = "
//         "%d, globalRank = %d\n",
//         zoneID, izone, localRank, globalRank);
//   }
// }

__global__ void kernelOutputRank()
{
    // if (d_globalLocalRank[0] == 0)
    // {
    //   cout << "d_globalRank = " << d_globalLocalRank[0] << " d_localRank = " <<
    //   d_globalLocalRank[1] << " globalZoneID = " << d_globalLocalRank[2] << "
    //   localZoneID = " << d_globalLocalRank[3] << endl;
    // }
    //  oss << "d_globalRank = " << d_globalLocalRank[0] << " d_localRank = " <<
    //  d_globalLocalRank[1] << " globalZoneID = " << d_globalLocalRank[2] << "
    //  localZoneID = " << d_globalLocalRank[3] << "\n"; PrintToWindow(oss);
    printf("d_globalRank = %d, d_localRank = %d, globalZoneID = %d, localZoneID = "
           "%d\n",
           d_globalLocalRank[0], d_globalLocalRank[1], d_globalLocalRank[2], d_globalLocalRank[3]);
}

// void stopDebugOutput() { fclose(debugOutput); }

// void outputHostValue(char const *functionName, char const *variableName,
//                      const int iter_step, const int rk_step,
//                      const double hostValue) {
//   fprintf(debugOutput, "%s: Iter %d, Rk %d, In %s, hostValue = %.30e\n",
//           variableName, iter_step, rk_step, functionName, hostValue);
//   printf("%s: Iter %d, Rk %d, In %s, hostValue = %.30e\n", variableName,
//          iter_step, rk_step, functionName, hostValue);
// }

// void outputMultiGPUHostValue(char const *functionName, char const
// *variableName,
//                              const int iter_step, const int rk_step,
//                              const double hostValue) {
//   // output hostValue
//   printf(
//       "%s: globalRank %d, localRank %d, Iter %d, Rk %d, In %s, hostValue = "
//       "%.30e\n",
//       variableName, globalRank, localRank, iter_step, rk_step, functionName,
//       hostValue);
// }
// /*
// void outputDeviceValue(char const *functionName, char const *variableName,
// const int iter_step, const int rk_step, const double deviceValue){
//         fprintf(debugOutput, "%s: Iter %d, Rk %d, In %s, deviceValue =
//         %.30e\n",
// variableName, iter_step, rk_step, functionName,  deviceValue);
// }
// */
// /*
// __device__ void outputDeviceValue(char const *functionName, char const
// *variableName, const int iter_step, const int rk_step, const double
// deviceValue){
// //    fprintf(debugOutput, "%s: Iter %d, Rk %d, In %s, deviceValue = %.30e\n",
// variableName, iter_step, rk_step, functionName,  deviceValue);
// }
// */
// __device__ void outputDeviceValue(char const *functionName,
//                                   char const *variableName,
//                                   const double deviceValue) {
//   printf("%s: In %s, deviceValue = %.30e\n", variableName, functionName,
//          deviceValue);
// }

// void outputMultiGPUHostValue(char const *functionName, char const
// *variableName,
//                              const int hostValue) {
//   printf(
//       "%s: in %s, globalRank %d, localRank %d, Iter %d, Rk %d hostValue =
//       %d\n", variableName, functionName, globalRank, localRank, d_iterStep,
//       d_rkStep, hostValue);
// }

// void outputMultiGPUHostValue(char const *functionName, char const
// *variableName,
//                              const double hostValue) {
//   printf(
//       "%s: in %s, globalRank %d, localRank %d, Iter %d, Rk %d hostValue = "
//       "%.30e\n",
//       variableName, functionName, globalRank, localRank, d_iterStep,
//       d_rkStep, hostValue);
// }

// __device__ void outputMultiGPUDeviceValue(char const *functionName,
//                                           char const *variableName,
//                                           const double deviceValue) {
//   printf("%s: In %s, in globalRank %d, in localRank %d, deviceValue =
//   %.30e\n",
//          variableName, functionName, d_globalLocalRank[0],
//          d_globalLocalRank[1], deviceValue);
// }

// __device__ void outputMultiGPUDeviceValue(char const *functionName,
//                                           char const *variableName,
//                                           const int deviceValue) {
//   printf("%s: In %s, in globalRank %d, in localRank %d, deviceValue = %d\n",
//          variableName, functionName, d_globalLocalRank[0],
//          d_globalLocalRank[1], deviceValue);
// }
// __device__ int getGPUGlobalRank() { return d_globalLocalRank[0]; }
// __device__ int getGPULocalRank() { return d_globalLocalRank[1]; }
// void outputHostGridValue() {
//   using namespace PHSPACE;  // namespace is used in PHengLEI. Otherwise,
//                             // function cannot be found!!
//   using namespace PHMPI;
//   using namespace GPUMemory;
//   using namespace GPUGeomVariables;
//   int nLocalZones = GetNumberofLocalZones();
//   for (int izone = 0; izone < nLocalZones; izone++) {
//     int zoneID = GetLocalZoneIDToGlobalZoneID(izone);
//     int level = 0;
//     UnstructGrid *grid = UnstructGridCast(GetGrid(zoneID, level));
//     RDouble *cpu_ycc = grid->GetCellCenterY();
//     printf("globalRank = %d, in host code, ycc[0] = %.30e\n", globalRank,
//            cpu_ycc[0]);
//     if (d_ycc) {
//       const int nTotalCell = grid->GetNTotalCell();
//       const int nBoundFace =
//           grid->GetNBoundFace();  // change the size of xcc ycc zcc
//       int nTotal = nTotalCell + nBoundFace;
//       RDouble *gpu_ycc = new RDouble[nTotal];
//       HANDLE_API_ERR(cudaMemcpy(gpu_ycc, d_ycc, nTotal * sizeof(RDouble),
//                                 cudaMemcpyDeviceToHost));
//       printf("globalRank = %d, in gpu code, ycc[0] = %.30e\n", globalRank,
//              gpu_ycc[0]);
//       delete[] gpu_ycc;
//     }
//   }
// }
// void outputCpuGpuValue() {
//   using namespace PHSPACE;  // namespace is used in PHengLEI. Otherwise,
//                             // function cannot be found!!
//   using namespace PHMPI;
//   using namespace GPUMemory;
//   using namespace GPUFlowVariables;
//   int nLocalZones = GetNumberofLocalZones();
//   for (int izone = 0; izone < nLocalZones; izone++) {
//     int zoneID = GetLocalZoneIDToGlobalZoneID(izone);
//     int level = 0;
//     UnstructGrid *grid = UnstructGridCast(GetGrid(zoneID, level));
//     RFloat **cpu_qNode = reinterpret_cast<RFloat
//     **>(grid->GetDataPtr("qnode")); printf("globalRank = %d, in host code,
//     qNode[0][5522] = %.30e\n",
//            globalRank, cpu_qNode[0][5522]);
//     if (d_qNode) {
//       int nl;
//       GlobalDataBase::GetData("nl", &nl, PHINT, 1);
//       int nchem;
//       GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);
//       int nEqn = nl + nchem;
//       int nTotalNode = grid->GetNTotalNode();
//       size_t nodeSize = nEqn * nTotalNode * sizeof(RFloat);
//       RFloat *gpu_qNode = new RFloat[nEqn * nTotalNode];
//       HANDLE_API_ERR(
//           cudaMemcpy(gpu_qNode, d_qNode, nodeSize, cudaMemcpyDeviceToHost));
//       printf("globalRank = %d, in device code, qNode[0][5522] = %.30e\n",
//              globalRank, gpu_qNode[0 * nTotalNode + 5522]);
//       delete[] gpu_qNode;
//     }
//   }
// }
// void outputCpuGpuQValue() {
//   using namespace PHSPACE;  // namespace is used in PHengLEI. Otherwise,
//                             // function cannot be found!!
//   using namespace PHMPI;
//   using namespace GPUMemory;
//   using namespace GPUFlowVariables;
//   int nLocalZones = GetNumberofLocalZones();
//   for (int izone = 0; izone < nLocalZones; izone++) {
//     int zoneID = GetLocalZoneIDToGlobalZoneID(izone);
//     int level = 0;
//     UnstructGrid *grid = UnstructGridCast(GetGrid(zoneID, level));
//     RFloat **cpu_q = reinterpret_cast<RFloat **>(grid->GetDataPtr("q"));
//     if (1 == globalRank)
//       printf("globalRank = %d, in host code, q[0][151726] = %.30e\n",
//              globalRank, cpu_q[0][151726]);
//     if (d_q_ns) {
//       int nl;
//       GlobalDataBase::GetData("nl", &nl, PHINT, 1);
//       int nchem;
//       GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);
//       int nEqn = nl + nchem;
//       int nBoundFace = grid->GetNBoundFace();
//       int nTotalCell = grid->GetNTotalCell();
//       int nTotal = nBoundFace + nTotalCell;
//       size_t sizeQ = nEqn * nTotal * sizeof(RFloat);
//       RFloat *gpu_q = new RFloat[nEqn * nTotal];
//       HANDLE_API_ERR(cudaMemcpy(gpu_q, d_q_ns, sizeQ,
//       cudaMemcpyDeviceToHost)); if (1 == globalRank)
//         printf("globalRank = %d, in device code, q[0][151726] = %.30e\n",
//                globalRank, gpu_q[0 * nTotal + 151726]);
//       delete[] gpu_q;
//     }
//   }
// }
// void outputCpuGpuqNodeValue() {
//   using namespace PHSPACE;  // namespace is used in PHengLEI. Otherwise,
//                             // function cannot be found!!
//   using namespace PHMPI;
//   using namespace GPUMemory;
//   using namespace GPUFlowVariables;
//   using namespace GPUGeomVariables;
//   int nLocalZones = GetNumberofLocalZones();
//   int globalRankCheckID = 2;
//   if (globalRankCheckID == globalRank) {
//     for (int izone = 0; izone < nLocalZones; izone++) {
//       int zoneID = GetLocalZoneIDToGlobalZoneID(izone);
//       int level = 0;
//       UnstructGrid *grid = UnstructGridCast(GetGrid(zoneID, level));
//       InterpointInfo *ipointinfo = grid->GetInterpointInfo();
//       int nIPoint = ipointinfo->GetNumberOfInterpoints();
//       int *interPoint2GlobalPoint = ipointinfo->GetInterPoint2GlobalPoint();
//       RFloat **cpu_qNode =
//           reinterpret_cast<RFloat **>(grid->GetDataPtr("qnode"));
//       RFloat **cpu_qInterPoint =
//           reinterpret_cast<RFloat **>(grid->GetDataPtr("qInterPoint"));
//       int globalPointCheckID = 1311;
//       int eqnCheckID = 0;
//       for (int localPoint = 0; localPoint < nIPoint; localPoint++) {
//         int globalPoint = interPoint2GlobalPoint[localPoint];
//         if (globalPointCheckID == globalPoint)
//           printf(
//               "globalRank = %d, in host code,  qInterPoint[%d][%d] =
//               %.30e\n", globalRank, eqnCheckID, localPoint,
//               cpu_qInterPoint[eqnCheckID][localPoint]);
//       }
//       // printf("globalRank = %d, in host code, localPoint = %d, globalPoint
//       =
//       // %d, qNode[0][500] = %.30e, qInterPoint[%d] = %.30e\n", globalRank,
//       // localPoint, globalPoint, cpu_qNode[0][500], localPoint,
//       // cpu_qInterPoint[0][localPoint]);
//       printf("globalRank = %d, in host code, qNode[%d][%d] = %.30e\n",
//              globalRank, eqnCheckID, globalPointCheckID,
//              cpu_qNode[eqnCheckID][globalPointCheckID]);
//       if ((d_qNode) && (d_qInterPoint) && (d_interPoint2GlobalPoint)) {
//         int nl;
//         GlobalDataBase::GetData("nl", &nl, PHINT, 1);
//         int nchem;
//         GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);
//         int neqn = nl + nchem;
//         int nTotalNode = grid->GetNTotalNode();
//         int nIPoint = ipointinfo->GetNumberOfInterpoints();
//         size_t sizeQNode = neqn * nTotalNode * sizeof(RFloat);
//         size_t sizeQIPoint = nIPoint * neqn * sizeof(RFloat);
//         size_t intPointSize = nIPoint * sizeof(int);
//         RFloat *gpu_qNode = new RFloat[neqn * nTotalNode];
//         RFloat *gpu_qInterPoint = new RFloat[nIPoint * neqn];
//         int *gpu_interPoint2GlobalPoint = new int[nIPoint];
//         HANDLE_API_ERR(
//             cudaMemcpy(gpu_qNode, d_qNode, sizeQNode,
//             cudaMemcpyDeviceToHost));
//         HANDLE_API_ERR(cudaMemcpy(gpu_qInterPoint, d_qInterPoint,
//         sizeQIPoint,
//                                   cudaMemcpyDeviceToHost));
//         HANDLE_API_ERR(cudaMemcpy(gpu_interPoint2GlobalPoint,
//                                   d_interPoint2GlobalPoint, intPointSize,
//                                   cudaMemcpyDeviceToHost));
//         for (int localPoint = 0; localPoint < nIPoint; localPoint++) {
//           int gpu_globalPoint = gpu_interPoint2GlobalPoint[localPoint];
//           if (globalPointCheckID == gpu_globalPoint)
//             printf(
//                 "globalRank = %d, in device code,  qInterPoint[%d][%d] = "
//                 "%.30e\n",
//                 globalRank, eqnCheckID, localPoint,
//                 gpu_qInterPoint[eqnCheckID * nIPoint + localPoint]);
//         }
//         // printf("globalRank = %d, in host code, localPoint = %d,
//         globalPoint =
//         // %d, qNode[0][500] = %.30e, qInterPoint[%d] = %.30e\n", globalRank,
//         // localPoint, globalPoint, cpu_qNode[0][500], localPoint,
//         // cpu_qInterPoint[0][localPoint]);
//         printf("globalRank = %d, in device code, qNode[%d][%d] = %.30e\n",
//                globalRank, eqnCheckID, globalPointCheckID,
//                gpu_qNode[eqnCheckID * nTotalNode + globalPointCheckID]);
//         delete[] gpu_qNode;
//         delete[] gpu_qInterPoint;
//         delete[] gpu_interPoint2GlobalPoint;
//       }
//     }
//   }  // end if (globalRankCheckID == globalRank)
// }
