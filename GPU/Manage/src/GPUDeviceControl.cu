#include <unistd.h>
#include <cstdlib>
#include <string>
#include "GPUDeviceControl.h"
#include "OutputDebug.h"
#include "TK_Exit.h"
#include "TK_Log.h"
#include "cudaErrorHandle.h"
#include "mpi.h"
#include "stdio.h"

using namespace std;
using namespace PHSPACE;

//! define the external variables of GPUManage namespace
int GPUNum = 0;

cudaDeviceProp GPUProp = cudaDeviceProp();

 //! get GPU count on the node
void GetGPUNum(int &GPUNum)
{
    HANDLE_API_ERR(cudaGetDeviceCount(&GPUNum));
    PrintToWindow("The numberOfGPUDevice = ", GPUNum, "for server processor", "\n");
    if (0 == GPUNum)
    {
        ostringstream oss;
        oss << "Error: No GPU Device Detected on This Compute Node." << endl
            << "Please Change the Platform or Turn off CUDA Compile Mode" << endl
            << "Exit!" << endl;
        TK_Exit::ExceptionExit(oss);
    }
}
 //! get detailed info about GPU card
void GPUInfoQuery(const int GPUId, cudaDeviceProp &GPUProp)
{
    HANDLE_API_ERR(cudaGetDeviceProperties(&GPUProp, GPUId));
}

 //! for multi-GPUs
void GPUDeviceMultiGPU(int numberOfProcessors, int &GPUNum, cudaDeviceProp &GPUProp)
{
    GetGPUNum(GPUNum);
     //! only one process
    if (1 == numberOfProcessors)
    {
        int GPUId = GPUNum - 1;
        GPUInfoQuery(GPUId, GPUProp);
        //! set the latest GPU as the working GPU
        cudaSetDevice(GPUId);
    }
     //! many processes: one process for one GPU
    else
    {
         //! get the global rank
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        int local_rank = -1;
        {
            MPI_Comm local_comm;
            MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &local_comm);
            MPI_Comm_rank(local_comm, &local_rank);
            MPI_Comm_free(&local_comm);
        }

        char hostname[32];
        if (gethostname(hostname, sizeof(hostname)))
        {
            TK_Exit::ExceptionExit("Error: gethostname calling error\n");
        }

        ostringstream oss;
        oss << "Program is running in MultiGPUDevice, "
            << "rank = " << rank << ", local_rank = " << local_rank << ", numberOfProcessors = " << numberOfProcessors
            << ", hostname = " << hostname << endl;
        PrintToWindow(oss);

        GPUInfoQuery(local_rank, GPUProp);
        cudaSetDevice(local_rank);
        getCodeID(numberOfProcessors);
    }
}

 //! report the occupancy for kernel
void ReportKernelPara(void *kernel, size_t dynamicSMem, int gridSize, int blockSize)
{
    //! int device;
    //! cudaDeviceProp prop;

    int    numBlocks;
    int    activeWarps;
    int    maxWarps;
    double occupancy;

    //! cout<<"warpSize from GPUProp: "<<GPUProp.warpSize<<endl;
    //! HANDLE_API_ERR(cudaGetDevice(&device));
    //! HANDLE_API_ERR(cudaGetDeviceProperties(&prop, device));

    //! cout<<"warpSize from prop: "<<prop.warpSize<<endl;
    HANDLE_API_ERR(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocks, kernel, blockSize, dynamicSMem));

    activeWarps = numBlocks * blockSize / GPUProp.warpSize;
    maxWarps    = GPUProp.maxThreadsPerMultiProcessor / GPUProp.warpSize;
    /*
      cout<<"gridSize from cudaOccupancyMaxActiveBlocksPerMultiprocessor:
      "<<numBlocks<<endl
          <<"gridSize from formula: "<<gridSize<<endl
          <<"blockSize: "<<blockSize<<endl
          <<"warpSize: "<<prop.warpSize<<endl
          <<"maxThreadsPerMultiProcessor:
          "<<prop.maxThreadsPerMultiProcessor<<endl;
          <<"activeWarps: "<<activeWarps<<endl
          <<"maxWarps"<<maxWarps<<endl;
  */
    occupancy = (double)activeWarps / maxWarps;
    cout << "For Kernel(formula): gridSize= " << gridSize << " blockSize= " << blockSize << endl;
    cout << " Occupancy : " << occupancy << endl;
}
/* calculate kernel parameter */
void KernelLaunchPara(void *kernel, int loopLen, size_t dynamicSMem, int &gridSize, int &blockSize, int blockSizeLimit)
{
    int minGridSize; //! minimal limitation for gridSize
    //! calculate the kernel launch parameter
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, kernel, dynamicSMem, blockSizeLimit);
    if (0 == loopLen)
    {
        gridSize = 1;
        return;
    }
    gridSize = (loopLen + blockSize - 1) / blockSize;
    if (gridSize > GPUProp.maxGridSize[0])
        //! reset gridSize value according the limitation
        gridSize = GPUProp.maxGridSize[0];
    //! test
    //! cout<<"minGridSize= "<<minGridSize<<endl;
}
void KernelLaunchParaDebug(void *kernel, int loopLen, size_t dynamicSMem, int &gridSize, int &blockSize,
                           int blockSizeLimit)
{
    int minGridSize; //! minimal limitation for gridSize
    //! calculate the kernel launch parameter
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, kernel, dynamicSMem, blockSizeLimit);
    if (0 == loopLen)
    {
        gridSize = 1;
        return;
    }
    gridSize = (loopLen + blockSize - 1) / blockSize;
    printf("In KernelLaunchParaDebug, minGridSize = %d, blockSize = %d, gridSize = "
           "% d, GPUProp.maxGridSize[0] = % d\n ",
           minGridSize, blockSize, gridSize, GPUProp.maxGridSize[0]);
    if (gridSize > GPUProp.maxGridSize[0])
        //! reset gridSize value according the limitation
        gridSize = GPUProp.maxGridSize[0];
    //! test
    //! cout<<"minGridSize= "<<minGridSize<<endl;
}
/* calculate the kernel's blockSize limit*/
int KernelBlockSizeLimit(int regsPerThread, int residentBlockNo, cudaDeviceProp &GPUProp)
{
    int blockLimit, blockLimitTmp;
    int regsPerSM = GPUProp.regsPerBlock;
    //! calculate according to the register volume on the SM
    blockLimitTmp = regsPerSM / (residentBlockNo * regsPerThread);
    blockLimitTmp /= 32;
    //! make the blockLimit is integral multiple of 32, a warpsize
    blockLimit = blockLimitTmp * 32;
    return blockLimit;
}
/* calculate the dynamic shared memory per block*/
size_t KerneldsMemPerBlock(int residentBlockNo, size_t dsMemPerThread, int &blockSizeLimit, cudaDeviceProp &GPUProp)
{
    size_t dsMemPerBlock;
    //! shared memory on a SM(B)
    size_t sharedMemPerSM = GPUProp.sharedMemPerBlock;
    int    blockSizeLimitTmp; //! blockSizeLimit according to shared  memory use
    if (dsMemPerThread == 0) return 0;
    blockSizeLimitTmp = sharedMemPerSM / (residentBlockNo * dsMemPerThread);
    blockSizeLimitTmp /= 32;
    //! make the blockSizeLimit is integral multiple of 32, a warpsize
    blockSizeLimitTmp *= 32;
    //! choose the smaller limit
    blockSizeLimit = (blockSizeLimitTmp < blockSizeLimit) ? blockSizeLimitTmp : blockSizeLimit;
    dsMemPerBlock  = blockSizeLimit * dsMemPerThread;
    return dsMemPerBlock;
}
/* calculate the kernel size*/
void KernelLaunchPara(void *kernel, int loopLen, int regsPerThread, int residentBlockNo, size_t dsMemPerThread,
                      size_t &dsMemPerBlock, int &gridSize, int &blockSize, cudaDeviceProp &GPUProp)
{
    int blockSizeLimit;
    //! calculate the blockSizeLimit by register use
    blockSizeLimit = KernelBlockSizeLimit(regsPerThread, residentBlockNo, GPUProp);
    //! adjust the blockSizeLimit by shared memory use, and return the dynamic
    //! shared memory of a threads block
    dsMemPerBlock = KerneldsMemPerBlock(residentBlockNo, dsMemPerThread, blockSizeLimit, GPUProp);
    //! calculate the gridSize and blockSize with blockSizeLimit and dsMemPerBlock
    int minGridSize; //! minimal limitation for gridSize
    //! calculate the kernel launch parameter
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, kernel, dsMemPerBlock, blockSizeLimit);
    if (0 == loopLen)
    {
        gridSize = 1;
        return;
    }
    gridSize = (loopLen + blockSize - 1) / blockSize;
    if (gridSize > GPUProp.maxGridSize[0])
        //! reset gridSize value according the girdSize limitation
        gridSize = GPUProp.maxGridSize[0];
}
