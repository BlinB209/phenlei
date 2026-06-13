//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      GPUDeviceControl.h
//! @brief     Device control.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <iostream>

using namespace std;

extern int            GPUNum;
extern cudaDeviceProp GPUProp;

void GPUDevice(int &GPUNum, cudaDeviceProp &GPUProp);
void GetGPUNum(int &GPUNum);
void GPUInfoQuery(const int GPUId, cudaDeviceProp &GPUProp);
void GPUDeviceMultiGPU(int numberOfProcessors, int &GPUNum, cudaDeviceProp &GPUProp);

//!  calculate kernel parameter
//!  void KernelLaunchPara(void *kernel, int loopLen, size_t dynamicSMem, int
//! &gridSize, int &blockSize);

//! calculate kernel parameter
void KernelLaunchPara(void *kernel, int loopLen, size_t dynamicSMem, int &gridSize, int &blockSize,
                      int blockSizeLimit = 0);

void KernelLaunchParaDebug(void *kernel, int loopLen, size_t dynamicSMem, int &gridSize, int &blockSize,
                           int blockSizeLimit = 0);

//! calculate the kernel occupancy
void ReportKernelPara(void *kernel, size_t dynamicSMem, int gridSize, int blockSize);

//! calculate the kernel's blockSize limit
int KernelBlockSizeLimit(int regsPerThread, int residentBlockNo, cudaDeviceProp &GPUProp);

//! calculate the dynamic shared memory per block
size_t KerneldsMemPerBlock(int residentBlockNo, size_t dsMemPerThread, int &blockLimit, cudaDeviceProp &GPUProp);

//! calculate the kernel size
void KernelLaunchPara(void *kernel, int loopLen, int regsPerThread, int residentBlockNo, size_t dsMemPerThread,
                      size_t &dsMemPerBlock, int &gridSize, int &blockSize, cudaDeviceProp &GPUProp);
