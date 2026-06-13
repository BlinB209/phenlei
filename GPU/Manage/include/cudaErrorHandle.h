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
//! @file      cudaErrorHandle.h
//! @brief     Check CUDA error.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <iostream>
#define HANDLE_API_ERR(err) (handleAPIErr(err, __FILE__, __LINE__))
#define HANDLE_KERNEL_ERR() (handleKernelErr(__FILE__, __LINE__))
/***
 * this eror handle cannot work under the environment:
 * CUDA/9.0.176 & intelcompiler/18(exist at the same time)
 ***/
using namespace std;
//! handle error for cuda api
void handleAPIErr(cudaError_t err, char const *file, const int line);
//! handle error for cuda kernel
void handleKernelErr(char const *file, const int line);
