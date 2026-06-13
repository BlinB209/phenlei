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
//! @file      TimeMeasure.h
//! @brief     GPU time measure.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include <ctime>
#include <cuda_runtime_api.h>

extern double dtime;

extern cudaEvent_t kernel_start;
extern cudaEvent_t kernel_stop;
//！ stage 0 means start, 1 means end
void TimeMeasure(int stage);
void GPUTimeMeasure(int stage);
void PrintAcmltTime(double dtime);
