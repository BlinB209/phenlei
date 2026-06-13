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
//! @file      OutputDebug.h
//! @brief     Debug function.
//! @author    Zhang Yong, Huang Can, Zhang Xi.

#include "stdio.h"
extern int   processNum;
extern int   d_iterStep;
extern int   d_rkStep;
extern char  hostName[60];
extern int   globalRank;
extern int   localRank;
extern int   globalZoneID;
extern int   localZoneID;
extern FILE *debugOutput;

void startDebugOutput();
void getCodeID(const int numberOfProcessors);
__global__ void kernelOutputRank();
