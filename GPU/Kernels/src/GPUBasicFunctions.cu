#include "GPUBasicFunctions.h"
//！ get the smaller one
__device__ RFloat GPUMIN(RFloat a, RFloat b) { return (a > b ? b : a); }
__device__ RFloat GPUMAX(RFloat a, RFloat b) { return (a > b ? a : b); }
__device__ RFloat GPUABS(RFloat a) { return (a < 0.0 ? -a : a); }
