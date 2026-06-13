#include "GPUBasicFunctionsPart2.h"
__device__ RFloat GPUSQR(RFloat a) { return (a * a); }

__device__ RFloat GPUSQR(RFloat a, RFloat b, RFloat c) { return (a * a + b * b + c * c); }
