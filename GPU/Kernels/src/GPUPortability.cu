#include "GPUPortability.h"

using namespace PHSPACE;
namespace GPUKernels
{
    __device__ void atomicMax(double *addr, double val)
    {
        //!__device__ void atomicMax(RFloat *addr, RFloat val){
        unsigned long long int *addr_as_ull = (unsigned long long int *)(addr);
        unsigned long long int  old         = *addr_as_ull, assumed;
        do
        {
            assumed = old;
            old     = atomicCAS(addr_as_ull, assumed,
                                __double_as_longlong(kernelMAXDOUBLE(val, __longlong_as_double(assumed))));
        } while (assumed != old);
    }

    __device__ void atomicMin(double *addr, double val)
    {
        //!__device__ void atomicMin(RFloat *addr, RFloat val){
        unsigned long long int *addr_as_ull = (unsigned long long int *)(addr);
        unsigned long long int  old         = *addr_as_ull, assumed;
        do
        {
            assumed = old;
            old     = atomicCAS(addr_as_ull, assumed,
                                __double_as_longlong(kernelMINDOUBLE(val, __longlong_as_double(assumed))));
        } while (assumed != old);
    }

#if __CUDA_ARCH__ < 600
    __device__ double atomicAddTest(double *address, double val)
    {
        unsigned long long int *address_as_ull = (unsigned long long int *)address;
        unsigned long long int  old            = *address_as_ull, assumed;
        do
        {
            assumed = old;
            old     = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
        } while (assumed != old);
        return __longlong_as_double(old);
    }
#endif

    __device__ double kernelMAXDOUBLE(double a, double b)
    {
        //!__device__ RFloat kernelMAXDOUBLE(RFloat a, RFloat b){
        return ((a > b) ? a : b);
    }

    __device__ double kernelMINDOUBLE(double a, double b)
    {
        //!__device__ RFloat kernelMINDOUBLE(RFloat a, RFloat b){
        return ((a < b) ? a : b);
    }
} //! namespace GPUKernels
