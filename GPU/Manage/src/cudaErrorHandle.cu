#include "cudaErrorHandle.h"

 //! handle error for cuda api
void handleAPIErr(cudaError_t err, char const *file, const int line)
{
    if (err != cudaSuccess)
    {
        cout << "Error: " << cudaGetErrorString(err) << " in codes file " << file << " at line:# " << line << endl;
        exit(EXIT_FAILURE);
    }
}
 //! handle error for cuda kernel
void handleKernelErr(char const *file, const int line)
{
    cudaError_t err;
    int         tline = line;  //! the handleKernelErr is called next to Kernel function
     //!err= cudaGetLastError();
    cudaDeviceSynchronize();
    err = cudaPeekAtLastError();
    if (err != cudaSuccess)
    {
        cout << "Error: Fail to Launch kernel in file " << file << " at line " << tline << endl
             << "Error Description: " << cudaGetErrorString(err) << endl;
        exit(EXIT_FAILURE);
    }
}
