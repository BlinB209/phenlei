#include "TimeMeasure.h"
#include <iostream>
#include <sys/time.h>
#include <stdio.h>
#include <iostream>

using namespace std;

cudaEvent_t kernel_start;
cudaEvent_t kernel_stop;

double dtime = 0.0;

/* The accurancy is millisecond */
void TimeMeasure(int stage)
{
    static double  start, stop;
    struct timeval mytime;
    double         epst;

    if (stage == 0) //! start
    {
      dtime = 0.0;
      cudaDeviceSynchronize(); //! it's not safe for multistream condition
      gettimeofday(&mytime, (struct timezone *)0);
      start = (double)(mytime.tv_sec + mytime.tv_usec * 1.0e-6);
      return;
    }
    else if (stage == 1) //! end
    {
      cudaDeviceSynchronize();
      gettimeofday(&mytime, (struct timezone *)0);
      stop = (double)(mytime.tv_sec + mytime.tv_usec * 1.0e-6);
      epst = stop - start;
      //!cout<<"Time elapsed: "<<epst<<" s"<<endl;

      dtime += epst;
      return;
    }
}
void GPUTimeMeasure(int stage)
{
    float          epst;
    struct timeval mytime;

    if (stage == 0) //! start
    {
      cudaDeviceSynchronize(); //! it's not safe for multistream condition
      cudaEventCreate(&kernel_start);
      cudaEventCreate(&kernel_stop);
      cudaEventRecord(kernel_start, 0);
      cudaEventSynchronize(kernel_start);

      return;
    }
    else if (stage == 1) //!end
    {
      cudaEventRecord(kernel_stop, 0);
      cudaEventSynchronize(kernel_stop);
      cudaEventElapsedTime(&epst, kernel_start, kernel_stop);
      epst /= 1000.0; //! from ms to s
      dtime += epst;

      cudaEventDestroy(kernel_start);
      cudaEventDestroy(kernel_stop);
      return;
    }
}
void PrintAcmltTime(double dtime)
{
    cout << "The kernel execution time: " << dtime << " s" << endl
        << "---------------------------------------------" << endl;
}
