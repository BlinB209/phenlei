#include "TK_Time.h"
#include "TK_Log.h"
#ifndef WIN32
#include <sys/time.h>
#include <unistd.h>
#include <sys/times.h>
#endif

namespace PHSPACE
{
namespace TIME_SPACE
{

LIB_EXPORT double GetWallTime()
{
    static time_t first_time[2];

    static int count = 1;

    time_t now_time[2];

    double GreenTime;

    double iCLK_TCK = 1.0;

#ifndef WIN32
    struct timeval va1;
    int ErrCode = 0;
    ErrCode = gettimeofday(&va1, (struct timezone *)0);

    now_time[0] = va1.tv_sec;
    now_time[1] = va1.tv_usec;
#else
    clock_t tnow;
    tnow = clock();
    iCLK_TCK = CLOCKS_PER_SEC;

    now_time[0] = tnow;
    now_time[1] = 0;
#endif

    if (count == 1)
    {
        first_time[0] = now_time[0];
        first_time[1] = now_time[1];

        count = 2;
    }

    now_time[0] -= first_time[0];
    now_time[1] -= first_time[1];

    GreenTime = now_time[0] + 0.000001 * now_time[1];
    GreenTime = GreenTime / iCLK_TCK;

    return GreenTime;
}

LIB_EXPORT double GetWallTimeForResidual()
{
    static int count = 1;

    double GreenTime;

    static double time0 = 0.0;

    if (count == 1)
    {
        time0 = (double)clock() / CLOCKS_PER_SEC;

        count = 2;
    }

    GreenTime = (double)clock() / CLOCKS_PER_SEC - time0;

    return GreenTime;
}

LIB_EXPORT int GetCPUTime(double *ptr_time)
{
    int ErrCode = 0;

    static clock_t first_time[2];

    static int count = 1;

    clock_t now_time[2];

    long iCLK_TCK;

    double time_array[2];

#ifndef WIN32
    struct tms tp;
    times(&tp);

    iCLK_TCK = sysconf(_SC_CLK_TCK);

    now_time[CPUTimeIndex] = tp.tms_utime;
    now_time[SYSTimeIndex] = tp.tms_stime;
#else
    iCLK_TCK = CLOCKS_PER_SEC;

    now_time[CPUTimeIndex] = clock();
    now_time[SYSTimeIndex] = 0;
#endif

    now_time[CPUTimeIndex] += now_time[SYSTimeIndex];

    if (count == 1)
    {
        first_time[CPUTimeIndex] = now_time[CPUTimeIndex];
        first_time[SYSTimeIndex] = now_time[SYSTimeIndex];

        count = 2;
    }

    now_time[CPUTimeIndex] -= first_time[CPUTimeIndex];
    now_time[SYSTimeIndex] -= first_time[SYSTimeIndex];

    time_array[CPUTimeIndex] = now_time[CPUTimeIndex] * 1.0 / iCLK_TCK;
    time_array[SYSTimeIndex] = now_time[SYSTimeIndex] * 1.0 / iCLK_TCK;

    if (ptr_time)
    {
        ptr_time[CPUTimeIndex] = time_array[CPUTimeIndex];
        ptr_time[SYSTimeIndex] = time_array[SYSTimeIndex];
    }

    return ErrCode;
}

LIB_EXPORT int WriStatTime()
{
    using namespace PHMPI;

    int ErrCode = 0;

    double cpu_time, sys_time, wall_time;

    double *time_array1 = new double[2]();

    ErrCode = GetCPUTime(time_array1);
    wall_time = GetWallTime();

    cpu_time = time_array1[CPUTimeIndex];
    sys_time = time_array1[SYSTimeIndex];

    delete [] time_array1;

    string cs;
    ostringstream oss;

    oss.setf(ios::fixed);
    oss.precision(3);

    oss << "\n    ----- Time statistics -----\n"
        << "          Wall time = " << wall_time << "\n"    //! Don't modify, ATP related.
        << "    ----- End time statistics -----\n";
    cs = oss.str();

    int root_proc           = 0;
    int myid                = GetCurrentProcessorID();
    int number_of_processor = GetNumberOfProcessor();

    double max_cpu_time  = -999.;
    double max_sys_time  = -999.;
    double max_wall_time = -999.;

    double min_cpu_time  = 999.;
    double min_sys_time  = 999.;
    double min_wall_time = 999.;

    double total_cpu_time  = 0.;
    double total_sys_time  = 0.;
    double total_wall_time = 0.;

    if (number_of_processor <= 1)
    {
        number_of_processor = 1;

        max_sys_time  = sys_time;
        max_cpu_time  = cpu_time;
        max_wall_time = wall_time;

        min_sys_time  = sys_time;
        min_cpu_time  = cpu_time;
        min_wall_time = wall_time;

        total_sys_time  = sys_time;
        total_cpu_time  = cpu_time;
        total_wall_time = wall_time;
    }
    else
    {
        int tag = 1;
        int proc;

        double *time_array = new double[3];

        if (myid == root_proc)
        {
            max_sys_time  = sys_time;
            max_cpu_time  = cpu_time;
            max_wall_time = wall_time;

            min_sys_time  = sys_time;
            min_cpu_time  = cpu_time;
            min_wall_time = wall_time;

            total_sys_time  = sys_time;
            total_cpu_time  = cpu_time;
            total_wall_time = wall_time;

            for (int i = 1; i < number_of_processor; ++ i)
            {
                proc = PH_Receive(time_array, 3*sizeof(double), tag);

                if (time_array[SYSTimeIndex]  > max_sys_time) max_sys_time  = time_array[SYSTimeIndex];
                if (time_array[CPUTimeIndex]  > max_cpu_time) max_cpu_time  = time_array[CPUTimeIndex];
                if (time_array[WallTimeIndex] > max_wall_time) max_wall_time = time_array[WallTimeIndex];

                if (time_array[SYSTimeIndex]  < min_sys_time) min_sys_time  = time_array[SYSTimeIndex];
                if (time_array[CPUTimeIndex]  < min_cpu_time) min_cpu_time  = time_array[CPUTimeIndex];
                if (time_array[WallTimeIndex] < min_wall_time) min_wall_time = time_array[WallTimeIndex];

                total_sys_time  += time_array[SYSTimeIndex];
                total_cpu_time  += time_array[CPUTimeIndex];
                total_wall_time += time_array[WallTimeIndex];
            }
        }
        else
        {
            time_array[SYSTimeIndex]  = sys_time;
            time_array[CPUTimeIndex]  = cpu_time;
            time_array[WallTimeIndex] = wall_time;

            PH_Send(time_array, 3 * sizeof(double), root_proc, tag);
        }

        delete [] time_array;
    }

    if (myid == root_proc)
    {
        cout << cs << endl;
    }

    (void)oss.flags();

    WriteLogFile(cs);

    return ErrCode;
}

LIB_EXPORT int WriStatTime(const char *title, double *last_time)
{
    int ErrCode = 0;

    double cpu_time, sys_time, wall_time;

    double *time_array = new double[2]();

    ErrCode   = GetCPUTime(time_array);
    wall_time = GetWallTime();

    cpu_time = time_array[CPUTimeIndex];
    sys_time = time_array[SYSTimeIndex];

    if (last_time)
    {
        time_array[0]           = sys_time;
        sys_time               -= last_time[SYSTimeIndex];
        last_time[SYSTimeIndex] = time_array[0];

        time_array[0]           = cpu_time;
        cpu_time               -= last_time[CPUTimeIndex];
        last_time[CPUTimeIndex] = time_array[0];

        time_array[0]            = wall_time;
        wall_time               -= last_time[WallTimeIndex];
        last_time[WallTimeIndex] = time_array[0];
    }

    delete [] time_array;

    string cs;
    ostringstream oss;

    oss.setf(ios::fixed);
    oss.precision(3);
    oss << title << " Time Statistics: SYS = " << sys_time
        << ", CPU=" << cpu_time << ", Wall = " << wall_time;
    cs = oss.str();
    (void)oss.flags();

    WriteLogFile(cs);

    return ErrCode;
}

LIB_EXPORT void StartTime(double *last_time)
{
    double *ptr_time = new double[2];

    TIME_SPACE::GetCPUTime(ptr_time);

    last_time[TIME_SPACE::WallTimeIndex] = TIME_SPACE::GetWallTime();
    last_time[TIME_SPACE::CPUTimeIndex]  = ptr_time[TIME_SPACE::CPUTimeIndex];
    last_time[TIME_SPACE::SYSTimeIndex]  = ptr_time[TIME_SPACE::SYSTimeIndex];

    delete [] ptr_time;
}

LIB_EXPORT void GetFileModifiedTime(string &fileName, RDouble &TimeFile)
{
    string filename = fileName;
    struct stat result;
    if (stat(filename.c_str(), &result) == 0 )
    {
        RDouble mod_time = static_cast<RDouble>(result.st_mtime);
        TimeFile = mod_time;
    }
}

}
#ifdef WIN32
LIB_EXPORT TimeSpan::TimeSpan()
{
    last_time = GetWallTime();
    now_time = last_time;
}
#else
TimeSpan::TimeSpan()
{
    gettimeofday(&last_time, NULL);
}
#endif
}