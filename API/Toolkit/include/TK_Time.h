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
//! @file      TK_Time.h
//! @brief     Time computing tools.
//! @author    Xu Qingxin, Bell.

#pragma once
#include "TK_Log.h"
#include <sys/types.h>
#include <sys/stat.h>
#ifndef WIN32
#include <sys/time.h>
#include <unistd.h>
#endif

#ifdef WIN32
#define stat _stat
#endif
using namespace std;

namespace PHSPACE
{
namespace TIME_SPACE
{

//! Set time type index value.
enum TimeTypeIndex
{
    CPUTimeIndex  = 0,
    SYSTimeIndex  = 1,
    WallTimeIndex = 2
};

//! Get Wall Time from first call.
LIB_EXPORT double GetWallTime();

//! Get Wall Time when dumping residuals.
LIB_EXPORT double GetWallTimeForResidual();

//! Get CPU Time from first call.
LIB_EXPORT int GetCPUTime(double *ptr_time);

//! Collect CPU time, system time, wall time, and write to log file.
LIB_EXPORT int WriStatTime();

//! Collect CPU time, system time, wall time from last_time, and write to log file.
//! @param[in] title    testing tag.
//! @param[in] last_time
LIB_EXPORT int WriStatTime(const char *title, double *last_time);

//! Write starting CPU time, system time, wall time to pri_time.
//! @param[in] pri_time    must be a vector of three element.
LIB_EXPORT void StartTime(double *ptr_time);

//! Get the last modification time of files.
LIB_EXPORT void GetFileModifiedTime(string &fileName, RDouble &TimeFile);
}

#ifdef WIN32
//! @brief The main interface class of time computing functions.
//! Windows.
class TimeSpan
{
public:
    LIB_EXPORT TimeSpan();

    LIB_EXPORT ~TimeSpan(){};

    //! Reset the starting testing time.
    void ResetTime();

    //! Get the elapsed time from ResetTime.
    //! @return (nowtime - lasttime).
    double GetTimeSpan();

    //! Print the elapsed time to screen.
    //! @param[in] title    testing tag.
    //! -# default is null("").
    void ShowTimeSpan(const string &title = "");

    //! Write the elapsed time to log file.
    //! @param[in] title    testing tag.
    //! -# default is null("").
    void ShowTimeSpanToLogFile(const string &title = "");

private:
    double last_time, now_time;

private:
    double GetWallTime();
    double GetSecondTime(const double &time_elapsed);
};

#else
//! Linux/Unix.
class TimeSpan
{
public:
    TimeSpan();

    ~TimeSpan(){};

    //! Reset the stating testing time.
    void ResetTime();

    //! Get the elapsed time from ResetTime.
    //! @return (nowtime - lasttime).
    double GetTimeSpan();

    //! Print the elapsed time to screen.
    //! @param[in] title    testing tag.
    //! -# default is null("").
    void ShowTimeSpan(const string &title = "");

    //! Write the elapsed time to log file.
    //! @param[in] title    testing tag.
    //! -# default is null("").
    void ShowTimeSpanToLogFile(const string &title = "");

private:
    struct timeval last_time, now_time;
    //! tv_sec : second.
    //! tv_usec: microsecond.
};
#endif

#include "TK_Time.hxx"
}