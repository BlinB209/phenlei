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
//! @file      TK_Log.h
//! @brief     some functions about Log file.
//! @author    He Xin.

#pragma once
#include "PHMpi.h"
using namespace std;
//#define LOG_DEBUG
#ifdef LOG_DEBUG
#include "Math_BasisFunction.h"
#endif

namespace PHSPACE
{
//! Open the log file.
//! @param[in] file_id    file ID.
//! @param[in] file       file stream flow.
void OpenLogFile(int file_id, fstream &file);

//! Close the log file.
//! @param[in] file    file stream flow.
void CloseLogFile(fstream &file);

//! Output information to log file.
//! @param[in] data    string data.
void WriteLogFile(const string &data);

//! Output information to log file.
//! @param[in] data    string data.
//! @param[in] writeEachProcessor    whether write each processor.
void WriteLogFile(const string &data, bool writeEachProcessor);

//! Output information to log file.
//! @param[in] oss    information string stream.
//! @param[in] writeEachProcessor    whether write each processor.
void WriteLogFile(const ostringstream &oss, bool writeEachProcessor = false);

//! Open warning information to log file.
//! @param[in] file_id    file ID.
//! @param[in] file       file stream flow.
void OpenWarningLogFile(int file_id, fstream &file);

//! Output warning information to log file.
//! @param[in] warningInfo    warning information.
void WriteWarningLogFile(const string &warningInfo);

//! Print information to window.
//! @param[in] information    information.
void PrintToWindow(const string &information);

//! Print information to window.
//! @param[in] information    information string stream.
void PrintToWindow(const ostringstream &information);

//! Create file directory.
//! @param[in] dir_name    file directory name.
void PHMakeDir(const string &dir_name);

template < typename T >
void WriteLogFile(const string &name, const T &data, bool writeEachProcessor = false)
{
    using namespace PHMPI;
#ifdef LOG_DEBUG
    int interval = GetProgressInterval((int)GetNumberOfProcessor(), 64);
    if (writeEachProcessor || GetCurrentProcessorID() % interval == 0 || GetCurrentProcessorID() == (GetNumberOfProcessor() - 1))
#else
    if (GetCurrentProcessorID() == server || writeEachProcessor)
#endif
    {
        fstream file;
        int myid = PHMPI::GetCurrentProcessorID();
        OpenLogFile(myid, file);
        file << name << data << "\n";
        CloseLogFile(file);
    }
}

template < typename T1, typename T2 >
void WriteLogFile(const string &name, const T1 &data1, const T2 &data2, bool writeEachProcessor = false)
{
    using namespace PHMPI;
#ifdef LOG_DEBUG
    int interval = GetProgressInterval((int)GetNumberOfProcessor(), 64);
    if (writeEachProcessor || GetCurrentProcessorID() % interval == 0 || GetCurrentProcessorID() == (GetNumberOfProcessor() - 1))
#else
    if (GetCurrentProcessorID() == server || writeEachProcessor)
#endif
    {
        fstream file;
        int myid = PHMPI::GetCurrentProcessorID();
        OpenLogFile(myid, file);
        file << name << data1 << " " << data2 << "\n";
        CloseLogFile(file);
    }
}

template < typename T1, typename T2, typename T3 >
void WriteLogFile(const string &name, const T1 &data1, const T2 &data2, const T3 &data3, bool writeEachProcessor = false)
{
    using namespace PHMPI;
#ifdef LOG_DEBUG
    int interval = GetProgressInterval((int)GetNumberOfProcessor(), 64);
    if (writeEachProcessor || GetCurrentProcessorID() % interval == 0 || GetCurrentProcessorID() == (GetNumberOfProcessor() - 1))
#else
    if (GetCurrentProcessorID() == server || writeEachProcessor)
#endif
    {
        fstream file;
        int myid = PHMPI::GetCurrentProcessorID();
        OpenLogFile(myid, file);
        file << name << data1 << " " << data2  << " " << data3 << "\n";
        CloseLogFile(file);
    }
}

template < typename T >
void PrintToWindow(const string &name, const T &data)
{
    if (PHMPI::GetCurrentProcessorID() != 0) return;

    ostringstream oss;
    oss << name << " " << data;
    PrintToWindow(oss.str());
}

template < typename T1, typename T2 >
void PrintToWindow(const string &name, const T1 &data1, const T2 &data2)
{
    if (PHMPI::GetCurrentProcessorID() != 0) return;

    ostringstream oss;
    oss << name << " " << data1 << " " << data2;
    PrintToWindow(oss.str());
}

template < typename T1, typename T2, typename T3 >
void PrintToWindow(const string &name, const T1 &data1, const T2 &data2, const T3 &data3)
{
    if (PHMPI::GetCurrentProcessorID() != 0) return;

    ostringstream oss;
    oss << name << " " << data1 << " " << data2 << " " << data3;
    PrintToWindow(oss.str());
}

template < typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9 >
void PrintToWindow(const string &name, const T1 &data1, const T2 &data2, const T3 &data3, const T4 &data4, const T5 &data5, const T6 &data6, const T7 &data7, const T8 &data8, const T9 &data9)
{
    if (PHMPI::GetCurrentProcessorID() != 0) return;

    ostringstream oss;
    oss << name << " " << data1 << " " << data2 << " " << data3 << " " << data4 << " " << data5 << " " << data6 << " " << data7 << " " << data8 << " " << data9;
    PrintToWindow(oss.str());
}

template < typename T1, typename T2, typename T3, typename T4>
void PrintToWindow(const string &name, const T1 &data1, const T2 &data2, const T3 &data3,  const T4 &data4)
{
    if (PHMPI::GetCurrentProcessorID() != 0) return;

    ostringstream oss;
    oss << name << " " << data1 << " " << data2 << " " << data3 << " " << data4 ;
    PrintToWindow(oss.str());
}

template < typename T1, typename T2, typename T3, typename T4, typename T5 >
void PrintToWindow(const string &name, const T1 &data1, const T2 &data2, const T3 &data3, const T4 &data4, const T5 &data5)
{
    if (PHMPI::GetCurrentProcessorID() != 0) return;

    ostringstream oss;
    oss << name << " " << data1 << " " << data2 << " " << data3 << " " << data4 << " " << data5;
    PrintToWindow(oss.str());
}

//! Monitor the progress to log files.
void ProgressMonitorToLogFile(const int &iNow, const int &nTotal, const string &monitorName);

//! Monitor the progress by printing to windows.
void ProgressMonitorToWindows(const int &iNow, const int &nTotal, const string &monitorName);

}