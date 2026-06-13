#include "TK_Log.h"
#include "TK_Exit.h"
#include <cmath>
#ifdef WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

using namespace std;

namespace PHSPACE
{

void OpenLogFile(int file_id, fstream &file)
{
    PHMakeDir("results");

    static int if_rewrite = 0;

    std::ostringstream oss;
    oss << "results/log" << file_id << ".log";

    if (if_rewrite == 0)
    {
        file.open(oss.str().c_str(), ios_base::out|ios_base::trunc);
    }
    else
    {
        file.open(oss.str().c_str(), ios_base::out|ios_base::app);
    }

    if_rewrite = 1;

    if (!file)
    {
        TK_Exit::FileOpenErrorExit(oss.str());
    }
}

void CloseLogFile(fstream &file)
{
    file.close();
    file.clear();
}

void WriteLogFile(const string &data)
{
    WriteLogFile(data, false);
}

void WriteLogFile(const string &data, bool writeEachProcessor)
{
    using namespace PHMPI;
    if (writeEachProcessor)
    {
        fstream file;
        OpenLogFile(GetCurrentProcessorID(), file);
        file << data << "\n";
        CloseLogFile(file);
    }
    else
    {
#ifndef LOG_DEBUG
        //! Only server write log.
        if (GetCurrentProcessorID() == server)
#endif
        {
            fstream file;
            OpenLogFile(GetCurrentProcessorID(), file);
            file << data << "\n";
            CloseLogFile(file);
        }
    }
}

void WriteLogFile(const ostringstream &oss, bool writeEachProcessor)
{
    if (writeEachProcessor == true)
    {
        fstream file;
        OpenLogFile(PHMPI::GetCurrentProcessorID(), file);
        file << oss.str() << "\n";
        CloseLogFile(file);
    }
    else
    {
        WriteLogFile(oss.str());
    }
}

void OpenWarningLogFile(int file_id, fstream &file)
{
    PHMakeDir("results");
    static int if_rewrite = 0;

    std::ostringstream oss;
    oss << "results/warning" << file_id << ".log";

    if (if_rewrite == 0)
    {
        file.open(oss.str().c_str(), ios_base::out|ios_base::trunc);
    }
    else
    {
        file.open(oss.str().c_str(), ios_base::out|ios_base::app);
    }

    if_rewrite = 1;

    if (!file)
    {
        TK_Exit::FileOpenErrorExit(oss.str());
    }
}

void WriteWarningLogFile(const string &warningInfo)
{
    if (PHMPI::GetCurrentProcessorID() == PHMPI::server)
    {
        fstream file;
        OpenWarningLogFile(PHMPI::GetCurrentProcessorID(), file);
        file << warningInfo << "\n";
        CloseLogFile(file);
    }
}

void PrintToWindow(const string &information)
{
    if (PHMPI::GetCurrentProcessorID() != PHMPI::GetServerProcessorID()) return;

    cout << information;
}

void PrintToWindow(const ostringstream &information)
{
    PrintToWindow(information.str());
}

void PHMakeDir(const string &dir_name)
{
    using namespace PHMPI;
    //if (GetCurrentProcessorID() == server)
    {
        int flag;

        #ifdef WIN32
            flag = _mkdir(dir_name.c_str());
        #else
            flag = mkdir(dir_name.c_str(), S_IRWXU);
        #endif
            if (flag == 0)
            {
                cout << dir_name<< " directory has been created successfully !\n";
            }
    }
}

void ProgressMonitorToLogFile(const int &iNow, const int &nTotal, const string &monitorName)
{
    ostringstream oss;
    oss.unsetf(ios::scientific);
    oss.setf(ios::fixed);
    oss.precision(0);
    oss << "     " << floor((iNow * 1.0) / nTotal * 100 + 0.5) << "% " << monitorName << " completed ...";
    WriteLogFile(oss);
}

void ProgressMonitorToWindows(const int &iNow, const int &nTotal, const string &monitorName)
{
    ostringstream oss;
    oss.unsetf(ios::scientific);
    oss.setf(ios::fixed);
    oss.precision(0);
    oss << "     " << floor((iNow * 1.0) / nTotal * 100 + 0.5) << "% " << monitorName << " completed ..." << endl;
    PrintToWindow(oss);
}

}