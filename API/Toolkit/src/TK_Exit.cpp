#include <stdarg.h>
#include "TK_Exit.h"
#include "TK_Time.h"

namespace PHSPACE
{
LIB_EXPORT void TK_Exit::ExitPHengLEI()
{
    PHMPI::FreeBlockData();

    TIME_SPACE::WriStatTime();
    PH_Finalize();
}

LIB_EXPORT void TK_Exit::ExceptionExit(const string &errorMessage, bool onlyServer)
{
    if (onlyServer)
    {
        int myid = PHMPI::GetCurrentProcessorID();
        int server = PHMPI::GetServerProcessorID();
        if (myid == server)
        {
            PHMPI::FreeBlockData();
            PrintToWindow(errorMessage);
            WriteLogFile(errorMessage, true);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    else
    {
        PHMPI::FreeBlockData();
        cout << errorMessage << endl;
        WriteLogFile(errorMessage, true);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}

LIB_EXPORT void TK_Exit::ExceptionExit(const ostringstream &errorMessage, bool onlyServer)
{
    ExceptionExit(errorMessage.str(), onlyServer);
}

/*LIB_EXPORT void TK_Exit::ExceptionExit(const char *format, ...)
{
    string var_str;
    va_list argList;
    va_start(argList, format);
    int len = vsnprintf(NULL, 0 ,format, argList);
    if (len > 0)
    {
        vector<char> buf(len+1);
        vsprintf(&buf.front(), format, argList);
        var_str.assign(buf.begin(), buf.end()-1);
    }
    va_end(argList);
    ExceptionExit(var_str);
}*/

LIB_EXPORT void TK_Exit::PrintDebugInfoExit(const string &debuginfo)
{
    PrintDebugInfoExit(debuginfo, __FILE__, __LINE__, __DATE__, __TIME__);
}

LIB_EXPORT void TK_Exit::PrintDebugInfoExit(const ostringstream &debuginfo)
{
    PrintDebugInfoExit(debuginfo.str());
}

LIB_EXPORT void TK_Exit::FileOpenErrorExit(const string &filename)
{
    ostringstream oss;
    oss << "Could not open " << filename << endl;
    ExceptionExit(oss.str());
}

void TK_Exit::PrintDebugInfoExit(const string &stopInformation, const string &fileName, const int &fileLine, const string &dateName, const string &timeName)
{
    ostringstream oss;
    oss << "\n";
    oss << "++++++++++++++++++  StopProgram Information  +++++++++++++++++++++++++++\n";
    oss <<  stopInformation << "\n";
    oss << " The stop filename is: " << fileName << "\n";
    oss << " at line " << fileLine << "\n";
    oss << " Compiled On line " << dateName << " at " << timeName << "\n";
    oss << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    ExceptionExit(oss.str());
}

}