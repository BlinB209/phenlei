#include "TK_Warning.h"
#include "TK_Log.h"

namespace PHSPACE
{

TK_Warning::TK_Warning() {}

LIB_EXPORT void TK_Warning::FunctionDeprecation(const string &oldFunctionName)
{
    ostringstream oss;
    oss << "Warning!!! " << oldFunctionName << " will be deprecated soon, please abandon it!" << endl;
    PrintToWindow(oss.str());
    WriteWarningLogFile(oss.str());
}

LIB_EXPORT void TK_Warning::FunctionDeprecation(const string &oldFunctionName, const string &newFunctionName)
{
    ostringstream oss;
    oss << "Warning!!! " << oldFunctionName << " will be deprecated soon, please replace it with " << newFunctionName << "." << endl;
    PrintToWindow(oss.str());
    WriteWarningLogFile(oss.str());
}

LIB_EXPORT void TK_Warning::Warning(const string &information)
{
    ostringstream oss;
    oss << "Warning: " << information << " !!!" << endl;
    PrintToWindow(oss.str());
    WriteWarningLogFile(oss.str());
}

}