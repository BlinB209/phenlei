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
//! @file      TK_Exit.h
//! @brief     Exit or abort the problem because of different reasons.
//! @author    Bell, Zhang Jian.

#pragma once
#include "PHMpi.h"
using namespace std;

namespace PHSPACE
{

//! @brief TK_Exit is a class that integrate several methods that terminate the program.
class TK_Exit
{
public:
    //! Normal Exit PHengLEI.
    //! @example: TK_Exit::ExitPHengLEI().
    LIB_EXPORT static void ExitPHengLEI();

    //! Abnormal Exit PHengLEI and dumping error information to Log file.
    //! @param[in] errorMessage    string type error information need to dump.
    //! @param[in] onlyServer      default false. If true, whether or not to abort program requires 
    //!                            to judge if the current process is server. Only server process call abort.
    //!                            Warning!!! Do not set true unless you ensure server process meets an exception.
    //! @example: TK_Exit::ExceptionExit("Front edge is not exist !");
    LIB_EXPORT static void ExceptionExit(const string &errorMessage, bool onlyServer = false);

    //! Abnormal Exit PHengLEI and dumping error information to Log file.
    //! @param[in] errorMessage    ostringstream type error information need to dump.
    //! @param[in] onlyServer      Default false. If true, whether or not to abort program requires 
    //!                            to judge if the current process is server. Only server process call abort.
    //!                            Warning!!! Do not set true unless you ensure server process meet an exception.
    //! @example: ostringstream oss;
    //!           oss << "Out of Range: pos " << pos << endl;
    //!           TK_Exit::ExceptionExit(oss);
    LIB_EXPORT static void ExceptionExit(const ostringstream &errorMessage, bool onlyServer = false);

    //! Print debug information to log file or screen.
    //! This function is used by developers to debug.
    //! Information include filename, line number, compile date and time.
    //! @param[in] debuginfo    what developers want to say.
    //! @example: TK_Exit::PrintDebugInfoExit("Error: numberOfTotalCompositeFaces != 0");
    LIB_EXPORT static void PrintDebugInfoExit(const string &debuginfo);

    //! Print debug information to log file or screen.
    //! This function is used by developers to debug.
    //! Information include filename, line number, compile date and time.
    //! @param[in] debuginfo    what developers want to say.
    //! @example: TK_Exit::PrintDebugInfoExit("Error: numberOfTotalCompositeFaces != 0");
    LIB_EXPORT static void PrintDebugInfoExit(const ostringstream &debuginfo);

    //! Print "Could not open [filename]" and exit when errors occurs during opening a file.
    //! @param[in] filename    file name.
    LIB_EXPORT static void FileOpenErrorExit(const string &filename);

    //! An unexpected situation occurs because of a variable's value is not right.
    //! It will cause unexpected error so that program will dump an information like 
    //! "Error: this situation has not been considered, for [variableName] = [value]" 
    //! and abort.
    //! @param[in] variableName    variable's name.
    //! @param[in] value           variable's value.
    //! @example: TK_Exit::UnexpectedVarValue("force_component", 0);
    template < typename T >
    static void UnexpectedVarValue(const string &variableName, const T &value);

private:
    //! Constructor private.
    TK_Exit();

    //! Called by PrintDebugInfoExit(string debuginfo).
    static void PrintDebugInfoExit(const string &stopInformation, const string &fileName, const int &fileLine, const string &dateName, const string &timeName);
};

#include "TK_Exit.hxx"
}