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
//! @file      TK_Warning.h
//! @brief     Write or print different warning information.
//! @author    Zhang Jian.
//#include <vld.h>
#pragma once
#include "TK_Log.h"
#include "LIB_Macro.h"
using namespace std;

namespace PHSPACE
{
class TK_Warning
{
public:
    //! Print warning information because of function will be deprecated.
    //! @param[in] oldFunctionName    old function's name which will be deprecated.
    LIB_EXPORT static void FunctionDeprecation(const string &oldFunctionName);

    //! Print warning information because of function will be deprecated.
    //! @param[in] oldFunctionName    old function's name which will be deprecated.
    //! @param[in] newFunctionName    new function's name.
    LIB_EXPORT static void FunctionDeprecation(const string &oldFunctionName, const string &newFunctionName);

    //! Print warning information to the window.
    //! @param[in] information    the warning information will be print to the window.
    LIB_EXPORT static void Warning(const string &information);

    //! Print warning information because of variable's value is not recommended.
    template < typename T >
    static void UnexpectedVarValue(const string &variableName, const T &value);

private:
    TK_Warning();
};

#include "TK_Warning.hxx"
}