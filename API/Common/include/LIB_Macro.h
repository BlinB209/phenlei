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
//! @file      LIB_Macro.h
//! @brief     This file is part of PHengLEI Technology software library.
//!            This file is intended to be the first file #included to any
//!            PHengLEI source. It defines platform-specific pre-processor 
//!            macros necessary for correct compilation of PHengLEI code
//!            It's used to define the DLL/Static library.
//! @author    Bell modified it from OCC.

#pragma once
#include <iostream>
using namespace std;

#if defined(__cplusplus) && (__cplusplus >= 201100L)
  //! Part of C++11 standard.
  #define Standard_OVERRIDE override
#elif defined(_MSC_VER) && (_MSC_VER >= 1700)
  //! MSVC extension since VS2012.
  #define Standard_OVERRIDE override
#else
  #define Standard_OVERRIDE
#endif

//======================================================
// Windows-specific definitions
//======================================================
# if defined(_WIN32) && !defined(HAVE_NO_DLL)

#  ifndef LIB_EXPORT
#   define LIB_EXPORT __declspec(dllexport)
//! For global variables :
#   define LIB_EXPORTEXTERN __declspec(dllexport) extern
#   define LIB_EXPORTEXTERNC extern "C" __declspec(dllexport)
#  endif  /* LIB_EXPORT */

#  ifndef LIB_IMPORT
#   define LIB_IMPORT __declspec(dllimport) extern
#   define LIB_IMPORTC extern "C" __declspec(dllimport)
#  endif  /* LIB_IMPORT */

# else  /* WNT */

//======================================================
// UNIX definitions
//======================================================

#  ifndef LIB_EXPORT
#   define LIB_EXPORT
//! For global variables :
#   define LIB_EXPORTEXTERN extern
#   define LIB_EXPORTEXTERNC extern "C"
#  endif  /* Standard_EXPORT */

#  ifndef LIB_IMPORT
#   define LIB_IMPORT extern
#   define LIB_IMPORTC extern "C"
#  endif  /* Standard_IMPORT */

//! Compatibility with old SUN compilers.

//! This preprocessor directive is a kludge to get around
//! a bug in the Sun Workshop 5.0 compiler, it keeps the
//! /usr/include/memory.h file from being #included
//! with an incompatible extern "C" definition of memchr
//! October 18, 2000  <rboehne@ricardo-us.com>.
#if __SUNPRO_CC_COMPAT == 5
#define _MEMORY_H
#endif

# endif  /* WNT */

//======================================================
// Other
//======================================================

# ifndef __Standard_API
//#  ifdef WNT
#   if !defined(_WIN32) || defined(__Standard_DLL) || defined(__FSD_DLL) || defined(__MMgt_DLL) || defined(__OSD_DLL) || defined(__Plugin_DLL) || defined(__Quantity_DLL) || defined(__Resource_DLL) || defined(__SortTools_DLL) || defined(__StdFail_DLL) || defined(__Storage_DLL) || defined(__TColStd_DLL) || defined(__TCollection_DLL) || defined(__TShort_DLL) || defined(__Units_DLL) || defined(__UnitsAPI_DLL) || defined(__Dico_DLL)
#    define __Standard_API Standard_EXPORT
#    define __Standard_APIEXTERN Standard_EXPORTEXTERN
#   else
#    define __Standard_API Standard_IMPORT
#    define __Standard_APIEXTERN Standard_IMPORT
#   endif  // __Standard_DLL
//#  else
//#   define __Standard_API
//#  endif  // WNT
# endif  // __Standard_API
 
