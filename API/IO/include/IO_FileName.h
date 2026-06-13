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
//! @file      IO_FileName.h
//! @brief     file name.
//! @author    Baka.

#pragma once
#include <sstream>
using namespace std;

namespace PHSPACE
{

//! Return the file name of the current processor.
//! The processor ID would be add to the input filename.
void GetFileNameofThisProcessor(string &filename);

//! Add digits to file name.
string AddDigitToFileName(const string &fname, const string &fext, const string &sep, int num);

string GetDefaultDirectory();

string AbsoluteFilename(const string &relativeFilename);

//! Get the file name's extension.
void GetFileNameExtension(const string &fullName, string &mainName, string &extensionName, const string &fileNameSeparator);

//! Change the file name's extension.
string ChangeExtensionOfFileName(const string &fileName, const string &newExtension);

//! Change the file's main name.
string ChangeFileMainName(const string &fileName, const string &newMainName);

int ExtractNumberOfPartFromFileName(const string &source);

template < typename T >
string AddSymbolToFileName(const string &fileName, const T &symbol)
{
    if (fileName == "") return "";

    string mainName, extensionName;
    GetFileNameExtension(fileName, mainName, extensionName, ".");

    ostringstream oss;
    oss << mainName << symbol << "." << extensionName;
    return oss.str();
}

template < typename T1, typename T2 >
string AddSymbolToFileName(const string &fileName, const T1 &v1, const T2 &v2)
{
    if (fileName == "") return "";

    string mainName, extensionName;
    GetFileNameExtension(fileName, mainName, extensionName, ".");

    ostringstream oss;
    oss << mainName << v1 << v2 << "." << extensionName;
    return oss.str();
}

template < typename T1, typename T2, typename T3 >
string AddSymbolToFileName(const string &fileName, const T1 &v1, const T2 &v2, const T3 &v3)
{
    if (fileName == "") return "";

    string mainName, extensionName;
    GetFileNameExtension(fileName, mainName, extensionName, ".");

    ostringstream oss;
    oss << mainName << v1 << v2 << v3 << "." << extensionName;
    return oss.str();
}

string GetProcessPath();

string GetCurrentWorkDir();

}