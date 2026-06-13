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
//! @file      TK_Parse.h
//! @brief     Parsing parameters from files to database.
//! @author    Zhang Jian.

#pragma once
#include <fstream>
#include "LIB_Macro.h"
#include <sstream>
using namespace std;

namespace PHSPACE
{
class Data_Param;

class TK_Parse
{
public:
    //! Read and parse data from file, then update to global database.
    //! @param[in] file        file stream, the file stores parameter data, e.g. cfd_para.hypara.
    //! @param[in] database    database that stores parameters.
    LIB_EXPORT static void ReadBasicData(fstream &file, Data_Param *database);

private:
    TK_Parse();
};

class GSStringToLowerCase
{
public:
    char operator()(char val)
    {
        return static_cast<char>(tolower(val));
    }
};

class GSStringToUpperCase
{
public:
    char operator()(char val)
    {
        return static_cast<char>(toupper(val));
    }
};

class ToLower
{
public:
    char operator()(char val)
    {
        return static_cast<char>(tolower(val));
    }
};

template < typename T >
inline T StringToDigit(const string &str, ios_base & (*f)(ios_base &))
{
    T value;
    istringstream stream(str);
    stream >> f >> value;
    return value;
}

//! Converts string to streamable value, and returns true on success and false otherwise.
template < typename T >
inline bool from_string(T &Value, const std::string &str, std::ios_base & (*f)(std::ios_base&))
{
    std::istringstream stream(str);
    stream >> f >> Value;
    return (!stream.fail()) && stream.get() == std::istringstream::traits_type::eof();
}

//! Converts value to string.
template < typename T >
inline std::string to_string(T t, std::ios_base & (*f)(std::ios_base&))
{
    std::ostringstream oss;
    oss << f << t;
    return oss.str();
}

template < typename T1, typename T2 >
string AggregateGeneralString(const T1 &v1, const T2 &v2)
{
    ostringstream oss;
    oss << v1 << v2;
    return oss.str();
}

template < typename T1, typename T2, typename T3 >
string AggregateGeneralString(const T1 &v1, const T2 &v2, const T3 &v3)
{
    ostringstream oss;
    oss << v1 << v2 << v3;
    return oss.str();
}

template < typename T1, typename T2, typename T3, typename T4 >
string AggregateGeneralString(const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4)
{
    ostringstream oss;
    oss << v1 << v2 << v3 << v4;
    return oss.str();
}

template < typename T1, typename T2, typename T3, typename T4, typename T5 >
string AggregateGeneralString(const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5)
{
    ostringstream oss;
    oss << v1 << v2 << v3 << v4 << v5;
    return oss.str();
}

template < typename T1, typename T2, typename T3, typename T4, typename T5, typename T6 >
string AggregateGeneralString(const T1 &v1, const T2 &v2, const T3 &v3, const T4 &v4, const T5 &v5, const T6 &v6)
{
    ostringstream oss;
    oss << v1 << v2 << v3 << v4 << v5 << v6;
    return oss.str();
}

void ToLowerCase(string &word);
void ToUpperCase(string &word);
void TrimBlanks(string &source);
void SkipLines(fstream &file, int nline);
void SkipWords(string &line, int nword);
bool IsEmptyLine(const string &line);
void ReadNextNonEmptyLine(fstream &file, string &line);
string FindNextWord(const string &source, string &word, const string &separator);
bool IsArrayParameter(const string &lineOfName);
void ReadNextLine(fstream &file, string &line);
bool IsCommentLine(const string &line);
string FindNextWord(string &source, string &word, const string &separator);
string FindNextWord(string &source, const string &separator);
void GetNameExt(const string &fullname, string &fname, string &fext, const string &sep);
void ReadNewLine(fstream &file, string &line);
string Dec2Hex(int data_in);

}