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
//! @file      IO_FileReader.h
//! @brief     file reader.
//! @author    Guo Yongheng.

#pragma once
#include "LIB_Macro.h"
#include "TK_Parse.h"

namespace PHSPACE
{
class IO_FileReader
{
private:
    std::string *line, *separator;
    fstream *file;
    int setfileFlag;

    std::string fileName;
    ios_base::openmode fileOpenMode;
    streamsize filePosition;
public:
    LIB_EXPORT  IO_FileReader();
    LIB_EXPORT ~IO_FileReader();

public:
    //! Open the file.
    //! @param[in] fileName        file name.
    //! @param[in] fileOpenMode    open mode.
    void OpenFile(const string &fileName, const ios_base::openmode &fileOpenMode);

    //! Close the file.
    void CloseFile();

    //!
    void MarkCurrentFilePosition();

    //!
    void MoveToPreviousFilePosition();

public:
    void SetDefaultFile(std::fstream *defaultFileIn);
    void SetDefaultSeparator(const std::string &separatorIn) { * this->separator = separatorIn; }

    std::string  * GetDefaultLine() { return line; }
    std::fstream * GetDefaultFile() { return file; }
    std::string  * GetDefaultSeparator() { return separator; }
public:
    //!
    void SetLineContent(const std::string &lineContent) { * this->line = lineContent; }

    //!
    void ShiftLineContent(int numberOfChars) { * this->line = (* this->line).substr(numberOfChars); }

    bool ReadNextMeaningfulLine();
    bool ReachTheEndOfFile();
    bool IsArrayParameter();

public:
    //! Skip lines.
    //! @param[in] numberOfLinesToSkip    the number of lines to skip.
    void SkipLines(int numberOfLinesToSkip);

    //!
    void ReadNextNonEmptyLine();

    //!
    void DumpLineContentToScreen();

    //!
    std::string ReadNextWord();

    //!
    std::string ReadNextWord(const std::string &separator);

    //!
    std::string ReadNextWordToLowerCase();

    //!
    std::string ReadNextWordToLowerCase(const std::string &separator);
public:
    //!
    void SkipReadSymbol(const string &stringSymbol);

    //!
    void SkipReadWholeBlock();

public:
    template < typename T >
    friend inline IO_FileReader & operator >> (IO_FileReader &ioStruct, T &value)
    {
        fstream &file = * ioStruct.GetDefaultFile();
        file >> value;
        return ioStruct;
    }

    template < typename T >
    T ReadNextDigit(ios_base & (* f)(ios_base &) = & std::dec)
    {
        std::string word = PHSPACE::FindNextWord(*line, *separator);
        T value = PHSPACE::StringToDigit< T >(word, f);
        return value;
    }
};

}