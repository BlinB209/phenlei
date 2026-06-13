#include "IO_FileReader.h"
#include "TK_Parse.h"
#include "PHIO.h"

namespace PHSPACE
{
LIB_EXPORT IO_FileReader::IO_FileReader()
{
    line        = new std::string;
    separator   = new std::string;
    file        = new fstream;
    setfileFlag = 0;
    string keyWordSeparator = " =\r\n\t#$,;\"";
    this->SetDefaultSeparator(keyWordSeparator);
}

LIB_EXPORT IO_FileReader::~IO_FileReader()
{
    delete line;
    delete separator;
    if (setfileFlag == 0)
    {
        delete file;
    }
}

void IO_FileReader::OpenFile(const string &fileName, const ios_base::openmode &fileOpenMode)
{
    this->fileName     = fileName;
    this->fileOpenMode = fileOpenMode;
    PHSPACE::OpenFile(*file, fileName, fileOpenMode);
}

void IO_FileReader::CloseFile()
{
    PHSPACE::CloseFile(*file);
}

void IO_FileReader::MarkCurrentFilePosition()
{
    filePosition = file->tellp();
}

void IO_FileReader::MoveToPreviousFilePosition()
{
    file->seekp(filePosition);
}

void IO_FileReader::SetDefaultFile(std::fstream *defaultFileIn)
{
    if (setfileFlag == 0)
    {
        delete file;
    }
    this->file = defaultFileIn;
    setfileFlag = 1;
}

bool IO_FileReader::ReadNextMeaningfulLine()
{
    while (!this->ReachTheEndOfFile())
    {
        PHSPACE::ReadNextLine(*file, *line);

        if (PHSPACE::IsEmptyLine(*line) || PHSPACE::IsCommentLine(*line))
        {
            continue;
        }
        return true;
    }
    return false;
}

bool IO_FileReader::ReachTheEndOfFile()
{
    if ((*file).eof())
    {
        return true;
    }
    return false;
}

bool IO_FileReader::IsArrayParameter()
{
    return PHSPACE::IsArrayParameter(* this->line);
}

void IO_FileReader::SkipLines(int numberOfLinesToSkip)
{
    PHSPACE::SkipLines(*file, numberOfLinesToSkip);
}

void IO_FileReader::ReadNextNonEmptyLine()
{
    PHSPACE::ReadNextNonEmptyLine(* this->file, * this->line);
}

void IO_FileReader::DumpLineContentToScreen()
{
    cout << *line << "\n";
}

std::string IO_FileReader::ReadNextWord()
{
    string word = PHSPACE::FindNextWord(* this->line, * this->separator);
    return word;
}

std::string IO_FileReader::ReadNextWord(const std::string &separator)
{
    string word = PHSPACE::FindNextWord(* this->line, separator);
    return word;
}

string IO_FileReader::ReadNextWordToLowerCase()
{
    string word = PHSPACE::FindNextWord(* this->line, * this->separator);

    PHSPACE::ToLowerCase(word);
    return word;
}

string IO_FileReader::ReadNextWordToLowerCase(const std::string &separator)
{
    string word = PHSPACE::FindNextWord(* this->line, separator);
    PHSPACE::ToLowerCase(word);
    return word;
}

void IO_FileReader::SkipReadSymbol(const string &stringSymbol)
{
    while (!this->ReachTheEndOfFile())
    {
        bool resultFlag = this->ReadNextMeaningfulLine();
        if (!resultFlag) break;

        string word = this->ReadNextWord();

        if (word == stringSymbol)
        {
            return;
        }
    }
}

void IO_FileReader::SkipReadWholeBlock()
{
    int countOfLeftBrackets  = 0;
    int countOfRightBrackets = 0;

    while (!this->ReachTheEndOfFile())
    {
        bool resultFlag = this->ReadNextMeaningfulLine();
        if (!resultFlag) break;

        string word = this->ReadNextWord();

        if (word == "{")
        {
            ++ countOfLeftBrackets;
        }
        else if (word == "}")
        {
            ++ countOfRightBrackets;
        }

        if (countOfLeftBrackets == countOfRightBrackets)
        {
            return;
        }
    }
}

}