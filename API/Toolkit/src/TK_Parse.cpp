#include <map>
#include "TK_Parse.h"
#include "Data_Param.h"
#include "PHHeader.h"
#include "TK_Exit.h"
#include "Constants.h"
#include "TK_Log.h"

namespace PHSPACE
{
LIB_EXPORT void TK_Parse::ReadBasicData(fstream &file, Data_Param *database)
{
    string line, keyword, name, word;
    //string *value;
    vector <string> values;
    //string separator  = " =\t\r\n#$,;\"";
    //! \t is the tab key.
    string separator = " =\r\n\t#$,;\"";
    string sep1 = "\r\n";
    string::size_type npos = string::npos;

    string errMsg = "error in parameter file";

    map <string, int> keywordmap;

    keywordmap.insert(pair<string, int>("int"   , PHINT));
    keywordmap.insert(pair<string, int>("float" , PHFLOAT));
    keywordmap.insert(pair<string, int>("double", PHDOUBLE));
    keywordmap.insert(pair<string, int>("string", PHSTRING));

    while (!file.eof())
    {
        getline(file, line);
        if (line == "") continue;
        FindNextWord(line, word, sep1);

        if (word == "{") continue;    //! Get into sub-domain.
        if (word == "}") break;       //! Get out of sub-domain.

        if (word.substr(0, 1) == "#" || word.substr(0, 2) == "//") continue;

        line = FindNextWord(line, keyword, separator);

        if (keyword == "") continue;

        line = FindNextWord(line, name, separator);

        int count = 0;
        if (name.find_first_of("[") != npos) ++ count;
        if (name.find_first_of("]") != npos) ++ count;
        int arraysize = 1;
        if (count == 2)
        {
            //! Array mode.
            string line1 = name;
            //string sep = " =\t\r\n#$,;\"[]";
            string sep = " =\r\n\t#$,;\"[]";
            line1 = FindNextWord(line1, name, sep);

            string value;
            line = FindNextWord(line, value, sep);
            while (value != "")
            {
                values.push_back(value);
                line = FindNextWord(line, value, sep);
            }

            arraysize = static_cast<int>(values.size());
        }
        else if (count == 0)
        {
            string value;
            if (name == "speciesName" || name == "initMassFraction")
            {
                value = line;
            }
            else
            {
                line = FindNextWord(line, value, separator);
            }
            values.push_back(value);
        }
        else
        {
            ostringstream oss;
            oss << errMsg << "\n";
            TK_Exit::ExceptionExit(oss.str(), true);
        }

        int type, size = arraysize;
        type = keywordmap[keyword];
        if (type == PHSTRING)
        {
            string *data = new string [size];
            for (int i = 0; i < size; ++ i)
            {
                data[i] = values[i];
            }
            database->UpdateData(name, data, type, size);
            delete [] data;
        }
        else if (type == PHINT)
        {
            int *data = new int [size];
            for (int i = 0; i < size; ++ i)
            {
                from_string<int>(data[i], values[i], std::dec);
            }

            database->UpdateData(name, data, type, size);
            delete [] data;
        }
        else if (type == PHFLOAT)
        {
            RFloat *data = new RFloat [size];
            for (int i = 0; i < size; ++ i)
            {
                from_string<RFloat>(data[i], values[i], std::dec);
            }
            database->UpdateData(name, data, type, size);
            delete [] data;
        }
        else if (type == PHDOUBLE)
        {
            RDouble *data = new RDouble [size];
            for (int i = 0; i < size; ++ i)
            {
                from_string<RDouble>(data[i], values[i], std::dec);
            }
            database->UpdateData(name, data, type, size);
            delete [] data;
        }
        //delete [] value;
        values.resize(0);
    }
}

void ToLowerCase(string &word)
{
    transform(word.begin(), word.end(), word.begin(), GSStringToLowerCase());
}

void ToUpperCase(string &word)
{
    transform(word.begin(), word.end(), word.begin(), GSStringToUpperCase());
}

void TrimBlanks(string &source)
{
    string::size_type firstindex, nextindex;
    firstindex = source.find_first_not_of(" ");
    nextindex  = source.find_last_not_of(" ");

    source = source.substr(firstindex, nextindex - firstindex + 1);
}

void SkipLines(fstream &file, int nline)
{
    string line;
    for (int i = 0; i < nline; ++ i)
    {
        getline(file, line);
        if (file.eof()) return;
    }
}

void SkipWords(string &line, int nword)
{
    string word;
    string separator  = " =\t\r\n#$,;\"'";
    for (int i = 0; i < nword; ++ i)
    {
        line = FindNextWord(line, word, separator);
    }
}

bool IsEmptyLine(const string &line)
{
    if (line == "")
    {
        return true;
    }
    else
    {
        const string notReadableSeparator = " \r\n\t";
        string word;
        PHSPACE::FindNextWord(line, word, notReadableSeparator);
        return word == "";
    }
}

void ReadNextNonEmptyLine(fstream &file, string &line)
{
    bool isSpaceLine = true;

    do
    {
        PHSPACE::ReadNextLine(file, line);

        for (string::iterator iter = line.begin(); iter != line.end(); ++ iter)
        {
            if (! isspace(* iter))
            {
                isSpaceLine = false;
                break;
            }
        }
        if (file.eof()) break;
    }
    while (isSpaceLine);

    return;
}

string FindNextWord(const string &source, string &word, const string &separator)
{
    int firstIndex, nextIndex, notFindIndex = -1;
    string emptyString = "";
    firstIndex = static_cast<int>(source.find_first_not_of(separator, 0));

    if (firstIndex == notFindIndex)
    {
        word = emptyString;
        return emptyString;
    }
    nextIndex = static_cast<int>(source.find_first_of(separator, firstIndex));
    if (nextIndex == notFindIndex)
    {
        word = source.substr(firstIndex);
        return emptyString;
    }
    else
    {
        word = source.substr(firstIndex, nextIndex - firstIndex);
        return source.substr(nextIndex);
    }
}

bool IsArrayParameter(const string &lineOfName)
{
    if (lineOfName.find_first_of("[") == string::npos)
    {
        return false;
    }

    if (lineOfName.find_first_of("]") == string::npos)
    {
        return false;
    }

    return true;
}

bool IsCommentLine(const string &line)
{
    const string notReadableSeparator = " \r\n\t";
    string word;
    PHSPACE::FindNextWord(line, word, notReadableSeparator);
    return (word.substr(0, 1) == "#" || word.substr(0, 2) == "//");
}

void ReadNextLine(fstream &file, string &line)
{
    std::getline(file, line);
}

string FindNextWord(string &source, string &word, const string &separator)
{
    string::size_type firstindex, nextindex, npos = string::npos;
    string nullstring = "";
    firstindex = source.find_first_not_of(separator, 0);

    if (firstindex == npos)
    {
        word = nullstring;
        return nullstring;
    }

    nextindex = source.find_first_of(separator, firstindex);
    if (nextindex == npos)
    {
        word = source.substr(firstindex);
        return nullstring;
    }
    else
    {
        word = source.substr(firstindex, nextindex - firstindex);
        return source.substr(nextindex);
    }
    return nullstring;
}

string FindNextWord(string &source, const string &separator)
{
    int firstIndex, nextIndex, notFindIndex = -1;
    string emptyString = "";
    firstIndex = static_cast<int>(source.find_first_not_of(separator, 0));

    if (firstIndex == notFindIndex)
    {
        return emptyString;
    }

    nextIndex = static_cast<int>(source.find_first_of(separator, firstIndex));
    if (nextIndex == notFindIndex)
    {
        string word = source.substr(firstIndex);
        source = emptyString;
        return word;
    }
    else
    {
        string word = source.substr(firstIndex, nextIndex - firstIndex);
        source = source.substr(nextIndex);
        return word;
    }
}

void GetNameExt(const string &fullname, string &fname, string &fext, const string &sep)
{
    basic_string <char>::size_type index;

    if (fullname.empty()) return;

    index = fullname.find_last_of(sep);
    fname = fullname.substr(0, index);
    fext  = fullname.substr(index+1, fullname.length() - index - 1);
}

void ReadNewLine(fstream &file, string &line)
{
    bool IsSpaceLine = true;

    do
    {
        getline(file, line);

        for (string::iterator iter = line.begin(); iter != line.end(); ++ iter)
        {
            if (!isspace(*iter))
            {
                IsSpaceLine = false;
                break;
            }
        }
        if (file.eof()) break;
    }
    while (IsSpaceLine);

    return;
}

string Dec2Hex(int data_in)
{
    stringstream oss;
    string temp;

    oss << setiosflags(ios::uppercase) << hex << data_in;
    oss >> temp;

    return temp;
}

}