#include <string>
#include <iostream>
#include <stdlib.h>
#include "Data_Util.h"
#include "Constants.h"

using namespace std;

namespace PHSPACE
{

void Data_Util::copy_data_sub(void *&target, void *source, int type, int size)
{
    if (type == PHSPACE::PHINT)
    {
        copy_void <int> (target, source, size);
    }
    else if (type == PHSPACE::PHFLOAT)
    {
        copy_void <RFloat> (target, source, size);
    }
    else if (type == PHSPACE::PHDOUBLE)
    {
        copy_void <RDouble> (target, source, size);
    }
    else if (type == PHSPACE::PHSTRING)
    {
        copy_void <std::string> (target, source, size);
    }
    else if (type == PHSPACE::PHBOOL)
    {
        copy_void <bool> (target, source, size);
    }
    else
    {
        cout << "No such type\n" << endl;
        abort();    //! It is temporary!!!
    }
}

void Data_Util::copy_data_exist(void *target, void *source, int type, int size)
{
    if (type == PHSPACE::PHINT)
    {
        copy_void_exist <int> (target, source, size);
    }
    else if (type == PHSPACE::PHFLOAT)
    {
        copy_void_exist <RFloat> (target, source, size);
    }
    else if (type == PHSPACE::PHDOUBLE)
    {
        copy_void_exist <RDouble> (target, source, size);
    }
    else if (type == PHSPACE::PHSTRING)
    {
        copy_void_exist <std::string> (target, source, size);
    }
    else if (type == PHSPACE::PHBOOL)
    {
        copy_void_exist <bool> (target, source, size);
    }
    else
    {
        cout << "No such type\n" << endl;
        abort();    //! It is temporary!!!
    }
}

}