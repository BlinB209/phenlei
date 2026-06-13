#include <iostream>
#include <stdlib.h>
#include "Data_SafeData.h"
#include "Data_Util.h"
#include "Constants.h"

namespace PHSPACE
{

Data_SafeData::Data_SafeData(const string &name, void *data, int type, int size)
{
    this->name = name;
    this->data = NULL;
    CopyData(data, type, size);
}

Data_SafeData::Data_SafeData(const Data_SafeData &rhs)
{
    name = rhs.name;
    this->data = 0;
    CopyData(rhs.data, rhs.type, rhs.size);
}

Data_SafeData::~Data_SafeData()
{
    this->DeleteData();
}

void Data_SafeData::OverWrite(void *data, int type, int size)
{
    //! Firstly, delete the old memory.
    DeleteData();
    //! Overwrite the old data.
    CopyData(data, type, size);
}

void Data_SafeData::ShowData()
{
    cout << "name = " << name << "\n";
    cout << "type = " << type << "\n";
    cout << "size = " << size << "\n";
    if (type == PHSPACE::PHSTRING)
    {
        string * mydata = reinterpret_cast<string *>(data);
        for (int i = 0; i < size; ++ i)
        {
            cout << "var[" << i << "] = " << mydata[i] << "\n";
        }
        cout << "\n";
    }
    else if (type == PHSPACE::PHINT)
    {
        int *mydata = reinterpret_cast<int *>(data);
        for (int i = 0; i < size; ++ i)
        {
            cout << "var[" << i << "] = " << mydata[i] << "\n";
        }
        cout << "\n";
    }
    else if (type == PHSPACE::PHFLOAT)
    {
        RFloat *mydata = reinterpret_cast<RFloat *>(data);
        for (int i = 0; i < size; ++ i)
        {
            cout << "var[" << i << "] = " << mydata[i] << "\n";
        }
        cout << "\n";
    }
    else if (type == PHSPACE::PHDOUBLE)
    {
        RDouble *mydata = reinterpret_cast<RDouble *>(data);
        for (int i = 0; i < size; ++ i)
        {
            cout << "var[" << i << "] = " << mydata[i] << "\n";
        }
        cout << "\n";
    }
}

bool Data_SafeData::operator < (const Data_SafeData &rhs) const
{
    return name < rhs.name;
}

Data_SafeData & Data_SafeData::operator = (const Data_SafeData &rhs)
{
    if (this == &rhs)
    {
        return *this;
    }
    this->name = rhs.name;

    OverWrite(rhs.data, rhs.type, rhs.size);
    return *this;
}

void Data_SafeData::CopyData(void *data, int type, int size)
{
    this->type = type;
    this->size = size;
    if (data == NULL) return;

    Data_Util::copy_data_sub(this->data, data, type, size);
}

void Data_SafeData::DeleteData()
{
    if (data == 0) return;

    if (type == PHSPACE::PHINT)
    {
        Data_Util::delete_void <int> (data);
    }
    else if (type == PHSPACE::PHFLOAT)
    {
        Data_Util::delete_void <RFloat> (data);
    }
    else if (type == PHSPACE::PHDOUBLE)
    {
        Data_Util::delete_void <RDouble> (data);
    }
    else if (type == PHSPACE::PHSTRING)
    {
        Data_Util::delete_void <string> (data);
    }
    else if (type == PHSPACE::PHBOOL)
    {
        Data_Util::delete_void <bool> (data);
    }
    else
    {
        cout << "No such type\n" << endl;
        abort();    //! It is temporary!!!
    }
    this->data = 0;
}

}