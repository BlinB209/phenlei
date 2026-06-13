#include "GlobalDataBase.h"
#include "Constants.h"
#include "DataContainer.h"

using namespace std;
#pragma warning (disable:913)
namespace PHSPACE
{
Data_Param * GlobalDataBase::globalParam = new Data_Param();
Data_Field * GlobalDataBase::globalField = new Data_Field();

void GlobalDataBase::UpdateData(const string &name, void *data, int type, int size)
{
    globalParam->UpdateData(name, data, type, size);
}

void GlobalDataBase::GetData(const string &name, void *data, int type, int size)
{
    globalParam->GetData(name, data, type, size);
}

int GlobalDataBase::GetSizeOfData(const string &name)
{
    return globalParam->GetSizeOfData(name);
}

void GlobalDataBase::UpdateDataPtr(const string &name, void *data)
{
    globalField->UpdateDataPtr(name, data);
};

void GlobalDataBase::DeleteDataPtr(const string &name)
{
    globalField->DeleteDataPtr(name);
};

void * GlobalDataBase::GetDataPtr(const string &name)
{
    return globalField->GetDataPtr(name);
};

set < Data_SafeData > * GlobalDataBase::GetBaseP()
{
    return globalParam->GetBaseP();
};

Data_Param * GlobalDataBase::GetDataPara() 
{
    return globalParam;
};

void GlobalDataBase::WriteString(DataContainer *cdata, const string &cs)
{
    int nlen = static_cast<int>(cs.length());
    cdata->Write(&nlen, sizeof(int));

    char *data = new char[nlen+1];

    cs.copy(data, nlen);
    data[nlen] = '\0';

    cdata->Write(data, nlen+1);

    delete [] data;
}

void GlobalDataBase::ReadString(DataContainer *cdata, string &cs)
{
    int nlen = 0;
    cdata->Read(&nlen, sizeof(int));

    char *data = new char[nlen+1];
    cdata->Read(data, nlen+1);

    cs = data;

    delete [] data;
}

void GlobalDataBase::WriteVoid(DataContainer *cdata, void *data, int type, int size)
{
    if (type == PHSTRING)
    {
        string *mydata = reinterpret_cast<string *>(data);
        for (int i = 0; i < size; ++ i)
        {
            WriteString(cdata, mydata[i]);
        }
    }
    else if (type == PHINT)
    {
        int *mydata = reinterpret_cast<int *>(data);
        cdata->Write(&mydata[0], size * sizeof(int));
    }
    else if (type == PHFLOAT)
    {
        RFloat *mydata = reinterpret_cast<RFloat *>(data);
        cdata->Write(&mydata[0], size * sizeof(RFloat));
    }
    else if (type == PHDOUBLE)
    {
        RDouble *mydata = reinterpret_cast<RDouble *>(data);
        cdata->Write(&mydata[0], size * sizeof(RDouble));
    }
}

void GlobalDataBase::UpdateData(DataContainer *cdata, const string &name, int type, int size)
{
    if (type == PHSTRING)
    {
        string *data = new string [size];
        for (int i = 0; i < size; ++ i)
        {
            ReadString(cdata, data[i]);
        }
        UpdateData(name, data, type, size);
        delete [] data;
    }
    else if (type == PHINT)
    {
        int *data = new int [size];
        for (int i = 0; i < size; ++ i)
        {
            cdata->Read(&data[i], sizeof(int));
        }
        UpdateData(name, data, type, size);
        delete [] data;
    }
    else if (type == PHFLOAT)
    {
        RFloat *data = new RFloat [size];
        for (int i = 0; i < size; ++ i)
        {
            cdata->Read(&data[i], sizeof(RFloat));
        }
        UpdateData(name, data, type, size);
        delete [] data;
    }
    else if (type == PHDOUBLE)
    {
        RDouble *data = new RDouble [size];
        for (int i = 0; i < size; ++ i)
        {
            cdata->Read(&data[i], sizeof(RDouble));
        }
        UpdateData(name, data, type, size);
        delete [] data;
    }
}

void GlobalDataBase::CompressData(DataContainer *&cdata)
{
    set < Data_SafeData > *basep = GetBaseP();
    set < Data_SafeData >::iterator iter;

    uint_t ndata = basep->size();

    cdata->Write(&ndata, sizeof(uint_t));

    for (iter = basep->begin(); iter != basep->end(); ++ iter)
    {
        int type = iter->GetType();
        int size = iter->GetSize();

        WriteString(cdata, iter->GetName());
        cdata->Write(&type, sizeof(int));
        cdata->Write(&size, sizeof(int));
        WriteVoid(cdata, iter->GetData(), type, size);
    }
}

void GlobalDataBase::DecompressData(DataContainer *cdata)
{
    //! It seems that move cdata to the beginning is very important.
    cdata->MoveToBegin();

    uint_t ndata = 0;
    cdata->Read(&ndata, sizeof(uint_t));

    for (uint_t i = 0; i < ndata; ++ i)
    {
        string name;
        int type;
        int size;

        ReadString(cdata, name);
        cdata->Read(&type, sizeof(int));
        cdata->Read(&size, sizeof(int));
        UpdateData(cdata, name, type, size);
    }
}

string GlobalDataBase::GetStrParaFromDB(const string &name)
{
    string data;
    GetData(name, &data, PHSTRING, 1);
    return data;
}

int GlobalDataBase::GetIntParaFromDB(const string &name)
{
    int data;
    GetData(name, &data, PHINT, 1);
    return data;
}

RDouble GlobalDataBase::GetDoubleParaFromDB(const string &name)
{
    RDouble data;
    GetData(name, &data, PHDOUBLE, 1);
    return data;
}

bool GlobalDataBase::IsExist(const string &name, int type, int size)
{
    return globalParam->IsExist(name, type, size);
}

RDouble * GlobalDataBase::GetRDoubleArrayFromDB(const string &name, int numberOfElements)
{
    RDouble *rRDoubleArray = new RDouble[numberOfElements];
    GetData(name, rRDoubleArray, PHDOUBLE, numberOfElements);
    return rRDoubleArray;
}

void GlobalDataBase::GetRDoubleVectorFromDB(vector< RDouble > &parameter, const string &nameOfRFloatArray, int numberOfElements)
{
    parameter.resize(numberOfElements);
    GetData(nameOfRFloatArray, &parameter[0], PHDOUBLE, numberOfElements);
}

LIB_EXPORT int GetTaskCode()
{
    int nsimutask = GlobalDataBase::GetIntParaFromDB("nsimutask");
    return nsimutask;
}

}