#include "Data_ParamFieldSuite.h"
#include "Data_Param.h"
#include "Data_Field.h"

namespace PHSPACE
{

LIB_EXPORT Data_ParamFieldSuite::Data_ParamFieldSuite(bool init)
{ 
    if (init)
    {
        param = new Data_Param();
        field = new Data_Field();
    }
    else
    {
        param = 0;
        field = 0;
    }
}

LIB_EXPORT Data_ParamFieldSuite::~Data_ParamFieldSuite()
{
    if (param)
    {
        delete param;
    }
    if (field)
    {
        delete field;
    }
}

LIB_EXPORT void Data_ParamFieldSuite::UpdateData(const string &name, void *data, int type, int size)
{
    param->UpdateData(name, data, type, size);
}

LIB_EXPORT void Data_ParamFieldSuite::GetData(const string &name, void *data, int type, int size)
{
    param->GetData(name, data, type, size);
}

LIB_EXPORT void Data_ParamFieldSuite::UpdateDataPtr(const string &name, void *data)
{
    field->UpdateDataPtr(name, data);
};

LIB_EXPORT void * Data_ParamFieldSuite::GetDataPtr(const string &name) const
{
    return field->GetDataPtr(name);
};

LIB_EXPORT void Data_ParamFieldSuite::DeleteDataPtr(const string &name)
{
    field->DeleteDataPtr(name);
}

}