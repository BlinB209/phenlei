#include "Data_IndexData.h"

namespace PHSPACE
{

Data_IndexData::Data_IndexData():data(0)
{
}

Data_IndexData::Data_IndexData(const string &name, void *data)
{
    this->name = name;
    this->data = data;
}

Data_IndexData::~Data_IndexData()
{
}

bool Data_IndexData::operator < (const Data_IndexData &rhs) const
{
    return name < rhs.name;
}

}