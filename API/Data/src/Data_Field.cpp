#include "Data_Field.h"

namespace PHSPACE
{

LIB_EXPORT Data_Field::Data_Field()
{
    dataSet = new set < Data_IndexData >;
}

LIB_EXPORT Data_Field::~Data_Field()
{
    delete dataSet;
}

LIB_EXPORT void Data_Field::UpdateDataPtr(const string &name, void *data)
{
    pair<set < Data_IndexData > ::iterator, bool> ret = dataSet->insert(Data_IndexData(name, data));
    if (!ret.second)
    {
        cout << "    Warning: UpdateDataPtr false, inserted repeatedly: " << name << " !!!" << endl;
    }
}

LIB_EXPORT void * Data_Field::GetDataPtr(const string &name) const
{
    Data_IndexData data(name, 0);
    set < Data_IndexData >::iterator it = dataSet->find(data);
    if (it != dataSet->end())
    {
        return it->GetData();
    }
    else
    {
        //cout << "Variable " << name << " is not in the database.\n";
        return 0;
    }
}

LIB_EXPORT void Data_Field::DeleteDataPtr(const string &name)
{
    Data_IndexData data(name, 0);
    set < Data_IndexData >::iterator it = dataSet->find(data);
    if (it != dataSet->end())
    {
        dataSet->erase(it);
    }
}

LIB_EXPORT int Data_Field::CheckDataExist(const string &name)
{
    int data;
    set < Data_IndexData >::iterator iter;
    Data_IndexData datas = Data_IndexData(name, &data);
    iter = dataSet->find(datas);
    if (iter != dataSet->end())
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

LIB_EXPORT set < Data_IndexData > * Data_Field::GetDataSet() const
{ 
    return dataSet;
}

}