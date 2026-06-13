#include <stdlib.h>
#include "Data_Param.h"
#include "Data_Util.h"

namespace PHSPACE
{

LIB_EXPORT Data_Param::Data_Param()
{
    dataSet = new set < Data_SafeData >;
}

LIB_EXPORT Data_Param::Data_Param(const Data_Param &rhs)
{
    dataSet = new set < Data_SafeData >;
    set < Data_SafeData >::iterator iter;
    for (iter = rhs.dataSet->begin(); iter != rhs.dataSet->end(); ++ iter)
    {
        dataSet->insert(*iter);
    }
}

LIB_EXPORT Data_Param::~Data_Param()
{
    delete dataSet;
    dataSet = NULL;
}

LIB_EXPORT void Data_Param::UpdateData(const string &name, void *data, int type, int size)
{
    Data_SafeData datas = Data_SafeData(name, data, type, size);
    //datas.ShowData();

    set < Data_SafeData >::iterator iter;
    iter = dataSet->find(datas);
    if (iter != dataSet->end())
    {
        dataSet->erase(iter);
        dataSet->insert(datas);
    }
    else
    {
        dataSet->insert(datas);
    }
}

LIB_EXPORT void Data_Param::GetData(const string &name, void *data, int type, int size)
{
    set < Data_SafeData >::iterator iter;
    Data_SafeData datas = Data_SafeData(name, data, type, size);
    iter = dataSet->find(datas);

    /*static set<string> bcNameSet;
    static int count = 0;*/

    if (iter != dataSet->end())
    {
        CopyDataValue(data, iter->GetData(), iter->GetType(), iter->GetSize());
 
        //bcNameSet.insert(name);
        //count ++;
    }
    else
    {
        cout << "This data does not exist: " << name.c_str() << endl;
        exit(0);    //! It is temporary!!!
    }

    /*int iCount = static_cast<int>(bcNameSet.size()) + 10000;
    if (iCount < count)
    {
        cout << "  Warning: this data has been read  more times: " << name.c_str() << endl;
    }*/
}

LIB_EXPORT int Data_Param::CheckDataExist(const string &name)
{
    int data;
    set < Data_SafeData >::iterator iter;
    Data_SafeData datas = Data_SafeData(name, &data, 1, 1);
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

LIB_EXPORT int Data_Param::GetSizeOfData(const string &name)
{
    set < Data_SafeData >::iterator iter;
    Data_SafeData datas = Data_SafeData(name, 0, 0, 0);
    iter = dataSet->find(datas);
    if (iter != dataSet->end())
    {
        return iter->GetSize();
    }
    else
    {
        cout << "This data does not exist: " << name.c_str() << endl;
        exit(0);    //! It is temporary!!!
    }
}

LIB_EXPORT bool Data_Param::IsExist(const string &name, int type, int size)
{
    set < Data_SafeData >::iterator iter;
    Data_SafeData datas = Data_SafeData(name, 0, type, size);
    iter = dataSet->find(datas);
    if (iter != dataSet->end())
    {
        return true;
    }
    else
    {
        return false;
    }
}

LIB_EXPORT set < Data_SafeData > * Data_Param::GetBaseP() const
{
    return dataSet;
};

LIB_EXPORT Data_Param & Data_Param::operator = (const Data_Param &rhs)
{
    //! Only do assignment if rhs is a different object from this.
    if (this != &rhs)
    {
        //! Deallocate, allocate new space, copy values.
        delete dataSet;
        dataSet = new set < Data_SafeData >;
        set < Data_SafeData >::iterator iter;
        for (iter = rhs.dataSet->begin(); iter != rhs.dataSet->end(); ++ iter)
        {
            dataSet->insert(*iter);
        }
    }

    return *this;
}

void Data_Param::CopyDataValue(void *dataout, void *data, int type, int size)
{
    if (data == 0) return;
    Data_Util::copy_data_exist(dataout, data, type, size);
}

}