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
//! @file      GlobalDataBase.h
//! @brief     GlobalDataBase class is similar to a database which stores global parameters or data pointer.
//!            Users can update or obtain the data by its name, type and size.
//!            It is different from Data_ParamFieldSuite class.
//!            It forbids users to create instance and only provides some interface function to users to access the data.
//!            These interface functions include UpdateData, GetData, UpdateDataPtr, GetDataPtr, etc.
//!            Example:
//!                 int nl = 5; GlobalDataBase::GetData("nl", &nl, PHINT, 1);
//!                 GlobalDataBase::UpdateData("dtmin", &dtmin, PHDOUBLE, 1);
//!                 RFloat *prim_inf = new RFloat[nl+nchem]; GlobalDataBase::UpdateDataPtr("prim_inf", prim_inf);
//!                 RFloat *prim_inf = reinterpret_cast<RFloat *>(GlobalDataBase::GetDataPtr("prim_inf")).
//! @author    Zhang Jian, He Xin.

#pragma once
#include "Data_Param.h"
#include "Data_Field.h"
#include "Precision.h"
#include <vector>

using namespace std;

namespace PHSPACE
{
class DataContainer;

class GlobalDataBase
{
private:
    //! Zonal control parameters.
    //! Data_Param is usually used to store the parameters from *.hypara file, such as refReNumber, viscousName, limitVariables, et. al.
    static Data_Param *globalParam;

    //! Zonal data pointers.
    static Data_Field *globalField;
public:
    //! Add or update a data(with name, type, size) to global database.
    //! If there was a data with same name, replace the old with the new.
    //! Obtain the data by its name, type and size.
    //! @param[in]  name    Source data's name.
    //! @param[in]  data    Source data's address.
    //! @param[in]  type    Source data's type.
    //! @param[in]  size    Source data's size(length of data array).
    static void UpdateData(const string &name, void *data, int type, int size);

    //! Obtain the data by its name, type and size from global database.
    //! @param[in]  name    Target data's name.
    //! @param[in]  data    Target data's address.
    //! @param[in]  type    Target data's type.
    //! @param[in]  size    Target data's size(length of data array).
    static void GetData(const string &name, void *data, int type, int size);

    //! Obtain the data's size (length of data array) by its name.
    static int GetSizeOfData(const string &name);

    //! Add or update a data(with name, type, size) to global database.
    //! If there was a data with same name, replace the old with the new.
    //! @param[in]  name    Source data's name.
    //! @param[in]  data    Pointer the Source data.
    static void UpdateDataPtr(const string &name, void *data);

    static void DeleteDataPtr(const string &name);

    //! Obtain data's pointer(address) from global database by name.
    //! If the dataSet doesn't has a same name data, returns null.
    //! @param[in]  name    Source data's name.
    static void * GetDataPtr(const string &name);

    //! Obtain the underlying data set which stores global parameters.
    static set < Data_SafeData > * GetBaseP();

    //! Obtain the private member globalParam that stores global parameters.
    static Data_Param * GetDataPara();

    //! Compress global parameters stored in globalParm into a DataContainer cdata.
    static void CompressData(DataContainer *&cdata);

    //! Decompress data from a DataContainer cdata, and update them to globalParam.
    static void DecompressData(DataContainer *cdata);

    //! UpdateData one data from a DataContainer cdata.
    //! @param[in]  name    data's name.
    //! @param[in]  type    data's type.
    //! @param[in]  size    length of data array.
    //! @param[out]  cdata    the DataContainer which stores one data.
    static void UpdateData(DataContainer *cdata, const string &name, int type, int size);

    //! Get one string type data from database by name.
    static string GetStrParaFromDB(const string &name);

    //! Get one integer type data from database by name.
    static int GetIntParaFromDB(const string &name);

    //! Get one double type data from database by name.
    static RDouble GetDoubleParaFromDB(const string &name);

    //! Return if the data named by 'name' exist or not.
    static bool IsExist(const string &name, int type, int size);

    //! Get double type vector from database by name.
    static RDouble * GetRDoubleArrayFromDB(const string &name, int numberOfElements);

    static void GetRDoubleVectorFromDB(vector< RDouble > &parameter, const string &nameOfRFloatArray, int numberOfElements);

private:
    //! Constructor is private, so instantiation of this class is forbidden.
    GlobalDataBase();

    //! WriteString function is used to write string object cs into cdata.
    //! @param[in]  cs       string type,source object.
    //! @param[in]  cdata    receive data from cs,target object.
    static void WriteString(DataContainer *cdata, const string &cs);

    //! ReadString function is used to write cdata into string object cs.
    //! @param[in] cs       string object,target object.
    //! @param[in] cdata    cope data to cs,source object.
    static void ReadString(DataContainer *cdata, string &cs);

    //! WriteVoid function is used to select suitable write function, restricted by type.
    //! @param[in]  cdata   receive and store data from data object.
    //! @param[in]  data    source object.
    //! @param[in]  type    the data type to be transformed.
    //! @param[in]  size    length of data array.
    static void WriteVoid(DataContainer *cdata, void *data, int type, int size);

    static void SetDefaultParamInfo();
};

//! Get task type.
LIB_EXPORT int GetTaskCode();
}
