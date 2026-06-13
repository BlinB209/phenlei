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
//! @file      Data_Param.h
//! @brief     Data_Param has a set of data which stores parameters by name.
//!            Users can update or obtain the data by its name, type and size.
//!            Data_Param is usually used to store the parameters from *.hypara file, such as refReNumber, viscousName, limitVariables, et. al.
//! @author    He Xin, Zhang Jian.

#pragma once
#include <set>
#include "Data_SafeData.h"
#include "LIB_Macro.h"

namespace PHSPACE
{

class Data_Param
{
private:
    //! A set that contains Data_SafeData type data.
    //! Data_SafeData has data and its' name, type and size.
    set < Data_SafeData > *dataSet;

public:
    //! Constructor.
    LIB_EXPORT Data_Param();

    //! Copy constructor.
    LIB_EXPORT Data_Param(const Data_Param &rhs);

    //! Destructor, free memory for dataSet.
    LIB_EXPORT ~Data_Param();

    //! Add or update a data(with name, type, size) to dataSet.
    //! If there was a data with same name, replace the old with the new.
    //! Obtain the data by its name, type and size.
    //! @param[in]  name    Source data's name.
    //! @param[in]  data    Source data's address.
    //! @param[in]  type    Source data's type.
    //! @param[in]  size    Source data's size(length of data array).
    LIB_EXPORT void UpdateData(const string &name, void *data, int type, int size);

    //! Obtain the data by its name, type and size.
    //! @param[in]  name    Target data's name.
    //! @param[in]  data    Target data's address.
    //! @param[in]  type    Target data's type.
    //! @param[in]  size    Target data's size(length of data array).
    LIB_EXPORT void GetData(const string &name, void *data, int type, int size);

    LIB_EXPORT int CheckDataExist(const string &name);

    //! Obtain the data's size (length of data array) by its name.
    LIB_EXPORT int GetSizeOfData(const string &name);

    //! Find one Data_Param object whether exist in dataSet.
    LIB_EXPORT bool IsExist(const string &name, int type, int size);

    //! Obtain a reference pointing to the underlying data set that stores the parameter data.
    LIB_EXPORT set < Data_SafeData > * GetBaseP() const;

    //! Override the operator =, assign another Data_Param instance to self.
    //! It will erase all the old data, and copy new data from right hand side.
    LIB_EXPORT Data_Param & operator = (const Data_Param &rhs);

private:
    //! Invoking function copy_data_exist.
    void CopyDataValue(void *dataout, void *data, int type, int size);
};

}