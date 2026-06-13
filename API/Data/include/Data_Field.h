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
//! @file      Data_Field.h
//! @brief     Data_Field has a set of data which stores field data by name.
//!            It only stores field data's first address, i.e. a pointer to data's memory.
//!            Users can update or obtain the data's pointer by its name.
//!            Usually, users use Data_Field to store flow field data, such as rho, u, p, and so on.
//! @author    He Xin, Zhang Jian.

#pragma once
#include <set>
#include "Data_IndexData.h"
#include "LIB_Macro.h"

namespace PHSPACE
{

class Data_Field
{
private:
    //! A set that stores Data_IndexData type data.
    //! Data_IndexData has data's name and address(pointer to the data).
    set < Data_IndexData > *dataSet;

public:
    //! Constructor, create Data_Field.
    LIB_EXPORT Data_Field();

    //! Destructor, free Data_Field.
    LIB_EXPORT ~Data_Field();

    //! Add or update a data(with name, type ,size) to dataSet.
    //! If there was a data with same name, replace the old with the new.
    //! @param[in]  name    Source data's name.
    //! @param[in]  data    Pointer the Source data.
    LIB_EXPORT void UpdateDataPtr(const string &name, void *data);

    //! Obtain data's address from dataSet by name.
    //! If the dataSet doesn't has a same name data, returns null.
    //! @param[in]  name    Source data's name.
    LIB_EXPORT void * GetDataPtr(const string &name) const;

    //! Delete one of the items in the dataSet by name.
    //! @param[in]  name    Item's name needs to be deleted.
    LIB_EXPORT void DeleteDataPtr(const string &name);

    LIB_EXPORT int CheckDataExist(const string &name);

    //! Obtain the underlying dataSet that stores all data.
    LIB_EXPORT set < Data_IndexData > * GetDataSet() const;
};

}