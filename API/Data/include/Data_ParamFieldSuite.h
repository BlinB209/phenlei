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
//! @file      Data_ParamFieldSuite.h
//! @brief     Data_ParamFieldSuite embeds a Data_Param type pointer and a Data_Field pointer into one class.
//!            It wraps Data_Param's public interface(UpdateData, GetData) and Data_Field's public interface(UpdateDataPtr, GetDataPtr).
//! @author    He Xin, Zhang Jian, reference from Dr. Wang.

#pragma once
#include <string>
#include "LIB_Macro.h"

using namespace std;

namespace PHSPACE
{

//! Forward declaration.
class Data_Param;
class Data_Field;

//! @brief Data_ParamFieldSuite wraps Data_Param's public interface(UpdateData, GetData) 
//!        and Data_Field's public interface(UpdateDataPtr, GetDataPtr).
class Data_ParamFieldSuite
{
private:
    //! Pointer to a Data_Param instance which stores parameters.
    Data_Param *param;

    //! Pointer to a Data_Field instance which stores field data's pointer.
    Data_Field *field;

public:
    //! Constructor.
    //! @param[in]  init    Create instance for param and field if true(default).
    LIB_EXPORT Data_ParamFieldSuite(bool init = true);

    //! Destructor, free param and field.
    LIB_EXPORT ~Data_ParamFieldSuite();

    //! Add or update a data(with name, type, size) for param.
    //! If there was a data with same name, replace the old with the new.
    //! Obtain the data by its name, type and size.
    //! @param[in]  name    Source data's name.
    //! @param[in]  data    Source data's address.
    //! @param[in]  type    Source data's type.
    //! @param[in]  size    Source data's size.
    LIB_EXPORT void UpdateData(const string &name, void *data, int type, int size);

    //! Obtain the data by its name, type and size.
    //! @param[in]  name    Target data's name.
    //! @param[in]  data    Target data's address.
    //! @param[in]  type    Target data's type.
    //! @param[in]  size    Target data's size.
    LIB_EXPORT void GetData(const string &name, void *data, int type, int size);

    //! Add or update a data(with name, type, size) to field.
    //! If there was a data with same name, replace the old with the new.
    //! @param[in]  name    Source data's name.
    //! @param[in]  data    Pointer the Source data.
    LIB_EXPORT void UpdateDataPtr(const string &name, void *data);

    //! Obtain data's address from dataSet by name.
    //! If the dataSet doesn't has a same name data, returns null.
    //! @param[in]  name    Source data's name.
    LIB_EXPORT void * GetDataPtr(const string &name) const;

    //! Remove data's address from dataSet by name.
    //! @param[in]  name    Source data's name.
    LIB_EXPORT void DeleteDataPtr(const string &name);
};

}