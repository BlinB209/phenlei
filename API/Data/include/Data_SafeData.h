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
//! @file      Data_SafeData.h
//! @brief     Data_SafeData is a basic class of data storage.
//! @author    He Xin, Zhang Jian, reference from Dr. Wang.

#pragma once
#include <string>

using namespace std;

namespace PHSPACE
{

//! @brief Data_SafeData is a basic class of data storage.
//! It stores the data's name, type, size, and allocate memory for the data.
//! There is a pointer to the data.
class Data_SafeData
{
private:
    string name;
    void  *data;
    int    type;
    int    size;

public:
    //! Construct function.
    Data_SafeData(const string &name, void *data, int type, int size);

    //! Copy Construct function.
    Data_SafeData(const Data_SafeData &rhs);

    //! Destructor, free memory.
    ~Data_SafeData();

    //! OverWrite old data by new type and size data.
    void OverWrite(void *data, int type, int size);

    //! Print Data to window.
    void ShowData();

    //! Return data's name.
    string GetName() const;

    //! Return data's type.
    int GetType() const;

    //! Return data's size.
    int GetSize() const;

    //! Return pointer to the data.
    void * GetData() const;

    //! Override operator < used by set<Data_SafeData>, sort by name.
    bool operator < (const Data_SafeData &rhs) const;

private:
    //! Override operator < used by set<Data_SafeData>.
    Data_SafeData & operator = (const Data_SafeData &rhs);

    //! Invoking function copy_data_sub.
    void CopyData(void *data, int type, int size);

    //! Invoking function delete_void.
    void DeleteData();
};

}

#include "Data_SafeData.hxx"