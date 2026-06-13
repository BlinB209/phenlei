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
//! @file      Data_IndexData.h
//! @brief     Data_IndexData is a basic class of data storage that only stores the pointer to the data.
//! @author    He Xin, Zhang Jian.

#pragma once
#include <string>

using namespace std;

namespace PHSPACE
{

//! @brief Data_IndexData stores the data's name and address(pointer to the data).
//! It doesn't allocate or free memory for the data.
class Data_IndexData
{
private:
    string name;
    void *data;
public:
    //! Default constructor.
    Data_IndexData();

    //! Constructor using name and data.
    Data_IndexData(const string &name, void *data);

    //! Destructor, do nothing.
    ~Data_IndexData();

    //! Return data's name.
    string GetName() const;

    //! Return pointer to the data.
    void * GetData() const;

    //! Override operator < used by set<Data_IndexData>, sort by name.
    bool operator < (const Data_IndexData &rhs) const;
};

}

#include "Data_IndexData.hxx"