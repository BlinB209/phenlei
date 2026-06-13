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
//! @file      DataStruct_Sort.h
//! @brief     DataStruct_Sort is usually used combined with set of STL, e.g. set<DataStruct_Sort<int>>.
//!            DataStruct_Sort objects stored in set is sorted by comparing its' value.
//! @author    Zhang Jian.

#pragma once
namespace PHSPACE
{
template < typename T >
class DataStruct_Sort
{
public:
    //! Value used to compare.
    T value;

    //! Index.
    int index;

public:
    //! Constructor.
    DataStruct_Sort();

    //! Constructor with value and index.
    DataStruct_Sort(T &value, int index);

    //! Overload operator < used to sort.
    bool operator < (const DataStruct_Sort &rhs) const;

    bool operator > (const DataStruct_Sort &rhs) const;
};

#include "DataStruct_Sort.hxx"
}
