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
//! @file      AMRDef.h
//! @brief     Explain this file briefly.
//! @author    Bell.

#pragma once
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <set>
#include <vector>
#include <algorithm>
#include <string>
using namespace std;

namespace PHSPACE
{
template < typename T >
class PHVector1D : public vector< T >
{
public:
    typedef typename vector< T >::size_type size_type;
    typedef typename vector< T >::value_type value_type;

public:
    PHVector1D() {};
    ~PHVector1D() {};

public:
    PHVector1D < T > (size_type count, const value_type &value) : vector< T >(count, value)
    {
    }

    explicit PHVector1D < T > (size_type count) : vector< T >(count)
    {
    }

    PHVector1D < T > (T *first, T *last) : vector< T >(first, last)
    {
    }

public:
    T * GetPointer() { return & (*this)[0]; }
};

template < typename T >
class PHVector2D
{
public:
    PHVector2D() {};
    ~PHVector2D() {};

public:
    explicit PHVector2D < T > (int count)
    {
        data2D.resize(count);
    }
    size_t size() const { return data2D.size(); }
    void push_back(const PHVector1D < T > &data1D) { data2D.push_back(data1D); }

protected:
    PHVector1D < PHVector1D < T > > data2D;

public:
    PHVector1D < T > & operator [] (const size_t &index) { return data2D[index]; }
    const PHVector1D < T > & operator [] (const size_t &index) const { return data2D[index]; }
    void resize(size_t nSize) { data2D.resize(nSize); }
};

}