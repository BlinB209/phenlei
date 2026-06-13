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
//! @file      Pointer.h
//! @brief     Pointer operation for multi-dimensional array.
//! @author    He Xin, Bell.

#pragma once
#include <vector>
#include <cstddef>
using namespace std;

namespace PHSPACE
{
template < typename T >
T * HeaderPointer(vector< T > &field)
{
    if (field.size() != 0)
    {
        return &field[0];
    }
    return 0;
}

template < typename T >
T * NewPointer(int length)
{
    T *data = new T [length];
    
    return data;
}

template < typename T >
T ** NewPointer2(int fistDim, int secondDim)
{
    T **data = new T * [fistDim];
    data[0] = new T[fistDim * secondDim];
    for (int m = 1; m < fistDim; ++ m)
    {
        data[m] = &data[m-1][secondDim];
    }
    return data;
}

template < typename T >
T ** NewPointer2(int nm, int *nlen)
{
    int count = 0;

    for (int m = 0; m < nm; ++ m)
    {
        count += nlen[m];
    }

    T **vv = new T * [nm];
    vv[0] = new T [count];
    for (int m = 1; m < nm; ++ m)
    {
        vv[m] = &vv[m-1][nlen[m-1]];
    }
    return vv;
}

template < typename T >
T *** NewPointer3(int fistDim, int secondDim, int thirdDim)
{
    T *** data = new T ** [fistDim];
    data[0]    = new T *  [fistDim * secondDim];
    data[0][0] = new T    [fistDim * secondDim * thirdDim];

    for (int i = 1; i < fistDim; ++ i)
    {
        data[i] = &data[i-1][secondDim];
    }

    for (int j = 1; j < secondDim; ++ j)
    {
        data[0][j] = &data[0][j-1][thirdDim];
    }

    for (int i = 1; i < fistDim; ++ i)
    {
        data[i][0] = &data[i-1][secondDim-1][thirdDim];
        for (int j = 1; j < secondDim; ++ j)
        {
            data[i][j] = &data[i][j-1][thirdDim];
        }
    }

    return data;
}

template < typename T >
T *** NewPointer3(int nCell, int nEqn, const int *nlen)
{
    int count = 0;

    for (int m = 0; m < nCell; ++ m)
    {
        count += nEqn * nlen[m];
    }

    T *** v3 = new T ** [nCell];
    v3[0]    = new T *  [nCell * nEqn];
    v3[0][0] = new T    [count];

    for (int i = 1; i < nCell; ++ i)
    {
        v3[i] = &v3[i-1][nEqn];
    }

    for (int j = 1; j < nEqn; ++ j)
    {
        v3[0][j] = &v3[0][j-1][nlen[0]];
    }

    for (int i = 1; i < nCell; ++ i)
    {
        v3[i][0] = &v3[i-1][nEqn-1][nlen[i-1]];
        for (int j = 1; j < nEqn; ++ j)
        {
            v3[i][j] = &v3[i][j-1][nlen[i]];
        }
    }

    return v3;
}

template < typename T >
T ** NewFastPointer(vector < T * > &vec)
{
    T **vv = new T * [vec.size()];
    for (int m = 0; m < vec.size(); ++ m)
    {
        vv[m] = vec[m];
    }
    return vv;
}

template < typename T >
void DelPointer(T *&v)
{
    if (!v) return;
    delete [] v;
    v = NULL;
}

template < typename T >
void DelPointer2(T **Data)
{
    if (!Data) return;
    delete [] Data[0];    Data[0] = NULL;
    delete [] Data;       Data = NULL;
}


template < typename T >
void DelPointer3(T ***v3)
{
    delete [] v3[0][0];    v3[0][0] = NULL;
    delete [] v3[0];       v3[0] = NULL;
    delete [] v3;          v3 = NULL;
}


template < typename T >
void DelFastPointer(T **v)
{
    delete [] v;    v = NULL;
}

template < typename T >
void FreePointer(T *&pointer)
{
    if (!pointer) return;
    delete pointer;    pointer = NULL;
}
}
