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
//! @file      PHHeader.h
//! @brief     Explain this file briefly.
//! @author    He Xin, Bell.

#pragma once
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <set>
#include <vector>
#include <algorithm>
#include "AMRDef.h"
#include "Precision.h"
using namespace std;

namespace PHSPACE
{

#define PHDEBUG
#ifdef PHDEBUG
#define ASSERT(x) \
    if (!(x))  \
{\
    cout << "ERROR: Assert " << #x << endl; \
    cout << "On line " << __LINE__ << endl; \
    cout << "On file " << __FILE__ << endl; \
    abort();                                \
}
#else
#define ASSERT(truth) (()0)
#endif

template < typename T >
void Set2Array(set<T> &data_set, T *data_array)
{
    typedef typename set<T>::iterator ITER;
    ITER iter;
    int icount = 0;
    for (iter = data_set.begin(); iter != data_set.end(); ++ iter)
    {
        data_array[icount++] = *iter;
    }
}

template < typename T >
void Array2Set(vector<T> &data_array, set<T> &data_set)
{
    for (std::size_t i = 0; i < data_array.size(); ++ i)
    {
        data_set.insert(data_array[i]);
    }
}

template < typename T >
void Array2Set(T *data_array, int size, set<T> &data_set)
{
    for (int i = 0; i < size; ++ i)
    {
        data_set.insert(data_array[i]);
    }
}

template < typename T1, typename T2 >
void GetPartOfArrayByCondition(vector<T1> &valueArray, const T1 &value, vector<T2> &referenceArray, vector<T2> &resultArray)
{
    resultArray.resize(0);
    for (std::size_t i = 0; i < valueArray.size(); ++ i)
    {
        if (valueArray[i] == value)
        {
            resultArray.push_back(referenceArray[i]);
        }
    }
}

template < typename T >
inline void SetField(vector< T > &field, const T &inputValue, int numberOfElements)
{
    for (int iElement = 0; iElement < numberOfElements; ++ iElement)
    {
        field[iElement] = inputValue;
    }
}

template < typename T >
inline void SetField(T *field, const T &inputValue, int numberOfElements)
{
    for (int iElement = 0; iElement < numberOfElements; ++ iElement)
    {
        field[iElement] = inputValue;
    }
}

template < typename T >
inline void SetField(T **field, int nDimension1, int nDimension2, const T &inputValue)
{
    for (int i = 0; i < nDimension1; ++ i)
    {
        for (int j = 0; j < nDimension2; ++ j)
        {
            field[i][j] = inputValue;
        }
    }
}

template < typename T1, typename T2 >
void SetField(PHVector1D <T1 *> &vectorA, const T2 &value, int numberOfElements)
{
    int numberOfEquations = vectorA.size();
    for (int iEquation = 0; iEquation < numberOfEquations; ++ iEquation)
    {
        for (int iElement = 0; iElement < numberOfElements; ++ iElement)
        {
            vectorA[iEquation][iElement] = value;
        }
    }
}

template < typename T1, typename T2 >
void SetField(PHVector1D < PHVector1D < T1 > > &vectorA, const T2 &value)
{
    int numberOfEquations = vectorA.size();

    for (int iEquation = 0; iEquation < numberOfEquations; ++ iEquation)
    {
        int numberOfElements = vectorA[iEquation].size();
        for (int iElement = 0; iElement < numberOfElements; ++ iElement)
        {
            vectorA[iEquation][iElement] = value;
        }
    }
}

template < typename T1, typename T2 >
void SetField(PHVector1D <T1 *> &vectorA, PHVector1D <T1 *> &vectorB, const T2 &coefficient, int numberOfElements)
{
    int numberOfEquations = vectorA.size();
    for (int iEquation = 0; iEquation < numberOfEquations; ++ iEquation)
    {
        for (int iElement = 0; iElement < numberOfElements; ++ iElement)
        {
            vectorA[iEquation][iElement] = coefficient * vectorB[iEquation][iElement];
        }
    }
}

template < typename T1, typename T2 >
inline void SetField(T1 *resultField, const T2 &inputValue, int numberOfElements)
{
    for (int iElement = 0; iElement < numberOfElements; ++ iElement)
    {
        resultField[iElement] = inputValue;
    }
}

template < typename T1, typename T2 >
inline void SetField(T1 *resultField, T2 *inputField, int numberOfElements)
{
    if (resultField == inputField) return;
    for (int iElement = 0; iElement < numberOfElements; ++ iElement)
    {
        resultField[iElement] = inputField[iElement];
    }
}

template < typename T1, typename T2 >
inline void SetField(PHVector1D < T1 > &resultField, const T2 &inputValue, int numberOfElements)
{
    if (resultField.size() == 0) return;
    SetField(&resultField[0], inputValue, numberOfElements);
}

template < typename T1, typename T2 >
void SetField(PHVector1D < PHVector1D < T1 > > &vectorA, PHVector1D < PHVector1D < T1 > > &vectorB, const T2 &coefficient)
{
    int numberOfEquations = vectorA.size();
    for (int iEquation = 0; iEquation < numberOfEquations; ++ iEquation)
    {
        int numberOfElements = vectorA[iEquation].size();
        for (int iElement = 0; iElement < numberOfElements; ++ iElement)
        {
            vectorA[iEquation][iElement] = coefficient * vectorB[iEquation][iElement];
        }
    }
}

template < typename T1, typename T2 >
inline void SetField(T1 *field, const T2 &inputvalue, int ist, int ied, int jst, int jed, int kst, int ked)
{
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                (*field)(i, j, k) = inputvalue;
            }
        }
    }
}

template < typename T >
void CopyField(T *It, const T *In, int size)
{
    if (!In || !It) return;

    for (int i = 0; i < size; ++ i)
    {
        It[i] = In[i];
    }
}

}