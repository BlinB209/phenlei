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
//! @file      DataStruct_Array.h
//! @brief     PHArray defines a multidimensional array class. It is a template class.
//!            PHArray<T,N> refers to a N-dimensional T type array
//!            We can construct a multidimensional array by specify every dimension's length or range.
//!            The way of underlying data's storage follow the example of fortran,
//!            e.g. we can construct an 3d RDouble type array by length, PHArray<RDouble,3> array3d = PH<RDouble,3>(10, 20, 30)
//!            and access the i,j,k index data by operator(), i.e. array3d(i,j,k).
//!            The underlying data's memory storage sequence is from innermost(i) layer to outermost layer(k).
//!            So, if we loop over all the data, we recommend you loop from k to i, i.e for(k,...){for(j,...){for(i,...)}}}
//!            PHArray class also implement some other array's scalar operatation +=, -=, *=, /=, e.g. array1 *= 0.1.
//!            Some other set operator like =, +, -, *, /. We can directly use these operators to PHArray object
//!            e.g. dql(I, J, K, M) = MUSCLCoefXb * qr(I, J, K, M)
//! @author    He Xin, Zhang Jian, referred from Blitz.
#pragma once
#include <cstddef>
#include <stdexcept>
#include <stdlib.h>
#include "PHSimpleArray.h"
#include "PHStorage.h"
#include "PHMemblock.h"
#include "PHProduct.h"
#include "PHETBase.h"
#include "PHSlice.h"    //! Sub-arrays and slicing.

namespace PHSPACE
{

template < typename T, int N >
class PHFastArrayIterator;

template < typename T, int N >
class PHArray : public MemoryBlockReference<T> , public PHETBase< PHArray<T,N> >
{
private:
    typedef MemoryBlockReference<T> T_base;
    using T_base::data_;
    using T_base::changeToNullBlock;
    using T_base::numReferences;
    GeneralArrayStorage<N> storage_;
    SimpleArray<int, N> length_;
    SimpleArray<int, N> stride_;
    int zeroOffset_;
public:
    //! T_numtype is the numeric type stored in the array.
    typedef T T_numtype;

    //! T_index is a vector type which can be used to access elements
    //! of many-dimensional arrays.
    typedef SimpleArray<int, N> T_index;
    
    //! T_array is the array type itself -- Array<T_numtype, N_rank>.
    typedef PHArray<T, N> T_array;

    //! T_iterator is a a fast iterator for the array, used for expression.
    //! Templates.
    typedef PHFastArrayIterator<T, N> T_iterator;

    //! This array's dimension.
    static const int ranks = N;

public:
    //! Create a reference of another array.
    PHArray(const PHArray<T, N> &array);

    //! Any missing length arguments will have their value taken from the
    //! last argument.  For example,
    //! Array<int,3> A(32,64);
    //! will create a 32x64x64 array.  This is handled by setupStorage().
    //! Construct an array by each dimension's length.
    PHArray(GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! Construct an array by each dimension's length.
    explicit PHArray(int length0, GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! Construct an array by each dimension's length.
    PHArray(int length0, int length1, GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! Construct an array by each dimension's length.
    PHArray(int length0, int length1, int length2,
            GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! Construct an array by each dimension's length.
    PHArray(int length0, int length1, int length2, int length3,
            GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! Construct an array by each dimension's length.
    PHArray(int length0, int length1, int length2, int length3, int length4,
            GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! This constructor takes an extent (length) vector and storage format.
    PHArray(const SimpleArray<int, N> &extent, 
            GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! This constructor takes a vector of bases (lbounds) and a vector of extents.
    PHArray(const SimpleArray<int, N> &lbounds,
            const SimpleArray<int, N> &extent,
            const GeneralArrayStorage<N> &storage = GeneralArrayStorage<N>());

    //! These constructors allow arbitrary bases (starting indices) to be set.
    //! e.g. Array<int,2> A(Range(10,20), Range(20,30))
    //! will create an 11x11 array whose indices are 10..20 and 20..30
    //! Construct an array by each dimension's range.
    PHArray(Range r0, 
            GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! Construct an array by each dimension's range.
    PHArray(Range r0, Range r1,
            GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! Construct an array by each dimension's range.
    PHArray(Range r0, Range r1, Range r2,
            GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! Construct an array by each dimension's range.
    PHArray(Range r0, Range r1, Range r2, Range r3,
            GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! Construct an array by each dimension's range.
    PHArray(Range r0, Range r1, Range r2, Range r3, Range r4,
            GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! Construct an array by each dimension's range.
    PHArray(Range r0, Range r1, Range r2, Range r3, Range r4, Range r5,
            GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! Construct an array from an existing block of memory.
    PHArray(T *dataFirst,
            Range r0,
            preexistingMemoryPolicy deletionPolicy,
            GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! Construct an array from an existing block of memory.
    PHArray(T *dataFirst,
            Range r0, Range r1,
            preexistingMemoryPolicy deletionPolicy,
            GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! Construct an array from an existing block of memory.
    PHArray(T *dataFirst,
            Range r0, Range r1, Range r2,
            preexistingMemoryPolicy deletionPolicy,
            GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! Construct an array from an existing block of memory.
    PHArray(T *dataFirst,
            Range r0, Range r1, Range r2, Range r3,
            preexistingMemoryPolicy deletionPolicy,
            GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! Construct an array from an existing block of memory.
    PHArray(T *dataFirst,
            Range r0, Range r1, Range r2, Range r3, Range r4,
            preexistingMemoryPolicy deletionPolicy,
            GeneralArrayStorage<N> storage = GeneralArrayStorage<N>());

    //! These constructors construct an sub-array from another array.
    //! Construct an sub-array from another array.
    PHArray(PHArray<T, N> &array_in, Range r0);

    //! Construct an sub-array from another array.
    PHArray(PHArray<T, N> &array_in, Range r0, Range r1);

    //! Construct an sub-array from another array.
    PHArray(PHArray<T, N> &array_in, Range r0, Range r1, Range r2);

    //! Construct an sub-array from another array.
    PHArray(PHArray<T, N> &array_in, Range r0, Range r1, Range r2, Range r3);

    //! Construct an sub-array from another array.
    PHArray(PHArray<T, N> &array_in, Range r0, Range r1, Range r2, Range r3, Range r4);

    //! Construct an sub-array from another array.
    PHArray(PHArray<T, N> &array_in, Range r0, Range r1, Range r2, Range r3, Range r4, Range r5);

    //! Construct an sub-array from another array.
    PHArray(PHArray<T, N> &array_in, Range r0, Range r1, Range r2, Range r3, Range r4, Range r5, Range r6);

    //! This constructor is invoked by the operator()'s which take
    //! a combination of integer and Range arguments.  
    //! It's not intended for end-user use.
    template < int N2, typename gasConstant, typename R1, typename R2, typename R3, typename R4,
        typename R5, typename R6, typename R7, typename R8, typename R9, typename R10 >
        PHArray(PHArray<T, N2> &array_in, gasConstant r0, R1 r1, R2 r2,
        R3 r3, R4 r4, R5 r5, R6 r6, R7 r7, R8 r8, R9 r9, R10 r10);

public:
    //! operator() functions:
    //! Return constant reference to array[i0].
    const T &operator()(int i0) const;

    //! Return reference to array[i0].
    T & operator()(int i0);

    //! Return constant reference to array[i0][i1].
    const T & operator()(int i0, int i1) const;

    //! Return reference to array[i0][i1].
    T & operator()(int i0, int i1);

    //! Return constant reference to array[i0][i1][i2].
    const T & operator()(int i0, int i1, int i2) const;

    //! Return reference to array[i0][i1][i2].
    T & operator()(int i0, int i1, int i2);

    //! Return constant reference to array[i0][i1][i2][i3].
    const T & operator()(int i0, int i1, int i2, int i3) const;

    //! Return reference to array[i0][i1][i2][i3].
    T & operator()(int i0, int i1, int i2, int i3);

    //! Return constant reference to array[i0][i1][i2][i3][i4].
    const T & operator()(int i0, int i1, int i2, int i3, int i4) const;

    //! Return reference to array[i0][i1][i2][i3][i4].
    T & operator()(int i0, int i1, int i2, int i3, int i4);
    
    //! UGKS.
    //! Return constant reference to array[i0][i1][i2][i3][i4][i5].
    const T & operator()(int i0, int i1, int i2, int i3, int i4, int i5) const;

    //! UGKS.
    //! Return reference to array[i0][i1][i2][i3][i4][i5]
    T & operator()(int i0, int i1, int i2, int i3, int i4,int i5);

    //! Return constant reference to underlying data storage at index i, 
    //! ie. data_[i].
    const T & operator[](int i) const;

    //! Return reference to underlying data storage at index i, 
    //! ie. data_[i].
    T & operator[](int i);

    //! Return const reference to underlying data stored by an index array.
    //! e.g index = (10, 20), calculate offset = stride_[0]*10+stride_[1]*20
    //! then returns data_[offset].
    template<int N_rank2>
    const T & operator()(const SimpleArray<int, N_rank2> &index) const;

    //! Return reference to underlying data stored by an index array.
    //! e.g index = (10, 20), calculate offset = stride_[0]*10+stride_[1]*20
    //! then returns data_[offset].
    template<int N_rank2>
    T & operator()(const SimpleArray<int, N_rank2> &index);

    //! These constructors construct an sub-array from this array self.
    //! Construct an sub-array from this array self by offering a range.
    T_array operator()(Range r0) const;

    //! Construct an sub-array from this array self by offering a range.
    T_array operator()(Range r0, Range r1) const;

    //! Construct an sub-array from this array self by offering a range.
    T_array operator()(Range r0, Range r1, Range r2) const;

    //! Construct an sub-array from this array self by offering a range.
    T_array operator()(Range r0, Range r1, Range r2, Range r3) const;

    //! Construct an sub-array from this array self by offering a range.
    T_array operator()(Range r0, Range r1, Range r2, Range r3, Range r4) const;

    //! Construct an sub-array from this array self by offering a range.
    T_array operator()(Range r0, Range r1, Range r2, Range r3, Range r4, Range r5) const;

    //! Construct an sub-array from this array self by offering a range.
    T_array operator()(Range r0, Range r1, Range r2, Range r3, Range r4, Range r5, Range r6) const;

    //! After calling slice(int rank, Range r), the array refers only to the
    //! range r of the original array.
    //! e.g. Array<int,1> x(100);
    //!      x.slice(firstRank, Range(25,50));
    //!      x = 0;       //! Sets elements 25..50 of the original array to 0.
    template < typename T1, typename T2 >
    typename PHSliceInfo<T, T1, T2>::T_slice operator()(T1 r1, T2 r2) const;

    template < typename T1, typename T2, typename T3 >
    typename PHSliceInfo<T, T1, T2, T3>::T_slice operator()(T1 r1, T2 r2, T3 r3) const;

    template < typename T1, typename T2, typename T3, typename T4 >
    typename PHSliceInfo<T, T1, T2, T3, T4>::T_slice operator()(T1 r1, T2 r2, T3 r3, T4 r4) const;

    template < typename T1, typename T2, typename T3, typename T4, typename T5 >
    typename PHSliceInfo<T, T1, T2, T3, T4, T5>::T_slice operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5) const;

    template < typename T1, typename T2, typename T3, typename T4, typename T5, typename T6 >
    typename PHSliceInfo<T, T1, T2, T3, T4, T5, T6>::T_slice operator()(T1 r1, T2 r2, T3 r3, T4 r4, T5 r5, T6 r6) const;

    //! Array expression operands.

    //! Scalar operand.
    //! NEEDS_WORK : need a precondition check on.
    //! isStorageContiguous when operator, is used.
    T_array & operator=(T x);

    //! operator = expr
    template<typename T_expr>
    T_array & operator = (const PHETBase<T_expr> &expr);

    //! operator = PHArray<T,N>
    T_array & operator = (const PHArray<T, N> &rhs);

    //! operator +=
    template<typename T_expr> T_array & operator += (const T_expr &expr);

    //! operator -=
    template<typename T_expr> T_array & operator -= (const T_expr &expr);

    //! operator *=
    template<typename T_expr> T_array & operator *= (const T_expr &expr);

    //! operator /=
    template<typename T_expr> T_array & operator /= (const T_expr &expr);

    //! Member functions.

    //! Return index relative to data_, ie. data_[index].
    //! index =  i0 * stride_[0] + i1 * stride_[1] + i2 * stride_[2] + i3 * stride_[3].
    int getindex(int i0, int i1, int i2, int i3) const;

    //! Return the rank dimension's lower bound, ie. range's first index, equals base_[rank].
    int lbound(int rank) const;

    //! Return the rank dimension's upper bound, ie. range's last index, equals base_[rank]+length_[rank]-1.
    int ubound(int rank) const;

    //! Return the array storing lbound, ie. base_.
    SimpleArray<int, N> lbound() const;

    //! Return the rank dimension's length, ig length_[rank]
    int length(int rank) const;

    //! Return the array storing length_.
    const SimpleArray<int, N> & length() const;

    //! Return the array storing base_.
    //! base_ is each dimension's lower bound.
    const SimpleArray<int, N> & base() const;

    //! Return the rank dimension's base, ie. base_[rank].
    int base(int rank) const;

    //! Return the rank dimension's stride, ie. stride_[rank].
    int stride(int rank) const;

    //! Total number of elements in this array.
    int numElements() const;

    //! data_ always refers to the point (0,0,...,0) which may
    //! not be in the array if the base is not zero in each rank.
    //! These data() routines return a pointer to the first
    //! element in the array (but note that it may not be
    //! stored first in memory if some ranks are stored descending).
    const T * data() const;

    //! data() is the location of the first data member,that is the location of the base.
    //! The base array starting subscript of each dimension,ri.first().
    //! It is especially important for slice.
    //! After slice, the first data point actually is (ifirst, jfirst, ...).
    //! The actual location is data_ + dataOffset().
    T * data();

    //! Judge whether the array's outer dimension's stride equals sum of inner dimension's stride.
    //! eg. if outerRank =3 and innerRank =2 , returns true if stride[3] == stride[2]*length[2], otherwise, returns false
    bool canCollapse(int outerRank, int innerRank) const;

    //! Return the PHFastArrayIterator of this array.
    //! It's not intended for end-user use.
    T_iterator beginFast() const;

    //! Return storage_.ordering(storageRankIndex).
    //! It's not intended for end-user use.
    int ordering(int storageRankIndex) const;

    //! Return storage_.ordering().
    //! It's not intended for end-user use.
    const SimpleArray<int, N> & ordering() const;

    //! True if storage_.ordering(rank) == 0.
    //! It's not intended for end-user use.
    bool isMajorRank(int rank) const;
    
    //! True if storage_.ordering(rank) != 0.
    //! It's not intended for end-user use.
    bool isMinorRank(int rank) const;

    //! Call storage_.isRankStoredAscending(rank).
    //! It's not intended for end-user use.
    bool isRankStoredAscending(int rank) const;
    //! The underlying data storage(data_) always refers to the point index at (0,0,...,0) which may
    //! not be in the array if the base is not zero in each rank.
    //! Return the pointer to data.
    T * getData() const { return data_; }

protected:
    //! Implementation routines.
    void computeStrides();

    void setupStorage(int rank);

    void setupStorage(int rank, T *data, preexistingMemoryPolicy deletionPolicy);

    void calculateZeroOffset();

    void Reference(const T_array&);

    T_array & noConst() const;

    void constructSubarray(PHArray<T, N> &array_in, Range r0);

    void constructSubarray(PHArray<T, N> &array_in, Range r0, Range r1);

    void constructSubarray(PHArray<T, N> &array_in, Range r0, Range r1, Range r2);

    void constructSubarray(PHArray<T, N> &array_in, Range r0, Range r1, Range r2, Range r3);

    void constructSubarray(PHArray<T, N> &array_in, Range r0, Range r1, Range r2, Range r3, Range r4);

    void constructSubarray(PHArray<T, N> &array_in, Range r0, Range r1, Range r2, Range r3, Range r4, Range r5);

    void constructSubarray(PHArray<T, N> &array_in, Range r0, Range r1, Range r2, Range r3, Range r4, Range r5, Range r6);

    template < int N2, typename gasConstant,
        typename R1, typename R2, typename R3, typename R4, typename R5, typename R6, typename R7,
        typename R8, typename R9, typename R10 >
        void constructSlice(PHArray<T, N2> &array_in,
        gasConstant r0, R1 r1, R2 r2, R3 r3, R4 r4, R5 r5, R6 r6, R7 r7, R8 r8, R9 r9, R10 r10);
    void slice(int rank, Range r);

    template<int N2>
    void slice(int &setRank, Range r, PHArray<T, N2> &array_in, SimpleArray<int, N2> &rankMap, int sourceRank);

    template<int N2>
    void slice(int &setRank, int i, PHArray<T, N2> &array_in, SimpleArray<int, N2> &rankMap, int sourceRank);

    template<int N2>
    void slice(int &setRank, PHnilArraySection nil, PHArray<T, N2> &array_in, SimpleArray<int, N2> &rankMap, int sourceRank);

    T_array & initialize(T x);

    int dataOffset() const;

    template<typename T_expr, typename T_update>
    T_array & evaluate(T_expr expr, T_update);

    template<typename T_expr, typename T_update>
    T_array & evaluateWithIndexTraversal1(T_expr expr, T_update);

    template<typename T_expr, typename T_update>
    T_array & evaluateWithIndexTraversalN(T_expr expr, T_update);

    template<typename T_expr, typename T_update>
    T_array & evaluateWithStackTraversal1(T_expr expr, T_update);

    template<typename T_expr, typename T_update>
    T_array & evaluateWithStackTraversalN(T_expr expr, T_update);

    //! #ifdef USE_SMARTARRY will construct an multi-array pointer to the underlying data_
    //! e.g. for a three dimensional T type array with range(r0, r1, r2).
    //! We can get and T*** type pointer called datap3 by GetSmartArrayP3() function
    //! and access the data at index i, j, k directly by datap3[k][j][i].
    //! Remember the way of array storage of range r0, r1, ... corresponding with innermost layer to outermost layer.
    //! Moreover, if use operator(i, j, k), will return datap3[k][j][i] instead of data_[i0 * stride_[0] + i1 * stride_[1]+ i2 * stride_[2]].
#ifdef USE_SMARTARRAY
public:
    typedef T      * TypeP1; //TypeP1 T *
    typedef TypeP1 * TypeP2; //TypeP2 T **
    typedef TypeP2 * TypeP3; //TypeP3 T ***
    typedef TypeP3 * TypeP4; //TypeP4 T ****
    typedef TypeP4 * TypeP5; //TypeP5 T *****
    typedef TypeP5 * TypeP6; //TypeP6 T ******
private:
    TypeP1 data1, datap1; //TypeP1 T *
    TypeP2 data2, datap2; //TypeP2 T **
    TypeP3 data3, datap3; //TypeP3 T ***
    TypeP4 data4, datap4; //TypeP4 T ****
    TypeP5 data5, datap5; //TypeP5 T *****
    TypeP6 data6, datap6; //TypeP6 T ******
private:
    void ZeroAllPointer();
    void AllocateRange(int dimension);
    void AllocateRange1();
    void AllocateRange2();
    void AllocateRange3();
    void AllocateRange4();
    void AllocateRange5();
    void AllocateRange6();
public:
    ~PHArray()
    {
        if (data1 != NULL) {                     data1 = NULL; }
        if (data2 != NULL) { delete [] data2;    data2 = NULL; }
        if (data3 != NULL) { delete [] data3;    data3 = NULL; }
        if (data4 != NULL) { delete [] data4;    data4 = NULL; }
        if (data5 != NULL) { delete [] data5;    data5 = NULL; }
        if (data6 != NULL) { delete [] data6;    data6 = NULL; }
    }

    TypeP1 GetSmartArrayP1() const
    {
        return datap1;
    }

    TypeP2 GetSmartArrayP2() const
    {
        return datap2;
    }

    TypeP3 GetSmartArrayP3() const
    {
        return datap3;
    }

    TypeP4 GetSmartArrayP4() const
    {
        return datap4;
    }

    TypeP5 GetSmartArrayP5() const
    {
        return datap5;
    }

    TypeP6 GetSmartArrayP6() const
    {
        return datap6;
    }
#endif
    //! End
};

}

const int firstRank    = 0;
const int secondRank   = 1;
const int thirdRank    = 2;
const int fourthRank   = 3;
const int fifthRank    = 4;
const int sixthRank    = 5;
const int seventhRank  = 6;
const int eighthRank   = 7;
const int ninthRank    = 8;
const int tenthRank    = 9;
const int eleventhRank = 10;

const int firstDim    = 0;
const int secondDim   = 1;
const int thirdDim    = 2;
const int fourthDim   = 3;
const int fifthDim    = 4;
const int sixthDim    = 5;
const int seventhDim  = 6;
const int eighthDim   = 7;
const int ninthDim    = 8;
const int tenthDim    = 9;
const int eleventhDim = 10;

#include "PHFastArrayIterator.h"    //! Fast Array iterators (for et).
#include "PHExpr.h"                 //! Array expression objects.
#include "PHExprWrap.h"
#include "PHET.h"                   //! Expression templates.
#include "PHOps.hxx"                //! Assignment operators.
#include "PHArray.hxx"
#include "PHSlice.hxx"              //! Slicing and sub-arrays.
