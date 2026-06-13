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
//! @file      PHSlice.h
//! @brief     Explain this file briefly.
//! @author    He Xin.
#pragma once

namespace PHSPACE
{
//! Forward declaration.
template<typename T, int N>
class PHArray;


class PHnilArraySection { };

template<typename T>
class PHArraySectionInfo
{
public:
    static const int isValidType = 0;
    static const int rank = 0;
    static const int isPick = 0;
};

template<>
class PHArraySectionInfo<Range> 
{
public:
    static const int isValidType = 1;
    static const int rank = 1;
    static const int isPick = 0;
};

template<>
class PHArraySectionInfo<int>
{
public:
    static const int isValidType = 1;
    static const int rank = 0;
    static const int isPick = 0;
};

template<>
class PHArraySectionInfo<PHnilArraySection>
{
public:
    static const int isValidType = 1;
    static const int rank = 0;
    static const int isPick = 0;
};

template<typename T,  typename T1,   typename T2 = PHnilArraySection, 
    typename T3 = PHnilArraySection, typename T4 = PHnilArraySection, 
    typename T5 = PHnilArraySection, typename T6 = PHnilArraySection, 
    typename T7 = PHnilArraySection, typename T8 = PHnilArraySection, 
    typename T9 = PHnilArraySection, typename T10 = PHnilArraySection, 
    typename T11 = PHnilArraySection>
class PHSliceInfo
{
public:
    static const int numValidTypes = PHArraySectionInfo<T1>::isValidType
                                   + PHArraySectionInfo<T2>::isValidType
                                   + PHArraySectionInfo<T3>::isValidType
                                   + PHArraySectionInfo<T4>::isValidType
                                   + PHArraySectionInfo<T5>::isValidType
                                   + PHArraySectionInfo<T6>::isValidType
                                   + PHArraySectionInfo<T7>::isValidType
                                   + PHArraySectionInfo<T8>::isValidType
                                   + PHArraySectionInfo<T9>::isValidType
                                   + PHArraySectionInfo<T10>::isValidType
                                   + PHArraySectionInfo<T11>::isValidType;

    static const int rank = PHArraySectionInfo<T1>::rank
                          + PHArraySectionInfo<T2>::rank
                          + PHArraySectionInfo<T3>::rank
                          + PHArraySectionInfo<T4>::rank
                          + PHArraySectionInfo<T5>::rank
                          + PHArraySectionInfo<T6>::rank
                          + PHArraySectionInfo<T7>::rank
                          + PHArraySectionInfo<T8>::rank
                          + PHArraySectionInfo<T9>::rank
                          + PHArraySectionInfo<T10>::rank
                          + PHArraySectionInfo<T11>::rank;

    static const int isPick = PHArraySectionInfo<T1>::isPick
                            + PHArraySectionInfo<T2>::isPick
                            + PHArraySectionInfo<T3>::isPick
                            + PHArraySectionInfo<T4>::isPick
                            + PHArraySectionInfo<T5>::isPick
                            + PHArraySectionInfo<T6>::isPick
                            + PHArraySectionInfo<T7>::isPick
                            + PHArraySectionInfo<T8>::isPick
                            + PHArraySectionInfo<T9>::isPick
                            + PHArraySectionInfo<T10>::isPick
                            + PHArraySectionInfo<T11>::isPick;

    typedef PHArray<T, numValidTypes> T_array;
    typedef PHArray<T, rank> T_slice;
};

}


