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
//! @file      Data_Util.h
//! @brief     This file defines some util function about copying data.
//! @author    He Xin, Zhang Jian.

#pragma once

namespace PHSPACE
{

struct Data_Util
{
    //! Function template, copy T type data from source to target.
    //! The target is a null pointer at first, this function will allocate memory internally.
    //! @param[in] size      size of the source object.
    //! @param[in] source    the source address.
    //! @param[in] target    the target address.
    template < typename T >
    static void copy_void(void *&target, void *source, int size)
    {
        target = new T[size];
        for (int i = 0; i < size; ++ i)
        {
            reinterpret_cast<T *>(target)[i] = reinterpret_cast<T *>(source)[i];
        }
    }

    //! Free memory for T type data.
    //! @param[in] data    the object to be deleted.
    template < typename T >
    static void delete_void(void *data)
    {
        delete [] reinterpret_cast<T *>(data);
    }

    //! Function template, copy T type data from source to target with same type and size.
    //! @param[in] size      size of the source object.
    //! @param[in] source    the source address.
    //! @param[in] target    the target address.
    template < typename T >
    static void copy_void_exist(void *target, void *source, int size)
    {
        for (int i = 0; i < size; ++ i)
        {
            reinterpret_cast<T *>(target)[i] = reinterpret_cast<T *>(source)[i];
        }
    }

    //! Copy data from source to target with the assumption that target hasn't allocated memory.
    //! This function will help to allocate memory internally.
    //! @param[out] target    The target address.
    //! @param[in]  source    The source address.
    //! @param[in]  type      The data's type.
    //! @param[in]  size      The data's size.
    static void copy_data_sub(void *&target, void *source, int type, int size);

    //! Copy data from source to target with the assumption that target has allocated memory.
    //! This function will replace the old data with new data.
    //! Users must guarantee the source has the same size and type with target.
    //! @param[out] target    The target address.
    //! @param[in]  source    The source address.
    //! @param[in]  type      The data's type.
    //! @param[in]  size      The data's size.
    static void copy_data_exist(void *target, void *source, int type, int size);
};

}