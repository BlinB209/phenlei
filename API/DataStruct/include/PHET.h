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
//! @file      PHET.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once

namespace PHSPACE
{
class PHUpdaterBase { };

template<typename X, typename Y>
class PHUpdate : public PHUpdaterBase
{
public:
    static inline void update(X &x, Y y)
    {
        x = (X)y; 
    }
};

#define PH_DEFINE_UPDATER(name, op)                         \
template<typename X, typename Y>                            \
class name : public PHUpdaterBase                           \
{                                                           \
public:                                                     \
    static inline void update(X &x, Y y)                    \
    {                                                       \
        x op y;                                             \
    }                                                       \
}

PH_DEFINE_UPDATER(fy_plus_update    , +=);
PH_DEFINE_UPDATER(fy_minus_update   , -=);
PH_DEFINE_UPDATER(fy_multiply_update, *=);
PH_DEFINE_UPDATER(fy_divide_update  , /=);

#define PH_DECLARE_ARRAY_UNARY(name, functor)                             \
                                                                          \
template <typename T1>                                                    \
inline                                                                    \
typename PHUnaryExprResult<functor, T1>::T_result                         \
name(const PHETBase<T1> &d1)                                              \
{                                                                         \
    typedef typename PHUnaryExprResult<functor, T1>::T_result result;     \
    return result(PHExprWrap<T1>::getExpr(d1.unwrap()));                  \
}

#define PH_DECLARE_ARRAY_BINARY(name, applic)                             \
                                                                          \
template <typename T1, typename T2>                                       \
inline                                                                    \
typename PHBinaryExprResult<applic, T1, T2>::T_result                     \
name(const PHETBase<T1> &d1,const PHETBase<T2> &d2)                       \
{                                                                         \
    typedef typename PHBinaryExprResult<applic, T1, T2>::T_result result; \
    return result(PHExprWrap<T1>::getExpr(d1.unwrap()),                   \
                  PHExprWrap<T2>::getExpr(d2.unwrap()));                  \
}

#define PH_DECLARE_ARRAY_BINARY_SCALAR(name, applic, scalar_type)                 \
template<typename T>                                                              \
inline                                                                            \
typename PHBinaryExprResult<applic, scalar_type,T>::T_result                      \
name(const scalar_type d1, const PHETBase<T> &d2)                                 \
{                                                                                 \
    typedef typename PHBinaryExprResult<applic, scalar_type, T>::T_result result; \
    return result(PHExprWrap<scalar_type >::getExpr(d1),                          \
                  PHExprWrap<T>::getExpr(d2.unwrap()));                           \
}                                                                                 \
                                                                                  \
template<typename T>                                                              \
inline                                                                            \
typename PHBinaryExprResult<applic, T, scalar_type >::T_result                    \
name(const PHETBase<T> &d1, const scalar_type d2)                                 \
{                                                                                 \
    typedef typename PHBinaryExprResult<applic, T, scalar_type >::T_result result;\
    return result(PHExprWrap<T>::getExpr(d1.unwrap()),                            \
                  PHExprWrap<scalar_type >::getExpr(d2));                         \
}

#define PH_DECLARE_ARRAY_SCALAR_OPERATIONS(scalar_type)                           \
                                                                                  \
PH_DECLARE_ARRAY_BINARY_SCALAR(operator +, Add     , scalar_type)                 \
PH_DECLARE_ARRAY_BINARY_SCALAR(operator -, Subtract, scalar_type)                 \
PH_DECLARE_ARRAY_BINARY_SCALAR(operator *, Multiply, scalar_type)                 \
PH_DECLARE_ARRAY_BINARY_SCALAR(operator /, Divide  , scalar_type)

PH_DECLARE_ARRAY_SCALAR_OPERATIONS(char)
PH_DECLARE_ARRAY_SCALAR_OPERATIONS(unsigned char)
PH_DECLARE_ARRAY_SCALAR_OPERATIONS(short)
PH_DECLARE_ARRAY_SCALAR_OPERATIONS(unsigned short)
PH_DECLARE_ARRAY_SCALAR_OPERATIONS(int)
PH_DECLARE_ARRAY_SCALAR_OPERATIONS(unsigned int)
PH_DECLARE_ARRAY_SCALAR_OPERATIONS(long)
PH_DECLARE_ARRAY_SCALAR_OPERATIONS(unsigned long)
PH_DECLARE_ARRAY_SCALAR_OPERATIONS(float)
PH_DECLARE_ARRAY_SCALAR_OPERATIONS(double)
PH_DECLARE_ARRAY_SCALAR_OPERATIONS(long double)

PH_DECLARE_ARRAY_BINARY(operator+, Add)
PH_DECLARE_ARRAY_BINARY(operator-, Subtract)
PH_DECLARE_ARRAY_BINARY(operator*, Multiply)
PH_DECLARE_ARRAY_BINARY(operator/, Divide)
// unary operators
PH_DECLARE_ARRAY_UNARY(operator+, UnaryPlus)
PH_DECLARE_ARRAY_UNARY(operator-, UnaryMinus)

}

