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
//! @file      Math_Limiter.h
//! @brief     Limiters' mathematic functions, such as vencat, vanleer, etc.
//! @author    Bell, He Xin.

#pragma once
#include "Precision.h"
#include "Math_BasisFunction.h"
#include <cmath>
#include <algorithm>
using namespace std;
 
namespace PHSPACE
{
//! CheckLimit.
//! @param[in]: value
//! @param[in]: v_min
//! @param[in]: v_max
inline bool CheckLimit(const RDouble &value, const RDouble &v_min, const RDouble &v_max)
{
    if (value < v_min || value > v_max) return false;
    return true;
}

//! NOLIMITER.
//! @param[in]: x
//! @param[in]: y
inline RDouble NOLIMITER(const RDouble &x)
{
    return x;
}

//! MINMOD.
//! @param[in]: x
//! @param[in]: y
inline RDouble MINMOD(const RDouble &x, const RDouble &y)
{
    return half * (SIGN(one, x) + SIGN(one, y)) * MIN(ABS(x), ABS(y));
}

//! VANLEER.
//! @param[in]: x
//! @param[in]: y
inline RDouble VANLEER(const RDouble &x, const RDouble &y)
{
    const RDouble eps = 1.0e-12;

    return (SIGN(one, x) + SIGN(one, y)) * x * y / (ABS(x + y) + eps);
}

//! MIX_MINMOD_VANLEER.
//! @param[in]: x
//! @param[in]: y
inline RDouble MIX_MINMOD_VANLEER(const RDouble &x, const RDouble &y)
{
    const RDouble eps = 1.0e-15;
    RDouble z, z2, z3, xy2, absx, absy, min_van;

    z = half * (SIGN(one, x) + SIGN(one, y));
    absx = abs(x);
    absy = abs(y);
    xy2 = (x + x) * y;
    z2 = xy2 / (x * x + y * y + eps);
    z3 = xy2 / (absx + absy + eps);
    min_van = z * ((one - z2) * min(absx , absy) + z2 * z3);
/*
    if (ABS(z) < eps)
    {
        min_van = 0.0;
    }
    else
    {
        absx = abs(x);
        absy = abs(y);
        xy2 = (x + x) * y;
        z2 = xy2 / (x * x + y * y + eps);
        z3 = xy2 / (absx + absy + eps);
        min_van = z * ((1.0 - z2) * min(absx, absy) + z2 * z3);
    }
*/
    return min_van;
}

//! VANALBADA.
//! @param[in]: x
//! @param[in]: y
inline RDouble VANALBADA(const RDouble &x, const RDouble &y)
{
    //const RDouble eps = 1.0e-7;
    //RDouble cx, cy, cd;
    //cx = x * x + eps;
    //cy = y * y + eps;
    //cd = 1.0 / (cx + cy);
    //return (cx * y + cy * x) * cd;
    const RDouble eps = 1.0e-6;
//    RDouble cx, cy, cd;
//    cx = x * x + eps;
//    cy = y * y + eps;
//    cd = 1.0 / (cx + cy);
    return 0.5 * (SIGN(one, x) + SIGN(one, y)) * ((x * y)* ABS(x + y)) / (x * x + y * y + eps);
}

//! VANALBADACLZ.
//! @param[in]: x
//! @param[in]: y
inline RDouble VANALBADACLZ(const RDouble &x, const RDouble &y)
{
    const RDouble eps = 1.0e-6;
    RDouble cx, cy, cd;
    cx = x * x + eps;
    cy = y * y + eps;
    cd = 1.0 / (cx + cy);
    return 0.5 * (SIGN(one, x) + SIGN(one, y)) * (cx * y + cy * x) * cd;
}

//! SMOOTH.
//! @param[in]: x
//! @param[in]: y
inline RDouble SMOOTH(const RDouble &x, const RDouble &y)
{
    const RDouble eps = 1.0e-6;
    RDouble Smooth= (2. * x * y + eps) / (x * x + y * y + eps);    //0.5 * (SIGN(1.0, x) + SIGN(1.0, y)) *
    return Smooth;
}

//! SMOOTHTURB.
//! @param[in]: x
//! @param[in]: y
inline RDouble SMOOTHTURB(const RDouble &x, const RDouble &y)
{
    const RDouble eps = 1.0e-6;
    RDouble Smooth = (2. * x * y) / (x * x + y * y + eps);    //0.5 * (SIGN(1.0, x) + SIGN(1.0, y)) *
    return Smooth;
}

//! THIRD_ORDER_SMOOTH.
//! @param[in]: x
//! @param[in]: y
inline RDouble THIRD_ORDER_SMOOTH(const RDouble &x, const RDouble &y)
{
    const RDouble eps2 = 1.0e-6;
    RDouble Smooth1 = x * (y * y + 2.0 * eps2) +  y * (2.0 * x * x + eps2);
    RDouble Smooth2 = 2.0 * x * x - x * y + 2.0 * y * y + 3.0 * eps2;
    RDouble Smooth = Smooth1 / Smooth2;
    return Smooth;
}

//! THIRD_ORDER_MINMOD_SMOOTH.
//! @param[in]: x
//! @param[in]: y
inline RDouble THIRD_ORDER_MINMOD_SMOOTH(const RDouble &x, const RDouble &y)
{
    const RDouble eps2 = 1.0e-6;
    RDouble Smooth1 = x * (y * y + 2.0 * eps2) +  y * (2.0 * x * x + eps2);
    RDouble Smooth2 = 2.0 * x * x - x * y + 2.0 * y * y + 3.0 * eps2;
    RDouble Smooth = 0.5 * (SIGN(one, x) + SIGN(one, y)) * Smooth1 / Smooth2;
    return Smooth;
}

//! VenFun.
//! @param[in]: x
//! @param[in]: y
//! @param[in]: eps
inline RDouble VenFun(const RDouble &x, const RDouble &y, const RDouble &eps)
{
    RDouble x2 = x * x;
    RDouble xy = x * y;
    RDouble y2 = y * y;
    return ((x2 + eps + 2.0 * xy) / (x2 + 2.0 * y2 + xy + eps));
}

PH_DEFINE_BINARY_FUNCTION(PH_MinMod, MINMOD)
PH_DEFINE_BINARY_FUNCTION(PH_Vanleer, VANLEER)
PH_DEFINE_BINARY_FUNCTION(PH_Vanalbada, VANALBADA)
PH_DEFINE_BINARY_FUNCTION(PH_VanalbadaCLZ, VANALBADACLZ)
PH_DEFINE_BINARY_FUNCTION(PH_Smooth, SMOOTH)
PH_DEFINE_BINARY_FUNCTION(PH_Third_Order_Smooth, THIRD_ORDER_SMOOTH)
PH_DEFINE_BINARY_FUNCTION(PH_Third_Order_Minmod_Smooth, THIRD_ORDER_MINMOD_SMOOTH)
PH_DEFINE_BINARY_FUNCTION(PH_Smooth_turb, SMOOTHTURB)
PH_DEFINE_BINARY_FUNCTION(PH_Min_Van, MIX_MINMOD_VANLEER)

PH_DECLARE_ARRAY_BINARY(MinMod, PH_MinMod)
PH_DECLARE_ARRAY_BINARY(Vanleer, PH_Vanleer)
PH_DECLARE_ARRAY_BINARY(Vanalbada, PH_Vanalbada)
PH_DECLARE_ARRAY_BINARY(VanalbadaCLZ, PH_VanalbadaCLZ)
PH_DECLARE_ARRAY_BINARY(Smooth, PH_Smooth)
PH_DECLARE_ARRAY_BINARY(Third_Order_SMOOTH, PH_Third_Order_Smooth)
PH_DECLARE_ARRAY_BINARY(Third_Order_Minmod_Smooth, PH_Third_Order_Minmod_Smooth)
PH_DECLARE_ARRAY_BINARY(Smoothturb, PH_Smooth_turb)

PH_DECLARE_ARRAY_BINARY(Minvan, PH_Min_Van)

}