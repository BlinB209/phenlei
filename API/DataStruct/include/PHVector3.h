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
//! @file      PHVector3D.h
//! @brief     Explain this file briefly.
//! @author    Wan YunBo.

#pragma once
#include <math.h>
#include <PHHeader.h>

using namespace std;

namespace PHSPACE
{
template <class T>
struct Vector3
{
    Vector3 (T x_=0, T y_=0, T z_=0) : x(x_), y(y_), z(z_) {}
 
    void set (T x_, T y_, T z_)
    { 
        x = x_;
        y = y_;
        z = z_;
    }
 
    Vector3 normalize() const
    {
        return ((*this) / norm());
    }
    RDouble norm () const 
    {
        return sqrt(normSquare());
    }
    T normSquare () const
    {
        return x*x + y*y + z*z;
    }
 
    bool operator == (const Vector3 &v) const
    {
        return x == v.x && y == v.y && z == v.z;
    }
    bool operator != (const Vector3 &v) const
    {
        return x != v.x || y != v.y || z != v.z;
    }
 
    Vector3 operator + (const Vector3 &v) const
    {
        return Vector3(x+v.x, y+v.y, z+v.z);
    }
    Vector3 & operator += (const Vector3 &v)
    {
        x+=v.x; y+=v.y; z+=v.z; return *this;
    }
    Vector3 operator - () const
    {
        return Vector3(-x, -y, -z);
    }
    Vector3 operator - (const Vector3 &v) const
    {
        return Vector3(x-v.x, y-v.y, z-v.z);
    }
    Vector3 & operator -= (const Vector3 &v)
    { 
        x-=v.x; y-=v.y; z-=v.z; return *this;
    }
    Vector3 operator * (T s) const
    {
        return Vector3(x*s, y*s, z*s);
    }
    T  operator * (const Vector3 &v) const
    { 
        return x*v.x + y*v.y + z*v.z;
    }
    //Vector3 operator * (const Vector3 &v) const
    //{
    //    return Vector3(*v.x, *v.y, *v.z);
    //}
    Vector3 & operator *= (float s)
    { 
        x*=s; y*=s; z*=s;
        return *this;
    }
    Vector3 operator / (float s) const
    { 
        ASSERT(s);
        return (*this)*(1/s);
    }
    Vector3& operator *= (RDouble s)
    {
        x *= s; y *= s; z *= s;
        return *this;
    }
    Vector3 operator / (RDouble s) const
    {
        ASSERT(s);
        return (*this) * (1 / s);
    }
    Vector3 & operator /= (float s)
    { 
        ASSERT(s); 
        return (*this)*=(1/s);
    }
    Vector3& operator /= (RDouble s) const
    {
        ASSERT(s);
        return (*this) *= (1 / s);
    }
    Vector3& operator /= (RDouble s)
    {
        ASSERT(s);
        return (*this) *= (1 / s);
    }
    T operator[](int i)
    {  
        if (i == 0) return this->x;
        if (i == 1) return this->y;
        if (i == 2) return this->z;
     }
    T operator[](int i) const // 声明为 const
    {
        if (i == 0) return this->x;
        if (i == 1) return this->y;
        if (i == 2) return this->z;
        // 可以添加越界检查
        //return T(); // 返回默认值
    }
    
    T x, y, z;
};
 
template <class T> inline
T Dot (const Vector3<T> &l, const Vector3<T> &r)
{
    return l.x*r.x + l.y*r.y + l.z*r.z;
}
 
template <class T> inline
Vector3<T> Cross (const Vector3<T> &l, const Vector3<T> &r)
{
    return Vector3<T>(l.y*r.z - l.z*r.y, l.z*r.x - l.x*r.z, l.x*r.y - l.y*r.x);
}
 
template <class T> inline
T BlendProduct (const Vector3<T> &l, const Vector3<T> &m, const Vector3<T> &r)
{
    return Dot(Cross(l, m), r);
}
 
typedef Vector3<char>    Vector3C;
typedef Vector3<int>     Vector3I;
typedef Vector3<float>   Vector3F;
typedef Vector3<RDouble> Vector3D;

}
