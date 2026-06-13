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
//! @file      vector3.hpp
//! @brief     The function of array operations.
//! @author    Huang Haibo, Zhang Jiaqian.

#ifndef VECTOR_H_FILE
#define VECTOR_H_FILE

template<class myType>
class vector3
{

public:
    vector3(myType x, myType y, myType z);
    vector3();
    myType dot(vector3<myType> other);
    vector3<myType> cross(vector3<myType> other);
    myType norm_square();
    myType x;
    myType y;
    myType z;
    myType get_x()    const
    {
        return x;
    }
    myType get_y()    const
    {
        return y;
    }
    myType get_z()    const
    {
        return z;
    }
    vector3<myType> operator- (const vector3<myType> &rhs) const
    {
        return vector3<myType>(x - rhs.x, y - rhs.y, z - rhs.z);
    }
    vector3<myType> operator+ (const vector3<myType> &rhs) const
    {
        return vector3<myType>(x + rhs.x, y + rhs.y, z + rhs.z);
    }
    vector3<myType> operator/ (const myType &rhs) const
    {
        return vector3<myType>(x/rhs, y/rhs, z/rhs);
    }
};

#endif
