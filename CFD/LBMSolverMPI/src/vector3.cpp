#include "vector3.hpp"

template<class myType>    //! 函数模板 function templates,任意类型的数据操作
vector3<myType>::vector3(myType x, myType y, myType z) : x(x), y(y), z(z)
{

}

template<class myType>
myType vector3<myType>::dot(vector3<myType> other)
{
    return (x * other.get_x() +
        y * other.get_y() +
        z * other.get_z());
}

template<class myType>
vector3<myType> vector3<myType>::cross(vector3<myType> other)
{
    return vector3<myType>(
        (y * other.get_z()) - (z * other.get_y()),
        (z * other.get_x()) - (x * other.get_z()),
        (x * other.get_y()) - (y * other.get_x()));
}

template<class myType>
vector3<myType>::vector3() : x(0), y(0), z(0)
{

}

template<class myType>
myType vector3<myType>::norm_square()
{
    return (x * x + y * y + z * z);
}
