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
//! @file      Math_BasisFunction.h
//! @brief     Basis mathematic function, such as min/max/sign, etc.
//! @author    Bell, He Xin.

#pragma once
#include <cmath>
#include <float.h>
#include <algorithm>
#include "Precision.h"
#include "TypeDefine.h"
using namespace std;

namespace PHSPACE
{
//! Get the interval: size of total is 'nTotal', it's divided into 'nInterval'.
int GetProgressInterval(int nTotal, int nInterval);

//! Return the minimum value of the two input values.
template < typename T >
inline T MIN(const T &a, const T &b)
{
    return (a < b) ? a : b;
}

template < typename T1, typename T2 >
inline T1 MIN(const T1 &a, const T2 &b)
{
    return static_cast<T1>((a < b) ? a : b);
}

//! Return the maximum value of the two input values.
template < typename T >
inline T MAX(const T &a, const T &b)
{
    return (a > b) ? a : b;
}

template < typename T1, typename T2 >
inline T1 MAX(const T1 &a, const T2 &b)
{
    return static_cast<T1>((a > b) ? a : b);
}

//! Return the absolute value of the input value.
template < typename T >
inline T ABS(const T &a)
{
    return (a < 0) ? -a : a;
}

//! Return the square value of the input value.
template < typename T >
inline T SQR(const T &a)
{
    return a * a;
}

//! Return the cube value of the input value.
template < typename T >
inline T POWER3(const T &a)
{
    return a * a * a;
}

//! Return the nearest round number value of the input value.
template < typename T >
inline T ROUND(const T &a)
{
    return (a > 0.0) ? floor(a + 0.5):ceil(a - 0.5);
}

//! Return the square distance of the three input values.
template < typename T >
inline T SQR(const T &a, const T &b, const T &c)
{
    return a * a + b * b + c * c;
}

inline RDouble SQR(RDouble *point1, RDouble *point2)
{
    return  (point1[0] - point2[0]) * (point1[0] - point2[0]) +
            (point1[1] - point2[1]) * (point1[1] - point2[1]) +
            (point1[2] - point2[2]) * (point1[2] - point2[2]);
}

//! Return the distance of the three input values.
template < typename T >
inline T DISTANCE(const T &a, const T &b, const T &c)
{
    return sqrt(SQR(a, b, c));
}

inline RDouble DISTANCE(RDouble *point1, RDouble *point2)
{
    return sqrt((point1[0] - point2[0]) * (point1[0] - point2[0]) +
                (point1[1] - point2[1]) * (point1[1] - point2[1]) +
                (point1[2] - point2[2]) * (point1[2] - point2[2]));
}

inline RDouble DistSQ(RDouble *a1, RDouble *a2, int dim)
{
    RDouble distSQ = 0.0;
    while (--dim >= 0)
    {
        RDouble diff = a1[dim] - a2[dim];
        distSQ += diff * diff;
    }
    return distSQ;
}

//! Swap the two input values.
template < typename T >
inline void SWAP(T &a, T &b)
{
    T c = a;
    a = b;
    b = c;
}

template < typename T >
bool isErrorData(T data)
{
#ifdef WIN32
    if (!_isnan(data) || ABS(data) < SMALL)
    {
        return false;
    }
#else
    if (ABS(data) < LARGE && (std::isnormal(data) || ABS(data) < SMALL))
    {
        return false;
    }
#endif
    return true;
}

//! sign(a, b): return |a| * sign(b). a is the same as b.
template < typename T >
inline T SIGN(const T &a, const T &b)
{
    int s = 1;
    if (b < 0) s = -1;
    return s * ABS(a);
}

//! sign(a, b): return |a| * sign(b). a and b are different types.
template < typename T1, typename T2 >
inline T1 SIGN(const T1 &a, const T2 &b)
{
    int s = 1;
    if (b < 0) s = -1;
    return static_cast <T1>(s * ABS(a));
}

template < typename T >
inline T SUMMARY(vector < T > &a, int n)
{
    T sum = 0.0;
    for (int i = 0; i < n; ++ i)
    {
        sum += a[i];
    }
    return sum;
}

//! Coordinate rotate with X axis:
//! x' = x.
//! y' = y*cos(degree) + z*sin(degree).
//! z' = z*cos(degree) - y*sin(degree).
template < typename T >
void RotateAxisX(const T *orginalR, T *newR, const RDouble &degree)
{
    newR[0] = orginalR[0];
    newR[1] = orginalR[1] * cos(degree) + orginalR[2] * sin(degree);
    newR[2] = orginalR[2] * cos(degree) - orginalR[1] * sin(degree);
}

//! Coordinate rotate with Y axis:
//! x' = x*cos(degree) - z*sin(degree).
//! y' = y.
//! z' = z*cos(degree) + x*sin(degree).
template < typename T >
void RotateAxisY(const T *orginalR, T *newR, const RDouble &degree)
{
    newR[0] = orginalR[0] * cos(degree) - orginalR[2] * sin(degree);
    newR[1] = orginalR[1];
    newR[2] = orginalR[2] * cos(degree) + orginalR[0] * sin(degree);
}

//! Coordinate rotate with Z axis:
//! x' = x*cos(degree) + y*sin(degree).
//! y' = y*cos(degree) - x*sin(degree).
//! z' = z.
template < typename T >
void RotateAxisZ(const T *orginalR, T *newR, const RDouble &degree)
{
    newR[0] = orginalR[0] * cos(degree) + orginalR[1] * sin(degree);
    newR[1] = orginalR[1] * cos(degree) - orginalR[0] * sin(degree);
    newR[2] = orginalR[2];
}

template < typename T1, typename T2, typename T3 >
void MatrixMultiply(T1 **a, T2 *b, T3 *c, int m, int n)
{
    for (int i = 0; i < m; ++ i)
    {
        c[i] = 0;
        for (int j = 0; j < n; ++ j)
        {
            c[i] += a[i][j] * b[j];
        }
    }
}

template < typename T1, typename T2, typename T3 >
void MatrixMultiply(T1 **a, vector< T2 > &b, vector< T3 > &c, int m, int n)
{
    MatrixMultiply(a, &b[0], &c[0], m, n);
}

template<typename T> 
string toString(const T &t)
{ 
    ostringstream oss;
    oss << t;
    return oss.str();
}

//! Compute the cross product of two vectors.
//! @param[in]  vec1        The 1st vector, 3D array.
//! @param[in]  vec2        The 2nd vector, 3D array.
//! @param[out] crossVec    The cross product, 3D array.
void CrossProduct(RDouble *vec1, RDouble *vec2, RDouble *crossVec);

template < typename T1, typename T2, typename T3 >
void CrossProduct(vector< T1 > &a, vector< T2 > &b, vector< T3 > &crossProduct)
{
    CrossProduct(&a[0], &b[0], &crossProduct[0]);
}

//! Compute the unit face (a line with two end-points) normal of 2D case.
//! @param[in] startPoint     Start point, 2D array.
//! @param[in] endPoint       End point, 2D array.
//! @param[out] unitNormal    The unit normal direction, 3D array, unitNormal[2] is the length of the original face.
void UnitFaceNormal2D(RDouble *startPoint, RDouble *endPoint, RDouble *unitNormal);

//! Compute the area of the quadrilateral that was composed by vector of r1, r2, r3, r4.
//! The area is called "volume" also in 2D mesh.
//! @param[in] r1       The 1st vector, 3D array.
//! @param[in] r2       The 2nd vector, 3D array.
//! @param[in] r3       The 3rd vector, 3D array.
//! @param[in] r4       The 4th vector, 3D array.
//! @param[out] area    The area of quadrilateral, that is called "volume" also in 2D mesh.
void Area2D(RDouble *r1, RDouble *r2, RDouble *r3, RDouble *r4, RDouble &area);

//! Compute the face normal of 3D case.
//! @param[in] r1    The 1st vector, 3D array.
//! @param[in] r2    The 2nd vector, 3D array.
//! @param[in] r3    The 3rd vector, 3D array.
//! @param[in] r4    The 4th vector, 3D array.
//! @param[out] faceNormal    The normal direction, 3D array.
void FaceNormal3D(RDouble *r1, RDouble *r2, RDouble *r3, RDouble *r4, RDouble *faceNormal);

//! Normalize the input vector for nDim-dimension vector.
//! @param[in]  vecIn     The input vector without normalization, with dimension of 'nDim'.
//! @param[out] vecOut    The output vector that has been normalized, with dimension of 'nDim + 1'.
//! @param[in]  nDim      dimension of the input vector that need to be normalized.
void NormalizeVector(RDouble *vecIn, RDouble *vecOut, int nDim);

//! Compute the six-times volume of a tetrahedron.
//! @param[in] r1,r2,r3,r4    The four edges of the tetrahedron.
//! @param[in] vol    The six-times of the real tetrahedron volume.
void Vol_Tetrahedron(RDouble *r1, RDouble *r2, RDouble *r3, RDouble *r4, RDouble &vol);

//! ComputDist.
//! @param[in] x1,y1,z1,x2,y2,z2
double ComputDist(double x1, double y1, double z1, double x2, double y2, double z2);

//! Determinant.
//! @param[in] A1,A2,A3,B1,B2,B3,C1,C2,C3
double Determinant(double A1, double A2, double A3,
                   double B1, double B2, double B3,
                   double C1, double C2, double C3);

void RestrictIndex(int &i, int &j, int &k, const int &maxi, const int &maxj, const int &maxk);

//! To solve the linear equation like the format A * x = b, where A is an matrix, while both x and b are vectors.
//! @param[in]: diagonalMatrix is the left-hand side matrix.
//! @param[in]: rhsMatrix is the right-hand side vector.
//! @param[in]: nEquation denotes the number of the equation.
//! @param[out]: deltaQ is solution vector of the linear equation.
void SolveLinearEquation(RDouble2D *diagonalMatrix, RDouble *rhsMatrix, int nEquation, RDouble *deltaQ);

}