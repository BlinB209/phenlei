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
//! @file      PHSuperMatrix.h
//! @brief     Explain this file briefly.
//! @author    xxx.

#pragma once
#include <vector>
#include "PHHeader.h"
#include "TK_Exit.h"
#include "Math_BasisFunction.h"

using namespace std;

namespace PHSPACE
{

template < typename T >
class PHMatrix
{
protected:
    int _nrow, _ncolumn, _nsize;
    PHVector1D< T > _data;
public:
    PHMatrix(int nrow, int ncolumn = 1)
    {
        _nrow    = nrow;
        _ncolumn = ncolumn;
        _nsize   = _nrow * _ncolumn;
        _data    = PHVector1D< T >(_nsize);
    }

    ~PHMatrix()
    {
    }
    int rows   () const { return _nrow; }
    int columns() const { return _ncolumn; }
    int sizes  () const { return _nsize; };
    PHVector1D< T > datas() { return _data; }

    const T & operator() (int row, int column = 1) const 
    {
        if ((row < 0 || row > _nrow) || (column < 0 || column > _ncolumn))
        {
            cout << endl << "Error : overflow in type PHMatrix ! ";
            TK_Exit::ExceptionExit("");
        }
        return  _data[(row - 1) * _ncolumn + (column - 1)] ;
    }

    T & operator()(int row, int column = 1)
    {
        if ((row < 0 || row > _nrow) || (column < 0 || column > _ncolumn))
        {
            TK_Exit::ExceptionExit("Error : overflow in type PHMatrix ! ");
        }
        return  _data[(row - 1) * _ncolumn + (column - 1)] ;
    }

};

template < typename T >
PHMatrix< T > operator * (const PHMatrix<T> & matrixA, const PHMatrix<T> & matrixB)
{
    int nrow = matrixA.rows();
    int ncol = matrixB.columns();

    if (matrixA.columns() != matrixB.rows())
    {
        TK_Exit::ExceptionExit("Error : Index not matched in PHMatrix<T> operator * ! ");
    }

    PHMatrix< T > result(nrow, ncol);

    for (int col = 1; col <= ncol; ++ col)
    {
        for (int row = 1; row <= nrow; ++ row)
        {
            T sum = static_cast< T >(0);
            for (int i = 1; i <= matrixA.columns(); ++ i)
            {
                sum += matrixA(row, i) * matrixB(i, col);
            }
            result(row, col) = sum;
        }
    }
    return result;
}

template < typename T >
PHMatrix< T > MATMUL(PHMatrix<T> & matrixA, PHMatrix<T> & matrixB)
{
    int nrow = matrixA.rows();
    int ncol = matrixB.columns();

    if (matrixA.columns() != matrixB.rows())
    {
        TK_Exit::ExceptionExit("Error : Index not matched in function MATMUL ! ");
    }

    PHMatrix<T> result(nrow, ncol);

    for (int row = 1; row <= nrow; ++ row)
    {
        for (int col = 1; col <= ncol; ++ col)
        {
            T sum = static_cast< T >(0);
            for (int i = 1; i <= matrixA.columns(); ++ i)
            {
                sum += matrixA(row, i) * matrixB(i, col);
            }
            result(row, col) = sum;
        }
    }
    return result;
}

template < typename T >
PHMatrix<T> TransposeMatrix(PHMatrix<T> & matrixA)
{
    int nrow = matrixA.rows();
    int ncol = matrixA.columns();

    if (nrow != ncol)
    {
        TK_Exit::ExceptionExit("Error : Index not matched in function TransposeMatrix ! ");
    }

    PHMatrix<T> matrixB(nrow, ncol);

    for (int row = 1; row <= nrow; ++ row)
    {
        for (int col = 1; col <= ncol; ++ col)
        {
            matrixB(row, col) = matrixA(col, row);
        }
    }
    return matrixB;
}

template < typename T >
PHMatrix<T>  Matrix3X3Inverse(PHMatrix<T> & matrix)
{
    PHMatrix<T> inverse(3, 3);

    if (matrix.rows() != 3 || matrix.columns() != 3)
    {
        TK_Exit::ExceptionExit("Error : Index must be 3 in function Matrix3X3Inverse ! ");
    }

    T detA = matrix(1, 1) * matrix(2, 2) * matrix(3, 3)
           + matrix(1, 2) * matrix(2, 3) * matrix(3, 1)
           + matrix(2, 1) * matrix(3, 2) * matrix(1, 3)
           - matrix(1, 1) * matrix(2, 3) * matrix(3, 2)
           - matrix(2, 2) * matrix(1, 3) * matrix(3, 1)
           - matrix(3, 3) * matrix(2, 1) * matrix(1, 2);

    if (detA == 0.0)
    {
        TK_Exit::ExceptionExit("Error : The Determinant of the Matrix is Zero£¬Singular Matrix ! ");
    }
    else
    {
        inverse(1, 1) = matrix(2, 2) * matrix(3, 3) - matrix(2, 3) * matrix(3, 2);
        inverse(1, 2) = matrix(1, 3) * matrix(3, 2) - matrix(1, 2) * matrix(3, 3);
        inverse(1, 3) = matrix(1, 2) * matrix(2, 3) - matrix(2, 2) * matrix(1, 3);

        inverse(2, 1) = matrix(2, 3) * matrix(3, 1) - matrix(2, 1) * matrix(3, 3);
        inverse(2, 2) = matrix(1, 1) * matrix(3, 3) - matrix(1, 3) * matrix(3, 1);
        inverse(2, 3) = matrix(2, 1) * matrix(1, 3) - matrix(1, 1) * matrix(2, 3);

        inverse(3, 1) = matrix(2, 1) * matrix(3, 2) - matrix(3, 1) * matrix(2, 2);
        inverse(3, 2) = matrix(1, 2) * matrix(3, 1) - matrix(1, 1) * matrix(3, 2);
        inverse(3, 3) = matrix(1, 1) * matrix(2, 2) - matrix(1, 2) * matrix(2, 1);

        for (int i = 1; i <= 3; ++ i)
        {
            for (int j = 1; j <= 3; ++ j)
            {
                inverse(i, j) /=  detA;
            }
        }
    }

    return inverse;

}

template < typename T >
PHMatrix<T> BaseMatrixOfRotate(const T & angle, string const & axisLabel)
{
    PHMatrix<T> matrix(3, 3);

    T cosine = cos(angle);
    T sine   = sin(angle);
    T oneTmp    = static_cast< T >(1);
    T zeroTmp   = static_cast< T >(0);

    if (axisLabel == "x" || axisLabel == "X")
    {
        matrix(1, 1) = oneTmp;
        matrix(1, 2) = zeroTmp;
        matrix(1, 3) = zeroTmp;

        matrix(2, 1) = zeroTmp;
        matrix(2, 2) = cosine;
        matrix(2, 3) = sine;

        matrix(3, 1) = zeroTmp;
        matrix(3, 2) =-sine;
        matrix(3, 3) = cosine;
    } 
    else if (axisLabel == "y" || axisLabel == "Y")
    {
        matrix(1, 1) = cosine;
        matrix(1, 2) = zeroTmp;
        matrix(1, 3) =-sine;

        matrix(2, 1) = zeroTmp;
        matrix(2, 2) = oneTmp;
        matrix(2, 3) = zeroTmp;

        matrix(3, 1) = sine;
        matrix(3, 2) = zeroTmp;
        matrix(3, 3) = cosine;
    }
    else if (axisLabel == "z" || axisLabel == "Z")
    {
        matrix(1, 1) = cosine;
        matrix(1, 2) = sine;
        matrix(1, 3) = zeroTmp;

        matrix(2, 1) =-sine;
        matrix(2, 2) = cosine;
        matrix(2, 3) = zeroTmp;

        matrix(3, 1) = zeroTmp;
        matrix(3, 2) = zeroTmp;
        matrix(3, 3) = oneTmp;
    }
    else
    {
        TK_Exit::ExceptionExit("Error : error axis label in BaseMatrixOfRotate ! ");
    }

    return matrix;
}

template < typename T >
PHMatrix<T> AngleRateTransformMatrixBodyToGlobal(PHMatrix<T> & eulerangle, string const & rotateorder)
{
    PHMatrix<T> matrix(3, 3);

    T pitchangle  = eulerangle(1);
    T rollangle   = eulerangle(3);

    T oneTmp    = static_cast< T >(1);
    T zeroTmp   = static_cast< T >(0);

    if (rotateorder == "y-z-x" || rotateorder == "Y-Z-X" || rotateorder == "2-3-1")
    {
        T cosineofpitchangle = cos(pitchangle);
        if (PHSPACE::ABS(cosineofpitchangle) < SMALL)
        {
            TK_Exit::ExceptionExit("Error : pitching angle is 90 DEG£¬singular ! ");
        }
        else
        {
            T sineofpitchangle  = sin(pitchangle);
            T cosineofrollangle = cos(rollangle);
            T sineofrollangle   = sin(rollangle);

            matrix(1, 1) = zeroTmp;
            matrix(1, 2) = sineofrollangle;
            matrix(1, 3) = cosineofrollangle;

            matrix(2, 1) = zeroTmp;
            matrix(2, 2) = cosineofrollangle / cosineofpitchangle;
            matrix(2, 3) =-sineofrollangle / cosineofpitchangle;

            matrix(3, 1) = oneTmp;
            matrix(3, 2) =-cosineofrollangle * sineofpitchangle / cosineofpitchangle;
            matrix(3, 3) = sineofrollangle * sineofpitchangle / cosineofpitchangle;
        }
    } 
    else
    {
        TK_Exit::ExceptionExit("Error : Undefined Euler Angle transformation order ! ");
    }

    return matrix;
}

template < typename T >
PHMatrix<T> CrossProductVector(PHMatrix<T> a, PHMatrix<T> b)
{
    PHMatrix<T> cp(a.rows());
    cp(1) = a(2) * b(3) - a(3) * b(2);
    cp(2) = a(3) * b(1) - a(1) * b(3);
    cp(3) = a(1) * b(2) - a(2) * b(1);
    return cp;
}

}