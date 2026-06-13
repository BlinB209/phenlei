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
//! @file      PHMatrix.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include <cstddef>
#include <vector>
#include <cmath>
#include "DataStruct_Array.h"
#include "Pointer.h"
#include "Constants.h"
using namespace std;

namespace PHSPACE
{
template < typename T >
class MAT3D
{
private:
    int ni, nj, nk;
    T ***elems;
    int Nelem;
public:
    MAT3D(int ni, int nj, int nk)
    {
        this->ni = ni;
        this->nj = nj;
        this->nk = nk;
        Nelem = ni * nj * nk;
        elems       = new T ** [ni];
        elems[0]    = new T *  [ni * nj];
        elems[0][0] = new T    [ni * nj * nk];
        for (int i = 1; i < ni; ++ i)
        {
            elems[i] = &elems[i-1][nj];
        }
        
        for (int j = 1; j < nj; ++ j)
        {
            elems[0][j] = &elems[0][j-1][nk];
        }

        for (int i = 1; i < ni; ++ i)
        {
            elems[i][0] = &elems[i-1][nj-1][nk];
            for (int j = 1; j < nj; ++ j)
            {
                elems[i][j] = &elems[i][j-1][nk];
            }
        }
    }

    ~MAT3D()
    {
        delete [] elems[0][0];
        delete [] elems[0];
        delete [] elems;
    }

    T **  GetT2D(int i) const { return elems[i]; }
    T *** GetT3D() const { return elems; }

    void SwitchSign()
    {
        T *mat = elems[0][0];
        for (int i = 0; i < Nelem; ++ i)
        {
            mat[i] = - mat[i];
        }
    }
};

template < typename T >
inline void SetValue(T *v, int ni, const T &x)
{
    for (int i = 0; i < ni; ++ i)
    {
        v[i] = x;
    }
}

template < typename T >
inline void SetValue(T **v, int ni, int nj, const T &x)
{
    for (int i = 0; i < ni; ++ i)
    {
        for (int j = 0; j < nj; ++ j)
        {
            v[i][j] = x;
        }
    }
}

template < typename T >
inline void SetValue(T **x, T **y, int ni, int nj, int coef)
{
    for (int i = 0; i < ni; ++ i)
    {
        for (int j = 0; j < nj; ++ j)
        {
            x[i][j] = coef * y[i][j];
        }
    }
}

template < typename T >
inline void SetValue(T **x, T **y, int ni, int nj)
{
    for (int i = 0; i < ni; ++ i)
    {
        for (int j = 0; j < nj; ++ j)
        {
            x[i][j] = y[i][j];
        }
    }
}


template < typename T >
void AllocateVector(vector < T * > &a, int length)
{
    for (int m = 0; m < a.size(); ++ m)
    {
        a[m] = new T [ length ];
    }
}

template < typename T >
void AllocateVector(vector < T * > *&a, int length, int nm)
{
    a = new vector < T * >(nm);
    for (int m = 0; m < a.size(); ++ m)
    {
        (*a)[m] = new T [length];
    }
}

template < typename T >
void AllocateVector(vector< vector< T > > &data, int numberOfEquations, int numberOfElements)
{
    data.resize(numberOfEquations);
    for (int iEquation = 0; iEquation < numberOfEquations; ++ iEquation)
    {
        data[iEquation].resize(numberOfElements);
    }
}

template < typename T >
void DeAllocateVector(vector < T * > &a)
{
    for (int m = 0; m < a.size(); ++ m)
    {
        delete [] a[m];
    }
}

template < typename T >
void DeAllocateVector(vector < T * > *a)
{
    for (int m = 0; m < a->size(); ++ m)
    {
        delete [] (*a)[m];
    }
    delete a;
}

// Compute the dot product of two vectors with length n
template < typename T >
T DotProduct(T *a, T *b, int n)
{
    T sum = 0.0;
     for (int i = 0; i < n; ++ i)
    {
        sum += a[i] * b[i];
    }
    return sum;
}

// Compute the dot product of two vectors with length n
template < typename T >
T GetNorm(T *a, int n)
{
    T sum = 0.0;
    for (int i = 0; i < n; ++ i)
    {
        sum += a[i] * a[i];
    }
    return sqrt(sum);
}

template < typename T >
T GetAngleOfTwoVector(T* x, T* y, int n)
{
    RDouble cosnum = DotProduct(x, y, n)/GetNorm(x, n)/GetNorm(y, n);
    if (cosnum > 1.0)
    {
        return 0;
    }
    else if (cosnum < -1.0)
    {
        return PI;
    }
    else
    {
        return acos(cosnum);
    }
}

template < typename T >
T ** MATMUL(T **matrixA, T **matrixB)
{
    int nrow = 3;
    int ncol = 3;

    T **result = NewPointer2 <T> (3, 3);

    for (int row = 0; row < 3; ++ row)
    {
        for (int col = 0; col < 3; ++ col)
        {
            T sum = static_cast< T >(0);
            for (int i = 0; i < 3; ++ i)
            {
                sum += matrixA[row][i] * matrixB[i][col];
            }
            result[row][col] = sum;
        }
    }
    return result;
}

template < typename T >
T * MATMUL(T **matrixA, T *matrixB)
{
    int nrow = 3;
    int ncol = 3;

    T *result = new T [3];

    for (int row = 0; row < 3; ++ row)
    {
        T sum = static_cast< T >(0);
        for (int i = 0; i <3; ++ i)
        {
            sum += matrixA[row][i] * matrixB[i];
        }
        result[row] = sum;
    }
    return result;
}

//solve the dq by LU decompound methods.
template < typename T >
void ScalarLUSolve(T **A, T *B, int N)
{
    // Forward substitution
    B[0] /= A[0][0];
    for (int i = 1; i < N; ++ i)
    {
        for (int k = 0; k < i; ++ k)
        {
            B[i] -= A[i][k] * B[k];
        }
        B[i] /= A[i][i];
    } 

    // Backward substitution
    for (int i = N - 1; i >= 0; -- i)
    {
        for (int k = i + 1; k < N; ++ k)
        {
            B[i] -= A[i][k] * B[k];
        }
    }
}

template < typename T >
void LUDCMPBlockDiag(MAT3D< T > *impmat, int len, int N)
{
    for (int iCell = 0; iCell < len; ++ iCell)
    {
        T **A = impmat->GetT2D(iCell);
        for (int j = 0; j < N; ++ j)
        {
            // Lower triangular matrix
            for (int i = j; i < N; ++ i)
            {
                for (int k = 0; k < j; ++ k)
                {
                    A[i][j] -= A[i][k] * A[k][j];
                }
            }

            // Upper triangular matrix
            for (int i = j + 1; i < N; ++ i)
            {
                for (int k = 0; k < j; ++ k)
                {
                    A[j][i] -= A[j][k] * A[k][i];
                }
                A[j][i] /= A[j][j];
            }
        }
    }
}

template < typename T >
void BlockLUInverse(T ***sub_matrix_A, T ***sub_matrix_B, T ***sub_matrix_C, int len, int N)
{
    for (int iCell = 0; iCell < len - 1; ++ iCell)
    {
        //A^-1 B
        //forward elimination
        for (int k = 0; k < N - 1; ++ k)
        {
            //scale row k by A[k][k]
            for (int j = k + 1; j < N; ++ j)
            {
                sub_matrix_A[iCell][k][j] /= sub_matrix_A[iCell][k][k];
            }
            for (int j = 0; j < N; ++ j)
            {
                sub_matrix_B[iCell][k][j] /= sub_matrix_A[iCell][k][k];
            }

            //loop over rows
            for (int i = k + 1; i < N; ++ i)
            {
                for (int j = k + 1; j < N; ++ j)
                {
                    sub_matrix_A[iCell][i][j] -= sub_matrix_A[iCell][i][k] * sub_matrix_A[iCell][k][j];
                }
                for (int j = 0; j < N; ++ j)
                {
                    sub_matrix_B[iCell][i][j] -= sub_matrix_A[iCell][i][k] * sub_matrix_B[iCell][k][j];
                }
            }
        } //end of k

        //scale by A[4][4]
        for (int j = 0; j < N; ++ j)
        {
            sub_matrix_B[iCell][N - 1][j] /= sub_matrix_A[iCell][N - 1][N - 1];
        }

        //end forward elimination && start backward substitution
        //loop over rows
        for (int k = N - 2; k >= 0; -- k)
        { ///3
            for (int i = k + 1; i < N; ++ i)
            {
                for (int j = 0; j < N; ++ j)
                {
                    sub_matrix_B[iCell][k][j] -= sub_matrix_A[iCell][k][i] * sub_matrix_B[iCell][i][j];
                }
            }
        } //end backward substitution

        // compute A[cell+1] := A[cell+1] - C[cell+1] * B[cell];
        for (int i = 0; i < N; ++ i)
        {
            for (int j = 0; j < N; ++ j)
            {
                for (int k = 0; k < N; ++ k)
                {
                    sub_matrix_A[iCell + 1][i][j] -= sub_matrix_C[iCell + 1][i][k] * sub_matrix_B[iCell][k][j];
                }
            }
        }
    } //end of jCell cell

    //forward elimination
    for (int k = 0; k < N - 1; ++ k)
    {
        //scale row k by A[k][k]
        for (int j = k + 1; j < N; ++ j)
        {
            sub_matrix_A[len - 1][k][j] /= sub_matrix_A[len - 1][k][k];
        }
        //loop over rows
        for (int i = k + 1; i < N; ++ i)
        {
            for (int j = k + 1; j < N; ++ j)
            {
                sub_matrix_A[len - 1][i][j] -= sub_matrix_A[len - 1][i][k] * sub_matrix_A[len - 1][k][j];
            }
        }
    } //end of k end forward elimination
}

//! Solve tri-diagonal block system
template < typename T >
void BlockLUSolve(T **DQ, T ***sub_matrix_A, T ***sub_matrix_B, T ***sub_matrix_C, int len, int N)
{
    //forward elimination for the tridiagonal block system
    for (int iCell = 0; iCell < len - 1; ++ iCell)
    {
        //A^-1 dq
        //forward elimination
        for (int k = 0; k < N - 1; ++ k)
        {
            //scale row k by A[k][k]
            DQ[iCell][k] /= sub_matrix_A[iCell][k][k];
            //loop over rows
            for (int i = k + 1; i < N; ++ i)
            {
                DQ[iCell][i] -= sub_matrix_A[iCell][i][k] * DQ[iCell][k];
            }
        } //end of k

        //scale by A[4][4]
        DQ[iCell][N - 1] /= sub_matrix_A[iCell][N - 1][N - 1];

        //end forward elimination && start backward substitution
        //loop over rows
        for (int k = N - 2; k >= 0; -- k)
        { ///3
            for (int i = k + 1; i < N; ++ i)
            {
                DQ[iCell][k] -= sub_matrix_A[iCell][k][i] * DQ[iCell][i];
            }
        } //end backward substitution

        // compute dq[cell+1] := dq[cell+1] - C[cell+1] * dq[cell];
        for (int i = 0; i < N; ++ i)
        {
            for (int k = 0; k < N; ++ k)
            {
                DQ[iCell + 1][i] -= sub_matrix_C[iCell + 1][i][k] * DQ[iCell][k];
            }
        }
    } //end of jCell cell

    //A^-1 dq at len-1
    //forward elimination
    for (int k = 0; k < N - 1; ++ k)
    {
        //scale row k by A[k][k]
        DQ[len - 1][k] /= sub_matrix_A[len - 1][k][k];
        //loop over rows
        for (int i = k + 1; i < N; ++ i)
        {
            DQ[len - 1][i] -= sub_matrix_A[len - 1][i][k] * DQ[len - 1][k];
        }
    } //end of k end forward elimination

    //scale by A[4][4]
    DQ[len - 1][N - 1] /= sub_matrix_A[len - 1][N - 1][N - 1];

    //end forward elimination && start backward substitution
    //loop over rows
    for (int k = N - 2; k >= 0; -- k)
    { ///3
        for (int i = k + 1; i < N; ++ i)
        {
            DQ[len - 1][k] -= sub_matrix_A[len - 1][k][i] * DQ[len - 1][i];
        }
    } //end backward substitution

    //backward substitution for the block tri-diagonal system
    //loop over block rows
    for (int iCell = len - 2; iCell >= 0; -- iCell)
    { ///-2
        //compute dq = dq- B dq[cell+1]
        for (int i = 0; i < N; ++ i)
        {
            for (int k = 0; k < N; ++ k)
            {
                DQ[iCell][i] -= sub_matrix_B[iCell][i][k] * DQ[iCell + 1][k];
            }
        }
    }
}

}
