#include "Math_BasisFunction.h"
#include "TypeDefine.h"

namespace PHSPACE
{

int GetProgressInterval(int nTotal, int nInterval)
{
    nInterval = MIN(nInterval, nTotal);
    int interval = MAX((int)((nTotal * 1.0) / (nInterval * 1.0)), 1);
    return interval;
}

void UnitFaceNormal2D(RDouble *startPoint, RDouble *endPoint, RDouble *unitNormal)
{
    RDouble dx, dy, ds, ods;

    dx = endPoint[0] - startPoint[0];
    dy = endPoint[1] - startPoint[1];
    ds = sqrt(dx * dx + dy * dy);
    const RDouble one = 1.0;
    ods = one / (ds + TINY);

    unitNormal[0] =   dy * ods;
    unitNormal[1] = - dx * ods;
    unitNormal[2] =   ds;
}

void Area2D(RDouble *r1, RDouble *r2, RDouble *r3, RDouble *r4, RDouble &area)
{
    RDouble v1[3], v2[3], p1[3], p2[3];

    for (int i = 0; i < 3; ++ i)
    {
        v1[i] = r2[i] - r1[i];
        v2[i] = r3[i] - r2[i];
    }

    CrossProduct(v1, v2, p1);

    for (int i = 0; i < 3; ++ i)
    {
        v1[i] = r4[i] - r3[i];
        v2[i] = r1[i] - r4[i];
    }

    CrossProduct(v1, v2, p2);

    const RDouble half = 0.5;
    area = half * (p1[2] + p2[2]);
}

void FaceNormal3D(RDouble *r1, RDouble *r2, RDouble *r3, RDouble *r4, RDouble *faceNormal)
{
    //! The following code is used to compute the face normal.
    RDouble v1[3], v2[3], p1[3], p2[3];
    for (int i = 0; i < 3; ++ i)
    {
        v1[i] = r2[i] - r1[i];
        v2[i] = r3[i] - r2[i];
    }

    CrossProduct(v1, v2, p1);

    for (int i = 0; i < 3; ++ i)
    {
        v1[i] = r4[i] - r3[i];
        v2[i] = r1[i] - r4[i];
    }

    CrossProduct(v1, v2, p2);

    const RDouble half = 0.5;
    for (int i = 0; i < 3; ++ i)
    {
        faceNormal[i] = half * (p1[i] + p2[i]) + SMALL;
    }
}

void NormalizeVector(RDouble *vecIn, RDouble *vecOut, int nDim)
{
    RDouble s = 0.0;
    for (int i = 0; i < nDim; ++ i)
    {
        s += vecIn[i] * vecIn[i];
    }
    vecOut[nDim] = sqrt(s);

    const RDouble one = 1.0;
    RDouble os = one / vecOut[nDim];

    for (int i = 0; i < nDim; ++ i)
    {
        vecOut[i] = vecIn[i] * os;
    }
}

void CrossProduct(RDouble *vec1, RDouble *vec2, RDouble *crossVec)
{
    //! The following code is used to compute cross product.
    crossVec[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    crossVec[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    crossVec[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

void Vol_Tetrahedron(RDouble *r1, RDouble *r2, RDouble *r3, RDouble *r4, RDouble &vol)
{
    //! The tetrahedron is six-times as much as the normal volume.
    RDouble a[3], b[3], c[3];
    for (int i = 0; i < 3; ++ i)
    {
        a[i] = r2[i] - r1[i];
        b[i] = r3[i] - r1[i];
        c[i] = r4[i] - r1[i];
    }
    vol = a[0] * b[1] * c[2] + a[1] * b[2] * c[0] + a[2] * b[0] * c[1]
        - c[0] * b[1] * a[2] - c[1] * b[2] * a[0] - c[2] * b[0] * a[1];
}

double ComputDist(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double dx, dy, dz;
    dx = x2 - x1;
    dy = y2 - y1;
    dz = z2 - z1;

    return sqrt(dx * dx + dy * dy + dz * dz);
}

double Determinant(double A1, double A2, double A3,
                   double B1, double B2, double B3,
                   double C1, double C2, double C3)
{
    double det;
    det = A1 * B2 * C3 + B1 * C2 * A3 + A2 * B3 * C1 
        - A3 * B2 * C1 - A2 * B1 * C3 - B3 * C2 * A1;

    return det;
}

void RestrictIndex(int &i, int &j, int &k, const int &maxi, const int &maxj, const int &maxk)
{
    i = MIN(MAX(1, i), MAX(1, maxi));
    j = MIN(MAX(1, j), MAX(1, maxj));
    k = MIN(MAX(1, k), MAX(1, maxk));
}

void SolveLinearEquation(RDouble2D *dMatrix, RDouble *rhsMatrix, int nEquation, RDouble *deltaQ)
{
    //! Allocate new matrix or vector in order to avoid changing the values of the original variables.
    Range II(0, nEquation - 1);
    RDouble2D *aMatrix = new RDouble2D(II, II, fortranArray);
    //! Transformation.
    RDouble2D &A = *reinterpret_cast<RDouble2D *>(aMatrix);
    //! Initialize the matrix of A.
    for (int m = 0; m < nEquation; ++ m)
    {
        for (int n = 0; n < nEquation; ++ n)
        {
            A(m, n) = (*dMatrix)(m, n);
        }
    }

    RDouble *B = new RDouble[nEquation];
    for (int m = 0; m < nEquation; ++ m)
    {
        B[m] = rhsMatrix[m];
    }

    //! To solve the linear equation called A * X = B using the Guass elimination method.
    for (int k = 0; k < nEquation - 1; ++ k)
    {
        //! Select the main element which is the maximum absolution value.
        RDouble maxElement = fabs(A(k, k));
        int p = k;    //! To record the row number of the main element.
        for (int i = k + 1; i < nEquation; ++ i)
        {
            if (fabs(A(i, k)) > maxElement)
            {
                p = i;
                maxElement = fabs(A(i, k));
            }
        }

        if (p != k)    //! Swap the elements between two rows.
        {
            RDouble temp = 0.0;
            for (int j = k; j < nEquation; ++ j)
            {
                temp = A(k, j);
                A(k, j) = A(p, j);
                A(p, j) = temp;
            }
            temp = B[k];
            B[k] = B[p];
            B[p] = temp;
        }

        //! To eliminate the elements.
        for (int i = k + 1; i < nEquation; ++ i)
        {
            RDouble mik = 0.0;
            if (fabs(A(k, k)) >= 1.0E-20)
            {
                mik = A(i, k) / A(k, k);
            }
            else
            {
                mik = A(i, k) / (A(k, k) + 1.0E-20);
            }
            for (int j = k + 1; j < nEquation; ++ j)
            {
                A(i, j) += - mik * A(k, j);
            }
            B[i] += - mik * B[k];
        }
    }

    //! The step of back substitution.
    if (fabs(A(nEquation - 1, nEquation - 1)) >= 1.0E-20)
    {
        deltaQ[nEquation - 1] = B[nEquation - 1] / A(nEquation - 1, nEquation - 1);
    }
    else
    {
        deltaQ[nEquation - 1] = B[nEquation - 1] / (A(nEquation - 1, nEquation - 1) + 1.0E-20);
    }

    for (int i = nEquation - 2; i >= 0; -- i)
    {
        RDouble fSum = 0.0;
        for (int j = i + 1; j < nEquation; ++ j)
        {
            fSum += A(i, j) * deltaQ[j];
        }
        if (fabs(A(i, i)) >= 1.0E-20)
        {
            deltaQ[i] = (B[i] - fSum) / A(i, i);
        }
        else
        {
            deltaQ[i] = (B[i] - fSum) / (A(i, i) + 1.0E-20);
        }
    }

    //! Deallocate the memories.
    delete aMatrix;
    delete [] B;    B = NULL;
}

}