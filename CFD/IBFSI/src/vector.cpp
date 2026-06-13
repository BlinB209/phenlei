#include "vector.hpp"
#include <iostream>

const double epsilon = 1e-8;

vector<int> creatmatrix_1D_int(int num, int value)
{
    vector<int> v;
    for (int i = 0; i < num; i++)
    {
        v.push_back(value);
    }
    return v;
}


/**
* creatmatrix_double    创建H列的一维矩阵(double)
*
* @return    vector<double>类型的数组
*/
vector<double> creatmatrix_1D(int num, double value)
{
    vector<double> v;
    for (int i = 0; i < num; i++)
    {
        v.push_back(value);
    }
    return v;
}


/**
* creatmatrix_int    创建H行L列的二维矩阵(int)
*
* @param H,L    行数，列数
* @param value    创建的二维数组的值；默认形参，默认值为0
*
* @return    vector<vector<double>>类型的矩阵(二维)
*/
vector<vector<int>> creatmatrix_int(int H, int L, int value)
{
    vector<vector<int>> v;
    for (int i = 0; i < H; i++)
    {
        vector<int>v1(L, value);
        v.push_back(v1);
    }
    return v;
}

/**
* creatmatrix_double    创建H行L列的二维矩阵(double)
* 
* @param H,L    行数，列数
* @param value    创建的二维数组的值；默认形参，默认值为0
* 
* @return    vector<vector<double>>类型的矩阵(二维)
*/
vector<vector<double>> creatmatrix_double(int H, int L, double value)
{
    vector<vector<double>> v;
    for (int i = 0; i < H; i++)
    {
        vector<double>v1(L, value);
        v.push_back(v1);
    }
    return v;
}

/**
* creatDiagMatrix_double    将一维矩阵转换为二维对角矩阵(double)
* 
* @return 所需对角阵
*/
vector<vector<double>> creatDiagMatrix_double(const vector<double> &A)
{
    int h_num = static_cast<int>(A.size());

    vector<vector<double>> A_D = creatmatrix_double(h_num, h_num);

    for (int i = 0; i < h_num; i++)
    {
        for (int j = 0; j < h_num; j++)
        {
            if (i == j)
            {
                A_D[i][j] = A[i];
            }
            else
            {
                A_D[i][j] = 0.;
            }
        }
    }
    return A_D;
}

/**
* add1D_num    一个数与一维矩阵相加(double + double)
*
* @param A    加法使用矩阵
* @param num    加法使用的数
*
* @return    vector<double>类型的矩阵(一维)
*/
vector<double> add1D_num(const vector<double> &A, double num)
{
    int A_size = static_cast<int>(A.size());
    vector<double> B(A_size, 0.);
    for (int i = 0; i < A_size; i++)
    {
        B[i] = num + A[i];
    }
    return B;
}


/**
* multiply1D_num    一个数乘以一维矩阵(double*double)
* 
* @param A    乘法使用矩阵
* @param num    乘法使用的乘数
* 
* @return    vector<double>类型的矩阵(一维)
*/
vector<double> multiply1D_num(const vector<double> &A, double num)
{
    int A_size = static_cast<int>(A.size());
    vector<double> B(A_size, 0.);
    for (int i = 0; i < A_size; i++)
    {
        B[i] = num * A[i];
    }
    return B;
}

/**
* multiply2D_num    一个数乘以一维矩阵(double*double)
* 
* @param A    乘法使用矩阵
* @param num    乘法使用的乘数
* 
* @return    vector<double>类型的矩阵(一维)
*/
vector<vector<double>> multiply2D_num(const vector<vector<double>> &A, double num)
{
    int A_H = static_cast<int>(A.size());
    int A_L = static_cast<int>(A[0].size());
    vector<vector<double>> B = creatmatrix_double(A_H, A_L);
    for (int i = 0; i < A_H; i++)
    {
        for (int j = 0; j < A_L; j++)
        {
            B[i][j] = num * A[i][j];
        }
    }
    return B;
}

/**
* add1D_matrix    两个一维矩阵相加(double + double),矩阵维数需相同
*
* @param A    加法使用矩阵
* @param B    加法使用另一个矩阵
*
* @return    vector<double>类型的矩阵(一维)
*/
vector<double> add1D_matrix(const vector<double> &A, const vector<double> &B)
{
    int A_size = static_cast<int>(A.size());
    int B_size = static_cast<int>(B.size());

    if (A_size != B_size)
    {
        throw "error add1D_matrix: The matrices that need to be summed have different dimensions!!";
    }

    vector<double> C(A_size, 0.);
    for (int i = 0; i < A_size; i++)
    {
        C[i] = A[i] + B[i];
    }
    return C;
}


/**
* minus1D_matrix    两个一维矩阵相减(double - double),矩阵维数需相同
*
* @param A    减法使用矩阵
* @param B    减法使用另一个矩阵
*
* @return    vector<double>类型的矩阵(一维)
*/
vector<double> minus1D_matrix(const vector<double> &A, const vector<double> &B)
{
    int A_size = static_cast<int>(A.size());
    int B_size = static_cast<int>(B.size());

    if (A_size != B_size)
    {
        throw "error add1D_matrix: The matrices that need to be summed have different dimensions!!";
    }

    vector<double> C(A_size, 0.);
    for (int i = 0; i < A_size; i++)
    {
        C[i] = A[i] - B[i];
    }
    return C;
}

/**
* trans_int    矩阵转置(int)
* 
* @param A    需转置的矩阵
* 
* @return    转置结果(int类型)
*/
vector<vector<int>> trans_int(const vector<vector<int>> &A)
{
    int h = static_cast<int>(A[0].size());
    int l = static_cast<int>(A.size());
    vector<vector<int>> AT = creatmatrix_int(h, l);
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < l; j++)
        {
            AT[i][j] = A[j][i];
        }
    }
    return AT;
}

/**
* trans_double    矩阵转置(double)
*
* @param A    需转置的矩阵
*
* @return    转置结果(double类型)
*/
vector<vector<double>> trans_double(const vector<vector<double>> &A)
{
    int h = static_cast<int>(A[0].size());
    int l = static_cast<int>(A.size());
    vector<vector<double>> AT = creatmatrix_double(h,l);
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < l; j++)
        {
            AT[i][j] = A[j][i];
        }
    }
    return AT;
}

/**
* multiply_double    二维矩阵A*一维矩阵B=一维矩阵C，并返回C
*/
vector<double> multiply_two_one(const vector<vector<double>> &A, const vector<double> &B)
{
    int A_h = static_cast<int>(A.size());
    int A_l = static_cast<int>(A[0].size());
    int B_h = static_cast<int>(B.size());

    if (A_l != B_h)
    {
        throw "error multiply_double: Two matrix dimensions cannot be multiplied!!";
    }

    vector<double> C(A_h, 0.0);
    for (int i = 0; i < A_h; ++i)
    {
        C[i] = 0;
        for (int k = 0; k < A_l; ++k)
        {
            C[i] += A[i][k] * B[k];
        }
        if (fabs(C[i]) < epsilon)
        {
            C[i] = 0.0;
        }
    }
    return C;
}

/**
* multiply_double    二维矩阵A*二维矩阵B=矩阵C，并返回C
*/
vector<vector<double>> multiply_double(const vector<vector<double>> &A, const vector<vector<double>> &B)
{
    int A_h = static_cast<int>(A.size());
    int A_l = static_cast<int>(A[0].size());
    int B_h = static_cast<int>(B.size());
    int B_l = static_cast<int>(B[0].size());

    if (A_l != B_h)
    {
        throw "error multiply_double: Two matrix dimensions cannot be multiplied!!";
    }

    vector<vector<double>> C = creatmatrix_double(A_h, B_l);
    for (int i = 0; i < A_h; i++)
    {
        for (int j = 0; j < B_l; j++)
        {
            C[i][j] = 0;
            for (int k = 0; k < A_l; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
            if (fabs(C[i][j]) < epsilon)
            {
                C[i][j] = 0.0;
            }
        }
    }
    return C;
}

/**
* show_matrix_int    打印显示矩阵（二维）(int)
*/
void show_matrix_int(const vector<vector<int>> &A)
{
    int h = static_cast<int>(A.size());
    int l = static_cast<int>(A[0].size());
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < l; j++)
        {
            cout << A[i][j] << "\t";
        }
        cout << endl;
    }
}

/**
* show_matrix_double    打印显示矩阵（二维）(double)
*/
void show_matrix_double(const vector<vector<double>> &A)
{
    int h = static_cast<int>(A.size());
    int l = static_cast<int>(A[0].size());
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < l; j++)
        {
            cout << A[i][j] << "\t";
        }
        cout << endl;
    }
}

/**
* tile_copy_double    将一个行向量向列方向复制N次
* 
* @param A    double类型的行向量
* @param N    复制次数，返回矩阵的行数
* 
* @return 复制结果,[N,***]
*/
vector<vector<double>> tile_copy_double(vector<double> A, int N)
{
    vector<vector<double>> A_k;
    for (int i = 0; i < N; i++)
    {
        A_k.push_back(A);
    }
    return A_k;
}

/**
* sum_Rows    对二维数组按行求和(double)
* 
* @return 大小为输入数组行数的一维数组
*/
vector<double> sum_Rows(vector<vector<double>> A)
{
    vector<double> sum_Results;

    int h = static_cast<int>(A.size());
    int l = static_cast<int>(A[0].size());
    for (int i = 0; i < h; i++)
    {
        double temp = 0.;
        for (int j = 0; j < l; j++)
        {
            temp += A[i][j];
        }
        sum_Results.push_back(temp);
    }
    return sum_Results;
}