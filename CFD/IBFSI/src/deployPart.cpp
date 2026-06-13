#include "IB.hpp"
#include <sstream>

//! 构造函数  实例对象时进行初始设置
IBFSI::IBFSI(int nx, int ny, int nz, double dT, vector3<double> vel0_, double density0, double elesizex, double elesizey, double elesizez, string input_IB, string meshFlieName, string outFileName)
{
    initialise(nx, ny, nz, elesizex, elesizey, elesizez, vel0_, outFileName);
    getInputFromFile(input_IB);    //! 读取IB参数 
    getLagFromMesh(meshFlieName);    //! 从Mesh中读取Lag
    setdenIn(density0);
    setCharacterQuantity();    //! 设置特征量
    initialBasicParams();    //! 初始化基础参数
    init_VIV();    //! 初始化涡激振动

    cout << "ds0 = " << ds0 << " ,  " << "dt = " << dt << endl;
    cout << "Lref = " << Lref << " ,  " << "Uref = " << Uref << endl;
}

/**
* velocity_Conversion_vector3_To_vector    将vector3的变量转化为vector类型用于数组计算
* 
* @param vel    从LBM获取的vector3类型的速度
*
*/
void IBFSI::velocity_Conversion_vector3_To_vector(vector3<double>* uvw_, double* rho_)
{
    if (Nz != 1)
    {
        cout << "error Nz != 1: At present, IB algorithm only supports two-dimensional computation!" << endl;
        exit(1);
    }

    int k = (Nz == 1) ? 0 : 0;
    //! 赋值到IBFSI的u、v、w速度中
    for (int j = 0; j < Ny; ++j)
    {
        for (int i = 0; i < Nx; ++i)
        {
            int index = i + j * Nx + k * Nx * Ny;
            u[j][i] = uvw_[index].x;
            v[j][i] = uvw_[index].y;
            w[j][i] = uvw_[index].z;
            rho[j][i] = rho_[index];
        }
    }
}


/**
* Conversion_vector_To_vecto3    将计算得到的vector类型体积力转换为vector3类型，用于输出
*
*/
 void IBFSI::Conversion_vector_To_vector3(vector3<double> **ff)
{
    int k = (Nz == 1) ? 0 : 0;

    //读取fxfyfz体积力输出
    //for (int k = 0; k < Nz; k++) {
    for (int j = 0; j < Ny; ++j)
    {
        for (int i = 0; i < Nx; ++i)
        {
            int ind = i + j * Nx + k * Nx * Ny;
            (*ff)[ind].x = fx[k][j][i];
            (*ff)[ind].y = fy[k][j][i];
            (*ff)[ind].z = fz[k][j][i];
        }
    }
}

/**
* getNumberFromString    提取字符串中整型数据
*
* @param a    输出的int数组，通过int*创建
*/
void IBFSI::IBgetNumberFromString(string s, int *a)
{
    stringstream str_strm;
    str_strm << s;
    string temp_str;
    int temp_int;
    int k = 0;
    while (!str_strm.eof())
    {
        str_strm >> temp_str;
        if (stringstream(temp_str) >> temp_int)
        {   //! try to convert string to int
            a[k] = temp_int;
            //cout << temp_int << endl;
            k++;
        }
        temp_str = "";
    }
}
