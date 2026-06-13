#include "LBMSolverMPI.hpp"
#include "vector3.cpp"
#include <iomanip>
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include<sstream>
#include<algorithm>
using namespace std;

void LBM::output_all(std::string filename)
{
    std::ofstream output_stream;
    output_stream.open(filename, std::ofstream::out);
    //! 若打开时有option：    | std::ofstream::app，则 app追加文件，不加为默认覆盖
    output_stream.precision(8);
    output_stream << "variables = x, y, z, rho, u, v, w" << '\n';
    output_stream << "zone i=" << total_NX << ",j=" << total_NY << ", k=" << total_NZ << ", f=point" << '\n';

    for (int z = 0; z < total_NZ; z++)
    {
        for (int y = 0; y < total_NY; y++)
        {
            for (int x = 0; x < total_NX; x++)
            {
                output_stream << x << " " << y << " " << z << " " <<
                rhoTotal[this->total_index(x, y, z)] << " " <<
                velTotal[this->total_index(x, y, z)].x << " " <<
                velTotal[this->total_index(x, y, z)].y << " " <<
                velTotal[this->total_index(x, y, z)].z << " " << '\n';
            }
        }
    }
    output_stream.close();
}

void LBM::output_all(std::string filename, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
{
    std::ofstream output_stream;
    output_stream.open(filename, std::ofstream::out);
    //! 若打开时有option：    | std::ofstream::app，则 app追加文件，不加为默认覆盖
    output_stream.precision(8);
    output_stream << "variables = x, y, z, rho, u, v, w" << '\n';
    output_stream << "zone i=" << total_NX + xmax - xmin << ",j=" << total_NY + ymax - ymin << ", k=" << total_NZ + zmax - zmin << ", f=point" << '\n';

    for (int z = zmin; z < total_NZ + zmax; z++)
    {
        for (int y = ymin; y < total_NY + ymax; y++)
        {
            for (int x = xmin; x < total_NX + xmax; x++)
            {
                output_stream << x << " " << y << " " << z << " " <<
                rhoTotal[this->total_index(x, y, z)] << " " <<
                velTotal[this->total_index(x, y, z)].x << " " <<
                velTotal[this->total_index(x, y, z)].y << " " <<
                velTotal[this->total_index(x, y, z)].z << " " << '\n';
            }
        }
    }
    output_stream.close();
}

void LBM::read_all(std::string filename)
{
    std::ifstream input_stream;
    int ab[3];
    char buffer[256];
    double temprho;
    double tempuvw[3];
    int x0, y0, z0;
    input_stream.open(filename, ios::in);
    input_stream.getline(buffer, 200);
    input_stream.getline(buffer, 200);

    for (int z = 0; z < total_NZ; z++)
    {
        for (int y = 0; y < total_NY; y++)
        {
            for (int x = 0; x < total_NX; x++)
            {
                input_stream.getline(buffer, 100);
                sscanf_s(buffer, "%d %d %d %lf %lf %lf %lf", &ab[0], &ab[1], &ab[2], &temprho, &tempuvw[0], &tempuvw[1], &tempuvw[2]);
                x0 = x - Xori;
                y0 = y - Yori;
                z0 = z - Zori;
                if (x0 >= 0 && x0 < NX && y0 >= 0 && y0 < NY && z0 >= 0 && z0 < NZ)
                {
                    rho[index(x0, y0, z0)] = (double)temprho;
                    vel[index(x0, y0, z0)].x = (double)tempuvw[0];
                    vel[index(x0, y0, z0)].y = (double)tempuvw[1];
                    vel[index(x0, y0, z0)].z = (double)tempuvw[2];
                }
            }
        }
    }
    input_stream.close();

    for (int x = 0; x < NX; x++)
    {
        for (int y = 0; y < NY; y++)
        {
            for (int z = 0; z < NZ; z++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                    pdf[index(x, y, z, i)] = calculate_feq(x, y, z, i);
                }
            }
        }
    }
    input_stream.close();
}

void LBM::read_Geo(std::string filename, int mpi_rank)
{
    int x0, y0, z0;
    force = new vector3<double>[1]();
    //for (int k = 0; k < direction_size; k++)
    //{
    //    cout << "(" << ei[k].x << " " << ei[k].y << "  " << ei[k].z << ")" << '\n';
    //}

    char buffer[256];
    int ele_point = 0, total_ibpoint = 0, *a;

    double x_min, x_max, y_min, y_max, z_min, z_max;
    double x_length, y_length, z_length;
    double *lagx, *lagy, *lagz;
    int *Ep0, *Ep1, *Ep2;

    a = new int[2]();

    ifstream infile(filename);
    //ofstream outfile("out.plt");    //! 输出文件
    if (!infile)
    {
        cout << "Unable to open input complex geometry (shape) file:mesh_3d.dat or mesh_2d.dat";
        exit(1);    //! terminate with error
    }
    //if (!outfile)
    //{
    //    cout << "Unable to open output file";
    //    exit(1); // terminate with error
    //}
    infile.getline(buffer, 60);    //! 跳过
    infile.getline(buffer, 60);    //! 跳过
    infile.getline(buffer, 90);    //! 读入结点个数、单元个数, !before numbers there should be at least a blank
    string str1 = buffer;
    //cout << str1 << '\n';
    //sscanf_s(buffer, "%*s, %*s %d", &total_ibpoint); // &total_ibpoint,&ele_point);&buf, _countof(buf)
    //cout<<total_ibpoint<<endl;
    getNumberFromString(str1, a, 2);
    total_ibpoint = a[0];
    ele_point = a[1];

    //cout << "total_ibpoint=" << total_ibpoint << "    ele_point=" << ele_point << '\n';

    lagx = new double[total_ibpoint]();
    lagy = new double[total_ibpoint]();
    lagz = new double[total_ibpoint]();

    Ep0 = new int[ele_point]();
    Ep1 = new int[ele_point]();
    Ep2 = new int[ele_point]();

    while (!infile.eof())    //! 逐行读入参数， 以文件结束处为标志，若文件最后多出一行空白行或一行未预料内容，则读入可能会有问题
    {
        for (int k = 0; k < total_ibpoint; k++)
        {
            infile.getline(buffer, 100);    //! 60太小，不行
            //string str1 = buffer;
            //cout << str1 << '\n';
            sscanf_s(buffer, "%lf  %lf  %lf", &lagx[k], &lagy[k], &lagz[k]);
            //cout << lagx[k] << "     " << lagy[k] << "  " << lagz[k] << '\n';
        }

        for (int k = 0; k < ele_point; k++)
        {
            infile.getline(buffer, 40);
            //                string str1 = buffer;
            //                cout << str1 << '\n';
            sscanf_s(buffer, "%d       %d       %d", &Ep0[k], &Ep1[k], &Ep2[k]);
            //cout << Ep0[k] << "     " << Ep1[k] << "  " << Ep2[k] << '\n';
        }
    }

    //! 为了检查一遍，重新输出Finite element面网格文件
    //std::ofstream output_stream;
    //output_stream.open("out.plt", std::ofstream::out);
    ////若打开时有option：    | std::ofstream::app，则 app追加文件，不加为默认覆盖
    //if (true) {
    //    output_stream << "Title = " << '"' << "finite-element data" << '"' << '\n';
    //    output_stream << "variables = x, y, z" << '\n';
    //    output_stream << "zone n=" << total_ibpoint << ", E=" << ele_point << ",f=FEpoint,et = triangle" << '\n';
    //}

    //for (int k = 0; k < total_ibpoint; k++) {
    //    output_stream << lagx[k] << "  " << lagy[k] << "  " << lagz[k] << '\n';
    //}
    //for (int k = 0; k < ele_point; k++) {
    //    output_stream << Ep0[k] << "  " << Ep1[k] << "  " << Ep2[k] << "  " << '\n';
    //}
    //output_stream.close();
    //-----------------------------------------------------
    x_max = *max_element(lagx, lagx + total_ibpoint);
    x_min = *min_element(lagx, lagx + total_ibpoint);
    y_max = *max_element(lagy, lagy + total_ibpoint);
    y_min = *min_element(lagy, lagy + total_ibpoint);
    z_max = *max_element(lagz, lagz + total_ibpoint);
    z_min = *min_element(lagz, lagz + total_ibpoint);

    //cout << "x_max=" << x_max << "  " << "xmin=" << x_min << '\n';
    //cout << "y_max=" << y_max << "  " << "ymin=" << y_min << '\n';
    //cout << "z_max=" << z_max << "  " << "zmin=" << z_min << '\n';
    x_length = x_max - x_min;
    y_length = y_max - y_min;
    z_length = z_max - z_min;
    //std::cout << "xlength=" << x_length << "   " << "ylength=" << y_length << "   " << "zlength=" << z_length << '\n';
    //std::cout << ac << "   " << "t=" << t << "   " << "u=" << u << "   " << "v=" << v << '\n';

    bool judge;
    vector3<double> v0, v1, v2;
    double t, u, v;
    int *t_rcd;
    t_rcd = new int[2]();
    int zmin0 = (int)floor(z_min);
    int zmax0 = (int)floor(z_max) + 2;    //! 3D
    if (direction_size == 9)    //! 2D
    {
        zmin0 = 0; zmax0 = 1;
    }

    for (int x = (int)floor(x_min); x < (int)floor(x_max) + 2; x++)
    {
        for (int y = (int)floor(y_min); y < (int)floor(y_max) + 2; y++)
        {
            //for (int z = (int)floor(z_min); z < (int)floor(z_max) + 2; z++) {
            for (int z = zmin0; z < zmax0; z++)
            {
                double origx = x; double origy = y; double origz = z;
                vector3<double> orig = { origx,origy,origz }, dir = { 0.0, 0.0, 1.0 };
                if (direction_size == 9) dir = { 0.0, 1.0, 0.0 };
                t_rcd[0] = 0;
                t_rcd[1] = 0;
                for (int k = 0; k < ele_point; k++)
                {
                    v0 = { lagx[Ep0[k] - 1],lagy[Ep0[k] - 1],lagz[Ep0[k] - 1] };
                    v1 = { lagx[Ep1[k] - 1],lagy[Ep1[k] - 1],lagz[Ep1[k] - 1] };
                    v2 = { lagx[Ep2[k] - 1],lagy[Ep2[k] - 1],lagz[Ep2[k] - 1] };
                    if (direction_size == 9)
                    {   //! 2D
                        v0 = { lagx[Ep0[k] - 1],lagy[Ep0[k] - 1],1.0 };
                        v1 = { lagx[Ep1[k] - 1],lagy[Ep1[k] - 1],0.0 };
                        v2 = { lagx[Ep2[k] - 1],lagy[Ep2[k] - 1], -1.0 };
                    }
                    //! 注意：每个三角单元的顶点（结点）编号Ep0[k]，Ep1[k]，Ep2[k]都是从1开始的，但是我们的结点坐标数组lagx[]是从0开始的。
                    //! 因此特别注意：正确得到Ep0[k]的x坐标值应该采用lagx[Ep0[k]-1]
                    judge = IntersectTriangle(orig, dir, v0, v1, v2, &t, &u, &v);
                    if (judge == true)
                    {   //! && u > 0.00001 && v > 0.00001
                        if (t > 0)
                        {
                            t_rcd[0] ++;    //! t为正值的次数
                        }
                        else
                        {
                            t_rcd[1] ++;    //! t为负值的次数
                        }
                    }

                    if (t_rcd[0] > 0 && t_rcd[1] > 0)
                    {
                        //std::cout << "p" << t_rcd[0] << "  n" << t_rcd[1] << '\n'; 
                        x0 = x - Xori;
                        y0 = y - Yori;
                        z0 = z - Zori;
                        if (x0 >= -2 && x0 < NX + 2 && y0 >= -2 && y0 < NY + 2 && z0 >= -2 && z0 < NZ + 2)
                        {
                            obst[index(x0, y0, z0, NX, NY)] = 1;
                        }
                    }
                }
            }
        }
    }

    //输出  obst=1  tecplot观看
    //std::ofstream output1_stream;
    //output1_stream.open("del" + std::to_string(mpi_rank) + ".plt", std::ofstream::out);
    ////若打开时有option：    | std::ofstream::app，则 app追加文件，不加为默认覆盖
    //if (true) {
    //    output1_stream << "variables = x, y, z, obst" << '\n';
    //    output1_stream << "zone i=" << NX << ",j=" << NY << ", k=" << NZ << ", f=point" << '\n';
    //}

    //for (int z = 0; z < NZ; z++) {
    //    for (int y = 0; y < NY; y++) {
    //        for (int x = 0; x < NX; x++) {
    //            output1_stream << x << " " << y << " " << z << " " <<
    //                obst[index(x, y, z, NX, NY)] << '\n';
    //        }
    //    }
    //}
    //output1_stream.close();

    int m = 0;
    //! 最靠近流体的固体点obst值改成2
    for (int x = -2; x < NX + 2; x++)
    {
        for (int y = -2; y < NY + 2; y++)
        {
            for (int z = -2; z < NZ + 2; z++)
            {
                if (obst[index(x, y, z, NX, NY)] == 1)
                {
                    int n = 0;
                    for (int k = 0; k < direction_size; k++)
                    {
                        int ind = index(x + ei[k].x, y + ei[k].y, z + ei[k].z, NX, NY);
                        if (ind >= 0 && ind <= index(NX + 1, NY + 1, NZ + 1))
                        {
                            if (obst[ind] == 0) { n++; }
                        }
                    }
                    if (n > 0) { obst[index(x, y, z, NX, NY)] = 2;   m++; }
                }
            }
        }
    }   //! 目的是得到obst2的结点数  以确定开数组的大小

    num_node2 = m;
    node2 = new int[num_node2]();
    delta = new double[num_node2 * direction_size]();

    m = 0;
    for (int x = -1; x < NX + 1; x++)
    {
        for (int y = -1; y < NY + 1; y++)
        {
            for (int z = -1; z < NZ + 1; z++)
            {
                if (obst[index(x, y, z, NX, NY)] == 2)
                {
                    node2[m] = index(x, y, z, NX, NY);
                    for (int k1 = 0; k1 < direction_size; k1++)
                    {
                        int ip = ei[k1].x, jp = ei[k1].y, kp = ei[k1].z;
                        int ind = index(x + ip, y + jp, z + kp, NX, NY);
                        bool judge;
                        double t, u, v;
                        vector3<double> v0, v1, v2;
                        vector3<double> orig = { (double)x,(double)y,(double)z },
                            dir = { (double)ip, (double)jp, (double)kp };
                        double b1 = 0., a1 = dir.norm_square();    //! 速度模型中速度长度 1, sqrt2, sqrt3...
                        if (ind >= 0 && ind <= index(NX + 1, NY + 1, NZ + 1))
                        {
                            if (obst[ind] == 0)
                            {
                                for (int k = 0; k < ele_point; k++)
                                {
                                    v0 = { lagx[Ep0[k] - 1],lagy[Ep0[k] - 1],lagz[Ep0[k] - 1] };
                                    v1 = { lagx[Ep1[k] - 1],lagy[Ep1[k] - 1],lagz[Ep1[k] - 1] };
                                    v2 = { lagx[Ep2[k] - 1],lagy[Ep2[k] - 1],lagz[Ep2[k] - 1] };
                                    if (direction_size == 9)
                                    {
                                        v0 = { lagx[Ep0[k] - 1],lagy[Ep0[k] - 1],1.0 };
                                        v1 = { lagx[Ep1[k] - 1],lagy[Ep1[k] - 1],0.0 };
                                        v2 = { lagx[Ep2[k] - 1],lagy[Ep2[k] - 1],-1.0 };
                                    }
                                    judge = IntersectTriangle(orig, dir, v0, v1, v2, &t, &u, &v);
                                    if (judge == true && t > 0)    //! && t<b1
                                    {
                                        b1 = t;
                                    }
                                }
                                if (b1 / a1 > 1) delta[m * direction_size + k1] = 0.;
                                else delta[m * direction_size + k1] = 1.0 - b1 / sqrt(a1);
                                //cout <<"q="<< 1.0- b1/a1<<"    " <<"a1=" <<a1 <<'\n';
                                //if (b1 / a1 > 1) cout <<"b1="<<b1<<" a1="<<a1<<" error!!!!!!"<< '\n';
                            }
                        }
                    }
                    m++;
                }
            }
        }
    }

    delete [] lagx;    lagx = nullptr;
    delete [] lagy;    lagy = nullptr;
    delete [] lagz;    lagz = nullptr;
    delete [] Ep0;    Ep0 = nullptr;
    delete [] Ep1;    Ep1 = nullptr;
    delete [] Ep2;    Ep2 = nullptr;
    delete [] t_rcd;    t_rcd = nullptr;
}

void LBM::getNumberFromString(string s, int *a, int n)
{
    stringstream str_strm;
    str_strm << s;    //! convert the string s into stringstream
    string temp_str;
    int temp_int;
    int k = 0;
    while (!str_strm.eof())
    {
        str_strm >> temp_str;    //! take words into temp_str one by one
        if (stringstream(temp_str) >> temp_int)
        {   //! try to convert string to int
            a[k] = temp_int;
            //cout << temp_int << endl;
            k++;
        }
        temp_str = "";    //! clear temp string
    }
}

bool LBM::IntersectTriangle(vector3<double> &orig, vector3<double> &dir,
    vector3<double> &v0, vector3<double> &v1, vector3<double> &v2,
    double *t, double *u, double *v)
{
    //! E1
    vector3<double> E1 = v1 - v0;

    //! E2
    vector3<double> E2 = v2 - v0;

    //! P
    vector3<double> P = dir.cross(E2);

    //! determinant
    double det = E1.dot(P);

    //! keep det > 0, modify T accordingly
    vector3<double> T;
    if (det > 0)
    {
        T = orig - v0;
    }
    else
    {
        T = v0 - orig;
        det = -det;
    }

    //! If determinant is near zero, ray lies in plane of triangle
    if (det < 0.0001f)
        return false;

    //! Calculate u and make sure u <= 1
    *u = T.dot(P);
    if (*u < 0.0f || *u > det)
        return false;

    //! Q
    vector3<double> Q = T.cross(E1);

    //! Calculate v and make sure u + v <= 1
    *v = dir.dot(Q);
    if (*v < 0.0f || *u + *v > det)
        return false;

    //! Calculate t, scale parameters, ray intersects triangle
    *t = E2.dot(Q);

    double fInvDet = 1.0f / det;
    *t *= fInvDet;
    *u *= fInvDet;
    *v *= fInvDet;

    return true;
}

void LBM::write_binary(std::string filename)    //! int im, int jm, int km
{
    int VarN = 7;    //! Number of variables(VarN) 
    //char filename[100];
    //sprintf_s(filename, 100, "%s", "P.plt");  //sprintf(filename, "P.plt");
    FILE *fp;
    errno_t err;
    if ((err = fopen_s(&fp, filename.c_str(), "wb")) != 0) std::cout << "error to open file: " << filename << '\n';
    int im = this->total_NX;
    int jm = this->total_NY;
    int km = this->total_NZ;

    char Title[] = "TubeFlow";
    char Vname[7][5] = { "X","Y","Z","u","v","w","p" };    //! "XYZp";

    char Zonename1[] = "Zone 001";
    float ZONEMARKER = 299.0;
    float EOHMARKER = 357.0;
    //! ==============Header Secontion =================//
    //! ------1.1 Magic number, Version number
    char MagicNumber[] = "#!TDV101";
    fwrite(MagicNumber, 8, 1, fp);

    //! ---- - 1.2.Integer value of 1.-----------------
    int IntegerValue = 1;
    fwrite(&IntegerValue, sizeof(IntegerValue), 1, fp);

    //! ---- - 1.3.Title and variable names.------------
    //! ---- - 1.3.1.The TITLE.
    dumpstring(Title, fp);
    fwrite(&VarN, sizeof(VarN), 1, fp);

    //! ------1.3.3 Variable names.
    for (int i = 0; i < VarN; i++)
    {
        dumpstring(Vname[i], fp);
    }

    //! ---- - 1.4.Zones----------------------------------
    //! --------Zone marker.Value = 299.0
    fwrite(&ZONEMARKER, 1, sizeof(ZONEMARKER), fp);

    //! --------Zone name.
    //! fwrite(Zonename1, sizeof(Zonename1), 1, fp);
    dumpstring(Zonename1, fp);

    int ZoneStuff[5] = { -1,0,1,0,0 };
    //! Zone color
    //! ZoneType
    //! DaraPacking 0=Block, 1=Point
    //! Specify Var Location. 0 = Don't specify, all c_str is located at the nodes. 1 = Specify
    //! Number of user defined face neighbor connections(value >= 0)
    for (int i = 0; i < 5; i++)
    {
        fwrite(&ZoneStuff[i], sizeof(int), 1, fp);
    }

    //! im, jm, km
    fwrite(&im, sizeof(im), 1, fp);
    fwrite(&jm, sizeof(jm), 1, fp);
    fwrite(&km, sizeof(km), 1, fp);

    //! -1 = Auxiliary name / value pair to follow 0 = No more Auxiliar name / value pairs.
    int AuxiliaryName = 0;
    fwrite(&AuxiliaryName, sizeof(AuxiliaryName), 1, fp);

    //! ----I HEADER OVER----------------------------------------------
    //! =============================Geometries section================
    //! =============================Text section======================
    //!  EOHMARKER, value = 357.0
    fwrite(&EOHMARKER, sizeof(EOHMARKER), 1, fp);

    //! ================II.Data section===============//
    //! ---- - 2.1 zone------------------------------------------------
    fwrite(&ZONEMARKER, sizeof(ZONEMARKER), 1, fp);

    //! variable format, 1 = Float, 2 = Double, 3 = LongInt, 4 = ShortInt, 5 = Byte, 6 = Bit
    int VariableType[7] = { 3,3,3,2,2,2,2 };
    for (int i = 0; i < VarN; i++)
    {
        fwrite(&VariableType[i], sizeof(int), 1, fp);
    }

    //! Has variable sharing 0 = no, 1 = yes.
    int HasVarSharing = 0;
    fwrite(&HasVarSharing, sizeof(HasVarSharing), 1, fp);

    //! Zone number to share connectivity list with(-1 = no sharing).
    int ZoneNumToShareConnectivity = -1;
    fwrite(&ZoneNumToShareConnectivity, sizeof(ZoneNumToShareConnectivity), 1, fp);
    //! Zone c_str.Each variable is in c_str format asspecified above.

    for (int k = 0; k < km; k++)
    {
        for (int j = 0; j < jm; j++)
        {
            for (int i = 0; i < im; i++)
            {
                int Va1 = i;
                int Va2 = j;
                int Va3 = k;
                fwrite(&Va1, sizeof(Va1), 1, fp);
                fwrite(&Va2, sizeof(Va2), 1, fp);
                fwrite(&Va3, sizeof(Va3), 1, fp);
                fwrite(&velTotal[this->total_index(i, j, k)].x, sizeof(double), 1, fp);
                fwrite(&velTotal[this->total_index(i, j, k)].y, sizeof(double), 1, fp);
                fwrite(&velTotal[this->total_index(i, j, k)].z, sizeof(double), 1, fp);
                fwrite(&rhoTotal[this->total_index(i, j, k)], sizeof(double), 1, fp);
            }
        }
    }
    fclose(fp);
    //getchar();
}

void dumpstring(char *str, FILE *file)
{
    int value = 0;
    while ((*str) != '\0')
    {
        value = (int)*str;
        fwrite(&value, sizeof(int), 1, file);
        str++;
    }
    char null_char[] = "";
    value = (int)*null_char;
    fwrite(&value, sizeof(int), 1, file);
}