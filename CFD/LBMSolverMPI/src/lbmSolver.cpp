#include "LBMSolverMPI.hpp"
#include "vector3.cpp"
#include <iomanip>
#include <iostream>
#include <cmath>
#include <mpi.h>
using namespace std;

LBM::LBM()
{

}

LBM::~LBM()
{
    delete [] rho;    rho = nullptr;
    delete [] pdf;    pdf = nullptr;
    delete [] pdf2;    pdf2 = nullptr;
    delete [] n_myid;    n_myid = nullptr;
    delete [] obst;    obst = nullptr;
    delete [] vel;    vel = nullptr;
    delete [] weights;    weights = nullptr;
    delete [] reverse_indexes;    reverse_indexes = nullptr;
    delete [] ei;    ei = nullptr;
    if (iMRT == 1)
    {
        delete [] M_MRT;    M_MRT = nullptr;
        delete [] M_MRTI;    M_MRTI = nullptr;
        delete [] MSM;    MSM = nullptr;
    }
    if (iGEO == 1)
    {
        delete [] force;    force = nullptr;
        delete [] node2;    node2 = nullptr;
        delete [] delta;    delta = nullptr;
    }
}

void LBM::ReadParameter(Params_Pakage* p, int mpi_rank)
{
    char buffer[256];
    std::string A;    //! 字符串
    char conti_f[20], geo_f[20];
    vector3<double> aa;
    double rho_specified;
    ifstream infile(p->in);    //! 输入文件
    ofstream outfile(p->out);    //! 输出文件

    if (!infile)
    {
        cout << "Unable to open input file";
        exit(1);    //! terminate with error
    }
    if (!outfile)
    {
        cout << "Unable to open output file";
        exit(1);    //! terminate with error
    }
    //! sscanf_s(tokenstring, "%s", s, _countof(s));
    //! sscanf_s(tokenstring, "%c", &c, sizeof(char));ggggg
    while (!infile.eof())    //! 逐行读入参数， 以文件结束处为标志，若文件最后多出一行空白行或一行未预料内容，则读入可能会有问题
    {
        infile.getline(buffer, 100);    //! 读入标志是否续算，是为1，否为0，续算则继续读入数据输入文件名
        sscanf_s(buffer, "%d", &this->iMRT);
        infile.getline(buffer, 120);    //! 读入标志是否续算，是为1，否为0，续算则继续读入数据输入文件名
        sscanf_s(buffer, "%d", &this->iContinue);
        if (this->iContinue == 1)
        {
            sscanf_s(buffer, "%d %s", &this->iContinue, &conti_f, (unsigned int)sizeof(conti_f));
            conti_file = conti_f;
        }
        infile.getline(buffer, 100);    //! 读入标志是否大涡模拟，是为1，否为0
        sscanf_s(buffer, "%d", &this->iLES);
        infile.getline(buffer, 120);    //! 读入标志是否需要输入复杂外形物体，是为1，否为0，需要则继续读入图形tecplot输入文件名
        sscanf_s(buffer, "%d", &this->iGEO);
        if (this->iGEO == 1)
        {
            sscanf_s(buffer, "%d %s", &this->iGEO, &geo_f, (unsigned int)sizeof(geo_f));
            geo_file = geo_f;
        }
        infile.getline(buffer, 100);    //! 读入标志是否多块网格，是为1，否为0
        sscanf_s(buffer, "%d", &this->FINE);
        infile.getline(buffer, 100);    //! 读入标志是否需要收敛判据，是为1，否为0，若需要则继续读入收敛阈值
        sscanf_s(buffer, "%d", &this->CONV);
        if (this->CONV == 1) sscanf_s(buffer, "%d %lf", &this->CONV, &this->conv_v);

        infile.getline(buffer, 60);    //! 读入计算域 X 方向网格数
        sscanf_s(buffer, "%d", &this->total_NX);
        infile.getline(buffer, 60);    //! 计算域 Y 方向网格数
        sscanf_s(buffer, "%d", &this->total_NY);
        infile.getline(buffer, 60);    //! 计算域 Z 方向网格数
        sscanf_s(buffer, "%d", &this->total_NZ);

        infile.getline(buffer, 60);    //! 读入密网格 X 方向网格数
        sscanf_s(buffer, "%d", &this->fine_NX);
        infile.getline(buffer, 60);    //! 密网格Y 方向网格数
        sscanf_s(buffer, "%d", &this->fine_NY);
        infile.getline(buffer, 60);    //! 密网格 Z 方向网格数
        sscanf_s(buffer, "%d", &this->fine_NZ);
        infile.getline(buffer, 60);    //! 密网格原点坐标
        sscanf_s(buffer, "%d", &this->LefLowx);
        infile.getline(buffer, 60);
        sscanf_s(buffer, "%d", &this->LefLowy);
        infile.getline(buffer, 60);
        sscanf_s(buffer, "%d", &this->LefLowz);

        infile.getline(buffer, 60);    //! 计算域 X 方向核数
        sscanf_s(buffer, "%d", &this->x_np);
        infile.getline(buffer, 60);    //! 计算域 Y 方向核数
        sscanf_s(buffer, "%d", &this->y_np);
        infile.getline(buffer, 60);    //! 计算域 Z 方向核数
        sscanf_s(buffer, "%d", &this->z_np);

        infile.getline(buffer, 60);    //! 每隔多少步输出一次文件
        sscanf_s(buffer, "%d", &this->FrameRate);
        infile.getline(buffer, 60);    //! 总共多少步终止程序
        sscanf_s(buffer, "%d", &this->ttstep);
        infile.getline(buffer, 60);    //! 初始密度
        sscanf_s(buffer, "%lf", &this->density0);
        infile.getline(buffer, 60);    //! 松弛因子，与运动粘性关系:\nu=c_s^2*(\tau-0.5)\delta t
        sscanf_s(buffer, "%lf", &this->tau);
        infile.getline(buffer, 100);    //! 体积力
        sscanf_s(buffer, "(%lf, %lf, %lf)", &this->gx, &this->gy, &this->gz);
        infile.getline(buffer, 60);    //! LBM速度模型
        sscanf_s(buffer, "%s", &p->velocity_set, (unsigned int)sizeof(p->velocity_set));
        this->velocity_set = p->velocity_set;
        //! 字符串提取时能默认以空格分割
        infile.getline(buffer, 60);    //! 边界条件类型
        sscanf_s(buffer, "%s", &p->bc, (unsigned int)sizeof(p->bc));
        this->boundary_condition = p->bc;
        infile.getline(buffer, 256);    //! 文件中的说明行，不取信息，跳过！
        infile.getline(buffer, 256);    //! 文件中的说明行，不取信息，跳过！
        infile.getline(buffer, 100);
        sscanf_s(buffer, "%s", &this->BC.xminFace, (unsigned int)sizeof(this->BC.xminFace));
        A = this->BC.xminFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            this->BC.vel_specified[0] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            this->BC.rho_specified[0] = rho_specified;
        }

        infile.getline(buffer, 100);
        sscanf_s(buffer, "%s", &this->BC.xmaxFace, (unsigned int)sizeof(this->BC.xmaxFace));
        A = this->BC.xmaxFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            this->BC.vel_specified[1] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            this->BC.rho_specified[1] = rho_specified;
        }

        infile.getline(buffer, 100);
        sscanf_s(buffer, "%s", &this->BC.yminFace, (unsigned int)sizeof(this->BC.yminFace));
        A = this->BC.yminFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            this->BC.vel_specified[2] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            this->BC.rho_specified[2] = rho_specified;
        }

        infile.getline(buffer, 256);
        sscanf_s(buffer, "%s", &this->BC.ymaxFace, (unsigned int)sizeof(this->BC.ymaxFace));
        A = this->BC.ymaxFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            this->BC.vel_specified[3] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            this->BC.rho_specified[3] = rho_specified;
        }

        infile.getline(buffer, 100);
        sscanf_s(buffer, "%s", &this->BC.zminFace, (unsigned int)sizeof(this->BC.zminFace));
        A = this->BC.zminFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            this->BC.vel_specified[4] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            this->BC.rho_specified[4] = rho_specified;
        }

        infile.getline(buffer, 100);
        sscanf_s(buffer, "%s", &this->BC.zmaxFace, (unsigned int)sizeof(this->BC.zmaxFace));
        A = this->BC.zmaxFace;
        if (A == "velocity") {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            this->BC.vel_specified[5] = aa;
        }
        if (A == "pressure") {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            this->BC.rho_specified[5] = rho_specified;
        }
        infile.getline(buffer, 100);//跳一行
        infile.getline(buffer, 100);
        sscanf_s(buffer, "(%lf, %lf, %lf)", &vel0.x, &vel0.y, &vel0.z);

        if (mpi_rank == 0)
        {
            cout << "Computational domain size:   " << total_NX << "  " << total_NY << "   " << total_NZ << endl;
            cout << "number of timesteps per frame: " << FrameRate << endl;
            cout << "The velocity model used:   " << velocity_set << endl;
            cout << "Tau: " << tau << '\n';
            cout << "G: (" << gx <<","<< gy <<"," << gz <<")" << '\n';
            cout << boundary_condition << '\n' << endl;
            //! 参数检查
            if (tau < 0.5)
            {
                std::cout << "Error: Tau must be greater than 0.5 for numerical stability." << '\n';
                exit(1);
            }
            if (velocity_set != "D3Q15" && velocity_set != "D3Q19" && velocity_set != "D3Q27" && velocity_set != "D2Q9")
            {
                std::cout << "Error: Please specify a valid velocity set such as D3Q15,D3Q27,D3Q19 or D2Q9." << '\n';
                exit(1);
            }
            if (total_NZ != 1 && velocity_set == "D2Q9")
            {
                std::cout << "Warning: NZ=1 for D2Q9.";
                exit(1);
            }
        }
    }

    infile.close();

    //! 读完后同时输出以上参数到另一文件，以验证读取的正确性
    outfile << "MRT:" << "   " << iMRT << endl;
    outfile << "CONTI:" << "   " << iContinue << endl;
    outfile << "LES:" << "   " << iLES << endl;
    outfile << "GEO:" << "   " << iGEO << endl;
    outfile << "MB:" << "   " << FINE << endl;
    outfile << "total_NX:" << "   " << total_NX << endl;
    outfile << "total_NY:" << "   " << total_NY << endl;
    outfile << "total_NZ:" << "   " << total_NZ << endl;
    outfile << "fine_NX:" << "   " << fine_NX << endl;
    outfile << "fine_NY:" << "   " << fine_NY << endl;
    outfile << "fine_NZ:" << "   " << fine_NZ << endl;
    outfile << "LefLowx:" << "   " << LefLowx << endl;
    outfile << "LefLowy:" << "   " << LefLowy << endl;
    outfile << "LefLowz:" << "   " << LefLowz << endl;
    outfile << "x_np:" << "   " << x_np << endl;
    outfile << "y_np:" << "   " << y_np << endl;
    outfile << "z_np:" << "   " << z_np << endl;
    outfile << "FrameRate:" << "   " << FrameRate << endl;
    outfile << "Total steps:" << " " << ttstep << endl;
    outfile << "LBM.initial density:" << " " << density0 << endl;
    outfile << "LBM.tau:" << " " << tau << endl;
    outfile << "LBM.gx:" << " " << gx << endl;
    outfile << "LBM.gy:" << " " << gy << endl;
    outfile << "LBM.gz:" << " " << gz << endl;
    outfile << "velocity model:" << " " << velocity_set << endl;
    outfile << "Boundary condition type:" << " " << boundary_condition << endl;
    outfile << "xmin Face BC:  " << " " << BC.xminFace << endl;
    outfile << "xmax Face BC:  " << " " << BC.xmaxFace << endl;
    outfile << "ymin Face BC:  " << " " << BC.yminFace << endl;
    outfile << "ymax Face BC:  " << " " << BC.ymaxFace << endl;
    outfile << "zmin Face BC:  " << " " << BC.zminFace << endl;
    outfile << "zmax Face BC:  " << " " << BC.zmaxFace << endl;
    outfile << "initial vel  " << " " << vel0.x << " " << vel0.y << " " << vel0.z << endl;
    outfile.close();
};

int LBM::initialise(Params_Pakage *p, int mpi_rank)
{
    ReadParameter(p, mpi_rank);
    mpi_ini(mpi_rank);

    set_velocity_set();
    iniFlow();

    if (iMRT == 1)
    {
        iniMRT();
    }

    if (iGEO == 1)
    {
        read_Geo(geo_file, mpi_rank);
    }

    iniobst();

    if (iContinue == 1)
    {
        read_all(conti_file);
    }
    else
    {
        //cout << "Loaded default initial conditions" << endl;
        double density = { 1 }, u_x = vel0.x, u_y = vel0.y, u_z = vel0.z;    //! 缺省初始化
        density = density0;
        for (int i = 0; i < NX; i++)
        {
            for (int j = 0; j < NY; j++)
            {
                for (int k = 0; k < NZ; k++)
                {
                    set_density(i, j, k, density);
                    set_velocity(i, j, k, u_x, u_y, u_z);
                }
            }
        }
        //std::cout << "If you wish to use your own initial conditions, please run the program but with command: generate_ic as a argument which will output ic.csv in format of p,u_x,u_y,u_z, assume indexes are incrementing i,j,k for i<NX,j<NY and k<NZ" << '\n';
    }

    //cout << "Kinematic shear viscosity: " << viscosity << '\n';
    //cout << "For velocity set D2Q9,D3Q15,D3Q19 and D3Q27, |u_max|<0.577\n";
    //output_tecplot("output/rank" + std::to_string(mpi_rank) + "_" + "0.plt");

    if (mpi_rank == 0)
    {
        string command;
        command = "mkdir output";
        system(command.c_str());
    }
    return 0;
}

void LBM::iniFlow()
{
    int NumNodes = (NX + 4) * (NY + 4) * (NZ + 4);
    int distributions_memo_size = NumNodes * direction_size;
    rho = new double[NumNodes];
    obst = new int[NumNodes];
    vel = new vector3<double>[NumNodes];
    pdf2 = new double[distributions_memo_size];
    pdf = new double[distributions_memo_size];
    for (int i = 0; i < NumNodes; i++)
    {
        rho[i] = 1.0;
        obst[i] = 0;
    }
    for (int x = -2; x < NX + 2; x++)
    {
        for (int y = -2; y < NY + 2; y++)
        {
            for (int z = -2; z < NZ + 2; z++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                    pdf2[index(x, y, z, i)] = weights[i];
                    pdf[index(x, y, z, i)] = weights[i];
                }
            }
        }
    }
}

void LBM::iniobst()
{
    //! 初始化obst
    string Astring;
    if (x_id == 0)
    {
        int x = 0;
        Astring = BC.xminFace;    //! 字符数组赋值给string，为字符串比较准备
        if (Astring == "nonslip")
        {
            for (int y = 0; y < NY; y++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    obst[index(x, y, z)] = 3;
                }
            }
        }
        else if (Astring == "velocity")
        {
            for (int y = 0; y < NY; y++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    obst[index(x, y, z)] = -1;
                }
            }
        }
        else if (Astring == "pressure")
        {
            for (int y = 0; y < NY; y++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    obst[index(x, y, z)] = -2;
                }
            }
        }
    }

    if (x_id == x_np - 1)
    {
        int x = NX - 1;
        Astring = BC.xmaxFace;
        if (Astring == "nonslip")
        {
            for (int y = 0; y < NY; y++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    obst[index(x, y, z)] = 3;
                }
            }
        }
        else if (Astring == "velocity")
        {
            for (int y = 0; y < NY; y++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    obst[index(x, y, z)] = -1;
                }
            }
        }
        else if (Astring == "pressure")
        {
            for (int y = 0; y < NY; y++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    obst[index(x, y, z)] = -2;
                }
            }
        }
    }

    if (y_id == 0)
    {
        int y = 0;
        Astring = BC.yminFace;
        if (Astring == "nonslip")
        {
            for (int x = 0; x < NX; x++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    obst[index(x, y, z)] = 3;
                }
            }
        }
        else if (Astring == "velocity")
        {
            for (int x = 0; x < NX; x++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    obst[index(x, y, z)] = -1;
                }
            }
        }
        else if (Astring == "pressure")
        {
            for (int x = 0; x < NX; x++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    obst[index(x, y, z)] = -2;
                }
            }
        }
    }

    if (y_id == y_np - 1)
    {
        int y = NY - 1;
        Astring = BC.ymaxFace;
        if (Astring == "nonslip")
        {
            for (int x = 0; x < NX; x++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    obst[index(x, y, z)] = 3;
                }
            }
        }
        else if (Astring == "velocity")
        {
            for (int x = 0; x < NX; x++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    obst[index(x, y, z)] = -1;
                }
            }
        }
        else if (Astring == "pressure")
        {
            for (int x = 0; x < NX; x++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    obst[index(x, y, z)] = -2;
                }
            }
        }
    }

    if (z_id == 0)
    {
        int z = 0;
        Astring = BC.zminFace;
        if (Astring == "nonslip")
        {
            for (int x = 0; x < NX; x++)
            {
                for (int y = 0; y < NY; y++)
                {
                    obst[index(x, y, z)] = 3;
                }
            }
        }
        else if (Astring == "velocity")
        {
            for (int x = 0; x < NX; x++)
            {
                for (int y = 0; y < NY; y++)
                {
                    obst[index(x, y, z)] = -1;
                }
            }
        }
        else if (Astring == "pressure")
        {
            for (int x = 0; x < NX; x++)
            {
                for (int y = 0; y < NY; y++)
                {
                    obst[index(x, y, z)] = -2;
                }
            }
        }
    }

    if (z_id == z_np - 1)
    {
        int z = NZ - 1;
        Astring = BC.zmaxFace;
        if (Astring == "nonslip")
        {
            for (int x = 0; x < NX; x++)
            {
                for (int y = 0; y < NY; y++)
                {
                    obst[index(x, y, z)] = 3;
                }
            }
        }
        else if (Astring == "velocity")
        {
            for (int x = 0; x < NX; x++)
            {
                for (int y = 0; y < NY; y++)
                {
                    obst[index(x, y, z)] = -1;
                }
            }
        }
        else if (Astring == "pressure")
        {
            for (int x = 0; x < NX; x++)
            {
                for (int y = 0; y < NY; y++)
                {
                    obst[index(x, y, z)] = -2;
                }
            }
        }
    }
}

void LBM::set_velocity_set()
{
    if (velocity_set == "D3Q15")
    {
        direction_size = 15;
        send_direction_size = 5;
        ei = new vector3<int>[15] { {0, 0, 0}, { 1,0,0 }, { -1,0,0 }, { 0,1,0 }, { 0,-1,0 }, { 0,0,1 }, { 0,0,-1 },
            { 1,1,1 }, { -1,-1,-1 }, { 1,1,-1 }, { -1,-1,1 }, { 1,-1,1 }, { -1,1,-1 }, { -1,1,1 }, { 1,-1,-1 } };
        weights = new double[15] { 2.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
            1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0  };
    }
    else if (velocity_set == "D3Q19")
    {
        direction_size = 19;
        send_direction_size = 5;
        ei = new vector3<int>[19] { {0, 0, 0}, { 1,0,0 }, { -1,0,0 }, { 0,1,0 }, { 0,-1,0 }, { 0,0,1 }, { 0,0,-1 },
            { 1,1,0 }, { 1,-1,0 }, { -1,1,0 }, { -1,-1,0 }, { 1,0,1 }, { -1,0,1 }, { 1,0,-1 }, { -1,0,-1 }, { 0,1,1 },
            { 0,1,-1 }, { 0,-1,1 }, { 0,-1,-1 } };

        weights = new double[19] {1.0 / 3.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0,
            1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
            1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
            1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
    }
    else if (velocity_set == "D3Q27")
    {
        direction_size = 27;
        send_direction_size = 9;
        ei = new vector3<int>[27] { {0, 0, 0}, { 1,0,0 }, { 0,1,0 }, { -1,0,0 }, { 0,-1,0 }, { 0,0,1 }, { 0,0,-1 },
            { 1,1,0 }, { -1,1,0 }, { -1,-1,0 }, { 1,-1,0 }, { 1,0,1 }, { 0,1,1 }, { -1,0,1 }, { 0,-1,1 }, { 1,0,-1 }, { 0,1,-1 }, { -1,0,-1 }, { 0,-1,-1 },
            { 1,1,1 }, { -1,1,1 }, { -1,-1,1 }, { 1,-1,1 }, { 1,1,-1 }, { -1,1,-1 }, { -1,-1,-1 }, { 1,-1,-1 }};

        weights = new double[27] {8.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0,
            2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0
            , 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,
            1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0
            , 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0};
    }
    else if (velocity_set == "D2Q9")
    {
        direction_size = 9;
        send_direction_size = 3;
        ei = new vector3<int>[9] { {1, 0, 0}, { 0,1,0 }, { -1,0,0 }, { 0,-1,0 }, { 1,1,0 }, { -1,1,0 }, { -1,-1,0 }, { 1,-1,0 }, { 0,0,0 }};
        weights = new double[9] {1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 4.0 / 9.0};
    }
    lookup_reverse();
}

void LBM::lookup_reverse()
{
    reverse_indexes = new int[direction_size];
    for (int i = 0; i < direction_size; i++)
    {
        for (int j = 0; j < direction_size; j++)
        {
            if (ei[i].x == -ei[j].x && ei[i].y == -ei[j].y && ei[i].z == -ei[j].z)
            {
                reverse_indexes[i] = j;
            }
        }
    }
}

double LBM::calculate_feq(int x, int y, int z, int i)
{
    double dot_product = (double)vel[index(x, y, z)].x * (double)ei[i].x + (double)vel[index(x, y, z)].y * (double)ei[i].y +
        (double)vel[index(x, y, z)].z * (double)ei[i].z;

    double feq = weights[i] * rho[index(x, y, z)] * (1.0 + dot_product / (c_s * c_s) + dot_product * dot_product / (2 * c_s * c_s * c_s * c_s) - vel[index(x, y, z)].norm_square() / (2 * c_s * c_s));
    return feq;
}

double LBM::calculate_feq(int x, int y, int z, int i, double u_le_x)
{
    double dot_product = (vel[index(x, y, z)].x + u_le_x) * ei[i].x + vel[index(x, y, z)].y * ei[i].y +
        vel[index(x, y, z)].z * ei[i].z;

    double feq = weights[i] * rho[index(x, y, z)] * (1.0 + dot_product / (c_s * c_s) + dot_product * dot_product / (2 * c_s * c_s * c_s * c_s) - vel[index(x, y, z)].norm_square() / (2 * c_s * c_s));
    return feq;
}

double LBM::calculate_fg(int x, int y, int z, int i, double feq)
{
    double dot_product = gx * ((double)ei[i].x - (double)vel[index(x, y, z)].x) + gy * ((double)ei[i].y - (double)vel[index(x, y, z)].y) +
        gz * ((double)ei[i].z - (double)vel[index(x, y, z)].z);

    double c_s_squ = 1.0 / 3.0;
    double fg = (1.0 - 1.0 / (2 * tau)) / c_s_squ * dot_product * feq;
    return fg;
}

double LBM::calculate_fneq(int x, int y, int z, int i, double feq)
{
    double fneq = pdf[index(x, y, z, i)] - feq;
    return fneq;
}

double LBM::calculate_pij(double *fneq)
{
    double p_xx = 0.0;
    double p_yy = 0.0;
    double p_zz = 0.0;
    double p_xy = 0.0;
    double p_yz = 0.0;
    double p_zx = 0.0;
    for (int i = 0; i < direction_size; i++)
    {
        p_xx += (double)ei[i].x * (double)ei[i].x * fneq[i];
        p_yy += (double)ei[i].y * (double)ei[i].y * fneq[i];
        p_zz += (double)ei[i].z * (double)ei[i].z * fneq[i];
        p_xy += (double)ei[i].x * (double)ei[i].y * fneq[i];
        p_yz += (double)ei[i].y * (double)ei[i].z * fneq[i];
        p_zx += (double)ei[i].z * (double)ei[i].x * fneq[i];
    }
    double p_ij = sqrt(p_xx * p_xx + p_yy * p_yy + p_zz * p_zz + 2 * p_xy * p_xy + 2 * p_zx * p_zx + 2 * p_yz * p_yz);
    return p_ij;
}

void LBM::set_velocity(int x_field, int y_field, int z_field, double u_x, double u_y, double u_z)
{
    vel[index(x_field, y_field, z_field)].x = u_x;
    vel[index(x_field, y_field, z_field)].y = u_y;
    vel[index(x_field, y_field, z_field)].z = u_z;
}

void LBM::set_density(int x_field, int y_field, int z_field, double density)
{
    rho[index(x_field, y_field, z_field)] = density;
}

void LBM::compute_density_momentum(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
{
    for (int x = xmin; x < xmax; x++)
    {
        for (int y = ymin; y < ymax; y++)
        {
            for (int z = zmin; z < zmax; z++)
            {
                if (obst[index(x, y, z)] <= 0)
                {
                    double new_density = 0;
                    vector3<double> u;
                    for (int i = 0; i < direction_size; i++)
                    {
                        new_density += pdf[index(x, y, z, i)];
                        u.x += pdf[index(x, y, z, i)] * ei[i].x;
                        u.y += pdf[index(x, y, z, i)] * ei[i].y;
                        u.z += pdf[index(x, y, z, i)] * ei[i].z;
                    }
                    rho[index(x, y, z)] = new_density;
                    vel[index(x, y, z)].x = (u.x + 0.5 * gx) / new_density;
                    vel[index(x, y, z)].y = (u.y + 0.5 * gy) / new_density;
                    vel[index(x, y, z)].z = (u.z + 0.5 * gz) / new_density;
                }
                else
                {
                    rho[index(x, y, z)] = 1.0;
                    vel[index(x, y, z)].x = 0.0;
                    vel[index(x, y, z)].y = 0.0;
                    vel[index(x, y, z)].z = 0.0;
                }
            }
        }
    }
}

void LBM::stream(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
{
    for (int x = xmin; x < xmax; x++)
    {
        for (int y = ymin; y < ymax; y++)
        {
            for (int z = zmin; z < zmax; z++)
            {
                if (obst[index(x, y, z)] <= 0)
                {
                    for (int i = 0; i < direction_size; i++)
                    {
                        //! 先缺省是周期性的
                        int xmd = x - (int)ei[i].x;
                        int ymd = y - (int)ei[i].y;
                        int zmd = z - (int)ei[i].z;
                        pdf[index(x, y, z, i)] = pdf2[index(xmd, ymd, zmd, i)];
                    }
                }
            }
        }
    }
}

void LBM::collision(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
   {//! Performs the collision step.
    double tauinv = 1.0 / tau;
    double omtauinv = 1.0 - tauinv;    //! 1 - 1/tau
    if (iMRT == 0)
    {
        if (iLES == 1)
        {
            double *fneq = new double[direction_size];
            double *feq = new double[direction_size];
            double *fg = new double[direction_size];
            for (int x = xmin; x < xmax; x++)
            {
                for (int y = ymin; y < ymax; y++)
                {
                    for (int z = zmin; z < zmax; z++)
                    {
                        for (int i = 0; i < direction_size; i++)
                        {
                            feq[i] = calculate_feq(x, y, z, i);
                            fg[i] = calculate_fg(x, y, z, i, feq[i]);
                            fneq[i] = pdf[index(x, y, z, i)] - feq[i];
                        }
                        double p_ij = calculate_pij(fneq);
                        const double tau_t = 0.5 * (sqrt(tau * tau + 2 * sqrt(2) * c_s * c_s * p_ij / (c_s * c_s * c_s * c_s)) - tau);
                        const double omega = 1.0 / (tau + tau_t);
                        for (int i = 0; i < direction_size; i++)
                        {
                            pdf2[index(x, y, z, i)] = (1.0 - omega) * pdf[index(x, y, z, i)] + omega * feq[i] + fg[i];
                        }
                    }
                }
            }
            delete [] fneq;    fneq = nullptr;
            delete [] feq;    feq = nullptr;
            delete [] fg;    fg = nullptr;
        }
        else
        {
            for (int x = xmin; x < xmax; x++)
            {
                for (int y = ymin; y < ymax; y++)
                {
                    for (int z = zmin; z < zmax; z++)
                    {
                        for (int i = 0; i < direction_size; i++)
                        {
                            double feq = calculate_feq(x, y, z, i);
                            double fg = calculate_fg(x, y, z, i, feq);
                            pdf2[index(x, y, z, i)] = omtauinv * pdf[index(x, y, z, i)] + tauinv * feq + fg;
                        }
                    }
                }
            }
        }
    }
    else {
        double *fneq = new double[direction_size];
        double *feq = new double[direction_size];
        double *fg = new double[direction_size];
        double *colli = new double[direction_size];
        if (iLES == 1)
        {
            for (int x = xmin; x < xmax; x++)
            {
                for (int y = ymin; y < ymax; y++)
                {
                    for (int z = zmin; z < zmax; z++)
                    {
                        for (int i = 0; i < direction_size; i++)
                        {
                            feq[i] = calculate_feq(x, y, z, i);
                            fg[i] = calculate_fg(x, y, z, i, feq[i]);
                            fneq[i] = pdf[index(x, y, z, i)] - feq[i];
                        }
                        double p_ij = calculate_pij(fneq);
                        const double tau_t = 0.5 * (sqrt(tau * tau + 2 * sqrt(2) * c_s * c_s * p_ij / (c_s * c_s * c_s * c_s)) - tau);
                        calculate_MSM(tau + tau_t);
                        for (int i = 0; i < direction_size; i++)
                        {
                            colli[i] = 0.0;
                            for (int k = 0; k < direction_size; k++)
                            {
                                colli[i] = colli[i] + MSM[i][k] * fneq[k];
                            }
                            pdf2[index(x, y, z, i)] = pdf[index(x, y, z, i)] + colli[i] + fg[i];
                        }
                    }
                }
            }
        }
        else
        {
            for (int x = xmin; x < xmax; x++)
            {
                for (int y = ymin; y < ymax; y++)
                {
                    for (int z = zmin; z < zmax; z++)
                    {
                        for (int i = 0; i < direction_size; i++)
                        {
                            feq[i] = calculate_feq(x, y, z, i);
                            fg[i] = calculate_fg(x, y, z, i, feq[i]);
                            fneq[i] = pdf[index(x, y, z, i)] - feq[i];
                        }
                        for (int i = 0; i < direction_size; i++)
                        {
                            colli[i] = 0.0;
                            for (int k = 0; k < direction_size; k++)
                            {
                                colli[i] = colli[i] + MSM[i][k] * fneq[k];
                            }
                            pdf2[index(x, y, z, i)] = pdf[index(x, y, z, i)] + colli[i] + fg[i];
                        }
                    }
                }
            }
        }
        delete [] colli;    colli = nullptr;
        delete [] fneq;    fneq = nullptr;
        delete [] feq;    feq = nullptr;
        delete [] fg;    fg = nullptr;
    }
}

void LBM::perform_timestep()
{
    time++;
    compute_density_momentum(0, NX, 0, NY, 0, NZ);
    collision(0, NX, 0, NY, 0, NZ);

    mpi_send();
    stream(0, NX, 0, NY, 0, NZ);
    wall_bounce_back();    //! 处理平直壁面无滑移（simple bounce-back） 

    if (iGEO == 1)
    {
        bouzidi();
    }

    velocity_boundary_condition();    //! 速度边界条件
    pressure_boundary_condition();    //! 压力边界条件
}

void LBM::output(int i,int mpi_rank)
{
    if ((i + 1) % FrameRate == 0)
    {
        if (mpi_rank == 0)
        {
            double percentage = (double)(i + 1) / (double)(ttstep) * 100.0;
            cout << "Saving data - " << (i + 1) << "/" << ttstep << " (" << percentage << "%)" << endl;
        }
        mpi_total(mpi_rank);
        if (mpi_rank ==x_np * y_np * z_np - 1)
        {
            write_binary("output/" + std::to_string(i + 1) + ".plt");
            if ((i + 1) % (10 * FrameRate) == 0)
            {
                output_all("continue.plt");
            }
        }
    }
}

void LBM::check_converge()
{
    double Err = 0;
    double R = 0;
    for (int z = 0; z < total_NZ; z++)
    {
        for (int y = 0; y < total_NY; y++)
        {
            for (int x = 0; x < total_NX; x++)
            {
                Err = Err + fabs(velTotal[total_index(x, y, z)].x - velTotalConv[total_index(x, y, z)].x)
                    + fabs(velTotal[total_index(x, y, z)].y - velTotalConv[total_index(x, y, z)].y)
                    + fabs(velTotal[total_index(x, y, z)].z - velTotalConv[total_index(x, y, z)].z);
                R = R + fabs(velTotal[total_index(x, y, z)].x)
                    + fabs(velTotal[total_index(x, y, z)].y)
                    + fabs(velTotal[total_index(x, y, z)].z);
                velTotalConv[total_index(x, y, z)].x = velTotal[total_index(x, y, z)].x;
                velTotalConv[total_index(x, y, z)].y = velTotal[total_index(x, y, z)].y;
                velTotalConv[total_index(x, y, z)].z = velTotal[total_index(x, y, z)].z;
            }
        }
    }
    Err = Err / R;
    if (Err <= conv_v)
    {
        output_all("output/convfinal.plt");
        cout << "CONV Condition Finish" << endl;
        exit (1);
    }
}

int LBM::get_time()
{
    return time;
}

int LBM::get_NX()
{
    return NX;
}

int LBM::get_NY()
{
    return NY;
}

int LBM::get_NZ()
{
    return NZ;
}

int LBM::get_MeshType()
{
    return MeshType;
}

int LBM::get_border_node_num()
{
    return border_node_num;
}

int LBM::get_total_steps()
{
    return ttstep;
}

int LBM::get_FrameRate()
{
    return FrameRate;
}

int LBM::getdirs()
{
    return send_direction_size;
}

int LBM::getalldirs()
{
    return direction_size;
}

vector3<int>* LBM::get_ei()
{
    return ei;
}
