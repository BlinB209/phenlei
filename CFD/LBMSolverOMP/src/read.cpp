#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include "LBMSolverOMP.hpp"
using namespace std;
//******************************************
// R E A D   P A R A M S (simulation parameters)
//
//  - Read the problem parameters from a file.
//

void read_params(Params_Pakage* p)
{
    char buffer[256];
    std::string A;    //! 字符串
    vector3<double> aa;
    double rho_specified;
    ifstream infile(p->in);     //! 输入文件
    ofstream outfile(p->out);   //! 输出文件

    if (!infile)
    {
        cout << "Unable to open input file";
        exit(1); // terminate with error
    }
    if (!outfile)
    {
        cout << "Unable to open output file";
        exit(1); // terminate with error
    }
    //sscanf_s(tokenstring, "%s", s, _countof(s));
    //sscanf_s(tokenstring, "%c", &c, sizeof(char));

    while (!infile.eof())    //! 逐行读入参数， 以文件结束处为标志，若文件最后多出一行空白行或一行未预料内容，则读入可能会有问题
    {
        infile.getline(buffer, 100);    //! 读入标志是否续算，是为1，否为0，续算则继续读入数据输入文件名
        sscanf_s(buffer, "%d", &p->MRT );
        infile.getline(buffer, 120);    //! 读入标志是否续算，是为1，否为0，续算则继续读入数据输入文件名
        sscanf_s(buffer, "%d", &p->CONTI);
        if(p->CONTI==1) sscanf_s(buffer, "%s", &p->conti_file, (unsigned int)sizeof(p->conti_file));
        infile.getline(buffer, 100);    //! 读入标志是否大涡模拟，是为1，否为0
        sscanf_s(buffer, "%d", &p->LES);
        infile.getline(buffer, 120);    //! 读入标志是否需要输入复杂外形物体，是为1，否为0，需要则继续读入图形tecplot输入文件名
        sscanf_s(buffer, "%d", &p->GEO);
        if (p->GEO == 1) sscanf_s(buffer, "%s", &p->geo_file, (unsigned int)sizeof(p->geo_file));
        infile.getline(buffer, 100);    //! 读入标志是否多块网格，是为1，否为0
        sscanf_s(buffer, "%d", &p->MB);
        infile.getline(buffer, 100);    //! 读入标志是否需要收敛判据，是为1，否为0，若需要则继续读入收敛阈值
        sscanf_s(buffer, "%d", &p->CONV);
        if (p->CONV == 1) sscanf_s(buffer, "%lf", &p->conv_v);

        infile.getline(buffer, 100);    //! 读入计算域 X 方向网格数
        sscanf_s(buffer, "%d", &p->LX);
        infile.getline(buffer, 100);    //! 计算域 Y 方向网格数
        sscanf_s(buffer, "%d", &p->LY);
        infile.getline(buffer, 100);    //! 计算域 Z 方向网格数
        sscanf_s(buffer, "%d", &p->LZ);
        infile.getline(buffer, 100);    //! 计算域 X 方向网格数
        sscanf_s(buffer, "%d", &p->LX2);
        infile.getline(buffer, 100);    //! 计算域 Y 方向网格数
        sscanf_s(buffer, "%d", &p->LY2);
        infile.getline(buffer, 100);    //! 计算域 Z 方向网格数
        sscanf_s(buffer, "%d", &p->LZ2);
        infile.getline(buffer, 100);    //
        sscanf_s(buffer, "%d", &p->LowX);
        infile.getline(buffer, 100);    //
        sscanf_s(buffer, "%d", &p->LowY);
        infile.getline(buffer, 100);    //
        sscanf_s(buffer, "%d", &p->LowZ);
        infile.getline(buffer, 100);    //! 每隔多少步输出一次文件

        sscanf_s(buffer, "%d", &p->FrameRate);
        infile.getline(buffer, 100);    //! 总共多少步终止程序
        sscanf_s(buffer, "%d", &p->ttstep);
        infile.getline(buffer, 100);    //! 初始密度
        sscanf_s(buffer, "%lf", &p->density0);
        infile.getline(buffer, 100);    //! 松弛因子，与运动粘性关系:\nu=c_s^2*(\tau-0.5)\delta t
        sscanf_s(buffer, "%lf", &p->tau);
        infile.getline(buffer, 120);    //! 体积力矢量
        sscanf_s(buffer, "(%lf, %lf, %lf)", &p->gf.x, &p->gf.y, &p->gf.z);
        infile.getline(buffer, 100);    //! LBM速度模型
        sscanf_s(buffer, "%s", &p->velocity_set, (unsigned int)sizeof(p->velocity_set));
        //! 字符串提取时能默认以空格分割
        infile.getline(buffer, 100);    //! 边界条件类型
        sscanf_s(buffer, "%s", &p->bc, (unsigned int)sizeof(p->bc));
        infile.getline(buffer, 256);    //! 文件中的说明行，不取信息，跳过！
        infile.getline(buffer, 256);    //! 文件中的说明行，不取信息，跳过！

        infile.getline(buffer, 100);
        sscanf_s(buffer, "%s", &p->BC.xminFace, (unsigned int)sizeof(p->BC.xminFace));
        A = p->BC.xminFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            p->BC.vel_specified[0] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            p->BC.rho_specified[0] = rho_specified;
        }

        infile.getline(buffer, 100);
        sscanf_s(buffer, "%s", &p->BC.xmaxFace, (unsigned int)sizeof(p->BC.xmaxFace));
        A = p->BC.xmaxFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            p->BC.vel_specified[1] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            p->BC.rho_specified[1] = rho_specified;
        }
        
        infile.getline(buffer, 100);
        sscanf_s(buffer, "%s", &p->BC.yminFace, (unsigned int)sizeof(p->BC.yminFace));
        A = p->BC.yminFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            p->BC.vel_specified[2] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            p->BC.rho_specified[2] = rho_specified;
        }

        infile.getline(buffer, 256);
        sscanf_s(buffer, "%s", &p->BC.ymaxFace, (unsigned int)sizeof(p->BC.ymaxFace));
        A = p->BC.ymaxFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            p->BC.vel_specified[3] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            p->BC.rho_specified[3] = rho_specified;
        }

        infile.getline(buffer, 100);
        sscanf_s(buffer, "%s", &p->BC.zminFace, (unsigned int)sizeof(p->BC.zminFace));
        A = p->BC.zminFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            p->BC.vel_specified[4] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            p->BC.rho_specified[4] = rho_specified;
        }

        infile.getline(buffer, 100);
        sscanf_s(buffer, "%s", &p->BC.zmaxFace, (unsigned int)sizeof(p->BC.zmaxFace));
        A = p->BC.zmaxFace;
        if (A == "velocity")
        {
            sscanf_s(buffer, "velocity (%lf, %lf, %lf)", &aa.x, &aa.y, &aa.z);
            p->BC.vel_specified[5] = aa;
        }
        if (A == "pressure")
        {
            sscanf_s(buffer, "pressure %lf", &rho_specified);
            p->BC.rho_specified[5] = rho_specified;
        }

        infile.getline(buffer, 100);    //! 跳一行
        infile.getline(buffer, 100);
        sscanf_s(buffer, "(%lf, %lf, %lf)", &p->vel0.x, &p->vel0.y, &p->vel0.z);

        infile.getline(buffer, 100);    //!IBFSI模块是否使用多个结构
        sscanf_s(buffer, "%d", &p->IB_Multi_Geo);
        infile.getline(buffer, 100);    //! 是否使用IBFSI模块
        sscanf_s(buffer, "%d", &p->IBFSI);

        cout <<"Computational domain size:   "<< p->LX  <<"  "<< p->LY<<"   "<< p->LZ<< endl;
        cout <<"The velocity model used:   "<< p->velocity_set <<'\n'<< endl;

        cout << p->bc << '\n' << endl;

    }

    infile.close();

    //! 读同时输出以上参数到另一文件，以验证读取的正确性
    outfile << "MRT:" << "   " << p->MRT << endl;
    outfile << "CONTI:" << "   " << p->CONTI << endl;
    outfile << "LES:" << "   " << p->LES << endl; 
    outfile << "GEO:" << "   " << p->GEO << endl;
    outfile << "MB:" << "   " << p->MB << endl;
    outfile << "LX:" << "   " << p->LX << endl;
    outfile << "LY:" << "   " << p->LY << endl;
    outfile << "LZ:" << "   " << p->LZ << endl;
    outfile << "LX2:" << "   " << p->LX2 << endl;
    outfile << "LY2:" << "   " << p->LY2<< endl;
    outfile << "LZ2:" << "   " << p->LZ2 << endl;
    outfile << "LowX:" << "   " << p->LowX << endl;
    outfile << "LowY:" << "   " << p->LowY << endl;
    outfile << "LowZ:" << "   " << p->LowZ << endl;
    outfile << "FrameRate:" << "   " << p->FrameRate<< endl;
    outfile << "Total steps:" << " " << p->ttstep << endl;
    outfile << "LBM.initial density:" << " " << p->density0 << endl;
    if (p->IBFSI) 
    {
        outfile << "LBM.tau:" << " Calculated by Re in IBFSI" << endl;
    }
    else 
    {
        outfile << "LBM.tau:" << " " << p->tau << endl;
    }
    outfile << "LBM.tau:" << " " << p->tau << endl;
    outfile << "LBM.gx:" << " " << p->gf.x << endl;
    outfile << "LBM.gy:" << " " << p->gf.y << endl;
    outfile << "LBM.gz:" << " " << p->gf.z << endl;
    outfile << "velocity model:" << " " << p->velocity_set << endl;
    outfile << "Boundary condition type:" << " " << p->bc << endl;
    outfile << "xmin Face BC:  " << " " << p->BC.xminFace << endl;
    outfile << "xmax Face BC:  " << " " << p->BC.xmaxFace << endl;
    outfile << "ymin Face BC:  " << " " << p->BC.yminFace << endl;
    outfile << "ymax Face BC:  " << " " << p->BC.ymaxFace << endl;
    outfile << "zmin Face BC:  " << " " << p->BC.zminFace << endl;
    outfile << "zmax Face BC:  " << " " << p->BC.zmaxFace << endl;
    outfile << "initial vel  " << " " << p->vel0.x << " " << p->vel0.y << " " << p->vel0.z << endl;
    outfile << "multiple geo number:" << " " << p->IB_Multi_Geo << endl;
    outfile << "UseIBPeskin or not:" << " " << p->IBFSI << endl;

    outfile.close();

//    return 0;
}





