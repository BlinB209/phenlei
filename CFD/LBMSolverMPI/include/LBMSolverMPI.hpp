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
//! @file      LBMSolverMPI.hpp
//! @brief     Lattice-Boltzmann solver in MPI parallel environment.
//! @author    Huang Haibo, Zhang Jiaqian.

#pragma once
#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <mpi.h>
#include <math.h>
#include "vector3.hpp"

void RunLBMSolverMPI();

//#define PARALLEL  1   //程序是否并行编译的开关
using namespace std;

struct Boundary_Condition
{   //! types :速度边界条件velocity;   压力边界条件pressure； 周期边界periodic
    char xminFace[20];
    char xmaxFace[20];
    char yminFace[20];
    char ymaxFace[20];
    char zminFace[20];
    char zmaxFace[20];
    vector3<double> vel_specified[6];    //! 若是速度边界：6个面（边界）上的速度矢量
    double rho_specified[6];    //! 若是压力边界:  6个边界上的密度（压力）
};

struct Params_Pakage
{
    char in[50];
    char out[50];
    char velocity_set[20];
    char bc[20];
};

class LBM
{

private:
    int FINE, iMRT, iLES, iContinue, iGEO, CONV;    //! 程序开关
    double conv_v;
    string conti_file, geo_file;
    int time = 0;
    int total_NX, total_NY, total_NZ;
    int NX, NY, NZ, x_np, y_np, z_np, x_id, y_id, z_id;

    //! MRT
    double** M_MRT, ** M_MRTI, ** MSM;

    //! fine
    int fine_NX, fine_NY, fine_NZ;
    int LefLowx, LefLowy, LefLowz;
    int Fdimx, Fdimy, Fdimz;
    int MeshType;
    int Xori, Yori, Zori;
    int finelocal[6];
    int border_node_num;

    int FrameRate, ttstep;
    double density0;
    double tau;
    double gx, gy, gz;
    const double c_s = 1.0 / sqrt(3);

    string boundary_condition;
    string velocity_set;
    Boundary_Condition BC;
    MPI_Comm *comm0;
    MPI_Comm *comm1;
    //
    double *rho;
    vector3<double> vel0;    //! ini_vel
    vector3<double> *vel;    //! bc_specified_vel;
    double *pdf2;
    double *pdf;
    int *obst;
    double *rhoTotal{};
    vector3<double> *velTotal{};
    vector3<double> *velTotalConv{};
    int *obstTotal{};
    int *n_myid;
    int direction_size;
    int send_direction_size;
    vector3<int> *ei;
    double *weights;
    int *reverse_indexes;

    vector3<double> *force;
    int *node2;
    double *delta;
    int num_node2;

    inline int index(int x, int y, int z) const
    {
        return ((z + 2) * (NX + 4) * (NY + 4) + (y + 2) * (NX + 4) + (x + 2));
    }
    inline int index(int x, int y, int z, int w) const
    {
        return ((x + 2) + (y + 2) * (NX + 4) + (z + 2) * (NX + 4) * (NY + 4) + w * (NX + 4) * (NY + 4) * (NZ + 4));
    }
    inline int index(int x, int y, int z, int NX, int NY)
    {
        return (z + 2) * (NX + 4) * (NY + 4) + (y + 2) * (NX + 4) + (x + 2);
    }
    inline int total_index(int x, int y, int z) const
    {
        return (z * total_NX * total_NY) + (y * total_NX) + x;
    }

public:
    LBM();
    ~LBM();
    int  initialise(Params_Pakage *p, int mpi_rank);    //! 总初始化计算域过程
    void perform_timestep();    //! 总计算步
    void output(int t, int mpi_rank);    //! 总输出步

    int get_time();
    int get_NX();
    int get_NY();
    int get_NZ();
    int get_MeshType();
    int get_border_node_num();
    int get_total_steps();
    int get_FrameRate();
    int getdirs();
    int getalldirs();
    vector3<int> *get_ei();

private:
    //! MPI
    void mpi_ini(int mpi_rank);    //! 并行计算初始化
    void mpi_send();    //! 传递边界pdf2
    void mpi_total(int mpi_rank);    //! 汇总各CPU数据
    //! 初始参数
    void ReadParameter(Params_Pakage *p, int mpi_rank);    //! 读入参数
    void iniFlow();    //! 流场区域初始化
    void iniobst();    //! 固体区域初始化
    //! output
    void write_binary(string filename);
    void read_all(std::string filename);    //! 读入流场数据，文件应该与输出匹配
    void output_all(std::string filename);    //! 输出流场数据
    void output_all(std::string filename, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);
    //! LBM iteration
    void compute_density_momentum(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);    //! 计算宏观量
    void stream(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);    //! Stream the current equilibrium distribution to the next distribution.
    void collision(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);    //! Perform the collision step. Assumes delta t / tau = 1.
    //! BC
    void bouzidi();
    void wall_bounce_back();    //! 处理平直壁面无滑移（simple bounce-back）
    void velocity_boundary_condition();    //! 处理速度边界条件
    void pressure_boundary_condition();    //! 处理压力边界条件
    void set_velocity_set();    //! Used to internally generate velocity_set.

    void set_velocity(int x_field, int y_field, int z_field, double u_x, double u_y, double u_z);    //! Set velocity at position in velocity field.
    void set_density(int x_field, int y_field, int z_field, double density);    //! Set density at position in density field.
    double calculate_feq(int i, int j, int k, int w);
    double calculate_feq(int i, int j, int k, int w, double u_le_x);
    //! 重载函数，该函数是为移动边界服务，f_i= f_i^eq（rho, (u+u_le_x), v, w）， 其中u_le_x为板滑移速度
    double calculate_fg(int i, int j, int k, int w, double feq);    //! 计算体积力
    double calculate_fneq(int i, int j, int k, int w, double feq);    //! 计算非平衡分布函数
    double calculate_pij(double *fneq);    //! 计算动量
    void lookup_reverse();    //! 求 i 方向的反方向
    //! MRT
    void iniMRT();    //! MRT计算初始化
    void calculate_MSM(double tau);    //! MRT矩阵计算
    void read_Geo(std::string filename, int mpi_rank);    //! 读入复杂外形
    void getNumberFromString(string s, int *a, int n);
    bool IntersectTriangle(vector3<double> &orig, vector3<double> &dir,    //! 判断固体点
        vector3<double> &v0, vector3<double> &v1, vector3<double> &v2,
        double *t, double *u, double *v);
    void check_converge();    //! 验证是否收敛
};

void RUNPHLBM(int mpi_rank);
void dumpstring(char *str, FILE *file);