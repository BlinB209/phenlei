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
//! @file      LBMSolverOMP.hpp
//! @brief     Lattice-Boltzmann solver in OpenMP parallel environment.
//! @author    Huang Haibo, Zhang Jiaqian.

#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include "vector3.hpp"
#ifndef LBM_HILE
#define LBM_HILE

#define PARALLEL  1    //! 程序是否并行编译的开关

int RunLBMSolverOMP();
struct Boundary_Condition
{   //!types :速度边界条件velocity;   压力边界条件pressure； 周期边界periodic
    char xminFace[20];
    char xmaxFace[20];
    char yminFace[20];
    char ymaxFace[20];
    char zminFace[20];
    char zmaxFace[20];
    vector3<double>  vel_specified[6];    //! 若是速度边界：6个面（边界）上的速度矢量
    double rho_specified[6];    //! 若是压力边界:  6个边界上的密度（压力）
};

class LBM
{

public:
//    LBM(int grid_size, std::string velocity_set, double c_s, double tau );
    LBM(int nx, int ny, int nz, std::string velocity_set, double c_s, double tau, vector3<double> gf, Boundary_Condition BC, vector3<double> vel0);
    ~LBM();
    void set_velocity(int x_field, int y_field, int z_field, double u_x, double u_y, double u_z);    //! Set velocity at position in velocity field.
    void set_density(int x_field, int y_field, int z_field, double density);    //! Set density at position in density field.
    double calculate_feq(int i, int j, int k, int w);
    void calculate_feq2(int i, int j, int k, double *feq);
    double calculate_feq(int i, int j, int k, int w, double u_le_x);
    //! 重载函数，该函数是为移动边界服务，f_i= f_i^eq（rho, (u+u_le_x), v, w）， 其中u_le_x为板滑移速度
    void calculate_fg(int i, int j, int k, double *feq, double *fg);    //! 计算体积力
    void calculate_fg_IB(int i, int j, int k, double* feq, double* fg);    //! 计算IB传递的体积力
    double calculate_fneq(int i, int j, int k, int w, double feq);    //!计算非平衡分布函数
    double calculate_pij(double *fneq);    //! 计算动量

    //! BC
    void bouzidi(int num_node2, int *node2,double *delta, vector3<double> *force);
    void wall_bounce_back();    //! 处理平直壁面无滑移（simple bounce-back）
    void velocity_boundary_condition();    //! 处理速度边界条件
    void pressure_boundary_condition();    //! 处理压力边界条件
    void lookup_reverse();    //! 求 i 方向的反方向
    //! LBM iteration
    void compute_density_momentum(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);
    void stream();    //! Stream the current equilibrium distribution to the next distribution.
    void collision(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);    //! Perform the collision step. Assumes delta t / tau = 1.
    void collision_BF(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);    //! include body force
    void collision_IB(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);    //! include body forc - IB
    void collision_LES(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);    //! including LES  SGS model
    void collision_zip(int MRT, int BF, int LES, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);    //! all collisons in a whole
    void perform_timestep(int MRT, int BF, int LES);    //! Delta t = 1 lattice unit.
    //! output
    void output_history(std::string filename, bool header, int i, vector3<double> *force);

    void output_velocity_profile_at_a_section();    //! 输出某截面上的速度剖面
    void output_indices_file();
    void output_tecplot(std::string filename, bool header = true);
    void output_all_fi(std::string filename);
    int  get_time();
    void write_binary(std::string filename);

    //! MPI
    void mpi_id(int mpi_rank);
    void mpi_send();
    //! MRT
    void iniMRT();
    void iniMRT2();
    void calculate_MSM(double tau);
    void D3Q19calculate_meq(int x, int y, int z, double *meq);
    //void D3Q19MRTcollision();
    void D2Q9calculate_meq(int x, int y, double *meq);
    //void D2Q9MRTcollision();
    void MRTcollision(int iLES, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);
    void MRTcollision_IB(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);
    //! C2F, F2C
    int getdirs();
    int getalldirs();
    void exchange(int id,int border_node_num, int *border, int *border_dir, double *f_border_all_dir, double *f_border_dir, int *Low, int *NN);
    void init_border(int *Low, int Fdimx, int Fdimy, int Fdimz, int &num, int *border, int *border_dir);
    void copy_data(int id, vector3<double> *vel, double *rho, int *Low, int *NN);    //! copy macrovariable data from fine to coarse
    vector3<int>* get_ei();

    //! IB
    //! 得到private速度宏观量vel值
    void getVel(vector3<double> **uvw_temp, double **rho_temp)
    {
        (*uvw_temp) = vel;
        (*rho_temp) = rho;
    }
    //! 将IB体积力加上private重力值
    void update_f(vector3<double> **f, int nx, int ny, int nz)
    {
        for (int z = 0; z < nz; z++)
        {
            for (int y = 0; y < ny; y++)
            {
                for (int x = 0; x < nx; x++)
                {
                    (*f)[index(x, y, z)].x += gf.x;
                    (*f)[index(x, y, z)].y += gf.y;
                    (*f)[index(x, y, z)].z += gf.z;
                }
            }
        }
    }
    //! 设置private是否使用IB以及传递IB体积力
    void setIBFSI(int ib)
    {
        IBFSI = ib;
    }
    void setIB_f(vector3<double> *f)
    {
        f_IB = f;
    }

private:
    int MRT, CONTI, LES, GEO, CONV;
    //! MRT
    double **M_MRT, **M_MRTI, **MSM;
    double *Inv_MS, *S_diag2, *Mtemp2, *F_diag2;
    int time = 0;
    int NX;
    int NY;
    int NZ;
    double c_s;
    double nu;
    double tau;
    vector3<double> gf;
    int IB_Multi_Geo;    //! Mutli_Geo numbers
    int IBFSI;    //! Judge whether to use IB
    vector3<double> *f_IB;    //! IB force
    int *obst;    //! obst 0: fluid nodes, 1: solid nodes. 2: solid nodes which is most closest to fluid 
    int *lab;    //! label the node type to detect the directions of streaming in the CF interface
    double *rho;
    int *border, *border_dir;
    double *f_border_dir;
    double *f_border_all_dir;
    int *node2;
    double *delta;
    vector3<double> *vel;    //! bc_specified_vel;
    vector3<double> *velF;
    double *rhoF;
    double *pdf2;
    double *pdf;
    vector3<double> *force;
    int x_id;
    int y_id;
    int z_id;
    int x_np=2;
    int y_np=1;
    int z_np=1;
    int n_myid[6];
    inline int index(int x, int y, int z) const
    {
        return (z * NX * NY) + (y * NX) + x;
    }
    inline int index(int x, int y, int z, int w) const
    {
        return (x + y * NX + z * NX * NY + w * NX * NY * NZ);
    }
    void output_array(double *array);    //! 目前没用上
    void set_velocity_set(std::string velocity_set);    //! Used to internally generate velocity_set.
    void initialise(vector3<double> vel0);
    //void initialise_re();
    //Lattice ei using D3DQ15. assumed speed of sound c_s = 1/sqrt(3).
    int direction_size;
    int send_direction_size;
    vector3<int> *ei;
    double *weights;
    int *reverse_indexes;
    std::string boundary_condition;
    //double gamma_dot;
    std::string velocity_set;
    Boundary_Condition BC;
};

struct Params_Pakage
{
    char in[50]; 
    char out[50];
    int MRT, CONTI, LES, GEO, MB, CONV;
    int IBFSI, IB_Multi_Geo;
    int LX, LY, LZ;
    int LX2, LY2, LZ2;
    int LowX, LowY, LowZ;
    int FrameRate;
    int ttstep;
    double density0, conv_v;
    vector3<double> vel0;
    double tau, Re_LBM;
    vector3<double> gf;
    char velocity_set[20];
    char conti_file[20],geo_file[20];
    char bc[20];
    Boundary_Condition BC;
};

void read_params( Params_Pakage *p);

#endif
