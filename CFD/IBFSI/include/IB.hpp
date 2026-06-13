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
//! @file      IB.hpp
//! @brief     The function of the immersed boundary method.
//! @author    Peng Zerui, Xu Pingyu. 

#include <iostream>
#include <string>
#include <fstream>
#include "..\..\LBMSolverOMP\include\vector3.hpp"
#include "vector.hpp"
#ifndef IB_HILE
#define IB_HILE


class IBFSI {

public:
    const double pi = 3.1415926535;

    IBFSI(int nx, int ny, int nz, double dT, vector3<double> vel0, double density0, double elex, double eley, double elez, string input_IB, string meshFlieName, string outFileName);
    void pIB(vector3<double>*, double*, int);
    void interpolate_to_lag_nodes();
    void get_lag_elements(vector<double>* , vector<double>* , vector<double>* , vector<double>* , vector<double>* , vector<double>*, vector<double>*);
    void interpolate_to_lag_elements(vector<double>*, vector<double>*, vector<double>, vector<double>);
    void calculate_lagrangian_element_forces(vector<double>*, vector<double>*, double, double, vector<double>, vector<double>,
        vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>);
    void calculate_body_forces_on_eular(vector<double>* , vector<double>* , vector<double>, vector<double>, vector<double>, vector<double>, vector<double>);
    void calculate_node_forces(vector<double>, vector<double>);
    void prescribed_motion(int, int);
    void filament_motion(int);
    void filament_solver(double dtfilm,
        vector<double> velfix,
        vector<double> angfix,
        vector<vector<double>> &xfo,
        vector<vector<double>> &xfn,
        vector<vector<double>> &ufn,
        vector<vector<double>> exf);
    void init_VIV();
    void VIV_motion();
    void solve_VIV(double deltat,
        vector<double> &X,
        vector<double> &Fn,
        vector<double> &Fn_1,
        vector<double> &Fn_2,
        vector<double> &Fn_3,
        vector<double> &Qn,
        vector<double> &Qn_1,
        vector<double> &Qn_2,
        vector<double> &Qn_3,
        vector<vector<double>> &A,
        vector<vector<double>> &B);
    void write_everything(int, string);

    void initialise(int nx, int ny, int nz, double elex, double eley, double elez, vector3<double> inflow_Velocity, string outFileName);
    void setCharacterQuantity();
    void initialBasicParams();

    void getInputFromFile(string file_str);
    void getLagFromMesh(string meshfilename);
    void output_tecplot(string filename, bool header = true);

    void IBgetNumberFromString(string, int*);

    void velocity_Conversion_vector3_To_vector(vector3<double>*, double*);
    void Conversion_vector_To_vector3(vector3<double>**);
    
    void setNxyz(int nx, int ny, int nz)
    {
        this->Nx = nx;
        this->Ny = ny;
        this->Nz = nz;
    }
    void setElExyz(double elesizex, double elesizey, double elesizez)
    {
        this->eleSizeX = elesizex;
        this->eleSizeY = elesizey;
        this->eleSizeZ = elesizez;
    }
    void setdenIn(double density0)
    {
        this->denIn = density0;
    }

    double getTau()
    {
        cout << "tau = " << 0.5 + 3.0 * Lref * Uref / (Re * dt) << endl;
        return 0.5 + 3.0 * Lref * Uref / (Re * dt);
    }
    int getNumomp()
    {
        return this->numsOMP;
    }

private:
    int Nx, Ny, Nz;
    //! ----------拉格朗日间距需改成每点对应数组-----------
    double dt;    //! dt计算步长；ds_Lag拉格朗日间距
    double eleSizeX, eleSizeY, eleSizeZ;    //! 流体域各方向

    vector<double> lagx, lagy, lagz;    //! 不含质量的每时间步实时拉格朗日坐标

    //! *********** ***********
    int multiGeo;    //! 结构数量
    int numsOMP;
    int iMotion, subMotionType1;
    int subMotionType2 = 2;
    int Nnode, Nelem;
    vector<double> Xnode, Ynode, Znode;
    vector<double> X0node, Y0node, Z0node;
    vector<double> X00node, Y00node, Z00node;
    vector<double> XIBnode, YIBnode, ZIBnode;
    vector<double> UIBnode, VIBnode, WIBnode, RhoIBnode;
    vector<int> Ielem, Jelem, Kelem;    //! 三个单元编号
    vector<double> xGrid, yGrid, zGrid;
    double ds0, dx, dy;
    vector3<double> vel0;
    double denIn;
    double Re;
    double Lref, Uref, Tref, Fref, KBref, KSref;    //! 特征量
    vector<double> Unode, Vnode, Wnode;
    vector<double> Axnode, Aynode, Aznode;
    vector<double> Fxnode, Fynode, Fznode;
    vector<double> extfulx, extfuly;
    vector<vector<double>> u, v, w, rho;    //! 流体域的速度

    vector<vector<vector<double>>> fx, fy, fz;
    double Cd_ib = 0.0;
    double Cl_ib = 0.0;
    string CdCl_name;

    //! flapping rigid plate
    int IsRotateOrnot;
    double iniAngle, iniPhase;
    double A_rotate, T_rotate;

    //! filament
    double kb_film, ks_film;
    double Ampl_heaving, Ampl_pitching;
    double Period_heaving, Period_pitching;
    double angfix_film1, angfix_film2;
    double velfix_film1, velfix_film2;
    double accfix_film1, accfix_film2;

    //! VIV
    double VIV_MassRotio;
    double VIV_NonDimNatFreq_y;    //! 两个方向无量纲固有频率
    double VIV_NonDimNatFreq_Ang;
    double VIV_DampRatio_y;    //! 阻尼比
    double VIV_DampRatio_Ang;
    double VIVMassCenter_X;    //! 质心位置
    double VIVMassCenter_Y;

    double VIV_Mass;    //! 质量
    double VIV_Inertia_Moment;

    vector<double> VIV_X;    //! 位移、弧度、加速度、角加速度
    vector<double> VIV_Fn, VIV_Fn_1, VIV_Fn_2, VIV_Fn_3;
    vector<double> VIV_Qn, VIV_Qn_1, VIV_Qn_2, VIV_Qn_3;
    vector<vector<double>> VIV_A, VIV_B;

    double VIV_K_y;    //!有量纲固有频率
    double VIV_K_Ang;
    double VIV_Damp_y;    //! 有量纲阻尼固有频率
    double VIV_Damp_Ang;

    double VIVForce_Fy;    //! 涡激气动力
    double VIVForceMoment_M;    //! 涡激气动力矩

    //! 基础参数
    double solvertime = 0.0;
    double solid_zero_positionX_percent, solid_zero_positionY_percent;
    double gravityX, gravityY, gravityZ;
    int supp;    //! 分布函数影响区域，必须为偶数
};

#endif