
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <streambuf>
//!    Now Linux only.
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <omp.h>
//#include <mpi.h>
#include "LBMSolverOMP.hpp"
#include "vector3.hpp"
#include "IB.hpp"

void read_Geo(int direction_size, vector3<int> *ei, int *obst, int *NN2, int &num_node2, int *node2, double *delta);

inline bool file_exists (const std::string& name)
{
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

inline int index(int x, int y, int z, int NX, int NY)
{
    return (z * NX * NY) + (y * NX) + x;
}

inline void mutilple_geo_add(vector3<double>** ff1, vector3<double>* ff2, int nx, int ny, int nz)
{
    int k = (nz == 1) ? 0 : 0;

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            int ind = i + j * nx + k * nx * ny;
            (*ff1)[ind].x = (*ff1)[ind].x + ff2[ind].x;
            (*ff1)[ind].y = (*ff1)[ind].y + ff2[ind].y;
            (*ff1)[ind].z = (*ff1)[ind].z + ff2[ind].z;
        }
    }
}

//int main(int argc, char** argv) {
int RunLBMSolverOMP() 
{
    std::cout.precision(4);    //! 数值8改成4了，要不然位数太多
    clock_t tStart = clock();    //! report timing information

    std::string command;
    command = "mkdir output";
    system(command.c_str());

    //! input filename and output filename
    struct Params_Pakage parameters = { "./bin/input.txt","./params_out.txt" };    //! 类部分初始化

/*    std::cout << "Do you want to clean the previous run? (1 - Yes, 0 - No): ";
    int choice;
    std::cin >> choice;
    if(choice == 1) {
        system("rm -rf output");
        system("mkdir output");
    }
*/

    omp_set_num_threads(1);
    read_params(&parameters);    //! 见read.cpp，从文件./input.txt读入模拟参数

    int NX = parameters.LX;    //! 计算域x方向网格数；
    int NY = parameters.LY;    //! 计算域y方向网格数；
    int NZ = parameters.LZ;    //! 计算域z方向网格数；
    int NX2 = parameters.LX2;    //! FINE计算域x方向网格数；
    int NY2 = parameters.LY2;    //! FINE计算域y方向网格数；
    int NZ2 = parameters.LZ2;    //! FINE计算域z方向网格数；
    int LowX = parameters.LowX;
    int LowY = parameters.LowY;
    int LowZ = parameters.LowZ;
    int Fdimx = (NX2) / 2;
    int Fdimy = (NY2) / 2;
    int Fdimz = (NZ2) / 2;
    if (Fdimz == 0) Fdimz = 1;    //! serve for 2D cases D2Q9

    double dT_lattice = 1.0;    //! 影响范围；时间步长
    double eleSizeX = 1.0, eleSizeY = 1.0, eleSizeZ = 1.0;    //! 三个方向的网格间距
    IBFSI *ibcontact = nullptr, *ibcontact2 = nullptr;

    //! 当不使用IBFSI模块时，下面的参数存在但不使用
    if (parameters.IBFSI == 1)
    {
        string input_IB = "input_IB.dat";
        string mesh_IBM = "mesh2d_IBM.dat";
        string CdCl_name = "output/CdCl_A.plt";
        ibcontact = new IBFSI(NX, NY, NZ, dT_lattice, parameters.BC.vel_specified[0], parameters.density0, eleSizeX, eleSizeY, eleSizeZ, input_IB, mesh_IBM, CdCl_name);
        if (parameters.IB_Multi_Geo == 2)
        {
            input_IB = "input_IB2.dat";
            mesh_IBM = "mesh2d_IBM2.dat";
            CdCl_name = "output/CdCl_B.plt";
        }
        ibcontact2 = new IBFSI(NX, NY, NZ, dT_lattice, parameters.BC.vel_specified[0], parameters.density0, eleSizeX, eleSizeY, eleSizeZ, input_IB, mesh_IBM, CdCl_name);
        parameters.tau = ibcontact->getTau();

        //! 目前IBM暂时没有添加多块网格MB参数与LES功能
        if (parameters.MB == 0 && parameters.LES == 0)
        {
            omp_set_num_threads(ibcontact->getNumomp());
            parameters.tau = ibcontact->getTau();
        }
        else
        {
            cout << "Currently, IBM is not adding multi-block grid MB parameters and LES functionality. Exit." << endl;
            exit(200);
        }
    }

    std::cout << "Grid size: " << NX << "x" << NY << "x" << NZ << '\n';
    if (parameters.MB == 1)
    {
        std::cout << "Fine Grid size: " << NX2 << "x" << NY2 << "x" << NZ2 << '\n';
    }
    int FrameRate = parameters.FrameRate;    //! FrameRate is the number of timesteps per frame.
    std::cout << "number of timesteps per frame: " << FrameRate << '\n';

    double c_s = 1.0 / sqrt(3);
    std::cout << "c_s (Speed of sound): " << c_s << '\n';

    double tau = parameters.tau;    //! initial tau value;
    double tau2 = (tau - 0.50)*2.0 + 0.50;
    std::cout << "Tau: " << tau << '\n';
    if (parameters.MB == 1)
    {
        std::cout << "TauFine: " << tau2 << '\n';
    }

    vector3<double> gf = parameters.gf;    //! initial Gx,Gy,Gz value;
    int BF = 0;    //! no body force
    if (fabs(gf.x) > 1.e-12 || fabs(gf.y) > 1.e-12 || fabs(gf.y) > 1.e-12)  BF = 1;    //! Body force 标志
    std::cout << "Body force   " << BF << '\n';

    if (tau < 0.5)
    {
        //! Section 3.5.5.1 with delta t=1.
        std::cout << "Error: Tau must be greater than 0.5 for numerical stability." << '\n';
        return -1;
    }

    std::string velocity_set = parameters.velocity_set;    //! 读入速度模型;

    //std::string boundary_conditions = parameters.bc;    //! 边界条件类型;
    if (velocity_set != "D3Q15" && velocity_set != "D3Q19" && velocity_set != "D3Q27" && velocity_set != "D2Q9")
    {
        std::cout << "Error: Please specify a valid velocity set such as D3Q15,D3Q27 or D2Q9." << '\n';
        return -1;
    }
    std::cout << "Velocity set: " << velocity_set << '\n';
    //if(boundary_conditions != "xxx") {
    //    std::cout << "Errors: boundary_conditions in *****)";
    //    return -1;
    //}可以加入一个检查输入文件一致性的这样的一个函数。比如取名叫check consistency。如果有什么不一致的地方，就提示出来。比如说本来是压力边界条件，然后输入文件中 Pressure后面又跟了速度矢量。

    if (NZ != 1 && velocity_set == "D2Q9")
    {
        std::cout << "Warning: Please set NZ=1 for D2Q9.";
        return -1;
    }
    int total_steps = parameters.ttstep;    //! total steps to terminate the simulation;

    LBM *solver = new LBM(NX, NY, NZ, velocity_set, c_s, tau, gf, parameters.BC, parameters.vel0);
    if (parameters.MB == 0) { NX2 = 1; NY2 = 1; NZ2 = 1; }
    LBM *solver2 = new LBM(NX2, NY2, NZ2, velocity_set, c_s, tau2, gf, parameters.BC, parameters.vel0);

    int Low[3] = { LowX,LowY,LowZ };
    int NN[3] = { NX,NY,NZ };
    //solver->mpi_id(mpi_rank);
    //serve for data exchange between solver (Coarse mesh) and solver2 (Fine mesh)
    int *border, *border_dir;
    double *f_border_dir;
    double *f_border_all_dir;
    int dirs = solver->getdirs();    //! get dirs
    int direction_size = solver->getalldirs();    //! get direction size

    //NX2 = 50; NY2 = 50; NZ2 = 50;//remove it!
    int NN2[3] = { NX2,NY2,NZ2 };    //! serve for "read_Geo"  only!
    int NH[3] = { NX2,NY2,NZ2 };    //! serve for "read_Geo" in single block
    vector3<int> *ei;
    ei = solver2->get_ei();
    int *obst;
    if (parameters.MB == 0) NH[0] = NX; NH[1] = NY; NH[2] = NZ;
    obst = new int[NH[0] * NH[1] * NH[2]]();    //! ()初始化为0
    int num_node2 = 10000;
    int *node2;
    double *delta;
    node2 = new int[num_node2]();
    delta = new double[num_node2 * direction_size]();    //! 初始化为0
    if (parameters.GEO == 1)
    {   //! && velocity_set != "D2Q9"
        if(parameters.MB==0) read_Geo(direction_size, ei, obst, NN, num_node2, node2, delta);
        if (parameters.MB == 1) read_Geo(direction_size, ei, obst, NN2, num_node2, node2, delta);
    }

    //if(parameters.GEO==1 && velocity_set == "D2Q9")    read_Geo(2, direction_size, ei, obst, NN2, num_node2, node2, delta);

    vector3<double> *velF;
    double *rhoF;
    velF = new vector3<double>[Fdimx*Fdimy*Fdimz]();
    rhoF = new double[Fdimx*Fdimy*Fdimz]();
    vector3<double> *force;
    force = new vector3<double>[1]();

    int border_node_num = 2 * (Fdimx* Fdimy + Fdimy * Fdimz + Fdimz * Fdimx)
        - 4 * (Fdimx + Fdimy + Fdimz) + 8;

    if (velocity_set == "D2Q9") border_node_num = 2 * (Fdimx + Fdimy) - 4;
    if (velocity_set == "D3Q19") border_node_num = border_node_num - 8;

    border = new int[border_node_num];
    border_dir = new int[border_node_num* dirs]();    //! ()初始化为0；model dependent,for D3Q19, at most outward directions: 5
    //actually dirs= send_direction_size
    f_border_dir = new double[border_node_num* dirs]();
    f_border_all_dir = new double[border_node_num*direction_size]();    //! direction_size
    int b_node_num = 1;

    int NZ2min = 2;    //! general 3D
    int NZ2max = NZ2 - 2;    //! general 3D
    if (velocity_set == "D2Q9")    //! 二维特殊处理下
    {
        NZ2min = 0;
        NZ2max = 1;
    }

    if (parameters.MB == 1) solver->init_border(Low, Fdimx, Fdimy, Fdimz, b_node_num, border, border_dir);

    std::cout << "actual num " << b_node_num << "  theory num " << border_node_num << '\n';

/*    for (int i = 0; i < argc; i++) {
        if (std::string(argv[i]) == "generate_ic") {
            solver->output_tecplot("ic.plt", true);
            std::cout << "Generated ic.plt" << '\n';
            return 0;
        }
    }*/
    if (file_exists("ic.plt"))
    {
        //! 若有初场文件，此处应读入
        int i, j, k, m = 0;
        double density, u_x , u_y , u_z, den_avg = 0.0;    //! 缺省初始化
        char buffer[256];
        std::ifstream infile("ic.plt");
        if (!infile)
        {
            std::cout << "Unable to open input file";
            exit(1);    //! terminate with error

        }
        infile.getline(buffer, 100);
        infile.getline(buffer, 100);
        while (!infile.eof())
        {
            m++;
            infile.getline(buffer, 100);
            sscanf_s(buffer, "%d %d %d %lf %lf %lf %lf", &i,&j,&k,&density,
                &u_x, &u_y,&u_z);
            solver->set_density(i, j, k, density);
            solver->set_velocity(i, j, k, u_x, u_y, u_z);
            den_avg = den_avg + density;
        }
        //solver->initialise_re();
        den_avg = den_avg / (double)m;
        std::cout << "Averged density:" <<den_avg << '\n';
        std::cout << "1Loaded initial conditions (read from file)" << '\n';
    }
    if (file_exists("ic2.plt"))
    {
        //! 若有初场文件，此处应读入
        int i, j, k, m = 0;
        double density, u_x, u_y, u_z;    //! 缺省初始化
        char buffer[256];
        std::ifstream infile("ic2.plt");
        if (!infile)
        {
            std::cout << "Unable to open input file";
            exit(1);    //! terminate with error
        }
        infile.getline(buffer, 100);
        infile.getline(buffer, 100);
        while (!infile.eof())
        {
            m++;
            infile.getline(buffer, 100);
            sscanf_s(buffer, "%d %d %d %lf %lf %lf %lf", &i, &j, &k, &density,
                &u_x, &u_y, &u_z);
            solver2->set_density(i, j, k, density);
            solver2->set_velocity(i, j, k, u_x, u_y, u_z);
         
        }
        //solver->initialise_re();

        std::cout << "2Loaded initial conditions (read from file)" << '\n';
    }
    if(!file_exists("ic2.plt")&& !file_exists("ic.plt"))
    {
        std::cout << "Loaded default initial conditions" << '\n';

        double density = { 1 }, u_x = parameters.vel0.x, u_y = parameters.vel0.y, u_z = parameters.vel0.z;    //! 缺省初始化
        std::cout<< u_x << " " << u_y <<"  "<<u_z << '\n';
        density = parameters.density0;
        for (int i = 0; i < NX; i++)
        {
            for (int j = 0; j < NY; j++)
            {
                for (int k = 0; k < NZ; k++)
                {
                    solver->set_density(i, j, k, density);
                    solver->set_velocity(i, j, k, u_x, u_y, u_z);
                }
            }
        }
        if (parameters.MB == 1)
        {
            for (int i = 0; i < NX2; i++)
            {
                for (int j = 0; j < NY2; j++)
                {
                    for (int k = 0; k < NZ2; k++)
                    {
                        solver2->set_density(i, j, k, density);
                        solver2->set_velocity(i, j, k, u_x, u_y, u_z);
                    }
                }
            }
        }
        //std::cout << "If you wish to use your own initial conditions, please run the program but with command: generate_ic as a argument which will output ic.csv in format of p,u_x,u_y,u_z, assume indexes are incrementing i,j,k for i<NX,j<NY and k<NZ" << '\n';
    }

    double viscosity = c_s * c_s * (tau - 0.5);
    std::cout << "Kinematic shear viscosity: " << viscosity << '\n';

    std::cout << "For velocity set D2Q9,D3Q15,D3Q19 and D3Q27, |u_max|<0.577\n";
    //solver->mpi_id(mpi_rank);
    solver->output_tecplot("0.plt");    //! output
    if (parameters.MB == 1) solver2->output_tecplot("00.plt");    //! output
    //solver->output_indicesile();    //! 输出output/indices.csv，这文件用处不大，不如输出ASCII tecplot文件

    int scale = 1;
    int runs = total_steps * scale * scale * scale;    //! 目前runs = total_steps

    if (parameters.MRT == 1)
    {
        solver->iniMRT();
        if(parameters.MB==1) solver2->iniMRT();
    }

    //! 当不使用IBFSI模块时，下面的参数存在但不使用
    vector3<double> *f_temp_ini = nullptr, *ff = nullptr, *ff2 = nullptr;
    vector3<double> *uvw_temp = nullptr; double *rho_temp = nullptr;
    if (parameters.IBFSI == 1)
    {
        f_temp_ini = new vector3<double>[NX * NY * NZ];
        for (int i = 0; i < NX * NY * NZ; i++)
        {
            f_temp_ini[i].x = 0;
            f_temp_ini[i].y = 0;
            f_temp_ini[i].z = 0;
        }
        ff = new vector3<double>[NX * NY * NZ]();
        ff2 = new vector3<double>[NX * NY * NZ]();
    }

    if (parameters.IBFSI == 1) 
    {
        solver->update_f(&f_temp_ini, NX, NY, NZ);    //! 1.将体积力加上LBM重力
        solver->setIB_f(f_temp_ini);    //! 2.加入LBM变量中
        solver->setIBFSI(parameters.IBFSI);
        for (int i = 0; i < runs; i = i + 1) 
        {
            solver->stream();
            solver->wall_bounce_back();
            solver->velocity_boundary_condition();
            solver->pressure_boundary_condition();

            solver->compute_density_momentum(0, NX, 0, NY, 0, NZ);    //! IB获取宏观量参数

            solver->getVel(&uvw_temp, &rho_temp);

            ibcontact->pIB(uvw_temp, rho_temp, i);    //! 计算体积力
            if (parameters.IB_Multi_Geo == 2) ibcontact2->pIB(uvw_temp, rho_temp, i);

            ibcontact->Conversion_vector_To_vector3(&ff);    //! 输出体积力
            if (parameters.IB_Multi_Geo == 2)
            {
                ibcontact2->Conversion_vector_To_vector3(&ff2);
                mutilple_geo_add(&ff, ff2, NX, NY, NZ);
            }

            solver->update_f(&ff, NX, NY, NZ);    //! 1.将体积力加上LBM重力
            solver->setIB_f(ff);    //! 2.加入LBM变量中

            solver->collision_zip(parameters.MRT, BF, parameters.LES, 0, NX, 0, NY, 0, NZ);

            if ((i + 1) % 10 == 0)
            {
                if (parameters.GEO == 1) solver2->output_history("output/force.plt", true, i, force);
            }

            if ((i + 1) % FrameRate == 0)
            {
#pragma warning(disable : 4996)    //! 获取时间，消除警告
                time_t now = time(0);
                char *dt = ctime(&now);
                std::cout << "present time：" << dt << std::endl;

                std::cout << "Executed in " << (clock() - tStart) << " clock cycles\n";
                tStart = clock();

                if (parameters.MB == 1) solver2->copy_data(2, velF, rhoF, Low, NN2);
                if (parameters.MB == 1) solver->copy_data(1, velF, rhoF, Low, NN2);
                double percentage = double(i + 1) / double(runs) * 100.0;

                std::cout << "Saving data - " << (i + 1) << "/" << runs << " (" << percentage << "%)" << '\n';

                //solver->output_tecplot("output/" +std::to_string(mpi_rank) + std::to_string(i + 1) + ".plt");    //! 输出tecplot
                solver->output_tecplot("output/solver" + std::to_string(i + 1) + "a.plt");    //! 输出tecplot
                ibcontact->output_tecplot("output/Lag_A" + std::to_string(i + 1) + ".plt");
                if (parameters.IB_Multi_Geo == 2) ibcontact2->output_tecplot("output/Lag_B" + std::to_string(i + 1) + ".plt");
                //solver->write_binary("output/solver" + std::to_string(i + 1) + "b.plt");
            }
        }
    }
    else
    {
        for (int i = 0; i < runs; i = i + 1)
        {
            if (parameters.MB == 0)
            {
                solver->compute_density_momentum(0, NX, 0, NY, 0, NZ);
                solver->collision_zip(parameters.MRT, BF, parameters.LES, 0, NX, 0, NY, 0, NZ);
                solver->stream();
                if (parameters.GEO == 1) solver->bouzidi(num_node2, node2, delta, force);
                solver->wall_bounce_back();    //! 处理平直壁面无滑移（simple bounce-back）
                solver->velocity_boundary_condition();    //! 速度边界条件
                solver->pressure_boundary_condition();    //! 压力边界条件

                //solver->perform_timestep(parameters.MRT, BF, parameters.LES); //同时执行 求宏观量、collision,和streaming 三个迭代步骤， 参见LBM.cpp
            }
            else
            {
                //solver->time++;
                solver->compute_density_momentum(0, NX, 0, NY, 0, NZ);
                solver->collision_zip(parameters.MRT, BF, parameters.LES, 0, NX, 0, NY, 0, NZ);
                solver->exchange(1, border_node_num, border, border_dir, f_border_all_dir, f_border_dir, Low, NN);    //! id=1, fetch data from coarse mesh

                solver2->compute_density_momentum(2, NX2 - 2, 2, NY2 - 2, NZ2min, NZ2max);
                solver2->collision_zip(parameters.MRT, BF, parameters.LES, 2, NX2 - 2, 2, NY2 - 2, NZ2min, NZ2max);
                solver2->exchange(3, border_node_num, border, border_dir, f_border_all_dir, f_border_dir, Low, NN);    //! id=3, C2F

                solver->stream();
                solver->wall_bounce_back();    //! 处理平直壁面无滑移（simple bounce-back）

                solver2->stream();
                solver2->bouzidi(num_node2, node2, delta, force);
                //solver2->wall_bounce_back();

                for (int i = 0; i < 1; i++)    //! fine mesh, one more iteration
                {
                    solver2->compute_density_momentum(2, NX2 - 2, 2, NY2 - 2, NZ2min, NZ2max);
                    solver2->collision_zip(parameters.MRT, BF, parameters.LES, 2, NX2 - 2, 2, NY2 - 2, NZ2min, NZ2max);
                    solver2->stream();
                    solver2->bouzidi(num_node2, node2, delta, force);
                }

                solver2->exchange(2, border_node_num, border, border_dir, f_border_all_dir, f_border_dir, Low, NN);    //! id=2, fetch data from fine mesh

                solver->exchange(4, border_node_num, border, border_dir, f_border_all_dir, f_border_dir, Low, NN);    //! id=4, F2C

                //solver2->wall_bounce_back();
                solver->velocity_boundary_condition();    //! 速度边界条件
                solver->pressure_boundary_condition();    //! 压力边界条件
                //solver2->velocity_boundary_condition();    //! 速度边界条件
                //solver2->pressure_boundary_condition();    //! 压力边界条件
            }
            if ((i + 1) % 10 == 0)
            {
                if (parameters.GEO == 1) solver2->output_history("output/force.plt", true, i, force);
            }
            //std::cout <<"i="<<i<< std::endl;
            if ((i + 1) % FrameRate == 0)
            {
#pragma warning(disable : 4996)    //! 获取时间，消除警告
                time_t now = time(0);
                char *dt = ctime(&now);
                std::cout << "present time：" << dt << std::endl;

                std::cout << "Executed in " << (clock() - tStart) << " clock cycles\n";
                //double ttt=(double)(clock() - tStart) / CLOCKS_PER_SEC;
                //std::cout << "Executed in " << ttt << " seconds\n";
                tStart = clock();

                if (parameters.MB == 1) solver2->copy_data(2, velF, rhoF, Low, NN2);
                if (parameters.MB == 1) solver->copy_data(1, velF, rhoF, Low, NN2);
                double percentage = (double)(i + 1) / (double)(runs) * 100.0;

                std::cout << "Saving data - " << (i + 1) << "/" << runs << " (" << percentage << "%)" << '\n';

                //solver->output_tecplot("output/" +std::to_string(mpi_rank) + std::to_string(i + 1) + ".plt");    //! 输出tecplot
                //solver->output_tecplot("output/solver" + std::to_string(i + 1) + "a.plt");    //! 输出tecplot
                solver->write_binary("output/solver" + std::to_string(i + 1) + "b.plt");
                //if (parameters.MB == 1) solver2->output_tecplot("output/solver2_" + std::to_string(i + 1) + "a.plt");    //! 输出tecplot
                if (parameters.MB == 1) solver2->write_binary("output/solver2_" + std::to_string(i + 1) + "b.plt");
                //solver->output_velocity_profile_at_a_section();
            }
        }
    }
    std::cout << std::endl;
    delete solver;    solver = nullptr;
    delete solver2;    solver2 = nullptr;
    if (parameters.IBFSI == 1)
    {
        delete ibcontact;
        delete ibcontact2;
        delete [] f_temp_ini;    f_temp_ini = nullptr;
        delete [] ff;    ff = nullptr;
        delete [] ff2;    ff2 = nullptr;
    }
    return 0;
}
