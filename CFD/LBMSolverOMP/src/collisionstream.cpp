#include "LBMSolverOMP.hpp"
#include "vector3.cpp"
#include <iomanip>
#include <cmath>
//#include <mpi.h>
/*LBM::LBM(int grid_size, std::string velocity_set, double m_c_s, double m_tau, std::string m_boundary_condition, double m_gamma_dot) : NX(grid_size),NY(grid_size),NZ(grid_size),
        c_s(m_c_s), tau(m_tau), boundary_condition(m_boundary_condition), gamma_dot(m_gamma_dot), velocity_set(velocity_set) {
    set_velocity_set(velocity_set);
    initialise();
}*/
//! 似乎用不到这么多(重)构造函数，上面这个构造函数没用，可以删除
LBM::LBM(int nx, int ny, int nz, std::string velocity_set, double m_c_s, double m_tau, vector3<double> m_gf, Boundary_Condition m_BC, vector3<double> m_vel0) : NX(nx), NY(ny), NZ(nz),
        c_s(m_c_s), tau(m_tau), gf(m_gf), BC(m_BC), velocity_set(velocity_set)
{
    set_velocity_set(velocity_set);
    initialise(m_vel0);
}

LBM::~LBM()
{
    delete [] rho;    rho = nullptr;
    delete [] pdf;    pdf = nullptr;
    delete [] pdf2;    pdf2 = nullptr;
    delete [] obst;    obst = nullptr;
    delete [] vel;    vel = nullptr;
    delete [] weights;    weights = nullptr;
    delete [] reverse_indexes;    reverse_indexes = nullptr;
    delete [] ei;    ei = nullptr;

    delete [] force;    force = nullptr;
    delete [] node2;    node2 = nullptr;
    delete [] delta;    delta = nullptr;

    delete [] border;    border = nullptr;
    delete [] border_dir;    border_dir = nullptr;
    delete [] f_border_dir;    f_border_dir = nullptr;
    delete [] f_border_all_dir;    f_border_all_dir = nullptr;

    delete [] velF;    velF = nullptr;
    delete [] rhoF;    rhoF = nullptr;

    if (MRT == 1)
    {
        delete [] M_MRT;    M_MRT = nullptr;
        delete [] M_MRTI;    M_MRTI = nullptr;
        delete [] MSM;    MSM = nullptr;
        delete [] Mtemp2;    Mtemp2 = nullptr;
        delete [] F_diag2;    F_diag2 = nullptr;
        delete [] S_diag2;    S_diag2 = nullptr;
        delete [] Inv_MS;    Inv_MS = nullptr;
    }
}

void LBM::initialise(vector3<double> vel0)
{
    int NumNodes = NX * NY * NZ;
    int distributions_memo_size = NumNodes * direction_size;
    lab = new int[NumNodes];
    obst= new int[NumNodes]();    //! initial obst
    rho = new double[NumNodes];
    vel = new vector3<double>[NumNodes];
    pdf2 = new double[distributions_memo_size];
    pdf = new double[distributions_memo_size];
    for(int i = 0; i < NX * NY * NZ; i++)
    {
        rho[i] = 1;
        vel[i] = vel0;
    }
    int y, z;
//#pragma omp parallel for private(y,z)
    for(int x = 0; x < NX; x++)
    {
        for(int y = 0; y < NY; y++)
        {
            for(int z = 0; z < NZ; z++)
            {
                for(int i = 0; i < direction_size; i++)
                {
                    pdf2[index(x, y, z, i)] =  weights[i];    //! calculate_feq(x,y,z,i);
                    pdf[index(x, y, z, i)] = weights[i];    //! calculate_feq(x, y, z, i);
                }
            }
        }
    }
}

void LBM::calculate_feq2(int x, int y, int z, double *feq)
{
    //double c_s_squ = 1.0 / 3.0;
    double a=vel[index(x, y, z)].norm_square();
//#pragma omp parallel for 
    for (int i = 0; i < direction_size; i++)
    {
        double dot_product = (double)vel[index(x, y, z)].x * (double)ei[i].x + (double)vel[index(x, y, z)].y * (double)ei[i].y +
            (double)vel[index(x, y, z)].z * (double)ei[i].z;
        feq[i] = weights[i] * rho[index(x, y, z)] * (1.0 + dot_product *3.0 
            + dot_product * dot_product*4.5 - a*1.5);
    }
    //return &feq;
}

double LBM::calculate_feq(int x, int y, int z, int i)
{
    double dot_product = (double)vel[index(x,y,z)].x * (double)ei[i].x + (double)vel[index(x,y,z)].y * (double)ei[i].y +
                         (double)vel[index(x,y,z)].z * (double)ei[i].z;

    double c_s_squ= 1.0/3.0;
    double feq = weights[i] * rho[index(x,y,z)] * (1.0 + dot_product / (c_s * c_s) + dot_product * dot_product / (2 * c_s * c_s * c_s * c_s) - vel[index(x,y,z)].norm_square() / (2 * c_s * c_s));
    return feq;
}

double LBM::calculate_feq(int x, int y, int z, int i, double u_le_x)
{
    double dot_product = (vel[index(x,y,z)].x + u_le_x) * ei[i].x + vel[index(x,y,z)].y * ei[i].y +
                          vel[index(x,y,z)].z * ei[i].z;
    double norm_square = (vel[index(x,y,z)].x + u_le_x) * (vel[index(x,y,z)].x + u_le_x) + vel[index(x,y,z)].y * ei[i].y +
                          vel[index(x,y,z)].z * ei[i].z;

    double feq = weights[i] * rho[index(x,y,z)] * (1.0 + dot_product / (c_s * c_s) + dot_product * dot_product / (2 * c_s * c_s * c_s * c_s) - norm_square / (2 * c_s * c_s));
    return feq;
}

void LBM::calculate_fg(int x, int y, int z, double *feq, double *fg)
{
    for (int i = 0; i < direction_size; i++)
    {
        double dot_product = gf.x * ((double)ei[i].x - (double)vel[index(x, y, z)].x) + gf.y * ((double)ei[i].y - (double)vel[index(x, y, z)].y) +
            gf.z * ((double)ei[i].z - (double)vel[index(x, y, z)].z);
        //double c_s_squ = 1.0 / 3.0;
    fg[i] = (1.0 - 1.0 / (2 * tau)) *3.0 * dot_product * feq[i];
    }
}

//! calculate_fg for IB
void LBM::calculate_fg_IB(int x, int y, int z, double *feq, double *fg)
{
    for (int i = 0; i < direction_size; i++)
    {
        double dot_product = f_IB[index(x, y, z)].x * ((double)ei[i].x - (double)vel[index(x, y, z)].x) + f_IB[index(x, y, z)].y * ((double)ei[i].y - (double)vel[index(x, y, z)].y) +
            f_IB[index(x, y, z)].z * ((double)ei[i].z - (double)vel[index(x, y, z)].z);
        //double c_s_squ = 1.0 / 3.0;
        fg[i] = (1.0 - 1.0 / (2 * tau)) * 3.0 * dot_product * feq[i];
    }
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


void LBM::set_velocity_set(std::string velocity_set)
{

    if(velocity_set == "D3Q15")
    {
//{ 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14 },
//{ 0,  1, -1,  0,  0 , 0 , 0 , 1, -1,  1 , -1, 1, -1,  1, -1 },
//{ 0 , 0 , 0 , 1 ,-1 , 0 , 0 , 1 , 1 ,-1 ,-1 , 1 , 1 ,-1 ,-1 },
//{ 0 , 0 , 0 , 0 , 0 , 1 ,-1 , 1 , 1 , 1 , 1 ,-1 ,-1 ,-1 ,-1 },
        direction_size = 15;
        send_direction_size = 5;
        ei = new vector3<int>[15]{ {0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
            {1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},{1,1,-1},{-1,1,-1},{1,-1,-1},{-1,-1,-1} };
        weights = new double[15] { 2.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,
            1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0  };
    }else if(velocity_set=="D3Q19")
    {
        direction_size = 19;
        send_direction_size = 5;
//        ei = new vector3<int>[19]{{0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
//            {1,1,0},{-1,-1,0},{1,0,1},{-1,0,-1},{0,1,1},{0,-1,-1},{1,-1,0},{-1,1,0},{1,0,-1},
//            {-1,0,1},{0,1,-1},{0,-1,1}};
        ei = new vector3<int>[19]{ {0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},
            {1,1,0},{1,-1,0},{-1,1,0},{-1,-1,0},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1},{0,1,1},
            {0,1,-1},{0,-1,1},{0,-1,-1} };
//        c   data ex / 1,-1, 0, 0, 0, 0, 1, 1,-1, -1, 1, -1, 1, -1, 0, 0,  0, 0 / ,
//        c     &  ey / 0, 0, 1,-1, 0, 0, 1,-1, 1, -1, 0, 0,  0,  0, 1, 1, -1, -1 / ,
//        c     &  ez / 0, 0, 0, 0, 1,-1, 0, 0, 0,  0, 1, 1, -1, -1, 1, -1, 1, -1 /
        weights = new double[19]{1.0/3.0,1.0/18.0,1.0/18.0,1.0/18.0,
            1.0/18.0,1.0/18.0,1.0/18.0,1.0/36.0,1.0/36.0,1.0/36.0,
            1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,
            1.0/36.0,1.0/36.0,1.0/36.0};
    }else if (velocity_set == "D3Q27")
    {
        direction_size = 27;
        send_direction_size = 9;
        ei = new vector3<int>[27]{{0,0,0}, {1,0,0},{0,1,0},{-1,0,0},{0,-1,0},{0,0,1},{0,0,-1},
        {1,1,0},{-1,1,0},{-1,-1,0},{1,-1,0},{1,0,1},{0,1,1},{-1,0,1},{0,-1,1},{1,0,-1},{0,1,-1},{-1,0,-1},{0,-1,-1},
        {1,1,1},{-1,1,1},{-1,-1,1},{1,-1,1},{1,1,-1},{-1,1,-1},{-1,-1,-1},{1,-1,-1}};
        weights = new double[27] {8.0/27.0,2.0/27.0,2.0/27.0,2.0/27.0,
            2.0/27.0,2.0/27.0,2.0/27.0, 1.0/54.0,1.0/54.0,1.0/54.0
            ,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,
            1.0/54.0,1.0/54.0,1.0/54.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0
            ,1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0};
    } else if(velocity_set == "D2Q9")
    {
        direction_size = 9;
        send_direction_size = 3;
        ei = new vector3<int>[9]{ {0,0,0},{1,0,0},{0,1,0},{-1,0,0},{0,-1,0},{1,1,0},{-1,1,0},{-1,-1,0},{1,-1,0}};
        weights = new double[9]{ 4.0 / 9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};
    }
    lookup_reverse();
}

void LBM::lookup_reverse() {
    reverse_indexes = new int[direction_size];
    for(int i = 0; i < direction_size; i++)
    {
        for(int j = 0; j < direction_size; j++)
        {
            if(ei[i].x == -ei[j].x && ei[i].y == -ei[j].y && ei[i].z == -ei[j].z)
            {
                reverse_indexes[i] = j;
            }
        }
    }
}

void LBM::output_array(double *array)
{
    std::cout << "x,y,z value" << std::endl;
    for(int x = 0; x < NX; x++)
    {
        for(int y = 0; y < NY; y++)
        {
            for(int z = 0; z < NZ; z++)
            {
                std::cout << x << "," << y << "," << z << ": " << array[this->index(x,y,z)] << std::endl;
            }
        }
    }
}

void LBM::output_velocity_profile_at_a_section()
{
    int z_index = 0;
    for(int x = 0; x < NX; x++)
    {
        for(int y = 0; y < NY; y++)
        {
            std::cout << vel[index(x,y,z_index)].x << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void LBM::set_velocity(int x_field, int y_field, int z_field, double u_x, double u_y, double u_z)
{
    vel[index(x_field, y_field, z_field)].x = u_x;
    vel[index(x_field, y_field, z_field)].y = u_y;
    vel[index(x_field, y_field, z_field)].z = u_z;

    for (int i = 0; i < direction_size; i++)
    {
        pdf2[index(x_field, y_field, z_field, i)] = calculate_feq(x_field, y_field, z_field, i);
        pdf[index(x_field, y_field, z_field, i)] = calculate_feq(x_field, y_field, z_field, i);
    }
}

void LBM::set_density(int x_field, int y_field, int z_field, double density)
{
    rho[index(x_field,y_field,z_field)] = density;
}

void LBM::compute_density_momentum(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax) 
{
    int y ; int z ; int i;
#pragma omp parallel for private(y,z,i)
    for (int x = xmin; x < xmax; x++)
    {
        for (int y = ymin; y < ymax; y++)
        {
            for (int z = zmin; z < zmax; z++)
            {
                double new_density = 0;
                //vector3<double> u;
                double a=0.0, b=0.0, c=0.0;
//#pragma omp parallel reduction(+:a,b,c) 
                for(int i = 0; i < direction_size; i++)
                {
                    new_density += pdf[index(x,y,z,i)];
                    a += pdf[index(x, y, z, i)] * ei[i].x;
                    b += pdf[index(x, y, z, i)] * ei[i].y;
                    c += pdf[index(x,y,z,i)] * ei[i].z;
                }

                //rho[index(x,y,z)] = new_density;
                //vel[index(x,y,z)].x = (a + 0.5 * gf.x) / new_density;
                //vel[index(x,y,z)].y = (b + 0.5 * gf.y) / new_density;
                //vel[index(x,y,z)].z = (c + 0.5 * gf.z) / new_density;
                vector3<double> gf_inner;
                if (IBFSI == 1)
                {
                    gf_inner = f_IB[index(x, y, z)];
                }
                else
                {
                    gf_inner = gf;
                }

                rho[index(x, y, z)] = new_density;
                vel[index(x, y, z)].x = (a + 0.5 * gf_inner.x) / new_density;
                vel[index(x, y, z)].y = (b + 0.5 * gf_inner.y) / new_density;
                vel[index(x, y, z)].z = (c + 0.5 * gf_inner.z) / new_density;
            }
        }
    }
}

void LBM::stream() 
{
    // int y, z ,i;
#pragma omp parallel for //private(y,z,i)
    for(int x = 0; x < NX; x++)
    {
        for(int y = 0; y < NY; y++)
        {
            for(int z = 0; z < NZ; z++)
            {
                for(int i = 0; i < direction_size; i++)
                {
                    //! 先缺省是周期性的
                    int xmd = (NX + x - (int)ei[i].x) % NX;
                    int ymd = (NY + y - (int)ei[i].y) % NY;
                    int zmd = (NZ + z - (int)ei[i].z) % NZ;
                    pdf[index(x,y,z,i)] = pdf2[index(xmd,ymd,zmd,i)];
                }
            }
        }
    }
}

void LBM::output_indices_file()
{
    std::ofstream output("output/indices.csv");
    output << "x,y,z" << '\n';
    for(int i = 0; i < NX; i++)
    {
        for(int j = 0; j < NY; j++)
        {
            for(int k = 0; k < NZ; k++)
            {
                output << i << "," << j << "," << k << '\n';
            }
        }
    }
    std::cout << std::endl;
    output.close();
}

void LBM::collision_zip(int MRT, int BF, int LES,int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
{
    if (IBFSI == 1 && MRT == 0) collision_IB(xmin, xmax, ymin, ymax, zmin, zmax);
    if (IBFSI == 1 && MRT == 1) MRTcollision_IB(xmin, xmax, ymin, ymax, zmin, zmax);

    if (LES == 0 && BF == 0 && MRT == 0)  collision(xmin, xmax, ymin, ymax, zmin, zmax);
    if (LES == 0 && BF == 1 && MRT == 0) collision_BF(xmin, xmax, ymin, ymax, zmin, zmax);
    if (LES == 1) collision_LES(xmin, xmax, ymin, ymax, zmin, zmax);

    if (MRT == 1) MRTcollision(LES, xmin, xmax, ymin, ymax, zmin, zmax);
}


void LBM::collision(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
{   //! Performs the collision step. //设置collision执行范围 12.17
    double tauinv = 1.0 / tau;
    double omtauinv = 1.0-tauinv;     // 1 - 1/tau

    double *feq = new double[direction_size]();

    //! firstprivate(feq,fg)//ordered //flush(feq,fg)//schedule(dynamic) //#pragma omp flush(a,b)  //firstprivate(feq)//schedule(dynamic)
    int i,j,y, z;
    double temp;
//#pragma omp parallel for ordered//private(i, feq)
#pragma omp parallel for //private(j,i)
    for(int x = xmin; x < xmax; x++)
    {
        for(int y = ymin; y < ymax; y++)
        {
            for(int z = zmin; z < zmax; z++)
            {
                //calculate_feq2(x, y, z, feq);
                double a = vel[index(x, y, z)].norm_square();

                for (int i = 0; i < direction_size; i++)
                {
                    double dot_product = (double)vel[index(x, y, z)].x * (double)ei[i].x + (double)vel[index(x, y, z)].y * (double)ei[i].y +
                        (double)vel[index(x, y, z)].z * (double)ei[i].z;
                    feq[i] = weights[i] * rho[index(x, y, z)] * (1.0 + dot_product * 3.0
                        + dot_product * dot_product*4.5 - a * 1.5);
                    pdf2[index(x, y, z, i)] = omtauinv * pdf[index(x, y, z, i)] + tauinv * feq[i];
                }
                
//                for (int j = 0; j < direction_size; j++) {
//                    pdf2[index(x, y, z, j)] = omtauinv * pdf[index(x, y, z, j)] + tauinv * feq[j];
//                }
                
            }
        }
    }
    delete [] feq; feq = nullptr;
}

void LBM::collision_BF(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
{   //! includes body force
    double tauinv = 1.0 / tau;
    double omtauinv = 1.0 - tauinv;    //! 1 - 1/tau
    double *fg = new double[direction_size]();
    double *feq = new double[direction_size]();
//#pragma omp parallel for
    for (int x = xmin; x < xmax; x++)
    {
        for (int y = ymin; y < ymax; y++)
        {
            for (int z = zmin; z < zmax; z++)
            {
                calculate_feq2(x, y, z, feq);
                calculate_fg(x, y, z, feq, fg);

                for (int i = 0; i < direction_size; i++)
                {
                    pdf2[index(x, y, z, i)] = omtauinv * pdf[index(x, y, z, i)] + tauinv * feq[i] + fg[i];
                }
            }
        }
    }
    delete [] feq;    feq = nullptr;
    delete [] fg;    fg = nullptr;
}

void LBM::collision_IB(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
{   //! includes body force
    double tauinv = 1.0 / tau;
    double omtauinv = 1.0 - tauinv;     //! 1 - 1/tau
    double *fg = new double[direction_size]();
    double *feq = new double[direction_size]();

    //int x, y, z, i;
#pragma omp parallel for //private(x,y,z,i)
    for (int x = xmin; x < xmax; x++)
    {
        for (int y = ymin; y < ymax; y++)
        {
            for (int z = zmin; z < zmax; z++)
            {
                calculate_feq2(x, y, z, feq);
                calculate_fg_IB(x, y, z, feq, fg);

                for (int i = 0; i < direction_size; i++)
                {
                    pdf2[index(x, y, z, i)] = omtauinv * pdf[index(x, y, z, i)] + tauinv * feq[i] + fg[i];
                }
            }
        }
    }
    delete [] feq;    feq = nullptr;
    delete [] fg;    fg = nullptr;
}

void LBM::collision_LES(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
{   //! includes LES SGS model, now "tau" is local
    double *fg = new double[direction_size]();
    double *feq = new double[direction_size]();
    double *fneq = new double[direction_size];

//#pragma omp for
    for (int x = xmin; x < xmax; x++)
    {
        for (int y = ymin; y < ymax; y++)
        {
            for (int z = zmin; z < zmax; z++)
            {
                calculate_feq2(x, y, z, feq);
                calculate_fg(x, y, z, feq, fg);
                int i = 0;
//#pragma omp parallel private(i)
                for (i = 0; i < direction_size; i++)
                {
                    fneq[i] = calculate_fneq(x, y, z, i, feq[i]);
                }
                double p_ij = calculate_pij(fneq);
                const double tau_t = 0.5 * (sqrt(tau * tau + 2 * sqrt(2) * c_s * c_s * p_ij / ( c_s * c_s * c_s * c_s)) - tau);
                const double omega = 1.0 / (tau + tau_t);
                for (int i = 0; i < direction_size; i++)
                {
                    pdf2[index(x, y, z, i)] = (1.0- omega) * pdf[index(x, y, z, i)] + omega * feq[i] + fg[i];
                }
            }
        }
    }
    delete [] feq;    feq = nullptr;
    delete [] fneq;    fneq = nullptr;
    delete [] fg;    fg = nullptr;
}

int LBM::get_time()
{
    return time;
}

void LBM::perform_timestep(int MRT, int BF, int LES)
{
    time++;
    compute_density_momentum(0, NX, 0, NY, 0, NZ);

    collision_zip(MRT, BF, LES, 0, NX, 0, NY, 0, NZ);
    //mpi_send();
    stream();
    //bouzidi(num_node2, node2, delta, force);    //
    wall_bounce_back();     //! 处理平直壁面无滑移（simple bounce-back）
    velocity_boundary_condition();    //! 速度边界条件
    pressure_boundary_condition();    //! 压力边界条件
}

void LBM::output_all_fi(std::string filename)
{
    std::ofstream output_stream;
    output_stream.open(filename, std::ofstream::out);
    for (int x = 0; x < NX; x++)
    {
        for (int y = 0; y < NY; y++)
        {
            for (int z = 0; z < NZ; z++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                    output_stream << pdf[index(x, y, z, i)]<<" ";
                }
            }
        }
    }
    output_stream.close();
}

void LBM::output_tecplot(std::string filename, bool header)
{
    std::ofstream output_stream;
    output_stream.open (filename, std::ofstream::out );
    //! 若打开时有option：    | std::ofstream::app，则 app追加文件，不加为默认覆盖
    if(header) 
    {
        output_stream << "variables = x, y, z, rho, u, v, w" << '\n'; //, psi, psi1
        output_stream << "zone i="<<NX<< ",j="<<NY<<", k="<<NZ<<", f=point" << '\n';
    }

    double *psi,*psi1;
    psi = new double[NX*NY];
    psi1 = new double[NX*NY];
// if (0 == 1) {
//     int z = 0;
//     for (int x = 0; x < NX; x++) {
//         psi[index(x, 0, 0)] = 0.0;
//         for (int y = 1; y < NY; y++) {
//             psi[index(x, y, z)] = psi[index(x, y - 1, z)] + 0.5*(vel[index(x, y - 1, z)].x + vel[index(x, y, z)].x);
//         }
//     }
//     for (int y = 0; y < NY; y++) {
//         psi1[index(0, y, 0)] = 0.0;
//         for (int x = 1; x < NX; x++) {
//             psi1[index(x, y, z)] = psi1[index(x - 1, y, z)] + 0.5*(vel[index(x - 1, y, z)].y + vel[index(x, y, z)].y);
//         }
//     }
// }
    for (int z = 0; z < NZ; z++) 
    {
        for(int y = 0; y < NY; y++) 
        {
             for (int x = 0; x < NX; x++) 
            {
                output_stream << x << " " << y << " " << z << " " <<
                    rho[this->index(x, y, z)] << " " <<
                    vel[this->index(x, y, z)].x << " " <<
                    vel[this->index(x, y, z)].y << " " <<
                    vel[this->index(x, y, z)].z << " " << '\n';
        //        psi[this->index(x, y, z)] << " " <<-psi1[this->index(x, y, z)] <<'\n';
            }
        }
    }
    output_stream.close();
}

void LBM::output_history(std::string filename, bool header, int i, vector3<double> *force)
{
    std::ofstream output_stream;
    output_stream.open(filename, std::ofstream::app);
    //若打开时有option：    | std::ofstream::app，则 app追加文件，不加为默认覆盖
    if (header) {
        output_stream << i << "  " << force[0].x << "  " << force[0].y << "  " << force[0].z << '\n';
    }
    output_stream.close();
}

//
//void LBM::setNX(int a) {
//    NX = a;
//}
//int LBM::getNX() {
//    return NX;
//}

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