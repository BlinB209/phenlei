#include <cmath>
#include <direct.h>
#include <sstream>
#include <omp.h>
#include <IB.hpp>
#include "vector.hpp"
//#include "vector.cpp"
//#include "..\..\LBMSolverOMP\src\vector3.cpp"

void IBFSI::pIB(vector3<double>* uvw, double* rho, int time)
{
#ifdef USE_OMP
    omp_set_num_threads(numsOMP);
#endif

    vector<vector<vector<double>>> fx_temp(Nz, vector<vector<double>>(Ny, vector<double>(Nx, 0)));
    vector<vector<vector<double>>> fy_temp(Nz, vector<vector<double>>(Ny, vector<double>(Nx, 0)));
    vector<vector<vector<double>>> fz_temp(Nz, vector<vector<double>>(Ny, vector<double>(Nx, 0)));    //! 初始化体积力fxfyfz
    fx = fx_temp; fy = fy_temp; fz = fz_temp;

    velocity_Conversion_vector3_To_vector(uvw, rho);

    interpolate_to_lag_nodes();

    vector<double> Xelem(Nelem, 0), Yelem(Nelem, 0), Zelem(Nelem, 0), areaElem(Nelem, 0);
    vector<double> Uelem(Nelem, 0), Velem(Nelem, 0), Welem(Nelem, 0);
    vector<double> XIBelem(Nelem, 0), YIBelem(Nelem, 0);
    vector<double> UIBelem(Nelem, 0), VIBelem(Nelem, 0);
    get_lag_elements(&Xelem, &Yelem, &XIBelem, &YIBelem, &Uelem, &Velem ,&areaElem);

    vector<double> Fxelem(Nelem, 0);
    vector<double> Fyelem(Nelem, 0);
    vector<double> Fxelem_temp(Nelem, 0);
    vector<double> Fyelem_temp(Nelem, 0);

    double MAXerror = 1.0;
    int iterNum = 0;
    int totalIter = 1;
    double dtolerror = 0.001;
    while (iterNum < totalIter && MAXerror > dtolerror)
    {
        double alpha = 0.0, beta = 2.0;

        interpolate_to_lag_elements(&UIBelem, &VIBelem, Xelem, Yelem);

        calculate_lagrangian_element_forces(&Fxelem_temp, &Fyelem_temp, alpha, beta, Xelem, Yelem, XIBelem, YIBelem, Uelem, Velem, UIBelem, VIBelem, areaElem);

        calculate_body_forces_on_eular(&Fxelem, &Fyelem, Xelem, Yelem, Fxelem_temp, Fyelem_temp, areaElem);

        double dtol_error = Uref * Nelem;

        MAXerror = 0.0;
        for (int i = 0; i < Nelem; i++)
        {
            MAXerror += sqrt((Uelem[i] - UIBelem[i]) * (Uelem[i] - UIBelem[i]) + (Velem[i] - VIBelem[i]) * (Velem[i] - VIBelem[i]));
        }
        MAXerror = MAXerror / dtol_error;
        iterNum = iterNum + 1;
    }
    calculate_node_forces(Fxelem, Fyelem);

    write_everything(time, CdCl_name);

    /*
    if (time % 50 == 0 && iMotion == 1) {
        cout << "IB_iterNum = " << iterNum << " ;  " << "IB_iterError = " << MAXerror << " ;  Cd = " << Cd_ib << " ;  Cl = " << Cl_ib << endl;
    }
    else if (time % 50 == 0 && iMotion == 2) {
        cout << "IB_iterNum = " << iterNum << " ;  " << "IB_iterError = " << MAXerror << " ;  Cd = " << Cd_ib << " ;  Cl = " << Cl_ib << " ;  xtail = " << Xnode[Nnode - 1] << " ;  ytail = " << Ynode[Nnode - 1] << endl;
    }
    else if(time % 50 == 0 && iMotion == 3) {
        cout << "IB_iterNum = " << iterNum << " ;  " << "IB_iterError = " << MAXerror << " ;  Cd = " << Cd_ib << " ;  Cl = " << Cl_ib << " ;  VIV_y = " << VIV_X[0] << " ;  VIV_theta = " << VIV_X[1] << endl;
    }
    else {
    }
    */

    int iternums = 200;
    switch (iMotion)
    {
        case 1:
            prescribed_motion(time, subMotionType1);
            break;
        case 2:
            filament_motion(iternums);
            break;
        case 3:
            VIV_motion();
            break;
        default:
            cout << "" << endl;
    }
}

void IBFSI::interpolate_to_lag_nodes()
{
    int x1, y1;
    if (iMotion == 1 && subMotionType1 == 999)
    {
        for (int iii = 0; iii < Nnode; ++iii)
        {
            /*RhoIBnode[iii] = 0.0;
            UIBnode[iii] = 0.0;
            VIBnode[iii] = 0.0;
            WIBnode[iii] = 0.0;*/
            XIBnode[iii] = Xnode[iii];
            YIBnode[iii] = Ynode[iii];
            ZIBnode[iii] = Znode[iii];
        }
    }

    for (int n = 0; n < Nnode; ++n)
    {
        RhoIBnode[n] = 0.0;
        UIBnode[n] = 0.0;
        VIBnode[n] = 0.0;
        WIBnode[n] = 0.0;

        x1 = (int)(fabs(Xnode[n] - xGrid[0]) / dx);
        y1 = (int)(fabs(Ynode[n] - yGrid[0]) / dy);

        double dist_x, dist_y;
        double wt_x, wt_y;
//#ifdef USE_OMP
//#pragma omp parallel for
//#endif
        int ss = (int)(supp/2);
        for (int j = y1 - ss; j <= y1 + ss + 1; ++j)
        {
            for (int i = x1 - ss; i <= x1 + ss + 1; ++i)
            {
                dist_x = fabs((Xnode[n] - xGrid[i]) / dx);
                dist_y = fabs((Ynode[n] - yGrid[j]) / dy);

                //! 求phi_x
                if (dist_x < 1.0)
                {
                    wt_x = (3.0 - 2.0 * dist_x + sqrt(1.0 + 4.0 * dist_x - 4.0 * dist_x * dist_x)) / 8.0;
                }
                else if (dist_x < 2.0)
                {
                    wt_x = (5.0 - 2.0 * dist_x - sqrt(-7.0 + 12.0 * dist_x - 4.0 * dist_x * dist_x)) / 8.0;
                }
                else
                {
                    wt_x = 0.0;
                }
                //! 求phi_y
                if (dist_y < 1.0)
                {
                    wt_y = (3.0 - 2.0 * dist_y + sqrt(1.0 + 4.0 * dist_y - 4.0 * dist_y * dist_y)) / 8.0;
                }
                else if (dist_y < 2.0)
                {
                    wt_y = (5.0 - 2.0 * dist_y - sqrt(-7.0 + 12.0 * dist_y - 4.0 * dist_y * dist_y)) / 8.0;
                }
                else
                {
                    wt_y = 0.0;
                }
                UIBnode[n] += (u[j][i] * wt_x * wt_y);
                VIBnode[n] += (v[j][i] * wt_x * wt_y);
                RhoIBnode[n] += (rho[j][i] * wt_x * wt_y);
            }
        }
        XIBnode[n] += (UIBnode[n] * dt);
        YIBnode[n] += (VIBnode[n] * dt);
    }
}

void IBFSI::get_lag_elements(vector<double> *Xelem, vector<double> *Yelem, vector<double> *XIBelem, vector<double> *YIBelem, vector<double> *Uelem, vector<double> *Velem, vector<double> *areaElem)
{
    int i, I, J;
//#ifdef USE_OMP
//#pragma omp parallel for private(i,I,J)
//#endif
    for (i = 0; i < Nelem; ++i)
    {
        double x1, y1;
        double x2, y2;

        I = Ielem[i];
        J = Jelem[i];

        x1 = Xnode[I];
        y1 = Ynode[I];

        x2 = Xnode[J];
        y2 = Ynode[J];
        
        (*Xelem)[i] = (x1 + x2) / 2.0;
        (*Yelem)[i] = (y1 + y2) / 2.0;

        (*Uelem)[i] = (Unode[I] + Unode[J]) / 2.0;
        (*Velem)[i] = (Vnode[I] + Vnode[J]) / 2.0;

        (*XIBelem)[i] = (XIBnode[I] + XIBnode[J]) / 2.0;
        (*YIBelem)[i] = (YIBnode[I] + YIBnode[J]) / 2.0;

        double ax, ay;
        ax = x1 - x2;
        ay = y1 - y2;
        (*areaElem)[i] = sqrt(ax * ax + ay * ay);
    }
}

void IBFSI::interpolate_to_lag_elements(vector<double> *UIBelem, vector<double> *VIBelem, vector<double> Xelem, vector<double> Yelem)
{
    int x1, y1;
    for (int n = 0; n < Nelem; ++n)
    {
        x1 = (int)(fabs(Xelem[n] - xGrid[0]) / dx);
        y1 = (int)(fabs(Yelem[n] - yGrid[0]) / dy);

        (*UIBelem)[n] = 0.0;
        (*VIBelem)[n] = 0.0;

        double dist_x, dist_y;
        double wt_x, wt_y;

    //#ifdef USE_OMP
    //#pragma omp parallel for
    //#endif
        int ss = (int)(supp/2);
        for (int j = y1 - ss; j <= y1 + ss + 1; ++j)
        {
            for (int i = x1 - ss; i <= x1 + ss + 1; ++i)
            {
                dist_x = fabs((Xelem[n] - xGrid[i]) / dx);
                dist_y = fabs((Yelem[n] - yGrid[j]) / dy);

                //! 求phi_x
                if (dist_x < 1.0)
                {
                    wt_x = (3.0 - 2.0 * dist_x + sqrt(1.0 + 4.0 * dist_x - 4.0 * dist_x * dist_x)) / 8.0;
                }
                else if (dist_x < 2.0)
                {
                    wt_x = (5.0 - 2.0 * dist_x - sqrt(-7.0 + 12.0 * dist_x - 4.0 * dist_x * dist_x)) / 8.0;
                }
                else
                {
                    wt_x = 0.0;
                }
                //! 求phi_y
                if (dist_y < 1.0)
                {
                    wt_y = (3.0 - 2.0 * dist_y + sqrt(1.0 + 4.0 * dist_y - 4.0 * dist_y * dist_y)) / 8.0;
                }
                else if (dist_y < 2.0)
                {
                    wt_y = (5.0 - 2.0 * dist_y - sqrt(-7.0 + 12.0 * dist_y - 4.0 * dist_y * dist_y)) / 8.0;
                }
                else
                {
                    wt_y = 0.0;
                }
                (*UIBelem)[n] += (u[j][i] * wt_x * wt_y);
                (*VIBelem)[n] += (v[j][i] * wt_x * wt_y);
            }
        }
    }
}

void IBFSI::calculate_lagrangian_element_forces(vector<double> *Fxelem_temp, vector<double> *Fyelem_temp, double alpha, double beta, 
    vector<double> Xelem, vector<double> Yelem, vector<double> XIBelem, vector<double> YIBelem, vector<double> Uelem,
    vector<double> Velem, vector<double> UIBelem, vector<double> VIBelem, vector<double> areaElem)
{
    for (int i = 0; i < Nelem; ++i)
    {
        (*Fxelem_temp)[i] = - 2.0 * alpha * areaElem[i] * (Xelem[i] - XIBelem[i])
            - 2.0 * beta * areaElem[i] * (Uelem[i] - UIBelem[i]);
        (*Fyelem_temp)[i] = - 2.0 * alpha * areaElem[i] * (Yelem[i] - YIBelem[i])
            - 2.0 * beta * areaElem[i] * (Velem[i] - VIBelem[i]);
    }
}

void IBFSI::calculate_body_forces_on_eular(vector<double> *Fxelem, vector<double> *Fyelem, vector<double> Xelem, vector<double> Yelem, vector<double> Fxelem_temp, vector<double> Fyelem_temp, vector<double> areaElem)
{
    vector<vector<vector<double>>> fx_temp(Nz, vector<vector<double>>(Ny, vector<double>(Nx, 0)));
    vector<vector<vector<double>>> fy_temp(Nz, vector<vector<double>>(Ny, vector<double>(Nx, 0)));    //! 初始化体积力fxfyfz

    int x1, y1;
    for (int n = 0; n < Nelem; ++n)
    {
        x1 = (int)(fabs(Xelem[n] - xGrid[0]) / dx);
        y1 = (int)(fabs(Yelem[n] - yGrid[0]) / dy);

        double dist_x, dist_y;
        double wt_x, wt_y;
    //#ifdef USE_OMP
    //#pragma omp parallel for
    //#endif
        int ss = (int)(supp/2);
        for (int j = y1 - ss; j <= y1 + ss + 1; ++j)
        {
            for (int i = x1 - ss; i <= x1 + ss + 1; ++i)
            {
                dist_x = fabs((Xelem[n] - xGrid[i]) / dx);
                dist_y = fabs((Yelem[n] - yGrid[j]) / dy);

                //! 求phi_x
                if (dist_x < 1.0)
                {
                    wt_x = (3.0 - 2.0 * dist_x + sqrt(1.0 + 4.0 * dist_x - 4.0 * dist_x * dist_x)) / 8.0;
                }
                else if (dist_x < 2.0)
                {
                    wt_x = (5.0 - 2.0 * dist_x - sqrt(-7.0 + 12.0 * dist_x - 4.0 * dist_x * dist_x)) / 8.0;
                }
                else
                {
                    wt_x = 0.0;
                }
                //! 求phi_y
                if (dist_y < 1.0)
                {
                    wt_y = (3.0 - 2.0 * dist_y + sqrt(1.0 + 4.0 * dist_y - 4.0 * dist_y * dist_y)) / 8.0;
                }
                else if (dist_y < 2.0)
                {
                    wt_y = (5.0 - 2.0 * dist_y - sqrt(-7.0 + 12.0 * dist_y - 4.0 * dist_y * dist_y)) / 8.0;
                }
                else
                {
                    wt_y = 0.0;
                }

                fx_temp[0][j][i] = fx_temp[0][j][i] - (Fxelem_temp[n] * wt_x * wt_y * areaElem[n] / (dx * dy));
                fy_temp[0][j][i] = fy_temp[0][j][i] - (Fyelem_temp[n] * wt_x * wt_y * areaElem[n] / (dx * dy));
                //fx_temp[0][j][i] = fx_temp[0][j][i] - ( Fxelem_temp[n] * wt_x * wt_y * ds0 / (dx * dy) );
                //fy_temp[0][j][i] = fy_temp[0][j][i] - ( Fyelem_temp[n] * wt_x * wt_y * ds0 / (dx * dy) );
            }
        }
    }

    //! 更新速度
    int j;
//#ifdef USE_OMP
//#pragma omp parallel for
//#endif
    for (j = 0; j < Ny; ++j)
    {
        int i;
        for (i = 0; i < Nx; ++i)
        {
            u[j][i] = u[j][i] + 0.5 * dt * fx_temp[0][j][i] / rho[j][i];
            v[j][i] = v[j][i] + 0.5 * dt * fy_temp[0][j][i] / rho[j][i];

            fx[0][j][i] = fx[0][j][i] + fx_temp[0][j][i];
            fy[0][j][i] = fy[0][j][i] + fy_temp[0][j][i];
        }
    }

//#ifdef USE_OMP
//#pragma omp parallel for
//#endif
    for (int i = 0; i < Nelem; ++i)
    {
        (*Fxelem)[i] += Fxelem_temp[i];
        (*Fyelem)[i] += Fyelem_temp[i];
    }
}

void IBFSI::calculate_node_forces(vector<double> Fxelem, vector<double> Fyelem)
{
    for (int iii = 0; iii < Nnode; ++iii)
    {
        Fxnode[iii] = 0.0;
        Fynode[iii] = 0.0;
        Fznode[iii] = 0.0;
    }

//#ifdef USE_OMP
//#pragma omp parallel for
//#endif
    for (int n = 0; n < Nelem; ++n)
    {
        int i = Ielem[n];
        int j = Jelem[n];

        Fxnode[i] = Fxnode[i] + (Fxelem[n] / 2.0);
        Fxnode[j] = Fxnode[j] + (Fxelem[n] / 2.0);
        Fynode[i] = Fynode[i] + (Fyelem[n] / 2.0);
        Fynode[j] = Fynode[j] + (Fyelem[n] / 2.0);
    }

    extfulx = Fxnode;
    extfuly = Fynode;
}

void IBFSI::write_everything(int time, string file_out)
{
    ofstream output_stream;
    output_stream.open(file_out, std::ofstream::app);

    //! 若打开时有option：    | std::ofstream::app，则 app追加文件，不加为默认覆盖
    Cd_ib = 0.0; Cl_ib = 0.0;
    for (int i = 0; i < Nnode; ++i)
    {
        Cd_ib += Fxnode[i] / Fref;
        Cl_ib += Fynode[i] / Fref;
    }

    if (time == 0)
    {
        if (iMotion == 3)
        {
            output_stream << "variables = " << "time " << "Cd  " << "Cl  " << "VIV_Y  " << "VIV_Theta" << '\n';
        }
        else if (iMotion == 2)
        {
            output_stream << "variables = " << "time " << "Cd  " << "Cl  " << "xtail  " << "ytail" << '\n';
        }
        else
        {
            output_stream << "variables = " << "time " << "Cd  " << "Cl" << '\n';
        }
    }

    if (iMotion == 3)
    {
        output_stream << time << "  " << Cd_ib << "  " << Cl_ib << "  " << VIV_X[0] << "  " << VIV_X[1] << '\n';
    }
    else if (iMotion == 2)
    {
        output_stream << time << "  " << Cd_ib << "  " << Cl_ib << "  " << Xnode[Nnode - 1] << "  " << Ynode[Nnode - 1] << '\n';
    }
    else
    {
        output_stream << time << "  " << Cd_ib << "  " << Cl_ib << '\n';
    }

    output_stream.close();
}

//! 初始化流体域，计算步长，分布函数影响区域，速度分布D..Q..
void IBFSI::initialise(int nx, int ny, int nz, double elesize0, double elesize1, double elesize2, vector3<double> inflow_Velocity, string outFileName)
{
    setNxyz(nx, ny, nz);
    setElExyz(elesize0, elesize1, elesize2);

    vel0 = inflow_Velocity;
    xGrid = creatmatrix_1D(Nx + 1, 0.0);
    yGrid = creatmatrix_1D(Ny + 1, 0.0);

    u = creatmatrix_double(Ny, Nx);
    v = creatmatrix_double(Ny, Nx);
    w = creatmatrix_double(Ny, Nx);
    rho = creatmatrix_double(Ny, Nx);

    vector<vector<vector<double>>> fx_t(Nz, vector<vector<double>>(Ny, vector<double>(Nx, 0)));
    vector<vector<vector<double>>> fy_t(Nz, vector<vector<double>>(Ny, vector<double>(Nx, 0)));
    vector<vector<vector<double>>> fz_t(Nz, vector<vector<double>>(Ny, vector<double>(Nx, 0)));
    fx = fx_t; fy = fy_t; fz = fz_t;    //! 初始化体积力fxfyfz

    CdCl_name = outFileName;
}

void IBFSI::setCharacterQuantity()
{
    double Lchod = 0.0;
    for (int i = 0; i < Nnode; i++)
    {
        for (int j = 0; j < Nnode; j++)
        {
            double lll = sqrt((Xnode[i] - Xnode[j]) * (Xnode[i] - Xnode[j]) + (Ynode[i] - Ynode[j]) * (Ynode[i] - Ynode[j]));
            if (lll > Lchod) Lchod = lll;
        }
    }
    if ((Lchod - 1.0) <= 0.01) Lchod = 1.0;
    //! 设置特征量
    Lref = Lchod;
    Uref = vel0.x;
    Tref = Lref / Uref;
    Fref = 0.5 * denIn * Uref * Uref * Lref;
    KBref = denIn * Uref * Uref * Lref * Lref * Lref;
    KSref = denIn * Uref * Uref * Lref;
    dt = ds0;
    dx = ds0; dy = ds0;
    //ds0 = Lref / Nelem;
}

void IBFSI::initialBasicParams()
{
    for (int i = 0; i < Nx + 1; i++)
    {
        xGrid[i] = (double(i) - double(Nx) * solid_zero_positionX_percent) * dx;
    }
    for (int j = 0; j < Ny + 1; j++)
    {
        yGrid[j] = (double(j) - double(Ny) * solid_zero_positionY_percent) * dy;
    }

    Unode = creatmatrix_1D(Nnode, 0.0);
    Vnode = creatmatrix_1D(Nnode, 0.0);
    Wnode = creatmatrix_1D(Nnode, 0.0);
    Axnode = creatmatrix_1D(Nnode, 0.0);
    Aynode = creatmatrix_1D(Nnode, 0.0);
    Aznode = creatmatrix_1D(Nnode, 0.0);
    Fxnode = creatmatrix_1D(Nnode, 0.0);
    Fynode = creatmatrix_1D(Nnode, 0.0);
    Fznode = creatmatrix_1D(Nnode, 0.0);
    RhoIBnode = creatmatrix_1D(Nnode, 0.0);
    UIBnode = creatmatrix_1D(Nnode, 0.0);
    VIBnode = creatmatrix_1D(Nnode, 0.0);
    WIBnode = creatmatrix_1D(Nnode, 0.0);

    if (iMotion == 2 && Ampl_heaving != 0)
    {
        for (int n = 0; n < Nnode; n++)
        {
            Ynode[n] -= Ampl_heaving;
        }
    }

    XIBnode = Xnode; YIBnode = Ynode; ZIBnode = Znode;
    X0node = Xnode; Y0node = Ynode; Z0node = Znode;
    X00node = Xnode; Y00node = Ynode; Z00node = Znode;
}

void IBFSI::getInputFromFile(string file_str)
{
    char buffer[256];
    ifstream infile(file_str);
    if (!infile)
    {
        cout << "Unable to open input complex geometry (shape) file:input_IB.dat";
        exit(1);
    }

    while (!infile.eof())
    {
        infile.getline(buffer, 60);    //! 跳过
        infile.getline(buffer, 60);    //! 跳过

        infile.getline(buffer, 100);
        sscanf_s(buffer, "%lf", &ds0);
        infile.getline(buffer, 100);
        sscanf_s(buffer, "(%lf  %lf)", &solid_zero_positionX_percent, &solid_zero_positionY_percent);
        infile.getline(buffer, 100);
        sscanf_s(buffer, "%d", &numsOMP);
        infile.getline(buffer, 100);
        sscanf_s(buffer, "%lf", &Re);
        cout << "Reynolds number : " << Re << endl;
        infile.getline(buffer, 100); 
        sscanf_s(buffer, "%d", &supp);
        infile.getline(buffer, 150);
        sscanf_s(buffer, "%d", &iMotion);
        infile.getline(buffer, 150);
        sscanf_s(buffer, "%d", &subMotionType1);
        infile.getline(buffer, 100);
        sscanf_s(buffer, "(%lf  %lf  %lf)", &gravityX, &gravityY, &gravityZ);

        infile.getline(buffer, 200);    //! 跳过

        infile.getline(buffer, 100);
        sscanf_s(buffer, "%d", &IsRotateOrnot);
        infile.getline(buffer, 200);
        sscanf_s(buffer, "%lf  %lf", &iniAngle, &iniPhase);
        infile.getline(buffer, 200);
        sscanf_s(buffer, "%lf", &A_rotate);
        infile.getline(buffer, 200);
        sscanf_s(buffer, "%lf", &T_rotate);

        infile.getline(buffer, 150);    //! 跳过

        infile.getline(buffer, 100);
        sscanf_s(buffer, "%lf  %lf", &kb_film, &ks_film);
        infile.getline(buffer, 150);
        sscanf_s(buffer, "%lf  %lf", &Ampl_heaving, &Ampl_pitching);
        infile.getline(buffer, 200);
        sscanf_s(buffer, "%lf  %lf", &Period_heaving, &Period_pitching);
        infile.getline(buffer, 150);
        sscanf_s(buffer, "(%lf  %lf)", &accfix_film1, &accfix_film2);

        infile.getline(buffer, 150);    //! 跳过

        infile.getline(buffer, 100);
        sscanf_s(buffer, "%lf", &VIV_MassRotio);
        infile.getline(buffer, 100);
        sscanf_s(buffer, "%lf  %lf", &VIV_NonDimNatFreq_y, &VIV_NonDimNatFreq_Ang);
        infile.getline(buffer, 150);
        sscanf_s(buffer, "%lf  %lf", &VIV_DampRatio_y, &VIV_DampRatio_Ang);
        infile.getline(buffer, 150);
        sscanf_s(buffer, "(%lf  %lf)", &VIVMassCenter_X, &VIVMassCenter_Y);
    }
    infile.close();
}

//! getLagFromMesh    从mesh文件中读取初始固体坐标
void IBFSI::getLagFromMesh(string meshFlieName)
{
    string line;
    /*char buffer[256];
    double buffer2[3] = { 0,0,0 };
    int buffer3[3] = { 0,0,0 };
    int* a;*/

    ifstream infile(meshFlieName);
    if (!infile)
    {
        cout << "Unable to open input complex geometry (shape) file:mesh_3d.dat or mesh_2d.dat";
        exit(1);
    }

    int pointNumber[3] = { 100,100,100 };

    getline(infile, line);    //! 跳过

    getline(infile, line);stringstream ss(line);
    int x; int nn = 0;
    while (ss >> x)
    {
        pointNumber[nn] = x;
        nn++;
    }
    Nnode = pointNumber[0]; Nelem = pointNumber[1];
    Xnode = creatmatrix_1D(Nnode, 0.0);
    Ynode = creatmatrix_1D(Nnode, 0.0);
    Znode = creatmatrix_1D(Nnode, 0.0);
    Ielem = creatmatrix_1D_int(Nelem, 0);
    Jelem = creatmatrix_1D_int(Nelem, 0);
    Kelem = creatmatrix_1D_int(Nelem, 0);

    getline(infile, line);    //! 跳过
    getline(infile, line);

    int iline = 0;
    double iXYZ[4] = { 100.,100.,100.,100. };
    while (iline < Nnode)
    {
        getline(infile, line);
        stringstream ss(line);
        nn = 0;
        double xx;
        while (ss >> xx)
        {
            iXYZ[nn] = xx;
            nn++;
        }
        Xnode[iline] = iXYZ[1];
        Ynode[iline] = iXYZ[2];
        Znode[iline] = iXYZ[3];
        iline++;
    }

    getline(infile, line);    //! 跳过
    getline(infile, line);

    iline = 0;
    int nIJK[6] = { 100,100,100,100,100,100 };
    while (iline < Nelem)
    {
        getline(infile, line);
        stringstream ss(line);
        x = 0; nn = 0;
        while (ss >> x)
        {
            nIJK[nn] = x;
            nn++;
        }
        Ielem[iline] = nIJK[1]-1;
        Jelem[iline] = nIJK[2]-1;
        Kelem[iline] = nIJK[3]-1;
        iline++;
    }

    infile.close();

    Xnode = add1D_num(Xnode, 0.0);
    Ynode = add1D_num(Ynode, 0.0);
}

//! 输出为asca码文件
void IBFSI::output_tecplot(string filename, bool header)
{
    std::ofstream output_stream;
    output_stream.open(filename, std::ofstream::out);
    //! 若打开时有option：    | std::ofstream::app，则 app追加文件，不加为默认覆盖
    if (header)
    {
        output_stream << "variables = x, y, z, uIB, vIB, wIB, Fx, Fy, Fz" << '\n';
    }
    for (int i = 0; i < Nnode; i++)
    {
        output_stream << Xnode[i] << " " << Ynode[i] << " " << Znode[i] << " " << UIBnode[i] << " " << VIBnode[i] << " " << WIBnode[i] << " " << Fxnode[i] << " " << Fynode[i] << " " << Fznode[i] << '\n';
        //        psi[this->index(x, y, z)] << " " <<-psi1[this->index(x, y, z)] <<'\n';
    }
    output_stream.close();
}
