#include "LBMSolverOMP.hpp"
//! reference: MRT_3Dsimplest/

inline int indexM(int x, int y,  int NX)
{
    return  (y * NX) + x;
}

void LBM::iniMRT()
{
    double *Mtemp2 = new double[direction_size*direction_size];
    double *F_diag2 = new double[direction_size];

    if (velocity_set == "D2Q9")
    {
        double Mtemp[9][9] =
        {
            {   1.,  1., 1., 1.,  1., 1., 1., 1., 1.},    //! M11, M12, M13, M14, M15
            {  -4., -1.,-1.,-1., -1., 2., 2., 2., 2.},
            {   4., -2.,-2.,-2., -2., 1., 1., 1., 1.},
            {   0.,  1., 0.,-1.,  0., 1.,-1.,-1., 1.},
            {   0., -2., 0., 2.,  0., 1.,-1.,-1., 1.},
            {   0.,  0., 1., 0., -1., 1., 1.,-1.,-1.},
            {   0.,  0.,-2., 0.,  2., 1., 1.,-1.,-1.},
            {   0.,  1.,-1., 1., -1., 0., 0., 0., 0.},
            {   0.,  0., 0., 0.,  0., 1.,-1., 1.,-1.}
        };
        double F_diag[9] = { 9.0, 36.0,36.0,6.0,12.0,6.0,12.0,4.0,4.0 };
        for (int i = 0; i < direction_size; i++)
        {
            F_diag2[i] = F_diag[i];
            for (int j = 0; j < direction_size; j++)
            {
                Mtemp2[indexM(i, j, direction_size)] = Mtemp[i][j];
            }
        }
    }
    else if (velocity_set == "D3Q15")
    {
        double Mtemp[15][15] =
        {
            //{  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14 },
            { 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1},
            {-2 ,-1, -1 ,-1 ,-1 ,-1 ,-1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1},
            {16 ,-4, -4 ,-4 ,-4 ,-4 ,-4 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1},
            { 0 , 1 ,-1 , 0 , 0 , 0 , 0 , 1 ,-1 , 1 ,-1 , 1 ,-1 , 1 ,-1},
            { 0 ,-4 , 4 , 0 , 0 , 0 , 0 , 1 ,-1 , 1 ,-1 , 1 ,-1 , 1 ,-1},
            { 0 , 0 , 0 , 1 ,-1 , 0 , 0 , 1 , 1 ,-1 ,-1 , 1 , 1 ,-1 ,-1},
            { 0 , 0 , 0 ,-4 , 4 , 0 , 0 , 1 , 1 ,-1 ,-1 , 1 , 1 ,-1 ,-1},
            { 0 , 0 , 0 , 0 , 0 , 1 ,-1 , 1 , 1 , 1 , 1 ,-1 ,-1 ,-1 ,-1},
            { 0 , 0 , 0 , 0 , 0 ,-4 , 4 , 1 , 1 , 1 , 1 ,-1 ,-1 ,-1 ,-1},
            { 0 , 2 , 2 ,-1 ,-1 ,-1 ,-1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0},
            { 0 , 0 , 0 , 1 , 1 ,-1 ,-1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0},
            { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 ,-1 ,-1 , 1 , 1 ,-1 ,-1 , 1},
            { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 1 ,-1 ,-1 ,-1 ,-1 , 1 , 1},
            { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 ,-1 , 1 ,-1 ,-1 , 1 ,-1 , 1},
            { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 ,-1 ,-1 , 1 ,-1 , 1 , 1 ,-1}
        };
        double F_diag[15] = { 15.,18.,360.,10.,40.,10.,40.,10.,40.,12.,4.,8.,8.,8.,8. };
        for (int i = 0; i < direction_size; i++)
        {
            F_diag2[i] = F_diag[i];
            for (int j = 0; j < direction_size; j++)
            {
                Mtemp2[indexM(i, j, direction_size)] = Mtemp[i][j];
            }
        }
    }
    else if (velocity_set == "D3Q19")
    {
        double Mtemp[19][19] =
        {
            //{  0   1   2   3   4   5   6  7   8  9 10 11 12 13 14 15 16 17 18},
            {  1,  1,  1,  1,  1,  1,  1, 1,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
            {-30,-11,-11,-11,-11,-11,-11, 8,  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8},
            { 12, -4, -4, -4, -4, -4, -4, 1,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
            {  0,  1, -1,  0,  0,  0,  0, 1,  1,-1,-1, 1,-1, 1,-1, 0, 0, 0, 0},    //! ex
            {  0, -4,  4,  0,  0,  0,  0, 1,  1,-1,-1, 1,-1, 1,-1, 0, 0, 0, 0},
            {  0,  0,  0,  1, -1,  0,  0, 1, -1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1},    //! ey
            {  0,  0,  0, -4,  4,  0,  0, 1, -1, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1},
            {  0,  0,  0,  0,  0,  1, -1, 0,  0, 0, 0, 1, 1,-1,-1, 1,-1, 1,-1},    //! ez
            {  0,  0,  0,  0,  0, -4,  4, 0,  0, 0, 0, 1, 1,-1,-1, 1,-1, 1,-1},
            {  0,  2,  2, -1, -1, -1, -1, 1,  1, 1, 1, 1, 1, 1, 1,-2,-2,-2,-2},
            {  0, -4, -4,  2,  2,  2,  2, 1,  1, 1, 1, 1, 1, 1, 1,-2,-2,-2,-2},
            {  0,  0,  0,  1,  1, -1, -1, 1,  1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0},
            {  0,  0,  0, -2, -2,  2,  2, 1,  1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0},
            {  0,  0,  0,  0,  0,  0,  0, 1, -1,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
            {  0,  0,  0,  0,  0,  0,  0, 0,  0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1},
            {  0,  0,  0,  0,  0,  0,  0, 0,  0, 0, 0, 1,-1,-1, 1, 0, 0, 0, 0},
            {  0,  0,  0,  0,  0,  0,  0, 1,  1,-1,-1,-1, 1,-1, 1, 0, 0, 0, 0},
            {  0,  0,  0,  0,  0,  0,  0,-1,  1,-1, 1, 0, 0, 0, 0, 1, 1,-1,-1},
            {  0,  0,  0,  0,  0,  0,  0, 0,  0, 0, 0, 1, 1,-1,-1,-1, 1,-1, 1}
        };

        double F_diag[19] = { 19.0,2394.0,252.0,10.0,40.0,10.0,40.0,10.0,40.0,36.0,72.0,12.0,24.0,4.0,4.0,4.0,8.0,8.0,8.0 };
        for (int i = 0; i < direction_size; i++)
        {
            F_diag2[i] = F_diag[i];
            for (int j = 0; j < direction_size; j++)
            {
                Mtemp2[indexM(i, j, direction_size)] = Mtemp[i][j];
            }
        }
    }
    else if (velocity_set == "D3Q27")
    {   //! A D3Q27 multiple-relaxation-time lattice Boltzmann method for turbulent flows
        double Mtemp[27][27] =
        {
            {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
            {0,1,0,-1,0,0,0,1,-1,-1,1,1,0,-1,0,1,0,-1,0,1,-1,-1,1,1,-1,-1,1},    //! ex
            {0,0,1,0,-1,0,0,1,1,-1,-1,0,1,0,-1,0,1,0,-1,1,1,-1,-1,1,1,-1,-1},    //! ey
            {0,0,0,0,0,1,-1,0,0,0,0,1,1,1,1,-1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1},    //! ez
            {-2,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1},
            {0,2,-1,2,-1,-1,-1,1,1,1,1,1,-2,1,-2,1,-2,1,-2,0,0,0,0,0,0,0,0},
            {0,0,1,0,1,-1,-1,1,1,1,1,-1,0,-1,0,-1,0,-1,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,1,-1,1,-1,0,0,0,0,0,0,0,0,1,-1,1,-1,1,-1,1,-1},
            {0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0,-1,0,1,1,1,-1,-1,-1,-1,1,1},
            {0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0,-1,0,1,0,1,-1,-1,1,-1,1,1,-1},
            {0,-4,0,4,0,0,0,-1,1,1,-1,-1,0,1,0,-1,0,1,0,2,-2,-2,2,2,-2,-2,2},
            {0,0,-4,0,4,0,0,-1,-1,1,1,0,-1,0,1,0,-1,0,1,2,2,-2,-2,2,2,-2,-2},
            {0,0,0,0,0,-4,4,0,0,0,0,-1,-1,-1,-1,1,1,1,1,2,2,2,2,-2,-2,-2,-2},
            {0,4,0,-4,0,0,0,-2,2,2,-2,-2,0,2,0,-2,0,2,0,1,-1,-1,1,1,-1,-1,1},
            {0,0,4,0,-4,0,0,-2,-2,2,2,0,-2,0,2,0,-2,0,2,1,1,-1,-1,1,1,-1,-1},
            {0,0,0,0,0,4,-4,0,0,0,0,-2,-2,-2,-2,2,2,2,2,1,1,1,1,-1,-1,-1,-1},
            {4,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1},
            {-8,4,4,4,4,4,4,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,1,1,1,1,1,1,1,1},
            {0,-4,2,-4,2,2,2,1,1,1,1,1,-2,1,-2,1,-2,1,-2,0,0,0,0,0,0,0,0},
            {0,0,-2,0,-2,2,2,1,1,1,1,-1,0,-1,0,-1,0,-1,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,-2,2,-2,2,0,0,0,0,0,0,0,0,1,-1,1,-1,1,-1,1,-1},
            {0,0,0,0,0,0,0,0,0,0,0,0,-2,0,2,0,2,0,-2,1,1,-1,-1,-1,-1,1,1},
            {0,0,0,0,0,0,0,0,0,0,0,-2,0,2,0,2,0,-2,0,1,-1,-1,1,-1,1,1,-1},
            {0,0,0,0,0,0,0,1,-1,-1,1,-1,0,1,0,-1,0,1,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,-1,-1,1,1,0,1,0,-1,0,1,0,-1,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,1,-1,1,-1,-1,1,-1,1,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,1,-1,-1,1,-1,1}
        };

        double F_diag[27] = { 27.,18.,18,18,18,36,12,12,12,12,72,72,72,72,72,72,36,216,72,24,24,24,24,8,8,8,8 };
        for (int i = 0; i < direction_size; i++)
        {
            F_diag2[i] = F_diag[i];
            for (int j = 0; j < direction_size; j++)
            {
                Mtemp2[indexM(i, j, direction_size)] = Mtemp[i][j];
            }
        }
    }
        M_MRT = new double*[direction_size];
        M_MRTI = new double*[direction_size];
        MSM = new double*[direction_size];
        for (int i = 0; i < direction_size; i++)
        {
            M_MRT[i] = new double[direction_size];
            M_MRTI[i] = new double[direction_size];
            MSM[i] = new double[direction_size];
            for (int j = 0; j < direction_size; j++)
            {
                M_MRT[i][j] = Mtemp2[indexM(i, j, direction_size)];
            }
        }
        for (int j = 0; j < direction_size; j++)
        {
            for (int i = 0; i < direction_size; i++)
            {
                M_MRTI[j][i] = M_MRT[i][j] / F_diag2[i];
            }
        }

    calculate_MSM(tau);
}

void LBM::calculate_MSM(double tau)
{
    double *Inv_MS = new double[direction_size*direction_size];
    double *S_diag2 = new double[direction_size];

    if (velocity_set == "D2Q9")
    {
        double S_diag[9] = { 1.64, 1.64, 1.54, 1.6, 1.9, 1.6, 1.9, 1., 1. };    //! s0, ...s8, s8 = s7, s2 = 1.54, s4 = s6 = 1.9
        S_diag[7] = 1.0 / tau;
        S_diag[8] = S_diag[7];
        for (int i = 0; i < direction_size; i++)
        {
            S_diag2[i] = S_diag[i];
        }
    }
    else if (velocity_set == "D3Q15")
    {
        double S_diag[15] = { 0.0, 1.6,1.2,0.0,1.6,0.0,1.6,0.0,1.6, 1.0,1.0,1.0,1.0,1.0,1.2};
        S_diag[9] = 1.0 / tau;
        S_diag[10] = S_diag[9];
        S_diag[11] = S_diag[9];
        S_diag[12] = S_diag[9];
        S_diag[13] = S_diag[9];
        for (int i = 0; i < direction_size; i++)
        {
            S_diag2[i] = S_diag[i];
        }
    }
    else if (velocity_set == "D3Q19")
    {
        double S_diag[19] = { 0.0, 1.19,1.4,0.0,1.2,0.0,1.2,0.0,1.2,1.0,1.4,1.0,1.4,1.0,1.0,1.0,1.98,1.98,1.98 };
        S_diag[9] = 1.0 / tau;
        S_diag[11] = S_diag[9];
        S_diag[13] = S_diag[9];
        S_diag[14] = S_diag[9];
        S_diag[15] = S_diag[9];
        for (int i = 0; i < direction_size; i++)
        {
            S_diag2[i] = S_diag[i];
        }
    }
    else if (velocity_set == "D3Q27")
    {
        double s4, s5, s7, s10, s13, s16, s17, s18, s20, s23, s26;
        s5= 1.0 / tau;
        s7 = 1.0 / tau;
        s4 = 1.54; s10 = 1.5; s13 = 1.83; s16 = 1.4; s17 = 1.61;
        s18 = s20 = 1.98; s23 = s26 = 1.74;
        double S_diag[27] = { 0, 0, 0, 0, s4, s5, s5, s7, s7, s7, s10, s10, s10, s13, s13, s13, s16, s17, s18, s18, s20, s20, s20, s23, s23, s23, s26 };

        for (int i = 0; i < direction_size; i++)
        {
            S_diag2[i] = S_diag[i];
        }
    }

    for (int i = 0; i < direction_size; i++)
    {
        for (int j = 0; j < direction_size; j++)
        {
            Inv_MS[indexM(i, j, direction_size)] = -M_MRTI[i][j] * S_diag2[j];
        }
    }
    for (int j = 0; j < direction_size; j++)
    {
        for (int i = 0; i < direction_size; i++)
        {
            MSM[i][j] = 0.0;
            for (int k = 0; k < direction_size; k++)
            {
                MSM[i][j] = MSM[i][j] + Inv_MS[indexM(i, k, direction_size)] * M_MRT[k][j];
            }
        }
    }
}

void LBM::MRTcollision(int iLES, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
{
    double tauinv = 1.0 / tau;
    double omtauinv = 1.0 - tauinv;    //! 1 - 1/tau
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
                        //fg[i] = calculate_fg(x, y, z, i, feq[i]);
                        fneq[i] = pdf[index(x, y, z, i)] - feq[i];
                    }
                    double p_ij = calculate_pij(fneq);
                    const double tau_t = 0.5 * (sqrt(tau * tau + 2 * sqrt(2) * c_s * c_s * p_ij / (c_s * c_s * c_s * c_s)) - tau);
                    const double omega = 1.0 / (tau + tau_t);
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
        int y, z, i, k;
//#pragma omp parallel for private(y, z, i, k)
        for (int x = xmin; x < xmax; x++)
        {
            for (int y = ymin; y < ymax; y++)
            {
                for (int z = zmin; z < zmax; z++)
                {
                    for (int i = 0; i < direction_size; i++)
                    {
                        feq[i] = calculate_feq(x, y, z, i);
                        //fg[i] = calculate_fg(x, y, z, i, feq[i]);
                        fneq[i] = pdf[index(x, y, z, i)] - feq[i];
                    }
                    for (int i = 0; i < direction_size; i++)
                    {
                        colli[i] = 0.0;
                        for (int k = 0; k < direction_size; k++)
                        {
                            colli[i] = colli[i] + MSM[i][k] * fneq[k];
                        }
                        pdf2[index(x, y, z, i)] = pdf[index(x, y, z, i)] + colli[i];    //! +fg[i];
                    }
                }
            }
        }
    }

    delete [] colli;    colli = nullptr;
    delete [] fneq;    fneq = nullptr;
    delete [] feq;    feq = nullptr;
    delete [] fg;    fg = nullptr;
    //double* m_eq = new double[direction_size];
    //double* Dm = new double[direction_size];
    /*    for (int x = 0; x < NX; x++) {
            for (int y = 0; y < NY; y++) {
                    D2Q9calculate_meq(x, y,  m_eq);
                    for (int i = 0; i < direction_size; i++)
                    {
                        mmt[i] = 0.0;
                        for (int k = 0; k < direction_size; k++)
                        {
                            mmt[i] = mmt[i] + Mm[i][k] * pdf[index(x, y, z, k)];
                        }
                    }
                    for (int i = 0; i < direction_size; i++)
                    {
                        Dm[i] = 0.0;
                        for (int k = 0; k < direction_size; k++)
                        {
                            Dm[i] = Dm[i] + Inv_MS[i][k] * (mmt[k] - m_eq[k]);
                        }
                    }
                    for (int i = 0; i < direction_size; i++)
                    {
                        pdf2[index(x, y, z, i)] = pdf[index(x, y, z, i)] + Dm[i];
                    }
            }
        }
        */
}

//    double S_diag[9] = { 1.0, 1.6, 1.5, 1.0, 1.2, 1.0, 1.2, 1., 1. };// !s0, ...s8, s8 = s7, s2 = 1.54, s4 = s6 = 1.9
void LBM::MRTcollision_IB(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
{
    double tauinv = 1.0 / tau;
    double omtauinv = 1.0 - tauinv;    //! 1 - 1/tau
    double *fneq = new double[direction_size];
    double *feq = new double[direction_size];
    double *fg = new double[direction_size];
    double *colli = new double[direction_size];

    int x, y, z, i, ii, k;
#pragma omp parallel for private(x, y, z, i, ii, k)
    for (x = xmin; x < xmax; x++)
    {
        for (y = ymin; y < ymax; y++)
        {
            for (z = zmin; z < zmax; z++)
            {
                calculate_feq2(x, y, z, feq);
                calculate_fg_IB(x, y, z, feq, fg);
                for (i = 0; i < direction_size; i++)
                {
                    fneq[i] = pdf[index(x, y, z, i)] - feq[i];
                }
                for (ii = 0; ii < direction_size; ii++)
                {
                    colli[ii] = 0.0;
                    for (k = 0; k < direction_size; k++)
                    {
                        colli[ii] = colli[ii] + MSM[ii][k] * fneq[k];
                    }
                    pdf2[index(x, y, z, ii)] = pdf[index(x, y, z, ii)] + colli[ii] + fg[ii];
                }
            }
        }
    }

    delete [] colli;    colli = nullptr;
    delete [] fneq;    fneq = nullptr;
    delete [] feq;    feq = nullptr;
    delete [] fg;    fg = nullptr;
}

void LBM::D2Q9calculate_meq(int x, int y, double *m_eq)
{
    int z = 0;
    double u_squ = vel[index(x, y, z)].x*vel[index(x, y, z)].x;
    double v_squ = vel[index(x, y, z)].y*vel[index(x, y, z)].y;
    double jx = vel[index(x, y, z)].x*rho[index(x, y, z)];
    double jy = vel[index(x, y, z)].y*rho[index(x, y, z)];

    m_eq[0] = rho[index(x, y, z)];
    m_eq[1] = -2.0*rho[index(x, y, z)] + 3.0*(jx*jx + jy * jy) / rho[index(x, y, z)];    //! energy
    m_eq[2] = rho[index(x, y, z)] - 3.0*(jx*jx + jy * jy) / rho[index(x, y, z)];    //! energy square
    m_eq[3] = jx;    //! momentum in x - dir
    m_eq[4] = -jx;    //! energy flux in x - dir
    m_eq[5] = jy;    //! momentum in y - dir
    m_eq[6] = -jy;    //! energy flux in y - dir
    m_eq[7] = (jx*jx - jy * jy) / rho[index(x, y, z)];    //! diagonal comp of stress tensor
    m_eq[8] = jx * jy / rho[index(x, y, z)];    //! off diagonal comp of stress tensor

    /*double jx = rho_[n] * field_.u[n][0];
    double jy = rho_[n] * field_.u[n][1];

    mEq_[n][0] = rho_[n];
    mEq_[n][1] = -2.0 * rho_[n] + 3.0 * (jx * jx + jy * jy);
    mEq_[n][2] = rho_[n] - 3.0 * (jx * jx + jy * jy);
    mEq_[n][3] = jx;
    mEq_[n][4] = -jx;
    mEq_[n][5] = jy;
    mEq_[n][6] = -jy;
    mEq_[n][7] = (jx * jx - jy * jy);
    mEq_[n][8] = jx * jy;
            S_eq(1) = 0.d0
        c    S_eq(2) = 6.d0*(ux(i, j)*Fx(i, j) + uy(i, j)*Fy(i, j))
        c    S_eq(3) = -S_eq(2)
        c    S_eq(4) = Fx(i, j)
        c    S_eq(5) = -Fx(i, j)
        c    S_eq(6) = Fy(i, j)
        c    S_eq(7) = -Fy(i, j)
        c    S_eq(8) = 2.d0*(ux(i, j)*Fx(i, j) - uy(i, j)*Fy(i, j))
        c    S_eq(9) = (ux(i, j)*Fy(i, j) + uy(i, j)*Fx(i, j))*/
}
/*        MRTS1 = 1.19d0
        MRTS2 = 1.4d0
        MRTS4 = 1.2d0
        MRTS9 = 1.d0 / tau
        MRTS10 = 1.4d0
        MRTS13 = MRTS9
        MRTS16 = 1.98d0       !correct
        */

void LBM::D3Q19calculate_meq(int x, int y, int z,  double *m_eq)
{
    double u_squ = vel[index(x, y, z)].x*vel[index(x, y, z)].x;
    double v_squ = vel[index(x, y, z)].y*vel[index(x, y, z)].y;
    double w_squ = vel[index(x, y, z)].z*vel[index(x, y, z)].z;
    //we = 3.d0  !0.d0
    //wej = -11.d0 / 2.d0 !- 475.d0 / 63.d0 !
    //wxx = -0.5d0     !0.d0 !
    //meq(2) = we * rho(x, y, z) + wej / density * (u_xsqu + u_ysqu + u_zsqu)
    //meq(10) = wxx * meq(9)  !!!!!
    //meq(11) = (u_ysqu - u_zsqu) / density
    //meq(12) = wxx * meq(11)

    m_eq[0] = rho[index(x, y, z)];
    m_eq[1] = (-11.0 * rho[index(x, y, z)] + 19.0 * (u_squ+ v_squ + w_squ));    //! density
    m_eq[2] = (-475.0 / 63.0) * (u_squ + v_squ + w_squ); 
    m_eq[3] = vel[index(x, y, z)].x;
    m_eq[4] = (-2.0 / 3.0) * vel[index(x, y, z)].x;
    m_eq[5] = vel[index(x, y, z)].y;
    m_eq[6] = (-2.0 / 3.0) * vel[index(x, y, z)].y;
    m_eq[7] = vel[index(x, y, z)].z;
    m_eq[8] = (-2.0 / 3.0) * vel[index(x, y, z)].z;
    m_eq[9] = 2.0 * u_squ -( v_squ + w_squ);    //! density
    m_eq[10] = 0.0;    //! wxx * meq(9)  !!!!!
    m_eq[11] = (v_squ - w_squ);    //! density
    m_eq[12] = 0.0;    //! wxx * meq(11)
    m_eq[13] = vel[index(x, y, z)].x * vel[index(x, y, z)].y;    //! *density
    m_eq[14] = vel[index(x, y, z)].y * vel[index(x, y, z)].z;    //! *density
    m_eq[15] = vel[index(x, y, z)].z * vel[index(x, y, z)].x;    //! *density
    m_eq[16] = 0.0;
    m_eq[17] = 0.0;
    m_eq[18] = 0.0;
}