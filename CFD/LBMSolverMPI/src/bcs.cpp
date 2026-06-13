#include "LBMSolverMPI.hpp"
#include "vector3.cpp"
#include <iomanip>
#include <iostream>

void LBM::wall_bounce_back()
{   //! nonslip壁面，bounce-back实现无滑移边界条件.
    for (int x = 0; x < NX; x++)
    {
        for (int y = 0; y < NY; y++)
        {
            for (int z = 0; z < NZ; z++)
            {
                if (obst[index(x, y, z)] == 3)
                {
                    for (int i = 0; i < direction_size; i++)
                    {
                        pdf[index(x, y, z, i)] = pdf[index(x, y, z, reverse_indexes[i])];
                    }
                }
            }
        }
    }
}

void LBM::velocity_boundary_condition()
{
    const double tauinv = 1.0 / tau;
    const double omtauinv = 1.0 - tauinv;
    std::string Astring;
    if (x_id == 0)
    {
        Astring = BC.xminFace;    //! 字符数组赋值给string，为字符串比较准备
        if (Astring == "velocity")
        {
            int x = 0;
            for (int y = 0; y < NY; y++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    vel[index(x, y, z)].x = BC.vel_specified[0].x;
                    vel[index(x, y, z)].y = BC.vel_specified[0].y;
                    vel[index(x, y, z)].z = BC.vel_specified[0].z;
                    rho[index(x, y, z)] = rho[index(x + 1, y, z)];    //! rho 由相邻流场内网格点(x+1, y, z)插值而来
                    for (int i = 0; i < direction_size; i++)
                    {
                        double feq = calculate_feq(x, y, z, i);
                        double feq1 = calculate_feq(x + 1, y, z, i);
                        pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x + 1, y, z, i)] - feq1);
                        //! 非平衡部分(pdf2[index(x+1, y, z, i)] - feq1) 由相邻流场内网格点(x+1, y, z)插值而来
                    }
                }
            }
        }
    }
    if (x_id == x_np - 1)
    {
        Astring = BC.xmaxFace;
        if (Astring == "velocity")
        {
            int x = NX - 1;
            for (int y = 0; y < NY; y++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    vel[index(x, y, z)].x = BC.vel_specified[1].x;
                    vel[index(x, y, z)].y = BC.vel_specified[1].y;
                    vel[index(x, y, z)].z = BC.vel_specified[1].z;
                    rho[index(x, y, z)] = rho[index(x - 1, y, z)];    //! rho 由相邻流场内网格点(x-1, y, z)插值而来
                    for (int i = 0; i < direction_size; i++) {
                        double feq = calculate_feq(x, y, z, i);
                        double feq1 = calculate_feq(x - 1, y, z, i);
                        pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x - 1, y, z, i)] - feq1);
                        //! 非平衡部分(pdf2[index(x-1, y, z, i)] - feq1) 由相邻流场内网格点(x-1, y, z)插值而来
                    }
                }
            }
        }
    }

    if (y_id == 0)
    {
        Astring = BC.yminFace;    //! 字符数组赋值给string，为字符串比较准备
        //! 以下给定这个面上的速度，可以是任意方向（包含垂直于该面的进口速度，也可以是滑动壁面（平行该面的速度））
        if (Astring == "velocity")
        {
            int y = 0;
            for (int x = 0; x < NX; x++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    vel[index(x, y, z)].x = BC.vel_specified[2].x;
                    vel[index(x, y, z)].y = BC.vel_specified[2].y;
                    vel[index(x, y, z)].z = BC.vel_specified[2].z;
                    rho[index(x, y, z)] = rho[index(x, y + 1, z)];    //! rho 由相邻流场内网格点(x, y+1, z)插值而来
                    for (int i = 0; i < direction_size; i++)
                    {
                        double feq = calculate_feq(x, y, z, i);
                        double feq1 = calculate_feq(x, y + 1, z, i);
                        pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x, y + 1, z, i)] - feq1);
                        //! 非平衡部分(pdf2[index(x, y+1, z, i)] - feq1) 由相邻流场内网格点(x, y+1, z)插值而来
                    }
                }
            }
        }
    }

    if (y_id == y_np - 1)
    {
        Astring = BC.ymaxFace;
        if (Astring == "velocity")
        {
            int y = NY - 1;
            for (int x = 0; x < NX; x++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    vel[index(x, y, z)].x = BC.vel_specified[3].x;
                    vel[index(x, y, z)].y = BC.vel_specified[3].y;
                    vel[index(x, y, z)].z = BC.vel_specified[3].z;
                    rho[index(x, y, z)] = rho[index(x, y - 1, z)];    //! rho 由相邻流场内网格点(x, y-1, z)插值而来
                    for (int i = 0; i < direction_size; i++)
                    {
                        double feq = calculate_feq(x, y, z, i);
                        double feq1 = calculate_feq(x, y - 1, z, i);
                        pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x, y - 1, z, i)] - feq1);
                        //! 非平衡部分(pdf2[index(x, y+1, z, i)] - feq1) 由相邻流场内网格点(x, y+1, z)插值而来
                    }
                }
            }
        }
    }
    if (direction_size != 9)
    {
        if (z_id == 0)
        {
            Astring = BC.zminFace;    //! 字符数组赋值给string，为字符串比较准备
            if (Astring == "velocity")
            {
                int z = 0;
                for (int x = 0; x < NX; x++)
                {
                    for (int y = 0; y < NY; y++)
                    {
                        vel[index(x, y, z)].x = BC.vel_specified[4].x;
                        vel[index(x, y, z)].y = BC.vel_specified[4].y;
                        vel[index(x, y, z)].z = BC.vel_specified[4].z;
                        rho[index(x, y, z)] = rho[index(x, y, z + 1)];    //! rho 由相邻流场内网格点(x, y, z+1)插值而来
                        for (int i = 0; i < direction_size; i++)
                        {
                            double feq = calculate_feq(x, y, z, i);
                            double feq1 = calculate_feq(x, y, z + 1, i);
                            pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x, y, z + 1, i)] - feq1);
                            //! 非平衡部分(pdf2[index(x, y, z+1, i)] - feq1) 由相邻流场内网格点(x, y, z+1)插值而来
                        }
                    }
                }
            }
        }

        if (z_id == z_np - 1)
        {
            Astring = BC.zmaxFace;    //! 字符数组赋值给string，为字符串比较准备
            if (Astring == "velocity")
            {
                int z = NZ - 1;
                for (int x = 0; x < NX; x++)
                {
                    for (int y = 0; y < NY; y++)
                    {
                        vel[index(x, y, z)].x = BC.vel_specified[5].x;
                        vel[index(x, y, z)].y = BC.vel_specified[5].y;
                        vel[index(x, y, z)].z = BC.vel_specified[5].z;
                        rho[index(x, y, z)] = rho[index(x, y, z - 1)];    //! rho 由相邻流场内网格点(x, y, z-1)插值而来
                        for (int i = 0; i < direction_size; i++)
                        {
                            double feq = calculate_feq(x, y, z, i);
                            double feq1 = calculate_feq(x, y, z - 1, i);
                            pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x, y, z - 1, i)] - feq1);
                            //! 非平衡部分(pdf2[index(x, y, z-1, i)] - feq1) 由相邻流场内网格点(x, y, z-1)插值而来
                        }
                    }
                }
            }
        }
    }
}

void LBM::pressure_boundary_condition()
{
    const double tauinv = 1.0 / tau;
    const double omtauinv = 1.0 - tauinv;
    std::string Astring;
    if (x_id == 0)
    {
        Astring = BC.xminFace;    //! 字符数组赋值给string，为字符串比较准备
        //! 以下给定面上的压强（密度）
        if (Astring == "pressure")
        {
            int x = 0;
            for (int y = 0; y < NY; y++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    vel[index(x, y, z)].x = vel[index(x + 1, y, z)].x;    //! 速度由相邻流场内网格点(x+1, y, z)插值而来
                    vel[index(x, y, z)].y = vel[index(x + 1, y, z)].y;
                    vel[index(x, y, z)].z = vel[index(x + 1, y, z)].z;
                    rho[index(x, y, z)] = BC.rho_specified[0];
                    for (int i = 0; i < direction_size; i++)
                    {
                        double feq = calculate_feq(x, y, z, i);
                        double feq1 = calculate_feq(x + 1, y, z, i);
                        pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x + 1, y, z, i)] - feq1);
                        //! 非平衡部分(pdf2[index(x+1, y, z, i)] - feq1) 由相邻流场内网格点(x+1, y, z)插值而来
                    }
                }
            }
        }
    }

    if (x_id == x_np - 1)
    {
        Astring = BC.xmaxFace;
        if (Astring == "pressure")
        {
            int x = NX - 1;
            for (int y = 0; y < NY; y++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    vel[index(x, y, z)].x = vel[index(x - 1, y, z)].x;
                    vel[index(x, y, z)].y = vel[index(x - 1, y, z)].y;
                    vel[index(x, y, z)].z = vel[index(x - 1, y, z)].z;
                    rho[index(x, y, z)] = BC.rho_specified[1];
                    for (int i = 0; i < direction_size; i++)
                    {
                        double feq = calculate_feq(x, y, z, i);
                        double feq1 = calculate_feq(x - 1, y, z, i);
                        pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x - 1, y, z, i)] - feq1);
                        //! 非平衡部分(pdf2[index(x-1, y, z, i)] - feq1) 由相邻流场内网格点(x-1, y, z)插值而来
                    }
                }
            }
        }
    }

    if (y_id == 0)
    {
        Astring = BC.yminFace;    //! 字符数组赋值给string，为字符串比较准备
        if (Astring == "pressure")
        {
            int y = 0;
            for (int x = 0; x < NX; x++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    vel[index(x, y, z)].x = vel[index(x, y + 1, z)].x;
                    vel[index(x, y, z)].y = vel[index(x, y + 1, z)].y;
                    vel[index(x, y, z)].z = vel[index(x, y + 1, z)].z;
                    rho[index(x, y, z)] = BC.rho_specified[2];
                    for (int i = 0; i < direction_size; i++)
                    {
                        double feq = calculate_feq(x, y, z, i);
                        double feq1 = calculate_feq(x, y + 1, z, i);
                        pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x, y + 1, z, i)] - feq1);
                        //! 非平衡部分(pdf2[index(x, y+1, z, i)] - feq1) 由相邻流场内网格点(x, y+1, z)插值而来
                    }
                }
            }
        }
    }

    if (y_id == y_np - 1)
    {
        Astring = BC.ymaxFace;
        if (Astring == "pressure")
        {
            int y = NY - 1;
            for (int x = 0; x < NX; x++)
            {
                for (int z = 0; z < NZ; z++)
                {
                    vel[index(x, y, z)].x = vel[index(x, y - 1, z)].x;
                    vel[index(x, y, z)].y = vel[index(x, y - 1, z)].y;
                    vel[index(x, y, z)].z = vel[index(x, y - 1, z)].z;
                    rho[index(x, y, z)] = BC.rho_specified[3];
                    for (int i = 0; i < direction_size; i++)
                    {
                        double feq = calculate_feq(x, y, z, i);
                        double feq1 = calculate_feq(x, y - 1, z, i);
                        pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x, y - 1, z, i)] - feq1);
                        //! 非平衡部分(pdf2[index(x, y+1, z, i)] - feq1) 由相邻流场内网格点(x, y+1, z)插值而来
                    }
                }
            }
        }
    }

    if (z_id == 0)
    {
        Astring = BC.zminFace;    //! 字符数组赋值给string，为字符串比较准备
        if (Astring == "pressure")
        {
            int z = 0;
            for (int x = 0; x < NX; x++)
            {
                for (int y = 0; y < NY; y++)
                {
                    vel[index(x, y, z)].x = vel[index(x, y, z + 1)].x;
                    vel[index(x, y, z)].y = vel[index(x, y, z + 1)].y;
                    vel[index(x, y, z)].z = vel[index(x, y, z + 1)].z;
                    rho[index(x, y, z)] = BC.rho_specified[4];
                    for (int i = 0; i < direction_size; i++)
                    {
                        double feq = calculate_feq(x, y, z, i);
                        double feq1 = calculate_feq(x, y, z + 1, i);
                        pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x, y, z + 1, i)] - feq1);
                        //! 非平衡部分(pdf2[index(x, y, z+1, i)] - feq1) 由相邻流场内网格点(x, y, z+1)插值而来
                    }
                }
            }
        }
    }

    if (z_id == z_np - 1)
    {
        Astring = BC.zmaxFace;    //! 字符数组赋值给string，为字符串比较准备
        if (Astring == "pressure")
        {
            int z = NZ - 1;
            for (int x = 0; x < NX; x++)
            {
                for (int y = 0; y < NY; y++)
                {
                    vel[index(x, y, z)].x = vel[index(x, y, z - 1)].x;    //! 由相邻流场内网格点(x, y, z - 1)插值而来
                    vel[index(x, y, z)].y = vel[index(x, y, z - 1)].y;
                    vel[index(x, y, z)].z = vel[index(x, y, z - 1)].z;
                    rho[index(x, y, z)] = BC.rho_specified[5];
                    for (int i = 0; i < direction_size; i++)
                    {
                        double feq = calculate_feq(x, y, z, i);
                        double feq1 = calculate_feq(x, y, z - 1, i);
                        pdf[index(x, y, z, i)] = feq + omtauinv * (pdf2[index(x, y, z - 1, i)] - feq1);
                        //! 非平衡部分(pdf2[index(x, y, z-1, i)] - feq1) 由相邻流场内网格点(x, y, z-1)插值而来
                    }
                }
            }
        }
    }
}

void  LBM::bouzidi()    //! Bouzidi 插值反弹，二阶空间精度
{   //! post-stream (after stream is performed)
    force[0].x = 0.0; force[0].y = 0.0; force[0].z = 0.0;
    for (int i = 0; i < num_node2; i++)
    {
        int j = node2[i];
        int z = (int)floor((double)(j) / (double)((NX+4) * (NY+4)));    //! relocate
        int y = (int)floor((double)(j % ((NX+4) * (NY+4))) / (double)(NX+4));
        int x = (j) % (NX+4);
        //#pragma omp parallel for
        for (int l = 0; l < direction_size; l++)
        {
            if (x >= 0 && x < NX && y >= 0 && y < NY && z >= 0 && z < NZ)
            {
                if (delta[l + i * direction_size] >= 0.00001)
                {
                    double q = delta[l + i * direction_size];
                    int ip = ei[l].x, jp = ei[l].y, kp = ei[l].z;
                    int m = reverse_indexes[l];

                    if (q < 0.5)
                    {
                        int n1 = index(x + 1 * ip, y + 1 * jp, z + 1 * kp);
                        int n2 = index(x + 2 * ip, y + 2 * jp, z + 2 * kp);

                        if (obst[n2] == 0)
                        {
                            pdf[index(x + ip, y + jp, z + kp, l)] = q * (1. + 2. * q) * pdf[index(x, y, z, m)]
                                + (1. - 4. * q * q) * pdf[index(x + ip, y + jp, z + kp, m)] - q * (1. - 2. * q) *
                                pdf[index(x + 2 * ip, y + 2 * jp, z + 2 * kp, m)];
                        }
                        else
                        {   //! the following is seldem used 以下一般用不上
                            if (obst[n1] == 0)
                            {
                                pdf[index(x + ip, y + jp, z + kp, l)] = 2. * q * pdf[index(x, y, z, m)]
                                    + (1. - 2. * q) * pdf[index(x + ip, y + jp, z + kp, m)];
                            }
                            else
                            {
                                pdf[index(x + ip, y + jp, z + kp, l)] = pdf[index(x, y, z, m)];
                            }
                        }
                    }
                    else
                    {
                        int n1 = index(x + 1 * ip, y + 1 * jp, z + 1 * kp);
                        //int n2 = index(x + 2 * ip, y + 2 * jp, z + 2 * kp);
                        //if (obst[n2] == 0) {
                        //    pdf[index(x + ip, y + jp, z + kp, l)] = 1.0 / (q * (2. * q + 1.0)) * pdf[index(x, y, z, m)]
                        //        + (2. * q - 1.) / q * pdf2[index(x + 1 * ip, y + 1 * jp, z + 1 * kp, l)] - (2. * q - 1.) /
                        //        (2. * q + 1.) * pdf2[index(x + 2 * ip, y + 2 * jp, z + 2 * kp, l)];
                        //}
                        //else //the following is seldem used 以下一般用不上
                        //{
                        if (obst[n1] == 0)
                        {
                            pdf[index(x + ip, y + jp, z + kp, l)] = 1. / (2. * q) * pdf[index(x, y, z, m)] +
                                (2. * q - 1.) / (2. * q) * pdf[index(x + 2 * ip, y + 2 * jp, z + 2 * kp, l)];
                        }
                        /*}*/
                    }
                    force[0].x = force[0].x + ei[m].x * (pdf[index(x, y, z, m)] + pdf[index(x + ip, y + jp, z + kp, l)]);
                    force[0].y = force[0].y + ei[m].y * (pdf[index(x, y, z, m)] + pdf[index(x + ip, y + jp, z + kp, l)]);
                    force[0].z = force[0].z + ei[m].z * (pdf[index(x, y, z, m)] + pdf[index(x + ip, y + jp, z + kp, l)]);
                }
            }
        }
    }
}