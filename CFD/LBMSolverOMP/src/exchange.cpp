#include "LBMSolverOMP.hpp"
#include "vector3.cpp"
#include <iomanip>
#include <cmath>
//#include <mpi.h>

void LBM::init_border(int *Low, int Fdimx, int Fdimy, int Fdimz, int &num, int *border, int *border_dir) // in solver (C)
{   //! function: 初始化记录被密网格盖住的最边缘粗网格点
    int dirs =  send_direction_size;
    int m ;

    int Zmin = Low[2];
    int Zmax = Low[2] + Fdimz - 1;    //! general 3D
    if (dirs == 3)    //!specifically for D2Q9 ! possible bug
    {
        Zmin = -1; Zmax = 1;
    }

    for (int x = 0; x < NX; x++)
    {
        for (int y = 0; y < NY; y++)
        {
            for (int z = 0; z < NZ; z++)
            {
                lab[index(x, y, z)] = 0;
                if (x > Low[0] &&  x<Low[0] + Fdimx-1 && y>Low[1] &&  y < Low[1] + Fdimy-1
                    && z>Zmin &&  z < Zmax)
                {
                    lab[index(x, y, z)] = 1;
                }
            }
        }
    }
    std::cout << "Lower Left corner("<<Low[0] << ", " << Low[1] << " ," << Low[2]  <<")"<< '\n';
    std::cout << Fdimx << " " << Fdimy << " " << Fdimz << '\n';
    m = 0;    //! initialize m
    for (int x = 0; x < NX; x++)
    {
        for (int y = 0; y < NY; y++)
        {
            for (int z = 0; z < NZ; z++)
            {
                int count = 0;    //! initialize direction count
                for (int k = 0; k < direction_size; k++)
                {
                    int xmd = (NX + x + ei[k].x) % NX;
                    int ymd = (NY + y + ei[k].y) % NY;
                    int zmd = (NZ + z + ei[k].z) % NZ;
                    if (lab[index(x, y, z)] == 0 &&
                        lab[index(xmd, ymd, zmd)] == 1)
                    {
                        border[m] = index(x, y, z);    //! record border node index记录被密网格盖住的最边缘粗网格点
                        border_dir[m*dirs + count] = k;    //! record associated directions记录有需求的方向
                        //std::cout << x <<" " <<y << " " <<z << " " <<k  << '\n';
                        count++;
                    }
                }
                if (count > 0) m++;
            }
        }
    }
    num = m ;    //! actual border_node_number
}

void LBM::exchange(int id, int border_node_num, int *border, int *border_dir, double *f_border_all_dir, double *f_border_dir, int *Low, int *NN)
{
    int dirs = send_direction_size;
    double temp = 0.0;

    if (id == 1)
    {   //! fetch data from coarse mesh, perform inside solver
        for (int m = 0; m < border_node_num; m++)
        {
            /*            for (int k = 0; k < dirs; k++) {
                int i = border_dir[m*dirs + k];
                if ( i!= 0) {
                    f_border_dir[dirs*m + k] = pdf2[border[m] + i * NX * NY * NZ];
                }
            }*/
           for (int k = 0; k < direction_size; k++)
           {
                f_border_all_dir[direction_size*m + k] = pdf2[border[m] + k * NX * NY * NZ];
           }
        }
    }

    if (id == 3)
    {   //! trasfer data from Coarse to Fine , perform inside solver2 (fine mesh),
        for (int m = 0; m < border_node_num; m++)
        {
            int z = (int)floor((double)(border[m]) / (double)((NN[0])*(NN[1])));    //! relocate grobal location in coarse mesh
            int y = (int)floor((double)((border[m]) % ((NN[0])*(NN[1]))) / (double)(NN[0]));
            int x = (border[m]) % (NN[0]);
            int x2 = 2 * (x - Low[0]);
            int y2 = 2 * (y - Low[1]);
            int z2 = 2 * (z - Low[2]);
            int z2p = z2 + 1;
            if (dirs == 3) z2p = z2;    //! D2Q9 model
/*        for (int k = 0; k < dirs; k++) {
                int i = border_dir[m*dirs + k];
                if(i!=0){
                {            
                temp = f_border_dir[m*dirs + k]; 
            pdf2[index(x2,     y2,    z2, i)] = temp;
            pdf2[index(x2 + 1, y2,    z2, i)] = temp;
            pdf2[index(x2,    y2 + 1, z2, i)] = temp;
            pdf2[index(x2,    y2,     z2p, i)] = temp;
            pdf2[index(x2 + 1, y2 + 1, z2, i)] = temp;
            pdf2[index(x2 + 1, y2,    z2p, i)] = temp;
            pdf2[index(x2,    y2 + 1, z2p, i)] = temp;
            pdf2[index(x2 + 1,y2 + 1, z2p, i)] = temp;
        }
    }
}*/
            for (int k = 0; k < direction_size; k++)
            {
                 temp = f_border_all_dir[m*direction_size + k];
                 pdf2[index(x2,  y2,  z2,   k)] = temp;
                 pdf2[index(x2+1,y2,  z2,   k)] = temp;
                 pdf2[index(x2,  y2+1,z2,   k)] = temp;
                 pdf2[index(x2,  y2,  z2p,  k)] = temp;
                 pdf2[index(x2+1,y2+1,z2,   k)] = temp;
                 pdf2[index(x2+1, y2, z2p,  k)] = temp;
                 pdf2[index(x2,  y2+1,z2p,  k)] = temp;
                 pdf2[index(x2+1,y2+1, z2p, k)] = temp;
            }
        }
    }

    if (id == 2)
    {   //! fetch data from fine mesh，perform inside fine mesh
        for (int m = 0; m < border_node_num; m++)
        {
            int z = (int)floor((double)(border[m]) / (double)((NN[0])*(NN[1])));    //! relocate grobal location in coarse mesh
            int y = (int)floor((double)((border[m]) % ((NN[0])*(NN[1]))) / (double)(NN[0]));
            int x = (border[m]) % (NN[0]);

            for (int k = 0; k < dirs; k++)
            {
                int re_dir = reverse_indexes[border_dir[dirs*m + k]];    //! 反方向

                if (re_dir != 0)
                {
                    int x2 = 2 * (x - Low[0]);
                    int y2 = 2 * (y - Low[1]);
                    int z2 = 2 * (z - Low[2]);
                    int z2p = z2 + 1;
                    if (dirs == 3) z2p = z2;    //! D2Q9 model
                    temp = pdf[index(x2, y2, z2, re_dir)] + pdf[index(x2 + 1, y2, z2, re_dir)] + pdf[index(x2, y2 + 1, z2, re_dir)];
                    temp += pdf[index(x2, y2, z2p, re_dir)] + pdf[index(x2 + 1, y2 + 1, z2, re_dir)] + pdf[index(x2, y2 + 1, z2p, re_dir)];
                    temp += pdf[index(x2 + 1, y2, z2p, re_dir)] + pdf[index(x2 + 1, y2 + 1, z2p, re_dir)];
                    temp = temp / 8.0;
                    f_border_dir[dirs*m + k] = temp;
                }
            }
        }
    }

    if (id == 4)
    {
        //! trandfer data from fine to coarse, perform inside Coarse mesh
        for (int m = 0; m < border_node_num; m++)
        {
            for (int k = 0; k < dirs; k++)
            {
                int re_dir = reverse_indexes[border_dir[dirs*m + k]];    //! 反方向
                if (re_dir != 0) { pdf[border[m] + re_dir * NX * NY * NZ] = f_border_dir[dirs*m + k]; }
            }
        }
    }
}

inline int indexa(int x, int y, int z, int NX, int NY)
{
    return (z * NX * NY) + (y * NX) + x;
}

void LBM::copy_data(int id, vector3<double> *vels,double *rhos, int *Low, int *NN)
{
    if (id == 1)
    {   //! copy data from fine to coarse
        if (NN[2] == 1)    //! 二维，NZ2=0
        {
            for (int x = 1; x < NN[0] / 2 - 1; x++)
            {
                for (int y = 1; y < NN[1] / 2 - 1; y++)
                {
                        vel[index(Low[0] + x, Low[1] + y, Low[2])] = vels[indexa(x, y, 0, NN[0] / 2, NN[1] / 2)];
                        rho[index(Low[0] + x, Low[1] + y, Low[2])] = rhos[indexa(x, y, 0, NN[0] / 2, NN[1] / 2)];
                }
            }
        }
        else
        {
            for (int x = 1; x < NN[0] / 2 - 1; x++)
            {
                for (int y = 1; y < NN[1] / 2 - 1; y++)
                {
                    for (int z = 1; z < NN[2] / 2 - 1; z++)
                    {
                        vel[index(Low[0] + x, Low[1] + y, Low[2] + z)] = vels[indexa(x, y, z, NN[0] / 2, NN[1] / 2)];
                        rho[index(Low[0] + x, Low[1] + y, Low[2] + z)] = rhos[indexa(x, y, z, NN[0] / 2, NN[1] / 2)];
                    }
                }
            }
        }
    }
    else
    {   //! id==2    //get data from fine mesh
        for (int x = 0; x < NN[0]; x=x+2)
        {
            for (int y = 0; y < NN[1]; y=y+2)
            {
                for (int z = 0; z < NN[2]; z=z+2)
                {
                    int zp = (z + 1)%NN[2];
                    vels[indexa(x / 2, y / 2, z / 2, NN[0] / 2, NN[1] / 2)] = (vel[index(x, y, z)] + vel[index(x + 1, y, z)]
                    +vel[index(x, y+1, z)] + vel[index(x, y, zp)]+vel[index(x+1, y+1, z)] + vel[index(x +1, y, zp)]
                    + vel[index(x, y+1, zp)] + vel[index(x + 1, y+1, zp)])/8.0;

                    rhos[indexa(x /2, y /2, z /2, NN[0]/2,NN[1]/2)] = (rho[index(x, y, z)] + rho[index(x + 1, y, z)]
                    + rho[index(x, y + 1, z)] + rho[index(x, y, zp)] + rho[index(x + 1, y + 1, z)] + rho[index(x + 1, y, zp)]
                    + rho[index(x, y + 1, zp)] + rho[index(x + 1, y + 1, zp)]) / 8.0;
                }
            }
        }
    }
}