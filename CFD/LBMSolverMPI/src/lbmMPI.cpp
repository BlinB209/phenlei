#include "LBMSolverMPI.hpp"
#include "vector3.cpp"
#include <mpi.h>

using namespace std;

void LBM::mpi_ini(int mpi_rank)
{
    n_myid = new int[6];
    NX = total_NX / x_np;
    NY = total_NY / y_np;
    NZ = total_NZ / z_np;
    z_id = int(mpi_rank / (x_np * y_np));
    y_id = int((mpi_rank % (x_np * y_np)) / x_np);
    x_id = (mpi_rank % (x_np * y_np)) % x_np;

    n_myid[0] = mpi_rank + 1;    //! x+1节点
    if (x_id == (x_np - 1))
    {
        n_myid[0] = mpi_rank - x_np + 1;
    }

    n_myid[1] = mpi_rank - 1;    //! x-1节点
    if (x_id == (0))
    {
        n_myid[1] = mpi_rank + x_np - 1;
    }

    n_myid[2] = mpi_rank + x_np;    //! y+1节点
    if (y_id == (y_np - 1))
    {
        n_myid[2] = mpi_rank - x_np * (y_np - 1);
    }

    n_myid[3] = mpi_rank - x_np;    //! y-1节点
    if (y_id == (0))
    {
        n_myid[3] = mpi_rank + x_np * (y_np - 1);
    }

    n_myid[4] = mpi_rank + y_np * x_np;    //! z+1节点
    if (z_id == (z_np - 1))
    {
        n_myid[4] = mpi_rank - x_np * y_np * (z_np - 1);
    }

    n_myid[5] = mpi_rank - y_np * x_np;    //! z-1节点
    if (z_id == (0))
    {
        n_myid[5] = mpi_rank + x_np * y_np * (z_np - 1);
    }
    Xori = x_id * NX;
    Yori = y_id * NY;
    Zori = z_id * NZ;

    if (mpi_rank == (x_np * y_np * z_np - 1))
    {
        int NumNodes = total_NX * total_NY * total_NZ;
        rhoTotal = new double[NumNodes] {};
        obstTotal = new int[NumNodes] {};
        velTotal = new vector3<double>[NumNodes] {};
        if (CONV == 1)
        {
            velTotalConv = new vector3<double>[NumNodes] {};
        }
    }
}

void LBM::mpi_send()
{
    int num;
    MPI_Request r1, r2;
    MPI_Status  status;
    //! MPI传值 x方向
    int send_numx = NZ * NY * direction_size;
    double *send_fx = new double[send_numx] {};
    double *recv_fx = new double[send_numx] {};

    int x = NX - 1;
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int y = 0; y < NY; y++)
        {
            for (int i = 0; i < direction_size; i++)
            {
                send_fx[num] = pdf2[index(x, y, z, i)];
                num++;
            }
        }
    }

    MPI_Isend(send_fx, send_numx, MPI_DOUBLE_PRECISION, n_myid[0], 200, MPI_COMM_WORLD, &r1);
    MPI_Irecv(recv_fx, send_numx, MPI_DOUBLE_PRECISION, n_myid[1], 200, MPI_COMM_WORLD, &r2);

    MPI_Wait(&r1, &status);
    MPI_Wait(&r2, &status);

    x = -1;
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int y = 0; y < NY; y++)
        {
            for (int i = 0; i < direction_size; i++)
            {
                pdf2[index(x, y, z, i)] = recv_fx[num];
                num++;
            }
        }
    }

    x = 0;
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int y = 0; y < NY; y++) {
            for (int i = 0; i < direction_size; i++)
            {
                send_fx[num] = pdf2[index(x, y, z, i)];
                num++;
            }
        }
    }

    MPI_Isend(send_fx, send_numx, MPI_DOUBLE_PRECISION, n_myid[1], 201, MPI_COMM_WORLD, &r1);
    MPI_Irecv(recv_fx, send_numx, MPI_DOUBLE_PRECISION, n_myid[0], 201, MPI_COMM_WORLD, &r2);

    MPI_Wait(&r1, &status);
    MPI_Wait(&r2, &status);

    x = NX;
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int y = 0; y < NY; y++)
        {
            for (int i = 0; i < direction_size; i++)
            {
                pdf2[index(x, y, z, i)] = recv_fx[num];
                num++;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    x = NX - 2;
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int y = 0; y < NY; y++)
        {
            for (int i = 0; i < direction_size; i++)
            {
                send_fx[num] = pdf2[index(x, y, z, i)];
                num++;
            }
        }
    }

    MPI_Isend(send_fx, send_numx, MPI_DOUBLE_PRECISION, n_myid[0], 200, MPI_COMM_WORLD, &r1);
    MPI_Irecv(recv_fx, send_numx, MPI_DOUBLE_PRECISION, n_myid[1], 200, MPI_COMM_WORLD, &r2);

    MPI_Wait(&r1, &status);
    MPI_Wait(&r2, &status);

    x = -2;
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int y = 0; y < NY; y++)
        {
            for (int i = 0; i < direction_size; i++)
            {
                pdf2[index(x, y, z, i)] = recv_fx[num];
                num++;
            }
        }
    }

    x = 1;
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int y = 0; y < NY; y++)
        {
            for (int i = 0; i < direction_size; i++)
            {
                send_fx[num] = pdf2[index(x, y, z, i)];
                num++;
            }
        }
    }

    MPI_Isend(send_fx, send_numx, MPI_DOUBLE_PRECISION, n_myid[1], 201, MPI_COMM_WORLD, &r1);
    MPI_Irecv(recv_fx, send_numx, MPI_DOUBLE_PRECISION, n_myid[0], 201, MPI_COMM_WORLD, &r2);

    MPI_Wait(&r1, &status);
    MPI_Wait(&r2, &status);

    x = NX + 1;
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int y = 0; y < NY; y++)
        {
            for (int i = 0; i < direction_size; i++)
            {
                pdf2[index(x, y, z, i)] = recv_fx[num];
                num++;
            }
        }
    }

    delete [] send_fx;    send_fx = nullptr;
    delete [] recv_fx;    recv_fx = nullptr;
    MPI_Barrier(MPI_COMM_WORLD);

    //! MPI传值 y方向
    int send_numy = (NX + 4) * NZ * direction_size;
    double *send_fy = new double[send_numy] {};
    double *recv_fy = new double[send_numy] {};

    int y = NY - 1;
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int x = -2; x < NX + 2; x++)
        {
            for (int i = 0; i < direction_size; i++)
            {
                send_fy[num] = pdf2[index(x, y, z, i)];
                num++;
            }
        }
    }

    MPI_Isend(send_fy, send_numy, MPI_DOUBLE_PRECISION, n_myid[2], 202, MPI_COMM_WORLD, &r1);
    MPI_Irecv(recv_fy, send_numy, MPI_DOUBLE_PRECISION, n_myid[3], 202, MPI_COMM_WORLD, &r2);

    MPI_Wait(&r1, &status);
    MPI_Wait(&r2, &status);

    y = -1;
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int x = -2; x < NX + 2; x++)
        {
            for (int i = 0; i < direction_size; i++)
            {
                pdf2[index(x, y, z, i)] = recv_fy[num];
                num++;
            }
        }
    }

    y = 0;
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int x = -2; x < NX + 2; x++)
        {
            for (int i = 0; i < direction_size; i++)
            {
                send_fy[num] = pdf2[index(x, y, z, i)];
                num++;
            }
        }
    }

    MPI_Isend(send_fy, send_numy, MPI_DOUBLE_PRECISION, n_myid[3], 203, MPI_COMM_WORLD, &r1);
    MPI_Irecv(recv_fy, send_numy, MPI_DOUBLE_PRECISION, n_myid[2], 203, MPI_COMM_WORLD, &r2);

    MPI_Wait(&r1, &status);
    MPI_Wait(&r2, &status);

    y = NY;
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int x = -2; x < NX + 2; x++)
        {
            for (int i = 0; i < direction_size; i++)
            {
                pdf2[index(x, y, z, i)] = recv_fy[num];
                num++;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    y = NY - 2;
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int x = -2; x < NX + 2; x++)
        {
            for (int i = 0; i < direction_size; i++)
            {
                send_fy[num] = pdf2[index(x, y, z, i)];
                num++;
            }
        }
    }

    MPI_Isend(send_fy, send_numy, MPI_DOUBLE_PRECISION, n_myid[2], 202, MPI_COMM_WORLD, &r1);
    MPI_Irecv(recv_fy, send_numy, MPI_DOUBLE_PRECISION, n_myid[3], 202, MPI_COMM_WORLD, &r2);

    MPI_Wait(&r1, &status);
    MPI_Wait(&r2, &status);

    y = -2;
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int x = -2; x < NX + 2; x++)
        {
            for (int i = 0; i < direction_size; i++)
            {
                pdf2[index(x, y, z, i)] = recv_fy[num];
                num++;
            }
        }
    }

    y = 1;
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int x = -2; x < NX + 2; x++)
        {
            for (int i = 0; i < direction_size; i++)
            {
                send_fy[num] = pdf2[index(x, y, z, i)];
                num++;
            }
        }
    }

    MPI_Isend(send_fy, send_numy, MPI_DOUBLE_PRECISION, n_myid[3], 203, MPI_COMM_WORLD, &r1);
    MPI_Irecv(recv_fy, send_numy, MPI_DOUBLE_PRECISION, n_myid[2], 203, MPI_COMM_WORLD, &r2);

    MPI_Wait(&r1, &status);
    MPI_Wait(&r2, &status);

    y = NY + 1;
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int x = -2; x < NX + 2; x++)
        {
            for (int i = 0; i < direction_size; i++)
            {
                pdf2[index(x, y, z, i)] = recv_fy[num];
                num++;
            }
        }
    }
    delete [] send_fy;    send_fy = nullptr;
    delete [] recv_fy;    recv_fy = nullptr;
    MPI_Barrier(MPI_COMM_WORLD);

    //! MPI传值 z方向

    if (direction_size != 9)
    {
        int send_numz = (NX + 4) * (NY + 4) * direction_size;
        double* send_fz = new double[send_numz] {};
        double* recv_fz = new double[send_numz] {};

        int z = NZ - 1;
        num = 0;
        for (int y = -2; y < NY + 2; y++)
        {
            for (int x = -2; x < NX + 2; x++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                    send_fz[num] = pdf2[index(x, y, z, i)];
                    num++;
                }
            }
        }

        MPI_Isend(send_fz, send_numz, MPI_DOUBLE_PRECISION, n_myid[4], 204, MPI_COMM_WORLD, &r1);
        MPI_Irecv(recv_fz, send_numz, MPI_DOUBLE_PRECISION, n_myid[5], 204, MPI_COMM_WORLD, &r2);

        MPI_Wait(&r1, &status);
        MPI_Wait(&r2, &status);

        z = -1;
        num = 0;
        for (int y = -2; y < NY + 2; y++)
        {
            for (int x = -2; x < NX + 2; x++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                    pdf2[index(x, y, z, i)] = recv_fz[num];
                    num++;
                }
            }
        }

        z = 0;
        num = 0;
        for (int y = -2; y < NY + 2; y++)
        {
            for (int x = -2; x < NX + 2; x++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                    send_fz[num] = pdf2[index(x, y, z, i)];
                    num++;
                }
            }
        }

        MPI_Isend(send_fz, send_numz, MPI_DOUBLE_PRECISION, n_myid[5], 205, MPI_COMM_WORLD, &r1);
        MPI_Irecv(recv_fz, send_numz, MPI_DOUBLE_PRECISION, n_myid[4], 205, MPI_COMM_WORLD, &r2);

        MPI_Wait(&r1, &status);
        MPI_Wait(&r2, &status);

        z = NZ;
        num = 0;
        for (int y = -2; y < NY + 2; y++)
        {
            for (int x = -2; x < NX + 2; x++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                    pdf2[index(x, y, z, i)] = recv_fz[num];
                    num++;
                }
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        z = NZ - 2;
        num = 0;
        for (int y = -2; y < NY + 2; y++)
        {
            for (int x = -2; x < NX + 2; x++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                    send_fz[num] = pdf2[index(x, y, z, i)];
                    num++;
                }
            }
        }

        MPI_Isend(send_fz, send_numz, MPI_DOUBLE_PRECISION, n_myid[4], 204, MPI_COMM_WORLD, &r1);
        MPI_Irecv(recv_fz, send_numz, MPI_DOUBLE_PRECISION, n_myid[5], 204, MPI_COMM_WORLD, &r2);

        MPI_Wait(&r1, &status);
        MPI_Wait(&r2, &status);

        z = -2;
        num = 0;
        for (int y = -2; y < NY + 2; y++)
        {
            for (int x = -2; x < NX + 2; x++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                    pdf2[index(x, y, z, i)] = recv_fz[num];
                    num++;
                }
            }
        }

        z = 1;
        num = 0;
        for (int y = -2; y < NY + 2; y++)
        {
            for (int x = -2; x < NX + 2; x++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                    send_fz[num] = pdf2[index(x, y, z, i)];
                    num++;
                }
            }
        }

        MPI_Isend(send_fz, send_numz, MPI_DOUBLE_PRECISION, n_myid[5], 205, MPI_COMM_WORLD, &r1);
        MPI_Irecv(recv_fz, send_numz, MPI_DOUBLE_PRECISION, n_myid[4], 205, MPI_COMM_WORLD, &r2);

        MPI_Wait(&r1, &status);
        MPI_Wait(&r2, &status);

        z = NZ + 1;
        num = 0;
        for (int y = -2; y < NY + 2; y++)
        {
            for (int x = -2; x < NX + 2; x++)
            {
                for (int i = 0; i < direction_size; i++)
                {
                    pdf2[index(x, y, z, i)] = recv_fz[num];
                    num++;
                }
            }
        }
        delete [] send_fz;    send_fz = nullptr;
        delete [] recv_fz;    recv_fz = nullptr;
        MPI_Barrier(MPI_COMM_WORLD);
    }

}

void LBM::mpi_total(int mpi_rank)
{
    int send_num = NX * NZ * NY * 5 + 6;
    int np = x_np * y_np * z_np;
    int recv_num;
    int num;
    double *send_all{}, *recv_all{};
    int *recv_Count{}, *displs{};

    if (mpi_rank == (np - 1))
    {
        recv_Count = new int[np] {};
        displs = new int[np] {};
    }

    MPI_Gather(&send_num, 1, MPI_INT, recv_Count, 1, MPI_INT, np - 1, MPI_COMM_WORLD);

    if (mpi_rank == (np - 1))
    {
        displs[0] = 0;
        for (int i = 0; i < np - 1; i++)
        {
            displs[i + 1] = displs[i] + recv_Count[i];
        }
        recv_num = displs[np - 1] + recv_Count[np - 1];
    }

    MPI_Barrier(MPI_COMM_WORLD);

    send_all = new double[send_num] {};
    num = 0;
    for (int z = 0; z < NZ; z++)
    {
        for (int y = 0; y < NY; y++)
        {
            for (int x = 0; x < NX; x++)
            {
                send_all[num] = vel[index(x, y, z)].x;
                send_all[num + NX * NZ * NY] = vel[index(x, y, z)].y;
                send_all[num + NX * NZ * NY * 2] = vel[index(x, y, z)].z;
                send_all[num + NX * NZ * NY * 3] = rho[index(x, y, z)];
                send_all[num + NX * NZ * NY * 4] = obst[index(x, y, z)];
                num = num + 1;
            }
        }
        send_all[NX * NZ * NY * 5 + 0] = static_cast<double>(Xori);
        send_all[NX * NZ * NY * 5 + 1] = static_cast<double>(Yori);
        send_all[NX * NZ * NY * 5 + 2] = static_cast<double>(Zori);
        send_all[NX * NZ * NY * 5 + 3] = static_cast<double>(NX);
        send_all[NX * NZ * NY * 5 + 4] = static_cast<double>(NY);
        send_all[NX * NZ * NY * 5 + 5] = static_cast<double>(NZ);
    }

    if (mpi_rank == (np - 1))
    {
        recv_all = new double[recv_num] {};
    }

    MPI_Gatherv(send_all, send_num, MPI_DOUBLE, recv_all, recv_Count, displs, MPI_DOUBLE, np - 1, MPI_COMM_WORLD);

    if (mpi_rank == (np - 1))
    {
        for (int i = 0; i < np; i++)
        {
            int num_n = displs[i];
            int num_m = recv_Count[i];
            int x0, y0, z0;
            int NX1, NY1, NZ1;
            x0 = static_cast<int>(recv_all[num_n + num_m - 6]);
            y0 = static_cast<int>(recv_all[num_n + num_m - 5]);
            z0 = static_cast<int>(recv_all[num_n + num_m - 4]);
            NX1 = static_cast<int>(recv_all[num_n + num_m - 3]);
            NY1 = static_cast<int>(recv_all[num_n + num_m - 2]);
            NZ1 = static_cast<int>(recv_all[num_n + num_m - 1]);
            num = 0;
            for (int z = 0; z < NZ1; z++)
            {
                for (int y = 0; y < NY1; y++)
                {
                    for (int x = 0; x < NX1; x++)
                    {
                        velTotal[total_index(x + x0, y + y0, z + z0)].x = recv_all[num_n + num];
                        velTotal[total_index(x + x0, y + y0, z + z0)].y = recv_all[num_n + num + NX1 * NZ1 * NY1];
                        velTotal[total_index(x + x0, y + y0, z + z0)].z = recv_all[num_n + num + NX1 * NZ1 * NY1 * 2];
                        rhoTotal[total_index(x + x0, y + y0, z + z0)] = recv_all[num_n + num + NX1 * NZ1 * NY1 * 3];
                        obstTotal[total_index(x + x0, y + y0, z + z0)] = static_cast<int>(recv_all[num_n + num + NX1 * NZ1 * NY1 * 4]);
                        num = num + 1;
                    }
                }
            }
        }
        delete [] recv_all;     recv_all = nullptr;
        delete [] recv_Count;    recv_Count = nullptr;
        delete [] displs;    displs = nullptr;
        if (CONV == 1)
        {
            check_converge();
        }
    }
    delete [] send_all;    send_all = nullptr;
    MPI_Barrier(MPI_COMM_WORLD);
}
