#include <streambuf>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <mpi.h>
#include "LBMSolverMPI.hpp"
#include "vector3.cpp"
//! 作者:  黄海波，张嘉骞，中国科学技术大学
using namespace std;

void RunLBMSolverMPI()
{   //! int argc, char** argv
    //! Initialize the MPI environment
    //MPI_Init(NULL, NULL);

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    RUNPHLBM(mpi_rank);
    // Finalize the MPI environment.
    MPI_Finalize();
}

void RUNPHLBM(int mpi_rank)
{
    LBM *solver = new LBM();
    struct Params_Pakage parameters = { "./bin/input.txt","./params_out.txt" };

    solver->initialise(&parameters, mpi_rank);

    int runs = solver->get_total_steps();
    for (int i = 0; i < runs; i = i + 1)
    {
        solver->perform_timestep();
        solver->output(i, mpi_rank);
    }
    delete solver;    solver = nullptr;
}