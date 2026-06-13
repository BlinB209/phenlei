#include "DemoSolver.h"
#include "Param_DemoSolver.h"
#include "TK_Log.h"

namespace PHSPACE
{

DemoSolver::DemoSolver()
{

}

DemoSolver::~DemoSolver()
{

}

LIB_EXPORT void DemoSolver::InitControlParameters()
{
    FreeControlParameters();
    controlParameters = new Param_DemoSolver();
    controlParameters->Init();
}

LIB_EXPORT Param_DemoSolver *DemoSolver::GetControlParameters() const
{
    return static_cast <Param_DemoSolver *> (controlParameters);
}

void DemoSolver::AllocateGlobalVar(Grid * grid)
{
    int nTotalNode = grid->GetNTotalNode();

    if (grid->GetLevel() == 0)
    {
        int *nodeBCType = new int[nTotalNode];
        PHSPACE::SetField(nodeBCType, 0, nTotalNode);
        grid->UpdateDataPtr("nodeBCType", nodeBCType);
    }
}

void DemoSolver::DeAllocateGlobalVar(Grid * grid)
{
    PrintToWindow("DemoSolver::DeAllocateGlobalVar\n");
}

void DemoSolver::ZeroResiduals(Grid * grid)
{
    PrintToWindow("DemoSolver::ZeroResiduals\n");
}

void DemoSolver::TimeStep(Grid * grid)
{
    PrintToWindow("DemoSolver::TimeStep\n");
}

void DemoSolver::ComputeResidual(Grid *grid_in, FieldProxy *rhs_proxy)
{
    PrintToWindow("DemoSolver::ComputeResidual\n");
}

void DemoSolver::SourceFlux(Grid *grid)
{
    PrintToWindow("DemoSolver::SourceFlux\n");
}

void DemoSolver::InviscidFlux(Grid *grid)
{
    PrintToWindow("DemoSolver::InviscidFlux\n");
}

void DemoSolver::ViscousFlux(Grid *grid)
{
    PrintToWindow("DemoSolver::ViscousFlux\n");
}

void DemoSolver::RungeKutta(Grid *grid, FieldProxy *rhs_proxy)
{
    PrintToWindow("DemoSolver::RungeKutta\n");
}

void DemoSolver::Post()
{
    PrintToWindow("DemoSolver::Post\n");
}

void DemoSolver::PostSolve(Grid * grid, int stage, int level)
{
    PrintToWindow("DemoSolver::PostSolve\n");
}

}