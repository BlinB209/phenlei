#include "PHMpi.h"
#include "Solver.h"

namespace PHSPACE
{
vector< vector< PHSolver *> *> * GlobalSolvers::globalSolversList = new vector< vector< PHSolver *> * >;

PHSolver *GlobalSolvers::GetSolver(int zoneID, int solverID)
{
    return (*(*globalSolversList)[zoneID])[solverID];
}

int GlobalSolvers::GetNumberOfSolvers(int iZone)
{
    return static_cast<int>((*globalSolversList)[iZone]->size());
}

void GlobalSolvers::AddSolverToGlobal(int iZone, vector < PHSolver * > *solversList)
{
    if (globalSolversList->size() == 0)
    {
        int nZones = PHMPI::GetNumberofGlobalZones();
        globalSolversList->resize(nZones, 0);
    }
    (*globalSolversList)[iZone] = solversList;
}

}
