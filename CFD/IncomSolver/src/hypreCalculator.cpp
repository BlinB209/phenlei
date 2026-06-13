#include <iostream>
#include <iomanip>
#include "hypreCalculator.h"
#include "Geo_UnstructGrid.h"
#include "TK_Time.h"
#include "TK_Exit.h"

namespace PHSPACE
{
hypreCalculator::hypreCalculator() : LinearEquationSystemCalculator()
{
#ifndef USE_HypreLib
  TK_Exit::ExceptionExit("hypre lib is not opened in CMAKEFile!");
#endif
}

hypreCalculator::~hypreCalculator()
{
    mathLibSolveFinalize();
}

void hypreCalculator::mathLibSolveInit(Grid *gridIn)
{
#ifdef USE_HypreLib
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int zoneID = grid->GetZoneID();

    int *globalCellNoOfGlobalZones = reinterpret_cast<int *>(GlobalDataBase::GetDataPtr("globalCellNoOfGlobalZones"));
    int *maxIterOfEachEq = reinterpret_cast<int *>(GlobalDataBase::GetDataPtr("maxIterOfEachEq"));
    RDouble *iterTolOfEachEq = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("iterTolOfEachEq"));

    localSize = nTotalCell;
    ilower = globalCellNoOfGlobalZones[zoneID];
    iupper = globalCellNoOfGlobalZones[zoneID + 1] - 1;
    rows = new int [localSize];
    for (int i = 0; i < localSize; i++)
    {
        rows[i] = ilower + i;
    }

    colValueOfOneRow = new RDouble [100];
    colNumberOfOneRow = new int [100];

    //HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A_hypre);
    //HYPRE_IJMatrixSetObjectType(A_hypre, HYPRE_PARCSR);
    //HYPRE_IJMatrixInitialize(A_hypre);

    ///* Create the rhs and solution */
    //HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b_hypre);
    //HYPRE_IJVectorSetObjectType(b_hypre, HYPRE_PARCSR);
    //HYPRE_IJVectorInitialize(b_hypre);

    //HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x_hypre);
    //HYPRE_IJVectorSetObjectType(x_hypre, HYPRE_PARCSR);
    //HYPRE_IJVectorInitialize(x_hypre);

    int linearSolverSize = 50;
    int *precondID = reinterpret_cast <int *>(GlobalDataBase::GetDataPtr("precondMethodOfLinEqSystem"));
    precondList = new HYPRE_Solver * [linearSolverSize];
    for (int iSolver = 0 ; iSolver < linearSolverSize; iSolver++)
    {
        if (maxIterOfEachEq[iSolver] > 0)
        {
            precondList[iSolver] = new HYPRE_Solver();
        }
    }
    HYPRE_PtrToSolverFcn *SolvePtr = new HYPRE_PtrToSolverFcn[linearSolverSize];
    HYPRE_PtrToSolverFcn *SetupPtr = new HYPRE_PtrToSolverFcn[linearSolverSize];
    for (int iSolver = 0 ; iSolver < linearSolverSize; iSolver++)
    {
        SolvePtr[iSolver] = NULL;
        SetupPtr[iSolver] = NULL;
    }

    for (int iSolver = 0 ; iSolver < linearSolverSize; iSolver++)
    {
        if (maxIterOfEachEq[iSolver] < 0)
        {
            continue;
        }

        int precondType = precondID[iSolver];
        HYPRE_Solver *precond = precondList[iSolver];
        if (precondType == iNone)
        {

        }
        else if (precondType == iILU)
        {
            HYPRE_ILUCreate(&(*precond));
            HYPRE_ILUSetTol(*precond, 0.0);
            HYPRE_ILUSetMaxIter(*precond, 1);
            HYPRE_ILUSetLogging(*precond, 0);
            
            SolvePtr[iSolver] = (HYPRE_PtrToSolverFcn) HYPRE_ILUSolve;
            SetupPtr[iSolver] = (HYPRE_PtrToSolverFcn) HYPRE_ILUSetup;
        }
        else if (precondType == iAMG)
        {
            HYPRE_BoomerAMGCreate(&(*precond));
            HYPRE_BoomerAMGSetTol            ((*precond), 0.0); // Only the number of iteration is important
            HYPRE_BoomerAMGSetMaxIter        ((*precond), 1);   // Do only one iteration
            HYPRE_BoomerAMGSetLogging(*precond, 0);

            SolvePtr[iSolver] = (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve;
            SetupPtr[iSolver] = (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup;
        }
        else if (precondType == iAMGDD)
        {
            HYPRE_BoomerAMGDDCreate(&(*precond));
            HYPRE_Solver amg_solver;
            HYPRE_BoomerAMGDDGetAMG((*precond), &amg_solver);

            HYPRE_Int amgdd_start_level = 0;
            HYPRE_Int amgdd_padding = 1;
            HYPRE_Int amgdd_fac_num_relax = 1;
            HYPRE_Int amgdd_num_comp_cycles = 2;
            HYPRE_Int amgdd_fac_relax_type = 3;
            HYPRE_Int amgdd_fac_cycle_type = 1;
            HYPRE_Int amgdd_num_ghost_layers = 1;

            HYPRE_BoomerAMGDDSetStartLevel((*precond), amgdd_start_level);
            HYPRE_BoomerAMGDDSetPadding((*precond), amgdd_padding);
            HYPRE_BoomerAMGDDSetFACNumRelax((*precond), amgdd_fac_num_relax);
            HYPRE_BoomerAMGDDSetFACNumCycles((*precond), amgdd_num_comp_cycles);
            HYPRE_BoomerAMGDDSetFACRelaxType((*precond), amgdd_fac_relax_type);
            HYPRE_BoomerAMGDDSetFACCycleType((*precond), amgdd_fac_cycle_type);
            HYPRE_BoomerAMGDDSetNumGhostLayers((*precond), amgdd_num_ghost_layers);

            HYPRE_BoomerAMGSetTol            (amg_solver, 0.0); // Only the number of iteration is important
            HYPRE_BoomerAMGSetMaxIter        (amg_solver, 1);   // Do only one iteration
            HYPRE_BoomerAMGSetLogging(amg_solver, 0);

            SolvePtr[iSolver] = (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGDDSolve;
            SetupPtr[iSolver] = (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGDDSetup;
        }
        else if (precondType == iEuclid)
        {
            HYPRE_EuclidCreate(MPI_COMM_WORLD, &(*precond));

            // Hardcoded parameters
            //HYPRE_EuclidSetSparseA ((*precond), 0.1);
            //HYPRE_EuclidSetRowScale((*precond), 1);
            SolvePtr[iSolver] = (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve;
            SetupPtr[iSolver] = (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup;
        }
        else if (precondType == DiagScale)
        {
            *precond = NULL;
            SolvePtr[iSolver] = (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScale;
            SetupPtr[iSolver] = (HYPRE_PtrToSolverFcn) HYPRE_ParCSRDiagScaleSetup;
        }
        else if (precondType == PILUT)
        {
            HYPRE_ParCSRPilutCreate( MPI_COMM_WORLD, &(*precond) );

            SolvePtr[iSolver] = (HYPRE_PtrToSolverFcn) HYPRE_ParCSRPilutSolve;
            SetupPtr[iSolver] = (HYPRE_PtrToSolverFcn) HYPRE_ParCSRPilutSetup;
        }
        else if (precondType == GSMG)
        {
            HYPRE_BoomerAMGCreate(&(*precond));
            HYPRE_BoomerAMGSetGSMG((*precond), 4);

            SolvePtr[iSolver] = (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve;
            SetupPtr[iSolver] = (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup;
        }
        else if (precondType == ParaSails)
        {
            HYPRE_ParaSailsCreate( MPI_COMM_WORLD, &(*precond) );
            SolvePtr[iSolver] = (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve;
            SetupPtr[iSolver] = (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup;
        }
    }

    int *solMethodOfLinEqSystem = reinterpret_cast <int *>(GlobalDataBase::GetDataPtr("solMethodOfLinEqSystem"));
    solverList = new HYPRE_Solver * [linearSolverSize];
    for (int iSolver = 0 ; iSolver < linearSolverSize; iSolver++)
    {
        if (maxIterOfEachEq[iSolver] < 0)
        {
            continue;
        }
        solverList[iSolver] = new HYPRE_Solver();
    }

    for (int iSolver = 0 ; iSolver < linearSolverSize; iSolver++)
    {
        if (maxIterOfEachEq[iSolver] < 0)
        {
            continue;
        }

        int linearSolver = solMethodOfLinEqSystem[iSolver];
        HYPRE_Solver *solver = solverList[iSolver];
        int maxIterStep = maxIterOfEachEq[iSolver];
        RDouble iterTol = iterTolOfEachEq[iSolver];
        
        if (linearSolver == GMRESMethod)
        {
            HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, solver);
            HYPRE_GMRESSetKDim((*solver), 10);
            HYPRE_GMRESSetMaxIter((*solver), maxIterStep);
            HYPRE_GMRESSetTol((*solver), iterTol);
            HYPRE_GMRESSetPrintLevel((*solver), 0);
            HYPRE_GMRESSetLogging((*solver), 0);
            if (precondID[iSolver] > 0)
            {
                HYPRE_GMRESSetPrecond((*solver), SolvePtr[iSolver], 
                            SetupPtr[iSolver], *(precondList[iSolver]));
            }
        }
        else if (linearSolver == BicgstabMethod)
        {
            HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &(*solver));
            HYPRE_ParCSRBiCGSTABSetMaxIter((*solver), maxIterStep);
            HYPRE_ParCSRBiCGSTABSetTol((*solver), iterTol);
            HYPRE_ParCSRBiCGSTABSetPrintLevel((*solver), 0);
            HYPRE_ParCSRBiCGSTABSetLogging((*solver), 0);

            if (precondID[iSolver] > 0)
            {
                HYPRE_BiCGSTABSetPrecond((*solver), SolvePtr[iSolver], 
                                        SetupPtr[iSolver], *(precondList[iSolver]));
            }
        }
        else if (linearSolver == AMGMethod)
        {
            HYPRE_BoomerAMGCreate(&(*solver));
            HYPRE_BoomerAMGSetTol((*solver), iterTol);
            HYPRE_BoomerAMGSetMaxIter((*solver), maxIterStep);
            HYPRE_BoomerAMGSetPrintLevel((*solver), 0);
            HYPRE_BoomerAMGSetLogging((*solver), 0);

            HYPRE_BoomerAMGSetAggNumLevels((*solver), 1);
            HYPRE_BoomerAMGSetRelaxType((*solver), 3);
            HYPRE_BoomerAMGSetCoarsenType((*solver), 10);
            HYPRE_BoomerAMGSetInterpType((*solver), 3);
        }
        else if (linearSolver == PCGMethod)
        {
            HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &(*solver));
            HYPRE_PCGSetMaxIter((*solver), maxIterStep); /* max iterations */
            HYPRE_PCGSetTol((*solver), iterTol); /* conv. tolerance */
            HYPRE_PCGSetTwoNorm((*solver), 1); /* use the two norm as the stopping criteria */
            HYPRE_PCGSetPrintLevel((*solver), 0); /* prints out the iteration info */
            HYPRE_PCGSetLogging((*solver), 0); /* needed to get run info later */
            if (precondID[iSolver] > 0)
            {
                HYPRE_PCGSetPrecond((*solver), SolvePtr[iSolver], 
                                    SetupPtr[iSolver], *(precondList[iSolver]));
            }
        }
        else if (linearSolver == FlexGMRESMethod)
        {
            HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, solver);
            HYPRE_FlexGMRESSetKDim((*solver), 5);
            HYPRE_FlexGMRESSetMaxIter((*solver), maxIterStep);
            HYPRE_FlexGMRESSetTol((*solver), iterTol);
            HYPRE_FlexGMRESSetPrintLevel((*solver), 0);
            HYPRE_FlexGMRESSetLogging((*solver), 0);
            if (precondID[iSolver] > 0)
            {
                HYPRE_FlexGMRESSetPrecond((*solver), SolvePtr[iSolver], 
                    SetupPtr[iSolver], *(precondList[iSolver]));
            }
        }
        else if (linearSolver == LGMRESMethod)
        {

            HYPRE_ParCSRLGMRESCreate(MPI_COMM_WORLD, solver);

            HYPRE_ParCSRLGMRESSetTol((*solver), iterTol);
            HYPRE_ParCSRLGMRESSetMaxIter((*solver), maxIterStep);

            if (precondID[iSolver] > 0)
            {
                HYPRE_LGMRESSetPrecond((*solver), SolvePtr[iSolver], 
                    SetupPtr[iSolver], *(precondList[iSolver]));
            }
        }
        else if (linearSolver == COGMRESMethod)
        {

            HYPRE_ParCSRCOGMRESCreate(MPI_COMM_WORLD, solver);
            HYPRE_ParCSRCOGMRESSetTol((*solver), iterTol);
            HYPRE_ParCSRCOGMRESSetMaxIter((*solver), maxIterStep);

            //if (precondID[iSolver] > 0)
            //{
            //    HYPRE_ParCSRCOGMRESSetPrecond((*solver), SolvePtr[iSolver], 
            //        SetupPtr[iSolver], *(precondList[iSolver]));
            //}
        }
        else if (linearSolver == CGNRMethod)
        {
            //HYPRE_ParCSRCGNRCreate(MPI_COMM_WORLD, solver);
            //HYPRE_ParCSRCGNRSetLogging(*solver,0);
            //if (precondID[iSolver] > 0)
            //{
            //    HYPRE_CGNRSetPrecond((*solver), SolvePtr[iSolver], SolvePtr[iSolver], 
            //        SetupPtr[iSolver], *(precondList[iSolver]));
            //}
        }
        else if (linearSolver == MGRMethod)
        {
            //! not complete.
            HYPRE_MGRCreate(solver);

            if (precondID[iSolver] > 0)
            {
                HYPRE_CGNRSetPrecond((*solver), SolvePtr[iSolver], SolvePtr[iSolver], 
                    SetupPtr[iSolver], *(precondList[iSolver]));
            }
        }
        else if (linearSolver == ILUMethod)
        {
            //! not complete.
            HYPRE_ILUCreate(solver);
            HYPRE_ILUSetMaxIter(*solver, maxIterStep );
            HYPRE_ILUSetTol(*solver, iterTol );
        }
        else if (linearSolver == AMGDDMethod)
        {
            HYPRE_BoomerAMGDDCreate(solver);
        }
        else if (linearSolver == HybridMethod)
        {
            //! not complete.
            HYPRE_ParCSRHybridCreate(solver);
            HYPRE_ParCSRHybridSetTol(*solver, iterTol);
            if (precondID[iSolver] > 0)
            {
                //HYPRE_ParCSRHybridSetPrecond((*solver), SolvePtr[iSolver], SolvePtr[iSolver], 
                //    SetupPtr[iSolver], *(precondList[iSolver]));
            }
        }
        else
        {
            HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &(*solver));
            HYPRE_GMRESSetKDim((*solver), 5);
            HYPRE_GMRESSetMaxIter((*solver), 100);
            HYPRE_GMRESSetTol((*solver), iterTol);
            HYPRE_GMRESSetPrintLevel((*solver), 0);
            HYPRE_GMRESSetLogging((*solver), 0);
        }
    }

    RDouble GlobalTotalCells = GlobalDataBase::GetDoubleParaFromDB("GlobalTotalCells");
    int *localCell2GlobalMap = grid->GetLocalCell2GlobalMap();
    UnstructBCSet **bcr = grid->GetBCRecord();

    ColVectors = new vector<int>[localSize];
    for (int iCell = 0; iCell < localSize; iCell++)
    {
        int globalCellNumber = localCell2GlobalMap[iCell];
        ColVectors[iCell].push_back(globalCellNumber);
    }

    for (int iFace = 0; iFace < nTotalFace; iFace++)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        if (re < localSize)
        {
            ColVectors[le].push_back(localCell2GlobalMap[re]);
            ColVectors[re].push_back(localCell2GlobalMap[le]);
        }
        else
        {
            int bcType = bcr[iFace]->GetKey();
            if (bcType != PHENGLEI::INTERFACE)
            {
                continue;
            }
            if (localCell2GlobalMap[re] < GlobalTotalCells)
            {
                if (localCell2GlobalMap[re] > localCell2GlobalMap[le])
                {
                    ColVectors[le].push_back(localCell2GlobalMap[re]);
                }
                else if (localCell2GlobalMap[re] < localCell2GlobalMap[le])
                {
                    ColVectors[le].push_back(localCell2GlobalMap[re]);
                }
            }
        }    
    }

    delete [] SolvePtr;
    delete [] SetupPtr;
#endif
}

void hypreCalculator::mathLibSolve(Grid *gridIn, int iEquation, RDouble *fi, RDouble *diagCoef, RDouble *b, RDouble *upNeighCoef, RDouble *dnNeighCoef)
{
#ifdef USE_HypreLib
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    
    int *localCell2GlobalMap = grid->GetLocalCell2GlobalMap();
    RDouble GlobalTotalCells = GlobalDataBase::GetDoubleParaFromDB("GlobalTotalCells");
    
    //int *constVarNumToSolveVarNum = reinterpret_cast<int *>(grid->GetDataPtr("constVarNumToSolveVarNum"));
    HYPRE_Solver *solver = solverList[iEquation];
    int *linearSolverID = reinterpret_cast <int *>(GlobalDataBase::GetDataPtr("solMethodOfLinEqSystem"));
    int linearSolver = linearSolverID[iEquation];
    
    int *cell2ZoneIDMap = grid->GetCell2ZoneIDMap();
    
    UnstructBCSet **bcr = grid->GetBCRecord();

    int *maxIterOfEachEq = reinterpret_cast<int *>(GlobalDataBase::GetDataPtr("maxIterOfEachEq"));
    RDouble *iterTolOfEachEq = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("iterTolOfEachEq"));

    int maxIterStep = maxIterOfEachEq[iEquation];
    RDouble iterTol = iterTolOfEachEq[iEquation];
    //int *ActualIterStepNum = reinterpret_cast <int * > (GlobalDataBase::GetDataPtr("ActualIterStepNum"));
    //RDouble *ActualIterTol = reinterpret_cast <RDouble * > (GlobalDataBase::GetDataPtr("ActualIterTol"));
    //RDouble *ActualSolveTime = reinterpret_cast <RDouble * > (GlobalDataBase::GetDataPtr("ActualSolveTime"));
    ////RDouble timeStartSolve = PHSPACE::TIME_SPACE::GetWallTime();

    int num_iterations;
    RDouble residual;

    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A_hypre);
    HYPRE_IJMatrixSetObjectType(A_hypre, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A_hypre);
    
    /* Create the rhs and solution */
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b_hypre);
    HYPRE_IJVectorSetObjectType(b_hypre, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b_hypre);
    
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x_hypre);
    HYPRE_IJVectorSetObjectType(x_hypre, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(x_hypre);
    
    
    ValuesVector = new vector<double>[localSize];
    for (int iCell = 0; iCell < localSize; iCell++)
    {
        ValuesVector[iCell].push_back(diagCoef[iCell]);
    }
    for (int iFace = 0; iFace < nTotalFace; iFace++)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        if (re < localSize)
        {
            ValuesVector[le].push_back(upNeighCoef[iFace]);
            ValuesVector[re].push_back(dnNeighCoef[iFace]);
        }
        else
        {
            int bcType = bcr[iFace]->GetKey();
            if (bcType != PHENGLEI::INTERFACE)
            {
                continue;
            }
            if (localCell2GlobalMap[re] < GlobalTotalCells)
            {
                if (localCell2GlobalMap[re] > localCell2GlobalMap[le])
                {
                    ValuesVector[le].push_back(upNeighCoef[iFace]);
                }
                else if (localCell2GlobalMap[re] < localCell2GlobalMap[le])
                {
                    ValuesVector[le].push_back(upNeighCoef[iFace]);
                }
            }
        }    
    }
    
    for (int i = 0; i < localSize; i++)
    {
        int colSize = ColVectors[i].size();
        int globalRowNumber = rows[i];
        for (int jCell = 0; jCell < ColVectors[i].size(); jCell++)
        {
            colNumberOfOneRow[jCell] = ColVectors[i][jCell];
        }
    
        for (int jCell = 0; jCell < ColVectors[i].size(); jCell++)
        {
            colValueOfOneRow[jCell] = ValuesVector[i][jCell];
        }
    
        HYPRE_IJMatrixSetValues(A_hypre, 1, &colSize, &globalRowNumber, colNumberOfOneRow, colValueOfOneRow);
    }
    
    HYPRE_IJMatrixAssemble(A_hypre);
    HYPRE_IJMatrixGetObject(A_hypre, (void **)&parcsr_A_hypre);  
    
    HYPRE_IJVectorSetValues(b_hypre, localSize, rows, b);
    HYPRE_IJVectorAssemble(b_hypre);
    HYPRE_IJVectorGetObject(b_hypre, (void **)&par_b_hypre);
    
    HYPRE_IJVectorSetValues(x_hypre, localSize, rows, fi);
    HYPRE_IJVectorAssemble(x_hypre);
    HYPRE_IJVectorGetObject(x_hypre, (void **)&par_x_hypre);
    
    int *precondID = reinterpret_cast <int *>(GlobalDataBase::GetDataPtr("precondID"));
    HYPRE_Solver *precond = precondList[iEquation];

    //if (precondID[iEquation] == iAMGDD)
    //{
    //    destoryTempDataInAMGDD(*precond);
    //}

    RDouble timeStartSolve = PHSPACE::TIME_SPACE::GetWallTime();
    if (linearSolver == GaussSeidelMethod) 
    {
        explicitSolve(gridIn, fi, diagCoef, b, upNeighCoef, dnNeighCoef);
    }
    else if (linearSolver == GMRESMethod) 
    {
        HYPRE_ParCSRGMRESSetup((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_ParCSRGMRESSolve((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_IJVectorGetValues(x_hypre, localSize, rows, fi);

        HYPRE_ParCSRGMRESGetNumIterations((*solver), &num_iterations);
        HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm((*solver), &residual);
    }
    else if (linearSolver == BicgstabMethod)
    {
        HYPRE_ParCSRBiCGSTABSetup((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_ParCSRBiCGSTABSolve((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_IJVectorGetValues(x_hypre, localSize, rows, fi);
        
        HYPRE_ParCSRBiCGSTABGetNumIterations((*solver), &num_iterations);
        HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm((*solver), &residual);
    }
    else if (linearSolver == AMGMethod)
    {
        HYPRE_BoomerAMGSetup((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_BoomerAMGSolve((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_IJVectorGetValues(x_hypre, localSize, rows, fi);
        HYPRE_BoomerAMGGetNumIterations((*solver), &num_iterations);
        HYPRE_BoomerAMGGetFinalRelativeResidualNorm((*solver), &residual);
    }
    else if (linearSolver == PCGMethod)
    {
        HYPRE_ParCSRPCGSetup((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_ParCSRPCGSolve((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_IJVectorGetValues(x_hypre, localSize, rows, fi);
        HYPRE_ParCSRPCGGetNumIterations((*solver), &num_iterations);
        HYPRE_ParCSRPCGGetFinalRelativeResidualNorm((*solver), &residual);
    }
    else if (linearSolver == FlexGMRESMethod) 
    {
        HYPRE_ParCSRFlexGMRESSetup((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_ParCSRFlexGMRESSolve((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_IJVectorGetValues(x_hypre, localSize, rows, fi);
        HYPRE_ParCSRFlexGMRESGetNumIterations((*solver), &num_iterations);
        HYPRE_ParCSRFlexGMRESGetFinalRelativeResidualNorm((*solver), &residual);
    }
    else if (linearSolver == LGMRESMethod)
    {
        HYPRE_ParCSRLGMRESSetup((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_ParCSRLGMRESSolve((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_IJVectorGetValues(x_hypre, localSize, rows, fi);
        HYPRE_ParCSRLGMRESGetNumIterations((*solver), &num_iterations);
        HYPRE_ParCSRLGMRESGetFinalRelativeResidualNorm((*solver), &residual);
    }
    else if (linearSolver == COGMRESMethod)
    {
        HYPRE_ParCSRCOGMRESSetup((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_ParCSRCOGMRESSolve((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_IJVectorGetValues(x_hypre, localSize, rows, fi);
        HYPRE_ParCSRCOGMRESGetNumIterations((*solver), &num_iterations);
        HYPRE_ParCSRCOGMRESGetFinalRelativeResidualNorm((*solver), &residual);
    }
    else if (linearSolver == CGNRMethod)
    {
        HYPRE_ParCSRCGNRCreate(MPI_COMM_WORLD, solver);
        HYPRE_ParCSRCGNRSetup((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_ParCSRCGNRSolve((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_IJVectorGetValues(x_hypre, localSize, rows, fi);
        HYPRE_ParCSRCGNRDestroy((*solver));
        HYPRE_ParCSRCGNRGetNumIterations((*solver), &num_iterations);
        HYPRE_ParCSRCGNRGetFinalRelativeResidualNorm((*solver), &residual);
    }
    else if (linearSolver == MGRMethod)
    {
        HYPRE_MGRSetup((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_MGRSolve((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_IJVectorGetValues(x_hypre, localSize, rows, fi);
        HYPRE_MGRGetNumIterations((*solver), &num_iterations);
        HYPRE_MGRGetFinalRelativeResidualNorm((*solver), &residual);
    }
    else if (linearSolver == ILUMethod)
    {
        HYPRE_ILUSetup((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_ILUSolve((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_IJVectorGetValues(x_hypre, localSize, rows, fi);
        HYPRE_ILUGetNumIterations((*solver), &num_iterations);
        HYPRE_ILUGetFinalRelativeResidualNorm((*solver), &residual);
    }
    else if (linearSolver == AMGDDMethod)
    {
        
        destoryTempDataInAMGDD((*solver));

        HYPRE_BoomerAMGDDSetup((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_BoomerAMGDDSolve((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_IJVectorGetValues(x_hypre, localSize, rows, fi);
        
        HYPRE_BoomerAMGDDGetNumIterations((*solver), &num_iterations);
        HYPRE_BoomerAMGDDGetFinalRelativeResidualNorm((*solver), &residual);
        
    }
    else if (linearSolver == HybridMethod)
    {
        HYPRE_ParCSRHybridSetup((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_ParCSRHybridSolve((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_IJVectorGetValues(x_hypre, localSize, rows, fi);
        HYPRE_ParCSRHybridGetNumIterations((*solver), &num_iterations);
        HYPRE_ParCSRHybridGetFinalRelativeResidualNorm((*solver), &residual);
    }
    else 
    {
        HYPRE_ParCSRGMRESSetup((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_ParCSRGMRESSolve((*solver), parcsr_A_hypre, par_b_hypre, par_x_hypre);
        HYPRE_IJVectorGetValues(x_hypre, localSize, rows, fi);
        HYPRE_ParCSRGMRESGetNumIterations((*solver), &num_iterations);
        HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm((*solver), &residual);
    }


    if (precondID[iEquation] == iAMGDD)
    {
        //HYPRE_BoomerAMGDDDestroy(AMGDDPrecond);
    }

    RDouble timeEndSolve = PHSPACE::TIME_SPACE::GetWallTime();
    //ActualIterStepNum[iEquation] = num_iterations;
    //ActualIterTol[iEquation] = residual;
    //ActualSolveTime[iEquation] = timeEndSolve - timeStartSolve;

    for (int iCell = 0; iCell < localSize; iCell++)
    {
        ValuesVector[iCell].clear();
    }
    delete [] ValuesVector;  
    HYPRE_IJMatrixDestroy(A_hypre);
    HYPRE_IJVectorDestroy(b_hypre);
    HYPRE_IJVectorDestroy(x_hypre);
#endif
}

void hypreCalculator::mathLibSolveFinalize()
{
#ifdef USE_HypreLib
    delete [] rows;
    for (int iCell = 0; iCell < localSize; iCell++)
    {
        ColVectors[iCell].clear();
    }
    delete [] ColVectors;
    delete [] colValueOfOneRow;
    delete [] colNumberOfOneRow;

    int linearSolverSize = GlobalDataBase::GetIntParaFromDB("linearSolverSize");

    int *linearSolverID = reinterpret_cast <int *>(GlobalDataBase::GetDataPtr("solMethodOfLinEqSystem"));
    for (int iSolver = 0 ; iSolver < linearSolverSize; iSolver++)
    {
        int linearSolver = linearSolverID[iSolver];
        HYPRE_Solver *solver = solverList[iSolver];

        if (linearSolver == GMRESMethod)
        {
            /* clean up */
            HYPRE_ParCSRGMRESDestroy((*solver));
        }
        else if (linearSolver == BicgstabMethod)
        {
            HYPRE_ParCSRBiCGSTABDestroy((*solver));
        }
        else if (linearSolver == AMGMethod)
        {
            /* Destroy solver */
            HYPRE_BoomerAMGDestroy((*solver));
        }
        else if (linearSolver == AMGDDMethod)
        {
            //HYPRE_BoomerAMGDDDestroy((*solver));
        }
        else if (linearSolver == PCGMethod)
        {
            /* Destroy solver */
            HYPRE_ParCSRPCGDestroy((*solver));
        }
        else if (linearSolver == FlexGMRESMethod)
        {
            HYPRE_ParCSRFlexGMRESDestroy((*solver));
        }
        else if (linearSolver == LGMRESMethod)
        {
            HYPRE_ParCSRLGMRESDestroy((*solver));
        }
        else if (linearSolver == COGMRESMethod)
        {
            HYPRE_ParCSRCOGMRESDestroy((*solver));
        }
        else if (linearSolver == CGNRMethod)
        {
            //HYPRE_ParCSRCGNRDestroy((*solver));
        }
        else if (linearSolver == MGRMethod)
        {
            HYPRE_MGRDestroy((*solver));
        }
        else if (linearSolver == ILUMethod)
        {
            HYPRE_ILUDestroy((*solver));
        }
        else
        {
            HYPRE_ParCSRGMRESDestroy((*solver));
        }
    }

    int *precondID = reinterpret_cast <int *>(GlobalDataBase::GetDataPtr("precondID"));
    
    for (int iSolver = 0 ; iSolver < linearSolverSize; iSolver++)
    {
        if (precondID[iSolver] == iAMG)
        {
            HYPRE_BoomerAMGDestroy(*(precondList[iSolver]));
        }
        else if (precondID[iSolver] == iEuclid)
        {
            HYPRE_EuclidDestroy(*(precondList[iSolver]));
        }
        else
        {
            delete precondList[iSolver];
        }
    }
    delete [] precondList;
#endif
}

#ifdef USE_HypreLib
void hypreCalculator::SetAMGDDPrecondParameters(HYPRE_Solver& AMGDDPrecond)
{

    HYPRE_Int amgdd_start_level = 0;
    HYPRE_Int amgdd_padding = 1;
    HYPRE_Int amgdd_fac_num_relax = 1;
    HYPRE_Int amgdd_num_comp_cycles = 2;
    HYPRE_Int amgdd_fac_relax_type = 3;
    HYPRE_Int amgdd_fac_cycle_type = 1;
    HYPRE_Int amgdd_num_ghost_layers = 1;

    HYPRE_BoomerAMGDDSetStartLevel(AMGDDPrecond, amgdd_start_level);
    HYPRE_BoomerAMGDDSetPadding(AMGDDPrecond, amgdd_padding);
    HYPRE_BoomerAMGDDSetFACNumRelax(AMGDDPrecond, amgdd_fac_num_relax);
    HYPRE_BoomerAMGDDSetFACNumCycles(AMGDDPrecond, amgdd_num_comp_cycles);
    HYPRE_BoomerAMGDDSetFACRelaxType(AMGDDPrecond, amgdd_fac_relax_type);
    HYPRE_BoomerAMGDDSetFACCycleType(AMGDDPrecond, amgdd_fac_cycle_type);
    HYPRE_BoomerAMGDDSetNumGhostLayers(AMGDDPrecond, amgdd_num_ghost_layers);

    HYPRE_Solver        amg_precond = NULL;
    HYPRE_BoomerAMGDDGetAMG(AMGDDPrecond, &amg_precond);

    /* for CGC BM Aug 25, 2006 */
    HYPRE_Int      cgcits = 1;
    /* for coordinate plotting BM Oct 24, 2006 */
    HYPRE_Int      plot_grids = 0;
    HYPRE_Int      coord_dim  = 3;
    float         *coordinates = NULL;

    /* parameters for ParaSAILS */
    HYPRE_Real   sai_threshold = 0.1;
    HYPRE_Real   sai_filter = 0.1;

    /* parameters for PILUT */
    HYPRE_Real   drop_tol = -1;
    HYPRE_Int    nonzeros_to_keep = -1;

    /* parameters for Euclid or ILU smoother in AMG */
    HYPRE_Real   eu_ilut = 0.0;
    HYPRE_Real   eu_sparse_A = 0.0;
    HYPRE_Int    eu_bj = 0;
    HYPRE_Int    eu_level = -1;
    HYPRE_Int    eu_stats = 0;
    HYPRE_Int    eu_mem = 0;
    HYPRE_Int    eu_row_scale = 0; /* Euclid only */

                                   /* parameters for GMRES */
                                   //HYPRE_Int    k_dim;
                                   /* parameters for COGMRES */
    HYPRE_Int    cgs = 1;
    HYPRE_Int    unroll = 0;
    /* parameters for LGMRES */
    //HYPRE_Int    aug_dim;
    /* parameters for GSMG */
    HYPRE_Int    gsmg_samples = 5;
    /* interpolation */
    HYPRE_Int    interp_type  = 6; /* default value */
    HYPRE_Int    post_interp_type  = 0; /* default value */
                                        /* RL: restriction */
    HYPRE_Int    restri_type = 0;
    /* aggressive coarsening */
    HYPRE_Int    agg_interp_type  = 4; /* default value */
    HYPRE_Int    agg_P_max_elmts  = 0; /* default value */
    HYPRE_Int    agg_P12_max_elmts  = 0; /* default value */
    HYPRE_Real   agg_trunc_factor  = 0; /* default value */
    HYPRE_Real   agg_P12_trunc_factor  = 0; /* default value */

    HYPRE_Int    print_system = 0;
    HYPRE_Int    rel_change = 0;

    /* begin lobpcg */
    HYPRE_Int    hybrid = 1;
    HYPRE_Int    num_sweep = 1;
    HYPRE_Int    relax_default = 3;

    HYPRE_Int  lobpcgFlag = 0;
    HYPRE_Int  lobpcgGen = 0;
    HYPRE_Int  constrained = 0;
    HYPRE_Int  vFromFileFlag = 0;
    HYPRE_Int  lobpcgSeed = 0;
    HYPRE_Int  blockSize = 1;
    HYPRE_Int  verbosity = 1;
    //HYPRE_Int  iterations;
    //HYPRE_Int  maxIterations = 100;
    HYPRE_Int  checkOrtho = 0;
    HYPRE_Int  printLevel = 0; /* also c.f. poutdat */
    HYPRE_Int  two_norm = 1;
    //HYPRE_Int  pcgIterations = 0;
    HYPRE_Int  pcgMode = 1;
    HYPRE_Real pcgTol = 1e-2;
    //HYPRE_Real nonOrthF;


    /* parameters for BoomerAMG */
    HYPRE_Real     A_drop_tol = 0.0;
    HYPRE_Int      A_drop_type = -1;
    HYPRE_Int      coarsen_cut_factor = 0;
    HYPRE_Real     strong_threshold;
    HYPRE_Real     strong_thresholdR;
    HYPRE_Real     filter_thresholdR;
    HYPRE_Real     trunc_factor;
    HYPRE_Real     jacobi_trunc_threshold;
    HYPRE_Real     S_commpkg_switch = 1.0;
    HYPRE_Real     CR_rate = 0.7;
    HYPRE_Real     CR_strong_th = 0.0;
    HYPRE_Int      CR_use_CG = 0;
    HYPRE_Int      P_max_elmts = 4;
    HYPRE_Int      cycle_type;
    HYPRE_Int      fcycle;
    HYPRE_Int      coarsen_type = 10;
    HYPRE_Int      measure_type = 0;
    HYPRE_Int      num_sweeps = 1;
    HYPRE_Int      IS_type = 1;
    HYPRE_Int      num_CR_relax_steps = 2;
    HYPRE_Int      relax_type = -1;
    HYPRE_Int      add_relax_type = 18;
    HYPRE_Int      relax_coarse = -1;
    HYPRE_Int      relax_up = -1;
    HYPRE_Int      relax_down = -1;
    HYPRE_Int      relax_order = 0;
    HYPRE_Int      level_w = -1;
    HYPRE_Int      level_ow = -1;
    /* HYPRE_Int    smooth_lev; */
    /* HYPRE_Int    smooth_rlx = 8; */
    HYPRE_Int      smooth_type = 6;
    HYPRE_Int      smooth_num_levels = 0;
    HYPRE_Int      smooth_num_sweeps = 1;
    HYPRE_Int      coarse_threshold = 9;
    HYPRE_Int      min_coarse_size = 0;
    /* redundant coarse grid solve */
    HYPRE_Int      seq_threshold = 0;
    HYPRE_Int      redundant = 0;
    /* additive versions */
    HYPRE_Int    additive = -1;
    HYPRE_Int    mult_add = -1;
    HYPRE_Int    simple = -1;
    HYPRE_Int    add_last_lvl = -1;
    HYPRE_Int    add_P_max_elmts = 0;
    HYPRE_Real   add_trunc_factor = 0;
    HYPRE_Int    rap2     = 0;
    HYPRE_Int    mod_rap2 = 0;
    HYPRE_Int    keepTranspose = 0;
#ifdef HYPRE_USING_DSUPERLU
    HYPRE_Int    dslu_threshold = -1;
#endif
    //HYPRE_Real   relax_wt;
    HYPRE_Real   add_relax_wt = 1.0;
    //HYPRE_Real   relax_wt_level;
    //HYPRE_Real   outer_wt;
    //HYPRE_Real   outer_wt_level;
    HYPRE_Real   tol = 1.e-8, pc_tol = 0.;
    HYPRE_Real   atol = 0.0;
    HYPRE_Real   max_row_sum = 1.;
    HYPRE_Int    converge_type = 0;

    HYPRE_Int  cheby_order = 2;
    HYPRE_Int  cheby_eig_est = 10;
    HYPRE_Int  cheby_variant = 0;
    HYPRE_Int  cheby_scale = 1;
    HYPRE_Real cheby_fraction = .3;

    HYPRE_Int           max_levels = 25;
    //HYPRE_Int           num_iterations;
    //HYPRE_Int           pcg_num_its, dscg_num_its;
    //HYPRE_Int           max_iter = 1000;
    //HYPRE_Int           mg_max_iter = 100;
    HYPRE_Int           nodal = 0;
    HYPRE_Int           nodal_diag = 0;
    HYPRE_Int           keep_same_sign = 0;
    HYPRE_Real          cf_tol = 0.9;

    /* coasening data */
    HYPRE_Int     num_cpt = 0;
    HYPRE_Int     num_fpt = 0;
    HYPRE_Int     num_isolated_fpt = 0;
    HYPRE_BigInt *cpt_index = NULL;
    HYPRE_BigInt *fpt_index = NULL;
    HYPRE_BigInt *isolated_fpt_index = NULL;

    HYPRE_BigInt *row_nums = NULL;
    HYPRE_Int *num_cols = NULL;
    HYPRE_BigInt *col_nums = NULL;
    HYPRE_Real *data = NULL;

    HYPRE_Int air = 0;
    HYPRE_Int **grid_relax_points = NULL;

    /* amg-dd options */
    //HYPRE_Int amgdd_start_level = 0;
    //HYPRE_Int amgdd_padding = 1;
    //HYPRE_Int amgdd_fac_num_relax = 1;
    //HYPRE_Int amgdd_num_comp_cycles = 2;
    //HYPRE_Int amgdd_fac_relax_type = 3;
    //HYPRE_Int amgdd_fac_cycle_type = 1;
    //HYPRE_Int amgdd_num_ghost_layers = 1;


    strong_threshold = 0.25;
    strong_thresholdR = 0.25;
    filter_thresholdR = 0.00;
    trunc_factor = 0.;
    jacobi_trunc_threshold = 0.01;
    cycle_type = 1;
    fcycle = 0;
    HYPRE_Real relax_wt = 1.;
    HYPRE_Real outer_wt = 1.;

    HYPRE_Int           build_matrix_type = 2;
    //build_matrix_arg_index = argc;
    HYPRE_Int           build_rhs_type = 2;
    //build_rhs_arg_index = argc;
    HYPRE_Int           build_src_type = -1;
    //build_src_arg_index = argc;
    HYPRE_Int           build_x0_type = -1;
    //build_x0_arg_index = argc;
    HYPRE_Int           build_funcs_type = 0;
    //build_funcs_arg_index = argc;
    HYPRE_Int           build_fpt_arg_index = 0;
    HYPRE_Int           build_sfpt_arg_index = 0;
    HYPRE_Int           build_cpt_arg_index = 0;
    //IS_type = 1;
    HYPRE_Int           debug_flag = 0;
    HYPRE_Int           solver_id = 0;
    HYPRE_Int           ioutdat = 3;
    HYPRE_Int           poutdat = 0;


    HYPRE_Int           num_functions = 1;
    HYPRE_Int           num_paths = 1;
    HYPRE_Int           agg_num_levels = 0;
    HYPRE_Int           ns_coarse = 1, ns_down = -1, ns_up = -1;


    /* defaults for Schwarz */
    HYPRE_Int           variant = 0;  /* multiplicative */
    HYPRE_Int           overlap = 1;  /* 1 layer overlap */
    HYPRE_Int           domain_type = 2; /* through agglomeration */
    HYPRE_Real          schwarz_rlx_weight = 1.;

    HYPRE_Int           use_nonsymm_schwarz = 0;

    /* default execution policy and memory space */
    HYPRE_ExecutionPolicy default_exec_policy = HYPRE_EXEC_DEVICE;
    HYPRE_MemoryLocation  memory_location     = HYPRE_MEMORY_DEVICE;

    HYPRE_BoomerAMGSetCGCIts(amg_precond, cgcits);
    HYPRE_BoomerAMGSetInterpType(amg_precond, interp_type);
    HYPRE_BoomerAMGSetRestriction(amg_precond, restri_type); /* 0: P^T, 1: AIR, 2: AIR-2 */
    HYPRE_BoomerAMGSetPostInterpType(amg_precond, post_interp_type);
    HYPRE_BoomerAMGSetNumSamples(amg_precond, gsmg_samples);
    HYPRE_BoomerAMGSetTol(amg_precond, pc_tol);
    HYPRE_BoomerAMGSetCoarsenType(amg_precond, coarsen_type);
    HYPRE_BoomerAMGSetCoarsenCutFactor(amg_precond, coarsen_cut_factor);
    HYPRE_BoomerAMGSetCPoints(amg_precond, max_levels, num_cpt, cpt_index);
    HYPRE_BoomerAMGSetFPoints(amg_precond, num_fpt, fpt_index);
    HYPRE_BoomerAMGSetIsolatedFPoints(amg_precond, num_isolated_fpt, isolated_fpt_index);
    HYPRE_BoomerAMGSetMeasureType(amg_precond, measure_type);
    HYPRE_BoomerAMGSetStrongThreshold(amg_precond, strong_threshold);
    HYPRE_BoomerAMGSetSeqThreshold(amg_precond, seq_threshold);
    HYPRE_BoomerAMGSetRedundant(amg_precond, redundant);
    HYPRE_BoomerAMGSetMaxCoarseSize(amg_precond, coarse_threshold);
    HYPRE_BoomerAMGSetMinCoarseSize(amg_precond, min_coarse_size);
    HYPRE_BoomerAMGSetTruncFactor(amg_precond, trunc_factor);
    HYPRE_BoomerAMGSetPMaxElmts(amg_precond, P_max_elmts);
    HYPRE_BoomerAMGSetJacobiTruncThreshold(amg_precond, jacobi_trunc_threshold);
    HYPRE_BoomerAMGSetSCommPkgSwitch(amg_precond, S_commpkg_switch);
    HYPRE_BoomerAMGSetPrintLevel(amg_precond, poutdat);
    HYPRE_BoomerAMGSetPrintFileName(amg_precond, "driver.out.log");
    HYPRE_BoomerAMGSetMaxIter(amg_precond, 1);
    HYPRE_BoomerAMGSetCycleType(amg_precond, cycle_type);
    HYPRE_BoomerAMGSetFCycle(amg_precond, fcycle);
    HYPRE_BoomerAMGSetNumSweeps(amg_precond, num_sweeps);
    HYPRE_BoomerAMGSetISType(amg_precond, IS_type);
    HYPRE_BoomerAMGSetNumCRRelaxSteps(amg_precond, num_CR_relax_steps);
    HYPRE_BoomerAMGSetCRRate(amg_precond, CR_rate);
    HYPRE_BoomerAMGSetCRStrongTh(amg_precond, CR_strong_th);
    HYPRE_BoomerAMGSetCRUseCG(amg_precond, CR_use_CG);


    HYPRE_BoomerAMGSetAddRelaxType(amg_precond, add_relax_type);
    HYPRE_BoomerAMGSetAddRelaxWt(amg_precond, add_relax_wt);
    HYPRE_BoomerAMGSetChebyOrder(amg_precond, cheby_order);
    HYPRE_BoomerAMGSetChebyFraction(amg_precond, cheby_fraction);
    HYPRE_BoomerAMGSetChebyEigEst(amg_precond, cheby_eig_est);
    HYPRE_BoomerAMGSetChebyVariant(amg_precond, cheby_variant);
    HYPRE_BoomerAMGSetChebyScale(amg_precond, cheby_scale);
    HYPRE_BoomerAMGSetRelaxOrder(amg_precond, relax_order);
    HYPRE_BoomerAMGSetRelaxWt(amg_precond, relax_wt);
    HYPRE_BoomerAMGSetOuterWt(amg_precond, outer_wt);

    HYPRE_BoomerAMGSetSmoothType(amg_precond, smooth_type);
    HYPRE_BoomerAMGSetSmoothNumLevels(amg_precond, smooth_num_levels);
    HYPRE_BoomerAMGSetSmoothNumSweeps(amg_precond, smooth_num_sweeps);
    HYPRE_BoomerAMGSetMaxLevels(amg_precond, max_levels);
    HYPRE_BoomerAMGSetMaxRowSum(amg_precond, max_row_sum);
    HYPRE_BoomerAMGSetDebugFlag(amg_precond, debug_flag);
    HYPRE_BoomerAMGSetNumFunctions(amg_precond, num_functions);
    HYPRE_BoomerAMGSetAggNumLevels(amg_precond, agg_num_levels);
    HYPRE_BoomerAMGSetAggInterpType(amg_precond, agg_interp_type);
    HYPRE_BoomerAMGSetAggTruncFactor(amg_precond, agg_trunc_factor);
    HYPRE_BoomerAMGSetAggP12TruncFactor(amg_precond, agg_P12_trunc_factor);
    HYPRE_BoomerAMGSetAggPMaxElmts(amg_precond, agg_P_max_elmts);
    HYPRE_BoomerAMGSetAggP12MaxElmts(amg_precond, agg_P12_max_elmts);
    HYPRE_BoomerAMGSetNumPaths(amg_precond, num_paths);
    HYPRE_BoomerAMGSetNodal(amg_precond, nodal);
    HYPRE_BoomerAMGSetNodalDiag(amg_precond, nodal_diag);
    HYPRE_BoomerAMGSetVariant(amg_precond, variant);
    HYPRE_BoomerAMGSetOverlap(amg_precond, overlap);
    HYPRE_BoomerAMGSetDomainType(amg_precond, domain_type);
    HYPRE_BoomerAMGSetSchwarzUseNonSymm(amg_precond, use_nonsymm_schwarz);
    HYPRE_BoomerAMGSetSchwarzRlxWeight(amg_precond, schwarz_rlx_weight);

    HYPRE_BoomerAMGSetEuLevel(amg_precond, eu_level);
    HYPRE_BoomerAMGSetEuBJ(amg_precond, eu_bj);
    HYPRE_BoomerAMGSetEuSparseA(amg_precond, eu_sparse_A);
    HYPRE_BoomerAMGSetCycleNumSweeps(amg_precond, ns_coarse, 3);


    HYPRE_BoomerAMGSetAdditive(amg_precond, additive);
    HYPRE_BoomerAMGSetMultAdditive(amg_precond, mult_add);
    HYPRE_BoomerAMGSetSimple(amg_precond, simple);
    HYPRE_BoomerAMGSetAddLastLvl(amg_precond, add_last_lvl);
    HYPRE_BoomerAMGSetMultAddPMaxElmts(amg_precond, add_P_max_elmts);
    HYPRE_BoomerAMGSetMultAddTruncFactor(amg_precond, add_trunc_factor);
    HYPRE_BoomerAMGSetRAP2(amg_precond, rap2);
    HYPRE_BoomerAMGSetModuleRAP2(amg_precond, mod_rap2);
    HYPRE_BoomerAMGSetKeepTranspose(amg_precond, keepTranspose);
}

void hypreCalculator::destoryTempDataInAMGDD(HYPRE_Solver& amgddsolver)
{
    HYPRE_Solver amg_solver;

    HYPRE_BoomerAMGDDGetAMG(amgddsolver, &amg_solver);

    hypre_ParAMGData* amg_data = (hypre_ParAMGData*)(amg_solver);

    HYPRE_Int     num_levels = hypre_ParAMGDataNumLevels(amg_data);

    if (hypre_ParAMGDataDofFuncArray(amg_data))
    {
        for (int i = 1; i < num_levels; i++)
        {
            hypre_IntArrayDestroy(hypre_ParAMGDataDofFuncArray(amg_data)[i]);
        }
        hypre_TFree(hypre_ParAMGDataDofFuncArray(amg_data), HYPRE_MEMORY_HOST);
        hypre_ParAMGDataDofFuncArray(amg_data) = NULL;
    }

    for (int i = 1; i < num_levels; i++)
    {
        hypre_ParVectorDestroy(hypre_ParAMGDataFArray(amg_data)[i]);
        hypre_ParVectorDestroy(hypre_ParAMGDataUArray(amg_data)[i]);

        if (hypre_ParAMGDataAArray(amg_data)[i])
        {
            hypre_ParCSRMatrixDestroy(hypre_ParAMGDataAArray(amg_data)[i]);
        }

        if (hypre_ParAMGDataPArray(amg_data)[i - 1])
        {
            hypre_ParCSRMatrixDestroy(hypre_ParAMGDataPArray(amg_data)[i - 1]);
        }

        if (hypre_ParAMGDataRestriction(amg_data))
        {
            if (hypre_ParAMGDataRArray(amg_data)[i - 1])
            {
                hypre_ParCSRMatrixDestroy(hypre_ParAMGDataRArray(amg_data)[i - 1]);
            }
        }

        hypre_IntArrayDestroy(hypre_ParAMGDataCFMarkerArray(amg_data)[i - 1]);

        /* get rid of any block structures */
        if (hypre_ParAMGDataABlockArray(amg_data)[i])
        {
            hypre_ParCSRBlockMatrixDestroy(hypre_ParAMGDataABlockArray(amg_data)[i]);
        }

        if (hypre_ParAMGDataPBlockArray(amg_data)[i - 1])
        {
            hypre_ParCSRBlockMatrixDestroy(hypre_ParAMGDataPBlockArray(amg_data)[i - 1]);
        }

        /* RL */
        if (hypre_ParAMGDataRestriction(amg_data))
        {
            if (hypre_ParAMGDataRBlockArray(amg_data)[i - 1])
            {
                hypre_ParCSRBlockMatrixDestroy(hypre_ParAMGDataRBlockArray(amg_data)[i - 1]);
            }
        }
    }

    hypre_ParAMGDataAArray(amg_data) = NULL;
    hypre_ParAMGDataPArray(amg_data) = NULL;
    hypre_ParAMGDataRArray(amg_data) = NULL;
    hypre_ParAMGDataCFMarkerArray(amg_data) = NULL;
    hypre_ParAMGDataFArray(amg_data) = NULL;
    hypre_ParAMGDataUArray(amg_data) = NULL;
    hypre_ParAMGDataDofFuncArray(amg_data) = NULL;

    hypre_ParAMGDataABlockArray(amg_data) = NULL;
    hypre_ParAMGDataPBlockArray(amg_data) = NULL;
    hypre_ParAMGDataRBlockArray(amg_data) = NULL;

}
#endif

}
