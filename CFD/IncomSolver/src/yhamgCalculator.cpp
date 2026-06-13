#include <iostream>
#include <iomanip>
#include "yhamgCalculator.h"
#include "TK_Time.h"
#include "TK_Exit.h"

namespace PHSPACE
{
    yhamgCalculator::yhamgCalculator() : LinearEquationSystemCalculator()
    {
#ifndef USE_YhamgLib
        TK_Exit::ExceptionExit("yhamg lib is not opened in CMAKEFile!");
#endif
    }

    yhamgCalculator::~yhamgCalculator()
    {
        mathLibSolveFinalize();
    }

    void yhamgCalculator::mathLibSolveInit(Grid *gridIn)
    {
#ifdef USE_YhamgLib
        UnstructGrid *grid = UnstructGridCast(gridIn);
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotalFace = grid->GetNTotalFace();
        int *leftCellOfFace = grid->GetLeftCellOfFace();
        int *rightCellOfFace = grid->GetRightCellOfFace();
        int zoneID = grid->GetZoneID();
        int *localCell2LocalMap = grid->GetLocalCell2LocalMap();
        int *globalCellNoOfGlobalZones = reinterpret_cast<int *>(GlobalDataBase::GetDataPtr("globalCellNoOfGlobalZones"));
        int *precondID = reinterpret_cast<int *>(GlobalDataBase::GetDataPtr("precondMethodOfLinEqSystem"));

        int linearSolverSize = 50;
        A = new ParCSRMatrix;
        PrecondList = new ParOperator* [linearSolverSize];
        for (int iSolver = 0 ; iSolver < linearSolverSize; iSolver++)
        {
            int precondType = precondID[iSolver];
            ParOperator *precond = PrecondList[iSolver];
            if (precondType == iNone)
            {
                precond = nullptr;
            }
            else if (precondType == iILU)
            {
                precond = new ParCSRPrecondILU0();
            }
            else if (precondType == iAMG)
            {
                precond = new ParCSRPrecondAMG();
                int    max_levels = 20;
                int    coarse_size = 8;
                RDouble strength_threshold = 0.25;
                int    aggressive_levels = 1;
                int    coarsen_type = 0; // 0: HMIS, 1: PMIS
                int    interp_type = 0;  // 0: Long-Range Interpolation, 1: Smoothed Aggregation
                int    interp_min_elements = 4;
                int    interp_max_elements = 6;
                RDouble truncation_factor = 0.1;
                RDouble sparsification_threshold = 0.01;
                int    cycle_type = 0; // 0: V-Cycle, 1: W-Cycle, 2: F-Cycle
                int    pre_sweeps = 2;
                int    post_sweeps = 2;
                int    coarse_sweeps = 2;
                int    smooth_type = 1; // 0: Jacobi, 1: SOR, 2: ILU, 3: Chebyshev

                RDouble jacobi_factor = 0.75;

                int    sor_type = 2; // 0: Forward, 1: Backword, 2: Symmetric
                RDouble sor_factor = 1.0;

                int    ilu_maxfil = 0;
                RDouble ilu_droptol = 0.01;

                int    chebyshev_order = 4;
                RDouble chebyshev_eigen_ratio = 0.3;
                bool   print_stats = 0;
                ((ParCSRPrecondAMG*)precond)->MaxLevels = max_levels;
                ((ParCSRPrecondAMG*)precond)->CoarseSize = coarse_size;
                ((ParCSRPrecondAMG*)precond)->StrengthThreshold = strength_threshold;
                ((ParCSRPrecondAMG*)precond)->AggressiveLevels = aggressive_levels;
                ((ParCSRPrecondAMG*)precond)->CoarsenType = coarsen_type;
                ((ParCSRPrecondAMG*)precond)->InterpType = interp_type;
                ((ParCSRPrecondAMG*)precond)->InterpMinElements = interp_min_elements;
                ((ParCSRPrecondAMG*)precond)->InterpMaxElements = interp_max_elements;
                ((ParCSRPrecondAMG*)precond)->TruncationFactor = truncation_factor;
                ((ParCSRPrecondAMG*)precond)->SparsificationThreshold = sparsification_threshold;
                ((ParCSRPrecondAMG*)precond)->CycleType = cycle_type;
                ((ParCSRPrecondAMG*)precond)->PreSweeps = pre_sweeps;
                ((ParCSRPrecondAMG*)precond)->PostSweeps = post_sweeps;
                ((ParCSRPrecondAMG*)precond)->CoarseSweeps = coarse_sweeps;
                ((ParCSRPrecondAMG*)precond)->SmoothType = smooth_type;
                ((ParCSRPrecondAMG*)precond)->JacobiFactor = jacobi_factor;
                ((ParCSRPrecondAMG*)precond)->SORType = sor_type;
                ((ParCSRPrecondAMG*)precond)->SORFactor = sor_factor;
                ((ParCSRPrecondAMG*)precond)->ILUMaxFillins = ilu_maxfil;
                ((ParCSRPrecondAMG*)precond)->ILUDropTolerance = ilu_droptol;
                ((ParCSRPrecondAMG*)precond)->ChebyshevEigenRatio = chebyshev_eigen_ratio;
                ((ParCSRPrecondAMG*)precond)->PrintStats = print_stats;


            }
            else if (precondType == iBJILU)
            {
                precond = new ParCSRPrecondBJILU();
            }
            else if (precondType == iBJSOR)
            {
                precond = new ParCSRPrecondBJSOR();
            }
            else if (precondType == iJacobi)
            {
                precond = new ParCSRPrecondJacobi();
            }
            else if (precondType == iSOR)
            {
                precond = new ParCSRPrecondSOR();
            }
            PrecondList[iSolver] = precond;
        }

        int * local_rowptr = new int[nTotalCell + 1];

        ColVectors = new vector<int>[nTotalCell];
        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            ColVectors[iCell].push_back(iCell);
        }

        for (int iFace = 0; iFace < nTotalFace; iFace++)
        {
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];
            if (re < nTotalCell)
            {
                ColVectors[le].push_back(re);
                ColVectors[re].push_back(le);
            }
        }

        int nonZeroNum = 0;
        local_rowptr[0] = nonZeroNum;
        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            nonZeroNum += ColVectors[iCell].size();
            local_rowptr[iCell + 1] = nonZeroNum;
        }

        int * local_colind = new int[nonZeroNum];
        int count = 0;
        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            for (int jCell = 0; jCell < ColVectors[iCell].size(); jCell++)
            {
                local_colind[count] = ColVectors[iCell][jCell];
                count += 1;
            }
        }

        RDouble *local_values = new RDouble[nonZeroNum];

        b = new ParVector(MPI_COMM_WORLD);
        b->Resize(nTotalCell);

        x = new ParVector(MPI_COMM_WORLD);
        x->Resize(nTotalCell);

        ParCSRMatrix* CSRMatrixA = reinterpret_cast <ParCSRMatrix*>(A);
        CSRMatrixA->Free();
        CSRMatrixA->comm = MPI_COMM_WORLD;
        CSRMatrixA->local.size[0] = nTotalCell;
        CSRMatrixA->local.size[1] = nTotalCell;
        CSRMatrixA->local.rowptr = local_rowptr;
        CSRMatrixA->local.colind = local_colind;
        CSRMatrixA->local.values = local_values;

        int nnb = 0;
        int *nbrank;
        int *recvptr;

        InterfaceInfo * interfaceInfo= grid->GetInterfaceInfo();
        if (interfaceInfo != NULL) 
        {
            int numberOfNeighbor = interfaceInfo->GetNumberOfNeighbor();
            nnb = numberOfNeighbor;
        }

        nbrank = new int[nnb];
        if (interfaceInfo != NULL) 
        {
            int numberOfNeighbor = interfaceInfo->GetNumberOfNeighbor();
            vector <InterfacePatch *> interfacePatch = interfaceInfo->GetInterfacePatch();
            for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; iNeighbor++)
            {
                int neighborID = interfacePatch[iNeighbor]->GetZoneIndexOfNeighbor();
                nbrank[iNeighbor] = neighborID;
            }
        }

        recvptr = new int[nnb + 1];
        int faceCount = 0;
        recvptr[0] = faceCount;
        if (interfaceInfo != NULL) 
        {
            int numberOfNeighbor = interfaceInfo->GetNumberOfNeighbor();
            vector <InterfacePatch *> interfacePatch = interfaceInfo->GetInterfacePatch();
            for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; iNeighbor++)
            {
                int nPatchFace = interfacePatch[iNeighbor]->GetNumberOfFace();
                faceCount += nPatchFace;
                recvptr[iNeighbor + 1] = faceCount;
            }
        }

        int * exter_colind = new int[faceCount];
        exter_colind2Face = new int[faceCount];
        RDouble *exter_values = new RDouble[faceCount];
        int * recvind = new int[faceCount];
        vector<int>* exterMatrixInfo = new vector<int>[nTotalCell];
        vector<int>* exterMatrixInfo2Face = new vector<int>[nTotalCell];

        faceCount = 0;
        if (interfaceInfo != NULL) 
        {
            int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();
            int numberOfNeighbor = interfaceInfo->GetNumberOfNeighbor();
            vector <InterfacePatch *> interfacePatch = interfaceInfo->GetInterfacePatch();
            for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; iNeighbor++)
            {
                int nPatchFace = interfacePatch[iNeighbor]->GetNumberOfFace();
                int neighborID = interfacePatch[iNeighbor]->GetZoneIndexOfNeighbor();

                for (int iPatchFace = 0; iPatchFace < nPatchFace; iPatchFace++)
                {
                    int *faceIndexForSend = interfacePatch[iNeighbor]->GetFaceIndexForSend();
                    int sendFaceIndex = faceIndexForSend[iPatchFace];
                    int faceIndex = interFace2BoundaryFace[sendFaceIndex];
                    int le = leftCellOfFace[faceIndex];
                    int re = rightCellOfFace[faceIndex];
                    recvind[faceCount] = localCell2LocalMap[re];
                    if (localCell2LocalMap[re] > -1) 
                    {
                        exterMatrixInfo[le].push_back(faceCount);
                        exterMatrixInfo2Face[le].push_back(faceIndex);
                    }
                    faceCount += 1;
                }
            }
        }

        int *exter_rowptr = new int[nTotalCell + 1];
        exter_rowptr[0] = 0;
        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            exter_rowptr[iCell + 1] = exter_rowptr[iCell] + exterMatrixInfo[iCell].size();
        }

        int cellCount = 0;
        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            for (int jCell = exter_rowptr[iCell]; jCell < exter_rowptr[iCell + 1]; jCell++)
            {
                exter_colind[cellCount] = exterMatrixInfo[iCell][jCell - exter_rowptr[iCell]];
                exter_colind2Face[cellCount] = exterMatrixInfo2Face[iCell][jCell - exter_rowptr[iCell]];
                cellCount += 1;
            }
        }

        delete [] exterMatrixInfo;
        delete [] exterMatrixInfo2Face;

        CSRMatrixA->exter.size[0] = nTotalCell;
        CSRMatrixA->exter.size[1] = recvptr[nnb];
        CSRMatrixA->exter.rowptr = exter_rowptr;
        CSRMatrixA->exter.colind = exter_colind;
        CSRMatrixA->exter.values = exter_values;
        CSRMatrixA->nnb = nnb;
        CSRMatrixA->nbrank = nbrank;
        CSRMatrixA->recvptr = recvptr;
        CSRMatrixA->recvind = recvind;
#endif
    }


    void yhamgCalculator::mathLibSolve(Grid *gridIn, int iEquation, RDouble *fi, RDouble *diagCoef, RDouble *bold, RDouble *upNeighCoef, RDouble *dnNeighCoef)
    {
#ifdef USE_YhamgLib
        UnstructGrid *grid = UnstructGridCast(gridIn);
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotalFace = grid->GetNTotalFace();
        int *leftCellOfFace = grid->GetLeftCellOfFace();
        int *rightCellOfFace = grid->GetRightCellOfFace();

        int *localCell2GlobalMap = grid->GetLocalCell2GlobalMap();
        RDouble GlobalTotalCells = GlobalDataBase::GetDoubleParaFromDB("GlobalTotalCells");
        int *linearSolverID = reinterpret_cast<int *>(GlobalDataBase::GetDataPtr("solMethodOfLinEqSystem"));
        int *precondID = reinterpret_cast<int *>(GlobalDataBase::GetDataPtr("precondMethodOfLinEqSystem"));
        int *maxIterOfEachEq = reinterpret_cast<int *>(GlobalDataBase::GetDataPtr("maxIterOfEachEq"));
        RDouble *iterTolOfEachEq = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("iterTolOfEachEq"));

        int maxIterStep = maxIterOfEachEq[iEquation];
        RDouble iterTol = iterTolOfEachEq[iEquation];
        int *ActualIterStepNum = reinterpret_cast <int * > (GlobalDataBase::GetDataPtr("ActualIterStepNum"));
        RDouble *ActualIterTol = reinterpret_cast <RDouble * > (GlobalDataBase::GetDataPtr("ActualIterTol"));
        RDouble *ActualSolveTime = reinterpret_cast <RDouble * > (GlobalDataBase::GetDataPtr("ActualSolveTime"));
        int linearSolver = linearSolverID[iEquation];
        RDouble timeStartSolve = PHSPACE::TIME_SPACE::GetWallTime();

        int    gmres_restart = 10;
        bool   print_stats = 0;
        int iter;
        RDouble relres;
        RDouble setup_time;
        RDouble solve_time;
        ParCSRMatrix* CSRMatrixA = reinterpret_cast <ParCSRMatrix*>(A);
        UnstructBCSet **bcr = grid->GetBCRecord();
        vector<RDouble>* ValuesVector = new vector<RDouble>[nTotalCell];
        RDouble *local_values = CSRMatrixA->local.values;

        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            ValuesVector[iCell].push_back(diagCoef[iCell]);
        }
        for (int iFace = 0; iFace < nTotalFace; iFace++)
        {
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];
            if (re < nTotalCell)
            {
                ValuesVector[le].push_back(upNeighCoef[iFace]);
                ValuesVector[re].push_back(dnNeighCoef[iFace]);
            }
        }

        int cellCount = 0;
        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            for (int jCell = 0; jCell < ColVectors[iCell].size(); jCell++)
            {
                local_values[cellCount] = ValuesVector[iCell][jCell];
                cellCount += 1;
            }
        }

        cellCount = 0;
        int *exter_rowptr = CSRMatrixA->exter.rowptr;
        RDouble *exter_values = CSRMatrixA->exter.values;
        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            for (int jCell = exter_rowptr[iCell]; jCell < exter_rowptr[iCell + 1]; jCell++)
            {
                int faceIndex = exter_colind2Face[cellCount];
                exter_values[cellCount] = upNeighCoef[faceIndex];
                cellCount += 1;
            }
        }

        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            (*x)[iCell] =  fi[iCell];
            (*b)[iCell] =  bold[iCell] + SMALL;
        }

        reinterpret_cast <ParCSRMatrix*>(A)->SetupHalo();

        ParOperator* precond = PrecondList[iEquation];
        int precondType = precondID[iEquation];

        if (linearSolver == GaussSeidelMethod) 
        {
            explicitSolve(gridIn, fi, diagCoef, bold, upNeighCoef, dnNeighCoef);
        }
        else if (linearSolver == GMRESMethod) 
        {
            if (precondType != 0)
            {
                ((ParCSRPrecond*)precond)->Setup(*(ParCSRMatrix*)A);
                ParSolverGMRES(gmres_restart, maxIterStep, iterTol, print_stats)(*A, *precond, *b, *x, iter, relres);
            }
            else if (precondType == 0)
            {
                ParSolverGMRES(gmres_restart, maxIterStep, iterTol, print_stats)(*A, *b, *x, iter, relres);
            }
        }
        else if (linearSolver == BicgstabMethod)
        {
            if (precondType != 0)
            {
                ((ParCSRPrecond*)precond)->Setup(*(ParCSRMatrix*)A);
                ParSolverBCGS(maxIterStep, iterTol, print_stats)(*A, *precond, *b, *x, iter, relres);
            }
            else
            {
                ParSolverBCGS(maxIterStep, iterTol, print_stats)(*A, *b, *x, iter, relres);
            }
        }
        else if (linearSolver == PCGMethod)
        {
            if (precondType != 0)
            {
                ((ParCSRPrecond*)precond)->Setup(*(ParCSRMatrix*)A);
                ParSolverCG(maxIterStep, iterTol, print_stats)(*A, *precond, *b, *x, iter, relres);
            }
            else
            {
                ParSolverCG(maxIterStep, iterTol, print_stats)(*A, *b, *x, iter, relres);
            }
        }
        else if (linearSolver == PipeBicgstabMethod)
        {
            if (precondType != 0)
            {
                ((ParCSRPrecond*)precond)->Setup(*(ParCSRMatrix*)A);
                ParSolverPipeBCGS(maxIterStep, iterTol, print_stats)(*A, *precond, *b, *x, iter, relres);
            }
            else
            {
                ParSolverPipeBCGS(maxIterStep, iterTol, print_stats)(*A, *b, *x, iter, relres);
            }
        }
        else if (linearSolver == PipeCGMethod)
        {
            if (precondType != 0)
            {
                ((ParCSRPrecond*)precond)->Setup(*(ParCSRMatrix*)A);
                ParSolverPipeCG(maxIterStep, iterTol, print_stats)(*A, *precond, *b, *x, iter, relres);
            }
            else
            {
                ParSolverPipeCG(maxIterStep, iterTol, print_stats)(*A, *b, *x, iter, relres);
            }
        }
        else if (linearSolver == PipeGMRESMethod) 
        {
            if (precondType != 0)
            {
                ((ParCSRPrecond*)precond)->Setup(*(ParCSRMatrix*)A);
                ParSolverPipeGMRES(gmres_restart, maxIterStep, iterTol, print_stats)(*A, *precond, *b, *x, iter, relres);
            }
            else
            {
                ParSolverPipeGMRES(gmres_restart, maxIterStep, iterTol, print_stats)(*A, *b, *x, iter, relres);
            }
        }

        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            fi[iCell] = (*x)[iCell];
        }

        //RDouble timeEndSolve = PHSPACE::TIME_SPACE::GetWallTime();
        //ActualIterStepNum[iEquation]] = iter;
        //ActualIterTol[constVarNumToSolveVarNum[iEquation]] = relres;
        //ActualSolveTime[constVarNumToSolveVarNum[iEquation]] = timeEndSolve - timeStartSolve;

        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            ValuesVector[iCell].clear();
        }
        delete [] ValuesVector; 
#endif
    }

    void yhamgCalculator::mathLibSolveFinalize()
    {
#ifdef USE_YhamgLib
        ParCSRMatrix* CSRMatrixA = reinterpret_cast <ParCSRMatrix*>(A);
        if (A)
        {
            CSRMatrixA->Free();
            delete [] CSRMatrixA->local.rowptr;
            delete [] CSRMatrixA->local.colind;
            delete [] CSRMatrixA->local.values;
            delete [] CSRMatrixA->exter.rowptr;
            delete [] CSRMatrixA->exter.colind;
            delete [] CSRMatrixA->exter.values;
            delete A;
        }
        int linearSolverSize = GlobalDataBase::GetIntParaFromDB("linearSolverSize");
        int *precondID = reinterpret_cast <int *>(GlobalDataBase::GetDataPtr("precondID"));
        for (int iSolver = 0 ; iSolver < linearSolverSize; iSolver++)
        {
            int precondType = precondID[iSolver];
            if (precondType > 0)
            {
                delete PrecondList[iSolver];
            }
        }

        delete [] PrecondList;
#endif
    }

}
