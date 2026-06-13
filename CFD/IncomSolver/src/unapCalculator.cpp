#include <iostream>
#include <iomanip>
#include "unapCalculator.h"
#include "TK_Time.h"
#include "TK_Exit.h"

namespace PHSPACE
{
unapCalculator::unapCalculator() : LinearEquationSystemCalculator()
{
#ifndef USE_UnapLib
    TK_Exit::ExceptionExit("unap lib is not opened in CMAKEFile!");
#endif
}

unapCalculator::~unapCalculator()
{
    mathLibSolveFinalize();
}

void unapCalculator::mathLibSolveInit(Grid *gridIn)
{
#ifdef USE_UnapLib
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    
    int nInternalFaces = nTotalFace - nBoundFace;
    
#if defined(REORDER)
#else
    initUtility();
#endif
    
    Communicator* comm = &COMM::getGlobalComm();
    
    using namespace UNAP;
    
    label nProcs = comm->getMySize();
    label myId = comm->getMyId();
    
    lduA = new LduMatrix(comm);
    upperAddr = new labelVector(nInternalFaces, comm);
    lowerAddr = new labelVector(nInternalFaces, comm);
    upper = new scalarVector(nInternalFaces, comm);
    lower = new scalarVector(nInternalFaces, comm);
    diag  = new scalarVector(nTotalCell, comm);
    b = new scalarVector(nTotalCell, comm);
    x = new scalarVector(nTotalCell, 1.0, comm);
    
    for (int iFace = nBoundFace; iFace < nTotalFace; iFace++)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        int jFace = iFace - nBoundFace;
    
        (*lowerAddr)[jFace] = le;
        (*upperAddr)[jFace] = re;
    }
    
    lduA->setMatrixTopology(nTotalCell, *upperAddr, *lowerAddr, false);

    int *localCell2GlobalMap = grid->GetLocalCell2GlobalMap();
    InterfaceInfo * interfaceInfo= grid->GetInterfaceInfo();
    if (interfaceInfo != NULL) 
    {
        nNeiProcs = interfaceInfo->GetNumberOfNeighbor();
        destRank = new int[nNeiProcs];
        int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();
        vector <InterfacePatch *> interfacePatch = interfaceInfo->GetInterfacePatch();
        for (int iNeighbor = 0; iNeighbor < nNeiProcs; iNeighbor++)
        {
            label neighborID = interfacePatch[iNeighbor]->GetZoneIndexOfNeighbor();
            destRank[iNeighbor] = neighborID;
        }

        int nIterFace = interfaceInfo->GetNIFace();
        offDiagRows = new int[nIterFace];
        offDiagCols = new int[nIterFace];
        offDiagCoeffs = new double[nIterFace];
        offDiagStarts = new int[nNeiProcs + 1];
        offDiagStarts[0] = 0;
        int InterFaceCount = 0;
        for (int iNeighbor = 0; iNeighbor < nNeiProcs; iNeighbor++)
        {
            label nPatchFace = interfacePatch[iNeighbor]->GetNumberOfFace();

            label neighborID = interfacePatch[iNeighbor]->GetZoneIndexOfNeighbor();

            for (int iPatchFace = 0; iPatchFace < nPatchFace; iPatchFace++)
            {
                int *faceIndexForSend = interfacePatch[iNeighbor]->GetFaceIndexForSend();
                int sendFaceIndex = faceIndexForSend[iPatchFace];
                int faceIndex = interFace2BoundaryFace[sendFaceIndex];
                int le = leftCellOfFace[faceIndex];
                int re = rightCellOfFace[faceIndex];

                offDiagRows[InterFaceCount] = le;
                offDiagCols[InterFaceCount] = localCell2GlobalMap[re];

                InterFaceCount += 1;
            }
            offDiagStarts[iNeighbor + 1] = InterFaceCount;
        }
        lduA->createInterfacesTopology(nNeiProcs, destRank, offDiagRows, offDiagCols, offDiagStarts);
    }
#endif
}

void unapCalculator::mathLibSolve(Grid *gridIn, int iEquation, RDouble *fi, RDouble *diagCoef, RDouble *binput, RDouble *upNeighCoef, RDouble *dnNeighCoef)
{
#ifdef USE_UnapLib
    UnstructGrid * grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int nInternalFaces = nTotalFace - nBoundFace;
    RDouble *area = grid->GetFaceArea();
    int *localCell2GlobalMap = grid->GetLocalCell2GlobalMap();

    Communicator *comm = &COMM::getGlobalComm();
    label nProcs = comm->getMySize();
    label myId = comm->getMyId();
    RDouble tol = 0.0;
    int minIter = 1;
    int *linearSolverID = reinterpret_cast<int *>(GlobalDataBase::GetDataPtr("solMethodOfLinEqSystem"));
    int *precondID = reinterpret_cast<int *>(GlobalDataBase::GetDataPtr("precondMethodOfLinEqSystem"));
    int *maxIterOfEachEq = reinterpret_cast<int *>(GlobalDataBase::GetDataPtr("maxIterOfEachEq"));
    RDouble *iterTolOfEachEq = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("iterTolOfEachEq"));

    int maxIterStep = maxIterOfEachEq[iEquation];
    RDouble iterTol = iterTolOfEachEq[iEquation];
    bool isPrintLinearSolver = false;
    //int *ActualIterStepNum = reinterpret_cast <int * > (GlobalDataBase::GetDataPtr("ActualIterStepNum"));
    //RDouble *ActualIterTol = reinterpret_cast <RDouble * > (GlobalDataBase::GetDataPtr("ActualIterTol"));
    //RDouble *ActualSolveTime = reinterpret_cast <RDouble * > (GlobalDataBase::GetDataPtr("ActualSolveTime"));
    RDouble timeStartSolve = PHSPACE::TIME_SPACE::GetWallTime();

    for (int iFace = nBoundFace; iFace < nTotalFace; iFace++)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        int jFace = iFace - nBoundFace;
        (*lower)[jFace] = dnNeighCoef[iFace]; 
        (*upper)[jFace] = upNeighCoef[iFace];
    }
    
    for (int iCell = 0; iCell < nTotalCell; iCell++)
    {
        (*diag)[iCell] = diagCoef[iCell];
        (*b)[iCell] = binput[iCell];
        (*x)[iCell] = 1.0;
    }
    
    lduA->setMatrixCoeffients(*diag, *upper, *lower, false);
    
    InterfaceInfo * interfaceInfo= grid->GetInterfaceInfo();
    if (interfaceInfo != NULL) 
    {
        nNeiProcs = interfaceInfo->GetNumberOfNeighbor();
        int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();
        vector <InterfacePatch *> interfacePatch = interfaceInfo->GetInterfacePatch();
        int InterFaceCount = 0;
        for (int iNeighbor = 0; iNeighbor < nNeiProcs; iNeighbor++)
        {
            label nPatchFace = interfacePatch[iNeighbor]->GetNumberOfFace();
            for (int iPatchFace = 0; iPatchFace < nPatchFace; iPatchFace++)
            {
                int *faceIndexForSend = interfacePatch[iNeighbor]->GetFaceIndexForSend();
                int sendFaceIndex = faceIndexForSend[iPatchFace];
                int faceIndex = interFace2BoundaryFace[sendFaceIndex];

                offDiagCoeffs[InterFaceCount] = upNeighCoef[faceIndex];
                InterFaceCount += 1;
            }
        }
        lduA->fillInterfacesCofficients(offDiagStarts, offDiagCoeffs);
    }

    if (precondID[iEquation] == iAMG) 
    {
        unapPrecond = new LduMGPrecond(*lduA);
    }
    else if (precondID[iEquation] == LduDiag) 
    {
        unapPrecond = new LduDiagPrecond(*lduA);
    }
    else if (precondID[iEquation] == LduDIC)
    {
        unapPrecond = new LduDICPrecond(*lduA);
    }
    else if (precondID[iEquation] == LduDILU)
    {
        unapPrecond = new LduDILUPrecond(*lduA);
    }
    else if (precondID[iEquation] == LduBJacobi)
    {
        //precond = new LduBJacobiPrecond(*lduA);
    }
    
    #if defined(DSWACC)
    lduA->constructRSSIterator();
    #endif

    if (linearSolverID[iEquation] == GaussSeidelMethod) 
    {
        explicitSolve(gridIn, fi, diagCoef, binput, upNeighCoef, dnNeighCoef);
    }
    else if (linearSolverID[iEquation] == GMRESMethod) 
    {
        unapSolver = new GMRES(*unapPrecond, 10);
        unapSolver->SET_minIter(minIter);
        unapSolver->SET_maxIter(maxIterStep);
        unapSolver->SET_ifPrint(isPrintLinearSolver);
        unapSolver->SET_relTol(iterTol);
        unapSolver->SET_tolerance(tol);
    
        solverPerf = static_cast<GMRES *>(unapSolver)->solve(*x, *lduA, *b);
        static_cast<GMRES*>(unapSolver)->~GMRES();
    }
    else if (linearSolverID[iEquation] == BicgstabMethod)
    {
        unapSolver = new PBiCGStab(*unapPrecond);
        unapSolver->SET_minIter(minIter);
        unapSolver->SET_maxIter(maxIterStep);
        unapSolver->SET_ifPrint(isPrintLinearSolver);
        unapSolver->SET_relTol(iterTol);
        unapSolver->SET_tolerance(tol);
    
        solverPerf = unapSolver->solve(*x, *lduA, *b);
        static_cast<PBiCGStab*>(unapSolver)->~PBiCGStab();
    }
    else if (linearSolverID[iEquation] == AMGMethod) 
    {
        DirectionalAgglomeration aggl(*lduA);
        aggl.SET_maxLevels(50);
        aggl.SET_rowSizeInCoarsestLevel(200);
        aggl.agglomerate();
        
        //unapSolver = new MGSolver(*lduA, aggl);

        unapSolver->SET_tolerance(tol);
        unapSolver->SET_relTol(iterTol);
        static_cast<MGSolver*>(unapSolver)->SET_nPreSweeps(2);
        unapSolver->SET_minIter(minIter);
        unapSolver->SET_maxIter(maxIterStep);
        unapSolver->SET_ifPrint(isPrintLinearSolver);
        //static_cast<MGSolver*>(unapSolver)->SET_smootherType("Jacobi");
        static_cast<MGSolver*>(unapSolver)->initSmoothers();
        static_cast<MGSolver*>(unapSolver)->setUp();
        solverPerf = static_cast<MGSolver*>(unapSolver)->solve(*x, *lduA, *b);
        
        static_cast<MGSolver*>(unapSolver)->~MGSolver();
    }
    else if (linearSolverID[iEquation] == PCGMethod) 
    {
        unapSolver = new PCG(*unapPrecond);
        unapSolver->SET_minIter(minIter);
        unapSolver->SET_maxIter(maxIterStep);
        unapSolver->SET_ifPrint(isPrintLinearSolver);
        unapSolver->SET_tolerance(tol);
        unapSolver->SET_relTol(iterTol);
        
        solverPerf = static_cast<PCG *>(unapSolver)->solve(*x, *lduA, *b);
        static_cast<PCG*>(unapSolver)->~PCG();
    }
    

    RDouble timeEndSolve = PHSPACE::TIME_SPACE::GetWallTime();
    //ActualIterStepNum[constVarNumToSolveVarNum[iEquation]] = solverPerf.nIterations();
    //ActualIterTol[constVarNumToSolveVarNum[iEquation]] = solverPerf.finalResidual();
    //ActualSolveTime[constVarNumToSolveVarNum[iEquation]] = timeEndSolve - timeStartSolve;

    if (precondID[iEquation] == iAMG) 
    {
        static_cast<LduMGPrecond*>(unapPrecond)->~LduMGPrecond();
    }
    else if (precondID[iEquation] == LduDiag) 
    {
        static_cast<LduDiagPrecond*>(unapPrecond)->~LduDiagPrecond();
    }
    else if (precondID[iEquation] == LduDIC)
    {
        static_cast<LduDICPrecond*>(unapPrecond)->~LduDICPrecond();
    }
    else if (precondID[iEquation] == LduDILU)
    {
        static_cast<LduDILUPrecond*>(unapPrecond)->~LduDILUPrecond();
    }
    else if (precondID[iEquation] == LduBJacobi)
    {
        //static_cast<LduBJacobiPrecond*>(precond)->~LduBJacobiPrecond();
    }

    for (int iCell = 0; iCell < nTotalCell; iCell++)
    {
        fi[iCell] = (*x)[iCell];
    }
#endif
}

void unapCalculator::mathLibSolveFinalize()
{
#ifdef USE_UnapLib

    DELETE_POINTERS(upperAddr);
    DELETE_POINTERS(lowerAddr);
    DELETE_POINTERS(upper);
    DELETE_POINTERS(lower);
    DELETE_POINTERS(diag);
    DELETE_POINTERS(b);
    DELETE_POINTERS(x);

    DELETE_POINTERS(destRank);
    DELETE_POINTERS(offDiagRows);
    DELETE_POINTERS(offDiagCols);
    DELETE_POINTERS(offDiagStarts);
    DELETE_POINTERS(offDiagCoeffs);

    lduA->~LduMatrix();

    delete unapSolver ;
    delete unapPrecond;

#endif
}

}
