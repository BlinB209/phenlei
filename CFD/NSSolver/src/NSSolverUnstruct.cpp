#include "NSSolverUnstruct.h"
#include "MultiGridOperation.h"
#include "GradientOperation.h"
#include "Math_Limiter.h"
#include "AerodynamicsBasicRelations.h"
#include "Glb_Dimension.h"
#include "Constants.h"
#include <list>
#include "Param_NSSolverUnstruct.h"
#include "Post_ForceMoment.h"
#include "Residual.h"
#include "Jacobian.h"
#include "FieldProxy.h"
#include "TK_Exit.h"
#include "TK_Log.h"
#include "IO_FileName.h"
#include "PHIO.h"
#include "Force.h"
#include "AleForceManager.h"
#include <cmath>
#include "Flux_RoeEntropyCorrection.h"
#include "PHMpi.h"
#include "TK_Time.h"
#include "Geo_GridIndex.h"
#ifdef USE_CUDA
#include "GPUKernelTest.h"
#include "GPUInviscidFlux.h"
#include "TemporaryOperations.h"
#include "BasicDeviceVariables.h"
#include "GPUNSSolver.h"
#include "GPUInvVisSpectrum.h"
#include "GPUMPICommunication.h"
#include "GPUTurbAhead.h"
#include "GPUComputeNodeValue.h"
#include "GPUFaceWeight.h"
#include "GPUGetVisFaceValue.h"
#include "GPUCompVisfluxTEST.h"
#include "GPULUSGS.h"
#include "GPUSpectrumRadius.h"
#include "GPULoadFlux.h"
#include "GPUCompGradientGGNode.h"
#include "GPUCompGradientGGCell.h"
#include "GPUCompGradientGGNodeNew.h"
#include "GPUCompGradientGGCell.h"
#include "TemporaryOperations.h"
#include "TemporaryOperationsPart2.h"
using namespace TemporaryOperations;
using namespace GPUNSSolverUnstruct;
using namespace GPUKernels;
#endif

#ifdef CUDAUNITTEST
#include "GPUKernelTestPart2.h"
#include "GPUTestFunctions.h"
#endif

#include "MixedGas.h"
#include "Gas.h"
#ifdef USE_GMRESSOLVER
 #include "petscksp.h"  
#endif
using namespace std;
#pragma warning(disable:6386)
#pragma warning(disable:6385)
#pragma warning(disable:26451)
#pragma warning (disable:913)

namespace PHSPACE
{
using namespace GAS_SPACE;
using namespace IDX;

NSSolverUnstruct::NSSolverUnstruct()
{
    gradientNSField = 0;
    gradientTemperature = 0;
}

NSSolverUnstruct::~NSSolverUnstruct()
{
    DeAllocateGlobalVariables();
    FreeControlParameters();
}

void NSSolverUnstruct::ReadParameter()
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    string uns_scheme_name = parameters->GetSpaceSecheme();

    int uns_scheme = GetSchemeID(uns_scheme_name);
    GlobalDataBase::UpdateData("uns_scheme", &uns_scheme, PHINT, 1);

    //! GMRES
    if (uns_scheme == ISCHEME_ROE || uns_scheme == ISCHEME_GMRES_ROE)
    {
        parameters->SetEntropyFixCoefficients();

        Grid *grid = GetGrid();
        int zoneID = grid->GetZoneID();
        RDouble entropyAcoustic = parameters->GetRoeEntropyFixCoef1();

        //! The output of Entropy fix coefficient for ROE scheme.
        if (zoneID == 0)    //! When grid->GetZoneID() = 0 and GetCurrentProcessorID() = 0.
        {
            PrintToWindow("Entropy fix coefficient for ROE scheme: ", entropyAcoustic, "\n");
        }
    }

    string uns_limiter_name;
    GlobalDataBase::GetData("uns_limiter_name", &uns_limiter_name, PHSTRING, 1);

    int limiterType = GetLimiterID(uns_limiter_name);
    parameters->SetLimiterType(limiterType);
    GlobalDataBase::UpdateData("uns_limiter", &limiterType, PHINT, 1);

    if (viscousType == INVISCID) return;

    string uns_vis_name = "default";
    GlobalDataBase::GetData("uns_vis_name", &uns_vis_name, PHSTRING, 1);

    int uns_vis;
    if (uns_vis_name.substr(0,4) == "test")
    {
        uns_vis = VIS_TEST;
    }
    else if (uns_vis_name.substr(0,4) == "aver")
    {
        uns_vis = VIS_AVER;
    }
    else if (uns_vis_name.substr(0,3) == "std")
    {
        uns_vis = VIS_STD;
    }
    else if (uns_vis_name.substr(0,4) == "new1")
    {
        uns_vis = VIS_NEW1;
    }
    else if (uns_vis_name.substr(0,4) == "new2")
    {
        uns_vis = VIS_NEW2;
    }
    else if (uns_vis_name.substr(0,4) == "eric")
    {
        uns_vis = VIS_ERIC;
    }
    else
    {
        uns_vis = VIS_STD;
    }

    GlobalDataBase::UpdateData("uns_vis", &uns_vis, PHINT, 1);
}

#ifdef USE_GMRESSOLVER
void NSSolverUnstruct::FindNeighborOfRoot(int root, int& neighbor, int nTotalFace, int nBoundFace, int* leftCellofFace, int* rightCellofFace, vector<int>& AI, vector<int>& AJ, int counter)
{
    if(counter == 0)
    {
        return;
    }

    for (int iFace = 0; iFace < nTotalFace; iFace++)
    {
        if (leftCellofFace[iFace] == neighbor)
        {
            int newcell = rightCellofFace[iFace];
            if (find(AJ.begin()+AI[root], AJ.begin()+AI[root + 1], newcell) == AJ.begin()+AI[root + 1])
            {
                AJ.push_back(newcell);
                AI[root + 1]++;
                FindNeighborOfRoot(root, newcell, nTotalFace, nBoundFace, leftCellofFace, rightCellofFace, AI, AJ, counter - 1);
            }
        }
        if (rightCellofFace[iFace] == neighbor)
        {
            int newcell = leftCellofFace[iFace];
            if (find(AJ.begin()+AI[root], AJ.begin()+AI[root + 1], newcell) == AJ.begin()+AI[root + 1])
            {
                AJ.push_back(newcell);
                AI[root + 1]++;
                FindNeighborOfRoot(root, newcell, nTotalFace, nBoundFace, leftCellofFace, rightCellofFace, AI, AJ, counter - 1);
            }
        }
    }
}

void NSSolverUnstruct::AllocateJacobianMatrix4GMRES(Grid *gridIn)
{
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");

    //! get limiter order 
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int limiterType = parameters->GetLimiterType();

    //! get flow type: inviscid or viscous
    int viscousType = parameters->GetViscousType();

    //! determine the order required for the Jacobian matrix components
    int order = 1;
    int order_inv = 1;
    int order_vis = 1;
    if (tscheme == GMRES || tscheme == LU_SGS)
    {
        //! inviscid part
        if(limiterType == ILMT_FIRST)
        {
            order_inv = 1;
        }
        else if (limiterType == ILMT_NOLIM)
        {
            order_inv = 2;
        }
        else
        {
            printf("GMRES solver: uns_limiter_name is fixed to [1st] and [nolim] so far\n");
            TK_Exit::ExceptionExit("If you implemented other methods, please change it in subroutine [FindJacobianIndex4GMRES]");
        }

        //! viscous part
        if(viscousType!=INVISCID)
        {
            string gradientName = parameters->GetGradientName();
            if(gradientName == "ggcell")
            {
                order_vis = 2;
            }
            else 
            {
                printf("GMRES solver: for NS/RANS, ggcell only so far\n");
                TK_Exit::ExceptionExit("If you implemented other methods, please change it in subroutine [FindJacobianIndex4GMRES]");
            }
        }
        else
        {
            order_vis = 1;
        }
        order = MAX(order_inv, order_vis);
    }

    //! get grid information
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int nEquation = GetNumberOfEquations();
    int nTotal = nTotalCell + nBoundFace;
    int nInterFace = nTotalFace - nBoundFace;
    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    vector<int> AI;
    vector<int> AJ;
    AI.insert(AI.begin(), nTotal + 1, 0);
    for (int iCell = 0; iCell < nTotal; iCell++)
    {
        AI[iCell + 1]++;
        AJ.push_back(iCell);
        FindNeighborOfRoot(iCell, iCell, nTotalFace, nBoundFace, leftCellofFace, rightCellofFace, AI, AJ, order);
        if(iCell < nTotal - 1) AI[iCell + 2] = AI[iCell + 1];
    }

    int TotalSize = AI[nTotal];
    RDouble **dRdq = NewPointer2<RDouble>(nEquation, nEquation * TotalSize);
    grid->UpdateDataPtr("dRdq", dRdq);
    PHSPACE::SetField(dRdq, nEquation, nEquation * TotalSize, 0.0);

    grid->SetJacobianAI4GMRES(AI);
    grid->SetJacobianAJ4GMRES(AJ);
}

void NSSolverUnstruct::AllocateJacobianMatrix4GMRES_Neighbor(Grid *gridIn)
{
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");

    //! get limiter order 
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int limiterType = parameters->GetLimiterType();

    //! get flow type: inviscid or viscous
    int viscousType = parameters->GetViscousType();

    //! determine the order required for the Jacobian matrix components
    int order = 1;
    int order_inv = 1;
    int order_vis = 1;
    //if (tscheme == GMRES)
    //{
    //! inviscid part
    if(limiterType == ILMT_FIRST)
    {
        order_inv = 1;
    }
    else if (limiterType == ILMT_NOLIM || limiterType == ILMT_VENCAT)
    {
        order_inv = 2;
    }
    else
    {
        printf("GMRES solver: uns_limiter_name is fixed to [1st] and [nolim] so far\n");
        TK_Exit::ExceptionExit("If you implemented other methods, please change it in subroutine [FindJacobianIndex4GMRES]");
    }
    //! viscous part
    if(viscousType!=INVISCID)
    {
        string gradientName = parameters->GetGradientName();
        if(gradientName == "ggcell")
        {
            order_vis = 2;
        }
        else 
        {
            printf("GMRES solver: for NS/RANS, ggcell only so far\n");
            TK_Exit::ExceptionExit("If you implemented other methods, please change it in subroutine [FindJacobianIndex4GMRES]");
        }
    }
    else
    {
        order_vis = 1;
    }
    order = MAX(order_inv, order_vis);
    //}
    //! get grid information
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int nEquation = GetNumberOfEquations();
    int nTotal = nTotalCell + nBoundFace;
    int nInterFace = nTotalFace - nBoundFace;
    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    // compute the neighbor information
    grid->ComputeNeighborinfo();
    vector<int>* neighborCells = grid->GMRESGetNeighborCells();
    //! GMRESParallel
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if( interfaceInfo != nullptr)
    {
        int* localInterfaceCellIndex = interfaceInfo->GetLocalInterfaceCellIndex();
    }
    //! GMRES CSR Format Jacobian matrix with respect to primitive variables
    //! #AI:The i-th element records the number of non-zero elements contained in the first i-1 row
    //! #AJ:The i-th element records the index of columns of the V [i] element, which V [i] is used to store the values of non zero elements in a matrix;
    vector<int> AI;
    vector<int> AJ;
    //! GMRESParallel
    vector<int> AK;    //! stores the number of ghost cell for interface bc at each row

    vector<int> AI1st;    //! GMRESJac1st  
    vector<int> AJ1st;    //! GMRESJac1st

    AI.insert(AI.begin(), nTotal + 1, 0);
    AK.insert(AK.begin(), nTotal, 0); 
    if (order == 2)
    {
        AI1st.insert(AI1st.begin(), nTotal + 1, 0);    //! GMRESJac1st
    }

    for (int iCell = 0; iCell < nTotal; iCell++)
    {
        AI[iCell + 1]++;
        AJ.push_back(iCell);

        //! GMRESJac1st
        if(order == 2)
        {
            AI1st[iCell + 1]++;
            AJ1st.push_back(iCell);
        }

        for(int index = 0; index < neighborCells[iCell].size(); index++)
        {
            int neighborIndex = neighborCells[iCell][index];
            AI[iCell + 1]++;
            AJ.push_back(neighborIndex);
            //! GMRESParallel
            if(neighborIndex >= nTotalCell && PHMPI::GetNumberofGlobalZones() != 1)
            {
                if(interfaceInfo->MatchedGlobalNeighborCellIndex(neighborIndex) != -1)
                {
                    AK[iCell]++;
                    //printf("AK %d iCell %d\n", AK[iCell], iCell);
                }
            }
            if(order == 2)
            {
                //! GMRESJac1st
                AI1st[iCell + 1]++;
                AJ1st.push_back(neighborIndex);
                
                for(int indexJ = 0; indexJ < neighborCells[neighborIndex].size(); indexJ++)
                {
                    int neighborIndexJ = neighborCells[neighborIndex][indexJ];
                    if(neighborIndexJ != iCell)
                    {
                        vector<int>::iterator result = find(AJ.begin()+AI[iCell], AJ.begin()+AI[iCell + 1], neighborIndexJ); // find qcell
                        if(result==AJ.begin()+AI[iCell+1])
                        {
                            AI[iCell + 1]++;
                            AJ.push_back(neighborIndexJ);
                        }
                    }
                }
            }
        }

        if(iCell < nTotal - 1) AI[iCell + 2]    = AI[iCell + 1];

        if( order == 2)
        {
            if(iCell < nTotal - 1) AI1st[iCell + 2] = AI1st[iCell + 1];    //! GMRESJac1st
        }
    }

     //! set AK
    if (interfaceInfo != nullptr && PHMPI::GetCurrentProcessorID() == 0)
    {
        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            //printf("!cell %d, expectValue  ", iCell + 1);
            int expectValue = 0;
            int rcell = nTotalCell;
            do
            {
                vector<int>::iterator result = find(AJ.begin() + AI[rcell], AJ.begin() + AI[rcell + 1], iCell); // find cell
                if(result != AJ.begin() + AI[rcell + 1])
                {
                    expectValue++;
                    int index = distance(AJ.begin(), result);
                    //printf("%d ", interfaceInfo->MatchedGlobalNeighborCellIndex(rcell) + 1);
                    // printf("%d ", AJ[index] + 1);
                }
                rcell++;
            } while (rcell < nTotal);
            //printf("\n");
            AK[iCell] = expectValue;
        }
    } 
    // TK_Exit::ExitPHengLEI();
 
    //! GMRESCoupled
    int blockSizeOfJacobianMatrix = nEquation;
    if(viscousType == ONE_EQU)
    {
        blockSizeOfJacobianMatrix++;
    }
    grid->SetBlockSizeOfJacobianMatrix(blockSizeOfJacobianMatrix);
    
    int TotalSize = AI[nTotal];
    RDouble **dRdq = NewPointer2<RDouble>(nEquation, nEquation * TotalSize);
    grid->UpdateDataPtr("dRdq", dRdq);
    PHSPACE::SetField(dRdq, nEquation, nEquation * TotalSize, 0.0);

    //! GMRESJac1st
    if( order == 2 )
    {
        int TotalSize1st    = AI1st[nTotal];
        RDouble **dRdq1st   = NewPointer2<RDouble>(nEquation, nEquation * TotalSize1st);
        grid->UpdateDataPtr("dRdq1st", dRdq1st);
        PHSPACE::SetField(dRdq1st, nEquation, nEquation * TotalSize1st, 0.0);
    }
    else
    {
        RDouble** dRdq1st = NewPointer2<RDouble>(1,1);
        grid->UpdateDataPtr("dRdq1st",dRdq1st);
    }

    grid->SetJacobianAI1st4GMRES(AI1st);
    grid->SetJacobianAJ1st4GMRES(AJ1st);
    grid->SetJacobianOrder(order);

    grid->SetJacobianAI4GMRES(AI);
    grid->SetJacobianAJ4GMRES(AJ);
    grid->SetJacobianAK4GMRES(AK);

    //! GMRESCoupled
    if( viscousType == ONE_EQU )
    {
        int TotalSize = AI[nTotal];
        RDouble **dRdqCoupledTerm = NewPointer2<RDouble>(nEquation, TotalSize);
        grid->UpdateDataPtr("dRdqCoupledTerm", dRdqCoupledTerm);
        PHSPACE::SetField(dRdqCoupledTerm, nEquation, TotalSize, 0.0);
    }

    RDouble **dRdqCoupledTerm_turb = NewPointer2<RDouble>(1, nEquation * TotalSize);
    grid->UpdateDataPtr("dRdqCoupledTerm_turb", dRdqCoupledTerm_turb);
    PHSPACE::SetField(dRdqCoupledTerm_turb, 1, nEquation * TotalSize, 0.0);

    //! GMRESCoupled1st
    if( viscousType == ONE_EQU && order == 2 )
    {
        int TotalSize1st    = AI1st[nTotal];
        RDouble **dRdqCoupledTerm1st = NewPointer2<RDouble>(nEquation, TotalSize1st);
        grid->UpdateDataPtr("dRdqCoupledTerm1st", dRdqCoupledTerm1st);
        PHSPACE::SetField(dRdqCoupledTerm1st, nEquation, TotalSize1st, 0.0);
    }
    else
    {
        RDouble **dRdqCoupledTerm1st = NewPointer2<RDouble>(1, 1);
        grid->UpdateDataPtr("dRdqCoupledTerm1st", dRdqCoupledTerm1st);
    }

    if( order==2 )
    {
        int TotalSize1st    = AI1st[nTotal];
        RDouble **dRdqCoupledTerm_turb1st = NewPointer2<RDouble>(1, nEquation*TotalSize1st);
        grid->UpdateDataPtr("dRdqCoupledTerm_turb1st", dRdqCoupledTerm_turb1st);
        PHSPACE::SetField(dRdqCoupledTerm_turb1st, 1, nEquation*TotalSize1st, 0.0);
    }
    else
    {
        RDouble **dRdqCoupledTerm_turb1st = NewPointer2<RDouble>(1, 1);
        grid->UpdateDataPtr("dRdqCoupledTerm_turb1st", dRdqCoupledTerm_turb1st);
    }

}
#endif

void NSSolverUnstruct::AllocateGlobalVar(Grid *gridIn)
{
    CFDSolver::AllocateGlobalVar(gridIn);

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **q   = NewPointer2<RDouble>(nEquation, nTotal);
    RDouble **res = NewPointer2<RDouble>(nEquation, nTotal);
    RDouble *diagonal = new RDouble [nTotal];

    grid->UpdateDataPtr("q", q);    //! Scalar flow field variable(rho/u/v/w/p).
    grid->UpdateDataPtr("res", res);    //! Residual or right-hand side.
    grid->UpdateDataPtr("diagonal", diagonal);    //! Sacrificing space to improve efficiency.

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");
    int nTurboZone = 0;
    if (referenceFrame)
    {
        nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
    }

#ifdef USE_GMRESSOLVER

    int originaltscheme = GlobalDataBase::GetIntParaFromDB("OriginalTscheme");
    if( originaltscheme == GMRES )
    {
        //! GMRES CSR Format Jacobian matrix with respect to primitive variables
        AllocateJacobianMatrix4GMRES_Neighbor(gridIn);
    }

    //! GMRESBoundary d(Dummy)/d(Physical)
    RDouble **dDdP = NewPointer2<RDouble>(nEquation, nEquation*nBoundFace);
    grid->UpdateDataPtr("dDdP", dDdP);
    PHSPACE::SetField(dDdP, nEquation, nEquation*nBoundFace, 0.0);
#endif

    RDouble **limit2D = NewPointer2<RDouble>(nEquation, nTotal);
    grid->UpdateDataPtr("limit2D", limit2D);
    PHSPACE::SetField(limit2D, nEquation, nTotal, zero);

    int nTotalNode = grid->GetNTotalNode();
    RDouble **qnode = NewPointer2<RDouble>(nEquation, nTotalNode);
    grid->UpdateDataPtr("qnode", qnode);    //! Scalar flow field variable(rho/u/v/w/p) at node.

    RDouble **tnode = NewPointer2<RDouble>(nTemperatureModel, nTotalNode);
    grid->UpdateDataPtr("tnode", tnode);    //! Static temperature at node.

    RDouble *nodeWeight = new RDouble[nTotalNode];
    grid->UpdateDataPtr("nodeWeight", nodeWeight);    //! Weight coefficient at node.

    RDouble *nodeWeightTrade = new RDouble[nTotalNode];
    PHSPACE::SetField(nodeWeightTrade, 0.0, nTotalNode);
    grid->UpdateDataPtr("nodeWeightTrade", nodeWeightTrade);

    if (grid->GetLevel() == 0)
    {
        int *nodeBCType = new int[nTotalNode];
        PHSPACE::SetField(nodeBCType, 0, nTotalNode);
        grid->UpdateDataPtr("nodeBCType", nodeBCType);
    }

    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation)
    {
        int numberOfInterpoints = interpointInformation->GetNumberOfInterpoints();
        RDouble **qInterPoint = NewPointer2<RDouble>(nEquation, numberOfInterpoints);
        grid->UpdateDataPtr("qInterPoint", qInterPoint);    //! Scalar flow field variable(rho/u/v/w/p) at interPoint.

        //RDouble **tInterPoint = NewPointer2<RDouble>(1, numberOfInterpoints);
        RDouble **tInterPoint = NewPointer2<RDouble>(nTemperatureModel, numberOfInterpoints);
        grid->UpdateDataPtr("tInterPoint", tInterPoint);    //! Static temperature at interPoint.
    }

    RDouble *rtem = new RDouble [nTotal];
    grid->UpdateDataPtr("rtem", rtem);    //! Pressure factor.

    RDouble *invSpectralRadius = new RDouble [nTotalCell];
    grid->UpdateDataPtr("invSpectralRadius", invSpectralRadius);    //! Inviscid spectral radius.

    RDouble *visSpectralRadius = new RDouble [nTotalCell];
    grid->UpdateDataPtr("visSpectralRadius", visSpectralRadius);    //! Viscous spectral radius.

    int nSpeciesNumber = parameters->GetNumberOfSpecies();
    int nChemicalRadius = parameters->GetNChemicalRadius();

    if (nChemical == 1 )
    {
        if(nChemicalRadius == 1)
        {
            RDouble **chemSpectralRadius = NewPointer2<RDouble>(nSpeciesNumber,nTotalCell);
            grid->UpdateDataPtr("chemSpectralRadius"  , chemSpectralRadius);    //! Chemical reaction source item.
        }
    }

    int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

    if (ifLowSpeedPrecon != 0)
    {
        RDouble *preconCoefficient = new RDouble [nTotal];
        grid->UpdateDataPtr("preconCoefficient", preconCoefficient);
        PHSPACE::SetField(preconCoefficient, 1.0, nTotal);
        RDouble *timeCoefficient = new RDouble [nTotal];
        grid->UpdateDataPtr("timeCoefficient", timeCoefficient);
        PHSPACE::SetField(timeCoefficient, 1.0, nTotal);
        RDouble *timeCoefficientInverse = new RDouble [nTotal];
        grid->UpdateDataPtr("timeCoefficientInverse", timeCoefficientInverse);
        PHSPACE::SetField(timeCoefficientInverse, 1.0, nTotal);
        int nElement = nEquation * nEquation;
        RDouble **preconMatrix = NewPointer2 <RDouble> (nElement, nTotal);
        grid->UpdateDataPtr("preconMatrix", preconMatrix);
    }

    int isUnsteady = parameters->GetIsUnsteady();

    if (isUnsteady)
    {
        RDouble **q_unsteady_n1   = NewPointer2<RDouble>(nEquation, nTotal);
        RDouble **q_unsteady_n2   = NewPointer2<RDouble>(nEquation, nTotal);
        RDouble **res_unsteady_n1 = NewPointer2<RDouble>(nEquation, nTotal);
        RDouble **res_unsteady_n2 = NewPointer2<RDouble>(nEquation, nTotal);

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble **res_unsteady_tmp = NewPointer2<RDouble>(nEquation, nTotal);

        grid->UpdateDataPtr("q_unsteady_n1"  , q_unsteady_n1 );
        grid->UpdateDataPtr("q_unsteady_n2"  , q_unsteady_n2 );
        grid->UpdateDataPtr("res_unsteady_n1", res_unsteady_n1);
        grid->UpdateDataPtr("res_unsteady_n2", res_unsteady_n2);
        grid->UpdateDataPtr("res_unsteady_tmp", res_unsteady_tmp);

        //! Statistical variables for unsteady simulation.
        RDouble **qAverage = NewPointer2<RDouble>(nEquation, nTotal); 
        grid->UpdateDataPtr("qAverage", qAverage);
        PHSPACE::SetField(qAverage, nEquation, nTotal, zero);

        RDouble **tauAverage = NewPointer2<RDouble>(6, nTotal); 
        grid->UpdateDataPtr("tauAverage", tauAverage);
        PHSPACE::SetField(tauAverage, 6, nTotal, zero);

        RDouble **q2Average = NewPointer2<RDouble>(6, nTotal); 
        grid->UpdateDataPtr("q2Average", q2Average);
        PHSPACE::SetField(q2Average, 6, nTotal, zero);
    }

    RDouble **t = NewPointer2<RDouble>(nTemperatureModel, nTotal);
    grid->UpdateDataPtr("t", t);    //! Static temperature.

    RDouble *dt  = new RDouble [nTotal];
    grid->UpdateDataPtr("dt", dt);    //! Time step.

    RDouble *gama  = new RDouble [nTotal];
    grid->UpdateDataPtr("gama", gama);    //! Ratio of specific heat coefficients at constant pressure and volume.

    bool isViscous = parameters->IsViscous();

    if (isViscous)
    {
        RDouble *visl = new RDouble [nTotal];
        RDouble *vist = new RDouble [nTotal];

        PHSPACE::SetField(visl, 1.0, nTotal);
        PHSPACE::SetField(vist, 0.0, nTotal);

        grid->UpdateDataPtr("visl", visl);    //! Laminar viscous coefficient.
        grid->UpdateDataPtr("vist", vist);    //! Turbulence viscous coefficient.
    }

    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
    if (iLES == LES_SOLVER)
    {
        RDouble *subgridScaleEnergy = new RDouble[nTotal];
        PHSPACE::SetField(subgridScaleEnergy, 0.0, nTotal);
        grid->UpdateDataPtr("subgridScaleEnergy", subgridScaleEnergy);
    }

    RDouble **gradPrimtiveVarX = NewPointer2<RDouble>(nEquation, nTotal);
    RDouble **gradPrimtiveVarY = NewPointer2<RDouble>(nEquation, nTotal);
    RDouble **gradPrimtiveVarZ = NewPointer2<RDouble>(nEquation, nTotal);
    PHSPACE::SetField(gradPrimtiveVarX, nEquation, nTotal, zero);
    PHSPACE::SetField(gradPrimtiveVarY, nEquation, nTotal, zero);
    PHSPACE::SetField(gradPrimtiveVarZ, nEquation, nTotal, zero);

    grid->UpdateDataPtr("gradPrimtiveVarX", gradPrimtiveVarX);    //! Gradient of scalar flow field variable at cell, for x direction.
    grid->UpdateDataPtr("gradPrimtiveVarY", gradPrimtiveVarY);    //! Gradient of scalar flow field variable at cell, for y direction.
    grid->UpdateDataPtr("gradPrimtiveVarZ", gradPrimtiveVarZ);    //! Gradient of scalar flow field variable at cell, for z direction.

    RDouble **gradTemperatureX = NewPointer2<RDouble>(nTemperatureModel, nTotal);
    RDouble **gradTemperatureY = NewPointer2<RDouble>(nTemperatureModel, nTotal);
    RDouble **gradTemperatureZ = NewPointer2<RDouble>(nTemperatureModel, nTotal);
    PHSPACE::SetField(gradTemperatureX, nTemperatureModel, nTotal, zero);
    PHSPACE::SetField(gradTemperatureY, nTemperatureModel, nTotal, zero);
    PHSPACE::SetField(gradTemperatureZ, nTemperatureModel, nTotal, zero);
    grid->UpdateDataPtr("gradTemperatureX", gradTemperatureX);    //! Gradient of static temperature at cell, for x direction.
    grid->UpdateDataPtr("gradTemperatureY", gradTemperatureY);    //! Gradient of static temperature at cell, for y direction.
    grid->UpdateDataPtr("gradTemperatureZ", gradTemperatureZ);    //! Gradient of static temperature at cell, for z direction.

    //! Temporary variable for periodic boundary condition.
    RDouble **rotNSgradValueX = NewPointer2<RDouble>(nEquation, nTotal);
    RDouble **rotNSgradValueY = NewPointer2<RDouble>(nEquation, nTotal);
    RDouble **rotNSgradValueZ = NewPointer2<RDouble>(nEquation, nTotal);
    PHSPACE::SetField(rotNSgradValueX, nEquation, nTotal, zero);
    PHSPACE::SetField(rotNSgradValueY, nEquation, nTotal, zero);
    PHSPACE::SetField(rotNSgradValueZ, nEquation, nTotal, zero);
    grid->UpdateDataPtr("rotNSgradValueX", rotNSgradValueX);
    grid->UpdateDataPtr("rotNSgradValueY", rotNSgradValueY);
    grid->UpdateDataPtr("rotNSgradValueZ", rotNSgradValueZ);

    //! allocate mixing plane variables.
    //! for multi-row turbomachinery.
    if (nTurboZone > 0)
    {
        int nSpanSection = PHSPACE::GlobalDataBase::GetIntParaFromDB("nSpanSection");
        int nType = 2;
        RDouble **q_Face = NewPointer2<RDouble>(nEquation, nTotal);
        RDouble **qCon_Face = NewPointer2<RDouble>(nEquation, nTotal);
        PHSPACE::SetField(q_Face, nEquation, nTotal, zero);
        PHSPACE::SetField(qCon_Face, nEquation, nTotal, zero);
        grid->UpdateDataPtr("q_Face", q_Face);
        grid->UpdateDataPtr("qCon_Face", qCon_Face);

        //! zone with MixingIn and MixingOut.
        int nFlag = 10 * nTurboZone + 2;
        RDouble **Radius_Face = NewPointer2<RDouble>(nFlag, nTotal);
        RDouble **Radius_Span = NewPointer2<RDouble>(nFlag, nSpanSection);
        PHSPACE::SetField(Radius_Face, nFlag, nTotal, zero);
        PHSPACE::SetField(Radius_Span, nFlag, nSpanSection, zero);
        grid->UpdateDataPtr("Radius_Face", Radius_Face);
        grid->UpdateDataPtr("Radius_Span", Radius_Span);

        RDouble ***SpanFlux      = NewPointer3<RDouble>(nFlag, nEquation, nSpanSection);
        RDouble ***q_Span        = NewPointer3<RDouble>(nFlag, nEquation, nSpanSection);
        RDouble ***qCon_Span     = NewPointer3<RDouble>(nFlag, nEquation, nSpanSection);
        RDouble ***qSpanNeighbor = NewPointer3<RDouble>(nFlag, nEquation, nSpanSection);
        RDouble ***dqSpanIn      = NewPointer3<RDouble>(nFlag, nEquation, nSpanSection);
        RDouble ***dqSpanEx      = NewPointer3<RDouble>(nFlag, nEquation, nSpanSection);
        RDouble ***dcSpan        = NewPointer3<RDouble>(nFlag, nEquation, nSpanSection);
        InitMixingPlane(SpanFlux, nFlag, nEquation, nSpanSection, zero);
        InitMixingPlane(q_Span, nFlag, nEquation, nSpanSection, zero);
        InitMixingPlane(qCon_Span, nFlag, nEquation, nSpanSection, zero);
        InitMixingPlane(qSpanNeighbor, nFlag, nEquation, nSpanSection, zero);
        InitMixingPlane(dqSpanIn, nFlag, nEquation, nSpanSection, zero);
        InitMixingPlane(dqSpanEx, nFlag, nEquation, nSpanSection, zero);
        InitMixingPlane(dcSpan, nFlag, nEquation, nSpanSection, zero);
        grid->UpdateDataPtr("SpanFlux", SpanFlux);
        grid->UpdateDataPtr("q_Span", q_Span);
        grid->UpdateDataPtr("qCon_Span", qCon_Span);
        grid->UpdateDataPtr("qSpanNeighbor", qSpanNeighbor);
        grid->UpdateDataPtr("dqSpanIn", dqSpanIn);
        grid->UpdateDataPtr("dqSpanEx", dqSpanEx);
        grid->UpdateDataPtr("dcSpan", dcSpan);

        RDouble **SpanArea = NewPointer2<RDouble>(nFlag, nSpanSection);
        PHSPACE::SetField(SpanArea, nFlag, nSpanSection, zero);
        grid->UpdateDataPtr("SpanArea", SpanArea);

        RDouble **nxsSpan = NewPointer2<RDouble>(nFlag, nSpanSection);
        RDouble **ntsSpan = NewPointer2<RDouble>(nFlag, nSpanSection);
        RDouble **nrsSpan = NewPointer2<RDouble>(nFlag, nSpanSection);
        PHSPACE::SetField(nxsSpan, nFlag, nSpanSection, zero);
        PHSPACE::SetField(ntsSpan, nFlag, nSpanSection, zero);
        PHSPACE::SetField(nrsSpan, nFlag, nSpanSection, zero);
        grid->UpdateDataPtr("nxsSpan", nxsSpan);
        grid->UpdateDataPtr("ntsSpan", ntsSpan);
        grid->UpdateDataPtr("nrsSpan", nrsSpan);
    }

    systemgridtype = GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if (MIXGRID == systemgridtype)
    {
        InitiallimiterOnlyForMixGrid(grid);

        InitialInviscidfaceproxyOnlyForMixGrid(grid);

        if (isViscous)
        {
            InitialViscousfaceproxyOnlyForMixGrid(grid);
        }
    }
    //DelPointer2(res);
    //DelPointer2(q);
    int level = grid->GetLevel();
    int zoneID = grid->GetZoneLocalID();
    if (level ==0 && zoneID == 0)
    {
        RDouble *monitorVariables = new RDouble [30];
        PHSPACE::SetField(monitorVariables, 0.0, 30);
        GlobalDataBase::UpdateDataPtr("monitorVariables", monitorVariables);
    }
}

void NSSolverUnstruct::DeAllocateGlobalVar(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nChemical = parameters->GetChemicalFlag();

    int nChemicalRadius = parameters->GetNChemicalRadius();

    RDouble **q    = reinterpret_cast <RDouble **> (grid->GetDataPtr("q"));
    RDouble **res  = reinterpret_cast <RDouble **> (grid->GetDataPtr("res"));
    RDouble *diagonal = reinterpret_cast <RDouble *> (grid->GetDataPtr("diagonal"));

    DelPointer2(q);
    DelPointer2(res);
    delete [] diagonal;    diagonal = nullptr;

    RDouble **limit2D = reinterpret_cast <RDouble **> (grid->GetDataPtr("limit2D"));
    DelPointer2(limit2D);

    RDouble **tnode = reinterpret_cast <RDouble **> (grid->GetDataPtr("tnode"));
    RDouble **qnode = reinterpret_cast <RDouble **> (grid->GetDataPtr("qnode"));
    DelPointer2(tnode);
    DelPointer2(qnode);

    RDouble *nodeWeight = reinterpret_cast <RDouble *> (grid->GetDataPtr("nodeWeight"));
    RDouble *nodeWeightTrade = reinterpret_cast <RDouble *> (grid->GetDataPtr("nodeWeightTrade"));
    delete [] nodeWeight;    nodeWeight = nullptr;
    delete [] nodeWeightTrade;    nodeWeightTrade = nullptr;

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");
    int nTurboZone = 0;
    if (referenceFrame)
    {
        nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
    }

#ifdef USE_GMRESSOLVER
    //! GMRESBoundary
    RDouble **dDdP = reinterpret_cast <RDouble**> (grid->GetDataPtr("dDdP"));
    DelPointer2(dDdP);

    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");

    if( tscheme == GMRES )
    {
        //! GMRES
        RDouble **dRdq = reinterpret_cast <RDouble**> (grid->GetDataPtr("dRdq"));
        //! GMRESJac1st
        RDouble **dRdq1st = reinterpret_cast <RDouble**> (grid->GetDataPtr("dRdq1st"));
        DelPointer2(dRdq);
        DelPointer2(dRdq1st);

        //! GMRESCoupled
        int GviscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
        if( GviscousType == ONE_EQU )
        {
            RDouble **dRdqCoupledTerm = reinterpret_cast <RDouble**> (grid->GetDataPtr("dRdqCoupledTerm"));

            DelPointer2(dRdqCoupledTerm);
        }

        RDouble **dRdqCoupledTerm_turb      = reinterpret_cast <RDouble**> (grid->GetDataPtr("dRdqCoupledTerm_turb"));
        RDouble **dRdqCoupledTerm1st        = reinterpret_cast <RDouble**> (grid->GetDataPtr("dRdqCoupledTerm1st"));
        RDouble **dRdqCoupledTerm_turb1st   = reinterpret_cast <RDouble**> (grid->GetDataPtr("dRdqCoupledTerm_turb1st"));
        DelPointer2(dRdqCoupledTerm_turb);
        DelPointer2(dRdqCoupledTerm1st);
        DelPointer2(dRdqCoupledTerm_turb1st);
    }
#endif
    if (grid->GetLevel() == 0)
    {
        int *nodeBCType = reinterpret_cast <int *> (grid->GetDataPtr("nodeBCType"));
        delete [] nodeBCType;    nodeBCType = nullptr;
    }

    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation)
    {
        RDouble **qInterPoint = reinterpret_cast <RDouble **> (grid->GetDataPtr("qInterPoint"));
        RDouble **tInterPoint = reinterpret_cast <RDouble **> (grid->GetDataPtr("tInterPoint"));
        DelPointer2(qInterPoint);
        DelPointer2(tInterPoint);
    }

    RDouble *rtem = reinterpret_cast <RDouble *> (grid->GetDataPtr("rtem"));
    delete [] rtem;    rtem = nullptr;
    RDouble *invSpectralRadius = reinterpret_cast <RDouble *> (grid->GetDataPtr("invSpectralRadius"));
    delete [] invSpectralRadius;    invSpectralRadius = nullptr;
    RDouble *visSpectralRadius = reinterpret_cast <RDouble *> (grid->GetDataPtr("visSpectralRadius"));
    delete [] visSpectralRadius;    visSpectralRadius = nullptr;

    RDouble* DiagValue = reinterpret_cast <RDouble*> (grid->GetDataPtr("DiagValue"));
    if (DiagValue != nullptr)
    {
        delete [] DiagValue;    DiagValue = nullptr;
    }

    RDouble ****DMatrixStartByLine = reinterpret_cast<RDouble****>(grid->GetDataPtr("DMatrixStartByLine"));
    RDouble ****UMatrixStartByLine = reinterpret_cast<RDouble****>(grid->GetDataPtr("UMatrixStartByLine"));
    RDouble ****LMatrixStartByLine = reinterpret_cast<RDouble****>(grid->GetDataPtr("LMatrixStartByLine"));

    if (DMatrixStartByLine != nullptr)
    {
        delete [] DMatrixStartByLine;    DMatrixStartByLine = nullptr;
        delete [] UMatrixStartByLine;    UMatrixStartByLine = nullptr;
        delete [] LMatrixStartByLine;    LMatrixStartByLine = nullptr;
    }

    MAT3D< RDouble > *diagMatrix = reinterpret_cast<MAT3D< RDouble > *>(grid->GetDataPtr("diagMatrix"));
    MAT3D< RDouble > *lowerMatrix = reinterpret_cast<MAT3D< RDouble > *>(grid->GetDataPtr("lowerMatrix"));
    MAT3D< RDouble > *upperMatrix = reinterpret_cast<MAT3D< RDouble > *>(grid->GetDataPtr("upperMatrix"));

    if (diagMatrix != nullptr)
    {
        delete diagMatrix;    diagMatrix = nullptr;
        delete upperMatrix;    upperMatrix = nullptr;
        delete lowerMatrix;    lowerMatrix = nullptr;
    }

    if (nChemical == 1)
    {
        if(nChemicalRadius == 1)
        {
            RDouble **chemSpectralRadius = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("chemSpectralRadius") );
            DelPointer2(chemSpectralRadius);
            //grid->DeleteDataPtr("chemSpectralRadius");
        }
    }

    int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();
    if (ifLowSpeedPrecon != 0)
    {
        RDouble *preconCoefficient = reinterpret_cast< RDouble * > (grid->GetDataPtr("preconCoefficient"));
        RDouble *timeCoefficient = reinterpret_cast< RDouble * > (grid->GetDataPtr("timeCoefficient"));
        RDouble *timeCoefficientInverse = reinterpret_cast< RDouble * > (grid->GetDataPtr("timeCoefficientInverse"));
        RDouble **preconMatrix = reinterpret_cast< RDouble ** > (grid->GetDataPtr("preconMatrix"));
        delete preconCoefficient;
        delete timeCoefficient;
        delete timeCoefficientInverse;
        DelPointer2(preconMatrix);
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble **q_unsteady_n1   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n1" ));
        RDouble **q_unsteady_n2   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n2" ));
        RDouble **res_unsteady_n1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n1"));
        RDouble **res_unsteady_n2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n2"));

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble **res_unsteady_tmp = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_tmp"));

        DelPointer2(q_unsteady_n1 );
        DelPointer2(q_unsteady_n2 );
        DelPointer2(res_unsteady_n1);
        DelPointer2(res_unsteady_n2);
        DelPointer2(res_unsteady_tmp);

        //! Statistical variables for unsteady simulation.
        RDouble **qAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("qAverage"));
        DelPointer2(qAverage);

        //! Statistical variables for unsteady simulation.
        RDouble **tauAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("tauAverage"));
        DelPointer2(tauAverage);

        RDouble **q2Average = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q2Average"));
        DelPointer2(q2Average);
    }

    RDouble **t   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t" ));
    DelPointer2(t );

    RDouble * dt  = reinterpret_cast< RDouble * > (grid->GetDataPtr("dt" ));
    delete [] dt;    dt = nullptr;

    RDouble * gama  = reinterpret_cast< RDouble * > (grid->GetDataPtr("gama" ));
    delete [] gama;    gama = nullptr;

    int viscousType = parameters->GetViscousType();

    if (viscousType > INVISCID)
    {
        RDouble * visl = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));
        RDouble * vist = reinterpret_cast< RDouble * > (grid->GetDataPtr("vist"));
        delete [] visl;    visl = nullptr;
        delete [] vist;    vist = nullptr;
    }

    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
    if(iLES == LES_SOLVER)
    {
        RDouble * subgridScaleEnergy = reinterpret_cast< RDouble * > (grid->GetDataPtr("subgridScaleEnergy"));
        delete [] subgridScaleEnergy;    subgridScaleEnergy = nullptr;
    }

    //! deallocate mixing plane array.
    if (nTurboZone > 0)
    {
        RDouble **q_Face    = reinterpret_cast <RDouble **> (grid->GetDataPtr("q_Face"));
        RDouble **qCon_Face = reinterpret_cast <RDouble **> (grid->GetDataPtr("qCon_Face"));
        RDouble **SpanArea  = reinterpret_cast <RDouble **> (grid->GetDataPtr("SpanArea"));

        DelPointer2(q_Face);
        DelPointer2(qCon_Face);
        DelPointer2(SpanArea);

        RDouble ***SpanFlux      = reinterpret_cast <RDouble ***> (grid->GetDataPtr("SpanFlux"));
        RDouble ***q_Span        = reinterpret_cast <RDouble ***> (grid->GetDataPtr("q_Span"));
        RDouble ***qCon_Span     = reinterpret_cast <RDouble ***> (grid->GetDataPtr("qCon_Span"));
        RDouble ***qSpanNeighbor = reinterpret_cast <RDouble ***> (grid->GetDataPtr("qSpanNeighbor"));
        RDouble ***dqSpanIn      = reinterpret_cast <RDouble ***> (grid->GetDataPtr("dqSpanIn"));
        RDouble ***dqSpanEx      = reinterpret_cast <RDouble ***> (grid->GetDataPtr("dqSpanEx"));
        RDouble ***dcSpan        = reinterpret_cast <RDouble ***> (grid->GetDataPtr("dcSpan"));

        DelPointer3(SpanFlux);
        DelPointer3(q_Span);
        DelPointer3(qCon_Span);
        DelPointer3(qSpanNeighbor);
        DelPointer3(dqSpanIn);
        DelPointer3(dqSpanEx);
        DelPointer3(dcSpan);

        RDouble **Radius_Face = reinterpret_cast <RDouble **> (grid->GetDataPtr("Radius_Face"));
        RDouble **Radius_Span = reinterpret_cast <RDouble **> (grid->GetDataPtr("Radius_Span"));

        DelPointer2(Radius_Face);
        DelPointer2(Radius_Span);

        RDouble **nxsSpan = reinterpret_cast <RDouble **> (grid->GetDataPtr("nxsSpan"));
        RDouble **ntsSpan = reinterpret_cast <RDouble **> (grid->GetDataPtr("ntsSpan"));
        RDouble **nrsSpan = reinterpret_cast <RDouble **> (grid->GetDataPtr("nrsSpan"));

        DelPointer2(nxsSpan);
        DelPointer2(ntsSpan);
        DelPointer2(nrsSpan);
    }

    RDouble **gradPrimtiveVarX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradPrimtiveVarY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradPrimtiveVarZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    DelPointer2(gradPrimtiveVarX);
    DelPointer2(gradPrimtiveVarY);
    DelPointer2(gradPrimtiveVarZ);

    RDouble **gradTemperatureX = reinterpret_cast< RDouble ** > (grid->GetDataPtr("gradTemperatureX" ));
    RDouble **gradTemperatureY = reinterpret_cast< RDouble ** > (grid->GetDataPtr("gradTemperatureY" ));
    RDouble **gradTemperatureZ = reinterpret_cast< RDouble ** > (grid->GetDataPtr("gradTemperatureZ" ));

    DelPointer2(gradTemperatureX);
    DelPointer2(gradTemperatureY);
    DelPointer2(gradTemperatureZ);

    RDouble **rotNSgradValueX = reinterpret_cast< RDouble ** > (grid->GetDataPtr("rotNSgradValueX"));
    RDouble **rotNSgradValueY = reinterpret_cast< RDouble ** > (grid->GetDataPtr("rotNSgradValueY"));
    RDouble **rotNSgradValueZ = reinterpret_cast< RDouble ** > (grid->GetDataPtr("rotNSgradValueZ"));

    DelPointer2(rotNSgradValueX);
    DelPointer2(rotNSgradValueY);
    DelPointer2(rotNSgradValueZ);

    int level = grid->GetLevel();
    int zoneID = grid->GetZoneLocalID();
    if (level ==0 && zoneID == 0)
    {
        RDouble *monitorVariables = reinterpret_cast< RDouble * > (GlobalDataBase::GetDataPtr("monitorVariables" ));
        delete [] monitorVariables;    monitorVariables = nullptr;
    }
}

bool NSSolverUnstruct::JudgeIfRestart()
{
    string restartNSFile = ".\results\flow.dat";
    GlobalDataBase::GetData("restartNSFile", &restartNSFile, PHSTRING, 1);

    if(PHMPI::IsParallelRun())
    {
        restartNSFile = PHSPACE::AddSymbolToFileName(restartNSFile, "_", 0);
    }

    bool restart_flag = false;

    ifstream infile(restartNSFile.c_str(), ios::in);
    if (infile)
    {
        restart_flag = true;
        infile.close();
        infile.clear();
    }
    if(restart_flag)
    {
         int originaltscheme = GlobalDataBase::GetIntParaFromDB("OriginalTscheme");
         if(GMRES == originaltscheme)
         {
            ostringstream oss;
            oss << "The array \"OriginalTscheme\" is error. The GMRES solver cannot support restart calculation from flow.dat." << endl;
            TK_Exit::ExceptionExit(oss);
         }
    }
    return restart_flag;
}

bool NSSolverUnstruct::JudgeIfInterpolate()
{
    string restartGridFile = "";
    GlobalDataBase::GetData("restartGridFile", &restartGridFile, PHSTRING, 1);

    string restartNSVarFile = "";
    GlobalDataBase::GetData("restartNSVarFile", &restartNSVarFile, PHSTRING, 1);

    if(PHMPI::IsParallelRun())
    {
        restartNSVarFile = PHSPACE::AddSymbolToFileName(restartNSVarFile, "_", 0);
    }

    bool restart_flag = false;

    ifstream inGridFile(restartGridFile.c_str(), ios::in);
    ifstream inFlowFile(restartNSVarFile.c_str(), ios::in);
    if (inGridFile && inFlowFile)
    {
        restart_flag = true;
        inGridFile.close();
        inGridFile.clear();
        inFlowFile.close();
        inFlowFile.clear();
    }

    return restart_flag;
}

bool NSSolverUnstruct::JudgeIfReadAverage()
{
    string restartNSFile = ".\results\flow.dat";
    GlobalDataBase::GetData("restartNSFile", &restartNSFile, PHSTRING, 1);
    string averageFlowFile = AddSymbolToFileName(restartNSFile, "_", "Average");

    if (PHMPI::IsParallelRun())
    {
        averageFlowFile = PHSPACE::AddSymbolToFileName(averageFlowFile, "_", 0);
    }

    bool restart_flag = false;

    ifstream infile(averageFlowFile.c_str(), ios::in);
    if (infile)
    {
        restart_flag = true;
        infile.close();
        infile.clear();
    }

    return restart_flag;
}

void NSSolverUnstruct::DumpRestart(ActionKey *actkey)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int level = actkey->level;
    UnstructGrid * grid = UnstructGridCast(GetGrid(level));

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    int nEquation = GetNumberOfEquations();

    int outnstep = 0;
    GlobalDataBase::GetData("outnstep", &outnstep, PHINT, 1);

    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

    DataContainer * cdata = actkey->GetData();
    cdata->MoveToBegin();
    cdata->Write(&outnstep,sizeof(int));

    for (int m = 0; m < nEquation; ++ m)
    {
        cdata->Write(q[m],sizeof(RDouble) * nTotal);
    }

    int isUnsteady = parameters->GetIsUnsteady();

    if (!isUnsteady) return;

    RDouble physicalTime = 0;
    GlobalDataBase::GetData("physicalTime",&physicalTime, PHDOUBLE, 1);

    cdata->Write(&physicalTime, sizeof(RDouble));

    RDouble **qn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n1"));
    RDouble **qn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n2"));

    for (int m = 0; m < nEquation; ++ m)
    {
        cdata->Write(qn1[m], sizeof(RDouble) * nTotal);
    }

    for (int m = 0; m < nEquation; ++ m)
    {
        cdata->Write(qn2[m],sizeof(RDouble) * nTotal);
    }

    RDouble **res   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"           ));
    RDouble **resn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n1"));
    RDouble **resn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n2"));

    for (int m = 0; m < nEquation; ++ m)
    {
        cdata->Write(res[m], sizeof(RDouble) * nTotal);
    }

    for (int m = 0; m < nEquation; ++ m)
    {
        cdata->Write(resn1[m], sizeof(RDouble) * nTotal);
    }

    for (int m = 0; m < nEquation; ++ m)
    {
        cdata->Write(resn2[m],sizeof(RDouble) * nTotal);
    }
}

void NSSolverUnstruct::DumpRestartH5(ActionKey *actkey)
{
    using namespace PHMPI;
    int currentProcessorID = GetCurrentProcessorID();

    int version = 1;
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady  = parameters->GetIsUnsteady();
    int outputStepInterval = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (currentProcessorID == GetServerProcessorID())
    {
        WriteData(actkey->filepos, &version, "Version");
        WriteData(actkey->filepos, &outputStepInterval, "outnstep");

        if (isUnsteady == 1)
        {
            RDouble physicalTime   = static_cast<RDouble>(GlobalDataBase::GetDoubleParaFromDB("physicalTime"));
            WriteData(actkey->filepos, &physicalTime, "physicalTime");

            if (IsNeedStatistics())
            {
                int nStatisticalStep = 0;
                GlobalDataBase::GetData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
                WriteData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
            }
        }
    }

    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    int nEquation = GetNumberOfEquations();

    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));
    double **primitiveVariables = NewPointer2<double>(nEquation, nTotal);

    for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            primitiveVariables[iEquation][iCell] = static_cast<double>(q[iEquation][iCell]);
        }
    }

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss<<"Group"<<gridID;
    grpName = oss.str();

    grploc = OpenGroup(actkey->filepos, grpName);

    WriteData(grploc, primitiveVariables[0], "q");
    WriteData(grploc, t[ITT], "t");
    WriteData(grploc, &nTotalCell, "nTotalCell");

    if (isUnsteady == 1)
    {
        double **qn1 = reinterpret_cast< double ** > (grid->GetDataPtr("q_unsteady_n1"));
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = static_cast<double>(qn1[iEquation][iCell]);
            }
        }
        WriteData(grploc, primitiveVariables[0], "q_unsteady_n1");

        RDouble **qn2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_unsteady_n2"));
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = static_cast<double>(qn2[iEquation][iCell]);
            }
        }
        WriteData(grploc, primitiveVariables[0], "q_unsteady_n2");

        RDouble **res   = reinterpret_cast<RDouble **> (grid->GetDataPtr("res"));
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = static_cast<double>(res[iEquation][iCell]);
            }
        }
        WriteData(grploc, primitiveVariables[0], "res");

        RDouble **resn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n1"));
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = static_cast<double>(resn1[iEquation][iCell]);
            }
        }
        WriteData(grploc, primitiveVariables[0], "res_unsteady_n1");

        RDouble **resn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n2"));
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = static_cast<double>(resn2[iEquation][iCell]);
            }
        }
        WriteData(grploc, primitiveVariables[0], "res_unsteady_n2");

        RDouble **qAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("qAverage"));
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = static_cast<double>(qAverage[iEquation][iCell]);
            }
        }
        WriteData(grploc, primitiveVariables[0], "qAverage");

        double **tauVariables = NewPointer2<double>(6, nTotal);

        RDouble **tauAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("tauAverage"));
        for (int m = 0; m < 6; ++ m)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                tauVariables[m][iCell] = static_cast<double>(tauAverage[m][iCell]);
            }
        }
        WriteData(grploc, tauVariables[0], "tauAverage");

        RDouble **q2Average = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q2Average"));
        for (int m = 0; m < 6; ++ m)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                tauVariables[m][iCell] = static_cast<double>(q2Average[m][iCell]);
            }
        }
        WriteData(grploc, tauVariables[0], "q2Average");
        DelPointer2(tauVariables);
    }

    H5Gclose(grploc);

    DelPointer2(primitiveVariables);
}

void NSSolverUnstruct::ReadRestartH5(ActionKey *actkey)
{
    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    int isRestartChangeInflow = parameters->GetIsRestartChangeInflow();
    int nEquation = GetNumberOfEquations();

    int outputStepInterval = 0;
    ReadData(actkey->filepos, &outputStepInterval, "outnstep");
    GlobalDataBase::UpdateData("outnstep", &outputStepInterval, PHINT, 1);

    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));

    double **primitiveVariables = NewPointer2<double>(nEquation, nTotal);

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss<<"Group"<<gridID;
    grpName = oss.str();

    grploc = OpenGroup(actkey->filepos, grpName);

    int nTotalCellRestart = 0;
    ReadData(grploc, &nTotalCellRestart, "nTotalCell");

    if(nTotalCellRestart != nTotalCell)
    {
        ostringstream erroeInfo;
        erroeInfo << " Error: the cell number in flow.dat is not equal to the cell number in grid file !" << endl;
        TK_Exit::ExceptionExit(erroeInfo.str());
    }

    ReadData(grploc, primitiveVariables[0], "q");
    for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            q[iEquation][iCell] = static_cast<RDouble>(primitiveVariables[iEquation][iCell]);
        }
    }
    ReadData(grploc, t[ITT], "t");

    if (isRestartChangeInflow)
    {
        RDouble AoA = parameters->GetAoA();
        RDouble attack = AoA * PI / 180.0;
        RDouble angleSlide = parameters->GetAngleOfSlide();
        RDouble sideslip = angleSlide * PI / 180.0;
        RDouble coefficientofStateEquation = gas->GetCoefficientOfStateEquation();
        for (int iCell = 0; iCell < nTotal; ++iCell)
        {
            RDouble velocity = sqrt(q[IU][iCell] * q[IU][iCell] + q[IV][iCell] * q[IV][iCell] + q[IW][iCell] * q[IW][iCell]);
            q[IU][iCell] = velocity * cos(attack) * cos(sideslip);
            q[IV][iCell] = velocity * sin(attack) * cos(sideslip);
            q[IW][iCell] = velocity * sin(sideslip);
            q[IP][iCell] = coefficientofStateEquation * q[IR][iCell] * t[ITT][iCell];
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady)
    {
        H5Gclose(grploc);
        DelPointer2(primitiveVariables);

        return;
    }

    RDouble **qn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n1"));
    RDouble **qn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n2"));
    RDouble **res   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"           ));
    RDouble **resn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n1"));
    RDouble **resn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n2"));
    RDouble **qAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("qAverage"));
    RDouble **tauAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("tauAverage"));
    RDouble **q2Average = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q2Average"));

    int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();
    double physicalTimeIn = 0.0;

    if (!CheckDataExist(actkey->filepos, "physicalTime"))
    {
        ifStartFromSteadyResults = 1;
    }

    if (ifStartFromSteadyResults)
    {
        PrintToWindow("Restart from steady NS flow field, reset outer step to be zero!\n");
        //! Start from steady flow field.
        //! Reset the outer step when start from steady flow.
        outputStepInterval = 0;
        GlobalDataBase::UpdateData("outnstep", &outputStepInterval, PHINT, 1);

        for (int m = 0; m < nEquation; ++ m)
        {
            for (int iCell = 0; iCell < nTotalCell; ++ iCell)
            {
                qn1[m][iCell] = q[m][iCell];
                qn2[m][iCell] = q[m][iCell];
                res[m][iCell]   = 0.0;
                resn1[m][iCell] = 0.0;
                resn2[m][iCell] = 0.0;
            }
        }
    }
    else
    {
        ReadData(actkey->filepos, &physicalTimeIn, "physicalTime");

        ReadData(grploc, primitiveVariables[0], "q_unsteady_n1");
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                qn1[iEquation][iCell] = static_cast<RDouble>(primitiveVariables[iEquation][iCell]);
            }
        }

        ReadData(grploc, primitiveVariables[0], "q_unsteady_n2");
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                qn2[iEquation][iCell] = static_cast<RDouble>(primitiveVariables[iEquation][iCell]);
            }
        }

        ReadData(grploc, primitiveVariables[0], "res");
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                res[iEquation][iCell] = static_cast<RDouble>(primitiveVariables[iEquation][iCell]);
            }
        }

        ReadData(grploc, primitiveVariables[0], "res_unsteady_n1");
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                resn1[iEquation][iCell] = static_cast<RDouble>(primitiveVariables[iEquation][iCell]);
            }
        }

        ReadData(grploc, primitiveVariables[0], "res_unsteady_n2");
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                resn2[iEquation][iCell] = static_cast<RDouble>(primitiveVariables[iEquation][iCell]);
            }
        }
    }

    int nStatisticalStep = 0;

    if(IsNeedStatistics())
    {
        if(ifStartFromSteadyResults)
        {
            for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < nTotalCell; ++ iCell)
                {
                    qAverage[iEquation][iCell] = q[iEquation][iCell];
                }
            }
        }
        else
        {
            ReadData(grploc, qAverage[0], "qAverage");

            ReadData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
        }
    }

    if(IsNeedReynoldsStressStatistics())
    {
        if(ifStartFromSteadyResults)
        {
            for (int m = 0; m < 6; ++ m)
            {
                for (int iCell = 0; iCell < nTotalCell; ++ iCell)
                {
                    tauAverage[m][iCell] = 0.0;
                    q2Average[m][iCell] = 0.0;
                }
            }
        }
        else
        {
            ReadData(grploc, tauAverage[0], "tauAverage");

            ReadData(grploc, q2Average[0], "q2Average");
        }
    }

    GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);

    RDouble physicalTime = static_cast<RDouble>(physicalTimeIn);
    GlobalDataBase::UpdateData("physicalTime",&physicalTime, PHDOUBLE, 1);

    H5Gclose(grploc);
    DelPointer2(primitiveVariables);
}

void NSSolverUnstruct::InterpolateFromRestartH5(ActionKey *actkey)
{
    UnstructGrid *grid = UnstructGridCast(GetGrid(actkey->level));
    int gridID = grid->GetZoneID();

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nEquation = GetNumberOfEquations();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));

    int nTemperatureModel = parameters->GetTemperatureModel();
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();

    grid->CreateInterpolateInfo(actkey);

    int *targetZone = reinterpret_cast<int*>(grid->GetDataPtr("targetZone"));
    int *targetCell = reinterpret_cast<int*>(grid->GetDataPtr("targetCell"));

    for (int iEquation = 0; iEquation < nEquation; iEquation++)
    {
        for (int iCell = 0; iCell < nTotal; iCell++)
        {
            q[iEquation][iCell] = primitiveVarFarfield[iEquation];
        }
    }

    SetField(t, nTemperatureModel, nTotal, 1.0);

    int outputStepInterval = 0;
    GlobalDataBase::UpdateData("outnstep", &outputStepInterval, PHINT, 1);


    vector <int> intpolateZone;
    intpolateZone.resize(0);
    for (int iCell = 0; iCell < nTotalCell; iCell++)
    {
        intpolateZone.push_back(targetZone[iCell]);
    }
    sort(intpolateZone.begin(), intpolateZone.end());
    intpolateZone.erase(unique(intpolateZone.begin(), intpolateZone.end()), intpolateZone.end());

    for (int iZone = 0; iZone < intpolateZone.size(); iZone++)
    {
        hid_t gridData;
        gridData = OpenGroup(actkey->gridfilepos, "Information");

        int *zoneType = ReadIntOneRow(gridData, "block_type");

        hid_t gridDataloc;
        string grpName;

        ostringstream oss;
        oss<<"Grid-"<<gridID;
        grpName = oss.str();
        gridDataloc = OpenGroup(actkey->gridfilepos, grpName);
        int *iDimensions = ReadIntOneRow(gridDataloc, "iDimensions");
        if (zoneType[iZone] == UNSTRUCTGRID)
        {
            int nTotalCellGrid = iDimensions[2];

            int nBoundFaceGrid = 0;
            hid_t faceData;
            faceData = OpenGroup(gridDataloc, "FaceTopology");
            ReadData(faceData, &nBoundFaceGrid, "nBoundFace");
            int nTotalGrid = nTotalCellGrid + nBoundFaceGrid;
            H5Gclose(faceData);

            RDouble **qIn =  NewPointer2<double>(nEquation, nTotalGrid);
            RDouble **tIn =  NewPointer2<double>(nTemperatureModel, nTotalGrid);

            hid_t flowDataloc;
            string grpFlowName;

            ostringstream ossFlow;
            ossFlow << "Group" << intpolateZone[iZone];
            grpFlowName = ossFlow.str();

            flowDataloc = OpenGroup(actkey->filepos, grpFlowName);

            ReadData(flowDataloc, qIn[0], "q");
            ReadData(flowDataloc, tIn[0], "t");

            for (int iCell = 0; iCell < nTotalCellGrid; iCell++)
            {
                if (intpolateZone[iZone] != targetZone[iCell])
                {
                    continue;
                }

                for (int iEquation = 0; iEquation < nEquation; iEquation++)
                {
                    for (int iCell = 0; iCell < nTotal; iCell++)
                    {
                        q[iEquation][iCell] = qIn[iEquation][targetCell[iCell]];
                    }
                }
                for (int iTemperatureModel = 0; iTemperatureModel < nTemperatureModel; iTemperatureModel++)
                {
                    for (int iCell = 0; iCell < nTotal; iCell++)
                    {
                        t[iTemperatureModel][iCell] = tIn[iTemperatureModel][targetCell[iCell]];
                    }
                }
            }
            H5Gclose(flowDataloc);
            DelPointer2(qIn);
            DelPointer2(tIn);
        }
        H5Gclose(gridDataloc);
    }
}


void NSSolverUnstruct::ReadRestart(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nEquation = GetNumberOfEquations();

    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();
    int outnstep = 0;
    cdata->Read(&outnstep,sizeof(int));

    GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);

    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    for (int m = 0; m < nEquation; ++ m)
    {
        cdata->Read(q[m], sizeof(RDouble) * nTotal);
    }

    //ReadInterface(actkey);

    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady) return;

    RDouble **qn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n1"));
    RDouble **qn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n2"));

    RDouble **res   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"           ));
    RDouble **resn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n1"));
    RDouble **resn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n2"));

    int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();

    RDouble physicalTime = 0.0;
    if (ifStartFromSteadyResults)
    {
        PrintToWindow("Restart from steady NS flow field, reset outer step to be zero!\n");
        //! Start from steady flow field.
        //! Reset the outer step when start from steady flow.
        outnstep = 0;
        GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);

        for (int m = 0; m < nEquation; ++ m)
        {
            for (int iCell = 0; iCell < nTotalCell; ++ iCell)
            {
                qn1[m][iCell] = q[m][iCell];
                qn2[m][iCell] = q[m][iCell];

                res[m][iCell]   = 0.0;
                resn1[m][iCell] = 0.0;
                resn2[m][iCell] = 0.0;
            }
        }
    }
    else
    {
        //! Start from unsteady flow field.
        cdata->Read(&physicalTime, sizeof(RDouble));

        for (int m = 0; m < nEquation; ++ m)
        {
            cdata->Read(qn1[m], sizeof(RDouble) * nTotal);
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            cdata->Read(qn2[m], sizeof(RDouble) * nTotal);
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            cdata->Read(res[m],sizeof(RDouble) * nTotal);
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            cdata->Read(resn1[m], sizeof(RDouble) * nTotal);
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            cdata->Read(resn2[m], sizeof(RDouble) * nTotal);
        }
    }

    GlobalDataBase::UpdateData("physicalTime",&physicalTime, PHDOUBLE, 1);
}

void NSSolverUnstruct::InitFlowAsRestart()
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    UnstructGrid *grid = UnstructGridCast(GetGrid(0));
    int outnstep = 0;
    GlobalDataBase::UpdateData("outnstep",&outnstep, PHINT, 1);

    RDouble physicalTime = 0.0;
    GlobalDataBase::UpdateData("physicalTime",&physicalTime, PHDOUBLE, 1);

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    int nEquation = GetNumberOfEquations();

    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            q[m][iCell] = primitiveVarFarfield[m];
        }
    }

    int sysGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    int flowInitMethod = GlobalDataBase::GetIntParaFromDB("flowInitMethod");
    RDouble attack = PI / 180.0 * GlobalDataBase::GetDoubleParaFromDB("attackd");
    RDouble sideslip = PI / 180.0 * GlobalDataBase::GetDoubleParaFromDB("angleSlide");
    RDouble direction_inlet[3] = { 0.0,0.0,0.0 };
    direction_inlet[0] = cos(attack) * cos(sideslip);
    direction_inlet[1] = sin(attack) * cos(sideslip);
    direction_inlet[2] = sin(sideslip);
    if (flowInitMethod == 1 && sysGridType == UNSTRUCTGRID)
    {
        RDouble maxBoundaryLayerWalldist = GlobalDataBase::GetDoubleParaFromDB("maxBoundaryLayerWalldist");
        RDouble boundaryHeight = maxBoundaryLayerWalldist;
        if (boundaryHeight > PHSPACE::SMALL)
        {
            for (int iCell = 0; iCell < nTotal; ++iCell)
            {
                RDouble *wallDistance = grid->GetWallDist();
                RDouble *nearestwallfacenormalx = grid->GetNearestWallFaceNormalX();
                RDouble *nearestwallfacenormaly = grid->GetNearestWallFaceNormalY();
                RDouble *nearestwallfacenormalz = grid->GetNearestWallFaceNormalZ();
                if (wallDistance[iCell] < boundaryHeight)
                {
                    //! iCell lay inside the boundary layer.
                    RDouble enta = MIN(wallDistance[iCell] / boundaryHeight, 1.0);
                    //! Normal velocity.
                    RDouble dot = direction_inlet[0] * nearestwallfacenormalx[iCell] + direction_inlet[1] * nearestwallfacenormaly[iCell] + direction_inlet[2] * nearestwallfacenormalz[iCell];
                    RDouble normalVelocity[3];
                    normalVelocity[0] = dot * nearestwallfacenormalx[iCell];
                    normalVelocity[1] = dot * nearestwallfacenormaly[iCell];
                    normalVelocity[2] = dot * nearestwallfacenormalz[iCell];
                    //! Tangential velocity.
                    RDouble tangentialVelocity[3];
                    tangentialVelocity[0] = direction_inlet[0] - normalVelocity[0];
                    tangentialVelocity[1] = direction_inlet[1] - normalVelocity[1];
                    tangentialVelocity[2] = direction_inlet[2] - normalVelocity[2];
                    //! Tangential velocity normalization.
                    RDouble tangentialVelocitymodulus = sqrt(tangentialVelocity[0] * tangentialVelocity[0] + tangentialVelocity[1] * tangentialVelocity[1] + tangentialVelocity[2] * tangentialVelocity[2]);
                    if (tangentialVelocitymodulus > PHSPACE::SMALL)
                    {
                        RDouble tangentialVelocitynormalization[3];
                        tangentialVelocitynormalization[0] = tangentialVelocity[0] / tangentialVelocitymodulus;
                        tangentialVelocitynormalization[1] = tangentialVelocity[1] / tangentialVelocitymodulus;
                        tangentialVelocitynormalization[2] = tangentialVelocity[2] / tangentialVelocitymodulus;
                        q[IU][iCell] = tangentialVelocitynormalization[0] * (2 * enta - 2 * enta * enta * enta + enta * enta * enta * enta);
                        q[IV][iCell] = tangentialVelocitynormalization[1] * (2 * enta - 2 * enta * enta * enta + enta * enta * enta * enta);
                        q[IW][iCell] = tangentialVelocitynormalization[2] * (2 * enta - 2 * enta * enta * enta + enta * enta * enta * enta);
                    }
                    else 
                    {
                        q[IU][iCell] = direction_inlet[0] * (2 * enta - 2 * enta * enta * enta + enta * enta * enta * enta);
                        q[IV][iCell] = direction_inlet[1] * (2 * enta - 2 * enta * enta * enta + enta * enta * enta * enta);
                        q[IW][iCell] = direction_inlet[2] * (2 * enta - 2 * enta * enta * enta + enta * enta * enta * enta);
                    }
                }
            }
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();

    if (!isUnsteady)
    {
        RDouble **res = reinterpret_cast <RDouble **> (grid->GetDataPtr("res"));
        for (int m = 0; m < nEquation; ++ m)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                res[m][iCell] = 0.0;
            }
        }
        return;
    }

    RDouble **qn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n1"));
    RDouble **qn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n2"));

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            qn1[m][iCell] = q[m][iCell];
        }
    }

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            qn2[m][iCell] = q[m][iCell];
        }
    }

    RDouble **res   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"           ));
    RDouble **resn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n1"));
    RDouble **resn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n2"));

    //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
    RDouble **resTmp = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_tmp"));

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            res  [m][iCell] = 0.0;
            resn1[m][iCell] = 0.0;
            resn2[m][iCell] = 0.0;
            resTmp[m][iCell] = 0.0;
        }
    }

    int nStatisticalStep = 0;
    GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);

}

void NSSolverUnstruct::UploadInterfaceData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo)
    {
        return;
    }

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    PHSPACE::UploadInterfaceValue(grid, q, "q",  nEquation);

    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));
    PHSPACE::UploadInterfaceValue(grid, t, "t", nTemperatureModel);
 }

void NSSolverUnstruct::DownloadInterfaceData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast< UnstructGrid * >(GetGrid(level));
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **q = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q"  ) );
    PHSPACE::DownloadInterfaceValue(grid, q, "q", nEquation);

    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t" ));
    PHSPACE::DownloadInterfaceValue(grid, t, "t", nTemperatureModel);
}

void NSSolverUnstruct::UploadInterpointData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast< UnstructGrid * >(GetGrid(level));
    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation == 0)
    {
        return;
    }
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **qnode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qnode"));
    PHSPACE::UploadInterpointValue(grid, qnode, "qnode",  nEquation);

    RDouble **tnode = reinterpret_cast<RDouble **>(grid->GetDataPtr("tnode"));
    PHSPACE::UploadInterpointValue(grid, tnode, "tnode", nTemperatureModel);
}

void NSSolverUnstruct::DownloadInterpointData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));
    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation == 0)
    {
        return;
    }
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **qnode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qnode"));
    PHSPACE::DownloadInterpointValue(grid, qnode, "qnode", nEquation);

    RDouble **tnode = reinterpret_cast<RDouble **> (grid->GetDataPtr("tnode"));
    PHSPACE::DownloadInterpointValue(grid, tnode, "tnode", nTemperatureModel);
}

void NSSolverUnstruct::CommunicationInterpointWeight(ActionKey *actkey)
{
    int level = actkey->level;

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int myid = GetCurrentProcessorID();
        int tag = iZone + level * nZones;
        int send_proc = GetZoneProcessorID(iZone);

        ZoneNeighbor *zoneNeighborForPoint = zoneConnectivityForPoint->GetZoneNeighborForPoint(iZone);
        std::size_t numberOfNeighbor = zoneNeighborForPoint->GetNumberOfNeighbor();

        for (std::size_t iNeighborZone = 0; iNeighborZone < numberOfNeighbor; ++ iNeighborZone)
        {
            int neighborZoneIndex = zoneNeighborForPoint->GetZoneIndexOfNeighbor(iNeighborZone);
            int recv_proc = GetZoneProcessorID(neighborZoneIndex);

            DataContainer *cdata = new DataContainer();

            if (myid == send_proc)
            {
                UnstructGrid *grid = UnstructGridCast(PHSPACE::GetGrid(iZone, level));
                InterpointInformation *interpointInformation = grid->GetInterpointInfo();
                RDouble *nodeWeight = reinterpret_cast<RDouble *>(grid->GetDataPtr("nodeWeight"));

                cdata->MoveToBegin();

                int iNeighbor = interpointInformation->FindIthNeighbor(neighborZoneIndex);
                int numberOfInterpointsForNeighbor = interpointInformation->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int *pointIndexForSend = interpointInformation->GetPointIndexForSend(iNeighbor);
                int *interPoint2GlobalPoint = interpointInformation->GetInterPoint2GlobalPoint();

                int iPoint;
                int globalPoint;
                RDouble valueWeight;
                RDouble *qtmp = new RDouble[numberOfInterpointsForNeighbor];
                for (int ipointLocal = 0; ipointLocal < numberOfInterpointsForNeighbor; ++ ipointLocal)
                {
                    iPoint = pointIndexForSend[ipointLocal];
                    globalPoint = interPoint2GlobalPoint[iPoint];
                    valueWeight = nodeWeight[globalPoint];

                    qtmp[ipointLocal] = valueWeight;
                }
                cdata->Write(qtmp, numberOfInterpointsForNeighbor * sizeof(RDouble));
                delete [] qtmp;     qtmp = nullptr;
            }

            PH_Trade(cdata, send_proc, recv_proc, tag);

            if (myid == recv_proc)
            {
                UnstructGrid *gridNeighbor = UnstructGridCast(PHSPACE::GetGrid(neighborZoneIndex, level));
                InterpointInformation *interpointInformationNeighbor = gridNeighbor->GetInterpointInfo();
                RDouble *nodeWeightTrade = reinterpret_cast<RDouble *>(gridNeighbor->GetDataPtr("nodeWeightTrade"));

                cdata->MoveToBegin();

                int iNeighbor = interpointInformationNeighbor->FindIthNeighbor(iZone);
                int numberOfInterpointsForNeighbor = interpointInformationNeighbor->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int *pointIndexForReceive = interpointInformationNeighbor->GetPointIndexForReceive(iNeighbor);
                int *interPoint2GlobalPoint = interpointInformationNeighbor->GetInterPoint2GlobalPoint();

                int iPoint;
                int globalPoint;
                RDouble valueWeight;
                for (int j = 0; j < numberOfInterpointsForNeighbor; ++ j)
                {
                    iPoint = pointIndexForReceive[j];
                    globalPoint = interPoint2GlobalPoint[iPoint];
                    cdata->Read(&valueWeight, sizeof(RDouble));

                    nodeWeightTrade[globalPoint] += valueWeight;
                }
            }
            delete cdata;    cdata = nullptr;
        }
    }
}

void NSSolverUnstruct::DownloadInterpointWeight(ActionKey *actkey)
{
    int level = actkey->level;

    using namespace PHMPI;

    int myid = GetCurrentProcessorID();
    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int currentZoneProc = GetZoneProcessorID(iZone);

        if (myid != currentZoneProc)
        {
            continue;
        }

        UnstructGrid *grid = UnstructGridCast(PHSPACE::GetGrid(iZone, level));
        int nTotalNode = grid->GetNTotalNode();
        RDouble *nodeWeight = reinterpret_cast<RDouble *>(grid->GetDataPtr("nodeWeight"));
        RDouble *nodeWeightTrade = reinterpret_cast<RDouble *>(grid->GetDataPtr("nodeWeightTrade"));

        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            nodeWeight[iNode] += nodeWeightTrade[iNode];
        }

        PHSPACE::SetField(nodeWeightTrade, 0.0, nTotalNode);
    }
}

void NSSolverUnstruct::UploadOversetData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast< UnstructGrid * >(GetGrid(level));

    int zoneIndex = grid->GetZoneID();
    OversetTopologyManager *oversetTopologyManager = PHSPACE::GetOversetTopologyManager();
    uint_t numberOfOversetNeighbors = oversetTopologyManager->GetNumberOfOversetNeighbors(zoneIndex);
    if (0 == numberOfOversetNeighbors)
    {
        return;
    }
   
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();
    int isSolve = GlobalDataBase::GetIntParaFromDB("isSolve");
    //int isNeedCommDqforOverset = GlobalDataBase::GetIntParaFromDB("isNeedCommDqforOverset");
    if (!isSolve)
    {
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q" ));
    PHSPACE::UploadOversetValue(grid, q, "q",  nEquation);

    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t" ));
    PHSPACE::UploadOversetValue(grid, t, "t",  nTemperatureModel);
}
    else
    {
        RDouble **gradPrimtiveVarX = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradPrimtiveVarX"));
        RDouble **gradPrimtiveVarY = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradPrimtiveVarY"));
        RDouble **gradPrimtiveVarZ = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradPrimtiveVarZ"));
        PHSPACE::UploadOversetValue(grid, gradPrimtiveVarX, "dqdx", nEquation);
        PHSPACE::UploadOversetValue(grid, gradPrimtiveVarY, "dqdy", nEquation);
        PHSPACE::UploadOversetValue(grid, gradPrimtiveVarZ, "dqdz", nEquation);

        RDouble **limit2D = reinterpret_cast <RDouble**> (grid->GetDataPtr("limit2D"));
        PHSPACE::UploadOversetValue(grid, limit2D, "limit2D", nEquation);

        RDouble **gradTemperatureX = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradTemperatureX"));
        RDouble **gradTemperatureY = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradTemperatureY"));
        RDouble **gradTemperatureZ = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradTemperatureZ"));

        PHSPACE::UploadOversetValue(grid, gradTemperatureX, "dtdx", nTemperatureModel);
        PHSPACE::UploadOversetValue(grid, gradTemperatureY, "dtdy", nTemperatureModel);
        PHSPACE::UploadOversetValue(grid, gradTemperatureZ, "dtdz", nTemperatureModel);
    }
}

void NSSolverUnstruct::DownloadOversetData(ActionKey *actkey)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast< UnstructGrid * >(GetGrid(level));

    int zoneIndex = grid->GetZoneID();
    OversetTopologyManager *oversetTopologyManager = PHSPACE::GetOversetTopologyManager();
    uint_t numberOfOversetNeighbors = oversetTopologyManager->GetNumberOfOversetNeighbors(zoneIndex);
    if (0 == numberOfOversetNeighbors)
    {
        return;
    }

    int nTemperatureModel = parameters->GetTemperatureModel();
    int nEquation = GetNumberOfEquations();

    int isSolve = GlobalDataBase::GetIntParaFromDB("isSolve");
    if (!isSolve)
    {
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q" ));
    PHSPACE::DownloadOversetValue(grid, q, "q", nEquation);

    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t" ));
    PHSPACE::DownloadOversetValue(grid, t, "t",  nTemperatureModel);
}
    else
    {
        RDouble **gradPrimtiveVarX = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradPrimtiveVarX"));
        RDouble **gradPrimtiveVarY = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradPrimtiveVarY"));
        RDouble **gradPrimtiveVarZ = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradPrimtiveVarZ"));
        PHSPACE::DownloadOversetValue(grid, gradPrimtiveVarX, "dqdx", nEquation);
        PHSPACE::DownloadOversetValue(grid, gradPrimtiveVarY, "dqdy", nEquation);
        PHSPACE::DownloadOversetValue(grid, gradPrimtiveVarZ, "dqdz", nEquation);
        
        RDouble **limit2D = reinterpret_cast <RDouble**> (grid->GetDataPtr("limit2D"));
        PHSPACE::DownloadOversetValue(grid, limit2D, "limit2D", nEquation);

        RDouble **gradTemperatureX = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradTemperatureX"));
        RDouble **gradTemperatureY = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradTemperatureY"));
        RDouble **gradTemperatureZ = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradTemperatureZ"));

        PHSPACE::DownloadOversetValue(grid, gradTemperatureX, "dtdx", nTemperatureModel);
        PHSPACE::DownloadOversetValue(grid, gradTemperatureY, "dtdy", nTemperatureModel);
        PHSPACE::DownloadOversetValue(grid, gradTemperatureZ, "dtdz", nTemperatureModel);
    }
}

void NSSolverUnstruct::InitMixingPlane(RDouble ***MixingPlaneVar, int Dim1, int Dim2, int Dim3, RDouble Value)
{

    for (int iDim1 = 0; iDim1 < Dim1; iDim1++)
    {
        for (int iDim2 = 0; iDim2 < Dim2; iDim2++)
        {
            for (int iDim3 = 0; iDim3 < Dim3; iDim3++)
            {
                MixingPlaneVar[iDim1][iDim2][iDim3] = Value;
            }
        }
    }
}

//! need to input q_average, which means flux averaged primitive variable of mixing plane.
//! this function should be called in boundary condition part.
void NSSolverUnstruct::AverageMixingPlane(Grid *grid)
{
    UnstructGrid *gridUnstruct = UnstructGridCast(grid);
    UnstructBCSet *unstructBCSet = gridUnstruct->GetUnstructBCSet();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refGama = parameters->GetRefGama();
    RDouble gama1 = refGama - 1.0;

    RDouble refDensity = parameters->GetRefDimensionalDensity();
    RDouble refPressure = parameters->GetRefDimensionalPressure();

    RDouble **q = reinterpret_cast< RDouble ** > (gridUnstruct->GetDataPtr("q"));

    RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));

    RDouble **r_ori    = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Face"));
    RDouble **r_target = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Span"));

    RDouble **SpanArea  = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("SpanArea"));

    RDouble **q_Face    = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("q_Face"));
    RDouble **qCon_Face = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("qCon_Face"));

    //! nSpan is the number of spanwise point to do data average
    //! also interpolate data on these points.
    int nSpanSection = PHSPACE::GlobalDataBase::GetIntParaFromDB("nSpanSection");
    RDouble nSpan = (RDouble)(nSpanSection);

    RDouble Area;

    RDouble *nxs, *nys, *nzs, *ns;

    nxs = gridUnstruct->GetFaceNormalX();
    nys = gridUnstruct->GetFaceNormalY();
    nzs = gridUnstruct->GetFaceNormalZ();
    ns = gridUnstruct->GetFaceArea();

    int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");

    string MixingPlane[100];
    GlobalDataBase::GetData("MixingPlane", &MixingPlane, PHSTRING, 2 * nTurboZone - 2);

    //! should define r_ori means different radius of each cell.
    //! also need to define r_target means the span radius that average data.
    //! r_min and r_max is used to find the maximum and minimum radius on mixing plane.
    //! use r_min and r_max to compute the range of interpolation.
    //! this also need compute during AverageMixingPlane.
    RDouble **nxsSpan = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("nxsSpan"));
    RDouble **ntsSpan = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("ntsSpan"));
    RDouble **nrsSpan = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("nrsSpan"));

    //! r_ori should be an array that contains coordinate of each node.
    //! the range of this array is equal to the face number of this boundary.
    RDouble r_min, r_max;

    r_min = 1e6;
    r_max = 0;

    RDouble *xfc = gridUnstruct->GetFaceCenterX();
    RDouble *yfc = gridUnstruct->GetFaceCenterY();
    RDouble *zfc = gridUnstruct->GetFaceCenterZ();

    RDouble *xcc = gridUnstruct->GetCellCenterX();
    RDouble *ycc = gridUnstruct->GetCellCenterY();
    RDouble *zcc = gridUnstruct->GetCellCenterZ();

    RDouble *vgn = gridUnstruct->GetFaceNormalVelocity();

    int **face2nodeArray = gridUnstruct->GetFace2NodeArray();
    int * node_number_of_each_face = gridUnstruct->GetNodeNumberOfEachFace();

    int *leftCellofFace = gridUnstruct->GetLeftCellOfFace();
    int *rightCellofFace = gridUnstruct->GetRightCellOfFace();
    int nTotalCell = gridUnstruct->GetNTotalCell();

    int nBCRegion = unstructBCSet->GetnBCRegion();

    RDouble PeriodicRotationAngle[100];
    GlobalDataBase::GetData("PeriodicRotationAngle", &PeriodicRotationAngle, PHDOUBLE, nTurboZone);

    using namespace IDX;

    RDouble ***qCon_Span  = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("qCon_Span"));
    RDouble ***SpanFlux   = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("SpanFlux"));

    int FlowDirection = -1;

    int nEquation = GetNumberOfEquations();

    RDouble *primi = new RDouble[nEquation]();
    RDouble *primo = new RDouble[nEquation]();

    //! this will remove later, iBCRegion should be input data.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        vector<int> *faceIndex = bcRegion->GetFaceIndex();

        int bcType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();

        //! find mixing plane and average data.
        //! use new bcType to define mixing plane.
        //! if (bcType == -1)

        //! Parallel
        int iTurboZone = gridUnstruct->GetOrdinaryGridIndex();
        //! Serial
        if (iTurboZone == -1)
        {
            iTurboZone = gridUnstruct->GetZoneID();
        }

        RDouble rotationAngle = PeriodicRotationAngle[iTurboZone] / 180.0 * PI;

        string MixingPlaneIn, MixingPlaneOut;
        int SourceFlag;

        if (iTurboZone == 0)
        {
            MixingPlaneOut = MixingPlane[0];
        }
        else if (iTurboZone == nTurboZone)
        {
            MixingPlaneIn = MixingPlane[2 * iTurboZone - 1];
        }
        else
        {
            MixingPlaneIn  = MixingPlane[2 * iTurboZone - 1];
            MixingPlaneOut = MixingPlane[2 * iTurboZone ];
        }

        if (bcName == MixingPlaneIn || bcName == MixingPlaneOut)
        {
            if (bcName == MixingPlaneIn)
            {
                FlowDirection = 0;
                SourceFlag    = 1 + (iTurboZone - 1) * 10;
            }
            if (bcName == MixingPlaneOut)
            {
                FlowDirection = 1;
                SourceFlag    = 0 + (iTurboZone + 1) * 10;
            }
            int MixingFlag = FlowDirection + iTurboZone * 10;

            //! preprocess of mixing plane.
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;

                int le = leftCellofFace[iFace];
                int re = rightCellofFace[iFace];

                //! need to compute total area, total flux of this bc region.
                //! also need to compute flux of each spanwise.
                //! spanwise can be defined based on mesh or manually.
                r_ori[MixingFlag][iFace] = sqrt(yfc[iFace] * yfc[iFace] + zfc[iFace] * zfc[iFace]);

                r_min = MIN(r_ori[MixingFlag][iFace], r_min);
                r_max = MAX(r_ori[MixingFlag][iFace], r_max);
            }

            //! get target radius.
            for (int iSpan = 0; iSpan <= nSpanSection; iSpan++)
            {
                r_target[MixingFlag][iSpan] = r_min + iSpan * (r_max - r_min) / nSpan;
            }

            //! compute average flux.
            for (int iSpan = 0; iSpan < nSpanSection; iSpan++)
            {
                RDouble Flux[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
                RDouble TotalArea = 0.0;

                nxsSpan[MixingFlag][iSpan] = 0.0;
                ntsSpan[MixingFlag][iSpan] = 0.0;
                nrsSpan[MixingFlag][iSpan] = 0.0;

                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    //! iFace is the face number in the set of faceIndex.
                    int iFace = *iter;

                    int le = leftCellofFace[iFace];
                    int re = rightCellofFace[iFace];

                    Area = ns[iFace];

                    for (int m = 0; m < 5; m++)
                    {
                        primi[m] = q[m][le];
                        primo[m] = q[m][re];
                    }

                    RDouble &rmi = primi[IR];
                    RDouble &umi = primi[IU];
                    RDouble &vmi = primi[IV];
                    RDouble &wmi = primi[IW];
                    RDouble &pmi = primi[IP];

                    RDouble &rmo = primo[IR];
                    RDouble &umo = primo[IU];
                    RDouble &vmo = primo[IV];
                    RDouble &wmo = primo[IW];
                    RDouble &pmo = primo[IP];

                    RDouble &vb = vgn[iFace];
                    //! need to compute velocity in radial and tangential direction.
                    //! 1. compute face normal in radial and tangential direction.
                    RDouble r_Face  = sqrt(yfc[iFace] * yfc[iFace] + zfc[iFace] * zfc[iFace]);
                    RDouble nts     = (nzs[iFace] * yfc[iFace] - nys[iFace] * zfc[iFace]) / r_Face;
                    RDouble nrs     = (nys[iFace] * yfc[iFace] + nzs[iFace] * zfc[iFace]) / r_Face;

                    //! 2. compute velocity in radial and tangential direction.
                    RDouble Velocity_ti = (vmi * zfc[iFace] - wmi * yfc[iFace]) / r_Face;
                    RDouble Velocity_ri = (vmi * yfc[iFace] + wmi * zfc[iFace]) / r_Face;

                    RDouble Velocity_to = (vmo * zfc[iFace] - wmo * yfc[iFace]) / r_Face;
                    RDouble Velocity_ro = (vmo * yfc[iFace] + wmo * zfc[iFace]) / r_Face;

                    RDouble vni = umi * nxs[iFace] + Velocity_ri * nrs - vb;
                    RDouble vno = umo * nxs[iFace] + Velocity_ro * nrs - vb;
                    //! find face center that locate in target span.
                    //! compute flux at each span.
                    if (r_Face > r_target[MixingFlag][iSpan] && r_Face <= r_target[MixingFlag][iSpan + 1])
                    {
                        //! direction
                        nxsSpan[MixingFlag][iSpan] += nxs[iFace] * Area;
                        ntsSpan[MixingFlag][iSpan] += nts * Area;
                        nrsSpan[MixingFlag][iSpan] += nrs * Area;
                        TotalArea += Area;

                        RDouble Hmi = (refGama / (refGama - 1)) * (pmi / rmi) + half * (umi * umi + vmi * vmi + wmi * wmi);
                        RDouble Hmo = (refGama / (refGama - 1)) * (pmo / rmo) + half * (umo * umo + vmo * vmo + wmo * wmo);

                        Flux[IR ] += (rmi * vni) * ns[iFace];
                        Flux[IRU] += (rmi * vni * umi + pmi * nxs[iFace]) * ns[iFace];
                        Flux[IRV] += (rmi * vni * Velocity_ti) * ns[iFace];
                        Flux[IRW] += (rmi * vni * Velocity_ri + pmi * nrs) * ns[iFace];
                        Flux[IRE] += (rmi * vni * Hmi + pmi * vgn[iFace])* ns[iFace];
                    }
                }
                nxsSpan[MixingFlag][iSpan] /= TotalArea;
                ntsSpan[MixingFlag][iSpan] /= TotalArea;
                nrsSpan[MixingFlag][iSpan] /= TotalArea;

                SpanArea[MixingFlag][iSpan] = TotalArea;

                for (int m = 0; m < 5; m++)
                {
                    qCon_Span[MixingFlag][m][iSpan] = Flux[m];
                }
            }

            for (int iSpan = 0; iSpan < nSpanSection; iSpan++)
            {
                //! SpanFlux will communicate between upstream and downstream of mixing plane.
                for (int m = 0; m < 5; m++)
                {
                    SpanFlux[MixingFlag][m][iSpan] = qCon_Span[MixingFlag][m][iSpan] / SpanArea[MixingFlag][iSpan];
                }
            }
        }
    }
}

//! SpanFlux from AverageMixingPlane should be input data.
//! Here we will obtain averaged primitive variable of each span from SpanFlux.
//! The averaged primitive variable will be treated as outlet data of upstream and inlet data of downsteam.
void NSSolverUnstruct::MixingPlaneDataTransfer(Grid *grid, Grid *NeighborGrid)
{
    UnstructGrid *gridUnstruct = UnstructGridCast(grid);
    UnstructGrid *gridSource   = UnstructGridCast(NeighborGrid);
    UnstructBCSet *unstructBCSet = gridUnstruct->GetUnstructBCSet();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refGama = parameters->GetRefGama();
    RDouble gama1 = refGama - 1;

    RDouble refDensity = parameters->GetRefDimensionalDensity();
    RDouble refPressure = parameters->GetRefDimensionalPressure();

    RDouble **q = reinterpret_cast< RDouble ** > (gridUnstruct->GetDataPtr("q"));

    RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));

    RDouble **r_target  = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Span"));

    RDouble **SpanArea = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("SpanArea"));

    //! nSpan is the number of spanwise point to do data average
    //! also interpolate data on these points.
    //! nSpanSection
    int nSpanSection = PHSPACE::GlobalDataBase::GetIntParaFromDB("nSpanSection");

    int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
    string MixingPlane[100];
    GlobalDataBase::GetData("MixingPlane", &MixingPlane, PHSTRING, 2 * nTurboZone - 2);

    RDouble *xfc = gridUnstruct->GetFaceCenterX();
    RDouble *yfc = gridUnstruct->GetFaceCenterY();
    RDouble *zfc = gridUnstruct->GetFaceCenterZ();

    int nBCRegion = unstructBCSet->GetnBCRegion();

    RDouble **nxsSpan = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("nxsSpan"));
    RDouble **ntsSpan = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("ntsSpan"));
    RDouble **nrsSpan = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("nrsSpan"));

    //! qSpanNeighbor: record flow variable in neighbor zone at each span.
    RDouble ***q_Span        = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("q_Span"));
    RDouble ***SpanFlux      = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("SpanFlux"));
    RDouble ***qSpanNeighbor = reinterpret_cast< RDouble *** > (gridSource->GetDataPtr("qSpanNeighbor"));

    using namespace IDX;

    int FlowDirection = -1;

    int nEquation = GetNumberOfEquations();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        vector<int> *faceIndex = bcRegion->GetFaceIndex();

        int bcType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();

        //! Parallel
        int iTurboZone = gridUnstruct->GetOrdinaryGridIndex();
        int SourceZone = gridSource->GetOrdinaryGridIndex();
        //! Serial
        if (iTurboZone == -1)
        {
            iTurboZone = gridUnstruct->GetZoneID();
            SourceZone = gridSource->GetZoneID();
        }

        //! use new bcType to define mixing plane.
        //! if (bcType == -1)
        string MixingPlaneIn, MixingPlaneOut;

        if (iTurboZone == 0)
        {
            MixingPlaneOut = MixingPlane[0];
        }
        else if (iTurboZone == nTurboZone)
        {
            MixingPlaneIn = MixingPlane[2 * iTurboZone - 1];
        }
        else
        {
            MixingPlaneIn = MixingPlane[2 * iTurboZone - 1];
            MixingPlaneOut = MixingPlane[2 * iTurboZone];
        }

        //! TargetMixingFlag is Mixingplane in current zone, which need neighbor zone data.
        //! SourceMixingFlag is Mixingplane in neighbor zone, which provide the data. 
        int MixingFlag, TargetFlag;
        if (bcName == MixingPlaneIn || bcName == MixingPlaneOut)
        {
            //! downstream to upstream: transfer data from MixingIn to MixingOut.
            if (SourceZone == iTurboZone + 1)
            {
                if (bcName == MixingPlaneIn) continue;
                //! upstream MixingIn.
                TargetFlag = 0 + SourceZone * 10;
                //! current MixingOut.
                MixingFlag = 1 + iTurboZone * 10;
            }
            //! upstream to downstream: transfer data from MixingOut to MixingIn.
            else if (SourceZone == iTurboZone - 1)
            {
                if (bcName == MixingPlaneOut) continue;
                TargetFlag = 1 + SourceZone * 10;
                MixingFlag = 0 + iTurboZone * 10;
            }

            //! compute flow variables.
            for (int iSpan = 0; iSpan < nSpanSection ; iSpan++)
            {
                RDouble a, b, c, discr;

                a = nxsSpan[MixingFlag][iSpan] * SpanFlux[MixingFlag][IRU][iSpan]
                    + nrsSpan[MixingFlag][iSpan] * SpanFlux[MixingFlag][IRW][iSpan];

                b = -nrsSpan[MixingFlag][iSpan] * SpanFlux[MixingFlag][IRU][iSpan]
                    +  nxsSpan[MixingFlag][iSpan] * SpanFlux[MixingFlag][IRW][iSpan];

                c = a * a+ b * b + SpanFlux[MixingFlag][IRV][iSpan] * SpanFlux[MixingFlag][IRV][iSpan];

                discr = c - 2 * SpanFlux[MixingFlag][IR][iSpan] * SpanFlux[MixingFlag][IP][iSpan];

                q_Span[MixingFlag][IP][iSpan] = (a + sqrt(a * a + (refGama*refGama - 1) * discr)) / (refGama + 1);

                RDouble vn = (a - q_Span[MixingFlag][IP][iSpan]) / SpanFlux[MixingFlag][IR][iSpan];
                RDouble vt = b / SpanFlux[MixingFlag][IR][iSpan];

                //! IU, IV, IW refers to vn, vtheta, vt
                q_Span[MixingFlag][IU][iSpan] = vn;
                q_Span[MixingFlag][IV][iSpan] = SpanFlux[MixingFlag][IRV][iSpan] / SpanFlux[MixingFlag][IR][iSpan];
                q_Span[MixingFlag][IW][iSpan] = vt;
                q_Span[MixingFlag][IR][iSpan] = SpanFlux[MixingFlag][IR][iSpan] / vn;

                for (int m = 0; m < nEquation; m++)
                {
                    qSpanNeighbor[TargetFlag][m][iSpan] = q_Span[MixingFlag][m][iSpan];
                }
                //! velocity direction is opposite, because flow in and flow out.
                qSpanNeighbor[TargetFlag][IU][iSpan] = -q_Span[MixingFlag][IU][iSpan];
                qSpanNeighbor[TargetFlag][IW][iSpan] = -q_Span[MixingFlag][IW][iSpan];
            }
        }
    }
}

void NSSolverUnstruct::NonReflective(Grid *grid, Grid *NeighborGrid)
{
    UnstructGrid *gridUnstruct = UnstructGridCast(grid);
    UnstructGrid *gridTarget = UnstructGridCast(NeighborGrid);
    UnstructBCSet *unstructBCSet = gridUnstruct->GetUnstructBCSet();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refPressure = parameters->GetRefDimensionalPressure();
    RDouble refGama = parameters->GetRefGama();
    RDouble gama1 = refGama - 1;

    RDouble **q = reinterpret_cast< RDouble ** > (gridUnstruct->GetDataPtr("q"));

    RDouble **r_ori    = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Face"));
    RDouble **r_target = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Span"));

    RDouble **SpanArea = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("SpanArea"));
    //! nSpan is the number of spanwise point to do data average
    //! also interpolate data on these points.
    int nSpanSection = PHSPACE::GlobalDataBase::GetIntParaFromDB("nSpanSection");

    int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");

    string MixingPlane[100];
    GlobalDataBase::GetData("MixingPlane", &MixingPlane, PHSTRING, 2 * nTurboZone - 2);

    int *leftCellofFace  = gridUnstruct->GetLeftCellOfFace();
    int *rightCellofFace = gridUnstruct->GetRightCellOfFace();
    int nTotalCell = gridUnstruct->GetNTotalCell();

    int nBCRegion = unstructBCSet->GetnBCRegion();

    RDouble *xcc = gridUnstruct->GetCellCenterX();
    RDouble *ycc = gridUnstruct->GetCellCenterY();
    RDouble *zcc = gridUnstruct->GetCellCenterZ();

    RDouble *xfc = gridUnstruct->GetFaceCenterX();
    RDouble *yfc = gridUnstruct->GetFaceCenterY();
    RDouble *zfc = gridUnstruct->GetFaceCenterZ();

    RDouble *nxs, *nys, *nzs, *ns;

    nxs = gridUnstruct->GetFaceNormalX();
    nys = gridUnstruct->GetFaceNormalY();
    nzs = gridUnstruct->GetFaceNormalZ();
    ns  = gridUnstruct->GetFaceArea();

    RDouble **nxsSpan = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("nxsSpan"));
    RDouble **ntsSpan = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("ntsSpan"));
    RDouble **nrsSpan = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("nrsSpan"));

    RDouble ***q_SpanCurrent = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("q_Span"));
    RDouble ***q_SpanTarget  = reinterpret_cast< RDouble *** > (gridTarget->GetDataPtr("q_Span"));
    RDouble ***qSpanNeighbor = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("qSpanNeighbor"));

    RDouble ***dqSpanIn = reinterpret_cast <RDouble ***> (gridUnstruct->GetDataPtr("dqSpanIn"));
    RDouble ***dqSpanEx = reinterpret_cast <RDouble ***> (gridUnstruct->GetDataPtr("dqSpanEx"));
    RDouble ***dcSpan   = reinterpret_cast <RDouble ***> (gridUnstruct->GetDataPtr("dcSpan"));

    //! angular velocity
    RDouble Omega[100];
    GlobalDataBase::GetData("Omega", &Omega, PHDOUBLE, nTurboZone);

    using namespace IDX;

    int FlowDirection = -1;

    int nEquation = GetNumberOfEquations();
    RDouble *ExtraDq = new RDouble[nEquation];
    for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        vector<int> *faceIndex = bcRegion->GetFaceIndex();

        int bcType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();

        //! Parallel
        int iTurboZone = gridUnstruct->GetOrdinaryGridIndex();
        int TargetZone = gridTarget->GetOrdinaryGridIndex();
        //! Serial
        if (iTurboZone == -1)
        {
            iTurboZone = gridUnstruct->GetZoneID();
            TargetZone = gridTarget->GetZoneID();
        }

        string MixingPlaneIn, MixingPlaneOut;

        if (iTurboZone == 0)
        {
            MixingPlaneOut = MixingPlane[0];
        }
        else if (iTurboZone == nTurboZone)
        {
            MixingPlaneIn = MixingPlane[2 * iTurboZone - 1];
        }
        else
        {
            MixingPlaneIn  = MixingPlane[2 * iTurboZone - 1];
            MixingPlaneOut = MixingPlane[2 * iTurboZone];
        }

        //! TargetMixingFlag is Mixingplane in current zone, which need neighbor zone data.
        //! SourceMixingFlag is Mixingplane in neighbor zone, which provide the data. 
        int MixingFlag, TargetFlag;
        if (bcName == MixingPlaneIn || bcName == MixingPlaneOut)
        {
            //! downstream to upstream: transfer data from MixingIn to MixingOut.
            if (TargetZone == iTurboZone + 1)
            {
                if (bcName == MixingPlaneIn) continue;
                //! upstream MixingIn.
                TargetFlag = 0 + TargetZone * 10;
                //! current MixingOut.
                MixingFlag = 1 + iTurboZone * 10;
            }
            //! upstream to downstream: transfer data from MixingOut to MixingIn.
            else if (TargetZone == iTurboZone - 1)
            {
                if (bcName == MixingPlaneOut) continue;
                TargetFlag = 1 + TargetZone * 10;
                MixingFlag = 0 + iTurboZone * 10;
            }
        }

        if (bcName == MixingPlaneIn || bcName == MixingPlaneOut)
        {
            //! compute variable perturbation
            //! when computing perturbation, notice that the direction of vn and vt are opposite.
            //! because flow in and flow out of the interface.
            for (int iSpan = 0; iSpan < nSpanSection; iSpan++)
            {
                for (int m = 0; m < nEquation; m++)
                {
                    ExtraDq[m] = 0.0;
                }

                //! interpolation
                for (int m = 0; m < nEquation; m++)
                {
                    dqSpanIn[MixingFlag][m][iSpan] = q_SpanTarget[TargetFlag][m][iSpan] - q_SpanCurrent[MixingFlag][m][iSpan];
                }

                dqSpanIn[MixingFlag][IU][iSpan] = -q_SpanTarget[TargetFlag][IU][iSpan] - q_SpanCurrent[MixingFlag][IU][iSpan];
                dqSpanIn[MixingFlag][IW][iSpan] = -q_SpanTarget[TargetFlag][IW][iSpan] - q_SpanCurrent[MixingFlag][IW][iSpan];

                //! extrapolation
                vector<int> *faceIndex = bcRegion->GetFaceIndex();
                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    int iFace = *iter;

                    int re = rightCellofFace[iFace];
                    int le = leftCellofFace[iFace];

                    RDouble r_Face = sqrt(yfc[iFace] * yfc[iFace] + zfc[iFace] * zfc[iFace]);

                    RDouble nts = (nzs[iFace] * yfc[iFace] - nys[iFace] * zfc[iFace]) / r_Face;
                    RDouble nrs = (nys[iFace] * yfc[iFace] + nzs[iFace] * zfc[iFace]) / r_Face;

                    RDouble Velocity_ti = (q[IV][le] * zfc[iFace] - q[IW][le] * yfc[iFace]) / r_Face;
                    RDouble Velocity_ri = (q[IV][le] * yfc[iFace] + q[IW][le] * zfc[iFace]) / r_Face;

                    RDouble Vni     =  q[IU][le] * nxs[iFace] + Velocity_ri * nrs;
                    RDouble Vthetai =  Velocity_ti;
                    RDouble Vti     = -q[IU][le] * nrs        + Velocity_ri * nxs[iFace];

                    //! extrapolating q[m][le] and q on the interface to obtain ExtraDq.
                    if (r_Face > r_target[MixingFlag][iSpan] && r_Face <= r_target[MixingFlag][iSpan + 1])
                    {
                        ExtraDq[IR] += (q_SpanCurrent[MixingFlag][IR][iSpan] - q[IR][le]) * ns[iFace];
                        ExtraDq[IU] += (q_SpanCurrent[MixingFlag][IU][iSpan] - Vni)       * ns[iFace];
                        ExtraDq[IV] += (q_SpanCurrent[MixingFlag][IV][iSpan] - Vthetai)   * ns[iFace];
                        ExtraDq[IW] += (q_SpanCurrent[MixingFlag][IW][iSpan] - Vti)       * ns[iFace];
                        ExtraDq[IP] += (q_SpanCurrent[MixingFlag][IP][iSpan] - q[IP][le]) * ns[iFace];
                    }
                }

                for (int m = 0; m < nEquation; m++)
                {
                    dqSpanEx[MixingFlag][m][iSpan] = ExtraDq[m] / SpanArea[MixingFlag][iSpan];
                }
            }

            for (int iSpan = 0; iSpan < nSpanSection; iSpan++)
            {
                RDouble deltaQ[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

                //! IU refers to vn.
                RDouble vn = q_SpanCurrent[MixingFlag][IU][iSpan];

                RDouble c = sqrt(refGama* q_SpanCurrent[MixingFlag][IP][iSpan] / q_SpanCurrent[MixingFlag][IR][iSpan]);

                RDouble avgRho = half * (q_SpanCurrent[MixingFlag][IR][iSpan] + q_SpanTarget[TargetFlag][IR][iSpan]);
                RDouble avgP   = half * (q_SpanCurrent[MixingFlag][IP][iSpan] + q_SpanTarget[TargetFlag][IP][iSpan]);

                RDouble avgc = sqrt(refGama * (avgP / avgRho));

                RDouble rhoc = half * (q_SpanCurrent[MixingFlag][IR][iSpan] + q_SpanTarget[TargetFlag][IR][iSpan]) * avgc;

                //! perturbation
                if (vn > 0)
                {
                    dcSpan[MixingFlag][0][iSpan] = -avgc * avgc * dqSpanEx[MixingFlag][IR][iSpan] + dqSpanEx[MixingFlag][IP][iSpan];
                    dcSpan[MixingFlag][2][iSpan] = rhoc * dqSpanEx[MixingFlag][IV][iSpan];
                    dcSpan[MixingFlag][3][iSpan] = rhoc * dqSpanEx[MixingFlag][IW][iSpan];
                }
                else
                {
                    dcSpan[MixingFlag][0][iSpan] = -avgc * avgc * dqSpanIn[MixingFlag][IR][iSpan] + dqSpanIn[MixingFlag][IP][iSpan];
                    dcSpan[MixingFlag][2][iSpan] = rhoc * dqSpanIn[MixingFlag][IV][iSpan];
                    dcSpan[MixingFlag][3][iSpan] = rhoc * dqSpanIn[MixingFlag][IW][iSpan];
                }

                if (vn + c > 0)
                {
                    dcSpan[MixingFlag][1][iSpan] = rhoc * dqSpanEx[MixingFlag][IU][iSpan] + dqSpanEx[MixingFlag][IP][iSpan];
                }
                else
                {
                    dcSpan[MixingFlag][1][iSpan] = rhoc * dqSpanIn[MixingFlag][IU][iSpan] + dqSpanIn[MixingFlag][IP][iSpan];
                }

                if (vn - c > 0)
                {
                    dcSpan[MixingFlag][4][iSpan] = -rhoc * dqSpanEx[MixingFlag][IU][iSpan] + dqSpanEx[MixingFlag][IP][iSpan];
                }
                else
                {
                    dcSpan[MixingFlag][4][iSpan] = -rhoc * dqSpanIn[MixingFlag][IU][iSpan] + dqSpanIn[MixingFlag][IP][iSpan];
                }

                //! IU, IV, IW refers to vn, vtheta, vt.
                deltaQ[IR] = (-dcSpan[MixingFlag][0][iSpan] + half * (dcSpan[MixingFlag][1][iSpan] + dcSpan[MixingFlag][4][iSpan])) / (avgc * avgc);
                deltaQ[IU] = half * (dcSpan[MixingFlag][1][iSpan] - dcSpan[MixingFlag][4][iSpan]) / rhoc;
                deltaQ[IV] = dcSpan[MixingFlag][2][iSpan] / rhoc;
                deltaQ[IW] = dcSpan[MixingFlag][3][iSpan] / rhoc;
                deltaQ[IP] = half * (dcSpan[MixingFlag][1][iSpan] + dcSpan[MixingFlag][4][iSpan]);

                for (int m = 0; m < nEquation; m++)
                {
                    //! qSpanNeighbor[MixingFlag][m][iSpan] += deltaQ[m];
                    qSpanNeighbor[MixingFlag][m][iSpan] = q_SpanCurrent[MixingFlag][m][iSpan] + deltaQ[m];
                }
            }
        }
    }
    delete [] ExtraDq;    ExtraDq = nullptr;
}

//! function to set radial profile of mixingin and mixingout.
//! treat mixingin as a kind of inlet and treat mixing as a kind of outlet.
void NSSolverUnstruct::SetMixingPlaneData(Grid *grid)
{
    UnstructGrid *gridUnstruct = UnstructGridCast(grid);
    UnstructBCSet *unstructBCSet = gridUnstruct->GetUnstructBCSet();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refGama = parameters->GetRefGama();
    RDouble gama1 = refGama - 1;

    RDouble **q = reinterpret_cast< RDouble ** > (gridUnstruct->GetDataPtr("q"));

    RDouble **r_ori    = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Face"));
    RDouble **r_target = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("Radius_Span"));

    //! nSpan is the number of spanwise point to do data average
    //! also interpolate data on these points.
    //! nSpan is set to 50 temporary.
    int nSpanSection = PHSPACE::GlobalDataBase::GetIntParaFromDB("nSpanSection");
    RDouble nSpan = (RDouble)(nSpanSection);

    RDouble *nxs, *nys, *nzs, *ns;

    nxs = gridUnstruct->GetFaceNormalX();
    nys = gridUnstruct->GetFaceNormalY();
    nzs = gridUnstruct->GetFaceNormalZ();
    ns = gridUnstruct->GetFaceArea();

    int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");

    string MixingPlane[100];
    GlobalDataBase::GetData("MixingPlane", &MixingPlane, PHSTRING, 2 * nTurboZone - 2);

    RDouble *xfc = gridUnstruct->GetFaceCenterX();
    RDouble *yfc = gridUnstruct->GetFaceCenterY();
    RDouble *zfc = gridUnstruct->GetFaceCenterZ();

    RDouble *xcc = gridUnstruct->GetCellCenterX();
    RDouble *ycc = gridUnstruct->GetCellCenterY();
    RDouble *zcc = gridUnstruct->GetCellCenterZ();

    RDouble refDensity  = parameters->GetRefDimensionalDensity();
    RDouble refPressure = parameters->GetRefDimensionalPressure();

    int **face2nodeArray = gridUnstruct->GetFace2NodeArray();
    int * node_number_of_each_face = gridUnstruct->GetNodeNumberOfEachFace();

    int *leftCellofFace = gridUnstruct->GetLeftCellOfFace();
    int *rightCellofFace = gridUnstruct->GetRightCellOfFace();
    int nTotalCell = gridUnstruct->GetNTotalCell();

    int nBCRegion = unstructBCSet->GetnBCRegion();

    int nEquation = GetNumberOfEquations();

    using namespace IDX;

    RDouble ***q_Span        = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("q_Span"));
    RDouble ***qSpanNeighbor = reinterpret_cast< RDouble *** > (gridUnstruct->GetDataPtr("qSpanNeighbor"));

    RDouble **nxsSpan = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("nxsSpan"));
    RDouble **ntsSpan = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("ntsSpan"));
    RDouble **nrsSpan = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("nrsSpan"));

    int FlowDirection = -1;

    for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        vector<int> *faceIndex = bcRegion->GetFaceIndex();

        int bcType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();

        //! Parallel
        int iTurboZone = gridUnstruct->GetOrdinaryGridIndex();
        //! Serial
        if (iTurboZone == -1)
        {
            iTurboZone = gridUnstruct->GetZoneID();
        }

        string MixingPlaneIn, MixingPlaneOut;

        if (iTurboZone == 0)
        {
            MixingPlaneOut = MixingPlane[0];
        }
        else if (iTurboZone == nTurboZone)
        {
            MixingPlaneIn = MixingPlane[2 * iTurboZone - 1];
        }
        else
        {
            MixingPlaneIn  = MixingPlane[2 * iTurboZone - 1];
            MixingPlaneOut = MixingPlane[2 * iTurboZone];
        }

        if (bcName == MixingPlaneIn || bcName == MixingPlaneOut)
        {
            if (bcName == MixingPlaneIn)
            {
                FlowDirection = 0;
            }
            if (bcName == MixingPlaneOut)
            {
                FlowDirection = 1;
            }
            int MixingFlag = FlowDirection + iTurboZone * 10;

            for (int iSpan = 0; iSpan < nSpanSection; iSpan++)
            {
                for (int m = 0; m < nEquation; m++)
                {
                    q_Span[MixingFlag][m][iSpan] = qSpanNeighbor[MixingFlag][m][iSpan];
                }

                //! return to vx, vtheta, vr.
                RDouble Temp_vx = q_Span[MixingFlag][IU][iSpan] * nxsSpan[MixingFlag][iSpan] - q_Span[MixingFlag][IW][iSpan] * nrsSpan[MixingFlag][iSpan];
                RDouble Temp_vr = q_Span[MixingFlag][IU][iSpan] * nrsSpan[MixingFlag][iSpan] + q_Span[MixingFlag][IW][iSpan] * nxsSpan[MixingFlag][iSpan];

                q_Span[MixingFlag][IU][iSpan] = Temp_vx;
                q_Span[MixingFlag][IW][iSpan] = Temp_vr;
            }

            for (int iSpan = 0; iSpan < nSpanSection; iSpan++)
            {

                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    //! iFace is the face number in the set of faceIndex.
                    int iFace = *iter;

                    int le = leftCellofFace[iFace];
                    int re = rightCellofFace[iFace];

                    RDouble Radius = sqrt(yfc[iFace] * yfc[iFace] + zfc[iFace] * zfc[iFace]);

                    if (Radius > r_target[MixingFlag][iSpan] && Radius <= r_target[MixingFlag][iSpan + 1])
                    {
                        RDouble Temp_vx = q_Span[MixingFlag][IU][iSpan];
                        RDouble Temp_vr = q_Span[MixingFlag][IW][iSpan];

                        //! vy and vz is different for each cell.
                        RDouble Temp_vy = (Temp_vr * yfc[iFace] + q_Span[MixingFlag][IV][iSpan] * zfc[iFace]) / Radius;
                        RDouble Temp_vz = (Temp_vr * zfc[iFace] - q_Span[MixingFlag][IV][iSpan] * yfc[iFace]) / Radius;

                        q[IR][re] = q_Span[MixingFlag][IR][iSpan];
                        q[IU][re] = Temp_vx;
                        q[IV][re] = Temp_vy;
                        q[IW][re] = Temp_vz;
                        q[IP][re] = q_Span[MixingFlag][IP][iSpan];
                    }
                }
            }
        }
    }
}

void NSSolverUnstruct::RotateVectorFromInterface(Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
{
    UnstructGrid *gridUnstruct = UnstructGridCast(gridIn);
    InterfaceInfo *interfaceInformation = gridUnstruct->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    UnstructBCSet *unstructBCSet = gridUnstruct->GetUnstructBCSet();

    RDouble **rotNSgradValueX = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("rotNSgradValueX"));
    RDouble **rotNSgradValueY = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("rotNSgradValueY"));
    RDouble **rotNSgradValueZ = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("rotNSgradValueZ"));

    int iNeighborZone = interfaceInformation->FindIthNeighbor(neighborZoneIndex);

    if (iNeighborZone == -1)
    {
        return;
    }

    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();
    int *interFace2BoundaryFace = interfaceInformation->GetInterFace2BoundaryFace();

    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");
    RDouble rotationAngle = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("rotationAngle");

    rotationAngle = rotationAngle / 180.0 * PI;

    RDouble *nxs, *nys, *nzs, *ns;

    nxs = gridUnstruct->GetFaceNormalX();
    nys = gridUnstruct->GetFaceNormalY();
    nzs = gridUnstruct->GetFaceNormalZ();
    ns  = gridUnstruct->GetFaceArea();

    RDouble *xfc = gridUnstruct->GetFaceCenterX();
    RDouble *yfc = gridUnstruct->GetFaceCenterY();
    RDouble *zfc = gridUnstruct->GetFaceCenterZ();

    RDouble *xcc = gridUnstruct->GetCellCenterX();
    RDouble *ycc = gridUnstruct->GetCellCenterY();
    RDouble *zcc = gridUnstruct->GetCellCenterZ();

    int *leftCellofFace = gridUnstruct->GetLeftCellOfFace();

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");

    if (referenceFrame == ROTATIONAL_FRAME)
    {
        int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
        RDouble PeriodicRotationAngle[100];
        GlobalDataBase::GetData("PeriodicRotationAngle", &PeriodicRotationAngle, PHDOUBLE, nTurboZone);

        //! Parallel
        int iTurboZone = gridUnstruct->GetOrdinaryGridIndex();
        //! Serial
        if (iTurboZone == -1)
        {
            iTurboZone = gridUnstruct->GetZoneID();
        }

        PeriodicRotationAngle[iTurboZone] = PeriodicRotationAngle[iTurboZone] * PI / 180.0;
        rotationAngle = PeriodicRotationAngle[iTurboZone];
    }

    if (nEquation >= 4)
    {
        RDouble **fieldRecvNSX = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarX"));
        RDouble **fieldRecvNSY = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarY"));
        RDouble **fieldRecvNSZ = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarZ"));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                int t1,le;
                le = leftCellofFace[iFace];
                gridUnstruct->GetTargetIndex(iFace, 1, t1);

                int iBFace = interFace2BoundaryFace[iFace];
                UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iBFace]);
                string bcName = bcRegion->GetBCName();

                RDouble theta;

                if (referenceFrame == ROTATIONAL_FRAME)
                {
                    int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
                    string Periodic_Name[100];
                    GlobalDataBase::GetData("Periodic_Name", &Periodic_Name, PHSTRING, 2 * nTurboZone);

                    //! Parallel
                    int iTurboZone = gridUnstruct->GetOrdinaryGridIndex();
                    //! Serial
                    if (iTurboZone == -1)
                    {
                        iTurboZone = gridUnstruct->GetZoneID();
                    }

                    if (bcName != Periodic_Name[2 * iTurboZone] && bcName != Periodic_Name[2 * iTurboZone + 1])
                    {
                        continue;
                    }

                    if (bcName == Periodic_Name[2 * iTurboZone])
                    {
                        theta = 2 * PI - rotationAngle;
                    }
                    else if (bcName == Periodic_Name[2 * iTurboZone + 1])
                    {
                        theta = rotationAngle;
                    }

                    //! for velocity gradient, both direction and nys,nzs direction need rotate.
                    rotNSgradValueX[IV][t1] = fieldRecvNSX[IV][t1] * cos(theta) - fieldRecvNSX[IW][t1] * sin(theta);
                    rotNSgradValueY[IV][t1] = fieldRecvNSY[IV][t1] * pow(cos(theta), 2) - fieldRecvNSZ[IV][t1] * sin(theta) * cos(theta)
                        - fieldRecvNSY[IW][t1] * sin(theta) * cos(theta) + fieldRecvNSZ[IW][t1] * pow(sin(theta), 2);
                    rotNSgradValueZ[IV][t1] = fieldRecvNSY[IV][t1] * sin(theta) * cos(theta) + fieldRecvNSZ[IV][t1] * pow(cos(theta), 2)
                        - fieldRecvNSY[IW][t1] * pow(sin(theta), 2) - fieldRecvNSZ[IW][t1] * sin(theta) * cos(theta);

                    rotNSgradValueX[IW][t1] = fieldRecvNSX[IV][t1] * sin(theta) + fieldRecvNSX[IW][t1] * cos(theta);
                    rotNSgradValueY[IW][t1] = fieldRecvNSY[IV][t1] * sin(theta) * cos(theta) - fieldRecvNSZ[IV][t1] * pow(sin(theta), 2)
                        + fieldRecvNSY[IW][t1] * pow(cos(theta), 2) - fieldRecvNSZ[IW][t1] * sin(theta) * cos(theta);
                    rotNSgradValueZ[IW][t1] = fieldRecvNSY[IV][t1] * pow(sin(theta), 2) + fieldRecvNSZ[IV][t1] * sin(theta) * cos(theta)
                        + fieldRecvNSY[IW][t1] * sin(theta) * cos(theta) + fieldRecvNSZ[IW][t1] * pow(cos(theta), 2);

                    fieldRecvNSX[IV][t1] = rotNSgradValueX[IV][t1];
                    fieldRecvNSY[IV][t1] = rotNSgradValueY[IV][t1];
                    fieldRecvNSZ[IV][t1] = rotNSgradValueZ[IV][t1];

                    fieldRecvNSX[IW][t1] = rotNSgradValueX[IW][t1];
                    fieldRecvNSY[IW][t1] = rotNSgradValueY[IW][t1];
                    fieldRecvNSZ[IW][t1] = rotNSgradValueZ[IW][t1];

                    for (int m = 0; m < nEquation; m++)
                    {
                        if (m == 2 || m == 3)
                        {
                            continue;
                        }

                        if (bcName == Periodic_Name[2 * iTurboZone])
                        {
                            rotNSgradValueY[m][t1] = fieldRecvNSY[m][t1] * cos(2 * PI - rotationAngle) - fieldRecvNSZ[m][t1] * sin(2 * PI - rotationAngle);
                            rotNSgradValueZ[m][t1] = fieldRecvNSY[m][t1] * sin(2 * PI - rotationAngle) + fieldRecvNSZ[m][t1] * cos(2 * PI - rotationAngle);
                        }
                        else if (bcName == Periodic_Name[2 * iTurboZone + 1])
                        {
                            rotNSgradValueY[m][t1] = fieldRecvNSY[m][t1] * cos(rotationAngle) - fieldRecvNSZ[m][t1] * sin(rotationAngle);
                            rotNSgradValueZ[m][t1] = fieldRecvNSY[m][t1] * sin(rotationAngle) + fieldRecvNSZ[m][t1] * cos(rotationAngle);
                        }

                        fieldRecvNSY[m][t1] = rotNSgradValueY[m][t1];
                        fieldRecvNSZ[m][t1] = rotNSgradValueZ[m][t1];
                    }
                }

                else
                {
                    if (bcName != "Periodic_up" && bcName != "Periodic_down")
                    {
                        continue;
                    }

                    if (bcName == "Periodic_up")
                    {
                        theta = 2 * PI - rotationAngle;
                    }
                    else if (bcName == "Periodic_down")
                    {
                        theta = rotationAngle;
                    }

                    //! for velocity gradient, both direction and nys,nzs direction need rotate.
                    rotNSgradValueX[IV][t1] = fieldRecvNSX[IV][t1] * cos(theta) - fieldRecvNSX[IW][t1] * sin(theta);
                    rotNSgradValueY[IV][t1] = fieldRecvNSY[IV][t1] * pow(cos(theta), 2) - fieldRecvNSZ[IV][t1] * sin(theta) * cos(theta)
                        - fieldRecvNSY[IW][t1] * sin(theta) * cos(theta) + fieldRecvNSZ[IW][t1] * pow(sin(theta), 2);
                    rotNSgradValueZ[IV][t1] = fieldRecvNSY[IV][t1] * sin(theta) * cos(theta) + fieldRecvNSZ[IV][t1] * pow(cos(theta), 2)
                        - fieldRecvNSY[IW][t1] * pow(sin(theta), 2) - fieldRecvNSZ[IW][t1] * sin(theta) * cos(theta);

                    rotNSgradValueX[IW][t1] = fieldRecvNSX[IV][t1] * sin(theta) + fieldRecvNSX[IW][t1] * cos(theta);
                    rotNSgradValueY[IW][t1] = fieldRecvNSY[IV][t1] * sin(theta) * cos(theta) - fieldRecvNSZ[IV][t1] * pow(sin(theta), 2)
                        + fieldRecvNSY[IW][t1] * pow(cos(theta), 2) - fieldRecvNSZ[IW][t1] * sin(theta) * cos(theta);
                    rotNSgradValueZ[IW][t1] = fieldRecvNSY[IV][t1] * pow(sin(theta), 2) + fieldRecvNSZ[IV][t1] * sin(theta) * cos(theta)
                        + fieldRecvNSY[IW][t1] * sin(theta) * cos(theta) + fieldRecvNSZ[IW][t1] * pow(cos(theta), 2);

                    fieldRecvNSX[IV][t1] = rotNSgradValueX[IV][t1];
                    fieldRecvNSY[IV][t1] = rotNSgradValueY[IV][t1];
                    fieldRecvNSZ[IV][t1] = rotNSgradValueZ[IV][t1];

                    fieldRecvNSX[IW][t1] = rotNSgradValueX[IW][t1];
                    fieldRecvNSY[IW][t1] = rotNSgradValueY[IW][t1];
                    fieldRecvNSZ[IW][t1] = rotNSgradValueZ[IW][t1];

                    for (int m = 0; m < nEquation; m++)
                    {
                        if (m == 2 || m == 3)
                        {
                            continue;
                        }

                        if (bcName == "Periodic_up")
                        {
                            rotNSgradValueY[m][t1] = fieldRecvNSY[m][t1] * cos(2 * PI - rotationAngle) - fieldRecvNSZ[m][t1] * sin(2 * PI - rotationAngle);
                            rotNSgradValueZ[m][t1] = fieldRecvNSY[m][t1] * sin(2 * PI - rotationAngle) + fieldRecvNSZ[m][t1] * cos(2 * PI - rotationAngle);
                        }
                        else if (bcName == "Periodic_down")
                        {
                            rotNSgradValueY[m][t1] = fieldRecvNSY[m][t1] * cos(rotationAngle) - fieldRecvNSZ[m][t1] * sin(rotationAngle);
                            rotNSgradValueZ[m][t1] = fieldRecvNSY[m][t1] * sin(rotationAngle) + fieldRecvNSZ[m][t1] * cos(rotationAngle);
                        }

                        fieldRecvNSY[m][t1] = rotNSgradValueY[m][t1];
                        fieldRecvNSZ[m][t1] = rotNSgradValueZ[m][t1];
                    }
                }
            }
        }
    }

    //! for temperature, only one equation.
    if (nEquation > 0 && nEquation < 4)
    {
        RDouble **fieldRecvTX = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTemperatureX"));
        RDouble **fieldRecvTY = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTemperatureY"));
        RDouble **fieldRecvTZ = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTemperatureZ"));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                int t1;
                gridUnstruct->GetTargetIndex(iFace, 1, t1);

                int iBFace = interFace2BoundaryFace[iFace];
                UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iBFace]);
                string bcName = bcRegion->GetBCName();

                for (int m = 0; m < nEquation; m++)
                {
                    RDouble Ty, Tz = 0;

                    if (referenceFrame == ROTATIONAL_FRAME)
                    {
                        int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
                        string Periodic_Name[100];
                        GlobalDataBase::GetData("Periodic_Name", &Periodic_Name, PHSTRING, 2 * nTurboZone);

                        //! Parallel
                        int iTurboZone = gridUnstruct->GetOrdinaryGridIndex();
                        //! Serial
                        if (iTurboZone == -1)
                        {
                            iTurboZone = gridUnstruct->GetZoneID();
                        }

                        if (bcName == Periodic_Name[2 * iTurboZone] || bcName == Periodic_Name[2 * iTurboZone + 1])
                        {
                            if (bcName == Periodic_Name[2 * iTurboZone])
                            {
                                Ty = fieldRecvTY[m][t1] * cos(2 * PI - rotationAngle) - fieldRecvTZ[m][t1] * sin(2 * PI - rotationAngle);
                                Tz = fieldRecvTY[m][t1] * sin(2 * PI - rotationAngle) + fieldRecvTZ[m][t1] * cos(2 * PI - rotationAngle);
                            }
                            else if (bcName == Periodic_Name[2 * iTurboZone + 1])
                            {
                                Ty = fieldRecvTY[m][t1] * cos(rotationAngle) - fieldRecvTZ[m][t1] * sin(rotationAngle);
                                Tz = fieldRecvTY[m][t1] * sin(rotationAngle) + fieldRecvTZ[m][t1] * cos(rotationAngle);
                            }
                            fieldRecvTY[m][t1] = Ty;
                            fieldRecvTZ[m][t1] = Tz;
                        }
                    }
                    else
                    {
                        if (bcName == "Periodic_up" || bcName == "Periodic_down")
                        {
                            if (bcName == "Periodic_up")
                            {
                                Ty = fieldRecvTY[m][t1] * cos(2 * PI - rotationAngle) - fieldRecvTZ[m][t1] * sin(2 * PI - rotationAngle);
                                Tz = fieldRecvTY[m][t1] * sin(2 * PI - rotationAngle) + fieldRecvTZ[m][t1] * cos(2 * PI - rotationAngle);
                            }
                            else if (bcName == "Periodic_down")
                            {
                                Ty = fieldRecvTY[m][t1] * cos(rotationAngle) - fieldRecvTZ[m][t1] * sin(rotationAngle);
                                Tz = fieldRecvTY[m][t1] * sin(rotationAngle) + fieldRecvTZ[m][t1] * cos(rotationAngle);
                            }
                            fieldRecvTY[m][t1] = Ty;
                            fieldRecvTZ[m][t1] = Tz;
                        }
                    }
                }
            }
        }
    }
}


FieldProxy * NSSolverUnstruct::CreateFieldProxy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    int nEquation = GetNumberOfEquations();

    RDouble **field = NewPointer2<RDouble>(nEquation, nTotal);

    FieldProxy *fieldProxy = new FieldProxy();

    fieldProxy->SetField_UNS(field, true);

    return fieldProxy;
}

FieldProxy * NSSolverUnstruct::GetFieldProxy(Grid *gridIn, const string &field_name)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **field = reinterpret_cast< RDouble ** > (grid->GetDataPtr(field_name));

    FieldProxy *field_proxy = new FieldProxy();

    field_proxy->SetField_UNS(field);

    return field_proxy;
}

FieldProxy * NSSolverUnstruct::GetResidualProxy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    FieldProxy *res_proxy = new FieldProxy();

    res_proxy->SetField_UNS(res);

    return res_proxy;
}

void NSSolverUnstruct::RecoverResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
#ifdef USE_CUDA
    int nEquation_gpu = GetNumberOfEquations();
    CallGPURecoverResidual(gridIn, nEquation_gpu);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    int nEquation = GetNumberOfEquations();

    RDouble **res = reinterpret_cast<RDouble **> (grid->GetDataPtr("res"));

    if (grid->GetLevel() != 0) 
    {
        RDouble **rhs = rhsProxy->GetField_UNS();
        for (int m = 0; m < nEquation; ++ m)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                res[m][iCell] = rhs[m][iCell];
            }
        }
    }
}

void NSSolverUnstruct::StoreRhsByResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
#ifdef USE_CUDA
    int nEquation_gpu = GetNumberOfEquations();
    CallGPUStoreRhsByResidual(gridIn, nEquation_gpu);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    int nEquation = GetNumberOfEquations();

    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));
    RDouble **rhs = rhsProxy->GetField_UNS();

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            rhs[m][iCell] = res[m][iCell];
        }
    }
}

void NSSolverUnstruct::InitResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    int nEquation = GetNumberOfEquations();

    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    RDouble **rhs = rhsProxy->GetField_UNS();
    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            res[m][iCell] = - rhs[m][iCell];
        }
    }
}

void NSSolverUnstruct::FreeGradientProxy(Grid *grid)
{

}

void NSSolverUnstruct::DumpLeakageInformation(Grid *gridIn)
{
    int nSpecies = gas->GetNumberOfSpecies();
    string *speciesNames = gas->GetNameOfSpecies();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nm = parameters->GetNSEquationNumber();

    int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep") ;
    int innerStep = GlobalDataBase::GetIntParaFromDB("innstep");
    RDouble physicalTimeStepDimensional = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStepDimensional");
    RDouble totalLeakageVolume = GlobalDataBase::GetDoubleParaFromDB("totalLeakageVolume");

    RDouble monitorThresholdValue = 0.02;
    if (GlobalDataBase::IsExist("monitorThresholdValue", PHDOUBLE, 1))
    {
        monitorThresholdValue = GlobalDataBase::GetDoubleParaFromDB("monitorThresholdValue") ;
    }

    int nEquation = GetNumberOfEquations();
    RDouble *primVars = new RDouble [nEquation];

    RDouble *tempVariables = new RDouble [30];
    RDouble *tempTotalVariables = new RDouble [30];
    for (int i = 0; i < 30; i++)
    {
        tempVariables[i] = 0.0;
        tempTotalVariables[i] = 0.0;
    }
    RDouble *monitorVariables = reinterpret_cast< RDouble * > (GlobalDataBase::GetDataPtr("monitorVariables" ));

    monitorVariables[0] = outnstep * 1.0;
    monitorVariables[1] = outnstep * physicalTimeStepDimensional;
    RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
    RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refVelocity = refDimensionalSonicSpeed * refMachNumber;
    using namespace PHMPI;
    if (innerStep == 1)
    {
        int nZones = GetNumberofGlobalZones();
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int zoneProcessorID = GetZoneProcessorID(iZone);
            int currentProcessorID = GetCurrentProcessorID();
            if (currentProcessorID == zoneProcessorID)
            {
                UnstructGrid *grid = UnstructGridCast(PHSPACE::GetGrid(iZone, 0));
                int nTotalCell = grid->GetNTotalCell();
                RDouble *vol   = grid->GetCellVolume();
                int *leftCellofFace = grid->GetLeftCellOfFace();
                RDouble *xfn  = grid->GetFaceNormalX();
                RDouble *yfn  = grid->GetFaceNormalY();
                RDouble *zfn  = grid->GetFaceNormalZ();
                RDouble *area = grid->GetFaceArea();

                RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

                UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
                int nBCRegion = unstructBCSet->GetnBCRegion();
                for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
                {
                    UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
                    int bcType = bcRegion->GetBCType();
                    if (bcType == PHENGLEI::PRESSURE_INLET)
                    {
                        Data_Param *bcData = bcRegion->GetBCParamDataBase();
                        RDouble totalPressure = 0.0;
                        RDouble totalTemperature = 0.0;
                        bcData->GetData("totalPressure", &totalPressure, PHDOUBLE, 1);
                        bcData->GetData("totalTemperature", &totalTemperature, PHDOUBLE, 1);

                        RDouble initMassFraction[MAX_SPECIES_NUM] = {0};
                        bcData->GetData("initMassFraction", &initMassFraction, PHDOUBLE, nSpecies);

                        for (int iEquation = 0; iEquation < nm; iEquation++)
                        {
                            primVars[iEquation] = 0.0;
                        }
                        for (int iEquation = nm; iEquation < nEquation; iEquation++)
                        {
                            primVars[iEquation] = initMassFraction[iEquation - nm];
                        }
                        RDouble oMass = 0.0, gasR = 0.0;
                        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
                        gasR = rjmk * oMass;

                        RDouble ro = totalPressure / (totalTemperature * gasR);

                        vector<int> *faceIndex = bcRegion->GetFaceIndex();
                        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                        {
                            // iFace is the face number in the set of faceIndex.
                            int iFace = *iter;
                            int le = leftCellofFace[iFace];

                            RDouble normalVelocity = refVelocity * ABS(q[IU][le] * xfn[iFace] + q[IV][le] * yfn[iFace] + q[IW][le] * zfn[iFace]);
                            tempVariables[2] += (refDimensionalDensity * q[IDX::IR][le] * normalVelocity * physicalTimeStepDimensional * area[iFace]) / ro;
                            tempVariables[4] += (refDimensionalDensity * q[IDX::IR][le] * normalVelocity * area[iFace]) / ro;
                        }
                    }
                }

                for (int s = 0; s < nSpecies; ++ s)
                {
                    monitorVariables[5 + s] = 0.0;
                    for (int iCell = 0; iCell < nTotalCell; iCell++)
                    {
                        if(q[nm + s][iCell] > monitorThresholdValue)
                        {
                            tempVariables[5 + s] += vol[iCell];
                        }
                    }
                }
            }
        }

        PH_AllReduce(tempVariables, tempTotalVariables, 30, PH_SUM);
        monitorVariables[2] += tempTotalVariables[2];
        monitorVariables[3] = totalLeakageVolume - monitorVariables[2];
        monitorVariables[4] = tempTotalVariables[4];

        RDouble maxLeakageTime = LARGE;
        if (GlobalDataBase::IsExist("maxLeakageTime", PHDOUBLE, 1))
        {
            maxLeakageTime = GlobalDataBase::GetDoubleParaFromDB("maxLeakageTime");
        }

        RDouble currentTime = outnstep * physicalTimeStepDimensional;
        if(monitorVariables[3] < 0 || maxLeakageTime < currentTime)
        {
            for (int iZone = 0; iZone < nZones; ++ iZone)
            {
                int zoneProcessorID = GetZoneProcessorID(iZone);
                int currentProcessorID = GetCurrentProcessorID();
                if (currentProcessorID == zoneProcessorID)
                {
                    UnstructGrid *grid = UnstructGridCast(PHSPACE::GetGrid(iZone, 0));
                    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
                    int nBCRegion = unstructBCSet->GetnBCRegion();
                    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
                    {
                        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
                        int bcType = bcRegion->GetBCType();
                        if (bcType == PHENGLEI::PRESSURE_INLET)
                        {
                            bcRegion->SetBCType(PHENGLEI::SOLID_SURFACE);
                        }
                    }
                }
            }
        }

        for (int i = 5; i < 30; i++)
        {
            monitorVariables[i] = tempTotalVariables[i];
        }

        if (GetCurrentProcessorID() == GetServerProcessorID())
        {
            ostringstream resData;
            resData << setiosflags(ios::left);
            resData << setprecision(5);
            resData << setiosflags(ios::scientific);
            resData << setiosflags(ios::showpoint);

            resData << setw(7) << outnstep << "    ";
            resData << monitorVariables[1] << "   ";
            resData << monitorVariables[2] << "   ";
            resData << monitorVariables[3] << "   ";
            resData << monitorVariables[4] << "   ";
            for (int s = 0; s < nSpecies; ++ s)
            {
                resData << monitorVariables[5 + s] << "   ";
            }
            resData << endl;

            //PrintToWindow(resData);
            string resFileName = "results/leakage.dat";
            ios_base::openmode openMode = ios_base::out | ios_base::app;

            fstream file;
            OpenFile(file, resFileName, openMode);
            ostringstream dumpFile;
            if (IfFileEmpty(file))
            {
                vector<string> title_tecplot;
                title_tecplot.push_back("Title=\"THE RESIDUAL\"");
                title_tecplot.push_back("Variables=");
                title_tecplot.push_back("\"iter\"");
                title_tecplot.push_back("\"physicalTime\"");
                title_tecplot.push_back("\"totalLeakage\"");
                title_tecplot.push_back("\"lastLeakageVolume\"");
                title_tecplot.push_back("\"MassLeakageRate\"");

                for (int s = 0; s < nSpecies; ++ s)
                {
                    title_tecplot.push_back("\"" +  speciesNames[s] + "\"");
                }

                for (std::size_t i = 0; i < title_tecplot.size(); ++i)
                {
                    dumpFile << title_tecplot[i] << "\n";
                }
            }
            dumpFile << resData.str();
            WriteASCIIFile(file, dumpFile.str());
            CloseFile(file);
        }
    }
    delete [] tempVariables;    tempVariables = nullptr;
    delete [] tempTotalVariables;    tempTotalVariables = nullptr;
    delete [] primVars;    primVars = nullptr;
}


//! Drive functions to calculate residuals for different orders and dimensions.
void NSSolverUnstruct::Turb_Sengy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    if (viscousType >= TWO_EQU)
    {
        int nrokplus = 1;
        GlobalDataBase::GetData("nrokplus", &nrokplus, PHINT, 1);

        if (nrokplus > 0)
        {
            RDouble **res  = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));
            RDouble *sengy = reinterpret_cast< RDouble * > (grid->GetDataPtr("sengy"  ));
            for (int iCell = 0; iCell < nTotalCell; ++ iCell)
            {
                res[4][iCell] += sengy[iCell];
            }
        }
    }
}

//! Load flow variables stored in grid to q.
void NSSolverUnstruct::LoadQ(Grid *gridIn, FieldProxy *qProxy)
{
#ifdef USE_CUDA
    Param_NSSolverUnstruct *parameters_gpu = GetControlParameters();
    CallGPULoadQ(gridIn, parameters_gpu);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = qProxy->GetField_UNS();

    RDouble **qold = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();

    for (int m = 0; m < nLaminar; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            q[m][iCell] = qold[m][iCell];
        }
    }
}

//! Restrict defect from fine to coarse grid.
void NSSolverUnstruct::RestrictDefect(Grid *fgrid_in, Grid *cgrid_in)
{
    UnstructGrid *fgrid = UnstructGridCast(fgrid_in);
    UnstructGrid *cgrid = UnstructGridCast(cgrid_in);

    int  fnTotalCell = fgrid->GetNTotalCell();
    int *cell2coarsegridcell    = fgrid->GetCell2CoarseGridCell();

    ZeroResiduals(cgrid);

    RDouble **cres = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("res"));
    RDouble **fres = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("res"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();

    for (int m = 0; m < nLaminar; ++ m)
    {
        for (int iCell = 0; iCell < fnTotalCell; ++ iCell)
        {
            cres[m][ cell2coarsegridcell[iCell] ] -= fres[m][iCell];
        }
    }
}

//! Put correction back on the coarse grid.
void NSSolverUnstruct::PutCorrectionBack(Grid *gridIn, FieldProxy *qProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = qProxy->GetField_UNS();

    RDouble **qold = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();

    for (int m = 0; m < nLaminar; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            qold[m][iCell] += q[m][iCell];
        }
    }
}

//! Correct variable in fine grid using correction in coarse grid.
void NSSolverUnstruct::CorrectFineGrid(Grid *fgrid_in, Grid *cgrid_in)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int mgProlongationType = parameters->GetMgProlongationType();

    if(mgProlongationType == 2)
    {
        CorrectFineGridHighOrder(fgrid_in, cgrid_in);
    }
    else
    {
        CorrectFineGridZeroOrder(fgrid_in, cgrid_in);
    }
}

void NSSolverUnstruct::CorrectFineGridZeroOrder(Grid *fgrid_in, Grid *cgrid_in)
{
    UnstructGrid *fgrid = UnstructGridCast(fgrid_in);
    UnstructGrid *cgrid = UnstructGridCast(cgrid_in);

    int fnTCell = fgrid->GetNTotalCell();

    int *cell2coarsegridcell = fgrid->GetCell2CoarseGridCell();

    RDouble **fq = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("q"));
    RDouble **cq = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("q"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();

    RDouble mgCorrectionLimit = parameters->GetMgCorrectionLimit();

    for (int m = 0; m < nLaminar; ++ m)
    {
        for (int iCell = 0; iCell < fnTCell; ++ iCell)
        {
            int coarseCell = cell2coarsegridcell[iCell];

            RDouble fgridValue = fq[m][iCell];
            RDouble correction = cq[m][coarseCell];
            RDouble DQ         = MIN(fabs(correction), fabs(fgridValue) * mgCorrectionLimit);
            fq[m][iCell]      += SIGN(DQ, correction);
        }
    }

    //FreeGradientProxy(cgrid);
}

void NSSolverUnstruct::CorrectFineGridHighOrder(Grid *fgrid_in, Grid *cgrid_in)
{
    UnstructGrid *fgrid = UnstructGridCast(fgrid_in);
    UnstructGrid *cgrid = UnstructGridCast(cgrid_in);

    int fnTNode = fgrid->GetNTotalNode();
    int fnTCell = fgrid->GetNTotalCell();
    int fnBFace = fgrid->GetNBoundFace();
    int fnTFace = fgrid->GetNTotalFace();

    int *cell2CoarseGridCell = fgrid->GetCell2CoarseGridCell();
    int *fleft_cell_of_face  = fgrid->GetLeftCellOfFace();
    int *fright_cell_of_face = fgrid->GetRightCellOfFace();

    int *cright_cell_of_face = cgrid->GetRightCellOfFace();

    int **face2nodeArray = fgrid->GetFace2NodeArray();
    int *fnode_number_of_each_face = fgrid->GetNodeNumberOfEachFace();

    RDouble **fq = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("q"));
    RDouble **cq = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("q"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();

    RDouble minRatio[5], maxRatio[5];
    for (int m = 0; m < 5; ++ m)
    {
        minRatio[m] = LARGE;
        maxRatio[m] = -LARGE;
    }

    RDouble mgCorrectionLimit = parameters->GetMgCorrectionLimit();

    for (int m = 0; m < nLaminar; ++ m)
    {
        //if(m == 1 || m == 2 || m == 3) continue;
        //! On coarse mesh, average the cell data onto nodes.
        RDouble *fqNode  = new RDouble [ fnTNode ];
        int *fnodeCount = new int [fnTNode];
        for (int iNode = 0; iNode < fnTNode; ++ iNode)
        {
            fqNode[iNode]     = 0.0;
            fnodeCount[iNode] = 0;
        }

        UnstructBCSet **bcRecord = fgrid->GetBCRecord();
        for (int iFace = 0; iFace < fnTFace; ++ iFace)
        {
            int bcType = 0;
            if(iFace < fnBFace)
            {
                bcType = bcRecord[iFace]->GetKey();
            }

            int le = fleft_cell_of_face[iFace];
            int re = fright_cell_of_face[iFace];
            int leftCoarseCell  = cell2CoarseGridCell[le];
            int rightCoarseCell = -1;
            if(iFace >= fnBFace)
            {
                rightCoarseCell = cell2CoarseGridCell[re];
            }
            else if(IsInterface(bcType))
            {
                int cFaceID = iFace;

                //! DCQ = CQ - CQ_old, CQ && CQ_old are obtained by communicating.
                rightCoarseCell = cright_cell_of_face[cFaceID];
            }  

            int nNode = fnode_number_of_each_face[iFace];
            for (int iNode = 0; iNode < nNode; ++ iNode)
            {
                int fnodeID = face2nodeArray[iFace][iNode];

                fqNode[fnodeID] += cq[m][leftCoarseCell];
                fnodeCount[fnodeID] ++;

                if(rightCoarseCell != -1)
                {
                    fqNode[fnodeID] += cq[m][rightCoarseCell];
                    fnodeCount[fnodeID] ++;
                }
            }
        }

        for (int iNode = 0; iNode < fnTNode; ++ iNode)
        {
            fqNode[iNode] /= fnodeCount[iNode];
        }
        delete [] fnodeCount;    fnodeCount = nullptr;

        //! On fine mesh, average the nodes data to cells.
        int *fcellCount   = new int [fnTCell];
        RDouble *fcellData = new RDouble [fnTCell];
        for (int iCell = 0; iCell < fnTCell; ++ iCell)
        {
            fcellCount[iCell] = 0;
            fcellData[iCell]  = 0.0;
        }
        
        for (int iFace = 0; iFace < fnBFace; ++ iFace)
        {
            int le = fleft_cell_of_face[iFace];

            int nNode = fnode_number_of_each_face[iFace];
            for (int iNode = 0; iNode < nNode; ++ iNode)
            {
                int fnodeID = face2nodeArray[iFace][iNode];

                fcellData[le] += fqNode[fnodeID];
                ++ fcellCount[le];
            }
        }

        for (int iFace = fnBFace; iFace < fnTFace; ++ iFace)
        {
            int le = fleft_cell_of_face[iFace];
            int re = fright_cell_of_face[iFace];

            int nNode = fnode_number_of_each_face[iFace];
            for (int iNode = 0; iNode < nNode; ++ iNode)
            {
                int fnodeID = face2nodeArray[iFace][iNode];

                fcellData[le] += fqNode[fnodeID];
                ++ fcellCount[le];

                fcellData[re] += fqNode[fnodeID];
                ++ fcellCount[re];
            }
        }
        delete [] fqNode;    fqNode = nullptr;

        //! Add the correction onto fine grid.
        for (int iCell = 0; iCell < fnTCell; ++ iCell)
        {
            RDouble fgridValue = fq[m][iCell];
            RDouble correction = fcellData[iCell] / fcellCount[iCell];
            RDouble DQ         = MIN(fabs(correction), fabs(fgridValue) * mgCorrectionLimit);
            fq[m][iCell]      += SIGN(DQ, correction);

            if(m == IDX::IP || m == IDX::IR)
            {
                if (fq[m][iCell] <= 1.0e-10)
                {
                    cout << " Negative found after correction, variable = " << m 
                         << ", find grid Value = " << fgridValue << ", correction = " << correction 
                         << ", fcellID = " << iCell << ", ccellID = " << cell2CoarseGridCell[iCell] 
                         << ", level = " << fgrid->GetLevel() << endl;
                }
            }
        }
        delete [] fcellCount;    fcellCount = nullptr;
        delete [] fcellData;    fcellData = nullptr;
    }
}

void NSSolverUnstruct::InterpolatFineGrid(Grid *fineGrid, Grid *coarseGrid)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nEquation = GetNumberOfEquations();

    //! Interpolate vist only if for turbulent simulation is run on coarse grid.
    //!   - 0: do not run turbulent simulation on coarse grid.
    //!   - 1: run turbulent simulation on coarse grid.
    int isRunTurbOnCoarseGrid = 0;

    RDouble **fineQ   = reinterpret_cast<RDouble **> (fineGrid->GetDataPtr("q"));
    RDouble **coarseQ = reinterpret_cast<RDouble **> (coarseGrid->GetDataPtr("q"));

    for (int m = 0; m < nEquation; ++ m)
    {
        InterpolatQ(fineGrid, fineQ[m], coarseGrid, coarseQ[m]);
    }

    RDouble *fineGama   = reinterpret_cast<RDouble *> (fineGrid->GetDataPtr("gama"));
    RDouble *coarseGama = reinterpret_cast<RDouble *> (coarseGrid->GetDataPtr("gama"));

    InterpolatQ(fineGrid, fineGama, coarseGrid, coarseGama);

    int viscousType = parameters->GetViscousType();
    if (viscousType > INVISCID)
    {
        RDouble *fineViscousLaminar   = reinterpret_cast<RDouble *> (fineGrid->GetDataPtr("visl"));
        RDouble *coarseViscousLaminar = reinterpret_cast<RDouble *> (coarseGrid->GetDataPtr("visl"));

        InterpolatQ(fineGrid, fineViscousLaminar, coarseGrid, coarseViscousLaminar);

        if (viscousType > LAMINAR && isRunTurbOnCoarseGrid)
        {
            RDouble *fineViscousTurbulence   = reinterpret_cast<RDouble *> (fineGrid->GetDataPtr("vist"));
            RDouble *coarseViscousTurbulence = reinterpret_cast<RDouble *> (coarseGrid->GetDataPtr("vist"));

            InterpolatQ(fineGrid, fineViscousTurbulence, coarseGrid, coarseViscousTurbulence);
        }
    }
}

//! Load rhs to residual res stored in grid.
void NSSolverUnstruct::LoadResiduals(Grid *gridIn, FieldProxy *rhsProxy)
{
#ifdef USE_CUDA
    int nEquation_gpu = GetNumberOfEquations();
    CallGPULoadResiduals(gridIn, nEquation_gpu);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **rhs = rhsProxy->GetField_UNS();

    int nEquation = GetNumberOfEquations();

    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            res[m][iCell] = - rhs[m][iCell];
        }
    }
}

RDouble NSSolverUnstruct::UnsteadyConvergence(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();

    if (!isUnsteady)
    {
        return zero;
    }

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **qn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n1"));
    RDouble **qn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n2"));
    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    int nLaminar  = parameters->GetLaminarNumber();

    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble refGama = parameters->GetRefGama();

    RDouble *prim0 = new RDouble[nEquation];
    RDouble *prim1 = new RDouble[nEquation];
    RDouble *prim2 = new RDouble[nEquation];

    RDouble *qcsv0 = new RDouble[nEquation];
    RDouble *qcsv1 = new RDouble[nEquation];
    RDouble *qcsv2 = new RDouble[nEquation];

    RDouble Ttr[3] = {0.0}, Tv[3] = {0.0}, Te[3] = {0.0};

    using namespace GAS_SPACE;

    RDouble sum1 = zero;
    RDouble sum2 = zero;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            prim0[m] = q  [m][iCell];
            prim1[m] = qn1[m][iCell];
            prim2[m] = qn2[m][iCell];
        }

        if (nTemperatureModel > 1)
        {
            gas->GetTemperature(prim0, Ttr[0], Tv[0], Te[0]);
            gas->GetTemperature(prim1, Ttr[1], Tv[1], Te[1]);
            gas->GetTemperature(prim2, Ttr[2], Tv[2], Te[2]);
        }

        //! Primitive2Conservative function does not need to input gama in fact.
        //! It will be modified later.
        gas->Primitive2Conservative(prim0, refGama, Tv[0], Te[0], qcsv0);
        gas->Primitive2Conservative(prim1, refGama, Tv[1], Te[1], qcsv1);
        gas->Primitive2Conservative(prim2, refGama, Tv[2], Te[2], qcsv2);

        //! We need nLaminar here.
        for (int m = 0; m < nLaminar; ++ m)
        {
            //! Now res has been converted to dq. It means that residual is not right hand side term,but dq.
            RDouble dq_p =  res[m][iCell];       //! qn+1,p+1 - qn+1,p
            RDouble dq_n =  qcsv0[m] - qcsv1[m]; //! qn+1,p+1 - qn
            sum1 += dq_p * dq_p;
            sum2 += dq_n * dq_n;
        }
    }

    delete [] prim0;    prim0 = nullptr;
    delete [] prim1;    prim1 = nullptr;
    delete [] prim2;    prim2 = nullptr;
    delete [] qcsv0;    qcsv0 = nullptr;
    delete [] qcsv1;    qcsv1 = nullptr;
    delete [] qcsv2;    qcsv2 = nullptr;

    RDouble cvg = sqrt(ABS(sum1 / (sum2 + SMALL)));
    return cvg;
}

void NSSolverUnstruct::UpdateUnsteadyFlow(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();

    if (!isUnsteady) return;

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **qn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n1"));
    RDouble **qn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n2"));

    RDouble **resn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n1"));
    RDouble **resn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n2"));
    RDouble **resTmp = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_tmp"));

    int nEquation = GetNumberOfEquations();

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            qn2[ m ][ iCell ] = qn1[ m ][ iCell ];
            qn1[ m ][ iCell ] = q  [ m ][ iCell ];

            resn2[ m ][ iCell ] = resn1 [ m ][ iCell];
            resn1[ m ][ iCell ] = resTmp[ m ][ iCell];    //! Here the current outnstep is over, the value of the stored resTmp should be assigned to resn1 for the next outnstep. It should be noticed that resTmp only contain the inviscid and viscous flux.
        }
    }

    //! Statistical variables for unsteady simulation.
    if (IsNeedStatistics())
    {
        int nStatisticalStep = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");
        RDouble statisticalTimePeriod = GlobalDataBase::GetDoubleParaFromDB("statisticalTimePeriod");
        RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

        RDouble c1 = 1.0 / (nStatisticalStep * 1.0);
        if(statisticalTimePeriod > 0.0)
        {
            c1 = MAX(c1, physicalTimeStep / statisticalTimePeriod);
        }
        RDouble c2 = 1.0 - c1;

        RDouble **qAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("qAverage"));
        for (int m = 0; m < nEquation; ++ m)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                qAverage[m][iCell] = c2 * qAverage[m][iCell] + c1 * q[m][iCell];
            }
        }
    }

    //! added by zzp 202108, Statistical Reynolds stress for unsteady simulation.
    if (IsNeedReynoldsStressStatistics())
    {
        int nStatisticalStep   = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");
        RDouble statisticalTimePeriod = GlobalDataBase::GetDoubleParaFromDB("statisticalTimePeriod");
        RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

        RDouble c1 = 1.0 / (nStatisticalStep * 1.0);
        if(statisticalTimePeriod > 0.0)
        {
            c1 = MAX(c1, physicalTimeStep / statisticalTimePeriod);
        }
        RDouble c2 = 1.0 - c1;

        RDouble **qAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("qAverage"));
        RDouble **tauAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("tauAverage"));
        RDouble **q2Average = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q2Average"));

        int statisticMethod = parameters->GetStatisticMethod();
        if(statisticMethod == 0)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                RDouble u = q[IU][iCell];
                RDouble v = q[IV][iCell];
                RDouble w = q[IW][iCell];

                q2Average[0][iCell] = c2 * q2Average[0][iCell] + c1 * u * u;
                q2Average[1][iCell] = c2 * q2Average[1][iCell] + c1 * v * v;
                q2Average[2][iCell] = c2 * q2Average[2][iCell] + c1 * w * w;
                q2Average[3][iCell] = c2 * q2Average[3][iCell] + c1 * u * v;
                q2Average[4][iCell] = c2 * q2Average[4][iCell] + c1 * u * w;
                q2Average[5][iCell] = c2 * q2Average[5][iCell] + c1 * v * w;
            }

            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                RDouble uavg = qAverage[IU][iCell];
                RDouble vavg = qAverage[IV][iCell];
                RDouble wavg = qAverage[IW][iCell];

                tauAverage[0][iCell] = q2Average[0][iCell] - uavg * uavg;
                tauAverage[1][iCell] = q2Average[1][iCell] - vavg * vavg;
                tauAverage[2][iCell] = q2Average[2][iCell] - wavg * wavg;
                tauAverage[3][iCell] = q2Average[3][iCell] - uavg * vavg;
                tauAverage[4][iCell] = q2Average[4][iCell] - uavg * wavg;
                tauAverage[5][iCell] = q2Average[5][iCell] - vavg * wavg;
            }
        }

        else
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                RDouble uprm = q[IU][iCell] - qAverage[IU][iCell];
                RDouble vprm = q[IV][iCell] - qAverage[IV][iCell];
                RDouble wprm = q[IW][iCell] - qAverage[IW][iCell];

                tauAverage[0][iCell] = c2 * tauAverage[0][iCell] + c1 * uprm * uprm;
                tauAverage[1][iCell] = c2 * tauAverage[1][iCell] + c1 * vprm * vprm;
                tauAverage[2][iCell] = c2 * tauAverage[2][iCell] + c1 * wprm * wprm;
                tauAverage[3][iCell] = c2 * tauAverage[3][iCell] + c1 * uprm * vprm;
                tauAverage[4][iCell] = c2 * tauAverage[4][iCell] + c1 * uprm * wprm;
                tauAverage[5][iCell] = c2 * tauAverage[5][iCell] + c1 * vprm * wprm;
            }
        }
    }
}

void NSSolverUnstruct::StoreHistoricalResidual(Grid *gridIn)
{
    ;
}

void NSSolverUnstruct::DualTimeSource(Grid *gridIn)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady)
    {
        return;
    }

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    RDouble **q     = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qn1   = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_unsteady_n1"));
    RDouble **qn2   = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_unsteady_n2"));

    RDouble **res   = reinterpret_cast<RDouble **> (grid->GetDataPtr("res"));
    RDouble **resn1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_unsteady_n1"));
    RDouble **resn2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_unsteady_n2"));
    RDouble **resTmp = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_unsteady_tmp"));

    RDouble *vol   = grid->GetCellVolume();    //! The volume at current timestep. If it is dualtime step method, it represents the volume at timestep of n+1.
    RDouble *voln1 = grid->GetCellVolume();    //! The volume at timestep of n  .
    RDouble *voln2 = grid->GetCellVolume();    //! the volume at timestep of n-1.

    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble refGama = parameters->GetRefGama();

    //! the primitive variables
    RDouble *prim0 = new RDouble[nEquation];
    RDouble *prim1 = new RDouble[nEquation];
    RDouble *prim2 = new RDouble[nEquation];

    //! the conservative variables
    RDouble *qcsv0 = new RDouble[nEquation];
    RDouble *qcsv1 = new RDouble[nEquation];
    RDouble *qcsv2 = new RDouble[nEquation];

    RDouble Ttr[3] = {0.0}, Tv[3] = {0.0}, Te[3] = {0.0};

    using namespace GAS_SPACE;

    //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                resTmp[m][iCell] = res [m][iCell];
            }
        }

    //! Computation of dualtime coefficients, including three coefficients for Residual of R(p), R(n) and R(n-1),
    //! and three coefficients for conservative variables of qcsv(p), qcsv(n), and qcsv(n-1),
    //! the fourth coefficient is also used for the real    time part of the diagonal data for the matrix of the left term when LUSGS is running,
    //! the last   coefficient is      used for the virtual time part of the diagonal data for the matrix of the left term when LUSGS is running.
    RDouble dualTimeCoefficient[7];
    const int methodOfDualTime = GlobalDataBase::GetIntParaFromDB("methodOfDualTime");
    ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);

    RDouble dualTimeResC1 = dualTimeCoefficient[0];
    RDouble dualTimeResC2 = dualTimeCoefficient[1];
    RDouble dualTimeResC3 = dualTimeCoefficient[2];
    RDouble dualTimeQC1   = dualTimeCoefficient[3];
    RDouble dualTimeQC2   = dualTimeCoefficient[4];
    RDouble dualTimeQC3   = dualTimeCoefficient[5];

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            prim0[m] = q  [m][iCell];
            prim1[m] = qn1[m][iCell];
            prim2[m] = qn2[m][iCell];
        }
        if (nTemperatureModel > 1)
        {
            gas->GetTemperature(prim0, Ttr[0], Tv[0], Te[0]);
            gas->GetTemperature(prim1, Ttr[1], Tv[1], Te[1]);
            gas->GetTemperature(prim2, Ttr[2], Tv[2], Te[2]);
        }
        //! In fact, there is no need to input the variable of refGama for the function of Primitive2Conservative, it is to be improved.
        gas->Primitive2Conservative(prim0, refGama, Tv[0], Te[0], qcsv0);
        gas->Primitive2Conservative(prim1, refGama, Tv[1], Te[1], qcsv1);
        gas->Primitive2Conservative(prim2, refGama, Tv[2], Te[2], qcsv2);

        for (int m = 0; m < nEquation; ++ m)
        {
            RDouble dualSrcRes = dualTimeResC1 * res  [ m ][ iCell ] +
                                dualTimeResC2 * resn1[ m ][ iCell ] +
                                dualTimeResC3 * resn2[ m ][ iCell ];

            RDouble dualSrcQ   = dualTimeQC1 * qcsv0[ m ] * vol  [ iCell ] + 
                                dualTimeQC2 * qcsv1[ m ] * voln1[ iCell ] + 
                                dualTimeQC3 * qcsv2[ m ] * voln2[ iCell ];

            res[m][iCell] = dualSrcRes + dualSrcQ;
        }
    }

    delete [] prim0;    prim0 = nullptr;
    delete [] prim1;    prim1 = nullptr;
    delete [] prim2;    prim2 = nullptr;
    delete [] qcsv0;    qcsv0 = nullptr;
    delete [] qcsv1;    qcsv1 = nullptr;
    delete [] qcsv2;    qcsv2 = nullptr;
}

void NSSolverUnstruct::ChemicalSource(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    RDouble *vol  = grid->GetCellVolume();

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));
    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nm = parameters->GetNSEquationNumber();
    int nLaminar  = parameters->GetLaminarNumber();

    RDouble **src = NewPointer2<RDouble>(nLaminar, nTotalCell);

    using namespace GAS_SPACE;
    using namespace IDX;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        gas->ChemicalSource(q, t, vol, src, iCell);
        for (int m = nm; m < nLaminar; ++ m)
        {
            res[m][iCell] += src[m - nm][iCell];
        }
    }

    DelPointer2(src);
}

void NSSolverUnstruct::gravitySource(Grid* gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int isGravity = 0;
    if (GlobalDataBase::IsExist("isGravity", PHINT, 1))
    {
        isGravity = GlobalDataBase::GetIntParaFromDB("isGravity");
    }
    if (isGravity == 0)
    {
        return;
    }

    int nTotalCell = grid->GetNTotalCell();
    RDouble *vol  = grid->GetCellVolume();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refVelocity = parameters->GetRefMachNumber() * parameters->GetRefDimensionalSonicSpeed();

    RDouble gx = 0.0;
    RDouble gy = -9.81;
    RDouble gz = 0.0;

    RDouble gxNonDim = gx / refVelocity / refVelocity;
    RDouble gyNonDim = gy / refVelocity / refVelocity;
    RDouble gzNonDim = gz / refVelocity / refVelocity;

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        res[IRU][iCell] += (q[IR][iCell] - 1.0 ) * vol[iCell] * gxNonDim;
        res[IRV][iCell] += (q[IR][iCell] - 1.0 ) * vol[iCell] * gyNonDim;
        res[IRW][iCell] += (q[IR][iCell] - 1.0 ) * vol[iCell] * gzNonDim;
        res[IRE][iCell] += (q[IR][iCell] - 1.0 ) * vol[iCell] * 
            (gxNonDim * q[IRU][iCell] + gyNonDim * q[IRV][iCell] + gzNonDim * q[IRW][iCell]);
    }
}

void NSSolverUnstruct::PorousMediumSource(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *vol = grid->GetCellVolume();

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int isPorousZone = parameters->GetIsPorousZone();

    if (isPorousZone == NON_POROUS)
    {
        return;
    }

    SimpleVC *volumeCondition = gridIn->GetVolumeConditionIn();
    int vcType = volumeCondition->GetVCType();
    if (vcType != 0)
    {
        return;
    }

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **res = reinterpret_cast<RDouble **> (grid->GetDataPtr("res"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));

    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble referenceLength = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");
    RDouble referenceLengthSquare = referenceLength * referenceLength;
    const RDouble *viscousResistanceCoeff = parameters->GetViscousResistanceCoeff();
    const RDouble *inertialResistanceCoeff = parameters->GetInertialResistanceCoeff();

    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        res[IRU][iCell] -= viscousResistanceCoeff[0] * referenceLengthSquare * viscousLaminar[iCell] * q[IU][iCell] * vol[iCell] / refReNumber;
        res[IRV][iCell] -= viscousResistanceCoeff[1] * referenceLengthSquare * viscousLaminar[iCell] * q[IV][iCell] * vol[iCell] / refReNumber;
        res[IRW][iCell] -= viscousResistanceCoeff[2] * referenceLengthSquare * viscousLaminar[iCell] * q[IW][iCell] * vol[iCell] / refReNumber;

    }

    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        res[IRU][iCell] -= half * inertialResistanceCoeff[0] * referenceLength * q[IR][iCell] * fabs(q[IU][iCell]) * q[IU][iCell] * vol[iCell];
        res[IRV][iCell] -= half * inertialResistanceCoeff[1] * referenceLength * q[IR][iCell] * fabs(q[IV][iCell]) * q[IV][iCell] * vol[iCell];
        res[IRW][iCell] -= half * inertialResistanceCoeff[2] * referenceLength * q[IR][iCell] * fabs(q[IW][iCell]) * q[IW][iCell] * vol[iCell];
    }

}

void NSSolverUnstruct::RotatingSource(Grid* gridIn)
{
    UnstructGrid *grid = UnstructGridCast (gridIn);
    int nTotalCell = grid->GetNTotalCell();

    RDouble *vol = grid->GetCellVolume();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refVelocity = parameters->GetRefMachNumber() * parameters->GetRefDimensionalSonicSpeed();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");

    //! angular velocity
    RDouble Omega[100];
    GlobalDataBase::GetData("Omega", &Omega, PHDOUBLE, nTurboZone);

    //! Parallel
    int iTurboZone = grid->GetOrdinaryGridIndex();
    //! Serial
    if (iTurboZone == -1)
    {
        iTurboZone = grid->GetZoneID();
    }

    Omega[iTurboZone] /= refVelocity;

    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        res[IRV][iCell] += q[IR][iCell] * Omega[iTurboZone] * q[IW][iCell] * vol[iCell];
        res[IRW][iCell] -= q[IR][iCell] * Omega[iTurboZone] * q[IV][iCell] * vol[iCell];
    }
}

void NSSolverUnstruct::ComputeFirstLayerGridHeight(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    using namespace IDX;

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            //! iFace is the face number in the set of faceIndex.
            if (!IsWall(bcType) && bcType != PHENGLEI::ABLATION_SURFACE)
            {
                continue;
            }

#ifdef USE_ALTERNATIVE_CODE
            RDouble cellVolume = gridCellVolume(iWall, jWall, kWall); //The volume of first cell.
                                                                      //RDouble deltaHeight = half * cellVolume / surfaceArea;
#else
            //! The projection of the dy called the distance between the center of the cell and that of the surface on normal vector.
            //RDouble deltaHeight = fabs((xcc - xfc) *xfn[iFace] + (ycc - yfc) * yfn[iFace] + (zcc - zfc) * zfn[iFace]);
#endif
        }
    }
}

void NSSolverUnstruct::ComputeGamaAndTemperature(Grid *gridIn)
{
#ifdef USE_CUDA
    Param_NSSolverUnstruct *parameters_gpu = GetControlParameters();
    CallGPUComputeGamaAndTemperature(gridIn, parameters_gpu);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    RDouble **t   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));
    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));

    RDouble refGama = parameters->GetRefGama();

    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int mTT = parameters->GetmTT();
    int mTV = parameters->GetmTV();
    int mTE = parameters->GetmTE();

    using namespace GAS_SPACE;
    using namespace IDX;

    if (nChemical == 1)
    {
        RDouble *prim = new RDouble[nEquation];
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                prim[m] = q[m][iCell];
            }
            gas->GetSpecificHeatRatioAndTemperatute(prim, gama[iCell], t[mTT][iCell], t[mTV][iCell], t[mTE][iCell]);
            if (nTemperatureModel > 1)
            {
                t[ITV][iCell] = t[mTE][iCell];
                if (nTemperatureModel == 3)
                {
                    t[ITE][iCell] = t[mTE][iCell];
                }
            }
        }
        delete [] prim;    prim = nullptr;
    }
    else
    {
        RDouble coefficientofStateEquation = gas->GetCoefficientOfStateEquation();

        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            RDouble &rm = q[IR][iCell];
            RDouble &pm = q[IP][iCell];
            gama[iCell] = refGama;
            t[ITT][iCell] = pm / (coefficientofStateEquation * rm);
        }
    }
}

void NSSolverUnstruct::ComputeNodeValue(Grid *gridIn)
{
#ifdef USE_CUDA
    int nEquation_gpu = GetNumberOfEquations();
    CallGPUComputeNodeValue(gridIn, nEquation_gpu);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalNode = grid->GetNTotalNode();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble **q = reinterpret_cast<RDouble **>(grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast<RDouble **>(grid->GetDataPtr("t"));

    RDouble     **qNode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qnode"));
    RDouble     **tNode = reinterpret_cast<RDouble **>(grid->GetDataPtr("tnode"));
    RDouble *nodeWeight = reinterpret_cast<RDouble *>(grid->GetDataPtr("nodeWeight"));

    int nEquation = GetNumberOfEquations();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int *nodeNumberOfEachFace = grid->GetNodeNumberOfEachFace();
    int **face2nodeArray      = grid->GetFace2NodeArray();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        nodeWeight[iNode] = 0.0;
        for(int m =0; m < nEquation; ++ m)
        {
            qNode[m][iNode] = 0; 
        }
        tNode[0][iNode] = 0;
    }

    //! Node BC type.
    int *nodeBC = NewPointer< int > (nTotalNode);
    PHSPACE::SetField(nodeBC, PHENGLEI::NO_BOUNDARY_CONDITION, nTotalNode);

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {            
                    int point = face2nodeArray[iFace][jNode];
                    nodeBC[point] = bcType;
                }
            }
        }
    }

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType == PHENGLEI::FARFIELD)
        {
            vector<int>* faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;

                //! Far field has the second highest priority.
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2nodeArray[iFace][jNode];
                    if(nodeBC[point] == PHENGLEI::SOLID_SURFACE)
                    {
                        continue;
                    }
                    nodeBC[point] = bcType;
                }
            }
        }
    }

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            if(bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::FARFIELD ||
                bcType == PHENGLEI::SYMMETRY || bcType == PHENGLEI::INTERFACE || bcType == PHENGLEI::OVERSET)
            {
                continue;
            }

            //! other BC, except symmetry and interface.
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                int point = face2nodeArray[iFace][jNode];
                if(nodeBC[point] == PHENGLEI::SOLID_SURFACE || nodeBC[point] == PHENGLEI::FARFIELD)
                {
                    continue;
                }
                nodeBC[point] = bcType;
            }
        }
    }

    //! Boundary faces.
    RDouble **qOnBCFace = NewPointer2< RDouble >(nEquation + 1, nBoundFace);    //! '1' is temperature.
    PHSPACE::SetField(qOnBCFace, nEquation + 1, nBoundFace, zero);

    //! Init q on BC faces.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            if (bcType == PHENGLEI::INTERFACE || bcType == PHENGLEI::SYMMETRY || bcType == PHENGLEI::OVERSET)
            {
                continue;
            }

            int le = leftCellOfFace [iFace];
            int re = rightCellOfFace[iFace];

            for (int m = 0; m < nEquation; ++ m)
            {
                qOnBCFace[m][iFace] = half * (q[m][le] + q[m][re]);
            }

            qOnBCFace[nEquation][iFace] = half * (t[0][le] + t[0][re]);
        }
    }

    //! Step1: solid surface.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2nodeArray[iFace][jNode];

                    RDouble dx = x[point] - xfc[iFace];
                    RDouble dy = y[point] - yfc[iFace];
                    RDouble dz = z[point] - zfc[iFace];
                    RDouble dist = DISTANCE(dx, dy, dz);
                    RDouble weightTemp = 1.0 / dist;

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        RDouble faceQtemp = qOnBCFace[m][iFace];
                        qNode[m][point] += faceQtemp * weightTemp;
                    }
                    tNode[0][point] += qOnBCFace[nEquation][iFace] * weightTemp;

                    nodeWeight[point] += weightTemp;
                }
            }
        }
    }

    //! Step2: Far-field.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType == PHENGLEI::FARFIELD)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2nodeArray[iFace][jNode];
                    if(nodeBC[point] == PHENGLEI::SOLID_SURFACE)
                    {
                        continue;
                    }

                    RDouble dx = x[point] - xfc[iFace];
                    RDouble dy = y[point] - yfc[iFace];
                    RDouble dz = z[point] - zfc[iFace];
                    RDouble dist = DISTANCE(dx, dy, dz);
                    RDouble weightTemp = 1.0 / dist;

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        RDouble faceQtemp = qOnBCFace[m][iFace];
                        qNode[m][point] += faceQtemp * weightTemp;
                    }
                    tNode[0][point] += qOnBCFace[nEquation][iFace] * weightTemp;

                    nodeWeight[point] += weightTemp;
                }
            }
        }
    }

    //! Step3: other BC except solid surface/far field/symmetry/Interface.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType != PHENGLEI::SOLID_SURFACE && bcType != PHENGLEI::FARFIELD && 
            bcType != PHENGLEI::SYMMETRY      && bcType != PHENGLEI::INTERFACE && bcType != PHENGLEI::OVERSET)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2nodeArray[iFace][jNode];
                    if(nodeBC[point] == PHENGLEI::SOLID_SURFACE || nodeBC[point] == PHENGLEI::FARFIELD)
                    {
                        continue;
                    }

                    RDouble dx = x[point] - xfc[iFace];
                    RDouble dy = y[point] - yfc[iFace];
                    RDouble dz = z[point] - zfc[iFace];
                    RDouble dist = DISTANCE(dx, dy, dz);
                    RDouble weightTemp = 1.0 / dist;

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        RDouble faceQtemp = qOnBCFace[m][iFace];
                        qNode[m][point] += faceQtemp * weightTemp;
                    }
                    tNode[0][point] += qOnBCFace[nEquation][iFace] * weightTemp;

                    nodeWeight[point] += weightTemp;
                }
            }
        }
    }

    //! Step4: Now the interior points.
    //!        Importantly, the symmetry points are used as interior points.
    int **cell2NodeArray      = grid->GetCell2NodeArray();
    int *nodeNumberOfEachCell = grid->GetNodeNumberOfEachCell();
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int nNode = nodeNumberOfEachCell[iCell];
        for (int jNode = 0; jNode < nNode; ++ jNode)
        {
            int point = cell2NodeArray[iCell][jNode];
            if(nodeBC[point] != PHENGLEI::NO_BOUNDARY_CONDITION)
            {
                continue;
            }

            //! NO_BC && Symmetry && Interface?
            RDouble dx = x[point] - xcc[iCell];
            RDouble dy = y[point] - ycc[iCell];
            RDouble dz = z[point] - zcc[iCell];
            RDouble dist = DISTANCE(dx, dy, dz);
            RDouble weightTemp = 1.0 / dist;

            for (int m = 0; m < nEquation; ++ m)
            {
                RDouble cellQtemp = q[m][iCell];
                qNode[m][point] += cellQtemp * weightTemp;
            }
            tNode[0][point] += t[0][iCell] * weightTemp;

            nodeWeight[point] += weightTemp; 
        }
    }

    //! The velocity modify has not been considered on the symmetry!
    //! Change the number of points divided by the interpoint to the number of points before partition,because it is a face loop,each point is calculated twice,so the number of points needs to be multiplied by 2.
    InterpointInformation *interPointInfor = grid->GetInterpointInfo();
    if (interPointInfor)
    {
        int numberOfInterpoints = interPointInfor->GetNumberOfInterpoints();
        int *interPoint2GlobalPoint = interPointInfor->GetInterPoint2GlobalPoint();
        int *labelOfInterPoint = interPointInfor->GetLabelOfInterPoint();

        for (int iPoint = 0; iPoint < numberOfInterpoints; ++ iPoint)
        {
            int globalPoint = interPoint2GlobalPoint[iPoint];
            if (labelOfInterPoint[iPoint] != 1)
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    qNode[m][globalPoint] = 0;
                }
                tNode[0][globalPoint] = 0;
            }
        }

        RDouble **qInterPoint = reinterpret_cast<RDouble **>(grid->GetDataPtr("qInterPoint"));
        RDouble **tInterPoint = reinterpret_cast<RDouble **>(grid->GetDataPtr("tInterPoint"));
        for (int iPoint = 0; iPoint < numberOfInterpoints; ++ iPoint)
        {
            for (int m = 0; m < nEquation; m++)
            {
                qInterPoint[m][iPoint] =0;
            }
            tInterPoint[0][iPoint] =0;
        }
    }

    DelPointer2(qOnBCFace);
    DelPointer(nodeBC);
}


void NSSolverUnstruct::ModifyNodeValue(Grid *grid_in)
{
#ifdef USE_CUDA
    int nEquation_gpu = GetNumberOfEquations();
    CallGPUModifyNodeValue(grid_in, nEquation_gpu);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(grid_in);
    int nTotalNode = grid->GetNTotalNode();

    int nEquation = GetNumberOfEquations();

    RDouble     **qnode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qnode"));
    RDouble     **tnode = reinterpret_cast<RDouble **>(grid->GetDataPtr("tnode"));
    RDouble *nodeWeight = reinterpret_cast<RDouble  *>(grid->GetDataPtr("nodeWeight"));

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (nodeWeight[iNode] > SMALL)
        {
            for (int m = 0; m < nEquation; m++)
            {
                qnode[m][iNode] /= nodeWeight[iNode];
            }
            tnode[0][iNode] /= nodeWeight[iNode];
        }
        else
        {
            for (int m = 0; m < nEquation; m++)
            {
                qnode[m][iNode] = 0.0;
            }
            tnode[0][iNode] = 0.0;
        }
    }

    string gradientName = GlobalDataBase::GetStrParaFromDB("gradientName");
    if (gradientName == "ggnodelaplacian")
    {
        CompNodeVarForGGNodeLaplacian(grid);
    }
}

void NSSolverUnstruct::CompNodeVarForGGNodeLaplacian(Grid *grid_in)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);
    int nTotalNode = grid->GetNTotalNode();

    int nEquation = GetNumberOfEquations();

    RDouble **qnode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qnode"));
    RDouble **tnode = reinterpret_cast<RDouble **>(grid->GetDataPtr("tnode"));

    RDouble **q = reinterpret_cast<RDouble **>(grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast<RDouble **>(grid->GetDataPtr("t"));

    //! Get necessary geometry information.
    int *nCPN = grid->GetCellNumberOfEachNode();
    int *n2c  = grid->GetNode2Cell();
    int inode, j, icell, count;
    int *knode;

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();
    RDouble *xn  = grid->GetX();
    RDouble *yn  = grid->GetY();
    RDouble *zn  = grid->GetZ();

    RDouble weight, weight1, weight2, *n_count;
    RDouble *lamdax = grid->GetLamdax();
    RDouble *lamday = grid->GetLamday();
    RDouble *lamdaz = grid->GetLamdaz();

    if (lamdax == nullptr || lamday == nullptr || lamdaz == nullptr)
    {
        grid->CalcLaplacianWeitht();
        lamdax = grid->GetLamdax();
        lamday = grid->GetLamday();
        lamdaz = grid->GetLamdaz();
    }

    knode = grid->Getknode();

    n_count = new RDouble [ nTotalNode ];
    for( inode = 0; inode < nTotalNode; ++ inode )
    {
        for (int m = 0; m < nEquation; m++)
        {
            qnode[m][inode] = 0.0;
        }
        tnode[0][inode] = 0.0;
        n_count[ inode ] = 0.0;
    }
    count = 0;
    for( inode = 0; inode < nTotalNode; ++ inode )
    {
        for( j = 0; j < nCPN[ inode ]; ++j )
        {
            icell = n2c[ count++ ];
            weight1 = 1 + lamdax[ inode ] * ( xcc[ icell ] - xn[ inode ] )
                + lamday[ inode ] * ( ycc[ icell ] - yn[ inode ] )
                + lamdaz[ inode ] * ( zcc[ icell ] - zn[ inode ] );
            weight2 = 1 / sqrt( ( xcc[ icell ] - xn[ inode ] ) * ( xcc[ icell ] - xn[ inode ] )
                + ( ycc[ icell ] - yn[ inode ] ) * ( ycc[ icell ] - yn[ inode ] )
                + ( zcc[ icell ] - zn[ inode ] ) * ( zcc[ icell ] - zn[ inode ] ) );
            weight  = weight1 * ( 1 - knode[ inode ] ) + weight2 * knode[ inode ];

            //q_n[ inode ] += weight * q[ icell ];

            for (int m = 0; m < nEquation; m++)
            {
                qnode[m][ inode ] += weight * q[m][ icell ];
            }
            tnode[0][ inode ] += weight * t[0][ icell ];
            n_count[ inode ] += weight;
        }
    }

    for( inode = 0; inode < nTotalNode; ++ inode ) 
    {
        for (int m = 0; m < nEquation; m++)
        {
            qnode[m][inode] /= n_count[ inode ] + TINY;
        }
        tnode[0][inode] /= n_count[ inode ] + TINY;
    }

    delete [] n_count;    n_count = nullptr;
}

RDouble* NSSolverUnstruct::CompSoundSpeed(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble * > (grid->GetDataPtr("gama"));
    RDouble *c = new RDouble[nTotal];

    using namespace IDX;

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        c[iCell] = sqrt(gama[iCell] * q[IP][iCell] / q[IR][iCell]);
    }

    return c;
}

RDouble** NSSolverUnstruct::ComputePrimitiveVariablesWithMoleFraction(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nm = parameters->GetNSEquationNumber();
    int nSpeciesNumber = parameters->GetNumberOfSpecies();

    RDouble **moleFraction = NewPointer2<RDouble>(nSpeciesNumber, nTotal);
    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        RDouble *massFractionTemp = new RDouble[nSpeciesNumber];
        RDouble *moleFractionTemp = new RDouble[nSpeciesNumber];
        for (int m = 0; m < nSpeciesNumber; ++ m)
        {
            massFractionTemp[m] = q[nm + m][iCell];
        }
        gas->ComputeMoleFractionByMassFraction(massFractionTemp,moleFractionTemp);

        for (int m = 0; m < nSpeciesNumber; ++ m)
        {
            moleFraction[m][iCell] = moleFractionTemp[m];
        }
        delete []massFractionTemp;    massFractionTemp = nullptr;
        delete []moleFractionTemp;    moleFractionTemp = nullptr;
    }
    return moleFraction;
}

RDouble** NSSolverUnstruct::ComputeDimensionalVariables(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refDensity = parameters->GetRefDimensionalDensity();
    RDouble refTemperature = parameters->GetRefDimensionalTemperature();
    RDouble refVelocity = parameters->GetRefMachNumber() * parameters->GetRefDimensionalSonicSpeed();
    RDouble refPressure = refDensity * refVelocity * refVelocity;

    RDouble velocityMagnitude;
    RDouble **dimensionalVariables   = NewPointer2<RDouble>(7, nTotal);

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        //! Compute the dimensional density.
        dimensionalVariables[0][iCell] = refDensity * q[IR][iCell];

        //! Compute the dimensional velocity.
        dimensionalVariables[1][iCell] = refVelocity * q[IU][iCell];
        dimensionalVariables[2][iCell] = refVelocity * q[IV][iCell];
        dimensionalVariables[3][iCell] = refVelocity * q[IW][iCell];
        velocityMagnitude = sqrt(q[IU][iCell] * q[IU][iCell] + q[IV][iCell] * q[IV][iCell] + q[IW][iCell] * q[IW][iCell]);
        dimensionalVariables[4][iCell] = refVelocity * velocityMagnitude;
        //! Compute the dimensional pressure.
        dimensionalVariables[5][iCell] = refPressure * q[IP][iCell];

        //! Compute the dimensional temperature.
        dimensionalVariables[6][iCell] = refTemperature * t[ITT][iCell];
    }

    return dimensionalVariables;
}

RDouble* NSSolverUnstruct::CompMachNumber(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble * > (grid->GetDataPtr("gama"));
    RDouble *mach = new RDouble[nTotal];

    using namespace IDX;

    RDouble rm,um,vm,wm,pm,c2,v2;
    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        rm = q[IR][iCell];
        um = q[IU][iCell];
        vm = q[IV][iCell];
        wm = q[IW][iCell];
        pm = q[IP][iCell];
        c2 = abs(gama[iCell] * pm / rm);
        v2 = um * um + vm * vm + wm * wm;

        mach[iCell] = sqrt(v2 / c2);
    }
    return mach;
}

RDouble * NSSolverUnstruct::ComputeVorticitybyQCriteria(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;
    RDouble * vorticity = new RDouble [ nTotal ];

    if(!WantVisualField(grid))
    {
        //! Even if do not visual field vorticity, 
        //! it would be dumped on boundary, to compatible with the field variable.
        PHSPACE::SetField(vorticity, -1.0e10, nTotal);
        return vorticity;
    }

    RDouble **dqdx = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarZ"));

    using namespace IDX;

    RDouble * dudxField = dqdx[ IU ];
    RDouble * dudyField = dqdy[ IU ];
    RDouble * dudzField = dqdz[ IU ];

    RDouble * dvdxField = dqdx[ IV ];
    RDouble * dvdyField = dqdy[ IV ];
    RDouble * dvdzField = dqdz[ IV ];

    RDouble * dwdxField = dqdx[ IW ];
    RDouble * dwdyField = dqdy[ IW ];
    RDouble * dwdzField = dqdz[ IW ];

    for (int iCell = 0; iCell < nTotal; ++ iCell) 
    {
        RDouble dudx = dudxField[ iCell ];
        RDouble dudy = dudyField[ iCell ];
        RDouble dudz = dudzField[ iCell ];

        RDouble dvdx = dvdxField[ iCell ];
        RDouble dvdy = dvdyField[ iCell ];
        RDouble dvdz = dvdzField[ iCell ];

        RDouble dwdx = dwdxField[ iCell ];
        RDouble dwdy = dwdyField[ iCell ];
        RDouble dwdz = dwdzField[ iCell ];

        RDouble q1 = dvdy * dwdz - dwdy * dvdz;
        RDouble q2 = dudx * dwdz - dwdx * dudz;
        RDouble q3 = dudx * dvdy - dvdx * dudy;

        vorticity[ iCell ] = q1 + q2 + q3; 
    }

    return vorticity;
}

void NSSolverUnstruct::ComputeVorticitybyQCriteria(Grid *gridIn, RDouble * vorticity_x, RDouble * vorticity_y, RDouble * vorticity_z, RDouble * vorticityMagnitude,  RDouble * strain_rate,  RDouble * Q_criteria)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **dqdx = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarZ"));
    using namespace IDX;

    RDouble * dudxField = dqdx[ IU ];
    RDouble * dudyField = dqdy[ IU ];
    RDouble * dudzField = dqdz[ IU ];

    RDouble * dvdxField = dqdx[ IV ];
    RDouble * dvdyField = dqdy[ IV ];
    RDouble * dvdzField = dqdz[ IV ];

    RDouble * dwdxField = dqdx[ IW ];
    RDouble * dwdyField = dqdy[ IW ];
    RDouble * dwdzField = dqdz[ IW ];

    RDouble Omx, Omy, Omz, Om;
    RDouble s11, s12, s13, s22, s23, s33, S;

    for (int iCell = 0; iCell < nTotal; ++ iCell) 
    {
        RDouble dudx = dudxField[ iCell ];
        RDouble dudy = dudyField[ iCell ];
        RDouble dudz = dudzField[ iCell ];

        RDouble dvdx = dvdxField[ iCell ];
        RDouble dvdy = dvdyField[ iCell ];
        RDouble dvdz = dvdzField[ iCell ];

        RDouble dwdx = dwdxField[ iCell ];
        RDouble dwdy = dwdyField[ iCell ];
        RDouble dwdz = dwdzField[ iCell ];


        s11 = dudx; s22 = dvdy; s33 = dwdz;
        s12 = half * (dudy + dvdx); s13 = half * (dudz + dwdx); s23 = half * (dvdz + dwdy);
        S = two * (s11 * s11 + s22 * s22 + s33 * s33 + two*(s12 * s12 + s13 * s13 + s23 * s23));
        S = sqrt (S);
        strain_rate[iCell] = S;

        Omx = dwdy - dvdz; 
        Omy = dudz - dwdx;
        Omz = dvdx - dudy;
        Om = DISTANCE(Omx, Omy, Omz);

        vorticity_x[iCell] = Omx; 
        vorticity_y[iCell] = Omy;
        vorticity_z[iCell] = Omz;
        vorticityMagnitude[iCell] = Om; 

        Q_criteria[iCell] = fourth * (Om * Om - S * S);
    }
}

RDouble * NSSolverUnstruct::ComputeCp(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();

    using namespace IDX;
    RDouble poo = primitiveVarFarfield[IP];

    //! Pressure cofficients: cp = [p - 1.0 / (gama * mach * mach)] / 2.0 = (p - poo) / 2.0.
    RDouble *cp = new RDouble [nTotal];
    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        cp[iCell] = 2.0 * (q[IP][iCell] - poo);
    }

    return cp;
}

RDouble * NSSolverUnstruct::ComputeCpAverage(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;
    RDouble **qAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("qAverage"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();

    using namespace IDX;
    RDouble poo = primitiveVarFarfield[IP];

    //! Pressure cofficients: cp = [p - 1.0 / (gama * mach * mach)] / 2.0 = (p - poo) / 2.0.
    RDouble * cp = new RDouble [nTotal];
    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        cp[iCell] = 2.0 * (qAverage[IP][iCell] - poo);
    }

    return cp;
}

void NSSolverUnstruct::GetQlQr(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    int nEquation = GetNumberOfEquations();

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();

    //! GMRESPassQC
    RDouble **qLC = faceProxy->GetQLC();
    RDouble **qRC = faceProxy->GetQRC();

    //! GMRES
    int* leftcellindexofFace = faceProxy->GetLeftCellIndexOfFace();
    int* rightcellindexofFace = faceProxy->GetRightCellIndexOfFace();

    //! iFace is the global ID, jFace is the local ID in this section.
    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];
        int jFace  = iFace - localStart;

        //! GMRES
        leftcellindexofFace[jFace] = le;
        rightcellindexofFace[jFace] = re;

        for (int m = 0; m < nEquation; ++ m)
        {
            qL[m][jFace] = q[m][le];
            qR[m][jFace] = q[m][re];

            //! GMRESPassQC
            qLC[m][jFace] = q[m][le];
            qRC[m][jFace] = q[m][re];
        }
    }
}

void NSSolverUnstruct::GetGamaLR(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble *gama = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));

    RDouble *gamal = faceProxy->GetGamaL();
    RDouble *gamar = faceProxy->GetGamaR();

    //! iFace is the global ID, jFace is the local ID in this section.
    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];
        int jFace  = iFace - localStart;

        gamal[jFace] = gama[le];
        gamar[jFace] = gama[re];
    }
}

void NSSolverUnstruct::GetPressureFactorLR(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble *rtem = reinterpret_cast <RDouble *> (grid->GetDataPtr("rtem"));

    RDouble *pressureCoeffL = faceProxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR = faceProxy->GetPressureCoefficientR();

    for (int iFace = localStart; iFace < localEnd; ++iFace)
    {
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];
        int jFace = iFace - localStart;

        pressureCoeffL[jFace] = rtem[le];
        pressureCoeffR[jFace] = rtem[re];
    }
}

void NSSolverUnstruct::GetTimeCoefficientLR(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble *timeCoefficientInverse = reinterpret_cast<RDouble *> (grid->GetDataPtr("timeCoefficientInverse"));

    RDouble *timeCoeffl = faceProxy->GetTimeCoefficientL();
    RDouble *timeCoeffr = faceProxy->GetTimeCoefficientR();

    //! iFace is the global ID, jFace is the local ID in this section.
    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];
        int jFace  = iFace - localStart;

        timeCoeffl[jFace] = timeCoefficientInverse[le];
        timeCoeffr[jFace] = timeCoefficientInverse[re];
    }
}

void NSSolverUnstruct::GetPreconCoefficientLR(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble *preconCoefficient = reinterpret_cast<RDouble *> (grid->GetDataPtr("preconCoefficient"));;

    RDouble *preconCoeffl = faceProxy->GetPreconCoefficientL();
    RDouble *preconCoeffr = faceProxy->GetPreconCoefficientR();

    //! iFace is the global ID, jFace is the local ID in this section.
    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];
        int jFace  = iFace - localStart;

        preconCoeffl[jFace] = preconCoefficient[le];
        preconCoeffr[jFace] = preconCoefficient[re];
    }
}

void NSSolverUnstruct::ZeroResiduals(Grid *gridIn)
{
#ifdef USE_CUDA
    int nEquation_gpu = GetNumberOfEquations();
    CallGPUZeroResiduals(gridIn, nEquation_gpu);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    int nEquation = GetNumberOfEquations();

    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    for (int m = 0; m < nEquation; ++ m)
    {
        PHSPACE::SetField(res[m], 0.0, nTotal);
    }

#ifdef USE_GMRESSOLVER
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");

    if( tscheme == GMRES )
    {
        //! GMRES zero residual
        RDouble **dRdq = reinterpret_cast<RDouble**> (grid->GetDataPtr("dRdq"));
        //! PHSPACE::SetField(dRdq, nEquation * nTotal, nEquation * nTotal, 0.0);
        //! GMRESCSR
        vector<int> AI = grid->GetJacobianAI4GMRES();
        PHSPACE::SetField(dRdq, nEquation, nEquation * AI[nTotal], 0.0);

        //! GMRESJac1st
        int jacOrder = grid->GetJacobianOrder();
        if(jacOrder == 2)
        {
            RDouble **dRdq1st = reinterpret_cast<RDouble**> (grid->GetDataPtr("dRdq1st"));
            vector<int> AI1st = grid->GetJacobianAI1st4GMRES();
            PHSPACE::SetField(dRdq1st, nEquation, nEquation * AI1st[nTotal], 0.0);
        }

        //! GMRESCoupled
        int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
        if( viscousType == ONE_EQU )
        {
            RDouble **dRdqCoupledTerm = reinterpret_cast <RDouble**> (grid->GetDataPtr("dRdqCoupledTerm"));
            PHSPACE::SetField(dRdqCoupledTerm, nEquation, AI[nTotal], 0.0);
        }
    }
#endif

}

void NSSolverUnstruct::ComputeMinTimeStep(Grid *gridIn, RDouble &minDt, RDouble &maxDt)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *dt = reinterpret_cast< RDouble *  > (grid->GetDataPtr("dt"));

    minDt = LARGE;
    maxDt = - LARGE;

#ifdef USE_CUDA
    CallGPUComputeMinTimeStep(nTotalCell);
    return;
#endif

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        minDt = MIN(minDt, dt[iCell]);
        maxDt = MAX(maxDt, dt[iCell]);
    }
}

void NSSolverUnstruct::ReduceMaxTimeStep(Grid *gridIn, RDouble globalMinDt)
{
#ifdef USE_CUDA
    CallGPUReduceMaxTimeStep(gridIn);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *volume = grid->GetCellVolume();
    RDouble *dt = reinterpret_cast< RDouble *  > (grid->GetDataPtr("dt"));

    RDouble ktmax = GlobalDataBase::GetDoubleParaFromDB("ktmax");

    if (ktmax > 0)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            dt[iCell] = MIN(dt[iCell], ktmax * globalMinDt);
        }
    }

    //! The time step is divided by the volume, so it is need not divided other else.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dt[iCell] /= volume[iCell];
    }

    grid->SetGhostCell(dt);
}

//! Inviscid spectral radius = 0.5 * (|V.n| + c) * S.
//! Referenced in R. F. Chen && Z. J. Wang, AIAA Journal 2000, 38(12).
void NSSolverUnstruct::SpectrumRadiusInviscid(Grid *gridIn)
{
#ifdef USE_CUDA
    CallGPUSpectrumRadiusInviscid(gridIn);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble *vgn  = grid->GetFaceNormalVelocity();

    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();

    RDouble **q   = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble *invSpectrum = reinterpret_cast<RDouble *> (grid->GetDataPtr("invSpectralRadius"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

    RDouble *preconCoefficient = nullptr;
    RDouble *timeCoefficientInverse = nullptr;

    if (ifLowSpeedPrecon != 0)
    {
        preconCoefficient = reinterpret_cast<RDouble *> (grid->GetDataPtr("preconCoefficient"));
        if (isUnsteady)
        {
            timeCoefficientInverse = reinterpret_cast<RDouble *> (grid->GetDataPtr("timeCoefficientInverse"));
        }
    }

    RDouble *rho, *u, *v, *w, *p;

    rho = q[IDX::IR];
    u   = q[IDX::IU];
    v   = q[IDX::IV];
    w   = q[IDX::IW];
    p   = q[IDX::IP];

    PHSPACE::SetField(invSpectrum, 0.0, nTotalCell);

    for(int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];

        RDouble rl  = rho[le];
        RDouble ul  = u[le];
        RDouble vl  = v[le];
        RDouble wl  = w[le];
        RDouble pl  = p[le];
        //RDouble vnl = xfn[iFace] * ul + yfn[iFace] * vl + zfn[iFace] * wl;
        RDouble vnl = xfn[iFace] * ul + yfn[iFace] * vl + zfn[iFace] * wl - vgn[iFace];

        RDouble rr  = rho[re];
        RDouble ur  = u[re];
        RDouble vr  = v[re];
        RDouble wr  = w[re];
        RDouble pr  = p[re];
        //RDouble vnr = xfn[iFace] * ur + yfn[iFace] * vr + zfn[iFace] * wr;
        RDouble vnr = xfn[iFace] * ur + yfn[iFace] * vr + zfn[iFace] * wr  - vgn[iFace];

        RDouble specificHeatRatio = half * (gama[le] + gama[re]);
        RDouble pm = half * (pl + pr);
        RDouble rm = half * (rl + rr);
        RDouble vn = half * (vnl + vnr);
        RDouble cm = sqrt(specificHeatRatio * pm / rm);

        if(ifLowSpeedPrecon != 0)
        {
            RDouble c2 = specificHeatRatio * pm / rm;
            RDouble preconCoeff = half * (preconCoefficient[le] + preconCoefficient[re]);

            cm = half * sqrt(((1.0 - preconCoeff) * vn) * ((1.0 - preconCoeff) * vn) + 4.0 * preconCoeff * c2);
            vn = half * vn * (1.0 + preconCoeff);

            if (isUnsteady)
            {
            RDouble timeCoeff = half * (timeCoefficientInverse[le] + timeCoefficientInverse[re]);
                cm *= timeCoeff;
                vn *= timeCoeff;
            }
        }

        RDouble inviscidSpectrumRadius = half * area[iFace] * (ABS(vn) + cm);

        invSpectrum[le] += inviscidSpectrumRadius;
        if(re < nTotalCell)
        {
            //! The re cell is the interior cell.
            invSpectrum[re] += inviscidSpectrumRadius;
        }
    }
}

void NSSolverUnstruct::SpectrumRadiusViscous(Grid *gridIn)
{
#ifdef USE_CUDA
    Param_NSSolverUnstruct *parameters_gpu = GetControlParameters();
    CallGPUSpectrumRadiusViscous(gridIn, parameters_gpu);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble *volume = grid->GetCellVolume();

    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();

    RDouble **q   = reinterpret_cast< RDouble **> (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble * > (grid->GetDataPtr("gama"));
    RDouble *visSpectrum = reinterpret_cast< RDouble *> (grid->GetDataPtr("visSpectralRadius"));

    RDouble *rho, *u, *v, *w, *p;
    rho = q[IDX::IR];
    u   = q[IDX::IU];
    v   = q[IDX::IV];
    w   = q[IDX::IW];
    p   = q[IDX::IP];

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *viscousLaminar = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast< RDouble * > (grid->GetDataPtr("vist"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();

    RDouble prandtlLaminar    = parameters->GetPrandtlLaminar();    //! default 0.72
    RDouble prandtlTurbulence = parameters->GetPrandtlTurbulence(); //! default 0.9

    const RDouble fourthThird = 4.0 / 3.0;
    RDouble oprl = 1.0 / prandtlLaminar;
    RDouble oprt = 1.0 / prandtlTurbulence;

    PHSPACE::SetField(visSpectrum, 0.0, nTotalCell);

    int visSpetrumMethod = 2;
    if(visSpetrumMethod == 1)
    {
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            int le = leftCellofFace [iFace];
            int re = rightCellofFace[iFace];

            RDouble visLaminar = half * (viscousLaminar[le] + viscousLaminar[re]);
            RDouble visTurb    = half * (viscousTurbulence[le] + viscousTurbulence[re]);
            RDouble density    = half * (rho[le] + rho[re]);

            RDouble specificHeatRatio = half * (gama[le] + gama[re]);
            RDouble coefOfSpectrum1 = fourthThird * (visLaminar + visTurb);
            RDouble coefOfSpectrum2 = specificHeatRatio * (visLaminar * oprl + visTurb * oprt);

            //! Viscous spectral radius type 2.
            //! Referenced in Wang Gang, Chinese Journal of Aeronautics, 2012 (25): 33-41, eq. 13.
            RDouble coefOfSpectrum  = two * PHSPACE::MAX(coefOfSpectrum1, coefOfSpectrum2) / (refReNumber * density);

            RDouble faceArea = area[iFace];
            RDouble faceAera2 = PHSPACE::SQR(faceArea);

            RDouble viscousSpectrumRadius = faceAera2 * coefOfSpectrum;

            visSpectrum[le] += viscousSpectrumRadius / volume[le];
            if(re < nTotalCell)
            {
                visSpectrum[re] += viscousSpectrumRadius / volume[re];
            }
        }
    }
    else if (visSpetrumMethod == 2)
    {
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            int le = leftCellofFace [iFace];
            int re = rightCellofFace[iFace];

            RDouble dx = xcc[re] - xcc[le];
            RDouble dy = ycc[re] - ycc[le];
            RDouble dz = zcc[re] - zcc[le];
            RDouble ds = ABS(xfn[iFace] * dx + yfn[iFace] * dy + zfn[iFace] * dz);

            RDouble visLaminar = half * (viscousLaminar[le] + viscousLaminar[re]);
            RDouble visTurb    = half * (viscousTurbulence[le] + viscousTurbulence[re]);
            RDouble density    = half * (rho[le] + rho[re]);
            RDouble viscosity  = visLaminar + visTurb;
            RDouble faceArea   = area[iFace];

            //! Viscous spectral radius = 2 * miu/(rho.|n.dr|) * S.
            //! Referenced in R. F. Chen && Z. J. Wang, AIAA Journal 2000, 38(12).
            RDouble viscousSpectrumRadius = 2.0 * viscosity / (density * ds * refReNumber + SMALL);
            viscousSpectrumRadius *= half * faceArea;

            visSpectrum[le] += viscousSpectrumRadius;
            if(re < nTotalCell)
            {
                visSpectrum[re] += viscousSpectrumRadius;
            }
        }
    }
}

void NSSolverUnstruct::SpectrumRadiusChemical(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    RDouble *vol = grid->GetCellVolume();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));
    RDouble **chemSpectrum = reinterpret_cast< RDouble ** > (grid->GetDataPtr("chemSpectralRadius"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nSpeciesNumber = parameters->GetNumberOfSpecies();

    PHSPACE::SetField(chemSpectrum, nSpeciesNumber, nTotalCell, zero);

    RDouble **srs = NewPointer2<RDouble>(nSpeciesNumber, nTotalCell);

    using namespace IDX;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        gas->ChemicalSpectrumRadius(q, t, vol, srs, iCell);

        for (int m = 0; m < nSpeciesNumber - 1; ++ m)
        {
            chemSpectrum[m][iCell] += srs[m][iCell];
        }
    }

    DelPointer2(srs);
}

void NSSolverUnstruct::LocalTimeStep(Grid *gridIn)
{
#ifdef USE_CUDA
    Param_NSSolverUnstruct *parameters_gpu = GetControlParameters();
    Param_NSSolver *parameters_nsgpu = GetControlParameters();
    CallGPULocalTimeStep(gridIn, parameters_gpu, parameters_nsgpu);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);
    RDouble *volume    = grid->GetCellVolume();

    int nTotalCell = grid->GetNTotalCell();

    RDouble *dt = reinterpret_cast< RDouble *  > (grid->GetDataPtr("dt"));
    RDouble *invSpectralRadius = reinterpret_cast<RDouble *> (grid->GetDataPtr("invSpectralRadius"));
    RDouble *visSpectralRadius = reinterpret_cast<RDouble *> (grid->GetDataPtr("visSpectralRadius"));

    RDouble *CFLCell = new RDouble [nTotalCell];

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    RDouble *spectralRadius = new RDouble [nTotalCell];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        spectralRadius[iCell] = invSpectralRadius[iCell];
    }

    bool isViscous = parameters->IsViscous();

    if (isViscous)
    {
        //! Consider the viscous time step.
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            spectralRadius[iCell] += visSpectralRadius[iCell];
        }    
    }

    LimitCFL(grid,CFLCell);

    //! Set dt to zero.
    PHSPACE::SetField(dt, 0.0, nTotalCell);

    //! time step.dt=CFL*volume/(invSpectralRadius + visSpectralRadius)
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dt[iCell] = CFLCell[iCell] * volume[iCell] / spectralRadius[iCell];
    }

    DelPointer(spectralRadius);
    DelPointer(CFLCell);
}

void NSSolverUnstruct::LimitCFL(Grid *gridIn, RDouble *CFLCell)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int level = grid->GetLevel();

    FieldProxy *qProxy = GetFieldProxy(gridIn, "q");
    RDouble **q = qProxy->GetField_UNS();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble CFL = ComputeCFL();

    int isInMGIniting = GlobalDataBase::GetIntParaFromDB("isInMGIniting");
    if (!grid->IsFinestGrid() && !isInMGIniting)
    {
        RDouble mgCFLScale = parameters->GetMgCFLScale();
        CFL *= mgCFLScale;
    }

    PHSPACE::SetField(CFLCell, CFL, nTotalCell);

    int nChemical = parameters->GetChemicalFlag();
    if(nChemical != 0)
    {
        return;
    }

    RDouble CFLStart = GlobalDataBase::GetDoubleParaFromDB("CFLStart");
    RDouble pTotal = GlobalDataBase::GetDoubleParaFromDB("pTotal");

    const RDouble LOWER = 0.001;
    const RDouble UPPER = 0.02;
    RDouble pLowerLimit = LOWER * pTotal;
    RDouble pUpperLimit = UPPER * pTotal;

    RDouble *cellSkewness = grid->GetCellSkewness();

    //! adjust CFLCell according to p in the cell.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if(q[IP][iCell] < pLowerLimit)
        {
            CFLCell[iCell] = CFLStart;
        }
        else if(q[IP][iCell] < pUpperLimit)
        {
            CFLCell[iCell] = (q[IP][iCell]-pUpperLimit)*(CFL-CFLStart)/(pUpperLimit-pLowerLimit) + CFL;
        }
    }

    if(0 == level)
    {
        //! adjust CFLCell according to skewness in the cell.
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            if(cellSkewness[iCell] < tenth)
            {
                CFLCell[iCell] *= tenth;
            }
            else if(cellSkewness[iCell] < ten)
            {
                CFLCell[iCell] *= (one-exp(zero-cellSkewness[iCell]));
            }
        }
    }
}

void NSSolverUnstruct::GlobalTimeStep(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble *vol = grid->GetCellVolume();
    RDouble *dt  = reinterpret_cast< RDouble * > (grid->GetDataPtr("dt" ));

    int nTotalCell = grid->GetNTotalCell();

    RDouble dtau;
    GlobalDataBase::GetData("dtau", &dtau, PHDOUBLE, 1);

    //! Set dt to zero.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dt[iCell] = dtau / vol[iCell];
    }

    grid->SetGhostCell(dt);
    RDouble dtmin = dtau;
    GlobalDataBase::UpdateData("dtmin", &dtmin, PHDOUBLE, 1);
}

void NSSolverUnstruct::LocalGlobalTimeStep(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble *vgn  = grid->GetFaceNormalVelocity();
    RDouble *vol  = grid->GetCellVolume();

    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));
    RDouble *dt   = reinterpret_cast< RDouble *  > (grid->GetDataPtr("dt"));

    RDouble cfl = ComputeCFL();

    //RDouble refGama = parameters->GetRefGama();

    int le, re;
    RDouble um, vm, wm, cm;

    RDouble *rho, *u, *v, *w, *p;

    using namespace IDX;

    rho = q[IR];
    u   = q[IU];
    v   = q[IV];
    w   = q[IW];
    p   = q[IP];

    //! Set dt to zero.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dt[iCell] = 0.0;
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellofFace[ iFace ];

        um      = u[le];
        vm      = v[le];
        wm      = w[le];
        cm      = sqrt(gama[le] * p[le] / rho[le]);
        dt[le] += (ABS(xfn[iFace] * um + yfn[iFace] * vm + zfn[iFace] * wm - vgn[iFace]) + cm) * area[iFace];
    }

    //! For interior faces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellofFace[ iFace ];
        re = rightCellofFace[ iFace ];

        um      = u[le];
        vm      = v[le];
        wm      = w[le];
        cm      = sqrt(gama[le] * p[le] / rho[le]);
        dt[le] += (ABS(xfn[iFace] * um + yfn[iFace] * vm + zfn[iFace] * wm - vgn[iFace]) + cm) * area[iFace];

        um      = u[re];
        vm      = v[re];
        wm      = w[re];
        cm      = sqrt(gama[re] * p[re] / rho[re]);
        dt[re] += (ABS(xfn[iFace] * um + yfn[iFace] * vm + zfn[iFace] * wm - vgn[iFace]) + cm) * area[iFace];
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dt[iCell] = cfl * vol[iCell] / dt[iCell];
    }

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    if (viscousType > INVISCID)
    {
        RDouble *dtv;
        dtv = new RDouble [ nTotalCell ];

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            dtv[iCell] = 0.0;
        }

        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            le = leftCellofFace[ iFace ];
            re = rightCellofFace[ iFace ];

            dtv[le] += area[iFace] * area[iFace];
        }

        //! For interior faces.
        for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
        {
            le = leftCellofFace[ iFace ];
            re = rightCellofFace[ iFace ];

            dtv[le] += area[iFace] * area[iFace];
            dtv[re] += area[iFace] * area[iFace];
        }

        RDouble *visl = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));
        RDouble *vist = reinterpret_cast< RDouble * > (grid->GetDataPtr("vist"));

        RDouble refReNumber  = parameters->GetRefReNumber();
        RDouble prandtlLaminar    = parameters->GetPrandtlLaminar();
        RDouble prandtlTurbulence = parameters->GetPrandtlTurbulence();

        RDouble coef1,coef2,coef,vis_l,vis_t;
        RDouble csrvist = 1.0;

        const RDouble foth = 4.0 / 3.0;

        RDouble oprl = 1.0 / prandtlLaminar;
        RDouble oprt = 1.0 / prandtlTurbulence;

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            vis_l = visl[iCell];
            vis_t = csrvist * vist[iCell];
            coef1 = foth * (vis_l + vis_t);
            coef2 = gama[iCell] * (vis_l * oprl + vis_t * oprt);
            coef  = MAX(coef1, coef2);
            dtv[iCell] = cfl * rho[iCell] * refReNumber * vol[iCell] * vol[iCell] / (coef * dtv[iCell]);
        }

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            dt[iCell] = dt[iCell] * dtv[iCell] / (dt[iCell] + dtv[iCell]);
        }

        delete [] dtv;    dtv = nullptr;
    }

    RDouble dtmin = LARGE;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dtmin = MIN(dtmin, static_cast<RDouble>(dt[iCell]));
    }

    GlobalDataBase::UpdateData("dtmin", &dtmin, PHDOUBLE, 1);

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dt[iCell] = dtmin / vol[iCell];
    }
}

void NSSolverUnstruct::FillField(Grid *gridIn, FieldProxy *field_proxy, RDouble value)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    int nEquation = GetNumberOfEquations();

    RDouble **field = field_proxy->GetField_UNS();

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            field[m][iCell] = value;
        }
    }
}

void NSSolverUnstruct::FillField(Grid *gridIn, FieldProxy *fieldTarget, FieldProxy *fieldSource)
{
#ifdef USE_CUDA
    int nEquation_gpu = GetNumberOfEquations();
    CallGPUFillField(gridIn, nEquation_gpu);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    int nEquation = GetNumberOfEquations();

    RDouble **fieldTar = fieldTarget->GetField_UNS();
    RDouble **fieldSrc = fieldSource->GetField_UNS();

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            fieldTar[m][iCell] = fieldSrc[m][iCell];
        }
    }
}

void NSSolverUnstruct::RungeKuttaResidual(Grid *gridIn, FieldProxy *dqProxy, RDouble coef)
{
#ifdef USE_CUDA
    int nEquation_gpu = GetNumberOfEquations();
    CallGPURungeKuttaResidual(gridIn, nEquation_gpu, coef);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int      nTotalCell = grid->GetNTotalCell();

    RDouble **dq = dqProxy->GetField_UNS();
    RDouble **res = reinterpret_cast <RDouble **> (grid->GetDataPtr("res"));

    int nEquation = GetNumberOfEquations();

    RDouble * dt = reinterpret_cast <RDouble *> (grid->GetDataPtr("dt"));

    for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            //! The volume has been divided in time step, so the dt is actually 'dt/vol'.
            dq[iEquation][iCell] = dt[iCell] * coef * res[iEquation][iCell];
        }
    }
}

//! Restrict from fine to coarse grid for all q.
void NSSolverUnstruct::RestrictAllQ(Grid *fgrid, Grid *cgrid)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nEquation = GetNumberOfEquations();

    RDouble **fq = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("q"));
    RDouble **cq = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("q"));

    for (int m = 0; m < nEquation; ++ m)
    {
        RestrictQ(fgrid, fq[m], cgrid, cq[m]);
    }

    RDouble *fgama = reinterpret_cast< RDouble * > (fgrid->GetDataPtr("gama"));
    RDouble *cgama = reinterpret_cast< RDouble * > (cgrid->GetDataPtr("gama"));

    RestrictQ(fgrid, fgama, cgrid, cgama);
    
    int viscousType = parameters->GetViscousType();
    if (viscousType > INVISCID)
    {
        RDouble * fvisl = reinterpret_cast< RDouble * > (fgrid->GetDataPtr("visl"));
        RDouble * cvisl = reinterpret_cast< RDouble * > (cgrid->GetDataPtr("visl"));
        RestrictQ(fgrid, fvisl, cgrid, cvisl);

        if (viscousType > LAMINAR)
        {
            RDouble * fvist = reinterpret_cast< RDouble * > (fgrid->GetDataPtr("vist"));
            RDouble * cvist = reinterpret_cast< RDouble * > (cgrid->GetDataPtr("vist"));
            RestrictQ(fgrid, fvist, cgrid, cvist);
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble **fqn1 = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("q_unsteady_n1"));
        RDouble **fqn2 = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("q_unsteady_n2"));
        RDouble **cqn1 = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("q_unsteady_n1"));
        RDouble **cqn2 = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("q_unsteady_n2"));

        RDouble **fres   = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("res"           ));
        RDouble **fresn1 = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("res_unsteady_n1"));
        RDouble **fresn2 = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("res_unsteady_n2"));
        RDouble **cres   = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("res"           ));
        RDouble **cresn1 = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("res_unsteady_n1"));
        RDouble **cresn2 = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("res_unsteady_n2"));

        for (int m = 0; m < nEquation; ++ m)
        {
            RestrictQ(fgrid, fqn1[m], cgrid, cqn1[m]);
            RestrictQ(fgrid, fqn2[m], cgrid, cqn2[m]);

            RestrictQ(fgrid, fres[m],   cgrid, cres[m]);
            RestrictQ(fgrid, fresn1[m], cgrid, cresn1[m]);
            RestrictQ(fgrid, fresn2[m], cgrid, cresn2[m]);
        }
    }
}

void NSSolverUnstruct::LoadFlux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
#ifdef USE_CUDA
    CallGPULoadFlux(gridIn, localStart, localEnd);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();

    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();

    RDouble **flux = faceProxy->GetFlux();

    //! Determine if there are boundary faces.
    int nMid  = localStart;
    if (localEnd <= nBoundFace)
    {
       //! If all boundary faces.
       nMid = localEnd;
    }
    else if (localStart < nBoundFace)
    {
       //! Part of them are boundary faces.
       nMid = nBoundFace;
    }

    //! 'rhs' := inviscid flux - viscid flux.
    //! 'res' is defined as '-rhs'.

    RDouble **res = reinterpret_cast<RDouble **> (grid->GetDataPtr("res"));

    //! For boundary faces, remember re is ghost cell.
    for (int iFace = localStart; iFace < nMid; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        int jFace = iFace - localStart;
        for (int m = 0; m < nLaminar; ++ m)
        {
            res[m][le] -= flux[m][jFace];
        }
    }

    //! Interior faces.
    for (int iFace = nMid; iFace < localEnd; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];
        int jFace = iFace - localStart;

        for (int m = 0; m < nLaminar; ++ m)
        {
            res[m][le] -= flux[m][jFace];
            res[m][re] += flux[m][jFace];
        }
    }
}


void NSSolverUnstruct::InitialInviscidfaceproxyOnlyForMixGrid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalFace = grid->GetNTotalFace();
    int nEquation = 5;

    GeomProxy *geomProxy = new GeomProxy();
    geomProxy->Create(nTotalFace);
    FillGeomProxy(grid, geomProxy, 0, nTotalFace);

    InviscidfaceProxyOnlyForMixGrid = new FaceProxy();
    InviscidfaceProxyOnlyForMixGrid->Create(nTotalFace, nEquation);
    InviscidfaceProxyOnlyForMixGrid->SetGeomProxy(geomProxy);

}

void NSSolverUnstruct::InitialViscousfaceproxyOnlyForMixGrid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalFace = grid->GetNTotalFace();
    int nEquation = 5;

    ViscousfaceProxyOnlyForMixGrid = new FaceProxy();
    ViscousfaceProxyOnlyForMixGrid->Create(nTotalFace, nEquation);

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int numberOfSpecies = parameters->GetNumberOfSpecies();
    
    NSFaceValue *facevar = new NSFaceValue(nEquation, numberOfSpecies, nTotalFace);
    ViscousfaceProxyOnlyForMixGrid->SetNSFaceValue(facevar);

    ComputeFaceWeight(grid, ViscousfaceProxyOnlyForMixGrid, 0, nTotalFace);
}


void NSSolverUnstruct::UpdateQlQrOnlyforMixGrid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalFace = grid->GetNTotalFace();

    int localStart = 0;
    int localEnd = nTotalFace;

    GetInviscidFaceValue(grid, InviscidfaceProxyOnlyForMixGrid, limiterOnlyForMixGrid, localStart, localEnd);
}

void NSSolverUnstruct::UpdateInviscidfluxOnlyforMixGrid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalFace = grid->GetNTotalFace();

    int localStart = 0;
    int localEnd = nTotalFace;

    ComputeInviscidFlux(grid, InviscidfaceProxyOnlyForMixGrid, localStart, localEnd);
}

void NSSolverUnstruct::UpdateViscousfluxOnlyforMixGrid(Grid *gridIn)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();
    if (viscousType == INVISCID) return;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalFace = grid->GetNTotalFace();

    int localStart = 0;
    int localEnd = nTotalFace;

    GetVisFaceValue(grid, ViscousfaceProxyOnlyForMixGrid, localStart, localEnd);

    ComputeVisflux(grid, ViscousfaceProxyOnlyForMixGrid, localStart, localEnd);
}

void NSSolverUnstruct::LoadGlobalfacefluxtoResOnlyforMixGrid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();

    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();

    RDouble **res = reinterpret_cast<RDouble **> (grid->GetDataPtr("res"));
    RDouble **Inviscidflux = InviscidfaceProxyOnlyForMixGrid->GetFlux();
    //! For boundary faces, remember re is ghost cell.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        for (int m = 0; m < nLaminar; ++ m)
        {
            res[m][le] -= Inviscidflux[m][iFace];
        }
    }

    //! Interior faces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

        for (int m = 0; m < nLaminar; ++ m)
        {
            res[m][le] -= Inviscidflux[m][iFace];
            res[m][re] += Inviscidflux[m][iFace];
        }
    }

    int viscousType = parameters->GetViscousType();
    if (viscousType == INVISCID) return;

    RDouble **Viscousflux = ViscousfaceProxyOnlyForMixGrid->GetFlux();
    //! For boundary faces, remember re is ghost cell.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        for (int m = 0; m < nLaminar; ++ m)
        {
            res[m][le] -= Viscousflux[m][iFace];
        }
    }

    //! Interior faces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

        for (int m = 0; m < nLaminar; ++ m)
        {
            res[m][le] -= Viscousflux[m][iFace];
            res[m][re] += Viscousflux[m][iFace];
        }
    }
}

void NSSolverUnstruct::ComputeViscousCoeff(Grid *gridIn)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nChemical = parameters->GetChemicalFlag();
    if (nChemical == 0)
    {
        ComputeViscousCoefficientWithoutChemical(gridIn);
    }
    else
    {
        ComputeViscousCoefficient(gridIn);
    }
}

void NSSolverUnstruct::ComputeViscousCoefficient(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble wallTemperature = parameters->GetWallTemperature();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));

    using namespace IDX;
    RDouble *primitiveVars  = new RDouble[nEquation];
    RDouble temperature, electronTemperature, viscosity;

    RDouble refTemperatureSutherland = GlobalDataBase::GetDoubleParaFromDB("tsuth");

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        temperature = t[ITT][iCell];
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            primitiveVars[iEquation] = q[iEquation][iCell];
        }
        electronTemperature = temperature;
        if (nTemperatureModel == 2)
        {
            electronTemperature = t[ITV][iCell];
        }
        else if (nTemperatureModel == 3)
        {
            electronTemperature = t[ITE][iCell];
        }
        gas->ComputeViscosityByPrimitive(primitiveVars, temperature, electronTemperature, refTemperatureSutherland, viscosity);
        viscousLaminar[iCell] = viscosity;
    }
    if (wallTemperature <= 0.0)
    {
        grid->SetGhostCell(viscousLaminar);
    }
    else
    {
        for (int iCell = nTotalCell; iCell < nTotal; ++ iCell)
        {
            temperature = t[ITT][iCell];
            for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
            {
                primitiveVars[iEquation] = q[iEquation][iCell];
            }
            electronTemperature = temperature;
            if (nTemperatureModel == 2)
            {
                electronTemperature = t[ITV][iCell];
            }
            else if (nTemperatureModel == 3)
            {
                electronTemperature = t[ITE][iCell];
            }
            gas->ComputeViscosityByPrimitive(primitiveVars, temperature, electronTemperature, refTemperatureSutherland, viscosity);
            viscousLaminar[iCell] = viscosity;
        }
    }

    delete [] primitiveVars;    primitiveVars = nullptr;
}

void NSSolverUnstruct::ComputeViscousCoefficientWithoutChemical(Grid *gridIn)
{
#ifdef USE_CUDA
    Param_NSSolverUnstruct *parameters_gpu = GetControlParameters();
    CallGPUComputeViscousCoefficientWithoutChemical(gridIn, parameters_gpu);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble **t = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();

    RDouble viscousLaminarMin;
    GlobalDataBase::GetData("visl_min", &viscousLaminarMin, PHDOUBLE, 1);

    using namespace IDX;

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        viscousLaminar[iCell] = t[ITT][iCell] * sqrt(t[ITT][iCell]) * (1.0 + nonDimensionalSutherlandTemperature) / (t[ITT][iCell] + nonDimensionalSutherlandTemperature);
    }
}

void NSSolverUnstruct::InitCGrid(Grid *fgrid_in, Grid *cgrid_in)
{
    UnstructGrid *fgrid = UnstructGridCast(fgrid_in);
    UnstructGrid *cgrid = UnstructGridCast(cgrid_in);
    int *cell2coarsegridcell = fgrid->GetCell2CoarseGridCell();
    int fnTotalCell = fgrid->GetNTotalCell();
    int cnTCell = cgrid->GetNTotalCell();
    RDouble *cvol = cgrid->GetCellVolume();
    RDouble *fvol = fgrid->GetCellVolume();

    RDouble **fq = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("q"));
    RDouble **cq = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("q"));

    int nEquation = GetNumberOfEquations();

    int cn = cnTCell + cgrid->GetNBoundFace();

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < cn; ++ iCell)
        {
            cq[m][iCell] = 0.0;
        }
    }

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < fnTotalCell; ++ iCell)
        {
            cq[m][ cell2coarsegridcell[iCell] ] += fvol[iCell] * fq[m][iCell];
        }
    }
    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < cnTCell; ++ iCell)
        {
            cq[m][iCell] /= cvol[iCell];
        }
    }
}

void NSSolverUnstruct::ComputePreconditionCoefficient(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gamma = reinterpret_cast< RDouble * > (grid->GetDataPtr("gama"));
    RDouble *preconCoefficient = reinterpret_cast< RDouble * > (grid->GetDataPtr("preconCoefficient"));

    RDouble *dt = reinterpret_cast< RDouble * > (grid->GetDataPtr("dt"));
    RDouble *timeCoefficient = reinterpret_cast< RDouble * > (grid->GetDataPtr("timeCoefficient"));
    RDouble *timeCoefficientInverse = reinterpret_cast< RDouble * > (grid->GetDataPtr("timeCoefficientInverse"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble kPreconCoeff = parameters->GetPreconCoefficient();

    int isUnsteady = parameters->GetIsUnsteady();
    RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");
    RDouble machinf2 = refMachNumber * refMachNumber;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble rm  = q[IR][iCell];
        RDouble um  = q[IU][iCell];
        RDouble vm  = q[IV][iCell];
        RDouble wm  = q[IW][iCell];
        RDouble pm  = q[IP][iCell];
        RDouble gama = gamma[iCell];

        RDouble c2 = gama * pm / rm;
        RDouble v2 = um * um + vm * vm + wm * wm;
        RDouble mach2 = v2 / c2;
        RDouble machPrec2 = MIN(MAX(mach2,kPreconCoeff * machinf2),1.0);

        preconCoefficient[iCell] = machPrec2;
    }

    grid->SetGhostCell(preconCoefficient);
}

void NSSolverUnstruct::ComputePreconditionCoefficientUnsteady(Grid *gridIn)
{
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

    RDouble *vol = grid->GetCellVolume();
    RDouble *cellVelocityX = grid->GetCellVelocityX();
    RDouble *cellVelocityY = grid->GetCellVelocityY();
    RDouble *cellVelocityZ = grid->GetCellVelocityZ();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gamma = reinterpret_cast< RDouble * > (grid->GetDataPtr("gama"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble kPreconCoeff = parameters->GetPreconCoefficient();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble machinf2Pre = kPreconCoeff * refMachNumber * refMachNumber;

    RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");
    RDouble *dt = reinterpret_cast< RDouble * > (grid->GetDataPtr("dt"));
    RDouble *timeCoefficient = reinterpret_cast< RDouble * > (grid->GetDataPtr("timeCoefficient"));
    RDouble *preconCoefficient = reinterpret_cast< RDouble * > (grid->GetDataPtr("preconCoefficient"));
    RDouble *timeCoefficientInverse = reinterpret_cast< RDouble * > (grid->GetDataPtr("timeCoefficientInverse"));
    RDouble timeAccuracyLowPre =  0.5; 
    RDouble machUsteady2 = GlobalDataBase::GetDoubleParaFromDB("preconCoeffConst");

    RDouble maxTimeCoeff = 0.1667;
    if (GlobalDataBase::IsExist("maxTimeCoeff", PHINT, 1))
    {
        GlobalDataBase::GetData("maxTimeCoeff", &maxTimeCoeff, PHINT, 1);
    }
    else
    {
    GlobalDataBase::UpdateData("maxTimeCoeff", &maxTimeCoeff, PHINT, 1);
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble rm  = q[IR][iCell];
        RDouble um  = q[IU][iCell];
        RDouble vm  = q[IV][iCell];
        RDouble wm  = q[IW][iCell];
        RDouble pm  = q[IP][iCell];
        RDouble gama = gamma[iCell];

        RDouble ucc = cellVelocityX[iCell];
        RDouble vcc = cellVelocityY[iCell];
        RDouble wcc = cellVelocityZ[iCell];

        RDouble c2 = gama * pm / rm;
        RDouble v2 = (um - ucc) * (um - ucc) + (vm - vcc) * (vm - vcc) + (wm - wcc) * (wm - wcc);
        RDouble mach2 = v2 / c2;

        RDouble machPrec2 = MIN(MAX(MAX(mach2, machUsteady2), machinf2Pre),1.0);
        preconCoefficient[iCell] = machPrec2;
        
        RDouble dtemp = dt[iCell] * vol[iCell];

        RDouble ratioDt = dtemp / physicalTimeStep;
        if (ratioDt > maxTimeCoeff)
        {
            RDouble dt_max_lim = maxTimeCoeff * physicalTimeStep;
            dtemp = MIN(dtemp, dt_max_lim);
            dt[iCell] = dtemp / (vol[iCell] + SMALL);
        }

        timeCoefficient[iCell] = (1.0 + timeAccuracyLowPre )* dtemp / physicalTimeStep ;
        timeCoefficientInverse[iCell] = 1.0 / (1.0 + timeCoefficient[iCell]);
        preconCoefficient[iCell] = preconCoefficient[iCell] / (timeCoefficientInverse[iCell] + (1.0 - timeCoefficientInverse[iCell]) * preconCoefficient[iCell]);
    }

    grid->SetGhostCellExceptInterface(timeCoefficient);
    grid->SetGhostCellExceptInterface(timeCoefficientInverse);
    grid->SetGhostCellExceptInterface(preconCoefficient);
}

//! Drive individual functions to calculate inviscid fluxes.
void NSSolverUnstruct::InviscidFlux(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalFace = grid->GetNTotalFace();

    Limiter *limiter = 0;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int limiterType = parameters->GetLimiterType();
    int isInIniting = GlobalDataBase::GetIntParaFromDB("isInIniting");
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
    bool isFirstOrder = (limiterType == ILMT_FIRST) || isInIniting;

     //! Face proxy.
    FaceProxy *faceProxy = CreateFaceProxy(grid);
    faceProxy->SetGeomProxy(CreateGeomProxy(grid));

#ifdef USE_GMRESSOLVER
    RDouble** dRdqfromGrid = reinterpret_cast<RDouble**> (grid->GetDataPtr("dRdq"));
    faceProxy->SetJacobianMatrix(dRdqfromGrid);
    // GMRESCSR
    vector<int> AI = grid->GetJacobianAI4GMRES();
    vector<int> AJ = grid->GetJacobianAJ4GMRES();
    faceProxy->SetJacobianAI4GMRES(AI);
    faceProxy->SetJacobianAJ4GMRES(AJ);

    //! GMRESJac1st
    int jacOrder = grid->GetJacobianOrder();
    vector<int> AI1st = grid->GetJacobianAI1st4GMRES();
    vector<int> AJ1st = grid->GetJacobianAJ1st4GMRES();
    faceProxy->SetJacobianOrder(jacOrder);
    faceProxy->SetJacobianAI1st4GMRES(AI1st);
    faceProxy->SetJacobianAJ1st4GMRES(AJ1st);
    RDouble** dRdq1stfromGrid = reinterpret_cast<RDouble**> (grid->GetDataPtr("dRdq1st"));
    faceProxy->SetJacobianMatrix1st(dRdq1stfromGrid);

    //! GMRESBoundary
    RDouble **dDdPfromGrid = reinterpret_cast<RDouble**> (grid->GetDataPtr("dDdP"));
    faceProxy->SetdDdPMatrix(dDdPfromGrid);
#endif

    int localStart, localEnd;
    localStart = 0;
    do
    {
        localEnd = localStart + SEGCTION_LENGTH;
        if (localEnd > nTotalFace) 
        {
            localEnd = nTotalFace;
        }

        GetInviscidFaceValue(grid, faceProxy, limiter, localStart, localEnd);

        ComputeInviscidFlux(grid, faceProxy, localStart, localEnd);

        LoadFlux(grid, faceProxy, localStart, localEnd);

        localStart = localEnd;
    } while (localStart < nTotalFace);

    delete faceProxy;    faceProxy = nullptr;
    delete limiter;    limiter = nullptr;
}

void NSSolverUnstruct::CalPressureFactor(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();

    int **cell2face = grid->GetCell2Face();
    vector<int> * c2c = grid->GetCell2Cell();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble minPressureRatio = 1.3;
    RDouble maxPressureRatio = 2.6;
    RDouble pi = 3.141592653589;
    
    int RoeEntropyFixMethod = this->GetControlParameters()->GetRoeEntropyFixMethod();

    RDouble **q   = reinterpret_cast <RDouble **> (grid->GetDataPtr("q"));
    RDouble *rtem = reinterpret_cast <RDouble *> (grid->GetDataPtr("rtem"));

    if (RoeEntropyFixMethod == 2)
    {
        for (int iCell = 0; iCell < nTotalCell; ++iCell)
        {
            RDouble LeftPressureOfCell = zero;
            RDouble RightPressureOfCell = zero;
            rtem[iCell] = zero;
            for (int j = 0; j < faceNumberOfEachCell[iCell]; ++j)
            {
                int face = cell2face[iCell][j];
                int le = leftCellofFace[face];
                int re = rightCellofFace[face];
                
                LeftPressureOfCell += q[PHSPACE::IDX::IP][le];
                RightPressureOfCell += q[PHSPACE::IDX::IP][re];
            }
            rtem[iCell] += ABS(RightPressureOfCell - LeftPressureOfCell) / ABS(RightPressureOfCell + LeftPressureOfCell);
            rtem[iCell] /= static_cast <RDouble> (faceNumberOfEachCell[iCell]);
        }
    }

    if (RoeEntropyFixMethod == 6)
    {
        for (int iCell = 0; iCell < nTotalCell; ++iCell)
        {
            int nNeighbor = static_cast<int>(c2c[iCell].size());

            RDouble minPressure = q[IP][iCell];
            RDouble maxPressure = q[IP][iCell];

            for (int j = 0; j < nNeighbor; ++j)
            {
                int neighborCell = c2c[iCell][j];
                if (q[IP][neighborCell] < minPressure)
                {
                    minPressure = q[IP][neighborCell];
                }
                if (q[IP][neighborCell] > maxPressure)
                {
                    maxPressure = q[IP][neighborCell];
                }
            }

            RDouble pressureRatio = maxPressure / minPressure;
            RDouble presFunc = min(one, max(zero, (maxPressureRatio - pressureRatio) / (maxPressureRatio - minPressureRatio)));
            RDouble pressCosFunc = cos(presFunc*pi);
            pressCosFunc = max(-one, min(pressCosFunc, one));

            rtem[iCell] = half * (one - pressCosFunc);
        }
    }

    grid->SetGhostCellExceptInterface(rtem);
}

void NSSolverUnstruct::GetHeatTransferCoeff(Grid *gridIn)
{

}

void NSSolverUnstruct::ViscousFlux(Grid *gridIn)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();
    if (viscousType == INVISCID) return;
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
    INVScheme inviscidScheme = GetINVScheme(tscheme);
    // GKS 本身含粘性效应，跳过额外粘性通量计算
    if (inviscidScheme == Gks_Scheme) return;
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalFace = grid->GetNTotalFace();

#ifdef USE_CUDA
    using namespace GPUMemory;
    using namespace GPUControlVariables;
    SetGPUVisFluxSolverControlVars(d_visflux_solver ,"NSSolverUnstruct");
#endif

    FaceProxy *faceProxy = CreateFaceProxy(grid);
    faceProxy->SetNSFaceValue(CreateNSFaceValue(grid));
    int localStart, localEnd;
    localStart = 0;
    do
    {
        localEnd = localStart + SEGCTION_LENGTH;
        if (localEnd > nTotalFace)
        {
            localEnd = nTotalFace;
        }

        ComputeFaceWeight(grid, faceProxy, localStart, localEnd);

        GetVisFaceValue(grid, faceProxy, localStart, localEnd);

        if( tscheme != GMRES)
        {
            ComputeVisflux(grid, faceProxy, localStart, localEnd);
        }
#ifdef USE_GMRESSOLVER
        else
        {
            ComputeVisflux_GMRES(grid, faceProxy, localStart, localEnd);
        }
#endif
        LoadFlux(grid, faceProxy, localStart, localEnd);

        localStart = localEnd;
    } while (localStart < nTotalFace);

#ifdef USE_CUDA
    using namespace GPUMemory;
    using namespace GPUControlVariables;
    SetGPUVisFluxSolverControlVars(d_visflux_solver ,"Empty");
#endif

    delete faceProxy;    faceProxy = nullptr;
}

void NSSolverUnstruct::ComputeFaceWeight(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble *deltL = faceProxy->GetWeightL();
    RDouble *deltR = faceProxy->GetWeightR();

    grid->FaceWeight(deltL, deltR, localStart, localEnd);
}

void NSSolverUnstruct::GetVisFaceValue(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
#ifdef USE_CUDA
    using namespace GPUKernels;
    int nEquation_gpu = GetNumberOfEquations();
    Param_NSSolverUnstruct *parameters_gpu = GetControlParameters();
    CallGPUGetVisFaceValue(gridIn, parameters_gpu, nEquation_gpu, localStart, localEnd);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();

    NSFaceValue *faceVariable = faceProxy->GetNSFaceValue();

    RDouble *faceVelocityX = grid->GetFaceVelocityX();
    RDouble *faceVelocityY = grid->GetFaceVelocityY();
    RDouble *faceVelocityZ = grid->GetFaceVelocityZ();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refGama = 1.4;
    GlobalDataBase::GetData("refGama", &refGama, PHDOUBLE, 1);
    //RDouble refGama = parameters->GetRefGama();

    RDouble oPrandtlLaminar    = parameters->GetoPrandtlLaminar();
    RDouble oPrandtlTurbulence = parameters->GetoPrandtlTurbulence();

    RDouble referenceMachNumber = parameters->GetRefMachNumber();

    RDouble tWall = parameters->GetWallTemperature();

    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();

    int numberOfSpecies = parameters->GetNumberOfSpecies();
    int nChemical = parameters->GetChemicalFlag();

    int isPorousZone = parameters->GetIsPorousZone();
    RDouble porosity, kSolid;

    using namespace PHENGLEI;
    int vcType = FLUID;
    if (isPorousZone != NON_POROUS)
    {
        porosity = parameters->GetPorosity();
        kSolid = parameters->GetHeatConductivitySolid();
        SimpleVC *volumeCondition = gridIn->GetVolumeConditionIn();
        vcType = volumeCondition->GetVCType();
    }

    int nEquation = GetNumberOfEquations();

    RDouble **rhoDiffusion = faceVariable->GetRhoDS()->GetField();
    RDouble **hintSpecies = faceVariable->GetHintS()->GetField();
    RDouble **primitiveVariable = faceVariable->GetPrim()->GetField();
    RDouble **tMid   = faceVariable->GetT()->GetField();

    RDouble *kCp = faceVariable->GetKCP();   //! Heat conductivity coefficient, usually computed by Cp.
    RDouble *viscousLaminarFace = faceVariable->GetMUL();
    RDouble *viscousTurbulenceFace = faceVariable->GetMUT();
    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");
    RDouble referenceMachNumberSquare = referenceMachNumber * referenceMachNumber;

    using namespace IDX;

    RDouble *qL                = new RDouble[nEquation];
    RDouble *qR                = new RDouble[nEquation];
    RDouble *primface          = new RDouble[nEquation];
    RDouble *rhoDiffusionFace  = new RDouble[numberOfSpecies];
    RDouble *hintSpeciesOfFace = new RDouble[numberOfSpecies];

    int le,re,jFace;
    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        jFace = iFace - localStart;

        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        if (iFace < nBoundFace) re = iFace + nTotalCell;

        viscousLaminarFace[jFace]   = half * (viscousLaminar[le] + viscousLaminar[re]);
        viscousTurbulenceFace[jFace] = half * (viscousTurbulence[le] + viscousTurbulence[re]);
    }

    int nMid = nBoundFace;
    if (localEnd < nBoundFace) nMid = localEnd;

    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        if (iFace < nBoundFace) re = iFace + nTotalCell;

        jFace = iFace - localStart;

        for (int m = 0; m < nEquation; ++ m)
        {
            primitiveVariable[m][jFace] = half * (q[m][le] + q[m][re]);
        }

        tMid[ITT][jFace] = half * (t[ITT][le] + t[ITT][re]);

        if (nChemical == 1)
        {
            using namespace GAS_SPACE;

            RDouble t1  = half * (t[ITT][le] + t[ITT][re]);

            for (int m = 0; m < nEquation; ++ m)
            {
                qL[m] = q[m][le];
                qR[m] = q[m][re];
                primface[m] = primitiveVariable[m][jFace];
            }

            gas->ComputeDensityDiffusionAndKCP(t1, primface, viscousLaminarFace[jFace], viscousTurbulenceFace[jFace], rhoDiffusionFace, hintSpeciesOfFace, kCp[jFace]);
            for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
            {
                rhoDiffusion[ispecies][jFace] = rhoDiffusionFace[ispecies];
            }

            for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
            {
                hintSpecies[ispecies][jFace] = hintSpeciesOfFace[ispecies];
            }
        }
        else
        {
            RDouble cp = 1.0 / ((refGama - 1.0) * referenceMachNumberSquare);
            kCp[jFace] = (viscousLaminarFace[jFace] * oPrandtlLaminar + viscousTurbulenceFace[jFace] * oPrandtlTurbulence) * cp;
            if (isPorousZone != NON_POROUS && vcType == UNDEFINED)
            {
                kCp[jFace] *= porosity;
                kCp[jFace] += kSolid * refDimensionalTemperature * (1.0 - porosity);
            }
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();
    int isAle = parameters->GetIsCodeOfAleModel();

    //! If there are some boundary faces, make modifications. 
    //! Check if there are boundary faces. If no, return.
    if (localStart < nBoundFace)
    {
        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

        RDouble tw = tWall / refDimensionalTemperature;
        for (int iFace = localStart; iFace < nMid; ++ iFace)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
            int bcType = bcRegion->GetBCType();
            jFace  = iFace - localStart;

            le = leftCellOfFace[iFace];
            re = rightCellOfFace[iFace];

            if (bcType == PHENGLEI::INTERFACE || bcType == PHENGLEI::OVERSET) continue;

            for (int m = 0; m < nEquation; ++ m)
            {
                RDouble tmp = half * (q[m][le] + q[m][re]);
                primitiveVariable[m][jFace] = tmp;
            }

            if (bcType == PHENGLEI::SOLID_SURFACE)
            {
                RDouble uWall = 0.0;
                RDouble vWall = 0.0;
                RDouble wWall = 0.0;
                Data_Param *bcData  = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace])->GetBCParamDataBase();
                if(bcData)
                {
                    if (bcData->IsExist("uWall", PHDOUBLE, 1))
                    {
                        bcData->GetData("uWall", &uWall, PHDOUBLE, 1);
                        bcData->GetData("vWall", &vWall, PHDOUBLE, 1);
                        bcData->GetData("wWall", &wWall, PHDOUBLE, 1);
                    }
                }

                RDouble velocityXWall = uWall;
                RDouble velocityYWall = vWall;
                RDouble velocityZWall = wWall;
                if (isUnsteady && isAle)
                {
                    velocityXWall = faceVelocityX[jFace] + uWall;
                    velocityYWall = faceVelocityY[jFace] + vWall;
                    velocityZWall = faceVelocityZ[jFace] + wWall;
                }

                //! set wall velocity for translating frame and rotating frame.
                if (referenceFrame == TRANSLATIONAL_FRAME || referenceFrame == ROTATIONAL_FRAME)
                {
                    int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");

                    string shroud[100];
                    GlobalDataBase::GetData("shroud", &shroud, PHSTRING, nTurboZone);

                    velocityXWall = faceVelocityX[jFace];
                    velocityYWall = faceVelocityY[jFace];
                    velocityZWall = faceVelocityZ[jFace];

                    for (int iTurboZone = 0; iTurboZone < nTurboZone; iTurboZone++)
                    {
                        if (bcRegion->GetBCName() == shroud[iTurboZone])
                        {
                            velocityXWall = 0;
                            velocityYWall = 0;
                            velocityZWall = 0;
                        }
                    }
                }

                primitiveVariable[IU][jFace] = velocityXWall;
                primitiveVariable[IV][jFace] = velocityYWall;
                primitiveVariable[IW][jFace] = velocityZWall;

                if (tw > 0.0)
                {
                    tMid[0][jFace] = tw;
                    RDouble cp = 1.0 / ((refGama - 1.0) * referenceMachNumberSquare);
                    RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();
                    viscousLaminarFace[jFace] = tMid[0][jFace] * sqrt(tMid[0][jFace]) * (1.0 + nonDimensionalSutherlandTemperature) / (tMid[0][jFace] + nonDimensionalSutherlandTemperature);
                    kCp[jFace] = (viscousLaminarFace[jFace] * oPrandtlLaminar + viscousTurbulenceFace[jFace] * oPrandtlTurbulence) * cp;
                    if (isPorousZone != NON_POROUS && vcType == UNDEFINED)
                    {
                        kCp[jFace] *= porosity;
                        kCp[jFace] += kSolid * refDimensionalTemperature * (1.0 - porosity);
                    }
                }
            }
        }
    }

    delete [] qL;    qL = nullptr;
    delete [] qR;    qR = nullptr;
    delete [] primface;    primface = nullptr;
    delete [] rhoDiffusionFace;    rhoDiffusionFace = nullptr;
    delete [] hintSpeciesOfFace;    hintSpeciesOfFace = nullptr;
}

#ifdef USE_GMRESSOLVER
// GMRESVis GMRESAD
void NSSolverUnstruct::ComputeVisflux_GMRES(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();

    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    RDouble **dqdx = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarZ"));

    RDouble **dtdx = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradTemperatureX"));
    RDouble **dtdy = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradTemperatureY"));
    RDouble **dtdz = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradTemperatureZ"));

    RDouble **flux = faceProxy->GetFlux();

    RDouble *deltaL = faceProxy->GetWeightL();
    RDouble *deltaR = faceProxy->GetWeightR();

    NSFaceValue *faceVariable = faceProxy->GetNSFaceValue();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble *vol  = grid->GetCellVolume();
    
    RDouble *xcc  = grid->GetCellCenterX();
    RDouble *ycc  = grid->GetCellCenterY();
    RDouble *zcc  = grid->GetCellCenterZ();

    RDouble *xfc  = grid->GetFaceCenterX();
    RDouble *yfc  = grid->GetFaceCenterY();
    RDouble *zfc  = grid->GetFaceCenterZ();

    RDouble **rhoDiffusion = faceVariable->GetRhoDS()->GetField();
    RDouble **hintSpecies = faceVariable->GetHintS()->GetField();
    RDouble **primitiveVariableFace = faceVariable->GetPrim()->GetField();
    RDouble **tm     = faceVariable->GetT()->GetField();

    RDouble *kCp     = faceVariable->GetKCP();
    RDouble *viscousLaminarFace   = faceVariable->GetMUL();
    RDouble *viscousTurbulenceFace = faceVariable->GetMUT();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int numberOfSpecies = parameters->GetNumberOfSpecies();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int nm = parameters->GetNSEquationNumber();

    RDouble refReNumber  = parameters->GetRefReNumber();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    int nolstress = GlobalDataBase::GetIntParaFromDB("nolstress");
    int nrokplus  = GlobalDataBase::GetIntParaFromDB("nrokplus");

    RDouble skewnessAngle = parameters->GetSkewnessAngle();

    RDouble **t                 = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));
    RDouble **primitiveVariable = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *gama               = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble** dRdq              = reinterpret_cast<RDouble **> (grid->GetDataPtr("dRdq"));
    RDouble** dDdP              = reinterpret_cast<RDouble **> (grid->GetDataPtr("dDdP"));
    vector<int> AI              = grid->GetJacobianAI4GMRES();
    vector<int> AJ              = grid->GetJacobianAJ4GMRES();
    // GMRESJac1st
    RDouble** dRdq1st              = reinterpret_cast<RDouble **> (grid->GetDataPtr("dRdq1st"));
    vector<int> AI1st              = grid->GetJacobianAI1st4GMRES();
    vector<int> AJ1st              = grid->GetJacobianAJ1st4GMRES();
    RDouble **qTurb = 0;
    RDouble **aniss  = 0;
    RDouble coefficientofStateEquation = gas->GetCoefficientOfStateEquation();

    if (nrokplus > 0)
    {
        qTurb = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
        aniss  = reinterpret_cast<RDouble **> (grid->GetDataPtr("aniss"));
    }

    //! GMRESCoupled
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    RDouble **dRdqCoupledTerm = 0;
    if( viscousType == ONE_EQU )
    {
        dRdqCoupledTerm = reinterpret_cast <RDouble**> (grid->GetDataPtr("dRdqCoupledTerm"));
    }
    RDouble **primitiveturb = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"  ));
    RDouble **dDdP_turb     = reinterpret_cast<RDouble **> (grid->GetDataPtr("dDdP_turb"));
    RDouble *dFluxdnutL     = new RDouble[nEquation]();
    RDouble *dFluxdnutR     = new RDouble[nEquation]();

    //! GMRES2ndCorrection
    vector<int> BCLeftCells     = grid->GetBCLeftCells();
    vector<int> BCRightCells    = grid->GetBCRightCells();
    vector<int> BCFaces         = grid->GetBCFaces();

    int nIndependentVars = nEquation * 2 + (nEquation + nTemperatureModel) * 6;
    int halfIndependentVars = nIndependentVars / 2;
    //! GMRESCoupled
    if( viscousType == ONE_EQU )
    {
        nIndependentVars += 2;    //! two extra terms include the derivative w.r.t nut_L, nut_R
    }

    int colidx, colidx2, colidx1st, colidx21st;

    ADReal *dqdxL = new ADReal[nEquation+nTemperatureModel];
    ADReal *dqdyL = new ADReal[nEquation+nTemperatureModel];
    ADReal *dqdzL = new ADReal[nEquation+nTemperatureModel];

    ADReal *dqdxR = new ADReal[nEquation+nTemperatureModel];
    ADReal *dqdyR = new ADReal[nEquation+nTemperatureModel];
    ADReal *dqdzR = new ADReal[nEquation+nTemperatureModel];
    //RDouble *dgrad = new RDouble[nEquation+nTemperatureModel];

    RDouble *primL = new RDouble[nEquation];
    RDouble *primR = new RDouble[nEquation];
    RDouble *fN    = new RDouble[nEquation];    //! for neighbors
    ADReal *fL     = new ADReal[nEquation];
    ADReal *fR     = new ADReal[nEquation];
    ADReal *fluxp  = new ADReal[nLaminar];
    ADReal fLturb;    //! GMRESCoupled
    ADReal fRturb;    //! GMRESCoupled

    RDouble t1x, t1y, t1z, t2x, t2y, t2z;
    RDouble txx, tyy, tzz;
    RDouble txy, txz, tyz;

    bool isFineGrid = grid->IsFinestGrid();

    using namespace IDX;
    RDouble dFluxdpL[nEquation][nEquation];
    RDouble dFluxdpR[nEquation][nEquation];
    RDouble dGradtmp[nEquation][nEquation];
    RDouble dGradtmpBC[nEquation][nEquation];
    RDouble dDdPlocal[nEquation][nEquation];
    RDouble dFluxdgradxL[nEquation][nEquation+nTemperatureModel];
    RDouble dFluxdgradyL[nEquation][nEquation+nTemperatureModel];
    RDouble dFluxdgradzL[nEquation][nEquation+nTemperatureModel];
    RDouble dFluxdgradxR[nEquation][nEquation+nTemperatureModel];
    RDouble dFluxdgradyR[nEquation][nEquation+nTemperatureModel];
    RDouble dFluxdgradzR[nEquation][nEquation+nTemperatureModel];
    RDouble dfpert;
    RDouble dgradpert;
    RDouble perturbScale = 1e-6;

    vector<int> *neighborCells = grid->GMRESGetNeighborCells();
    vector<int> *neighborFaces = grid->GMRESGetNeighborFaces();
    vector<int> *neighborLR = grid->GMRESGetNeighborLR();
    RDouble **dqdcvL= new RDouble*[nEquation];
    RDouble **dqdcvR= new RDouble*[nEquation];
    RDouble **dqdcvN= new RDouble*[nEquation];    //! dqdcv for neighbors
    RDouble **dgraddqx = new RDouble*[nEquation+nTemperatureModel];
    RDouble **dgraddqy = new RDouble*[nEquation+nTemperatureModel];
    RDouble **dgraddqz = new RDouble*[nEquation+nTemperatureModel];

    for(int indexI=0; indexI < nEquation; indexI++)
    {
        dqdcvL[indexI]   = new RDouble[nEquation];
        dqdcvR[indexI]   = new RDouble[nEquation];
        dqdcvN[indexI]   = new RDouble[nEquation];
        dgraddqx[indexI] = new RDouble[nEquation];
        dgraddqy[indexI] = new RDouble[nEquation];
        dgraddqz[indexI] = new RDouble[nEquation];
        
        for(int indexJ =0; indexJ < nEquation; indexJ ++)
        {
            dqdcvL[indexI][indexJ] = 0.0;
            dqdcvR[indexI][indexJ] = 0.0;
            dqdcvN[indexI][indexJ] = 0.0;
        }
    }
    dgraddqx[nEquation]    = new RDouble[nEquation];
    dgraddqy[nEquation]    = new RDouble[nEquation];
    dgraddqz[nEquation]    = new RDouble[nEquation];

    int nMid = 0;

    //! GMRESBoundary
     if (localStart >= nBoundFace)
     {
         nMid = localStart;
     }
     else if (localEnd <= nBoundFace) 
     {
        //! If they are all boundary faces.
         nMid = localEnd;
     }
     else
     {
         //! Part of them are boundary faces.
         nMid = nBoundFace;
     }

    //! boundary faces
    for (int iFace = localStart; iFace < nMid; ++ iFace)
    {
        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];
        int jFace  = iFace - localStart;

        //! init the value
        for (int m = 0; m < nEquation; ++ m)
        {
            fL[m]    = primitiveVariable[m][le];
            fR[m]    = primitiveVariable[m][re];
            primL[m] = primitiveVariable[m][le];
            primR[m] = primitiveVariable[m][re];
            dqdxL[m] = dqdx[m][le];
            dqdyL[m] = dqdy[m][le];
            dqdzL[m] = dqdz[m][le];
            dqdxR[m] = dqdx[m][re];
            dqdyR[m] = dqdy[m][re];
            dqdzR[m] = dqdz[m][re];
        }
        dqdxL[nEquation] = dtdx[ITT][le];
        dqdyL[nEquation] = dtdy[ITT][le];
        dqdzL[nEquation] = dtdz[ITT][le];
        dqdxR[nEquation] = dtdx[ITT][re];
        dqdyR[nEquation] = dtdy[ITT][re];
        dqdzR[nEquation] = dtdz[ITT][re];

        //! GMRESCoupled
        if( viscousType == ONE_EQU )
        {
            fLturb = primitiveturb[0][le];
            fRturb = primitiveturb[0][re];
        }

        //! define the sequence of the independent variables
        for (int m = 0; m < nEquation; m++)
        {
            fL[m].diff(m, nIndependentVars);
            fR[m].diff(halfIndependentVars + m, nIndependentVars);
        }
        for (int m = 0; m < nEquation + nTemperatureModel; m++)
        {
            dqdxL[m].diff(nEquation + m, nIndependentVars);
            dqdyL[m].diff(nEquation + (nEquation + nTemperatureModel) + m, nIndependentVars);
            dqdzL[m].diff(nEquation + 2 * (nEquation + nTemperatureModel) + m, nIndependentVars);
            dqdxR[m].diff(halfIndependentVars + nEquation + m, nIndependentVars);
            dqdyR[m].diff(halfIndependentVars + nEquation + (nEquation + nTemperatureModel) + m, nIndependentVars);
            dqdzR[m].diff(halfIndependentVars + nEquation + 2 * (nEquation + nTemperatureModel) + m, nIndependentVars);
        }

        //! GMRESCoupled
        if( viscousType == ONE_EQU )
        {
            fLturb.diff(nIndependentVars-2,nIndependentVars);
            fRturb.diff(nIndependentVars-1,nIndependentVars);
        }

        if( viscousType == LAMINAR )
        {
            Cal_GMRES_Visflux_AD(fL, fR, dqdxL, dqdyL, dqdzL, 
                                dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                                gridIn, faceProxy, parameters); 
        }
        else if( viscousType == ONE_EQU )
        {
            Cal_GMRES_Visflux_AD_Coupled(fL, fR, fLturb,fRturb,dqdxL, dqdyL, dqdzL, 
                                dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                                gridIn, faceProxy, parameters); 
        }

        for (int m = 0; m < nLaminar; ++ m)
        {
            //! value information
            flux[m][jFace] = fluxp[m].val();
            //! AD information
            for (int n = 0; n < nLaminar; ++ n)
            {
                dFluxdpL[m][n] = fluxp[m].dx(n);
                dFluxdpR[m][n] = fluxp[m].dx(halfIndependentVars + n);
            }
            for (int n = 0; n < nEquation + nTemperatureModel; n++)
            {
                dFluxdgradxL[m][n] = fluxp[m].dx(nEquation + n);
                dFluxdgradyL[m][n] = fluxp[m].dx(nEquation + nEquation + nTemperatureModel + n);
                dFluxdgradzL[m][n] = fluxp[m].dx(nEquation + 2 * (nEquation + nTemperatureModel) + n);
                dFluxdgradxR[m][n] = fluxp[m].dx(halfIndependentVars + nEquation + n);
                dFluxdgradyR[m][n] = fluxp[m].dx(halfIndependentVars + nEquation + nEquation + nTemperatureModel + n);
                dFluxdgradzR[m][n] = fluxp[m].dx(halfIndependentVars + nEquation + 2 * (nEquation + nTemperatureModel) + n);
            }
        }

        if( viscousType == ONE_EQU )
        {
            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdnutL[m] = fluxp[m].dx(nIndependentVars-2);
                dFluxdnutR[m] = fluxp[m].dx(nIndependentVars-1);
            }

            int colidx1      = GetIndexOfBlockJacobianMatrix(AI, AJ, 1, re, le);
            int colidx2      = GetIndexOfBlockJacobianMatrix(AI, AJ, 1, le, le);
            int colidx3      = GetIndexOfBlockJacobianMatrix(AI, AJ, 1, re, re);
            int colidx4      = GetIndexOfBlockJacobianMatrix(AI, AJ, 1, le, re);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                dRdqCoupledTerm[indexI][colidx1] -= dFluxdnutL[indexI];
                dRdqCoupledTerm[indexI][colidx2] += dFluxdnutL[indexI];
                dRdqCoupledTerm[indexI][colidx3] -= dFluxdnutR[indexI];
                dRdqCoupledTerm[indexI][colidx4] += dFluxdnutR[indexI];
                dRdqCoupledTerm[indexI][colidx2] += dFluxdnutR[indexI]*dDdP_turb[0][re-nTotalCell];
            }
        }
        //! assemble the matrix
        //! contributions from the bc
        int indexre = (re - nTotalCell) * nEquation;
        for (int indexI = 0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
            }
        }

        //! GMRESPV , obtain the dqdcv for the left cell
        gas->dPrimitive2dConservative(primL,gama[le],dqdcvL);

        //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
        //! then minus to the corresponding location in dRdq.
        int indexf = re * nEquation;    //! first index
        int indexs = le * nEquation;    //! second index
        colidx     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
        colidx1st  = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, re, le);    //! GMRESJac1st
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx + indexJ]   -= dFluxdpL[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];

                    dRdq1st[indexI][colidx1st + indexJ] -= dFluxdpL[indexI][indexK]*
                                                                dqdcvL[indexK][indexJ];    //! GMRESJac1st
                }
            }
        }

        //! contributions from the bc, 
        //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dFluxdpL[indexI][indexJ]  += dFluxdpR[indexI][indexK]*
                                                dDdPlocal[indexK][indexJ];
                }
            }
        }

        //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = le * nEquation;    //! second index
        colidx      = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
        colidx1st   = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, le);    //! GMRESJac1st
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx+indexJ] += dFluxdpL[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];

                    dRdq1st[indexI][colidx1st+indexJ] += dFluxdpL[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ]; // GMRESJac1st
                }
            }
        }

        //! considering the perturbation of the gradient in the flux calculation
        RDouble nxs    = xfn[iFace];
        RDouble nys    = yfn[iFace];
        RDouble nzs    = zfn[iFace];
        RDouble ns     = area[iFace];
        RDouble volume = vol[re];
        //dgradtmp = dfluxdgradxR*dgradRdqLx + dfluxdgradyR*dgradRdqLy + dfluxdgradzR*dgraddqLz
        //dRdq[le][le] += dgradtmp*dqdcvL
        gas->dGradient2dPrimitive(primL,-1,dgraddqx, 'x', nxs, nys,nzs,ns,volume);
        gas->dGradient2dPrimitive(primL,-1,dgraddqy, 'y', nxs, nys,nzs,ns,volume);
        gas->dGradient2dPrimitive(primL,-1,dgraddqz, 'z', nxs, nys,nzs,ns,volume);

        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dGradtmp[indexI][indexJ] = 0;
                
                for(int indexK=0; indexK < (nEquation+1); indexK++)
                {
                    dGradtmp[indexI][indexJ] += dFluxdgradxR[indexI][indexK]*
                                                             dgraddqx[indexK][indexJ];
                    dGradtmp[indexI][indexJ] += dFluxdgradyR[indexI][indexK]*
                                                             dgraddqy[indexK][indexJ];
                    dGradtmp[indexI][indexJ] += dFluxdgradzR[indexI][indexK]*
                                                             dgraddqz[indexK][indexJ];
                }
            }
        }

        //! convert dGradtmp by right multiplying the dqdcvL,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = le * nEquation;    //! second index
        colidx      = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
        colidx1st   = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, le);    //! GMRESJac1st
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx+indexJ] += dGradtmp[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];

                    dRdq1st[indexI][colidx1st+indexJ] += dGradtmp[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ]; // GMRESJac1st
                }
            }
        }

        indexf = re * nEquation;    //! first index
        indexs = le * nEquation;    //! second index
        colidx      = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
        colidx1st   = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, re, le);    //! GMRESJac1st
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx+indexJ] -= dGradtmp[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];

                    dRdq1st[indexI][colidx1st+indexJ] -= dGradtmp[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ]; // GMRESJac1st
                }
            }
        }

        nxs    = xfn[iFace];
        nys    = yfn[iFace];
        nzs    = zfn[iFace];
        ns     = area[iFace];
        volume = vol[le];
        // dgradtmp = dfluxdgradxL*dgradLdqRx + dfluxdgradyL*dgradLdqRy + dfluxdgradzL*dgradLdqRz
        // dRdq[le][re] += dgradtmp*dqdcvL
        gas->dGradient2dPrimitive(primR,1,dgraddqx, 'x', nxs, nys,nzs,ns,volume);
        gas->dGradient2dPrimitive(primR,1,dgraddqy, 'y', nxs, nys,nzs,ns,volume);
        gas->dGradient2dPrimitive(primR,1,dgraddqz, 'z', nxs, nys,nzs,ns,volume);

        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dGradtmp[indexI][indexJ] = 0;
                
                for(int indexK=0; indexK < (nEquation+1); indexK++)
                {
                    dGradtmp[indexI][indexJ] += dFluxdgradxL[indexI][indexK]*
                                                             dgraddqx[indexK][indexJ];
                    dGradtmp[indexI][indexJ] += dFluxdgradyL[indexI][indexK]*
                                                             dgraddqy[indexK][indexJ];
                    dGradtmp[indexI][indexJ] += dFluxdgradzL[indexI][indexK]*
                                                             dgraddqz[indexK][indexJ];
                }
            }
        }

        //! contributions from the bc, 
        //! dFluxdPR*dDdPlocal
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dGradtmpBC[indexI][indexJ] = 0.0;
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dGradtmpBC[indexI][indexJ]  += dGradtmp[indexI][indexK]*
                                                    dDdPlocal[indexK][indexJ];
                }
            }
        }

        //! convert dGradtmpBC by right multiplying the dqdcvL,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation; // first index
        indexs = le * nEquation; // second index
        colidx      = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
        colidx1st   = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, le);    //! GMRESJac1st
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx+indexJ] += dGradtmpBC[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];

                    dRdq1st[indexI][colidx1st+indexJ] += dGradtmpBC[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];  // GMRESJac1st
                }
            }
        }

        //! the neighbors of the left cells
        for(int index = 0; index < neighborCells[le].size(); index++)
        {
            int neighborIndex = neighborCells[le][index];

            if(neighborIndex != re)
            {
                int Faceindex   = neighborFaces[le][index];
                //! judge whether it is the boundary face and obtain its boundary type
                bool isGhostCell = false;
                UnstructBC *neighborBCRegion = nullptr;
                int neighborBCType;
                if (neighborIndex >= nTotalCell)
                {
                    isGhostCell = true;
                    neighborBCRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[Faceindex]);
                    neighborBCType = neighborBCRegion->GetBCType();
                }

                int sign       = neighborLR[le][index];
                RDouble nxs    = xfn[Faceindex];
                RDouble nys    = yfn[Faceindex];
                RDouble nzs    = zfn[Faceindex];
                RDouble ns     = area[Faceindex];
                RDouble volume = vol[le];
                for(int m =0; m < nEquation; m++)
                {
                    fN[m] = primitiveVariable[m][neighborIndex];
                }

                gas->dPrimitive2dConservative(fN,gama[neighborIndex],dqdcvN);

                gas->dGradient2dPrimitive(fN,sign, dgraddqx, 'x', nxs, nys,nzs,ns,volume);
                gas->dGradient2dPrimitive(fN,sign, dgraddqy, 'y', nxs, nys,nzs,ns,volume);
                gas->dGradient2dPrimitive(fN,sign, dgraddqz, 'z', nxs, nys,nzs,ns,volume);

                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dGradtmp[indexI][indexJ] = 0;

                        for(int indexK=0; indexK < (nEquation+1); indexK++)
                        {
                            dGradtmp[indexI][indexJ] += dFluxdgradxL[indexI][indexK]*
                                                                     dgraddqx[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradyL[indexI][indexK]*
                                                                     dgraddqy[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradzL[indexI][indexK]*
                                                                     dgraddqz[indexK][indexJ];
                        }
                    }
                }

                int indexf  = le * nEquation;    //! first index
                int indexs  = neighborIndex * nEquation;    //! second index
                int indexf2 = re * nEquation; 
                colidx  = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {

                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            // dRL/dqLNeighbors += dRL/dGradL*dGradL/dqLNeighbors ===> dFluxdgradL*dgraddq
                            dRdq[indexI][colidx+indexJ] += dGradtmp[indexI][indexK]*
                                                                     dqdcvN[indexK][indexJ];
                            
                            // dRR/dqLNeighbors -= dRR/dGradL*dGradL/dqLNeighbors ===> dFluxdgradL*dgraddq
                            dRdq[indexI][colidx2+indexJ] -= dGradtmp[indexI][indexK]*
                                                                     dqdcvN[indexK][indexJ];
                        }
                    }
                } 

                if (isGhostCell && neighborBCType != PHENGLEI::INTERFACE)
                {
                    int indexre = (neighborIndex - nTotalCell) * nEquation;
                    for (int indexI = 0; indexI < nEquation; indexI++)
                    {
                        for (int indexJ = 0; indexJ < nEquation; indexJ++)
                        {
                            dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre + indexJ];
                        }
                    }

                    for (int indexI = 0; indexI < nEquation; indexI++)
                    {
                        for (int indexJ = 0; indexJ < nEquation; indexJ++)
                        {
                            dGradtmpBC[indexI][indexJ] = 0.0;
                            for (int indexK = 0; indexK < nEquation; indexK++)
                            {
                                dGradtmpBC[indexI][indexJ] += dGradtmp[indexI][indexK] * dDdPlocal[indexK][indexJ];
                            }
                        }
                    }

                    gas->dPrimitive2dConservative(primL, gama[le], dqdcvL);
                    int colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                    int colidx1st = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, le);
                    int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                    int colidx21st = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, re, le);
                    for (int indexI = 0; indexI < nEquation; indexI++)
                    {
                        for (int indexJ = 0; indexJ < nEquation; indexJ++)
                        {
                            for (int indexK = 0; indexK < nEquation; indexK++)
                            {
                                RDouble tmp = dGradtmpBC[indexI][indexK] * dqdcvL[indexK][indexJ];
                                dRdq[indexI][colidx + indexJ] += tmp;
                                dRdq[indexI][colidx2 + indexJ] -= tmp;
                                dRdq1st[indexI][colidx1st + indexJ] += tmp;
                                dRdq1st[indexI][colidx21st + indexJ] -= tmp;
                            }
                        }
                    }
                }
            }
        }
    }

    //! interior faces
    for (int iFace = nMid; iFace < localEnd; ++ iFace)
    {
        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];
        int jFace  = iFace - localStart;

        //! init the value
        for (int m = 0; m < nEquation; ++ m)
        {
            fL[m]    = primitiveVariable[m][le];
            fR[m]    = primitiveVariable[m][re];
            primL[m] = primitiveVariable[m][le];
            primR[m] = primitiveVariable[m][re];
            dqdxL[m] = dqdx[m][le];
            dqdyL[m] = dqdy[m][le];
            dqdzL[m] = dqdz[m][le];
            dqdxR[m] = dqdx[m][re];
            dqdyR[m] = dqdy[m][re];
            dqdzR[m] = dqdz[m][re];
        }
        dqdxL[nEquation] = dtdx[ITT][le];
        dqdyL[nEquation] = dtdy[ITT][le];
        dqdzL[nEquation] = dtdz[ITT][le];
        dqdxR[nEquation] = dtdx[ITT][re];
        dqdyR[nEquation] = dtdy[ITT][re];
        dqdzR[nEquation] = dtdz[ITT][re];

        //! GMRESCoupled
        if( viscousType == ONE_EQU )
        {
            fLturb = primitiveturb[0][le];
            fRturb = primitiveturb[0][re];
        }

        //! define the sequence of the independent variables
        for (int m = 0; m < nEquation; m++)
        {
            fL[m].diff(m, nIndependentVars);
            fR[m].diff(halfIndependentVars + m, nIndependentVars);
        }
        for (int m = 0; m < nEquation + nTemperatureModel; m++)
        {
            dqdxL[m].diff(nEquation + m, nIndependentVars);
            dqdyL[m].diff(nEquation + (nEquation + nTemperatureModel) + m, nIndependentVars);
            dqdzL[m].diff(nEquation + 2 * (nEquation + nTemperatureModel) + m, nIndependentVars);
            dqdxR[m].diff(halfIndependentVars + nEquation + m, nIndependentVars);
            dqdyR[m].diff(halfIndependentVars + nEquation + (nEquation + nTemperatureModel) + m, nIndependentVars);
            dqdzR[m].diff(halfIndependentVars + nEquation + 2 * (nEquation + nTemperatureModel) + m, nIndependentVars);
        }

        //! GMRESCoupled
        if( viscousType == ONE_EQU )
        {
            fLturb.diff(nIndependentVars-2,nIndependentVars);
            fRturb.diff(nIndependentVars-1,nIndependentVars);
        }

        // std::cout << "before AD " << iFace << "\n";

        if( viscousType == LAMINAR )
        {
            Cal_GMRES_Visflux_AD(fL, fR, dqdxL, dqdyL, dqdzL, 
                                dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                                gridIn, faceProxy, parameters); 
        }
        else if( viscousType == ONE_EQU )
        {
            Cal_GMRES_Visflux_AD_Coupled(fL, fR, fLturb,fRturb,dqdxL, dqdyL, dqdzL, 
                                dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                                gridIn, faceProxy, parameters); 
        }


        // std::cout << "after AD " << iFace << "\n";

        for (int m = 0; m < nLaminar; ++ m)
        {
            //! value information
            flux[m][jFace] = fluxp[m].val();
            //! AD information
            for (int n = 0; n < nLaminar; ++ n)
            {
                dFluxdpL[m][n] = fluxp[m].dx(n);
                dFluxdpR[m][n] = fluxp[m].dx(halfIndependentVars + n);
            }
            for (int n = 0; n < nEquation + nTemperatureModel; n++)
            {
                dFluxdgradxL[m][n] = fluxp[m].dx(nEquation + n);
                dFluxdgradyL[m][n] = fluxp[m].dx(nEquation + nEquation + nTemperatureModel + n);
                dFluxdgradzL[m][n] = fluxp[m].dx(nEquation + 2 * (nEquation + nTemperatureModel) + n);
                dFluxdgradxR[m][n] = fluxp[m].dx(halfIndependentVars + nEquation + n);
                dFluxdgradyR[m][n] = fluxp[m].dx(halfIndependentVars + nEquation + nEquation + nTemperatureModel + n);
                dFluxdgradzR[m][n] = fluxp[m].dx(halfIndependentVars + nEquation + 2 * (nEquation + nTemperatureModel) + n);
            }
        }

        if( viscousType == ONE_EQU )
        {
            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdnutL[m] = fluxp[m].dx(nIndependentVars-2);
                dFluxdnutR[m] = fluxp[m].dx(nIndependentVars-1);
            }

            int colidx1      = GetIndexOfBlockJacobianMatrix(AI, AJ, 1, re, le);
            int colidx2      = GetIndexOfBlockJacobianMatrix(AI, AJ, 1, le, le);
            int colidx3      = GetIndexOfBlockJacobianMatrix(AI, AJ, 1, re, re);
            int colidx4      = GetIndexOfBlockJacobianMatrix(AI, AJ, 1, le, re);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                dRdqCoupledTerm[indexI][colidx1] -= dFluxdnutL[indexI];
                dRdqCoupledTerm[indexI][colidx2] += dFluxdnutL[indexI];
                dRdqCoupledTerm[indexI][colidx3] -= dFluxdnutR[indexI];
                dRdqCoupledTerm[indexI][colidx4] += dFluxdnutR[indexI];
            }
        }

        //! GMRESPV , obtain the dqdcv for the left cell
        gas->dPrimitive2dConservative(primL,gama[le],dqdcvL);

        //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
        //! then minus to the corresponding location in dRdq.
        int indexf = re * nEquation;    //! first index
        int indexs = le * nEquation;    //! second index
        colidx      = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
        colidx1st   = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, re, le);    //! GMRESJac1st
        for (int indexI = 0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx+indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];

                    dRdq1st[indexI][colidx1st+indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];    //! GMRESJac1st
                }
            }
        }

        //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = le * nEquation;    //! second index
        colidx      = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
        colidx1st   = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, le);    //! GMRESJac1st
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx+indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];

                    dRdq1st[indexI][colidx1st+indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];    //! GMRESJac1st
                }
            }
        }

         //! GMRESPV , obtain the dqdcv for the right cell
        gas->dPrimitive2dConservative(primR,gama[re],dqdcvR);

        //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
        //! then minus to the corresponding location in dRdq.
        indexf = re * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        colidx    = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
        colidx1st = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, re, re);    //! GMRESJac1st
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx+indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];

                    dRdq1st[indexI][colidx1st+indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ]; // GMRESJac1st
                }
            }
        }

        //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        colidx      = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
        colidx1st   = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, re);    //! GMRESJac1st
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx+indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];

                    dRdq1st[indexI][colidx1st+indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];    //! GMRESJac1st
                }
            }
        }

        //! considering the perturbation of the gradient in the flux calculation
        RDouble nxs         = xfn[iFace];
        RDouble nys         = yfn[iFace];
        RDouble nzs         = zfn[iFace];
        RDouble ns          = area[iFace];
        RDouble volume      = vol[re];
        // dgradtmp = dfluxdgradxR*dgradRdqLx + dfluxdgradyR*dgradRdqLy + dfluxdgradzR*dgraddqLz
        // dRdq[le][le] += dgradtmp*dqdcvL
        gas->dGradient2dPrimitive(primL,-1,dgraddqx, 'x', nxs, nys,nzs,ns,volume);
        gas->dGradient2dPrimitive(primL,-1,dgraddqy, 'y', nxs, nys,nzs,ns,volume);
        gas->dGradient2dPrimitive(primL,-1,dgraddqz, 'z', nxs, nys,nzs,ns,volume);

        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dGradtmp[indexI][indexJ] = 0;
                for(int indexK=0; indexK < (nEquation+1); indexK++)
                {
                    dGradtmp[indexI][indexJ] += dFluxdgradxR[indexI][indexK] * dgraddqx[indexK][indexJ];
                    dGradtmp[indexI][indexJ] += dFluxdgradyR[indexI][indexK] * dgraddqy[indexK][indexJ];
                    dGradtmp[indexI][indexJ] += dFluxdgradzR[indexI][indexK] * dgraddqz[indexK][indexJ];
                }
            }
        }

        //! convert dGradtmp by right multiplying the dqdcvL,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = le * nEquation;    //! second index
        colidx      = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
        colidx1st   = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, le);    //! GMRESJac1st
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx+indexJ] += dGradtmp[indexI][indexK] * dqdcvL[indexK][indexJ];

                    dRdq1st[indexI][colidx1st+indexJ] += dGradtmp[indexI][indexK] * dqdcvL[indexK][indexJ];    //! GMRESJac1st
                }
            }
        }

        //! dgradtmp = dfluxdgradxR*dgradRdqLx + dfluxdgradyR*dgradRdqLy + dfluxdgradzR*dgraddqLz
        //! dRdq[re][le] += dgradtmp*dqdcvL
        indexf = re * nEquation;    //! first index
        indexs = le * nEquation;    //! second index
        colidx      = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
        colidx1st   = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, re, le);    //! GMRESJac1st
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx+indexJ] -= dGradtmp[indexI][indexK] * dqdcvL[indexK][indexJ];

                    dRdq1st[indexI][colidx1st+indexJ] -= dGradtmp[indexI][indexK] * dqdcvL[indexK][indexJ]; // GMRESJac1st
                }
            }
        }

        nxs         = xfn[iFace];
        nys         = yfn[iFace];
        nzs         = zfn[iFace];
        ns          = area[iFace];
        volume      = vol[le];
        //! dgradtmp = dfluxdgradxL*dgradLdqRx + dfluxdgradyL*dgradLdqRy + dfluxdgradzL*dgradLdqRz
        //! dRdq[le][re] += dgradtmp*dqdcvL
        gas->dGradient2dPrimitive(primR,1,dgraddqx, 'x', nxs, nys,nzs,ns,volume);
        gas->dGradient2dPrimitive(primR,1,dgraddqy, 'y', nxs, nys,nzs,ns,volume);
        gas->dGradient2dPrimitive(primR,1,dgraddqz, 'z', nxs, nys,nzs,ns,volume);

        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dGradtmp[indexI][indexJ] = 0;
                
                for(int indexK=0; indexK < (nEquation+1); indexK++)
                {
                    dGradtmp[indexI][indexJ] += dFluxdgradxL[indexI][indexK]*
                                                             dgraddqx[indexK][indexJ];
                    dGradtmp[indexI][indexJ] += dFluxdgradyL[indexI][indexK]*
                                                             dgraddqy[indexK][indexJ];
                    dGradtmp[indexI][indexJ] += dFluxdgradzL[indexI][indexK]*
                                                             dgraddqz[indexK][indexJ];
                }
            }
        }

        //! convert dGradtmp by right multiplying the dqdcvL,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        colidx      = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
        colidx1st   = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, re);    //! GMRESJac1st
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx+indexJ] += dGradtmp[indexI][indexK]*
                                                             dqdcvR[indexK][indexJ];

                    dRdq1st[indexI][colidx1st+indexJ] += dGradtmp[indexI][indexK]*
                                                             dqdcvR[indexK][indexJ];    //! GMRESJac1st
                }
            }
        }

        // dgradtmp = dfluxdgradxL*dgradLdqRx + dfluxdgradyL*dgradLdqRy + dfluxdgradzL*dgradLdqRz
        // dRdq[re][re] += dgradtmp*dqdcvL
        indexf = re * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        colidx      = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
        colidx1st   = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, re, re);    //! GMRESJac1st
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx+indexJ] -= dGradtmp[indexI][indexK]*
                                                             dqdcvR[indexK][indexJ];

                    dRdq1st[indexI][colidx1st+indexJ] -= dGradtmp[indexI][indexK]*
                                                             dqdcvR[indexK][indexJ];    //! GMRESJac1st
                }
            }
        }

        //! the neighbors of the left cells
        for(int index = 0; index < neighborCells[le].size(); index++)
        {
            int neighborIndex = neighborCells[le][index];

            if(neighborIndex != re)
            {
                int Faceindex   = neighborFaces[le][index];
                //! judge whether it is the boundary face and obtain its boundary type
                bool isGhostCell = false;
                UnstructBC *neighborBCRegion = nullptr;
                int neighborBCType;
                if (neighborIndex >= nTotalCell)
                {
                    isGhostCell = true;
                    neighborBCRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[Faceindex]);
                    neighborBCType = neighborBCRegion->GetBCType();
                }
                
                int sign       = neighborLR[le][index];
                RDouble nxs    = xfn[Faceindex];
                RDouble nys    = yfn[Faceindex];
                RDouble nzs    = zfn[Faceindex];
                RDouble ns     = area[Faceindex];
                RDouble volume = vol[le];

                for(int m =0; m < nEquation; m++)
                {
                    fN[m] = primitiveVariable[m][neighborIndex];
                }

                gas->dPrimitive2dConservative(fN,gama[neighborIndex],dqdcvN);

                gas->dGradient2dPrimitive(fN,sign,dgraddqx, 'x', nxs, nys,nzs,ns,volume);
                gas->dGradient2dPrimitive(fN,sign,dgraddqy, 'y', nxs, nys,nzs,ns,volume);
                gas->dGradient2dPrimitive(fN,sign,dgraddqz, 'z', nxs, nys,nzs,ns,volume);

                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dGradtmp[indexI][indexJ] = 0;

                        for(int indexK=0; indexK < (nEquation+1); indexK++)
                        {
                            dGradtmp[indexI][indexJ] += dFluxdgradxL[indexI][indexK]*
                                                                     dgraddqx[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradyL[indexI][indexK]*
                                                                     dgraddqy[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradzL[indexI][indexK]*
                                                                     dgraddqz[indexK][indexJ];
                        }
                    }
                }

                int indexf  = le * nEquation;    //! first index
                int indexs  = neighborIndex * nEquation;    //! second index
                int indexf2 = re * nEquation; 
                colidx  = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {

                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            // dRL/dqLNeighbors += dRL/dGradL*dGradL/dqLNeighbors ===> dFluxdgradL*dgraddq
                            dRdq[indexI][colidx+indexJ] += dGradtmp[indexI][indexK]*
                                                                     dqdcvN[indexK][indexJ];

                            // dRR/dqLNeighbors -= dRR/dGradL*dGradL/dqLNeighbors ===> dFluxdgradL*dgraddq
                            dRdq[indexI][colidx2+indexJ] -= dGradtmp[indexI][indexK]*
                                                                     dqdcvN[indexK][indexJ];
                        }
                    }
                }

                //! boundary condition except interface
                if (isGhostCell && neighborBCType != PHENGLEI::INTERFACE)
                {
                    int indexre = (neighborIndex - nTotalCell) * nEquation;
                    for (int indexI = 0; indexI < nEquation; indexI++)
                    {
                        for (int indexJ = 0; indexJ < nEquation; indexJ++)
                        {
                            dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre + indexJ];
                        }
                    }

                    for (int indexI = 0; indexI < nEquation; indexI++)
                    {
                        for (int indexJ = 0; indexJ < nEquation; indexJ++)
                        {
                            dGradtmpBC[indexI][indexJ] = 0.0;
                            for (int indexK = 0; indexK < nEquation; indexK++)
                            {
                                dGradtmpBC[indexI][indexJ] += dGradtmp[indexI][indexK] * dDdPlocal[indexK][indexJ];
                            }
                        }
                    }

                    gas->dPrimitive2dConservative(primL, gama[le], dqdcvL);
                    int colidx     = GetIndexOfBlockJacobianMatrix(AI   , AJ   , nEquation, le, le);
                    int colidx1st  = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, le);
                    int colidx2    = GetIndexOfBlockJacobianMatrix(AI,    AJ,    nEquation, re, le);
                    int colidx21st = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, re, le);
                    for (int indexI = 0; indexI < nEquation; indexI++)
                    {
                        for (int indexJ = 0; indexJ < nEquation; indexJ++)
                        {
                            for (int indexK = 0; indexK < nEquation; indexK++)
                            {
                                RDouble tmp = dGradtmpBC[indexI][indexK] * dqdcvL[indexK][indexJ];
                                dRdq[indexI][colidx + indexJ] += tmp;
                                dRdq[indexI][colidx2 + indexJ] -= tmp;
                                dRdq1st[indexI][colidx1st + indexJ] += tmp;
                                dRdq1st[indexI][colidx21st + indexJ] -= tmp;
                            }
                        }
                    }
                }
            }
        }

        //! the neighbors of the right cells
        for(int index = 0; index < neighborCells[re].size(); index++)
        {
            int neighborIndex = neighborCells[re][index];

            if(neighborIndex != le)
            {
                int Faceindex       = neighborFaces[re][index];
                 //! judge whether it is the boundary face and obtain its boundary type
                bool isGhostCell = false;
                UnstructBC *neighborBCRegion = nullptr;
                int neighborBCType;
                if (neighborIndex >= nTotalCell)
                {
                    isGhostCell = true;
                    neighborBCRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[Faceindex]);
                    neighborBCType = neighborBCRegion->GetBCType();
                }

                int sign = neighborLR[re][index];
                RDouble nxs    = xfn[Faceindex];
                RDouble nys    = yfn[Faceindex];
                RDouble nzs    = zfn[Faceindex];
                RDouble ns     = area[Faceindex];
                RDouble volume = vol[re];

                for(int m=0; m < nEquation; m++)
                {
                    fN[m] = primitiveVariable[m][neighborIndex];
                }

                gas->dPrimitive2dConservative(fN,gama[neighborIndex],dqdcvN);

                gas->dGradient2dPrimitive(fN,sign,dgraddqx, 'x', nxs, nys,nzs,ns,volume);
                gas->dGradient2dPrimitive(fN,sign,dgraddqy, 'y', nxs, nys,nzs,ns,volume);
                gas->dGradient2dPrimitive(fN,sign,dgraddqz, 'z', nxs, nys,nzs,ns,volume);

                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dGradtmp[indexI][indexJ] = 0;
                        for(int indexK=0; indexK < (nEquation+1); indexK++)
                        {
                            dGradtmp[indexI][indexJ] += dFluxdgradxR[indexI][indexK]*
                                                                     dgraddqx[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradyR[indexI][indexK]*
                                                                     dgraddqy[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradzR[indexI][indexK]*
                                                                     dgraddqz[indexK][indexJ];
                        }
                    }
                }

                int indexf  = le * nEquation;    //! first index
                int indexs  = neighborIndex * nEquation;    //! second index
                int indexf2 = re * nEquation; 
                colidx  = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            // dRL/dqRNeighbors += dRL/dGradR*dGradR/dqRNeighbors ===> dFluxdgradR*dgraddq
                            dRdq[indexI][colidx+indexJ] += dGradtmp[indexI][indexK]*
                                                                     dqdcvN[indexK][indexJ];


                            // dRR/dqRNeighbors -= dRR/dGradR*dGradR/dqRNeighbors ===> dFluxdgradR*dgraddq
                            dRdq[indexI][colidx2+indexJ] -= dGradtmp[indexI][indexK]*
                                                                     dqdcvN[indexK][indexJ];
                        }
                    }
                }

                if (isGhostCell && neighborBCType != PHENGLEI::INTERFACE)
                {
                    int indexre = (neighborIndex - nTotalCell) * nEquation;
                    for (int indexI = 0; indexI < nEquation; indexI++)
                    {
                        for (int indexJ = 0; indexJ < nEquation; indexJ++)
                        {
                            dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre + indexJ];
                        }
                    }

                    for (int indexI = 0; indexI < nEquation; indexI++)
                    {
                        for (int indexJ = 0; indexJ < nEquation; indexJ++)
                        {
                            dGradtmpBC[indexI][indexJ] = 0.0;
                            for (int indexK = 0; indexK < nEquation; indexK++)
                            {
                                dGradtmpBC[indexI][indexJ] += dGradtmp[indexI][indexK] * dDdPlocal[indexK][indexJ];
                            }
                        }
                    }

                    gas->dPrimitive2dConservative(primR, gama[re], dqdcvR);
                    int colidx     = GetIndexOfBlockJacobianMatrix(AI   , AJ   , nEquation, le, re);
                    int colidx1st  = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, re);
                    int colidx2    = GetIndexOfBlockJacobianMatrix(AI   , AJ   , nEquation, re, re);
                    int colidx21st = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, re, re);
                    for (int indexI = 0; indexI < nEquation; indexI++)
                    {
                        for (int indexJ = 0; indexJ < nEquation; indexJ++)
                        {
                            for (int indexK = 0; indexK < nEquation; indexK++)
                            {
                                RDouble tmp = dGradtmp[indexI][indexK] * dqdcvR[indexK][indexJ];
                                dRdq[indexI][colidx + indexJ] += tmp;
                                dRdq[indexI][colidx2 + indexJ] -= tmp;
                                dRdq1st[indexI][colidx1st + indexJ] += tmp;
                                dRdq1st[indexI][colidx21st + indexJ] -= tmp;
                            }
                        }
                    }
                }
            }
        }
    }
    delete [] dqdxL;    dqdxL = nullptr;
    delete [] dqdyL;    dqdyL = nullptr;
    delete [] dqdzL;    dqdzL = nullptr;

    delete [] dqdxR;    dqdxR = nullptr;
    delete [] dqdyR;    dqdyR = nullptr;
    delete [] dqdzR;    dqdzR = nullptr;

    delete [] fL;       fL = nullptr;
    delete [] fR;       fR = nullptr;
    delete [] fN;       fN = nullptr;
    delete [] fluxp;    fluxp = nullptr;
    delete [] primL;    primL = nullptr;
    delete [] primR;    primR = nullptr;

    for(int index = 0; index < nEquation; index++)
    {
        delete [] dqdcvL[index];
        delete [] dqdcvR[index];
        delete [] dqdcvN[index];
        delete [] dgraddqx[index];
        delete [] dgraddqy[index];
        delete [] dgraddqz[index];
    }
    delete [] dgraddqx[nEquation];
    delete [] dgraddqy[nEquation];
    delete [] dgraddqz[nEquation];
    delete [] dqdcvL;    dqdcvL = nullptr;
    delete [] dqdcvR;    dqdcvR = nullptr;
    delete [] dqdcvN;    dqdcvN = nullptr;
    delete [] dgraddqx;    dgraddqx = nullptr;
    delete [] dgraddqy;    dgraddqy = nullptr;
    delete [] dgraddqz;    dgraddqz = nullptr;

    delete [] dFluxdnutL;    dFluxdnutL = nullptr;
    delete [] dFluxdnutR;    dFluxdnutR = nullptr;
}

//! GMRESVis
void NSSolverUnstruct::ComputeVisflux_FD(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();

    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    RDouble **dqdx = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarZ"));

    RDouble **dtdx = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradTemperatureX"));
    RDouble **dtdy = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradTemperatureY"));
    RDouble **dtdz = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradTemperatureZ"));

    RDouble  **flux = faceProxy->GetFlux();

    RDouble *deltaL = faceProxy->GetWeightL();
    RDouble *deltaR = faceProxy->GetWeightR();

    NSFaceValue *faceVariable = faceProxy->GetNSFaceValue();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble *vol  = grid->GetCellVolume();
    
    RDouble *xcc  = grid->GetCellCenterX();
    RDouble *ycc  = grid->GetCellCenterY();
    RDouble *zcc  = grid->GetCellCenterZ();

    RDouble *xfc  = grid->GetFaceCenterX();
    RDouble *yfc  = grid->GetFaceCenterY();
    RDouble *zfc  = grid->GetFaceCenterZ();

    RDouble **rhoDiffusion = faceVariable->GetRhoDS()->GetField();
    RDouble **hintSpecies = faceVariable->GetHintS()->GetField();
    RDouble **primitiveVariableFace = faceVariable->GetPrim()->GetField();
    RDouble **tm     = faceVariable->GetT()->GetField();

    RDouble *kCp     = faceVariable->GetKCP();
    RDouble *viscousLaminarFace   = faceVariable->GetMUL();
    RDouble *viscousTurbulenceFace = faceVariable->GetMUT();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int numberOfSpecies = parameters->GetNumberOfSpecies();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int nm = parameters->GetNSEquationNumber();

    RDouble refReNumber  = parameters->GetRefReNumber();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    int nolstress = GlobalDataBase::GetIntParaFromDB("nolstress");
    int nrokplus  = GlobalDataBase::GetIntParaFromDB("nrokplus");

    RDouble skewnessAngle = parameters->GetSkewnessAngle();

    RDouble **t                 = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));
    RDouble **primitiveVariable = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *gama               = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble** dRdq              = reinterpret_cast<RDouble **> (grid->GetDataPtr("dRdq"));
    RDouble** dDdP              = reinterpret_cast<RDouble **> (grid->GetDataPtr("dDdP"));

    RDouble **qTurb = 0;
    RDouble **aniss  = 0;
    RDouble coefficientofStateEquation = gas->GetCoefficientOfStateEquation();

    if (nrokplus > 0)
    {
        qTurb = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
        aniss  = reinterpret_cast<RDouble **> (grid->GetDataPtr("aniss"));
    }

    RDouble *dqdxL = new RDouble[nEquation+nTemperatureModel];
    RDouble *dqdyL = new RDouble[nEquation+nTemperatureModel];
    RDouble *dqdzL = new RDouble[nEquation+nTemperatureModel];

    RDouble *dqdxR = new RDouble[nEquation+nTemperatureModel];
    RDouble *dqdyR = new RDouble[nEquation+nTemperatureModel];
    RDouble *dqdzR = new RDouble[nEquation+nTemperatureModel];
    RDouble *dgrad = new RDouble[nEquation+nTemperatureModel];

    RDouble *fL       = new RDouble[nEquation+nTemperatureModel];
    RDouble *fLpert   = new RDouble[nEquation+nTemperatureModel];
    RDouble *fR       = new RDouble[nEquation+nTemperatureModel];
    RDouble *fRpert   = new RDouble[nEquation+nTemperatureModel];
    RDouble *fN       = new RDouble[nEquation]; // for neighbors
    RDouble *fMid     = new RDouble[nEquation+nTemperatureModel];
    RDouble *fMidpert = new RDouble[nEquation+nTemperatureModel];
    RDouble *fluxp    = new RDouble[nLaminar];

    RDouble t1x, t1y, t1z, t2x, t2y, t2z;
    RDouble txx,tyy,tzz;
    RDouble txy,txz,tyz;

    bool isFineGrid = grid->IsFinestGrid();

    using namespace IDX;
    RDouble  dFluxdpL[nEquation][nEquation];
    RDouble  dFluxdpR[nEquation][nEquation];
    RDouble  dGradtmp[nEquation][nEquation];
    RDouble  dGradtmpBC[nEquation][nEquation];
    RDouble  dDdPlocal[nEquation][nEquation];
    RDouble  dFluxdgradxL[nEquation][nEquation+nTemperatureModel];
    RDouble  dFluxdgradyL[nEquation][nEquation+nTemperatureModel];
    RDouble  dFluxdgradzL[nEquation][nEquation+nTemperatureModel];
    RDouble  dFluxdgradxR[nEquation][nEquation+nTemperatureModel];
    RDouble  dFluxdgradyR[nEquation][nEquation+nTemperatureModel];
    RDouble  dFluxdgradzR[nEquation][nEquation+nTemperatureModel];
    RDouble  dfpert;
    RDouble  dgradpert;
    RDouble  perturbScale = 1e-6;

    vector<int> *neighborCells = grid->GMRESGetNeighborCells();
    vector<int> *neighborFaces = grid->GMRESGetNeighborFaces();
    vector<int> *neighborLR = grid->GMRESGetNeighborLR();

    RDouble** dqdcvL= new RDouble*[nEquation];
    RDouble** dqdcvR= new RDouble*[nEquation];
    RDouble** dqdcvN= new RDouble*[nEquation];    //! dqdcv for neighbors
    RDouble** dgraddqx = new RDouble*[nEquation+nTemperatureModel];
    RDouble** dgraddqy = new RDouble*[nEquation+nTemperatureModel];
    RDouble** dgraddqz = new RDouble*[nEquation+nTemperatureModel];

    for(int indexI=0; indexI < nEquation; indexI++)
    {
        dqdcvL[indexI]      = new RDouble[nEquation];
        dqdcvR[indexI]      = new RDouble[nEquation];
        dqdcvN[indexI]      = new RDouble[nEquation];
        dgraddqx[indexI]    = new RDouble[nEquation];
        dgraddqy[indexI]    = new RDouble[nEquation];
        dgraddqz[indexI]    = new RDouble[nEquation];
        
        for(int indexJ =0; indexJ < nEquation; indexJ ++)
        {
            dqdcvL[indexI][indexJ] = 0.0;
            dqdcvR[indexI][indexJ] = 0.0;
            dqdcvN[indexI][indexJ] = 0.0;
        }
    }
    dgraddqx[nEquation] = new RDouble[nEquation];
    dgraddqy[nEquation] = new RDouble[nEquation];
    dgraddqz[nEquation] = new RDouble[nEquation];

    int nMid;
        //! GMRESBoundary
    if (localStart >= nBoundFace)
    {
        nMid = localStart;    //! a bug 01.07
    }
    else if (localEnd <= nBoundFace)    //! add else 
    {
        //! If they are all boundary faces.
        nMid = localEnd;
    }
    else
    {
        //! Part of them are boundary faces.
        nMid = nBoundFace;
    }

    //! boundary faces
    for (int iFace = localStart; iFace < nMid; ++ iFace)
    {
        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];
        int jFace  = iFace - localStart;

        for (int m = 0; m < nEquation; ++ m)
        {
            fL[m]       = primitiveVariable[m][le];
            fR[m]       = primitiveVariable[m][re];
            fMid[m]     = primitiveVariableFace[m][jFace];
            dqdxL[m]    = dqdx[m][le];
            dqdyL[m]    = dqdy[m][le];
            dqdzL[m]    = dqdz[m][le];
            dqdxR[m]    = dqdx[m][re];
            dqdyR[m]    = dqdy[m][re];
            dqdzR[m]    = dqdz[m][re];
        }
        fL[nEquation]           = t[ITT][le];
        fLpert[nEquation]       = t[ITT][le];
        fR[nEquation]           = t[ITT][re];
        fRpert[nEquation]       = t[ITT][re];
        fMid[nEquation]         = tm[ITT][jFace];
        fMidpert[nEquation]     = tm[ITT][jFace];
        dqdxL[nEquation]        = dtdx[ITT][le];
        dqdyL[nEquation]        = dtdy[ITT][le];
        dqdzL[nEquation]        = dtdz[ITT][le];
        dqdxR[nEquation]        = dtdx[ITT][re];
        dqdyR[nEquation]        = dtdy[ITT][re];
        dqdzR[nEquation]        = dtdz[ITT][re];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            flux[m][jFace] = fluxp[m];
        }

        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        //! for solid walls, primitiveVariableFace[IU,IV,IW][jFace] = 0
        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            //! perturb the density of the left cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fLpert[index]   = fL[index];
                fMidpert[index] = fMid[index];  
            }
            fLpert[IR]          *=  (1.0 + perturbScale);
            fLpert[nEquation]   =   fLpert[IP]/(coefficientofStateEquation * fLpert[IR]);
            dfpert              =   fLpert[IR] - fL[IR];
            fMidpert[IR]        =   half * (fLpert[IR]+fR[IR]);
            fMidpert[nEquation] =   half * (fLpert[nEquation]+fR[nEquation]);

            Cal_GMRES_Visflux(fLpert, fR, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpL[m][IR] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            //! perturb the u of the left cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fLpert[index]   = fL[index];
                fMidpert[index] = fMid[index];  
            }
            fLpert[nEquation]   =   fL[nEquation];
            fMidpert[nEquation] =   half * (fLpert[nEquation]+fR[nEquation]);
            fLpert[IU]          +=  perturbScale;
            dfpert              =   fLpert[IU] - fL[IU];
            // fMidpert[IU]        =   half * (fLpert[IU]+fR[IU]);

            Cal_GMRES_Visflux(fLpert, fR, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpL[m][IU] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            //! perturb the v of the left cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fLpert[index]   = fL[index];
                fMidpert[index] = fMid[index];  
            }
            fLpert[nEquation]   =   fL[nEquation];
            fMidpert[nEquation] =   half * (fLpert[nEquation]+fR[nEquation]);
            fLpert[IV]          +=  perturbScale;
            dfpert              =   fLpert[IV] - fL[IV];
            // fMidpert[IV]        =   half * (fLpert[IV]+fR[IV]);

            Cal_GMRES_Visflux(fLpert, fR, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpL[m][IV] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            //! perturb the w of the left cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fLpert[index]   = fL[index];
                fMidpert[index] = fMid[index];  
            }
            fLpert[nEquation]   =   fL[nEquation];
            fMidpert[nEquation] =   half * (fLpert[nEquation]+fR[nEquation]);
            fLpert[IW]          +=  perturbScale;
            dfpert              =   fLpert[IW] - fL[IW];
            // fMidpert[IW]        =   half * (fLpert[IW]+fR[IW]);

            Cal_GMRES_Visflux(fLpert, fR, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpL[m][IW] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            //! perturb the pressure of the left cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fLpert[index]   = fL[index];
                fMidpert[index] = fMid[index];  
            }
            fLpert[IP]          *=  (1.0 + perturbScale);
            fLpert[nEquation]   =   fLpert[IP]/(coefficientofStateEquation * fLpert[IR]);
            dfpert              =   fLpert[IP] - fL[IP];
            fMidpert[IP]        =   half * (fLpert[IP]+fR[IP]);
            fMidpert[nEquation] =   half * (fLpert[nEquation]+fR[nEquation]);

            Cal_GMRES_Visflux(fLpert, fR, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpL[m][IP] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            //! perturb the density of the right cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fRpert[index]   = fR[index];
                fMidpert[index] = fMid[index];  
            }
            fRpert[IR]          *=  (1.0 + perturbScale);
            fRpert[nEquation]   =   fRpert[IP]/(coefficientofStateEquation * fRpert[IR]);
            dfpert              =   fRpert[IR] - fR[IR];
            fMidpert[IR]        =   half * (fL[IR]+fRpert[IR]);
            fMidpert[nEquation] =   half * (fL[nEquation]+fRpert[nEquation]);

            Cal_GMRES_Visflux(fL, fRpert, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpR[m][IR] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            //! perturb the u of the right cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fRpert[index]   = fR[index];
                fMidpert[index] = fMid[index];  
            }
            fRpert[nEquation]   =   fR[nEquation];
            fMidpert[nEquation] =   half * (fL[nEquation]+fRpert[nEquation]);
            fRpert[IU]          +=  perturbScale;
            dfpert              =   fRpert[IU] - fR[IU];
            // fMidpert[IU]        =   half * (fL[IU]+fRpert[IU]);

            Cal_GMRES_Visflux(fL, fRpert, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpR[m][IU] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            // perturb the v of the right cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fRpert[index]   = fR[index];
                fMidpert[index] = fMid[index];  
            }
            fRpert[nEquation]   =   fR[nEquation];
            fMidpert[nEquation] =   half * (fL[nEquation]+fRpert[nEquation]);
            fRpert[IV]          +=  perturbScale;
            dfpert              =   fRpert[IV] - fR[IV];
            // fMidpert[IV]        =   half * (fL[IV]+fRpert[IV]);

            Cal_GMRES_Visflux(fL, fRpert, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpR[m][IV] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            // perturb the w of the right cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fRpert[index]   = fR[index];
                fMidpert[index] = fMid[index];  
            }
            fRpert[nEquation]   =   fR[nEquation];
            fMidpert[nEquation] =   half * (fL[nEquation]+fRpert[nEquation]);
            fRpert[IW]          +=  perturbScale;
            dfpert              =   fRpert[IW] - fR[IW];
            // fMidpert[IW]        =   half * (fL[IW]+fRpert[IW]);

            Cal_GMRES_Visflux(fL, fRpert, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpR[m][IW] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            //! perturb the pressure of the right cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fRpert[index]   = fR[index];
                fMidpert[index] = fMid[index];  
            }
            fRpert[IP]          *=  (1.0 + perturbScale);
            fRpert[nEquation]   =   fRpert[IP]/(coefficientofStateEquation * fRpert[IR]);
            dfpert              =   fRpert[IP] - fR[IP];
            fMidpert[IP]        =   half * (fL[IP]+fRpert[IP]);
            fMidpert[nEquation] =   half * (fL[nEquation]+fRpert[nEquation]);

            Cal_GMRES_Visflux(fL, fRpert, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpR[m][IP] = (fluxp[m] - flux[m][jFace])/dfpert;
            }
        }
        else
        {
            //! perturb the density of the left cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fLpert[index]   = fL[index];
                fMidpert[index] = fMid[index];  
            }
            fLpert[IR]          *=  (1.0 + perturbScale);
            fLpert[nEquation]   =   fLpert[IP]/(coefficientofStateEquation * fLpert[IR]);
            dfpert              =   fLpert[IR] - fL[IR];
            fMidpert[IR]        =   half * (fLpert[IR]+fR[IR]);
            fMidpert[nEquation] =   half * (fLpert[nEquation]+fR[nEquation]);

            Cal_GMRES_Visflux(fLpert, fR, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpL[m][IR] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            //! perturb the u of the left cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fLpert[index]   = fL[index];
                fMidpert[index] = fMid[index];  
            }
            fLpert[nEquation]   =   fL[nEquation];
            fMidpert[nEquation] =   half * (fLpert[nEquation]+fR[nEquation]);
            fLpert[IU]          +=  perturbScale;
            dfpert              =   fLpert[IU] - fL[IU];
            fMidpert[IU]        =   half * (fLpert[IU]+fR[IU]);

            Cal_GMRES_Visflux(fLpert, fR, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpL[m][IU] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            //! perturb the v of the left cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fLpert[index]   = fL[index];
                fMidpert[index] = fMid[index];  
            }
            fLpert[nEquation]   =   fL[nEquation];
            fMidpert[nEquation] =   half * (fLpert[nEquation]+fR[nEquation]);
            fLpert[IV]          +=  perturbScale;
            dfpert              =   fLpert[IV] - fL[IV];
            fMidpert[IV]        =   half * (fLpert[IV]+fR[IV]);

            Cal_GMRES_Visflux(fLpert, fR, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpL[m][IV] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            //! perturb the w of the left cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fLpert[index]   = fL[index];
                fMidpert[index] = fMid[index];  
            }
            fLpert[nEquation]   =   fL[nEquation];
            fMidpert[nEquation] =   half * (fLpert[nEquation]+fR[nEquation]);
            fLpert[IW]          +=  perturbScale;
            dfpert              =   fLpert[IW] - fL[IW];
            fMidpert[IW]        =   half * (fLpert[IW]+fR[IW]);

            Cal_GMRES_Visflux(fLpert, fR, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpL[m][IW] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            //! perturb the pressure of the left cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fLpert[index]   = fL[index];
                fMidpert[index] = fMid[index];  
            }
            fLpert[IP]          *=  (1.0 + perturbScale);
            fLpert[nEquation]   =   fLpert[IP]/(coefficientofStateEquation * fLpert[IR]);
            dfpert              =   fLpert[IP] - fL[IP];
            fMidpert[IP]        =   half * (fLpert[IP]+fR[IP]);
            fMidpert[nEquation] =   half * (fLpert[nEquation]+fR[nEquation]);

            Cal_GMRES_Visflux(fLpert, fR, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpL[m][IP] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            //! perturb the density of the right cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fRpert[index]   = fR[index];
                fMidpert[index] = fMid[index];  
            }
            fRpert[IR]          *=  (1.0 + perturbScale);
            fRpert[nEquation]   =   fRpert[IP]/(coefficientofStateEquation * fRpert[IR]);
            dfpert              =   fRpert[IR] - fR[IR];
            fMidpert[IR]        =   half * (fL[IR]+fRpert[IR]);
            fMidpert[nEquation] =   half * (fL[nEquation]+fRpert[nEquation]);

            Cal_GMRES_Visflux(fL, fRpert, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpR[m][IR] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            //! perturb the u of the right cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fRpert[index]   = fR[index];
                fMidpert[index] = fMid[index];  
            }
            fRpert[nEquation]   =   fR[nEquation];
            fMidpert[nEquation] =   half * (fL[nEquation]+fRpert[nEquation]);
            fRpert[IU]          +=  perturbScale;
            dfpert              =   fRpert[IU] - fR[IU];
            fMidpert[IU]        =   half * (fL[IU]+fRpert[IU]);

            Cal_GMRES_Visflux(fL, fRpert, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpR[m][IU] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            //! perturb the v of the right cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fRpert[index]   = fR[index];
                fMidpert[index] = fMid[index];  
            }
            fRpert[nEquation]   =   fR[nEquation];
            fMidpert[nEquation] =   half * (fL[nEquation]+fRpert[nEquation]);
            fRpert[IV]          +=  perturbScale;
            dfpert              =   fRpert[IV] - fR[IV];
            fMidpert[IV]        =   half * (fL[IV]+fRpert[IV]);

            Cal_GMRES_Visflux(fL, fRpert, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpR[m][IV] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            //! perturb the w of the right cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fRpert[index]   = fR[index];
                fMidpert[index] = fMid[index];  
            }
            fRpert[nEquation]   =   fR[nEquation];
            fMidpert[nEquation] =   half * (fL[nEquation]+fRpert[nEquation]);
            fRpert[IW]          +=  perturbScale;
            dfpert              =   fRpert[IW] - fR[IW];
            fMidpert[IW]        =   half * (fL[IW]+fRpert[IW]);

            Cal_GMRES_Visflux(fL, fRpert, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpR[m][IW] = (fluxp[m] - flux[m][jFace])/dfpert;
            }

            //! perturb the pressure of the right cell of the face
            for(int index = 0; index < nEquation; index++)
            {
                fRpert[index]   = fR[index];
                fMidpert[index] = fMid[index];
            }
            fRpert[IP]          *= (1.0 + perturbScale);
            fRpert[nEquation]   =  fRpert[IP] / (coefficientofStateEquation * fRpert[IR]);
            dfpert              =  fRpert[IP] - fR[IP];
            fMidpert[IP]        =  half * (fL[IP] + fRpert[IP]);
            fMidpert[nEquation] =  half * (fL[nEquation] + fRpert[nEquation]);

            Cal_GMRES_Visflux(fL, fRpert, fMidpert, dqdxL, dqdyL, dqdzL, 
                              dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                              gridIn, faceProxy);

            for (int m = 0; m < nLaminar; ++ m)
            {
                dFluxdpR[m][IP] = (fluxp[m] - flux[m][jFace])/dfpert;
            }
        }

        //! assemble the matrix
        //! contributions from the bc
        int indexre = (re-nTotalCell)*nEquation;
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
            }
        }

        //! GMRESPV , obtain the dqdcv for the left cell
        gas->dPrimitive2dConservative(fL,gama[le],dqdcvL);

        //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
        //! then minus to the corresponding location in dRdq.
        int indexf = re * nEquation;    //! first index
        int indexs = le * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] -= dFluxdpL[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];
                }
            }
        }

        //! contributions from the bc, 
        //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dFluxdpL[indexI][indexJ]  += dFluxdpR[indexI][indexK]*
                                                dDdPlocal[indexK][indexJ];
                }
            }
        }

        //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = le * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] += dFluxdpL[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];
                }
            }
        }

        //! GMRESPV , obtain the dqdcv for the right cell
        gas->dPrimitive2dConservative(fR,gama[re],dqdcvR);

        //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
        //! then minus to the corresponding location in dRdq.
        indexf = re * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] -= dFluxdpR[indexI][indexK]*
                                                             dqdcvR[indexK][indexJ];
                }
            }
        }

        //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] += dFluxdpR[indexI][indexK]*
                                                             dqdcvR[indexK][indexJ];
                }
            }
        }

        //! considering the perturbation of the gradient in the flux calculation
        //! perturb the drhodx of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxL[index];
        }
        dgrad[nEquation] =  dqdxL[nEquation];
        dgrad[IR]        += perturbScale;
        dgradpert        =  dgrad[IR] - dqdxL[IR];

        Cal_GMRES_Visflux(fL, fR, fMid, dgrad, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxL[m][IR] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dudx of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxL[index];
        }
        dgrad[nEquation] =  dqdxL[nEquation];
        dgrad[IU]        += perturbScale;
        dgradpert        =  dgrad[IU] - dqdxL[IU];

        Cal_GMRES_Visflux(fL, fR, fMid, dgrad, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxL[m][IU] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dvdx of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxL[index];
        }
        dgrad[nEquation] =  dqdxL[nEquation];
        dgrad[IV]        += perturbScale;
        dgradpert        =  dgrad[IV] - dqdxL[IV];

        Cal_GMRES_Visflux(fL, fR, fMid, dgrad, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxL[m][IV] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dwdx of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxL[index];
        }
        dgrad[nEquation] =  dqdxL[nEquation];
        dgrad[IW]        += perturbScale;
        dgradpert        =  dgrad[IW] - dqdxL[IW];

        Cal_GMRES_Visflux(fL, fR, fMid, dgrad, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxL[m][IW] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dPdx of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxL[index];
        }
        dgrad[nEquation] =   dqdxL[nEquation];
        dgrad[IP]        +=  perturbScale;
        dgradpert        =   dgrad[IP] - dqdxL[IP];

        Cal_GMRES_Visflux(fL, fR, fMid, dgrad, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxL[m][IP] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dTdx of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxL[index];
        }
        dgrad[nEquation] =  dqdxL[nEquation];
        dgrad[nEquation] += perturbScale;
        dgradpert        =  dgrad[nEquation] - dqdxL[nEquation];

        Cal_GMRES_Visflux(fL, fR, fMid, dgrad, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxL[m][nEquation] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the drhody of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyL[index];
        }
        dgrad[nEquation] =  dqdyL[nEquation];
        dgrad[IR]        += perturbScale;
        dgradpert        =  dgrad[IR] - dqdyL[IR];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dgrad, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyL[m][IR] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dudy of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index]   = dqdyL[index];
        }
        dgrad[nEquation]    =   dqdyL[nEquation];
        dgrad[IU]           +=  perturbScale;
        dgradpert           =   dgrad[IU] - dqdyL[IU];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dgrad, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyL[m][IU] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dvdy of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyL[index];
        }
        dgrad[nEquation] =  dqdyL[nEquation];
        dgrad[IV]        += perturbScale;
        dgradpert        =  dgrad[IV] - dqdyL[IV];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dgrad, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyL[m][IV] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dwdy of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyL[index];
        }
        dgrad[nEquation] =  dqdyL[nEquation];
        dgrad[IW]        += perturbScale;
        dgradpert        =  dgrad[IW] - dqdyL[IW];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dgrad, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyL[m][IW] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dPdy of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyL[index];
        }
        dgrad[nEquation] =  dqdyL[nEquation];
        dgrad[IP]        += perturbScale;
        dgradpert        =  dgrad[IP] - dqdyL[IP];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dgrad, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyL[m][IP] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dTdy of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyL[index];
        }
        dgrad[nEquation] =  dqdyL[nEquation];
        dgrad[nEquation] += perturbScale;
        dgradpert        =  dgrad[nEquation] - dqdyL[nEquation];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dgrad, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyL[m][nEquation] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the drhodz of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzL[index];
        }
        dgrad[nEquation] =  dqdzL[nEquation];
        dgrad[IR]        += perturbScale;
        dgradpert        =  dgrad[IR] - dqdzL[IR];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dgrad, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzL[m][IR] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dudz of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzL[index];
        }
        dgrad[nEquation] =  dqdzL[nEquation];
        dgrad[IU]        += perturbScale;
        dgradpert        =  dgrad[IU] - dqdzL[IU];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dgrad, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzL[m][IU] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dvdz of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzL[index];
        }
        dgrad[nEquation] =  dqdzL[nEquation];
        dgrad[IV]        += perturbScale;
        dgradpert        =  dgrad[IV] - dqdzL[IV];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dgrad, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzL[m][IV] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dwdz of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzL[index];
        }
        dgrad[nEquation] =  dqdzL[nEquation];
        dgrad[IW]        += perturbScale;
        dgradpert        =  dgrad[IW] - dqdzL[IW];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dgrad, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzL[m][IW] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dPdz of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzL[index];
        }
        dgrad[nEquation] =  dqdzL[nEquation];
        dgrad[IP]        += perturbScale;
        dgradpert        =  dgrad[IP] - dqdzL[IP];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dgrad, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzL[m][IP] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dTdz of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index]   = dqdzL[index];
        }
        dgrad[nEquation]    =   dqdzL[nEquation];
        dgrad[nEquation]    +=  perturbScale;
        dgradpert           =   dgrad[nEquation] - dqdzL[nEquation];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dgrad, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzL[m][nEquation] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the drhodx of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index]   = dqdxR[index];
        }
        dgrad[nEquation]    =   dqdxR[nEquation];
        dgrad[IR]           +=  perturbScale;
        dgradpert           =   dgrad[IR] - dqdxR[IR];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dgrad, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxR[m][IR] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dudx of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index]   = dqdxR[index];
        }
        dgrad[nEquation]    =   dqdxR[nEquation];
        dgrad[IU]           +=  perturbScale;
        dgradpert           =   dgrad[IU] - dqdxR[IU];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dgrad, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxR[m][IU] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dvdx of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index]   = dqdxR[index];
        }
        dgrad[nEquation]    =   dqdxR[nEquation];
        dgrad[IV]           +=  perturbScale;
        dgradpert           =   dgrad[IV] - dqdxR[IV];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dgrad, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxR[m][IV] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dwdx of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxR[index];
        }
        dgrad[nEquation] =  dqdxR[nEquation];
        dgrad[IW]        += perturbScale;
        dgradpert        =  dgrad[IW] - dqdxR[IW];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dgrad, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxR[m][IW] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dPdx of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxR[index];
        }
        dgrad[nEquation] =  dqdxR[nEquation];
        dgrad[IP]        += perturbScale;
        dgradpert        =  dgrad[IP] - dqdxR[IP];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dgrad, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxR[m][IP] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dTdx of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxR[index];
        }
        dgrad[nEquation] =  dqdxR[nEquation];
        dgrad[nEquation] += perturbScale;
        dgradpert        =  dgrad[nEquation] - dqdxR[nEquation];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dgrad, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxR[m][nEquation] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //!  perturb the drhody of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyR[index];
        }
        dgrad[nEquation] =  dqdyR[nEquation];
        dgrad[IR]        += perturbScale;
        dgradpert        =  dgrad[IR] - dqdyR[IR];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dgrad, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyR[m][IR] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dudy of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyR[index];
        }
        dgrad[nEquation] =  dqdyR[nEquation];
        dgrad[IU]        += perturbScale;
        dgradpert        =  dgrad[IU] - dqdyR[IU];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dgrad, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyR[m][IU] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dvdy of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyR[index];
        }
        dgrad[nEquation] =  dqdyR[nEquation];
        dgrad[IV]        += perturbScale;
        dgradpert        =  dgrad[IV] - dqdyR[IV];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dgrad, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyR[m][IV] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dwdy of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyR[index];
        }
        dgrad[nEquation] =  dqdyR[nEquation];
        dgrad[IW]        += perturbScale;
        dgradpert        =  dgrad[IW] - dqdyR[IW];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dgrad, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyR[m][IW] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dPdy of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyR[index];
        }
        dgrad[nEquation] =   dqdyR[nEquation];
        dgrad[IP]        +=  perturbScale;
        dgradpert        =   dgrad[IP] - dqdyR[IP];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dgrad, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyR[m][IP] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dTdy of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyR[index];
        }
        dgrad[nEquation] =  dqdyR[nEquation];
        dgrad[nEquation] += perturbScale;
        dgradpert        =  dgrad[nEquation] - dqdyR[nEquation];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dgrad, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyR[m][nEquation] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the drhodz of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzR[index];
        }
        dgrad[nEquation] =  dqdzR[nEquation];
        dgrad[IR]        += perturbScale;
        dgradpert        =  dgrad[IR] - dqdzR[IR];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dgrad, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzR[m][IR] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //!S perturb the dudz of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzR[index];
        }
        dgrad[nEquation] =  dqdzR[nEquation];
        dgrad[IU]        += perturbScale;
        dgradpert        =  dgrad[IU] - dqdzR[IU];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dgrad, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzR[m][IU] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dvdz of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzR[index];
        }
        dgrad[nEquation] =  dqdzR[nEquation];
        dgrad[IV]        += perturbScale;
        dgradpert        =  dgrad[IV] - dqdzR[IV];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dgrad, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzR[m][IV] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dwdz of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzR[index];
        }
        dgrad[nEquation] =  dqdzR[nEquation];
        dgrad[IW]        += perturbScale;
        dgradpert        =  dgrad[IW] - dqdzR[IW];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dgrad, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzR[m][IW] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dPdz of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzR[index];
        }
        dgrad[nEquation] =  dqdzR[nEquation];
        dgrad[IP]        += perturbScale;
        dgradpert        =  dgrad[IP] - dqdzR[IP];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dgrad, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzR[m][IP] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dTdz of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzR[index];
        }
        dgrad[nEquation] =  dqdzR[nEquation];
        dgrad[nEquation] += perturbScale;
        dgradpert        =  dgrad[nEquation] - dqdzR[nEquation];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dgrad, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzR[m][nEquation] = (fluxp[m] - flux[m][jFace])/dgradpert;
        }

        RDouble nxs    = xfn[iFace];
        RDouble nys    = yfn[iFace];
        RDouble nzs    = zfn[iFace];
        RDouble ns     = area[iFace];
        RDouble volume = vol[re];
        // dgradtmp = dfluxdgradxR*dgradRdqLx + dfluxdgradyR*dgradRdqLy + dfluxdgradzR*dgraddqLz
        // dRdq[le][le] += dgradtmp*dqdcvL
        gas->dGradient2dPrimitive(fL,-1,dgraddqx, 'x', nxs, nys,nzs,ns,volume);
        gas->dGradient2dPrimitive(fL,-1,dgraddqy, 'y', nxs, nys,nzs,ns,volume);
        gas->dGradient2dPrimitive(fL,-1,dgraddqz, 'z', nxs, nys,nzs,ns,volume);

        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dGradtmp[indexI][indexJ] = 0;
                
                for(int indexK=0; indexK < (nEquation+1); indexK++)
                {
                    dGradtmp[indexI][indexJ] += dFluxdgradxR[indexI][indexK]*
                                                             dgraddqx[indexK][indexJ];
                    dGradtmp[indexI][indexJ] += dFluxdgradyR[indexI][indexK]*
                                                             dgraddqy[indexK][indexJ];
                    dGradtmp[indexI][indexJ] += dFluxdgradzR[indexI][indexK]*
                                                             dgraddqz[indexK][indexJ];
                }
            }
        }

        //! convert dGradtmp by right multiplying the dqdcvL,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = le * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] += dGradtmp[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];
                }
            }
        }

        // dgradtmp = dfluxdgradxR*dgradRdqLx + dfluxdgradyR*dgradRdqLy + dfluxdgradzR*dgraddqLz
        // dRdq[re][le] += dgradtmp*dqdcvL
        indexf = re * nEquation; // first index
        indexs = le * nEquation; // second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] -= dGradtmp[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];
                }
            }
        }

        nxs         = xfn[iFace];
        nys         = yfn[iFace];
        nzs         = zfn[iFace];
        ns          = area[iFace];
        volume      = vol[le];
        // dgradtmp = dfluxdgradxL*dgradLdqRx + dfluxdgradyL*dgradLdqRy + dfluxdgradzL*dgradLdqRz
        // dRdq[le][re] += dgradtmp*dqdcvL
        gas->dGradient2dPrimitive(fR,1,dgraddqx, 'x', nxs, nys,nzs,ns,volume);
        gas->dGradient2dPrimitive(fR,1,dgraddqy, 'y', nxs, nys,nzs,ns,volume);
        gas->dGradient2dPrimitive(fR,1,dgraddqz, 'z', nxs, nys,nzs,ns,volume);

        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dGradtmp[indexI][indexJ] = 0;
                
                for(int indexK=0; indexK < (nEquation+1); indexK++)
                {
                    dGradtmp[indexI][indexJ] += dFluxdgradxL[indexI][indexK]*
                                                             dgraddqx[indexK][indexJ];
                    dGradtmp[indexI][indexJ] += dFluxdgradyL[indexI][indexK]*
                                                             dgraddqy[indexK][indexJ];
                    dGradtmp[indexI][indexJ] += dFluxdgradzL[indexI][indexK]*
                                                             dgraddqz[indexK][indexJ];
                }
            }
        }

        //! convert dGradtmp by right multiplying the dqdcvL,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] += dGradtmp[indexI][indexK]*
                                                             dqdcvR[indexK][indexJ];
                }
            }
        }

        // dgradtmp = dfluxdgradxL*dgradLdqRx + dfluxdgradyL*dgradLdqRy + dfluxdgradzL*dgradLdqRz
        // dRdq[re][re] += dgradtmp*dqdcvL
        indexf = re * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] -= dGradtmp[indexI][indexK]*
                                                             dqdcvR[indexK][indexJ];
                }
            }
        }

        //! contributions from the bc, 
        //! dFluxdPR*dDdPlocal
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dGradtmpBC[indexI][indexJ] = 0.0;
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dGradtmpBC[indexI][indexJ]  += dGradtmp[indexI][indexK]*
                                                dDdPlocal[indexK][indexJ];
                }
            }
        }

        //! convert dGradtmpBC by right multiplying the dqdcvL,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = le * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] += dGradtmpBC[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];
                }
            }
        }

        //! the neighbors of the left cells
        for(int index = 0; index < neighborCells[le].size(); index++)
        {
            int neighborIndex = neighborCells[le][index];

            if(neighborIndex != re)
            {
                int Faceindex  = neighborFaces[le][index];
                int sign       = neighborLR[le][index];
                RDouble nxs    = xfn[Faceindex];
                RDouble nys    = yfn[Faceindex];
                RDouble nzs    = zfn[Faceindex];
                RDouble ns     = area[Faceindex];
                RDouble volume = vol[le];
                for(int m =0; m < nEquation; m++)
                {
                    fN[m] = primitiveVariable[m][neighborIndex];
                }

                gas->dPrimitive2dConservative(fN,gama[neighborIndex],dqdcvN);

                gas->dGradient2dPrimitive(fN,sign, dgraddqx, 'x', nxs, nys,nzs,ns,volume);
                gas->dGradient2dPrimitive(fN,sign, dgraddqy, 'y', nxs, nys,nzs,ns,volume);
                gas->dGradient2dPrimitive(fN,sign, dgraddqz, 'z', nxs, nys,nzs,ns,volume);

                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dGradtmp[indexI][indexJ] = 0;

                        for(int indexK=0; indexK < (nEquation+1); indexK++)
                        {
                            dGradtmp[indexI][indexJ] += dFluxdgradxL[indexI][indexK]*
                                                                     dgraddqx[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradyL[indexI][indexK]*
                                                                     dgraddqy[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradzL[indexI][indexK]*
                                                                     dgraddqz[indexK][indexJ];
                        }
                    }
                }

                int indexf  = le * nEquation; // first index
                int indexs  = neighborIndex * nEquation; // second index
                int indexf2 = re * nEquation; 
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {

                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            // dRL/dqLNeighbors += dRL/dGradL*dGradL/dqLNeighbors ===> dFluxdgradL*dgraddq
                            dRdq[indexf+indexI][indexs+indexJ] += dGradtmp[indexI][indexK]*
                                                                     dqdcvN[indexK][indexJ];
                            
                            // dRR/dqLNeighbors -= dRR/dGradL*dGradL/dqLNeighbors ===> dFluxdgradL*dgraddq
                            dRdq[indexf2+indexI][indexs+indexJ] -= dGradtmp[indexI][indexK]*
                                                                     dqdcvN[indexK][indexJ];
                        }
                    }
                }
            }
        }

        //! the neighbors of the right cells
        for(int index = 0; index < neighborCells[re].size(); index ++)
        {
            int neighborIndex = neighborCells[re][index];

            if(neighborIndex != le)
            {
                int Faceindex  = neighborFaces[re][index];
                int sign       = neighborLR[re][index];
                RDouble nxs    = xfn[Faceindex];
                RDouble nys    = yfn[Faceindex];
                RDouble nzs    = zfn[Faceindex];
                RDouble ns     = area[Faceindex];
                RDouble volume = vol[re];

                for(int m=0; m < nEquation; m++)
                {
                    fN[m] = primitiveVariable[m][neighborIndex];
                }

                gas->dPrimitive2dConservative(fN,gama[neighborIndex],dqdcvN);
    
                gas->dGradient2dPrimitive(fN,sign,dgraddqx, 'x', nxs, nys,nzs,ns,volume);
                gas->dGradient2dPrimitive(fN,sign,dgraddqy, 'y', nxs, nys,nzs,ns,volume);
                gas->dGradient2dPrimitive(fN,sign,dgraddqz, 'z', nxs, nys,nzs,ns,volume);

                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dGradtmp[indexI][indexJ] = 0;
                        for(int indexK=0; indexK < (nEquation+1); indexK++)
                        {
                            dGradtmp[indexI][indexJ] += dFluxdgradxR[indexI][indexK]*
                                                                     dgraddqx[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradyR[indexI][indexK]*
                                                                     dgraddqy[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradzR[indexI][indexK]*
                                                                     dgraddqz[indexK][indexJ];
                        }
                    }
                }
    
                int indexf  = le * nEquation;    //! first index
                int indexs  = neighborIndex * nEquation;    //! second index
                int indexf2 = re * nEquation; 
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                    
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            // dRL/dqRNeighbors += dRL/dGradR*dGradR/dqRNeighbors ===> dFluxdgradR*dgraddq
                            dRdq[indexf+indexI][indexs+indexJ] += dGradtmp[indexI][indexK]*
                                                                     dqdcvN[indexK][indexJ];

                            // dRR/dqRNeighbors -= dRR/dGradR*dGradR/dqRNeighbors ===> dFluxdgradR*dgraddq
                            dRdq[indexf2+indexI][indexs+indexJ] -= dGradtmp[indexI][indexK]*
                                                                     dqdcvN[indexK][indexJ];
                        }
                    }
                }
            }   
        }
    }

    //! interior faces
    for (int iFace = nMid; iFace < localEnd; ++ iFace)
    {
        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];
        int jFace  = iFace - localStart;

        for (int m = 0; m < nEquation; ++ m)
        {
            fL[m]    = primitiveVariable[m][le];
            fR[m]    = primitiveVariable[m][re];
            fMid[m]  = primitiveVariableFace[m][jFace];
            dqdxL[m] = dqdx[m][le];
            dqdyL[m] = dqdy[m][le];
            dqdzL[m] = dqdz[m][le];
            dqdxR[m] = dqdx[m][re];
            dqdyR[m] = dqdy[m][re];
            dqdzR[m] = dqdz[m][re];
        }
        fL[nEquation]       = t[ITT][le];
        fLpert[nEquation]   = t[ITT][le];
        fR[nEquation]       = t[ITT][re];
        fRpert[nEquation]   = t[ITT][re];
        fMid[nEquation]     = tm[ITT][jFace];
        fMidpert[nEquation] = tm[ITT][jFace];
        dqdxL[nEquation]    = dtdx[ITT][le];
        dqdyL[nEquation]    = dtdy[ITT][le];
        dqdzL[nEquation]    = dtdz[ITT][le];
        dqdxR[nEquation]    = dtdx[ITT][re];
        dqdyR[nEquation]    = dtdy[ITT][re];
        dqdzR[nEquation]    = dtdz[ITT][re];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            flux[m][jFace] = fluxp[m];
        }

        //! perturb the density of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            fLpert[index]   = fL[index];
            fMidpert[index] = fMid[index];  
        }
        fLpert[IR]          *= (1.0 + perturbScale);
        fLpert[nEquation]   =  fLpert[IP]/(coefficientofStateEquation * fLpert[IR]);
        dfpert              =  fLpert[IR] - fL[IR];
        fMidpert[IR]        =  half * (fLpert[IR]+fR[IR]);
        fMidpert[nEquation] =  half * (fLpert[nEquation]+fR[nEquation]);
        
        Cal_GMRES_Visflux(fLpert, fR, fMidpert, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdpL[m][IR] = (fluxp[m] - flux[m][jFace]) / dfpert;
        }

        //! perturb the u of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            fLpert[index]   = fL[index];
            fMidpert[index] = fMid[index];  
        }
        fLpert[nEquation]   =  fL[nEquation];
        fMidpert[nEquation] =  half * (fLpert[nEquation]+fR[nEquation]);
        fLpert[IU]          += perturbScale;
        dfpert              =  fLpert[IU] - fL[IU];
        fMidpert[IU]        =  half * (fLpert[IU]+fR[IU]);

        Cal_GMRES_Visflux(fLpert, fR, fMidpert, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdpL[m][IU] = (fluxp[m] - flux[m][jFace]) / dfpert;
        }

        //! perturb the v of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            fLpert[index]   = fL[index];
            fMidpert[index] = fMid[index];  
        }
        fLpert[nEquation]   =   fL[nEquation];
        fMidpert[nEquation] =   half * (fLpert[nEquation]+fR[nEquation]);
        fLpert[IV]          +=  perturbScale;
        dfpert              =   fLpert[IV] - fL[IV];
        fMidpert[IV]        =   half * (fLpert[IV]+fR[IV]);

        Cal_GMRES_Visflux(fLpert, fR, fMidpert, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdpL[m][IV] = (fluxp[m] - flux[m][jFace]) / dfpert;
        }

        //! perturb the w of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            fLpert[index]   = fL[index];
            fMidpert[index] = fMid[index];  
        }
        fLpert[nEquation]   =  fL[nEquation];
        fMidpert[nEquation] =  half * (fLpert[nEquation]+fR[nEquation]);
        fLpert[IW]          += perturbScale;
        dfpert              =  fLpert[IW] - fL[IW];
        fMidpert[IW]        =  half * (fLpert[IW]+fR[IW]);
        
        Cal_GMRES_Visflux(fLpert, fR, fMidpert, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdpL[m][IW] = (fluxp[m] - flux[m][jFace]) / dfpert;
        }

        //! perturb the pressure of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            fLpert[index]   = fL[index];
            fMidpert[index] = fMid[index];  
        }
        fLpert[IP]          *= (1.0 + perturbScale);
        fLpert[nEquation]   =  fLpert[IP]/(coefficientofStateEquation * fLpert[IR]);
        dfpert              =  fLpert[IP] - fL[IP];
        fMidpert[IP]        =  half * (fLpert[IP]+fR[IP]);
        fMidpert[nEquation] =  half * (fLpert[nEquation]+fR[nEquation]);

        Cal_GMRES_Visflux(fLpert, fR, fMidpert, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdpL[m][IP] = (fluxp[m] - flux[m][jFace]) / dfpert;
        }

        //! perturb the density of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            fRpert[index]   = fR[index];
            fMidpert[index] = fMid[index];  
        }
        fRpert[IR]          *= (1.0 + perturbScale);
        fRpert[nEquation]   =  fRpert[IP]/(coefficientofStateEquation * fRpert[IR]);
        dfpert              =  fRpert[IR] - fR[IR];
        fMidpert[IR]        =  half * (fL[IR]+fRpert[IR]);
        fMidpert[nEquation] =  half * (fL[nEquation]+fRpert[nEquation]);
        
        Cal_GMRES_Visflux(fL, fRpert, fMidpert, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdpR[m][IR] = (fluxp[m] - flux[m][jFace]) / dfpert;
        }

        //! perturb the u of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            fRpert[index]   = fR[index];
            fMidpert[index] = fMid[index];  
        }
        fRpert[nEquation]   =  fR[nEquation];
        fMidpert[nEquation] =  half * (fL[nEquation]+fRpert[nEquation]);
        fRpert[IU]          += perturbScale;
        dfpert              =  fRpert[IU] - fR[IU];
        fMidpert[IU]        =  half * (fL[IU]+fRpert[IU]);
        
        Cal_GMRES_Visflux(fL, fRpert, fMidpert, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdpR[m][IU] = (fluxp[m] - flux[m][jFace]) / dfpert;
        }

        //! perturb the v of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            fRpert[index]   = fR[index];
            fMidpert[index] = fMid[index];  
        }
        fRpert[nEquation]   =  fR[nEquation];
        fMidpert[nEquation] =  half * (fL[nEquation]+fRpert[nEquation]);
        fRpert[IV]          += perturbScale;
        dfpert              =  fRpert[IV] - fR[IV];
        fMidpert[IV]        =  half * (fL[IV]+fRpert[IV]);
        
        Cal_GMRES_Visflux(fL, fRpert, fMidpert, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdpR[m][IV] = (fluxp[m] - flux[m][jFace])/dfpert;
        }

        //! perturb the w of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            fRpert[index]   = fR[index];
            fMidpert[index] = fMid[index];  
        }
        fRpert[nEquation]   =  fR[nEquation];
        fMidpert[nEquation] =  half * (fL[nEquation]+fRpert[nEquation]);
        fRpert[IW]          += perturbScale;
        dfpert              =  fRpert[IW] - fR[IW];
        fMidpert[IW]        =  half * (fL[IW]+fRpert[IW]);

        Cal_GMRES_Visflux(fL, fRpert, fMidpert, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdpR[m][IW] = (fluxp[m] - flux[m][jFace]) / dfpert;
        }

        //! perturb the pressure of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            fRpert[index]   = fR[index];
            fMidpert[index] = fMid[index];  
        }
        fRpert[IP]          *= (1.0 + perturbScale);
        fRpert[nEquation]   =  fRpert[IP] / (coefficientofStateEquation * fRpert[IR]);
        dfpert              =  fRpert[IP] - fR[IP];
        fMidpert[IP]        =  half * (fL[IP]+fRpert[IP]);
        fMidpert[nEquation] =  half * (fL[nEquation]+fRpert[nEquation]);

        Cal_GMRES_Visflux(fL, fRpert, fMidpert, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdpR[m][IP] = (fluxp[m] - flux[m][jFace])/dfpert;
        }

        //! GMRESPV , obtain the dqdcv for the left cell
        gas->dPrimitive2dConservative(fL,gama[le],dqdcvL);

        //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
        //! then minus to the corresponding location in dRdq.
        int indexf = re * nEquation;    //! first index
        int indexs = le * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] -= dFluxdpL[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];
                }
            }
        }

        //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = le * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] += dFluxdpL[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];
                }
            }
        }

         //! GMRESPV , obtain the dqdcv for the right cell
        gas->dPrimitive2dConservative(fR,gama[re],dqdcvR);

        //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
        //! then minus to the corresponding location in dRdq.
        indexf = re * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] -= dFluxdpR[indexI][indexK]*
                                                             dqdcvR[indexK][indexJ];
                }
            }
        }

        //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation; // first index
        indexs = re * nEquation; // second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] += dFluxdpR[indexI][indexK]*
                                                             dqdcvR[indexK][indexJ];
                }
            }
        }

        //! considering the perturbation of the gradient in the flux calculation
        //! perturb the drhodx of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxL[index];
        }
        dgrad[nEquation] =  dqdxL[nEquation];
        dgrad[IR]        += perturbScale;
        dgradpert        =  dgrad[IR] - dqdxL[IR];

        Cal_GMRES_Visflux(fL, fR, fMid, dgrad, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxL[m][IR] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dudx of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxL[index];
        }
        dgrad[nEquation] =  dqdxL[nEquation];
        dgrad[IU]        += perturbScale;
        dgradpert        =  dgrad[IU] - dqdxL[IU];

        Cal_GMRES_Visflux(fL, fR, fMid, dgrad, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxL[m][IU] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dvdx of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxL[index];
        }
        dgrad[nEquation] =  dqdxL[nEquation];
        dgrad[IV]        += perturbScale;
        dgradpert        =  dgrad[IV] - dqdxL[IV];

        Cal_GMRES_Visflux(fL, fR, fMid, dgrad, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxL[m][IV] = (fluxp[m] - flux[m][jFace])/dgradpert;
        }

        //! perturb the dwdx of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxL[index];
        }
        dgrad[nEquation] =  dqdxL[nEquation];
        dgrad[IW]        += perturbScale;
        dgradpert        =  dgrad[IW] - dqdxL[IW];

        Cal_GMRES_Visflux(fL, fR, fMid, dgrad, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxL[m][IW] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dPdx of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxL[index];
        }
        dgrad[nEquation] =  dqdxL[nEquation];
        dgrad[IP]        += perturbScale;
        dgradpert        =  dgrad[IP] - dqdxL[IP];

        Cal_GMRES_Visflux(fL, fR, fMid, dgrad, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxL[m][IP] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dTdx of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxL[index];
        }
        dgrad[nEquation] =  dqdxL[nEquation];
        dgrad[nEquation] += perturbScale;
        dgradpert        =  dgrad[nEquation] - dqdxL[nEquation];

        Cal_GMRES_Visflux(fL, fR, fMid, dgrad, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxL[m][nEquation] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the drhody of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyL[index];
        }
        dgrad[nEquation] =  dqdyL[nEquation];
        dgrad[IR]        += perturbScale;
        dgradpert        =  dgrad[IR] - dqdyL[IR];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dgrad, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyL[m][IR] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dudy of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyL[index];
        }
        dgrad[nEquation] =  dqdyL[nEquation];
        dgrad[IU]        += perturbScale;
        dgradpert        =  dgrad[IU] - dqdyL[IU];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dgrad, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyL[m][IU] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dvdy of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyL[index];
        }
        dgrad[nEquation] =  dqdyL[nEquation];
        dgrad[IV]        += perturbScale;
        dgradpert        =  dgrad[IV] - dqdyL[IV];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dgrad, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyL[m][IV] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dwdy of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyL[index];
        }
        dgrad[nEquation] =  dqdyL[nEquation];
        dgrad[IW]        += perturbScale;
        dgradpert        =  dgrad[IW] - dqdyL[IW];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dgrad, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyL[m][IW] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dPdy of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyL[index];
        }
        dgrad[nEquation] =  dqdyL[nEquation];
        dgrad[IP]        += perturbScale;
        dgradpert        =  dgrad[IP] - dqdyL[IP];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dgrad, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyL[m][IP] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dTdy of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyL[index];
        }
        dgrad[nEquation] =  dqdyL[nEquation];
        dgrad[nEquation] += perturbScale;
        dgradpert        =  dgrad[nEquation] - dqdyL[nEquation];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dgrad, dqdzL, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyL[m][nEquation] = (fluxp[m] - flux[m][jFace])/dgradpert;
        }

        //! perturb the drhodz of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzL[index];
        }
        dgrad[nEquation] =  dqdzL[nEquation];
        dgrad[IR]        += perturbScale;
        dgradpert        =  dgrad[IR] - dqdzL[IR];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dgrad, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzL[m][IR] = (fluxp[m] - flux[m][jFace])/dgradpert;
        }

        //! perturb the dudz of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index]   = dqdzL[index];
        }
        dgrad[nEquation]    =   dqdzL[nEquation];
        dgrad[IU]           +=  perturbScale;
        dgradpert           =   dgrad[IU] - dqdzL[IU];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dgrad, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzL[m][IU] = (fluxp[m] - flux[m][jFace])/dgradpert;
        }

        //! perturb the dvdz of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzL[index];
        }
        dgrad[nEquation] =  dqdzL[nEquation];
        dgrad[IV]        += perturbScale;
        dgradpert        =  dgrad[IV] - dqdzL[IV];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dgrad, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzL[m][IV] = (fluxp[m] - flux[m][jFace])/dgradpert;
        }

        //! perturb the dwdz of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzL[index];
        }
        dgrad[nEquation] =  dqdzL[nEquation];
        dgrad[IW]        += perturbScale;
        dgradpert        =  dgrad[IW] - dqdzL[IW];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dgrad, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzL[m][IW] = (fluxp[m] - flux[m][jFace])/dgradpert;
        }

        //! perturb the dPdz of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzL[index];
        }
        dgrad[nEquation] =  dqdzL[nEquation];
        dgrad[IP]        += perturbScale;
        dgradpert        =  dgrad[IP] - dqdzL[IP];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dgrad, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzL[m][IP] = (fluxp[m] - flux[m][jFace])/dgradpert;
        }

        //! perturb the dTdz of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzL[index];
        }
        dgrad[nEquation] =  dqdzL[nEquation];
        dgrad[nEquation] += perturbScale;
        dgradpert        =  dgrad[nEquation] - dqdzL[nEquation];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dgrad, 
                          dqdxR, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzL[m][nEquation] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the drhodx of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxR[index];
        }
        dgrad[nEquation] =  dqdxR[nEquation];
        dgrad[IR]        += perturbScale;
        dgradpert        =  dgrad[IR] - dqdxR[IR];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dgrad, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxR[m][IR] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dudx of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxR[index];
        }
        dgrad[nEquation] =  dqdxR[nEquation];
        dgrad[IU]        += perturbScale;
        dgradpert        =  dgrad[IU] - dqdxR[IU];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dgrad, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxR[m][IU] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dvdx of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index]   = dqdxR[index];
        }
        dgrad[nEquation]    =   dqdxR[nEquation];
        dgrad[IV]           +=  perturbScale;
        dgradpert           =   dgrad[IV] - dqdxR[IV];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dgrad, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxR[m][IV] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dwdx of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxR[index];
        }
        dgrad[nEquation] =  dqdxR[nEquation];
        dgrad[IW]        += perturbScale;
        dgradpert        =  dgrad[IW] - dqdxR[IW];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dgrad, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxR[m][IW] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dPdx of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxR[index];
        }
        dgrad[nEquation] =  dqdxR[nEquation];
        dgrad[IP]        += perturbScale;
        dgradpert        =  dgrad[IP] - dqdxR[IP];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dgrad, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxR[m][IP] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dTdx of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdxR[index];
        }
        dgrad[nEquation] =  dqdxR[nEquation];
        dgrad[nEquation] += perturbScale;
        dgradpert        =  dgrad[nEquation] - dqdxR[nEquation];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dgrad, dqdyR, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradxR[m][nEquation] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the drhody of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyR[index];
        }
        dgrad[nEquation] =  dqdyR[nEquation];
        dgrad[IR]        += perturbScale;
        dgradpert        =  dgrad[IR] - dqdyR[IR];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dgrad, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyR[m][IR] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dudy of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyR[index];
        }
        dgrad[nEquation] =  dqdyR[nEquation];
        dgrad[IU]        += perturbScale;
        dgradpert        =  dgrad[IU] - dqdyR[IU];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dgrad, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyR[m][IU] = (fluxp[m] - flux[m][jFace])/dgradpert;
        }

        //! perturb the dvdy of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyR[index];
        }
        dgrad[nEquation] =  dqdyR[nEquation];
        dgrad[IV]        += perturbScale;
        dgradpert        =  dgrad[IV] - dqdyR[IV];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dgrad, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyR[m][IV] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dwdy of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyR[index];
        }
        dgrad[nEquation] =  dqdyR[nEquation];
        dgrad[IW]        += perturbScale;
        dgradpert        =  dgrad[IW] - dqdyR[IW];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dgrad, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyR[m][IW] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dPdy of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyR[index];
        }
        dgrad[nEquation] =  dqdyR[nEquation];
        dgrad[IP]        += perturbScale;
        dgradpert        =  dgrad[IP] - dqdyR[IP];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dgrad, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyR[m][IP] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dTdy of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdyR[index];
        }
        dgrad[nEquation] =  dqdyR[nEquation];
        dgrad[nEquation] += perturbScale;
        dgradpert        =  dgrad[nEquation] - dqdyR[nEquation];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dgrad, dqdzR, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradyR[m][nEquation] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the drhodz of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzR[index];
        }
        dgrad[nEquation] =  dqdzR[nEquation];
        dgrad[IR]        += perturbScale;
        dgradpert        =  dgrad[IR] - dqdzR[IR];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dgrad, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzR[m][IR] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dudz of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzR[index];
        }
        dgrad[nEquation] =  dqdzR[nEquation];
        dgrad[IU]        += perturbScale;
        dgradpert        =  dgrad[IU] - dqdzR[IU];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dgrad, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzR[m][IU] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dvdz of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzR[index];
        }
        dgrad[nEquation] =  dqdzR[nEquation];
        dgrad[IV]        += perturbScale;
        dgradpert        =  dgrad[IV] - dqdzR[IV];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dgrad, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzR[m][IV] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dwdz of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzR[index];
        }
        dgrad[nEquation] =  dqdzR[nEquation];
        dgrad[IW]        += perturbScale;
        dgradpert        =  dgrad[IW] - dqdzR[IW];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dgrad, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzR[m][IW] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dPdz of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzR[index];
        }
        dgrad[nEquation] =  dqdzR[nEquation];
        dgrad[IP]        += perturbScale;
        dgradpert        =  dgrad[IP] - dqdzR[IP];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dgrad, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzR[m][IP] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        //! perturb the dTdz of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            dgrad[index] = dqdzR[index];
        }
        dgrad[nEquation] =  dqdzR[nEquation];
        dgrad[nEquation] += perturbScale;
        dgradpert        =  dgrad[nEquation] - dqdzR[nEquation];

        Cal_GMRES_Visflux(fL, fR, fMid, dqdxL, dqdyL, dqdzL, 
                          dqdxR, dqdyR, dgrad, iFace, jFace, fluxp, 
                          gridIn, faceProxy);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dFluxdgradzR[m][nEquation] = (fluxp[m] - flux[m][jFace]) / dgradpert;
        }

        RDouble nxs    = xfn[iFace];
        RDouble nys    = yfn[iFace];
        RDouble nzs    = zfn[iFace];
        RDouble ns     = area[iFace];
        RDouble volume = vol[re];
        // dgradtmp = dfluxdgradxR*dgradRdqLx + dfluxdgradyR*dgradRdqLy + dfluxdgradzR*dgraddqLz
        // dRdq[le][le] += dgradtmp*dqdcvL
        gas->dGradient2dPrimitive(fL,-1,dgraddqx, 'x', nxs, nys,nzs,ns,volume);
        gas->dGradient2dPrimitive(fL,-1,dgraddqy, 'y', nxs, nys,nzs,ns,volume);
        gas->dGradient2dPrimitive(fL,-1,dgraddqz, 'z', nxs, nys,nzs,ns,volume);

        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dGradtmp[indexI][indexJ] = 0;

                for(int indexK=0; indexK < (nEquation+1); indexK++)
                {
                    dGradtmp[indexI][indexJ] += dFluxdgradxR[indexI][indexK]*
                                                             dgraddqx[indexK][indexJ];
                    dGradtmp[indexI][indexJ] += dFluxdgradyR[indexI][indexK]*
                                                             dgraddqy[indexK][indexJ];
                    dGradtmp[indexI][indexJ] += dFluxdgradzR[indexI][indexK]*
                                                             dgraddqz[indexK][indexJ];
                }
            }
        }

        //! convert dGradtmp by right multiplying the dqdcvL,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = le * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] += dGradtmp[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];
                }
            }
        }

        // dgradtmp = dfluxdgradxR*dgradRdqLx + dfluxdgradyR*dgradRdqLy + dfluxdgradzR*dgraddqLz
        // dRdq[re][le] += dgradtmp*dqdcvL
        indexf = re * nEquation;    //! first index
        indexs = le * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] -= dGradtmp[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];
                }
            }
        }


        nxs    = xfn[iFace];
        nys    = yfn[iFace];
        nzs    = zfn[iFace];
        ns     = area[iFace];
        volume = vol[le];
        // dgradtmp = dfluxdgradxL*dgradLdqRx + dfluxdgradyL*dgradLdqRy + dfluxdgradzL*dgradLdqRz
        // dRdq[le][re] += dgradtmp*dqdcvL
        gas->dGradient2dPrimitive(fR,1,dgraddqx, 'x', nxs, nys,nzs,ns,volume);
        gas->dGradient2dPrimitive(fR,1,dgraddqy, 'y', nxs, nys,nzs,ns,volume);
        gas->dGradient2dPrimitive(fR,1,dgraddqz, 'z', nxs, nys,nzs,ns,volume);

        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dGradtmp[indexI][indexJ] = 0;

                for(int indexK=0; indexK < (nEquation+1); indexK++)
                {
                    dGradtmp[indexI][indexJ] += dFluxdgradxL[indexI][indexK]*
                                                             dgraddqx[indexK][indexJ];
                    dGradtmp[indexI][indexJ] += dFluxdgradyL[indexI][indexK]*
                                                             dgraddqy[indexK][indexJ];
                    dGradtmp[indexI][indexJ] += dFluxdgradzL[indexI][indexK]*
                                                             dgraddqz[indexK][indexJ];
                }
            }
        }

        //! convert dGradtmp by right multiplying the dqdcvL,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] += dGradtmp[indexI][indexK]*
                                                             dqdcvR[indexK][indexJ];
                }
            }
        }

        // dgradtmp = dfluxdgradxL*dgradLdqRx + dfluxdgradyL*dgradLdqRy + dfluxdgradzL*dgradLdqRz
        // dRdq[re][re] += dgradtmp*dqdcvL
        indexf = re * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] -= dGradtmp[indexI][indexK]*
                                                             dqdcvR[indexK][indexJ];
                }
            }
        }

        //! the neighbors of the left cells
        for(int index = 0; index < neighborCells[le].size(); index++)
        {
            int neighborIndex = neighborCells[le][index];

            if(neighborIndex != re)
            {
                int Faceindex  = neighborFaces[le][index];
                int sign       = neighborLR[le][index];
                RDouble nxs    = xfn[Faceindex];
                RDouble nys    = yfn[Faceindex];
                RDouble nzs    = zfn[Faceindex];
                RDouble ns     = area[Faceindex];
                RDouble volume = vol[le];

                for(int m =0; m < nEquation; m++)
                {
                    fN[m] = primitiveVariable[m][neighborIndex];
                }

                gas->dPrimitive2dConservative(fN,gama[neighborIndex],dqdcvN);

                gas->dGradient2dPrimitive(fN,sign,dgraddqx, 'x', nxs, nys,nzs,ns,volume);
                gas->dGradient2dPrimitive(fN,sign,dgraddqy, 'y', nxs, nys,nzs,ns,volume);
                gas->dGradient2dPrimitive(fN,sign,dgraddqz, 'z', nxs, nys,nzs,ns,volume);

                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dGradtmp[indexI][indexJ] = 0;

                        for(int indexK=0; indexK < (nEquation+1); indexK++)
                        {
                            dGradtmp[indexI][indexJ] += dFluxdgradxL[indexI][indexK]*
                                                                     dgraddqx[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradyL[indexI][indexK]*
                                                                     dgraddqy[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradzL[indexI][indexK]*
                                                                     dgraddqz[indexK][indexJ];
                        }
                    }
                }

                int indexf  = le * nEquation;    //! first index
                int indexs  = neighborIndex * nEquation;    //! second index
                int indexf2 = re * nEquation; 
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            // dRL/dqLNeighbors += dRL/dGradL*dGradL/dqLNeighbors ===> dFluxdgradL*dgraddq
                            dRdq[indexf+indexI][indexs+indexJ] += dGradtmp[indexI][indexK]*
                                                                     dqdcvN[indexK][indexJ];
                            
                            // dRR/dqLNeighbors -= dRR/dGradL*dGradL/dqLNeighbors ===> dFluxdgradL*dgraddq
                            dRdq[indexf2+indexI][indexs+indexJ] -= dGradtmp[indexI][indexK]*
                                                                     dqdcvN[indexK][indexJ];
                        }
                    }
                }
            }
        }

        //! the neighbors of the right cells
        for(int index = 0; index < neighborCells[re].size(); index++)
        {
            int neighborIndex = neighborCells[re][index];

            if(neighborIndex != le)
            {
                int Faceindex       = neighborFaces[re][index];
                int sign            = neighborLR[re][index];
                RDouble nxs         = xfn[Faceindex];
                RDouble nys         = yfn[Faceindex];
                RDouble nzs         = zfn[Faceindex];
                RDouble ns          = area[Faceindex];
                RDouble volume      = vol[re];

                for(int m=0; m < nEquation; m++)
                {
                    fN[m] = primitiveVariable[m][neighborIndex];
                }

                gas->dPrimitive2dConservative(fN,gama[neighborIndex],dqdcvN);

                gas->dGradient2dPrimitive(fN,sign,dgraddqx, 'x', nxs, nys,nzs,ns,volume);
                gas->dGradient2dPrimitive(fN,sign,dgraddqy, 'y', nxs, nys,nzs,ns,volume);
                gas->dGradient2dPrimitive(fN,sign,dgraddqz, 'z', nxs, nys,nzs,ns,volume);

                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dGradtmp[indexI][indexJ] = 0;
                        for(int indexK=0; indexK < (nEquation+1); indexK++)
                        {
                            dGradtmp[indexI][indexJ] += dFluxdgradxR[indexI][indexK]*
                                                                     dgraddqx[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradyR[indexI][indexK]*
                                                                     dgraddqy[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradzR[indexI][indexK]*
                                                                     dgraddqz[indexK][indexJ];
                        }
                    }
                }

                int indexf  = le * nEquation;    //! first index
                int indexs  = neighborIndex * nEquation;    //! second index
                int indexf2 = re * nEquation; 
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            // dRL/dqRNeighbors += dRL/dGradR*dGradR/dqRNeighbors ===> dFluxdgradR*dgraddq
                            dRdq[indexf+indexI][indexs+indexJ] += dGradtmp[indexI][indexK]*
                                                                     dqdcvN[indexK][indexJ];

                            // dRR/dqRNeighbors -= dRR/dGradR*dGradR/dqRNeighbors ===> dFluxdgradR*dgraddq
                            dRdq[indexf2+indexI][indexs+indexJ] -= dGradtmp[indexI][indexK]*
                                                                     dqdcvN[indexK][indexJ];
                        }
                    }
                }
            }
        }

        // output convert matrix
        // std::cout << "===================dgraddqx===========================\n";
        // for(int indexI=0; indexI < nEquation+1; indexI++)
        // {
        //     
        //     for(int indexJ =0; indexJ < nEquation; indexJ ++)
        //     {
        //        std::cout << dgraddqx[indexI][indexJ] << " ";
        //     }
        //     std::cout << "\n";
        // }
        // std::cout << "===================dgraddqy===========================\n";
        // for(int indexI=0; indexI < nEquation+1; indexI++)
        // {
        //     
        //     for(int indexJ =0; indexJ < nEquation; indexJ ++)
        //     {
        //        std::cout << dgraddqy[indexI][indexJ] << " ";
        //     }
        //     std::cout << "\n";
        // } 
        // std::cout << "===================dgraddqz===========================\n";
        // for(int indexI=0; indexI < nEquation+1; indexI++)
        // {
        //     
        //     for(int indexJ =0; indexJ < nEquation; indexJ ++)
        //     {
        //        std::cout << dgraddqz[indexI][indexJ] << " ";
        //     }
        //     std::cout << "\n";
        // }

        // std::cout << "the left cell is: " << le << " "<< "the neighboring cells are : " ;
        // for(int index = 0; index < neighborCells[le].size(); index++)
        // {
        //     std::cout << neighborCells[le][index] << " ";
        // }
        // std::cout << "\n";

        // std::cout << "the left cell2 is: " << le << " "<< "the neighboring cells are : " ;
        // for(auto iter = neighborCells[le].begin(); iter != neighborCells[le].end(); iter++)
        // {
        //     std::cout << *iter << " ";
        // }
        // std::cout << "\n";
  
        // std::cout << "the right cell is: " << re << " "<< "the neighboring cells are : " ;
        // for(int index = 0; index < neighborCells[re].size(); index++)
        // {
        //     std::cout << neighborCells[re][index] << " ";
        // }
        // std::cout << "\n";

    }

    delete [] dqdxL;    dqdxL = nullptr;
    delete [] dqdyL;    dqdyL = nullptr;
    delete [] dqdzL;    dqdzL = nullptr;

    delete [] dqdxR;    dqdxR = nullptr;
    delete [] dqdyR;    dqdyR = nullptr;
    delete [] dqdzR;    dqdzR = nullptr;
    delete [] dgrad;    dgrad = nullptr;

    delete [] fL;       fL = nullptr;
    delete [] fR;       fR = nullptr;
    delete [] fN;       fN = nullptr;
    delete [] fMid;     fMid = nullptr;
    delete [] fluxp;    fluxp = nullptr;

    delete [] fLpert;   fLpert = nullptr;
    delete [] fRpert;   fRpert = nullptr;
    delete [] fMidpert; fMidpert = nullptr;

    for(int index = 0; index < nEquation; index++)
    {
        delete [] dqdcvL[index];
        delete [] dqdcvR[index];
        delete [] dqdcvN[index];
        delete [] dgraddqx[index];
        delete [] dgraddqy[index];
        delete [] dgraddqz[index];
    }
    delete [] dgraddqx[nEquation];
    delete [] dgraddqy[nEquation];
    delete [] dgraddqz[nEquation];
    delete [] dqdcvL;    dqdcvL = nullptr;
    delete [] dqdcvR;    dqdcvR = nullptr;
    delete [] dqdcvN;    dqdcvN = nullptr;
    delete [] dgraddqx;    dgraddqx = nullptr;
    delete [] dgraddqy;    dgraddqy = nullptr;
    delete [] dgraddqz;    dgraddqz = nullptr;
}

//! GMRESCoupled
template <typename T>
void NSSolverUnstruct::Cal_GMRES_Visflux_AD_Coupled(T* fL, T* fR, T fLturb, T fRturb, T* dqdxL,T* dqdyL, T* dqdzL, 
                                                    T* dqdxR, T* dqdyR, T* dqdzR, int iFace,
                                                    int jFace, T* flux, Grid *gridIn, 
                                                    FaceProxy *faceProxy, Param_NSSolverUnstruct *parameters)
{
        UnstructGrid *grid = UnstructGridCast(gridIn);
        int nBoundFace = grid->GetNBoundFace();

        int *leftCellofFace  = grid->GetLeftCellOfFace();
        int *rightCellofFace = grid->GetRightCellOfFace();
        int nTotalCell       = grid->GetNTotalCell();

        RDouble **dDdP_turb = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP_turb"));

        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

        RDouble *xfn  = grid->GetFaceNormalX();
        RDouble *yfn  = grid->GetFaceNormalY();
        RDouble *zfn  = grid->GetFaceNormalZ();
        RDouble *area = grid->GetFaceArea();
    
        RDouble *xcc  = grid->GetCellCenterX();
        RDouble *ycc  = grid->GetCellCenterY();
        RDouble *zcc  = grid->GetCellCenterZ();

        RDouble *xfc  = grid->GetFaceCenterX();
        RDouble *yfc  = grid->GetFaceCenterY();
        RDouble *zfc  = grid->GetFaceCenterZ();

        // Param_NSSolverUnstruct *parameters = GetControlParameters();
        int numberOfSpecies = parameters->GetNumberOfSpecies();

        int nLaminar  = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();

        int nEquation = nLaminar;

        int nm  = parameters->GetNSEquationNumber();
        RDouble skewnessAngle = parameters->GetSkewnessAngle();
        bool isFineGrid = grid->IsFinestGrid();
        RDouble *deltaL = faceProxy->GetWeightL();
        RDouble *deltaR = faceProxy->GetWeightR();
        NSFaceValue *faceVariable = faceProxy->GetNSFaceValue();
        // RDouble *kCp     = faceVariable->GetKCP();
        // RDouble *viscousLaminarFace   = faceVariable->GetMUL();
        // RDouble *viscousTurbulenceFace = faceVariable->GetMUT();

        RDouble refReNumber  = parameters->GetRefReNumber();
        RDouble oRefReNumber = parameters->GetoRefReNumber();
        RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();

        int nolstress = GlobalDataBase::GetIntParaFromDB("nolstress");
        int nrokplus  = GlobalDataBase::GetIntParaFromDB("nrokplus");

        RDouble **qTurb = 0;
        RDouble ** aniss  = 0;

        if (nrokplus > 0)
        {
            qTurb = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
            aniss  = reinterpret_cast<RDouble **> (grid->GetDataPtr("aniss"));
        }

        T *fvis = new T[nLaminar];
        T *dfdx = new T[nEquation];
        T *dfdy = new T[nEquation];
        T *dfdz = new T[nEquation];

        T *dfd1  = new T[nEquation];
        T *dfd2  = new T[nEquation];
        T *dfdn  = new T[nEquation];
        T *dfdt1 = new T[nEquation];
        T *dfdt2 = new T[nEquation];

        T *f1    = new T[nEquation];
        T *f2    = new T[nEquation];
        T *fMid  = new T[nEquation];
        T *primitiveVariableFace = new T[nEquation+nTemperatureModel];

        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];

        RDouble nxs = xfn[iFace];
        RDouble nys = yfn[iFace];
        RDouble nzs = zfn[iFace];

        RDouble t1x, t1y, t1z, t2x, t2y, t2z;
        T txx, tyy, tzz;
        T txy, txz, tyz;

        T viscousLaminarMin,viscousLaminarL,viscousLaminarR,viscousLaminarFace;
        T viscousTurbulenceL, viscousTurbulenceR, viscousTurbulenceFace,ld,ld3,fv1;
        GlobalDataBase::GetData("visl_min", &viscousLaminarMin, PHDOUBLE, 1);

        RDouble refGama = 1.4;
        RDouble referenceMachNumber         = parameters->GetRefMachNumber();
        RDouble referenceMachNumberSquare   = referenceMachNumber * referenceMachNumber;
        RDouble oPrandtlLaminar             = parameters->GetoPrandtlLaminar();
        RDouble oPrandtlTurbulence          = parameters->GetoPrandtlTurbulence();
        T kCpFace;
        RDouble coefficientofStateEquation = gas->GetCoefficientOfStateEquation();

        //! Get first tangential vector on the face.
        if (abs(nxs) > SMALL)
        {
            t1x =   nys;
            t1y = - nxs;
            t1z =   0.0;
        }
        else if (abs(nys) > SMALL)
        {
            t1x = - nys;
            t1y =   nxs;
            t1z =   0.0;
        }
        else if (abs(nzs) > SMALL)
        {
            t1x =   0.0;
            t1y = - nzs;
            t1z =   nys;
        }
        else
        {
            for (int m = 0; m < nLaminar; ++ m)
            {
                flux[m] = 0.0;
            }
            return;
        }

        //! Normalize the tangential vector.
        RDouble oNormal = 1.0 / DISTANCE(t1x, t1y, t1z);
        t1x *= oNormal;
        t1y *= oNormal;
        t1z *= oNormal;

        //! Get the second tangential vector by cross dot t1 to normal.
        t2x = nys * t1z - nzs * t1y;
        t2y = nzs * t1x - nxs * t1z;
        t2z = nxs * t1y - nys * t1x;

        RDouble dxL = xcc[le] - xfc[iFace];
        RDouble dyL = ycc[le] - yfc[iFace];
        RDouble dzL = zcc[le] - zfc[iFace];

        RDouble dxR = xcc[re] - xfc[iFace];
        RDouble dyR = ycc[re] - yfc[iFace];
        RDouble dzR = zcc[re] - zfc[iFace];

        RDouble dL = nxs * dxL + nys * dyL + nzs * dzL;
        RDouble dR = nxs * dxR + nys * dyR + nzs * dzR;

        RDouble distTemp = - dL / (DISTANCE(dxL, dyL, dzL) + SMALL);
        distTemp = MIN(distTemp,  1.0);
        distTemp = MAX(distTemp, -1.0);
        RDouble angle1 = asin(distTemp) * 180.0 / PI;
     
        distTemp = dR / (DISTANCE(dxR, dyR, dzR) + SMALL);
        distTemp = MIN(distTemp,  1.0);
        distTemp = MAX(distTemp, -1.0);
        RDouble angle2 = asin(distTemp) * 180.0 / PI;

        RDouble dxnL  = nxs * dL - dxL;
        RDouble dynL  = nys * dL - dyL;
        RDouble dznL  = nzs * dL - dzL;

        RDouble dxnR  = nxs * dR - dxR;
        RDouble dynR  = nys * dR - dyR;
        RDouble dznR  = nzs * dR - dzR;
        //! Quantities at points 1 and 2.
        for (int m = 0; m < nEquation; ++ m)
        {
            f1[m] = fL[m];
            f2[m] = fR[m];
        }

        T tL    = fL[IP]/(coefficientofStateEquation*fL[IR]);
        T tR    = fR[IP]/(coefficientofStateEquation*fR[IR]);
        T dtdxL = dqdxL[nEquation];
        T dtdyL = dqdyL[nEquation];
        T dtdzL = dqdzL[nEquation];
        T dtdxR = dqdxR[nEquation];
        T dtdyR = dqdyR[nEquation];
        T dtdzR = dqdzR[nEquation];

        for (int m = 0; m < nEquation; ++ m)
        {
            fMid[m] = half * (f1[m] + f2[m]);
            primitiveVariableFace[m] = fMid[m];
        }
        T tMid = half * (tL + tR);
        primitiveVariableFace[nEquation] = tMid;

         //! calculate the viscousLaminar here
        viscousLaminarL     = tL * sqrt(tL) * (1.0 + nonDimensionalSutherlandTemperature) / (tL + nonDimensionalSutherlandTemperature);
        viscousLaminarL     = MAX(viscousLaminarMin, viscousLaminarL);
        viscousLaminarR     = tR * sqrt(tR) * (1.0 + nonDimensionalSutherlandTemperature) / (tR + nonDimensionalSutherlandTemperature);
        viscousLaminarR     = MAX(viscousLaminarMin, viscousLaminarR);

        //! calculate the viscousTurbulence here
        RDouble nuoo  = 1.0;
        RDouble nuMin = 1.0e-5 * nuoo;
        RDouble nuMax = 1.0e10 * nuoo;
        RDouble SA_cv1_cube = GlobalDataBase::GetDoubleParaFromDB("SA_cv1_cube");

        ld  = fLturb * fL[0] / (viscousLaminarL + SMALL);
        ld3 = ld * ld * ld;
        fv1 = ld3 / (ld3 + SA_cv1_cube);
        viscousTurbulenceL = fL[0] * fv1 * fLturb;
        viscousTurbulenceL = min(nuMax, viscousTurbulenceL);
        viscousTurbulenceL = max(nuMin, viscousTurbulenceL);

        ld      = fRturb * fR[0] / (viscousLaminarR + SMALL);
        ld3     = ld * ld * ld;
        fv1     = ld3 / (ld3 + SA_cv1_cube);
        viscousTurbulenceR = fR[0] * fv1 * fRturb;
        viscousTurbulenceR = min(nuMax, viscousTurbulenceR);
        viscousTurbulenceR = max(nuMin, viscousTurbulenceR);

        if(iFace < nBoundFace)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
            int bcType = bcRegion->GetBCType();

            //! for solid walls, primitiveVariableFace[IU,IV,IW][jFace] = 0
            if (bcType == PHENGLEI::SOLID_SURFACE)
            {
                RDouble uWall = 0.0;
                RDouble vWall = 0.0;
                RDouble wWall = 0.0;
                Data_Param *bcData  = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace])->GetBCParamDataBase();
                if(bcData)
                {
                    if (bcData->IsExist("uWall", PHDOUBLE, 1))
                    {
                        bcData->GetData("uWall", &uWall, PHDOUBLE, 1);
                        bcData->GetData("vWall", &vWall, PHDOUBLE, 1);
                        bcData->GetData("wWall", &wWall, PHDOUBLE, 1);
                    }
                }

                T velocityXWall = uWall;
                T velocityYWall = vWall;
                T velocityZWall = wWall;
                primitiveVariableFace[IU] = velocityXWall;
                primitiveVariableFace[IV] = velocityYWall;
                primitiveVariableFace[IW] = velocityZWall;

                viscousTurbulenceR = -viscousTurbulenceL;
            }
            else if(bcType == PHENGLEI::FARFIELD)
            {
                viscousTurbulenceR = dDdP_turb[0][re-nTotalCell]*viscousTurbulenceL;
            }
        }

         //! calculate the viscousLaminarFace 
        viscousLaminarFace    = half*(viscousLaminarL    + viscousLaminarR);
        viscousTurbulenceFace = half*(viscousTurbulenceL + viscousTurbulenceR);

        RDouble cp = 1.0 / ((refGama - 1.0) * referenceMachNumberSquare);
        kCpFace = (viscousLaminarFace * oPrandtlLaminar + viscousTurbulenceFace * oPrandtlTurbulence) * cp;
        //! More accurate for skewness cells.
        //! If the cell degenerates too much, then interpolate the cell center value
        //! to the left and right points that vertical to the face normal.
        //! However, it will effect the robustness, so the skewnessAngle is suggested to be larger than 60 degree
        //! skewnessAngle, more large, more robust.
        //! For robust reason, this interpolation is NOT done on coarse grids.
        //! angle1/2 is defined as the complement of the angle between face normal and face-center/cell-center line.
        //! meaning that, (90 - angle1/2) is the angle between the face normal and the face-cell center line.
        //! angle1, angle2 = 0, represent completely orthogonal right angle, perfect quality.
        if (isFineGrid && angle1 > skewnessAngle && angle2 > skewnessAngle)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                f1[m] += dqdxL[m] * dxnL + dqdyL[m]* dynL + dqdzL[m]* dznL;
                f2[m] += dqdxR[m] * dxnR + dqdyR[m]* dynR + dqdzR[m]* dznR;
            }

            //! Only consider the rho limiter, may consider p limiter in future.
            tL += dtdxL * dxnL + dtdyL * dynL + dtdzL * dznL;
            tR += dtdxR * dxnR + dtdyR * dynR + dtdzR * dznR;

            if (tL < SMALL)
            {
                tL = fL[IP]/(coefficientofStateEquation*fL[IR]);
            }
            if (tR < SMALL)
            {
                tR = fR[IP]/(coefficientofStateEquation*fR[IR]);
            }

            //! Quantities at the face.
            for (int m = 0; m < nEquation; ++ m)
            {
                fMid[m] = primitiveVariableFace[m];
            }

            tMid = primitiveVariableFace[nEquation];
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            dfdn[m] = 0.0;
        }
        T dtdn = 0.0;

        if (angle1 > 0.0 && angle2 > 0.0 && ABS(dL) > TINY && ABS(dR) > TINY)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                dfd1[m] = (f1[m] - fMid[m]) / dL;
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                dfd2[m] = (f2[m] - fMid[m]) / dR;
            }

            T dtd1 = (tL - tMid) / dL;
            T dtd2 = (tR - tMid) / dR;

            RDouble dtmp = dL * dL + dR * dR;
            RDouble weightL = dL * dL / dtmp;
            RDouble weightR = 1.0 - weightL;

            if (iFace < nBoundFace)
            {
                UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
                int bcType = bcRegion->GetBCType();

                if (bcType != PHENGLEI::INTERFACE)
                {
                    weightL = 1.0;
                    weightR = 0.0;
                }
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                dfdn[m] = weightL * dfd1[m] + weightR * dfd2[m];
            }

            dtdn = dtd1 * weightL + dtd2 * weightR;
        }
        
        RDouble weightL = deltaL[jFace];
        RDouble weightR = deltaR[jFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            dfdx[m] = weightL * dqdxL[m] + weightR * dqdxR[m];
            dfdy[m] = weightL * dqdyL[m] + weightR * dqdyR[m];
            dfdz[m] = weightL * dqdzL[m] + weightR * dqdzR[m];
        }

        for (int m = IU; m <= IW; ++ m)
        {
            dfdt1[m] = t1x * dfdx[m] + t1y * dfdy[m] + t1z * dfdz[m];
            dfdt2[m] = t2x * dfdx[m] + t2y * dfdy[m] + t2z * dfdz[m];
        }

        //! Now true gradients. 
        for (int m = IU; m <= IW; ++ m)
        {
            dfdx[m] = nxs * dfdn[m] + t1x * dfdt1[m] + t2x * dfdt2[m];
            dfdy[m] = nys * dfdn[m] + t1y * dfdt1[m] + t2y * dfdt2[m];
            dfdz[m] = nzs * dfdn[m] + t1z * dfdt1[m] + t2z * dfdt2[m];
        }

        T dudx = dfdx[IU];
        T dudy = dfdy[IU];
        T dudz = dfdz[IU];

        T dvdx = dfdx[IV];
        T dvdy = dfdy[IV];
        T dvdz = dfdz[IV];

        T dwdx = dfdx[IW];
        T dwdy = dfdy[IW];
        T dwdz = dfdz[IW];

        T qNorm = 0.0;
        // if (nChemical == 1)
        // {
        //     for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
        //     {
        //         qNorm += rhoDiffusion[iSpecies][jFace] * hintSpecies[iSpecies][jFace] * dfdn[iSpecies + nm];
        //     }
// 
        //     for (int m = nm; m < nLaminar; ++ m)
        //     {
        //         fvis[m] = rhoDiffusion[m - nm][jFace] * dfdn[m];
        //     }
        // }
        
        qNorm += kCpFace * dtdn;
        T divv2p3 = two3rd * (dudx + dvdy + dwdz);

        T vis = viscousLaminarFace + viscousTurbulenceFace;
        T um  = fMid[1];
        T vm  = fMid[2];
        T wm  = fMid[3];

        //! Stress components.
        txx = vis * (two * dudx - divv2p3);
        tyy = vis * (two * dvdy - divv2p3);
        tzz = vis * (two * dwdz - divv2p3);
        txy = vis * (dudy + dvdx);
        txz = vis * (dudz + dwdx);
        tyz = vis * (dvdz + dwdy);

        if (nrokplus > 0)
        {
            T rkl  = fL[IR] * qTurb[0][le];
            T rkr  = fR[IR] * qTurb[0][re];
            T rhok = half * (rkl + rkr) * refReNumber;
            if (nolstress > 0)
            {
                RDouble b11 = half * (aniss[0][le] + aniss[0][re]);
                RDouble b22 = half * (aniss[1][le] + aniss[1][re]);
                RDouble b12 = half * (aniss[2][le] + aniss[2][re]);
                RDouble b13 = half * (aniss[3][le] + aniss[3][re]);
                RDouble b23 = half * (aniss[4][le] + aniss[4][re]);
                RDouble b33 = - b11 - b22;
                txx = txx + b11 - two3rd * rhok;
                tyy = tyy + b22 - two3rd * rhok;
                tzz = tzz + b33 - two3rd * rhok;
                txy = txy + b12;
                txz = txz + b13;
                tyz = tyz + b23;
            }
            else
            {
                txx = txx - two3rd * rhok;
                tyy = tyy - two3rd * rhok;
                tzz = tzz - two3rd * rhok;
            }
        }

        fvis[IR ] = 0.0;
        fvis[IRU] = nxs * txx + nys * txy + nzs * txz;
        fvis[IRV] = nxs * txy + nys * tyy + nzs * tyz;
        fvis[IRW] = nxs * txz + nys * tyz + nzs * tzz;
        fvis[IRE] = um * fvis[IRU] + vm * fvis[IRV] + wm * fvis[IRW] + qNorm;

        for (int m = 0; m < nLaminar; ++ m)
        {
            flux[m] = - oRefReNumber * area[iFace] * fvis[m];
        }

         delete [] fvis;    fvis = nullptr;
         delete [] f1;    f1 = nullptr;
         delete [] f2;    f2 = nullptr;
         delete [] fMid;    fMid = nullptr;
         delete [] dfdx;    dfdx = nullptr;
         delete [] dfdy;    dfdy = nullptr;
         delete [] dfdz;    dfdz = nullptr;
     
         delete [] dfd1;    dfd1 = nullptr;
         delete [] dfd2;    dfd2 = nullptr;
         delete [] dfdn;    dfdn = nullptr;
         delete [] dfdt1;    dfdt1 = nullptr;
         delete [] dfdt2;    dfdt2 = nullptr;
         delete [] primitiveVariableFace;    primitiveVariableFace = nullptr;
}

//! GMRESAD 
template <typename T>
void NSSolverUnstruct::Cal_GMRES_Visflux_AD(T* fL, T* fR, T* dqdxL,T* dqdyL, T* dqdzL, 
                          T* dqdxR, T* dqdyR, T* dqdzR, int iFace,
                          int jFace, T* flux, Grid *gridIn, 
                          FaceProxy *faceProxy, Param_NSSolverUnstruct *parameters)
{
        UnstructGrid *grid = UnstructGridCast(gridIn);
        int nBoundFace = grid->GetNBoundFace();

        int *leftCellofFace = grid->GetLeftCellOfFace();
        int *rightCellofFace = grid->GetRightCellOfFace();

        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

        RDouble *xfn  = grid->GetFaceNormalX();
        RDouble *yfn  = grid->GetFaceNormalY();
        RDouble *zfn  = grid->GetFaceNormalZ();
        RDouble *area = grid->GetFaceArea();

        RDouble *xcc  = grid->GetCellCenterX();
        RDouble *ycc  = grid->GetCellCenterY();
        RDouble *zcc  = grid->GetCellCenterZ();

        RDouble *xfc  = grid->GetFaceCenterX();
        RDouble *yfc  = grid->GetFaceCenterY();
        RDouble *zfc  = grid->GetFaceCenterZ();

        // Param_NSSolverUnstruct *parameters = GetControlParameters();
        int numberOfSpecies = parameters->GetNumberOfSpecies();

        int nLaminar  = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();

        int nEquation = nLaminar;

        int nm  = parameters->GetNSEquationNumber();
        RDouble skewnessAngle = parameters->GetSkewnessAngle();
        bool isFineGrid = grid->IsFinestGrid();
        RDouble *deltaL = faceProxy->GetWeightL();
        RDouble *deltaR = faceProxy->GetWeightR();
        NSFaceValue *faceVariable = faceProxy->GetNSFaceValue();
        // RDouble *kCp     = faceVariable->GetKCP();
        // RDouble *viscousLaminarFace   = faceVariable->GetMUL();
        RDouble *viscousTurbulenceFace = faceVariable->GetMUT();

        RDouble refReNumber  = parameters->GetRefReNumber();
        RDouble oRefReNumber = parameters->GetoRefReNumber();
        RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();

        int nolstress = GlobalDataBase::GetIntParaFromDB("nolstress");
        int nrokplus  = GlobalDataBase::GetIntParaFromDB("nrokplus");
        RDouble wallTemperature             = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
        RDouble refDimensionalTemperature   = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

        RDouble **qTurb = 0;
        RDouble ** aniss  = 0;

        if (nrokplus > 0)
        {
            qTurb = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
            aniss  = reinterpret_cast<RDouble **> (grid->GetDataPtr("aniss"));
        }

        T *fvis = new T[nLaminar];
        T *dfdx = new T[nEquation];
        T *dfdy = new T[nEquation];
        T *dfdz = new T[nEquation];

        T *dfd1  = new T[nEquation];
        T *dfd2  = new T[nEquation];
        T *dfdn  = new T[nEquation];
        T *dfdt1 = new T[nEquation];
        T *dfdt2 = new T[nEquation];

        T *f1    = new T[nEquation];
        T *f2    = new T[nEquation];
        T *fMid  = new T[nEquation];
        T *primitiveVariableFace = new T[nEquation+nTemperatureModel];

        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];

        RDouble nxs = xfn[iFace];
        RDouble nys = yfn[iFace];
        RDouble nzs = zfn[iFace];

        RDouble t1x, t1y, t1z, t2x, t2y, t2z;
        T txx, tyy, tzz;
        T txy, txz, tyz;

        T viscousLaminarMin,viscousLaminarL,viscousLaminarR,viscousLaminarFace;
        GlobalDataBase::GetData("visl_min", &viscousLaminarMin, PHDOUBLE, 1);

        RDouble refGama = 1.4;
        RDouble referenceMachNumber       = parameters->GetRefMachNumber();
        RDouble referenceMachNumberSquare = referenceMachNumber * referenceMachNumber;
        RDouble oPrandtlLaminar           = parameters->GetoPrandtlLaminar();
        RDouble oPrandtlTurbulence        = parameters->GetoPrandtlTurbulence();
        T kCpFace;
        RDouble coefficientofStateEquation = gas->GetCoefficientOfStateEquation();

        //! Get first tangential vector on the face.
        if (abs(nxs) > SMALL)
        {
            t1x =   nys;
            t1y = - nxs;
            t1z =   0.0;
        }
        else if (abs(nys) > SMALL)
        {
            t1x = - nys;
            t1y =   nxs;
            t1z =   0.0;
        }
        else if (abs(nzs) > SMALL)
        {
            t1x =   0.0;
            t1y = - nzs;
            t1z =   nys;
        }
        else
        {
            for (int m = 0; m < nLaminar; ++ m)
            {
                flux[m] = 0.0;
            }
            return;
        }

        //! Normalize the tangential vector.
        RDouble oNormal = 1.0 / DISTANCE(t1x, t1y, t1z);
        t1x *= oNormal;
        t1y *= oNormal;
        t1z *= oNormal;

        //! Get the second tangential vector by cross dot t1 to normal.
        t2x = nys * t1z - nzs * t1y;
        t2y = nzs * t1x - nxs * t1z;
        t2z = nxs * t1y - nys * t1x;

        RDouble dxL = xcc[le] - xfc[iFace];
        RDouble dyL = ycc[le] - yfc[iFace];
        RDouble dzL = zcc[le] - zfc[iFace];

        RDouble dxR = xcc[re] - xfc[iFace];
        RDouble dyR = ycc[re] - yfc[iFace];
        RDouble dzR = zcc[re] - zfc[iFace];

        RDouble dL = nxs * dxL + nys * dyL + nzs * dzL;
        RDouble dR = nxs * dxR + nys * dyR + nzs * dzR;

        RDouble distTemp = - dL / (DISTANCE(dxL, dyL, dzL) + SMALL);
        distTemp = MIN(distTemp,  1.0);
        distTemp = MAX(distTemp, -1.0);
        RDouble angle1 = asin(distTemp) * 180.0 / PI;
     
        distTemp = dR / (DISTANCE(dxR, dyR, dzR) + SMALL);
        distTemp = MIN(distTemp,  1.0);
        distTemp = MAX(distTemp, -1.0);
        RDouble angle2 = asin(distTemp) * 180.0 / PI;

        RDouble dxnL  = nxs * dL - dxL;
        RDouble dynL  = nys * dL - dyL;
        RDouble dznL  = nzs * dL - dzL;

        RDouble dxnR  = nxs * dR - dxR;
        RDouble dynR  = nys * dR - dyR;
        RDouble dznR  = nzs * dR - dzR;
        //! Quantities at points 1 and 2.
        for (int m = 0; m < nEquation; ++ m)
        {
            f1[m] = fL[m];
            f2[m] = fR[m];
        }

        T tL    = fL[IP] / (coefficientofStateEquation*fL[IR]);
        T tR    = fR[IP] / (coefficientofStateEquation*fR[IR]);
        T dtdxL = dqdxL[nEquation];
        T dtdyL = dqdyL[nEquation];
        T dtdzL = dqdzL[nEquation];
        T dtdxR = dqdxR[nEquation];
        T dtdyR = dqdyR[nEquation];
        T dtdzR = dqdzR[nEquation];
        //! calculate the viscousLaminar here
        viscousLaminarL    = tL * sqrt(tL) * (1.0 + nonDimensionalSutherlandTemperature) / (tL + nonDimensionalSutherlandTemperature);
        viscousLaminarL    = MAX(viscousLaminarMin, viscousLaminarL);
        viscousLaminarR    = tR * sqrt(tR) * (1.0 + nonDimensionalSutherlandTemperature) / (tR + nonDimensionalSutherlandTemperature);
        viscousLaminarR    = MAX(viscousLaminarMin, viscousLaminarR);
        viscousLaminarFace = half*(viscousLaminarL + viscousLaminarR);

        RDouble cp = 1.0 / ((refGama - 1.0) * referenceMachNumberSquare);
        kCpFace = (viscousLaminarFace * oPrandtlLaminar + viscousTurbulenceFace[jFace] * oPrandtlTurbulence) * cp;

        for (int m = 0; m < nEquation; ++ m)
        {
            fMid[m] = half * (f1[m] + f2[m]);
            primitiveVariableFace[m] = fMid[m];
        }
        T tMid = half * (tL + tR);
        primitiveVariableFace[nEquation] = tMid;

        if(iFace < nBoundFace)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
            int bcType = bcRegion->GetBCType();

            //! for solid walls, primitiveVariableFace[IU,IV,IW][jFace] = 0
            if (bcType == PHENGLEI::SOLID_SURFACE)
            {
                RDouble uWall = 0.0;
                RDouble vWall = 0.0;
                RDouble wWall = 0.0;
                Data_Param *bcData  = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace])->GetBCParamDataBase();
                if(bcData)
                {
                    if (bcData->IsExist("uWall", PHDOUBLE, 1))
                    {
                        bcData->GetData("uWall", &uWall, PHDOUBLE, 1);
                        bcData->GetData("vWall", &vWall, PHDOUBLE, 1);
                        bcData->GetData("wWall", &wWall, PHDOUBLE, 1);
                    }
                }

                T velocityXWall = uWall;
                T velocityYWall = vWall;
                T velocityZWall = wWall;
                primitiveVariableFace[IU] = velocityXWall;
                primitiveVariableFace[IV] = velocityYWall;
                primitiveVariableFace[IW] = velocityZWall;

                //! GMRES3D
                if( wallTemperature > 0.0 )
                {
                    tMid                             = wallTemperature/refDimensionalTemperature;
                    primitiveVariableFace[nEquation] = tMid;
                }
            }
        }
        //! More accurate for skewness cells.
        //! If the cell degenerates too much, then interpolate the cell center value
        //! to the left and right points that vertical to the face normal.
        //! However, it will effect the robustness, so the skewnessAngle is suggested to be larger than 60 degree
        //! skewnessAngle, more large, more robust.
        //! For robust reason, this interpolation is NOT done on coarse grids.
        //! angle1/2 is defined as the complement of the angle between face normal and face-center/cell-center line.
        //! meaning that, (90 - angle1/2) is the angle between the face normal and the face-cell center line.
        //! angle1, angle2 = 0, represent completely orthogonal right angle, perfect quality.
        if (isFineGrid && angle1 > skewnessAngle && angle2 > skewnessAngle)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                f1[m] += dqdxL[m] * dxnL + dqdyL[m]* dynL + dqdzL[m]* dznL;
                f2[m] += dqdxR[m] * dxnR + dqdyR[m]* dynR + dqdzR[m]* dznR;
            }

            //! Only consider the rho limiter, may consider p limiter in future.
            tL += dtdxL * dxnL + dtdyL * dynL + dtdzL * dznL;
            tR += dtdxR * dxnR + dtdyR * dynR + dtdzR * dznR;

            if (tL < SMALL)
            {
                tL = fL[IP]/(coefficientofStateEquation*fL[IR]);
            }
            if (tR < SMALL)
            {
                tR = fR[IP]/(coefficientofStateEquation*fR[IR]);
            }

            //! Quantities at the face.
            for (int m = 0; m < nEquation; ++ m)
            {
                fMid[m] = primitiveVariableFace[m];
            }

            tMid = primitiveVariableFace[nEquation];
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            dfdn[m] = 0.0;
        }
        T dtdn = 0.0;

        if (angle1 > 0.0 && angle2 > 0.0 && ABS(dL) > TINY && ABS(dR) > TINY)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                dfd1[m] = (f1[m] - fMid[m]) / dL;
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                dfd2[m] = (f2[m] - fMid[m]) / dR;
            }

            T dtd1 = (tL - tMid) / dL;
            T dtd2 = (tR - tMid) / dR;

            RDouble dtmp = dL * dL + dR * dR;
            RDouble weightL = dL * dL / dtmp;
            RDouble weightR = 1.0 - weightL;

            if (iFace < nBoundFace)
            {
                UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
                int bcType = bcRegion->GetBCType();

                if (bcType != PHENGLEI::INTERFACE)
                {
                    weightL = 1.0;
                    weightR = 0.0;
                }
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                dfdn[m] = weightL * dfd1[m] + weightR * dfd2[m];
            }

            dtdn = dtd1 * weightL + dtd2 * weightR;
        }
        
        RDouble weightL = deltaL[jFace];
        RDouble weightR = deltaR[jFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            dfdx[m] = weightL * dqdxL[m] + weightR * dqdxR[m];
            dfdy[m] = weightL * dqdyL[m] + weightR * dqdyR[m];
            dfdz[m] = weightL * dqdzL[m] + weightR * dqdzR[m];
        }

        for (int m = IU; m <= IW; ++ m)
        {
            dfdt1[m] = t1x * dfdx[m] + t1y * dfdy[m] + t1z * dfdz[m];
            dfdt2[m] = t2x * dfdx[m] + t2y * dfdy[m] + t2z * dfdz[m];
        }

        //! Now true gradients. 
        for (int m = IU; m <= IW; ++ m)
        {
            dfdx[m] = nxs * dfdn[m] + t1x * dfdt1[m] + t2x * dfdt2[m];
            dfdy[m] = nys * dfdn[m] + t1y * dfdt1[m] + t2y * dfdt2[m];
            dfdz[m] = nzs * dfdn[m] + t1z * dfdt1[m] + t2z * dfdt2[m];
        }
        
        T dudx = dfdx[IU];
        T dudy = dfdy[IU];
        T dudz = dfdz[IU];

        T dvdx = dfdx[IV];
        T dvdy = dfdy[IV];
        T dvdz = dfdz[IV];

        T dwdx = dfdx[IW];
        T dwdy = dfdy[IW];
        T dwdz = dfdz[IW];

        T qNorm = 0.0;
        // if (nChemical == 1)
        // {
        //     for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
        //     {
        //         qNorm += rhoDiffusion[iSpecies][jFace] * hintSpecies[iSpecies][jFace] * dfdn[iSpecies + nm];
        //     }
// 
        //     for (int m = nm; m < nLaminar; ++ m)
        //     {
        //         fvis[m] = rhoDiffusion[m - nm][jFace] * dfdn[m];
        //     }
        // }

        qNorm += kCpFace * dtdn;
        T divv2p3 = two3rd * (dudx + dvdy + dwdz);

        T vis = viscousLaminarFace + viscousTurbulenceFace[jFace];
        T um  = fMid[1];
        T vm  = fMid[2];
        T wm  = fMid[3];

        //! Stress components.
        txx = vis * (two * dudx - divv2p3);
        tyy = vis * (two * dvdy - divv2p3);
        tzz = vis * (two * dwdz - divv2p3);
        txy = vis * (dudy + dvdx);
        txz = vis * (dudz + dwdx);
        tyz = vis * (dvdz + dwdy);

        if (nrokplus > 0)
        {
            T rkl  = fL[IR] * qTurb[0][le];
            T rkr  = fR[IR] * qTurb[0][re];
            T rhok = half * (rkl + rkr) * refReNumber;
            if (nolstress > 0)
            {
                RDouble b11 = half * (aniss[0][le] + aniss[0][re]);
                RDouble b22 = half * (aniss[1][le] + aniss[1][re]);
                RDouble b12 = half * (aniss[2][le] + aniss[2][re]);
                RDouble b13 = half * (aniss[3][le] + aniss[3][re]);
                RDouble b23 = half * (aniss[4][le] + aniss[4][re]);
                RDouble b33 = - b11 - b22;
                txx = txx + b11 - two3rd * rhok;
                tyy = tyy + b22 - two3rd * rhok;
                tzz = tzz + b33 - two3rd * rhok;
                txy = txy + b12;
                txz = txz + b13;
                tyz = tyz + b23;
            }
            else
            {
                txx = txx - two3rd * rhok;
                tyy = tyy - two3rd * rhok;
                tzz = tzz - two3rd * rhok;
            }
        }

        fvis[IR ] = 0.0;
        fvis[IRU] = nxs * txx + nys * txy + nzs * txz;
        fvis[IRV] = nxs * txy + nys * tyy + nzs * tyz;
        fvis[IRW] = nxs * txz + nys * tyz + nzs * tzz;
        fvis[IRE] = um * fvis[IRU] + vm * fvis[IRV] + wm * fvis[IRW] + qNorm;

        for (int m = 0; m < nLaminar; ++ m)
        {
            flux[m] = - oRefReNumber * area[iFace] * fvis[m];
        }

         delete [] fvis;    fvis = nullptr;
         delete [] f1;    f1 = nullptr;
         delete [] f2;    f2 = nullptr;
         delete [] fMid;    fMid = nullptr;
         delete [] dfdx;    dfdx = nullptr;
         delete [] dfdy;    dfdy = nullptr;
         delete [] dfdz;    dfdz = nullptr;
     
         delete [] dfd1;    dfd1 = nullptr;
         delete [] dfd2;    dfd2 = nullptr;
         delete [] dfdn;    dfdn = nullptr;
         delete [] dfdt1;    dfdt1 = nullptr;
         delete [] dfdt2;    dfdt2 = nullptr;
         delete [] primitiveVariableFace;    primitiveVariableFace = nullptr;
}

//! GMRESVis  GMRESVisFD
void NSSolverUnstruct::Cal_GMRES_Visflux(RDouble* fL, RDouble* fR, RDouble* primitiveVariableFace, RDouble* dqdxL,
                           RDouble* dqdyL, RDouble* dqdzL, RDouble* dqdxR, RDouble* dqdyR, RDouble* dqdzR, int iFace,
                           int jFace, RDouble* flux, Grid *gridIn, FaceProxy *faceProxy)
{
        UnstructGrid *grid = UnstructGridCast(gridIn);
        int nBoundFace = grid->GetNBoundFace();

        int *leftCellofFace = grid->GetLeftCellOfFace();
        int *rightCellofFace = grid->GetRightCellOfFace();

        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

        RDouble *xfn  = grid->GetFaceNormalX();
        RDouble *yfn  = grid->GetFaceNormalY();
        RDouble *zfn  = grid->GetFaceNormalZ();
        RDouble *area = grid->GetFaceArea();
    
        RDouble *xcc  = grid->GetCellCenterX();
        RDouble *ycc  = grid->GetCellCenterY();
        RDouble *zcc  = grid->GetCellCenterZ();

        RDouble *xfc  = grid->GetFaceCenterX();
        RDouble *yfc  = grid->GetFaceCenterY();
        RDouble *zfc  = grid->GetFaceCenterZ();

        Param_NSSolverUnstruct *parameters = GetControlParameters();
        int numberOfSpecies = parameters->GetNumberOfSpecies();

        int nLaminar  = parameters->GetLaminarNumber();
        int nChemical = parameters->GetChemicalFlag();
        int nTemperatureModel = parameters->GetTemperatureModel();

        int nEquation = GetNumberOfEquations();

        int nm  = parameters->GetNSEquationNumber();
        RDouble skewnessAngle = parameters->GetSkewnessAngle();
        bool isFineGrid = grid->IsFinestGrid();
        RDouble *deltaL = faceProxy->GetWeightL();
        RDouble *deltaR = faceProxy->GetWeightR();
        NSFaceValue *faceVariable = faceProxy->GetNSFaceValue();
        // RDouble *kCp     = faceVariable->GetKCP();
        // RDouble *viscousLaminarFace   = faceVariable->GetMUL();
        RDouble *viscousTurbulenceFace = faceVariable->GetMUT();

        RDouble refReNumber  = parameters->GetRefReNumber();
        RDouble oRefReNumber = parameters->GetoRefReNumber();
        RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();

        int nolstress = GlobalDataBase::GetIntParaFromDB("nolstress");
        int nrokplus  = GlobalDataBase::GetIntParaFromDB("nrokplus");

        RDouble **qTurb = 0;
        RDouble ** aniss  = 0;

        if (nrokplus > 0)
        {
            qTurb = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
            aniss  = reinterpret_cast<RDouble **> (grid->GetDataPtr("aniss"));
        }

        RDouble *fvis = new RDouble[nLaminar];
        RDouble *dfdx = new RDouble[nEquation];
        RDouble *dfdy = new RDouble[nEquation];
        RDouble *dfdz = new RDouble[nEquation];

        RDouble *dfd1  = new RDouble[nEquation];
        RDouble *dfd2  = new RDouble[nEquation];
        RDouble *dfdn  = new RDouble[nEquation];
        RDouble *dfdt1 = new RDouble[nEquation];
        RDouble *dfdt2 = new RDouble[nEquation];

        RDouble *f1    = new RDouble[nEquation];
        RDouble *f2    = new RDouble[nEquation];
        RDouble *fMid  = new RDouble[nEquation];

        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];

        RDouble nxs = xfn[iFace];
        RDouble nys = yfn[iFace];
        RDouble nzs = zfn[iFace];

        RDouble t1x, t1y, t1z, t2x, t2y, t2z;
        RDouble txx, tyy, tzz;
        RDouble txy, txz, tyz;

        RDouble viscousLaminarMin,viscousLaminarL,viscousLaminarR,viscousLaminarFace;
        GlobalDataBase::GetData("visl_min", &viscousLaminarMin, PHDOUBLE, 1);

        RDouble refGama = 1.4;
        RDouble referenceMachNumber       = parameters->GetRefMachNumber();
        RDouble referenceMachNumberSquare = referenceMachNumber * referenceMachNumber;
        RDouble oPrandtlLaminar           = parameters->GetoPrandtlLaminar();
        RDouble oPrandtlTurbulence        = parameters->GetoPrandtlTurbulence();
        RDouble kCpFace;

        //! Get first tangential vector on the face.
        if (abs(nxs) > SMALL)
        {
            t1x =   nys;
            t1y = - nxs;
            t1z =   0.0;
        }
        else if (abs(nys) > SMALL)
        {
            t1x = - nys;
            t1y =   nxs;
            t1z =   0.0;
        }
        else if (abs(nzs) > SMALL)
        {
            t1x =   0.0;
            t1y = - nzs;
            t1z =   nys;
        }
        else
        {
            for (int m = 0; m < nLaminar; ++ m)
            {
                flux[m] = 0.0;
            }
            return;
        }

        //! Normalize the tangential vector.
        RDouble oNormal = 1.0 / DISTANCE(t1x, t1y, t1z);
        t1x *= oNormal;
        t1y *= oNormal;
        t1z *= oNormal;

        //! Get the second tangential vector by cross dot t1 to normal.
        t2x = nys * t1z - nzs * t1y;
        t2y = nzs * t1x - nxs * t1z;
        t2z = nxs * t1y - nys * t1x;

        RDouble dxL = xcc[le] - xfc[iFace];
        RDouble dyL = ycc[le] - yfc[iFace];
        RDouble dzL = zcc[le] - zfc[iFace];

        RDouble dxR = xcc[re] - xfc[iFace];
        RDouble dyR = ycc[re] - yfc[iFace];
        RDouble dzR = zcc[re] - zfc[iFace];

        RDouble dL = nxs * dxL + nys * dyL + nzs * dzL;
        RDouble dR = nxs * dxR + nys * dyR + nzs * dzR;

        RDouble distTemp = - dL / (DISTANCE(dxL, dyL, dzL) + SMALL);
        distTemp = MIN(distTemp,  1.0);
        distTemp = MAX(distTemp, -1.0);
        RDouble angle1 = asin(distTemp) * 180.0 / PI;
     
        distTemp = dR / (DISTANCE(dxR, dyR, dzR) + SMALL);
        distTemp = MIN(distTemp,  1.0);
        distTemp = MAX(distTemp, -1.0);
        RDouble angle2 = asin(distTemp) * 180.0 / PI;

        RDouble dxnL  = nxs * dL - dxL;
        RDouble dynL  = nys * dL - dyL;
        RDouble dznL  = nzs * dL - dzL;

        RDouble dxnR  = nxs * dR - dxR;
        RDouble dynR  = nys * dR - dyR;
        RDouble dznR  = nzs * dR - dzR;

        //! Quantities at points 1 and 2.
        for (int m = 0; m < nEquation; ++ m)
        {
            f1[m] = fL[m];
            f2[m] = fR[m];
        }

        RDouble tL    = fL[nEquation];
        RDouble tR    = fR[nEquation];
        RDouble dtdxL = dqdxL[nEquation];
        RDouble dtdyL = dqdyL[nEquation];
        RDouble dtdzL = dqdzL[nEquation];
        RDouble dtdxR = dqdxR[nEquation];
        RDouble dtdyR = dqdyR[nEquation];
        RDouble dtdzR = dqdzR[nEquation];

        //! calculate the viscousLaminar here
        viscousLaminarL    = tL * sqrt(tL) * (1.0 + nonDimensionalSutherlandTemperature) / (tL + nonDimensionalSutherlandTemperature);
        viscousLaminarL    = MAX(viscousLaminarMin, viscousLaminarL);
        viscousLaminarR    = tR * sqrt(tR) * (1.0 + nonDimensionalSutherlandTemperature) / (tR + nonDimensionalSutherlandTemperature);
        viscousLaminarR    = MAX(viscousLaminarMin, viscousLaminarR);
        viscousLaminarFace = half*(viscousLaminarL + viscousLaminarR);

        RDouble cp = 1.0 / ((refGama - 1.0) * referenceMachNumberSquare);
        kCpFace = (viscousLaminarFace * oPrandtlLaminar + viscousTurbulenceFace[jFace] * oPrandtlTurbulence) * cp;

        for (int m = 0; m < nEquation; ++ m)
        {
            fMid[m] = half * (f1[m] + f2[m]);
        }
        RDouble tMid = half * (tL + tR);

        //! More accurate for skewness cells.
        //! If the cell degenerates too much, then interpolate the cell center value
        //! to the left and right points that vertical to the face normal.
        //! However, it will effect the robustness, so the skewnessAngle is suggested to be larger than 60 degree
        //! skewnessAngle, more large, more robust.
        //! For robust reason, this interpolation is NOT done on coarse grids.
        //! angle1/2 is defined as the complement of the angle between face normal and face-center/cell-center line.
        //! meaning that, (90 - angle1/2) is the angle between the face normal and the face-cell center line.
        //! angle1, angle2 = 0, represent completely orthogonal right angle, perfect quality.
        if (isFineGrid && angle1 > skewnessAngle && angle2 > skewnessAngle)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                f1[m] += dqdxL[m] * dxnL + dqdyL[m]* dynL + dqdzL[m]* dznL;
                f2[m] += dqdxR[m] * dxnR + dqdyR[m]* dynR + dqdzR[m]* dznR;
            }

            //! Only consider the rho limiter, may consider p limiter in future.
            tL += dtdxL * dxnL + dtdyL * dynL + dtdzL * dznL;
            tR += dtdxR * dxnR + dtdyR * dynR + dtdzR * dznR;

            if (tL < SMALL)
            {
                tL = fL[nEquation];
            }
            if (tR < SMALL)
            {
                tR = fR[nEquation];
            }

            //! Quantities at the face.
            for (int m = 0; m < nEquation; ++ m)
            {
                fMid[m] = primitiveVariableFace[m];
            }

            tMid = primitiveVariableFace[nEquation];
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            dfdn[m] = 0.0;
        }
        RDouble dtdn = 0.0;

        if (angle1 > 0.0 && angle2 > 0.0 && ABS(dL) > TINY && ABS(dR) > TINY)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                dfd1[m] = (f1[m] - fMid[m]) / dL;
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                dfd2[m] = (f2[m] - fMid[m]) / dR;
            }

            RDouble dtd1 = (tL - tMid) / dL;
            RDouble dtd2 = (tR - tMid) / dR;

            RDouble dtmp = dL * dL + dR * dR;
            RDouble weightL = dL * dL / dtmp;
            RDouble weightR = 1.0 - weightL;

            if (iFace < nBoundFace)
            {
                UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
                int bcType = bcRegion->GetBCType();

                if (bcType != PHENGLEI::INTERFACE)
                {
                    weightL = 1.0;
                    weightR = 0.0;
                }
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                dfdn[m] = weightL * dfd1[m] + weightR * dfd2[m];
            }

            dtdn = dtd1 * weightL + dtd2 * weightR;
        }

        RDouble weightL = deltaL[jFace];
        RDouble weightR = deltaR[jFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            dfdx[m] = weightL * dqdxL[m] + weightR * dqdxR[m];
            dfdy[m] = weightL * dqdyL[m] + weightR * dqdyR[m];
            dfdz[m] = weightL * dqdzL[m] + weightR * dqdzR[m];
        }

        for (int m = IU; m <= IW; ++ m)
        {
            dfdt1[m] = t1x * dfdx[m] + t1y * dfdy[m] + t1z * dfdz[m];
            dfdt2[m] = t2x * dfdx[m] + t2y * dfdy[m] + t2z * dfdz[m];
        }

        //! Now true gradients. 
        for (int m = IU; m <= IW; ++ m)
        {
            dfdx[m] = nxs * dfdn[m] + t1x * dfdt1[m] + t2x * dfdt2[m];
            dfdy[m] = nys * dfdn[m] + t1y * dfdt1[m] + t2y * dfdt2[m];
            dfdz[m] = nzs * dfdn[m] + t1z * dfdt1[m] + t2z * dfdt2[m];
        }

        RDouble dudx = dfdx[IU];
        RDouble dudy = dfdy[IU];
        RDouble dudz = dfdz[IU];

        RDouble dvdx = dfdx[IV];
        RDouble dvdy = dfdy[IV];
        RDouble dvdz = dfdz[IV];

        RDouble dwdx = dfdx[IW];
        RDouble dwdy = dfdy[IW];
        RDouble dwdz = dfdz[IW];

        RDouble qNorm = 0.0;
        // if (nChemical == 1)
        // {
        //     for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
        //     {
        //         qNorm += rhoDiffusion[iSpecies][jFace] * hintSpecies[iSpecies][jFace] * dfdn[iSpecies + nm];
        //     }
// 
        //     for (int m = nm; m < nLaminar; ++ m)
        //     {
        //         fvis[m] = rhoDiffusion[m - nm][jFace] * dfdn[m];
        //     }
        // }

        qNorm += kCpFace * dtdn;
        RDouble divv2p3 = two3rd * (dudx + dvdy + dwdz);

        RDouble vis = viscousLaminarFace + viscousTurbulenceFace[jFace];
        RDouble um  = fMid[1];
        RDouble vm  = fMid[2];
        RDouble wm  = fMid[3];

        //! Stress components.
        txx = vis * (two * dudx - divv2p3);
        tyy = vis * (two * dvdy - divv2p3);
        tzz = vis * (two * dwdz - divv2p3);
        txy = vis * (dudy + dvdx);
        txz = vis * (dudz + dwdx);
        tyz = vis * (dvdz + dwdy);

        if (nrokplus > 0)
        {
            RDouble rkl  = fL[IR] * qTurb[0][le];
            RDouble rkr  = fR[IR] * qTurb[0][re];
            RDouble rhok = half * (rkl + rkr) * refReNumber;
            if (nolstress > 0)
            {
                RDouble b11 = half * (aniss[0][le] + aniss[0][re]);
                RDouble b22 = half * (aniss[1][le] + aniss[1][re]);
                RDouble b12 = half * (aniss[2][le] + aniss[2][re]);
                RDouble b13 = half * (aniss[3][le] + aniss[3][re]);
                RDouble b23 = half * (aniss[4][le] + aniss[4][re]);
                RDouble b33 = - b11 - b22;
                txx = txx + b11 - two3rd * rhok;
                tyy = tyy + b22 - two3rd * rhok;
                tzz = tzz + b33 - two3rd * rhok;
                txy = txy + b12;
                txz = txz + b13;
                tyz = tyz + b23;
            }
            else
            {
                txx = txx - two3rd * rhok;
                tyy = tyy - two3rd * rhok;
                tzz = tzz - two3rd * rhok;
            }
        }

        fvis[IR ] = 0.0;
        fvis[IRU] = nxs * txx + nys * txy + nzs * txz;
        fvis[IRV] = nxs * txy + nys * tyy + nzs * tyz;
        fvis[IRW] = nxs * txz + nys * tyz + nzs * tzz;
        fvis[IRE] = um * fvis[IRU] + vm * fvis[IRV] + wm * fvis[IRW] + qNorm;

        for (int m = 0; m < nLaminar; ++ m)
        {
            flux[m] = - oRefReNumber * area[iFace] * fvis[m];
        }

         delete [] fvis;    fvis = nullptr;
         delete [] f1;    f1 = nullptr;
         delete [] f2;    f2 = nullptr;
         delete [] fMid;    fMid = nullptr;
         delete [] dfdx;    dfdx = nullptr;
         delete [] dfdy;    dfdy = nullptr;
         delete [] dfdz;    dfdz = nullptr;

         delete [] dfd1;    dfd1 = nullptr;
         delete [] dfd2;    dfd2 = nullptr;
         delete [] dfdn;    dfdn = nullptr;
         delete [] dfdt1;    dfdt1 = nullptr;
         delete [] dfdt2;    dfdt2 = nullptr;
}
#endif

void NSSolverUnstruct::ComputeVisflux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
#ifdef USE_CUDA
    Param_NSSolverUnstruct *parameters_gpu = GetControlParameters();
    int nEquation_gpu = GetNumberOfEquations();
    CallGPUComputeVisflux(gridIn, parameters_gpu, nEquation_gpu, localStart, localEnd);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();

    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    RDouble **dqdx = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarZ"));

    RDouble **dtdx = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradTemperatureX"));
    RDouble **dtdy = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradTemperatureY"));
    RDouble **dtdz = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradTemperatureZ"));

    RDouble  **flux = faceProxy->GetFlux();

    RDouble *deltaL = faceProxy->GetWeightL();
    RDouble *deltaR = faceProxy->GetWeightR();

    NSFaceValue *faceVariable = faceProxy->GetNSFaceValue();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    RDouble **rhoDiffusion = faceVariable->GetRhoDS()->GetField();
    RDouble **hintSpecies = faceVariable->GetHintS()->GetField();
    RDouble **primitiveVariableFace = faceVariable->GetPrim()->GetField();
    RDouble **tm = faceVariable->GetT()->GetField();

    RDouble *kCp = faceVariable->GetKCP();
    RDouble *viscousLaminarFace    = faceVariable->GetMUL();
    RDouble *viscousTurbulenceFace = faceVariable->GetMUT();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int numberOfSpecies = parameters->GetNumberOfSpecies();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();

    int nEquation = GetNumberOfEquations();

    int nm = parameters->GetNSEquationNumber();

    RDouble refReNumber  = parameters->GetRefReNumber();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    int nolstress = GlobalDataBase::GetIntParaFromDB("nolstress");
    int nrokplus  = GlobalDataBase::GetIntParaFromDB("nrokplus");

    RDouble skewnessAngle = parameters->GetSkewnessAngle();

    RDouble **t = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));
    RDouble **primitiveVariable = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble **qTurb = 0;
    RDouble **aniss  = 0;

    if (nrokplus > 0)
    {
        qTurb = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
        aniss  = reinterpret_cast<RDouble **> (grid->GetDataPtr("aniss"));
    }

    RDouble *fvis = new RDouble[nLaminar];

    RDouble *dfdx = new RDouble[nEquation];
    RDouble *dfdy = new RDouble[nEquation];
    RDouble *dfdz = new RDouble[nEquation];

    RDouble *dfd1  = new RDouble[nEquation];
    RDouble *dfd2  = new RDouble[nEquation];
    RDouble *dfdn  = new RDouble[nEquation];
    RDouble *dfdt1 = new RDouble[nEquation];
    RDouble *dfdt2 = new RDouble[nEquation];

    RDouble *f1   = new RDouble[nEquation];
    RDouble *f2   = new RDouble[nEquation];
    RDouble *fMid = new RDouble[nEquation];

    RDouble t1x, t1y, t1z, t2x, t2y, t2z;
    RDouble txx, tyy, tzz;
    RDouble txy, txz, tyz;

    bool isFineGrid = grid->IsFinestGrid();

    using namespace GAS_SPACE;
    RDouble coefficientofstateEquation = gas->GetCoefficientOfStateEquation();

    using namespace IDX;
    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];
        int jFace  = iFace - localStart;

        RDouble nxs = xfn[iFace];
        RDouble nys = yfn[iFace];
        RDouble nzs = zfn[iFace];

        //! Get first tangential vector on the face.
        if (abs(nxs) > SMALL)
        {
            t1x =   nys;
            t1y = - nxs;
            t1z =   0.0;
        }
        else if (abs(nys) > SMALL)
        {
            t1x = - nys;
            t1y =   nxs;
            t1z =   0.0;
        }
        else if (abs(nzs) > SMALL)
        {
            t1x =   0.0;
            t1y = - nzs;
            t1z =   nys;
        }
        else
        {
            for (int m = 0; m < nLaminar; ++ m)
            {
                flux[m][jFace] = 0.0;
            }
            continue;
        }

        //! Normalize the tangential vector.
        RDouble oNormal = 1.0 / DISTANCE(t1x, t1y, t1z);
        t1x *= oNormal;
        t1y *= oNormal;
        t1z *= oNormal;

        //! Get the second tangential vector by cross dot t1 to normal.
        t2x = nys * t1z - nzs * t1y;
        t2y = nzs * t1x - nxs * t1z;
        t2z = nxs * t1y - nys * t1x;

        RDouble dxL = xcc[le] - xfc[iFace];
        RDouble dyL = ycc[le] - yfc[iFace];
        RDouble dzL = zcc[le] - zfc[iFace];

        RDouble dxR = xcc[re] - xfc[iFace];
        RDouble dyR = ycc[re] - yfc[iFace];
        RDouble dzR = zcc[re] - zfc[iFace];

        RDouble dL = nxs * dxL + nys * dyL + nzs * dzL;
        RDouble dR = nxs * dxR + nys * dyR + nzs * dzR;

        RDouble distTemp = - dL / (DISTANCE(dxL, dyL, dzL) + SMALL);
        distTemp = MIN(distTemp,  1.0);
        distTemp = MAX(distTemp, -1.0);
        RDouble angle1 = asin(distTemp) * 180.0 / PI;

        distTemp = dR / (DISTANCE(dxR, dyR, dzR) + SMALL);
        distTemp = MIN(distTemp,  1.0);
        distTemp = MAX(distTemp, -1.0);
        RDouble angle2 = asin(distTemp) * 180.0 / PI;

        RDouble dxnL  = nxs * dL - dxL;
        RDouble dynL  = nys * dL - dyL;
        RDouble dznL  = nzs * dL - dzL;

        RDouble dxnR  = nxs * dR - dxR;
        RDouble dynR  = nys * dR - dyR;
        RDouble dznR  = nzs * dR - dzR;

        //! Quantities at points 1 and 2.
        for (int m = 0; m < nEquation; ++ m)
        {
            f1[m] = primitiveVariable[m][le];
            f2[m] = primitiveVariable[m][re];
        }

        RDouble tL = t[ITT][le];
        RDouble tR = t[ITT][re];

        for (int m = 0; m < nEquation; ++ m)
        {
            fMid[m] = half * (f1[m] + f2[m]);
        }
        RDouble tMid = half * (tL + tR);

        if (iFace < nBoundFace)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
            int bcType = bcRegion->GetBCType();
            if (IsWall(bcType))
            {
                tMid = tm[ITT][jFace];
                fMid[IDX::IR] = fMid[IDX::IP] / (coefficientofstateEquation * tMid);
            }
        }

        //! More accurate for skewness cells.
        //! If the celLoadResidualsl degenerates too much, then interpolate the cell center value
        //! to the left and right points that vertical to the face normal.
        //! However, it will effect the robustness, so the skewnessAngle is suggested to be larger than 60 degree
        //! skewnessAngle, more large, more robust.
        //! For robust reason, this interpolation is NOT done on coarse grids.
        //! angle1/2 is defined as the complement of the angle between face normal and face-center/cell-center line.
        //! meaning that, (90 - angle1/2) is the angle between the face normal and the face-cell center line.
        //! angle1, angle2 = 0, represent completely orthogonal right angle, perfect quality.
        if (isFineGrid && angle1 > skewnessAngle && angle2 > skewnessAngle)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                f1[m] += dqdx[m][le] * dxnL + dqdy[m][le] * dynL + dqdz[m][le] * dznL;
                f2[m] += dqdx[m][re] * dxnR + dqdy[m][re] * dynR + dqdz[m][re] * dznR;
            }

            //! Only consider the rho limiter, may consider p limiter in future.
            tL += dtdx[ITT][le] * dxnL + dtdy[ITT][le] * dynL + dtdz[ITT][le] * dznL;
            tR += dtdx[ITT][re] * dxnR + dtdy[ITT][re] * dynR + dtdz[ITT][re] * dznR;

            if (tL < SMALL)
            {
                tL = t[ITT][le];
            }
            if (tR < SMALL)
            {
                tR = t[ITT][re];
            }

            //! Quantities at the face.
            for (int m = 0; m < nEquation; ++ m)
            {
                fMid[m] = primitiveVariableFace[m][jFace];
            }

            tMid = tm[ITT][jFace];
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            dfdn[m] = 0.0;
        }
        RDouble dtdn = 0.0;

        if (angle1 > 0.0 && angle2 > 0.0 && ABS(dL) > TINY && ABS(dR) > TINY)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                dfd1[m] = (f1[m] - fMid[m]) / dL;
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                dfd2[m] = (f2[m] - fMid[m]) / dR;
            }

            RDouble dtd1 = (tL - tMid) / dL;
            RDouble dtd2 = (tR - tMid) / dR;

            RDouble dtmp = dL * dL + dR * dR;
            RDouble weightL = dL * dL / dtmp;
            RDouble weightR = 1.0 - weightL;

            if (iFace < nBoundFace)
            {
                UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
                int bcType = bcRegion->GetBCType();

                if (bcType != PHENGLEI::INTERFACE && bcType != PHENGLEI::OVERSET)
                {
                    weightL = 1.0;
                    weightR = 0.0;
                }
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                dfdn[m] = weightL * dfd1[m] + weightR * dfd2[m];
            }

            dtdn = dtd1 * weightL + dtd2 * weightR;
        }

        RDouble weightL = deltaL[jFace];
        RDouble weightR = deltaR[jFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            dfdx[m] = weightL * dqdx[m][le] + weightR * dqdx[m][re];
            dfdy[m] = weightL * dqdy[m][le] + weightR * dqdy[m][re];
            dfdz[m] = weightL * dqdz[m][le] + weightR * dqdz[m][re];
        }

        for (int m = IU; m <= IW; ++ m)
        {
            dfdt1[m] = t1x * dfdx[m] + t1y * dfdy[m] + t1z * dfdz[m];
            dfdt2[m] = t2x * dfdx[m] + t2y * dfdy[m] + t2z * dfdz[m];
        }

        //! Now true gradients. 
        for (int m = IU; m <= IW; ++ m)
        {
            dfdx[m] = nxs * dfdn[m] + t1x * dfdt1[m] + t2x * dfdt2[m];
            dfdy[m] = nys * dfdn[m] + t1y * dfdt1[m] + t2y * dfdt2[m];
            dfdz[m] = nzs * dfdn[m] + t1z * dfdt1[m] + t2z * dfdt2[m];
        }

        RDouble dudx = dfdx[IU];
        RDouble dudy = dfdy[IU];
        RDouble dudz = dfdz[IU];

        RDouble dvdx = dfdx[IV];
        RDouble dvdy = dfdy[IV];
        RDouble dvdz = dfdz[IV];

        RDouble dwdx = dfdx[IW];
        RDouble dwdy = dfdy[IW];
        RDouble dwdz = dfdz[IW];

        RDouble qNorm = 0.0;
        if (nChemical == 1)
        {
            for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
            {
                qNorm += rhoDiffusion[iSpecies][jFace] * hintSpecies[iSpecies][jFace] * dfdn[iSpecies + nm];
            }

            for (int m = nm; m < nLaminar; ++ m)
            {
                fvis[m] = rhoDiffusion[m - nm][jFace] * dfdn[m];
            }
        }

        qNorm += kCp[jFace] * dtdn;
        RDouble divv2p3 = two3rd * (dudx + dvdy + dwdz);

        RDouble vis = viscousLaminarFace[jFace] + viscousTurbulenceFace[jFace];
        RDouble um  = fMid[1];
        RDouble vm  = fMid[2];
        RDouble wm  = fMid[3];

        //! Stress components.
        txx = vis * (two * dudx - divv2p3);
        tyy = vis * (two * dvdy - divv2p3);
        tzz = vis * (two * dwdz - divv2p3);
        txy = vis * (dudy + dvdx);
        txz = vis * (dudz + dwdx);
        tyz = vis * (dvdz + dwdy);

        if (nrokplus > 0)
        {
            RDouble rkl  = primitiveVariable[IR][le] * qTurb[0][le];
            RDouble rkr  = primitiveVariable[IR][re] * qTurb[0][re];
            RDouble rhok = half * (rkl + rkr) * refReNumber;
            if (nolstress > 0)
            {
                RDouble b11 = half * (aniss[0][le] + aniss[0][re]);
                RDouble b22 = half * (aniss[1][le] + aniss[1][re]);
                RDouble b12 = half * (aniss[2][le] + aniss[2][re]);
                RDouble b13 = half * (aniss[3][le] + aniss[3][re]);
                RDouble b23 = half * (aniss[4][le] + aniss[4][re]);
                RDouble b33 = - b11 - b22;
                txx = txx + b11 - two3rd * rhok;
                tyy = tyy + b22 - two3rd * rhok;
                tzz = tzz + b33 - two3rd * rhok;
                txy = txy + b12;
                txz = txz + b13;
                tyz = tyz + b23;
            }
            else
            {
                txx = txx - two3rd * rhok;
                tyy = tyy - two3rd * rhok;
                tzz = tzz - two3rd * rhok;
            }
        }

        fvis[IR ] = 0.0;
        fvis[IRU] = nxs * txx + nys * txy + nzs * txz;
        fvis[IRV] = nxs * txy + nys * tyy + nzs * tyz;
        fvis[IRW] = nxs * txz + nys * tyz + nzs * tzz;
        fvis[IRE] = um * fvis[IRU] + vm * fvis[IRV] + wm * fvis[IRW] + qNorm;

        for (int m = 0; m < nLaminar; ++ m)
        {
            flux[m][jFace] = - oRefReNumber * area[iFace] * fvis[m];
        }
    }

    delete [] fvis;    fvis = nullptr;

    delete [] dfdx;    dfdx = nullptr;
    delete [] dfdy;    dfdy = nullptr;
    delete [] dfdz;    dfdz = nullptr;

    delete [] dfd1;    dfd1 = nullptr;
    delete [] dfd2;    dfd2 = nullptr;
    delete [] dfdn;    dfdn = nullptr;
    delete [] dfdt1;    dfdt1 = nullptr;
    delete [] dfdt2;    dfdt2 = nullptr;

    delete [] f1;    f1 = nullptr;
    delete [] f2;    f2 = nullptr;
    delete [] fMid;    fMid = nullptr;
}

void NSSolverUnstruct::CalculateBoundaryData()
{
    CalculateMassFluxRatio();
}

void NSSolverUnstruct::Boundary(Grid *gridIn)
{
#ifdef USE_CUDA
    using namespace GPUMemory;
    using namespace GPUControlVariables;
    Param_NSSolverUnstruct *parameters_gpu = GetControlParameters();
    CallGPUBoundary(gridIn, parameters_gpu);
    return;
#endif   
    using namespace PHENGLEI;
    UnstructGrid *grid = UnstructGridCast(gridIn);

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble wallTemperature = parameters->GetWallTemperature();
    int isOverset = parameters->GetIsOverLapping();

    if (isOverset)
    {
        ReSetOversetBoundary(grid);
    }

    CalculateBoundaryData();

    bool isViscous = parameters->IsViscous();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");
    //! Parallel
    int currentZone = grid->GetOrdinaryGridIndex();
    //! Serial
    if (currentZone == -1)
    {
        currentZone = grid->GetZoneID();
    }

    int MixingPlaneNumber = 0;
    string MixingIn, MixingOut;
    if (referenceFrame == ROTATIONAL_FRAME)
    {
        int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
        string MixingPlane[100];
        GlobalDataBase::GetData("MixingPlane", &MixingPlane, PHSTRING, 2 * nTurboZone - 2);

        if (currentZone == 0)
        {
            MixingOut = MixingPlane[2 * currentZone];
        }
        else if (currentZone == nTurboZone)
        {
            MixingIn = MixingPlane[2 * currentZone - 1];
        }
        else
        {
            MixingIn  = MixingPlane[2 * currentZone - 1];
            MixingOut = MixingPlane[2 * currentZone];
        }
    }

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();
        if (bcName == MixingIn || bcName == MixingOut)
        {
            continue;
        }

        if (IsInterface(bcType)|| bcType == OVERSET)
        {
            continue;
        }
        else if (bcType == EXTRAPOLATION)
        {
            OutflowBCRegion(grid, bcRegion);
        }
        else if (IsWall(bcType))
        {
            if (!isViscous)
            {
                SymmetryBCRegion(grid, bcRegion);
            }
            else if (wallTemperature <= 0.0)
            {
                ViscousAdiabaticWallBCRegion(gridIn, bcRegion);
            }
            else
            {
                ViscousIsotropicWallBCRegion(grid, bcRegion);
            }
        }
        else if (bcType == SYMMETRY)
        {
            SymmetryBCRegion(grid, bcRegion);
        }
        else if (bcType == FARFIELD)
        {
            FarfieldBCRegion(grid, bcRegion);
        }
        else if (bcType == INFLOW)
        {
            InflowBCRegion(gridIn, bcRegion);
        }
        else if (bcType == OUTFLOW)
        {
            OutflowBCRegion(gridIn, bcRegion);
        }
        else if (bcType == PRESSURE_OUTLET)
        {
            PressureOutletBCRegion(gridIn, bcRegion);
        }
        else if (bcType == PRESSURE_INLET)
        {
            PressureInletBCRiemannRegion(gridIn, bcRegion);
        }
        else if (bcType == MASS_FLOW_INLET)
        {
            MassFlowInletBCRegion(gridIn, bcRegion);
        }
        else if (bcType == MASS_FLOW_OUTLET)
        {
            MassFlowOutletBCRegion(gridIn, bcRegion);
        }
        else
        {
            TK_Exit::ExceptionExit("Error: this boundary type does not exist!\n");
        }
    }
}

//! this is from Geo_UnstructGrid
void NSSolverUnstruct::CornerPoint(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **qnode = reinterpret_cast <RDouble **> (grid->GetDataPtr("qnode"));

    int nTotalNode = grid->GetNTotalNode();
    int nBoundFace = grid->GetNBoundFace();
    int *face2node = grid->GetFace2Node();
    int nTotalCell = grid->GetNTotalCell();
    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();
    //! corner point of periodic boundary and solid wall boundary.
    int * corner_point = new int[nTotalNode];

    for (int ipoint = 0; ipoint < nTotalNode; ++ipoint)
    {
        corner_point[ipoint] = 0;
    }

    int nodepos = 0;
    int pt;

    //! find corner point.
    for (int iFace = 0; iFace < nBoundFace; iFace++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        if (bcType == -1)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++j)
            {
                pt = face2node[nodepos + j];
                corner_point[pt] = 1;
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }
    nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; iFace++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();
        if (bcType == 2)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++j)
            {
                pt = face2node[nodepos + j];
                if (corner_point[pt] == 1)
                {
                    corner_point[pt] = 2;
                }
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

    int le, re, bcType;
    nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++iFace)
    {
        le = leftCellOfFace[iFace];
        re = iFace + nTotalCell;

        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        bcType = bcRegion->GetBCType();

        if (bcType == 2)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++j)
            {
                pt = face2node[nodepos + j];

                //! for corner point
                if (corner_point[pt] == 2)
                {
                    for (int m = 1; m <= 3; m++)
                    {
                        q[m][iFace]=0;
                    }
                }
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

    delete [] corner_point;    corner_point = nullptr;
}

void NSSolverUnstruct::ComputeMassFlow(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    using namespace PHENGLEI;

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        
        if (bcType == INFLOW)
        {
            ComputeMassInFlow(gridIn, bcRegion);
        }
        else if (bcType == OUTFLOW)
        {
            ComputeMassOutFlow(gridIn, bcRegion);
        }
    }
}

void NSSolverUnstruct::ComputeMassInFlow(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    RDouble massIn  = 0.0;
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        //! iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

        RDouble rmi, umi, vmi, wmi, vni;
        RDouble rmo, umo, vmo, wmo, vno;

        rmi = q[IR][le];
        umi = q[IU][le];
        vmi = q[IV][le];
        wmi = q[IW][le];
        vni = xfn[iFace] * umi + yfn[iFace] * vmi + zfn[iFace] * wmi;

        rmo = q[IR][re];
        umo = q[IU][re];
        vmo = q[IV][re];
        wmo = q[IW][re];
        vno = xfn[iFace] * umo + yfn[iFace] * vmo + zfn[iFace] * wmo;

        massIn -= half * (rmi * vni + rmo * vno);
    }
    GlobalDataBase::UpdateData("massIn", &massIn, PHDOUBLE, 1);
}

void NSSolverUnstruct::ComputeMassOutFlow(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    RDouble massOut = 0.0;
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        //! iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellofFace[iFace];

        RDouble rm, um, vm, wm, vn;
        rm = q[IR][le];
        um = q[IU][le];
        vm = q[IV][le];
        wm = q[IW][le];
        vn = xfn[iFace] * um + yfn[iFace] * vm + zfn[iFace] * wm;

        massOut += rm * vn;
    }
    GlobalDataBase::UpdateData("massOut", &massOut, PHDOUBLE, 1);
}

void NSSolverUnstruct::BoundaryQlQrFix(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();
    int nMid;
    int nTotalCell = grid->GetNTotalCell();

    //! GMRESBoundary
    faceProxy->SetlocalStart(localStart);
    faceProxy->SetlocalEnd(localEnd);
    faceProxy->SetnBoundFace(nBoundFace);
    faceProxy->SetnTotalCell(nTotalCell);

    //! Check if there are boundary faces. If no, return.
    if (localStart >= nBoundFace)
    {
        return;
    }

    if (localEnd <= nBoundFace)
    {
        //! If they are all boundary faces.
        nMid = localEnd;
    }
    else
    {
        //! Part of them are boundary faces.
        nMid = nBoundFace;
    }

    //! Get grid information.
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *faceVelocityX = grid->GetFaceVelocityX();
    RDouble *faceVelocityY = grid->GetFaceVelocityY();
    RDouble *faceVelocityZ = grid->GetFaceVelocityZ();
    RDouble *vgn           = grid->GetFaceNormalVelocity();
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nChemical = parameters->GetChemicalFlag();

    int nEquation = GetNumberOfEquations();

    bool isViscous = parameters->IsViscous();
    RDouble tWall = parameters->GetWallTemperature();
    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
    RDouble coefficientOfStateEquation = gas->GetCoefficientOfStateEquation();

    int isUnsteady = parameters->GetIsUnsteady();
    int isAle = parameters->GetIsCodeOfAleModel();

    using namespace PHSPACE;
    using namespace PHENGLEI;

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();

    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble *primitiveVariable = new RDouble[nEquation];
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");    //! GMRES3D

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");

    int *nodeBCType = reinterpret_cast <int *> (grid->GetDataPtr("nodeBCType"));
    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();
    int *face2node = grid->GetFace2Node();

    using namespace IDX;

#ifdef USE_GMRESSOLVER
    // GMRESBCorrection
    if(localStart < nMid)
    {
        int nBCRegion = unstructBCSet->GetnBCRegion();
        for(int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion++)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
            int bcType = bcRegion->GetBCType();
            if(IsWall(bcType))
            {
                vector<int>* faceIndex = bcRegion->GetFaceIndex();
                faceProxy->SetWallType(isViscous);
                faceProxy->SetWallFaceIndex(faceIndex);
            }
        }
    }
#endif

    for (int iFace = localStart; iFace < nMid; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();
        int jFace = iFace - localStart;
        int nodepos = 0;

        if (bcType != PHENGLEI::SOLID_SURFACE && bcType != PHENGLEI::SYMMETRY)
        {
            nodepos += node_number_of_each_face[iFace];
            continue;
        }

        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            RFloat temp =  half * (q[m][le] + q[m][re]);
            qL[m][jFace] = temp;
            qR[m][jFace] = temp;
        }

        if( tscheme == GMRES )
        {
            continue;
        }

        if (bcType == PHENGLEI::SYMMETRY)
        {
            RDouble vn = 2.0 * (xfn[iFace] * qL[IU][jFace] + yfn[iFace] * qL[IV][jFace] + zfn[iFace] * qL[IW][jFace] - vgn[iFace]);

            qR[IR][jFace] = qL[IR][jFace];
            qR[IU][jFace] = qL[IU][jFace] - vn * xfn[iFace];
            qR[IV][jFace] = qL[IV][jFace] - vn * yfn[iFace];
            qR[IW][jFace] = qL[IW][jFace] - vn * zfn[iFace];
            qR[IP][jFace] = qL[IP][jFace];
        }

        if (bcType == PHENGLEI::SOLID_SURFACE && isViscous)
        {
            RDouble uWall = 0.0;
            RDouble vWall = 0.0;
            RDouble wWall = 0.0;
            Data_Param *bcData  = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace])->GetBCParamDataBase();
            if(bcData)
            {
                if (bcData->IsExist("uWall", PHDOUBLE, 1))
                {
                    bcData->GetData("uWall", &uWall, PHDOUBLE, 1);
                    bcData->GetData("vWall", &vWall, PHDOUBLE, 1);
                    bcData->GetData("wWall", &wWall, PHDOUBLE, 1);
                }
            }

            RDouble velocityXWall = uWall;
            RDouble velocityYWall = vWall;
            RDouble velocityZWall = wWall;
            if (isUnsteady && isAle)
            {
                velocityXWall = faceVelocityX[iFace] + uWall;
                velocityYWall = faceVelocityY[iFace] + vWall;
                velocityZWall = faceVelocityZ[iFace] + wWall;
            }

            //! set wall velocity for translating frame and rotating frame.
            if (referenceFrame == TRANSLATIONAL_FRAME || referenceFrame == ROTATIONAL_FRAME)
            {
                int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
                string shroud[100];
                GlobalDataBase::GetData("shroud", &shroud, PHSTRING, nTurboZone);

                //! Parallel
                int iTurboZone = grid->GetOrdinaryGridIndex();
                //! Serial
                if (iTurboZone == -1)
                {
                    iTurboZone = grid->GetZoneID();
                }

                velocityXWall = faceVelocityX[iFace];
                velocityYWall = faceVelocityY[iFace];
                velocityZWall = faceVelocityZ[iFace];

                if (bcRegion->GetBCName() == shroud[iTurboZone])
                {
                    velocityXWall = 0;
                    velocityYWall = 0;
                    velocityZWall = 0;
                }
            }

            qL[IU][jFace] = velocityXWall;
            qL[IV][jFace] = velocityYWall;
            qL[IW][jFace] = velocityZWall;

            qR[IU][jFace] = velocityXWall;
            qR[IV][jFace] = velocityYWall;
            qR[IW][jFace] = velocityZWall;

            if (tWall <= 0.0)
            {
                //! Viscous WALL, adiabatic.
            }
            else
            {
                //! iso-thermal wall.
                RDouble tw = tWall / refDimensionalTemperature;

                for (int m = 0; m < nEquation; ++ m)
                {
                    primitiveVariable[m] = qL[m][jFace];
                }

                RDouble omav = one;
                using namespace GAS_SPACE;
                if (nChemical == 1)
                {
                    omav = gas->ComputeMolecularWeightReciprocal(primitiveVariable);
                }

                RDouble pressureWall = primitiveVariable[IP];
                RDouble rhoWall  = pressureWall / (coefficientOfStateEquation * tw * omav);

                qL[IR][jFace] = rhoWall;
                qR[IR][jFace] = rhoWall;
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

    delete [] primitiveVariable;    primitiveVariable = nullptr;
}

void NSSolverUnstruct::GetInviscidFaceValue(Grid *gridIn, FaceProxy *faceProxy, Limiter * limiter, int localStart, int localEnd)
{
#ifdef USE_CUDA
    Param_NSSolverUnstruct *parameters_gpu = GetControlParameters();
    CallGPUGetQlQr(gridIn, localStart, localEnd, SEGCTION_LENGTH);
    CallGPUGetGamaLR(localStart, localEnd);
    //! PressureFactorLR is not use in GPU code at present. 
    CallGPUReConstructFaceValue(gridIn, limiter, localStart, localEnd, SEGCTION_LENGTH, parameters_gpu);
    CallGPUBoundaryQlQrFix(gridIn, faceProxy, parameters_gpu, localStart, localEnd, SEGCTION_LENGTH);
    return;
#endif

    UnstructGrid *grid = UnstructGridCast(gridIn);
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

    GetQlQr(grid, faceProxy, localStart, localEnd);

    GetGamaLR(grid, faceProxy, localStart, localEnd);

    GetPressureFactorLR(grid, faceProxy, localStart, localEnd);

    ReConstructFaceValue(grid, faceProxy, limiter, localStart, localEnd);

    BoundaryQlQrFix(grid, faceProxy, localStart, localEnd);

    if(ifLowSpeedPrecon)
    {
        if (isUnsteady)
        {
        GetTimeCoefficientLR(grid, faceProxy, localStart, localEnd);
        }
        GetPreconCoefficientLR(grid, faceProxy, localStart, localEnd);
    }
}

//! Face values reconstruction by considering scalar limiters.
void NSSolverUnstruct::ReConstructFaceValue(Grid *gridIn, FaceProxy *faceProxy, Limiter * limiter, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    if(!grid->IsFinestGrid())
    {
        //! Reconstruction is not necessary for coarse grid.
        return;
    }

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int limiterType = parameters->GetLimiterType();
    int isInIniting = GlobalDataBase::GetIntParaFromDB("isInIniting");
    int flowInitMethod = GlobalDataBase::GetIntParaFromDB("flowInitMethod");
    int tscheme     = GlobalDataBase::GetIntParaFromDB("tscheme");
    bool isFirstOrder = (limiterType == ILMT_FIRST) || isInIniting;
    if(isFirstOrder && flowInitMethod != 3)
    {
        return;
    }

    if ( tscheme == GMRES && limiterType == ILMT_VENCAT )
    {
        return;
    }

    //! GMRESnolim GMRESSign
    int *qlsign = faceProxy->GetQLSign();
    int *qrsign = faceProxy->GetQRSign();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    bool isStructuredLimiter = (limiterType >= ILMT_STRUCT);
    if (isStructuredLimiter)
    { 
        //! Limiter for structured limiter.
        ReConstructQlQr_STR(gridIn, faceProxy, nEquation, localStart, localEnd);
        return;
    }

    RDouble **limit2D = reinterpret_cast <RDouble **> (grid->GetDataPtr("limit2D"));;

    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();

    RDouble *qTry = new RDouble[nEquation];

    int reconstructMethod = GlobalDataBase::GetIntParaFromDB("reconmeth");
    const int RECONSTRUCT_USING_RESPECTIVE_LIMITER = 0;

    if (reconstructMethod == RECONSTRUCT_USING_RESPECTIVE_LIMITER)
    {
        //! Q+, Q- reconstructing using their self limiters.
        for (int i = localStart; i < localEnd; ++ i)
        {
            int le,re,j;
            j  = i - localStart;
            le = leftCellofFace[i];
            re = rightCellofFace[i];

            RDouble dx = xfc[i] - xcc[le];
            RDouble dy = yfc[i] - ycc[le];
            RDouble dz = zfc[i] - zcc[le];

            //! GMRESSign
            qlsign[j] = -1;
            qrsign[j] = -1;

            for (int m = 0; m < nEquation; ++ m)
            {
                qTry[m] = qL[m][j];
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                RDouble *limit =  limit2D[m];
                qTry[m] += limit[le] * (dqdx[m][le] * dx + dqdy[m][le] * dy + dqdz[m][le] * dz);
            }

            if (PositiveCheckForDensityAndPressure(qTry))
            {
                //! GMRESSign
                qlsign[j] = 1;
                for (int m = 0; m < nEquation; ++ m)
                {
                    qL[m][j] = qTry[m];
                }
            }

            dx = xfc[i] - xcc[re];
            dy = yfc[i] - ycc[re];
            dz = zfc[i] - zcc[re];

            for (int m = 0; m < nEquation; ++ m)
            {
                qTry[m] = qR[m][j];
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                RDouble *limit =  limit2D[m];
                qTry[m] += limit[re] * (dqdx[m][re] * dx + dqdy[m][re] * dy + dqdz[m][re] * dz);
            }

            if (PositiveCheckForDensityAndPressure(qTry))
            {
                //! GMRESSign
                qrsign[j] = 1;
                for (int m = 0; m < nEquation; ++ m)
                {
                    qR[m][j] = qTry[m];
                }
            }
        }
    }
    else
    {
        //! Q+, Q- use the min limiter coefficients of left and right cell.
        for (int i = localStart; i < localEnd; ++ i)
        {
            int le,re,j;
            j  = i - localStart;
            le = leftCellofFace[i];
            re = rightCellofFace[i];

            //! GMRESSign
            qlsign[j] = -1;
            qrsign[j] = -1;

            RDouble dx = xfc[i] - xcc[le];
            RDouble dy = yfc[i] - ycc[le];
            RDouble dz = zfc[i] - zcc[le];

            for (int m = 0; m < nEquation; ++ m)
            {
                qTry[m] = qL[m][j];
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                RDouble *limit =  limit2D[m];
                qTry[m] += MIN(limit[le], limit[re]) * (dqdx[m][le] * dx + dqdy[m][le] * dy + dqdz[m][le] * dz);
            }

            if (PositiveCheckForDensityAndPressure(qTry))
            {
                //! GMRESSign
                qlsign[j] = 1;

                for (int m = 0; m < nEquation; ++ m)
                {
                    qL[m][j] = qTry[m];
                }
            }

            dx = xfc[i] - xcc[re];
            dy = yfc[i] - ycc[re];
            dz = zfc[i] - zcc[re];

            for (int m = 0; m < nEquation; ++ m)
            {
                qTry[m] = qR[m][j];
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                RDouble *limit =  limit2D[m];
                qTry[m] += MIN(limit[le], limit[re]) * (dqdx[m][re] * dx + dqdy[m][re] * dy + dqdz[m][re] * dz);
            }

            if (PositiveCheckForDensityAndPressure(qTry))
            {
                //! GMRESSign
                qrsign[j] = 1;
                for (int m = 0; m < nEquation; ++ m)
                {
                    qR[m][j] = qTry[m];
                }
            }
        }
    }

    delete [] qTry;    qTry = nullptr;
}

void NSSolverUnstruct::ComputeInviscidFlux(Grid *grid, FaceProxy *faceProxy, int localStart, int localEnd)
{
    int unstructScheme = GlobalDataBase::GetIntParaFromDB("uns_scheme");
    int isInMGIniting = GlobalDataBase::GetIntParaFromDB("isInMGIniting");
    int isInIniting = GlobalDataBase::GetIntParaFromDB("isInIniting");
    int flowInitMethod = GlobalDataBase::GetIntParaFromDB("flowInitMethod");
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

    INVScheme inviscidScheme = GetINVScheme(unstructScheme);

#ifdef USE_GMRESSOLVER
    //1 GMRESInit
    if( 3 == flowInitMethod )
    {
        int originaltscheme       = GlobalDataBase::GetIntParaFromDB("OriginalTscheme");
        int newIterStep           = GlobalDataBase::GetIntParaFromDB("newnstep");
        int flowInitStep          = GlobalDataBase::GetIntParaFromDB("flowInitStep");
        // string uns_scheme_name    = GlobalDataBase::GetStrParaFromDB("uns_scheme_name");

        if( originaltscheme == GMRES && newIterStep <= flowInitStep )
        {
            if( unstructScheme == ISCHEME_GMRES_ROE )
            {
                inviscidScheme = Roe_Scheme_ConservativeForm;
            }
            else if( unstructScheme == ISCHEME_GMRES_Steger )
            {
                inviscidScheme = Steger_Scheme;
            }
            
        }
    }
#endif

    if ((!grid->IsFinestGrid() && isInMGIniting && !ifLowSpeedPrecon) ||
        (grid->IsFinestGrid() && isInIniting   && flowInitMethod != 3 && !ifLowSpeedPrecon))
    {
        //! The precondition has not been considered for vanleer method.
        //! so, only in the compressible multi-grid initialization state, run Vanleer.
        inviscidScheme = Vanleer_Scheme;
    }

    InviscidFluxWrap(grid, faceProxy, inviscidScheme, localStart, localEnd);
}

void NSSolverUnstruct::FillGeomProxy(Grid *gridIn, GeomProxy *geom_proxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *vgn  = grid->GetFaceNormalVelocity();
    RDouble *area = grid->GetFaceArea();

    std::copy(xfn  + localStart, xfn  + localEnd, geom_proxy->GetFaceNormalX());
    std::copy(yfn  + localStart, yfn  + localEnd, geom_proxy->GetFaceNormalY());
    std::copy(zfn  + localStart, zfn  + localEnd, geom_proxy->GetFaceNormalZ());
    std::copy(vgn  + localStart, vgn  + localEnd, geom_proxy->GetFaceVelocity());
    std::copy(area + localStart, area + localEnd, geom_proxy->GetFaceArea()   );
}

void NSSolverUnstruct::GetResidual(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));

#ifdef USE_CUDA
    GPUResCopyBack(grid);
#endif

    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    int nEquation = GetNumberOfEquations();

    PHSPACE::ComputeResidualonGrid(grid, actkey, res, nEquation);
}

void NSSolverUnstruct::InviscidFluxWrap(Grid* gridIn, FaceProxy* faceProxy, INVScheme inviscidScheme, int localStart, int localEnd)
{
    UnstructGrid* grid = UnstructGridCast(gridIn);
    StructGrid* grid0 = StructGridCast(gridIn);

    Param_NSSolverUnstruct* parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();
    int nFaceLength = localEnd - localStart;
    faceProxy->setsize(nFaceLength);

    RDouble** flux = faceProxy->GetFlux();
    GeomProxy* geomProxy = faceProxy->GetGeomProxy();

    if (inviscidScheme == Gks_Scheme)
    {
        GlobalTimeStep(gridIn);
    }

    if (systemgridtype == UNSTRUCTGRID)
    {
        FillGeomProxy(grid, geomProxy, localStart, localEnd);
    }
    else
    {
        //! mixgrid
        //! geomProxy has been set up once for ever.
    }


    //! Inviscid flux scheme parameters proxy.
    InviscidSchemeParameter inviscidSchemeParameterProxy;

#ifdef USE_GMRESSOLVER 
    // GMRESnolim
    int limiterType = parameters->GetLimiterType();
    int limitVariables = parameters->GetLimitVariables();
    int limitVector = parameters->GetLimitVector();
    inviscidSchemeParameterProxy.SetLimiterType(limiterType);
    inviscidSchemeParameterProxy.SetLimitModel(limitVariables);    //! GMRESVenkat
    inviscidSchemeParameterProxy.SetLimitVector(limitVector);

    // GMRESnolim
    RDouble* xfc = grid->GetFaceCenterX();
    RDouble* yfc = grid->GetFaceCenterY();
    RDouble* zfc = grid->GetFaceCenterZ();

    RDouble* xcc = grid->GetCellCenterX();
    RDouble* ycc = grid->GetCellCenterY();
    RDouble* zcc = grid->GetCellCenterZ();

    RDouble** dqdx = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble** dqdy = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble** dqdz = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble* nsx = grid->GetFaceNormalX();
    RDouble* nsy = grid->GetFaceNormalY();
    RDouble* nsz = grid->GetFaceNormalZ();
    RDouble* ns = grid->GetFaceArea();
    RDouble* vol = grid->GetCellVolume();
    RDouble** prim = reinterpret_cast<RDouble**> (grid->GetDataPtr("q"));
    RDouble* gama = reinterpret_cast<RDouble*> (grid->GetDataPtr("gama"));
    RDouble averageVolume = grid->GetAverageVolume();    //! GMRESVenkat

    vector<int>* neighborCells = grid->GMRESGetNeighborCells();
    vector<int>* neighborFaces = grid->GMRESGetNeighborFaces();
    vector<int>* neighborLR = grid->GMRESGetNeighborLR();

    //! GMRES2ndCorrection
    vector<int> BCLeftCells = grid->GetBCLeftCells();
    vector<int> BCRightCells = grid->GetBCRightCells();
    vector<int> BCFaces = grid->GetBCFaces();

    //! GMRES3D
    UnstructBCSet* unstructBCSet = grid->GetUnstructBCSet();

    inviscidSchemeParameterProxy.SetFaceCenter(xfc, yfc, zfc);
    inviscidSchemeParameterProxy.SetCellCenter(xcc, ycc, zcc);
    inviscidSchemeParameterProxy.SetGradient(dqdx, dqdy, dqdz);
    inviscidSchemeParameterProxy.SetFaceNormalAbs(nsx, nsy, nsz);
    inviscidSchemeParameterProxy.SetFaceAreaAbs(ns);
    inviscidSchemeParameterProxy.SetVolume(vol);
    inviscidSchemeParameterProxy.SetPrimitive(prim);
    inviscidSchemeParameterProxy.SetGama(gama);
    inviscidSchemeParameterProxy.SetNeighborCells(neighborCells);
    inviscidSchemeParameterProxy.SetNeighborFaces(neighborFaces);
    inviscidSchemeParameterProxy.SetNeighborLR(neighborLR);

    inviscidSchemeParameterProxy.SetAverageVolume(averageVolume);    //! GMRESVenkat

    //! GMRES2ndCorrection
    inviscidSchemeParameterProxy.SetBCLeftCells(BCLeftCells);
    inviscidSchemeParameterProxy.SetBCRightCells(BCRightCells);
    inviscidSchemeParameterProxy.SetBCFaces(BCFaces);

    //! GMRES3D
    inviscidSchemeParameterProxy.SetUnstructBCSet(unstructBCSet);
#endif
    if (inviscidScheme == Gks_Scheme)
    {
        //int limiterType = parameters->GetLimiterType();
        //int limitVariables = parameters->GetLimitVariables();
        //int limitVector = parameters->GetLimitVector();
        //inviscidSchemeParameterProxy.SetLimiterType(limiterType);
        //inviscidSchemeParameterProxy.SetLimitModel(limitVariables);    //! GMRESVenkat
        //inviscidSchemeParameterProxy.SetLimitVector(limitVector);
        //faceProxy* faceProxy = CreateFaceProxy(grid);
        //faceProxy->SetNSFaceValue(CreateNSFaceValue(grid));
        //GetVisFaceValue(grid, faceProxy, localStart, localEnd);
        //RDouble* gama = reinterpret_cast<RDouble*> (grid->GetDataPtr("gama"));
        //RDouble* xcc = grid->GetCellCenterX();
        //RDouble* ycc = grid->GetCellCenterY();
        //RDouble* zcc = grid->GetCellCenterZ();
        //RDouble* xfc = grid->GetFaceCenterX();
        //RDouble* yfc = grid->GetFaceCenterY();
        //RDouble* zfc = grid->GetFaceCenterZ();
        RDouble** prim = reinterpret_cast<RDouble**> (grid->GetDataPtr("q"));
        int* node_number_of_each_face = grid->GetNodeNumberOfEachFace();  // ÿ����Ķ�����
        // ��ȡ�����������
        /*RDouble* x = grid->GetX();
        RDouble* y = grid->GetY();
        RDouble* z = grid->GetZ();
        int* face2node = grid->GetFace2Node();  // ��-��������ӳ��
        int nTotalFace = grid->GetNTotalFace();
        

        // ��ʼ�����������飨�����ѷ����ڴ棩
        RDouble* xfv = new RDouble[nTotalFace];
        RDouble* yfv = new RDouble[nTotalFace];
        RDouble* zfv = new RDouble[nTotalFace];
        int nodepos = 0;  // ����������face2node�е���ʼλ��
        for (int iFace = 0; iFace < nTotalFace; ++iFace)
        {
            // ��ȡ��ǰ��Ķ�������
            int nNodes = node_number_of_each_face[iFace];
            if (nNodes < 2)  // ������Ҫ2��������ܹ�������
            {
                cerr << "iFace=" << iFace << "ʱ��" << "nNodes < 2" << endl;
                xfv[iFace] = 0.0;
                yfv[iFace] = 0.0;
                zfv[iFace] = 0.0;
                nodepos += nNodes;
                continue;
            }

            // ȡ��ĵ�1���͵�2�����ڶ��㣨����0��1��
            int p1 = face2node[nodepos + 0];  // ��һ����������
            int p2 = face2node[nodepos + 1];  // �ڶ�����������

            // ������������P2 - P1
            xfv[iFace] = x[p2] - x[p1];
            yfv[iFace] = y[p2] - y[p1];
            zfv[iFace] = z[p2] - z[p1];
            //cerr << "xfv[" << iFace << "]=" << xfv[iFace] << "    " << "x[" << p2 << "]=" << x[p2] << "    " << "x[" << p1 << "]=" << x[p1]<<endl;
            //cerr << "yfv[" << iFace << "]=" << yfv[iFace] << "    " << "y[" << p2 << "]=" << y[p2] << "    " << "y[" << p1 << "]=" << y[p1]<<endl;
            //cerr << "zfv[" << iFace << "]=" << zfv[iFace] << "    " << "z[" << p2 << "]=" << z[p2] << "    " << "z[" << p1 << "]=" << z[p1]<<endl;
            // ��λ���������������Ҫ��λ������
            RDouble len = sqrt(xfv[iFace] * xfv[iFace] + yfv[iFace] * yfv[iFace] + zfv[iFace] * zfv[iFace]);
            xfv[iFace] /= len+SMALL;
            yfv[iFace] /= len+SMALL;
            zfv[iFace] /= len+SMALL;
            // ���¶�������λ�ã��ƶ�����һ����Ķ��㣩
            nodepos += nNodes;
        }*/
        //cerr << "xfv[0]=" << xfv[0] << ";yfv[0]=" << yfv[0] << ";zfv[0]=" << zfv[0] << endl;
        //cerr << "xfv[1]=" << xfv[1] << ";yfv[1]=" << yfv[1] << ";zfv[1]=" << zfv[1] << endl;
		RDouble dtmin = GlobalDataBase::GetDoubleParaFromDB("dtmin");
        UnstructBCSet* unstructBCSet = grid->GetUnstructBCSet();
        //inviscidSchemeParameterProxy.SetCellCenter(xcc, ycc, zcc);
        //inviscidSchemeParameterProxy.SetFaceCenter(xfc, yfc, zfc);
        RDouble** dqdx = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradPrimtiveVarX"));
        RDouble** dqdy = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradPrimtiveVarY"));
        RDouble** dqdz = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradPrimtiveVarZ"));
        int* leftcellindexofFace = faceProxy->GetLeftCellIndexOfFace();
        int* rightcellindexofFace = faceProxy->GetRightCellIndexOfFace();
		//int* leftCellOfFace = grid->GetLeftCellOfFace();  // 未使用
		//int* rightCellOfFace = grid->GetRightCellOfFace();
        RDouble* area = geomProxy->GetFaceArea();
		//inviscidSchemeParameterProxy.SetFaceVector(xfv, yfv, zfv, nTotalFace);
        inviscidSchemeParameterProxy.SetGradient(dqdx, dqdy, dqdz);
        inviscidSchemeParameterProxy.SetPrimitive(prim);
        //inviscidSchemeParameterProxy.SetGama(gama);
		inviscidSchemeParameterProxy.SetNodeNumOfEachFace(node_number_of_each_face);
		inviscidSchemeParameterProxy.SetTimeStep(dtmin);
        inviscidSchemeParameterProxy.SetUnstructBCSet(unstructBCSet);
        inviscidSchemeParameterProxy.SetFaceArea(area);
        //delete[] xfv;    xfv = NULL;
        //delete[] yfv;    yfv = NULL;
        //delete[] zfv;    zfv = NULL;
	}

    SetInviscidSchemeParameters(&inviscidSchemeParameterProxy, faceProxy, parameters);

#ifdef USE_CUDA
    CallGPUinviscidScheme(faceProxy, &inviscidSchemeParameterProxy, localStart, localEnd);
    CallGPUInviscidFluxWrapLoop1(nLaminar, localStart, localEnd, SEGCTION_LENGTH);
    return;
#endif

    //! The flux is actually computed here!
    inviscidScheme(faceProxy, &inviscidSchemeParameterProxy);

    RDouble* area = geomProxy->GetFaceArea();
    //! iFace is the global ID, jFace is the local ID in this section.
    if (inviscidScheme != Gks_Scheme)
    {
        //! The flux is computed in the faceProxy, so we need to multiply the area.
        for (int iFace = localStart; iFace < localEnd; ++ iFace)
        {
            int jFace = iFace - localStart;
            RDouble &areas = area[jFace];
            for (int m = 0; m < nLaminar; ++ m)
            {
                flux[m][jFace] *= areas;
            }
		}
    }
}

void NSSolverUnstruct::RegisterCFDSolverInterfaceField()
{
    NSSolver::RegisterCFDSolverInterfaceField();

    NSSolver::RegisterCFDSolverInterfaceFieldUnstruct();
}

void NSSolverUnstruct::RegisterCFDSolverInterpointField()
{
    NSSolver::RegisterCFDSolverInterpointField();
}

void NSSolverUnstruct::RegisterOversetField()
{
    NSSolver::RegisterOversetField();
}

void NSSolverUnstruct::SolutionFix(Grid *gridIn, RDouble *primitiveVariable, int iCell)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    int nEquation = GetNumberOfEquations();

    int **cell2face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    using namespace GAS_SPACE;
    using namespace IDX;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    RDouble roo = primitiveVarFarfield[IR];
    RDouble poo = primitiveVarFarfield[IP];
    RDouble densityLimitLoc = 1.0e-6 * roo;
    RDouble pressureLimitLoc = 1.0e-6 * poo;

    for (int m = 0; m < nEquation; ++ m)
    {
        primitiveVariable[m] = zero;
    }

    int iFace, le, re, neighborCell;
    RDouble volumeSum, ovolumeSum, vsur;
    volumeSum = 0.0;
    for (int icell2face = 0; icell2face < faceNumberOfEachCell[iCell]; ++ icell2face)
    {
        iFace = cell2face[iCell][icell2face];
        le    = leftCellofFace[iFace];
        re    = rightCellofFace[iFace];

        neighborCell = le;
        if (iCell == le) neighborCell = re;

        if (q[IR][neighborCell] < densityLimitLoc || q[IP][neighborCell] < pressureLimitLoc)
        {
            continue;
        }

        vsur = one;

        volumeSum += vsur;

        for (int m = 0; m < nEquation; ++ m)
        {
            RDouble f = q[m][neighborCell];
            if (m == IR || m == IP)
            {
                f = ABS(f);
            }
            primitiveVariable[m] += f * vsur;
        }
    }

    if (volumeSum > zero)
    {
        ovolumeSum = 1.0 / volumeSum;
        for (int m = 0; m < nEquation; ++ m)
        {
            primitiveVariable[m] *= ovolumeSum;
        }
    }
    else
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            primitiveVariable[m] = primitiveVarFarfield[m];
        }
    }
}

void NSSolverUnstruct::UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy)
{
#ifdef USE_CUDA
    Param_NSSolverUnstruct *parameters_gpu = GetControlParameters();
    int nEquation_gpu = GetNumberOfEquations();
    CallGPUUpdateFlowField(gridIn, parameters_gpu, nEquation_gpu);
    return;
#endif
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    RDouble **q  = qProxy->GetField_UNS();
    RDouble **dq = dqProxy->GetField_UNS();

    RDouble **qnew = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *gamma = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble **t    = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();
    int *iBlank = grid->GetBlankIndex();

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    //! Three modes for temperature, including the translational-rotational temperature, vibrational temperature and electron temperature.
    RDouble rm, pm, temperature[3];

    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    RDouble gama;
    using namespace GAS_SPACE;
    using namespace IDX;

    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    RDouble roo = primitiveVarFarfield[IR];
    RDouble poo = primitiveVarFarfield[IP];
    RDouble densityLimit  = 1.0e-6 * roo;
    RDouble pressureLimit = 1.0e-6 * poo;

    RDouble *primitiveVariable = new RDouble[nEquation];
    RDouble *qtry = new RDouble[nEquation];

    int *negativeCellList = new int [nTotalCell];
    PHSPACE::SetField(negativeCellList, 0, nTotalCell);

    int nNegativeCell    = 0;
    int mostNegativeCell = -1;
    RDouble mostNegativePressure = PHSPACE::LARGE;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (iBlank[iCell] != ACTIVE)
        {
            continue;
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            primitiveVariable[m] = q[m][iCell];
        }

        gama = gamma[iCell];
        //! Obtain the temperature.
        temperature[ITT] = t[ITT][iCell];
        if (nTemperatureModel == 1)
        {
            temperature[ITV] = temperature[ITT];
            temperature[ITE] = temperature[ITT];
        }
        else
        {
            temperature[ITV] = t[ITV][iCell];
            temperature[ITE] = temperature[ITV];
            if (nTemperatureModel == 3)
            {
                temperature[ITE] = t[ITE][iCell];
            }
        }
        gas->Primitive2Conservative(primitiveVariable, gama, temperature[ITV], temperature[ITE], qtry);

        for (int m = 0; m < nLaminar; ++ m)
        {
            RDouble dqTemp = dq[m][iCell];
            qtry[m] += dqTemp;
        }

        for (int m = 0; m < nTemperatureModel-1; ++ m)
        {
            qtry[nLaminar + nChemical + m] += dq[nLaminar + nChemical + m][iCell];
        }

        //! As input parameter, temperature is the initial temperature of the next iteration step.
        //! As returned parameter, temperature saves the temperatures of the next time-advancing step.
        gas->Conservative2Primitive(qtry, gama, primitiveVariable, temperature);

        rm = primitiveVariable[IR];
        pm = primitiveVariable[IP];

        if (rm < densityLimit || pm < pressureLimit)
        {
            if (pm < mostNegativePressure)
            {
                negativeCellList[iCell] = 1;
                mostNegativePressure = pm;
                mostNegativeCell     = iCell;
            }

            //cout << "Negative rho or pressure found, level = " << grid->GetLevel() << ", cell = " << iCell << endl;
            //TK_Exit::ExceptionExit("Negative!");

            nNegativeCell += 1;
            SolutionFix(grid, primitiveVariable, iCell);
        }

        if (nChemical == 1)
        {
            gas->NormalizePrimitive(primitiveVariable);
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            qnew[m][iCell] = primitiveVariable[m];
        }
    }

    if (nNegativeCell > 0)
    {
        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

        RDouble negativeRatio = (RDouble)nNegativeCell / (RDouble)nTotalCell;
        cout.setf(ios::scientific);
        cout.precision(4);
        cout << "      Warning: negative pressure appears in " << nNegativeCell << " cells";
        cout << ", ratio = " << negativeRatio;
        cout << ", iter = " << GlobalDataBase::GetIntParaFromDB("outnstep") ;
        cout << ", level = " << grid->GetLevel();
        cout << ", zone  = " << grid->GetZoneID();
        cout << ", cellID = " << mostNegativeCell;
        cout << ", the minimum pressure appears at: (" 
             << xcc[mostNegativeCell] << ", " << ycc[mostNegativeCell] << ", " << zcc[mostNegativeCell] << "). ";

        int nBoundFace = grid->GetNBoundFace();
        int **cell2face = grid->GetCell2Face();
        int isInteriorCell = 1;
        for (int jFace = 0; jFace < faceNumberOfEachCell[mostNegativeCell]; ++ jFace)
        {
            int faceID = cell2face[mostNegativeCell][jFace];
            if(faceID < nBoundFace)
            {
                if(isInteriorCell)
                {
                    cout << endl;
                }

                isInteriorCell = 0;
                UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[faceID]);
                int bcType = bcRegion->GetBCType();
                cout << "            face bcType = " << bcType << endl;
            }
        }

        if(isInteriorCell)
        {
            cout << " is interior cell ..." << endl;
        }
    }

    DelPointer(negativeCellList);
    DelPointer(primitiveVariable);
    DelPointer(qtry);
}

void NSSolverUnstruct::SetGhostDQLUSGS(Grid *gridIn, RDouble **dq, RDouble **q)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *rightCellofFace = grid->GetRightCellOfFace();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    bool isViscous = parameters->IsViscous();

    using namespace PHSPACE;
    using namespace GAS_SPACE; 

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        if (bcType == PHENGLEI::SYMMETRY || (IsWall(bcType) && !isViscous))
        {
            SetSymmBCGhostDQLUSGS(grid, bcRegion, dq, q);
        }
        else if (IsWall(bcType) && isViscous)
        {
            SetWallBCGhostDQLUSGS(grid, bcRegion, dq, q);
        }
        else if (bcType == PHENGLEI::EXTRAPOLATION|| bcType == PHENGLEI::OUTFLOW)
        {
            SetExtraBCGhostDQLUSGS(grid, bcRegion, dq, q);
        }
        else
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                int re = rightCellofFace[iFace];

                //! FIX_ALL set DQ = 0.
                //! Far_FIELD, PROFILE and FIX_V have not been implemented yet.
                for (int m = 0; m < nLaminar; ++ m)
                {
                    dq[m][re] = 0.0;
                }
                for (int m = 0; m < nTemperatureModel - 1; m++)
                {
                    dq[nLaminar + nChemical + m][re] = 0.0;
                }
            }
        }
    }
}

void NSSolverUnstruct::SetSymmBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq, RDouble **q)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    //! Get grid metrics.
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();

    RDouble *gama = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble **t   = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));
    Param_NSSolverUnstruct *parameters = GetControlParameters();    
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble *primitiveVariable = new RDouble[nEquation];
    RDouble *Q = new RDouble[nEquation];

    RDouble gam1,tv = 0.0,te = 0.0;
    RDouble rm, um, vm, wm, pm;
    RDouble rvn, ru, rv, rw, p_new;
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        //! iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            primitiveVariable[m] = q[m][le];
        }

        tv = t[ITT][le];
        te = tv;
        if (nTemperatureModel > 1)    //! Consider Two-Temperature model and Three-Temperature model.
        {
            tv = t[ITV][le];
            te = tv;
            if (nTemperatureModel == 3)
            {
                te = t[ITE][le];
            }
        }
        gas->Primitive2Conservative(primitiveVariable, gama[le], tv, te, Q);

        for (int m = 0; m < nLaminar; ++ m)
        {
            Q[m] += dq[m][le];
        }
        for (int m = 0; m < nTemperatureModel - 1; m++)
        {
            Q[nLaminar + nChemical + m] += dq[nLaminar + nChemical + m][le];
        }

        //! Non-conservative variables in re in level n.
        rm = q[IR][re];
        um = q[IU][re];
        vm = q[IV][re];
        wm = q[IW][re];
        pm = q[IP][re];

        for (int m = 0; m < nLaminar; ++ m)
        {
            dq[m][re] = dq[m][le];
        }
        for (int m = 0; m < nTemperatureModel - 1; m++)
        {
            dq[nLaminar + nChemical + m][re] = dq[nLaminar + nChemical + m][le];
        }

        //! Inviscid wall or symmetry boundary.
        dq[IR][re] = dq[IR][le];
        rvn = 2.0 * (xfn[iFace]*dq[IRU][le] + yfn[iFace]*dq[IRV][le] + zfn[iFace]*dq[IRW][le]);
        dq[IRU][re] = dq[IRU][le] - xfn[iFace] * rvn;
        dq[IRV][re] = dq[IRV][le] - yfn[iFace] * rvn;
        dq[IRW][re] = dq[IRW][le] - zfn[iFace] * rvn;
        rvn   = 2.0 * (xfn[iFace]*Q[IRU] + yfn[iFace]*Q[IRV] + zfn[iFace]*Q[IRW]);
        ru    = Q[IRU] - xfn[iFace] * rvn;
        rv    = Q[IRV] - yfn[iFace] * rvn;
        rw    = Q[IRW] - zfn[iFace] * rvn;
        p_new = Q[IRE] - 0.5 * (Q[IRU]*Q[IRU] + Q[IRV]*Q[IRV] + Q[IRW]*Q[IRW]) / Q[IR];

        gam1 = gama[le] - 1.0;
        dq[IRE][re] = p_new - pm/gam1 + half*((ru*ru + rv*rv + rw*rw)/Q[IR] - rm*(um*um + vm*vm + wm*wm));
    }
    delete [] primitiveVariable;    primitiveVariable = nullptr;
    delete [] Q;    Q = nullptr;
}

void NSSolverUnstruct::SetWallBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq, RDouble **q)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    //! Get grid metrics.

    RDouble *gama = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble **t   = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));
    Param_NSSolverUnstruct *parameters = GetControlParameters();    
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble *primitiveVariable = new RDouble[nEquation];
    RDouble *Q = new RDouble[nEquation];

    RDouble tv = 0.0,te = 0.0;
    RDouble rm, um, vm, wm, pm, ub, vb, wb;

    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        //! iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

        //! Conservative Variables on level n+1.
        for (int m = 0; m < nEquation; ++ m)
        {
            primitiveVariable[m] = q[m][le];
        }
    
        tv = t[ITT][le];
        te = tv;
        if (nTemperatureModel > 1)    //! Consider Two-Temperature model and Three-Temperature model.
        {
            tv = t[ITV][le];
            te = tv;
            if (nTemperatureModel == 3)
            {
                te = t[ITE][le];
            }
        }
        gas->Primitive2Conservative(primitiveVariable, gama[le], tv, te, Q);
    
        for (int m = 0; m < nLaminar; ++ m)
        {
            Q[m] += dq[m][le];
        }
        for (int m = 0; m < nTemperatureModel - 1; m++)
        {
            Q[nLaminar + nChemical + m] += dq[nLaminar + nChemical + m][le];
        }

        //! Non-conservative variables in re in level n.
        rm = q[IR][re];
        um = q[IU][re];
        vm = q[IV][re];
        wm = q[IW][re];
        pm = q[IP][re];

        for (int m = 0; m < nLaminar; ++ m)
        {
            dq[m][re] = dq[m][le];
        }
        for (int m = 0; m < nTemperatureModel - 1; m++)
        {
            dq[nLaminar + nChemical + m][re] = dq[nLaminar + nChemical + m][le];
        }

        ub = 0.0;
        vb = 0.0;
        wb = 0.0;
        dq[IR ][re] =   dq[IR ][le];
        dq[IRU][re] = - dq[IRU][le] + two * ub * dq[IR][le];
        dq[IRV][re] = - dq[IRV][le] + two * vb * dq[IR][le];
        dq[IRW][re] = - dq[IRW][le] + two * wb * dq[IR][le];
        dq[IRE][re] =   dq[IRE][le] + ub * (dq[IRU][re] - dq[IRU][le])
            + vb * (dq[IRV][re] - dq[IRV][le])
            + wb * (dq[IRW][re] - dq[IRW][le]);
    }
    delete [] primitiveVariable;    primitiveVariable = nullptr;
    delete [] Q;    Q = nullptr;
}

void NSSolverUnstruct::SetExtraBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq, RDouble **q)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    //! Get grid metrics.

    RDouble *gama = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble **t   = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));
    Param_NSSolverUnstruct *parameters = GetControlParameters();    
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble *primitiveVariable = new RDouble[nEquation];
    RDouble *Q = new RDouble[nEquation];

    RDouble tv = 0.0,te = 0.0;
    RDouble rm, um, vm, wm, pm;

    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        //! iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

        //! Conservative Variables on level n+1.
        for (int m = 0; m < nEquation; ++ m)
        {
            primitiveVariable[m] = q[m][le];
        }

        tv = t[ITT][le];
        te = tv;
        if (nTemperatureModel > 1)    //! Consider Two-Temperature model and Three-Temperature model.
        {
            tv = t[ITV][le];
            te = tv;
            if (nTemperatureModel == 3)
            {
                te = t[ITE][le];
            }
        }
        gas->Primitive2Conservative(primitiveVariable, gama[le], tv, te, Q);

        for (int m = 0; m < nLaminar; ++ m)
        {
            Q[m] += dq[m][le];
        }
        for (int m = 0; m < nTemperatureModel - 1; m++)
        {
            Q[nLaminar + nChemical + m] += dq[nLaminar + nChemical + m][le];
        }

        //! Non-conservative variables in re in level n.
        rm = q[IR][re];
        um = q[IU][re];
        vm = q[IV][re];
        wm = q[IW][re];
        pm = q[IP][re];

        for (int m = 0; m < nLaminar; ++ m)
        {
            dq[m][re] = dq[m][le];
        }
        for (int m = 0; m < nTemperatureModel - 1; m++)
        {
            dq[nLaminar + nChemical + m][re] = dq[nLaminar + nChemical + m][le];
        }

        for (int m = 0; m < nLaminar; ++ m)
        {
            dq[m][re] =  dq[m][le];
        }
        for (int m = 0; m < nTemperatureModel - 1; m++)
        {
            dq[nLaminar + nChemical + m][re] = dq[nLaminar + nChemical + m][le];
        }
    }
    delete [] primitiveVariable;    primitiveVariable = nullptr;
    delete [] Q;    Q = nullptr;
}

RDouble NSSolverUnstruct::GetRad(RDouble *prim, RDouble *nxyz, RDouble gama)
{
    using namespace IDX;

    RDouble &xfn  = nxyz[0];
    RDouble &yfn  = nxyz[1];
    RDouble &zfn  = nxyz[2];
    RDouble &area = nxyz[3];
    RDouble &vgn  = nxyz[4];

    RDouble &rm = prim[IR];
    RDouble &um = prim[IU];
    RDouble &vm = prim[IV];
    RDouble &wm = prim[IW];
    RDouble &pm = prim[IP];

    RDouble vn = xfn * um + yfn * vm + zfn * wm - vgn;
    RDouble eigen = ABS(vn) + sqrt(gama * pm / rm);
    return eigen * area;
}

void NSSolverUnstruct::SolveLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble **dq = dqProxy->GetField_UNS();

    //! Get grid metrics.
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble *vgn  = grid->GetFaceNormalVelocity();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int *iBlank = grid->GetBlankIndex();

    //! Get flow variables.
    RDouble **q = reinterpret_cast <RDouble **> (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast <RDouble **> (grid->GetDataPtr("t"));
    RDouble *gamma = reinterpret_cast <RDouble *> (grid->GetDataPtr("gama"));
    RDouble *diagonal = reinterpret_cast <RDouble *> (grid->GetDataPtr("diagonal"));
    RDouble **res = reinterpret_cast <RDouble **> (grid->GetDataPtr("res"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int nm = parameters->GetNSEquationNumber();

    int nElectronIndex = parameters->GetIndexOfElectron();

    RDouble refReNumber = parameters->GetRefReNumber();

    int isUnsteady = parameters->GetIsUnsteady();
    int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

    bool isViscous = parameters->IsViscous();

    RDouble *viscousLaminar = nullptr, *viscousTurbulence = nullptr;
    if (isViscous)
    {
        viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
        viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
    }

    //! Some intermediate variables.
    RDouble nxyz[5];
    RDouble *qNeighbor  = new RDouble[nEquation];
    RDouble *dqNeighbor = new RDouble[nEquation];
    RDouble *dqTry = new RDouble[nEquation];
    RDouble *dqOld = new RDouble [nEquation];
    RDouble *rhs0  = new RDouble[nEquation];
    RDouble *df    = new RDouble[nEquation];
    RDouble *facePrimitiveVariable = new RDouble [nEquation];
    RDouble *temperature = new RDouble[nTemperatureModel];
    RDouble *Ux = new RDouble [nEquation];

    RDouble *preconCoefficient = nullptr;
    RDouble *timeCoefficientInverse = nullptr;
    RDouble **preconMatrix = nullptr;
    int nElement = 0;
    RDouble *preconMatrixNeighbor = nullptr;
    RDouble *preconMatrixTemp = nullptr;

    if (ifLowSpeedPrecon != 0)
    {
        preconCoefficient = reinterpret_cast< RDouble * > (grid->GetDataPtr("preconCoefficient"));
        preconMatrix = reinterpret_cast< RDouble ** > (grid->GetDataPtr("preconMatrix"));
        nElement = nEquation * nEquation;
        preconMatrixNeighbor = new RDouble [nElement];
        preconMatrixTemp = new RDouble [nElement];

        if (isUnsteady)
        {
            timeCoefficientInverse = reinterpret_cast< RDouble * > (grid->GetDataPtr("timeCoefficientInverse"));
        }
    }
    using namespace IDX;

    //SpectrumRadius(grid, diagonal);

    int **cell2face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    RDouble **deltaFlux = LUplusDQ->GetField_UNS();

    //! Now the Forward Sweep.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (iBlank[iCell] != ACTIVE)
        {
            for (int m = 0; m < nLaminar; ++ m)
            {
                dq[m][iCell] = 0.0;
            }
            continue;
        }

        for (int m = 0; m < nLaminar; ++ m)
        {
            rhs0[m] = 0.0;

            //! Back up the old dq, to compute the convergence.
            //! it is initialized by zero or res in the first sweep,
            //! then it is updated by backward sweep, that is dq*.
            dqOld[m] = dq[m][iCell];

            //! Here, the deltaFlux is computed in Upper backward sweep!
            //! res: b      deltaFlux: Ux, which is computed in backward sweep.
            //! the dq is not the real dq, now.
            //! Although the 'dq' changed can not the right one, it dosen't matter, since 
            //! the following only using the lower neighbor cells, whose 'dq' has been updated.
            //! It is convenient for precondition transform.
            dq[m][iCell] = res[m][iCell];

            Ux[m] = deltaFlux[m][iCell];

            //! Then reset it to zero to store the Lower forward sweep!
            deltaFlux[m][iCell] = 0.0;
        }
        //! Compute the new Ux (deltaFlux).
        for (int jFace = 0; jFace < faceNumberOfEachCell[iCell]; ++ jFace)
        {
            int face = cell2face[iCell][jFace];
            int le   = leftCellofFace[face];
            int re   = rightCellofFace[face];
            if (iBlank[le] == INTERPOLATION || iBlank[re] == INTERPOLATION)
            {
                continue;
            }
            //! One of le and re must be cell itself.
            if (le > iCell || re > iCell) continue;

            nxyz[0] = xfn [face];
            nxyz[1] = yfn [face];
            nxyz[2] = zfn [face];
            nxyz[3] = area[face];
            nxyz[4] = vgn [face];

            if (re == iCell)
            {
                //! If the re is the iCell itself, then swap with 'le'.
                //! Now its neighboring cell 're' belongs to lower triangular.
                SWAP(le, re);

                //! Swap the normal and the moving face velocity.
                nxyz[0] = - nxyz[0];
                nxyz[1] = - nxyz[1];
                nxyz[2] = - nxyz[2];
                
                nxyz[4] = - nxyz[4];
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                qNeighbor[m] = q[m][re];
            }

            for (int m = 0; m < nLaminar; ++ m)
            {
                dqNeighbor[m] = dq[m][re];
            }

            if (ifLowSpeedPrecon != 0)
            {
                for (int m = 0; m < nElement; ++ m)
                {
                    preconMatrixNeighbor[m] = preconMatrix[m][re];
                }
            }
            //! Compute the ROE averaged face primitive variables.
            RDouble rhoRoe, uRoe, vRoe, wRoe, pRoe;
            RoeAveragedVariables(q[IDX::IR][le], q[IDX::IU][le], q[IDX::IV][le], q[IDX::IW][le], q[IDX::IP][le], gamma[le],
                                 q[IDX::IR][re], q[IDX::IU][re], q[IDX::IV][re], q[IDX::IW][re], q[IDX::IP][re], gamma[re],
                                 rhoRoe, uRoe, vRoe, wRoe, pRoe);

            facePrimitiveVariable[IDX::IR] = rhoRoe;
            facePrimitiveVariable[IDX::IU] = uRoe;
            facePrimitiveVariable[IDX::IV] = vRoe;
            facePrimitiveVariable[IDX::IW] = wRoe;
            facePrimitiveVariable[IDX::IP] = pRoe;

            //! The inviscid Delta Flux.
            RDouble gama = gamma[re];
            if (ifLowSpeedPrecon == 0)
            {
                RDouble rad = GetRad(facePrimitiveVariable, nxyz, gama);

                if (nChemical > 0)
                {
                    for (int m = 0; m < nTemperatureModel; ++ m)
                    {
                        temperature[m] = t[m][re];
                    }
                    GAS_SPACE::ChemicalMXDQ(qNeighbor,temperature ,nxyz[0], nxyz[1], nxyz[2], nxyz[3], nxyz[4], dqNeighbor, df, nm, nLaminar, nTemperatureModel, rad, -1,nElectronIndex);
                }
                else
                {
                    MXDQ_STD(qNeighbor, nxyz[0], nxyz[1], nxyz[2], nxyz[3], nxyz[4], gama, dqNeighbor, df, nm, nLaminar, rad, -1);
                }
            }
            else
            {
                //! Compute the inviscid difference flux (df).
                if (!isUnsteady)
                {
                RDouble preconCoeff = preconCoefficient[re];
                    PreconMatrixFreeDF(qNeighbor, dqNeighbor, nxyz, preconMatrixNeighbor, gama, preconCoeff, df, nm, nLaminar, -1);
                }
                else
                {
                    RDouble preconCoeff = preconCoefficient[re];
                RDouble timeCoeff = timeCoefficientInverse[re];
                PreconMatrixFreeDF(qNeighbor, dqNeighbor, nxyz, preconMatrixNeighbor, gama, preconCoeff, timeCoeff, df, nm, nLaminar, -1);
            }
            }

            if (isViscous)
            {
                //! The viscous Delta Flux.
                RDouble distance;
                distance = ABS(xfn[face] * (xcc[re] - xcc[le])
                             + yfn[face] * (ycc[re] - ycc[le])
                             + zfn[face] * (zcc[re] - zcc[le]));

                RDouble viscousLaminarFace = half * (viscousLaminar[le] + viscousLaminar[re]);
                RDouble viscousTurbulenceFace = half * (viscousTurbulence[le] + viscousTurbulence[re]);
                RDouble density = half * (q[IDX::IR][le] + q[IDX::IR][re]);

                RDouble viscosity = viscousLaminarFace + viscousTurbulenceFace;
                RDouble faceArea  = area[face];
                RDouble viscousSpectrumRadius = 2.0 * viscosity / (density * distance * refReNumber + SMALL);
                viscousSpectrumRadius *= half * faceArea;

                for (int m = 0; m < nLaminar; ++ m)
                {
                    RDouble deltaVisFlux = viscousSpectrumRadius * dqNeighbor[m];
                    df[m] -= deltaVisFlux;
                }
            }

            //! Add Flux together.
            for (int m = 0; m < nLaminar; ++ m)
            {
                rhs0[m] += df[m];
            }
        }

        if(ifLowSpeedPrecon == 1)
        {
            //! Here, the 'dq' is 'res' actually.
            //! The reason why operate on 'dq' other than 'res', is that the 'dq' is a temporary variable.
            //! and the real 'res' is forbidden to be changed.

            for (int m = 0; m < nLaminar; ++ m)
            {
                dqTry[m] = dq[m][iCell];
            }

            for (int m = 0; m < nElement; ++ m)
            {
                preconMatrixTemp[m] = preconMatrix[m][iCell];
            }

            if(!isUnsteady)
            {
                ConservativePreconMatrixVariables(dqTry, preconMatrixTemp, nLaminar);
            }
            else
            {
            RDouble timeCoeffTemp = timeCoefficientInverse[iCell];
            ConservativePreconMatrixVariables(dqTry, preconMatrixTemp, timeCoeffTemp, nLaminar);
            }

            for (int m = 0; m < nLaminar; ++ m)
            {
                dq[m][iCell] = dqTry[m];
            }
        }

        for (int m = 0; m < nLaminar; ++ m)
        {
            //! dq = { (b - Ux) - Lx } / D.
            //! Note: the 'dq' has actually initialized by the residual at the beginning of this function.
            //! the 'rhs0' is the total Delta Flux of the neighbors.
            //! rhs0: Lx    diagonal: D.
            dqTry[m] = (dq[m][iCell] - Ux[m] - rhs0[m]) / diagonal[iCell];

            //! Store the lower forward sweep delta-flux, which will be used in the backward sweep.
            deltaFlux[m][iCell] += rhs0[m];
        }

        for (int m = 0; m < nLaminar; ++ m)
        {
            dq[m][iCell] = dqTry[m];

            sweepNormal += SQR(dq[m][iCell] - dqOld[m]);
        }
    }   //! iCell loop.

    SetGhostDQLUSGS(grid, dq, q);

    delete [] facePrimitiveVariable;    facePrimitiveVariable = nullptr;
    delete [] qNeighbor;    qNeighbor = nullptr;
    delete [] dqNeighbor;    dqNeighbor = nullptr;
    delete [] temperature;    temperature = nullptr;
    delete [] dqTry;    dqTry = nullptr;
    delete [] dqOld;    dqOld = nullptr;
    delete [] rhs0;    rhs0 = nullptr;
    delete [] df;    df = nullptr;
    delete [] Ux;    Ux = nullptr;
    if (ifLowSpeedPrecon != 0)
    {
        delete [] preconMatrixNeighbor;    preconMatrixNeighbor = nullptr;
        delete [] preconMatrixTemp;    preconMatrixTemp = nullptr;
    }
}

void NSSolverUnstruct::SolveLUSGSBackward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    //! Here, the dq is the dq-star, which has been initialized in the forward sweep stage.
    //! that it is solve by the lower-triangle equations.
    RDouble **dq = dqProxy->GetField_UNS();

    using namespace IDX;

    //! Get grid metrics.
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble *vgn  = grid->GetFaceNormalVelocity();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int *iBlank = grid->GetBlankIndex();

    //! Get flow variables.
    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));
    RDouble *gamma = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble *diagonal = reinterpret_cast <RDouble *> (grid->GetDataPtr("diagonal"));
    RDouble **res = reinterpret_cast <RDouble **> (grid->GetDataPtr("res"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int nm = parameters->GetNSEquationNumber();

    int nElectronIndex = parameters->GetIndexOfElectron();

    RDouble refReNumber = parameters->GetRefReNumber();

    int isUnsteady = parameters->GetIsUnsteady();
    int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

    bool isViscous = parameters->IsViscous();
    RDouble *viscousLaminar = nullptr, *viscousTurbulence = nullptr;
    if (isViscous)
    {
        viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
        viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
    }

    //! Some intermediate variables.
    RDouble nxyz[5];
    RDouble *qNeighbor  = new RDouble [nEquation];
    RDouble *dqNeighbor = new RDouble [nEquation];
    RDouble *dqTry = new RDouble [nEquation];
    RDouble *dqOld = new RDouble [nEquation];
    RDouble *rhs0  = new RDouble [nEquation];
    RDouble *df    = new RDouble [nEquation];
    RDouble *facePrimitiveVariable = new RDouble [nEquation];
    RDouble *temperature = new RDouble[nTemperatureModel];
    RDouble *Lx = new RDouble [nEquation];

    RDouble *preconCoefficient = nullptr;
    RDouble *timeCoefficientInverse = nullptr;
    RDouble **preconMatrix = nullptr;
    int nElement = 0;
    RDouble *preconMatrixNeighbor = nullptr;
    RDouble *preconMatrixTemp = nullptr;

    if (ifLowSpeedPrecon != 0)
    {
        preconCoefficient = reinterpret_cast <RDouble *> (grid->GetDataPtr("preconCoefficient"));
        preconMatrix = reinterpret_cast< RDouble ** > (grid->GetDataPtr("preconMatrix"));
        nElement = nEquation * nEquation;
        preconMatrixNeighbor = new RDouble [nElement];
        preconMatrixTemp = new RDouble [nElement];

        if (isUnsteady)
        {
            timeCoefficientInverse = reinterpret_cast< RDouble * > (grid->GetDataPtr("timeCoefficientInverse"));
        }
    }

    int **cell2face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    RDouble **deltaFlux = LUplusDQ->GetField_UNS();

    //! Backward sweep.
    for (int iCell = nTotalCell - 1; iCell >= 0; -- iCell)
    {
        if (iBlank[iCell] != ACTIVE)
        {
            for (int m = 0; m < nLaminar; ++ m)
            {
                dq[m][iCell] = 0.0;
            }
            continue;
        }

        for (int m = 0; m < nLaminar; ++ m)
        {
            rhs0[m] = 0.0;

            //! Back up the old dq, to compute the convergence.
            //! the 'dq' is dq*, which has been updated in forward.
            dqOld[m] = dq[m][iCell];

            //! the dq is not the real dq, now.
            //! it is convenient for precondition transform.
            dq[m][iCell] = res[m][iCell];

            Lx[m] = deltaFlux[m][iCell];

            deltaFlux[m][iCell] = 0.0;
        }

        for (int jFace = 0; jFace < faceNumberOfEachCell[iCell]; ++ jFace)
        {
            int face  = cell2face[iCell][jFace];
            int le    = leftCellofFace[face];
            int re    = rightCellofFace[face];

            if (iBlank[le] == INTERPOLATION || iBlank[re] == INTERPOLATION) continue; 

            //! One of le and re must be cell itself.
            if (le < iCell || re < iCell) continue;

            nxyz[0] = xfn [face];
            nxyz[1] = yfn [face];
            nxyz[2] = zfn [face];
            nxyz[3] = area[face];
            nxyz[4] = vgn [face];

            if (re == iCell)
            {
                //! If the re is the iCell itself, then swap with 'le'.
                //! Now its neighboring cell 're' belongs to upper triangular.
                SWAP(le, re);

                //! Swap the normal and the moving face velocity.
                nxyz[0] = - nxyz[0];
                nxyz[1] = - nxyz[1];
                nxyz[2] = - nxyz[2];

                nxyz[4] = - nxyz[4];
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                qNeighbor [m] = q[m][re];
            }

            for (int m = 0; m < nLaminar; ++ m)
            {
                dqNeighbor[m] = dq[m][re];
            }

            if (ifLowSpeedPrecon != 0)
            {
                for (int m = 0; m < nElement; ++ m)
                {
                    preconMatrixNeighbor[m] = preconMatrix[m][re];
                }
            }

            //! Compute the ROE averaged face primitive variables.
            RDouble rhoRoe, uRoe, vRoe, wRoe, pRoe;
            RoeAveragedVariables(q[IDX::IR][le], q[IDX::IU][le], q[IDX::IV][le], q[IDX::IW][le], q[IDX::IP][le], gamma[le],
                q[IDX::IR][re], q[IDX::IU][re], q[IDX::IV][re], q[IDX::IW][re], q[IDX::IP][re], gamma[re],
                rhoRoe, uRoe, vRoe, wRoe, pRoe);

            facePrimitiveVariable[IDX::IR] = rhoRoe;
            facePrimitiveVariable[IDX::IU] = uRoe;
            facePrimitiveVariable[IDX::IV] = vRoe;
            facePrimitiveVariable[IDX::IW] = wRoe;
            facePrimitiveVariable[IDX::IP] = pRoe;

            RDouble gama = gamma[re];
            if (ifLowSpeedPrecon  == 0 )
            {
                RDouble rad = GetRad(facePrimitiveVariable,nxyz,gama);

                if (nChemical > 0)
                {
                    for (int m = 0; m < nTemperatureModel; ++ m)
                    {
                        temperature[m] = t[m][re];
                    }
                    GAS_SPACE::ChemicalMXDQ(qNeighbor,temperature ,nxyz[0], nxyz[1], nxyz[2], nxyz[3], nxyz[4], dqNeighbor, df, nm, nLaminar, nTemperatureModel, rad, -1,nElectronIndex);
                }
                else
                {
                    MXDQ_STD(qNeighbor, nxyz[0], nxyz[1], nxyz[2], nxyz[3], nxyz[4], gama, dqNeighbor, df, nm, nLaminar, rad, -1);
                }
            }
            else
            {
                if (!isUnsteady)
                {
                RDouble preconCoeff = preconCoefficient[re];
                    PreconMatrixFreeDF(qNeighbor, dqNeighbor, nxyz, preconMatrixNeighbor, gama, preconCoeff, df, nm, nLaminar, -1);
                }
                else
                {
                    RDouble preconCoeff = preconCoefficient[re];
                RDouble timeCoeff = timeCoefficientInverse[re];
                PreconMatrixFreeDF(qNeighbor, dqNeighbor, nxyz, preconMatrixNeighbor, gama, preconCoeff, timeCoeff, df, nm, nLaminar, -1);
                }
            }

            if (isViscous)
            {
                RDouble dist;
                dist = ABS( xfn[face] * (xcc[re] - xcc[le])
                           + yfn[face] * (ycc[re] - ycc[le])
                           + zfn[face] * (zcc[re] - zcc[le]));

                RDouble vis_l = half * (viscousLaminar[le] + viscousLaminar[re]);
                RDouble vis_t = half * (viscousTurbulence[le] + viscousTurbulence[re]);
                RDouble density = half * (q[IDX::IR][le] + q[IDX::IR][re]);

                RDouble viscosity = vis_l + vis_t;
                RDouble faceArea  = area[face];
                RDouble viscousSpectrumRadius = 2.0 * viscosity / (density * dist * refReNumber + SMALL);
                viscousSpectrumRadius *= half * faceArea;

                for (int m = 0; m < nLaminar; ++ m)
                {
                    RDouble deltaVisFlux = viscousSpectrumRadius * dqNeighbor[m];
                    df[m] -= deltaVisFlux;
                }
            }

            //! Add Flux together.
            for (int m = 0; m < nLaminar; ++ m)
            {
                rhs0[m] += df[m];
            }
        }

        if(ifLowSpeedPrecon == 1)
        {
            //! Here, the 'dq' is 'res' actually.
            //! The reason why operate on 'dq' other than 'res', is that the 'dq' is a temporary variable.
            //! and the real 'res' is forbidden to be changed.

            for (int m = 0; m < nLaminar; ++ m)
            {
                dqTry[m] = dq[m][iCell];
            }

            for (int m = 0; m < nElement; ++ m)
            {
                preconMatrixTemp[m] = preconMatrix[m][iCell];
            }

            if(!isUnsteady)
            {
                ConservativePreconMatrixVariables(dqTry, preconMatrixTemp, nLaminar);
            }
            else
            {
            RDouble timeCoeffTemp = timeCoefficientInverse[iCell];
            ConservativePreconMatrixVariables(dqTry, preconMatrixTemp, timeCoeffTemp, nLaminar);
            }

            for (int m = 0; m < nLaminar; ++ m)
            {
                dq[m][iCell] = dqTry[m];
            }
        }

        for (int m = 0; m < nLaminar; ++ m)
        {
            //! Note: the 'dq' has been updated by the forward sweep.
            //! the 'rhs0' is the total Delta Flux of the neighbors.
            //! x = {(b - LX) - Ux} / D.
            //! rhs0: Ux.    diagonal: D.
            dqTry[m] = (dq[m][iCell] - Lx[m] - rhs0[m]) / diagonal[iCell];

            //! Store the upper backward sweep delta-flux, which will be used in the forward sweep.
            deltaFlux[m][iCell] += rhs0[m];
        }

        //RestrictDQ(dqTry, nLaminar, dqlim);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dq[m][iCell] = dqTry[m];

            sweepNormal += SQR(dq[m][iCell] - dqOld[m]);
        }
    }   //! iCell loop.

    delete [] facePrimitiveVariable;    facePrimitiveVariable = nullptr;
    delete [] qNeighbor;    qNeighbor = nullptr;
    delete [] dqNeighbor;    dqNeighbor = nullptr;
    delete [] temperature;    temperature = nullptr;
    delete [] dqTry;    dqTry = nullptr;
    delete [] dqOld;    dqOld = nullptr;
    delete [] rhs0;    rhs0 = nullptr;
    delete [] df;    df = nullptr;
    delete [] Lx;    Lx = nullptr;
    if (ifLowSpeedPrecon != 0)
    {
        delete [] preconMatrixNeighbor;    preconMatrixNeighbor = nullptr;
        delete [] preconMatrixTemp;    preconMatrixTemp = nullptr;
    }
}


void NSSolverUnstruct::ZeroResidualOfSpecialCells(Grid *gridIn)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isOverset = parameters->GetIsOverLapping();
    if (!isOverset)
    {
        return;
    }

    UnstructGrid *grid = PHSPACE::UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    RDouble **residual = reinterpret_cast<RDouble **> (grid->GetDataPtr("res"));
    int nEquation = GetNumberOfEquations();
    int *iBlank = grid->GetBlankIndex();

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (iBlank[iCell] != ACTIVE)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                residual[m][iCell] = 0.0;
            }
        }
    }
    return;
}

void NSSolverUnstruct::SpectrumRadius(Grid *gridIn)
{
    SpectrumRadiusInviscid(gridIn);
    //! If the flow is viscous, it need to count the contribution from viscosity.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isViscous = parameters->GetViscousType();
    if (isViscous)
    {
        SpectrumRadiusViscous(gridIn);
    }

    int nChemical = parameters->GetChemicalFlag();
    int nChemicalRadius = parameters->GetNChemicalRadius();

    if (nChemical == 1 && nChemicalRadius == 1)
    {
        SpectrumRadiusChemical(UnstructGridCast(gridIn));
    }
}

void NSSolverUnstruct::ConservativePreconMatrix(Grid *gridIn, int sign)
{
    UnstructGrid *grid = UnstructGridCast( gridIn );
    int nTotalCell = grid->GetNTotalCell();

    //! Get flow variables.
    RDouble **q = reinterpret_cast <RDouble **> (grid->GetDataPtr("q"));
    RDouble *gamma = reinterpret_cast <RDouble *> (grid->GetDataPtr("gama"));
    RDouble **preconMatrix = reinterpret_cast <RDouble **> (grid->GetDataPtr("preconMatrix"));
    RDouble *preconCoefficient = reinterpret_cast< RDouble * > (grid->GetDataPtr("preconCoefficient"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nNSEquation = parameters->GetNSEquationNumber();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nEquation = nLaminar + nChemical;

    int nElement = nEquation * nEquation;
    RDouble *qPrecondition = new RDouble [nEquation];
    RDouble **MG   = NewPointer2 <RDouble> (nEquation, nEquation);

    using namespace IDX;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble gama = gamma[iCell];
        RDouble preconCoeff = preconCoefficient[iCell];

        for (int m = 0; m < nEquation; ++ m)
        {
            qPrecondition[m] = q[m][iCell];
        }

        ComputeConservativePreconMatrix(sign, MG, qPrecondition, gama, preconCoeff, nNSEquation, nEquation);

        for (int row = 0; row < nEquation; ++ row)
        {
            for (int col = 0; col < nEquation; ++ col)
            {
                int index = row * nEquation + col;
                preconMatrix[index][iCell] = MG[row][col];
            }
        }
    }

    grid->SetGhostCell(preconMatrix,nElement);

    delete [] qPrecondition;    qPrecondition = nullptr;
    delete [] MG;    MG = nullptr;
}

void NSSolverUnstruct::Diagonal(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble *vol = grid->GetCellVolume();

    RDouble *dt = reinterpret_cast<RDouble *> (grid->GetDataPtr("dt"));
    RDouble *invSpectralRadius = reinterpret_cast<RDouble *> (grid->GetDataPtr("invSpectralRadius"));
    RDouble *visSpectralRadius = reinterpret_cast<RDouble *> (grid->GetDataPtr("visSpectralRadius"));
    RDouble *diagonal = reinterpret_cast <RDouble *> (grid->GetDataPtr("diagonal"));

    using namespace IDX;

    //! If flow is unsteady, need to count the contribution of unsteady.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();
    if (!isUnsteady)
    {
        //! Part 1: Diagonal = volume / dt.
        //! Notice: the dt[iCell] has divided vol, so the diagonal[iCell] has multiplied the vol by default.
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            diagonal[iCell] = 1.0 / dt[iCell];    //! Notice: the dt[iCell] has divided vol, so the diagonal[iCell] has multiplied the vol by default.
        }
    }
    else
    {
        RDouble dualTimeCoefficient[7];
        const int methodOfDualTime = GlobalDataBase::GetIntParaFromDB("methodOfDualTime");
        ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);
        RDouble dualTimeSpectrumC1 = - dualTimeCoefficient[3];
        RDouble dualTimeSpectrumC2 =   dualTimeCoefficient[6];

        if (ifLowSpeedPrecon)
        {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
                diagonal[iCell] = dualTimeSpectrumC2 / dt[iCell];
            }
        }
        else
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                diagonal[iCell] = dualTimeSpectrumC1 * vol[iCell] + (dualTimeSpectrumC2 / dt[iCell]);    //! Notice: the dt[iCell] has divided vol, as same as above.
            }
        }
    }

    //! Part 2: Diagonal += (|V.n| + c). (|V.n| + c) is inviscid spectral radius.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        diagonal[iCell] += invSpectralRadius[iCell];
    }

    //! If the flow is viscous, it need to count the contribution from viscosity.
    bool isViscous = parameters->IsViscous();
    if (isViscous)
    {
        //! Part 3: Diagonal += 2*miu/(rho.|n.dr|). 2*miu/(rho.|n.dr|) is viscous spectral radius.
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            diagonal[iCell] += visSpectralRadius[iCell];
        }
    }

    int nm        = parameters->GetNSEquationNumber();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nChemicalRadius = parameters->GetNChemicalRadius();
    int nChemicalSource  = parameters->GetNChemicalSource();

    RDouble rad_chem;

    if (nChemical == 1 && nChemicalSource == 1 && nChemicalRadius == 1)
    {
        //! Part 4: Diagonal += chemical radius.
        RDouble **chemSpectralRadius = reinterpret_cast<RDouble **> (grid->GetDataPtr("chemSpectralRadius"));
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            rad_chem = zero;
            for (int m = nm; m < nLaminar; ++ m)
            {
                rad_chem = MAX(nChemicalRadius * chemSpectralRadius[ m - nm ][ iCell ], rad_chem);
            }
            diagonal[iCell] += rad_chem;    //! Referring to the structure part, the dt[iCell] contained in the diagonal[iCell] can be improved for more stability here.
        }
    }
    //grid->SetGhostCell(diagonal, nLaminar);
}

void NSSolverUnstruct::RestrictDQ(RDouble *dqtry, int nm, RDouble lim_max)
{
    for (int m = 0; m < nm; ++ m)
    {
        dqtry[m] = MAX(MIN(dqtry[m], lim_max),-lim_max);
    }
}

void NSSolverUnstruct::ModifyDensity(RDouble *dq, int nNSEquation,int nLaminar,int nEquation,RDouble &densityRelax)
{
    int nChemcalSourceModified = gas->GetnChemcalSourceModified();
    if(nChemcalSourceModified==0)
    {
        for (int m = 1; m < nNSEquation; ++ m)
        {
            dq[m] = 1.0;
        }
        for (int m = nNSEquation; m <= nLaminar; ++ m)
        {
            dq[m] = densityRelax;
        }
        for (int m = nLaminar+1; m < nEquation; ++ m)
        {
            dq[m] = 1.0;
        }
    }
    else
    {
        for (int m = 1; m < nNSEquation - 1; ++ m)
        {
            dq[m] = 1.0;
        }
        for (int m = nNSEquation - 1; m < nEquation; ++ m)
        {
            dq[m] = densityRelax;
        }
    }
}

void NSSolverUnstruct::AirForceCoef(ActionKey *actkey)
{
#ifdef USE_CUDA
    // conduct after main loop
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));
    CallGPUFlowVariablesCopyBack(grid);
#endif
    const int ALL_PART = -1;
    AirForceCoefParts(actkey, ALL_PART);

    int Coordinate = LocalCoordinate;
    uint_t nTBC = GlobalBoundaryCondition::GetNumberOfBC();
    for (int iPart = 0; iPart < nTBC; ++iPart)
    {
        if (!GlobalBoundaryCondition::IsSolidWallBC(iPart)) continue;

        //! Compute local part force.
        AirForceCoefParts(actkey, iPart, Coordinate);

        //! Compute global part force.
        AirForceCoefParts(actkey, iPart);
    }
}

void NSSolverUnstruct::AirForceCoefParts(ActionKey *actkey, int partID, int Coordinate)
{
    int level = actkey->level;

    using namespace PHMPI;
    //fstream file;

    UnstructGrid * grid = UnstructGridCast(GetGrid(level));

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *visl = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    int viscousType = parameters->GetViscousType();

    const int ALL_PART = -1;
    SimpleBC *boundaryCondition = nullptr;
    string partBCName = "";
    RDouble TorqueRefX = 0.0, TorqueRefY = 0.0, TorqueRefZ = 0.0;

    //! Used to compute hinge moment:
    int dumpHingeMoment = 0;
    RDouble localCoordAxis0[3];
    RDouble localCoordAxis1[3];
    RDouble hingeMoment = 0.0;

    if(partID == ALL_PART)
    {
        TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
        TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
        TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");
    }
    else if ((partID != ALL_PART) && (Coordinate == GlobalCoordinate))
    {
        boundaryCondition = GlobalBoundaryCondition::GetBC(partID);
        partBCName = boundaryCondition->GetBCName();

        TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
        TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
        TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");
    }
    else if ((partID != ALL_PART) && (Coordinate == LocalCoordinate))
    {
        boundaryCondition = GlobalBoundaryCondition::GetBC(partID);
        partBCName = boundaryCondition->GetBCName();

        boundaryCondition->GetParamData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
        boundaryCondition->GetParamData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
        boundaryCondition->GetParamData("TorqueRefZ", &TorqueRefZ, PHDOUBLE, 1);
        
        boundaryCondition->GetParamData("dumpHingeMoment", &dumpHingeMoment, PHINT, 1);
        boundaryCondition->GetParamData("localCoordAxis0", localCoordAxis0, PHDOUBLE, 3);
        boundaryCondition->GetParamData("localCoordAxis1", localCoordAxis1, PHDOUBLE, 3);
    }

    int vis_run = 0;
    if (viscousType != INVISCID)
    {
        vis_run = 1;
        //vis_run = 0;
    }

    using namespace IDX;
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    RDouble poo = primitiveVarFarfield[IP];

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea() ;

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble cpx,cpy,cpz;
    cpx = 0.0;
    cpy = 0.0;
    cpz = 0.0;

    RDouble CA_f = 0.0;
    RDouble CA_p = 0.0;
    RDouble CN_f = 0.0;
    RDouble CN_p = 0.0;
    RDouble CZ_f = 0.0;
    RDouble CZ_p = 0.0;
    RDouble Cl_f = 0.0;
    RDouble Cl_p = 0.0;
    RDouble Cn_f = 0.0;
    RDouble Cn_p = 0.0;
    RDouble Cm_f = 0.0;
    RDouble Cm_p = 0.0;

    RDouble pw, cp;
    RDouble dfx, dfy, dfz, dpx, dpy, dpz, fvsx, fvsy, fvsz;
    RDouble nx, ny, nz;
    RDouble vis, txx, txy, txz, tyx, tyy, tyz, tzx, tzy, tzz;
    RDouble dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz;
    RDouble xc, yc, zc, dx, dy, dz, ods;

    RDouble **gradPrimtiveVarX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradPrimtiveVarY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradPrimtiveVarZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    using namespace IDX;

    //! Get the bcRegion information of unstruct grid.
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    //! Cycle all bcregions and find the solid surface type of boundary.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();

        //! If the bc type is wall, calculate the wall coefficient.
        if(bcType == PHENGLEI::SOLID_SURFACE)
        {
            if (partID != ALL_PART && bcName != partBCName)
            {
                continue;
            }
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            int iFace, le, re;
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                iFace = *iter;
                le = leftCellofFace [ iFace ];
                re = rightCellofFace[ iFace ];

                nx  = xfn[iFace] * area[iFace];
                ny  = yfn[iFace] * area[iFace];
                nz  = zfn[iFace] * area[iFace];

                dx  = xfc[iFace] - xcc[le];
                dy  = yfc[iFace] - ycc[le];
                dz  = zfc[iFace] - zcc[le];

                //! Pressure drag
                pw = q[IP][le];
                cp    = two * (pw - poo);

                dpx = nx * cp;
                dpy = ny * cp;
                dpz = nz * cp;

                CA_p += dpx;
                CN_p += dpy;
                CZ_p += dpz;

                dfx = dpx;
                dfy = dpy;
                dfz = dpz;

                xc = xfc[iFace];
                yc = yfc[iFace];
                zc = zfc[iFace];

                if (viscousType > INVISCID)
                {
                    dudx = gradPrimtiveVarX[IU][le];
                    dudy = gradPrimtiveVarY[IU][le];
                    dudz = gradPrimtiveVarZ[IU][le];

                    dvdx = gradPrimtiveVarX[IV][le];
                    dvdy = gradPrimtiveVarY[IV][le];
                    dvdz = gradPrimtiveVarZ[IV][le];

                    dwdx = gradPrimtiveVarX[IW][le];
                    dwdy = gradPrimtiveVarY[IW][le];
                    dwdz = gradPrimtiveVarZ[IW][le];

                    //! Gradient correction.
                    dx  = xcc[re] - xcc[le];
                    dy  = ycc[re] - ycc[le];
                    dz  = zcc[re] - zcc[le];
                    ods = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
                    dx  = dx * ods;
                    dy  = dy * ods;
                    dz  = dz * ods;

                    CorrectGradient(q[IU][le],q[IU][re],dudx,dudy,dudz,dx,dy,dz,ods);
                    CorrectGradient(q[IV][le],q[IV][re],dvdx,dvdy,dvdz,dx,dy,dz,ods);
                    CorrectGradient(q[IW][le],q[IW][re],dwdx,dwdy,dwdz,dx,dy,dz,ods);

                    //! It is unnecessary to consider because the turbulence viscosity coefficient is 0.
                    vis = visl[le];

                    const RDouble c23 = 2.0 / 3.0;
                    txx = vis * c23 * (2.0 * dudx - dvdy - dwdz);
                    tyy = vis * c23 * (2.0 * dvdy - dwdz - dudx);
                    tzz = vis * c23 * (2.0 * dwdz - dudx - dvdy);
                    txy = vis * (dudy + dvdx);
                    txz = vis * (dudz + dwdx);
                    tyz = vis * (dvdz + dwdy);
                    tyx = txy;
                    tzx = txz;
                    tzy = tyz;

                    fvsx = - two * (nx * txx + ny * tyx + nz * tzx) * oRefReNumber;
                    fvsy = - two * (nx * txy + ny * tyy + nz * tzy) * oRefReNumber;
                    fvsz = - two * (nx * txz + ny * tyz + nz * tzz) * oRefReNumber;

                    CA_f += fvsx;
                    CN_f += fvsy;
                    CZ_f += fvsz;

                    dfx += fvsx;
                    dfy += fvsy;
                    dfz += fvsz;

                    Cl_f += (yc - TorqueRefY) * fvsz - (zc - TorqueRefZ) * fvsy;
                    Cn_f += (zc - TorqueRefZ) * fvsx - (xc - TorqueRefX) * fvsz;
                    Cm_f += (xc - TorqueRefX) * fvsy - (yc - TorqueRefY) * fvsx;
                }

                Cl_p += (yc - TorqueRefY) * dpz - (zc - TorqueRefZ) * dpy;
                Cn_p += (zc - TorqueRefZ) * dpx - (xc - TorqueRefX) * dpz;
                Cm_p += (xc - TorqueRefX) * dpy - (yc - TorqueRefY) * dpx;

                cpx += dpx;
                cpy += dpy;
                cpz += dpz;

                if(dumpHingeMoment)
                {
                    Post_ForceMoment forceMoment;
                    RDouble point[3] = {xfc[iFace], yfc[iFace], zfc[iFace]};
                    RDouble faceForce[3] = {dfx, dfy, dfz};
                    hingeMoment += forceMoment.ComputeMoment(point, localCoordAxis0, localCoordAxis1, faceForce);
                }
            }
        }
    }

    DataContainer *cdata = actkey->GetData();

    cdata->Write(&CA_f,sizeof(RDouble));
    cdata->Write(&CA_p,sizeof(RDouble));
    cdata->Write(&CN_f,sizeof(RDouble));
    cdata->Write(&CN_p,sizeof(RDouble));
    cdata->Write(&CZ_f,sizeof(RDouble));
    cdata->Write(&CZ_p,sizeof(RDouble));

    cdata->Write(&cpx,sizeof(RDouble));
    cdata->Write(&cpy,sizeof(RDouble));
    cdata->Write(&cpz,sizeof(RDouble));

    cdata->Write(&Cl_f,sizeof(RDouble));
    cdata->Write(&Cl_p,sizeof(RDouble));
    cdata->Write(&Cn_f,sizeof(RDouble));
    cdata->Write(&Cn_p,sizeof(RDouble));
    cdata->Write(&Cm_f,sizeof(RDouble));
    cdata->Write(&Cm_p,sizeof(RDouble));

    cdata->Write(&hingeMoment, sizeof(RDouble));
}

void NSSolverUnstruct::AirForceCoefBodies(ActionKey *actkey)
{
    int level = actkey->level;

    using namespace PHMPI;
    //fstream file;

    UnstructGrid *grid = UnstructGridCast(GetGrid(level));

    RDouble** q = reinterpret_cast<RDouble**> (grid->GetDataPtr("q"));
    RDouble* visl = reinterpret_cast<RDouble*> (grid->GetDataPtr("visl"));
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    int viscousType = parameters->GetViscousType();

    RDouble TorqueRefXGlobal = 0.0, TorqueRefYGlobal = 0.0, TorqueRefZGlobal = 0.0;

    TorqueRefXGlobal = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    TorqueRefYGlobal = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    TorqueRefZGlobal = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

    //! Used to compute hinge moment.

    int vis_run = 0;
    if (viscousType != INVISCID)
    {
        vis_run = 1;
    }

    using namespace IDX;
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    RDouble poo = primitiveVarFarfield[IP];

    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble cpx, cpy, cpz;
    cpx = 0.0;
    cpy = 0.0;
    cpz = 0.0;

    RDouble pw, cp;
    RDouble dfx, dfy, dfz, dpx, dpy, dpz, fvsx, fvsy, fvsz;
    RDouble nx, ny, nz;
    RDouble vis, txx, txy, txz, tyx, tyy, tyz, tzx, tzy, tzz;
    RDouble dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz;
    RDouble xc, yc, zc, dx, dy, dz, ods;

    RDouble **gradPrimtiveVarX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradPrimtiveVarY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradPrimtiveVarZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    using namespace IDX;

    //! Get the BCRegion information of unstruct grid.
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    int numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    vector< BasicAerodynamicForce * > *basicAerodynamicForceVector = CreateBasicAerodynamicForceVector();

    //! Cycle all bcregions and find the solid surface type of boundary.
    for (int iBody = 0; iBody < numberOfMovingBodies; ++ iBody)
    {
        for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegion; ++ iBCRegionUnstruct)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
            int bcType = bcRegion->GetBCType();

            BasicAerodynamicForce tmpBasicAerodynamicForce;
            //! if the bc type is wall, calculate the wall coefficient.
            if (bcType == PHENGLEI::SOLID_SURFACE)
            {
                string bodyName = bcRegion->GetBodyName();
                if (FaceInAerodynamicForceBodies(bodyName, iBody))
                {
                    tmpBasicAerodynamicForce.SetMomentReferencePoint(TorqueRefXGlobal, TorqueRefYGlobal, TorqueRefZGlobal);

                    vector<int> *faceIndex = bcRegion->GetFaceIndex();
                    int iFace, le, re;
                    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); ++iter)
                    {
                        iFace = *iter;
                        le = leftCellofFace[iFace];
                        re = rightCellofFace[iFace];

                        nx = xfn[iFace] * area[iFace];
                        ny = yfn[iFace] * area[iFace];
                        nz = zfn[iFace] * area[iFace];

                        dx = xfc[iFace] - xcc[le];
                        dy = yfc[iFace] - ycc[le];
                        dz = zfc[iFace] - zcc[le];

                        //! pressure drag
                        pw = q[IP][le];
                        cp = two * (pw - poo);

                        dpx = nx * cp;
                        dpy = ny * cp;
                        dpz = nz * cp;

                        tmpBasicAerodynamicForce.SetPressureAerodynamicForce(dpx, dpy, dpz);

                        dfx = dpx;
                        dfy = dpy;
                        dfz = dpz;

                        xc = xfc[iFace];
                        yc = yfc[iFace];
                        zc = zfc[iFace];

                        tmpBasicAerodynamicForce.SetPressureAerodynamicMoment(xc, yc, zc);

                        if (viscousType > INVISCID)
                        {
                            dudx = gradPrimtiveVarX[IU][le];
                            dudy = gradPrimtiveVarY[IU][le];
                            dudz = gradPrimtiveVarZ[IU][le];

                            dvdx = gradPrimtiveVarX[IV][le];
                            dvdy = gradPrimtiveVarY[IV][le];
                            dvdz = gradPrimtiveVarZ[IV][le];

                            dwdx = gradPrimtiveVarX[IW][le];
                            dwdy = gradPrimtiveVarY[IW][le];
                            dwdz = gradPrimtiveVarZ[IW][le];

                            //! Gradient correction.
                            dx = xcc[re] - xcc[le];
                            dy = ycc[re] - ycc[le];
                            dz = zcc[re] - zcc[le];
                            ods = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
                            dx = dx * ods;
                            dy = dy * ods;
                            dz = dz * ods;

                            CorrectGradient(q[IU][le], q[IU][re], dudx, dudy, dudz, dx, dy, dz, ods);
                            CorrectGradient(q[IV][le], q[IV][re], dvdx, dvdy, dvdz, dx, dy, dz, ods);
                            CorrectGradient(q[IW][le], q[IW][re], dwdx, dwdy, dwdz, dx, dy, dz, ods);

                            //! It is unnecessary to consider because the turbulence viscosity coefficient is 0.
                            vis = visl[le];

                            const RDouble c23 = 2.0 / 3.0;
                            txx = vis * c23 * (2.0 * dudx - dvdy - dwdz);
                            tyy = vis * c23 * (2.0 * dvdy - dwdz - dudx);
                            tzz = vis * c23 * (2.0 * dwdz - dudx - dvdy);
                            txy = vis * (dudy + dvdx);
                            txz = vis * (dudz + dwdx);
                            tyz = vis * (dvdz + dwdy);
                            tyx = txy;
                            tzx = txz;
                            tzy = tyz;

                            fvsx = -two * (nx * txx + ny * tyx + nz * tzx) * oRefReNumber;
                            fvsy = -two * (nx * txy + ny * tyy + nz * tzy) * oRefReNumber;
                            fvsz = -two * (nx * txz + ny * tyz + nz * tzz) * oRefReNumber;

                            tmpBasicAerodynamicForce.SetViscousAerodynamicForce(fvsx, fvsy, fvsz);

                            tmpBasicAerodynamicForce.SetViscousAerodynamicMoment(xc, yc, zc);
                        }
                        tmpBasicAerodynamicForce.ComputeResultantAerodynamicForce();

                        tmpBasicAerodynamicForce.ComputeResultantAerodynamicMoment();

                        (*basicAerodynamicForceVector)[iBody]->AddAerodynamicForce(&tmpBasicAerodynamicForce);
                    }
                }
            }
        }
    }

    PHSPACE::WriteAerodynamicForceToActionKey(actkey, basicAerodynamicForceVector);
    PHSPACE::FreeBasicAerodynamicForceVector(basicAerodynamicForceVector);
}

void NSSolverUnstruct::CpDistriCoef(ActionKey *actkey)
{
    int level = actkey->level;
    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();

    UnstructGrid *grid = UnstructGridCast(GetGrid(level));

    int oriGridIndex = grid->GetOrdinaryGridIndex();
    PHWrite(cdata, oriGridIndex);

    int GridID = grid->GetGridID()->GetIndex();
    PHWrite(cdata, GridID);

    int nsolid_surface = grid->GetNumberOfWallCell();
    if ((nsolid_surface == 0) || (level != 0))
    {
        return;
    }

    using namespace PHMPI;
    using namespace PHENGLEI;
    using namespace IDX;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *visl = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));
    RDouble *vist = reinterpret_cast< RDouble* >  (grid->GetDataPtr("vist"));

    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    RDouble poo = primitiveVarFarfield[IP];

    RDouble oRefReNumber = parameters->GetoRefReNumber();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble refMachNumber = parameters->GetRefMachNumber();

    int wallFunctionType = GlobalDataBase::GetIntParaFromDB("wallFunctionType");
    int viscousType = parameters->GetViscousType();

    //! To obtain the reference values.
    RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
    RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
    RDouble refVelocity = refDimensionalSonicSpeed * refMachNumber;
    RDouble refTotalEnthalpy = gas->ComputeReferenceTotalEnthalpy();

    RDouble TorqueRefX, TorqueRefY, TorqueRefZ;
    TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea() ;

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    int *node_number_of_each_face  = grid->GetNodeNumberOfEachFace();
    long long int *nodePosi  = grid->GetFace2NodeSubscript();
    int *face2node = grid->GetFace2Node();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    int nTotalCell = grid->GetNTotalCell();
    int nTotalNode = grid->GetNTotalNode();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble cfx, cfy, cfz, cpx, cpy, cpz, cmx, cmy, cmz;
    cfx = 0.0;
    cfy = 0.0;
    cfz = 0.0;
    cpx = 0.0;
    cpy = 0.0;
    cpz = 0.0;
    cmx = 0.0;
    cmy = 0.0;
    cmz = 0.0;

    RDouble AoA = 0, attack = 0, angleSlide = 0, angleSlided = 0;
    AoA = parameters->GetAoA();
    attack   = AoA * PI / 180.0;
    angleSlided = parameters->GetAngleOfSlide();
    angleSlide = angleSlided * PI / 180.0;

    RDouble sina = 0, cosa = 0, sinb = 0, cosb = 0;
    sina = sin(attack);
    cosa = cos(attack);
    sinb = sin(angleSlide);
    cosb = cos(angleSlide);

    RDouble **gradPrimtiveVarX = reinterpret_cast< RDouble ** > (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradPrimtiveVarY = reinterpret_cast< RDouble ** > (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradPrimtiveVarZ = reinterpret_cast< RDouble ** > (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble pw = 0, cp = 0, pw_Dim = 0;
    RDouble dfx = 0, dfy = 0, dfz = 0, dpx = 0, dpy = 0, dpz = 0, dmx = 0, dmy = 0, dmz = 0, fvsx = 0, fvsy = 0, fvsz = 0;
    RDouble nx = 0, ny = 0, nz = 0;
    RDouble normalComponetX = 0, normalComponetY = 0, normalComponetZ = 0;
    RDouble visLaminar = 0.0, visTurbulence = 0.0, txx = 0, txy = 0, txz = 0, tyx = 0, tyy = 0, tyz = 0, tzx = 0, tzy = 0, tzz = 0;
    RDouble dudx = 0, dudy = 0, dudz = 0, dvdx = 0, dvdy = 0, dvdz = 0, dwdx = 0, dwdy = 0, dwdz = 0;
    RDouble xc = 0, yc = 0, zc = 0, dx = 0, dy = 0, dz = 0, ods = 0;
    RDouble StantonNumber = 0.0, wallTotalEnthalpy = 0.0, deltaH = 0.0, absVelocity = 0.0;

    int nNSEquation = parameters->GetNSEquationNumber();
    RDouble *wallVariables = new RDouble[nNSEquation];

    RDouble *wallDistance = grid->GetWallDist();

    FluidParameter refParam;
    gas->GetReferenceParameters(refParam);
    RDouble refVicosity = refParam.GetViscosity();
    RDouble refDimensionalMolecular = refParam.GetAverageMolecularWeight();

    RDouble wallTemperature = parameters->GetWallTemperature();
    RDouble knudsenNumber = 0.0, meanFreePath = 0.0, massReciprocal = 1.0;
    RDouble knLength = parameters->GetCharacteristicLength();

    int *visualVariablesType = postVisualWall->GetVisualVariablesType();
    int nWallVariables = postVisualWall->GetVisualVariablesNumber();

    RDouble visualVariables[50];    //! Store the values of variables.
    //! the order of variables in the computational array, eg. 0 denotes the coefficient of pressure,
    //! and the value of array indicates the index of the variable in the computation array.
    int variablesOrder[50];
    for (int m = 0; m < 50; ++m)
    {
        variablesOrder[m] = m;
    }

    RDouble *faceCP    = new RDouble[nBoundFace]();
    RDouble *faceCF    = new RDouble[nBoundFace]();
    RDouble *faceQ     = new RDouble[nBoundFace]();
    RDouble *faceYPlus = new RDouble[nBoundFace]();
    RDouble *facePW    = new RDouble[nBoundFace]();
    RDouble *faceST    = new RDouble[nBoundFace]();
    RDouble *faceKN    = new RDouble[nBoundFace]();

    //! Compute cp and cf on wall face.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (IsWall(bcType))
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                int le = leftCellofFace[iFace];
                int re = rightCellofFace[iFace];

                normalComponetX = xfn[iFace];
                normalComponetY = yfn[iFace];
                normalComponetZ = zfn[iFace];

                nx = normalComponetX * area[iFace];
                ny = normalComponetY * area[iFace];
                nz = normalComponetZ * area[iFace];

                //! Pressure drag.
                dx = xfc[iFace] - xcc[le];
                dy = yfc[iFace] - ycc[le];
                dz = zfc[iFace] - zcc[le];
                pw = q[IP][le];
                cp = two * (pw - poo);

                //! Change the nondimensional pressure into dimensional pressure(Pa).
                pw_Dim = pw * refDimensionalDensity * pow(refVelocity, 2);

                //! Note that the definition of cp is right here and now, but dpx and dfx have multiplied by area.
                dpx = nx * cp;
                dpy = ny * cp;
                dpz = nz * cp;

                dfx = dpx;
                dfy = dpy;
                dfz = dpz;

                fvsx = 0.0;
                fvsy = 0.0;
                fvsz = 0.0;
                RDouble heatCoeff = 0.0, vis = 0.0;

                if (viscousType > INVISCID)
                {
                    dudx = gradPrimtiveVarX[IU][le];
                    dudy = gradPrimtiveVarY[IU][le];
                    dudz = gradPrimtiveVarZ[IU][le];

                    dvdx = gradPrimtiveVarX[IV][le];
                    dvdy = gradPrimtiveVarY[IV][le];
                    dvdz = gradPrimtiveVarZ[IV][le];

                    dwdx = gradPrimtiveVarX[IW][le];
                    dwdy = gradPrimtiveVarY[IW][le];
                    dwdz = gradPrimtiveVarZ[IW][le];

                    //! Gradient correction.
                    dx  = xcc[re] - xcc[le];
                    dy  = ycc[re] - ycc[le];
                    dz  = zcc[re] - zcc[le];
                    ods = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
                    dx  = dx * ods;
                    dy  = dy * ods;
                    dz  = dz * ods;

                    CorrectGradient(q[IU][le],q[IU][re],dudx,dudy,dudz,dx,dy,dz,ods);
                    CorrectGradient(q[IV][le],q[IV][re],dvdx,dvdy,dvdz,dx,dy,dz,ods);
                    CorrectGradient(q[IW][le],q[IW][re],dwdx,dwdy,dwdz,dx,dy,dz,ods);

                    //! It is unnecessary to consider because the turbulence viscosity coefficient is 0.
                    visLaminar = half * (visl[le] + visl[re]);
                    visTurbulence = 0.0;
                    if (wallFunctionType != WALLFUNCTION::NONE)
                    {
                        visTurbulence = half * (vist[le] + vist[re]);      
                    }
                    vis = visLaminar + visTurbulence;

                    txx = vis * two3rd * (2.0 * dudx - dvdy - dwdz);
                    tyy = vis * two3rd * (2.0 * dvdy - dwdz - dudx);
                    tzz = vis * two3rd * (2.0 * dwdz - dudx - dvdy);
                    txy = vis * (dudy + dvdx);
                    txz = vis * (dudz + dwdx);
                    tyz = vis * (dvdz + dwdy);
                    tyx = txy;
                    tzx = txz;
                    tzy = tyz;

                    fvsx = - two * (nx * txx + ny * tyx + nz * tzx) * oRefReNumber;
                    fvsy = - two * (nx * txy + ny * tyy + nz * tzy) * oRefReNumber;
                    fvsz = - two * (nx * txz + ny * tyz + nz * tzz) * oRefReNumber;

                    dfx += fvsx;
                    dfy += fvsy;
                    dfz += fvsz;
                }
                xc = xfc[iFace];
                yc = yfc[iFace];
                zc = zfc[iFace];

                dmx = (yc - TorqueRefY) * dfz - (zc - TorqueRefZ) * dfy;
                dmy = (zc - TorqueRefZ) * dfx - (xc - TorqueRefX) * dfz;
                dmz = (xc - TorqueRefX) * dfy - (yc - TorqueRefY) * dfx;

                cpx += dpx;
                cpy += dpy;
                cpz += dpz;

                cfx += dfx;
                cfy += dfy;
                cfz += dfz;

                cmx += dmx;
                cmy += dmy;
                cmz += dmz;

                RDouble dir = cosa * fvsx + sina * fvsy;

                RDouble density = q[IR][le];

                //! Compute taoWall.
                RDouble cfn  = fvsx * normalComponetX + fvsy * normalComponetY + fvsz * normalComponetZ;
                RDouble cfxt = fvsx - cfn * normalComponetX;
                RDouble cfyt = fvsy - cfn * normalComponetY;
                RDouble cfzt = fvsz - cfn * normalComponetZ;

                RDouble cft = sqrt(cfxt * cfxt + cfyt * cfyt + cfzt * cfzt);

                RDouble cf = cft / area[iFace] * SIGN(1.0,dir);
                RDouble taoWall = half * cft / area[iFace];

                //! Compute yplus.
                RDouble yPlus = 0.0;
                if (viscousType > INVISCID)
                {
                    yPlus = refReNumber * sqrt(taoWall * density) * wallDistance[le] / visLaminar;

                    //! Compute the mean free path.
                    meanFreePath = ((vis * refVicosity) / pw_Dim) * sqrt(0.5 * PI * rjmk * massReciprocal * wallTemperature / refDimensionalMolecular);
                    knudsenNumber = meanFreePath / knLength;
                }

                faceCP   [iFace] = cp;
                faceCF   [iFace] = cf;
                faceYPlus[iFace] = yPlus;
                facePW   [iFace] = pw_Dim;
                faceKN   [iFace] = knudsenNumber;
            }
        }
    }

    //! Compute Q on wall face.
    RDouble dimensionalQ = 0.0;
    if (wallTemperature >= 0)
    {
        RDouble refGama = parameters->GetRefGama();
        RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();

        RDouble velocityInflow = refMachNumber * refDimensionalSonicSpeed;
        dimensionalQ = refDimensionalDensity * velocityInflow * velocityInflow * velocityInflow * 0.001;    //! 0.001: w --> kw

        RDouble prandtlLaminar = parameters->GetPrandtlLaminar();
        RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();

        RDouble Tw  = wallTemperature / refDimensionalTemperature;
        RDouble Cv  = 1.0 / refGama / (refGama - 1.0) / refMachNumber / refMachNumber;
        RDouble tmu = pow(Tw, 1.5) * (1.0 + nonDimensionalSutherlandTemperature) / (Tw + nonDimensionalSutherlandTemperature);
        //! Twall stores the translation-rotation temperature, vibration temperature and the electron temperature.
        RDouble Twall[3] = {Tw, Tw, Tw};

        RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
            int bcType = bcRegion->GetBCType();

            if (IsWall(bcType))
            {
                bool wallTempArrayExist = false;
                RDouble *wallTempArray = nullptr;
                int countFace = 0;
                if (bcRegion->CheckFieldData("wallTempArray"))
                {
                    wallTempArrayExist = true;
                    wallTempArray = reinterpret_cast <RDouble *> (bcRegion->GetFieldDataPtr("wallTempArray"));
                }
                vector<int> *faceIndex = bcRegion->GetFaceIndex();
                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    //! iFace is the face number in the set of faceIndex.
                    int iFace = *iter;
                    int le = leftCellofFace[iFace];
                    int re = rightCellofFace[iFace];

                    if (wallTempArrayExist)
                    {
                        Tw = wallTempArray[countFace] / refDimensionalTemperature;
                        countFace ++;
                    }

                    if (viscousType > INVISCID)
                    {
                        wallVariables[IU] = 0.0;
                        wallVariables[IV] = 0.0;
                        wallVariables[IW] = 0.0;
                    }
                    else
                    {
                        //! The velocity is used to extract the wall streamline, so this velocity is not the real velocity on wall.
                        wallVariables[IU] = q[IU][le];
                        wallVariables[IV] = q[IV][le];
                        wallVariables[IW] = q[IW][le];
                    }

                    wallVariables[IR] = q[IR][le];
                    wallVariables[IP] = q[IP][le];

                    absVelocity = wallVariables[IU] * wallVariables[IU] + wallVariables[IV] * wallVariables[IV] + wallVariables[IW] * wallVariables[IW];
                    wallTotalEnthalpy = 0.0;
                    gas->ComputeEnthalpyByPrimitive(wallVariables, refGama, wallTotalEnthalpy, Twall);
                    wallTotalEnthalpy += 0.5 * absVelocity;
                    wallTotalEnthalpy = wallTotalEnthalpy * pow(refVelocity, 2);
                    deltaH = refTotalEnthalpy - wallTotalEnthalpy;

                    RDouble deltT = t[ITT][le] - Tw; 
                    dx = xcc[re] - xcc[le];
                    dy = ycc[re] - ycc[le];
                    dz = zcc[re] - zcc[le];
                    RDouble dd = dx * dx + dy * dy + dz * dz;

                    dd = sqrt(dd)/2.0;
                    deltT /= dd;

                    faceQ[iFace] = Cv * tmu * deltT * refGama * oRefReNumber / prandtlLaminar;

                    StantonNumber = (faceQ[iFace] * dimensionalQ * 1000) / (refDimensionalDensity * refVelocity * deltaH);
                    faceST[iFace] = StantonNumber;
                }
            }
        }
    }

    //! Communicate data on interface.
    int *lableCell = new int[nTotalCell];
    std::fill_n(lableCell, nTotalCell, 0);

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        if (IsWall(bcType))
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                int iFace = *iter;
                int le = leftCellofFace[iFace];
                lableCell[le] = 1;
            }
        }
        if (bcType == INTERFACE || bcType == PHENGLEI::OVERSET)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                int iFace = *iter;
                int cL = leftCellofFace[iFace];
                if (lableCell[cL] == 1)
                {
                    lableCell[cL] = 2;
                }
            }
        }
    }

    RDouble **cellQ2D = NewPointer2<RDouble>(1, nTotal);
    RDouble * cellQ   = cellQ2D[0];
    std::fill_n(cellQ, nTotal, 0);

    RDouble **cellCP2D = NewPointer2<RDouble>(1, nTotal);
    RDouble * cellCP   = cellCP2D[0];
    std::fill_n(cellCP, nTotal, 0);

    RDouble **cellCF2D = NewPointer2<RDouble>(1, nTotal);
    RDouble * cellCF   = cellCF2D[0];
    std::fill_n(cellCF, nTotal, 0);

    RDouble **cellYPlus2D = NewPointer2<RDouble>(1, nTotal);
    RDouble * cellYPlus   = cellYPlus2D[0];
    std::fill_n(cellYPlus, nTotal, 0);

    RDouble **cellPW2D = NewPointer2<RDouble>(1, nTotal);
    RDouble * cellPW   = cellPW2D[0];
    std::fill_n(cellPW, nTotal, 0);

    RDouble **cellST2D = NewPointer2<RDouble>(1, nTotal);
    RDouble * cellST   = cellST2D[0];
    std::fill_n(cellST, nTotal, 0);

    RDouble **cellKN2D = NewPointer2<RDouble>(1, nTotal);
    RDouble * cellKN   = cellKN2D[0];
    std::fill_n(cellKN, nTotal, 0);

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        if (IsWall(bcType))
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                int iFace = *iter;
                int le = leftCellofFace[iFace];
                if(lableCell[le] == 2)
                {
                    cellCF[le] = faceCF[iFace];
                    cellCP[le] = faceCP[iFace];
                    cellYPlus[le] = faceYPlus[iFace];
                    cellPW[le]    = facePW[iFace];
                    cellKN[le]    = faceKN[iFace];
                    if (wallTemperature >= 0)
                    {
                        cellQ[le] = faceQ[iFace];
                        cellST[le] = faceST[iFace];
                    }
                }
            }
        }
    }

    delete [] lableCell;    lableCell = nullptr;

    PHSPACE::CommunicateInterfaceValue(grid, cellCP2D, "cellCP", 1);
    PHSPACE::CommunicateInterfaceValue(grid, cellCF2D, "cellCF", 1);
    PHSPACE::CommunicateInterfaceValue(grid, cellYPlus2D, "cellYPlus", 1);
    PHSPACE::CommunicateInterfaceValue(grid, cellPW2D, "cellPW", 1);
    PHSPACE::CommunicateInterfaceValue(grid, cellKN2D, "cellKN", 1);
    if (wallTemperature >= 0)
    {
        PHSPACE::CommunicateInterfaceValue(grid, cellQ2D, "cellQ", 1);
        PHSPACE::CommunicateInterfaceValue(grid, cellST2D, "cellST", 1);
    }

    //! Compute data on wall node.
    RDouble *nodeWeight = new RDouble[nTotalNode]();
    std::fill_n(nodeWeight, nTotalNode, 0.0);

    RDouble *nodeQ = new RDouble[nTotalNode]();
    std::fill_n(nodeQ, nTotalNode, 0.0);

    RDouble *nodeCP = new RDouble[nTotalNode]();
    std::fill_n(nodeCP, nTotalNode, 0.0);

    RDouble *nodeCF = new RDouble[nTotalNode]();
    std::fill_n(nodeCF, nTotalNode, 0.0);

    RDouble *nodeYPlus = new RDouble[nTotalNode]();
    std::fill_n(nodeYPlus, nTotalNode, 0.0);

    RDouble *nodePW = new RDouble[nTotalNode]();
    std::fill_n(nodePW, nTotalNode, 0.0);

    RDouble *nodeST = new RDouble[nTotalNode]();
    std::fill_n(nodeST, nTotalNode, 0.0);

    RDouble *nodeKN = new RDouble[nTotalNode]();
    std::fill_n(nodeKN, nTotalNode, 0.0);

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        if (IsWall(bcType))
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                int iFace = *iter;
                int leftCellIndex = leftCellofFace[iFace];
                for (int iNode = static_cast<int>(nodePosi[iFace]); iNode < nodePosi[iFace+1]; ++ iNode)
                {
                    int nodeIndex = face2node[iNode];
                    dx = x[nodeIndex] - xcc[leftCellIndex];
                    dy = y[nodeIndex] - ycc[leftCellIndex];
                    dz = z[nodeIndex] - zcc[leftCellIndex];
                    RDouble dist = DISTANCE(dx, dy, dz);
                    RDouble weightTemp = 1.0 / dist;

                    nodeCP   [nodeIndex] += faceCP[iFace] * weightTemp;
                    nodeCF   [nodeIndex] += faceCF[iFace] * weightTemp;
                    nodeYPlus[nodeIndex] += faceYPlus[iFace] * weightTemp;
                    nodePW   [nodeIndex] += facePW[iFace] * weightTemp;
                    nodeKN   [nodeIndex] += faceKN[iFace] * weightTemp;

                    if (wallTemperature >= 0)
                    {
                        nodeQ[nodeIndex] += faceQ[iFace] * weightTemp;
                        nodeST[nodeIndex] += faceST[iFace] * weightTemp;
                    }
                    nodeWeight[nodeIndex] += weightTemp;
                }
            }
        }
    }

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (nodeWeight[iNode] > SMALL)
        {
            nodeYPlus[iNode] /= nodeWeight[iNode];
            nodePW   [iNode] /= nodeWeight[iNode];
            nodeKN   [iNode] /= nodeWeight[iNode];
        }
    }

    //! For interface.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        if (bcType == INTERFACE || bcType == PHENGLEI::OVERSET)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                int iFace = *iter;
                int cR = rightCellofFace[iFace];
                if (ABS(cellQ[cR]) > SMALL)
                {
                    for (int iNode = static_cast<int>(nodePosi[iFace]); iNode < nodePosi[iFace+1]; ++ iNode)
                    {
                        int nodeIndex = face2node[iNode];
                        if (nodeWeight[nodeIndex] > SMALL)
                        {
                            dx = x[nodeIndex] - xcc[cR];
                            dy = y[nodeIndex] - ycc[cR];
                            dz = z[nodeIndex] - zcc[cR];
                            RDouble dist = DISTANCE(dx, dy, dz);
                            RDouble weightTemp = 1.0 / dist;

                            nodeCP[nodeIndex] += cellCP[cR] * weightTemp;
                            nodeCF[nodeIndex] += cellCF[cR] * weightTemp;
                            if (wallTemperature >= 0)
                            {
                                nodeQ[nodeIndex] += cellQ[cR] * weightTemp;
                            }
                            nodeWeight[nodeIndex] += weightTemp;
                        }
                    }
                }
            }
        }
    }

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (nodeWeight[iNode] > SMALL)
        {
            nodeCP[iNode] /= nodeWeight[iNode];
            nodeCF[iNode] /= nodeWeight[iNode];
            if (wallTemperature >= 0)
            {
                nodeST[iNode] /= nodeWeight[iNode];
                nodeQ[iNode] /= nodeWeight[iNode];
            }
        }
    }

    delete [] faceCP;     faceCP = nullptr;
    delete [] faceCF;     faceCF = nullptr;
    delete [] faceQ;      faceQ = nullptr;
    delete [] nodeWeight; nodeWeight = nullptr;
    delete [] faceYPlus;  faceYPlus = nullptr;
    delete [] facePW;     facePW = nullptr;
    delete [] faceST;     faceST = nullptr;
    delete [] faceKN;     faceKN = nullptr;
    DelPointer2(cellQ2D);
    DelPointer2(cellCP2D);
    DelPointer2(cellCF2D);
    DelPointer2(cellYPlus2D);
    DelPointer2(cellPW2D);
    DelPointer2(cellST2D);
    DelPointer2(cellKN2D);

    //! Compute data on wall face center.
    RDouble *cpFaceCenter    = new RDouble[nBoundFace]();
    RDouble *cfFaceCenter    = new RDouble[nBoundFace]();
    RDouble *qFaceCenter     = new RDouble[nBoundFace]();
    RDouble *yPlusFaceCenter = new RDouble[nBoundFace]();
    RDouble *PWFaceCenter    = new RDouble[nBoundFace]();
    RDouble *STFaceCenter    = new RDouble[nBoundFace]();
    RDouble *KNFaceCenter    = new RDouble[nBoundFace]();

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        if (face2node[iFace] > 4)
        {
            RDouble valueCP    = 0.0;
            RDouble valueCF    = 0.0;
            RDouble valueQ     = 0.0;
            RDouble valueYPlus = 0.0;
            RDouble valuePW    = 0.0;
            RDouble valueST    = 0.0;
            RDouble valueKN    = 0.0;

            int nodeIndex;
            for (int iNode = static_cast<int>(nodePosi[iFace]); iNode < nodePosi[iFace+1]; ++ iNode)
            {
                nodeIndex   = face2node[iNode];
                valueCP    += nodeCP   [nodeIndex];
                valueCF    += nodeCF   [nodeIndex];
                valueYPlus += nodeYPlus[nodeIndex];
                valuePW    += nodePW   [nodeIndex];
                valueKN    += nodeKN   [nodeIndex];

                if (wallTemperature >= 0)
                {
                    valueST += nodeST[nodeIndex];
                    valueQ += nodeQ[nodeIndex];
                }
            }
            cpFaceCenter[iFace]    = valueCP    / node_number_of_each_face[iFace];
            cfFaceCenter[iFace]    = valueCF    / node_number_of_each_face[iFace];
            qFaceCenter [iFace]    = valueQ     / node_number_of_each_face[iFace];
            yPlusFaceCenter[iFace] = valueYPlus / node_number_of_each_face[iFace];
            PWFaceCenter[iFace]    = valuePW / node_number_of_each_face[iFace];
            STFaceCenter[iFace]    = valueST    / node_number_of_each_face[iFace];
            KNFaceCenter[iFace]    = valueKN    / node_number_of_each_face[iFace];
        }
    }

    set<pair<int, string> > bcnameMap;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();

        bcnameMap.insert(pair<int, string>(bcType, bcName));
    }

    //! Reconstruct grid mesh of wall.
    for (set<pair<int, string> >::iterator iter = bcnameMap.begin(); iter != bcnameMap.end(); ++ iter)
    {
        if ((*iter).first != PHENGLEI::SOLID_SURFACE) continue;

    vector < vector < int > > face2nodelist(nBoundFace);
    vector < int > linkmap;
        GetFace2NodeList(grid, *iter, linkmap, face2nodelist);

        int NumPts            = static_cast<int>(linkmap.size());
        int NumElements       = static_cast<int>(face2nodelist.size());
        int NumFaces          = NumElements * 4;
        int TotalNumFaceNodes = NumElements * 4;
        int *cell2node        = new int[TotalNumFaceNodes]();
        RDouble **dataStore   = NewPointer2<RDouble>(nWallVariables, linkmap.size());

        for (int iVar = 0; iVar < nWallVariables; ++iVar)
        {
            for (std::size_t iNode = 0; iNode < linkmap.size(); ++iNode)
            {
                dataStore[iVar][iNode] = 0.0;
            }
        }

    int count = 0;
    for (std::size_t iFace = 0; iFace < face2nodelist.size(); ++ iFace)
    {
        uint_t nodes  = face2nodelist[iFace].size();
        int nodeIndex = face2nodelist[iFace][0];
        for (int iNode = 0; iNode < nodes; ++ iNode)
        {
            cell2node[count] = face2nodelist[iFace][iNode] + 1;
            count ++;
        }

        for (uint_t ip = nodes; ip < 4; ++ ip)
        {
            cell2node[count] = nodeIndex + 1;
            count ++;
        }
    }

    RDouble *coorX     = new RDouble[NumPts]();
    RDouble *coorY     = new RDouble[NumPts]();
    RDouble *coorZ     = new RDouble[NumPts]();
    RDouble *dumpCP    = new RDouble[NumPts]();
    RDouble *dumpCF    = new RDouble[NumPts]();
    RDouble *dumpQ     = new RDouble[NumPts]();
    RDouble *dumpQNon  = new RDouble[NumPts]();
    RDouble *dumpYPlus = new RDouble[NumPts]();
    RDouble *dumpPW    = new RDouble[NumPts]();
    RDouble *dumpST    = new RDouble[NumPts]();
    RDouble *dumpKN    = new RDouble[NumPts]();

    for (std::size_t iNode = 0; iNode < linkmap.size(); ++ iNode)
    {
        int nodeIndex = linkmap[iNode];

        if (nodeIndex < nTotalNode)
        {
            coorX[iNode] = x[nodeIndex];
            coorY[iNode] = y[nodeIndex];
            coorZ[iNode] = z[nodeIndex];

            dumpCP   [iNode] = nodeCP   [nodeIndex];
            dumpCF   [iNode] = nodeCF   [nodeIndex];
            dumpYPlus[iNode] = nodeYPlus[nodeIndex];
            dumpPW   [iNode] = nodePW   [nodeIndex];
            dumpQ    [iNode] = nodeQ    [nodeIndex] * dimensionalQ;
            dumpQNon [iNode] = nodeQ    [nodeIndex];
            dumpST   [iNode] = nodeST   [nodeIndex];
            dumpKN   [iNode] = nodeKN   [nodeIndex];
        }
        else
        {
            int faceIndex = nodeIndex - nTotalNode;
            coorX[iNode] = xfc[faceIndex];
            coorY[iNode] = yfc[faceIndex];
            coorZ[iNode] = zfc[faceIndex];

            dumpCP   [iNode] = cpFaceCenter   [faceIndex];
            dumpCF   [iNode] = cfFaceCenter   [faceIndex];
            dumpYPlus[iNode] = yPlusFaceCenter[faceIndex];
            dumpPW   [iNode] = PWFaceCenter   [faceIndex];
            dumpQ    [iNode] = qFaceCenter    [faceIndex] * dimensionalQ;
            dumpQNon [iNode] = qFaceCenter    [faceIndex];
            dumpST   [iNode] = STFaceCenter   [faceIndex];
            dumpKN   [iNode] = KNFaceCenter   [faceIndex];
        }

        visualVariables[VISUAL_WALL_CP]      = dumpCP   [iNode];
        visualVariables[VISUAL_WALL_CF]      = dumpCF   [iNode];
        visualVariables[VISUAL_WALL_YPLUS]   = dumpYPlus[iNode];
        visualVariables[VISUAL_WALL_QNONDIM] = dumpQNon [iNode];
        visualVariables[VISUAL_WALL_QDIM]    = dumpQ    [iNode];
        visualVariables[VISUAL_WALL_PW]      = dumpPW   [iNode];
        visualVariables[VISUAL_WALL_ST]      = dumpST   [iNode];
        visualVariables[VISUAL_WALL_KN]      = dumpKN   [iNode];

        for (int iVar = 0; iVar < nWallVariables; ++iVar)
        {
            int varType = visualVariablesType[iVar];
            int varIndex = variablesOrder[varType];

            dataStore[iVar][iNode] = visualVariables[varIndex];
        }
    }

    int TecioMission = 1;
    PHWrite(cdata, &TecioMission, 1);

    string bcName = (*iter).second;
        cdata->WriteString(bcName);

    PHWrite(cdata, &NumPts, 1);
    PHWrite(cdata, &NumElements, 1);
    PHWrite(cdata, &NumFaces, 1);

    PHWrite(cdata, NumPts);
    PHWrite(cdata, coorX, NumPts);
    PHWrite(cdata, coorY, NumPts);
    PHWrite(cdata, coorZ, NumPts);

    int ValueLocation = 1;

        for (int iData = 0; iData < nWallVariables; ++iData)
        {
            PHWrite(cdata, NumPts);
            PHWrite(cdata, ValueLocation);
            PHWrite(cdata, dataStore[iData], NumPts);
        }

    PHWrite(cdata, &TotalNumFaceNodes, 1);
    PHWrite(cdata, cell2node, TotalNumFaceNodes);

    delete [] coorX;     coorX = nullptr;
    delete [] coorY;     coorY = nullptr;
    delete [] coorZ;     coorZ = nullptr;
    delete [] dumpCP;    dumpCP = nullptr;
    delete [] dumpCF;    dumpCF = nullptr;
    delete [] dumpQ;     dumpQ = nullptr;
    delete [] dumpQNon;  dumpQNon = nullptr;
    delete [] dumpYPlus; dumpYPlus = nullptr;
    delete [] dumpPW;    dumpPW = nullptr;
    delete [] dumpST;    dumpST = nullptr;
    delete [] dumpKN;    dumpKN = nullptr;
    delete [] cell2node; cell2node = nullptr;
    DelPointer2(dataStore);
    }

    delete [] nodeQ;           nodeQ = nullptr;
    delete [] nodeCP;          nodeCP = nullptr;
    delete [] nodeCF;          nodeCF = nullptr;
    delete [] nodeYPlus;       nodeYPlus = nullptr;
    delete [] nodePW;          nodePW = nullptr;
    delete [] nodeST;          nodeST = nullptr;
    delete [] nodeKN;          nodeKN = nullptr;
    delete [] cpFaceCenter;    cpFaceCenter = nullptr;
    delete [] cfFaceCenter;    cfFaceCenter = nullptr;
    delete [] qFaceCenter;     qFaceCenter = nullptr;
    delete [] yPlusFaceCenter; yPlusFaceCenter = nullptr;
    delete [] PWFaceCenter;    PWFaceCenter = nullptr;
    delete [] STFaceCenter;    STFaceCenter = nullptr;
    delete [] KNFaceCenter;    KNFaceCenter = nullptr;
    delete [] wallVariables;   wallVariables = nullptr;
}

//! Dump out some flow field variables, for the case of turbulent flat plate of NASA Turbulence Center.
//! Inflow condition: Ma = 0.2, Re = 5.0e6.
void NSSolverUnstruct::Turbulence_Flat_Plate_Output_NASA(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalNode = grid->GetNTotalNode();
    RDouble **qVist;

    //! 100.000 kPa = 1 bar = 14.504 psi
    //! Condition for calculation
    //! Static Pressure (psia) 6.0 = 100000 / 14.504 * 6 Pa = 41367.8986 Pa
    //! Temperature(R) 700.0 = (5/9)*700 = 388.88K

    //! Dump out the Velocity.
    vector<string> title_tecplot;
    title_tecplot.push_back("title=\"Flow Fields of PHengLEI\"");
    title_tecplot.push_back("variables=");
    title_tecplot.push_back("\"y\"");
    title_tecplot.push_back("\"u\"");
    title_tecplot.push_back("\"y+\"");
    title_tecplot.push_back("\"u+\"");    
    int nvarplot;
    nvarplot = 2;

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nEquation = GetNumberOfEquations();

    using namespace PHSPACE;
    RDouble refGama = parameters->GetRefGama();
    RDouble reference_mach_number = parameters->GetRefMachNumber();
    RDouble refReNumber = 1.0e5;
    refReNumber = parameters->GetRefReNumber();
    RDouble tsuth;
    GlobalDataBase::GetData("tsuth", &tsuth, PHDOUBLE, 1);

    fstream file;
    PHSPACE::OpenFile(file, "results/turblent_flat_plate.dat", ios_base::out);

    int wordwidth = 20;

    RDouble *vist = reinterpret_cast< RDouble * > (grid->GetDataPtr("vist"));

    int nline = 2;

    vector<RDouble> xPositionList(nline);
    xPositionList[0] = 0.97008;
    xPositionList[1] = 1.90334;

    for (std::size_t i = 0; i < title_tecplot.size(); ++ i)
    {
        file << title_tecplot[i] << "\n";
    }

    RDouble **nodevar = NewPointer2<RDouble>(nEquation + 1, nTotalNode);
    CreateNodeVariable(grid, nodevar);

    for (int m = 0; m < nline; ++ m)
    {
        PHCutPlane *cutpl = new PHCutPlane();
        CutPlane(grid, nodevar, nEquation, xPositionList[m], cutpl, X_DIR);

        cutpl->SortData(cutpl->y);

        int nT = cutpl->nTotalNode;
        file << "zone  T = \"x = " << xPositionList[m] << "\" \n";
        //
        vector<RDouble> ylist(nT);
        for (int i = 0; i < nT; ++ i)
        {
            ylist[i] = cutpl->y[i];
        }
        sort(ylist.begin(),ylist.end());

        int yindex = 1;
        for (int i = 0; i < nT; ++ i)
        {
            if (ylist[i])
            {
                yindex = i;
                break;
            }
        }

        for (int i = 0; i < nT; ++ i)
        {
            if (ylist[yindex] == cutpl->y[i])
            {
                yindex = i;
                break;
            }
        }

        for (int i = 0; i < nT; ++ i)
        {
            RDouble rm,um,vm,wm,pm,tm,xm,ym,zm,mu,muw,dudy,tauw,utau,up,yp,rw;
            RDouble uw,vw,ww,pw,yw,tw;
            xm = cutpl->x[i];
            ym = cutpl->y[i];
            zm = cutpl->z[i];

            rm = cutpl->data[0][i];
            um = cutpl->data[1][i];
            vm = cutpl->data[2][i];
            wm = cutpl->data[3][i];
            pm = cutpl->data[4][i];
            tm = refGama * reference_mach_number * reference_mach_number * pm / rm;
            mu = tm * sqrt(tm) * (1.0 + tsuth) / (tm + tsuth);

            rw = cutpl->data[0][yindex];
            uw = cutpl->data[1][yindex];
            vw = cutpl->data[2][yindex];
            ww = cutpl->data[3][yindex];
            pw = cutpl->data[4][yindex];
            tw = refGama * reference_mach_number * reference_mach_number * pm / rm;

            muw  = tm * sqrt(tm) * (1.0 + tsuth) / (tm + tsuth);
            yw   = cutpl->y[yindex];

            dudy = uw / yw;
            tauw = muw * dudy;

            utau = sqrt(tauw/(rw * refReNumber));

            up   = um / utau;
            yp   = utau * ym * refReNumber / (mu / rm);

            file << setiosflags(ios::left);
            file << setiosflags(ios::scientific);
            file << setprecision(10);
            file << setw(wordwidth) << ym 
                 << setw(wordwidth) << um
                 << setw(wordwidth) << yp
                 << setw(wordwidth) << up;
            file << "\n";
        }
        delete cutpl;
    }
    file.close();
    file.clear();
    DelPointer2(nodevar);

    //! Dump out the friction.
    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q" ));
    PHCutPlane * cutpl = new PHCutPlane();
    CutPlaneCellCenter(grid, q, nEquation, 2.0e-6, cutpl, Y_DIR);
    cutpl->SortData(cutpl->x);

    int nT = cutpl->nTotalNode;
    if(nT > 0)
    {
        ostringstream frictionFileName;
        frictionFileName << "results/turblent_flat_plate_nasa_cf_" << grid->GetZoneID() + 1<< ".dat";
        PHSPACE::OpenFile(file, frictionFileName.str().c_str(), ios_base::out);

        file << "title=\"Flow Fields of PHengLEI\"\n";
        file << "variables = \"x\" \"cf\" \n";
        file << "zone  i = " << nT  << " \n";
        for (int i = 0; i < nT; ++ i)
        {
            RDouble rm, um, vm, wm, pm, tm, xm, ym, zm, mu, dudy;
            xm = cutpl->x[i];
            ym = cutpl->y[i];
            zm = cutpl->z[i];

            rm   = cutpl->data[0][i];
            um   = cutpl->data[1][i];
            vm   = cutpl->data[2][i];
            wm   = cutpl->data[3][i];
            pm   = cutpl->data[4][i];
            tm   = refGama * reference_mach_number * reference_mach_number * pm / rm;
            mu   = tm * sqrt(tm) * (1.0 + tsuth) / (tm + tsuth);
            dudy = um / ym;

            file << setiosflags(ios::left);
            file << setiosflags(ios::scientific);
            file << setprecision(10);
            file << setw(wordwidth) << xm
                << setw(wordwidth) << 2 * mu * dudy / refReNumber
                << "\n";
        }

        frictionFileName.clear();

        PHSPACE::CloseFile(file);
    }

    delete cutpl;    cutpl = nullptr;

    //! Dump out the turbulent viscosity of x = 0.97
    qVist = new RDouble *[1];    //! Tempal 2D array for vist, to dump both x and y value.
    qVist[0] = vist;
    cutpl = new PHCutPlane();
    CutPlaneCellCenter(grid,qVist,1,0.97,cutpl,X_DIR);
    cutpl->SortData(cutpl->y);
    nT = cutpl->nTotalNode;
    if(nT > 0)
    {
        ostringstream xSliceFileName;
        xSliceFileName << "results/turblent_flat_plate_nasa_vist_x0d97_" << grid->GetZoneID() + 1 << ".dat";
        PHSPACE::OpenFile(file, xSliceFileName.str().c_str(), ios_base::out);

        file << "title=\"turbulent viscosity at x = 0.97\"\n";
        file << "variables = \"vist\"    \"Y\"\n";
        file << "zone  i = " << nT  << " \n";
        for (int i = 0; i < nT; ++ i)
        {
            RDouble ym   = cutpl->y[i];
            RDouble vist1 = cutpl->data[0][i];

            file << setiosflags(ios::left);
            file << setiosflags(ios::scientific);
            file << setprecision(10);
            file << setw(wordwidth) << vist1
                << setw(wordwidth) << ym
                << "\n";
        }

        xSliceFileName.clear();

        PHSPACE::CloseFile(file);
    }

    delete cutpl;    cutpl = nullptr;
    delete [] qVist;    qVist = nullptr;

    //! Dump out the the MAXIMUM eddy viscosity as a function of x
    using namespace PHMPI;
    if (GetNumberofGlobalZones() > 1)
    {
        return;
    }
    string maxVistOfX = "results/turblent_flat_plate_nasa_vist_maximum.dat";
    PHSPACE::OpenFile(file, maxVistOfX.c_str(), ios_base::out);

    file << "title=\"the maximum turbulent viscosity\"\n";
    file << "variables = \"X\"    \"vist\"\n";

    qVist = new RDouble *[1];    //! Tempal 2D array for vist, to dump both x and y value.
    qVist[0] = vist;

    //! Get the number of x stand position.
    int nXStand = 0;
    RDouble* xfc = grid->GetFaceCenterX();
   
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();
    
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        vector<int>* faceIndex = bcRegion->GetFaceIndex();
        if (bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::SYMMETRY)
        {
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                ++ nXStand;
            }

        }
    }

    //! Get the x stand position.
    RDouble *xOfStand = new RDouble [nXStand];
    nXStand = 0;

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        if (bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::SYMMETRY)
        {
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                int iFace = *iter;
                xOfStand[nXStand++] = xfc[iFace];
            }

        }
    }

    vector< vector < RDouble > > var_list;
    vector < RDouble > x_list;
    vector < RDouble > var(1);
    cutpl = new PHCutPlane();    //! Store the maximum viscosity of each xcut stand.

    RDouble *maxVistOfStandX = new RDouble [nXStand];
    for (int iStand = 0; iStand < nXStand; ++ iStand)
    {
        RDouble xCut = xOfStand[iStand];
        PHCutPlane *cutpX = new PHCutPlane();    //! Store the viscosity of ONE x stand, x = xcut.

        RDouble maxVist = -LARGE;
        CutPlaneCellCenter(grid,qVist,1,xCut,cutpX,X_DIR);
        nT = cutpX->nTotalNode;
        for (int i = 0; i < nT; ++ i)
        {
            RDouble vist2 = cutpX->data[0][i];
            maxVist = MAX(maxVist, vist2);
        }

        var[0] = maxVist;
        var_list.push_back(var);
        x_list.push_back(xCut);

        delete cutpX;    cutpX = nullptr;
    }

    cutpl->neqn   = 1;
    cutpl->nTotalNode = static_cast<int>(var_list.size());
    cutpl->data   = new RDouble * [ 1 ];
    for (int m = 0; m < 1; ++ m)
    {
        cutpl->data[m] = new RDouble [ cutpl->nTotalNode ];
    }
    cutpl->x = new RDouble [ cutpl->nTotalNode ];
    cutpl->y = new RDouble [ cutpl->nTotalNode ];
    cutpl->z = new RDouble [ cutpl->nTotalNode ];
    for (int iNode = 0; iNode < cutpl->nTotalNode; ++ iNode)
    {
        cutpl->x[iNode] = x_list[iNode];
        cutpl->y[iNode] = 0;
        cutpl->z[iNode] = 0;
    }
    for (int iNode = 0; iNode < cutpl->nTotalNode; ++ iNode)
    {
        for (int m = 0; m < 1; ++ m)
        {
            cutpl->data[m][iNode] = var_list[iNode][m];
        }
    }

    cutpl->SortData(cutpl->x);

    nT = cutpl->nTotalNode;
    file << "zone  i = " << nT  << " \n";
    for (int i = 0; i < nT; ++ i)
    {
        RDouble xm   = cutpl->x[i];
        RDouble vist3 = cutpl->data[0][i];

        file << setiosflags(ios::left);
        file << setiosflags(ios::scientific);
        file << setprecision(10);
        file << setw(wordwidth) << xm
            << setw(wordwidth) << vist3
            << "\n";
    }

    PHSPACE::CloseFile(file);
    delete cutpl;    cutpl = nullptr;
    delete [] xOfStand;    xOfStand = nullptr;
    delete [] maxVistOfStandX;    maxVistOfStandX = nullptr;
    delete [] qVist;    qVist = nullptr;
}

void NSSolverUnstruct::CreateNodeVariable(Grid *gridIn, RDouble **nodevar)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    int nEquation = GetNumberOfEquations();

    ComputeNodeVariable(grid, nodevar, q, nEquation);
}

void NSSolverUnstruct::ComputeAverageVisualNodeField(Grid *gridIn, RDouble **qn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **qAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("qAverage"));

    int varCount = 0;

    using namespace IDX;
    CompNodeVar(grid, qn[varCount ++], qAverage[IR]);
    CompNodeVarForVisual(grid, qn[varCount ++], qAverage[IU]);
    CompNodeVarForVisual(grid, qn[varCount ++], qAverage[IV]);
    CompNodeVarForVisual(grid, qn[varCount ++], qAverage[IW]);
    CompNodeVar(grid, qn[varCount ++], qAverage[IP]);

    RDouble *cp = ComputeCpAverage(gridIn);
    CompNodeVar(grid, qn[varCount ++], cp);
    delete [] cp;    cp = nullptr;
}

void NSSolverUnstruct::SliceVisualization(Grid *gridIn, ActionKey *actkey, RDouble **qn, int nvarplot)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    DataContainer *cdata = actkey->GetTecData();

    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();
    int *face2node = grid->GetFace2Node();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    vector < RDouble > xSlice;
    vector < RDouble > ySlice;
    vector < RDouble > zSlice;
    vector < vector <RDouble> > qnSlice;
    qnSlice.resize(nvarplot);

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int sliceAxis = parameters->GetSliceAxis();
    RDouble slicePostion = parameters->GetSlicePostion();

    RDouble *xyz = 0;
    if (sliceAxis == X_DIR)
    {
        xyz = x;
    }
    else if (sliceAxis == Y_DIR)
    {
        xyz = y;
    }
    else
    {
        xyz = z;
    }

    RDouble Small = 1.0e-10;
    int nCount = 0;

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            if (bcType != PHENGLEI::SOLID_SURFACE)
            {
                nCount += node_number_of_each_face[iFace];
                continue;
            }

            RDouble coor_min =   LARGE;
            RDouble coor_max = - LARGE;

            for (int iNode = 0; iNode < node_number_of_each_face[ iFace ]; ++ iNode)
            {
                int index = face2node[ nCount + iNode ];
                coor_min = MIN(xyz[ index ],coor_min);
                coor_max = MAX(xyz[ index ],coor_max);
            }

            if ((coor_min - slicePostion) * (coor_max - slicePostion) > 0.0)
            {
                nCount += node_number_of_each_face[ iFace ];
                continue;
            }

            for (int iNode = 0; iNode < node_number_of_each_face[ iFace ]; ++ iNode)
            {
                int Node_index0 = face2node[ nCount + iNode % node_number_of_each_face[ iFace ] ];
                int Node_index1 = face2node[ nCount + (iNode + 1) % node_number_of_each_face[ iFace ] ];

                if ((xyz[ Node_index0 ] - slicePostion) * (xyz[ Node_index1 ] - slicePostion) > 0.0) continue;

                if (fabs(xyz[ Node_index1 ] - xyz[ Node_index0 ]) < Small)
                {
                    xSlice.push_back(x[ Node_index0 ]);
                    xSlice.push_back(x[ Node_index1 ]);

                    ySlice.push_back(x[ Node_index0 ]);
                    ySlice.push_back(x[ Node_index1 ]);

                    zSlice.push_back(x[ Node_index0 ]);
                    zSlice.push_back(x[ Node_index1 ]);

                    continue;
                }

                if (sliceAxis == X_DIR)
                {
                    xSlice.push_back(slicePostion);
                }
                else if (fabs(x[ Node_index0 ] - x[ Node_index1 ]) < Small)
                {
                    xSlice.push_back(x[ Node_index0]);
                }
                else
                {
                    RDouble Node_x = (slicePostion - xyz[ Node_index0]) * (x[ Node_index1 ] - x[ Node_index0 ]) / (xyz[ Node_index1 ] - xyz[ Node_index0 ]) + x[ Node_index0 ];
                    xSlice.push_back(Node_x);
                }

                if (sliceAxis == Y_DIR)
                {
                    ySlice.push_back(slicePostion);
                }
                else if (fabs(y[ Node_index0 ] - y[ Node_index1 ]) < Small)
                {
                    ySlice.push_back(y[ Node_index0 ]);
                }
                else
                {
                    RDouble Node_y = (slicePostion - xyz[ Node_index0 ]) * (y[ Node_index1 ] - y[ Node_index0 ]) / (xyz[ Node_index1 ] - xyz[ Node_index0 ]) + y[ Node_index0 ];
                    ySlice.push_back(Node_y);
                }

                if (sliceAxis == Z_DIR)
                {
                    zSlice.push_back(slicePostion);
                }
                else if (fabs(z[ Node_index0 ] - z[ Node_index1 ]) < Small)
                {
                    zSlice.push_back(z[ Node_index0 ]);
                }
                else
                {
                    RDouble Node_z = (slicePostion - xyz[ Node_index0 ]) * (z[ Node_index1 ] - z[ Node_index0 ]) / (xyz[ Node_index1 ] - xyz[ Node_index0 ]) + z[ Node_index0 ];
                    zSlice.push_back(Node_z);
                }

                for (int iVarplot = 0; iVarplot < nvarplot; ++ iVarplot)
                {
                    if (fabs(qn[ iVarplot ][ Node_index0 ] - qn[ iVarplot ][ Node_index1 ]) < Small)
                    {
                        RDouble Node_q = qn[ iVarplot ][ Node_index0 ];
                        qnSlice[ iVarplot ].push_back(Node_q);
                    }
                    else
                    {
                        RDouble Node_q = (slicePostion - xyz[ Node_index0 ]) * (qn[ iVarplot ][ Node_index1 ] - qn[ iVarplot ][ Node_index0 ]) / (xyz[ Node_index1 ] - xyz[ Node_index0 ]) + qn[ iVarplot ][ Node_index0 ];
                        qnSlice[ iVarplot ].push_back(Node_q);
                    }
                }
            }

            nCount += node_number_of_each_face[ iFace ];
        }
    }

    size_t NumPts = xSlice.size();
    PHWrite(cdata, &NumPts, 1);

    if (NumPts == 0)  return;

    PHWrite(cdata, &xSlice[ 0 ], NumPts);
    PHWrite(cdata, &ySlice[ 0 ], NumPts);
    PHWrite(cdata, &zSlice[ 0 ], NumPts);
    for (int iVarplot = 0; iVarplot < nvarplot; ++ iVarplot)
    {
        PHWrite(cdata, &qnSlice[ iVarplot ][ 0 ], NumPts);
    }
}

void NSSolverUnstruct::BoundaryVisualization(Grid *gridIn, ActionKey *actkey, RDouble **qn, int nvarplot)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    DataContainer *cdata = actkey->GetTecData();

    int nBoundFace = grid->GetNBoundFace();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    using namespace PHENGLEI;
    vector < vector < int > > face2nodelist(nBoundFace);
    vector < int > linkmap;
    set<int> bclist;

    set< pair<int, string> > bcnameMap;

    UnstructBCSet **bcr = grid->GetBCRecord();
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcr[iFace]->GetKey();
        bclist.insert(bcType);
        string bcname = bcr[iFace]->GetBCName();

        bcnameMap.insert(pair<int,string>(bcType, bcname));
    }
    
    int TecioMission;
    set< pair<int, string> >::iterator iter;
    for (iter = bcnameMap.begin(); iter != bcnameMap.end(); ++ iter)
    {
        if ((*iter).first == PHENGLEI::INTERFACE || (*iter).first == PHENGLEI::OVERSET || (*iter).first < 0) continue;

        face2nodelist.resize(0);
        linkmap.resize(0);

        GetFace2NodeList(grid, *iter, linkmap, face2nodelist);

        if(linkmap.size() == 0)
        {
            continue;
        }

        TecioMission = WriteBoundary;
        PHWrite(cdata, &TecioMission, 1);

        int BCtype = (*iter).first;
        PHWrite(cdata, &BCtype, 1);

        size_t NumPts = linkmap.size();
        size_t NumElements = face2nodelist.size();
        size_t NumFaces = NumElements * 4;
        PHWrite(cdata, &NumPts, 1);
        PHWrite(cdata, &NumElements, 1);
        PHWrite(cdata, &NumFaces, 1);

        size_t TotalNumFaceNodes_Rect = NumElements * 4;
        int *cell2node = new int [TotalNumFaceNodes_Rect];
        int *node_number_of_each_cell = new int [NumElements];

        int count = 0;
        for (std::size_t i = 0; i < face2nodelist.size(); ++ i)
        {
            node_number_of_each_cell[i] = 4;
            int np     = static_cast<int>(face2nodelist[i].size());
            int index0 = face2nodelist[i][0];
            for (int ip = 0; ip < np; ++ ip)
            {
                cell2node[count] = face2nodelist[i][ip];
                count ++;
            }

            for (int ip = np; ip < 4; ++ ip)
            {
                cell2node[count] = index0;
                count ++;
            }
        }

        PHWrite(cdata, &TotalNumFaceNodes_Rect, 1);
        PHWrite(cdata, node_number_of_each_cell, NumElements);

        PHWrite(cdata, &NumPts, 1);

        RDouble *xx = new RDouble [NumPts];
        RDouble *yy = new RDouble [NumPts];
        RDouble *zz = new RDouble [NumPts];

        for (std::size_t i = 0; i < linkmap.size(); ++ i)
        {
            int it = linkmap[i];
            xx[i] = x[it];
            yy[i] = y[it];
            zz[i] = z[it];
        }

        PHWrite(cdata, xx, NumPts);
        PHWrite(cdata, yy, NumPts);
        PHWrite(cdata, zz, NumPts);

        for (int m = 0; m < nvarplot; ++ m)
        {
            RDouble *qq = new RDouble [NumPts];

            for (std::size_t i = 0; i < linkmap.size(); ++ i)
            {
                int it = linkmap[i];
                qq[i] = qn[m][it];
            }

            PHWrite(cdata, qq, NumPts);
            delete [] qq;    qq = nullptr;
        }

        PHWrite(cdata, cell2node, TotalNumFaceNodes_Rect);

        delete [] cell2node;    cell2node = nullptr;
        delete [] node_number_of_each_cell;    node_number_of_each_cell = nullptr;
        delete [] xx;    xx = nullptr;
        delete [] yy;    yy = nullptr;
        delete [] zz;    zz = nullptr;
    }
}

void NSSolverUnstruct::VisualizationAverageFlow(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));
    int nTotalNode = grid->GetNTotalNode();

    //! Average flow of DES statistics.
    int nvarplot = 6;

    RDouble **qn = NewPointer2<RDouble>(nvarplot, nTotalNode);
    ComputeAverageVisualNodeField(grid, qn);

    int gridID = grid->GetGridID()->GetIndex();
    PHWrite(actkey->GetTecData(), gridID);

    int dimension = grid->GetDim();
    PHWrite(actkey->GetTecData(), dimension);

    int type = grid->Type();
    PHWrite(actkey->GetTecData(), &type, 1);

    BoundaryVisualization(grid, actkey, qn, nvarplot);

    if (WantVisualField(grid))
    {
        SaveDataForTecio(grid, actkey, qn, nvarplot);
    }

    int TecioMission = Writecomplete;
    PHWrite(actkey->GetTecData(), &TecioMission, 1);

    int visualSlice = GlobalDataBase::GetIntParaFromDB("visualSlice");
    if (grid->GetDim() == THREE_D && visualSlice == 1)
    {
        SliceVisualization(grid, actkey, qn, nvarplot);
    }

    DelPointer2(qn);
}

Limiter *NSSolverUnstruct::CreateLimiter(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int limitVariables = parameters->GetLimitVariables();
    int limitVector = parameters->GetLimitVector();
    int limiterType = parameters->GetLimiterType();

    int nEquation = GetNumberOfEquations();

    int nVariables = 1;
    if (limitVector != 0 || limitVariables != 0)
    {
        nVariables = nEquation;
    }

    RDouble **q = reinterpret_cast<RDouble **>(grid->GetDataPtr("q"));

    Limiter *limiter = new Limiter(grid, q, nVariables, "limit");
    limiter->SetLimiterTypeID(limiterType);
    limiter->SetLimiterModel(limitVariables);
    limiter->SetLimiterVector(limitVector);
    limiter->InitLimitData();
    if (limiterType == ILMT_VENCAT)
    {
        int ivencat = GlobalDataBase::GetIntParaFromDB("ivencat");
        limiter->SetVencatLimiterType(ivencat);
    }

    limiter->Calculation(gridIn);

    return limiter;
}

FaceProxy * NSSolverUnstruct::CreateFaceProxy(Grid *gridIn)
{
    int nEquation = GetNumberOfEquations();

    FaceProxy *faceProxy = new FaceProxy();
    faceProxy->Create(SEGCTION_LENGTH, nEquation);
    return faceProxy;
}

GeomProxy * NSSolverUnstruct::CreateGeomProxy(Grid *gridIn)
{
    GeomProxy *geomProxy = new GeomProxy();
    geomProxy->Create(SEGCTION_LENGTH);
    return geomProxy;
}

NSFaceValue *NSSolverUnstruct::CreateNSFaceValue(Grid *gridIn)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nEquation = GetNumberOfEquations();

    int numberOfSpecies = parameters->GetNumberOfSpecies();

    NSFaceValue *facevar = new NSFaceValue(nEquation, numberOfSpecies, SEGCTION_LENGTH);
    return facevar;
}

RDouble NSSolverUnstruct::P0Wall(RDouble p1, RDouble gama, RDouble mach)
{
    RDouble p2, m2;
    if (mach > one)
    {
        p2 = p1 * ShockRelations::P2OverP1(gama, mach);
        m2 = ShockRelations::MachBehindShock(gama, mach);
        return p2  * IsentropicRelations::P0OverP(gama, m2);
    }
    else
    {
        return p1 * IsentropicRelations::P0OverP(gama, mach);
    }
}

RDouble NSSolverUnstruct::R0Wall(RDouble r1, RDouble gama, RDouble mach)
{
    RDouble r2, m2;
    if (mach > one)
    {
        r2 = r1 * ShockRelations::Rho2OverRho1(gama, mach);
        m2 = ShockRelations::MachBehindShock(gama, mach);
        return r2  * IsentropicRelations::Rho0OverRho(gama, m2);
    }
    else
    {
        return r1 * IsentropicRelations::Rho0OverRho(gama, mach);
    }
}

void NSSolverUnstruct::ComputeGradient(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    bool isViscous = parameters->IsViscous();

    //Gradient *gradientNS = 0;
    if(grid->IsFinestGrid() || isViscous)
    {
        GetGradientField(gridIn);
    }
    else
    {
        return;
    }
}

void NSSolverUnstruct::ComputeLimiter(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    Param_NSSolverUnstruct *parameters = GetControlParameters();

    RDouble **limit2D = reinterpret_cast <RDouble **> (grid->GetDataPtr("limit2D"));

    int nEquation = GetNumberOfEquations();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal     = nTotalCell + nBoundFace;

    Limiter *limiter = 0;
    int limiterType = parameters->GetLimiterType();
    int isInIniting = GlobalDataBase::GetIntParaFromDB("isInIniting");
    bool isFirstOrder = (limiterType == ILMT_FIRST) || isInIniting;
    for (int m = 0; m < nEquation; ++ m)
    {
        SetField(limit2D[m], 0.0, nTotal);
    }
    if(limiterType == ILMT_NOLIM)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            SetField(limit2D[m], 1.0, nTotal);
        }
        return;
    }

    if (grid->IsFinestGrid() && !isFirstOrder)
    {
        limiter = CreateLimiter(grid);

        int limitVector = parameters->GetLimitVector();
        bool usingVectorLimiter = (limitVector == 1);
        if (usingVectorLimiter)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                SetField(limit2D[m], limiter->GetLimiter(m), nTotal);
            }
        }
        else
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                SetField(limit2D[m], limiter->GetLimiter(0), nTotal);
            }
        }
        delete limiter;    limiter = nullptr;
    }

    //PHSPACE::CommunicateInterfaceValue(grid, limit2D, "limit2D", nEquation);
}

void NSSolverUnstruct::GetGradientField(Grid *gridIn)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    UnstructGrid *gridUnstruct = UnstructGridCast(gridIn);
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **gradPrimtiveVarX = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradPrimtiveVarY = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradPrimtiveVarZ = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarZ"));
    FieldProxy *qProxy = GetFieldProxy(gridIn, "q");
    RDouble **q     = qProxy->GetField_UNS();
    RDouble **qnode = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("qnode"));
    string gradientName = parameters->GetGradientName();

    for (int m = 0; m < nEquation; ++ m)
    {
        if (gradientName == "ggnode" || gradientName == "ggnodelaplacian")
        {
            gridUnstruct->CompGradientGGNode_NEW(q[m], qnode[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "ggcell")
        {
            gridUnstruct->CompGradientGGCell(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m], 0);
        }
        else if (gradientName == "lsq")
        {
            gridUnstruct->CompGradientLSQ(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }    
        else if (gradientName == "ggnode_weight")
        {
            gridUnstruct->CompGradientGGNodeWeight(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "gg_m2")
        {
            gridUnstruct->CompGradientGGModified2(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "ggcellnew")
        {
            gridUnstruct->CompGradientGGCellNew(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "ggcellw")
        {
            gridUnstruct->CompGradientGGCellW(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else
        {
            TK_Exit::ExceptionExit("No reconstruction method has been choosed ! /n");
        }
    }

    RDouble **gradTemperatureX = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTemperatureX"));
    RDouble **gradTemperatureY = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTemperatureY"));
    RDouble **gradTemperatureZ = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTemperatureZ"));
    RDouble **t     = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("t"));
    RDouble **tnode = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("tnode"));

    for (int m = 0; m < nTemperatureModel; ++ m)
    {
        if (gradientName == "ggnode" || gradientName == "ggnodelaplacian")
        {
            gridUnstruct->CompGradientGGNode_NEW(t[m], tnode[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else if (gradientName == "ggcell")
        {
            gridUnstruct->CompGradientGGCell(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else if (gradientName == "lsq")
        {
            gridUnstruct->CompGradientLSQ(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }    
        else if (gradientName == "ggnode_weight")
        {
            gridUnstruct->CompGradientGGNodeWeight(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else if (gradientName == "gg_m2")
        {
            gridUnstruct->CompGradientGGModified2(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else if (gradientName == "ggcellnew")
        {
            gridUnstruct->CompGradientGGCellNew(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else if (gradientName == "ggcellw")
        {
            gridUnstruct->CompGradientGGCellW(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else
        {
            TK_Exit::ExceptionExit("No reconstruction method has been choosed ! /n");
        }
    }

    delete qProxy;    qProxy = nullptr;
}

void NSSolverUnstruct::InitiallimiterOnlyForMixGrid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    RDouble **q = reinterpret_cast<RDouble **>(grid->GetDataPtr("q"));
    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int limitVector = parameters->GetLimitVector();
    int limiterType = parameters->GetLimiterType();

    int nVariables = 1;
    if (limitVector == 1)
    {
        nVariables = 5;
    }

    limiterOnlyForMixGrid = new Limiter(grid, q, nVariables, "limit");

    limiterOnlyForMixGrid->InitLimitData();

    if (limiterType == ILMT_VENCAT)
    {
        int ivencat = GlobalDataBase::GetIntParaFromDB("ivencat");
        limiterOnlyForMixGrid->SetVencatLimiterType(ivencat);
    }
}

void NSSolverUnstruct::ComputelimiterOnlyForMixGrid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    RDouble **q = reinterpret_cast<RDouble **>(grid->GetDataPtr("q"));
    Param_NSSolverUnstruct *parameters = GetControlParameters();

    RDouble **limit2D = reinterpret_cast <RDouble **> (grid->GetDataPtr("limit2D"));

    int nEquation = GetNumberOfEquations();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

    int limitVariables = parameters->GetLimitVariables();
    int limitVector = parameters->GetLimitVector();
    int limiterType = parameters->GetLimiterType();
    
    if (limiterType >= ILMT_STRUCT || limiterType == ILMT_NOLIM || limiterType == ILMT_FIRST)
    {
        return;
    }

    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    const int SCALAR_LIMITER    = 0;
    const int VECTOR_LIMITER    = 1;
    
    const int LIMIT_BY_RHO_P    = 0;
    const int LIMIT_BY_EACH_VAR = 1;

    RDouble *limit_tmp = 0;
    if (limiterType == ILMT_BARTH)
    {
        if(limitVector == SCALAR_LIMITER)
        {
            //! Each variables has the same limiter value.
            limit_tmp = limiterOnlyForMixGrid->GetLimiter(0);
            if (limitVariables == LIMIT_BY_RHO_P)
            {
                limiterOnlyForMixGrid->BarthLimiter(limit_tmp, q[IR], dqdx[IR], dqdy[IR], dqdz[IR]);
                limiterOnlyForMixGrid->BarthLimiter(limit_tmp, q[IP], dqdx[IP], dqdy[IP], dqdz[IP]);
            }
            else if (limitVariables == LIMIT_BY_EACH_VAR)
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    limiterOnlyForMixGrid->BarthLimiter(limit_tmp, q[m], dqdx[m], dqdy[m], dqdz[m]);
                }
            }
        }
        else if(limitVector == VECTOR_LIMITER)
        {
            //! Limiter value of each variable is computed by itself.
            for (int m = 0; m < nEquation; ++ m)
            {
                limit_tmp = limiterOnlyForMixGrid->GetLimiter(m);
                limiterOnlyForMixGrid->BarthLimiter(limit_tmp, q[m], dqdx[m], dqdy[m], dqdz[m]);
            }
        }
    }
    else if (limiterType == ILMT_VENCAT)
    {
        if (limitVector == SCALAR_LIMITER)
        {
            //! Each variables has the same limiter value.
            limit_tmp = limiterOnlyForMixGrid->GetLimiter(0);
            if (limitVariables == LIMIT_BY_RHO_P)
            {
                limiterOnlyForMixGrid->VencatLimiter(limit_tmp, q[IR], dqdx[IR], dqdy[IR], dqdz[IR], IR);
                limiterOnlyForMixGrid->VencatLimiter(limit_tmp, q[IP], dqdx[IP], dqdy[IP], dqdz[IP], IP);
            }
            else if (limitVariables == LIMIT_BY_EACH_VAR)
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    limiterOnlyForMixGrid->VencatLimiter(limit_tmp, q[m], dqdx[m], dqdy[m], dqdz[m], m);
                }
            }
        }
        else if(limitVector == VECTOR_LIMITER)
        {
            //! Limiter value of each variable is computed by itself.
            for (int m = 0; m < nEquation; ++ m)
            {
                limit_tmp = limiterOnlyForMixGrid->GetLimiter(m);
                limiterOnlyForMixGrid->VencatLimiter(limit_tmp, q[m], dqdx[m], dqdy[m], dqdz[m], m);
            }
        }
    }
    else
    {
        TK_Exit::ExceptionExit("Error: this unstructured limiter is not exist !\n", true);
    }
    
    for (int m = 0; m < nEquation; ++m)
    {
        SetField(limit2D[m], 0.0, nTotal);
    }
    
    if (limitVector == SCALAR_LIMITER)
    {
        limit_tmp = limiterOnlyForMixGrid->GetLimiter(0);
        SetLimitBoundary(grid, limit_tmp);
        for (int m = 0; m < nEquation; ++m)
        {
            SetField(limit2D[m], limit_tmp, nTotal);
        }
    }
    else if(limitVector == VECTOR_LIMITER)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            limit_tmp = limiterOnlyForMixGrid->GetLimiter(m);
            SetLimitBoundary(grid, limit_tmp);
            SetField(limit2D[m], limit_tmp, nTotal);
        }
    }
}

LIB_EXPORT void NSSolverUnstruct::ComputePostVisualVariables(Post_Visual * postVisualization)
{
    UnstructGrid *grid = UnstructGridCast(postVisualization->GetGrid());

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));

    if (postVisualization->GetFlowType() == AverageFlow)
    {
        RDouble **qAverage = reinterpret_cast <RDouble **> (grid->GetDataPtr("qAverage"));

        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DENSITY);
        postVisualization->UpdateVisualNodeVarPtr(varName, qAverage[IDX::IR]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_U);
        postVisualization->UpdateVisualNodeVarPtr(varName, qAverage[IDX::IU]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_V);
        postVisualization->UpdateVisualNodeVarPtr(varName, qAverage[IDX::IV]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_W);
        postVisualization->UpdateVisualNodeVarPtr(varName, qAverage[IDX::IW]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_PRESSURE);
        postVisualization->UpdateVisualNodeVarPtr(varName, qAverage[IDX::IP]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_CP);
        RDouble *cp = ComputeCpAverage(grid);
        postVisualization->UpdateVisualNodeVarPtr(varName, cp);
        delete [] cp;    cp = nullptr;

        return;
    }

    //! added by zzp 202108, for ReynoldsStress output
    if (postVisualization->GetFlowType() == AverageReynoldsStress)
    {
        RDouble **tauAverage = reinterpret_cast <RDouble **> (grid->GetDataPtr("tauAverage"));

        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_XX);
        postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage[0]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_YY);
        postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage[1]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_ZZ);
        postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage[2]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_XY);
        postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage[3]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_XZ);
        postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage[4]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_YZ);
        postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage[5]);

        return;
    }

    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_IBLANK))
    {
        int    *iBlank = grid->GetBlankIndex();
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_IBLANK);
        postVisualization->UpdateVisualNodeNotInterpolation(varName, iBlank);
    }

    using namespace IDX;
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DENSITY))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DENSITY);
        postVisualization->UpdateVisualNodeVarPtr(varName, q[IDX::IR]);
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_U))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_U);
        postVisualization->UpdateVisualNodeVarPtr(varName, q[IDX::IU]);
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_V))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_V);
        postVisualization->UpdateVisualNodeVarPtr(varName, q[IDX::IV]);
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_W))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_W);
        postVisualization->UpdateVisualNodeVarPtr(varName, q[IDX::IW]);
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_PRESSURE))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_PRESSURE);
        postVisualization->UpdateVisualNodeVarPtr(varName, q[IDX::IP]);
    }

    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_TEMPERATURE))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TEMPERATURE);
        postVisualization->UpdateVisualNodeVarPtr(varName, t[IDX::ITT]);
    }    

    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_MACH))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_MACH);
        RDouble *mach  = CompMachNumber(grid);
        postVisualization->UpdateVisualNodeVarPtr(varName, mach);
        delete [] mach;    mach = nullptr;
    }

    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_STREAMLINE_U))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_STREAMLINE_U);
        postVisualization->UpdateVisualNodeVarPtr(varName, q[IDX::IU]);
    }

    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_STREAMLINE_V))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_STREAMLINE_V);
        postVisualization->UpdateVisualNodeVarPtr(varName, q[IDX::IV]);
    }

    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_STREAMLINE_W))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_STREAMLINE_W);
        postVisualization->UpdateVisualNodeVarPtr(varName, q[IDX::IW]);
    }

    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_TIME_STEP))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TIME_STEP);
        RDouble *dt = reinterpret_cast< RDouble * > (grid->GetDataPtr("dt" ));
        postVisualization->UpdateVisualNodeVarPtr(varName, dt);
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VOLUME))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VOLUME);
        RDouble *vol = grid->GetCellVolume();
        postVisualization->UpdateVisualNodeVarPtr(varName, vol);
    }
    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_WALL_DIST))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_WALL_DIST);
        RDouble *wallDistance = grid->GetWallDist();
        postVisualization->UpdateVisualNodeVarPtr(varName, wallDistance);
    }

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();
    string viscousName = parameters->GetViscousName();
    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");

    if (viscousType == LAMINAR && iLES == NOLES_SOLVER)
    {
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VISCOSITY_LAMINAR))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VISCOSITY_LAMINAR);
            RDouble *visl = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));
            postVisualization->UpdateVisualNodeVarPtr(varName, visl);
        }
    }
    else if (viscousType > LAMINAR || iLES == LES_SOLVER)
    {
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VISCOSITY_LAMINAR))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VISCOSITY_LAMINAR);
            RDouble *visl = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));
            postVisualization->UpdateVisualNodeVarPtr(varName, visl);
        }

        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VISCOSITY_TURBULENT))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VISCOSITY_TURBULENT);
            RDouble *vist = reinterpret_cast< RDouble * > (grid->GetDataPtr("vist"));
            postVisualization->UpdateVisualNodeVarPtr(varName, vist);
        }

        if (viscousName.substr(0,6) == "2eq-kw")
        {
            RDouble **q_turb = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_turb" ));
            if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_MODELED_TKE))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_MODELED_TKE);
                postVisualization->UpdateVisualNodeVarPtr(varName, q_turb[IDX::IKE]);
            }
            if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_MODELED_DISSIPATION))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_MODELED_DISSIPATION);
                postVisualization->UpdateVisualNodeVarPtr(varName, q_turb[IDX::IKW]);
            }
            if (viscousName.substr(0,17) == "2eq-kw-menter-sst")    //! For blending function F1 and F2
            {
                RDouble **q_transition = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_transition" ));
                if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_SST_F1))
                {
                    RDouble *blend = reinterpret_cast< RDouble * > (grid->GetDataPtr("blend" ));    //! For blengding function F1
                    string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_SST_F1);
                    postVisualization->UpdateVisualNodeVarPtr(varName, blend);
                }
                if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_SST_F2))
                {
                    RDouble *SST_F2 = reinterpret_cast< RDouble * > (grid->GetDataPtr("SST_F2" ));   //! For blengding function F2
                    string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_SST_F2);
                    postVisualization->UpdateVisualNodeVarPtr(varName, SST_F2);
                }
                int transitionType = GlobalDataBase::GetIntParaFromDB("transitionType");
                if (transitionType == IREGAMA)
                {
                    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_INTERMITTENCY))
                    {
                        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_INTERMITTENCY);
                        postVisualization->UpdateVisualNodeVarPtr(varName, q_transition[IDX::IGAMA]);
                    }
                    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_MOMENTUMTHICKREYNOLDS))
                    {
                        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_MOMENTUMTHICKREYNOLDS);
                        postVisualization->UpdateVisualNodeVarPtr(varName, q_transition[IDX::IRECT]);
                    }
                }
            }
        }
    }

    int visualField = parameters->GetPlotFieldType();
    if(WantVisualField(grid) || visualField == TEC_SPACE::BlockVisual)
    {
        RDouble *vorticity_x, *vorticity_y, *vorticity_z, *vorticityMagnitude, *strain_rate, *Q_criteria;

        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal = nTotalCell + nBoundFace;

        vorticity_x = new RDouble [nTotal];
        vorticity_y = new RDouble [nTotal];
        vorticity_z = new RDouble [nTotal];
        vorticityMagnitude = new RDouble [nTotal];
        strain_rate = new RDouble [nTotal];
        Q_criteria = new RDouble [nTotal];
        ComputeVorticitybyQCriteria(grid, vorticity_x, vorticity_y, vorticity_z, vorticityMagnitude, strain_rate, Q_criteria);

        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_X))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VORTICITY_X);
            postVisualization->UpdateVisualNodeVarPtr(varName, vorticity_x);
        }

        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_Y))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VORTICITY_Y);
            postVisualization->UpdateVisualNodeVarPtr(varName, vorticity_y);
        }

        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_Z))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VORTICITY_Z);
            postVisualization->UpdateVisualNodeVarPtr(varName, vorticity_z);
        }

        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_MAGNITUDE))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VORTICITY_MAGNITUDE);
            postVisualization->UpdateVisualNodeVarPtr(varName, vorticityMagnitude);
        }

        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_STRAIN_RATE))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_STRAIN_RATE);
            postVisualization->UpdateVisualNodeVarPtr(varName, strain_rate);
        }

        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_Q_CRITERIA))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_Q_CRITERIA);
            postVisualization->UpdateVisualNodeVarPtr(varName, Q_criteria);
        }

        delete [] vorticity_x;    vorticity_x = nullptr;
        delete [] vorticity_y;    vorticity_y = nullptr;
        delete [] vorticity_z;    vorticity_z = nullptr;
        delete [] vorticityMagnitude;    vorticityMagnitude = nullptr;
        delete [] strain_rate;    strain_rate = nullptr;
        delete [] Q_criteria;    Q_criteria = nullptr;
    }

    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_CP))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_CP);
        RDouble *cp = ComputeCp(grid);
        postVisualization->UpdateVisualNodeVarPtr(varName, cp);
        delete [] cp;    cp = nullptr;
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_GRADIENT_UX))
    {
        RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));        
        RDouble *dudx = dqdx[IU];
        postVisualization->UpdateVisualNodeVarPtr("gradientUx", dudx);       
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_GRADIENT_UY))
    {
        RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));        
        RDouble *dudy = dqdy[IU];
        postVisualization->UpdateVisualNodeVarPtr("gradientUy", dudy);       
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_GRADIENT_VX))
    {
        RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));        
        RDouble *dvdx = dqdx[IV];
        postVisualization->UpdateVisualNodeVarPtr("gradientVx", dvdx);       
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_GRADIENT_VY))
    {
        RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));        
        RDouble *dvdy = dqdy[IV];
        postVisualization->UpdateVisualNodeVarPtr("gradientVy", dvdy);       
    }

    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_DENSITY)
     || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_U)
     || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_V)
     || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_W)
     || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VELOCITY_MAGNITUDE)
     || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_PRESSURE)
     || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_TEMPERATURE))
    {
        RDouble **dimensionalVariables = ComputeDimensionalVariables(grid);

        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_DENSITY))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_DENSITY);
            postVisualization->UpdateVisualNodeVarPtr(varName, dimensionalVariables[0]);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_U))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_U);
            postVisualization->UpdateVisualNodeVarPtr(varName, dimensionalVariables[1]);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_V))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_V);
            postVisualization->UpdateVisualNodeVarPtr(varName, dimensionalVariables[2]);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_W))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_W);
            postVisualization->UpdateVisualNodeVarPtr(varName, dimensionalVariables[3]);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VELOCITY_MAGNITUDE))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VELOCITY_MAGNITUDE);
            postVisualization->UpdateVisualNodeVarPtr(varName, dimensionalVariables[4]);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_PRESSURE))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_PRESSURE);
            postVisualization->UpdateVisualNodeVarPtr(varName, dimensionalVariables[5]);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_TEMPERATURE))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_TEMPERATURE);
            postVisualization->UpdateVisualNodeVarPtr(varName, dimensionalVariables[6]);
        }

        DelPointer2(dimensionalVariables);
    }

    //! Chemical species must be the last.
    int nChemical = parameters->GetChemicalFlag();
    if(nChemical)
    {
        using namespace GAS_SPACE;

        int nm = parameters->GetNSEquationNumber();
        int nEquation = GetNumberOfEquations();
        int nSpeciesNumber = parameters->GetNumberOfSpecies();

        string *varname = gas->GetNameOfSpecies();
        for (int m = nm; m < nEquation; ++ m)
        {
            string varName = "massfraction-" + varname[m-nm];
            postVisualization->UpdateVisualNodeVarPtr(varName, q[m]);
        }

        RDouble **moleFraction = ComputePrimitiveVariablesWithMoleFraction(grid);
        for (int m = 0; m < nSpeciesNumber; ++ m)
        {
            string varName = "molefraction-" + varname[m];
            postVisualization->UpdateVisualNodeVarPtr(varName, moleFraction[m]);
    }
        DelPointer2(moleFraction);
}
}

LIB_EXPORT void NSSolverUnstruct::ComputePostProbesVariables(Post_Probes * postProbesVar)
{
    UnstructGrid *grid = UnstructGridCast(postProbesVar->GetGrid());

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));
    Param_NSSolverUnstruct *parameters = GetControlParameters();

    using namespace IDX;
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_DENSITY))
    {
        string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_DENSITY);
        postProbesVar->UpdateProbesVarPtr(varName, q[IDX::IR]);
    }

    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_U))
    {
        string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_U);
        postProbesVar->UpdateProbesVarPtr(varName, q[IDX::IU]);
    }
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_V))
    {
        string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_V);
        postProbesVar->UpdateProbesVarPtr(varName, q[IDX::IV]);
    }
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_W))
    {
        string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_W);
        postProbesVar->UpdateProbesVarPtr(varName, q[IDX::IW]);
    }
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_PRESSURE))
    {
        string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_PRESSURE);
        postProbesVar->UpdateProbesVarPtr(varName, q[IDX::IP]);
    }
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_TEMPERATURE))
    {
        string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_TEMPERATURE);
        postProbesVar->UpdateProbesVarPtr(varName, t[IDX::ITT]);
    }
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_MACH))
    {
        string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_MACH);
        RDouble *mach  = CompMachNumber(grid);
        postProbesVar->UpdateProbesVarPtr(varName, mach);
        delete [] mach;    mach = nullptr;
    }
    if (postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_DIMENSIONAL_DENSITY)
        || postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_DIMENSIONAL_U)
        || postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_DIMENSIONAL_V)
        || postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_DIMENSIONAL_W)
        || postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_VELOCITY_MAGNITUDE)
        || postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_DIMENSIONAL_PRESSURE)
        || postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_DIMENSIONAL_TEMPERATURE))
    {
        RDouble **dimensionalVariables = ComputeDimensionalVariables(grid);

        if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_DIMENSIONAL_DENSITY))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_DIMENSIONAL_DENSITY);
            postProbesVar->UpdateProbesVarPtr(varName, dimensionalVariables[0]);
}
        if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_DIMENSIONAL_U))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_DIMENSIONAL_U);
            postProbesVar->UpdateProbesVarPtr(varName, dimensionalVariables[1]);
        }
        if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_DIMENSIONAL_V))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_DIMENSIONAL_V);
            postProbesVar->UpdateProbesVarPtr(varName, dimensionalVariables[2]);
        }
        if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_DIMENSIONAL_W))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_DIMENSIONAL_W);
            postProbesVar->UpdateProbesVarPtr(varName, dimensionalVariables[3]);
        }
        if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_VELOCITY_MAGNITUDE))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_VELOCITY_MAGNITUDE);
            postProbesVar->UpdateProbesVarPtr(varName, dimensionalVariables[4]);
        }
        if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_DIMENSIONAL_PRESSURE))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_DIMENSIONAL_PRESSURE);
            postProbesVar->UpdateProbesVarPtr(varName, dimensionalVariables[5]);
        }
        if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_DIMENSIONAL_TEMPERATURE))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_DIMENSIONAL_TEMPERATURE);
            postProbesVar->UpdateProbesVarPtr(varName, dimensionalVariables[6]);
        }

        DelPointer2(dimensionalVariables);
    }
    //! Chemical species must be the last.postVisualization->IsNeedVisualization
    int nChemical = parameters->GetChemicalFlag();
    if(nChemical)
    {
        using namespace GAS_SPACE;

        int nm = parameters->GetNSEquationNumber();
        int nEquation = GetNumberOfEquations();
        int nSpeciesNumber = parameters->GetNumberOfSpecies();

        string *varname = gas->GetNameOfSpecies();
        for (int m = nm; m < nEquation; ++ m)
        {
            string varName = "massfraction-" + varname[m-nm];
            postProbesVar->UpdateProbesVarPtr(varName, q[m]);
        }
        RDouble **moleFraction = ComputePrimitiveVariablesWithMoleFraction(grid);

        for (int m = 0; m < nSpeciesNumber; ++ m)
        {
            string varName = "molefraction-" + varname[m];
            postProbesVar->UpdateProbesVarPtr(varName, moleFraction[m]);
        }
    }

    int viscousType = parameters->GetViscousType();
    string viscousName = parameters->GetViscousName();
    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");

    if (viscousType == LAMINAR && iLES == NOLES_SOLVER)
    {
        if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_VISCOSITY_LAMINAR))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_VISCOSITY_LAMINAR);
            RDouble *visl = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));
            postProbesVar->UpdateProbesVarPtr(varName, visl);
        }
    }
    else if (viscousType > LAMINAR || iLES == LES_SOLVER)
    {
        if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_VISCOSITY_LAMINAR))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_VISCOSITY_LAMINAR);
            RDouble *visl = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));
            postProbesVar->UpdateProbesVarPtr(varName, visl);
        }

        if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_VISCOSITY_TURBULENT))
        {
            string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_VISCOSITY_TURBULENT);
            RDouble *vist = reinterpret_cast< RDouble * > (grid->GetDataPtr("vist"));
            postProbesVar->UpdateProbesVarPtr(varName, vist);
        }

        if (viscousName.substr(0,6) == "2eq-kw")
        {
            RDouble **q_turb = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_turb" ));
            if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_MODELED_TKE))
            {
                string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_MODELED_TKE);
                postProbesVar->UpdateProbesVarPtr(varName, q_turb[IDX::IKE]);
            }
            if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_MODELED_DISSIPATION))
            {
                string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_MODELED_DISSIPATION);
                postProbesVar->UpdateProbesVarPtr(varName, q_turb[IDX::IKW]);
            }
        }
    }
}

LIB_EXPORT void NSSolverUnstruct::InitControlParameters()
{
    FreeControlParameters();
    controlParameters = new Param_NSSolverUnstruct();
    controlParameters->Init();
}

LIB_EXPORT Param_NSSolverUnstruct * NSSolverUnstruct::GetControlParameters() const
{
    return static_cast< Param_NSSolverUnstruct * > (controlParameters);
}

void NSSolverUnstruct::CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex, const int &nEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation) return;

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
    dataContainer->MoveToBegin();

    RDouble **dq = fieldProxy->GetField_UNS();

    int nlen = GetNumberOfGhostCellLayers() * interfaceNumberBetweenTwoNeighboringZone * nEquation;
    RDouble *dqTmp = new RDouble[nlen];
    int count = 0;

    for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
    {
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int sourceCell;
        int iFace = interfaceIndexContainerForSend[ iLocalFace ];
            grid->GetSourceIndex(iFace, iGhostLayer + 1, sourceCell);
        for (int m = 0; m < nEquation; ++ m)
        {
                dqTmp[count++] = dq[m][sourceCell];
        }
    }
    }
    PHWrite(dataContainer, dqTmp, nlen);
    delete [] dqTmp;    dqTmp = nullptr;
}

void NSSolverUnstruct::DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex, const int &nEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation) return;

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);

    int *interFace2BoundaryFace = interfaceInformation->GetInterFace2BoundaryFace();

    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");
    RDouble rotationAngle = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("rotationAngle");
    rotationAngle = rotationAngle * PI / 180.0;

    dataContainer->MoveToBegin();
    RDouble **dq = fieldProxy->GetField_UNS();


    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBoundFace = grid->GetNBoundFace();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");
    if (referenceFrame == ROTATIONAL_FRAME)
    {
        int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
        RDouble PeriodicRotationAngle[100];
        GlobalDataBase::GetData("PeriodicRotationAngle", &PeriodicRotationAngle, PHDOUBLE, nTurboZone);

        //! Parallel
        int iTurboZone = grid->GetOrdinaryGridIndex();
        //! Serial
        if (iTurboZone == -1)
        {
            iTurboZone = grid->GetZoneID();
        }

        PeriodicRotationAngle[iTurboZone] = PeriodicRotationAngle[iTurboZone] * PI / 180.0;
        rotationAngle = PeriodicRotationAngle[iTurboZone];
    }

    for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
    {
        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int iFace = interfaceIndexContainerForReceive[ iLocalFace ];
            int targetCell;
            grid->GetTargetIndex(iFace, iGhostLayer + 1, targetCell);

            for (int m = 0; m < nEquation; ++ m)
            {
                PHRead(dataContainer, dq[m][targetCell]);
            }

            int iBFace = interFace2BoundaryFace[iFace];
            int le = leftCellofFace[iBFace];
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iBFace]);
            string bcName = bcRegion->GetBCName();

            if (referenceFrame == ROTATIONAL_FRAME)
            {
                int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
                string Periodic_Name[100];
                GlobalDataBase::GetData("Periodic_Name", &Periodic_Name, PHSTRING, 2 * nTurboZone);
                RDouble PeriodicRotationAngle[100];
                GlobalDataBase::GetData("PeriodicRotationAngle", &PeriodicRotationAngle, PHDOUBLE, nTurboZone);

                //! Parallel
                int iTurboZone = grid->GetOrdinaryGridIndex();
                //! Serial
                if (iTurboZone == -1)
                {
                    iTurboZone = grid->GetZoneID();
                }

                if (bcName == Periodic_Name[2 * iTurboZone] || bcName == Periodic_Name[2 * iTurboZone + 1])
                {
                    RDouble rotfg[3] = { 0, 0, 0 };
                    if (bcName == Periodic_Name[2 * iTurboZone])
                    {
                        rotfg[0] = dq[1][targetCell];
                        rotfg[1] = dq[2][targetCell] * cos(2.0 * PI - rotationAngle) - dq[3][targetCell] * sin(2.0 * PI - rotationAngle);
                        rotfg[2] = dq[2][targetCell] * sin(2.0 * PI - rotationAngle) + dq[3][targetCell] * cos(2.0 * PI - rotationAngle);
                    }
                    else if (bcName == Periodic_Name[2 * iTurboZone + 1])
                    {
                        rotfg[0] = dq[1][targetCell];
                        rotfg[1] = dq[2][targetCell] * cos(rotationAngle) - dq[3][targetCell] * sin(rotationAngle);
                        rotfg[2] = dq[2][targetCell] * sin(rotationAngle) + dq[3][targetCell] * cos(rotationAngle);
                    }

                    for (int m = 1; m <= 3; ++m)
                    {
                        dq[m][targetCell] = rotfg[m - 1];
                    }
                }
            }
            else
            {
                if (periodicType == ROTATIONAL_PERIODICITY)
                {
                    if (bcName == "Periodic_up" || bcName == "Periodic_down")
                    {
                        RDouble rotfg[3] = { 0, 0, 0 };
                        if (bcName == "Periodic_up")
                        {
                            rotfg[0] = dq[1][targetCell];
                            rotfg[1] = dq[2][targetCell] * cos(2.0 * PI - rotationAngle) - dq[3][targetCell] * sin(2.0 * PI - rotationAngle);
                            rotfg[2] = dq[2][targetCell] * sin(2.0 * PI - rotationAngle) + dq[3][targetCell] * cos(2.0 * PI - rotationAngle);
                        }
                        else if (bcName == "Periodic_down")
                        {
                            rotfg[0] = dq[1][targetCell];
                            rotfg[1] = dq[2][targetCell] * cos(rotationAngle) - dq[3][targetCell] * sin(rotationAngle);
                            rotfg[2] = dq[2][targetCell] * sin(rotationAngle) + dq[3][targetCell] * cos(rotationAngle);
                        }

                        for (int m = 1; m <= 3; ++m)
                        {
                            dq[m][targetCell] = rotfg[m - 1];
                        }
                    }
                }
            }
        }
    }
}

void NSSolverUnstruct::CompressSpecificArrayToInterface(DataContainer *&dataContainer, const string &fieldName, Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
    dataContainer->MoveToBegin();
    if (nEquation == 0)
    {
        RDouble *fieldSend = reinterpret_cast <RDouble *> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int sourceCell;
                int iFace = interfaceIndexContainerForSend[iLocalFace];
                grid->GetSourceIndex(iFace, iGhostLayer + 1, sourceCell);
                PHWrite(dataContainer, fieldSend[sourceCell]);
            }
        }
    }
    else
    {
        RDouble **fieldSend = reinterpret_cast <RDouble **> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int sourceCell;
                int iFace = interfaceIndexContainerForSend[iLocalFace];
                grid->GetSourceIndex(iFace, iGhostLayer + 1, sourceCell);
                for (int m = 0; m < nEquation; ++ m)
                {
                    PHWrite(dataContainer, fieldSend[m][sourceCell]);
                }
            }
        }
    }
}

void NSSolverUnstruct::DecompressArrayFromInterface(DataContainer *&dataContainer, const string &fieldName, Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);
    dataContainer->MoveToBegin();
    if (nEquation == 0)
    {
        RDouble *fieldRecv = reinterpret_cast <RDouble *> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int targetCell;
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                grid->GetTargetIndex(iFace, iGhostLayer + 1, targetCell);
                PHRead(dataContainer, fieldRecv[targetCell]);
            }
        }
    }
    else
    {
        RDouble **fieldRecv = reinterpret_cast <RDouble **> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                int targetCell;
                grid->GetTargetIndex(iFace, iGhostLayer + 1, targetCell);
                for (int m = 0; m < nEquation; ++ m)
                {
                    PHRead(dataContainer, fieldRecv[m][targetCell]);
                }
            }
        }
    }
}

void NSSolverUnstruct::OutflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    //! Get the number of equations.
    int nTotalCell = grid->GetNTotalCell();

    int nEquation = GetNumberOfEquations();

    int iFace, le, re;
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

#ifdef USE_GMRESSOLVER
    //! GMRESBoundary GMRES3D
    RDouble **dDdP = reinterpret_cast <RDouble**> (grid->GetDataPtr("dDdP"));
#endif

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        //! iFace is the face number in the set of faceIndex.
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];
        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = q[m][le];
        }
#ifdef USE_GMRESSOLVER
        //! GMRESBoundary
        int idx = (re - nTotalCell) * nEquation;
        for (int m = 0; m < nEquation; ++m)
        {
            dDdP[m][idx + m] = 1.0;
        }
#endif
    }
}

void NSSolverUnstruct::ReSetOversetBoundary(Grid* gridIn)
{
    int isOversetSlip = PHSPACE::GlobalDataBase::GetIntParaFromDB("isOversetSlip");
    if (isOversetSlip)
    {
        return;
    }
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    vector<int> * c2c = grid->GetCell2Cell();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));
    int *iBlank = grid->GetBlankIndex();

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nTemperatureModel = parameters->GetTemperatureModel();
    int nEquation = GetNumberOfEquations();

    for (int iCell = 0; iCell < nTotalCell; iCell++)
    {
        if (iBlank[iCell] == INTERPOLATION)
        {
            int nNeighbor = static_cast<int>(c2c[iCell].size());
            for (int j = 0; j < nNeighbor; ++j)
            {
                int neighborCell = c2c[iCell][j];
                if (neighborCell < nTotalCell && iBlank[neighborCell] == INACTIVE)
                {
                    for (int m = 0; m < nEquation; ++m)
                    {
                        q[m][neighborCell] = q[m][iCell];
                    }
                    for (int m = 0; m < nTemperatureModel; ++m)
                    {
                        t[m][neighborCell] = t[m][iCell];
                    }
                }
            }
        }
    }
}

void NSSolverUnstruct::SymmetryBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme"); 
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *xfn         = grid->GetFaceNormalX();
    RDouble *yfn         = grid->GetFaceNormalY();
    RDouble *zfn         = grid->GetFaceNormalZ();
    RDouble *vgn         = grid->GetFaceNormalVelocity();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **dDdP = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP"));

    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    //! Get the number of equations.
    int nEquation = GetNumberOfEquations();
    int nTotalCell = grid->GetNTotalCell();

    int iFace = 0, le = 0, re = 0;
    if (tscheme != GMRES)
    {
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            iFace = *iter;
            le = leftCellofFace[iFace];
            re = rightCellofFace[iFace];

            for (int m = 0; m < nEquation; ++ m)
            {
                q[m][re] = q[m][le];
            }

            RDouble uL = q[IU][le];
            RDouble vL = q[IV][le];
            RDouble wL = q[IW][le];
            RDouble vn = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL - vgn[iFace];

            q[IU][re] = uL - two * xfn[iFace] * vn;
            q[IV][re] = vL - two * yfn[iFace] * vn;
            q[IW][re] = wL - two * zfn[iFace] * vn;
        }
    }
    else
    {
#ifdef USE_GMRESSOLVER
        //! Sacado
        ADReal *ql = new ADReal[nEquation]();
        ADReal *qr = new ADReal[nEquation]();

        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            iFace = *iter;
            le = leftCellofFace[iFace];
            re = rightCellofFace[iFace];

            //! init and define the sequence order of independent variables
            for (int m = 0; m < nEquation; m++)
            {
                ql[m] = q[m][le];
                ql[m].diff(m, nEquation);
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                qr[m] = ql[m];
            }

            ADReal uL = ql[IU];
            ADReal vL = ql[IV];
            ADReal wL = ql[IW];
            ADReal vn = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL - vgn[iFace];

            qr[IU] = uL - two * xfn[iFace] * vn;
            qr[IV] = vL - two * yfn[iFace] * vn;
            qr[IW] = wL - two * zfn[iFace] * vn;

            int idx = (re - nTotalCell) * nEquation;
            for (int m = 0; m < nEquation; m++)
            {
                q[m][re] = qr[m].val();
                for (int n = 0; n < nEquation; n++)
                {
                    dDdP[m][idx + n] = qr[m].dx(n);
                }
            }
        }

        delete [] ql;    ql = nullptr;
        delete [] qr;    qr = nullptr;
#endif
    }
}

#ifdef USE_GMRESSOLVER
//! GMRESBoundary AD   GMRESAD
void NSSolverUnstruct::GMRES_SymmetryBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble   *xfn       = grid->GetFaceNormalX();
    RDouble   *yfn       = grid->GetFaceNormalY();
    RDouble   *zfn       = grid->GetFaceNormalZ();
    RDouble   *vgn       = grid->GetFaceNormalVelocity();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **dDdP = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();
    int nTotalCell = grid->GetNTotalCell();

    int iFace, le, re;

    //! Sacado
    ADReal *ql = new ADReal[nEquation]();
    ADReal *qr = new ADReal[nEquation]();

    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        //! init and define the sequence order of independent variables
        for (int m = 0; m < nEquation; m++)
        {
            ql[m] = q[m][le];
            ql[m].diff(m, nEquation);
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            qr[m] = ql[m];
        }

        ADReal uL = ql[IU];
        ADReal vL = ql[IV];
        ADReal wL = ql[IW];
        ADReal vn = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL - vgn[iFace];

        qr[IU] = uL - two * xfn[iFace] * vn;
        qr[IV] = vL - two * yfn[iFace] * vn;
        qr[IW] = wL - two * zfn[iFace] * vn;

        int idx = (re - nTotalCell) * nEquation;
        for (int m = 0; m < nEquation; m++)
        {
            q[m][re] = qr[m].val();
            for (int n = 0; n < nEquation; n++)
            {
                dDdP[m][idx + n] = qr[m].dx(n);
            }
        }
    }

    delete [] ql; ql = nullptr;
    delete [] qr; qr = nullptr;
}

//! GMRESPV GMRESBoundary
void NSSolverUnstruct::GMRES_SymmetryBCRegion_FD(Grid * gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **dDdP = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();
    int nTotalCells = grid->GetNTotalCell();

    RDouble dqpert;
    RDouble *primL = new RDouble[nEquation]();
    RDouble *primR = new RDouble[nEquation]();
    RDouble temperaturesL[3] = { 0 }, temperaturesR[3] = { 0 };
    RDouble perturbScale = 0.000001;

    int iFace, le, re;

    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        //! firstly, update the boundary, q[m][le] and q[m][re] are the standard
        for (int m = 0; m < nEquation; m++)
        {
            primL[m] = q[m][le];
        } 
        Cal_GMRES_SymmetryBCRegion(gridIn, primL, primR, iFace);
        for (int m = 0; m < nEquation; m++)
        {
            q[m][re] = primR[m];
        }

        //! then perturb the primative variables at the left cell
        //! perturb rho
        int idx = (re - nTotalCells) * nEquation;

        for (int m = 0; m < nEquation; m++)
        {
            primL[m] = q[m][le];
        }
        primL[IR] *= (1.0+perturbScale);
        dqpert     = primL[IR] - q[IR][le];
        Cal_GMRES_SymmetryBCRegion(gridIn,primL,primR,iFace);
        for (int m = 0; m < nEquation; m++)
        {
            dDdP[m][idx] = (primR[m] - q[m][re]) / dqpert;
        }

        //! perturb u
        idx++;
        for (int m = 0; m < nEquation; m++)
        {
            primL[m] = q[m][le];
        }
        primL[IU] += (perturbScale);
        dqpert     = primL[IU] - q[IU][le];
        Cal_GMRES_SymmetryBCRegion(gridIn,primL,primR,iFace);
        for (int m = 0; m < nEquation; m++)
        {
            dDdP[m][idx] = (primR[m] - q[m][re]) / dqpert;
        }

        //! perturb v
        idx++;
        for (int m = 0; m < nEquation; m++)
        {
            primL[m] = q[m][le];
        }
        primL[IV] += (perturbScale);
        dqpert     = primL[IV] - q[IV][le];
        Cal_GMRES_SymmetryBCRegion(gridIn,primL,primR,iFace);
        for (int m = 0; m < nEquation; m++)
        {
            dDdP[m][idx] = (primR[m] - q[m][re]) / dqpert;
        }

        //! perturb w
        idx++;
        for (int m = 0; m < nEquation; m++)
        {
            primL[m] = q[m][le];
        }
        primL[IW] += (1.0+perturbScale);
        dqpert     = primL[IW] - q[IW][le];
        Cal_GMRES_SymmetryBCRegion(gridIn,primL,primR,iFace);
        for (int m = 0; m < nEquation; m++)
        {
            dDdP[m][idx] = (primR[m] - q[m][re]) / dqpert;
        }

        //! perturb p
        idx++;
        for (int m = 0; m < nEquation; m++)
        {
            primL[m] = q[m][le];
        }
        primL[IP] *= (1.0+perturbScale);
        dqpert     = primL[IP] - q[IP][le];
        Cal_GMRES_SymmetryBCRegion(gridIn,primL,primR,iFace);
        for (int m = 0; m < nEquation; m++)
        {
            dDdP[m][idx] = (primR[m] - q[m][re]) / dqpert;
        }
    }
    delete [] primL;    primL = nullptr;
    delete [] primR;    primR    = nullptr;
}

//! GMRESBoundary
void NSSolverUnstruct::GMRES_SymmetryBCRegion_CV(Grid * gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *xfn         = grid->GetFaceNormalX();
    RDouble *yfn         = grid->GetFaceNormalY();
    RDouble *zfn         = grid->GetFaceNormalZ();
    RDouble *vgn         = grid->GetFaceNormalVelocity();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **dDdP = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP"));
    RDouble *gama = reinterpret_cast<RDouble *>(grid->GetDataPtr("gama"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();
    int nTotalCells = grid->GetNTotalCell();

    RDouble dCvpert;
    RDouble *primL    = new RDouble[nEquation]();
    RDouble *primR    = new RDouble[nEquation]();
    RDouble *qCvL     = new RDouble[nEquation]();
    RDouble *qCvR     = new RDouble[nEquation]();
    RDouble *qCvLpert = new RDouble[nEquation]();
    RDouble *qCvRpert = new RDouble[nEquation]();
    RDouble temperaturesL[3] = { 0 }, temperaturesR[3] = { 0 };
    RDouble perturbScale = 0.001;

    int iFace, le, re;

    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        RDouble gmL = gama[le];
        RDouble gmR = gama[re];

        //! firstly, update the boundary
        for (int m = 0; m < nEquation; m++)
        {
            primL[m] = q[m][le];
        }

        gas->Primitive2Conservative(primL, gmL, temperaturesL[ITV], temperaturesL[ITE], qCvL);
        // gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);

        Cal_GMRES_SymmetryBCRegion(gridIn, primL, primR, iFace);
        for (int m = 0; m < nEquation; m++)
        {
            q[m][re] = primR[m];
        }
        gas->Primitive2Conservative(primR, gmR, temperaturesR[ITV], temperaturesR[ITE], qCvR);    //! Set standard

        //! then perturb the conservative variables at the left cell
        //! perturb rho
        int idx = (re - nTotalCells) * nEquation;

        for (int m = 0; m < nEquation; m++)
        {
            qCvLpert[m] = qCvL[m];
        }
        qCvLpert[IR] *= (1.0+perturbScale);
        dCvpert       = qCvLpert[IR] - qCvL[IR];
        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);    //! now [primL] becomes a new value
        Cal_GMRES_SymmetryBCRegion(gridIn,primL,primR,iFace);
        gas->Primitive2Conservative(primR, gmR, temperaturesR[ITV], temperaturesR[ITE], qCvRpert);
        for (int m = 0; m < nEquation; m++)
        {
            dDdP[m][idx] = (qCvRpert[m] - qCvR[m]) / dCvpert;
        }

        //! perturb rhou
        idx++;
        for (int m = 0; m < nEquation; m++)
        {
            qCvLpert[m] = qCvL[m];
        }
        qCvLpert[IU] += (perturbScale);
        dCvpert       = qCvLpert[IU] - qCvL[IU];
        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);    //! now [primL] becomes a new value
        Cal_GMRES_SymmetryBCRegion(gridIn,primL,primR,iFace);
        gas->Primitive2Conservative(primR, gmR, temperaturesR[ITV], temperaturesR[ITE], qCvRpert);
        for (int m = 0; m < nEquation; m++)
        {
            dDdP[m][idx] = (qCvRpert[m] - qCvR[m]) / dCvpert;
        }

        //! perturb rhov
        idx++;
        for (int m = 0; m < nEquation; m++)
        {
            qCvLpert[m] = qCvL[m];
        }
        qCvLpert[IV] += (perturbScale);
        dCvpert       = qCvLpert[IV] - qCvL[IV];
        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);    //! now [primL] becomes a new value
        Cal_GMRES_SymmetryBCRegion(gridIn,primL,primR,iFace);
        gas->Primitive2Conservative(primR, gmR, temperaturesR[ITV], temperaturesR[ITE], qCvRpert);
        for (int m = 0; m < nEquation; m++)
        {
            dDdP[m][idx] = (qCvRpert[m] - qCvR[m]) / dCvpert;
        }

        //! perturb rhow
        idx++;
        for (int m = 0; m < nEquation; m++)
        {
            qCvLpert[m] = qCvL[m];
        }
        qCvLpert[IW] += (perturbScale);
        dCvpert       = qCvLpert[IW] - qCvL[IW];
        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);    //! now [primL] becomes a new value
        Cal_GMRES_SymmetryBCRegion(gridIn,primL,primR,iFace);
        gas->Primitive2Conservative(primR, gmR, temperaturesR[ITV], temperaturesR[ITE], qCvRpert);
        for (int m = 0; m < nEquation; m++)
        {
            dDdP[m][idx] = (qCvRpert[m] - qCvR[m]) / dCvpert;
        }

        //! perturb rhoe
        idx++;
        for (int m = 0; m < nEquation; m++)
        {
            qCvLpert[m] = qCvL[m];
        }
        qCvLpert[IP] *= (1.0+perturbScale);
        dCvpert       = qCvLpert[IP] - qCvL[IP];
        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);    //! now [primL] becomes a new value
        Cal_GMRES_SymmetryBCRegion(gridIn,primL,primR,iFace);
        gas->Primitive2Conservative(primR, gmR, temperaturesR[ITV], temperaturesR[ITE], qCvRpert);
        for (int m = 0; m < nEquation; m++)
        {
            dDdP[m][idx] = (qCvRpert[m] - qCvR[m]) / dCvpert;
        }
    }

    delete [] primL;       primL    = nullptr;
    delete [] primR;       primR    = nullptr;
    delete [] qCvL;        qCvL     = nullptr;
    delete [] qCvR;        qCvR     = nullptr;
    delete [] qCvLpert;    qCvLpert = nullptr;
    delete [] qCvRpert;    qCvRpert = nullptr;
}

void NSSolverUnstruct::Cal_GMRES_SymmetryBCRegion(Grid *gridIn, RDouble *PrimL, RDouble *PrimR, const int &iFace)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid = UnstructGridCast(gridIn);
    RDouble *xfn       = grid->GetFaceNormalX();
    RDouble *yfn       = grid->GetFaceNormalY();
    RDouble *zfn       = grid->GetFaceNormalZ();
    RDouble *vgn       = grid->GetFaceNormalVelocity();

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    for (int m = 0; m < nEquation; ++m)
    {
        PrimR[m] = PrimL[m];
    }

    RDouble uL = PrimL[IU];
    RDouble vL = PrimL[IV];
    RDouble wL = PrimL[IW];
    RDouble vn = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL - vgn[iFace];

    PrimR[IU] = uL - two * xfn[iFace] * vn;
    PrimR[IV] = vL - two * yfn[iFace] * vn;
    PrimR[IW] = wL - two * zfn[iFace] * vn;
}

void NSSolverUnstruct::GMRES_ViscousIsotropicWallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *faceVelocityX = grid->GetFaceVelocityX();
    RDouble *faceVelocityY = grid->GetFaceVelocityY();
    RDouble *faceVelocityZ = grid->GetFaceVelocityZ();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **dDdP = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();
    int nTotalCells = grid->GetNTotalCell();

    int nEquation = GetNumberOfEquations();

    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();

    //! Get the variables defined in boundarycondition.hypara file.
    Data_Param *bcParamDB = bcRegionUnstruct->GetBCParamDataBase();
    RDouble wallTemperature = 0.0;

    if (bcParamDB)
    {
        if (!bcParamDB->IsExist("wallTemperature", PHDOUBLE, 1))
        {
            wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
            bcParamDB->UpdateData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
        }
        else
        {
            bcParamDB->GetData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
        }
    }
    else
    {
        wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
    }
    wallTemperature /= refDimensionalTemperature;

    bool wallTempArrayExist = false;
    RDouble *wallTempArray;
    int countFace = 0;
    if (bcRegionUnstruct->CheckFieldData("wallTempArray"))
    {
        wallTempArrayExist = true;
        wallTempArray = reinterpret_cast <RDouble *> (bcRegionUnstruct->GetFieldDataPtr("wallTempArray"));
    }

    int iFace, le, re;
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    RDouble molecularWeight = one;
    using namespace GAS_SPACE;
    RDouble coefficientofstateEquation = gas->GetCoefficientOfStateEquation();

    ADReal rGhost = 0.0, pGhost = 0.0, tGhost = 0.0;    //! rho/pressure/temprature on ghost cell.
    RDouble *qL = nullptr;
    if (nChemical == 1)
    {
        // qL = new RDouble [nEquation];
    }

    ADReal * qleft  = new ADReal [nEquation]();
    ADReal * qright = new ADReal [nEquation]();

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            qleft[m] = q[m][le];
            qleft[m].diff(m,nEquation);
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            qright[m] = qleft[m];
        }

        qright[IU] = - qleft[IU] + two * faceVelocityX[iFace];
        qright[IV] = - qleft[IV] + two * faceVelocityY[iFace];
        qright[IW] = - qleft[IW] + two * faceVelocityZ[iFace];

        if (nChemical == 1)
        {
            // for (int m = 0; m < nEquation; ++ m)
            // {
            //     qL[m] = q[m][le];
            // }
            // molecularWeight = gas->ComputeMolecularWeightReciprocal(qL);
        }

        ADReal tL = (1.0 / coefficientofstateEquation) * qleft[IP] / qleft[IR];

        //! Pressure keep the same on the wall normal direction.
        pGhost = qleft[IP];

        if (wallTempArrayExist)
        {
            wallTemperature = wallTempArray[countFace] / refDimensionalTemperature;
            countFace ++;
        }

        //! Extrapolate temperature from inside to outside.
        tGhost = two * wallTemperature - tL;
        if (tGhost <= 0.0) tGhost = wallTemperature;

        if (nChemical == 1)
        {
            // tGhost = wallTemperature;
        }

        //! Then compute the density from the extrapolated temperature.
        rGhost = pGhost / (coefficientofstateEquation * tGhost * molecularWeight);

        qright[IR] = rGhost;
        qright[IP] = pGhost;

        int idx = (re - nTotalCells)*nEquation;

        for(int m = 0; m < nEquation; m++)
        {
            q[m][re] = qright[m].val();
            for(int n = 0; n < nEquation; n++)
            {
                dDdP[m][idx+n] = qright[m].dx(n);
            }
        }
    }

    if (nChemical == 1)
    {
        // delete [] qL;
    }
    delete [] qleft;    qleft = nullptr;
    delete [] qright;    qright = nullptr;
}

//! GMRESVis GMRESBoundary
void NSSolverUnstruct::GMRES_ViscousAdiabaticWallBCRegion_FD(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    int nTotalCells = grid->GetNTotalCell();
    RDouble *faceVelocityX = grid->GetFaceVelocityX();
    RDouble *faceVelocityY = grid->GetFaceVelocityY();
    RDouble *faceVelocityZ = grid->GetFaceVelocityZ();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));
    RDouble **dDdP = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();

    RDouble uWall = 0.0;
    RDouble vWall = 0.0;
    RDouble wWall = 0.0;

    if(bcData)
    {
        if (bcData->IsExist("uWall", PHDOUBLE, 1))
        {
            bcData->GetData("uWall", &uWall, PHDOUBLE, 1);
            bcData->GetData("vWall", &vWall, PHDOUBLE, 1);
            bcData->GetData("wWall", &wWall, PHDOUBLE, 1);
        }
    }    

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    int iFace, le, re;
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    RDouble molecularWeight = one;
    using namespace GAS_SPACE;
    RDouble coefficientofstateEquation = gas->GetCoefficientOfStateEquation();
    RDouble rGhost = 0.0, pGhost = 0.0, tGhost = 0.0;
    RDouble *qL = nullptr;
    if (nChemical == 1)
    {
        qL = new RDouble [nEquation];
    }
    RDouble *primL = new RDouble [nEquation];
    RDouble *primR = new RDouble [nEquation];
    RDouble dqpert;
    RDouble perturbScale = 1e-6;

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            primL[m] = q[m][le];
        }

        Cal_GMRES_ViscousAdiabaticWallBCRegion(gridIn, primL, primR, iFace, uWall, vWall, wWall);

        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = primR[m];
        }

        //! then perturb the primative variables at the left cell
        //! perturb rho
        int idx = (re - nTotalCells) * nEquation;

        for (int m = 0; m < nEquation; m++)
        {
            primL[m] = q[m][le];
        }
        primL[IR] *= (1.0+perturbScale);
        dqpert     = primL[IR] - q[IR][le];
        Cal_GMRES_ViscousAdiabaticWallBCRegion(gridIn, primL, primR, iFace, uWall, vWall, wWall);
        for (int m = 0; m < nEquation; m++)
        {
            dDdP[m][idx] = (primR[m] - q[m][re]) / dqpert;
        }

        //! perturb u
        idx++;
        for (int m = 0; m < nEquation; m++)
        {
            primL[m] = q[m][le];
        }
        primL[IU] += (perturbScale);
        dqpert     = primL[IU] - q[IU][le];
        Cal_GMRES_ViscousAdiabaticWallBCRegion(gridIn, primL, primR, iFace, uWall, vWall, wWall);
        for (int m = 0; m < nEquation; m++)
        {
            dDdP[m][idx] = (primR[m] - q[m][re]) / dqpert;
        }

        //! perturb v
        idx++;
        for (int m = 0; m < nEquation; m++)
        {
            primL[m] = q[m][le];
        }
        primL[IV] += (perturbScale);
        dqpert     = primL[IV] - q[IV][le];
        Cal_GMRES_ViscousAdiabaticWallBCRegion(gridIn, primL, primR, iFace, uWall, vWall, wWall);
        for (int m = 0; m < nEquation; m++)
        {
            dDdP[m][idx] = (primR[m] - q[m][re]) / dqpert;
        }

        //! perturb w
        idx++;
        for (int m = 0; m < nEquation; m++)
        {
            primL[m] = q[m][le];
        }
        primL[IW] += (1.0+perturbScale);
        dqpert     = primL[IW] - q[IW][le];
        Cal_GMRES_ViscousAdiabaticWallBCRegion(gridIn, primL, primR, iFace, uWall, vWall, wWall);
        for (int m = 0; m < nEquation; m++)
        {
            dDdP[m][idx] = (primR[m] - q[m][re]) / dqpert;
        }

        //! perturb p
        idx++;
        for (int m = 0; m < nEquation; m++)
        {
            primL[m] = q[m][le];
        }
        primL[IP] *= (1.0+perturbScale);
        dqpert     = primL[IP] - q[IP][le];
        Cal_GMRES_ViscousAdiabaticWallBCRegion(gridIn, primL, primR, iFace, uWall, vWall, wWall);
        for (int m = 0; m < nEquation; m++)
        {
            dDdP[m][idx] = (primR[m] - q[m][re]) / dqpert;
        }
    }
    if (nChemical == 1)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            qL[m] = q[m][le];
        }
        molecularWeight = gas->ComputeMolecularWeightReciprocal(qL);

        tGhost = t[ITT][le];
        pGhost = q[IP][le];
        //! Then compute the density from the extrapolated temperature.
        rGhost = pGhost / (coefficientofstateEquation * tGhost * molecularWeight);

        q[IR][re] = rGhost;
        q[IP][re] = pGhost;
    }

    delete [] primL;  primL = nullptr;
    delete [] primR;  primR = nullptr;
}
//! GMRESVis GMRESBoundary
void NSSolverUnstruct::Cal_GMRES_ViscousAdiabaticWallBCRegion(Grid *gridIn, RDouble *PrimL, RDouble *PrimR, int iFace, RDouble uWall, RDouble vWall, RDouble wWall)
{
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    RDouble *faceVelocityX = grid->GetFaceVelocityX();
    RDouble *faceVelocityY = grid->GetFaceVelocityY();
    RDouble *faceVelocityZ = grid->GetFaceVelocityZ();
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");
    int nEquation = GetNumberOfEquations();

    for (int m = 0; m < nEquation; ++ m)
    {
        PrimR[m] = PrimL[m];
    }

    RDouble velocityXWall = uWall;
    RDouble velocityYWall = vWall;
    RDouble velocityZWall = wWall;

    if (isUnsteady && isAle)
    {
        velocityXWall = faceVelocityX[iFace] + uWall;
        velocityYWall = faceVelocityY[iFace] + vWall;
        velocityZWall = faceVelocityZ[iFace] + wWall;
    }

    PrimR[IU] = - PrimL[IU] + two * velocityXWall;
    PrimR[IV] = - PrimL[IV] + two * velocityYWall;
    PrimR[IW] = - PrimL[IW] + two * velocityZWall;
}
//! AD  GMRESAD
void NSSolverUnstruct::GMRES_ViscousAdiabaticWallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    int nTotalCell         = grid->GetNTotalCell();
    RDouble *faceVelocityX = grid->GetFaceVelocityX();
    RDouble *faceVelocityY = grid->GetFaceVelocityY();
    RDouble *faceVelocityZ = grid->GetFaceVelocityZ();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));
    RDouble **dDdP = reinterpret_cast <RDouble**> (grid->GetDataPtr("dDdP"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();

    RDouble uWall = 0.0;
    RDouble vWall = 0.0;
    RDouble wWall = 0.0;

    if(bcData)
    {
        if (bcData->IsExist("uWall", PHDOUBLE, 1))
        {
            bcData->GetData("uWall", &uWall, PHDOUBLE, 1);
            bcData->GetData("vWall", &vWall, PHDOUBLE, 1);
            bcData->GetData("wWall", &wWall, PHDOUBLE, 1);
        }
    }

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    int iFace, le, re;
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    RDouble molecularWeight = one;
    using namespace GAS_SPACE;
    RDouble coefficientofstateEquation = gas->GetCoefficientOfStateEquation();
    RDouble rGhost = 0.0, pGhost = 0.0, tGhost = 0.0;
    ADReal *qL = new ADReal[nEquation]();
    ADReal *qR = new ADReal[nEquation]();

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        //! init for q[m][le] and define the independent variabels
        for (int m = 0; m < nEquation; m++)
        {
            qL[m] = q[m][le];
            qL[m].diff(m, nEquation);
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            qR[m] = qL[m];
        }

        ADReal velocityXWall = uWall;
        ADReal velocityYWall = vWall;
        ADReal velocityZWall = wWall;
        if (isUnsteady && isAle)
        {
            velocityXWall = faceVelocityX[iFace] + uWall;
            velocityYWall = faceVelocityY[iFace] + vWall;
            velocityZWall = faceVelocityZ[iFace] + wWall;
        }

        qR[IU] = - qL[IU] + two * velocityXWall;
        qR[IV] = - qL[IV] + two * velocityYWall;
        qR[IW] = - qL[IW] + two * velocityZWall;

        //! assign to qr and dqr/dql
        int idx = (re - nTotalCell) * nEquation;
        for (int m = 0; m < nEquation; m++)
        {
            q[m][re] = qR[m].val();
            for (int n = 0; n < nEquation; ++n)
            {
                dDdP[m][idx + n] = qR[m].dx(n);
            }
        }
    }

    if (nChemical == 1)
    {

    }
    delete [] qL;    qL = nullptr;
    delete [] qR;    qR = nullptr;
}
//! GMRESPV GMRESBoundary FD
void NSSolverUnstruct::GMRES_FarfieldBCRegion_FD(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    int nTotalCell = grid->GetNTotalCell();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));
    RDouble **dDdP = reinterpret_cast <RDouble**> (grid->GetDataPtr("dDdP"));
    int nEquation = GetNumberOfEquations();
    int iFace, le, re;
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();


    using namespace GAS_SPACE;

    RDouble *qL = new RDouble[nEquation]();
    RDouble *qR = new RDouble[nEquation]();
    RDouble dqpert;
    RDouble temperaturesL[3] = { 0 }, temperaturesR[3] = { 0 };
    RDouble perturbScale = 0.000001;

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            qL[m] = q[m][le];
        }
        Cal_GMRES_FarfieldBCRegion(qL, qR, le, iFace, gridIn, bcRegionUnstruct);
        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = qR[m];
        }

        //! perturb rho of the left cell
        qL[IR] = q[IR][le];
        qL[IU] = q[IU][le];
        qL[IV] = q[IV][le];
        qL[IW] = q[IW][le];
        qL[IP] = q[IP][le];
        qL[IR] *= (1.0+perturbScale);
        dqpert = qL[IR] - q[IR][le];
        Cal_GMRES_FarfieldBCRegion(qL, qR, le, iFace, gridIn, bcRegionUnstruct);
        dDdP[IR][nEquation * (re - nTotalCell) + IR] = (qR[IR] - q[IR][re]) / dqpert;
        dDdP[IU][nEquation * (re - nTotalCell) + IR] = (qR[IU] - q[IU][re]) / dqpert;
        dDdP[IV][nEquation * (re - nTotalCell) + IR] = (qR[IV] - q[IV][re]) / dqpert;
        dDdP[IW][nEquation * (re - nTotalCell) + IR] = (qR[IW] - q[IW][re]) / dqpert;
        dDdP[IP][nEquation * (re - nTotalCell) + IR] = (qR[IP] - q[IP][re]) / dqpert;

        //! perturb u of the left cell
        qL[IR] = q[IR][le];
        qL[IU] = q[IU][le];
        qL[IV] = q[IV][le];
        qL[IW] = q[IW][le];
        qL[IP] = q[IP][le];
        qL[IU] += (perturbScale);
        dqpert = qL[IU] - q[IU][le];
        Cal_GMRES_FarfieldBCRegion(qL, qR, le, iFace, gridIn, bcRegionUnstruct);
        dDdP[IR][nEquation * (re - nTotalCell) + IU] = (qR[IR] - q[IR][re]) / dqpert;
        dDdP[IU][nEquation * (re - nTotalCell) + IU] = (qR[IU] - q[IU][re]) / dqpert;
        dDdP[IV][nEquation * (re - nTotalCell) + IU] = (qR[IV] - q[IV][re]) / dqpert;
        dDdP[IW][nEquation * (re - nTotalCell) + IU] = (qR[IW] - q[IW][re]) / dqpert;
        dDdP[IP][nEquation * (re - nTotalCell) + IU] = (qR[IP] - q[IP][re]) / dqpert;

        //! perturb v of the left cell
        qL[IR] = q[IR][le];
        qL[IU] = q[IU][le];
        qL[IV] = q[IV][le];
        qL[IW] = q[IW][le];
        qL[IP] = q[IP][le];
        qL[IV] += (perturbScale);
        dqpert = qL[IV] - q[IV][le];
        Cal_GMRES_FarfieldBCRegion(qL, qR, le, iFace, gridIn, bcRegionUnstruct);
        dDdP[IR][nEquation * (re - nTotalCell) + IV] = (qR[IR] - q[IR][re]) / dqpert;
        dDdP[IU][nEquation * (re - nTotalCell) + IV] = (qR[IU] - q[IU][re]) / dqpert;
        dDdP[IV][nEquation * (re - nTotalCell) + IV] = (qR[IV] - q[IV][re]) / dqpert;
        dDdP[IW][nEquation * (re - nTotalCell) + IV] = (qR[IW] - q[IW][re]) / dqpert;
        dDdP[IP][nEquation * (re - nTotalCell) + IV] = (qR[IP] - q[IP][re]) / dqpert;

        //! perturb w of the left cell
        qL[IR] = q[IR][le];
        qL[IU] = q[IU][le];
        qL[IV] = q[IV][le];
        qL[IW] = q[IW][le];
        qL[IP] = q[IP][le];
        qL[IW] += (perturbScale);
        dqpert = qL[IW] - q[IW][le];
        Cal_GMRES_FarfieldBCRegion(qL, qR, le, iFace, gridIn, bcRegionUnstruct);
        dDdP[IR][nEquation * (re - nTotalCell) + IW] = (qR[IR] - q[IR][re]) / dqpert;
        dDdP[IU][nEquation * (re - nTotalCell) + IW] = (qR[IU] - q[IU][re]) / dqpert;
        dDdP[IV][nEquation * (re - nTotalCell) + IW] = (qR[IV] - q[IV][re]) / dqpert;
        dDdP[IW][nEquation * (re - nTotalCell) + IW] = (qR[IW] - q[IW][re]) / dqpert;
        dDdP[IP][nEquation * (re - nTotalCell) + IW] = (qR[IP] - q[IP][re]) / dqpert;

        //! perturb p of the left cell
        qL[IR] = q[IR][le];
        qL[IU] = q[IU][le];
        qL[IV] = q[IV][le];
        qL[IW] = q[IW][le];
        qL[IP] = q[IP][le];
        qL[IP] *= (1.0 + perturbScale);
        dqpert = qL[IP] - q[IP][le];
        Cal_GMRES_FarfieldBCRegion(qL, qR, le, iFace, gridIn, bcRegionUnstruct);
        dDdP[IR][nEquation * (re - nTotalCell) + IP] = (qR[IR] - q[IR][re]) / dqpert;
        dDdP[IU][nEquation * (re - nTotalCell) + IP] = (qR[IU] - q[IU][re]) / dqpert;
        dDdP[IV][nEquation * (re - nTotalCell) + IP] = (qR[IV] - q[IV][re]) / dqpert;
        dDdP[IW][nEquation * (re - nTotalCell) + IP] = (qR[IW] - q[IW][re]) / dqpert;
        dDdP[IP][nEquation * (re - nTotalCell) + IP] = (qR[IP] - q[IP][re]) / dqpert;     
    }
    delete [] qL;    qL = nullptr;
    delete [] qR;    qR = nullptr;
}

//! GMRESBoundary
void NSSolverUnstruct::GMRES_FarfieldBCRegion_CV(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    int nTotalCell = grid->GetNTotalCell();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));
    RDouble **dDdP = reinterpret_cast <RDouble**> (grid->GetDataPtr("dDdP"));
    int nEquation = GetNumberOfEquations();
    int iFace, le, re;
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    using namespace GAS_SPACE;

    RDouble *qL       = new RDouble[nEquation]();
    RDouble *qCvL     = new RDouble[nEquation]();
    RDouble *qCvLpert = new RDouble[nEquation]();
    RDouble *qR       = new RDouble[nEquation]();
    RDouble *qCvR     = new RDouble[nEquation]();
    RDouble *qCvRpert = new RDouble[nEquation]();
    RDouble dCvpert;
    RDouble temperaturesL[3] = { 0 }, temperaturesR[3] = { 0 };
    RDouble perturbScale = 0.001;

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            qL[m]= q[m][le];
        }

        gas->Primitive2Conservative(qL, gama[le], temperaturesL[ITV], temperaturesL[ITE], qCvL);

        Cal_GMRES_FarfieldBCRegion(qL,qR, le,iFace, gridIn, bcRegionUnstruct);

        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re]= qR[m];
        }

        gas->Primitive2Conservative(qR, gama[le], temperaturesR[ITV], temperaturesR[ITE], qCvR);

        //! perturb the density of the left cell
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IR] *= (1.0+perturbScale);
        dCvpert      = qCvLpert[IR] - qCvL[IR];
        gas->Conservative2Primitive(qCvLpert, gama[le], qL, temperaturesL);
        Cal_GMRES_FarfieldBCRegion(qL,qR, le,iFace, gridIn, bcRegionUnstruct);
        gas->Primitive2Conservative(qR, gama[le], temperaturesR[ITV], temperaturesR[ITE], qCvRpert);
        dDdP[IR][nEquation*(re-nTotalCell)+IR] = (qCvRpert[IR] - qCvR[IR])/dCvpert;
        dDdP[IU][nEquation*(re-nTotalCell)+IR] = (qCvRpert[IU] - qCvR[IU])/dCvpert;
        dDdP[IV][nEquation*(re-nTotalCell)+IR] = (qCvRpert[IV] - qCvR[IV])/dCvpert;
        dDdP[IW][nEquation*(re-nTotalCell)+IR] = (qCvRpert[IW] - qCvR[IW])/dCvpert;
        dDdP[IP][nEquation*(re-nTotalCell)+IR] = (qCvRpert[IP] - qCvR[IP])/dCvpert;

        //! perturb the momentum rhou of the left cell
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IU] += perturbScale;
        dCvpert      = qCvLpert[IU] - qCvL[IU];
        gas->Conservative2Primitive(qCvLpert, gama[le], qL, temperaturesL);
        Cal_GMRES_FarfieldBCRegion(qL,qR, le,iFace, gridIn, bcRegionUnstruct);
        gas->Primitive2Conservative(qR, gama[le], temperaturesR[ITV], temperaturesR[ITE], qCvRpert);
        dDdP[IR][nEquation*(re-nTotalCell)+IU] = (qCvRpert[IR] - qCvR[IR])/dCvpert;
        dDdP[IU][nEquation*(re-nTotalCell)+IU] = (qCvRpert[IU] - qCvR[IU])/dCvpert;
        dDdP[IV][nEquation*(re-nTotalCell)+IU] = (qCvRpert[IV] - qCvR[IV])/dCvpert;
        dDdP[IW][nEquation*(re-nTotalCell)+IU] = (qCvRpert[IW] - qCvR[IW])/dCvpert;
        dDdP[IP][nEquation*(re-nTotalCell)+IU] = (qCvRpert[IP] - qCvR[IP])/dCvpert;

        //! perturb the momentum rhov of the left cell
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IV] += perturbScale;
        dCvpert      = qCvLpert[IV] - qCvL[IV];
        gas->Conservative2Primitive(qCvLpert, gama[le], qL, temperaturesL);
        Cal_GMRES_FarfieldBCRegion(qL,qR, le,iFace, gridIn, bcRegionUnstruct);
        gas->Primitive2Conservative(qR, gama[le], temperaturesR[ITV], temperaturesR[ITE], qCvRpert);
        dDdP[IR][nEquation*(re-nTotalCell)+IV] = (qCvRpert[IR] - qCvR[IR])/dCvpert;
        dDdP[IU][nEquation*(re-nTotalCell)+IV] = (qCvRpert[IU] - qCvR[IU])/dCvpert;
        dDdP[IV][nEquation*(re-nTotalCell)+IV] = (qCvRpert[IV] - qCvR[IV])/dCvpert;
        dDdP[IW][nEquation*(re-nTotalCell)+IV] = (qCvRpert[IW] - qCvR[IW])/dCvpert;
        dDdP[IP][nEquation*(re-nTotalCell)+IV] = (qCvRpert[IP] - qCvR[IP])/dCvpert;

        //! perturb the momentum rhow of the left cell
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IW] += perturbScale;
        dCvpert      = qCvLpert[IW] - qCvL[IW];
        gas->Conservative2Primitive(qCvLpert, gama[le], qL, temperaturesL);
        Cal_GMRES_FarfieldBCRegion(qL,qR, le,iFace, gridIn, bcRegionUnstruct);
        gas->Primitive2Conservative(qR, gama[le], temperaturesR[ITV], temperaturesR[ITE], qCvRpert);
        dDdP[IR][nEquation*(re-nTotalCell)+IW] = (qCvRpert[IR] - qCvR[IR])/dCvpert;
        dDdP[IU][nEquation*(re-nTotalCell)+IW] = (qCvRpert[IU] - qCvR[IU])/dCvpert;
        dDdP[IV][nEquation*(re-nTotalCell)+IW] = (qCvRpert[IV] - qCvR[IV])/dCvpert;
        dDdP[IW][nEquation*(re-nTotalCell)+IW] = (qCvRpert[IW] - qCvR[IW])/dCvpert;
        dDdP[IP][nEquation*(re-nTotalCell)+IW] = (qCvRpert[IP] - qCvR[IP])/dCvpert;

        //! perturb the momentum rhoE of the left cell
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IP] *= (1.0+perturbScale);
        dCvpert      = qCvLpert[IP] - qCvL[IP];
        gas->Conservative2Primitive(qCvLpert, gama[le], qL, temperaturesL);
        Cal_GMRES_FarfieldBCRegion(qL,qR, le,iFace, gridIn, bcRegionUnstruct);
        gas->Primitive2Conservative(qR, gama[le], temperaturesR[ITV], temperaturesR[ITE], qCvRpert);
        dDdP[IR][nEquation*(re-nTotalCell)+IP] = (qCvRpert[IR] - qCvR[IR])/dCvpert;
        dDdP[IU][nEquation*(re-nTotalCell)+IP] = (qCvRpert[IU] - qCvR[IU])/dCvpert;
        dDdP[IV][nEquation*(re-nTotalCell)+IP] = (qCvRpert[IV] - qCvR[IV])/dCvpert;
        dDdP[IW][nEquation*(re-nTotalCell)+IP] = (qCvRpert[IW] - qCvR[IW])/dCvpert;
        dDdP[IP][nEquation*(re-nTotalCell)+IP] = (qCvRpert[IP] - qCvR[IP])/dCvpert;

    }
    delete [] qL;        qL       = nullptr; 
    delete [] qCvL;      qCvL     = nullptr;
    delete [] qCvLpert;  qCvLpert = nullptr;
    delete [] qR;        qR       = nullptr;
    delete [] qCvR;      qCvR     = nullptr;
    delete [] qCvRpert;  qCvRpert = nullptr;
}

// GMRESBoundary
template <typename T>
void NSSolverUnstruct::Cal_GMRES_FarfieldBCRegion(T* qL,T* qR,int le,int iFace, Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCelCal_GMRES_FarfieldBCRegionlofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *xfn         = grid->GetFaceNormalX();
    RDouble *yfn         = grid->GetFaceNormalY();
    RDouble *zfn         = grid->GetFaceNormalZ();
    RDouble *vgn         = grid->GetFaceNormalVelocity();
    RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));
    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    RDouble refGama = parameters->GetRefGama();

    int nEquation = GetNumberOfEquations();

    RDouble *primitiveVarFarfield = new RDouble [nEquation];
    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();
    bcData->GetData("primitiveVarFarfield", primitiveVarFarfield, PHDOUBLE, nEquation);

    T rInside, uInside, vInside, wInside, pInside, cInside, vNormalInside, vNormalOutside;
    T gamm1, velocityInside;
    T vTangentX, vTangentY, vTangentZ, entropyS;
    T R_Minus, R_Plus;
    T rBC, uBC, vBC, wBC, pBC, vNormalBC, cBC;
    T gamale = gama[le];

    //! Infinite variables.
    T rOutside = primitiveVarFarfield[IR];
    T uOutside = primitiveVarFarfield[IU];
    T vOutside = primitiveVarFarfield[IV];
    T wOutside = primitiveVarFarfield[IW];
    T pOutside = primitiveVarFarfield[IP];
    T cOutside = sqrt(abs(refGama * pOutside / rOutside));

    for (int m = 0; m < nEquation; ++ m)
    {
        qR[m]= qL[m];
    }
    //! Inner point variables.
    rInside = qL[IR];
    uInside = qL[IU];
    vInside = qL[IV];
    wInside = qL[IW];
    pInside = qL[IP];
    vNormalOutside = xfn[iFace] * uOutside + yfn[iFace] * vOutside + zfn[iFace] * wOutside - vgn[iFace];
    vNormalInside  = xfn[iFace] * uInside  + yfn[iFace] * vInside  + zfn[iFace] * wInside  - vgn[iFace];
    cInside = sqrt(abs(gamale * pInside / rInside));
    gamm1   = gamale - one;
    velocityInside = sqrt(uInside * uInside + vInside * vInside + wInside * wInside);

    //! Supersonic.
    if (velocityInside > cInside)
    {
        if (vNormalInside >= 0.0)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                qR[m] = qL[m];
            }
        }
        else
        {
            //! Inlet.
            qR[IR] = rOutside;
            qR[IU] = uOutside;
            qR[IV] = vOutside;
            qR[IW] = wOutside;
            qR[IP] = pOutside;
        }
        return;
    }
    //! Subsonic.
    R_Plus  = vNormalInside + 2.0 * cInside / gamm1;
    R_Minus = vNormalOutside - 2.0 * cOutside / gamm1;
    vNormalBC =  half * (R_Plus + R_Minus);
    cBC =  0.25 * (R_Plus - R_Minus) * gamm1;

    //! GMRESFarFieldCorrection
    T bf,vNth;
    vNth = 1e-4;
    bf  = cos((vNormalBC - vNth) / (2.0 * vNth)*PI);
    bf  = pow((bf+1)/2.0,2);

    if (vNormalBC > vNth)
    {
        entropyS   = pInside / pow(rInside, gama[le]);
        vTangentX  = uInside - xfn[iFace] * vNormalInside;
        vTangentY  = vInside - yfn[iFace] * vNormalInside;
        vTangentZ  = wInside - zfn[iFace] * vNormalInside;
    }
    else if (vNormalBC < -vNth)
    {
        //! Inlet.
        entropyS   = pOutside / pow(rOutside, gama[le]);
        vTangentX  = uOutside - xfn[iFace] * vNormalOutside;
        vTangentY  = vOutside - yfn[iFace] * vNormalOutside;
        vTangentZ  = wOutside - zfn[iFace] * vNormalOutside;
    }
    else
    {
        entropyS  = pInside / pow(rInside, gama[le]);
        vTangentX = uInside - xfn[iFace] * vNormalInside;
        vTangentY = vInside - yfn[iFace] * vNormalInside;
        vTangentZ = wInside - zfn[iFace] * vNormalInside;

        entropyS  = entropyS  *bf + (1 - bf)*(pOutside / pow(rOutside, gama[le]));
        vTangentX = vTangentX *bf + (1 - bf)*(uOutside - xfn[iFace] * vNormalOutside);
        vTangentY = vTangentY *bf + (1 - bf)*(vOutside - yfn[iFace] * vNormalOutside);
        vTangentZ = vTangentZ *bf + (1 - bf)*(wOutside - zfn[iFace] * vNormalOutside);
    }

    rBC = pow((cBC * cBC / (entropyS * gamale)), one / gamm1);
    uBC = vTangentX + xfn[iFace] * vNormalBC;
    vBC = vTangentY + yfn[iFace] * vNormalBC;
    wBC = vTangentZ + zfn[iFace] * vNormalBC;
    pBC = cBC * cBC * rBC / gamale;
    qR[IR] = rBC;
    qR[IU] = uBC;
    qR[IV] = vBC;
    qR[IW] = wBC;
    qR[IP] = pBC;
}

//! GMRESBoundary GMRESAD
void NSSolverUnstruct::GMRES_FarfieldBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    int nTotalCell = grid->GetNTotalCell();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));
    RDouble **dDdP = reinterpret_cast <RDouble**> (grid->GetDataPtr("dDdP"));
    int nEquation = GetNumberOfEquations();
    int iFace, le, re;
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    using namespace GAS_SPACE;

    ADReal *qL = new ADReal[nEquation]();
    ADReal *qR = new ADReal[nEquation]();

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            qL[m] = q[m][le];
            qL[m].diff(m, nEquation);
        }
        Cal_GMRES_FarfieldBCRegion(qL, qR, le, iFace, gridIn, bcRegionUnstruct);

        int idx = (re - nTotalCell) * nEquation;
        for (int m = 0; m < nEquation; ++m)
        {
            q[m][re] = qR[m].val();
            for (int n = 0; n < nEquation; ++n)
            {
                dDdP[m][idx + n] = qR[m].dx(n);
            }
        }
    }
    delete [] qL;  qL = nullptr;
    delete [] qR;  qR = nullptr;
}
#endif

void NSSolverUnstruct::ViscousIsotropicWallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme"); 
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *faceVelocityX = grid->GetFaceVelocityX();
    RDouble *faceVelocityY = grid->GetFaceVelocityY();
    RDouble *faceVelocityZ = grid->GetFaceVelocityZ();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **dDdP = reinterpret_cast<RDouble **>(grid->GetDataPtr("dDdP"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();
    int nTotalCells = grid->GetNTotalCell();

    int nEquation = GetNumberOfEquations();

    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();

    //! Get the variables defined in boundarycondition.hypara file.
    Data_Param *bcParamDB = bcRegionUnstruct->GetBCParamDataBase();
    RDouble wallTemperature = 0.0;

    if (bcParamDB)
    {
        if (!bcParamDB->IsExist("wallTemperature", PHDOUBLE, 1))
        {
            wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
            bcParamDB->UpdateData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
        }
        else
        {
            bcParamDB->GetData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
        }
    }
    else
    {
        wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
    }
    wallTemperature /= refDimensionalTemperature;

    bool wallTempArrayExist = false;
    RDouble *wallTempArray = nullptr;
    int countFace = 0;
    if (bcRegionUnstruct->CheckFieldData("wallTempArray"))
    {
        wallTempArrayExist = true;
        wallTempArray = reinterpret_cast <RDouble *> (bcRegionUnstruct->GetFieldDataPtr("wallTempArray"));
    }

    int iFace, le, re;
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    RDouble molecularWeight = one;
    using namespace GAS_SPACE;
    RDouble coefficientofstateEquation = gas->GetCoefficientOfStateEquation();

    if (tscheme != GMRES)
    {
        RDouble rGhost = 0.0, pGhost = 0.0, tGhost = 0.0;    //! rho/pressure/temprature on ghost cell.
    RDouble *qL = nullptr;
        if (nChemical == 1)
        {
            qL = new RDouble [nEquation];
        }

        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            iFace = *iter;
            le = leftCellofFace[iFace];
            re = rightCellofFace[iFace];

            for (int m = 0; m < nEquation; ++ m)
            {
                q[m][re] = q[m][le];
            }

            q[IU][re] = - q[IU][le] + two * faceVelocityX[iFace];
            q[IV][re] = - q[IV][le] + two * faceVelocityY[iFace];
            q[IW][re] = - q[IW][le] + two * faceVelocityZ[iFace];

            if (nChemical == 1)
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    qL[m] = q[m][le];
                }
                molecularWeight = gas->ComputeMolecularWeightReciprocal(qL);
            }

            RDouble tL = (1.0 / coefficientofstateEquation) * q[IP][le] / q[IR][le];

            //! Pressure keep the same on the wall normal direction.
            pGhost = q[IP][le];

            if (wallTempArrayExist)
            {
                wallTemperature = wallTempArray[countFace] / refDimensionalTemperature;
                countFace ++;
            }

            //! Extrapolate temperature from inside to outside.
            tGhost = two * wallTemperature - tL;
            if (tGhost <= wallTemperature*0.1) tGhost = wallTemperature*0.1;

            if (nChemical == 1)
            {
                tGhost = wallTemperature;
            }

            //! Then compute the density from the extrapolated temperature.
            rGhost = pGhost / (coefficientofstateEquation * tGhost * molecularWeight);

            q[IR][re] = rGhost;
            q[IP][re] = pGhost;
        }

        if (nChemical == 1)
        {
        delete [] qL;    qL = nullptr;
        }
    }
    else
    {
#ifdef USE_GMRESSOLVER

        ADReal rGhost  = 0.0, pGhost = 0.0, tGhost = 0.0;    //! rho/pressure/temprature on ghost cell.
        ADReal *qleft  = new ADReal [nEquation]();
        ADReal *qright = new ADReal [nEquation]();

        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            iFace = *iter;
            le = leftCellofFace[iFace];
            re = rightCellofFace[iFace];

            for (int m = 0; m < nEquation; ++ m)
            {
                qleft[m] = q[m][le];
                qleft[m].diff(m,nEquation);
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                qright[m] = qleft[m];
            }

            qright[IU] = - qleft[IU] + two * faceVelocityX[iFace];
            qright[IV] = - qleft[IV] + two * faceVelocityY[iFace];
            qright[IW] = - qleft[IW] + two * faceVelocityZ[iFace];


            ADReal tL = (1.0 / coefficientofstateEquation) * qleft[IP] / qleft[IR];

            //! Pressure keep the same on the wall normal direction.
            pGhost = qleft[IP];

            if (wallTempArrayExist)
            {
                wallTemperature = wallTempArray[countFace] / refDimensionalTemperature;
                countFace ++;
            }

            //! Extrapolate temperature from inside to outside.
            tGhost = two * wallTemperature - tL;
            if (tGhost <= 0.0) tGhost = wallTemperature;

            //! Then compute the density from the extrapolated temperature.
            rGhost = pGhost / (coefficientofstateEquation * tGhost * molecularWeight);

            qright[IR] = rGhost;
            qright[IP] = pGhost;

            int idx = (re - nTotalCells)*nEquation;

            for(int m = 0; m < nEquation; m++)
            {
                q[m][re] = qright[m].val();
                for(int n = 0; n < nEquation; n++)
                {
                    dDdP[m][idx+n] = qright[m].dx(n);
                }
            }
        }

        delete [] qleft;    qleft = nullptr;
        delete [] qright;    qright = nullptr;
#endif
    }
}

void NSSolverUnstruct::ViscousAdiabaticWallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme"); 
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    int nTotalCell         = grid->GetNTotalCell();
    RDouble *faceVelocityX = grid->GetFaceVelocityX();
    RDouble *faceVelocityY = grid->GetFaceVelocityY();
    RDouble *faceVelocityZ = grid->GetFaceVelocityZ();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));
    RDouble **dDdP = reinterpret_cast <RDouble**> (grid->GetDataPtr("dDdP"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();

    RDouble uWall = 0.0;
    RDouble vWall = 0.0;
    RDouble wWall = 0.0;

    if(bcData)
    {
        if (bcData->IsExist("uWall", PHDOUBLE, 1))
        {
            bcData->GetData("uWall", &uWall, PHDOUBLE, 1);
            bcData->GetData("vWall", &vWall, PHDOUBLE, 1);
            bcData->GetData("wWall", &wWall, PHDOUBLE, 1);
        }
    }

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    int iFace, le, re;
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    RDouble molecularWeight = one;
    using namespace GAS_SPACE;
    int nDensityForWallMethod = gas->GetnDensityForWallMethod();
    RDouble coefficientofstateEquation = gas->GetCoefficientOfStateEquation();
    RDouble rGhost = 0.0, pGhost = 0.0, tGhost = 0.0;

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");

    RDouble refGama = parameters->GetRefGama();
    RDouble gama1 = refGama - 1;

    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble coefficientOfStateEquation = 1.0 / (refGama * refMachNumber * refMachNumber);

    if (tscheme != GMRES)
    {
        RDouble *qL = nullptr;
        if (nChemical == 1)
        {
            qL = new RDouble [nEquation];
        }

        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            iFace = *iter;
            le = leftCellofFace[iFace];
            re = rightCellofFace[iFace];

            for (int m = 0; m < nEquation; ++ m)
            {
                q[m][re] = q[m][le];
            }

            RDouble velocityXWall = uWall;
            RDouble velocityYWall = vWall;
            RDouble velocityZWall = wWall;
            if (isUnsteady && isAle)
            {
                velocityXWall = faceVelocityX[iFace] + uWall;
                velocityYWall = faceVelocityY[iFace] + vWall;
                velocityZWall = faceVelocityZ[iFace] + wWall;
            }

            //! set wall velocity for translating frame and rotating frame.
            if (referenceFrame == TRANSLATIONAL_FRAME || referenceFrame == ROTATIONAL_FRAME)
            {
                int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");

                string shroud[100];
                GlobalDataBase::GetData("shroud", &shroud, PHSTRING, nTurboZone);

                velocityXWall = faceVelocityX[iFace];
                velocityYWall = faceVelocityY[iFace];
                velocityZWall = faceVelocityZ[iFace];

                //! Parallel
                int iTurboZone = grid->GetOrdinaryGridIndex();
                //! Serial
                if (iTurboZone == -1)
                {
                    iTurboZone = grid->GetZoneID();
                }

                if (bcRegionUnstruct->GetBCName() == shroud[iTurboZone])
                {
                    velocityXWall = 0;
                    velocityYWall = 0;
                    velocityZWall = 0;
                }
            }

            q[IU][re] = - q[IU][le] + two * velocityXWall;
            q[IV][re] = - q[IV][le] + two * velocityYWall;
            q[IW][re] = - q[IW][le] + two * velocityZWall;

            if (nChemical == 1)
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    qL[m] = q[m][le];
                }
                molecularWeight = gas->ComputeMolecularWeightReciprocal(qL);

                tGhost = t[ITT][le];
                pGhost = q[IP][le];
                //! Then compute the density from the extrapolated temperature.
                rGhost = pGhost / (coefficientofstateEquation * tGhost * molecularWeight);
                if (nDensityForWallMethod == 0)
                {
                    rGhost = MAX(0.5 * rGhost,2.0 * rGhost - q[IR][le]);
                }
                q[IR][re] = rGhost;
                q[IU][re] = - q[IU][le];
                q[IV][re] = - q[IV][le];
                q[IW][re] = - q[IW][le];
                q[IP][re] = pGhost;
            }
        }

        if (nChemical == 1)
        {
        delete [] qL;    qL = nullptr;
        }

    }
    else
    {
#ifdef USE_GMRESSOLVER
        ADReal *qL = new ADReal[nEquation]();
        ADReal *qR = new ADReal[nEquation]();

        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            iFace = *iter;
            le = leftCellofFace[iFace];
            re = rightCellofFace[iFace];

            //! init for q[m][le] and define the independent variabels
            for (int m = 0; m < nEquation; m++)
            {
                qL[m] = q[m][le];
                qL[m].diff(m, nEquation);
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                qR[m] = qL[m];
            }

            ADReal velocityXWall = uWall;
            ADReal velocityYWall = vWall;
            ADReal velocityZWall = wWall;
            if (isUnsteady && isAle)
            {
                velocityXWall = faceVelocityX[iFace] + uWall;
                velocityYWall = faceVelocityY[iFace] + vWall;
                velocityZWall = faceVelocityZ[iFace] + wWall;
            }

            qR[IU] = - qL[IU] + two * velocityXWall;
            qR[IV] = - qL[IV] + two * velocityYWall;
            qR[IW] = - qL[IW] + two * velocityZWall;

            //! assign to qr and dqr/dql
            int idx = (re - nTotalCell) * nEquation;
            for (int m = 0; m < nEquation; m++)
            {
                q[m][re] = qR[m].val();
                for (int n = 0; n < nEquation; ++n)
                {
                    dDdP[m][idx + n] = qR[m].dx(n);
                }
            }
        }

        delete [] qL;    qL = nullptr;
        delete [] qR;    qR = nullptr;
#endif
    }
}

void NSSolverUnstruct::FarfieldBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme"); 
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid    = UnstructGridCast(gridIn);
    int *leftCellofFace   = grid->GetLeftCellOfFace();
    int *rightCellofFace  = grid->GetRightCellOfFace();
    RDouble *wallDistance = grid->GetWallDist();
    int nTotalCell        = grid->GetNTotalCell();
    RDouble *xfn          = grid->GetFaceNormalX();
    RDouble *yfn          = grid->GetFaceNormalY();
    RDouble *zfn          = grid->GetFaceNormalZ();
    RDouble *vgn          = grid->GetFaceNormalVelocity();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));
    RDouble **dDdP = reinterpret_cast <RDouble**> (grid->GetDataPtr("dDdP"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nm       = parameters->GetNSEquationNumber();
    int nLaminar = parameters->GetLaminarNumber();

    RDouble refGama = parameters->GetRefGama();

    int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

    RDouble *preconCoefficient = nullptr;

    if (ifLowSpeedPrecon != 0)
    {
        preconCoefficient = reinterpret_cast<RDouble *> (grid->GetDataPtr("preconCoefficient"));
    }

    int nEquation = GetNumberOfEquations();

    RDouble *primitiveVarFarfield = new RDouble [nEquation];
    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();
    bcData->GetData("primitiveVarFarfield", primitiveVarFarfield, PHDOUBLE, nEquation);

    int iFace, le, re;
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    RDouble rInside, uInside, vInside, wInside, pInside, cInside, vNormalInside, vNormalOutside;
    RDouble gamm1, velocityInside;
    RDouble vTangentX, vTangentY, vTangentZ, entropyS;
    RDouble R_Minus, R_Plus;
    RDouble rBC, uBC, vBC, wBC, pBC, vNormalBC, cBC;

    //! Infinite variables.
    RDouble rOutside = primitiveVarFarfield[IR];
    RDouble uOutside = primitiveVarFarfield[IU];
    RDouble vOutside = primitiveVarFarfield[IV];
    RDouble wOutside = primitiveVarFarfield[IW];
    RDouble pOutside = primitiveVarFarfield[IP];
    RDouble cOutside = sqrt(ABS(refGama * pOutside / rOutside));

    if (tscheme != GMRES)
    {
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            iFace = *iter;
            le = leftCellofFace[iFace];
            re = rightCellofFace[iFace];

            for (int m = 0; m < nEquation; ++ m)
            {
                q[m][re] = q[m][le];
            }

            //! Inner point variables.
            rInside = q[IR][le];
            uInside = q[IU][le];
            vInside = q[IV][le];
            wInside = q[IW][le];
            pInside = q[IP][le];

            vNormalOutside = xfn[iFace] * uOutside + yfn[iFace] * vOutside + zfn[iFace] * wOutside - vgn[iFace];
            vNormalInside  = xfn[iFace] * uInside  + yfn[iFace] * vInside  + zfn[iFace] * wInside  - vgn[iFace];

            cInside = sqrt(ABS(gama[le] * pInside / rInside));
            gamm1   = gama[le] - one;

            velocityInside = sqrt(uInside * uInside + vInside * vInside + wInside * wInside);

            //! Supersonic.
            if (velocityInside > cInside)
            {
                if (vNormalInside >= 0.0)
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        q[m][re] = q[m][le];
                    }
                }
                else
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        q[m][re] = primitiveVarFarfield[m];
                    }
                }
            }
            else
            {
                //! Subsonic.
            if (ifLowSpeedPrecon == 0)
                {
                    R_Plus  = vNormalInside + 2.0 * cInside / gamm1;
                    R_Minus = vNormalOutside - 2.0 * cOutside / gamm1;
                    vNormalBC =  half * (R_Plus + R_Minus);
                    cBC =  0.25 * (R_Plus - R_Minus) * gamm1;
                    if (vNormalBC >= 0.0)
                    {
                        entropyS = pInside / pow(rInside, gama[le]);
                        vTangentX  = uInside - xfn[iFace] * vNormalInside;
                        vTangentY  = vInside - yfn[iFace] * vNormalInside;
                        vTangentZ  = wInside - zfn[iFace] * vNormalInside;
                    }
                    else
                    {
                        //! Inlet.
                        entropyS = pOutside / pow(rOutside, gama[le]);
                        vTangentX  = uOutside - xfn[iFace] * vNormalOutside;
                        vTangentY  = vOutside - yfn[iFace] * vNormalOutside;
                        vTangentZ  = wOutside - zfn[iFace] * vNormalOutside;
                    }

                    rBC = pow((cBC * cBC / (entropyS * gama[le])), one / gamm1);
                    uBC = vTangentX + xfn[iFace] * vNormalBC;
                    vBC = vTangentY + yfn[iFace] * vNormalBC;
                    wBC = vTangentZ + zfn[iFace] * vNormalBC;
                    pBC = cBC * cBC * rBC / gama[le];

                    q[IR][re] = rBC;
                    q[IU][re] = uBC;
                    q[IV][re] = vBC;
                    q[IW][re] = wBC;
                    q[IP][re] = pBC;

                    for (int m = nLaminar; m < nEquation; ++m)
                    {
                        q[m][re] = q[m][le];
                    }
                }
                else
                {
                    int preconFarfieldBCMethod = parameters->GetPreconFarFieldBCMethod();

                    if (preconFarfieldBCMethod == 0)
                    {
                        RDouble preconCoeff, B0, A0, cPrecondition;

                        preconCoeff = preconCoefficient[le];
                        cPrecondition = half * sqrt(((preconCoeff - 1) * vNormalInside) * ((preconCoeff - 1) * vNormalInside) + 4 * preconCoeff * cInside * cInside);

                        A0 = 2.0 * cInside/((-preconCoeff * vNormalInside * 0.5) + cPrecondition + 0.5 * vNormalInside) / gamm1;
                        B0 = 2.0 * cInside/((-preconCoeff * vNormalInside * 0.5) - cPrecondition + 0.5 * vNormalInside) / gamm1;
                        R_Plus = A0 * cInside + vNormalInside;
                        R_Minus = B0 * cOutside + vNormalOutside;
                        cBC  = (R_Plus - R_Minus) / (A0 - B0);
                        vNormalBC = R_Plus - A0 * cBC;

                        if (vNormalBC >= 0.0)
                        {
                            entropyS = pInside / pow(rInside, gama[le]);
                            vTangentX  = uInside - xfn[iFace] * vNormalInside;
                            vTangentY  = vInside - yfn[iFace] * vNormalInside;
                            vTangentZ  = wInside - zfn[iFace] * vNormalInside;
                        }
                        else
                        {
                            //! Inlet.
                            entropyS = pOutside / pow(rOutside, gama[le]);
                            vTangentX  = uOutside - xfn[iFace] * vNormalOutside;
                            vTangentY  = vOutside - yfn[iFace] * vNormalOutside;
                            vTangentZ  = wOutside - zfn[iFace] * vNormalOutside;
                        }

                        rBC = pow((cBC * cBC / (entropyS * gama[le])), one / gamm1);
                        uBC = vTangentX + xfn[iFace] * vNormalBC;
                        vBC = vTangentY + yfn[iFace] * vNormalBC;
                        wBC = vTangentZ + zfn[iFace] * vNormalBC;
                        pBC = cBC * cBC * rBC / gama[le];

                        q[IR][re] = rBC;
                        q[IU][re] = uBC;
                        q[IV][re] = vBC;
                        q[IW][re] = wBC;
                        q[IP][re] = pBC;
                    }
                    else if (preconFarfieldBCMethod == 1)
                    {
                        R_Plus  = vNormalInside + 2.0 * cInside / gamm1;
                        R_Minus = vNormalOutside - 2.0 * cOutside / gamm1;
                        vNormalBC =  half * (R_Plus + R_Minus);

                        if (vNormalBC >= 0.0)
                        {
                            q[IR][re] = rInside;
                            q[IU][re] = uInside;
                            q[IV][re] = vInside;
                            q[IW][re] = wInside;
                            q[IP][re] = pOutside;
                        }
                        else
                        {
                            //! Inlet.
                            q[IR][re] = rOutside;
                            q[IU][re] = uOutside;
                            q[IV][re] = vOutside;
                            q[IW][re] = wOutside;
                            q[IP][re] = pInside;
                        }
                    }
                }
                for (int m = nm; m < nEquation; ++ m)
                {
                    q[m][re] = q[m][le];
                }
            }
        }
    }
    else
    {
#ifdef USE_GMRESSOLVER
    using namespace GAS_SPACE;

    ADReal *qL = new ADReal[nEquation]();
    ADReal *qR = new ADReal[nEquation]();

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            qL[m] = q[m][le];
            qL[m].diff(m, nEquation);
        }
        Cal_GMRES_FarfieldBCRegion(qL, qR, le, iFace, gridIn, bcRegionUnstruct);

        int idx = (re - nTotalCell) * nEquation;
        for (int m = 0; m < nEquation; ++m)
        {
            q[m][re] = qR[m].val();
            for (int n = 0; n < nEquation; ++n)
            {
                dDdP[m][idx + n] = qR[m].dx(n);
            }
        }
    }
    delete [] qL;    qL = nullptr;
    delete [] qR;    qR = nullptr;
#endif
    }
    delete [] primitiveVarFarfield;    primitiveVarFarfield = nullptr;
}

void NSSolverUnstruct::InflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *wallDistance = grid->GetWallDist();

    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nEquation = GetNumberOfEquations();

    RDouble *primitiveVarFarfield = new RDouble [nEquation];
    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();
    bcData->GetData("primitiveVarFarfield", primitiveVarFarfield, PHDOUBLE, nEquation);

    int directionMethod = 0;
    if (bcData->IsExist("directionMethod", PHINT, 1))
    {
        bcData->GetData("directionMethod", &directionMethod, PHINT, 1);
    }

    RDouble localNonDimVelocity = 0.0;
    if (directionMethod == 2)
    {
        RDouble refDimensionalVelocity = parameters->GetReferenceDimensionalVelocity();
        RDouble localDimensionalVelocity;
        bcData->GetData("refDimensionalVelocity", &localDimensionalVelocity, PHINT, 1);
        localNonDimVelocity = localDimensionalVelocity / refDimensionalVelocity;
    }

    int iFace, le, re;
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];
        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = primitiveVarFarfield[m];
        }
        if (directionMethod == 2)
        {
            q[IU][re] = -1.0 * localNonDimVelocity * xfn[iFace];
            q[IV][re] = -1.0 * localNonDimVelocity * yfn[iFace];
            q[IW][re] = -1.0 * localNonDimVelocity * zfn[iFace];
        }
    }

    delete [] primitiveVarFarfield;    primitiveVarFarfield = nullptr;
}

void NSSolverUnstruct::PressureOutletBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble * > (grid->GetDataPtr("gama"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nEquation = GetNumberOfEquations();

    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();

    RDouble staticPressure;
    bcData->GetData("staticPressure", &staticPressure, PHDOUBLE, 1);

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");

    //! rotating frame, will move to other place next step.
    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");

    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();

    RDouble *faceVelocityX = grid->GetFaceVelocityX();
    RDouble *faceVelocityY = grid->GetFaceVelocityY();
    RDouble *faceVelocityZ = grid->GetFaceVelocityZ();

    //! Dimensional.
    RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
    RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble velocity = refMachNumber * refDimensionalSonicSpeed;
    RDouble dynamic_pressure = refDimensionalDensity * velocity * velocity;

    staticPressure /= dynamic_pressure;

    int iFace, le, re;
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];
        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = q[m][le];
        }
        q[IP][re] = staticPressure;
    }
}

void NSSolverUnstruct::MassFlowOutletBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    using namespace IDX;
    //! Dimensional.
    Param_NSSolverUnstruct *parameters = GetControlParameters();  
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
    RDouble refGama = parameters->GetRefGama();

    int nEquation = GetNumberOfEquations();

    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();
    RDouble massFluxRatio, massFlow, BCTotalArea;
    bcData->GetData("massFluxRatio", &massFluxRatio, PHDOUBLE, 1);
    bcData->GetData("massFlow", &massFlow, PHDOUBLE, 1);
    bcData->GetData("BCTotalArea", &BCTotalArea, PHDOUBLE, 1);

    //! Nondimensionalization.
    RDouble coefficientOfStateEquation = 1.0/(refGama * refMachNumber * refMachNumber);

    //! fluxDensity : rho * u;
    RDouble fluxDensity = massFlow/(refDimensionalDensity * refDimensionalDensity * refMachNumber)/BCTotalArea;

    RDouble mfTbCoeff = 0.5 * (refGama - 1.0) / refGama / coefficientOfStateEquation * 
        (1.0 - massFluxRatio * massFluxRatio);

    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    int iFace, le, re;

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = q[m][le];
        }

        //! calculate the temperature of the inside cell.
        RDouble staticTemperature = q[IP][le] / q[IR][le] / coefficientOfStateEquation;
        RDouble tb   = staticTemperature + mfTbCoeff * (q[IU][le] * q[IU][le] +
            q[IV][le] * q[IV][le] + q[IW][le] * q[IW][le]);

        if(tb < TINY) 
        {
            tb = staticTemperature;
        }

        RDouble pb = q[IP][le] * pow(tb / staticTemperature, refGama / (refGama - 1.0));
        RDouble rb = pb / tb / coefficientOfStateEquation;

        RDouble ub = q[IU][le] * massFluxRatio;
        RDouble vb = q[IV][le] * massFluxRatio;
        RDouble wb = q[IW][le] * massFluxRatio;

        //! judge the unusual situation, Vn < 0.
        RDouble Vn  = ub * xfn[iFace] + vb * yfn[iFace] + wb * zfn[iFace];        
        if(Vn < 0.0)
        {
            Vn = fluxDensity / rb;
            ub = Vn * xfn[iFace];
            vb = Vn * yfn[iFace];
            wb = Vn * zfn[iFace];
        }

        //! if the outlet is supersonic flow.
        RDouble cb = sqrt(refGama * coefficientOfStateEquation * tb);
        if(Vn >= cb)
        {
            rb = q[IR][le];
            ub = q[IU][le];
            vb = q[IV][le];
            wb = q[IW][le];
            pb = q[IP][le];
        }

        //! allocate value to the ghost cell.
        q[IR][re] = rb;
        q[IU][re] = ub;
        q[IV][re] = vb;
        q[IW][re] = wb;
        q[IP][re] = pb;
    }
}

void NSSolverUnstruct::calculateBCArea(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();

    //! BC total area does not change with calculating, so it only calculates once.
    if (bcData->IsExist("BCTotalArea", PHDOUBLE, 1))
    {
        return;
    }

    using namespace PHMPI;
    string bcName = bcRegionUnstruct->GetBCName();
    RDouble BCTotalArea = 0.0;

    //! One processor may owe some zones, so calculate all BC area in one processor.
    int nLocalZones = GetNumberofLocalZones();
    for (int iZone = 0; iZone < nLocalZones; iZone++)
    {
        int iZoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        Grid * iGrid = PHSPACE::GetGrid(iZoneID);
        UnstructGrid *grid = UnstructGridCast(iGrid);      
        
        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int nBCRegion = unstructBCSet->GetnBCRegion();

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
            string jBCName = bcRegion->GetBCName();

            //! Judge the same BC region according the BC name.
            if (jBCName == bcName)
            {
                int iFace;
                vector<int> *faceIndex = bcRegion->GetFaceIndex();
                RDouble   *area    = grid->GetFaceArea();
                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    iFace = *iter;
                    BCTotalArea += area[iFace];
                }
            }
        }
    }

    //! Get the total area in all processors with the same BC name.
    RDouble GLB_BCTotalArea;
    PH_AllReduce(&BCTotalArea, &GLB_BCTotalArea, 1, MPI_SUM);
    BCTotalArea = GLB_BCTotalArea;

    for (int iZone = 0; iZone < nLocalZones; iZone++)
    {
        int iZoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        Grid * iGrid = PHSPACE::GetGrid(iZoneID);
        UnstructGrid *grid = UnstructGridCast(iGrid); 

        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int nBCRegion = unstructBCSet->GetnBCRegion();

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
            string jBCName = bcRegion->GetBCName();

            if (jBCName == bcName)
            {
                Data_Param *bcDataTmp = bcRegion->GetBCParamDataBase();
                bcDataTmp->UpdateData("BCTotalArea", &BCTotalArea, PHDOUBLE, 1);
            }
        }
    }
}

void NSSolverUnstruct::CalculateMassFluxRatio()
{
    vector<SimpleBC*> *globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();
    if (globalBCList == 0)
    {
        return;
    }

    for (int iBC = 0; iBC < globalBCList->size(); ++iBC)
    {
        SimpleBC *boundaryCondition = (*globalBCList)[iBC];
        int bcType = boundaryCondition->GetBCType();

        if (bcType == PHENGLEI::MASS_FLOW_OUTLET)
        {
            RDouble massFlux = 0.0;
            string bcName = boundaryCondition->GetBCName();

            int nLocalZones = PHMPI::GetNumberofLocalZones();
            for (int iZone = 0; iZone < nLocalZones; iZone++)
            {
                int iZoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);
                Grid *iGrid = PHSPACE::GetGrid(iZoneID);
                UnstructGrid *grid = UnstructGridCast(iGrid);
                RDouble **q = reinterpret_cast<RDouble**> (grid->GetDataPtr("q"));

                UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
                int nBCRegion = unstructBCSet->GetnBCRegion();
                for (int iBCRegion = 0; iBCRegion < nBCRegion; ++iBCRegion)
                {
                    UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
                    string jBCName = bcRegion->GetBCName();

                    if (jBCName == bcName)
                    {
                        int iFace;
                        vector<int> *faceIndex = bcRegion->GetFaceIndex();
                        RDouble *area = grid->GetFaceArea();
                        int *leftCellofFace = grid->GetLeftCellOfFace();
                        RDouble *xfn = grid->GetFaceNormalX();
                        RDouble *yfn = grid->GetFaceNormalY();
                        RDouble *zfn = grid->GetFaceNormalZ();
                        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                        {
                            iFace = *iter;
                            int le = leftCellofFace[iFace];

                            RDouble vn_c1 = q[IU][le] * xfn[iFace] + q[IV][le] * yfn[iFace] + q[IW][le] * zfn[iFace];
                            massFlux += q[IR][le] * vn_c1 * area[iFace];
                        }
                    }
                }
            }

            RDouble GLB_MassFlux;
            PH_AllReduce(&massFlux, &GLB_MassFlux, 1, MPI_SUM);

            Param_NSSolverUnstruct *parameters = GetControlParameters();
            RDouble refMachNumber = parameters->GetRefMachNumber();
            RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
            RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();

            RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

            //! Calculte the dimensional mass flux in this BC name region.
            GLB_MassFlux *= refDimensionalDensity * refDimensionalSonicSpeed * refMachNumber * gridScaleFactor * gridScaleFactor;

            for (int iZone = 0; iZone < nLocalZones; iZone++)
            {
                int iZoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);
                Grid *iGrid = PHSPACE::GetGrid(iZoneID);
                UnstructGrid *grid = UnstructGridCast(iGrid);

                UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
                int nBCRegion = unstructBCSet->GetnBCRegion();

                for (int iBCRegion = 0; iBCRegion < nBCRegion; ++iBCRegion)
                {
                    UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
                    string jBCName = bcRegion->GetBCName();

                    if (jBCName == bcName)
                    {
                        Data_Param *bcData = bcRegion->GetBCParamDataBase();
                        RDouble massFlow = 0.0;
                        bcData->GetData("massFlow", &massFlow, PHDOUBLE, 1);

                        RDouble massFluxRatio = massFlow / (GLB_MassFlux + TINY);
                        massFluxRatio = MAX(0.8, MIN(massFluxRatio, 1.2));

                        bcData->UpdateData("massFluxRatio", &massFluxRatio, PHDOUBLE, 1);
                    }
                }
            }
        }
    }
}

void NSSolverUnstruct::MassFlowInletBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! get the basic mesh and flow field information.
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();

    int nEquation = GetNumberOfEquations();

    int iFace = 0, le = 0, re = 0;
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();

    //! Get the basic variables for mass flow in boundary.
    RDouble direction_inlet[3] = {0.0}, massFlow = 0.0, totalTemperature = 0.0;
    bcData->GetData("massFlow", &massFlow, PHDOUBLE, 1);
    bcData->GetData("totalTemperature", &totalTemperature, PHDOUBLE, 1);
    bcData->GetData("direction_inlet", direction_inlet, PHDOUBLE, 3);

    RDouble BCTotalArea = 0.0;
    if (bcData->IsExist("BCTotalArea", PHDOUBLE, 1))
    {
        bcData->GetData("BCTotalArea", &BCTotalArea, PHDOUBLE, 1);
    }else
    {
        TK_Exit::ExceptionExit("BCTotalArea data does not exist!\n");
    }

    massFlow /= BCTotalArea;

    //! Dimensional.
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
    RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
    RDouble refGama = parameters->GetRefGama();
    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();

    //! Nondimensionalization.
    totalTemperature /= refDimensionalTemperature;
    RDouble coefficientOfStateEquation = 1.0/(refGama * refMachNumber * refMachNumber);

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = primitiveVarFarfield[m];
        }

        q[IR][re]  = q[IR][le];
        RDouble Vnb = massFlow/ (refDimensionalDensity*refDimensionalSonicSpeed*refMachNumber) / q[IR][le];    //! massFlow

        q[IU][re] = Vnb * direction_inlet[0];
        q[IV][re] = Vnb * direction_inlet[1];
        q[IW][re] = Vnb * direction_inlet[2];

        //! square of the stagnation sound speed.
        RDouble c2   = refGama * coefficientOfStateEquation * totalTemperature; 
        //! square of the critical sound speed.
        RDouble c2_ct = 2.0 * c2 / (refGama+1.0); 
        //! square of the Laval number.
        RDouble lamd2 = Vnb * Vnb / c2_ct;
        RDouble tb    = totalTemperature * (1.0 - (refGama - 1.0)/(refGama + 1.0) * lamd2);
        q[IP][re] = q[IR][le] * coefficientOfStateEquation * tb;
    }
}

void NSSolverUnstruct::PressureInletBCRiemannRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble **q  = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nEquation = GetNumberOfEquations();

    int iFace = 0, le = 0, re = 0;
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();

    RDouble direction_inlet[3], totalPressure, totalTemperature;
    bcData->GetData("totalPressure", &totalPressure, PHDOUBLE, 1);
    bcData->GetData("totalTemperature", &totalTemperature, PHDOUBLE, 1);
    bcData->GetData("direction_inlet", direction_inlet, PHDOUBLE, 3);

    RDouble *primitiveVarFarfield = new RDouble [nEquation];
    bcData->GetData("primitiveVarFarfield", primitiveVarFarfield, PHDOUBLE, nEquation);

    int directionMethod;
    bcData->GetData("directionMethod", &directionMethod, PHINT, 1);

    //! Dimensional.
    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
    RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
    RDouble velocity = refMachNumber * refDimensionalSonicSpeed;
    RDouble dynamic_pressure = refDimensionalDensity * velocity * velocity;
    RDouble refGama = parameters->GetRefGama();
    RDouble gama1 = refGama - 1;

    //! Nondimensionalization.
    totalPressure /= dynamic_pressure;
    totalTemperature /= refDimensionalTemperature;
    RDouble coefficientOfStateEquation = 1.0/(refGama * refMachNumber * refMachNumber);

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = primitiveVarFarfield[m];
        }

        //! Get inner cell flow field information.
        RDouble rInside = q[IR][le];
        RDouble V2Inside = q[IU][le] * q[IU][le] + q[IV][le] * q[IV][le]  + q[IW][le] * q[IW][le];
        RDouble pInside = q[IP][le];

        //! Calculate the total energy, total enthalpy and sound speed.
        RDouble totalEnergy = pInside / (rInside * gama1) + 0.5 * V2Inside;
        RDouble totalEnthalpy = (refGama * coefficientOfStateEquation / gama1) * totalTemperature;
        RDouble c2 = refGama * pInside / rInside;

        //! Calculate the riemann invariants R+.
        RDouble Riemann = 2.0 * sqrt(c2) / gama1;
        Riemann += q[IU][le] * xfn[iFace] + q[IV][le] * yfn[iFace] + q[IW][le] * zfn[iFace];

        RDouble theta = xfn[iFace] * direction_inlet[0] + yfn[iFace] * direction_inlet[1] + 
            zfn[iFace] * direction_inlet[2];

        RDouble totalc2 = gama1 * (totalEnthalpy - (totalEnergy + pInside / rInside)+0.5 * V2Inside) + c2;

        RDouble a =  1.0 + 0.5 * gama1 * theta * theta;
        RDouble b = -1.0 * gama1 * theta * Riemann;
        RDouble c =  0.5 * gama1 * Riemann * Riemann - 2.0 * totalc2 / gama1;

        RDouble d = b * b - 4.0 * a * c;
        d = sqrt(max(zero, d));
        RDouble Vin   = (-b + d)/(2.0*a);
        Vin   = max(zero, Vin);
        V2Inside = Vin * Vin;

        c2 = totalc2 - 0.5 * gama1 * V2Inside;

        //! Fix the mach number in case of mach > 1.0.
        RDouble Mach2 = V2Inside / c2;
        Mach2 = min(one, Mach2);
        V2Inside = Mach2 * c2;
        Vin = sqrt(V2Inside);
        c2 = totalc2 - 0.5 * gama1 * V2Inside;

        RDouble tInside = c2 / (refGama * coefficientOfStateEquation);
        pInside = totalPressure * pow((tInside / totalTemperature), refGama / gama1);
        rInside = pInside/(coefficientOfStateEquation * tInside);

        q[IR][re] = rInside;
        q[IU][re] = Vin * direction_inlet[0];
        q[IV][re] = Vin * direction_inlet[1];
        q[IW][re] = Vin * direction_inlet[2];
        q[IP][re] = pInside;

        if (directionMethod == 2)
        {
            q[IU][re] = -1.0 * Vin * xfn[iFace];
            q[IV][re] = -1.0 * Vin * yfn[iFace];
            q[IW][re] = -1.0 * Vin * zfn[iFace];
        }
    }

    delete [] primitiveVarFarfield;    primitiveVarFarfield = nullptr;
}

void NSSolverUnstruct::SwapNeighborsDQNonblocking(Grid *gridIn, FieldProxy *dqProxy)
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int nZones = PHMPI::GetNumberofGlobalZones();
    vector< PH_Request > requestContainer;
    vector< vector< DataContainer * > > receivedDataBuffer;
    vector< vector< DataContainer * > > sendDataBuffer;
    receivedDataBuffer.resize(nZones);
    sendDataBuffer.resize(nZones);
    int level = gridIn->GetLevel();

    //! Step 0: Compressing data firstly.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcessor = PHMPI::GetZoneProcessorID(iZone);

        ZoneNeighbor * globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
        if(currentProcessor == sendProcessor)
        {
            sendDataBuffer[iZone].resize(numberOfNeighbor);
        }

        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if(currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate and compress the buffers for sending.
            DataContainer * sendBuffer = 0;
            if(currentProcessor == sendProcessor)
            {
                sendBuffer = new DataContainer;
                sendDataBuffer[ iZone ][iNeighbor] = sendBuffer;

                //! Compress the send information into the actkey.
                UnstructGrid * grid = UnstructGridCast(PHSPACE::GetGrid(iZone, level));
                CompressDQ(sendBuffer, dqProxy, grid, iZone, neighborZone, 5);
            }
        }
    }

    //! Step 1: Communication.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor * globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        //! Communicating.
        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessor    = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if(currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate the buffers for sending and receiving.
            DataContainer * receivedBuffer = 0;
            DataContainer * sendBuffer = 0;
            if(currentProcessor == receiveProcessor)
            {
                receivedBuffer = new DataContainer;
                receivedDataBuffer[ iZone ].push_back(receivedBuffer);
            }
            if(currentProcessor == sendProcessor)
            {
                sendBuffer = sendDataBuffer[iZone][iNeighbor];
            }

            int tag = iZone;

            //if (sendProcessor == receiveProcessor) continue;

            //! Communication.
            if (currentProcessor == sendProcessor)
            {
                if (sendProcessor != receiveProcessor)
                {
                    streamsize nlen = sendBuffer->Size();

                    //! Send the data to neighbors.
                    send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
                }
            }
            if (currentProcessor == receiveProcessor)
            {
                //! Get the neighbor order: which order that 'iZone' on neighbors of neighborZone.
                int neighborOrer = -1;
                ZoneNeighbor * globalNeighborZonesTemp = zoneConnectivity->GetZoneNeighbor(neighborZone);
                std::size_t numberOfNeighborTemp = globalNeighborZonesTemp->GetNumberOfNeighbor();

                for (std::size_t neighborID = 0; neighborID < numberOfNeighborTemp; ++ neighborID)
                {
                    int neighborZoneTemp = globalNeighborZonesTemp->GetZoneIndexOfNeighbor(neighborID);
                    if(neighborZoneTemp == iZone)
                    {
                        neighborOrer = static_cast<int>(neighborID);
                        break;
                    }
                }

                ASSERT(neighborOrer != -1);

                if (sendProcessor != receiveProcessor)
                {
                    //! The data length of the received data is same to the length that send to the neighbor.
                    CharVecSizeType nlen = sendDataBuffer[neighborZone][neighborOrer]->Size();

                    //! Receive data from neighbors.
                    receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
                }
                else
                {
                    SWAP(receivedDataBuffer[ iZone ].back(), sendDataBuffer[iZone][iNeighbor]);
                }
            }
        }
    }

    //! Step2: MPI waiting for Non-blocking communications.
    if(PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_Wait(static_cast<int>(requestContainer.size()), &(requestContainer[0]));
    }

    //! Step3: Translate the data container.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor * globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor == receiveProcessor)
            {
                DataContainer * receiveData = receivedDataBuffer[ iZone ][ count ];

                //! Because get into here when currentProcessor==receiveProcessor,
                //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!

                //! Decompress the interface data from data-container.
                UnstructGrid * grid = UnstructGridCast(PHSPACE::GetGrid(neighborZone, level));
                DecompressDQ(receiveData, dqProxy, grid, neighborZone, iZone, 5);

                ++ count;
            }
        }
    }

    //! Step4: Free the buffers.
    for (std::size_t iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < receivedDataBuffer[ iDim ].size(); ++ jDim)
        {
            delete receivedDataBuffer[ iDim ][ jDim ];
        }
    }
    receivedDataBuffer.clear();
    for (std::size_t iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < sendDataBuffer[ iDim ].size(); ++ jDim)
        {
            delete sendDataBuffer[ iDim ][ jDim ];
        }
    }
    sendDataBuffer.clear();
    requestContainer.clear();
}

void NSSolverUnstruct::CompressGradientandLimitToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **limit2D = reinterpret_cast <RDouble **> (grid->GetDataPtr("limit2D"));

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
    dataContainer->MoveToBegin();
    
    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble **dtdx  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureX"));
    RDouble **dtdy  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureY"));
    RDouble **dtdz  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureZ"));
    
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int sourceCell;
        int iFace = interfaceIndexContainerForSend[iLocalFace];
        grid->GetSourceIndex(iFace, 1, sourceCell);
        for (int m = 0; m <= 4; ++ m)
        {
            PHWrite(dataContainer, dqdx[m][sourceCell]);
            PHWrite(dataContainer, dqdy[m][sourceCell]);
            PHWrite(dataContainer, dqdz[m][sourceCell]);
        }

        PHWrite(dataContainer, dtdx[0][sourceCell]);
        PHWrite(dataContainer, dtdy[0][sourceCell]);
        PHWrite(dataContainer, dtdz[0][sourceCell]);

        PHWrite(dataContainer, limit2D[0][sourceCell]);
        PHWrite(dataContainer, limit2D[1][sourceCell]);
        PHWrite(dataContainer, limit2D[2][sourceCell]);
        PHWrite(dataContainer, limit2D[3][sourceCell]);
        PHWrite(dataContainer, limit2D[4][sourceCell]);
    }
}

void NSSolverUnstruct::DecompressGradientandLimitToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **limit2D = reinterpret_cast <RDouble **> (grid->GetDataPtr("limit2D"));

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);
    dataContainer->MoveToBegin();
    
    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble **dtdx  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureX"));
    RDouble **dtdy  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureY"));
    RDouble **dtdz  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureZ"));

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];
        int targetCell;
        grid->GetTargetIndex(iFace, 1, targetCell);
        for (int m = 0; m <= 4; ++ m)
        {
            PHRead(dataContainer, dqdx[m][targetCell]);
            PHRead(dataContainer, dqdy[m][targetCell]);
            PHRead(dataContainer, dqdz[m][targetCell]);
        }

        PHRead(dataContainer, dtdx[0][targetCell]);
        PHRead(dataContainer, dtdy[0][targetCell]);
        PHRead(dataContainer, dtdz[0][targetCell]);

        PHRead(dataContainer, limit2D[0][targetCell]);
        PHRead(dataContainer, limit2D[1][targetCell]);
        PHRead(dataContainer, limit2D[2][targetCell]);
        PHRead(dataContainer, limit2D[3][targetCell]);
        PHRead(dataContainer, limit2D[4][targetCell]);
    }
}

void NSSolverUnstruct::CompressQlQrToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
    int * interFace2BoundaryFace = interfaceInformation->GetInterFace2BoundaryFace();
    dataContainer->MoveToBegin();

    RDouble **qL = InviscidfaceProxyOnlyForMixGrid->GetQL();

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForSend[iLocalFace];
        int FaceID = interFace2BoundaryFace[iFace];
        
        for (int m = 0; m < 5; ++ m)
        {
            PHWrite(dataContainer, qL[m][FaceID]);
        }
    }
}

void NSSolverUnstruct::DecompressQlQrToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);
    int * interFace2BoundaryFace = interfaceInformation->GetInterFace2BoundaryFace();
    dataContainer->MoveToBegin();

    RDouble **qR = InviscidfaceProxyOnlyForMixGrid->GetQR();

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];
        int FaceID = interFace2BoundaryFace[iFace];
        
        for (int m = 0; m < 5; ++ m)
        {
            PHRead(dataContainer, qR[m][FaceID]);
        }
    }
}

void NSSolverUnstruct::CompressFacefluxToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int gridtypeofcurrentZone = grid->Type();
    int gridtypeofneighborZone = PHMPI::GetZoneGridType(neighborZoneIndex);
    if(gridtypeofcurrentZone == gridtypeofneighborZone)
    {
        RDouble uploadData = 0.0;
        dataContainer->MoveToBegin();
        PHWrite(dataContainer, uploadData);
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
    int * interFace2BoundaryFace = interfaceInformation->GetInterFace2BoundaryFace();
    dataContainer->MoveToBegin();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    if(viscousType <= INVISCID)
    {
        //! Euler
        RDouble **Inviscidflux = InviscidfaceProxyOnlyForMixGrid->GetFlux();

        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int iFace = interfaceIndexContainerForSend[iLocalFace];
            int FaceID = interFace2BoundaryFace[iFace];
            
            for (int m = 0; m < 5; ++ m)
            {
                PHWrite(dataContainer, Inviscidflux[m][FaceID]);
            }
        }
    }
    else
    {
        //! laminar
        RDouble **Inviscidflux = InviscidfaceProxyOnlyForMixGrid->GetFlux();
        RDouble **Viscousflux = ViscousfaceProxyOnlyForMixGrid->GetFlux();
        
        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int iFace = interfaceIndexContainerForSend[iLocalFace];
            int FaceID = interFace2BoundaryFace[iFace];
            
            for (int m = 0; m < 5; ++ m)
            {
                PHWrite(dataContainer, Inviscidflux[m][FaceID]);
                PHWrite(dataContainer, Viscousflux[m][FaceID]);
            }
        }
    }
}

void NSSolverUnstruct::DecompressFacefluxToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int gridtypeofcurrentZone = grid->Type();
    int gridtypeofneighborZone = PHMPI::GetZoneGridType(neighborZoneIndex);
    if(gridtypeofcurrentZone == gridtypeofneighborZone)
    {
        RDouble downdData = 0.0;
        dataContainer->MoveToBegin();
        PHRead(dataContainer, downdData);
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);
    int * interFace2BoundaryFace = interfaceInformation->GetInterFace2BoundaryFace();
    dataContainer->MoveToBegin();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    RDouble coefofstrflux = GlobalDataBase::GetDoubleParaFromDB("coefofstrflux");
    RDouble coefofunstrflux = 1.0 - coefofstrflux;

    if(viscousType <= INVISCID)
    {
        //! Euler
        RDouble **Inviscidflux = InviscidfaceProxyOnlyForMixGrid->GetFlux();

        RDouble Inviscidflux_str[5] = {0.0};
        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int iFace = interfaceIndexContainerForReceive[iLocalFace];
            int FaceID = interFace2BoundaryFace[iFace];
            
            for (int m = 0; m < 5; ++ m)
            {
                PHRead(dataContainer, Inviscidflux_str[m]);

                Inviscidflux[m][FaceID] = (coefofstrflux * Inviscidflux_str[m] + coefofunstrflux * Inviscidflux[m][FaceID]);
            }
        }
    }
    else
    {
        //! laminar
        RDouble **Inviscidflux = InviscidfaceProxyOnlyForMixGrid->GetFlux();
        RDouble **Viscousflux = ViscousfaceProxyOnlyForMixGrid->GetFlux();

        RDouble Inviscidflux_str[5];
        RDouble Viscousflux_str[5];
        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int iFace = interfaceIndexContainerForReceive[iLocalFace];
            int FaceID = interFace2BoundaryFace[iFace];
            
            for (int m = 0; m < 5; ++ m)
            {
                PHRead(dataContainer, Inviscidflux_str[m]);
                PHRead(dataContainer, Viscousflux_str[m]);

                Inviscidflux[m][FaceID] = (coefofstrflux * Inviscidflux_str[m] + coefofunstrflux * Inviscidflux[m][FaceID]);
                Viscousflux[m][FaceID]  = (coefofstrflux * Viscousflux_str[m]  + coefofunstrflux * Viscousflux[m][FaceID]);
            }
        }
    }
}

#ifdef USE_GMRESSOLVER
//! GMRES solver -- a linear system solver
void NSSolverUnstruct::GMRESSolverSingle(Grid *gridIn, FieldProxy *dqProxy)
{
    //! GMRESCoupled
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    if ( viscousType == ONE_EQU )
    {
        return;
    }

    UnstructGrid *grid = UnstructGridCast(gridIn);
    //! Memory allocating
    RDouble **dRdq = reinterpret_cast<RDouble**>(grid->GetDataPtr("dRdq"));
    RDouble **dDdP = reinterpret_cast<RDouble**>(grid->GetDataPtr("dDdP"));

    int nTotalCells = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nEquations = GetNumberOfEquations(); 
        
    RDouble **res = reinterpret_cast<RDouble **>(grid->GetDataPtr("res"));
    RDouble *vol = grid->GetCellVolume();
    RDouble *dt = reinterpret_cast<RDouble*> (grid->GetDataPtr("dt"));
    vector<int> AI = grid->GetJacobianAI4GMRES();    //! GMRESSparse GMRESCSR
    vector<int> AJ = grid->GetJacobianAJ4GMRES();    //! GMRESSparse GMRESCSR

    //! GMRESCSR
    for (int iCell = 0; iCell < nTotalCells; iCell++)
    {
        vector<int>::iterator result = find(AJ.begin() + AI[iCell], AJ.begin() + AI[iCell + 1], iCell);    //! find qcell
        int index = distance(AJ.begin(), result);
        index *= nEquations;
        for (int iEquation = 0; iEquation < nEquations; iEquation++)
        {
             dRdq[iEquation][index + iEquation] += 1.0 / dt[iCell];
        }
    }

    //! GMRESJac1st
    vector<int> AI1st = grid->GetJacobianAI1st4GMRES(); 
    vector<int> AJ1st = grid->GetJacobianAJ1st4GMRES();
    RDouble **dRdq1st = reinterpret_cast<RDouble**>(grid->GetDataPtr("dRdq1st"));
    int jacOrder    = grid->GetJacobianOrder();
    if(jacOrder == 2)
    {
        for (int iCell = 0; iCell < nTotalCells; iCell++)
        {
            vector<int>::iterator result = find(AJ1st.begin() + AI1st[iCell], AJ1st.begin() + AI1st[iCell + 1], iCell); // find qcell
            int index = distance(AJ1st.begin(), result);
            index *= nEquations;
            for (int iEquation = 0; iEquation < nEquations; iEquation++)
            {
                 dRdq1st[iEquation][index + iEquation] += 1.0 / dt[iCell];
            }
        }
    }

    //! Entry to PETSc for NS equations.
    Vec x, b;
    Mat A;
    Mat A1st;                       //! GMRESJac1st
    KSP ksp;                        //! define solver
    PC pc;                          //! define matrix
    PetscInt nTotalSize;            //! dimension size of matrix and vector
    PetscInt nLevels = 3;           //! define number of levels used by preconditioner ILU
    PetscInt maxStepsRestart = 500; //! define maximum steps required by [restart]
    PetscInt maxSteps = 500;        //! define maximum steps
    PetscReal rTolerance = 1.e-1;   //! define tolerance of relative error // -1
    PetscReal dropTolerance = 1.e5; //! define tolerance of drop error

    nTotalSize = nEquations * nTotalCells;

    //! initialize PETSc
    //! PetscFunctionBeginUser;
    PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);

    //! create vectors;
    VecCreate(PETSC_COMM_SELF, &x);
    PetscObjectSetName((PetscObject)x, "solution");
    VecSetSizes(x, PETSC_DECIDE, nTotalSize);
    VecSetUp(x);
    VecDuplicate(x, &b);

    //! create matrix
    //! GMRESCSR
    PetscInt nnz[nTotalCells];
    for (int iCell = 0; iCell < nTotalCells; iCell++)
    {
        int count = 0;
        for (int j = AI[iCell]; j < AI[iCell + 1]; j++)
        {
            if(AJ[j] < nTotalCells)
            {
                count++;
            }
        }
        nnz[iCell] = count;
    }
    MatCreate(PETSC_COMM_SELF, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, nTotalSize, nTotalSize);
    MatSetType(A, MATBAIJ);
    MatSetBlockSize(A, nEquations);
    MatSeqBAIJSetPreallocation(A, nEquations, PETSC_DEFAULT, nnz);
    MatSetUp(A);

    //! GMRESJac1st
    if( jacOrder == 2)
    {
        PetscInt nnz1st[nTotalCells];
        for (int iCell = 0; iCell < nTotalCells; iCell++)
        {
            int count = 0;
            for (int j = AI1st[iCell]; j < AI1st[iCell + 1]; j++)
            {
                if(AJ1st[j] < nTotalCells)
                {
                    count++;
                }
            }
            nnz1st[iCell] = count;
        }
        MatCreate(PETSC_COMM_SELF, &A1st);
        MatSetSizes(A1st, PETSC_DECIDE, PETSC_DECIDE, nTotalSize, nTotalSize);
        MatSetType(A1st, MATBAIJ);
        MatSetBlockSize(A1st, nEquations);
        MatSeqBAIJSetPreallocation(A1st, nEquations, PETSC_DEFAULT, nnz1st);
        MatSetUp(A1st);
    }

    //! improved approach to assemble matrix CSR
    PetscInt petsccol[nEquations];
    PetscInt petscrow;
    int jmax = AI[nTotalCells];
    int rcell = 0;
    int j = AI[rcell];
    do
    {
        int rowidx = rcell * nEquations;
        for (j = AI[rcell]; j < AI[rcell + 1]; j++)
        {
            int qcell = AJ[j];    //! get qcell (colume index of matrix)
            //! start to insert values: limited to dimension of nTotalSize
            if(qcell < nTotalCells)
            {
                int colidx = qcell * nEquations;

                for (int irow = 0; irow < nEquations; irow++)
                {
                    petscrow = rowidx + irow;
                   
                    for (int n = 0; n < nEquations; n++)
                    {
                        petsccol[n] = colidx + n;
                    }
                    MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdq[irow][j * nEquations], INSERT_VALUES);
                }
            }
        }
        rcell++; //! update rcell, go to the next row index of matrix
    } while (j < jmax);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    if(jacOrder == 2)
    {
        int jmax = AI1st[nTotalCells];
        int rcell = 0;
        int j = AI1st[rcell];
        do
        {
            int rowidx = rcell * nEquations;
            for (j = AI1st[rcell]; j < AI1st[rcell + 1]; j++)
            {
                int qcell = AJ1st[j];    //! get qcell (colume index of matrix)
                //! start to insert values: limited to dimension of nTotalSize
                if(qcell < nTotalCells)
                {
                    int colidx = qcell * nEquations;

                    for (int irow = 0; irow < nEquations; irow++)
                    {
                        petscrow = rowidx + irow;
                       
                        for (int n = 0; n < nEquations; n++)
                        {
                            petsccol[n] = colidx + n;
                        }
                        MatSetValues(A1st, 1, &petscrow, nEquations, petsccol, &dRdq1st[irow][j * nEquations], INSERT_VALUES);
                    }
                }
            }
            rcell++;    //! update rcell, go to the next row index of matrix
        } while (j < jmax);
        MatAssemblyBegin(A1st, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(A1st, MAT_FINAL_ASSEMBLY);
    }

    //! assemble vector
    VecSet(b, zero);
    for (PetscInt iCell = 0; iCell < nTotalCells; iCell++)
    {
        for (PetscInt iEquation = 0; iEquation < nEquations; iEquation++)
        {
            PetscInt isize = iCell * nEquations + iEquation;
            VecSetValues(b, 1, &isize, &res[iEquation][iCell], ADD_VALUES);
        }
    }
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    //! define solver
    KSPCreate(PETSC_COMM_SELF, &ksp);
    KSPSetType(ksp, KSPGMRES);
    if(jacOrder == 2)
    {
        KSPSetOperators(ksp, A, A1st);
    }
    else
    {
        KSPSetOperators(ksp, A, A);
    }
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCILU);
    KSPSetPCSide(ksp,PC_RIGHT);
    PCFactorSetLevels(pc, nLevels);
    KSPGMRESSetRestart(ksp, maxStepsRestart);
    KSPGMRESSetOrthogonalization(ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);
    KSPSetTolerances(ksp, rTolerance, 1.e-20, dropTolerance, maxSteps);

    //! solution
    KSPSolve(ksp, b, x);

    //! GMRESResidual
    RDouble **dq = dqProxy->GetField_UNS();

    //! convert x to dqProxy
    for (PetscInt iCell = 0; iCell < nTotalCells; iCell++)
    {
        for (PetscInt iEquation = 0; iEquation < nEquations; iEquation++)
        {
            PetscInt isize = iCell * nEquations + iEquation;
            PetscReal value;
            VecGetValues(x, 1, &isize, &dq[iEquation][iCell]);
        }
    }

    //! free work space
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    if( jacOrder == 2)
    {
        MatDestroy(&A1st);
    }
    KSPDestroy(&ksp);
    PetscFinalize();
    return;
}

//! GMRESParallel
//! GMRES solver -- a linear system solver for parallel
void NSSolverUnstruct::GMRESSolverRow(Grid *gridIn, FieldProxy *dqProxy)
{
    //! GMRESCoupled
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    if ( viscousType == ONE_EQU )
    {
        return;
    }

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int globalCellIndexShift = grid->GetGlobalCellIndexShift();
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    //! Memory allocating
    RDouble **dRdq = reinterpret_cast<RDouble**>(grid->GetDataPtr("dRdq"));
    RDouble **dDdP = reinterpret_cast<RDouble**>(grid->GetDataPtr("dDdP"));

    // RDouble **dRdqCoupledTerm = reinterpret_cast<RDouble**>(grid->GetDataPtr("dRdqCoupledTerm"));

    int nTotalCells = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nEquations = GetNumberOfEquations(); 

    RDouble **res = reinterpret_cast<RDouble **>(grid->GetDataPtr("res"));
    RDouble *vol = grid->GetCellVolume();
    RDouble *dt = reinterpret_cast<RDouble*> (grid->GetDataPtr("dt"));
    vector<int> AI = grid->GetJacobianAI4GMRES();    //! GMRESSparse GMRESCSR
    vector<int> AJ = grid->GetJacobianAJ4GMRES();    //! GMRESSparse GMRESCSR
    vector<int> AK = grid->GetJacobianAK4GMRES();    //! GMRES Parallel

    for (int iCell = 0; iCell < nTotalCells; iCell++)
    {
        vector<int>::iterator result = find(AJ.begin() + AI[iCell], AJ.begin() + AI[iCell + 1], iCell);    //! find qcell
        int index = distance(AJ.begin(), result);
        index *= nEquations;
        for (int iEquation = 0; iEquation < nEquations; iEquation++){
             dRdq[iEquation][index + iEquation] += 1.0 / dt[iCell];
        }
    }

    //! GMRESJac1st
    vector<int> AI1st = grid->GetJacobianAI1st4GMRES(); 
    vector<int> AJ1st = grid->GetJacobianAJ1st4GMRES();
    RDouble **dRdq1st = reinterpret_cast<RDouble**>(grid->GetDataPtr("dRdq1st"));
    int jacOrder    = grid->GetJacobianOrder();
    if(jacOrder == 2)
    {
        for (int iCell = 0; iCell < nTotalCells; iCell++)
        {
            vector<int>::iterator result = find(AJ1st.begin() + AI1st[iCell], AJ1st.begin() + AI1st[iCell + 1], iCell); // find qcell
            int index = distance(AJ1st.begin(), result);
            index *= nEquations;
            for (int iEquation = 0; iEquation < nEquations; iEquation++)
            {
                 dRdq1st[iEquation][index + iEquation] += 1.0 / dt[iCell];
            }
        }
    }
    //! entry to PETSc for NS equations
    Vec x, b;
    Mat A;
    Mat APC;    //! GMRES Parallel
    KSP ksp;                        //! define solver
    PC pc;                          //! define matrix
    PetscInt nTotalSize;            //! dimension size of matrix and vector
    PetscInt nLevels = 3;           //! define number of levels used by preconditioner ILU
    PetscInt maxStepsRestart = 500; //! define maximum steps required by [restart]
    PetscInt maxSteps = 500;        //! define maximum steps
    PetscReal rTolerance = 1.e-6;   //! define tolerance of relative error // -1
    PetscReal dropTolerance = 1.e5; //! define tolerance of drop error    
    PetscMPIInt size;

    int globalTotalCells;
    globalTotalCells = (int)GlobalDataBase::GetDoubleParaFromDB("GlobalTotalCells");
    int nInterfaceCells;
    if (interfaceInfo != nullptr)
    {
        nInterfaceCells = interfaceInfo->GetNIFace();
    }
    nTotalSize = nEquations * globalTotalCells;

    //! initialize PETSc
    //! PetscFunctionBeginUser;
    PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);    
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    //! create vectors;
    VecCreateMPI(PETSC_COMM_WORLD, nTotalCells * nEquations, nTotalSize, &x);
    //! VecCreate(PETSC_COMM_WORLD, &x);
    PetscObjectSetName((PetscObject)x, "solution");
    //! VecSetSizes(x, nTotalCells * nEquations, nTotalSize);
    VecSetUp(x);
    VecDuplicate(x, &b);

    //! create matrix
    //! GMRESCSR
    PetscInt d_nnz[nTotalCells];
    PetscInt o_nnz[nTotalCells];
    for (int iCell = 0; iCell < nTotalCells; iCell++)
    {
        int count = 0;
        for (int j = AI[iCell]; j < AI[iCell + 1]; j++)
        {
            if(AJ[j] < nTotalCells)
            {
                count++;
            }
            
        }
        d_nnz[iCell] = count;
        o_nnz[iCell] = AK[iCell];
    }
    PH_Barrier();
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, nTotalCells * nEquations, (nTotalCells) * nEquations, nTotalSize, nTotalSize);
    MatSetType(A, MATMPIBAIJ);
    MatSetBlockSize(A, nEquations);
    MatMPIBAIJSetPreallocation(A, nEquations, PETSC_DEFAULT, d_nnz, PETSC_DEFAULT, o_nnz);
    MatSeqBAIJSetPreallocation(A, nEquations, PETSC_DEFAULT, d_nnz);
    PetscInt firstRow = globalTotalCells * nEquations;
    PetscInt lastOneMoreRow = (globalTotalCells + nTotalCells) * nEquations;
    MatGetOwnershipRange(A, &firstRow, &lastOneMoreRow);
    MatSetUp(A);
    MatZeroEntries(A);

    MatCreate(PETSC_COMM_WORLD, &APC);
    MatSetSizes(APC, nTotalCells * nEquations, (nTotalCells) * nEquations, nTotalSize, nTotalSize);
    MatSetType(APC, MATMPIBAIJ);
    MatSetBlockSize(APC, nEquations);

    //! GMRES Precond GMRESJac1st
    if( jacOrder == 2)
    {
        PetscInt nnz1st[nTotalCells];
        for (int iCell = 0; iCell < nTotalCells; iCell++)
        {
            int count = 0;
            for (int j = AI1st[iCell]; j < AI1st[iCell + 1]; j++)
            {
                if(AJ1st[j] < nTotalCells)
                {
                    count++;
                }
            }
            nnz1st[iCell] = count;
        }

        MatMPIBAIJSetPreallocation(APC, nEquations, PETSC_DEFAULT, nnz1st, 0, PETSC_NULL);
    }
    else 
    {
        MatMPIBAIJSetPreallocation(APC, nEquations, PETSC_DEFAULT, d_nnz, 0, PETSC_NULL);
    }
    MatGetOwnershipRange(APC, &firstRow, &lastOneMoreRow);
    MatSetUp(APC);
    MatZeroEntries(APC);

    //! improved approach to assemble matrix CSR
    PetscInt petsccol[nEquations];
    PetscInt petscrow;
    int jmax = AI[nTotalCells];
    int rcell = 0;
    int j = AI[rcell];

    //! consider the second order template for 2nd order Jacobian matrix for domain decomposition strategy
    // BUG
    if(interfaceInfo != nullptr)
    {
        vector<int> *neighborCells = grid->GMRESGetNeighborCells();
        int *localInterfacePhysicalCellIndex;
        localInterfacePhysicalCellIndex = interfaceInfo->GetLocalInterfacePhysicalCellIndex();
        rcell = nTotalCells;
        jmax  = AI[nTotalCells + nBoundFace];
        j = AI[rcell];
        do
        {
            int rowidx = interfaceInfo->MatchedGlobalNeighborCellIndex(rcell) * nEquations;
            if (rowidx > 0)    //! rcell is the ghost cell, rowidx > 0 means this cell locates on the interface boundary
            {
                int leftcell = interfaceInfo->MatchedLocalPhysicalCellIndex(rcell);
                for (int index = 0; index < neighborCells[leftcell].size(); index++)
                {
                    int qcell = neighborCells[leftcell][index];
                    if(qcell != rcell && qcell < nTotalCells)
                    {
                        printf("%d ", qcell + 1);
                    }
                }
            }
            rcell++;
        }while(rcell < nTotalCells + nBoundFace);
        printf("Finish writting\n");
    }
    PH_Barrier();
    TK_Exit::ExceptionExit(">");


    if(jacOrder == 2)
    {
        int jmax = AI1st[nTotalCells];
        int rcell = 0;
        int j = AI1st[rcell];
        do
        {
            int rowidx = (rcell + globalCellIndexShift) * nEquations;
            for (j = AI1st[rcell]; j < AI1st[rcell + 1]; j++)
            {
                int qcell = AJ1st[j];    //! get qcell (colume index of matrix)
                //! start to insert values: limited to dimension of nTotalSize
                if(qcell < nTotalCells)
                {
                    int colidx = (qcell + globalCellIndexShift) * nEquations;

                    for (int irow = 0; irow < nEquations; irow++)
                    {
                        petscrow = rowidx + irow;
                        for (int n = 0; n < nEquations; n++)
                        {
                            petsccol[n] = colidx + n;
                        }
                        MatSetValues(APC, 1, &petscrow, nEquations, petsccol, &dRdq1st[irow][j * nEquations], INSERT_VALUES);
                    }
                }
            }
            rcell++;    //! update rcell, go to the next row index of matrix
        } while (j < jmax);
    }
    else
    {
        int jmax = AI[nTotalCells];
        int rcell = 0;
        int j = AI[rcell];
        do
        {
            int rowidx = (rcell + globalCellIndexShift) * nEquations;
            int interfaceProcess = 0;
            for (j = AI[rcell]; j < AI[rcell + 1]; j++)
            {
                int qcell = AJ[j];    //! get qcell (colume index of matrix)
                //! start to insert values: limited to dimension of nTotalSize
                if(qcell < nTotalCells)
                {
                    int colidx = (qcell + globalCellIndexShift) * nEquations;

                    for (int irow = 0; irow < nEquations; irow++)
                    {
                        petscrow = rowidx + irow;
                        for (int n = 0; n < nEquations; n++)
                        {
                            petsccol[n] = colidx + n;
                        }
                        MatSetValues(APC, 1, &petscrow, nEquations, petsccol, &dRdq[irow][j * nEquations], INSERT_VALUES);
                    }
                }
            }
            rcell++;    //! update rcell, go to the next row index of matrix
        } while (j < jmax);
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(APC, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(APC, MAT_FINAL_ASSEMBLY);

    //! assemble vector
    VecSet(b, zero);
    for (PetscInt iCell = 0; iCell < nTotalCells; iCell++)
    {
        for (PetscInt iEquation = 0; iEquation < nEquations; iEquation++)
        {
            PetscInt isize = (iCell + globalCellIndexShift) * nEquations + iEquation;
            VecSetValues(b, 1, &isize, &res[iEquation][iCell], ADD_VALUES);
        }
    }
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    //! define solver
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);

    KSPSetOperators(ksp, A, APC);

    KSPSetPCSide(ksp, PC_RIGHT);
    KSPSetUp(ksp);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCBJACOBI);
    //! set sub-ksp and sub-pc
    KSP *subksp;
    PC  subpc;
    int localTotalSize = nTotalCells * nEquations;
    PCBJacobiSetLocalBlocks(pc, 1, &localTotalSize);
    int nLocal = 1;
    int firstLocal = grid->GetZoneID();
    PCBJacobiGetSubKSP(pc, &nLocal, &firstLocal, &subksp);
    KSPGetPC(subksp[0], &subpc);
    PCSetType(subpc, PCILU);
    PCFactorSetLevels(subpc, nLevels);

    KSPGMRESSetRestart(ksp, maxStepsRestart);
    KSPGMRESSetOrthogonalization(ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);
    KSPSetTolerances(ksp, rTolerance, 1.e-20, dropTolerance, maxSteps);

    //! solution
    KSPSolve(ksp, b, x);

    //! GMRESResidual
    RDouble **dq = dqProxy->GetField_UNS();

    // convert x to dqProxy
    for (PetscInt iCell = 0; iCell < nTotalCells; iCell++)
    {
        for (PetscInt iEquation = 0; iEquation < nEquations; iEquation++)
        {
            PetscInt isize = (iCell + globalCellIndexShift) * nEquations + iEquation;
            VecGetValues(x, 1, &isize, &dq[iEquation][iCell]);
        }
    }

    //! free work space
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    MatDestroy(&APC);
    KSPDestroy(&ksp);
    PetscFinalize();
    return;
}

void NSSolverUnstruct::GMRESSolver(Grid *gridIn, FieldProxy *dqProxy)
{
    //! GMRESCoupled
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    if (viscousType == ONE_EQU)
    {
        return;
    }

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int globalCellIndexShift = grid->GetGlobalCellIndexShift();
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();

    //! Memory allocating
    RDouble **dRdq = reinterpret_cast<RDouble**>(grid->GetDataPtr("dRdq"));
    RDouble **dDdP = reinterpret_cast<RDouble**>(grid->GetDataPtr("dDdP"));

    int nTotalCells = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nEquations = GetNumberOfEquations();

    RDouble **res = reinterpret_cast<RDouble**>(grid->GetDataPtr("res"));
    RDouble *vol = grid->GetCellVolume();
    RDouble *dt = reinterpret_cast<RDouble*> (grid->GetDataPtr("dt"));

    vector<int> AI = grid->GetJacobianAI4GMRES();    //! GMRESSparse GMRESCSR
    vector<int> AJ = grid->GetJacobianAJ4GMRES();    //! GMRESSparse GMRESCSR
    vector<int> AK = grid->GetJacobianAK4GMRES();    //! GMRES Parallel

    //! GMRESCSR
    for (int iCell = 0; iCell < nTotalCells; iCell++)
    {
        vector<int>::iterator result = find(AJ.begin() + AI[iCell], AJ.begin() + AI[iCell + 1], iCell);    //! find qcell
        int index = distance(AJ.begin(), result);
        index *= nEquations;
        for (int iEquation = 0; iEquation < nEquations; iEquation++)
        {
            dRdq[iEquation][index + iEquation] += 1.0 / dt[iCell];
        }
    }

    //! GMRESJac1st
    vector<int> AI1st = grid->GetJacobianAI1st4GMRES();
    vector<int> AJ1st = grid->GetJacobianAJ1st4GMRES();
    RDouble **dRdq1st = reinterpret_cast<RDouble**>(grid->GetDataPtr("dRdq1st"));
    int jacOrder = grid->GetJacobianOrder();
    if (jacOrder == 2)
    {
        for (int iCell = 0; iCell < nTotalCells; iCell++)
        {
            vector<int>::iterator result = find(AJ1st.begin() + AI1st[iCell], AJ1st.begin() + AI1st[iCell + 1], iCell); // find qcell
            int index = distance(AJ1st.begin(), result);
            index *= nEquations;
            for (int iEquation = 0; iEquation < nEquations; iEquation++)
            {
                dRdq1st[iEquation][index + iEquation] += 1.0 / dt[iCell];
            }
        }
    }

    //! entry to PETSc solve Ax = b equation for NS equations.
    Vec x, b;
    Mat A;
    Mat APC;                        //! GMRES Parallel
    KSP ksp;                        //! define solver
    PC pc;                          //! define matrix
    PetscInt nTotalSize;            //! dimension size of matrix and vector
    PetscInt nLevels = 3;           //! define number of levels used by preconditioner ILU
    PetscInt maxStepsRestart = 500; //! define maximum steps required by [restart]
    PetscInt maxSteps = 500;        //! define maximum steps
    PetscReal rTolerance = 1.e-1;   //! define tolerance of relative error // -1
    PetscReal dropTolerance = 1.e5; //! define tolerance of drop error    
    PetscMPIInt size;

    int globalTotalCells = (int)GlobalDataBase::GetDoubleParaFromDB("GlobalTotalCells");
    int nInterfaceCells = 0;
    if (interfaceInfo != nullptr)
    {
        nInterfaceCells = interfaceInfo->GetNIFace();
    }
    nTotalSize = nEquations * globalTotalCells;

    //! initialize PETSc
    PetscInitialize(PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    //! create vectors;
    VecCreateMPI(PETSC_COMM_WORLD, nTotalCells * nEquations, nTotalSize, &x);
    //! VecCreate(PETSC_COMM_WORLD, &x);
    PetscObjectSetName((PetscObject)x, "solution");
    //! VecSetSizes(x, nTotalCells * nEquations, nTotalSize);
    VecSetUp(x);
    VecDuplicate(x, &b);

    //! create matrix
    //! GMRESCSR
    PetscInt d_nnz[nTotalCells];
    PetscInt o_nnz[nTotalCells];
    for (int iCell = 0; iCell < nTotalCells; iCell++)
    {
        int count = 0;
        for (int j = AI[iCell]; j < AI[iCell + 1]; j++)
        {
            if (AJ[j] < nTotalCells)
            {
                count++;
            }
        }
        d_nnz[iCell] = count;
        o_nnz[iCell] = AK[iCell];
    }
    PH_Barrier();
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, nTotalCells * nEquations, (nTotalCells)*nEquations, nTotalSize, nTotalSize);
    MatSetType(A, MATMPIBAIJ);
    MatSetBlockSize(A, nEquations);
    MatMPIBAIJSetPreallocation(A, nEquations, PETSC_DEFAULT, d_nnz, 20, PETSC_NULL);
    PetscInt firstRow = globalTotalCells * nEquations;
    PetscInt lastOneMoreRow = (globalTotalCells + nTotalCells) * nEquations;
    MatGetOwnershipRange(A, &firstRow, &lastOneMoreRow);
    MatSetUp(A);
    MatZeroEntries(A);

    MatCreate(PETSC_COMM_WORLD, &APC);
    MatSetSizes(APC, nTotalCells * nEquations, (nTotalCells)*nEquations, nTotalSize, nTotalSize);
    MatSetType(APC, MATMPIBAIJ);
    MatSetBlockSize(APC, nEquations);

    //! GMRES Precond GMRESJac1st
    if (jacOrder == 2)
    {
        PetscInt nnz1st[nTotalCells];
        for (int iCell = 0; iCell < nTotalCells; iCell++)
        {
            int count = 0;
            for (int j = AI1st[iCell]; j < AI1st[iCell + 1]; j++)
            {
                if (AJ1st[j] < nTotalCells)
                {
                    count++;
                }
            }
            nnz1st[iCell] = count;
        }

        MatMPIBAIJSetPreallocation(APC, nEquations, PETSC_DEFAULT, nnz1st, 0, PETSC_NULL);
    }
    else
    {
        MatMPIBAIJSetPreallocation(APC, nEquations, PETSC_DEFAULT, d_nnz, 0, PETSC_NULL);
    }
    MatGetOwnershipRange(APC, &firstRow, &lastOneMoreRow);
    MatSetUp(APC);
    MatZeroEntries(APC);

    //! improved approach to assemble matrix CSR
    PetscInt petsccol[nEquations];
    PetscInt petscrow;
    int jmax = AI[nTotalCells];
    int rcell = 0;
    int j = AI[rcell];
    do
    {
        int rowidx = (rcell + globalCellIndexShift) * nEquations;
        int interfaceProcess = 0;
        for (j = AI[rcell]; j < AI[rcell + 1]; j++)
        {
            int qcell = AJ[j];    //! get qcell (colume index of matrix)
            //! start to insert values: limited to dimension of nTotalSize
            if (qcell < nTotalCells)
            {
                int colidx = (qcell + globalCellIndexShift) * nEquations;

                for (int irow = 0; irow < nEquations; irow++)
                {
                    petscrow = rowidx + irow;
                    for (int n = 0; n < nEquations; n++)
                    {
                        petsccol[n] = colidx + n;
                    }
                    MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdq[irow][j * nEquations], ADD_VALUES);
                }
            }
        }
        rcell++;    //! update rcell, go to the next row index of matrix
    } while (j < jmax);

    //! consider the second order template for 2nd order Jacobian matrix for domain decomposition strategy
    if (interfaceInfo != nullptr && jacOrder == 2)
    {
        vector<int> *neighborCells = grid->GMRESGetNeighborCells();
        int *localInterfacePhysicalCellIndex;
        // localInterfacePhysicalCellIndex = interfaceInfo->GetLocalInterfacePhysicalCellIndex();
        rcell = nTotalCells;
        jmax = AI[nTotalCells + nBoundFace];
        j = AI[rcell];
        do
        {
            int rowidx = interfaceInfo->MatchedGlobalNeighborCellIndex(rcell) * nEquations;
            if (rowidx >= 0)    //! rcell is the ghost cell, rowidx > 0 means this cell locates on the interface boundary
            {
                int leftcell = interfaceInfo->MatchedLocalPhysicalCellIndex(rcell); // return its corresponding physical cell index
                int colidx = (leftcell + globalCellIndexShift) * nEquations;
                vector<int>::iterator result = find(AJ.begin() + AI[rcell], AJ.begin() + AI[rcell + 1], leftcell); // find leftcell
                int index = distance(AJ.begin(), result);
                index *= nEquations;
                if (result == AJ.begin() + AI[rcell + 1])
                {

                    PrintToWindow(" rcell ", rcell, " cannot find qcell ", leftcell, " from ", "\n");
                    for (int i = AI[rcell]; i < AI[rcell + 1]; i++)
                    {
                        PrintToWindow("AJ ", AJ[i]);
                    }
                    TK_Exit::ExceptionExit("Sparse matrix index is wrong");

                }

                for (int irow = 0; irow < nEquations; irow++)
                {
                    petscrow = rowidx + irow;
                    for (int n = 0; n < nEquations; n++)
                    {
                        petsccol[n] = colidx + n;
                    }
                    MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdq[irow][index], ADD_VALUES);
                }

                //! obtain the list of the neighbor physical cell for the left cell
                for (int id = 0; id < neighborCells[leftcell].size(); id++)
                {
                    int neighborIndex = neighborCells[leftcell][id];
                    if (neighborIndex >= nTotalCells)
                    {
                        continue;
                    }

                    int colidx = (neighborIndex + globalCellIndexShift) * nEquations;
                    vector<int>::iterator result = find(AJ.begin() + AI[rcell], AJ.begin() + AI[rcell + 1], neighborIndex); // find neighborIndex
                    int index = distance(AJ.begin(), result);
                    index *= nEquations;
                    if (result == AJ.begin() + AI[rcell + 1])
                    {
                        PrintToWindow(" rcell ", rcell, " cannot find qcell ", neighborIndex, " from ", "\n");
                        for (int i = AI[rcell]; i < AI[rcell + 1]; i++)
                        {
                            PrintToWindow("AJ ", AJ[i]);
                        }
                        TK_Exit::ExceptionExit("Sparse matrix index is wrong");
                    }

                    for (int irow = 0; irow < nEquations; irow++)
                    {
                        petscrow = rowidx + irow;
                        for (int n = 0; n < nEquations; n++)
                        {
                            petsccol[n] = colidx + n;
                        }
                        MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdq[irow][index], ADD_VALUES);
                    }
                }

                //! obtain the R_physical/q_ghost caused by the neighbor cell at interface bc
                result = find(AJ.begin() + AI[leftcell], AJ.begin() + AI[leftcell + 1], rcell); // find rcell
                index = distance(AJ.begin(), result);
                index *= nEquations;
                if (result == AJ.begin() + AI[leftcell + 1])
                {
                    PrintToWindow(" rcell ", leftcell, " cannot find qcell ", rcell, " from ", "\n");
                    for (int i = AI[leftcell]; i < AI[leftcell + 1]; i++)
                    {
                        PrintToWindow("AJ ", AJ[i]);
                    }
                    TK_Exit::ExceptionExit("Sparse matrix index is wrong");
                }

                for (int irow = 0; irow < nEquations; irow++)
                {
                    petscrow = (leftcell + globalCellIndexShift) * nEquations + irow;
                    for (int n = 0; n < nEquations; n++)
                    {
                        petsccol[n] = rowidx + n;
                    }
                    MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdq[irow][index], ADD_VALUES);
                }

                //! obtain the list of the neighbor physical cell for the left cell
                for (int id = 0; id < neighborCells[leftcell].size(); id++)
                {
                    int neighborIndex = neighborCells[leftcell][id];
                    if (neighborIndex >= nTotalCells)
                    {
                        continue;
                    }

                    vector<int>::iterator result = find(AJ.begin() + AI[neighborIndex], AJ.begin() + AI[neighborIndex + 1], rcell); // find rcell
                    int index = distance(AJ.begin(), result);
                    index *= nEquations;
                    if (result == AJ.begin() + AI[neighborIndex + 1])
                    {
                        PrintToWindow(" rcell ", neighborIndex, " cannot find qcell ", rcell, " from ", "\n");
                        for (int i = AI[neighborIndex]; i < AI[neighborIndex + 1]; i++)
                        {
                            PrintToWindow("AJ ", AJ[i]);
                        }
                        TK_Exit::ExceptionExit("Sparse matrix index is wrong");
                    }

                    for (int irow = 0; irow < nEquations; irow++)
                    {
                        petscrow = (neighborIndex + globalCellIndexShift) * nEquations + irow;
                        for (int n = 0; n < nEquations; n++)
                        {
                            petsccol[n] = rowidx + n;
                        }
                        MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdq[irow][index], ADD_VALUES);
                    }
                }
            }
            rcell++;
        } while (rcell < nTotalCells + nBoundFace);
    }
    else if (interfaceInfo != nullptr && jacOrder != 2)
    {
        vector<int> *neighborCells = grid->GMRESGetNeighborCells();
        int *localInterfacePhysicalCellIndex;
        localInterfacePhysicalCellIndex = interfaceInfo->GetLocalInterfacePhysicalCellIndex();
        rcell = nTotalCells;
        jmax = AI[nTotalCells + nBoundFace];
        j = AI[rcell];

        do
        {
            int rowidx = interfaceInfo->MatchedGlobalNeighborCellIndex(rcell);
            if (rowidx >= 0) //! rcell is the ghost cell, rowidx > 0 means this cell locates on the interface boundary
            {
                int leftcell = interfaceInfo->MatchedLocalPhysicalCellIndex(rcell); // return its corresponding physical cell index
                int colidx = (leftcell + globalCellIndexShift) * nEquations;

                vector<int>::iterator result = find(AJ.begin() + AI[rcell], AJ.begin() + AI[rcell + 1], leftcell); // find leftcell
                int index = distance(AJ.begin(), result);
                index *= nEquations;
                if (result == AJ.begin() + AI[rcell + 1])
                {
                    PrintToWindow(" rcell ", rcell, " cannot find qcell ", leftcell, " from ", "\n");
                    for (int i = AI[rcell]; i < AI[rcell + 1]; i++)
                    {
                        PrintToWindow("AJ ", AJ[i]);
                    }
                    TK_Exit::ExceptionExit("Sparse matrix index is wrong");
                }

                rowidx *= nEquations;
                for (int irow = 0; irow < nEquations; irow++)
                {
                    petscrow = rowidx + irow;
                    for (int n = 0; n < nEquations; n++)
                    {
                        petsccol[n] = colidx + n;
                    }
                    MatSetValues(A, 1, &petscrow, nEquations, petsccol, &dRdq[irow][index], ADD_VALUES);
                }
            }
            rcell++;
        } while (rcell < nTotalCells + nBoundFace);
    }

    if (jacOrder == 2)
    {
        int jmax = AI1st[nTotalCells];
        int rcell = 0;
        int j = AI1st[rcell];
        do
        {
            int rowidx = (rcell + globalCellIndexShift) * nEquations;
            for (j = AI1st[rcell]; j < AI1st[rcell + 1]; j++)
            {
                int qcell = AJ1st[j];    //! get qcell (colume index of matrix)
                //! start to insert values: limited to dimension of nTotalSize
                if (qcell < nTotalCells)
                {
                    int colidx = (qcell + globalCellIndexShift) * nEquations;

                    for (int irow = 0; irow < nEquations; irow++)
                    {
                        petscrow = rowidx + irow;
                        for (int n = 0; n < nEquations; n++)
                        {
                            petsccol[n] = colidx + n;
                        }
                        MatSetValues(APC, 1, &petscrow, nEquations, petsccol, &dRdq1st[irow][j * nEquations], INSERT_VALUES);
                    }
                }
            }
            rcell++; // update rcell, go to the next row index of matrix
        } while (j < jmax);
    }
    else
    {
        int jmax = AI[nTotalCells];
        int rcell = 0;
        int j = AI[rcell];
        do
        {
            int rowidx = (rcell + globalCellIndexShift) * nEquations;
            int interfaceProcess = 0;
            for (j = AI[rcell]; j < AI[rcell + 1]; j++)
            {
                int qcell = AJ[j];    //! get qcell (colume index of matrix)
                //! start to insert values: limited to dimension of nTotalSize
                if (qcell < nTotalCells)
                {
                    int colidx = (qcell + globalCellIndexShift) * nEquations;

                    for (int irow = 0; irow < nEquations; irow++)
                    {
                        petscrow = rowidx + irow;
                        for (int n = 0; n < nEquations; n++)
                        {
                            petsccol[n] = colidx + n;
                        }
                        MatSetValues(APC, 1, &petscrow, nEquations, petsccol, &dRdq[irow][j * nEquations], INSERT_VALUES);
                    }
                }
            }
            rcell++;    //! update rcell, go to the next row index of matrix
        } while (j < jmax);
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(APC, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(APC, MAT_FINAL_ASSEMBLY);

    //! assemble vector
    VecSet(b, zero);
    for (PetscInt iCell = 0; iCell < nTotalCells; iCell++)
    {
        for (PetscInt iEquation = 0; iEquation < nEquations; iEquation++)
        {
            PetscInt isize = (iCell + globalCellIndexShift) * nEquations + iEquation;
            VecSetValues(b, 1, &isize, &res[iEquation][iCell], ADD_VALUES);
        }
    }
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    //! define solver
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPGMRES);

    KSPSetOperators(ksp, A, APC);

    KSPSetPCSide(ksp, PC_RIGHT);
    KSPSetUp(ksp);
    KSPGetPC(ksp, &pc);
    PCSetType(pc, PCBJACOBI);
    // set sub-ksp and sub-pc
    KSP* subksp;
    PC  subpc;
    int localTotalSize = nTotalCells * nEquations;
    PCBJacobiSetLocalBlocks(pc, 1, &localTotalSize);
    int nLocal = 1;
    int firstLocal = grid->GetZoneID();
    PCBJacobiGetSubKSP(pc, &nLocal, &firstLocal, &subksp);
    KSPGetPC(subksp[0], &subpc);
    PCSetType(subpc, PCILU);
    PCFactorSetLevels(subpc, nLevels);

    KSPGMRESSetRestart(ksp, maxStepsRestart);
    KSPGMRESSetOrthogonalization(ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);
    KSPSetTolerances(ksp, rTolerance, 1.e-20, dropTolerance, maxSteps);

    //! solution
    KSPSolve(ksp, b, x);

    //! GMRESResidual
    RDouble** dq = dqProxy->GetField_UNS();

    //! convert x to dqProxy
    for (PetscInt iCell = 0; iCell < nTotalCells; iCell++)
    {
        for (PetscInt iEquation = 0; iEquation < nEquations; iEquation++)
        {
            PetscInt isize = (iCell + globalCellIndexShift) * nEquations + iEquation;
            VecGetValues(x, 1, &isize, &dq[iEquation][iCell]);
        }
    }

    //! free work space
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    MatDestroy(&APC);
    KSPDestroy(&ksp);
    PetscFinalize();
    return;
}
#endif
}


