using namespace std;

#include "HODefine.h"

#include "ActionKey.h"
#include "Math_BasisFunction.h"
#include "MultiGridOperation.h"
#include "Pointer.h"
#include "Precision.h"
#include "TK_Exit.h"
#include "HOSolverUnstruct.h"
#include "HOGeometryStructure.h"

#include "HOSolverUnstructParam.h"

#include "HOBasisFunction.h"
#include "HOShapeFunction.h"
#include "HOStandardElement.h"
#include "HOGaussJacobiQuadrature.h"
#include "Geo_UnstructBC.h"
#include "HOBoundary.h"
#include "HOInviscidFlux.h"

#include "Residual.h"
#include "Post_Visual.h"
#include "IO_FileName.h"

namespace PHSPACE
{

HOSolverUnstruct::HOSolverUnstruct()
{
#ifdef HO_UNIT_TEST

    //TestOrthogoMatrix();

    //TestBasisFunction();

    //TestShapeFunction();

    //TestGaussLegendrePoints1D();
#endif
}

HOSolverUnstruct::~HOSolverUnstruct()
{

}

LIB_EXPORT void HOSolverUnstruct::InitControlParameters()
{
    if (HODebug) cout << "InitControlParameters: init. params, check params" << endl;

    controlParameters = new HOSolverUnstructParam();
    controlParameters->Init();
}

void HOSolverUnstruct::ReadParameter()
{
    if (HODebug) cout << "ReadParameter" << endl;
    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    string uns_scheme_name = parameters->GetSpaceSecheme();

    int uns_scheme = GetSchemeID(uns_scheme_name);
    GlobalDataBase::UpdateData("uns_scheme", &uns_scheme, PHINT, 1);
}

void HOSolverUnstruct::AllocateGlobalVar(Grid * gridIn)
{
    CFDSolver::AllocateGlobalVar(gridIn);
    
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int hLevel = grid->GetLevel();
    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);

    if (hLevel == 0) highOrderGrids.resize(parameters->GetNMGLevel());
 
    int nTotalCell = grid->GetNTotalCell();
    //int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;
    int nTotalNode = grid->GetNTotalNode();

    //int nm = 5;
    //! GlobalDataBase::GetData("nm", &nm, PHINT, 1);

    int nl = 5;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);
    int nchem = 0;
    GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);
    int neqn = nl + nchem;

    if (hLevel == 0) highOrderGrids.resize(parameters->GetNMGLevel());
    HighOrderGrid & highOrderGrid = highOrderGrids[hLevel];

    //! Cell
    //std::vector< HighOrderCell > & highOrderCells = highOrderGrids[hLevel].highOrderCells;
    highOrderGrid.CopyGeometryInfoToHighOrderCell(grid);

    //! Face
    //std::vector< HighOrderFace > & highOrderFaces = highOrderGrids[hLevel].highOrderFaces;
    highOrderGrid.CopyGeometryInfoToHighOrderFace(grid);

    //! Node
    //std::vector< HighOrderNode > & highOrderNodes = highOrderGrids[hLevel].highOrderNodes;
    highOrderGrid.CopyGeometryInfoToHighOrderNode(grid);

    HighOrderGrid::CellSolType & highOrderCellsSol = highOrderGrid.highOrderCellsSol;
    HighOrderGrid::FaceSolType & highOrderFacesSol = highOrderGrid.highOrderFacesSol;

    int pNLevels   = parameters->pMultiGrid;
    int dgSolOrder = parameters->dgSolOrder;
    const int isUnsteady  = parameters->GetIsUnsteady();
    const int viscousType = parameters->GetViscousType();

    if (hLevel == 0)
    {
        if (pNLevels <= 0)
        {
            pNLevels = 1;
        }

        if ((pNLevels != (dgSolOrder + 1)) && (pNLevels > 1))
        {
            cout << " The P_multigrid of DG_parm is wrong, please checking  " << endl;
        }

        //! totalFaceGaussPoint
        highOrderGrids[hLevel].totalFaceGaussPoint.resize(pNLevels);
        highOrderCellsSol.resize(pNLevels);
        highOrderFacesSol.resize(pNLevels);

        for (int img = 0; img < pNLevels; ++ img)
        {
            highOrderGrid.InitHighOrderCellSol(grid, dgSolOrder, img, isUnsteady, viscousType);
            highOrderGrid.InitHighOrderFaceSol(grid, img, viscousType);
        }
    }
    else
    {
        //! totalFaceGaussPoint
        highOrderGrids[hLevel].totalFaceGaussPoint.resize(1);
        highOrderCellsSol.resize(1);
        highOrderFacesSol.resize(1);
        highOrderGrid.InitHighOrderCellSol(grid, 0, 0, isUnsteady, viscousType);
        highOrderGrid.InitHighOrderFaceSol(grid, 0, viscousType);
    }    

    RDouble *dt = new RDouble [ nTotal ];
    grid->UpdateDataPtr("dt" , dt);

    int *dofsOfCells = new int[nTotalCell];
    highOrderGrid.GetNumberOfDOFs(nTotalCell, dofsOfCells);

    RDouble ***q   = NewPointer3<RDouble>(nTotalCell, neqn, dofsOfCells);
    RDouble ***res = NewPointer3<RDouble>(nTotalCell, neqn, dofsOfCells);

    delete [] dofsOfCells;
    
    grid->UpdateDataPtr("q"  , q );
    grid->UpdateDataPtr("res", res);

    if (grid->GetLevel() == 0)
    {
        int *nodeBCType = new int[nTotalNode];
        PHSPACE::SetField(nodeBCType, 0, nTotalNode);
        grid->UpdateDataPtr("nodeBCType", nodeBCType);
    }
}

void HOSolverUnstruct::DeAllocateGlobalVar(Grid * gridIn)
{
    if (HODebug) cout << "DeAllocateGlobalVar" << endl;
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble *dt = reinterpret_cast<RDouble *> (gridIn->GetDataPtr("dt"));
    delete [] dt;

    RDouble ***q   = reinterpret_cast< RDouble *** > (gridIn->GetDataPtr("q" ));
    RDouble ***res = reinterpret_cast< RDouble *** > (gridIn->GetDataPtr("res"));

    DelPointer3(q );
    DelPointer3(res);

    if (grid->GetLevel() == 0)
    {
        int *nodeBCType = reinterpret_cast <int *> (grid->GetDataPtr("nodeBCType"));
        delete [] nodeBCType;
    }
}

void HOSolverUnstruct::RegisterCFDSolverInterfaceField()
{
    if (HODebug) cout << "RegisterCFDSolverInterfaceField" << endl;
}

void HOSolverUnstruct::InitMemory()
{
    InitControlParameters();
    ReadParameter();
    AllocateGlobalVariables();
}

void HOSolverUnstruct::ReleaseMemory()
{
    DeAllocateGlobalVariables();
}

void HOSolverUnstruct::AllocateGlobalVariables()
{
    Grid *grid = GetGrid();
    while (grid)
    {
        AllocateGlobalVar(grid);
        grid = grid->GetCoarseGrid();
    }
}

bool HOSolverUnstruct::JudgeIfReadAverage()
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

bool HOSolverUnstruct::JudgeIfRestart()
{
    string restartNSFile = ".\results\flow.dat";
    GlobalDataBase::GetData("restartNSFile", &restartNSFile, PHSTRING, 1);

    if (PHMPI::IsParallelRun())
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

    return restart_flag;
}

void HOSolverUnstruct::InitFlow()
{
    int nstart;

    Param_CFDSolver *parameters = GetControlParameters();

    //! This method would be wrong when turbulent computational from laminar flow field.
    //! But this situation breaks the law, so it is not considered.
    bool readFlowFile = false;
    if (JudgeIfRestart())
    {
        readFlowFile = true;

        nstart = 1;
        GlobalDataBase::UpdateData("nstart", &nstart, PHINT, 1);
    }

    if (!readFlowFile)
    {
        int outnstep = 0;
        GlobalDataBase::UpdateData("outnstep",&outnstep, PHINT, 1);

        nstart = 0;
        GlobalDataBase::UpdateData("nstart", &nstart, PHINT, 1);

        InitFlowAsRestart();
    }
    else
    {
        //! To read data in restart file for continual simulation.
        ActionKey *actkeyReadRestartData = new ActionKey();
        FillActionKey(actkeyReadRestartData, READ_RESTART, 0);
        InitFlowAsReadingRestart(actkeyReadRestartData);
        delete actkeyReadRestartData;
    }

    RDouble globalMinTimeStep = -LARGE;
    GlobalDataBase::UpdateData("globalMinTimeStep",&globalMinTimeStep, PHDOUBLE, 1);

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady)
    {
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
        int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");
        bool isReadAverageFlow = ifStaticsFlowField && (outnstep >= startStatisticStep);
        if (isReadAverageFlow && readFlowFile && JudgeIfReadAverage())
        {
            //! To dump restart data for continual simulation.
            ActionKey *actkeyReadAverageFlow = new ActionKey();
            FillActionKey(actkeyReadAverageFlow, READ_AVERAGE_FLOW, 0);
            //ReadAverageFlow(actkeyReadAverageFlow);
            delete actkeyReadAverageFlow;
        }
        else
        {
            int nStatisticalStep = 0;
            GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
        }
    }

    //! Compute gama, temperature, viscousity, and set boundary condition once.
    InitDependentVariables();
}

void HOSolverUnstruct::InitFlowAsRestart()
{
    UnstructGrid * grid = UnstructGridCast(GetGrid(0));
    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);

    int outnstep = 0;
    GlobalDataBase::UpdateData("outnstep",&outnstep, PHINT, 1);

    RDouble physicalTime = 0.0;
    GlobalDataBase::UpdateData("physicalTime",&physicalTime, PHFLOAT, 1);

    int nTotalCell = grid->GetNTotalCell();

    int nl = 5;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);
    int nchem = 0;
    GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);
    int neqn = nl + nchem;

    RDouble *prim_inf = reinterpret_cast<RDouble *> (GlobalDataBase::GetDataPtr("prim_inf" ));
    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");

    vector < RDouble > conserValue(5,0.0);

    conserValue[0] = prim_inf[0];
    conserValue[1] = prim_inf[0] * prim_inf[1];
    conserValue[2] = prim_inf[0] * prim_inf[2];
    conserValue[3] = prim_inf[0] * prim_inf[3];
    conserValue[4] = prim_inf[4]/(refGama - 1) + 0.5 * prim_inf[0] * (pow(prim_inf[1],2) + pow(prim_inf[2],2) + pow(prim_inf[3],2));

    //! Init dof of cell and initialize value of gausspoint
    highOrderGrids[0].InitDofOfCell(grid, conserValue, neqn);

    //! Initialize value of gausspoint for face
    highOrderGrids[0].InitValueOfFace(grid, conserValue, neqn);

    int hLevel = grid->GetLevel();

    RDouble *** q = reinterpret_cast< RDouble *** > (grid->GetDataPtr("q"));

    int ndof = 0;
    for (int ieqn = 0; ieqn < neqn; ++ ieqn)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            ndof = highOrderGrids[hLevel].highOrderCellsSol[0][iCell].nDOFsSol;

            for (int idof = 0; idof < ndof; ++ idof)
            {
                q[iCell][ieqn][idof] = highOrderGrids[hLevel].highOrderCellsSol[0][iCell].dofQ[ieqn][idof];
            }
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();

    if (!isUnsteady) return;
    //! For unsteady problem, the initial condition must be high order
    cout << " the initial condition must be high order for unsteady problem------> InitFlowAsRestart()" << endl;

    highOrderGrids[0].InitDofOfCellUnsteady(grid, conserValue, neqn);
}

void HOSolverUnstruct::InitDependentVariables()
{
    Grid *fgrid = GetGrid(0);
    Boundary(fgrid);
    
    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    int viscousType = parameters->GetViscousType();
    if (viscousType > INVISCID)
    {
        CompViscousCoef(fgrid);
        InitViscousCoefTurb(fgrid);
    }
}

void HOSolverUnstruct::Boundary(Grid *gridIn)
{
    if (HODebug) cout << "Boundary: set the ghost cell q£¬DG needs to implement the special BC accroding to the BC type." << endl;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *left_cell_of_face  = grid->GetLeftCellOfFace();
    int *right_cell_of_face = grid->GetRightCellOfFace();

    UnstructBCSet **bcr = grid->GetBCRecord();

    int nBoundFace = grid->GetNBoundFace();
    
    HOSolverUnstructParam *parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    int hLevel = grid->GetLevel(); 
    int pLevel = parameters->pMultiGrid;
    
    if (pLevel == 1)
    {
        pLevel = 0;
    }

    if (hLevel > 1 && pLevel > 1)
    {
        cout << "the multigrid para is wrong for unstructure high order in HOSolverUnstruct::Boundary(Grid *gridIn)" << endl;
    }

    int nm = 5;
    GlobalDataBase::GetData("nm", &nm, PHINT, 1);

    int nl;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);

    int nchem;
    GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);
    int nEquation = nl + nchem;

    double refGama = parameters->GetRefGama();

    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble wallTemperature = parameters->GetWallTemperature();

    double refDimensionalTemperature = parameters->GetRefDimensionalTemperature();

    double tw = wallTemperature / refDimensionalTemperature;

    int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

    RDouble *prims, *primt;
    prims = new RDouble[nl+nchem];
    primt = new RDouble[nl+nchem];

    RDouble *primInflow = reinterpret_cast< RDouble * >(GlobalDataBase::GetDataPtr("prim_inf"));
    bool massFlow = false;
    RDouble massAlpha = 1.0;
    double massIn  = 1.0;
    double massOut = 1.0;

    if (massFlow)
    {
        //ComputeMassFlow(grid);
        GlobalDataBase::GetData("massIn", &massIn, PHDOUBLE, 1);
        GlobalDataBase::GetData("massOut", &massOut, PHDOUBLE, 1);
        massAlpha = massIn / massOut;
        cout << " massAlpha = " << massAlpha << "\n";
    }

    using namespace IDX;
    using namespace PHENGLEI;

    int le, re, bcType, numberGaussPoint, rGausspoint;
    RDouble xNormal, yNormal, zNormal, vgn;
    RDouble xtn = 0.0, ytn = 0.0, ztn = 0.0;
    //vector< vector< RDouble > > conservationQ;

    bool isViscous = parameters->IsViscous();

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = left_cell_of_face [ iFace ];
        re = right_cell_of_face[ iFace ];
        bcType = bcr[iFace]->GetKey();
        HighOrderFaceSol & highOrderFaceSol = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace];
        numberGaussPoint = highOrderFaceSol.integNumberOfGaussPoint;

        if (IsInterface(bcType)) //! Interface
        {
            continue;
        }
        else if (bcType == SYMMETRY)
        {
            for (int iGaussPoint = 0; iGaussPoint < numberGaussPoint; ++iGaussPoint)
            {
                xNormal = highOrderFaceSol.metricNormalsFace[iGaussPoint][0];
                yNormal = highOrderFaceSol.metricNormalsFace[iGaussPoint][1];
                zNormal = highOrderFaceSol.metricNormalsFace[iGaussPoint][2];

                vector< vector< RDouble > > & conservationQ = highOrderFaceSol.q;

                vgn = 0.0;

                for (int m = 0; m < nEquation; ++ m)
                {
                    prims[m] = conservationQ[m][iGaussPoint];
                }

                prims[1] = prims[1] / prims[0];
                prims[2] = prims[2] / prims[0];
                prims[3] = prims[3] / prims[0];
                prims[4] = HO_GAMAM1 * (prims[4] - 0.5 * prims[0] * (pow(prims[1],2) + pow(prims[2],2) + pow(prims[3],2)));

                SymmetryBC(prims, primt, xNormal, yNormal, zNormal, vgn, nl, nchem);

                int rGaussPoint = iGaussPoint + numberGaussPoint;
                conservationQ[0][rGaussPoint] = primt[0];
                conservationQ[1][rGaussPoint] = primt[0] * primt[1];
                conservationQ[2][rGaussPoint] = primt[0] * primt[2];
                conservationQ[3][rGaussPoint] = primt[0] * primt[3];
                conservationQ[4][rGaussPoint] = primt[4]/HO_GAMAM1 + 0.5 * primt[0] * (pow(primt[1],2) + pow(primt[2],2) + pow(primt[3],2));
            }

            continue;
        }
        else if (IsWall(bcType))
        {
            if (!isViscous)
            {
                for (int iGausspoint = 0; iGausspoint < numberGaussPoint; ++iGausspoint)
                {
                    xNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][0];
                    yNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][1];
                    zNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][2];
                    vector< vector< RDouble > > & conservationQ = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].q;

                    vgn = 0.0;

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prims[m] = conservationQ[m][iGausspoint];
                    }
                    prims[1] = prims[1] / prims[0];
                    prims[2] = prims[2] / prims[0];
                    prims[3] = prims[3] / prims[0];
                    prims[4] = HO_GAMAM1 * (prims[4] - 0.5 * prims[0] * (pow(prims[1],2) + pow(prims[2],2) + pow(prims[3],2)));

                    SymmetryBC(prims, primt, xNormal, yNormal, zNormal, vgn, nl, nchem);

                    rGausspoint = iGausspoint + numberGaussPoint;
                    conservationQ[0][rGausspoint] = primt[0];
                    conservationQ[1][rGausspoint] = primt[0] * primt[1];
                    conservationQ[2][rGausspoint] = primt[0] * primt[2];
                    conservationQ[3][rGausspoint] = primt[0] * primt[3];
                    conservationQ[4][rGausspoint] = primt[4]/HO_GAMAM1 + 0.5 * primt[0] * (pow(primt[1],2) + pow(primt[2],2) + pow(primt[3],2));
                }
            }
            else if (wallTemperature <= 0.0)
            {
                //! Viscous adiabatic wall
                for (int iGausspoint = 0; iGausspoint < numberGaussPoint; ++iGausspoint)
                {
                    xNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][0];
                    yNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][1];
                    zNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][2];
                    vector< vector< RDouble > > & conservationQ = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].q;                    

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prims[m] = conservationQ[m][iGausspoint];
                    }
                    prims[1] = prims[1] / prims[0];
                    prims[2] = prims[2] / prims[0];
                    prims[3] = prims[3] / prims[0];
                    prims[4] = HO_GAMAM1 * (prims[4] - 0.5 * prims[0] * (pow(prims[1],2) + pow(prims[2],2) + pow(prims[3],2)));

                    ViscousAdiabaticWallBC(prims, primt, xNormal, yNormal, zNormal, xtn, ytn, ztn, refMachNumber, nl, nchem);
                    
                    rGausspoint = iGausspoint + numberGaussPoint;
                    conservationQ[0][rGausspoint]  =   primt[0];
                    conservationQ[1][rGausspoint]  =   primt[0] * primt[1];
                    conservationQ[2][rGausspoint]  =   primt[0] * primt[2];
                    conservationQ[3][rGausspoint]  =   primt[0] * primt[3];
                    conservationQ[4][rGausspoint]  =   primt[4]/HO_GAMAM1 + 0.5 * primt[0] * (pow(primt[1],2) + pow(primt[2],2) + pow(primt[3],2));
                }
            }
            else
            {
                //! Viscous iso-thermal wall
                for (int iGausspoint = 0; iGausspoint < numberGaussPoint; ++iGausspoint)
                {
                    xNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][0];
                    yNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][1];
                    zNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][2];
                    vector< vector< RDouble > > & conservationQ = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].q;

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prims[m] = conservationQ[m][iGausspoint];
                    }
                    prims[1] = prims[1] / prims[0];
                    prims[2] = prims[2] / prims[0];
                    prims[3] = prims[3] / prims[0];
                    prims[4] = HO_GAMAM1 * (prims[4] - 0.5 * prims[0] * (pow(prims[1],2) + pow(prims[2],2) + pow(prims[3],2)));

                    ViscousIsotropicWallBC(prims, primt, xNormal, yNormal, zNormal, xtn, ytn, ztn, refMachNumber, tw, nl, nchem);
                    
                    rGausspoint = iGausspoint + numberGaussPoint;
                    conservationQ[0][rGausspoint]  =   primt[0];
                    conservationQ[1][rGausspoint]  =   primt[0] * primt[1];
                    conservationQ[2][rGausspoint]  =   primt[0] * primt[2];
                    conservationQ[3][rGausspoint]  =   primt[0] * primt[3];
                    conservationQ[4][rGausspoint]  =   primt[4]/HO_GAMAM1 + 0.5 * primt[0] * (pow(primt[1],2) + pow(primt[2],2) + pow(primt[3],2));
                }
            }
            
            continue;
        }
        else if (bcType == FARFIELD)
        {
            if (ifLowSpeedPrecon == 0)
            {
                for (int iGausspoint = 0; iGausspoint < numberGaussPoint; ++ iGausspoint)
                {
                    xNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][0];
                    yNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][1];
                    zNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][2];
                    vector< vector< RDouble > > & conservationQ = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].q;
                    vgn = 0.0;

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        prims[m] = conservationQ[m][iGausspoint];
                    }

                    prims[1] = prims[1] / prims[0];
                    prims[2] = prims[2] / prims[0];
                    prims[3] = prims[3] / prims[0];

                    prims[4] = HO_GAMAM1 * (prims[4] - 0.5 * prims[0] * (pow(prims[1],2) + pow(prims[2],2) + pow(prims[3],2)));

                    //FarfieldBC(prims, primInflow, primt, refGama, HO_GAMA, xNormal, yNormal, zNormal, vgn, nl, nchem);
                    FarFieldRiemannInvariants(prims, primInflow, primt, refGama, HO_GAMA, xNormal, yNormal, zNormal, vgn, nl, nchem);

                    rGausspoint = iGausspoint + numberGaussPoint;
                    conservationQ[0][rGausspoint] = primt[0];
                    conservationQ[1][rGausspoint] = primt[0] * primt[1];
                    conservationQ[2][rGausspoint] = primt[0] * primt[2];
                    conservationQ[3][rGausspoint] = primt[0] * primt[3];
                    conservationQ[4][rGausspoint] = primt[4] / HO_GAMAM1 + 0.5 * primt[0] * (pow(primt[1],2) + pow(primt[2],2) + pow(primt[3],2));
                }
            }
            else
            {
            }
        }
        else if (bcType == INFLOW)
        {
            for (int iGausspoint = 0; iGausspoint < numberGaussPoint; ++ iGausspoint)
            {
                xNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][0];
                yNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][1];
                zNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][2];
                vector< vector< RDouble > > & conservationQ = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].q;
                vgn = 0.0;

                InflowBC(primInflow, primt, nl, nchem);

                rGausspoint = iGausspoint + numberGaussPoint;
                conservationQ[0][rGausspoint] = primt[0];
                conservationQ[1][rGausspoint] = primt[0] * primt[1];
                conservationQ[2][rGausspoint] = primt[0] * primt[2];
                conservationQ[3][rGausspoint] = primt[0] * primt[3];
                conservationQ[4][rGausspoint] = primt[4] / HO_GAMAM1 + 0.5 * primt[0] * (pow(primt[1],2) + pow(primt[2],2) + pow(primt[3],2));
            } 

            continue;
        }
        else if (bcType == OUTFLOW)
        { 
            for (int iGausspoint = 0; iGausspoint < numberGaussPoint; ++ iGausspoint)
            {
                xNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][0];
                yNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][1];
                zNormal = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].metricNormalsFace[iGausspoint][2];
                vector< vector< RDouble > > & conservationQ = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].q;
                vgn = 0.0;

                for (int m = 0; m < nEquation; ++ m)
                {
                    prims[m] = conservationQ[m][iGausspoint];
                }
                prims[1] = prims[1] / prims[0];
                prims[2] = prims[2] / prims[0];
                prims[3] = prims[3] / prims[0];
                prims[4] = HO_GAMAM1 * (prims[4] - 0.5 * prims[0] * (pow(prims[1],2) + pow(prims[2],2) + pow(prims[3],2)));

                OutflowBC(prims, primt, nl, nchem);

                rGausspoint = iGausspoint + numberGaussPoint;
                conservationQ[0][rGausspoint]  =   primt[0];
                conservationQ[1][rGausspoint]  =   primt[0] * primt[1];
                conservationQ[2][rGausspoint]  =   primt[0] * primt[2];
                conservationQ[3][rGausspoint]  =   primt[0] * primt[3];
                conservationQ[4][rGausspoint]  =   primt[4]/HO_GAMAM1 + 0.5 * primt[0] * (pow(primt[1],2) + pow(primt[2],2) + pow(primt[3],2));
            }

            continue;
        }
        else
        {
            TK_Exit::ExceptionExit("Error: this boundary type does not exist!\n");
        }
    }
    delete [] prims;
    delete [] primt;
}

void HOSolverUnstruct::CompViscousCoef(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int      nTotalCell = grid->GetNTotalCell();
    int      nTotalFace = grid->GetNTotalFace();
    int      nBoundFace = grid->GetNBoundFace();
    int      nTotal = nTotalCell + nBoundFace;

    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    int hLevel = grid->GetLevel(); 
    int pLevel     = parameters->pMultiGrid;
    if (pLevel == 1)
    {
        pLevel = 0;
    }

    RDouble tsuth;
    GlobalDataBase::GetData("tsuth", &tsuth, PHDOUBLE, 1);

    RDouble visl_min;
    GlobalDataBase::GetData("visl_min", &visl_min, PHDOUBLE, 1);

    RDouble R = parameters->GetR();

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        //! Compute gauss point temperature of volume
        int numberGaussPoint = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell].integNumberOfGaussPoint;
        for (int iGaussPoint = 0; iGaussPoint < numberGaussPoint; ++ iGaussPoint)
        {
            RDouble conserValue[5];
            //! Get conservation variable
            conserValue[0] = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell].q[0][iGaussPoint];
            conserValue[1] = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell].q[1][iGaussPoint];
            conserValue[2] = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell].q[2][iGaussPoint];
            conserValue[3] = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell].q[3][iGaussPoint];
            conserValue[4] = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell].q[4][iGaussPoint];

            RDouble pre = ComputePressure(conserValue);
            RDouble T = pre / conserValue[0] / R;
            highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell].visl[iGaussPoint] = pow(T,1.5) * (1.0 + tsuth) / (T + tsuth);
        }
    }

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        //! Compute gauss point temperature of volume
        int numberGaussPoint = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].integNumberOfGaussPoint;
        for (int iGaussPoint = 0; iGaussPoint < numberGaussPoint; ++ iGaussPoint)
        {
            RDouble conserValue[5];
            //! Get conservation variable
            conserValue[0] = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].q[0][iGaussPoint];
            conserValue[1] = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].q[1][iGaussPoint];
            conserValue[2] = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].q[2][iGaussPoint];
            conserValue[3] = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].q[3][iGaussPoint];
            conserValue[4] = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].q[4][iGaussPoint];

            RDouble pre = ComputePressure(conserValue);
            RDouble T = pre / conserValue[0] / R;

            highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].visl[iGaussPoint] = pow(T,1.5) * (1.0 + tsuth) / (T + tsuth);
        }
    }
}

void HOSolverUnstruct::ComputeGaussPointQ(Grid * gridIn)
{
    UnstructGrid * grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int * left_cell_of_face  = grid->GetLeftCellOfFace();
    int * right_cell_of_face = grid->GetRightCellOfFace();

    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    int hLevel = grid->GetLevel();
    int pLevel = parameters->pMultiGrid;

    if (pLevel == 1)
    {
        pLevel = 0;
    }

    RDouble tsuth;
    GlobalDataBase::GetData("tsuth", &tsuth, PHDOUBLE, 1);

    RDouble visl_min;
    GlobalDataBase::GetData("visl_min", &visl_min, PHDOUBLE, 1);

    int nl = 5;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);
    int nchem;
    GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);
    int nEquation = nl + nchem;

    int numberGaussPoint = 0, rGaussPoint = 0, ndof = 0;
    int leftCell = 0, rightCell = 0;

    RDouble qTmp = 0.0, qLTmp, qRTmp;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        //! Compute gauss point temperature of volume
        HighOrderCellSol & uHOCellSol = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell];
        numberGaussPoint = uHOCellSol.integNumberOfGaussPoint;
        ndof = uHOCellSol.nDOFsSol;        
        for (int iGaussPoint = 0; iGaussPoint < numberGaussPoint; ++ iGaussPoint)
        {
            //! Compute conservation variable
            for (int ieqn = 0; ieqn < nEquation; ++ ieqn)
            {
                //uHOCellSol.q[ieqn][iGaussPoint] = 0.0;
                qTmp = 0.0;
                for (int idof =0; idof < ndof; ++ idof)
                {
                    qTmp += uHOCellSol.basisFuncCellIntegration[iGaussPoint][idof] * uHOCellSol.dofQ[ieqn][idof];
                }

                uHOCellSol.q[ieqn][iGaussPoint] = qTmp;
            }
        }
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        HighOrderFaceSol & uHOFaceSol = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace];
        //! Compute gauss point temperature of face
        numberGaussPoint = uHOFaceSol.integNumberOfGaussPoint;
        leftCell = left_cell_of_face[iFace];
        HighOrderCellSol & uHOLeftCellSol = highOrderGrids[hLevel].highOrderCellsSol[pLevel][leftCell];

        for (int iGaussPoint = 0; iGaussPoint < numberGaussPoint; ++ iGaussPoint)
        {
            //! Compute conservation variable
            rGaussPoint = iGaussPoint + numberGaussPoint;
            for (int ieqn = 0; ieqn < nEquation; ++ ieqn)
            {
                qLTmp = 0.0;
                for (int idof =0; idof < ndof; ++ idof)
                {
                    qLTmp += uHOFaceSol.leftBasisFuncFaceIntegration[iGaussPoint][idof] * uHOLeftCellSol.dofQ[ieqn][idof];
                }
                uHOFaceSol.q[ieqn][iGaussPoint] = qLTmp;
                uHOFaceSol.q[ieqn][rGaussPoint] = 0.0;
            }
        }
    }
    
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        HighOrderFaceSol & uHOFaceSol = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace];
        //! Compute gauss point temperature of face
        numberGaussPoint = uHOFaceSol.integNumberOfGaussPoint;
        leftCell  = left_cell_of_face[iFace];
        rightCell = right_cell_of_face[iFace];
        HighOrderCellSol & uHOLeftCellSol = highOrderGrids[hLevel].highOrderCellsSol[pLevel][leftCell];
        HighOrderCellSol & uHORigtCellSol = highOrderGrids[hLevel].highOrderCellsSol[pLevel][rightCell];
        for (int iGaussPoint = 0; iGaussPoint < numberGaussPoint; ++ iGaussPoint)
        {
            //! Compute conservation variable
            rGaussPoint = iGaussPoint + numberGaussPoint;
            for (int ieqn = 0; ieqn < nEquation; ++ ieqn)
            {
                qLTmp = 0.0;
                qRTmp = 0.0;
                for (int idof =0; idof < ndof; ++ idof)
                {
                    qLTmp += uHOFaceSol.leftBasisFuncFaceIntegration[iGaussPoint][idof] * uHOLeftCellSol.dofQ[ieqn][idof];
                    qRTmp += uHOFaceSol.rightBasisFuncFaceIntegration[iGaussPoint][idof] * uHORigtCellSol.dofQ[ieqn][idof];
                }
                uHOFaceSol.q[ieqn][iGaussPoint] = qLTmp;
                uHOFaceSol.q[ieqn][rGaussPoint] = qRTmp;
            }
        }
    }
}

void HOSolverUnstruct::ComputeQAverage(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int      nTotalCell = grid->GetNTotalCell();

    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    int hLevel = grid->GetLevel(); 
    int pLevel     = parameters->pMultiGrid;
    bool isViscous = parameters->IsViscous();
    if (pLevel == 1)
    {
        pLevel = 0;
    }

    RDouble tsuth;
    GlobalDataBase::GetData("tsuth", &tsuth, PHDOUBLE, 1);

    RDouble visl_min;
    GlobalDataBase::GetData("visl_min", &visl_min, PHDOUBLE, 1);

    int nl = 5;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);

    int nchem;
    GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);
    int nEquation = nl + nchem;
    int numberGaussPoint = 0; 

    const RDouble *weightCoeff;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        //! Compute QAverage of volume
        HighOrderCellSol & uHOCellSol = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell];
        HighOrderCell    & uHOCell    = highOrderGrids[hLevel].highOrderCells[iCell];
        numberGaussPoint = uHOCellSol.integNumberOfGaussPoint;
        for (int ieqn = 0; ieqn < nEquation; ++ ieqn)
        {
            uHOCellSol.qAverage[ieqn] = 0.0;
        }

        weightCoeff = highOrderGrids[hLevel].standardCellElements[uHOCell.standardElementIndex].GetWtIntegration();

        for (int iGaussPoint = 0; iGaussPoint < numberGaussPoint; ++ iGaussPoint)
        {
            //! Compute conservation variable
            for (int ieqn = 0; ieqn < nEquation; ++ ieqn)
            {
                uHOCellSol.qAverage[ieqn] += weightCoeff[iGaussPoint] * uHOCellSol.q[ieqn][iGaussPoint] * uHOCellSol.JacDetCellIntegration[iGaussPoint];
            }
        }
        for (int ieqn = 0; ieqn < nEquation; ++ ieqn)
        {
            uHOCellSol.qAverage[ieqn] /= uHOCell.volume;
        }
    }

    if (isViscous == 1)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            //! Compute QAverage of volume
            HighOrderCellSol & uHOCellSol = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell];
            HighOrderCell    & uHOCell    = highOrderGrids[hLevel].highOrderCells[iCell];
            numberGaussPoint = uHOCellSol.integNumberOfGaussPoint;
            
            uHOCellSol.vislAverage = 0.0;
            weightCoeff   =   highOrderGrids[hLevel].standardCellElements[uHOCell.standardElementIndex].GetWtIntegration();
            for (int iGaussPoint = 0; iGaussPoint < numberGaussPoint; ++ iGaussPoint)
            {
                //! Compute conservation variable.
                uHOCellSol.vislAverage += weightCoeff[iGaussPoint] * uHOCellSol.visl[iGaussPoint] * uHOCellSol.JacDetCellIntegration[iGaussPoint];
            }
            uHOCellSol.vislAverage /= uHOCell.volume;
        }
    }
    else if (isViscous > 1)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            //! Compute QAverage of volume
            HighOrderCellSol & uHOCellSol = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell];
            HighOrderCell    & uHOCell    = highOrderGrids[hLevel].highOrderCells[iCell];
            numberGaussPoint = uHOCellSol.integNumberOfGaussPoint;

            uHOCellSol.vistAverage = 0.0;
            weightCoeff = highOrderGrids[hLevel].standardCellElements[uHOCell.standardElementIndex].GetWtIntegration();
            for (int iGaussPoint = 0; iGaussPoint < numberGaussPoint; ++ iGaussPoint)
            {
                //! Compute conservation variable
                uHOCellSol.vistAverage += weightCoeff[iGaussPoint] * uHOCellSol.vist[iGaussPoint] * uHOCellSol.JacDetCellIntegration[iGaussPoint];
            }
            uHOCellSol.vistAverage /= uHOCell.volume;
        }
    }
}

void HOSolverUnstruct::InitViscousCoefTurb(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int      nTotalCell = grid->GetNTotalCell();
    int      nTotalFace = grid->GetNTotalFace();
    int      nBoundFace = grid->GetNBoundFace();
    int      nTotal = nTotalCell + nBoundFace;

    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    int hLevel = grid->GetLevel(); 
    int pLevel     = parameters->pMultiGrid;
    if (pLevel == 1)
    {
        pLevel = 0;
    }
    int numberGaussPoint = 0;

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        //! Compute gauss point temperature of volume
        numberGaussPoint = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell].integNumberOfGaussPoint;
        for (int iGaussPoint = 0; iGaussPoint < numberGaussPoint; ++ iGaussPoint)
        {
            highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell].vist[iGaussPoint] = 0.0;
        }        
    }

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        //! Compute gauss point temperature of volume
        numberGaussPoint = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].integNumberOfGaussPoint;
        for (int iGaussPoint = 0; iGaussPoint < numberGaussPoint; ++ iGaussPoint)
        {
            highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace].vist[iGaussPoint] = 0.0;
        }        
    }
}

void HOSolverUnstruct::CompGamaAndTField(Grid *gridIn)
{
    int nl = 5;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);

    int nchem = 0;
    GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);
}

void HOSolverUnstruct::InitCoarseGridsFlow()
{
    if (HODebug) cout << "InitCoarseGridsFlow: Restrict from fine to coarse grid for all q set ghost cell q" << endl;

    Grid *fgrid = GetGrid(0);

    Grid *cgrid = fgrid->GetCoarseGrid();
    while (cgrid)
    {
        RestrictAllQ(fgrid, cgrid);
        Boundary(cgrid);

        fgrid = cgrid;
        cgrid = cgrid->GetCoarseGrid();
    }
}

//! Restrict from fine to coarse grid for all q
void HOSolverUnstruct::RestrictAllQ(Grid *fgrid, Grid *cgrid)
{
    if (HODebug) cout << "RestrictAllQ: Restrict from fine to coarse grid for all q" << endl;
}

void HOSolverUnstruct::ActionReflect(ActionTag *acttag)
{
}

void HOSolverUnstruct::UploadInterfaceData(ActionKey *actkey)
{
}

void HOSolverUnstruct::DownloadInterfaceData(ActionKey *actkey)
{
}

void HOSolverUnstruct::ZeroResiduals(Grid * gridIn)
{
    if (HODebug) cout << "ZeroResiduals:  \"res\" = 0" << endl;
    
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    int nl = 5;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);

    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    int hLevel = grid->GetLevel(); 
    int pLevel = parameters->pMultiGrid;
    
    if (pLevel == 1)
    {
        pLevel = 0;
    }

    int ndof = 0;
    int neqn = nl;

    RDouble *** res = reinterpret_cast< RDouble *** > (grid->GetDataPtr("res"));

    for (int ieqn = 0; ieqn < neqn; ++ ieqn)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            ndof = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell].nDOFsSol;

            for (int idof = 0; idof < ndof; ++ idof)
            {
                highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell].res[ieqn][idof] = 0.0;

                res[iCell][ieqn][idof] = 0.0;
            }
        }
    }
}

FieldProxy * HOSolverUnstruct::GetResidualProxy(Grid *gridIn)
{
    if (HODebug) cout << "GetResidualProxy: return proxy of \"res\"" << endl;

    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble ***res = reinterpret_cast< RDouble *** > (grid->GetDataPtr("res"));

    FieldProxy *res_proxy = new FieldProxy();

    res_proxy->SetField_UHO(res);

    return res_proxy;
}

FieldProxy * HOSolverUnstruct::CreateFieldProxy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int hLevel = grid->GetLevel();
    int nTotalCell = grid->GetNTotalCell();

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int nEquation  = nl + nchem;

    HighOrderGrid & highOrderGrid = highOrderGrids[hLevel];
    
    int * dofsOfCells = new int[nTotalCell];
    highOrderGrid.GetNumberOfDOFs(nTotalCell, dofsOfCells);

    RDouble *** field = NewPointer3<RDouble>(nTotalCell, nEquation, dofsOfCells);

    delete [] dofsOfCells;

    FieldProxy *fieldProxy = new FieldProxy();

    fieldProxy->SetField_UHO(field, true);

    return fieldProxy;
}

//! Load flow variables stored in grid to q
void HOSolverUnstruct::LoadQ(Grid *gridIn, FieldProxy *qProxy)
{
    if (HODebug) cout << "LoadQ: qProxy = \"q\"" << endl;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int hLevel = grid->GetLevel();
    int nTotalCell = grid->GetNTotalCell();
    //int nBoundFace = grid->GetNBoundFace();
    //int nTotal = nTotalCell + nBoundFace;

    RDouble ***qPro = qProxy->GetField_UHO();
    RDouble ***qold = reinterpret_cast< RDouble *** > (grid->GetDataPtr("q"));

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int nEquation  = nl + nchem;

    HighOrderGrid & highOrderGrid = highOrderGrids[hLevel];

    int * dofsOfCells = new int[nTotalCell];
    highOrderGrid.GetNumberOfDOFs(nTotalCell, dofsOfCells);

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            for (int iDof = 0; iDof < dofsOfCells[iCell]; ++ iDof)
            {
                qPro[iCell][m][iDof] = qold[iCell][m][iDof];
            }
        }
    }

    delete [] dofsOfCells;
}

void HOSolverUnstruct::Lhs(Grid *gridIn, FieldProxy *dqProxy, double coef)
{
    if (HODebug) cout << "Lhs: dqProxy *= dt * coef" << endl;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int hLevel = grid->GetLevel();
    int nTotalCell = grid->GetNTotalCell();

    RDouble ***dq = dqProxy->GetField_UHO();

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int nEquation  = nl + nchem;

    RDouble * dt = reinterpret_cast< RDouble * > (grid->GetDataPtr("dt" ));

    HighOrderGrid & highOrderGrid = highOrderGrids[hLevel];

    int * dofsOfCells = new int[nTotalCell];
    highOrderGrid.GetNumberOfDOFs(nTotalCell, dofsOfCells);
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            for (int iDof = 0; iDof < dofsOfCells[iCell]; ++ iDof)
            {
                dq[iCell][m][iDof] *= dt[iCell] * coef;
            }
        }
    }

    delete [] dofsOfCells;
}

void HOSolverUnstruct::StoreRhsByResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
    if (HODebug) cout << "StoreRhsByResidual:  rhsProxy = res" << endl;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int hLevel = grid->GetLevel();
    int nTotalCell = grid->GetNTotalCell();
    //int nBoundFace = grid->GetNBoundFace();
   // int nTotal = nTotalCell + nBoundFace;

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int nEquation  = nl + nchem;

    RDouble ***res = reinterpret_cast< RDouble *** > (grid->GetDataPtr("res"));
    RDouble ***rhs = rhsProxy->GetField_UHO();

    HighOrderGrid & highOrderGrid = highOrderGrids[hLevel];

    int * dofsOfCells = new int[nTotalCell];
    highOrderGrid.GetNumberOfDOFs(nTotalCell, dofsOfCells);

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            for (int iDof = 0; iDof < dofsOfCells[iCell]; ++ iDof)
            {
                rhs[iCell][m][iDof] = res[iCell][m][iDof];
            }
        }
    }

    delete [] dofsOfCells;
}

void HOSolverUnstruct::TimeStep(Grid *gridIn)
{
    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    int pLevel     = parameters->pMultiGrid;
    if (pLevel == 1)
    {
        pLevel = 0;
    }

    //! Compute varible in GaussPoint 
    ComputeGaussPointQ(gridIn);
    //! Compute varible of cellAverage 
    ComputeQAverage(gridIn);

    int ifLocalTimeStep = parameters->GetIfLocalTimeStep();

    if (ifLocalTimeStep == 0)
    {
        LocalTimeStep(gridIn);
    }
    else if (ifLocalTimeStep == 1)
    {
        GlobalTimeStep(gridIn);
    }
    else
    {
        LocalGlobalTimeStep(gridIn);
    }

}

void HOSolverUnstruct::ComputeMinTimeStep(Grid *gridIn, RDouble &minDt, RDouble &maxDt)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *dt = reinterpret_cast< RDouble *  > (grid->GetDataPtr("dt"));

    minDt = LARGE;
    maxDt = - LARGE;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        minDt = MIN(minDt, dt[iCell]);
        maxDt = MAX(maxDt, dt[iCell]);
    }
}

void HOSolverUnstruct::ReduceMaxTimeStep(Grid *grid, RDouble globalMinDt)
{
    if (HODebug) cout << "ReduceMaxTimeStep" << endl;
}

void HOSolverUnstruct::LoadResiduals(Grid *gridIn, FieldProxy *rhsProxy)
{
    if (HODebug) cout << "LoadResiduals  \"res\" = - rhsProxy" << endl;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int hLevel = grid->GetLevel();
    int nTotalCell = grid->GetNTotalCell();
    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    int pLevel = parameters->pMultiGrid;
    if (pLevel == 1)
    {
        pLevel = 0;
    }

    RDouble ***rhs = rhsProxy->GetField_UHO();

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int nEquation  = nl + nchem;

    RDouble ***res = reinterpret_cast< RDouble *** > (grid->GetDataPtr("res"));

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        HighOrderCellSol & uHOCellSol = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell];
        for (int m = 0; m < nEquation; ++ m)
        {
            for (int iDof = 0; iDof < uHOCellSol.nDOFsSol; ++ iDof)
            {
                res[iCell][m][iDof] = -rhs[iCell][m][iDof];
            }
        }
    }
}

//! Drive functions to calculate residuals for different orders and dimensions
void HOSolverUnstruct::UpdateResiduals(Grid *grid)
{
    if (HODebug) cout << "UpdateResiduals" << endl;

    ComputeGaussPointQ(grid);
    
    ComputeQAverage(grid);

    Boundary(grid);

    //CompGamaAndTField(grid);

    RightHandSide(grid);
}

void HOSolverUnstruct::RightHandSide(Grid *grid)
{
    if (HODebug) cout << "RightHandSide" << endl;

    InviscidFlux(grid);

    ViscousFlux(grid);

    SourceFlux(grid);

    //ZeroResidualOfSpecialCells(grid); //Added By Guo Yongheng 20160825
    //delete simu_ts_RightHandSide;

    ResDivideVol(grid);
}

//! Drive individaul functions to calculate inviscid fluxes
void HOSolverUnstruct::InviscidFlux(Grid * gridIn)
{
    if (HODebug) cout << "InviscidFlux" << endl;

    InviscidFluxForFace(gridIn);

    InviscidFluxForCell(gridIn);
}

void HOSolverUnstruct::ViscousFlux(Grid * gridIn)
{
}

void HOSolverUnstruct::SourceFlux(Grid * gridIn)
{
}

void HOSolverUnstruct::InviscidFluxForFace(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int      nTotalFace = grid->GetNTotalFace();
    int      nBoundFace = grid->GetNBoundFace();
    int      *leftCellOfFace    = grid->GetLeftCellOfFace();
    int      *rigtCellOfFace    = grid->GetRightCellOfFace();

    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    int hLevel = grid->GetLevel(); 
    int pLevel = parameters->pMultiGrid;
    if (pLevel == 1)
    {
        pLevel = 0;
    }

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int nEquation  = nl + nchem;

    RDouble ***res = reinterpret_cast< RDouble *** > (grid->GetDataPtr("res"));
    
    RDouble **flux  = new RDouble *[nEquation];
    RDouble **primL = new RDouble *[nEquation];
    RDouble **primR = new RDouble *[nEquation];
    RDouble **fluxLeftTotal = new RDouble * [nEquation];
    RDouble **fluxRigtTotal = new RDouble * [nEquation];
    
    int maxDof= 100; //! A arbitrarily number of dof related to accuracy!!!
                     //! The number need to be modified of different accuracy!!!
    fluxLeftTotal[0] = new RDouble[nEquation*maxDof];
    fluxRigtTotal[0] = new RDouble[nEquation*maxDof];
    primL[0]         = new RDouble[nEquation * (highOrderGrids[hLevel].totalFaceGaussPoint[pLevel])];
    primR[0]         = new RDouble[nEquation * (highOrderGrids[hLevel].totalFaceGaussPoint[pLevel])];
    flux[0]          = new RDouble[nEquation * (highOrderGrids[hLevel].totalFaceGaussPoint[pLevel])];

    for (int i=1; i<nEquation; i++)
    {
        fluxLeftTotal[i] = &(fluxLeftTotal[i-1][maxDof]);
        fluxRigtTotal[i] = &(fluxRigtTotal[i-1][maxDof]);
        primL[i]         = &(primL[i-1][highOrderGrids[hLevel].totalFaceGaussPoint[pLevel]]);
        primR[i]         = &(primR[i-1][highOrderGrids[hLevel].totalFaceGaussPoint[pLevel]]);
        flux[i]          = &(flux[i-1][highOrderGrids[hLevel].totalFaceGaussPoint[pLevel]]);
    }

    int numberGaussPoint = 0, rGaussPoint = 0, nLeftDof = 0, nRigtDof = 0, curGaussPoint = 0, offsetGaussPoint = 0;
    int leftCell = 0, rigtCell = 0, faceVTKType = 0;
    
    int uns_scheme = GlobalDataBase::GetIntParaFromDB("uns_scheme");

    RDouble * xNormal = new RDouble [highOrderGrids[hLevel].totalFaceGaussPoint[pLevel]];
    RDouble * yNormal = new RDouble [highOrderGrids[hLevel].totalFaceGaussPoint[pLevel]];
    RDouble * zNormal = new RDouble [highOrderGrids[hLevel].totalFaceGaussPoint[pLevel]];

    curGaussPoint = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        HighOrderFaceSol & uHOFaceSol = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace];
        HighOrderFace & uHOFace       = highOrderGrids[hLevel].highOrderFaces[iFace];
        numberGaussPoint = uHOFaceSol.integNumberOfGaussPoint;
        offsetGaussPoint = uHOFaceSol.offsetGaussPoint;
        faceVTKType      = uHOFace.VTK_Type;
        
        for (int iGaussPoint = 0; iGaussPoint < numberGaussPoint; ++ iGaussPoint)
        {
            //! Compute conservation variable
            curGaussPoint = offsetGaussPoint + iGaussPoint;
            primL[0][curGaussPoint] = uHOFaceSol.q[0][iGaussPoint];
            primL[1][curGaussPoint] = uHOFaceSol.q[1][iGaussPoint]/primL[0][curGaussPoint];
            primL[2][curGaussPoint] = uHOFaceSol.q[2][iGaussPoint]/primL[0][curGaussPoint];
            primL[3][curGaussPoint] = uHOFaceSol.q[3][iGaussPoint]/primL[0][curGaussPoint];
            primL[4][curGaussPoint] = uHOFaceSol.q[4][iGaussPoint];
            primL[4][curGaussPoint] = HO_GAMAM1 * (primL[4][curGaussPoint] - 0.5 * primL[0][curGaussPoint] * 
                (pow(primL[1][curGaussPoint],2) + pow(primL[2][curGaussPoint],2) + pow(primL[3][curGaussPoint],2)));

            if (primL[4][curGaussPoint] < 0)
            {
                cout << "Negative Pressure of Face:  " << iFace << " P L =  " << primL[4][curGaussPoint] << endl;
            }

            rGaussPoint = iGaussPoint + numberGaussPoint;
            primR[0][curGaussPoint] = uHOFaceSol.q[0][rGaussPoint];
            primR[1][curGaussPoint] = uHOFaceSol.q[1][rGaussPoint]/primR[0][curGaussPoint];
            primR[2][curGaussPoint] = uHOFaceSol.q[2][rGaussPoint]/primR[0][curGaussPoint];
            primR[3][curGaussPoint] = uHOFaceSol.q[3][rGaussPoint]/primR[0][curGaussPoint];
            primR[4][curGaussPoint] = uHOFaceSol.q[4][rGaussPoint];
            primR[4][curGaussPoint] = HO_GAMAM1 * (primR[4][curGaussPoint] - 0.5 * primR[0][curGaussPoint] * 
                (pow(primR[1][curGaussPoint],2) + pow(primR[2][curGaussPoint],2) + pow(primR[3][curGaussPoint],2)));

            if (primR[4][curGaussPoint]<0)
            {
                cout << "Negative Pressure of Face:  " << iFace << " P R =  " << primR[4][curGaussPoint] << endl;
            }

            xNormal[curGaussPoint] = uHOFaceSol.metricNormalsFace[iGaussPoint][0];
            yNormal[curGaussPoint] = uHOFaceSol.metricNormalsFace[iGaussPoint][1];
            zNormal[curGaussPoint] = uHOFaceSol.metricNormalsFace[iGaussPoint][2];
        }
    }

    switch(uns_scheme)
    {
        case ISCHEME_ROE:
            break;
        case ISCHEME_VANLEER:
            Vanleer_Scheme(highOrderGrids[hLevel].totalFaceGaussPoint[pLevel], primL, primR, xNormal, yNormal, zNormal, HO_GAMA, flux);
            break;
        case ISCHEME_HLLE:
            break;
        default:
            break;
    }

    RDouble fluxTmp = 0.0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        HighOrderFaceSol & uHOFaceSol = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace];
        HighOrderFace & uHOFace       = highOrderGrids[hLevel].highOrderFaces[iFace];        
        numberGaussPoint = uHOFaceSol.integNumberOfGaussPoint;
        offsetGaussPoint = uHOFaceSol.offsetGaussPoint;
        faceVTKType      = uHOFace.VTK_Type;
        //const RDouble *weightCoeff;        
        leftCell = leftCellOfFace[iFace];
        HighOrderCellSol & uHOCellSol = highOrderGrids[hLevel].highOrderCellsSol[pLevel][leftCell];
        nLeftDof    =   uHOCellSol.nDOFsSol; 
        const RDouble *weightCoeff = highOrderGrids[hLevel].standardFaceElements[uHOFace.standardElementIndex].GetWtIntegration();
        for (int ieqn = 0; ieqn < nEquation; ++ ieqn)
        {
            for (int idof = 0; idof < nLeftDof; ++ idof)
            {
                fluxLeftTotal[ieqn][idof] = 0.0;
                for (int iGaussPoint = 0; iGaussPoint < numberGaussPoint; ++ iGaussPoint)
                {
                    curGaussPoint = offsetGaussPoint + iGaussPoint;
                    fluxTmp = flux[ieqn][curGaussPoint] * weightCoeff[iGaussPoint] * uHOFaceSol.JacDetFaceIntegration[iGaussPoint];
                    fluxLeftTotal[ieqn][idof] += fluxTmp * uHOFaceSol.leftBasisFuncFaceIntegration[iGaussPoint][idof];    
                }
                res[leftCell][ieqn][idof] -= fluxLeftTotal[ieqn][idof];

                if ((res[leftCell][ieqn][idof] != res[leftCell][ieqn][idof]) || fabs(res[leftCell][ieqn][idof]) > 1.0e20)
                {
                    cout << "iface res " << iFace << " = " << res[leftCell][ieqn][idof] << endl;
                }
            }
        }
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        HighOrderFaceSol & uHOFaceSol = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace];
        HighOrderFace & uHOFace       = highOrderGrids[hLevel].highOrderFaces[iFace];
        numberGaussPoint =   uHOFaceSol.integNumberOfGaussPoint;
        offsetGaussPoint =   uHOFaceSol.offsetGaussPoint;
        faceVTKType      =   uHOFace.VTK_Type;
        //HighOrderFaceSoldµÄtype
        leftCell    =   leftCellOfFace[iFace];
        HighOrderCellSol & uHOLeftCellSol       =   highOrderGrids[hLevel].highOrderCellsSol[pLevel][leftCell];
        nLeftDof    =   uHOLeftCellSol.nDOFsSol;

        rigtCell    =   rigtCellOfFace[iFace];
        HighOrderCellSol & uHORigtCellSol       =   highOrderGrids[hLevel].highOrderCellsSol[pLevel][rigtCell];
        nRigtDof    =   uHORigtCellSol.nDOFsSol;
        const RDouble *weightCoeff  = highOrderGrids[hLevel].standardFaceElements[uHOFace.standardElementIndex].GetWtIntegration();
        for (int ieqn = 0; ieqn < nEquation; ++ ieqn)
        {
            for (int idof = 0; idof < nLeftDof; ++ idof)
            {
                fluxLeftTotal[ieqn][idof] = 0.0;
                for (int iGaussPoint = 0; iGaussPoint < numberGaussPoint; ++ iGaussPoint)
                {
                    curGaussPoint = offsetGaussPoint + iGaussPoint;  
                    fluxTmp = flux[ieqn][curGaussPoint] * weightCoeff[iGaussPoint] * uHOFaceSol.JacDetFaceIntegration[iGaussPoint];
                    fluxLeftTotal[ieqn][idof] += fluxTmp * uHOFaceSol.leftBasisFuncFaceIntegration[iGaussPoint][idof];
                }
                res[leftCell][ieqn][idof] -= fluxLeftTotal[ieqn][idof];
            }
        }
        for (int ieqn = 0; ieqn < nEquation; ++ ieqn)
        {
            for (int idof = 0; idof < nRigtDof; ++ idof)
            {
                fluxRigtTotal[ieqn][idof] = 0.0;
                for (int iGaussPoint = 0; iGaussPoint < numberGaussPoint; ++ iGaussPoint)
                {
                    curGaussPoint = offsetGaussPoint + iGaussPoint;
                    fluxTmp = flux[ieqn][curGaussPoint] * weightCoeff[iGaussPoint] * uHOFaceSol.JacDetFaceIntegration[iGaussPoint];
                    fluxRigtTotal[ieqn][idof] += fluxTmp * uHOFaceSol.rightBasisFuncFaceIntegration[iGaussPoint][idof];
                }
                res[rigtCell][ieqn][idof] += fluxRigtTotal[ieqn][idof];

                if ((res[rigtCell][ieqn][idof] != res[rigtCell][ieqn][idof]) || fabs(res[rigtCell][ieqn][idof]) > 1.0e20)
                {
                    cout << "iface res " << iFace << " = " << res[rigtCell][ieqn][idof] << endl;
                }
            }
        }
    }

    delete [] flux[0];
    delete [] primL[0];
    delete [] primR[0];
    delete [] flux;
    delete [] primL;
    delete [] primR;
    delete [] fluxLeftTotal[0];
    delete [] fluxLeftTotal;
    delete [] fluxRigtTotal[0];
    delete [] fluxRigtTotal;
    delete [] xNormal;
    delete [] yNormal;
    delete [] zNormal;
}

void HOSolverUnstruct::InviscidFluxForCell(Grid *grid_in)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);
    int      nTotalCell = grid->GetNTotalCell();
    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    
    int hLevel = grid->GetLevel();
    int pLevel     = parameters->pMultiGrid;
    if (pLevel == 1)
    {
        pLevel = 0;
    }

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int nEquation  = nl + nchem;

    RDouble ***res = reinterpret_cast< RDouble *** > (grid->GetDataPtr("res"));

    RDouble *fluxTotal = new RDouble [nEquation];

    int numberGaussPoint, cellVTKType, ndof;
    RDouble flux[5][3], gradFia[3], ro, u, v, w, p, e, rou, rov, row, fluxTmp;
    
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        HighOrderCellSol & uHOCellSol = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell];
        HighOrderCell & uHOCell       = highOrderGrids[hLevel].highOrderCells[iCell];
        numberGaussPoint              = uHOCellSol.integNumberOfGaussPoint;
        cellVTKType                   = uHOCell.VTK_Type;
        const RDouble * weightCoeff = highOrderGrids[hLevel].standardCellElements[uHOCell.standardElementIndex].GetWtIntegration();
        ndof = uHOCellSol.nDOFsSol;
        for (int idof = 0; idof < ndof; ++ idof)
        {
            fluxTotal[0] = 0.0;
            fluxTotal[1] = 0.0;
            fluxTotal[2] = 0.0;
            fluxTotal[3] = 0.0;
            fluxTotal[4] = 0.0;
            for (int iGaussPoint = 0; iGaussPoint < numberGaussPoint; ++ iGaussPoint)
            {
                gradFia[0] = uHOCellSol.dxBasisFuncCellIntegration[iGaussPoint][idof];
                gradFia[1] = uHOCellSol.dyBasisFuncCellIntegration[iGaussPoint][idof];
                gradFia[2] = uHOCellSol.dzBasisFuncCellIntegration[iGaussPoint][idof];

                ro  = uHOCellSol.q[0][iGaussPoint];
                rou = uHOCellSol.q[1][iGaussPoint];
                rov = uHOCellSol.q[2][iGaussPoint];
                row = uHOCellSol.q[3][iGaussPoint];
                e   = uHOCellSol.q[4][iGaussPoint];

                u   = rou/ro;
                v   = rov/ro;
                w   = row/ro;
                p   = HO_GAMAM1 * (e - 0.5*(rou*rou + rov*rov + row*row)/ro);

                flux[0][0] = rou;
                flux[0][1] = rov;
                flux[0][2] = row;

                flux[1][0] = rou*rou/ro + p;
                flux[1][1] = rou*rov/ro;
                flux[1][2] = rou*row/ro;
                    
                flux[2][0] = rov*rou/ro;
                flux[2][1] = rov*rov/ro + p;
                flux[2][2] = rov*row/ro;
                    
                flux[3][0] = row*rou/ro;
                flux[3][1] = row*rov/ro;
                flux[3][2] = row*row/ro + p;

                flux[4][0] = (e + p)*u;
                flux[4][1] = (e + p)*v;
                flux[4][2] = (e + p)*w;

                fluxTmp = weightCoeff[iGaussPoint] * uHOCellSol.JacDetCellIntegration[iGaussPoint];

                for (int ieqn = 0; ieqn < nEquation; ++ ieqn)
                {
                    fluxTotal[ieqn] += fluxTmp * (gradFia[0] * flux[ieqn][0] + gradFia[1] * flux[ieqn][1] + gradFia[2] * flux[ieqn][2]);
                }
            }

            for (int ieqn = 0; ieqn < nEquation; ++ ieqn)
            {
                res[iCell][ieqn][idof] += fluxTotal[ieqn];
            }
        }
    }

    delete [] fluxTotal;
}

void HOSolverUnstruct::ResDivideVol(Grid *grid_in)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);
    int      nTotalCell = grid->GetNTotalCell();
    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);

    int hLevel = grid->GetLevel();
    int pLevel     = parameters->pMultiGrid;
    if (pLevel == 1)
    {
        pLevel = 0;
    }

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int nEquation  = nl + nchem;

    RDouble ***res = reinterpret_cast< RDouble *** > (grid->GetDataPtr("res"));

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        HighOrderCellSol & uHOCellSol = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell];
        HighOrderCell & uHOCell       = highOrderGrids[hLevel].highOrderCells[iCell];
        int ndof = uHOCellSol.nDOFsSol;
        for (int iDof = 0; iDof < ndof; ++ iDof)
        {
            for (int ieqn = 0; ieqn < nEquation; ++ ieqn)
            {
                res[iCell][ieqn][iDof] /= uHOCell.volume;
            }
        }
    }
}

void HOSolverUnstruct::ForwardLUSGS(Grid *grid_in)
{

}

void HOSolverUnstruct::SolveLUSGS(Grid *gridIn, FieldProxy *dq_proxy)
{
}


void HOSolverUnstruct::SolveLUSGS(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *rhsProxy, int sweeps, double epsilon)
{
    if (HODebug) cout << "SolveLUSGS solve dq" << endl;
}

//! Explicit multi-stage time marching method
void HOSolverUnstruct::RungeKutta(Grid *grid, FieldProxy *rhsProxy)
{
    int RKStage = 1;
    GlobalDataBase::GetData("RKStage", &RKStage, PHINT, 1);

    FieldProxy *resProxy = GetResidualProxy(grid);
    FieldProxy *qProxy = CreateFieldProxy(grid);
    FieldProxy *qProxyTemp = CreateFieldProxy(grid);

    if (grid->GetLevel() == 0)
    {
        RDouble *lamda = new RDouble[RKStage];
        GlobalDataBase::GetData("lamda", lamda, PHDOUBLE, RKStage);

        LoadQ(grid, qProxy);
        FillField(grid, qProxyTemp, qProxy);
        for (int istage = 0; istage < RKStage; ++ istage)
        {
            LoadResiduals(grid, rhsProxy);            
            UpdateResiduals(grid);
            Lhs(grid, resProxy, lamda[istage]);
            FillField(grid, qProxy, qProxyTemp);
            UpdateFlowField(grid, qProxy, resProxy);
        }

        delete [] lamda;
    }
    else
    {
        LoadQ(grid, qProxy);
        LoadResiduals(grid, rhsProxy);
        UpdateResiduals(grid);
        Lhs(grid,resProxy,1.0);
        UpdateFlowField(grid, qProxy, resProxy);
    }

    delete qProxy;
    delete qProxyTemp;
    delete resProxy;
}

void HOSolverUnstruct::FillField(Grid *gridIn, FieldProxy *field1Proxy, FieldProxy *field2Proxy)
{
    if (HODebug) cout << "FillField: field1Proxy = field2Proxy" << endl;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int hLevel = grid->GetLevel();
    int nTotalCell = grid->GetNTotalCell();

    RDouble ***field1 = field1Proxy->GetField_UHO();
    RDouble ***field2 = field2Proxy->GetField_UHO();

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int nEquation  = nl + nchem;

    HighOrderGrid & highOrderGrid = highOrderGrids[hLevel];

    int * dofsOfCells = new int[nTotalCell];
    highOrderGrid.GetNumberOfDOFs(nTotalCell, dofsOfCells);

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            for (int iDof = 0; iDof < dofsOfCells[iCell]; ++ iDof)
            {
                field1[iCell][m][iDof] = field2[iCell][m][iDof];
            }
        }
    }

    delete [] dofsOfCells;
}

void HOSolverUnstruct::FillField(Grid *gridIn, FieldProxy *fieldProxy, RDouble value)
{
    if (HODebug) cout << "FillField: fieldProxy = value" << endl;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int hLevel = grid->GetLevel();
    int nTotalCell = grid->GetNTotalCell();

    RDouble ***field = fieldProxy->GetField_UHO();

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int nEquation  = nl + nchem;

    HighOrderGrid & highOrderGrid = highOrderGrids[hLevel];

    int * dofsOfCells = new int[nTotalCell];
    highOrderGrid.GetNumberOfDOFs(nTotalCell, dofsOfCells);

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            for (int iDof = 0; iDof < dofsOfCells[iCell]; ++ iDof)
            {
                field[iCell][m][iDof] = value;
            }
        }
    }

    delete [] dofsOfCells;
}

void HOSolverUnstruct::UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy)
{
    if (HODebug) cout << "UpdateFlowField solve dq" << endl;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int hLevel = grid->GetLevel();
    int nTotalCell = grid->GetNTotalCell();
    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    int pLevel = parameters->pMultiGrid;
    if (pLevel == 1)
    {
        pLevel = 0;
    }

    RDouble ***q = reinterpret_cast< RDouble *** > (grid->GetDataPtr("q"));

    RDouble *** dq = dqProxy->GetField_UHO();

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int nEquation  = nl + nchem;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
         HighOrderCellSol & uHOCellSol = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell];
        for (int ieqn = 0; ieqn < nEquation; ++ ieqn)
        {
            for (int iDof = 0; iDof < uHOCellSol.nDOFsSol; ++ iDof)
            {
                q[iCell][ieqn][iDof] += dq[iCell][ieqn][iDof];
                highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell].dofQ[ieqn][iDof] = q[iCell][ieqn][iDof];
            }
        }
    }
}

void HOSolverUnstruct::RecoverResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
}

int HOSolverUnstruct::GetNPostSolve()
{
    return 3;
}

extern ResidualGlobal residual;
extern MaxResidual *maxResidual;

void HOSolverUnstruct::GetResidual(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));

    RDouble ***res = reinterpret_cast< RDouble *** > (grid->GetDataPtr("res"));

    int nl = 5;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);

    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int ntmodel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    int neqn = nl + nchem + ntmodel - 1;

    int nTotalCell = grid->GetNTotalCell();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    vector< RDouble > residualNormal(neqn);

    //! Initialization.
    for (int m = 0; m < neqn; ++ m)
    {
        residualNormal[m] = 0.0;
    }

    RDouble maximumResidual = maxResidual->GetMaxRes();

    RDouble maxResLocal = maximumResidual;
    RDouble maxResCoorXLocal, maxResCoorYLocal, maxResCoorZLocal;
    int maxResIndexMLocal;

    bool maxResChange = false;

    int * cellIBlank = grid->GetBlankIndex();

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (cellIBlank[ iCell ] != 1) continue;

        for (int m = 0; m < neqn; ++ m)
        {
            RDouble res2 = res[iCell][m][0] * res[iCell][m][0];
            residualNormal[m] += res2;

            if (res2 > maxResLocal)
            {
                maxResChange = true;

                maxResLocal = res2;

                maxResCoorXLocal = xcc[iCell];
                maxResCoorYLocal = ycc[iCell];
                maxResCoorZLocal = zcc[iCell];

                maxResIndexMLocal = m;
            }
        }
    }

    if (maxResChange)
    {
        maxResidual->SetMaxRes(maxResLocal);
        maxResidual->SetMaxResCoorX(maxResCoorXLocal);
        maxResidual->SetMaxResCoorY(maxResCoorYLocal);
        maxResidual->SetMaxResCoorZ(maxResCoorZLocal);
        maxResidual->setMaxResVariableIndex(maxResIndexMLocal);
    }

    DataContainer *cdata = actkey->GetData();
    PHWrite(cdata, neqn);
    PHWrite(cdata, residualNormal, neqn);
}

LIB_EXPORT void HOSolverUnstruct::ComputePostVisualVariables(Post_Visual * postVisualization)
{
    UnstructGrid * grid = UnstructGridCast(postVisualization->GetGrid());

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    int nl = 5;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);
    int nchem = 0;
    GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);
    int neqn = nl + nchem;

    int hLevel = grid->GetLevel();

    ComputeQAverage(grid);

    RDouble ** qToOut = NewPointer2<RDouble>(neqn, nTotal);
    
    for (int ieqn = 0; ieqn < neqn; ++ ieqn)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            qToOut[ieqn][iCell] = highOrderGrids[hLevel].highOrderCellsSol[0][iCell].qAverage[ieqn];
        }
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble conserValue[5];

        for (int ieqn = 0; ieqn < neqn; ++ ieqn)
        {
            conserValue[ieqn] = qToOut[ieqn][iCell];
        }

        qToOut[1][iCell] /= qToOut[0][iCell];
        qToOut[2][iCell] /= qToOut[0][iCell];
        qToOut[3][iCell] /= qToOut[0][iCell];

        qToOut[4][iCell] = ComputePressure(conserValue);
    }

    int *left_cell_of_face  = grid->GetLeftCellOfFace();
    int *right_cell_of_face = grid->GetRightCellOfFace();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    UnstructBCSet **bcr = grid->GetBCRecord();

    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    int pLevel = parameters->pMultiGrid;

    if (pLevel == 1)
    {
        pLevel = 0;
    }

    if (hLevel > 1 && pLevel > 1)
    {
        cout << "the multigrid para is wrong for unstructure high order in HOSolverUnstruct::Boundary(Grid *gridIn)" << endl;
    }

    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble wallTemperature = parameters->GetWallTemperature();

    double refDimensionalTemperature = parameters->GetRefDimensionalTemperature();

    double tw = wallTemperature / refDimensionalTemperature;

    RDouble *prims, *primt;
    prims = new RDouble[nl+nchem];
    primt = new RDouble[nl+nchem];

    bool massFlow = false;
    RDouble massAlpha = 1.0;
    double massIn  = 1.0;
    double massOut = 1.0;

    if (massFlow)
    {
        //ComputeMassFlow(grid);
        GlobalDataBase::GetData("massIn", &massIn, PHDOUBLE, 1);
        GlobalDataBase::GetData("massOut", &massOut, PHDOUBLE, 1);
        massAlpha = massIn / massOut;
        cout << " massAlpha = " << massAlpha << "\n";
    }

    using namespace IDX;
    using namespace PHENGLEI;

    int le, re, bcType, numberGaussPoint;
    RDouble vgn = 0.0;
    RDouble xtn = 0.0, ytn = 0.0, ztn = 0.0;

    bool isViscous = parameters->IsViscous();

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = left_cell_of_face [ iFace ];
        re = right_cell_of_face[ iFace ];
        bcType = bcr[iFace]->GetKey();
        HighOrderFaceSol & highOrderFaceSol = highOrderGrids[hLevel].highOrderFacesSol[pLevel][iFace];
        numberGaussPoint = highOrderFaceSol.integNumberOfGaussPoint;

        if (IsInterface(bcType)) //! Interface
        {
            continue;
        }
        else if (bcType == SYMMETRY)
        {
            using namespace IDX;

            qToOut[IR][re] = qToOut[IR][le];
            qToOut[IP][re] = qToOut[IP][le];

            RDouble um = qToOut[IU][le];
            RDouble vm = qToOut[IV][le];
            RDouble wm = qToOut[IW][le];

            RDouble vn = xfn[iFace] * um + yfn[iFace] * vm + zfn[iFace] * wm - vgn;

            qToOut[IU][re] = um - two * xfn[iFace] * vn;
            qToOut[IV][re] = vm - two * yfn[iFace] * vn;
            qToOut[IW][re] = wm - two * zfn[iFace] * vn;

            continue;
        }
        else if (IsWall(bcType))
        {
            if (!isViscous)
            {
                using namespace IDX;

                qToOut[IR][re] = qToOut[IR][le];
                qToOut[IP][re] = qToOut[IP][le];

                RDouble um = qToOut[IU][le];
                RDouble vm = qToOut[IV][le];
                RDouble wm = qToOut[IW][le];

                RDouble vn = xfn[iFace] * um + yfn[iFace] * vm + zfn[iFace] * wm - vgn;

                qToOut[IU][re] = um - two * xfn[iFace] * vn;
                qToOut[IV][re] = vm - two * yfn[iFace] * vn;
                qToOut[IW][re] = wm - two * zfn[iFace] * vn;
            }
            else if (wallTemperature <= 0.0) //! Viscous adiabatic wall.
            {
                using namespace IDX;

                qToOut[IR][re] = qToOut[IR][le];
                qToOut[IP][re] = qToOut[IP][le];

                qToOut[IU][re] = - qToOut[IU][le] + two * xtn;
                qToOut[IV][re] = - qToOut[IV][le] + two * ytn;
                qToOut[IW][re] = - qToOut[IW][le] + two * ztn;
            }
            else //! Viscous iso-thermal wall
            {
                using namespace IDX;

                qToOut[IR][re] = qToOut[IR][le];
                qToOut[IP][re] = qToOut[IP][le];

                qToOut[IU][re] = - qToOut[IU][le] + two * xtn;
                qToOut[IV][re] = - qToOut[IV][le] + two * ytn;
                qToOut[IW][re] = - qToOut[IW][le] + two * ztn;

                RDouble omav = one;
                RDouble coefficientOfStateEquation = 1.0 / (HO_GAMA * refMachNumber * refMachNumber);

                RDouble tm = (1.0 / coefficientOfStateEquation) * qToOut[IP][le] / qToOut[IR][le];

                RDouble rg, pg, tg; //! g for ghost

                pg = qToOut[IP][le];
                tg = two * tw - tm;
                if (tg <= 0.0) tg = tw;

                rg = pg / (coefficientOfStateEquation * tg * omav);

                qToOut[IR][re] = rg;
                qToOut[IP][re] = pg;
            }

            continue;
        }
        else if (bcType == FARFIELD)
        {
            using namespace IDX;

            qToOut[IR][re] = qToOut[IR][le];
            qToOut[IP][re] = qToOut[IP][le];
            qToOut[IU][re] = qToOut[IU][le];
            qToOut[IV][re] = qToOut[IV][le];
            qToOut[IW][re] = qToOut[IW][le];
        }
        else if (bcType == INFLOW)
        {
            qToOut[IR][re] = qToOut[IR][le];
            qToOut[IP][re] = qToOut[IP][le];
            qToOut[IU][re] = qToOut[IU][le];
            qToOut[IV][re] = qToOut[IV][le];
            qToOut[IW][re] = qToOut[IW][le];

            continue;
        }
        else if (bcType == OUTFLOW)
        { 
            qToOut[IR][re] = qToOut[IR][le];
            qToOut[IP][re] = qToOut[IP][le];
            qToOut[IU][re] = qToOut[IU][le];
            qToOut[IV][re] = qToOut[IV][le];
            qToOut[IW][re] = qToOut[IW][le];

            continue;
        }
        else //if (bcType == PRESSURE_OUTLET)
        {
            qToOut[IR][re] = qToOut[IR][le];
            qToOut[IP][re] = qToOut[IP][le];
            qToOut[IU][re] = qToOut[IU][le];
            qToOut[IV][re] = qToOut[IV][le];
            qToOut[IW][re] = qToOut[IW][le];

            continue;
        }
    }
    delete [] prims;
    delete [] primt;


    using namespace IDX;
    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DENSITY))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DENSITY);
        postVisualization->UpdateVisualNodeVarPtr(varName, qToOut[IDX::IR]);
    }
    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_U))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_U);
        postVisualization->UpdateVisualNodeVarPtr(varName, qToOut[IDX::IU]);
    }
    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_V))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_V);
        postVisualization->UpdateVisualNodeVarPtr(varName, qToOut[IDX::IV]);
    }
    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_W))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_W);
        postVisualization->UpdateVisualNodeVarPtr(varName, qToOut[IDX::IW]);
    }
    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_PRESSURE))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_PRESSURE);
        postVisualization->UpdateVisualNodeVarPtr(varName, qToOut[IDX::IP]);
    }

    DelPointer2(qToOut);
}

void HOSolverUnstruct::PostSolve(Grid *gridIn, int stage, int level)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int hLevel = grid->GetLevel();
    int nTotalCell = grid->GetNTotalCell();

    HighOrderGrid & highOrderGrid = highOrderGrids[hLevel];

    int * dofsOfCells = new int[nTotalCell];
    highOrderGrid.GetNumberOfDOFs(nTotalCell, dofsOfCells);

    delete [] dofsOfCells;

    if (stage == 0)
    {
        CommunicationInterfaceData();
        return;
    }

    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    int intervalStepFlow = parameters->GetIntervalStepFlow();

    int intervalStepPlot = parameters->GetIntervalStepPlot();
    int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");

    bool isSubIterationDump = true;
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady == HO_UNSTEADY)
    {
        GlobalDataBase::GetData("isSubIterationDump", &isSubIterationDump, PHBOOL, 1);
    }

    if (stage == 1)
    {
        return;
    }

    if (outnstep % intervalStepFlow == 0 && level == 0 && isSubIterationDump)
    {
        //! To dump restart data for continual simulation.
        ActionKey *actkeyDumpRestartData = new ActionKey();
        FillActionKey(actkeyDumpRestartData, DUMP_RESTART, level);
        DumpRestartData(actkeyDumpRestartData);
        delete actkeyDumpRestartData;
    }

    if (outnstep % intervalStepPlot == 0 && isSubIterationDump)
    {
    }

    int intervalStepRes = parameters->GetIntervalStepRes();
    if (outnstep % intervalStepRes == 0 && (!isUnsteady || !isSubIterationDump))
    {
        //! To dump residual data.
        ActionKey *actkeyDumpResidual = new ActionKey();
        FillActionKey(actkeyDumpResidual, DUMP_RESIDUAL, level);
        DumpResidual(actkeyDumpResidual);
        delete actkeyDumpResidual;
    }

    int intervalStepForce = parameters->GetIntervalStepForce();
    if (intervalStepForce > 0)
    {
        if (outnstep % intervalStepForce == 0 && isSubIterationDump)
        {
            if (level == 0)
            {
                //! To dump airforce coeference data.
                ActionKey *actkeyDumpAirForceCoef = new ActionKey();
                FillActionKey(actkeyDumpAirForceCoef, DUMP_AIR_FORCE_COEF, level);
                DumpAirForceCoef(actkeyDumpAirForceCoef);
                delete actkeyDumpAirForceCoef;

                //! To dump wall airforce coeference data.
                ActionKey *actkeyDumpWallAircoef = new ActionKey();
                FillActionKey(actkeyDumpWallAircoef, DUMP_CP_DISTRI, 0);
                DumpCpDistri(actkeyDumpWallAircoef);
                delete actkeyDumpWallAircoef;

                //! To dump surface information data.
                ActionKey *actkeyDumpSurfaceInfo = new ActionKey();
                FillActionKey(actkeyDumpSurfaceInfo, DUMP_SURFACE_INFO, 0);
                DumpSurfaceInfo(actkeyDumpSurfaceInfo);
                delete actkeyDumpSurfaceInfo;
            }
        }
    }
}

int HOSolverUnstruct::GetSchemeID(const string &scheme_name)
{
    int scheme_id = ISCHEME_ROE;
    if (scheme_name.substr(0,3) == "roe")
    {
        scheme_id = ISCHEME_ROE;
    }
    else if (scheme_name.substr(0,7) == "ausmpw+")
    {
        scheme_id = ISCHEME_AUSMPW_PLUS;
    }
    else if (scheme_name.substr(0,6) == "ausmpw")
    {
        scheme_id = ISCHEME_AUSMPW;
    }
    else if (scheme_name.substr(0,6) == "ausm+w")
    {
        scheme_id = ISCHEME_AUSM_W;
    }
    else if (scheme_name.substr(0,7) == "ausm+up")
    {
        scheme_id = ISCHEME_AUSMPUP;
    }
    else if (scheme_name.substr(0,5) == "ausm+")
    {    
        scheme_id = ISCHEME_AUSMP;
    }
    else if (scheme_name.substr(0,6) == "ausmdv")
    {
        scheme_id = ISCHEME_AUSMDV;
    }
    else if (scheme_name.substr(0,14) == "vanleer_Rotate")
    {
        scheme_id = ISCHEME_Rotate;
    }
    else if (scheme_name.substr(0,7) == "vanleer")
    {
        scheme_id = ISCHEME_VANLEER;
    }
    else if (scheme_name.substr(0,6) == "steger")
    {
        scheme_id = ISCHEME_STEGER;
    }
    else if (scheme_name.substr(0,4) == "hlle")
    {
        scheme_id = ISCHEME_HLLE;
    }
    else if (scheme_name.substr(0,5) == "lax_f")
    {
        scheme_id = ISCHEME_LAX_FRIEDRICHS;
    }
    //! Bell 20130418 add
    else
    {
        TK_Exit::ExceptionExit("Error: this inv-scheme is not exist !\n", true);
    }
    
    return scheme_id;
}

void HOSolverUnstruct::LocalTimeStep(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble *volume  = grid->GetCellVolume();

    int nTotalCell = grid->GetNTotalCell();
    RDouble *dt   = reinterpret_cast< RDouble *  > (grid->GetDataPtr("dt"));

    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    int hLevel = gridIn->GetLevel(); 
    int pLevel     = parameters->pMultiGrid;
    if (pLevel == 1)
    {
        pLevel = 0;
    }

    double cfl = ComputeCFL(grid->GetLevel());

    //! Set dt to zero
    PHSPACE::SetField(dt, 0.0, nTotalCell);

    //! Inviscid time step.
    RDouble *invSpectrumRadius = new RDouble [nTotalCell];
    InvSpectrumRadius(gridIn, invSpectrumRadius);

    int nPolySol = 1;

    RDouble dtMin = 1.0e30;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        nPolySol  = highOrderGrids[hLevel].highOrderCellsSol[pLevel][iCell].nPolySol; //polynomial order
        nPolySol  = 2 * nPolySol + 1;
        dt[iCell] = cfl * volume[iCell] / invSpectrumRadius[iCell]/nPolySol;
        dtMin = min(dtMin, dt[iCell]);
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dt[iCell] = dtMin;
    }

    delete [] invSpectrumRadius;

    bool isViscous = parameters->IsViscous();
    if (isViscous)
    {
        //! Viscous time step.
        RDouble *visSpectrumRadius = new RDouble [nTotalCell];
        PHSPACE::SetField(visSpectrumRadius, 0.0, nTotalCell);

        VisSpectrumRadius(gridIn, visSpectrumRadius);

        RDouble *dtViscous = new RDouble [nTotalCell];
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            dtViscous[iCell] = cfl * volume[iCell] / visSpectrumRadius[iCell];
        }

        delete [] visSpectrumRadius;


        //! Consider the viscous time step.
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            dt[iCell] = dt[iCell] * dtViscous[iCell] / (dt[iCell] + dtViscous[iCell]);
        }

        delete [] dtViscous;
    }
}

void HOSolverUnstruct::GlobalTimeStep(Grid *gridIn)
{
}

void HOSolverUnstruct::LocalGlobalTimeStep(Grid *gridIn)
{
}

RDouble HOSolverUnstruct::ComputeCFL(int level)
{
    int iterationStep = 0;

    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);

    bool isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        //! Unsteady.
        GlobalDataBase::GetData("innstep", &iterationStep, PHINT, 1);
    }
    else
    {
        int nMGLevel = parameters->GetNMGLevel();
        //! Steady
        if (nMGLevel <= 1)
        {
            //! Finest grid.
            GlobalDataBase::GetData("outnstep", &iterationStep, PHINT, 1);
        }
        else
        {
            //! Coarse grid.
            GlobalDataBase::GetData("newnstep", &iterationStep, PHINT, 1);
        }
    }

    int cflNstep     = parameters->GetCFLVaryStep();
    RDouble cflStart = parameters->GetCFLStart();
    RDouble cflEnd   = parameters->GetCFLEnd();

    RDouble cfl = cflEnd;
    if (iterationStep < cflNstep)
    {
        cfl  = cflStart * (cflNstep - iterationStep) + cflEnd * iterationStep;
        cfl /= cflNstep;
    }

    return cfl;
}

void HOSolverUnstruct::InvSpectrumRadius(Grid *gridIn, RDouble *invSpectrum)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();

    RDouble rl, ul, vl, wl, pl, el, cl, vnl, rr, ur, vr, wr, pr, er, cr, vnr;
    RDouble specificHeatRatio, inviscidSpectrumRadius, pm, rm, vn, cm;
    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    int hLevel = gridIn->GetLevel(); 
    int pLevel     = parameters->pMultiGrid;
    if (pLevel == 1)
    {
        pLevel = 0;
    }
    int leftCell, rigtCell;

    PHSPACE::SetField(invSpectrum, 0.0, nTotalCell);

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        leftCell = leftCellofFace [iFace];
        HighOrderCellSol & uHOLeftCellSol  =   highOrderGrids[hLevel].highOrderCellsSol[pLevel][leftCell];     
        rl  = uHOLeftCellSol.qAverage[0];
        ul  = uHOLeftCellSol.qAverage[1]/uHOLeftCellSol.qAverage[0];
        vl  = uHOLeftCellSol.qAverage[2]/uHOLeftCellSol.qAverage[0];
        wl  = uHOLeftCellSol.qAverage[3]/uHOLeftCellSol.qAverage[0];
        el  = uHOLeftCellSol.qAverage[4];  
        pl  = HO_GAMAM1 * (el - 0.5 * rl *(ul * ul + vl * vl + wl * wl));
        cl  = sqrt(HO_GAMA * pl / rl);
        vnl = xfn[iFace] * ul + yfn[iFace] * vl + zfn[iFace] * wl;

        specificHeatRatio = HO_GAMA;
        pm = pl;
        rm = rl;
        vn = vnl;
        cm = sqrt(specificHeatRatio * pm / rm);

        inviscidSpectrumRadius = half * area[iFace] * (ABS(vn) + cm);

        invSpectrum[leftCell] += inviscidSpectrumRadius;
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        leftCell = leftCellofFace [iFace];
        rigtCell = rightCellofFace[iFace];
        HighOrderCellSol & uHOLeftCellSol  =   highOrderGrids[hLevel].highOrderCellsSol[pLevel][leftCell];
        HighOrderCellSol & uHORigtCellSol  =   highOrderGrids[hLevel].highOrderCellsSol[pLevel][rigtCell];      
        rl  = uHOLeftCellSol.qAverage[0];
        ul  = uHOLeftCellSol.qAverage[1]/uHOLeftCellSol.qAverage[0];
        vl  = uHOLeftCellSol.qAverage[2]/uHOLeftCellSol.qAverage[0];
        wl  = uHOLeftCellSol.qAverage[3]/uHOLeftCellSol.qAverage[0];
        el  = uHOLeftCellSol.qAverage[4];  
        pl  = HO_GAMAM1 * (el - 0.5 * rl *(ul * ul + vl * vl + wl * wl));
        cl  = sqrt(HO_GAMA * pl / rl);
        vnl = xfn[iFace] * ul + yfn[iFace] * vl + zfn[iFace] * wl;

        rr  = uHORigtCellSol.qAverage[0];
        ur  = uHORigtCellSol.qAverage[1]/uHORigtCellSol.qAverage[0];
        vr  = uHORigtCellSol.qAverage[2]/uHORigtCellSol.qAverage[0];
        wr  = uHORigtCellSol.qAverage[3]/uHORigtCellSol.qAverage[0];
        er  = uHORigtCellSol.qAverage[4];  
        pr  = HO_GAMAM1 * (er - 0.5 * rr *(ur * ur + vr * vr + wr * wr));
        cr  = sqrt(HO_GAMA * pr / rr);
        vnr = xfn[iFace] * ur + yfn[iFace] * vr + zfn[iFace] * wr;

        specificHeatRatio = HO_GAMA;
        pm = half * (pl + pr);
        rm = half * (rl + rr);
        vn = half * (vnl + vnr);
        cm = sqrt(specificHeatRatio * pm / rm);

        inviscidSpectrumRadius = half * area[iFace] * (ABS(vn) + cm);

        invSpectrum[leftCell] += inviscidSpectrumRadius;
        if (rigtCell < nTotalCell)
        {
            //! The re cell is the interior cell.
            invSpectrum[rigtCell] += inviscidSpectrumRadius;
        }
    }
}

void HOSolverUnstruct::VisSpectrumRadius(Grid *gridIn, RDouble *visSpectrum)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();

    RDouble *xfn    = grid->GetFaceNormalX();
    RDouble *yfn    = grid->GetFaceNormalY();
    RDouble *zfn    = grid->GetFaceNormalZ();
    RDouble *xcc    = grid->GetCellCenterX();
    RDouble *ycc    = grid->GetCellCenterY();
    RDouble *zcc    = grid->GetCellCenterZ();
    RDouble *area   = grid->GetFaceArea();
    RDouble *volume = grid->GetCellVolume();

    HOSolverUnstructParam * parameters = dynamic_cast< HOSolverUnstructParam * > (controlParameters);
    int hLevel = gridIn->GetLevel(); 
    int pLevel     = parameters->pMultiGrid;
    if (pLevel == 1)
    {
        pLevel = 0;
    }
    int leftCell, rigtCell;

    PHSPACE::SetField(visSpectrum, 0.0, nTotalCell);

    RDouble refReNumber = GetControlParameters()->GetRefReNumber();

    RDouble prl = GlobalDataBase::GetDoubleParaFromDB("prl"); //! default 0.72
    RDouble prt = GlobalDataBase::GetDoubleParaFromDB("prt"); //! default 0.9

    const RDouble fourthThird = 4.0 / 3.0;
    RDouble oprl = 1.0 / prl;
    RDouble oprt = 1.0 / prt;

    RDouble visLaminar, visTurb, density, viscosity, specificHeatRatio;
    RDouble coefOfSpectrum1, coefOfSpectrum2, coefOfSpectrum; 
    RDouble faceArea, faceAera2, viscousSpectrumRadius;
    RDouble dx, dy, dz, ds;
    
    int visSpetrumMethod = 2;
    if (visSpetrumMethod == 1)
    {
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            leftCell = leftCellofFace [iFace];
            rigtCell = rightCellofFace[iFace];

            HighOrderCellSol & uHOLeftCellSol  =   highOrderGrids[hLevel].highOrderCellsSol[pLevel][leftCell];
            HighOrderCellSol & uHORigtCellSol  =   highOrderGrids[hLevel].highOrderCellsSol[pLevel][rigtCell];      

            visLaminar = half * (uHOLeftCellSol.vislAverage + uHORigtCellSol.vislAverage);
            visTurb    = half * (uHOLeftCellSol.vistAverage + uHORigtCellSol.vistAverage);
            density    = half * (uHOLeftCellSol.qAverage[0] + uHORigtCellSol.qAverage[0]);
            viscosity  = visLaminar + visTurb;

            specificHeatRatio = HO_GAMA;
            coefOfSpectrum1 = fourthThird * (visLaminar + visTurb);
            coefOfSpectrum2 = specificHeatRatio * (visLaminar * oprl + visTurb * oprt);
            coefOfSpectrum  = two * PHSPACE::MAX(coefOfSpectrum1, coefOfSpectrum2) / (refReNumber * density);

            faceArea = area[iFace];
            faceAera2 = PHSPACE::SQR(faceArea);

            viscousSpectrumRadius = faceAera2 * coefOfSpectrum;

            visSpectrum[leftCell] += viscousSpectrumRadius / volume[leftCell];
            if (rigtCell < nTotalCell)
            {
                visSpectrum[rigtCell] += viscousSpectrumRadius / volume[rigtCell];
            }
        }
    }
    else if (visSpetrumMethod == 2)
    {
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            leftCell = leftCellofFace [iFace];
            rigtCell = rightCellofFace[iFace];

            HighOrderCellSol & uHOLeftCellSol  =   highOrderGrids[hLevel].highOrderCellsSol[pLevel][leftCell];
            HighOrderCellSol & uHORigtCellSol  =   highOrderGrids[hLevel].highOrderCellsSol[pLevel][rigtCell];      

            dx  = xcc[rigtCell] - xcc[leftCell];
            dy  = ycc[rigtCell] - ycc[leftCell];
            dz  = zcc[rigtCell] - zcc[leftCell];
            ds  = ABS(xfn[iFace] * dx + yfn[iFace] * dy + zfn[iFace] * dz);

            visLaminar = half * (uHOLeftCellSol.vislAverage + uHORigtCellSol.vislAverage);
            visTurb    = half * (uHOLeftCellSol.vistAverage + uHORigtCellSol.vistAverage);
            density    = half * (uHOLeftCellSol.qAverage[0] + uHORigtCellSol.qAverage[0]);
            viscosity  = visLaminar + visTurb;

            faceArea  = area[iFace];
            viscousSpectrumRadius  = 2.0 * viscosity / (density * ds * refReNumber + SMALL);
            viscousSpectrumRadius *= half * faceArea;

            visSpectrum[leftCell] += viscousSpectrumRadius;
            if (rigtCell < nTotalCell)
            {
                visSpectrum[rigtCell] += viscousSpectrumRadius;
            }
        }
    }
}

}