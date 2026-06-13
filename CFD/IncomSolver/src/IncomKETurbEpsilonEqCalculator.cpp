#include <iostream>
#include "IncomKETurbEpsilonEqCalculator.h"




using namespace std;

namespace PHSPACE
{
IncomKETurbEpsilonEqCalculator::IncomKETurbEpsilonEqCalculator():IncomScalarEqCalculator()
{

}

IncomKETurbEpsilonEqCalculator::~IncomKETurbEpsilonEqCalculator()
{

}

void IncomKETurbEpsilonEqCalculator::SetDiffusionCoeff(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    RDouble prandtl = 1.3;
    RDouble **DiffusionCoeff = reinterpret_cast <RDouble **> (GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    int solverIndex = GetSolverIndex();

    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        DiffusionCoeff[solverIndex][iCell] = visl[iCell] + vist[iCell] / prandtl;
    }

    for (int iFace = 0; iFace < nBoundFace; ++iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        DiffusionCoeff[solverIndex][re] = DiffusionCoeff[solverIndex][le];
    }

    CommunicateAnInterfaceVar(DiffusionCoeff[solverIndex]);
}

void IncomKETurbEpsilonEqCalculator::InitFlowAsRestart(Grid *gridIn)
{
    IncomScalarEqCalculator::InitFlowAsRestart(gridIn);
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    RDouble initPhi = GlobalDataBase::GetDoubleParaFromDB("init" + varNameIncom[solverIndex]);
    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(varNameIncom[solverIndex]));
    PHSPACE::SetField(phi, initPhi, nTotal);
}

void IncomKETurbEpsilonEqCalculator::solveScalarEquation(Grid *gridIn, int iEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int solverIndex = GetSolverIndex();
    CalcGrad(grid);
    SetDiffusionCoeff(grid);

    constructMatrixACoeff(grid, solverIndex);
    constructBCoeff(grid, solverIndex);

    WallModification(grid);

    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("upperMatrixCoeff"));
    RDouble *lowerMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("lowerMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));

    string *varNameIncom = reinterpret_cast <string *>(GlobalDataBase::GetDataPtr("varNameIncom"));
    string varName = varNameIncom[solverIndex];
    RDouble *q = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
    calculateLinearEquation(grid, solverIndex, q, diagMatrixCoeff, bCoeff, upperMatrixCoeff , lowerMatrixCoeff);

    UpdateBCValue(grid);
    UpdateAfterIterloop(grid);
    UpdateProperties(grid);
}

void IncomKETurbEpsilonEqCalculator::calWallBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff)
{
    
}

void IncomKETurbEpsilonEqCalculator::WallModification(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCells = grid->GetNTotalCell();
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFaces = grid->GetNBoundFace() - grid->GetNIFace();
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();
    int **cell2face  = grid->GetCell2Face();

    RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
    RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("upperMatrixCoeff"));
    RDouble *lowerMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("lowerMatrixCoeff"));
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();

        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            for (int i = 0; i < faceIndex->size(); ++ i)
            {
                int bface = (*faceIndex)[i];
                int le = leftCellOfFace[bface];

                bCoeff[le] = Epsilon[le];
                diagMatrixCoeff[le] = 1.0;

                for (int j = 0; j < grid->GetFaceNumberOfEachCell()[le]; ++j)
                {
                    int iFace = grid->GetCell2Face()[le][j];

                    if (iFace < nBoundFaces)
                    {
                        continue;
                    }
                    else
                    {
                        if (grid->GetLeftCellOfFace()[iFace] == le)
                        {
                            upperMatrixCoeff[iFace] = 0.0;
                        }
                        else
                        {
                            lowerMatrixCoeff[iFace] = 0.0;
                        }
                    }
                }
            }
        }
    }
    
}

void IncomKETurbEpsilonEqCalculator::CalcOtherMatrixACoeff(Grid * gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
    RDouble *vol = reinterpret_cast<RDouble *>(grid->GetCellVolume());
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));
    RDouble c1m_ = 1.44;
    RDouble c2m_ = 1.92;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        diagMatrixCoeff[iCell] += c2m_ * rho[iCell] * vol[iCell] * Epsilon[iCell] / (k[iCell] + 1e-20);
    }
}

void IncomKETurbEpsilonEqCalculator::CalcOtherbCoeff(Grid *gridIn, int iEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
    RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
    RDouble *vol = reinterpret_cast<RDouble *>(grid->GetCellVolume());
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));
    RDouble c1m_ = 1.44;
    RDouble c2m_ = 1.92;
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    CalcGenTerm(grid);

    UnstructBCSet* unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        if (2 == bcType)
        {
            vector<int>* faceIndex = bcRegion->GetFaceIndex();

            for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];

                gen[le] = 0.0;
                Epsilon[le] = 0.0;
            }
        }
    }

    UpdateWallGen(grid);

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble ke = vol[iCell] * Epsilon[iCell] / (k[iCell] + 1e-20);
        bCoeff[iCell] += c1m_ * gen[iCell] * ke;
        bCoeff[iCell] = max(bCoeff[iCell], 0.0);
    }
}

void IncomKETurbEpsilonEqCalculator::GetResidual(Grid *gridIn, vector<RDouble>& res)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble resNow = 0.0;

    grid->GetData("EpsilonResNow", &resNow, PHDOUBLE, 1);

    res.push_back(resNow);
}


void IncomKETurbEpsilonEqCalculator::InitialUnsteadyVar(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    PHString1D phiNameList;
    phiNameList.push_back("Epsilon");

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady == 0)
    {
        return;
    }

    string TranCalcMethod = GlobalDataBase::GetStrParaFromDB("TranCalcMethod");
    if (TranCalcMethod == "IMPLICIT_EULER")
    {
        ImplicitEuler_ReInitTimeVar(grid, phiNameList);
    }
    else if (TranCalcMethod == "IMPLICIT_2ND_ORDER")
    {
        Implicit2ndOrder_ReInitTimeVar(grid, phiNameList);
    }

    PHString1D().swap(phiNameList);
}

void IncomKETurbEpsilonEqCalculator::UpdateUnsteadyVariable(Grid *grid)
{
    UnstructGrid *gridIn = UnstructGridCast(grid); 
    std::vector<std::string> phiNameList;
    phiNameList.push_back("Epsilon");

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady == 0)
    {
        return;
    }

    string TranCalcMethod = GlobalDataBase::GetStrParaFromDB("TranCalcMethod");
    if (TranCalcMethod == "IMPLICIT_EULER")
    {
        ImplicitEuler_SaveOldTimeValue(gridIn, phiNameList);
    }
    else if (TranCalcMethod == "IMPLICIT_2ND_ORDER")
    {
        Implicit2ndOrder_SaveOldTimeValue(gridIn, phiNameList);
    }

}

void IncomKETurbEpsilonEqCalculator::UpdateBCValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
    RDouble *wDist = reinterpret_cast<RDouble *>(grid->GetWallDist());
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *yplus = reinterpret_cast<RDouble *>(grid->GetDataPtr("yplus")) + nTotalCell;
    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            for (int i = 0; i < faceIndex->size(); ++i)
            {
                int bface = (*faceIndex)[i];
                int le = grid->GetLeftCellOfFace()[bface];
                int re = rightCellOfFace[bface];
                phi[re] = 2 * phi[le];
            }
        }
        else if (bcType == PHENGLEI::SYMMETRY)
        {
            for (size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                phi[re] = phi[le];
            }
        }
        else if (bcType == PHENGLEI::FARFIELD)
        {
            RDouble epsilonb = 0.0;

            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                const int iFacelocal = (*faceIndex)[iFace];

                if (bcData->IsExist("initEpsilon", PHDOUBLE, 1))
                {
                    bcData->GetData("initEpsilon", &epsilonb, PHDOUBLE, 1);
                }
                else
                {
                    std::cout << "Please assign a value of Epsilon on inlet" << std::endl;;
                }
            }
            for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
                if (FaceFlux[iFacelocal] < 0)
                {

                    phi[re] = epsilonb;
                }
                else
                {

                    phi[re] = phi[le];
                }
            }
        }
        else if (bcType == PHENGLEI::INFLOW)
        {
            RDouble epsilonb = 0.0;

            RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                const int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                if (bcData->IsExist("initEpsilon", PHDOUBLE, 1))
                {
                    bcData->GetData("initEpsilon", &epsilonb, PHDOUBLE, 1);
                }
                else
                {
                    std::cout << "Please assign a value of Epsilon on inlet" << std::endl;;
                }
                if (FaceFlux[iFacelocal] < 0)
                {
                    phi[re] = epsilonb;
                }
                else
                {
                    phi[re] = phi[le];
                }
            }
        }
        else if (bcType == PHENGLEI::OUTFLOW)
        {
            RDouble epsilonb = 0.0;
            RDouble *FaceFlux= reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
            RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
            string varName = "Epsilon";
            RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dx"));
            RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dy"));
            RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dz"));
            
            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                if (FaceFlux[iFacelocal] < 0.0)
                {
                    if (bcData->IsExist("initEpsilon", PHDOUBLE, 1))
                    {
                        bcData->GetData("initEpsilon", &epsilonb, PHDOUBLE, 1);
                        phi[re] = RDouble(epsilonb);
                    }
                    else
                    {
                        std::cout << " Please assign a value " << std::endl;
                    }
                }
                else
                {
                    phi[re] = phi[le];
                }
            }
        }
        else if (bcType == PHENGLEI::PRESSURE_INLET)
        {
            RDouble epsilonb = 0.0;

            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                const int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                if (bcData->IsExist("initEpsilon", PHDOUBLE, 1))
                {
                    bcData->GetData("initEpsilon", &epsilonb, PHDOUBLE, 1);
                }
                else
                {
                    std::cout << "Please assign a value of Epsilon on inlet" << std::endl;;
                }
                if (FaceFlux[iFacelocal] < 0)
                {
                    phi[re] = epsilonb;
                }
                else
                {
                    phi[re] = phi[le];
                }
            }
        }
        else if (bcType == PHENGLEI::PRESSURE_OUTLET)
        {
            RDouble epsilonb = 0.0;
            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];

                if (FaceFlux[iFacelocal] < 0.0)
                {
                    if (bcData->IsExist("initEpsilon", PHDOUBLE, 1))
                    {
                        bcData->GetData("initEpsilon", &epsilonb, PHDOUBLE, 1);
                        phi[re] = RDouble(epsilonb);
                    }
                    else
                    {
                        std::cout << " Please assign a value " << std::endl;
                    }
                }
                else
                {
                    phi[re] = phi[le];
                }
            }
        }
        else 
        {
            for (size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                phi[re] = phi[le];
            }
        }
    }
    CommunicateAnInterfaceVar(phi);
}

void IncomKETurbEpsilonEqCalculator::UpdateProperties(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));

    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        if (fabs(Epsilon[iCell]) < 1e-20)
        {
            Epsilon[iCell] = 1e-20;
        }
        else if (Epsilon[iCell] < 0.0)
        {
            Epsilon[iCell] = -Epsilon[iCell];
        }
    }
}

void IncomKETurbEpsilonEqCalculator::CalcGenTerm(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *dudx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdx"));
    RDouble *dudy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdy"));
    RDouble *dudz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdz"));
    RDouble *dvdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdx"));
    RDouble *dvdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdy"));
    RDouble *dvdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdz"));
    RDouble *dwdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdx"));
    RDouble *dwdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdy"));
    RDouble *dwdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdz"));
    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
    RDouble *yplus = reinterpret_cast<RDouble *>(grid->GetDataPtr("yplus"));
    RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    
    GradientCalculation(grid, "U", "dUdx", "dUdy", "dUdz");
    GradientCalculation(grid, "V", "dVdx", "dVdy", "dVdz");
    GradientCalculation(grid, "W", "dWdx", "dWdy", "dWdz");

    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        RDouble _dudx = dudx[iCell];
        RDouble _dudy = dudy[iCell];
        RDouble _dudz = dudz[iCell];
        RDouble _dvdx = dvdx[iCell];
        RDouble _dvdy = dvdy[iCell];
        RDouble _dvdz = dvdz[iCell];
        RDouble _dwdx = dwdx[iCell];
        RDouble _dwdy = dwdy[iCell];
        RDouble _dwdz = dwdz[iCell];

        gen[iCell] = (vist[iCell] + visl[iCell]) * (2.0 * (pow(_dudx, 2) + pow(_dvdy, 2) + pow(_dwdz, 2)) +
            pow((_dvdx + _dudy), 2) + pow((_dwdx + _dudz), 2) + pow((_dwdy + _dvdz), 2));

        gen[iCell] = min(gen[iCell], 10.0 * rho[iCell] * Epsilon[iCell]);
    }

}

void IncomKETurbEpsilonEqCalculator::UpdateWallGen(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    UnstructBCSet* unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int>* faceIndex = bcRegion->GetFaceIndex();
        Data_Param* bcData = bcRegion->GetBCParamDataBase();

        if (bcType == PHENGLEI::SOLID_SURFACE)

        {
            int nTotalCell = grid->GetNTotalCell();
            int *leftCellOfFace = grid->GetLeftCellOfFace();
            int *rightCellOfFace = grid->GetRightCellOfFace();
            int nBoundFace = grid->GetNBoundFace();
            RDouble *xfn = grid->GetFaceNormalX();
            RDouble *yfn = grid->GetFaceNormalY();
            RDouble *zfn = grid->GetFaceNormalZ();
            RDouble *WD = reinterpret_cast<RDouble *>(grid->GetDataPtr("wd"));
            RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
            RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
            RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
            RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
            RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
            RDouble *vistb = vist + nTotalCell;
            RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
            RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
            RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
            RDouble *mu = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
            RDouble *yplus = reinterpret_cast<RDouble *>(grid->GetDataPtr("yplus")) + nTotalCell;
            RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
            RDouble *nu_eff_w = reinterpret_cast<RDouble *>(grid->GetDataPtr("nu_eff")) + nTotalCell;
            RDouble *wallFlag = reinterpret_cast<RDouble *>(grid->GetDataPtr("wallFlag"));

            RDouble cmu75_ = 0.1643;
            RDouble cappa_ = 0.4187;
            RDouble cmu = 0.09;
            RDouble cmu25_ = 0.5477;
            RDouble yptr_ = 11.225;
            RDouble yvstar_ = 20.0;
            RDouble cl_ = 2.55;
            RDouble elog_ = 9.793;

            for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];

                RDouble wd = WD[re];
                RDouble cmute = cmu25_ * sqrt(k[le]);
                yplus[iFacelocal] = rho[le] * wd * cmute / (visl[le]);

                RDouble du = u[re] - u[le];
                RDouble dv = v[re] - v[le];
                RDouble dw = w[re] - w[le];
                RDouble ax = grid->GetFaceNormalX()[iFacelocal];
                RDouble ay = grid->GetFaceNormalY()[iFacelocal];
                RDouble az = grid->GetFaceNormalZ()[iFacelocal];
                RDouble dun = (du * ax + dv * ay + dw * az);
                RDouble Ux = du - dun * ax;
                RDouble Uy = dv - dun * ay;
                RDouble Uz = dw - dun * az;
                RDouble Uwall = std::sqrt(Ux * Ux + Uy * Uy + Uz * Uz);//壁面切向速度

                yplus[iFacelocal] = max(yplus[iFacelocal], yptr_);
                double tau_w = Uwall * cmute * rho[le] / (std::log(elog_ * yplus[iFacelocal]) / cappa_);
                double d = yplus[iFacelocal] * visl[le] / (rho[le] * cmute);
                nu_eff_w[iFacelocal] = tau_w * wd / Uwall;
                gen[le] += (tau_w * tau_w / (cappa_ * rho[le] * cmute * d)) / wallFlag[le];
                Epsilon[le] += (cmu75_ * std::pow(k[le], 1.5) / (cappa_ * d)) / wallFlag[le];
            }
        }
    }
}

}