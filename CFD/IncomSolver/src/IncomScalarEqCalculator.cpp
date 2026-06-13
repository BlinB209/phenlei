#include "IncomScalarEqCalculator.h"



namespace PHSPACE
{
IncomScalarEqCalculator::IncomScalarEqCalculator():IncomCalculator()
{

}

IncomScalarEqCalculator::~IncomScalarEqCalculator()
{

}

void IncomScalarEqCalculator::GetResidual(Grid *gridIn, vector<RDouble>& res)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
}

void IncomScalarEqCalculator::InitFlowAsRestart(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    RDouble **DiffusionCoeff = reinterpret_cast <RDouble **> (GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    int solverIndex = GetSolverIndex();

    RDouble initMu = GlobalDataBase::GetDoubleParaFromDB("initMu");
    PHSPACE::SetField(DiffusionCoeff[solverIndex], initMu, nTotal);


    RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varNameIncom[solverIndex] + "dx"));
    RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varNameIncom[solverIndex] + "dy"));
    RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varNameIncom[solverIndex] + "dz"));

    PHSPACE::SetField(dphidx, 0.0, nTotal);
    PHSPACE::SetField(dphidy, 0.0, nTotal);
    PHSPACE::SetField(dphidz, 0.0, nTotal);
}

void IncomScalarEqCalculator::AllocateGlobalVar(Grid *gridIn)
{     
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    RDouble **DiffusionCoeff = reinterpret_cast <RDouble **> (GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    int solverIndex = GetSolverIndex();

    DiffusionCoeff[solverIndex] = NewPointer<RDouble>(nTotal);
    
    RDouble *phi = NewPointer<RDouble>(nTotal);
    grid->UpdateDataPtr(varNameIncom[solverIndex], phi);

    RDouble *dphidx = NewPointer<RDouble>(nTotal);
    RDouble *dphidy = NewPointer<RDouble>(nTotal);
    RDouble *dphidz = NewPointer<RDouble>(nTotal);

    grid->UpdateDataPtr("d" + varNameIncom[solverIndex] + "dx", dphidx);
    grid->UpdateDataPtr("d" + varNameIncom[solverIndex] + "dy", dphidy);
    grid->UpdateDataPtr("d" + varNameIncom[solverIndex] + "dz", dphidz);

    AllocateUnsteadyVar(grid);
}

void IncomScalarEqCalculator::AllocateUnsteadyVar(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    string *varNameIncom = reinterpret_cast<string *>(GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    string varName = varNameIncom[solverIndex];

    PHString1D phiNameList;
    phiNameList.push_back(varName);

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady == 0)
    {
        return;
    }

    string TranCalcMethod = GlobalDataBase::GetStrParaFromDB("TranCalcMethod");
    if (TranCalcMethod == "IMPLICIT_EULER")
    {
        ImplicitEuler_AllocateMemory(grid, phiNameList);
    }
    else if (TranCalcMethod == "IMPLICIT_2ND_ORDER")
    {
        Implicit2ndOrder_AllocateMemory(grid, phiNameList);
    }

    PHString1D().swap(phiNameList);
}

void IncomScalarEqCalculator::IncompressibleInitial(Grid * gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    UpdateBCValue(grid);

    InitialUnsteadyVar(grid);
}

void IncomScalarEqCalculator::UpdateUnsteadyVariable(Grid *grid)
{
}


void IncomScalarEqCalculator::solveScalarEquation(Grid *gridIn, int iEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();

    RDouble *diagMatrixCoeff  = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("upperMatrixCoeff"));
    RDouble *lowerMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("lowerMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));

    string *varNameIncom = reinterpret_cast <string *>(GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    string varName = varNameIncom[solverIndex];

    CalcGrad(grid);
    SetDiffusionCoeff(grid);

    constructMatrixACoeff(grid, solverIndex);

    constructBCoeff(grid, solverIndex);

    RDouble *q = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
    calculateLinearEquation(grid, solverIndex, q, diagMatrixCoeff, bCoeff, upperMatrixCoeff , lowerMatrixCoeff);

    UpdateBCValue(grid);
    UpdateAfterIterloop(grid);
    UpdateProperties(grid);
}

void IncomScalarEqCalculator::constructMatrixACoeff(Grid *gridIn, int iEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 

    int solverIndex = GetSolverIndex(iEquation);

    InitializeMatrixACoeff(grid);

    treatBC(grid, solverIndex, 1.0);

    FirstUpwindInvisTerm(grid);

    VisMatrixTerm(grid, solverIndex);

    calcTransMatrixTerm(grid);

    CalcOtherMatrixACoeff(grid);

    relaxMatrixTerm(grid, solverIndex);
}

void IncomScalarEqCalculator::constructBCoeff(Grid *gridIn, int iEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 

    int solverIndex = GetSolverIndex(iEquation);

    InitializeBCoeff(grid);

    treatBC(grid, solverIndex, 0.0);

    calcInvisSourceTerm(grid, solverIndex);

    VisSourceTerm(grid, solverIndex);

    calcTransSourceTerm(grid, solverIndex);

    CalcOtherbCoeff(grid, solverIndex);

    relaxSourceTerm(grid, solverIndex);
}

void IncomScalarEqCalculator::UpdateBCValue(Grid *gridIn)
{

}

void IncomScalarEqCalculator::CalcGrad(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    string varName = varNameIncom[solverIndex];
    RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dx"));
    RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dy"));
    RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dz"));

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        dphidx[iCell] = 0.0;
        dphidy[iCell] = 0.0;
        dphidz[iCell] = 0.0;
    }

    GradientCalculation(grid, varName, "d" + varName + "dx", "d" + varName + "dy",
        "d" + varName + "dz");
}


void IncomScalarEqCalculator::SetDiffusionCoeff(Grid *gridIn) 
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int solverIndex = this->GetSolverIndex();
    RDouble **DiffusionCoeff = reinterpret_cast<RDouble **>(GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    CommunicateAnInterfaceVar(DiffusionCoeff[solverIndex]);
}

void IncomScalarEqCalculator::treatBC(Grid *gridIn, int iVariable, RDouble matrixCoeff)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    RDouble **DiffusionCoeff = reinterpret_cast <RDouble **> (GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("bCoeff"));
    int solverIndex = GetSolverIndex(iVariable);

    RDouble *mu = DiffusionCoeff[solverIndex];
    string varName = varNameIncom[solverIndex];
    RDouble *q = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
    RDouble *dqdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dx"));
    RDouble *dqdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dy"));
    RDouble *dqdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dz"));

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            calWallBC(grid, faceIndex, bcData, diagMatrixCoeff, bCoeff, q, dqdx, dqdy, dqdz, mu, matrixCoeff);
        }
        else if (bcType == PHENGLEI::SYMMETRY)
        {
            if (iVariable == IDX::S_IU || iVariable == IDX::S_IV || iVariable == IDX::S_IW)
            {
                calSymmetryBC(grid, faceIndex, bcData, diagMatrixCoeff, bCoeff, q, dqdx, dqdy, dqdz, mu, matrixCoeff);
            }
        }
        else if (bcType == PHENGLEI::FARFIELD)
        {
            calFarfieldBC(grid, faceIndex, bcData, diagMatrixCoeff, bCoeff, q, dqdx, dqdy, dqdz, mu, matrixCoeff);
        }
        else if (bcType == PHENGLEI::INFLOW)
        {
            calVelocityInletBC(grid, faceIndex, bcData, diagMatrixCoeff, bCoeff, q, dqdx, dqdy, dqdz, mu, matrixCoeff);
        }
        else if (bcType == PHENGLEI::OUTFLOW)
        {
            calPressureOutletBC(grid, faceIndex, bcData, diagMatrixCoeff, bCoeff, q, dqdx, dqdy, dqdz, mu, matrixCoeff);
        }
        else if (bcType == PHENGLEI::PRESSURE_INLET)
        {
            calPressureInletBC(grid, faceIndex, bcData, diagMatrixCoeff, bCoeff, q, dqdx, dqdy, dqdz, mu, matrixCoeff);
        }
        else if (bcType == PHENGLEI::PRESSURE_OUTLET)
        {
            calPressureOutletBC(grid, faceIndex, bcData, diagMatrixCoeff, bCoeff, q, dqdx, dqdy, dqdz, mu, matrixCoeff);
        }
        else if (bcType == PHENGLEI::INTERFACE)
        {
            calInterfaceBC(grid, faceIndex, bcData, diagMatrixCoeff, bCoeff, q, dqdx, dqdy, dqdz, mu, matrixCoeff);
        }
    }
}

void IncomScalarEqCalculator::calFarfieldBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *area = grid->GetFaceArea();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));

    for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        int re = rightCellOfFace[iFacelocal];

        RDouble muLocal = mu[le];

        RDouble dx = xfc[iFacelocal] - xcc[le];
        RDouble dy = yfc[iFacelocal] - ycc[le];
        RDouble dz = zfc[iFacelocal] - zcc[le];
        RDouble area2 = muLocal * 
            area[iFacelocal] / (xfn[iFacelocal] * dx + yfn[iFacelocal] * dy + zfn[iFacelocal] * dz);
        RDouble ax = muLocal * xfn[iFacelocal] * area[iFacelocal] - dx * area2;
        RDouble ay = muLocal * yfn[iFacelocal] * area[iFacelocal] - dy * area2;
        RDouble az = muLocal * zfn[iFacelocal] * area[iFacelocal] - dz * area2;

        diagMatrixCoeff[le] = diagMatrixCoeff[le] + area2 * matrixCoeff;
        bCoeff[le] = bCoeff[le] + area2 * q[re] + dqdx[le] * ax + dqdy[le] * ay + dqdz[le] * az;

        if (FaceFlux[iFacelocal] <0)
        {
            RDouble a12 = MAX(FaceFlux[iFacelocal], 0.0);
            RDouble a21 = a12 - FaceFlux[iFacelocal];
            diagMatrixCoeff[le] = diagMatrixCoeff[le] + a21 * matrixCoeff;
            bCoeff[le] = bCoeff[le] + a21 * q[re];
        }
    }
}

void IncomScalarEqCalculator::calPressureInletBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *area = grid->GetFaceArea();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));

    for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        int re = rightCellOfFace[iFacelocal];

        RDouble muLocal = mu[le];

        RDouble dx = xfc[iFacelocal] - xcc[le];
        RDouble dy = yfc[iFacelocal] - ycc[le];
        RDouble dz = zfc[iFacelocal] - zcc[le];

        RDouble area2 = muLocal * 
            area[iFacelocal] / (xfn[iFacelocal] * dx + yfn[iFacelocal] * dy + zfn[iFacelocal] * dz);
        RDouble ax = muLocal * xfn[iFacelocal] * area[iFacelocal] - dx * area2;
        RDouble ay = muLocal * yfn[iFacelocal] * area[iFacelocal] - dy * area2;
        RDouble az = muLocal * zfn[iFacelocal] * area[iFacelocal] - dz * area2;

        diagMatrixCoeff[le] = diagMatrixCoeff[le] + area2 * matrixCoeff;
        bCoeff[le] = bCoeff[le] + area2 * q[re] + dqdx[le] * ax + dqdy[le] * ay + dqdz[le] * az;

        if (FaceFlux[iFacelocal] < 0)
        {
            RDouble a12 = std::max(FaceFlux[iFacelocal], 0.0);
            RDouble a21 = a12 - FaceFlux[iFacelocal];

            diagMatrixCoeff[le] = diagMatrixCoeff[le] + a21 * matrixCoeff;
            bCoeff[le] = bCoeff[le] + a21 * q[re];
        }
    }
}


void IncomScalarEqCalculator::calPressureOutletBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *area = grid->GetFaceArea();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));

    for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        int re = rightCellOfFace[iFacelocal];

        RDouble muLocal = mu[le];
        RDouble dx = xfc[iFacelocal] - xcc[le];
        RDouble dy = yfc[iFacelocal] - ycc[le];
        RDouble dz = zfc[iFacelocal] - zcc[le];

        RDouble area2 = muLocal * 
            area[iFacelocal] / (xfn[iFacelocal] * dx + yfn[iFacelocal] * dy + zfn[iFacelocal] * dz);
        RDouble ax = muLocal * xfn[iFacelocal] * area[iFacelocal] - dx * area2;
        RDouble ay = muLocal * yfn[iFacelocal] * area[iFacelocal] - dy * area2;
        RDouble az = muLocal * zfn[iFacelocal] * area[iFacelocal] - dz * area2;


        diagMatrixCoeff[le] = diagMatrixCoeff[le] + area2 * matrixCoeff;
        bCoeff[le] = bCoeff[le] + area2 * q[re] + dqdx[le] * ax + dqdy[le] * ay + dqdz[le] * az;

        if (FaceFlux[iFacelocal] < 0)
        {
            RDouble a12 = MAX(FaceFlux[iFacelocal], 0.0);
            RDouble a21 = a12 - FaceFlux[iFacelocal];

            diagMatrixCoeff[le] = diagMatrixCoeff[le] + a21 * matrixCoeff;
            bCoeff[le] = bCoeff[le] + a21 * q[re];
        }
    }
}

void IncomScalarEqCalculator::calInterfaceBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *area = grid->GetFaceArea();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();
    RDouble *faceWeightOfLeftCell = reinterpret_cast<RDouble *>(grid->GetDataPtr("faceWeightOfLeftCell"));
    RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("upperMatrixCoeff"));
    RDouble *lowerMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("lowerMatrixCoeff"));
    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));

    for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        int re = rightCellOfFace[iFacelocal];

        RDouble muLocal = faceWeightOfLeftCell[iFace] * mu[le] + (1.0 - faceWeightOfLeftCell[iFace]) * mu[re];
        RDouble dx = xcc[re] - xcc[le];
        RDouble dy = ycc[re] - ycc[le];
        RDouble dz = zcc[re] - zcc[le];
        RDouble area2 = muLocal * 
            area[iFacelocal] / (xfn[iFacelocal] * dx + yfn[iFacelocal] * dy + zfn[iFacelocal] * dz);
        RDouble ax = muLocal * xfn[iFacelocal] * area[iFacelocal] - dx * area2;
        RDouble ay = muLocal * yfn[iFacelocal] * area[iFacelocal] - dy * area2;
        RDouble az = muLocal * zfn[iFacelocal] * area[iFacelocal] - dz * area2;


        diagMatrixCoeff[le] = diagMatrixCoeff[le] + area2 * matrixCoeff;
        upperMatrixCoeff[iFacelocal] = upperMatrixCoeff[iFacelocal] - area2 * matrixCoeff;

        RDouble a12 = MAX(FaceFlux[iFacelocal], 0.0);
        RDouble a21 = a12 - FaceFlux[iFacelocal];

        diagMatrixCoeff[le] = diagMatrixCoeff[le] + a21 * matrixCoeff;
        upperMatrixCoeff[iFacelocal] = upperMatrixCoeff[iFacelocal] - a21 * matrixCoeff;

        RDouble sLR = ax * (faceWeightOfLeftCell[iFace] * dqdx[le] + (1.0 - faceWeightOfLeftCell[iFace]) * dqdx[re]) +
            ay * (faceWeightOfLeftCell[iFace] * dqdy[le] + (1.0 - faceWeightOfLeftCell[iFace]) * dqdy[re]) +
            az * (faceWeightOfLeftCell[iFace] * dqdz[le] + (1.0 - faceWeightOfLeftCell[iFace]) * dqdz[re]);

        bCoeff[le] += sLR;

    }
}

void IncomScalarEqCalculator::calSymmetryBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *area = grid->GetFaceArea();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();
    int nTotalCell = grid->GetNTotalCell(); 

    for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        int re = rightCellOfFace[iFacelocal];
        RDouble muLocal = mu[le];

        RDouble dx = xfc[iFacelocal] - xcc[le];
        RDouble dy = yfc[iFacelocal] - ycc[le];
        RDouble dz = zfc[iFacelocal] - zcc[le];

        RDouble area2 = muLocal * 
            area[iFacelocal] / (xfn[iFacelocal] * dx + yfn[iFacelocal] * dy + zfn[iFacelocal] * dz);
        RDouble ax = muLocal * xfn[iFacelocal] * area[iFacelocal] - dx * area2;
        RDouble ay = muLocal * yfn[iFacelocal] * area[iFacelocal] - dy * area2;
        RDouble az = muLocal * zfn[iFacelocal] * area[iFacelocal] - dz * area2;

        diagMatrixCoeff[le] = diagMatrixCoeff[le] + area2 * matrixCoeff;
        bCoeff[le] = bCoeff[le] + area2 * q[re] + dqdx[le] * ax + dqdy[le] * ay + dqdz[le] * az;
    }
}

void IncomScalarEqCalculator::calVelocityInletBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *area    = grid->GetFaceArea();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));

    for (size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        int re = rightCellOfFace[iFacelocal];

        RDouble dx = xfc[iFacelocal] - xcc[le];
        RDouble dy = yfc[iFacelocal] - ycc[le];
        RDouble dz = zfc[iFacelocal] - zcc[le];

        RDouble area2 = mu[le] * 
            area[iFacelocal] / (xfn[iFacelocal] * dx + yfn[iFacelocal] * dy + zfn[iFacelocal] * dz);
        RDouble ax = mu[le] * xfn[iFacelocal] * area[iFacelocal] - dx * area2;
        RDouble ay = mu[le] * yfn[iFacelocal] * area[iFacelocal] - dy * area2;
        RDouble az = mu[le] * zfn[iFacelocal] * area[iFacelocal] - dz * area2;

        diagMatrixCoeff[le] = diagMatrixCoeff[le] + area2 * matrixCoeff;
        bCoeff[le] = bCoeff[le] + area2 * q[re] + dqdx[le] * ax + dqdy[le] * ay + dqdz[le] * az;

        if (FaceFlux[iFacelocal] < 0)
        {
            RDouble a12 = MAX(FaceFlux[iFacelocal], 0.0);
            RDouble a21 = a12 - FaceFlux[iFacelocal];
            diagMatrixCoeff[le] = diagMatrixCoeff[le] + a21 * matrixCoeff;
            bCoeff[le] = bCoeff[le] + a21 * q[re];
        }
    }
}


void IncomScalarEqCalculator::calWallBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *area = grid->GetFaceArea();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        int re = rightCellOfFace[iFacelocal];

        RDouble muLocal = mu[re];

        RDouble dx = xfc[iFacelocal] - xcc[le];
        RDouble dy = yfc[iFacelocal] - ycc[le];
        RDouble dz = zfc[iFacelocal] - zcc[le];

        RDouble area2 = muLocal * 
            area[iFacelocal] / (xfn[iFacelocal] * dx + yfn[iFacelocal] * dy + zfn[iFacelocal] * dz);
        RDouble ax = muLocal * xfn[iFacelocal] * area[iFacelocal] - dx * area2;
        RDouble ay = muLocal * yfn[iFacelocal] * area[iFacelocal] - dy * area2;
        RDouble az = muLocal * zfn[iFacelocal] * area[iFacelocal] - dz * area2;

        diagMatrixCoeff[le] = diagMatrixCoeff[le] + area2 * matrixCoeff;
        bCoeff[le] = bCoeff[le] + area2 * q[re] + dqdx[le] * ax + dqdy[le] * ay + dqdz[le] * az;
    }
}


void IncomScalarEqCalculator::calcInvisSourceTerm(Grid *grid, int solverIndex) 
{
    string ConvCalcMethod = GlobalDataBase::GetStrParaFromDB("ConvCalcMethod");
    if (ConvCalcMethod == "UPWIND")
    {
    } 
    else if (ConvCalcMethod == "CDS")
    {
        SecondCentralSourceTerm(grid, solverIndex);
    }
    else if (ConvCalcMethod == "QUICK")
    {
        QUICKSourceTerm(grid, solverIndex);
    }
    else if (ConvCalcMethod == "SUDS")
    {
        SecondUpWindSourceTerm(grid, solverIndex);
    }
}


void IncomScalarEqCalculator::FirstUpwindInvisTerm(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("upperMatrixCoeff"));
    RDouble *lowerMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("lowerMatrixCoeff"));
    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        RDouble a12 = MAX(FaceFlux[iFace], 0.0);
        RDouble a21 = a12 - FaceFlux[iFace];

        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];

        upperMatrixCoeff[iFace] -= a21;
        diagMatrixCoeff[le]  += a21;

        lowerMatrixCoeff[iFace] -= a12;
        diagMatrixCoeff[re]  += a12;
    }
}

void IncomScalarEqCalculator::QUICKSourceTerm(Grid *gridIn, int solverIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    
    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
    RDouble *bCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("bCoeff"));

    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    string varName = varNameIncom[solverIndex];
    RDouble *q = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
    RDouble *dqdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dx"));
    RDouble *dqdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dy"));
    RDouble *dqdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dz"));

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];
        int ic = le;
        int id = re;
        if (FaceFlux[iFace] < 0)
        {
            ic = re;
            id = le;
        }
        RDouble dx = xcc[id] - xcc[ic];
        RDouble dy = ycc[id] - ycc[ic];
        RDouble dz = zcc[id] - zcc[ic];
        RDouble cq = q[id] - 2 * (dqdx[ic] * dx + dqdy[ic] * dy + dqdz[ic] * dz);
        RDouble cq_r = MAX((q[ic] - cq) / (q[id] - q[ic] + SMALL),0);
        RDouble b = MIN(2 * cq_r, 2);
        RDouble a12 = fabs(FaceFlux[iFace]) * (0.5 * MIN((3 + cq_r) / 4, b) * (q[ic] - q[id]));
        bCoeff[ic] += a12;
        bCoeff[id] -= a12;
    }
}


void IncomScalarEqCalculator::SecondUpWindSourceTerm(Grid *gridIn, int solverIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
    RDouble *bCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("bCoeff"));

    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    string varName = varNameIncom[solverIndex];
    RDouble *q = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
    RDouble *dqdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dx"));
    RDouble *dqdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dy"));
    RDouble *dqdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dz"));

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];
        int ic = le;
        int id = re;
        if (FaceFlux[iFace] < 0)
        {
            ic = re;
            id = le;
        }
        RDouble dx = xcc[id] - xcc[ic];
        RDouble dy = ycc[id] - ycc[ic];
        RDouble dz = zcc[id] - zcc[ic];
        RDouble cu = q[id] - 2 * (dqdx[ic] * dx + dqdy[ic] * dy + dqdz[ic] * dz);
        RDouble cu_r = MAX((q[ic] - cu) / (q[id] - q[ic] + SMALL),0);

        RDouble a12 = fabs(FaceFlux[iFace]) * (0.5 * MIN(cu_r, 2)) * (q[ic] - q[id]);
        bCoeff[ic] += a12;
        bCoeff[id] -= a12;
    }

}

void IncomScalarEqCalculator::SecondCentralSourceTerm(Grid *gridIn, int solverIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
    RDouble *bCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("bCoeff"));

    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    string varName = varNameIncom[solverIndex];
    RDouble *q = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
    RDouble *dqdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dx"));
    RDouble *dqdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dy"));
    RDouble *dqdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dz"));

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];
        int ic = le;
        int id = re;
        if (FaceFlux[iFace] < 0)
        {
            ic = re;
            id = le;
        }
        RDouble dx = xcc[id] - xcc[ic];
        RDouble dy = ycc[id] - ycc[ic];
        RDouble dz = zcc[id] - zcc[ic];
        RDouble cq = q[id] - 2 * (dqdx[ic] * dx + dqdy[ic] * dy + dqdz[ic] * dz);
        RDouble cq_r = MAX((q[ic] - cq) / (q[id] - q[ic] + SMALL),0);
        RDouble b = MIN(2 * cq_r, 2);
        RDouble a12 = fabs(FaceFlux[iFace]) * (0.5 * MIN((3 + cq_r) / 4, b) * (q[ic] - q[id]));
        bCoeff[ic] += a12;
        bCoeff[id] -= a12;
    }
}

void IncomScalarEqCalculator::VisMatrixTerm(Grid *gridIn, int solverIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *area = grid->GetFaceArea();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *faceWeightOfLeftCell = reinterpret_cast<RDouble *>(grid->GetDataPtr("faceWeightOfLeftCell"));

    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("upperMatrixCoeff"));
    RDouble *lowerMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("lowerMatrixCoeff"));

    RDouble *mu;
    RDouble *tempmu;
    using namespace IDX;
    RDouble **DiffusionCoeff = reinterpret_cast <RDouble **> (GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    
    mu = DiffusionCoeff[solverIndex];
    if (solverIndex == S_ITEMP) 
    {
        int nTotalCell = grid->GetNTotalCell();
        RDouble *cp = reinterpret_cast<RDouble *>(grid->GetDataPtr("cp"));
        tempmu = new RDouble [nTotalCell];
        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            tempmu[iCell] = DiffusionCoeff[solverIndex][iCell] / cp[iCell];
        }
        mu = tempmu;
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        RDouble ax = area[iFace] * xfn[iFace];
        RDouble ay = area[iFace] * yfn[iFace];
        RDouble az = area[iFace] * zfn[iFace];

        RDouble ex = xcc[re] - xcc[le];
        RDouble ey = ycc[re] - ycc[le];
        RDouble ez = zcc[re] - zcc[le];

        RDouble diffCoeF = faceWeightOfLeftCell[iFace] * mu[le] + (1.0 - faceWeightOfLeftCell[iFace]) * mu[re];

        RDouble aLR = diffCoeF * (ax * ax + ay * ay + az * az) / (ax * ex + ay * ey + az * ez);

        upperMatrixCoeff[iFace] -= aLR;
        diagMatrixCoeff[le]  += aLR;

        lowerMatrixCoeff[iFace] -= aLR;
        diagMatrixCoeff[re]  += aLR;
    }

    if (solverIndex == S_ITEMP)
    {
        delete [] tempmu;
    }
}


void IncomScalarEqCalculator::VisSourceTerm(Grid *gridIn, int solverIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *area = grid->GetFaceArea();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *faceWeightOfLeftCell = reinterpret_cast<RDouble *>(grid->GetDataPtr("faceWeightOfLeftCell"));

    RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));
    string *varNameIncom = reinterpret_cast<string *>(GlobalDataBase::GetDataPtr("varNameIncom"));
    string varName = varNameIncom[solverIndex];
    RDouble *q = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
    RDouble *dqdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dx"));
    RDouble *dqdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dy"));
    RDouble *dqdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dz"));

    RDouble **DiffusionCoeff = reinterpret_cast<RDouble **>(GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    RDouble *mu = DiffusionCoeff[solverIndex];

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        RDouble ax = area[iFace] * xfn[iFace];
        RDouble ay = area[iFace] * yfn[iFace];
        RDouble az = area[iFace] * zfn[iFace];

        RDouble ex = xcc[re] - xcc[le];
        RDouble ey = ycc[re] - ycc[le];
        RDouble ez = zcc[re] - zcc[le];

        RDouble diffCoeF = faceWeightOfLeftCell[iFace] * mu[le] + (1.0 - faceWeightOfLeftCell[iFace]) * mu[re];

        RDouble aLR = diffCoeF * (ax * ax + ay * ay + az * az) / (ax * ex + ay * ey + az * ez);

        ax = diffCoeF * ax - ex * aLR;
        ay = diffCoeF * ay - ey * aLR;
        az = diffCoeF * az - ez * aLR;

        RDouble sLR = ax * (faceWeightOfLeftCell[iFace] * dqdx[le] + (1.0 - faceWeightOfLeftCell[iFace]) * dqdx[re]) +
                      ay * (faceWeightOfLeftCell[iFace] * dqdy[le] + (1.0 - faceWeightOfLeftCell[iFace]) * dqdy[re]) +
                      az * (faceWeightOfLeftCell[iFace] * dqdz[le] + (1.0 - faceWeightOfLeftCell[iFace]) * dqdz[re]);

        bCoeff[le] += sLR;
        bCoeff[re] -= sLR;
    }
}

void IncomScalarEqCalculator::calcTransMatrixTerm(Grid *grid)
{
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady == 1)
    {
        string TranCalcMethod = GlobalDataBase::GetStrParaFromDB("TranCalcMethod");
        if (TranCalcMethod == "IMPLICIT_EULER")
        {
            TransMatrixTerm_1st(grid);
        }
        else if (TranCalcMethod == "IMPLICIT_2ND_ORDER")
        {
            TransMatrixTerm_2nd(grid);
        }
    }
}

void IncomScalarEqCalculator::TransMatrixTerm_1st(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    RDouble dt = GlobalDataBase::GetDoubleParaFromDB("dt");
    RDouble *cvOld = grid->GetCellVolume();
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        diagMatrixCoeff[iCell] += rho[iCell] * cvOld[iCell] / dt;
    }
}

void IncomScalarEqCalculator::TransMatrixTerm_2nd(Grid *gridIn) 
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble dt = GlobalDataBase::GetDoubleParaFromDB("dt");

    int timeStepNow;
    GlobalDataBase::GetData("timeStepNow", &timeStepNow, PHINT, 1);

    RDouble *cvOld = grid->GetCellVolume();
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));

    if (timeStepNow == 0)
    {
        for (int iCell = 0; iCell < nTotalCell; ++iCell)
        {
            diagMatrixCoeff[iCell] += rho[iCell] * cvOld[iCell] / dt;
        }
    }
    else
    {
        for (int iCell = 0; iCell < nTotalCell; ++iCell)
        {
            diagMatrixCoeff[iCell] += 3 * rho[iCell] * cvOld[iCell] / (2 * dt);
        }
    }
}

void IncomScalarEqCalculator::calcTransSourceTerm(Grid *grid, int iVariable)
{
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady == 1)
    {
        string TranCalcMethod = GlobalDataBase::GetStrParaFromDB("TranCalcMethod");
        if (TranCalcMethod == "IMPLICIT_EULER")
        {
            TransSourceTerm_1st(grid, iVariable);
        }
        else if (TranCalcMethod == "IMPLICIT_2ND_ORDER")
        {
            TransSourceTerm_2nd(grid, iVariable);
        }
    }
}

void IncomScalarEqCalculator::TransSourceTerm_1st(Grid *gridIn, int iVariable)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    RDouble *cvOld = grid->GetCellVolume();
    RDouble dt = GlobalDataBase::GetDoubleParaFromDB("dt");
    RDouble *rhoOld = reinterpret_cast<RDouble *>(grid->GetDataPtr("rhoOld"));
    RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));
    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    string varName = varNameIncom[iVariable];
    RDouble *qOld = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName + "Old"));

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        bCoeff[iCell]  += rhoOld[iCell] * cvOld[iCell] / dt * qOld[iCell];
    }
}

void IncomScalarEqCalculator::TransSourceTerm_2nd(Grid *gridIn, int iVariable)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    RDouble *cvOld = grid->GetCellVolume();
    RDouble dt = GlobalDataBase::GetDoubleParaFromDB("dt");

    int timeStepNow;
    GlobalDataBase::GetData("timeStepNow", &timeStepNow, PHINT, 1);
    RDouble *rhoOld = reinterpret_cast<RDouble *>(grid->GetDataPtr("rhoOld"));
    RDouble *rhoOldOld = reinterpret_cast<RDouble *>(grid->GetDataPtr("rhoOldOld"));
    RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));
    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    string varName = varNameIncom[iVariable];
    RDouble *qOld = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName + "Old"));
    RDouble *qOldOld = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName + "OldOld"));

    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        RDouble tmpOld = 2 * rhoOld[iCell] * cvOld[iCell] / dt;
        RDouble tmpOldOld = -rhoOldOld[iCell] * cvOld[iCell] / (2 * dt);
        bCoeff[iCell] += tmpOld * qOld[iCell] + tmpOldOld * qOldOld[iCell];
    }
}

void IncomScalarEqCalculator::relaxMatrixTerm(Grid *grid, int iVariable)
{
    int nTotalCell = grid->GetNTotalCell();
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *soluRelaxCoeff = reinterpret_cast <RDouble *>(GlobalDataBase::GetDataPtr("soluRelaxCoeff"));

    RDouble RelaxCoeff = soluRelaxCoeff[iVariable];

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        diagMatrixCoeff[iCell] =  diagMatrixCoeff[iCell] / RelaxCoeff;
    }
}

void IncomScalarEqCalculator::relaxSourceTerm(Grid *grid, int iVariable)
{
    int nTotalCell = grid->GetNTotalCell();
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));

    string *varNameIncom = reinterpret_cast <string *>(GlobalDataBase::GetDataPtr("varNameIncom"));
    string varName = varNameIncom[iVariable];
    RDouble *q = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));

    RDouble *soluRelaxCoeff = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("soluRelaxCoeff"));
    RDouble RelaxCoeff = soluRelaxCoeff[iVariable];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        bCoeff[iCell] = bCoeff[iCell] + (1 - RelaxCoeff) * diagMatrixCoeff[iCell] * q[iCell];
    }
}

void IncomScalarEqCalculator::ImplicitEuler_AllocateMemory(Grid *grid, vector<string> &phiName)
{
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    for (vector<string>::iterator iter = phiName.begin(); iter != phiName.end(); ++ iter)
    {
        RDouble *phi = NewPointer<RDouble>(nTotal);
        grid->UpdateDataPtr((*iter) + "Old", phi);
    }
}

void IncomScalarEqCalculator::ImplicitEuler_ReInitTimeVar(Grid *grid, vector<string> &phiName)
{
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    for (vector<string>::iterator iter = phiName.begin(); iter != phiName.end(); ++iter)
    {
        RDouble *initPhi = reinterpret_cast<RDouble *>(grid->GetDataPtr(*iter));
        RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr((*iter) + "Old"));

        for (int iCell = 0; iCell < nTotal; ++iCell)
        {
            phi[iCell] = initPhi[iCell];
        }
    }
}

void IncomScalarEqCalculator::ImplicitEuler_SaveOldTimeValue(Grid *gridIn, vector<string> &phiName)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int nTotalCell = grid->GetNTotalCell();
    vector<RDouble *> phiNow;
    vector<RDouble *> phiOld;
    for (vector<string>::iterator iter = phiName.begin(); iter != phiName.end(); ++ iter)
    {
        RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr((*iter)));
        RDouble *phiOldPtr = reinterpret_cast<RDouble *>(grid->GetDataPtr((*iter) + "Old"));
        phiNow.push_back(phi);
        phiOld.push_back(phiOldPtr);
    }

    for (int index = 0; index < phiName.size(); ++ index)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            phiOld[index][iCell] = phiNow[index][iCell];

        }
    }

    vector<RDouble *>().swap(phiNow);
    vector<RDouble *>().swap(phiOld);
}

void IncomScalarEqCalculator::ImplicitEuler_InitRestartTimeVar(Grid *grid, vector<string> &nameLists)
{
    int nTotalCell = grid->GetNTotalCell();
    for (vector<string>::iterator iter = nameLists.begin(); iter != nameLists.end(); ++iter)
    {
        RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(*iter));
        RDouble *phiOldPtr = reinterpret_cast<RDouble *>(grid->GetDataPtr((*iter) + "Old"));

        for (int iCell = 0; iCell < nTotalCell; ++iCell)
        {
            phiOldPtr[iCell] = phi[iCell];
        }
    }
}

void IncomScalarEqCalculator::ImplicitEuler_ReadRestartFile(Grid * grid, vector<string>&nameLists, hid_t grploc)
{
    for (vector<string>::iterator iter = nameLists.begin(); iter != nameLists.end(); ++iter)
    {
        RDouble *phiOldPtr = reinterpret_cast<RDouble *>(grid->GetDataPtr((*iter) + "Old"));
        if (phiOldPtr != nullptr)
        {
            ReadData(grploc, &phiOldPtr[0], (*iter));
        }
    }
}



void IncomScalarEqCalculator::Implicit2ndOrder_AllocateMemory(Grid *grid, vector<string> &phiName)
{
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;
    RDouble *phiOld;
    RDouble *phiOldOld;

    for (vector<string>::iterator iter = phiName.begin(); iter != phiName.end(); ++ iter)
    {
        phiOld = NewPointer<RDouble>(nTotal);
        grid->UpdateDataPtr((*iter) + "Old", phiOld);

        phiOldOld = NewPointer<RDouble>(nTotal);
        grid->UpdateDataPtr((*iter) + "OldOld", phiOldOld);
    }
}

void IncomScalarEqCalculator::Implicit2ndOrder_ReInitTimeVar(Grid *grid, vector<string> &phiName)
{
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    for (vector<string>::iterator iter = phiName.begin(); iter != phiName.end(); ++iter)
    {
        RDouble *initPhi = reinterpret_cast<RDouble *>(grid->GetDataPtr(*iter));
        RDouble *phiOld = reinterpret_cast<RDouble *>(grid->GetDataPtr((*iter) + "Old"));
        RDouble *phiOldOld = reinterpret_cast<RDouble *>(grid->GetDataPtr((*iter) + "OldOld"));

        for (int iCell = 0; iCell < nTotal; ++iCell)
        {
            phiOld[iCell] = initPhi[iCell];
            phiOldOld[iCell] = initPhi[iCell];
        }
    }
}

void IncomScalarEqCalculator::Implicit2ndOrder_SaveOldTimeValue(Grid *gridIn, vector<string> &phiName)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int nTotalCell = grid->GetNTotalCell();
    vector<RDouble *> phiNow;
    vector<RDouble *> phiOld;
    vector<RDouble *> phiOldOld;

    for (vector<string>::iterator iter = phiName.begin(); iter != phiName.end(); ++ iter)
    {
        RDouble *phiNowPtr = reinterpret_cast<RDouble *>(grid->GetDataPtr((*iter)));
        RDouble *phiOldPtr = reinterpret_cast<RDouble *>(grid->GetDataPtr((*iter) + "Old"));
        RDouble *phiOldOldPtr = reinterpret_cast<RDouble *>(grid->GetDataPtr((*iter) + "OldOld"));
        phiNow.push_back(phiNowPtr);
        phiOld.push_back(phiOldPtr);
        phiOldOld.push_back(phiOldOldPtr);
    }

    for (int index = 0; index < phiName.size(); ++ index)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            phiOldOld[index][iCell] = phiOld[index][iCell];
            phiOld[index][iCell]    = phiNow[index][iCell];
        }
    }

    vector<RDouble *>().swap(phiNow);
    vector<RDouble *>().swap(phiOld);
    vector<RDouble *>().swap(phiOldOld);
}

void IncomScalarEqCalculator::Implicit2ndOrder_InitRestartTimeVar(Grid *grid, vector<string> &nameLists)
{
    int nTotalCell = grid->GetNTotalCell();
    for (vector<string>::iterator iter = nameLists.begin(); iter != nameLists.end(); ++iter)
    {
        RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(*iter));
        RDouble *phiOldPtr = reinterpret_cast<RDouble *>(grid->GetDataPtr((*iter) + "Old"));
        RDouble *phiOldOldPtr = reinterpret_cast<RDouble *>(grid->GetDataPtr((*iter) + "OldOld"));

        for (int iCell = 0; iCell < nTotalCell; ++iCell)
        {
            phiOldPtr[iCell] = phi[iCell];
            phiOldOldPtr[iCell] = phi[iCell];
        }
    }
}

void IncomScalarEqCalculator::Implicit2ndOrder_ReadRestartFile(Grid *grid, vector<string> &nameLists, hid_t grploc)
{
    for (vector<string>::iterator iter = nameLists.begin(); iter != nameLists.end(); ++iter)
    {
        RDouble *phiOldPtr = reinterpret_cast<RDouble *>(grid->GetDataPtr((*iter) + "Old"));
        if (phiOldPtr != nullptr)
        {
            ReadData(grploc, &phiOldPtr[0], (*iter));
        }

        RDouble *phiOldOldPtr = reinterpret_cast<RDouble *>(grid->GetDataPtr((*iter) + "OldOld"));
        if (phiOldOldPtr != nullptr)
        {
            ReadData(grploc, &phiOldOldPtr[0], (*iter));
        }
    }
}


}
