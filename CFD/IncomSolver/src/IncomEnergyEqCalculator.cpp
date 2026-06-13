#include "IncomEnergyEqCalculator.h"


namespace PHSPACE
{
IncomEnergyEqCalculator::IncomEnergyEqCalculator():IncomScalarEqCalculator()
{

}

IncomEnergyEqCalculator::~IncomEnergyEqCalculator()
{
}

void IncomEnergyEqCalculator::CalcGrad(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    string varName = varNameIncom[solverIndex];

    GradientCalculation(grid, "T", "d" + varName + "dx", "d" + varName + "dy", "d" + varName + "dz");
}

void IncomEnergyEqCalculator::calWallBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();

    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    string *varNameIncom = reinterpret_cast<string *>(GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    string varName = varNameIncom[solverIndex];

    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr("T"));
    RDouble *cp = reinterpret_cast<RDouble *>(grid->GetDataPtr("cp"));
    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dx"));
    RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dy"));
    RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dz"));
    RDouble **DiffusionCoeff = reinterpret_cast<RDouble **>(GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    
    for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
    {
        int iFacelocal = (*faceIndex)[iFace]; //! The boundary element serial number is arranged in front of the element face.
        int le = leftCellOfFace[iFacelocal];
        int re = grid->GetRightCellOfFace()[iFacelocal];

        RDouble dx = xfc[iFacelocal] - xcc[le];
        RDouble dy = yfc[iFacelocal] - ycc[le];
        RDouble dz = zfc[iFacelocal] - zcc[le];

        RDouble area2 = DiffusionCoeff[solverIndex][re] * area[iFacelocal] / (xfn[iFacelocal] * dx + yfn[iFacelocal] * dy + zfn[iFacelocal] * dz);
        RDouble ax = DiffusionCoeff[solverIndex][re] * xfn[iFacelocal] * area[iFacelocal] - dx * area2;
        RDouble ay = DiffusionCoeff[solverIndex][re] * yfn[iFacelocal] * area[iFacelocal] - dy * area2;
        RDouble az = DiffusionCoeff[solverIndex][re] * zfn[iFacelocal] * area[iFacelocal] - dz * area2;

        diagMatrixCoeff[le] += area2 / cp[le] * matrixCoeff;
        bCoeff[le] += area2 * phi[re] + dphidx[le] * ax + dphidy[le] * ay + dphidz[le] * az;
    }
}

void IncomEnergyEqCalculator::calVelocityInletBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell =grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();

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

    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    string varName = varNameIncom[solverIndex];

    RDouble *Enthalpy = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
    RDouble *phi= reinterpret_cast<RDouble *>(grid->GetDataPtr("T"));
    RDouble *cp = reinterpret_cast<RDouble *>(grid->GetDataPtr("cp"));
    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dx"));
    RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dy"));
    RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dz"));
    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
    RDouble **DiffusionCoeff = reinterpret_cast <RDouble **> (GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    
    for (int iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        int iFacelocal = (*faceIndex)[iFace]; //! The boundary element serial number is arranged in front of the element face.
        int le = leftCellOfFace[iFacelocal];
        int re = grid->GetRightCellOfFace()[iFacelocal];

        RDouble dx = xfc[iFacelocal] - xcc[le];
        RDouble dy = yfc[iFacelocal] - ycc[le];
        RDouble dz = zfc[iFacelocal] - zcc[le];

        RDouble area2 = DiffusionCoeff[solverIndex][re] * 
            area[iFacelocal] / (xfn[iFacelocal] * dx + yfn[iFacelocal] * dy + zfn[iFacelocal] * dz);
        RDouble ax = DiffusionCoeff[solverIndex][re] * xfn[iFacelocal] * area[iFacelocal] - dx * area2;
        RDouble ay = DiffusionCoeff[solverIndex][re] * yfn[iFacelocal] * area[iFacelocal] - dy * area2;
        RDouble az = DiffusionCoeff[solverIndex][re] * zfn[iFacelocal] * area[iFacelocal] - dz * area2;

        diagMatrixCoeff[le] += area2 / cp[re] * matrixCoeff;
        bCoeff[le] += area2 * phi[re] + dphidx[le] * ax + dphidy[le] * ay + dphidz[le] * az;


        if (FaceFlux[iFacelocal] < 0)
        {
            diagMatrixCoeff[le] -= FaceFlux[iFacelocal] * matrixCoeff;
            bCoeff[le] -= FaceFlux[iFacelocal] * Enthalpy[re];
        }
    }
}

void IncomEnergyEqCalculator::calPressureOutletBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();

    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    string varName = varNameIncom[solverIndex];

    RDouble *Enthalpy = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr("T"));
    
    RDouble *cp = reinterpret_cast<RDouble *>(grid->GetDataPtr("cp"));
    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dx"));
    RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dy"));
    RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dz"));
    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
    RDouble **DiffusionCoeff = reinterpret_cast<RDouble **>(GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    
    for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
    {
        int iFacelocal = (*faceIndex)[iFace]; //! The boundary element serial number is arranged in front of the element face.
        int le = leftCellOfFace[iFacelocal];
        int re = grid->GetRightCellOfFace()[iFacelocal];

        RDouble dx = xfc[iFacelocal] - xcc[le];
        RDouble dy = yfc[iFacelocal] - ycc[le];
        RDouble dz = zfc[iFacelocal] - zcc[le];

        RDouble area2 = DiffusionCoeff[solverIndex][re] * 
            area[iFacelocal] / (xfn[iFacelocal] * dx + yfn[iFacelocal] * dy + zfn[iFacelocal] * dz);
        RDouble ax = DiffusionCoeff[solverIndex][re] * xfn[iFacelocal] * area[iFacelocal] - dx * area2;
        RDouble ay = DiffusionCoeff[solverIndex][re] * yfn[iFacelocal] * area[iFacelocal] - dy * area2;
        RDouble az = DiffusionCoeff[solverIndex][re] * zfn[iFacelocal] * area[iFacelocal] - dz * area2;

        diagMatrixCoeff[le] += area2 / cp[re] * matrixCoeff;
        bCoeff[le] += area2 * phi[re] + dphidx[le] * ax + dphidy[le] * ay + dphidz[le] * az;

        if (FaceFlux[iFacelocal] < 0)
        {
            diagMatrixCoeff[le] -= FaceFlux[iFacelocal] * matrixCoeff;
            bCoeff[le] -= FaceFlux[iFacelocal] * Enthalpy[re];
        }
    }
}

void IncomEnergyEqCalculator::calPressureInletBC(Grid *grid, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff)
{
    calVelocityInletBC(grid, faceIndex, bcData, diagMatrixCoeff, bCoeff, q, dqdx, dqdy, dqdz, mu, matrixCoeff);
}

void IncomEnergyEqCalculator::calFarfieldBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();

    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    string varName = varNameIncom[solverIndex];

    RDouble *Enthalpy = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr("T"));
    
    RDouble *cp = reinterpret_cast<RDouble *>(grid->GetDataPtr("cp"));
    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dx"));
    RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dy"));
    RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dz"));
    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
    RDouble **DiffusionCoeff = reinterpret_cast<RDouble **>(GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    
    for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
    {
        int iFacelocal = (*faceIndex)[iFace]; //! The boundary element serial number is arranged in front of the element face.
        int le = leftCellOfFace[iFacelocal];
        int re = grid->GetRightCellOfFace()[iFacelocal];

        RDouble dx = xfc[iFacelocal] - xcc[le];
        RDouble dy = yfc[iFacelocal] - ycc[le];
        RDouble dz = zfc[iFacelocal] - zcc[le];

        RDouble area2 = DiffusionCoeff[solverIndex][re] * 
            area[iFacelocal] / (xfn[iFacelocal] * dx + yfn[iFacelocal] * dy + zfn[iFacelocal] * dz);

        RDouble ax = DiffusionCoeff[solverIndex][re] * xfn[iFacelocal] * area[iFacelocal] - dx * area2;
        RDouble ay = DiffusionCoeff[solverIndex][re] * yfn[iFacelocal] * area[iFacelocal] - dy * area2;
        RDouble az = DiffusionCoeff[solverIndex][re] * zfn[iFacelocal] * area[iFacelocal] - dz * area2;

        diagMatrixCoeff[le] += area2 / cp[re] * matrixCoeff;
        bCoeff[le] += area2 * phi[re] + dphidx[le] * ax + dphidy[le] * ay + dphidz[le] * az;

        if (FaceFlux[iFacelocal] < 0)
        {
            diagMatrixCoeff[le] -= FaceFlux[iFacelocal] * matrixCoeff;
            bCoeff[le] -= FaceFlux[iFacelocal] * Enthalpy[re];
        }
    }
}

void IncomEnergyEqCalculator::calInterfaceBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
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
    RDouble *cp = reinterpret_cast<RDouble *>(grid->GetDataPtr("cp"));

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


        diagMatrixCoeff[le] = diagMatrixCoeff[le] + area2 / cp[re] * matrixCoeff;
        upperMatrixCoeff[iFacelocal] = upperMatrixCoeff[iFacelocal] - area2 / cp[re] * matrixCoeff;

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

void IncomEnergyEqCalculator::CalcOtherbCoeff(Grid *gridIn, int iEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nInteriorFace = nTotalFace - nBoundFace;
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *area = grid->GetFaceArea();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *mu= reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
    RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));
    RDouble *k     = reinterpret_cast<RDouble *>(grid->GetDataPtr("k"));
    RDouble *dudx  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdx"));
    RDouble *dudy  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdy"));
    RDouble *dudz  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdz"));
    RDouble *dvdx  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdx"));
    RDouble *dvdy  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdy"));
    RDouble *dvdz  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdz"));
    RDouble *dwdx  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdx"));
    RDouble *dwdy  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdy"));
    RDouble *dwdz  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdz"));
    RDouble *faceWeightOfLeftCell  = reinterpret_cast<RDouble *>(grid->GetDataPtr("faceWeightOfLeftCell"));
    RDouble **tau1 = NewPointer2<RDouble> (3, 3);
    RDouble **tau2 = NewPointer2<RDouble> (3, 3);
    RDouble *tauu1 = new double [3];
    RDouble *tauu2 = new double [3];
    RDouble *sf = new double [3];

    RDouble *vol = grid->GetCellVolume();
    RDouble *dpdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdx"));
    RDouble *dpdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdy"));
    RDouble *dpdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdz"));
    RDouble *p = reinterpret_cast<RDouble *>(grid->GetDataPtr("P"));
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    UnstructBCSet **bcr = grid->GetBCRecord();
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcr[iFace]->GetKey();
        if (bcType != PHENGLEI::INTERFACE)
        {
            continue;
        }
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        RDouble muf = faceWeightOfLeftCell[iFace] * mu[le] + (1.0 - faceWeightOfLeftCell[iFace]) * mu[re];
        tau1[0][0] = 2 * dudx[le];
        tau1[1][1] = 2 * dvdy[le];
        tau1[2][2] = 2 * dwdz[le];
        tau1[0][1] = dudy[le] + dvdx[le];
        tau1[1][0] = tau1[0][1];
        tau1[0][2] = dudz[le] + dwdx[le];
        tau1[2][0] = tau1[0][2];
        tau1[1][2] = dvdz[le] + dwdy[le];
        tau1[2][1] = tau1[1][2];
        tauu1[0] = tau1[0][0] * u[le] + tau1[0][1] * v[le] + tau1[0][2] * w[le];
        tauu1[1] = tau1[1][0] * u[le] + tau1[1][1] * v[le] + tau1[1][2] * w[le];
        tauu1[2] = tau1[1][0] * u[le] + tau1[2][1] * v[le] + tau1[2][2] * w[le];

        tau2[0][0] = 2 * dudx[re];
        tau2[1][1] = 2 * dvdy[re];
        tau2[2][2] = 2 * dwdz[re];
        tau2[0][1] = dudy[re] + dvdx[re];
        tau2[1][0] = tau2[0][1];
        tau2[0][2] = dudz[re] + dwdx[re];
        tau2[2][0] = tau2[0][2];
        tau2[1][2] = dvdz[re] + dwdy[re];
        tau2[2][1] = tau2[1][2];
        tauu2[0] = tau2[0][0] * u[re] + tau2[0][1] * v[re] + tau2[0][2] * w[re];
        tauu2[1] = tau2[1][0] * u[re] + tau2[1][1] * v[re] + tau2[1][2] * w[re];
        tauu2[2] = tau2[1][0] * u[re] + tau2[2][1] * v[re] + tau2[2][2] * w[re];

        sf[0] = faceWeightOfLeftCell[iFace] * tauu1[0] + (1.0 - faceWeightOfLeftCell[iFace]) * tauu2[0];
        sf[1] = faceWeightOfLeftCell[iFace] * tauu1[1] + (1.0 - faceWeightOfLeftCell[iFace]) * tauu2[1];
        sf[2] = faceWeightOfLeftCell[iFace] * tauu1[2] + (1.0 - faceWeightOfLeftCell[iFace]) * tauu2[2];

        RDouble s1 = muf * (xfn[iFace] * sf[0] + yfn[iFace] * sf[1] + zfn[iFace] * sf[2]) * area[iFace];
        bCoeff[le] += s1;
    }


    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        RDouble muf = faceWeightOfLeftCell[iFace] * mu[le] + (1.0 - faceWeightOfLeftCell[iFace]) * mu[re];
        tau1[0][0] = 2 * dudx[le];
        tau1[1][1] = 2 * dvdy[le];
        tau1[2][2] = 2 * dwdz[le];
        tau1[0][1] = dudy[le] + dvdx[le];
        tau1[1][0] = tau1[0][1];
        tau1[0][2] = dudz[le] + dwdx[le];
        tau1[2][0] = tau1[0][2];
        tau1[1][2] = dvdz[le] + dwdy[le];
        tau1[2][1] = tau1[1][2];
        tauu1[0] = tau1[0][0] * u[le] + tau1[0][1] * v[le] + tau1[0][2] * w[le];
        tauu1[1] = tau1[1][0] * u[le] + tau1[1][1] * v[le] + tau1[1][2] * w[le];
        tauu1[2] = tau1[1][0] * u[le] + tau1[2][1] * v[le] + tau1[2][2] * w[le];

        tau2[0][0] = 2 * dudx[re];
        tau2[1][1] = 2 * dvdy[re];
        tau2[2][2] = 2 * dwdz[re];
        tau2[0][1] = dudy[re] + dvdx[re];
        tau2[1][0] = tau2[0][1];
        tau2[0][2] = dudz[re] + dwdx[re];
        tau2[2][0] = tau2[0][2];
        tau2[1][2] = dvdz[re] + dwdy[re];
        tau2[2][1] = tau2[1][2];
        tauu2[0] = tau2[0][0] * u[re] + tau2[0][1] * v[re] + tau2[0][2] * w[re];
        tauu2[1] = tau2[1][0] * u[re] + tau2[1][1] * v[re] + tau2[1][2] * w[re];
        tauu2[2] = tau2[1][0] * u[re] + tau2[2][1] * v[re] + tau2[2][2] * w[re];

        sf[0] = faceWeightOfLeftCell[iFace] * tauu1[0] + (1.0 - faceWeightOfLeftCell[iFace]) * tauu2[0];
        sf[1] = faceWeightOfLeftCell[iFace] * tauu1[1] + (1.0 - faceWeightOfLeftCell[iFace]) * tauu2[1];
        sf[2] = faceWeightOfLeftCell[iFace] * tauu1[2] + (1.0 - faceWeightOfLeftCell[iFace]) * tauu2[2];

        RDouble s1 = muf * (xfn[iFace] * sf[0] + yfn[iFace] * sf[1] + zfn[iFace] * sf[2]) * area[iFace];
        bCoeff[le] += s1;
        bCoeff[re] -= s1;
    }
    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        bCoeff[iCell] += (u[iCell] * dpdx[iCell] + v[iCell] * dpdy[iCell] + w[iCell] * dpdz[iCell]) * vol[iCell];
    }
    if (isUnsteady==1)
    {
        RDouble dt = GlobalDataBase::GetDoubleParaFromDB("dt");

        RDouble *pOld = reinterpret_cast<RDouble *>(grid->GetDataPtr("POld"));
        for (int iCell = 0; iCell < nTotalCell; ++iCell)
        {
            bCoeff[iCell] += ((p[iCell]-pOld[iCell])/dt) * vol[iCell];
        }
    }
    DelPointer2(tau1);
    DelPointer2(tau2);
    delete [] tauu1;
    delete [] tauu2;
    delete [] sf;
}

void IncomEnergyEqCalculator::GetResidual(Grid *gridIn, vector<RDouble>& res)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble resNow = 0.0;

    grid->GetData("EnthalpyResNow", &resNow, PHDOUBLE, 1);

    res.push_back(resNow);
}

void IncomEnergyEqCalculator::InitialUnsteadyVar(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    PHString1D phiNameList;
    phiNameList.push_back("Enthalpy");

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

void IncomEnergyEqCalculator::UpdateUnsteadyVariable(Grid *grid)
{
    UnstructGrid *gridIn = UnstructGridCast(grid); 
    std::vector<std::string> phiNameList;
    phiNameList.push_back("Enthalpy");

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
    else
    {
        cout << "we not support Tran Method : " << TranCalcMethod  << endl;
    }

}

void IncomEnergyEqCalculator::UpdateBCValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    string *varNameIncom = reinterpret_cast<string *>(GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    string varName = varNameIncom[solverIndex];

    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr("T"));
    RDouble *Enthalpy = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
    RDouble *cp = reinterpret_cast<RDouble *>(grid->GetDataPtr("cp"));
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble * faceCenterX = grid->GetFaceCenterX();
    RDouble * faceCenterY = grid->GetFaceCenterY();
    RDouble * faceCenterZ = grid->GetFaceCenterZ();
    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
    
    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            RDouble spBound;

            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                spBound = phi[le];
                if (bcData->IsExist("initT", PHDOUBLE, 1))
                {
                    bcData->GetData("initT", &spBound, PHDOUBLE, 1);
                }
                else if (bcData->IsExist("TPolynomial", PHDOUBLE, 4))
                {
                    RDouble *TPolynomial = new RDouble[4];
                    bcData->GetData("TPolynomial", TPolynomial, PHDOUBLE, 4);
                    spBound = TPolynomial[0] + TPolynomial[1] * faceCenterX[iFacelocal] + TPolynomial[2] * faceCenterY[iFacelocal]
                            + TPolynomial[3] * faceCenterZ[iFacelocal];
                }

                phi[re] = spBound;
                Enthalpy[re] = cp[re] * spBound;
            }

        }
        else if (bcType == PHENGLEI::SYMMETRY)
        {
            //! not this boundary!
        }
        else if (bcType == PHENGLEI::FARFIELD)
        {
            RDouble spBound;
            if (bcData->IsExist("initT", PHDOUBLE, 1))
            {
                bcData->GetData("initT", &spBound, PHDOUBLE, 1);
            }
            else
            {
                cout << " No initial value has assigned to the boundary " << endl;
            }

            for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                if (FaceFlux[iFacelocal] < 0)
                {
                    phi[re] = spBound;
                    Enthalpy[re] = cp[re] * spBound;
                }
                else
                {
                    phi[re] = phi[le];
                    Enthalpy[re] =  cp[iFacelocal] * phi[le];
                }
            }
        }
        else if (bcType == PHENGLEI::INFLOW)
        {
            RDouble spBound;

            if (bcData->IsExist("initT", PHDOUBLE, 1))
            {
                bcData->GetData("initT", &spBound, PHDOUBLE, 1);

                for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int re = rightCellOfFace[iFacelocal];
                    phi[re] = spBound;
                    Enthalpy[re] = cp[re] * spBound;
                }
            }
            else if (bcData->IsExist("TPolynomial", PHDOUBLE, 3))
            {
                RDouble *TPolynomial = new RDouble[3];
                bcData->GetData("TPolynomial", TPolynomial, PHDOUBLE, 3);
                for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int re = rightCellOfFace[iFacelocal];
                    spBound = TPolynomial[0] * faceCenterX[iFacelocal] + TPolynomial[1] * faceCenterY[iFacelocal] 
                        + TPolynomial[2] * faceCenterZ[iFacelocal];

                    phi[re] = spBound;
                    Enthalpy[re] = cp[re] * spBound;
                }
            }
            else
            {
                cout << " No initial value has assigned to the boundary " << endl;
            }
        }
        else if (bcType == PHENGLEI::OUTFLOW)
        {
            RDouble spBound;

            if (bcData->IsExist("initT", PHDOUBLE, 1))
            {
                bcData->GetData("initT", &spBound, PHDOUBLE, 1);
            }

            for (int iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                phi[re] = phi[le];
                Enthalpy[re] = cp[re] * phi[re];
                if (FaceFlux[iFacelocal] < 0.0)
                {            
                    phi[re] = spBound;
                    Enthalpy[re] = cp[re] * spBound;
                }
            }
        }
        else if (bcType == PHENGLEI::PRESSURE_INLET)
        {
            RDouble spBound;
            if (bcData->IsExist("initT", PHDOUBLE, 1))
            {
                bcData->GetData("initT", &spBound, PHDOUBLE, 1);

                for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int re = rightCellOfFace[iFacelocal];
                    phi[re] = spBound;
                    Enthalpy[re] = cp[re] * spBound;
                }
            }
            else if (bcData->IsExist("TPolynomial", PHDOUBLE, 3))
            {
                RDouble *TPolynomial = new RDouble[3];
                bcData->GetData("TPolynomial", TPolynomial, PHDOUBLE, 3);
                for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int re = rightCellOfFace[iFacelocal];
                    spBound = TPolynomial[0] * faceCenterX[iFacelocal] + TPolynomial[1] * faceCenterY[iFacelocal] 
                        + TPolynomial[2] * faceCenterZ[iFacelocal];

                    phi[re] = spBound;
                    Enthalpy[re] = cp[re] * spBound;
                }
            }
            else
            {
                cout << " No initial value has assigned to the boundary " << endl;
            }
        }
        else if (bcType == PHENGLEI::PRESSURE_OUTLET)
        {
            RDouble spBound;

            if (bcData->IsExist("initT", PHDOUBLE, 1))
            {
                bcData->GetData("initT", &spBound, PHDOUBLE, 1);
            }

            for (int iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                phi[re] = phi[le];
                Enthalpy[re] = cp[re] * phi[re];
                if (FaceFlux[iFacelocal] < 0.0)
                {            
                    phi[re] = spBound;
                    Enthalpy[re] = cp[re] * spBound;
                }
            }
        }
    }

    CommunicateAnInterfaceVar(phi);
    CommunicateAnInterfaceVar(Enthalpy);
}

void IncomEnergyEqCalculator::InitFlowAsRestart(Grid *gridIn)
{
    IncomScalarEqCalculator::InitFlowAsRestart(gridIn);
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();

    RDouble initCPg = GlobalDataBase::GetDoubleParaFromDB("initCPg");
    RDouble initT = GlobalDataBase::GetDoubleParaFromDB("initT");
    RDouble initPhi = initCPg * initT;
    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(varNameIncom[solverIndex]));
    PHSPACE::SetField(phi, initPhi, nTotal);

    GAS_SPACE::gas->SetGasCp("gasCp");
    GAS_SPACE::gas->SetGasK("gasK");
    GAS_SPACE::gas->InitParameterForK();

    RDouble *t = reinterpret_cast<RDouble *>(grid->GetDataPtr("T"));
    PHSPACE::SetField(t, initT, nTotal);

    RDouble *cpg = reinterpret_cast<RDouble *>(grid->GetDataPtr("cp"));
    PHSPACE::SetField(cpg, initCPg, nTotal);

    RDouble initk = GlobalDataBase::GetDoubleParaFromDB("initK");
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("k"));
    PHSPACE::SetField(k, initk, nTotal);
}

void IncomEnergyEqCalculator::AllocateGlobalVar(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    IncomScalarEqCalculator::AllocateGlobalVar(grid);

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble *t = NewPointer<RDouble>(nTotal);
    grid->UpdateDataPtr("T", t);

    RDouble *cpg = NewPointer<RDouble>(nTotal);
    grid->UpdateDataPtr("cp", cpg);

    RDouble *k = NewPointer<RDouble>(nTotal);
    grid->UpdateDataPtr("k", k);

    RDouble *wallq = NewPointer<RDouble>(nBoundFace);
    grid->UpdateDataPtr("wallq", wallq);

    RDouble *heatFlux = NewPointer<RDouble>(nTotal);
    grid->UpdateDataPtr("heatFlux", heatFlux);
}

void IncomEnergyEqCalculator::IncompressibleInitial(Grid * gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    RDouble usq;
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int nTotalInnerFace = nTotalFace - nBoundFace;
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int nTotal = nTotalCell + nBoundFace;
    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *t = reinterpret_cast<RDouble *>(grid->GetDataPtr("T"));
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("k"));
    RDouble *h = reinterpret_cast<RDouble *>(grid->GetDataPtr("Enthalpy"));
    RDouble *cpg = reinterpret_cast<RDouble *>(grid->GetDataPtr("cp"));

    GAS_SPACE::gas->UpdateK(grid, k);
    GAS_SPACE::gas->UpdateCp(grid, cpg);
    RDouble **DiffusionCoeff = reinterpret_cast <RDouble **> (GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    int solverIndex = GetSolverIndex();
    
    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        DiffusionCoeff[solverIndex][iCell] = k[iCell];

        usq = (u[iCell] * u[iCell] + v[iCell] * v[iCell] + w[iCell] * w[iCell]) / 2;
        h[iCell] = cpg[iCell] * t[iCell] + usq;
    }

    UpdateBCValue(grid);

    InitialUnsteadyVar(grid);
}

void IncomEnergyEqCalculator::UpdateProperties(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *t = reinterpret_cast<RDouble *>(grid->GetDataPtr("T"));
    RDouble *h = reinterpret_cast<RDouble *>(grid->GetDataPtr("Enthalpy"));
    RDouble *cp = reinterpret_cast<RDouble *>(grid->GetDataPtr("cp"));
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("k"));

    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nIFace = grid->GetNIFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble tMin = 1.0;
    RDouble tMax = 5000.0;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        t[iCell] = h[iCell] / cp[iCell];

        if (t[iCell] < tMin)
        {
            t[iCell] = tMin;
        }
        else if (t[iCell] > tMax)
        {
            t[iCell] = tMax;
        }
    }

    UpdateBCValue(grid);

    GAS_SPACE::gas->UpdateK(grid, k);
    GAS_SPACE::gas->UpdateCp(grid, cp);

    CommunicateAnInterfaceVar(h);
    CommunicateAnInterfaceVar(t);
    CommunicateAnInterfaceVar(k);
    CommunicateAnInterfaceVar(cp);

}

}