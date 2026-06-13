#include <iostream>
#include <memory>
#include <cmath>
#include "IncomSATurbKEqCalculator.h"



using namespace std;

namespace PHSPACE
{
IncomSATurbKEqCalculator::IncomSATurbKEqCalculator():IncomScalarEqCalculator()
{
    
}

IncomSATurbKEqCalculator::~IncomSATurbKEqCalculator()
{

}

void IncomSATurbKEqCalculator::InitFlowAsRestart(Grid *gridIn)
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

    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    PHSPACE::SetField(visl, 0.0, nTotal);

    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    PHSPACE::SetField(vist, 0.0, nTotal);

    RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
    PHSPACE::SetField(gen, 0.0, nTotal);

    RDouble *des = reinterpret_cast<RDouble *>(grid->GetDataPtr("des"));
    PHSPACE::SetField(des, 0.0, nTotal);

    RDouble *gra = reinterpret_cast<RDouble *>(grid->GetDataPtr("gra"));
    PHSPACE::SetField(gra, 0.0, nTotal);

    RDouble *fv1 = reinterpret_cast<RDouble *>(grid->GetDataPtr("fv1"));
    PHSPACE::SetField(fv1, 0.0, nTotal);

    RDouble *dKineticdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dKineticdx"));
    RDouble *dKineticdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dKineticdy"));
    RDouble *dKineticdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dKineticdz"));
    PHSPACE::SetField(dKineticdx, 0.0, nTotal);
    PHSPACE::SetField(dKineticdy, 0.0, nTotal);
    PHSPACE::SetField(dKineticdz, 0.0, nTotal);

    RDouble *WD = reinterpret_cast<RDouble *>(grid->GetDataPtr("WD"));
    PHSPACE::SetField(WD, 0.0, nTotal);
}

void IncomSATurbKEqCalculator::AllocateGlobalVar(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    IncomScalarEqCalculator::AllocateGlobalVar(grid);

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble *visl = NewPointer<RDouble>(nTotal);
    grid->UpdateDataPtr("visl", visl);

    RDouble *vist = NewPointer<RDouble>(nTotal);
    grid->UpdateDataPtr("vist", vist);

    RDouble *gen = NewPointer<RDouble>(nTotal);
    grid->UpdateDataPtr("gen", gen);

    RDouble *des = NewPointer<RDouble>(nTotal);
    grid->UpdateDataPtr("des", des);

    RDouble *gra = NewPointer<RDouble>(nTotal);
    grid->UpdateDataPtr("gra", gra);

    RDouble *fv1 = NewPointer<RDouble>(nTotal);
    grid->UpdateDataPtr("fv1", fv1);

    RDouble *WD = NewPointer<RDouble>(nTotal);
    grid->UpdateDataPtr("WD", WD);
}

//! The initialization of each solver depends on the associated variables of other solvers.
void IncomSATurbKEqCalculator::IncompressibleInitial(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    RDouble *fv1  = reinterpret_cast<RDouble *>(grid->GetDataPtr("fv1"));
    RDouble *rho  = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));

    RDouble delatv = 0.6666;
    RDouble Cb2    = 0.622;
    RDouble Cv1    = 7.1;
    RDouble Cb1    = 0.1355;
    RDouble kappa  = 0.4187;
    RDouble Cw2    = 0.3;
    RDouble Cw3    = 2;
    RDouble *mu    = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble chi = rho[iCell] * k[iCell] / mu[iCell];
        fv1[iCell] = (pow(chi, 3)) / ((pow(chi, 3)) + (pow(Cv1, 3)));
        vist[iCell] = rho[iCell] * k[iCell] * fv1[iCell];
        visl[iCell] = mu[iCell];
        mu[iCell] = visl[iCell] + vist[iCell];
        
        if (vist[iCell] > 1e5 * visl[iCell])
        {
            vist[iCell] = visl[iCell] * 1e5;
        }
    }
    UpdateBCValue(grid);
}

void IncomSATurbKEqCalculator::CalcOtherMatrixACoeff(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
    RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));
    RDouble *gra = reinterpret_cast<RDouble *>(grid->GetDataPtr("gra"));
    RDouble *des = reinterpret_cast<RDouble *>(grid->GetDataPtr("des"));
    RDouble *vol = reinterpret_cast<RDouble *>(grid->GetCellVolume());
    RDouble *dudx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdx"));
    RDouble *dudy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdy"));
    RDouble *dudz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdz"));
    RDouble *dvdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdx"));
    RDouble *dvdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdy"));
    RDouble *dvdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdz"));
    RDouble *dwdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdx"));
    RDouble *dwdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdy"));
    RDouble *dwdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdz"));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *fv1 = reinterpret_cast<RDouble *>(grid->GetDataPtr("fv1"));
    RDouble *wDist = reinterpret_cast<RDouble *>(grid->GetWallDist());
    RDouble *WD = reinterpret_cast<RDouble *>(grid->GetDataPtr("WD"));

    RDouble delatv = 0.6666;
    RDouble Cb2    = 0.622;
    RDouble Cv1    = 7.1;
    RDouble Cb1    = 0.1355;
    RDouble kappa  = 0.4187;
    RDouble Cw2    = 0.3;
    RDouble Cw3    = 2;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {

        RDouble wd = wDist[iCell];
        WD[iCell] = wd;

        RDouble _dudy = dudy[iCell];
        RDouble _dudz = dudz[iCell];
        RDouble _dvdx = dvdx[iCell];
        RDouble _dvdz = dvdz[iCell];
        RDouble _dwdx = dwdx[iCell];
        RDouble _dwdy = dwdy[iCell];

        RDouble chi = rho[iCell] * k[iCell] / visl[iCell];
        fv1[iCell] = (pow(chi, 3)) / ((pow(chi, 3)) + (pow(Cv1, 3)));

        RDouble fv2 = 1 - chi / (1 + chi * fv1[iCell]);

        RDouble S = sqrt(pow((_dwdy - _dvdz), 2) + pow((_dudz - _dwdx), 2) + pow((_dvdx - _dudy), 2));
        RDouble S_0 = S + k[iCell] * fv2 / (kappa * kappa * wd * wd);
        RDouble S_ = MAX(S_0, 0.3 * S);

        RDouble Cw1 = Cb1 / (kappa * kappa) + (1 + Cb2) / delatv;
        RDouble r = MIN((k[iCell]) / (S_ * kappa * kappa * wd * wd), 10.0);

        RDouble g = r + Cw2 * (pow(r, 6) - r);
        RDouble fw = g * pow(((1+Cw3)/(pow(g,6)+pow(Cw3,6))),1.0/6.0);

        des[iCell] = Cw1 * rho[iCell] * fw * (k[iCell] / wd) * (k[iCell] / wd);

        diagMatrixCoeff[iCell] +=  des[iCell] / (k[iCell] + 1e-20) * vol[iCell];
    }
}

void IncomSATurbKEqCalculator::CalcOtherbCoeff(Grid *gridIn, int iEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
    RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));
    RDouble *gra = reinterpret_cast<RDouble *>(grid->GetDataPtr("gra"));
    RDouble *des = reinterpret_cast<RDouble *>(grid->GetDataPtr("des"));
    RDouble *vol = reinterpret_cast<RDouble *>(grid->GetCellVolume());

    RDouble *dKineticdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dKineticdx"));
    RDouble *dKineticdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dKineticdy"));
    RDouble *dKineticdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dKineticdz"));


    RDouble *WD = reinterpret_cast<RDouble *>(grid->GetDataPtr("WD"));
    RDouble *dudx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdx"));
    RDouble *dudy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdy"));
    RDouble *dudz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdz"));
    RDouble *dvdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdx"));
    RDouble *dvdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdy"));
    RDouble *dvdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdz"));
    RDouble *dwdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdx"));
    RDouble *dwdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdy"));
    RDouble *dwdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdz"));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *fv1 = reinterpret_cast<RDouble *>(grid->GetDataPtr("fv1"));
    RDouble *wDist = reinterpret_cast<RDouble *>(grid->GetWallDist());

    RDouble delatv = 0.6666;
    RDouble Cb2    = 0.622;
    RDouble Cv1    = 7.1;
    RDouble Cb1    = 0.1355;
    RDouble kappa  = 0.4187;
    RDouble Cw2    = 0.3;
    RDouble Cw3    = 2;

    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {

        RDouble wd = wDist[iCell];
        WD[iCell] = wd;

        RDouble _dudy = dudy[iCell];
        RDouble _dudz = dudz[iCell];
        RDouble _dvdx = dvdx[iCell];
        RDouble _dvdz = dvdz[iCell];
        RDouble _dwdx = dwdx[iCell];
        RDouble _dwdy = dwdy[iCell];

        RDouble chi = rho[iCell] * k[iCell] / visl[iCell];
        fv1[iCell] = (pow(chi, 3)) / ((pow(chi, 3)) + (pow(Cv1, 3)));

        RDouble fv2 = 1 - chi / (1 + chi * fv1[iCell]);
        RDouble S = sqrt(pow((_dwdy - _dvdz), 2) + pow((_dudz - _dwdx), 2) + pow((_dvdx - _dudy), 2));
        RDouble S_0 = S + k[iCell] * fv2 / (kappa * kappa * wd * wd);
        RDouble S_ = MAX(S_0, 0.3 * S);

        gen[iCell] = Cb1 * rho[iCell] * S_ * k[iCell];
        gen[iCell] = MAX(gen[iCell], 0.0);

        RDouble _dkdx = dKineticdx[iCell];
        RDouble _dkdy = dKineticdy[iCell];
        RDouble _dkdz = dKineticdz[iCell];
        gra[iCell] = Cb2 * rho[iCell] * (_dkdx * _dkdx + _dkdy * _dkdy + _dkdz * _dkdz) / delatv;

    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        bCoeff[iCell] += (gen[iCell] + gra[iCell] ) * vol[iCell];
    }
}

void IncomSATurbKEqCalculator::calVelocityInletBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();

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
    int *rightCellOfFace = grid->GetRightCellOfFace();
    string *varNameIncom = reinterpret_cast<string *>(GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    string varName = varNameIncom[solverIndex];

    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
    
    RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dx"));
    RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dy"));
    RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dz"));
    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
    RDouble **DiffusionCoeff = reinterpret_cast<RDouble **>(GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    
    for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        if (bcData->IsExist("init" + varName, PHDOUBLE, 1))
        {
            int iFacelocal = (*faceIndex)[iFace];
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

            diagMatrixCoeff[le] += area2 * matrixCoeff;
            bCoeff[le] += area2 * phi[re] + dphidx[le] * ax + dphidy[le] * ay + dphidz[le] * az;

            if (FaceFlux[iFacelocal] < 0)
            {
                diagMatrixCoeff[le] -= FaceFlux[iFacelocal] * matrixCoeff;
                bCoeff[le] -= FaceFlux[iFacelocal] * phi[re];
            }
        }
        else
        {
            int iFacelocal = (*faceIndex)[iFace];
            int le = leftCellOfFace[iFacelocal];
            int re = rightCellOfFace[iFacelocal];
            if (FaceFlux[iFacelocal] < 0)
            {
                diagMatrixCoeff[le] -= FaceFlux[iFacelocal] * matrixCoeff;
                bCoeff[le] -= FaceFlux[iFacelocal] * phi[re];
            }
        }
    }
}

void IncomSATurbKEqCalculator::calPressureOutletBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();

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
    int *rightCellOfFace = grid->GetRightCellOfFace();
    string *varNameIncom = reinterpret_cast<string *>(GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    string varName = varNameIncom[solverIndex];

    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
    
    RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dx"));
    RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dy"));
    RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dz"));
    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
    RDouble **DiffusionCoeff = reinterpret_cast<RDouble **>(GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    
    for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        if (bcData->IsExist("init" + varName, PHDOUBLE, 1))
        {
            int iFacelocal = (*faceIndex)[iFace];  
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

            diagMatrixCoeff[le] += area2 * matrixCoeff;
            bCoeff[le] += area2 * phi[re] + dphidx[le] * ax + dphidy[le] * ay + dphidz[le] * az;

            if (FaceFlux[iFacelocal] < 0)
            {
                diagMatrixCoeff[le] -= FaceFlux[iFacelocal] * matrixCoeff;
                bCoeff[le] -= FaceFlux[iFacelocal] * phi[re];
            }
        }
        else
        {
            int iFacelocal = (*faceIndex)[iFace];
            int le = grid->GetLeftCellOfFace()[iFacelocal];
            int re = rightCellOfFace[iFacelocal];
            if (FaceFlux[iFacelocal] < 0)
            {
                diagMatrixCoeff[le] -= FaceFlux[iFacelocal] * matrixCoeff;
                bCoeff[le] -= FaceFlux[iFacelocal] * phi[re];
            }
        }
    }
}


void IncomSATurbKEqCalculator::calWallBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();

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
    string *varNameIncom = reinterpret_cast<string *>(GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    string varName = varNameIncom[solverIndex];

    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dx"));
    RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dy"));
    RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dz"));
    RDouble **DiffusionCoeff = reinterpret_cast<RDouble **>(GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    
    for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        int iFacelocal = (*faceIndex)[iFace]; 
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

        diagMatrixCoeff[le] += area2  * matrixCoeff;
        bCoeff[le] += dphidx[le] * ax + dphidy[le] * ay + dphidz[le] * az;
    }
}

void IncomSATurbKEqCalculator::GetResidual(Grid *gridIn, vector<RDouble> &res)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    RDouble resNow = 0.0;

    grid->GetData("KineticResNow", &resNow, PHDOUBLE, 1);

    res.push_back(resNow);
}

void IncomSATurbKEqCalculator::InitialUnsteadyVar(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    PHString1D phiNameList;
    phiNameList.push_back("Kinetic");

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

void IncomSATurbKEqCalculator::UpdateUnsteadyVariable(Grid *grid)
{
    UnstructGrid *gridIn = UnstructGridCast(grid); 
    std::vector<std::string> phiNameList;
    phiNameList.push_back("Kinetic");

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

void IncomSATurbKEqCalculator::UpdateBCValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    string varName = varNameIncom[solverIndex];
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *area = grid->GetFaceArea();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();
    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
    RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dx"));
    RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dy"));
    RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr("d" + varName + "dz"));

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        if (bcType == PHENGLEI::INFLOW)
        {
            RDouble kb = 0.0;
            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                kb = phi[le];
                if (bcData->IsExist("init" + varName, PHDOUBLE, 1))
                {
                    bcData->GetData("init" + varName, &kb, PHDOUBLE, 1);
                }
                phi[re] = kb;
            }
        }
        else if (bcType == PHENGLEI::OUTFLOW)
        {
            RDouble kb = 0.0;
            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                RDouble dx = xfc[iFacelocal] - xcc[le];
                RDouble dy = yfc[iFacelocal] - ycc[le];
                RDouble dz = zfc[iFacelocal] - zcc[le];
                phi[re] = phi[le] + dphidx[le] * dx + dphidy[le] * dy + dphidz[le] * dz;
                if (FaceFlux[iFacelocal] < 0.0)
                {
                    if (bcData->IsExist("init" + varName, PHDOUBLE, 1))
                    {
                        bcData->GetData("init" + varName, &kb, PHDOUBLE, 1);
                        phi[re] = RDouble(kb);
                    }
                    else
                    {
                        std::cout << " Please assign a value " << std::endl;
                    }
                }
            }
        }
        else if (bcType == PHENGLEI::PRESSURE_INLET)
        {
            RDouble kb = 0.0;
            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                kb = phi[le];
                if (bcData->IsExist("init" + varName, PHDOUBLE, 1))
                {
                    bcData->GetData("init" + varName, &kb, PHDOUBLE, 1);
                }
                phi[re] = kb;
            }
        }
        else if (bcType == PHENGLEI::PRESSURE_OUTLET)
        {
            RDouble kb = 0.0;
            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];

                RDouble dx = xfc[iFacelocal] - xcc[le];
                RDouble dy = yfc[iFacelocal] - ycc[le];
                RDouble dz = zfc[iFacelocal] - zcc[le];
                phi[re] = phi[le] + dphidx[le] * dx + dphidy[le] * dy + dphidz[le] * dz;
                if (FaceFlux[iFacelocal] < 0.0)
                {
                    if (bcData->IsExist("init" + varName, PHDOUBLE, 1))
                    {
                        bcData->GetData("init" + varName, &kb, PHDOUBLE, 1);
                        phi[re] = RDouble(kb);
                    }
                    else
                    {
                        std::cout << " Please assign a value " << std::endl;
                    }
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

void IncomSATurbKEqCalculator::SetDiffusionCoeff(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nIFace = grid->GetNIFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *rhoFace = reinterpret_cast<RDouble *>(grid->GetDataPtr("rhoFace"));

    RDouble delatv = 0.6666;
    RDouble Cb2    = 0.622;
    RDouble Cv1    = 7.1;
    RDouble Cb1    = 0.1355;
    RDouble kappa  = 0.4187;
    RDouble Cw2    = 0.3;
    RDouble Cw3    = 2;
    RDouble **DiffusionCoeff = reinterpret_cast <RDouble **> (GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    int solverIndex = GetSolverIndex();
    
    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        DiffusionCoeff[solverIndex][iCell] = (visl[iCell] + rho[iCell] * k[iCell]) / delatv;
    }

    for (int iFace = 0; iFace < nBoundFace; ++iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        DiffusionCoeff[solverIndex][re] = (visl[le] + rhoFace[iFace] * k[le]) / delatv;
    }
    CommunicateAnInterfaceVar(DiffusionCoeff[solverIndex]);
}

void IncomSATurbKEqCalculator::CalcGrad(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    GradientCalculation(grid, "U", "dUdx", "dUdy", "dUdz");
    GradientCalculation(grid, "V", "dVdx", "dVdy", "dVdz");
    GradientCalculation(grid, "W", "dWdx", "dWdy", "dWdz");

    RDouble *dudx  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdx"));
    RDouble *dudy  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdy"));
    RDouble *dudz  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdz"));
    RDouble *dvdx  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdx"));
    RDouble *dvdy  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdy"));
    RDouble *dvdz  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdz"));
    RDouble *dwdx  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdx"));
    RDouble *dwdy  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdy"));
    RDouble *dwdz  = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdz"));
    CommunicateAnInterfaceVar(dudx);
    CommunicateAnInterfaceVar(dvdx);
    CommunicateAnInterfaceVar(dwdx);
    CommunicateAnInterfaceVar(dudy);
    CommunicateAnInterfaceVar(dvdy);
    CommunicateAnInterfaceVar(dwdy);
    CommunicateAnInterfaceVar(dudz);
    CommunicateAnInterfaceVar(dvdz);
    CommunicateAnInterfaceVar(dwdz);

    IncomScalarEqCalculator::CalcGrad(grid);
}

void IncomSATurbKEqCalculator::UpdateAfterIterloop(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nINFace = grid->GetNIFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
    RDouble *fv1 = reinterpret_cast<RDouble *>(grid->GetDataPtr("fv1"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *rhoFace = reinterpret_cast<RDouble *>(grid->GetDataPtr("rhoFace"));
    RDouble *mu   = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
    
    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        RDouble vistold = vist[iCell];
        vist[iCell] = vistold + (rho[iCell] * k[iCell] * fv1[iCell] - vistold);
        mu[iCell] = visl[iCell] + vist[iCell];
        
        if (vist[iCell] > 1e5 * visl[iCell])
        {
            mu[iCell] = visl[iCell] * 1e5;
            vist[iCell] = visl[iCell] * 1e5;
        }
    }

    for (int iFace = 0; iFace < nBoundFace - nINFace; ++iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        vist[re] = vist[le];
        mu[re] = mu[le];
    }

    UnstructBCSet* unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        if (bcType == 2)
        {
            vector<int>* faceIndex = bcRegion->GetFaceIndex();
            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                mu[re] = visl[le];
            }
        }
    }

    CommunicateAnInterfaceVar(mu);
    CommunicateAnInterfaceVar(visl);
    CommunicateAnInterfaceVar(vist);

}

}