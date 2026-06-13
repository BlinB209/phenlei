#include "MomEqCalculator.h"
#include "TK_Log.h"
#include "TK_Time.h"
#include "Post_Probes.h"
#include "Post_Visual.h"
#include "Geo_UnstructBC.h"
#include "Data_Param.h"

#include "Post_WriteTecplot.h"
#include "PHMpi.h"


namespace PHSPACE
{
    
    MomEqCalculator::MomEqCalculator():IncomScalarEqCalculator()
    {
    }

    void MomEqCalculator::solveMomentumEquations(Grid *gridIn)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);
        CalcGrad(grid);
        
        constructMatrixACoeff(grid, IDX::S_IU);
        CoefPrepForPres(grid);

        solveScalarEquation(grid, IDX::S_IU);
        solveScalarEquation(grid, IDX::S_IV);
        solveScalarEquation(grid, IDX::S_IW);

        updateBCValue(grid);
    }

    void MomEqCalculator::CalcGrad(Grid *gridIn)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);
        GradientCalculation(grid, "P","dpdx", "dpdy", "dpdz");
        GradientCalculation(grid, "U", "dUdx", "dUdy", "dUdz");
        GradientCalculation(grid, "V", "dVdx", "dVdy", "dVdz");
        GradientCalculation(grid, "W", "dWdx", "dWdy", "dWdz");
    }

    void MomEqCalculator::solveScalarEquation(Grid* gridIn, int iEquation)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);
        int solverIndex = GetSolverIndex(iEquation);
        RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("diagMatrixCoeff"));
        RDouble *bCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("bCoeff"));
        RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("upperMatrixCoeff"));
        RDouble *lowerMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("lowerMatrixCoeff"));
        string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));

        constructBCoeff(grid, solverIndex);
        RDouble *q = reinterpret_cast<RDouble *> (grid->GetDataPtr(varNameIncom[solverIndex]));
        calculateLinearEquation(grid, solverIndex, q, diagMatrixCoeff, bCoeff, upperMatrixCoeff , lowerMatrixCoeff);
    }

    void MomEqCalculator::calcFaceFlux(Grid*  gridIn)
    {
        UnstructGrid* grid = UnstructGridCast(gridIn);
        int nTotalCell = grid->GetNTotalCell();
        int nIFace = grid->GetNIFace();
        int nTotalFace = grid->GetNTotalFace();
        int nBoundFace = grid->GetNBoundFace();
        int nTotalInFace = nTotalFace - nBoundFace;
        RDouble *pc = reinterpret_cast<RDouble *> (grid->GetDataPtr("pc"));
        RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("diagMatrixCoeff"));
        RDouble *bCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("bCoeff"));

        int *leftCellOfFace = grid->GetLeftCellOfFace();
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
        RDouble FaceFluxRelaxCoeff = GlobalDataBase::GetDoubleParaFromDB("FaceFluxRelaxCoeff");
        RDouble *u = reinterpret_cast<RDouble *> (grid->GetDataPtr("U"));
        RDouble *v = reinterpret_cast<RDouble *> (grid->GetDataPtr("V"));
        RDouble *w = reinterpret_cast<RDouble *> (grid->GetDataPtr("W"));
        RDouble *p = reinterpret_cast<RDouble *> (grid->GetDataPtr("P"));
        RDouble *InverseAp = reinterpret_cast<RDouble *> (grid->GetDataPtr("InverseAp"));
        RDouble *bfU = reinterpret_cast<RDouble *> (grid->GetDataPtr("bfU"));
        RDouble *bfV = reinterpret_cast<RDouble *> (grid->GetDataPtr("bfV"));
        RDouble *bfW = reinterpret_cast<RDouble *> (grid->GetDataPtr("bfW"));
        RDouble *dpdx = reinterpret_cast<RDouble *> (grid->GetDataPtr("dpdx"));
        RDouble *dpdy = reinterpret_cast<RDouble *> (grid->GetDataPtr("dpdy"));
        RDouble *dpdz = reinterpret_cast<RDouble *> (grid->GetDataPtr("dpdz"));
        RDouble* rho  = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
        RDouble* mu   = reinterpret_cast<RDouble*>(grid->GetDataPtr("mu"));
        RDouble* FaceFlux = reinterpret_cast<RDouble*>(grid->GetDataPtr("FaceFlux"));
        RDouble* dun  = reinterpret_cast<RDouble*>(grid->GetDataPtr("dun"));
        RDouble *rhoFace = reinterpret_cast<RDouble*>(grid->GetDataPtr("rhoFace"));
        RDouble *faceWeightOfLeftCell = reinterpret_cast<RDouble *>(grid->GetDataPtr("faceWeightOfLeftCell"));
        RDouble* cRho = reinterpret_cast<RDouble*>(grid->GetDataPtr("cRho"));
        RDouble* ptotal = reinterpret_cast<RDouble*>(grid->GetDataPtr("Ptotal"));

        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

        for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
        {
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];

            RDouble ax = area[iFace] * xfn[iFace];
            RDouble ay = area[iFace] * yfn[iFace];
            RDouble az = area[iFace] * zfn[iFace];

            RDouble dx = xcc[re] - xcc[le];
            RDouble dy = ycc[re] - ycc[le];
            RDouble dz = zcc[re] - zcc[le];

            RDouble ad = ax * dx + ay * dy + az * dz;
            RDouble a2 = ax * ax + ay * ay + az * az;

            RDouble ubar_a = (faceWeightOfLeftCell[iFace] * u[le] + (1.0 - faceWeightOfLeftCell[iFace]) * u[re]) * ax 
                + (faceWeightOfLeftCell[iFace] * v[le] + (1.0 - faceWeightOfLeftCell[iFace]) * v[re]) * ay 
                + (faceWeightOfLeftCell[iFace] * w[le] + (1.0 - faceWeightOfLeftCell[iFace]) * w[re]) * az;

            RDouble pDiff;
            RDouble pdf;

            int bodyForceFlag = GlobalDataBase::GetIntParaFromDB("bodyForceFlag");
            if (bodyForceFlag != 1)
            {
                pDiff = (InverseAp[le] * faceWeightOfLeftCell[iFace] + InverseAp[re] * (1.0 - faceWeightOfLeftCell[iFace])) *( p[re] - p[le]) ;
                pdf = ((InverseAp[le] * faceWeightOfLeftCell[iFace] * dpdx[le] + InverseAp[re] * (1.0 - faceWeightOfLeftCell[iFace]) * dpdx[re]) * dx +
                    (InverseAp[le] * faceWeightOfLeftCell[iFace] * dpdy[le] + InverseAp[re] * (1.0 - faceWeightOfLeftCell[iFace]) * dpdy[re]) * dy+
                    (InverseAp[le] * faceWeightOfLeftCell[iFace] * dpdz[le] + InverseAp[re] * (1.0 - faceWeightOfLeftCell[iFace]) * dpdz[re]) * dz);
            }
            else
            {
                RDouble dxfL = xfc[iFace] - xcc[le];
                RDouble dyfL = yfc[iFace] - ycc[le];
                RDouble dzfL = zfc[iFace] - zcc[le];
                RDouble dxfR = xfc[iFace] - xcc[re];
                RDouble dyfR = yfc[iFace] - ycc[re];
                RDouble dzfR = zfc[iFace] - zcc[re];

                pDiff = (InverseAp[le] * faceWeightOfLeftCell[iFace] + InverseAp[re] * (1.0 - faceWeightOfLeftCell[iFace])) * (p[re] - p[le] +
                    (bfU[re] * dxfR + bfV[re] * dyfR + bfW[re] * dzfR)- (bfU[le] * dxfL + bfV[le] * dyfL + bfW[le] * dzfL));
                pdf = ((InverseAp[le] * faceWeightOfLeftCell[iFace] * (dpdx[le] - bfU[le]) + InverseAp[re] * (1.0 - faceWeightOfLeftCell[iFace]) * (dpdx[re] - bfU[re])) * dx +
                    (InverseAp[le] * faceWeightOfLeftCell[iFace] * (dpdy[le] - bfV[le]) + InverseAp[re] * (1.0 - faceWeightOfLeftCell[iFace]) * (dpdy[re] - bfV[re])) * dy +
                       (InverseAp[le] * faceWeightOfLeftCell[iFace] * (dpdz[le] - bfW[le]) + InverseAp[re] * (1.0 - faceWeightOfLeftCell[iFace]) * (dpdz[re] - bfW[re])) * dz);
            };

            pDiff *= a2/ad;
            pdf*= a2 / ad;

            FaceFlux[iFace] = rhoFace[iFace] * (ubar_a - pDiff + pdf + FaceFluxRelaxCoeff * dun[iFace]);
        }

        for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
            int bcType = bcRegion->GetBCType();
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            Data_Param* bcData = bcRegion->GetBCParamDataBase();

            if (bcType == PHENGLEI::SOLID_SURFACE)
            {
            }
            else if (bcType == PHENGLEI::SYMMETRY)
            {
            }
            else if (bcType == PHENGLEI::FARFIELD)
            {
                RDouble refflowU;
                RDouble refflowV;
                RDouble refflowW;
                bcData->GetData("initU", &refflowU, PHDOUBLE, 1);
                bcData->GetData("initV", &refflowV, PHDOUBLE, 1);
                bcData->GetData("initW", &refflowW, PHDOUBLE, 1);

                for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];
                    if (FaceFlux[iFacelocal] < 0)
                    {
                        u[re] = refflowU;
                        v[re] = refflowV;
                        w[re] = refflowW;
                        RDouble muLocal = mu[le];

                        RDouble fluxbValue = rho[re] * area[iFacelocal] * (u[re] * xfn[iFacelocal] + v[re] * yfn[iFacelocal] + w[re] * zfn[iFacelocal]);
                        RDouble fluxbCorrectValue = FaceFlux[iFacelocal] / rho[le] * cRho[re] * pc[re];
                        FaceFlux[iFacelocal] = fluxbValue + fluxbCorrectValue;
                    }
                    else
                    {
                        RDouble muLocal = mu[le];

                        RDouble dx = xfc[iFacelocal] - xcc[le];
                        RDouble dy = yfc[iFacelocal] - ycc[le];
                        RDouble dz = zfc[iFacelocal] - zcc[le];

                        RDouble ad = xfn[iFacelocal] * dx + yfn[iFacelocal] * dy + zfn[iFacelocal] * dz;

                        RDouble ub = u[le] + InverseAp[le] * xfn[iFacelocal]  / ad *
                            (p[re] - p[le] - dpdx[le] * dx - dpdy[le] * dy - dpdz[le] * dz);
                        RDouble vb = v[le] + InverseAp[le] * yfn[iFacelocal]  / ad *
                            (p[re] - p[le] - dpdx[le] * dx - dpdy[le] * dy - dpdz[le] * dz);
                        RDouble wb = w[le] + InverseAp[le] * zfn[iFacelocal]  / ad *
                            (p[re] - p[le] - dpdx[le] * dx - dpdy[le] * dy - dpdz[le] * dz);

                        RDouble fluxbValue = rho[le] * 
                            (area[iFacelocal] * (ub * xfn[iFacelocal]  + vb * yfn[iFacelocal] + wb * zfn[iFacelocal]) 
                                + FaceFluxRelaxCoeff * dun[iFacelocal]);
                        RDouble fluxbCorrectValue = FaceFlux[iFacelocal] / rho[le] * cRho[re] * pc[le];
                        FaceFlux[iFacelocal] = fluxbValue + fluxbCorrectValue;
                    }
                }
            }
            else if (bcType == PHENGLEI::INFLOW)
            {
                for (size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];
                    FaceFlux[iFacelocal] = rho[re] * area[iFacelocal] * (u[re] * xfn[iFacelocal] + v[re] * yfn[iFacelocal] + w[re] * zfn[iFacelocal]);
                }

                if (cRho != nullptr)
                {
                    for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                    {
                        int iFacelocal = (*faceIndex)[iFace];
                        int re = rightCellOfFace[iFacelocal];
                        FaceFlux[iFacelocal] += FaceFlux[iFacelocal] / rho[re] * cRho[re] * pc[re];
                    }
                }
            }
            else if (bcType == PHENGLEI::OUTFLOW)
            {

                for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];

                    RDouble muLocal = mu[le];

                    RDouble dx = xfc[iFacelocal] - xcc[le];
                    RDouble dy = yfc[iFacelocal] - ycc[le];
                    RDouble dz = zfc[iFacelocal] - zcc[le];

                    RDouble ad = xfn[iFacelocal] * dx + yfn[iFacelocal] * dy + zfn[iFacelocal] * dz;

                    RDouble ub = u[le] + InverseAp[le] * xfn[iFacelocal] / ad *
                        ( p[re] - p[le] - dpdx[le] * dx - dpdy[le] * dy - dpdz[le] * dz );
                    RDouble vb = v[le] + InverseAp[le] * yfn[iFacelocal]  / ad *
                        ( p[re] - p[le] - dpdx[le] * dx - dpdy[le] * dy - dpdz[le] * dz );
                    RDouble wb = w[le] + InverseAp[le] * zfn[iFacelocal]  / ad *
                        ( p[re] - p[le] - dpdx[le] * dx - dpdy[le] * dy - dpdz[le] * dz );

                    FaceFlux[iFacelocal] = rho[re] * area[iFacelocal] * 
                        ((ub * xfn[iFacelocal]  + vb * yfn[iFacelocal] + wb * zfn[iFacelocal])+ FaceFluxRelaxCoeff * dun[iFacelocal]);
                }

                if (cRho != nullptr)
                {
                    for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                    {
                        int iFacelocal = (*faceIndex)[iFace];
                        int re = rightCellOfFace[iFacelocal];
                        FaceFlux[iFacelocal] += FaceFlux[iFacelocal] / rho[re] * cRho[re] * pc[re];
                    }
                }
            }
            else if (bcType == PHENGLEI::PRESSURE_INLET)
            {
                for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];

                    RDouble dx = xfc[iFacelocal] - xcc[le];
                    RDouble dy = yfc[iFacelocal] - ycc[le];
                    RDouble dz = zfc[iFacelocal] - zcc[le];

                    if (p[re] > ptotal[re])
                    {
                        p[re] = ptotal[re];
                    }

                    u[re] = u[le];
                    v[re] = v[le];
                    w[re] = w[le];
                    FaceFlux[iFacelocal] = rho[re] * area[iFacelocal] * (u[re] * xfn[iFacelocal] + v[re] * yfn[iFacelocal] + w[re] * zfn[iFacelocal]);
                }
            }
            else if (bcType == PHENGLEI::PRESSURE_OUTLET)
            {
                for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];

                    RDouble dx = xfc[iFacelocal] - xcc[le];
                    RDouble dy = yfc[iFacelocal] - ycc[le];
                    RDouble dz = zfc[iFacelocal] - zcc[le];

                    RDouble ad = xfn[iFacelocal] * area[iFacelocal] * dx + yfn[iFacelocal] * area[iFacelocal] * dy + zfn[iFacelocal] * area[iFacelocal] * dz;

                    RDouble ub = u[le] + InverseAp[le] * xfn[iFacelocal] * area[iFacelocal] / ad *
                        ( p[re] - p[le] - dpdx[le] * dx - dpdy[le] * dy - dpdz[le] * dz );
                    RDouble vb = v[le] + InverseAp[le] * yfn[iFacelocal] * area[iFacelocal] / ad *
                        ( p[re] - p[le] - dpdx[le] * dx - dpdy[le] * dy - dpdz[le] * dz );
                    RDouble wb = w[le] + InverseAp[le] * zfn[iFacelocal] * area[iFacelocal] / ad *
                        ( p[re] - p[le] - dpdx[le] * dx - dpdy[le] * dy - dpdz[le] * dz );

                    FaceFlux[iFacelocal] = rho[re] * 
                        ((ub * xfn[iFacelocal] * area[iFacelocal] + vb * yfn[iFacelocal] * area[iFacelocal] + wb * zfn[iFacelocal] * area[iFacelocal])
                            + FaceFluxRelaxCoeff * dun[iFacelocal]);
                }

                if (cRho != nullptr)
                {
                    for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                    {
                        int iFacelocal = (*faceIndex)[iFace];
                        int le = leftCellOfFace[iFacelocal];
                        int re = rightCellOfFace[iFacelocal];
                        FaceFlux[iFacelocal] += FaceFlux[iFacelocal] / rho[re] * cRho[re] * pc[re];
                    }
                }
            }
            else if (bcType == PHENGLEI::INTERFACE)
            {
                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    int iFace = *iter;
                    int le = leftCellOfFace[iFace];
                    int re = rightCellOfFace[iFace];

                    RDouble ax = area[iFace] * xfn[iFace];
                    RDouble ay = area[iFace] * yfn[iFace];
                    RDouble az = area[iFace] * zfn[iFace];

                    RDouble dx = xcc[re] - xcc[le];
                    RDouble dy = ycc[re] - ycc[le];
                    RDouble dz = zcc[re] - zcc[le];

                    RDouble ad = ax * dx + ay * dy + az * dz;
                    RDouble a2 = ax * ax + ay * ay + az * az;

                    RDouble ubar_a = (faceWeightOfLeftCell[iFace] * u[le] + (1.0 - faceWeightOfLeftCell[iFace]) * u[re]) * ax 
                        + (faceWeightOfLeftCell[iFace] * v[le] + (1.0 - faceWeightOfLeftCell[iFace]) * v[re]) * ay 
                        + (faceWeightOfLeftCell[iFace] * w[le] + (1.0 - faceWeightOfLeftCell[iFace]) * w[re]) * az;

                    RDouble pDiff;
                    RDouble pdf;

                    int bodyForceFlag = GlobalDataBase::GetIntParaFromDB("bodyForceFlag");
                    if (bodyForceFlag != 1)
                    {
                        pDiff = (InverseAp[le] * faceWeightOfLeftCell[iFace] + InverseAp[re] * (1.0 - faceWeightOfLeftCell[iFace])) *( p[re] - p[le]) ;
                        pdf = ((InverseAp[le] * faceWeightOfLeftCell[iFace] * dpdx[le] + InverseAp[re] * (1.0 - faceWeightOfLeftCell[iFace]) * dpdx[re]) * dx +
                            (InverseAp[le] * faceWeightOfLeftCell[iFace] * dpdy[le] + InverseAp[re] * (1.0 - faceWeightOfLeftCell[iFace]) * dpdy[re]) * dy+
                            (InverseAp[le] * faceWeightOfLeftCell[iFace] * dpdz[le] + InverseAp[re] * (1.0 - faceWeightOfLeftCell[iFace]) * dpdz[re]) * dz);
                    }
                    else
                    {
                        RDouble dxfL = xfc[iFace] - xcc[le];
                        RDouble dyfL = yfc[iFace] - ycc[le];
                        RDouble dzfL = zfc[iFace] - zcc[le];
                        RDouble dxfR = xfc[iFace] - xcc[re];
                        RDouble dyfR = yfc[iFace] - ycc[re];
                        RDouble dzfR = zfc[iFace] - zcc[re];

                        pDiff = (InverseAp[le] * faceWeightOfLeftCell[iFace] + InverseAp[re] * (1.0 - faceWeightOfLeftCell[iFace])) * (p[re] - p[le] +
                            (bfU[re] * dxfR + bfV[re] * dyfR + bfW[re] * dzfR)- (bfU[le] * dxfL + bfV[le] * dyfL + bfW[le] * dzfL));
                        pdf = ((InverseAp[le] * faceWeightOfLeftCell[iFace] * (dpdx[le] - bfU[le]) + InverseAp[re] * (1.0 - faceWeightOfLeftCell[iFace]) * (dpdx[re] - bfU[re])) * dx +
                            (InverseAp[le] * faceWeightOfLeftCell[iFace] * (dpdy[le] - bfV[le]) + InverseAp[re] * (1.0 - faceWeightOfLeftCell[iFace]) * (dpdy[re] - bfV[re])) * dy +
                            (InverseAp[le] * faceWeightOfLeftCell[iFace] * (dpdz[le] - bfW[le]) + InverseAp[re] * (1.0 - faceWeightOfLeftCell[iFace]) * (dpdz[re] - bfW[re]) * dz));

                    };

                    pDiff *= a2/ad;
                    pdf*= a2 / ad;

                    FaceFlux[iFace] = rhoFace[iFace] * (ubar_a - pDiff + pdf + FaceFluxRelaxCoeff * dun[iFace]);
                }
            }
        }
    }

    void MomEqCalculator::updateBCValue(Grid*  gridIn)
    {
        UnstructGrid* grid = UnstructGridCast(gridIn);
        int nTotalCell = grid->GetNTotalCell();
        int nIFace = grid->GetNIFace();
        int nTotalFace = grid->GetNTotalFace();
        int nBoundFace = grid->GetNBoundFace();
        int nTotalInFace = nTotalFace - nBoundFace;
        RDouble *pc = reinterpret_cast<RDouble *> (grid->GetDataPtr("pc"));
        RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("diagMatrixCoeff"));
        RDouble *bCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("bCoeff"));

        int *leftCellOfFace = grid->GetLeftCellOfFace();
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
        RDouble *u = reinterpret_cast<RDouble *> (grid->GetDataPtr("U"));
        RDouble *v = reinterpret_cast<RDouble *> (grid->GetDataPtr("V"));
        RDouble *w = reinterpret_cast<RDouble *> (grid->GetDataPtr("W"));
        RDouble *p = reinterpret_cast<RDouble *> (grid->GetDataPtr("P"));
        RDouble *InverseAp = reinterpret_cast<RDouble *> (grid->GetDataPtr("InverseAp"));
        RDouble *bfU = reinterpret_cast<RDouble *> (grid->GetDataPtr("bfU"));
        RDouble *bfV = reinterpret_cast<RDouble *> (grid->GetDataPtr("bfV"));
        RDouble *bfW = reinterpret_cast<RDouble *> (grid->GetDataPtr("bfW"));
        RDouble *dpdx = reinterpret_cast<RDouble *> (grid->GetDataPtr("dpdx"));
        RDouble *dpdy = reinterpret_cast<RDouble *> (grid->GetDataPtr("dpdy"));
        RDouble *dpdz = reinterpret_cast<RDouble *> (grid->GetDataPtr("dpdz"));
        RDouble* rho  = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
        RDouble* mu   = reinterpret_cast<RDouble*>(grid->GetDataPtr("mu"));
        RDouble* FaceFlux = reinterpret_cast<RDouble*>(grid->GetDataPtr("FaceFlux"));
        RDouble* dun  = reinterpret_cast<RDouble*>(grid->GetDataPtr("dun"));
        RDouble *faceWeightOfLeftCell = reinterpret_cast<RDouble *>(grid->GetDataPtr("faceWeightOfLeftCell"));

        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

        for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
        {
            UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
            int bcType = bcRegion->GetBCType();
            vector<int>* faceIndex = bcRegion->GetFaceIndex();
            Data_Param* bcData = bcRegion->GetBCParamDataBase();

            if (bcType == PHENGLEI::FARFIELD)
            {
                RDouble refflowU;
                RDouble refflowV;
                RDouble refflowW;
                bcData->GetData("initU", &refflowU, PHDOUBLE, 1);
                bcData->GetData("initV", &refflowV, PHDOUBLE, 1);
                bcData->GetData("initW", &refflowW, PHDOUBLE, 1);

                for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];
                    if (FaceFlux[iFacelocal] < 0)
                    {
                        u[re] = refflowU;
                        v[re] = refflowV;
                        w[re] = refflowW;
                    }
                    else
                    {
                        RDouble muLocal = mu[le];

                        RDouble dx = xfc[iFacelocal] - xcc[le];
                        RDouble dy = yfc[iFacelocal] - ycc[le];
                        RDouble dz = zfc[iFacelocal] - zcc[le];

                        RDouble ad = xfn[iFacelocal] * area[iFacelocal] * dx + yfn[iFacelocal] * area[iFacelocal] * dy + zfn[iFacelocal] * area[iFacelocal] * dz;

                        RDouble tempTerm = (p[re] - p[le] - dpdx[le] * dx - dpdy[le] * dy - dpdz[le] * dz);
                        u[re] = u[le] + InverseAp[le] * xfn[iFacelocal] * area[iFacelocal] / ad * tempTerm;
                        v[re] = v[le] + InverseAp[le] * yfn[iFacelocal] * area[iFacelocal] / ad * tempTerm;
                        w[re] = w[le] + InverseAp[le] * zfn[iFacelocal] * area[iFacelocal] / ad * tempTerm;
                    }
                }
            }
            else if (bcType == PHENGLEI::INFLOW)
            {
            }
            else if (bcType == PHENGLEI::SOLID_SURFACE)
            {
                RDouble ub = 0.0;
                if (bcData->IsExist("initU", PHDOUBLE, 1))
                {
                    bcData->GetData("initU", &ub, PHDOUBLE, 1);
                }

                RDouble vb = 0.0;
                if (bcData->IsExist("initV", PHDOUBLE, 1))
                {
                    bcData->GetData("initV", &vb, PHDOUBLE, 1);
                }

                RDouble wb = 0.0;
                if (bcData->IsExist("initW", PHDOUBLE, 1))
                {
                    bcData->GetData("initW", &wb, PHDOUBLE, 1);
                }

                for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];

                    u[re] = RDouble(ub);
                    v[re] = RDouble(vb);
                    w[re] = RDouble(wb);
                }
            }
            else if (bcType == PHENGLEI::OUTFLOW)
            {
                for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];

                    RDouble muLocal = mu[le];

                    RDouble dx = xfc[iFacelocal] - xcc[le];
                    RDouble dy = yfc[iFacelocal] - ycc[le];
                    RDouble dz = zfc[iFacelocal] - zcc[le];

                    RDouble ad = xfn[iFacelocal] * area[iFacelocal] * dx + yfn[iFacelocal] * area[iFacelocal] * dy + zfn[iFacelocal] * area[iFacelocal] * dz;

                    RDouble tempTerm = (p[re] - p[le] - dpdx[le] * dx - dpdy[le] * dy - dpdz[le] * dz);
                    u[re] = u[le] + InverseAp[le] * xfn[iFacelocal] * area[iFacelocal] / ad * tempTerm;
                    v[re] = v[le] + InverseAp[le] * yfn[iFacelocal] * area[iFacelocal] / ad * tempTerm;
                    w[re] = w[le] + InverseAp[le] * zfn[iFacelocal] * area[iFacelocal] / ad * tempTerm;
                }
            }
            else if (bcType == PHENGLEI::PRESSURE_INLET)
            {
                for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];

                    u[re] = u[le];
                    v[re] = v[le];
                    w[re] = w[le];

                }
            }
            else if (bcType == PHENGLEI::PRESSURE_OUTLET)
            {
                for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];

                    RDouble muLocal = mu[le];

                    RDouble dx = xfc[iFacelocal] - xcc[le];
                    RDouble dy = yfc[iFacelocal] - ycc[le];
                    RDouble dz = zfc[iFacelocal] - zcc[le];

                    RDouble ad = xfn[iFacelocal] * area[iFacelocal] * dx + yfn[iFacelocal] * area[iFacelocal] * dy + zfn[iFacelocal] * area[iFacelocal] * dz;

                    RDouble tempTerm = (p[re] - p[le] - dpdx[le] * dx - dpdy[le] * dy - dpdz[le] * dz);
                    u[re] = u[le] + InverseAp[le] * xfn[iFacelocal] * area[iFacelocal] / ad * tempTerm;
                    v[re] = v[le] + InverseAp[le] * yfn[iFacelocal] * area[iFacelocal] / ad * tempTerm;
                    w[re] = w[le] + InverseAp[le] * zfn[iFacelocal] * area[iFacelocal] / ad * tempTerm;
                }
            }
        }

        CommunicateAnInterfaceVar(u);
        CommunicateAnInterfaceVar(v);
        CommunicateAnInterfaceVar(w);
    }

    void MomEqCalculator::CalcVelocityRes(UnstructGrid *grid)
    {
        int nTotalCell = grid->GetNTotalCell();
        RDouble *scU = reinterpret_cast<RDouble *> (grid->GetDataPtr("scU"));
        RDouble *scV = reinterpret_cast<RDouble *> (grid->GetDataPtr("scV"));
        RDouble *scW = reinterpret_cast<RDouble *> (grid->GetDataPtr("scW"));

        RDouble uResNow = 0.0;
        RDouble vResNow = 0.0;
        RDouble wResNow = 0.0;

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            uResNow += std::abs(scU[iCell]);
            vResNow += std::abs(scV[iCell]);
            wResNow += std::abs(scW[iCell]);
        }

        RDouble temp = uResNow;
        PH_AllreduceSepMode(&temp, &uResNow, 1, MPI_DOUBLE, MPI_SUM);

        temp = vResNow;
        PH_AllreduceSepMode(&temp, &vResNow, 1, MPI_DOUBLE, MPI_SUM);

        temp = wResNow;
        PH_AllreduceSepMode(&temp, &wResNow, 1, MPI_DOUBLE, MPI_SUM);

        grid->UpdateData("UResNow", &uResNow, PHDOUBLE, 1);
        grid->UpdateData("VResNow", &vResNow, PHDOUBLE, 1);
        grid->UpdateData("WResNow", &wResNow, PHDOUBLE, 1);
    }

    void MomEqCalculator::CoefPrepForPres(Grid*  gridIn)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotalFace = grid->GetNTotalFace();
        int *leftCellOfFace = grid->GetLeftCellOfFace();
        int *rightCellOfFace = grid->GetRightCellOfFace();
        RDouble* vol = grid->GetCellVolume();
        RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("diagMatrixCoeff"));
        RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("upperMatrixCoeff"));
        RDouble *lowerMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("lowerMatrixCoeff"));
        RDouble *hu = reinterpret_cast<RDouble *> (grid->GetDataPtr("hu"));
        RDouble *InverseAp = reinterpret_cast<RDouble *> (grid->GetDataPtr("InverseAp"));

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            hu[iCell] = diagMatrixCoeff[iCell];
        }

        int isSIMPLEC = GlobalDataBase::GetIntParaFromDB("SIMPLEC");
        if (isSIMPLEC)
        {
            for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
            {
                int le = leftCellOfFace[iFace];
                int re = rightCellOfFace[iFace];
                hu[le] += upperMatrixCoeff[iFace];
                hu[re] += lowerMatrixCoeff[iFace];
            }
            UnstructBCSet **bcr = grid->GetBCRecord();
            for (int iFace = 0; iFace < nBoundFace; ++iFace)
            {
                int le = leftCellOfFace[iFace];
                int re = rightCellOfFace[iFace];
                int bcType = bcr[iFace]->GetKey();
                if (bcType == PHENGLEI::INTERFACE)
                {
                    hu[le] += upperMatrixCoeff[iFace];
                }
            }
        }

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            hu[iCell] = vol[iCell] / hu[iCell];
            InverseAp[iCell] = vol[iCell] / diagMatrixCoeff[iCell];
        }

        string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
        if (flowSolverName == "PISO")
        {
            SetField(hu, InverseAp, nTotalCell);
        }

        CommunicateAnInterfaceVar(hu);
        CommunicateAnInterfaceVar(InverseAp);

    }

    void MomEqCalculator:: CorrectFaceFlux(Grid*  gridIn)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotalFace = grid->GetNTotalFace();
        int *leftCellOfFace = grid->GetLeftCellOfFace();
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
        RDouble *u = reinterpret_cast<RDouble *> (grid->GetDataPtr("U"));
        RDouble *v = reinterpret_cast<RDouble *> (grid->GetDataPtr("V"));
        RDouble *w = reinterpret_cast<RDouble *> (grid->GetDataPtr("W"));
        RDouble *p = reinterpret_cast<RDouble *> (grid->GetDataPtr("P"));
        RDouble *pc = reinterpret_cast<RDouble *> (grid->GetDataPtr("pc"));
        RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("upperMatrixCoeff"));
        RDouble* rho  = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
        RDouble* mu   = reinterpret_cast<RDouble*>(grid->GetDataPtr("mu"));
        RDouble* FaceFlux = reinterpret_cast<RDouble*>(grid->GetDataPtr("FaceFlux"));
        RDouble* dun  = reinterpret_cast<RDouble*>(grid->GetDataPtr("dun"));
        RDouble *rhoFace = reinterpret_cast<RDouble*>(grid->GetDataPtr("rhoFace"));
        RDouble *faceWeightOfLeftCell = reinterpret_cast<RDouble *>(grid->GetDataPtr("faceWeightOfLeftCell"));

        RDouble* cellCenterX = grid->GetCellCenterX();
        RDouble* cellCenterY = grid->GetCellCenterY();
        RDouble* cellCenterZ = grid->GetCellCenterZ();
        RDouble *hu = reinterpret_cast<RDouble *> (grid->GetDataPtr("hu"));
        RDouble *cRho = reinterpret_cast<RDouble*>(grid->GetDataPtr("cRho"));

        string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
        if (flowSolverName == "CompressibleSIMPLE")
        {
            for (int iFace = nBoundFace; iFace < nTotalFace; ++iFace)
            {
                int le = leftCellOfFace[iFace];
                int re = rightCellOfFace[iFace];

                RDouble ax = area[iFace] * xfn[iFace];
                RDouble ay = area[iFace] * yfn[iFace];
                RDouble az = area[iFace] * zfn[iFace];
                RDouble cRhof = faceWeightOfLeftCell[iFace] * cRho[le] + (1.0 - faceWeightOfLeftCell[iFace]) * cRho[re];
                RDouble pcf = faceWeightOfLeftCell[iFace] * pc[le] + (1.0 - faceWeightOfLeftCell[iFace]) * pc[re];
                RDouble dx = cellCenterX[re] - cellCenterX[le];
                RDouble dy = cellCenterY[re] - cellCenterY[le];
                RDouble dz = cellCenterZ[re] - cellCenterZ[le];
                RDouble ad = ax * dx + ay * dy + az * dz;
                RDouble a2 = ax * ax + ay * ay + az * az;
                RDouble uf = faceWeightOfLeftCell[iFace] * u[le] + (1.0 - faceWeightOfLeftCell[iFace]) * u[re];
                RDouble vf = faceWeightOfLeftCell[iFace] * v[le] + (1.0 - faceWeightOfLeftCell[iFace]) * v[re];
                RDouble wf = faceWeightOfLeftCell[iFace] * w[le] + (1.0 - faceWeightOfLeftCell[iFace]) * w[re];
                RDouble un = uf * ax + vf * ay + wf * az;

                if (FaceFlux[iFace] > 0)
                {
                    rhoFace[iFace] = rho[le];
                }
                else
                {
                    rhoFace[iFace] = rho[re];
                }

                FaceFlux[iFace] += FaceFlux[iFace] * cRhof * pcf / rhoFace[iFace];

                FaceFlux[iFace] -= rhoFace[iFace] * 
                    ((hu[le] * faceWeightOfLeftCell[iFace] + hu[re] * (1.0 - faceWeightOfLeftCell[iFace])) * a2) / ad * (pc[re] - pc[le]);
                dun[iFace] = FaceFlux[iFace] / rhoFace[iFace] - un;
            }
        }
        else
        {
            for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
            {
                int le = leftCellOfFace[iFace];
                int re = rightCellOfFace[iFace];

                RDouble ax = area[iFace] * xfn[iFace];
                RDouble ay = area[iFace] * yfn[iFace];
                RDouble az = area[iFace] * zfn[iFace];
                RDouble uf = faceWeightOfLeftCell[iFace] * u[le] + (1.0 - faceWeightOfLeftCell[iFace]) * u[re];
                RDouble vf = faceWeightOfLeftCell[iFace] * v[le] + (1.0 - faceWeightOfLeftCell[iFace]) * v[re];
                RDouble wf = faceWeightOfLeftCell[iFace] * w[le] + (1.0 - faceWeightOfLeftCell[iFace]) * w[re];
                RDouble un = uf * ax + vf * ay + wf * az;

                FaceFlux[iFace] += upperMatrixCoeff[iFace] * (pc[re] - pc[le]);
                dun[iFace] = FaceFlux[iFace] / rhoFace[iFace] - un;
            }
        }

        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

        for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
        {
            UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
            int bcType = bcRegion->GetBCType();
            vector<int>* faceIndex = bcRegion->GetFaceIndex();
            Data_Param* bcData = bcRegion->GetBCParamDataBase();


            if (bcType == PHENGLEI::FARFIELD)
            {
                RDouble un;

                RDouble refflowU;
                RDouble refflowV;
                RDouble refflowW;
                bcData->GetData("initU", &refflowU, PHDOUBLE, 1);
                bcData->GetData("initV", &refflowV, PHDOUBLE, 1);
                bcData->GetData("initW", &refflowW, PHDOUBLE, 1);

                for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];

                    if (FaceFlux[iFacelocal] < 0)
                    {
                        FaceFlux[iFacelocal] = 0;
                        RDouble fluxbCorrectValue = rho[re] * (u[re] * xfn[iFacelocal] * area[iFacelocal] + v[re] * yfn[iFacelocal] * area[iFacelocal] + w[re] * zfn[iFacelocal] * area[iFacelocal]);
                        if (fluxbCorrectValue < 0)
                        {
                            FaceFlux[iFacelocal] = fluxbCorrectValue;
                        }
                    }
                    else
                    {
                        RDouble fluxbValue = rho[re] * area[iFacelocal] * (u[re] * xfn[iFacelocal] + v[re] * yfn[iFacelocal] + w[re] * zfn[iFacelocal]);
                        FaceFlux[iFacelocal] = fluxbValue;
                        un = u[re] * xfn[iFacelocal] * area[iFacelocal] + v[re] * yfn[iFacelocal] * area[iFacelocal] + w[re] * zfn[iFacelocal] * area[iFacelocal];
                        dun[iFacelocal] = fluxbValue / rho[le] - un;
                        FaceFlux[iFacelocal] = fluxbValue;
                        un = u[re] * xfn[iFacelocal] * area[iFacelocal] + v[re] * yfn[iFacelocal] * area[iFacelocal] + w[re] * zfn[iFacelocal] * area[iFacelocal];
                        dun[iFacelocal] = fluxbValue / rho[le] - un;
                        RDouble fluxbCorrectValue = FaceFlux[iFacelocal] / rho[le] * cRho[re] * pc[le];
                        FaceFlux[iFacelocal] += fluxbCorrectValue;
                        if (FaceFlux[iFacelocal] < 0)
                        {
                            FaceFlux[iFacelocal] = 0;

                            RDouble fluxbCorrectValue = rho[re] * area[iFacelocal] * (u[re] * xfn[iFacelocal]  + v[re] * yfn[iFacelocal] + w[re] * zfn[iFacelocal]);

                            if (fluxbCorrectValue < 0)
                            {
                                FaceFlux[iFacelocal] = fluxbCorrectValue;
                            }
                        }
                    }
                }
            }
            else if (bcType == PHENGLEI::INFLOW)
            {
                if (cRho != nullptr)
                {
                    for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                    {
                        int iFacelocal = (*faceIndex)[iFace];
                        int re = rightCellOfFace[iFacelocal];
                        FaceFlux[iFacelocal] += FaceFlux[iFacelocal] / rho[re] * cRho[re] * pc[re];
                    }
                }
            }
            else if (bcType == PHENGLEI::OUTFLOW)
            {
                for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];


                    FaceFlux[iFacelocal] = rho[re] * (u[re] * xfn[iFacelocal] * area[iFacelocal] + v[re] * yfn[iFacelocal] * area[iFacelocal] + w[re] * zfn[iFacelocal] * area[iFacelocal]);

                    RDouble un = u[re] * xfn[iFacelocal] * area[iFacelocal] + v[re] * yfn[iFacelocal] * area[iFacelocal] + w[re] * zfn[iFacelocal] * area[iFacelocal];
                    dun[iFacelocal] = rho[re] * (u[re] * xfn[iFacelocal] * area[iFacelocal] + v[re] * yfn[iFacelocal] * area[iFacelocal] + w[re] * zfn[iFacelocal] * area[iFacelocal]) / rho[re] - un;
                }

                RDouble* cRho = reinterpret_cast<RDouble*>(grid->GetDataPtr("cRho"));
                if (cRho != nullptr)
                {
                    RDouble* pc = reinterpret_cast<RDouble*>(grid->GetDataPtr("pc"));
                    for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                    {
                        int iFacelocal = (*faceIndex)[iFace];
                        int re = rightCellOfFace[iFacelocal];
                        FaceFlux[iFacelocal] += FaceFlux[iFacelocal] / rho[re] * cRho[re] * pc[re];
                    }
                }
            }
            else if (bcType == PHENGLEI::PRESSURE_INLET)
            {
                for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];

                    FaceFlux[iFacelocal] = rho[re] * (u[re] * xfn[iFacelocal] * area[iFacelocal] + v[re] * yfn[iFacelocal] * area[iFacelocal] + w[re] * zfn[iFacelocal] * area[iFacelocal]);
                }
            }
            else if (bcType == PHENGLEI::PRESSURE_OUTLET)
            {
                for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];

                    RDouble fluxbValue = rho[re] * area[iFacelocal] * (u[re] * xfn[iFacelocal] + v[re] * yfn[iFacelocal] + w[re] * zfn[iFacelocal]);

                    FaceFlux[iFacelocal] = fluxbValue;

                    RDouble un = u[re] * xfn[iFacelocal] * area[iFacelocal] + v[re] * yfn[iFacelocal] * area[iFacelocal] + w[re] * zfn[iFacelocal] * area[iFacelocal];
                    dun[iFacelocal] = fluxbValue / rho[re] - un;
                }

                RDouble* cRho = reinterpret_cast<RDouble*>(grid->GetDataPtr("cRho"));
                if (cRho != nullptr)
                {
                    RDouble* pc = reinterpret_cast<RDouble*>(grid->GetDataPtr("pc"));
                    for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                    {
                        int iFacelocal = (*faceIndex)[iFace];
                        int le = leftCellOfFace[iFacelocal];
                        int re = rightCellOfFace[iFacelocal];
                        FaceFlux[iFacelocal] += FaceFlux[iFacelocal] / rho[re] * cRho[re] * pc[re];
                    }
                }
            }
            else if (bcType == PHENGLEI::INTERFACE)
            {
                for (std::size_t iFaceLocal = 0; iFaceLocal < faceIndex->size(); ++iFaceLocal)
                {
                    int iFace = (*faceIndex)[iFaceLocal];
                    int le = leftCellOfFace[iFace];
                    int re = rightCellOfFace[iFace];
                    RDouble ax = area[iFace] * xfn[iFace];
                    RDouble ay = area[iFace] * yfn[iFace];
                    RDouble az = area[iFace] * zfn[iFace];
                    RDouble uf = faceWeightOfLeftCell[iFace] * u[le] + (1.0 - faceWeightOfLeftCell[iFace]) * u[re];
                    RDouble vf = faceWeightOfLeftCell[iFace] * v[le] + (1.0 - faceWeightOfLeftCell[iFace]) * v[re];
                    RDouble wf = faceWeightOfLeftCell[iFace] * w[le] + (1.0 - faceWeightOfLeftCell[iFace]) * w[re];
                    RDouble un = uf * ax + vf * ay + wf * az;
                    if (flowSolverName == "CompressibleSIMPLE")
                    {
                        RDouble cRhof = faceWeightOfLeftCell[iFace] * cRho[le] + (1.0 - faceWeightOfLeftCell[iFace]) * cRho[re];
                        RDouble pcf = faceWeightOfLeftCell[iFace] * pc[le] + (1.0 - faceWeightOfLeftCell[iFace]) * pc[re];
                        RDouble dx = cellCenterX[re] - cellCenterX[le];
                        RDouble dy = cellCenterY[re] - cellCenterY[le];
                        RDouble dz = cellCenterZ[re] - cellCenterZ[le];
                        RDouble ad = ax * dx + ay * dy + az * dz;
                        RDouble a2 = ax * ax + ay * ay + az * az;

                        if (FaceFlux[iFace] > 0)
                        {
                            rhoFace[iFace] = rho[le];
                        }
                        else
                        {
                            rhoFace[iFace] = rho[re];
                        }

                        FaceFlux[iFace] += FaceFlux[iFace] * cRhof * pcf / rhoFace[iFace];

                        FaceFlux[iFace] -= rhoFace[iFace] *
                            ((hu[le] * faceWeightOfLeftCell[iFace] + hu[re] * (1.0 - faceWeightOfLeftCell[iFace])) * a2) / ad * (pc[re] - pc[le]);
                        dun[iFace] = FaceFlux[iFace] / rhoFace[iFace] - un;
                    }
                    else
                    {
                        FaceFlux[iFace] += upperMatrixCoeff[iFace] * (pc[re] - pc[le]);
                        dun[iFace] = FaceFlux[iFace] / rhoFace[iFace] - un;
                    }
                }
            }
        }
    }

    void MomEqCalculator:: pressureSourceTerm(Grid*  gridIn, int iVariable)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);

        int nTotalCell = grid->GetNTotalCell();
        RDouble *vol = grid->GetCellVolume();
        RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));
        RDouble *dpdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdx"));
        RDouble *dpdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdy"));
        RDouble *dpdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdz"));
        RDouble* gradp;

        switch (iVariable)
        {
            case 0:
                gradp = dpdx;
                break;
            case 1:
                gradp = dpdy;
                break;
            case 2:
                gradp = dpdz;
                break;
        }

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            bCoeff[iCell] -= gradp[iCell] * vol[iCell];
        }
    }

    void MomEqCalculator:: BulkVisSourceTerm(Grid*  gridIn, int iVariable)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);
        int nTotalCell = grid->GetNTotalCell();
        int nTotalFace = grid->GetNTotalFace();
        int nBoundFace = grid->GetNBoundFace();

        int* leftCellOfFace = grid->GetLeftCellOfFace();
        int* rightCellOfFace = grid->GetRightCellOfFace();
        RDouble* area = grid->GetFaceArea();
        RDouble* faceNormalX = grid->GetFaceNormalX();
        RDouble* faceNormalY = grid->GetFaceNormalY();
        RDouble* faceNormalZ = grid->GetFaceNormalZ();

        RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));
        RDouble* faceWeightOfLeftCell = reinterpret_cast<RDouble*>(grid->GetDataPtr("faceWeightOfLeftCell"));
        RDouble* mu = reinterpret_cast<RDouble*>(grid->GetDataPtr("mu"));
        RDouble *dudx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdx"));
        RDouble *dvdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdy"));
        RDouble *dwdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdz"));

        RDouble* gradq;
        RDouble* faceNormal;

        switch (iVariable)
        {
            case IDX::S_IU:
                gradq = dudx;
                faceNormal = faceNormalX;
                break;
            case IDX::S_IV:
                gradq = dvdy;
                faceNormal = faceNormalY;
                break;
            case IDX::S_IW:
                gradq = dwdz;
                faceNormal = faceNormalZ;
                break;
        }

        for (int iFace = nBoundFace; iFace < nTotalFace; ++iFace)
        {
            int icL = leftCellOfFace[iFace];
            int icR = rightCellOfFace[iFace];
            RDouble ax = area[iFace] * faceNormal[iFace];
            RDouble muF = faceWeightOfLeftCell[iFace] * mu[icL] + (1.0 - faceWeightOfLeftCell[iFace]) * mu[icR];
            RDouble dudxF = faceWeightOfLeftCell[iFace] * gradq[icL] + (1.0 - faceWeightOfLeftCell[iFace]) * gradq[icR];
            bCoeff[icL] -= (2 / 3) * muF * dudxF * ax;
            bCoeff[icR] -= (2 / 3) * muF * dudxF * ax;
        }

        for (int ibFace = 0; ibFace < nBoundFace; ++ibFace)
        {
            int icL = leftCellOfFace[ibFace];
            int icR = rightCellOfFace[ibFace];
            RDouble muF = mu[icR];
            RDouble ax = area[ibFace] * faceNormal[ibFace];
            RDouble dudxF = gradq[icL];

            bCoeff[icL] -= (2 / 3) * muF * dudxF * ax;
        }
    }

    void MomEqCalculator:: GravitySourceTerm(Grid*  gridIn, int iVariable)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);

        int nTotalCell = grid->GetNTotalCell();
        RDouble *vol = grid->GetCellVolume();
        RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
        RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));

        int isBoussinesqApproximation = GlobalDataBase::GetIntParaFromDB("isBoussinesqApproximation");

        RDouble gx = GlobalDataBase::GetDoubleParaFromDB("gravityX");
        RDouble gy = GlobalDataBase::GetDoubleParaFromDB("gravityY");
        RDouble gz = GlobalDataBase::GetDoubleParaFromDB("gravityZ");

        RDouble gDir;
        if (iVariable ==IDX::S_IU)
        {
            gDir = gx;
        }
        else if (iVariable ==IDX::S_IV)
        {
            gDir = gy;
        }
        else if (iVariable ==IDX::S_IW)
        {
            gDir = gz;
        }

        for (int iCell = 0; iCell < nTotalCell; ++iCell)
        {
            RDouble bouCoe = 1.0;
            if (isBoussinesqApproximation)
            {
                RDouble thermalExpansionCoeff = GlobalDataBase::GetDoubleParaFromDB("thermalExpansionCoeff");  //! coefficient of thermal expansion.
                RDouble refT = GlobalDataBase::GetDoubleParaFromDB("refT");
                RDouble *T = reinterpret_cast<RDouble *>(grid->GetDataPtr("T"));

                bouCoe = -thermalExpansionCoeff * (T[iCell] - refT);
            }
            bCoeff[iCell] += gDir * rho[iCell] * vol[iCell] * bouCoe;
        }
    }

    void MomEqCalculator::CalcOtherbCoeff(Grid* gridIn, int iVariable)
    {
        UnstructGrid* grid = UnstructGridCast(gridIn);

        pressureSourceTerm(grid, iVariable);

        BulkVisSourceTerm(grid, iVariable);

        int isCalGravityVisSource = GlobalDataBase::GetIntParaFromDB("isCalGravityVisSource"); 
        if (isCalGravityVisSource == 1)
        {
            GravitySourceTerm(grid, iVariable);
        }
    }


    void MomEqCalculator::SolveExplicitVelocityEquations(Grid*  gridIn)
    {  
        UnstructGrid *grid = UnstructGridCast(gridIn);
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotalFace = grid->GetNTotalFace();
        int nTotalInnerFace = nTotalFace - nBoundFace;

        int* leftCellOfFace = grid->GetLeftCellOfFace();
        int* rightCellOfFace = grid->GetRightCellOfFace();
        RDouble* vol = grid->GetCellVolume();
        RDouble *uc = reinterpret_cast<RDouble *> (grid->GetDataPtr("uc"));
        RDouble *vc = reinterpret_cast<RDouble *> (grid->GetDataPtr("vc"));
        RDouble *wc = reinterpret_cast<RDouble *> (grid->GetDataPtr("wc"));
        RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("diagMatrixCoeff"));
        RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("upperMatrixCoeff"));
        RDouble *lowerMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("lowerMatrixCoeff"));
        RDouble *dpcdx = reinterpret_cast<RDouble *> (grid->GetDataPtr("dpcdx"));
        RDouble *dpcdy = reinterpret_cast<RDouble *> (grid->GetDataPtr("dpcdy"));
        RDouble *dpcdz = reinterpret_cast<RDouble *> (grid->GetDataPtr("dpcdz"));
        RDouble PPEqRelaxCoeff = GlobalDataBase::GetDoubleParaFromDB("PPEqRelaxCoeff");
        RDouble *u = reinterpret_cast<RDouble *> (grid->GetDataPtr("U"));
        RDouble *v = reinterpret_cast<RDouble *> (grid->GetDataPtr("V"));
        RDouble *w = reinterpret_cast<RDouble *> (grid->GetDataPtr("W"));
        RDouble *p = reinterpret_cast<RDouble *> (grid->GetDataPtr("P"));
        RDouble *InverseAp = reinterpret_cast<RDouble *> (grid->GetDataPtr("InverseAp"));

        CalcGrad(grid);

        constructMatrixACoeff(grid, IDX::S_IU);
        CoefPrepForPres(grid);

        for (int iFace = nBoundFace; iFace < nTotalFace; ++iFace)
        {
            int le = leftCellOfFace [iFace];
            int re = rightCellOfFace[iFace];

            u[le] += (InverseAp[le]  * dpcdx[le] * PPEqRelaxCoeff * upperMatrixCoeff[iFace]) / diagMatrixCoeff[le];
            v[le] += (InverseAp[le]  * dpcdy[le] * PPEqRelaxCoeff * upperMatrixCoeff[iFace]) / diagMatrixCoeff[le];
            w[le] += (InverseAp[le]  * dpcdz[le] * PPEqRelaxCoeff * upperMatrixCoeff[iFace]) / diagMatrixCoeff[le];

            u[re] += (InverseAp[re]  * dpcdx[re] * PPEqRelaxCoeff * lowerMatrixCoeff[iFace]) / diagMatrixCoeff[re];
            v[re] += (InverseAp[re]  * dpcdy[re] * PPEqRelaxCoeff * lowerMatrixCoeff[iFace]) / diagMatrixCoeff[re];
            w[re] += (InverseAp[re]  * dpcdz[re] * PPEqRelaxCoeff * lowerMatrixCoeff[iFace]) / diagMatrixCoeff[re];
        }

        updateBCValue(grid);
    }

    void MomEqCalculator::UpdateUnsteadyFlux(Grid*  gridIn)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);

        string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
        if (flowSolverName == "CompressibleSIMPLE")
        {
            int nTotalCell = grid->GetNTotalCell();
            int nBoundFace = grid->GetNBoundFace();
            int nTotalFace = grid->GetNTotalFace();
            int* leftCellOfFace = grid->GetLeftCellOfFace();
            int* rightCellOfFace = grid->GetRightCellOfFace();
            RDouble* area = grid->GetFaceArea();
            RDouble* xfn = grid->GetFaceNormalX();
            RDouble* yfn = grid->GetFaceNormalY();
            RDouble* zfn = grid->GetFaceNormalZ();
            RDouble* cellCenterX = grid->GetCellCenterX();
            RDouble* cellCenterY = grid->GetCellCenterY();
            RDouble* cellCenterZ = grid->GetCellCenterZ();
            unsigned int nTotal = nTotalCell + nBoundFace;
            RDouble *u = reinterpret_cast<RDouble *> (grid->GetDataPtr("U"));
            RDouble *v = reinterpret_cast<RDouble *> (grid->GetDataPtr("V"));
            RDouble *w = reinterpret_cast<RDouble *> (grid->GetDataPtr("W"));
            RDouble *p = reinterpret_cast<RDouble *> (grid->GetDataPtr("P"));
            RDouble *pc = reinterpret_cast<RDouble *> (grid->GetDataPtr("pc"));
            RDouble *hu = reinterpret_cast<RDouble *> (grid->GetDataPtr("hu"));
            RDouble* rho  = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
            RDouble* mu   = reinterpret_cast<RDouble*>(grid->GetDataPtr("mu"));
            RDouble* FaceFlux = reinterpret_cast<RDouble*>(grid->GetDataPtr("FaceFlux"));
            RDouble* dun  = reinterpret_cast<RDouble*>(grid->GetDataPtr("dun"));
            RDouble *rhoFace = reinterpret_cast<RDouble*>(grid->GetDataPtr("rhoFace"));
            RDouble *cRho = reinterpret_cast<RDouble*>(grid->GetDataPtr("cRho"));
            RDouble *faceWeightOfLeftCell = reinterpret_cast<RDouble*>(grid->GetDataPtr("faceWeightOfLeftCell"));

            for (int iFace = nBoundFace; iFace < nTotalFace; ++iFace)
            {
                int le = leftCellOfFace[iFace];
                int re = rightCellOfFace[iFace];

                RDouble ax = area[iFace] * xfn[iFace];
                RDouble ay = area[iFace] * yfn[iFace];
                RDouble az = area[iFace] * zfn[iFace];

                RDouble cRhof = faceWeightOfLeftCell[iFace] * cRho[le] + (1.0 - faceWeightOfLeftCell[iFace]) * cRho[re];
                RDouble pcf = faceWeightOfLeftCell[iFace] * pc[le] + (1.0 - faceWeightOfLeftCell[iFace]) * pc[re];
                RDouble dx = cellCenterX[re] - cellCenterX[le];
                RDouble dy = cellCenterY[re] - cellCenterY[le];
                RDouble dz = cellCenterZ[re] - cellCenterZ[le];
                RDouble ad = ax * dx + ay * dy + az * dz;
                RDouble a2 = ax * ax + ay * ay + az * az;
                RDouble uFace = faceWeightOfLeftCell[iFace] * u[le] + (1.0 - faceWeightOfLeftCell[iFace]) * u[re];
                RDouble vFace = faceWeightOfLeftCell[iFace] * v[le] + (1.0 - faceWeightOfLeftCell[iFace]) * v[re];
                RDouble wFace = faceWeightOfLeftCell[iFace] * w[le] + (1.0 - faceWeightOfLeftCell[iFace]) * w[re];
                RDouble un = uFace * ax + vFace * ay + wFace * az;

                if (FaceFlux[iFace] > 0)
                {
                    rhoFace[iFace] = rho[le];
                }
                else
                {
                    rhoFace[iFace] = rho[re];
                }

                FaceFlux[iFace] += FaceFlux[iFace] * cRhof * pcf / rhoFace[iFace];

                RDouble aj = rhoFace[iFace] * ((hu[le] * faceWeightOfLeftCell[iFace] + hu[re] * (1.0 - faceWeightOfLeftCell[iFace])) * a2) / ad;
                FaceFlux[iFace] -= aj * (pc[re] - pc[le]);
                dun[iFace] = FaceFlux[iFace] / rhoFace[iFace] - un;
            }

            InterfaceInfo* iinfo = grid->GetInterfaceInfo();
            int nIFace = grid->GetNIFace();

            for (int iFace = 0; iFace < nIFace; ++iFace)
            {
                int* interFace2BoundaryFace = iinfo->GetInterFace2BoundaryFace();
                int boundFace = interFace2BoundaryFace[iFace];
                int le = leftCellOfFace[boundFace];
                int re = rightCellOfFace[boundFace];

                RDouble ax = area[boundFace] * xfn[boundFace];
                RDouble ay = area[boundFace] * yfn[boundFace];
                RDouble az = area[boundFace] * zfn[boundFace];
                RDouble cRhof = faceWeightOfLeftCell[boundFace] * cRho[le] + (1.0 - faceWeightOfLeftCell[boundFace]) * cRho[re];
                RDouble pcf = faceWeightOfLeftCell[boundFace] * pc[le] + (1.0 - faceWeightOfLeftCell[boundFace]) * pc[re];
                RDouble dx = cellCenterX[re] - cellCenterX[le];
                RDouble dy = cellCenterY[re] - cellCenterY[le];
                RDouble dz = cellCenterZ[re] - cellCenterZ[le];
                RDouble ad = ax * dx + ay * dy + az * dz;
                RDouble a2 = ax * ax + ay * ay + az * az;
                RDouble uFace = faceWeightOfLeftCell[boundFace] * u[le] + (1.0 - faceWeightOfLeftCell[boundFace]) * u[re];
                RDouble vFace = faceWeightOfLeftCell[boundFace] * v[le] + (1.0 - faceWeightOfLeftCell[boundFace]) * v[re];
                RDouble wFace = faceWeightOfLeftCell[boundFace] * w[le] + (1.0 - faceWeightOfLeftCell[boundFace]) * w[re];
                RDouble un = uFace * ax + vFace * ay + wFace * az;

                if (FaceFlux[boundFace] > 0)
                {
                    rhoFace[boundFace] = rho[le];
                }
                else
                {
                    rhoFace[boundFace] = rho[re];
                }

                FaceFlux[boundFace] += FaceFlux[boundFace] * cRhof * pcf / rhoFace[boundFace];

                RDouble aj = rhoFace[boundFace] * ((hu[le] * faceWeightOfLeftCell[boundFace] + hu[re] * (1.0 - faceWeightOfLeftCell[boundFace])) * a2) / ad;
                FaceFlux[boundFace] -= aj * (pc[re] - pc[le]);
                dun[boundFace] = FaceFlux[boundFace] / rhoFace[boundFace] - un;
            }
        }
    }

    void MomEqCalculator::IncompressibleInitial(Grid*  gridIn)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);

        string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
        if (flowSolverName == "CompressibleSIMPLE")
        {
            RDouble* Rg = reinterpret_cast<RDouble*>(grid->GetDataPtr("Rg"));
            GAS_SPACE::gas->UpdateRg(grid, Rg);
        }

        RDouble* rho  = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
        RDouble* mu   = reinterpret_cast<RDouble*>(grid->GetDataPtr("mu"));
        RDouble* FaceFlux = reinterpret_cast<RDouble*>(grid->GetDataPtr("FaceFlux"));
        RDouble* dun  = reinterpret_cast<RDouble*>(grid->GetDataPtr("dun"));
        GAS_SPACE::gas->UpdateRho(grid, rho);
        GAS_SPACE::gas->UpdateMu(grid, mu);

        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotalFace = grid->GetNTotalFace();
        int *leftCellOfFace = grid->GetLeftCellOfFace();
        int *rightCellOfFace = grid->GetRightCellOfFace();
        RDouble *area    = grid->GetFaceArea();
        RDouble *xfn = grid->GetFaceNormalX();
        RDouble *yfn = grid->GetFaceNormalY();
        RDouble *zfn = grid->GetFaceNormalZ();
        RDouble *u = reinterpret_cast<RDouble *> (grid->GetDataPtr("U"));
        RDouble *v = reinterpret_cast<RDouble *> (grid->GetDataPtr("V"));
        RDouble *w = reinterpret_cast<RDouble *> (grid->GetDataPtr("W"));
        RDouble *p = reinterpret_cast<RDouble *> (grid->GetDataPtr("P"));
        RDouble *bfU = reinterpret_cast<RDouble *> (grid->GetDataPtr("bfU"));
        RDouble *bfV = reinterpret_cast<RDouble *> (grid->GetDataPtr("bfV"));
        RDouble *bfW = reinterpret_cast<RDouble *> (grid->GetDataPtr("bfW"));

        RDouble *dpdx = reinterpret_cast<RDouble *> (grid->GetDataPtr("dpdx"));
        RDouble *dpdy = reinterpret_cast<RDouble *> (grid->GetDataPtr("dpdy"));
        RDouble *dpdz = reinterpret_cast<RDouble *> (grid->GetDataPtr("dpdz"));
        RDouble *rhoFace = reinterpret_cast<RDouble*>(grid->GetDataPtr("rhoFace"));
        RDouble *faceWeightOfLeftCell = reinterpret_cast<RDouble *>(grid->GetDataPtr("faceWeightOfLeftCell"));

        for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
        {
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];
            rhoFace[iFace] = faceWeightOfLeftCell[iFace] * rho[le] + (1.0 - faceWeightOfLeftCell[iFace]) * rho[re];
            RDouble uFace = faceWeightOfLeftCell[iFace] * u[le] + (1.0 - faceWeightOfLeftCell[iFace]) * u[re];
            RDouble vFace = faceWeightOfLeftCell[iFace] * v[le] + (1.0 - faceWeightOfLeftCell[iFace]) * v[re];
            RDouble wFace = faceWeightOfLeftCell[iFace] * w[le] + (1.0 - faceWeightOfLeftCell[iFace]) * w[re];

            RDouble ax = area[iFace] * xfn[iFace];
            RDouble ay = area[iFace] * yfn[iFace];
            RDouble az = area[iFace] * zfn[iFace];
            FaceFlux[iFace] = rhoFace[iFace] * (uFace * ax + vFace * ay + wFace * az);
            dun[iFace] = 0;
        }

        int bodyForceFlag = GlobalDataBase::GetIntParaFromDB("bodyForceFlag");
        if (bodyForceFlag == 1)
        {
            for (int iCell = 0; iCell < nTotalCell; ++ iCell)
            {
                bfU[iCell] = 0.0;
                bfV[iCell] = 0.0;
                bfW[iCell] = 0.0;
                dpdx[iCell] = 0.0;
                dpdy[iCell] = 0.0;
                dpdz[iCell] = 0.0;
            }
        }

        UpdateBoundaryPressureAndVelocity(grid);

        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];
            RDouble uFace = u[re];
            RDouble vFace = v[re];
            RDouble wFace = w[re];

            RDouble ax = area[iFace] * xfn[iFace];
            RDouble ay = area[iFace] * yfn[iFace];
            RDouble az = area[iFace] * zfn[iFace];
            rhoFace[iFace] = rho[le];
            FaceFlux[iFace] = rho[le] * (uFace * ax + vFace * ay + wFace * az);
            dun[iFace] = 0;
        }

        InterfaceInfo* iinfo = grid->GetInterfaceInfo();
        int nIFace = grid->GetNIFace();

        for (int iFace = 0; iFace < nIFace; ++iFace)
        {
            int* interFace2BoundaryFace = iinfo->GetInterFace2BoundaryFace();
            int boundFace = interFace2BoundaryFace[iFace];

            int le = leftCellOfFace[boundFace];
            int re = rightCellOfFace[boundFace];
            rhoFace[boundFace] = faceWeightOfLeftCell[boundFace] * rho[le] + (1.0 - faceWeightOfLeftCell[boundFace]) * rho[re];
            RDouble uFace = faceWeightOfLeftCell[boundFace] * u[le] + (1.0 - faceWeightOfLeftCell[boundFace]) * u[re];
            RDouble vFace = faceWeightOfLeftCell[boundFace] * v[le] + (1.0 - faceWeightOfLeftCell[boundFace]) * v[re];
            RDouble wFace = faceWeightOfLeftCell[boundFace] * w[le] + (1.0 - faceWeightOfLeftCell[boundFace]) * w[re];
            RDouble ax = area[boundFace] * xfn[boundFace];
            RDouble ay = area[boundFace] * yfn[boundFace];
            RDouble az = area[boundFace] * zfn[boundFace];
            FaceFlux[boundFace] = rhoFace[boundFace] * (uFace * ax + vFace * ay + wFace * az);
            dun[boundFace] = 0;
        }

        InitialUnsteadyVar(grid);

        CalcWalldistOnWall(grid);
    }

    void MomEqCalculator::CalcWalldistOnWall(Grid*  gridIn)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);

        RDouble *xfn = grid->GetFaceNormalX();
        RDouble *yfn = grid->GetFaceNormalY();
        RDouble *zfn = grid->GetFaceNormalZ();
        RDouble *xfc = grid->GetFaceCenterX();
        RDouble *yfc = grid->GetFaceCenterY();
        RDouble *zfc = grid->GetFaceCenterZ();
        RDouble *xcc = grid->GetCellCenterX();
        RDouble *ycc = grid->GetCellCenterY();
        RDouble *zcc = grid->GetCellCenterZ();
        int* leftCellOfFace = grid->GetLeftCellOfFace();
        int* rightCellOfFace = grid->GetRightCellOfFace();
        RDouble *wd = reinterpret_cast<RDouble *> (grid->GetDataPtr("wd"));

        UnstructBCSet* unstructBCSet = grid->GetUnstructBCSet();
        int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
        for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
        {
            UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
            int bcType = bcRegion->GetBCType();
            if (bcType == 2)
            {
                vector<int>* faceIndex = bcRegion->GetFaceIndex();

                for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];

                    double FaceNormalX = -xfn[iFacelocal];
                    double FaceNormalY = -yfn[iFacelocal];
                    double FaceNormalZ = -zfn[iFacelocal];

                    double DisCelltoFaceX = xcc[le] - xfc[iFacelocal];
                    double DisCelltoFaceY = ycc[le] - yfc[iFacelocal];
                    double DisCelltoFaceZ = zcc[le] - zfc[iFacelocal];
                    wd[re] = (FaceNormalX * DisCelltoFaceX + FaceNormalY * DisCelltoFaceY + FaceNormalZ * DisCelltoFaceZ)
                        / std::sqrt(FaceNormalX * FaceNormalX + FaceNormalY * FaceNormalY + FaceNormalZ * FaceNormalZ);
                }
            }
        }
    }

    void MomEqCalculator::UpdateBoundaryPressureAndVelocity(Grid*  gridIn)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);
        RDouble *area    = grid->GetFaceArea();
        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

        for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
            int bcType = bcRegion->GetBCType();
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            Data_Param *bcData = bcRegion->GetBCParamDataBase();

            const int SubsonicInletBC       = 100001;
            const int SupersonicInletBC       = 100002;
            const int SupersonicOutletBC       = 100003;

            if (bcType == PHENGLEI::SOLID_SURFACE)
            {
                RDouble ub, vb, wb;
                int *leftCellOfFace = grid->GetLeftCellOfFace();
                int *rightCellOfFace = grid->GetRightCellOfFace();
                RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
                RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
                RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));

                for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];

                    int re = rightCellOfFace[iFacelocal];

                    if (bcData)
                    {
                        if (bcData->IsExist("initU", PHDOUBLE, 1))
                        {
                            bcData->GetData("initU", &ub, PHDOUBLE, 1);
                        }
                        else
                        {
                            ub = 0.0;
                        }

                        if (bcData->IsExist("initV", PHDOUBLE, 1))
                        {
                            bcData->GetData("initV", &vb, PHDOUBLE, 1);
                        }
                        else
                        {
                            vb = 0.0;
                        }

                        if (bcData->IsExist("initW", PHDOUBLE, 1))
                        {
                            bcData->GetData("initW", &wb, PHDOUBLE, 1);
                        }
                        else
                        {
                            wb = 0.0;
                        }

                        u[re] = RDouble(ub);
                        v[re] = RDouble(vb);
                        w[re] = RDouble(wb);
                    }
                }
            }
            else if (bcType == PHENGLEI::SYMMETRY)
            {
                int* leftCellOfFace = grid->GetLeftCellOfFace();
                int* rightCellOfFace = grid->GetRightCellOfFace();
                int nTotalCell = grid->GetNTotalCell();
                RDouble* u = reinterpret_cast<RDouble*>(grid->GetDataPtr("U"));
                RDouble* v = reinterpret_cast<RDouble*>(grid->GetDataPtr("V"));
                RDouble* w = reinterpret_cast<RDouble*>(grid->GetDataPtr("W"));
                RDouble* mu = reinterpret_cast<RDouble*>(grid->GetDataPtr("mu"));
                RDouble* xfn = grid->GetFaceNormalX();
                RDouble* yfn = grid->GetFaceNormalY();
                RDouble* zfn = grid->GetFaceNormalZ();

                for (std::size_t i = 0; i < faceIndex->size(); ++i)
                {
                    int iFace = (*faceIndex)[i];
                    int le = leftCellOfFace[iFace];
                    int re = rightCellOfFace[iFace];
                    RDouble muLocal = mu[le];

                    u[re] = u[le] - u[le] * xfn[iFace];
                    v[re] = v[le] - v[le] * yfn[iFace];
                    w[re] = w[le] - w[le] * zfn[iFace];
                }
            }
            else if (bcType == PHENGLEI::FARFIELD)
            {
                FlowFarField_GetBoundrayData(grid, bcData);

                RDouble* u = reinterpret_cast<RDouble*>(grid->GetDataPtr("U"));
                RDouble* v = reinterpret_cast<RDouble*>(grid->GetDataPtr("V"));
                RDouble* w = reinterpret_cast<RDouble*>(grid->GetDataPtr("W"));
                RDouble* FaceFlux = reinterpret_cast<RDouble*>(grid->GetDataPtr("FaceFlux"));        
                RDouble* rho = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
                RDouble* p = reinterpret_cast<RDouble*>(grid->GetDataPtr("P"));
                RDouble* T = reinterpret_cast<RDouble*>(grid->GetDataPtr("T"));
                RDouble refflowU;
                RDouble refflowV;
                RDouble refflowW;
                bcData->GetData("initU", &refflowU, PHDOUBLE, 1);
                bcData->GetData("initV", &refflowV, PHDOUBLE, 1);
                bcData->GetData("initW", &refflowW, PHDOUBLE, 1);

                RDouble refflowT;
                RDouble refflowP;
                bcData->GetData("initT", &refflowT, PHDOUBLE, 1);
                bcData->GetData("initP", &refflowP, PHDOUBLE, 1);

                int* leftCellOfFace = grid->GetLeftCellOfFace();
                int* rightCellOfFace = grid->GetRightCellOfFace();
                RDouble* xfn = grid->GetFaceNormalX();
                RDouble* yfn = grid->GetFaceNormalY();
                RDouble* zfn = grid->GetFaceNormalZ();
                for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];
                    if (FaceFlux[iFacelocal] < 0)
                    {
                        u[re] = refflowU;
                        v[re] = refflowV;
                        w[re] = refflowW;
                        p[re] = p[le];

                        FaceFlux[iFacelocal] = rho[re]* (u[re] * xfn[iFacelocal] + v[re] * yfn[iFacelocal] + w[re] * zfn[iFacelocal]);
                    }
                    else
                    {
                        u[re] = u[le];
                        v[re] = v[le];
                        w[re] = w[le];
                        p[re] = refflowP;
                        FaceFlux[iFacelocal] = rho[le] * (u[le] * xfn[iFacelocal] + v[le] * yfn[iFacelocal] + w[le] * zfn[iFacelocal]);
                    }
                }
            }
            else if (bcType == PHENGLEI::INFLOW)
            {
                RDouble ub, vb, wb;
                int *leftCellOfFace = grid->GetLeftCellOfFace();
                int *rightCellOfFace = grid->GetRightCellOfFace();
                RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
                RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
                RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));

                RDouble* xfn = grid->GetFaceNormalX();
                RDouble* yfn = grid->GetFaceNormalY();
                RDouble* zfn = grid->GetFaceNormalZ();

                for (size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];

                    int re = rightCellOfFace[iFacelocal];

                    if (bcData)
                    {
                        if (bcData->IsExist("initU", PHDOUBLE, 1))
                        {
                            bcData->GetData("initU", &ub, PHDOUBLE, 1);
                        }
                        else
                        {
                            ub = 0.0;
                        }

                        if (bcData->IsExist("initV", PHDOUBLE, 1))
                        {
                            bcData->GetData("initV", &vb, PHDOUBLE, 1);
                        }
                        else
                        {
                            vb = 0.0;
                        }

                        if (bcData->IsExist("initW", PHDOUBLE, 1))
                        {
                            bcData->GetData("initW", &wb, PHDOUBLE, 1);
                        }
                        else
                        {
                            wb = 0.0;
                        }

                        if (bcData->IsExist("totalVelocity", PHDOUBLE, 1))
                        {
                            RDouble totalVelocity;
                            bcData->GetData("totalVelocity", &totalVelocity, PHDOUBLE, 1);
                            ub = -totalVelocity * xfn[iFacelocal];
                            vb = -totalVelocity * yfn[iFacelocal];
                            wb = -totalVelocity * zfn[iFacelocal];
                        }

                        u[re] = ub;
                        v[re] = vb;
                        w[re] = wb;
                    }
                }
            }
            else if (bcType == PHENGLEI::OUTFLOW)
            {
            }
            else if (bcType == PHENGLEI::PRESSURE_INLET)
            {
                int* leftCellOfFace = grid->GetLeftCellOfFace();
                int* rightCellOfFace = grid->GetRightCellOfFace();
                RDouble* u = reinterpret_cast<RDouble*>(grid->GetDataPtr("U"));
                RDouble* v = reinterpret_cast<RDouble*>(grid->GetDataPtr("V"));
                RDouble* w = reinterpret_cast<RDouble*>(grid->GetDataPtr("W"));
                RDouble* area = grid->GetFaceArea();
                RDouble* xfn = grid->GetFaceNormalX();
                RDouble* yfn = grid->GetFaceNormalY();
                RDouble* zfn = grid->GetFaceNormalZ();
                RDouble* xfc = grid->GetFaceCenterX();
                RDouble* yfc = grid->GetFaceCenterY();
                RDouble* zfc = grid->GetFaceCenterZ();
                RDouble* xcc = grid->GetCellCenterX();
                RDouble* ycc = grid->GetCellCenterY();
                RDouble* zcc = grid->GetCellCenterZ();
                RDouble* rho = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
                RDouble* FaceFlux = reinterpret_cast<RDouble*>(grid->GetDataPtr("FaceFlux"));
                RDouble* p = reinterpret_cast<RDouble*>(grid->GetDataPtr("P"));
                RDouble* ptotal = reinterpret_cast<RDouble*>(grid->GetDataPtr("Ptotal"));
                for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];
                    RDouble ptotalb;
                    if (bcData)
                    {
                        if (bcData->IsExist("totalPressure", PHDOUBLE, 1))
                        {
                            bcData->GetData("totalPressure", &ptotalb, PHDOUBLE, 1);
                        }
                        else
                        {
                            ptotalb = 0.0;
                        }
            
                        ptotal[re] = ptotalb;
                        p[re] = 1;
                    }
                    RDouble vtotal = pow(((ptotal[re] - p[re]) * 2 / rho[re]), 0.5);
                    u[re] = -vtotal * xfn[iFacelocal];
                    v[re] = -vtotal * yfn[iFacelocal];
                    w[re] = -vtotal * zfn[iFacelocal];
                    FaceFlux[iFacelocal] = rho[re] * (u[re] * xfn[iFacelocal] * area[iFacelocal] + v[re] * yfn[iFacelocal] * area[iFacelocal] + w[re] * zfn[iFacelocal] * area[iFacelocal]);
                }
            }
            else if (bcType == PHENGLEI::PRESSURE_OUTLET)
            {
            }
            else 
            {
            }

            if (bcType == PHENGLEI::SOLID_SURFACE)

            {
                int *leftCellOfFace = grid->GetLeftCellOfFace();
                int *rightCellOfFace = grid->GetRightCellOfFace();

                RDouble *xfc = grid->GetFaceCenterX();
                RDouble *yfc = grid->GetFaceCenterY();
                RDouble *zfc = grid->GetFaceCenterZ();

                RDouble *xcc = grid->GetCellCenterX();
                RDouble *ycc = grid->GetCellCenterY();
                RDouble *zcc = grid->GetCellCenterZ();

                RDouble *p = reinterpret_cast<RDouble *>(grid->GetDataPtr("P"));
                RDouble *dpdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdx"));
                RDouble *dpdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdy"));
                RDouble *dpdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdz"));
                if (GlobalDataBase::GetIntParaFromDB("bodyForceFlag"))
                {
                    RDouble *bfU = reinterpret_cast<RDouble *>(grid->GetDataPtr("bfU"));
                    RDouble *bfV = reinterpret_cast<RDouble *>(grid->GetDataPtr("bfV"));
                    RDouble *bfW = reinterpret_cast<RDouble *>(grid->GetDataPtr("bfW"));
                    RDouble *cellVolume = grid->GetCellVolume();

                    for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                    {
                        int iFacelocal = (*faceIndex)[iFace];
                        int le = leftCellOfFace[iFacelocal];
                        int re = rightCellOfFace[iFacelocal];
                        RDouble faceCenterX = xfc[iFacelocal];
                        RDouble faceCenterY = yfc[iFacelocal];
                        RDouble faceCenterZ = zfc[iFacelocal];
                        RDouble cellCenterXle = xcc[le];
                        RDouble cellCenterYle = ycc[le];
                        RDouble cellCenterZle = zcc[le];

                        RDouble dx = faceCenterX - cellCenterXle;
                        RDouble dy = faceCenterY - cellCenterYle;
                        RDouble dz = faceCenterZ - cellCenterZle;
                        p[re] = p[le] + bfU[le] * dx + bfV[le] * dy + bfW[le] * dz;
                    }
                }
                else
                {
                    for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                    {
                        int iFacelocal = (*faceIndex)[iFace];
                        int le = leftCellOfFace[iFacelocal];
                        int re = rightCellOfFace[iFacelocal];
                        RDouble faceCenterX = xfc[iFacelocal];
                        RDouble faceCenterY = yfc[iFacelocal];
                        RDouble faceCenterZ = zfc[iFacelocal];
                        RDouble cellCenterXle = xcc[le];
                        RDouble cellCenterYle = ycc[le];
                        RDouble cellCenterZle = zcc[le];
                        RDouble dx = faceCenterX - cellCenterXle;
                        RDouble dy = faceCenterY - cellCenterYle;
                        RDouble dz = faceCenterZ - cellCenterZle;
                        p[re] = p[le];
                    }
                }
            }
            else if (bcType == PHENGLEI::SYMMETRY)
            {
                int* leftCellOfFace = grid->GetLeftCellOfFace();
                int* rightCellOfFace = grid->GetRightCellOfFace();
                RDouble* p = reinterpret_cast<RDouble*>(grid->GetDataPtr("P"));
                RDouble* dpdx = (RDouble*)grid->GetDataPtr("dpdx");
                RDouble* dpdy = (RDouble*)grid->GetDataPtr("dpdy");
                RDouble* dpdz = (RDouble*)grid->GetDataPtr("dpdz");
                RDouble* xfn = grid->GetFaceNormalX();
                RDouble* yfn = grid->GetFaceNormalY();
                RDouble* zfn = grid->GetFaceNormalZ();

                for (std::size_t i = 0; i < faceIndex->size(); ++i)
                {
                    int iFace = (*faceIndex)[i];
                    int re = rightCellOfFace[iFace];
                    int le = leftCellOfFace[iFace];

                    RDouble dpbdx = dpdx[le] - (dpdx[le] * xfn[iFace] + dpdy[le] * yfn[iFace] + dpdz[le] * zfn[iFace]) * xfn[iFace];
                    RDouble dpbdy = dpdy[le] - (dpdx[le] * xfn[iFace] + dpdy[le] * yfn[iFace] + dpdz[le] * zfn[iFace]) * yfn[iFace];
                    RDouble dpbdz = dpdz[le] - (dpdx[le] * xfn[iFace] + dpdy[le] * yfn[iFace] + dpdz[le] * zfn[iFace]) * zfn[iFace];
                    p[re] = p[le];
                }
            }
            else if (bcType == PHENGLEI::FARFIELD)
            {
                RDouble refMachNumber = 0.0;
                int* leftCellOfFace = grid->GetLeftCellOfFace();
                int* rightCellOfFace = grid->GetRightCellOfFace();
                RDouble* p = reinterpret_cast<RDouble*>(grid->GetDataPtr("P"));
                RDouble* FaceFlux = reinterpret_cast<RDouble*>(grid->GetDataPtr("FaceFlux"));
                RDouble* u = reinterpret_cast<RDouble*>(grid->GetDataPtr("U"));
                RDouble* v = reinterpret_cast<RDouble*>(grid->GetDataPtr("V"));
                RDouble* w = reinterpret_cast<RDouble*>(grid->GetDataPtr("W"));
                RDouble* xfn = grid->GetFaceNormalX();
                RDouble* yfn = grid->GetFaceNormalY();
                RDouble* zfn = grid->GetFaceNormalZ();

                RDouble* rho = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
                RDouble initRg = GlobalDataBase::GetDoubleParaFromDB("initRg");
                RDouble refflowU;
                RDouble refflowV;
                RDouble refflowW;
                bcData->GetData("initU", &refflowU, PHDOUBLE, 1);
                bcData->GetData("initV", &refflowV, PHDOUBLE, 1);
                bcData->GetData("initW", &refflowW, PHDOUBLE, 1);
                RDouble refflowT;
                RDouble refflowP;
                bcData->GetData("initT", &refflowT, PHDOUBLE, 1);
                bcData->GetData("initP", &refflowP, PHDOUBLE, 1);
                for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];
                    if (FaceFlux[iFacelocal] < 0)
                    {
                        u[re] = refflowU;
                        v[re] = refflowV;
                        w[re] = refflowW;
                        p[re] = p[le];
                        RDouble fluxbValue = rho[re] * area[iFacelocal] * (u[re] * xfn[iFacelocal] + v[re] * yfn[iFacelocal] + w[re] * zfn[iFacelocal]);
                        if (fluxbValue < 0)
                        {
                            FaceFlux[iFacelocal] = fluxbValue;
                        }
                        else
                        {
                            FaceFlux[iFacelocal] = 0;
                        }
                    }
                    else
                    {
                        u[re] = u[le];
                        v[re] = v[le];
                        w[re] = w[le];
                        p[re] = refflowP;
                        RDouble fluxbValue = rho[re] * area[iFacelocal] * (u[re] * xfn[iFacelocal] + v[re] * yfn[iFacelocal] + w[re] * zfn[iFacelocal]);

                        if (fluxbValue >= 0)
                        {
                            FaceFlux[iFacelocal] = fluxbValue;
                        }
                        else
                        {
                            FaceFlux[iFacelocal] = 0;
                        }
                    }
                }
            }
            else if (bcType == PHENGLEI::INFLOW)

            {
                int *leftCellOfFace = grid->GetLeftCellOfFace();
                int *rightCellOfFace = grid->GetRightCellOfFace();
                RDouble *p = reinterpret_cast<RDouble *>(grid->GetDataPtr("P"));
                RDouble *pc = reinterpret_cast<RDouble *>(grid->GetDataPtr("pc"));
                for (size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];
                    p[re] = p[le];
                }
            }
            else if (bcType == PHENGLEI::OUTFLOW)

            {
                RDouble pb;
                int *rightCellOfFace = grid->GetRightCellOfFace();
                RDouble *p = reinterpret_cast<RDouble *>(grid->GetDataPtr("P"));
                for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
                {
                    const int iFacelocal = (*faceIndex)[iFace];
                    int re = rightCellOfFace[iFacelocal];

                    if (bcData)
                    {
                        if (bcData->IsExist("initP", PHDOUBLE, 1))
                        {
                            bcData->GetData("initP", &pb, PHDOUBLE, 1);
                        }
                        else
                        {
                            pb = 0.0;
                        }

                        p[re] = pb;
                    }
                }
            }
            else if (bcType == PHENGLEI::PRESSURE_INLET)
            {
                RDouble ptotalb;
                int* leftCellOfFace = grid->GetLeftCellOfFace();
                int* rightCellOfFace = grid->GetRightCellOfFace();
            
                RDouble* p = reinterpret_cast<RDouble*>(grid->GetDataPtr("P"));
                RDouble* ptotal = reinterpret_cast<RDouble*>(grid->GetDataPtr("Ptotal"));
                RDouble* rho = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
                for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                {
                    int iFacelocal = (*faceIndex)[iFace];
                    int le = leftCellOfFace[iFacelocal];
                    int re = rightCellOfFace[iFacelocal];
                    RDouble* u = reinterpret_cast<RDouble*>(grid->GetDataPtr("U"));
                    RDouble* v = reinterpret_cast<RDouble*>(grid->GetDataPtr("V"));
                    RDouble* w = reinterpret_cast<RDouble*>(grid->GetDataPtr("W"));
                    if (bcData)
                    {
                        if (bcData->IsExist("totalPressure", PHDOUBLE, 1))
                        {
                            bcData->GetData("totalPressure", &ptotalb, PHDOUBLE, 1);
                        }
                        else
                        {
                            ptotalb = 0.0;
                        }
            
                        ptotal[re] = ptotalb;
                        p[re] = ptotal[re]-0.5*rho[le]*(pow(u[re],2)+pow(v[re], 2)+pow(w[re], 2));
                    }
                }
            }
            else if (bcType == PHENGLEI::PRESSURE_OUTLET)
            
            {
                RDouble pb;
                int *rightCellOfFace = grid->GetRightCellOfFace();
                RDouble *p = reinterpret_cast<RDouble *>(grid->GetDataPtr("P"));
                for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
                {
                    const int iFacelocal = (*faceIndex)[iFace];
                    int re = rightCellOfFace[iFacelocal];
            
                    if (bcData)
                    {
                        if (bcData->IsExist("initP", PHDOUBLE, 1))
                        {
                            bcData->GetData("initP", &pb, PHDOUBLE, 1);
                        }
                        else
                        {
                            pb = 0.0;
                        }
            
                        p[re] = pb;
                    }
                }
            }
            else if (bcType == SubsonicInletBC)
            {
            }
            else if (bcType == SupersonicInletBC)
            {
            }
            else if (bcType == SupersonicOutletBC)
            
            {
                RDouble pb;
                int* leftCellOfFace = grid->GetLeftCellOfFace();
                int* rightCellOfFace = grid->GetRightCellOfFace();
                RDouble* p = reinterpret_cast<RDouble*>(grid->GetDataPtr("P"));
            
                int isSupersonic = GlobalDataBase::GetIntParaFromDB("isSupersonic");
                if (isSupersonic)
                {
                    for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                    {
                        int iFacelocal = (*faceIndex)[iFace];
                        int le = leftCellOfFace[iFacelocal];
                        int re = rightCellOfFace[iFacelocal];
                        p[re] = p[le];
                    }
                }
                else
                {
                    for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
                    {
                        int iFacelocal = (*faceIndex)[iFace];
                        int le = leftCellOfFace[iFacelocal];
                        int re = rightCellOfFace[iFacelocal];
            
                        if (bcData)
                        {
                            if (bcData->IsExist("initP", PHDOUBLE, 1))
                            {
                                bcData->GetData("initP", &pb, PHDOUBLE, 1);
                            }
                            else
                            {
                                pb = 0.0;
                            }
            
                            p[re] = pb;
                        }
                    }
                }
            }
        }
    }

    void MomEqCalculator::AllocateGlobalVar(Grid*  gridIn)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);
        int nTotalFace = grid->GetNTotalFace();
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nIFace = grid->GetNIFace();
        int nTotal = nTotalCell + nBoundFace;

        RDouble *u = NewPointer<RDouble>(nTotal);
        grid->UpdateDataPtr("U", u);

        RDouble *v = NewPointer<RDouble>(nTotal);
        grid->UpdateDataPtr("V", v);

        RDouble *w = NewPointer<RDouble>(nTotal);
        grid->UpdateDataPtr("W", w);

        RDouble *p = NewPointer<RDouble>(nTotal);
        grid->UpdateDataPtr("P", p);

        RDouble *pc = NewPointer<RDouble>(nTotal);
        grid->UpdateDataPtr("pc", pc);

        RDouble *ptotal = NewPointer<RDouble>(nTotal);
        grid->UpdateDataPtr("Ptotal", ptotal);

        RDouble *diagMatrixCoeff = NewPointer<RDouble>(nTotal);
        grid->UpdateDataPtr("diagMatrixCoeff", diagMatrixCoeff);

        RDouble *bCoeff = NewPointer<RDouble>(nTotal);
        grid->UpdateDataPtr("bCoeff", bCoeff);

        RDouble *upperMatrixCoeff = NewPointer<RDouble>(nTotalFace);
        grid->UpdateDataPtr("upperMatrixCoeff", upperMatrixCoeff);

        RDouble *lowerMatrixCoeff = NewPointer<RDouble>(nTotalFace);
        grid->UpdateDataPtr("lowerMatrixCoeff", lowerMatrixCoeff);

        RDouble *FaceFlux = NewPointer<RDouble>(nTotalFace);
        grid->UpdateDataPtr("FaceFlux", FaceFlux);

        int bodyForceFlag = GlobalDataBase::GetIntParaFromDB("bodyForceFlag");
        if (bodyForceFlag == 1)
        {
            RDouble *bfU = NewPointer<RDouble>(nTotal);
            grid->UpdateDataPtr("bfU", bfU);

            RDouble *bfV = NewPointer<RDouble>(nTotal);
            grid->UpdateDataPtr("bfV", bfV);

            RDouble *bfW = NewPointer<RDouble>(nTotal);
            grid->UpdateDataPtr("bfW", bfW);

            RDouble* dbfdx = NewPointer<RDouble>(nTotalCell);
            RDouble* dbfdy = NewPointer<RDouble>(nTotalCell);
            RDouble* dbfdz = NewPointer<RDouble>(nTotalCell);

            grid->UpdateDataPtr("dbfdx", dbfdx);
            grid->UpdateDataPtr("dbfdy", dbfdy);
            grid->UpdateDataPtr("dbfdz", dbfdz);
        }

        RDouble *dun = NewPointer<RDouble>(nTotalFace);
        grid->UpdateDataPtr("dun", dun);

        RDouble *rho = NewPointer<RDouble>(nTotal);
        grid->UpdateDataPtr("rho", rho);

        RDouble *rhoFace = NewPointer<RDouble>(nTotalFace);
        grid->UpdateDataPtr("rhoFace", rhoFace);

        RDouble *mu = NewPointer<RDouble>(nTotal);
        grid->UpdateDataPtr("mu", mu);

        RDouble *hu = NewPointer<RDouble>(nTotal);
        grid->UpdateDataPtr("hu", hu);

        RDouble *InverseAp = NewPointer<RDouble>(nTotal);
        grid->UpdateDataPtr("InverseAp", InverseAp);

        RDouble *dudx = NewPointer<RDouble>(nTotal);
        RDouble *dudy = NewPointer<RDouble>(nTotal);
        RDouble *dudz = NewPointer<RDouble>(nTotal);
        RDouble *dvdx = NewPointer<RDouble>(nTotal);
        RDouble *dvdy = NewPointer<RDouble>(nTotal);
        RDouble *dvdz = NewPointer<RDouble>(nTotal);
        RDouble *dwdx = NewPointer<RDouble>(nTotal);
        RDouble *dwdy = NewPointer<RDouble>(nTotal);
        RDouble *dwdz = NewPointer<RDouble>(nTotal);
        RDouble *dpdx = NewPointer<RDouble>(nTotal);
        RDouble *dpdy = NewPointer<RDouble>(nTotal);
        RDouble *dpdz = NewPointer<RDouble>(nTotal);
        RDouble *dpcdx = NewPointer<RDouble>(nTotal);
        RDouble *dpcdy = NewPointer<RDouble>(nTotal);
        RDouble *dpcdz = NewPointer<RDouble>(nTotal);

        grid->UpdateDataPtr("dUdx", dudx);
        grid->UpdateDataPtr("dUdy", dudy);
        grid->UpdateDataPtr("dUdz", dudz);
        grid->UpdateDataPtr("dVdx", dvdx);
        grid->UpdateDataPtr("dVdy", dvdy);
        grid->UpdateDataPtr("dVdz", dvdz);
        grid->UpdateDataPtr("dWdx", dwdx);
        grid->UpdateDataPtr("dWdy", dwdy);
        grid->UpdateDataPtr("dWdz", dwdz);
        grid->UpdateDataPtr("dpdx", dpdx);
        grid->UpdateDataPtr("dpdy", dpdy);
        grid->UpdateDataPtr("dpdz", dpdz);
        grid->UpdateDataPtr("dpcdx", dpcdx);
        grid->UpdateDataPtr("dpcdy", dpcdy);
        grid->UpdateDataPtr("dpcdz", dpcdz);

        RDouble* faceWeightOfLeftCell = NewPointer<RDouble>(nTotalFace);
        grid->UpdateDataPtr("faceWeightOfLeftCell", faceWeightOfLeftCell);
        
        RDouble* wd = NewPointer<RDouble>(nTotal);
        grid->UpdateDataPtr("wd", wd);

        AllocateUnsteadyVar(grid);

        string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
        if (flowSolverName == "CompressibleSIMPLE")
        {
            RDouble* cRho = NewPointer<RDouble>(nTotal);
            grid->UpdateDataPtr("cRho", cRho);
            PHSPACE::SetField(cRho, 0.0, nTotal);

            RDouble initRg = GlobalDataBase::GetDoubleParaFromDB("initRg");
            RDouble* Rg = NewPointer<RDouble>(nTotal);
            grid->UpdateDataPtr("Rg", Rg);
            PHSPACE::SetField(Rg, initRg, nTotal);
        }
    }

    void MomEqCalculator::InitFlowAsRestart(Grid*  grid)
    {
        GAS_SPACE::gas->SetGasName("speciesNameIncom");
        GAS_SPACE::gas->SetGasMolarMass("gasMolarMass");
        GAS_SPACE::gas->SetGasRho("gasRho");
        GAS_SPACE::gas->SetGasMu("gasMu");

        int iunsteady =  GlobalDataBase::GetIntParaFromDB("iunsteady");
        if (iunsteady)
        {
            int timeStepNow = 0;
            GlobalDataBase::UpdateData("timeStepNow", &timeStepNow, PHINT, 1);
            RDouble startTime = GlobalDataBase::GetDoubleParaFromDB("startTime");
            GlobalDataBase::UpdateData("currentTime", &startTime, PHDOUBLE, 1);
        }

        int nTotalFace = grid->GetNTotalFace();
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal = nTotalCell + nBoundFace;

        RDouble* u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
        RDouble initU = GlobalDataBase::GetDoubleParaFromDB("initU");
        PHSPACE::SetField(u, initU, nTotal);

        RDouble* v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
        RDouble initV = GlobalDataBase::GetDoubleParaFromDB("initV");
        PHSPACE::SetField(v, initV, nTotal);

        RDouble* w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
        RDouble initW = GlobalDataBase::GetDoubleParaFromDB("initW");
        PHSPACE::SetField(w, initW, nTotal);

        RDouble *p = reinterpret_cast<RDouble *> (grid->GetDataPtr("P"));
        RDouble initP = GlobalDataBase::GetDoubleParaFromDB("initP");
        PHSPACE::SetField(p, initP, nTotal);

        RDouble* ptotal = reinterpret_cast<RDouble*>(grid->GetDataPtr("Ptotal"));
        PHSPACE::SetField(ptotal, initP, nTotal);

        RDouble *pc = reinterpret_cast<RDouble *> (grid->GetDataPtr("pc"));
        PHSPACE::SetField(pc, 0.0, nTotal);

        RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("diagMatrixCoeff"));
        PHSPACE::SetField(diagMatrixCoeff, 0.0, nTotal);

        RDouble *bCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("bCoeff"));
        PHSPACE::SetField(bCoeff, 0.0, nTotal);

        RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("upperMatrixCoeff"));
        PHSPACE::SetField(upperMatrixCoeff, 0.0, nTotalFace);

        RDouble *lowerMatrixCoeff = reinterpret_cast<RDouble *> (grid->GetDataPtr("lowerMatrixCoeff"));
        PHSPACE::SetField(lowerMatrixCoeff, 0.0, nTotalFace);

        RDouble* FaceFlux = reinterpret_cast<RDouble*>(grid->GetDataPtr("FaceFlux"));
        PHSPACE::SetField(FaceFlux, 0.0, nTotalFace);

        RDouble* dun  = reinterpret_cast<RDouble*>(grid->GetDataPtr("dun"));
        PHSPACE::SetField(dun, 0.0, nTotalFace);

        RDouble* rho  = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
        RDouble initRho = GlobalDataBase::GetDoubleParaFromDB("initRho");
        PHSPACE::SetField(rho, initRho, nTotal);

        RDouble *rhoFace = reinterpret_cast<RDouble*>(grid->GetDataPtr("rhoFace"));
        PHSPACE::SetField(rhoFace, 0.0, nTotalFace);

        RDouble* mu   = reinterpret_cast<RDouble*>(grid->GetDataPtr("mu"));
        RDouble initMu = GlobalDataBase::GetDoubleParaFromDB("initMu");
        PHSPACE::SetField(mu, initMu, nTotal);

        RDouble **DiffusionCoeff = reinterpret_cast <RDouble **> (GlobalDataBase::GetDataPtr("DiffusionCoeff"));
        DiffusionCoeff[IDX::S_IU] = mu;
        DiffusionCoeff[IDX::S_IV] = mu;
        DiffusionCoeff[IDX::S_IW] = mu;

        calcFaceWeight(grid);

        RDouble *hu = reinterpret_cast<RDouble *> (grid->GetDataPtr("hu"));
        PHSPACE::SetField(hu, 0.0, nTotal);

        RDouble *InverseAp = reinterpret_cast<RDouble *> (grid->GetDataPtr("InverseAp"));
        PHSPACE::SetField(InverseAp, 0.0, nTotal);

        RDouble *wd = reinterpret_cast<RDouble *> (grid->GetDataPtr("wd"));
        PHSPACE::SetField(wd, 0.0, nTotal);


        RDouble *dudx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdx"));
        RDouble *dudy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdy"));
        RDouble *dudz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdz"));
        RDouble *dvdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdx"));
        RDouble *dvdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdy"));
        RDouble *dvdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdz"));
        RDouble *dwdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdx"));
        RDouble *dwdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdy"));
        RDouble *dwdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdz"));
        RDouble* dpdx = reinterpret_cast<RDouble*>(grid->GetDataPtr("dpdx"));
        RDouble* dpdy = reinterpret_cast<RDouble*>(grid->GetDataPtr("dpdy"));
        RDouble* dpdz = reinterpret_cast<RDouble*>(grid->GetDataPtr("dpdz"));
        RDouble *dpcdx = reinterpret_cast<RDouble *> (grid->GetDataPtr("dpcdx"));
        RDouble *dpcdy = reinterpret_cast<RDouble *> (grid->GetDataPtr("dpcdy"));
        RDouble *dpcdz = reinterpret_cast<RDouble *> (grid->GetDataPtr("dpcdz"));
        PHSPACE::SetField(dudx, 0.0, nTotal);
        PHSPACE::SetField(dudy, 0.0, nTotal);
        PHSPACE::SetField(dudz, 0.0, nTotal);
        PHSPACE::SetField(dvdx, 0.0, nTotal);
        PHSPACE::SetField(dvdy, 0.0, nTotal);
        PHSPACE::SetField(dvdz, 0.0, nTotal);
        PHSPACE::SetField(dwdx, 0.0, nTotal);
        PHSPACE::SetField(dwdy, 0.0, nTotal);
        PHSPACE::SetField(dwdz, 0.0, nTotal);
        PHSPACE::SetField(dpdx, 0.0, nTotal);
        PHSPACE::SetField(dpdy, 0.0, nTotal);
        PHSPACE::SetField(dpdz, 0.0, nTotal);
        PHSPACE::SetField(dpcdx, 0.0, nTotal);
        PHSPACE::SetField(dpcdy, 0.0, nTotal);
        PHSPACE::SetField(dpcdz, 0.0, nTotal);

        string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
        if (flowSolverName == "CompressibleSIMPLE")
        {
            RDouble* cRho = reinterpret_cast<RDouble*>(grid->GetDataPtr("cRho"));
            PHSPACE::SetField(cRho, 0.0, nTotal);

            RDouble* Rg = reinterpret_cast<RDouble*>(grid->GetDataPtr("Rg"));
            RDouble initRg = GlobalDataBase::GetDoubleParaFromDB("initRg");
            PHSPACE::SetField(Rg, initRg, nTotal);
        }

        int bodyForceFlag = GlobalDataBase::GetIntParaFromDB("bodyForceFlag");
        if (bodyForceFlag == 1)
        {
            RDouble *bfU = reinterpret_cast<RDouble *> (grid->GetDataPtr("bfU"));
            RDouble *bfV = reinterpret_cast<RDouble *> (grid->GetDataPtr("bfV"));
            RDouble *bfW = reinterpret_cast<RDouble *> (grid->GetDataPtr("bfW"));
            PHSPACE::SetField(bfU, 0.0, nTotal);
            PHSPACE::SetField(bfV, 0.0, nTotal);
            PHSPACE::SetField(bfW, 0.0, nTotal);

            RDouble *dbfdx = reinterpret_cast<RDouble *> (grid->GetDataPtr("dbfdx"));
            RDouble *dbfdy = reinterpret_cast<RDouble *> (grid->GetDataPtr("dbfdy"));
            RDouble *dbfdz = reinterpret_cast<RDouble *> (grid->GetDataPtr("dbfdz"));
            PHSPACE::SetField(dbfdx, 0.0, nTotalCell);
            PHSPACE::SetField(dbfdy, 0.0, nTotalCell);
            PHSPACE::SetField(dbfdz, 0.0, nTotalCell);
        }
    }

    void MomEqCalculator::AllocateUnsteadyVar(Grid*  grid)
    {
        PHString1D phiNameList;
        phiNameList.resize(4);

        phiNameList[0] = "U";
        phiNameList[1] = "V";
        phiNameList[2] = "W";
        phiNameList[3] = "P";
        phiNameList.push_back("rho");

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

    void MomEqCalculator::InitialUnsteadyVar(Grid*  grid)
    {
        PHString1D phiNameList;
        phiNameList.resize(4);

        phiNameList[0] = "U";
        phiNameList[1] = "V";
        phiNameList[2] = "W";
        phiNameList[3] = "P";


        string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
        if (flowSolverName == "CompressibleSIMPLE")
        {
            phiNameList.push_back("rho");
        }

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

    void MomEqCalculator::UpdateUnsteadyProperties(Grid* gridIn)
    {
        UnstructGrid * grid = UnstructGridCast(gridIn);
        UpdateProperties(grid);
    }

    void MomEqCalculator::UpdateProperties(Grid* gridIn)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);
        string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
        if (flowSolverName == "CompressibleSIMPLE")
        {
            int isStableUnsteadyMethod = GlobalDataBase::GetIntParaFromDB("isStableUnsteadyMethod");
            if (isStableUnsteadyMethod)
            {
                return;
            }
        }

        RDouble* rho  = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
        RDouble* mu   = reinterpret_cast<RDouble*>(grid->GetDataPtr("mu"));
        GAS_SPACE::gas->UpdateRho(grid, rho);
        GAS_SPACE::gas->UpdateMu(grid, mu);

        CommunicateAnInterfaceVar(rho);
        CommunicateAnInterfaceVar(mu);

        RDouble* cRho = reinterpret_cast<RDouble*>(grid->GetDataPtr("cRho"));

        if (cRho != nullptr)
        {
            CommunicateAnInterfaceVar(cRho);
        }
    }

    void MomEqCalculator::GetResidual(Grid *gridIn, vector<RDouble>& res)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);

        RDouble resNow = 0.0;

        grid->GetData("UResNow", &resNow, PHDOUBLE, 1);
        res.push_back(resNow);
        grid->GetData("VResNow", &resNow, PHDOUBLE, 1);
        res.push_back(resNow);
        grid->GetData("WResNow", &resNow, PHDOUBLE, 1);
        res.push_back(resNow);
        grid->GetData("PResNow", &resNow, PHDOUBLE, 1);
        res.push_back(resNow);
    }

    void MomEqCalculator::UpdateUnsteadyVariable(Grid *grid)
    {
        UnstructGrid *gridIn = UnstructGridCast(grid);

        vector<string> phiNameList;
        phiNameList.push_back("U");
        phiNameList.push_back("V");
        phiNameList.push_back("W");
        phiNameList.push_back("P");

        string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
        if (flowSolverName == "CompressibleSIMPLE")
        {
            phiNameList.push_back("rho");
        }

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
        vector<string>().swap(phiNameList);
    }

    void FlowFarField_GetBoundrayData(Grid*  grid, Data_Param* bcData)
    {
        RDouble refMachNumber = 0.0;
        RDouble nx;
        RDouble ny;
        RDouble nz;
        RDouble refflowT;
        RDouble refflowP;
        if (bcData)
        {
            if (bcData->IsExist("refMachNumber", PHDOUBLE, 1))
            {
                bcData->GetData("refMachNumber", &refMachNumber, PHDOUBLE, 1);
                bcData->GetData("nx", &nx, PHDOUBLE, 1);
                bcData->GetData("ny", &ny, PHDOUBLE, 1);
                bcData->GetData("nz", &nz, PHDOUBLE, 1);
                bcData->GetData("initT", &refflowT, PHDOUBLE, 1);
                bcData->GetData("initP", &refflowP, PHDOUBLE, 1);
            }
        }
        RDouble gama = 1.4;
        RDouble initRg = GlobalDataBase::GetDoubleParaFromDB("initRg");
        RDouble totalVelocity = sqrt(gama * initRg * refflowT) * refMachNumber;
        RDouble refflowU = totalVelocity * nx;
        RDouble refflowV = totalVelocity * ny;
        RDouble refflowW = totalVelocity * nz;

        bcData->UpdateData("initU", &refflowU, PHDOUBLE, 1);
        bcData->UpdateData("initV", &refflowV, PHDOUBLE, 1);
        bcData->UpdateData("initW", &refflowW, PHDOUBLE, 1);
        bcData->UpdateData("initT", &refflowT, PHDOUBLE, 1);
        bcData->UpdateData("initP", &refflowP, PHDOUBLE, 1);
    }


}