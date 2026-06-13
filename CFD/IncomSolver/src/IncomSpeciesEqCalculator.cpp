#include "IncomSpeciesEqCalculator.h"




namespace PHSPACE
{
IncomSpeciesEqCalculator::IncomSpeciesEqCalculator():IncomScalarEqCalculator()
{

}

IncomSpeciesEqCalculator::~IncomSpeciesEqCalculator()
{

}

void IncomSpeciesEqCalculator::InitFlowAsRestart(Grid *gridIn)
{
    IncomScalarEqCalculator::InitFlowAsRestart(gridIn);
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    int *localSolverIndex = reinterpret_cast <int *> (GlobalDataBase::GetDataPtr("localSolverIndex"));
    int localIndex = localSolverIndex[solverIndex];

    int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
    RDouble *initMassFractionIncom = new RDouble[numberOfSpeciesIncom];
    GlobalDataBase::GetData("initMassFractionIncom", initMassFractionIncom, PHDOUBLE, numberOfSpeciesIncom);
    RDouble initPhi = initMassFractionIncom[localIndex];
    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(varNameIncom[solverIndex]));
    PHSPACE::SetField(phi, initPhi, nTotal);

    GAS_SPACE::gas->InitParameterForMassdiff();

    delete [] initMassFractionIncom;
}

void IncomSpeciesEqCalculator::IncompressibleInitial(Grid * gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int *localSolverIndex = reinterpret_cast <int *> (GlobalDataBase::GetDataPtr("localSolverIndex"));
    RDouble **DiffusionCoeff = reinterpret_cast <RDouble **> (GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    int solverIndex = GetSolverIndex();

    GAS_SPACE::gas->UpdateMassDiff(grid, localSolverIndex[solverIndex], DiffusionCoeff[solverIndex]);

    UpdateBCValue(grid);

    InitialUnsteadyVar(grid);
}


void IncomSpeciesEqCalculator::GetResidual(Grid *gridIn, vector<RDouble>& res)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    string varName = varNameIncom[solverIndex];
    RDouble resNow = 0.0;
    grid->GetData(varName + "ResNow", &resNow, PHDOUBLE, 1);

    res.push_back(resNow);
}

void IncomSpeciesEqCalculator::InitialUnsteadyVar(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
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
        ImplicitEuler_ReInitTimeVar(grid, phiNameList);
    }
    else if (TranCalcMethod == "IMPLICIT_2ND_ORDER")
    {
        Implicit2ndOrder_ReInitTimeVar(grid, phiNameList);
    }

    PHString1D().swap(phiNameList);
}

void IncomSpeciesEqCalculator::UpdateUnsteadyVariable(Grid *grid)
{
    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    string varName = varNameIncom[solverIndex];
    UnstructGrid *gridIn = UnstructGridCast(grid); 
    std::vector<std::string> phiNameList;
    string *scalarNameList = GAS_SPACE::gas->GetGasName();
    phiNameList.push_back(varName);

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

void IncomSpeciesEqCalculator::UpdateBCValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    string varName = varNameIncom[solverIndex];
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
    
    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int>* faceIndex = bcRegion->GetFaceIndex();
        Data_Param* bcData = bcRegion->GetBCParamDataBase();

        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            for (std::size_t i = 0; i < faceIndex->size(); ++i)
            {
                const int iface = (*faceIndex)[i];
                int le = leftCellOfFace[iface];
                int re = rightCellOfFace[iface];
                RDouble  spBound = phi[le];

                int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
                if (bcData->IsExist("initMassFractionIncom", PHDOUBLE, numberOfSpeciesIncom))
                {
                    RDouble *initMassFractionIncom = new RDouble[numberOfSpeciesIncom];
                    bcData->GetData("initMassFractionIncom", initMassFractionIncom, PHDOUBLE, numberOfSpeciesIncom);

                    int *localSolverIndex = reinterpret_cast <int *> (GlobalDataBase::GetDataPtr("localSolverIndex"));
                    int localIndex = localSolverIndex[solverIndex];
                    spBound = initMassFractionIncom[localIndex];

                    delete [] initMassFractionIncom;
                }
                phi[re] = RDouble(spBound);
            }

        }
        else if (bcType == PHENGLEI::INFLOW)
        {
            for (std::size_t i = 0; i < faceIndex->size(); ++i)
            {
                const int iface = (*faceIndex)[i];
                int re = rightCellOfFace[iface];
                RDouble spBound;
                int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
                if (bcData->IsExist("initMassFractionIncom", PHDOUBLE, numberOfSpeciesIncom))
                {
                    RDouble *initMassFractionIncom = new RDouble[numberOfSpeciesIncom];
                    bcData->GetData("initMassFractionIncom", initMassFractionIncom, PHDOUBLE, numberOfSpeciesIncom);

                    int *localSolverIndex = reinterpret_cast <int *> (GlobalDataBase::GetDataPtr("localSolverIndex"));
                    int localIndex = localSolverIndex[solverIndex];
                    spBound = initMassFractionIncom[localIndex];

                    delete [] initMassFractionIncom;
                }
                else
                {
                    std::cout << " No initial value has assigned to the boundary " << std::endl;
                }
                if (re == 3908)
                {
                    int mm = 0;
                }
                phi[re] = RDouble(spBound);
            }

        }
        else if (bcType == PHENGLEI::OUTFLOW)
        {
            for (std::size_t i = 0; i < faceIndex->size(); ++i)
            {
                const int iface = (*faceIndex)[i];
                int le = leftCellOfFace[iface];
                int re = rightCellOfFace[iface];
                RDouble spBound = phi[le];

                if (FaceFlux[iface] < 0.0)
                {
                    int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
                    if (bcData->IsExist("initMassFractionIncom", PHDOUBLE, numberOfSpeciesIncom))
                    {
                        RDouble *initMassFractionIncom = new RDouble[numberOfSpeciesIncom];
                        bcData->GetData("initMassFractionIncom", initMassFractionIncom, PHDOUBLE, numberOfSpeciesIncom);

                        int *localSolverIndex = reinterpret_cast <int *> (GlobalDataBase::GetDataPtr("localSolverIndex"));
                        int localIndex = localSolverIndex[solverIndex];
                        spBound = initMassFractionIncom[localIndex];

                        delete [] initMassFractionIncom;
                    }
                }

                phi[re] = spBound;
            }
        }
        else if (bcType == PHENGLEI::PRESSURE_INLET)
        {
            for (std::size_t i = 0; i < faceIndex->size(); ++i)
            {
                RDouble spBound;

                const int iface = (*faceIndex)[i];
                int re = rightCellOfFace[iface];
                int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
                if (bcData->IsExist("initMassFractionIncom", PHDOUBLE, numberOfSpeciesIncom))
                {
                    RDouble *initMassFractionIncom = new RDouble[numberOfSpeciesIncom];
                    bcData->GetData("initMassFractionIncom", initMassFractionIncom, PHDOUBLE, numberOfSpeciesIncom);

                    int *localSolverIndex = reinterpret_cast <int *> (GlobalDataBase::GetDataPtr("localSolverIndex"));
                    int localIndex = localSolverIndex[solverIndex];
                    spBound = initMassFractionIncom[localIndex];

                    delete [] initMassFractionIncom;
                }
                else
                {
                    std::cout << " No initial value has assigned to the boundary " << std::endl;
                }

                phi[re] = spBound;
            }

        }
        else if (bcType == PHENGLEI::PRESSURE_OUTLET)
        {
            for (std::size_t i = 0; i < faceIndex->size(); ++i)
            {
                const int iface = (*faceIndex)[i];
                int le = leftCellOfFace[iface];
                int re = rightCellOfFace[iface];

                RDouble spBound = phi[le];

                if (FaceFlux[iface] < 0.0)
                {
                    int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
                    if (bcData->IsExist("initMassFractionIncom", PHDOUBLE, numberOfSpeciesIncom))
                    {
                        RDouble *initMassFractionIncom = new RDouble[numberOfSpeciesIncom];
                        bcData->GetData("initMassFractionIncom", initMassFractionIncom, PHDOUBLE, numberOfSpeciesIncom);

                        int *localSolverIndex = reinterpret_cast <int *> (GlobalDataBase::GetDataPtr("localSolverIndex"));
                        int localIndex = localSolverIndex[solverIndex];
                        spBound = initMassFractionIncom[localIndex];

                        delete [] initMassFractionIncom;
                    }
                }
                phi[re] = spBound;
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

void IncomSpeciesEqCalculator::UpdateAfterIterloop(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
    int *localSolverIndex = reinterpret_cast <int *> (GlobalDataBase::GetDataPtr("localSolverIndex"));
    int solverIndex = GetSolverIndex();
    if (localSolverIndex[solverIndex] != numberOfSpeciesIncom - 1)
    {
        return;
    }

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int nTotal = nTotalCell + nBoundFace;

    int *leftCellOfFace = grid->GetLeftCellOfFace();

    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *T = reinterpret_cast<RDouble *>(grid->GetDataPtr("T"));
    RDouble *p = reinterpret_cast<RDouble *>(grid->GetDataPtr("P"));

    RDouble *Rg = reinterpret_cast<RDouble *>(grid->GetDataPtr("Rg"));
    RDouble *cp = reinterpret_cast<RDouble *>(grid->GetDataPtr("cp"));
    RDouble *cRho = reinterpret_cast<RDouble *>(grid->GetDataPtr("cRho"));

    RDouble *spSum = new RDouble[nTotalCell];
    string *scalarNameList = GAS_SPACE::gas->GetGasName();
    RDouble *sp = (RDouble *)grid->GetDataPtr(scalarNameList[0]);

    SetField(spSum, 0.0, nTotalCell);
    
    int numOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
    for (int iSp = 0; iSp < numOfSpecies; ++iSp)
    {
        sp = (RDouble *)grid->GetDataPtr(scalarNameList[iSp]);
        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            if (sp[iCell] < 0.0)
            {
                sp[iCell] = 0.0;
            }
    
            spSum[iCell] += sp[iCell];
        }
    }
    
    for (int iSp = 0; iSp < numOfSpecies; iSp++)
    {
        sp = (RDouble *)grid->GetDataPtr(scalarNameList[iSp]);
        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            sp[iCell] /= spSum[iCell];
        }
    }

    UpdateSpBoundary(grid);

    delete [] spSum;

    for (int iSp = 0; iSp < numOfSpecies; iSp++)
    {
        sp = (RDouble *)grid->GetDataPtr(scalarNameList[iSp]);
        CommunicateAnInterfaceVar(sp);
    }
}

void IncomSpeciesEqCalculator::UpdateSpBoundary(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
    string *scalarNameList = GAS_SPACE::gas->GetGasName();
    int numOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
    for (int iSp = 0; iSp < numOfSpecies; ++iSp)
    {
        RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(scalarNameList[iSp]));
        for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
            int bcType = bcRegion->GetBCType();
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            Data_Param *bcData = bcRegion->GetBCParamDataBase();

            if (bcType == PHENGLEI::SOLID_SURFACE)
            {
                for (std::size_t i = 0; i < faceIndex->size(); ++i)
                {
                    const int iface = (*faceIndex)[i];
                    int le = leftCellOfFace[iface];
                    int re = rightCellOfFace[iface];
                    RDouble spBound = phi[le];

                    int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
                    if (bcData->IsExist("initMassFractionIncom", PHDOUBLE, numberOfSpeciesIncom))
                    {
                        RDouble *initMassFractionIncom = new RDouble[numberOfSpeciesIncom];
                        bcData->GetData("initMassFractionIncom", initMassFractionIncom, PHDOUBLE, numberOfSpeciesIncom);
                        spBound = initMassFractionIncom[iSp];
                        delete [] initMassFractionIncom;
                    }

                    phi[re] = spBound;
                }
            }
            else if (bcType == PHENGLEI::INFLOW)
            {
                for (std::size_t i = 0; i < faceIndex->size(); ++i)
                {
                    RDouble spBound;

                    const int iface = (*faceIndex)[i];
                    int re = rightCellOfFace[iface];
                    int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
                    if (bcData->IsExist("initMassFractionIncom", PHDOUBLE, numberOfSpeciesIncom))
                    {
                        RDouble *initMassFractionIncom = new RDouble[numberOfSpeciesIncom];
                        bcData->GetData("initMassFractionIncom", initMassFractionIncom, PHDOUBLE, numberOfSpeciesIncom);
                        spBound = initMassFractionIncom[iSp];
                        delete [] initMassFractionIncom;
                    }
                    else
                    {
                        std::cout << " No initial value has assigned to the boundary " << std::endl;
                    }

                    phi[re] = spBound;
                }

            }
            else if (bcType == PHENGLEI::OUTFLOW)
            {
                for (std::size_t i = 0; i < faceIndex->size(); ++i)
                {
                    const int iface = (*faceIndex)[i];
                    int le = leftCellOfFace[iface];
                    int re = rightCellOfFace[iface];

                    RDouble spBound = phi[le];
                    if (FaceFlux[iface] < 0.0)
                    {
                        int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
                        if (bcData->IsExist("initMassFractionIncom", PHDOUBLE, numberOfSpeciesIncom))
                        {
                            RDouble *initMassFractionIncom = new RDouble[numberOfSpeciesIncom];
                            bcData->GetData("initMassFractionIncom", initMassFractionIncom, PHDOUBLE, numberOfSpeciesIncom);
                            spBound = initMassFractionIncom[iSp];
                            delete [] initMassFractionIncom;
                        }
                    }
                    phi[re] = spBound;
                }
            }
            else if (bcType == PHENGLEI::PRESSURE_INLET)
            {
                for (std::size_t i = 0; i < faceIndex->size(); ++i)
                {
                    RDouble spBound;

                    const int iface = (*faceIndex)[i];
                    int re = rightCellOfFace[iface];

                    int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
                    if (bcData->IsExist("initMassFractionIncom", PHDOUBLE, numberOfSpeciesIncom))
                    {
                        RDouble *initMassFractionIncom = new RDouble[numberOfSpeciesIncom];
                        bcData->GetData("initMassFractionIncom", initMassFractionIncom, PHDOUBLE, numberOfSpeciesIncom);
                        spBound = initMassFractionIncom[iSp];
                        delete [] initMassFractionIncom;
                    }
                    else
                    {
                        std::cout << " No initial value has assigned to the boundary " << std::endl;
                    }
                    phi[re] = spBound;
                }
            }
            else if (bcType == PHENGLEI::PRESSURE_OUTLET)
            {
                for (std::size_t i = 0; i < faceIndex->size(); ++i)
                {
                    const int iface = (*faceIndex)[i];
                    int le = leftCellOfFace[iface];
                    int re = rightCellOfFace[iface];
                    RDouble spBound = phi[le];
                    
                    if (FaceFlux[iface] < 0.0)
                    {
                        int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
                        if (bcData->IsExist("initMassFractionIncom", PHDOUBLE, numberOfSpeciesIncom))
                        {
                            RDouble *initMassFractionIncom = new RDouble[numberOfSpeciesIncom];
                            bcData->GetData("initMassFractionIncom", initMassFractionIncom, PHDOUBLE, numberOfSpeciesIncom);
                            spBound = initMassFractionIncom[iSp];
                            delete [] initMassFractionIncom;
                        }
                    }
                    phi[re] = spBound;
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
    }
}

void IncomSpeciesEqCalculator::UpdateProperties(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
    int *localSolverIndex = reinterpret_cast <int *> (GlobalDataBase::GetDataPtr("localSolverIndex"));
    int solverIndex = GetSolverIndex();
    if (localSolverIndex[solverIndex] != numberOfSpeciesIncom - 1)
    {
        return;
    }

    int isStableUnsteadyMethod = GlobalDataBase::GetIntParaFromDB("isStableUnsteadyMethod");
    if (isStableUnsteadyMethod)
    {
        return;
    }

    RDouble *Rg = reinterpret_cast<RDouble *>(grid->GetDataPtr("Rg"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    GAS_SPACE::gas->UpdateRg(grid, Rg);
    GAS_SPACE::gas->UpdateRho(grid, rho);

    CommunicateAnInterfaceVar(Rg);
    CommunicateAnInterfaceVar(rho);
}

void IncomSpeciesEqCalculator::SetDiffusionCoeff(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int *localSolverIndex = reinterpret_cast<int *>(GlobalDataBase::GetDataPtr("localSolverIndex"));
    RDouble **DiffusionCoeff = reinterpret_cast<RDouble **>(GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    int solverIndex = GetSolverIndex();
    GAS_SPACE::gas->UpdateMassDiff(grid, localSolverIndex[solverIndex], DiffusionCoeff[solverIndex]);
    CommunicateAnInterfaceVar(DiffusionCoeff[solverIndex]);
}

}