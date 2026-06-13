#pragma once
#include"IncomUnstructSolver.h"
#include"TK_Time.h"
#include"unapCalculator.h"
#include"yhamgCalculator.h"
#include"Post_ForceMoment.h"
#include"GradientOperation.h"
#include "PHIO.h"


namespace PHSPACE
{
    IncomUnstructSolver::IncomUnstructSolver()
    {
        FreeControlParameters();
        controlParameters = new Param_INCompSolverUnstruct();

        
        int mathLibType = GlobalDataBase::GetIntParaFromDB("mathLibType");
        if (mathLibType == HYPRELib)
        {
            AxEqualbCalculator = new hypreCalculator();
        }
        else if (mathLibType == UNAPLib)
        {
            AxEqualbCalculator = new unapCalculator();
        }
        else if (mathLibType == YHAMGLib)
        {
            AxEqualbCalculator = new yhamgCalculator();
        }

        momEqCalculator = new MomEqCalculator();

        PPEquaCalculator = new PPEqCalculator();

        RDouble *soluRelaxCoeff = reinterpret_cast <RDouble *> (GlobalDataBase::GetDataPtr("soluRelaxCoeff"));
        RDouble MomEqRelaxCoeff = GlobalDataBase::GetDoubleParaFromDB("MomEqRelaxCoeff");
        RDouble PPEqRelaxCoeff = GlobalDataBase::GetDoubleParaFromDB("PPEqRelaxCoeff");

        soluRelaxCoeff[IDX::S_IU] = MomEqRelaxCoeff;
        soluRelaxCoeff[IDX::S_IV] = MomEqRelaxCoeff;
        soluRelaxCoeff[IDX::S_IW] = MomEqRelaxCoeff;
        soluRelaxCoeff[IDX::S_IP] = PPEqRelaxCoeff;

        int *localSolverIndex = reinterpret_cast <int *> (GlobalDataBase::GetDataPtr("localSolverIndex"));
        localSolverIndex[IDX::S_IU] = 0;
        localSolverIndex[IDX::S_IV] = 1;
        localSolverIndex[IDX::S_IW] = 2;
        localSolverIndex[IDX::S_IP] = 3;

        int *solMethodOfLinEqSystem = reinterpret_cast <int *> (GlobalDataBase::GetDataPtr("solMethodOfLinEqSystem"));
        int MomEqSolMethod = GlobalDataBase::GetIntParaFromDB("MomEqSolMethod");
        int PPEqSolMethod = GlobalDataBase::GetIntParaFromDB("PPEqSolMethod");
        solMethodOfLinEqSystem[IDX::S_IU] = MomEqSolMethod;
        solMethodOfLinEqSystem[IDX::S_IV] = MomEqSolMethod;
        solMethodOfLinEqSystem[IDX::S_IW] = MomEqSolMethod;
        solMethodOfLinEqSystem[IDX::S_IP] = PPEqSolMethod;

        int *precondMethodOfLinEqSystem = reinterpret_cast <int *> (GlobalDataBase::GetDataPtr("precondMethodOfLinEqSystem"));
        int MomEqPrecondMethod = GlobalDataBase::GetIntParaFromDB("MomEqPrecondMethod");
        int PPEqPrecondMethod = GlobalDataBase::GetIntParaFromDB("PPEqPrecondMethod");
        precondMethodOfLinEqSystem[IDX::S_IU] = MomEqPrecondMethod;
        precondMethodOfLinEqSystem[IDX::S_IV] = MomEqPrecondMethod;
        precondMethodOfLinEqSystem[IDX::S_IW] = MomEqPrecondMethod;
        precondMethodOfLinEqSystem[IDX::S_IP] = PPEqPrecondMethod;
        
        int *maxIterOfEachEq = reinterpret_cast <int *> (GlobalDataBase::GetDataPtr("maxIterOfEachEq"));
        int MomEqMaxSweep = GlobalDataBase::GetIntParaFromDB("MomEqMaxSweep");
        int PPEqMaxSweep = GlobalDataBase::GetIntParaFromDB("PPEqMaxSweep");
        maxIterOfEachEq[IDX::S_IU] = MomEqMaxSweep;
        maxIterOfEachEq[IDX::S_IV] = MomEqMaxSweep;
        maxIterOfEachEq[IDX::S_IW] = MomEqMaxSweep;
        maxIterOfEachEq[IDX::S_IP] = PPEqMaxSweep;

        RDouble *iterTolOfEachEq = reinterpret_cast <RDouble *> (GlobalDataBase::GetDataPtr("iterTolOfEachEq"));
        RDouble MomEqIterSolvTol = GlobalDataBase::GetDoubleParaFromDB("MomEqIterSolvTol");
        RDouble PPEqIterSolvTol = GlobalDataBase::GetDoubleParaFromDB("PPEqIterSolvTol");
        iterTolOfEachEq[IDX::S_IU] = MomEqIterSolvTol;
        iterTolOfEachEq[IDX::S_IV] = MomEqIterSolvTol;
        iterTolOfEachEq[IDX::S_IW] = MomEqIterSolvTol;
        iterTolOfEachEq[IDX::S_IP] = PPEqIterSolvTol;

        int isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
        if (isSolveTurb == 1)
        {
            string TurbEqSolverName = GlobalDataBase::GetStrParaFromDB("TurbEqSolverName");
            if (TurbEqSolverName == "KE")
            {
                RDouble TurbEqRelaxCoeff = GlobalDataBase::GetDoubleParaFromDB("TurbEqRelaxCoeff");

                IncomTurbKSolver = new IncomKETurbKEqCalculator();
                IncomTurbKSolver->SetSolverIndex(IDX::S_ITURBK);
                soluRelaxCoeff[IDX::S_ITURBK] = TurbEqRelaxCoeff;
                localSolverIndex[IDX::S_ITURBK] = 0;

                int TurbEqSolMethod = GlobalDataBase::GetIntParaFromDB("TurbEqSolMethod");
                solMethodOfLinEqSystem[IDX::S_ITURBK] = TurbEqSolMethod;
                int TurbEqPrecondMethod = GlobalDataBase::GetIntParaFromDB("TurbEqPrecondMethod");
                precondMethodOfLinEqSystem[IDX::S_ITURBK] = TurbEqPrecondMethod;
                int TurbEqMaxSweep = GlobalDataBase::GetIntParaFromDB("TurbEqMaxSweep");
                maxIterOfEachEq[IDX::S_ITURBK] = TurbEqMaxSweep;
                RDouble TurbEqIterSolvTol = GlobalDataBase::GetDoubleParaFromDB("TurbEqIterSolvTol");
                iterTolOfEachEq[IDX::S_ITURBK] = TurbEqIterSolvTol;

                IncomTurbEplisionSolver = new IncomKETurbEpsilonEqCalculator();
                IncomTurbEplisionSolver->SetSolverIndex(IDX::S_IEPSILON);
                soluRelaxCoeff[IDX::S_IEPSILON] = TurbEqRelaxCoeff;
                localSolverIndex[IDX::S_IEPSILON] = 1;
                solMethodOfLinEqSystem[IDX::S_IEPSILON] = TurbEqSolMethod;
                precondMethodOfLinEqSystem[IDX::S_IEPSILON] = TurbEqPrecondMethod;
                maxIterOfEachEq[IDX::S_IEPSILON] = TurbEqMaxSweep;
                iterTolOfEachEq[IDX::S_IEPSILON] = TurbEqIterSolvTol;
            }
            else if (TurbEqSolverName == "SA")
            {
                IncomTurbSASolver = new IncomSATurbKEqCalculator();
                IncomTurbSASolver->SetSolverIndex(IDX::S_ISA);

                RDouble TurbEqRelaxCoeff = GlobalDataBase::GetDoubleParaFromDB("TurbEqRelaxCoeff");
                soluRelaxCoeff[IDX::S_ISA] = TurbEqRelaxCoeff;

                localSolverIndex[IDX::S_ISA] = 0;

                int TurbEqSolMethod = GlobalDataBase::GetIntParaFromDB("TurbEqSolMethod");
                solMethodOfLinEqSystem[IDX::S_ISA] = TurbEqSolMethod;

                int TurbEqPrecondMethod = GlobalDataBase::GetIntParaFromDB("TurbEqPrecondMethod");
                precondMethodOfLinEqSystem[IDX::S_ISA] = TurbEqPrecondMethod;

                int TurbEqMaxSweep = GlobalDataBase::GetIntParaFromDB("TurbEqMaxSweep");
                maxIterOfEachEq[IDX::S_ISA] = TurbEqMaxSweep;

                RDouble TurbEqIterSolvTol = GlobalDataBase::GetDoubleParaFromDB("TurbEqIterSolvTol");
                iterTolOfEachEq[IDX::S_ISA] = TurbEqIterSolvTol;
            }
        }

        int isSolveEnergy = GlobalDataBase::GetIntParaFromDB("isSolveEnergy");
        if (isSolveEnergy == 1)
        {
            IncomEnergySolver = new IncomEnergyEqCalculator();
            IncomEnergySolver->SetSolverIndex(IDX::S_ITEMP);
            RDouble EnergyEqRelaxCoeff = GlobalDataBase::GetDoubleParaFromDB("EnergyEqRelaxCoeff");
            soluRelaxCoeff[IDX::S_ITEMP] = EnergyEqRelaxCoeff;

            localSolverIndex[IDX::S_ITEMP] = 0;

            int EnergyEqSolMethod = GlobalDataBase::GetIntParaFromDB("EnergyEqSolMethod");
            solMethodOfLinEqSystem[IDX::S_ITEMP] = EnergyEqSolMethod;

            int EnergyEqPrecondMethod = GlobalDataBase::GetIntParaFromDB("EnergyEqPrecondMethod");
            precondMethodOfLinEqSystem[IDX::S_ITEMP] = EnergyEqPrecondMethod;

            int EnergyEqMaxSweep = GlobalDataBase::GetIntParaFromDB("EnergyEqMaxSweep");
            maxIterOfEachEq[IDX::S_ITEMP] = EnergyEqMaxSweep;

            RDouble EnergyEqIterSolvTol = GlobalDataBase::GetDoubleParaFromDB("EnergyEqIterSolvTol");
            iterTolOfEachEq[IDX::S_ITEMP] = EnergyEqIterSolvTol;
        }

        int isSolveSpecies = GlobalDataBase::GetIntParaFromDB("isSolveSpecies");
        if (isSolveSpecies == 1)
        {

            RDouble SpeciesEqRelaxCoeff = GlobalDataBase::GetDoubleParaFromDB("SpeciesEqRelaxCoeff");
            int SpeciesEqSolMethod = GlobalDataBase::GetIntParaFromDB("SpeciesEqSolMethod");
            int SpeciesEqPrecondMethod = GlobalDataBase::GetIntParaFromDB("SpeciesEqPrecondMethod");
            int SpeciesEqMaxSweep = GlobalDataBase::GetIntParaFromDB("SpeciesEqMaxSweep");
            RDouble SpeciesEqIterSolvTol = GlobalDataBase::GetDoubleParaFromDB("SpeciesEqIterSolvTol");

            int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
            IncomSpeciesSolver = new IncomScalarEqCalculator * [numberOfSpeciesIncom];
            string *gasNameList = new string[numberOfSpeciesIncom];
            GlobalDataBase::GetData("speciesNameIncom", gasNameList, PHSTRING, numberOfSpeciesIncom);
            for (int iGas = 0; iGas < numberOfSpeciesIncom; iGas++)
            {
                IncomSpeciesSolver[iGas] = new IncomSpeciesEqCalculator();

                string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
                for (int iIndex = 0; iIndex < 50; iIndex++)
                {
                    if (varNameIncom[iIndex] == gasNameList[iGas])
                    {
                        IncomSpeciesSolver[iGas]->SetSolverIndex(iIndex);
                        soluRelaxCoeff[iIndex] = SpeciesEqRelaxCoeff;
                        localSolverIndex[iIndex] = iGas;
                        
                        solMethodOfLinEqSystem[iIndex] = SpeciesEqSolMethod;
                        precondMethodOfLinEqSystem[iIndex] = SpeciesEqPrecondMethod;
                        maxIterOfEachEq[iIndex] = SpeciesEqMaxSweep;
                        iterTolOfEachEq[iIndex] = SpeciesEqIterSolvTol;
                    }
                }
            }
            delete [] gasNameList;
        }

    }

    IncomUnstructSolver::~IncomUnstructSolver()
    {

    }

    void IncomUnstructSolver::DumpResultFile(Grid *grid, int level)
    {
        int outnstep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
        int intervalStepRes = GlobalDataBase::GetIntParaFromDB("intervalStepRes");
        int intervalStepFlow = GlobalDataBase::GetIntParaFromDB("intervalStepFlow");
        int intervalStepPlot = GlobalDataBase::GetIntParaFromDB("intervalStepPlot");
        int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");

        if (outnstep % intervalStepRes == 0)
        {
            DumpResidual(grid);
        }

        bool isSubIterationDump = true;
        if (isUnsteady)
        {
            isSubIterationDump = GlobalDataBase::GetIntParaFromDB("isSubIterationDump");
        }

        if (outnstep % intervalStepFlow == 0 && level == 0 && isSubIterationDump)
        {
            //! To dump restart data for continual simulation.
            ActionKey *actkeyDumpRestartData = new ActionKey();
            FillActionKey(actkeyDumpRestartData, DUMP_RESTART, 0);
            DumpRestartData(actkeyDumpRestartData);
            delete actkeyDumpRestartData;
        }

        int intervalStepForce = GlobalDataBase::GetIntParaFromDB("intervalStepForce");
        if (intervalStepForce > 0)
        {
            if (outnstep % intervalStepForce == 0 && isSubIterationDump || outnstep == 1)
            {
                if (level == 0)
                {
                    //! To dump airforce coefficient data.
                    ActionKey *actkeyDumpAirForceCoef = new ActionKey();
                    FillActionKey(actkeyDumpAirForceCoef, DUMP_AIR_FORCE_COEF, level);
                    DumpAirForceCoef(actkeyDumpAirForceCoef);
                    delete actkeyDumpAirForceCoef;

                    //! To dump wall airforce coefficient data.
                    ActionKey *actkeyDumpWallAircoef = new ActionKey();
                    FillActionKey(actkeyDumpWallAircoef, DUMP_CP_DISTRI, 0);
                    //DumpCpDistri(actkeyDumpWallAircoef);
                    delete actkeyDumpWallAircoef;
                }
            }
        }
    }

    void IncomUnstructSolver::AirForceCoef(ActionKey *actkey)
    {
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

    void IncomUnstructSolver::AirForceCoefParts(ActionKey *actkey, int partID, int Coordinate)
    {
        int level = actkey->level;

        using namespace PHMPI;
        //fstream file;

        UnstructGrid * grid = UnstructGridCast(GetGrid(level));

        const int ALL_PART = -1;
        SimpleBC * boundaryCondition = NULL;
        string partBCName = "";
        RDouble TorqueRefX = 0.0, TorqueRefY = 0.0, TorqueRefZ = 0.0;

        //! Used to compute hinge moment:
        int dumpHingeMoment = 0;
        RDouble localCoordAxis0[3];
        RDouble localCoordAxis1[3];
        RDouble hingeMoment = 0.0;

        if (partID == ALL_PART)
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

        using namespace IDX;

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

        RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
        RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
        RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
        RDouble *p = reinterpret_cast<RDouble *>(grid->GetDataPtr("P"));

        RDouble *dUdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdx"));
        RDouble *dUdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdy"));
        RDouble *dUdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdz"));
        RDouble *dVdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdx"));
        RDouble *dVdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdy"));
        RDouble *dVdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdz"));
        RDouble *dWdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdx"));
        RDouble *dWdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdy"));
        RDouble *dWdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdz"));
        RDouble poo = GlobalDataBase::GetDoubleParaFromDB("initP");
        RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));

        RDouble initRho = GlobalDataBase::GetDoubleParaFromDB("initRho");
        RDouble initU = GlobalDataBase::GetDoubleParaFromDB("initU");
        RDouble initV = GlobalDataBase::GetDoubleParaFromDB("initV");
        RDouble initW = GlobalDataBase::GetDoubleParaFromDB("initW");
        RDouble refVelocity = sqrt(initU * initU + initV * initV + initW * initW);
        if (refVelocity < 1.0e-10) refVelocity = 1.0e-10;

        RDouble refDynamicPressure = initRho * refVelocity * refVelocity;

        RDouble cpx, cpy, cpz;
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

        RDouble pw,cp;
        RDouble dfx,dfy,dfz,dpx,dpy,dpz,fvsx,fvsy,fvsz;
        RDouble nx,ny,nz;
        RDouble vis,txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz;
        RDouble dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz;
        RDouble xc,yc,zc,dx,dy,dz,ods;

        using namespace IDX;

        // get the bcRegion information of unstruct grid.
        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int nBCRegion = unstructBCSet->GetnBCRegion();

        // cycle all bcregions and find the solid surface type of boundary.
        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
            int bcType = bcRegion->GetBCType();
            string bcName = bcRegion->GetBCName();

            // if the bc type is wall, calculate the wall coefficient.
            if (bcType == PHENGLEI::SOLID_SURFACE)
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

                    //pressure drag
                    pw = p[le];
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

                    {
                        dudx = dUdx[le];
                        dudy = dUdy[le];
                        dudz = dUdz[le];

                        dvdx = dVdx[le];
                        dvdy = dVdy[le];
                        dvdz = dVdz[le];

                        dwdx = dWdx[le];
                        dwdy = dWdy[le];
                        dwdz = dWdz[le];

                        //! Gradient correction.
                        dx  = xcc[re] - xcc[le];
                        dy  = ycc[re] - ycc[le];
                        dz  = zcc[re] - zcc[le];
                        ods = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
                        dx  = dx * ods;
                        dy  = dy * ods;
                        dz  = dz * ods;

                        CorrectGradient(u[le], u[re],dudx,dudy,dudz,dx,dy,dz,ods);
                        CorrectGradient(v[le], v[re],dvdx,dvdy,dvdz,dx,dy,dz,ods);
                        CorrectGradient(w[le], w[re],dwdx,dwdy,dwdz,dx,dy,dz,ods);

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

                        fvsx = - two * (nx * txx + ny * tyx + nz * tzx);
                        fvsy = - two * (nx * txy + ny * tyy + nz * tzy);
                        fvsz = - two * (nx * txz + ny * tyz + nz * tzz);

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

                    if (dumpHingeMoment)
                    {
                        Post_ForceMoment forceMoment;
                        RDouble point[3] = {xfc[iFace], yfc[iFace], zfc[iFace]};
                        RDouble faceForce[3] = {dfx, dfy, dfz};
                        hingeMoment += forceMoment.ComputeMoment(point, localCoordAxis0, localCoordAxis1, faceForce);
                    }
                }
            }
        }

        CA_f   /= refDynamicPressure;
        CA_p   /= refDynamicPressure;
        CN_f   /= refDynamicPressure;
        CN_p   /= refDynamicPressure;
        CZ_f   /= refDynamicPressure;
        CZ_p   /= refDynamicPressure;
        cpx    /= refDynamicPressure;
        cpy    /= refDynamicPressure;
        cpz    /= refDynamicPressure;
        Cl_f   /= refDynamicPressure;
        Cl_p   /= refDynamicPressure;
        Cn_f   /= refDynamicPressure;
        Cn_p   /= refDynamicPressure;
        Cm_f   /= refDynamicPressure;
        Cm_p   /= refDynamicPressure;

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

    void IncomUnstructSolver::InitFlowAsReadingRestart(ActionKey *actkey)
    {
        InitFlowAsRestart();
        actkey->filename = GetRestartFileName();
        if (actkey->filename == "")
        {
            return;
        }

        if (PHMPI::IsParallelRun())
        {
            actkey->filename = PHSPACE::AddSymbolToFileName(actkey->filename, "_", 0);
        }

        hid_t file;
        file = OpenHDF5File(actkey->filename);
        actkey->filepos = file;

        ReadRestartH5(actkey);

        H5Fclose(file);
        actkey->filepos = 0;
    }

    void IncomUnstructSolver::ReadRestartH5(ActionKey *actkey)
    {
        Grid *grid = GetGrid(actkey->level);
        int gridID = grid->GetZoneID();

        int outnstep = 0;
        ReadData(actkey->filepos, &outnstep, "outnstep");
        GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);

        RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
        RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
        RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
        RDouble *p = reinterpret_cast<RDouble *>(grid->GetDataPtr("P"));
        RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));

        int nBoundFace = grid->GetNBoundFace();
        int nTotalCell = grid->GetNTotalCell();
        int nTotal = nTotalCell + nBoundFace;

        hid_t grploc;
        string grpName;

        ostringstream oss;
        oss << "Group" << gridID;
        grpName = oss.str();

        grploc = OpenGroup(actkey->filepos, grpName);

        int nTotalCellRestart = 0;
        ReadData(grploc, &nTotalCellRestart, "nTotalCell");

        if (nTotalCellRestart != nTotalCell)
        {
            ostringstream erroeInfo;
            erroeInfo << " Error: the cell number in flow.dat is not equal to the cell number in grid file !" << endl;
            TK_Exit::ExceptionExit(erroeInfo.str());
        }

        ReadData(grploc, u, "U");
        ReadData(grploc, v, "V");
        ReadData(grploc, w, "W");
        ReadData(grploc, p, "P");
        ReadData(grploc, rho, "rho");
        string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
        if (flowSolverName == "CompressibleSIMPLE")
        {
            RDouble *cRho = reinterpret_cast<RDouble *>(grid->GetDataPtr("cRho"));
            ReadData(grploc, cRho, "cRho");
        }

        int isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
        if (isSolveTurb == 1)
        {
            RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
            RDouble *e = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));

            string TurbEqSolverName = GlobalDataBase::GetStrParaFromDB("TurbEqSolverName");
            if (TurbEqSolverName == "KE")
            {
                ReadData(grploc, k, "Kinetic");
                ReadData(grploc, e, "Epsilon");
            }
            else if (TurbEqSolverName == "SA")
            {
                ReadData(grploc, k, "Kinetic");
            }
        }

        int isSolveEnergy = GlobalDataBase::GetIntParaFromDB("isSolveEnergy");
        if (isSolveEnergy == 1)
        {
            RDouble *T = reinterpret_cast<RDouble *> (grid->GetDataPtr("T"));
            ReadData(grploc, T, "T");
        }

        int isSolveSpecies = GlobalDataBase::GetIntParaFromDB("isSolveSpecies");
        if (isSolveSpecies == 1)
        {
            string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
            int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
            for (int iGas = 0; iGas < numberOfSpeciesIncom; iGas++)
            {
                IncomSpeciesEqCalculator* SpeciesSolver = static_cast<IncomSpeciesEqCalculator *>(IncomSpeciesSolver[iGas]);
                int solverIndex = SpeciesSolver->GetSolverIndex();
                string varName = varNameIncom[solverIndex];
                RDouble *q = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
                ReadData(grploc, q, varName);
            }
        }

        H5Gclose(grploc);
    }

    void IncomUnstructSolver::DumpRestartH5(ActionKey *actkey)
    {
        using namespace PHMPI;
        int currentProcessorID = GetCurrentProcessorID();
        int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
        int version = 1;
        if (currentProcessorID == GetServerProcessorID())
        {
            WriteData(actkey->filepos, &version, "Version");
            WriteData(actkey->filepos, &outnstep, "outnstep");
        }

        Grid *grid = GetGrid(actkey->level);
        int gridID = grid->GetZoneID();

        int nEquation = GetNumberOfEquations();

        int nBoundFace = grid->GetNBoundFace();
        int nTotalCell = grid->GetNTotalCell();
        int nTotal = nTotalCell + nBoundFace;

        RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
        RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
        RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
        RDouble *p = reinterpret_cast<RDouble *>(grid->GetDataPtr("P"));
        RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));

        hid_t grploc;
        string grpName;

        ostringstream oss;
        oss << "Group" << gridID;
        grpName = oss.str();

        grploc = OpenGroup(actkey->filepos, grpName);

        WriteData(grploc, u, "U");
        WriteData(grploc, v, "V");
        WriteData(grploc, w, "W");
        WriteData(grploc, p, "P");
        WriteData(grploc, rho, "rho");

        string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
        if (flowSolverName == "CompressibleSIMPLE")
        {
            RDouble *cRho = reinterpret_cast<RDouble *>(grid->GetDataPtr("cRho"));
            WriteData(grploc, cRho, "cRho");
        }

        int isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
        if (isSolveTurb == 1)
        {
            RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
            RDouble *e = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));

            string TurbEqSolverName = GlobalDataBase::GetStrParaFromDB("TurbEqSolverName");
            if (TurbEqSolverName == "KE")
            {
                WriteData(grploc, k, "Kinetic");
                WriteData(grploc, e, "Epsilon");
            }
            else if (TurbEqSolverName == "SA")
            {
                WriteData(grploc, k, "Kinetic");
            }
        }

        int isSolveEnergy = GlobalDataBase::GetIntParaFromDB("isSolveEnergy");
        if (isSolveEnergy == 1)
        {
            RDouble *T = reinterpret_cast<RDouble *>(grid->GetDataPtr("T"));
            WriteData(grploc, T, "T");
        }

        int isSolveSpecies = GlobalDataBase::GetIntParaFromDB("isSolveSpecies");
        if (isSolveSpecies == 1)
        {
            string *varNameIncom = reinterpret_cast <string *>(GlobalDataBase::GetDataPtr("varNameIncom"));
            int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
            for (int iGas = 0; iGas < numberOfSpeciesIncom; iGas++)
            {
                IncomSpeciesEqCalculator *SpeciesSolver = static_cast<IncomSpeciesEqCalculator *>(IncomSpeciesSolver[iGas]);
                int solverIndex = SpeciesSolver->GetSolverIndex();
                string varName = varNameIncom[solverIndex];
                RDouble *q = reinterpret_cast<RDouble *> (grid->GetDataPtr(varName));
                WriteData(grploc, q, varName);
            }
        }

        WriteData(grploc, &nTotalCell, "nTotalCell");

        H5Gclose(grploc);
    }

    void IncomUnstructSolver::CreateH5RestartFile(ActionKey *actkey)
    {
        using namespace PHMPI;
        int nZones = GetNumberofGlobalZones();
        int currentProcessorID = GetCurrentProcessorID();
        int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            ostringstream oss;
            oss << "Group" << iZone;

            string grpName;

            grpName = oss.str();
            CreateGroup(actkey->filepos, grpName);

            int length = 0;
            int faceLength = 0;
            int nWallBC = 0, nCurMaxCellNum = 0;
            int sendProcessID = GetZoneProcessorID(iZone);
            if (currentProcessorID == sendProcessID)
            {
                int gridTypeofLocalZone = GetZoneGridType(iZone);
                Grid *grid = PHSPACE::GetGrid(iZone, actkey->level);
                int nBoundFace = grid->GetNBoundFace();
                int nTotalCell = grid->GetNTotalCell();
                faceLength = grid->GetNTotalFace();
                length = nBoundFace + nTotalCell;
            }

            int getherLength;
            MPI_Allreduce(&length, &getherLength, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            int getherFaceLength;
            MPI_Allreduce(&faceLength, &getherFaceLength, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            hid_t  grpData;
            grpData = OpenGroup(actkey->filepos, grpName);

            CreateEmptyData(grpData, "U", 1, getherLength, PHDOUBLE);
            CreateEmptyData(grpData, "V", 1, getherLength, PHDOUBLE);
            CreateEmptyData(grpData, "W", 1, getherLength, PHDOUBLE);
            CreateEmptyData(grpData, "P", 1, getherLength, PHDOUBLE);
            CreateEmptyData(grpData, "rho", 1, getherLength, PHDOUBLE);

            string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
            if (flowSolverName == "CompressibleSIMPLE")
            {
                CreateEmptyData(grpData, "cRho", 1, getherLength, PHDOUBLE);
            }

            int isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
            if (isSolveTurb == 1)
            {
                string TurbEqSolverName = GlobalDataBase::GetStrParaFromDB("TurbEqSolverName");
                if (TurbEqSolverName == "KE")
                {
                    CreateEmptyData(grpData, "Kinetic", 1, getherLength, PHDOUBLE);
                    CreateEmptyData(grpData, "Epsilon", 1, getherLength, PHDOUBLE);
                }
                else if (TurbEqSolverName == "SA")
                {
                    CreateEmptyData(grpData, "Kinetic", 1, getherLength, PHDOUBLE);
                }
            }

            int isSolveEnergy = GlobalDataBase::GetIntParaFromDB("isSolveEnergy");
            if (isSolveEnergy == 1)
            {
                CreateEmptyData(grpData, "T", 1, getherLength, PHDOUBLE);
            }

            int isSolveSpecies = GlobalDataBase::GetIntParaFromDB("isSolveSpecies");
            if (isSolveSpecies == 1)
            {
                string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
                int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
                for (int iGas = 0; iGas < numberOfSpeciesIncom; iGas++)
                {
                    IncomSpeciesEqCalculator* SpeciesSolver = static_cast<IncomSpeciesEqCalculator *>(IncomSpeciesSolver[iGas]);
                    int solverIndex = SpeciesSolver->GetSolverIndex();
                    string varName = varNameIncom[solverIndex];
                    CreateEmptyData(grpData, varName, 1, getherLength, PHDOUBLE);
                }
            }

            CreateEmptyData(grpData, "nTotalCell", 1, 1, PHINT);

            H5Gclose(grpData);
        }

        CreateEmptyData(actkey->filepos, "Version", 1, 1, PHINT);
        CreateEmptyData(actkey->filepos, "outnstep", 1, 1, PHINT);
    }

    void IncomUnstructSolver::CpDistriCoef(ActionKey *actkey)
    {
        int level = actkey->level;
        DataContainer *cdata = actkey->GetData();
        cdata->MoveToBegin();

        UnstructGrid *grid = UnstructGridCast(GetGrid(level));
        int GridID = grid->GetGridID()->GetIndex();
        PHWrite(cdata, GridID);

        using namespace PHMPI;
        using namespace PHENGLEI;
        using namespace IDX;

        RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
        RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
        RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
        RDouble *p = reinterpret_cast<RDouble *>(grid->GetDataPtr("P"));
        RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));

        RDouble *dUdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdx"));
        RDouble *dUdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdy"));
        RDouble *dUdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdz"));
        RDouble *dVdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdx"));
        RDouble *dVdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdy"));
        RDouble *dVdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdz"));
        RDouble *dWdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdx"));
        RDouble *dWdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdy"));
        RDouble *dWdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdz"));

        RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("initRho");
        RDouble initU = GlobalDataBase::GetDoubleParaFromDB("initU");
        RDouble initV = GlobalDataBase::GetDoubleParaFromDB("initV");
        RDouble initW = GlobalDataBase::GetDoubleParaFromDB("initW");
        RDouble refVelocity = sqrt(initU * initU + initV * initV + initW * initW );
        RDouble attack= atan(initV/initU);
        RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
        RDouble poo = GlobalDataBase::GetDoubleParaFromDB("initP");

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


        RDouble sina = 0, cosa = 0;
        sina = sin(attack);
        cosa = cos(attack);


        RDouble pw = 0, cp = 0, pw_Dim = 0;
        RDouble dfx = 0, dfy = 0, dfz = 0, dpx = 0, dpy = 0, dpz = 0, dmx = 0, dmy = 0, dmz = 0, fvsx = 0, fvsy = 0, fvsz = 0;
        RDouble nx = 0, ny = 0, nz = 0;
        RDouble normalComponetX = 0, normalComponetY = 0, normalComponetZ = 0;
        RDouble vis = 0, txx = 0, txy = 0, txz = 0, tyx = 0, tyy = 0, tyz = 0, tzx = 0, tzy = 0, tzz = 0;
        RDouble dudx = 0, dudy = 0, dudz = 0, dvdx = 0, dvdy = 0, dvdz = 0, dwdx = 0, dwdy = 0, dwdz = 0;
        RDouble xc = 0, yc = 0, zc = 0, dx = 0, dy = 0, dz = 0, ods = 0;

        RDouble *wallDistance = grid->GetWallDist();

        RDouble *faceCP    = new RDouble[nBoundFace]();
        RDouble *faceCF    = new RDouble[nBoundFace]();
        RDouble *faceQ     = new RDouble[nBoundFace]();
        RDouble *faceYPlus = new RDouble[nBoundFace]();
        RDouble *facePW    = new RDouble[nBoundFace]();

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
                    // iFace is the face number in the set of faceIndex.
                    int iFace = *iter;
                    int le = leftCellofFace[iFace];
                    int re = rightCellofFace[iFace];

                    normalComponetX = xfn[iFace];
                    normalComponetY = yfn[iFace];
                    normalComponetZ = zfn[iFace];

                    nx  = normalComponetX * area[iFace];
                    ny  = normalComponetY * area[iFace];
                    nz  = normalComponetZ * area[iFace];

                    //! Pressure drag.
                    dx  = xfc[iFace] - xcc[le];
                    dy  = yfc[iFace] - ycc[le];
                    dz  = zfc[iFace] - zcc[le];
                    pw = p[le];
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

                    //if (viscousType > INVISCID)
                    {
                        dudx = dUdx[le];
                        dudy = dUdy[le];
                        dudz = dUdz[le];

                        dvdx = dVdx[le];
                        dvdy = dVdy[le];
                        dvdz = dVdz[le];

                        dwdx = dWdx[le];
                        dwdy = dWdy[le];
                        dwdz = dWdz[le];

                        //! Gradient correction.
                        dx  = xcc[re] - xcc[le];
                        dy  = ycc[re] - ycc[le];
                        dz  = zcc[re] - zcc[le];
                        ods = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
                        dx  = dx * ods;
                        dy  = dy * ods;
                        dz  = dz * ods;

                        CorrectGradient(u[le],u[re],dudx,dudy,dudz,dx,dy,dz,ods);
                        CorrectGradient(v[le],v[re],dvdx,dvdy,dvdz,dx,dy,dz,ods);
                        CorrectGradient(w[le],w[re],dwdx,dwdy,dwdz,dx,dy,dz,ods);

                        //! It is unnecessary to consider because the turbulence viscosity coefficient is 0.
                        vis = visl[le];

                        txx = vis * two3rd * (2.0 * dudx - dvdy - dwdz);
                        tyy = vis * two3rd * (2.0 * dvdy - dwdz - dudx);
                        tzz = vis * two3rd * (2.0 * dwdz - dudx - dvdy);
                        txy = vis * (dudy + dvdx);
                        txz = vis * (dudz + dwdx);
                        tyz = vis * (dvdz + dwdy);
                        tyx = txy;
                        tzx = txz;
                        tzy = tyz;

                        fvsx = - two * (nx * txx + ny * tyx + nz * tzx);
                        fvsy = - two * (nx * txy + ny * tyy + nz * tzy);
                        fvsz = - two * (nx * txz + ny * tyz + nz * tzz);

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

                    RDouble density = rho[le];

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
                    //if (viscousType > INVISCID)
                    {
                        yPlus = sqrt(taoWall * density) * wallDistance[le] / vis;
                    }

                    faceCP   [iFace] = cp;
                    faceCF   [iFace] = cf;
                    faceYPlus[iFace] = yPlus;
                    facePW   [iFace] = pw_Dim;
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
            if (bcType == INTERFACE)
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
                    if (lableCell[le] == 2)
                    {
                        cellCF[le] = faceCF[iFace];
                        cellCP[le] = faceCP[iFace];
                    }
                }
            }
        }

        delete [] lableCell;    lableCell = NULL;

        momEqCalculator->CommunicateAnInterfaceVar(cellCP2D[0]);
        momEqCalculator->CommunicateAnInterfaceVar(cellCF2D[0]);

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
                    for (int iNode = nodePosi[iFace]; iNode < nodePosi[iFace+1]; ++ iNode)
                    {
                        int nodeIndex = face2node[iNode];
                        RDouble dx = x[nodeIndex] - xcc[leftCellIndex];
                        RDouble dy = y[nodeIndex] - ycc[leftCellIndex];
                        RDouble dz = z[nodeIndex] - zcc[leftCellIndex];
                        RDouble dist = DISTANCE(dx, dy, dz);
                        RDouble weightTemp = 1.0 / dist;

                        nodeCP   [nodeIndex] += faceCP[iFace] * weightTemp;
                        nodeCF   [nodeIndex] += faceCF[iFace] * weightTemp;
                        nodeYPlus[nodeIndex] += faceYPlus[iFace] * weightTemp;
                        nodePW   [nodeIndex] += facePW[iFace] * weightTemp;

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
            }
        }

        //! For interface.
        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
            int bcType = bcRegion->GetBCType();

            if (bcType == INTERFACE)
            {
                vector<int> *faceIndex = bcRegion->GetFaceIndex();
                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    int iFace = *iter;
                    int cR = rightCellofFace[iFace];
                    if (ABS(cellQ[cR]) > SMALL)
                    {
                        for (int iNode = nodePosi[iFace]; iNode < nodePosi[iFace+1]; ++ iNode)
                        {
                            int nodeIndex = face2node[iNode];
                            if (nodeWeight[nodeIndex] > SMALL)
                            {
                                RDouble dx = x[nodeIndex] - xcc[cR];
                                RDouble dy = y[nodeIndex] - ycc[cR];
                                RDouble dz = z[nodeIndex] - zcc[cR];
                                RDouble dist = DISTANCE(dx, dy, dz);
                                RDouble weightTemp = 1.0 / dist;

                                nodeCP[nodeIndex] += cellCP[cR] * weightTemp;
                                nodeCF[nodeIndex] += cellCF[cR] * weightTemp;
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
            }
        }

        delete [] faceCP;     faceCP = NULL;
        delete [] faceCF;     faceCF = NULL;
        delete [] faceQ;      faceQ = NULL;
        delete [] nodeWeight; nodeWeight = NULL;
        delete [] faceYPlus;  faceYPlus = NULL;
        delete [] facePW;     facePW = NULL;
        DelPointer2(cellQ2D);
        DelPointer2(cellCP2D);
        DelPointer2(cellCF2D);

        //! Compute data on wall face center.
        RDouble *cpFaceCenter    = new RDouble[nBoundFace]();
        RDouble *cfFaceCenter    = new RDouble[nBoundFace]();
        RDouble *qFaceCenter     = new RDouble[nBoundFace]();
        RDouble *yPlusFaceCenter = new RDouble[nBoundFace]();
        RDouble *PWFaceCenter    = new RDouble[nBoundFace]();

        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            if (face2node[iFace] > 4)
            {
                RDouble valueCP    = 0.0;
                RDouble valueCF    = 0.0;
                RDouble valueQ     = 0.0;
                RDouble valueYPlus = 0.0;
                RDouble valuePW    = 0.0;

                int nodeIndex;
                for (int iNode = nodePosi[iFace]; iNode < nodePosi[iFace+1]; ++ iNode)
                {
                    nodeIndex   = face2node[iNode];
                    valueCP    += nodeCP   [nodeIndex];
                    valueCF    += nodeCF   [nodeIndex];
                    valueYPlus += nodeYPlus[nodeIndex];
                    valuePW    += nodePW   [nodeIndex];
                }
                cpFaceCenter[iFace]    = valueCP    / node_number_of_each_face[iFace];
                cfFaceCenter[iFace]    = valueCF    / node_number_of_each_face[iFace];
                qFaceCenter [iFace]    = valueQ     / node_number_of_each_face[iFace];
                yPlusFaceCenter[iFace] = valueYPlus / node_number_of_each_face[iFace];
                PWFaceCenter[iFace] = valuePW / node_number_of_each_face[iFace];
            }
        }

        //! Reconstruct grid mesh of wall.
        vector < vector < int > > face2nodelist(nBoundFace);
        vector < int > linkmap;
        GetFace2NodeList(grid, SOLID_SURFACE, linkmap, face2nodelist);

        int NumPts      = static_cast<int>(linkmap.size());
        int NumElements = static_cast<int>(face2nodelist.size());
        int NumFaces    = NumElements * 4;
        int TotalNumFaceNodes = NumElements * 4;
        int *cell2node        = new int[TotalNumFaceNodes]();

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
        RDouble *dumpYPlus = new RDouble[NumPts]();
        RDouble *dumpPW    = new RDouble[NumPts]();

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
            }
        }

        delete [] nodeQ;    nodeQ = NULL;
        delete [] nodeCP;    nodeCP = NULL;
        delete [] nodeCF;    nodeCF = NULL;
        delete [] nodeYPlus;    nodeYPlus = NULL;
        delete [] nodePW;    nodePW = NULL;
        delete [] cpFaceCenter;    cpFaceCenter = NULL;
        delete [] cfFaceCenter;    cfFaceCenter = NULL;
        delete [] qFaceCenter;    qFaceCenter = NULL;
        delete [] yPlusFaceCenter;    yPlusFaceCenter = NULL;
        delete [] PWFaceCenter;    PWFaceCenter = NULL;

        int TecioMission = 1;
        PHWrite(cdata, &TecioMission, 1);

        PHWrite(cdata, &NumPts, 1);
        PHWrite(cdata, &NumElements, 1);
        PHWrite(cdata, &NumFaces, 1);

        PHWrite(cdata, NumPts);
        PHWrite(cdata, coorX, NumPts);
        PHWrite(cdata, coorY, NumPts);
        PHWrite(cdata, coorZ, NumPts);

        int ValueLocation = 1;

        PHWrite(cdata, NumPts);
        PHWrite(cdata, ValueLocation);
        PHWrite(cdata, dumpCP, NumPts);

        PHWrite(cdata, NumPts);
        PHWrite(cdata, ValueLocation);
        PHWrite(cdata, dumpCF, NumPts);

        PHWrite(cdata, NumPts);
        PHWrite(cdata, ValueLocation);
        PHWrite(cdata, dumpYPlus, NumPts);

        ValueLocation = 1;

        PHWrite(cdata, NumPts);
        PHWrite(cdata, ValueLocation);

        PHWrite(cdata, NumPts);
        PHWrite(cdata, ValueLocation);

        PHWrite(cdata, NumPts);
        PHWrite(cdata, ValueLocation);
        PHWrite(cdata, dumpPW, NumPts);

        PHWrite(cdata, &TotalNumFaceNodes, 1);
        PHWrite(cdata, cell2node, TotalNumFaceNodes);

        delete [] coorX;    coorX = NULL;
        delete [] coorY;    coorY = NULL;
        delete [] coorZ;    coorZ = NULL;
        delete [] dumpCP;    dumpCP = NULL;
        delete [] dumpCF;    dumpCF = NULL;
        delete [] dumpYPlus;    dumpYPlus = NULL;
        delete [] dumpPW;    dumpPW = NULL;
        delete [] cell2node;    cell2node = NULL;
    }

    void IncomUnstructSolver::IncompressibleInitial(Grid * gridIn)
    {
        UnstructGrid * grid = UnstructGridCast(gridIn);
        momEqCalculator->IncompressibleInitial(grid);

        int isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
        if (isSolveTurb == 1)
        {
            string TurbEqSolverName = GlobalDataBase::GetStrParaFromDB("TurbEqSolverName");
            if (TurbEqSolverName == "KE")
            {
                IncomTurbKSolver->IncompressibleInitial(grid);
                IncomTurbEplisionSolver->IncompressibleInitial(grid);
            }
            else if (TurbEqSolverName == "SA")
            {
                IncomTurbSASolver->IncompressibleInitial(grid);
            }
        }

        int isSolveEnergy = GlobalDataBase::GetIntParaFromDB("isSolveEnergy");
        if (isSolveEnergy == 1)
        {
            IncomEnergySolver->IncompressibleInitial(grid);
        }

        int isSolveSpecies = GlobalDataBase::GetIntParaFromDB("isSolveSpecies");
        if (isSolveSpecies == 1)
        {
            int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
            for (int iGas = 0; iGas < numberOfSpeciesIncom; iGas++)
            {
                IncomSpeciesEqCalculator* SpeciesSolver = static_cast<IncomSpeciesEqCalculator *>(IncomSpeciesSolver[iGas]);
                SpeciesSolver->IncompressibleInitial(grid);
            }
        }
    }

    void IncomUnstructSolver::UpdateUnsteadyProperties(Grid *gridIn)
    {
        UnstructGrid * grid = UnstructGridCast(gridIn);
        momEqCalculator->UpdateUnsteadyProperties(grid);

        int isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
        if (isSolveTurb == 1)
        {
            string TurbEqSolverName = GlobalDataBase::GetStrParaFromDB("TurbEqSolverName");
            if (TurbEqSolverName == "KE")
            {
                IncomTurbKSolver->UpdateUnsteadyProperties(grid);
                IncomTurbEplisionSolver->UpdateUnsteadyProperties(grid);
            }
            else if (TurbEqSolverName == "SA")
            {
                IncomTurbSASolver->UpdateUnsteadyProperties(grid);
            }
        }

        int isSolveEnergy = GlobalDataBase::GetIntParaFromDB("isSolveEnergy");
        if (isSolveEnergy == 1)
        {
            IncomEnergySolver->UpdateUnsteadyProperties(grid);
        }

        int isSolveSpecies = GlobalDataBase::GetIntParaFromDB("isSolveSpecies");
        if (isSolveSpecies == 1)
        {
            int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
            for (int iGas = 0; iGas < numberOfSpeciesIncom; iGas++)
            {
                IncomSpeciesEqCalculator* SpeciesSolver = static_cast<IncomSpeciesEqCalculator *>(IncomSpeciesSolver[iGas]);
                SpeciesSolver->UpdateUnsteadyProperties(grid);
            }
        }
    }
    
    void IncomUnstructSolver::UpdateUnsteadyFlux(Grid *gridIn)
    {
        UnstructGrid * grid = UnstructGridCast(gridIn);
        int isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
        int isSolveEnergy = GlobalDataBase::GetIntParaFromDB("isSolveEnergy");
        int isSolveSpecies = GlobalDataBase::GetIntParaFromDB("isSolveSpecies");
        momEqCalculator->UpdateUnsteadyFlux(grid);

        if (isSolveTurb == 1)
        {
            string TurbEqSolverName = GlobalDataBase::GetStrParaFromDB("TurbEqSolverName");
            if (TurbEqSolverName == "KE")
            {
                IncomTurbKSolver->UpdateUnsteadyFlux(grid);
                IncomTurbEplisionSolver->UpdateUnsteadyFlux(grid);
            }
            else if (TurbEqSolverName == "SA")
            {
                IncomTurbSASolver->UpdateUnsteadyFlux(grid);
            }
        }

        if (isSolveEnergy == 1)
        {
            IncomEnergySolver->UpdateUnsteadyFlux(grid);
        }

        if (isSolveSpecies == 1)
        {
            int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
            for (int iGas = 0; iGas < numberOfSpeciesIncom; iGas++)
            {
                IncomSpeciesEqCalculator* SpeciesSolver = static_cast<IncomSpeciesEqCalculator *>(IncomSpeciesSolver[iGas]);
                SpeciesSolver->UpdateUnsteadyFlux(grid);
            }
        }
    }

    void IncomUnstructSolver::UpdateUnsteadyVariable(Grid *grid)
    {
        int isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
        int isSolveEnergy = GlobalDataBase::GetIntParaFromDB("isSolveEnergy");
        int isSolveSpecies = GlobalDataBase::GetIntParaFromDB("isSolveSpecies");
        momEqCalculator->UpdateUnsteadyVariable(grid);

        if (isSolveTurb == 1)
        {
            string TurbEqSolverName = GlobalDataBase::GetStrParaFromDB("TurbEqSolverName");
            if (TurbEqSolverName == "KE")
            {
                IncomTurbKSolver->UpdateUnsteadyVariable(grid);
                IncomTurbEplisionSolver->UpdateUnsteadyVariable(grid);
            }
            else if (TurbEqSolverName == "SA")
            {
                IncomTurbSASolver->UpdateUnsteadyVariable(grid);
            }
        }

        if (isSolveEnergy == 1)
        {
            IncomEnergySolver->UpdateUnsteadyVariable(grid);
        }

        if (isSolveSpecies == 1)
        {
            int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
            for (int iGas = 0; iGas < numberOfSpeciesIncom; iGas++)
            {
                IncomSpeciesEqCalculator* SpeciesSolver = static_cast<IncomSpeciesEqCalculator *>(IncomSpeciesSolver[iGas]);
                SpeciesSolver->UpdateUnsteadyVariable(grid);
            }
        }
    }

    void IncomUnstructSolver::GetResidual(Grid *gridIn, vector<RDouble>& res)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);

        int isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
        int isSolveEnergy = GlobalDataBase::GetIntParaFromDB("isSolveEnergy");
        int isSolveSpecies = GlobalDataBase::GetIntParaFromDB("isSolveSpecies");
        momEqCalculator->GetResidual(grid, res);

        if (isSolveTurb == 1)
        {
            string TurbEqSolverName = GlobalDataBase::GetStrParaFromDB("TurbEqSolverName");
            if (TurbEqSolverName == "KE")
            {
                IncomTurbKSolver->GetResidual(grid, res);
                IncomTurbEplisionSolver->GetResidual(grid, res);
            }
            else if (TurbEqSolverName == "SA")
            {
                IncomTurbSASolver->GetResidual(grid, res);
            }
        }

        if (isSolveEnergy == 1)
        {
            IncomEnergySolver->GetResidual(grid, res);
        }

        if (isSolveSpecies == 1)
        {
            int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
            for (int iGas = 0; iGas < numberOfSpeciesIncom; iGas++)
            {
                IncomSpeciesEqCalculator* SpeciesSolver = static_cast<IncomSpeciesEqCalculator *>(IncomSpeciesSolver[iGas]);
                SpeciesSolver->GetResidual(grid, res);
            }
        }
    }

    void IncomUnstructSolver::GetPhiName(vector<string>& phiName)
    {
        phiName.push_back("U");
        phiName.push_back("V");
        phiName.push_back("W");
        phiName.push_back("P");

        int isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
        if (isSolveTurb == 1)
        {
            string TurbEqSolverName = GlobalDataBase::GetStrParaFromDB("TurbEqSolverName");
            if (TurbEqSolverName == "KE")
            {
                phiName.push_back("TE");
                phiName.push_back("ED");
            }
            else if (TurbEqSolverName == "SA")
            {
                phiName.push_back("TE");
            }
        }

        int isSolveEnergy = GlobalDataBase::GetIntParaFromDB("isSolveEnergy");
        if (isSolveEnergy == 1)
        {
            phiName.push_back("H");
        }

        int isSolveSpecies = GlobalDataBase::GetIntParaFromDB("isSolveSpecies");
        if (isSolveSpecies == 1)
        {
            int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
            for (int iGas = 0; iGas < numberOfSpeciesIncom; iGas++)
            {
                IncomSpeciesEqCalculator* SpeciesSolver = static_cast<IncomSpeciesEqCalculator *>(IncomSpeciesSolver[iGas]);
                string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
                int solverIndex = SpeciesSolver->GetSolverIndex();
                string varName = varNameIncom[solverIndex];

                phiName.push_back(varName);
            }
        }
    }


    void IncomUnstructSolver::SolveIncomSteadyField()
    {
        UnstructGrid *grid = UnstructGridCast(GetGeometry()->GetGrid());

        int loopTimes = 1;
        string flowSolverName = GlobalDataBase::GetStrParaFromDB("flowSolverName");
        if (flowSolverName == "PISO")
        {
            loopTimes = 4;
        }

        for (int i = 0; i < loopTimes; i++)
        {
            if (i == 0)
            {
                momEqCalculator->solveMomentumEquations(grid);
            }
            else
            {
                momEqCalculator->SolveExplicitVelocityEquations(grid);
            }

            momEqCalculator->calcFaceFlux(grid);

            PPEquaCalculator->solvePPEquation(grid);

            momEqCalculator->UpdateProperties(grid);
            momEqCalculator->CorrectFaceFlux(grid);
        }

        SolveTurbulenceEquation(grid);

        SolveEnergyEquation(grid);

        SolveSpeciesEquation(grid);
    }

    void IncomUnstructSolver::SolveTurbulenceEquation(Grid *gridIn) 
    {
        int isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
        if (isSolveTurb != 1)
        {
            return;
        }
        UnstructGrid *grid = UnstructGridCast(gridIn);
        string *varNameIncom = reinterpret_cast<string *>(GlobalDataBase::GetDataPtr("varNameIncom"));
        
        string TurbEqSolverName = GlobalDataBase::GetStrParaFromDB("TurbEqSolverName");
        if (TurbEqSolverName == "KE")
        {
            IncomTurbKSolver->solveScalarEquation(grid, IDX::S_ITURBK);

            IncomTurbEplisionSolver->solveScalarEquation(grid, IDX::S_IEPSILON);
        }
        else if (TurbEqSolverName == "SA")
        {
            IncomTurbSASolver->solveScalarEquation(grid, IDX::S_ISA);
        }
    }

    void IncomUnstructSolver::SolveEnergyEquation(Grid *gridIn) 
    {
        int isSolveEnergy = GlobalDataBase::GetIntParaFromDB("isSolveEnergy");
        if (isSolveEnergy != 1)
        {
            return;
        }
        UnstructGrid *grid = UnstructGridCast(gridIn);
        IncomEnergySolver->solveScalarEquation(grid, IDX::S_ITEMP);
    }

    void IncomUnstructSolver::SolveSpeciesEquation(Grid *gridIn) 
    {
        int isSolveSpecies = GlobalDataBase::GetIntParaFromDB("isSolveSpecies");
        if (isSolveSpecies != 1)
        {
            return;
        }
        UnstructGrid *grid = UnstructGridCast(gridIn);
        string *varNameIncom = reinterpret_cast<string *>(GlobalDataBase::GetDataPtr("varNameIncom"));
        int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
        for (int iGas = 0; iGas < numberOfSpeciesIncom; iGas++)
        {
            IncomSpeciesEqCalculator* SpeciesSolver = static_cast<IncomSpeciesEqCalculator *>(IncomSpeciesSolver[iGas]);
            int *localSolverIndex = reinterpret_cast <int *> (GlobalDataBase::GetDataPtr("localSolverIndex"));
            int solverIndex = SpeciesSolver->GetSolverIndex();
            SpeciesSolver->solveScalarEquation(grid, solverIndex);
        }
    }

    void IncomUnstructSolver::SolutionIsConverged(Grid *gridIn, bool *flag)
    {
        UnstructGrid * grid = UnstructGridCast(gridIn);
        int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
        RDouble uResNow;
        RDouble vResNow;
        RDouble wResNow;
        RDouble pResNow;
        grid->GetData("UResNow", &uResNow, PHDOUBLE, 1);
        grid->GetData("VResNow", &vResNow, PHDOUBLE, 1);
        grid->GetData("WResNow", &wResNow, PHDOUBLE, 1);
        grid->GetData("PResNow", &pResNow, PHDOUBLE, 1);
        RDouble resU = GlobalDataBase::GetDoubleParaFromDB("resU");
        RDouble resV = GlobalDataBase::GetDoubleParaFromDB("resV");
        RDouble resW = GlobalDataBase::GetDoubleParaFromDB("resW");
        RDouble resP = GlobalDataBase::GetDoubleParaFromDB("resP");

        if (outnstep > 0)
        {
            if ((uResNow < resU) && (vResNow < resV) && (wResNow < resW) && (pResNow < resP))
            {
                *flag = true;
            }
            else
            {
                *flag = false;
            }
        }
        else
        {
            *flag = false;
        }
    }

    void IncomUnstructSolver::AllocateGlobalVar(Grid * gridIn)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);
        CFDSolver::AllocateGlobalVar(grid);

        int isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
        int isSolveEnergy = GlobalDataBase::GetIntParaFromDB("isSolveEnergy");
        int isSolveSpecies = GlobalDataBase::GetIntParaFromDB("isSolveSpecies");
        momEqCalculator->AllocateGlobalVar(grid);

        if (isSolveTurb == 1)
        {
            string TurbEqSolverName = GlobalDataBase::GetStrParaFromDB("TurbEqSolverName");
            if (TurbEqSolverName == "KE")
            {
                IncomTurbKSolver->AllocateGlobalVar(grid);
                IncomTurbEplisionSolver->AllocateGlobalVar(grid);
            }
            else if (TurbEqSolverName == "SA")
            {
                IncomTurbSASolver->AllocateGlobalVar(grid);
            }
        }

        if (isSolveEnergy == 1)
        {
            IncomEnergySolver->AllocateGlobalVar(grid);
        }

        if (isSolveSpecies == 1)
        {
            int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
            for (int iGas = 0; iGas < numberOfSpeciesIncom; iGas++)
            {
                IncomSpeciesEqCalculator* SpeciesSolver = static_cast<IncomSpeciesEqCalculator *>(IncomSpeciesSolver[iGas]);
                SpeciesSolver->AllocateGlobalVar(grid);
            }
        }
    }

    void IncomUnstructSolver::InitControlParameters()
    {
        controlParameters->Init();
    }

    LIB_EXPORT Param_INCompSolverUnstruct * IncomUnstructSolver::GetControlParameters() const
    {
        return static_cast< Param_INCompSolverUnstruct * > (controlParameters);
    }

    void IncomUnstructSolver::ComputePostVisualVariables(Post_Visual* postVisualization)
    {
        UnstructGrid *grid = UnstructGridCast(postVisualization->GetGrid());
        RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
        RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
        RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
        RDouble *p = reinterpret_cast<RDouble *>(grid->GetDataPtr("P"));
        RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
        RDouble *mu = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));

        postVisualization->UpdateVisualNodeVarPtr("U", u);
        postVisualization->UpdateVisualNodeVarPtr("V", v);
        postVisualization->UpdateVisualNodeVarPtr("W", w);
        postVisualization->UpdateVisualNodeVarPtr("P", p);
        postVisualization->UpdateVisualNodeVarPtr("DEN", rho);
        postVisualization->UpdateVisualNodeVarPtr("MU", mu);

        int isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
        if (isSolveTurb == 1)
        {
            string TurbEqSolverName = GlobalDataBase::GetStrParaFromDB("TurbEqSolverName");
            if (TurbEqSolverName == "KE")
            {
                RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
                RDouble *e = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
                RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
                RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
                RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
                RDouble *yplus = reinterpret_cast<RDouble *>(grid->GetDataPtr("yplus"));
                RDouble *wDist = reinterpret_cast<RDouble *>(grid->GetWallDist());

                postVisualization->UpdateVisualNodeVarPtr("VIS", vist);
                postVisualization->UpdateVisualNodeVarPtr("TE", k);
                postVisualization->UpdateVisualNodeVarPtr("ED", e);
            }
            else if (TurbEqSolverName == "SA")
            {
                RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
                RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));

                postVisualization->UpdateVisualNodeVarPtr("VIS", vist);
                postVisualization->UpdateVisualNodeVarPtr("TE", k);
            }
        }

        int isSolveEnergy = GlobalDataBase::GetIntParaFromDB("isSolveEnergy");
        if (isSolveEnergy == 1)
        {
            RDouble *Enthalpy = reinterpret_cast<RDouble *>(grid->GetDataPtr("Enthalpy"));
            RDouble *T = reinterpret_cast<RDouble *>(grid->GetDataPtr("T"));
            RDouble *cp = reinterpret_cast<RDouble *>(grid->GetDataPtr("cp"));
            RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("k"));

            postVisualization->UpdateVisualNodeVarPtr("Enthalpy", Enthalpy);
            postVisualization->UpdateVisualNodeVarPtr("T", T);   
            postVisualization->UpdateVisualNodeVarPtr("CP", cp);   
            postVisualization->UpdateVisualNodeVarPtr("K", k);   
        }

        int isSolveSpecies = GlobalDataBase::GetIntParaFromDB("isSolveSpecies");
        if (isSolveSpecies == 1)
        {
            int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
            for (int iGas = 0; iGas < numberOfSpeciesIncom; iGas++)
            {
                IncomSpeciesEqCalculator *SpeciesSolver = static_cast<IncomSpeciesEqCalculator *>(IncomSpeciesSolver[iGas]);
                string *varNameIncom = reinterpret_cast<string *>(GlobalDataBase::GetDataPtr("varNameIncom"));
                int solverIndex = SpeciesSolver->GetSolverIndex();
                RDouble *specieDataPtr = reinterpret_cast<RDouble *>(grid->GetDataPtr( varNameIncom[solverIndex]));
                postVisualization->UpdateVisualNodeVarPtr(varNameIncom[solverIndex], specieDataPtr);
            }
        }
    }

    void IncomUnstructSolver::InitFlowAsRestart() 
    {
        Grid *gridIn = this->GetGeometry()->GetGrid();
        UnstructGrid *grid = UnstructGridCast(gridIn);
        int nTotalCell = grid->GetNTotalCell();
        int nTotalFace = grid->GetNTotalFace();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal = nTotalCell + nBoundFace;

        int isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
        int isSolveEnergy = GlobalDataBase::GetIntParaFromDB("isSolveEnergy");
        int isSolveSpecies = GlobalDataBase::GetIntParaFromDB("isSolveSpecies");
        momEqCalculator->InitFlowAsRestart(grid);

        if (isSolveTurb == 1)
        {
            string TurbEqSolverName = GlobalDataBase::GetStrParaFromDB("TurbEqSolverName");
            if (TurbEqSolverName == "KE")
            {
                IncomTurbKSolver->InitFlowAsRestart(grid);
                IncomTurbEplisionSolver->InitFlowAsRestart(grid);
            }
            else if (TurbEqSolverName == "SA")
            {
                IncomTurbSASolver->InitFlowAsRestart(grid);
            }
        }

        if (isSolveEnergy == 1)
        {
            IncomEnergySolver->InitFlowAsRestart(grid);
        }

        if (isSolveSpecies == 1)
        {
            int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
            for (int iGas = 0; iGas < numberOfSpeciesIncom; iGas++)
            {
                IncomSpeciesEqCalculator* SpeciesSolver = static_cast<IncomSpeciesEqCalculator *>(IncomSpeciesSolver[iGas]);
                SpeciesSolver->InitFlowAsRestart(grid);
            }
        }
    }


    bool IncomUnstructSolver::JudgeIfRestart()
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


    void IncomUnstructSolver::InitDependentVariables()
    {
        Grid *gridIn = this->GetGeometry()->GetGrid();
        UnstructGrid *grid = UnstructGridCast(gridIn);

        ConstructLocalCell2GlobalCellMap(grid);

        SetAxEqbSolver(grid);

        IncompressibleInitial(grid);
    }

    void IncomUnstructSolver::ConstructLocalCell2GlobalCellMap(Grid *gridIn)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotalFace = grid->GetNTotalFace();
        int nTotal = nTotalCell + nBoundFace;
        int *leftCellOfFace = grid->GetLeftCellOfFace();
        int *rightCellOfFace = grid->GetRightCellOfFace();

        RDouble *neighborCellNumber = NewPointer<RDouble>(nTotal);
        RDouble *zoneNumberOfNeighborCell = NewPointer<RDouble>(nTotal);

        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            neighborCellNumber[iCell] = iCell;
        }

        for (int iFace = 0; iFace < nBoundFace; iFace++)
        {
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];
            neighborCellNumber[re] = le;
        }

        int zoneID = grid->GetZoneID();
        for (int iCell = 0; iCell < nTotal; iCell++)
        {
            zoneNumberOfNeighborCell[iCell] = zoneID;
        }

        momEqCalculator->CommunicateAnInterfaceVar(neighborCellNumber);

        momEqCalculator->CommunicateAnInterfaceVar(zoneNumberOfNeighborCell);

        int *localCell2GlobalMap = grid->GetLocalCell2GlobalMap();
        if (localCell2GlobalMap == NULL)
        {
            localCell2GlobalMap = new int[nTotal];
            grid->SetLocalCell2GlobalMap(localCell2GlobalMap);
        }

        int *localCell2LocalMap = grid->GetLocalCell2LocalMap();
        if (localCell2LocalMap == NULL)
        {
            localCell2LocalMap = new int[nTotal];
            grid->SetLocalCell2LocalMap(localCell2LocalMap);
        }

        int * cell2ZoneIDMap = grid->GetCell2ZoneIDMap();
        if (cell2ZoneIDMap == NULL)
        {
            cell2ZoneIDMap = new int[nTotal];
            grid->SetCell2ZoneIDMap(cell2ZoneIDMap);
        }

        int *globalCellNoOfGlobalZones = reinterpret_cast<int *>(GlobalDataBase::GetDataPtr("globalCellNoOfGlobalZones"));

        for (int iCell = 0; iCell < nTotal; iCell++)
        {
            int zoneID = zoneNumberOfNeighborCell[iCell];
            int cellStartOfThisZone = globalCellNoOfGlobalZones[zoneID];
            localCell2GlobalMap[iCell] = cellStartOfThisZone + neighborCellNumber[iCell];
            localCell2LocalMap[iCell] = neighborCellNumber[iCell];
            cell2ZoneIDMap[iCell] = zoneNumberOfNeighborCell[iCell];
        }

        for (int iCell = 0; iCell < nTotalCell; iCell++)
        {
            localCell2LocalMap[iCell] = -2;
        }

        UnstructBCSet **bcr = grid->GetBCRecord();
        for (int iFace = 0; iFace < nBoundFace; iFace++)
        {
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];
            int bcType = bcr[iFace]->GetKey();
            if (bcType != PHENGLEI::INTERFACE)
            {
                localCell2GlobalMap[re] = -1;
                localCell2LocalMap[re] = -1;
            }
        }

        DelPointer(neighborCellNumber);
        DelPointer(zoneNumberOfNeighborCell);
    }

    void IncomUnstructSolver::SetAxEqbSolver(Grid *gridIn)
    {
        UnstructGrid * grid = UnstructGridCast(gridIn);
        int isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
        int isSolveEnergy = GlobalDataBase::GetIntParaFromDB("isSolveEnergy");
        int isSolveSpecies = GlobalDataBase::GetIntParaFromDB("isSolveSpecies");

        AxEqualbCalculator->mathLibSolveInit(grid);

        momEqCalculator->SetAxEqualbCalculator(AxEqualbCalculator);
        PPEquaCalculator->SetAxEqualbCalculator(AxEqualbCalculator);


        if (isSolveTurb == 1)
        {
            string TurbEqSolverName = GlobalDataBase::GetStrParaFromDB("TurbEqSolverName");
            if (TurbEqSolverName == "KE")
            {
                IncomTurbKSolver->SetAxEqualbCalculator(AxEqualbCalculator);
                IncomTurbEplisionSolver->SetAxEqualbCalculator(AxEqualbCalculator);
            }
            else if (TurbEqSolverName == "SA")
            {
                IncomTurbSASolver->SetAxEqualbCalculator(AxEqualbCalculator);
            }
        }

        if (isSolveEnergy == 1)
        {
            IncomEnergySolver->SetAxEqualbCalculator(AxEqualbCalculator);
        }

        if (isSolveSpecies == 1)
        {
            int numberOfSpeciesIncom = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
            for (int iGas = 0; iGas < numberOfSpeciesIncom; iGas++)
            {
                IncomSpeciesEqCalculator* SpeciesSolver = static_cast<IncomSpeciesEqCalculator *>(IncomSpeciesSolver[iGas]);
                SpeciesSolver->SetAxEqualbCalculator(AxEqualbCalculator);
            }
        }
    }

    void IncomUnstructSolver::DumpResidual(Grid *gridIn)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);
        using namespace PHMPI;
        int currentProcessor = GetCurrentProcessorID();
        if (currentProcessor != server)
        {
            return;
        }

        int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
        int maxSimuStep = GlobalDataBase::GetIntParaFromDB("maxSimuStep");
        int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
        int intervalStepRes = GlobalDataBase::GetIntParaFromDB("intervalStepRes");
        
        bool isSubIterationDump = true;
        if (isUnsteady)
        {
            GlobalDataBase::GetData("isSubIterationDump", &isSubIterationDump, PHBOOL, 1);
        }
        
        if ((outnstep == 1 || outnstep % intervalStepRes == 0 || outnstep == maxSimuStep) && (!isUnsteady || !isSubIterationDump))
        {
            string fileName = "res.dat";
            string outputFileName = string("results/") + fileName;

            if ((outnstep == 1 && !isUnsteady))
            {
                remove(&outputFileName[0]);
            }
            else
            {
                int innstep = 0;
                GlobalDataBase::GetData("innstep", &innstep, PHINT, 1);
                if (outnstep == 1 && innstep == 1)
                {
                    remove(&outputFileName[0]);
                }
            }

            ios_base::openmode openmode = ios_base::out | ios_base::app;

            ofstream file(outputFileName.c_str(), openmode);
            ostringstream oss;

            if (!isUnsteady && outnstep == 1)
            {
                file << "Title=\"THE RESIDUAL\"\n";
                file << "Variables=\n";
                file << "\"iter\"\n";
                oss << "iter" << "    ";

                vector<string> phiName;
                this->GetPhiName(phiName);

                int numOfVar = phiName.size();
                for (int iVar = 0; iVar < numOfVar; ++iVar)
                {
                    file << "\"" << phiName[iVar] + "Res" << "\"\n";
                    oss << setw(7) << phiName[iVar] + "Res" << "    ";
                }
                file << "\"" << "WallTime" << "\"\n";
                oss << setw(7) << "WallTime" << std::endl;
            }

            if (isUnsteady)
            {
                int innstep = 0;
                GlobalDataBase::GetData("innstep", &innstep, PHINT, 1);
                if (outnstep == 1 && innstep == 1)
                {
                    file << "Title=\"THE RESIDUAL\"\n";
                    file << "Variables=\n";
                    file << "\"iter\"\n";
                    file << "\"Sub-iter\"\n";
                    oss  << "iter" << "    " << "Sub-iter" << "    ";

                    vector<string> phiName;

                    this->GetPhiName(phiName);

                    int numOfVar = phiName.size();
                    for (int iVar = 0; iVar < numOfVar; ++iVar)
                    {
                        file << "\"" << phiName[iVar] + "Res" << "\"\n";
                        oss << setw(7) << phiName[iVar] + "Res" << "    ";
                    }
                    file << "\"" << "WallTime" << "\"\n";
                    oss << setw(10) << "WallTime" << std::endl;
                }
            }

            file << setiosflags(std::ios::left);
            file << setprecision(5);
            file << setiosflags(std::ios::scientific);
            file << setiosflags(ios::showpoint);

            oss << setiosflags(std::ios::left);
            oss << setprecision(5);
            oss << setiosflags(std::ios::scientific);
            oss << setiosflags(ios::showpoint);

            if (isUnsteady)
            {
                int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
                int innerStep = GlobalDataBase::GetIntParaFromDB("innstep");
                file << setw(5) << outnstep << "    " << innerStep << "    ";
                oss << setw(5) << outnstep << "    " << innerStep << "    ";
            }
            else
            {
                file << setw(5) << outnstep << "    ";
                oss << setw(5) << outnstep << "    ";
            }

            vector<RDouble> res;

            this->GetResidual(grid, res);

            int numOfSolver = res.size();
            for (int iVar = 0; iVar < numOfSolver; ++iVar)
            {
                file << setw(1) << res[iVar] << "  ";
                oss << setw(1) << res[iVar] << "  ";
            }

            RDouble step_time = PHSPACE::TIME_SPACE::GetWallTime();
            file << setw(5) << setprecision(2) << setiosflags(ios::hex) << step_time << endl;
            oss << setw(2) << setprecision(2) << setiosflags(ios::hex) << step_time << endl;
            cout << oss.str();
            vector<RDouble>().swap(res);

        }
    }


    void IncomUnstructSolver::UpdateUnsteadyFlow(Grid *grid)
    {
        UpdateUnsteadyProperties(grid);
        UpdateUnsteadyFlux(grid);
        UpdateUnsteadyVariable(grid);
    }

}
