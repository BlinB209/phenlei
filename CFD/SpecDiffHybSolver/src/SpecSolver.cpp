#include "SpecSolver.h"
#include "SpecGrid.h"
#include "Spectrum.h"
#include "RescaleField.h"
#include "Statistics.h"
#include "LES.h"
#include "PoissonSolver.h"
#include "PHMpi.h"
#include "GlobalDataBase.h"
#include "Constants.h"
#include "IO_FileName.h"
#include "TK_Log.h"
#include "TK_Exit.h"
#include "PHIO.h"
#include "Math_BasisFunction.h"
#include "PHMatrix.h"
#include "p3dfft.h"
#include <stdio.h>

namespace PHSPACE
{

SpecSolver::SpecSolver()
{
    controlParameters = 0;
}

SpecSolver::~SpecSolver()
{
    FreeControlParameters();
}

void SpecSolver::Run()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    Init();

    Param_SpecSolver *parameters = GetControlParameters();

    int maxSimuStep = parameters->GetMaxSimuStep();
    int startStatisticStep = parameters->GetStartStatisticStep();
    int intervalStepRes = parameters->GetIntervalStepRes();
    int intervalStepStatistic = parameters->GetIntervalStepStatistic();
    int intervalStepFlow = parameters->GetIntervalStepFlow();

    int nstep = GlobalDataBase::GetIntParaFromDB("nstep");
    int nstep0 = nstep;
    do {
        ++ nstep;
        GlobalDataBase::UpdateData("nstep", &nstep, PHINT, 1);

        //if( (myid == 0) && (istep%intervalStepRes == 0) )
        //{
        //    cout << setw(20) << istep;
        //    cout << setiosflags(ios::right);
        //    cout << setprecision(10);
        //    cout << setiosflags(ios::scientific);
        //    cout << setiosflags(ios::showpoint);
        //    cout << setw(20) << real((*complexU)(2, 2, 2)) << setw(20) << imag((*complexU)(2, 2, 2)) << "\n" << endl;
        //}

        GetRHS();

        Solve();

        int nstart = GlobalDataBase::GetIntParaFromDB("nstart");
        int timeOrderOfPressureTerm = GlobalDataBase::GetIntParaFromDB("timeOrderOfPressureTerm");
        if( (!nstart) && ((nstep - nstep0) == timeOrderOfPressureTerm) )
        {
            SetTimeDataCoef();
        }
        else if( IsConvectiveTermFailed && ((nstep - nstep0) == timeOrderOfPressureTerm) )
        {
            SetTimeDataCoef();
            IsConvectiveTermFailed = false;
            //IsPressureGradientFailed = false;
        }

        //RescaleVelocityField(nstep);

        if(nstep >= startStatisticStep)
        {
            GetStatistics();
        }

        if(nstep % intervalStepRes == 0)
        {
            OutputResidual();
        }

        if(nstep % intervalStepStatistic == 0)
        {
            OutputStatisticVisual();
        }

        if(nstep % intervalStepFlow == 0)
        {
            OutputRestartData();
        }

    } while(nstep < maxSimuStep);

    WriteLogFile("Iteration Finished!\n");

    DeAllocateGlobalVariables();
}

void SpecSolver::Init()
{
    bool readFlowFile = false;
    int nstart = 0;
    if(JudgeIfRestart())
    {
        readFlowFile = true;

        nstart = 1;
    }
    GlobalDataBase::UpdateData("nstart", &nstart, PHINT, 1);

    if(nstart == 0)
    {
        WriteLogFile("Initialize Flow!\n");
    }
    else
    {
        WriteLogFile("Reading Restart File As InitFlow!\n");
    }

    InitControlParameters();

    //SpecturmData = new Spectrum();
    //SpecturmData -> InitSpectrumData();
    InitSpectrumData();

    InitGridData();

    InitTimeDependentData();

    InitGlobal();

    InitRescaleData();

    InitStatisticsData();

    InitFieldData();

    if(nstart)
    {
        InputOldTimeData();
    }

    InitLESData();

    OutputFieldName();
}

void SpecSolver::InitSpectrumData()
{
    Spectrum *SpectrumData = GetSpectrumData();

    SpectrumData = new Spectrum();

    SpectrumData -> InitSpectrumData();

    SetSpectrumData(SpectrumData);

    WriteLogFile("Spectrum data initialized!\n");
    //SpecturmData -> InitSpectrumData();
}

void SpecSolver::InitGridData()
{
    SpecGrid *GridData = GetGrid();

    GridData = new SpecGrid();

    GridData -> InitGridData();

    SetGrid(GridData);

    WriteLogFile("Grid data initialized!\n");
}

void SpecSolver::InitRescaleData()
{
    SpecGrid *GridData = GetGrid();

    RescaleFieldHIT *RescaleFieldData = GetRescaleFieldData();

    RescaleFieldData = new RescaleFieldHIT();

    RescaleFieldData -> InitRescaleData();

    SetRescaleFieldData(RescaleFieldData);

    WriteLogFile("Rescale data initialized!\n");
}

void SpecSolver::InitStatisticsData()
{
    SpecGrid *GridData = GetGrid();

    StatisticsHIT *StatisticData = GetStatisticData();

    StatisticData = new StatisticsHIT();

    StatisticData -> InitStatisticsData(GridData);

    SetStatisticData(StatisticData);

    WriteLogFile("Statistics data initialized!\n");
}

void SpecSolver::InitLESData()
{
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");

    if(viscousType > 1)
    {
        LESHIT *LESData = GetLESData();

        LESData = new LESHIT();

        SpecGrid *GridData = GetGrid();

        LESData -> InitLESData(GridData);

        SetLESData(LESData);

        WriteLogFile("LES data initialized!\n");
    }
}

void SpecSolver::InputOldTimeData()
{
    ReadConvectiveTerm();

    ReadTimeDependentInformation();

    if(IsConvectiveTermFailed)
    {
        SetInitTimeDataCoef();
    }
}

void SpecSolver::ReadConvectiveTerm()
{
    SpecGrid *GridData = GetGrid();

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    fstream file;
    ios_base::openmode openMode = ios_base::in | ios_base::binary;
    PHSPACE::OpenFile(file, nonlinearTermFile.c_str(), openMode);

    int inttmp[5];
    PHRead(file, inttmp, 5);

    int nout = 3 * (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);

    PHRead(file, &(*convectiveTerm)(istf, jstf, kstf, 1), nout);

    PHSPACE::CloseFile(file);
}

void SpecSolver::ReadTimeDependentInformation()
{
    SpecGrid *GridData = GetGrid();

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    int nstep = GlobalDataBase::GetIntParaFromDB("nstep");

    fstream file;        
    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();
    //if(myid == procIDWithKx0Ky0)
    //{
        file.open(timeInformationRestartFile.c_str(), ios_base::in);
        if((!file))
        {
            TK_Exit::FileOpenErrorExit(timeInformationRestartFile);
        }
        else
        {
            file >> nstep;
        }
        file.close();
        file.clear();

        GlobalDataBase::UpdateData("nstep", &nstep, PHINT, 1);
    //}

    //GlobalDataBase::UpdateData("nstep", &nstep, PHINT, 1);
//#ifdef PH_PARALLEL
//    MPI_Bcast(&nstep, 1, MPI_INTEGER, procIDWithKx0Ky0, MPI_COMM_WORLD); 
//#endif
}

void SpecSolver::GetRHS()
{
    GetVelocityGradient();
    GetNonlinearTerm();

    Param_SpecSolver * parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    if(viscousType > 0)
    {
        //GetViscousTerm();
    }	
}

void SpecSolver::GetVelocityGradient()
{
    SpecGrid *GridData = GetGrid();

    int ist = GridData -> iStartFourierIndex;
    int jst = GridData -> jStartFourierIndex;
    int kst = GridData -> kStartFourierIndex;
    int ied = GridData -> iEndFourierIndex;
    int jed = GridData -> jEndFourierIndex;
    int ked = GridData -> kEndFourierIndex;

    GetGradientVector( &(*complexU)(ist, jst, kst), &(*complexVelocityGradient)(ist, jst, kst, 1), &(*complexVelocityGradient)(ist, jst, kst, 2), &(*complexVelocityGradient)(ist, jst, kst, 3) );
    GetGradientVector( &(*complexV)(ist, jst, kst), &(*complexVelocityGradient)(ist, jst, kst, 4), &(*complexVelocityGradient)(ist, jst, kst, 5), &(*complexVelocityGradient)(ist, jst, kst, 6) );
    GetGradientVector( &(*complexW)(ist, jst, kst), &(*complexVelocityGradient)(ist, jst, kst, 7), &(*complexVelocityGradient)(ist, jst, kst, 8), &(*complexVelocityGradient)(ist, jst, kst, 9) );
}

void SpecSolver::GetGradientVector(PHComplex *var, PHComplex *dVardX, PHComplex *dVardY, PHComplex *dVardZ)
{
    SpecGrid *GridData = GetGrid();

    int ist = GridData -> iStartFourierIndex;
    int jst = GridData -> jStartFourierIndex;
    int kst = GridData -> kStartFourierIndex;
    int ied = GridData -> iEndFourierIndex;
    int jed = GridData -> jEndFourierIndex;
    int ked = GridData -> kEndFourierIndex;

    Complex1D *kx = GridData -> complexKX;
    Complex1D *ky = GridData -> complexKY;
    Complex1D *kz = GridData -> complexKZ;

    Range I(ist, ied);
    Range J(jst, jed);
    Range K(kst, ked);

    Complex3D var3D(var, I, J, K, neverDeleteData, fortranArray);

    Complex3D dVardXVec(dVardX, I, J, K, neverDeleteData, fortranArray);
    Complex3D dVardYVec(dVardY, I, J, K, neverDeleteData, fortranArray);
    Complex3D dVardZVec(dVardZ, I, J, K, neverDeleteData, fortranArray);

    for(int iz = ist; iz <= ied; ++iz)
    {
        for(int iy = jst; iy <= jed; ++iy)
        {
            for(int ix = kst; ix <= ked; ++ix)
            {
                dVardXVec(iz, iy, ix) = (*kx)(ix) * var3D(iz, iy, ix);
                dVardYVec(iz, iy, ix) = (*ky)(iy) * var3D(iz, iy, ix);
                dVardZVec(iz, iy, ix) = (*kz)(iy) * var3D(iz, iy, ix);
            }
        }
    }
}

void SpecSolver::GetNonlinearTerm()
{
    Param_SpecSolver *parameters = GetControlParameters();

    RDouble reynolds = parameters->GetRefReNumber();
    int viscousType = parameters->GetViscousType();
    string viscousName = parameters->GetViscousName();

    SpecGrid *GridData = GetGrid();

    int istf = GridData -> iStartFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    Range IF(istf, iedf);
    Range JF(jstf, jedf);
    Range KF(kstf, kedf);

    int istp = GridData -> iStartPhysicalIndex;
    int jstp = GridData -> jStartPhysicalIndex;
    int kstp = GridData -> kStartPhysicalIndex;
    int iedp = GridData -> iEndPhysicalIndex;
    int jedp = GridData -> jEndPhysicalIndex;
    int kedp = GridData -> kEndPhysicalIndex;

    Range IP(istp, iedp);
    Range JP(jstp, jedp);
    Range KP(kstp, kedp);

    int nXYZAll = GridData -> nXYZAll;
    RDouble dnXYZAll = 1.0 / static_cast<double>(nXYZAll);

    (*complexVelocityGradient)(IF, JF, KF, 10) = (*complexU)(IF, JF, KF);
    (*complexVelocityGradient)(IF, JF, KF, 11) = (*complexV)(IF, JF, KF);
    (*complexVelocityGradient)(IF, JF, KF, 12) = (*complexW)(IF, JF, KF);

    int nXYZFourier = GridData -> nXYZFourier;
    int nXYZPhysical = GridData -> nXYZPhysical;

    unsigned char opc2r[] = "tff";
    Cp3dfft_btran_c2r_many(reinterpret_cast<double*>(&(*complexVelocityGradient)(istf, jstf, kstf, 1)), nXYZFourier, &(*realVelocityGradient)(istp, jstp, kstp, 1), nXYZPhysical, 12, opc2r);

    (*realVelocitySquare)(IP, JP, KP, 1) = (*realVelocityGradient)(IP, JP, KP, 10) * (*realVelocityGradient)(IP, JP, KP, 10); // u*u
    (*realVelocitySquare)(IP, JP, KP, 2) = (*realVelocityGradient)(IP, JP, KP, 11) * (*realVelocityGradient)(IP, JP, KP, 10); // v*u
    (*realVelocitySquare)(IP, JP, KP, 3) = (*realVelocityGradient)(IP, JP, KP, 12) * (*realVelocityGradient)(IP, JP, KP, 10); // w*u

    (*realVelocitySquare)(IP, JP, KP, 4) = (*realVelocityGradient)(IP, JP, KP, 10) * (*realVelocityGradient)(IP, JP, KP, 11); // u*v
    (*realVelocitySquare)(IP, JP, KP, 5) = (*realVelocityGradient)(IP, JP, KP, 11) * (*realVelocityGradient)(IP, JP, KP, 11); // v*v
    (*realVelocitySquare)(IP, JP, KP, 6) = (*realVelocityGradient)(IP, JP, KP, 12) * (*realVelocityGradient)(IP, JP, KP, 11); // w*v

    (*realVelocitySquare)(IP, JP, KP, 7) = (*realVelocityGradient)(IP, JP, KP, 10) * (*realVelocityGradient)(IP, JP, KP, 12); // u*w
    (*realVelocitySquare)(IP, JP, KP, 8) = (*realVelocityGradient)(IP, JP, KP, 11) * (*realVelocityGradient)(IP, JP, KP, 12); // v*w
    (*realVelocitySquare)(IP, JP, KP, 9) = (*realVelocityGradient)(IP, JP, KP, 12) * (*realVelocityGradient)(IP, JP, KP, 12); // w*w

    (*realVelocitySquare)(IP, JP, KP, 10) = (*realVelocityGradient)(IP, JP, KP, 10) * (*realVelocityGradient)(IP, JP, KP, 1)
        + (*realVelocityGradient)(IP, JP, KP, 11) * (*realVelocityGradient)(IP, JP, KP, 2)
        + (*realVelocityGradient)(IP, JP, KP, 12) * (*realVelocityGradient)(IP, JP, KP, 3); // u*du/dx + v*du/dy +w*du/dz
    (*realVelocitySquare)(IP, JP, KP, 11) = (*realVelocityGradient)(IP, JP, KP, 10) * (*realVelocityGradient)(IP, JP, KP, 4)
        + (*realVelocityGradient)(IP, JP, KP, 11) * (*realVelocityGradient)(IP, JP, KP, 5)
        + (*realVelocityGradient)(IP, JP, KP, 12) * (*realVelocityGradient)(IP, JP, KP, 6); // u*dv/dx + v*dv/dy +w*dv/dz
    (*realVelocitySquare)(IP, JP, KP, 12) = (*realVelocityGradient)(IP, JP, KP, 10) * (*realVelocityGradient)(IP, JP, KP, 7)
        + (*realVelocityGradient)(IP, JP, KP, 11) * (*realVelocityGradient)(IP, JP, KP, 8)
        + (*realVelocityGradient)(IP, JP, KP, 12) * (*realVelocityGradient)(IP, JP, KP, 9); // u*dw/dx + v*dw/dy +w*dw/dz

    unsigned char opr2c[] = "fft";
    Cp3dfft_ftran_r2c_many(&(*realVelocitySquare)(istp, jstp, kstp, 1), nXYZPhysical, reinterpret_cast<double*>(&(*complexVelocitySquare)(istf, jstf, kstf, 1)), nXYZFourier, 12, opr2c);

    (*complexVelocitySquare)(IF, JF, KF, Range(1, 12)) = (*complexVelocitySquare)(IF, JF, KF, Range(1, 12)) * dnXYZAll;

    Complex3D *complexDivergenceOfVelocitySquare = new Complex3D(IF, JF, KF, fortranArray);

    Complex4D *convectiveTermTemp = new Complex4D(IF, JF, KF, Range(1, 3), fortranArray);

    GetDivergence(&(*complexVelocitySquare)(istf, jstf, kstf, 1), &(*complexVelocitySquare)(istf, jstf, kstf, 2), &(*complexVelocitySquare)(istf, jstf, kstf, 3), complexDivergenceOfVelocitySquare);
    (*convectiveTermTemp)(IF, JF, KF, 1) = 0.5 * ( (*complexVelocitySquare)(IF, JF, KF, 10) + (*complexDivergenceOfVelocitySquare)(IF, JF, KF) );

    GetDivergence(&(*complexVelocitySquare)(istf, jstf, kstf, 4), &(*complexVelocitySquare)(istf, jstf, kstf, 5), &(*complexVelocitySquare)(istf, jstf, kstf, 6), complexDivergenceOfVelocitySquare);
    (*convectiveTermTemp)(IF, JF, KF, 2) = 0.5 * ( (*complexVelocitySquare)(IF, JF, KF, 11) + (*complexDivergenceOfVelocitySquare)(IF, JF, KF) );

    GetDivergence(&(*complexVelocitySquare)(istf, jstf, kstf, 7), &(*complexVelocitySquare)(istf, jstf, kstf, 8), &(*complexVelocitySquare)(istf, jstf, kstf, 9), complexDivergenceOfVelocitySquare);
    (*convectiveTermTemp)(IF, JF, KF, 3) = 0.5 * ( (*complexVelocitySquare)(IF, JF, KF, 12) + (*complexDivergenceOfVelocitySquare)(IF, JF, KF) );

    if (viscousType <= 1)
    {
        realnut = 0.0;
    }
    else
    {
    }

    if(viscousType > 1)
    {
    }

    (*complexNonlinearTerm)(IF, JF, KF, Range(1, 3)) = (*nonlinearTermCoef)(0) * (*convectiveTermTemp)(IF, JF, KF, Range(1, 3)) + (*nonlinearTermCoef)(-1) * (*convectiveTerm)(IF, JF, KF, Range(1, 3));

    (*convectiveTerm)(IF, JF, KF, Range(1, 3)) = (*convectiveTermTemp)(IF, JF, KF, Range(1, 3));

    delete complexDivergenceOfVelocitySquare; complexDivergenceOfVelocitySquare = NULL;
    delete convectiveTermTemp; convectiveTermTemp = NULL;
}

void SpecSolver::GetDivergence(PHComplex *vecX, PHComplex *vecY, PHComplex *vecZ, Complex3D *divVar)
{
    int ist = GridData -> iStartFourierIndex;
    int jst = GridData -> jStartFourierIndex;
    int kst = GridData -> kStartFourierIndex;
    int ied = GridData -> iEndFourierIndex;
    int jed = GridData -> jEndFourierIndex;
    int ked = GridData -> kEndFourierIndex;

    Range I(ist, ied);
    Range J(jst, jed);
    Range K(kst, ked);

    Complex1D *kx = GridData -> complexKX;
    Complex1D *ky = GridData -> complexKY;
    Complex1D *kz = GridData -> complexKZ;

    Complex3D vecVarX(vecX, I, J, K, neverDeleteData, fortranArray);
    Complex3D vecVarY(vecY, I, J, K, neverDeleteData, fortranArray);
    Complex3D vecVarZ(vecZ, I, J, K, neverDeleteData, fortranArray);

    (*divVar)(I, J, K) = PHComplex(0.0, 0.0);

    for( int ix = kst; ix <= ked; ++ix)
    {
        for( int iy = jst; iy <= jed; ++iy)
        {
            for( int iz = ist; iz <= ied; ++iz)
            {
                (*divVar)(iz, iy, ix) = (*divVar)(iz, iy, ix) + (*kx)(ix) * vecVarX(iz, iy, ix) + (*ky)(iy) * vecVarY(iz, iy, ix) + (*kz)(iz) * vecVarZ(iz, iy, ix);
            }
        }
    }
}

void SpecSolver::Solve()
{
    Param_SpecSolver *parameters = GetControlParameters();
    RDouble reynolds = parameters->GetRefReNumber();
    RDouble timeStep = parameters->GetTimeStep();

    RDouble phi = 2.0 * reynolds / (1.0 + realnut);

    SpecGrid *GridData = GetGrid();

    int istf = GridData -> iStartFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    RDouble1D &realKX = *reinterpret_cast <RDouble1D *> (GridData -> realKX);
    RDouble1D &realKY = *reinterpret_cast <RDouble1D *> (GridData -> realKY);
    RDouble1D &realKZ = *reinterpret_cast <RDouble1D *> (GridData -> realKZ);

    Complex1D &complexKX = *reinterpret_cast <Complex1D *> (GridData -> complexKX);
    Complex1D &complexKY = *reinterpret_cast <Complex1D *> (GridData -> complexKY);
    Complex1D &complexKZ = *reinterpret_cast <Complex1D *> (GridData -> complexKZ);

    RDouble realKNyquistX = GridData -> realKNyquistX;
    RDouble realKNyquistY = GridData -> realKNyquistY;
    RDouble realKNyquistZ = GridData -> realKNyquistZ;

    const PHComplex cZero = PHComplex(0.0, 0.0);

    for(int ix = kstf; ix <= kedf; ++ix)
    {
        for(int iy = jstf; iy <= jedf; ++iy)
        {
            for(int iz = istf; iz <= iedf; ++iz)
            {
                RDouble kx = realKX(ix);
                RDouble ky = realKY(iy);
                RDouble kz = realKZ(iz);

                RDouble k2  = kx * kx + ky * ky + kz * kz;

                if(k2 < TINY)
                {
                    (*complexP)(iz, iy, ix) = cZero;
                    (*complexU)(iz, iy ,ix) = cZero;
                    (*complexV)(iz, iy, ix) = cZero;
                    (*complexW)(iz, iy, ix) = cZero;

                    continue;
                }

                if( (abs(abs(ky) - abs(realKNyquistY)) < TINY) || (abs(abs(kx) - abs(realKNyquistX)) < TINY) || (abs(abs(kz) - abs(realKNyquistZ)) < TINY) )
                {
                    (*complexP)(iz, iy, ix) = cZero;
                    (*complexU)(iz, iy ,ix) = cZero;
                    (*complexV)(iz, iy, ix) = cZero;
                    (*complexW)(iz, iy, ix) = cZero;

                    continue;
                }

                PHComplex ckx = complexKX(ix);
                PHComplex cky = complexKY(iy);
                PHComplex ckz = complexKZ(iz);

                PHComplex pLocal = (*complexP)(iz, iy, ix);

                PHComplex uLocal = -(*complexNonlinearTerm)(iz, iy, ix, 1) - ckx * pLocal + (*complexViscousTerm)(iz, iy ,ix, 1);
                PHComplex vLocal = -(*complexNonlinearTerm)(iz, iy, ix, 2) - cky * pLocal + (*complexViscousTerm)(iz, iy ,ix, 2);
                PHComplex wLocal = -(*complexNonlinearTerm)(iz, iy, ix, 3) - ckz * pLocal + (*complexViscousTerm)(iz, iy ,ix, 3);

                PoissonSolver::PoissonSolverVelocity(uLocal, timeStep, k2, phi);
                PoissonSolver::PoissonSolverVelocity(vLocal, timeStep, k2, phi);
                PoissonSolver::PoissonSolverVelocity(wLocal, timeStep, k2, phi);

                uLocal = (*complexU)(iz, iy, ix) + uLocal;
                vLocal = (*complexV)(iz, iy, ix) + vLocal;
                wLocal = (*complexW)(iz, iy, ix) + wLocal;

                PHComplex divergenceOfUStar = ckx * uLocal + cky * vLocal + ckz * wLocal;

                RDouble psi = 1.0 / (-1.0 * timeStep * k2);
                PHComplex deltaPressure = divergenceOfUStar * psi;

                // velocity projected into divergence-free space
                uLocal = uLocal - timeStep * ckx * deltaPressure;
                vLocal = vLocal - timeStep * cky * deltaPressure;
                wLocal = wLocal - timeStep * ckz * deltaPressure;

                // pressure Update
                pLocal = pLocal + deltaPressure - 1.0 / phi * divergenceOfUStar;

                // update fields
                (*complexU)(iz, iy, ix) = uLocal;
                (*complexV)(iz, iy, ix) = vLocal;
                (*complexW)(iz, iy, ix) = wLocal;
                (*complexP)(iz, iy, ix) = pLocal;
            }
        }
    }
}

void SpecSolver::InitTimeDependentData()
{
    InitTimeDependentFileName();

    AllocTimeDependentData();

    int nstart = GlobalDataBase::GetIntParaFromDB("nstart");
    if(!nstart)
    {
        SetInitTimeDataCoef();
    }
    else
    {
        SetTimeDataCoef();
    }

    WriteLogFile("Time dependent data initialized!\n");
}

void SpecSolver::InitTimeDependentFileName()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    stringstream strmMyid;
    strmMyid << myid;

    string outputdir = "./results";
    outputdir = GlobalDataBase::GetStrParaFromDB("outputdir");

    nonlinearTermFile = outputdir + "/NonlinearTerm_" + strmMyid.str() + ".rsta"; 
    nonlinearTermBackupFile = outputdir + "/NonlinearTerm_" + strmMyid.str() + "_bk.rsta";

    pressureGradientTermFile = outputdir + "/PressureGradientTerm_" + strmMyid.str() + ".rsta";
    pressureGradientTermBackupFile = outputdir + "/PressureGradientTerm_" + strmMyid.str() + "_bk.rsta";

    laplacianVelocityTermFile = outputdir + "/LaplacianVelocity_" + strmMyid.str() + ".rsta";
    laplacianVelocityTermBackupFile = outputdir + "/LaplacianVelocity_" + strmMyid.str() + "_bk.rsta";

    timeInformationRestartFile = outputdir + "/Timeinfo.rsta";
}

void SpecSolver::AllocTimeDependentData()
{
    SpecGrid *GridData = GetGrid();

    int ist = GridData -> iStartFourierIndex;
    int jst = GridData -> jStartFourierIndex;
    int kst = GridData -> kStartFourierIndex;
    int ied = GridData -> iEndFourierIndex;
    int jed = GridData -> jEndFourierIndex;
    int ked = GridData -> kEndFourierIndex;

    int timeOrderOfConvectiveTerm = GlobalDataBase::GetIntParaFromDB("timeOrderOfConvectiveTerm"); 
    Range I( -(timeOrderOfConvectiveTerm-1), 0 );
    nonlinearTermCoef = new RDouble1D( I, fortranArray );

    int timeOrderOfPressureTerm = GlobalDataBase::GetIntParaFromDB("timeOrderOfPressureTerm");
    Range J( -timeOrderOfPressureTerm, 0 );
    Range K( -timeOrderOfPressureTerm, -1);
    pressureGradientTermCoef = new RDouble1D( J, fortranArray );

    Range II( ist , ied );
    Range JJ( jst , jed );
    Range KK( kst , ked );
    Range NN( 1   , 3   );
    convectiveTerm = new Complex4D( II, JJ, KK, NN , fortranArray );
    complexLaplacianVelocityTerm = new Complex4D( II, JJ, KK, NN, fortranArray );
    complexPressureGradientTerm = new Complex5D( II, JJ, KK, NN, K, fortranArray );
}

void SpecSolver::SetInitTimeDataCoef()
{
    (*nonlinearTermCoef)( 0) = 1.0;
    (*nonlinearTermCoef)(-1) = 0.0;

    int timeOrderOfPressureTerm = GlobalDataBase::GetIntParaFromDB("timeOrderOfPressureTerm");
    for ( int i = -timeOrderOfPressureTerm ; i <= 0 ; ++i )
    {
        (*pressureGradientTermCoef)(i) = 0.0;
    }
}

void SpecSolver::SetTimeDataCoef()
{
    (*nonlinearTermCoef)( 0) = 1.5;
    (*nonlinearTermCoef)(-1) = -0.5;

    int timeOrderOfPressureTerm = GlobalDataBase::GetIntParaFromDB("timeOrderOfPressureTerm");
    if ( timeOrderOfPressureTerm == 1 )
    {
        (*pressureGradientTermCoef)( 0) = 0.0;
        (*pressureGradientTermCoef)(-1) = 1.0;
    }
    else if ( timeOrderOfPressureTerm == 2 )
    {
        (*pressureGradientTermCoef)( 0) = 0.0;
        (*pressureGradientTermCoef)(-1) = 2.0;
        (*pressureGradientTermCoef)(-2) = -1.0;
    }
    else
    {
        TK_Exit::UnexpectedVarValue( "timeOrderOfPressureTerm", timeOrderOfPressureTerm );
    }
}

void SpecSolver::InitGlobal()
{
    SpecGrid *GridData = GetGrid();

    int istf = GridData -> iStartFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kedf = GridData -> kEndFourierIndex;
    Range IF(istf, iedf);
    Range JF(jstf, jedf);
    Range KF(kstf, kedf);

    int istp = GridData -> iStartPhysicalIndex;
    int jstp = GridData -> jStartPhysicalIndex;
    int kstp = GridData -> kStartPhysicalIndex;
    int iedp = GridData -> iEndPhysicalIndex;
    int jedp = GridData -> jEndPhysicalIndex;
    int kedp = GridData -> kEndPhysicalIndex;
    Range IP(istp, iedp);
    Range JP(jstp, jedp);
    Range KP(kstp, kedp);

    int nstep = 0;
    GlobalDataBase::UpdateData("nstep", &nstep, PHINT, 1);

    //int nstart = GlobalDataBase::GetIntParaFromDB("nstart");
    //if(!nstart)
    //{
    //    IsInitialization = true;
    //}
    //else
    //{
    //    IsInitialization = false;
    //}

    IsConvectiveTermFailed = false;
    //IsPressureGradientFailed = false;

    complexVelocityGradient = new Complex4D(IF, JF, KF, Range(1, 12), fortranArray);
    realVelocityGradient = new RDouble4D(IP, JP, KP, Range(1, 12), fortranArray);
    (*complexVelocityGradient)(IF, JF, KF, Range(1, 12)) = PHComplex(0.0, 0.0);
    (*realVelocityGradient)(IP, JP, KP, Range(1, 12)) = 0.0;

    complexVelocitySquare = new Complex4D(IF, JF, KF, Range(1, 12), fortranArray);
    realVelocitySquare = new RDouble4D(IP, JP, KP, Range(1, 12), fortranArray);
    (*complexVelocitySquare)(IF, JF, KF, Range(1, 12)) = PHComplex(0.0, 0.0);
    (*realVelocitySquare)(IP, JP, KP, Range(1, 12)) = 0.0; 

    complexViscousTerm = new Complex4D(IF, JF, KF, Range(1, 3), fortranArray);
    (*complexViscousTerm)(IF, JF, KF, Range(1, 3)) = PHComplex(0.0, 0.0);
    complexNonlinearTerm = new Complex4D(IF, JF, KF, Range(1, 3), fortranArray);
    (*complexNonlinearTerm)(IF, JF, KF, Range(1, 3)) = PHComplex(0.0, 0.0);

    WriteLogFile("Global data initialized!\n");
}

void SpecSolver::InitFieldData()
{
    InitFieldFileName();

    AllocFieldData();

    GetInitField();

    WriteLogFile("Field data initialized!\n");
}

void SpecSolver::InitFieldFileName()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    stringstream strmMyid;
    strmMyid << myid;

    string outputdir = "./results";
    outputdir = GlobalDataBase::GetStrParaFromDB("outputdir");

    velocityRestartFile = outputdir + "/Velocity_" + strmMyid.str() + ".rsta";
    velocityBackupFile = outputdir + "/Velocity_" + strmMyid.str() + "_bk.rsta";
    pressureRestartFile = outputdir + "/Pressure_" + strmMyid.str() + ".rsta";
    pressureBackupFile = outputdir + "/Pressure_" + strmMyid.str() + "_bk.rsta";
    velocityNameFile = outputdir + "/Velocity.nam";
    fieldNameFile = outputdir + "/field.nam";
}

void SpecSolver::AllocFieldData()
{
    SpecGrid *GridData = GetGrid();

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex  ;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex  ;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex  ;

    int istp = GridData -> iStartPhysicalIndex;
    int jstp = GridData -> jStartPhysicalIndex;
    int kstp = GridData -> kStartPhysicalIndex;
    int iedp = GridData -> iEndPhysicalIndex;
    int jedp = GridData -> jEndPhysicalIndex;
    int kedp = GridData -> kEndPhysicalIndex;

    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;
    int kEndFourierIndexGlobal = GridData -> kEndFourierIndexGlobal;

    Range IF(istf, iedf);
    Range JF(jstf, jedf);
    Range KF(kstf, kedf);

    Range IP(istp, iedp);
    Range JP(jstp, jedp);
    Range KP(kstp, kedp);

    complexU = new Complex3D(IF, JF, KF, fortranArray);
    complexV = new Complex3D(IF, JF, KF, fortranArray);
    complexW = new Complex3D(IF, JF, KF, fortranArray);
    complexP = new Complex3D(IF, JF, KF, fortranArray);

    complexDivergenceU = new Complex3D(IF, JF, KF, fortranArray);

    (*complexU)(IF, JF, KF) = PHComplex(0.0, 0.0);
    (*complexV)(IF, JF, KF) = PHComplex(0.0, 0.0);
    (*complexW)(IF, JF, KF) = PHComplex(0.0, 0.0);
    (*complexP)(IF, JF, KF) = PHComplex(0.0, 0.0);
    (*complexDivergenceU)(IF, JF, KF) = PHComplex(0.0, 0.0);

    realnut3d = new RDouble3D(IP, JP, KP, fortranArray);
    (*realnut3d)(IP, JP, KP) = 0.0;

    realnut = 0.0;

    int nx = kEndFourierIndexGlobal - kStartFourierIndexGlobal + 1;
    int nMax = int(sqrt(3.0 * nx * nx) + 0.5);

    ek3DInit = new RDouble1D(Range(1, nMax), fortranArray);
    (*ek3DInit)(Range(1, nMax)) = 0.0;

    ek3D = new RDouble1D(Range(1, nMax), fortranArray);
    (*ek3D)(Range(1, nMax)) = 0.0;
}

void SpecSolver::GetInitField()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    if(myid == server)
    {
        WriteVelocityName();
    }

    int nstart = GlobalDataBase::GetIntParaFromDB("nstart");
    if(!nstart)
    {
        InitProfile();
    }
    else
    {
        InitFromOldField();

        int ifRescaleField = GlobalDataBase::GetIntParaFromDB("ifRescaleField");
        if(ifRescaleField)
        {
            GetInitSpectrum();
        }
    }

    OutputFieldRestartFile();
}

void SpecSolver::InitProfile()
{
    SpecGrid *GridData = GetGrid();

    int ist = GridData -> iStartFourierIndex;
    int jst = GridData -> jStartFourierIndex;
    int kst = GridData -> kStartFourierIndex;
    int ied = GridData -> iEndFourierIndex;
    int jed = GridData -> jEndFourierIndex;
    int ked = GridData -> kEndFourierIndex;

    Range I(ist, ied);
    Range J(jst, jed);
    Range K(kst, ked);

    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;
    int kEndFourierIndexGlobal = GridData -> kEndFourierIndexGlobal;

    RDouble1D &realKX = *reinterpret_cast <RDouble1D *> (GridData -> realKX);
    RDouble1D &realKY = *reinterpret_cast <RDouble1D *> (GridData -> realKY);
    RDouble1D &realKZ = *reinterpret_cast <RDouble1D *> (GridData -> realKZ);

    RDouble XL = GlobalDataBase::GetDoubleParaFromDB("XL");
    RDouble YL = GlobalDataBase::GetDoubleParaFromDB("YL");
    RDouble ZL = GlobalDataBase::GetDoubleParaFromDB("ZL");

    const double pi = 2.0 * acos(0.0);

    int nx = kEndFourierIndexGlobal - kStartFourierIndexGlobal + 1;
    int nMax = int(sqrt(3.0 * nx * nx) + 0.5);

    Range N(1, nMax);
    RDouble1D *kAll = new RDouble1D(N, fortranArray);
    RDouble1D *ekAll = new RDouble1D (N, fortranArray);

    (*kAll)(N) = 0.0;
    (*ekAll)(N) = 0.0;

    RDouble kMax = static_cast<double> (nMax - kStartFourierIndexGlobal) / XL;

    for(int i = 1; i <= nMax; ++i)
    {
        (*kAll)(i) = static_cast<double> (i - kStartFourierIndexGlobal) / XL;
    }

    (*complexP)(I, J, K) = PHComplex(0.0, 0.0);
    (*complexU)(I, J, K) = PHComplex(0.0, 0.0);
    (*complexV)(I, J, K) = PHComplex(0.0, 0.0);
    (*complexW)(I, J, K) = PHComplex(0.0, 0.0);

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    srand((unsigned int)time(NULL) + (unsigned int)myid);

    RDouble cp[3], vecr10[3], vecr20[3];

    for(int ix = kst; ix <= ked ;++ix)
    {
        for(int iy = jst; iy <= jed; ++iy)
        {
            for(int iz = ist; iz <= ied; ++iz)
            {
                RDouble kx = realKX(ix);
                RDouble ky = realKY(iy);
                RDouble kz = realKZ(iz);

                RDouble k2  = kx * kx + ky * ky + kz * kz;
                RDouble modk = sqrt(k2);

                if(k2 < TINY) continue;

                RDouble vecK[3] = {kx, ky, kz};
                RDouble vecK0[3] = {kx/modk, ky/modk, kz/modk};

                // to find a vector vec0 which is not parallel to vecK0
                RDouble vec0[3] = {1.0, 0.0, 0.0};
                RDouble vec1[3] = {0.0, 1.0, 0.0};
                RDouble vec2[3] = {0.0, 0.0, 1.0};

                if(abs(DotProduct(vec1, vecK0, 3)) < abs(DotProduct(vec0, vecK0, 3)))
                {
                    vec0[0] = 0.0;
                    vec0[1] = 1.0;
                    vec0[2] = 0.0;
                }

                if(abs(DotProduct(vec2, vecK0, 3)) < abs(DotProduct(vec0, vecK0, 3)))
                {
                    vec0[0] = 0.0;
                    vec0[1] = 0.0;
                    vec0[2] = 1.0;
                }

                // to find vecr10 and vecr20 ,so that vecr10,vecr20,and vecK0 form a local rectangular coordinate system
                CrossProduct(vecK0, vec0, cp);
                NormVector(cp, vecr10, 3);

                CrossProduct(vecK0, vecr10, vecr20);

                RDouble phi = rand()/RDouble(RAND_MAX);
                phi = 2.0 * pi * phi;

                RDouble theta1 = rand()/RDouble(RAND_MAX);
                theta1 = 2.0 * pi * theta1;

                PHComplex alpha1 = PHComplex(cos(theta1), sin(theta1)) * cos(phi);
                PHComplex vecUa[3] = {vecr10[0] * alpha1, vecr10[1] * alpha1, vecr10[2] * alpha1};

                RDouble theta2 = rand()/RDouble(RAND_MAX);
                theta2 = 2.0 * pi * theta2;

                PHComplex alpha2 = PHComplex(cos(theta2), sin(theta2)) * cos(phi);
                PHComplex vecUb[3] = {vecr20[0] * alpha2, vecr20[1] * alpha2, vecr20[2] * alpha2};

                PHComplex vecU[3] = {vecUa[0]+vecUb[0], vecUa[1]+vecUb[1], vecUa[2]+vecUb[2]};

                RDouble ek = 0.0;
                GetAssumedSpectrum(modk, ek);

                RDouble amplitude = sqrt(ek / (2.0 * pi * k2));

                (*complexU)(iz, iy, ix) = amplitude * vecU[0];
                (*complexV)(iz, iy, ix) = amplitude * vecU[1];
                (*complexW)(iz, iy, ix) = amplitude * vecU[2];
            }
        }
    }

    InitPressure();

    GetEkAll(ekAll, nMax);

    if(myid == 0)
    {
        for(int ix = 1; ix <= nMax; ++ix)
        {
            RDouble ek = 0.0;
            RDouble modk = (*kAll)(ix);

            GetAssumedSpectrum(modk, ek);

            int k = int(modk + 0.5);

            if(k < 1 || k > nMax) continue;

            (*ek3DInit)(k) = ek;
        }
    }
    

#ifdef PH_PARALLEL
    PH_Bcast(&(*ek3DInit)(1), nMax * sizeof(RDouble), 0);
#endif

    if(myid == 0)
    {
        std::ostringstream oss;

        oss << "variables = k, Ek, Ek3d_init" << endl;

        oss << setiosflags(ios::left);
        oss << setprecision(5);
        oss << setiosflags(ios::scientific);
        oss << setiosflags(ios::showpoint);

        for(int ix = 2; ix <= nMax; ++ix)
        {
            RDouble modk = (*kAll)(ix);

            int k = int(modk + 0.5);

            oss << modk << "    ";
            oss << (*ekAll)(k) << "    ";
            oss << (*ek3DInit)(k) << "    ";
            oss << "\n";
        }
        oss << endl;

        string fileName = "initEk.dat";

        fstream file;

        ios_base::openmode openmode = ios_base::out|ios_base::trunc;

        file.open(fileName.c_str(), openmode);

        WriteASCIIFile(file, oss.str());

        PHSPACE::CloseFile(file);
    }

    delete kAll; kAll = NULL;
    delete ekAll; ekAll = NULL;
}

void SpecSolver::GetAssumedSpectrum(RDouble modk, RDouble &ek)
{
    Param_SpecSolver *parameters = GetControlParameters();
    int spectrumType = parameters->GetSpectrumType();

    if(spectrumType == 3)
    {
        GetCBCSpectrum(modk, ek);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("spectrumType", spectrumType);
    }
}

void SpecSolver::GetCBCSpectrum(RDouble modk, RDouble &ek)
{
    if(abs(modk) < 10 * TINY)
    {
        ek = 0.0;

        return;
    }

    RDouble logmodk = log(modk);

    RDouble logEk = 0.0;

    RDouble *logK_CBC = SpectrumData -> logK_CBC;

    RDouble *logEK_CBC = SpectrumData -> logEK_CBC;

    RDouble *d2logEKdK2_CBC = SpectrumData -> d2logEKdK2_CBC;

    int nPoints = SpectrumData -> nPoints;

    Splint(logK_CBC, logEK_CBC, d2logEKdK2_CBC, nPoints, logmodk, logEk);

    ek = exp(logEk);
}

void SpecSolver::InitPressure()
{
    SpecGrid *GridData = GetGrid();

    RDouble1D &realKX = *reinterpret_cast <RDouble1D *> (GridData -> realKX);
    RDouble1D &realKY = *reinterpret_cast <RDouble1D *> (GridData -> realKY);
    RDouble1D &realKZ = *reinterpret_cast <RDouble1D *> (GridData -> realKZ);

    int nXYZFourier = GridData -> nXYZFourier;
    int nXYZPhysical = GridData -> nXYZPhysical;
    int nXYZAll = GridData -> nXYZAll;

    RDouble dnXYZAll = 1.0 / static_cast<double>(nXYZAll);

    int istf = GridData -> iStartFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    Range IF(istf, iedf);
    Range JF(jstf, jedf);
    Range KF(kstf, kedf);

    int istp = GridData -> iStartPhysicalIndex;
    int jstp = GridData -> jStartPhysicalIndex;
    int kstp = GridData -> kStartPhysicalIndex;
    int iedp = GridData -> iEndPhysicalIndex;
    int jedp = GridData -> jEndPhysicalIndex;
    int kedp = GridData -> kEndPhysicalIndex;

    Range IP(istp, iedp);
    Range JP(jstp, jedp);
    Range KP(kstp, kedp);

    (*complexVelocityGradient)(IF, JF, KF, 1) = (*complexU)(IF, JF, KF);
    (*complexVelocityGradient)(IF, JF, KF, 2) = (*complexV)(IF, JF, KF);
    (*complexVelocityGradient)(IF, JF, KF, 3) = (*complexW)(IF, JF, KF);

    unsigned char opc2r[] = "tff";
    Cp3dfft_btran_c2r_many(reinterpret_cast<double*>(&(*complexVelocityGradient)(istf, jstf, kstf, 1)), nXYZFourier, &(*realVelocityGradient)(istp, jstp, kstp, 1), nXYZPhysical, 3, opc2r);

    for(int ix = istp; ix <= iedp; ++ix)
    {
        for(int iy = jstp; iy <= jedp; ++iy)
        {
            for(int iz = kstp; iz <= kedp; ++iz)
            {
                RDouble u = (*realVelocityGradient)(ix, iy, iz, 1);
                RDouble v = (*realVelocityGradient)(ix, iy, iz, 2);
                RDouble w = (*realVelocityGradient)(ix, iy, iz, 3);

                (*realVelocityGradient)(ix, iy, iz, 1) = u * u;
                (*realVelocityGradient)(ix, iy, iz, 2) = u * v;
                (*realVelocityGradient)(ix, iy, iz, 3) = u * w;

                (*realVelocityGradient)(ix, iy, iz, 4) = v * u;
                (*realVelocityGradient)(ix, iy, iz, 5) = v * v;
                (*realVelocityGradient)(ix, iy, iz, 6) = v * w;

                (*realVelocityGradient)(ix, iy, iz, 7) = w * u;
                (*realVelocityGradient)(ix, iy, iz, 8) = w * v;
                (*realVelocityGradient)(ix, iy, iz, 9) = w * w;
            }
        }
    }

    unsigned char opr2c[] = "fft";
    Cp3dfft_ftran_r2c_many(&(*realVelocityGradient)(istp, jstp, kstp, 1), nXYZPhysical, reinterpret_cast<double*>(&(*complexVelocityGradient)(istf, jstf, kstf, 1)), nXYZFourier, 9, opr2c);

    (*complexVelocityGradient)(IF, JF, KF, Range(1, 9)) = (*complexVelocityGradient)(IF, JF, KF, Range(1, 9)) * dnXYZAll;

    for(int ix = kstf; ix <= kedf; ++ix)
    {
        for(int iy = jstf; iy <= jedf; ++iy)
        {
            for(int iz = istf; iz <= iedf; ++iz)
            {
                RDouble kx = realKX(ix);
                RDouble ky = realKY(iy);
                RDouble kz = realKZ(iz);

                RDouble k2  = kx * kx + ky * ky + kz * kz;

                if(k2 < TINY)
                {
                    (*complexP)(iz, iy, ix) = 0.0;
                    continue;
                }

                PHComplex uu = (*complexVelocityGradient)(iz, iy, ix, 1);
                PHComplex uv = (*complexVelocityGradient)(iz, iy, ix, 2);
                PHComplex uw = (*complexVelocityGradient)(iz, iy, ix, 3);
                PHComplex vu = (*complexVelocityGradient)(iz, iy, ix, 4);
                PHComplex vv = (*complexVelocityGradient)(iz, iy, ix, 5);
                PHComplex vw = (*complexVelocityGradient)(iz, iy, ix, 6);
                PHComplex wu = (*complexVelocityGradient)(iz, iy, ix, 7);
                PHComplex wv = (*complexVelocityGradient)(iz, iy, ix, 8);
                PHComplex ww = (*complexVelocityGradient)(iz, iy, ix, 9);

                PHComplex laplacianU = -kx * kx * uu - kx * ky * uv - kx * kz * uw - ky * kx * vu - ky * ky * vv - ky * kz * vw - kz * kx * wu - kz * ky * wv - kz * kz * ww;

                (*complexP)(iz, iy, ix) = laplacianU;
            }
        }
    }
}

void SpecSolver::GetEkAll(RDouble1D *ek3D, int nMax)
{
    SpecGrid *GridData = GetGrid();

    int ist = GridData -> iStartFourierIndex;
    int jst = GridData -> jStartFourierIndex;
    int kst = GridData -> kStartFourierIndex;
    int ied = GridData -> iEndFourierIndex;
    int jed = GridData -> jEndFourierIndex;
    int ked = GridData -> kEndFourierIndex;

    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;

    RDouble1D &realKX = *reinterpret_cast <RDouble1D *> (GridData -> realKX);
    RDouble1D &realKY = *reinterpret_cast <RDouble1D *> (GridData -> realKY);
    RDouble1D &realKZ = *reinterpret_cast <RDouble1D *> (GridData -> realKZ);

    RDouble1D *eK3DLocal = new RDouble1D(Range(1, nMax));
    Int1D *incunt = new Int1D(Range(1, nMax));
    Int1D *incuntLocal = new Int1D(Range(1, nMax));

    (*eK3DLocal)(Range(1, nMax)) = 0.0;
    (*incunt)(Range(1, nMax)) = 0;
    (*incuntLocal)(Range(1, nMax)) = 0;

    RDouble XL = GlobalDataBase::GetDoubleParaFromDB("XL");

    const RDouble pi = 2.0 * acos(0.0);

    RDouble kMax = static_cast<double> (nMax - 1) / XL;

    for( int ix = kst; ix <= ked ;++ix )
    {
        for( int iy = jst; iy <= jed; ++iy )
        {
            for( int iz = ist; iz <= ied; ++iz)
            {
                RDouble kx = realKX(ix);
                RDouble ky = realKY(iy);
                RDouble kz = realKZ(iz);

                RDouble k2  = kx * kx + ky * ky + kz * kz;
                RDouble modk = sqrt(k2);

                int kShell = int(modk + 0.5);

                if(kShell > kMax || kShell == 0) continue;

                int ib = 2;
                if(kStartFourierIndexGlobal == kst && ix == kst)
                {
                    ib = 1;
                }

                RDouble uu = real( (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix)) );
                RDouble vv = real( (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix)) );
                RDouble ww = real( (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix)) );

                // Sum between (kShell-0.5, kShell+0.5)
                (*incuntLocal)(kShell) = (*incuntLocal)(kShell) + ib;

                (*eK3DLocal)(kShell) = (*eK3DLocal)(kShell) + ib * (uu + vv + ww) * (2.0 * pi * modk * modk);
            }
        }
    }

    //MPI_Reduce(&incuntLocal(1), &incunt(1), nMax, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);

#ifdef PH_PARALLEL
    PH_Reduce(&(*incuntLocal)(1), &(*incunt)(1), nMax, PH_SUM);
    PH_Reduce(&(*eK3DLocal)(1), &(*ek3D)(1), nMax, PH_SUM);
#endif

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    if(myid == 0)
    {
        for(int i = 1; i <= nMax; ++i)
        {
            (*incunt)(i) = max(1, (*incunt)(i));
        }

        for(int kShell = 1; kShell <= nMax; ++kShell)
        {
            RDouble ek = (*ek3D)(kShell);
            RDouble ncunt = static_cast<RDouble> ((*incunt)(kShell));

            (*ek3D)(kShell) = ek / ncunt;
        }
    }

    PH_Bcast(&(*ek3D)(1), nMax * sizeof(RDouble), 0);

    delete eK3DLocal; eK3DLocal = NULL;
    delete incunt; incunt = NULL;
    delete incuntLocal; incuntLocal = NULL;
}

void SpecSolver::InitFromOldField()
{
    ReadinVelocityRestartFile();

    ReadinPressureRestartFile();
}

void SpecSolver::ReadinVelocityRestartFile()
{
    SpecGrid *GridData = GetGrid();

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    fstream velocityFile;
    ios_base::openmode openMode = ios_base::in | ios_base::binary;
    PHSPACE::OpenFile(velocityFile, velocityRestartFile.c_str(), openMode);

    int inttmp[5];
    PHRead(velocityFile, inttmp, 5);

    int nout = (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);

    PHRead(velocityFile, &(*complexU)(istf, jstf, kstf), nout);
    PHRead(velocityFile, &(*complexV)(istf, jstf, kstf), nout);
    PHRead(velocityFile, &(*complexW)(istf, jstf, kstf), nout);

    PHSPACE::CloseFile(velocityFile);
}

void SpecSolver::ReadinPressureRestartFile()
{
    SpecGrid *GridData = GetGrid();

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    fstream pressureFile;
    ios_base::openmode openMode = ios_base::in | ios_base::binary;
    PHSPACE::OpenFile(pressureFile, pressureRestartFile.c_str(), openMode);

    int inttmp[5];
    PHRead(pressureFile, inttmp, 5);

    int nout = (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);

    PHRead(pressureFile, &(*complexP)(istf, jstf, kstf), nout);

    PHSPACE::CloseFile(pressureFile);
}

void SpecSolver::GetInitSpectrum()
{
    SpecGrid *GridData = GetGrid();

    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;
    int kEndFourierIndexGlobal = GridData -> kEndFourierIndexGlobal;

    int nx = kEndFourierIndexGlobal - kStartFourierIndexGlobal + 1;
    int nMax = int(sqrt(3.0 * nx * nx) + 0.5);

    Range N(1, nMax);

    RDouble1D &kAll = *new RDouble1D(N, fortranArray);

    RDouble XL = GlobalDataBase::GetDoubleParaFromDB("XL");

    for(int i = 1; i <= nMax; ++i)
    {
        kAll(i) = static_cast<double> (i - kStartFourierIndexGlobal) / XL;
    }

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    if(myid == 0)
    {
        for(int ix = 1; ix <= nMax; ++ix)
        {
            RDouble modk = kAll(ix);
            RDouble ek = 0.0;

            GetAssumedSpectrum(modk, ek);

            (*ek3DInit)(int(modk + 0.5)) = ek;
        }
    }

    PH_Bcast(&(*ek3DInit)(1), nMax * sizeof(RDouble), 0);

    if(myid == 0)
    {
        std::ostringstream oss;

        oss << "variables = k, Ek3d_init" << endl;

        oss << setiosflags(ios::left);
        oss << setprecision(5);
        oss << setiosflags(ios::scientific);
        oss << setiosflags(ios::showpoint);

        for(int ix = 2; ix <= nMax; ++ix)
        {
            RDouble modk = kAll(ix);

            oss << modk << "    ";
            oss << (*ek3DInit)(int(modk + 0.5)) << "    ";
            oss << "\n";
        }
        oss << endl;

        string fileName = "initEk_bk.dat";

        fstream file;

        ios_base::openmode openmode = ios_base::out|ios_base::trunc;

        file.open(fileName.c_str(), openmode);

        WriteASCIIFile(file, oss.str());

        PHSPACE::CloseFile(file);
    }

    delete &kAll;
}

void SpecSolver::OutputFieldRestartFile()
{
    Param_SpecSolver *parameters = GetControlParameters();
    int intervalStepFlow = parameters->GetIntervalStepFlow();

    int nstep = GlobalDataBase::GetIntParaFromDB("nstep");
    
    int nTwrest = nstep / intervalStepFlow;

    string velocityFileName = "Velocity.rsta";
    string pressureFileName = "Pressure.rsta";
    if( nTwrest % 2 == 0 )
    {
        velocityFileName = velocityBackupFile;
        pressureFileName = pressureBackupFile;
    }
    else
    {
        velocityFileName = velocityRestartFile;
        pressureFileName = pressureRestartFile;
    }

    WriteVelocityRestartFile(velocityFileName);

    WritePressureRestartFile(pressureFileName);
}

void SpecSolver::WriteVelocityRestartFile(string velocityFileName)
{
    SpecGrid *GridData = GetGrid();

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    Int1D *fourierIndexSize = GridData -> fourierIndexSize;

    fstream velocityFile;
    ios_base::openmode openMode = ios_base::out | ios_base::binary;
    PHSPACE::OpenFile(velocityFile, velocityFileName.c_str(), openMode);

    int inttmp = 1;
    PHWrite(velocityFile, inttmp);

    PHWrite(velocityFile, &(*fourierIndexSize)(1), 3);

    inttmp = 3;
    PHWrite(velocityFile, inttmp);

    int nout = (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);

    PHWrite(velocityFile, &(*complexU)(istf, jstf, kstf), nout);
    PHWrite(velocityFile, &(*complexV)(istf, jstf, kstf), nout);
    PHWrite(velocityFile, &(*complexW)(istf, jstf, kstf), nout);

    PHSPACE::CloseFile(velocityFile);
}

void SpecSolver::WritePressureRestartFile(string pressureFileName)
{
    SpecGrid *GridData = GetGrid();

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    Int1D *fourierIndexSize = GridData -> fourierIndexSize;

    fstream pressureFile;
    ios_base::openmode openMode = ios_base::out | ios_base::binary;
    PHSPACE::OpenFile(pressureFile, pressureFileName.c_str(), openMode);

    int inttmp = 1;
    PHWrite(pressureFile, inttmp);

    PHWrite(pressureFile, &(*fourierIndexSize)(1), 3);

    PHWrite(pressureFile, inttmp);

    int nout = (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);

    PHWrite(pressureFile, &(*complexP)(istf, jstf, kstf), nout);

    PHSPACE::CloseFile(pressureFile);
}

void SpecSolver::Splint(RDouble *xa, RDouble *ya, RDouble *y2a, int n, RDouble x, RDouble &y)
{
    int klo = 0;
    int khi = n - 1;

    while((khi - klo) > 1)
    {
        int k = (khi + klo) / 2;

        if(xa[k] > x)
        {
            khi = k;
        }
        else
        {
            klo = k;
        }
    }

    RDouble h = xa[khi] - xa[klo];
    if(abs(h) < TINY)
    {
        cout << "bad xa input in splint!" << endl;

        TK_Exit::ExceptionExit( "bad xa input in splint! \n" );
    }

    RDouble a = (xa[khi] - x) / h;
    RDouble b = (x - xa[klo]) / h;

    y = a * ya[klo] + b * ya[khi] + ((pow(a, 3.0) - a) * y2a[klo] + (pow(b, 3.0) - b) * y2a[khi]) * (h * h) / 6.0;
}

/*
RDouble SpecSolver:: DotProduct(RDouble *vec1, RDouble *vec2)
{
    RDouble product;

    product = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];

    return product;
}
*/

void SpecSolver::WriteVelocityName()
{
    fstream file;
    file.open(velocityNameFile.c_str(), ios_base::out);
    if ( !file )
    {
        TK_Exit::ExceptionExit("could not open velocityNameFile\n");
    }
    file << "u;velocity \n";
    file << "v \n";
    file << "w \n";
    file.close();
    file.clear();
}

void SpecSolver::OutputFieldName()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    if(myid != 0)
    {
        return;
    }

    fstream file;
    file.open(fieldNameFile.c_str(), ios_base::out);
    if(!file)
    {
        TK_Exit::FileOpenErrorExit(fieldNameFile);
    }
    file << "p\n"
         << "u;velocity\n"
         << "v\n"
         << "w\n";
    file.close();
    file.clear();
}

const string SpecSolver::GetResidualFileName()
{
    string resSaveFile = "res.dat";
    GlobalDataBase::GetData("resSaveFile", &resSaveFile, PHSTRING, 1);
    return resSaveFile;
}

void SpecSolver::RescaleVelocityField(const int nstep)
{
    int nRescale = GlobalDataBase::GetIntParaFromDB("nRescale");

    RescaleFieldHIT *RescaleFieldData = GetRescaleFieldData();
}

void SpecSolver::GetStatistics()
{
    Param_SpecSolver *parameters = GetControlParameters();
    int statisticMethod = parameters->GetStatisticMethod();
    int viscousType = parameters->GetViscousType();

    int nStatisticalStep = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");
    ++ nStatisticalStep;
    GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);

    SpecGrid *GridData = GetGrid();

    int nY = GridData -> nY;
    int nZ = GridData -> nZ;

    int istf = GridData -> iStartFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;
    int jStartFourierIndexGlobal = GridData -> jStartFourierIndexGlobal;
    int iStartFourierIndexGlobal = GridData -> iStartFourierIndexGlobal;
    int kEndFourierIndexGlobal = GridData -> kEndFourierIndexGlobal;
    int jEndFourierIndexGlobal = GridData -> jEndFourierIndexGlobal;
    int iEndFourierIndexGlobal = GridData -> iEndFourierIndexGlobal;

    Int2D *startAndEndFourierIndexOfAll = GridData -> startAndEndFourierIndexOfAll;

    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();

    Range I(istf, iedf);
    Range J(jstf, jedf);
    Range K(kstf, kedf);
    Range XN(kStartFourierIndexGlobal, kEndFourierIndexGlobal);
    Range YN(jStartFourierIndexGlobal, jEndFourierIndexGlobal);
    Range ZN(iStartFourierIndexGlobal, iEndFourierIndexGlobal);
    Range YNHalf(jStartFourierIndexGlobal, nY/2+1);
    Range ZNHalf(kStartFourierIndexGlobal, nZ/2+1);
    Range N3(1, 3);
    Range N19(1, 19);

    StatisticsHIT *StatisticData = GetStatisticData();

    RDouble2D *energyOfUUKxAll = StatisticData -> energyOfUUKxAll;
    RDouble2D *energyOfUUKyAll = StatisticData -> energyOfUUKyAll;
    RDouble2D *energyOfUUKzAll = StatisticData -> energyOfUUKzAll;

    RDouble2D *energyOfUUKxLocal = StatisticData -> energyOfUUKxLocal;
    RDouble2D *energyOfUUKyLocal = StatisticData -> energyOfUUKyLocal;
    RDouble2D *energyOfUUKzLocal = StatisticData -> energyOfUUKzLocal;

    RDouble2D *energyOfUUKx = StatisticData -> energyOfUUKx;
    RDouble2D *energyOfUUKy = StatisticData -> energyOfUUKy;
    RDouble2D *energyOfUUKz = StatisticData -> energyOfUUKz;

    StatisticData -> GetStatisticsLocal(GridData, complexU, complexV, complexW);

    RDouble1D *statisticLocal = StatisticData -> statisticLocal;
    RDouble1D *statisticAll = StatisticData -> statisticAll;
    RDouble1D *statistic = StatisticData -> statistic;

    int nstat = 19;
    PH_Reduce(&(*statisticLocal)(1), &(*statisticAll)(1), nstat, PH_SUM);

    int nx = kEndFourierIndexGlobal - kStartFourierIndexGlobal + 1;
    int nMax = int(sqrt(3.0 * nx * nx) + 0.5);
    GetEkAll(ek3D, nMax);

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int nTProcessor = GetNumberOfProcessor();

    // !energyOfUUKxAll
#ifdef PH_PARALLEL
    if(myid == 0)
    {
#endif
        (*energyOfUUKxAll)(XN, N3) = 0.0;
        (*energyOfUUKxAll)(K, N3) = (*energyOfUUKxLocal)(K, N3);

#ifdef PH_PARALLEL
        RDouble2D *energyOfUUKxTmp = new RDouble2D(XN, N3, fortranArray);
        for(int ip = 1; ip <= (nTProcessor - 1); ++ip)
        {
            int kstf_ip = (*startAndEndFourierIndexOfAll)(3, ip);
            int kedf_ip = (*startAndEndFourierIndexOfAll)(6, ip);
            int nReceive = kedf_ip - kstf_ip + 1;

            (*energyOfUUKxTmp)(XN, N3) = 0.0;

            for(int in = 1; in <= 3; ++in)
            {
                MPI_Status status;
                MPI_Recv(&(*energyOfUUKxTmp)(kstf_ip, in), nReceive, MPI_DOUBLE, ip, ip, MPI_COMM_WORLD, &status);
                (*energyOfUUKxAll)(Range(kstf_ip, kedf_ip), in) = (*energyOfUUKxAll)(Range(kstf_ip, kedf_ip), in) + (*energyOfUUKxTmp)(Range(kstf_ip, kedf_ip), in);
            }
        }

        delete energyOfUUKxTmp; energyOfUUKxTmp = NULL;
    }
    else
    {
        int nSend = kedf - kstf + 1;
        for(int in = 1; in <= 3; ++in)
        {
            MPI_Send(&(*energyOfUUKxLocal)(kstf, in), nSend,  MPI_DOUBLE, 0, myid, MPI_COMM_WORLD);
        }
    }
#endif

    // !energyOfUUKyAll
#ifdef PH_PARALLEL
    if(myid == 0)
    {
#endif
        (*energyOfUUKyAll)(YN, N3) = 0.0;
        (*energyOfUUKyAll)(J, N3) = (*energyOfUUKyLocal)(J, N3);

#ifdef PH_PARALLEL
        RDouble2D *energyOfUUKyTmp = new RDouble2D(YN, N3, fortranArray);
        for(int ip = 1; ip <= (nTProcessor - 1); ++ip)
        {
            int jstf_ip = (*startAndEndFourierIndexOfAll)(2, ip);
            int jedf_ip = (*startAndEndFourierIndexOfAll)(5, ip);
            int nReceive = jedf_ip - jstf_ip + 1;

            (*energyOfUUKyTmp)(YN, N3) = 0.0;

            for(int in = 1; in <= 3; ++in)
            {
                MPI_Status status;
                MPI_Recv(&(*energyOfUUKyTmp)(jstf_ip, in), nReceive, MPI_DOUBLE, ip, ip, MPI_COMM_WORLD, &status);
                (*energyOfUUKyAll)(Range(jstf_ip, jedf_ip), in) = (*energyOfUUKyAll)(Range(jstf_ip, jedf_ip), in) + (*energyOfUUKyTmp)(Range(jstf_ip, jedf_ip), in);
            }
        }

        delete energyOfUUKyTmp; energyOfUUKyTmp = NULL;
    }
    else
    {
        int nSend = jedf - jstf + 1;
        for(int in = 1; in <= 3; ++in)
        {
            MPI_Send(&(*energyOfUUKyLocal)(jstf, in), nSend,  MPI_DOUBLE, 0, myid, MPI_COMM_WORLD);
        }
    }
#endif

    // !energyOfUUKzAll
#ifdef PH_PARALLEL
    if(myid == 0)
    {
#endif
        (*energyOfUUKzAll)(ZN, N3) = 0.0;
        (*energyOfUUKzAll)(I, N3) = (*energyOfUUKzLocal)(I, N3);

#ifdef PH_PARALLEL
        RDouble2D *energyOfUUKzTmp = new RDouble2D(ZN, N3, fortranArray);
        for(int ip = 1; ip <= (nTProcessor - 1); ++ip)
        {
            int istf_ip = (*startAndEndFourierIndexOfAll)(1, ip);
            int iedf_ip = (*startAndEndFourierIndexOfAll)(4, ip);
            int nReceive = iedf_ip - istf_ip + 1;

            (*energyOfUUKzTmp)(ZN, N3) = 0.0;

            for(int in = 1; in <= 3; ++in)
            {
                MPI_Status status;
                MPI_Recv(&(*energyOfUUKzTmp)(istf_ip, in), nReceive, MPI_DOUBLE, ip, ip, MPI_COMM_WORLD, &status);
                (*energyOfUUKzAll)(Range(istf_ip, iedf_ip), in) = (*energyOfUUKzAll)(Range(istf_ip, iedf_ip), in) + (*energyOfUUKzTmp)(Range(istf_ip, iedf_ip), in);
            }
        }

        delete energyOfUUKzTmp; energyOfUUKzTmp = NULL;
    }
    else
    {
        int nSend = iedf - istf + 1;
        for(int in = 1; in <= 3; ++in)
        {
            MPI_Send(&(*energyOfUUKzLocal)(istf, in), nSend,  MPI_DOUBLE, 0, myid, MPI_COMM_WORLD);
        }
    }
#endif

    if(myid != procIDWithKx0Ky0)
    {
        return;
    }

    // !energyOfUUKyT
    RDouble2D *energyOfUUKyT = new RDouble2D(YNHalf, N3, fortranArray);
    (*energyOfUUKyT)(YNHalf, N3) = 0.0;

    (*energyOfUUKyT)(jStartFourierIndexGlobal, N3) = (*energyOfUUKyAll)(jStartFourierIndexGlobal, N3);
    (*energyOfUUKyT)(nY/2+1, N3) = (*energyOfUUKyAll)(nY/2+1, N3);

    for(int iy = jStartFourierIndexGlobal+1; iy <= nY/2; ++iy)
    {
        (*energyOfUUKyT)(iy, N3) = (*energyOfUUKyAll)(iy, N3) + (*energyOfUUKyAll)(nY - iy + 2, N3);
    }

    // !energyOfUUKzT
    RDouble2D *energyOfUUKzT = new RDouble2D(ZNHalf, N3, fortranArray);
    (*energyOfUUKzT)(ZNHalf, N3) = 0.0;

    (*energyOfUUKzT)(kStartFourierIndexGlobal, N3) = (*energyOfUUKzAll)(kStartFourierIndexGlobal, N3);
    (*energyOfUUKzT)(nZ/2+1, N3) = (*energyOfUUKzAll)(nZ/2+1, N3);

    for(int ix = kStartFourierIndexGlobal+1; ix <= nZ/2; ++ix)
    {
        (*energyOfUUKzT)(ix, N3) = (*energyOfUUKzAll)(ix, N3) + (*energyOfUUKzAll)(nZ - ix + 2, N3);
    }

    if(statisticMethod == 1)
    {
        RDouble c1 = 1.0 / (nStatisticalStep * 1.0);
        RDouble c2 = 1.0 - c1;

        (*statistic)(N19) = c2 * (*statistic)(N19) + c1 * (*statisticAll)(N19);

        (*energyOfUUKx)(XN, N3) = c2 * (*energyOfUUKx)(XN, N3) + c1 * (*energyOfUUKxAll)(XN, N3);
        (*energyOfUUKy)(YNHalf, N3) = c2 * (*energyOfUUKy)(YNHalf, N3) + c1 * (*energyOfUUKyT)(YNHalf, N3);
        (*energyOfUUKz)(ZNHalf, N3) = c2 * (*energyOfUUKz)(ZNHalf, N3) + c1 * (*energyOfUUKzT)(ZNHalf, N3);

        RDouble nutAverage = StatisticData -> nutAverage;

        if(viscousType > 1)
        {
            nutAverage = c2 * nutAverage + c1 * realnut;
        }
    }
    else
    {
        TK_Exit::UnexpectedVarValue("statisticMethod", statisticMethod);
    }

    delete energyOfUUKyT; energyOfUUKyT = NULL;
    delete energyOfUUKzT; energyOfUUKzT = NULL;
}

void SpecSolver::OutputResidual()
{
    ActionKey *actkeyDumpResidual = new ActionKey();
    FillActionKey(actkeyDumpResidual, DUMP_RESIDUAL, 0);
    DumpRes(actkeyDumpResidual);
    delete actkeyDumpResidual;
}

void SpecSolver::DumpRes(ActionKey *actkey)
{
    SpecGrid *GridData = GetGrid();

    string formatRes;
    PrepareFormatedResidual(formatRes);

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();

    if(myid != procIDWithKx0Ky0)
    {
        return;
    }

    string resFile = GetResidualFileName();
    if (resFile == "")
    {
        return;
    }

    actkey->filename = resFile;
    ios_base::openmode openmode = ios_base::out|ios_base::app;
    actkey->openmode = openmode;

    if (!actkey->IsNeedOpenFile())
    {
        return;
    }

    ParallelOpenFile(actkey);

    std::ostringstream oss;
    fstream &file = *(actkey->file);

    if (IfFileEmpty(file))
    {
        oss << "Title=\"RES\"" << endl;
        oss << "Variables=" << endl;
        oss << "nstep" << endl;
        oss << "time" << endl;
        oss << "Ek" << endl;
        oss << "Re_EddyLength" << endl;
        oss << "Re_TaylorScale" << endl;
        oss << "eddyLength" << endl;
        oss << "integralLength" << endl;
        oss << "KolmogorovLength" << endl;
        oss << "lambda" << endl;
        oss << "epsilon" << endl;
        oss << "skewness" << endl;
    }

    oss << formatRes;

    WriteASCIIFile(file, oss.str());

    ParallelCloseFile(actkey);
}

void SpecSolver::PrepareFormatedResidual(string &res)
{
    Param_SpecSolver *parameters = GetControlParameters();
    RDouble reynolds = parameters->GetRefReNumber();
    RDouble timeStep = parameters->GetTimeStep();

    int istf = GridData -> iStartFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    int istp = GridData -> iStartPhysicalIndex;
    int jstp = GridData -> jStartPhysicalIndex;
    int kstp = GridData -> kStartPhysicalIndex;
    int iedp = GridData -> iEndPhysicalIndex;
    int jedp = GridData -> jEndPhysicalIndex;
    int kedp = GridData -> kEndPhysicalIndex;

    int nX = GridData -> nX;

    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();

    int nX2 = GridData -> nX2;
    int nY2 = GridData -> nY2;
    int nZ2 = GridData -> nZ2;

    Range I(istf, iedf);
    Range J(jstf, jedf);
    Range K(kstf, kedf);

    Range IP(istp, iedp);
    Range JP(jstp, jedp);
    Range KP(kstp, kedp);

    RDouble1D *statisticAll = StatisticData -> statisticAll;

    (*complexDivergenceU)(I, J, K) = (*complexVelocityGradient)(I, J, K, 1) + (*complexVelocityGradient)(I, J, K, 5) + (*complexVelocityGradient)(I, J, K, 9);

    RDouble sumDivU = 0.0;
    for(int ix = kstf; ix <= kedf; ++ix)
    {
        for(int iy = jstf; iy <= jedf; ++iy)
        {
            for(int iz = istf; iz <= iedf; ++iz)
            {
                sumDivU = sumDivU + real( (*complexDivergenceU)(iz, iy, ix) * conj((*complexDivergenceU)(iz, iy, ix)) );
            }
        }
    }

    RDouble divavg = sqrt(2.0 * sumDivU);

    RDouble dudx2 = 0.0;
    RDouble dudx3 = 0.0;  

    for(int iz = kstf; iz <= kedf; ++iz)
    {
        for(int iy = jstf; iy <= jedf; ++iy)
        {
            for(int ix = istf; ix <= iedf; ++ix)
            {
                RDouble dudx = (*realVelocityGradient)(ix, iy, iz, 1);

                dudx2 = dudx2 + dudx * dudx;
                dudx3 = dudx3 + dudx * dudx * dudx;
            }
        }
    }

    RDouble *wmpiloc = new RDouble[3];
    RDouble *wmpi = new RDouble[3];

    wmpiloc[0] = divavg;
    wmpiloc[1] = dudx2;
    wmpiloc[2] = dudx3;

    wmpi[0] = divavg;
    wmpi[1] = dudx2;
    wmpi[2] = dudx3;

    PH_Reduce(wmpiloc, wmpi, 3, PH_SUM);

    divavg = wmpi[0];
    dudx2 = wmpi[1] / (nX2 * nY2 * nZ2);
    dudx3 = wmpi[2] / (nX2 * nY2 * nZ2);

    RDouble skewness = dudx3 / (sqrt(dudx2 * dudx2 * dudx2) + TINY);

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    if(myid == procIDWithKx0Ky0)
    {
        RDouble avgGradu2 = ((*statisticAll)(11) + (*statisticAll)(15) + (*statisticAll)(19)) / 3.0;

        RDouble eps = (*statisticAll)(10) / reynolds;

        RDouble kineticEnergy = 0.5 * ((*statisticAll)(4) + (*statisticAll)(5) + (*statisticAll)(6));

        RDouble uprm = sqrt(2.0 / 3.0 * kineticEnergy);

        RDouble eddyLength = sqrt(kineticEnergy * kineticEnergy * kineticEnergy) / eps;

        RDouble reynoldsEddyLength = sqrt(kineticEnergy) * eddyLength * reynolds;

        RDouble transverseTaylorScale = uprm / sqrt(avgGradu2);

        RDouble reynoldsTaylorScale = uprm * transverseTaylorScale * reynolds;

        RDouble KolmogorovLength = sqrt(sqrt(1.0 / (eps * reynolds * reynolds *reynolds)));

        RDouble intEdk = 0.0;
        for(int kshell = 1; kshell <= nX / 2; ++kshell)
        {
            intEdk = intEdk + (*ek3D)(kshell) / static_cast<RDouble> (kshell);
        }

        const double pi = 2.0 * acos(0.0);

        RDouble integralLength = pi / 2.0 * intEdk / (uprm + TINY);

        int nstep = GlobalDataBase::GetIntParaFromDB("nstep");

        ostringstream oss;

        oss << setiosflags(ios::left);
        oss << setprecision(5);
        oss << setiosflags(ios::scientific);
        oss << setiosflags(ios::showpoint);

        oss << nstep << "    ";
        oss << static_cast<double>(nstep) * timeStep << "    ";
        oss << kineticEnergy << "    ";
        oss << reynoldsEddyLength << "    ";
        oss << reynoldsTaylorScale << "    ";
        oss << eddyLength << "    ";
        oss << integralLength << "    ";
        oss << KolmogorovLength << "    ";
        oss << transverseTaylorScale << "    ";
        oss << eps << "    ";
        oss << skewness << endl;

        res = oss.str();
    }

    delete []wmpiloc;
    delete []wmpi;
}

void SpecSolver::OutputStatisticVisual()
{
    StatisticData -> OutputStatisticVisual(GridData);
}

void SpecSolver::OutputRestartData()
{
    OutputFieldRestartFile();

    OutputOldTimeData();

    StatisticData -> OutputStatisticRestart(GridData);
}

void SpecSolver::OutputOldTimeData()
{
    Param_SpecSolver *parameters = GetControlParameters();
    int intervalStepFlow = parameters->GetIntervalStepFlow();

    int nstep = GlobalDataBase::GetIntParaFromDB("nstep");

    int inout = nstep / intervalStepFlow;

    string nonlinearTermFileName, pressureGradientTermFileName, laplacianVelocityTermFileName;
    if(inout % 2  == 0)
    {
        nonlinearTermFileName = nonlinearTermBackupFile;
        pressureGradientTermFileName = pressureGradientTermBackupFile;
        laplacianVelocityTermFileName = laplacianVelocityTermBackupFile;
    }
    else
    {
        nonlinearTermFileName = nonlinearTermFile;
        pressureGradientTermFileName = pressureGradientTermFile;
        laplacianVelocityTermFileName = laplacianVelocityTermFile;
    }

    WriteConvectiveTerm(nonlinearTermFileName);

    WriteTimeDependentInformation();
}

void SpecSolver::WriteConvectiveTerm(const string nonlinearTermFileName)
{
    //SpecDiffHybGrid *GridData = GetGrid();

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    Int1D *fourierIndexSize = GridData -> fourierIndexSize;

    fstream file;
    ios_base::openmode openMode = ios_base::out | ios_base::binary;
    PHSPACE::OpenFile(file, nonlinearTermFileName.c_str(), openMode);

    int inttmp = 1;
    PHWrite(file, inttmp);

    PHWrite(file, &(*fourierIndexSize)(1), 3);

    inttmp = 3;
    PHWrite(file, inttmp);

    int nout = 3 * (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);

    PHWrite(file, &(*convectiveTerm)(istf, jstf, kstf, 1), nout);

    PHSPACE::CloseFile(file);
}

void SpecSolver::WriteTimeDependentInformation()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();
    if(myid != procIDWithKx0Ky0)
    {
        return;
    }

    int nstep = GlobalDataBase::GetIntParaFromDB("nstep");
    fstream file;
    file.open(timeInformationRestartFile.c_str(), ios_base::out);
    if(!file)
    {
        TK_Exit::FileOpenErrorExit(timeInformationRestartFile);
    }
    file << setw(20) << nstep;
    file.close();
    file.clear();
}

void SpecSolver::FillActionKey(ActionKey *actkey, int action, int level)
{
    actkey->action     = action;
    actkey->solver     = 1;
    actkey->solverID   = 0;
    actkey->kind       = SOLVER_BASED;
    actkey->level      = level;
}

bool SpecSolver::JudgeIfRestart()
{
    string outputdir = "./results";
    outputdir = GlobalDataBase::GetStrParaFromDB("outputdir");

    string restartNSFile = outputdir + "/Velocity_0.rsta";
    //GlobalDataBase::GetData("restartNSFile", &restartNSFile, PHSTRING, 1);

    //if(PHMPI::IsParallelRun())
    //{
    //    restartNSFile = PHSPACE::AddSymbolToFileName(restartNSFile, "_", 0);
    //}

    bool restart_flag = false;

    ifstream infile(restartNSFile.c_str(), ios::in);
    if ( infile )
    {
        restart_flag = true;
        infile.close();
        infile.clear();
    }

    return restart_flag;
}

void SpecSolver::NormVector(RDouble *vecIn, RDouble *vecOut, int nDim)
{
    RDouble s = 0.0;
    for (int i = 0; i < nDim; ++ i)
    {
        s += vecIn[i] * vecIn[i];
    }

    const RDouble one = 1.0;
    RDouble os = one / (s + SMALL);

    for (int i = 0; i < nDim; ++ i)
    {
        vecOut[i] = vecIn[i] * os;
    }
}

void SpecSolver::DeAllocateGlobalVariables()
{
    delete nonlinearTermCoef; nonlinearTermCoef = NULL;
    delete pressureGradientTermCoef; pressureGradientTermCoef = NULL;
    delete convectiveTerm; convectiveTerm = NULL;
    delete complexLaplacianVelocityTerm; complexLaplacianVelocityTerm = NULL;
    delete complexPressureGradientTerm; complexPressureGradientTerm = NULL;

    delete complexViscousTerm; complexViscousTerm = NULL;
    delete complexNonlinearTerm; complexNonlinearTerm = NULL;
    delete complexVelocityGradient; complexVelocityGradient = NULL;
    delete realVelocityGradient; realVelocityGradient = NULL;
    delete complexVelocitySquare; complexVelocitySquare = NULL;
    delete realVelocitySquare; realVelocitySquare = NULL;

    delete complexU; complexU = NULL;
    delete complexV; complexV = NULL;
    delete complexW; complexW = NULL;
    delete complexP; complexP = NULL;
    delete complexDivergenceU; complexDivergenceU = NULL;
    delete realnut3d; realnut3d = NULL;

    delete ek3DInit; ek3DInit = NULL;
    delete ek3D; ek3D = NULL;

    delete SpectrumData; SpectrumData = NULL;

    delete GridData; GridData = NULL;
    
    delete RescaleFieldData; RescaleFieldData = NULL;
    
    delete StatisticData; StatisticData = NULL;

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    if(viscousType > 1)
    {
        delete LESData; LESData = NULL;
    }
}

void SpecSolver::InitControlParameters()
{
    FreeControlParameters();
    controlParameters = new Param_SpecSolver();
    controlParameters->Init();

    WriteLogFile("Control parameters initialized!\n");
}

LIB_EXPORT Param_SpecSolver *SpecSolver::GetControlParameters() const
{
    return controlParameters;
}

void SpecSolver::FreeControlParameters()
{
    if (controlParameters)
    {
        delete controlParameters;
        controlParameters = NULL;
    }
}

}