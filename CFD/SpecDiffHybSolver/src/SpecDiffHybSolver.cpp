#include "SpecDiffHybSolver.h"
#include "SpecDiffHybGrid.h"
#include "ExplicitDifferenceBoundary.h"
#include "CompactDifferenceFirstDerivative.h"
#include "CompactDifferenceSecondDerivative.h"
#include "PoissonSolver.h"
#include "CorrectWallDivergence.h"
#include "RescaleField.h"
#include "Statistics.h"
#include "LES.h"
#include "Param_SpecDiffHybSolver.h"
#include "p3dfft.h"
#include "TK_Time.h"
#include "TK_Exit.h"
#include "TK_Warning.h"
#include "PHMpi.h"
#include "PHHeader.h"
#include "GlobalDataBase.h"
#include "Constants.h"
#include "PHIO.h"
#include "IO_HDF5File.h"
#include <iostream>
#include <sstream>
#include <cstdlib>

namespace PHSPACE
{

SpecDiffHybSolver::SpecDiffHybSolver()
{
}

SpecDiffHybSolver::~SpecDiffHybSolver()
{
    FreeControlParameters();
}

void SpecDiffHybSolver::Run()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    Init();

    Param_SpecDiffHybSolver *parameters = GetControlParameters();
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

        GetRHS();

        int iVortProjection = GlobalDataBase::GetIntParaFromDB("iVortProjection");
        if( iVortProjection == 0 ) 
        {
            Solve();
        }
        else
        {
            TK_Exit::UnexpectedVarValue( "iVortProjection", iVortProjection );
        }

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

        RescaleVelocityField(nstep);

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

        //IsInitialization = false;
    } while(nstep < maxSimuStep);

    WriteLogFile("Iteration Finished!\n");

    //DeAllocateGlobalVariables();
}

void SpecDiffHybSolver::Init()
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

    InitGridData();

    InitTimeDependentData();

    InitBoundaryData();

    InitCompactDifferenceFirstDerivativeData();

    InitCompactDifferenceSecondDerivativeData();

    InitGlobal();

    InitRescaleData();

    InitFieldData();

    InitPoissonSolverData();

    if(nstart)
    {
        InputOldTimeData();
    }

    InitCorrectWallDivergenceData();

    if (!nstart)
    {
        VanishingDIVproj();
    }

    InitStatisticsData();

    InitLESData();

    OutputFieldName();
}

/*
void SpecDiffHybSolver::Init()
{
    bool readFlowFile = false;
    int nstart = 0;
    if(JudgeIfRestart())
    {
        readFlowFile = true;

        nstart = 1;
    }
    GlobalDataBase::UpdateData("nstart", &nstart, PHINT, 1);

    InitControlParameters();
    WriteLogFile("Control parameters initialized!\n");

    GridData = new SpecDiffHybGrid();
    GridData -> InitGridData();
    WriteLogFile("Grid data initialized!\n");

    InitTimeDependentData();
    WriteLogFile("Time dependent data initialized!\n");

    DiffBoundaryData = new ExplicitDifferenceBoundary();
    DiffBoundaryData -> GetSchemeCoefDiffBoundaryData();
    WriteLogFile("Boundary data initialized!\n");

    CompactDifferenceFirstDerivativeData = new CompactDifferenceFirstDerivative();
    CompactDifferenceFirstDerivativeData -> InitCompactDifferenceFirstDerivativeData(GridData, DiffBoundaryData);
    WriteLogFile("Compact scheme data for the first derivative initialized!\n");

    CompactDifferenceSecondDerivativeData = new CompactDifferenceSecondDerivative();
    CompactDifferenceSecondDerivativeData -> InitCompactDifferenceSecondDerivativeData(GridData);
    WriteLogFile("Compact scheme data for the second derivative initialized!\n");

    InitGlobal();
    WriteLogFile("Global data initialized!\n");

    CompactDifferenceFirstDerivativeData -> GetDetaDz(GridData);

    RescaleFieldData = new RescaleField();
    RescaleFieldData -> InitRescaleData( GridData );
    WriteLogFile("Rescale data initialized!\n");

    InitFieldData();
    WriteLogFile("Field data initialized!\n");

    PoissonSolverData = new PoissonSolver();
    PoissonSolverData -> InitPoissonSolverData( GridData, DiffBoundaryData, CompactDifferenceSecondDerivativeData, realnut, realnut_00 );
    WriteLogFile("Poisson solver data initialized!\n");

    if(nstart)
    {
        InputOldTimeData();
    }

    CorrectWallDivergenceData = new CorrectWallDivergence();
    CorrectWallDivergenceData -> InitCorrectWallDivergenceData( GridData, CompactDifferenceFirstDerivativeData, CompactDifferenceSecondDerivativeData, PoissonSolverData, realnut );
    WriteLogFile("Correct wall divergence data initialized!\n");

    if (!nstart)
    {
        VanishingDIVproj();
    }

    StatisticsData = new Statistics();
    StatisticsData -> InitStatisticsData(GridData);
    WriteLogFile("Statistics data initialized!\n");

    Param_SpecDiffHybSolver * parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    if(viscousType > 1)
    {
        LESData = new LES();

        LESData -> InitLESData(GridData);

        WriteLogFile("LES data initialized!\n");
    }

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();;
    if(myid == procIDWithKx0Ky0)
    {
        OutputFieldName();
    }
}
*/

void SpecDiffHybSolver::InitGridData()
{
    SpecDiffHybGrid *GridData = GetGrid();

    GridData = new SpecDiffHybGrid();

    GridData -> InitGridData();

    SetGrid(GridData);

    WriteLogFile("Grid data initialized!\n");
}

void SpecDiffHybSolver::InitTimeDependentData()
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

/*
void SpecDiffHybSolver::InitTimeDependentFileName()
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
*/

void SpecDiffHybSolver::AllocTimeDependentData()
{
    SpecDiffHybGrid *GridData = GetGrid();

    int ist = GridData -> iStartFourierIndex;
    int ied = GridData -> iEndFourierIndex  ;
    int jst = GridData -> jStartFourierIndex;
    int jed = GridData -> jEndFourierIndex  ;
    int kst = GridData -> kStartFourierIndex;
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

/*
void SpecDiffHybSolver::SetInitTimeDataCoef()
{
    (*nonlinearTermCoef)( 0) = 1.0;
    (*nonlinearTermCoef)(-1) = 0.0;

    int timeOrderOfPressureTerm = GlobalDataBase::GetIntParaFromDB("timeOrderOfPressureTerm");
    for ( int i = -timeOrderOfPressureTerm ; i <= 0 ; ++i )
    {
        (*pressureGradientTermCoef)(i) = 0.0;
    }
}

void SpecDiffHybSolver::SetTimeDataCoef()
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
*/

void SpecDiffHybSolver::InitBoundaryData()
{
    DiffBoundaryData = new ExplicitDifferenceBoundary();

    DiffBoundaryData -> GetSchemeCoefDiffBoundaryData();

    WriteLogFile("Boundary data initialized!\n");
}

void SpecDiffHybSolver::InitCompactDifferenceFirstDerivativeData()
{
    SpecDiffHybGrid *GridData = GetGrid();

    CompactDifferenceFirstDerivativeData = new CompactDifferenceFirstDerivative();

    CompactDifferenceFirstDerivativeData -> InitCompactDifferenceFirstDerivativeData(GridData, DiffBoundaryData);

    CompactDifferenceFirstDerivativeData -> GetDetaDz(GridData);

    WriteLogFile("Compact scheme data for the first derivative initialized!\n");
}

void SpecDiffHybSolver::InitCompactDifferenceSecondDerivativeData()
{
    SpecDiffHybGrid *GridData = GetGrid();

    CompactDifferenceSecondDerivativeData = new CompactDifferenceSecondDerivative();

    CompactDifferenceSecondDerivativeData -> InitCompactDifferenceSecondDerivativeData(GridData);

    WriteLogFile("Compact scheme data for the second derivative initialized!\n");
}

void SpecDiffHybSolver::InitGlobal()
{
    SpecDiffHybGrid *GridData = GetGrid();

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

    screenOutFile = "screen.out";

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

void SpecDiffHybSolver::InitRescaleData()
{
    SpecDiffHybGrid *GridData = GetGrid();

    RescaleField *RescaleFieldData = GetRescaleFieldData();

    RescaleFieldData = new RescaleField();

    RescaleFieldData -> InitRescaleData(GridData);

    SetRescaleFieldData(RescaleFieldData);

    WriteLogFile("Rescale data initialized!\n");
}

void SpecDiffHybSolver::InitStatisticsData()
{
    SpecDiffHybGrid *GridData = GetGrid();

    Statistics *StatisticData = GetStatisticData();

    StatisticData = new Statistics();

    StatisticData -> InitStatisticsData(GridData);

    SetStatisticData(StatisticData);

    WriteLogFile("Statistics data initialized!\n");
}

void SpecDiffHybSolver::InitLESData()
{
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");

    if(viscousType > 1)
    {
        LES *LESData = GetLESData();

        LESData = new LES();

        SpecDiffHybGrid *GridData = GetGrid();

        LESData -> InitLESData(GridData);

        SetLESData(LESData);

        WriteLogFile("LES data initialized!\n");
    }
}

/*
void SpecDiffHybSolver::InitFieldData()
{
    InitFieldFileName();

    AllocFieldData();

    GetInitField();

    WriteLogFile("Field data initialized!\n");
}

void SpecDiffHybSolver::InitFieldFileName()
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
*/

void SpecDiffHybSolver::AllocFieldData()
{
    SpecDiffHybGrid *GridData = GetGrid();

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

    Range IF(istf, iedf);
    Range JF(jstf, jedf);
    Range KF(kstf, kedf);

    Range IP(istp, iedp);
    Range JP(jstp, jedp);
    Range KP(kstp, kedp);

    complexU           = new Complex3D(IF, JF, KF, fortranArray);
    complexV           = new Complex3D(IF, JF, KF, fortranArray);
    complexW           = new Complex3D(IF, JF, KF, fortranArray);
    complexP           = new Complex3D(IF, JF, KF, fortranArray);
    complexDivergenceU = new Complex3D(IF, JF, KF, fortranArray);

    (*complexU)(IF, JF, KF) = PHComplex(0.0, 0.0);
    (*complexV)(IF, JF, KF) = PHComplex(0.0, 0.0);
    (*complexW)(IF, JF, KF) = PHComplex(0.0, 0.0);
    (*complexP)(IF, JF, KF) = PHComplex(0.0, 0.0);
    (*complexDivergenceU)(IF, JF, KF) = PHComplex(0.0, 0.0);

    realnut = new RDouble1D(IF, fortranArray);
    realnut_00 = new RDouble1D(IF, fortranArray);
    realnut3d = new RDouble3D(IP, JP, KP, fortranArray);
    (*realnut)(IF) = 0.0;
    (*realnut_00)(IF) = 0.0;
    (*realnut3d)(IP, JP, KP) = 0.0;
}

void SpecDiffHybSolver::GetInitField()
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
    }

    OutputFieldRestartFile();
}

void SpecDiffHybSolver::InitPoissonSolverData()
{
    SpecDiffHybGrid *GridData = GetGrid();

    PoissonSolverData = new PoissonSolver();

    PoissonSolverData -> InitPoissonSolverData(GridData, DiffBoundaryData, CompactDifferenceSecondDerivativeData, realnut, realnut_00);

    WriteLogFile("Poisson solver data initialized!\n");
}

void SpecDiffHybSolver::InitCorrectWallDivergenceData()
{
    SpecDiffHybGrid *GridData = GetGrid();

    CorrectWallDivergenceData = new CorrectWallDivergence();

    CorrectWallDivergenceData -> InitCorrectWallDivergenceData(GridData, CompactDifferenceFirstDerivativeData, CompactDifferenceSecondDerivativeData, PoissonSolverData, realnut);

    WriteLogFile("Correct wall divergence data initialized!\n");
}

/*
void SpecDiffHybSolver::WriteVelocityName()
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

void SpecDiffHybSolver::InputOldTimeData()
{
    ReadConvectiveTerm();

    ReadTimeDependentInformation();

    if(IsConvectiveTermFailed)
    {
        SetInitTimeDataCoef();
    }
}
*/

void SpecDiffHybSolver::ReadConvectiveTerm()
{
    SpecDiffHybGrid *GridData = GetGrid();

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

    /*
    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;
    Range I(istf, iedf);
    Range J(jstf, jedf);
    Range K(kstf, kedf);
    Range N(1, 3);

    fstream file;
    file.open(nonlinearTermFile.c_str(), ios_base::in|ios_base::binary);
    if(!file)
    {
        IsConvectiveTermFailed = true;
        (*convectiveTerm)(I, J, K, N) = PHComplex(0.0, 0.0);     
        TK_Warning::Warning("ReadConvectiveTerm: NonlinearTerm.rsta not exist. ConvectiveTerm will be set zero.");
    }
    else
    {
        int inttmp[5];
        int nin = 3 * (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);
        file.read(reinterpret_cast<char*>(inttmp), 5*sizeof(int));
        file.read(reinterpret_cast<char*>(&(*convectiveTerm)(istf, jstf, kstf, 1)), nin*sizeof(PHComplex));
    }
    file.close();
    file.clear();
    */
}

void SpecDiffHybSolver::ReadTimeDependentInformation()
{
    SpecDiffHybGrid *GridData = GetGrid();

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

void SpecDiffHybSolver::InitProfile()
{
    SpecDiffHybGrid *GridData = GetGrid();

    int ist = GridData -> iStartFourierIndex;
    int jst = GridData -> jStartFourierIndex;
    int kst = GridData -> kStartFourierIndex;
    int ied = GridData -> iEndFourierIndex;
    int jed = GridData -> jEndFourierIndex;
    int ked = GridData -> kEndFourierIndex;

    Range I(ist, ied);
    Range J(jst, jed);
    Range K(kst, ked);

    RDouble1D *uAvg = new RDouble1D(I, fortranArray);
    Complex1D *cfu = new Complex1D(I, fortranArray);
    Complex1D *cfv = new Complex1D(I, fortranArray);
    Complex1D *cfw = new Complex1D(I, fortranArray);
    Complex1D *cfp = new Complex1D(I, fortranArray);

    const double pi = 2.0 * acos(0.0);
    const int nXZ = (jed - jst + 1) * (ked - kst + 1);

    (*complexP)(I, J, K) = PHComplex(0.0, 0.0);
    (*complexU)(I, J, K) = PHComplex(0.0, 0.0);
    (*complexV)(I, J, K) = PHComplex(0.0, 0.0);
    (*complexW)(I, J, K) = PHComplex(0.0, 0.0);

    RDouble1D *rZ = GridData -> realZ;
    for( int iz = ist; iz <= ied; ++iz )
    {
        (*uAvg)(iz) = 1.5 * (1.0 - (*rZ)(iz) * (*rZ)(iz));
    }

    Param_SpecDiffHybSolver *parameters = GetControlParameters();
    RDouble reynolds = parameters->GetRefReNumber();

    RDouble uavgDivideByUtau = 2.960975610 + 2.439024390 * log(reynolds);
    (*uAvg)(I) = (*uAvg)(I) * uavgDivideByUtau;

    RDouble1D *realKX = GridData -> realKX;
    RDouble1D *realKY = GridData -> realKY;
    RDouble kNyquistX = GridData -> realKNyquistX;
    RDouble kNyquistY = GridData -> realKNyquistY;
    RDouble rrrmax = -1.0e9;
    RDouble ekx = 1.0;
    RDouble eky = 1.0;
    RDouble ek = 1.0;

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    srand( (unsigned int)time(NULL) + (unsigned int)myid );
    for( int ix = kst; ix <= ked ;++ix )
    {
        for( int iy = jst; iy <= jed; ++iy )
        {
            for( int iz = ist; iz <= ied; ++iz)
            {
                RDouble kx = (*realKX)(ix);
                RDouble ky = (*realKY)(iy);
                RDouble k2  = kx * kx + ky * ky;

                RDouble rd = abs( (*rZ)(iz) - (*rZ)(ist) );
                rd = min( rd, abs( (*rZ)(iz) - (*rZ)(ied) ) );
                RDouble rdplus = rd * reynolds;

                // !ekx
                RDouble kt = abs(kx);
                if( kt > TINY )
                {
                    if( rdplus > 100.0)
                    {
                        GetAssumedSpectrum(kt, rd, ekx);
                    }
                    else if( rd > TINY )
                    {
                        GetAssumedSpectrumNearWall(kt, rd, ekx);
                    }
                    else
                    {
                        ekx = 1.0;
                    }
                }
                else
                {
                    ekx = 1.0;
                }

                // !eky
                kt = abs(ky);
                if( kt > TINY )
                {
                    if( rdplus > 100.0)
                    {
                        GetAssumedSpectrum(kt, rd, eky);
                    }
                    else if( rd > TINY )
                    {
                        GetAssumedSpectrumNearWall(kt, rd, eky);
                    }
                    else
                    {
                        eky = 1.0;
                    }
                }
                else
                {
                    eky = 1.0;
                }

                ek = ekx * eky;
                RDouble rrr = ek;

                if( rrr > rrrmax )
                {
                    rrrmax = rrr;
                }

                if( (abs(kx) < kNyquistX) && (abs(ky) < kNyquistY) )
                {
                    RDouble theta = rand()/RDouble(RAND_MAX);
                    theta = 2.0 * pi * theta;
                    RDouble zz = rand()/RDouble(RAND_MAX);
                    zz = 2.0 * zz - 1.0;
                    (*complexU)(iz, iy, ix) = rrr * PHComplex(cos(theta), sin(theta)) * zz;

                    theta = rand()/RDouble(RAND_MAX);
                    theta = 2.0 * pi * theta;
                    zz = rand()/RDouble(RAND_MAX);
                    zz = 2.0 * zz - 1.0;
                    (*complexV)(iz, iy, ix) = rrr * PHComplex(cos(theta), sin(theta)) * zz;

                    theta = rand()/RDouble(RAND_MAX);
                    theta = 2.0 * pi * theta;
                    zz = rand()/RDouble(RAND_MAX);
                    zz = 2.0 * zz - 1.0;
                    (*complexW)(iz, iy, ix) = rrr * PHComplex(cos(theta), sin(theta)) * zz;
                }
                else
                {
                    continue;
                }
            }// !for iz
        }// !for iy
    }// !for ix

    if( rrrmax > 1.0/nXZ )
    {
        (*complexU)(I, J, K) = (*complexU)(I, J, K) / rrrmax * 1.0 / static_cast<RDouble>(nXZ);
        (*complexV)(I, J, K) = (*complexV)(I, J, K) / rrrmax * 1.0 / static_cast<RDouble>(nXZ);
        (*complexW)(I, J, K) = (*complexW)(I, J, K) / rrrmax * 1.0 / static_cast<RDouble>(nXZ);
    }

    // !filtering the fields:
    for(int ix = kst; ix <= ked; ++ix)
    {
        for(int iy = jst; iy <= jed; ++iy)
        {
            (*cfu)(I) = (*complexU)(I, iy, ix);
            (*cfv)(I) = (*complexV)(I, iy, ix);
            (*cfw)(I) = (*complexW)(I, iy, ix);
            (*cfp)(I) = (*complexP)(I, iy, ix);
            for(int iz = ist + 1; iz <= ied - 1; ++iz)
            {
                (*complexU)(iz, iy, ix) = ( (*cfu)(iz - 1) + (*cfu)(iz + 1) ) / 6.0 + 2.0 * (*cfu)(iz) / 3.0;
                (*complexV)(iz, iy, ix) = ( (*cfv)(iz - 1) + (*cfv)(iz + 1) ) / 6.0 + 2.0 * (*cfv)(iz) / 3.0;
                (*complexW)(iz, iy, ix) = ( (*cfw)(iz - 1) + (*cfw)(iz + 1) ) / 6.0 + 2.0 * (*cfw)(iz) / 3.0;
                (*complexP)(iz, iy, ix) = ( (*cfp)(iz - 1) + (*cfp)(iz + 1) ) / 6.0 + 2.0 * (*cfp)(iz) / 3.0;
            }
        }
    }

    //! clearing (0,0) mode:
    //int x0Fall = GridData -> x0Fall;
    //int y0Fall = GridData -> y0Fall;
    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;
    int jStartFourierIndexGlobal = GridData -> jStartFourierIndexGlobal;
    if( (kst == kStartFourierIndexGlobal) && (jst == jStartFourierIndexGlobal) )
    {
        (*complexU)(I, jst, kst) = PHComplex(0.0, 0.0);
        (*complexV)(I, jst, kst) = PHComplex(0.0, 0.0);
        (*complexW)(I, jst, kst) = PHComplex(0.0, 0.0);
        (*complexP)(I, jst, kst) = PHComplex(0.0, 0.0);
    }

    // !vanishing velocity on walls:
    (*complexU)(ist, J, K) = PHComplex(0.0, 0.0);
    (*complexV)(ist, J, K) = PHComplex(0.0, 0.0);
    (*complexW)(ist, J, K) = PHComplex(0.0, 0.0);

    (*complexU)(ied, J, K) = PHComplex(0.0, 0.0);
    (*complexV)(ied, J, K) = PHComplex(0.0, 0.0);
    (*complexW)(ied, J, K) = PHComplex(0.0, 0.0);

    // !cu(iz,jst,kst) = cmplx(uavg,0)
    if( (jst == jStartFourierIndexGlobal) && (kst == kStartFourierIndexGlobal) )
    {
        for( int iz = ist; iz <= ied; ++iz )
        {
            (*complexP)(iz, jst, kst) = PHComplex(0.0, 0.0);
            (*complexU)(iz, jst, kst) = PHComplex((*uAvg)(iz), 0.0);
        }
    }

    RescaleVelocityField(0);

    delete uAvg; uAvg = NULL;
    delete cfu; cfu = NULL;
    delete cfv; cfv = NULL;
    delete cfw; cfw = NULL;
    delete cfp; cfp = NULL;
}

void SpecDiffHybSolver::InitFromOldField()
{
    ReadinVelocityRestartFile();

    int iVortProjection = GlobalDataBase::GetIntParaFromDB("iVortProjection");
    if( iVortProjection == 0 )
    {
        ReadinPressureRestartFile();
    }
    else
    {
        return;
    }
}

void SpecDiffHybSolver::OutputFieldRestartFile()
{
    Param_SpecDiffHybSolver *parameters = GetControlParameters();
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

    int iVortProjection = GlobalDataBase::GetIntParaFromDB("iVortProjection");
    if(iVortProjection == 0)
    {
        WritePressureRestartFile(pressureFileName);
    }
}

void SpecDiffHybSolver::GetAssumedSpectrum( const double kt, const double rd, double &ek )
{
    const RDouble p0 = 2.0;
    const RDouble cL = 6.78;
    const RDouble beta = 5.2;
    const RDouble ceta = 0.4;

    RDouble kd = kt*rd;
    RDouble kL = kt * 3.0;

    Param_SpecDiffHybSolver *parameters = GetControlParameters();
    RDouble reynolds = parameters->GetRefReNumber();

    RDouble keta = kt / (0.1 * reynolds);

    RDouble ekInertia = pow(kd, -5.0/3.0);
    RDouble ekL = pow(kL / sqrt(kL * kL + cL), 5.0/3.0 + p0);
    RDouble ekEpsilon = exp( -beta * ( pow(pow(keta, 4.0) + pow(ceta, 4.0), 0.25) - ceta ) );

    ek = ekEpsilon * ekL * ekInertia;
}

void SpecDiffHybSolver::GetAssumedSpectrumNearWall( const double kt, const double rd, double &ek )
{
    const RDouble beta = 5.2;
    
    Param_SpecDiffHybSolver *parameters = GetControlParameters();
    RDouble reynolds = parameters->GetRefReNumber();

    RDouble rdplus = rd * reynolds;
    RDouble keta = kt / (0.1 * rdplus); 

    RDouble ekEpsilon = exp( -beta * keta );
    ek = pow(0.5, ekEpsilon);
}

void SpecDiffHybSolver::ReadinVelocityRestartFile()
{
    SpecDiffHybGrid *GridData = GetGrid();

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
    /*
    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    fstream file;
    file.open(velocityRestartFile.c_str(),ios_base::in|ios_base::binary);
    if ( !file )
    {
        TK_Exit::FileOpenErrorExit(velocityRestartFile);
    }

    int inttmp[5];
    file.read(reinterpret_cast<char*>(inttmp), 5*sizeof(int));
    int nin = (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);
    file.read(reinterpret_cast<char*>(&(*complexU)(istf, jstf, kstf)), nin*sizeof(PHComplex));
    file.read(reinterpret_cast<char*>(&(*complexV)(istf, jstf, kstf)), nin*sizeof(PHComplex));
    file.read(reinterpret_cast<char*>(&(*complexW)(istf, jstf, kstf)), nin*sizeof(PHComplex));
    file.close();
    file.clear();
    */
}

void SpecDiffHybSolver::ReadinPressureRestartFile()
{
    SpecDiffHybGrid *GridData = GetGrid();

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

    /*
    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    fstream file;
    file.open(pressureRestartFile.c_str(),ios_base::in|ios_base::binary);
    if ( !file )
    {
        TK_Exit::FileOpenErrorExit(pressureRestartFile);
    }

    int inttmp[5];
    file.read(reinterpret_cast<char*>(inttmp), 5*sizeof(int));
    int nin = (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);
    file.read(reinterpret_cast<char*>(&(*complexP)(istf, jstf, kstf)), nin*sizeof(PHComplex));
    file.close();
    file.clear();
    */
}

void SpecDiffHybSolver::WriteVelocityRestartFile(string velocityFileName)
{
    SpecDiffHybGrid *GridData = GetGrid();

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

    /*
    file.open(velocityFileName.c_str(), ios_base::out|ios_base::binary);
    if ( !file )
    {
        TK_Exit::FileOpenErrorExit(velocityFileName);
    }

    int inttmp = 1;
    file.write(reinterpret_cast<char*>(&inttmp), sizeof(int));
    inttmp = 3;
    Int1D *fourierIndexSize = GridData -> fourierIndexSize;
    file.write(reinterpret_cast<char*>(&(*fourierIndexSize)(1)), 3*sizeof(int));
    file.write(reinterpret_cast<char*>(&inttmp), sizeof(int));
    int nout = (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);
    file.write(reinterpret_cast<char*>(&(*complexU)(istf, jstf, kstf)), nout*sizeof(PHComplex));
    file.write(reinterpret_cast<char*>(&(*complexV)(istf, jstf, kstf)), nout*sizeof(PHComplex));
    file.write(reinterpret_cast<char*>(&(*complexW)(istf, jstf, kstf)), nout*sizeof(PHComplex));
    file.close();
    file.clear();
    */
}

void SpecDiffHybSolver::WriteVelocityRestartFile_HDF5(string velocityFileName)
{
    ActionKey *actkeyDumpVelocityRestartData = new ActionKey();

    FillActionKey(actkeyDumpVelocityRestartData, DUMP_RESTART, 0);

    actkeyDumpVelocityRestartData->filename = velocityFileName;
    if (actkeyDumpVelocityRestartData->filename == "")
    {
        return;
    }

    hid_t file;
    file = CreateHDF5File(actkeyDumpVelocityRestartData->filename);
    actkeyDumpVelocityRestartData->filepos = file;

    CreateH5VelocityRestartFile(actkeyDumpVelocityRestartData);

    DumpVelocityRestartH5(actkeyDumpVelocityRestartData);

    H5Fclose(file);
    actkeyDumpVelocityRestartData->filepos = 0;

    delete actkeyDumpVelocityRestartData;
}

void SpecDiffHybSolver::ReadVelocityRestartFile_HDF5(string velocityFileName)
{
    ActionKey *actkeyReadVelocityRestartData = new ActionKey();

    FillActionKey(actkeyReadVelocityRestartData, DUMP_RESTART, 0);

    actkeyReadVelocityRestartData->filename = velocityFileName;
    if (actkeyReadVelocityRestartData->filename == "")
    {
        return;
    }

    hid_t file;
    file = CreateHDF5File(actkeyReadVelocityRestartData->filename);
    actkeyReadVelocityRestartData->filepos = file;

    ReadVelocityRestartH5(actkeyReadVelocityRestartData);

    H5Fclose(file);
    actkeyReadVelocityRestartData->filepos = 0;

    delete actkeyReadVelocityRestartData;
}

void SpecDiffHybSolver::ReadVelocityRestartH5(ActionKey *actkey)
{
    PHComplex a = PHComplex(0.0, 0.0);

    ReadData(actkey->filepos, &a, "a");
}

void SpecDiffHybSolver::WritePressureRestartFile(string pressureFileName)
{
    SpecDiffHybGrid *GridData = GetGrid();

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

    /*
    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    fstream file;
    file.open(pressureFileName.c_str(), ios_base::out|ios_base::binary);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(pressureFileName);
    }

    int inttmp = 1;
    file.write(reinterpret_cast<char*>(&inttmp), sizeof(int));
    Int1D *fourierIndexSize = GridData -> fourierIndexSize;
    file.write(reinterpret_cast<char*>(&(*fourierIndexSize)(1)), 3*sizeof(int));
    file.write(reinterpret_cast<char*>(&inttmp), sizeof(int));
    int nout = (iedf - istf + 1) * (jedf - jstf + 1) * (kedf - kstf + 1);
    file.write(reinterpret_cast<char*>(&(*complexP)(istf, jstf, kstf)), nout*sizeof(PHComplex));
    file.close();
    file.clear();
    */
}

void SpecDiffHybSolver::CreateH5VelocityRestartFile(ActionKey *actkey)
{
}

void SpecDiffHybSolver::DumpVelocityRestartH5(ActionKey *actkey)
{
}

void SpecDiffHybSolver::GetRHS()
{
    GetVelocityGradient();
    GetNonlinearTerm();

    int iVortProjection = GlobalDataBase::GetIntParaFromDB("iVortProjection");

    Param_SpecDiffHybSolver * parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    if( (viscousType > 0) || (iVortProjection !=0 ) )
    {
        GetViscousTerm();
    }		
}

void SpecDiffHybSolver::Solve()
{
    SpecDiffHybGrid *GridData = GetGrid();

    Param_SpecDiffHybSolver *parameters = GetControlParameters();
    RDouble timeStep = parameters->GetTimeStep();
    RDouble pressureGradient = parameters->GetPressureGradient();
    RDouble reynolds = parameters->GetRefReNumber();

    const PHComplex cZero = PHComplex(0.0, 0.0);
    const PHComplex cOne = PHComplex(1.0, 0.0);
    const PHComplex cOneI = PHComplex(0.0, 1.0);

    RDouble2D *coef_LHS = CompactDifferenceSecondDerivativeData -> lHSCoefMatrix;
    RDouble2D *coef_RHS = CompactDifferenceSecondDerivativeData -> rHSCoefMatrix;
    Int1D *lowerBound = CompactDifferenceSecondDerivativeData -> lowerBoundofRHSCoefMatrix;
    Int1D *upperBound = CompactDifferenceSecondDerivativeData -> upperBoundofRHSCoefMatrix;
    RDouble4D *coefVelocity = PoissonSolverData -> coefMatrixVelocityDirichlet;

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex  ;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex  ;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex  ;
    int nSolve = iedf - istf;
    Range IF(istf, iedf);
    
    RDouble1D *realPHI = new RDouble1D(IF , fortranArray);
    (*realPHI)(IF) = 2.0 * reynolds/( 1.0 + (*realnut)(IF) );

    Complex1D *complexULocal = new Complex1D(IF , fortranArray);
    Complex1D *complexVLocal = new Complex1D(IF , fortranArray);
    Complex1D *complexWLocal = new Complex1D(IF , fortranArray);
    Complex1D *complexPLocal = new Complex1D(IF , fortranArray);
    Complex1D *cR = new Complex1D(IF, fortranArray);
    RDouble2D *rC = new RDouble2D(IF, Range(-2, 2), fortranArray);
    RDouble1D *realKX = GridData -> realKX;
    RDouble1D *realKY = GridData -> realKY;
    Complex1D *complexKX = GridData -> complexKX;
    Complex1D *complexKY = GridData -> complexKY;
    //int x0Fall = GridData -> x0Fall;
    //int y0Fall = GridData -> y0Fall;
    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;
    int jStartFourierIndexGlobal = GridData -> jStartFourierIndexGlobal;

    if( (kstf == kStartFourierIndexGlobal) && (jstf == jStartFourierIndexGlobal) )
    {
        int imode00fix = GlobalDataBase::GetIntParaFromDB("imode00fix");
        if( imode00fix > 0 )
        {
            (*realPHI)(IF) = 2.0 * reynolds/( 1.0 + (*realnut_00)(IF) );
        }

        int ix = kstf;
        int iy = jstf;
        RDouble kxx = (*realKX)(ix);
        RDouble kyy = (*realKY)(iy);
        RDouble k2  = kxx * kxx + kyy * kyy;

        for( int iz = istf; iz <= iedf; ++iz )
        {
            (*complexULocal)(iz) = -(*complexNonlinearTerm)(iz, iy, ix, 1) - PHComplex(pressureGradient, 0.0) + (*complexViscousTerm)(iz, iy, ix, 1);
        }
        (*complexP)(IF, iy, ix) = cZero;

        PoissonSolver::PoissonSolverVelocity( nSolve, &(*coef_LHS)(istf+1, -2), &(*complexULocal)(istf), &(*lowerBound)(istf+1), &(*upperBound)(istf+1), cZero, cZero, &(*realPHI)(istf), 
            &(*coefVelocity)(istf, -2, iy, ix), &(*cR)(istf) );

        (*complexULocal)(IF) = (*complexU)(IF, iy, ix) + (*complexULocal)(IF);

        (*complexU)(IF, iy ,ix) = (*complexULocal)(IF);
        (*complexV)(IF, iy ,ix) = cZero;
        (*complexW)(IF, iy ,ix) = cZero;
    }

    (*realPHI)(IF) = 2.0 * reynolds/( 1.0 + (*realnut)(IF) );

    RDouble1D *detaDz = GridData -> realDetaDz;

    int nBoundary_Neumann = DiffBoundaryData -> nBoundary_Neumann;

    RDouble4D *coefMatrixPressureNeumann = PoissonSolverData -> coefMatrixPressureNeumann;
    RDouble4D *boundaryMatrixPressureNeumann = PoissonSolverData -> boundaryMatrixPressureNeumann;

    Complex3D *complexDeltaU1 = CorrectWallDivergenceData -> complexDeltaU1;
    Complex3D *complexDeltaU2 = CorrectWallDivergenceData -> complexDeltaU2;
    Complex3D *complexDeltaV1 = CorrectWallDivergenceData -> complexDeltaV1;
    Complex3D *complexDeltaV2 = CorrectWallDivergenceData -> complexDeltaV2;
    Complex3D *complexDeltaW1 = CorrectWallDivergenceData -> complexDeltaW1;
    Complex3D *complexDeltaW2 = CorrectWallDivergenceData -> complexDeltaW2;
    Complex3D *complexDeltaP1 = CorrectWallDivergenceData -> complexDeltaP1;
    Complex3D *complexDeltaP2 = CorrectWallDivergenceData -> complexDeltaP2;
    Complex2D *complexDivDeltaU1LowerBound = CorrectWallDivergenceData -> complexDivDeltaU1LowerBound;
    Complex2D *complexDivDeltaU1UpperBound = CorrectWallDivergenceData -> complexDivDeltaU1UpperBound;
    Complex2D *complexDivDeltaU2LowerBound = CorrectWallDivergenceData -> complexDivDeltaU2LowerBound;
    Complex2D *complexDivDeltaU2UpperBound = CorrectWallDivergenceData -> complexDivDeltaU2UpperBound;

    Complex1D *dVardZ = new Complex1D(IF, fortranArray);
    Complex1D *omegaZ = new Complex1D(IF, fortranArray);
    Complex1D *deltaPressure = new Complex1D(IF, fortranArray);
    Complex1D *divergenceOfUStar = new Complex1D(IF, fortranArray);        
    for( int ix = kstf; ix <= kedf; ++ix )
    {
        for( int iy = jstf; iy <= jedf; ++iy )
        {
            RDouble kxx = (*realKX)(ix);
            RDouble kyy = (*realKY)(iy);
            RDouble k2  = kxx * kxx + kyy * kyy;
            PHComplex ckx = (*complexKX)(ix);
            PHComplex cky = (*complexKY)(iy);

            if( k2 < TINY )
            {
                continue;
            }

            // ! Nyquist wavenumbers:
            RDouble realKNyquistY = GridData -> realKNyquistY;
            RDouble realKNyquistX = GridData -> realKNyquistX;
            if( (abs(kyy) == abs(realKNyquistY)) || (abs(kxx) == abs(realKNyquistX)) )
            {
                (*complexP)(IF, iy, ix) = cZero;
                (*complexU)(IF, iy ,ix) = cZero;
                (*complexV)(IF, iy, ix) = cZero;
                (*complexW)(IF, iy, ix) = cZero;
                continue;
            }

            (*complexPLocal)(IF) = (*complexP)(IF, iy ,ix);
            CompactDifferenceFirstDerivativeData -> GetDvarDz( istf, iedf, &(*detaDz)(istf), &(*complexPLocal)(istf), &(*dVardZ)(istf), &(*cR)(istf), &(*rC)(istf, -2) );

            for( int iz = istf; iz <= iedf; ++iz)
            {
                (*complexULocal)(iz) = -(*complexNonlinearTerm)(iz, iy, ix, 1) - ckx * (*complexPLocal)(iz) + (*complexViscousTerm)(iz, iy ,ix, 1);
                (*complexVLocal)(iz) = -(*complexNonlinearTerm)(iz, iy, ix, 2) - cky * (*complexPLocal)(iz) + (*complexViscousTerm)(iz, iy ,ix, 2);
                (*complexWLocal)(iz) = -(*complexNonlinearTerm)(iz, iy, ix, 3) - (*dVardZ)(iz) + (*complexViscousTerm)(iz, iy ,ix, 3);
            }

            PoissonSolver::PoissonSolverVelocity( nSolve, &(*coef_LHS)(istf+1, -2), &(*complexULocal)(istf), &(*lowerBound)(istf+1), &(*upperBound)(istf+1), cZero, cZero, &(*realPHI)(istf), 
                &(*coefVelocity)(istf, -2, iy, ix), &(*cR)(istf) );
            PoissonSolver::PoissonSolverVelocity( nSolve, &(*coef_LHS)(istf+1, -2), &(*complexVLocal)(istf), &(*lowerBound)(istf+1), &(*upperBound)(istf+1), cZero, cZero, &(*realPHI)(istf), 
                &(*coefVelocity)(istf, -2, iy, ix), &(*cR)(istf) );
            PoissonSolver::PoissonSolverVelocity( nSolve, &(*coef_LHS)(istf+1, -2), &(*complexWLocal)(istf), &(*lowerBound)(istf+1), &(*upperBound)(istf+1), cZero, cZero, &(*realPHI)(istf), 
                &(*coefVelocity)(istf, -2, iy, ix), &(*cR)(istf) );

            (*complexULocal)(IF) = (*complexU)(IF, iy, ix) + (*complexULocal)(IF);
            (*complexVLocal)(IF) = (*complexV)(IF, iy, ix) + (*complexVLocal)(IF);
            (*complexWLocal)(IF) = (*complexW)(IF, iy, ix) + (*complexWLocal)(IF);

            for( int iz = istf; iz <= iedf; ++iz )
            {
                (*divergenceOfUStar)(iz) = ckx * (*complexULocal)(iz) + cky * (*complexVLocal)(iz);
            }
            CompactDifferenceFirstDerivativeData -> GetDvarDz( istf, iedf, &(*detaDz)(istf), &(*complexWLocal)(istf), &(*dVardZ)(istf), &(*cR)(istf), &(*rC)(istf, -2) );
            (*divergenceOfUStar)(IF) = (*divergenceOfUStar)(IF) + (*dVardZ)(IF);
            (*deltaPressure)(IF) = (*divergenceOfUStar)(IF) / timeStep;

            PoissonSolverData -> PoissonSolverNeumannBC( nSolve, &(*coef_LHS)(istf+1, -2), &(*deltaPressure)(istf), &(*lowerBound)(istf+1), &(*upperBound)(istf+1), nBoundary_Neumann, cZero, cZero, 
                &(*cR)(istf), &(*coefMatrixPressureNeumann)(istf, -2, iy, ix), &(*boundaryMatrixPressureNeumann)(1, 1, iy, ix) );

            CompactDifferenceFirstDerivativeData -> GetDvarDz( istf, iedf, &(*detaDz)(istf), &(*deltaPressure)(istf), &(*dVardZ)(istf), &(*cR)(istf), &(*rC)(istf, -2) );
            for( int iz = istf; iz <= iedf; ++iz )
            {
                (*complexWLocal)(iz) = (*complexWLocal)(iz) - timeStep * (*dVardZ)(iz);
            }

            (*complexWLocal)(istf) = cZero;
            (*complexWLocal)(iedf) = cZero;

            CompactDifferenceFirstDerivativeData -> GetDvarDz( istf, iedf, &(*detaDz)(istf), &(*complexWLocal)(istf), &(*dVardZ)(istf), &(*cR)(istf), &(*rC)(istf, -2) );

            for( int iz = istf; iz <= iedf; ++iz )
            {
                (*omegaZ)(iz) = kyy * (*complexULocal)(iz) - kxx * (*complexVLocal)(iz);
                (*complexULocal)(iz) = ( cOneI * kxx * (*dVardZ)(iz) + kyy * (*omegaZ)(iz) ) / k2;
                (*complexVLocal)(iz) = ( cOneI * kyy * (*dVardZ)(iz) - kxx * (*omegaZ)(iz) ) / k2;
            }
            (*complexULocal)(istf) = cZero;
            (*complexULocal)(iedf) = cZero;
            (*complexVLocal)(istf) = cZero;
            (*complexVLocal)(iedf) = cZero;

            (*complexPLocal)(IF) = (*complexPLocal)(IF) + (*deltaPressure)(IF) - 1.0 / (*realPHI)(IF) * (*divergenceOfUStar)(IF);

            PHComplex dWdZLowerBound = (*dVardZ)(istf);
            PHComplex dWdZUpperBound = (*dVardZ)(iedf);
            PHComplex divOfDeltaU1LowerBound = (*complexDivDeltaU1LowerBound)(iy ,ix);
            PHComplex divOfDeltaU1UpperBound = (*complexDivDeltaU1UpperBound)(iy ,ix);
            PHComplex divOfDeltaU2LowerBound = (*complexDivDeltaU2LowerBound)(iy ,ix);
            PHComplex divOfDeltaU2UpperBound = (*complexDivDeltaU2UpperBound)(iy ,ix);

            PHComplex c1 = (dWdZUpperBound * divOfDeltaU2LowerBound - dWdZLowerBound * divOfDeltaU2UpperBound) 
                / (divOfDeltaU2UpperBound * divOfDeltaU1LowerBound - divOfDeltaU2LowerBound * divOfDeltaU1UpperBound);
            PHComplex c2 = (dWdZLowerBound * divOfDeltaU1UpperBound - dWdZUpperBound * divOfDeltaU1LowerBound) 
                / (divOfDeltaU2UpperBound * divOfDeltaU1LowerBound - divOfDeltaU2LowerBound * divOfDeltaU1UpperBound);

            for( int iz = istf; iz <= iedf; ++iz )
            {
                (*complexULocal)(iz) = (*complexULocal)(iz) + c1 * (*complexDeltaU1)(iz, iy ,ix) + c2 * (*complexDeltaU2)(iz, iy ,ix);
                (*complexVLocal)(iz) = (*complexVLocal)(iz) + c1 * (*complexDeltaV1)(iz, iy ,ix) + c2 * (*complexDeltaV2)(iz, iy ,ix);
                (*complexWLocal)(iz) = (*complexWLocal)(iz) + c1 * (*complexDeltaW1)(iz, iy ,ix) + c2 * (*complexDeltaW2)(iz, iy ,ix);
                (*complexPLocal)(iz) = (*complexPLocal)(iz) + c1 * (*complexDeltaP1)(iz, iy ,ix) + c2 * (*complexDeltaP2)(iz, iy ,ix);
            }

            (*complexU)(IF, iy, ix) = (*complexULocal)(IF);
            (*complexV)(IF, iy, ix) = (*complexVLocal)(IF);
            (*complexW)(IF, iy, ix) = (*complexWLocal)(IF);
            (*complexP)(IF, iy, ix) = (*complexPLocal)(IF);                            
        }
    }

    delete realPHI; realPHI = NULL;
    delete complexULocal; complexULocal = NULL;
    delete complexVLocal; complexVLocal = NULL;
    delete complexWLocal; complexWLocal = NULL;
    delete complexPLocal; complexPLocal = NULL;
    delete cR; cR = NULL;
    delete rC; rC = NULL;
    delete dVardZ; dVardZ = NULL;
    delete deltaPressure; deltaPressure = NULL;
    delete divergenceOfUStar; divergenceOfUStar = NULL;
    delete omegaZ; omegaZ = NULL;
}

void SpecDiffHybSolver::GetVelocityGradient()
{
    SpecDiffHybGrid *GridData = GetGrid();

    int ist = GridData -> iStartFourierIndex;
    int jst = GridData -> jStartFourierIndex;
    int kst = GridData -> kStartFourierIndex;
    int ied = GridData -> iEndFourierIndex;
    int jed = GridData -> jEndFourierIndex;
    int ked = GridData -> kEndFourierIndex;

    Complex1D *varTemp = new Complex1D( Range(ist, ied), fortranArray );
    Complex1D *dvarTemp = new Complex1D( Range(ist, ied), fortranArray );
    Complex1D *workSpace1 = new Complex1D( Range(ist, ied), fortranArray );
    RDouble1D *workSpace2 = new RDouble1D( Range(ist, ied), fortranArray );

    CompactDifferenceFirstDerivativeData -> GetGradientVector( &(*complexU)(ist, jst, kst), &(*complexVelocityGradient)(ist, jst, kst, 1), &(*complexVelocityGradient)(ist, jst, kst, 2), 
        &(*complexVelocityGradient)(ist, jst, kst, 3), varTemp, dvarTemp, workSpace1, workSpace2, GridData );
    CompactDifferenceFirstDerivativeData -> GetGradientVector( &(*complexV)(ist, jst, kst), &(*complexVelocityGradient)(ist, jst, kst, 4), &(*complexVelocityGradient)(ist, jst, kst, 5), 
        &(*complexVelocityGradient)(ist, jst, kst, 6), varTemp, dvarTemp, workSpace1, workSpace2, GridData );
    CompactDifferenceFirstDerivativeData -> GetGradientVector( &(*complexW)(ist, jst, kst), &(*complexVelocityGradient)(ist, jst, kst, 7), &(*complexVelocityGradient)(ist, jst, kst, 8), 
        &(*complexVelocityGradient)(ist, jst, kst, 9), varTemp, dvarTemp, workSpace1, workSpace2, GridData );

    delete varTemp; varTemp = NULL;
    delete dvarTemp; dvarTemp = NULL;
    delete workSpace1; workSpace1 = NULL;
    delete workSpace2; workSpace2 = NULL;
}

void SpecDiffHybSolver::GetNonlinearTerm()
{
    SpecDiffHybGrid *GridData = GetGrid();
    LES *LESData = GetLESData();

    Param_SpecDiffHybSolver *parameters = GetControlParameters();
    RDouble reynolds = parameters->GetRefReNumber();
    int viscousType = parameters->GetViscousType();
    string viscousName = parameters->GetViscousName();

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex  ;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex  ;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;
    Range IF(istf, iedf);
    Range JF(jstf, jedf);
    Range KF(kstf, kedf);

    int istp = GridData -> iStartPhysicalIndex;
    int iedp = GridData -> iEndPhysicalIndex  ;
    int jstp = GridData -> jStartPhysicalIndex;
    int jedp = GridData -> jEndPhysicalIndex  ;
    int kstp = GridData -> kStartPhysicalIndex;
    int kedp = GridData -> kEndPhysicalIndex  ;
    Range IP(istp, iedp);
    Range JP(jstp, jedp);
    Range KP(kstp, kedp);

    int nXYAll = GridData -> nXYAll;
    RDouble dnXYAll = 1.0 / static_cast<double>(nXYAll);

    (*complexVelocityGradient)(IF, JF, KF, 10) = (*complexU)(IF, JF, KF);
    (*complexVelocityGradient)(IF, JF, KF, 11) = (*complexV)(IF, JF, KF);
    (*complexVelocityGradient)(IF, JF, KF, 12) = (*complexW)(IF, JF, KF);
    int nXYZFourier = GridData -> nXYZFourier;
    int nXYZPhysical = GridData -> nXYZPhysical;

    unsigned char opc2r[] = "nff";
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

    unsigned char opr2c[] = "ffn";
    Cp3dfft_ftran_r2c_many(&(*realVelocitySquare)(istp, jstp, kstp, 1), nXYZPhysical, reinterpret_cast<double*>(&(*complexVelocitySquare)(istf, jstf, kstf, 1)), nXYZFourier, 12, opr2c);
    (*complexVelocitySquare)(IF, JF, KF, Range(1, 12)) = (*complexVelocitySquare)(IF, JF, KF, Range(1, 12)) * dnXYAll;

    Complex3D *complexDivergenceOfVelocitySquare = new Complex3D(IF, JF, KF, fortranArray);
    Complex1D *varTemp = new Complex1D(IF, fortranArray);
    Complex1D *dVarTemp = new Complex1D(IF, fortranArray);
    Complex3D *dVardZ = new Complex3D(IF, JF, KF, fortranArray); 
    Complex1D *workSpace1 = new Complex1D(IF, fortranArray);
    RDouble1D *workSpace2 = new RDouble1D(IF, fortranArray);
    Complex4D *convectiveTermTemp = new Complex4D(IF, JF, KF, Range(1, 3), fortranArray);

    // ! 0.5*[ vec{u}\dot\grad(u) + \nabla\dot(\vec{u}u) ]
    CompactDifferenceFirstDerivativeData -> GetDivergence( &(*complexVelocitySquare)(istf, jstf, kstf, 1), &(*complexVelocitySquare)(istf, jstf, kstf, 2), &(*complexVelocitySquare)(istf, jstf, kstf, 3),
        complexDivergenceOfVelocitySquare, varTemp, dVarTemp, dVardZ, workSpace1, workSpace2, GridData );
    (*convectiveTermTemp)(IF, JF, KF, 1) = 0.5 * ( (*complexVelocitySquare)(IF, JF, KF, 10) + (*complexDivergenceOfVelocitySquare)(IF, JF, KF) );

    CompactDifferenceFirstDerivativeData -> GetDivergence( &(*complexVelocitySquare)(istf, jstf, kstf, 4), &(*complexVelocitySquare)(istf, jstf, kstf, 5), &(*complexVelocitySquare)(istf, jstf, kstf, 6),
        complexDivergenceOfVelocitySquare, varTemp, dVarTemp, dVardZ, workSpace1, workSpace2, GridData );
    (*convectiveTermTemp)(IF, JF, KF, 2) = 0.5 * ( (*complexVelocitySquare)(IF, JF, KF, 11) + (*complexDivergenceOfVelocitySquare)(IF, JF, KF) );

    CompactDifferenceFirstDerivativeData -> GetDivergence( &(*complexVelocitySquare)(istf, jstf, kstf, 7), &(*complexVelocitySquare)(istf, jstf, kstf, 8), &(*complexVelocitySquare)(istf, jstf, kstf, 9),
        complexDivergenceOfVelocitySquare, varTemp, dVarTemp, dVardZ, workSpace1, workSpace2, GridData );
    (*convectiveTermTemp)(IF, JF, KF, 3) = 0.5 * ( (*complexVelocitySquare)(IF, JF, KF, 12) + (*complexDivergenceOfVelocitySquare)(IF, JF, KF) );

    if (viscousType <= 1)
    {
        (*realnut)(IF) = 0.0;
        (*realnut_00)(IF) = 0.0;
    }
    else
    {
        if(viscousName.substr(0, 11) == "Smagorinsky")
        {
            LESData -> Smagorinsky(GridData, realnut, realVelocityGradient);
            PoissonSolverData -> GetCoefMatrixVelocityDirichlet(GridData, DiffBoundaryData, CompactDifferenceSecondDerivativeData, realnut, realnut_00);
        }
        else if(viscousName.substr(0, 11) == "DynamicSmag")
        {
            LESData -> DynamicSmagorinsky(GridData, CompactDifferenceFirstDerivativeData, realnut, realnut3d, realVelocityGradient, complexU, complexV, complexW);
            PoissonSolverData -> GetCoefMatrixVelocityDirichlet(GridData, DiffBoundaryData, CompactDifferenceSecondDerivativeData, realnut, realnut_00);
        }
        else if(viscousName.substr(0, 5) == "Sigma")
        {
            LESData -> Sigma(GridData, realnut, realVelocityGradient);
            PoissonSolverData -> GetCoefMatrixVelocityDirichlet(GridData, DiffBoundaryData, CompactDifferenceSecondDerivativeData, realnut, realnut_00);
        }
        else
        {
            TK_Exit::UnexpectedVarValue("viscousName", viscousName);
        }
    }

    mode00fix();

    if(viscousType > 1)
    {
        int inutxyavg = GlobalDataBase::GetIntParaFromDB("inutxyavg");
        if(inutxyavg ==0 || inutxyavg ==1)
        {
            Complex4D *complexDiffTau = LESData -> complexDiffTau;
            (*convectiveTermTemp)(IF, JF, KF, Range(1, 3)) = (*convectiveTermTemp)(IF, JF, KF, Range(1, 3)) + (*complexDiffTau)(IF, JF, KF, Range(1, 3));
        }

        if(inutxyavg == 1 || inutxyavg ==2)
        {
            RDouble1D *detaDz = GridData -> realDetaDz;
            RDouble1D *dnutdz = new RDouble1D(IF, fortranArray);
            RDouble1D *workSpace3 = new RDouble1D(IF, fortranArray);
            CompactDifferenceFirstDerivativeData -> GetDvarDz(istf, iedf, &(*detaDz)(istf), &(*realnut)(istf), &(*dnutdz)(istf), &(*workSpace2)(istf), &(*workSpace3)(istf));

            for(int iz = istf; iz <= iedf; ++iz)
            {
                for(int iy = jstf; iy <= jedf; ++iy)
                {
                    for(int ix = kstf; ix <= kedf; ++ix)
                    {
                        (*convectiveTermTemp)(iz, iy, ix, 1) = (*convectiveTermTemp)(iz, iy, ix, 1) 
                            - (*dnutdz)(iz) / reynolds * ( (*complexVelocityGradient)(iz, iy, ix, 3) + (*complexVelocityGradient)(iz, iy, ix, 7) );
                        (*convectiveTermTemp)(iz, iy, ix, 2) = (*convectiveTermTemp)(iz, iy, ix, 2) 
                            - (*dnutdz)(iz) / reynolds * ( (*complexVelocityGradient)(iz, iy, ix, 6) + (*complexVelocityGradient)(iz, iy, ix, 8) );
                        (*convectiveTermTemp)(iz, iy, ix, 3) = (*convectiveTermTemp)(iz, iy, ix, 3) 
                            - (*dnutdz)(iz) / reynolds * ( 2.0 * (*complexVelocityGradient)(iz, iy, ix, 9) );
                    }
                }
            }

            int imode00fix = GlobalDataBase::GetIntParaFromDB("imode00fix");
            using namespace PHMPI;
            int myid = GetCurrentProcessorID();
            if((imode00fix > 0) && (myid ==0))
            {
                CompactDifferenceFirstDerivativeData -> GetDvarDz(istf, iedf, &(*detaDz)(istf), &(*realnut_00)(istf), &(*dnutdz)(istf), &(*workSpace2)(istf), &(*workSpace3)(istf));

                int iy = jstf;
                int ix = kstf;
                for(int iz = istf; iz <= iedf; ++iz)
                {
                    (*convectiveTermTemp)(iz, iy, ix, 1) = (*convectiveTermTemp)(iz, iy, ix, 1) 
                        - (*dnutdz)(iz) / reynolds * ( (*complexVelocityGradient)(iz, iy, ix, 3) + (*complexVelocityGradient)(iz, iy, ix, 7) );
                    (*convectiveTermTemp)(iz, iy, ix, 2) = (*convectiveTermTemp)(iz, iy, ix, 2) 
                        - (*dnutdz)(iz) / reynolds * ( (*complexVelocityGradient)(iz, iy, ix, 6) + (*complexVelocityGradient)(iz, iy, ix, 8) );
                    (*convectiveTermTemp)(iz, iy, ix, 3) = (*convectiveTermTemp)(iz, iy, ix, 3) 
                        - (*dnutdz)(iz) / reynolds * ( 2.0 * (*complexVelocityGradient)(iz, iy, ix, 9) );
                }
            }

            delete dnutdz; dnutdz = NULL;
            delete workSpace3; workSpace3 = NULL;
        }
    }

    (*complexNonlinearTerm)(IF, JF, KF, Range(1, 3)) = (*nonlinearTermCoef)(0) * (*convectiveTermTemp)(IF, JF, KF, Range(1, 3)) + (*nonlinearTermCoef)(-1) * (*convectiveTerm)(IF, JF, KF, Range(1, 3));
    (*convectiveTerm)(IF, JF, KF, Range(1, 3)) = (*convectiveTermTemp)(IF, JF, KF, Range(1, 3));

    delete complexDivergenceOfVelocitySquare; complexDivergenceOfVelocitySquare = NULL;
    delete varTemp; varTemp = NULL;
    delete dVarTemp; dVarTemp = NULL;
    delete dVardZ; dVardZ = NULL;
    delete workSpace1; workSpace1 = NULL;
    delete workSpace2; workSpace2 = NULL;
    delete convectiveTermTemp; convectiveTermTemp = NULL;
}

void SpecDiffHybSolver::mode00fix()
{
    SpecDiffHybGrid *GridData = GetGrid();

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    Range IF(istf, iedf);

    int imode00fix = GlobalDataBase::GetIntParaFromDB("imode00fix");
    if(imode00fix == 0)
    {
        (*realnut_00)(IF) = (*realnut)(IF);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("imode00fix", imode00fix);
    }
}

/*
void SpecDiffHybSolver::OutputFieldName()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    SpecDiffHybGrid *GridData = GetGrid();
    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();
    if(myid != procIDWithKx0Ky0)
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
*/

void SpecDiffHybSolver::GetViscousTerm()
{
    SpecDiffHybGrid *GridData = GetGrid();

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;
    Range IF(istf, iedf);
    Range JF(jstf, jedf);
    Range KF(kstf, kedf);

    Complex1D *c1tmp = new Complex1D(IF, fortranArray);
    Complex1D *c2tmp = new Complex1D(IF, fortranArray);
    Complex3D *c3tmp = new Complex3D(IF, JF, KF, fortranArray);
    Complex3D *cLap_u = new Complex3D(IF, JF, KF, fortranArray);
    Complex3D *cLap_v = new Complex3D(IF, JF, KF, fortranArray);
    Complex3D *cLap_w = new Complex3D(IF, JF, KF, fortranArray);

    CompactDifferenceSecondDerivativeData -> GetLaplacian(GridData, complexU, cLap_u, c1tmp, c2tmp, c3tmp);
    CompactDifferenceSecondDerivativeData -> GetLaplacian(GridData, complexV, cLap_v, c1tmp, c2tmp, c3tmp);
    CompactDifferenceSecondDerivativeData -> GetLaplacian(GridData, complexW, cLap_w, c1tmp, c2tmp, c3tmp);

    for(int iz = istf; iz <= iedf; ++iz)
    {
        (*complexViscousTerm)(iz, JF, KF, 1) = (1.0 + (*realnut)(iz)) * (*cLap_u)(iz, JF, KF);
        (*complexViscousTerm)(iz, JF, KF, 2) = (1.0 + (*realnut)(iz)) * (*cLap_v)(iz, JF, KF);
        (*complexViscousTerm)(iz, JF, KF, 3) = (1.0 + (*realnut)(iz)) * (*cLap_w)(iz, JF, KF);
    }

    int imode00fix = GlobalDataBase::GetIntParaFromDB("imode00fix");
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if((imode00fix > 0) && (myid ==0))
    {
        for(int iz = istf; iz <= iedf; ++iz)
        {
            (*complexViscousTerm)(iz, jstf, kstf, 1) = (1.0 + (*realnut)(iz)) * (*cLap_u)(iz, jstf, kstf);
            (*complexViscousTerm)(iz, jstf, kstf, 2) = (1.0 + (*realnut)(iz)) * (*cLap_v)(iz, jstf, kstf);
            (*complexViscousTerm)(iz, jstf, kstf, 3) = (1.0 + (*realnut)(iz)) * (*cLap_w)(iz, jstf, kstf);
        }
    }

    Param_SpecDiffHybSolver *parameters = GetControlParameters();
    RDouble reynolds = parameters->GetRefReNumber();

    (*complexViscousTerm)(IF, JF, KF, Range(1, 3)) = (*complexViscousTerm)(IF, JF, KF, Range(1, 3))/reynolds;

    delete c1tmp; c1tmp = NULL;
    delete c2tmp; c2tmp = NULL;
    delete c3tmp; c3tmp = NULL;
    delete cLap_u; cLap_u = NULL;
    delete cLap_v; cLap_v = NULL;
    delete cLap_w; cLap_w = NULL;
}

void SpecDiffHybSolver::VanishingDIVproj()
{
    SpecDiffHybGrid *GridData = GetGrid();

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;
    Range IF(istf, iedf);

    RDouble2D *rC = new RDouble2D(IF, Range(-2, 2), fortranArray);
    Complex1D *c1 = new Complex1D(IF, fortranArray);
    Complex1D *c2 = new Complex1D(IF, fortranArray);
    Complex1D *cu_star = new Complex1D(IF, fortranArray);
    Complex1D *cv_star = new Complex1D(IF, fortranArray);
    Complex1D *cw_star = new Complex1D(IF, fortranArray);
    Complex1D *cu_new = new Complex1D(IF, fortranArray);
    Complex1D *cv_new = new Complex1D(IF, fortranArray);
    Complex1D *cw_new = new Complex1D(IF, fortranArray);
    Complex1D *cdudz_star = new Complex1D(IF, fortranArray);
    Complex1D *cdvdz_star = new Complex1D(IF, fortranArray);
    Complex1D *cdwdz = new Complex1D(IF, fortranArray);
    Complex1D *comega_x = new Complex1D(IF, fortranArray);
    Complex1D *comega_y = new Complex1D(IF, fortranArray);
    Complex1D *cR = new Complex1D(IF, fortranArray);

    //OutputFieldRestartFile(0);

    const PHComplex cZero = PHComplex(0.0, 0.0);
    const PHComplex cOne = PHComplex(1.0, 0.0);
    int nSolve = iedf - istf;

    RDouble1D *realKX = GridData -> realKX;
    RDouble1D *realKY = GridData -> realKY;
    Complex1D *complexKX = GridData -> complexKX;
    Complex1D *complexKY = GridData -> complexKY;
    RDouble1D *detaDz = GridData -> realDetaDz;

    RDouble2D *coef_LHS = CompactDifferenceSecondDerivativeData -> lHSCoefMatrix;
    RDouble2D *coef_RHS = CompactDifferenceSecondDerivativeData -> rHSCoefMatrix;
    Int1D *lowerBound = CompactDifferenceSecondDerivativeData -> lowerBoundofRHSCoefMatrix;
    Int1D *upperBound = CompactDifferenceSecondDerivativeData -> upperBoundofRHSCoefMatrix;

    Complex2D *complexDivDeltaU1LowerBound = CorrectWallDivergenceData -> complexDivDeltaU1LowerBound;
    Complex2D *complexDivDeltaU1UpperBound = CorrectWallDivergenceData -> complexDivDeltaU1UpperBound;
    Complex2D *complexDivDeltaU2LowerBound = CorrectWallDivergenceData -> complexDivDeltaU2LowerBound;
    Complex2D *complexDivDeltaU2UpperBound = CorrectWallDivergenceData -> complexDivDeltaU2UpperBound;
    Complex3D *complexDeltaU1 = CorrectWallDivergenceData -> complexDeltaU1;
    Complex3D *complexDeltaV1 = CorrectWallDivergenceData -> complexDeltaV1;
    Complex3D *complexDeltaW1 = CorrectWallDivergenceData -> complexDeltaW1;
    Complex3D *complexDeltaU2 = CorrectWallDivergenceData -> complexDeltaU2;
    Complex3D *complexDeltaV2 = CorrectWallDivergenceData -> complexDeltaV2;
    Complex3D *complexDeltaW2 = CorrectWallDivergenceData -> complexDeltaW2;

    for(int ix = kstf; ix <= kedf; ++ix)
    {
        for(int iy = jstf; iy <= jedf; ++iy)
        {
            double rkx = (*realKX)(ix);
            double rky = (*realKY)(iy);
            PHComplex ckx = (*complexKX)(ix);
            PHComplex cky = (*complexKY)(iy);
            double k2  = rkx * rkx + rky * rky;

            if(k2 < TINY)
            {
                continue;
            }

            (*cu_star)(IF) = (*complexU)(IF, iy, ix);
            (*cv_star)(IF) = (*complexV)(IF, iy, ix);
            (*cw_star)(IF) = (*complexW)(IF, iy, ix);
            CompactDifferenceFirstDerivativeData -> GetDvarDz(istf, iedf, &(*detaDz)(istf), &(*cu_star)(istf), &(*cdudz_star)(istf), &(*cR)(istf), &(*rC)(istf, -2));
            CompactDifferenceFirstDerivativeData -> GetDvarDz(istf, iedf, &(*detaDz)(istf), &(*cv_star)(istf), &(*cdvdz_star)(istf), &(*cR)(istf), &(*rC)(istf, -2));
            for(int iz = istf; iz <= iedf; ++iz)
            {
                (*comega_x)(iz) = cky * (*cw_star)(iz) - (*cdvdz_star)(iz);
                (*comega_y)(iz) = -ckx * (*cw_star)(iz) + (*cdudz_star)(iz);
                (*cw_new)(iz) = cky * (*comega_x)(iz) - ckx * (*comega_y)(iz);
            }

            PoissonSolver::PoissonSolverDirichletBC(nSolve, &(*coef_LHS)(istf+1, -2), &(*coef_RHS)(istf+1, -2), &(*cw_new)(istf), k2, &(*lowerBound)(istf+1), &(*upperBound)(istf+1),
                cZero, cZero, &(*rC)(istf, -2), &(*cR)(istf));

            (*cw_new)(istf) = cZero;
            (*cw_new)(iedf) = cZero;
            CorrectWallDivergenceData -> rotProjection();

            (*cu_new)(istf) = cZero;
            (*cv_new)(istf) = cZero;
            (*cu_new)(iedf) = cZero;
            (*cv_new)(iedf) = cZero;

            PHComplex cmDIVp = -(*cdwdz)(iedf);
            PHComplex cmDIVm = -(*cdwdz)(istf);
            PHComplex cdd = (*complexDivDeltaU1UpperBound)(iy, ix) * (*complexDivDeltaU2LowerBound)(iy, ix) - (*complexDivDeltaU1LowerBound)(iy, ix) * (*complexDivDeltaU2UpperBound)(iy, ix);
            PHComplex cb1 = ( cmDIVp * (*complexDivDeltaU2LowerBound)(iy, ix) - cmDIVm * (*complexDivDeltaU2UpperBound)(iy, ix)) / cdd;
            PHComplex cb2 = (-cmDIVp * (*complexDivDeltaU1LowerBound)(iy, ix) + cmDIVm * (*complexDivDeltaU1UpperBound)(iy, ix)) / cdd;

            for(int iz = istf; iz <= iedf; ++iz)
            {
                (*cu_new)(iz) = (*cu_new)(iz) + cb1 * (*complexDeltaU1)(iz, iy, ix) + cb2 * (*complexDeltaU2)(iz, iy, ix);
                (*cv_new)(iz) = (*cv_new)(iz) + cb1 * (*complexDeltaV1)(iz, iy, ix) + cb2 * (*complexDeltaV2)(iz, iy, ix);
                (*cw_new)(iz) = (*cw_new)(iz) + cb1 * (*complexDeltaW1)(iz, iy, ix) + cb2 * (*complexDeltaW2)(iz, iy, ix);
            }
        }
    }

    //OutputFieldRestartFile(10);

    delete rC; rC = NULL;
    delete c1; c1 = NULL;
    delete c2; c2 = NULL;
    delete cu_star; cu_star = NULL;
    delete cv_star; cv_star = NULL;
    delete cw_star; cw_star = NULL;
    delete cu_new; cu_new = NULL;
    delete cv_new; cv_new = NULL;
    delete cw_new; cw_new = NULL;
    delete cdudz_star; cdudz_star = NULL;
    delete cdvdz_star; cdvdz_star = NULL;
    delete cdwdz; cdwdz = NULL;
    delete comega_x; comega_x = NULL;
    delete comega_y; comega_y = NULL;
    delete cR; cR = NULL;
}

void SpecDiffHybSolver::RescaleVelocityField(const int nstep)
{
    SpecDiffHybGrid *GridData = GetGrid();
    RescaleField *RescaleFieldData = GetRescaleFieldData();

    int nRescaleFluctuation = GlobalDataBase::GetIntParaFromDB("nRescaleFluctuation");
    int nRescaleUAverage = GlobalDataBase::GetIntParaFromDB("nRescaleUAverage");
    int nFixUAverage = GlobalDataBase::GetIntParaFromDB("nFixUAverage");

    if( (nstep < nRescaleFluctuation) || (nstep < nRescaleUAverage) )
    {
        RescaleFieldData -> GetUVWMagnitude(GridData, complexU, complexV, complexW);
    }

    if(nstep < nRescaleFluctuation)
    {
        RescaleFieldData -> RescaleFluctuation(GridData, complexU, complexV, complexW);
    }

    if(nstep < nRescaleUAverage)
    {
        RescaleFieldData -> RescaleUAverage(GridData, complexU);
    }

    if( (nFixUAverage < 0) || (nstep < nFixUAverage) )
    {
        RescaleFieldData -> FixUAverage(GridData, complexU);
    }
}

void SpecDiffHybSolver::GetStatistics()
{
    Param_SpecDiffHybSolver *parameters = GetControlParameters();
    RDouble reynolds = parameters->GetRefReNumber();
    int statisticMethod = parameters->GetStatisticMethod();
    int viscousType = parameters->GetViscousType();
    RDouble timeStep = parameters->GetTimeStep();

    int nStatisticalStep = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");
    ++ nStatisticalStep;
    GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);

    SpecDiffHybGrid *GridData = GetGrid();
    Statistics *StatisticData = GetStatisticData();

    StatisticData -> GetStatisticsLocal(GridData, complexU, complexV, complexW, complexVelocityGradient);

    int istf = GridData -> iStartFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;
    int jStartFourierIndexGlobal = GridData -> jStartFourierIndexGlobal;
    int kEndFourierIndexGlobal = GridData -> kEndFourierIndexGlobal;
    int jEndFourierIndexGlobal = GridData -> jEndFourierIndexGlobal;

    Range I(istf, iedf);
    Range J(jstf, jedf);
    Range K(kstf, kedf);
    Range XN(kStartFourierIndexGlobal, kEndFourierIndexGlobal);
    Range YN(jStartFourierIndexGlobal, jEndFourierIndexGlobal);
    Range N19(1, 19);
    Range N9(1, 9);
    Range N3(1, 3);

    RDouble2D *statisticLocal = StatisticData -> statisticLocal;
    RDouble2D *statisticAll = StatisticData -> statisticAll;

    (*statisticAll)(I, N19) = (*statisticLocal)(I, N19);

    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();

#ifdef PH_PARALLEL
    int nstat = 19 * (iedf - istf + 1);
    MPI_Reduce(&(*statisticLocal)(istf, 1), &(*statisticAll)(istf, 1), nstat, MPI_DOUBLE, MPI_SUM, procIDWithKx0Ky0, MPI_COMM_WORLD);
#endif

    RDouble3D *energyOfUUKxAll = StatisticData -> energyOfUUKxAll;
    RDouble3D *energyOfUUKyAll = StatisticData -> energyOfUUKyAll;
    RDouble3D *energyOfUUKxLocal = StatisticData -> energyOfUUKxLocal;
    RDouble3D *energyOfUUKyLocal = StatisticData -> energyOfUUKyLocal;
    Int2D *startAndEndFourierIndexOfAll = GridData -> startAndEndFourierIndexOfAll;

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int nTProcessor = GetNumberOfProcessor();   

#ifdef PH_PARALLEL
    if(myid == 0)
    {
#endif
        (*energyOfUUKxAll)(I, XN, N3) = 0.0;
        (*energyOfUUKxAll)(I, K, N3) = (*energyOfUUKxLocal)(I, K, N3);

#ifdef PH_PARALLEL
        RDouble3D *energyOfUUKxTmp = new RDouble3D(I, XN, N3, fortranArray);
        for(int ip = 1; ip <= (nTProcessor - 1); ++ip)
        {
            int kstf_ip = (*startAndEndFourierIndexOfAll)(3, ip);
            int kedf_ip = (*startAndEndFourierIndexOfAll)(6, ip);
            int nenergy = (kedf_ip - kstf_ip + 1) * (iedf - istf + 1);
            (*energyOfUUKxTmp)(I, XN, N3) = 0.0;

            for(int in = 1; in <= 3; ++in)
            {
                MPI_Status status;
                MPI_Recv(&(*energyOfUUKxTmp)(istf, kstf_ip, in), nenergy, MPI_DOUBLE, ip, ip, MPI_COMM_WORLD, &status);
                (*energyOfUUKxAll)(I, Range(kstf_ip, kedf_ip), in) = (*energyOfUUKxAll)(I, Range(kstf_ip, kedf_ip), in) + (*energyOfUUKxTmp)(I, Range(kstf_ip, kedf_ip), in);
            }   
        }

        delete energyOfUUKxTmp; energyOfUUKxTmp = NULL;
    }
    else
    {
        int nenergy = (kedf - kstf + 1) * (iedf - istf + 1);
        for(int in = 1; in <= 3; ++in)
        {
            MPI_Send(&(*energyOfUUKxLocal)(istf, kstf, in), nenergy,  MPI_DOUBLE, 0, myid, MPI_COMM_WORLD);
        }
    }
#endif

#ifdef PH_PARALLEL
    if(myid == 0)
    {
#endif
        (*energyOfUUKyAll)(I, YN, N3) = 0.0;
        (*energyOfUUKyAll)(I, J, N3) = (*energyOfUUKyLocal)(I, J, N3);

#ifdef PH_PARALLEL
        RDouble3D *energyOfUUKyTmp = new RDouble3D(I, YN, N3, fortranArray);
        for(int ip = 1; ip <= (nTProcessor - 1); ++ip)
        {
            int jstf_ip = (*startAndEndFourierIndexOfAll)(2, ip);
            int jedf_ip = (*startAndEndFourierIndexOfAll)(5, ip);
            int nenergy = (jedf_ip - jstf_ip + 1) * (iedf - istf + 1);
            (*energyOfUUKyTmp)(I, YN, N3) = 0.0;

            for(int in = 1; in <= 3; ++in)
            {
                MPI_Status status;
                MPI_Recv(&(*energyOfUUKyTmp)(istf, jstf_ip, in), nenergy, MPI_DOUBLE, ip, ip, MPI_COMM_WORLD, &status);
                (*energyOfUUKyAll)(I, Range(jstf_ip, jedf_ip), in) = (*energyOfUUKyAll)(I, Range(jstf_ip, jedf_ip), in) + (*energyOfUUKyTmp)(I, Range(jstf_ip, jedf_ip), in);
            }      
        }

        delete energyOfUUKyTmp; energyOfUUKyTmp = NULL;
    }
    else
    {
        int nenergy = (jedf - jstf + 1) * (iedf - istf + 1);
        for(int in = 1; in <= 3; ++in)
        {
            MPI_Send(&(*energyOfUUKyLocal)(istf, jstf, in), nenergy,  MPI_DOUBLE, 0, myid, MPI_COMM_WORLD);
        }   
    }
#endif

    if(myid != procIDWithKx0Ky0)
    {
        return;
    }

    int innz = iedf - istf + 1;
    int iz_half = StatisticData -> iz_half;
    Range IHalf(istf, iz_half);
    //RDouble2D *epsilonXYZTime = new RDouble2D(IHalf, N9, fortranArray);
    RDouble2D *statisticTime = new RDouble2D(IHalf, N19, fortranArray);
    RDouble3D *energyOfUUKxT = new RDouble3D(IHalf, XN, N3, fortranArray);
    RDouble3D *energyOfUUKyTT = new RDouble3D(IHalf, YN, N3, fortranArray);
    //(*epsilonXYZTime)(IHalf, N9) = 0.0;
    (*statisticTime)(IHalf, N19) = 0.0;
    (*energyOfUUKxT)(IHalf, XN, N3) = 0.0;
    (*energyOfUUKyTT)(IHalf, YN, N3) = 0.0;

    for(int iz = istf; iz <= iz_half; ++iz)
    {
        int izup = innz - iz + istf;

        (*statisticTime)(iz, N19) = 0.5 * ((*statisticAll)(iz, N19) + (*statisticAll)(izup, N19));

        (*energyOfUUKxT)(iz, XN, N3) = 0.5 * ((*energyOfUUKxAll)(iz, XN, N3) + (*energyOfUUKxAll)(izup, XN, N3));
        (*energyOfUUKyTT)(iz, YN, N3) = 0.5 * ((*energyOfUUKyAll)(iz, YN, N3) + (*energyOfUUKyAll)(izup, YN, N3));
    }

    int nY = GridData -> nY;
    Range YNHalf(jStartFourierIndexGlobal, nY/2+1);
    RDouble3D *energyOfUUKyT = new RDouble3D(IHalf, YNHalf, N3, fortranArray);
    (*energyOfUUKyT)(IHalf, YNHalf, N3) = 0.0;
    (*energyOfUUKyT)(IHalf, jStartFourierIndexGlobal, N3) = (*energyOfUUKyTT)(IHalf, jStartFourierIndexGlobal, N3);
    (*energyOfUUKyT)(IHalf, nY/2+1, N3) = (*energyOfUUKyTT)(IHalf, nY/2+1, N3);
    for(int iy = jStartFourierIndexGlobal+1; iy <= nY/2; ++iy)
    {
        (*energyOfUUKyT)(IHalf, iy, N3) = (*energyOfUUKyTT)(IHalf, iy, N3) + (*energyOfUUKyTT)(IHalf, nY - iy + 2, N3);
    }

    RDouble tavg = GlobalDataBase::GetDoubleParaFromDB("tavg");
    bool ifInitStatistic = StatisticData -> GetIfInitStatistic();
    RDouble2D *statistic = StatisticData -> statistic;

    RDouble3D *energyOfUUKx = StatisticData -> energyOfUUKx;
    RDouble3D *energyOfUUKy = StatisticData -> energyOfUUKy;

    if(statisticMethod == 0)
    {
        RDouble dtdtavg = timeStep / tavg;
        if(ifInitStatistic)
        {
            dtdtavg = 1.0;
        }

        (*statistic)(IHalf, N19) = (1.0 - dtdtavg) * (*statistic)(IHalf, N19) + dtdtavg * (*statisticTime)(IHalf, N19);

        (*energyOfUUKx)(IHalf, XN, N3) = (1.0 - dtdtavg) * (*energyOfUUKx)(IHalf, XN, N3) + dtdtavg * (*energyOfUUKxT)(IHalf, XN, N3);
        (*energyOfUUKy)(IHalf, YNHalf, N3) = (1.0 - dtdtavg) * (*energyOfUUKy)(IHalf, YNHalf, N3) + dtdtavg * (*energyOfUUKyT)(IHalf, YNHalf, N3);

        RDouble dUdZLowerBound = real((*complexVelocityGradient)(istf, jstf, kstf, 3));
        RDouble dUdZUpperBound = real((*complexVelocityGradient)(iedf, jstf, kstf, 3));
        RDouble tauWAverage = StatisticData -> GetTauWAverage();
        RDouble tauW = 0.5 * ( abs(dUdZLowerBound) +abs(dUdZUpperBound) ) / reynolds;
        tauWAverage = (1.0 - dtdtavg) * tauWAverage + dtdtavg * tauW;
        StatisticData->SetTauWAverage(tauWAverage);

        RDouble1D *nutAverage = StatisticData -> nutAverage;

        //int inutavg = StatisticData -> inutavg;
        if(viscousType > 1)
        {
            //if(inutavg == 0 || ifInitStatistic)
            //{
            //    (*nutAverage)(I) = (*realnut)(I);
            //}
            //else
            //{
                (*nutAverage)(I) = (1.0 - dtdtavg) * (*nutAverage)(I) + dtdtavg * (*realnut)(I);
            //}
            //inutavg = 1;
        }

        //if(ifInitStatistic)
        //{
        //    ifInitStatistic = false;
        //    StatisticData->SetIfInitStatistic(ifInitStatistic);
        //}
    }
    else if(statisticMethod == 1)
    {           
        //int nStatisticalStep = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");

        RDouble c1 = 1.0 / (nStatisticalStep * 1.0);
        RDouble c2 = 1.0 - c1;

        (*statistic)(IHalf, N19) = c2 * (*statistic)(IHalf, N19) + c1 * (*statisticTime)(IHalf, N19);

        (*energyOfUUKx)(IHalf, XN, N3) = c2 * (*energyOfUUKx)(IHalf, XN, N3) + c1 * (*energyOfUUKxT)(IHalf, XN, N3);
        (*energyOfUUKy)(IHalf, YNHalf, N3) = c2 * (*energyOfUUKy)(IHalf, YNHalf, N3) + c1 * (*energyOfUUKyT)(IHalf, YNHalf, N3);

        RDouble dUdZLowerBound = real((*complexVelocityGradient)(istf, jstf, kstf, 3));
        RDouble dUdZUpperBound = real((*complexVelocityGradient)(iedf, jstf, kstf, 3));

        RDouble tauWAverage = StatisticData -> GetTauWAverage();

        RDouble tauW = 0.5 * ( abs(dUdZLowerBound) +abs(dUdZUpperBound) ) / reynolds;

        tauWAverage = c2 * tauWAverage + c1 * tauW;

        StatisticData->SetTauWAverage(tauWAverage);

        RDouble1D *nutAverage = StatisticData -> nutAverage;

        //int inutavg = StatisticData -> inutavg;
        if(viscousType > 1)
        {
            //if(inutavg == 0 || ifInitStatistic)
            //{
            //    (*nutAverage)(I) = (*realnut)(I);
            //}
            //else
            //{
                (*nutAverage)(I) = c2 * (*nutAverage)(I) + c1 * (*realnut)(I);
            //}
            //inutavg = 1;
        }

        //if(ifInitStatistic)
        //{
        //    ifInitStatistic = false;
        //    StatisticData->SetIfInitStatistic(ifInitStatistic);
        //}
    }
    else
    {
        TK_Exit::UnexpectedVarValue("statisticMethod", statisticMethod);
    }

    delete statisticTime; statisticTime = NULL;
    delete energyOfUUKxT; energyOfUUKxT = NULL;
    delete energyOfUUKyTT; energyOfUUKyTT = NULL;
    delete energyOfUUKyT; energyOfUUKyT = NULL;
}

/*
void SpecDiffHybSolver::OutputResidual()
{
    //SpecDiffHybGrid *GridData = GetGrid();

    //PrepareResidual(res);

    //string outputdir = "./results";
    //outputdir = GlobalDataBase::GetStrParaFromDB("outputdir");
    //string residualFile = outputdir + "/monitor.dat";
    //string nutFile = outputdir + "/nut.dat";

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex  ;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex  ;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex; 
    Range I(istf, iedf);
    Range J(jstf, jedf);
    Range K(kstf, kedf);
    (*complexDivergenceU)(I, J, K) = (*complexVelocityGradient)(I, J, K, 1) + (*complexVelocityGradient)(I, J, K, 5) + (*complexVelocityGradient)(I, J, K, 9);

    RDouble1D *rZ = GridData -> realZ;
    RDouble ub = 0.0;
    RDouble tauW = 0.0;
    int procIDWithKx0Ky0 = GridData -> procIDWithKx0Ky0;
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if(myid == procIDWithKx0Ky0)
    {
        RDouble1D *uAverage = new RDouble1D(I, fortranArray);
        for(int iz = istf; iz <= iedf; ++iz)
        {
            (*uAverage)(iz) = real((*complexU)(iz, jstf, kstf));
        }

        RDouble dUdZLowerBound = real((*complexVelocityGradient)(istf, jstf, kstf, 3));
        RDouble dUdZUpperBound = real((*complexVelocityGradient)(iedf, jstf, kstf, 3));

        tauW = 0.5 * ( abs(dUdZLowerBound) +abs(dUdZUpperBound) ) / reynolds;

        ub = 0.0;
        for(int iz = (istf+1); iz <= (iedf-1); ++iz)
        {
            ub = ub + 0.5 * ((*rZ)(iz+1) - (*rZ)(iz-1)) * (*uAverage)(iz);
        }
        ub = ub / ((*rZ)(iedf) - (*rZ)(istf));

        delete uAverage; uAverage = NULL;
    }

    RDouble1D *kx = GridData -> realKX;
    RDouble1D *ky = GridData -> realKY;
    RDouble ek = 0.0;
    for(int ix = kstf; ix <= kedf; ++ix)
    {
        for(int iy = jstf; iy <= jedf; ++iy)
        {
            RDouble k2  = (*kx)(ix) * (*kx)(ix) + (*ky)(iy) * (*ky)(iy);
            if( k2 < TINY)
            {
                continue;
            }
            for(int iz = (istf+1); iz <= (iedf-1); ++iz)
            {
                ek = ek + 2.0 * real( (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix))
                    + (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix))
                    + (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix))
                    );
            }
        }// !for iy
    }// !for ix
    ek = ek / (iedf - istf - 1);

    RDouble divavg = 0.0;
    for(int iz = istf; iz <= iedf; ++iz)
    {
        RDouble sumDivU = 0.0;
        for(int ix = kstf; ix <= kedf; ++ix)
        {
            for(int iy = jstf; iy <= jedf; ++iy)
            {
                sumDivU = sumDivU + real( (*complexDivergenceU)(iz, iy, ix) * conj((*complexDivergenceU)(iz, iy, ix)) );
            }
        }
        divavg = divavg + sqrt(2.0 * sumDivU);
    }
    divavg = divavg / (iedf - istf + 1);

#ifdef PH_PARALLEL
    RDouble1D *wmpiloc = new RDouble1D(Range(1, 2), fortranArray);
    RDouble1D *wmpi = new RDouble1D(Range(1, 2), fortranArray);
    (*wmpiloc)(1) = ek;
    (*wmpiloc)(2) = divavg;
    MPI_Reduce(&(*wmpiloc)(1), &(*wmpi)(1), 2, MPI_DOUBLE, MPI_SUM, procIDWithKx0Ky0, MPI_COMM_WORLD);
    ek = (*wmpi)(1);
    divavg = (*wmpi)(2);
    delete wmpiloc; wmpiloc = NULL;
    delete wmpi; wmpi = NULL;
#endif


    //using namespace PHMPI;
    //int myid = GetCurrentProcessorID();
    //int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();

    ActionKey *actkeyDumpResidual = new ActionKey();
    FillActionKey(actkeyDumpResidual, DUMP_RESIDUAL, 0);
    DumpRes(actkeyDumpResidual);
    delete actkeyDumpResidual;

        /*
        fstream file;
        file.open(screenOutFile.c_str(), ios_base::out|ios_base::app);
        file << "\n";
        file << "-----------------------------------------------------------------------------------------------------\n";
        file << "nstep = " << nstep << "\n";
        file << "DIVavg = " << divavg << "\n";
        file.close();
        file.clear();

        RDouble tauWAverage = StatisticsData -> tauWAverage;

        int nstart = GlobalDataBase::GetIntParaFromDB("nstart");
        file.open(residualFile.c_str(), ios_base::in);
        if ((!file) || (!nstart))
        {
            file.open(residualFile.c_str(), ios_base::out|ios_base::trunc);
            file << "variables= \n";
            file << "       nstep       time        rDIVavg          rUb            rEk            Tauw            rTauw_avg \n";
            file << setw(20) << nstep;
            file << setiosflags(ios::right);
            file << setprecision(10);
            file << setiosflags(ios::scientific);
            file << setiosflags(ios::showpoint);
            file << setw(20) << static_cast<double>(nstep) * timeStep 
                << setw(20) << divavg 
                << setw(20) << ub 
                << setw(20) << ek 
                << setw(20) << tauW 
                << setw(20) << tauWAverage 
                << "\n";
            file.close();
            file.clear();
        }
        else
        {
            file.close();
            file.clear();
            file.open(residualFile.c_str(), ios_base::out|ios_base::app);
            file << setw(20) << nstep;
            file << setiosflags(ios::right);
            file << setprecision(10);
            file << setiosflags(ios::scientific);
            file << setiosflags(ios::showpoint);
            file << setw(20) << static_cast<double>(nstep) * timeStep 
                << setw(20) << divavg 
                << setw(20) << ub 
                << setw(20) << ek 
                << setw(20) << tauW 
                << setw(20) << tauWAverage 
                << "\n";
            file.close();
            file.clear();
        }

        int viscousType = parameters->GetViscousType();
        if(viscousType > 1)
        {
            RDouble1D *nutAverage = StatisticsData -> nutAverage;
            file.open(nutFile.c_str(), ios_base::out);
            if ( !file )
            {
                TK_Exit::ExceptionExit("could not open nutFile\n");
            }
            file << "variables= z, nut, rnut_avg, rnut_00 \n";
            for(int iz = istf; iz <= iedf; ++iz)
            {
                file << setiosflags(ios::right);
                file << setprecision(10);
                file << setiosflags(ios::scientific);
                file << setiosflags(ios::showpoint);
                file << setw(20) << (*rZ)(iz) << setw(20) << (*realnut)(iz) << setw(20) << (*nutAverage)(iz) << setw(20) << (*realnut_00)(iz) << "\n";
            }
        }
        
}
*/

void SpecDiffHybSolver::DumpRes(ActionKey *actkey)
{
    SpecDiffHybGrid *GridData = GetGrid();

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
        oss << "DIVavg" << endl;
        oss << "U<sub>b</sub>" << endl;
        oss << "E<sub>k</sub>" << endl;
        oss << "<greek>t</greek><sub>w</sub>" << endl;
        oss << "<greek>t</greek><sub>w_avg</sub>" << endl;
    }

    oss << formatRes;

    WriteASCIIFile(file, oss.str());

    ParallelCloseFile(actkey);
}

void SpecDiffHybSolver::PrepareFormatedResidual(string &res)
{
    SpecDiffHybGrid *GridData = GetGrid();
    Statistics *StatisticData = GetStatisticData();

    Param_SpecDiffHybSolver *parameters = GetControlParameters();
    RDouble reynolds = parameters->GetRefReNumber();
    RDouble timeStep = parameters->GetTimeStep();

    int nstep = GlobalDataBase::GetIntParaFromDB("nstep");

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex  ;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex  ;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex; 
    Range I(istf, iedf);
    Range J(jstf, jedf);
    Range K(kstf, kedf);
    (*complexDivergenceU)(I, J, K) = (*complexVelocityGradient)(I, J, K, 1) + (*complexVelocityGradient)(I, J, K, 5) + (*complexVelocityGradient)(I, J, K, 9);

    RDouble1D *rZ = GridData -> realZ;
    RDouble ub = 0.0;
    RDouble tauW = 0.0;
    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if(myid == procIDWithKx0Ky0)
    {
        RDouble1D *uAverage = new RDouble1D(I, fortranArray);
        for(int iz = istf; iz <= iedf; ++iz)
        {
            (*uAverage)(iz) = real((*complexU)(iz, jstf, kstf));
        }

        RDouble dUdZLowerBound = real((*complexVelocityGradient)(istf, jstf, kstf, 3));
        RDouble dUdZUpperBound = real((*complexVelocityGradient)(iedf, jstf, kstf, 3));

        tauW = 0.5 * ( abs(dUdZLowerBound) +abs(dUdZUpperBound) ) / reynolds;

        ub = 0.0;
        for(int iz = (istf+1); iz <= (iedf-1); ++iz)
        {
            ub = ub + 0.5 * ((*rZ)(iz+1) - (*rZ)(iz-1)) * (*uAverage)(iz);
        }
        ub = ub / ((*rZ)(iedf) - (*rZ)(istf));

        delete uAverage; uAverage = NULL;
    }

    RDouble1D *kx = GridData -> realKX;
    RDouble1D *ky = GridData -> realKY;
    RDouble ek = 0.0;
    for(int ix = kstf; ix <= kedf; ++ix)
    {
        for(int iy = jstf; iy <= jedf; ++iy)
        {
            RDouble k2  = (*kx)(ix) * (*kx)(ix) + (*ky)(iy) * (*ky)(iy);
            if( k2 < TINY)
            {
                continue;
            }
            for(int iz = (istf+1); iz <= (iedf-1); ++iz)
            {
                ek = ek + 2.0 * real( (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix))
                    + (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix))
                    + (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix))
                    );
            }
        }// !for iy
    }// !for ix
    ek = ek / (iedf - istf - 1);

    RDouble divavg = 0.0;
    for(int iz = istf; iz <= iedf; ++iz)
    {
        RDouble sumDivU = 0.0;
        for(int ix = kstf; ix <= kedf; ++ix)
        {
            for(int iy = jstf; iy <= jedf; ++iy)
            {
                sumDivU = sumDivU + real( (*complexDivergenceU)(iz, iy, ix) * conj((*complexDivergenceU)(iz, iy, ix)) );
            }
        }
        divavg = divavg + sqrt(2.0 * sumDivU);
    }
    divavg = divavg / (iedf - istf + 1);

#ifdef PH_PARALLEL
    RDouble1D *wmpiloc = new RDouble1D(Range(1, 2), fortranArray);
    RDouble1D *wmpi = new RDouble1D(Range(1, 2), fortranArray);
    (*wmpiloc)(1) = ek;
    (*wmpiloc)(2) = divavg;
    MPI_Reduce(&(*wmpiloc)(1), &(*wmpi)(1), 2, MPI_DOUBLE, MPI_SUM, procIDWithKx0Ky0, MPI_COMM_WORLD);
    ek = (*wmpi)(1);
    divavg = (*wmpi)(2);
    delete wmpiloc; wmpiloc = NULL;
    delete wmpi; wmpi = NULL;
#endif

    if(myid == procIDWithKx0Ky0)
    {
        RDouble tauWAverage = StatisticData -> GetTauWAverage();

        ostringstream oss;

        oss << setiosflags(ios::left);
        oss << setprecision(10);
        oss << setiosflags(ios::scientific);
        oss << setiosflags(ios::showpoint);

        oss << nstep << "    ";
        oss << static_cast<double>(nstep) * timeStep << "    ";
        oss << divavg << "    ";
        oss << ub << "    ";
        oss << ek << "    ";
        oss << tauW << "    ";
        oss << tauWAverage << endl;

        res = oss.str();
    }
}

void SpecDiffHybSolver::OutputRestartData()
{
    OutputFieldRestartFile();

    OutputOldTimeData();

    StatisticData -> OutputStatisticRestart(GridData);
}

void SpecDiffHybSolver::OutputStatisticVisual()
{
    StatisticData -> OutputStatisticVisual(GridData);
}

void SpecDiffHybSolver::OutputOldTimeData()
{
    Param_SpecDiffHybSolver *parameters = GetControlParameters();
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

    int iVortProjection = GlobalDataBase::GetIntParaFromDB("iVortProjection");
    if(iVortProjection != 0)
    {
        TK_Exit::UnexpectedVarValue( "iVortProjection", iVortProjection );
    }

    WriteTimeDependentInformation();
}

void SpecDiffHybSolver::WriteConvectiveTerm(const string nonlinearTermFileName)
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

void SpecDiffHybSolver::WriteTimeDependentInformation()
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

void SpecDiffHybSolver::DeAllocateGlobalVariables()
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

    delete realnut; realnut = NULL;
    delete realnut_00; realnut_00 = NULL;

    delete GridData; GridData = NULL;

    //delete RescaleFieldData; RescaleFieldData = NULL;

    //delete StatisticData; StatisticData = NULL;

    delete DiffBoundaryData; DiffBoundaryData = NULL;

    delete CompactDifferenceFirstDerivativeData; CompactDifferenceFirstDerivativeData = NULL;

    delete CompactDifferenceSecondDerivativeData; CompactDifferenceSecondDerivativeData = NULL;

    delete PoissonSolverData; PoissonSolverData = NULL;

    delete CorrectWallDivergenceData; CorrectWallDivergenceData = NULL;

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    if(viscousType > 1)
    {
        delete LESData; LESData = NULL;
    }
}

LIB_EXPORT void SpecDiffHybSolver::InitControlParameters()
{
    FreeControlParameters();
    controlParameters = new Param_SpecDiffHybSolver();
    controlParameters->Init();

    WriteLogFile("Control parameters initialized!\n");
}

LIB_EXPORT Param_SpecDiffHybSolver * SpecDiffHybSolver::GetControlParameters() const
{
    return static_cast <Param_SpecDiffHybSolver *> (controlParameters);
}

}