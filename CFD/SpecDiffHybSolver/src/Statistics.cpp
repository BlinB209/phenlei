#include "Statistics.h"
#include "GlobalDataBase.h"
#include "PHMpi.h"
#include "TK_Exit.h"
#include "PHHeader.h"
#include "TK_Parse.h"
#include "Constants.h"
#include "SpecDiffHybSolver.h"
//#include "SpecSolver.h"
#include "PHIO.h"

using namespace std;

namespace PHSPACE
{
StatisticsHIT::StatisticsHIT()
{
}

StatisticsHIT::~StatisticsHIT()
{
    delete statisticAll; statisticAll = NULL;
    delete statisticLocal; statisticLocal = NULL;
    delete statistic; statistic = NULL;

    delete energyOfUUKx; energyOfUUKx = NULL;
    delete energyOfUUKxAll; energyOfUUKxAll = NULL;
    delete energyOfUUKxLocal; energyOfUUKxLocal = NULL;

    delete energyOfUUKy; energyOfUUKy = NULL;
    delete energyOfUUKyAll; energyOfUUKyAll = NULL;
    delete energyOfUUKyLocal; energyOfUUKyLocal = NULL;
}

Statistics::Statistics()
{
}

Statistics::~Statistics()
{
    delete statisticAll; statisticAll = NULL;
    delete statisticLocal; statisticLocal = NULL;
    delete statistic; statistic = NULL;

    delete nutAverage; nutAverage = NULL;

    delete energyOfUUKx; energyOfUUKx = NULL;
    delete energyOfUUKxAll; energyOfUUKxAll = NULL;
    delete energyOfUUKxLocal; energyOfUUKxLocal = NULL;

    delete energyOfUUKy; energyOfUUKy = NULL;
    delete energyOfUUKyAll; energyOfUUKyAll = NULL;
    delete energyOfUUKyLocal; energyOfUUKyLocal = NULL;
}

void StatisticsHIT::InitStatisticsData(SpecGrid *GridData)
{
    int nStatisticalStep = 0;
    GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);

    int istf = GridData->iStartFourierIndex;
    int iedf = GridData->iEndFourierIndex;
    iz_half = istf + (iedf - istf + 1) / 2;

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();

    if (myid == procIDWithKx0Ky0)
    {
        InitStatisticsFileName();
    }

    AllocStatisticsdata(GridData);

    InitStatisticsDataValue(GridData);

    int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");
    if (myid == procIDWithKx0Ky0)
    {
        //if (startStatisticStep < 0)
        if(IfStatisticExist())
        {
            InputStatistics(GridData);
        }
        else
        {
            ifInitStatistic = true;
        }

        WriteStatisticsName();
    }
}

void Statistics::InitStatisticsData(SpecDiffHybGrid *GridData)
{
    int nStatisticalStep = 0;
    GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);

    //inutavg = 0;

    tauWAverage = 0.0;

    int istf = GridData->iStartFourierIndex;
    int iedf = GridData->iEndFourierIndex;
    iz_half = istf + (iedf - istf + 1) / 2;

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();

    if (myid == procIDWithKx0Ky0)
    {
        InitStatisticsFileName();
        OutputPlot3dGrid(GridData);
    }

    AllocStatisticsdata(GridData);
    InitStatisticsDataValue(GridData);

    int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");
    if (myid == procIDWithKx0Ky0)
    {
        //if (startStatisticStep < 0)
        if(IfStatisticExist())
        {
            InputStatistics(GridData);
        }
        else
        {
            ifInitStatistic = true;
        }

        WriteStatisticsName();
    }
}

void StatisticsHIT::InitStatisticsFileName()
{
    string outputdir = "./results";
    outputdir = GlobalDataBase::GetStrParaFromDB("outputdir");

    zGridPlot3dFile = outputdir + "/Zgrid.x";
    zHalfGridPlot3dFile = outputdir + "/Zgrid_half.x";
    statisticsNameFile = outputdir + "/Statistic.nam";
    statisticsRestartFile = outputdir + "/Statistic.rsta";
    statisticsRestartTimeFile = outputdir + "/Statistic_t.rsta";
    statisticsInformationFile = outputdir + "/Statistic.info";
}

void Statistics::OutputPlot3dGrid(SpecDiffHybGrid *GridData)
{
    int istf = GridData->iStartFourierIndex;
    int iedf = GridData->iEndFourierIndex;
    Int1D *fourierIndexSize = GridData->fourierIndexSize;
    RDouble1D *rZ = GridData->realZ;

    int wordwidth = 20;
    fstream file;
    file.open(zGridPlot3dFile.c_str(), ios_base::out);
    if (!file)
    {
        TK_Exit::ExceptionExit("could not open zGridPlot3dFile\n");
    }
    file << setw(wordwidth) << 1 << "\n";
    file << setw(wordwidth) << 1 << setw(wordwidth) << 1 << setw(wordwidth) << (*fourierIndexSize)(1) << "\n";
    for (int iz = istf; iz <= iedf; ++ iz)
    {
        file << setiosflags(ios::left);
        file << setiosflags(ios::scientific);
        file << setprecision(10);
        file << setw(wordwidth);
        file << 0.0;
    }
    for (int iz = istf; iz <= iedf; ++ iz)
    {
        file << setiosflags(ios::left);
        file << setiosflags(ios::scientific);
        file << setprecision(10);
        file << setw(wordwidth);
        file << 0.0;
    }
    for (int iz = istf; iz <= iedf; ++ iz)
    {
        file << setiosflags(ios::left);
        file << setiosflags(ios::scientific);
        file << setprecision(10);
        file << setw(wordwidth);
        file << (*rZ)(iz);
    }
    file.close();
    file.clear();

    int nout = iz_half - istf + 1;

    file.open(zHalfGridPlot3dFile.c_str(), ios_base::out);
    if (!file)
    {
        TK_Exit::ExceptionExit("could not open zHalfGridPlot3dFile\n");
    }
    file << setw(wordwidth) << 1 << "\n";
    file << setw(wordwidth) << 1 << setw(wordwidth) << 1 << setw(wordwidth) << nout << "\n";
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setiosflags(ios::left);
        file << setiosflags(ios::scientific);
        file << setprecision(10);
        file << setw(wordwidth);
        file << 0.0;
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setiosflags(ios::left);
        file << setiosflags(ios::scientific);
        file << setprecision(10);
        file << setw(wordwidth);
        file << 0.0;
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setiosflags(ios::left);
        file << setiosflags(ios::scientific);
        file << setprecision(10);
        file << setw(wordwidth);
        file << (*rZ)(iz);
    }
    file.close();
    file.clear();
}

void StatisticsHIT::AllocStatisticsdata(SpecGrid *GridData)
{
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

    int nY = GridData -> nY;
    int nZ = GridData -> nZ;

    Range I(istf, iedf);
    Range J(jstf, jedf);
    Range K(kstf, kedf);

    Range XN(kStartFourierIndexGlobal, kEndFourierIndexGlobal);
    Range YN(jStartFourierIndexGlobal, jEndFourierIndexGlobal);
    Range ZN(iStartFourierIndexGlobal, iEndFourierIndexGlobal);

    Range YNHalf(jStartFourierIndexGlobal, nY/2+1);
    Range ZNHalf(iStartFourierIndexGlobal, nZ/2+1);

    Range N19(1, 19);
    Range N9(1, 9);
    Range N3(1, 3);

    statisticAll = new RDouble1D(N19, fortranArray);
    statisticLocal = new RDouble1D(N19, fortranArray);
    statistic = new RDouble1D(N19, fortranArray);

    energyOfUUKx = new RDouble2D(XN, N3, fortranArray);            
    energyOfUUKxAll = new RDouble2D(XN, N3, fortranArray);       
    energyOfUUKxLocal = new RDouble2D(K, N3, fortranArray); 
    
    energyOfUUKy = new RDouble2D(YNHalf, N3, fortranArray);
    energyOfUUKyAll = new RDouble2D(YN, N3, fortranArray);
    energyOfUUKyLocal = new RDouble2D(J, N3, fortranArray); 

    energyOfUUKz = new RDouble2D(ZNHalf, N3, fortranArray);
    energyOfUUKzAll = new RDouble2D(ZN, N3, fortranArray);
    energyOfUUKzLocal = new RDouble2D(I, N3, fortranArray);
}

void Statistics::AllocStatisticsdata(SpecDiffHybGrid *GridData)
{
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

    int nY = GridData -> nY;

    Range I(istf, iedf);
    Range J(jstf, jedf);
    Range K(kstf, kedf);

    Range IHalf(istf, iz_half);

    Range XN(kStartFourierIndexGlobal, kEndFourierIndexGlobal);
    Range YN(jStartFourierIndexGlobal, jEndFourierIndexGlobal);

    Range YNHalf(jStartFourierIndexGlobal, nY/2+1);

    Range N19(1, 19);
    Range N9(1, 9);
    Range N3(1, 3);

    statisticAll = new RDouble2D(I, N19, fortranArray);
    statisticLocal = new RDouble2D(I, N19, fortranArray);
    statistic = new RDouble2D(IHalf, N19, fortranArray);

    nutAverage = new RDouble1D(I, fortranArray);

    energyOfUUKx = new RDouble3D(IHalf, XN, N3, fortranArray);            
    energyOfUUKxAll = new RDouble3D(I, XN, N3, fortranArray);       
    energyOfUUKxLocal = new RDouble3D(I, K, N3, fortranArray); 
    
    energyOfUUKy = new RDouble3D(IHalf, YNHalf, N3, fortranArray);
    energyOfUUKyAll = new RDouble3D(I, YN, N3, fortranArray);
    energyOfUUKyLocal = new RDouble3D(I, J, N3, fortranArray); 
}

void StatisticsHIT::InitStatisticsDataValue(SpecGrid *GridData)
{
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

    int nY = GridData -> nY;
    int nZ = GridData -> nZ;

    Range I(istf, iedf);
    Range J(jstf, jedf);
    Range K(kstf, kedf);

    Range IHalf(istf, iz_half);

    Range XN(kStartFourierIndexGlobal, kEndFourierIndexGlobal);
    Range YN(jStartFourierIndexGlobal, jEndFourierIndexGlobal);
    Range ZN(iStartFourierIndexGlobal, iEndFourierIndexGlobal);

    Range YNHalf(jStartFourierIndexGlobal, nY/2+1);
    Range ZNHalf(iStartFourierIndexGlobal, nZ/2+1);

    Range N19(1, 19);
    Range N9(1, 9);
    Range N3(1, 3);

    (*statisticAll)(N19) = 0.0;
    (*statisticLocal)(N19) = 0.0;
    (*statistic)(N19) = 0.0;

    (*energyOfUUKx)(XN, N3) = 0.0;         
    (*energyOfUUKxAll)(XN, N3) = 0.0;       
    (*energyOfUUKxLocal)(K, N3) = 0.0; 

    (*energyOfUUKy)(YNHalf, N3) = 0.0;
    (*energyOfUUKyAll)(YN, N3) = 0.0;
    (*energyOfUUKyLocal)(J, N3) = 0.0;

    (*energyOfUUKz)(ZNHalf, N3) = 0.0;
    (*energyOfUUKzAll)(ZN, N3) = 0.0;
    (*energyOfUUKzLocal)(I, N3) = 0.0;
}

void Statistics::InitStatisticsDataValue(SpecDiffHybGrid *GridData)
{
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

    int nY = GridData -> nY;

    Range I(istf, iedf);
    Range J(jstf, jedf);
    Range K(kstf, kedf);

    Range IHalf(istf, iz_half);

    Range XN(kStartFourierIndexGlobal, kEndFourierIndexGlobal);
    Range YN(jStartFourierIndexGlobal, jEndFourierIndexGlobal);

    Range YNHalf(jStartFourierIndexGlobal, nY/2+1);
    Range N19(1, 19);
    Range N9(1, 9);
    Range N3(1, 3);

    (*statisticAll)(I, N19) = 0.0;
    (*statisticLocal)(I, N19) = 0.0;
    (*statistic)(IHalf, N19) = 0.0;

    (*nutAverage)(I) = 0.0;

    (*energyOfUUKx)(IHalf, XN, N3) = 0.0;         
    (*energyOfUUKxAll)(I, XN, N3) = 0.0;       
    (*energyOfUUKxLocal)(I, K, N3) = 0.0; 
    
    (*energyOfUUKy)(IHalf, YNHalf, N3) = 0.0;
    (*energyOfUUKyAll)(I, YN, N3) = 0.0;
    (*energyOfUUKyLocal)(I, J, N3) = 0.0;  
}

void StatisticsHIT::InputStatistics(SpecGrid *GridData)
{
    int nstart = GlobalDataBase::GetIntParaFromDB("nstart");
    if (!nstart)
    {
        ifInitStatistic = true;
        return;
    }

    ReadStatisticsInformation();

    ReadStatisticRestart(GridData);
}

void Statistics::InputStatistics(SpecDiffHybGrid *GridData)
{
    int nstart = GlobalDataBase::GetIntParaFromDB("nstart");
    if (!nstart)
    {
        ifInitStatistic = true;
        return;
    }

    ReadStatisticsInformation();

    ReadStatisticRestart(GridData);
}

void StatisticsHIT::WriteStatisticsName()
{
    fstream file;
    file.open(statisticsNameFile.c_str(), ios_base::out);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(statisticsNameFile);
    }

    file << "Uavg; mean velocity\n";
    file << "Vavg\n";
    file << "Wavg\n";
    file << "<uu>\n";
    file << "<vv>\n";
    file << "<ww>\n";
    file << "<uv>\n";
    file << "<uw>\n";
    file << "<vw>\n";
    file << "dissipationRate\n";

    file.close();
    file.clear();
}

void StatisticsHIT::ReadStatisticsInformation()
{
    int statisticMethod = GlobalDataBase::GetIntParaFromDB("statisticMethod");

    int istatMethTmp = 0;
    int nStatisticalStep = 0;

    string line,word;
    string separator  = " =\t\r\n#$,;\"'";
    fstream file;
    file.open(statisticsInformationFile.c_str(), ios_base::in);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(statisticsInformationFile);
    }

    ReadNewLine(file, line);
    line = FindNextWord(line, word, separator);
    from_string<int>(istatMethTmp, word, std::dec);
    line = FindNextWord(line, word, separator);
    from_string<int>(nStatisticalStep, word, std::dec);

    if (istatMethTmp != statisticMethod)
    {
        TK_Exit::ExceptionExit("istatMeth in statisticsInformationFile is NOT equal to that in para.inp!");
    }

    GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);

    file.close();
    file.clear();
}

void StatisticsHIT::ReadStatisticRestart(SpecGrid *GridData)
{
    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;
    int jStartFourierIndexGlobal = GridData -> jStartFourierIndexGlobal;
    int iStartFourierIndexGlobal = GridData -> iStartFourierIndexGlobal;

    int kEndFourierIndexGlobal = GridData -> kEndFourierIndexGlobal;

    int nY = GridData -> nY;
    int nZ = GridData -> nZ;

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;

    //int noutput1 = iz_half - istf + 1;
    int nOutputX = 3 * (kEndFourierIndexGlobal - kStartFourierIndexGlobal + 1);
    int nOutputY = 3 * (nY/2 + 1 - jStartFourierIndexGlobal + 1);
    int noutputZ = 3 * (nZ/2 + 1 - iStartFourierIndexGlobal + 1);

    fstream file;
    ios_base::openmode openMode = ios_base::in | ios_base::binary;
    PHSPACE::OpenFile(file, statisticsRestartFile.c_str(), openMode);

    PHRead(file, &(*statistic)(1), 19);

    PHRead(file, &(*energyOfUUKx)(kStartFourierIndexGlobal, 1), nOutputX);

    PHRead(file, &(*energyOfUUKy)(jStartFourierIndexGlobal, 1), nOutputY);

    PHRead(file, &(*energyOfUUKz)(iStartFourierIndexGlobal, 1), noutputZ);

    PHSPACE::CloseFile(file);
}

void Statistics::ReadStatisticRestart(SpecDiffHybGrid *GridData)
{
    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;
    int jStartFourierIndexGlobal = GridData -> jStartFourierIndexGlobal;
    int kEndFourierIndexGlobal = GridData -> kEndFourierIndexGlobal;

    int nY = GridData -> nY;

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;

    int noutput1 = iz_half - istf + 1;
    int noutput2 = (iz_half - istf + 1) * (kEndFourierIndexGlobal - kStartFourierIndexGlobal + 1);
    int noutput3 = (iz_half - istf + 1) * (nY/2 + 1 - jStartFourierIndexGlobal + 1);

    fstream file;
    ios_base::openmode openMode = ios_base::in | ios_base::binary;
    PHSPACE::OpenFile(file, statisticsRestartFile.c_str(), openMode);

    PHRead(file, &(*statistic)(istf, 1), 19 * noutput1);

    PHRead(file, &(*energyOfUUKx)(istf, kStartFourierIndexGlobal, 1), 3 * noutput2);

    PHRead(file, &(*energyOfUUKy)(istf, jStartFourierIndexGlobal, 1), 3 * noutput3);

    PHRead(file, &tauWAverage, 1);

    PHRead(file, &(*nutAverage)(istf), iedf - istf + 1);

    PHSPACE::CloseFile(file);
}

/*
void Statistics::ReadStatisticsRestart(SpecDiffHybGrid *GridData)
{
    int istf = GridData -> iStartFourierIndex;

    string line,word;
    string separator = " =\t\r\n#$,;\"'";
    fstream file;
    file.open(statisticsRestartFile.c_str(), ios_base::in);
    if (!file)
    {
        cout << "Warning!!!" << statisticsRestartFile << "does not exist or end during reading!" << "The statisic will begin from scratch!" << "\n" << endl;
        ifInitStatistic = true;
        file.close();
        file.clear();
    }
    else
    {
        ReadNewLine(file, line);
        ReadNewLine(file, line);

        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*uAverage)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*vAverage)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*wAverage)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*uSquare)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*vSquare)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*wSquare)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*uMultiplyV)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*uMultiplyW)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*vMultiplyW)(iz);
        }
        for (int iz = istf; iz <= iz_half; ++ iz)
        {
            file >> (*epsilon)(iz);
        }
        file.close();
        file.clear();

        ReadStatisticsSelected(GridData);
    }
}

void Statistics::ReadStatisticsSelected(SpecDiffHybGrid *GridData)
{
    int istf = GridData -> iStartFourierIndex;
    int x0Fall = GridData -> x0Fall;
    int xNFall = GridData -> xNFall;
    int y0Fall = GridData -> y0Fall;
    int nY = GridData -> nY;
    Range IHalf(istf, iz_half);
    Range XN(x0Fall, xNFall);
    Range YNHalf(y0Fall, nY/2+1);
    Range N9(1, 9);

    Param_SpecDiffHybSolver *parameters = GetControlParameters();
    RDouble reynolds = parameters->GetRefReNumber();
    RDouble dreynolds2 = 1.0 / reynolds / reynolds;

    int itmp = 0;
    RDouble rtmp = 0.0;

    string line,word;
    string separator  = " =\t\r\n#$,;\"'";
    fstream file;
    file.open("statistic_selected.dat", ios_base::in);
    if (!file)
    {
        TK_Exit::ExceptionExit("could not open statistic_selected.dat\n");
    }

    ReadNewLine(file, line);
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        ReadNewLine(file, line);

        line = FindNextWord(line, word, separator);
        line = FindNextWord(line, word, separator);
        line = FindNextWord(line, word, separator);
        line = FindNextWord(line, word, separator);
        line = FindNextWord(line, word, separator);
        line = FindNextWord(line, word, separator);
        line = FindNextWord(line, word, separator);

        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 1), word, std::dec );
        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 2), word, std::dec );
        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 3), word, std::dec );
        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 4), word, std::dec );
        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 5), word, std::dec );
        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 6), word, std::dec );
        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 7), word, std::dec );
        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 8), word, std::dec );
        line = FindNextWord(line, word, separator);
        from_string< double >( (*epsilonXYZ)(iz, 9), word, std::dec );
    }
    (*epsilonXYZ)(IHalf, N9) = (*epsilonXYZ)(IHalf, N9) / dreynolds2;
    file.close();
    file.clear();

    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        // !Read rEuu_kx
        string cheuuf = "./spectra/rEuu_kx_";
        stringstream strmMyid;
        strmMyid << iz;
        cheuuf = cheuuf +  strmMyid.str() + ".dat";
        file.open(cheuuf.c_str(), ios_base::in);
        if (!file)
        {
            cout << "Error in reading in " << cheuuf << "\n";
            cout << "rEuu_kx=rEuu_ky=0.0 will be used !!!\n" << endl;
            (*energyOfUUKx)(IHalf, XN) = 0.0;
            (*energyOfUUKy)(IHalf, YNHalf) = 0.0;
        }

        ReadNewLine(file, line);
        ReadNewLine(file, line);
        for (int ix = x0Fall; ix <= xNFall; ++ ix)
        {
            ReadNewLine(file, line);

            line = FindNextWord(line, word, separator);

            line = FindNextWord(line, word, separator);
            from_string< double >( (*energyOfUUKx)(iz, ix, 1), word, std::dec );
            line = FindNextWord(line, word, separator);
            from_string< double >( (*energyOfUUKx)(iz, ix, 2), word, std::dec );
            line = FindNextWord(line, word, separator);
            from_string< double >( (*energyOfUUKx)(iz, ix, 3), word, std::dec );
        }

        file.close();
        file.clear();

        //! Read rEuu_ky.
        cheuuf = "./spectra/rEuu_ky_";
        cheuuf = cheuuf +  strmMyid.str() + ".dat";
        file.open(cheuuf.c_str(), ios_base::in);
        if (!file)
        {
            cout << "Error in reading in " << cheuuf << "\n";
            cout << "rEuu_kx=rEuu_ky=0.0 will be used !!!\n" << endl;
            (*energyOfUUKx)(IHalf, XN) = 0.0;
            (*energyOfUUKy)(IHalf, YNHalf) = 0.0;
        }

        ReadNewLine(file, line);
        ReadNewLine(file, line);
        for (int iy = y0Fall; iy <= (nY/2+1); ++ iy)
        {
            ReadNewLine(file, line);

            line = FindNextWord(line, word, separator);

            line = FindNextWord(line, word, separator);
            from_string< double >( (*energyOfUUKy)(iz, iy, 1), word, std::dec );
            line = FindNextWord(line, word, separator);
            from_string< double >( (*energyOfUUKy)(iz, iy, 2), word, std::dec );
            line = FindNextWord(line, word, separator);
            from_string< double >( (*energyOfUUKy)(iz, iy, 3), word, std::dec );
        }

        file.close();
        file.clear();
    }
}
*/

void StatisticsHIT::GetStatisticsLocal(SpecGrid *GridData, Complex3D *complexU, Complex3D *complexV, Complex3D *complexW)
{
    int istf = GridData -> iStartFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    int iStartFourierIndexGlobal = GridData -> iStartFourierIndexGlobal;
    int jStartFourierIndexGlobal = GridData -> jStartFourierIndexGlobal;
    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;

    RDouble1D *kx = GridData -> realKX;
    RDouble1D *ky = GridData -> realKY;
    RDouble1D *kz = GridData -> realKZ;

    if (kstf == kStartFourierIndexGlobal)
    {
        PHComplex sumUU = PHComplex(0.0, 0.0);
        PHComplex sumVV = PHComplex(0.0, 0.0);
        PHComplex sumWW = PHComplex(0.0, 0.0);
        PHComplex sumUV = PHComplex(0.0, 0.0);
        PHComplex sumUW = PHComplex(0.0, 0.0);
        PHComplex sumVW = PHComplex(0.0, 0.0);

        for(int iy = jstf; iy <= jedf; ++ iy)
        {
            for(int iz = istf; iz <= iedf; ++iz)
            {
                sumUU = sumUU + (*complexU)(iz, iy, kstf) * conj((*complexU)(iz, iy, kstf));
                sumVV = sumVV + (*complexV)(iz, iy, kstf) * conj((*complexV)(iz, iy, kstf));
                sumWW = sumWW + (*complexW)(iz, iy, kstf) * conj((*complexW)(iz, iy, kstf));
                sumUV = sumUV + (*complexU)(iz, iy, kstf) * conj((*complexV)(iz, iy, kstf));
                sumUW = sumUW + (*complexU)(iz, iy, kstf) * conj((*complexW)(iz, iy, kstf));
                sumVW = sumVW + (*complexV)(iz, iy, kstf) * conj((*complexW)(iz, iy, kstf));
            }
        }

        for(int ix = kstf + 1; ix <= kedf; ++ ix)
        {
            for(int iy = jstf; iy <= jedf; ++ iy)
            {
                for(int iz = istf; iz <= iedf; ++iz)
                {
                    sumUU = sumUU + 2.0 * ( (*complexU)(iz, iy, kstf) * conj((*complexU)(iz, iy, kstf)) );
                    sumVV = sumVV + 2.0 * ( (*complexV)(iz, iy, kstf) * conj((*complexV)(iz, iy, kstf)) );
                    sumWW = sumWW + 2.0 * ( (*complexW)(iz, iy, kstf) * conj((*complexW)(iz, iy, kstf)) );
                    sumUV = sumUV + 2.0 * ( (*complexU)(iz, iy, kstf) * conj((*complexV)(iz, iy, kstf)) );
                    sumUW = sumUW + 2.0 * ( (*complexU)(iz, iy, kstf) * conj((*complexW)(iz, iy, kstf)) );
                    sumVW = sumVW + 2.0 * ( (*complexV)(iz, iy, kstf) * conj((*complexW)(iz, iy, kstf)) );
                }
            }
        }

        (*statisticLocal)(4) = real(sumUU);
        (*statisticLocal)(5) = real(sumVV);
        (*statisticLocal)(6) = real(sumWW);
        (*statisticLocal)(7) = real(sumUV);
        (*statisticLocal)(8) = real(sumUW);
        (*statisticLocal)(9) = real(sumVW);

        if (jstf == jStartFourierIndexGlobal && istf == iStartFourierIndexGlobal)
        {
            (*statisticLocal)(1) = real((*complexU)(istf, jstf, kstf));

            (*statisticLocal)(4) = (*statisticLocal)(4) - real((*complexU)(istf, jstf, kstf) * conj((*complexU)(istf, jstf, kstf)));
            (*statisticLocal)(5) = (*statisticLocal)(5) - real((*complexV)(istf, jstf, kstf) * conj((*complexV)(istf, jstf, kstf)));
            (*statisticLocal)(6) = (*statisticLocal)(6) - real((*complexW)(istf, jstf, kstf) * conj((*complexW)(istf, jstf, kstf)));
            (*statisticLocal)(7) = (*statisticLocal)(7) - real((*complexU)(istf, jstf, kstf) * conj((*complexV)(istf, jstf, kstf)));
            (*statisticLocal)(8) = (*statisticLocal)(8) - real((*complexU)(istf, jstf, kstf) * conj((*complexW)(istf, jstf, kstf)));
            (*statisticLocal)(9) = (*statisticLocal)(9) - real((*complexV)(istf, jstf, kstf) * conj((*complexW)(istf, jstf, kstf)));
        }

        // !dissipation
        int ix0 = kstf;
        for(int iy = jstf; iy <= jedf; ++iy)
        {
            for(int iz = istf; iz <= iedf; ++iz)
            {
                if(jstf == jStartFourierIndexGlobal && istf == iStartFourierIndexGlobal) continue;

                RDouble kx2 = (*kx)(ix0) * (*kx)(ix0);
                RDouble ky2 = (*ky)(iy) * (*ky)(iy);
                RDouble kz2 = (*kz)(iz) * (*ky)(iz);

                RDouble k2 = kx2 + ky2 + kz2;

                RDouble uu = real( (*complexU)(iz, iy, ix0) * conj((*complexU)(iz, iy, ix0)) );
                RDouble vv = real( (*complexV)(iz, iy, ix0) * conj((*complexV)(iz, iy, ix0)) );
                RDouble ww = real( (*complexW)(iz, iy, ix0) * conj((*complexW)(iz, iy, ix0)) );

                RDouble tmp = k2 * (uu + vv + ww);

                (*statisticLocal)(10) = (*statisticLocal)(10) + tmp;

                (*energyOfUUKxLocal)(ix0, 1) = (*energyOfUUKxLocal)(ix0, 1) + uu;
                (*statisticLocal)(11) = (*statisticLocal)(11) + kx2 * uu;

                (*energyOfUUKxLocal)(ix0, 2) = (*energyOfUUKxLocal)(ix0, 2) + vv;
                (*statisticLocal)(12) = (*statisticLocal)(12) + kx2 * vv;

                (*energyOfUUKxLocal)(ix0, 3) = (*energyOfUUKxLocal)(ix0, 3) + ww;
                (*statisticLocal)(13) = (*statisticLocal)(13) + kx2 * ww;

                (*energyOfUUKyLocal)(iy, 1) = (*energyOfUUKyLocal)(iy, 1) + uu;
                (*statisticLocal)(14) = (*statisticLocal)(14) + ky2 * uu;

                (*energyOfUUKyLocal)(iy, 2) = (*energyOfUUKyLocal)(iy, 2) + vv;
                (*statisticLocal)(15) = (*statisticLocal)(15) + ky2 * vv;

                (*energyOfUUKyLocal)(iy, 3) = (*energyOfUUKyLocal)(iy, 3) + ww;
                (*statisticLocal)(16) = (*statisticLocal)(16) + ky2 * ww;

                (*energyOfUUKzLocal)(iz, 1) = (*energyOfUUKzLocal)(iz, 1) + uu;
                (*statisticLocal)(17) = (*statisticLocal)(17) + kz2 * uu;

                (*energyOfUUKzLocal)(iz, 2) = (*energyOfUUKzLocal)(iz, 2) + vv;
                (*statisticLocal)(18) = (*statisticLocal)(18) + kz2 * vv;

                (*energyOfUUKzLocal)(iz, 3) = (*energyOfUUKzLocal)(iz, 3) + ww;
                (*statisticLocal)(19) = (*statisticLocal)(19) + kz2 * ww;
            }
        }

        for(int ix = kstf+1; ix <= kedf; ++ix)
        {
            for(int iy = jstf; iy <= jedf; ++iy)
            {
                for(int iz = istf; iz <= iedf; ++iz)
                {
                    RDouble kx2 = (*kx)(ix) * (*kx)(ix);
                    RDouble ky2 = (*ky)(iy) * (*ky)(iy);
                    RDouble kz2 = (*kz)(iz) * (*ky)(iz);

                    RDouble k2 = kx2 + ky2 + kz2;

                    RDouble uu = real( (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix)) );
                    RDouble vv = real( (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix)) );
                    RDouble ww = real( (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix)) );

                    RDouble tmp = k2 * (uu + vv + ww);

                    (*statisticLocal)(10) = (*statisticLocal)(10) + 2.0 * tmp;

                    (*energyOfUUKxLocal)(ix, 1) = (*energyOfUUKxLocal)(ix, 1) + 2.0 * uu;
                    (*statisticLocal)(11) = (*statisticLocal)(11) + 2.0 * kx2 * uu;

                    (*energyOfUUKxLocal)(ix, 2) = (*energyOfUUKxLocal)(ix, 2) + 2.0 * vv;
                    (*statisticLocal)(12) = (*statisticLocal)(12) + 2.0 * kx2 * vv;

                    (*energyOfUUKxLocal)(ix, 3) = (*energyOfUUKxLocal)(ix, 3) + 2.0 * ww;
                    (*statisticLocal)(13) = (*statisticLocal)(13) + 2.0 * kx2 * ww;

                    (*energyOfUUKyLocal)(iy, 1) = (*energyOfUUKyLocal)(iy, 1) + 2.0 * uu;
                    (*statisticLocal)(14) = (*statisticLocal)(14) + 2.0 * ky2 * uu;

                    (*energyOfUUKyLocal)(iy, 2) = (*energyOfUUKyLocal)(iy, 2) + 2.0 * vv;
                    (*statisticLocal)(15) = (*statisticLocal)(15) + 2.0 * ky2 * vv;

                    (*energyOfUUKyLocal)(iy, 3) = (*energyOfUUKyLocal)(iy, 3) + 2.0 * ww;
                    (*statisticLocal)(16) = (*statisticLocal)(16) + 2.0 * ky2 * ww;

                    (*energyOfUUKzLocal)(iz, 1) = (*energyOfUUKzLocal)(iz, 1) + 2.0 * uu;
                    (*statisticLocal)(17) = (*statisticLocal)(17) + 2.0 * kz2 * uu;

                    (*energyOfUUKzLocal)(iz, 2) = (*energyOfUUKzLocal)(iz, 2) + 2.0 * vv;
                    (*statisticLocal)(18) = (*statisticLocal)(18) + 2.0 * kz2 * vv;

                    (*energyOfUUKzLocal)(iz, 3) = (*energyOfUUKzLocal)(iz, 3) + 2.0 * ww;
                    (*statisticLocal)(19) = (*statisticLocal)(19) + 2.0 * kz2 * ww;
                }
            }
        }
    }

    else
    {
        PHComplex sumUU = PHComplex(0.0, 0.0);
        PHComplex sumVV = PHComplex(0.0, 0.0);
        PHComplex sumWW = PHComplex(0.0, 0.0);
        PHComplex sumUV = PHComplex(0.0, 0.0);
        PHComplex sumUW = PHComplex(0.0, 0.0);
        PHComplex sumVW = PHComplex(0.0, 0.0);

        for (int ix = kstf; ix <= kedf; ++ ix)
        {
            for (int iy = jstf; iy <= jedf; ++ iy)
            {
                for(int iz = istf; iz <= iedf; ++iz)
                {
                    sumUU = sumUU + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix));
                    sumVV = sumVV + 2.0 * (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix));
                    sumWW = sumWW + 2.0 * (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                    sumUV = sumUV + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexV)(iz, iy, ix));
                    sumUW = sumUW + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                    sumVW = sumVW + 2.0 * (*complexV)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                }
            }
        }

        (*statisticLocal)(4) = real(sumUU);
        (*statisticLocal)(5) = real(sumVV);
        (*statisticLocal)(6) = real(sumWW);
        (*statisticLocal)(7) = real(sumUV);
        (*statisticLocal)(8) = real(sumUW);
        (*statisticLocal)(9) = real(sumVW);

        for (int ix = kstf; ix <= kedf; ++ ix)
        {
            for (int iy = jstf; iy <= jedf; ++ iy)
            {
                for(int iz = istf; iz <= iedf; ++iz)
                {
                    RDouble kx2 = (*kx)(ix) * (*kx)(ix);
                    RDouble ky2 = (*ky)(iy) * (*ky)(iy);
                    RDouble kz2 = (*kz)(iz) * (*ky)(iz);

                    RDouble k2 = kx2 + ky2 + kz2;

                    RDouble uu = real( (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix)) );
                    RDouble vv = real( (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix)) );
                    RDouble ww = real( (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix)) );

                    RDouble tmp = k2 * (uu + vv + ww);

                    (*statisticLocal)(10) = (*statisticLocal)(10) + 2.0 * tmp;

                    (*energyOfUUKxLocal)(ix, 1) = (*energyOfUUKxLocal)(ix, 1) + 2.0 * uu;
                    (*statisticLocal)(11) = (*statisticLocal)(11) + 2.0 * kx2 * uu;

                    (*energyOfUUKxLocal)(ix, 2) = (*energyOfUUKxLocal)(ix, 2) + 2.0 * vv;
                    (*statisticLocal)(12) = (*statisticLocal)(12) + 2.0 * kx2 * vv;

                    (*energyOfUUKxLocal)(ix, 3) = (*energyOfUUKxLocal)(ix, 3) + 2.0 * ww;
                    (*statisticLocal)(13) = (*statisticLocal)(13) + 2.0 * kx2 * ww;

                    (*energyOfUUKyLocal)(iy, 1) = (*energyOfUUKyLocal)(iy, 1) + 2.0 * uu;
                    (*statisticLocal)(14) = (*statisticLocal)(14) + 2.0 * ky2 * uu;

                    (*energyOfUUKyLocal)(iy, 2) = (*energyOfUUKyLocal)(iy, 2) + 2.0 * vv;
                    (*statisticLocal)(15) = (*statisticLocal)(15) + 2.0 * ky2 * vv;

                    (*energyOfUUKyLocal)(iy, 3) = (*energyOfUUKyLocal)(iy, 3) + 2.0 * ww;
                    (*statisticLocal)(16) = (*statisticLocal)(16) + 2.0 * ky2 * ww;
                }
            }
        }
    }
}

void Statistics::GetStatisticsLocal(SpecDiffHybGrid *GridData, Complex3D *complexU, Complex3D *complexV, Complex3D *complexW, Complex4D *complexVelocityGradient)
{
    int istf = GridData -> iStartFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;
    int jStartFourierIndexGlobal = GridData -> jStartFourierIndexGlobal;
    RDouble1D *kx = GridData -> realKX;
    RDouble1D *ky = GridData -> realKY;

    Range I(istf, iedf);
    Range IHalf(istf, iz_half);
    Range J(jstf, jedf);
    Range K(kstf, kedf);
    Range N19(1, 19);
    Range N9(1, 9);
    Range N3(1, 3);

    (*statisticLocal)(I, N19) = 0.0;
    //(*epsXYZLocal)(I, N9) = 0.0;
    (*energyOfUUKxLocal)(I, K, N3) = 0.0;
    (*energyOfUUKyLocal)(I, J, N3) = 0.0;

    for (int iz = istf; iz <= iedf; ++ iz)
    {
        if (kstf == kStartFourierIndexGlobal)
        {
            // !process with (0,:)modes
            int iy0 = jstf;
            if (jstf == jStartFourierIndexGlobal)
            {
                // !process with (0,0) mode
                iy0 = jstf + 1;
                (*statisticLocal)(iz, 1) = real((*complexU)(iz, jstf, kstf));
            }

            // !Reynolds stress
            PHComplex sumUU = PHComplex(0.0, 0.0);
            PHComplex sumVV = PHComplex(0.0, 0.0);
            PHComplex sumWW = PHComplex(0.0, 0.0);
            PHComplex sumUV = PHComplex(0.0, 0.0);
            PHComplex sumUW = PHComplex(0.0, 0.0);
            PHComplex sumVW = PHComplex(0.0, 0.0);
            for (int iy = iy0; iy <= jedf; ++ iy)
            {
                sumUU = sumUU + (*complexU)(iz, iy, kstf) * conj((*complexU)(iz, iy, kstf));
                sumVV = sumVV + (*complexV)(iz, iy, kstf) * conj((*complexV)(iz, iy, kstf));
                sumWW = sumWW + (*complexW)(iz, iy, kstf) * conj((*complexW)(iz, iy, kstf));
                sumUV = sumUV + (*complexU)(iz, iy, kstf) * conj((*complexV)(iz, iy, kstf));
                sumUW = sumUW + (*complexU)(iz, iy, kstf) * conj((*complexW)(iz, iy, kstf));
                sumVW = sumVW + (*complexV)(iz, iy, kstf) * conj((*complexW)(iz, iy, kstf));
            }
            for (int ix = kstf + 1; ix <= kedf; ++ ix)
            {
                for (int iy = jstf; iy <= jedf; ++ iy)
                {
                    sumUU = sumUU + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix));
                    sumVV = sumVV + 2.0 * (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix));
                    sumWW = sumWW + 2.0 * (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                    sumUV = sumUV + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexV)(iz, iy, ix));
                    sumUW = sumUW + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                    sumVW = sumVW + 2.0 * (*complexV)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                }
            }
            (*statisticLocal)(iz, 4) = real(sumUU);
            (*statisticLocal)(iz, 5) = real(sumVV);
            (*statisticLocal)(iz, 6) = real(sumWW);
            (*statisticLocal)(iz, 7) = real(sumUV);
            (*statisticLocal)(iz, 8) = real(sumUW);
            (*statisticLocal)(iz, 9) = real(sumVW);

            // !dissipation
            int ix0 = kstf;
            for (int iy = iy0; iy <= jedf; ++ iy)
            {
                RDouble kx2 = (*kx)(ix0) * (*kx)(ix0);
                RDouble ky2 = (*ky)(iy) * (*ky)(iy);
                RDouble k2 = kx2 + ky2;

                RDouble uu = real( (*complexU)(iz, iy, ix0) * conj((*complexU)(iz, iy, ix0)) );
                RDouble vv = real( (*complexV)(iz, iy, ix0) * conj((*complexV)(iz, iy, ix0)) );
                RDouble ww = real( (*complexW)(iz, iy, ix0) * conj((*complexW)(iz, iy, ix0)) );

                PHComplex complexTmp = k2 * (uu + vv + ww)
                                     + (*complexVelocityGradient)(iz, iy, ix0, 3) * conj((*complexVelocityGradient)(iz, iy, ix0, 3))
                                     + (*complexVelocityGradient)(iz, iy, ix0, 6) * conj((*complexVelocityGradient)(iz, iy, ix0, 6))
                                     + (*complexVelocityGradient)(iz, iy, ix0, 9) * conj((*complexVelocityGradient)(iz, iy, ix0, 9));
                (*statisticLocal)(iz, 10) = (*statisticLocal)(iz, 10) + real(complexTmp);

                (*energyOfUUKxLocal)(iz, ix0, 1) = (*energyOfUUKxLocal)(iz, ix0, 1) + uu;
                (*statisticLocal)(iz, 11) = (*statisticLocal)(iz, 11) + kx2 * uu;

                (*energyOfUUKxLocal)(iz, ix0, 2) = (*energyOfUUKxLocal)(iz, ix0, 2) + vv;
                (*statisticLocal)(iz, 12) = (*statisticLocal)(iz, 12) + kx2 * vv;

                (*energyOfUUKxLocal)(iz, ix0, 3) = (*energyOfUUKxLocal)(iz, ix0, 3) + ww;
                (*statisticLocal)(iz, 13) = (*statisticLocal)(iz, 13) + kx2 * ww;

                (*energyOfUUKyLocal)(iz, iy, 1) = (*energyOfUUKyLocal)(iz, iy, 1) + uu;
                (*statisticLocal)(iz, 14) = (*statisticLocal)(iz, 14) + ky2 * uu;

                (*energyOfUUKyLocal)(iz, iy, 2) = (*energyOfUUKyLocal)(iz, iy, 2) + vv;
                (*statisticLocal)(iz, 15) = (*statisticLocal)(iz, 15) + ky2 * vv;

                (*energyOfUUKyLocal)(iz, iy, 3) = (*energyOfUUKyLocal)(iz, iy, 3) + ww;
                (*statisticLocal)(iz, 16) = (*statisticLocal)(iz, 16) + ky2 * ww;

                complexTmp = (*complexVelocityGradient)(iz, iy, ix0, 3) * conj((*complexVelocityGradient)(iz, iy, ix0, 3));
                (*statisticLocal)(iz, 17) = (*statisticLocal)(iz, 17) + real(complexTmp);

                complexTmp = (*complexVelocityGradient)(iz, iy, ix0, 6) * conj((*complexVelocityGradient)(iz, iy, ix0, 6));
                (*statisticLocal)(iz, 18) = (*statisticLocal)(iz, 18) + real(complexTmp);

                complexTmp = (*complexVelocityGradient)(iz, iy, ix0, 9) * conj((*complexVelocityGradient)(iz, iy, ix0, 9));
                (*statisticLocal)(iz, 19) = (*statisticLocal)(iz, 19) + real(complexTmp);
            }

            for (int ix = kstf+1; ix <= kedf; ++ ix)
            {
                for (int iy = jstf; iy <= jedf; ++ iy)
                {
                    RDouble kx2 = (*kx)(ix) * (*kx)(ix);
                    RDouble ky2 = (*ky)(iy) * (*ky)(iy);
                    RDouble k2 = kx2 + ky2;

                    RDouble uu = real( (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix)) );
                    RDouble vv = real( (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix)) );
                    RDouble ww = real( (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix)) );

                    PHComplex complexTmp = k2 * (uu + vv + ww)
                        + (*complexVelocityGradient)(iz, iy, ix, 3) * conj((*complexVelocityGradient)(iz, iy, ix, 3))
                        + (*complexVelocityGradient)(iz, iy, ix, 6) * conj((*complexVelocityGradient)(iz, iy, ix, 6))
                        + (*complexVelocityGradient)(iz, iy, ix, 9) * conj((*complexVelocityGradient)(iz, iy, ix, 9));
                    (*statisticLocal)(iz, 10) = (*statisticLocal)(iz, 10) + 2.0 * real(complexTmp);

                    (*energyOfUUKxLocal)(iz, ix, 1) = (*energyOfUUKxLocal)(iz, ix, 1) + 2.0 * uu;
                    (*statisticLocal)(iz, 11) = (*statisticLocal)(iz, 11) + 2.0 * kx2 * uu;

                    (*energyOfUUKxLocal)(iz, ix, 2) = (*energyOfUUKxLocal)(iz, ix, 2) + 2.0 * vv;
                    (*statisticLocal)(iz, 12) = (*statisticLocal)(iz, 12) + 2.0 * kx2 * vv;

                    (*energyOfUUKxLocal)(iz, ix, 3) = (*energyOfUUKxLocal)(iz, ix, 3) + 2.0 * ww;
                    (*statisticLocal)(iz, 13) = (*statisticLocal)(iz, 13) + 2.0 * kx2 * ww;

                    (*energyOfUUKyLocal)(iz, iy, 1) = (*energyOfUUKyLocal)(iz, iy, 1) + 2.0 * uu;
                    (*statisticLocal)(iz, 14) = (*statisticLocal)(iz, 14) + 2.0 * ky2 * uu;

                    (*energyOfUUKyLocal)(iz, iy, 2) = (*energyOfUUKyLocal)(iz, iy, 2) + 2.0 * vv;
                    (*statisticLocal)(iz, 15) = (*statisticLocal)(iz, 15) + 2.0 * ky2 * vv;

                    (*energyOfUUKyLocal)(iz, iy, 3) = (*energyOfUUKyLocal)(iz, iy, 3) + 2.0 * ww;
                    (*statisticLocal)(iz, 16) = (*statisticLocal)(iz, 16) + 2.0 * ky2 * ww;

                    complexTmp = (*complexVelocityGradient)(iz, iy, ix, 3) * conj((*complexVelocityGradient)(iz, iy, ix, 3));
                    (*statisticLocal)(iz, 17) = (*statisticLocal)(iz, 17) + 2.0 * real(complexTmp);

                    complexTmp = (*complexVelocityGradient)(iz, iy, ix, 6) * conj((*complexVelocityGradient)(iz, iy, ix, 6));
                    (*statisticLocal)(iz, 18) = (*statisticLocal)(iz, 18) + 2.0 * real(complexTmp);

                    complexTmp = (*complexVelocityGradient)(iz, iy, ix, 9) * conj((*complexVelocityGradient)(iz, iy, ix, 9));
                    (*statisticLocal)(iz, 19) = (*statisticLocal)(iz, 19) + 2.0 * real(complexTmp);
                }// !for iy
            }// !for ix             
        }// !if(kstf == x0Fall)
        else
        {
            // !process with kx!=0
            // !Reynolds stress
            PHComplex sumUU = PHComplex(0.0, 0.0);
            PHComplex sumVV = PHComplex(0.0, 0.0);
            PHComplex sumWW = PHComplex(0.0, 0.0);
            PHComplex sumUV = PHComplex(0.0, 0.0);
            PHComplex sumUW = PHComplex(0.0, 0.0);
            PHComplex sumVW = PHComplex(0.0, 0.0);
            for (int ix = kstf; ix <= kedf; ++ ix)
            {
                for (int iy = jstf; iy <= jedf; ++ iy)
                {
                    sumUU = sumUU + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix));
                    sumVV = sumVV + 2.0 * (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix));
                    sumWW = sumWW + 2.0 * (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                    sumUV = sumUV + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexV)(iz, iy, ix));
                    sumUW = sumUW + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                    sumVW = sumVW + 2.0 * (*complexV)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                }
            }
            (*statisticLocal)(iz, 4) = real(sumUU);
            (*statisticLocal)(iz, 5) = real(sumVV);
            (*statisticLocal)(iz, 6) = real(sumWW);
            (*statisticLocal)(iz, 7) = real(sumUV);
            (*statisticLocal)(iz, 8) = real(sumUW);
            (*statisticLocal)(iz, 9) = real(sumVW);

            // !dissipation
            // !epsilon = -(kx^2+ky^2)*( cu*conjg(cu) +cv*conjg(cv) cw*conjg(cw) ) + cdudz*conjg(cdudz) +cdvdz*conjg(cdvdz) +cdwdz*conjg(cdwdz)
            for (int ix = kstf; ix <= kedf; ++ ix)
            {
                for (int iy = jstf; iy <= jedf; ++ iy)
                {
                    RDouble kx2 = (*kx)(ix) * (*kx)(ix);
                    RDouble ky2 = (*ky)(iy) * (*ky)(iy);
                    RDouble k2 = kx2 + ky2;

                    RDouble uu = real( (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix)) );
                    RDouble vv = real( (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix)) );
                    RDouble ww = real( (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix)) );

                    PHComplex complexTmp = k2 * (uu + vv + ww)
                                         + (*complexVelocityGradient)(iz, iy, ix, 3) * conj((*complexVelocityGradient)(iz, iy, ix, 3))
                                         + (*complexVelocityGradient)(iz, iy, ix, 6) * conj((*complexVelocityGradient)(iz, iy, ix, 6))
                                         + (*complexVelocityGradient)(iz, iy, ix, 9) * conj((*complexVelocityGradient)(iz, iy, ix, 9));
                    (*statisticLocal)(iz, 10) = (*statisticLocal)(iz, 10) + 2.0 * real(complexTmp);

                    (*energyOfUUKxLocal)(iz, ix, 1) = (*energyOfUUKxLocal)(iz, ix, 1) + 2.0 * uu;
                    (*statisticLocal)(iz, 11) = (*statisticLocal)(iz, 11) + 2.0 * kx2 * uu;

                    (*energyOfUUKxLocal)(iz, ix, 2) = (*energyOfUUKxLocal)(iz, ix, 2) + 2.0 * vv;
                    (*statisticLocal)(iz, 12) = (*statisticLocal)(iz, 12) + 2.0 * kx2 * vv;

                    (*energyOfUUKxLocal)(iz, ix, 3) = (*energyOfUUKxLocal)(iz, ix, 3) + 2.0 * ww;
                    (*statisticLocal)(iz, 13) = (*statisticLocal)(iz, 13) + 2.0 * kx2 * ww;

                    (*energyOfUUKyLocal)(iz, iy, 1) = (*energyOfUUKyLocal)(iz, iy, 1) + 2.0 * uu;
                    (*statisticLocal)(iz, 14) = (*statisticLocal)(iz, 14) + 2.0 * ky2 * uu;

                    (*energyOfUUKyLocal)(iz, iy, 2) = (*energyOfUUKyLocal)(iz, iy, 2) + 2.0 * vv;
                    (*statisticLocal)(iz, 15) = (*statisticLocal)(iz, 15) + 2.0 * ky2 * vv;

                    (*energyOfUUKyLocal)(iz, iy, 3) = (*energyOfUUKyLocal)(iz, iy, 3) + 2.0 * ww;
                    (*statisticLocal)(iz, 16) = (*statisticLocal)(iz, 16) + 2.0 * ky2 * ww;

                    complexTmp = (*complexVelocityGradient)(iz, iy, ix, 3) * conj((*complexVelocityGradient)(iz, iy, ix, 3));
                    (*statisticLocal)(iz, 17) = (*statisticLocal)(iz, 17) + 2.0 * real(complexTmp);

                    complexTmp = (*complexVelocityGradient)(iz, iy, ix, 6) * conj((*complexVelocityGradient)(iz, iy, ix, 6));
                    (*statisticLocal)(iz, 18) = (*statisticLocal)(iz, 18) + 2.0 * real(complexTmp);

                    complexTmp = (*complexVelocityGradient)(iz, iy, ix, 9) * conj((*complexVelocityGradient)(iz, iy, ix, 9));
                    (*statisticLocal)(iz, 19) = (*statisticLocal)(iz, 19) + 2.0 * real(complexTmp);
                }//! for iy
            }//! for ix
        }//! else
    }//! for iz
}

void StatisticsHIT::OutputStatisticRestart(SpecGrid *GridData)
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();

    if (myid != procIDWithKx0Ky0)
    {
        return;
    }
    else
    {
        WriteStatisticInformation();

        WriteStatisticRestart(GridData);

        //WriteStatisticVisual(GridData);

        //WriteSpectraVisual(GridData);
    }
}

void Statistics::OutputStatisticRestart(SpecDiffHybGrid *GridData)
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();

    if (myid != procIDWithKx0Ky0)
    {
        return;
    }
    else
    {
        WriteStatisticInformation();

        WriteStatisticRestart(GridData);

        //WriteStatisticVisual(GridData);

        //WriteSpectraVisual(GridData);
    }
}

void StatisticsHIT::OutputStatisticVisual(SpecGrid *GridData)
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();

    if (myid != procIDWithKx0Ky0)
    {
        return;
    }
    else
    {
        WriteStatisticVisual(GridData);

        WriteSpectraVisual(GridData);
    }
}

void Statistics::OutputStatisticVisual(SpecDiffHybGrid *GridData)
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();

    if (myid != procIDWithKx0Ky0)
    {
        return;
    }
    else
    {
        WriteStatisticVisual(GridData);

        WriteSpectraVisual(GridData);
    }
}

void StatisticsHIT::WriteStatisticInformation()
{
    int statisticMethod = GlobalDataBase::GetIntParaFromDB("statisticMethod");

    int nStatisticalStep = 0;
    GlobalDataBase::GetData("nStatisticalStep", &nStatisticalStep, PHINT, 1);

    fstream file;
    file.open(statisticsInformationFile.c_str(), ios_base::out);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(statisticsInformationFile);
    }
    file << setw(20) << statisticMethod << setw(20) << nStatisticalStep << "\n";
    file.close();
    file.clear();
}

void StatisticsHIT::WriteStatisticRestart(SpecGrid *GridData)
{
    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;
    int jStartFourierIndexGlobal = GridData -> jStartFourierIndexGlobal;
    int iStartFourierIndexGlobal = GridData -> iStartFourierIndexGlobal;

    int kEndFourierIndexGlobal = GridData -> kEndFourierIndexGlobal;

    int nY = GridData -> nY;
    int nZ = GridData -> nY;

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;

    //int noutput1 = iz_half - istf + 1;
    int nOutputX = 3 * (kEndFourierIndexGlobal - kStartFourierIndexGlobal + 1);
    int noutputY = 3 * (nY/2 + 1 - jStartFourierIndexGlobal + 1);
    int noutputZ = 3 * (nZ/2 + 1 - iStartFourierIndexGlobal + 1);

    fstream file;
    ios_base::openmode openMode = ios_base::out | ios_base::binary;
    PHSPACE::OpenFile(file, statisticsRestartFile.c_str(), openMode);

    PHWrite(file, &(*statistic)(1), 19);

    PHWrite(file, &(*energyOfUUKx)(kStartFourierIndexGlobal, 1), nOutputX);

    PHWrite(file, &(*energyOfUUKy)(jStartFourierIndexGlobal, 1), noutputY);

    PHWrite(file, &(*energyOfUUKz)(iStartFourierIndexGlobal, 1), noutputZ);

    PHSPACE::CloseFile(file);
}

void Statistics::WriteStatisticRestart(SpecDiffHybGrid *GridData)
{
    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;
    int jStartFourierIndexGlobal = GridData -> jStartFourierIndexGlobal;
    int kEndFourierIndexGlobal = GridData -> kEndFourierIndexGlobal;

    int nY = GridData -> nY;

    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;

    int noutput1 = iz_half - istf + 1;
    int noutput2 = (iz_half - istf + 1) * (kEndFourierIndexGlobal - kStartFourierIndexGlobal + 1);
    int noutput3 = (iz_half - istf + 1) * (nY/2 + 1 - jStartFourierIndexGlobal + 1);

    fstream file;
    ios_base::openmode openMode = ios_base::out | ios_base::binary;
    PHSPACE::OpenFile(file, statisticsRestartFile.c_str(), openMode);

    PHWrite(file, &(*statistic)(istf, 1), 19 * noutput1);

    PHWrite(file, &(*energyOfUUKx)(istf, kStartFourierIndexGlobal, 1), 3 * noutput2);

    PHWrite(file, &(*energyOfUUKy)(istf, jStartFourierIndexGlobal, 1), 3 * noutput3);

    PHWrite(file, &tauWAverage, 1);

    PHWrite(file, &(*nutAverage)(istf), iedf - istf + 1);

    PHSPACE::CloseFile(file);
}

/*
void Statistics::WriteStatisticsRestart(SpecDiffHybGrid *GridData)
{
    int istf = GridData -> iStartFourierIndex;
    int noutput = iz_half - istf + 1;

    //ActionKey *actkeyDumpStatistic = new ActionKey();
    //FillActionKey(actkeyDumpStatistic, DUMP_RESIDUAL, 0);
    //string formatStatistic;
    //PrepareFormatedStatistic(formatStatistic);
    //delete actkeyDumpStatistic;

    fstream file;
    file.open(statisticsRestartFile.c_str(), ios_base::out);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(statisticsRestartFile);
    }
    file << setw(20) << 1 << "\n";
    file << setw(20) << 1 << setw(20) << 1 << setw(20) << noutput << setw(20) << 10 << "\n";
    file << setiosflags(ios::right);
    file << setprecision(10);
    file << setiosflags(ios::scientific);
    file << setiosflags(ios::showpoint);
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*uAverage)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*vAverage)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*wAverage)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*uSquare)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*vSquare)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*wSquare)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*uMultiplyV)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*uMultiplyW)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*vMultiplyW)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*epsilon)(iz);
    }
    file.close();
    file.clear();

    file.open(statisticsRestartTimeFile.c_str(), ios_base::out);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(statisticsRestartTimeFile);
    }
    file << setw(20) << 1 << "\n";
    file << setw(20) << 1 << setw(20) << 1 << setw(20) << noutput << setw(20) << 10 << "\n";
    file << setiosflags(ios::right);
    file << setprecision(10);
    file << setiosflags(ios::scientific);
    file << setiosflags(ios::showpoint);
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*uAverageTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*vAverageTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*wAverageTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*uSquareTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*vSquareTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*wSquareTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*uMultiplyVTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*uMultiplyWTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*vMultiplyWTime)(iz);
    }
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(20) << (*epsilonTime)(iz);
    }
    file.close();
    file.clear();
}
*/

void StatisticsHIT::WriteStatisticVisual(SpecGrid *GridData)
{
    ActionKey *actkeyDumpStatistic = new ActionKey();
    FillActionKey(actkeyDumpStatistic, VISUALIZATION_AVERAGE_FLOW, 0);

    string formatStatistic;
    PrepareFormatedStatistic(GridData, formatStatistic);

    string statisticFileName = GetStatisticVisualFileName();
    if (statisticFileName == "")
    {
        return;
    }

    actkeyDumpStatistic->filename = statisticFileName;
    actkeyDumpStatistic->openmode = ios_base::out|ios_base::trunc;

    if (!actkeyDumpStatistic->IsNeedOpenFile())
    {
        return;
    }

    ParallelOpenFile(actkeyDumpStatistic);

    std::ostringstream oss;
    fstream &file = *(actkeyDumpStatistic->file);

    if (IfFileEmpty(file))
    {
        oss << "Title=\"STATISTIC\"" << endl;
        oss << "Variables=" << endl;
        oss << "<U>" << endl;
        oss << "<V>" << endl;
        oss << "<W>" << endl;
        oss << "<uu>" << endl;
        oss << "<vv>" << endl;
        oss << "<ww>" << endl;
        oss << "<uv>" << endl;
        oss << "<uw>" << endl;
        oss << "<vw>" << endl;
        oss << "epsilon" << endl;
        oss << "eps_ux" << endl;
        oss << "eps_vx" << endl;
        oss << "eps_wx" << endl;
        oss << "eps_uy" << endl;
        oss << "eps_vy" << endl;
        oss << "eps_wy" << endl;
        oss << "eps_uz" << endl;
        oss << "eps_vz" << endl;
        oss << "eps_wz" << endl;
    }

    oss << formatStatistic;

    WriteASCIIFile(file, oss.str());

    ParallelCloseFile(actkeyDumpStatistic);

    delete actkeyDumpStatistic;
}

void Statistics::WriteStatisticVisual(SpecDiffHybGrid *GridData)
{
    ActionKey *actkeyDumpStatistic = new ActionKey();
    FillActionKey(actkeyDumpStatistic, VISUALIZATION_AVERAGE_FLOW, 0);

    string formatStatistic;
    PrepareFormatedStatistic(GridData, formatStatistic);

    string statisticFileName = GetStatisticVisualFileName();
    if (statisticFileName == "")
    {
        return;
    }

    actkeyDumpStatistic->filename = statisticFileName;
    actkeyDumpStatistic->openmode = ios_base::out|ios_base::trunc;

    if (!actkeyDumpStatistic->IsNeedOpenFile())
    {
        return;
    }

    ParallelOpenFile(actkeyDumpStatistic);

    std::ostringstream oss;
    fstream &file = *(actkeyDumpStatistic->file);

    if (IfFileEmpty(file))
    {
        oss << "Title=\"STATISTIC\"" << endl;
        oss << "Variables=" << endl;
        oss << "i" << endl;
        oss << "Z" << endl;
        oss << "Z<sup>+</sup>" << endl;
        oss << "<U>" << endl;
        oss << "<V>" << endl;
        oss << "<W>" << endl;
        oss << "<uu>" << endl;
        oss << "<vv>" << endl;
        oss << "<ww>" << endl;
        oss << "<uv>" << endl;
        oss << "<uw>" << endl;
        oss << "<vw>" << endl;
        oss << "epsilon" << endl;
        oss << "eps_ux" << endl;
        oss << "eps_vx" << endl;
        oss << "eps_wx" << endl;
        oss << "eps_uy" << endl;
        oss << "eps_vy" << endl;
        oss << "eps_wy" << endl;
        oss << "eps_uz" << endl;
        oss << "eps_vz" << endl;
        oss << "eps_wz" << endl;
        oss << "<greek>n</greek><sub>t</sub>" << endl;
    }

    oss << formatStatistic;

    WriteASCIIFile(file, oss.str());

    ParallelCloseFile(actkeyDumpStatistic);

    delete actkeyDumpStatistic;
}

const string StatisticsHIT::GetStatisticVisualFileName()
{
    string statisticVisualFile = "statisticVisual.dat";
    GlobalDataBase::GetData("statisticVisualFile", &statisticVisualFile, PHSTRING, 1);
    return statisticVisualFile;
}

void StatisticsHIT::PrepareFormatedStatistic(SpecGrid *GridData, string &formatStatistic)
{
    RDouble reynolds = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    RDouble dreynolds2 = 1.0 / reynolds / reynolds;

    int istf = GridData -> iStartFourierIndex;

    ostringstream oss;

    oss << setiosflags(ios::left);
    oss << setprecision(5);
    oss << setiosflags(ios::scientific);
    oss << setiosflags(ios::showpoint);

    oss << (*statistic)(1) << "    ";
    oss << (*statistic)(2) << "    ";
    oss << (*statistic)(3) << "    ";
    oss << (*statistic)(4) << "    ";
    oss << (*statistic)(5) << "    ";
    oss << (*statistic)(6) << "    ";
    oss << (*statistic)(7) << "    ";
    oss << (*statistic)(8) << "    ";
    oss << (*statistic)(9) << "    ";
    oss << (*statistic)(10) * dreynolds2 << "    ";
    oss << (*statistic)(11) * dreynolds2 << "    ";
    oss << (*statistic)(12) * dreynolds2 << "    ";
    oss << (*statistic)(13) * dreynolds2 << "    ";
    oss << (*statistic)(14) * dreynolds2 << "    ";
    oss << (*statistic)(15) * dreynolds2 << "    ";
    oss << (*statistic)(16) * dreynolds2 << "    ";
    oss << (*statistic)(17) * dreynolds2 << "    ";
    oss << (*statistic)(18) * dreynolds2 << "    ";
    oss << (*statistic)(19) * dreynolds2 << "    ";
    oss << "\n";

    oss << endl;

    formatStatistic = oss.str();
}

void Statistics::PrepareFormatedStatistic(SpecDiffHybGrid *GridData, string &formatStatistic)
{
    RDouble reynolds = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    RDouble dreynolds2 = 1.0 / reynolds / reynolds;

    RDouble1D *realZ = GridData -> realZ;
    int istf = GridData -> iStartFourierIndex;

    ostringstream oss;

    oss << setiosflags(ios::left);
    oss << setprecision(5);
    oss << setiosflags(ios::scientific);
    oss << setiosflags(ios::showpoint);

    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        oss << iz << "    ";
        oss << (*realZ)(iz) << "    ";
        oss << (1.0 + (*realZ)(iz)) * reynolds << "    ";
        oss << (*statistic)(iz, 1) << "    ";
        oss << (*statistic)(iz, 2) << "    ";
        oss << (*statistic)(iz, 3) << "    ";
        oss << (*statistic)(iz, 4) << "    ";
        oss << (*statistic)(iz, 5) << "    ";
        oss << (*statistic)(iz, 6) << "    ";
        oss << (*statistic)(iz, 7) << "    ";
        oss << (*statistic)(iz, 8) << "    ";
        oss << (*statistic)(iz, 9) << "    ";
        oss << (*statistic)(iz, 10) * dreynolds2 << "    ";
        oss << (*statistic)(iz, 11) * dreynolds2 << "    ";
        oss << (*statistic)(iz, 12) * dreynolds2 << "    ";
        oss << (*statistic)(iz, 13) * dreynolds2 << "    ";
        oss << (*statistic)(iz, 14) * dreynolds2 << "    ";
        oss << (*statistic)(iz, 15) * dreynolds2 << "    ";
        oss << (*statistic)(iz, 16) * dreynolds2 << "    ";
        oss << (*statistic)(iz, 17) * dreynolds2 << "    ";
        oss << (*statistic)(iz, 18) * dreynolds2 << "    ";
        oss << (*statistic)(iz, 19) * dreynolds2 << "    ";
        oss << (*nutAverage)(iz) << "    ";
        oss << "\n";
    }
    oss << endl;

    formatStatistic = oss.str();
}

void StatisticsHIT::WriteSpectraVisual(SpecGrid *GridData)
{
    WriteSpectraKxVisual(GridData);

    WriteSpectraKyVisual(GridData);
}

void Statistics::WriteSpectraVisual(SpecDiffHybGrid *GridData)
{
    WriteSpectraKxVisual(GridData);

    WriteSpectraKyVisual(GridData);
}

void StatisticsHIT::WriteSpectraKxVisual(SpecGrid *GridData)
{
    ActionKey *actkeyDumpSpectra = new ActionKey();
    FillActionKey(actkeyDumpSpectra, VISUALIZATION_AVERAGE_FLOW, 0);

    string spectraFileName = "results/Euu_kx";

    actkeyDumpSpectra->filename = spectraFileName;
    actkeyDumpSpectra->openmode = ios_base::out|ios_base::trunc;

    if (!actkeyDumpSpectra->IsNeedOpenFile())
    {
        return;
    }

    ParallelOpenFile(actkeyDumpSpectra);

    std::ostringstream oss;
    fstream &file = *(actkeyDumpSpectra->file);

    if (IfFileEmpty(file))
    {
        oss << "Title=\"SPECTRA\"" << endl;
        oss << "Variables=" << endl;
        oss << "kx" << endl;
        oss << "Euu" << endl;
        oss << "Evv" << endl;
        oss << "Eww" << endl;
    }

    string formatSpectra;
    PrepareFormatedSpectraKx(GridData, formatSpectra);

    oss << formatSpectra;

    WriteASCIIFile(file, oss.str());

    ParallelCloseFile(actkeyDumpSpectra);

    delete actkeyDumpSpectra;
}

void Statistics::WriteSpectraKxVisual(SpecDiffHybGrid *GridData)
{
    ActionKey *actkeyDumpSpectra = new ActionKey();
    FillActionKey(actkeyDumpSpectra, VISUALIZATION_AVERAGE_FLOW, 0);

    string spectraFileName = "results/Euu_kx";

    actkeyDumpSpectra->filename = spectraFileName;
    actkeyDumpSpectra->openmode = ios_base::out|ios_base::trunc;

    if (!actkeyDumpSpectra->IsNeedOpenFile())
    {
        return;
    }

    ParallelOpenFile(actkeyDumpSpectra);

    std::ostringstream oss;
    fstream &file = *(actkeyDumpSpectra->file);

    if (IfFileEmpty(file))
    {
        oss << "Title=\"SPECTRA\"" << endl;
        oss << "Variables=" << endl;
        oss << "kx" << endl;
        oss << "Euu" << endl;
        oss << "Evv" << endl;
        oss << "Eww" << endl;
    }

    string formatSpectra;
    PrepareFormatedSpectraKx(GridData, formatSpectra);

    oss << formatSpectra;

    WriteASCIIFile(file, oss.str());

    ParallelCloseFile(actkeyDumpSpectra);

    delete actkeyDumpSpectra;
}

void StatisticsHIT::PrepareFormatedSpectraKx(SpecGrid *GridData, string &formatSpectra)
{
    RDouble reynolds = GlobalDataBase::GetDoubleParaFromDB("refReNumber");

    RDouble realXL = GlobalDataBase::GetDoubleParaFromDB("XL");       

    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;
    int kEndFourierIndexGlobal = GridData -> kEndFourierIndexGlobal;

    ostringstream oss;

    oss << setiosflags(ios::left);
    oss << setprecision(5);
    oss << setiosflags(ios::scientific);
    oss << setiosflags(ios::showpoint);

    for (int ix = kStartFourierIndexGlobal; ix <= kEndFourierIndexGlobal; ++ ix)
    {
        RDouble rk_x = static_cast<double>(ix - kStartFourierIndexGlobal) / realXL;
        oss << rk_x << "    ";
        oss << (*energyOfUUKx)(ix, 1) << "    ";
        oss << (*energyOfUUKx)(ix, 2) << "    ";
        oss << (*energyOfUUKx)(ix, 3) << "    ";
        oss << "\n";
    }

    oss << endl;

    formatSpectra = oss.str();
}

void Statistics::PrepareFormatedSpectraKx(SpecDiffHybGrid *GridData, string &formatSpectra)
{
    RDouble reynolds = GlobalDataBase::GetDoubleParaFromDB("refReNumber");

    RDouble realXL = GlobalDataBase::GetDoubleParaFromDB("XL");       

    RDouble1D *realZ = GridData -> realZ;

    int istf = GridData -> iStartFourierIndex;

    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;
    int kEndFourierIndexGlobal = GridData -> kEndFourierIndexGlobal;

    ostringstream oss;

    oss << setiosflags(ios::left);
    oss << setprecision(5);
    oss << setiosflags(ios::scientific);
    oss << setiosflags(ios::showpoint);

    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        oss << "zone t = \"yplus = " << (1.0 + (*realZ)(iz)) * reynolds << "\"\n";

        for (int ix = kStartFourierIndexGlobal; ix <= kEndFourierIndexGlobal; ++ ix)
        {
            RDouble rk_x = static_cast<double>(ix - kStartFourierIndexGlobal) / realXL;
            oss << rk_x << "    ";
            oss << (*energyOfUUKx)(iz, ix, 1) << "    ";
            oss << (*energyOfUUKx)(iz, ix, 2) << "    ";
            oss << (*energyOfUUKx)(iz, ix, 3) << "    ";
            oss << "\n";
        }
    }

    oss << endl;

    formatSpectra = oss.str();
}

void StatisticsHIT::WriteSpectraKyVisual(SpecGrid *GridData)
{
    ActionKey *actkeyDumpSpectra = new ActionKey();
    FillActionKey(actkeyDumpSpectra, VISUALIZATION_AVERAGE_FLOW, 0);

    string spectraFileName = "results/Euu_ky";

    actkeyDumpSpectra->filename = spectraFileName;
    actkeyDumpSpectra->openmode = ios_base::out|ios_base::trunc;

    if (!actkeyDumpSpectra->IsNeedOpenFile())
    {
        return;
    }

    ParallelOpenFile(actkeyDumpSpectra);

    std::ostringstream oss;
    fstream &file = *(actkeyDumpSpectra->file);

    if (IfFileEmpty(file))
    {
        oss << "Title=\"SPECTRA\"" << endl;
        oss << "Variables=" << endl;
        oss << "ky" << endl;
        oss << "Euu" << endl;
        oss << "Evv" << endl;
        oss << "Eww" << endl;
    }

    string formatSpectra;
    PrepareFormatedSpectraKy(GridData, formatSpectra);

    oss << formatSpectra;

    WriteASCIIFile(file, oss.str());

    ParallelCloseFile(actkeyDumpSpectra);

    delete actkeyDumpSpectra;
}

void Statistics::WriteSpectraKyVisual(SpecDiffHybGrid *GridData)
{
    ActionKey *actkeyDumpSpectra = new ActionKey();
    FillActionKey(actkeyDumpSpectra, VISUALIZATION_AVERAGE_FLOW, 0);

    string spectraFileName = "results/Euu_ky";

    actkeyDumpSpectra->filename = spectraFileName;
    actkeyDumpSpectra->openmode = ios_base::out|ios_base::trunc;

    if (!actkeyDumpSpectra->IsNeedOpenFile())
    {
        return;
    }

    ParallelOpenFile(actkeyDumpSpectra);

    std::ostringstream oss;
    fstream &file = *(actkeyDumpSpectra->file);

    if (IfFileEmpty(file))
    {
        oss << "Title=\"SPECTRA\"" << endl;
        oss << "Variables=" << endl;
        oss << "ky" << endl;
        oss << "Euu" << endl;
        oss << "Evv" << endl;
        oss << "Eww" << endl;
    }

    string formatSpectra;
    PrepareFormatedSpectraKy(GridData, formatSpectra);

    oss << formatSpectra;

    WriteASCIIFile(file, oss.str());

    ParallelCloseFile(actkeyDumpSpectra);

    delete actkeyDumpSpectra;
}

void StatisticsHIT::PrepareFormatedSpectraKy(SpecGrid *GridData, string &formatSpectra)
{
    RDouble reynolds = GlobalDataBase::GetDoubleParaFromDB("refReNumber");

    RDouble realYL = GlobalDataBase::GetDoubleParaFromDB("YL");

    int jStartFourierIndexGlobal = GridData -> jStartFourierIndexGlobal;

    int nY = GridData -> nY;

    ostringstream oss;

    oss << setiosflags(ios::left);
    oss << setprecision(5);
    oss << setiosflags(ios::scientific);
    oss << setiosflags(ios::showpoint);

    for (int iy = jStartFourierIndexGlobal; iy <= (nY/2+1); ++ iy)
    {
        RDouble rk_y = static_cast<double>(iy - jStartFourierIndexGlobal) / realYL;
        oss << rk_y << "    ";
        oss << (*energyOfUUKy)(iy, 1) << "    ";
        oss << (*energyOfUUKy)(iy, 2) << "    ";
        oss << (*energyOfUUKy)(iy, 3) << "    ";
        oss << "\n";
    }

    oss << endl;

    formatSpectra = oss.str();
}

void Statistics::PrepareFormatedSpectraKy(SpecDiffHybGrid *GridData, string &formatSpectra)
{
    RDouble reynolds = GlobalDataBase::GetDoubleParaFromDB("refReNumber");

    RDouble realYL = GlobalDataBase::GetDoubleParaFromDB("YL");

    RDouble1D *realZ = GridData -> realZ;

    int istf = GridData -> iStartFourierIndex;

    int jStartFourierIndexGlobal = GridData -> jStartFourierIndexGlobal;

    int nY = GridData -> nY;

    ostringstream oss;

    oss << setiosflags(ios::left);
    oss << setprecision(5);
    oss << setiosflags(ios::scientific);
    oss << setiosflags(ios::showpoint);

    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        oss << "zone t = \"yplus = " << (1.0 + (*realZ)(iz)) * reynolds << "\"\n";

        for (int iy = jStartFourierIndexGlobal; iy <= (nY/2+1); ++ iy)
        {
            RDouble rk_y = static_cast<double>(iy - jStartFourierIndexGlobal) / realYL;
            oss << rk_y << "    ";
            oss << (*energyOfUUKy)(iz, iy, 1) << "    ";
            oss << (*energyOfUUKy)(iz, iy, 2) << "    ";
            oss << (*energyOfUUKy)(iz, iy, 3) << "    ";
            oss << "\n";
        }
    }

    oss << endl;

    formatSpectra = oss.str();
}

void StatisticsHIT::WriteSpectraKzVisual(SpecGrid *GridData)
{
    ActionKey *actkeyDumpSpectra = new ActionKey();
    FillActionKey(actkeyDumpSpectra, VISUALIZATION_AVERAGE_FLOW, 0);

    string spectraFileName = "results/Euu_kz";

    actkeyDumpSpectra->filename = spectraFileName;
    actkeyDumpSpectra->openmode = ios_base::out|ios_base::trunc;

    if (!actkeyDumpSpectra->IsNeedOpenFile())
    {
        return;
    }

    ParallelOpenFile(actkeyDumpSpectra);

    std::ostringstream oss;
    fstream &file = *(actkeyDumpSpectra->file);

    if (IfFileEmpty(file))
    {
        oss << "Title=\"SPECTRA\"" << endl;
        oss << "Variables=" << endl;
        oss << "kz" << endl;
        oss << "Euu" << endl;
        oss << "Evv" << endl;
        oss << "Eww" << endl;
    }

    string formatSpectra;
    PrepareFormatedSpectraKz(GridData, formatSpectra);

    oss << formatSpectra;

    WriteASCIIFile(file, oss.str());

    ParallelCloseFile(actkeyDumpSpectra);

    delete actkeyDumpSpectra;
}

void StatisticsHIT::PrepareFormatedSpectraKz(SpecGrid *GridData, string &formatSpectra)
{
    RDouble reynolds = GlobalDataBase::GetDoubleParaFromDB("refReNumber");

    RDouble realZL = GlobalDataBase::GetDoubleParaFromDB("ZL");

    int iStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;

    int nZ = GridData -> nZ;

    ostringstream oss;

    oss << setiosflags(ios::left);
    oss << setprecision(5);
    oss << setiosflags(ios::scientific);
    oss << setiosflags(ios::showpoint);

    for (int iz = iStartFourierIndexGlobal; iz <= (nZ/2+1); ++ iz)
    {
        RDouble kz = static_cast<double>(iz - iStartFourierIndexGlobal) / realZL;
        oss << kz << "    ";
        oss << (*energyOfUUKz)(iz, 1) << "    ";
        oss << (*energyOfUUKz)(iz, 2) << "    ";
        oss << (*energyOfUUKz)(iz, 3) << "    ";
        oss << "\n";
    }

    oss << endl;

    formatSpectra = oss.str();
}

/*
void Statistics::WriteStatisticVisual(SpecDiffHybGrid *GridData)
{
    RDouble1D *realZ = GridData -> realZ;
    int istf = GridData -> iStartFourierIndex;
    int x0Fall = GridData -> x0Fall;
    int xNFall = GridData -> xNFall;
    int y0Fall = GridData -> y0Fall;
    int nY = GridData -> nY;
    Range IHalf(istf, iz_half);
    Range XN(x0Fall, xNFall);
    Range YNHalf(y0Fall, nY/2+1);
    Range N9(1, 9);

    Param_SpecDiffHybSolver *parameters = GetControlParameters();
    RDouble reynolds = parameters->GetRefReNumber();
    RDouble dreynolds2 = 1.0 / reynolds / reynolds;

    fstream file;
    file.open("statistic_selected.dat", ios_base::out);
    if (!file)
    {
        TK_Exit::ExceptionExit("could not open statistic_selected.dat\n");
    }
    file << "variables = i, z, zplus, uavg, vavg, wavg, ruu, rvv, rww, epsilon, eps_ux, eps_vx, eps_wx, eps_uy, eps_vy, eps_wy, eps_uz, eps_vz, eps_wz \n";
    int wordwidth = 20;
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        file << setw(wordwidth) << iz;
        file << setiosflags(ios::right);
        file << setprecision(10);
        file << setiosflags(ios::scientific);
        file << setiosflags(ios::showpoint);
        file << setw(wordwidth) << (*realZ)(iz)
             << setw(wordwidth) << (1.0 + (*realZ)(iz)) * reynolds
             << setw(wordwidth) << (*uAverage)(iz)
             << setw(wordwidth) << (*vAverage)(iz)
             << setw(wordwidth) << (*wAverage)(iz)
             << setw(wordwidth) << (*uSquare)(iz)
             << setw(wordwidth) << (*vSquare)(iz)
             << setw(wordwidth) << (*wSquare)(iz)
             << setw(wordwidth) << (*epsilon)(iz) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 1) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 2) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 3) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 4) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 5) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 6) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 7) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 8) * dreynolds2
             << setw(wordwidth) << (*epsilonXYZ)(iz, 9) * dreynolds2
             << "\n";
    }
    file.close();
    file.clear();

    double realXL = GlobalDataBase::GetDoubleParaFromDB("XL");       
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        //! write rEuu_kx.
        string cheuufx = "./spectra/rEuu_kx_";
        stringstream strmiz;
        strmiz << iz;
        cheuufx = cheuufx + strmiz.str() + ".dat";
        file.open(cheuufx.c_str(), ios_base::out);
        if (!file)
        {
            TK_Exit::ExceptionExit("could not open rEuu_kx.dat\n");
        }
        file << "variables = kx, Euu, Evv, Eww \n";
        file << "zone t = \" yplus = " << setw(wordwidth) << (1.0 + (*realZ)(iz)) * reynolds << "\"" << "\n";
        for (int ix = x0Fall; ix <= xNFall; ++ ix)
        {
            double rk_x = static_cast<double>(ix - x0Fall) / realXL;
            file << setiosflags(ios::right);
            file << setprecision(10);
            file << setiosflags(ios::scientific);
            file << setiosflags(ios::showpoint);
            file << setw(wordwidth) << rk_x
                 << setw(wordwidth) << (*energyOfUUKx)(iz, ix, 1)
                 << setw(wordwidth) << (*energyOfUUKx)(iz, ix, 2)
                 << setw(wordwidth) << (*energyOfUUKx)(iz, ix, 3)
                 << "\n";
        }
        file.close();
        file.clear();
    }

    double realYL = GlobalDataBase::GetDoubleParaFromDB("YL");
    for (int iz = istf; iz <= iz_half; ++ iz)
    {
        // !write rEuu_ky
        string cheuufy = "./spectra/rEuu_ky_";
        stringstream strmiz;
        strmiz << iz;
        cheuufy = cheuufy +  strmiz.str() + ".dat";
        file.open(cheuufy.c_str(), ios_base::out);
        if (!file)
        {
            TK_Exit::ExceptionExit("could not open rEuu_ky.dat\n");
        }
        file << "variables = ky, Euu, Evv, Eww \n";
        file << "zone t = \"yplus = " << setw(wordwidth) << (1.0 + (*realZ)(iz)) * reynolds << "\"" << "\n";
        for (int iy = y0Fall; iy <= (nY/2+1); ++ iy)
        {
            double rk_y = static_cast<double>(iy - y0Fall) / realYL;
            file << setiosflags(ios::right);
            file << setprecision(10);
            file << setiosflags(ios::scientific);
            file << setiosflags(ios::showpoint);
            file << setw(wordwidth) << rk_y
                 << setw(wordwidth) << (*energyOfUUKy)(iz, iy, 1)
                 << setw(wordwidth) << (*energyOfUUKy)(iz, iy, 2)
                 << setw(wordwidth) << (*energyOfUUKy)(iz, iy, 3)
                 << "\n";
        }
        file.close();
        file.clear();
    }// !for iz
}
*/

void StatisticsHIT::FillActionKey(ActionKey *actkey, int action, int level)
{
    actkey->action     = action;
    actkey->solver     = 1;
    actkey->solverID   = 0;
    actkey->kind       = SOLVER_BASED;
    actkey->level      = level;
}

bool StatisticsHIT::IfStatisticExist()
{
    string outputdir = "./results";
    outputdir = GlobalDataBase::GetStrParaFromDB("outputdir");

    bool statisticExist_flag = false;

    ifstream infile(statisticsRestartFile.c_str(), ios::in);
    if (infile)
    {
        statisticExist_flag = true;
        infile.close();
        infile.clear();
    }

    return statisticExist_flag;
}

RDouble Statistics::GetTauWAverage()
{
    return this->tauWAverage;
}

void Statistics::SetTauWAverage(const RDouble &tauWAverageIn)
{
    this->tauWAverage = tauWAverageIn;
}

bool Statistics::GetIfInitStatistic()
{
    return this->ifInitStatistic;
}

void Statistics::SetIfInitStatistic(const bool &ifInitStatisticIn)
{
    this->ifInitStatistic = ifInitStatisticIn;
}

/*
LIB_EXPORT void Statistics::InitControlParameters()
{
    FreeControlParameters();
    controlParameters = new Param_SpecSolver();
    controlParameters->Init();
}

LIB_EXPORT Param_SpecSolver *Statistics::GetControlParameters() const
{
    return static_cast <Param_SpecSolver *> (controlParameters);
}
*/
}