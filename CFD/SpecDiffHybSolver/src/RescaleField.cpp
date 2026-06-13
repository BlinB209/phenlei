#include "RescaleField.h"
#include "SpecGrid.h"
//#include "SpecSolver.h"
//#include "SpecDiffHybSolver.h"
#include "GlobalDataBase.h"
#include "PHMpi.h"
#include "TK_Exit.h"
#include "PHHeader.h"
#include "TK_Parse.h"

using namespace std;

namespace PHSPACE
{
RescaleFieldHIT::RescaleFieldHIT()
{
}

RescaleFieldHIT::~RescaleFieldHIT()
{
    delete rescaleStatistic; rescaleStatistic = NULL;
    delete rescaleStatisticLocal; rescaleStatisticLocal = NULL;
}

RescaleField::RescaleField()
{
}

RescaleField::~RescaleField()
{
    delete uAverage; uAverage = NULL;
    delete uFluctuation; uFluctuation = NULL;
    delete vFluctuation; vFluctuation = NULL;
    delete wFluctuation; wFluctuation = NULL;
    delete uAverageOfReference; uAverageOfReference = NULL;
    delete uFluctuationOfReference; uFluctuationOfReference = NULL;
    delete vFluctuationOfReference; vFluctuationOfReference = NULL;
    delete wFluctuationOfReference; wFluctuationOfReference = NULL;
    //delete rescaleStatistic; rescaleStatistic = NULL;
    //delete rescaleStatisticLocal; rescaleStatisticLocal = NULL;

    //FreeControlParameters();
}

void RescaleFieldHIT::InitRescaleData()
{
    //InitControlParameters();

    AllocRescaleData();
}

void RescaleField::InitRescaleData(SpecDiffHybGrid *GridData)
{
    //InitControlParameters();

    AllocRescaleData(GridData);

    SetUavg();

    GetUVWMagnitudeOfReference(GridData);
}

void RescaleFieldHIT::AllocRescaleData()
{
    rescaleStatistic = new RDouble1D(Range(1, 4), fortranArray);
    rescaleStatisticLocal = new RDouble1D(Range(1, 4), fortranArray);
    (*rescaleStatistic)(Range(1, 4)) = 0.0;
    (*rescaleStatisticLocal)(Range(1, 4)) = 0.0;
}

void RescaleField::AllocRescaleData( SpecDiffHybGrid *GridData )
{
    int ist = GridData -> iStartFourierIndex;
    int ied = GridData -> iEndFourierIndex;
    Range I(ist, ied);

    uAverage = new RDouble1D(I, fortranArray);
    uFluctuation = new RDouble1D(I, fortranArray);
    vFluctuation = new RDouble1D(I, fortranArray);
    wFluctuation = new RDouble1D(I, fortranArray);
    (*uAverage)(I) = 0.0;
    (*uFluctuation)(I) = 0.0;
    (*vFluctuation)(I) = 0.0;
    (*wFluctuation)(I) = 0.0;

    uAverageOfReference = new RDouble1D(I, fortranArray);
    uFluctuationOfReference = new RDouble1D(I, fortranArray);
    vFluctuationOfReference = new RDouble1D(I, fortranArray);
    wFluctuationOfReference = new RDouble1D(I, fortranArray);
    (*uAverageOfReference)(I) = 0.0;
    (*uFluctuationOfReference)(I) = 0.0;
    (*vFluctuationOfReference)(I) = 0.0;
    (*wFluctuationOfReference)(I) = 0.0;

    rescaleStatistic = new RDouble2D(I, Range(1, 4), fortranArray);
    rescaleStatisticLocal = new RDouble2D(I, Range(1, 4), fortranArray);
    (*rescaleStatistic)(I, Range(1, 4)) = 0.0;
    (*rescaleStatisticLocal)(I, Range(1, 4)) = 0.0;
}

void RescaleField::SetUavg()
{
    RDouble reynolds = GlobalDataBase::GetDoubleParaFromDB("refReNumber");

    uavgDivideByUtau = 2.960975610 + 2.439024390 * log(reynolds);
}

void RescaleField::GetUVWMagnitudeOfReference( SpecDiffHybGrid *GridData )
{
    int ist = GridData -> iStartFourierIndex;
    int ied = GridData -> iEndFourierIndex;
    Range I(ist, ied);

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();
    if( myid == procIDWithKx0Ky0 )
    {
        InputUVWMagnitudeOfReference( GridData );
        (*rescaleStatistic)(I, 1) = (*uAverageOfReference)(I);
        (*rescaleStatistic)(I, 2) = (*uFluctuationOfReference)(I);
        (*rescaleStatistic)(I, 3) = (*vFluctuationOfReference)(I);
        (*rescaleStatistic)(I, 4) = (*wFluctuationOfReference)(I);
    }

    int n = 4 * (ied - ist + 1);
#ifdef PH_PARALLEL
    MPI_Bcast(&(*rescaleStatistic)(ist, 1) , n, MPI_DOUBLE, server, MPI_COMM_WORLD);
#endif
    (*uAverageOfReference)(I) = (*rescaleStatistic)(I, 1);
    (*uFluctuationOfReference)(I) = (*rescaleStatistic)(I, 2);
    (*vFluctuationOfReference)(I) = (*rescaleStatistic)(I, 3);
    (*wFluctuationOfReference)(I) = (*rescaleStatistic)(I, 4);
}

void RescaleField::GetUVWMagnitude(SpecDiffHybGrid *GridData, Complex3D *complexU, Complex3D *complexV, Complex3D *complexW)
{
    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int n = 4 * (iedf - istf + 1);
    Range I(istf, iedf);

    GetUVWLocal(GridData, complexU, complexV, complexW);

    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();
    (*rescaleStatistic)(I, Range(1, 4)) = (*rescaleStatisticLocal)(I, Range(1, 4));

#ifdef PH_PARALLEL
    MPI_Reduce( &(*rescaleStatisticLocal)(istf, 1), &(*rescaleStatistic)(istf, 1), n, MPI_DOUBLE, MPI_SUM, procIDWithKx0Ky0, MPI_COMM_WORLD );
    MPI_Bcast(&(*rescaleStatistic)(istf, 1), n, MPI_DOUBLE, procIDWithKx0Ky0, MPI_COMM_WORLD);
#endif

    for(int iz = istf; iz <= iedf; ++iz)
    {
        (*uAverage)(iz) = (*rescaleStatistic)(iz, 1);
        (*uFluctuation)(iz) = sqrt((*rescaleStatistic)(iz, 2));
        (*vFluctuation)(iz) = sqrt((*rescaleStatistic)(iz, 3));
        (*wFluctuation)(iz) = sqrt((*rescaleStatistic)(iz, 4));
    }  
}

void RescaleField::InputUVWMagnitudeOfReference(SpecDiffHybGrid *GridData)
{
    int ist = GridData -> iStartFourierIndex;
    int ied = GridData -> iEndFourierIndex;
    int innz = ied - ist + 1;

    string filename = "DNS.dat";
    filename = GlobalDataBase::GetStrParaFromDB("DNSDataFile");
    fstream file;
    file.open(filename.c_str(),ios_base::in);
    if ( !file )
    {
        TK_Exit::FileOpenErrorExit(filename);
    }

    string line,word;
    string separator  = " =\t\r\n#$,;\"'";

    for(int n = 1; n <= 3; ++n)
    {
        ReadNewLine(file, line);
    }
    int innz_dns;
    ReadNewLine(file, line);
    line = FindNextWord(line,word,separator);
    from_string< int >( innz_dns, word, std::dec );

    RDouble1D *yDNS = new RDouble1D( Range(1, innz_dns), fortranArray);
    RDouble1D *uavgDNS = new RDouble1D( Range(1, innz_dns), fortranArray);
    RDouble1D *ufDNS = new RDouble1D( Range(1, innz_dns), fortranArray);
    RDouble1D *vfDNS = new RDouble1D( Range(1, innz_dns), fortranArray);
    RDouble1D *wfDNS = new RDouble1D( Range(1, innz_dns), fortranArray);
    RDouble1D *yy = new RDouble1D( Range(ist, ied), fortranArray);
    double yplustmp;
    for( int n = 1; n <= innz_dns; ++n )
    {
        ReadNewLine(file, line);

        line = FindNextWord(line, word, separator);
        from_string< double >( (*yDNS)(n), word, std::dec );

        line = FindNextWord(line, word, separator);
        from_string< double >( yplustmp, word, std::dec );

        line = FindNextWord(line, word, separator);
        from_string< double >( (*uavgDNS)(n), word, std::dec );

        line = FindNextWord(line, word, separator);
        from_string< double >( (*ufDNS)(n), word, std::dec );

        line = FindNextWord(line, word, separator);
        from_string< double >( (*wfDNS)(n), word, std::dec );

        line = FindNextWord(line, word, separator);
        from_string< double >( (*vfDNS)(n), word, std::dec );
    }
    file.close();
    file.clear();

    RDouble1D *rZ = GridData -> realZ;
    (*yy)( Range(ist, ied) ) = 1.0 + (*rZ)( Range(ist, ied) );
    Interpolation( &(*yy)(ist), &(*uAverageOfReference)(ist), innz/2+1, &(*yDNS)(1), &(*uavgDNS)(1), innz_dns );
    Interpolation( &(*yy)(ist), &(*uFluctuationOfReference)(ist), innz/2+1, &(*yDNS)(1), &(*ufDNS)(1), innz_dns );
    Interpolation( &(*yy)(ist), &(*vFluctuationOfReference)(ist), innz/2+1, &(*yDNS)(1), &(*vfDNS)(1), innz_dns );
    Interpolation( &(*yy)(ist), &(*wFluctuationOfReference)(ist), innz/2+1, &(*yDNS)(1), &(*wfDNS)(1), innz_dns );

    int iz_half = ist + innz/2;
    for(int iz = ist; iz <= iz_half; ++iz)
    {
        int izup = innz - iz + ist;
        (*uAverageOfReference)(izup) = (*uAverageOfReference)(iz);
        (*uFluctuationOfReference)(izup) = (*uFluctuationOfReference)(iz);
        (*vFluctuationOfReference)(izup) = (*vFluctuationOfReference)(iz);
        (*wFluctuationOfReference)(izup) = (*wFluctuationOfReference)(iz);
    }

    int wordwidth = 20;
    filename = "DNSinterp.dat";
    file.open(filename.c_str(), ios_base::out);
    if ( !file )
    {
        TK_Exit::ExceptionExit("could not open DNSinterp.dat \n");
    }
    file << "variables = y, Uav+, u+, v+, w+ \n";
    for(int iz = ist; iz <= ied; ++iz)
    {
        file << setiosflags(ios::left);
        file << setiosflags(ios::scientific);
        file << setprecision(10);
        file << setw(wordwidth) << (*rZ)(iz) 
            << setw(wordwidth) << (*uAverageOfReference)(iz) 
            << setw(wordwidth) << (*uFluctuationOfReference)(iz) 
            << setw(wordwidth) << (*vFluctuationOfReference)(iz) 
            << setw(wordwidth) << (*wFluctuationOfReference)(iz) 
            << "\n";
    }

    delete yDNS; yDNS = NULL;
    delete uavgDNS; uavgDNS = NULL;
    delete ufDNS; ufDNS = NULL;
    delete vfDNS; vfDNS = NULL;
    delete wfDNS; wfDNS = NULL;
    delete yy; yy = NULL;
}

void RescaleField::FixUAverage(SpecDiffHybGrid *GridData, Complex3D *complexU)
{ 
    RDouble1D *rZ = GridData -> realZ;
    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;
    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();
    Range I(istf, iedf);

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if(myid != procIDWithKx0Ky0)
    {
        return;
    }

    RDouble1D *uav = new RDouble1D(I, fortranArray);
    for(int iz = istf; iz <= iedf; ++iz)
    {
        (*uav)(iz) = real((*complexU)(iz, jstf, kstf));
    }

    double ub = 0.0;
    for(int iz = istf + 1; iz <= iedf - 1; ++iz)
    {
        ub = ub + 0.5 * ((*rZ)(iz+1) - (*rZ)(iz-1)) * (*uav)(iz);
    }
    ub = ub / ((*rZ)(iedf) - (*rZ)(istf));

    (*complexU)(I, jstf, kstf) = (*complexU)(I, jstf, kstf) * uavgDivideByUtau / ub;

    delete uav; uav = NULL;
}

void RescaleField::GetUVWLocal(SpecDiffHybGrid *GridData, Complex3D *complexU, Complex3D *complexV, Complex3D *complexW)
{
    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;
    //int x0Fall = GridData -> x0Fall;
    //int y0Fall = GridData -> y0Fall;
    int kStartFourierIndexGlobal = GridData -> kStartFourierIndexGlobal;
    int jStartFourierIndexGlobal = GridData -> jStartFourierIndexGlobal;

    (*rescaleStatisticLocal)(Range(istf, iedf), Range(1, 4)) = 0.0;

    for(int iz = istf; iz <= iedf; ++iz)
    {
        PHComplex sumU = PHComplex(0.0, 0.0);
        PHComplex sumV = PHComplex(0.0, 0.0);
        PHComplex sumW = PHComplex(0.0, 0.0);
        if(kstf == kStartFourierIndexGlobal)
        {
            int iy0 = jstf;
            if(jstf == jStartFourierIndexGlobal)
            {
                iy0 = jstf + 1;
                (*rescaleStatisticLocal)(iz, 1) = real( (*complexU)(iz, jstf, kstf) );
            }

            for(int iy = iy0; iy <= jedf; ++iy)
            {
                sumU = sumU + (*complexU)(iz, iy, kstf) * conj((*complexU)(iz, iy, kstf));
                sumV = sumV + (*complexV)(iz, iy, kstf) * conj((*complexV)(iz, iy, kstf));
                sumW = sumW + (*complexW)(iz, iy, kstf) * conj((*complexW)(iz, iy, kstf));
            }
            for(int ix = kstf + 1; ix <= kedf; ++ix)
            {
                for(int iy = jstf; iy <= jedf; ++iy)
                {
                    sumU = sumU + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix));
                    sumV = sumV + 2.0 * (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix));
                    sumW = sumW + 2.0 * (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                }
            }
        }
        else
        {
            for(int ix = kstf; ix <= kedf; ++ix)
            {
                for(int iy = jstf; iy <= jedf; ++iy)
                {
                    sumU = sumU + 2.0 * (*complexU)(iz, iy, ix) * conj((*complexU)(iz, iy, ix));
                    sumV = sumV + 2.0 * (*complexV)(iz, iy, ix) * conj((*complexV)(iz, iy, ix));
                    sumW = sumW + 2.0 * (*complexW)(iz, iy, ix) * conj((*complexW)(iz, iy, ix));
                }
            }
        }
        (*rescaleStatisticLocal)(iz, 2) = real(sumU);
        (*rescaleStatisticLocal)(iz, 3) = real(sumV);
        (*rescaleStatisticLocal)(iz, 4) = real(sumW);
    }
}

void RescaleField::RescaleUAverage(SpecDiffHybGrid *GridData, Complex3D *complexU)
{      
    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    int procIDWithKx0Ky0 = GridData -> GetProcIDWithKx0Ky0();

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    if(myid != procIDWithKx0Ky0)
    {
        return;
    }
    else
    {
        for(int iz = istf; iz <= iedf; ++iz)
        {
            (*complexU)(iz, jstf, kstf) = PHComplex((*uAverageOfReference)(iz), 0.0);
        }
    }  
}

void RescaleField::RescaleFluctuation(SpecDiffHybGrid *GridData, Complex3D *complexU, Complex3D *complexV, Complex3D *complexW)
{
    int istf = GridData -> iStartFourierIndex;
    int iedf = GridData -> iEndFourierIndex;
    int jstf = GridData -> jStartFourierIndex;
    int jedf = GridData -> jEndFourierIndex;
    int kstf = GridData -> kStartFourierIndex;
    int kedf = GridData -> kEndFourierIndex;

    for(int iz = istf + 1; iz <= iedf - 1; ++iz)
    {
        (*complexU)(iz, Range(jstf+1, jedf), Range(kstf, kedf)) = (*complexU)(iz, Range(jstf+1, jedf), Range(kstf, kedf)) / (*uFluctuation)(iz) * (*uFluctuationOfReference)(iz);
        (*complexV)(iz, Range(jstf+1, jedf), Range(kstf, kedf)) = (*complexV)(iz, Range(jstf+1, jedf), Range(kstf, kedf)) / (*vFluctuation)(iz) * (*vFluctuationOfReference)(iz);
        (*complexW)(iz, Range(jstf+1, jedf), Range(kstf, kedf)) = (*complexW)(iz, Range(jstf+1, jedf), Range(kstf, kedf)) / (*wFluctuation)(iz) * (*wFluctuationOfReference)(iz);

        (*complexU)(iz, jstf, Range(kstf+1, kedf)) = (*complexU)(iz, jstf, Range(kstf+1, kedf)) / (*uFluctuation)(iz) * (*uFluctuationOfReference)(iz);
        (*complexV)(iz, jstf, Range(kstf+1, kedf)) = (*complexV)(iz, jstf, Range(kstf+1, kedf)) / (*vFluctuation)(iz) * (*vFluctuationOfReference)(iz);
        (*complexW)(iz, jstf, Range(kstf+1, kedf)) = (*complexW)(iz, jstf, Range(kstf+1, kedf)) / (*wFluctuation)(iz) * (*wFluctuationOfReference)(iz);
    }
}

void RescaleField::Interpolation( RDouble *yout, RDouble *uout, int n, RDouble *yin, RDouble *uin, int n0 )
{
    RDouble1D y( yout, Range(1, n), neverDeleteData, fortranArray );
    RDouble1D u( uout, Range(1, n), neverDeleteData, fortranArray );
    RDouble1D y0( yin, Range(1, n0), neverDeleteData, fortranArray );
    RDouble1D u0( uin, Range(1, n0), neverDeleteData, fortranArray );

    u( Range(1, n) ) = 0.0;
    u(1) = 0.0;
    int i0 = 1;
    double yl = 0.0;
    double yr = 0.0;
    double ul = 0.0;
    double ur = 0.0;

    for(int i = 2; i <= n; ++i)
    {
        while( (y0(i0) < y(i)) && (i0 < n0) )
        {
            i0 = i0 + 1;
        }
        if( i0 == n0 )
        {
            u( Range(i, n) ) = u0(n0);
            break;
        }
        else if( i0 == 1 )
        {
            yl = 0.0;
            yr = y0(i0);
            ul = 0.0;
            ur = u0(i0);
        }
        else
        {
            yl = y0(i0 - 1);
            yr = y0(i0);
            ul = u0(i0 - 1);
            ur = u0(i0);
        }
        u(i) = ul + (y(i) - yl) * (ur - ul) / (yr - yl);
    }
}

/*
LIB_EXPORT void RescaleField::InitControlParameters()
{
    FreeControlParameters();
    controlParameters = new Param_SpecSolver();
    controlParameters->Init();
}

LIB_EXPORT Param_SpecSolver *RescaleField::GetControlParameters() const
{
    return static_cast <Param_SpecSolver *> (controlParameters);
}
*/
}

