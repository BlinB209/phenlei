#include "Spectrum.h"
#include "Param_Spectrum.h"
#include "SpecDiffHybGrid.h"
#include "PHHeader.h"
#include "GlobalDataBase.h"
#include "TK_Exit.h"
#include "TK_Parse.h"

using namespace std;

namespace PHSPACE
{
Spectrum::Spectrum()
{
}

Spectrum::~Spectrum()
{
    delete []k_CBC;
    delete []eK_CBC;
    delete []logK_CBC;
    delete []logEK_CBC;
    delete []d2logEKdK2_CBC;
}

void Spectrum::InitSpectrumData()
{
    //Param_Spectrum *parameters = GetControlParameters();
    spectrumType = GlobalDataBase::GetIntParaFromDB("spectrumType");
    nPoints = 20;

    AllocSpectrumData(nPoints);

    if(spectrumType == 3)
    {
        InitCBCSpectrum();
    }
    else
    {
        TK_Exit::UnexpectedVarValue("spectrumType", spectrumType);
    }
}

void Spectrum::AllocSpectrumData(int nPoints)
{
    k_CBC = new RDouble [nPoints];
    eK_CBC = new RDouble [nPoints];
    logK_CBC = new RDouble [nPoints];
    logEK_CBC = new RDouble [nPoints];
    d2logEKdK2_CBC = new RDouble [nPoints];

    PHSPACE::SetField(k_CBC, 0.0, nPoints);
    PHSPACE::SetField(eK_CBC, 0.0, nPoints);
    PHSPACE::SetField(logK_CBC, 0.0, nPoints);
    PHSPACE::SetField(logEK_CBC, 0.0, nPoints);
    PHSPACE::SetField(d2logEKdK2_CBC, 0.0, nPoints);
}

void Spectrum::InitCBCSpectrum()
{
    const RDouble pi = 2.0 * acos(0.0);
    const RDouble m = 0.0508;
    const RDouble u0 = 10.0;
    const RDouble re = 34000.0;
    const RDouble lRef = 10.0 * m / (2.0 * pi);
    const RDouble tRef = 64.0 * m / u0;
    const RDouble uRef = lRef / tRef;

    RDouble dlogEKdK1_CBC = 0.0;
    RDouble dlogEKdKn_CBC = 0.0;

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    //int procIDWithKx0Ky0 = GridData -> procIDWithKx0Ky0;
    if(myid == 0)
    {
        ReadCBCSpectrum();

        RDouble kScale = 100.0 * lRef;
        RDouble eKScale = 1.0e-6 * tRef * tRef / (lRef * lRef * lRef);

        for(int n = 0; n < nPoints; ++n)
        {
            k_CBC[n] = k_CBC[n] * kScale;
            eK_CBC[n] = eK_CBC[n] * eKScale;

            logK_CBC[n] = log(k_CBC[n]);
            logEK_CBC[n] = log(eK_CBC[n]);
        }

        WriteCBCSpectrum();

        dlogEKdK1_CBC = (logEK_CBC[1] - logEK_CBC[0]) / (logK_CBC[1] - logK_CBC[0]);
        dlogEKdKn_CBC = (logEK_CBC[nPoints - 1] - logEK_CBC[nPoints - 2]) / (logK_CBC[nPoints - 1] - logK_CBC[nPoints - 2]);

        Spline(k_CBC, eK_CBC, nPoints, dlogEKdK1_CBC, dlogEKdKn_CBC, d2logEKdK2_CBC);
    }
    
    MPI_Bcast(&nPoints, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
    MPI_Bcast(logK_CBC, nPoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(logEK_CBC, nPoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(d2logEKdK2_CBC, nPoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dlogEKdK1_CBC, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dlogEKdKn_CBC, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void Spectrum::ReadCBCSpectrum()
{
    string filename = "ek.dat";
    filename = GlobalDataBase::GetStrParaFromDB("ekCBCfilename");
    fstream file;
    file.open(filename.c_str(),ios_base::in);
    if ( !file )
    {
        TK_Exit::FileOpenErrorExit(filename);
    }

    string line,word;
    string separator  = " =\t\r\n#$,;\"'";

    for(int n = 0; n < nPoints; ++n)
    {
        ReadNewLine(file, line);

        line = FindNextWord(line, word, separator);
        from_string< double >( k_CBC[n], word, std::dec );

        line = FindNextWord(line, word, separator);
        from_string< double >( eK_CBC[n], word, std::dec );
    }
    file.close();
    file.clear();
}

void Spectrum::WriteCBCSpectrum()
{
    string filename = "ek_rescale.dat";

    fstream file;
    file.open(filename.c_str(),ios_base::out);
    if ( !file )
    {
        TK_Exit::FileOpenErrorExit(filename);
    }
    file << "variables = k, Ek \n";
    for (int n = 0; n < nPoints; ++ n)
    {
        file << setiosflags(ios::right);
        file << setprecision(10);
        file << setiosflags(ios::scientific);
        file << setiosflags(ios::showpoint);
        file << setw(20) << k_CBC[n]
             << setw(20) << eK_CBC[n]
             << "\n";
    }
    file.close();
    file.clear();
}

void Spectrum::Spline(RDouble *x, RDouble *y, int n, RDouble yp1, RDouble ypn, RDouble *y2)
{
    RDouble *u = new RDouble [nPoints];
    PHSPACE::SetField(u, 0.0, nPoints);

    if(yp1 > LARGE)
    {
        y2[0] = 0.0;
        u[0] = 0.0;
    }
    else
    {
        y2[0] = -0.5;
        u[0] = 3.0 / (x[1] - x[0]) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
    }

    for(int i = 1; i < n-1; ++i)
    {
        RDouble sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
        RDouble p = sig * y2[i-1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        u[i] = (6.0 * ((y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1]) / (x[i] - x[i-1])) / (x[i+1] - x[i-1]) - sig * u[i-1]) / p;
    }

    RDouble qn = 0.0;
    RDouble un = 0.0;
    if(ypn < LARGE)
    {
        qn = 0.5;
        un = 3.0 / (x[n-1] - x[n-2]) * (ypn - (y[n-1] - y[n-2]) / (x[n-1] - x[n-2]));
    }

    y2[n-1] = (un - qn * u[n-2]) / (qn * y2[n-2] + 1.0);

    for(int i = n-2; i >= 0; --i)
    {
        y2[i] = y2[i] * y2[i+1] + u[i];
    }

    delete []u;
}

}
