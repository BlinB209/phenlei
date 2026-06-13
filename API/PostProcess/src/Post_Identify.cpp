#include "Post_Identify.h"
#include "PHIO.h"
#include "TK_Exit.h"
#include "TK_Log.h"
#include "GlobalDataBase.h"
#include "Constants.h"
#include "AleManager.h"
#include <cmath>

using namespace std;
namespace PHSPACE
{
    Identify::Identify()
    {
        amplitude = -PHSPACE::GetRDoubleParameterFromDataBase(0, "amplitude") * PI / 180.0;
        reduceFrequency  = PHSPACE::GetRDoubleParameterFromDataBase(0, "reduceFrequency");
        lenthReference   = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
        physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

        reduceFrequency  = reduceFrequency * lenthReference;
        physicalTimeStep = physicalTimeStep / lenthReference;

        nForceVar = 15;
        nIdentifyVar = 8;

        c0 = new RDouble[nIdentifyVar];
        c1 = new RDouble[nIdentifyVar];
        c2 = new RDouble[nIdentifyVar];

        nameList.push_back("cfx");
        nameList.push_back("cfy");
        nameList.push_back("cfz");
        nameList.push_back("cmx");
        nameList.push_back("cmy");
        nameList.push_back("cmz");
        nameList.push_back("cl");
        nameList.push_back("cd");
    }

    Identify::~Identify()
    {
        delete [] c0;    c0 = NULL;
        delete [] c1;    c1 = NULL;
        delete [] c2;    c2 = NULL; 
    }

    void Identify::Run()
    {
        ReadForce();
        ForceArrange();
        DumpHysteresisData();

        if (isIdentify)
        {
            IdentifyProcess();
            DumpResults();
        }
    }

    void Identify::ReadForce()
    {
        isIdentify = true;
        int nBlankLine = 17;
        fstream file;
        string forceFile = GlobalDataBase::GetStrParaFromDB("aircoeffile");
        file.open(forceFile.c_str(), ios_base::in);
        vector <RDouble> force(nForceVar);

        string blankLine, line;
        for (int iLine = 0; iLine < nBlankLine; ++ iLine)
        {
            getline(file, blankLine);
        }

        nDataLine = 1;
        while (!file.eof())
        {
            for (int iVar = 0; iVar < nForceVar; ++ iVar)
            {
                file >> force[iVar];
            }

            vector <RDouble> last = force;
            if (!aeroForce.empty())
            {
                last = aeroForce.back();
            }

            nwStep = max(int(force[0] - last[0]), nwStep);

            //! remove data line in aircoef.dat which is repeated!
           if (aeroForce.empty() || int(force[0] - nwStep) == int(last[0]))
            {
                aeroForce.push_back(force);
                ++ nDataLine;
            }
        }
        nDataLine -= 1;

        if (nDataLine < 3)
        {
            TK_Exit::ExceptionExit("Error: The data number is less than one cycle! \n");
        }

        RDouble period = PI;
        nLineOneCycle = static_cast<int>(period / (reduceFrequency * physicalTimeStep));
        if (nDataLine < nLineOneCycle)
        {
            TK_Exit::ExceptionExit("Error: The data number is less than one cycle! \n");
        }
        file.close();
        file.clear();
    }

    void Identify::ForceArrange()
    {
        for (int iLine = 0; iLine < nDataLine; ++ iLine)
        {
            vector <RDouble> &force = aeroForce[iLine];
            //! arrange force[1~8] by cfx,cfy,cfz,cmx,cmy,cmz,cl,cd
            RDouble cl  = force[1 ];
            RDouble cd  = force[2 ];
            RDouble cfx = force[8 ];
            RDouble cfy = force[9 ];
            RDouble cfz = force[10];
            RDouble cmx = force[11];
            RDouble cmy = force[12];
            RDouble cmz = force[13];
            force[1] = -cfx;
            force[2] =  cfy;
            force[3] = -cfz;
            force[4] = -cmx;
            force[5] =  cmy;
            force[6] = -cmz;
            force[7] =   cl;
            force[8] =   cd;
        }
    }

    void Identify::IdentifyProcess()
    {
        GetNumberOfCycle();
        for (int iVar = 1; iVar <= nIdentifyVar; ++ iVar)
        {
            getIntegral(iVar);
        }
    }

    void Identify::getIntegral(int iVar)
    {
        RDouble *f = new RDouble[nLineEnd - nLinestart + 1]();
        getIntegrand(iVar, f);

        RDouble sum1 = 0.0, sum2 = 0.0;
        int integralOrder = 4;
        if (GlobalDataBase::IsExist("derivativeFileName", PHSTRING, 1))
        {
            integralOrder = GlobalDataBase::GetIntParaFromDB("integralOrder");
        }

        if (integralOrder == 2)
        {
            GetIntegral2ndOrder(iVar, f, sum1, sum2);
        }
        else if (integralOrder == 4)
        {
            GetIntegral4thOrder(iVar, f, sum1, sum2);
        }
        else
        {
            TK_Exit::ExceptionExit("No such integral Order! \n");
        }

        sum1 = sum1 / nCycle;
        sum2 = sum2 / nCycle;

        RDouble amp = amplitude * PI;
        sum1 = 2.0 * reduceFrequency * sum1 / amp;
        sum2 = 2.0 * sum2 / amp;

        c1[iVar - 1] = sum1;
        c2[iVar - 1] = sum2;

        delete [] f;    f = NULL;
    }

    void Identify::getIntegrand(int iVar, RDouble *f)
    {
        RDouble fSum = 0.0;
        for (int iLine = static_cast<int>(nLinestart); iLine <= nLineEnd; ++ iLine)
        {
            fSum += aeroForce[iLine][iVar];
        }
        RDouble f0 = fSum / (nLineEnd - nLinestart + 1.0);

        c0[iVar - 1] = f0;

        for (int iLine = 0; iLine <= nLineEnd - nLinestart; ++ iLine)
        {
            f[iLine] = aeroForce[iLine + nLinestart][iVar] - f0;
        }
        f[0] = half * (f[0] + f[nLineEnd - nLinestart]);
        f[nLineEnd - nLinestart] = f[0];
    }

    void Identify::GetIntegral2ndOrder(int iVar, RDouble *f, RDouble &sum1, RDouble &sum2)
    {
        sum1 = 0.0;
        sum2 = 0.0;
        int nbeg = 0;
        RDouble *f1= new RDouble[nLineEnd - nLinestart + 1]();
        RDouble *f2= new RDouble[nLineEnd - nLinestart + 1]();

        RDouble ddt = nwStep * physicalTimeStep;

        for (int m = 0; m <= nLineEnd - nLinestart; m ++)
        {
            int n = int(aeroForce[m + nLinestart][0]) + nbeg;
            RDouble kt = 2.0 * n * reduceFrequency * physicalTimeStep;

            f1[m] = f[m] * sin(kt);
            f2[m] = f[m] * cos(kt);
        }

        for (int m = 1; m < nLineEnd - nLinestart; m ++)
        {
            sum1 += f1[m];
            sum2 += f2[m];
        }

        sum1 = ddt * (half * (f1[0] + f1[nLineEnd - nLinestart]) + sum1);
        sum2 = ddt * (half * (f2[0] + f2[nLineEnd - nLinestart]) + sum2);
        delete [] f1;    f1 = NULL;
        delete [] f2;    f2 = NULL;
    }

    void Identify::GetIntegral4thOrder(int iVar, RDouble *f, RDouble &sum1, RDouble &sum2)
    {
        sum1 = 0.0;
        sum2 = 0.0;
        int nbeg = 0;
        RDouble *f1= new RDouble[nLineEnd - nLinestart + 1]();
        RDouble *f2= new RDouble[nLineEnd - nLinestart + 1]();

        RDouble ddt = nwStep * physicalTimeStep;
        int nMiddle = static_cast<int>((nLineEnd - nLinestart)) / 2;

        for (int m = 0; m < nMiddle; ++ m)
        {
            int j = 2 * m + 1;
            int n = int(aeroForce[j + nLinestart][0]) + nbeg;
            RDouble kt = 2.0 * n * reduceFrequency * physicalTimeStep;

            f1[j] = 4.0 * f[j] * sin(kt);
            f2[j] = 4.0 * f[j] * cos(kt);

            j = j + 1;
            n = int(aeroForce[j + nLinestart][0]) + nbeg;
            kt = 2.0 * n * reduceFrequency * physicalTimeStep;

            f1[j] = 2.0 * f[j] * sin(kt);
            f2[j] = 2.0 * f[j] * cos(kt);
        }

        for (int m = 0; m <= nLineEnd - nLinestart; m += static_cast<int>(nLineEnd - nLinestart))
        {
            int n = int(aeroForce[m + nLinestart][0]) + nbeg;
            RDouble kt = 2.0 * n * reduceFrequency * physicalTimeStep;

            f1[m] = f[m] * sin(kt);
            f2[m] = f[m] * cos(kt);
        }

        for (int m = 0; m <= nLineEnd - nLinestart; m ++)
        {
            sum1 += f1[m];
            sum2 += f2[m];
        }

        sum1 = ddt/3.0*sum1;
        sum2 = ddt/3.0*sum2;
        delete [] f1;    f1 = NULL;
        delete [] f2;    f2 = NULL;
    }

    void Identify::GetNumberOfCycle()
    {
        int totalCycle = nDataLine / nLineOneCycle;
        int startCycle = 0;
        if (totalCycle > 2)
        {
            nCycle = totalCycle - 1;
            startCycle = 1;
        }
        else if (totalCycle < 1)
        {
            nCycle = 0;
            TK_Exit::ExceptionExit("Error: The data number is less than one cycle! \n");
        }
        else
        {
            nCycle = 1;
        }
    
        nLinestart = nDataLine - nCycle * nLineOneCycle;
        nLineEnd   = nDataLine - 1;
    }

    void Identify::DumpResults()
    {
        string derivativeFileName = "results/identify.dat";
        if (GlobalDataBase::IsExist("derivativeFileName", PHSTRING, 1))
        {
            derivativeFileName = GlobalDataBase::GetStrParaFromDB("derivativeFileName");
        }
        int integralOrder = 4;
        if (GlobalDataBase::IsExist("derivativeFileName", PHSTRING, 1))
        {
            integralOrder = GlobalDataBase::GetIntParaFromDB("integralOrder");
        }
        string forceFile = GlobalDataBase::GetStrParaFromDB("aircoeffile");
        amplitude = PHSPACE::GetRDoubleParameterFromDataBase(0, "amplitude");
        reduceFrequency = PHSPACE::GetRDoubleParameterFromDataBase(0, "reduceFrequency");
        lenthReference  = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
        physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

        fstream file;
        file.open(derivativeFileName.c_str(), ios_base::out | ios_base::trunc);

        int wordWidth = 16;
        int floatNumber = 4;
        ostringstream oss;
        oss << setprecision(floatNumber);
        oss << setiosflags(ios::right);
        oss << setiosflags(ios::scientific);

        //if (PHSPACE::IfFileEmpty(file))
        {
            oss << "Static and damping derivatives identified result." << "\n";
            oss << "\n";
            oss << setw(wordWidth) << "Data file:" << "\t" << forceFile.c_str() << "\n";

            oss << setw(wordWidth) << "Frequency:"  << "\t" << reduceFrequency << "\n";
            oss << setw(wordWidth) << "Amplitude:"  << "\t" << amplitude << "\n";
            oss << setw(wordWidth) << "Timestep:"   << "\t" << physicalTimeStep << "\n";
            oss << setw(wordWidth) << "LenthRef:"   << "\t" << lenthReference << "\n";
            oss << setw(wordWidth) << "Cycles used:" << "\t" << nCycle << "\n";
            oss << setw(wordWidth) << "Integral Order:" << "\t" << integralOrder << "\n";
            oss << "\n";

            oss << setw(wordWidth) << "Parameter";
            oss << setw(wordWidth) << "Average";
            oss << setw(wordWidth) << "dynamic_zero";
            oss << setw(wordWidth) << "dynamic_1st";
            oss << "\n";
            for (int iVar = 0; iVar < nIdentifyVar; ++ iVar)
            {
                oss << setw(wordWidth) << nameList[iVar];
                oss << setw(wordWidth) << c0[iVar];
                oss << setw(wordWidth) << c1[iVar];
                oss << setw(wordWidth) << c2[iVar];
                oss << "\n";
            }

            PHSPACE::WriteASCIIFile(file, oss.str());

            file.close();
            file.clear();

        }
    }

    void Identify::DumpHysteresisData()
    {
        string hysteresisFileName = "results/force_beta.plt";
        if (GlobalDataBase::IsExist("hysteresisFileName", PHSTRING, 1))
        {
            hysteresisFileName = GlobalDataBase::GetStrParaFromDB("hysteresisFileName");
        }
        RDouble amp = PHSPACE::GetRDoubleParameterFromDataBase(0, "amplitude");
        RDouble frequency = PHSPACE::GetRDoubleParameterFromDataBase(0, "reduceFrequency");
        RDouble timeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

        fstream file;
        file.open(hysteresisFileName.c_str(), ios_base::out | ios_base::trunc);

        int wordWidth = 16;
        int floatNumber = 6;
        ostringstream oss;
        oss << setprecision(floatNumber);
        oss << setiosflags(ios::right);
        oss << setiosflags(ios::scientific);

        oss << "variables = beta, cfx, cfy, cfz, cmx, cmy, cmz" << "\n";
        for (int iLine = 0; iLine < nDataLine; ++ iLine)
        {
            RDouble beta = -amp * sin(2.0 * aeroForce[iLine][0] * frequency * timeStep);
            oss << setw(wordWidth) << beta;
            for (int iVar = 1; iVar <= 6; ++ iVar)
            {
                oss << setw(wordWidth) << aeroForce[iLine][iVar];
            }
            oss << "\n";
        }

        PHSPACE::WriteASCIIFile(file, oss.str());

        file.close();
        file.clear();
    }
}
