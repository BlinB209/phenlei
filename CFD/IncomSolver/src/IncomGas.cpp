#include "IncomGas.h"
#include "GlobalDataBase.h"
#include "PHIO.h"
#include <algorithm>

using namespace std;

namespace PHSPACE
{
namespace GAS_SPACE
{
    IncomGas::IncomGas()
    {
        
        componentLib = {"AIR", "H2", "H2_l", "NH3_l", "CH4", "C2H4", "C2H6", "C3H8", "NG", "LNG", "H2S", "Cl2", "NH3"};
        SetPropertyLib();
    }

    IncomGas::~IncomGas()
    {
        int muType = GlobalDataBase::GetIntParaFromDB("muType");
        if (muType == muIdealGasMixingLaw)
        {
            DelPointer2(phiij);
        }
        int massdiffType = GlobalDataBase::GetIntParaFromDB("massdiffType");
        if (massdiffType == massDiffDiluteApprox)
        {
            delete[] dilute;
        }
    }

    void IncomGas::InitCommonParameter()
    {
        PerfectGas::InitCommonParameter();

           //! Operating Pressure
           //! Operating Pressure
        
    }

    void IncomGas::InitParameterForMu()
    {
        RDouble refT = GlobalDataBase::GetDoubleParaFromDB("refT");
        int muType = GlobalDataBase::GetIntParaFromDB("muType");
        if (muType == muIdealGasMixingLaw)
        {
            ComputePhiIJ();
        }
    }

    void IncomGas::InitParameterForK()
    {
        int kType = GlobalDataBase::GetIntParaFromDB("kType");
        if (kType == kIdealGasMixingLaw)
        {
            ComputePhiIJ();
        }
    }

    void IncomGas::InitParameterForMassdiff()
    {
        int massdiffType = GlobalDataBase::GetIntParaFromDB("massdiffType");
        int numOfGas = GlobalDataBase::GetIntParaFromDB("numOfGas");
        if (massdiffType == massDiffConstantDiluteApprox)
        {
            dilute = new RDouble[1];
            if (GlobalDataBase::IsExist("dilute", PHINT, 1))
            {
                GlobalDataBase::GetData("dilute", dilute, PHDOUBLE, 1);
            }
            else
            {
                dilute[0] = 2.88e-5;
            }
        }
        else if (massdiffType == massDiffDiluteApprox)
        {
            dilute = new RDouble[numOfGas];
            GlobalDataBase::GetData("dilute", dilute, PHDOUBLE, numOfGas);
        }
        else if (massdiffType == massDiffMulticomponent)
        {

        }
    }

    void IncomGas::SetGasName(string gasName)
    {
        int numOfGas = GlobalDataBase::GetIntParaFromDB("numOfGas");
        gasNameList = new string[numOfGas];
        GlobalDataBase::GetData(gasName, gasNameList, PHSTRING, numOfGas);
    }

    void IncomGas::SetGasMolarMass(string gasMolarMass)
    {
        int numOfGas = GlobalDataBase::GetIntParaFromDB("numOfGas");
        gasMolarMassList = new RDouble[numOfGas];
        for (int iSp = 0; iSp < numOfGas; ++iSp)
        {
            auto iter = find(componentLib.begin(), componentLib.end(), gasNameList[iSp]);
            if (iter != componentLib.end())
            {
                gasMolarMassList[iSp] = propertyLib[gasNameList[iSp]].molarMass;
            }
        }
    }

    void IncomGas::SetGasRho(string gasRho)
    {
        int numOfGas = GlobalDataBase::GetIntParaFromDB("numOfGas");
        gasRhoList = new RDouble[numOfGas];
        for (int iSp = 0; iSp < numOfGas; ++iSp)
        {
            auto iter = find(componentLib.begin(), componentLib.end(), gasNameList[iSp]);
            if (iter != componentLib.end())
            {
                gasRhoList[iSp] = propertyLib[gasNameList[iSp]].rho;
            }

        }
    }

    void IncomGas::SetGasMu(string gasMu)
    {
        int numOfGas = GlobalDataBase::GetIntParaFromDB("numOfGas");
        gasMuList = new RDouble[numOfGas];
        for (int iSp = 0; iSp < numOfGas; ++iSp)
        {
            auto iter = find(componentLib.begin(), componentLib.end(), gasNameList[iSp]);
            if (iter != componentLib.end())
            {
                gasMuList[iSp] = propertyLib[gasNameList[iSp]].mu;
            }
        }
    }

    void IncomGas::SetGasCp(string gasCp)
    {
        int numOfGas = GlobalDataBase::GetIntParaFromDB("numOfGas");
        gasCpList = new RDouble[numOfGas];
        for (int iSp = 0; iSp < numOfGas; ++iSp)
        {
            auto iter = find(componentLib.begin(), componentLib.end(), gasNameList[iSp]);
            if (iter != componentLib.end())
            {
                gasCpList[iSp] = propertyLib[gasNameList[iSp]].cp;
            }
        }
    }

    void IncomGas::SetGasK(string gasK)
    {
        int numOfGas = GlobalDataBase::GetIntParaFromDB("numOfGas");
        gasKList = new RDouble[numOfGas];
        for (int iSp = 0; iSp < numOfGas; ++iSp)
        {
            auto iter = find(componentLib.begin(), componentLib.end(), gasNameList[iSp]);
            if (iter != componentLib.end())
            {
                gasKList[iSp] = propertyLib[gasNameList[iSp]].k;
            }
        }
    }

    void IncomGas::UpdateRg(UnstructGrid* grid, RDouble* Rg)
    {
        int nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();
        int numOfGas = GlobalDataBase::GetIntParaFromDB("numOfGas");
        string gasName = gasNameList[0];
        RDouble* massFracOfGas = reinterpret_cast<RDouble*>(grid->GetDataPtr(gasName));
        RDouble gasR = 8.314;
        if (massFracOfGas == NULL)
        {
            for (int iCell = 0; iCell < nTotal; ++iCell)
            {
                Rg[iCell] = gasR / gasMolarMassList[0];
            }
        }
        else
        {
            SetField(Rg, 0.0 , nTotal);
            for (int iGas = 0; iGas < numOfGas; ++iGas)
            {
                massFracOfGas = reinterpret_cast<RDouble*>(grid->GetDataPtr(gasNameList[iGas]));
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    Rg[iCell] += massFracOfGas[iCell] * gasR / gasMolarMassList[iGas];
                }
            }
        }        
    }

    void IncomGas::UpdateRho(UnstructGrid* grid, RDouble* rho)
    {
        RDouble refT = GlobalDataBase::GetDoubleParaFromDB("refT");
        RDouble refP = GlobalDataBase::GetDoubleParaFromDB("refP");
        int rhoType = GlobalDataBase::GetIntParaFromDB("rhoType");
        switch (rhoType)
        {
        case rhoConstant:
        {
            break;
        }

        case rhoIncompressibleIdealGas:
        {
            int nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();
            RDouble* T = reinterpret_cast<RDouble*>(grid->GetDataPtr("T"));
            RDouble* Rg = reinterpret_cast<RDouble*>(grid->GetDataPtr("Rg"));
            RDouble* cRho = reinterpret_cast<RDouble*>(grid->GetDataPtr("cRho"));
            
            if (T == NULL)
            {
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    cRho[iCell] = 1 / (Rg[iCell] * refT);
                    rho[iCell] = refP * cRho[iCell];
                }
            }
            else
            {
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    cRho[iCell] = 1 / (Rg[iCell] * T[iCell]);
                    rho[iCell] = refP * cRho[iCell];
                }
            }        

            break;
        }

        case rhoIdealGas:
        {
            int nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();
            RDouble* T = reinterpret_cast<RDouble*>(grid->GetDataPtr("T"));
            RDouble* p = reinterpret_cast<RDouble*>(grid->GetDataPtr("P"));
            RDouble* Rg = reinterpret_cast<RDouble*>(grid->GetDataPtr("Rg"));
            RDouble* cRho = reinterpret_cast<RDouble*>(grid->GetDataPtr("cRho"));

            if (T == NULL)
            {
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    cRho[iCell] = 1 / (Rg[iCell] * refT);
                    rho[iCell] = (refP + p[iCell]) * cRho[iCell];
                }
            }
            else
            {
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    cRho[iCell] = 1 / (Rg[iCell] * T[iCell]);
                    rho[iCell] = (refP + p[iCell]) * cRho[iCell];
                }
            }        

            break;
        }

        case rhoVolumnWeightedMixingLaw:
        {
            int nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();
            string gasName = gasNameList[0];
            RDouble* massFracOfGas = reinterpret_cast<RDouble*>(grid->GetDataPtr(gasName));
            int numOfGas = GlobalDataBase::GetIntParaFromDB("numOfGas");
            for (int iCell = 0; iCell < nTotal; ++iCell)
            {
                rho[iCell] = massFracOfGas[iCell] / gasRhoList[0];
            }

            for (int iGas = 1; iGas < numOfGas; ++iGas)
            {
                gasName = gasNameList[iGas];
                massFracOfGas = reinterpret_cast<RDouble*>(grid->GetDataPtr(gasName));

                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    rho[iCell] += massFracOfGas[iCell] / gasRhoList[iGas];
                }
            }

            for (int iCell = 0; iCell < nTotal; ++iCell)
            {
                rho[iCell] = 1 / rho[iCell];
            }

            break;
        }
        }
    }

    void IncomGas::UpdateMu(UnstructGrid* grid, RDouble* mu)
    {
        RDouble refT = GlobalDataBase::GetDoubleParaFromDB("refT");
        int muType = GlobalDataBase::GetIntParaFromDB("muType");
        switch (muType)
        {
        case muConstant:
        {
            break;
        }

        case muSutherLand:
        {
            int nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();
            RDouble* T = reinterpret_cast<RDouble*>(grid->GetDataPtr("T"));
            int coeOfSutherland = GlobalDataBase::GetIntParaFromDB("numOfCoeSutherland");

            if (coeOfSutherland == 2)
            {
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    mu[iCell] = 1.458e-6 * pow(T[iCell], 3/2) / (T[iCell] + 110.4);
                }
            }
            else if (coeOfSutherland == 3)
            {
                RDouble referenceValue = refT + 110.56;
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    mu[iCell] = 1.716e-5 * pow(T[iCell] / refT, 3 / 2) * referenceValue / (T[iCell] + 110.56);
                }
            }
            
            break;
        }

        case muIdealGasMixingLaw:
        {
            int numOfGas = GlobalDataBase::GetIntParaFromDB("numOfGas");
            int nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();
            RDouble* sumOfAmount = NewPointer<RDouble>(nTotal);
            RDouble** amountOfSubstance = NewPointer2<RDouble>(numOfGas, nTotal);
            RDouble** amountFracOfSubstance = NewPointer2<RDouble>(numOfGas, nTotal);
            ComputeAmountOfSubstance(grid, amountOfSubstance, amountFracOfSubstance);

            for (int iCell = 0; iCell < nTotal; ++iCell)
            {
                sumOfAmount[iCell] = amountFracOfSubstance[0][iCell] * phiij[0][0];
            }

            for (int iGas = 1; iGas < numOfGas; ++iGas)
            {
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    sumOfAmount[iCell] += amountFracOfSubstance[iGas][iCell] * phiij[0][iGas];
                }
            }

            for (int iCell = 0; iCell < nTotal; ++iCell)
            {
                mu[iCell] = amountFracOfSubstance[0][iCell] * gasMuList[0] / sumOfAmount[iCell];
            }

            for (int iGas = 1; iGas < numOfGas; ++iGas)
            {
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    mu[iCell] += amountFracOfSubstance[iGas][iCell] * gasMuList[iGas] / sumOfAmount[iCell];
                }
            }

            DelPointer(sumOfAmount);
            DelPointer2(amountOfSubstance);
            DelPointer2(amountFracOfSubstance);
            break;
        }

        case muMassWeightedMixingLaw:
        {
            int numOfGas = GlobalDataBase::GetIntParaFromDB("numOfGas");
            int nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();
            string gasName = gasNameList[0];
            RDouble* massFracOfGas = reinterpret_cast<RDouble*>(grid->GetDataPtr(gasName));

            for (int iCell = 0; iCell < nTotal; ++iCell)
            {
                mu[iCell] = massFracOfGas[iCell] * gasMuList[0];
            }

            for (int iGas = 1; iGas < numOfGas; ++iGas)
            {
                gasName = gasNameList[iGas];
                massFracOfGas = reinterpret_cast<RDouble*>(grid->GetDataPtr(gasName));
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    mu[iCell] += massFracOfGas[iCell] * gasMuList[iGas];
                }
            }

            break;
        }

        case muPiecewiseLinear:
        {

            break;
        }

        case muPolynomial:
        {
            break;
        }

        case muPowerLaw:
        {
            break;
        }
        }
    }

    void IncomGas::UpdateCp(UnstructGrid* grid, RDouble* cp)
    {
        int cpType = GlobalDataBase::GetIntParaFromDB("cpType");
        switch (cpType)
        {
        case cpConstant:
        {
            break;
        }

        case cpMixing:
        {
            int numOfGas = GlobalDataBase::GetIntParaFromDB("numOfGas");
            int nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();
            string gasName = gasNameList[0];
            RDouble* massFracOfGas = reinterpret_cast<RDouble*>(grid->GetDataPtr(gasName));

            for (int iCell = 0; iCell < nTotal; ++iCell)
            {
                cp[iCell] = massFracOfGas[iCell] * gasCpList[0];
            }

            for (int iGas = 1; iGas < numOfGas; ++iGas)
            {
                gasName = gasNameList[iGas];
                massFracOfGas = reinterpret_cast<RDouble*>(grid->GetDataPtr(gasName));
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    cp[iCell] += massFracOfGas[iCell] * gasCpList[iGas];
                }
            }

            break;
        }

        case cpPiecewiseLinear:
        {
            break;
        }

        case cpPiecewisePolynomial:
        {
            break;
        }

        case cpPolynomial:
        {
            break;
        }

        }
    }

    void IncomGas::UpdateK(UnstructGrid* grid, RDouble* k)
    {
        int kType = GlobalDataBase::GetIntParaFromDB("kType");
        switch (kType)
        {
        case kConstant:
        {
            break;
        }

        case kMassWeightedMixingLaw:
        {
            int numOfGas = GlobalDataBase::GetIntParaFromDB("numOfGas");
            int nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();
            string gasName = gasNameList[0];
            RDouble* massFracOfGas = reinterpret_cast<RDouble*>(grid->GetDataPtr(gasName));

            for (int iCell = 0; iCell < nTotal; ++iCell)
            {
                k[iCell] = massFracOfGas[iCell] * gasKList[0];
            }

            for (int iGas = 1; iGas < numOfGas; ++iGas)
            {
                gasName = gasNameList[iGas];
                massFracOfGas = reinterpret_cast<RDouble*>(grid->GetDataPtr(gasName));
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    k[iCell] += massFracOfGas[iCell] * gasKList[iGas];
                }
            }

            break;
        }

        case kIdealGasMixingLaw:
        {
            int numOfGas = GlobalDataBase::GetIntParaFromDB("numOfGas");
            int nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();
            RDouble* sumOfAmount = NewPointer<RDouble>(nTotal);
            RDouble** amountOfSubstance = NewPointer2<RDouble>(numOfGas, nTotal);
            RDouble** amountFracOfSubstance = NewPointer2<RDouble>(numOfGas, nTotal);
            ComputeAmountOfSubstance(grid, amountOfSubstance, amountFracOfSubstance);

            for (int iCell = 0; iCell < nTotal; ++iCell)
            {
                sumOfAmount[iCell] = amountFracOfSubstance[0][iCell] * phiij[0][0];
            }

            for (int iGas = 1; iGas < numOfGas; ++iGas)
            {
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    sumOfAmount[iCell] += amountFracOfSubstance[iGas][iCell] * phiij[0][iGas];
                }
            }

            for (int iCell = 0; iCell < nTotal; ++iCell)
            {
                k[iCell] = amountFracOfSubstance[0][iCell] * gasKList[0] / sumOfAmount[iCell];
            }

            for (int iGas = 1; iGas < numOfGas; ++iGas)
            {
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    k[iCell] += amountFracOfSubstance[iGas][iCell] * gasKList[iGas] / sumOfAmount[iCell];
                }
            }

            DelPointer(sumOfAmount);
            DelPointer2(amountOfSubstance);
            DelPointer2(amountFracOfSubstance);

            break;
        }

        case kPiecewiseLinear:
        {
            break;
        }

        case kPiecewisePolynomial:
        {
            break;
        }

        case kPolynomial:
        {
            break;
        }

        }
    }

    void IncomGas::UpdateMassDiff(UnstructGrid* grid, int iGas, RDouble* massdiff)
    {
        RDouble sct = 0.7;
        int isSolveTurb = GlobalDataBase::GetIntParaFromDB("isSolveTurb");
        int massdiffType = GlobalDataBase::GetIntParaFromDB("massdiffType");
        switch (massdiffType)
        {
        case massDiffConstantDiluteApprox:
        {
            int nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();
            RDouble* rho = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));

            if (isSolveTurb == 0)
            {
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    massdiff[iCell] = rho[iCell] * dilute[0];
                }
            }
            else
            {
                RDouble* vist = reinterpret_cast<RDouble*>(grid->GetDataPtr("vist"));
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    massdiff[iCell] = rho[iCell] * dilute[0];
                    massdiff[iCell] += vist[iCell] / sct;
                }
            }
            break;
        }

        case massDiffDiluteApprox:
        {
            int nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();
            RDouble* rho = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));

            if (isSolveTurb == 0)
            {
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    massdiff[iCell] = rho[iCell] * dilute[iGas];
                }
            }
            else
            {
                RDouble* vist = reinterpret_cast<RDouble*>(grid->GetDataPtr("vist"));
                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    massdiff[iCell] = rho[iCell] * dilute[iGas];
                    massdiff[iCell] += vist[iCell] / sct;
                }
            }
            break;
        }

        case massDiffMulticomponent:
        {
            break;
        }

        case massDiffUnityLewisNumber:
        {
            int nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();
            RDouble* k = reinterpret_cast<RDouble*>(grid->GetDataPtr("k"));
            RDouble* cp = reinterpret_cast<RDouble*>(grid->GetDataPtr("cp"));
            RDouble* rho = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
            for (int iCell = 0; iCell < nTotal; ++iCell)
            {
                massdiff[iCell] = k[iCell] / (rho[iCell] * cp[iCell]);
            }
            break;
        }
        }
    }

    void IncomGas::ComputePhiIJ()
    {
        int numOfGas = GlobalDataBase::GetIntParaFromDB("numOfGas");
        phiij = NewPointer2<RDouble>(numOfGas, numOfGas);

        for (int iGas = 0; iGas < numOfGas; ++iGas)
        {
            for (int jGas = 0; jGas < numOfGas; ++jGas)
            {
                RDouble molarMassI = gasMolarMassList[iGas];
                RDouble molarMassJ = gasMolarMassList[jGas];
                RDouble muI = gasMuList[iGas];
                RDouble muJ = gasMuList[jGas];

                RDouble molarMassRatio = molarMassJ / molarMassI;
                RDouble muRatio = muI / muJ;

                phiij[iGas][jGas] = pow(1 + pow(muRatio, 1 / 2) * pow(molarMassRatio, 1 / 4), 2) / pow((1 + molarMassRatio) * 8, 1 / 2);
            }
        }
    }

    void IncomGas::ComputeAmountOfSubstance(UnstructGrid* grid, RDouble** amountOfSubstance, RDouble** amountFracOfSubstance)
    {
        int numOfGas = GlobalDataBase::GetIntParaFromDB("numOfGas");
        int nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();
        RDouble* vol = grid->GetCellVolume();
        RDouble* rho = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
        RDouble* totalAmountOfSubstance = new RDouble[nTotal];
        for (int iCell = 0; iCell < nTotal; ++iCell)
        {
            totalAmountOfSubstance[iCell] = 0.0;
        }

        if (numOfGas == 1)
        {
            string gasName = gasNameList[0];
            RDouble molarMass = gasMolarMassList[0];

            for (int iCell = 0; iCell < nTotal; ++iCell)
            {
                amountOfSubstance[0][iCell] = rho[iCell] * vol[iCell] / molarMass;
                totalAmountOfSubstance[iCell] += amountOfSubstance[0][iCell];
            }
        }
        else
        {
            for (int iGas = 0; iGas < numOfGas; ++iGas)
            {
                string gasName = gasNameList[iGas];
                RDouble molarMass = gasMolarMassList[iGas];
                RDouble* massFracOfGas = reinterpret_cast<RDouble*>(grid->GetDataPtr(gasName));

                for (int iCell = 0; iCell < nTotal; ++iCell)
                {
                    amountOfSubstance[iGas][iCell] = massFracOfGas[iCell] * rho[iCell] * vol[iCell] / molarMass;
                    totalAmountOfSubstance[iCell] += amountOfSubstance[iGas][iCell];
                }
            }
        }        

        for (int iGas = 0; iGas < numOfGas; ++iGas)
        {
            for (int iCell = 0; iCell < nTotal; ++iCell)
            {
                amountFracOfSubstance[iGas][iCell] = amountOfSubstance[iGas][iCell] / totalAmountOfSubstance[iCell];
            }
        }

        delete[] totalAmountOfSubstance;
    }

    void IncomGas::SetPropertyLib()
    {
        propertySp p_AIR;
        propertySp p_H2;
        propertySp p_H2_l;
        propertySp p_NH3_l;
        propertySp p_CH4;
        propertySp p_C2H4;
        propertySp p_C3H8;
        propertySp p_C2H6;
        propertySp p_NG;
        propertySp p_LNG;
        propertySp p_H2S;
        propertySp p_Cl2;
        propertySp p_NH3;        

        p_AIR.rho = 1.225;
        p_AIR.mu  = 1.7894e-5;
        p_AIR.k   = 0.0242;
        p_AIR.cp  = 1006.43;
        p_AIR.molarMass  = 0.028;

        p_H2.rho = 0.08189;
        p_H2.mu  = 8.411e-06;
        p_H2.k   = 0.1672;
        p_H2.cp  = 14283;
        p_H2.molarMass  = 0.00201594;

        p_H2_l.rho = 70.85;
        p_H2_l.mu  = 1.332e-5;
        p_H2_l.k   = 0.10382;
        p_H2_l.cp  = 9772.2;
        p_H2_l.molarMass = 0.00201594;

        p_NH3_l.rho = 610;
        p_NH3_l.mu  = 0.000152;
        p_NH3_l.k   = 0.493;
        p_NH3_l.cp  = 4758;
        p_NH3_l.molarMass = 0.01703061;

        p_CH4.rho = 0.6679;
        p_CH4.mu  = 1.087e-05;
        p_CH4.k   = 0.0332;
        p_CH4.cp  = 2222;
        p_CH4.molarMass = 0.01604303;

        p_C2H4.rho = 1.137;
        p_C2H4.mu  = 1.03e-05;
        p_C2H4.k   = 0.0214;
        p_C2H4.cp  = 2233;
        p_C2H4.molarMass = 0.02805418;

        p_C3H8.rho = 1.91;
        p_C3H8.mu  = 7.95e-06;
        p_C3H8.k   = 0.0177;
        p_C3H8.cp  = 1549;
        p_C3H8.molarMass = 0.04409;


        p_C2H6.rho = 1.263;
        p_C2H6.mu  = 9.29e-6;
        p_C2H6.k   = 0.0207;
        p_C2H6.cp  = 1731;
        p_C2H6.molarMass = 0.03;

        p_NG.rho = 0.753;
        p_NG.mu  = 1.07e-5;
        p_NG.k   = 0.0318;
        p_NG.cp  = 2188.9;
        p_NG.molarMass = 0.018046;

        p_LNG.rho = 434.5566;
        p_LNG.mu  = 0.14331;
        p_LNG.k = 0.019549;
        p_LNG.cp  = 3425.3;
        p_LNG.molarMass = 0.018046;

        p_H2S.rho = 1.46;
        p_H2S.mu  = 1.2e-05;
        p_H2S.k   = 0.0134;
        p_H2S.cp  = 1170;
        p_H2S.molarMass = 0.03407994;

        p_Cl2.rho = 2.95;
        p_Cl2.mu  = 1.33e-05;
        p_Cl2.k   = 0.00882;
        p_Cl2.cp  = 650;
        p_Cl2.molarMass = 0.070906;

        p_NH3.rho = 0.6894;
        p_NH3.mu  = 1.015e-05;
        p_NH3.k   = 0.0247;
        p_NH3.cp  = 2158;
        p_NH3.molarMass = 0.01703061;

        propertyLib = { {"AIR", p_AIR}, {"H2", p_H2}, {"H2_l", p_H2_l}, {"NH3_l", p_NH3_l}, {"CH4", p_CH4}, 
            {"C2H4",p_C2H4}, {"C3H8", p_C3H8},  {"C2H6", p_C2H6}, {"NG", p_NG}, {"LNG", p_LNG}, 
            {"H2S", p_H2S}, {"Cl2", p_Cl2}, {"NH3", p_NH3}};
    }

}
}