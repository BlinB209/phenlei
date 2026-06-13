#include "MixedGas.h"
#include "Constants.h"
#include "PHMatrix.h"
#include "TK_Exit.h"
#include "TK_Parse.h"
#include "Glb_Dimension.h"
#include "Geo_SimpleBC.h"
#include "TK_Log.h"
#include "PHIO.h"
#pragma warning(disable:6386)
#pragma warning(disable:4100)
#pragma warning(disable:4701)
#pragma warning(disable:26451)

namespace PHSPACE
{
namespace GAS_SPACE
{

ThermodynamicFunction::ThermodynamicFunction()
{
    number_of_polynomial_coefficient = 0;
    polynomial_coefficient = 0;
    number_of_temperature_interval = 0;
    enthalpyFitCoefficient = 0;
    polynomial_coefficientT2 = 0;
}

ThermodynamicFunction::~ThermodynamicFunction()
{
    delete [] enthalpyFitCoefficient;
    DelPointer2(polynomial_coefficient);
    DelPointer2(polynomial_coefficientT2);
}

void ThermodynamicFunction::Init(int number_of_temperature_interval, int number_of_polynomial_coefficient)
{
    this->number_of_temperature_interval   = number_of_temperature_interval;
    this->number_of_polynomial_coefficient = number_of_polynomial_coefficient;
    polynomial_coefficient = NewPointer2<RDouble>(number_of_temperature_interval, number_of_polynomial_coefficient);
    enthalpyFitCoefficient = new RDouble[number_of_polynomial_coefficient];
    polynomial_coefficientT2 = NewPointer2<RDouble>(number_of_temperature_interval + 1, number_of_polynomial_coefficient);
}

void ThermodynamicFunction::ReadPolynomialCoefficient(fstream &file)
{
    string line, word;
    string separator  = " =\t\r\n#$,;\"'";

    for (int i_interval = 0; i_interval < number_of_temperature_interval; ++ i_interval)
    {
        getline(file, line);

        RDouble *coef = this->GetPolynomialCoefficient(i_interval);

        for (int icoef = 0; icoef < number_of_polynomial_coefficient; ++ icoef)
        {
            line = FindNextWord(line, word, separator);
            from_string< RDouble >(coef[icoef], word, std::dec);
        }
    }
}

void ThermodynamicFunction::CreatePolynomialCoefficient(RDouble **curveFitData, RDouble **curveFitDataT2, RDouble *enthalpyFitData)
{
    for (int m = 0; m < number_of_temperature_interval; ++ m)
    {
        RDouble *coef = this->GetPolynomialCoefficient(m);

        for (int n = 0; n < number_of_polynomial_coefficient; ++ n)
        {
            coef[n] = curveFitData[m][n];
        }
    }
    //The curve fit data for the two-temperature model.
    for (int m = 0; m < number_of_temperature_interval + 1; ++ m)
    {
        RDouble *coef = this->GetPolynomialCoefficientT2(m);

        for (int n = 0; n < number_of_polynomial_coefficient; ++ n)
        {
            coef[n] = curveFitDataT2[m][n];
        }
    }

    for (int n = 0; n < number_of_polynomial_coefficient; ++ n)
    {
        enthalpyFitCoefficient[n] = enthalpyFitData[n];
    }
}

void ThermodynamicFunction::ReadPolynomialCoefficient(DataContainer *cdata)
{
    for (int i_interval = 0; i_interval < number_of_temperature_interval; ++ i_interval)
    {
        RDouble *coef = this->GetPolynomialCoefficient(i_interval);

        for (int icoef = 0; icoef < number_of_polynomial_coefficient; ++ icoef)
        {
            cdata->Read(&coef[icoef], sizeof(RDouble));
        }
    }
    //The curve fit data for the two-temperature model.
    for (int i_interval = 0; i_interval < number_of_temperature_interval + 1; ++ i_interval)
    {
        RDouble *coef = this->GetPolynomialCoefficientT2(i_interval);

        for (int icoef = 0; icoef < number_of_polynomial_coefficient; ++ icoef)
        {
            cdata->Read(&coef[icoef], sizeof(RDouble));
        }
    }

    for (int icoef = 0; icoef < number_of_polynomial_coefficient; ++ icoef)
    {
        cdata->Read(&enthalpyFitCoefficient[icoef], sizeof(RDouble));
    }
}

void ThermodynamicFunction::WritePolynomialCoefficient(DataContainer *cdata)
{
    for (int i_interval = 0; i_interval < number_of_temperature_interval; ++ i_interval)
    {
        RDouble *coef = this->GetPolynomialCoefficient(i_interval);

        for (int icoef = 0; icoef < number_of_polynomial_coefficient; ++ icoef)
        {
            cdata->Append(&coef[icoef], sizeof(RDouble));
        }
    }
    //The curve fit data for the two-temperature model.
    for (int i_interval = 0; i_interval < number_of_temperature_interval + 1; ++ i_interval)
    {
        RDouble *coef = this->GetPolynomialCoefficientT2(i_interval);

        for (int icoef = 0; icoef < number_of_polynomial_coefficient; ++ icoef)
        {
            cdata->Append(&coef[icoef], sizeof(RDouble));
        }
    }

    for (int icoef = 0; icoef < number_of_polynomial_coefficient; ++ icoef)
    {
        cdata->Append(&enthalpyFitCoefficient[icoef], sizeof(RDouble));
    }
}

ThermodynamicManager::ThermodynamicManager()
{
    thermoDynamicFunction = 0;
    temperatureRange = 0;
    temperatureStart = 0;
}

ThermodynamicManager::~ThermodynamicManager()
{
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        delete thermoDynamicFunction[ispecies];
    }
    delete [] thermoDynamicFunction;
    delete [] temperatureRange;
    delete [] temperatureStart;
}

void ThermodynamicManager::Init(int numberOfSpecies)
{
    this->numberOfSpecies = numberOfSpecies;
    thermoDynamicFunction = new ThermodynamicFunction * [numberOfSpecies];
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        thermoDynamicFunction[ispecies] = new ThermodynamicFunction();
    }
}

void ThermodynamicManager::Read(fstream &file)
{
    string line, word;
    string separator = " =\t\r\n#$,;\"'";

    //! Read the number of thermodynamic temperature polynomial coefficients of species.
    SkipLines(file, 3);
    getline(file, line);

    //! numberofpolynomialcoefficient is the number of polynomial coefficients(a1-a7,or a1-a6).

    line = FindNextWord(line, word, separator);
    from_string< int >(numberOfTemperatureInterval, word, std::dec);

    line = FindNextWord(line, word, separator);
    from_string< int >(numberofpolynomialcoefficient, word, std::dec);

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        thermoDynamicFunction[ispecies]->Init(numberOfTemperatureInterval, numberofpolynomialcoefficient);
    }

    temperatureRange = new RDouble[numberOfTemperatureInterval+1];

    //! Read the thermodynamic temperature polynomial coefficients,temperature range and so on of species.
    SkipLines(file, 4);
    getline(file, line);

    for (int m = 0; m < numberOfTemperatureInterval + 1; ++ m)
    {
        line = FindNextWord(line, word, separator);
        from_string< RDouble >(temperatureRange[m], word, std::dec);
    }

    SkipLines(file, 2);

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        SkipLines(file, 3);
        thermoDynamicFunction[ispecies]->ReadPolynomialCoefficient(file);
    }
}

void ThermodynamicManager::GenerateData(vector<string> &namelist)
{
    ChemkinFitData myChemkinData;
    ChemkinFitData *currentData = new ChemkinFitData();

    myChemkinData.GenerateFullChemkinData();
    myChemkinData.CreateChemkinFitData(namelist, currentData);

    //! To obtain the basic information of interpolation.
    int nSpecies, nInterval, nCoefNum;
    currentData->GetInterpolationInfo(nSpecies, nInterval, nCoefNum);
    numberOfSpecies = nSpecies;
    numberOfTemperatureInterval = nInterval;
    numberofpolynomialcoefficient = nCoefNum;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        thermoDynamicFunction[ispecies]->Init(numberOfTemperatureInterval, numberofpolynomialcoefficient);
    }

    //! To obtain the temperature range of interpolation.
    RDouble *temperatures = currentData->GetTemperatureRange();
    temperatureRange = new RDouble[numberOfTemperatureInterval + 1];
    for (int m = 0; m < numberOfTemperatureInterval + 1; ++ m)
    {
        temperatureRange[m] = temperatures[m];
    }

    //! To obtain the curve fitting data of interpolation of each species.
    RDouble ***curveFittingData = currentData->GetCurveFittingData();
    RDouble ***curveFittingDataT2 = currentData->GetCurveFittingDataT2();
    RDouble **curveFitData = NewPointer2 <RDouble> (numberOfTemperatureInterval, numberofpolynomialcoefficient);
    RDouble **curveFitDataT2 = NewPointer2 <RDouble> (numberOfTemperatureInterval + 1, numberofpolynomialcoefficient);
    RDouble *tmpFitData = new RDouble[numberofpolynomialcoefficient];
    RDouble **enthalpyFitData = currentData->GetEnthalpyFitData();
    RDouble *benckmarkT2 = currentData->GetStartTemperature();
    temperatureStart = new RDouble [numberOfSpecies];

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        temperatureStart[ispecies] = benckmarkT2[ispecies];

        for (int m = 0; m < numberOfTemperatureInterval; ++ m)
        {
            for (int n = 0; n < numberofpolynomialcoefficient; ++ n)
            {
                curveFitData[m][n] = curveFittingData[ispecies][m][n];
            }
        }

        for (int n = 0; n < numberofpolynomialcoefficient; ++ n)
        {
            tmpFitData[n] = enthalpyFitData[ispecies][n];
        }

        for (int m = 0; m < numberOfTemperatureInterval + 1; ++ m)
        {
            for (int n = 0; n < numberofpolynomialcoefficient; ++ n)
            {
                curveFitDataT2[m][n] = curveFittingDataT2[ispecies][m][n];
            }
        }
        thermoDynamicFunction[ispecies]->CreatePolynomialCoefficient(curveFitData, curveFitDataT2, tmpFitData);
    }

    delete [] tmpFitData;
    DelPointer2(curveFitData);
    DelPointer2(curveFitDataT2);
    delete currentData;
}

void ThermodynamicManager::Read(DataContainer *cdata)
{
    cdata->Read(&numberOfTemperatureInterval, sizeof(int));
    cdata->Read(&numberofpolynomialcoefficient, sizeof(int));

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        thermoDynamicFunction[ispecies]->Init(numberOfTemperatureInterval, numberofpolynomialcoefficient);
    }

    temperatureStart = new RDouble [numberOfSpecies];
    cdata->Read(temperatureStart, numberOfSpecies * sizeof(RDouble));

    temperatureRange = new RDouble[numberOfTemperatureInterval + 1];
    cdata->Read(temperatureRange, (numberOfTemperatureInterval + 1) * sizeof(RDouble));

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        thermoDynamicFunction[ispecies]->ReadPolynomialCoefficient(cdata);
    }
}

void ThermodynamicManager::Write(DataContainer *cdata)
{
    cdata->Append(&numberOfTemperatureInterval, sizeof(int));
    cdata->Append(&numberofpolynomialcoefficient, sizeof(int));
    cdata->Append(temperatureStart, numberOfSpecies * sizeof(RDouble));
    cdata->Append(temperatureRange, (numberOfTemperatureInterval + 1) * sizeof(RDouble));
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        thermoDynamicFunction[ispecies]->WritePolynomialCoefficient(cdata);
    }
}

RDouble * ThermodynamicManager::GetPolynomialCoefficient(int iSpecies, int iTemperature)
{
    return thermoDynamicFunction[iSpecies]->GetPolynomialCoefficient(iTemperature);
}

RDouble * ThermodynamicManager::GetPolynomialCoefficientT2(int iSpecies, int iTemperature)
{
    return thermoDynamicFunction[iSpecies]->GetPolynomialCoefficientT2(iTemperature);
}

RDouble * ThermodynamicManager::GetEnthalpyFitCoefficient(int iSpecies)
{
    return thermoDynamicFunction[iSpecies]->GetEnthalpyFitCoefficient();
}

void ThermodynamicManager::GetTemperatureRangeIndex(const RDouble &temperature_dimensional, int &indexOfTemperatureRange)
{
    //! Avoid abnormal situations.
    indexOfTemperatureRange = 0;
    if (temperatureRange[0] >= temperature_dimensional)
    {
        indexOfTemperatureRange = 0;
    }
    else if (temperatureRange[numberOfTemperatureInterval] <= temperature_dimensional)
    {
        indexOfTemperatureRange = numberOfTemperatureInterval - 1;
    }
    else
    {
        for (int i = 0; i < numberOfTemperatureInterval; ++ i)
        {
            if ((temperatureRange[i] < temperature_dimensional) && (temperature_dimensional <= temperatureRange[i+1]))
            {
                indexOfTemperatureRange = i;
                return;
            }
        }
    }
    return;
}

RDouble ThermodynamicManager::GetTemperatureMin()
{
    return temperatureRange[0];
}

RDouble ThermodynamicManager::GetTemperatureMax()
{
    return temperatureRange[numberOfTemperatureInterval];
}

int ThermodynamicManager::GetNumberofPolynomialCoefficient()
{
    return this->numberofpolynomialcoefficient;
}

RDouble ThermodynamicManager::GetBenchmarkTemperature(int iSpecies)
{
    return temperatureStart[iSpecies];
}

ReactionRate::ReactionRate()
{
    numberOfReaction = 0;
    afr = 0;
    bfr = 0;
    cfr = 0;

    abr = 0;
    bbr = 0;
    cbr = 0;
}

ReactionRate::~ReactionRate()
{
    delete [] afr;
    delete [] bfr;
    delete [] cfr;

    delete [] abr;
    delete [] bbr;
    delete [] cbr;
}

void ReactionRate::Init(int numberOfReaction)
{
    this->numberOfReaction = numberOfReaction;

    afr = new RDouble [numberOfReaction];
    bfr = new RDouble [numberOfReaction];
    cfr = new RDouble [numberOfReaction];

    abr = new RDouble [numberOfReaction];
    bbr = new RDouble [numberOfReaction];
    cbr = new RDouble [numberOfReaction];
}

void ReactionRate::Read(fstream &file)
{
    string line, word;
    string separator = " =\t\r\n#$,;\"'";

    if (numberOfReaction <= 0) return;

    //! Read the forward and reverse reaction rate coefficients of each chemical reaction.
    SkipLines(file, 3);
    for (int irection = 0; irection < numberOfReaction; ++ irection)
    {
        getline(file, line);
        //! Read irtmp.
        line = FindNextWord(line, word, separator);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (afr[irection], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (bfr[irection], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (cfr[irection], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (abr[irection], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (bbr[irection], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (cbr[irection], word, std::dec);
    }
}

void ReactionRate::Create(RDouble **reationCoef)
{
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        afr[r] = reationCoef[r][0];
        bfr[r] = reationCoef[r][1];
        cfr[r] = reationCoef[r][2];

        abr[r] = reationCoef[r][3];
        bbr[r] = reationCoef[r][4];
        cbr[r] = reationCoef[r][5];
    }
}

void ReactionRate::Read(DataContainer *cdata)
{
    cdata->Read(afr, numberOfReaction * sizeof(RDouble));
    cdata->Read(bfr, numberOfReaction * sizeof(RDouble));
    cdata->Read(cfr, numberOfReaction * sizeof(RDouble));
                                        
    cdata->Read(abr, numberOfReaction * sizeof(RDouble));
    cdata->Read(bbr, numberOfReaction * sizeof(RDouble));
    cdata->Read(cbr, numberOfReaction * sizeof(RDouble));
}

void ReactionRate::Write(DataContainer *cdata)
{
    cdata->Append(afr, numberOfReaction * sizeof(RDouble));
    cdata->Append(bfr, numberOfReaction * sizeof(RDouble));
    cdata->Append(cfr, numberOfReaction * sizeof(RDouble));

    cdata->Append(abr, numberOfReaction * sizeof(RDouble));
    cdata->Append(bbr, numberOfReaction * sizeof(RDouble));
    cdata->Append(cbr, numberOfReaction * sizeof(RDouble));
}

StoichiometricEquation::StoichiometricEquation()
{
    cvp = 0;
    cvn = 0;
    cpz = 0;
    numberOfReaction = 0;
    numberOfSpecies = 0;

    reactionType = 0;
    isCollisionType = 0;
}

StoichiometricEquation::~StoichiometricEquation()
{
    DelPointer2(cvp);
    DelPointer2(cvn);
    DelPointer2(cpz);

    delete [] reactionType;
    delete [] isCollisionType;
}

void StoichiometricEquation::Init(int numberOfReaction, int numberOfSpecies)
{
    this->numberOfReaction = numberOfReaction;
    this->numberOfSpecies  = numberOfSpecies ;

    reactionType = new int[numberOfReaction];
    isCollisionType = new bool[numberOfReaction];

    cvp = NewPointer2<int> (numberOfReaction, numberOfSpecies);
    cvn = NewPointer2<int> (numberOfReaction, numberOfSpecies);
    cpz = NewPointer2<RDouble> (numberOfReaction, numberOfSpecies);
}

void StoichiometricEquation::Read(fstream &file)
{
    string line, word;
    string separator = " =\t\r\n#$,;\"'";

    if (numberOfReaction <= 0) return;

    SkipLines(file, 4);

    //! Read chemical reaction equation coefficients of species.
    for (int irection = 0; irection < numberOfReaction; ++ irection)
    {
        getline(file, line);
        //! Read irtmp.
        line = FindNextWord(line, word, separator);

        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            line = FindNextWord(line, word, separator);
            from_string< int >(cvp[irection][ispecies], word, std::dec);
        }
    }

    SkipLines(file, 4);
    //! Read chemical reaction equation coefficients of species.
    for (int irection = 0; irection < numberOfReaction; ++ irection)
    {
        getline(file, line);
        //! Read irtmp.
        line = FindNextWord(line, word, separator);

        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            line = FindNextWord(line, word, separator);
            from_string< int >(cvn[irection][ispecies], word, std::dec);
        }
    }

    SkipLines(file, 4);
    //! Read chemical reaction equation coefficients of each collision body.
    for (int irection = 0; irection < numberOfReaction; ++ irection)
    {
        getline(file, line);
        //! Read irtmp.
        line = FindNextWord(line, word, separator);

        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            line = FindNextWord(line, word, separator);
            from_string< RDouble >(cpz[irection][ispecies], word, std::dec);
        }
    }

    //Read the types of chemical reactions.
    SkipLines(file, 3);
    getline(file, line);
    for (int irection = 0; irection < numberOfReaction; ++ irection)
    {
        line = FindNextWord(line, word, separator);
        from_string<int>(reactionType[irection], word, std::dec);

        isCollisionType[irection] = 0;
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            //Determine whether the third collision body exists.
            if (cpz[irection][ispecies] > 0.0) 
            {
                isCollisionType[irection] = 1;
                break;
            }
        }
    }
}

void StoichiometricEquation::Create(int **forwardCoef, int **backwardCoef, RDouble **collisionCoef, int **reactionFlag)
{
    //! To create the coefficient of forward reaction equation.
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            cvp[r][s] = forwardCoef[r][s];
        }
    }

    //! To create the coefficient of backward reaction equation.
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            cvn[r][s] = backwardCoef[r][s];
        }
    }

    //! To create the coefficient of the third collision body in the reaction equation.
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            cpz[r][s] = collisionCoef[r][s];
        }
    }

    //! To generate the types of chemical reactions.
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        reactionType[r] = reactionFlag[r][1];
        if (reactionFlag[r][0] == 1)
            isCollisionType[r] = true;
        else
            isCollisionType[r] = false;
    }
}

void StoichiometricEquation::Read(DataContainer *cdata)
{
    for (int irection = 0; irection < numberOfReaction; ++ irection)
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            cdata->Read(&cvp[irection][ispecies], sizeof(int));
        }
    }

    for (int irection = 0; irection < numberOfReaction; ++ irection)
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            cdata->Read(&cvn[irection][ispecies], sizeof(int));
        }
    }

    for (int irection = 0; irection < numberOfReaction; ++ irection)
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            cdata->Read(&cpz[irection][ispecies], sizeof(RDouble));
        }
    }

    cdata->Read(reactionType, numberOfReaction * sizeof(int));
    cdata->Read(isCollisionType, numberOfReaction * sizeof(bool));
}

void StoichiometricEquation::Write(DataContainer *cdata)
{
    for (int irection = 0; irection < numberOfReaction; ++ irection)
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            cdata->Append(&cvp[irection][ispecies], sizeof(int));
        }
    }

    for (int irection = 0; irection < numberOfReaction; ++ irection)
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            cdata->Append(&cvn[irection][ispecies], sizeof(int));
        }
    }

    for (int irection = 0; irection < numberOfReaction; ++ irection)
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            cdata->Append(&cpz[irection][ispecies], sizeof(RDouble));
        }
    }

    cdata->Append(reactionType, numberOfReaction * sizeof(int));
    cdata->Append(isCollisionType, numberOfReaction * sizeof(bool));
}

CurveFitsMethod::CurveFitsMethod()
{
    blotter_ai = 0;
    blotter_bi = 0;
    blotter_ci = 0;
    blotter_di = 0;
    blotter_ei = 0;

    omega_ai = 0;
    omega_bi = 0;
    omega_ci = 0;
    omega_di = 0;

    curvefit_av = 0;
    curvefit_bv = 0; 
    curvefit_cv = 0;

    curvefit_aes = 0;
    curvefit_bes = 0; 
    curvefit_ces = 0;

    numberOfSpecies = 0;
}

CurveFitsMethod::~CurveFitsMethod()
{
    if (numberOfSpecies > 0)
    {
        delete [] blotter_ai;
        delete [] blotter_bi;
        delete [] blotter_ci;
        delete [] blotter_di;
        delete [] blotter_ei;

        delete [] omega_ai;
        delete [] omega_bi;
        delete [] omega_ci;
        delete [] omega_di;

        delete [] curvefit_av;
        delete [] curvefit_bv; 
        delete [] curvefit_cv;

        delete [] curvefit_aes;
        delete [] curvefit_bes;
        delete [] curvefit_ces;

        numberOfSpecies = 0;
    }
}

void CurveFitsMethod::Init(int numberOfSpecies, int nGasModel/* = 0*/)
{
    if (numberOfSpecies <= 0)
    {
        return;
    }
    this->numberOfSpecies = numberOfSpecies;

    blotter_ai = new RDouble[numberOfSpecies];
    blotter_bi = new RDouble[numberOfSpecies];
    blotter_ci = new RDouble[numberOfSpecies];
    blotter_di = new RDouble[numberOfSpecies];
    blotter_ei = new RDouble[numberOfSpecies];

    omega_ai = new RDouble[numberOfSpecies];
    omega_bi = new RDouble[numberOfSpecies];
    omega_ci = new RDouble[numberOfSpecies];
    omega_di = new RDouble[numberOfSpecies];

    curvefit_aes = new RDouble[numberOfSpecies];
    curvefit_bes = new RDouble[numberOfSpecies];
    curvefit_ces = new RDouble[numberOfSpecies];

    //! Set curve fit coefficients of vibrational excitation rate of electron collision for Nitrogen.
    SetNitrogenVibrationalExcitationRate();

    //! Set the curve fit coefficients of effective collision cross-sectional area between the electron and neutral heavy particle.
    SetEffectiveCollisionCrossSectionalAreaCurveFits(nGasModel);
}

//! Set the curve fit coefficients of vibrational excitation rate of electron collision-N2.
void CurveFitsMethod::SetNitrogenVibrationalExcitationRate()
{
    RDouble av[10] = {8.034, 7.924, 7.876, 7.626, 7.326, 4.9, 2.457, 1.119, 0.4681, 0.1837};
    RDouble bv[10] = {-2.227, -2.235, -2.257, -2.334, -2.454, -2.556, -2.702, -2.865, -3.042, -3.223};
    RDouble cv[10] = {2.005, 1.479, 1.054, 0.6499, 0.2049, 7.448e-3, 2.952e-3, 1.133e-3, 4.312e-3, 2.219e-4};

    int nLevel = 10;
    curvefit_av = new RDouble[nLevel];
    curvefit_bv = new RDouble[nLevel];
    curvefit_cv = new RDouble[nLevel];
    for (int i = 0; i < nLevel; i ++)
    {
        curvefit_av[i] = av[i];
        curvefit_bv[i] = bv[i];
        curvefit_cv[i] = cv[i];
    }
}

//! Set the curve fit coefficients of effective collision cross-sectional area between the electron and neutral heavy particle.
void CurveFitsMethod::SetEffectiveCollisionCrossSectionalAreaCurveFits(int nGasModel /*= 0*/)
{
    //Curve fitting data for electronic collision area.
    RDouble aes[8] = {1.2e-20, 2.0e-20, 1.0e-19, 5.0e-20, 7.5e-20, 1.0e-20, 1.0e-20, 1.0e-20};
    RDouble bes[8] = {1.7e-24, 6.0e-24, 0.0, 0.0, 5.5e-24, 0.0, 0.0, 0.0};
    RDouble ces[8] = {-2.0e-29, 0.0, 0.0, 0.0, -1.0e-28, 0.0, 0.0, 0.0};

    if (nGasModel == 0) //Earth Gas.
    {
        int nStart = MIN(4, numberOfSpecies);
        for (int i = 0; i < nStart; ++ i)
        {
            curvefit_aes[i] = aes[i];
            curvefit_bes[i] = bes[i];
            curvefit_ces[i] = ces[i];
        }
        for (int i = nStart; i < numberOfSpecies; ++ i)
        {
            curvefit_aes[i] = 0.0;
            curvefit_bes[i] = 0.0;
            curvefit_ces[i] = 0.0;
        }
        if (numberOfSpecies <= 5)
        {
            curvefit_aes[numberOfSpecies - 1] = aes[4];
            curvefit_bes[numberOfSpecies - 1] = bes[4];
            curvefit_ces[numberOfSpecies - 1] = ces[4];
        }
        else
        {
            curvefit_aes[numberOfSpecies - 2] = aes[4];
            curvefit_bes[numberOfSpecies - 2] = bes[4];
            curvefit_ces[numberOfSpecies - 2] = ces[4];
        }
    }
    else
    {
        for (int i = 0; i < numberOfSpecies; ++ i)
        {
            curvefit_aes[i] = aes[i];
            curvefit_bes[i] = bes[i];
            curvefit_ces[i] = ces[i];
        }
    }
}

void CurveFitsMethod::SetEffectiveCollisionCrossSectionalAreaCurveFits(vector<string> &namelist)
{
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        if (namelist[s] == "O")
        {
            curvefit_aes[s] = 1.2e-20;
            curvefit_bes[s] = 1.7e-24;
            curvefit_ces[s] = -2.0e-29;
        }
        else if (namelist[s] == "O2")
        {
            curvefit_aes[s] = 2.0e-20;
            curvefit_bes[s] = 6.0e-24;
            curvefit_ces[s] = 0.0;
        }
        else if (namelist[s] == "NO")
        {
            curvefit_aes[s] = 1.0e-19;
            curvefit_bes[s] = 0.0;
            curvefit_ces[s] = 0.0;
        }
        else if (namelist[s] == "N")
        {
            curvefit_aes[s] = 5.0e-20;
            curvefit_bes[s] = 0.0;
            curvefit_ces[s] = 0.0;
        }
        else if (namelist[s] == "N2")
        {
            curvefit_aes[s] = 7.5e-20;
            curvefit_bes[s] = 5.5e-24;
            curvefit_ces[s] = -1.0e-28;
        }
        else if (namelist[s] == "C")
        {
            curvefit_aes[s] = 1.0e-20;
            curvefit_bes[s] = 0.0;
            curvefit_ces[s] = 0.0;
        }
        else if (namelist[s] == "CO")
        {
            curvefit_aes[s] = 1.0e-20;
            curvefit_bes[s] = 0.0;
            curvefit_ces[s] = 0.0;
        }
        else if (namelist[s] == "CO2")
        {
            curvefit_aes[s] = 1.0e-20;
            curvefit_bes[s] = 0.0;
            curvefit_ces[s] = 0.0;
        }
        else
        {
            curvefit_aes[s] = 0.0;
            curvefit_bes[s] = 0.0;
            curvefit_ces[s] = 0.0;
        }
    }
}

void CurveFitsMethod::Read(fstream &file)
{
    string line, word;
    string separator = " =\t\r\n#$,;\"'";

    //! Obtain the curve fit coefficients of Gupta and Blottner methods for computing the viscosities of species.
    SkipLines(file, 2);
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        SkipLines(file, 3);
        getline(file, line);

        line = FindNextWord(line, word, separator);
        from_string< RDouble >(blotter_ai[ispecies], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string< RDouble >(blotter_bi[ispecies], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string< RDouble >(blotter_ci[ispecies], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string< RDouble >(blotter_di[ispecies], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string< RDouble >(blotter_ei[ispecies], word, std::dec);
    }

    //! Obtain the curve fit coefficients of computing the average collision areas for species couples called PI*Omega(2,2).
    SkipLines(file, 2);
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        SkipLines(file, 3);
        getline(file, line);

        line = FindNextWord(line, word, separator);
        from_string< RDouble >(omega_ai[ispecies], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string< RDouble >(omega_bi[ispecies], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string< RDouble >(omega_ci[ispecies], word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string< RDouble >(omega_di[ispecies], word, std::dec);
    }
}

void CurveFitsMethod::GenerateData(vector<string> &namelist)
{
    TransportFitData myCurveFitData;
    TransportFitData *currentData = new TransportFitData();

    myCurveFitData.GenerateFullCurveFittingData();
    myCurveFitData.CreateTransportFitData(namelist, currentData);

    RDouble **viscosityFitData = currentData->GetViscosityCurveFittingData();
    RDouble **collisionAreaFitData = currentData->GetCollisionAreaCurveFittingData();
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        blotter_ai[s] = viscosityFitData[s][0];
        blotter_bi[s] = viscosityFitData[s][1];
        blotter_ci[s] = viscosityFitData[s][2];
        blotter_di[s] = viscosityFitData[s][3];
        blotter_ei[s] = viscosityFitData[s][4];

        omega_ai[s] = collisionAreaFitData[s][0];
        omega_bi[s] = collisionAreaFitData[s][1];
        omega_ci[s] = collisionAreaFitData[s][2];
        omega_di[s] = collisionAreaFitData[s][3];
    }

    SetEffectiveCollisionCrossSectionalAreaCurveFits(namelist);
    delete currentData;
}

void CurveFitsMethod::Read(DataContainer *cdata)
{
    cdata->Read(blotter_ai, numberOfSpecies * sizeof(RDouble));
    cdata->Read(blotter_bi, numberOfSpecies * sizeof(RDouble));
    cdata->Read(blotter_ci, numberOfSpecies * sizeof(RDouble));
    cdata->Read(blotter_di, numberOfSpecies * sizeof(RDouble));
    cdata->Read(blotter_ei, numberOfSpecies * sizeof(RDouble));

    cdata->Read(omega_ai, numberOfSpecies * sizeof(RDouble));
    cdata->Read(omega_bi, numberOfSpecies * sizeof(RDouble));
    cdata->Read(omega_ci, numberOfSpecies * sizeof(RDouble));
    cdata->Read(omega_di, numberOfSpecies * sizeof(RDouble));

    cdata->Read(curvefit_av, 10 * sizeof(RDouble));
    cdata->Read(curvefit_bv, 10 * sizeof(RDouble));
    cdata->Read(curvefit_cv, 10 * sizeof(RDouble));

    cdata->Read(curvefit_aes, numberOfSpecies * sizeof(RDouble));
    cdata->Read(curvefit_bes, numberOfSpecies * sizeof(RDouble));
    cdata->Read(curvefit_ces, numberOfSpecies * sizeof(RDouble));
}

void CurveFitsMethod::Write(DataContainer *cdata)
{
    cdata->Append(blotter_ai, numberOfSpecies * sizeof(RDouble));
    cdata->Append(blotter_bi, numberOfSpecies * sizeof(RDouble));
    cdata->Append(blotter_ci, numberOfSpecies * sizeof(RDouble));
    cdata->Append(blotter_di, numberOfSpecies * sizeof(RDouble));
    cdata->Append(blotter_ei, numberOfSpecies * sizeof(RDouble));

    cdata->Append(omega_ai, numberOfSpecies * sizeof(RDouble));
    cdata->Append(omega_bi, numberOfSpecies * sizeof(RDouble));
    cdata->Append(omega_ci, numberOfSpecies * sizeof(RDouble));
    cdata->Append(omega_di, numberOfSpecies * sizeof(RDouble));

    cdata->Append(curvefit_av, 10 * sizeof(RDouble));
    cdata->Append(curvefit_bv, 10 * sizeof(RDouble));
    cdata->Append(curvefit_cv, 10 * sizeof(RDouble));

    cdata->Append(curvefit_aes, numberOfSpecies * sizeof(RDouble));
    cdata->Append(curvefit_bes, numberOfSpecies * sizeof(RDouble));
    cdata->Append(curvefit_ces, numberOfSpecies * sizeof(RDouble));
}

SchmidtNumber::SchmidtNumber()
{
    scl  = 0;
    sct  = 0;
    oscl = 0;
    osct = 0;
    numberOfSpecies = 0;
}

SchmidtNumber::~SchmidtNumber()
{
    delete [] scl ;
    delete [] sct ;
    delete [] oscl;
    delete [] osct;
}

void SchmidtNumber::Init(int numberOfSpecies)
{
    this->numberOfSpecies = numberOfSpecies;

    scl  = new RDouble[numberOfSpecies];
    sct  = new RDouble[numberOfSpecies];
    oscl = new RDouble[numberOfSpecies];
    osct = new RDouble[numberOfSpecies];
}

void SchmidtNumber::ComputeSpeciesSchmidtNumber(int *ionTypeOfSpecies)
{
    //! Calculate the Schmidt number of laminar or turbulent flow of species.

    RDouble sc_l;
    GlobalDataBase::GetData("sc_l", &sc_l, PHDOUBLE, 1);
    RDouble sc_t;
    GlobalDataBase::GetData("sc_t", &sc_t, PHDOUBLE, 1);

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        if (ionTypeOfSpecies[ispecies] != 0 && ABS(ionTypeOfSpecies[ispecies]) < 100)
        {
            scl[ispecies] = half * sc_l;
            sct[ispecies] = half * sc_t;
        }
        else
        {
            scl[ispecies] = 1.0 * sc_l;
            sct[ispecies] = 1.0 * sc_t;
        }

        oscl[ispecies] = 1.0 / scl[ispecies];
        osct[ispecies] = 1.0 / sct[ispecies];
    }
}

void SchmidtNumber::Read(DataContainer *cdata)
{
    cdata->Read(scl, numberOfSpecies * sizeof(RDouble));
    cdata->Read(sct, numberOfSpecies * sizeof(RDouble));
                                        
    cdata->Read(oscl, numberOfSpecies * sizeof(RDouble));
    cdata->Read(osct, numberOfSpecies * sizeof(RDouble));
}

void SchmidtNumber::Write(DataContainer *cdata)
{
    cdata->Append(scl, numberOfSpecies * sizeof(RDouble));
    cdata->Append(sct, numberOfSpecies * sizeof(RDouble));
                                          
    cdata->Append(oscl, numberOfSpecies * sizeof(RDouble));
    cdata->Append(osct, numberOfSpecies * sizeof(RDouble));
}

MolecularProperty::MolecularProperty()
{
    schmidt_number              = 0;
    molecularWeightDimensional  = 0;
    molecularWeight             = 0;
    oMolecularWeightDimensional = 0;
    oMolecularWeight            = 0;

    characteristicTemperatureOfSpecies = 0;
    mass_fraction_reference = 0;
    name_of_species         = 0;
    ionTypeOfSpecies        = 0;
    collisionCrossSection   = 0;

    thermodynamicProperties = 0;
    numberOfSpecies = 0;
}

MolecularProperty::~MolecularProperty()
{
    delete schmidt_number;
    delete [] molecularWeightDimensional ;
    delete [] molecularWeight            ;
    delete [] oMolecularWeightDimensional;
    delete [] oMolecularWeight           ;

    delete [] characteristicTemperatureOfSpecies;
    delete [] mass_fraction_reference;
    delete [] name_of_species        ;
    delete [] ionTypeOfSpecies     ;
    delete [] collisionCrossSection;
    delete [] thermodynamicProperties;
}

void MolecularProperty::Init(int numberOfSpecies)
{
    this->numberOfSpecies = numberOfSpecies;

    molecularWeightDimensional  = new RDouble[numberOfSpecies];
    molecularWeight             = new RDouble[numberOfSpecies];
    oMolecularWeightDimensional = new RDouble[numberOfSpecies];
    oMolecularWeight            = new RDouble[numberOfSpecies];

    characteristicTemperatureOfSpecies = new RDouble[numberOfSpecies];
    name_of_species = new string[numberOfSpecies];

    //! ionTypeOfSpecies is the ion type of species.
    ionTypeOfSpecies        = new int[numberOfSpecies];
    mass_fraction_reference = new RDouble[numberOfSpecies];
    collisionCrossSection   = new RDouble[numberOfSpecies];

    schmidt_number = new SchmidtNumber();
    schmidt_number->Init(numberOfSpecies);

    //! To set the thermodynamic properties of species under the multi-temperature model.
    thermodynamicProperties = new Thermo_Param[numberOfSpecies];
}

void MolecularProperty::Read(fstream &file)
{
    string line, word;
    string separator = " =\t\r\n#$,;\"'";

    SkipLines(file, 3);
    //! Obtain the name of species.
    getline(file, line);
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        line = FindNextWord(line, name_of_species[ispecies], separator);
    }

    SkipLines(file, 3);
    //! Obtain the ion type of species.
    getline(file, line);
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        line = FindNextWord(line, word, separator);
        from_string< int >(ionTypeOfSpecies[ispecies], word, std::dec);
    }

    SkipLines(file, 3);
    //! Obtain the molecular weight of species.
    getline(file, line);
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        line = FindNextWord(line, word, separator);
        from_string< RDouble >(molecularWeightDimensional[ispecies], word, std::dec);
    }

    SkipLines(file, 3);
    //! Obtain the mass ratio of species.
    getline(file, line);
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        line = FindNextWord(line, word, separator);
        from_string< RDouble >(mass_fraction_reference[ispecies], word, std::dec);
    }

    SkipLines(file, 3);
    //! Obtain the collision cross section of species.
    getline(file, line);
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        line = FindNextWord(line, word, separator);
        from_string< RDouble >(collisionCrossSection[ispecies], word, std::dec);
    }

    SkipLines(file, 3);
    //! Obtain the characteristic temperature of species.
    getline(file, line);
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        line = FindNextWord(line, word, separator);
        from_string< RDouble >(characteristicTemperatureOfSpecies[ispecies], word, std::dec);
    }

    schmidt_number->ComputeSpeciesSchmidtNumber(ionTypeOfSpecies);

    //! The thermodynamic properties under the multi-temperature model.
    int nParam = 6;     //! To add 6 parameters including the degeneracy, characteristic temperature etc.
    for (int iparam = 0; iparam < nParam; iparam ++)
    {
        SkipLines(file, 3);
        //! Obtain the characteristic temperature of species.
        getline(file, line);
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            line = FindNextWord(line, word, separator);
            if (iparam == 0)             //! To set the type of species.
            {
                from_string< int >(thermodynamicProperties[ispecies].nType, word, std::dec);
            }else if (iparam == 1)       //! To set the first degeneracy of species.
            {
                from_string< int >(thermodynamicProperties[ispecies].gs0, word, std::dec);
            }else if (iparam == 2)       //! To set the second degeneracy of species.
            {
                from_string< int >(thermodynamicProperties[ispecies].gs1, word, std::dec);
            }else if (iparam == 3)       //! To set the characteristic temperature of species in the vibration mode.
            {
                from_string< RDouble >(thermodynamicProperties[ispecies].Tvs.Tv[0], word, std::dec);
            }else if (iparam == 4)       //! To set the characteristic temperature of species in the electron mode.
            {
                from_string< RDouble >(thermodynamicProperties[ispecies].Tes, word, std::dec);
            }else if (iparam == 5)       //! To set the formation enthalpies of species.
            {
                from_string< RDouble >(thermodynamicProperties[ispecies].hs0, word, std::dec);
            }
        }
    }         //! End of reading the parameters to species.
}
void MolecularProperty::GenerateData(vector<string> &namelist, vector<RDouble> &initMassFractions)
{
    AtmosphereData *baseData = new AtmosphereData();
    AtmosphereData *currentData = new AtmosphereData();

    baseData->GenerateFullPhysicalChemicalData();
    baseData->CreateAtmosphereData(namelist, currentData);

    //Obtain the new data according to the name list.
    int *speciesType = currentData->GetSpeciesType();
    int *speciesCharge = currentData->GetSpeciesCharge();
    int *speciesG0 = currentData->GetSpeciesDegeneracy0();
    int *speciesG1 = currentData->GetSpeciesDegeneracy1();
    RDouble *speciesMass = currentData->GetSpeciesMolecularWeight();
    RDouble *speciesTtr = currentData->GetSpeciesTemperature();
    RDouble *speciesTes = currentData->GetSpeciesElectronTemperature();
    RDouble *speciesSigma = currentData->GetSpeciesCollisionArea();
    RDouble *speciesEnthalpy = currentData->GetSpeciesFormationEnthalpy();
    VibrationModeData *speciesVibrationModeData = currentData->GetSpeciesVibrationalModeData();

    for (int s = 0; s < this->numberOfSpecies; ++ s)
    {
        name_of_species[s] = namelist[s];
        mass_fraction_reference[s] = initMassFractions[s];
        //! Generate the ion type of species.
        ionTypeOfSpecies[s] = speciesCharge[s];
        //! Generate the molecular weight of species.
        molecularWeightDimensional[s] = speciesMass[s];
        //! Generate the collision cross section of species.
        collisionCrossSection[s] = speciesSigma[s];
        //! Generate the characteristic temperature of species.
        characteristicTemperatureOfSpecies[s] = speciesTtr[s];

        thermodynamicProperties[s].nType = speciesType[s];
        thermodynamicProperties[s].gs0 = speciesG0[s];
        thermodynamicProperties[s].gs1 = speciesG1[s];
        thermodynamicProperties[s].Tvs = speciesVibrationModeData[s];
        thermodynamicProperties[s].Tes = speciesTes[s];
        thermodynamicProperties[s].hs0 = speciesEnthalpy[s];
    }
    schmidt_number->ComputeSpeciesSchmidtNumber(ionTypeOfSpecies);

    delete baseData;
    delete currentData;
}

void MolecularProperty::Read(DataContainer *cdata)
{
    schmidt_number->Read(cdata);
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        cdata->ReadString(name_of_species[ispecies]);
    }

    cdata->Read(ionTypeOfSpecies, numberOfSpecies * sizeof(int));
    cdata->Read(mass_fraction_reference, numberOfSpecies * sizeof(RDouble));

    cdata->Read(collisionCrossSection, numberOfSpecies * sizeof(RDouble));
    cdata->Read(characteristicTemperatureOfSpecies, numberOfSpecies * sizeof(RDouble));

    cdata->Read(molecularWeightDimensional, numberOfSpecies * sizeof(RDouble));
    cdata->Read(molecularWeight, numberOfSpecies * sizeof(RDouble));

    cdata->Read(oMolecularWeightDimensional, numberOfSpecies * sizeof(RDouble));
    cdata->Read(oMolecularWeight, numberOfSpecies * sizeof(RDouble));

    cdata->Read(thermodynamicProperties, numberOfSpecies * sizeof(Thermo_Param));
}

void MolecularProperty::Write(DataContainer *cdata)
{
    schmidt_number->Write(cdata);
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        cdata->AppendString(name_of_species[ispecies]);
    }

    cdata->Append(ionTypeOfSpecies, numberOfSpecies * sizeof(int));
    cdata->Append(mass_fraction_reference, numberOfSpecies * sizeof(RDouble));

    cdata->Append(collisionCrossSection, numberOfSpecies * sizeof(RDouble));
    cdata->Append(characteristicTemperatureOfSpecies, numberOfSpecies * sizeof(RDouble));

    cdata->Append(molecularWeightDimensional, numberOfSpecies * sizeof(RDouble));
    cdata->Append(molecularWeight, numberOfSpecies * sizeof(RDouble));

    cdata->Append(oMolecularWeightDimensional, numberOfSpecies * sizeof(RDouble));
    cdata->Append(oMolecularWeight, numberOfSpecies * sizeof(RDouble));

    cdata->Append(thermodynamicProperties, numberOfSpecies * sizeof(Thermo_Param));
}

RDouble * MolecularProperty::GetLaminarSchmidtNumberReciprocal() const
{ 
    return schmidt_number->GetLaminarSchmidtNumberReciprocal();
}

RDouble * MolecularProperty::GetTurbulentSchmidtNumberReciprocal() const
{
    return schmidt_number->GetTurbulentSchmidtNumberReciprocal();
}

MixedGas::MixedGas()
{
    molecularProperty      = 0;
    thermodynamicManager   = 0;
    reactionRate           = 0;
    stoichiometricEquation = 0;
    CurveFits = 0;
    nGasModel = 0; //nGasModel = 0 indicates the gas of Earth, nGasModel = 1 is for gas of Mars.
    parkVDPower = 0.6;
    nTEnergyModel = 0; //The conventional method to computing temperature energy terms.
    nIsChemicalFreeze = 0;
    nIsSuperCatalytic = 1;
    polynomialCoefT2 = 0;
}

MixedGas::~MixedGas()
{
    if (nchem > 0)
    {
        delete molecularProperty;
        delete thermodynamicManager;
        delete reactionRate;
        delete stoichiometricEquation;
        delete CurveFits;

        DeAllocateWorkingSpace();
        RDouble *prim_inf = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("prim_inf"));
        delete [] prim_inf;
        RDouble *catalyticMassFraction = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("catalyticMassFraction"));
        delete [] catalyticMassFraction;
    }
}

RDouble MixedGas::ComputeReferenceTotalEnthalpy()
{
    RDouble totalEnthalpy = 0.0, staticEnthalpy = 0.0;
    RDouble *initMassFraction = this->GetInitMassFraction();
    RDouble initVibrationTemperature = 1.0;
    if (GlobalDataBase::IsExist("freestream_vibration_temperature", PHDOUBLE, 1))
    {
        initVibrationTemperature = GlobalDataBase::GetDoubleParaFromDB("freestream_vibration_temperature");
        initVibrationTemperature = initVibrationTemperature / referenceParameterDimensional.GetTemperature();
    }
    RDouble velocity = referenceParameterDimensional.GetVelocity();
    if (this->ntmodel == 1)
    {
        staticEnthalpy = this->GetMixtureGasEnthalpy(initMassFraction, 1.0);
    }
    else
    {
        
    }
    totalEnthalpy = (staticEnthalpy + 0.5) * velocity * velocity;

    return totalEnthalpy;
}

void MixedGas::AllocateWorkingSpace()
{
    moleFractions = new RDouble[numberOfSpecies];
    massFractions = new RDouble[numberOfSpecies];

    speciesViscosity = new RDouble[numberOfSpecies];
    speciesWeight = new RDouble[numberOfSpecies];

    speciesCp = new RDouble[numberOfSpecies];
    speciesCv = new RDouble[numberOfSpecies];
    speciesEnthalpy = new RDouble[numberOfSpecies];

    workSpecies = new RDouble[numberOfSpecies];
    speciesEnergyTr = new RDouble[numberOfSpecies];
    speciesEnergyVb = new RDouble[numberOfSpecies];
    speciesEnergyEe = new RDouble[numberOfSpecies];
    speciesConductivity = new RDouble[numberOfSpecies];

    dfsdx = new RDouble[nl + nchem];
    dfsdy = new RDouble[nl + nchem];
    dfsdz = new RDouble[nl + nchem];
    fvis  = new RDouble[nl + nchem + ntmodel - 1];

    maximumSpecies = new RDouble[numberOfSpecies];

    Thermo_Energy_temparay.Cps = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Cvs = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Cvtrs = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Cvvs = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Cves = new RDouble[numberOfSpecies];

    Thermo_Energy_temparay.Hs = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Es = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Etrs = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Evs = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Ees = new RDouble[numberOfSpecies];
    Thermo_Energy_temparay.Rms = new RDouble[numberOfSpecies];

    Thermo_Energy_temparay.forwardRates = new RDouble [numberOfReaction];
    Thermo_Energy_temparay.backwardRates = new RDouble [numberOfReaction];
    Thermo_Energy_temparay.collisionRates = new RDouble [numberOfReaction];
    Thermo_Energy_temparay.OmegaDerivative = NewPointer2 <RDouble> (numberOfSpecies, ntmodel - 1);

    Thermo_Energy_temparay.partialValues = NewPointer2 <RDouble> (numberOfReaction, nEquation);
    Thermo_Energy_temparay.partialValues2 = NewPointer2 <RDouble> (numberOfReaction, nEquation);

    //Obtain the interpolation coefficients for the polynomial.
    int n = this->thermodynamicManager->GetNumberofPolynomialCoefficient();
    polynomialCoefT2 = NewPointer2<RDouble>(numberOfSpecies, n);
}

void MixedGas::DeAllocateWorkingSpace()
{
    delete [] moleFractions;
    delete [] massFractions;
    delete [] speciesViscosity;
    delete [] speciesWeight;

    delete [] speciesCp;
    delete [] speciesCv;
    delete [] speciesEnthalpy;
    delete [] workSpecies;

    delete [] speciesEnergyTr;
    delete [] speciesEnergyVb;
    delete [] speciesEnergyEe;
    delete [] speciesConductivity;

    delete [] dfsdx;
    delete [] dfsdy;
    delete [] dfsdz;
    delete [] fvis;
    
    delete [] maximumSpecies;

    delete [] Thermo_Energy_temparay.Cps;
    delete [] Thermo_Energy_temparay.Cvs;
    delete [] Thermo_Energy_temparay.Cvtrs ;
    delete [] Thermo_Energy_temparay.Cvvs;
    delete [] Thermo_Energy_temparay.Cves;

    delete [] Thermo_Energy_temparay.Hs ;
    delete [] Thermo_Energy_temparay.Es ;
    delete [] Thermo_Energy_temparay.Etrs;
    delete [] Thermo_Energy_temparay.Evs;
    delete [] Thermo_Energy_temparay.Ees;
    delete [] Thermo_Energy_temparay.Rms;

    delete [] Thermo_Energy_temparay.forwardRates;
    delete [] Thermo_Energy_temparay.backwardRates;
    delete [] Thermo_Energy_temparay.collisionRates;
    DelPointer2(Thermo_Energy_temparay.OmegaDerivative);
    DelPointer2(Thermo_Energy_temparay.partialValues);
    DelPointer2(Thermo_Energy_temparay.partialValues2);
    DelPointer2(polynomialCoefT2);
}

void MixedGas::Primitive2Conservative(RDouble *prim, RDouble gama, RDouble Tv, RDouble Te, RDouble *q)
{
    using namespace IDX;
    RDouble em;
    RDouble &density = prim[IR];
    RDouble &um = prim[IU];
    RDouble &vm = prim[IV];
    RDouble &wm = prim[IW];
    
    //! Obtain the mass fraction of the last species whose label in the collection is ns-1.
    if (nchem == 1)
    {
        NormalizePrimitive(prim);
    }

    //! Obtain the total internal energy Em.
    ComputeInternalEnergy(prim, gama, Tv, Te, em);

    q[IR ] = density;
    q[IRU] = density * um;
    q[IRV] = density * vm;
    q[IRW] = density * wm;
    q[IRE] = density * em;

    for (int m = nm; m < nl; ++ m)
    {
        q[m] = density * prim[m];
    }
    if (nchem == 1)
    {
        q[nl] = density * prim[nl];
    }
    for (int i = 0; i < ntmodel - 1; i ++)    //! multi-temperature model.
    {
        q[nl + nchem + i] = density * prim[nl + nchem + i];
    }
}

void MixedGas::Primitive2ConservativeR(RDouble *prim, RDouble gama, RDouble Tv, RDouble Te, RDouble staticE, RDouble *q)
{
    using namespace IDX;
    RDouble em;
    RDouble &density = prim[IR];
    RDouble &um = prim[IU];
    RDouble &vm = prim[IV];
    RDouble &wm = prim[IW];

    RDouble squareVelocity;
    squareVelocity = um * um + vm * vm + wm * wm;
    em = staticE + half * squareVelocity;

    q[IR ] = density;
    q[IRU] = density * um;
    q[IRV] = density * vm;
    q[IRW] = density * wm;
    q[IRE] = density * em;

    for (int m = nm; m < nEquation; ++ m)
    {
        q[m] = density * prim[m];
    }
}

void MixedGas::Conservative2Primitive(RDouble *q, RDouble gama, RDouble *prim, RDouble *temperature)
{
    using namespace IDX;
    RDouble density, oDensity, um, vm, wm, pressure, totalInternalEnergy, squareVelocity, internalEnergy;
    density = q[IR];
    if (density <= zero)
    {
        prim[IR] = density;
        return;
    }

    oDensity = 1.0 / density;
    um = q[IU] * oDensity;
    vm = q[IV] * oDensity;
    wm = q[IW] * oDensity;
    totalInternalEnergy  = q[IP];    ///den*H
    squareVelocity = um * um + vm * vm + wm * wm;
    internalEnergy = totalInternalEnergy - half * density * squareVelocity;

    prim[IR] = density;
    prim[IU] = um;
    prim[IV] = vm;
    prim[IW] = wm;

    if (nchem == 1)
    {
        for (int m = nm; m < nl; ++m)
        {
            prim[m] = MIN(static_cast<RDouble>(one), MAX(SMALL, q[m] * oDensity));
        }
        NormalizePrimitive(prim);
    }

    for (int i = 0; i < ntmodel - 1; ++ i)    //! multi-temperature model.
    {
        prim[nl + nchem + i] = MAX(q[nl + nchem + i] * oDensity, 1.0e-30);
    }

    int nIdealState = GlobalDataBase::GetIntParaFromDB("nIdealState");
    if (nIdealState == 1)
    {
        prim[IP] = (gama - one) * internalEnergy;
        //! Compute temperature.
        //RDouble gasConstant = this->GetCoefficientOfStateEquation();
        //temperature[ITT] = prim[IP] / (prim[IR] * gasConstant);
        //temperature[ITV] = temperature[ITT];
        //temperature[ITE] = temperature[ITT];
        GetTemperature(prim, temperature[ITT], temperature[ITV], temperature[ITE]);
    }
    else
    {
        //GetPressure(prim, gama, internalEnergy, temperature[0], pressure);
        this->GetTemperatureAndPressure(prim, internalEnergy, temperature, pressure);
        prim[IP] = MAX(SMALL, pressure);
    }
}

void MixedGas::Conservative2PrimitiveR(RDouble *q, RDouble gama, RDouble *prim, RDouble *temperature)
{
    using namespace IDX;
    RDouble gasConstant = this->GetUniversalGasConstant();
    RDouble *speciesMass = molecularProperty->GetMolecularWeight();

    RDouble oDensity, squareVelocity, internalEnergy;
    RDouble massReciprocal, ceDivideMe;
    RDouble gasMixtureCv, gasMixtureEnergy, TkLimit;
    int nCount;

    RDouble &density = prim[IR];
    RDouble &um = prim[IU];
    RDouble &vm = prim[IV];
    RDouble &wm = prim[IW];
    RDouble &pressure = prim[IP];

    density = q[IR];    
    oDensity = 1.0 / density;
    um = q[IRU] * oDensity;
    vm = q[IRV] * oDensity;
    wm = q[IRW] * oDensity;
 
    squareVelocity = um * um + vm * vm + wm * wm;
    internalEnergy = q[IRE] * oDensity - half * squareVelocity;

    for (int m = nm; m < nl; ++ m)
    {
        prim[m] = MIN(static_cast<RDouble>(one), MAX(SMALL, q[m] * oDensity));
    }

    NormalizePrimitive(prim);

    for (int i = 1; i < ntmodel; ++ i)    //! multi-temperature model.
    {
        prim[nl + i] = MAX(q[nl + i] * oDensity, 1.0e-30);
    }

    //! Obtain the average molecular weight.
    massReciprocal = 0.0;
    for (int i = 0; i < numberOfSpecies; ++ i)
    {
        massFractions[i] = prim[nm + i];
        massReciprocal += massFractions[i] / speciesMass[i];
    }
    
    if (nElectronIndex >= 0)
    {
        ceDivideMe = massFractions[nElectronIndex] / speciesMass[nElectronIndex];
    }
    else
    {
        ceDivideMe = 0.0;
    }

    //compute T and P;

    int nMaxStepTemperature = 5;
    RDouble precisionTemperature = 1.0e-3;

    if (ntmodel == 1)
    {
        RDouble &trTemperature = temperature[mTT];
        nCount = 0;
        while (nCount < nMaxStepTemperature)
        {
            nCount++;
            ComputeMixtureCvAndEnergy(massFractions, trTemperature, gasMixtureCv, gasMixtureEnergy);
           
            TkLimit = (gasMixtureEnergy - internalEnergy) / max(gasMixtureCv, SMALL);

            if (TkLimit > 0.9 * trTemperature)
            {
                TkLimit = 0.9 * trTemperature;
            }
            else if (TkLimit < -9.0 * trTemperature)
            {
                TkLimit = -9.0 * trTemperature;
            }
            trTemperature = MAX(trTemperature - TkLimit, trTemperatureMinNonDimensional);

            if (fabs(TkLimit / trTemperature) <= precisionTemperature)
            {
                break;
            }
        }
        //trTemperature = MAX(trTemperature, SMALL);
        pressure = density * gasConstant * trTemperature * massReciprocal;
    }
    else if (ntmodel == 2)
    {

    }
    else
    {

    }
}

RDouble * MixedGas::GetLaminarSchmidtNumberReciprocal() const
{ 
    return molecularProperty->GetLaminarSchmidtNumberReciprocal();
}

RDouble * MixedGas::GetTurbulentSchmidtNumberReciprocal() const
{
    return molecularProperty->GetTurbulentSchmidtNumberReciprocal();
}

//! To obtain the laminar Schmidt number of each component.
void MixedGas::GetLaminarSchmidtNumber(RDouble *Scl) const
{
    if (Scl == NULL)
    {
        return;
    }
    else
    {
        RDouble *pSchmidt = molecularProperty->GetLaminarSchmidtNumberReciprocal();
        for (int i = 0; i < numberOfSpecies; ++ i)
        {
            Scl[i] = 1.0 / pSchmidt[i];
        }
    }
}

//! To obtain the turbulent Schmidt number of each component.
void MixedGas::GetTurbulentSchmidtNumber(RDouble *Sct) const
{
    if (Sct == NULL)
    {
        return;
    }
    else
    {
        RDouble *pSchmidt = molecularProperty->GetTurbulentSchmidtNumberReciprocal();
        for (int i = 0; i < numberOfSpecies; ++ i)
        {
            Sct[i] = 1.0 / pSchmidt[i];
        }
    }
}

void MixedGas::GetSpecificHeatRatio(RDouble *prim, RDouble &gama)
{
    if (nchem == 0)      //! Perfect gas.
    {
        gama = 1.4;
    }
    else                 //! The chemical non-equilibrium flow.
    {
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            massFractions[i] = prim[nm + i];
        }
        RDouble cp, cv;
        RDouble R = this->GetUniversalGasConstant();
        RDouble M = this->GetMixedGasMolecularWeight(massFractions);
        if (ntmodel >= 2)     //! Two-Temperature and Three-Temperature model.
        {

        }
        else                 //! The default setting is One-Temperature model.
        {
            RDouble Ttr = prim[IDX::IP] * M / (prim[IDX::IR] * R);
            this->ComputeConstantPressureSpecificHeatByMassFraction(massFractions, Ttr, cp);
            cv = cp - R/M;
        }

        gama = cp / cv;
    }
}

void MixedGas::GetSpecificHeatRatioAndTemperatute(RDouble *primitiveVariables, RDouble &gama, RDouble &transRotationTemperature, RDouble &vibrationTemperature, RDouble &electronTemperature)
  {
  
    GetTemperature(primitiveVariables, transRotationTemperature, vibrationTemperature, electronTemperature);

    if (nchem == 0)       //! Perfect gas.
    {
        gama = 1.4;
    }
    else                  //! The chemical non-equilibrium flow.
    {
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            massFractions[i] = primitiveVariables[nm + i];
        }
        RDouble cp, cv;
        RDouble gasConstant = this->GetUniversalGasConstant();
        RDouble mass = this->GetMixedGasMolecularWeight(massFractions);
        //! The NormalizePrimitive() function has been used in the function of GetTemperature().
        if (ntmodel >= 2)        //! Two-Temperature and Three-Temperature model.
        {

        }
        else                    //! The default setting is One-Temperature model.
        {
            this->ComputeConstantPressureSpecificHeatByMassFraction(massFractions, transRotationTemperature, cp);
            cv = cp - gasConstant / mass;
        }
        gama = cp / cv;
    }
}

void MixedGas::ReadGasModel()
{
    string gasfile;
    GlobalDataBase::GetData("gasfile", &gasfile, PHSTRING, 1);

    string chemicalModel = "";
    int nSpecies = 0, nReaction = 0;
    if (gasfile == "Gu5")       //! The Gupta 5-species-6-reactions model.
    {
        chemicalModel = "Gupta";
        nSpecies = 5;
        nReaction = 6;
    }
    else if (gasfile == "Gu7")  //! The Gupta 7-species-9-reactions model.
    {
        chemicalModel = "Gupta";
        nSpecies = 7;
        nReaction = 9;
    }
    else if (gasfile == "Gu11") //! The Gupta 11-species-20-reactions model.
    {
        chemicalModel = "Gupta";
        nSpecies = 11;
        nReaction = 20;
    }
    else if (gasfile == "Pa5")  //! The Park 5-species-17-reactions model.
    {
        chemicalModel = "Park";
        nSpecies = 5;
        nReaction = 17;
    }
    else if (gasfile == "Pa7")  //! The Park 7-species-22-reactions model.
    {
        chemicalModel = "Park";
        nSpecies = 7;
        nReaction = 22;
    }
    else if (gasfile == "Pa11") //! The Park 11-species-48-reactions model.
    {
        chemicalModel = "Park";
        nSpecies = 11;
        nReaction = 48;
    }
    else if (gasfile == "DK5")  //! The Dunn-Kang 5-species-11-reactions model.
    {
        chemicalModel = "Dunn-Kang";
        nSpecies = 5;
        nReaction = 11;
    }
    else if (gasfile == "DK7")  //! The Dunn-Kang 7-species-15-reactions model.
    {
        chemicalModel = "Dunn-Kang";
        nSpecies = 7;
        nReaction = 15;
    }
    else if (gasfile == "DK11") //! The Dunn-Kang 11-species-26-reactions model.
    {
        chemicalModel = "Dunn-Kang";
        nSpecies = 11;
        nReaction = 26;
    }
    else if (gasfile == "Mars-Pa8") //! The Park 8-species-12-reactions model for Mars.
    {
        chemicalModel = "Mars-Park";
        nSpecies = 8;
        nReaction = 12;
    }
    else if (gasfile == "Mars-Mc8") //! The McKenzie 8-species-12-reactions model for Mars.
    {
        chemicalModel = "Mars-McKenzie";
        nSpecies = 8;
        nReaction = 12;
    }
    else if (gasfile == "Combustion-12")
    {
        chemicalModel = "CombustionGas"; //The 12-species-20-reactions model of Combustion chamber Gas which reference AEROPH.
        nSpecies = 12;
        nReaction = 20;
    }
    else if (gasfile == "Gas-Mixture")
    {
        chemicalModel = "Gas-Mixture";
        nSpecies = 2;
        nReaction = 1;
    }
    else
    {
        if (gasfile.find("Gu") != gasfile.npos)
        {
            chemicalModel = "Gupta";
        }
        else if (gasfile.find("Pa") != gasfile.npos)
        {
            if (gasfile.find("Mars-Pa") != gasfile.npos)
            {
                chemicalModel = "Mars-Park";
            }
            else
            {
                chemicalModel = "Park";
            }
        }
        else if (gasfile.find("DK") != gasfile.npos)
        {
            chemicalModel = "Dunn-Kang";
        }
        else if (gasfile.find("Mars-Mc") != gasfile.npos)
        {
            chemicalModel = "Mars-McKenzie";
        }
    }

    if (chemicalModel == "")    //! Obtain the gas model from file.
    {
        fstream file;
        file.open(gasfile.c_str(), ios_base::in);

        ReadGasModel(file);

        file.close();
        file.clear();
    }
    else    //! Generate the gas model in the code.
    {
        GenerateGasModel(chemicalModel);
    }
}

void MixedGas::InitGasModel()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    if (nchem <= 0) return;

    if (myid == server)
    {
        ReadGasModel();
    }

    int number_of_processor = GetNumberOfProcessor();

    if (number_of_processor <= 1) return;

    int tag = 0;
    PH_BcastComposite(ChemCompressData, ChemDecompressData, server, tag);
}

void MixedGas::ComputeReferenceViscosityDimensionalChemical()
{
    RDouble reference_viscosity_dimensional;
    RDouble reference_tempereature = 1.0;

    RDouble *mass_fraction_reference = this->molecularProperty->GetInitialMassFraction();
    ComputeMoleFractionByMassFraction(mass_fraction_reference, moleFractions);

    if (ntmodel > 1)    //! Multi-Temperature Model.
    {

    }
    else         //! The default model: one-temperature model.
    { 
        ComputeOneTemperatureModelSpeciesViscosityDimensional(reference_tempereature, speciesViscosity);
        ComputeMixtureCoefficientByWilkeFormula(moleFractions, speciesViscosity, speciesWeight);
        ComputeMixtureByWilkeFormula(moleFractions, speciesViscosity, speciesWeight, reference_viscosity_dimensional);
    }

    referenceParameterDimensional.SetViscosity(reference_viscosity_dimensional);
}

void MixedGas::ComputeReferenceViscosityDimensionalNoChemical()
{
    RDouble t0i  = 288.15;
    RDouble tsi  = 110.4;
    RDouble mu0i = 1.7894E-05;

    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();
    RDouble reference_viscosity_dimensional = (t0i + tsi) / (refDimensionalTemperature + tsi) * pow(refDimensionalTemperature / t0i, 1.5) * mu0i;
    referenceParameterDimensional.SetViscosity(reference_viscosity_dimensional);
}

void MixedGas::ComputeReferenceTsuth()
{
    RDouble tsi = 110.4;

    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();

    RDouble tsuth = tsi / refDimensionalTemperature;
    GlobalDataBase::UpdateData("tsuth", &tsuth, PHDOUBLE, 1);
}

void MixedGas::ComputeReferenceViscosityDimensional()
{
    if (nchem > 0)
    {
        ComputeReferenceViscosityDimensionalChemical();
    }
    else
    {
        ComputeReferenceViscosityDimensionalNoChemical();
    }

    ComputeReferenceTsuth();
}

void MixedGas::ComputeReferenceSpecificHeatRatio()
{
    if (nchem > 0)
    {
        ComputeReferenceSpecificHeatRatioWithChemical();
    }
    else
    {
        RDouble refGama = 1.4;
        refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
        referenceParameterDimensional.SetAverageSpecificHeatRatio(refGama);
    }
}

void MixedGas::ComputeReferenceSpecificHeatRatioWithChemical()
{
    RDouble reference_temperature = 1.0;
    RDouble refGama = 1.4;
    RDouble *mass_fraction_reference = this->molecularProperty->GetInitialMassFraction();

    if (ntmodel > 1)       //! multi-temperature model.
    {

    }
    else                  //! one-temperature model.
    {
        RDouble cp_dimensional, cv_dimensional;
        ComputeSpeciesConstantPressureSpecificHeatDimensional(reference_temperature, speciesCp);
        ComputeMixtureByMassFraction(mass_fraction_reference, speciesCp, cp_dimensional);

        RDouble reference_average_molecular_weight_dimensional = referenceParameterDimensional.GetAverageMolecularWeight();
        RDouble o_reference_average_molecular_weight_dimensional = 1.0 / reference_average_molecular_weight_dimensional;

        cv_dimensional = cp_dimensional - rjmk * o_reference_average_molecular_weight_dimensional;
        refGama = cp_dimensional / cv_dimensional;
    }

    referenceParameterDimensional.SetAverageSpecificHeatRatio(refGama);
    GlobalDataBase::UpdateData("refGama", &refGama, PHDOUBLE, 1);
}

void MixedGas::ComputeReferenceGeneralGasConstant()
{
    RDouble reference_average_molecular_weight_dimensional = referenceParameterDimensional.GetAverageMolecularWeight();
    RDouble reference_average_general_gas_constant = rjmk / reference_average_molecular_weight_dimensional;

    referenceParameterDimensional.SetAverageGeneralGasConstant(reference_average_general_gas_constant);
    GlobalDataBase::UpdateData("reference_average_general_gas_constant", &reference_average_general_gas_constant, PHDOUBLE, 1);
}

void MixedGas::SetAirInformationByDataBase()
{
    RDouble height = GlobalDataBase::GetDoubleParaFromDB("height");

    RDouble refDimensionalTemperature;
    RDouble refDimensionalPressure;
    RDouble refDimensionalDensity;
    RDouble refDimensionalSonicSpeed;

    GetAirInfo(height, refDimensionalTemperature, refDimensionalPressure, refDimensionalDensity, refDimensionalSonicSpeed);

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    referenceParameterDimensional.SetPressure(refDimensionalPressure);
    referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
    referenceParameterDimensional.SetSoundVelocity(refDimensionalSonicSpeed);

    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalTemperature", &refDimensionalTemperature, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalSonicSpeed", &refDimensionalSonicSpeed, PHDOUBLE, 1);
}

void MixedGas::SetAirInformationByExpData()
{
    RDouble refDimensionalTemperature;
    RDouble refDimensionalPressure;
    RDouble refDimensionalDensity;
    RDouble refDimensionalSonicSpeed;

    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    RDouble T0 = referenceParameterDimensional.GetTemperature();
    RDouble P0 = referenceParameterDimensional.GetPressure();

    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble R = 287.053;

    RDouble mi = refGama / (refGama - 1.0);
    RDouble coef = 1.0 + (refGama - 1.0) * refMachNumber * refMachNumber / 2.0;

    RDouble T, P, rho, c;
    T = T0 / coef;
    P = P0 / pow(coef, mi);
    rho = P / R / T;
    c = sqrt(refGama * R * T);

    refDimensionalDensity = rho;
    refDimensionalPressure = P;
    refDimensionalTemperature = T;
    refDimensionalSonicSpeed = c;

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    referenceParameterDimensional.SetPressure(refDimensionalPressure);
    referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
    referenceParameterDimensional.SetSoundVelocity(refDimensionalSonicSpeed);

    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalTemperature", &refDimensionalTemperature, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalSonicSpeed", &refDimensionalSonicSpeed, PHDOUBLE, 1);
}

void MixedGas::ComputeReferenceGasInformation()
{
    int inflowParaType = GlobalDataBase::GetIntParaFromDB("inflowParaType");

    if (inflowParaType == NONDIMENSIONCONDITION)    //! Obtain the reference density and reference pressure for non-dimension condition.
    {
        ComputeDensityWithReynoldsNumber();
        ComputePressureInGasStateEquation();
    }
    else if (inflowParaType == FLIGHTCONDITION || inflowParaType == EXPERIMENTCONDITION)      //! Obtain the reference density for non-dimension, and the Re.
    {
        NormalizeAirInformation();
        ComputeReferenceReynoldsNumber();
    }
    else if (inflowParaType == TEMPERATURE_DENSITY)
    {
        RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();
        RDouble refDimensionalDensity = referenceParameterDimensional.GetDensity();
        RDouble gasConstant = referenceParameterDimensional.GetAverageGeneralGasConstant();
        
        RDouble refDimensionalPressure = refDimensionalDensity * gasConstant * refDimensionalTemperature;
        referenceParameterDimensional.SetPressure(refDimensionalPressure);
        GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);

        ComputeReferenceReynoldsNumber();
    }
    else if (inflowParaType == TEMPERATURE_PRESSURE)
    {
        NormalizeAirInformation();
        ComputeReferenceReynoldsNumber();
    }
    else
    {
        TK_Exit::UnexpectedVarValue("inflowParaType", inflowParaType);
    }
}

void MixedGas::ComputeReferenceSoundVelocity()
{
    RDouble gama0 = referenceParameterDimensional.GetAverageSpecificHeatRatio();
    RDouble reference_general_gas_constant = referenceParameterDimensional.GetAverageGeneralGasConstant();
    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();

    //! Unit: m/s.
    RDouble refDimensionalSonicSpeed = sqrt(gama0 * reference_general_gas_constant * refDimensionalTemperature);

    referenceParameterDimensional.SetSoundVelocity(refDimensionalSonicSpeed);
    GlobalDataBase::UpdateData("refDimensionalSonicSpeed", &refDimensionalSonicSpeed, PHDOUBLE, 1);
}

void MixedGas::ComputeReferenceInitMassFractionInformation()
{
    RDouble *mass_fraction_reference = this->molecularProperty->GetInitialMassFraction();

    MassFractionConversion(mass_fraction_reference);

}

void MixedGas::ComputeReferenceMolecularInformation()
{
    if (nchem > 0)
    {
        ComputeReferenceMolecularInformationWithChemical();
    }
    else
    {
        ComputeReferenceMolecularInformationWithoutChemical();
    }

    RDouble refMass = referenceParameterDimensional.GetAverageMolecularWeight();
    this->electronMassDimension = 5.486e-7;     //! Unit: kg/m3.
    this->electronMass = 5.486e-7 / refMass;    //! non-dimensional value.
}

void MixedGas::ComputeReferenceMolecularInformationWithoutChemical()
{
    //! Unit: kg/mol.
    RDouble reference_average_molecular_weight_dimensional = 28.9644 * 1.0e-3;

    referenceParameterDimensional.SetAverageMolecularWeight(reference_average_molecular_weight_dimensional);
}

void MixedGas::ComputeReferenceMolecularInformationWithChemical()
{
    RDouble *molecularWeight  = molecularProperty->GetMolecularWeight();
    RDouble *oMolecularWeight = molecularProperty->GetMolecularWeightReciprocal();

    RDouble *molecularWeightDimensional  = molecularProperty->GetMolecularWeightDimensional();
    RDouble *oMolecularWeightDimensional = molecularProperty->GetMolecularWeightReciprocalDimensional();

    RDouble *mass_fraction_reference = this->molecularProperty->GetInitialMassFraction();

    //! The molecular weight of the species.
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        molecularWeightDimensional[ispecies] *= 1.0e-3;
    }

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        oMolecularWeightDimensional[ispecies] = 1.0 / molecularWeightDimensional[ispecies];
    }

    RDouble o_reference_average_molecular_weight_dimensional = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        o_reference_average_molecular_weight_dimensional += mass_fraction_reference[ispecies] * oMolecularWeightDimensional[ispecies];
    }

    RDouble reference_average_molecular_weight_dimensional = 1.0 / o_reference_average_molecular_weight_dimensional;

    referenceParameterDimensional.SetAverageMolecularWeight(reference_average_molecular_weight_dimensional);

    //! The non-dimensionalization of the molecular weight of species.
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        molecularWeight[ispecies] = molecularWeightDimensional[ispecies] * o_reference_average_molecular_weight_dimensional;
    }

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        oMolecularWeight[ispecies] = 1.0 / molecularWeight[ispecies];
    }
}

void MixedGas::ComputeReferenceVelocity()
{
    RDouble refDimensionalSonicSpeed = referenceParameterDimensional.GetSoundVelocity();

    int inflowParaType = GlobalDataBase::GetIntParaFromDB("inflowParaType");

    if (inflowParaType == TEMPERATURE_DENSITY || inflowParaType == TEMPERATURE_PRESSURE)
    {
        RDouble refDimensionalVelocity = referenceParameterDimensional.GetVelocity();
        RDouble refMachNumber = refDimensionalVelocity / refDimensionalSonicSpeed;
        referenceParameterDimensional.SetMachNumber(refMachNumber);
        GlobalDataBase::UpdateData("refMachNumber", &refMachNumber, PHDOUBLE, 1);
    }
    else
    {
        RDouble refMachNumber = referenceParameterDimensional.GetMachNumber();
        RDouble refDimensionalVelocity = refDimensionalSonicSpeed * refMachNumber;
        referenceParameterDimensional.SetVelocity(refDimensionalVelocity);
        GlobalDataBase::UpdateData("refDimensionalVelocity", &refDimensionalVelocity, PHDOUBLE, 1);
    }
}

void MixedGas::ComputeReferenceReynoldsNumber()
{
    RDouble refDimensionalDensity = referenceParameterDimensional.GetDensity();
    RDouble reference_velocity_dimensional  = referenceParameterDimensional.GetVelocity();
    RDouble reference_viscosity_dimensional = referenceParameterDimensional.GetViscosity();

    RDouble refReNumber = refDimensionalDensity * reference_velocity_dimensional / reference_viscosity_dimensional;
    RDouble oreynolds = 1.0 / refReNumber;

    referenceParameterDimensional.SetReynoldsNumber(refReNumber);
    GlobalDataBase::UpdateData("refReNumber", &refReNumber, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("oreynolds", &oreynolds, PHDOUBLE, 1);
}

void MixedGas::ComputeProperReynoldsNumberForGrid()
{
    RDouble refReNumber = referenceParameterDimensional.GetReynoldsNumber();
    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

    refReNumber *= reynoldsReferenceLengthDimensional;
    RDouble oreynolds = 1.0 / refReNumber;

    referenceParameterDimensional.SetReynoldsNumber(refReNumber);
    GlobalDataBase::UpdateData("refReNumber", &refReNumber, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("oreynolds", &oreynolds, PHDOUBLE, 1);
}

void MixedGas::ComputeOtherProperParameterForGrid()
{
    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");
    RDouble forceReferenceLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
    RDouble forceReferenceLengthSpanWise = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");
    RDouble forceReferenceArea = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");

    RDouble TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    RDouble TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    RDouble TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

    if (GetDim() == THREE_D)
    {
        forceReferenceArea /= reynoldsReferenceLengthDimensional * reynoldsReferenceLengthDimensional;
    }
    else if (GetDim() == TWO_D)
    {
        forceReferenceArea /= reynoldsReferenceLengthDimensional;
    }
    forceReferenceLength /= reynoldsReferenceLengthDimensional;
    forceReferenceLengthSpanWise /= reynoldsReferenceLengthDimensional;
    TorqueRefX /= reynoldsReferenceLengthDimensional;
    TorqueRefY /= reynoldsReferenceLengthDimensional;
    TorqueRefZ /= reynoldsReferenceLengthDimensional;

    GlobalDataBase::UpdateData("forceReferenceLength", &forceReferenceLength, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("forceReferenceLengthSpanWise", &forceReferenceLengthSpanWise, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("forceReferenceArea", &forceReferenceArea, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("TorqueRefZ", &TorqueRefZ, PHDOUBLE, 1);

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady == 1)
    {
        PrintToWindow("Unsteady Informations" , "\n");
        RDouble reference_velocity_dimensional = referenceParameterDimensional.GetVelocity();
        RDouble unit = reynoldsReferenceLengthDimensional / reference_velocity_dimensional;
        RDouble physicalTimeStepDimensional = 0.0;
        RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");
        if (GlobalDataBase::IsExist("physicalTimeStepDimensional", PHDOUBLE, 1))
        {
            physicalTimeStepDimensional = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStepDimensional");
            if (physicalTimeStepDimensional >= 0.0)
            {
                physicalTimeStep = physicalTimeStepDimensional / unit;
            }
        }
        physicalTimeStepDimensional = physicalTimeStep * unit;
        GlobalDataBase::UpdateData("physicalTimeStepDimensional", &physicalTimeStepDimensional, PHDOUBLE, 1);
        //! The non-dimensional physical time step based on grid scale.
        physicalTimeStep /= reynoldsReferenceLengthDimensional;
        GlobalDataBase::UpdateData("physicalTimeStep", &physicalTimeStep, PHDOUBLE, 1);

        PrintToWindow("  Physical Time Step       :    ", physicalTimeStepDimensional , "s", "\n");
    }
}

void MixedGas::ComputeDensityWithReynoldsNumber()
{
    RDouble refReNumber = referenceParameterDimensional.GetReynoldsNumber();
    RDouble reference_viscosity_dimensional = referenceParameterDimensional.GetViscosity();
    RDouble reference_velocity_dimensional = referenceParameterDimensional.GetVelocity();

    RDouble refDimensionalDensity = refReNumber * reference_viscosity_dimensional / reference_velocity_dimensional;

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
}

void MixedGas::ComputePressureInGasStateEquation()
{
    RDouble refDimensionalDensity = referenceParameterDimensional.GetDensity();
    RDouble reference_average_generanl_gas_constant = referenceParameterDimensional.GetAverageGeneralGasConstant();
    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();

    RDouble refDimensionalPressure = refDimensionalDensity * reference_average_generanl_gas_constant * refDimensionalTemperature;

    referenceParameterDimensional.SetPressure(refDimensionalPressure);
    GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);
}

void MixedGas::InitNumberOfSpecies()
{
    if (nchem == 0)
    {
        numberOfSpecies = 0;
        nl = nm;
        nEquation = nm;
    }
    else
    {
        nl = nm + numberOfSpecies - 1;
        nEquation = nm + numberOfSpecies + ntmodel - 1;
    }
    GlobalDataBase::UpdateData("nl", &nl, PHINT, 1);
    GlobalDataBase::UpdateData("numberOfSpecies", &numberOfSpecies, PHINT, 1);

    if (numberOfSpecies > 0)
    {
        this->nElectronIndex = this->GetSpeciesIndex("e-");
        nNitrogenIndex = this->GetSpeciesIndex("N2");
        if (this->nElectronIndex >= 0)
        {
            this->nSpeciesEquation = numberOfSpecies - 2;
            this->nLastSpeciesIndex = numberOfSpecies - 2;
        }
        else
        {
            this->nSpeciesEquation = numberOfSpecies - 1;
            this->nLastSpeciesIndex = numberOfSpecies - 1;
        }
    }
    else
    {
        this->nSpeciesEquation = 0;
        this->nElectronIndex = -1;
        this->nNitrogenIndex = -1;
        this->nLastSpeciesIndex = -1;
    }

    //! This variable is used to indicate the index of the species with the maximum mass fraction in the array of mixture.
    //for Earth air, Nitrogen is the species with the maximum mass fraction, so it is employed as the variable name.
    GlobalDataBase::UpdateData("nNitrogenIndex", &this->nLastSpeciesIndex, PHINT, 1);
    GlobalDataBase::UpdateData("nElectronIndex", &this->nElectronIndex, PHINT, 1);

    if (nchem > 0)
    {
        AllocateWorkingSpace();
    }
}

RDouble MixedGas::ComputeReferenceReynoldsNumberFromTotalPT(void)
{
    RDouble Re0, refMachNumber;
    RDouble T0 = referenceParameterDimensional.GetTemperature();
    RDouble P0 = referenceParameterDimensional.GetPressure();
    GlobalDataBase::GetData("refMachNumber", &refMachNumber, PHDOUBLE, 1);

    RDouble refGama, c;
    GlobalDataBase::GetData("refGama", &refGama, PHDOUBLE, 1);
    c = referenceParameterDimensional.GetSoundVelocity();
    RDouble R = 287.0;
    RDouble coef, T, P, rho, miuInf;
    RDouble mi = refGama / (refGama - 1.0);
    coef = 1.0 + (refGama - 1.0) * refMachNumber * refMachNumber / 2.0;
    T = T0 / coef;
    P = P0 / pow(coef, mi);
    rho = P / R / T;

    RDouble miu0 = 1.7894e-5;
    RDouble T1 = 288.15, C = 110.4;
    miuInf = pow(T / T1, 1.5) * (T1 + C) / (T + C) * miu0;

    Re0 = rho * refMachNumber * c / miuInf;
    return Re0;
}

void MixedGas::Print()
{
    ostringstream oss;
    string *speciesName = GetNameOfSpecies();
    Thermo_Param *param = molecularProperty->GetThermodynamicProperties();
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        oss << speciesName[s] << endl;
        oss << param[s].nType << ", " << param[s].Tvs.Tv[0] << ", " << param[s].Tes << ", " << param[s].gs0 << ", " << param[s].gs1 << ", " << param[s].hs0 << endl;
    }

    RDouble *blotterA = CurveFits->GetBlotterCurveFittingCoefficientA();
    RDouble *blotterB = CurveFits->GetBlotterCurveFittingCoefficientB();
    RDouble *blotterC = CurveFits->GetBlotterCurveFittingCoefficientC();
    RDouble *blotterD = CurveFits->GetBlotterCurveFittingCoefficientD();
    RDouble *blotterE = CurveFits->GetBlotterCurveFittingCoefficientE();
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        oss << speciesName[s] << endl;
        oss << blotterA[s] << ", " << blotterB[s] << ", " << blotterC[s] << ", " << blotterD[s] << ", " << blotterE[s] << endl;
    }

    int nTemp = 5, nCoef = 0;
    nCoef = thermodynamicManager->GetNumberofPolynomialCoefficient();
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        oss << speciesName[s] << endl;
        for (int i = 0; i < nTemp; ++ i)
        {
            RDouble *coef = thermodynamicManager->GetPolynomialCoefficient(s, i);
            for (int j = 0; j < nCoef; ++ j)
            {
                oss << coef[j] << ", ";
            }
            oss << endl;
        }
    }
    oss << endl;
    int **forward = stoichiometricEquation->GetForwardReactionStoichiometricCoefficient();
    int **backward = stoichiometricEquation->GetBackwardReactionStoichiometricCoefficient();
    RDouble **thirdBody = stoichiometricEquation->GetReactionStoichiometricThirdBodyCoefficient();
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            oss << forward[r][s] << ", ";
        }
        oss << endl;
    }
    oss << endl;
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            oss << backward[r][s] << ", ";
        }
        oss << endl;
    }
    oss << endl;
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            oss << thirdBody[r][s] << ", ";
        }
        oss << endl;
    }
    oss << endl;

    RDouble *f1 = reactionRate->GetForwardReactionRateCoefficientA();
    RDouble *f2 = reactionRate->GetForwardReactionRateCoefficientB();
    RDouble *f3 = reactionRate->GetForwardReactionRateCoefficientC();
    RDouble *b1 = reactionRate->GetBackwardReactionRateCoefficientA();
    RDouble *b2 = reactionRate->GetBackwardReactionRateCoefficientB();
    RDouble *b3 = reactionRate->GetBackwardReactionRateCoefficientC();
    for (int r = 0; r < numberOfReaction; ++ r)
    {
        oss << f1[r] << ", " << f2[r] << ", " << f3[r] << ", " << b1[r] << ", " << b2[r] << ", " << b3[r] << endl;
    }
    oss << endl;

    WriteLogFile(oss);
}

void MixedGas::InitReferenceParameter()
{
    InitGasModel();
    InitNumberOfSpecies();

    int inflowParaType = GlobalDataBase::GetIntParaFromDB("inflowParaType");

    if (inflowParaType == NONDIMENSIONCONDITION)
    {
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

        referenceParameterDimensional.SetMachNumber(refMachNumber);
        referenceParameterDimensional.SetReynoldsNumber(refReNumber);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
    }
    else if (inflowParaType == FLIGHTCONDITION)
    {
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");

        referenceParameterDimensional.SetMachNumber(refMachNumber);

        SetAirInformationByDataBase();
    }
    else if (inflowParaType == EXPERIMENTCONDITION)
    {
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble refDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");

        referenceParameterDimensional.SetMachNumber(refMachNumber);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
        referenceParameterDimensional.SetPressure(refDimensionalPressure);

        SetAirInformationByExpData();
    }
    else if (inflowParaType == TEMPERATURE_DENSITY)    //! where temperature and density are fixed.
    {
        RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");

        referenceParameterDimensional.SetVelocity(refDimensionalVelocity);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
        referenceParameterDimensional.SetDensity(refDimensionalDensity);
    }
    else if (inflowParaType == TEMPERATURE_PRESSURE)    //! where temperature and pressure are fixed.
    {
        RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble refDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");

        referenceParameterDimensional.SetVelocity(refDimensionalVelocity);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
        referenceParameterDimensional.SetPressure(refDimensionalPressure);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("inflowParaType", inflowParaType);
    }

    RDouble prl = GlobalDataBase::GetDoubleParaFromDB("prl");
    RDouble prt = GlobalDataBase::GetDoubleParaFromDB("prt");

    RDouble oprl = 1.0 / prl;
    GlobalDataBase::UpdateData("oprl", &oprl, PHDOUBLE, 1);
    RDouble oprt = 1.0 / prt;
    GlobalDataBase::UpdateData("oprt", &oprt, PHDOUBLE, 1);

    this->oprl = oprl;
    this->oprt = oprt;
}

void MixedGas::ComputeCoefficientOfStateEquation()
{
    //! No matter whether chemical reactions exist or not.
    RDouble gama0 = referenceParameterDimensional.GetAverageSpecificHeatRatio();
    RDouble refMachNumber = referenceParameterDimensional.GetMachNumber();

    //! The variable is the non-dimensional value of 8.314/M,M denotes the non-dimensional molecular weight of the mixture gas.
    coefficientOfStateEquation = 1.0 / (gama0 * refMachNumber * refMachNumber); //! Universal gas constant.
    GlobalDataBase::UpdateData("coefficientOfStateEquation", &coefficientOfStateEquation, PHDOUBLE, 1);
}

void MixedGas::NormalizeAirInformation()
{
    RDouble refDimensionalPressure = referenceParameterDimensional.GetPressure();
    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();
    RDouble reference_average_generanl_gas_constant = referenceParameterDimensional.GetAverageGeneralGasConstant();

    RDouble refDimensionalDensity = refDimensionalPressure / reference_average_generanl_gas_constant / refDimensionalTemperature;

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
}

void MixedGas::InitCommonParameter()
{
    nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    nm    = GlobalDataBase::GetIntParaFromDB("nm");
    ntmodel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    nGasModel = GlobalDataBase::GetIntParaFromDB("nGasModel");
    nTEnergyModel = GlobalDataBase::GetIntParaFromDB("nTEnergyModel");
    parkVDPower = GlobalDataBase::GetDoubleParaFromDB("parkVDPower");
    chemicalRelaxCorf = GlobalDataBase::GetDoubleParaFromDB("chemicalRelaxCorf");
    nSpeciesLimit = GlobalDataBase::GetIntParaFromDB("nSpeciesLimit");
    nEnergyRecycle = GlobalDataBase::GetIntParaFromDB("nEnergyRecycle") * nchem;
    veTemperatureMin = GlobalDataBase::GetDoubleParaFromDB("veTemperatureMin");
    nChemcalSourceModified = GlobalDataBase::GetIntParaFromDB("nChemcalSourceModified");
    nChemcalSourceEsMethod = GlobalDataBase::GetIntParaFromDB("nChemcalSourceEsMethod");

    maxViscous = log(1.0);
    lnMaximum =  log(1.0e36);
    if (GlobalDataBase::IsExist("maxViscous", PHDOUBLE, 1))
    {
        maxViscous = log(GlobalDataBase::GetDoubleParaFromDB("maxViscous")*10.0);
    }

    trTemperatureMin = 10.0;
    if (GlobalDataBase::IsExist("trTemperatureMin", PHDOUBLE, 1))
    {
        trTemperatureMin = GlobalDataBase::GetDoubleParaFromDB("trTemperatureMin");
    }

    nDensityForWallMethod = 0;
    if (GlobalDataBase::IsExist("nDensityForWallMethod", PHINT, 1))
    {
        nDensityForWallMethod = GlobalDataBase::GetIntParaFromDB("nDensityForWallMethod");
    }

    if (GlobalDataBase::IsExist("nIsChemicalFreeze", PHINT, 1))
    {
        nIsChemicalFreeze = GlobalDataBase::GetIntParaFromDB("nIsChemicalFreeze");
    }
    if (GlobalDataBase::IsExist("nIsSuperCatalytic", PHINT, 1))
    {
        nIsSuperCatalytic = GlobalDataBase::GetIntParaFromDB("nIsSuperCatalytic");
        if (nIsSuperCatalytic != 0 && nIsSuperCatalytic != 1)
        {
            //cout<<"The parameter \"nIsSuperCatalytic\" does not equal to zero or one"<<endl;
            nIsSuperCatalytic = 1;
            GlobalDataBase::UpdateData("nIsSuperCatalytic", &nIsSuperCatalytic, PHINT, 1);
        }
    }

    if (parkVDPower < 0.0 || parkVDPower > 1.0)
    {
        parkVDPower = 0.6;
        GlobalDataBase::UpdateData("parkVDPower", &parkVDPower, PHDOUBLE, 1);
    }

    reynoldsReferenceLengthDimensional = 1.0;

    RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

    reynoldsReferenceLengthDimensional *= gridScaleFactor;

    GlobalDataBase::UpdateData("reynoldsReferenceLengthDimensional", &reynoldsReferenceLengthDimensional, PHDOUBLE, 1);

    mTT = 0;
    if (ntmodel == 3)
    {

    }
    else if (ntmodel == 2)
    {

    }
    else
    {
        mTV = 0;
        mTE = 0;
    }
}

void MixedGas::InitOtherParameter()
{
    int nolstress, nrokplus;
    nolstress = -1;
    nrokplus  = -1;
    GlobalDataBase::UpdateData("nolstress", &nolstress, PHINT, 1);
    GlobalDataBase::UpdateData("nrokplus" , &nrokplus , PHINT, 1);

    int outnstep = 0;
    GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);
    int innstep = 0;
    GlobalDataBase::UpdateData("innstep", &innstep, PHINT, 1);

    RDouble physicalTime = 0.0;
    GlobalDataBase::UpdateData("physicalTime", &physicalTime, PHDOUBLE, 1);

    int nIndexTT = 0;
    int nIndexTV = 0;
    int nIndexTE = 0;
    if (ntmodel == 2)
    {

    }
    else if (ntmodel == 3)
    {

    }
    else
    {
        nIndexTT = 0;
        nIndexTV = 0;
        nIndexTE = 0;
    }
    GlobalDataBase::UpdateData("indexOfTT", &nIndexTT, PHINT, 1);
    GlobalDataBase::UpdateData("indexOfTV", &nIndexTV, PHINT, 1);
    GlobalDataBase::UpdateData("indexOfTE", &nIndexTE, PHINT, 1);
}

void MixedGas::ComputeReferenceParameter()
{
    //! Convert to mass fraction, if the initial setting is other fractions convert other type fractions.
    ComputeReferenceInitMassFractionInformation();

    //! The first two variables is temperature independent, 
    //! depend on molecular weight and mass fraction of species if chemical reactions exist.
    ComputeReferenceMolecularInformation();

    ComputeReferenceGeneralGasConstant();

    //! The following five variables depend on temperature, 
    //! and mass fraction of species if chemical reactions exist.
    ComputeReferenceViscosityDimensional();

    ComputeReferenceSpecificHeatRatio();

    ComputeReferenceSoundVelocity();

    ComputeReferenceVelocity();

    ComputeCoefficientOfStateEquation();
    //! The above seven variables should be computed firstly.

    //! The variables in this function are computed by different methods according to "inflowParaType", 
    //! including density, pressure or reynolds.
    ComputeReferenceGasInformation();

    ComputeReferencePrimitive();
}

void MixedGas::InitGlobalParameterOfNSEquation()
{
    //! Obtain the variable of "reynoldsReferenceLengthDimensional" according to "gridScaleFactor".
    InitCommonParameter();

    InitReferenceParameter();

    ComputeReferenceParameter();

    PrintInflowConditionToWindow();

    ComputeProperReynoldsNumberForGrid();

    //! The five reference variables for the whole aerodynamic coefficients integration, which are from the file of "cfd_para.hypara", should be updated with "reynoldsReferenceLengthDimensional".
    ComputeOtherProperParameterForGrid();

    InitOtherParameter();

    //! Read parameters from file of "boundary_condition.hypara". Notice that the five reference variables for the component aerodynamic coefficients integration should be updated with "reynoldsReferenceLengthDimensional".
    GlobalBoundaryCondition::ReadGlobalBoundaryCondition();

    GlobalBoundaryCondition::SetGlobalBCByGlobalDataBase();

    GlobalBoundaryCondition::InitMassFlowBoundary();

    GlobalBoundaryCondition::SetBCDataBaseByGlobalBC();
}

void MixedGas::ComputeReferencePrimitive()
{
    RDouble refMachNumber = referenceParameterDimensional.GetMachNumber();
    RDouble reference_sound_velocity, reference_pressure;
    
    RDouble massoo_ref = 1.0;

    RDouble reference_density     = 1.0;
    RDouble reference_temperature = 1.0;

    reference_sound_velocity = 1.0 / refMachNumber;
    //! coefficientOfStateEquation = 1.0 / (gama0 * refMachNumber * refMachNumber);
    reference_pressure = coefficientOfStateEquation * 
                         reference_density * 
                         reference_temperature / massoo_ref;

    //! Set the global reference value of energy.
    RDouble refVelocity = referenceParameterDimensional.GetVelocity();
    RDouble refTemperature = referenceParameterDimensional.GetTemperature();
    referenceEnergy = refVelocity * refVelocity;
    referenceSpecificHeat = refTemperature / referenceEnergy;

    veTemperatureMinNonDimensional = veTemperatureMin / refTemperature;
    trTemperatureMinNonDimensional = trTemperatureMin / refTemperature;

    if (nchem == 1)
    {
        if (ntmodel > 1)
        {

        }
        else       //! One-temperature model.
        {
            reference_pressure = coefficientOfStateEquation * 
                                 reference_density * 
                                 reference_temperature / massoo_ref;
        }
    }

    RDouble attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    RDouble angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");

    RDouble attack   = attackd * PI / 180.0;
    RDouble sideslip = angleSlide * PI / 180.0;

    GlobalDataBase::UpdateData("attack", &attack, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("sideslip", &sideslip, PHDOUBLE, 1);

    using namespace IDX;
    RDouble *prim_inf = new RDouble[nl + nchem + ntmodel - 1];
    GlobalDataBase::UpdateDataPtr("prim_inf", prim_inf);

    prim_inf[IR] = reference_density;
    prim_inf[IU] = 1.0 * cos(attack) * cos(sideslip);
    prim_inf[IV] = 1.0 * sin(attack) * cos(sideslip);
    prim_inf[IW] = 1.0 *               sin(sideslip);
    prim_inf[IP] = reference_pressure;

    if (nchem > 0)
    {
        RDouble * mass_fraction_reference = this->molecularProperty->GetInitialMassFraction();

        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            prim_inf[nm + ispecies] = mass_fraction_reference[ispecies];
        }

        //! Set the vibration energy and electron energy.
        if (ntmodel == 2)    //! Tve = Tv = Te = Tinf.
        {

        }
        else if (ntmodel == 3)
        {

        }

        RDouble *catalyticMassFraction = new RDouble[numberOfSpecies];
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            catalyticMassFraction[s] = prim_inf[nm + s];
        }

        if (nIsSuperCatalytic == 1)
        {
            if (this->nGasModel == 0) //Gas of Earth.
            {
                GetEarthFullyCatalyticMassFraction(mass_fraction_reference, catalyticMassFraction);
            }
            else if (this->nGasModel == 1) //Gas of Mars.
            {
                GetMarsFullyCatalyticMassFraction(mass_fraction_reference, catalyticMassFraction);
            }
        }

        SetMaximumSpecies(mass_fraction_reference);

        GlobalDataBase::UpdateDataPtr("catalyticMassFraction", catalyticMassFraction);
    }
}

void MixedGas::CreateAllDataClass(int numberOfSpecies, int numberOfReaction)
{
    molecularProperty = new MolecularProperty();
    molecularProperty->Init(numberOfSpecies);

    reactionRate = new ReactionRate();
    reactionRate->Init(numberOfReaction);

    stoichiometricEquation = new StoichiometricEquation();
    stoichiometricEquation->Init(numberOfReaction, numberOfSpecies);

    CurveFits = new CurveFitsMethod();
    CurveFits->Init(numberOfSpecies);

    thermodynamicManager = new ThermodynamicManager();
    thermodynamicManager->Init(numberOfSpecies);
}

void MixedGas::Read(fstream &file)
{
    molecularProperty->Read(file);
    reactionRate->Read(file);
    thermodynamicManager->Read(file);
    CurveFits->Read(file);
    stoichiometricEquation->Read(file);
}

void MixedGas::ReadGasModel(fstream &file)
{
    string line, word;
    string separator = " =\t\r\n#$,;\"'";

    SkipLines(file, 2);

    getline(file, line);
    line = FindNextWord(line, word, separator);

    from_string< int >(numberOfSpecies, word, std::dec);
    line = FindNextWord(line, word, separator);
    from_string< int >(numberOfReaction, word, std::dec);

    CreateAllDataClass(numberOfSpecies, numberOfReaction);
    Read(file);
}

void MixedGas::GenerateGasModel(string chemicalModel)
{
    int nGasType = 0;
    int Ablation = GlobalDataBase::GetIntParaFromDB("nAblation");
    if (chemicalModel == "Gupta")
    {
        if (Ablation == 0)
        {
            nGasType = 0;
        }
        else
        {
            nGasType = 5;    //! Carbon material ablation.
        }
    }
    else if (chemicalModel == "Park")
    {
        if (Ablation == 0)
        {
            nGasType = 1;
        }
        else
        {
            nGasType = 6;    //! Carbon material ablation.
        }
    }
    else if (chemicalModel == "Dunn-Kang")    //! The Dunn-Kang Model.
    {
        if (Ablation == 0)
        {
            nGasType = 2;
        }
        else
        {
            nGasType = 7;    //! Carbon material ablation.
        }
    }
    else if (chemicalModel == "Mars-Park")    //! The Park Model for Mars.
    {
        nGasType = 3;
    }
    else if (chemicalModel == "Mars-McKenzie")    //! The McKenzie Model for Mars.
    {
        nGasType = 4;
    }
    else if (chemicalModel == "CombustionGas")    //! The Combustion chamber Gas Model.
    {
        nGasType = 8;
    }
    else if (chemicalModel == "Gas-Mixture")    //! The Mixed Gas Model.
    {
        nGasType = 9;
    }

    //Obtain the initial data.
    vector<string> namelist;
    vector<RDouble> initMassFraction;
    this->GetNameList(namelist, initMassFraction);

    //! Generate the chemical reaction model.
    ChemicalReactions myChemModel;
    ChemicalReactions *currentReaction = new ChemicalReactions();
    myChemModel.CreateBasicChemicalReactionModel(nGasType);
    myChemModel.CreateChemicalReactionModel(namelist, currentReaction);

    int ns = 0, nr = 0;
    currentReaction->GetChemicalModelInfo(ns, nr);
    if (ns <= 0 || nr <= 0)
    {
        cout << "Failed to create the gas model!" << endl;
        delete currentReaction;
        return;
    }

    //! The number of species and reaction are imported from the database.
    numberOfSpecies = ns;
    numberOfReaction = nr;
    //! Allocate memories.
    CreateAllDataClass(numberOfSpecies, numberOfReaction);

    //-------------------------------------------------------------------
    //! Generate the physical and chemical data of species.
    molecularProperty->GenerateData(namelist, initMassFraction);

    //! Generate the curve fitting data of chemical reaction rate.
    RDouble **reactionCoef = currentReaction->GetReactionRateCurveFittingData();
    reactionRate->Create(reactionCoef);

    //! Generate the curve fitting data of chemkin model for computing the specific heat and enthalpy in the one-temperature model.
    thermodynamicManager->GenerateData(namelist);

    //! Generate the curve fitting data of transport coefficients, such as viscosity.
    CurveFits->GenerateData(namelist);

    //! Generate the chemical reaction equations.
    int **reactionFlag = currentReaction->GetReactionType();
    int **forwardCoef = currentReaction->GetForwardStoichiometricCoefficient();
    int **backwardCoef = currentReaction->GetBackwardStoichiometricCoefficient();
    RDouble **collisionCoef = currentReaction->GetThirdBodyCollisionCoefficient();
    stoichiometricEquation->Create(forwardCoef, backwardCoef, collisionCoef, reactionFlag);

    delete currentReaction;
}

void MixedGas::GetNameList(vector<string> &namelist, vector<RDouble> &initMassFraction)
{
    string speciesName[30], line1, line2;
    RDouble massFraction[30];
    GlobalDataBase::GetData("speciesName", &line1, PHSTRING, 1);
    GlobalDataBase::GetData("initMassFraction", &line2, PHSTRING, 1);
    string word;
    string separator  = " =\t\r\n#$,;\"'";
    int nSpecies = 0;
    for (int s = 0; s < 30; ++ s)
    {
        line1 = FindNextWord(line1, word, separator);
        if (word == "")
        {
            break;
        }
        ++ nSpecies;
        speciesName[s] = word;

        line2 = FindNextWord(line2, word, separator);
        from_string<RDouble>(massFraction[s], word, std::dec);
    }
    namelist.clear();
    initMassFraction.clear();
    for (int s = 0; s < nSpecies; ++ s)
    {
        namelist.push_back(speciesName[s]);
        initMassFraction.push_back(massFraction[s]);
    }
}

void MixedGas::Read(DataContainer *cdata)
{
    molecularProperty->Read(cdata);
    reactionRate->Read(cdata);
    stoichiometricEquation->Read(cdata);
    CurveFits->Read(cdata);
    thermodynamicManager->Read(cdata);
}

void MixedGas::Write(DataContainer *cdata)
{
    molecularProperty->Write(cdata);
    reactionRate->Write(cdata);
    stoichiometricEquation->Write(cdata);
    CurveFits->Write(cdata);
    thermodynamicManager->Write(cdata);
}

void MixedGas::CompressData(DataContainer *&cdata)
{
    cdata->Append(&numberOfSpecies, sizeof(int));
    cdata->Append(&numberOfReaction, sizeof(int));

    Write(cdata);
}

void MixedGas::DecompressData(DataContainer *cdata)
{
    cdata->MoveToBegin();

    cdata->Read(&numberOfSpecies, sizeof(int));
    cdata->Read(&numberOfReaction, sizeof(int));

    CreateAllDataClass(numberOfSpecies, numberOfReaction);
    Read(cdata);
}

void MixedGas::GetTemperatureRangeIndex(const RDouble &temperature_dimensional, int &indexOfTemperatureRange)
{
    this->thermodynamicManager->GetTemperatureRangeIndex(temperature_dimensional, indexOfTemperatureRange);
}

void MixedGas::ComputeSpeciesConstantPressureSpecificHeatDimensional(RDouble temperature, RDouble *cpSpeciesDimensional)
{
    int ntype = 1;    //! 0 stands for applying the old method, 1 for the new method.

    if (ntype == 0)   //! The primary algorithm.
    {
        RDouble t1, t2, t3, t4;

        RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    
        t1 = temperature * referenceTemperature;
        t2 = t1 * t1;
        t3 = t1 * t2;
        t4 = t1 * t3;

        int itemperature;
        GetTemperatureRangeIndex(t1, itemperature);

        RDouble *oMolecularWeightDimensional = molecularProperty->GetMolecularWeightReciprocalDimensional();

        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            RDouble *coef = this->thermodynamicManager->GetPolynomialCoefficient(ispecies, itemperature);

            cpSpeciesDimensional[ispecies] = coef[0]
                                           + coef[1] * t1
                                           + coef[2] * t2
                                           + coef[3] * t3
                                           + coef[4] * t4;
            cpSpeciesDimensional[ispecies] *= rjmk * oMolecularWeightDimensional[ispecies];
        }
        return;
    }
    else      //! The new method that modified the coefficients in the temperature boundary range.
    {
        RDouble gasConstant = rjmk;      //! The universal gas constant equals to 8.31434 J/mol-K.
        //! The following variables are dimensional values.
         RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
         RDouble *speciesMass = this->molecularProperty->GetMolecularWeightDimensional(); //The mole weight kg/mol.

        int n = this->thermodynamicManager->GetNumberofPolynomialCoefficient();
        if (n < 1)
        {
            return;
        }

        //! Allocate the memory space.
        RDouble coef[7];
        RDouble T1, T2, T3, T4;
        T1 = temperature * referenceTemperature;
        T1 = MAX(T1, 50.0);
        T2 = T1 * T1;
        T3 = T1 * T2;
        T4 = T1 * T3;

        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            //! Apply the polynomial by linearly averaging the coefficients on the temperature boundary.
            GetLinearlyAveragedPolynomialCoefficients(ispecies, T1, coef);

            cpSpeciesDimensional[ispecies] = coef[0] + coef[1] * T1 + coef[2] * T2 + coef[3] * T3 + coef[4] * T4;
            cpSpeciesDimensional[ispecies] *= gasConstant / speciesMass[ispecies];
        }
    }
}

RDouble MixedGas::GetMixedGasMolecularWeight(RDouble *massFraction)
{
    RDouble mixedMass = 0.0;
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();
    for (int i = 0; i < numberOfSpecies; ++ i)
    {
        mixedMass += massFraction[i] / speciesMass[i];
    }
    mixedMass = 1.0 / mixedMass;

    //! Return the non-dimensional molecular weight of the mixture gas.
    return mixedMass;
}

RDouble MixedGas::GetMixedGasMolecularWeightReciprocalWithoutElectron(RDouble *massFraction, RDouble &ceDivideMe)
{
    RDouble *speciesMass = molecularProperty->GetMolecularWeight();
    //! Obtain the average molecular weight.
    RDouble mixedMass = GetMixedGasMolecularWeight(massFraction);
    ceDivideMe = 0.0;
    if (nElectronIndex >= 0)
    {
        ceDivideMe = massFraction[nElectronIndex] / speciesMass[nElectronIndex];
    }

    return (1.0 / mixedMass - ceDivideMe);
}

RDouble MixedGas::GetUniversalGasConstant()
{
    //! To obtain the non-dimensional universal gas constant called R = 8.314J/(mol*K).
    //RDouble gasConstant = rjmk;//! unit: J/(mol*K).
    //RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    //RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    //RDouble referenceMass = this->referenceParameterDimensional.GetAverageMolecularWeight();
    //gasConstant = gasConstant * referenceTemperature / (referenceMass * referenceVelocity * referenceVelocity);

    //// Return the non-dimensional universal gas constant.
    //return gasConstant;
    return coefficientOfStateEquation;
}

RDouble MixedGas::GetMixtureGasEnthalpy(RDouble *massFraction, RDouble temperature)
{
    //! Obtain the enthalpy of each species.
    this->ComputeSpeciesEnthalpy(temperature, speciesEnthalpy);
    //! Integral.
    RDouble mixedGasEnthalpy = 0.0;
    this->ComputeMixtureByMassFraction(massFraction, speciesEnthalpy, mixedGasEnthalpy);

    return mixedGasEnthalpy;
}

void MixedGas::ComputeConstantVolumeSpecificHeatByMassFraction(RDouble *massFraction, RDouble temperature, RDouble &mixedGasCv)
{
    RDouble mixedGasCp;       //! Specific heat at constant pressure.
    this->ComputeConstantPressureSpecificHeatByMassFraction(massFraction, temperature, mixedGasCp);

    //! Obtain universal gas constant.
    RDouble gasConstant = this->GetUniversalGasConstant();
    //! Obtain molecular weight of the mixture gas.
    RDouble mixedGasMass = this->GetMixedGasMolecularWeight(massFraction);

    //! Compute the constant volume specific heat of the mixture gas.
    mixedGasCv = mixedGasCp - gasConstant / mixedGasMass;
}

void MixedGas::ComputeSpeciesConstantVolumeSpecificHeat(RDouble temperature, RDouble *cvSpecies)
{
    //! To obtain the non-dimensional constant volume specific heat called cv of each species by fitting formula.
    //! temperature is the non-dimensional temperature of the mixture gas.
    //! cvSpecies is an array saving the non-dimensional constant volume specific heat of each species.

    //! Obtain the dimensional value of specific heat at constant volume.
    ComputeSpeciesConstantVolumeSpecificHeatDimensional(temperature, cvSpecies);

    //! Obtain the reference values.
    RDouble referenceVelocity = this->referenceParameterDimensional.GetVelocity();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceValue = referenceTemperature / (referenceVelocity * referenceVelocity);

    //! Translate the dimensional value to the non-dimensional value.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        cvSpecies[i] = cvSpecies[i] * referenceValue;
    }
}

void MixedGas::ComputeSpeciesConstantVolumeSpecificHeatDimensional(RDouble temperature, RDouble *cvSpeciesDimensional)
{
    //! To obtain the dimensional constant volume specific heat called cv[J/(kg*K)] of each species by fitting formula.
    //! temperature is the non-dimensional temperature of the mixture gas.
    //! cvSpeciesDimensional is an array saving the dimensional constant volume specific heat of each species.
    //! The following variables are the non-dimensional values.

    //! Obtain the dimensional value of constant pressure specific heat of each species.
    this->ComputeSpeciesConstantPressureSpecificHeatDimensional(temperature, cvSpeciesDimensional);

    RDouble gasConstant = rjmk;// dimensional universal gas constant.
    //! Obtain molecular weight of each species.
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();

    //! Compute the specific heat at constant volume.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        cvSpeciesDimensional[i] = cvSpeciesDimensional[i] - gasConstant / speciesMass[i];
    }
}

void MixedGas::ComputeConstantPressureSpecificHeatByMassFraction(RDouble *massFraction, RDouble temperature, RDouble &mixedcp)
{
    ComputeSpeciesConstantPressureSpecificHeat(temperature, speciesCp);
    ComputeMixtureByMassFraction(massFraction, speciesCp, mixedcp);
}

void MixedGas::ComputeSpeciesConstantPressureSpecificHeat(RDouble temperature, RDouble *cpSpecies)
{
    //! Obtain the dimensional value of species specific heat at constant pressure.
    this->ComputeSpeciesConstantPressureSpecificHeatDimensional(temperature, cpSpecies);

    //! Change to the non-dimensional values.
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        cpSpecies[ispecies] = cpSpecies[ispecies] * referenceSpecificHeat;
    }
}

void MixedGas::ComputeSpeciesEnthalpy(RDouble temperature, RDouble *speciesEnthaly)
{
    //! Obtain the dimensional values of species enthalpies.
    this->ComputeSpeciesEnthalpyDimensional(temperature, speciesEnthaly);

    //! Change to the non-dimensional values.
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        speciesEnthaly[s] = speciesEnthaly[s] / referenceEnergy;
    }
}

void MixedGas::ComputeSpeciesEnthalpyDimensional(RDouble temperature, RDouble *speciesEnthaly)
{
    int ntype = 1;     //! 0 stands for applying the old method, 1 for the new method.

    if (ntype == 0)    //! The primary algorithm.
    {
        RDouble t1, t2, t3, t4, ot1;

        RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();

        t1 = temperature * referenceTemperature;
        t2 = t1 * t1;
        t3 = t1 * t2;
        t4 = t1 * t3;
        ot1 = 1.0 / t1;

        int itemperature;
        GetTemperatureRangeIndex(t1, itemperature);

        RDouble *oMolecularWeightDimensional = molecularProperty->GetMolecularWeightReciprocalDimensional();

        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            RDouble *coef = this->thermodynamicManager->GetPolynomialCoefficient(ispecies, itemperature);
            speciesEnthaly[ispecies] = coef[0]            +
                                       coef[1] * t1 / 2.0 +
                                       coef[2] * t2 / 3.0 +
                                       coef[3] * t3 / 4.0 +
                                       coef[4] * t4 / 5.0 +
                                       coef[5] * ot1;
            speciesEnthaly[ispecies] *= rjmk * oMolecularWeightDimensional[ispecies] * t1;
        }
    }
    else      //! The new method that modified the coefficients in the temperature boundary range.
    {
        RDouble gasConstant = rjmk;     //! The universal gas constant equals to 8.31434 J/mol-K.
        //! The following variables are dimensional values.
        RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
        RDouble *speciesMass = this->molecularProperty->GetMolecularWeightDimensional();      //! The mole weight kg/mol.

        int n = this->thermodynamicManager->GetNumberofPolynomialCoefficient();
        if (n < 1)
        {
            return;
        }

        RDouble coef[7];
        RDouble T1, T2, T3, T4;
        T1 = temperature * referenceTemperature;
        T1 = MAX(T1, 50.0);
        T2 = T1 * T1;
        T3 = T1 * T2;
        T4 = T1 * T3;

        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            //! Apply the polynomial by linearly averaging the coefficients on the temperature boundary.
            GetLinearlyAveragedPolynomialCoefficients(ispecies, T1, coef);

            speciesEnthaly[ispecies] = coef[0] + coef[1] * T1/2.0 + coef[2] * T2/3.0 + coef[3] * T3/4.0 + coef[4] * T4/5.0 + coef[5]/T1;
            speciesEnthaly[ispecies] *= T1 * gasConstant / speciesMass[ispecies];
        }
    }
}

void MixedGas::ComputeSpeciesEnergy(RDouble temperature, RDouble *speciesEnergy)
{
}

void MixedGas::ComputeSpeciesEnergyDimensional(RDouble temperature, RDouble *speciesEnergy)
{
}

void MixedGas::ComputeMixtureByMassFraction(RDouble *massFraction, RDouble *mixtureOfSpecies, RDouble &mixture)
{
    mixture = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        mixture += massFraction[ispecies] * mixtureOfSpecies[ispecies];
    }
}

void MixedGas::ComputeMixtureByPrimitive(RDouble *prim, RDouble *mixtureOfSpecies, RDouble &mixture)
{
    mixture = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        mixture += prim[nm + ispecies] * mixtureOfSpecies[ispecies];
    }
}

void MixedGas::MassFractionConversion(RDouble *fractionOfSpecies)
{
    int nFraction = GlobalDataBase::GetIntParaFromDB("nFraction");

    if (nFraction == 1)
    {
        RDouble *molecularWeightDimensional  = molecularProperty->GetMolecularWeightDimensional();

        RDouble wiSum = 0.0;

        for (int ispecies = 0; ispecies < numberOfSpecies; ++ispecies)
        {
            wiSum += fractionOfSpecies[ispecies] * molecularWeightDimensional[ispecies];
        }

        RDouble owiSum = 1.0 / wiSum;

        for (int ispecies = 0; ispecies < numberOfSpecies; ++ispecies)
        {
            fractionOfSpecies[ispecies] = fractionOfSpecies[ispecies] * molecularWeightDimensional[ispecies] * owiSum;
        }
    }
}

void MixedGas::ComputeMoleFractionByMassFraction(RDouble *massFractionOfSpecies, RDouble *moleFractionOfSpecies)
{
    RDouble *oMolecularWeight = molecularProperty->GetMolecularWeightReciprocal();     //! The reciprocal value of mass fraction.

    RDouble xiSum, oxiSum;

    xiSum = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        xiSum += massFractionOfSpecies[ispecies] * oMolecularWeight[ispecies];
    }
    oxiSum = 1.0 / xiSum;

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
       moleFractionOfSpecies[ispecies] = massFractionOfSpecies[ispecies] * oMolecularWeight[ispecies] * oxiSum;
    }
}

void MixedGas::ComputeMixtureCoefficientByWilkeFormula(RDouble *moleFractionOfSpecies, RDouble *chem_var, RDouble *chem_phi)
{
//! Arguments: moleFractionOfSpecies saves the mole fraction of the species,chem_var stores the viscosity of the given species,
//! return the value phi stored in the variable chem_phi.

    RDouble *molecularWeight = molecularProperty->GetMolecularWeight();

    RDouble tmp, tmp1, tmp2;

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        chem_phi[ispecies] = 0.0;
        for (int jspecies = 0; jspecies < numberOfSpecies; ++ jspecies)
        {
            tmp = sqrt(molecularWeight[jspecies] / molecularWeight[ispecies]);
            tmp1 = 1.0 + sqrt(chem_var[ispecies] / chem_var[jspecies]) * sqrt(tmp);
            tmp2 = sqrt(8.0 * (1.0 + molecularWeight[ispecies] / molecularWeight[jspecies]));
            chem_phi[ispecies] += moleFractionOfSpecies[jspecies] * SQR(tmp1) / tmp2;
        }
    }
}

void MixedGas::ComputeMixtureByWilkeFormula(RDouble *moleFractionOfSpecies, RDouble *mixtureOfSpecies, RDouble *phiOfSpecies, RDouble &mixture)
{
    mixture = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        mixture += moleFractionOfSpecies[ispecies] * mixtureOfSpecies[ispecies] / phiOfSpecies[ispecies];
    }
}

void MixedGas::ComputeTemperature(RDouble *primitiveVariables, RDouble &temperature)
{
}

void MixedGas::GetTemperature(RDouble *primitiveVariables, RDouble &transRotationTemperature, RDouble &vibrationTemperature, RDouble &electronTemperature)
{
    using namespace IDX;
    RDouble density, pressure, gasConstant;
    density = primitiveVariables[IR];
    pressure = primitiveVariables[IP];
    gasConstant = GetUniversalGasConstant();;

    if (nchem == 0)     //! Perfect gas.
    {
        transRotationTemperature = pressure / (density * gasConstant);
    }
    else                //! Chemical non-equilibrium flow.
    {
        //! Compute the mass fraction of the last species.
        NormalizePrimitive(primitiveVariables);

        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            massFractions[i] = primitiveVariables[nm + i];
        }

        RDouble massReciprocal = 0.0, ceDivideMe = 0.0;
        massReciprocal = GetMixedGasMolecularWeightReciprocalWithoutElectron(massFractions, ceDivideMe);

        if (ntmodel == 1) //! The default setting is One-Temperature model.
        {
            transRotationTemperature = pressure / (density * gasConstant * (massReciprocal + ceDivideMe));
        }    
        else if (ntmodel == 2)      //! Two-Temperature model.
        {

        }
        else  //! Three-Temperature model.
        {

        }
    }
}

void MixedGas::GetTemperatureR(RDouble *primitiveVariables, RDouble &transRotationTemperature, RDouble &vibrationTemperature, RDouble &electronTemperature)
{
    using namespace IDX;
    RDouble density, pressure, gasConstant;
    density = primitiveVariables[IR];
    pressure = primitiveVariables[IP];
    gasConstant = GetUniversalGasConstant();

   // Compute the mass fraction of the last species.
    NormalizePrimitive(primitiveVariables);

    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        massFractions[i] = primitiveVariables[nm + i];
    }

    RDouble massReciprocal = 0.0, ceDivideMe = 0.0;
    massReciprocal = GetMixedGasMolecularWeightReciprocalWithoutElectron(massFractions, ceDivideMe);

    if (ntmodel == 1)
    {
        transRotationTemperature = pressure / (density * gasConstant * (massReciprocal + ceDivideMe));
    }
    else if (ntmodel == 2)
    {

    }
    else
    {

    }

    return;
}

RDouble MixedGas::GetTemperatureMin()
{
    RDouble refDimensionalTemperature = this->referenceParameterDimensional.GetTemperature();
    return this->thermodynamicManager->GetTemperatureMin() / refDimensionalTemperature;
}

RDouble MixedGas::GetTemperatureMax()
{
    RDouble refDimensionalTemperature = this->referenceParameterDimensional.GetTemperature();
    return this->thermodynamicManager->GetTemperatureMax() / refDimensionalTemperature;
}

string * MixedGas::GetNameOfSpecies() const 
{
    if (!molecularProperty) return 0;
    return molecularProperty->GetNameOfSpecies(); 
};

void MixedGas::GetElectronMolecularWeightAndMassFraction(RDouble *massFraction, RDouble &ce, RDouble &Me)
{
}

void MixedGas::GetPressure(RDouble *primitiveVars, RDouble gama, RDouble internalEnergy, RDouble temperature, RDouble &pressure)
{
    if (nchem != 1)      //! Perfect gas.
    {
        pressure = (gama - one) * internalEnergy;
        return;
    }

    using namespace IDX;
    RDouble density, specificInternalEnergy, gasConstant;
    gasConstant = this->GetUniversalGasConstant();
    density = primitiveVars[IR];
    specificInternalEnergy = internalEnergy / density;

    //! Obtain the mass fractions of species.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        massFractions[i] = primitiveVars[nm + i];
    }

    RDouble transRotationTemperature, massReciprocal, ceDivideMe, Tk;
    massReciprocal = GetMixedGasMolecularWeightReciprocalWithoutElectron(massFractions, ceDivideMe);
    Tk = temperature;
    if (ntmodel == 2)     //! Two-Temperature model.
    {

    }
    else if (ntmodel == 3)     //! Three-Temperature model.
    {

    }
    else       //! The default is One-Temperature model.
    {
        transRotationTemperature = ComputeNonequilibrumTemperature(massFractions, specificInternalEnergy, Tk);
        pressure = density * gasConstant * transRotationTemperature * (massReciprocal + ceDivideMe);
    }
}

void MixedGas::GetTemperatureAndPressure(RDouble *primitiveVars, RDouble internalEnergy, RDouble *temperature, RDouble &pressure)
{
    using namespace IDX;
    RDouble density, specificInternalEnergy, gasConstant;
    gasConstant = this->GetUniversalGasConstant();
    density = primitiveVars[IR];
    specificInternalEnergy = internalEnergy / density;

    //! Obtain the mass fractions of species.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        massFractions[i] = primitiveVars[nm + i];
    }

    RDouble transRotationTemperature, massReciprocal, ceDivideMe, Tk;
    massReciprocal = GetMixedGasMolecularWeightReciprocalWithoutElectron(massFractions, ceDivideMe);
    
    if (ntmodel == 2)     //! Two-Temperature model.
    {

    }
    else if (ntmodel == 3)     //! Three-Temperature model.
    {

    }
    else       //! The default is One-Temperature model.
    {
        Tk = temperature[ITT];
        transRotationTemperature = ComputeNonequilibrumTemperature(massFractions, specificInternalEnergy, Tk);
        pressure = density * gasConstant * transRotationTemperature * (massReciprocal + ceDivideMe);

        //Save the temperature of the next time-advancing step.
        temperature[ITT] = transRotationTemperature;
        temperature[ITV] = transRotationTemperature;
        temperature[ITE] = transRotationTemperature;
    }
}

void MixedGas::NormalizePrimitive(RDouble *primitiveVariables)
{
    RDouble fSum = 0.0;
    for (int s = 0; s < this->numberOfSpecies - 1; ++s)
    {
        fSum += primitiveVariables[nm + s];
    }

    //! Compute mass fraction of the last species.
    primitiveVariables[nm + numberOfSpecies - 1] = MAX(0.0, 1.0 - fSum);
    fSum += primitiveVariables[nm + numberOfSpecies - 1];

    //! Normalize the mass fractions.
    for (int s = 0; s < this->numberOfSpecies; ++ s)
    {
        primitiveVariables[nm + s] /= fSum;
    }
}

void MixedGas::ComputeMolecularWeightReciprocal(RDouble *primitiveVariables, RDouble &massReciprocal)
{
    RDouble *oMolecularWeight = molecularProperty->GetMolecularWeightReciprocal();

    massReciprocal = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        massReciprocal += primitiveVariables[ispecies + nm] * oMolecularWeight[ispecies];
    }
}

void MixedGas::ComputeMolecularWeightReciprocalDimensional(RDouble *prim, RDouble &omavDimensional)
{
    RDouble *oMolecularWeightDimensional = molecularProperty->GetMolecularWeightReciprocalDimensional();

    omavDimensional = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        omavDimensional += prim[ispecies + nm] * oMolecularWeightDimensional[ispecies];
    }
}

void MixedGas::ComputeOneTemperatureModelSpeciesViscosity(RDouble temperature, RDouble *viscosityOfSpecies)
{
    //! Obtain the dimensional values of each species viscosity.
    ComputeOneTemperatureModelSpeciesViscosityDimensional(temperature, speciesViscosity);

    RDouble referenceViscosity = this->referenceParameterDimensional.GetViscosity();

    //! Translate to the non-dimensional value.
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        viscosityOfSpecies[ispecies] = speciesViscosity[ispecies] / referenceViscosity;
    }
}

void MixedGas::ComputeOneTemperatureModelSpeciesViscosityDimensional(RDouble temperature, RDouble *viscosityOfSpeciesDimensional)
{
    //! To obtain the dimensional viscosity of each species according to the Blotter fitting formula.
    //! temperature is the non-dimensional temperature of the mixture gas.
    //! viscosityOfSpeciesDimensional is an array saving the dimensional viscosity of each species [N*s/m2].

    int nType = 0;     //! 0 stands for selecting the Blotter curve fits mode, 1 is for Sutherland relation.
    nType = GlobalDataBase::GetIntParaFromDB("nSutherland");
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble T = temperature * referenceTemperature;

    if (nType == 0)     //! Compute species viscosity using the Blotter curve fits for single-temperature model.
    {
        RDouble *blotterAi = CurveFits->GetBlotterCurveFittingCoefficientA();
        RDouble *blotterBi = CurveFits->GetBlotterCurveFittingCoefficientB();
        RDouble *blotterCi = CurveFits->GetBlotterCurveFittingCoefficientC();
        RDouble *blotterDi = CurveFits->GetBlotterCurveFittingCoefficientD();
        RDouble *blotterEi = CurveFits->GetBlotterCurveFittingCoefficientE();
   
        RDouble lnt1, lnt2, lnt3, lnt4, tmp;
        lnt1 = log(T);
        lnt2 = lnt1 * lnt1;
        lnt3 = lnt1 * lnt2;
        lnt4 = lnt1 * lnt3;

        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            tmp = blotterAi[ispecies] * lnt4
                + blotterBi[ispecies] * lnt3
                + blotterCi[ispecies] * lnt2
                + blotterDi[ispecies] * lnt1
                + blotterEi[ispecies];

            //! The Blotter curve fits.
            viscosityOfSpeciesDimensional[ispecies] = 0.10 * exp(tmp) + SMALL;
        }
    }
    else       //! Compute species viscosity using the Sutherland relation.
    { 
        RDouble miu0, T0, C;
        miu0 = 1.7894E-5;     //! Unit: kg/(m*s).
        T0 = 288.15;          //! Unit: K.
        C  = 110.4;           //! Unit: K.
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            viscosityOfSpeciesDimensional[ispecies] = miu0 * pow(T / T0, 1.5) * (T0 + C) / (T + C);
        }
    }
}

void MixedGas::ComputeSpeciesEnthalpy(RDouble temperature, RDouble *speciesCv, RDouble *speciesE, RDouble *speciesH)
{
    ComputeSpeciesConstantVolumeSpecificHeat(temperature, speciesCv);

    ComputeSpeciesConstantPressureSpecificHeat(temperature, speciesH);

    RDouble gasConstant = this->GetCoefficientOfStateEquation();
    RDouble *speciesMass = molecularProperty->GetMolecularWeight();
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        speciesE[s] = speciesH[s] - gasConstant / speciesMass[s];
    }
}

void MixedGas::ComputeSpeciesConstantPressureSpecificHeatByEuckenFormula(RDouble *cpSpecies)
{
}

void MixedGas::ComputeEnthalpyByPrimitive(RDouble *primitiveVariables, RDouble gama, RDouble &enthalpy, RDouble *temperatures)
{
    using namespace IDX;

    RDouble &density  = primitiveVariables[IR];
    RDouble &pressure = primitiveVariables[IP];
    int nIdealState = GlobalDataBase::GetIntParaFromDB("nIdealState");
    if (nchem == 0 || nIdealState == 1)
    //if (nchem == 0)
    {
        enthalpy = (gama / (gama - one)) * (pressure / density);
    }
    else
    {
        enthalpy = this->ComputeMixedGasEnthalpy(primitiveVariables, temperatures);
    }
}

void MixedGas::ComputeInternalEnergy(RDouble *primitiveVariables, RDouble gama, RDouble vibrationTemperature, RDouble electronTemperature, RDouble &totalEnergy)
{
    using namespace IDX;
    RDouble temperature = 0.0, squareVelocity = 0.0, massReciprocal = 0.0, staticEnthalpy = 0.0;

    RDouble &um = primitiveVariables[IU];
    RDouble &vm = primitiveVariables[IV];
    RDouble &wm = primitiveVariables[IW];
    RDouble &density = primitiveVariables[IR];
    RDouble &pressure = primitiveVariables[IP];
    squareVelocity = um * um + vm * vm + wm * wm;
    int nIdealState = GlobalDataBase::GetIntParaFromDB("nIdealState");
    if (nchem == 0 || nIdealState == 1)//! Perfect gas.
    //if (nchem == 0)
    {
        totalEnergy = (pressure / density) / (gama - one) + half * squareVelocity;
    }
    else
    {
        if (ntmodel == 1)     //! Single-Temperature model.
        {
            ComputeMolecularWeightReciprocal(primitiveVariables, massReciprocal);
            temperature = pressure / (coefficientOfStateEquation * density * massReciprocal);
            ComputeSpeciesEnthalpy(temperature, speciesEnthalpy);
            ComputeMixtureByPrimitive(primitiveVariables, speciesEnthalpy, staticEnthalpy);
            totalEnergy = staticEnthalpy - pressure / density + half * squareVelocity;
        }
        else      //! Multi-temperature model.
        {

        }
    }
}

void MixedGas::ComputeDensityDiffusionAndKCP(RDouble temperature, RDouble *primface, RDouble mul, RDouble mut, RDouble *rho_ds_face, RDouble *hintSpeciesOfFace, RDouble &kcp)
{
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        massFractions[ispecies] = primface[nm + ispecies];
    }

    ComputeMoleFractionByMassFraction(massFractions, moleFractions);

    RDouble *oscl = GetLaminarSchmidtNumberReciprocal();
    RDouble *osct = GetTurbulentSchmidtNumberReciprocal();

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        RDouble tmp = 1.0 - moleFractions[ispecies];
        if (tmp <= 0.0)
        {
            rho_ds_face[ispecies] = mul * oscl[ispecies] + mut * osct[ispecies];
        }
        else
        {
            rho_ds_face[ispecies] = (mul * oscl[ispecies] + mut * osct[ispecies]) * (1.0 - massFractions[ispecies]) / tmp;
        }
    }


    ComputeSpeciesEnthalpy(temperature, speciesEnthalpy);

    RDouble cp = -1.0;
    ComputeConstantPressureSpecificHeatByMassFraction(massFractions, temperature, cp);
    kcp = (mul * oprl + mut * oprt) * cp;

    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        hintSpeciesOfFace[ispecies] = speciesEnthalpy[ispecies];
    }
}

void MixedGas::ComputeDHAndTotalEnthalpy(RDouble *primitiveVariables, RDouble &gama, RDouble *deltaQ, RDouble &deltaEnthalpy, RDouble &totalEnthalpy)
{
    using namespace IDX;

    RDouble squareVelocity, ae, af;
    RDouble enthalpy;
    //! Obtain the primitive variables.
    RDouble &um = primitiveVariables[IU];
    RDouble &vm = primitiveVariables[IV];
    RDouble &wm = primitiveVariables[IW];
    RDouble &density = primitiveVariables[IR];
    RDouble &pressure = primitiveVariables[IP];

    squareVelocity = um * um + vm * vm + wm * wm;
    ae = gama - one;
    af = half * ae * squareVelocity;
    deltaEnthalpy = -ae * (um * deltaQ[IRU] + vm * deltaQ[IRV] + wm * deltaQ[IRW] - deltaQ[IRE]);

    int nIdealState = GlobalDataBase::GetIntParaFromDB("nIdealState");
    if (nchem == 0 || nIdealState ==1)
    //if (nchem == 0)
    {
        //! h = Cp*T = (gamma*Rg/gamma-1)*T = (gamma/gamma-1)*(p/rho).
        enthalpy = (gama / (gama - one)) * (pressure / density);
    }
    else
    {
        //! The obtained value is the reciprocal value of molecular weight.
        RDouble *oMolecularWeight = molecularProperty->GetMolecularWeightReciprocal();

        RDouble temperature, massReciprocal, as;
        ComputeMolecularWeightReciprocal(primitiveVariables, massReciprocal);
        temperature = pressure / (coefficientOfStateEquation * density * massReciprocal);
        ComputeSpeciesEnthalpy(temperature, speciesEnthalpy);
        ComputeMixtureByPrimitive(primitiveVariables, speciesEnthalpy, enthalpy);

        for (int ispecies = 0; ispecies < numberOfSpecies - 1; ++ ispecies)
        {
            as = -ae * (speciesEnthalpy[ispecies] - speciesEnthalpy[numberOfSpecies - 1])
               + coefficientOfStateEquation * gama * temperature * (oMolecularWeight[ispecies] - oMolecularWeight[numberOfSpecies - 1]);
            deltaEnthalpy += as * deltaQ[nm + ispecies];
        }
        af += coefficientOfStateEquation * gama * temperature * oMolecularWeight[numberOfSpecies - 1] - ae * speciesEnthalpy[numberOfSpecies - 1];
    }

    deltaEnthalpy += af * deltaQ[IR];
    // ht = h+v*v/2
    totalEnthalpy = enthalpy + half * squareVelocity;
}

void MixedGas::ComputeTotalEnthalpyAndDH(RDouble *primitiveVariables, RDouble &gama, RDouble *deltaQ, RDouble &totalEnthalpy, RDouble &deltaEnthalpy)
{
    //! This function is used for perfect gas and single temperature model.
    //! To compute the total enthalpy and the variable dh in the function MXDQ_STD(). dh=b2, b2 denotes the coefficient \n
    //! of the vector M*dQ which can be referred to the forumla (A.7) and (A.8) in the appendix A of the PHengLEI Theory manual.
    using namespace IDX;

    //! Obtain the primitive variables.
    RDouble &um = primitiveVariables[IU];
    RDouble &vm = primitiveVariables[IV];
    RDouble &wm = primitiveVariables[IW];
    RDouble &density = primitiveVariables[IR];
    RDouble &pressure = primitiveVariables[IP];

    //! c1 is product coefficient of the variable dq[IR],c2 = R*T/Mns - (gama-1)*ens,gama1 = gama-1.
    RDouble squareVelocity, c1, c2, gama1;
    squareVelocity = um * um + vm * vm + wm * wm;
    gama1 = gama - 1.0;
    c1 = half * gama1 * squareVelocity;
    deltaEnthalpy = -gama1 * (um * deltaQ[IRU] + vm * deltaQ[IRV] + wm * deltaQ[IRW] - deltaQ[IRE]);

    //! enthalpy is static enthalpy.
    RDouble enthalpy;
    int nIdealState = GlobalDataBase::GetIntParaFromDB("nIdealState");
    if (nchem == 0 || nIdealState == 1)     //! Compute the static enthalpy.
    //if (nchem == 0)
    {
        // h = cp*T = [gama/(gama-1)]*(R/M)*T = [gama/(gama-1)]*(p/rho).
        enthalpy = (gama / gama1) * (pressure / density);
    }
    else      //! Compute the static enthalpy of the multi-species model.
    {
        //! Obtain the non-dimensional molecular weight of each species saving in the array Ms.
        RDouble *speciesMass = molecularProperty->GetMolecularWeight();
        RDouble gasConstant, temperature, mass;
        for (int i = 0; i < numberOfSpecies; ++ i)
        {
            massFractions[i] = primitiveVariables[nm + i];
        }

        gasConstant = GetUniversalGasConstant();
        mass = GetMixedGasMolecularWeight(massFractions);
        //! Compute temperature using the state equation.
        temperature = (pressure / density) * (mass / gasConstant);

        //! Obtain constant volume specific heat of each species.
        ComputeSpeciesConstantVolumeSpecificHeat(temperature, this->speciesCv);

        // c2 = R*T/Mns - (gama-1)*ens
        c2 = gasConstant * temperature / speciesMass[numberOfSpecies - 1] - gama1 * speciesCv[numberOfSpecies - 1] * temperature;
        c1 += c2;     //! Added to the product coefficient of the variable dq[IR].

        //! fSum saves the summation of each term.
        RDouble fSum = 0.0;
        for (int i = 0; i < numberOfSpecies - 1; ++ i)
        {
            fSum += deltaQ[nm + i];
        }
        deltaEnthalpy += -c2 * fSum;

        fSum = 0.0;
        for (int i = 0; i < numberOfSpecies - 1; ++ i)
        {
            fSum += deltaQ[nm + i] / speciesMass[i];
        }
        deltaEnthalpy += gasConstant * temperature * fSum;

        fSum = 0.0;
        for (int i = 0; i < numberOfSpecies - 1; ++ i)
        {
            fSum += deltaQ[nm + i] * speciesCv[i] * temperature;
        }
        deltaEnthalpy += -gama1 * fSum;

        this->ComputeSpeciesEnthalpy(temperature, speciesEnthalpy);
        enthalpy = 0.0;     //! Compute the static enthalpy.
        for (int i = 0; i < numberOfSpecies; ++ i)
        {
            enthalpy += massFractions[i] * speciesEnthalpy[i];
        }
    }

    deltaEnthalpy += c1 * deltaQ[IR];
    // h0 = h+v*v/2, total enthalpy equals to static enthalpy plus half of the squared velocity.
    totalEnthalpy = enthalpy + half * squareVelocity;
}

void MixedGas::ComputeViscosityByPrimitive(RDouble *primitiveVariables, RDouble &transRotationTemperature, RDouble &electronTemperature, RDouble &suthTemperature, RDouble &viscosity)
{
    using namespace IDX;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        massFractions[ispecies] = primitiveVariables[nm + ispecies];
    }

    ComputeViscosity(primitiveVariables[IR], massFractions, transRotationTemperature, electronTemperature, suthTemperature, viscosity);
}

void MixedGas::ComputeTransportCoefficients(RDouble *primitiveVariables, RDouble *temperature, RDouble visTurbulence, RDouble &viscosity, RDouble *heatConductivity, RDouble *rhoDs)
{
    using namespace IDX;
    //RDouble rho = primitiveVariables[IR];
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        massFractions[s] = primitiveVariables[nm + s];
    }
    // Compute mole fractions using the mass fraction.
    ComputeMoleFractionByMassFraction(massFractions, moleFractions);

    RDouble *speciesMass = molecularProperty->GetMolecularWeight();

    if (ntmodel > 1) //Multi-temperature model.
    {

    }
    else//One-temperature model.
    {
        ComputeOneTemperatureModelSpeciesViscosity(temperature[ITT], speciesViscosity);
        GetPartitionFunctionValues(speciesViscosity, speciesMass, moleFractions, speciesWeight);

        ComputeSpeciesConstantPressureSpecificHeat(temperature[ITT], speciesCp);
        ComputeSpeciesHeatConductivityByEuckenFormula(speciesViscosity, speciesCp, speciesConductivity);

        // Compute heat conductivity.
        ComputeMixtureByWilkeFormula(moleFractions, speciesConductivity, speciesWeight, heatConductivity[ITT]);
    }

    // Compute viscosity.
    ComputeMixtureByWilkeFormula(moleFractions, speciesViscosity, speciesWeight, viscosity);

    // Compute mass diffusion coefficients.
    ComputeSpeciesMassDiffusionCoefficient(massFractions, viscosity, visTurbulence, rhoDs);

}

void MixedGas::ComputeTransportCoefficientsR(RDouble *massFractions, RDouble *temperature, RDouble *speciesCps, RDouble *speciesCvvs, RDouble *speciesCves, RDouble visTurbulence, RDouble &viscosity, RDouble *heatConductivity, RDouble *rhoDs)
{
    using namespace IDX;

    RDouble *oMolecularWeight = molecularProperty->GetMolecularWeightReciprocal();     //! The reciprocal value of mass fraction.
    RDouble *speciesMass = molecularProperty->GetMolecularWeight();
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referenceViscosity = this->referenceParameterDimensional.GetViscosity();
    int *ionType = this->molecularProperty->GetIonTypeOfSpecies();
    RDouble *oSchmidtLaminar = GetLaminarSchmidtNumberReciprocal();
    RDouble *oSchmidtTurbulence = GetTurbulentSchmidtNumberReciprocal();
    RDouble gasConstant = this->GetUniversalGasConstant(); 

    RDouble *blotterAi = CurveFits->GetBlotterCurveFittingCoefficientA();
    RDouble *blotterBi = CurveFits->GetBlotterCurveFittingCoefficientB();
    RDouble *blotterCi = CurveFits->GetBlotterCurveFittingCoefficientC();
    RDouble *blotterDi = CurveFits->GetBlotterCurveFittingCoefficientD();
    RDouble *blotterEi = CurveFits->GetBlotterCurveFittingCoefficientE();

    RDouble Ttr = temperature[mTT] * referenceTemperature;

    RDouble xiSum, oxiSum;
    xiSum = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        moleFractions[ispecies] = massFractions[ispecies] * oMolecularWeight[ispecies];
        xiSum += moleFractions[ispecies];
    }

    oxiSum = 1.0 / xiSum;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
       moleFractions[ispecies] *= oxiSum;
    }

    RDouble lnT = log(static_cast<RealData>(Ttr));
    for (int ispecies = 0; ispecies <= nLastSpeciesIndex; ++ ispecies)
    {
        oxiSum = min(maxViscous, (((blotterAi[ispecies] * lnT + blotterBi[ispecies]) * lnT + blotterCi[ispecies]) * lnT + blotterDi[ispecies]) * lnT + blotterEi[ispecies]);
        speciesViscosity[ispecies] = 0.10 * exp(static_cast<RealData>(oxiSum)) / referenceViscosity + SMALL;
    }

    if (nElectronIndex >= 0)
    {
        RDouble Te = temperature[mTE] * referenceTemperature;
        RDouble lnTe = log(static_cast<RealData>(Te));
        oxiSum = min(maxViscous, (((blotterAi[nElectronIndex] * lnTe + blotterBi[nElectronIndex]) * lnTe + blotterCi[nElectronIndex]) * lnTe + blotterDi[nElectronIndex]) * lnTe + blotterEi[nElectronIndex]);
        speciesViscosity[nElectronIndex] = 0.10 * exp(static_cast<RealData>(oxiSum)) / referenceViscosity + SMALL;
    }
    
    RDouble tmp0, tmp1, phiTmp;
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        phiTmp = 0.0;
        for (int i = 0; i < numberOfSpecies; ++ i)
        {
            tmp0 = speciesMass[i] / speciesMass[s];
            tmp1 = 1.0 + sqrt(static_cast<RealData>(speciesViscosity[s] / speciesViscosity[i])) * sqrt((sqrt(static_cast<RealData>(tmp0))));
            phiTmp += moleFractions[i] * tmp1 * tmp1 / sqrt(static_cast<RealData>(8.0 + 8.0 / tmp0));
        }
        speciesWeight[s] = phiTmp;
    }

    //! Compute the viscosity using Eucken formula.
    viscosity = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        speciesConductivity[ispecies] = moleFractions[ispecies] * speciesViscosity[ispecies] / speciesWeight[ispecies];
        viscosity += speciesConductivity[ispecies];
    }

    //! Compute mass diffusion coefficient for Species
    tmp1 = 0;
    phiTmp = 0;
    for (int i = 0; i <= this->nLastSpeciesIndex; ++ i)
    {
        if (1.0 - moleFractions[i] <= 1.0e-20)
        {
            rhoDs[i] = 0.0;
        }
        else
        {
            rhoDs[i] = (1.0 - massFractions[i]) / (1.0 - moleFractions[i]) * (viscosity * oSchmidtLaminar[i] + visTurbulence * oSchmidtTurbulence[i]);
        }

        if (ionType[i] != 0)
        {
            tmp0  = massFractions[i] * oMolecularWeight[i];
            tmp1 += rhoDs[i] * tmp0;
            phiTmp += tmp0;
        }
    }

    //! Compute the mass diffusion coefficient of electron.
    if (nElectronIndex >= 0)
    {
        rhoDs[nElectronIndex] = tmp1 / MAX(phiTmp, 1.0e-20);
    }

    //! Compute the heat conductivity using Eucken formula.
    heatConductivity[ITT] = 0.0;
    for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
    {
        heatConductivity[ITT] += speciesConductivity[ispecies] * (speciesCps[ispecies] + 1.25 * gasConstant * oMolecularWeight[ispecies]);
    }

    heatConductivity[ITV] = 0.0;
    heatConductivity[ITE] = 0.0;
    if (ntmodel == 2)
    {

    }
    else if (ntmodel == 3)
    {

    }
}

void MixedGas::ComputeViscosity(RDouble density, RDouble *massFractionOfSpecies, RDouble &transRotationTemperature, RDouble &electronTemperature, RDouble &suthTemperature, RDouble &viscosity)
{
    if (nchem == 0)     //! Perfect gas.
    {
        viscosity = transRotationTemperature * sqrt(transRotationTemperature) * (1.0 + suthTemperature) / (transRotationTemperature + suthTemperature);
    }
    else      //! Chemical non-equilibrium flow.
    {
        RDouble *speciesMass = molecularProperty->GetMolecularWeight();
        ComputeMoleFractionByMassFraction(massFractionOfSpecies, moleFractions);

        if (ntmodel > 1)     //! Multi-temperature model.
        {

        }
        else     //! One-temperature model.
        {
            ComputeOneTemperatureModelSpeciesViscosity(transRotationTemperature, speciesViscosity);
        }

        GetPartitionFunctionValues(speciesViscosity, speciesMass, moleFractions, speciesWeight);

        viscosity = GetMixedGasViscosityWithWilkeFormula(moleFractions, speciesViscosity, speciesWeight);
    }
}

void MixedGas::ComputeSpeciesHeatConductivityByEuckenFormula(RDouble *viscosity, RDouble *cps, RDouble *conductivity)
{
    //! Obtain the universal gas constant.
    RDouble gasConstant = this->GetUniversalGasConstant(); 
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();

    //! Compute the heat conductivity using Eucken formula.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        conductivity[i] = viscosity[i] * cps[i];
        conductivity[i] += 1.25 * viscosity[i] * gasConstant / speciesMass[i];
    }
}

RDouble MixedGas::ComputeMixtureGasHeatConductivityByWassilewaFormula(RDouble *primitiveVariables, RDouble temperature)
{
    for (int i = 0; i < numberOfSpecies; ++ i)
    {
        massFractions[i] = primitiveVariables[nm + i];
    }

    //! Obtain the mole fraction of each species.
    this->ComputeMoleFractionByMassFraction(massFractions, moleFractions);

    //! Obtain the molecular weight of each species.
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();

    //! Obtain the viscosity of each species.
    this->ComputeOneTemperatureModelSpeciesViscosity(temperature, speciesViscosity);

    //! Obtain the species specific heat ratio at constant pressure.
    this->ComputeSpeciesConstantPressureSpecificHeat(temperature, speciesCp);

    //! Compute the heat conductivity of each species.
    this->ComputeSpeciesHeatConductivityByEuckenFormula(speciesViscosity, speciesCp, speciesConductivity);

    //! Compute species partition function values.
    this->GetPartitionFunctionValues(speciesViscosity, speciesMass, moleFractions, this->workSpecies);

    RDouble lambda = 0.0;
    for (int i = 0; i < numberOfSpecies; ++ i)
    {
        lambda += speciesConductivity[i] * moleFractions[i] / workSpecies[i];
    }

    return lambda;
}

void MixedGas::GetPartitionFunctionValues(RDouble *target, RDouble *speciesMass, RDouble *moleFraction, RDouble *phi)
{
    //! To obtain the partition function values based on the targets, e.g conductivity or viscosity.
    RDouble tmp = 0.0;
    RDouble tmp0, tmp1, phiTmp;
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        phiTmp = 0.0;
        for (int i = 0; i < numberOfSpecies; ++ i)
        {
            tmp0 = speciesMass[i] / speciesMass[s];
            tmp1 = 1.0 + sqrt(target[s] / target[i]) * sqrt(static_cast<float>(sqrt(static_cast<float>(tmp0))));
            tmp = moleFraction[i] * tmp1 * tmp1;
            phiTmp += tmp / sqrt(static_cast<float>(8.0 + 8.0 / tmp0));
        }
        phi[s] = phiTmp;
    }
}

void MixedGas::ChemicalSource(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, int iCell)
{
}

void MixedGas::EnergySource(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble **energySourceTerms, int nCellNumber)
{
}

RDouble MixedGas::ComputeVibrationAndElectronEnergySourceTerm(RDouble *primitiveVariables, RDouble *omega, RDouble *temperature)
{
    return 0.0;
}

RDouble MixedGas::ComputeVibrationEnergySourceTerm(RDouble *primitiveVariables, RDouble *omega, RDouble *temperature)
{
    return 0.0;
}

void MixedGas::GetTranslationAndVibrationRelaxationTime(Thermo_Param *thermodynamicProperties, RDouble *massFraction,
                                                        RDouble *speciesMass, RDouble temperature, RDouble pressure, RDouble density, RDouble *tvs)
{
}

void MixedGas::GetTranslationAndVibrationRelaxationTimeR(Thermo_Param *thermodynamicProperties, RDouble *massFraction,
                                                        RDouble *speciesMass, RDouble temperature, RDouble pressure, RDouble density, RDouble *tvs)
{
}

void MixedGas::GetElectronAndVibrationRelaxationTime(Thermo_Param *thermodynamicProperties, RDouble electronPressure, RDouble electronTemperature, RDouble *tes)
{
}

RDouble MixedGas::ComputeElectronEnergySourceTerm(RDouble *primitiveVariables, RDouble *omega, RDouble *temperature)
{
    return 0.0;
}

void MixedGas::GetEffectiveCollisionCrossSectionalArea(Thermo_Param *thermodynamicProperties, RDouble electronTemperature, RDouble electronNumberDensity, RDouble *sigma)
{
}

void MixedGas::ChemicalSpectrumRadius(RDouble **q, RDouble **t, RDouble *vol, RDouble **srs, int iCell)
{
    using namespace IDX;
    RDouble *molecularWeightDimensional  = molecularProperty->GetMolecularWeightDimensional();
    RDouble *oMolecularWeightDimensional = molecularProperty->GetMolecularWeightReciprocalDimensional();

    RDouble refDimensionalTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble reference_velocity_dimensional = this->referenceParameterDimensional.GetVelocity();
    RDouble refDimensionalDensity = this->referenceParameterDimensional.GetDensity();

    RDouble temperatureDimensional, densityDimensional, oTemperatureDimensional, oTemperatureDimensional2;
    RDouble diffp, diffn, posr, negr, comm, pztt;
    RDouble xis, srst;
    RDouble density, temperature;
    RDouble odensityDimensionalJS, odensityDimensionalNS;
    RDouble cpDimensional, cvDimensional, omav, coef_tr;

    RDouble lov = reynoldsReferenceLengthDimensional / reference_velocity_dimensional;

    int neqn = nl + nchem;

    RDouble *afr = reactionRate->GetForwardReactionRateCoefficientA();
    RDouble *bfr = reactionRate->GetForwardReactionRateCoefficientB();
    RDouble *cfr = reactionRate->GetForwardReactionRateCoefficientC();

    RDouble *abr = reactionRate->GetBackwardReactionRateCoefficientA();
    RDouble *bbr = reactionRate->GetBackwardReactionRateCoefficientB();
    RDouble *cbr = reactionRate->GetBackwardReactionRateCoefficientC();

    int **cvp = stoichiometricEquation->GetForwardReactionStoichiometricCoefficient();
    int **cvn = stoichiometricEquation->GetBackwardReactionStoichiometricCoefficient();
    RDouble **cpz = stoichiometricEquation->GetReactionStoichiometricThirdBodyCoefficient();

    RDouble *cq = new RDouble [numberOfReaction];
    RDouble *kf = new RDouble [numberOfReaction];
    RDouble *kb = new RDouble [numberOfReaction];

    RDouble *pos = new RDouble[numberOfReaction];
    RDouble *neg = new RDouble[numberOfReaction];
    RDouble *pzt = new RDouble[numberOfReaction];
    int     *npz = new int    [numberOfReaction];

    RDouble *c_a1 = new RDouble[numberOfSpecies];
    RDouble *c_a2 = new RDouble[numberOfSpecies];
    RDouble *c_bk = new RDouble[numberOfSpecies];
    RDouble *c_tr = new RDouble[numberOfSpecies];

    RDouble *prim = new RDouble[neqn];

    density     = q[IR ][iCell];
    temperature = t[ITT][iCell];

    densityDimensional = density * refDimensionalDensity;
    temperatureDimensional = temperature * refDimensionalTemperature;
    oTemperatureDimensional = one / (temperatureDimensional + SMALL);
    oTemperatureDimensional2 = SQR(oTemperatureDimensional);

    for (int m = 0; m < neqn; ++ m)
    {
        prim[m] = q[m][iCell];
    }

    for (int ireaction = 0; ireaction < numberOfReaction; ++ ireaction)
    {
        pos[ireaction] = 1.0;
        neg[ireaction] = 1.0;
        pzt[ireaction] = 0.0;
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            xis = densityDimensional * (prim[nm+ispecies] + SMALL) * oMolecularWeightDimensional[ispecies];

            for (int nexp = 0; nexp < cvp[ireaction][ispecies]; ++ nexp)
            {
                pos[ireaction] *= xis;
            }

            for (int nexp = 0; nexp < cvn[ireaction][ispecies]; ++ nexp)
            {
                neg[ireaction] *= xis;
            }

            pzt[ireaction] += xis * cpz[ireaction][ispecies];
        }

        npz[ireaction] = 0;
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            npz[ireaction] += static_cast<int>(cpz[ireaction][ispecies] + 0.1);
        }

        if (npz[ireaction] == 0)
        {
            pzt[ireaction] = one;
        }

        kf[ireaction] = afr[ireaction] * pow(temperatureDimensional, bfr[ireaction]) * exp(-cfr[ireaction] * oTemperatureDimensional);
        kb[ireaction] = abr[ireaction] * pow(temperatureDimensional, bbr[ireaction]) * exp(-cbr[ireaction] * oTemperatureDimensional);
    }

    for (int ispecies = 0; ispecies < numberOfSpecies - 1; ++ ispecies)
    {
        c_a2[ispecies] = 0.0;
        for (int ireaction = 0; ireaction < numberOfReaction; ++ ireaction)
        {
            diffp = bfr[ireaction] * oTemperatureDimensional + cfr[ireaction] * oTemperatureDimensional2;
            diffn = bbr[ireaction] * oTemperatureDimensional + cbr[ireaction] * oTemperatureDimensional2;
            posr = pos[ireaction] * pzt[ireaction] * kf[ireaction] * diffp;
            negr = neg[ireaction] * pzt[ireaction] * kb[ireaction] * diffn;
            comm = (cvn[ireaction][ispecies] - cvp[ireaction][ispecies]);
            c_a2[ispecies] += comm * (posr - negr);
        }
        c_a2[ispecies] *= molecularWeightDimensional[ispecies];
    }

    ComputeSpeciesEnthalpyDimensional(temperature, speciesEnthalpy);
    ComputeSpeciesConstantPressureSpecificHeatDimensional(temperature, speciesCp);
    ComputeMixtureByPrimitive(prim, speciesCp, cpDimensional);
    ComputeMixtureByPrimitive(prim, oMolecularWeightDimensional, omav);

    cvDimensional = cpDimensional - rjmk * omav;
    coef_tr = - 1.0 / (densityDimensional * cvDimensional);

    for (int ispecies = 0; ispecies < numberOfSpecies - 1; ++ ispecies)
    {
        c_tr[ispecies] = coef_tr * (speciesEnthalpy[ispecies] - speciesEnthalpy[numberOfSpecies-1] - rjmk * temperatureDimensional * 
            (oMolecularWeightDimensional[ispecies] - oMolecularWeightDimensional[numberOfSpecies-1]));
    }

    for (int ispecies = 0; ispecies < numberOfSpecies - 1; ++ ispecies)
    {
        for (int jspecies = 0; jspecies < numberOfSpecies - 1; ++ jspecies)
        {
            c_bk[jspecies] = 0.0;
            odensityDimensionalJS = 1.0 / (densityDimensional * (prim[nm+jspecies         ] + SMALL));
            odensityDimensionalNS = 1.0 / (densityDimensional * (prim[nm+numberOfSpecies-1] + SMALL));
            for (int ireaction = 0; ireaction < numberOfReaction; ++ ireaction)
            {
                diffp = cvp[ireaction][jspecies] * odensityDimensionalJS - cvp[ireaction][numberOfSpecies-1] * odensityDimensionalNS;
                diffn = cvn[ireaction][jspecies] * odensityDimensionalJS - cvn[ireaction][numberOfSpecies-1] * odensityDimensionalNS;
                posr  = pos[ireaction] * kf[ireaction] * diffp;
                negr  = neg[ireaction] * kb[ireaction] * diffn;
                comm  = (cvn[ireaction][ispecies] - cvp[ireaction][ispecies]);
                c_bk[jspecies] += comm * (posr - negr) * pzt[ireaction];
                if (npz[ireaction] != 0)
                {
                    pztt = cpz[ireaction][jspecies] * oMolecularWeightDimensional[jspecies]
                    - cpz[ireaction][numberOfSpecies-1] * oMolecularWeightDimensional[numberOfSpecies-1];
                    c_bk[jspecies] += comm * (kf[ireaction] * pos[ireaction] - kb[ireaction] * neg[ireaction]) * pztt;
                }
            }
            c_bk[jspecies] = molecularWeightDimensional[ispecies] * c_bk[jspecies];
        }
        srst = 0.0;
        for (int jspecies = 0; jspecies < numberOfSpecies - 1; ++ jspecies)
        {
            srst += ABS(c_bk[jspecies] + c_a2[ispecies] * c_tr[jspecies]);
        }
        srs[ispecies][iCell] = srst * lov * vol[iCell];
    }

    delete [] cq;
    delete [] kf;
    delete [] kb;

    delete [] pos;
    delete [] neg;
    delete [] pzt;
    delete [] npz;

    delete [] c_a1;
    delete [] c_a2;
    delete [] c_bk;
    delete [] c_tr;

    delete [] prim;
}

void MixedGas::ComputeChemicalSourceDerivatives(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber)
{
}

void MixedGas::ComputeChemicalSourceJacobian(RDouble *primitiveVariables, RDouble *temperatures, RDouble cellVolume, RDouble *chemicalSourceTerms, RDouble **sourceDerivatives)
{
}

void MixedGas::ComputeChemicalSourceAndJacobian(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber)
{
}

void MixedGas::ComputeChemicalSourceAndJacobianD3(int *Inzone, RDouble4D &primitiveVariables, RDouble4D &temperatures, RDouble4D &speciesEs, RDouble4D &speciesCvs, RDouble3D &totalCv, RDouble3D &cellVolume, RDouble4D &chemicalSourceTerms, RDouble5D &sourceDerivatives)
{
}

void MixedGas::ComputeChemicalSourceDerivatives1(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber)
{
}

void MixedGas::ComputeChemicalSourceDerivatives2(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber)
{
}

void MixedGas::ComputeEnergySourceDerivatives1(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
                                               RDouble ***chemicalSourceDerivatives, RDouble ***energySourceDerivatives, int nCellNumber)
{

}

void MixedGas::ComputeEnergySourceDiagonalDerivatives1(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
                                                       RDouble ***chemicalSourceDerivatives, RDouble **energySourceDerivatives, int nCellNumber)
{
 
}

void MixedGas::ComputeEnergySourceDiagonalDerivatives2(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
                                                       RDouble ***chemicalSourceDerivatives, RDouble **energySourceDerivatives, int nCellNumber)
{

}

void MixedGas::ComputeEnergySourceJacobianT2(RDouble *primitiveVariables, RDouble *temperatures, RDouble cellVolume, RDouble *chemicalSourceTerms, RDouble **sourceDerivatives)
{
 
}

void MixedGas::ComputeEnergySourceJacobianT3(RDouble *primitiveVariables, RDouble *temperatures, RDouble cellVolume, RDouble *chemicalSourceTerms, RDouble **sourceDerivatives)
{

}

void MixedGas::ComputeEnergySourceAndJacobianT2(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber)
{

}

void MixedGas::ComputeEnergySourceAndJacobianT2D3(int *Inzone, RDouble4D &primitiveVariables, RDouble4D &temperatures, RDouble4D &speciesCvvs, RDouble4D &speciesCves, RDouble4D &speciesEvs, RDouble4D &speciesEes, RDouble3D &totalCvtr, RDouble3D &totalCvv, RDouble3D &totalCve, RDouble3D &cellVolume, RDouble4D &chemicalSourceTerms, RDouble5D &sourceDerivatives)
{                                              

}

void MixedGas::ComputeEnergySourceAndJacobianT3(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber)
{

}

void MixedGas::ComputeEnergySourceAndJacobianT3D3(int *Inzone, RDouble4D &primitiveVariables, RDouble4D &temperatures, RDouble4D &speciesCvvs, RDouble4D &speciesCves, RDouble4D &speciesEvs, RDouble4D &speciesEes, RDouble3D &totalCvtr, RDouble3D &totalCvv, RDouble3D &totalCve, RDouble3D &cellVolume, RDouble4D &chemicalSourceTerms, RDouble5D &sourceDerivatives)
{
 
}

void MixedGas::GetLinearlyAveragedPolynomialCoefficients(int nSpecies, RDouble temperature, RDouble *polynomialCoef)
{
    int nIndex1, nIndex2;      //! Record the index of temperature range.
    RDouble coef1, coef2;      //! The coefficients for linear polynomial on the temperature boundary.
    RDouble T = temperature;   //! temperature is the dimensional value.

    //! Select the temperature range and compute the coefficients.
    if (T > 15500.0)
    {
        nIndex1 = 3;
        nIndex2 = 4;
        if (T > 25500.0)     //! Range: (25500, 30000).
        {
            coef1 = 0.0;
            coef2 = 1.0;
        }
        else if (T <= 24500.0)     //! Range: (15500, 24500).
        {
            coef1 = 1.0;
            coef2 = 0.0;
        }
        else     //! Range: (24500, 25500).
        {
            coef2 = 0.001 * (T - 24500.0);
            coef1 = 1.0 - coef2;
        }
    }
    else if (T > 6500.0)
    {
        nIndex1 = 2;
        nIndex2 = 3;
        if (T <= 14500.0)     //! Range: (6500, 14500).
        {
            coef1 = 1.0;
            coef2 = 0.0;
        }
        else      //! Range: (14500, 15500).
        {
            coef2 = 0.001 * (T - 14500.0);
            coef1 = 1.0 - coef2;
        }
    }
    else if (T > 1200.0)
    {
        nIndex1 = 1;
        nIndex2 = 2;
        if (T <= 5500.0)     //! Range: (1200, 5500).
        {
            coef1 = 1.0;
            coef2 = 0.0;
        }
        else      //! Range: (5500, 6500).
        {
            coef2 = 0.001 * (T - 5500.0);
            coef1 = 1.0 - coef2;
        }
    }
    else
    {
        nIndex1 = 0;
        nIndex2 = 1;
        if (T <= 800.0)      //! Range: (0, 800).
        {
            coef1 = 1.0;
            coef2 = 0.0;
        }
        else      //! Range: (800, 1200).
        {
            coef2 = (T - 800.0)/400.0;
            coef1 = 1.0 - coef2;
        }
    }

    RDouble *a = this->thermodynamicManager->GetPolynomialCoefficient(nSpecies, nIndex1);
    RDouble *b = this->thermodynamicManager->GetPolynomialCoefficient(nSpecies, nIndex2);
    int n = thermodynamicManager->GetNumberofPolynomialCoefficient();
    for (int i = 0; i < n; ++ i)
    {
        polynomialCoef[i] = coef1 * a[i] + coef2 * b[i];
    }
}

void MixedGas::GetLinearlyAveragedPolynomialCoefficientsT2(RDouble temperature, RDouble **polynomialCoef)
{
    int nIndex1, nIndex2;      //! Record the index of temperature range.
    RDouble coef1, coef2;      //! The coefficients for linear polynomial on the temperature boundary.
    RDouble T = temperature;   //! temperature is the dimensional value.

    //! Select the temperature range and compute the coefficients.
    if (T > 15500.0)
    {
        nIndex1 = 4;
        nIndex2 = 5;
        if (T > 25500.0)     //! Range: (25500, 30000).
        {
            coef1 = 0.0;
            coef2 = 1.0;
        }
        else if (T <= 24500.0)     //! Range: (15500, 24500).
        {
            coef1 = 1.0;
            coef2 = 0.0;
        }
        else     //! Range: (24500, 25500).
        {
            coef2 = 0.001 * (T - 24500.0);
            coef1 = 1.0 - coef2;
        }
    }
    else if (T > 6500.0)
    {
        nIndex1 = 3;
        nIndex2 = 4;
        if (T <= 14500.0)     //! Range: (6500, 14500).
        {
            coef1 = 1.0;
            coef2 = 0.0;
        }
        else      //! Range: (14500, 15500).
        {
            coef2 = 0.001 * (T - 14500.0);
            coef1 = 1.0 - coef2;
        }
    }
    else if (T > 1200.0)
    {
        nIndex1 = 2;
        nIndex2 = 3;
        if (T <= 5500.0)     //! Range: (1200, 5500).
        {
            coef1 = 1.0;
            coef2 = 0.0;
        }
        else      //! Range: (5500, 6500).
        {
            coef2 = 0.001 * (T - 5500.0);
            coef1 = 1.0 - coef2;
        }
    }
    else
    {
        nIndex1 = 1;
        nIndex2 = 2;
        if (T <= 800.0)      //! Range: (0, 800).
        {
            coef1 = 1.0;
            coef2 = 0.0;
        }
        else      //! Range: (800, 1200).
        {
            coef2 = (T - 800.0) / 400.0;
            coef1 = 1.0 - coef2;
        }
    }

    RDouble benchmarkTemp = 200, T1 = 0.0, T2 = 0.0;
    RDouble *polynomial1 = NULL, *polynomial2 = NULL;

    int n = thermodynamicManager->GetNumberofPolynomialCoefficient();
    if (ABS(coef1 - 1) <= EPSILON)
    {
        for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
        {
            polynomial1 = this->thermodynamicManager->GetPolynomialCoefficientT2(iSpecies, nIndex1);
            for (int i = 0; i < n; ++ i)
            {
                polynomialCoef[iSpecies][i] = polynomial1[i];
            }
        }
    }
    else if (ABS(coef2 - 1) <= EPSILON)
    {
        for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
        {
            polynomial2 = this->thermodynamicManager->GetPolynomialCoefficientT2(iSpecies, nIndex2);
            for (int i = 0; i < n; ++ i)
            {
                polynomialCoef[iSpecies][i] = polynomial2[i];
            }
        }
    }
    else
    {
        for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
        {
            polynomial1 = this->thermodynamicManager->GetPolynomialCoefficientT2(iSpecies, nIndex1);
            polynomial2 = this->thermodynamicManager->GetPolynomialCoefficientT2(iSpecies, nIndex2);
            for (int i = 0; i < n; ++ i)
            {
                polynomialCoef[iSpecies][i] = coef1 * polynomial1[i] + coef2 * polynomial2[i];
            }
        }
    }

    //Modify the coefficient when the temperature is low.
    bool isSmooth = false;
    for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
    {
        benchmarkTemp = this->thermodynamicManager->GetBenchmarkTemperature(iSpecies);
        if (isSmooth)
        {
            T1 = benchmarkTemp - 100;
            T2 = benchmarkTemp + 100;
            if (T < T1)
            {
                polynomial1 = this->thermodynamicManager->GetPolynomialCoefficientT2(iSpecies, 0);
                for (int i = 0; i < n; ++ i)
                {
                    polynomialCoef[iSpecies][i] = polynomial1[i];
                }
            }
            else if (T <= T2)
            {
                polynomial1 = this->thermodynamicManager->GetPolynomialCoefficientT2(iSpecies, 0);
                polynomial2 = this->thermodynamicManager->GetPolynomialCoefficientT2(iSpecies, 1);
                coef2 = (T - T1) / 200.0;
                coef1 = 1.0 - coef2;
                for (int i = 0; i < n; ++ i)
                {
                    polynomialCoef[iSpecies][i] = coef1 * polynomial1[i] + coef2 * polynomial2[i];
                }
            }
        }
        else //No smooth.
        {
            if (T < benchmarkTemp)
            {
                polynomial1 = this->thermodynamicManager->GetPolynomialCoefficientT2(iSpecies, 0);
                for (int i = 0; i < n; ++ i)
                {
                    polynomialCoef[iSpecies][i] = polynomial1[i];
                }
            }
        }
    }
}

void MixedGas::GetLinearlyAveragedPolynomialCoefficientsIndex(RDouble temperature, RDouble &coef, int &nIndex1, int &nIndex2)
{
}

void MixedGas::GetLinearlyAveragedPolynomialCoefficientsByIndex(int nSpecies, RDouble *polynomialCoef, RDouble coef, int nIndex1, int nIndex2)
{
}

void MixedGas::GetAirInfo(const RDouble &h1, RDouble &t, RDouble &p, RDouble &den, RDouble &a)
{
    //! USSA 76 Model for atmosphere.
    RDouble h      = h1 * 1000.0;
    RDouble r      = 287.053;
    RDouble g0     = 9.80665;
    RDouble rp     = 6.37111e6;
    RDouble g      = (rp/(rp+h))*(rp/(rp+h))*g0;
    RDouble t0     = 288.15;
    RDouble p0     = 10.1325e2;
    RDouble rho0   = 1.225;
    RDouble t11    = 216.65;
    RDouble p11    = 2.2632e2;
    RDouble rho11  = 3.6392e-1;
    RDouble t20    = t11;
    RDouble p20    = 5.4747e1;
    RDouble rho20  = 8.8035e-2;
    RDouble t32    = 228.65;
    RDouble p32    = 8.6789;
    RDouble rho32  = 1.3225e-2;
    RDouble t47    = 270.65;
    RDouble p47    = 1.1090;
    RDouble rho47  = 1.4275e-3;
    RDouble t52    = t47;
    RDouble p52    = 5.8997e-1;
    RDouble rho52  = 7.5943e-4;
    RDouble t61    = 252.65;
    RDouble p61    = 1.8209e-1;
    RDouble rho61  = 2.5109e-4;
    RDouble t79    = 180.65;
    RDouble p79    = 1.0376e-2;
    RDouble rho79  = 2.0010e-5;
    RDouble t90    = t79;
    RDouble p90    = 1.6437e-3;
    RDouble rho90  = 3.4165e-6;
    RDouble t100   = 210.02;
    RDouble p100   = 3.0070e-4;
    RDouble rho100 = 5.6044e-7;
    RDouble t110   = 257.00;
    RDouble p110   = 7.3527e-5;
    RDouble rho110 = 9.7081e-8;
    RDouble t120   = 349.49;
    RDouble p120   = 2.5209e-5;
    RDouble rho120 = 2.2222e-8;
    RDouble t150   = 892.79;
    RDouble p150   = 5.0599e-6;
    RDouble rho150 = 2.0752e-9;
    RDouble t160   = 1022.20;
    RDouble p160   = 3.6929e-6;
    RDouble rho160 = 1.2336e-9;
    RDouble t170   = 1103.40;
    RDouble p170   = 2.7915e-6;
    RDouble rho170 = 7.8155e-10;
    RDouble t190   = 1205.40;
    RDouble p190   = 1.6845e-6;
    RDouble rho190 = 3.5807e-10;
    RDouble t230   = 1322.30;
    RDouble p230   = 6.7138e-7;
    RDouble rho230 = 1.5640e-10;
    RDouble t300   = 1432.10;
    RDouble p300   = 1.8828e-7;
    RDouble rho300 = 1.9159e-11;
    RDouble t400   = 1487.40;
    RDouble p400   = 4.0278e-8;
    RDouble rho400 = 2.8028e-12;
    RDouble t500   = 1499.20;
    RDouble p500   = 1.0949e-8;
    RDouble rho500 = 5.2148e-13;
    RDouble t600   = 1506.10;
    RDouble p600   = 3.4475e-9;
    RDouble rho600 = 1.1367e-13;
    RDouble t700   = 1507.60;

    RDouble rho = 0.0;

    if (h <= 11019.0) 
    {
        RDouble al1 = (t11 - t0) / 11019.0;
        t   = t0 + al1 * h;
        p   = p0 * pow(t/t0, -g / (r * al1));
        rho = rho0 * pow(t/t0, -1.0 - g / (r*al1));
    }
    else if (h <= 20063.0)
    {
        t   = t11;
        p   = p11   * exp(-g*(h-11019.0)/(r*t11));
        rho = rho11 * exp(-g*(h-11019.0)/(r*t11));
    }
    else if (h <= 32162.0)
    {
        RDouble al2 = (t32-t20) / (32162.0-20063.0);
        t   = t11 + al2*(h-20063.0);
        p   = p20 * pow(t/t11, -g/(r*al2));
        rho = rho20 * pow(t/t11, -1.0-g/(r*al2));
    }
    else if (h <= 47350.0)
    {
        RDouble al3 = (t47-t32) / (47350.0-32162.0);
        t   = t32 + al3 * (h-32162.0);
        p   = p32 * pow(t/t32, -g/(r*al3));
        rho = rho32 * pow(t/t32, -1.0-g/(r*al3));
    }

    else if (h <= 52429.0)
    {
        t   = t47;
        p   = p47*exp(-g*(h-47350.0)/(r*t47));
        rho = rho47*exp(-g*(h-47350.0)/(r*t47));
    }
    else if (h <= 61591.0)
    {
        RDouble al4 =(t61-t52)/(61591.0-52429.0);
        t = t47 + al4*(h-52429.0);
        p = p52 * pow(t/t47, -g/(r*al4));
        rho = rho52* pow(t/t47, -1.0-g/(r*al4));
    }
    else if (h <= 79994.0)
    {
        RDouble al5 = (t79-t61) / (79994.0-61591.0);
        t = t61 + al5 * (h-61591.0);
        p = p61 * pow(t/t61, -g/(r*al5));
        rho = rho61 * pow(t/t61, -1.0-g/(r*al5));
    }
    else if (h <= 90000.0)
    {
        t = t79;
        p = p79 * exp(-g*(h-79994.0)/(r*t79));
        rho = rho79 * exp(-g*(h-79994.0)/(r*t79));
    }
    else if (h <= 100000.0)
    {
        RDouble al6 = (t100-t90)/10000.0;
        t = t79 + al6*(h-90000.0);
        p = p90 * pow(t/t79, -g/(r*al6));
        rho = rho90 * pow(t/t79, -1.0-g/(r*al6));
    }
    else if (h <= 110000.0)
    {
        RDouble al7 = (t110-t100)/10000.0;
        t = t100 + al7*(h-100000.0);
        p = p100 * pow(t/t100, -g/(r*al7));
        rho = rho100 * (t/t100, -1.0-g/(r*al7));
    }
    else if (h <= 120000.0)
    {
        RDouble al8 = (t120-t110)/10000.0;
        t = t110 + al8*(h-110000.0);
        p = p110 * pow(t/t110, -g/(r*al8));
        rho = rho110 * pow(t/t110, -1.0-g/(r*al8));
    }
    else if (h <= 150000.0)
    {
        RDouble al9 = (t150-t120)/30000.0;
        t = t120 + al9*(h-120000.0);
        p = p120 * pow(t/t120, -g/(r*al9));
        rho = rho120 * pow(t/t120, -1.0-g/(r*al9));
    }
    else if (h <= 160000.0)
    {
        RDouble al10 = (t160-t150)/10000.0;
        t = t150 + al10*(h-150000.0);
        p = p150 * pow(t/t150, -g/(r*al10));
        rho = rho150 * pow(t/t150, -1.0-g/(r*al10));
    }
    else if (h <= 170000.0)
    {
        RDouble al11 = (t170-t160)/10000.0;
        t = t160 + al11*(h-160000.0);
        p = p160 * pow(t/t160, -g/(r*al11));
        rho = rho160 * pow(t/t160, -1.0-g/(r*al11));
    }
    else if (h <= 190000.0)
    {
        RDouble al12 = (t190-t170)/20000.0;
        t = t170 + al12*(h-170000.0);
        p = p170 * pow(t/t170, -g/(r*al12));
        rho = rho170 * pow(t/t170, -1.0-g/(r*al12));
    }
    else if (h <= 230000.0)
    {
        RDouble al13 = (t230-t190)/40000.0;
        t = t190 + al13*(h-190000.0);
        p = p190 * pow(t/t190, -g/(r*al13));
        rho = rho190 * pow(t/t190, -1.0-g/(r*al13));
    }
    else if (h <= 300000.0)
    {
        RDouble al14 = (t300-t230)/70000.0;
        t = t230 + al14*(h-230000.0);
        p = p230 * pow(t/t230, -g/(r*al14));
        rho = rho230 * pow(t/t230, -1.0-g/(r*al14));
    }
    else if (h <= 400000.0)
    {
        RDouble al15 = (t400-t300)/100000.0;
        t = t300 + al15*(h-300000.0);
        p = p300 * pow(t/t300, -g/(r*al15));
        rho = rho300 * pow(t/t300, -1.0-g/(r*al15));
    }
    else if (h <= 500000.0)
    {
        RDouble al16 = (t500-t400)/100000.0;
        t = t400 + al16*(h-400000.0);
        p = p400 * pow(t/t400, -g/(r*al16));
        rho = rho400 * pow(t/t400, -1.0-g/(r*al16));
    }
    else if (h <= 600000.0)
    {
        RDouble al17 = (t600-t500)/100000.0;
        t = t500 + al17*(h-500000.0);
        p = p500 * pow(t/t500, -g/(r*al17));
        rho = rho500 * pow(t/t500, -1.0-g/(r*al17));
    }
    else if (h <= 700000.0)
    {
        RDouble al18 = (t700-t600)/100000.0;
        t = t600 + al18*(h-600000.0);
        p = p600 * pow(t/t600, -g/(r*al18));
        rho = rho600 * pow(t/t600, -1.0-g/(r*al18));
    }
    else
    {
        TK_Exit::ExceptionExit("Error: height must be less than 700km!!!\n", true);
    }

    a   = sqrt(1.4*r*t);
    p   = p * 100.0;
    den = rho;
}

void MixedGas::ComputeCoefficientInMXDQ(RDouble *primitiveVars, RDouble &alpha, RDouble &beta, RDouble *speciesBeta)
{
    //! The arguments alpha, beta and beta_s are referred to the formular (A.5) in the appendix A of the PHengLEI Theory manual.
    //! alpha = R/(M*cv)
    //! beta = R/(M*cv)*(0.5*V^2 - ens) + R*T/Mns
    //! betas =  R/(M*cv) * (ens - es) + R*T(1/Ms - 1/Mns)
    int nSpecies = this->GetNumberOfSpecies();
    for (int s = 0; s < nSpecies; ++ s)
    {
        massFractions[s] = primitiveVars[this->nm + s];
    }

    //! Obtain the primitive values.
    using namespace IDX;
    RDouble &density   = primitiveVars[IR];
    RDouble &uVelocity = primitiveVars[IU];
    RDouble &vVelocity = primitiveVars[IV];
    RDouble &wVelocity = primitiveVars[IW];
    RDouble &pressure  = primitiveVars[IP];
    RDouble squareVelocity = uVelocity * uVelocity + vVelocity * vVelocity + wVelocity * wVelocity;
    //! Obtain the non-dimensional universal gas constant.
    RDouble gasConstant = this->GetUniversalGasConstant();
    //! Obtain the molecular weight of mixture gas.
    RDouble mixedGasMass = this->GetMixedGasMolecularWeight(massFractions);
    //! Compute the temperature using the state equation of mixture gas.
    RDouble temperature = pressure * mixedGasMass / (density * gasConstant);

    //! Compute the internal energy of last species with the label ns-1.
    this->ComputeSpeciesConstantVolumeSpecificHeat(temperature, speciesCv);
    //this->ComputeSpeciesEnthalpy(temperature, speciesCv, speciesEnergyTr, speciesEnthalpy);
    RDouble mixedGasCv = 0.0;
    this->ComputeMixtureByMassFraction(massFractions, speciesCv, mixedGasCv);
    //! Compute the first coefficient called alpha.
    alpha = gasConstant / (mixedGasMass * mixedGasCv);

    //! Obtain the molecular weight of each species.
    RDouble *speciesMassReciprocal = this->molecularProperty->GetMolecularWeightReciprocal();
    RDouble nitrogenEnergy = speciesCv[nLastSpeciesIndex] * temperature;//speciesEnergyTr[nLastSpeciesIndex];
    RDouble electronEnergy = 0.0, electronMassReciprocal = 1.0 / electronMass;
    if (nElectronIndex >= 0)
    {
        electronEnergy = speciesCv[nElectronIndex] * temperature;
        //electronEnergy = speciesEnergyTr[nElectronIndex];
        electronMassReciprocal = speciesMassReciprocal[nElectronIndex];
    }

    //! Compute the second coefficient called beta.
    RDouble alpha1 = 0.5 * squareVelocity - nitrogenEnergy;
    RDouble alpha2 = speciesMassReciprocal[nLastSpeciesIndex];
    beta = alpha * alpha1 + gasConstant * temperature * alpha2;

    //! Compute the other coefficient stored in speciesBeta.
    RDouble speciesGama, speciesEta, tempCoef;
    int *ionType = molecularProperty->GetIonTypeOfSpecies();
    speciesBeta[nLastSpeciesIndex] = 0.0;
    for (int s = 0; s < nLastSpeciesIndex; ++ s)
    {
        tempCoef = ionType[s] * speciesMassReciprocal[s] / electronMassReciprocal;

        speciesGama = speciesMassReciprocal[s] - speciesMassReciprocal[nLastSpeciesIndex];
        speciesGama += tempCoef * (electronMassReciprocal - speciesMassReciprocal[nLastSpeciesIndex]);

        speciesEta = nitrogenEnergy - speciesCv[s] * temperature;
        //speciesEta = nitrogenEnergy - speciesEnergyTr[s];
        speciesEta += tempCoef * (nitrogenEnergy - electronEnergy);

        speciesBeta[s] = alpha * speciesEta + gasConstant * temperature * speciesGama;
    }
}

void MixedGas::ComputeCoefficientInMXDQR(RDouble *primitiveVars, RDouble trTemperature, RDouble squareVelocity, RDouble *speciesCvs, RDouble totalCv, RDouble &alpha, RDouble &beta, RDouble *speciesBeta)
{
    //! The arguments alpha, beta and beta_s are referred to the formular (A.5) in the appendix A of the PHengLEI Theory manual.
    //! alpha = R/(M*cv)
    //! beta = R/(M*cv)*(0.5*V^2 - ens) + R*T/Mns
    //! betas =  R/(M*cv) * (ens - es) + R*T(1/Ms - 1/Mns)
    //! Obtain the primitive values.

    RDouble gasConstant = this->GetUniversalGasConstant();
    RDouble *speciesMassReciprocal = this->molecularProperty->GetMolecularWeightReciprocal();
    int *ionType = molecularProperty->GetIonTypeOfSpecies();

    RDouble mixedGasMass = 0.0;
    for (int i = 0; i < numberOfSpecies; ++ i)
    {
        massFractions[i] = primitiveVars[this->nm + i];
        mixedGasMass += massFractions[i] * speciesMassReciprocal[i];
    }

    alpha = gasConstant * mixedGasMass / totalCv;

    //! Obtain the molecular weight of each species.

    RDouble nitrogenEnergy = speciesCvs[nLastSpeciesIndex] * trTemperature;   //speciesEnergyTr[nLastSpeciesIndex];

    RDouble electronEnergy = 0.0, electronMassReciprocal = 1.0 / electronMass;
    if (nElectronIndex >= 0)
    {
        electronEnergy = speciesCvs[nElectronIndex] * trTemperature;
        //electronEnergy = speciesEnergyTr[nElectronIndex];
        electronMassReciprocal = speciesMassReciprocal[nElectronIndex];
    }

    //! Compute the second coefficient called beta.
    RDouble alpha1 = 0.5 * squareVelocity - nitrogenEnergy;
    RDouble alpha2 = speciesMassReciprocal[nLastSpeciesIndex];
    beta = alpha * alpha1 + gasConstant * trTemperature * alpha2;

    //! Compute the other coefficient stored in speciesBeta.
    RDouble speciesGama, speciesEta, tempCoef;

    for (int s = 0; s < nLastSpeciesIndex; ++ s)
    {
        tempCoef = ionType[s] * speciesMassReciprocal[s] / electronMassReciprocal;

        speciesGama = speciesMassReciprocal[s] - speciesMassReciprocal[nLastSpeciesIndex];
        speciesGama += tempCoef * (electronMassReciprocal - speciesMassReciprocal[nLastSpeciesIndex]);

        speciesEta = nitrogenEnergy - speciesCvs[s] * trTemperature;
        //speciesEta = nitrogenEnergy - speciesEnergyTr[s];
        speciesEta += tempCoef * (nitrogenEnergy - electronEnergy);

        speciesBeta[s] = alpha * speciesEta + gasConstant * trTemperature * speciesGama;
    }
    speciesBeta[nLastSpeciesIndex] = 0.0;
}

//! To obtain the coefficients of computing the M*dQ in the function MVXDQ_STD_LP(), added by LiPeng on Jan. 21, 2019.
void MixedGas::ComputeCoefficientInMVXDQ(RDouble *prim, RDouble K, RDouble rD, RDouble *theta_s, RDouble *phi_s, RDouble *Ds)
{
    //! K=visl/Prl + vist/Prt denotes the thermal conductivity coefficient, rD=visl/Scl + vist/Sct is the species diffusion coefficient,\n
    //! Ds is an array of saving the diffusion coefficient of each species.

    int ns = this->GetNumberOfSpecies();
    //! To obtain the mass fraction of each species.
    RDouble *fs = new RDouble[ns];
    RDouble *xs = new RDouble[ns];
    for (int i = 0; i < ns; i++)
    {
        fs[i] = prim[i + this->nm];
    }

    //! Obtain the mole fraction of each species.
    this->ComputeMoleFractionByMassFraction(fs, xs);

    //! Obtain the density and pressure.
    RDouble &ro = prim[0];
    RDouble &pm = prim[4];
    RDouble R = this->GetUniversalGasConstant();    //! The non-dimensional universal gas constant.
    RDouble M = this->GetMixedGasMolecularWeight(fs);     //! The molecular weight of mixture gas.
    RDouble tm = pm * M / (ro * R);     //! Compute the temperature using the state equation of mixture gas.

    RDouble cp;    //! Save the value of constant pressure specific heat.
    this->ComputeConstantPressureSpecificHeatByMassFraction(fs, tm, cp);
    RDouble lambda = cp * K;      //! K = visl/Prl + vist/Prt denotes the thermal conductivity coefficient.

    //! Compute the diffusion coefficient of each species.
    for (int i = 0; i < ns; i++)
    {
        Ds[i] = (rD/ro) * (1.0 - fs[i]) / (1.0 - xs[i]);     //! rD=visl/Scl + vist/Sct.
    }

    RDouble cv;//!Save the value of constant volume specific heat.
    this->ComputeConstantVolumeSpecificHeatByMassFraction(fs, tm, cv);
    RDouble e = cv * tm;//!the total internal energy of mixture gas

    //! Compute the internal energy of last species with label ns-1 in the collection.
    RDouble *cvs = new RDouble[ns];    //! Save the value of constant volume specific heat of each species.
    this->ComputeSpeciesConstantVolumeSpecificHeat(tm, cvs);
    RDouble ens = cvs[ns-1] * tm;
    RDouble *Ms = this->molecularProperty->GetMolecularWeight();    //! Obtain the molecular weight of each species.

    //! Compute the coefficient called theta_s.
    RDouble ei = 0.0;
    for (int i = 0; i < ns; i++)
    {
        ei = cvs[i] * tm;
        theta_s[i] = (ei - ens) - M * (1.0/Ms[i] - 1.0/Ms[ns-1]) * e;
    }

    //! Compute the coefficient called phi_s.
    RDouble *cps = new RDouble[ns];    //! Save the value of constant pressure specific heat of each species.
    this->ComputeSpeciesConstantPressureSpecificHeat(tm, cps);
    RDouble roDns = ro * Ds[ns-1] * cps[ns-1] * tm;
    RDouble hi = 0.0;

    //! Compute the coefficient called phi_s of each species.
    for (int i = 0; i < ns; i ++)
    {//!lambda=cp*(visl/Prl + vist/Prt)
        hi = cps[i] * tm;
        phi_s[i] = -lambda * M * (1.0/Ms[i] - 1.0/Ms[ns-1]) * tm + (ro * Ds[i] * hi  - roDns);
    }

    //! Remove the dynamic memory space.
    delete [] fs;
    delete [] xs;
    delete [] cvs;
    delete [] cps;
}

void MixedGas::GetEarthFullyCatalyticMassFraction(RDouble *fs, RDouble *fcw)
{
    RDouble *Ms = this->molecularProperty->GetMolecularWeightDimensional();
    RDouble M = 0.0;
    //! Compute the total molecular weight.
    for (int i = 0; i < this->numberOfSpecies; i ++)
    {
        fcw[i] = 0.0; //Initialization.
        M += fs[i] / Ms[i];
    }
    M = 1.0 / M;

    RDouble *xs = new RDouble[this->numberOfSpecies];
    //! Compute mole fraction of each species.
    for (int i = 0; i < this->numberOfSpecies; i ++)
    {
        xs[i] = fs[i] * M / Ms[i];
    }

    //! To compute the mass of oxygen.
    RDouble co2 = 0.0, cn2 = 0.0;
    int coef = 0;
    RDouble Mo = 15.9994E-3;     //! The molecular weight.
    string *species_name = this->molecularProperty->GetNameOfSpecies();
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        if (FindChemicalElement(species_name[i], 'O', coef))
        {
            co2 += xs[i] * coef * Mo;
        }
    }
    delete [] xs;
    //Compute the mass fraction of Oxygen.
    co2 = co2 / M;
    //Compute the mass fraction of Nitrogen.
    cn2 = MAX(1.0 - co2, 0.0);

    int nIndex = this->GetSpeciesIndex("O2");
    if (nIndex >= 0)
    {
        fcw[nIndex] = co2;
    }
    nIndex = this->GetSpeciesIndex("N2");
    if (nIndex >= 0)
    {
        fcw[nIndex] = cn2;
    }
}

void MixedGas::GetMarsFullyCatalyticMassFraction(RDouble *fs, RDouble *fcw)
{
    RDouble *Ms = this->molecularProperty->GetMolecularWeightDimensional();
    RDouble M = 0.0;
    //! Compute the total molecular weight.
    for (int i = 0; i < this->numberOfSpecies; i++)
    {
        fcw[i] = 0.0; //Initialization.
        M += fs[i] / Ms[i];
    }
    M = 1.0 / M;

    RDouble *xs = new RDouble[this->numberOfSpecies];
    //! Compute mole fraction of each species.
    for (int i = 0; i < this->numberOfSpecies; i ++)
    {
        xs[i] = fs[i] * M / Ms[i];
    }

    //! To compute the mass of Oxygen, Nitrogen and Carbon dioxide.
    RDouble cO2 = 0.0, cN2 = 0.0, cCO2 = 0.0;
    int coef = 0;
    RDouble atomicNitrogenMass = 1.400674E-2, carbonDioxideMass = 4.401E-2; //! The molecular weight(kg/mol).
    string *species_name = this->molecularProperty->GetNameOfSpecies();

    RDouble moleFracationCO2 = 0.0;
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        //Compute the total mass of Nitrogen.
        if (FindChemicalElement(species_name[i], 'N', coef))
        {
            cN2 += xs[i] * coef * atomicNitrogenMass;
        }
        //Compute the total mass of Carbon.
        if (FindChemicalElement(species_name[i], 'C', coef))
        {
            moleFracationCO2 += xs[i] * coef; //Obtain the total mole fractions of CO2;
        }
    }
    delete [] xs;

    //Compute the mass fraction of Nitrogen.
    cN2 = cN2 / M;
    //Compute the mass fraction of Carbon dioxide.
    cCO2 = moleFracationCO2 * carbonDioxideMass / M;
    //Compute the mass fraction of Oxygen.
    cO2 = MAX(1.0 - cN2 - cCO2, 0.0);

    if (cCO2 > 1.0 && cN2 == 0 && cO2 == 0)
    {
        cCO2 = 1.0;
    }

    int nIndex = this->GetSpeciesIndex("O2");
    if (nIndex >= 0)
    {
        fcw[nIndex] = cO2;
    }
    nIndex = this->GetSpeciesIndex("N2");
    if (nIndex >= 0)
    {
        fcw[nIndex] = cN2;
    }
    nIndex = this->GetSpeciesIndex("CO2");
    if (nIndex >= 0)
    {
        fcw[nIndex] = cCO2;
    }
}

bool MixedGas::FindChemicalElement(string species_name, char element_name, int &coef)
{
    bool bflag = false;
    coef = 0;
    int i = static_cast<int> (species_name.find(element_name));
    if (i >= 0)
    {
        coef = 1;
        bflag = true;
        int n = static_cast<int> (species_name.length());
        if (i + 1 < n)
        {
            if (species_name[i + 1] >= '1' && species_name[i + 1] <= '9')
            {
                coef = atoi(species_name.substr(i + 1, 1).c_str());
            }

            if (element_name == 'N' && species_name[i + 1] == 'a')  //Na
            {
                coef = 0;
                bflag = false;   //这个地方没考虑Na的氮化物
            }

            if (element_name == 'S' && species_name[i + 1] == 'i')  //Si
            {
                coef = 0;     //
                bflag = false;   //这个地方没考虑Si的硫化物
            }
        }
    }

    return bflag;
}

int MixedGas::GetSpeciesIndex(string species_name)
{
    string *name_array = this->molecularProperty->GetNameOfSpecies();
    for (int i = 0; i < this->numberOfSpecies; i ++)
    {
        if (name_array[i] == species_name)
        {
            return i;
        }
    }
    return -1;
}

void MixedGas::ComputeTranslationAndRotationSpecificHeatAtConstantVolume(RDouble *speciesTransRotationCv)
{
    //! The dimensional value of translation-rotation specific heat at constant volume.
    ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesTransRotationCv);

    //! Translate to the non-dimensional values.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        speciesTransRotationCv[i] *= referenceSpecificHeat;
    }
}

void MixedGas::ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(RDouble *speciesTransRotationCv)
{
    RDouble gasConstant = rjmk; //Unit: J/(kg*K)
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = molecularProperty->GetThermodynamicProperties();
    RDouble coef = 1.0;
    for (int i = 0; i <= this->nLastSpeciesIndex; ++ i)
    {
        speciesTransRotationCv[i] = 1.5 * gasConstant / speciesMass[i];
        coef = (RDouble)(thermodynamicProperties[i].nType - 1);
        //! Polyatom species.
        speciesTransRotationCv[i] += coef * gasConstant / speciesMass[i];
    }
    if (this->nElectronIndex >= 0)
    {
        speciesTransRotationCv[this->nElectronIndex] = 0.0;
    }
}

void MixedGas::ComputeTranslationAndRotationEnergy(RDouble transRotationTemperature, RDouble *speciesTransRotationEnergy)
{
    //! Obtain the translation-rotation specific heat at constant volume.
    ComputeTranslationAndRotationSpecificHeatAtConstantVolume(speciesTransRotationEnergy);

    //! Obtain the thermodynamic properties of species.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();

    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        speciesTransRotationEnergy[i]  = speciesTransRotationEnergy[i] * transRotationTemperature;
        speciesTransRotationEnergy[i] += thermodynamicProperties[i].hs0 / referenceEnergy;
    }
}

void MixedGas::ComputeVibrationSpecificHeatAtConstantVolume(RDouble vibrationTemperature, RDouble *speciesVibrationCv)
{
}

void MixedGas::ComputeDimensionalVibrationSpecificHeatAtConstantVolume(RDouble vibrationTemperature, RDouble *speciesVibrationCv)
{
    RDouble gasConstant = rjmk; //Unit: J/(kg*K)
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();
    //! The reference temperature and the dimensional vibration temperature.
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble Tv = vibrationTemperature * referenceTemperature;

    if (nTEnergyModel >= 1)
    {
        RDouble speciesCvtr[MAX_SPECIES_NUM], speciesCvv[MAX_SPECIES_NUM];
        if (nTEnergyModel == 1)
        {
            this->ComputeSpecificHeatByFitMethod(vibrationTemperature, speciesCvv, true);
        }
        else
        {
            this->ComputeSpecificHeatByPolynomialFitMethod(vibrationTemperature, speciesCvv, true);
        }
        this->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCvtr);
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType < 2)
            {
                speciesVibrationCv[i] = 0.0;
                continue;
            }
            //! Polyatom species.
            speciesVibrationCv[i] = speciesCvv[i] - speciesCvtr[i];
        }
    }
    else
    {
        RDouble theta = 0.0, expValue = 0.0, powValue = 0.0;
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType < 2)
            {
                speciesVibrationCv[i] = 0.0;
                continue;
            }
            //! Polyatom species.
            speciesVibrationCv[i] = 0.0;
            for (int m = 0; m < thermodynamicProperties[i].Tvs.nMode; ++ m)
            {
                theta = thermodynamicProperties[i].Tvs.Tv[m] / (Tv + SMALL);

                expValue = MIN(exp(theta), LARGE);
                powValue = theta / ((expValue - 1.0) + SMALL);
                powValue = powValue * powValue;

                speciesVibrationCv[i] += thermodynamicProperties[i].Tvs.g[m] * (gasConstant / speciesMass[i]) * expValue * powValue;
            }
        }
    }
}

void MixedGas::ComputeDimensionalVibrationEnergy(RDouble vibrationTemperature, RDouble *speciesVibrationEnergy)
{
    RDouble gasConstant = rjmk;     //! Unit: J/(kg*K).
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();

    //! Obtain the reference temperature and velocity.
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble Tv = vibrationTemperature * referenceTemperature;

    if (nTEnergyModel >= 1)
    {
        RDouble speciesCvtr[MAX_SPECIES_NUM], speciesEnthaly[MAX_SPECIES_NUM];
        if (nTEnergyModel == 1)
        {
            this->ComputeEnthalpyByFitMethod(vibrationTemperature, speciesEnthaly, true);
        }
        else
        {
            this->ComputeEnthalpyByPolynomialFitMethod(vibrationTemperature, speciesEnthaly, true);
        }
        this->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCvtr);
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType < 2)
            {
                speciesVibrationEnergy[i] = 0.0;
                continue;
            }
            //! Polyatom species.
            speciesVibrationEnergy[i] = speciesEnthaly[i] - gasConstant * Tv / speciesMass[i];
            speciesVibrationEnergy[i] -= speciesCvtr[i] * Tv + thermodynamicProperties[i].hs0;
        }
    }
    else
    {
        RDouble theta = 0.0, expValue = 0.0, powValue = 0.0;
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType < 2)
            {
                speciesVibrationEnergy[i] = 0.0;
                continue;
            }
            //! Polyatom species.
            speciesVibrationEnergy[i] = 0.0;
            for (int m = 0; m < thermodynamicProperties[i].Tvs.nMode; ++ m)
            {
                theta = thermodynamicProperties[i].Tvs.Tv[m] / (Tv + SMALL);
                expValue = MIN(exp(theta), LARGE);
                powValue = thermodynamicProperties[i].Tvs.Tv[m] / ((expValue - 1.0) + SMALL);

                speciesVibrationEnergy[i] += thermodynamicProperties[i].Tvs.g[m] * (gasConstant / speciesMass[i]) * powValue;
            }
        }
    }
}

void MixedGas::ComputeVibrationEnergy(RDouble vibrationTemperature, RDouble *speciesVibrationEnergy)
{
    //! Obtain the dimensional value of vibration energy.
    ComputeDimensionalVibrationEnergy(vibrationTemperature, speciesVibrationEnergy);

    //! Translate to non-dimensional value.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        speciesVibrationEnergy[i] /= referenceEnergy;
    }
}

void MixedGas::ComputeElectronSpecificHeatAtConstantVolume(RDouble electronTemperature, RDouble *speciesElectronCv)
{
}

void MixedGas::ComputeDimensionalElectronSpecificHeatAtConstantVolume(RDouble electronTemperature, RDouble *speciesElectronCv)
{
    RDouble gasConstant = rjmk;      //! Unit: J/(kg*K).
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();

    //! Obtain the reference temperature and compute the dimensional value of electron temperature.
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble Te = electronTemperature * referenceTemperature;

    if (nTEnergyModel >= 1)
    {
        RDouble speciesCvtr[MAX_SPECIES_NUM], speciesCvv[MAX_SPECIES_NUM];
        if (nTEnergyModel == 1)
        {
            this->ComputeSpecificHeatByFitMethod(electronTemperature, speciesCvv, false);
        }
        else
        {
            this->ComputeSpecificHeatByPolynomialFitMethod(electronTemperature, speciesCvv, false);
        }
        this->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCvtr);
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType >= 2)
            {
                speciesElectronCv[i] = 0.0;
                continue;
            }
            speciesElectronCv[i] = speciesCvv[i] - speciesCvtr[i];
        }
    }
    else
    {
        RDouble theta = 0.0, expValue = 0.0, gs0 = 0.0, gs1 = 0.0, powValue = 0.0;
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType == 0) //electron
            {
                speciesElectronCv[i] = 1.5 * gasConstant / speciesMass[i];
                continue;
            }

            //! The degeneracy of species.
            gs0 = thermodynamicProperties[i].gs0;
            gs1 = thermodynamicProperties[i].gs1;

            //! The monoaton species and polyatom species.
            theta = thermodynamicProperties[i].Tes / (Te + SMALL);
            expValue = MIN(exp(theta), LARGE);
            powValue = theta / (gs0 * expValue + gs1);
            powValue = powValue * powValue;

            speciesElectronCv[i] = (gasConstant / speciesMass[i]) * (gs0 * gs1) * expValue * powValue;
        }
    }
}

void MixedGas::ComputeDimensionalElectronEnergy(RDouble electronTemperature, RDouble *speciesElectronEnergy)
{
    RDouble gasConstant = rjmk;     //! Unit: J/(kg*K).
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();

    //! Obtain the reference temperature and compute the dimensional value of electron temperature.
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble Te = electronTemperature * referenceTemperature;

    if (nTEnergyModel >= 1)
    {
        RDouble speciesCvtr[MAX_SPECIES_NUM], speciesEnthaly[MAX_SPECIES_NUM];
        if (nTEnergyModel == 1)
        {
            this->ComputeEnthalpyByFitMethod(electronTemperature, speciesEnthaly, false);
        }
        else
        {
            this->ComputeEnthalpyByPolynomialFitMethod(electronTemperature, speciesEnthaly, false);
        }
        this->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCvtr);
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType >= 2)
            {
                speciesElectronEnergy[i] = 0.0;
                continue;
            }
            speciesElectronEnergy[i] = speciesEnthaly[i] - gasConstant * Te / speciesMass[i];
            speciesElectronEnergy[i] -= speciesCvtr[i] * Te + thermodynamicProperties[i].hs0;
        }
    }
    else
    {
        RDouble theta = 0.0, expValue = 0.0, gs0 = 0.0, gs1 = 0.0;
        for (int i = 0; i < this->numberOfSpecies; ++ i)
        {
            if (thermodynamicProperties[i].nType == 0) //electron
            {
                speciesElectronEnergy[i] = 1.5 * gasConstant * Te / speciesMass[i];
                continue;
            }

            //! The monoaton species and polyatom species.
            theta = thermodynamicProperties[i].Tes / (Te + SMALL);
            //! The degeneracy of species.
            gs0 = thermodynamicProperties[i].gs0;
            gs1 = thermodynamicProperties[i].gs1;
            expValue = MIN(exp(theta), LARGE);

            speciesElectronEnergy[i] = (gasConstant / speciesMass[i]) * thermodynamicProperties[i].Tes * (gs1 / (gs0 * expValue + gs1));
        }
    }
}

void MixedGas::ComputeElectronEnergy(RDouble electronTemperature, RDouble *speciesElectronEnergy)
{
    //! Obtain the dimensional value of electron energy.
    ComputeDimensionalElectronEnergy(electronTemperature, speciesElectronEnergy);

    //! Translate to non-dimensional value.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        speciesElectronEnergy[i] /= referenceEnergy;
    }
}

void MixedGas::ComputeSpeciesEnthalpy(RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature, RDouble *speciesEnthalpy)
{
    //! Compute energies.
    ComputeTranslationAndRotationEnergy(transRotationTemperature, this->speciesEnergyTr);
    ComputeVibrationEnergy(vibrationTemperature, this->speciesEnergyVb);
    ComputeElectronEnergy(electronTemperature, this->speciesEnergyEe);

    RDouble gasConstant = GetUniversalGasConstant();
    RDouble *speciesMass = molecularProperty->GetMolecularWeight();
    if (this->nElectronIndex >= 0)
    {
        //! Electron.
        speciesEnthalpy[nElectronIndex] = 2.5 * gasConstant * electronTemperature / speciesMass[nElectronIndex];
    }
    for (int i = 0; i <= nLastSpeciesIndex; ++ i)
    {
        //! The monoaton species and polyatom species.
        speciesEnthalpy[i] = speciesEnergyTr[i] + speciesEnergyVb[i] + speciesEnergyEe[i];
        speciesEnthalpy[i] += gasConstant * transRotationTemperature / speciesMass[i];
    }
}

RDouble MixedGas::GetMixedGasTranslationAndRotationSpecificHeatAtConstantVolume(RDouble *massFraction)
{
    //! Get dimensional value of translation-rotation specific heat at constant volume.
    RDouble transRotationCv = GetMixedGasDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(massFraction);

    //! Translate to the non-dimensional values.
    return (transRotationCv * referenceSpecificHeat);
}

RDouble MixedGas::GetMixedGasDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(RDouble *massFraction)
{
    ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCv);

    RDouble transRotationCv = 0.0;     //! The translation-rotation specific heat at constant volume.
    ComputeMixtureByMassFraction(massFraction, speciesCv, transRotationCv);

    return transRotationCv;
}

RDouble MixedGas::GetMixedGasVibrationSpecificHeatAtConstantVolume(RDouble vibrationTemperature, RDouble *massFraction)
{
    //! Get dimensional value of vibration specific heat at constant volume.
    RDouble vibrationCv = GetMixedGasDimensionalVibrationSpecificHeatAtConstantVolume(vibrationTemperature, massFraction);

    // Translate to the non-dimensional values.
    return (vibrationCv * referenceSpecificHeat);
}

RDouble MixedGas::GetMixedGasDimensionalVibrationSpecificHeatAtConstantVolume(RDouble vibrationTemperature, RDouble *massFraction)
{
    ComputeDimensionalVibrationSpecificHeatAtConstantVolume(vibrationTemperature, speciesCv);

    //! The vibration specific heat at constant volume.
    RDouble vibrationCv = 0.0;
    ComputeMixtureByMassFraction(massFraction, speciesCv, vibrationCv);

    return vibrationCv;
}

RDouble MixedGas::GetMixedGasElectronSpecificHeatAtConstantVolume(RDouble electronTemperature, RDouble *massFraction)
{
    //! Get dimensional value of electron specific heat at constant volume.
    RDouble electronCv = GetMixedGasDimensionalElectronSpecificHeatAtConstantVolume(electronTemperature, massFraction);

    //! Translate to the non-dimensional values.
    return (electronCv * referenceSpecificHeat);
}

RDouble MixedGas::GetMixedGasDimensionalElectronSpecificHeatAtConstantVolume(RDouble electronTemperature, RDouble *massFraction)
{
    this->ComputeDimensionalElectronSpecificHeatAtConstantVolume(electronTemperature, speciesCv);

    RDouble electronCv = 0.0;     //! The electron specific heat at constant volume.
    ComputeMixtureByMassFraction(massFraction, speciesCv, electronCv);

    return electronCv;
}

RDouble MixedGas::GetMixedGasTranslationAndRotationEnergy(RDouble transRotationTemperature, RDouble *massFraction)
{
    //! To obtain the translation-rotation energy.
    ComputeTranslationAndRotationEnergy(transRotationTemperature, speciesEnergyTr);

    RDouble transRotationEnergy = 0.0;
    ComputeMixtureByMassFraction(massFraction, speciesEnergyTr, transRotationEnergy);

    return transRotationEnergy;
}

RDouble MixedGas::GetMixedGasDimensionalVibrationEnergy(RDouble vibrationTemperature, RDouble *massFraction)
{
    //! To obtain the vibration energy.
    this->ComputeDimensionalVibrationEnergy(vibrationTemperature, speciesEnergyVb);

    RDouble vibrationEnergy = 0.0;
    ComputeMixtureByMassFraction(massFraction, speciesEnergyVb, vibrationEnergy);

    return vibrationEnergy;
}

RDouble MixedGas::GetMixedGasVibrationEnergy(RDouble vibrationTemperature, RDouble *massFraction)
{
    RDouble vibrationEnergy = GetMixedGasDimensionalVibrationEnergy(vibrationTemperature, massFraction);
    //! Translate to non-dimensional value.
    vibrationEnergy = vibrationEnergy / referenceEnergy;

    return vibrationEnergy;
}

RDouble MixedGas::GetMixedGasDimensionalElectronEnergy(RDouble electronTemperature, RDouble *massFraction)
{
    //! To obtain the electron energy.
    this->ComputeDimensionalElectronEnergy(electronTemperature, speciesEnergyEe);

    RDouble electronEnergy = 0.0;
    ComputeMixtureByMassFraction(massFraction, speciesEnergyEe, electronEnergy);

    return electronEnergy;
}

RDouble MixedGas::GetMixedGasElectronEnergy(RDouble electronTemperature, RDouble *massFraction)
{
    RDouble electronEnergy = GetMixedGasDimensionalElectronEnergy(electronTemperature, massFraction);
    //! Translate to non-dimensional value.
    electronEnergy = electronEnergy / referenceEnergy;

    return electronEnergy;
}

RDouble MixedGas::ComputeMixedGasEnthalpy(RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature, RDouble *massFraction)
{
    //! To obtain the enthalpies of species.
    ComputeSpeciesEnthalpy(transRotationTemperature, vibrationTemperature, electronTemperature, speciesEnthalpy);

    RDouble enthalpy = 0.0;
    ComputeMixtureByMassFraction(massFraction, speciesEnthalpy, enthalpy);

    return enthalpy;
}

RDouble MixedGas::ComputeMixedGasEnthalpy(RDouble *primitiveVars, RDouble *temperatures)
{
    using namespace IDX;
    for (int s = 0; s < numberOfSpecies; ++ s)
    {
        massFractions[s] = primitiveVars[this->nm + s];
    }
    RDouble staticEnthalpy = 0.0;
    if (this->ntmodel == 1)
    {
        staticEnthalpy = this->GetMixtureGasEnthalpy(massFractions, temperatures[ITT]);
     }
    else
    {
    }
    return staticEnthalpy;
}

RDouble MixedGas::GetTotalFormationEnthalpy(RDouble *massFraction)
{
    return 0.0;
}

RDouble MixedGas::ComputeTranslationAndRotationTemperatureByPrimitiveVariables(RDouble *primitiveVariables, RDouble electronPressure/* = 0.0*/)
{
    return 0.0;
}

RDouble MixedGas::ComputeTranslationAndRotationTemperature(RDouble *massFraction, RDouble transRotationEnergy)
{
    return 0.0;
}

RDouble MixedGas::ComputeVibrationTemperature(RDouble *massFraction, RDouble vibrationEnergy, RDouble initTemperature/* = 0.0*/)
{
    return 0.0;
}

RDouble MixedGas::ComputeElectronTemperature(RDouble *massFraction, RDouble electronEnergy, RDouble initTemperature)
{
    return 0.0;
}

RDouble MixedGas::ComputeVibrationAndElectronTemperature(RDouble *massFraction, RDouble vibrationElectronEnergy, RDouble initTemperature)
{
    return 0.0;
}

RDouble MixedGas::ComputeOneTemperatureModelTemperature(RDouble *massFraction, RDouble internalEnergy)
{
    return 0.0;
}

RDouble MixedGas::ComputeNonequilibrumTemperatureViaBisectionMethod(RDouble *massFraction, RDouble internalEnergy)
{
    return 0.0;
}

RDouble MixedGas::ComputeNonequilibrumTemperature(RDouble *massFraction, RDouble internalEnergy, RDouble initTemperature)
{
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble e1 = 0.0, e2 = 0.0, T1 = 100.0 / referenceTemperature, T2 = 30000.0 / referenceTemperature;
    RDouble ft1, ft2, Tk = 0.0, Tk1 = 0.0, precision = 1.0e-12;
    RDouble gasConstant = this->GetUniversalGasConstant();
    RDouble mass = this->GetMixedGasMolecularWeight(massFraction);
    RDouble gasR = gasConstant / mass;

    //! To determine the region where the root exists.
    int nCount = 0, nMaxStep = 500;
    if (initTemperature <= 0.0)
    {
        e1 = GetMixtureGasEnthalpy(massFraction, T1) - gasR * T1;
        e2 = GetMixtureGasEnthalpy(massFraction, T2) - gasR * T2;
        if (internalEnergy <= e1)
        {
            return T1;
        }
        else if (internalEnergy >= e2)
        {
            return T2;
        }

        //! To determine the region where the root exists.
        ft1 = e1 - internalEnergy;
        RDouble deltaT = 3000.0 / referenceTemperature;
        while (nCount <= 10)
        {
            nCount ++;
            T2 = nCount * deltaT;
            ft2 = GetMixtureGasEnthalpy(massFraction, T2) - gasR * T2 - internalEnergy;
            if (fabs(ft2) <= precision)
            {
                return T2;
            }
            else if (ft1 * ft2 <= 0.0)
            {
                break;     //! Find the region where the root exists.
            }
            T1 = T2;
            ft1 = ft2;
        }
        if (nCount > 10)     //! Failed to find the region, set the minimum temperature.
        {
            return T2;
        }
        Tk = (T1 + T2) / 2.0;
    }
    else
    {
        nMaxStep = 5;
        Tk = initTemperature;
    }

    //! To compute the vibration-electron temperature using the Newton iteration method.
    nCount = 0;
    RDouble gasMixtureCv = 0.0, deltaEnergy = 0.0;
    while (nCount < nMaxStep)
    {
        nCount++;
        ComputeConstantVolumeSpecificHeatByMassFraction(massFraction, Tk, gasMixtureCv);
        deltaEnergy = GetMixtureGasEnthalpy(massFraction, Tk) - gasR * Tk - internalEnergy;
        Tk1 = Tk - deltaEnergy / gasMixtureCv;
        if (fabs(Tk1 - Tk) <= precision)
        {
            break;
        }
        Tk = Tk1;
    }
    //printf("Number of iteration is %d\n", nCount);
    return Tk1;
}

RDouble MixedGas::GetElectronPressure(RDouble *primitiveVariables, RDouble electronTemperature)
{
    RDouble density = primitiveVariables[IDX::IR];
    RDouble gasConstant = this->GetUniversalGasConstant();
    //RDouble referenceMass = referenceParameterDimensional.GetAverageMolecularWeight();
    //! To obtain the molecular weight of species.

    //! Search the index of the electron.
    RDouble electronMassFraction = 0.0;
    //RDouble electronMass = 5.486e-7 / referenceMass;
    if (nElectronIndex >= 0)
    {
        electronMassFraction = primitiveVariables[nm + nElectronIndex];
    }

    //! Compute electron pressure.
    RDouble electronPressure = (density * gasConstant * electronTemperature) * electronMassFraction / electronMass;
    return electronPressure;
}

void MixedGas::ComputeSpeciesTranslationAndRotationHeatConductivity(RDouble *viscosity, RDouble *conductivity)
{

}

void MixedGas::ComputeSpeciesVibrationHeatConductivity(RDouble *viscosity, RDouble *speciesVibrationCv, RDouble *conductivity)
{

}

void MixedGas::ComputeSpeciesElectronHeatConductivity(RDouble *viscosity, RDouble *speciesElectronCv, RDouble *conductivity)
{

}

RDouble MixedGas::GetMixedGasTranslationAndRotationHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble electronTemperature)
{
    return 0.0;
}

RDouble MixedGas::GetMixedGasVibrationHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature)
{
    return 0.0;
}

RDouble MixedGas::GetMixedGasElectronHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble electronTemperature)
{
    return 0.0;
}

void MixedGas::ComputeMixedGasTranslationVibrationElectronHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature, RDouble &transRotationConductivity, RDouble &vibrationConductivity, RDouble &electronConductivity)
{

}

RDouble MixedGas::GetMixedGasHeatConductivityWithWilkeFormula(RDouble *moleFraction, RDouble *conductivity, RDouble *phi)
{
    RDouble lambda = 0.0;
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        lambda += moleFraction[i] * conductivity[i] / phi[i];
    }
    return lambda;
}

void MixedGas::ComputeSpeciesViscosityWithCurveFitMethod(RDouble transRotationTemperature, RDouble electronTemperature, RDouble *speciesViscosity)
{
    //! Obtain the reference viscosity.
    RDouble referenceViscosity = this->referenceParameterDimensional.GetViscosity();

    //! Obtain the dimensional value of viscosity.
    ComputeSpeciesDimensionalViscosityWithCurveFitMethod(transRotationTemperature, electronTemperature, speciesViscosity);

    //! Translate to the non-dimensional value.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        speciesViscosity[i] /= referenceViscosity;
    }
}

void MixedGas::ComputeSpeciesDimensionalViscosityWithCurveFitMethod(RDouble transRotationTemperature, RDouble electronTemperature, RDouble *speciesViscosity)
{
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble dimensionalTransRotationTemperature, dimensionalElectronTemperature;
    dimensionalTransRotationTemperature = transRotationTemperature * referenceTemperature;
    dimensionalElectronTemperature = electronTemperature * referenceTemperature;
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.

    //! Obtain the coefficients of the curve fitting.
    RDouble *coefA = CurveFits->GetBlotterCurveFittingCoefficientA();
    RDouble *coefB = CurveFits->GetBlotterCurveFittingCoefficientB();
    RDouble *coefC = CurveFits->GetBlotterCurveFittingCoefficientC();
    RDouble *coefD = CurveFits->GetBlotterCurveFittingCoefficientD();
    RDouble *coefE = CurveFits->GetBlotterCurveFittingCoefficientE();

    RDouble lnt1, lnt2, lnt3, lnt4, tmp;
    lnt1 = log(dimensionalTransRotationTemperature);     //! Polyatom and monoatom.
    lnt2 = lnt1 * lnt1;
    lnt3 = lnt1 * lnt2;
    lnt4 = lnt1 * lnt3;
    for (int i = 0; i <= nLastSpeciesIndex; ++ i)
    {
        //! To compute viscosity using the curve fitting method.
        tmp = coefA[i] * lnt4 + coefB[i] * lnt3 + coefC[i] * lnt2 + coefD[i] * lnt1 + coefE[i];

        //! Compute the dimensional value of viscosity.
        speciesViscosity[i] = 0.1 * exp(tmp);
    }
    if (nElectronIndex >= 0)// electron.
    {
        lnt1 = log(dimensionalElectronTemperature);
        lnt2 = lnt1 * lnt1;
        lnt3 = lnt1 * lnt2;
        lnt4 = lnt1 * lnt3;
        tmp = coefA[nElectronIndex] * lnt4 + coefB[nElectronIndex] * lnt3 + coefC[nElectronIndex] * lnt2 + coefD[nElectronIndex] * lnt1 + coefE[nElectronIndex];
        speciesViscosity[nElectronIndex] = 0.1 * exp(tmp);
    }
}

void MixedGas::ComputeSpeciesViscosityWithLennardJonesMethod(RDouble transRotationTemperature, RDouble electronPressure, RDouble *speciesViscosity)
{
    //! Obtain the reference viscosity.
    RDouble referenceViscosity = this->referenceParameterDimensional.GetViscosity();

    //! Obtain the dimensional value of viscosity.
    ComputeSpeciesDimensionalViscosityWithLennardJonesMethod(transRotationTemperature, electronPressure, speciesViscosity);

    //! Translate to the non-dimensional value.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        speciesViscosity[i] /= referenceViscosity;
    }
}

void MixedGas::ComputeSpeciesDimensionalViscosityWithLennardJonesMethod(RDouble transRotationTemperature, RDouble electronPressure, RDouble *speciesViscosity)
{
    //! Obtain the reference temperature.
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble dimensionalTransRotationTemperature = transRotationTemperature * referenceTemperature;

    //! The average area of collision between species.
    RDouble *collisionArea = new RDouble[this->numberOfSpecies];
    ComputeAverageCollisionAreaOmega22(transRotationTemperature, electronPressure, collisionArea);

    //! Obtain the molecular weights of species.
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();
    //! Translate to the dimensional value.
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        speciesViscosity[i] = 8.4411e-8 * PI * sqrt(speciesMass[i] * dimensionalTransRotationTemperature) / collisionArea[i];
    }
    delete [] collisionArea;
}

void MixedGas::ComputeAverageCollisionAreaOmega22(RDouble transRotationTemperature, RDouble electronPressure, RDouble *piOmega)
{
    //! Obtain the reference temperature and reference pressure.
    RDouble referenceTemperature = this->referenceParameterDimensional.GetTemperature();
    RDouble referencePressure = this->referenceParameterDimensional.GetPressure();
    RDouble dimensionalTransRotationTemperature, temperature, lnt, tmp;
    dimensionalTransRotationTemperature = transRotationTemperature * referenceTemperature;
    temperature = dimensionalTransRotationTemperature / 1000.0;
    lnt = log(dimensionalTransRotationTemperature);

    RDouble p0, pe, pem, pressure, lnpe;
    pe  = electronPressure * referencePressure;
    p0  = 101325.0; //Pa
    pem = p0 * 0.0975 * pow(temperature, 4);
    pressure = p0 / pe;

    int *ionType = this->molecularProperty->GetIonTypeOfSpecies();
    RDouble *coefAss, *coefBss, *coefCss, *coefDss;
    coefAss = this->CurveFits->GetOmegaCurveFittingCoefficientA();
    coefBss = this->CurveFits->GetOmegaCurveFittingCoefficientB();
    coefCss = this->CurveFits->GetOmegaCurveFittingCoefficientC();
    coefDss = this->CurveFits->GetOmegaCurveFittingCoefficientD();

    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        if (ionType[i] == 0)     //! The atom or molecular.
        {
            lnpe = 1.0;
        }
        else
        {
            if (fabs(pe - pem) <= 1.0e-8)
            {
                lnpe = 1.0;
            }
            else     //! Modify the coefficient.
            {
                lnpe = 0.5 * log(2.09e-2 * pow(temperature, 4) * pressure + 1.52 * pow(temperature, 8.0 / 3.0) * pow(pressure, 2.0 / 3.0));
            }
        }
        tmp = coefAss[i] * pow(lnt, 3) + coefBss[i] * pow(lnt, 2) + coefCss[i] * lnt + coefDss[i];
        piOmega[i] = 1.0e-20 * lnpe * exp(tmp);
    }
}

RDouble MixedGas::GetMixedGasViscosityWithWilkeFormula(RDouble *moleFraction, RDouble *viscosity, RDouble *phi)
{
    RDouble mixedViscosity = 0.0;
    for (int i = 0; i < this->numberOfSpecies; ++ i)
    {
        mixedViscosity += moleFraction[i] * viscosity[i] / phi[i];
    }
    return mixedViscosity;
}

void MixedGas::ComputeSpeciesMassDiffusionCoefficient(RDouble *massFraction, RDouble viscosityLaminar, RDouble viscosityTurbulence, RDouble *speciesDiffusionCoef)
{
    //! The function returns the diffusion coefficients of each species which saved in rDs = (1 - Cs)/(1- Xs) * (mul/Scl + mut/Sct).
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();
    int *ionType = this->molecularProperty->GetIonTypeOfSpecies();

    //! To compute the mole fractions of species.
    RDouble *speciesMass = this->molecularProperty->GetMolecularWeight();
    this->ComputeMoleFractionByMassFraction(massFraction, moleFractions);

    //! Obtain the Schmidt number.
    RDouble *oSchmidtLaminar = GetLaminarSchmidtNumberReciprocal();
    RDouble *oSchmidtTurbulence = GetTurbulentSchmidtNumberReciprocal();

    RDouble tmp1 = 0.0, tmp2 = 0.0;
    int nIndex = -1;
    for(int i = 0; i < this->numberOfSpecies; ++ i)
    {
        if (thermodynamicProperties[i].nType == 0)
        {
            nIndex = i;      //! To save the index of the electron.
            continue;      //! Ignore the electron.
        }

        if (1.0 - moleFractions[i] <= 1.0e-20)
        {
            speciesDiffusionCoef[i] = 0.0;
        }
        else
        {
            speciesDiffusionCoef[i] = (1.0 - massFraction[i]) / (1.0 - moleFractions[i]) * (viscosityLaminar * oSchmidtLaminar[i] + viscosityTurbulence * oSchmidtTurbulence[i]);
        }

        if (ionType[i] == 0)
        {
            continue;     //! Ignore the atom and molecular species.
        }
        else
        {
            tmp1 += speciesDiffusionCoef[i] * massFraction[i] / speciesMass[i];
            tmp2 += massFraction[i] / speciesMass[i];
        }
    }

    //! Compute the mass diffusion coefficient of electron.
    if (nIndex >= 0)
    {
        speciesDiffusionCoef[nIndex] = tmp1 / MAX(tmp2, 1.0e-20);
    }
}

RDouble MixedGas::GetMassDiffusionCoefficientWithCLNModel(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble electronTemperature)
{
    return 0.0;
}

void MixedGas::GetElectronPressurePartialDerivatives(RDouble *primitiveVars, RDouble electronTemperature, RDouble &gama, RDouble *speciesGama)
{

}

void MixedGas::GetPressurePartialDerivatives(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble electronTemperature,
                                             RDouble &alpha, RDouble &beta, RDouble *speciesBeta)
{
}

void MixedGas::GetMultiTemperatureModelPartialDerivatives(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble electronTemperature,
                                                          RDouble &alpha, RDouble &beta, RDouble &gama, RDouble *speciesBeta, RDouble *speciesGama)
{
}

void MixedGas::GetMultiTemperatureModelPartialDerivativesR(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble electronTemperature,
                                                           RDouble &alpha, RDouble &beta, RDouble &gama, RDouble *speciesBeta, RDouble *speciesGama,
                                                           RDouble squareVelocity, RDouble *speciesEtrs, RDouble totalCvtr, RDouble totalCvv, RDouble totalCve)
{

}

RDouble MixedGas::ComputePressureDerivativeToTotalDensity(RDouble *primitiveVars)
{
    return 0.0;
}

RDouble MixedGas::GetMixedGasEnthalpy(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature)
{
    int nChemical = GetChemicalType();
    RDouble staticEnthaly = 0.0;
    if (nChemical == 0)
    {
        staticEnthaly = 3.5 * primitiveVars[IDX::IP] / primitiveVars[IDX::IR];
    }
    else
    {
        RDouble temperatures[3] = {transRotationTemperature, vibrationTemperature, electronTemperature};
        staticEnthaly = ComputeMixedGasEnthalpy(primitiveVars, temperatures);
    }

    //delete []massFraction;
    return staticEnthaly;
}

RDouble MixedGas::ComputeMolecularWeightReciprocal(RDouble *prim)
{
    RDouble omav = 1.0;
    ComputeMolecularWeightReciprocal(prim, omav);
    return omav;
}

RDouble MixedGas::ComputeMolecularWeight(RDouble *prim)
{
    return 0.0;
}

void MixedGas::ComputeMixtureGasHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature,
                                                 RDouble &transRotationConductivity, RDouble &vibrationConductivity, RDouble &electronConductivity)
{
    //int nTemperatureModel = GetTemperatureModel();
    if (ntmodel > 1)
    {
    }
    else     //! One-temperature model.
    {
        //! Compute heat conductivity of mixture gas.
        transRotationConductivity = ComputeMixtureGasHeatConductivityByWassilewaFormula(primitiveVariables, transRotationTemperature);
        vibrationConductivity = 0.0;
        electronConductivity = 0.0;
    }
}

void MixedGas::ComputeMixtureGasSpecificHeatAtConstantVolume(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature,
                                                             RDouble &cvtr, RDouble &cvv, RDouble &cve)
{

}

void MixedGas::ComputeEnergySourceDerivatives(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
                                              RDouble ***chemicalSourceDerivatives, RDouble ***energySourceDerivatives, int nCellNumber)
{
}

void MixedGas::ComputeEnergySourceDiagonalDerivatives(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
                                                      RDouble ***chemicalSourceDerivatives, RDouble **energySourceDerivatives, int nCellNumber)
{
}

RDouble MixedGas::GetMaximumSpeciesMassDiffusionCoef(RDouble *primitiveVars, RDouble viscosityLaminar, RDouble viscosityTurbulence)
{
    return 0.0;     //! Returned value.
}

void MixedGas::GetSpeciesMassDiffusionCoef(RDouble *primitiveVars, RDouble viscosityLaminar, RDouble viscosityTurbulence, RDouble *speciesDiffusionCoef)
{
    //! The function returns the diffusion coefficients of each species which saved in rDs = (1 - Cs)/(1- Xs) * (mul/Scl + mut/Sct).
    for (int i = 0; i < numberOfSpecies; ++ i)
    {
        massFractions[i] = primitiveVars[nm + i];
    }

    ComputeSpeciesMassDiffusionCoefficient(massFractions, viscosityLaminar, viscosityTurbulence, speciesDiffusionCoef);
}

void MixedGas::GetTemperatureAndCP(RDouble *prim, RDouble &Tm, RDouble &cp)
{
    int nchem = GetChemicalType();
    int nm = GetNSEquationNumber();

    //! Specific heat ratio
    RDouble gama = 1.4;
    RDouble &ro = prim[0];
    RDouble &pm = prim[4];
    //! The non-dimensional universal gas constant.
    RDouble R = GetCoefficientOfStateEquation();
    //! The constant pressure specific heat of perfect gas.
    cp = gama / (gama - 1.0) * R;
    //! Compute the temperature using the state equation of perfect gas.
    Tm = pm / (ro * R);

    //! Mixture gas condition.
    if (nchem == 1)
    {
        //! Save the number of species
        int ns = GetNumberOfSpecies();
        RDouble *fs = new RDouble[ns];
        for (int i = 0; i < ns; i ++)
        {
            fs[i] = prim[i + nm];     //! Obtain the mass fraction of each species.
        }

        //! Obtain the average molecular weight of mixture gas.
        RDouble M = GetMixedGasMolecularWeight(fs);
        R = GetUniversalGasConstant();     //! The non-dimensional universal gas constant.
        //! Compute the temperature using the state equation of mixture gas.
        Tm = pm * M / (ro * R);
        
        //! Compute constant pressure specific heat of each species.
        ComputeConstantPressureSpecificHeatByMassFraction(fs, Tm, cp);

        //! Remove the memory space.
        delete [] fs;    fs = nullptr;
    }
}

RDouble MixedGas::GetTemperature(RDouble *prim)
{
    int nchem = GetChemicalType();
    int nm = GetNSEquationNumber();
    //! Obtain the density and pressure.
    RDouble &ro = prim[0];
    RDouble &pm = prim[4];
    RDouble R = GetCoefficientOfStateEquation();
    //! Compute the temperature using the state equation of mixture gas.
    RDouble Tm = pm / (ro * R);

    if (nchem == 1)
    {
        //! Save the number of species.
        int ns = GetNumberOfSpecies();
        RDouble *fs = new RDouble[ns];
        //! Obtain the mass fraction of each species.
        for (int i = 0; i < ns; i ++)
        {
            fs[i] = prim[i + nm];
        }

        //! Obtain the average molecular weight of mixture gas.
        RDouble M = GetMixedGasMolecularWeight(fs);
        R = GetUniversalGasConstant();     //! The non-dimensional universal gas constant.
        //! Compute the temperature using the state equation of mixture gas.
        Tm = pm * M / (ro * R);
        //! Remove the memory space.
        delete [] fs;    fs = nullptr;
    }
    return Tm;
}

void MixedGas::GetEverySpeciesEnthalpy(RDouble trTemperature, RDouble vTemperature, RDouble eTemperature, RDouble *hs)
{
    if (ntmodel > 1)      //! Multi-temperature model.
    {
    }
    else      //! Single temperature model.
    {
        ComputeSpeciesEnthalpy(trTemperature, hs);
    }
}

RDouble MixedGas::ComputeMolecularWeightReciprocalWithoutElectron(RDouble *prim, RDouble &ce_div_me)
{
    RDouble massFraction[MAX_SPECIES_NUM];
    for (int i = 0; i < numberOfSpecies; i ++)
    {
        massFraction[i] = prim[i + nm];
    }

    RDouble M = GetMixedGasMolecularWeightReciprocalWithoutElectron(massFraction, ce_div_me);
    return M;
}

RDouble MixedGas::GetMixedGasTranslationAndRotationEnergy(RDouble *prim, RDouble Ttr)
{
    RDouble massFraction[MAX_SPECIES_NUM];
    for (int i = 0; i < numberOfSpecies; i ++)
    {
        massFraction[i] = prim[nm + i];
    }

    RDouble Etr = GetMixedGasTranslationAndRotationEnergy(Ttr, massFraction);
    return Etr;
}

RDouble MixedGas::GetMixedGasVibrationEnergy(RDouble *prim, RDouble Tv)
{
    RDouble massFraction[MAX_SPECIES_NUM];
    for (int i = 0; i < numberOfSpecies; ++ i)
    {
        massFraction[i] = prim[nm + i];
    }

    RDouble vibrationEnergy = GetMixedGasVibrationEnergy(Tv, massFraction);
    return vibrationEnergy;
}

RDouble MixedGas::GetMixedGasElectronEnergy(RDouble *prim, RDouble Te)
{
    RDouble massFraction[MAX_SPECIES_NUM];
    for (int i = 0; i < numberOfSpecies; ++ i)
    {
        massFraction[i] = prim[nm + i];
    }

    RDouble electronEnergy = GetMixedGasElectronEnergy(Te, massFraction);
    return electronEnergy;
}

RDouble MixedGas::GetMixedGasTransAndRotatSpecHeatAtConstVolume(RDouble *prim)
{
    return 0.0;
}

RDouble MixedGas::GetMixedGasVibrationSpecificHeatAtConstantVolume(RDouble *prim, RDouble Tv)
{
    RDouble massFraction[MAX_SPECIES_NUM];
    for (int i = 0; i < numberOfSpecies; i ++)
    {
        massFraction[i] = prim[nm + i];
    }

    RDouble cvv = GetMixedGasVibrationSpecificHeatAtConstantVolume(Tv, massFraction);
    return cvv;
}

RDouble MixedGas::GetMixedGasElectronSpecificHeatAtConstantVolume(RDouble *prim, RDouble Te)
{
    RDouble massFraction[MAX_SPECIES_NUM];
    for (int i = 0; i < numberOfSpecies; i ++)
    {
        massFraction[i] = prim[nm + i];
    }

    RDouble cve = GetMixedGasElectronSpecificHeatAtConstantVolume(Te, massFraction);
    return cve;
}

RDouble MixedGas::ComputeFrozenSoundSpeed(RDouble *primitiveVars)
{
    return 0.0;
}

void MixedGas::ComputeSpecificHeat(RDouble *primitiveVars, RDouble *temperatures, RDouble *spciesCvtr, RDouble *spciesCvv, RDouble &totalCvtr, RDouble &totalCvv)
{

}

void MixedGas::ComputeEnergy(RDouble *primitiveVars, RDouble *temperatures, RDouble *spciesEtr, RDouble *spciesEv, RDouble *spciesEnthalpy, RDouble &totalEtr, RDouble &totalEv, RDouble &totalEnthalpy)
{
}

void MixedGas::ComputePressureDerivatives(RDouble *primitiveVars, RDouble *temperatures, RDouble &alpha, RDouble &beta, RDouble &eta, RDouble *speciesBeta)
{
}

void MixedGas::ComputeSpecificHeatByFitMethod(RDouble temperature, RDouble *cvSpecies, bool isVibrationModel/* = true*/)
{
    RDouble gasConstant = rjmk;     //! Unit: J/(kg*K).
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();
    RDouble referenceTemperature = referenceParameterDimensional.GetTemperature();
    RDouble *coef, cp;
    RDouble T1, T2, T3, T4;
    T1 = temperature * referenceTemperature;
    //T1 = MAX(T1, 50.0);
    T2 = T1 * T1;
    T3 = T1 * T2;
    T4 = T1 * T3;

    if (isVibrationModel)
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            if (thermodynamicProperties[ispecies].nType < 2)
            {
                cvSpecies[ispecies] = 0.0;
                continue;
            }
            coef = thermodynamicManager->GetEnthalpyFitCoefficient(ispecies);
            cp = coef[1] + 2 * coef[2] * T1 + 3 * coef[3] * T2 + 4 * coef[4] * T3 + 5 * coef[5] * T4;
            cvSpecies[ispecies] = (cp - 1.0) * gasConstant / speciesMass[ispecies];
        }
    }
    else
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            if (thermodynamicProperties[ispecies].nType >= 2)
            {
                cvSpecies[ispecies] = 0.0;
                continue;
            }
            coef = thermodynamicManager->GetEnthalpyFitCoefficient(ispecies);
            cp = coef[1] + 2 * coef[2] * T1 + 3 * coef[3] * T2 + 4 * coef[4] * T3 + 5 * coef[5] * T4;
            cvSpecies[ispecies] = (cp - 1.0) * gasConstant / speciesMass[ispecies];
        }
    }
}

void MixedGas::ComputeEnthalpyByFitMethod(RDouble temperature, RDouble *hSpecies, bool isVibrationModel/* = true*/)
{
    RDouble gasConstant = rjmk;     //! Unit: J/(kg*K).
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();
    RDouble referenceTemperature = referenceParameterDimensional.GetTemperature();
    RDouble *coef, hs;
    RDouble T1, T2, T3, T4, T5;
    T1 = temperature * referenceTemperature;
    //T1 = MAX(T1, 50.0);
    T2 = T1 * T1;
    T3 = T1 * T2;
    T4 = T1 * T3;
    T5 = T1 * T4;

    if (isVibrationModel)
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            if (thermodynamicProperties[ispecies].nType < 2)
            {
                hSpecies[ispecies] = 0.0;
                continue;
            }
            coef = thermodynamicManager->GetEnthalpyFitCoefficient(ispecies);
            hs = coef[0] + coef[1] * T1 + coef[2] * T2 + coef[3] * T3 + coef[4] * T4 + coef[5] * T5;
            hSpecies[ispecies] = hs * gasConstant / speciesMass[ispecies];
        }
    }
    else
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            if (thermodynamicProperties[ispecies].nType >= 2)
            {
                hSpecies[ispecies] = 0.0;
                continue;
            }
            coef = thermodynamicManager->GetEnthalpyFitCoefficient(ispecies);
            hs = coef[0] + coef[1] * T1 + coef[2] * T2 + coef[3] * T3 + coef[4] * T4 + coef[5] * T5;
            hSpecies[ispecies] = hs * gasConstant / speciesMass[ispecies];
        }
    }
}

void MixedGas::ComputeSpecificHeatByPolynomialFitMethod(RDouble temperature, RDouble *cvSpecies, bool isVibrationModel/* = true*/)
{
    RDouble gasConstant = rjmk;     //! Unit: J/(kg*K).
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();
    RDouble referenceTemperature = referenceParameterDimensional.GetTemperature();
    RDouble cp, T1, T2, T3, T4;
    T1 = temperature * referenceTemperature;
    //T1 = MAX(T1, 1.0);
    T2 = T1 * T1;
    T3 = T1 * T2;
    T4 = T1 * T3;

    //! Obtain the interpolation coefficients for the polynomial.
    this->GetLinearlyAveragedPolynomialCoefficientsT2(T1, this->polynomialCoefT2);

    if (isVibrationModel)
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            if (thermodynamicProperties[ispecies].nType < 2)
            {
                cvSpecies[ispecies] = 0.0;
                continue;
            }
            cp = polynomialCoefT2[ispecies][1] + 2 * polynomialCoefT2[ispecies][2] * T1 + 3 * polynomialCoefT2[ispecies][3] * T2 + 4 * polynomialCoefT2[ispecies][4] * T3 + 5 * polynomialCoefT2[ispecies][5] * T4;
            cvSpecies[ispecies] = (cp - 1.0) * gasConstant / speciesMass[ispecies];
        }
    }
    else
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            if (thermodynamicProperties[ispecies].nType >= 2)
            {
                cvSpecies[ispecies] = 0.0;
                continue;
            }
            cp = polynomialCoefT2[ispecies][1] + 2 * polynomialCoefT2[ispecies][2] * T1 + 3 * polynomialCoefT2[ispecies][3] * T2 + 4 * polynomialCoefT2[ispecies][4] * T3 + 5 * polynomialCoefT2[ispecies][5] * T4;
            cvSpecies[ispecies] = (cp - 1.0) * gasConstant / speciesMass[ispecies];
        }
    }
}

void MixedGas::ComputeEnthalpyByPolynomialFitMethod(RDouble temperature, RDouble *hSpecies, bool isVibrationModel/* = true*/)
{
    RDouble gasConstant = rjmk;     //! Unit: J/(kg*K).
    RDouble *speciesMass = molecularProperty->GetMolecularWeightDimensional();
    //! nTypes denotes the types of species which is used to identity the electron, monoatom and polyatom.
    Thermo_Param *thermodynamicProperties = this->molecularProperty->GetThermodynamicProperties();
    RDouble referenceTemperature = referenceParameterDimensional.GetTemperature();
    RDouble hs, T1, T2, T3, T4, T5;
    T1 = temperature * referenceTemperature;
    //T1 = MAX(T1, 1.0);
    T2 = T1 * T1;
    T3 = T1 * T2;
    T4 = T1 * T3;
    T5 = T1 * T4;

    //! Obtain the interpolation coefficients for the polynomial.
    this->GetLinearlyAveragedPolynomialCoefficientsT2(T1, this->polynomialCoefT2);

    if (isVibrationModel)
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            if (thermodynamicProperties[ispecies].nType < 2)
            {
                hSpecies[ispecies] = 0.0;
                continue;
            }
            hs = polynomialCoefT2[ispecies][0] + polynomialCoefT2[ispecies][1] * T1 + polynomialCoefT2[ispecies][2] * T2 + polynomialCoefT2[ispecies][3] * T3 + polynomialCoefT2[ispecies][4] * T4 + polynomialCoefT2[ispecies][5] * T5;
            hSpecies[ispecies] = hs * gasConstant / speciesMass[ispecies];
        }
    }
    else
    {
        for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
        {
            if (thermodynamicProperties[ispecies].nType >= 2)
            {
                hSpecies[ispecies] = 0.0;
                continue;
            }
            hs = polynomialCoefT2[ispecies][0] + polynomialCoefT2[ispecies][1] * T1 + polynomialCoefT2[ispecies][2] * T2 + polynomialCoefT2[ispecies][3] * T3 + polynomialCoefT2[ispecies][4] * T4 + polynomialCoefT2[ispecies][5] * T5;
            hSpecies[ispecies] = hs * gasConstant / speciesMass[ispecies];
        }
    }
}

void MixedGas::ComputeTemperatureFromTotalTemperature(RDouble *massF, RDouble totalTemperature, RDouble mach, RDouble &Temperature, RDouble &gama, int Dimensional)
{
    RDouble Tmin = 10.0;
    RDouble Tmax = 50000.0;
    RDouble mach2 = mach * mach;
    RDouble totalTemperatureLocal;
    RDouble referenceTemperature = referenceParameterDimensional.GetTemperature();
    if (Dimensional == 0)
    {
        Tmin = 10.0 / referenceTemperature;
        Tmax = 50000.0 / referenceTemperature;
    }
    int count = 0;
    while (count < 100)
    {
        Temperature = 0.5 * (Tmin + Tmax);
        ComputeGama(massF, Temperature, gama, Dimensional);
        totalTemperatureLocal = Temperature * (1.0 + 0.5 * (gama - 1.0) * mach2);
        if(ABS(totalTemperatureLocal - totalTemperature) / totalTemperature < 1.0e-8) return;
        if (totalTemperatureLocal > totalTemperature)
        {
            Tmax = Temperature;
        }
        else
        {
            Tmin = Temperature;
        }
        count++;
    }
}

void ChemCompressData(DataContainer *&cdata)
{
    gas->CompressData(cdata);
}

void ChemDecompressData(DataContainer *cdata)
{
    gas->DecompressData(cdata);
}

void MixedGas::SetMaximumSpecies(RDouble *fs)
{
    RDouble *Ms = this->molecularProperty->GetMolecularWeightDimensional();
    string *species_name = this->molecularProperty->GetNameOfSpecies();

    char    elementName[4] = {'O', 'N', 'C', 'H'}; 
    RDouble elementWeight[4] = {1.599940E-2, 1.400674E-2, 1.400674E-2, 1.0E-3};
    RDouble elementMaximum[4] = {0.0, 0.0, 0.0, 0.0};

    int coef = 0;
    RDouble elementMaximumC1;

    RDouble *xs = new RDouble[this->numberOfSpecies]();

    for (int i = 0; i < this->numberOfSpecies; i ++)
    {
        maximumSpecies[i] = 1.0;    //ini
        xs[i] = fs[i] / Ms[i];
    }

    if (nSpeciesLimit == 0) return;

    //! Compute the maximum of each element.
    for (int n = 0; n < 4; n ++)
    {
        for (int i = 0; i < this->numberOfSpecies; i ++)
        {
            if (FindChemicalElement(species_name[i], elementName[n], coef))
            {
                elementMaximum[n] += xs[i] * coef * elementWeight[n];
            }
        }
    }

    //! Compute the maximum of each Specie.
    for (int i = 0; i < this->numberOfSpecies; i ++)
    {
        for (int n = 0; n < 4; n ++)
        {
            if (FindChemicalElement(species_name[i], elementName[n], coef))
            {
                elementMaximumC1 = Ms[i] * elementMaximum[n] / (elementWeight[n] * coef);
                maximumSpecies[i] = MIN(maximumSpecies[i], elementMaximumC1);
            }
        }
    }
    delete [] xs;    xs = NULL;
}

void ChemicalMXDQ(RDouble *primitiveVars, RDouble *temperature, RDouble nx, RDouble ny, RDouble nz, RDouble faceArea, RDouble moveVelocity,
                  RDouble *deltaQ, RDouble *flux, int nNSEquation, int nLaminar, int nTemperatureModel, RDouble radius, int iSign, int nElectronIndex)
{
    using namespace IDX;
    using namespace GAS_SPACE;
    //! The number of species.
    int nSpecies = nLaminar - nNSEquation + 1;
    RDouble alpha = 0.0, beta = 0.0, gama = 0.0, *speciesBeta;
    speciesBeta = new RDouble[nSpecies];

    RDouble uVelocity = primitiveVars[IU];
    RDouble vVelocity = primitiveVars[IV];
    RDouble wVelocity = primitiveVars[IW];
    RDouble squareVelocity, absoluteVelocity, relativeVelocity;

    //! The absolute velocity.
    absoluteVelocity = nx * uVelocity + ny * vVelocity + nz * wVelocity;
    //! The relative velocity.
    relativeVelocity = absoluteVelocity - moveVelocity;
    //! Compute squared velocity.
    squareVelocity = uVelocity * uVelocity + vVelocity * vVelocity + wVelocity * wVelocity;

    //! Obtain the temperatures.
    RDouble transRotationTemperature = temperature[ITT];
    RDouble vibrationTemperature = transRotationTemperature;
    RDouble electronTemperature = transRotationTemperature;
    if (nTemperatureModel > 1)
    {
        vibrationTemperature = temperature[ITV];
        electronTemperature = vibrationTemperature;
        if (nTemperatureModel == 3)
        {
            electronTemperature = temperature[ITE];
        }
    }

    //! Compute the total enthalpy.
    RDouble totalEnthalpy = gas->GetMixedGasEnthalpy(primitiveVars, transRotationTemperature, vibrationTemperature, electronTemperature);
    totalEnthalpy += 0.5 * squareVelocity;

    //! Compute the coefficient b1 and b2 (see formula A.5 in PHengLEI Theory manual).
    RDouble b1 = nx * deltaQ[IRU] + ny * deltaQ[IRV] + nz * deltaQ[IRW] - absoluteVelocity * deltaQ[IR];
    RDouble dh = uVelocity * deltaQ[IRU] + vVelocity * deltaQ[IRV] + wVelocity * deltaQ[IRW] - deltaQ[IRE];
    RDouble b2 = 0.0, b3 = 0.0, electronPressure = 0.0;

    if (nTemperatureModel == 1)
    {
        //! Compute the coefficients in the chemical reaction.
        gas->ComputeCoefficientInMXDQ(primitiveVars, alpha, beta, speciesBeta);

        b2 = beta * deltaQ[IR] - alpha * dh;
        for (int s = 0; s < nSpecies - 1; ++ s)
        {
            b2 += speciesBeta[s] * deltaQ[nNSEquation + s];
        }
    }
    else    //! Multi-temperature model.
    {
        //! Compute the pressure of electron.
        electronPressure = gas->GetElectronPressure(primitiveVars, electronTemperature);
        RDouble *speciesGama = new RDouble[nSpecies];

        //! Compute partial derivatives of pressure and electron pressure to conservative variables.
        gas->GetMultiTemperatureModelPartialDerivatives(primitiveVars, transRotationTemperature, electronTemperature,
            alpha, beta, gama, speciesBeta, speciesGama);

        if (nTemperatureModel == 2)
        {
            dh += deltaQ[nLaminar + 1];
            b2 += gama * deltaQ[nLaminar + 1];
            b3 += gama * deltaQ[nLaminar + 1];
        }
        else if (nTemperatureModel == 3)
        {
            b2 += gama * deltaQ[nLaminar + 1];
            b3 += gama * deltaQ[nLaminar + 2];
            dh += deltaQ[nLaminar + 1] + deltaQ[nLaminar + 2];
        }
        b2 += beta * deltaQ[IR] - alpha * dh;

        for (int iSpices = 0; iSpices < nSpecies - 1; ++ iSpices)
        {
            b2 += speciesBeta[iSpices] * deltaQ[nNSEquation + iSpices];
            b3 += speciesGama[iSpices] * deltaQ[nNSEquation + iSpices];
        }

        delete [] speciesGama;    speciesGama = nullptr;
    }

    //! Comoute the fluxes called M*dQ.
    flux[IR ] = relativeVelocity * deltaQ[IR ] + b1;
    flux[IRU] = relativeVelocity * deltaQ[IRU] + b1 * uVelocity     + b2 * nx;
    flux[IRV] = relativeVelocity * deltaQ[IRV] + b1 * vVelocity     + b2 * ny;
    flux[IRW] = relativeVelocity * deltaQ[IRW] + b1 * wVelocity     + b2 * nz;
    flux[IRE] = relativeVelocity * deltaQ[IRE] + b1 * totalEnthalpy + b2 * absoluteVelocity;
    for (int m = nNSEquation; m < nLaminar; ++ m)
    {
        flux[m] = relativeVelocity * deltaQ[m] + b1 * primitiveVars[m];
    }

    //! Modify the fluxes of multi-temperature model.
    if (nTemperatureModel == 2) 
    {
        flux[nLaminar + 1] = relativeVelocity * deltaQ[nLaminar + 1] + b1 * (primitiveVars[nLaminar + 1] + electronPressure / primitiveVars[IR]) + b3 * absoluteVelocity;
    }
    else if (nTemperatureModel == 3)
    {
        flux[nLaminar + 1] = relativeVelocity * deltaQ[nLaminar + 1] + b1 * primitiveVars[nLaminar + 1];
        flux[nLaminar + 2] = relativeVelocity * deltaQ[nLaminar + 2] + b1 * (primitiveVars[nLaminar + 2] + electronPressure / primitiveVars[IR]) + b3 * absoluteVelocity;
    }

    //! To add rM * dQ that is the product of the invscid spectrum radius rM and increment dQ.
    RDouble spectrumRadiusTerm = iSign * radius / faceArea;
    for (int m = 0; m < nLaminar; ++ m)
    {
        flux[m] += spectrumRadiusTerm * deltaQ[m];
        flux[m] *= 0.5 * faceArea;
    }

    flux[nLaminar] = 0.0;
    if (nElectronIndex >= 0)
    {
        flux[nLaminar - 1] = 0.0;
    }

    //! Modify the fluxes of vibration and electron energy terms.
    for (int m = 0; m < nTemperatureModel - 1; ++ m)
    {
        flux[nLaminar + 1 + m] += spectrumRadiusTerm * deltaQ[nLaminar + 1 + m];
        flux[nLaminar + 1 + m] *= 0.5 * faceArea;
    }
    delete [] speciesBeta;    speciesBeta = nullptr;
}

void ChemicalMXDQR(RDouble *primitiveVars, RDouble *temperature, RDouble nx, RDouble ny, RDouble nz, RDouble faceArea, RDouble moveVelocity,
                   RDouble *deltaQ, RDouble *flux, int nNSEquation, int nLaminar, int nTemperatureModel, RDouble radius, int iSign, int nElectronIndex,
                   RDouble *speciesCvs, RDouble *speciesEtrs, RDouble totalH,  RDouble totalCv, RDouble totalCvtr, RDouble totalCvv, RDouble totalCve)
{
    using namespace IDX;
    using namespace GAS_SPACE;
    //! The number of species.
    int nSpecies = nLaminar - nNSEquation + 1;
    //! int nElectronIndex = GlobalDataBase::GetIntParaFromDB("nElectronIndex");
    RDouble alpha = 0.0, beta = 0.0, gama = 0.0, *speciesBeta;

    int mTT = gas->GetmTT();
    int mTE = gas->GetmTE();
    Thermo_Energy  *Thermo_Energy_temparay = gas->GetThermo_Energy_temparay();
    speciesBeta = Thermo_Energy_temparay->Rms;    //! be careful. by dms

    RDouble uVelocity = primitiveVars[IU];
    RDouble vVelocity = primitiveVars[IV];
    RDouble wVelocity = primitiveVars[IW];
    RDouble squareVelocity, absoluteVelocity, relativeVelocity;

    //! The absolute velocity.
    absoluteVelocity = nx * uVelocity + ny * vVelocity + nz * wVelocity;
    //! The relative velocity.
    relativeVelocity = absoluteVelocity - moveVelocity;
    //! Compute squared velocity.
    squareVelocity = uVelocity * uVelocity + vVelocity * vVelocity + wVelocity * wVelocity;

    //! Obtain the temperatures.
    RDouble transRotationTemperature = temperature[mTT];
    RDouble electronTemperature = temperature[mTE];

    //! Compute the total enthalpy.
    RDouble totalEnthalpy = totalH + 0.5 * squareVelocity;

    //! Compute the coefficient b1 and b2 (see formula A.5 in PHengLEI Theory manual).
    RDouble b1 = nx * deltaQ[IRU] + ny * deltaQ[IRV] + nz * deltaQ[IRW] - absoluteVelocity * deltaQ[IR];
    RDouble dh = uVelocity * deltaQ[IRU] + vVelocity * deltaQ[IRV] + wVelocity * deltaQ[IRW] - deltaQ[IRE];
    RDouble b2 = 0.0, b3 = 0.0, electronPressure = 0.0;

    if (nTemperatureModel == 1)
    {
        //! Compute the coefficients in the chemical reaction.
        //gas->ComputeCoefficientInMXDQ(primitiveVars, alpha, beta, speciesBeta);
        gas->ComputeCoefficientInMXDQR(primitiveVars, transRotationTemperature, squareVelocity, speciesCvs, totalCv, alpha, beta, speciesBeta);

        b2 = beta * deltaQ[IR] - alpha * dh;
        for (int s = 0; s < nSpecies - 1; ++ s)
        {
            b2 += speciesBeta[s] * deltaQ[nNSEquation + s];
        }
    }
    else    //! Multi-temperature model.
    {
        //! Compute the pressure of electron.
        electronPressure = gas->GetElectronPressure(primitiveVars, electronTemperature);
        RDouble *speciesGama = Thermo_Energy_temparay->Ees;    //! be careful. by dms


        gas->GetMultiTemperatureModelPartialDerivativesR(primitiveVars, transRotationTemperature, electronTemperature,alpha, beta, gama, speciesBeta, speciesGama,
            squareVelocity, speciesEtrs, totalCvtr, totalCvv, totalCve);

        if (nTemperatureModel == 2)
        {
            dh += deltaQ[nLaminar + 1];
            b2 += gama * deltaQ[nLaminar + 1];
            b3 += gama * deltaQ[nLaminar + 1];
        }
        else if (nTemperatureModel == 3)
        {
            b2 += gama * deltaQ[nLaminar + 1];
            b3 += gama * deltaQ[nLaminar + 2];
            dh += deltaQ[nLaminar + 1] + deltaQ[nLaminar + 2];
        }
        b2 += beta * deltaQ[IR] - alpha * dh;

        for (int iSpices = 0; iSpices < nSpecies - 1; ++ iSpices)
        {
            b2 += speciesBeta[iSpices] * deltaQ[nNSEquation + iSpices];
            b3 += speciesGama[iSpices] * deltaQ[nNSEquation + iSpices];
        }
    }

    //! Comoute the fluxes called M*dQ.
    flux[IR ] = relativeVelocity * deltaQ[IR ] + b1;
    flux[IRU] = relativeVelocity * deltaQ[IRU] + b1 * uVelocity     + b2 * nx;
    flux[IRV] = relativeVelocity * deltaQ[IRV] + b1 * vVelocity     + b2 * ny;
    flux[IRW] = relativeVelocity * deltaQ[IRW] + b1 * wVelocity     + b2 * nz;
    flux[IRE] = relativeVelocity * deltaQ[IRE] + b1 * totalEnthalpy + b2 * absoluteVelocity;
    for (int m = nNSEquation; m < nLaminar; ++ m)
    {
        flux[m] = relativeVelocity * deltaQ[m] + b1 * primitiveVars[m];
    }

    //! Modify the fluxes of multi-temperature model.
    if (nTemperatureModel == 2) 
    {
        flux[nLaminar + 1] = relativeVelocity * deltaQ[nLaminar + 1] + b1 * (primitiveVars[nLaminar + 1] + electronPressure / primitiveVars[IR]) + b3 * absoluteVelocity;
    }
    else if (nTemperatureModel == 3)
    {
        flux[nLaminar + 1] = relativeVelocity * deltaQ[nLaminar + 1] + b1 * primitiveVars[nLaminar + 1];
        flux[nLaminar + 2] = relativeVelocity * deltaQ[nLaminar + 2] + b1 * (primitiveVars[nLaminar + 2] + electronPressure / primitiveVars[IR]) + b3 * absoluteVelocity;
    }

    //! To add rM * dQ that is the product of the invscid spectrum radius rM and increment dQ.
    RDouble spectrumRadiusTerm = iSign * radius / faceArea;
    for (int m = 0; m < nLaminar; ++ m)
    {
        flux[m] += spectrumRadiusTerm * deltaQ[m];
        flux[m] *= 0.5 * faceArea;
    }

    flux[nLaminar] = 0.0;
    if (nElectronIndex >= 0)
    {
        flux[nLaminar - 1] = 0.0;
    }

    //! Modify the fluxes of vibration and electron energy terms.
    for (int m = 0; m < nTemperatureModel - 1; ++ m)
    {
        flux[nLaminar + 1 + m] += spectrumRadiusTerm * deltaQ[nLaminar + 1 + m];
        flux[nLaminar + 1 + m] *= 0.5 * faceArea;
    }
}

}
}
