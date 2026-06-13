#include <cmath>
#include "Gas.h"
#include "GlobalDataBase.h"
#include "Constants.h"
#include "TK_Exit.h"


namespace PHSPACE
{
namespace GAS_SPACE
{
    Gas *gas;

    int GetIndexBySpeciesName(string speciesName)
    {
        int nIndex = -1;
        if (speciesName == "O")
        {
            nIndex = 0;
        }
        else if (speciesName == "O2")
        {
            nIndex = 1;
        }
        else if (speciesName == "NO")
        {
            nIndex = 2;
        }
        else if (speciesName == "N")
        {
            nIndex = 3;
        }
        else if (speciesName == "N2")
        {
            nIndex = 4;
        }
        else if (speciesName == "O+")
        {
            nIndex = 5;
        }
        else if (speciesName == "O2+")
        {
            nIndex = 6;
        }
        else if (speciesName == "NO+")
        {
            nIndex = 7;
        }
        else if (speciesName == "N+")
        {
            nIndex = 8;
        }
        else if (speciesName == "N2+")
        {
            nIndex = 9;
        }
        else if (speciesName == "e-")
        {
            nIndex = 10;
        }
        else if (speciesName == "C")
        {
            nIndex = 11;
        }
        else if (speciesName == "C2")
        {
            nIndex = 12;
        }
        else if (speciesName == "CO")
        {
            nIndex = 13;
        }
        else if (speciesName == "CO2")
        {
            nIndex = 14;
        }
        else if (speciesName == "CN")
        {
            nIndex = 15;
        }
        else
        {
            nIndex = -1;
        }
        return nIndex;
    }

    string GetSpeciesNameByIndex(int nIndex)
    {
        string speciesName = "";
        if (nIndex == 0)
        {
            speciesName = "O";
        }
        else if (nIndex == 1)
        {
            speciesName = "O2";
        }
        else if (nIndex == 2)
        {
            speciesName = "NO";
        }
        else if (nIndex == 3)
        {
            speciesName = "N";
        }
        else if (nIndex == 4)
        {
            speciesName = "N2";
        }
        else if (nIndex == 5)
        {
            speciesName = "O+";
        }
        else if (nIndex == 6)
        {
            speciesName = "O2+";
        }
        else if (nIndex == 7)
        {
            speciesName = "NO+";
        }
        else if (nIndex == 8)
        {
            speciesName = "N+";
        }
        else if (nIndex == 9)
        {
            speciesName = "N2+";
        }
        else if (nIndex == 10)
        {
            speciesName = "e-";
        }
        else if (nIndex == 11)
        {
            speciesName = "C";
        }
        else if (nIndex == 12)
        {
            speciesName = "C2";
        }
        else if (nIndex == 13)
        {
            speciesName = "CO";
        }
        else if (nIndex == 14)
        {
            speciesName = "CO2";
        }
        else if (nIndex == 15)
        {
            speciesName = "CN";
        }
        else
        {
            speciesName = "";
        }
        return speciesName;
    }

    void SpeciesNameToInteger(int numberOfSpecies, string *nameList, Int1D *speciesOrder)
    {
        //Initialization.
        for (int n = 0; n < numberOfSpecies; ++ n)
        {
            (*speciesOrder)(n) = -1;
        }
        for (int n = 0; n < numberOfSpecies; ++ n)
        {
            int nIndex = GetIndexBySpeciesName(nameList[n]);
            if (nIndex >= 0)
            {
                (*speciesOrder)(nIndex) = n;
            }
        }
    }
    void SpeciesNameToInteger(int numberOfSpecies, string *nameList, int *speciesOrder)
    {
        //Initialization.
        for (int n = 0; n < numberOfSpecies; ++ n)
        {
            speciesOrder[n] = -1;
        }
        for (int n = 0; n < numberOfSpecies; ++ n)
        {
            int nIndex = GetIndexBySpeciesName(nameList[n]);
            if (nIndex >= 0)
            {
                speciesOrder[nIndex] = n;
            }
        }
    }

    bool IsEqualGasModel(Int1D *speciesOrder1, Int1D *speciesOrder2)
    {
        //the index of the array are O, O2, NO, N, N2, O+, O2+, NO+, N+, N2+, e-, C, C2, CO, CO2, CN, respectively.
        //the elements indicate the order of the species in the gas model, -1 indicates the species does not exist.
        int totalNumber = 16;
        for (int n = 0; n < totalNumber; ++ n)
        {
            if ((*speciesOrder1)(n) != (*speciesOrder2)(n))
            {
                return false;
            }
            //if (((*speciesOrder1)(n) >= 0 && (*speciesOrder2)(n) < 0) ||
            //    ((*speciesOrder1)(n) < 0 && (*speciesOrder2)(n) >= 0))
            //{
            //    return false;
            //}
        }
        return true;
    }

    Gas::Gas()
    {
        wallMultiTemperature = 0;
        if (GlobalDataBase::IsExist("wallMultiTemperature", PHINT, 1))
        {
            wallMultiTemperature = GlobalDataBase::GetIntParaFromDB("wallMultiTemperature");
        }

        nDiagonalModified = 0;
       
        nDiagonalModifiedTurb = 0;
        
        trTemperatureMinNonDim = 10.0 ;
        if (GlobalDataBase::IsExist("trTemperatureMin", PHDOUBLE, 1))
        {
            trTemperatureMinNonDim = GlobalDataBase::GetDoubleParaFromDB("trTemperatureMin") ;
        }
        trTemperatureMinNonDim /= GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    }



Thermo_Energy_Data::Thermo_Energy_Data(int numberOfSpecies)
{
    Cps = new RDouble[numberOfSpecies];
    Cvs = new RDouble[numberOfSpecies];
    Cvtrs = new RDouble[numberOfSpecies];
    Cvvs = new RDouble[numberOfSpecies];
    Cves = new RDouble[numberOfSpecies];

    Hs = new RDouble[numberOfSpecies];
    Es = new RDouble[numberOfSpecies];
    Etrs = new RDouble[numberOfSpecies];
    Evs = new RDouble[numberOfSpecies];
    Ees = new RDouble[numberOfSpecies];
    Rms = new RDouble[numberOfSpecies];
}

Thermo_Energy_Data::Thermo_Energy_Data()
{
    Cps = NULL;
    Cvs = NULL;
    Cvtrs = NULL;
    Cvvs = NULL;
    Cves = NULL;

    Hs = NULL;
    Es = NULL;
    Etrs = NULL;
    Evs = NULL;
    Ees = NULL;
    Rms = NULL;
}

Thermo_Energy_Data::~Thermo_Energy_Data()
{
    if (Cps != NULL) delete[]Cps;
    if (Cvs != NULL) delete[]Cvs;
    if (Cvtrs != NULL) delete[]Cvtrs;
    if (Cvvs != NULL) delete[]Cvvs;
    if (Cves != NULL) delete[]Cves;

    if (Hs != NULL) delete[]Hs;
    if (Es != NULL) delete[]Es;
    if (Etrs != NULL) delete[]Etrs;
    if (Evs != NULL) delete[]Evs;
    if (Ees != NULL) delete[]Ees;
    if (Rms != NULL) delete[]Rms;
}

void Thermo_Energy_Data::Init(int numberOfSpecies, int ntmodel)
{
    Cps = new RDouble[numberOfSpecies];
    Cvs = new RDouble[numberOfSpecies];
    Cvtrs = new RDouble[numberOfSpecies];
    Cvvs = new RDouble[numberOfSpecies];
    Cves = new RDouble[numberOfSpecies];

    Hs = new RDouble[numberOfSpecies];
    Es = new RDouble[numberOfSpecies];
    Etrs = new RDouble[numberOfSpecies];
    Evs = new RDouble[numberOfSpecies];
    Ees = new RDouble[numberOfSpecies];
    Rms = new RDouble[numberOfSpecies];
}

void Thermo_Energy_Data::freeData()
{
    if (Cps != NULL) delete[]Cps;
    if (Cvs != NULL) delete[]Cvs;
    if (Cvtrs != NULL) delete[]Cvtrs;
    if (Cvvs != NULL) delete[]Cvvs;
    if (Cves != NULL) delete[]Cves;

    if (Hs != NULL) delete[]Hs;
    if (Es != NULL) delete[]Es;
    if (Etrs != NULL) delete[]Etrs;
    if (Evs != NULL) delete[]Evs;
    if (Ees != NULL) delete[]Ees;
    if (Rms != NULL) delete[]Rms;
}


Thermo_Energy_DB ::Thermo_Energy_DB(int imin, int imax ,int jmin, int jmax, int kmin, int kmax, int numberSpecies, int ntmodel)
{
    ni = imax - imin + 1;
    nj = jmax - jmin + 1;
    nk = kmax - kmin + 1;

    istart = imin;
    jstart = jmin;
    kstart = kmin;

    nmax = ni * nj * nk;
    numberOfSpecies = numberSpecies;

    Data = new Thermo_Energy_Data[nmax];
    for (int m = 0; m < nmax; ++ m)
    {
        Data[m].Init(numberOfSpecies, ntmodel);
    }
}

Thermo_Energy_DB ::~Thermo_Energy_DB()
{
    if (Data != NULL)
    {
        if(numberOfSpecies > 0)
        {
            for (int m = 0; m < nmax; ++ m)
            {
                Data[m].freeData();
            }
        }
        delete[]Data;
    }
}

Thermo_Energy_Data * Thermo_Energy_DB ::GetThermo_Energy_Data(int i, int j ,int k)
{
    return  &Data[((k - kstart) * nj + (j - jstart)) * ni + i];
}

}

}