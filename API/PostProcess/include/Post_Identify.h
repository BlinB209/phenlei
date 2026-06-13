//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      Post_Identify.h
//! @brief     Identify the motive derivative
//! @author    He Kun.

#pragma once
#include "PHHeader.h"

using namespace std;

namespace PHSPACE
{

class Identify
{
public:
    bool isIdentify;
    int nForceVar;
    int nIdentifyVar;
    int nDataLine;
    int nwStep;
    size_t nLinestart, nLineEnd;
    int nLineOneCycle, nCycle;
    int nFrequency, nOrder;
    string identifyMethod;
    string motionType;
    RDouble reduceFrequency, amplitude;
    RDouble physicalTimeStep;
    RDouble lenthReference;
    vector <string> nameList;

    vector < vector <RDouble> > aeroForce;

    RDouble *c0, *c1, *c2;

public:
    Identify();
    ~Identify();

    void Run();
    void ReadForce();
    void ForceArrange();
    void IdentifyProcess();
    void GetNumberOfCycle();
    void getIntegral(int iVar);
    void getIntegrand(int iVar, RDouble *f);
    void GetIntegral4thOrder(int iVar, RDouble *f, RDouble &sum1, RDouble &sum2);
    void GetIntegral2ndOrder(int iVar, RDouble *f, RDouble &sum1, RDouble &sum2);

    void DumpHysteresisData();
    void DumpResults();
};

}