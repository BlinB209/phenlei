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
//! @file      PostProcess.h
//! @brief     Flow result post-processing.
//! @author    Baka.

#pragma once
#include "IO_FlowReader.h"
#include "PostProcessData.h"

using namespace std;

namespace PHSPACE
{

class PostProcess
{
public:
    Region *region;

public:
    PostProcess(Region *region_in);
    ~PostProcess();

public:
    void Run();

private:
    void ReadGridFile();
    void ReadFlowFile();
    void ComputeMetrics();

    void InitGlobalValues();
    void ComputeCoefficientOfStateEquation();
    void ComputeReferenceTsuth();
    void ComputeReferenceLengthDimensional();
    void ReadBCFile();
    void CompAirForceCoef();

private:

};

const int COMPAIRFORCE = 0;





}