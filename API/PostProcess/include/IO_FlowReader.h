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
//! @file      IO_FlowReader.h
//! @brief     Flow file raed.
//! @author    Baka.

#pragma once
#include "Region.h"
#include "IO_HDF5File.h"

using namespace std;

namespace PHSPACE
{

class IO_FlowReader
{
public:
    Region *region;

public:
    IO_FlowReader(Region *region_in);
    ~IO_FlowReader();

    void Run();

private:
    void Init();
    void InitMemoryStruct(Grid *gridIn);
    void InitMemoryUnstruct(Grid *gridIn);

    void ReadFlowFile();
    void ReadFlowFileStr(hid_t loc_id, Grid *gridIn);
    void ReadFlowFileUnStr(hid_t loc_id, Grid *gridIn);

    void ReadTurbFile();
    void ReadTurbFileStr(hid_t loc_id, Grid *gridIn);
    void ReadTurbFileUnStr(hid_t loc_id, Grid *gridIn);

    void ReadTransitionFile();
    void ReadTransitionFileStr(hid_t loc_id, Grid *gridIn);
    void ReadTransitionFileUnStr(hid_t loc_id, Grid *gridIn);

private:
    string restartNSFile;
    string turbFile;
    string transitionFile;

    int nEquation;
    int nTurbulenceEquation;
    int nTemperatureModel;
    int nTransitionEquation;
};

}