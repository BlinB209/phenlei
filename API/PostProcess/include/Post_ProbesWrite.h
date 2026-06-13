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
//! @file      Post_ProbesWrite.h
//! @brief     Write the monitored probes flow field variables into file.
//! @author    Meng Liyuan.

#pragma once
#include "LIB_Macro.h"
#include "TypeDefine.h"

namespace PHSPACE
{
class ActionKey;
class Grid;
class DataContainer;
class Post_Probes;

class Post_ProbesWrite
{
public:
    LIB_EXPORT Post_ProbesWrite(int flowType = 0);
    LIB_EXPORT ~Post_ProbesWrite();

public:
    //LIB_EXPORT void Run(bool dumpAllLevel = false);
    LIB_EXPORT void Run();

private:
    vector<DataContainer *> datalist;
    int flowType;

private:
    //! Communicate probes variables with each processor.
    void ServerProbesWrite(ActionKey *actkey);

    //! Store the probes variables.
    void StoreProbesVarariables(ActionKey *actkey, Post_Probes *postProbesVar);

    //! Monitored probes post-processing.
    void PostProbesWrite(ActionKey *actkey);

    void Clear();
};

}