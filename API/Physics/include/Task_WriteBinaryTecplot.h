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
//! @file      Task_WriteBinaryTecplot.h
//! @brief     Write flow field data into tecio film.
//! @author    Xu Gang, Zhang Jian, He Xin.

#pragma once
#include "Geo_UnstructGrid.h"
#include "Task.h"
using namespace std;

namespace PHSPACE
{
//! @brief Task_WriteBinaryTecplot class defines the method about Writing flow field data into binary tecplot film.\n
//!
//! Task_WriteBinaryTecplot
class Task_WriteBinaryTecplot : public Task
{
public:
    LIB_EXPORT Task_WriteBinaryTecplot(int flowType = 0);
    LIB_EXPORT ~Task_WriteBinaryTecplot();
    LIB_EXPORT void ServerWrite(ActionKey *actkey);
    void PreTask(ActionKey *actkey);
    void MainTask(ActionKey *actkey);
    void PostTask(ActionKey *actkey);

private:
    vector < DataContainer * > datalist;
    int flowType;
};





}