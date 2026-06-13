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
//! @file      Geo_MultiGridInfo_Struct.h
//! @brief     It defines the multi-grid step of the structured grid,
//!            such as StepI, StepJ, StepK.
//! @author    Zhang Jian.

#pragma once
namespace PHSPACE
{

class Geo_MultiGridInfo_Struct
{
private:
    int istp, jstp, kstp;    //! Multigrid step.
public:
    Geo_MultiGridInfo_Struct();
    ~Geo_MultiGridInfo_Struct();
public:
    //!
    int GetMultigridStepI() const;

    //!
    int GetMultigridStepJ() const;

    //!
    int GetMultigridStepK() const;

    //!
    void SetMultigridStepI(int istp);

    //!
    void SetMultigridStepJ(int jstp);

    //!
    void SetMultigridStepK(int kstp);
};

#include "Geo_MultiGridInfo_struct.hxx"
}