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
//! @file      CFDSimulation.h
//! @brief     The entrance of CFD simulation.
//! @author    He Xin.

#pragma once

namespace PHSPACE
{
//! Entrance of PHengLEI.
void RunPHengLEI();

//! Show welcome information.
void ShowWelcomeTitle();

//! After reading control parameters from files, a preliminary self-check will be carried out to avoid the easy-neglected mistakes.
static void ControlParametersPreliminarySelfCheck();

//! According to the functions of PHengLEI, the control parameters self-check process includes six parts listed as follows.
static void ParametersSelfCheckForCFD();
}
