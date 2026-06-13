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
//! @file      SolversList.h
//! @brief     This file include all solvers' head file, if you want develop new solver, include here.
//! @author    Bell, He Xin.

#pragma once

#include "NSSolver.h"                 //! Base NS solver.
#include "NSSolverStructFD.h"         //! Structured NS solver, Finite Difference Method.
#include "NSSolverStructFV.h"         //! Structured NS solver, Finite Volume Method.
#include "NSSolverUnstruct.h"         //! Unstructured NS solver, Finite Volume Method.

#include "TurbSolver.h"               //! Base turbulent solver.
#include "TurbSolverStruct.h"         //! Structured turbulent solver.
#include "TurbSolverUnstr.h"          //! Unstructured turbulent solver.
#include "TurbSolverStructFD.h"       //! Structured FD turbulent solver.

#include "TransitionSolver.h"               //! Base transition solver.
#include "TransitionSolverStruct.h"         //! Structured transition solver.
#include "TransitionSolverUnstr.h"          //! Unstructured transition solver.

#ifdef USE_DEMOSOLVER
#include "DemoSolver.h"
#endif
