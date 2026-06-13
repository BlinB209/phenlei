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
//! @file      Mesh_DeformationSPRING.h
//! @brief     Grid deformation by SPRING method.
//! @author    Bell, Baka.

//! REFERENCES
//! [1] Zhao Z, et al. Numerical simulation of unsteady flows on bird-like flapping 
//!     wing[J]. Transactions of Nanjing University of Aeronautics & Astronautics, 
//!     2013, 30(sup):93-100.

#pragma once
#include "Mesh_Deformation.h"

namespace PHSPACE
{

//! @brief Mesh_DeformationSPRING class achieve unstructured grid deform by SPRING method.
class Mesh_DeformationSPRING : public Mesh_Deformation
{
public:
    //! @param[in] nZonesIn         Number of grid.
    //! @param[in] stationalGridIn  Original grid.
    Mesh_DeformationSPRING(int nZonesIn, Grid **stationalGridIn);
    ~Mesh_DeformationSPRING();

private:
    void Deforming();
    void SurfaceGridMove(int iStep);
    void SurfaceGridMoveByImput();
    void TransformGrid();
    void MatchControlPoints();
    void PostTreat();
};

}
