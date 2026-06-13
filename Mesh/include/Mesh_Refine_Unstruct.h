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
//! @file      Mesh_Refine_Unstruct.h
//! @brief     Parallel unstructured grid refine.
//! @author    Baka, Bell.

//! REFERENCES
//! [1] Zhao Z et al. A large-scale parallel hybrid grid generation technique for 
//!     realistic complex geometry[J]. Int J Numer Meth Fluids, 2020, 92(10): 1235-1255

#pragma once
#include <stdio.h>
#include <string>
#include "Region.h"
#include "Mesh_Refine.h"
#include "Mesh_RefineController.h"

namespace PHSPACE
{

//! @brief Mesh_Refine_Unstruct achieve unstructured grid refine.
class Mesh_Refine_Unstruct : public Mesh_Refine
{
public:
    Mesh_Refine_Unstruct(const string &gridFileNameIn);
    ~Mesh_Refine_Unstruct();

private:
    //! Initialize: allocate memory.
    void AllocateMemory();

    //! Construct grid topology information.
    void ConstructGridTopo();
    void ConstructGridTopo(Grid *grid);
    void ConstructCellTopo(Grid *grid);
    void ConstructNodeTopo(Grid *grid);
    void ConstructFaceTopo(Grid *grid);
    void ConstructMetrics (Grid *grid);
    void ConstructInterFaceTopo(Grid *grid);

    //! Set Refine property of each cell.
    void BuildRefineProperty();
    void SetDefaultRefineProperty(int iZone);
    void ComputeRefineProperty(int iZone);
    void PrepareDataStructure(int iZone);
    void SetCellModifiedStatusToNoChange(int iZone);
    void SetCellComputationalStatusToHiddenState(int iZone);
    void SetCurrentCellModifiedStatusByOldCellModifiedStatus(int iZone);
    void ResetCellProperty(int iZone);

    //! Set current cell computational status to ON, if old cell modified status is NOCHANGE;
    void SetCurrentCellComputationalStatus(int iZone);

    //! Set refine type of each cell.
    void FixAnisoRefineType();
    void InitAnisoRefineType();
    void SwapInterfaceRefineType();
    bool IsNeedInfectInterfaceIsotropicRefineType();

    //! Grid refine.
    void RefineGrid();

    //! Generate refined grid, and then dump grid.
    void GenerateAndDumpComputationalGrid();
    void AllocateComputationalGrid();
    void GenerateComputationalGrid();
    void GenerateGlobalInterfaceLink();

private:
    Mesh_RefineController **gridRefineControllers;

};

}
