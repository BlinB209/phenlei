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
//! @file      Mesh_Refine.h
//! @brief     Grid refine.
//! @author    Bell, Baka.

#pragma once
#include "Region.h"

namespace PHSPACE
{

//! @brief RefineParameter class storage the parameter of grid refine.
class RefineParameter
{
private:
    //! Mesh refine type, ISOTROPICREFINE or ANISOTROPICREFINE.
    int anisoRefine;

public:
    void SetAnisoRefineType(int anisoRefineIn);
    int GetAnisoRefineType();

};

//! @brief Mesh_Refine achieve grid refine.
class Mesh_Refine
{
private:
    string gridFileName;

public:
    Mesh_Refine(const string &gridFileNameIn);
    ~Mesh_Refine();

public:
    void Run();
    void SetRefineParameter(RefineParameter *refineParameterIn);

protected:
    //! Read original grid.
    void ReadGrid();

    //! Initialize: allocate memory.
    virtual void AllocateMemory() {};

    //! Construct grid information.
    virtual void ConstructGridTopo() {};

    //! Set Refine property of each cell.
    virtual void BuildRefineProperty() {};

    //! Set refine type of each cell.
    virtual void FixAnisoRefineType() {};

    //! Grid refine.
    virtual void RefineGrid() {};

    //! Generate refined grid, and then dump grid.
    virtual void GenerateAndDumpComputationalGrid() {};
    void DumpComputationalGrid();

protected:
    //! Number of grid in this process.
    int numberOfZones;

    //! Grid refine parameter.
    RefineParameter *refineParameter;

    Region *region;
    Grid **refinedGrids;
};

}
