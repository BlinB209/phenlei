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
//! @file      Pre_GridConversion.h
//! @brief     Grid conversion from other format to PHengLEI.
//! @author    Bell.

#pragma once
#include "GridType.h"
using namespace std;

namespace PHSPACE
{
//! @brief Pre_GridConversion class defines the method of 
//! Grid conversion from other format to PHengLEI.
class Pre_GridConversion
{
protected:
    //! Grid file name of original grid that need to be converted.
    string gridFileName;

    //! The number of blocks in the grid file.
    int nBlocks;

    //! The final converted computational grids.
    Grid **grids;

public:
    //! @param[in] Grid file name of original grid that need to be converted.
    LIB_EXPORT Pre_GridConversion(const string &gridFileName);

    LIB_EXPORT virtual ~Pre_GridConversion();

public:
    //! Run to convert field view grid.
    LIB_EXPORT void Run();

    //! Return the converted computational grids.
    Grid ** GetGrids();

    //! Return the number of grid blocks.
    int GetNumberofBlocks();

    //! Assign the number of grid blocks.
    void SetNumberofBlocks(int dataIn);

    //! Clear the memory after run, which is the converted computational grids.
    //! Warning: if this function is executed, the converted computational grids would be deleted.
    LIB_EXPORT void Clear();

    void WriteVolumeName(const string &out_grid_file, Grid **grids_in);

private:
    virtual void ReadGrid() = 0;
    virtual void Conversion() = 0;

protected:
    void WriteFaceBC(const string &out_grid_file, Grid **grids_in);
    void WriteFaceBoundaryName(const string &out_grid_file, Grid **grids_in);
    void WriteCellToNode(const string &targetGridFileName, Grid **grids_in);
};

//! Dump the converted grid to file.
LIB_EXPORT void DumpGrid(const string &gridFileName, Grid **grids, int nBlocks);

LIB_EXPORT void DumpOldGrid(const string &gridFileName, Grid **grids, int nBlocks);
//!
void GetGridsBoundingBox(Grid **grids, int nBlocks, RDouble *pmin, RDouble *pmax);

//! Get the minimum and maximum edge length of grid.
void GetGridsMinMaxDS(Grid **grids, int nBlocks, RDouble &mindis, RDouble &maxdis);

//!
void MatchInterface(Grid *grid, LinkStruct *link);

//!
void ConstructGlobalInterfaceLink(Grid **grids,int nBlocks);

//!
LIB_EXPORT void DumpAdditionalData(const string &gridFileName, Grid **grids, int nBlocks);
//!
LIB_EXPORT void DumpFaceBC(const string &gridFileName, Grid **grids, int nBlocks);

//! Dump Characteristic Boundary to file.
void DumpCharacteristicBoundary(const string &gridFileName, Grid **grids, int nBlocks);
void WriteVolumeInfo(Grid **grids_in, int nBlocks);
void WriteBoundaryInfo(Grid **grids_in, int nBlocks);
void DumpBoundaryInfo(const set< pair<int, string> > &bcNamePairSet);

void WriteBCNameInfo(const string &gridFileName, Grid **grids_in, int nBlocks);
void DumpBCNameInfo(const set< pair<int, string> > &bcNamePairSet, const string &gridFileName);

void DumpVCNameFile(const set< pair<int, string> > &vcNamePairSet, const string &gridFileName);
void DumpVolumeConditionFile(const set< pair<int, string> > &vcNamePairSet);

void ComputeCharacteristicBoundary(Grid *grid_in, DataContainer *data);
void ComputeCharacteristicBoundaryUns(Grid *gridIn, DataContainer *data);
void ComputeCharacteristicBoundaryStr(Grid *grid_in, DataContainer *data);

#include "Pre_GridConversion.hxx"

}