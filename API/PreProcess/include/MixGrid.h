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
//! @file      MixGrid.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include "Geo_Grid.h"
#include "Pre_GridBase.h"
#include "Geo_StructGrid.h"
#include "Geo_UnstructGrid.h"

namespace PHSPACE
{
class Grid;
class StructGrid;
class UnstructGrid;

const int MULTI_TYPE     = -1;
const int PHENGLEI_TYPE  = 1;
const int CGNS_TYPE      = 2;
const int PLOT3D_TYPE    = 3;
const int FIELDVIEW_TYPE = 4;
const int FLUENT_TYPE    = 5;
const int USTAR_TYPE     = 6;
const int MIXGRID_TYPE   = 7;
const int GMSH_TYPE      = 8;
const int GRIDGEN_TYPE   = 9;

typedef void (*CREATE_SINGLE_GRID)(Grid *comm_grid);

void SetBcRegion(StructGrid *grid, VInt &imin, VInt &imax, VInt &jmin, VInt &jmax, VInt &kmin, VInt &kmax, VInt &bctype);

void Fantasy2Ustar();
void Fantasy2Ustar2D();
void Fantasy2Ustar3D();

void ReProcessBCInfo(Grid **grids, int nBlocks);
void ReadLnkFile(const string &fileName_in, Grid **grids, int nBlocks);
void ReadLnkInfor(const string &fileName_in, Grid **grids, int nBlocks);
void FantasyWriteToGrd(const string &cel_file, Grid **grids, int nBlocks);
void FantasyWriteToInp(const string &bnd_file, Grid **grids, int nBlocks);

void FantasyWriteToFluent(const string &gridFileName, Grid **grids);

void ReadUstarGrid(const string &gridfile, UnstructGrid *grid);
void ReadUstarGrid3D_Binary(const string &gridfile, UnstructGrid *grid);

void StrGridToUnsGrid(StructGrid *grid, UnstructGrid *uns_grid);
void StrGridToUnsGrid(Grid **grids, int nBlocks, UnstructGrid *uns_grid);

void WriteLnkFile(const string &gridFileName, Grid **grids_in, int nBlocks);
void CheckMeshMultigrid(Grid **grids_in, int nBlocks);
void WriteBcFileFromFTS(const string &gridFileName, Grid **grids_in, int nBlocks);

void GetUnsIndex(StructGrid *grid, DataStruct_AdtTree<int, RDouble> *adtTree, int &count, RDouble tol, Int3D &unsindex, RDouble *xuns, RDouble *yuns, RDouble *zuns);

void FillSection(Grid **grids, int nBlocks, CGNSBase *base_cgns, BaseElement *base_elem, Int3D **unsindexlist);

class Ustar2Fantasy
{
private:
    string mapfilename;
    map<int, int> bcmap;
public:
    Ustar2Fantasy(string mapfile_in);
    ~Ustar2Fantasy();
    int TranslateBC(int bctype);
};

}