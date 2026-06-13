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
//! @file      GridGroup.h
//! @brief     Explain this file briefly.
//! @author    xxx.

#pragma once
#include "Pointer.h"
#include "Geo_StructBC.h"
#include "TK_Exit.h"
#include "MultigridFactory.h"
#include "MixGrid.h"

using namespace std;

namespace PHSPACE
{
class Region;

class GridGroup
{
public:
    GridGroup(int zoneStart = 0);
    ~GridGroup();
private:
    int nzones;
    int *block_proc_dump;    //! Bell 20131124 add
    int * block_proc;
    int * block_proc_grid;
    int * block_type;
    int * block_idx ;
    int * block_fileProc;
    vector< Grid *> grids;
    int iBlock;
    int zoneStart;
    Region *region;
public:
    void AddGrid(Grid *grid);

public:
    void SetIBlock(int iBlock) { this->iBlock = iBlock; }
    bool IsZoneLayoutInitialized();
    void SetZoneHandler(Region *region) { this->region = region; }
    Region * GetRegion() { return region; }
    void SetNZones(int nZonesIn) { this->nzones = nZonesIn; }
    int GetNZones() const { return nzones; }
    int GetNumberofGrid() const;
    Grid * GetGrid(int iZone) const { return grids[iZone]; }

    //void ReadSingleGrid(fstream &file, int iZone);
    void InitZoneLayoutInformation(const string &filename);

    void SetBlockProcDump(int * block_proc_dump) { this->block_proc_dump = block_proc_dump; }
    void SetBlockProc(int * block_proc) { this->block_proc = block_proc; }
    void SetBlockProcGrid(int * block_proc) { this->block_proc_grid = block_proc; }
    void SetBlockType(int * block_type) { this->block_type = block_type; }
    void SetBlockIndex(int * block_idx) { this->block_idx = block_idx; }
    void SetBockFileProc(int * block_fileProc) { this->block_fileProc = block_fileProc; }

    int * GetBlockProcDump() { return this->block_proc_dump; }
    int * GetBlockProc()     { return this->block_proc; }
    int * GetBlockProcGrid() { return this->block_proc_grid; }
    int * GetBlockType()     { return this->block_type; }
    int * GetBlockIndex()    { return this->block_idx; }
    int * GetBlockFileProc() { return this->block_fileProc; }
    int   GetZoneStart()     { return this->zoneStart; }

    void SetGlobalZoneLayout();
private:
    void InitZoneLayout(fstream &file);
    
};
}