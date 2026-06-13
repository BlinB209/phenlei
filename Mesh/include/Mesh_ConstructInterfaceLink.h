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
//! @file      Mesh_ConstructInterfaceLink.h
//! @brief     Construct grid interface information, support for parallel operation. 
//! @author    Baka.

#pragma once
#include "Geo_Grid.h"
#include "ComputationalGrid.h"

namespace PHSPACE
{

class Mesh_ConstructInterfaceLink
{
private:
    //! Number of grid in this process.
    int nZonesInCurrentProcess;

    //! Grids in this process.
    Grid **gridContainer;

    //! The number of global zones.
    int nZonesGlobal;

    //! 
    int *zoneProc;

    int *localZoneIndex;

    //! Interface information
    vector< ZoneInterface * > zoneInterfacesTopologyList;

public:
    Mesh_ConstructInterfaceLink(Grid **gridIn, int nZonesInCurrentProcessIn);
    ~Mesh_ConstructInterfaceLink();

public:
    void Run();

private:
    void ReorderGridIndex();
    void ChangeNoBoundaryCondition();

    void BuildZoneInterfacesInfo();
    void CollectUnsGridInterfaceInfo(Grid *grid, DataContainer *cdata);
    void DecodeInterfaceInfo(ZoneInterface *zoneInterface, DataContainer *cdata);
    void ConstructInterfaceLink();
    void BuildGridLink(int iZone, LinkStruct *link);

    void GetBoundingBox(RDouble *pmin, RDouble *pmax);
    void GetMinMaxDS(RDouble &mindis, RDouble &maxdis);


};




}