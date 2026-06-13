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
//! @file      Zone.h
//! @brief     Explain this file briefly.
//! @author    He Xin, Bell.

#pragma once
#include <vector>
#include "Data_Field.h"
#include "Data_Param.h"
#include <sstream>
using namespace std;

namespace PHSPACE
{
class Grid;
class UnstructGrid;
class PHSolver;
class PHGeometry;
class Data_Param;
class Data_Field;
class ZoneNeighbor;
class Region;

class Zone
{
private:
    //! Zone Index.
    int index;

    //! Computational grid.
    PHGeometry *geom;

    //! Solvers list, load on the current zone.
    vector <PHSolver *> *solvers;

    //! Control parameters, usually read from file.
    Data_Param *zPara;

    //! Zonal data pointers.
    Data_Field *zField;

    //! Number of solvers on this zone.
    int nsolver;

public:
    Zone(int index);
    ~Zone();

public:
    int GetIndex(void) { return index; }
    int GetNumberOfSolvers() { return nsolver; }
    void SetNSolver(int n) { nsolver = n; }
    void ComputeWeight();
    void CoarseGrids();
    void ShowCoarseGrids();

    PHGeometry * GetGeometry() { return geom; }
    void SetGeometry(PHGeometry *In) { this->geom = In; }

    void AddSolver(PHSolver *solver);
    void AddGridAndCopyZonePara(Grid *grid);

    void AddGrid(Grid *grid);
    void AddGridToGlobal(int iZone);
    void UpdateAllData();
    void UpdateData(const string &name, void *data, int type, int size)
    {
        zPara->UpdateData(name, data, type, size);
    }

    void GetData(const string &name, void *data, int type, int size)
    {
        zPara->GetData(name, data, type, size);
    }

    void UpdateDataPtr(const string &name, void *data)
    {
        zField->UpdateDataPtr(name, data);
    }

    void * GetDataPtr(const string &name) const
    {
        return zField->GetDataPtr(name);
    }

    void SetInterface(Grid *base_grid, Grid *cgrid);
};

extern vector <Zone *> *GlobalZones;

void InitGlobalValuesOfAllZones();

void InitGlobalValuesGrid();
void InitGlobalValuesSolver();

//! To release the allocated memories for the object of chemical class.
void CleanGlobalValuesSolver();

}
