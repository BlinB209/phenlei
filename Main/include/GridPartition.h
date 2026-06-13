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
//! @file      GridPartition.h
//! @brief     Grid partition: Structured and Unstructured.
//! @author    He Xin, Bell, Baka.

#pragma once
#include <fstream>
#include <vector>
#include "Math_BasisFunction.h"
using namespace std;

namespace PHSPACE
{
class Grid;
class StructGrid;
class InterfaceInfo;
class DataContainer;
class GridRefine;

void PartitionGrid();
void PartitionUnstructGrid();
void RefineStructGrid();
void PartitionStructGrid();

void MultiGridPartition();
void GetPartitionPara(PHVectorString1D &gridFileNameList, PHVectorString1D &outFileNameList, PHVectorInt1D &gridTypeList, PHVectorInt1D &maxProcList);

void RefineStructGrid(const string &coordinateFile, const string &layoutFile);

void ReadGridgenBC(const string &BCFile, DataContainer *cData);
void ReadGridgenBC(fstream &file, DataContainer *cData);

class BlockBCPatchRefine
{
public:
    BlockBCPatchRefine() {};
    ~BlockBCPatchRefine() {};

public:
    void SetContent(const int &nb, const int &imin, const int &imax, const int &jmin, const int &jmax, const int &kmin, const int &kmax);
    void RefineBC();

private:
    int st[3], ed[3];
    int nd[3];
    int nb;
};

class BlockRefine;

class BlockBCRefine
{
public:
    BlockBCRefine() {next = 0;};
    ~BlockBCRefine() {delete next;};

public:
    void SetContent(const int &BCType, const int &imin, const int &imax, const int &jmin, const int &jmax, const int &kmin, const int &kmax);
    void SetPatchContent(const int &nb, const int &imin, const int &imax, const int &jmin, const int &jmax, const int &kmin, const int &kmax);
    void Refine();
    void RefineBC();

private:
    int st[3],ed[3];
    int nd[3];

    BlockBCPatchRefine *next;
    int BCType;
    int nb;
};

class BlockRefine
{
public:
    BlockRefine() {};
    ~BlockRefine() {};

public:
    int size() {return dim[0] * dim[1] * dim[2];}
    int cellsize() {return (dim[0] - 1) * (dim[1] - 1) * MAX(dim[2] - 1, 1);}
    void SetID(const int &id) {this->id = id;}
    int GetID() const {return id;}

    void SetNI(const int &ni) {this->dim[0] = ni;}
    void SetNJ(const int &nj) {this->dim[1] = nj;}
    void SetNK(const int &nk) {this->dim[2] = nk;}

    int  GetNI() const {return dim[0];}
    int  GetNJ() const {return dim[1];}
    int  GetNK() const {return dim[2];}

    vector <BlockBCRefine *> &GetBCList() {return BCList;}

public:
    void Refine();
    void SetNewDimension();
    void DumpGrid(fstream &file, Grid *gridIn);
    void RefineCoordinate(RDouble3D &x, RDouble3D &xold, const int &ni, const int &nj, const int &nk);
    StructGrid *CreateStructGrid(const int &ni, const int &nj, const int &nk);

private:
    int dim[3], oldDim[3];
    int id;
    vector <BlockBCRefine *> BCList;
};

class GridRefine
{
private:
    vector <BlockRefine *> refblock;
    int numberOfTotalBlocks;
    int numberOfBlocks;
public:
    GridRefine() {};
    ~GridRefine() {};
public:
    void Init(DataContainer *cData);
    void Refine();
    void DumpGrid();
    void DumpLayout(DataContainer *cData);
    void RefineStructGrid();
private:
    void SetNewDimension();
    int  MapBC(int i);
};

void GetIJKRegion(int *dim, int &ist, int &ied, int &jst, int &jed, int &kst, int &ked, int m);
void GetIJKRegion(int *st, int *ed, int &ist, int &ied, int &jst, int &jed, int &kst, int &ked, int m);
}
