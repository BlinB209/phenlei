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
//! @file      Pre_GridTranslator.h
//! @brief     Explain this file briefly.
//! @author    xxx

#pragma once
#include "Region.h"

using namespace std;

namespace PHSPACE
{

class GridTranslator
{
private:
    int numberOfZones;
    int isOverset;
    vector<StructGrid *> *gridGroup;
    vector<StructFace *> *structFaceGroup;
public:
     GridTranslator();
     GridTranslator(Grid **&strgrids, int &nBlocks);
    ~GridTranslator();
public:
    void Run();
private:
    void ReadGridgenFiles();
    void ReadGridgenCoordinateFile();
    void ReadGridgenBoundaryFile();
    void ProcessBoundaryCondition();
    void GenerateStructFaceGroup();
    void ColorateStructFaceGroup();
    void CollectInterfaceInformation();
    void OutputMultiZoneComputationalGrids();
    void OutputGrid();
    void CheckGrid();
};

}