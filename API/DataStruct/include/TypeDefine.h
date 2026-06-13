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
//! @file      TypeDefine.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include "Precision.h"
#include "DataStruct_Array.h"
#include "DataStruct_AdtTree.h"
#include "DataStruct_KDTree.h"
#include "AMRDef.h"
#include <string>
#include <vector>
#include <map>
#include <numeric>
using namespace std;

namespace PHSPACE
{
    typedef PHArray<int, 1> Int1D;
    typedef PHArray<int, 2> Int2D;
    typedef PHArray<int, 3> Int3D;
    typedef PHArray<int, 4> Int4D;
    typedef PHArray<int, 5> Int5D;
    typedef PHArray<int, 6> Int6D;

    typedef PHArray<RFloat, 1> RFloat1D;
    typedef PHArray<RFloat, 2> RFloat2D;
    typedef PHArray<RFloat, 3> RFloat3D;
    typedef PHArray<RFloat, 4> RFloat4D;
    typedef PHArray<RFloat, 5> RFloat5D;
    typedef PHArray<RFloat, 6> RFloat6D;

    typedef PHArray<RDouble, 1> RDouble1D;
    typedef PHArray<RDouble, 2> RDouble2D;
    typedef PHArray<RDouble, 3> RDouble3D;
    typedef PHArray<RDouble, 4> RDouble4D;
    typedef PHArray<RDouble, 5> RDouble5D;
    typedef PHArray<RDouble, 6> RDouble6D;

    typedef vector < RDouble * > VDoubleP;

    typedef vector < int  > VInt;
    typedef vector < VInt > VVInt;

    typedef map<string, int    > IntDBase;
    typedef map<string, RDouble> RDoubleDBase;
    typedef map<string, RFloat > RFloatDBase;
    typedef map<string, void * > PointerDBase;

    typedef streamsize PHLongLong;

    typedef vector<int> PHInt1D;
    typedef vector<PHInt1D> PHInt2D;
    typedef vector<RDouble> PHDouble1D;
    typedef vector<PHDouble1D> PHDouble2D;
    typedef vector<PHDouble2D> PHDouble3D;
    typedef vector<string> PHString1D;

    typedef DataStruct_AdtNode<int, RDouble> PHProcessorNode;
    typedef DataStruct_AdtTree<int, RDouble> PHProcessorTree;

    typedef PHVector1D<int    > PHVectorInt1D;
    typedef PHVector1D<string > PHVectorString1D;
    typedef PHVector1D<RDouble> PHVectorRDouble1D;
    typedef PHVector2D<int    > PHVectorInt2D;

}

