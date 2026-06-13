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
//! @file      Geo_StructMerge.h
//! @brief     It is the class 'OriginalStructGridMerge', which is type of Structured geometry grid operation for grid merging.
//! @author    Xu Fengding, Xu Gang .

#pragma once
#include <iostream>
#include "Mesh_Deformation.h"
#include "Mesh_DeformationSPRING.h"
#include "Mesh_DeformationRBF.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "GeometryUnit.h"
#include "Geo_Grid.h"
#include "Geo_StructBC.h"
using namespace std;
namespace PHSPACE
{
class OriginalStructGridMerge
{
private:
    int                               nTotalblock;
    int                               nTotalZone;
    StructGrid                        **mergegrid;
    vector < vector < StructGrid* > > originalgridofBlock;
    int *ni, *nj, *nk;
    int                             *nTotalNode;
    int                             *nTotalFace;
    int                             *nTotalCell;
    RDouble                         ***coordinates;
    int                             *nBCRegion;
    vector < vector < StructBC* > > dumpOriginalBCRegion;
    StructBCSet                     **compositeBCRegion;
    int                             *nIFace;
    InterfaceInfo                   **interfaceInfo;
    int                             nDim;
public:
    OriginalStructGridMerge();

    ~OriginalStructGridMerge();

    //! Initial all information need to merge.
    void Initial(vector < vector < Grid* > > originalgridofBlockIn, int nTotalblock, int nTotalZone, int nDim);

    //! Compute the mergegrid.
    void Run();

    //! Return the mergegrid.
    void GetMergeGrid(StructGrid **mergeGridOut);

    //! Return nTotalblock.
    int GetnTotalblock() { return nTotalblock; };
};

//! The struct of Informationofboundary is used to store bctype¡¢bodyname and boundaryname of boundary.
struct Informationofboundary
{
    int    bcType;
    string bodyName;
    string boundaryName;
};

//! Certify the type of the two Informationofboundary or the two StructBC is the same.
int WhetherInformationofBoundarytheSame(Informationofboundary s1, Informationofboundary s2);

int WhetherInformationofBoundarytheSame(StructBC *s1, StructBC *s2);

//! WhethercanMergeI/J/KDirect is used to certify whether the boundaries whose face direction is i/j/k direction can be merged.
//! MergeBCI/J/K is used to merge all boundaries whose face direction is i/j/k direction.
//! the effective length of merged boundaries is stored in *nBC.
int WhethercanMergeIDirect(vector < StructBC* > bc, int nBC);

void MergeBCI(vector < StructBC* > &bc, int *nBC);

int WhethercanMergeJDirect(vector < StructBC* > bc, int nBC);

void MergeBCJ(vector < StructBC* > &bc, int *nBC);

int WhethercanMergeKDirect(vector < StructBC* > bc, int nBC);

void MergeBCK(vector < StructBC* > &bc, int *nBC);

//! Select all start/end faces of i/j/k direction ,divide them by type and merge them
void SelectDivideMerge(int iDirection, int iFacePosition,
    vector < StructBC* > &OriginalBCRegion, int &nOriginalBCRegion,
    vector < vector < StructBC* > > &bcRegionMerged, int &nBCType, vector < int > &numberofEveryType);

}