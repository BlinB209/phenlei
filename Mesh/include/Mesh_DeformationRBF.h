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
//! @file      Mesh_DeformationRBF.h
//! @brief     Grid deformation by RBF method.
//! @author    Bell, Baka.

//! REFERENCES
//! [1] Zhao Z, et al. Numerical simulation of unsteady flows on bird-like flapping
//!     wing[J]. Transactions of Nanjing University of Aeronautics & Astronautics, 
//!     2013, 30(sup):93-100.
//! [2] Zhao Z et al. An efficient large-scale mesh deformation method based on 
//!     MPI/OpenMP hybrid parallel radial basis function interpolation[J]. 
//!     Chin J Aeronaut, 2020, 33(5): 1392-1404.

#pragma once
#include "Mesh_Deformation.h"

namespace PHSPACE
{

typedef PHVector1D< Point3D > PointArray;
typedef DataStruct_AdtTree<int, RDouble> CPTree;
typedef DataStruct_AdtNode<int, RDouble> CPNode;

//! @brief Mesh_DeformationRBF class achieve unstructured grid deform by RBF method.
class Mesh_DeformationRBF : public Mesh_Deformation
{
private:
    //! Number of control point.
    int nControlPoints;

    //! The coordinate of each control point.
    PointArray controlPointCoordinates;

    //! The most displacement of each control point.
    PointArray controlPointVectors;

    //! The zone index and node index each control point belong
    PHVectorInt2D controlPointToGlobalPointIndex;

    //! The RBF influence radius parameter: multiply it to the maximum boundary deformation distance
    //! is the influence radius used in RBF method (influenceRadius in following).
    //! the points that out of the radius will not been deformed.
    //! 10.0-50.0 is recommended.
    RDouble influencePara;

    //! The ADTree of control point
    CPTree *controlPointTree;

    //! Number of reference Control Points.
    int numberOfReferenceCP;

    RDouble wallTolerance;
    PHVectorInt1D secondSegment;
    PHVectorInt1D rotateNodeIndex;

public:
    //! @param[in] nZonesIn         Number of grid.
    //! @param[in] stationalGridIn  Original grid.
    Mesh_DeformationRBF(int nZonesIn, Grid **stationalGridIn);
    ~Mesh_DeformationRBF();

    void SetNumberOfReferenceCP(int numberOfReferenceCPIn);
    void SetInfluencePara(RDouble influenceParaIn);

private:
    void TransformGrid();
    void MatchControlPoints();
    void PostTreat();
    void Deforming();
    void SurfaceGridMove(int iStep);
    void BuildControlPointTree();
    void FindOutReferenceCP(DYNode point, PHVectorInt1D &referenceCP);
    void RBFDeformingKernel(DynamicGrid *grid, int nodeIndex, PHVectorInt1D &referenceCP);
};

RDouble Wendland  (RDouble distance, RDouble radius);
RDouble WendlandC0(RDouble distance, RDouble radius);
RDouble WendlandC2(RDouble distance, RDouble radius);
RDouble WendlandC4(RDouble distance, RDouble radius);
//void MatrixInverse(vector< vector< RDouble > > &maxtrix, const int size);

}
