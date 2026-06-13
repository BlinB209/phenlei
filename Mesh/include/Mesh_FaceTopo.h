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
//! @file      Mesh_FaceTopo.h
//! @brief     Store face information of grid, used in unstructured grid refining.
//! @author    Baka, Bell.

#pragma once
#include "TypeDefine.h"
#include "Pre_GridBase.h"

namespace PHSPACE
{
//! @brief Mesh_FaceTopo store face information of grid, used in unstructured grid refining.
class Mesh_FaceTopo
{
private:
    //! Face type, such as TRI_3, QUAD_4.
    PHVectorInt1D faceType;

    //! Face topo information of grid.
    PHVectorInt1D leftCellIndex;
    PHVectorInt1D rightCellIndex;
    PHVectorInt2D face2Node;

    //! Face to node information, after refined.
    PHVectorInt2D compositeFace2Node;

    //! Parent face index of each face.
    PHVectorInt1D parentFaceIndex;

    //! Children face index of each face.
    PHVectorInt2D childrenFaceIndex;

    //! Face index map between refined face and global face.
    PHVectorInt1D computationalToGlobalFaceIndex;
    PHVectorInt1D globalToComputationalFaceIndex;

    //! Face type face refined.
    PHVectorInt1D compositeFaceType;

    //! Face index map between children face and refined face.
    PHVectorInt1D childrenToCompositeFaceIndex;
    PHVectorInt1D compositeToChildrenFaceIndex;

    //! Middle node of each edge.
    PHVectorInt1D middleNodeOfEdge;
    set< MultiIndex > edgesWithMiddleNode;

    //! Face metrics information.
    PHVectorRDouble1D faceArea;
    PHVectorRDouble1D faceNormalX;
    PHVectorRDouble1D faceNormalY;
    PHVectorRDouble1D faceNormalZ;

    //! All face of grid, used to face search.
    set< MultiIndex > faceList;
    set< MultiIndex > compositeFacesList;

    //! Boundary condition of each grid.
    vector< SimpleBC * > faceBoundaryCondition;

public:
    Mesh_FaceTopo();
    ~Mesh_FaceTopo();

public:
    PHVectorInt1D & GetFaceType() { return faceType; }
    PHVectorInt1D & GetLeftCellIndex() { return leftCellIndex; }
    PHVectorInt1D & GetRightCellIndex() { return rightCellIndex; }
    PHVectorInt2D & GetFace2Node() { return face2Node; }
    PHVectorInt2D & GetCompositeFace2Node() { return compositeFace2Node; }
    PHVectorInt2D & GetIntermediateFace2Node() { return intermediateFace2Node; }

    PHVectorInt1D & GetParentFaceIndex() { return parentFaceIndex; }
    PHVectorInt2D & GetChildrenFaceIndex() { return childrenFaceIndex; }

    PHVectorInt1D & GetComputationalToGlobalFaceIndex() { return computationalToGlobalFaceIndex; }
    PHVectorInt1D & GetGlobalToComputationalFaceIndex() { return globalToComputationalFaceIndex; }

    PHVectorInt1D & GetCompositeFaceType() { return compositeFaceType; }
    PHVectorInt1D & GetChildrenToCompositeFaceIndex() { return childrenToCompositeFaceIndex; }
    PHVectorInt1D & GetCompositeToChildrenFaceIndex() { return compositeToChildrenFaceIndex; }
    PHVectorInt1D & GetMiddleNodeOfEdge() { return middleNodeOfEdge; }

    PHVectorRDouble1D & GetFaceArea() { return faceArea; }
    PHVectorRDouble1D & GetFaceNormalX() { return faceNormalX; }
    PHVectorRDouble1D & GetFaceNormalY() { return faceNormalY; }
    PHVectorRDouble1D & GetFaceNormalZ() { return faceNormalZ; }

    int GetNTotalCompositeFaces() { return static_cast<int>(compositeFace2Node.size()); }

    set < MultiIndex > & GetFaceList() { return faceList; }
    vector< SimpleBC * > & GetFaceBoundaryCondition() { return faceBoundaryCondition; }

    int GetNTotalFaces();
    int GetNBoundaryFaces();
    int GetNumberOfTotalEdgesWithMiddleNode();

    //! Find if face is exist.
    int FindFace(MultiIndex face);
    int FindFace(PHVectorInt1D face2Node);
    int FindCompositeFace(MultiIndex compositeFace);
    int FindEdgeWithMiddleNode(MultiIndex edgeWithMiddleNode);

    //! Insert face.
    void InsertFace(MultiIndex face);
    void InsertCompositeFace(MultiIndex compositeface);
    void InsertEdgeWithMiddleNode(MultiIndex edgeWithMiddleNode);

    //! Add face.
    int AddFace(PHVectorInt1D &face2NodeIn, int faceTypeIn, int parentFaceIndex);
    void ResizeFaceNumber(int faceNumber);

    //! Set default information of face, after refined.
    void SetDefaultFaceInformation(int faceIndex, int faceTypeIn, PHVectorInt1D &face2NodeIn, int parentFaceIndex);
    void AddFaceHierachical(int faceIndex, int parentFaceIndexIn);
    void ScanChildFace(PHVectorInt1D cell2Node, int cellType, int cellIndex);

private:
    PHVectorInt2D intermediateFace2Node;

};















}