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
//! @file      InterfaceLinkStructure.h
//! @brief     Store interface information of grid, used in unstructured grid refining.
//! @author    Baka, Bell.

#pragma once
#include "Geo_Point.h"
#include "DataStruct_Sort.h"
using namespace std;

namespace PHSPACE
{
class Grid;
class ZoneInterface;

//! @brief FaceSearching achieve face search.
class FaceSearching
{
public:
    FaceSearching(const PHVectorInt1D &nodeIndexes, int faceIndex = 0);
    ~FaceSearching();

public:
    //! Face index.
    int faceIndex;

    //! Node index on face.
    PHVectorInt1D nodeIndexes;

    //! Node index on face, after reorder.
    PHVectorInt1D sortedNodeIndexes;

public:
    int GetFaceIndex() { return faceIndex; }
    uint_t GetNumberOfNodes() { return nodeIndexes.size(); }
    const PHVectorInt1D & GetNodeIndexes() { return nodeIndexes; }

};

//! @brief CompareFaceByMethod1 used to construct face list.
class CompareFaceByMethod1
{
public:
    bool operator() (const FaceSearching *lhs, const FaceSearching *rhs) const
    {
        return lhs->sortedNodeIndexes < rhs->sortedNodeIndexes;
    }
};

//! @brief FaceSearchingManager achieve face search.
class FaceSearchingManager
{
public:
    FaceSearchingManager();
    ~FaceSearchingManager();

protected:
    //! Face list, used to face search.
    PHVector1D < FaceSearching * > referenceFacesArray;
    set< FaceSearching *, CompareFaceByMethod1 > referenceFacesSet;

    //! Face status.
    PHVectorInt1D statusContainerOfFaces;

    //! Children face of each face.
    PHVectorInt2D childFaceIndexes;

    //! Node index map between child face and parent face.
    PHVectorInt2D relativeChildFaceNodeIndexes;

public:
    //! Add face.
    void AddFace(const PHVectorInt1D &faceNodeIndexes);

    //! Split QUAD_4 face to TRI_3 face.
    void ComputeNewFaceIndex();

protected:
    //! Find if face is exist.
    set< FaceSearching *, CompareFaceByMethod1 >::iterator FindFace(FaceSearching *faceSearching);

    //! Split QUAD_4 face to TRI_3 face.
    void SplitQuadrilateralToTriangle(FaceSearching *parentFaceSearching);

    //! Get face index and coordinates.
    void GetTriangleIndexes(FaceSearching *parentFaceSearching, PHVectorInt2D &localTriangleIndexesContainer, PHVectorInt2D &triangleIndexesContainer);
    void GetLocalTriangleIndexes(PHVectorInt2D &localTriangleIndexesContainer);

public:
    PHVectorInt2D & GetGlobalChildFaceIndexes() { return childFaceIndexes; }
    PHVectorInt2D & GetRelativeChildFaceNodeIndexes() { return relativeChildFaceNodeIndexes; }

};

//! @brief InterfaceLinkStructure store interface information of grid.
class InterfaceLinkStructure
{
public:
    typedef DataStruct_AdtTree < int, RDouble > AdtTree;

protected:
    //! Node coordinate of interface.
    AdtTree coordinateTree;

    //! FaceList, used to face search.
    set < DataStruct_Sort< PHVectorInt1D > > referenceInterfaceListForSearching;
    set < DataStruct_Sort< set< int > > > referenceInterfaceSetForSearching;

    //! Interface index map between local interface and global interface.
    PHVectorInt2D globalInterfaceIndexToZoneIndexesMapping;
    PHVectorInt2D globalInterfaceIndexToLocalInterfaceIndexesMapping;
    PHVectorInt2D localInterfaceIndexToGlobalInterfaceIndexMapping;

    //! Interface index map between local interface and global interface finally.
    PHVectorInt2D newGlobalInterfaceIndexToZoneIndexesMapping;
    PHVectorInt2D newGlobalInterfaceIndexToLocalInterfaceIndexesMapping;

    //! Global faces status.
    PHVectorInt1D statusContainerOfGlobalFaces;

    //! Face search.
    FaceSearchingManager *faceSearchingManager;

public:
    InterfaceLinkStructure(RDouble *pmin, RDouble *pmax, RDouble tolerance, int nZones);
    ~InterfaceLinkStructure();

public:
    //! Create interfacelink information
    void CreateLinkInformation(PHVectorInt1D &facePointIndexes, int zoneIndex, int iLocalInterfaceCount);

    //! Search neighbor grid index and interface index of in neighbor grid of each interface in current grid.
    void MatchInterfaceTopology(Grid *grid);

    //! Reconstruct interface topology, only used to QUAD_4 face.
    void ReconstructInterfaceTopology();

    //! Initialize size by number of interfaces of the grid.
    void InitializeInterface(Grid *grid);
    void InitializeInterface(ZoneInterface *zoneInterfaceTopology);

    //! Init new interface map by old interface map.
    void InitializeNewLocalGlobalInterfaceMapping();

    //! Update interface index map between local interface and global interface.
    void UpdateLocalGlobalInterfaceMapping();

    //! Init face status.
    void InitializeStatusContainerOfGlobalFaces();

    //! Add face.
    void AddFace(const PHVectorInt1D &facePointIndexes);

    //! Create interfaceLink information.
    void CreateLinkInformationBySet(PHVectorInt1D &facePointIndexes, int zoneIndex, int iLocalInterfaceCount);

public:
    int GetNumberOfInterfacePoints() { return coordinateTree.GetNodeNum(); }
    AdtTree & GetCoordinateTree() { return coordinateTree; }
    RDouble GetTolerance() const { return tolerance; }
    PHVectorInt1D & GetStatusContainerOfGlobalFaces() { return statusContainerOfGlobalFaces; }

    set < DataStruct_Sort< set< int > > >    & GetReferenceInterfaceSetForSearching () { return referenceInterfaceSetForSearching; }
    set < DataStruct_Sort< PHVectorInt1D > > & GetReferenceInterfaceListForSearching() { return referenceInterfaceListForSearching; }

    PHVectorInt2D & GetGlobalInterfaceIndexToZoneIndexesMapping             () { return globalInterfaceIndexToZoneIndexesMapping; }
    PHVectorInt2D & GetGlobalInterfaceIndexToLocalInterfaceIndexesMapping   () { return globalInterfaceIndexToLocalInterfaceIndexesMapping; }
    PHVectorInt2D & GetLocalInterfaceIndexToGlobalInterfaceIndexMapping     () { return localInterfaceIndexToGlobalInterfaceIndexMapping; }
    PHVectorInt2D & GetNewGlobalInterfaceIndexToZoneIndexesMapping          () { return newGlobalInterfaceIndexToZoneIndexesMapping; }
    PHVectorInt2D & GetNewGlobalInterfaceIndexToLocalInterfaceIndexesMapping() { return newGlobalInterfaceIndexToLocalInterfaceIndexesMapping; }

    FaceSearchingManager * GetFaceSearchingManager() { return faceSearchingManager; }

private:
    RDouble tolerance;
};

//! Create interfacelink information
void CreateLinkInformation(PHVectorInt1D &facePointIndexes, int zoneIndex, int iLocalInterfaceCount, InterfaceLinkStructure *interfaceLinkStructure);

}
