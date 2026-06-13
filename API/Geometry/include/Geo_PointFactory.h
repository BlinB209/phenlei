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
//! @file      Geo_PointFactory.h
//! @brief     Explain this file briefly.
//! @author    Bell, Zhang Yang.

#pragma once
#include "PHHeader.h"
#include "GRHash.h"
using namespace GRIDER_SPACE;

namespace PHSPACE
{

class Geo_PointFactory
{
public:
    Geo_PointFactory ();
    ~Geo_PointFactory();

public:
    typedef RDouble value_type;
    typedef Point3D point_type;
    typedef set< point_type >::iterator iterator;

protected:
    PHVectorRDouble1D pointRBFWeight;
    
    PHVectorRDouble1D pointNormalX;
    PHVectorRDouble1D pointNormalY;
    PHVectorRDouble1D pointNormalZ;

    PointHash pointList;
    PHVector1D < point_type > pointArray;

    PHVectorInt1D * computationalToGlobalNodeIndexMapping;
    PHVectorInt1D * globalToComputationalNodeIndexMapping;

    PHVector1D < set< int > > pointBoundaryConditionType;

    int oldNumberOfComputationalNodes;

    RDouble minDistance;

public:
    //!
    PHVectorRDouble1D & GetPointRBFWeight() { return pointRBFWeight; }

    //!
    PHVectorRDouble1D & GetPointNormalX() { return pointNormalX; }

    //!
    PHVectorRDouble1D & GetPointNormalY() { return pointNormalY; }

    //!
    PHVectorRDouble1D & GetPointNormalZ() { return pointNormalZ; }

    //!
    void SetMinDistance(RDouble dataIn) { this->minDistance = dataIn; }

    //! Get the minimum distance between two points.
    RDouble GetMinDistance() { return this->minDistance; }

    //!
    point_type & GetPoint(int id) { return pointArray[id]; }

    //!
    PHVectorInt1D & GetComputationalToGlobalNodeIndexMapping(){ return *computationalToGlobalNodeIndexMapping; }

    //!
    PHVectorInt1D & GetGlobalToComputationalNodeIndexMapping(){ return *globalToComputationalNodeIndexMapping; }

    //!
    void SetPoint(int id, point_type &point) { pointArray[id] = point; }

    //!
    bool FindPoint(point_type &point, iterator &iter);

    //!
    void AddPointDirectly(point_type &point, int boundaryConditionType = 0);

    //!
    int AddPoint(point_type *&point, bool &isNodeExist);

    //!
    int AddPointTest(point_type &point, int iNode, int numberOfNodes);

    //!
    PointHash & GetPointList() { return pointList; }

    //!
    PHVector1D < point_type > & GetPointArray() { return pointArray; }

    //!
    uint_t GetNumberOfComputationalNodes() { return (*computationalToGlobalNodeIndexMapping).size(); }    //! Green.

    //!
    uint_t GetNumberOfPoints() { return pointArray.size(); }

    //!
    Geo_PointFactory::point_type * GetCenterPoint(PHVectorInt1D &element, vector< int > &pointList);

    //!
    Geo_PointFactory::point_type CirclePoint(PHVectorInt1D &element, PHVectorInt1D &pointList);

    //!
    void SetBoundaryPoint(set< int > &bcVertex, int boundaryConditionType);

    //!
    bool PointWithBoundary(int index, int boundaryConditionType);

    //!
    void InsertPointBoundaryConditionType(int index, int boundaryConditionType);

    //!
    void SetOldNumberOfComputationalNodes(int oldNumberOfComputationalNodes);

    //!
    int GetOldNumberOfComputationalNodes();
};

Geo_PointFactory::point_type GetCenterPoint(PHVectorInt1D &element, PHVectorInt1D &pointList);
}
