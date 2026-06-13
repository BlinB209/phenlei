#include "Geo_PointFactory.h"
#include "Constants.h"
using namespace std;

namespace PHSPACE
{

Geo_PointFactory::Geo_PointFactory()
{
    computationalToGlobalNodeIndexMapping = new PHVectorInt1D;
    globalToComputationalNodeIndexMapping = new PHVectorInt1D;
    minDistance = 1.0e-10;
    oldNumberOfComputationalNodes = 0;
}

Geo_PointFactory::~Geo_PointFactory()
{
    delete computationalToGlobalNodeIndexMapping;
    delete globalToComputationalNodeIndexMapping;

    pointList.clear();
}

Geo_PointFactory::point_type Geo_PointFactory::CirclePoint(PHVectorInt1D & element, PHVectorInt1D & pointList)
{
    RDouble xc = zero;
    RDouble yc = zero;
    RDouble zc = zero;

    RDouble radius = one;

    uint_t nodeNumber = pointList.size();

    for (int iNode = 0; iNode < nodeNumber; ++ iNode)
    {
        Geo_PointFactory::point_type pt = GetPoint(element[ pointList[ iNode ] ]);
        xc += pt.X();
        yc += pt.Y();
        zc += pt.Z();
    }

    RDouble oNodeNumber = one / (nodeNumber + SMALL);

    xc *= oNodeNumber;
    yc *= oNodeNumber;
    zc *= oNodeNumber;

    RDouble osc = one / (DISTANCE(xc, yc, zc) + SMALL);

    RDouble x0 = xc * osc * radius;
    RDouble y0 = yc * osc * radius;
    RDouble z0 = zc * osc * radius;

    Geo_PointFactory::point_type center_point(x0, y0, z0);
    return center_point;
}

void Geo_PointFactory::SetBoundaryPoint(set< int > & bcVertex, int boundaryConditionType)
{
    for (set< int >::iterator iter = bcVertex.begin(); iter != bcVertex.end(); ++ iter)
    {
        pointBoundaryConditionType[* iter ].insert(boundaryConditionType);  //! Are the  point boundary condition type allowed to be repeated?
    }
}

bool Geo_PointFactory::PointWithBoundary(int index, int boundaryConditionType)
{
    set< int >::iterator iter = pointBoundaryConditionType[ index ].find(boundaryConditionType);
    if (iter != pointBoundaryConditionType[ index ].end())
    {
        return true;
    }
    return false;
}

void Geo_PointFactory::InsertPointBoundaryConditionType(int index, int boundaryConditionType)
{
    pointBoundaryConditionType[ index ].insert(boundaryConditionType);
}


Geo_PointFactory::point_type * Geo_PointFactory::GetCenterPoint(PHVectorInt1D & element, vector< int > & pointList)
{
    RDouble xc = zero;
    RDouble yc = zero;
    RDouble zc = zero;

    uint_t nodeNumber = pointList.size(); //! The node number of the center point.

    //! New
    //! Reorder the point, so that the point start from lower to super.
    set<int> orderedPoints;
    for (int iNode = 0; iNode < nodeNumber; ++ iNode)
    {
        int point = element[ pointList[ iNode ] ];
        orderedPoints.insert(point);
    }
    for (set<int>::iterator iter = orderedPoints.begin(); iter != orderedPoints.end(); ++ iter)
    {
        int point = * iter;
        Geo_PointFactory::point_type pt = GetPoint(point); //! Find out the vertex.
        xc += pt.X();
        yc += pt.Y();
        zc += pt.Z();
    }
    //*/

    RDouble oNodeNumber = one / (nodeNumber + SMALL);

    xc *= oNodeNumber;    //! Average.
    yc *= oNodeNumber;
    zc *= oNodeNumber;

    Geo_PointFactory::point_type * center_point = new Geo_PointFactory::point_type(xc, yc, zc);
    return center_point;
}

bool Geo_PointFactory::FindPoint(point_type & point, iterator & iter)
{
    return pointList.search(&point) == 0 ? false : true;
}

void Geo_PointFactory::AddPointDirectly(point_type & point, int boundaryConditionType)
{
    pointArray.push_back(point);
}

int Geo_PointFactory::AddPoint(point_type *& point, bool & isNodeExist)
{
    //! Insert the node into the tetNode, which is a hash table.
    int localIndex = pointList.insert(point, isNodeExist) - 1;

    //! Push the node into node list of GRAdvanceFront.
    if (isNodeExist)
    {
        //! Exist old node.
        delete point;
        point = pointList[localIndex];
    }else
    {
        //! The nodeIn is a new node.
        value_type tolerance = this->GetMinDistance();
        tolerance = MIN(tolerance / 1.0e10, 1.0e-13);
        point->SetTolerance(tolerance);

        //int index = this->GetNumberOfPoints();
        //point->SetID(index);

        point->SetID(localIndex);

        AddPointDirectly(*point);
    }

    return localIndex;
}

int Geo_PointFactory::AddPointTest(point_type & point, int iNode, int numberOfNodes)
{
    Geo_PointFactory::iterator iter;   //! Firstly,check whether the point exists, add in if it does not exist in order to avoid repetition.
    int resultsID;
    if (FindPoint(point, iter))
    {
        //cout << " iNode = " << iNode << " numberOfNodes = " << numberOfNodes << "\n";
        //cout << " Point Has Existed ! Point ID = " << iter->ID() << " " << point.ID() << "\n";
        resultsID = iter->ID();
        //return iter->ID();
    }
    else
    {
        int index = static_cast<int>(this->GetNumberOfPoints());
        point.SetID(index);
        AddPointDirectly(point);
        resultsID = index;
        //return index;
    }

    return resultsID;
}

void Geo_PointFactory::SetOldNumberOfComputationalNodes(int oldNumberOfComputationalNodes)
{
    this->oldNumberOfComputationalNodes = oldNumberOfComputationalNodes;
}

int Geo_PointFactory::GetOldNumberOfComputationalNodes()
{
    return oldNumberOfComputationalNodes;
}

}

