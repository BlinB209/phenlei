#include "Mesh_DeformationRBF.h"
#include "Mesh_DeformationSPRING.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "GeometryUnit.h"
#include "TK_Log.h"
#include "Geo_Element.h"
//#include "Eigen/Eigen"
#include "MatrixLUSGS.h"

//using namespace Eigen;

namespace PHSPACE
{
//typedef Matrix< RDouble, Dynamic, Dynamic > RBFMatrix;

Mesh_DeformationRBF::Mesh_DeformationRBF(int nZonesIn, Grid **stationalGridIn) : 
    Mesh_Deformation(nZonesIn, stationalGridIn)
{
    nControlPoints = 0;
    controlPointVectors.resize(0);
    controlPointCoordinates.resize(0);
    controlPointToGlobalPointIndex.resize(0);

    secondSegment.resize(0);
    rotateNodeIndex.resize(0);

    influencePara       = 0;
    controlPointTree    = 0;
    numberOfReferenceCP = 0;
}

Mesh_DeformationRBF::~Mesh_DeformationRBF()
{
    delete controlPointTree;
}

void Mesh_DeformationRBF::SetNumberOfReferenceCP(int numberOfReferenceCPIn)
{
    this->numberOfReferenceCP = numberOfReferenceCPIn;
}

void Mesh_DeformationRBF::SetInfluencePara(RDouble influenceParaIn)
{
    this->influencePara = influenceParaIn;
}

void Mesh_DeformationRBF::TransformGrid()
{
    deformedGrid = new DynamicGrid *[numberOfZones];

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        deformedGrid[iZone] = new DynamicGrid(stationalGrid[iZone]);
        deformedGrid[iZone]->StationalData2DynamicData();
        SetNodeAttribute(deformedGrid[iZone]);
    }

    SetSymmetryToZero();
}

void Mesh_DeformationRBF::SurfaceGridMove(int iStep)
{
    RDouble flapAngle = deformParameter->GetFlapAngle();
    int nDeformStep = deformParameter->GetNDeformStep();
    RDouble rotatePostionZ = deformParameter->GetRotatePostionZ();
    RDouble rotatePostionY = deformParameter->GetRotatePostionY();

    RDouble xVector, yVector, zVector;
    RDouble xNew, yNew, zNew;
    RDouble tNow = (iStep + 1) * 1.0 / nDeformStep;
    RDouble tOld = iStep * 1.0 / nDeformStep;
    RDouble omega = 2.0 * PI;
    RDouble amplitude = flapAngle;

    RDouble angle_step = amplitude * sin(omega * tNow);
    angle_step        -= amplitude * sin(omega * tOld);

    for (int iNode = 0; iNode < nControlPoints; ++ iNode)
    {
        RDouble currentRotatePosition;
        if (controlPointCoordinates[iNode].Z() > 0)
        {
            currentRotatePosition = rotatePostionZ;
        }
        else
        {
            currentRotatePosition = -rotatePostionZ;
        }

        xNew =  controlPointCoordinates[iNode].X();
        yNew = -sin(angle_step) * (abs(controlPointCoordinates[iNode].Z() - currentRotatePosition)) + cos(angle_step) * (controlPointCoordinates[iNode].Y() - rotatePostionY);
        zNew =  cos(angle_step) * (abs(controlPointCoordinates[iNode].Z() - currentRotatePosition)) + sin(angle_step) * (controlPointCoordinates[iNode].Y() - rotatePostionY);

        yNew += rotatePostionY;
        if (controlPointCoordinates[iNode].Z() > 0)
        {
            zNew = currentRotatePosition + zNew;
        }
        else
        {
            zNew = currentRotatePosition - zNew;
        }

        xVector = xNew - controlPointCoordinates[iNode].X();
        yVector = yNew - controlPointCoordinates[iNode].Y();
        zVector = zNew - controlPointCoordinates[iNode].Z();

        Point3D node;
        node.SetX(xVector);
        node.SetY(yVector);
        node.SetZ(zVector);

        controlPointVectors[iNode] = node;

        PHVectorInt1D globalIndex = controlPointToGlobalPointIndex[iNode];
        int zoneIndex = globalIndex[0];
        int nodeIndex = globalIndex[1];

        DynamicGrid *grid = deformedGrid[zoneIndex];
        DYNode *nodeArray = grid->GetNodeArray();

        nodeArray[nodeIndex].xnew = xNew;
        nodeArray[nodeIndex].ynew = yNew;
        nodeArray[nodeIndex].znew = zNew;
    }

    //! Second segment.
    RDouble phase = - PI * 3 / 4.0;
    angle_step    = amplitude * sin(omega * tNow + phase);
    angle_step   -= amplitude * sin(omega * tOld + phase);

    //! The rotate angle times of segment1.
    double times = 1.0;

    for (int iNode = 0; iNode < nControlPoints; ++ iNode)
    {
        if (secondSegment[iNode] != 1)
        {
            continue;
        }

        int rotateNodeIndexCurrent = rotateNodeIndex[iNode];

        RDouble secondSegmentRotatePostionY = controlPointCoordinates[rotateNodeIndexCurrent].Y() + controlPointVectors[rotateNodeIndexCurrent].Y();
        RDouble secondSegmentRotatePostionZ = controlPointCoordinates[rotateNodeIndexCurrent].Z() + controlPointVectors[rotateNodeIndexCurrent].Z();

        RDouble xCoor = controlPointCoordinates[iNode].X() + controlPointVectors[iNode].X();
        RDouble yCoor = controlPointCoordinates[iNode].Y() + controlPointVectors[iNode].Y();
        RDouble zCoor = controlPointCoordinates[iNode].Z() + controlPointVectors[iNode].Z();

        RDouble currentRotatePosition;
        if (controlPointCoordinates[iNode].Z() > 0)
        {
            currentRotatePosition = secondSegmentRotatePostionZ;
        }
        else
        {
            currentRotatePosition = - secondSegmentRotatePostionZ;
        }

        xNew =  xCoor;
        yNew = -sin(angle_step * times) * (abs(zCoor - currentRotatePosition)) + cos(angle_step * times) * (yCoor - secondSegmentRotatePostionY);
        zNew =  cos(angle_step * times) * (abs(zCoor - currentRotatePosition)) + sin(angle_step * times) * (yCoor - secondSegmentRotatePostionY);

        yNew += secondSegmentRotatePostionY;
        if (controlPointCoordinates[iNode].Z() > 0)
        {
            zNew = currentRotatePosition + zNew;
        }
        else
        {
            zNew = currentRotatePosition - zNew;
        }

        xVector = xNew - controlPointCoordinates[iNode].X();
        yVector = yNew - controlPointCoordinates[iNode].Y();
        zVector = zNew - controlPointCoordinates[iNode].Z();

        Point3D node;
        node.SetX(xVector);
        node.SetY(yVector);
        node.SetZ(zVector);

        controlPointVectors[iNode] = node;

        PHVectorInt1D globalIndex = controlPointToGlobalPointIndex[iNode];
        int zoneIndex = globalIndex[0];
        int nodeIndex = globalIndex[1];

        DynamicGrid *grid = deformedGrid[zoneIndex];
        DYNode *nodeArray = grid->GetNodeArray();

        nodeArray[nodeIndex].xnew = xNew;
        nodeArray[nodeIndex].ynew = yNew;
        nodeArray[nodeIndex].znew = zNew;
    }
}

void Mesh_DeformationRBF::Deforming()
{
    string visualFileName = GlobalDataBase::GetStrParaFromDB("visualFileName");
    PostVisual(visualFileName, 0);

    PrintToWindow("  Start Deforming ... ", "\n");
    WriteLogFile("  Start Deforming ... ", "\n");

    //! RBF deform
    int nDeformStep = deformParameter->GetNDeformStep();
    for (int iStep = 0; iStep < nDeformStep; ++ iStep)
    {
        BuildControlPointTree();

        PrintToWindow("  iStep =  ", iStep + 1, "\n");
        SurfaceGridMove(iStep);

        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            DynamicGrid *grid = deformedGrid[iZone];
            DYNode *nodeArray = grid->GetNodeArray();
            int nTotalNode = grid->GetNTotalNode();

            int nNeedDeformPoints = 0;
            for (int iNode = 0; iNode < nTotalNode; ++ iNode)
            {
                if (nodeArray[iNode].moveType != DYNAMIC_POINT_FIELD) continue;
                ++ nNeedDeformPoints;
            }

            int count = 0;
            int interval = GetProgressInterval(nNeedDeformPoints, 20);
            for (int iNode = 0; iNode < nTotalNode; ++ iNode)
            {
                if (nodeArray[iNode].moveType != DYNAMIC_POINT_FIELD) continue;

                count ++;
                if (count > 0 && count % interval == 0)
                {
                    ProgressMonitorToWindows(count, nNeedDeformPoints, "Field point deformation");
                    ProgressMonitorToLogFile(count, nNeedDeformPoints, "Field point deformation");
                }

                //! Fine reference control point.
                PHVectorInt1D referenceCP;
                FindOutReferenceCP(nodeArray[iNode], referenceCP);
                if (referenceCP.size() == 0) continue;

                RBFDeformingKernel(grid, iNode, referenceCP);
            }
        }

        PostTreat();
        PostVisual(visualFileName, iStep + 1);
    }

    PrintToWindow("  End Deforming ... ", "\n");
    WriteLogFile("  End Deforming ... ", "\n");
}

void Mesh_DeformationRBF::BuildControlPointTree()
{
    // Compute the boundary.
    RDouble wallBound[6];
    for (int iDim = 0; iDim < 3; ++ iDim)
    {
        wallBound[iDim]     =  PHSPACE::LARGE;
        wallBound[iDim + 3] = -PHSPACE::LARGE;
    }

    for (int iPoint = 0; iPoint < nControlPoints; ++ iPoint)
    {
        RDouble x0 = controlPointCoordinates[iPoint].X();
        RDouble y0 = controlPointCoordinates[iPoint].Y();
        RDouble z0 = controlPointCoordinates[iPoint].Z();

        wallBound[0] = MIN(x0, wallBound[0]);
        wallBound[1] = MIN(y0, wallBound[1]);
        wallBound[2] = MIN(z0, wallBound[2]);

        wallBound[3] = MAX(x0, wallBound[3]);
        wallBound[4] = MAX(y0, wallBound[4]);
        wallBound[5] = MAX(z0, wallBound[5]);
    }

    wallTolerance = 0;
    for (int iDim = 0; iDim < 3; ++ iDim)
    {
        wallTolerance = MAX(wallTolerance, wallBound[iDim + 3] - wallBound[iDim]);
    }
    wallTolerance *= 0.01;

    for (int iDim = 0; iDim < 3; ++ iDim)
    {
        wallBound[iDim]     -= wallTolerance;
        wallBound[iDim + 3] += wallTolerance;
    }

    //! Build Tree.
    RDouble wallMin[6], wallMax[6];
    for (int iDim = 0; iDim < 3; ++ iDim)
    {
        wallMin[iDim] = wallBound[iDim];
        wallMax[iDim] = wallBound[iDim + 3];
    }
    for (int iDim = 3; iDim < 6; ++ iDim)
    {
        wallMin[iDim] = wallBound[iDim - 3];
        wallMax[iDim] = wallBound[iDim];
    }

    if (controlPointTree)
    {
        delete controlPointTree;
    }
    controlPointTree = new CPTree(6, wallMin, wallMax);

    RDouble triBox[6];

    for (int iPoint = 0; iPoint < nControlPoints; ++ iPoint)
    {
        for (int iDim = 0; iDim < 3; ++ iDim)
        {
            triBox[iDim]     =  PHSPACE::LARGE;
            triBox[iDim + 3] = -PHSPACE::LARGE;
        }

        RDouble xCoor = controlPointCoordinates[iPoint].X();
        RDouble yCoor = controlPointCoordinates[iPoint].Y();
        RDouble zCoor = controlPointCoordinates[iPoint].Z();

        triBox[0] = MIN(triBox[0], xCoor);
        triBox[1] = MIN(triBox[1], yCoor);
        triBox[2] = MIN(triBox[2], zCoor);

        triBox[3] = MAX(triBox[3], xCoor);
        triBox[4] = MAX(triBox[4], yCoor);
        triBox[5] = MAX(triBox[5], zCoor);

        CPNode *cpNode = new CPNode(6, triBox, iPoint);
        controlPointTree->AddNode(cpNode);
    }
}

void Mesh_DeformationRBF::FindOutReferenceCP(DYNode point, PHVectorInt1D &referenceCP)
{
    RDouble xCoor = point.x;
    RDouble yCoor = point.y;
    RDouble zCoor = point.z;

    RDouble *rmin, *rmax;
    rmin = controlPointTree->GetMin();
    rmax = controlPointTree->GetMax();
    RDouble pmin[6], pmax[6];
    RDouble dis = 0.1 * wallTolerance;

    // The maximum wanted wall points that in the box,
    // it is set to be sqrt(nTPoints), usually.
    int nCPLimitedMax = MAX(1 * (int) sqrt(nControlPoints * 1.0), 1);
    nCPLimitedMax     = MIN(nCPLimitedMax, nControlPoints);
    nCPLimitedMax     = MIN(nCPLimitedMax, 1000);
    int nCPLimitedMin = this->numberOfReferenceCP;

    // Find out by control points tree.
    bool finalFound = false;
    int cpMin = nCPLimitedMin;
    int cpMax = nCPLimitedMax;
    cpMax = MAX(cpMin, cpMax);

    CPTree::AdtNodeList nList;

    // Enlarge time, the number of probe time.
    // If lose once, enlarge the search scope then, by multiple 2 times.
    int nEnlargeTime = 1;
    int iEnlarge = 0;
    while (iEnlarge < nEnlargeTime)
    {
        ++ iEnlarge;

        pmin[0] = rmin[0];
        pmin[1] = rmin[1];
        pmin[2] = rmin[2];
        pmax[3] = rmax[3];
        pmax[4] = rmax[4];
        pmax[5] = rmax[5];

        pmin[3] = xCoor - dis;
        pmin[4] = yCoor - dis;
        pmin[5] = zCoor - dis;

        pmax[0] = xCoor + dis;
        pmax[1] = yCoor + dis;
        pmax[2] = zCoor + dis;

        RDouble size_min = 0;
        RDouble size_max = dis * 2.0;
        RDouble radius;

        int iSearch = 0;
        bool isFoundNearestCP = false;
        int nSearchTimes = 100;
        while (iSearch < nSearchTimes)
        {
            ++ iSearch;

            radius = 0.5 * (size_min + size_max);
            pmin[3] = xCoor - radius;
            pmin[4] = yCoor - radius;
            pmin[5] = zCoor - radius;

            pmax[0] = xCoor + radius;
            pmax[1] = yCoor + radius;
            pmax[2] = zCoor + radius;

            nList.clear();
            nList.reserve(cpMax + 1);
            controlPointTree->FindNodesInRegion(pmin, pmax, nList, cpMax + 1);

            int nNearCP = static_cast<int>(nList.size());
            if ((nNearCP > cpMin && nNearCP <= cpMax))
            {
                // case 1: the number of control points in the box is according with the anticipate.
                // case 2: the worst case, to ensure the success, accept this case, although 
                //         the control point found are to much.
                PHVectorInt1D     referenceCPTemp  (numberOfReferenceCP, PHSPACE::INVALID_INDEX);
                PHVectorRDouble1D referenceDistance(numberOfReferenceCP, PHSPACE::LARGE);

                for (int iPoint = 0; iPoint < nNearCP; ++ iPoint)
                {
                    int cpIndex = nList[iPoint]->GetData();

                    RDouble dx = xCoor - controlPointCoordinates[cpIndex].X();
                    RDouble dy = yCoor - controlPointCoordinates[cpIndex].Y();
                    RDouble dz = zCoor - controlPointCoordinates[cpIndex].Z();

                    // Second, find the minimum distance.
                    RDouble dist = SQR(dx, dy, dz);

                    // Find out the maximum distance and the points in the reference CP queue.
                    // In order to compare it with the others.
                    int     maxPoint    = PHSPACE::INVALID_INDEX;
                    RDouble maxDistance = PHSPACE::SMALL;
                    for (int iReference = 0; iReference < numberOfReferenceCP; ++ iReference)
                    {
                        if (referenceDistance[iReference] > maxDistance)
                        {
                            maxDistance = referenceDistance[iReference];
                            maxPoint    = iReference;
                        }
                    }

                    if (dist < maxDistance)
                    {
                        referenceCPTemp  [maxPoint] = cpIndex;
                        referenceDistance[maxPoint] = dist;
                    }
                }

                // Find out the real reference Control Points.
                int numberOfRealCP = 0;
                for (int iReference = 0; iReference < numberOfReferenceCP; ++ iReference)
                {
                    if (referenceCPTemp[iReference] != PHSPACE::INVALID_INDEX)
                    {
                        ++ numberOfRealCP;
                    }
                }

                referenceCP.clear();
                referenceCP.resize(numberOfRealCP);
                int count = 0;
                for (int iReference = 0; iReference < numberOfReferenceCP; ++ iReference)
                {
                    if (referenceCPTemp[iReference] != PHSPACE::INVALID_INDEX)
                    {
                        referenceCP[count] = referenceCPTemp[iReference];
                        ++ count;
                    }
                }

                if (numberOfRealCP == 0) TK_Exit::ExceptionExit("Error: The real Control Points do not exist!");

                isFoundNearestCP = true;
                break;
            } 
            else if (nNearCP > cpMax && (size_max - size_min) > 0.01 * dis)
            {
                // The box is too large and so lead to too large number of wall faces,
                // in order to reduce the wall faces, decrease the search radius.
                size_max = radius;
                continue;
            }
            else
            {
                // The box is too small and so lead to nearly non wall faces,
                // in order to get more the wall faces, increase the search radius.
                size_min = radius;
                size_max = size_max * 2.0;
                continue;
            }
        }

        if (isFoundNearestCP)
        {
            finalFound = true;
            break;
        }
        else
        {
            // Enlarge the limitation.
            cpMax *= 2;
            cpMin = MAX(1 * (int)sqrt(cpMax * 1.0), 1);
            continue;
        }
    }

    if (!finalFound)
    {
        referenceCP.clear();
    }
}

void Mesh_DeformationRBF::MatchControlPoints()
{
    RDouble rotatePostionZ = this->deformParameter->GetRotatePostionZ();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        DynamicGrid *grid = deformedGrid[iZone];
        DYFace *faceArray = grid->GetFaceArray();
        DYNode *nodeArray = grid->GetNodeArray();

        int nBoundFace = grid->GetNBoundFace();
        int nTotalNode = grid->GetNTotalNode();

        int *nodeMark = new int[nTotalNode];
        SetField(nodeMark, 0, nTotalNode);

        int *secondSegmentCurrent = grid->GetSecondSegment();
        int rotateNodeIndexCurrent = grid->GetRotateNodeIndex();

        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            int bctype = faceArray[iFace].bcType;
            if (bctype != PHENGLEI::SOLID_SURFACE)
            {
                continue;
            }

            for (int iNode = 0; iNode < faceArray[iFace].nodeNumber; ++ iNode)
            {
                int nodeIndex = faceArray[iFace].face2node[iNode];
                if (nodeMark[nodeIndex] == 1)
                {
                    continue;
                }

                if (fabs(nodeArray[nodeIndex].z) < rotatePostionZ)
                {
                    continue;
                }
                nodeMark[nodeIndex] = 1;

                nodeArray[nodeIndex].moveType = DYNAMIC_POINT_WALL;

                Point3D node;
                node.SetX(nodeArray[nodeIndex].x);
                node.SetY(nodeArray[nodeIndex].y);
                node.SetZ(nodeArray[nodeIndex].z);
                controlPointCoordinates.push_back(node);

                secondSegment.push_back(secondSegmentCurrent[nodeIndex]);

                //rotateNodeIndex.push_back(rotateNodeIndexCurrent);
                if (rotateNodeIndexCurrent == nodeIndex)
                {
                    rotateNodeIndexCurrent = static_cast<int>(controlPointCoordinates.size()) - 1;
                }


                PHVectorInt1D globalIndex;
                globalIndex.push_back(iZone);
                globalIndex.push_back(nodeIndex);

                controlPointToGlobalPointIndex.push_back(globalIndex);
            }
        }

        for (int iControlPoint = 0; iControlPoint < controlPointCoordinates.size(); ++ iControlPoint)
        {
            rotateNodeIndex.push_back(rotateNodeIndexCurrent);
        }

        delete [] nodeMark;
    }

    nControlPoints = static_cast<int>(controlPointCoordinates.size());
    controlPointVectors.resize(nControlPoints);
}

void Mesh_DeformationRBF::PostTreat()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        DynamicGrid *grid = deformedGrid[iZone];
        DYNode *nodeArray = grid->GetNodeArray();

        int nTotalNode = grid->GetNTotalNode();
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            if (nodeArray[iNode].moveType == STATIC_POINT)
            {
                continue;
            }

            nodeArray[iNode].x = nodeArray[iNode].xnew;
            nodeArray[iNode].y = nodeArray[iNode].ynew;
            nodeArray[iNode].z = nodeArray[iNode].znew;
        }
    }

    for (int iNode = 0; iNode < nControlPoints; ++ iNode)
    {
        PHVectorInt1D globalIndex = controlPointToGlobalPointIndex[iNode];
        int zoneIndex = globalIndex[0];
        int nodeIndex = globalIndex[1];

        DynamicGrid *grid = deformedGrid[zoneIndex];
        DYNode *nodeArray = grid->GetNodeArray();

        controlPointCoordinates[iNode].SetX(nodeArray[nodeIndex].x);
        controlPointCoordinates[iNode].SetY(nodeArray[nodeIndex].y);
        controlPointCoordinates[iNode].SetZ(nodeArray[nodeIndex].z);
    }

    SetSymmetryToZero();
}

void Mesh_DeformationRBF::RBFDeformingKernel(DynamicGrid *grid, int nodeIndex, PHVectorInt1D &referenceCP)
{
    DYNode *nodeArray = grid->GetNodeArray();
    int numberOfRealCP = static_cast<int>(referenceCP.size());

    Range M(0, numberOfRealCP - 1);
    RDouble2D *controlMatrix = new RDouble2D(M, M, fortranArray);
    RDouble2D &controlMatrix2D = *reinterpret_cast<RDouble2D *> (controlMatrix);
    RDouble2D *inverseMatrix = new RDouble2D(M, M, fortranArray);
    RDouble2D &inverseMatrix2D = *reinterpret_cast<RDouble2D *> (inverseMatrix);

    /*vector< vector< RDouble > > controlMatrix;
    controlMatrix.resize(numberOfRealCP);
    for (int iPoint = 0; iPoint < numberOfRealCP; ++ iPoint)
    {
        controlMatrix[iPoint].resize(numberOfRealCP);
    }*/

    // Compute the influence Radius.
    RDouble maxProjectDist = PHSPACE::SMALL;
    for (int iPoint = 0; iPoint < numberOfRealCP; ++ iPoint)
    {
        int nodeIndexTmp = referenceCP[iPoint];
        RDouble dx = controlPointVectors[nodeIndexTmp].X();
        RDouble dy = controlPointVectors[nodeIndexTmp].Y();
        RDouble dz = controlPointVectors[nodeIndexTmp].Z();
        RDouble dist = DISTANCE(dx, dy, dz);
        maxProjectDist = MAX(maxProjectDist, dist);
    }
    RDouble influenceRadius = maxProjectDist * this->influencePara;

    for (int iPoint = 0; iPoint < numberOfRealCP; ++ iPoint)
    {
        int iPointIndex = referenceCP[iPoint];
        for (int jPoint = 0; jPoint < numberOfRealCP; ++ jPoint)
        {
            int jPointIndex = referenceCP[jPoint];

            RDouble dx = controlPointCoordinates[iPointIndex].X() - controlPointCoordinates[jPointIndex].X();
            RDouble dy = controlPointCoordinates[iPointIndex].Y() - controlPointCoordinates[jPointIndex].Y();
            RDouble dz = controlPointCoordinates[iPointIndex].Z() - controlPointCoordinates[jPointIndex].Z();
            RDouble distance = DISTANCE(dx, dy, dz);

            if (iPoint != jPoint && distance <= 1.0e-10)
            {
                ostringstream oss;
                oss << "Error: Two overlapping control points are found! " << endl;
                TK_Exit::ExceptionExit(oss.str());
            }

            controlMatrix2D(iPoint, jPoint) = Wendland(distance, influenceRadius);
            //controlMatrix[iPoint][jPoint] = Wendland(distance, influenceRadius);
        }
    }

    //MatrixInverse(controlMatrix, numberOfRealCP);
    ObtainInverseDiagonalMatrix(controlMatrix, numberOfRealCP, inverseMatrix);

    // Compute Coefficient Vector.
    PHVectorRDouble1D coefficientX(numberOfRealCP);
    PHVectorRDouble1D coefficientY(numberOfRealCP);
    PHVectorRDouble1D coefficientZ(numberOfRealCP);

    for (int iPoint = 0; iPoint < numberOfRealCP; ++ iPoint)
    {
        coefficientX[iPoint] = 0;
        coefficientY[iPoint] = 0;
        coefficientZ[iPoint] = 0;

        for (int jPoint = 0; jPoint < numberOfRealCP; ++ jPoint)
        {
            int CPIndex = referenceCP[jPoint];

            RDouble dx = controlPointVectors[CPIndex].X();
            RDouble dy = controlPointVectors[CPIndex].Y();
            RDouble dz = controlPointVectors[CPIndex].Z();

            coefficientX[iPoint] += inverseMatrix2D(iPoint, jPoint) * dx;
            coefficientY[iPoint] += inverseMatrix2D(iPoint, jPoint) * dy;
            coefficientZ[iPoint] += inverseMatrix2D(iPoint, jPoint) * dz;
        }
    }

    // Interpolating from boundary to field.
    RDouble volumDx = 0;
    RDouble volumDy = 0;
    RDouble volumDz = 0;

    RDouble xCoor = nodeArray[nodeIndex].x;
    RDouble yCoor = nodeArray[nodeIndex].y;
    RDouble zCoor = nodeArray[nodeIndex].z;

    for (int iPoint = 0; iPoint < numberOfRealCP; ++ iPoint)
    {
        int nodeIndexOfCP = referenceCP[iPoint];

        RDouble dx = controlPointCoordinates[nodeIndexOfCP].X() - xCoor;
        RDouble dy = controlPointCoordinates[nodeIndexOfCP].Y() - yCoor;
        RDouble dz = controlPointCoordinates[nodeIndexOfCP].Z() - zCoor;
        RDouble distance = DISTANCE(dx, dy, dz);
        RDouble fai = Wendland(distance, influenceRadius);

        volumDx += fai * coefficientX[iPoint];
        volumDy += fai * coefficientY[iPoint];
        volumDz += fai * coefficientZ[iPoint];
    }

    RDouble xNew = xCoor + volumDx;
    RDouble yNew = yCoor + volumDy;
    RDouble zNew = zCoor + volumDz;

    nodeArray[nodeIndex].xnew = xNew;
    nodeArray[nodeIndex].ynew = yNew;
    nodeArray[nodeIndex].znew = zNew;

    delete controlMatrix;
    delete inverseMatrix;
}

RDouble Wendland(RDouble distance, RDouble radius)
{
    //return WendlandC0(distance, radius);
    return WendlandC2(distance, radius);
    //return WendlandC4(distance, radius);
}

RDouble WendlandC0(RDouble distance, RDouble radius)
{
    if (distance > radius)
    {
        return 0;
    }
    else
    {
        return (1 - distance / radius) * (1 - distance / radius);
    }
}

RDouble WendlandC2(RDouble distance, RDouble radius)
{
    if (distance > radius)
    {
        return 0;
    }
    else
    {
        RDouble enta = distance / (radius + 1.0e-12);
        return pow((1 - enta), 4) * (4 * enta + 1.0);
    }
}

RDouble WendlandC4(RDouble distance, RDouble radius)
{
    if (distance > radius)
    {
        return 0;
    }
    else
    {
        RDouble enta = distance / (radius + 1.0e-12);
        return pow((1 - enta), 6.0) * (35. * enta + 18. * enta + 3.0);
    }
}

//void MatrixInverse(vector< vector< RDouble > > &maxtrix, const int size)
//{
//    RBFMatrix maxtrixOfEigenForm(size, size);
//
//    for (int i = 0; i < size; ++ i)
//    {
//        for (int j = 0; j < size; ++ j)
//        {
//            maxtrixOfEigenForm(i, j) = maxtrix[i][j];
//        }
//    }
//
//    RBFMatrix inversedMatrixOfEigenForm = maxtrixOfEigenForm.inverse();
//
//    for (int i = 0; i < size; ++ i)
//    {
//        for (int j = 0; j < size; ++ j)
//        {
//            maxtrix[i][j] = inversedMatrixOfEigenForm(i, j);
//        }
//    }
//}

}