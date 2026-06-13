#include "Mesh_DeformationSPRING.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "GeometryUnit.h"
#include "TK_Log.h"
#include "Geo_Element.h"

namespace PHSPACE
{

Mesh_DeformationSPRING::Mesh_DeformationSPRING(int nZonesIn, Grid **stationalGridIn) : 
    Mesh_Deformation(nZonesIn, stationalGridIn)
{

}

Mesh_DeformationSPRING::~Mesh_DeformationSPRING()
{

}

void Mesh_DeformationSPRING::SurfaceGridMoveByImput()
{
    vector < vector < RDouble > > nodeDisplacement = deformParameter->GetNodeDisplacement();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        DynamicGrid *grid      = deformedGrid[iZone];
        DYFace      *faceArray = grid->GetFaceArray();
        DYNode      *nodeArray = grid->GetNodeArray();

        int nBoundFace = grid->GetNBoundFace();
        int nTotalNode = grid->GetNTotalNode();

        int *nodeMark = new int[nTotalNode];
        SetField(nodeMark, 0, nTotalNode);

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

                for (int iControlPoint = 0; iControlPoint < nodeDisplacement.size(); ++ iControlPoint)
                {
                    RDouble disX = nodeArray[nodeIndex].x - nodeDisplacement[iControlPoint][0];
                    RDouble disY = nodeArray[nodeIndex].y - nodeDisplacement[iControlPoint][1];
                    RDouble disZ = nodeArray[nodeIndex].z - nodeDisplacement[iControlPoint][2];

                    RDouble dis = DISTANCE(disX, disY, disZ);
                    if (dis < EPSILON)
                    {
                        nodeMark[nodeIndex] = 1;

                        nodeArray[nodeIndex].xnew = nodeArray[nodeIndex].x + nodeDisplacement[iControlPoint][3];
                        nodeArray[nodeIndex].ynew = nodeArray[nodeIndex].y + nodeDisplacement[iControlPoint][4];
                        nodeArray[nodeIndex].znew = nodeArray[nodeIndex].z + nodeDisplacement[iControlPoint][5];

                        break;
                    }
                }
            }
        }

        delete [] nodeMark;    nodeMark = nullptr;
    }
}

void Mesh_DeformationSPRING::SurfaceGridMove(int iStep)
{
    bool imputNodeDisplacement = deformParameter->IsImputNodeDisplacement();
    if (imputNodeDisplacement)
    {
        SurfaceGridMoveByImput();
        return;
    }

    RDouble flapAngle = deformParameter->GetFlapAngle();
    int nDeformStep = deformParameter->GetNDeformStep();
    RDouble rotatePostionZ = deformParameter->GetRotatePostionZ();
    RDouble rotatePostionY = deformParameter->GetRotatePostionY();

    RDouble xNew, yNew, zNew;
    RDouble tNow = (iStep + 1) * 1.0 / nDeformStep;
    RDouble tOld = iStep * 1.0 / nDeformStep;
    RDouble omega = 2.0 * PI;
    RDouble amplitude = flapAngle;

    RDouble angle_step = amplitude * sin(omega * tNow);
    angle_step        -= amplitude * sin(omega * tOld);

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        DynamicGrid *grid = deformedGrid[iZone];
        DYNode *nodeArray = grid->GetNodeArray();

        int nTotalNode = grid->GetNTotalNode();
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            if (nodeArray[iNode].moveType != DYNAMIC_POINT_WALL)
            {
                continue;
            }

            RDouble currentRotatePosition;
            if (nodeArray[iNode].z > 0)
            {
                currentRotatePosition = rotatePostionZ;
            }
            else
            {
                currentRotatePosition = - rotatePostionZ;
            }

            xNew = nodeArray[iNode].x;
            yNew = -sin(angle_step) * (abs(nodeArray[iNode].z - currentRotatePosition)) + cos(angle_step) * (nodeArray[iNode].y - rotatePostionY);
            zNew =  cos(angle_step) * (abs(nodeArray[iNode].z - currentRotatePosition)) + sin(angle_step) * (nodeArray[iNode].y - rotatePostionY);

            yNew += rotatePostionY;
            if (nodeArray[iNode].z > 0)
            {
                zNew = currentRotatePosition + zNew;
            }
            else
            {
                zNew = currentRotatePosition - zNew;
            }

            nodeArray[iNode].xnew = xNew;
            nodeArray[iNode].ynew = yNew;
            nodeArray[iNode].znew = zNew;
        }

        //! Second segment.
        RDouble phase = - PI * 3 / 4.0;
        angle_step    = amplitude * sin(omega * tNow + phase);
        angle_step   -= amplitude * sin(omega * tOld + phase);

        int rotateNodeIndex = grid->GetRotateNodeIndex();
        if (rotateNodeIndex == -1)
        {
            return;
        }

        RDouble secondSegmentRotatePostionY = nodeArray[rotateNodeIndex].ynew;
        RDouble secondSegmentRotatePostionZ = nodeArray[rotateNodeIndex].znew;

        //! The rotate angle times of segment1.
        double times = 1.0;
        int *secondSegment = grid->GetSecondSegment();

        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            if (nodeArray[iNode].moveType != DYNAMIC_POINT_WALL)
            {
                continue;
            }

            if (secondSegment[iNode] != 1)
            {
                continue;
            }

            RDouble currentRotatePosition;
            if (nodeArray[iNode].z > 0)
            {
                currentRotatePosition = secondSegmentRotatePostionZ;
            }
            else
            {
                currentRotatePosition = - secondSegmentRotatePostionZ;
            }

            xNew = nodeArray[iNode].xnew;
            yNew = -sin(angle_step * times) * (abs(nodeArray[iNode].znew - currentRotatePosition)) + cos(angle_step * times) * (nodeArray[iNode].ynew - secondSegmentRotatePostionY);
            zNew =  cos(angle_step * times) * (abs(nodeArray[iNode].znew - currentRotatePosition)) + sin(angle_step * times) * (nodeArray[iNode].ynew - secondSegmentRotatePostionY);

            yNew += secondSegmentRotatePostionY;
            if (nodeArray[iNode].z > 0)
            {
                zNew = currentRotatePosition + zNew;
            }
            else
            {
                zNew = currentRotatePosition - zNew;
            }

            nodeArray[iNode].xnew = xNew;
            nodeArray[iNode].ynew = yNew;
            nodeArray[iNode].znew = zNew;
        }
    }
}

void Mesh_DeformationSPRING::TransformGrid()
{
    deformedGrid = new DynamicGrid *[numberOfZones];

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        deformedGrid[iZone] = new DynamicGrid(stationalGrid[iZone]);
        deformedGrid[iZone]->StationalData2DynamicData();
        deformedGrid[iZone]->ReconstrcutImplicitGeomeInfor();
        SetNodeAttribute(deformedGrid[iZone]);
    }

    SetSymmetryToZero();
}

void Mesh_DeformationSPRING::Deforming()
{
    string visualFileName = GlobalDataBase::GetStrParaFromDB("visualFileName");
    PostVisual(visualFileName, 0);

    PrintToWindow("  Start Deforming ... ", "\n");
    WriteLogFile("  Start Deforming ... ", "\n");

    int nDeformStep = deformParameter->GetNDeformStep();
    for (int iStep = 0; iStep < nDeformStep; ++ iStep)
    {
        PrintToWindow("  iStep =  ", iStep + 1, "\n");
        SurfaceGridMove(iStep);

        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            DynamicGrid *grid = deformedGrid[iZone];
            DYNode *nodeArray = grid->GetNodeArray();
            int nTotalNode = grid->GetNTotalNode();

            RDouble norm, norm1=1.0;
            RDouble normx, normy, normz;
            RDouble dxn, dyn, dzn;
            RDouble rn;
            RDouble dx, dy, dz;

            int springOrder = 2;
            RDouble epsilon = 1.0e-12;
            int numberOfSpringRelax = 4000;
            for (int iter = 1; iter < numberOfSpringRelax; ++ iter)
            {
                normx = 0;
                normy = 0;
                normz = 0;
                norm  = 0;

                for (int iNode = 0; iNode < nTotalNode; ++ iNode)
                {
                    if (nodeArray[iNode].moveType == DYNAMIC_POINT_FIELD)
                    {
                        dxn = nodeArray[iNode].xnew - nodeArray[iNode].x;
                        dyn = nodeArray[iNode].ynew - nodeArray[iNode].y;
                        dzn = nodeArray[iNode].znew - nodeArray[iNode].z;

                        dx = 0;
                        dy = 0;
                        dz = 0;
                        rn = 0;

                        for (int jNode = 0; jNode < nodeArray[iNode].nodeNumberAround; ++ jNode)
                        {
                            int nodeIndex   = nodeArray[iNode].node2node[jNode];
                            RDouble springK = nodeArray[iNode].node2nodeK[jNode];
                            springK         = pow(springK, springOrder * 1.0);

                            dx = dx + springK * (nodeArray[nodeIndex].xnew - nodeArray[nodeIndex].x);
                            dy = dy + springK * (nodeArray[nodeIndex].ynew - nodeArray[nodeIndex].y);
                            dz = dz + springK * (nodeArray[nodeIndex].znew - nodeArray[nodeIndex].z);
                            rn = rn + springK;
                        }

                        dx = dx / rn;
                        dy = dy / rn;
                        dz = dz / rn;

                        normx = normx + fabs(dx - dxn);
                        normy = normy + fabs(dy - dyn);
                        normz = normz + fabs(dz - dzn);
                        norm  = norm  + sqrt((dx - dxn) * (dx - dxn)
                              + (dy - dyn) * (dy - dyn)
                              + (dz - dzn) * (dz - dzn));

                        nodeArray[iNode].xnew = nodeArray[iNode].x + dx;
                        nodeArray[iNode].ynew = nodeArray[iNode].y + dy;
                        nodeArray[iNode].znew = nodeArray[iNode].z + dz;
                    }
                }

                if (iter == 1)
                {
                    cout << "#iter       normx     normy     normz    norm\n";
                    cout << "iter = " << iter << "	" << normx << "	" << normy << "	" << normz << "	" << norm << endl;
                    norm1 = norm;
                }

                if (iter % 100 == 0)
                {
                    cout << "iter = " << iter << "	" << normx << "	" << normy << "	" << normz << "	" << norm << endl;
                }

                if (norm / norm1 < epsilon)
                {
                    break;
                }
            }
        }

        PostTreat();
        PostVisual(visualFileName, iStep + 1);
    }

    PrintToWindow("  End Deforming ... ", "\n");
    WriteLogFile("  End Deforming ... ", "\n");
}

void Mesh_DeformationSPRING::MatchControlPoints()
{
    RDouble rotatePostion = this->deformParameter->GetRotatePostionZ();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        DynamicGrid *grid = deformedGrid[iZone];
        DYFace *faceArray = grid->GetFaceArray();
        DYNode *nodeArray = grid->GetNodeArray();

        int nBoundFace = grid->GetNBoundFace();
        int nTotalNode = grid->GetNTotalNode();

        int *nodeMark = new int[nTotalNode];
        SetField(nodeMark, 0, nTotalNode);

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

                if (fabs(nodeArray[nodeIndex].z) < rotatePostion)
                {
                    continue;
                }
                nodeMark[nodeIndex] = 1;

                nodeArray[nodeIndex].moveType = DYNAMIC_POINT_WALL;
            }
        }

        delete [] nodeMark;
    }
}

void Mesh_DeformationSPRING::PostTreat()
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

    SetSymmetryToZero();
}

}