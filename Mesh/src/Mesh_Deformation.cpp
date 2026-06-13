#include "Mesh_Deformation.h"
#include "Mesh_DeformationSPRING.h"
#include "Mesh_DeformationRBF.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "GeometryUnit.h"
#include "TK_Log.h"

namespace PHSPACE
{

Mesh_Deformation::Mesh_Deformation(int nZonesIn, Grid **stationalGridIn)
{
    this->numberOfZones = nZonesIn;
    this->stationalGrid = stationalGridIn;

    this->deformedGrid = 0;
    this->deformParameter = 0;
}

Mesh_Deformation::~Mesh_Deformation()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        delete deformedGrid[iZone];
    }
    delete [] deformedGrid;
}

void Mesh_Deformation::SetDeformParameter(DeformParameter *deformParameterIn)
{
    this->deformParameter = deformParameterIn;
}

void Mesh_Deformation::Run()
{
    TransformGrid();
    MatchDynamicPoints();
    Deforming();
    ModifyOriginalGrid();
}

void Mesh_Deformation::SetNodeAttribute(DynamicGrid *grid)
{
    int nTotalNode = grid->GetNTotalNode();
    int nBoundFace = grid->GetNBoundFace();

    DYFace *faceArray = grid->GetFaceArray();
    DYNode *nodeArray = grid->GetNodeArray();

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        nodeArray[iNode].moveType = DYNAMIC_POINT_FIELD;
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bctype = faceArray[iFace].bcType;
        if (bctype != PHENGLEI::SYMMETRY)
        {
            int *face2node = faceArray[iFace].face2node;
            for (int iNode = 0; iNode < faceArray[iFace].nodeNumber; ++ iNode)
            {
                int nodeIndex = face2node[iNode];
                nodeArray[nodeIndex].moveType = STATIC_POINT;
            }
        }
    }
}

void Mesh_Deformation::SetSymmetryToZero()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        DynamicGrid *grid = deformedGrid[iZone];
        grid->SetSymmetryToZero();
    }
}

void Mesh_Deformation::MatchDynamicPoints()
{
    bool imputNodeDisplacement = deformParameter->IsImputNodeDisplacement();
    if (!imputNodeDisplacement)
    {
        MatchControlPoints();
        return;
    }
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
                        nodeMark [nodeIndex]          = 1;
                        nodeArray[nodeIndex].moveType = DYNAMIC_POINT_WALL;

                        break;
                    }
                }
            }
        }

        delete [] nodeMark; nodeMark = nullptr;
    }
}

void Mesh_Deformation::PostVisual(string visualFileName, int iStep)
{
    string fileName = AddSymbolToFileName(visualFileName, "_", iStep);
    fstream visualFile;
    OpenFile(visualFile, fileName, ios::out);
    visualFile.setf(ios::scientific);
    visualFile.precision(6);

    visualFile << "TITLE     = \"3D Grid\"\n";
    visualFile << "VARIABLES = \"X\"\n";
    visualFile << "            \"Y\"\n";
    visualFile << "            \"Z\"\n";

    const int OUT   = 1;
    const int UNOUT = 2;

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        DynamicGrid *grid = deformedGrid[iZone];

        int nTotalNode = grid->GetNTotalNode();
        int nTotalFace = grid->GetNTotalFace();
        int nTotalCell = grid->GetNTotalCell();

        DYNode *nodeArray = grid->GetNodeArray();
        DYFace *faceArray = grid->GetFaceArray();

        int nNodesOut, nFacesOut;

        int *nodeMark       = new int[nTotalNode];
        int *faceMark       = new int[nTotalFace];
        int *cellMark       = new int[nTotalCell];
        int *localNodeIndex = new int[nTotalNode];

        //! Dump wall
        for (int iNode = 0 ; iNode < nTotalNode; ++ iNode)
        {
            nodeMark[iNode] = UNOUT;
            localNodeIndex[iNode] = -1;
        }

        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            faceMark[iFace] = UNOUT;
        }

        nFacesOut = 0;
        int *face2node;
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            if (faceArray[iFace].bcType == PHENGLEI::SOLID_SURFACE)
            {
                ++ nFacesOut;
                face2node = faceArray[iFace].face2node;
                for (int iNode = 0; iNode < faceArray[iFace].nodeNumber; ++ iNode)
                {
                    int nodeIndex = face2node[iNode];
                    nodeMark[nodeIndex] = OUT;
                }
            }
        }

        nNodesOut = 0;
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            if (nodeMark[iNode] == OUT)
            {
                localNodeIndex[iNode] = nNodesOut;
                nNodesOut ++;
            }
        }

        if (nNodesOut)
        {
            visualFile << "ZONE T=\" Grid" << iZone << " SolidSurface \"\n";
            visualFile << "N= " << nNodesOut << ", E= " << nFacesOut << ", F=FEPOINT, ET=Quadrilateral\n";
            visualFile << "DT=(SINGLE SINGLE SINGLE)\n";

            for (int iNode = 0; iNode < nTotalNode; ++ iNode)
            {
                if (nodeMark[iNode] == OUT)
                {
                    visualFile << nodeArray[iNode].x << "	" << nodeArray[iNode].y << "	" << nodeArray[iNode].z << endl;
                }
            }

            for (int iFace = 0; iFace < nTotalFace; ++ iFace)
            {
                if (faceArray[iFace].bcType == PHENGLEI::SOLID_SURFACE)
                {
                    face2node = faceArray[iFace].face2node;
                    for (int iNode = 0; iNode < faceArray[iFace].nodeNumber; ++ iNode)
                    {
                        int nodeIndex = face2node[iNode];

                        nodeIndex = localNodeIndex[nodeIndex] + 1;
                        visualFile << nodeIndex << "    ";
                    }

                    for (int iNode = faceArray[iFace].nodeNumber; iNode < 4; ++ iNode)
                    {
                        visualFile << localNodeIndex[face2node[0]] + 1;
                        visualFile << "    ";
                    }
                    visualFile << endl;
                }
            }
        }

        int gridSlice = deformParameter->GetGridSlice();
        if (gridSlice)
        {
            int sliceAxis = deformParameter->GetSliceAxis();
            RDouble slicePosition = deformParameter->GetSlicePosition();

            for (int iNode = 0; iNode < nTotalNode; ++ iNode)
            {
                nodeMark[iNode] = UNOUT;
                localNodeIndex[iNode] = -1;
            }

            for (int iFace = 0; iFace < nTotalFace; ++ iFace)
            {
                faceMark[iFace] = UNOUT;
            }

            for (int iCell = 0; iCell < nTotalCell; ++ iCell)
            {
                cellMark[iCell] = UNOUT;
            }

            RDouble rmin, rmax;
            for (int iFace = 0; iFace < nTotalFace; ++ iFace)
            {
                rmin =  1.0e30;
                rmax = -1.0e30;

                face2node = faceArray[iFace].face2node;
                for (int iNode = 0; iNode < faceArray[iFace].nodeNumber; ++ iNode)
                {
                    int nodeIndex = face2node[iNode];
                    if (sliceAxis == X_DIR)
                    {
                        if (nodeArray[nodeIndex].x > rmax)
                        {
                            rmax = nodeArray[nodeIndex].x;
                        }

                        if (nodeArray[nodeIndex].x < rmin)
                        {
                            rmin = nodeArray[nodeIndex].x;
                        }
                    }
                    else if (sliceAxis == Y_DIR)
                    {
                        if (nodeArray[nodeIndex].y > rmax)
                        {
                            rmax = nodeArray[nodeIndex].y;
                        }

                        if (nodeArray[nodeIndex].y < rmin)
                        {
                            rmin = nodeArray[nodeIndex].y;
                        }
                    }
                    else
                    {
                        if (nodeArray[nodeIndex].z > rmax)
                        {
                            rmax = nodeArray[nodeIndex].z;
                        }

                        if (nodeArray[nodeIndex].z < rmin)
                        {
                            rmin = nodeArray[nodeIndex].z;
                        }
                    }
                }

                if ((rmax - slicePosition) * (rmin - slicePosition) <= 0.0)
                {
                    int leftCell = faceArray[iFace].leftCell;
                    cellMark[leftCell] = OUT;

                    int rightCell = faceArray[iFace].rightCell;
                    if (rightCell >= 0 && rightCell < nTotalCell)
                    {
                        cellMark[rightCell] = OUT;
                    }
                }
            }

            nFacesOut  = 0;
            for (int iFace = 0; iFace < nTotalFace; ++ iFace)
            {
                if (faceArray[iFace].bcType == PHENGLEI::INTERIOR)
                {
                    int leftCell  = faceArray[iFace].leftCell;
                    int rightCell = faceArray[iFace].rightCell;
                    if (cellMark[leftCell] == OUT || cellMark[rightCell] == OUT)
                    {
                        ++ nFacesOut;
                        faceMark[iFace] = OUT;
                        face2node = faceArray[iFace].face2node;
                        for (int iNode = 0; iNode < faceArray[iFace].nodeNumber; ++ iNode)
                        {
                            int nodeIndex = face2node[iNode];
                            nodeMark[nodeIndex] = OUT;
                        }
                    }
                }
            }

            nNodesOut = 0;
            for (int iNode = 0; iNode < nTotalNode; ++ iNode)
            {
                if (nodeMark[iNode] == OUT)
                {
                    localNodeIndex[iNode] = nNodesOut;
                    nNodesOut ++ ;
                }
            }

            if (nNodesOut)
            {
                visualFile << "ZONE T=\" Grid" << iZone << " SliceGrid \"\n";
                visualFile << "N= " << nNodesOut << ", E= " << nFacesOut << ", F=FEPOINT, ET=Quadrilateral\n";
                visualFile << "DT=(SINGLE SINGLE SINGLE)\n";

                for (int iNode = 0; iNode < nTotalNode; ++ iNode)
                {
                    if (nodeMark[iNode] == OUT)
                    {
                        visualFile << nodeArray[iNode].x << "	" << nodeArray[iNode].y << "	" << nodeArray[iNode].z << endl;
                    }
                }

                for (int iFace = 0; iFace < nTotalFace; ++ iFace)
                {
                    if (faceMark[iFace] == OUT)
                    {
                        face2node = faceArray[iFace].face2node;
                        for (int iNode = 0; iNode < faceArray[iFace].nodeNumber; ++ iNode)
                        {
                            int nodeIndex = face2node[iNode];

                            nodeIndex = localNodeIndex[nodeIndex] + 1;
                            visualFile << nodeIndex << "    ";
                        }

                        if (faceArray[iFace].nodeNumber == 3)
                        {
                            visualFile << localNodeIndex [face2node[2]] + 1;
                        }
                        visualFile << endl;
                    }
                }
            }
        }

        delete [] nodeMark;
        delete [] faceMark;
        delete [] cellMark;
        delete [] localNodeIndex;
    }

    CloseFile(visualFile);
}

void Mesh_Deformation::ModifyOriginalGrid()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        DynamicGrid *grid      = deformedGrid[iZone];
        DYNode      *nodeArray = grid->GetNodeArray();

        Grid    *originalGrid = stationalGrid[iZone];
        RDouble *x            = originalGrid->GetX();
        RDouble *y            = originalGrid->GetY();
        RDouble *z            = originalGrid->GetZ();

        int nTotalNode = grid->GetNTotalNode();
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            if (nodeArray[iNode].moveType == STATIC_POINT)
            {
                continue;
            }

            x[iNode] = nodeArray[iNode].x;
            y[iNode] = nodeArray[iNode].y;
            z[iNode] = nodeArray[iNode].z;
        }
    }
}


DeformParameter::DeformParameter()
{
    imputNodeDisplacement = false;
    nodeDisplacement.resize(0);
}

DeformParameter::~DeformParameter()
{

}

void DeformParameter::SetNDeformStep(int nDeformStepIn)
{
    this->nDeformStep = nDeformStepIn;
}

int DeformParameter::GetNDeformStep()
{
    return this->nDeformStep;
}

void DeformParameter::SetFlapAngle(RDouble flapAngleIn)
{
    this->flapAngle = flapAngleIn * PI / 180.0;
}

RDouble DeformParameter::GetFlapAngle()
{
    return this->flapAngle;
}

void DeformParameter::SetRotatePostionZ(RDouble rotatePostionZIn)
{
    this->rotatePostionZ = rotatePostionZIn;
}

RDouble DeformParameter::GetRotatePostionZ()
{
    return this->rotatePostionZ;
}

void DeformParameter::SetRotatePostionY(RDouble rotatePostionYIn)
{
    this->rotatePostionY = rotatePostionYIn;
}

RDouble DeformParameter::GetRotatePostionY()
{
    return this->rotatePostionY;
}

void DeformParameter::SetGridSlice(int gridSliceIn)
{
    this->gridSlice = gridSliceIn;
}

int DeformParameter::GetGridSlice()
{
    return this->gridSlice;
}

void DeformParameter::SetSliceAxis(int sliceAxisIn)
{
    this->sliceAxis = sliceAxisIn;
}

int DeformParameter::GetSliceAxis()
{
    return this->sliceAxis;
}

void DeformParameter::SetSlicePosition(RDouble slicePositionIn)
{
    this->slicePosition = slicePositionIn;
}

RDouble DeformParameter::GetSlicePosition()
{
    return this->slicePosition;
}

void DeformParameter::SetNodeDisplacement(vector < vector < RDouble > > nodeDisplacementIn)
{
    this->nodeDisplacement      = nodeDisplacementIn;
    this->imputNodeDisplacement = true;
}

vector < vector < RDouble > > DeformParameter::GetNodeDisplacement()
{
    return this->nodeDisplacement;
}

bool DeformParameter::IsImputNodeDisplacement()
{
    return this->imputNodeDisplacement;
}


void DeformGrid()
{
    using namespace PHMPI;

    //! Read stational grid.
    string grid_file;
    GlobalDataBase::GetData("stationalGridFile", &grid_file, PHSTRING, 1);
    GlobalDataBase::UpdateData("gridfile", &grid_file, PHSTRING, 1);

    Region *gridHandler = new Region();
    gridHandler->ReadGrid();

    GlobalBoundaryCondition::ReadGlobalBoundaryCondition();
    gridHandler->ReSetBoundaryConditionByGlobalBC();

    int nZones = GetNumberofGlobalZones();
    Grid **stationalGrid = new Grid *[nZones];

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        stationalGrid[iZone] = GetGrid(iZone, 0);
    }

    //! Get deform parameter.
    int nDeformStep = GlobalDataBase::GetIntParaFromDB("nDeformStep");
    RDouble flapAngle = GlobalDataBase::GetDoubleParaFromDB("flapAngle");
    RDouble rotatePostionZ = GlobalDataBase::GetDoubleParaFromDB("rotatePostionZ");
    RDouble rotatePostionY = GlobalDataBase::GetDoubleParaFromDB("rotatePostionY");
    int gridSlice = GlobalDataBase::GetIntParaFromDB("gridSlice");
    int sliceAxis = GlobalDataBase::GetIntParaFromDB("sliceAxis");
    RDouble slicePosition = GlobalDataBase::GetDoubleParaFromDB("slicePosition");

    DeformParameter *deformParameter = new DeformParameter();
    deformParameter->SetNDeformStep(nDeformStep);
    deformParameter->SetFlapAngle(flapAngle);
    deformParameter->SetRotatePostionZ(rotatePostionZ);
    deformParameter->SetRotatePostionY(rotatePostionY);
    deformParameter->SetGridSlice(gridSlice);
    deformParameter->SetSliceAxis(sliceAxis);
    deformParameter->SetSlicePosition(slicePosition);

    //! Grid deformation.
    Mesh_Deformation *gridDeformation = 0;
    int deformationMethod = GlobalDataBase::GetIntParaFromDB("deformationMethod");
    if (deformationMethod == SPRING)
    {
        gridDeformation = new Mesh_DeformationSPRING(nZones, stationalGrid);
    }
    else
    {
        int numberOfReferenceCP = GlobalDataBase::GetIntParaFromDB("numberOfReferenceCP");
        RDouble influencePara = GlobalDataBase::GetDoubleParaFromDB("influencePara");

        gridDeformation = new Mesh_DeformationRBF(nZones, stationalGrid);
        gridDeformation->SetNumberOfReferenceCP(numberOfReferenceCP);
        gridDeformation->SetInfluencePara(influencePara);
    }

    gridDeformation->SetDeformParameter(deformParameter);
    gridDeformation->Run();

    delete deformParameter;
    delete gridDeformation;
    delete gridHandler;
    delete [] stationalGrid;
}

}