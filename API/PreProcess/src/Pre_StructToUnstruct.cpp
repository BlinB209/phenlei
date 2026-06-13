#include "Pre_StructToUnstruct.h"
#include "MixGrid.h"
#include "GridType.h"
#include "Glb_Dimension.h"
#include "Pre_GridConversion.h"
#include "TK_Log.h"

using namespace std;
namespace PHSPACE
{
Pre_StructToUnstruct::Pre_StructToUnstruct(Grid **structGridIn, int nZonesIn)
{
    this->structGrid = structGridIn;
    this->numberOfZones = nZonesIn;

    structToUnstruct = 0;
}

Pre_StructToUnstruct::~Pre_StructToUnstruct()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        delete structToUnstruct[iZone];
    }
    delete [] structToUnstruct;
}

void Pre_StructToUnstruct::Run()
{
    InitMemory();
    CopyCoordinate();
    CopyIBlock();
    ConstructBCType();
    ConstructFaceTopo();
    CopyInterface();
    ConstructCellMapStrToUns();
}

void Pre_StructToUnstruct::ConstructCellMapStrToUns()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        if (structGrid[iZone]->Type() != STRUCTGRID) continue;
        StructGrid   *strGrid = StructGridCast  (structGrid[iZone]);
        int nI = strGrid->GetNI();
        int nJ = strGrid->GetNJ();
        int nK = strGrid->GetNK();

        Range I, J, K;
        GetRange(nI, nJ, nK, -1, 0, I, J, K);
        Int3D *cellMapStrToUns = new Int3D(I, J, K, fortranArray);
        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*cellMapStrToUns)(i, j, k) = -1;
                }
            }
        }

        int nTotalCell = strGrid->GetNTotalCell();

        //! Gost cell.
        StructBCSet *structBCSet = strGrid->GetStructBCSet();
        int nBCRegion = structBCSet->GetnBCRegion();
        int iBoundFace = 0;
        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);
            string bcName = structBC->GetBCName();

            int ist1, ied1, jst1, jed1, kst1, ked1;
            structBC->GetIJKRegion(ist1, ied1, jst1, jed1, kst1, ked1);

            int iSurface         = structBC->GetFaceDirection() + 1;
            int leftOrRightIndex = structBC->GetFaceLeftOrRightIndex();

            for (int k = kst1; k <= ked1; ++ k)
            {
                for (int j = jst1; j <= jed1; ++ j)
                {
                    for (int i = ist1; i <= ied1; ++ i)
                    {
                        int iGost = i;
                        int jGost = j;
                        int kGost = k;

                        if (iSurface == X_DIR)
                        {
                            iGost = iGost + leftOrRightIndex;
                        }
                        else if (iSurface == Y_DIR)
                        {
                            jGost = jGost + leftOrRightIndex;
                        }
                        else if (iSurface == Z_DIR)
                        {
                            kGost = kGost + leftOrRightIndex;
                        }

                        (*cellMapStrToUns)(iGost, jGost, kGost) = nTotalCell + iBoundFace;
                        iBoundFace ++;

                    }
                }
            }
        }

        //! Interior cell.
        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        strGrid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd);

        for (int kCell = kCellStart; kCell <= kCellEnd; ++ kCell)
        {
            for (int jCell = jCellStart; jCell <= jCellEnd; ++ jCell)
            {
                for (int iCell = iCellStart; iCell <= iCellEnd; ++ iCell)
                {
                    int unsCellIndex = GetCellIndex(strGrid, iCell, jCell, kCell);
                    (*cellMapStrToUns)(iCell, jCell, kCell) = unsCellIndex;
                }
            }
        }

        strGrid->UpdateDataPtr("cellMapStrToUns", cellMapStrToUns);
    }

}

void Pre_StructToUnstruct::CopyInterface()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        if (structGrid[iZone]->Type() != STRUCTGRID) continue;
        UnstructGrid *unsGrid = UnstructGridCast(structToUnstruct[iZone]);
        InterfaceInfo *iifo = structGrid[iZone]->GetInterfaceInfo();
        unsGrid->SetInterfaceInfo(iifo);
        InterfaceFields *interfaceFields = structGrid[iZone]->GetInterfaceFields();
        unsGrid->SetInterfaceFields(interfaceFields);
    }
}

Grid **Pre_StructToUnstruct::GetUnstructGrids()
{
    return this->structToUnstruct;
}

void Pre_StructToUnstruct::InitMemory()
{
    structToUnstruct = new Grid *[numberOfZones];
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        if (structGrid[iZone]->Type() != STRUCTGRID) continue;
        int gridIndex = structGrid[iZone]->GetZoneID();
        GridID *index = new GridID(gridIndex);
        Grid *grid;
        grid = CreateGridGeneral(UNSTRUCTGRID, index, 0, GetDim());
        structToUnstruct[iZone] = grid;
    }
}

void Pre_StructToUnstruct::CopyIBlock()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        if (structGrid[iZone]->Type() != STRUCTGRID) continue;
        Grid *strGrid = structGrid[iZone];
        Grid *unsGrid = structToUnstruct[iZone];
        int iBlock = strGrid->GetIBlock();
        unsGrid->SetIBlock(iBlock);
    }
}

void Pre_StructToUnstruct::CopyCoordinate()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        if (structGrid[iZone]->Type() != STRUCTGRID) continue;
        StructGrid   *strGrid = StructGridCast  (structGrid[iZone]);
        UnstructGrid *unsGrid = UnstructGridCast(structToUnstruct[iZone]);

        int nTotalNode = strGrid->GetNTotalNode();
        int nTotalFace = strGrid->GetNTotalFace();
        int nTotalCell = strGrid->GetNTotalCell();

        unsGrid->SetNTotalNode(nTotalNode);
        unsGrid->SetNTotalFace(nTotalFace);
        unsGrid->SetNTotalCell(nTotalCell);

        RDouble *x = new RDouble[nTotalNode];
        RDouble *y = new RDouble[nTotalNode];
        RDouble *z = new RDouble[nTotalNode];
        unsGrid->SetX(x);
        unsGrid->SetY(y);
        unsGrid->SetZ(z);

        RDouble *xx = strGrid->GetX();
        RDouble *yy = strGrid->GetY();
        RDouble *zz = strGrid->GetZ();

        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            x[iNode] = xx[iNode];
            y[iNode] = yy[iNode];
            z[iNode] = zz[iNode];
        }
    }
}

void Pre_StructToUnstruct::ConstructBCType()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        if (structGrid[iZone]->Type() != STRUCTGRID) continue;
        StructGrid *strGrid = StructGridCast(structGrid[iZone]);
        UnstructGrid *unsGrid = UnstructGridCast(structToUnstruct[iZone]);

        int nBoundFace = strGrid->GetNBoundFace();
        unsGrid->SetNBoundFace(nBoundFace);
        UnstructBCSet **bcrs = new UnstructBCSet * [nBoundFace];

        int boundFaceIndex = 0;
        StructBCSet *structBCSet = strGrid->GetStructBCSet();
        int nBCRegion = structBCSet->GetnBCRegion();

        int nBCRegionUnstruct = static_cast<int>(nBCRegion);
        unsGrid->CreateUnstructBCSet(nBCRegionUnstruct);
        UnstructBCSet *unstructBCSet = unsGrid->GetUnstructBCSet();
        int *bcRegionIDofBCFace = new int[nBoundFace];

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);
            int BCType = structBC->GetBCType();
            string bcName = structBC->GetBCName();

            UnstructBC *unstructBC = new UnstructBC(iBCRegion);
            unstructBCSet->SetBCRegion(iBCRegion, unstructBC);

            unstructBC->SetBCType(BCType);
            unstructBC->SetBCName(bcName);
            int ist, ied, jst, jed, kst, ked;
            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        bcRegionIDofBCFace[boundFaceIndex] = iBCRegion;
                        vector<int>* faceIndex = unstructBC->GetFaceIndex();
                        faceIndex->push_back(boundFaceIndex);

                        bcrs[boundFaceIndex] = new UnstructBCSet();
                        bcrs[boundFaceIndex]->SetBCName(bcName);
                        bcrs[boundFaceIndex]->SetKey(BCType);

                        boundFaceIndex ++;
                    }
                }
            }
        }
        unstructBCSet->SetBCFaceInfo(bcRegionIDofBCFace);
        unsGrid->SetBCRecord(bcrs);
    }
}

void Pre_StructToUnstruct::ConstructFaceTopo()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        if (structGrid[iZone]->Type() != STRUCTGRID) continue;
        StructGrid *strGrid = StructGridCast(structGrid[iZone]);
        UnstructGrid *unsGrid = UnstructGridCast(structToUnstruct[iZone]);

        int nTotalFace = unsGrid->GetNTotalFace();
        int nTotalCell = unsGrid->GetNTotalCell();

        int *node_number_of_each_face = new int[nTotalFace];

        int gridDim = GetDim();
        int face2NodeSize = 0;
        if (gridDim == TWO_D)
        {
            face2NodeSize = nTotalFace * 2;
            PHSPACE::SetField(node_number_of_each_face, 2, nTotalFace);
        }
        else
        {
            face2NodeSize = nTotalFace * 4;
            PHSPACE::SetField(node_number_of_each_face, 4, nTotalFace);
        }

        unsGrid->SetNodeNumberOfEachFace(node_number_of_each_face);

        int *face2Node = new int[face2NodeSize];
        unsGrid->SetFace2Node(face2Node);

        int *leftCellOfFace  = new int[nTotalFace];
        int *rightCellOfFace = new int[nTotalFace];
        unsGrid->SetLeftCellOfFace (leftCellOfFace);
        unsGrid->SetRightCellOfFace(rightCellOfFace);

        //! Boundary face
        int faceIndex = 0;
        int nodeCount = 0;
        StructBCSet *structBCSet = strGrid->GetStructBCSet();
        int nBCRegion = structBCSet->GetnBCRegion();
        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);

            int ist, ied, jst, jed, kst, ked;
            structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);
            int iSurface = structBC->GetFaceDirection() + 1;
            int leftOrRightIndex = structBC->GetFaceLeftOrRightIndex();

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        int leftCellIndex = GetCellIndex(strGrid, i, j, k);
                        leftCellOfFace [faceIndex] = leftCellIndex;
                        rightCellOfFace[faceIndex] = nTotalCell + faceIndex;

                        int iNodeIndex1 = 1, iNodeIndex2 = 1, iNodeIndex3, iNodeIndex4;
                        int jNodeIndex1 = 1, jNodeIndex2 = 1, jNodeIndex3, jNodeIndex4;
                        int kNodeIndex1, kNodeIndex2, kNodeIndex3, kNodeIndex4;

                        int iWall, jWall, kWall;
                        structBC->GetBoundaryFaceIndex(i, j, k, iWall, jWall, kWall);

                        if (gridDim == TWO_D)
                        {
                            if (iSurface == 1)
                            {
                                if (leftOrRightIndex == -1)
                                {
                                    iNodeIndex2 = iWall;
                                    jNodeIndex2 = jWall;

                                    iNodeIndex1 = iWall;
                                    jNodeIndex1 = jWall + 1;
                                }
                                else
                                {
                                    iNodeIndex1 = iWall;
                                    jNodeIndex1 = jWall;

                                    iNodeIndex2 = iWall;
                                    jNodeIndex2 = jWall + 1;
                                }
                            }
                            else if (iSurface == 2)
                            {
                                if (leftOrRightIndex == -1)
                                {
                                    iNodeIndex1 = iWall;
                                    jNodeIndex1 = jWall;
                                    iNodeIndex2 = iWall + 1;
                                    jNodeIndex2 = jWall;
                                }
                                else
                                {
                                    iNodeIndex2 = iWall;
                                    jNodeIndex2 = jWall;
                                    iNodeIndex1 = iWall + 1;
                                    jNodeIndex1 = jWall;
                                }
                            }

                            kNodeIndex1 = 1;
                            kNodeIndex2 = 1;

                            int nodeIndex1, nodeIndex2;
                            nodeIndex1 = GetNodeIndex(strGrid, iNodeIndex1, jNodeIndex1, kNodeIndex1);
                            nodeIndex2 = GetNodeIndex(strGrid, iNodeIndex2, jNodeIndex2, kNodeIndex2);

                            face2Node[nodeCount++] = nodeIndex1;
                            face2Node[nodeCount++] = nodeIndex2;
                        }
                        else
                        {
                            if (iSurface == 1)
                            {
                                iNodeIndex1 = iWall;
                                jNodeIndex1 = jWall;
                                kNodeIndex1 = kWall;

                                iNodeIndex2 = iWall;
                                jNodeIndex2 = jWall + 1;
                                kNodeIndex2 = kWall;

                                iNodeIndex3 = iWall;
                                jNodeIndex3 = jWall + 1;
                                kNodeIndex3 = kWall + 1;

                                iNodeIndex4 = iWall;
                                jNodeIndex4 = jWall;
                                kNodeIndex4 = kWall + 1;
                            }
                            else if (iSurface == 2)
                            {
                                iNodeIndex1 = iWall;
                                jNodeIndex1 = jWall;
                                kNodeIndex1 = kWall;

                                iNodeIndex2 = iWall;
                                jNodeIndex2 = jWall;
                                kNodeIndex2 = kWall + 1;

                                iNodeIndex3 = iWall + 1;
                                jNodeIndex3 = jWall;
                                kNodeIndex3 = kWall + 1;

                                iNodeIndex4 = iWall + 1;
                                jNodeIndex4 = jWall;
                                kNodeIndex4 = kWall;
                            }
                            else
                            {
                                iNodeIndex1 = iWall;
                                jNodeIndex1 = jWall;
                                kNodeIndex1 = kWall;

                                iNodeIndex2 = iWall + 1;
                                jNodeIndex2 = jWall;
                                kNodeIndex2 = kWall;

                                iNodeIndex3 = iWall + 1;
                                jNodeIndex3 = jWall + 1;
                                kNodeIndex3 = kWall;

                                iNodeIndex4 = iWall;
                                jNodeIndex4 = jWall + 1;
                                kNodeIndex4 = kWall;
                            }

                            int nodeIndex1, nodeIndex2, nodeIndex3, nodeIndex4;
                            nodeIndex1 = GetNodeIndex(strGrid, iNodeIndex1, jNodeIndex1, kNodeIndex1);
                            nodeIndex2 = GetNodeIndex(strGrid, iNodeIndex2, jNodeIndex2, kNodeIndex2);
                            nodeIndex3 = GetNodeIndex(strGrid, iNodeIndex3, jNodeIndex3, kNodeIndex3);
                            nodeIndex4 = GetNodeIndex(strGrid, iNodeIndex4, jNodeIndex4, kNodeIndex4);

                            if (leftOrRightIndex == 1)
                            {
                                face2Node[nodeCount++] = nodeIndex1;
                                face2Node[nodeCount++] = nodeIndex2;
                                face2Node[nodeCount++] = nodeIndex3;
                                face2Node[nodeCount++] = nodeIndex4;
                            }
                            else
                            {
                                face2Node[nodeCount++] = nodeIndex1;
                                face2Node[nodeCount++] = nodeIndex4;
                                face2Node[nodeCount++] = nodeIndex3;
                                face2Node[nodeCount++] = nodeIndex2;
                            }

                        }
                        faceIndex ++;
                    }
                }
            }
        }

        //! Interior face
        int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
        strGrid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd);

        for (int kCell = kCellStart; kCell <= kCellEnd; ++ kCell)
        {
            for (int jCell = jCellStart; jCell <= jCellEnd; ++ jCell)
            {
                for (int iCell = iCellStart; iCell <= iCellEnd; ++ iCell)
                {
                    int iNodeIndex1, iNodeIndex2, iNodeIndex3, iNodeIndex4;
                    int jNodeIndex1, jNodeIndex2, jNodeIndex3, jNodeIndex4;
                    int kNodeIndex1, kNodeIndex2, kNodeIndex3, kNodeIndex4;
                    int nodeIndex1, nodeIndex2, nodeIndex3, nodeIndex4;

                    if (gridDim == TWO_D)
                    {
                        if (iCell != iCellEnd)
                        {
                            iNodeIndex1 = iCell + 1;
                            jNodeIndex1 = jCell;

                            iNodeIndex2 = iCell + 1;
                            jNodeIndex2 = jCell + 1;

                            kNodeIndex1 = 1;
                            kNodeIndex2 = 1;

                            nodeIndex1 = GetNodeIndex(strGrid, iNodeIndex1, jNodeIndex1, kNodeIndex1);
                            nodeIndex2 = GetNodeIndex(strGrid, iNodeIndex2, jNodeIndex2, kNodeIndex2);

                            face2Node[nodeCount++] = nodeIndex1;
                            face2Node[nodeCount++] = nodeIndex2;

                            int leftCellIndex  = GetCellIndex(strGrid, iCell    , jCell, kCell);
                            int rightCellIndex = GetCellIndex(strGrid, iCell + 1, jCell, kCell);
                            leftCellOfFace [faceIndex] = leftCellIndex;
                            rightCellOfFace[faceIndex] = rightCellIndex;

                            faceIndex ++;
                        }

                        if (jCell != jCellEnd)
                        {
                            iNodeIndex2 = iCell;
                            jNodeIndex2 = jCell + 1;

                            iNodeIndex1 = iCell + 1;
                            jNodeIndex1 = jCell + 1;

                            kNodeIndex1 = 1;
                            kNodeIndex2 = 1;

                            nodeIndex1 = GetNodeIndex(strGrid, iNodeIndex1, jNodeIndex1, kNodeIndex1);
                            nodeIndex2 = GetNodeIndex(strGrid, iNodeIndex2, jNodeIndex2, kNodeIndex2);

                            face2Node[nodeCount++] = nodeIndex1;
                            face2Node[nodeCount++] = nodeIndex2;

                            int leftCellIndex  = GetCellIndex(strGrid, iCell, jCell    , kCell);
                            int rightCellIndex = GetCellIndex(strGrid, iCell, jCell + 1, kCell);
                            leftCellOfFace [faceIndex] = leftCellIndex;
                            rightCellOfFace[faceIndex] = rightCellIndex;

                            faceIndex ++;
                        }
                    }
                    else
                    {
                        //! I direction
                        if (iCell != iCellEnd)
                        {
                            iNodeIndex1 = iCell + 1;
                            jNodeIndex1 = jCell;
                            kNodeIndex1 = kCell;

                            iNodeIndex2 = iCell + 1;
                            jNodeIndex2 = jCell;
                            kNodeIndex2 = kCell + 1;

                            iNodeIndex3 = iCell + 1;
                            jNodeIndex3 = jCell + 1;
                            kNodeIndex3 = kCell + 1;

                            iNodeIndex4 = iCell + 1;
                            jNodeIndex4 = jCell + 1;
                            kNodeIndex4 = kCell;

                            nodeIndex1 = GetNodeIndex(strGrid, iNodeIndex1, jNodeIndex1, kNodeIndex1);
                            nodeIndex2 = GetNodeIndex(strGrid, iNodeIndex2, jNodeIndex2, kNodeIndex2);
                            nodeIndex3 = GetNodeIndex(strGrid, iNodeIndex3, jNodeIndex3, kNodeIndex3);
                            nodeIndex4 = GetNodeIndex(strGrid, iNodeIndex4, jNodeIndex4, kNodeIndex4);

                            face2Node[nodeCount++] = nodeIndex1;
                            face2Node[nodeCount++] = nodeIndex4;
                            face2Node[nodeCount++] = nodeIndex3;
                            face2Node[nodeCount++] = nodeIndex2;

                            int leftCellIndex  = GetCellIndex(strGrid, iCell    , jCell, kCell);
                            int rightCellIndex = GetCellIndex(strGrid, iCell + 1, jCell, kCell);
                            leftCellOfFace [faceIndex] = leftCellIndex;
                            rightCellOfFace[faceIndex] = rightCellIndex;

                            faceIndex ++;
                        }

                        //! J direction
                        if (jCell != jCellEnd)
                        {
                            iNodeIndex1 = iCell;
                            jNodeIndex1 = jCell + 1;
                            kNodeIndex1 = kCell;

                            iNodeIndex2 = iCell + 1;
                            jNodeIndex2 = jCell + 1;
                            kNodeIndex2 = kCell;

                            iNodeIndex3 = iCell + 1;
                            jNodeIndex3 = jCell + 1;
                            kNodeIndex3 = kCell + 1;

                            iNodeIndex4 = iCell;
                            jNodeIndex4 = jCell + 1;
                            kNodeIndex4 = kCell + 1;

                            nodeIndex1 = GetNodeIndex(strGrid, iNodeIndex1, jNodeIndex1, kNodeIndex1);
                            nodeIndex2 = GetNodeIndex(strGrid, iNodeIndex2, jNodeIndex2, kNodeIndex2);
                            nodeIndex3 = GetNodeIndex(strGrid, iNodeIndex3, jNodeIndex3, kNodeIndex3);
                            nodeIndex4 = GetNodeIndex(strGrid, iNodeIndex4, jNodeIndex4, kNodeIndex4);

                            face2Node[nodeCount++] = nodeIndex1;
                            face2Node[nodeCount++] = nodeIndex4;
                            face2Node[nodeCount++] = nodeIndex3;
                            face2Node[nodeCount++] = nodeIndex2;

                            int leftCellIndex  = GetCellIndex(strGrid, iCell, jCell    , kCell);
                            int rightCellIndex = GetCellIndex(strGrid, iCell, jCell + 1, kCell);
                            leftCellOfFace [faceIndex] = leftCellIndex;
                            rightCellOfFace[faceIndex] = rightCellIndex;

                            faceIndex ++;
                        }

                        //! K direction
                        if (kCell != kCellEnd)
                        {
                            iNodeIndex1 = iCell;
                            jNodeIndex1 = jCell;
                            kNodeIndex1 = kCell + 1;

                            iNodeIndex2 = iCell;
                            jNodeIndex2 = jCell + 1;
                            kNodeIndex2 = kCell + 1;

                            iNodeIndex3 = iCell + 1;
                            jNodeIndex3 = jCell + 1;
                            kNodeIndex3 = kCell + 1;

                            iNodeIndex4 = iCell + 1;
                            jNodeIndex4 = jCell;
                            kNodeIndex4 = kCell + 1;

                            nodeIndex1 = GetNodeIndex(strGrid, iNodeIndex1, jNodeIndex1, kNodeIndex1);
                            nodeIndex2 = GetNodeIndex(strGrid, iNodeIndex2, jNodeIndex2, kNodeIndex2);
                            nodeIndex3 = GetNodeIndex(strGrid, iNodeIndex3, jNodeIndex3, kNodeIndex3);
                            nodeIndex4 = GetNodeIndex(strGrid, iNodeIndex4, jNodeIndex4, kNodeIndex4);

                            face2Node[nodeCount++] = nodeIndex1;
                            face2Node[nodeCount++] = nodeIndex4;
                            face2Node[nodeCount++] = nodeIndex3;
                            face2Node[nodeCount++] = nodeIndex2;

                            int leftCellIndex  = GetCellIndex(strGrid, iCell, jCell, kCell);
                            int rightCellIndex = GetCellIndex(strGrid, iCell, jCell, kCell + 1);
                            leftCellOfFace [faceIndex] = leftCellIndex;
                            rightCellOfFace[faceIndex] = rightCellIndex;

                            faceIndex ++;
                        }
                    }
                }
            }
        }
    }
}

int Pre_StructToUnstruct::GetCellIndex(StructGrid *strGrid, int iCell, int jCell, int kCell)
{
    int nI = strGrid->GetNI();
    int nJ = strGrid->GetNJ();

    int faceSize = (nI - 1) * (nJ - 1);

    int cellIndex = faceSize * (kCell - 1) + (nI - 1) * (jCell - 1) + (iCell - 1);
    return cellIndex;
}

int Pre_StructToUnstruct::GetNodeIndex(StructGrid *strGrid, int iNode, int jNode, int kNode)
{
    int nI = strGrid->GetNI();
    int nJ = strGrid->GetNJ();

    int nodeSize = nI * nJ;

    int nodeIndex = nodeSize * (kNode - 1) + nI * (jNode - 1) + (iNode - 1);
    return nodeIndex;
}

}
