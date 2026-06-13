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
#pragma once
#include "PHMpi.h"
#include "Constants.h"
#include "ActionKey.h"
#include "Geo_UnstructGrid.h"
#include "Geo_GridManager.h"
#include "Glb_Dimension.h"
#include "Geo_FaceMetrics_Unstruct.h"
#include "Geo_CellMetrics_Unstruct.h"
#include "IO_FileName.h"
#include "PHIO.h"
#include "TK_Parse.h"
#include "TK_Exit.h"

namespace PHSPACE
{

GridManager::GridManager(UnstructGrid *unstructGridIn)
{
    oldFaceNormalX = NULL;
    oldFaceNormalY = NULL;
    oldFaceNormalZ = NULL;
    oldFaceArea = NULL;

    oldFaceCenterX = NULL;
    oldFaceCenterY = NULL;
    oldFaceCenterZ = NULL;

    cellVolumeN1 = NULL;
    cellVolumeN2 = NULL;

    oldFaceNormalVelocity = NULL;
    oldFaceVelocityX = NULL;
    oldFaceVelocityY = NULL;
    oldFaceVelocityZ = NULL;

    newFaceNormalVelocity = NULL;
    newFaceVelocityX = NULL;
    newFaceVelocityY = NULL;
    newFaceVelocityZ = NULL;

    blankIndex = NULL;

    localToGlobalNodeIndexMapping = NULL;
    globalCoordinateNew = NULL;

    this->unstructGrid = unstructGridIn;
    cellCenterX = unstructGrid->GetCellCenterX();
    cellCenterY = unstructGrid->GetCellCenterY();
    cellCenterZ = unstructGrid->GetCellCenterZ();

    faceNormalX = unstructGrid->GetFaceNormalX();
    faceNormalY = unstructGrid->GetFaceNormalY();
    faceNormalZ = unstructGrid->GetFaceNormalZ();

    cellVolume  = unstructGrid->GetCellVolume();
    faceCenterX = unstructGrid->GetFaceCenterX();
    faceCenterY = unstructGrid->GetFaceCenterY();
    faceCenterZ = unstructGrid->GetFaceCenterZ();
    faceArea    = unstructGrid->GetFaceArea();

    faceNormalVelocity = unstructGrid->GetFaceNormalVelocity();
    faceVelocityX = unstructGrid->GetFaceVelocityX();
    faceVelocityY = unstructGrid->GetFaceVelocityY();
    faceVelocityZ = unstructGrid->GetFaceVelocityZ();

    blankIndex = unstructGrid->GetBlankIndex();
}

GridManager::~GridManager()
{
    DeAllocateMetrics();

    DeAllocateALEMetrics();

    delete [] localToGlobalNodeIndexMapping;    localToGlobalNodeIndexMapping = NULL;
}


RDouble *GridManager::GetCellVolume(int iStep)
{
    //! set this here in temporary.
    if (iStep == 0)
    {
        return cellVolume;
    }
    else if (iStep == 1)
    {
        return cellVolumeN1;
    }
    else if (iStep == 2)
    {
        return cellVolumeN2;
    }

    return cellVolume;
}

void GridManager::AllocateMetrics()
{
    int numberOfFaces = this->unstructGrid->GetNTotalFace();

    if (IsNotAllocated(oldFaceNormalX))
    {
        oldFaceNormalX = new RDouble[numberOfFaces];
    }
    if (IsNotAllocated(oldFaceNormalY))
    {
        oldFaceNormalY = new RDouble[numberOfFaces];
    }
    if (IsNotAllocated(oldFaceNormalZ))
    {
        oldFaceNormalZ = new RDouble[numberOfFaces];
    }
}

void GridManager::DeAllocateMetrics()
{
    if (oldFaceNormalX)
    {
        delete [] oldFaceNormalX;
    }
    if (oldFaceNormalY)
    {
        delete [] oldFaceNormalY;
    }
    if (oldFaceNormalZ)
    {
        delete [] oldFaceNormalZ;
    }

    if (oldFaceCenterX)
    {
        delete [] oldFaceCenterX;
    }
    if (oldFaceCenterY)
    {
        delete [] oldFaceCenterY;
    }
    if (oldFaceCenterZ)
    {
        delete [] oldFaceCenterZ;
    }

    if (oldFaceArea)
    {
        delete [] oldFaceArea;
    }
}

void GridManager::AllocateALEMetrics()
{
    int numberOfCells = this->unstructGrid->GetNTotalCell();
    int numberOfBoundaryFaces = this->unstructGrid->GetNBoundFace();
    int numberOfFaces = this->unstructGrid->GetNTotalFace();
    int numberOfTotalCells = numberOfCells + numberOfBoundaryFaces;

    //ALE and other solver all need the data of faceNormalVelocity.
    //if (faceNormalVelocity == 0) faceNormalVelocity = new RDouble[numberOfFaces];

    //Open up the meory for cellVolumeN1Container\cellVolumeN2Container\oldFaceNormalVelocity only when the unsteady moving meshes conditions.
    if (IsNotAllocated(cellVolumeN1))
    {
        cellVolumeN1 = new RDouble[numberOfTotalCells];
    }
    if (IsNotAllocated(cellVolumeN2))
    {
        cellVolumeN2 = new RDouble[numberOfTotalCells];
    }

    if (IsNotAllocated(oldFaceVelocityX))
    {
        oldFaceVelocityX = new RDouble[numberOfFaces];
    }
    if (IsNotAllocated(oldFaceVelocityY))
    {
        oldFaceVelocityY = new RDouble[numberOfFaces];
    }
    if (IsNotAllocated(oldFaceVelocityZ))
    {
        oldFaceVelocityZ = new RDouble[numberOfFaces];
    }
    if (IsNotAllocated(oldFaceNormalVelocity))
    {
        oldFaceNormalVelocity = new RDouble[numberOfFaces];
    }

    if (IsNotAllocated(newFaceVelocityX))
    {
        newFaceVelocityX = new RDouble[numberOfFaces];
    }
    if (IsNotAllocated(newFaceVelocityY))
    {
        newFaceVelocityY = new RDouble[numberOfFaces];
    }
    if (IsNotAllocated(newFaceVelocityZ))
    {
        newFaceVelocityZ = new RDouble[numberOfFaces];
    }
    if (IsNotAllocated(newFaceNormalVelocity))
    {
        newFaceNormalVelocity = new RDouble[numberOfFaces];
    }
}

void GridManager::DeAllocateALEMetrics()
{
    if (cellVolumeN1)
    {
        delete [] cellVolumeN1;
    }
    if (cellVolumeN2)
    {
        delete [] cellVolumeN2;
    }

    if (oldFaceVelocityX)
    {
        delete [] oldFaceVelocityX;
    }
    if (oldFaceVelocityY)
    {
        delete [] oldFaceVelocityY;
    }
    if (oldFaceVelocityZ)
    {
        delete [] oldFaceVelocityZ;
    }
    if (oldFaceNormalVelocity)
    {
        delete [] oldFaceNormalVelocity;
    }

    if (newFaceVelocityX)
    {
        delete [] newFaceVelocityX;
    }
    if (newFaceVelocityY)
    {
        delete [] newFaceVelocityY;
    }
    if (newFaceVelocityZ)
    {
        delete [] newFaceVelocityZ;
    }
    if (newFaceNormalVelocity)
    {
        delete [] newFaceNormalVelocity;
    }
}

void GridManager::ComputeMetrics(ActionKey *actkey)
{
    ActionKey *actionKey = PHSPACE::GetCurrentActionKey();
    this->unstructGrid->AllocateMetrics(actionKey);

    if (PHSPACE::GetDim() == ONE_D)
    {
    }
    else if (PHSPACE::GetDim() == TWO_D)
    {
        this->ComputeFaceNormal2D();
        this->ComputeFaceCenter2D();
        this->ComputeCellCenterVol2D(actkey);
    }
    else
    {
        this->ComputeFaceNormal3D();
        this->ComputeFaceCenter3D();
        this->ComputeCellCenterVol3D(actkey);
    }

    RDouble *faceNormalX = unstructGrid->GetFaceNormalX();
    RDouble *faceNormalY = unstructGrid->GetFaceNormalY();
    RDouble *faceNormalZ = unstructGrid->GetFaceNormalZ();
    RDouble *faceArea    = unstructGrid->GetFaceArea();

    this->ClosureCheck(faceNormalX, faceArea);
    this->ClosureCheck(faceNormalY, faceArea);
    this->ClosureCheck(faceNormalZ, faceArea);
}

void GridManager::ClosureCheck(RDouble *faceNormalX, RDouble *faceArea)
{
    int numberOfCells = unstructGrid->GetNTotalCell();
    int numberOfFaces = unstructGrid->GetNTotalFace();

    int *leftCellIndexContainer = unstructGrid->GetLeftCellOfFace();
    int *rightCellIndexContainer = unstructGrid->GetRightCellOfFace();

    RDouble *geometricSummation = new RDouble[numberOfCells];

    PHSPACE::ZeroField(geometricSummation, numberOfCells);

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int leftCellIndex = leftCellIndexContainer[iFace];
        int rightCellIndex = rightCellIndexContainer[iFace];

        geometricSummation[leftCellIndex] += faceNormalX[iFace] * faceArea[iFace];

        if (rightCellIndex < 0 || rightCellIndex >= numberOfCells)
        {
            continue;
        }
        geometricSummation[rightCellIndex] -= faceNormalX[iFace] * faceArea[iFace];
    }

    RDouble totalGeometricSummation = 0.0;
    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        totalGeometricSummation += PHSPACE::ABS(geometricSummation[iCell]);
        if (PHSPACE::ABS(geometricSummation[iCell]) > 1.0e-8)
        {
            cout << "faceArea sum for cell " << iCell << " is " << geometricSummation[iCell] << "\n";
        }
    }
    delete [] geometricSummation;
    //cout << "zone " << this->GetZoneIndex() << " faceArea sums is " << total << "\n";
}

void GridManager::ComputeFaceNormal2D()
{
    RDouble *x = unstructGrid->GetX();
    RDouble *y = unstructGrid->GetY();

    int numberOfFaces = unstructGrid->GetNTotalFace();

    RDouble *faceNormalX = unstructGrid->GetFaceNormalX();
    RDouble *faceNormalY = unstructGrid->GetFaceNormalY();
    RDouble *faceNormalZ = unstructGrid->GetFaceNormalZ();
    RDouble *oldFaceNormalX = this->GetOldFaceNormalX();
    RDouble *oldFaceNormalY = this->GetOldFaceNormalY();
    RDouble *oldFaceNormalZ = this->GetOldFaceNormalZ();
    RDouble *faceArea = unstructGrid->GetFaceArea();

    int *faceNodeIndexContainer = unstructGrid->GetFace2Node();

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        oldFaceNormalX[iFace] = faceNormalX[iFace];
        oldFaceNormalY[iFace] = faceNormalY[iFace];
        oldFaceNormalZ[iFace] = faceNormalZ[iFace];

        int p1 = faceNodeIndexContainer[2 * iFace];
        int p2 = faceNodeIndexContainer[2 * iFace + 1];

        faceNormalX[iFace] = y[p2] - y[p1];
        faceNormalY[iFace] = x[p1] - x[p2];
        faceNormalZ[iFace] = 0.0;

        faceArea[iFace] = PHSPACE::DISTANCE(faceNormalX[iFace], faceNormalY[iFace], faceNormalZ[iFace]);
        RDouble reciprocalArea = 1.0 / (faceArea[iFace] + SMALL);
        faceNormalX[iFace] *= reciprocalArea;
        faceNormalY[iFace] *= reciprocalArea;
        faceNormalZ[iFace] *= reciprocalArea;
    }
}

void GridManager::ComputeFaceCenter2D()
{
    RDouble *x = unstructGrid->GetX();
    RDouble *y = unstructGrid->GetY();
    RDouble *z = unstructGrid->GetZ();

    int *faceNodeIndexContainer = unstructGrid->GetFace2Node();

    int numberOfFaces = unstructGrid->GetNTotalFace();

    RDouble *faceCenterX = unstructGrid->GetFaceCenterX();
    RDouble *faceCenterY = unstructGrid->GetFaceCenterY();
    RDouble *faceCenterZ = unstructGrid->GetFaceCenterZ();

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int p1 = faceNodeIndexContainer[2 * iFace];
        int p2 = faceNodeIndexContainer[2 * iFace + 1];
        faceCenterX[iFace] = half * (x[p1] + x[p2]);
        faceCenterY[iFace] = half * (y[p1] + y[p2]);
        faceCenterZ[iFace] = half * (z[p1] + z[p2]);
    }
}

void GridManager::ComputeCellCenterVol2D(ActionKey *actkey)
{
    int numberOfCells = this->unstructGrid->GetNTotalCell();
    int numberOfBoundaryFaces = this->unstructGrid->GetNBoundFace();
    int numberOfFaces = this->unstructGrid->GetNTotalFace();
    int numberOfTotalCells = numberOfCells + numberOfBoundaryFaces;

    int *leftCellIndexContainer = unstructGrid->GetLeftCellOfFace();
    int *rightCellIndexContainer = unstructGrid->GetRightCellOfFace();

    RDouble *cellCenterX = unstructGrid->GetCellCenterX();
    RDouble *cellCenterY = unstructGrid->GetCellCenterY();
    RDouble *cellCenterZ = unstructGrid->GetCellCenterZ();

    RDouble *faceCenterX = unstructGrid->GetFaceCenterX();
    RDouble *faceCenterY = unstructGrid->GetFaceCenterY();
    RDouble *faceCenterZ = unstructGrid->GetFaceCenterZ();

    RDouble *faceNormalX = unstructGrid->GetFaceNormalX();
    RDouble *faceNormalY = unstructGrid->GetFaceNormalY();
    RDouble *faceNormalZ = unstructGrid->GetFaceNormalZ();
    RDouble *faceArea = unstructGrid->GetFaceArea();
    RDouble *cellVolume = unstructGrid->GetCellVolume();

    PHSPACE::ZeroField(cellCenterX, numberOfTotalCells);
    PHSPACE::ZeroField(cellCenterY, numberOfTotalCells);
    PHSPACE::ZeroField(cellCenterZ, numberOfTotalCells);
    PHSPACE::ZeroField(cellVolume, numberOfTotalCells);

    // now the cell center and cellVolume
    // For each face, boundary faces first
    // For boundary faces
    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++iFace)
    {
        int leftCellIndex = leftCellIndexContainer[iFace];
        RDouble dot = (faceCenterX[iFace] * faceNormalX[iFace] +
                       faceCenterY[iFace] * faceNormalY[iFace] +
                       faceCenterZ[iFace] * faceNormalZ[iFace]) * faceArea[iFace];
        cellCenterX[leftCellIndex] += faceCenterX[iFace] * dot;
        cellCenterY[leftCellIndex] += faceCenterY[iFace] * dot;
        cellCenterZ[leftCellIndex] += faceCenterZ[iFace] * dot;
        cellVolume[leftCellIndex] += dot;
    }

    // For interior cell faces
    for (int iFace = numberOfBoundaryFaces; iFace < numberOfFaces; ++iFace)
    {
        int leftCellIndex = leftCellIndexContainer[iFace];
        int rightCellIndex = rightCellIndexContainer[iFace];
        RDouble dot = (faceCenterX[iFace] * faceNormalX[iFace] +
                       faceCenterY[iFace] * faceNormalY[iFace] +
                       faceCenterZ[iFace] * faceNormalZ[iFace]) * faceArea[iFace];
        cellCenterX[leftCellIndex] += faceCenterX[iFace] * dot;
        cellCenterY[leftCellIndex] += faceCenterY[iFace] * dot;
        cellCenterZ[leftCellIndex] += faceCenterZ[iFace] * dot;
        cellVolume[leftCellIndex] += dot;

        cellCenterX[rightCellIndex] -= faceCenterX[iFace] * dot;
        cellCenterY[rightCellIndex] -= faceCenterY[iFace] * dot;
        cellCenterZ[rightCellIndex] -= faceCenterZ[iFace] * dot;
        cellVolume[rightCellIndex] -= dot;
    }

    ostringstream oss;

    // Don't forget the coefficients
    int numberOfCellsHaveNegativeVolumes = 0;
    RDouble minvol = LARGE, maxvol = 0.0;
    int indexMinv = 0, indexMaxv = 0;
    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        RDouble tmp = 1.0 / (1.5 * cellVolume[iCell] + SMALL);
        cellCenterX[iCell] *= tmp;
        cellCenterY[iCell] *= tmp;
        cellCenterZ[iCell] *= tmp;
        cellVolume[iCell] *= half;

        if (minvol > cellVolume[iCell])
        {
            minvol = cellVolume[iCell];
            indexMinv = iCell;
        }
        if (maxvol < cellVolume[iCell])
        {
            maxvol = cellVolume[iCell];
            indexMaxv = iCell;
        }

        if (cellVolume[iCell] <= 0.0)
        {
            cellVolume[iCell] = -cellVolume[iCell];
            ++numberOfCellsHaveNegativeVolumes;
            oss << "  Warning: negative volume cell index " << iCell << ", (x, y) " << cellCenterX[iCell] << " " << cellCenterY[iCell] << " !\n";
        }
    }
    if (numberOfCellsHaveNegativeVolumes) oss << numberOfCellsHaveNegativeVolumes << " cells have negative volumes \n";

    // For ghost cells
    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++iFace)
    {
        int leftCellIndex = leftCellIndexContainer[iFace];
        int rightCellIndex = iFace + numberOfCells;
        if (faceArea[iFace] > SMALL)
        {
            RDouble tmp = 2.0 * ((cellCenterX[leftCellIndex] - faceCenterX[iFace]) * faceNormalX[iFace]
                               + (cellCenterY[leftCellIndex] - faceCenterY[iFace]) * faceNormalY[iFace]
                               + (cellCenterZ[leftCellIndex] - faceCenterZ[iFace]) * faceNormalZ[iFace]);
            cellCenterX[iFace + numberOfCells] = cellCenterX[leftCellIndex] - faceNormalX[iFace] * tmp;
            cellCenterY[iFace + numberOfCells] = cellCenterY[leftCellIndex] - faceNormalY[iFace] * tmp;
            cellCenterZ[iFace + numberOfCells] = cellCenterZ[leftCellIndex] - faceNormalZ[iFace] * tmp;
        }
        else
        {
            // Degenerated faces
            cellCenterX[iFace + numberOfCells] = -cellCenterX[leftCellIndex] + 2.0 * faceCenterX[iFace];
            cellCenterY[iFace + numberOfCells] = -cellCenterY[leftCellIndex] + 2.0 * faceCenterY[iFace];
            cellCenterZ[iFace + numberOfCells] = -cellCenterZ[leftCellIndex] + 2.0 * faceCenterZ[iFace];
        }
        cellVolume[rightCellIndex] = cellVolume[leftCellIndex];
    }

    PHSPACE::StreamToActionKey(actkey, oss);
}

void GridManager::ComputeFaceNormal3D()
{
    RDouble *x = unstructGrid->GetX();
    RDouble *y = unstructGrid->GetY();
    RDouble *z = unstructGrid->GetZ();

    int numberOfFaces = this->unstructGrid->GetNTotalFace();

    RDouble *faceNormalX = unstructGrid->GetFaceNormalX();
    RDouble *faceNormalY = unstructGrid->GetFaceNormalY();
    RDouble *faceNormalZ = unstructGrid->GetFaceNormalZ();
    RDouble *oldFaceNormalX = this->GetOldFaceNormalX();
    RDouble *oldFaceNormalY = this->GetOldFaceNormalY();
    RDouble *oldFaceNormalZ = this->GetOldFaceNormalZ();
    RDouble *faceArea = unstructGrid->GetFaceArea();

    int *faceNodeIndexContainer = unstructGrid->GetFace2Node();
    int *faceNodeNumberContainer = unstructGrid->GetNodeNumberOfEachFace();

    int nodePosition = 0;

    PHSPACE::SetField(oldFaceNormalX, faceNormalX, numberOfFaces);
    PHSPACE::SetField(oldFaceNormalY, faceNormalY, numberOfFaces);
    PHSPACE::SetField(oldFaceNormalZ, faceNormalZ, numberOfFaces);

    PHSPACE::ZeroField(faceNormalX, numberOfFaces);
    PHSPACE::ZeroField(faceNormalY, numberOfFaces);
    PHSPACE::ZeroField(faceNormalZ, numberOfFaces);

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int faceNodeNumber = faceNodeNumberContainer[iFace];
        for (int iNodeInFace = 0; iNodeInFace < faceNodeNumber; ++iNodeInFace)
        {
            int index1 = iNodeInFace;
            int index2 = (iNodeInFace + 1) % faceNodeNumberContainer[iFace];
            int p1 = faceNodeIndexContainer[nodePosition + index1];
            int p2 = faceNodeIndexContainer[nodePosition + index2];

            RDouble dx1 = x[p1];
            RDouble dy1 = y[p1];
            RDouble dz1 = z[p1];

            RDouble dx2 = x[p2];
            RDouble dy2 = y[p2];
            RDouble dz2 = z[p2];

            faceNormalX[iFace] += half * (dy1 * dz2 - dy2 * dz1);
            faceNormalY[iFace] += half * (dz1 * dx2 - dz2 * dx1);
            faceNormalZ[iFace] += half * (dx1 * dy2 - dx2 * dy1);
        }
        faceArea[iFace] = PHSPACE::DISTANCE(faceNormalX[iFace], faceNormalY[iFace], faceNormalZ[iFace]);
        nodePosition += faceNodeNumber;
    }

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        RDouble reciprocalArea = 1.0 / (faceArea[iFace] + SMALL);
        faceNormalX[iFace] *= reciprocalArea;
        faceNormalY[iFace] *= reciprocalArea;
        faceNormalZ[iFace] *= reciprocalArea;
    }
}

void GridManager::ComputeFaceCenter3D()
{
    RDouble *x = unstructGrid->GetX();
    RDouble *y = unstructGrid->GetY();
    RDouble *z = unstructGrid->GetZ();

    int *faceNodeIndexContainer = unstructGrid->GetFace2Node();
    int *faceNodeNumberContainer = unstructGrid->GetNodeNumberOfEachFace();

    int numberOfFaces = this->unstructGrid->GetNTotalFace();

    RDouble *faceCenterX = unstructGrid->GetFaceCenterX();
    RDouble *faceCenterY = unstructGrid->GetFaceCenterY();
    RDouble *faceCenterZ = unstructGrid->GetFaceCenterZ();

    RDouble *faceNormalX = unstructGrid->GetFaceNormalX();
    RDouble *faceNormalY = unstructGrid->GetFaceNormalY();
    RDouble *faceNormalZ = unstructGrid->GetFaceNormalZ();

    int nodePosition = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        RDouble x0 = 0.0;
        RDouble y0 = 0.0;
        RDouble z0 = 0.0;

        int faceNodeNumber = faceNodeNumberContainer[iFace];
        for (int iNodeInFace = 0; iNodeInFace < faceNodeNumber; ++iNodeInFace)
        {
            int index = faceNodeIndexContainer[nodePosition + iNodeInFace];
            x0 += x[index];
            y0 += y[index];
            z0 += z[index];
        }

        RDouble coef = 1.0 / faceNodeNumberContainer[iFace];

        x0 *= coef;
        y0 *= coef;
        z0 *= coef;

        faceCenterX[iFace] = 0.0;
        faceCenterY[iFace] = 0.0;
        faceCenterZ[iFace] = 0.0;
        RDouble sarea = 0.0;

        for (int iNodeInFace = 0; iNodeInFace < faceNodeNumber; ++iNodeInFace)
        {
            int index1 = iNodeInFace;
            int index2 = (iNodeInFace + 1) % faceNodeNumberContainer[iFace];
            int p1 = faceNodeIndexContainer[nodePosition + index1];
            int p2 = faceNodeIndexContainer[nodePosition + index2];

            RDouble dx1 = x[p1] - x0;
            RDouble dy1 = y[p1] - y0;
            RDouble dz1 = z[p1] - z0;

            RDouble dx2 = x[p2] - x0;
            RDouble dy2 = y[p2] - y0;
            RDouble dz2 = z[p2] - z0;

            RDouble x00 = third * (x[p1] + x[p2] + x0);
            RDouble y00 = third * (y[p1] + y[p2] + y0);
            RDouble z00 = third * (z[p1] + z[p2] + z0);

            RDouble anx = dy1 * dz2 - dy2 * dz1;
            RDouble any = dz1 * dx2 - dz2 * dx1;
            RDouble anz = dx1 * dy2 - dx2 * dy1;

            RDouble faceArea = half * DISTANCE(anx, any, anz);
            RDouble norm = anx * faceNormalX[iFace] + any * faceNormalY[iFace] + anz * faceNormalZ[iFace];

            if (norm < 0.0) faceArea = -faceArea;
            sarea += faceArea;

            faceCenterX[iFace] += faceArea * x00;
            faceCenterY[iFace] += faceArea * y00;
            faceCenterZ[iFace] += faceArea * z00;
        }

        RDouble osarea = 1.0 / (sarea + SMALL);

        faceCenterX[iFace] *= osarea;
        faceCenterY[iFace] *= osarea;
        faceCenterZ[iFace] *= osarea;

        nodePosition += faceNodeNumber;
    }
}

void GridManager::ComputeCellCenterVol3D(ActionKey *actkey)
{
    int numberOfCells = this->unstructGrid->GetNTotalCell();
    int numberOfBoundaryFaces = this->unstructGrid->GetNBoundFace();
    int numberOfFaces = this->unstructGrid->GetNTotalFace();
    int numberOfTotalCells = numberOfCells + numberOfBoundaryFaces;

    int *leftCellIndexContainer = unstructGrid->GetLeftCellOfFace();
    int *rightCellIndexContainer = unstructGrid->GetRightCellOfFace();

    int *faceNodeIndexContainer = unstructGrid->GetFace2Node();
    int *faceNodeNumberContainer = unstructGrid->GetNodeNumberOfEachFace();

    RDouble *x = unstructGrid->GetX();
    RDouble *y = unstructGrid->GetY();
    RDouble *z = unstructGrid->GetZ();

    RDouble *cellCenterX = unstructGrid->GetCellCenterX();
    RDouble *cellCenterY = unstructGrid->GetCellCenterY();
    RDouble *cellCenterZ = unstructGrid->GetCellCenterZ();

    RDouble *faceCenterX = unstructGrid->GetFaceCenterX();
    RDouble *faceCenterY = unstructGrid->GetFaceCenterY();
    RDouble *faceCenterZ = unstructGrid->GetFaceCenterZ();

    RDouble *faceNormalX = unstructGrid->GetFaceNormalX();
    RDouble *faceNormalY = unstructGrid->GetFaceNormalY();
    RDouble *faceNormalZ = unstructGrid->GetFaceNormalZ();

    RDouble *faceArea = unstructGrid->GetFaceArea();
    RDouble *cellVolume = unstructGrid->GetCellVolume();

    PHSPACE::ZeroField(cellCenterX, numberOfTotalCells);
    PHSPACE::ZeroField(cellCenterY, numberOfTotalCells);
    PHSPACE::ZeroField(cellCenterZ, numberOfTotalCells);
    PHSPACE::ZeroField(cellVolume,  numberOfTotalCells);

    int nodePosition = 0;

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int leftCellIndex = leftCellIndexContainer[iFace];
        int rightCellIndex = rightCellIndexContainer[iFace];

        int faceNodeNumber = faceNodeNumberContainer[iFace];

        for (int iNode = 0; iNode < faceNodeNumber; ++iNode)
        {
            int index1 = iNode;
            int index2 = (iNode + 1) % faceNodeNumber;
            int p2 = faceNodeIndexContainer[nodePosition + index1];
            int p3 = faceNodeIndexContainer[nodePosition + index2];

            RDouble x21 = x[p2] - faceCenterX[iFace];
            RDouble y21 = y[p2] - faceCenterY[iFace];
            RDouble z21 = z[p2] - faceCenterZ[iFace];

            RDouble x31 = x[p3] - faceCenterX[iFace];
            RDouble y31 = y[p3] - faceCenterY[iFace];
            RDouble z31 = z[p3] - faceCenterZ[iFace];

            RDouble nx = y21 * z31 - y31 * z21;
            RDouble ny = z21 * x31 - z31 * x21;
            RDouble nz = x21 * y31 - x31 * y21;

            RDouble fcX = faceCenterX[iFace] + x[p2] + x[p3];
            RDouble fcY = faceCenterY[iFace] + y[p2] + y[p3];
            RDouble fcZ = faceCenterZ[iFace] + z[p2] + z[p3];

            // Cell Center and Volume
            RDouble tmp = nx * fcX + ny * fcY + nz * fcZ;
            fcX *= tmp;
            fcY *= tmp;
            fcZ *= tmp;

            cellCenterX[leftCellIndex] += fcX;
            cellCenterY[leftCellIndex] += fcY;
            cellCenterZ[leftCellIndex] += fcZ;
            cellVolume[leftCellIndex] += tmp;

            cellCenterX[rightCellIndex] -= fcX;
            cellCenterY[rightCellIndex] -= fcY;
            cellCenterZ[rightCellIndex] -= fcZ;
            cellVolume[rightCellIndex] -= tmp;
        }
        nodePosition += faceNodeNumber;
    }

    ostringstream oss;

    // Don't forget the coefficients
    int numberOfErrorCells = 0;
    RDouble minvol = LARGE, maxvol = 0.0;
    int indexMinv = 0, indexMaxv = 0;
    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        RDouble tmp = 1.0 / (4.0 * cellVolume[iCell] + SMALL);
        cellCenterX[iCell] *= tmp;
        cellCenterY[iCell] *= tmp;
        cellCenterZ[iCell] *= tmp;
        cellVolume[iCell] /= 18.0;
        if (minvol > cellVolume[iCell])
        {
            minvol = cellVolume[iCell];
            indexMinv = iCell;
        }
        if (maxvol < cellVolume[iCell])
        {
            maxvol = cellVolume[iCell];
            indexMaxv = iCell;
        }

        if (cellVolume[iCell] <= 0.0)
        {
            cellVolume[iCell] = -cellVolume[iCell];
            ++numberOfErrorCells;
            oss << "  Warning: negative volume cell index " << iCell << ", (x, y, z) " << cellCenterX[iCell] << " " << cellCenterY[iCell] << " " << cellCenterZ[iCell] << " !\n";
        }
    }

    if (numberOfErrorCells) oss << numberOfErrorCells << " cells have negative vols \n";

    // For ghost cells
    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++iFace)
    {
        int leftCellIndex = leftCellIndexContainer[iFace];
        int rightCellIndex = iFace + numberOfCells;

        if (faceArea[iFace] > SMALL)
        {
            RDouble tmp = 2.0 * ((cellCenterX[leftCellIndex] - faceCenterX[iFace]) * faceNormalX[iFace]
                + (cellCenterY[leftCellIndex] - faceCenterY[iFace]) * faceNormalY[iFace]
                + (cellCenterZ[leftCellIndex] - faceCenterZ[iFace]) * faceNormalZ[iFace]);
            cellCenterX[rightCellIndex] = cellCenterX[leftCellIndex] - faceNormalX[iFace] * tmp;
            cellCenterY[rightCellIndex] = cellCenterY[leftCellIndex] - faceNormalY[iFace] * tmp;
            cellCenterZ[rightCellIndex] = cellCenterZ[leftCellIndex] - faceNormalZ[iFace] * tmp;

        }
        else
        {
            // Degenerated faces
            cellCenterX[rightCellIndex] = -cellCenterX[leftCellIndex] + 2.0 * faceCenterX[iFace];
            cellCenterY[rightCellIndex] = -cellCenterY[leftCellIndex] + 2.0 * faceCenterY[iFace];
            cellCenterZ[rightCellIndex] = -cellCenterZ[leftCellIndex] + 2.0 * faceCenterZ[iFace];
        }
        cellVolume[rightCellIndex] = cellVolume[leftCellIndex];
    }

    PHSPACE::StreamToActionKey(actkey, oss);
}

void GridManager::ComputeFaceNormal3DNew()
{
    RDouble *x = unstructGrid->GetX();
    RDouble *y = unstructGrid->GetY();
    RDouble *z = unstructGrid->GetZ();

    int numberOfFaces = this->unstructGrid->GetNTotalFace();

    RDouble *faceNormalX = this->GetFaceNormalX();
    RDouble *faceNormalY = this->GetFaceNormalY();
    RDouble *faceNormalZ = this->GetFaceNormalZ();
    RDouble *faceArea = this->GetFaceArea();

    int *faceNodeIndexContainer = unstructGrid->GetFace2Node();
    int *faceNodeNumberContainer = unstructGrid->GetNodeNumberOfEachFace();

    int nodePosition = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        RDouble xct = 0.0;
        RDouble yct = 0.0;
        RDouble zct = 0.0;

        int faceNodeNumber = faceNodeNumberContainer[iFace];
        for (int iNodeInFace = 0; iNodeInFace < faceNodeNumber; ++iNodeInFace)
        {
            int nodeIndex = faceNodeIndexContainer[nodePosition + iNodeInFace];
            xct += x[nodeIndex];
            yct += y[nodeIndex];
            zct += z[nodeIndex];
        }
        xct /= faceNodeNumberContainer[iFace];
        yct /= faceNodeNumberContainer[iFace];
        zct /= faceNodeNumberContainer[iFace];

        faceNormalX[iFace] = 0.0;
        faceNormalY[iFace] = 0.0;
        faceNormalZ[iFace] = 0.0;

        for (int iNodeInFace = 0; iNodeInFace < faceNodeNumber; ++iNodeInFace)
        {
            int p1 = faceNodeIndexContainer[nodePosition + iNodeInFace];
            int p2 = faceNodeIndexContainer[nodePosition + (iNodeInFace + 1) % faceNodeNumber];

            RDouble dx1 = xct - x[p1];
            RDouble dy1 = yct - y[p1];
            RDouble dz1 = zct - z[p1];

            RDouble dx2 = xct - x[p2];
            RDouble dy2 = yct - y[p2];
            RDouble dz2 = zct - z[p2];

            faceNormalX[iFace] += dy1 * dz2 - dy2 * dz1;
            faceNormalY[iFace] += dz1 * dx2 - dz2 * dx1;
            faceNormalZ[iFace] += dx1 * dy2 - dx2 * dy1;
        }

        faceNormalX[iFace] *= half;
        faceNormalY[iFace] *= half;
        faceNormalZ[iFace] *= half;

        faceArea[iFace] = PHSPACE::DISTANCE(faceNormalX[iFace], faceNormalY[iFace], faceNormalZ[iFace]);
        nodePosition += faceNodeNumber;
    }

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        RDouble reciprocalArea = 1.0 / (faceArea[iFace] + SMALL);
        faceNormalX[iFace] *= reciprocalArea;
        faceNormalY[iFace] *= reciprocalArea;
        faceNormalZ[iFace] *= reciprocalArea;
    }
}

void GridManager::ComputeFaceCenter3DNew()
{
    RDouble *x = unstructGrid->GetX();
    RDouble *y = unstructGrid->GetY();
    RDouble *z = unstructGrid->GetZ();

    int numberOfFaces = unstructGrid->GetNTotalFace();

    RDouble *faceCenterX = this->GetFaceCenterX();
    RDouble *faceCenterY = this->GetFaceCenterY();
    RDouble *faceCenterZ = this->GetFaceCenterZ();

    RDouble *faceNormalX = this->GetFaceNormalX();
    RDouble *faceNormalY = this->GetFaceNormalY();
    RDouble *faceNormalZ = this->GetFaceNormalZ();

    int *faceNodeIndexContainer  = unstructGrid->GetFace2Node();
    int *faceNodeNumberContainer = unstructGrid->GetNodeNumberOfEachFace();

    int nodePosition = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        RDouble xct = 0.0;
        RDouble yct = 0.0;
        RDouble zct = 0.0;

        int faceNodeNumber = faceNodeNumberContainer[iFace];
        for (int iNodeInFace = 0; iNodeInFace < faceNodeNumber; ++iNodeInFace)
        {
            int faceNodeIndex = faceNodeIndexContainer[nodePosition + iNodeInFace];
            xct += x[faceNodeIndex];
            yct += y[faceNodeIndex];
            zct += z[faceNodeIndex];
        }

        RDouble tmp = 1.0 / faceNodeNumberContainer[iFace];

        xct *= tmp;
        yct *= tmp;
        zct *= tmp;

        faceCenterX[iFace] = 0.0;
        faceCenterY[iFace] = 0.0;
        faceCenterZ[iFace] = 0.0;

        RDouble suma = 0.0;

        for (int iNodeInFace = 0; iNodeInFace < faceNodeNumber; ++iNodeInFace)
        {
            int p1 = faceNodeIndexContainer[nodePosition + iNodeInFace];
            int p2 = faceNodeIndexContainer[nodePosition + (iNodeInFace + 1) % faceNodeNumber];

            RDouble dx1 = xct - x[p1];
            RDouble dy1 = yct - y[p1];
            RDouble dz1 = zct - z[p1];

            RDouble dx2 = xct - x[p2];
            RDouble dy2 = yct - y[p2];
            RDouble dz2 = zct - z[p2];

            RDouble x00 = x[p1] + x[p2] + xct;
            RDouble y00 = y[p1] + y[p2] + yct;
            RDouble z00 = z[p1] + z[p2] + zct;

            RDouble anx = dy1 * dz2 - dy2 * dz1;
            RDouble any = dz1 * dx2 - dz2 * dx1;
            RDouble anz = dx1 * dy2 - dx2 * dy1;

            RDouble sss = PHSPACE::DISTANCE(anx, any, anz);
            RDouble norm = anx * faceNormalX[iFace] + any * faceNormalY[iFace] + anz * faceNormalZ[iFace];

            if (norm < 0.0) sss = -sss;
            suma += sss;

            faceCenterX[iFace] += sss * x00;
            faceCenterY[iFace] += sss * y00;
            faceCenterZ[iFace] += sss * z00;
        }

        tmp = third / (suma + SMALL);

        faceCenterX[iFace] *= tmp;
        faceCenterY[iFace] *= tmp;
        faceCenterZ[iFace] *= tmp;

        nodePosition += faceNodeNumber;
    }
}

void GridManager::ComputeCellCenterVol3DNew()
{
    int numberOfCells = this->unstructGrid->GetNTotalCell();
    int numberOfBoundaryFaces = this->unstructGrid->GetNBoundFace();
    int numberOfFaces = this->unstructGrid->GetNTotalFace();
    int numberOfTotalCells = numberOfCells + numberOfBoundaryFaces;

    int *leftCellIndexContainer = unstructGrid->GetLeftCellOfFace();
    int *rightCellIndexContainer = unstructGrid->GetRightCellOfFace();

    int *faceNodeIndexContainer = unstructGrid->GetFace2Node();
    int *faceNodeNumberContainer = unstructGrid->GetNodeNumberOfEachFace();

    RDouble *x = unstructGrid->GetX();
    RDouble *y = unstructGrid->GetY();
    RDouble *z = unstructGrid->GetZ();

    RDouble *cellCenterX = this->GetCellCenterX();
    RDouble *cellCenterY = this->GetCellCenterY();
    RDouble *cellCenterZ = this->GetCellCenterZ();

    RDouble *faceCenterX = this->GetFaceCenterX();
    RDouble *faceCenterY = this->GetFaceCenterY();
    RDouble *faceCenterZ = this->GetFaceCenterZ();

    RDouble *faceNormalX = this->GetFaceNormalX();
    RDouble *faceNormalY = this->GetFaceNormalY();
    RDouble *faceNormalZ = this->GetFaceNormalZ();
    RDouble *faceArea = this->GetFaceArea();
    RDouble *cellVolume = this->GetCellVolume();

    RDouble c34 = 3.0 / 4.0;

    PHSPACE::ZeroField(cellCenterX, numberOfTotalCells);
    PHSPACE::ZeroField(cellCenterY, numberOfTotalCells);
    PHSPACE::ZeroField(cellCenterZ, numberOfTotalCells);
    PHSPACE::ZeroField(cellVolume,  numberOfTotalCells);

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int leftCellIndex = leftCellIndexContainer[iFace];
        int rightCellIndex = rightCellIndexContainer[iFace];
        if (iFace < numberOfBoundaryFaces) rightCellIndex = iFace + numberOfCells;

        RDouble xc = faceCenterX[iFace];
        RDouble yc = faceCenterY[iFace];
        RDouble zc = faceCenterZ[iFace];

        RDouble vvol = third * (faceNormalX[iFace] * xc + faceNormalY[iFace] * yc + faceNormalZ[iFace] * zc) * faceArea[iFace];
        xc *= c34 * vvol;
        yc *= c34 * vvol;
        zc *= c34 * vvol;

        cellCenterX[leftCellIndex] += xc;
        cellCenterY[leftCellIndex] += yc;
        cellCenterZ[leftCellIndex] += zc;
        cellVolume[leftCellIndex] += vvol;

        cellCenterX[rightCellIndex] -= xc;
        cellCenterY[rightCellIndex] -= yc;
        cellCenterZ[rightCellIndex] -= zc;
        cellVolume[rightCellIndex] -= vvol;
    }

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        RDouble tmp = 1.0 / (cellVolume[iCell] + SMALL);
        cellCenterX[iCell] *= tmp;
        cellCenterY[iCell] *= tmp;
        cellCenterZ[iCell] *= tmp;
    }

    RDouble *txc = new RDouble[numberOfTotalCells];
    RDouble *tyc = new RDouble[numberOfTotalCells];
    RDouble *tzc = new RDouble[numberOfTotalCells];

    PHSPACE::SetField(txc, cellCenterX, numberOfCells);
    PHSPACE::SetField(tyc, cellCenterY, numberOfCells);
    PHSPACE::SetField(tzc, cellCenterZ, numberOfCells);

    PHSPACE::ZeroField(cellCenterX, numberOfTotalCells);
    PHSPACE::ZeroField(cellCenterY, numberOfTotalCells);
    PHSPACE::ZeroField(cellCenterZ, numberOfTotalCells);
    PHSPACE::ZeroField(cellVolume, numberOfTotalCells);

    int nodePosition = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int leftCellIndex = leftCellIndexContainer[iFace];
        int rightCellIndex = rightCellIndexContainer[iFace];
        if (iFace < numberOfBoundaryFaces) rightCellIndex = iFace + numberOfCells;

        int faceNodeNumber = faceNodeNumberContainer[iFace];
        for (int iNodeInFace = 0; iNodeInFace < faceNodeNumber; ++iNodeInFace)
        {
            int p2 = faceNodeIndexContainer[nodePosition + iNodeInFace];
            int p3 = faceNodeIndexContainer[nodePosition + (iNodeInFace + 1) % faceNodeNumber];

            RDouble x21 = x[p2] - faceCenterX[iFace];
            RDouble y21 = y[p2] - faceCenterY[iFace];
            RDouble z21 = z[p2] - faceCenterZ[iFace];

            RDouble x31 = x[p3] - faceCenterX[iFace];
            RDouble y31 = y[p3] - faceCenterY[iFace];
            RDouble z31 = z[p3] - faceCenterZ[iFace];

            RDouble nx = half * (y21 * z31 - y31 * z21);
            RDouble ny = half * (z21 * x31 - z31 * x21);
            RDouble nz = half * (x21 * y31 - x31 * y21);

            RDouble xc = third * (faceCenterX[iFace] + x[p2] + x[p3]);
            RDouble yc = third * (faceCenterY[iFace] + y[p2] + y[p3]);
            RDouble zc = third * (faceCenterZ[iFace] + z[p2] + z[p3]);

            RDouble hxcl = xc - txc[leftCellIndex];
            RDouble hycl = yc - tyc[leftCellIndex];
            RDouble hzcl = zc - tzc[leftCellIndex];

            RDouble hxcr = xc - txc[rightCellIndex];
            RDouble hycr = yc - tyc[rightCellIndex];
            RDouble hzcr = zc - tzc[rightCellIndex];

            // Cell Center and Volume
            RDouble vvoll = third * (nx * hxcl + ny * hycl + nz * hzcl);
            RDouble vvolr = third * (nx * hxcr + ny * hycr + nz * hzcr);

            RDouble xcl = xc - fourth * hxcl;
            RDouble ycl = yc - fourth * hycl;
            RDouble zcl = zc - fourth * hzcl;
            xcl *= vvoll;
            ycl *= vvoll;
            zcl *= vvoll;

            cellCenterX[leftCellIndex] += xcl;
            cellCenterY[leftCellIndex] += ycl;
            cellCenterZ[leftCellIndex] += zcl;
            cellVolume[leftCellIndex] += vvoll;

            RDouble xcr = xc - fourth * hxcr;
            RDouble ycr = yc - fourth * hycr;
            RDouble zcr = zc - fourth * hzcr;
            xcr *= vvolr;
            ycr *= vvolr;
            zcr *= vvolr;

            cellCenterX[rightCellIndex] -= xcr;
            cellCenterY[rightCellIndex] -= ycr;
            cellCenterZ[rightCellIndex] -= zcr;
            cellVolume[rightCellIndex] -= vvolr;
        }
        nodePosition += faceNodeNumber;
    }

    ostringstream oss;

    // Don't forget the coefficients
    int  numberOfCellsHaveNegativeVolumes = 0;
    RDouble minvol = LARGE, maxvol = 0.0;
    int indexMinv = 0, indexMaxv = 0;
    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        RDouble tmp = 1.0 / (cellVolume[iCell] + SMALL);
        cellCenterX[iCell] *= tmp;
        cellCenterY[iCell] *= tmp;
        cellCenterZ[iCell] *= tmp;

        if (minvol > cellVolume[iCell])
        {
            minvol = cellVolume[iCell];
            indexMinv = iCell;
        }
        if (maxvol < cellVolume[iCell])
        {
            maxvol = cellVolume[iCell];
            indexMaxv = iCell;
        }

        if (cellVolume[iCell] <= 0.0)
        {
            cellVolume[iCell] = -cellVolume[iCell];
            ++numberOfCellsHaveNegativeVolumes;
            oss << "  Warning: negative volume cell index " << iCell << ", (x, y, z) " << cellCenterX[iCell] << " " << cellCenterY[iCell] << " " << cellCenterZ[iCell] << " !\n";
        }
    }

    if (numberOfCellsHaveNegativeVolumes) oss << numberOfCellsHaveNegativeVolumes << " cells have negative vols \n";

    // For ghost cells
    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++iFace)
    {
        int leftCellIndex = leftCellIndexContainer[iFace];
        int rightCellIndex = rightCellIndexContainer[iFace];

        if (faceArea[iFace] > SMALL)
        {
            RDouble tmp = 2.0 * ((cellCenterX[leftCellIndex] - faceCenterX[iFace]) * faceNormalX[iFace] +
                (cellCenterY[leftCellIndex] - faceCenterY[iFace]) * faceNormalY[iFace] +
                (cellCenterZ[leftCellIndex] - faceCenterZ[iFace]) * faceNormalZ[iFace]);
            cellCenterX[rightCellIndex] = cellCenterX[leftCellIndex] - faceNormalX[iFace] * tmp;
            cellCenterY[rightCellIndex] = cellCenterY[leftCellIndex] - faceNormalY[iFace] * tmp;
            cellCenterZ[rightCellIndex] = cellCenterZ[leftCellIndex] - faceNormalZ[iFace] * tmp;
        }
        else
        {
            // Degenerated faces
            cellCenterX[rightCellIndex] = -cellCenterX[leftCellIndex] + 2.0 * faceCenterX[iFace];
            cellCenterY[rightCellIndex] = -cellCenterY[leftCellIndex] + 2.0 * faceCenterY[iFace];
            cellCenterZ[rightCellIndex] = -cellCenterZ[leftCellIndex] + 2.0 * faceCenterZ[iFace];
        }
        cellVolume[rightCellIndex] = cellVolume[leftCellIndex];
    }

    delete [] txc;    txc = NULL;
    delete [] tyc;    tyc = NULL;
    delete [] tzc;    tzc = NULL;

    ActionKey *actionKey = PHSPACE::GetCurrentActionKey();
    PHSPACE::StreamToActionKey(actionKey, oss);
}

void GridManager::BackUpOldGrid()
{
    SimpleGrid *oldGrid = this->unstructGrid->GetOldGrid();
    if (! oldGrid)
    {
        oldGrid = new SimpleGrid();
        this->unstructGrid->SetOldGrid(oldGrid);
    }
    *oldGrid = * this->unstructGrid;

    int numberOfFaces         = this->unstructGrid->GetNTotalFace();
    int numberOfCells         = this->unstructGrid->GetNTotalCell();
    int numberOfBoundaryFaces = this->unstructGrid->GetNBoundFace();
    int numberOfTotalCells    = numberOfCells + numberOfBoundaryFaces;

    if (IsNotAllocated(oldFaceCenterX)) oldFaceCenterX = new RDouble[numberOfFaces]();
    if (IsNotAllocated(oldFaceCenterY)) oldFaceCenterY = new RDouble[numberOfFaces]();
    if (IsNotAllocated(oldFaceCenterZ)) oldFaceCenterZ = new RDouble [numberOfFaces]();
    if (IsNotAllocated(oldFaceArea))  oldFaceArea  = new RDouble[numberOfFaces]();
    if (IsNotAllocated(cellVolumeN1)) cellVolumeN1 = new RDouble[numberOfTotalCells]();
    if (IsNotAllocated(cellVolumeN2)) cellVolumeN2 = new RDouble [numberOfTotalCells]();

    PHSPACE::SetField(oldFaceCenterX, unstructGrid->GetFaceCenterX(), numberOfFaces);
    PHSPACE::SetField(oldFaceCenterY, unstructGrid->GetFaceCenterY(), numberOfFaces);
    PHSPACE::SetField(oldFaceCenterZ, unstructGrid->GetFaceCenterZ(), numberOfFaces);
    PHSPACE::SetField(oldFaceArea, unstructGrid->GetFaceArea(), numberOfFaces);

    PHSPACE::SetField(cellVolumeN2, cellVolumeN1, numberOfTotalCells);
    PHSPACE::SetField(cellVolumeN1, unstructGrid->GetCellVolume(), numberOfTotalCells);

    PHSPACE::SetField(oldFaceNormalVelocity, newFaceNormalVelocity, numberOfFaces);
    PHSPACE::SetField(oldFaceVelocityX, newFaceVelocityX, numberOfFaces);
    PHSPACE::SetField(oldFaceVelocityY, newFaceVelocityY, numberOfFaces);
    PHSPACE::SetField(oldFaceVelocityZ, newFaceVelocityZ, numberOfFaces);
}

void GridManager::ComputeMassCharacterOfObject()
{
    //0- prepare the data.
    int numberOfCells   = unstructGrid->GetNTotalCell();
    RDouble *cellVolume = unstructGrid->GetCellVolume();

    //0.1- compute the mass of every cell.
    RDouble totalMass   = 10000.0;
    RDouble totalVolume = 0.0;
    RDouble *cellMass = new RDouble [numberOfCells];
    
    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        totalVolume += cellVolume[ iCell ];
    }

    ofstream Volume ("totalVolume.dat", ios_base::out|ios_base::app);
    Volume << totalVolume << "\n";
    Volume.close();
    Volume.clear();

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        cellMass[ iCell ] = totalMass / totalVolume * cellVolume[ iCell ];
    }

    //1-calculate the mass center coordinate.
    RDouble massCenterX = 0.0;
    RDouble massCenterY = 0.0;
    RDouble massCenterZ = 0.0;

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        //! minus the origin position later!
        RDouble x = cellCenterX[ iCell ] - 0.0;
        RDouble y = cellCenterY[ iCell ] - 0.0;
        RDouble z = cellCenterZ[ iCell ] - 0.0;
        
        massCenterX += x * cellMass[ iCell ];
        massCenterY += y * cellMass[ iCell ];
        massCenterZ += z * cellMass[ iCell ];
    }

    massCenterX /= totalMass;
    massCenterY /= totalMass;
    massCenterZ /= totalMass;

    ofstream massCenter ("massCenter.dat", ios_base::out|ios_base::app);
    massCenter << massCenterX << '\t' << massCenterY << '\t' << massCenterZ << "\n";
    massCenter.close();
    massCenter.clear();

    //2. calculate moment of inertia.
    RDouble iXX = 0.0;
    RDouble iYY = 0.0;
    RDouble iZZ = 0.0;
    RDouble iXY = 0.0;
    RDouble iYZ = 0.0;
    RDouble iZX = 0.0;

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        //! minus the mass cell center position later!
        RDouble x = cellCenterX[ iCell ] - 2.4880824;
        RDouble y = cellCenterY[ iCell ] - 0.0;
        RDouble z = cellCenterZ[ iCell ] - 0.0;
        
        iXX += cellMass[ iCell ] * (PHSPACE::SQR(y) + PHSPACE::SQR(z));
        iYY += cellMass[ iCell ] * (PHSPACE::SQR(z) + PHSPACE::SQR(x));
        iZZ += cellMass[ iCell ] * (PHSPACE::SQR(x) + PHSPACE::SQR(y));
        iXY += cellMass[ iCell ] * (x * y);
        iYZ += cellMass[ iCell ] * (y * z);
        iZX += cellMass[ iCell ] * (z * x);
    }

    ofstream massCharacter ("massCharacter.dat", ios_base::out|ios_base::app);
    massCharacter << iXX << '\t' << iYY << '\t' << iZZ << '\t' << iXY << '\t' << iYZ<< '\t' << iZX << "\n";
    massCharacter.close();
    massCharacter.clear();
    delete [] cellMass;    cellMass = NULL;
}

void GridManager::DumpCellCenterFile()
{
    fstream cellCenterFile;
    cellCenterFile.open("./grid/cellCenterFile.ccf", ios::out);

    int numberOfBoundaryFaces = unstructGrid->GetNBoundFace();
    int numberOfCells = unstructGrid->GetNTotalCell();
    int numberOfTotalCells = numberOfCells + numberOfBoundaryFaces;
    RDouble *cellCenterX = unstructGrid->GetCellCenterX();
    RDouble *cellCenterY = unstructGrid->GetCellCenterY();
    RDouble *cellCenterZ = unstructGrid->GetCellCenterZ();

    cellCenterFile << numberOfTotalCells << endl;

    for (int iCell = 0; iCell < numberOfTotalCells; ++ iCell)
    {
        cellCenterFile << cellCenterX[ iCell ] << " "
                       << cellCenterY[ iCell ] << " "
                       << cellCenterZ[ iCell ] << " "
                       << endl;
    }

    cellCenterFile.close();
    cellCenterFile.clear();
}

void GridManager::ComputeGridFaceVelocity()
{
    int strategyFaceNormalVelocity = GlobalDataBase::GetIntParaFromDB("strategyForFaceNormalVelocity");
    if (strategyFaceNormalVelocity == 0)
    {
        this->ComputeGridFaceVelocitySatisfyingGCL();
    }
    else if (strategyFaceNormalVelocity == 1)
    {
        this->ComputeGridFaceVelocityDirectly();
    }
    else if (strategyFaceNormalVelocity == 2)
    {
        this->ComputeGridFaceVelocityDirectly();
    }
    else
    {
        //StopProgram("ÎÞ´ËstrategyForFaceNormalVelocity!");
    }
}

void GridManager::ComputeGridFaceVelocityDirectly()
{
    RDouble *newFaceCenterX = unstructGrid->GetFaceCenterX();
    RDouble *newFaceCenterY = unstructGrid->GetFaceCenterY();
    RDouble *newFaceCenterZ = unstructGrid->GetFaceCenterZ();

    RDouble *oldFaceCenterX = this->GetOldFaceCenterX();
    RDouble *oldFaceCenterY = this->GetOldFaceCenterY();
    RDouble *oldFaceCenterZ = this->GetOldFaceCenterZ();

    int numberOfFaces = unstructGrid->GetNTotalFace();

    RDouble *newFaceNormalVelocity = this->GetNewFaceNormalVelocity();
    RDouble *newFaceVelocityX      = this->GetNewFaceVelocityX();
    RDouble *newFaceVelocityY      = this->GetNewFaceVelocityY();
    RDouble *newFaceVelocityZ      = this->GetNewFaceVelocityZ();

    RDouble *faceNormalVelocity = unstructGrid->GetFaceNormalVelocity();
    RDouble *faceVelocityX      = unstructGrid->GetFaceVelocityX();
    RDouble *faceVelocityY      = unstructGrid->GetFaceVelocityY();
    RDouble *faceVelocityZ      = unstructGrid->GetFaceVelocityZ();

    RDouble *oldFaceVelocityX = this->GetOldFaceVelocityX();
    RDouble *oldFaceVelocityY = this->GetOldFaceVelocityY();
    RDouble *oldFaceVelocityZ = this->GetOldFaceVelocityZ();

    RDouble *newFaceNormalX = unstructGrid->GetFaceNormalX();
    RDouble *newFaceNormalY = unstructGrid->GetFaceNormalY();
    RDouble *newFaceNormalZ = unstructGrid->GetFaceNormalZ();

    RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        RDouble newFcX = newFaceCenterX[iFace];
        RDouble newFcY = newFaceCenterY[iFace];
        RDouble newFcZ = newFaceCenterZ[iFace];
        RDouble oldFcX = oldFaceCenterX[iFace];
        RDouble oldFcY = oldFaceCenterY[iFace];
        RDouble oldFcZ = oldFaceCenterZ[iFace];
    
        newFaceVelocityX[iFace] = (newFcX - oldFcX) / physicalTimeStep;
        newFaceVelocityY[iFace] = (newFcY - oldFcY) / physicalTimeStep;
        newFaceVelocityZ[iFace] = (newFcZ - oldFcZ) / physicalTimeStep;

        newFaceNormalVelocity[iFace] = newFaceVelocityX[iFace] * newFaceNormalX[iFace]
                                     + newFaceVelocityY[iFace] * newFaceNormalY[iFace]
                                     + newFaceVelocityZ[iFace] * newFaceNormalZ[iFace];
    }

    RFloat a = 1.0;
    RFloat b = 0.0;

    int strategyFaceNormalVelocity = GlobalDataBase::GetIntParaFromDB("strategyForFaceNormalVelocity");
    if (strategyFaceNormalVelocity == 2)
    {
        a = 1.5;
        b = 0.5;
    }

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        //! faceVelocity is combined by the one order face velocities of the two times.
        faceVelocityX[ iFace ] = a * newFaceVelocityX[ iFace ] - b * oldFaceVelocityX[ iFace ];
        faceVelocityY[ iFace ] = a * newFaceVelocityY[ iFace ] - b * oldFaceVelocityY[ iFace ];
        faceVelocityZ[ iFace ] = a * newFaceVelocityZ[ iFace ] - b * oldFaceVelocityZ[ iFace ];

        faceNormalVelocity[ iFace ] = faceVelocityX[ iFace ] * newFaceNormalX[ iFace ]
                                    + faceVelocityY[ iFace ] * newFaceNormalY[ iFace ]
                                    + faceVelocityZ[ iFace ] * newFaceNormalZ[ iFace ];
    }
}

void GridManager::ComputeGridFaceVelocitySatisfyingGCL()
{
    int geometricDimension = GlobalDataBase::GetIntParaFromDB("ndim");

    if (geometricDimension == TWO_D)
    {
        this->ComputeGridFaceVelocity2DSatisfyingGCL();
    }
    else
    {
        this->ComputeGridFaceVelocity3DSatisfyingGCL();
    }
}

// this function needs to fit all the 2D situation in the arbitrary plane, not only in xy plane.
// the direction shoule be the same as others when computing the volume rate of change.
void GridManager::ComputeGridFaceVelocity2DSatisfyingGCL()
{
    SimpleGrid *oldGrid = unstructGrid->GetOldGrid();

    RDouble *x = unstructGrid->GetX();
    RDouble *y = unstructGrid->GetY();
    RDouble *z = unstructGrid->GetZ();

    RDouble *oldX = oldGrid->GetX();
    RDouble *oldY = oldGrid->GetY();
    RDouble *oldZ = oldGrid->GetZ();

    RDouble *faceCenterX = unstructGrid->GetFaceCenterX();
    RDouble *faceCenterY = unstructGrid->GetFaceCenterY();
    RDouble *faceCenterZ = unstructGrid->GetFaceCenterZ();

    RDouble *oldFaceCenterX = this->GetOldFaceCenterX();
    RDouble *oldFaceCenterY = this->GetOldFaceCenterY();
    RDouble *oldFaceCenterZ = this->GetOldFaceCenterZ();

    RDouble *faceArea    = unstructGrid->GetFaceArea();
    RDouble *oldFaceArea = this->GetOldFaceArea();

    int *faceNodeIndex = unstructGrid->GetFace2Node();
    int numberOfFaces  = unstructGrid->GetNTotalFace();

    RDouble *newFaceNormalVelocity = this->GetNewFaceNormalVelocity();
    RDouble *newFaceVelocityX      = this->GetNewFaceVelocityX();
    RDouble *newFaceVelocityY      = this->GetNewFaceVelocityY();
    RDouble *newFaceVelocityZ      = this->GetNewFaceVelocityZ();

    RDouble *faceNormalVelocity = unstructGrid->GetFaceNormalVelocity();
    RDouble *faceVelocityX      = unstructGrid->GetFaceVelocityX();
    RDouble *faceVelocityY      = unstructGrid->GetFaceVelocityY();
    RDouble *faceVelocityZ      = unstructGrid->GetFaceVelocityZ();

    RDouble *oldFaceNormalVelocity = this->GetOldFaceNormalVelocity();
    RDouble *oldFaceVelocityX      = this->GetOldFaceVelocityX();
    RDouble *oldFaceVelocityY      = this->GetOldFaceVelocityY();
    RDouble *oldFaceVelocityZ      = this->GetOldFaceVelocityZ();

    RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        RDouble darea  = 0.0;

        RDouble newFcX = faceCenterX[ iFace ];
        RDouble newFcY = faceCenterY[ iFace ];
        RDouble newFcZ = faceCenterZ[ iFace ];
        RDouble oldFcX = oldFaceCenterX[ iFace ];
        RDouble oldFcY = oldFaceCenterY[ iFace ];
        RDouble oldFcZ = oldFaceCenterZ[ iFace ];

        int p1 = faceNodeIndex[ 2 * iFace     ];
        int p2 = faceNodeIndex[ 2 * iFace + 1 ];
        
        RDouble oldX1 = oldX[ p1 ];
        RDouble oldY1 = oldY[ p1 ];
        RDouble oldZ1 = oldZ[ p1 ];

        RDouble oldX2 = oldX[ p2 ];
        RDouble oldY2 = oldY[ p2 ];
        RDouble oldZ2 = oldZ[ p2 ];

        RDouble newX1 = x[ p1 ];
        RDouble newY1 = y[ p1 ];
        RDouble newZ1 = z[ p1 ];

        RDouble newX2 = x[ p2 ];
        RDouble newY2 = y[ p2 ];
        RDouble newZ2 = z[ p2 ];

        RDouble dX1 = newX2 - oldX1;
        RDouble dY1 = newY2 - oldY1;
        RDouble dZ1 = newZ2 - oldZ1;
        
        RDouble dX2 = oldX2 - newX1;
        RDouble dY2 = oldY2 - newY1;
        RDouble dZ2 = oldZ2 - newZ1;
        //! the four points must in the same plane, so the area can be computed by the cross product of the diagonal.
        //! and it is right for arbitrary quadrilateral!
        darea = half * (dY1 * dZ2 - dZ1 * dY2 + dZ1 * dX2 - dX1 * dZ2 + dX1 * dY2 - dY1 * dX2);
    
        newFaceNormalVelocity[ iFace ] = darea / physicalTimeStep / (faceArea[iFace] + TINY);
        newFaceVelocityX[ iFace ] = (newFcX - oldFcX) / physicalTimeStep;
        newFaceVelocityY[ iFace ] = (newFcY - oldFcY) / physicalTimeStep;
        newFaceVelocityZ[ iFace ] = (newFcZ - oldFcZ) / physicalTimeStep;
    }

    RDouble a = 1.0;
    RDouble b = 0.0;

    int linearTwoStepMethods = PHSPACE::GlobalDataBase::GetIntParaFromDB("linearTwoStepMethods");

    if (linearTwoStepMethods == 1)
    {
        //linearTwoStepMethods == 1
        //Backward Euler: Order 1
        a = 1.0;
        b = 0.0;
    }
    else if (linearTwoStepMethods == 3)
    {
        //linearTwoStepMethods == 3
        //Backward differentiation: Order 2
        a = 1.5;
        b = 0.5;
    }

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        //! faceVelocity is combined by the one order face velocities of the two times(old and new).
        faceNormalVelocity[ iFace ] = a * newFaceNormalVelocity[ iFace ] - b * oldFaceNormalVelocity[ iFace ] * (oldFaceArea[ iFace ] / faceArea[ iFace ]);
        faceVelocityX     [ iFace ] = a * newFaceVelocityX     [ iFace ] - b * oldFaceVelocityX     [ iFace ];
        faceVelocityY     [ iFace ] = a * newFaceVelocityY     [ iFace ] - b * oldFaceVelocityY     [ iFace ];
        faceVelocityZ     [ iFace ] = a * newFaceVelocityZ     [ iFace ] - b * oldFaceVelocityZ     [ iFace ];
    }

    int ifLowSpeedPrecon = GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");

    if (ifLowSpeedPrecon == 1)
    {
        int numberOfNodes      = unstructGrid->GetNTotalNode();
        int numberOfCells      = unstructGrid->GetNTotalCell();
        RDouble *nodeVelocityX = new RDouble[numberOfNodes];
        RDouble *nodeVelocityY = new RDouble[numberOfNodes];
        RDouble *nodeVelocityZ = new RDouble[numberOfNodes];

        RDouble *cellVelocityX = unstructGrid->GetCellVelocityX();
        RDouble *cellVelocityY = unstructGrid->GetCellVelocityY();
        RDouble *cellVelocityZ = unstructGrid->GetCellVelocityZ();

        int *cell2node         = unstructGrid->GetCell2Node();
        int **cell2NodeArray   = unstructGrid->GetCell2NodeArray();
        int *nodeNumberOfEachCell = unstructGrid->GetNodeNumberOfEachCell();
        for (int iNode = 0; iNode < numberOfNodes; ++iNode)
        {
            nodeVelocityX[iNode] = (x[iNode] - oldX[iNode]) / physicalTimeStep;
            nodeVelocityY[iNode] = (y[iNode] - oldY[iNode]) / physicalTimeStep;
            nodeVelocityZ[iNode] = (z[iNode] - oldZ[iNode]) / physicalTimeStep;
        }

        for (int iCell = 0; iCell < numberOfCells; iCell++)
        {
            for (int jNode = 0; jNode < nodeNumberOfEachCell[iCell]; jNode++)
            {
                int nodeIndex = cell2NodeArray[iCell][jNode];

                cellVelocityX[iCell] += nodeVelocityX[nodeIndex];
                cellVelocityY[iCell] += nodeVelocityY[nodeIndex];
                cellVelocityZ[iCell] += nodeVelocityZ[nodeIndex];

            }

            cellVelocityX[iCell] /= nodeNumberOfEachCell[iCell];
            cellVelocityY[iCell] /= nodeNumberOfEachCell[iCell];
            cellVelocityZ[iCell] /= nodeNumberOfEachCell[iCell];
        }

        DelPointer(nodeVelocityX);
        DelPointer(nodeVelocityY);
        DelPointer(nodeVelocityZ);
    }
}

void GridManager::ComputeGridFaceVelocity3DSatisfyingGCL()
{
    SimpleGrid *oldGrid = unstructGrid->GetOldGrid();

    RDouble *x = unstructGrid->GetX();
    RDouble *y = unstructGrid->GetY();
    RDouble *z = unstructGrid->GetZ();

    RDouble *oldX = oldGrid->GetX();
    RDouble *oldY = oldGrid->GetY();
    RDouble *oldZ = oldGrid->GetZ();

    RDouble *faceCenterX = unstructGrid->GetFaceCenterX();
    RDouble *faceCenterY = unstructGrid->GetFaceCenterY();
    RDouble *faceCenterZ = unstructGrid->GetFaceCenterZ();

    RDouble *oldFaceCenterX = this->GetOldFaceCenterX();
    RDouble *oldFaceCenterY = this->GetOldFaceCenterY();
    RDouble *oldFaceCenterZ = this->GetOldFaceCenterZ();

    RDouble *faceArea    = unstructGrid->GetFaceArea();
    RDouble *oldFaceArea = this->GetOldFaceArea();

    int *faceNodeIndex  = unstructGrid->GetFace2Node();
    int *faceNodeNumber = unstructGrid->GetNodeNumberOfEachFace();

    int numberOfFaces = unstructGrid->GetNTotalFace();
    RDouble *newFaceNormalVelocity = this->GetNewFaceNormalVelocity();
    RDouble *newFaceVelocityX      = this->GetNewFaceVelocityX();
    RDouble *newFaceVelocityY      = this->GetNewFaceVelocityY();
    RDouble *newFaceVelocityZ      = this->GetNewFaceVelocityZ();

    RDouble *faceNormalVelocity = unstructGrid->GetFaceNormalVelocity();
    RDouble *faceVelocityX      = unstructGrid->GetFaceVelocityX();
    RDouble *faceVelocityY      = unstructGrid->GetFaceVelocityY();
    RDouble *faceVelocityZ      = unstructGrid->GetFaceVelocityZ();

    RDouble *oldFaceNormalVelocity = this->GetOldFaceNormalVelocity();
    RDouble *oldFaceVelocityX      = this->GetOldFaceVelocityX();
    RDouble *oldFaceVelocityY      = this->GetOldFaceVelocityY();
    RDouble *oldFaceVelocityZ      = this->GetOldFaceVelocityZ();

    RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

    int nodePosition = 0;

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        RDouble dvol = 0.0;
        int p0 = faceNodeIndex[ nodePosition ];
        int nodeNumberEachFace = faceNodeNumber[iFace];

        RDouble newFcX = faceCenterX[ iFace ];
        RDouble newFcY = faceCenterY[ iFace ];
        RDouble newFcZ = faceCenterZ[ iFace ];
        RDouble oldFcX = oldFaceCenterX[ iFace ];
        RDouble oldFcY = oldFaceCenterY[ iFace ];
        RDouble oldFcZ = oldFaceCenterZ[ iFace ];

        for (int iNodeInFace = 0; iNodeInFace < nodeNumberEachFace; ++ iNodeInFace)
        {
            int p1, p2;
            if (iNodeInFace == nodeNumberEachFace - 1)
            {
                p1 = faceNodeIndex[ nodePosition ++ ];
                p2 = p0;
            }
            else
            {
                p1 = faceNodeIndex[ nodePosition ++ ];
                p2 = faceNodeIndex[ nodePosition    ];
            }

            RDouble oldX1 = oldX[ p1 ];
            RDouble oldY1 = oldY[ p1 ];
            RDouble oldZ1 = oldZ[ p1 ];

            RDouble oldX2 = oldX[ p2 ];
            RDouble oldY2 = oldY[ p2 ];
            RDouble oldZ2 = oldZ[ p2 ];

            RDouble newX1 = x[ p1 ];
            RDouble newY1 = y[ p1 ];
            RDouble newZ1 = z[ p1 ];

            RDouble newX2 = x[ p2 ];
            RDouble newY2 = y[ p2 ];
            RDouble newZ2 = z[ p2 ];
            RDouble xm    = 0.25 * (oldX1 + oldX2 + newX1 + newX2);
            RDouble ym    = 0.25 * (oldY1 + oldY2 + newY1 + newY2);
            RDouble zm    = 0.25 * (oldZ1 + oldZ2 + newZ1 + newZ2);

            dvol += ComputeVolume4P(newFcX, newFcY, newFcZ, xm, ym, zm, oldX1, oldY1, oldZ1, oldX2, oldY2, oldZ2)
                  + ComputeVolume4P(newFcX, newFcY, newFcZ, xm, ym, zm, newX1, newY1, newZ1, oldX1, oldY1, oldZ1)
                  + ComputeVolume4P(newFcX, newFcY, newFcZ, xm, ym, zm, newX2, newY2, newZ2, newX1, newY1, newZ1)
                  + ComputeVolume4P(newFcX, newFcY, newFcZ, xm, ym, zm, oldX2, oldY2, oldZ2, newX2, newY2, newZ2)
                  + ComputeVolume4P(oldFcX, oldFcY, oldFcZ, newFcX, newFcY, newFcZ, oldX1, oldY1, oldZ1, oldX2, oldY2, oldZ2);
        }
        newFaceNormalVelocity[ iFace ] = dvol / physicalTimeStep / (faceArea[iFace] + TINY);
        newFaceVelocityX[ iFace ] = (newFcX - oldFcX) / physicalTimeStep;
        newFaceVelocityY[ iFace ] = (newFcY - oldFcY) / physicalTimeStep;
        newFaceVelocityZ[ iFace ] = (newFcZ - oldFcZ) / physicalTimeStep;
    }

    RDouble a = 1.0;
    RDouble b = 0.0;

    int linearTwoStepMethods = PHSPACE::GlobalDataBase::GetIntParaFromDB("linearTwoStepMethods");
    if (linearTwoStepMethods == 1)
    {
        //linearTwoStepMethods == 1
        //Backward Euler: Order 1
        a = 1.0;
        b = 0.0;
    }
    else if (linearTwoStepMethods == 3)
    {
        //linearTwoStepMethods == 3
        //Backward differentiation: Order 2
        a = 1.5;
        b = 0.5;
    }

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        //! faceVelocity is combined by the one order face velocities of the two times(old and new).
        faceNormalVelocity[ iFace ] = a * newFaceNormalVelocity[ iFace ] - b * oldFaceNormalVelocity[ iFace ] * (oldFaceArea[ iFace ] / (faceArea[ iFace ] + TINY));
        faceVelocityX     [ iFace ] = a * newFaceVelocityX     [ iFace ] - b * oldFaceVelocityX     [ iFace ];
        faceVelocityY     [ iFace ] = a * newFaceVelocityY     [ iFace ] - b * oldFaceVelocityY     [ iFace ];
        faceVelocityZ     [ iFace ] = a * newFaceVelocityZ     [ iFace ] - b * oldFaceVelocityZ     [ iFace ];
    }

    int ifLowSpeedPrecon =  GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");

    if (ifLowSpeedPrecon == 1)
    {
        int numberOfNodes         = unstructGrid->GetNTotalNode();
        int numberOfCells         = unstructGrid->GetNTotalCell();
        RDouble *nodeVelocityX = new RDouble[numberOfNodes];
        RDouble *nodeVelocityY = new RDouble[numberOfNodes];
        RDouble *nodeVelocityZ = new RDouble[numberOfNodes];

        RDouble *cellVelocityX = unstructGrid->GetCellVelocityX();
        RDouble *cellVelocityY = unstructGrid->GetCellVelocityY();
        RDouble *cellVelocityZ = unstructGrid->GetCellVelocityZ();

        int *cell2node = unstructGrid->GetCell2Node();
        int **cell2NodeArray      = unstructGrid->GetCell2NodeArray();
        int *nodeNumberOfEachCell = unstructGrid->GetNodeNumberOfEachCell();
        for (int iNode = 0; iNode < numberOfNodes; ++iNode)
        {
            nodeVelocityX[iNode] = (x[iNode] - oldX[iNode]) / physicalTimeStep;
            nodeVelocityY[iNode] = (y[iNode] - oldY[iNode]) / physicalTimeStep;
            nodeVelocityZ[iNode] = (z[iNode] - oldZ[iNode]) / physicalTimeStep;
        }

        for (int iCell = 0; iCell < numberOfCells; iCell++)
        {
            for (int jNode = 0; jNode < nodeNumberOfEachCell[iCell]; jNode++)
            {
                int nodeIndex = cell2NodeArray[iCell][jNode];

                cellVelocityX[iCell] += nodeVelocityX[nodeIndex];
                cellVelocityY[iCell] += nodeVelocityY[nodeIndex];
                cellVelocityZ[iCell] += nodeVelocityZ[nodeIndex];

            }

            cellVelocityX[iCell] /= nodeNumberOfEachCell[iCell];
            cellVelocityY[iCell] /= nodeNumberOfEachCell[iCell];
            cellVelocityZ[iCell] /= nodeNumberOfEachCell[iCell];
        }

        DelPointer(nodeVelocityX);
        DelPointer(nodeVelocityY);
        DelPointer(nodeVelocityZ);
    }
}

void GridManager::ComputeLocalToGlobalNodeIndexMappingDirectly()
{
    int numberOfNodes = unstructGrid->GetNTotalNode();

    int *localToGlobalNodeIndexMapping = this->GetLocalToGlobalNodeIndexMapping();

    localToGlobalNodeIndexMapping = new int[numberOfNodes]();

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        localToGlobalNodeIndexMapping[ iNode ] = iNode;
    }
}

void GridManager::SkipReadLayoutFile(fstream & file, int targetZone)
{
    PHSPACE::PHSkipRead< int >(file, 3);

    for (int iZone = 0; iZone < targetZone; ++ iZone)
    {
        int numberOfNodes;
        int numberOfFaces;
        int numberOfCells;

        PHSPACE::PHRead(file, numberOfNodes); 
        PHSPACE::PHRead(file, numberOfFaces);
        PHSPACE::PHRead(file, numberOfCells);

        PHSPACE::PHSkipRead< int >(file, numberOfNodes);
        PHSPACE::PHSkipRead< int >(file, numberOfFaces);
        PHSPACE::PHSkipRead< int >(file, numberOfCells);
    }
}

void GridManager::ComputeLocalToGlobalNodeIndexMappingByReadingLayoutFile(int targetZone)
{
    fstream infile;
    string originalGridLayoutFileName = GlobalDataBase::GetStrParaFromDB("originalGridLayoutFileName");

    PHSPACE::OpenFile(infile, originalGridLayoutFileName, ios_base::in|ios_base::binary);

    int numberOfNodes = unstructGrid->GetNTotalNode();
    int numberOfCells = unstructGrid->GetNTotalCell();
    int numberOfFaces = unstructGrid->GetNTotalFace();

    SkipReadLayoutFile(infile, targetZone);

    int numberOfNodesRead;
    int numberOfFacesRead;
    int numberOfCellsRead;

    PHSPACE::PHRead(infile, numberOfNodesRead);
    PHSPACE::PHRead(infile, numberOfFacesRead);
    PHSPACE::PHRead(infile, numberOfCellsRead);

    if ((numberOfNodesRead != numberOfNodes) ||
        (numberOfFacesRead != numberOfFaces) ||
        (numberOfCellsRead != numberOfCells))
    {
        TK_Exit::ExceptionExit("The number of nodes/faces/cells is not patch! \n");
    }

    int *localToGlobalNodeIndexMapping = this->GetLocalToGlobalNodeIndexMapping();
    if (IsNotAllocated(localToGlobalNodeIndexMapping)) localToGlobalNodeIndexMapping = new int [numberOfNodes]();

    PHSPACE::PHRead(infile, localToGlobalNodeIndexMapping, numberOfNodes);

    if (localToGlobalNodeIndexMapping != NULL)
    {
        delete [] localToGlobalNodeIndexMapping;    localToGlobalNodeIndexMapping = NULL;
    }
    PHSPACE::CloseFile(infile);

}

void GridManager::ReadNewCoordinate()
{
    using namespace PHMPI;
    int numberOfDymamicGridStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("numberOfDymamicGridStep");

    int iStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
    if (iStep > numberOfDymamicGridStep)
    {
        iStep = (iStep - 1) % (numberOfDymamicGridStep) + 1;
    }

    string existentMovingGridFileFolder = PHSPACE::GlobalDataBase::GetStrParaFromDB("existentMovingGridFileFolder");
    string newNodeFileName = PHSPACE::AggregateGeneralString(existentMovingGridFileFolder, "step", iStep, "/mmgrid.in");
    fstream infile;
    PHSPACE::OpenFile(infile, newNodeFileName, ios_base::in|ios_base::binary);

    int numberOfGlobalNodes;
    PHSPACE::PHRead(infile, numberOfGlobalNodes);

    int count = 3 * numberOfGlobalNodes;

    globalCoordinateNew = new RDouble[count];
    PHSPACE::PHRead(infile, globalCoordinateNew, count);

    PHSPACE::CloseFile(infile);
}

void GridManager::UpdateGridNode()
{
    SimpleGrid *oldGrid = this->unstructGrid->GetOldGrid();

    RDouble *newX = this->unstructGrid->GetX();
    RDouble *newY = this->unstructGrid->GetY();
    RDouble *newZ = this->unstructGrid->GetZ();

    RDouble *oldX = oldGrid->GetX();
    RDouble *oldY = oldGrid->GetY();
    RDouble *oldZ = oldGrid->GetZ();

    int numberOfNodes = this->unstructGrid->GetNTotalNode();
    int *localToGlobalNodeIndexMapping = this->GetLocalToGlobalNodeIndexMapping();

    RDouble maxDx = 0;
    RDouble maxDY = 0;
    RDouble maxDZ = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        int globalNodeIndex = localToGlobalNodeIndexMapping[ iNode ];

        newX[ iNode ] = globalCoordinateNew[ globalNodeIndex * 3     ];
        newY[ iNode ] = globalCoordinateNew[ globalNodeIndex * 3 + 1 ];
        newZ[ iNode ] = globalCoordinateNew[ globalNodeIndex * 3 + 2 ];

        maxDx = PHSPACE::MAX(maxDx, PHSPACE::ABS(newX[ iNode ] - oldX[ iNode ]));
        maxDY = PHSPACE::MAX(maxDY, PHSPACE::ABS(newY[ iNode ] - oldY[ iNode ]));
        maxDZ = PHSPACE::MAX(maxDZ, PHSPACE::ABS(newZ[ iNode ] - oldZ[ iNode ]));
    }

    cout << "<oldGrid, newGrid> max dx, dy, dz = " << maxDx << "	" << maxDY << "	" << maxDZ << endl;
}

void GridManager::UpdateGridCoordinateByReadingFile(int targetZone)
{
    int innerIterationSteps = GlobalDataBase::GetIntParaFromDB("innstep");
    if (innerIterationSteps != 1) return;

    if (PHMPI::GetNumberofGlobalZones() == 1)
    {
        this->ComputeLocalToGlobalNodeIndexMappingDirectly();
    }
    else
    {
        this->ComputeLocalToGlobalNodeIndexMappingByReadingLayoutFile(targetZone);
    }

    this->ReadNewCoordinate();

    this->UpdateGridNode();
}

void GridManager::InitializeMovingGrids()
{
    this->InitializeOldFaceNormalVelocity();
    this->InitializeOldCellVolume();
}

void GridManager::InitializeOldFaceNormalVelocity()
{
    int numberOfFaces = this->unstructGrid->GetNTotalFace();

    GridManager *gridManager = this->unstructGrid->GetGridManager();

    RDouble *faceNormalVelocity    = gridManager-> GetFaceNormalVelocity();
    RDouble *oldFaceNormalVelocity = gridManager-> GetOldFaceNormalVelocity();

    PHSPACE::SetField(oldFaceNormalVelocity, faceNormalVelocity, numberOfFaces);
}

void GridManager::InitializeOldCellVolume()
{
    int numberOfCells = this->unstructGrid->GetNTotalCell();
    int numberOfBoundaryFaces = this->unstructGrid->GetNBoundFace();
    int numberOfTotalCells = numberOfCells + numberOfBoundaryFaces;

    GridManager * gridManager = this->unstructGrid->GetGridManager();

    RDouble *cellVolume   = gridManager->GetCellVolume();
    RDouble *cellVolumeN1 = gridManager->GetCellVolume(1);
    RDouble *cellVolumeN2 = gridManager->GetCellVolume(2);

    PHSPACE::SetField(cellVolumeN2, cellVolume, numberOfTotalCells);
    PHSPACE::SetField(cellVolumeN1, cellVolume, numberOfTotalCells);
}

void GridManager::ReadMovingGrid()
{
    string fileName = PHSPACE::AggregateGeneralString("moveGrid", this->unstructGrid->GetZoneID(), ".dat");

    RDouble *x = this->unstructGrid->GetX();
    RDouble *y = this->unstructGrid->GetY();
    RDouble *z = this->unstructGrid->GetZ();

    RDouble *oldCellVolumeN1Container = this->GetOldCellVolume();

    int numberOfNodes = unstructGrid->GetNTotalNode();
    int numberOfCells = this->unstructGrid->GetNTotalCell();
    int numberOfBoundaryFaces = this->unstructGrid->GetNBoundFace();
    int numberOfTotalCells = numberOfCells + numberOfBoundaryFaces;
    fstream file;
    OpenFile(file, fileName, ios_base::in|ios_base::binary);

    PHMPI::PH_Read(file, & x[ 0 ], numberOfNodes);
    PHMPI::PH_Read(file, & y[ 0 ], numberOfNodes);
    PHMPI::PH_Read(file, & z[ 0 ], numberOfNodes);
    PHMPI::PH_Read(file, oldCellVolumeN1Container, numberOfTotalCells);

    CloseFile(file);
}

void GridManager::UpdateOldCellVolume()
{  
    int numberOfCells = this->unstructGrid->GetNTotalCell();
    int numberOfBoundaryFaces = this->unstructGrid->GetNBoundFace();
    int numberOfTotalCells = numberOfCells + numberOfBoundaryFaces;

    GridManager *gridManager = this->unstructGrid->GetGridManager();

    RDouble *cellVolume = gridManager->GetCellVolume();
    RDouble *cellVolumeN1 = gridManager->GetCellVolume(1);
    RDouble *cellVolumeN2 = gridManager->GetCellVolume(2);

    PHSPACE::SetField(cellVolumeN2, cellVolumeN1, numberOfTotalCells);
    PHSPACE::SetField(cellVolumeN1, cellVolume  , numberOfTotalCells);
}

void GridManager::ConvergenceNorm(FieldProxy *q1Proxy, FieldProxy *q0Proxy, int numberOfEquations, RFloat &norm)
{
    int numberOfCells = this->unstructGrid->GetNTotalCell();

    RDouble **q1 = q1Proxy->GetField_UNS();
    RDouble **q0 = q0Proxy->GetField_UNS();

    norm = zero;
    for (int iEquation = 0; iEquation < numberOfEquations; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfCells; ++ iCell)
        {
            RFloat diff = q1[ iEquation ][ iCell ] - q0[ iEquation ][ iCell ];
            norm += PHSPACE::SQR(diff);
        }
    }
    
    norm = sqrt(norm / (numberOfCells * numberOfEquations));
}

RDouble ComputeVolume4P(RDouble x1, RDouble y1, RDouble z1, RDouble x2, RDouble y2, RDouble z2, RDouble x3, RDouble y3, RDouble z3, RDouble x4, RDouble y4, RDouble z4)
{
    RDouble xx1, yy1, zz1, xx2, yy2, zz2, xx3, yy3, zz3, vol;

    xx1 = x2 - x1;
    yy1 = y2 - y1;
    zz1 = z2 - z1;
    xx2 = x3 - x1;
    yy2 = y3 - y1;
    zz2 = z3 - z1;
    xx3 = x4 - x1;
    yy3 = y4 - y1;
    zz3 = z4 - z1;

    vol = xx1 * yy2 * zz3 + xx2 * yy3 * zz1 + xx3 * yy1 * zz2 - xx3 * yy2 * zz1 - xx2 * yy1 * zz3 - xx1 * yy3 * zz2;
    vol /= 6.0;

    return vol;
}

bool IsNotAllocated(void *pointer)
{
    if (pointer) return false;
    return true;
}

}

