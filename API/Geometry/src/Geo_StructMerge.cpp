#include <iostream>
#include <Geo_StructMerge.h>
#include <iostream>
#include "Mesh_Deformation.h"
#include "Mesh_DeformationSPRING.h"
#include "Mesh_DeformationRBF.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "GeometryUnit.h"
#include "Geo_Grid.h"
#include "PHMpi.h"
#include "Glb_Dimension.h"
#include "Pre_HDF5File.h"
#include <algorithm> 
#include <string.h>
#include <Geo_SimpleBC.h>
#include "GridType.h"
#include <Region.h>
#include <Post_WriteTecplot.h>
using namespace std;

namespace PHSPACE
{
using namespace PHMPI;

//! Certify the type of the two Informationofboundary or the two StructBC is the same.
int WhetherInformationofBoundarytheSame(Informationofboundary s1, Informationofboundary s2)
{
    int    bcType1, bcType2;
    string bodyName1, bodyName2, bcName1, bcName2;
    bcType1   = s1.bcType;
    bodyName1 = s1.bodyName;
    bcName1   = s1.boundaryName;
    bcType2   = s2.bcType;
    bodyName2 = s2.bodyName;
    bcName2   = s2.boundaryName;

    if ((bcType1   == bcType2)
     && (bodyName1 == bodyName2)
     && (bcName1   == bcName2))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int WhetherInformationofBoundarytheSame(StructBC *s1, StructBC *s2)
{
    int    bcType1, bcType2;
    string bodyName1, bodyName2,bcName1, bcName2;
    bcType1   = s1->GetBCType();
    bodyName1 = s1->GetBodyName();
    bcName1   = s1->GetBCName();
    bcType2   = s2->GetBCType();
    bodyName2 = s2->GetBodyName();
    bcName2   = s2->GetBCName();

    if ((bcType1   == bcType2)
     && (bodyName1 == bodyName2)
     && (bcName1   == bcName2))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

//! WhethercanMergeI/J/KDirect is used to certify whether the boundaries whose face direction is i/j/k direction can be merged.
//! MergeBCI/J/K is used to merge all boundaries whose face direction is i/j/k direction.
//! the effective length of merged boundaries is stored in *nBC.
int WhethercanMergeIDirect(vector < StructBC* > bc, int nBC)
{
    if (nBC > 1)
    {
        int iBCOrignalJst, iBCOrignalJed, iBCOrignalKst, iBCOrignalKed,
            jBCOrignalJst, jBCOrignalJed, jBCOrignalKst, jBCOrignalKed;

        for (int iBC = 1; iBC < nBC; ++ iBC)
        {
            for (int jBC = 0; jBC < iBC; ++ jBC)
            {
                bc[iBC]->ComputeOriginal();
                bc[jBC]->ComputeOriginal();
                iBCOrignalJst = bc[iBC]->GetStartOrignal(1);
                iBCOrignalKst = bc[iBC]->GetStartOrignal(2);
                iBCOrignalJed = bc[iBC]->GetEndOrignal(1);
                iBCOrignalKed = bc[iBC]->GetEndOrignal(2);
                jBCOrignalJst = bc[jBC]->GetStartOrignal(1);
                jBCOrignalKst = bc[jBC]->GetStartOrignal(2);
                jBCOrignalJed = bc[jBC]->GetEndOrignal(1);
                jBCOrignalKed = bc[jBC]->GetEndOrignal(2);

                if (((iBCOrignalJst == jBCOrignalJed) && ((iBCOrignalKst == jBCOrignalKst) && (iBCOrignalKed == jBCOrignalKed)))
                 || ((iBCOrignalJed == jBCOrignalJst) && ((iBCOrignalKst == jBCOrignalKst) && (iBCOrignalKed == jBCOrignalKed)))
                 || ((iBCOrignalKst == jBCOrignalKed) && ((iBCOrignalJst == jBCOrignalJst) && (iBCOrignalJed == jBCOrignalJed)))
                 || ((iBCOrignalKed == jBCOrignalKst) && ((iBCOrignalJst == jBCOrignalJst) && (iBCOrignalJed == jBCOrignalJed))))
                {
                    return 1;
                }
            }
        }
    }
    else
    {
        return 0;
    }
}

void MergeBCI(vector < StructBC* > &bc, int *nBC)
{
    if (WhethercanMergeIDirect(bc, *nBC))
    {
        int iBCOrignalJst, iBCOrignalJed, iBCOrignalKst, iBCOrignalKed,
            jBCOrignalJst, jBCOrignalJed, jBCOrignalKst, jBCOrignalKed;
        int iBCFind = 1, jBCFind = 0;
        StructGrid *jBCGrid;
        int jBCGridIndex;
        int whethertheSame = 0;

        for (int iBC = 1; iBC < *nBC; ++ iBC)
        {
            for (int jBC = 0; jBC < iBC; ++ jBC)
            {
                bc[iBC]->ComputeOriginal();
                bc[jBC]->ComputeOriginal();
                iBCOrignalJst = bc[iBC]->GetStartOrignal(1);
                iBCOrignalKst = bc[iBC]->GetStartOrignal(2);
                iBCOrignalJed = bc[iBC]->GetEndOrignal(1);
                iBCOrignalKed = bc[iBC]->GetEndOrignal(2);
                jBCOrignalJst = bc[jBC]->GetStartOrignal(1);
                jBCOrignalKst = bc[jBC]->GetStartOrignal(2);
                jBCOrignalJed = bc[jBC]->GetEndOrignal(1);
                jBCOrignalKed = bc[jBC]->GetEndOrignal(2);

                if (((iBCOrignalJst == jBCOrignalJed) && ((iBCOrignalKst == jBCOrignalKst) && (iBCOrignalKed == jBCOrignalKed)))
                 || ((iBCOrignalJed == jBCOrignalJst) && ((iBCOrignalKst == jBCOrignalKst) && (iBCOrignalKed == jBCOrignalKed)))
                 || ((iBCOrignalKst == jBCOrignalKed) && ((iBCOrignalJst == jBCOrignalJst) && (iBCOrignalJed == jBCOrignalJed)))
                 || ((iBCOrignalKed == jBCOrignalKst) && ((iBCOrignalJst == jBCOrignalJst) && (iBCOrignalJed == jBCOrignalJed))))
                {
                    whethertheSame = 1;
                    iBCFind        = iBC;
                    jBCFind        = jBC;
                    break;
                }
            }
            if (1 == whethertheSame)
            {
                break;
            }
        }

        int    iBCType   = bc[jBCFind]->GetBCType();
        string iBCName   = bc[jBCFind]->GetBCName();
        string iBodyName = bc[jBCFind]->GetBodyName();
        int bcDirectionNew, targetBCDirectionNew;
        int leftorRightNew, targetLeftorRightNew;
        int currentZoneIndex, targetCurrentZoneIndex;
        int stOrignalNew[3], edOrignalNew[3];
        bcDirectionNew       = bc[jBCFind]->GetFaceDirection();
        targetBCDirectionNew = bc[jBCFind]->GetFaceDirectionOfTargetBlock();
        leftorRightNew       = bc[jBCFind]->GetFaceLeftOrRightIndex();
        targetLeftorRightNew = bc[jBCFind]->GetFaceLeftOrRightIndexOfTargetBlock();
        if (MERGED != (bc[jBCFind]->GetWhetherMerge()))
        {
            jBCGridIndex           = bc[jBCFind]->GetRegionBlock();
            jBCGrid                = StructGridCast(GetGrid(jBCGridIndex, 0));
            currentZoneIndex       = jBCGrid->GetOrdinaryGridIndex();
            jBCGridIndex           = bc[jBCFind]->GetTargetRegionBlock();
            jBCGrid                = StructGridCast(GetGrid(jBCGridIndex, 0));
            targetCurrentZoneIndex = jBCGrid->GetOrdinaryGridIndex();
        }
        else
        {
            currentZoneIndex       = bc[jBCFind]->GetRegionBlock();
            targetCurrentZoneIndex = bc[jBCFind]->GetTargetRegionBlock();
        }

        if ((iBCOrignalJst == jBCOrignalJed) 
         && (iBCOrignalKst == jBCOrignalKst) 
         && (iBCOrignalKed == jBCOrignalKed))
        {
            stOrignalNew[0] = bc[jBCFind]->GetStartOrignal(0);
            stOrignalNew[1] = min(jBCOrignalJst, iBCOrignalJst);
            stOrignalNew[2] = jBCOrignalKst;
            edOrignalNew[0] = bc[jBCFind]->GetEndOrignal(0);
            edOrignalNew[1] = max(iBCOrignalJed, jBCOrignalJed);
            edOrignalNew[2] = jBCOrignalKed;
        }
        else if ((iBCOrignalJed == jBCOrignalJst) 
              && (iBCOrignalKst == jBCOrignalKst) 
              && (iBCOrignalKed == jBCOrignalKed))
        {
            stOrignalNew[0] = bc[jBCFind]->GetStartOrignal(0);
            stOrignalNew[1] = min(jBCOrignalJst, iBCOrignalJst);
            stOrignalNew[2] = jBCOrignalKst;
            edOrignalNew[0] = bc[jBCFind]->GetEndOrignal(0);
            edOrignalNew[1] = max(iBCOrignalJed, jBCOrignalJed);
            edOrignalNew[2] = jBCOrignalKed;
        }
        else if ((iBCOrignalKst == jBCOrignalKed) 
              && (iBCOrignalJst == jBCOrignalJst) 
              && (iBCOrignalJed == jBCOrignalJed))
        {
            stOrignalNew[0] = bc[jBCFind]->GetStartOrignal(0);
            stOrignalNew[1] = jBCOrignalJst;
            stOrignalNew[2] = min(jBCOrignalKst, iBCOrignalKst);
            edOrignalNew[0] = bc[jBCFind]->GetEndOrignal(0);
            edOrignalNew[1] = jBCOrignalJed;
            edOrignalNew[2] = max(iBCOrignalKed, jBCOrignalKed);
        }
        else if ((iBCOrignalKed == jBCOrignalKst) 
              && (iBCOrignalJst == jBCOrignalJst) 
              && (iBCOrignalJed == jBCOrignalJed))
        {
            stOrignalNew[0] = bc[jBCFind]->GetStartOrignal(0);
            stOrignalNew[1] = jBCOrignalJst;
            stOrignalNew[2] = min(jBCOrignalKst, iBCOrignalKst);
            edOrignalNew[0] = bc[jBCFind]->GetEndOrignal(0);
            edOrignalNew[1] = jBCOrignalJed;
            edOrignalNew[2] = max(iBCOrignalKed, jBCOrignalKed);
        }

        bc[jBCFind]->SetCurrentZoneIndex(currentZoneIndex);
        bc[jBCFind]->SetFaceDirection(bcDirectionNew);
        bc[jBCFind]->SetFaceDirectionOfTargetBlock(targetBCDirectionNew);
        bc[jBCFind]->SetFaceLeftOrRightIndex(leftorRightNew);
        bc[jBCFind]->SetFaceLeftOrRightIndexOfTargetBlock(targetLeftorRightNew);
        bc[jBCFind]->SetTargetRegionBlock(targetCurrentZoneIndex);
        bc[jBCFind]->SetStartOrignal(stOrignalNew);
        bc[jBCFind]->SetEndOrignal(edOrignalNew);
        bc[jBCFind]->SetWhetherMerge(MERGED);
        bc[jBCFind]->SetBCType(iBCType);
        bc[jBCFind]->SetBCName(iBCName);
        bc[jBCFind]->SetBodyName(iBodyName);

        for (int mBC = iBCFind; mBC < (*nBC - 1); ++ mBC)
        {
            bc[mBC]->SetBCType(bc[mBC + 1]->GetBCType());
            bc[mBC]->SetBCName(bc[mBC + 1]->GetBCName());
            bc[mBC]->SetBodyName(bc[mBC + 1]->GetBodyName());
            bc[mBC] = bc[mBC + 1];
        }
        *nBC = *nBC - 1;
        MergeBCI(bc, nBC);
    }
    else
    {
        if (MERGED != (bc[0]->GetWhetherMerge()))
        {
            int currentZoneIndex, targetCurrentZoneIndex;
            int bcGridIndex;
            StructGrid *bcGrid;
            bcGridIndex            = bc[0]->GetRegionBlock();
            bcGrid                 = StructGridCast(GetGrid(bcGridIndex, 0));
            currentZoneIndex       = bcGrid->GetOrdinaryGridIndex();
            bcGridIndex            = bc[0]->GetTargetRegionBlock();
            bcGrid                 = StructGridCast(GetGrid(bcGridIndex, 0));
            targetCurrentZoneIndex = bcGrid->GetOrdinaryGridIndex();
            bc[0]->SetCurrentZoneIndex(currentZoneIndex);
            bc[0]->SetTargetRegionBlock(targetCurrentZoneIndex);
            bc[0]->SetWhetherMerge(MERGED);
        }
    }
}

int WhethercanMergeJDirect(vector < StructBC* > bc, int nBC)
{
    if (nBC > 1)
    {
        int iBCOrignalIst, iBCOrignalIed, iBCOrignalKst, iBCOrignalKed,
            jBCOrignalIst, jBCOrignalIed, jBCOrignalKst, jBCOrignalKed;

        for (int iBC = 1; iBC < nBC; ++ iBC)
        {
            for (int jBC = 0; jBC < iBC; ++ jBC)
            {
                bc[iBC]->ComputeOriginal();
                bc[jBC]->ComputeOriginal();
                iBCOrignalIst = bc[iBC]->GetStartOrignal(0);
                iBCOrignalKst = bc[iBC]->GetStartOrignal(2);
                iBCOrignalIed = bc[iBC]->GetEndOrignal(0);
                iBCOrignalKed = bc[iBC]->GetEndOrignal(2);
                jBCOrignalIst = bc[jBC]->GetStartOrignal(0);
                jBCOrignalKst = bc[jBC]->GetStartOrignal(2);
                jBCOrignalIed = bc[jBC]->GetEndOrignal(0);
                jBCOrignalKed = bc[jBC]->GetEndOrignal(2);

                if (((iBCOrignalIst == jBCOrignalIed) && ((iBCOrignalKst == jBCOrignalKst) && (iBCOrignalKed == jBCOrignalKed)))
                 || ((iBCOrignalIed == jBCOrignalIst) && ((iBCOrignalKst == jBCOrignalKst) && (iBCOrignalKed == jBCOrignalKed)))
                 || ((iBCOrignalKst == jBCOrignalKed) && ((iBCOrignalIst == jBCOrignalIst) && (iBCOrignalIed == jBCOrignalIed)))
                 || ((iBCOrignalKed == jBCOrignalKst) && ((iBCOrignalIst == jBCOrignalIst) && (iBCOrignalIed == jBCOrignalIed))))
                {
                    return 1;
                }
            }
        }
    }
    else
    {
        return 0;
    }
}

void MergeBCJ(vector < StructBC* > &bc, int *nBC)
{
    if (WhethercanMergeJDirect(bc, *nBC))
    {
        int iBCOrignalIst, iBCOrignalIed, iBCOrignalKst, iBCOrignalKed,
            jBCOrignalIst, jBCOrignalIed, jBCOrignalKst, jBCOrignalKed;
        int iBCFind = 1, jBCFind = 0;
        StructGrid *jBCGrid;
        int jBCGridIndex;
        int whethertheSame = 0;

        for (int iBC = 1; iBC < *nBC; ++ iBC)
        {
            for (int jBC = 0; jBC < iBC; ++ jBC)
            {
                bc[iBC]->ComputeOriginal();
                bc[jBC]->ComputeOriginal();
                iBCOrignalIst = bc[iBC]->GetStartOrignal(0);
                iBCOrignalKst = bc[iBC]->GetStartOrignal(2);
                iBCOrignalIed = bc[iBC]->GetEndOrignal(0);
                iBCOrignalKed = bc[iBC]->GetEndOrignal(2);
                jBCOrignalIst = bc[jBC]->GetStartOrignal(0);
                jBCOrignalKst = bc[jBC]->GetStartOrignal(2);
                jBCOrignalIed = bc[jBC]->GetEndOrignal(0);
                jBCOrignalKed = bc[jBC]->GetEndOrignal(2);

                if (((iBCOrignalIst == jBCOrignalIed) && ((iBCOrignalKst == jBCOrignalKst) && (iBCOrignalKed == jBCOrignalKed)))
                 || ((iBCOrignalIed == jBCOrignalIst) && ((iBCOrignalKst == jBCOrignalKst) && (iBCOrignalKed == jBCOrignalKed)))
                 || ((iBCOrignalKst == jBCOrignalKed) && ((iBCOrignalIst == jBCOrignalIst) && (iBCOrignalIed == jBCOrignalIed)))
                 || ((iBCOrignalKed == jBCOrignalKst) && ((iBCOrignalIst == jBCOrignalIst) && (iBCOrignalIed == jBCOrignalIed))))
                {
                    whethertheSame = 1;
                    iBCFind        = iBC;
                    jBCFind        = jBC;
                    break;
                }

            }
            if (1 == whethertheSame)
            {
                break;
            }
        }

        int    iBCType   = bc[jBCFind]->GetBCType();
        string iBCName   = bc[jBCFind]->GetBCName();
        string iBodyName = bc[jBCFind]->GetBodyName();
        int bcDirectionNew, targetBCDirectionNew;
        int leftorRightNew, targetLeftorRightNew;
        int currentZoneIndex, targetCurrentZoneIndex;
        int stOrignalNew[3], edOrignalNew[3];
        bcDirectionNew       = bc[jBCFind]->GetFaceDirection();
        targetBCDirectionNew = bc[jBCFind]->GetFaceDirectionOfTargetBlock();
        leftorRightNew       = bc[jBCFind]->GetFaceLeftOrRightIndex();
        targetLeftorRightNew = bc[jBCFind]->GetFaceLeftOrRightIndexOfTargetBlock();
        if ((bc[jBCFind]->GetWhetherMerge()) != MERGED)
        {
            jBCGridIndex           = bc[jBCFind]->GetRegionBlock();
            jBCGrid                = StructGridCast(GetGrid(jBCGridIndex, 0));
            currentZoneIndex       = jBCGrid->GetOrdinaryGridIndex();
            jBCGridIndex           = bc[jBCFind]->GetTargetRegionBlock();
            jBCGrid                = StructGridCast(GetGrid(jBCGridIndex, 0));
            targetCurrentZoneIndex = jBCGrid->GetOrdinaryGridIndex();
        }
        else
        {
            currentZoneIndex       = bc[jBCFind]->GetRegionBlock();
            targetCurrentZoneIndex = bc[jBCFind]->GetTargetRegionBlock();
        }

        if ((iBCOrignalIst == jBCOrignalIed)
         && (iBCOrignalKst == jBCOrignalKst)
         && (iBCOrignalKed == jBCOrignalKed))
        {
            stOrignalNew[1] = bc[jBCFind]->GetStartOrignal(1);
            stOrignalNew[0] = min(jBCOrignalIst, iBCOrignalIst);
            stOrignalNew[2] = jBCOrignalKst;
            edOrignalNew[1] = bc[jBCFind]->GetEndOrignal(1);
            edOrignalNew[0] = max(iBCOrignalIed, jBCOrignalIed);
            edOrignalNew[2] = jBCOrignalKed;
        }
        else if ((iBCOrignalIed == jBCOrignalIst)
              && (iBCOrignalKst == jBCOrignalKst)
              && (iBCOrignalKed == jBCOrignalKed))
        {
            stOrignalNew[1] = bc[jBCFind]->GetStartOrignal(1);
            stOrignalNew[0] = min(jBCOrignalIst, iBCOrignalIst);
            stOrignalNew[2] = jBCOrignalKst;
            edOrignalNew[1] = bc[jBCFind]->GetEndOrignal(1);
            edOrignalNew[0] = max(iBCOrignalIed, jBCOrignalIed);
            edOrignalNew[2] = jBCOrignalKed;
        }
        else if ((iBCOrignalKst == jBCOrignalKed)
              && (iBCOrignalIst == jBCOrignalIst)
              && (iBCOrignalIed == jBCOrignalIed))
        {
            stOrignalNew[1] = bc[jBCFind]->GetStartOrignal(1);
            stOrignalNew[0] = jBCOrignalIst;
            stOrignalNew[2] = min(jBCOrignalKst, iBCOrignalKst);
            edOrignalNew[1] = bc[jBCFind]->GetEndOrignal(1);
            edOrignalNew[0] = jBCOrignalIed;
            edOrignalNew[2] = max(iBCOrignalKed, jBCOrignalKed);
        }
        else if ((iBCOrignalKed == jBCOrignalKst)
              && (iBCOrignalIst == jBCOrignalIst)
              && (iBCOrignalIed == jBCOrignalIed))
        {
            stOrignalNew[1] = bc[jBCFind]->GetStartOrignal(1);
            stOrignalNew[0] = jBCOrignalIst;
            stOrignalNew[2] = min(jBCOrignalKst, iBCOrignalKst);
            edOrignalNew[1] = bc[jBCFind]->GetEndOrignal(1);
            edOrignalNew[0] = jBCOrignalIed;
            edOrignalNew[2] = max(iBCOrignalKed, jBCOrignalKed);
        }

        bc[jBCFind]->SetCurrentZoneIndex(currentZoneIndex);
        bc[jBCFind]->SetFaceDirection(bcDirectionNew);
        bc[jBCFind]->SetFaceDirectionOfTargetBlock(targetBCDirectionNew);
        bc[jBCFind]->SetFaceLeftOrRightIndex(leftorRightNew);
        bc[jBCFind]->SetFaceLeftOrRightIndexOfTargetBlock(targetLeftorRightNew);
        bc[jBCFind]->SetTargetRegionBlock(targetCurrentZoneIndex);
        bc[jBCFind]->SetStartOrignal(stOrignalNew);
        bc[jBCFind]->SetEndOrignal(edOrignalNew);
        bc[jBCFind]->SetWhetherMerge(MERGED);
        bc[jBCFind]->SetBCType(iBCType);
        bc[jBCFind]->SetBCName(iBCName);
        bc[jBCFind]->SetBodyName(iBodyName);

        for (int mBC = iBCFind; mBC < (*nBC - 1); ++ mBC)
        {
            bc[mBC]->SetBCType(bc[mBC + 1]->GetBCType());
            bc[mBC]->SetBCName(bc[mBC + 1]->GetBCName());
            bc[mBC]->SetBodyName(bc[mBC + 1]->GetBodyName());
            bc[mBC] = bc[mBC + 1];
        }
        *nBC = *nBC - 1;
        MergeBCJ(bc, nBC);
    }
    else
    {
        if ((bc[0]->GetWhetherMerge()) != MERGED)
        {
            int currentZoneIndex, targetCurrentZoneIndex;
            int bcGridIndex;
            StructGrid *bcGrid;
            bcGridIndex            = bc[0]->GetRegionBlock();
            bcGrid                 = StructGridCast(GetGrid(bcGridIndex, 0));
            currentZoneIndex       = bcGrid->GetOrdinaryGridIndex();
            bcGridIndex            = bc[0]->GetTargetRegionBlock();
            bcGrid                 = StructGridCast(GetGrid(bcGridIndex, 0));
            targetCurrentZoneIndex = bcGrid->GetOrdinaryGridIndex();
            bc[0]->SetCurrentZoneIndex(currentZoneIndex);
            bc[0]->SetTargetRegionBlock(targetCurrentZoneIndex);
            bc[0]->SetWhetherMerge(MERGED);
        }
    }
}

int WhethercanMergeKDirect(vector < StructBC* > bc, int nBC)
{
    if (nBC > 1)
    {
        int iBCOrignalIst, iBCOrignalIed, iBCOrignalJst, iBCOrignalJed,
            jBCOrignalIst, jBCOrignalIed, jBCOrignalJst, jBCOrignalJed;

        for (int iBC = 1; iBC < nBC; ++ iBC)
        {
            for (int jBC = 0; jBC < iBC; ++ jBC)
            {
                bc[iBC]->ComputeOriginal();
                bc[jBC]->ComputeOriginal();
                iBCOrignalIst = bc[iBC]->GetStartOrignal(0);
                iBCOrignalJst = bc[iBC]->GetStartOrignal(1);
                iBCOrignalIed = bc[iBC]->GetEndOrignal(0);
                iBCOrignalJed = bc[iBC]->GetEndOrignal(1);
                jBCOrignalIst = bc[jBC]->GetStartOrignal(0);
                jBCOrignalJst = bc[jBC]->GetStartOrignal(1);
                jBCOrignalIed = bc[jBC]->GetEndOrignal(0);
                jBCOrignalJed = bc[jBC]->GetEndOrignal(1);

                if (((iBCOrignalIst == jBCOrignalIed) && ((iBCOrignalJst == jBCOrignalJst) && (iBCOrignalJed == jBCOrignalJed)))
                 || ((iBCOrignalIed == jBCOrignalIst) && ((iBCOrignalJst == jBCOrignalJst) && (iBCOrignalJed == jBCOrignalJed)))
                 || ((iBCOrignalJst == jBCOrignalJed) && ((iBCOrignalIst == jBCOrignalIst) && (iBCOrignalIed == jBCOrignalIed)))
                 || ((iBCOrignalJed == jBCOrignalJst) && ((iBCOrignalIst == jBCOrignalIst) && (iBCOrignalIed == jBCOrignalIed))))
                {
                    return 1;
                }
            }
        }
    }
    else
    {
        return 0;
    }
}

void MergeBCK(vector < StructBC* > &bc, int *nBC)
{
    if (WhethercanMergeKDirect(bc, *nBC))
    {
        int iBCOrignalIst, iBCOrignalIed, iBCOrignalJst, iBCOrignalJed,
            jBCOrignalIst, jBCOrignalIed, jBCOrignalJst, jBCOrignalJed;
        int iBCFind = 1, jBCFind = 0;
        StructGrid *jBCGrid;
        int jBCGridIndex;
        int whethertheSame = 0;

        for (int iBC = 1; iBC < *nBC; ++ iBC)
        {
            for (int jBC = 0; jBC < iBC; ++ jBC)
            {
                bc[iBC]->ComputeOriginal();
                bc[jBC]->ComputeOriginal();
                iBCOrignalIst = bc[iBC]->GetStartOrignal(0);
                iBCOrignalJst = bc[iBC]->GetStartOrignal(1);
                iBCOrignalIed = bc[iBC]->GetEndOrignal(0);
                iBCOrignalJed = bc[iBC]->GetEndOrignal(1);
                jBCOrignalIst = bc[jBC]->GetStartOrignal(0);
                jBCOrignalJst = bc[jBC]->GetStartOrignal(1);
                jBCOrignalIed = bc[jBC]->GetEndOrignal(0);
                jBCOrignalJed = bc[jBC]->GetEndOrignal(1);

                if (((iBCOrignalIst == jBCOrignalIed) && ((iBCOrignalJst == jBCOrignalJst) && (iBCOrignalJed == jBCOrignalJed)))
                 || ((iBCOrignalIed == jBCOrignalIst) && ((iBCOrignalJst == jBCOrignalJst) && (iBCOrignalJed == jBCOrignalJed)))
                 || ((iBCOrignalJst == jBCOrignalJed) && ((iBCOrignalIst == jBCOrignalIst) && (iBCOrignalIed == jBCOrignalIed)))
                 || ((iBCOrignalJed == jBCOrignalJst) && ((iBCOrignalIst == jBCOrignalIst) && (iBCOrignalIed == jBCOrignalIed))))
                {
                    whethertheSame = 1;
                    iBCFind        = iBC;
                    jBCFind        = jBC;
                    break;
                }
            }
            if (1 == whethertheSame)
            {
                break;
            }
        }

        int    iBCType   = bc[jBCFind]->GetBCType();
        string iBCName   = bc[jBCFind]->GetBCName();
        string iBodyName = bc[jBCFind]->GetBodyName();
        int bcDirectionNew, targetBCDirectionNew;
        int leftorRightNew, targetLeftorRightNew;
        int currentZoneIndex, targetCurrentZoneIndex;
        int stOrignalNew[3], edOrignalNew[3];
        bcDirectionNew       = bc[jBCFind]->GetFaceDirection();
        targetBCDirectionNew = bc[jBCFind]->GetFaceDirectionOfTargetBlock();
        leftorRightNew       = bc[jBCFind]->GetFaceLeftOrRightIndex();
        targetLeftorRightNew = bc[jBCFind]->GetFaceLeftOrRightIndexOfTargetBlock();
        if ((bc[jBCFind]->GetWhetherMerge()) != MERGED)
        {
            jBCGridIndex           = bc[jBCFind]->GetRegionBlock();
            jBCGrid                = StructGridCast(GetGrid(jBCGridIndex, 0));
            jBCGridIndex           = bc[jBCFind]->GetTargetRegionBlock();
            jBCGrid                = StructGridCast(GetGrid(jBCGridIndex, 0));
            currentZoneIndex       = jBCGrid->GetOrdinaryGridIndex();
            targetCurrentZoneIndex = jBCGrid->GetOrdinaryGridIndex();
        }
        else
        {
            currentZoneIndex = bc[jBCFind]->GetRegionBlock();
            targetCurrentZoneIndex = bc[jBCFind]->GetTargetRegionBlock();
        }

        if ((iBCOrignalIst == jBCOrignalIed)
         && (iBCOrignalJst == jBCOrignalJst)
         && (iBCOrignalJed == jBCOrignalJed))
        {
            stOrignalNew[2] = bc[jBCFind]->GetStartOrignal(2);
            stOrignalNew[0] = min(jBCOrignalIst, iBCOrignalIst);
            stOrignalNew[1] = jBCOrignalJst;
            edOrignalNew[2] = bc[jBCFind]->GetEndOrignal(2);
            edOrignalNew[0] = max(iBCOrignalIed, jBCOrignalIed);
            edOrignalNew[1] = jBCOrignalJed;
        }
        else if ((iBCOrignalIed == jBCOrignalIst)
              && (iBCOrignalJst == jBCOrignalJst)
              && (iBCOrignalJed == jBCOrignalJed))
        {
            stOrignalNew[2] = bc[jBCFind]->GetStartOrignal(2);
            stOrignalNew[0] = min(jBCOrignalIst, iBCOrignalIst);
            stOrignalNew[1] = jBCOrignalJst;
            edOrignalNew[2] = bc[jBCFind]->GetEndOrignal(2);
            edOrignalNew[0] = max(iBCOrignalIed, jBCOrignalIed);
            edOrignalNew[1] = jBCOrignalJed;
        }
        else if ((iBCOrignalJst == jBCOrignalJed)
              && (iBCOrignalIst == jBCOrignalIst)
              && (iBCOrignalIed == jBCOrignalIed))
        {
            stOrignalNew[2] = bc[jBCFind]->GetStartOrignal(2);
            stOrignalNew[0] = jBCOrignalIst;
            stOrignalNew[1] = min(jBCOrignalJst, iBCOrignalJst);
            edOrignalNew[2] = bc[jBCFind]->GetEndOrignal(2);
            edOrignalNew[0] = jBCOrignalIed;
            edOrignalNew[1] = max(iBCOrignalJed, jBCOrignalJed);
        }
        else if ((iBCOrignalJed == jBCOrignalJst)
              && (iBCOrignalIst == jBCOrignalIst)
              && (iBCOrignalIed == jBCOrignalIed))
        {
            stOrignalNew[2] = bc[jBCFind]->GetStartOrignal(2);
            stOrignalNew[0] = jBCOrignalIst;
            stOrignalNew[1] = min(jBCOrignalJst, iBCOrignalJst);
            edOrignalNew[2] = bc[jBCFind]->GetEndOrignal(2);
            edOrignalNew[0] = jBCOrignalIed;
            edOrignalNew[1] = max(iBCOrignalJed, jBCOrignalJed);
        }

        bc[jBCFind]->SetCurrentZoneIndex(currentZoneIndex);
        bc[jBCFind]->SetFaceDirection(bcDirectionNew);
        bc[jBCFind]->SetFaceDirectionOfTargetBlock(targetBCDirectionNew);
        bc[jBCFind]->SetFaceLeftOrRightIndex(leftorRightNew);
        bc[jBCFind]->SetFaceLeftOrRightIndexOfTargetBlock(targetLeftorRightNew);
        bc[jBCFind]->SetTargetRegionBlock(targetCurrentZoneIndex);
        bc[jBCFind]->SetStartOrignal(stOrignalNew);
        bc[jBCFind]->SetEndOrignal(edOrignalNew);
        bc[jBCFind]->SetWhetherMerge(MERGED);
        bc[jBCFind]->SetBCType(iBCType);
        bc[jBCFind]->SetBCName(iBCName);
        bc[jBCFind]->SetBodyName(iBodyName);

        for (int mBC = iBCFind; mBC < (*nBC - 1); ++ mBC)
        {
            bc[mBC]->SetBCType(bc[mBC + 1]->GetBCType());
            bc[mBC]->SetBCName(bc[mBC + 1]->GetBCName());
            bc[mBC]->SetBodyName(bc[mBC + 1]->GetBodyName());
            bc[mBC] = bc[mBC + 1];
        }
        *nBC = *nBC - 1;
        MergeBCK(bc, nBC);
    }
    else
    {
        if ((bc[0]->GetWhetherMerge()) != MERGED)
        {
            int currentZoneIndex, targetCurrentZoneIndex;
            int bcGridIndex;
            StructGrid *bcGrid;
            bcGridIndex            = bc[0]->GetRegionBlock();
            bcGrid                 = StructGridCast(GetGrid(bcGridIndex, 0));
            currentZoneIndex       = bcGrid->GetOrdinaryGridIndex();
            bcGridIndex            = bc[0]->GetTargetRegionBlock();
            bcGrid                 = StructGridCast(GetGrid(bcGridIndex, 0));
            targetCurrentZoneIndex = bcGrid->GetOrdinaryGridIndex();
            bc[0]->SetCurrentZoneIndex(currentZoneIndex);
            bc[0]->SetTargetRegionBlock(targetCurrentZoneIndex);
            bc[0]->SetWhetherMerge(MERGED);
        }
    }
}

//! Select all start/end faces of i/j/k direction ,divide them by type and merge them
void SelectDivideMerge(int iDirection, int iFacePosition,
    vector < StructBC* > &OriginalBCRegion, int &nOriginalBCRegion,
    vector < vector < StructBC* > > &bcRegionMerged, int &nBCType, vector < int > &numberofEveryType)
{
    //! Extract all the boundaries whose directions are i and that are on the left.
    int nOriginalBCRegionSelected = 0;
    vector < StructBC* > originalBCRegionSelected;
    vector < Informationofboundary > originalBCRegionSelectedInformation;
    vector < vector < StructBC* > > originalBCRegionDivided;

    for (int iRegion = 0; iRegion < nOriginalBCRegion; ++ iRegion)
    {
        int mTest, nTest;
        mTest = OriginalBCRegion[iRegion]->GetFaceDirection();
        OriginalBCRegion[iRegion]->ComputeOriginal();
        nTest = OriginalBCRegion[iRegion]->GetStartOrignal(iDirection);
        if ((mTest == iDirection) && (nTest == iFacePosition))
        {
            nOriginalBCRegionSelected ++;
            originalBCRegionSelected.push_back(OriginalBCRegion[iRegion]);
        }
    }

    //! Divide all these boundries whose directions are i and that are on the left by boundry condition.
    originalBCRegionSelectedInformation.resize(nOriginalBCRegionSelected);
    for (int iregion = 0; iregion < nOriginalBCRegionSelected; ++ iregion)
    {
        originalBCRegionSelectedInformation[iregion].bcType       = originalBCRegionSelected[iregion]->GetBCType();
        originalBCRegionSelectedInformation[iregion].bodyName     = originalBCRegionSelected[iregion]->GetBodyName();
        originalBCRegionSelectedInformation[iregion].boundaryName = originalBCRegionSelected[iregion]->GetBCName();
    }

    int bcTypeCount = 0;
    originalBCRegionDivided.resize(1);
    originalBCRegionDivided[bcTypeCount].resize(1);
    originalBCRegionDivided[bcTypeCount][0] = originalBCRegionSelected[0];
    ++ bcTypeCount;
    for (int iFace = 1; iFace < nOriginalBCRegionSelected; ++ iFace)
    {
        for (int jFace = 0; jFace < iFace; ++ jFace)
        {
            if (!WhetherInformationofBoundarytheSame(originalBCRegionSelectedInformation[iFace], originalBCRegionSelectedInformation[jFace]))
            {
                if (jFace == (iFace - 1))
                {
                    ++ bcTypeCount;
                }
            }
            else
            {
                break;
            }
        }
    }
    nBCType = bcTypeCount;
    numberofEveryType.resize(nBCType);
    for (int iBCtype = 0; iBCtype < nBCType; ++ iBCtype)
    {
        numberofEveryType[iBCtype] = 0;
    }
    bcTypeCount = 0;
    originalBCRegionDivided[bcTypeCount][0] = originalBCRegionSelected[0];
    ++ (numberofEveryType[bcTypeCount]);
    ++ bcTypeCount;
    originalBCRegionDivided.resize(nBCType);

    for (int iFace = 1; iFace < nOriginalBCRegionSelected; ++ iFace)
    {
        for (int jFace = 0; jFace < iFace; ++ jFace)
        {
            if (WhetherInformationofBoundarytheSame(originalBCRegionSelectedInformation[iFace], originalBCRegionSelectedInformation[jFace]))
            {
                for (int iBCType = 0; iBCType < bcTypeCount;++ iBCType)
                {
                    if (WhetherInformationofBoundarytheSame((originalBCRegionSelected[iFace]), (originalBCRegionDivided[iBCType][0])))
                    {
                        originalBCRegionDivided[iBCType].push_back(originalBCRegionSelected[iFace]);
                        ++ (numberofEveryType[iBCType]);
                    }
                }
                break;
            }
            else
            {
                if (!WhetherInformationofBoundarytheSame(originalBCRegionSelectedInformation[iFace], originalBCRegionSelectedInformation[jFace]))
                {
                    if (jFace == (iFace - 1))
                    {
                        originalBCRegionDivided[bcTypeCount].resize(1);
                        originalBCRegionDivided[bcTypeCount][numberofEveryType[bcTypeCount]] = originalBCRegionSelected[iFace];
                        ++ (numberofEveryType[bcTypeCount]);
                        ++ bcTypeCount;
                    }
                }
            }
        }
    }

    //! Merge all these boundries whose directions are i and that are on the left.
    for (int iBCType = 0; iBCType < nBCType; ++ iBCType)
    {
        if (iDirection == 0)
        {
            MergeBCI(originalBCRegionDivided[iBCType], &(numberofEveryType[iBCType]));
        }
        else if (iDirection == 1)
        {
            MergeBCJ(originalBCRegionDivided[iBCType], &(numberofEveryType[iBCType]));
        }
        else if (iDirection == 2)
        {
            MergeBCK(originalBCRegionDivided[iBCType], &(numberofEveryType[iBCType]));
        }
    }

    bcRegionMerged.resize(nBCType);
    for (int iBCType=0; iBCType < nBCType; ++ iBCType)
    {
        bcRegionMerged[iBCType].resize(numberofEveryType[iBCType]);
        for (int iBC=0; iBC < numberofEveryType[iBCType]; ++ iBC)
        {
            bcRegionMerged[iBCType][iBC] = originalBCRegionDivided[iBCType][iBC];
        }
    }
}

OriginalStructGridMerge::OriginalStructGridMerge()
{
    nTotalZone        = 0;
    mergegrid         = 0;
    ni                = 0;
    nj                = 0;
    nk                = 0;
    nTotalNode        = 0;
    nTotalFace        = 0;
    nTotalCell        = 0;
    coordinates       = 0;
    nBCRegion         = 0;
    compositeBCRegion = 0;
    nIFace            = 0;
    interfaceInfo     = 0;
    nTotalblock       = 0;
    nDim              = 0;
}

OriginalStructGridMerge::~OriginalStructGridMerge()
{
    nTotalZone = 0;
    if (mergegrid)
    {
        DelPointer(mergegrid);
        DelPointer(ni);
        DelPointer(nj);
        DelPointer(nk);
        DelPointer(nTotalNode);
        DelPointer(nTotalFace);
        DelPointer(nTotalCell);

        for (int iBlock = 0; iBlock < nTotalblock; ++ iBlock)
        {
            DelPointer(coordinates[iBlock]);
        }

        DelPointer(coordinates);
        DelPointer(nBCRegion);
        DelPointer(compositeBCRegion);
        DelPointer(nIFace);
        DelPointer(interfaceInfo);
    }
    nTotalblock = 0;
    nDim        = 0;
}

void OriginalStructGridMerge::Initial(vector < vector < Grid* > > originalgridofBlockIn, int nTotalblock, int nTotalZone, int nDim)
{
    this->nTotalblock         = nTotalblock;
    this->nTotalZone          = nTotalZone;
    originalgridofBlock.resize(nTotalblock);
    this->nDim = nDim;

    vector < vector < int > > oriEndIndexIDir;
    oriEndIndexIDir.resize(nTotalblock);
    vector < vector < int > > oriEndIndexJDir;
    oriEndIndexJDir.resize(nTotalblock);
    vector < vector < int > > oriEndIndexKDir;
    oriEndIndexKDir.resize(nTotalblock);

    vector < vector < int > > oriStartIndexIDir;
    oriStartIndexIDir.resize(nTotalblock);
    vector < vector < int > > oriStartIndexJDir;
    oriStartIndexJDir.resize(nTotalblock);
    vector < vector < int > > oriStartIndexKDir;
    oriStartIndexKDir.resize(nTotalblock);

    vector < vector < int > > oriNI;
    oriNI.resize(nTotalblock);
    vector < vector < int > > oriNJ;
    oriNJ.resize(nTotalblock);
    vector < vector < int > > oriNK;
    oriNK.resize(nTotalblock);
    int ZERO = 0;
    int ONE  = 1;
    int TWO  = 2;

    vector < vector < int > > nDumpOriginalBCRegion;
    nDumpOriginalBCRegion.resize(nTotalblock);
    dumpOriginalBCRegion.resize(nTotalblock);
    int *nOriginalBCRegion = new int[nTotalblock];

    vector < vector < vector < StructBC* > > > bcRegionMergedIDirectionStart;
    bcRegionMergedIDirectionStart.resize(nTotalblock);
    int *nBCTypeIDirectionStart = new int [nTotalblock];
    vector < vector < int > > numberofEveryTypeIStart;
    numberofEveryTypeIStart.resize(nTotalblock);

    vector < vector < vector < StructBC* > > > bcRegionMergedIDirectionEnd;
    bcRegionMergedIDirectionEnd.resize(nTotalblock);
    int *nBCTypeIDirectionEnd = new int [nTotalblock];
    vector < vector < int > > numberofEveryTypeIEnd;
    numberofEveryTypeIEnd.resize(nTotalblock);

    vector < vector < vector < StructBC* > > > bcRegionMergedJDirectionStart;
    bcRegionMergedJDirectionStart.resize(nTotalblock);
    int *nBCTypeJDirectionStart = new int[nTotalblock];
    vector < vector < int > > numberofEveryTypeJStart;
    numberofEveryTypeJStart.resize(nTotalblock);

    vector < vector < vector < StructBC* > > > bcRegionMergedJDirectionEnd;
    bcRegionMergedJDirectionEnd.resize(nTotalblock);
    int *nBCTypeJDirectionEnd = new int[nTotalblock];
    vector < vector < int > > numberofEveryTypeJEnd;
    numberofEveryTypeJEnd.resize(nTotalblock);

    vector < vector < vector < StructBC* > > > bcRegionMergedKDirectionStart;
    bcRegionMergedKDirectionStart.resize(nTotalblock);
    int *nBCTypeKDirectionStart = new int[nTotalblock];
    vector < vector < int > > numberofEveryTypeKStart;
    numberofEveryTypeKStart.resize(nTotalblock);

    vector < vector < vector < StructBC* > > > bcRegionMergedKDirectionEnd;
    bcRegionMergedKDirectionEnd.resize(nTotalblock);
    int *nBCTypeKDirectionEnd = new int[nTotalblock];
    vector < vector < int > > numberofEveryTypeKEnd;
    numberofEveryTypeKEnd.resize(nTotalblock);

    ni                = new int [nTotalblock];
    nj                = new int [nTotalblock];
    nk                = new int [nTotalblock];
    nTotalNode        = new int [nTotalblock];
    coordinates       = new RDouble **[nTotalblock];
    nBCRegion         = new int [nTotalblock];
    compositeBCRegion = new StructBCSet *[nTotalblock];
    nIFace            = new int [nTotalblock];
    mergegrid         = new StructGrid *[nTotalblock];
    for (int iBlock = 0; iBlock < nTotalblock; ++ iBlock)
    {
        coordinates[iBlock] = new RDouble *[3];
    }

    for (int iBlock = 0; iBlock < nTotalblock; ++ iBlock)
    {
        if(!originalgridofBlockIn[iBlock].size())
        {
            continue;
        }

        for (int iZone = 0; iZone < originalgridofBlockIn[iBlock].size(); ++ iZone)
        {
            originalgridofBlock[iBlock].push_back(StructGridCast(originalgridofBlockIn[iBlock][iZone]));
        }

        //! Calculaete ni¡¢nj¡¢nk of every block.
        for (int iZone = 0; iZone < originalgridofBlock[iBlock].size(); ++ iZone)
        {
            oriEndIndexIDir[iBlock].push_back(*(originalgridofBlock[iBlock][iZone]->GetOrdinaryDimEndIndex() + ZERO));
            oriEndIndexJDir[iBlock].push_back(*(originalgridofBlock[iBlock][iZone]->GetOrdinaryDimEndIndex() + ONE));
            oriEndIndexKDir[iBlock].push_back(*(originalgridofBlock[iBlock][iZone]->GetOrdinaryDimEndIndex() + TWO));
            oriStartIndexIDir[iBlock].push_back(*(originalgridofBlock[iBlock][iZone]->GetOrdinaryDimStartIndex() + ZERO));
            oriStartIndexJDir[iBlock].push_back(*(originalgridofBlock[iBlock][iZone]->GetOrdinaryDimStartIndex() + ONE));
            oriStartIndexKDir[iBlock].push_back(*(originalgridofBlock[iBlock][iZone]->GetOrdinaryDimStartIndex() + TWO));
            oriNI[iBlock].push_back(originalgridofBlock[iBlock][iZone]->GetNI());
            oriNJ[iBlock].push_back(originalgridofBlock[iBlock][iZone]->GetNJ());
            oriNK[iBlock].push_back(originalgridofBlock[iBlock][iZone]->GetNK());
        }
        ni[iBlock] = 0;
        nj[iBlock] = 0;
        nk[iBlock] = 0;
        for (int iZone = 0; iZone < originalgridofBlock[iBlock].size(); ++ iZone)
        {
            ni[iBlock] = max(oriEndIndexIDir[iBlock][iZone], ni[iBlock]);
            nj[iBlock] = max(oriEndIndexJDir[iBlock][iZone], nj[iBlock]);
            nk[iBlock] = max(oriEndIndexKDir[iBlock][iZone], nk[iBlock]);
        }

        //! Calculate the coordinates of the point of every block.
        nTotalNode[iBlock] = ni[iBlock] * nj[iBlock] * nk[iBlock];
        RDouble *mergegridX = new RDouble [nTotalNode[iBlock]];
        RDouble *mergegridY = new RDouble [nTotalNode[iBlock]];
        RDouble *mergegridZ = new RDouble [nTotalNode[iBlock]];
        for (int iZone = 0; iZone < originalgridofBlock[iBlock].size(); ++ iZone)
        {
            RDouble3D &structx = *originalgridofBlock[iBlock][iZone]->GetStructX();
            RDouble3D &structy = *originalgridofBlock[iBlock][iZone]->GetStructY();
            RDouble3D &structz = *originalgridofBlock[iBlock][iZone]->GetStructZ();

            for (int iz = 1; iz <= oriNK[iBlock][iZone]; ++ iz)
            {
                for (int iy = 1; iy <= oriNJ[iBlock][iZone]; ++ iy)
                {
                    for (int ix = 1; ix <= oriNI[iBlock][iZone]; ++ ix)
                    {
                        int IX, IY, IZ, ix0, iy0, iz0;
                        ix0 = oriStartIndexIDir[iBlock][iZone];
                        iy0 = oriStartIndexJDir[iBlock][iZone];
                        iz0 = oriStartIndexKDir[iBlock][iZone];
                        int Dx, Dy, Dz;
                        Dx                                      = ni[iBlock];
                        Dy                                      = nj[iBlock];
                        Dz                                      = nk[iBlock];
                        IX                                      = ix + ix0 - 2;
                        IY                                      = iy + iy0 - 2;
                        IZ                                      = iz + iz0 - 2;
                        mergegridX[IX + IY * Dx + IZ * Dx * Dy] = structx(ix, iy, iz);
                        mergegridY[IX + IY * Dx + IZ * Dx * Dy] = structy(ix, iy, iz);
                        mergegridZ[IX + IY * Dx + IZ * Dx * Dy] = structz(ix, iy, iz);
                    }
                }
            }
        }

        coordinates[iBlock][0] = new RDouble [nTotalNode[iBlock]];
        coordinates[iBlock][1] = new RDouble [nTotalNode[iBlock]];
        coordinates[iBlock][2] = new RDouble [nTotalNode[iBlock]];

        for (int iNode = 0; iNode < nTotalNode[iBlock]; ++ iNode)
        {
            coordinates[iBlock][0][iNode] = mergegridX[iNode];
            coordinates[iBlock][1][iNode] = mergegridY[iNode];
            coordinates[iBlock][2][iNode] = mergegridZ[iNode];
        }

        DelPointer(mergegridX);
        DelPointer(mergegridY);
        DelPointer(mergegridZ);

        nBCRegion[iBlock]         = 0;
        nOriginalBCRegion[iBlock] = 0;
        for (int iZone = 0; iZone < originalgridofBlock[iBlock].size(); ++ iZone)
        {
            nDumpOriginalBCRegion[iBlock].push_back(originalgridofBlock[iBlock][iZone]->GetStructBCSet()->GetnBCRegion());
            nOriginalBCRegion[iBlock] = nOriginalBCRegion[iBlock] + nDumpOriginalBCRegion[iBlock][iZone];
        }

        //! Divide all the boundaries by block.
        for (int iZone = 0; iZone < originalgridofBlock[iBlock].size(); ++ iZone)
        {
            originalgridofBlock[iBlock][iZone]->ComputeBCOriginal();
            for (int iRegion = 0; iRegion < nDumpOriginalBCRegion[iBlock][iZone]; ++iRegion)
            {
                StructBC *newBC = new StructBC(originalgridofBlock[iBlock][iZone]->GetZoneID(), iRegion);
                StructBC *oldBC = originalgridofBlock[iBlock][iZone]->GetStructBCSet()->GetBCRegion(iRegion);
                newBC->CopyBC(oldBC);
                dumpOriginalBCRegion[iBlock].push_back(newBC);
            }
        }

        SelectDivideMerge(0, 1, dumpOriginalBCRegion[iBlock], nOriginalBCRegion[iBlock], bcRegionMergedIDirectionStart[iBlock], nBCTypeIDirectionStart[iBlock], numberofEveryTypeIStart[iBlock]);
        SelectDivideMerge(0, ni[iBlock], dumpOriginalBCRegion[iBlock], nOriginalBCRegion[iBlock], bcRegionMergedIDirectionEnd[iBlock], nBCTypeIDirectionEnd[iBlock], numberofEveryTypeIEnd[iBlock]);
        SelectDivideMerge(1, 1, dumpOriginalBCRegion[iBlock], nOriginalBCRegion[iBlock], bcRegionMergedJDirectionStart[iBlock], nBCTypeJDirectionStart[iBlock], numberofEveryTypeJStart[iBlock]);
        SelectDivideMerge(1, nj[iBlock], dumpOriginalBCRegion[iBlock], nOriginalBCRegion[iBlock], bcRegionMergedJDirectionEnd[iBlock], nBCTypeJDirectionEnd[iBlock], numberofEveryTypeJEnd[iBlock]);
        if (nDim != TWO_D)
        {
            SelectDivideMerge(2, 1, dumpOriginalBCRegion[iBlock], nOriginalBCRegion[iBlock], bcRegionMergedKDirectionStart[iBlock], nBCTypeKDirectionStart[iBlock], numberofEveryTypeKStart[iBlock]);
            SelectDivideMerge(2, nk[iBlock], dumpOriginalBCRegion[iBlock], nOriginalBCRegion[iBlock], bcRegionMergedKDirectionEnd[iBlock], nBCTypeKDirectionEnd[iBlock], numberofEveryTypeKEnd[iBlock]);
        }

        //! Calculate the number of total of boundaries of every block.
        nBCRegion[iBlock] = 0;
        for (int iBctype = 0; iBctype < nBCTypeIDirectionStart[iBlock]; ++ iBctype)
        {
            nBCRegion[iBlock] = nBCRegion[iBlock] + numberofEveryTypeIStart[iBlock][iBctype];
        }
        for (int iBctype = 0; iBctype < nBCTypeIDirectionEnd[iBlock]; ++ iBctype)
        {
            nBCRegion[iBlock] = nBCRegion[iBlock] + numberofEveryTypeIEnd[iBlock][iBctype];
        }
        for (int iBctype = 0; iBctype < nBCTypeJDirectionStart[iBlock]; ++ iBctype)
        {
            nBCRegion[iBlock] = nBCRegion[iBlock] + numberofEveryTypeJStart[iBlock][iBctype];
        }
        for (int iBctype = 0; iBctype < nBCTypeJDirectionEnd[iBlock]; ++ iBctype)
        {
            nBCRegion[iBlock] = nBCRegion[iBlock] + numberofEveryTypeJEnd[iBlock][iBctype];
        }
        if (nDim != TWO_D)
        {
            for (int iBctype = 0; iBctype < nBCTypeKDirectionStart[iBlock]; ++iBctype)
            {
                nBCRegion[iBlock] = nBCRegion[iBlock] + numberofEveryTypeKStart[iBlock][iBctype];
            }
            for (int iBctype = 0; iBctype < nBCTypeKDirectionEnd[iBlock]; ++iBctype)
            {
                nBCRegion[iBlock] = nBCRegion[iBlock] + numberofEveryTypeKEnd[iBlock][iBctype];
            }
        }

        //! Construct of StructBCset
        compositeBCRegion[iBlock] = new StructBCSet(iBlock);
        compositeBCRegion[iBlock]->CreateBCRegion(nBCRegion[iBlock]);
        int iBCCount = 0;
        for (int iBCType = 0; iBCType < nBCTypeIDirectionStart[iBlock]; ++ iBCType)
        {
            for (int iBC = 0; iBC < numberofEveryTypeIStart[iBlock][iBCType]; ++ iBC)
            {
                compositeBCRegion[iBlock]->SetBCRegion(iBCCount, bcRegionMergedIDirectionStart[iBlock][iBCType][iBC]);
                ++ iBCCount;
            }
        }
        for (int iBCType = 0; iBCType < nBCTypeIDirectionEnd[iBlock]; ++ iBCType)
        {
            for (int iBC = 0; iBC < numberofEveryTypeIEnd[iBlock][iBCType]; ++ iBC)
            {
                compositeBCRegion[iBlock]->SetBCRegion(iBCCount, bcRegionMergedIDirectionEnd[iBlock][iBCType][iBC]);
                ++ iBCCount;
            }
        }
        for (int iBCType = 0; iBCType < nBCTypeJDirectionStart[iBlock]; ++ iBCType)
        {
            for (int iBC = 0; iBC < numberofEveryTypeJStart[iBlock][iBCType]; ++ iBC)
            {
                compositeBCRegion[iBlock]->SetBCRegion(iBCCount, bcRegionMergedJDirectionStart[iBlock][iBCType][iBC]);
                ++ iBCCount;
            }
        }
        for (int iBCType = 0; iBCType < nBCTypeJDirectionEnd[iBlock]; ++ iBCType)
        {
            for (int iBC = 0; iBC < numberofEveryTypeJEnd[iBlock][iBCType]; ++ iBC)
            {
                compositeBCRegion[iBlock]->SetBCRegion(iBCCount, bcRegionMergedJDirectionEnd[iBlock][iBCType][iBC]);
                ++ iBCCount;
            }
        }
        if (nDim != TWO_D)
        {
            for (int iBCType = 0; iBCType < nBCTypeKDirectionStart[iBlock]; ++iBCType)
            {
                for (int iBC = 0; iBC < numberofEveryTypeKStart[iBlock][iBCType]; ++iBC)
                {
                    compositeBCRegion[iBlock]->SetBCRegion(iBCCount, bcRegionMergedKDirectionStart[iBlock][iBCType][iBC]);
                    ++iBCCount;
                }
            }
            for (int iBCType = 0; iBCType < nBCTypeKDirectionEnd[iBlock]; ++iBCType)
            {
                for (int iBC = 0; iBC < numberofEveryTypeKEnd[iBlock][iBCType]; ++iBC)
                {
                    compositeBCRegion[iBlock]->SetBCRegion(iBCCount, bcRegionMergedKDirectionEnd[iBlock][iBCType][iBC]);
                    ++iBCCount;
                }
            }
        }
    }

    DelPointer(nOriginalBCRegion);
    DelPointer(nBCTypeIDirectionStart);
    DelPointer(nBCTypeIDirectionEnd);
    DelPointer(nBCTypeJDirectionStart);
    DelPointer(nBCTypeJDirectionEnd);
    DelPointer(nBCTypeKDirectionStart);
    DelPointer(nBCTypeKDirectionEnd);
}

void OriginalStructGridMerge::Run()
{
    for (int iBlock = 0; iBlock < nTotalblock; ++ iBlock)
    {
        if (0 == dumpOriginalBCRegion[iBlock].size())
        {
            mergegrid[iBlock] = 0;
            continue;
        }

        GridID *index = new GridID(iBlock);
        mergegrid[iBlock] = StructGridCast(CreateGridGeneral(STRUCTGRID, index, 0, nDim));
        mergegrid[iBlock]->SetNI(ni[iBlock]);
        mergegrid[iBlock]->SetNJ(nj[iBlock]);
        mergegrid[iBlock]->SetNK(nk[iBlock]);
        mergegrid[iBlock]->SetBasicDimension();
        mergegrid[iBlock]->SetX(coordinates[iBlock][0]);
        mergegrid[iBlock]->SetY(coordinates[iBlock][1]);
        mergegrid[iBlock]->SetZ(coordinates[iBlock][2]);
        mergegrid[iBlock]->RotateAboutAxis();
        mergegrid[iBlock]->ComputeMinMaxBox();
        mergegrid[iBlock]->SetArrayLayout();
        mergegrid[iBlock]->CreateCompositeBCRegion(nBCRegion[iBlock]);
        mergegrid[iBlock]->CopyStructBCSet(compositeBCRegion[iBlock]);
        mergegrid[iBlock]->SetNIFace(nIFace[iBlock]);
        mergegrid[iBlock]->SetOrdinaryGridIndex(iBlock);
    }
}

void OriginalStructGridMerge::GetMergeGrid(StructGrid ** mergeGridOut)
{
    for (int iBlock = 0; iBlock < nTotalblock; ++ iBlock)
    {
        mergeGridOut[iBlock] = mergegrid[iBlock];
    }
}

}