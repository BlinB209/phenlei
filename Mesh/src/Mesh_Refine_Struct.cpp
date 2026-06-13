#include "Mesh_Refine_Struct.h"
#include "Glb_Dimension.h"
#include "Pre_GridConversion.h"
#include "Zone.h"

using namespace std;

namespace PHSPACE
{
LIB_EXPORT Mesh_Refine_Struct::Mesh_Refine_Struct(const string &gridFileNameIn) :
    Mesh_Refine(gridFileNameIn)
{
    OrdinaryGrid.resize(0);
}

LIB_EXPORT Mesh_Refine_Struct::~Mesh_Refine_Struct()
{

}

void Mesh_Refine_Struct::AllocateMemory()
{
    numberOfZones = PHMPI::GetNumberofLocalZones();
    refinedGrids = new Grid * [numberOfZones];

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Zone *zone = region->GetZone(iZone);
        if (!zone) continue;

        Grid *grid = zone->GetGeometry()->GetGrid();
        OrdinaryGrid.push_back(grid);
    }
}

void Mesh_Refine_Struct::RefineGrid()
{
    for (int iZone = 0; iZone < numberOfZones; iZone ++)
    {
        StructGrid *grid = StructGridCast(OrdinaryGrid[iZone]);
        RefineGrid(iZone, grid);
    }
}

void Mesh_Refine_Struct::RefineGrid(int iZone, StructGrid *gridIn)
{
    int newNodeI, newNodeJ, newNodeK;

    GridID *indexNew = new GridID(iZone);
    Grid *gridNew = CreateGridGeneral(STRUCTGRID, indexNew, 0, GetDim());
    refinedGrids[iZone] = gridNew;
    StructGrid *strGridNew = StructGridCast(gridNew);

    int ni = gridIn->GetNI();
    int nj = gridIn->GetNJ();
    int nk = gridIn->GetNK();

    strGridNew->SetNI((ni-1)*2+1);
    strGridNew->SetNJ((nj-1)*2+1);
    strGridNew->SetNK((nk-1)*2+1);
    strGridNew->SetBasicDimension();

    int nTotalNode = strGridNew->GetNTotalNode();
    RDouble *nx = new RDouble [nTotalNode];
    RDouble *ny = new RDouble [nTotalNode];
    RDouble *nz = new RDouble [nTotalNode];
    strGridNew->SetX(nx);
    strGridNew->SetY(ny);
    strGridNew->SetZ(nz);
    strGridNew->SetArrayLayout();

    int ist, ied, jst, jed, kst, ked;
    gridIn->GetNodeIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble3D &xx = *gridIn->GetStructX();
    RDouble3D &yy = *gridIn->GetStructY();
    RDouble3D &zz = *gridIn->GetStructZ();

    RDouble3D &nxx = *strGridNew->GetStructX();
    RDouble3D &nyy = *strGridNew->GetStructY();
    RDouble3D &nzz = *strGridNew->GetStructZ();

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                newNodeI = (i-1) * 2 + 1;
                newNodeJ = (j-1) * 2 + 1;
                newNodeK = (k-1) * 2 + 1;
                nxx(newNodeI, newNodeJ, newNodeK) = xx(i, j, k);
                nyy(newNodeI, newNodeJ, newNodeK) = yy(i, j, k);
                nzz(newNodeI, newNodeJ, newNodeK) = zz(i, j, k);
            }
        }
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied-1; ++ i)
            {
                newNodeI = (i-1) * 2 + 2;
                newNodeJ = (j-1) * 2 + 1;
                newNodeK = (k-1) * 2 + 1;
                nxx(newNodeI, newNodeJ, newNodeK) = 0.5 * (nxx(newNodeI-1, newNodeJ, newNodeK) + nxx(newNodeI+1, newNodeJ, newNodeK));
                nyy(newNodeI, newNodeJ, newNodeK) = 0.5 * (nyy(newNodeI-1, newNodeJ, newNodeK) + nyy(newNodeI+1, newNodeJ, newNodeK));
                nzz(newNodeI, newNodeJ, newNodeK) = 0.5 * (nzz(newNodeI-1, newNodeJ, newNodeK) + nzz(newNodeI+1, newNodeJ, newNodeK));
            }
        }
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed-1; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                newNodeI = (i-1) * 2 + 1;
                newNodeJ = (j-1) * 2 + 2;
                newNodeK = (k-1) * 2 + 1;
                nxx(newNodeI, newNodeJ, newNodeK) = 0.5 * (nxx(newNodeI, newNodeJ-1, newNodeK) + nxx(newNodeI, newNodeJ+1, newNodeK));
                nyy(newNodeI, newNodeJ, newNodeK) = 0.5 * (nyy(newNodeI, newNodeJ-1, newNodeK) + nyy(newNodeI, newNodeJ+1, newNodeK));
                nzz(newNodeI, newNodeJ, newNodeK) = 0.5 * (nzz(newNodeI, newNodeJ-1, newNodeK) + nzz(newNodeI, newNodeJ+1, newNodeK));
            }
        }
    }

    for (int k = kst; k <= ked-1; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                newNodeI = (i-1) * 2 + 1;
                newNodeJ = (j-1) * 2 + 1;
                newNodeK = (k-1) * 2 + 2;
                nxx(newNodeI, newNodeJ, newNodeK) = 0.5 * (nxx(newNodeI, newNodeJ, newNodeK-1) + nxx(newNodeI, newNodeJ, newNodeK+1));
                nyy(newNodeI, newNodeJ, newNodeK) = 0.5 * (nyy(newNodeI, newNodeJ, newNodeK-1) + nyy(newNodeI, newNodeJ, newNodeK+1));
                nzz(newNodeI, newNodeJ, newNodeK) = 0.5 * (nzz(newNodeI, newNodeJ, newNodeK-1) + nzz(newNodeI, newNodeJ, newNodeK+1));
            }
        }
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed-1; ++ j)
        {
            for (int i = ist; i <= ied-1; ++ i)
            {
                newNodeI = (i-1) * 2 + 2;
                newNodeJ = (j-1) * 2 + 2;
                newNodeK = (k-1) * 2 + 1;
                nxx(newNodeI, newNodeJ, newNodeK) = 0.25 * (nxx(newNodeI-1, newNodeJ, newNodeK) + nxx(newNodeI+1, newNodeJ, newNodeK) + nxx(newNodeI, newNodeJ-1, newNodeK) + nxx(newNodeI, newNodeJ+1, newNodeK));
                nyy(newNodeI, newNodeJ, newNodeK) = 0.25 * (nyy(newNodeI-1, newNodeJ, newNodeK) + nyy(newNodeI+1, newNodeJ, newNodeK) + nyy(newNodeI, newNodeJ-1, newNodeK) + nyy(newNodeI, newNodeJ+1, newNodeK));
                nzz(newNodeI, newNodeJ, newNodeK) = 0.25 * (nzz(newNodeI-1, newNodeJ, newNodeK) + nzz(newNodeI+1, newNodeJ, newNodeK) + nzz(newNodeI, newNodeJ-1, newNodeK) + nzz(newNodeI, newNodeJ+1, newNodeK));
            }
        }
    }

    for (int k = kst; k <= ked-1; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied-1; ++ i)
            {
                newNodeI = (i-1) * 2 + 2;
                newNodeJ = (j-1) * 2 + 1;
                newNodeK = (k-1) * 2 + 2;
                nxx(newNodeI, newNodeJ, newNodeK) = 0.25 * (nxx(newNodeI-1, newNodeJ, newNodeK) + nxx(newNodeI+1, newNodeJ, newNodeK) + nxx(newNodeI, newNodeJ, newNodeK-1) + nxx(newNodeI, newNodeJ, newNodeK+1));
                nyy(newNodeI, newNodeJ, newNodeK) = 0.25 * (nyy(newNodeI-1, newNodeJ, newNodeK) + nyy(newNodeI+1, newNodeJ, newNodeK) + nyy(newNodeI, newNodeJ, newNodeK-1) + nyy(newNodeI, newNodeJ, newNodeK+1));
                nzz(newNodeI, newNodeJ, newNodeK) = 0.25 * (nzz(newNodeI-1, newNodeJ, newNodeK) + nzz(newNodeI+1, newNodeJ, newNodeK) + nzz(newNodeI, newNodeJ, newNodeK-1) + nzz(newNodeI, newNodeJ, newNodeK+1));
            }
        }
    }

    for (int k = kst; k <= ked-1; ++ k)
    {
        for (int j = jst; j <= jed-1; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                newNodeI = (i-1) * 2 + 1;
                newNodeJ = (j-1) * 2 + 2;
                newNodeK = (k-1) * 2 + 2;
                nxx(newNodeI, newNodeJ, newNodeK) = 0.25 * (nxx(newNodeI, newNodeJ-1, newNodeK) + nxx(newNodeI, newNodeJ+1, newNodeK) + nxx(newNodeI, newNodeJ, newNodeK-1) + nxx(newNodeI, newNodeJ, newNodeK+1));
                nyy(newNodeI, newNodeJ, newNodeK) = 0.25 * (nyy(newNodeI, newNodeJ-1, newNodeK) + nyy(newNodeI, newNodeJ+1, newNodeK) + nyy(newNodeI, newNodeJ, newNodeK-1) + nyy(newNodeI, newNodeJ, newNodeK+1));
                nzz(newNodeI, newNodeJ, newNodeK) = 0.25 * (nzz(newNodeI, newNodeJ-1, newNodeK) + nzz(newNodeI, newNodeJ+1, newNodeK) + nzz(newNodeI, newNodeJ, newNodeK-1) + nzz(newNodeI, newNodeJ, newNodeK+1));
            }
        }
    }

    for (int k = kst; k <= ked-1; ++ k)
    {
        for (int j = jst; j <= jed-1; ++ j)
        {
            for (int i = ist; i <= ied-1; ++ i)
            {
                newNodeI = (i-1) * 2 + 2;
                newNodeJ = (j-1) * 2 + 2;
                newNodeK = (k-1) * 2 + 2;
                nxx(newNodeI, newNodeJ, newNodeK) = 0.125 * (nxx(newNodeI-1, newNodeJ-1, newNodeK-1) + nxx(newNodeI+1, newNodeJ-1, newNodeK-1) + nxx(newNodeI-1, newNodeJ+1, newNodeK-1) + nxx(newNodeI+1, newNodeJ+1, newNodeK-1)
                                                  + nxx(newNodeI-1, newNodeJ-1, newNodeK+1) + nxx(newNodeI+1, newNodeJ-1, newNodeK+1) + nxx(newNodeI-1, newNodeJ+1, newNodeK+1) + nxx(newNodeI+1, newNodeJ+1, newNodeK+1));
                nyy(newNodeI, newNodeJ, newNodeK) = 0.125 * (nyy(newNodeI-1, newNodeJ-1, newNodeK-1) + nyy(newNodeI+1, newNodeJ-1, newNodeK-1) + nyy(newNodeI-1, newNodeJ+1, newNodeK-1) + nyy(newNodeI+1, newNodeJ+1, newNodeK-1)
                                                  + nyy(newNodeI-1, newNodeJ-1, newNodeK+1) + nyy(newNodeI+1, newNodeJ-1, newNodeK+1) + nyy(newNodeI-1, newNodeJ+1, newNodeK+1) + nyy(newNodeI+1, newNodeJ+1, newNodeK+1));
                nzz(newNodeI, newNodeJ, newNodeK) = 0.125 * (nzz(newNodeI-1, newNodeJ-1, newNodeK-1) + nzz(newNodeI+1, newNodeJ-1, newNodeK-1) + nzz(newNodeI-1, newNodeJ+1, newNodeK-1) + nzz(newNodeI+1, newNodeJ+1, newNodeK-1)
                                                  + nzz(newNodeI-1, newNodeJ-1, newNodeK+1) + nzz(newNodeI+1, newNodeJ-1, newNodeK+1) + nzz(newNodeI-1, newNodeJ+1, newNodeK+1) + nzz(newNodeI+1, newNodeJ+1, newNodeK+1));
            }
        }
    }

    strGridNew->ComputeMinMaxBox();
}

void Mesh_Refine_Struct::GenerateAndDumpComputationalGrid()
{
    GenerateBCRegion();

    ConstructGlobalInterfaceLink(refinedGrids, numberOfZones);

    DumpComputationalGrid();
}

void Mesh_Refine_Struct::GenerateBCRegion()
{
    for (int iZone = 0; iZone < numberOfZones; iZone ++)
    {
        StructGrid *grid       = StructGridCast(OrdinaryGrid[iZone]);
        StructGrid *strGridNew = StructGridCast(refinedGrids[iZone]);

        strGridNew->CopyStructBCSet(grid->GetStructBCSet());
        StructBCSet *rhs = strGridNew->GetStructBCSet();

        int numberOfBoundaryFaces = rhs->GetnBCRegion();
        for (int iBCRegion = 0; iBCRegion < numberOfBoundaryFaces; iBCRegion ++)
        {
            StructBC *bcregion = rhs->GetBCRegion(iBCRegion);

            int bctype = bcregion->GetBCType();

            int *s_st   = bcregion->GetStartPoint();
            int *s_ed   = bcregion->GetEndPoint();
            int *t_st   = bcregion->GetTargetStart();
            int *t_ed   = bcregion->GetTargetEnd();
            int *s_lr3d = bcregion->GetFaceDirectionIndex();
            int *t_lr3d = bcregion->GetFaceDirectionIndexOfTargetBlock();

            for (int m = 0; m < GetDim(); ++ m)
            {
                if (s_st[m] == s_ed[m])
                {
                    if (s_st[m] > 1 || s_lr3d[m] == 1)
                    {
                        s_st[m] += 1;
                        s_ed[m] += 1;
                    }
                }
                else
                {
                    s_ed[m] += 1;
                }
            }

            for (int m = 0; m < GetDim(); ++ m)
            {
                s_st[m] = (s_st[m]-1) * 2 + 1;
                s_ed[m] = (s_ed[m]-1) * 2 + 1;
            }

            if (bctype == PHENGLEI::INTERFACE)
            {
                for (int m = 0; m < GetDim(); ++ m)
                {
                    if (t_st[m] == t_ed[m])
                    {
                        if (t_st[m] > 1 || t_lr3d[m] == 1)
                        {
                            t_st[m] += 1;
                            t_ed[m] += 1;
                        }
                    }
                    else
                    {
                        t_ed[m] += 1;
                    }
                }

                for (int m = 0; m < GetDim(); ++ m)
                {
                    t_st[m] = (t_st[m]-1) * 2 + 1;
                    t_ed[m] = (t_ed[m]-1) * 2 + 1;
                }
            }
        }

        strGridNew->ProcessBCInfo();
    }
}

}
