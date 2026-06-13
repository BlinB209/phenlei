#include "Pre_CGNSToPHengLEI_Struct.h"
#include "Geo_StructGrid.h"
#include "Glb_Dimension.h"
#include "TK_Log.h"

namespace PHSPACE
{

LIB_EXPORT Pre_CGNSToPHengLEI_Struct::Pre_CGNSToPHengLEI_Struct(const string &gridFileName) : 
        Pre_CGNSConversion_Struct(gridFileName)
{

}

LIB_EXPORT Pre_CGNSToPHengLEI_Struct::~Pre_CGNSToPHengLEI_Struct()
{

}

void Pre_CGNSToPHengLEI_Struct::Conversion()
{
    ResetGridScaleAndTranslate();

    grids = new Grid* [nBlocks];

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        GridID *index = new GridID(iZone);
        StructGrid *grid = StructGridCast(CreateGridGeneral(STRUCTGRID, index, 0, GetDim()));
        grids[iZone] = grid;
    }

    SetVolumeInfo();
    SetCoordinate();
    SetBCInfo();

    CheckMeshMultigrid();

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        Grid *grid = grids[iZone];
        grid->SetLnkInfo();
        grid->ProcessBCInfo();
    }
}

void Pre_CGNSToPHengLEI_Struct::SetCoordinate()
{
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        BaseData *base_cgns = Str_CGNS_Data->GetBaseData(iZone);

        int nCoords =base_cgns->GetNCoords();
        int idim = base_cgns->GetIdim();
        int jdim = base_cgns->GetJdim();
        int kdim = base_cgns->GetKdim();

        RDouble ***x = base_cgns->GetX();
        RDouble ***y = base_cgns->GetY();
        RDouble ***z = base_cgns->GetZ();

        StructGrid *grid = StructGridCast(grids[iZone]);

        grid->SetNI(idim);
        grid->SetNJ(jdim);
        grid->SetNK(kdim);

        grid->SetBasicDimension();

        int nTotalNode = grid->GetNTotalNode();
        RDouble *xCoor = new RDouble[nTotalNode];
        RDouble *yCoor = new RDouble[nTotalNode];
        RDouble *zCoor = new RDouble[nTotalNode];
        grid->SetX(xCoor);
        grid->SetY(yCoor);
        grid->SetZ(zCoor);
        grid->SetArrayLayout();

        RDouble3D &xx = * grid->GetStructX();
        RDouble3D &yy = * grid->GetStructY();
        RDouble3D &zz = * grid->GetStructZ();

        for (int k = 0; k < kdim; ++ k)
        {
            for (int j = 0; j < jdim; ++ j)
            {
                for (int i = 0; i < idim; ++ i)
                {
                    xx(i+1, j+1, k+1) = x[k][j][i];
                    yy(i+1, j+1, k+1) = y[k][j][i];

                    if (nCoords == 2) z[k][j][i] = 0.0;
                    zz(i+1, j+1, k+1) = z[k][j][i];
                }
            }
        }

        grid->ComputeMinMaxBox();
    }
}
void Pre_CGNSToPHengLEI_Struct::SetVolumeInfo()
{
    for (int iZone = 0; iZone < nBlocks; ++iZone)
    {
        BaseData *base_cgns = Str_CGNS_Data->GetBaseData(iZone);
        StructGrid *grid = StructGridCast(grids[iZone]);

        SimpleVC *volumeCondition = new SimpleVC();
        string vcName = base_cgns->GetVCName();
        volumeCondition->SetVCName(vcName);
        int vcType = base_cgns->GetVCType();
        volumeCondition->SetVCType(vcType);
        grid->SetVolumeCondition(volumeCondition);
    }

}
void Pre_CGNSToPHengLEI_Struct::SetBCInfo()
{
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        BaseData *base_cgns = Str_CGNS_Data->GetBaseData(iZone);

        int    nCoords       = base_cgns->GetNCoords();
        int    n1to1         = base_cgns->GetN1to1();
        int    nBCRegions    = base_cgns->GetNBCRegions();
        int    *lab1         = base_cgns->GetLab1();
        string *bcName       = base_cgns->GetBcName();
        cgsize_t_60000 *pnts = base_cgns->GetPnts();

        StructGrid *grid = StructGridCast(grids[iZone]);
        grid->CreateCompositeBCRegion(nBCRegions + n1to1);
        StructBCSet *structBCSet = grid->GetStructBCSet();

        for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
        {
            StructBC *BCRegions = new StructBC(iZone, iBCRegion);
            structBCSet->SetBCRegion(iBCRegion, BCRegions);

            string bcNameCurrent;
            int iMin, iMax, jMin, jMax, kMin, kMax, BCType;

            if (nCoords == 2)
            {
                iMin = static_cast<int>(pnts[iBCRegion][0]);
                iMax = static_cast<int>(pnts[iBCRegion][2]);

                jMin = static_cast<int>(pnts[iBCRegion][1]);
                jMax = static_cast<int>(pnts[iBCRegion][3]);

                kMin = 1;
                kMax = 1;
            }
            else
            {
                iMin = static_cast<int>(pnts[iBCRegion][0]);
                iMax = static_cast<int>(pnts[iBCRegion][3]);

                jMin = static_cast<int>(pnts[iBCRegion][1]);
                jMax = static_cast<int>(pnts[iBCRegion][4]);

                kMin = static_cast<int>(pnts[iBCRegion][2]);
                kMax = static_cast<int>(pnts[iBCRegion][5]);
            }

            BCType = lab1[iBCRegion];
            bcNameCurrent = bcName[iBCRegion];

            BCRegions->SetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);
            BCRegions->SetBCType(BCType);
            BCRegions->SetBCName(bcNameCurrent);
        }

        int ier = base_cgns->GetIer();
        if (ier == CG_OK)
        {
            n1to1 = base_cgns->GetN1to1();
            for (int one21 = 0; one21 < n1to1; ++ one21)
            {
                int *b = base_cgns->GetB(one21);
                int *vt = base_cgns->GetVt(one21);
                cgsize_t *sranges      = base_cgns->GetSrange(one21);
                cgsize_t *donor_ranges = base_cgns->GetDonor_range(one21);

                StructBC *BCRegions = new StructBC(iZone, nBCRegions + one21);
                structBCSet->SetBCRegion(nBCRegions + one21, BCRegions);

                string bcNameCurrent;
                int iMin, iMax, jMin, jMax, kMin, kMax, BCType;

                iMin = static_cast<int>(sranges[0]) * b[0];
                iMax = static_cast<int>(sranges[0 + nCoords]) * b[0];
                jMin = static_cast<int>(sranges[1]) * b[1];
                jMax = static_cast<int>(sranges[1 + nCoords]) * b[1];

                if (nCoords == 3)
                {
                    kMin = static_cast<int>(sranges[2]) * b[2];
                    kMax = static_cast<int>(sranges[2 + nCoords]) * b[2];
                }
                else
                {
                    kMin = 1;
                    kMax = 1;
                }

                BCType = -1;
                bcNameCurrent = "Connection";

                BCRegions->SetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);
                BCRegions->SetBCType(BCType);
                BCRegions->SetBCName(bcNameCurrent);

                iMin = static_cast<int>(donor_ranges[0]) * vt[0];
                iMax = static_cast<int>(donor_ranges[0 + nCoords]) * vt[0];
                jMin = static_cast<int>(donor_ranges[1]) * vt[1];
                jMax = static_cast<int>(donor_ranges[1 + nCoords]) * vt[1];

                if (nCoords == 3)
                {
                    kMin = static_cast<int>(donor_ranges[2]) * vt[2];
                    kMax = static_cast<int>(donor_ranges[2 + nCoords]) * vt[2];
                }
                else
                {
                    kMin = 1;
                    kMax = 1;
                }

                int nbt = base_cgns->Getlab2(one21) - 1;

                BCRegions->SetTargetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);
                BCRegions->SetTargetRegionBlock(nbt);

                if (nCoords == THREE_D)
                {
                    BCRegions->ComputeRelativeParameters();
                }
            }
        }
    }
}

void Pre_CGNSToPHengLEI_Struct::CheckMeshMultigrid()
{
    int km, NN;
    int km_grid = 1024;

    for (int iZone = 0; iZone < nBlocks; iZone ++)
    {
        StructGrid *grid = StructGridCast(grids[iZone]);

        StructBCSet *structBCSet = grid->GetStructBCSet();

        int nBCRegion = structBCSet->GetnBCRegion();
        for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion ++)
        {
            km = 1;
            NN = 2;

            StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);
            int *s_st = structBC->GetStartPoint();
            int *s_ed = structBC->GetEndPoint();

            while ((((abs(s_st[0])-1)%NN)==0) && (((abs(s_ed[0])-1)%NN)==0) && (((abs(s_st[1])-1)%NN)==0) && (((abs(s_ed[1])-1)%NN)==0) && (((abs(s_st[2])-1)%NN)==0) && (((abs(s_ed[2])-1)%NN)==0))
            {
                NN = NN * 2;
                km = km + 1;
            }
            km_grid = min(km_grid, km);
        }
    }

    PrintToWindow("The most valid multi-grid level is: ", km_grid, "\n");
}

}