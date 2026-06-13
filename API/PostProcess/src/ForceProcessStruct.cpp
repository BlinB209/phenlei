#include "ForceProcessStruct.h"
#include "Geo_StructBC.h"
#include "Glb_Dimension.h"
#include "GradientOperation.h"
#include "Post_ForceMoment.h"
#include <cmath>
#include "Geo_StructGrid.h"

namespace PHSPACE
{

ForceProcessStruct::ForceProcessStruct()
{
    CollectData = new DataContainer();
}

ForceProcessStruct::~ForceProcessStruct()
{
    delete CollectData;
}

void ForceProcessStruct::Run()
{
    CompValueField();
    InitForce();
    CompForce();
    CollectionForce();
    DumpForce();
}

void ForceProcessStruct::InitForce()
{
    RDouble attack, sideslip;
    RDouble attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    RDouble angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");

    attack   = attackd * PI / 180.0;
    sideslip = angleSlide * PI / 180.0;

    RDouble forceReferenceArea   = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
    RDouble forceReferenceLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
    RDouble forceReferenceLengthSpanWise = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");

    globalAerodynamicForce.Init(attack, sideslip, forceReferenceArea, forceReferenceLength, forceReferenceLengthSpanWise);
}

void ForceProcessStruct::CollectionForce()
{
    globalAerodynamicForce.CollectionForce(CollectData);
}

void ForceProcessStruct::AirForceCoef(int iZone)
{
    const int ALL_PART = -1;
    AirForceCoefParts(iZone, ALL_PART);

    uint_t nTBC = GlobalBoundaryCondition::GetNumberOfBC();
    for (int iPart = 0; iPart < nTBC; ++ iPart)
    {
        if (!GlobalBoundaryCondition::IsSolidWallBC(iPart)) continue;

        AirForceCoefParts(iZone, iPart, LocalCoordinate);

        AirForceCoefParts(iZone, iPart);
    }
}

void ForceProcessStruct::AirForceCoefParts(int iZone, int partID, int Coordinate)
{
    StructGrid *grid = StructGridCast(GetGrid(iZone));

    RDouble4D &primitiveVars = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &viscousLaminar = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D &viscousTurbulence = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("vist"));

    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    RDouble oRefReNumber = 1.0 / refReNumber;

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");

    const int ALL_PART = -1;
    SimpleBC *boundaryCondition = NULL;
    string partBCName = "";
    RDouble TorqueRefX = 0.0, TorqueRefY = 0.0, TorqueRefZ = 0.0;

    //! Compute hinge moment.
    int dumpHingeMoment = 0;
    RDouble localCoordAxis0[3] = { 0 };
    RDouble localCoordAxis1[3] = { 0 };

    if (partID == ALL_PART)
    {
        TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
        TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
        TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");
    }
    else if ((partID != ALL_PART) && (Coordinate == GlobalCoordinate))
    {
        boundaryCondition = GlobalBoundaryCondition::GetBC(partID);
        partBCName = boundaryCondition->GetBCName();

        TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
        TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
        TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");
    }
    else if ((partID != ALL_PART) && (Coordinate == LocalCoordinate))
    {
        boundaryCondition = GlobalBoundaryCondition::GetBC(partID);
        partBCName = boundaryCondition->GetBCName();
        boundaryCondition->GetParamData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
        boundaryCondition->GetParamData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
        boundaryCondition->GetParamData("TorqueRefZ", &TorqueRefZ, PHDOUBLE, 1);

        boundaryCondition->GetParamData("dumpHingeMoment", &dumpHingeMoment, PHINT, 1);
        boundaryCondition->GetParamData("localCoordAxis0", localCoordAxis0, PHDOUBLE, 3);
        boundaryCondition->GetParamData("localCoordAxis1", localCoordAxis1, PHDOUBLE, 3);
    }

    using namespace IDX;
    RDouble *primitiveVarFarfield = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("prim_inf"));
    RDouble pressureFarfield = primitiveVarFarfield[IP];

    RDouble4D &xfv = *(grid->GetFaceVectorX());
    RDouble4D &yfv = *(grid->GetFaceVectorY());
    RDouble4D &zfv = *(grid->GetFaceVectorZ());
    RDouble3D &vol = *(grid->GetCellVolume());

    RDouble cpx = 0.0;
    RDouble cpy = 0.0;
    RDouble cpz = 0.0;
    RDouble hingeMoment = 0.0;

    RDouble CA_f = 0.0;
    RDouble CA_p = 0.0;
    RDouble CN_f = 0.0;
    RDouble CN_p = 0.0;
    RDouble CZ_f = 0.0;
    RDouble CZ_p = 0.0;
    RDouble Cl_f = 0.0;
    RDouble Cl_p = 0.0;
    RDouble Cn_f = 0.0;
    RDouble Cn_p = 0.0;
    RDouble Cm_f = 0.0;
    RDouble Cm_p = 0.0;

    int ndim = GetDim();

    using namespace IDX;

    int ix = 0;

    StructBCSet *structBCSet = grid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();
    
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcRegion = structBCSet->GetBCRegion(iBCRegion);
        int BCType = bcRegion->GetBCType();

        int ist = 0, ied = 0, jst = 0, jed = 0, kst = 0, ked = 0;
        bcRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        bool isNeedComputeAirForce = (BCType == PHENGLEI::SOLID_SURFACE || BCType == PHENGLEI::ABLATION_SURFACE);
        if (isNeedComputeAirForce && partID != ALL_PART)
        {
            const string &bcName = bcRegion->GetBCName();
            isNeedComputeAirForce = isNeedComputeAirForce && (bcName == partBCName);
        }

        if (!isNeedComputeAirForce)
        {
            ix += (ied - ist + 1) * (jed - jst + 1) * (ked - kst + 1);
            continue;
        }

        int *s_lr3d = bcRegion->GetFaceDirectionIndex();
        int nSurface = bcRegion->GetFaceDirection() + 1;
        int leftOrRightIndex = bcRegion->GetFaceLeftOrRightIndex();

        int iWall = 0, jWall = 0, kWall = 0;

        RDouble dfx = 0, dfy = 0, dfz = 0, dpx = 0, dpy = 0, dpz = 0;

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    //! i, j, k - is first cell index near SOLID_SURFACE.
                    bcRegion->GetBoundaryFaceIndex(i, j, k, iWall, jWall, kWall);

                    //! The surface direction of SOLID_SURFACE.
                    RDouble nx = leftOrRightIndex * xfv(iWall, jWall, kWall, nSurface);
                    RDouble ny = leftOrRightIndex * yfv(iWall, jWall, kWall, nSurface);
                    RDouble nz = leftOrRightIndex * zfv(iWall, jWall, kWall, nSurface);

                    //! The index of  first ghost cell layer near SOLID_SURFACE.
                    int iCellGhostLayer = i + s_lr3d[0];
                    int jCellGhostLayer = j + s_lr3d[1];
                    int kCellGhostLayer = k + s_lr3d[2];

                    RDouble pCellfirstLayer = primitiveVars(i, j, k, IP);

                    RDouble pressureWall = pCellfirstLayer;
                    RDouble pressureCoef = two * (pressureWall - pressureFarfield);

                    dpx = nx * pressureCoef;
                    dpy = ny * pressureCoef;
                    dpz = nz * pressureCoef;

                    CA_p += dpx;
                    CN_p += dpy;
                    CZ_p += dpz;

                    dfx = dpx;
                    dfy = dpy;
                    dfz = dpz;

                    RDouble xc = 0, yc = 0, zc = 0;
                    grid->FaceCoor(iWall, jWall, kWall, nSurface, xc, yc, zc);

                    RDouble visLaminar = 0, visTurbulence = 0;
                    RDouble dudx = 0, dudy = 0, dudz = 0, dvdx = 0, dvdy = 0, dvdz = 0, dwdx = 0, dwdy = 0, dwdz = 0;
                    if (viscousType > INVISCID)
                    {
                        Geometry3D *geometry3D = new Geometry3D;

                        if (nSurface == 1)
                        {
                            //! Get i-th face vector on the wall surface.
                            Get_GEO_I(i, j, k, iWall, iWall, ndim, geometry3D, xfv, yfv, zfv, vol);
                            //! Get derivative to the wall in i-direction of u,v,w.
                            DXDYDZ_I(primitiveVars, geometry3D, i, j, k, ndim, IU, dudx, dudy, dudz, leftOrRightIndex);
                            DXDYDZ_I(primitiveVars, geometry3D, i, j, k, ndim, IV, dvdx, dvdy, dvdz, leftOrRightIndex);
                            DXDYDZ_I(primitiveVars, geometry3D, i, j, k, ndim, IW, dwdx, dwdy, dwdz, leftOrRightIndex);
                            visLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(iCellGhostLayer, j, k));
                            visTurbulence = half * (viscousTurbulence(i, j, k) + viscousTurbulence(iCellGhostLayer, j, k));
                        }
                        else if (nSurface == 2)
                        {
                            //! Get j-th face vector on the wall surface.
                            Get_GEO_J(i, j, k, jWall, jWall, ndim, geometry3D, xfv, yfv, zfv, vol);
                            //! Get derivative to the wall in j-direction of u,v,w.
                            DXDYDZ_J(primitiveVars, geometry3D, i, j, k, ndim, IU, dudx, dudy, dudz, leftOrRightIndex);
                            DXDYDZ_J(primitiveVars, geometry3D, i, j, k, ndim, IV, dvdx, dvdy, dvdz, leftOrRightIndex);
                            DXDYDZ_J(primitiveVars, geometry3D, i, j, k, ndim, IW, dwdx, dwdy, dwdz, leftOrRightIndex);
                            visLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(i, jCellGhostLayer, k));
                            visTurbulence = half * (viscousTurbulence(i, j, k) + viscousTurbulence(i, jCellGhostLayer, k));
                        }
                        else
                        {
                            //! Get k-th face vector on the wall surface.
                            Get_GEO_K(i, j, k, kWall, kWall, geometry3D, xfv, yfv, zfv, vol);
                            //! Get derivative to the wall in k-direction of u,v,w.
                            DXDYDZ_K(primitiveVars, geometry3D, i, j, k, IU, dudx, dudy, dudz, leftOrRightIndex);
                            DXDYDZ_K(primitiveVars, geometry3D, i, j, k, IV, dvdx, dvdy, dvdz, leftOrRightIndex);
                            DXDYDZ_K(primitiveVars, geometry3D, i, j, k, IW, dwdx, dwdy, dwdz, leftOrRightIndex);
                            visLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(i, j, kCellGhostLayer));
                            visTurbulence = half * (viscousTurbulence(i, j, k) + viscousTurbulence(i, j, kCellGhostLayer));
                        }

                        //! It is unnecessary to consider because the turbulence viscosity coefficient is 0.
                        //visTurbulence = 0.0;
                        RDouble viscousCoef = visLaminar + visTurbulence;

                        RDouble divv2p3 = two3rd * (dudx + dvdy + dwdz);

                        //! Stress components.
                        RDouble txx = viscousCoef * (two * dudx - divv2p3);
                        RDouble tyy = viscousCoef * (two * dvdy - divv2p3);
                        RDouble tzz = viscousCoef * (two * dwdz - divv2p3);
                        RDouble txy = viscousCoef * (dudy + dvdx);
                        RDouble txz = viscousCoef * (dudz + dwdx);
                        RDouble tyz = viscousCoef * (dvdz + dwdy);

                        RDouble fvsx = -two * (nx * txx + ny * txy + nz * txz) * oRefReNumber;
                        RDouble fvsy = -two * (nx * txy + ny * tyy + nz * tyz) * oRefReNumber;
                        RDouble fvsz = -two * (nx * txz + ny * tyz + nz * tzz) * oRefReNumber;

                        dfx += fvsx;
                        dfy += fvsy;
                        dfz += fvsz;

                        CA_f += fvsx;
                        CN_f += fvsy;
                        CZ_f += fvsz;

                        Cl_f += (yc - TorqueRefY) * fvsz - (zc - TorqueRefZ) * fvsy;
                        Cn_f += (zc - TorqueRefZ) * fvsx - (xc - TorqueRefX) * fvsz;
                        Cm_f += (xc - TorqueRefX) * fvsy - (yc - TorqueRefY) * fvsx;

                        delete geometry3D;
                    }

                    Cl_p += (yc - TorqueRefY) * dpz - (zc - TorqueRefZ) * dpy;
                    Cn_p += (zc - TorqueRefZ) * dpx - (xc - TorqueRefX) * dpz;
                    Cm_p += (xc - TorqueRefX) * dpy - (yc - TorqueRefY) * dpx;

                    cpx += dpx;
                    cpy += dpy;
                    cpz += dpz;

                    if (dumpHingeMoment)
                    {
                        Post_ForceMoment forceMoment;
                        RDouble point[3] = { xc, yc, zc };
                        RDouble faceForce[3] = { dfx, dfy, dfz };
                        hingeMoment += forceMoment.ComputeMoment(point, localCoordAxis0, localCoordAxis1, faceForce);
                    }

                    ++ ix;
                }
            }
        }
    }

    this->CollectData->Write(&CA_f, sizeof(RDouble));
    this->CollectData->Write(&CA_p, sizeof(RDouble));
    this->CollectData->Write(&CN_f, sizeof(RDouble));
    this->CollectData->Write(&CN_p, sizeof(RDouble));
    this->CollectData->Write(&CZ_f, sizeof(RDouble));
    this->CollectData->Write(&CZ_p, sizeof(RDouble));

    this->CollectData->Write(&cpx, sizeof(RDouble));
    this->CollectData->Write(&cpy, sizeof(RDouble));
    this->CollectData->Write(&cpz, sizeof(RDouble));

    this->CollectData->Write(&Cl_f, sizeof(RDouble));
    this->CollectData->Write(&Cl_p, sizeof(RDouble));
    this->CollectData->Write(&Cn_f, sizeof(RDouble));
    this->CollectData->Write(&Cn_p, sizeof(RDouble));
    this->CollectData->Write(&Cm_f, sizeof(RDouble));
    this->CollectData->Write(&Cm_p, sizeof(RDouble));

    this->CollectData->Write(&hingeMoment, sizeof(RDouble));
}

void ForceProcessStruct::CompForce()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int zone_proc = GetZoneProcessorID(iZone);
        int myid = GetCurrentProcessorID();

        if (myid == zone_proc)
        {
            AirForceCoef(iZone);
        }
    }
}

void ForceProcessStruct::DumpForce()
{
    ActionKey *actkeyDumpAirForceCoef = new ActionKey();
    delete actkeyDumpAirForceCoef->data;
    actkeyDumpAirForceCoef->data = this->CollectData;

    string airCoefFile = "aircoef.dat";
    GlobalDataBase::GetData("aircoeffile", &airCoefFile, PHSTRING, 1);

    actkeyDumpAirForceCoef->filename = airCoefFile;
    ios_base::openmode openmode = ios_base::out|ios_base::app;
    actkeyDumpAirForceCoef->openmode = openmode;

    int iterationStep;
    GlobalDataBase::GetData("outnstep", &iterationStep, PHINT, 1);

    globalAerodynamicForce.DumpForce(actkeyDumpAirForceCoef, iterationStep);

    DataContainer *dataNew = new DataContainer();
    actkeyDumpAirForceCoef->data = dataNew;

    delete actkeyDumpAirForceCoef;
}

void ForceProcessStruct::CompValueField()
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    #ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Grid *grid = GetGrid(zoneID,0);

        ComputeGamaAndTemperature(grid);
        ComputeViscousLaminarCoefficient(grid);
    }
}

void ForceProcessStruct::ComputeGamaAndTemperature(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &q = * reinterpret_cast<RDouble4D *>(grid->GetDataPtr("q"));
    RDouble4D &t = * reinterpret_cast<RDouble4D *>(grid->GetDataPtr("t"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    using namespace IDX;

    RDouble omav = one;
    RDouble coefficientOfStateEquation = GlobalDataBase::GetDoubleParaFromDB("coefficientOfStateEquation");

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble &rm = q(i, j, k, IR);
                RDouble &pm = q(i, j, k, IP);
                t(i, j, k, ITT)  = pm / (coefficientOfStateEquation * rm * omav);
            }
        }
    }
}

void ForceProcessStruct::ComputeViscousLaminarCoefficient(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble3D &visl = * reinterpret_cast<RDouble3D *>(grid->GetDataPtr("visl"));
    RDouble4D &t    = * reinterpret_cast<RDouble4D *>(grid->GetDataPtr("t"));

    RDouble tsuth;
    GlobalDataBase::GetData("tsuth", &tsuth, PHDOUBLE, 1);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    using namespace IDX;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble tm = t(i, j, k, ITT);
                visl(i, j, k) = tm * sqrt(tm) * (1.0 + tsuth) / (tm + tsuth);
            }
        }
    }

    RDouble wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");

    if (wallTemperature <= 0.0)
    {
        GhostCell3D(visl, ni, nj, nk);
    }
    else
    {
        int idx_bc[2];
        idx_bc[0] = 0;

        idx_bc[1] = ni;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int m = 0; m <= 1; ++ m)
                {
                    int i = idx_bc[m];
                    RDouble tm = t(i, j, k, ITT);
                    visl(i, j, k) = tm * sqrt(tm) * (1.0 + tsuth) / (tm + tsuth);
                }
            }
        }

        idx_bc[1] = nj;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m <= 1; ++ m)
                {
                    int j = idx_bc[m];
                    RDouble tm = t(i, j, k, ITT);
                    visl(i, j, k) = tm * sqrt(tm) * (1.0 + tsuth) / (tm + tsuth);
                }
            }
        }

        if (GetDim() == THREE_D)
        {
            idx_bc[1] = nk;
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    for (int m = 0; m <= 1; ++ m)
                    {
                        int k = idx_bc[m];
                        RDouble tm = t(i, j, k, ITT);
                        visl(i, j, k) = tm * sqrt(tm) * (1.0 + tsuth) / (tm + tsuth);
                    }
                }
            }
        }
    }
}

}