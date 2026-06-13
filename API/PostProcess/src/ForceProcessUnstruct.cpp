#include "ForceProcessUnstruct.h"
#include "Post_ForceMoment.h"
#include "GradientOperation.h"
#include "Geo_UnstructGrid.h"
#include "PHMpi.h"
#include "TK_Exit.h"

namespace PHSPACE
{

ForceProcessUnstruct::ForceProcessUnstruct()
{
    CollectData = new DataContainer();
}

ForceProcessUnstruct::~ForceProcessUnstruct()
{
    delete CollectData;
}

void ForceProcessUnstruct::Run()
{
    InitValueField();
    CompValueField();
    InitForce();
    CompForce();
    CollectionForce();
    DumpForce();
}

void ForceProcessUnstruct::InitValueField()
{
    using namespace PHMPI;

    int nLocalZones = GetNumberofLocalZones();
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        UnstructGrid *grid = UnstructGridCast(GetGrid(zoneID));
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal = nTotalCell + nBoundFace;

        int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
        RDouble **t = NewPointer2<RDouble>(nTemperatureModel, nTotal);
        grid->UpdateDataPtr("t", t);

        RDouble *visl = new RDouble [nTotal];
        grid->UpdateDataPtr("visl",visl);
    }
}

void ForceProcessUnstruct::InitForce()
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

void ForceProcessUnstruct::CollectionForce()
{
    globalAerodynamicForce.CollectionForce(CollectData);
}

void ForceProcessUnstruct::AirForceCoef(int iZone)
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

void ForceProcessUnstruct::AirForceCoefParts(int iZone, int partID, int Coordinate)
{
    using namespace PHMPI;

    UnstructGrid * grid = UnstructGridCast(GetGrid(iZone));

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *visl = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));

    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    RDouble oRefReNumber = 1.0 / refReNumber;

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");

    const int ALL_PART = -1;
    SimpleBC * boundaryCondition = NULL;
    string partBCName = "";
    RDouble TorqueRefX = 0.0, TorqueRefY = 0.0, TorqueRefZ = 0.0;

    //! Used to compute hinge moment:
    int dumpHingeMoment = 0;
    RDouble localCoordAxis0[3];
    RDouble localCoordAxis1[3];
    RDouble hingeMoment = 0.0;

    if(partID == ALL_PART)
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

    int vis_run = 0;
    if (viscousType != INVISCID)
    {
        vis_run = 1;
    }

    using namespace IDX;
    RDouble *primitiveVarFarfield = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("prim_inf"));
    RDouble poo = primitiveVarFarfield[IP];

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea() ;

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int * leftCellofFace = grid->GetLeftCellOfFace();
    int * rightCellofFace = grid->GetRightCellOfFace();

    RDouble cpx,cpy,cpz;
    cpx = 0.0;
    cpy = 0.0;
    cpz = 0.0;

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

    RDouble pw,cp;
    RDouble dfx,dfy,dfz,dpx,dpy,dpz,fvsx,fvsy,fvsz;
    RDouble nx,ny,nz;
    RDouble vis,txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz;
    RDouble dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz;
    RDouble xc,yc,zc,dx,dy,dz,ods;

    RDouble **gradPrimtiveVarX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradPrimtiveVarY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradPrimtiveVarZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    using namespace IDX;

    // get the bcRegion information of unstruct grid.
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();


    // cycle all bcregions and find the solid surface type of boundary.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        string bcName = bcRegion->GetBCName();
        int bcType = bcRegion->GetBCType();

        // if the bc type is wall, calculate the wall coefficient.
        if(bcType == PHENGLEI::SOLID_SURFACE)
        {
            if (partID != ALL_PART && bcName != partBCName)
            {
                continue;
            }
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            int iFace, le, re;
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                iFace = *iter;
                le = leftCellofFace [ iFace ];
                re = rightCellofFace[ iFace ];

                nx  = xfn[iFace] * area[iFace];
                ny  = yfn[iFace] * area[iFace];
                nz  = zfn[iFace] * area[iFace];

                dx  = xfc[iFace] - xcc[le];
                dy  = yfc[iFace] - ycc[le];
                dz  = zfc[iFace] - zcc[le];

                //pressure drag
                pw = q[IP][le];
                cp    = two * (pw - poo);

                dpx = nx * cp;
                dpy = ny * cp;
                dpz = nz * cp;

                CA_p += dpx;
                CN_p += dpy;
                CZ_p += dpz;

                dfx = dpx;
                dfy = dpy;
                dfz = dpz;

                xc = xfc[iFace];
                yc = yfc[iFace];
                zc = zfc[iFace];

                if (viscousType > INVISCID)
                {
                    dudx = gradPrimtiveVarX[IU][le];
                    dudy = gradPrimtiveVarY[IU][le];
                    dudz = gradPrimtiveVarZ[IU][le];

                    dvdx = gradPrimtiveVarX[IV][le];
                    dvdy = gradPrimtiveVarY[IV][le];
                    dvdz = gradPrimtiveVarZ[IV][le];

                    dwdx = gradPrimtiveVarX[IW][le];
                    dwdy = gradPrimtiveVarY[IW][le];
                    dwdz = gradPrimtiveVarZ[IW][le];

                    //! Gradient correction.
                    dx  = xcc[re] - xcc[le];
                    dy  = ycc[re] - ycc[le];
                    dz  = zcc[re] - zcc[le];
                    ods = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
                    dx  = dx * ods;
                    dy  = dy * ods;
                    dz  = dz * ods;

                    CorrectGradient(q[IU][le],q[IU][re],dudx,dudy,dudz,dx,dy,dz,ods);
                    CorrectGradient(q[IV][le],q[IV][re],dvdx,dvdy,dvdz,dx,dy,dz,ods);
                    CorrectGradient(q[IW][le],q[IW][re],dwdx,dwdy,dwdz,dx,dy,dz,ods);

                    //! It is unnecessary to consider because the turbulence viscosity coefficient is 0.
                    vis = visl[le];

                    const RDouble c23 = 2.0 / 3.0;
                    txx = vis * c23 * (2.0 * dudx - dvdy - dwdz);
                    tyy = vis * c23 * (2.0 * dvdy - dwdz - dudx);
                    tzz = vis * c23 * (2.0 * dwdz - dudx - dvdy);
                    txy = vis * (dudy + dvdx);
                    txz = vis * (dudz + dwdx);
                    tyz = vis * (dvdz + dwdy);
                    tyx = txy;
                    tzx = txz;
                    tzy = tyz;

                    fvsx = - two * (nx * txx + ny * tyx + nz * tzx) * oRefReNumber;
                    fvsy = - two * (nx * txy + ny * tyy + nz * tzy) * oRefReNumber;
                    fvsz = - two * (nx * txz + ny * tyz + nz * tzz) * oRefReNumber;

                    CA_f += fvsx;
                    CN_f += fvsy;
                    CZ_f += fvsz;

                    dfx += fvsx;
                    dfy += fvsy;
                    dfz += fvsz;

                    Cl_f += (yc - TorqueRefY) * fvsz - (zc - TorqueRefZ) * fvsy;
                    Cn_f += (zc - TorqueRefZ) * fvsx - (xc - TorqueRefX) * fvsz;
                    Cm_f += (xc - TorqueRefX) * fvsy - (yc - TorqueRefY) * fvsx;
                }

                Cl_p += (yc - TorqueRefY) * dpz - (zc - TorqueRefZ) * dpy;
                Cn_p += (zc - TorqueRefZ) * dpx - (xc - TorqueRefX) * dpz;
                Cm_p += (xc - TorqueRefX) * dpy - (yc - TorqueRefY) * dpx;

                cpx += dpx;
                cpy += dpy;
                cpz += dpz;

                if(dumpHingeMoment)
                {
                    Post_ForceMoment forceMoment;
                    RDouble point[3] = {xfc[iFace], yfc[iFace], zfc[iFace]};
                    RDouble faceForce[3] = {dfx, dfy, dfz};
                    hingeMoment += forceMoment.ComputeMoment(point, localCoordAxis0, localCoordAxis1, faceForce);
                }
            }
        }
    }

    this->CollectData->Write(&CA_f,sizeof(RDouble));
    this->CollectData->Write(&CA_p,sizeof(RDouble));
    this->CollectData->Write(&CN_f,sizeof(RDouble));
    this->CollectData->Write(&CN_p,sizeof(RDouble));
    this->CollectData->Write(&CZ_f,sizeof(RDouble));
    this->CollectData->Write(&CZ_p,sizeof(RDouble));

    this->CollectData->Write(&cpx,sizeof(RDouble));
    this->CollectData->Write(&cpy,sizeof(RDouble));
    this->CollectData->Write(&cpz,sizeof(RDouble));

    this->CollectData->Write(&Cl_f,sizeof(RDouble));
    this->CollectData->Write(&Cl_p,sizeof(RDouble));
    this->CollectData->Write(&Cn_f,sizeof(RDouble));
    this->CollectData->Write(&Cn_p,sizeof(RDouble));
    this->CollectData->Write(&Cm_f,sizeof(RDouble));
    this->CollectData->Write(&Cm_p,sizeof(RDouble));

    this->CollectData->Write(&hingeMoment, sizeof(RDouble));
}

void ForceProcessUnstruct::CompForce()
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

void ForceProcessUnstruct::DumpForce()
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

void ForceProcessUnstruct::CompValueField()
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
        GetGradientField(grid);
    }
}

void ForceProcessUnstruct::ComputeGamaAndTemperature(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **t   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("t") );
    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    using namespace IDX;

    RDouble coefficientOfStateEquation = GlobalDataBase::GetDoubleParaFromDB("coefficientOfStateEquation");

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        RDouble &rm = q[IR][iCell];
        RDouble &pm = q[IP][iCell];
        t[ITT][iCell] = pm / (coefficientOfStateEquation * rm);
    }
}

void ForceProcessUnstruct::GetGradientField(Grid *gridIn)
{
    UnstructGrid *gridUnstruct = UnstructGridCast(gridIn);

    int nChemical           = GlobalDataBase::GetIntParaFromDB("nchem");
    int nNSEquation         = GlobalDataBase::GetIntParaFromDB("nm");
    int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    int nEquation = nNSEquation + nChemical + nTemperatureModel - 1;

    RDouble **gradPrimtiveVarX = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradPrimtiveVarY = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradPrimtiveVarZ = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarZ"));
    RDouble **q     = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("q"));

    string gradientName = GlobalDataBase::GetStrParaFromDB("gradientName");

    for (int m = 0; m < nEquation; ++ m)
    {
        if (gradientName == "ggnode" || gradientName == "ggnodelaplacian")
        {
            gridUnstruct->CompGradientGGCell(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "ggcell")
        {
            gridUnstruct->CompGradientGGCell(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "lsq")
        {
            gridUnstruct->CompGradientLSQ(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }    
        else if (gradientName == "ggnode_weight")
        {
            gridUnstruct->CompGradientGGNodeWeight(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "gg_m2")
        {
            gridUnstruct->CompGradientGGModified2(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "ggcellnew")
        {
            gridUnstruct->CompGradientGGCellNew(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "ggcellw")
        {
            gridUnstruct->CompGradientGGCellW(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else
        {
            TK_Exit::ExceptionExit("No reconstruction method has been choosed ! /n");
        }
    }
}

void ForceProcessUnstruct::ComputeViscousLaminarCoefficient(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble **t = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));

    RDouble tsuth;
    GlobalDataBase::GetData("tsuth", &tsuth, PHDOUBLE, 1);

    using namespace IDX;

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        viscousLaminar[iCell] = t[ITT][iCell] * sqrt(t[ITT][iCell]) * (1.0 + tsuth) / (t[ITT][iCell] + tsuth);
    }
}

}