#include "Force.h"
#include "Math_BasisFunction.h"
#include "Constants.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "Geo_SimpleBC.h"
#include "TK_Time.h"
#include "TK_Parse.h"
#include "TK_Exit.h"
#include "Glb_Dimension.h"
#include "GlobalDataBase.h"
#include "PHHeader.h"
#include "PHMpi.h"
#include "AleManager.h"

using namespace std;

namespace PHSPACE
{

Force globalAerodynamicForce;

void InitForce()
{

    RDouble attack, sideslip;
    RDouble attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    RDouble angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");

    attack = attackd * PI / 180.0;
    sideslip = angleSlide * PI / 180.0;

    RDouble forceReferenceArea = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
    RDouble forceReferenceLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
    RDouble forceReferenceLengthSpanWise = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");

    globalAerodynamicForce.Init(attack, sideslip, forceReferenceArea, forceReferenceLength, forceReferenceLengthSpanWise);
}

void CollectionForce(ActionKey *actkey)
{
    globalAerodynamicForce.CollectionForce(actkey->GetData());
}

void DumpForce(ActionKey *actkey)
{
    int iterationStep;
    GlobalDataBase::GetData("outnstep", &iterationStep, PHINT, 1);

    globalAerodynamicForce.DumpForce(actkey, iterationStep);
}

Force::Force()
{
    partBCID = 0;

    cpx = 0;
    cpy = 0;
    cpz = 0;

    cl_tot = 0;
    cd_tot = 0;
    cd_pr = 0;
    cd_sf = 0;
    cd_cl2pa = 0;
    xcp = 0;
    ycp = 0;
    zcp = 0;
    side = 0;

    CA_f = 0;
    CA_p = 0;
    CN_f = 0;
    CN_p = 0;
    CZ_f = 0;
    CZ_p = 0;
    Cl_f = 0;
    Cl_p = 0;
    Cn_f = 0;
    Cn_p = 0;
    Cm_f = 0;
    Cm_p = 0;

    CA_f_body = 0;
    CN_f_body = 0;
    CZ_f_body = 0;
    Cl_f_body = 0;
    Cn_f_body = 0;
    Cm_f_body = 0;
    CA_p_body = 0;
    CN_p_body = 0;
    CZ_p_body = 0;
    Cl_p_body = 0;
    Cn_p_body = 0;
    Cm_p_body = 0;

    numberOfBody = 0;
    wallNameOfEachBody.resize(0);
}

Force::~Force()
{
    delete [] partBCID;    partBCID = NULL;

    delete [] cpx;    cpx = NULL;
    delete [] cpy;    cpy = NULL;
    delete [] cpz;    cpz = NULL;

    delete [] cl_tot;    cl_tot = NULL;
    delete [] cd_tot;    cd_tot = NULL;
    delete [] cd_pr;    cd_pr = NULL;
    delete [] cd_sf;    cd_sf = NULL;
    delete [] cd_cl2pa;    cd_cl2pa = NULL;
    delete [] xcp;    xcp = NULL;
    delete [] ycp;    ycp = NULL;
    delete [] zcp;    zcp = NULL;
    delete [] side;    side = NULL;

    delete [] CA_f;    CA_f = NULL;
    delete [] CA_p;    CA_p = NULL;
    delete [] CN_f;    CN_f = NULL;
    delete [] CN_p;    CN_p = NULL;
    delete [] CZ_f;    CZ_f = NULL;
    delete [] CZ_p;    CZ_p = NULL;
    delete [] Cl_f;    Cl_f = NULL;
    delete [] Cl_p;    Cl_p = NULL;
    delete [] Cn_f;    Cn_f = NULL;
    delete [] Cn_p;    Cn_p = NULL;
    delete [] Cm_f;    Cm_f = NULL;
    delete [] Cm_p;    Cm_p = NULL;

    delete [] hingeMoment;    hingeMoment = NULL;
   
    delete [] CA_f_body;    CA_f_body = NULL;
    delete [] CN_f_body;    CN_f_body = NULL;
    delete [] CZ_f_body;    CZ_f_body = NULL;
    delete [] Cl_f_body;    Cl_f_body = NULL;
    delete [] Cn_f_body;    Cn_f_body = NULL;
    delete [] Cm_f_body;    Cm_f_body = NULL;
    delete [] CA_p_body;    CA_p_body = NULL;
    delete [] CN_p_body;    CN_p_body = NULL;
    delete [] CZ_p_body;    CZ_p_body = NULL;
    delete [] Cl_p_body;    Cl_p_body = NULL;
    delete [] Cn_p_body;    Cn_p_body = NULL;
    delete [] Cm_p_body;    Cm_p_body = NULL;
}

void Force::Init(RDouble attack, RDouble sideslip, RDouble forceReferenceArea, RDouble forceReferenceLength, RDouble forceReferenceLengthSpanWise)
{
    this->attack = attack;
    this->sideslip = sideslip;
    this->forceReferenceArea = forceReferenceArea;
    this->forceReferenceLength = forceReferenceLength;
    this->forceReferenceLengthSpanWise = forceReferenceLengthSpanWise;

    isPartExist = (GlobalBoundaryCondition::GetGlobalBoundaryConditionList() != 0);
    numberofSolidWallPart = 0;
    if (isPartExist)
    {
        numberofSolidWallPart = GlobalBoundaryCondition::GetNumberOfSolidWallPart();
    }
    //! Because of global and local part force, numberofSolidWallPart * 2
    numberofTotalPart = numberofSolidWallPart * 2 + 1;

    int myid = PHMPI::GetCurrentProcessorID();
    int server = PHMPI::GetServerProcessorID();
    if ((myid == server) && (!cpx))
    {
        cpx = new RDouble[numberofTotalPart];
        cpy = new RDouble[numberofTotalPart];
        cpz = new RDouble[numberofTotalPart];
        hingeMoment = new RDouble[numberofTotalPart];

        CA_f = new RDouble[numberofTotalPart];
        CA_p = new RDouble[numberofTotalPart];
        CN_f = new RDouble[numberofTotalPart];
        CN_p = new RDouble[numberofTotalPart];
        CZ_f = new RDouble[numberofTotalPart];
        CZ_p = new RDouble[numberofTotalPart];
        Cl_f = new RDouble[numberofTotalPart];
        Cl_p = new RDouble[numberofTotalPart];
        Cn_f = new RDouble[numberofTotalPart];
        Cn_p = new RDouble[numberofTotalPart];
        Cm_f = new RDouble[numberofTotalPart];
        Cm_p = new RDouble[numberofTotalPart];

        PHSPACE::SetField(cpx, 0.0, numberofTotalPart);
        PHSPACE::SetField(cpy, 0.0, numberofTotalPart);
        PHSPACE::SetField(cpz, 0.0, numberofTotalPart);
        PHSPACE::SetField(hingeMoment, 0.0, numberofTotalPart);

        PHSPACE::SetField(CA_f, 0.0, numberofTotalPart);
        PHSPACE::SetField(CA_p, 0.0, numberofTotalPart);
        PHSPACE::SetField(CN_f, 0.0, numberofTotalPart);
        PHSPACE::SetField(CN_p, 0.0, numberofTotalPart);
        PHSPACE::SetField(CZ_f, 0.0, numberofTotalPart);
        PHSPACE::SetField(CZ_p, 0.0, numberofTotalPart);
        PHSPACE::SetField(Cl_f, 0.0, numberofTotalPart);
        PHSPACE::SetField(Cl_p, 0.0, numberofTotalPart);
        PHSPACE::SetField(Cn_f, 0.0, numberofTotalPart);
        PHSPACE::SetField(Cn_p, 0.0, numberofTotalPart);
        PHSPACE::SetField(Cm_f, 0.0, numberofTotalPart);
        PHSPACE::SetField(Cm_p, 0.0, numberofTotalPart);
    }

    numberOfBody = static_cast<int>(GlobalBoundaryCondition::GetNumberOfBody());
    if ((myid == server) && (!CA_f_body) && (numberOfBody))
    {
        CA_f_body = new RDouble[numberOfBody];
        CN_f_body = new RDouble[numberOfBody];
        CZ_f_body = new RDouble[numberOfBody];
        Cl_f_body = new RDouble[numberOfBody];
        Cn_f_body = new RDouble[numberOfBody];
        Cm_f_body = new RDouble[numberOfBody];
        CA_p_body = new RDouble[numberOfBody];
        CN_p_body = new RDouble[numberOfBody];
        CZ_p_body = new RDouble[numberOfBody];
        Cl_p_body = new RDouble[numberOfBody];
        Cn_p_body = new RDouble[numberOfBody];
        Cm_p_body = new RDouble[numberOfBody];

        set < string > *bodyList = GlobalBoundaryCondition::GetBodyNameList();
        set < string >::iterator bodyListIter;
        wallNameOfEachBody.resize(numberOfBody);
        int bodyIndex = 0;
        for (bodyListIter = bodyList->begin(); bodyListIter != bodyList->end(); ++ bodyListIter)
        {
            string bodyNameInList = *bodyListIter;
            wallNameOfEachBody[bodyIndex].push_back(bodyNameInList);
            bodyIndex ++;
        }
    }

    if ((myid == server) && (numberOfBody))
    {
        PHSPACE::SetField(CA_f_body, 0.0, numberOfBody);
        PHSPACE::SetField(CN_f_body, 0.0, numberOfBody);
        PHSPACE::SetField(CZ_f_body, 0.0, numberOfBody);
        PHSPACE::SetField(Cl_f_body, 0.0, numberOfBody);
        PHSPACE::SetField(Cn_f_body, 0.0, numberOfBody);
        PHSPACE::SetField(Cm_f_body, 0.0, numberOfBody);
        PHSPACE::SetField(CA_p_body, 0.0, numberOfBody);
        PHSPACE::SetField(CN_p_body, 0.0, numberOfBody);
        PHSPACE::SetField(CZ_p_body, 0.0, numberOfBody);
        PHSPACE::SetField(Cl_p_body, 0.0, numberOfBody);
        PHSPACE::SetField(Cn_p_body, 0.0, numberOfBody);
        PHSPACE::SetField(Cm_p_body, 0.0, numberOfBody);
    }

    if ((!partBCID) && (myid == server))
    {
        partBCID = new int[numberofTotalPart];

        const int ALL_SOLID_WALL = -1;
        PHSPACE::SetField(partBCID, ALL_SOLID_WALL, numberofTotalPart);

        int partCount = 1;
        uint_t nTBC = GlobalBoundaryCondition::GetNumberOfBC();
        for (int iBC = 0; iBC < nTBC; ++ iBC)
        {
            if ((!GlobalBoundaryCondition::IsSolidWallBC(iBC)) && ((GlobalBoundaryCondition::GetBC(iBC)->GetBCType() != PHENGLEI::ABLATION_SURFACE))) continue;

            //! Because of global and local part force, partBCID set twice
            partBCID[partCount ++] = iBC;
            partBCID[partCount ++] = iBC;

            if (numberOfBody == 0)
            {
                continue;
            }
            string wallName = GlobalBoundaryCondition::GetBoundaryName(iBC);
            int bodyIndex = static_cast<int>(GlobalBoundaryCondition::GetBodyIndex(iBC));
            wallNameOfEachBody[bodyIndex].push_back(wallName);
        }
    }
}

void Force::CollectionForce(DataContainer *data_in)
{
    data_in->MoveToBegin();
    int nLocalZones = PHMPI::GetNumberofLocalZones();

    int nTotalVarDim = nForceVar * numberofTotalPart;
    double *localForceVar = new double[nTotalVarDim];
    PHSPACE::SetField(localForceVar, 0.0, nTotalVarDim);

    int nTBC = static_cast<int>(GlobalBoundaryCondition::GetNumberOfBC());
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int varCount = 0;

        //! '-1'     , represent total wall force;\n
        //! 0 - nWALL, represent each wall part force.
        for (int iPart = ALL_SOLID_WALL; iPart < nTBC; ++ iPart)
        {
            bool ifRead = (iPart == ALL_SOLID_WALL); // total wall force.
            if (iPart >= 0)
            {
                ifRead = ((GlobalBoundaryCondition::IsSolidWallBC(iPart)) || ((GlobalBoundaryCondition::GetBC(iPart)->GetBCType() == PHENGLEI::ABLATION_SURFACE)));
            }
            if (!ifRead) continue;

            RDouble dCA_f, dCA_p, dCN_f, dCN_p, dCZ_f, dCZ_p, dcpx, dcpy, dcpz, dCl_f, dCl_p, dCn_f, dCn_p, dCm_f, dCm_p, hingeMoment;

            data_in->Read(&dCA_f, sizeof(RDouble));
            data_in->Read(&dCA_p, sizeof(RDouble));
            data_in->Read(&dCN_f, sizeof(RDouble));
            data_in->Read(&dCN_p, sizeof(RDouble));
            data_in->Read(&dCZ_f, sizeof(RDouble));
            data_in->Read(&dCZ_p, sizeof(RDouble));

            data_in->Read(&dcpx, sizeof(RDouble));
            data_in->Read(&dcpy, sizeof(RDouble));
            data_in->Read(&dcpz, sizeof(RDouble));

            data_in->Read(&dCl_f, sizeof(RDouble));
            data_in->Read(&dCl_p, sizeof(RDouble));
            data_in->Read(&dCn_f, sizeof(RDouble));
            data_in->Read(&dCn_p, sizeof(RDouble));
            data_in->Read(&dCm_f, sizeof(RDouble));
            data_in->Read(&dCm_p, sizeof(RDouble));

            data_in->Read(&hingeMoment, sizeof(RDouble));

            localForceVar[varCount ++] += static_cast<double>(dCA_f);
            localForceVar[varCount ++] += static_cast<double>(dCA_p);
            localForceVar[varCount ++] += static_cast<double>(dCN_f);
            localForceVar[varCount ++] += static_cast<double>(dCN_p);
            localForceVar[varCount ++] += static_cast<double>(dCZ_f);
            localForceVar[varCount ++] += static_cast<double>(dCZ_p);

            localForceVar[varCount ++] += static_cast<double>(dcpx);
            localForceVar[varCount ++] += static_cast<double>(dcpy);
            localForceVar[varCount ++] += static_cast<double>(dcpz);

            localForceVar[varCount ++] += static_cast<double>(dCl_f);
            localForceVar[varCount ++] += static_cast<double>(dCl_p);
            localForceVar[varCount ++] += static_cast<double>(dCn_f);
            localForceVar[varCount ++] += static_cast<double>(dCn_p);
            localForceVar[varCount ++] += static_cast<double>(dCm_f);
            localForceVar[varCount ++] += static_cast<double>(dCm_p);

            localForceVar[varCount ++] += static_cast<double>(hingeMoment);

            //! Because of global and local part force, when (iPart > ALL_SOLID_WALL), read second.
            if (iPart > ALL_SOLID_WALL)
            {
                data_in->Read(&dCA_f, sizeof(RDouble));
                data_in->Read(&dCA_p, sizeof(RDouble));
                data_in->Read(&dCN_f, sizeof(RDouble));
                data_in->Read(&dCN_p, sizeof(RDouble));
                data_in->Read(&dCZ_f, sizeof(RDouble));
                data_in->Read(&dCZ_p, sizeof(RDouble));

                data_in->Read(&dcpx, sizeof(RDouble));
                data_in->Read(&dcpy, sizeof(RDouble));
                data_in->Read(&dcpz, sizeof(RDouble));

                data_in->Read(&dCl_f, sizeof(RDouble));
                data_in->Read(&dCl_p, sizeof(RDouble));
                data_in->Read(&dCn_f, sizeof(RDouble));
                data_in->Read(&dCn_p, sizeof(RDouble));
                data_in->Read(&dCm_f, sizeof(RDouble));
                data_in->Read(&dCm_p, sizeof(RDouble));

                data_in->Read(&hingeMoment, sizeof(RDouble));

                localForceVar[varCount ++] += static_cast<double>(dCA_f);
                localForceVar[varCount ++] += static_cast<double>(dCA_p);
                localForceVar[varCount ++] += static_cast<double>(dCN_f);
                localForceVar[varCount ++] += static_cast<double>(dCN_p);
                localForceVar[varCount ++] += static_cast<double>(dCZ_f);
                localForceVar[varCount ++] += static_cast<double>(dCZ_p);

                localForceVar[varCount ++] += static_cast<double>(dcpx);
                localForceVar[varCount ++] += static_cast<double>(dcpy);
                localForceVar[varCount ++] += static_cast<double>(dcpz);

                localForceVar[varCount ++] += static_cast<double>(dCl_f);
                localForceVar[varCount ++] += static_cast<double>(dCl_p);
                localForceVar[varCount ++] += static_cast<double>(dCn_f);
                localForceVar[varCount ++] += static_cast<double>(dCn_p);
                localForceVar[varCount ++] += static_cast<double>(dCm_f);
                localForceVar[varCount ++] += static_cast<double>(dCm_p);

                localForceVar[varCount ++] += static_cast<double>(hingeMoment);
            }
        }
    }

    //! Collection.
    int myid = PHMPI::GetCurrentProcessorID();
    int server = PHMPI::GetServerProcessorID();

    double *globalForceVar = new double[nTotalVarDim];
    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_AllreduceSepMode(localForceVar, globalForceVar, nTotalVarDim, MPI_DOUBLE, PH_SUM);
    }
    else
    {
        for (int iVar = 0; iVar < nTotalVarDim; ++ iVar)
        {
            globalForceVar[iVar] = localForceVar[iVar];
        }
    }

    if (myid == server)
    {
        int varCount = 0;
        for (int iPart = 0; iPart < numberofTotalPart; ++ iPart)
        {
            CA_f[iPart] = static_cast<RDouble>(globalForceVar[varCount ++]);
            CA_p[iPart] = static_cast<RDouble>(globalForceVar[varCount ++]);
            CN_f[iPart] = static_cast<RDouble>(globalForceVar[varCount ++]);
            CN_p[iPart] = static_cast<RDouble>(globalForceVar[varCount ++]);
            CZ_f[iPart] = static_cast<RDouble>(globalForceVar[varCount ++]);
            CZ_p[iPart] = static_cast<RDouble>(globalForceVar[varCount ++]);

            cpx[iPart] = static_cast<RDouble>(globalForceVar[varCount ++]);
            cpy[iPart] = static_cast<RDouble>(globalForceVar[varCount ++]);
            cpz[iPart] = static_cast<RDouble>(globalForceVar[varCount ++]);

            Cl_f[iPart] = static_cast<RDouble>(globalForceVar[varCount ++]);
            Cl_p[iPart] = static_cast<RDouble>(globalForceVar[varCount ++]);
            Cn_f[iPart] = static_cast<RDouble>(globalForceVar[varCount ++]);
            Cn_p[iPart] = static_cast<RDouble>(globalForceVar[varCount ++]);
            Cm_f[iPart] = static_cast<RDouble>(globalForceVar[varCount ++]);
            Cm_p[iPart] = static_cast<RDouble>(globalForceVar[varCount ++]);

            hingeMoment[iPart] = static_cast<RDouble>(globalForceVar[varCount ++]);
        }
    }

    delete [] localForceVar;
    delete [] globalForceVar;
}

void Force::CompForce()
{
    RDouble sina, cosa, sinb, cosb;

    sina = sin(attack);
    cosa = cos(attack);
    sinb = sin(sideslip);
    cosb = cos(sideslip);

    bool isGlobalCoordinate = false;
    for (int iPart = 0; iPart < numberofTotalPart; ++ iPart)
    {
        RDouble scf, vcm, Sb, forceReferenceAreaPart, forceReferenceLengthPart, forceReferenceLengthPart_B;
        int bcID = partBCID[iPart];
        if ((bcID == ALL_SOLID_WALL) || (isGlobalCoordinate == true))
        {
            forceReferenceAreaPart = this->forceReferenceArea;
            forceReferenceLengthPart = this->forceReferenceLength;
            forceReferenceLengthPart_B = this->forceReferenceLengthSpanWise;
            isGlobalCoordinate = false;
        }
        else
        {
            SimpleBC *boundaryCondition = GlobalBoundaryCondition::GetBC(bcID);

            boundaryCondition->GetParamData("forceReferenceArea", &forceReferenceAreaPart, PHDOUBLE, 1);
            boundaryCondition->GetParamData("forceReferenceLength", &forceReferenceLengthPart, PHDOUBLE, 1);
            boundaryCondition->GetParamData("forceReferenceLengthSpanWise", &forceReferenceLengthPart_B, PHDOUBLE, 1);
            isGlobalCoordinate = true;
        }

        scf = 1.0 / forceReferenceAreaPart;
        vcm = 1.0 / (forceReferenceAreaPart * forceReferenceLengthPart);
        Sb = 1.0 / (forceReferenceAreaPart * forceReferenceLengthPart_B);

        CA_f[iPart] *= scf;
        CA_p[iPart] *= scf;
        CN_f[iPart] *= scf;
        CN_p[iPart] *= scf;
        CZ_f[iPart] *= scf;
        CZ_p[iPart] *= scf;

        cpx[iPart] *= scf;
        cpy[iPart] *= scf;
        cpz[iPart] *= scf;

        Cl_f[iPart] *= Sb;
        Cl_p[iPart] *= Sb;
        Cn_f[iPart] *= Sb;
        Cn_p[iPart] *= Sb;
        Cm_f[iPart] *= vcm;
        Cm_p[iPart] *= vcm;

        hingeMoment[iPart] *= vcm;
    }

    //* GRID_SIZE (F10.0) Number of grid nodes or cells
    //* MACH      (F10.2) Mach Number.
    //* ALPHA     (F10.3) Angle of Attack in Degrees.
    //* CL_TOT    (F10.3) Total Lift Coefficient.
    //* CD_TOT    (F10.5) Total Drag Coefficient.
    //* CD_PR     (F10.5) Surface-Pressure Integrated Drag Coefficient.
    //* CD_SF     (F10.5) Skin-Friction Integrated Drag Coefficient.
    //* CD-CL2/PA (F10.5) CD_TOT - (CL_TOT*CL_TOT)/(PI*AR), Note: AR=9.5
    //* CM_TOT    (F10.4) Total Moment Coefficient.

    if (!cl_tot || !cd_tot || !cd_pr || !cd_sf || !cd_cl2pa || !xcp || !side)
    {
        cl_tot = new RDouble[numberofTotalPart];
        cd_tot = new RDouble[numberofTotalPart];
        cd_pr = new RDouble[numberofTotalPart];
        cd_sf = new RDouble[numberofTotalPart];
        cd_cl2pa = new RDouble[numberofTotalPart];
        xcp = new RDouble[numberofTotalPart];
        ycp = new RDouble[numberofTotalPart];
        zcp = new RDouble[numberofTotalPart];
        side = new RDouble[numberofTotalPart];
    }

    const RDouble AR = 9.5;
    for (int iPart = 0; iPart < numberofTotalPart; ++ iPart)
    {
        RDouble cfx, cfy, cfz;
        cfx = CA_f[iPart] + CA_p[iPart];
        cfy = CN_f[iPart] + CN_p[iPart];
        cfz = CZ_f[iPart] + CZ_p[iPart];

        RDouble cmx, cmy, cmz;
        cmx = Cl_f[iPart] + Cl_p[iPart];
        cmy = Cn_f[iPart] + Cn_p[iPart];
        cmz = Cm_f[iPart] + Cm_p[iPart];

        cl_tot[iPart] = - sina * cfx + cosa * cfy;
        cd_tot[iPart] = cosa * cosb * cfx + sina * cosb * cfy + sinb * cfz;
        side[iPart] = -cosa * sinb * cfx - sina * sinb * cfy + cosb * cfz;
        cd_pr[iPart] = cosa * cosb * cpx[iPart] + sina * cosb * cpy[iPart] + sinb * cpz[iPart];
        cd_sf[iPart] = cd_tot[iPart] - cd_pr[iPart];
        cd_cl2pa[iPart] = cd_tot[iPart] - (cl_tot[iPart] * cl_tot[iPart]) / (PI * AR);

        xcp[iPart] = SIGN(one, cfy) * cmz / (ABS(cfy) + SMALL);
    }

    if (numberOfBody == 0)
    {
        return;
    }

    isGlobalCoordinate = false;
    for (int iPart = 0; iPart < numberofTotalPart; ++ iPart)
    {
        int bcID = partBCID[iPart];
        if (bcID == ALL_SOLID_WALL)
        {
            continue;
        }

        if (isGlobalCoordinate == false)
        {
            isGlobalCoordinate = true;
            continue;
        }

        int bodyIndex = static_cast<int>(GlobalBoundaryCondition::GetBodyIndex(bcID));

        CA_f_body[bodyIndex] += CA_f[iPart];
        CN_f_body[bodyIndex] += CN_f[iPart];
        CZ_f_body[bodyIndex] += CZ_f[iPart];
        Cl_f_body[bodyIndex] += Cl_f[iPart];
        Cn_f_body[bodyIndex] += Cn_f[iPart];
        Cm_f_body[bodyIndex] += Cm_f[iPart];
        CA_p_body[bodyIndex] += CA_p[iPart];
        CN_p_body[bodyIndex] += CN_p[iPart];
        CZ_p_body[bodyIndex] += CZ_p[iPart];
        Cl_p_body[bodyIndex] += Cl_p[iPart];
        Cn_p_body[bodyIndex] += Cn_p[iPart];
        Cm_p_body[bodyIndex] += Cm_p[iPart];

        isGlobalCoordinate = false;
    }
}

//! Bell modify 20120605.
void Force::DumpForce(int outnstep, std::ostringstream &oss)
{
    oss << setiosflags(ios::left);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    const int ALL = 0;
    oss << outnstep << "	";
    oss << cl_tot[ALL] << "	";
    oss << cd_tot[ALL] << "	";
    oss << side[ALL] << "	";
    oss << cd_pr[ALL] << "	";
    oss << cd_sf[ALL] << "	";
    oss << cd_cl2pa[ALL] << "	";
    oss << xcp[ALL] << "	";
    oss << CA_f[ALL] + CA_p[ALL] << "	";
    oss << CN_f[ALL] + CN_p[ALL] << "	";
    oss << CZ_f[ALL] + CZ_p[ALL] << "	";
    oss << Cl_f[ALL] + Cl_p[ALL] << "	";
    oss << Cn_f[ALL] + Cn_p[ALL] << "	";
    oss << Cm_f[ALL] + Cm_p[ALL] << "	";

    RDouble wallTime = PHSPACE::TIME_SPACE::GetWallTime();
    oss << wallTime;

    oss << "\n";
}

void Force::DumpPartForce(int outnstep, int partID, std::ostringstream &oss)
{
    oss << setiosflags(ios::left);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    oss << outnstep << "	";

    oss << cl_tot[partID] << "	";
    oss << cd_tot[partID] << "	";
    oss << side[partID] << "	";
    oss << cd_pr[partID] << "	";
    oss << cd_sf[partID] << "	";
    oss << cd_cl2pa[partID] << "	";
    oss << xcp[partID] << "	";
    oss << CA_f[partID] + CA_p[partID] << "	";
    oss << CN_f[partID] + CN_p[partID] << "	";
    oss << CZ_f[partID] + CZ_p[partID] << "	";
    oss << Cl_f[partID] + Cl_p[partID] << "	";
    oss << Cn_f[partID] + Cn_p[partID] << "	";
    oss << Cm_f[partID] + Cm_p[partID] << "	";

    oss << "\n";
}


void Force::DumpForce(ActionKey *actkey, int outnstep)
{
    CompForce();

    if (!actkey->IsNeedOpenFile())
    {
        return;
    }

    //! For total wall force.
    ParallelOpenFile(actkey);

    std::ostringstream oss;
    fstream &file = *(actkey->file);

    if (IfFileEmpty(file))
    {
        oss << "Title=\"THE AIRCOEF\"" << endl;
        oss << "Variables=" << endl;
        oss << "iter" << endl;
        oss << "<i>C<sub>L</sub></i>" << endl;
        oss << "<i>C<sub>D</sub></i>" << endl;
        oss << "<i>C<sub>Z</sub></i>" << endl;
        oss << "<i>C<sub>D_p</sub></i>" << endl;
        oss << "<i>C<sub>D_f</sub></i>" << endl;
        oss << "<i>C<sub>D</sub></i>-<i>C<sub>L</sub></i><sup>2</sup>/(<greek>p</greek><math>4</math>A<sub>r</sub>)</i>" << endl;
        oss << "<i>X<sub>cp</sub></i>" << endl;
        oss << "<i>C<sub>A</sub></i>" << endl;
        oss << "<i>C<sub>N</sub></i>" << endl;
        oss << "<i>C<sub>Z1</sub></i>" << endl;
        oss << "<i>C<sub>ml</sub></i>" << endl;
        oss << "<i>C<sub>mn</sub></i>" << endl;
        oss << "<i>C<sub>m</sub></i>" << endl;

        oss << "\"Walltime\"" << endl;
    }

    DumpForce(outnstep, oss);

    WriteASCIIFile(file, oss.str());

    ParallelCloseFile(actkey);

    //! For part force.
    bool isGlobalCoordinate = false;
    for (int iPart = 1; iPart < numberofTotalPart; ++ iPart)
    {
        int bcID = partBCID[iPart];

        SimpleBC *boundaryCondition = GlobalBoundaryCondition::GetBC(bcID);
        const string &bcName = boundaryCondition->GetBCName();

        string totalForceFileName = actkey->filename;
        string partAircoefFileName = PHSPACE::AddSymbolToFileName(totalForceFileName, "_", bcName);

        if (isGlobalCoordinate == false)
        {
            isGlobalCoordinate = true;
            continue;
        }
        else
        {
            isGlobalCoordinate = false;
        }

        fstream fileTmp;
        ParallelOpenFile(fileTmp, partAircoefFileName, ios_base::out | ios_base::app);

        ostringstream ossTmp;
        if (IfFileEmpty(fileTmp))
        {
            ossTmp << "Title=\"Aerodynamic Force of " << bcName << "\"" << endl;
            ossTmp << "Variables=" << endl;
            ossTmp << "iter" << endl;

            ossTmp << "<i>C<sub>L</sub></i>" << endl;
            ossTmp << "<i>C<sub>D</sub></i>" << endl;
            ossTmp << "<i>C<sub>Z</sub></i>" << endl;
            ossTmp << "<i>C<sub>D_p</sub></i>" << endl;
            ossTmp << "<i>C<sub>D_f</sub></i>" << endl;
            ossTmp << "<i>C<sub>D</sub></i>-<i>C<sub>L</sub></i><sup>2</sup>/(<greek>p</greek><math>4</math>A<sub>r</sub>)</i>" << endl;
            ossTmp << "<i>X<sub>cp</sub></i>" << endl;
            ossTmp << "<i>C<sub>A</sub></i>" << endl;
            ossTmp << "<i>C<sub>N</sub></i>" << endl;
            ossTmp << "<i>C<sub>Z1</sub></i>" << endl;
            ossTmp << "<i>C<sub>ml</sub></i>" << endl;
            ossTmp << "<i>C<sub>mn</sub></i>" << endl;
            ossTmp << "<i>C<sub>m</sub></i>" << endl;
        }
        DumpPartForce(outnstep, iPart, ossTmp);
        WriteASCIIFile(fileTmp, ossTmp.str());

        ParallelCloseFile(fileTmp);
    }

    DumpAllForceData(actkey, outnstep);
}

void Force::DumpAllForceData(ActionKey *actkey, int outnstep)
{
    string fileName = actkey->filename;

    string mainName, fileDir, extensionName;
    string globalFileName, localFileName;
    GetFileNameExtension(fileName, mainName, extensionName, ".");
    GetFileNameExtension(mainName, fileDir, extensionName, "/");

    ostringstream ossGlobalFileName, ossLocalFileName, ossGlobalData, ossLocalData;
    ossGlobalFileName << fileDir << "/force_global.dat";
    ossLocalFileName << fileDir << "/force_local.dat";
    globalFileName = ossGlobalFileName.str();
    localFileName = ossLocalFileName.str();

    string gridFile = GlobalDataBase::GetStrParaFromDB("gridfile");
    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    RDouble attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    RDouble angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");
    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    RDouble TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    RDouble TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    RDouble TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");
    RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");
    RDouble wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");
    RDouble height = 0.0;
    int inflowParaType = GlobalDataBase::GetIntParaFromDB("inflowParaType");
    if (inflowParaType == 1)
    {
        height = GlobalDataBase::GetDoubleParaFromDB("height");
    }

    //! Write global force data file title
    ossGlobalData << "#     MeshFile                     :       " << gridFile << endl;
    if (inflowParaType == 1)
    {
        ossGlobalData << "#     Height                       :       " << height << "  km" << endl;
    }

    ossGlobalData << setiosflags(ios::scientific);
    ossGlobalData << setprecision(6);
    ossGlobalData << "#     Reynolds                     :       " << refReNumber / reynoldsReferenceLengthDimensional << "  1/m" << endl;
    ossGlobalData << resetiosflags(ios::scientific);

    ossGlobalData << "#     Mach                         :       " << refMachNumber << endl;
    ossGlobalData << "#     Alpha                        :       " << attackd << "  degree" << endl;
    ossGlobalData << "#     Beta                         :       " << angleSlide << "  degree" << endl;

    if (wallTemperature < 0.0)
    {
        ossGlobalData << "#     Wall Temperature             :       " << "Adiabatic Wall" << endl;
    }
    else
    {
        ossGlobalData << "#     Wall Temperature             :       " << wallTemperature << "  K" << endl;
    }

    ossGlobalData << "#     gridScaleFactor              :       " << gridScaleFactor << endl;
    ossGlobalData << "#     forceReferenceLength         :       " << this->forceReferenceLength * reynoldsReferenceLengthDimensional << "  m" << endl;
    ossGlobalData << "#     forceReferenceLengthSpanWise :       " << this->forceReferenceLengthSpanWise * reynoldsReferenceLengthDimensional << "  m" << endl;

    if (GetDim() == THREE_D)
    {
        ossGlobalData << "#     forceReferenceArea           :       " << this->forceReferenceArea * reynoldsReferenceLengthDimensional * reynoldsReferenceLengthDimensional << "  m2" << endl;
    }
    else
    {
        ossGlobalData << "#     forceReferenceArea           :       " << this->forceReferenceArea * reynoldsReferenceLengthDimensional << "  m^2" << endl;
    }

    ossGlobalData << "#     TorqueRefX                   :       " << TorqueRefX * reynoldsReferenceLengthDimensional << endl;
    ossGlobalData << "#     TorqueRefY                   :       " << TorqueRefY * reynoldsReferenceLengthDimensional << endl;
    ossGlobalData << "#     TorqueRefZ                   :       " << TorqueRefZ * reynoldsReferenceLengthDimensional << endl;
    ossGlobalData << "#     iter                         :       " << outnstep << endl << endl;
    ossGlobalData << "#        CA             CN            CZ1           Cml             Cmn            Cm" << endl;
    ossGlobalData << "#        CA_p           CN_p          CZ1_p         Cml_p           Cmn_p          Cm_p" << endl;
    ossGlobalData << "#        CA_f           CN_f          CZ1_f         Cml_f           Cmn_f          Cm_f" << endl;

    //! Write local force data file title
    ossLocalData << "#     MeshFile             :       " << gridFile << endl;
    if (inflowParaType == 1)
    {
        ossLocalData << "#     Height               :       " << height << "  km" << endl;
    }

    ossLocalData << setiosflags(ios::scientific);
    ossLocalData << setprecision(6);
    ossLocalData << "#     Reynolds             :       " << refReNumber / reynoldsReferenceLengthDimensional << "  1/m" << endl;
    ossLocalData << resetiosflags(ios::scientific);

    ossLocalData << "#     Mach                 :       " << refMachNumber << endl;
    ossLocalData << "#     Alpha                :       " << attackd << "  degree" << endl;
    ossLocalData << "#     Beta                 :       " << angleSlide << "  degree" << endl;

    if (wallTemperature < 0.0)
    {
        ossLocalData << "#     Wall Temperature     :       " << "Adiabatic Wall" << endl;
    }
    else
    {
        ossLocalData << "#     Wall Temperature     :       " << wallTemperature << "  K" << endl;
    }


    ossLocalData << "#     gridScaleFactor      :       " << gridScaleFactor << endl;
    ossLocalData << "#     iter                 :       " << outnstep << endl << endl;
    ossLocalData << "#        CA             CN            CZ1           Cml             Cmn            Cm" << endl;
    ossLocalData << "#        CA_p           CN_p          CZ1_p         Cml_p           Cmn_p          Cm_p" << endl;
    ossLocalData << "#        CA_f           CN_f          CZ1_f         Cml_f           Cmn_f          Cm_f" << endl;
    ossLocalData << "#        CMh" << endl;

    bool isGlobalData = true;
    //! Write force data file
    for (int iPart = 0; iPart < numberofTotalPart; ++ iPart)
    {
        int bcID = partBCID[iPart];

        string partName;
        if (bcID == -1)
        {
            partName = "whole";
        }
        else
        {
            SimpleBC *boundaryCondition = GlobalBoundaryCondition::GetBC(bcID);
            const string &bcName = boundaryCondition->GetBCName();
            partName = bcName;
        }

        if (isGlobalData == true)
        {
            DumpAllGlobalForceData(partName, iPart, ossGlobalData);
            isGlobalData = false;
        }
        else
        {
            DumpAllLocalForceData(partName, iPart, ossLocalData);
            isGlobalData = true;
        }
    }

    for (int iBody = 0; iBody < numberOfBody; ++ iBody)
    {
        DumpBodyForceData(iBody, ossGlobalData);
    }

    //! Write global force data file
    fstream globalFile;
    ParallelOpenFile(globalFile, globalFileName, ios_base::out | ios_base::trunc);
    WriteASCIIFile(globalFile, ossGlobalData.str());
    ParallelCloseFile(globalFile);

    //! Write local force data file
    fstream localFile;
    ParallelOpenFile(localFile, localFileName, ios_base::out | ios_base::trunc);
    WriteASCIIFile(localFile, ossLocalData.str());
    ParallelCloseFile(localFile);
}

void Force::DumpAllLocalForceData(const string &partName, int partID, std::ostringstream &oss)
{
    RDouble TorqueRefX, TorqueRefY, TorqueRefZ, forceReferenceAreaPart, forceReferenceLengthPart, forceReferenceLengthSpanWisePart;
    int bcID = partBCID[partID];
    SimpleBC *boundaryCondition = GlobalBoundaryCondition::GetBC(bcID);

    boundaryCondition->GetParamData("forceReferenceArea", &forceReferenceAreaPart, PHDOUBLE, 1);
    boundaryCondition->GetParamData("forceReferenceLength", &forceReferenceLengthPart, PHDOUBLE, 1);
    boundaryCondition->GetParamData("forceReferenceLengthSpanWise", &forceReferenceLengthSpanWisePart, PHDOUBLE, 1);
    boundaryCondition->GetParamData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
    boundaryCondition->GetParamData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
    boundaryCondition->GetParamData("TorqueRefZ", &TorqueRefZ, PHDOUBLE, 1);

    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

    oss << resetiosflags(ios::scientific);

    oss << endl;
    oss << "*" << partName << endl;
    oss << "#forceReferenceLength         :    " << forceReferenceLengthPart * reynoldsReferenceLengthDimensional << "  m" << endl;
    oss << "#forceReferenceLengthSpanWise :    " << forceReferenceLengthSpanWisePart * reynoldsReferenceLengthDimensional << "  m" << endl;

    if (GetDim() == THREE_D)
    {
        oss << "#forceReferenceArea           :    " << forceReferenceAreaPart * reynoldsReferenceLengthDimensional * reynoldsReferenceLengthDimensional << "  m2" << endl;
    }
    else
    {
        oss << "#forceReferenceArea           :    " << forceReferenceAreaPart * reynoldsReferenceLengthDimensional << "  m^2" << endl;
    }

    oss << "#TorqueRefX                   :    " << TorqueRefX * reynoldsReferenceLengthDimensional << endl;
    oss << "#TorqueRefY                   :    " << TorqueRefY * reynoldsReferenceLengthDimensional << endl;
    oss << "#TorqueRefZ                   :    " << TorqueRefZ * reynoldsReferenceLengthDimensional << endl;

    oss << setiosflags(ios::right);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    oss << CA_f[partID] + CA_p[partID] << "	";
    oss << CN_f[partID] + CN_p[partID] << "	";
    oss << CZ_f[partID] + CZ_p[partID] << "	";
    oss << Cl_f[partID] + Cl_p[partID] << "	";
    oss << Cn_f[partID] + Cn_p[partID] << "	";
    oss << Cm_f[partID] + Cm_p[partID] << endl;

    oss << CA_p[partID] << "	";
    oss << CN_p[partID] << "	";
    oss << CZ_p[partID] << "	";
    oss << Cl_p[partID] << "	";
    oss << Cn_p[partID] << "	";
    oss << Cm_p[partID] << endl;

    oss << CA_f[partID] << "	";
    oss << CN_f[partID] << "	";
    oss << CZ_f[partID] << "	";
    oss << Cl_f[partID] << "	";
    oss << Cn_f[partID] << "	";
    oss << Cm_f[partID] << endl;

    oss << hingeMoment[partID] << endl;
}

void Force::DumpAllGlobalForceData(const string &partName, int partID, std::ostringstream &oss)
{
    oss << setiosflags(ios::right);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    oss << endl;
    oss << "*" << partName << endl;

    oss << CA_f[partID] + CA_p[partID] << "	";
    oss << CN_f[partID] + CN_p[partID] << "	";
    oss << CZ_f[partID] + CZ_p[partID] << "	";
    oss << Cl_f[partID] + Cl_p[partID] << "	";
    oss << Cn_f[partID] + Cn_p[partID] << "	";
    oss << Cm_f[partID] + Cm_p[partID] << endl;

    oss << CA_p[partID] << "	";
    oss << CN_p[partID] << "	";
    oss << CZ_p[partID] << "	";
    oss << Cl_p[partID] << "	";
    oss << Cn_p[partID] << "	";
    oss << Cm_p[partID] << endl;

    oss << CA_f[partID] << "	";
    oss << CN_f[partID] << "	";
    oss << CZ_f[partID] << "	";
    oss << Cl_f[partID] << "	";
    oss << Cn_f[partID] << "	";
    oss << Cm_f[partID] << endl;
}

void Force::DumpBodyForceData(int bodyID, std::ostringstream &oss)
{
    oss << setiosflags(ios::right);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    oss << endl;
    if (wallNameOfEachBody[bodyID].size() == 1)
    {
        return;
    }
    oss << "*" << wallNameOfEachBody[bodyID][0] << "(";
    oss << wallNameOfEachBody[bodyID][1];

    for (int ibc = 2; ibc < wallNameOfEachBody[bodyID].size(); ++ ibc)
    {
        oss << " + " << wallNameOfEachBody[bodyID][ibc];
    }
    oss << ")" << endl;

    oss << CA_f_body[bodyID] + CA_p_body[bodyID] << "	";
    oss << CN_f_body[bodyID] + CN_p_body[bodyID] << "	";
    oss << CZ_f_body[bodyID] + CZ_p_body[bodyID] << "	";
    oss << Cl_f_body[bodyID] + Cl_p_body[bodyID] << "	";
    oss << Cn_f_body[bodyID] + Cn_p_body[bodyID] << "	";
    oss << Cm_f_body[bodyID] + Cm_p_body[bodyID] << endl;

    oss << CA_p_body[bodyID] << "	";
    oss << CN_p_body[bodyID] << "	";
    oss << CZ_p_body[bodyID] << "	";
    oss << Cl_p_body[bodyID] << "	";
    oss << Cn_p_body[bodyID] << "	";
    oss << Cm_p_body[bodyID] << endl;

    oss << CA_f_body[bodyID] << "	";
    oss << CN_f_body[bodyID] << "	";
    oss << CZ_f_body[bodyID] << "	";
    oss << Cl_f_body[bodyID] << "	";
    oss << Cn_f_body[bodyID] << "	";
    oss << Cm_f_body[bodyID] << endl;
}

void Force::DumpData(DataContainer *data)
{
    data->MoveToBegin();

    PHWrite(data, numberofSolidWallPart);
    PHWrite(data, numberofTotalPart);

    PHWrite(data, partBCID, numberofTotalPart);

    PHWrite(data, cpx, numberofTotalPart);
    PHWrite(data, cpy, numberofTotalPart);
    PHWrite(data, cpz, numberofTotalPart);

    PHWrite(data, cl_tot, numberofTotalPart);
    PHWrite(data, cd_tot, numberofTotalPart);
    PHWrite(data, cd_pr, numberofTotalPart);
    PHWrite(data, cd_sf, numberofTotalPart);
    PHWrite(data, cd_cl2pa, numberofTotalPart);
    PHWrite(data, xcp, numberofTotalPart);
    PHWrite(data, side, numberofTotalPart);
    PHWrite(data, hingeMoment, numberofTotalPart);

}

bool AerodynamicForceComponent::initialize = false;
vector< std::set< string > >  AerodynamicForceComponent::components;

AerodynamicForceComponent::AerodynamicForceComponent()
{
}

AerodynamicForceComponent::~AerodynamicForceComponent()
{
}

void AerodynamicForceComponent::Initialize()
{
    if (!AerodynamicForceComponent::initialize)
    {
        int numberOfAerodynamicForceComponents = GlobalDataBase::GetIntParaFromDB("numberOfAerodynamicForceComponents");
        AerodynamicForceComponent::components.resize(numberOfAerodynamicForceComponents);
        int nTBCWall = GlobalBoundaryCondition::GetNumberOfSolidWallPart();
        SimpleBC *boundaryCondition = 0;
        int iComponent = 0;
        int nTBC = static_cast<int>(GlobalBoundaryCondition::GetNumberOfBC());
        if (numberOfAerodynamicForceComponents != nTBCWall)
        {
            TK_Exit::ExceptionExit("Error in AerodynamicForceComponent::initialize : numberOfAerodynamicForceComponents != nTBCWall");

        }
        for (int iBC = 0; iBC < nTBC; ++ iBC)
        {
            if (GlobalBoundaryCondition::IsSolidWallBC(iBC))
            {
                boundaryCondition = GlobalBoundaryCondition::GetBC(iBC);
                string bcName = boundaryCondition->GetBCName();
                AerodynamicForceComponent::components[iComponent++].insert(bcName);
            }
        }

        AerodynamicForceComponent::initialize = true;
    }
}

bool AerodynamicForceComponent::FaceInAerodynamicForceComponents(string componentName, int aerodynamicForceComponentIndex)
{
    AerodynamicForceComponent::Initialize();
    set< string >::iterator iter = AerodynamicForceComponent::components[aerodynamicForceComponentIndex].find(componentName);
    return iter != AerodynamicForceComponent::components[aerodynamicForceComponentIndex].end();
}

bool FaceInAerodynamicForceComponents(string componentName, int aerodynamicForceComponentIndex)
{
    return AerodynamicForceComponent::FaceInAerodynamicForceComponents(componentName, aerodynamicForceComponentIndex);
}

}

