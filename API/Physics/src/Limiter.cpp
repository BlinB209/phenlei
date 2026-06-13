#include "Limiter.h"
#include "Glb_Dimension.h"
#include "Geo_StructBC.h"
#include "Geo_UnstructBC.h"
#include "TK_Exit.h"
#include "Math_Limiter.h"
using namespace std;

namespace PHSPACE
{

LIB_EXPORT Limiter::Limiter(UnstructGrid * gridin, RDouble ** q, int nVariable, const string &limiterVarName)
{
    this->grid              = gridin;
    this->q                 = q;
    this->nLimiterVariables = nVariable;
    this->limiterVarName    = limiterVarName;

    limit = new Data2D<RDouble>();

    limitModel     = 0;
    limitVector   = 0;
    limiterTypeID  = 1;

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal     = nTotalCell + nBoundFace;

    this->CreateLimiter(nTotal, nLimiterVariables);
}

LIB_EXPORT Limiter::~Limiter()
{
    delete limit;    limit = NULL;
}

void Limiter::CreateLimiter(int nsize, int neqn)
{
    limit->CreateData(nsize, neqn);
}

LIB_EXPORT Limiter::iterator Limiter::GetLimiter(int m)
{
    return limit->GetData(m);
}

LIB_EXPORT void Limiter::Calculation(Grid *gridIn)
{
    if (limiterTypeID >= ILMT_STRUCT || limiterTypeID == ILMT_NOLIM || limiterTypeID == ILMT_FIRST)
    {
        return;
    }

    UnstructGrid *gridUnstruct = UnstructGridCast(gridIn);

    RDouble **dqdx = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarZ"));
    using namespace IDX;

    const int SCALAR_LIMITER    = 0;
    const int VECTOR_LIMITER    = 1;

    const int LIMIT_BY_RHO_P    = 0;
    const int LIMIT_BY_EACH_VAR = 1;

    RDouble *limit = 0;
    if (limiterTypeID == ILMT_BARTH)
    {
        if (limitVector == SCALAR_LIMITER)
        {
            //! Each variables has the same limiter value.
            limit = this->GetLimiter(0);
            if (limitModel == LIMIT_BY_RHO_P)
            {
                BarthLimiter(limit, q[IR], dqdx[IR], dqdy[IR], dqdz[IR]);
                BarthLimiter(limit, q[IP], dqdx[IP], dqdy[IP], dqdz[IP]);
            }
            else if (limitModel == LIMIT_BY_EACH_VAR)
            {
                for (int iVariable = 0; iVariable < nLimiterVariables; ++ iVariable)
                {
                    BarthLimiter(limit, q[iVariable], dqdx[iVariable], dqdy[iVariable], dqdz[iVariable]);
                }
            }
        }
        else if (limitVector == VECTOR_LIMITER)
        {
            //! Limiter value of each variable is computed by itself.
            for (int iVariable = 0; iVariable < nLimiterVariables; ++ iVariable)
            {
                limit = this->GetLimiter(iVariable);
                BarthLimiter(limit, q[iVariable], dqdx[iVariable], dqdy[iVariable], dqdz[iVariable]);
            }
        }
    }
    else if (limiterTypeID == ILMT_VENCAT)
    {
        if (limitVector == SCALAR_LIMITER)
        {
            //! The limiter values of each variable are same.
            limit = this->GetLimiter(0);
            if (limitModel == LIMIT_BY_RHO_P)
            {
                VencatLimiter(limit, q[IR], dqdx[IR], dqdy[IR], dqdz[IR], IR);
                VencatLimiter(limit, q[IP], dqdx[IP], dqdy[IP], dqdz[IP], IP);
            }
            else if (limitModel == LIMIT_BY_EACH_VAR)
            {
                for (int m = 0; m < nLimiterVariables; ++ m)
                {
                    VencatLimiter(limit, q[m], dqdx[m], dqdy[m], dqdz[m], m);
                }
            }
        }
        else if (limitVector == VECTOR_LIMITER)
        {
            //! The limiter values of each variable are different.
            for (int m = 0; m < nLimiterVariables; ++ m)
            {
                limit = this->GetLimiter(m);
                VencatLimiter(limit, q[m], dqdx[m], dqdy[m], dqdz[m], m);
            }
        }
    }
    else
    {
        TK_Exit::ExceptionExit("Error: this unstructured limiter is not exist !\n", true);
    }

    //! Set limiter at boundary.
    if (limitVector == SCALAR_LIMITER)
    {
        limit = this->GetLimiter(0);
        SetLimitBoundary(gridUnstruct, limit);
    }
    else if (limitVector == VECTOR_LIMITER)
    {
        for (int m = 0; m < nLimiterVariables; ++ m)
        {
            limit = this->GetLimiter(m);
            SetLimitBoundary(gridUnstruct, limit);
        }
    }
}

LIB_EXPORT void Limiter::InitLimitData()
{
    limit->SetData(one);

    if (limiterTypeID == ILMT_FIRST)
    {
        limit->SetData(zero);
    }
}

void Limiter::BarthLimiter(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *xfc  = grid->GetFaceCenterX();
    RDouble *yfc  = grid->GetFaceCenterY();
    RDouble *zfc  = grid->GetFaceCenterZ();
    RDouble *xcc  = grid->GetCellCenterX();
    RDouble *ycc  = grid->GetCellCenterY();
    RDouble *zcc  = grid->GetCellCenterZ();
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();

    //! Need temporary arrays for the differences between the value in every cell and
    //! maximum/minimum in the neighboring cells.
    RDouble *dmin = new RDouble[nTotal];
    RDouble *dmax = new RDouble[nTotal];

    //! Find the the differences for q.
    MinMaxDiff(dmin, dmax, q);

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le     = leftCellOfFace [iFace];

        RDouble dx     = xfc[iFace] - xcc[le];
        RDouble dy     = yfc[iFace] - ycc[le];
        RDouble dz     = zfc[iFace] - zcc[le];
        RDouble dqFace = dqdx[le] * dx + dqdy[le] * dy + dqdz[le] * dz;

        RDouble ds     = sqrt(dx * dx + dy * dy + dz * dz);
        RDouble dot    = (xfn[iFace] * dx + yfn[iFace] * dy + zfn[iFace] * dz) / (ds + SMALL);
        
        if (dot < 0.0)
        {
            limit[le] = 0.0;
        }
        else if (dqFace > dmax[le])
        {
            limit[le] = MIN(limit[le], dmax[le] / (dqFace + SMALL));
        }
        else if (dqFace < dmin[le])
        {
            limit[le] = MIN(limit[le], dmin[le] / (dqFace + SMALL));
        }
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le     = leftCellOfFace [iFace];
        int re     = rightCellOfFace[iFace];

        RDouble dx     = xfc[iFace] - xcc[le];
        RDouble dy     = yfc[iFace] - ycc[le];
        RDouble dz     = zfc[iFace] - zcc[le];
        RDouble dqFace = dqdx[le] * dx + dqdy[le] * dy + dqdz[le] * dz;

        RDouble ds     = sqrt(dx * dx + dy * dy + dz * dz);
        RDouble dot    = (xfn[iFace] * dx + yfn[iFace] * dy + zfn[iFace] * dz) / (ds + SMALL);

        if (dot < 0.0)
        {
            limit[le] = 0.0;
        }
        else if (dqFace > dmax[le])
        {
            limit[le] = MIN(limit[le], dmax[le] / (dqFace + SMALL));
        }
        else if (dqFace < dmin[le])
        {
            limit[le] = MIN(limit[le], dmin[le] / (dqFace + SMALL));
        }

        dx     = xfc[iFace] - xcc[re];
        dy     = yfc[iFace] - ycc[re];
        dz     = zfc[iFace] - zcc[re];

        dqFace = dqdx[re] * dx + dqdy[re] * dy + dqdz[re] * dz;

        ds     = sqrt(dx * dx + dy * dy + dz * dz);
        dot    = - (xfn[iFace] * dx + yfn[iFace] * dy + zfn[iFace] * dz) / (ds + SMALL);

        if (dot < 0.0)
        {
            limit[re] = 0.0;
        }
        else if (dqFace > dmax[re])
        {
            limit[re] = MIN(limit[re], dmax[re] / (dqFace + SMALL));
        }
        else if (dqFace < dmin[re])
        {
            limit[re] = MIN(limit[re], dmin[re] / (dqFace + SMALL));
        }
    }

    delete [] dmin;    dmin = NULL;
    delete [] dmax;    dmax = NULL;
}

/************************************************************************/
/*         Control switch of Vencat limiter, Bell 20111226              */
/************************************************************************/
void Limiter::VencatLimiter(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, int iVariable)
{
    switch (vencatLimiterType)
    {
    case 0:
        VencatLimiter0(limit, q, dqdx, dqdy, dqdz);
        break;
    case 1:
        VencatLimiter1(limit, q, dqdx, dqdy, dqdz);
        break;
    case 2:
        VencatLimiter2(limit, q, dqdx, dqdy, dqdz, iVariable);
        break;
    case 3:
        VencatLimiter3(limit, q, dqdx, dqdy, dqdz, iVariable);
        break;
    case 4:
        VencatLimiter4(limit, q, dqdx, dqdy, dqdz);
        break;
    case 5:
        VencatLimiter5(limit, q, dqdx, dqdy, dqdz);
        break;
    case 6:
        VencatLimiter6(limit, q, dqdx, dqdy, dqdz);
        break;
    case 7:
        VencatLimiter7(limit, q, dqdx, dqdy, dqdz, iVariable);
        break;
    default:
        TK_Exit::UnexpectedVarValue("vencat limiter type", vencatLimiterType);
        break;
    }
}
/************************************************************************
Dr. Wang improved vencat limiter, independent of grid unit.
This limiter is used in Fluent.
Reference: AIAA-96-2091
A fast nested multigrid viscous flow solver for adaptive Cartesian Quad grids
venkatCoeff = [0.01, 0.2]
But, the dmin and dmax should be the whole grid rather than one partition grid.
************************************************************************/
void Limiter::VencatLimiter0(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal     = nTotalCell + nBoundFace;
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *xfc   = grid->GetFaceCenterX();
    RDouble *yfc   = grid->GetFaceCenterY();
    RDouble *zfc   = grid->GetFaceCenterZ();
    RDouble *xcc   = grid->GetCellCenterX();
    RDouble *ycc   = grid->GetCellCenterY();
    RDouble *zcc   = grid->GetCellCenterZ();

    RDouble dx, dy, dz;
    RDouble dqFace, tmp;

    RDouble *dmin = new RDouble[nTotal];
    RDouble *dmax = new RDouble[nTotal];

    //! Find the the differences for q.
    MinMaxDiff(dmin, dmax, q);

    //! Find the maximum and minimum of dq in the field to set epsilons.
    RDouble dqMin = dmin[0], dqMax = dmax[0];

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dqMin = MIN(dqMin, dmin[iCell]);
        dqMax = MAX(dqMax, dmax[iCell]);
    }

    grid->CompareMaxMinValue(dqMax, 1);
    grid->CompareMaxMinValue(dqMin, 2);

    RDouble venkatCoeff = 1.0e-5;
    GlobalDataBase::GetData("venkatCoeff", &venkatCoeff, PHDOUBLE, 1);

    RDouble eps;
    eps = venkatCoeff * (dqMax - dqMin);
    eps = eps * eps + SMALL;

    int le, re;

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        dx      = xfc[iFace] - xcc[le];
        dy      = yfc[iFace] - ycc[le];
        dz      = zfc[iFace] - zcc[le];
        dqFace  = dqdx[le] * dx + dqdy[le] * dy + dqdz[le] * dz;

        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        dx     = xfc[iFace] - xcc[le];
        dy     = yfc[iFace] - ycc[le];
        dz     = zfc[iFace] - zcc[le];
        dqFace = dqdx[le] * dx + dqdy[le] * dy + dqdz[le] * dz;

        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }

        dx     = xfc[iFace] - xcc[re];
        dy     = yfc[iFace] - ycc[re];
        dz     = zfc[iFace] - zcc[re];
        dqFace = dqdx[re] * dx + dqdy[re] * dy + dqdz[re] * dz;

        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[re], dqFace, eps);
            limit[re] = MIN(limit[re], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[re], dqFace, eps);
            limit[re] = MIN(limit[re], tmp);
        }
    }

    delete [] dmin;    dmin = NULL;
    delete [] dmax;    dmax = NULL;
}

void Limiter::VencatLimiter1(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *xfc   = grid->GetFaceCenterX();
    RDouble *yfc   = grid->GetFaceCenterY();
    RDouble *zfc   = grid->GetFaceCenterZ();
    RDouble *xcc   = grid->GetCellCenterX();
    RDouble *ycc   = grid->GetCellCenterY();
    RDouble *zcc   = grid->GetCellCenterZ();

    RDouble dx, dy, dz, ds;
    RDouble dqFace, tmp;

    RDouble *dmin = new RDouble[nTotal];
    RDouble *dmax = new RDouble[nTotal];

    //! Find the the differences for q.
    MinMaxDiff(dmin, dmax, q);

    //! Find the maximum and minimum of dq in the field to set epsilons.
    RDouble dqMin = dmin[0], dqMax = dmax[0];

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dqMin = MIN(dqMin, dmin[iCell]);
        dqMax = MAX(dqMax, dmax[iCell]);
    }

    RDouble venkatCoeff = 1.0e-5;
    GlobalDataBase::GetData("venkatCoeff", &venkatCoeff, PHDOUBLE, 1);

    RDouble eps, eps_ref;
    eps_ref = venkatCoeff * (dqMax - dqMin);
    eps = eps_ref + SMALL;
    //eps = eps * eps + SMALL;

    int le, re;

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        dx = xcc[re] - xcc[le];
        dy = ycc[re] - ycc[le];
        dz = zcc[re] - zcc[le];
        
        ds  = sqrt(dx * dx + dy * dy + dz * dz);
        eps = ds;
        eps = eps_ref * eps * eps * eps;

        dx      = xfc[iFace] - xcc[le];
        dy      = yfc[iFace] - ycc[le];
        dz      = zfc[iFace] - zcc[le];
        dqFace  = dqdx[le] * dx + dqdy[le] * dy + dqdz[le] * dz;

        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        dx = xcc[re] - xcc[le];
        dy = ycc[re] - ycc[le];
        dz = zcc[re] - zcc[le];
        
        ds  = sqrt(dx * dx + dy * dy + dz * dz);
        eps = ds;
        eps = eps_ref * eps * eps * eps;

        dx     = xfc[iFace] - xcc[le];
        dy     = yfc[iFace] - ycc[le];
        dz     = zfc[iFace] - zcc[le];
        dqFace = dqdx[le] * dx + dqdy[le] * dy + dqdz[le] * dz;

        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }

        dx     = xfc[iFace] - xcc[re];
        dy     = yfc[iFace] - ycc[re];
        dz     = zfc[iFace] - zcc[re];
        dqFace = dqdx[re] * dx + dqdy[re] * dy + dqdz[re] * dz;

        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[re], dqFace, eps);
            limit[re] = MIN(limit[re], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[re], dqFace, eps);
            limit[re] = MIN(limit[re], tmp);
        }
    }

    delete [] dmin;    dmin = NULL;
    delete [] dmax;    dmax = NULL;
}

//! This limiter type was proposed by CXH.
//! Which has nothing to do with mesh size, mach number.
//! But it need to be V&V in more numerical experiments.
void Limiter::VencatLimiter2(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, int iVariable)
{
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;      
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *xfc   = grid->GetFaceCenterX();
    RDouble *yfc   = grid->GetFaceCenterY();
    RDouble *zfc   = grid->GetFaceCenterZ();
    RDouble *xcc   = grid->GetCellCenterX();
    RDouble *ycc   = grid->GetCellCenterY();
    RDouble *zcc   = grid->GetCellCenterZ();

    RDouble dx, dy, dz;
    RDouble dqFace, tmp;

    RDouble *dmin = new RDouble[nTotal];
    RDouble *dmax = new RDouble[nTotal];

    //! Find the the differences for q.
    MinMaxDiff(dmin, dmax, q);

    RDouble venkatCoeff = GlobalDataBase::GetDoubleParaFromDB("venkatCoeff");
    RDouble globalTotalCells = GlobalDataBase::GetDoubleParaFromDB("UnstrGlobalTotalCells");

    //! The following 'eps' is the eps^2 actually in the mathematic formation.
    //! And the 'venkatCoeff' is the coefficient of 'K' in the mathematic formation.
    int dimension = GetDim();
    RDouble n = 3.0;
    RDouble eps;
    eps = pow((venkatCoeff * pow((1.0 / globalTotalCells), (1.0 / dimension))), n);
    if (iVariable == IDX::IP)
    {
        RDouble * primInflow = reinterpret_cast< RDouble * >(GlobalDataBase::GetDataPtr("prim_inf"));
        eps *= primInflow[IDX::IP] * primInflow[IDX::IP];
    }

    int le, re;

    //! The following loop is to calculate the limiter coefficient of each cell on the boundary face.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        //! Get the left and right cell index of iFace.
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        //! Compute the vector from face center to left cell center.
        dx = xfc[iFace] - xcc[le];
        dy = yfc[iFace] - ycc[le];
        dz = zfc[iFace] - zcc[le];

        //! Compute the dq from left cell center to face center.
        dqFace  = dqdx[le] * dx + dqdy[le] * dy + dqdz[le] * dz;
        if (dqFace > SMALL)
        {
            //! According to the vencat limiter of Theory manual.
            tmp       = VenFun(dmax[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }
        else if (dqFace < - SMALL)
        {
            //! According to the vencat limiter of Theory manual.
            tmp       = VenFun(dmin[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }
    }

    //! The following loop is to calculate the limiter coefficient of each interior cell.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        dx     = xfc[iFace] - xcc[le];
        dy     = yfc[iFace] - ycc[le];
        dz     = zfc[iFace] - zcc[le];
        dqFace = dqdx[le] * dx + dqdy[le] * dy + dqdz[le] * dz;
        if (dqFace > SMALL)
        {
            tmp      = VenFun(dmax[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp      = VenFun(dmin[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }

        dx     = xfc[iFace] - xcc[re];
        dy     = yfc[iFace] - ycc[re];
        dz     = zfc[iFace] - zcc[re];
        dqFace = dqdx[re] * dx + dqdy[re] * dy + dqdz[re] * dz;
        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[re], dqFace, eps);
            limit[re] = MIN(limit[re], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[re], dqFace, eps);
            limit[re] = MIN(limit[re], tmp);
        }
    }

    delete [] dmin;    dmin = NULL;
    delete [] dmax;    dmax = NULL;
}

//! This limiter type was modified by NeverMore.
//! Which the limiter coefficient of eps^2 is influenced by cell volume,
//! other than distance between cell center and face center, that was used in VencatLimiter4.
void Limiter::VencatLimiter3(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, int iVariable)
{
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *xfc    = grid->GetFaceCenterX();
    RDouble *yfc    = grid->GetFaceCenterY();
    RDouble *zfc    = grid->GetFaceCenterZ();
    RDouble *xcc    = grid->GetCellCenterX();
    RDouble *ycc    = grid->GetCellCenterY();
    RDouble *zcc    = grid->GetCellCenterZ();
    RDouble *volume = grid->GetCellVolume();

    RDouble dx, dy, dz;
    RDouble dqFace, tmp;
    int le, re;

    RDouble *dmin = new RDouble[nTotal];
    RDouble *dmax = new RDouble[nTotal];

    //! Find the the differences for q.
    MinMaxDiff(dmin, dmax, q);

    RDouble globalTotalVolume;
    GlobalDataBase::GetData("UnstrGlobalTotalVolume", &globalTotalVolume, PHDOUBLE, 1);

    //! The following 'eps' is the eps^2 actually in the mathematic formation.
    //! And the 'venkatCoeff' is the coefficient of 'K' in the mathematic formation.
    RDouble venkatCoeff = 1.0e-5;
    GlobalDataBase::GetData("venkatCoeff", &venkatCoeff, PHDOUBLE, 1);
    RDouble *primInflow = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("prim_inf"));
        
    RDouble n = 3.0;
    int dimension = GetDim();
    RDouble dim1 = 1.0 / (dimension * 1.0);

    //! Calculate eps for each cell.
    RDouble squareOfReferenceVariable = primInflow[IDX::IP] * primInflow[IDX::IP];
    RDouble *eps = new RDouble [nTotalCell];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        eps[iCell] = pow(venkatCoeff * pow(volume[iCell] / globalTotalVolume, dim1), n);
        if (iVariable == IDX::IP)
        {
            eps[iCell]  *= squareOfReferenceVariable;
        }
    }

    //! The following loop is to calculate the limiter coefficient of each cell on the boundary face.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        dx      = xfc[iFace] - xcc[le];
        dy      = yfc[iFace] - ycc[le];
        dz      = zfc[iFace] - zcc[le];

        dqFace  = dqdx[le] * dx + dqdy[le] * dy + dqdz[le] * dz;
        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[le], dqFace, eps[le]);
            limit[le] = MIN(limit[le], tmp);

        }
        else if (dqFace < - SMALL)
        { 
            tmp       = VenFun(dmin[le], dqFace, eps[le]);
            limit[le] = MIN(limit[le], tmp);
        }
    }

    //! The following loop is to calculate the limiter coefficient of each interior cell.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        dx     = xfc[iFace] - xcc[le];
        dy     = yfc[iFace] - ycc[le];
        dz     = zfc[iFace] - zcc[le];
        dqFace = dqdx[le] * dx + dqdy[le] * dy + dqdz[le] * dz;
        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[le], dqFace, eps[le]);
            limit[le] = MIN(limit[le], tmp);

        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[le], dqFace, eps[le]);
            limit[le] = MIN(limit[le], tmp);
        }

        dx     = xfc[iFace] - xcc[re];
        dy     = yfc[iFace] - ycc[re];
        dz     = zfc[iFace] - zcc[re];
        dqFace = dqdx[re] * dx + dqdy[re] * dy + dqdz[re] * dz;
        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[re], dqFace, eps[re]);
            limit[re] = MIN(limit[re], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[re], dqFace, eps[re]);
            limit[re] = MIN(limit[re], tmp);
        }
    }

    delete [] dmin;    dmin = NULL;
    delete [] dmax;    dmax = NULL;
    delete [] eps;     eps = NULL;
}

/************************************************************************/
/*  According to the literature formula,related to the grid unit        */
/*                        Bell 20121028 add                             */
/************************************************************************/
void Limiter::VencatLimiter4(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *xfc   = grid->GetFaceCenterX();
    RDouble *yfc   = grid->GetFaceCenterY();
    RDouble *zfc   = grid->GetFaceCenterZ();
    RDouble *xcc   = grid->GetCellCenterX();
    RDouble *ycc   = grid->GetCellCenterY();
    RDouble *zcc   = grid->GetCellCenterZ();

    RDouble dx, dy, dz, ds;
    RDouble dqFace, tmp;

    RDouble *dmin = new RDouble[nTotal];
    RDouble *dmax = new RDouble[nTotal];

    //! Find the the differences for q.
    MinMaxDiff(dmin, dmax, q);

    //! Find the maximum and minimum of dq in the field to set epsilons.
    RDouble dqMin = dmin[0], dqMax = dmax[0];

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dqMin = MIN(dqMin, dmin[iCell]);
        dqMax = MAX(dqMax, dmax[iCell]);
    }

    //! Bell 20131123 mod.
    //! Bell 20130307 add.
    grid->CompareMaxMinValue(dqMax, 1);
    grid->CompareMaxMinValue(dqMin, 2);
    //*/

    RDouble venkatCoeff = 1.0e-5;
    GlobalDataBase::GetData("venkatCoeff", &venkatCoeff, PHDOUBLE, 1);

    RDouble eps;

    int le, re;

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        dx = xcc[re] - xcc[le];
        dy = ycc[re] - ycc[le];
        dz = zcc[re] - zcc[le];

        ds  = sqrt(dx * dx + dy * dy + dz * dz);
        eps = venkatCoeff * venkatCoeff * venkatCoeff * ds * ds * ds;

        dx      = xfc[iFace] - xcc[le];
        dy      = yfc[iFace] - ycc[le];
        dz      = zfc[iFace] - zcc[le];
        dqFace  = dqdx[le] * dx + dqdy[le] * dy + dqdz[le] * dz;

        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        dx = xcc[re] - xcc[le];
        dy = ycc[re] - ycc[le];
        dz = zcc[re] - zcc[le];

        ds  = sqrt(dx * dx + dy * dy + dz * dz);
        eps = venkatCoeff * venkatCoeff * venkatCoeff * ds * ds * ds;

        dx     = xfc[iFace] - xcc[le];
        dy     = yfc[iFace] - ycc[le];
        dz     = zfc[iFace] - zcc[le];
        dqFace = dqdx[le] * dx + dqdy[le] * dy + dqdz[le] * dz;

        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }

        dx     = xfc[iFace] - xcc[re];
        dy     = yfc[iFace] - ycc[re];
        dz     = zfc[iFace] - zcc[re];
        dqFace = dqdx[re] * dx + dqdy[re] * dy + dqdz[re] * dz;

        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[re], dqFace, eps);
            limit[re] = MIN(limit[re], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[re], dqFace, eps);
            limit[re] = MIN(limit[re], tmp);
        }
    }

    delete [] dmin;    dmin = NULL;
    delete [] dmax;    dmax = NULL;
}

/************************************************************************/
/*       According to the grid volume,related to the grid unit          */
/*                          Bell 20161114 add                           */
/************************************************************************/
void Limiter::VencatLimiter5(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *xfc    = grid->GetFaceCenterX();
    RDouble *yfc    = grid->GetFaceCenterY();
    RDouble *zfc    = grid->GetFaceCenterZ();
    RDouble *xcc    = grid->GetCellCenterX();
    RDouble *ycc    = grid->GetCellCenterY();
    RDouble *zcc    = grid->GetCellCenterZ();
    RDouble *volume = grid->GetCellVolume();

    RDouble dx, dy, dz;
    RDouble dqFace, tmp;

    RDouble *dmin = new RDouble[nTotal];
    RDouble *dmax = new RDouble[nTotal];

    //! Find the the differences for q.
    MinMaxDiff(dmin, dmax, q);

    RDouble venkatCoeff = 1.0e-5;
    GlobalDataBase::GetData("venkatCoeff", &venkatCoeff, PHDOUBLE, 1);

    RDouble eps;

    int le, re;

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];

        eps = venkatCoeff * venkatCoeff * venkatCoeff * volume[le];

        dx      = xfc[iFace] - xcc[le];
        dy      = yfc[iFace] - ycc[le];
        dz      = zfc[iFace] - zcc[le];
        dqFace  = dqdx[le] * dx + dqdy[le] * dy + dqdz[le] * dz;

        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        eps = venkatCoeff * venkatCoeff * venkatCoeff * volume[le];
        dx     = xfc[iFace] - xcc[le];
        dy     = yfc[iFace] - ycc[le];
        dz     = zfc[iFace] - zcc[le];
        dqFace = dqdx[le] * dx + dqdy[le] * dy + dqdz[le] * dz;

        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[le], dqFace, eps);
            limit[le] = MIN(limit[le], tmp);
        }

        eps = venkatCoeff * venkatCoeff * venkatCoeff * volume[re];
        dx     = xfc[iFace] - xcc[re];
        dy     = yfc[iFace] - ycc[re];
        dz     = zfc[iFace] - zcc[re];
        dqFace = dqdx[re] * dx + dqdy[re] * dy + dqdz[re] * dz;

        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[re], dqFace, eps);
            limit[re] = MIN(limit[re], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[re], dqFace, eps);
            limit[re] = MIN(limit[re], tmp);
        }
    }

    delete [] dmin;    dmin = NULL;
    delete [] dmax;    dmax = NULL;
}

//! Multidimensional limiting Process.
void Limiter::VencatLimiter6(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;
    int iNode,j;
    int nTotalNode = grid->GetNTotalNode();
    int *nCPN = grid->GetCellNumberOfEachNode();
    int *npc  = grid->GetNodeNumberOfEachCell();
    int *n2c  = grid->GetNode2Cell();
    int *c2n  = grid->GetCell2Node();
    RDouble *xn = grid->GetX();
    RDouble *yn = grid->GetY();
    RDouble *zn = grid->GetZ();

    RDouble *xcc    = grid->GetCellCenterX();
    RDouble *ycc    = grid->GetCellCenterY();
    RDouble *zcc    = grid->GetCellCenterZ();
    RDouble *volume = grid->GetCellVolume();

    RDouble globalTotalVolume = 0.0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        globalTotalVolume += volume[iCell];  
    }
    RDouble **qAll = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    using namespace IDX;

    RDouble dx, dy, dz, ds;
    RDouble dqNode, tmp;
    RDouble dqmax, dqmin;

    RDouble venkatCoeff = 1.0e-5;
    GlobalDataBase::GetData("venkatCoeff", &venkatCoeff, PHDOUBLE, 1);

    RDouble eps;

    RDouble *qmax = new RDouble[nTotalNode];
    RDouble *qmin = new RDouble[nTotalNode];

    int count = 0;
    for (iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        qmax[iNode] = -1;
        qmin[iNode] = 100000;
        for (j = 0; j < nCPN[ iNode ]; ++ j)
        {
            qmax[iNode] = max(qmax[iNode], q[n2c[count]]);
            qmin[iNode] = min(qmin[iNode], q[n2c[count]]);
            count ++;
        }

    }

    count = 0;
    int nodeno;
    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        for (j = 0; j < npc[iCell]; ++ j)
        {
            nodeno = c2n[count ++];
            dqmax  = qmax[nodeno] - q[iCell];
            dqmin  = qmin[nodeno] - q[iCell];
            dx     =  xn[nodeno] - xcc[iCell];
            dy     =  yn[nodeno] - ycc[iCell];
            dz     =  zn[nodeno] - zcc[iCell];
            ds     =  sqrt(dx * dx + dy * dy + dz * dz);
            dqNode = dqdx[iCell] * dx + dqdy[iCell] * dy + dqdz[iCell] * dz;

            eps = venkatCoeff * venkatCoeff * venkatCoeff * ds * ds * ds * qAll[IP][iCell] * qAll[IP][iCell] / globalTotalVolume;
            if (dqNode > SMALL)
            {
                tmp          = VenFun(dqmax, dqNode, eps);
                limit[iCell] = MIN(limit[iCell], tmp);
            }
            else if (dqNode < - SMALL)
            {
                tmp          = VenFun(dqmin, dqNode, eps);
                limit[iCell] = MIN(limit[iCell], tmp);
            }
        }
    }

    delete [] qmax;    qmax = NULL;
    delete [] qmin;    qmin = NULL;
}
#pragma warning(disable:4100)
void Limiter::VencatLimiter7(RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, int iVariable)
{
#pragma warning(default:4100)
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal     = nTotalCell + nBoundFace;      
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *xfc   = grid->GetFaceCenterX();
    RDouble *yfc   = grid->GetFaceCenterY();
    RDouble *zfc   = grid->GetFaceCenterZ();
    RDouble *xcc   = grid->GetCellCenterX();
    RDouble *ycc   = grid->GetCellCenterY();
    RDouble *zcc   = grid->GetCellCenterZ();
    RDouble *volume   = grid->GetCellVolume();

    int le, re;
    RDouble dx, dy, dz;
    RDouble dqFace, tmp;

    RDouble *dmin = new RDouble[nTotal];
    RDouble *dmax = new RDouble[nTotal];
    //! Find the the differences for q.
    MinMaxDiff(dmin, dmax, q);

    RDouble venkatCoeff = GlobalDataBase::GetDoubleParaFromDB("venkatCoeff");
    RDouble averageVolume = grid->GetAverageVolume();

    //! The following 'eps' is the eps^2 actually in the mathematic formation.
    //! And the 'venkatCoeff' is the coefficient of 'K' in the mathematic formation.
    int dimension = GetDim();
    RDouble n = three;
    RDouble *epsCell = new RDouble[nTotalCell];
    RDouble eps_tmp = venkatCoeff * venkatCoeff * venkatCoeff;
    eps_tmp *= one/pow(pow(averageVolume, one/dimension), n);

    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        epsCell[iCell]  = eps_tmp * q[iCell] * q[iCell];
        epsCell[iCell] *= pow(pow(volume[iCell], one/dimension), n);
    }
    //! The following loop is to calculate the limiter coefficient of each cell on the boundary face.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        //! Get the left and right cell index of iFace.
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        //! Compute the vector from face to left cell center.
        dx = xfc[iFace] - xcc[le];
        dy = yfc[iFace] - ycc[le];
        dz = zfc[iFace] - zcc[le];

        //! Compute the dq from left cell center to face center.
        dqFace  = dqdx[le] * dx + dqdy[le] * dy + dqdz[le] * dz;
        if (dqFace > SMALL)
        {
            //! According to the vencat limiter of Theory manual.
            tmp       = VenFun(dmax[le], dqFace, epsCell[le]);
            limit[le] = MIN(limit[le], tmp);
        }
        else if (dqFace < - SMALL)
        {
            //! According to the vencat limiter of Theory manual.
            tmp       = VenFun(dmin[le], dqFace, epsCell[le]);
            limit[le] = MIN(limit[le], tmp);
        }
    }

    //! The following loop is to calculate the limiter coefficient of each interior cell.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        dx     = xfc[iFace] - xcc[le];
        dy     = yfc[iFace] - ycc[le];
        dz     = zfc[iFace] - zcc[le];
        dqFace = dqdx[le] * dx + dqdy[le] * dy + dqdz[le] * dz;
        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[le], dqFace, epsCell[le]);
            limit[le] = MIN(limit[le], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[le], dqFace, epsCell[le]);
            limit[le] = MIN(limit[le], tmp);
        }

        dx     = xfc[iFace] - xcc[re];
        dy     = yfc[iFace] - ycc[re];
        dz     = zfc[iFace] - zcc[re];
        dqFace = dqdx[re] * dx + dqdy[re] * dy + dqdz[re] * dz;
        if (dqFace > SMALL)
        {
            tmp       = VenFun(dmax[re], dqFace, epsCell[re]);
            limit[re] = MIN(limit[re], tmp);
        }
        else if (dqFace < - SMALL)
        {
            tmp       = VenFun(dmin[re], dqFace, epsCell[re]);
            limit[re] = MIN(limit[re], tmp);
        }
    }

    delete [] dmin;    dmin = NULL;
    delete [] dmax;    dmax = NULL;
    delete [] epsCell;    epsCell = NULL;
}
void Limiter::MinMaxDiff(RDouble *dmin, RDouble *dmax, RDouble *q)
{
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    //! Find the maximum and minimum in the neighbor of each cell.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dmin[iCell] = q[iCell];
        dmax[iCell] = q[iCell];
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int le       = leftCellOfFace[iFace];
        int re       = rightCellOfFace[iFace];

        int bcType = bcRegion->GetBCType();
        if (!IsInterface(bcType) && bcType != PHENGLEI::OVERSET)
        {
            continue;
        }

        dmin[le] = MIN(dmin[le], q[re]);
        dmax[le] = MAX(dmax[le], q[re]);
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le       = leftCellOfFace [iFace];
        int re       = rightCellOfFace[iFace];
        dmin[le] = MIN(dmin[le], q[re]);
        dmax[le] = MAX(dmax[le], q[re]);

        dmin[re] = MIN(dmin[re], q[le]);
        dmax[re] = MAX(dmax[re], q[le]);
    }

    //! Get the maximum and the minimum difference.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dmin[iCell] -= q[iCell];
        dmax[iCell] -= q[iCell];
    }
}

void Init_Limiters(Grid *grid_in, RDouble *limit)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    int uns_limiter;
    GlobalDataBase::GetData("uns_limiter", &uns_limiter, PHINT, 1);

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        limit[iCell] = one;
    }
    
    if (uns_limiter == ILMT_FIRST)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            limit[iCell] = zero;
        }
    }
}

int GetLimiterID(const string &limiter_name)
{
    int limiter_id = ILMT_NOLIM;
    if (limiter_name.substr(0,5) == "nolim")
    {
        limiter_id = ILMT_NOLIM;
    }
    else if (limiter_name.substr(0,3) == "1st")
    {
        limiter_id = ILMT_FIRST;
    }
    else if (limiter_name.substr(0,5) == "barth")
    {
        limiter_id = ILMT_BARTH;
    }
    else if (limiter_name.substr(0,6) == "vencat")
    {
        limiter_id = ILMT_VENCAT;
    }
    else if (limiter_name.substr(0,6) == "minmod")
    {
        limiter_id = ILMT_MINMOD;
    }
    else if (limiter_name.substr(0,7) == "vanleer")
    {
        limiter_id = ILMT_VANLEER;
    }
    else if (limiter_name.substr(0,10) == "minmodtest")
    {
        limiter_id = ILMT_MINMODTEST;
    }
    else if (limiter_name.substr(0,9) == "vanalbada")
    {
        limiter_id = ILMT_VAN_ALBADA;
        if (limiter_name.substr(0,13) == "vanalbada_clz")
        {
            limiter_id = ILMT_VAN_ALBADA_CLZ;
        }
    }
    else if (limiter_name.substr(0,6) == "smooth")
    {
        limiter_id = ILMT_SMOOTH;
    }
    else if (limiter_name.substr(0,6) == "minvan")
    {
        limiter_id = ILMT_MIN_VAN;
    }
    else if (limiter_name.substr(0,9) == "3rdsmooth")
    {
        limiter_id = ILMT_3rd_SMOOTH;
    }
    else if (limiter_name.substr(0,17) == "3rd_minmod_smooth")
    {
        limiter_id = ILMT_3rd_Minmod_SMOOTH;
    }
    else if (limiter_name.substr(0, 8) == "weno3_js")
    {
        limiter_id = ILMT_WENO3_JS;
    }
    else if (limiter_name.substr(0, 12) == "wenn3_prm211")
    {
        limiter_id = ILMT_WENN3_PRM211;
    }
    else if (limiter_name.substr(0, 8) == "wenn3_zm")
    {
        limiter_id = ILMT_WENN3_ZM;
    }
    else if (limiter_name.substr(0, 10) == "wenn3_zes2")
    {
        limiter_id = ILMT_WENN3_ZES2;
    }
    else if (limiter_name.substr(0, 10) == "wenn3_zes3")
    {
        limiter_id = ILMT_WENN3_ZES3;
    }
    else
    {
        TK_Exit::ExceptionExit("this limiter does'nt exist\n");
    }
    return limiter_id;
}

}

