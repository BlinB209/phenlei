#include "PHMpi.h"
#include "AleModel.h"
#include "PHIO.h"
#include "Constants.h"
#include "Math_BasisFunction.h"
#include "PHMatrix.h"
#include "Geo_Grid.h"
using namespace std;

namespace PHSPACE
{
void lmsmcoef(int lmsmflag, int outnstep, RDouble *outerdt, RDouble &xi, RDouble &theta, RDouble &phi)
{
    //! This subroutine calculates the coefficients for the
    //! Dual Time Stepping Method with a variable outer time step.

    //! Note that it may be easily adapted to a three or more time
    //! stepping method, ie Linear Multi Step Methods.
    //! WARNING: these higher order schemes are not available but
    //! they could be implemented with little effort, given examples show 
    //! that.
    //! PLEASE: verify relations linking xi theta and phi before
    //! using such schemes.
    
    //! outerdt(1)= dt for the time step of the outer loop at n.
    //! outerdt(2)= outer time-step at n-1.
    //! outerdt(3)= outer time-step at n-2 (not used yet).
    //! outnstep  = current time-step for the outer loop.
    //! outtimesc = time stepping scheme for the outer loop.
    
    //! xi, theta, phi, coefficients for the dual time stepping.
    
    int order;

    if (lmsmflag == LMSM_IMPEUL)
    {
        //! FIRST ORDER SCHEME.
        //! Backward euler - A stable.
        theta = one;
        phi   = zero;
        xi    = zero;
        order = 1;
    }
    else if (lmsmflag == LMSM_IMPTRZ)
    {
        //! SECOND ORDER SCHEMES.
        //! One-step trapezoidal - A stable.
        theta = half;
        phi   = zero;
        order = 2;
        ordcoef(order,xi,theta,phi,outerdt,2);
    }
    else if (lmsmflag == LMSM_IMPBD2)
    {
        //! Backward differentiation - A stable.
        if (outnstep == 1)
        {
            theta = one;
            phi   = zero;
            xi    = zero;
            order = 1;
        }
        else
        {
            theta = one;
            phi   = zero;
            order = 2;
            ordcoef(order,xi,theta,phi,outerdt,2);
        }
    }
    else if (lmsmflag == LMSM_IMPTHR)
    {
        //! THIRD ORDER SCHEMES / NOT YET IMPLEMENTED.
        //! Third-order implicit.
        if (outnstep == 1)
        {
            theta = one;
            phi   = zero;
            xi    = zero;
            order = 1;
        }
        else if (outnstep == 2)
        {
            theta = one;
            phi   = zero;
            order = 2;
            ordcoef(order,xi,theta,phi,outerdt,2);
        }
        else
        {
            phi   = zero;
            order = 3;
            ordcoef(order,xi,theta,phi,outerdt,3);
        }
    }
    else if (lmsmflag == LMSM_IMPFOR)
    {
        //! FOURTH ORDER SCHEME / NOT YET IMPLEMENTED.
        //! Milne's scheme.
        if (outnstep == 1)
        {
            theta = one;
            phi   = zero;
            xi    = zero;
        }
        else if (outnstep == 2)
        {
            theta = one;
            phi   = zero;
            order = 2;
            ordcoef(order,xi,theta,phi,outerdt,2);
        }
        else if (outnstep == 3)
        {
            phi   = zero;
            order = 3;
            ordcoef(order,xi,theta,phi,outerdt,3);
        }
        else
        {
            order = 4;
            ordcoef(order,xi,theta,phi,outerdt,1);
        }
    }
    else if (lmsmflag == LMSM_NODUAL)
    {
        //! No linear multi step method.
        xi    = - one;
        theta =   one;
        phi   =   zero;
    }
}

void ordcoef(int order, RDouble &xi, RDouble &theta, RDouble &phi, RDouble *outerdt, int iact)
{
    //! This routine returns the value of the missing coefficient to 
    //! ensure the required precision of the time integration scheme.
    RDouble a;
    
    a = outerdt[1] / (outerdt[0] + SMALL);

    if (order == 2)
    {
        if (iact == 1)
        {
            //! xi,theta given.
            phi   = (half + half * xi * (1 + a) - theta) / (a + SMALL);
        }
        else if (iact == 2)
        {
            //! theta,phi given.
            xi    = (two * theta + two * a * phi - one) / (one + a);
        }
        else if (iact == 3)
        {
            //! phi,xi given.
            theta = half + half * xi * (one + a) - a * phi;
        }
    }
    else if (order == 3)
    {
        if (iact == 1)
        {
            //! xi given.
            phi   = (one + (one + three * a + two * a * a) * xi) / 
                    (six * a * (one + a) + SMALL);
            theta = half + half * (one + a) * xi - a * phi;
        }
        else if (iact == 2)
        {
            //! theta given.
            phi = ((one + three * a + two * a * a) * theta - a - a * a) / 
                  (a * (two + three * a + a * a) + SMALL);
            xi  = (two * theta + two * a * phi - one) / (one + a);
        }
        else if (iact == 3)
        {
            //! phi given.
            xi    = (- one + six * a * (one + a) * phi) / 
                    (one + three * a + two * a * a);
            theta = half + half * (one + a) * xi - a * phi;
        }
    }
    else if (order == 4)
    {
        if (iact == 1)
        {
            phi   = - half / (one + a + a * a);
            xi    = two * (one + two * a) / (one + a);
            theta = a * (two + three * a) * phi - a * (one + a) * xi;
        }
    }
}

void chckmovmt(int *ialepar, RDouble **angle0, RDouble **coord0, RDouble **angle, RDouble **coord,RDouble *outerdt, 
               int lmsmflag, int outnstep)
{
    //! time      = current time for unsteady calculations.
    //! outerdt   = dt for the time step of the outer loop.
    //! outnstep  = current time-step for the outer loop.
    //! outtimesc = time stepping scheme for the outer loop.
    //! ialepar   = array of integer.
    //! ralepar   = array of real.
    //! u0,v0,w0  = free stream velocity components.

    //! angle     = current Euler's angles at time n, n-1, n-2.
    //! coord     = current position of the reference point at time n, n-1, n-2.

    //! angle0    = current Euler's angles and their time derivatives.
    //! coord0    = current position of the reference point and
    //!             its time derivatives.

    //! Global definitions.
    int itrans = ialepar[0];
    int irotat = ialepar[1];
    int ipdfm  = ialepar[2];

    RDouble dt     = outerdt[0];
    RDouble dtold  = outerdt[1];
    
    //! Compute Euler's angles and the reference position
    //! coordinates and their time derivatives.
    //! These expressions are second order accurate only
    //! for outnstep greater than 1.

    RDouble xi,theta,phi;
    RDouble fac1,fac2,fac3,fac4,fac5,fac6;
    
    lmsmcoef(lmsmflag,outnstep,outerdt,xi,theta,phi);

    fac3 =        xi  / (dtold + SMALL);
    fac1 = (one + xi) / (dt    + SMALL);
    fac2 =    - (fac1 + fac3);
    fac4 =        two / (dt * dtold + SMALL);
    fac6 =         dt / (dt + dtold + SMALL);
    fac5 =        one - fac6;
    
    if (outnstep == 1) fac4 = 0.0;

    //! This is the case of a predefined movement.
    if (ipdfm >= 1 && ipdfm <= 3)
    {
        return;
    }
    else if (ipdfm == 4)
    {
        return;
    }
    else if (ipdfm == 5)
    {
        return;
    }
    else if (ipdfm == 6)
    {
        return;
    }
    
    //! Translation : a constant translation velocity is considered.
    if (itrans == 1)
    {
        for (int m = 0; m < 3; ++ m)
        {
            coord0[m][0] = coord0[m][0] + coord0[m][1] * dt;
            coord0[m][2] = 0.0;
            coord [m][0] = coord0[m][0];
        }
    }
    //! Translation : the reference point coordinates are prescribed.
    else if (itrans == 2)
    {
        //deftrans(coord,time,outerdt,ialepar,ralepar);
        for (int m = 0; m < 3; ++ m)
        {
            coord0[m][0] = coord[m][0];
            coord0[m][1] = fac1 * coord[m][0] + fac2 * coord[m][1] + fac3 * coord[m][2];
            coord0[m][2] = fac4 * (fac5 * coord[m][0] - coord[m][1] + fac6 * coord[m][2]);
        }
    }
    //! Translation : the movement is enforced in another way.
    else if (itrans == 3)
    {
        //! Not yet implemented!
        cout << "this option is not yet implemented\n";
        return;
    }

    //! Rotation : a constant angular velocity is considered.
    if (irotat == 1)
    {
        for (int m = 0; m < 3; ++ m)
        {
            angle0[m][0] = angle0[m][0] + angle0[m][1] * dt;
            angle0[m][2] = 0.0;
            angle [m][0] = angle0[m][0];
        }
    }
    //! Rotation : the values of Euler's angles are analytically prescribed.
    else if (irotat == 2)
    {
        //defrotat(angle,time,outerdt,ialepar,ralepar);
        for (int m = 0; m < 3; ++ m)
        {
            angle0[m][0] = angle[m][0];
            angle0[m][1] = fac1 * angle[m][0] + fac2 * angle[m][1] + fac3 * angle[m][2];
            angle0[m][2] = fac4 * (fac5 * angle[m][0] - angle[m][1] + fac6 * angle[m][2]);
        }
    }
    //! Rotation : the movement is enforced in another way.
    else if (irotat == 3)
    {
        //! Not yet implemented!
        cout << "this option is not yet implemented\n";
        return;
    }
}

//void deftrans(RDouble **coord, RDouble time, RDouble *outerdt, int *ialepar, RDouble *ralepar)
//{
//    //! This routine enables to write the reference position coordinates
//    //! (coord(.,1)) as a function of time or outerdt.
//    
//    //! coord(1,1) = ....
//    //! coord(2,1) = ....
//    //! coord(3,1) = ....
//    return;
//}

void chckdeform()
{
    return;
}

RDouble * AleModel::euler_angle  = 0;
RDouble **AleModel::trans_ref    = 0;
RDouble **AleModel::xyz_ref      = 0;
RDouble **AleModel::rotate_ref   = 0;
RDouble **AleModel::angle_ref    = 0;
RDouble **AleModel::matrix       = 0;
AleParameter *AleModel::param = 0;

class AleInitializingClass
{
public:
    AleInitializingClass();
    ~AleInitializingClass();
};

AleInitializingClass::AleInitializingClass()
{
    AleModel::Allocate();
}

AleInitializingClass::~AleInitializingClass()
{
    AleModel::Free();
}

AleInitializingClass ale_local_initialize;

AleParameter::AleParameter()
{
    outerdt = new RDouble[3];
}

AleParameter::~AleParameter()
{
    delete [] outerdt;
}

void AleParameter::SetLMSMFlag()
{
    string outtimesc;
    GlobalDataBase::GetData("outtimesc", &outtimesc, PHSTRING, 1);

    if (outtimesc.substr(0,6) == "impeul")
    {
        lmsmflag = LMSM_IMPEUL;
    }
    else if (outtimesc.substr(0,6) == "imptrz")
    {
        lmsmflag = LMSM_IMPTRZ;
    }
    else if (outtimesc.substr(0,6) == "impbd2")
    {
        lmsmflag = LMSM_IMPBD2;
    }
    else if (outtimesc.substr(0,6) == "impthr")
    {
        lmsmflag = LMSM_IMPTHR;
    }
    else if (outtimesc.substr(0,6) == "impfor")
    {
        lmsmflag = LMSM_IMPFOR;
    }
    else if (outtimesc.substr(0,6) == "nodual")
    {
        lmsmflag = LMSM_NODUAL;
    }
    else
    {
        lmsmflag = LMSM_NODUAL;
    }
}

void AleParameter::Init()
{
    //! At this time,we can not get the value of outnstep ,temporarily set to 0.

    int outnstep = 0;
    RDouble dtau = GlobalDataBase::GetDoubleParaFromDB("dtau");
    RDouble timemax = GlobalDataBase::GetDoubleParaFromDB("timemax");
    RDouble dtsave = GlobalDataBase::GetDoubleParaFromDB("dtsave");

    RDouble simutime = zero;

    this->outnstep = outnstep;
    this->dtau     = dtau;
    this->simutime = simutime;
    this->dtsave   = dtsave;
    this->timemax  = timemax;

    this->outerdt[0] = zero;
    this->outerdt[1] = zero;
    this->outerdt[2] = zero;

    SetLMSMFlag();
}

void AleModel::Allocate()
{
    trans_ref  = CreateMatrix();
    xyz_ref    = CreateMatrix();
    rotate_ref = CreateMatrix();
    angle_ref  = CreateMatrix();
    matrix     = CreateMatrix();
    param      = new AleParameter();
}

void AleModel::Free()
{
    DestroyMatrix(trans_ref);
    DestroyMatrix(xyz_ref  );
    DestroyMatrix(rotate_ref);
    DestroyMatrix(angle_ref);
    DestroyMatrix(matrix   );
    delete param;
}

void AleModel::SetMatrix(RDouble ** matrix_x, RDouble ** matrix_y)
{
    for (int i = 0; i < 3; ++ i)
    {
        for (int j = 0; j < 3; ++ j)
        {
            matrix_x[i][j] = matrix_y[i][j];
        }
    }
}

void AleModel::TransposeMatrix(RDouble ** matrix)
{
    for (int i = 0; i < 3; ++ i)
    {
        for (int j = 0; j < i; ++ j)
        {
            SWAP(matrix[i][j],matrix[j][i]);
        }
    }
}

void AleModel::ComputeRelativeRotateMatrix(RDouble ** matrix_now, RDouble ** matrix_old)
{
    //! Create matrix_now and matrix_old.
    RDouble ** matrix_t = CreateMatrix();

    SetMatrix(matrix_t, matrix_old);

    TransposeMatrix(matrix_old);
    MatrixMultiply(matrix_now, matrix_t, matrix);

    DestroyMatrix(matrix_t);
}

void AleModel::ComputeRelativeRotateMatrix()
{
    //! Create matrix_now and matrix_old.
    RDouble ** matrix_now = CreateMatrix();
    RDouble ** matrix_old = CreateMatrix();

    MatrixRotate(angle_ref[0], matrix_now);
    MatrixRotate(angle_ref[1], matrix_old);

    ComputeRelativeRotateMatrix(matrix_now, matrix_old);

    //! Free matrix_now and matrix_old.
    DestroyMatrix(matrix_now);
    DestroyMatrix(matrix_old);
}

void AleModel::RotateTranslation()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    ComputeRelativeRotateMatrix();
    int level = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Grid *grid = GetGrid(iZone, level);
        if (!grid) continue;
        RotateTranslation(grid);
    }
}

void AleModel::RotateTranslation(Grid *grid)
{
    int nTotalNode = grid->GetNTotalNode();
    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble xref_new = xyz_ref[0][0];
    RDouble yref_new = xyz_ref[0][1];
    RDouble zref_new = xyz_ref[0][2];

    RDouble xref_old = xyz_ref[1][0];
    RDouble yref_old = xyz_ref[1][1];
    RDouble zref_old = xyz_ref[1][2];

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        RDouble dx = x[iNode] - xref_old;
        RDouble dy = y[iNode] - yref_old;
        RDouble dz = z[iNode] - zref_old;
        
        x[iNode] = matrix[0][0] * dx + matrix[0][1] * dy + matrix[0][2] * dz + xref_new;
        y[iNode] = matrix[1][0] * dx + matrix[1][1] * dy + matrix[1][2] * dz + yref_new;
        z[iNode] = matrix[2][0] * dx + matrix[2][1] * dy + matrix[2][2] * dz + zref_new;
    }
}

void AleModel::EulerAngle()
{
    MatrixRotate(euler_angle, matrix);
}

void AleModel::MatrixMultiply(RDouble **x, RDouble **y, RDouble **z)
{
    for (int i = 0; i < 3; ++ i)
    {
        for (int j = 0; j < 3; ++ j)
        {
            z[i][j] = zero;
            for (int k = 0; k < 3; ++ k)
            {
                z[i][j] += x[i][k] * y[k][j];
            }
        }
    }
}

void AleModel::MatrixRotate(RDouble *angle, RDouble **matrix)
{
    //! Create matrix1£¬matrix2£¬matrix3£¬matrix4.
    RDouble ** matrix1 = CreateMatrix();
    RDouble ** matrix2 = CreateMatrix();
    RDouble ** matrix3 = CreateMatrix();
    RDouble ** matrix4 = CreateMatrix();

    //! y-z-x.
    MatrixRotateY(angle[0], matrix1);
    MatrixRotateZ(angle[1], matrix2);
    MatrixRotateX(angle[2], matrix3);

    MatrixMultiply(matrix2, matrix1, matrix4);
    MatrixMultiply(matrix3, matrix4, matrix);

    //! Free matrix1£¬matrix2£¬matrix3£¬matrix4.
    DestroyMatrix(matrix1);
    DestroyMatrix(matrix2);
    DestroyMatrix(matrix3);
    DestroyMatrix(matrix4);
}

void AleModel::MatrixRotateX(RDouble angle, RDouble **matrix)
{
    RDouble cosa = cos(angle);
    RDouble sina = sin(angle);

    matrix[0][0] = one ; matrix[0][1] = zero; matrix[0][2] =   zero;
    matrix[1][0] = zero; matrix[1][1] = cosa; matrix[1][2] = - sina;
    matrix[2][0] = zero; matrix[2][1] = sina; matrix[2][2] =   cosa;
}

void AleModel::MatrixRotateY(RDouble angle, RDouble **matrix)
{
    RDouble cosa = cos(angle);
    RDouble sina = sin(angle);

    matrix[0][0] =   cosa; matrix[0][1] = zero; matrix[0][2] =  sina;
    matrix[1][0] =   zero; matrix[1][1] = one ; matrix[1][2] =  zero;
    matrix[2][0] = - sina; matrix[2][1] = zero; matrix[2][2] =  cosa;
}

void AleModel::MatrixRotateZ(RDouble angle, RDouble **matrix)
{
    RDouble cosa = cos(angle);
    RDouble sina = sin(angle);

    matrix[0][0] = cosa; matrix[0][1] = - sina; matrix[0][2] = zero;
    matrix[1][0] = sina; matrix[1][1] =   cosa; matrix[1][2] = zero;
    matrix[2][0] = zero; matrix[2][1] =   zero; matrix[2][2] = one ;
}

RDouble ** AleModel::CreateMatrix()
{
    return NewPointer2<RDouble>(3, 3);
}

void AleModel::DestroyMatrix(RDouble **matrix)
{
    DelPointer2(matrix);
}

RDouble * AleModel::CreateVector()
{
    return new RDouble[3];
}

void AleModel::DestroyVector(RDouble *v)
{
    delete [] v;
}

void AleModel::CrossMultiply(RDouble *a, RDouble *b, RDouble *c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

void AleModel::UpdateGridInfo()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    int level = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Grid *grid = GetGrid(iZone, level);
        if (!grid) continue;
        grid->ComputeMinMaxBox();
        grid->ComputeMetrics();
    }
}

void AleModel::DynamicGrid()
{
    BackUpOldGrid();
    GenerateDynamicGrid();
    ComputeGridVelocity();
    UpdateGridInfo();
}

void AleModel::GenerateDynamicGrid()
{
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    if (isAle)
    {
        //! Rigid body movement.
        CheckMovement();
        RotateTranslation();
    }
    else
    {
        //! Arbitrary movement.
        ReadDynamicGrid();
    }
}

void AleModel::GetGridFileName(string &filename)
{
    filename = "movegrid.fts";
}

void AleModel::BackUpOldGrid()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    int level = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Grid *grid = GetGrid(iZone, level);
        if (!grid) continue;
        grid->BackUpOldGrid();
    }
}

void AleModel::ComputeGridVelocity()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    int level = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Grid *grid = GetGrid(iZone, level);
        if (!grid) continue;

        int nTotalNode = grid->GetNTotalNode();
        RDouble *xt = new RDouble[nTotalNode];
        RDouble *yt = new RDouble[nTotalNode];
        RDouble *zt = new RDouble[nTotalNode];

        GridVerticeVelocity(grid, xt, yt, zt);
        GridSurfaceVelocity(grid, xt, yt, zt);

        delete [] xt;
        delete [] yt;
        delete [] zt;
    }
}

void AleModel::ReadDynamicGrid()
{
    fstream file;
    ios_base::openmode openmode = ios_base::in|ios_base::binary;

    string filename;

    GetGridFileName(filename);

    ParallelOpenFile(file, filename, openmode);

    ReadDynamicGrid(file);

    ParallelCloseFile(file);
}

void AleModel::ReadDynamicGrid(fstream &file)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ReadGridCoordinate(file, iZone);
    }
}

void AleModel::ReadGridCoordinate(fstream &file, int iZone)
{
    using namespace PHMPI;

    int send_proc = 0;
    int recv_proc = 0;

    GetSendRecvProc(iZone, send_proc, recv_proc);

    int myid = GetCurrentProcessorID();

    //! If the process neither sends nor receives, it will return.
    if (myid != send_proc && myid != recv_proc) return;

    DataContainer *cdata = new DataContainer();

    //! Read grid file ,and send to the corresponding process.
    ReadAbstractData(file, cdata, send_proc, recv_proc);

    if (myid == recv_proc)
    {
        Grid *grid = GetGrid(iZone, 0);
        grid->UpdateCoordinate(cdata);
    }

    delete cdata;    cdata = nullptr;
}

void AleModel::GridSurfaceVelocity(Grid *grid, RDouble *xt, RDouble *yt, RDouble *zt)
{
    grid->GridSurfaceVelocity(xt, yt, zt);
}

void AleModel::GridVerticeVelocity(Grid *grid, RDouble *xt, RDouble *yt, RDouble *zt)
{
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    if (isAle)
    {
        //! Rigid body movement.
        GridVerticeVelocityRigid(grid, xt, yt, zt);
    }
    else
    {
        //! Arbitrary movement.
        GridVerticeVelocityGeneral(grid, xt, yt, zt);
    }
}

void AleModel::GridVerticeVelocityGeneral(Grid *grid, RDouble *xt, RDouble *yt, RDouble *zt)
{
    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    SimpleGrid * oldGrid = grid->GetOldGrid();

    RDouble *xold = oldGrid->GetX();
    RDouble *yold = oldGrid->GetY();
    RDouble *zold = oldGrid->GetZ();

    int nTotalNode = grid->GetNTotalNode();

    RDouble dtau;
    GlobalDataBase::GetData("dtau", &dtau, PHDOUBLE, 1);

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        xt[iNode] = (x[iNode] - xold[iNode]) / dtau;
        yt[iNode] = (y[iNode] - yold[iNode]) / dtau;
        zt[iNode] = (z[iNode] - zold[iNode]) / dtau;
    }
}

void AleModel::GridVerticeVelocityRigid(Grid *grid, RDouble *xt, RDouble *yt, RDouble *zt)
{
    //! This subroutine calculates the grid velocity at the cell vertices.
    int nTotalNode = grid->GetNTotalNode();
    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    //! Compute the components of the rotation vector.
    //! Omega (psip,thetap,phip) in the frame (k0,jh,i) where
    //! k0 = -sin theta * i + cos theta sin phi * j + cos theta cos phi * k
    //! jh = cos phi * j - sin phi * k
    //! Omega(p,q,r) in the frame (i,j,k) angle0(1,.) = (psi,theta,phi)
    //! psi is not set since it is not needed.

    RDouble xref_old = xyz_ref[1][0];
    RDouble yref_old = xyz_ref[1][1];
    RDouble zref_old = xyz_ref[1][2];

    RDouble *dr    = CreateVector();
    RDouble *omega = CreateVector();
    RDouble *vel   = CreateVector();

    //! Calculate the grid velocity with the formula dOM =dOP+Omega*PM.
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        RDouble dx = x[iNode] - xref_old;
        RDouble dy = y[iNode] - yref_old;
        RDouble dz = z[iNode] - zref_old;

        dr[0] = dx;
        dr[1] = dy;
        dr[2] = dz;

        CrossMultiply(omega, dr, vel);

        //xt[iNode] = trans_ref[1][0] + vel[0];
        //yt[iNode] = trans_ref[1][1] + vel[1];
        //zt[iNode] = trans_ref[1][2] + vel[2];
        xt[iNode] = zero;
        yt[iNode] = zero;
        zt[iNode] = zero;
    }

    DestroyVector(dr  );
    DestroyVector(omega);
    DestroyVector(vel );
}

void AleModel::CheckMovement()
{
    //! time      = current time for unsteady calculations.
    //! outerdt   = dt for the time step of the outer loop.
    //! outnstep  = current time-step for the outer loop.
    //! outtimesc = time stepping scheme for the outer loop.
    //! ialepar   = array of integer.
    //! ralepar   = array of real.

    //! angle     = current Euler's angles at time n, n-1, n-2.
    //! coord     = current position of the reference point at time n, n-1, n-2.

    //! angle0    = current Euler's angles and their time derivatives.
    //! coord0    = current position of the reference point and
    //!             its time derivatives.

    //! Global definitions.
    int itrans = 1;
    int irotat = 0;

    //! Compute Euler's angles and the reference position
    //! coordinates and their time derivatives.
    //! These expressions are second order accurate only
    //! for outnstep greater than 1.

    RDouble dt = param->dtau;
    
    if (itrans == 1)
    {
        //! Translation : a constant translation velocity is considered.
        for (int m = 0; m < 3; ++ m)
        {
            trans_ref[0][m] += trans_ref[1][m] * dt;
            trans_ref[2][m]  = zero;
        }

        for (int m = 0; m < 3; ++ m)
        {
            xyz_ref[2][m] = xyz_ref[1][m];
            xyz_ref[1][m] = xyz_ref[0][m];
            xyz_ref[0][m] = trans_ref[0][m];
        }
    }
    else if (itrans == 2)
    {
        //! Translation : the reference point coordinates are prescribed.
    }
    else if (itrans == 3)
    {
        //! Translation : the movement is enforced in another way.
        //! Not yet implemented!
        cout << "this option is not yet implemented\n";
        return;
    }

    //! Rotation : a constant angular velocity is considered.
    if (irotat == 1)
    {
        for (int m = 0; m < 3; ++ m)
        {
            rotate_ref[0][m] = rotate_ref[0][m] + rotate_ref[1][m] * dt;
            rotate_ref[2][m] = 0.0;
        }

        for (int m = 0; m < 3; ++ m)
        {
            angle_ref[2][m] = angle_ref[1][m];
            angle_ref[1][m] = angle_ref[0][m];
            angle_ref[0][m] = rotate_ref[0][m];
        }
    }
    //! Rotation : the values of Euler's angles are analytically prescribed.
    else if (irotat == 2)
    {
    }
    //! Rotation : the movement is enforced in another way.
    else if (irotat == 3)
    {
        //! Not yet implemented!
        cout << "this option is not yet implemented\n";
        return;
    }
}

void AleModel::Deform()
{
}

void AleModel::Init()
{
    param->Init();

    //! Return judgment isn't necessary here.
    for (int i = 0; i < 3; ++ i)
    {
        for (int j = 0; j < 3; ++ j)
        {
            xyz_ref   [i][j] = 0.0;
            trans_ref [i][j] = 0.0;
            angle_ref [i][j] = 0.0;
            rotate_ref[i][j] = 0.0;
        }
    }

    for (int i = 0; i < 3; ++ i)
    {
        for (int j = 0; j < 3; ++ j)
        {
            xyz_ref   [i][j] = 0.0;
            trans_ref [i][j] = 0.0;
            angle_ref [i][j] = 0.0;
            rotate_ref[i][j] = 0.0;
        }
    }

    trans_ref[1][0] = - 1.0;
    trans_ref[1][1] =   0.0;
    trans_ref[1][2] =   0.0;
}

void AleModel::MoveGrid()
{
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");

    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    int ialetype = 0;
    GlobalDataBase::GetData("ialetype", &ialetype, PHINT, 1);

    if (!isUnsteady || !isAle) return;

    if (ialetype == 2)
    {
        int ipdfm = 0;
        //ipdfm  = ialepar(3);
        if (ipdfm == 5)
        {
            //simu_compute_force;
        }
    }

    //! codeOfAleModel = 0: no ALE method;
    //!                  1: ALE method for non-moving grids (relative frame);
    //!                  2: ALE method for moving grids (absolute frame);
    //!                  3: ALE method for deforming grids.
    
    //! coord     = current position of the reference point at time n, n-1, n-2.
    //! angle     = current Euler's angles at time n, n-1, n-2.

    
    //! psi, theta and phi angles.
    //! angle0(1,1),angle0(2,1),angle0(3,1)
    //! phi, theta and psi velocity.
    //! angle0(1,2),angle0(2,2),angle0(3,2)
    //! phi, theta and psi acceleration.
    //! angle0(1,3),angle0(2,3),angle0(3,3)
    
    //! x,y,z position.
    //! coord0(1,1),coord0(2,1),coord0(3,1)
    //! x,y,z velocity.
    //! coord0(1,2),coord0(2,2),coord0(3,2)
    //! x,y,z accelerations.
    //! coord0(1,3),coord0(2,3),coord0(3,3)


    int lmsmflag = 1;
    GlobalDataBase::GetData("lmsmflag", &lmsmflag, PHINT, 1);

    int outnstep = 1;
    GlobalDataBase::GetData("outnstep", &outnstep, PHINT, 1);

    int    * ialepar = reinterpret_cast< int    *  > (GlobalDataBase::GetDataPtr("ialepar"));
    RDouble * outerdt = reinterpret_cast< RDouble *  > (GlobalDataBase::GetDataPtr("outerdt"));
    RDouble **angle0  = reinterpret_cast< RDouble ** > (GlobalDataBase::GetDataPtr("angle0"));
    RDouble **coord0  = reinterpret_cast< RDouble ** > (GlobalDataBase::GetDataPtr("coord0"));
    RDouble **angle   = reinterpret_cast< RDouble ** > (GlobalDataBase::GetDataPtr("angle"));
    RDouble **coord   = reinterpret_cast< RDouble ** > (GlobalDataBase::GetDataPtr("coord"));

    if (ialetype == 2)
    {
        //! Check how the movement is enforced.
        chckmovmt(ialepar,angle0,coord0,angle,coord,
            outerdt,lmsmflag,outnstep);
    }
    else if (ialetype >= 3)
    {
        //! Check how the deformation is enforced.
        //! WARNING: chckdeform does nothing!
        chckdeform();
    }

    RDouble *dangle = new RDouble[3];
    RDouble *tcoord = new RDouble[3];

    //! Update position parameters by specifying their temporal increment.
    for (int m = 0; m < 3; ++ m)
    {
        dangle[m] = angle[m][0] - angle[m][1];
        tcoord[m] = coord[m][0] - coord[m][1];
    }

    delete [] dangle;
    delete [] tcoord;

    //! Update Euler's angles and the reference point position.
    for (int m = 0; m < 3; ++ m)
    {
        angle[m][2] = angle[m][1];
        angle[m][1] = angle[m][0];
        
        coord[m][2] = coord[m][1];
        coord[m][1] = coord[m][0];
    }
}

void AleModel::LinearMultiStepMethodsCoeffcient(AleParameter *param)
{
    //! This subroutine calculates the coefficients for the
    //! Dual Time Stepping Method with a variable outer time step.

    //! Note that it may be easily adapted to a three or more time
    //! stepping method, ie Linear Multi Step Methods.
    //! WARNING : these higher order schemes are not available but
    //! they could be implemented with little effort, given examples show 
    //! that.
    //! PLEASE : verify relations linking xi theta and phi before
    //! using such schemes.
    
    //! outerdt(1) = dt for the time step of the outer loop at n.
    //! outerdt(2) = outer time-step at n-1.
    //! outerdt(3) = outer time-step at n-2 (not used yet).
    //! outnstep   = current time-step for the outer loop.
    //! outtimesc  = time stepping scheme for the outer loop.
    
    //! xi, theta, phi, coefficients for the dual time stepping.
    
    if (param->lmsmflag == LMSM_IMPEUL)
    {
        //! FIRST ORDER SCHEME.
        //! Backward euler - A stable.
        param->theta = one;
        param->phi   = zero;
        param->xi    = zero;
        param->order = 1;
    }
    else if (param->lmsmflag == LMSM_IMPTRZ)
    {
        //! SECOND ORDER SCHEMES.
        //! One-step trapezoidal - A stable.
        param->theta = half;
        param->phi   = zero;
        param->xi    = zero;
        param->order = 2;
    }
    else if (param->lmsmflag == LMSM_IMPBD2)
    {
        //! Backward differentiation - A stable.
        if (param->outnstep == 1)
        {
            param->theta = one;
            param->phi   = zero;
            param->xi    = zero;
            param->order = 1;
        }
        else
        {
            param->theta = one;
            param->phi   = zero;
            param->xi    = half;
            param->order = 2;
        }
    }
    else if (param->lmsmflag == LMSM_NODUAL)
    {
        //! No linear multi step method.
        param->theta =   one;
        param->phi   =   zero;
        param->xi    = - one;
    }
}

void AleModel::GetOuterdt(AleParameter *param)
{
    //! This subroutine computes the size of the real time step for the
    //! uter loop, outerdt, when using dual time stepping. It is also
    //! updating the new time level, time.

    // Extract Dual Time Stepping coefficients.
    LinearMultiStepMethodsCoeffcient(param);

    //! Use constant outer time step.
    param->outerdt[0] = param->dtau;
    
    //! Update the "physical" time.
    param->simutime += param->outerdt[0];

    param->outerdt[2] = param->outerdt[1];
    param->outerdt[1] = param->outerdt[0];

}

void AleModel::UnsteadyPreSolve()
{
    //! Pre-work for unsteady.
    int innstep = 0;
    GlobalDataBase::UpdateData("innstep", &innstep, PHINT, 1);

    RDouble physicalTime = 0;
    GlobalDataBase::GetData("physicalTime", &physicalTime, PHDOUBLE, 1);

    RDouble physicalTimeStep = 0;
    GlobalDataBase::GetData("physicalTimeStep", &physicalTimeStep, PHDOUBLE, 1);

    physicalTime += physicalTimeStep;
    GlobalDataBase::UpdateData("physicalTime", &physicalTime, PHDOUBLE, 1);

    //! Pre-work for ALE.
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    if (!isAle) return;

    DynamicGrid();
}

void AleModel::GetOuterSrc()
{
    GetOuterdt(param);
}

}

