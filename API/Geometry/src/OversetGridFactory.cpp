#include <cmath>
#include "OversetGridFactory.h"
#include "Geo_Interface.h"
#include "Geo_Interpoint.h"
#include "Geo_MultiGridInfo_Struct.h"
#include "PHIO.h"
#include "Geo_NodeTopo_Struct.h"
#include "Geo_FaceMetrics_Struct.h"
#include "Geo_CellMetrics_Struct.h"
#include "Geo_DynamicGridMetrics_Struct.h"
#include "Geo_OversetGridTopo_Struct.h"
using namespace std;

namespace PHSPACE
{
int GetCodeOfDigHoles()
{
    return GlobalDataBase::GetIntParaFromDB("codeOfDigHoles");
}

int GetCodeOfTurbulenceModel()
{
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");

    int codeOfTurbulenceModel = 0;
    if (viscousType != INVISCID && viscousType != LAMINAR)
    {
        codeOfTurbulenceModel = 1;
    }
    
    return codeOfTurbulenceModel;
}

//---------------------------------------------------------------------------------
vector< StructGrid * > * structGridTable = NULL;

LIB_EXPORT void GenerateStructGridTable()
{
    using namespace PHMPI;

    int nZones = PHMPI::GetNumberofGlobalZones();
    structGridTable = new vector< StructGrid * >(nZones, NULL);

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Grid *grid = PHSPACE::GetGrid(iZone, 0);
        if (! grid) continue;
        (* structGridTable)[iZone] = PHSPACE::StructGridCast(grid);
    }

    return;
}

LIB_EXPORT void DeleteStructGridTable()
{
    using namespace PHMPI;

    int nZones = PHMPI::GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        (* structGridTable)[iZone] = NULL;
    }

    delete structGridTable; structGridTable = NULL;

    return;
}

LIB_EXPORT StructGrid * GetStructGrid(int i)
{
    return (* structGridTable)[i];
}


LIB_EXPORT HoleGrid::HoleGrid()
{
    nsi = NULL;
    nsj = NULL;
    nsp = NULL;

    x = NULL;
    y = NULL;
    z = NULL;

    box = NULL;

    capsule = NULL;
    zone_st = 0;
    zone_ed = 0;
    zone_at = 0;
    ns = 0;
    np = 0;
}

LIB_EXPORT HoleGrid::~HoleGrid()
{
    delete x; x = NULL;
    delete y; y = NULL;
    delete z; z = NULL;

    delete box; box = NULL;

    delete nsi; nsi = NULL;
    delete nsj; nsj = NULL;
    delete nsp; nsp = NULL;

    uint_t size = capsule->size();
    for (int i = 0; i < size; ++ i)
    {
        delete (* capsule)[i]; (* capsule)[i] = NULL;
    }
    delete capsule; capsule = NULL;
}
LIB_EXPORT void HoleGrid::DigHoles()
{
    GenerateCapsualSpace();

    MeasureBoxAndCapsual();

    ProbeCellsInHole();

    ProbeBorderCell(IN_HOLE_COLOR, SECOND_COLOR);

    ProbeBorderCell(SECOND_COLOR, OVERSET_COLOR);

    return;
}

LIB_EXPORT void HoleGrid::GenerateSurfaceSpace()
{
    nsi = new vector< int >(ns);
    nsj = new vector< int >(ns);
    nsp = new vector< int >(ns);

    return;
}

LIB_EXPORT void HoleGrid::GeneratePointSpace()
{
    x = new vector< RDouble >(np);
    y = new vector< RDouble >(np);
    z = new vector< RDouble >(np);

    box = new vector< RDouble >(6);

    return;
}

void HoleGrid::GenerateCapsualSpace()
{
    int numberOfCapsules = 0;
    for (int is = 0; is < ns; ++ is)
    {
        numberOfCapsules += ((* nsi)[is] - 1) * ((* nsj)[is] - 1);
    }

    capsule = new vector< vector< RDouble > * >(numberOfCapsules, NULL);

    for (int ic = 0; ic < numberOfCapsules; ++ ic)
    {
        (* capsule)[ic] = new vector< RDouble >(6);
    }

    return;
}

void HoleGrid::MeasureBoxAndCapsual()
{
    using namespace PHMPI;
    vector< vector< RDouble > > holeLocalFrame(3, vector< RDouble >(4));
    phy.resize(3, vector< RDouble >(np));

    int iZone = PHSPACE::GetZoneInverseIndex(zone_at);
    int iProc = PHMPI::GetZoneProcessorID(iZone);
    int jProc = PHMPI::GetCurrentProcessorID();

    if (iProc == jProc)
    {
        StructGrid * grid = PHSPACE::GetStructGrid(iZone);
        holeLocalFrame = grid->GetLocalFrame();
    }

    for (int i = 0; i < 3; ++i)
    {
        PH_Bcast(& holeLocalFrame[i][0], 4 * sizeof(RDouble), iProc, 0);
    }

    vector< RDouble > ss(3);
    for (int i = 0; i < np; ++ i)
    {
        PHSPACE::ComputePhysicalCoordinate(holeLocalFrame, (* x)[i], (* y)[i], (* z)[i], ss);
        phy[0][i] = ss[0]; phy[1][i] = ss[1]; phy[2][i] = ss[2];
    }

    for (int i = 0; i < 3; ++ i)
    {
        (* box)[i*2  ] = * min_element(phy[i].begin(), phy[i].end());
        (* box)[i*2+1] = * max_element(phy[i].begin(), phy[i].end());
    }

    vector< int > nd(4);
    vector< vector< RDouble > > nodeCoordinate(3, vector< RDouble >(4));
    int iCapsule = 0;
    for (int m = 0; m < ns; ++ m)
    {
        int spst = (* nsp)[m];
        int js = (* nsj)[m];
        for (int j = 0; j < js - 1; ++ j)
        {
            int is = (* nsi)[m];
            for (int i = 0; i < is - 1; ++ i)
            {
                nd[0] = spst + is * j + i;
                nd[1] = spst + is * j + (i + 1);
                nd[2] = spst + is * (j + 1) + i;
                nd[3] = spst + is * (j + 1) + (i + 1);

                CopyNodeCoordinate(nodeCoordinate, phy, nd);

                ProduceCasing(nodeCoordinate, iCapsule);

                iCapsule += 1;
            }
        }
    }

    return;
}

void HoleGrid::ProbeCellsInHole()
{
    using namespace PHMPI;

    for (int iSmile = zone_st; iSmile <= zone_ed; ++ iSmile)
    {
        int iZone = PHSPACE::GetZoneInverseIndex(iSmile);

        StructGrid * grid = PHSPACE::GetStructGrid(iZone);

        if (! grid) continue;

        grid->ProbeHoleCell(this);
    }
}

void HoleGrid::ProbeBorderCell(int iType, int jType)
{
    using namespace PHMPI;
    int numberOfProcessors = PHMPI::GetNumberOfProcessor();
    vector< int > * friendPoint = new vector< int > [numberOfProcessors];
    vector< int > * friendZone  = new vector< int > [numberOfProcessors];

    vector< int > myPoint, myZone;

    for (int iSmile = zone_st; iSmile <= zone_ed; ++ iSmile)
    {
        int iZone = PHSPACE::GetZoneInverseIndex(iSmile);
        StructGrid * grid = PHSPACE::GetStructGrid(iZone);
        if (! grid) continue;
        grid->ProbeBorderCell(iType, jType, friendPoint, friendZone);
    }
        
    int myid = PHMPI::GetCurrentProcessorID();

    //std_mm
    for (int ip = 0; ip < numberOfProcessors; ++ ip)
    {
        if (ip != myid)
        {
            int nfps = static_cast<int>(friendPoint[ip].size());
            PH_Send(& nfps, sizeof(int), ip, 0);

            if (nfps == 0) continue;

            PH_Send(& friendPoint[ip][0], nfps * sizeof(int), ip, 0);
            PH_Send(& friendZone[ip][0],  nfps * sizeof(int), ip, 0);
        }
        else
        {
            int nfps = 0;
            for (int jp = 0; jp < numberOfProcessors; ++ jp)
            {
                if (jp == myid) continue;
                PH_Receive(& nfps, sizeof(int), jp, 0);

                if (nfps == 0) continue;
                vector< int > carrier(nfps);

                PH_Receive(& carrier[0], nfps * sizeof(int), jp, 0);
                myPoint.insert(myPoint.end(), carrier.begin(), carrier.end());

                PH_Receive(& carrier[0], nfps * sizeof(int), jp, 0);
                myZone.insert(myZone.end(), carrier.begin(), carrier.end());
            }
        }
    }

    delete [] friendPoint; friendPoint = NULL;
    delete [] friendZone ; friendZone  = NULL;

    uint_t numberOfPoints = myPoint.size();
    for (int pCounter = 0; pCounter < numberOfPoints; ++ pCounter)
    {
        int zoneIndex = myZone[pCounter];
        int globalPointLabel = myPoint[pCounter];
        StructGrid * grid = PHSPACE::GetStructGrid(zoneIndex);
        grid->ColorOverlapCell(globalPointLabel, jType);
    }

    return;
}

LIB_EXPORT int HoleGrid::GetCellColor(RDouble x, RDouble y, RDouble z, RDouble far)
{
    vector< RDouble > pst(3), ped(3), xqd(4), yqd(4), zqd(4), pp_ist(6);
    vector< RDouble > xi, yi, zi, xs, ys, zs;

    vector< RDouble > & box = * this->box;
    if (x < box[0] || x > box[1] || y < box[2] || y > box[3] || z < box[4] || z > box[5])
    {
        return GENERIC_COLOR;
    }

    pst[0] = x; ped[0] = x;
    pst[1] = y; ped[1] = far;
    pst[2] = z; ped[2] = z;

    vector< int > nd(4);
    int i_cap = -1;
    for (int m = 0; m < ns; ++ m)
    {
        int spst = (* nsp)[m];

        int js = (* nsj)[m];
        for (int j = 0; j < js - 1; ++ j)
        {
            int is = (* nsi)[m];
            for (int i = 0; i < is - 1; ++ i)
            {
                nd[0] = spst + is * j + i;
                nd[1] = spst + is * j + (i + 1);
                nd[2] = spst + is * (j + 1) + i;
                nd[3] = spst + is * (j + 1) + (i + 1);

                i_cap += 1; //very important part!

                RDouble xqmin = (* (* capsule)[i_cap])[0];
                RDouble xqmax = (* (* capsule)[i_cap])[1];
                if (x < xqmin || x > xqmax) continue;

                RDouble zqmin = (* (* capsule)[i_cap])[4];
                RDouble zqmax = (* (* capsule)[i_cap])[5];
                if (z < zqmin || z > zqmax) continue;

                for (int k = 0; k < 4; ++ k)
                {
                    xqd[k] = phy[0][nd[k]];
                    yqd[k] = phy[1][nd[k]];
                    zqd[k] = phy[2][nd[k]];
                }

                int pCounter = 0;
                Radial_Triangle_Intersect(pst, ped, xqd, yqd, zqd, pp_ist, pCounter);

                for (int ie = 0; ie < pCounter; ++ ie)
                {
                    xi.push_back(pp_ist[ie * 3 + 0]);
                    yi.push_back(pp_ist[ie * 3 + 1]); 
                    zi.push_back(pp_ist[ie * 3 + 2]);
                }
            }
        }
    }

    uint_t microPoints = xi.size();
    if (microPoints != 0)
    {
        xs.push_back(xi[0]); ys.push_back(yi[0]); zs.push_back(zi[0]);

        for (int i = 1; i < microPoints; ++ i)
        {
            uint_t search_size = xs.size();
            bool ifpushback = true;
            for (int j = 0; j < search_size; ++ j)
            {
                RDouble distance = abs(xs[j] - xi[i]) + abs(ys[j] - yi[i]) + abs(zs[j] - zi[i]);
                if (distance < 1.e-13)
                {
                    ifpushback = false;
                    break;
                }
            }

            if (ifpushback) 
            {
                xs.push_back(xi[i]); ys.push_back(yi[i]); zs.push_back(zi[i]);
            }
        }
    }

    uint_t numberOfViewPoints = xs.size();

    if (numberOfViewPoints % 2 == 0)
    {
        return GENERIC_COLOR;
    }
    else
    {
        return IN_HOLE_COLOR;
    }
}

void HoleGrid::CopyNodeCoordinate(vector< vector< RDouble > > & r, vector< vector< RDouble > > & s, vector< int > & node)
{
    for (int j = 0; j < 3; ++ j)
    {
        for (int i = 0; i < 4; ++ i)
        {
            int node_i = node[i];

            r[j][i] = s[j][node_i];
        }
    }
    return;
}

void HoleGrid::ProduceCasing(vector< vector< RDouble > > & nodeCoordinate, int i)
{
    vector< RDouble > & ss = * (* capsule)[i];

    for (int m = 0; m < 3; ++ m)
    {
        ss[2*m  ] = * min_element(nodeCoordinate[m].begin(), nodeCoordinate[m].end());
        ss[2*m+1] = * max_element(nodeCoordinate[m].begin(), nodeCoordinate[m].end());
    }

    return;
}

void HoleGrid::Radial_Triangle_Intersect(vector<RDouble> & pst, vector<RDouble> & ped, vector<RDouble> & xqd, vector<RDouble> & yqd, vector<RDouble> & zqd, vector<RDouble> & pp_ist, int & pCounter)
{
    const int numberOfTriangles = 2;
    int flag = 0;
    vector<RDouble> direct(3), xx(4), bb(4);
    vector< vector<RDouble> > amer(4, vector<RDouble>(4));

    RDouble accumulator = 0.0;
    for (int i = 0; i < 3; ++ i)
    {
        direct[i] = ped[i] - pst[i];
        accumulator += direct[i] * direct[i];
    }
    accumulator = 1.0 / sqrt(accumulator);
    for (int i = 0; i < 3; ++i)
    {
        direct[i] *= accumulator;
    }

    for (int iTriangle = 0; iTriangle < numberOfTriangles; ++iTriangle)
    {
        for (int j = 0; j < 3; ++j)
        {
            amer[0][j] = xqd[j+iTriangle];
            amer[1][j] = yqd[j+iTriangle];
            amer[2][j] = zqd[j+iTriangle];
            amer[3][j] = 1.0;
        }

        for (int i = 0; i < 3; ++i)
        {
            amer[i][3] = -direct[i];
            bb[i] = pst[i];
        }
        amer[3][3] = 0.0;
        bb[3] = 1.0;

        GaussElimination(amer, bb, 4, xx, flag);
        if (flag == 0) continue;

        if (xx[0] > 0.0 && xx[1] > 0.0 && xx[2] > 0.0 && xx[3] > 0.0)
        {
            for (int i = 0; i < 3; ++i)
            {
                pp_ist[3 * pCounter + i] = pst[i] + direct[i] * xx[3];
            }
            pCounter += 1;
        }
    }

    return;
}

void HoleGrid::GaussElimination(vector< vector< RDouble > > & a, vector< RDouble > & b, int n, vector< RDouble > & x, int & flag)
{
    int is = 0;
    flag = 1;
    vector< int > js(n);

    for (int k = 0; k < n-1; ++k)
    {
        RDouble d = 0.0;
        for (int i = k; i < n; ++i)
        {
            for (int j = k; j < n; ++j)
            {
                if (abs(a[i][j]) > d)
                {
                    d = abs(a[i][j]);
                    js[k] = j;
                    is = i;
                }
            }
        }

        if (d < SMALL)
        {
            flag = 0;
        }
        else
        {
            if (js[k] != k)
            {
                for (int i = 0; i < n; ++i)
                {
                    RDouble t = a[i][k];
                    a[i][k] = a[i][js[k]];
                    a[i][js[k]] = t;
                }
            }

            if (is != k)
            {
                for (int j = k; j < n; ++j)
                {
                    RDouble t = a[k][j];
                    a[k][j] = a[is][j];
                    a[is][j] = t;
                }

                RDouble t = b[k];
                b[k] = b[is];
                b[is] = t;
            }
        }

        if (flag == 0) return;

        RDouble sfn = 1.0/a[k][k];
        for (int j = k+1; j < n; ++j)
        {
            a[k][j] *= sfn;
        }

        b[k] *= sfn;

        for (int i = k+1; i < n; ++i)
        {
            for (int j = k+1; j < n; ++j)
            {
                a[i][j] -= a[i][k]*a[k][j];
            }
            b[i] -= a[i][k]*b[k];
        }
    }

    if (abs(a[n-1][n-1]) < SMALL)
    {
        flag = 0;
        return;
    }

    x[n-1] = b[n-1]/a[n-1][n-1];
    for (int i = n-2; i >= 0; --i)
    {
        RDouble t = 0;
        for (int j = i+1; j < n; ++j)
        {
            t += a[i][j]*x[j];
        }
        x[i] = b[i]-t;
    }

    js[n-1] = n-1;
    for (int k = n-1; k >= 0; --k)
    {
        if (js[k] != k)
        {
            RDouble t = x[k];
            x[k] = x[js[k]];
            x[js[k]] = t;
        }
    }

    return;
}

LIB_EXPORT Stainer::Stainer()
{
    holeGridContainer = NULL;
    numberOfHoles = 0;
}

LIB_EXPORT Stainer::~Stainer()
{
    for (int iHole = 0; iHole < numberOfHoles; ++ iHole)
    {
        delete (* holeGridContainer)[iHole]; (* holeGridContainer)[iHole] = NULL;
    }
    delete holeGridContainer; holeGridContainer = NULL;
}

LIB_EXPORT void Stainer::ReadHoleData()
{
    using namespace PHMPI;
    if (! PHSPACE::GetCodeOfDigHoles()) return;

    int serverTmp = GetServerProcessorID();
    fstream file;
    string fileName = GlobalDataBase::GetStrParaFromDB("holeFullFileName");
    ParallelOpenFile(file, fileName, ios_base::in|ios_base::binary);

    PH_Read_Bcast(file, & numberOfHoles, sizeof(int), serverTmp);
    holeGridContainer = new vector< HoleGrid* >(numberOfHoles, NULL);

    vector< int > holeZoneStart(numberOfHoles), holeZoneEnd(numberOfHoles), holeZoneAttached(numberOfHoles), numberOfHoleSurfaces(numberOfHoles), numberOfHolePoints(numberOfHoles);
    PH_Read_Bcast(file, & holeZoneStart[0],        numberOfHoles * sizeof(int), serverTmp);
    PH_Read_Bcast(file, & holeZoneEnd[0],          numberOfHoles * sizeof(int), serverTmp);
    PH_Read_Bcast(file, & holeZoneAttached[0],     numberOfHoles * sizeof(int), serverTmp);
    PH_Read_Bcast(file, & numberOfHoleSurfaces[0], numberOfHoles * sizeof(int), serverTmp);
    PH_Read_Bcast(file, & numberOfHolePoints[0],   numberOfHoles * sizeof(int), serverTmp);

    for (int iHole = 0; iHole < numberOfHoles; ++ iHole)
    {
        HoleGrid* holeGrid = new HoleGrid();
        holeGrid->SetZoneStart      (holeZoneStart[iHole]       );
        holeGrid->SetZoneEnd        (holeZoneEnd[iHole]         );
        holeGrid->SetZoneAttached   (holeZoneAttached[iHole]    );
        holeGrid->SetNumberOfSurface(numberOfHoleSurfaces[iHole]);
        holeGrid->SetNumberOfPoints (numberOfHolePoints[iHole]  );

        holeGrid->GenerateSurfaceSpace();
        PH_Read_Bcast(file, holeGrid->GetNSI(), numberOfHoleSurfaces[iHole] * sizeof(int), serverTmp);
        PH_Read_Bcast(file, holeGrid->GetNSJ(), numberOfHoleSurfaces[iHole] * sizeof(int), serverTmp);
        PH_Read_Bcast(file, holeGrid->GetNSP(), numberOfHoleSurfaces[iHole] * sizeof(int), serverTmp);

        holeGrid->GeneratePointSpace();
        PH_Read_Bcast(file, holeGrid->GetX(), numberOfHolePoints[iHole] * sizeof(RDouble), serverTmp);
        PH_Read_Bcast(file, holeGrid->GetY(), numberOfHolePoints[iHole] * sizeof(RDouble), serverTmp);
        PH_Read_Bcast(file, holeGrid->GetZ(), numberOfHolePoints[iHole] * sizeof(RDouble), serverTmp);

        (* holeGridContainer)[iHole] = holeGrid;
    }

    ParallelCloseFile(file);
}

LIB_EXPORT void Stainer::ProbeOversetCells()
{
    PaintGenericColor();

    ProbeInternalOverlapCells();

    ProbeExternalOverlapCells();
}

void Stainer::PaintGenericColor()
{
    using namespace PHMPI;

    int nZones = PHMPI::GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        StructGrid * grid = (* structGridTable)[iZone];
        
        if (! grid) continue;

        int * iBlank = grid->GetCellTypeContainer();
        int numberOfCells       = grid->GetNTotalCell();
        for (int iCell = 0; iCell < numberOfCells; ++iCell)
        {
            iBlank[iCell] = GENERIC_COLOR;
        }
    }

    return;
}

void Stainer::ProbeInternalOverlapCells()
{
    using namespace PHMPI;

    if (! PHSPACE::GetCodeOfDigHoles()) return;
    for (int iHole = 0; iHole < numberOfHoles; ++ iHole)
    {
        HoleGrid* holeGrid = (* holeGridContainer)[iHole];
        holeGrid->DigHoles();
    }
    EraseSecondColor();

    int numberOfInternalOversetCells = 0;
    int numberOfCellsInHole          = 0;

    int nZones = PHMPI::GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        StructGrid * grid = PHSPACE::GetStructGrid(iZone);
        if (! grid) continue;
        grid->GatherCellsByColor(numberOfCellsInHole, IN_HOLE_COLOR);
        grid->GatherCellsByColor(numberOfInternalOversetCells, OVERSET_COLOR);
    }

    cout << "number of internal overset cells = " << numberOfInternalOversetCells << endl;
    cout << "number of cells in hole          = " << numberOfCellsInHole << endl;

    return;
}

void Stainer::ProbeExternalOverlapCells()
{
    using namespace PHMPI;

    int borderCellCounter = 0;
    int nZones = PHMPI::GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        StructGrid * grid = PHSPACE::GetStructGrid(iZone);

        if (! grid) continue;

        grid->ProbeExternalOverlapCells(borderCellCounter);
    }

    cout << "number of external overset cells = " << borderCellCounter << endl;

    return;
}

void Stainer::EraseSecondColor()
{
    using namespace PHMPI;
    
    int nZones = PHMPI::GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        StructGrid * grid = PHSPACE::GetStructGrid(iZone);
        if (! grid) continue;
        grid->ResetCellType();
    }

    return;
}


LIB_EXPORT BackgroundTree :: BackgroundTree ()
{
    unitCellCapsule = NULL;

    inf_plane = NULL;
    sup_plane = NULL;

    zone_location = NULL;
    i_location    = NULL;
    j_location    = NULL;
    k_location    = NULL;

    fruitLevel = NULL;
    fruitInheritMapping = NULL;

    isEmptyBinaryTree = false;
    zone_st = 0;
    zone_ed = 0;
    numberOfUnitCells = 0;
}

LIB_EXPORT BackgroundTree ::~BackgroundTree ()
{
    delete [] unitCellCapsule; unitCellCapsule = NULL;

    delete [] inf_plane; inf_plane = NULL;
    delete [] sup_plane; sup_plane = NULL;

    delete [] zone_location; zone_location = NULL;
    delete [] i_location; i_location = NULL;
    delete [] j_location; j_location = NULL;
    delete [] k_location; k_location = NULL;

    delete [] fruitLevel; fruitLevel = NULL;
    delete [] fruitInheritMapping; fruitInheritMapping = NULL;
}

LIB_EXPORT void BackgroundTree ::Generate(int zone_st, int zone_ed)
{
    SetZoneViewRange(zone_st, zone_ed);

    ComputeNumberOfUnitCells();

    if (IsEmpty()) return;

    GenerateSourceSpace();

    ComputeUnitCellCapsule();

    InitFruitLevelAndInheritMapping();

    ComputeRootSpaceParameter();

    AppendNewFruit();

    return;
}

void BackgroundTree ::ComputeNumberOfUnitCells()
{
    using namespace PHMPI;

    int counter = 0;
    for (int iSmile = zone_st; iSmile <= zone_ed; ++iSmile)
    {
        int iZone = PHSPACE::GetZoneInverseIndex(iSmile);

        StructGrid * grid = PHSPACE::GetStructGrid(iZone);

        if (! grid) continue;

        counter += grid->GetNTotalNode();
    }

    numberOfUnitCells = counter;// unitCell(i,j,k)<-->node(i,j,k)

    if (numberOfUnitCells == 0) isEmptyBinaryTree = true; //very important part!

    return;
}

void BackgroundTree ::GenerateSourceSpace()
{
    zone_location = new int [numberOfUnitCells];
    i_location    = new int [numberOfUnitCells];
    j_location    = new int [numberOfUnitCells];
    k_location    = new int [numberOfUnitCells];

    unitCellCapsule = new vector< RDouble >[numberOfUnitCells];
    inf_plane = new vector< RDouble >[numberOfUnitCells];
    sup_plane = new vector< RDouble >[numberOfUnitCells];

    fruitLevel    = new int [numberOfUnitCells];
    fruitInheritMapping = new vector< int >[numberOfUnitCells];

    for (int iUnitCell = 0; iUnitCell < numberOfUnitCells; ++ iUnitCell)
    {
        unitCellCapsule[iUnitCell].resize(6);
        inf_plane[iUnitCell].resize(6);
        sup_plane[iUnitCell].resize(6);
        fruitInheritMapping[iUnitCell].resize(3);
    }
    return;
}

void BackgroundTree ::ComputeUnitCellCapsule()
{
    using namespace PHMPI;
    int iUnitCell = 0;
    for (int iSmile = zone_st; iSmile <= zone_ed; ++ iSmile)
    {
        int iZone = PHSPACE::GetZoneInverseIndex(iSmile);

        StructGrid * grid = PHSPACE::GetStructGrid(iZone);

        if (! grid) continue;

        grid->ComputeUnitCellCapsule(this, iUnitCell);
    }

    return;
}

LIB_EXPORT void BackgroundTree ::SetUnitCellParameter(int iUnitCell, vector< RDouble > & xx, vector< RDouble > & yy, vector< RDouble > & zz)
{
    unitCellCapsule[iUnitCell][0] = * min_element(xx.begin(), xx.end());
    unitCellCapsule[iUnitCell][1] = * min_element(yy.begin(), yy.end());
    unitCellCapsule[iUnitCell][2] = * min_element(zz.begin(), zz.end());

    unitCellCapsule[iUnitCell][3] = * max_element(xx.begin(), xx.end());
    unitCellCapsule[iUnitCell][4] = * max_element(yy.begin(), yy.end());
    unitCellCapsule[iUnitCell][5] = * max_element(zz.begin(), zz.end());

    return;
}

void BackgroundTree ::InitFruitLevelAndInheritMapping()
{
    for (int iUnitCell = 0; iUnitCell < numberOfUnitCells; ++iUnitCell)
    {
        fruitLevel[iUnitCell] = 0;

        for (int j = 0; j < 3; ++j)
        {
            fruitInheritMapping[iUnitCell][j] = 0;
        }
    }

    return;
}

void BackgroundTree ::ComputeRootSpaceParameter()
{
    const int ROOT_NODE_INDEX = 0;

    vector< RDouble > infimumOfRootSpace(3), supremumOfRootSpace(3);
    vector< RDouble > carrier(numberOfUnitCells);

    for (int m = 0; m < 3; ++ m)
    {
        for (int iUnitCell = 0; iUnitCell < numberOfUnitCells; ++ iUnitCell)
        {
            carrier[iUnitCell] = unitCellCapsule[iUnitCell][m];
        }
        infimumOfRootSpace[m] = * min_element(carrier.begin(), carrier.end());

        for (int iUnitCell = 0; iUnitCell < numberOfUnitCells; ++ iUnitCell)
        {
            carrier[iUnitCell] = unitCellCapsule[iUnitCell][m + 3];
        }
        supremumOfRootSpace[m] = * max_element(carrier.begin(), carrier.end());
    }

    for (int m = 0; m < 3; ++ m)
    {
        inf_plane[ROOT_NODE_INDEX][m] = infimumOfRootSpace[m];
        sup_plane[ROOT_NODE_INDEX][m] = supremumOfRootSpace[m];

        inf_plane[ROOT_NODE_INDEX][m + 3] = infimumOfRootSpace[m];
        sup_plane[ROOT_NODE_INDEX][m + 3] = supremumOfRootSpace[m];
    }

    cout << "Root node coordinate: " << endl;
    cout << "(xmin, ymin, zmin) = " << "(" << infimumOfRootSpace[0]  << ", " << infimumOfRootSpace[1]  << ", " << infimumOfRootSpace[2]  << ")" << endl;
    cout << "(xmax, ymax, zmax) = " << "(" << supremumOfRootSpace[0] << ", " << supremumOfRootSpace[1] << ", " << supremumOfRootSpace[2] << ")" << endl;

    return;
}

bool BackgroundTree ::UnitCellCapsuleBelongsToNewNode(vector< RDouble > & capsule, vector< RDouble > & inf_plane, vector< RDouble > & sup_plane)
{
    if (inf_plane[0] <= capsule[0] && capsule[0] <= sup_plane[0] && 
        inf_plane[1] <= capsule[1] && capsule[1] <= sup_plane[1] &&
        inf_plane[2] <= capsule[2] && capsule[2] <= sup_plane[2] &&
        inf_plane[3] <= capsule[3] && capsule[3] <= sup_plane[3] &&
        inf_plane[4] <= capsule[4] && capsule[4] <= sup_plane[4] &&
        inf_plane[5] <= capsule[5] && capsule[5] <= sup_plane[5])
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool BackgroundTree ::PointBelongsToUnitCellCapsule(RDouble xx, RDouble yy, RDouble zz, vector< RDouble > & capsule)
{
    const RDouble EPS = 1.e-20;
    if (capsule[0] <= xx+EPS && xx <= capsule[3]+EPS &&
         capsule[1] <= yy+EPS && yy <= capsule[4]+EPS &&
         capsule[2] <= zz+EPS && zz <= capsule[5]+EPS)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool BackgroundTree ::PointBelongsToNewNode(RDouble xx, RDouble yy, RDouble zz, int kNode)
{
    if (kNode == 0) return false;

    if (inf_plane[kNode][0] <= xx && xx <= sup_plane[kNode][3] &&
         inf_plane[kNode][1] <= yy && yy <= sup_plane[kNode][4] &&
         inf_plane[kNode][2] <= zz && zz <= sup_plane[kNode][5])
    {
        return true;
    }
    else
    {
        return false;
    }
}

void BackgroundTree ::ProduceNewLeftNode(int kNode, vector< RDouble > & inf_left, vector< RDouble > & sup_left)
{
    int j = (fruitLevel[kNode] + 1) % 6;

    for (int i = 0; i < 6; ++ i)
    {
        inf_left[i] = inf_plane[kNode][i];
        if (i == j)
        {
            sup_left[i] = 0.5 * (inf_plane[kNode][i] + sup_plane[kNode][i]);
        }
        else
        {
            sup_left[i] = sup_plane[kNode][i];
        }
    }

    return;
}

void BackgroundTree ::ProduceNewRightNode(int kNode, vector< RDouble > & inf_right, vector< RDouble > & sup_right)
{
    int j = (fruitLevel[kNode] + 1) % 6;
    for (int i = 0; i < 6; ++ i)
    {
        if (i == j)
        {
            inf_right[i] = 0.5 * (inf_plane[kNode][i] + sup_plane[kNode][i]);
        }
        else
        {
            inf_right[i] = inf_plane[kNode][i];
        }
        sup_right[i] = sup_plane[kNode][i];
    }

    return;
}

void BackgroundTree ::AppendNewLeftNode(int iUnitCell, int kNode, vector< RDouble > & inf_new_node, vector< RDouble > & sup_new_node)
{
    const int LEFT = 0, PARENT = 1;

    fruitInheritMapping[kNode][LEFT] = iUnitCell;
    fruitInheritMapping[iUnitCell][PARENT] = kNode;
    fruitLevel[iUnitCell] = fruitLevel[kNode] + 1;
    inf_plane[iUnitCell] = inf_new_node;
    sup_plane[iUnitCell] = sup_new_node;

    return;
}

void BackgroundTree ::AppendNewRightNode(int iUnitCell, int kNode, vector< RDouble > & inf_new_node, vector< RDouble > & sup_new_node)
{
    const int PARENT = 1, RIGHT = 2;

    fruitInheritMapping[kNode][RIGHT] = iUnitCell;
    fruitInheritMapping[iUnitCell][PARENT] = kNode;
    fruitLevel[iUnitCell] = fruitLevel[kNode] + 1;
    inf_plane[iUnitCell] = inf_new_node;
    sup_plane[iUnitCell] = sup_new_node;

    return;
}

void BackgroundTree ::ComputeLocalCoordinateByNewtonMethod(RDouble xx, RDouble yy, RDouble zz, vector< RDouble > & ss, vector< int > & node)
{
    uint_t node_size = fruitIndexStorage.size();

    for (int i = 0; i < node_size; ++ i)
    {
        int nodeIndex = fruitIndexStorage[i];
        int iZone = zone_location[nodeIndex];
        int iCore = i_location[nodeIndex];
        int jCore = j_location[nodeIndex];
        int kCore = k_location[nodeIndex];

        StructGrid * grid = PHSPACE::GetStructGrid(iZone);
        grid->ComputeLocalCoordinateByNewtonMethod(iCore, jCore, kCore, xx, yy, zz, ss, node);

        if (! ss.empty()) break;
    }

    return;
}

LIB_EXPORT void BackgroundTree ::ProbeNearestPoint(RDouble xx, RDouble yy, RDouble zz, RDouble & distance, int & index)
{
    vector< RDouble > ss;
    vector< int > globalPointIndex;

    for (int iSmile = zone_st; iSmile <= zone_ed; ++ iSmile)
    {
        int iZone = PHSPACE::GetZoneInverseIndex(iSmile);

        StructGrid * grid = PHSPACE::GetStructGrid(iZone);

        if (! grid) continue;

        grid->ProbeNearestPoint(xx, yy, zz, ss, globalPointIndex);
    }

    distance = ss[0];
    index = globalPointIndex[0];

    uint_t numberOfPoints = globalPointIndex.size();
    for (int i = 1; i < numberOfPoints; ++ i)
    {
        if (ss[i] < distance)
        {
            distance = ss[i];
            index = globalPointIndex[i];
        }
    }

    return;
}

void BackgroundTree ::AppendNewFruit()
{
    for (int iUnitCell = 0; iUnitCell < numberOfUnitCells; ++ iUnitCell)
    {
        AppendNewFruit(iUnitCell);
    }

    return;
}

LIB_EXPORT void BackgroundTree ::SearchNodes(OversetCellCollector  * oversetCellCollector, RDouble xx, RDouble yy, RDouble zz, int ip, int iElement)
{
    fruitIndexStorage.clear();

    SearchNodes(xx, yy, zz, 0);

    if (fruitIndexStorage.empty())
    {
        RDouble distance = 0.0;
        int index = -1;
        ProbeNearestPoint(xx, yy, zz, distance, index);
        oversetCellCollector->SetNearestPointParameter(ip, iElement, distance, index);
    }
    else
    {
        vector< int > node;
        vector< RDouble > ss;
        ComputeLocalCoordinateByNewtonMethod(xx, yy, zz, ss, node);
        if (! ss.empty())
        {
            oversetCellCollector->SetOptimalPointParameter(ip, iElement, ss, node);
        }
        else
        {
            RDouble distance = 0.0;
            int index = -1;
            ProbeNearestPoint(xx, yy, zz, distance, index);
            oversetCellCollector->SetNearestPointParameter(ip, iElement, distance, index);
        }
    }

    return;
}

void BackgroundTree ::SearchNodes(RDouble xx, RDouble yy, RDouble zz, int n)
{
    const int LEFT = 0, RIGHT = 2;

    vector< RDouble > & capsule = unitCellCapsule[n];
    if (PointBelongsToUnitCellCapsule(xx, yy, zz, capsule))
    {
        fruitIndexStorage.push_back(n);
    }

    int leftNodeIndex = fruitInheritMapping[n][LEFT];
    if (PointBelongsToNewNode(xx, yy, zz, leftNodeIndex))
    {
        SearchNodes(xx, yy, zz, leftNodeIndex);
    }

    int rightNodeIndex = fruitInheritMapping[n][RIGHT];
    if (PointBelongsToNewNode(xx, yy, zz, rightNodeIndex))
    {
        SearchNodes(xx, yy, zz, rightNodeIndex);
    }

    return;
}

void BackgroundTree ::AppendNewFruit(int iUnitCell)
{
    const int ROOT_NODE_INDEX = 0;
    const int LEFT = 0, RIGHT = 2;

    vector< RDouble > inf_new_node(6), sup_new_node(6);

    if (iUnitCell == ROOT_NODE_INDEX) return;

    int kNode = ROOT_NODE_INDEX;
    vector< RDouble > & capsule = unitCellCapsule[iUnitCell];

    while (true)
    {
        int kLeftNode = fruitInheritMapping[kNode][LEFT];
        if (kLeftNode != 0)
        {
            if (UnitCellCapsuleBelongsToNewNode(capsule, inf_plane[kLeftNode], sup_plane[kLeftNode]))
            {
                kNode = kLeftNode;
                continue;
            }
        }
        else
        {
            ProduceNewLeftNode(kNode, inf_new_node, sup_new_node);

            if (UnitCellCapsuleBelongsToNewNode(capsule, inf_new_node, sup_new_node))
            {
                AppendNewLeftNode(iUnitCell, kNode, inf_new_node, sup_new_node);
                
                break;
            }
        }

        int kRightNode = fruitInheritMapping[kNode][RIGHT];
        if (kRightNode != 0)
        {
            if (UnitCellCapsuleBelongsToNewNode(capsule, inf_plane[kRightNode], sup_plane[kRightNode]))
            {
                kNode = kRightNode;
                continue;
            }
        }
        else
        {
            ProduceNewRightNode(kNode, inf_new_node, sup_new_node);

            if (UnitCellCapsuleBelongsToNewNode(capsule, inf_new_node, sup_new_node))
            {
                AppendNewRightNode(iUnitCell, kNode, inf_new_node, sup_new_node);
                
                break;
            }
        }
    }

    return;
}

LIB_EXPORT OversetCellCollector ::OversetCellCollector ()
{
    localCellIndex = NULL;
    zoneIndex      = NULL;
    xCenter        = NULL;
    yCenter        = NULL;
    zCenter        = NULL;

    spaceLattice   = NULL;
    xLocal         = NULL;
    yLocal         = NULL;
    zLocal         = NULL;
    radius         = NULL;
    zone_st        = 0;
    zone_ed        = 0;
}

LIB_EXPORT OversetCellCollector ::~OversetCellCollector ()
{
    delete [] localCellIndex; localCellIndex = NULL;
    delete [] zoneIndex; zoneIndex = NULL;
    delete [] xCenter; xCenter = NULL;
    delete [] yCenter; yCenter = NULL;
    delete [] zCenter; zCenter = NULL;

    delete [] spaceLattice; spaceLattice = NULL;
    delete [] xLocal; xLocal = NULL;
    delete [] yLocal; yLocal = NULL;
    delete [] zLocal; zLocal = NULL;
    delete [] radius; radius = NULL;
}

void OversetCellCollector ::GenerateServiceSpace()
{
    using namespace PHMPI;
    const int NODE_NUM = 8;
    int numberOfProcessors = PHMPI::GetNumberOfProcessor();

    spaceLattice = new vector< int > [numberOfProcessors];
    xLocal       = new vector< RDouble > [numberOfProcessors];
    yLocal       = new vector< RDouble > [numberOfProcessors];
    zLocal       = new vector< RDouble > [numberOfProcessors];
    radius       = new vector< RDouble > [numberOfProcessors];

    for (int ip = 0; ip < numberOfProcessors; ++ ip)
    {
        int numberOfOversetCells = elementNumberRecorder[ip];

        spaceLattice[ip].resize(numberOfOversetCells * NODE_NUM);

        xLocal[ip].resize(numberOfOversetCells);
        yLocal[ip].resize(numberOfOversetCells);
        zLocal[ip].resize(numberOfOversetCells);
        radius[ip].resize(numberOfOversetCells);
    }
    
    return;
}

LIB_EXPORT void OversetCellCollector ::SearchOptimalUnitCell(BackgroundTree  * backgroundTree)
{
    using namespace PHMPI;

    if (! backgroundTree->IsEmpty())
    {
        int numberOfProcessors = PHMPI::GetNumberOfProcessor();
        for (int ip = 0; ip < numberOfProcessors; ++ ip)
        {
            int numberOfElements = elementNumberRecorder[ip];
            for (int iElement = 0; iElement < numberOfElements; ++ iElement)
            {
                RDouble xx = xCenter[ip][iElement];
                RDouble yy = yCenter[ip][iElement];
                RDouble zz = zCenter[ip][iElement];

                backgroundTree->SearchNodes(this, xx, yy, zz, ip, iElement);
            }
        }
    }
    else
    {
        SetVirtualPointParameter();
    }

    return;
}

void OversetCellCollector ::GenerateTargetSpace()
{
    using namespace PHMPI;

    int numberOfProcessors = PHMPI::GetNumberOfProcessor();

    localCellIndex = new vector< int > [numberOfProcessors];
    zoneIndex      = new vector< int > [numberOfProcessors];

    xCenter = new vector< RDouble > [numberOfProcessors];
    yCenter = new vector< RDouble > [numberOfProcessors];
    zCenter = new vector< RDouble > [numberOfProcessors];

    elementNumberRecorder = vector< int >(numberOfProcessors, 0);

    return;
}

LIB_EXPORT void OversetCellCollector ::RunGeometricAnalysis(int zone_st, int zone_ed)
{
    SetZoneViewRange(zone_st, zone_ed);

    GenerateTargetSpace();

    CollectOversetCells();

    BcastOversetCells();

    GenerateServiceSpace();

    return;
}

void OversetCellCollector ::BcastOversetCells()
{
    using namespace PHMPI;

    int numberOfProcessors = PHMPI::GetNumberOfProcessor();
    for (int ip = 0; ip < numberOfProcessors; ++ ip)
    {
        int numberOfOversetCells = static_cast<int>(localCellIndex[ip].size());
        PH_Bcast(& numberOfOversetCells, sizeof(int), ip, 0);
        elementNumberRecorder[ip] = numberOfOversetCells;

        if (numberOfOversetCells == 0) continue;

        localCellIndex[ip].resize(numberOfOversetCells);
        PH_Bcast(& localCellIndex[ip][0], numberOfOversetCells * sizeof(int), ip, 0);

        zoneIndex[ip].resize(numberOfOversetCells);
        PH_Bcast(& zoneIndex[ip][0], numberOfOversetCells * sizeof(int), ip, 0);

        xCenter[ip].resize(numberOfOversetCells);
        PH_Bcast(& xCenter[ip][0], numberOfOversetCells * sizeof(RDouble), ip, 0);

        yCenter[ip].resize(numberOfOversetCells);
        PH_Bcast(& yCenter[ip][0], numberOfOversetCells * sizeof(RDouble), ip, 0);

        zCenter[ip].resize(numberOfOversetCells);
        PH_Bcast(& zCenter[ip][0], numberOfOversetCells * sizeof(RDouble), ip, 0);
    }

    return;
}

void OversetCellCollector ::CollectOversetCells()
{
    using namespace PHMPI;

    int myid = PHMPI::GetCurrentProcessorID();

    for (int iSmile = zone_st; iSmile <= zone_ed; ++ iSmile)
    {
        int iZone = PHSPACE::GetZoneInverseIndex(iSmile);

        StructGrid * grid = PHSPACE::GetStructGrid(iZone);

        if (! grid) continue;

        grid->CollectOversetCells(this, myid);
    }

    return;
}

void OversetCellCollector ::ProbeNearestPoint(BackgroundTree  * backgroundTree, RDouble xx, RDouble yy, RDouble zz, int ip, int iElement)
{
    const int NODE_NUM = 8;
    RDouble distance = 0.0;
    int globalPointIndex = -1;

    backgroundTree->ProbeNearestPoint(xx, yy, zz, distance, globalPointIndex);

    xLocal[ip][iElement] = 0.0;
    yLocal[ip][iElement] = 0.0;
    zLocal[ip][iElement] = 0.0;

    radius[ip][iElement] = distance;
    
    for (int iCell = 0; iCell < NODE_NUM; ++ iCell)
    {
        spaceLattice[ip][iElement * NODE_NUM + iCell] = globalPointIndex;
    }

    return;
}

LIB_EXPORT void OversetCellCollector ::ResizeSourceIndexContainer()
{
    using namespace PHMPI;

    int myid = PHMPI::GetCurrentProcessorID();

    int numberOfOversetCells = elementNumberRecorder[myid];

    sourceProcessorIndex.resize(numberOfOversetCells);

    int numberOfProcessors = PHMPI::GetNumberOfProcessor();

    vector< RDouble > * trans = new vector< RDouble > [numberOfProcessors];
    for (int ip = 0; ip < numberOfProcessors; ++ ip)
    {
        trans[ip].resize(numberOfOversetCells);
    }

    trans[myid] = radius[myid];

    for (int ip = 0; ip < numberOfProcessors; ++ ip)
    {
        if (ip != myid)
        {
            int numberOfElements = elementNumberRecorder[ip];
            PH_Send(& numberOfElements, sizeof(int), ip, 0);

            if (numberOfElements == 0) continue;
            PH_Send(&radius[ip][0], numberOfElements * sizeof(RDouble), ip, 0);
        }
        else
        {
            int numberOfElements = 0;
            for (int jp = 0; jp < numberOfProcessors; ++ jp)
            {
                if (jp == myid) continue;
                PH_Receive(& numberOfElements, sizeof(int), jp, 0);

                if (numberOfElements == 0) continue;

                PH_Receive(& trans[jp][0], numberOfElements * sizeof(RDouble), jp, 0);
            }
        }
    }

    for (int iElement = 0; iElement < numberOfOversetCells; ++ iElement)
    {
        RDouble distance = trans[0][iElement];
        for (int ip = 1; ip < numberOfProcessors; ++ ip)
        {
            if (trans[ip][iElement] < distance)
            {
                distance = trans[ip][iElement];

                sourceProcessorIndex[iElement] = ip;
            }
        }
    }

    delete [] trans; trans = NULL;

    return;
}


LIB_EXPORT void OversetCellCollector ::CollectOversetCells(int iZone, vector< int > & cell, vector< RDouble > & xx, vector< RDouble > & yy, vector< RDouble > & zz)
{
    using namespace PHMPI;

    int myid = PHMPI::GetCurrentProcessorID();
    uint_t numberOfOversetCells = cell.size();

    zoneIndex[myid].insert(zoneIndex[myid].end(), numberOfOversetCells, iZone);
    localCellIndex[myid].insert(localCellIndex[myid].end(), cell.begin(), cell.end());
    xCenter[myid].insert(xCenter[myid].end(), xx.begin(), xx.end());
    yCenter[myid].insert(yCenter[myid].end(), yy.begin(), yy.end());
    zCenter[myid].insert(zCenter[myid].end(), zz.begin(), zz.end());

    return;
}

LIB_EXPORT void OversetCellCollector ::SetNearestPointParameter(int ip, int iElement, RDouble distance, int index)
{
    const int NODE_NUM = 8;

    xLocal[ip][iElement] = 0.0;
    yLocal[ip][iElement] = 0.0;
    zLocal[ip][iElement] = 0.0;

    radius[ip][iElement] = distance;
    
    for (int iCell = 0; iCell < NODE_NUM; ++ iCell)
    {
        spaceLattice[ip][iElement * NODE_NUM + iCell] = index;
    }

    return;
}

LIB_EXPORT void OversetCellCollector ::SetOptimalPointParameter(int ip, int iElement, vector< RDouble > & ss, vector< int > & node)
{
    const int NODE_NUM = 8;

    xLocal[ip][iElement] = ss[0];
    yLocal[ip][iElement] = ss[1];
    zLocal[ip][iElement] = ss[2];

    //WriteLogFile("ksai = ",ss[0],ss[1],ss[2]);

    radius[ip][iElement] = 0.0;

    for (int iCell = 0; iCell < NODE_NUM; ++ iCell)
    {
        spaceLattice[ip][iElement * NODE_NUM + iCell] = node[iCell];
    }

    return;
}

void OversetCellCollector ::SetVirtualPointParameter()
{
    using namespace PHMPI;

    const int NODE_NUM = 8, VIRTUAL_POINT_LABEL = -1;
    const RDouble FAR = 1.e20, RELATIVE_LOCATION = 0.0;

    int numberOfProcessors = PHMPI::GetNumberOfProcessor();
    for (int ip = 0; ip < numberOfProcessors; ++ ip)
    {
        int numberOfElements = elementNumberRecorder[ip];
        for (int iElement = 0; iElement < numberOfElements; ++ iElement)
        {
            xLocal[ip][iElement] = RELATIVE_LOCATION;
            yLocal[ip][iElement] = RELATIVE_LOCATION;
            zLocal[ip][iElement] = RELATIVE_LOCATION;

            radius[ip][iElement] = FAR;

            for (int iCell = 0; iCell < NODE_NUM; ++ iCell)
            {
                spaceLattice[ip][iElement * NODE_NUM + iCell] = VIRTUAL_POINT_LABEL;
            }
        }
    }

    return;
}

LIB_EXPORT void OversetCellCollector ::InsertOversetParameter(vector< int > * zoneContainer, vector< int > * cellContainer)
{
    using namespace PHMPI;

    int myid = PHMPI::GetCurrentProcessorID();

    zoneContainer->insert(zoneContainer->end(), zoneIndex[myid].begin(),      zoneIndex[myid].end()     );
    cellContainer->insert(cellContainer->end(), localCellIndex[myid].begin(), localCellIndex[myid].end());

    return;
}

LIB_EXPORT void OversetCellCollector ::InsertServiceParameter(vector< RDouble > * ksai, vector< int > * waiter)
{
    using namespace PHMPI;

    const int NODE_NUM = 8;

    int numberOfProcessors = PHMPI::GetNumberOfProcessor();
    int myid = PHMPI::GetCurrentProcessorID();
    int numberOfOversetCells = elementNumberRecorder[myid];

    vector< RDouble > xgroup(numberOfOversetCells), ygroup(numberOfOversetCells), zgroup(numberOfOversetCells);
    vector< int > node_g(numberOfOversetCells * NODE_NUM);

    vector< RDouble > * x_trans = new vector< RDouble > [numberOfProcessors];
    vector< RDouble > * y_trans = new vector< RDouble > [numberOfProcessors];
    vector< RDouble > * z_trans = new vector< RDouble > [numberOfProcessors];
    vector< int > * node_trans = new vector< int >    [numberOfProcessors];

    
    for (int ip = 0; ip < numberOfProcessors; ++ ip)
    {
        x_trans[ip].resize(numberOfOversetCells);
        y_trans[ip].resize(numberOfOversetCells);
        z_trans[ip].resize(numberOfOversetCells);

        node_trans[ip].resize(numberOfOversetCells * NODE_NUM);
    }

    x_trans[myid] = xLocal[myid];
    y_trans[myid] = yLocal[myid];
    z_trans[myid] = zLocal[myid];

    node_trans[myid] = spaceLattice[myid];

    for (int ip = 0; ip < numberOfProcessors; ++ ip)
    {
        if (ip != myid)
        {
            int numberOfElements = elementNumberRecorder[ip];
            PH_Send(& numberOfElements, sizeof(int), ip, 0);

            if (numberOfElements == 0) continue;

            PH_Send(&xLocal[ip][0], numberOfElements * sizeof(RDouble), ip, 0);
            PH_Send(&yLocal[ip][0], numberOfElements * sizeof(RDouble), ip, 0);
            PH_Send(&zLocal[ip][0], numberOfElements * sizeof(RDouble), ip, 0);

            PH_Send(&spaceLattice[ip][0], numberOfElements * NODE_NUM * sizeof(int), ip, 0);
        }
        else
        {
            int numberOfElements = 0;
            for (int jp = 0; jp < numberOfProcessors; ++ jp)
            {
                if (jp == myid) continue;
                PH_Receive(& numberOfElements, sizeof(int), jp, 0);

                if (numberOfElements == 0) continue;

                PH_Receive(& x_trans[jp][0], numberOfElements * sizeof(RDouble), jp, 0);
                PH_Receive(& y_trans[jp][0], numberOfElements * sizeof(RDouble), jp, 0);
                PH_Receive(& z_trans[jp][0], numberOfElements * sizeof(RDouble), jp, 0);

                PH_Receive(& node_trans[jp][0], numberOfElements * NODE_NUM * sizeof(int), jp, 0);
            }
        }
    }


    for (int iElement = 0; iElement < numberOfOversetCells; ++ iElement)
    {
        int ip = sourceProcessorIndex[iElement];

        xgroup[iElement] = x_trans[ip][iElement];
        ygroup[iElement] = y_trans[ip][iElement];
        zgroup[iElement] = z_trans[ip][iElement];

        for (int iCell = 0; iCell < NODE_NUM; ++ iCell)
        {
            node_g[iElement * NODE_NUM + iCell] = node_trans[ip][iElement * NODE_NUM + iCell];
        }
    }

    delete [] x_trans; x_trans = NULL;
    delete [] y_trans; y_trans = NULL;
    delete [] z_trans; z_trans = NULL;
    delete [] node_trans; node_trans = NULL;

    ksai[0].insert(ksai[0].end(), xgroup.begin(), xgroup.end());
    ksai[1].insert(ksai[1].end(), ygroup.begin(), ygroup.end());
    ksai[2].insert(ksai[2].end(), zgroup.begin(), zgroup.end());

    waiter->insert(waiter->end(), node_g.begin(), node_g.end());

    return;
}


LIB_EXPORT OversetStructGrid::OversetStructGrid()
{
    zoneContainer = NULL;
    cellContainer = NULL;

    ksai          = NULL;
    waiter        = NULL;

    nodeSend      = NULL;
    nodeRecv      = NULL;

    field         = NULL;
    dataSend      = NULL;
    dataRecv      = NULL;
    searchTimes   = 0;
}

LIB_EXPORT OversetStructGrid::~OversetStructGrid()
{
    delete zoneContainer; zoneContainer = NULL;
    delete cellContainer; cellContainer = NULL;

    delete [] ksai; ksai   = NULL;
    delete waiter;  waiter = NULL;

    delete [] nodeSend; nodeSend = NULL;
    delete [] nodeRecv; nodeRecv = NULL;

    delete [] field;    field    = NULL;
    delete [] dataSend; dataSend = NULL;
    delete [] dataRecv; dataRecv = NULL;
}

void OversetStructGrid::ProbeOversetCells()
{
    Stainer *stainer = new Stainer();

    stainer->ReadHoleData();

    stainer->ProbeOversetCells();

    delete stainer;

    return;
}

LIB_EXPORT void OversetStructGrid::RunGeometricAnalysis()
{
    ReadBasicParameter();

    ProbeOversetCells();

    GenerateSourceAndTargetSpace();

    ProbeRelativeOversetLocation();
    
    FillNodeIndexContainer();

    return;
}

void OversetStructGrid::ReadBasicParameter()
{
    using namespace PHMPI;
    int serverTmp = PHMPI::GetServerProcessorID();
    int thisip = PHMPI::GetCurrentProcessorID();

    fstream file0;
    string fileName0 = PHSPACE::GlobalDataBase::GetStrParaFromDB("masterFileName");
    ParallelOpenFile(file0, fileName0, ios_base::in);
    if (thisip == serverTmp) file0 >> searchTimes;
    PH_Bcast(& searchTimes, sizeof(int), serverTmp);

    sourceZoneStart.resize(searchTimes); sourceZoneEnd.resize(searchTimes);
    targetZoneStart.resize(searchTimes); targetZoneEnd.resize(searchTimes);

    if (thisip == serverTmp)
    {
        for (int iSearch = 0; iSearch < searchTimes; ++iSearch)
        {
            file0 >> sourceZoneStart[iSearch] >> sourceZoneEnd[iSearch] >> targetZoneStart[iSearch] >> targetZoneEnd[iSearch];
        }
    }

    PH_Bcast(& sourceZoneStart[0], searchTimes * sizeof(int), serverTmp);
    PH_Bcast(& sourceZoneEnd[0],   searchTimes * sizeof(int), serverTmp);
    PH_Bcast(& targetZoneStart[0], searchTimes * sizeof(int), serverTmp);
    PH_Bcast(& targetZoneEnd[0],   searchTimes * sizeof(int), serverTmp);
    ParallelCloseFile(file0);

    int nZones = PHMPI::GetNumberofGlobalZones();
    zoneInverseIndex.resize(nZones);

    fstream file1;
    string fileName1 = PHSPACE::GlobalDataBase::GetStrParaFromDB("zoneInverseFileName");
    ParallelOpenFile(file1, fileName1, ios_base::in|ios_base::binary);
    PH_Read_Bcast(file1, & zoneInverseIndex[0], nZones * sizeof(int), serverTmp);
    ParallelCloseFile(file1);

    return;
}

void OversetStructGrid::ProbeRelativeOversetLocation()
{
    for (int iSearch = 0; iSearch < searchTimes; ++ iSearch)
    {
        BackgroundTree  * backgroundTree = new BackgroundTree ();

        backgroundTree->Generate(sourceZoneStart[iSearch], sourceZoneEnd[iSearch]);

        OversetCellCollector  * oversetCellCollector = new OversetCellCollector ();

        oversetCellCollector->RunGeometricAnalysis(targetZoneStart[iSearch], targetZoneEnd[iSearch]);

        oversetCellCollector->SearchOptimalUnitCell(backgroundTree);

        InsertOversetAndServiceParameter(oversetCellCollector);

        delete oversetCellCollector; oversetCellCollector = NULL;

        delete backgroundTree; backgroundTree = NULL;
    }

    //cout << "ksai[0].size() = " << ksai[0].size() << endl;
    //cout << "ksai[1].size() = " << ksai[1].size() << endl;
    //cout << "ksai[2].size() = " << ksai[2].size() << endl;

    //cout << "waiter->size() = " << waiter->size() << endl;

    //int numberOfOversetCells = ksai[0].size();
    //for (int i = 0; i < numberOfOversetCells; ++ i)
    //{
    //    WriteLogFile("ksai = ", ksai[0][i], ksai[1][i], ksai[2][i]);
    //}

    //for (int i = 0; i < 8 * numberOfOversetCells; ++ i)
    //{
    //    WriteLogFile("waitress = ", (* waiter)[i]);
    //}

    return;
}

void OversetStructGrid::InsertOversetAndServiceParameter(OversetCellCollector  * oversetCellCollector)
{
    oversetCellCollector->ResizeSourceIndexContainer();

    oversetCellCollector->InsertOversetParameter(zoneContainer, cellContainer);

    oversetCellCollector->InsertServiceParameter(ksai, waiter);

    return;
}

void OversetStructGrid::LinearInsertValue(int numberOfEquations, vector< RDouble > & backdrop, RDouble xx, RDouble yy, RDouble zz, vector< RDouble > & newCellValue)
{
    const int NODE_NUM = 8;
    vector< RDouble > a(NODE_NUM), dq(NODE_NUM);

    for (int m = 0; m < numberOfEquations; ++ m)
    {
        for (int iCell = 0; iCell < NODE_NUM; ++ iCell)
        {
            dq[iCell] = backdrop[m * NODE_NUM + iCell ];
        }

        a[0] =  dq[0];
        a[1] = -dq[0] + dq[1];
        a[2] = -dq[0] + dq[2];
        a[3] = -dq[0] + dq[3];
        a[4] =  dq[0] - dq[1] + dq[4] - dq[2];
        a[5] =  dq[0] - dq[1] - dq[3] + dq[5];
        a[6] =  dq[0] - dq[2] - dq[3] + dq[6];
        a[7] = -dq[0] + dq[1] - dq[4] + dq[2] + dq[3] - dq[5] + dq[7] - dq[6];
        
        newCellValue[m] = a[0] + a[1] * xx + a[2] * yy + a[3] * zz + a[4] * xx * yy + a[5] * xx * zz + a[6] * yy * zz + a[7] * xx * yy * zz;
    }

    return;
}

void OversetStructGrid::GenerateSourceAndTargetSpace()
{
    using namespace PHMPI;

    int numberOfProcessors = PHMPI::GetNumberOfProcessor(); 

    zoneContainer = new vector< int >();
    cellContainer = new vector< int >();

    ksai = new vector< RDouble >[3];
    waiter = new vector< int >();

    nodeSend = new vector< int > [numberOfProcessors];
    nodeRecv = new vector< int > [numberOfProcessors];

    dataSend = new vector< RDouble > * [numberOfProcessors];
    dataRecv = new vector< RDouble > * [numberOfProcessors];

    for (int ip = 0; ip < numberOfProcessors; ++ ip)
    {
        dataSend[ip] = NULL;
        dataRecv[ip] = NULL;
    }

    return;
}

void OversetStructGrid::FillNodeIndexContainer()
{
    using namespace PHMPI;

    uint_t numberOfPorts = waiter->size();
    indexGroup.resize(numberOfPorts);
    processGroup.resize(numberOfPorts);

    for (int iPort = 0; iPort < numberOfPorts; ++ iPort)
    {
        int globalNodeIndex = (* waiter)[iPort];

        int zoneIndex = PHMPI::GetZoneIndexAccordingToGlobalPointLabel(globalNodeIndex);

        int processorIndex = PHMPI::GetZoneProcessorID(zoneIndex);

        nodeSend[processorIndex].push_back(globalNodeIndex);

        indexGroup[iPort] = static_cast<int>(nodeSend[processorIndex].size()) - 1;

        processGroup[iPort] = processorIndex;
    }

    int numberOfProcessors = PHMPI::GetNumberOfProcessor();
    int myid = PHMPI::GetCurrentProcessorID();

    nodeRecv[myid] = nodeSend[myid];

    for (int ip = 0; ip < numberOfProcessors; ++ ip)
    {
        if (ip != myid)
        {
            int numberOfElements = static_cast<int>(nodeSend[ip].size());

            PH_Send(& numberOfElements, sizeof(int), ip, 0);

            if (numberOfElements == 0) continue;

            PH_Send(& nodeSend[ip][0], numberOfElements * sizeof(int), ip, 0);
        }
        else
        {
            for (int jp = 0; jp < numberOfProcessors; ++ jp)
            {
                if (jp == myid) continue;

                int numberOfElements;

                PH_Receive(& numberOfElements, sizeof(int), jp, 0);

                if (numberOfElements == 0) continue;

                nodeRecv[jp].resize(numberOfElements);

                PH_Receive(& nodeRecv[jp][0], numberOfElements * sizeof(int), jp, 0);
            }
        }
    }

    return;
}

LIB_EXPORT void OversetStructGrid::FillSolverContainer()
{
    solverNameContainer.clear();
    equationNumberContainer.clear();

    solverNameContainer.push_back("q");

    int nl = PHSPACE::GlobalDataBase::GetIntParaFromDB("nl");
    equationNumberContainer.push_back(nl);

    if (PHSPACE::GetCodeOfTurbulenceModel())
    {
        solverNameContainer.push_back("q_turb");

        int n_turb = PHSPACE::GlobalDataBase::GetIntParaFromDB("n_turb");
        equationNumberContainer.push_back(n_turb);
    }

    return;
}

LIB_EXPORT void OversetStructGrid::InitFieldNew()
{
    using namespace PHMPI;
    
    uint_t numberOfSolvers = solverNameContainer.size();

    int nZones = PHMPI::GetNumberofGlobalZones();

    field = new vector< RDouble4D * > [nZones];

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        field[iZone].resize(numberOfSolvers, NULL);

        StructGrid * grid = PHSPACE::GetStructGrid(iZone);
        if (! grid) continue;

        for (int iSolver = 0; iSolver < numberOfSolvers; ++ iSolver)
        {
            string name = solverNameContainer[iSolver];
            field[iZone][iSolver] = reinterpret_cast< RDouble4D * >(grid->GetDataPtr(name));
        }
    }
    return;
}

LIB_EXPORT void OversetStructGrid::InitField()
{
    using namespace PHMPI;
    
    uint_t numberOfSolvers = solverNameContainer.size();

    int nZones = PHMPI::GetNumberofGlobalZones();

    field = new vector< RDouble4D * > [nZones];

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        field[iZone].resize(numberOfSolvers, NULL);

        Grid *grid = PHSPACE::GetGrid(iZone);

        if (! grid) continue;

        for (int iSolver = 0; iSolver < numberOfSolvers; ++ iSolver)
        {
            string name = solverNameContainer[iSolver];
            field[iZone][iSolver] = static_cast< RDouble4D * >(grid->GetDataPtr(name));
        }
    }
    return;
}

LIB_EXPORT void OversetStructGrid::CommunicateOversetField(int iSolver)
{
    using namespace PHMPI;
    int iLocal, jLocal, kLocal;

    int numberOfEquations = equationNumberContainer[iSolver];
    int numberOfProcessors = PHMPI::GetNumberOfProcessor();

    for (int ip = 0; ip < numberOfProcessors; ++ ip)
    {
        dataSend[ip] = new vector< RDouble >();
        dataRecv[ip] = new vector< RDouble >();
    }

    for (int ip = 0; ip < numberOfProcessors; ++ ip)
    {
        uint_t numberOfElements = nodeRecv[ip].size();
        for (int iElement = 0; iElement < numberOfElements; ++iElement)
        {
            int globalNodeLabel = nodeRecv[ip][iElement];

            int zoneIndex = PHMPI::GetZoneIndexAccordingToGlobalPointLabel(globalNodeLabel);

            StructGrid * grid = PHSPACE::GetStructGrid(zoneIndex);

            grid->GetLocalIndex(globalNodeLabel, iLocal, jLocal, kLocal);

            RDouble4D & myField = * (field[zoneIndex][iSolver]);

            for (int iEquation = 0; iEquation < numberOfEquations; ++ iEquation)
            {
                (* dataSend[ip]).push_back(myField(iLocal, jLocal, kLocal, iEquation));
            }
        }
    }

    int myid = PHMPI::GetCurrentProcessorID();

    * dataRecv[myid] = * dataSend[myid];

    for (int ip = 0; ip < numberOfProcessors; ++ ip)
    {
        if (ip != myid)
        {
            int numberOfElements = static_cast<int>((* dataSend[ip]).size());
            PH_Send(& numberOfElements, sizeof(int), ip, 0);

            if (numberOfElements == 0) continue;

            PH_Send(& (* dataSend[ip])[0], numberOfElements * sizeof(RDouble), ip, 0);
        }
        else
        {
            for (int jp = 0; jp < numberOfProcessors; ++ jp)
            {
                if (jp == myid) continue;

                int numberOfElements;
                PH_Receive(& numberOfElements, sizeof(int), jp, 0);

                if (numberOfElements == 0) continue;

                dataRecv[jp]->resize(numberOfElements);

                PH_Receive(& (* dataRecv[jp])[0], numberOfElements * sizeof(RDouble), jp, 0);
            }
        }
    }

    return;
}

LIB_EXPORT void OversetStructGrid::CorrectOversetField(int iSolver)
{
    using namespace PHMPI;

    int ic, jc, kc;
    const int NODE_NUM = 8;
    int numberOfEquations = equationNumberContainer[iSolver];

    vector< RDouble > backdrop(numberOfEquations * NODE_NUM), newCellValue(numberOfEquations);

    uint_t numberOfOversetCells = cellContainer->size();
    for (int iOversetCell = 0; iOversetCell < numberOfOversetCells; ++ iOversetCell)
    {
        int iZone = (* zoneContainer)[iOversetCell];

        StructGrid * grid = PHSPACE::GetStructGrid(iZone);

        if (! grid) continue; //Maybe no use!

        RDouble4D &myField = * field[iZone][iSolver];

        grid->GetLocalCenter((* cellContainer)[iOversetCell], ic, jc, kc);

        for (int iNeighbor = 0; iNeighbor < NODE_NUM; ++ iNeighbor)
        {
            int neighborCellLabel = NODE_NUM * iOversetCell + iNeighbor;

            int iProcessor = processGroup[neighborCellLabel];

            int iElement = indexGroup[neighborCellLabel];

            for (int iEquation = 0; iEquation < numberOfEquations; ++iEquation)
            {
                backdrop[NODE_NUM * iEquation + iNeighbor] = (* dataRecv[iProcessor])[numberOfEquations * iElement + iEquation];
            }
        }

        LinearInsertValue(numberOfEquations, backdrop, ksai[0][iOversetCell], ksai[1][iOversetCell], ksai[2][iOversetCell], newCellValue);

        for (int iEquation = 0; iEquation < numberOfEquations; ++ iEquation)
        {
            myField(ic + 1, jc + 1, kc + 1, iEquation) = newCellValue[iEquation];
        }
    }

    int numberOfProcessors = PHMPI::GetNumberOfProcessor();
    for (int ip = 0; ip < numberOfProcessors; ++ ip)
    {
        delete dataSend[ip]; dataSend[ip] = NULL;
        delete dataRecv[ip]; dataRecv[ip] = NULL;
    }

    return;
}

OversetStructGrid * oversetStructGrid = NULL;

LIB_EXPORT OversetStructGrid * GetOversetStructGrid()
{
    return oversetStructGrid;
}

LIB_EXPORT void CreateOversetStructGrid()
{
    oversetStructGrid = new OversetStructGrid();
}

LIB_EXPORT void DeleteOversetStructGrid()
{
    delete oversetStructGrid; oversetStructGrid = NULL;
}

LIB_EXPORT int GetZoneInverseIndex(int i_view)
{
    return oversetStructGrid->GetZoneInverseIndex(i_view);
}

LIB_EXPORT void ComputePhysicalCoordinate(vector< vector< RDouble > > & frame, RDouble x, RDouble y, RDouble z, vector< RDouble > & r)
{
    for (int i = 0; i < 3; ++ i)
    {
        r[i] = frame[i][0] + x * frame[i][1] + y * frame[i][2] + z * frame[i][3];
    }
    return;
}

}