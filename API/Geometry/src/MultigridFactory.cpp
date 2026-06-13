#include <cmath>
#include "OversetGridFactory.h"
#include "Geo_MultiGridInfo_Struct.h"
#include "Geo_SimpleBC.h"
#include "MultigridFactory.h"
#include "Geo_NodeTopo_Struct.h"
#include "Geo_FaceMetrics_Struct.h"
#include "Geo_CellMetrics_Struct.h"
#include "Geo_DynamicGridMetrics_Struct.h"
#include "Geo_OversetGridTopo_Struct.h"
#include "TK_Exit.h"
#include "Geo_StructBC.h"
#include "Geo_Interface.h"
#include "Constants.h"

using namespace std;

namespace PHSPACE
{
PHSurface::PHSurface(Range I, Range J, Range K)
{
    ist = I.first();
    jst = J.first();
    kst = K.first();

    ni = I.length();
    nj = J.length();
    nk = K.length();

    np = ni * nj * nk;
    data = new int [np];
}

PHSurface::~PHSurface()
{
    delete [] data;
    data = NULL;
}

int & PHSurface::operator()(int i, int j, int k)
{
    return data[ni * nj * (k - kst) + ni * (j - jst) + (i - ist)];
}

int * PHSurface::GetData()
{
    return data;
}

int PHSurface::GetNumberOfFaces()
{
    return np;
}

StructFace::StructFace(StructGrid * structGrid)
{
    this->structGrid                         = structGrid;

    neighborZoneIndexContainer               = new vector< PHSurface * >(6, NULL);
    neighborZoneLocalInterfaceIndexContainer = new vector< PHSurface * >(6, NULL);
    localInterfaceIndexContainer             = new vector< PHSurface * >(6, NULL);
    localBoundaryFaceIndexContainer          = new vector< PHSurface * >(6, NULL);

    numberOfFaceElements.resize(6);

    Init();

    structFaceGroup = NULL;
}

StructFace::~StructFace()
{
    for (int m = 0; m < 6; ++ m)
    {
        delete (* neighborZoneIndexContainer)[m];
        delete (* neighborZoneLocalInterfaceIndexContainer)[m];
        delete (* localInterfaceIndexContainer)[m];
        delete (* localBoundaryFaceIndexContainer)[m];
    }
    delete neighborZoneIndexContainer;               neighborZoneIndexContainer               = NULL;
    delete neighborZoneLocalInterfaceIndexContainer; neighborZoneLocalInterfaceIndexContainer = NULL;
    delete localInterfaceIndexContainer;             localInterfaceIndexContainer             = NULL;
    delete localBoundaryFaceIndexContainer;          localBoundaryFaceIndexContainer          = NULL;
}

void StructFace::PaintSelf()
{
    const int GHOST_INDEX = -1;
    StructBCSet *structBCSet = structGrid->GetStructBCSet();
    
    int iBoundaryFace = 0, iInterface = 0;
    int nBCRegion = structBCSet->GetnBCRegion();
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC * bcRegion = structBCSet->GetBCRegion(iBCRegion);
        
        int fds = bcRegion->GetFaceDirection();
        int fns = bcRegion->GetFaceLeftOrRightIndex();

        int delta = (fns + 1) / 2;
        int mfds = fds * 2 + delta;

        PHSurface &localBoundaryFace = this->GetLocalBoundaryFaceIndexContainer(mfds);
        PHSurface &localInterface    = this->GetLocalInterfaceIndexContainer(mfds);
        vector< int > sst(3), sed(3);
        for (int m = 0; m < 3; ++ m)
        {
            sst[m] = bcRegion->GetStartPoint(m);
            sed[m] = bcRegion->GetEndPoint(m);
        }

        sst[fds] += delta;
        sed[fds] += delta;

        for (int k = sst[2]; k <= sed[2]; ++ k)
        {
            for (int j = sst[1]; j <= sed[1]; ++ j)
            {
                for (int i = sst[0]; i <= sed[0]; ++ i)
                {
                    localBoundaryFace(i, j, k) = iBoundaryFace;
                    localInterface(i, j, k) = GHOST_INDEX;
                    iBoundaryFace += 1;
                }
            }
        }

        int bc_type = bcRegion->GetBCType();

        if (! IsInterface(bc_type)) continue;

        int targetZoneIndex = bcRegion->GetTargetRegionBlock();
        PHSurface & neighborZone = this->GetNeighborZoneIndexContainer(mfds);

        for (int k = sst[2]; k <= sed[2]; ++ k)
        {
            for (int j = sst[1]; j <= sed[1]; ++ j)
            {
                for (int i = sst[0]; i <= sed[0]; ++ i)
                {
                    localInterface(i, j, k) = iInterface;
                    neighborZone(i, j, k) = targetZoneIndex;
                    iInterface += 1;
                }
            }
        }
    }

    return;
}

void StructFace::PaintNeighbor()
{
    vector< int > ti(3), di(3);
    StructBCSet * structBCSet = structGrid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC * bcRegion = structBCSet->GetBCRegion(iBCRegion);
        int bc_type = bcRegion->GetBCType();
        if (! IsInterface(bc_type)) continue;

        int *st_map = bcRegion->GetStam();
        int *rate   = bcRegion->GetRate();

        int fds = bcRegion->GetFaceDirection();
        int fns = bcRegion->GetFaceLeftOrRightIndex();

        int ds = (fns + 1) / 2;
        int ms = fds * 2 + ds;

        vector< int > sst(3), sed(3);
        for (int m = 0; m < 3; ++ m)
        {
            sst[m] = bcRegion->GetStartPoint(m);
            sed[m] = bcRegion->GetEndPoint(m);
        }

        sst[fds] += ds;
        sed[fds] += ds;

        int fdt = bcRegion->GetFaceDirectionOfTargetBlock();
        int fnt = bcRegion->GetFaceLeftOrRightIndexOfTargetBlock();

        int dt = (fnt + 1) / 2;
        int mt = fdt * 2 + dt;

        vector< int > tst(3), ted(3);
        for (int m = 0; m < 3; ++ m)
        {
            tst[m] = bcRegion->GetTargetStart(m);
            ted[m] = bcRegion->GetTargetEnd(m);
        }

        tst[fdt] += dt;
        ted[fdt] += dt;

        int targetZoneIndex = bcRegion->GetTargetRegionBlock();
        StructFace * neighbor = GetStructFaceGroup(targetZoneIndex);
        PHSurface & lid = neighbor->GetLocalInterfaceIndexContainer(mt);
        PHSurface & nid = this->GetNeighborZoneLocalInterfaceIndexContainer(ms);
        
        for (int k = sst[2]; k <= sed[2]; ++ k)
        {
            di[2] = k - sst[2];
            for (int j = sst[1]; j <= sed[1]; ++ j)
            {
                di[1] = j - sst[1];
                for (int i = sst[0]; i <= sed[0]; ++ i)
                {
                    di[0] = i - sst[0];
                    for (int m = 0; m < 3; ++ m)
                    {
                        int n = st_map[m];
                        ti[n] = (tst[n] + ted[n] + rate[m] * (tst[n] - ted[n])) / 2 + rate[m] * di[m];
                    }
                    nid(i, j, k) = lid(ti[0], ti[1], ti[2]);
                }
            }
        }
    }

    return;
}

int StructFace::ComputeNumberOfInterfaces()
{
    const int GHOST_INDEX = -1;

    int numberOfInterfaces = 0;
    for (int m = 0; m < 6; ++ m)
    {
        PHSurface & localInterface = GetLocalInterfaceIndexContainer(m);

        int * data = localInterface.GetData();

        int numberOfElements = localInterface.GetNumberOfFaces();
 
        for (int iElement = 0; iElement < numberOfElements; ++ iElement)
        {
            if (data[iElement] == GHOST_INDEX) continue;

            numberOfInterfaces += 1;
        }
    }

    return numberOfInterfaces;
}

void StructFace::CollectInterfaceInformation(InterfaceInfo * interfaceInformation)
{
    if (!interfaceInformation) return;

    const int GHOST_INDEX = -1;

    int * neighborZoneIndexContainer                 = interfaceInformation->GetInterFace2ZoneID();
    int * neighborZoneLocalInterfaceIndexContainer   = interfaceInformation->GetInterFace2InterFaceID();
    int * localInterfaceToBoundaryFaceIndexesMapping = interfaceInformation->GetInterFace2BoundaryFace();

    for (int m = 0; m < 6; ++ m)
    {
        int * sa = this->GetNeighborZoneIndex(m);
        int * sb = this->GetNeighborZoneLocalInterfaceIndex(m);
        int * sc = this->GetLocalInterfaceIndex(m);
        int * sd = this->GetLocalBoundaryFaceIndex(m);

        int numberOfElements = numberOfFaceElements[m];
        for (int iElement = 0; iElement < numberOfElements; ++ iElement)
        {
            if (sc[iElement] == GHOST_INDEX) continue;
            
            int jZone           = sa[iElement];
            int jLocalInterface = sb[iElement];
            int iLocalInterface = sc[iElement];
            int iBoundaryFace   = sd[iElement];
            
            neighborZoneIndexContainer[iLocalInterface]                 = jZone;
            neighborZoneLocalInterfaceIndexContainer[iLocalInterface]   = jLocalInterface;
            localInterfaceToBoundaryFaceIndexesMapping[iLocalInterface] = iBoundaryFace;
        }
    }

    return;
}

void StructFace::Init()
{
    this->ni = structGrid->GetNI();
    this->nj = structGrid->GetNJ();
    this->nk = structGrid->GetNK();

    numberOfFaceElements[0] = (nj - 1) * (nk - 1); numberOfFaceElements[1] = (nj - 1) * (nk - 1);
    numberOfFaceElements[2] = (ni - 1) * (nk - 1); numberOfFaceElements[3] = (ni - 1) * (nk - 1);
    numberOfFaceElements[4] = (ni - 1) * (nj - 1); numberOfFaceElements[5] = (ni - 1) * (nj - 1);

    IST.setRange(1, 1); ISE.setRange(1, ni - 1); IED.setRange(ni, ni);
    JST.setRange(1, 1); JSE.setRange(1, nj - 1); JED.setRange(nj, nj);
    KST.setRange(1, 1); KSE.setRange(1, nk - 1); KED.setRange(nk, nk);

    GenerateSubspace(neighborZoneIndexContainer);
    GenerateSubspace(neighborZoneLocalInterfaceIndexContainer);
    GenerateSubspace(localInterfaceIndexContainer);
    GenerateSubspace(localBoundaryFaceIndexContainer);

    return;
}

void StructFace::GenerateSubspace(vector< PHSurface * > * micro)
{
    vector< PHSurface * > & smile = * micro;

    smile[0] = new PHSurface(IST, JSE, KSE);
    smile[1] = new PHSurface(IED, JSE, KSE);

    smile[2] = new PHSurface(ISE, JST, KSE);
    smile[3] = new PHSurface(ISE, JED, KSE);

    smile[4] = new PHSurface(ISE, JSE, KST);
    smile[5] = new PHSurface(ISE, JSE, KED);

    return;
}

MultigridManager::MultigridManager(Grid ** OrdinaryGrid, int isOverset, int globalBlockNumber)
{
    this->OrdinaryGrid      = OrdinaryGrid;
    this->isOverset         = isOverset;
    this->globalBlockNumber = globalBlockNumber;
}

MultigridManager::~MultigridManager()
{

}

void MultigridManager::Run()
{
    if (!isOverset) return;

    GenerateLinkStructure();

    GenerateHoleStructure();

    return;
}

void MultigridManager::GenerateGlobalBlockSpace()
{
    startCenterLabel.resize(globalBlockNumber);
    iCenterDimension.resize(globalBlockNumber);
    jCenterDimension.resize(globalBlockNumber);
    kCenterDimension.resize(globalBlockNumber);
}

void MultigridManager::GenerateGlobalCellLinkStructure()
{
    ReadGridCoordinate();
    GridExtend();
}

void MultigridManager::GridExtend()
{
    ComputeCellCenter();
    FaceExtend();
    PrismLineExtend();
    CornerExtend();
}

void MultigridManager::ComputeCellCenter()
{
    const RDouble SMALLTmp = -1.e9, EIGHTH = 0.125;
    int np0, np1, np2, np3, np4, np5, np6, np7;

    for (int iCenter = 0; iCenter < globalCenterNumber; ++ iCenter)
    {
        cellCenterX[ iCenter ] = SMALLTmp;
        cellCenterY[ iCenter ] = SMALLTmp;
        cellCenterZ[ iCenter ] = SMALLTmp;
    }

    for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        for (int k = 1; k < kDimension[ iBlock ]; ++ k)
        {
            for (int j = 1; j < jDimension[ iBlock ]; ++ j)
            {
                for (int i = 1; i < iDimension[ iBlock ]; ++ i)
                {
                    np0 = startPointLabel[ iBlock ]+jDimension[ iBlock ]*iDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;
                    np1 = startPointLabel[ iBlock ]+jDimension[ iBlock ]*iDimension[ iBlock ]*k+iDimension[ iBlock ]*j+(i-1);
                    np2 = startPointLabel[ iBlock ]+jDimension[ iBlock ]*iDimension[ iBlock ]*k+iDimension[ iBlock ]*(j-1)+i;
                    np3 = startPointLabel[ iBlock ]+jDimension[ iBlock ]*iDimension[ iBlock ]*(k-1)+iDimension[ iBlock ]*j+i;
                    np4 = startPointLabel[ iBlock ]+jDimension[ iBlock ]*iDimension[ iBlock ]*k+iDimension[ iBlock ]*(j-1)+(i-1);
                    np5 = startPointLabel[ iBlock ]+jDimension[ iBlock ]*iDimension[ iBlock ]*(k-1)+iDimension[ iBlock ]*j+(i-1);
                    np6 = startPointLabel[ iBlock ]+jDimension[ iBlock ]*iDimension[ iBlock ]*(k-1)+iDimension[ iBlock ]*(j-1)+i;
                    np7 = startPointLabel[ iBlock ]+jDimension[ iBlock ]*iDimension[ iBlock ]*(k-1)+iDimension[ iBlock ]*(j-1)+(i-1);
                    
                    int iCenter = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
                    cellCenterX[ iCenter ] = EIGHTH*(coordinateX[ np0 ]+coordinateX[ np1 ]+coordinateX[ np2 ]+coordinateX[ np3 ]+coordinateX[ np4 ]+coordinateX[ np5 ]+coordinateX[ np6 ]+coordinateX[ np7 ]);
                    cellCenterY[ iCenter ] = EIGHTH*(coordinateY[ np0 ]+coordinateY[ np1 ]+coordinateY[ np2 ]+coordinateY[ np3 ]+coordinateY[ np4 ]+coordinateY[ np5 ]+coordinateY[ np6 ]+coordinateY[ np7 ]);
                    cellCenterZ[ iCenter ] = EIGHTH*(coordinateZ[ np0 ]+coordinateZ[ np1 ]+coordinateZ[ np2 ]+coordinateZ[ np3 ]+coordinateZ[ np4 ]+coordinateZ[ np5 ]+coordinateZ[ np6 ]+coordinateZ[ np7 ]);
                
                    characterPoint[ iCenter ] = np0;
                }
            }
        }
    }
}

void MultigridManager::FaceExtend()
{
    int iBlock, i, j, k, ic, jc, kc, nPointExtend, nPointCorner;
    int nPointPatch, nBlockPatch, iPatch, jPatch, kPatch, nPointTarget;

    for (iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        for (k = 1; k < kDimension[ iBlock ]; ++ k)
        {
            for (j = 1; j < jDimension[ iBlock ]; ++ j)
            {
                i = 1;
                ic = 0;
                nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+ic;
                nPointCorner = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;
                if (patchIPoint[ nPointCorner ] >= 0)
                {
                    nPointPatch = patchIPoint[ nPointCorner ];
                    nBlockPatch = blockLocation[ nPointPatch ];
                    iPatch = iLocation[ nPointPatch ];
                    jPatch = jLocation[ nPointPatch ];
                    kPatch = kLocation[ nPointPatch ];
                    nPointTarget = startCenterLabel[ nBlockPatch ]+iCenterDimension[ nBlockPatch ]*jCenterDimension[ nBlockPatch ]*kPatch+iCenterDimension[ nBlockPatch ]*jPatch+iPatch;
                    characterPoint[ nPointExtend ] = patchIPoint[ nPointCorner ];
                }
                else
                {
                    nPointTarget = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
                    characterPoint[ nPointExtend ] = nPointCorner;
                }
                cellCenterX[ nPointExtend ] = cellCenterX[ nPointTarget ];
                cellCenterY[ nPointExtend ] = cellCenterY[ nPointTarget ];
                cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointTarget ];

                i  = iCenterDimension[ iBlock ]-2;
                ic = iCenterDimension[ iBlock ]-1;
                
                nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+ic;
                nPointCorner = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;
                if (patchIPoint[ nPointCorner ] >= 0)
                {
                    nPointPatch = patchIPoint[ nPointCorner ];
                    nBlockPatch = blockLocation[ nPointPatch ];
                    iPatch = iLocation[ nPointPatch ];
                    jPatch = jLocation[ nPointPatch ];
                    kPatch = kLocation[ nPointPatch ];
                    nPointTarget = startCenterLabel[ nBlockPatch ]+iCenterDimension[ nBlockPatch ]*jCenterDimension[ nBlockPatch ]*kPatch+iCenterDimension[ nBlockPatch ]*jPatch+iPatch;
                    characterPoint[ nPointExtend ] = patchIPoint[ nPointCorner ];
                }
                else
                {
                    nPointTarget = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
                    characterPoint[ nPointExtend ] = nPointCorner;
                }
                cellCenterX[ nPointExtend ] = cellCenterX[ nPointTarget ];
                cellCenterY[ nPointExtend ] = cellCenterY[ nPointTarget ];
                cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointTarget ];
            }
        }

        for (i = 1; i < iDimension[ iBlock ]; ++ i)
        {
            for (k = 1; k < kDimension[ iBlock ]; ++ k)
            {
                j = 1;
                jc = 0;

                nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*jc+i;
                nPointCorner = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;

                if (patchJPoint[ nPointCorner ] >= 0)
                {
                    nPointPatch = patchJPoint[ nPointCorner ];
                    nBlockPatch = blockLocation[ nPointPatch ];
                    iPatch = iLocation[ nPointPatch ];
                    jPatch = jLocation[ nPointPatch ];
                    kPatch = kLocation[ nPointPatch ];
                    nPointTarget = startCenterLabel[ nBlockPatch ]+iCenterDimension[ nBlockPatch ]*jCenterDimension[ nBlockPatch ]*kPatch+iCenterDimension[ nBlockPatch ]*jPatch+iPatch;
                    characterPoint[ nPointExtend ] = patchJPoint[ nPointCorner ];
                }
                else
                {
                    nPointTarget = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
                    characterPoint[ nPointExtend ] = nPointCorner;
                }
                cellCenterX[ nPointExtend ] = cellCenterX[ nPointTarget ];
                cellCenterY[ nPointExtend ] = cellCenterY[ nPointTarget ];
                cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointTarget ];

                j  = jCenterDimension[ iBlock ]-2;
                jc = jCenterDimension[ iBlock ]-1;

                nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*jc+i;
                nPointCorner = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;

                if (patchJPoint[ nPointCorner ] >= 0)
                {
                    nPointPatch = patchJPoint[ nPointCorner ];
                    nBlockPatch = blockLocation[ nPointPatch ];
                    iPatch = iLocation[ nPointPatch ];
                    jPatch = jLocation[ nPointPatch ];
                    kPatch = kLocation[ nPointPatch ];
                    nPointTarget = startCenterLabel[ nBlockPatch ]+iCenterDimension[ nBlockPatch ]*jCenterDimension[ nBlockPatch ]*kPatch+iCenterDimension[ nBlockPatch ]*jPatch+iPatch;
                    characterPoint[ nPointExtend ] = patchJPoint[ nPointCorner ];
                }
                else
                {
                    nPointTarget = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
                    characterPoint[ nPointExtend ] = nPointCorner;
                }
                cellCenterX[ nPointExtend ] = cellCenterX[ nPointTarget ];
                cellCenterY[ nPointExtend ] = cellCenterY[ nPointTarget ];
                cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointTarget ];
            }
        }

        for (j = 1; j < jDimension[ iBlock ]; ++ j)
        {
            for (i = 1; i < iDimension[ iBlock ]; ++ i)
            {
                k = 1;
                kc = 0;

                nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*kc+iCenterDimension[ iBlock ]*j+i;
                nPointCorner = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;

                if (patchKPoint[ nPointCorner ] >= 0)
                {
                    nPointPatch = patchKPoint[ nPointCorner ];
                    nBlockPatch = blockLocation[ nPointPatch ];
                    iPatch = iLocation[ nPointPatch ];
                    jPatch = jLocation[ nPointPatch ];
                    kPatch = kLocation[ nPointPatch ];
                    nPointTarget = startCenterLabel[ nBlockPatch ]+iCenterDimension[ nBlockPatch ]*jCenterDimension[ nBlockPatch ]*kPatch+iCenterDimension[ nBlockPatch ]*jPatch+iPatch;
                    characterPoint[ nPointExtend ] = patchKPoint[ nPointCorner ];
                }
                else
                {
                    nPointTarget = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
                    characterPoint[ nPointExtend ] = nPointCorner;
                }
                cellCenterX[ nPointExtend ] = cellCenterX[ nPointTarget ];
                cellCenterY[ nPointExtend ] = cellCenterY[ nPointTarget ];
                cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointTarget ];

                k  = kCenterDimension[ iBlock ]-2;
                kc = kCenterDimension[ iBlock ]-1;

                nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*kc+iCenterDimension[ iBlock ]*j+i;
                nPointCorner = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;

                if (patchKPoint[ nPointCorner ] >= 0)
                {
                    nPointPatch = patchKPoint[ nPointCorner ];
                    nBlockPatch = blockLocation[ nPointPatch ];
                    iPatch = iLocation[ nPointPatch ];
                    jPatch = jLocation[ nPointPatch ];
                    kPatch = kLocation[ nPointPatch ];
                    nPointTarget = startCenterLabel[ nBlockPatch ]+iCenterDimension[ nBlockPatch ]*jCenterDimension[ nBlockPatch ]*kPatch+iCenterDimension[ nBlockPatch ]*jPatch+iPatch;
                    characterPoint[ nPointExtend ] = patchKPoint[ nPointCorner ];
                }
                else
                {
                    nPointTarget = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
                    characterPoint[ nPointExtend ] = nPointCorner;
                }
                cellCenterX[ nPointExtend ] = cellCenterX[ nPointTarget ];
                cellCenterY[ nPointExtend ] = cellCenterY[ nPointTarget ];
                cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointTarget ];
            }
        }
    }
}

void MultigridManager::PrismLineExtend()
{
    RDouble xx,yy,zz;
    int iBlock, i, j, k, nPointer, iC, jC, kC, nPointExtend, nPointCorner;
    int nFlap, nPointBy, nBlock, iPoint, jPoint, kPoint, nPointG;
    int iCounter = 0;

    for (iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        //! 1.
        for (k = 1; k < kDimension[ iBlock ]; ++ k)
        {
            j = 1;
            i = 1;
            nPointer = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;
            iC = 0;
            jC = 0;
            kC = k;
            nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*kC+iCenterDimension[ iBlock ]*jC+iC;
            nPointCorner = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
            cellCenterX[ nPointExtend ] = cellCenterX[ nPointCorner ];
            cellCenterY[ nPointExtend ] = cellCenterY[ nPointCorner ];
            cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointCorner ];
            characterPoint[ nPointExtend ] = characterPoint[ nPointCorner ];
            if (patchIPoint[ nPointer ] >= 0 && patchJPoint[ nPointer ] >= 0)
            {
                nFlap = 0;
                GetCrossPoint(patchIPoint[ nPointer ], patchJPoint[ nPointer ], cellCenterX[ nPointCorner ], cellCenterY[ nPointCorner ], cellCenterZ[ nPointCorner ], xx, yy, zz, nFlap, nPointBy, nBlock, iPoint, jPoint, kPoint);
                if (nFlap == 1)
                {
                    cellCenterX[ nPointExtend ] = xx;
                    cellCenterY[ nPointExtend ] = yy;
                    cellCenterZ[ nPointExtend ] = zz;
                    characterPoint[ nPointExtend ] = characterPoint[ nPointBy ];
                }
                else
                {
                    nPointG = startCenterLabel[ nBlock ]+iCenterDimension[ nBlock ]*jCenterDimension[ nBlock ]*kPoint+iCenterDimension[ nBlock ]*jPoint+iPoint;
                    cellCenterX[ nPointExtend ] = cellCenterX[ nPointG ];
                    cellCenterY[ nPointExtend ] = cellCenterY[ nPointG ];
                    cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointG ];
                    characterPoint[ nPointExtend ] = characterPoint[ nPointG ];
                    iCounter += 1;
                }
            }
        }

        //! 2.
        for (k = 1; k < kDimension[ iBlock ]; ++ k)
        {
            i = 1;
            j = jDimension[ iBlock ]-1;
            nPointer = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;
            iC = 0;
            jC = jCenterDimension[ iBlock ]-1;
            kC = k;
            nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*kC+iCenterDimension[ iBlock ]*jC+iC;
            nPointCorner = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
            cellCenterX[ nPointExtend ] = cellCenterX[ nPointCorner ];
            cellCenterY[ nPointExtend ] = cellCenterY[ nPointCorner ];
            cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointCorner ];
            characterPoint[ nPointExtend ] = characterPoint[ nPointCorner ];
            if (patchIPoint[ nPointer ] >= 0 && patchJPoint[ nPointer ] >= 0)
            {
                nFlap = 0;
                GetCrossPoint(patchIPoint[ nPointer ], patchJPoint[ nPointer ], cellCenterX[ nPointCorner ], cellCenterY[ nPointCorner ], cellCenterZ[ nPointCorner ], xx, yy, zz, nFlap, nPointBy, nBlock, iPoint, jPoint, kPoint);
                if (nFlap == 1)
                {
                    cellCenterX[ nPointExtend ] = xx;
                    cellCenterY[ nPointExtend ] = yy;
                    cellCenterZ[ nPointExtend ] = zz;
                    characterPoint[ nPointExtend ] = characterPoint[ nPointBy ];
                }
                else
                {
                    nPointG = startCenterLabel[ nBlock ]+iCenterDimension[ nBlock ]*jCenterDimension[ nBlock ]*kPoint+iCenterDimension[ nBlock ]*jPoint+iPoint;
                    cellCenterX[ nPointExtend ] = cellCenterX[ nPointG ];
                    cellCenterY[ nPointExtend ] = cellCenterY[ nPointG ];
                    cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointG ];
                    characterPoint[ nPointExtend ] = characterPoint[ nPointG ];
                    iCounter += 1;
                }
            }
        }

        //! 3.
        for (k = 1; k < kDimension[ iBlock ]; ++ k)
        {
            i = iDimension[ iBlock ]-1;
            j = 1;
            nPointer = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;
            iC = iCenterDimension[ iBlock ]-1;
            jC = 0;
            kC = k;
            nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*kC+iCenterDimension[ iBlock ]*jC+iC;
            nPointCorner = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
            cellCenterX[ nPointExtend ] = cellCenterX[ nPointCorner ];
            cellCenterY[ nPointExtend ] = cellCenterY[ nPointCorner ];
            cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointCorner ];
            characterPoint[ nPointExtend ] = characterPoint[ nPointCorner ];
            if (patchIPoint[ nPointer ] >= 0 && patchJPoint[ nPointer ] >= 0)
            {
                nFlap = 0;
                GetCrossPoint(patchIPoint[ nPointer ], patchJPoint[ nPointer ], cellCenterX[ nPointCorner ], cellCenterY[ nPointCorner ], cellCenterZ[ nPointCorner ], xx, yy, zz, nFlap, nPointBy, nBlock, iPoint, jPoint, kPoint);
                if (nFlap == 1)
                {
                    cellCenterX[ nPointExtend ] = xx;
                    cellCenterY[ nPointExtend ] = yy;
                    cellCenterZ[ nPointExtend ] = zz;
                    characterPoint[ nPointExtend ] = characterPoint[ nPointBy ];
                }
                else
                {
                    nPointG = startCenterLabel[ nBlock ]+iCenterDimension[ nBlock ]*jCenterDimension[ nBlock ]*kPoint+iCenterDimension[ nBlock ]*jPoint+iPoint;
                    cellCenterX[ nPointExtend ] = cellCenterX[ nPointG ];
                    cellCenterY[ nPointExtend ] = cellCenterY[ nPointG ];
                    cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointG ];
                    characterPoint[ nPointExtend ] = characterPoint[ nPointG ];
                    iCounter += 1;
                }
            }
        }

        //! 4.
        for (k = 1; k < kDimension[ iBlock ]; ++ k)
        {
            i = iDimension[ iBlock ]-1;
            j = jDimension[ iBlock ]-1;
            nPointer = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;
            iC = iCenterDimension[ iBlock ]-1;
            jC = jCenterDimension[ iBlock ]-1;
            kC = k;
            nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*kC+iCenterDimension[ iBlock ]*jC+iC;
            nPointCorner = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
            cellCenterX[ nPointExtend ] = cellCenterX[ nPointCorner ];
            cellCenterY[ nPointExtend ] = cellCenterY[ nPointCorner ];
            cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointCorner ];
            characterPoint[ nPointExtend ] = characterPoint[ nPointCorner ];
            if (patchIPoint[ nPointer ] >= 0 && patchJPoint[ nPointer ] >= 0)
            {
                nFlap = 0;
                GetCrossPoint(patchIPoint[ nPointer ], patchJPoint[ nPointer ], cellCenterX[ nPointCorner ], cellCenterY[ nPointCorner ], cellCenterZ[ nPointCorner ], xx, yy, zz, nFlap, nPointBy, nBlock, iPoint, jPoint, kPoint);
                if (nFlap == 1)
                {
                    cellCenterX[ nPointExtend ] = xx;
                    cellCenterY[ nPointExtend ] = yy;
                    cellCenterZ[ nPointExtend ] = zz;
                    characterPoint[ nPointExtend ] = characterPoint[ nPointBy ];
                }
                else
                {
                    nPointG = startCenterLabel[ nBlock ]+iCenterDimension[ nBlock ]*jCenterDimension[ nBlock ]*kPoint+iCenterDimension[ nBlock ]*jPoint+iPoint;
                    cellCenterX[ nPointExtend ] = cellCenterX[ nPointG ];
                    cellCenterY[ nPointExtend ] = cellCenterY[ nPointG ];
                    cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointG ];
                    characterPoint[ nPointExtend ] = characterPoint[ nPointG ];
                    iCounter += 1;
                }
            }
        }

        //! 5.
        for (j = 1; j < jDimension[ iBlock ]; ++ j)
        {
            i = 1;
            k = 1;
            nPointer = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;
            iC = 0;
            jC = j;
            kC = 0;
            nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*kC+iCenterDimension[ iBlock ]*jC+iC;
            nPointCorner = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
            cellCenterX[ nPointExtend ] = cellCenterX[ nPointCorner ];
            cellCenterY[ nPointExtend ] = cellCenterY[ nPointCorner ];
            cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointCorner ];
            characterPoint[ nPointExtend ] = characterPoint[ nPointCorner ];
            if (patchIPoint[ nPointer ] >= 0 && patchKPoint[ nPointer ] >= 0)
            {
                nFlap = 0;
                GetCrossPoint(patchIPoint[ nPointer ], patchKPoint[ nPointer ], cellCenterX[ nPointCorner ], cellCenterY[ nPointCorner ], cellCenterZ[ nPointCorner ], xx, yy, zz, nFlap, nPointBy, nBlock, iPoint, jPoint, kPoint);
                if (nFlap == 1)
                {
                    cellCenterX[ nPointExtend ] = xx;
                    cellCenterY[ nPointExtend ] = yy;
                    cellCenterZ[ nPointExtend ] = zz;
                    characterPoint[ nPointExtend ] = characterPoint[ nPointBy ];
                }
                else
                {
                    nPointG = startCenterLabel[ nBlock ]+iCenterDimension[ nBlock ]*jCenterDimension[ nBlock ]*kPoint+iCenterDimension[ nBlock ]*jPoint+iPoint;
                    cellCenterX[ nPointExtend ] = cellCenterX[ nPointG ];
                    cellCenterY[ nPointExtend ] = cellCenterY[ nPointG ];
                    cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointG ];
                    characterPoint[ nPointExtend ] = characterPoint[ nPointG ];
                    iCounter += 1;
                }
            }
        }

        //! 6.
        for (j = 1; j < jDimension[ iBlock ]; ++ j)
        {
            i = 1;
            k = kDimension[ iBlock ]-1;
            nPointer = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;
            iC = 0;
            jC = j;
            kC = kCenterDimension[ iBlock ]-1;
            nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*kC+iCenterDimension[ iBlock ]*jC+iC;
            nPointCorner = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
            cellCenterX[ nPointExtend ] = cellCenterX[ nPointCorner ];
            cellCenterY[ nPointExtend ] = cellCenterY[ nPointCorner ];
            cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointCorner ];
            characterPoint[ nPointExtend ] = characterPoint[ nPointCorner ];
            if (patchIPoint[ nPointer ] >= 0 && patchKPoint[ nPointer ] >= 0)
            {
                nFlap = 0;
                GetCrossPoint(patchIPoint[ nPointer ], patchKPoint[ nPointer ], cellCenterX[ nPointCorner ], cellCenterY[ nPointCorner ], cellCenterZ[ nPointCorner ], xx, yy, zz, nFlap, nPointBy, nBlock, iPoint, jPoint, kPoint);
                if (nFlap == 1)
                {
                    cellCenterX[ nPointExtend ] = xx;
                    cellCenterY[ nPointExtend ] = yy;
                    cellCenterZ[ nPointExtend ] = zz;
                    characterPoint[ nPointExtend ] = characterPoint[ nPointBy ];
                }
                else
                {
                    nPointG = startCenterLabel[ nBlock ]+iCenterDimension[ nBlock ]*jCenterDimension[ nBlock ]*kPoint+iCenterDimension[ nBlock ]*jPoint+iPoint;
                    cellCenterX[ nPointExtend ] = cellCenterX[ nPointG ];
                    cellCenterY[ nPointExtend ] = cellCenterY[ nPointG ];
                    cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointG ];
                    characterPoint[ nPointExtend ] = characterPoint[ nPointG ];
                    iCounter += 1;
                }
            }
        }

        //! 7.
        for (j = 1; j < jDimension[ iBlock ]; ++ j)
        {
            k = 1;
            i = iDimension[ iBlock ]-1;
            nPointer = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;
            iC = iCenterDimension[ iBlock ]-1;
            jC = j;
            kC = 0;
            nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*kC+iCenterDimension[ iBlock ]*jC+iC;
            nPointCorner = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
            cellCenterX[ nPointExtend ] = cellCenterX[ nPointCorner ];
            cellCenterY[ nPointExtend ] = cellCenterY[ nPointCorner ];
            cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointCorner ];
            characterPoint[ nPointExtend ] = characterPoint[ nPointCorner ];
            if (patchIPoint[ nPointer ] >= 0 && patchKPoint[ nPointer ] >= 0)
            {
                nFlap = 0;
                GetCrossPoint(patchIPoint[ nPointer ], patchKPoint[ nPointer ], cellCenterX[ nPointCorner ], cellCenterY[ nPointCorner ], cellCenterZ[ nPointCorner ], xx, yy, zz, nFlap, nPointBy, nBlock, iPoint, jPoint, kPoint);
                if (nFlap == 1)
                {
                    cellCenterX[ nPointExtend ] = xx;
                    cellCenterY[ nPointExtend ] = yy;
                    cellCenterZ[ nPointExtend ] = zz;
                    characterPoint[ nPointExtend ] = characterPoint[ nPointBy ];
                }
                else
                {
                    nPointG = startCenterLabel[ nBlock ]+iCenterDimension[ nBlock ]*jCenterDimension[ nBlock ]*kPoint+iCenterDimension[ nBlock ]*jPoint+iPoint;
                    cellCenterX[ nPointExtend ] = cellCenterX[ nPointG ];
                    cellCenterY[ nPointExtend ] = cellCenterY[ nPointG ];
                    cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointG ];
                    characterPoint[ nPointExtend ] = characterPoint[ nPointG ];
                    iCounter += 1;
                }
            }
        }
        //! 8.
        for (j = 1; j < jDimension[ iBlock ]; ++ j)
        {
            k = kDimension[ iBlock ]-1;
            i = iDimension[ iBlock ]-1;
            nPointer = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;
            iC = iCenterDimension[ iBlock ]-1;
            jC = j;
            kC = kCenterDimension[ iBlock ]-1;
            nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*kC+iCenterDimension[ iBlock ]*jC+iC;
            nPointCorner = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
            cellCenterX[ nPointExtend ] = cellCenterX[ nPointCorner ];
            cellCenterY[ nPointExtend ] = cellCenterY[ nPointCorner ];
            cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointCorner ];
            characterPoint[ nPointExtend ] = characterPoint[ nPointCorner ];
            if (patchIPoint[ nPointer ] >= 0 && patchKPoint[ nPointer ] >= 0)
            {
                nFlap = 0;
                GetCrossPoint(patchIPoint[ nPointer ], patchKPoint[ nPointer ], cellCenterX[ nPointCorner ], cellCenterY[ nPointCorner ], cellCenterZ[ nPointCorner ], xx, yy, zz, nFlap, nPointBy, nBlock, iPoint, jPoint, kPoint);
                if (nFlap == 1)
                {
                    cellCenterX[ nPointExtend ] = xx;
                    cellCenterY[ nPointExtend ] = yy;
                    cellCenterZ[ nPointExtend ] = zz;
                    characterPoint[ nPointExtend ] = characterPoint[ nPointBy ];
                }
                else
                {
                    nPointG = startCenterLabel[ nBlock ]+iCenterDimension[ nBlock ]*jCenterDimension[ nBlock ]*kPoint+iCenterDimension[ nBlock ]*jPoint+iPoint;
                    cellCenterX[ nPointExtend ] = cellCenterX[ nPointG ];
                    cellCenterY[ nPointExtend ] = cellCenterY[ nPointG ];
                    cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointG ];
                    characterPoint[ nPointExtend ] = characterPoint[ nPointG ];
                    iCounter += 1;
                }
            }
        }

        //! 9.
        for (i = 1; i < iDimension[ iBlock ]; ++ i)
        {
            k = kDimension[ iBlock ]-1;
            j = jDimension[ iBlock ]-1;
            nPointer = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;
            iC = i;
            jC = jCenterDimension[ iBlock ]-1;
            kC = kCenterDimension[ iBlock ]-1;
            nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*kC+iCenterDimension[ iBlock ]*jC+iC;
            nPointCorner = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
            cellCenterX[ nPointExtend ] = cellCenterX[ nPointCorner ];
            cellCenterY[ nPointExtend ] = cellCenterY[ nPointCorner ];
            cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointCorner ];
            characterPoint[ nPointExtend ] = characterPoint[ nPointCorner ];
            if (patchJPoint[ nPointer ] >= 0 && patchKPoint[ nPointer ] >= 0)
            {
                nFlap = 0;
                GetCrossPoint(patchJPoint[ nPointer ], patchKPoint[ nPointer ], cellCenterX[ nPointCorner ], cellCenterY[ nPointCorner ], cellCenterZ[ nPointCorner ], xx, yy, zz, nFlap, nPointBy, nBlock, iPoint, jPoint, kPoint);
                if (nFlap == 1)
                {
                    cellCenterX[ nPointExtend ] = xx;
                    cellCenterY[ nPointExtend ] = yy;
                    cellCenterZ[ nPointExtend ] = zz;
                    characterPoint[ nPointExtend ] = characterPoint[ nPointBy ];
                }
                else
                {
                    nPointG = startCenterLabel[ nBlock ]+iCenterDimension[ nBlock ]*jCenterDimension[ nBlock ]*kPoint+iCenterDimension[ nBlock ]*jPoint+iPoint;
                    cellCenterX[ nPointExtend ] = cellCenterX[ nPointG ];
                    cellCenterY[ nPointExtend ] = cellCenterY[ nPointG ];
                    cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointG ];
                    characterPoint[ nPointExtend ] = characterPoint[ nPointG ];
                    iCounter += 1;
                }
            }
        }

        //! 10.
        for (i = 1; i < iDimension[ iBlock ]; ++ i)
        {
            k = 1;
            j = jDimension[ iBlock ]-1;
            nPointer = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;
            iC = i;
            jC = jCenterDimension[ iBlock ]-1;
            kC = 0;
            nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*kC+iCenterDimension[ iBlock ]*jC+iC;
            nPointCorner = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
            cellCenterX[ nPointExtend ] = cellCenterX[ nPointCorner ];
            cellCenterY[ nPointExtend ] = cellCenterY[ nPointCorner ];
            cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointCorner ];
            characterPoint[ nPointExtend ] = characterPoint[ nPointCorner ];
            if (patchJPoint[ nPointer ] >= 0 && patchKPoint[ nPointer ] >= 0)
            {
                nFlap = 0;
                GetCrossPoint(patchJPoint[ nPointer ], patchKPoint[ nPointer ], cellCenterX[ nPointCorner ], cellCenterY[ nPointCorner ], cellCenterZ[ nPointCorner ], xx, yy, zz, nFlap, nPointBy, nBlock, iPoint, jPoint, kPoint);
                if (nFlap == 1)
                {
                    cellCenterX[ nPointExtend ] = xx;
                    cellCenterY[ nPointExtend ] = yy;
                    cellCenterZ[ nPointExtend ] = zz;
                    characterPoint[ nPointExtend ] = characterPoint[ nPointBy ];
                }
                else
                {
                    nPointG = startCenterLabel[ nBlock ]+iCenterDimension[ nBlock ]*jCenterDimension[ nBlock ]*kPoint+iCenterDimension[ nBlock ]*jPoint+iPoint;
                    cellCenterX[ nPointExtend ] = cellCenterX[ nPointG ];
                    cellCenterY[ nPointExtend ] = cellCenterY[ nPointG ];
                    cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointG ];
                    characterPoint[ nPointExtend ] = characterPoint[ nPointG ];
                    iCounter += 1;
                }
            }
        }

        //! 11.
        for (i = 1; i < iDimension[ iBlock ]; ++ i)
        {
            j = 1;
            k = kDimension[ iBlock ]-1;
            nPointer = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;
            iC = i;
            jC = 0;
            kC = kCenterDimension[ iBlock ]-1;
            nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*kC+iCenterDimension[ iBlock ]*jC+iC;
            nPointCorner = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
            cellCenterX[ nPointExtend ] = cellCenterX[ nPointCorner ];
            cellCenterY[ nPointExtend ] = cellCenterY[ nPointCorner ];
            cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointCorner ];
            characterPoint[ nPointExtend ] = characterPoint[ nPointCorner ];
            if (patchJPoint[ nPointer ] >= 0 && patchKPoint[ nPointer ] >= 0)
            {
                nFlap = 0;
                GetCrossPoint(patchJPoint[ nPointer ], patchKPoint[ nPointer ], cellCenterX[ nPointCorner ], cellCenterY[ nPointCorner ], cellCenterZ[ nPointCorner ], xx, yy, zz, nFlap, nPointBy, nBlock, iPoint, jPoint, kPoint);
                if (nFlap == 1)
                {
                    cellCenterX[ nPointExtend ] = xx;
                    cellCenterY[ nPointExtend ] = yy;
                    cellCenterZ[ nPointExtend ] = zz;
                    characterPoint[ nPointExtend ] = characterPoint[ nPointBy ];
                }
                else
                {
                    nPointG = startCenterLabel[ nBlock ]+iCenterDimension[ nBlock ]*jCenterDimension[ nBlock ]*kPoint+iCenterDimension[ nBlock ]*jPoint+iPoint;
                    cellCenterX[ nPointExtend ] = cellCenterX[ nPointG ];
                    cellCenterY[ nPointExtend ] = cellCenterY[ nPointG ];
                    cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointG ];
                    characterPoint[ nPointExtend ] = characterPoint[ nPointG ];
                    iCounter += 1;
                }
            }
        }

        //! 12.
        for (i = 1; i < iDimension[ iBlock ]; ++ i)
        {
            j = 1;
            k = 1;
            nPointer = startPointLabel[ iBlock ]+iDimension[ iBlock ]*jDimension[ iBlock ]*k+iDimension[ iBlock ]*j+i;
            iC = i;
            jC = 0;
            kC = 0;
            nPointExtend = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*kC+iCenterDimension[ iBlock ]*jC+iC;
            nPointCorner = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
            cellCenterX[ nPointExtend ] = cellCenterX[ nPointCorner ];
            cellCenterY[ nPointExtend ] = cellCenterY[ nPointCorner ];
            cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointCorner ];
            characterPoint[ nPointExtend ] = characterPoint[ nPointCorner ];
            if (patchJPoint[ nPointer ] >= 0 && patchKPoint[ nPointer ] >= 0)
            {
                nFlap = 0;
                GetCrossPoint(patchJPoint[ nPointer ], patchKPoint[ nPointer ], cellCenterX[ nPointCorner ], cellCenterY[ nPointCorner ], cellCenterZ[ nPointCorner ], xx, yy, zz, nFlap, nPointBy, nBlock, iPoint, jPoint, kPoint);
                if (nFlap == 1)
                {
                    cellCenterX[ nPointExtend ] = xx;
                    cellCenterY[ nPointExtend ] = yy;
                    cellCenterZ[ nPointExtend ] = zz;
                    characterPoint[ nPointExtend ] = characterPoint[ nPointBy ];
                }
                else
                {
                    nPointG = startCenterLabel[ nBlock ]+iCenterDimension[ nBlock ]*jCenterDimension[ nBlock ]*kPoint+iCenterDimension[ nBlock ]*jPoint+iPoint;
                    cellCenterX[ nPointExtend ] = cellCenterX[ nPointG ];
                    cellCenterY[ nPointExtend ] = cellCenterY[ nPointG ];
                    cellCenterZ[ nPointExtend ] = cellCenterZ[ nPointG ];
                    characterPoint[ nPointExtend ] = characterPoint[ nPointG ];
                    iCounter += 1;
                }
            }
        }
    }
}

void MultigridManager::GetCrossPoint(int mPatch, int nPatch, RDouble x, RDouble y, RDouble z, RDouble & xx, RDouble & yy, RDouble & zz, int & nFlap, int & nPointBy, int & nBlock, int & iPoint, int & jPoint, int & kPoint)
{
    const RDouble EPS = 1.e-4;
    int iBym[ 6 ], jBym[ 6 ], kBym[ 6 ], nPointBym[ 6 ], iByn[ 6 ], jByn[ 6 ], kByn[ 6 ], nPointByn[ 6 ];
    RDouble xBym[ 6 ], yBym[ 6 ], zBym[ 6 ], xByn[ 6 ], yByn[ 6 ], zByn[ 6 ];
    RDouble dx, dy, dz, p, q, r, distance;
    int mBlock, ipm, jpm, kpm;
    mBlock = blockLocation[ mPatch ];
    ipm = iLocation[ mPatch ];
    jpm = jLocation[ mPatch ];
    kpm = kLocation[ mPatch ];
    nBlock = blockLocation[ nPatch ];
    iPoint = iLocation[ nPatch ];
    jPoint = jLocation[ nPatch ];
    kPoint = kLocation[ nPatch ];

    iBym[ 0 ] = ipm-1; iBym[ 1 ] = ipm+1;  iBym[ 2 ] = ipm;
    jBym[ 0 ] = jpm;   jBym[ 1 ] = jpm;    jBym[ 2 ] = jpm-1;
    kBym[ 0 ] = kpm;   kBym[ 1 ] = kpm;    kBym[ 2 ] = kpm;

    iBym[ 3 ] = ipm;   iBym[ 4 ] = ipm;    iBym[ 5 ] = ipm;
    jBym[ 3 ] = jpm+1; jBym[ 4 ] = jpm;    jBym[ 5 ] = jpm;
    kBym[ 3 ] = kpm;   kBym[ 4 ] = kpm-1;  kBym[ 5 ] = kpm+1;

    iByn[ 0 ] = iPoint-1;iByn[ 1 ] = iPoint+1;iByn[ 2 ] = iPoint;
    jByn[ 0 ] = jPoint;jByn[ 1 ] = jPoint;jByn[ 2 ] = jPoint-1;
    kByn[ 0 ] = kPoint;kByn[ 1 ] = kPoint;kByn[ 2 ] = kPoint;

    iByn[ 3 ] = iPoint;iByn[ 4 ] = iPoint;iByn[ 5 ] = iPoint;
    jByn[ 3 ] = jPoint+1;jByn[ 4 ] = jPoint;jByn[ 5 ] = jPoint;
    kByn[ 3 ] = kPoint;kByn[ 4 ] = kPoint-1;kByn[ 5 ] = kPoint+1;

    for (int i = 0; i < 6; ++ i)
    {
        nPointBym[ i ] = startCenterLabel[ mBlock ]+iCenterDimension[ mBlock ]*jCenterDimension[ mBlock ]*kBym[ i ]+iCenterDimension[ mBlock ]*jBym[ i ]+iBym[ i ];
        xBym[ i ] = cellCenterX[ nPointBym[ i ] ];
        yBym[ i ] = cellCenterY[ nPointBym[ i ] ];
        zBym[ i ] = cellCenterZ[ nPointBym[ i ] ];

        nPointByn[ i ] = startCenterLabel[ nBlock ]+iCenterDimension[ nBlock ]*jCenterDimension[ nBlock ]*kByn[ i ]+iCenterDimension[ nBlock ]*jByn[ i ]+iByn[ i ];
        xByn[ i ] = cellCenterX[ nPointByn[ i ] ];
        yByn[ i ] = cellCenterY[ nPointByn[ i ] ];
        zByn[ i ] = cellCenterZ[ nPointByn[ i ] ];
    }

    for (int m = 0; m < 6; ++ m)
    {
        for (int n = 0; n < 6; ++ n)
        {
            dx = abs(xBym[ m ]-xByn[ n ]);
            dy = abs(yBym[ m ]-yByn[ n ]);
            dz = abs(zBym[ m ]-zByn[ n ]);
            if (dx < EPS && dy < EPS && dz < EPS)
            {
                p = x-xByn[ n ];
                q = y-yByn[ n ];
                r = z-zByn[ n ];
                distance = sqrt(p*p+q*q+r*r);
                if (distance > EPS)
                {
                    xx = xByn[ n ];
                    yy = yByn[ n ];
                    zz = zByn[ n ];
                    nPointBy = nPointByn[ n ];
                    nFlap = 1;
                    return;
                }
            }
        }
    }

    nFlap = 0;
    return;
}

void MultigridManager::CornerExtend()
{
    int iBlock, i, j, k, mPoint, nPoint;
    for (iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        //! p1.
        i = 0;
        j = 0;
        k = 0;
        mPoint = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*(k+1)+iCenterDimension[ iBlock ]*j+i;
        nPoint = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
        cellCenterX[ nPoint ] = cellCenterX[ mPoint ];
        cellCenterY[ nPoint ] = cellCenterY[ mPoint ];
        cellCenterZ[ nPoint ] = cellCenterZ[ mPoint ];
        characterPoint[ nPoint ] = characterPoint[ mPoint ];

        //! p2.
        i = iCenterDimension[ iBlock ]-1;
        j = 0;
        k = 0;
        mPoint = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*(k+1)+iCenterDimension[ iBlock ]*j+i;
        nPoint = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
        cellCenterX[ nPoint ] = cellCenterX[ mPoint ];
        cellCenterY[ nPoint ] = cellCenterY[ mPoint ];
        cellCenterZ[ nPoint ] = cellCenterZ[ mPoint ];
        characterPoint[ nPoint ] = characterPoint[ mPoint ];

        //! p3.
        i = 0;
        j = jCenterDimension[ iBlock ]-1;
        k = 0;
        mPoint = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*(k+1)+iCenterDimension[ iBlock ]*j+i;
        nPoint = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
        cellCenterX[ nPoint ] = cellCenterX[ mPoint ];
        cellCenterY[ nPoint ] = cellCenterY[ mPoint ];
        cellCenterZ[ nPoint ] = cellCenterZ[ mPoint ];
        characterPoint[ nPoint ] = characterPoint[ mPoint ];

        //! p4.
        i = 0;
        j = 0;
        k = kCenterDimension[ iBlock ]-1;
        mPoint = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*(k-1)+iCenterDimension[ iBlock ]*j+i;
        nPoint = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
        cellCenterX[ nPoint ] = cellCenterX[ mPoint ];
        cellCenterY[ nPoint ] = cellCenterY[ mPoint ];
        cellCenterZ[ nPoint ] = cellCenterZ[ mPoint ];
        characterPoint[ nPoint ] = characterPoint[ mPoint ];

        //! p5.
        i = iCenterDimension[ iBlock ]-1;
        j = jCenterDimension[ iBlock ]-1;
        k = 0;
        mPoint = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*(k+1)+iCenterDimension[ iBlock ]*j+i;
        nPoint = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
        cellCenterX[ nPoint ] = cellCenterX[ mPoint ];
        cellCenterY[ nPoint ] = cellCenterY[ mPoint ];
        cellCenterZ[ nPoint ] = cellCenterZ[ mPoint ];
        characterPoint[ nPoint ] = characterPoint[ mPoint ];

        //! p6.
        i = iCenterDimension[ iBlock ]-1;
        j = 0;
        k = kCenterDimension[ iBlock ]-1;
        mPoint = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*(k-1)+iCenterDimension[ iBlock ]*j+i;
        nPoint = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
        cellCenterX[ nPoint ] = cellCenterX[ mPoint ];
        cellCenterY[ nPoint ] = cellCenterY[ mPoint ];
        cellCenterZ[ nPoint ] = cellCenterZ[ mPoint ];
        characterPoint[ nPoint ] = characterPoint[ mPoint ];

        //! p7.
        i = 0;
        j = jCenterDimension[ iBlock ]-1;
        k = kCenterDimension[ iBlock ]-1;
        mPoint = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*(k-1)+iCenterDimension[ iBlock ]*j+i;
        nPoint = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
        cellCenterX[ nPoint ] = cellCenterX[ mPoint ];
        cellCenterY[ nPoint ] = cellCenterY[ mPoint ];
        cellCenterZ[ nPoint ] = cellCenterZ[ mPoint ];
        characterPoint[ nPoint ] = characterPoint[ mPoint ];

        //! p8.
        i = iCenterDimension[ iBlock ]-1;
        j = jCenterDimension[ iBlock ]-1;
        k = kCenterDimension[ iBlock ]-1;
        mPoint = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*(k-1)+iCenterDimension[ iBlock ]*j+i;
        nPoint = startCenterLabel[ iBlock ]+iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*k+iCenterDimension[ iBlock ]*j+i;
        cellCenterX[ nPoint ] = cellCenterX[ mPoint ];
        cellCenterY[ nPoint ] = cellCenterY[ mPoint ];
        cellCenterZ[ nPoint ] = cellCenterZ[ mPoint ];
        characterPoint[ nPoint ] = characterPoint[ mPoint ];
    }
}

void MultigridManager::GenerateGlobalPointSpace()
{
    coordinateX.resize(globalPointNumber);
    coordinateY.resize(globalPointNumber);
    coordinateZ.resize(globalPointNumber);
}

void MultigridManager::GenerateGlobalCenterSpace()
{
    characterPoint.resize(globalCenterNumber);
    cellCenterX.resize(globalCenterNumber);
    cellCenterY.resize(globalCenterNumber);
    cellCenterZ.resize(globalCenterNumber);
}

void MultigridManager::GenerateLocalPointSpace()
{
    blockLocation.resize(globalPointNumber);
    iLocation.resize(globalPointNumber);
    jLocation.resize(globalPointNumber);
    kLocation.resize(globalPointNumber);
}

void MultigridManager::ReadGridCoordinate()
{
    GenerateGlobalBlockSpace();

    int pCounter = 0;
    for (int iBlock = 0; iBlock < globalBlockNumber; ++iBlock)
    {
        startPointLabel[iBlock] = pCounter;
        pCounter += OrdinaryGrid[iBlock]->GetNTotalNode();
    }
    globalPointNumber = pCounter;
    GenerateGlobalPointSpace();

    int iCounter = 0;
    for (int iBlock = 0; iBlock < globalBlockNumber; ++ iBlock)
    {
        startCenterLabel[ iBlock ] = iCounter;
        iCenterDimension[ iBlock ] = iDimension[ iBlock ]+1;
        jCenterDimension[ iBlock ] = jDimension[ iBlock ]+1;
        kCenterDimension[ iBlock ] = kDimension[ iBlock ]+1;
        iCounter += iCenterDimension[ iBlock ]*jCenterDimension[ iBlock ]*kCenterDimension[ iBlock ];
    }
    globalCenterNumber = iCounter;
    GenerateGlobalCenterSpace();

    for (int iBlock = 0; iBlock < globalBlockNumber; ++iBlock)
    {
        RDouble *x = OrdinaryGrid[iBlock]->GetX();
        RDouble *y = OrdinaryGrid[iBlock]->GetY();
        RDouble *z = OrdinaryGrid[iBlock]->GetZ();

        int localPointNumber = OrdinaryGrid[iBlock]->GetNTotalNode();
        for (pCounter = 0; pCounter < localPointNumber; ++ pCounter)
        {
            int globalPointLabel = startPointLabel[iBlock] + pCounter;
            coordinateX[globalPointLabel] = x[pCounter];
            coordinateY[globalPointLabel] = y[pCounter];
            coordinateZ[globalPointLabel] = z[pCounter];

        }
    }

    GenerateLocalPointSpace();

    int globalPointLabel = 0;
    for (int iBlock = 0; iBlock < globalBlockNumber; ++iBlock)
    {
        for (int k = 0; k < kDimension[iBlock]; ++k)
        {
            for (int j = 0; j < jDimension[iBlock]; ++j)
            {
                for (int i = 0; i < iDimension[iBlock]; ++i)
                {
                    blockLocation[globalPointLabel] = iBlock;
                    iLocation[globalPointLabel] = i;
                    jLocation[globalPointLabel] = j;
                    kLocation[globalPointLabel] = k;
                    globalPointLabel += 1;
                }
            }
        }
    }
}

void MultigridManager::GenerateGlobalNodeLinkStructure()
{
    ReadBoundaryCondition();
    InitialNodeLinksData();
    SetNewPatchModel();
}

void MultigridManager::ReadBoundaryCondition()
{
    GenerateGlobalFaceSpace();

    maxLocalFaceNumber = 0;
    specialBlockLabel = -1;

    for (int iBlock = 0; iBlock < globalBlockNumber; ++iBlock)
    {
        StructGrid * grid = StructGridCast(OrdinaryGrid[iBlock]);

        StructBCSet * structBCSet = grid->GetStructBCSet();
        overallFaceNumber[iBlock] = structBCSet->GetnBCRegion();
        
        int patchCounter = 0, innerCounter = 0, ComparePatchInner;
        int boundaryLabel;
        
        iDimension[iBlock] = grid->GetNI();
        jDimension[iBlock] = grid->GetNJ();
        kDimension[iBlock] = grid->GetNK();
        
        for (int iBCRegion = 0; iBCRegion < overallFaceNumber[iBlock]; ++ iBCRegion)
        {
            StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);

            boundaryLabel = structBC->GetBCType();
            if (boundaryLabel < 0)
            {
                patchCounter += 1;
            }
            else
            {
                innerCounter += 1;
            }
        }

        ComparePatchInner = max(patchCounter, innerCounter);
        if (ComparePatchInner > maxLocalFaceNumber)
        {
            maxLocalFaceNumber = ComparePatchInner;
            specialBlockLabel = iBlock;
        }
    }

    cout << " Grid Information :" << endl;
    cout << " maxLocalFaceNumber = " << maxLocalFaceNumber << endl;
    cout << " specialblockLabel  = " << specialBlockLabel << endl;

    GenerateLocalFaceSpace();
    
    for (int iBlock = 0; iBlock < globalBlockNumber; ++iBlock)
    {
        StructGrid * grid = StructGridCast(OrdinaryGrid[iBlock]);

        StructBCSet * structBCSet = grid->GetStructBCSet();

        int patchLabel = 0, innerLabel = 0;
        int imin, imax, jmin, jmax, kmin, kmax, boundaryLabel, targetLabel;

        for (int iBCRegion = 0; iBCRegion < overallFaceNumber[iBlock]; ++ iBCRegion)
        {
            StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);

            int * s_st = structBC->GetStartPoint();
            int * s_ed = structBC->GetEndPoint();

            imin = s_st[0];
            jmin = s_st[1];
            kmin = s_st[2];

            imax = s_ed[0];
            jmax = s_ed[1];
            kmax = s_ed[2];

            boundaryLabel = structBC->GetBCType();
            if (boundaryLabel < 0)
            {
                iminMaster[iBlock][patchLabel] = imin;
                imaxMaster[iBlock][patchLabel] = imax;
                jminMaster[iBlock][patchLabel] = jmin;
                jmaxMaster[iBlock][patchLabel] = jmax;
                kminMaster[iBlock][patchLabel] = kmin;
                kmaxMaster[iBlock][patchLabel] = kmax;

                int * t_st = structBC->GetTargetStart();
                int * t_ed = structBC->GetTargetEnd();

                imin = t_st[0];
                jmin = t_st[1];
                kmin = t_st[2];

                imax = t_ed[0];
                jmax = t_ed[1];
                kmax = t_ed[2];

                targetLabel = structBC->GetTargetRegionBlock();
                
                iminSlave[iBlock][patchLabel] = imin;
                imaxSlave[iBlock][patchLabel] = imax;
                jminSlave[iBlock][patchLabel] = jmin;
                jmaxSlave[iBlock][patchLabel] = jmax;
                kminSlave[iBlock][patchLabel] = kmin;
                kmaxSlave[iBlock][patchLabel] = kmax;
                tarPatch[iBlock][patchLabel] = targetLabel;
                patchLabel += 1;
            }
            else
            {
                iminInner[iBlock][innerLabel] = min(imin, imax);
                imaxInner[iBlock][innerLabel] = max(imin, imax);
                jminInner[iBlock][innerLabel] = min(jmin, jmax);
                jmaxInner[iBlock][innerLabel] = max(jmin, jmax);
                kminInner[iBlock][innerLabel] = min(kmin, kmax);
                kmaxInner[iBlock][innerLabel] = max(kmin, kmax);
                faceType[iBlock][innerLabel] = boundaryLabel;
                innerLabel += 1;
            }
        }
        externalFaceNumber[iBlock] = patchLabel;
        internalFaceNumber[iBlock] = innerLabel;
    }

    return;
}

void MultigridManager::GenerateGlobalFaceSpace()
{
    startPointLabel.resize(globalBlockNumber);
    iDimension.resize(globalBlockNumber);
    jDimension.resize(globalBlockNumber);
    kDimension.resize(globalBlockNumber);
    overallFaceNumber.resize(globalBlockNumber);
    externalFaceNumber.resize(globalBlockNumber);
    internalFaceNumber.resize(globalBlockNumber);
    processLabel.resize(globalBlockNumber);
}

void MultigridManager::InitialNodeLinksData()
{
    int iCounter = 0;
    for (int iBlock = 0; iBlock < globalBlockNumber; ++iBlock)
    {
        startPointLabel[iBlock] = iCounter;
        iCounter += iDimension[iBlock]*jDimension[iBlock]*kDimension[iBlock];
    }
    globalPointNumber = iCounter;

    for (int iPoint = 0; iPoint < globalPointNumber; ++iPoint)
    {
        patchIPoint.push_back(-1);
        patchJPoint.push_back(-1);
        patchKPoint.push_back(-1);
    }   //! Different from the initialization of PMB3D.
}

void MultigridManager::GenerateLocalFaceSpace()
{
    InitializeLocalFaceSpace(iminMaster, 0);
    InitializeLocalFaceSpace(imaxMaster, 0);
    InitializeLocalFaceSpace(jminMaster, 0);
    InitializeLocalFaceSpace(jmaxMaster, 0);
    InitializeLocalFaceSpace(kminMaster, 0);
    InitializeLocalFaceSpace(kmaxMaster, 0);

    InitializeLocalFaceSpace(iminSlave, 0);
    InitializeLocalFaceSpace(imaxSlave, 0);
    InitializeLocalFaceSpace(jminSlave, 0);
    InitializeLocalFaceSpace(jmaxSlave, 0);
    InitializeLocalFaceSpace(kminSlave, 0);
    InitializeLocalFaceSpace(kmaxSlave, 0);

    InitializeLocalFaceSpace(iminInner, 0);
    InitializeLocalFaceSpace(imaxInner, 0);
    InitializeLocalFaceSpace(jminInner, 0);
    InitializeLocalFaceSpace(jmaxInner, 0);
    InitializeLocalFaceSpace(kminInner, 0);
    InitializeLocalFaceSpace(kmaxInner, 0);

    InitializeLocalFaceSpace(faceType, -1000);
    InitializeLocalFaceSpace(tarPatch, -1000);
}

void MultigridManager::InitializeLocalFaceSpace(PHInt2D & myInt2D, int value)
{
    myInt2D.resize(globalBlockNumber);
    for (int iBlock = 0; iBlock < globalBlockNumber; ++iBlock)
    {
        myInt2D[iBlock].resize(maxLocalFaceNumber);
        for (int iFace = 0; iFace < maxLocalFaceNumber; ++iFace)
        {
            myInt2D[iBlock][iFace] = value;
        }
    }
}

void MultigridManager::SetNewPatchModel()
{
    vector< int > sourceToLocalAspectMapping(3), localToTargetAspectMapping(3), sourceToTargetAspectMapping(3);
    vector< int > relativeRate(3), colorSource(3), colorTarget(3), sst(3), sed(3), tst(3), ted(3), maps(3);
    int sside = 0, tside = 0;
    const int NEGATIVE_UNIT = -1, POSITIVE_UNIT = 1;

    for (int iBlock = 0; iBlock < globalBlockNumber; ++iBlock)
    {
        for (int iPatch = 0; iPatch < externalFaceNumber[iBlock]; ++iPatch)
        {
            sst[0] = iminMaster[iBlock][iPatch];
            sed[0] = imaxMaster[iBlock][iPatch];
            sst[1] = jminMaster[iBlock][iPatch];
            sed[1] = jmaxMaster[iBlock][iPatch];
            sst[2] = kminMaster[iBlock][iPatch];
            sed[2] = kmaxMaster[iBlock][iPatch];

            tst[0] = iminSlave[iBlock][iPatch];
            ted[0] = imaxSlave[iBlock][iPatch];
            tst[1] = jminSlave[iBlock][iPatch];
            ted[1] = jmaxSlave[iBlock][iPatch];
            tst[2] = kminSlave[iBlock][iPatch];
            ted[2] = kmaxSlave[iBlock][iPatch];

            //! Local analysis of the source grid block.
            for (int m = 0; m < 3; ++ m)
            {
                if (sst[m] == sed[m])
                {
                    sourceToLocalAspectMapping[m] = 0; //! The m-th axis is the relative main axis of the source grid block.

                    if (sst[m] == 1)
                    {
                        sside = 2 * m;
                        colorSource[m] = POSITIVE_UNIT; //! The patching face is the start point of the source grid block relative main axis. 
                    }
                    else
                    {
                        sside = 2 * m + 1;
                        colorSource[m] = NEGATIVE_UNIT; //! The patching face is the end point of the source grid block relative main axis. 
                    }
                }
                else if (sst[m] < 0)
                {
                    sourceToLocalAspectMapping[m] = 1; //! The m-th axis is the relative negative axis of the source grid block.

                    if (sst[m] < sed[m])
                    {
                        colorSource[m] = NEGATIVE_UNIT; //! The index on the the relative negative axis of the source grid block is from small to large.
                    }
                    else
                    {
                        colorSource[m] = POSITIVE_UNIT; //! The index on the the relative negative axis of the source grid block is from large to small.
                    }
                }
                else
                {
                    sourceToLocalAspectMapping[m] = 2; //! The m-th axis is the relative positive axis of the source grid block.

                    if (sst[m] > sed[m])
                    {
                        colorSource[m] = NEGATIVE_UNIT; //! The index on the the relative positive axis of the source grid block is from large to small.
                    }
                    else
                    {
                        colorSource[m] = POSITIVE_UNIT; //! The index on the the relative positive axis of the source grid block is from small to large.
                    }
                }

                //! Local analysis of the target grid block.
                if (tst[m] == ted[m])
                {
                    localToTargetAspectMapping[0] = m; //! The m-th axis is the relative main axis of the target grid block.

                    if (tst[m] == 1)
                    {
                        tside = 2 * m;
                        colorTarget[m] = NEGATIVE_UNIT; //! The patching face is the start point of the target grid block relative main axis.
                    }
                    else
                    {
                        tside = 2 * m + 1;
                        colorTarget[m] = POSITIVE_UNIT; //! The patching face is the end point of the target grid block relative main axis.
                    }
                }
                else if (tst[m] < 0)
                {
                    localToTargetAspectMapping[1] = m; //! The m-th axis is the relative negative axis of the target grid block.

                    if (tst[m] < ted[m])
                    {
                        colorTarget[m] = NEGATIVE_UNIT; //! The index on the the relative negative axis of the target grid block is from small to large.
                    }
                    else
                    {
                        colorTarget[m] = POSITIVE_UNIT; //! The index on the the relative negative axis of the target grid block is from large to small.
                    }
                }
                else
                {
                    localToTargetAspectMapping[2] = m; //! The m-th axis is the relative positive axis of the target grid block.

                    if (tst[m] > ted[m])
                    {
                        colorTarget[m] = NEGATIVE_UNIT; //! The index on the the relative positive axis of the target grid block is from large to small.
                    }
                    else
                    {
                        colorTarget[m] = POSITIVE_UNIT; //! The index on the the relative positive axis of the target grid block is from small to large.
                    }
                }
            }

            //! The analysis from source to target.
            for (int m = 0; m < 3; ++ m)
            {
                int i = sourceToLocalAspectMapping[m];
                int j = localToTargetAspectMapping[i];

                sourceToTargetAspectMapping[m] = j;
        
                //! Main axis analysis takes (s start-t end) and (s end-t start)as true, (s start-t start) and (s end-t end) as false.
                //! Negative axis analysis takes (s increase-t increase) and (s decrease-t decrease)as true, (s increase-t decrease) and (s decrease-t increase) as false.
                //! Positive axis analysis takes (s increase-t increase) and (s decrease-t decrease)as true, (s increase-t decrease) and (s decrease-t increase) as false.
                //! It can be concluded that taking colorSource equals colorTarget as true, colorSource unequals colorTarget as false.
                if (colorSource[m] == colorTarget[j])
                {
                    relativeRate[m] = POSITIVE_UNIT;
                }
                else
                {
                    relativeRate[m] = NEGATIVE_UNIT;
                }
            }

            PHSPACE::Standardize(sst, sed, sside);
            PHSPACE::Standardize(tst, ted, tside);

            int nbpat = tarPatch[iBlock][iPatch];
            int mk = sourceToTargetAspectMapping[2];
            int mj = sourceToTargetAspectMapping[1];
            int mi = sourceToTargetAspectMapping[0];

            for (int k = sst[2]; k <= sed[2]; ++ k)
            {
                int ks = k + 1;
                maps[mk] = (1 + relativeRate[2]) / 2 * tst[mk] + (1 - relativeRate[2]) / 2 * ted[mk] + relativeRate[2] * (k - sst[2]);
                for (int j = sst[1]; j <= sed[1]; ++ j)
                {
                    int js = j + 1; 
                    maps[mj] = (1 + relativeRate[1]) / 2 * tst[mj] + (1 - relativeRate[1]) / 2 * ted[mj] + relativeRate[1] * (j - sst[1]);
                    for (int i = sst[0]; i <= sed[0]; ++ i)
                    {
                        int is = i + 1;
                        maps[mi] = (1 + relativeRate[0]) / 2 * tst[mi] + (1 - relativeRate[0]) / 2 * ted[mi] + relativeRate[0] * (i - sst[0]);

                        int kt = maps[2] + 1;
                        int jt = maps[1] + 1;
                        int it = maps[0] + 1;

                        int ps = startPointLabel[iBlock] + iDimension[iBlock] * jDimension[iBlock] * ks + iDimension[iBlock] * js + is;
                        int pt = startPointLabel[nbpat]  + iDimension[nbpat]  * jDimension[nbpat]  * kt + iDimension[nbpat]  * jt + it;

                        int id = sside / 2;
                        if (id == 0)
                        {
                            patchIPoint[ps] = pt;
                        }
                        else if (id == 1)
                        {
                            patchJPoint[ps] = pt;
                        }
                        else
                        {
                            patchKPoint[ps] = pt;
                        }
                    }
                }
            }
        }
    }

    return;
}

void MultigridManager::GetOtherBlockDirections(int id1, int & id2, int & id3)
{
    if (id1 == 1)
    {
        id2 = 2;
        id3 = 3;
    } 

    if (id1 == 2)
    {
        id2 = 3;
        id3 = 1;
    } 

    if (id1 == 3)
    {
        id2 = 1;
        id3 = 2;
    }
}

void MultigridManager::ComputeContravariantComponent()
{
    PHDouble2D amer(3), aa;
    PHDouble1D xx(3), bb(3);
    for (int i = 0; i < 3; ++i)
    {
        amer[i].resize(3);
    }

    for (int iBlock = 0; iBlock < globalBlockNumber; ++iBlock)
    {
        int iHarf = iDimension[iBlock]/2;
        int jHarf = jDimension[iBlock]/2;
        int kHarf = kDimension[iBlock]/2;

        int originalLabel = startPointLabel[iBlock];
        int iDirectLabel = originalLabel+iHarf;
        int jDirectLabel = originalLabel+iDimension[iBlock]*jHarf;
        int kDirectLabel = originalLabel+iDimension[iBlock]*jDimension[iBlock]*kHarf;

        amer[0][0] = coordinateX[iDirectLabel]-coordinateX[originalLabel];
        amer[0][1] = coordinateX[jDirectLabel]-coordinateX[originalLabel];
        amer[0][2] = coordinateX[kDirectLabel]-coordinateX[originalLabel];

        amer[1][0] = coordinateY[iDirectLabel]-coordinateY[originalLabel];
        amer[1][1] = coordinateY[jDirectLabel]-coordinateY[originalLabel];
        amer[1][2] = coordinateY[kDirectLabel]-coordinateY[originalLabel];

        amer[2][0] = coordinateZ[iDirectLabel]-coordinateZ[originalLabel];
        amer[2][1] = coordinateZ[jDirectLabel]-coordinateZ[originalLabel];
        amer[2][2] = coordinateZ[kDirectLabel]-coordinateZ[originalLabel];

        for (int k = 0; k < kCenterDimension[iBlock]; ++k)
        {
            for (int j = 0; j < jCenterDimension[iBlock]; ++j)
            {
                for (int i = 0; i < iCenterDimension[iBlock]; ++i)
                {
                    int core = startCenterLabel[iBlock]+iCenterDimension[iBlock]*jCenterDimension[iBlock]*k+iCenterDimension[iBlock]*j+i;
                    bb[0] = cellCenterX[core]-coordinateX[originalLabel];
                    bb[1] = cellCenterY[core]-coordinateY[originalLabel];
                    bb[2] = cellCenterZ[core]-coordinateZ[originalLabel];
                    aa = amer;
                    int flag;
                    GaussElimination(aa, bb, 3, xx, flag);
                    if (flag == 0)
                    {
                        TK_Exit::ExceptionExit("Linear correlation !\n");
                    }

                    cellCenterX[core] = xx[0];
                    cellCenterY[core] = xx[1];
                    cellCenterZ[core] = xx[2];
                }
            }
        }
    }

    return;
}

void MultigridManager::OutputLinkData()
{
    string linkFileName = GetLinkFileName();
    ofstream outfile(linkFileName.c_str(), ios_base::binary);
    outfile.write((char *)&globalBlockNumber, sizeof(int));

    for (int iBlock = 0; iBlock < globalBlockNumber; ++iBlock)
    {
        outfile.write((char *)&startPointLabel[iBlock], sizeof(int));
    }

    outfile.write((char *)&globalPointNumber, sizeof(int));
    
    for (int ip = 0; ip < globalPointNumber; ++ip)
    {
        outfile.write((char *)&patchIPoint[ip], sizeof(int));
    }

    for (int ip = 0; ip < globalPointNumber; ++ip)
    {
        outfile.write((char *)&patchJPoint[ip], sizeof(int));
    }

    for (int ip = 0; ip < globalPointNumber; ++ip)
    {
        outfile.write((char *)&patchKPoint[ip], sizeof(int));
    }

    for (int iBlock = 0; iBlock < globalBlockNumber; ++iBlock)
    {
        outfile.write((char *)&startCenterLabel[iBlock], sizeof(int));
    }

    outfile.write((char *)&globalCenterNumber, sizeof(int));

    for (int iCenter = 0; iCenter < globalCenterNumber; ++iCenter)
    {
        outfile.write((char *)&characterPoint[iCenter], sizeof(int));
    }

    for (int iCenter = 0; iCenter < globalCenterNumber; ++iCenter)
    {
        outfile.write((char *)&cellCenterX[iCenter], sizeof(RDouble));
    }

    for (int iCenter = 0; iCenter < globalCenterNumber; ++iCenter)
    {
        outfile.write((char *)&cellCenterY[iCenter], sizeof(RDouble));
    }

    for (int iCenter = 0; iCenter < globalCenterNumber; ++iCenter)
    {
        outfile.write((char *)&cellCenterZ[iCenter], sizeof(RDouble));
    }

    outfile.close();
}

string MultigridManager::GetLinkFileName()
{
    string linkFileName;

    GlobalDataBase::GetData("linkFileName", &linkFileName, PHSTRING,1);

    return linkFileName;
}

void MultigridManager::GenerateLinkStructure()
{
    GenerateGlobalNodeLinkStructure();

    GenerateGlobalCellLinkStructure();

    ComputeContravariantComponent();    //! Very important part Guo Yongheng 20151204. Contravariant Component Transformation.

    OutputLinkData();

    return;
}

void MultigridManager::GenerateHoleStructure()
{
    int codeOfDigHoles;
    GlobalDataBase::GetData("codeOfDigHoles", &codeOfDigHoles, PHINT,1);
    if (codeOfDigHoles == 0) return;

    ReadZoneInverseMapping();

    ReadOriginalHoleFile();

    ComputeLocalCoordinate();

    OutputNewHoleFile();

    return;
}

void MultigridManager::GenerateHoleGridStructure()
{
    if (! GlobalDataBase::GetIntParaFromDB("codeOfDigHoles")) return;

    ReadZoneInverseMapping();

    ReadOriginalHoleGrid();

    ComputeLocalHoleGrid();

    OutputNewHoleGrid();

    return;
}

void MultigridManager::OutputNewHoleGrid()
{
    string holeFullFileName = GlobalDataBase::GetStrParaFromDB("holeFullFileName");
    fstream file(holeFullFileName.c_str(), ios_base::out|ios_base::binary);

    file.write((char *)&numberOfHoles, sizeof(int));

    file.close();
    file.clear();
    return;
}

void MultigridManager::ComputeLocalHoleGrid()
{

}

void MultigridManager::ReadOriginalHoleGrid()
{

}

void MultigridManager::ReadOriginalHoleFile()
{
    string holeGridFileName;
    string holeBasicFileName = PHSPACE::GetHoleBasicFileName();
    fstream basicFile(holeBasicFileName.c_str(), ios_base::in);
    basicFile>>numberOfHoles;

    holeZoneStart.resize(numberOfHoles);
    holeZoneEnd.resize(numberOfHoles);
    holeZoneAttached.resize(numberOfHoles);
    numberOfHoleSurfaces.resize(numberOfHoles);
    nsi.resize(numberOfHoles);
    nsj.resize(numberOfHoles);
    spst.resize(numberOfHoles);
    numberOfHolePoints.resize(numberOfHoles);
    xsurf.resize(numberOfHoles);
    ysurf.resize(numberOfHoles);
    zsurf.resize(numberOfHoles);

    for (int iHole = 0; iHole < numberOfHoles; ++iHole)
    {
        basicFile>>holeGridFileName;
        basicFile>>holeZoneStart[iHole]>>holeZoneEnd[iHole]>>holeZoneAttached[iHole];

        fstream holeFile(holeGridFileName.c_str(), ios_base::in|ios_base::binary);
        holeFile.read((char *)&numberOfHoleSurfaces[iHole], sizeof(int));
        nsi[iHole].resize(numberOfHoleSurfaces[iHole]);
        nsj[iHole].resize(numberOfHoleSurfaces[iHole]);
        spst[iHole].resize(numberOfHoleSurfaces[iHole]);
        
        int kSurfaceDimension, pCounter = 0;
        for (int iSurf = 0; iSurf < numberOfHoleSurfaces[iHole]; ++iSurf)
        {
            holeFile.read((char *)&nsi[iHole][iSurf], sizeof(int));
            holeFile.read((char *)&nsj[iHole][iSurf], sizeof(int));
            holeFile.read((char *)&kSurfaceDimension, sizeof(int));
            spst[iHole][iSurf] = pCounter;
            pCounter += nsi[iHole][iSurf] * nsj[iHole][iSurf];
        }

        numberOfHolePoints[iHole] = pCounter;
        xsurf[iHole].resize(numberOfHolePoints[iHole]);
        ysurf[iHole].resize(numberOfHolePoints[iHole]);
        zsurf[iHole].resize(numberOfHolePoints[iHole]);

        for (int iSurf = 0; iSurf < numberOfHoleSurfaces[iHole]; ++iSurf)
        {
            int surfNodeNumber = nsi[iHole][iSurf]*nsj[iHole][iSurf];
            int startNode = spst[iHole][iSurf];
            for (int ip = 0; ip < surfNodeNumber; ++ip)
            {
                holeFile.read((char *)&xsurf[iHole][startNode+ip], sizeof(RDouble));
            }

            for (int ip = 0; ip < surfNodeNumber; ++ip)
            {
                holeFile.read((char *)&ysurf[iHole][startNode+ip], sizeof(RDouble));
            }

            for (int ip = 0; ip < surfNodeNumber; ++ip)
            {
                holeFile.read((char *)&zsurf[iHole][startNode+ip], sizeof(RDouble));
            }
        }
        holeFile.close();
    }
    basicFile.close();

    return;
}

void MultigridManager::ReadZoneInverseMapping()
{
    int nZones = globalBlockNumber;
    zoneInverseIndex.resize(nZones);
    string zoneInverseFileName = GlobalDataBase::GetStrParaFromDB("zoneInverseFileName");
    fstream file(zoneInverseFileName.c_str(), ios_base::in|ios_base::binary);
    file.read((char *)&zoneInverseIndex[0], nZones * sizeof(int));
    file.close();
    file.clear();
    return;
}

void MultigridManager::ComputeLocalCoordinate()
{
    PHDouble1D bb(3), xx(3);
    PHDouble2D amer(3), aa;
    for (int i = 0; i < 3; ++i)
    {
        amer[i].resize(3);
    }

    for (int iHole = 0; iHole < numberOfHoles; ++iHole)
    {
        int iZone = zoneInverseIndex[holeZoneAttached[iHole]];
        Grid * gridIn = OrdinaryGrid[iZone];
        StructGrid * grid = PHSPACE::StructGridCast(gridIn);
        PHDouble2D & localFrame = grid->GetLocalFrame();

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                amer[i][j] = localFrame[i][j+1];
            }
        }

        for (int ip = 0; ip < numberOfHolePoints[iHole]; ++ip)
        {
            bb[0] = xsurf[iHole][ip]-localFrame[0][0];
            bb[1] = ysurf[iHole][ip]-localFrame[1][0];
            bb[2] = zsurf[iHole][ip]-localFrame[2][0];

            aa = amer;
            int flag;
            GaussElimination(aa, bb, 3, xx, flag);
            if (flag == 0)
            {
                TK_Exit::ExceptionExit("Linear correlation !\n");
            }

            xsurf[iHole][ip] = xx[0];
            ysurf[iHole][ip] = xx[1];
            zsurf[iHole][ip] = xx[2];
        }
    }
    return;
}

void MultigridManager::OutputNewHoleFile()
{
    string holeFullFileName;
    GlobalDataBase::GetData("holeFullFileName", &holeFullFileName, PHSTRING, 1);
    fstream file(holeFullFileName.c_str(), ios_base::out|ios_base::binary);
    file.write((char *)&numberOfHoles, sizeof(int));

    for (int iHole = 0; iHole < numberOfHoles; ++iHole)
    {
        file.write((char *)&holeZoneStart[iHole], sizeof(int));
    }

    for (int iHole = 0; iHole < numberOfHoles; ++iHole)
    {
        file.write((char *)&holeZoneEnd[iHole], sizeof(int));
    }

    for (int iHole = 0; iHole < numberOfHoles; ++iHole)
    {
        file.write((char *)&holeZoneAttached[iHole], sizeof(int));
    }

    for (int iHole = 0; iHole < numberOfHoles; ++iHole)
    {
        file.write((char *)&numberOfHoleSurfaces[iHole], sizeof(int));
    }

    for (int iHole = 0; iHole < numberOfHoles; ++iHole)
    {
        file.write((char *)&numberOfHolePoints[iHole], sizeof(int));
    }

    for (int iHole = 0; iHole < numberOfHoles; ++iHole)
    {
        for (int is = 0; is < numberOfHoleSurfaces[iHole]; ++is)
        {
            file.write((char *)&nsi[iHole][is], sizeof(int));
        }

        for (int is = 0; is < numberOfHoleSurfaces[iHole]; ++is)
        {
            file.write((char *)&nsj[iHole][is], sizeof(int));
        }

        for (int is = 0; is < numberOfHoleSurfaces[iHole]; ++is)
        {
            file.write((char *)&spst[iHole][is], sizeof(int));
        }

        for (int ip = 0; ip < numberOfHolePoints[iHole]; ++ip)
        {
            file.write((char *)&xsurf[iHole][ip], sizeof(RDouble));
        }

        for (int ip = 0; ip < numberOfHolePoints[iHole]; ++ip)
        {
            file.write((char *)&ysurf[iHole][ip], sizeof(RDouble));
        }

        for (int ip = 0; ip < numberOfHolePoints[iHole]; ++ip)
        {
            file.write((char *)&zsurf[iHole][ip], sizeof(RDouble));
        }
    }
    file.close();
    return;
}

void Standardize(vector< int > & st, vector< int > & ed, int side)
{
    for (int i = 0; i < 3; ++ i)
    {
        st[i] = abs(st[i]) - 1;
        ed[i] = abs(ed[i]) - 1;

        if (st[i] > ed[i])
        {
            int se = st[i];
            st[i]  = ed[i];
            ed[i]  = se;
        }
        ed[i] -= 1;
    }

    int id = side / 2;
    int ic = side % 2;

    if (ic == 0)
    {
        ed[id] += 1;
    }
    else
    {
        st[id] -= 1;
    }

    return;
}

string GetHoleBasicFileName()
{
    string holeBasicFileName;
    GlobalDataBase::GetData("holeBasicFileName", &holeBasicFileName, PHSTRING,1);
    return holeBasicFileName;
}

void GaussElimination(PHDouble2D & a, PHDouble1D & b, int n, PHDouble1D & x, int & flag)
{
    int is = 0;
    flag = 1;
    PHInt1D js(n);

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

        if (flag == 0)
        {
            return;
        }

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

}


