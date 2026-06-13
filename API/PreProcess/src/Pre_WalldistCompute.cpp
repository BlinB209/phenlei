#include "Pre_WalldistCompute.h"
#include "Geo_StructBC.h"
#include "Geo_UnstructBC.h"
#include "Glb_Dimension.h"
#include "Math_BasisFunction.h"
#include "TK_Exit.h"
#include "Geo_StructGrid.h"
#include "Geo_Grid.h"
#include "Geo_UnstructGrid.h"
#include "Pre_HDF5File.h"
#include "AleManager.h"
#ifdef USE_TecplotLib
#include "TECXXX.h"
#endif
#pragma warning(disable:6386)
#pragma warning(disable:6385)
#pragma warning(disable:26451)

using namespace std;

namespace PHSPACE
{
LIB_EXPORT Pre_WalldistCompute::Pre_WalldistCompute(int walldistComputeMethod, int level)
{
    this->wallstructureList = new vector<WallStructure *>;
    this->walldistComputeMethod = walldistComputeMethod;
    this->level = level;
    this->myid = PHMPI::GetCurrentProcessorID();
    nTWFace = 0;
    wallFaceTree = 0;
    wallFaceKDTree = 0;
    wallNodeKDTree = 0;
    rotateDegree = 33.33 * PI / 180.0;
    presentBlock = 0;
    *wallMin = 0;
    *wallMax = 0;
    wallTolerance = 0.0;
    nTWallFace = 0;
    rotateAxis = 0;

    cellMethodOrNodeMethod = GlobalDataBase::GetIntParaFromDB("cellMethodOrNodeMethod");
}

LIB_EXPORT Pre_WalldistCompute::~Pre_WalldistCompute()
{
    for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
    {
        WallStructure *walldist = (*wallstructureList)[iWall];
        delete walldist;
    }
    wallstructureList->clear();
    delete wallstructureList; wallstructureList = NULL;

    WriteLogFile("  Wall face structure deleting ...\n");
    delete wallFaceTree;
}

LIB_EXPORT void Pre_WalldistCompute::Run()
{
    WriteLogFile("  Start Wall Distance Computing ...");

    FillWallStructureNew();

    InitWallDistance();

    int systemGridType = GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    int dimension = PHSPACE::GetDim();
    if (dimension == TWO_D)
    {
        walldistComputeMethod = ACCURATE_METHOD;
    }

    TimeSpan *timeSpan = new TimeSpan();

    int parallelAxis = JudgeInitialRotateAxis();
    this->rotateAxis = 1;
    while (!IsWalldistComputOver())
    {
        WriteLogFile("Wall dist computing for axis: ", this->rotateAxis);
        PrintToWindow("Wall dist computing for axis: ", this->rotateAxis, "\n");

        if (this->rotateAxis == 5) break;

        if (this->rotateAxis == parallelAxis)
        {
            ++ this->rotateAxis;
            continue;
        }

        if (walldistComputeMethod == FAST_METHOD || walldistComputeMethod == SUPER_FAST_METHOD)
        {
            BuildWallStructTree();
        }
        else if (walldistComputeMethod == KDTREE_METHOD && systemGridType == STRUCTGRID)
        {
            BuildWallStructKDTree();
        }

        ComputeWallDistance();

        if (walldistComputeMethod != FAST_METHOD) break;

        ++ this->rotateAxis;
        delete wallFaceTree;    wallFaceTree = nullptr;

        if (walldistComputeMethod == KDTREE_METHOD && systemGridType == STRUCTGRID)
        {
            delete wallFaceKDTree;    wallFaceKDTree = nullptr;
        }

        if (this->rotateAxis == 4)
        {
            //! The last search time.
            walldistComputeMethod = ACCURATE_METHOD;
        }

        if (this->rotateAxis == 5) break;
    }

    PrintToWindow("    Wall distance computed over, waiting for other processors ...\n");

    CheckWallDistance();

    PostWallDistance();

    timeSpan->ShowTimeSpanToLogFile("Wall distance computing");
    WriteLogFile("  End Wall Distance Computing ...");
    delete timeSpan;    timeSpan = nullptr;
}

LIB_EXPORT void Pre_WalldistCompute::ComputeNodeDistance()
{
    FillWallStructureNew();

    ComputeNodeDistanceByKDTreeMethod();

    PrintToWindow("    Node Wall distance computed over, waiting for other processors ...\n");

    PH_Barrier();

    WriteLogFile("  End Grid Node Distance To Wall Compute ...");
}

LIB_EXPORT void Pre_WalldistCompute::ComputeNodeDistanceInOtherBlock()
{
    FillWallStructureNew();

    ComputeNodeDistanceByKDTreeMethodInOtherBlock();
}

void Pre_WalldistCompute::FillWallStructureNew()
{
    TimeSpan *timeSpan = new TimeSpan();

    ServerCollection();
    timeSpan->ShowTimeSpanToLogFile("Wall structure collection to server");

    ServerBcast();
    timeSpan->ShowTimeSpanToLogFile("Wall structure BCasting from server");

    PrintToWindow("Barrier for wall structure filling ...\n");
    PrintToWindow("End for wall structure filling ...\n");
    delete timeSpan;    timeSpan = nullptr;
}

void Pre_WalldistCompute::FillWallStructure()
{
    WriteLogFile("  Start to fill wall structure ...");

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int tag  = GetSendRecvTag(GRID_BASED, 0, iZone);
        int proc = GetZoneProcessorID(iZone);

        DataContainer *cdata = new DataContainer();

        if (myid == proc)
        {
            Grid *grid = GetGrid(iZone, this->level);
            CompressWallStructure(cdata, grid);
        }

        PH_Bcast(cdata, proc, tag);

        DeCompressWallStructure(cdata, true);

        delete cdata;    cdata = nullptr;
    }

    WriteLogFile("  End filling wall structure ...");
}

void Pre_WalldistCompute::InitWallDistance()
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Grid *grid = GetGrid(zoneID, this->level);
        InitWallDistance(grid);
    }
}

void Pre_WalldistCompute::InitWallDistance(Grid *grid)
{
    if (grid->Type() == PHSPACE::STRUCTGRID)
    {
        StructGrid *structGrid = StructGridCast(grid);
        InitWallDistance(structGrid);
    }
    else
    {
        UnstructGrid *unstructGrid = UnstructGridCast(grid);
        InitWallDistance(unstructGrid);
    }
}

void Pre_WalldistCompute::InitWallDistance(StructGrid *grid)
{
    RDouble3D &walldist = *grid->GetWallDist();
    RDouble3D &nearestwallfacenormalx = *grid->GetNearestWallFaceNormalX();
    RDouble3D &nearestwallfacenormaly = *grid->GetNearestWallFaceNormalY();
    RDouble3D &nearestwallfacenormalz = *grid->GetNearestWallFaceNormalZ();

    //! Initialization.
    walldist = LARGE;
    nearestwallfacenormalx = LARGE;
    nearestwallfacenormaly = LARGE;
    nearestwallfacenormalz = LARGE;
}

void Pre_WalldistCompute::InitWallDistance(UnstructGrid *grid)
{
    RDouble *walldist = grid->GetWallDist();
    RDouble *wallDistNode = grid->GetWallDistNode();
    RDouble *nearestwallfacenormalx = grid->GetNearestWallFaceNormalX();
    RDouble *nearestwallfacenormaly = grid->GetNearestWallFaceNormalY();
    RDouble *nearestwallfacenormalz = grid->GetNearestWallFaceNormalZ();

    int nTotalNode = grid->GetNTotalNode();
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    //! Initialization.
    PHSPACE::SetField(walldist, LARGE, nTotal);
    PHSPACE::SetField(wallDistNode, LARGE, nTotalNode);
    PHSPACE::SetField(nearestwallfacenormalx, LARGE, nTotal);
    PHSPACE::SetField(nearestwallfacenormaly, LARGE, nTotal);
    PHSPACE::SetField(nearestwallfacenormalz, LARGE, nTotal);

    int **faceNodeIndex = grid->GetFace2NodeArray();
    int *faceNodeNumber = grid->GetNodeNumberOfEachFace();

    UnstructBCSet **bcr = grid->GetBCRecord();
    
    for (int iFace = 0; iFace < nBoundFace; ++iFace)
    {
        int bcType = bcr[iFace]->GetKey();
        if (IsWall(bcType))
        {
            int nNode = faceNodeNumber[iFace];
            for (int iNode = 0; iNode < nNode; ++iNode)
            {
                int nodeID = faceNodeIndex[iFace][iNode];
                wallDistNode[nodeID] = 0.0;
            }
        }
    }
}

void Pre_WalldistCompute::ComputeWallDistance()
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    PrintToWindow("  Number of wall structure:", wallstructureList->size(), "\n");

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Grid *grid = GetGrid(zoneID, this->level);

        ComputeWallDistance(grid);
    }
}

void Pre_WalldistCompute::ComputeWallDistance(Grid *grid)
{
    if (grid->Type() == PHSPACE::STRUCTGRID)
    {
        StructGrid *structGrid = StructGridCast(grid);
        if (walldistComputeMethod == Pre_WalldistCompute::KDTREE_METHOD)
        {
            ComputeWallDistanceKDTree(structGrid);
        }
        else if (walldistComputeMethod == Pre_WalldistCompute::FAST_METHOD)
        {
            ComputeWallDist3DByProjectingFast(structGrid);
        }
        else
        {
            ComputeWallDistance(structGrid);
        }
    }
    else
    {
        UnstructGrid *unstructGrid = UnstructGridCast(grid);
        ComputeWallDistance(unstructGrid);
    }
}

void Pre_WalldistCompute::ComputeWallDistance(UnstructGrid *grid)
{
    if (grid->GetDim() == PHSPACE::TWO_D)
    {
        ComputeWallDist2DByBinarySearching(grid);
    }
    else
    {
        if (walldistComputeMethod == Pre_WalldistCompute::ACCURATE_METHOD)
        {
            if (cellMethodOrNodeMethod == CELL_METHOD)
            {
                ComputeWallDist3DByProjectingAccurate(grid);
            }
            else
            {
                ComputeWallDist3DByProjectingAccurateAndNodeMethod(grid);
            }
        }
        else
        {
            if (cellMethodOrNodeMethod == CELL_METHOD)
            {
                ComputeWallDist3DByProjectingFast(grid);
            }
            else
            {
                ComputeWallDist3DByProjectingFastAndNodeMethod(grid);
            }
        }
    }
}

void Pre_WalldistCompute::ComputeWallDist2DByBinarySearching(UnstructGrid *grid)
{
    using namespace PHMPI;

    RDouble *walldist = grid->GetWallDist();

    RDouble *nearestwallfacenormalx = grid->GetNearestWallFaceNormalX();
    RDouble *nearestwallfacenormaly = grid->GetNearestWallFaceNormalY();
    RDouble *nearestwallfacenormalz = grid->GetNearestWallFaceNormalZ();

    int nTotalCell = grid->GetNTotalCell();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
    int *wallCellLabel = grid->GetWallCellLabel();

    cout.unsetf(ios::scientific);
    cout.setf(ios::fixed);
    cout.precision(0);

    int interval = GetProgressInterval(nTotalCell, 10);

#ifdef AI_RandomForestRegressor
    RDouble faceNormal[3];
    vector<RDouble> icellToWallNormal(3);
    vector<vector<RDouble>> allCellToWallNormal;
#endif

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (viscousType <= LAMINAR && wallCellLabel[iCell] != 1 && iLES == NOLES_SOLVER) continue;

        if (myid == 0 && iCell > 0 && iCell % interval == 0)
        {
            ProgressMonitorToWindows(iCell, nTotalCell, "wall distance");
        }

        RDouble xcc0 = xcc[iCell];
        RDouble ycc0 = ycc[iCell];
        RDouble zcc0 = zcc[iCell];

        for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
        {
            WallStructure *cwalldist = (*wallstructureList)[iWall];

            uint_t nWallFaces = cwalldist->GetNumberOfWallFaces();
            WallStructure::value_type &xfc = cwalldist->GetXFaceCenter();
            WallStructure::value_type &yfc = cwalldist->GetYFaceCenter();
            WallStructure::value_type &zfc = cwalldist->GetZFaceCenter();
            WallStructure::value_type &x = cwalldist->GetX();
            WallStructure::value_type &y = cwalldist->GetY();
            WallStructure::value_type &z = cwalldist->GetZ();

#ifdef AI_RandomForestRegressor
            WallStructure::value_type &xfnAI = cwalldist->GetXFaceNormal();
            WallStructure::value_type &yfnAI = cwalldist->GetYFaceNormal();
            WallStructure::value_type &zfnAI = cwalldist->GetZFaceNormal();
#endif           
            int *nPointPerFace = cwalldist->GetnPointPerFace();
            int *wallFace2Node = cwalldist->GetWallFace2Node();

            vector<vector<int> > face2node;
            ConstructFace2Node(nPointPerFace, wallFace2Node, nWallFaces, face2node);

            int markFace = -1;
            for (uint_t iFace = 0; iFace < nWallFaces; ++ iFace)
            {
                RDouble &xMid = xfc[iFace];
                RDouble &yMid = yfc[iFace];
                RDouble &zMid = zfc[iFace];

#ifdef AI_RandomForestRegressor
                RDouble &xNormal = xfnAI[iFace];
                RDouble &yNormal = yfnAI[iFace];
                RDouble &zNormal = zfnAI[iFace];
#endif

                //! The following computing would be executed two times for 2d case,
                //! but it's regardless, because the time consuming is little.
                uint_t numberofPointsOfTheFace = face2node[iFace].size();
                for (int index = 0; index < numberofPointsOfTheFace; ++ index)
                {
                    int iStart = index;
                    int iEnd = (index + 1) % numberofPointsOfTheFace;

                    int startPoint = face2node[iFace][iStart];
                    int endPoint   = face2node[iFace][iEnd];

                    const RDouble xStart = x[startPoint];
                    const RDouble yStart = y[startPoint];
                    const RDouble zStart = z[startPoint];

                    const RDouble xEnd = x[endPoint];
                    const RDouble yEnd = y[endPoint];
                    const RDouble zEnd = z[endPoint];

                    RDouble dx = xEnd - xStart;
                    RDouble dy = yEnd - yStart;
                    RDouble dz = zEnd - zStart;
                    RDouble lengthOfTwoPoints = SQR(dx, dy, dz);
                    if (lengthOfTwoPoints <= 1.0e-16)
                    {
                        //! Two points are overlapped.
                        continue;
                    }

                    RDouble wdst = ComputeDistanceOfANodeToSegment(xcc0, ycc0, zcc0, xStart, yStart, zStart,
                                                                   xMid, yMid, zMid, xEnd, yEnd, zEnd);

                    if (walldist[iCell] > wdst)
                    {
                        walldist[iCell] = wdst;
                        markFace = static_cast<int>(iFace);

#ifdef AI_RandomForestRegressor
                        icellToWallNormal[0] = xNormal;
                        icellToWallNormal[1] = yNormal;
                        icellToWallNormal[2] = zNormal;
#endif
                    }
                }
            }

#ifdef AI_RandomForestRegressor
            allCellToWallNormal.push_back(icellToWallNormal);
#endif

            if (markFace <= -1)
            {
                continue;
            }

            WallStructure::value_type& xfn = cwalldist->GetXFaceNormal();
            WallStructure::value_type& yfn = cwalldist->GetYFaceNormal();
            WallStructure::value_type& zfn = cwalldist->GetZFaceNormal();

            nearestwallfacenormalx[iCell] = xfn[markFace];
            nearestwallfacenormaly[iCell] = yfn[markFace];
            nearestwallfacenormalz[iCell] = zfn[markFace];
        }
    }
#ifdef AI_RandomForestRegressor
    DumpCellToWallNormal(allCellToWallNormal);
#endif
}

void Pre_WalldistCompute::ComputeWallDist3DByProjectingAccurateAndNodeMethod(UnstructGrid *grid)
{
    using namespace PHMPI;

    RDouble *wallDistNode = grid->GetWallDistNode();

    int nTotalNode = grid->GetNTotalNode();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble *nearestwallfacenormalx = grid->GetNearestWallFaceNormalX();
    RDouble *nearestwallfacenormaly = grid->GetNearestWallFaceNormalY();
    RDouble *nearestwallfacenormalz = grid->GetNearestWallFaceNormalZ();

    RDouble nodeCoord[3], faceCenter[3], faceNormal[3];

    int interval = GetProgressInterval(nTotalNode, 10);
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (iNode > 0 && iNode % interval == 0)
        {
            ProgressMonitorToWindows(iNode, nTotalNode, "wall distance");
            ProgressMonitorToLogFile(iNode, nTotalNode, "wall distance");

            if (rotateAxis >= 3 || timeSpan.GetTimeSpan() > 3600.0)
            {
                ostringstream oss;
                oss.unsetf(ios::scientific);
                oss.setf(ios::fixed);
                oss.precision(0);
                oss << "     rotateAxis = " << this->rotateAxis << ", " << floor((iNode * 1.0) / nTotalNode * 100 + 0.5) << "% " << "wall distance" << " completed ...";
                WriteLogFile(oss, false);
            }
        }

        if (wallDistNode[iNode] < LARGE * 0.9)
        {
            continue;
        }

        nodeCoord[0] = x[iNode];
        nodeCoord[1] = y[iNode];
        nodeCoord[2] = z[iNode];

        for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
        {
            WallStructure *cwalldist = (*wallstructureList)[iWall];

            uint_t nWallFaces              = cwalldist->GetNumberOfWallFaces();
            WallStructure::value_type &xfc = cwalldist->GetXFaceCenter();
            WallStructure::value_type &yfc = cwalldist->GetYFaceCenter();
            WallStructure::value_type &zfc = cwalldist->GetZFaceCenter();
            WallStructure::value_type &xW   = cwalldist->GetX();
            WallStructure::value_type &yW   = cwalldist->GetY();
            WallStructure::value_type &zW   = cwalldist->GetZ();
            WallStructure::value_type &xfn = cwalldist->GetXFaceNormal();
            WallStructure::value_type &yfn = cwalldist->GetYFaceNormal();
            WallStructure::value_type &zfn = cwalldist->GetZFaceNormal();

            int *nPointPerFace = cwalldist->GetnPointPerFace();
            int *wallFace2Node = cwalldist->GetWallFace2Node();

            int nodeCount = 0;
            int markFace = -1;
            for (int iFace = 0; iFace < nWallFaces; ++ iFace)
            {
                faceCenter[0] = xfc[iFace];
                faceCenter[1] = yfc[iFace];
                faceCenter[2] = zfc[iFace];

                faceNormal[0] = xfn[iFace];
                faceNormal[1] = yfn[iFace];
                faceNormal[2] = zfn[iFace];

                vector<int> face2node;
                face2node.resize(nPointPerFace[iFace]);
                for (int jNode = 0; jNode < nPointPerFace[iFace]; ++jNode)
                {
                    face2node[jNode] = wallFace2Node[nodeCount++];
                }

                RDouble wdst = ComputeDistanceOfANodeToAFace(nodeCoord, faceCenter, faceNormal, face2node, xW, yW, zW);
                if (wallDistNode[iNode] > wdst)
                {
                    wallDistNode[iNode] = wdst;
                    markFace = iFace;
                }
            }

            if (markFace <= -1)
            {
                return;
            }

            nearestwallfacenormalx[iNode] = xfn[markFace];
            nearestwallfacenormaly[iNode] = yfn[markFace];
            nearestwallfacenormalz[iNode] = zfn[markFace];

        }
    }
}

void Pre_WalldistCompute::ComputeWallDist3DByProjectingAccurate(UnstructGrid *grid)
{
    using namespace PHMPI;

    RDouble *walldist = grid->GetWallDist();

    RDouble *nearestwallfacenormalx = grid->GetNearestWallFaceNormalX();
    RDouble *nearestwallfacenormaly = grid->GetNearestWallFaceNormalY();
    RDouble *nearestwallfacenormalz = grid->GetNearestWallFaceNormalZ();

    int nTotalCell = grid->GetNTotalCell();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
    int *wallCellLabel = grid->GetWallCellLabel();

    cout.unsetf(ios::scientific);
    cout.setf(ios::fixed);
    cout.precision(0);

    RDouble cellCenter[3], faceCenter[3], faceNormal[3];

    //int interval = GetProgressInterval(nTotal, 10);
    int interval = GetProgressInterval(nTotalCell, 10);
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (viscousType <= LAMINAR && wallCellLabel[iCell] != 1 && iLES == NOLES_SOLVER) continue;

        if (iCell > 0 && iCell % interval == 0)
        {
            //ProgressMonitorToWindows(iCell, nTotal, "wall distance");
            //ProgressMonitorToLogFile(iCell, nTotal, "wall distance");
            ProgressMonitorToWindows(iCell, nTotalCell, "wall distance");
            ProgressMonitorToLogFile(iCell, nTotalCell, "wall distance");

            if (rotateAxis >= 3 || timeSpan.GetTimeSpan() > 3600.0)
            {
                ostringstream oss;
                oss.unsetf(ios::scientific);
                oss.setf(ios::fixed);
                oss.precision(0);
                oss << "     rotateAxis = " << this->rotateAxis << ", " << floor((iCell * 1.0) / nTotalCell * 100 + 0.5) << "% " << "wall distance" << " completed ...";
                WriteLogFile(oss, false);
            }
        }

        if (walldist[iCell] < LARGE * 0.9)
        {
            continue;
        }

        cellCenter[0] = xcc[iCell];
        cellCenter[1] = ycc[iCell];
        cellCenter[2] = zcc[iCell];

        for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
        {
            WallStructure *cwalldist = (*wallstructureList)[iWall];

            uint_t nWallFaces = cwalldist->GetNumberOfWallFaces();
            WallStructure::value_type &xfc = cwalldist->GetXFaceCenter();
            WallStructure::value_type &yfc = cwalldist->GetYFaceCenter();
            WallStructure::value_type &zfc = cwalldist->GetZFaceCenter();
            WallStructure::value_type &x = cwalldist->GetX();
            WallStructure::value_type &y = cwalldist->GetY();
            WallStructure::value_type &z = cwalldist->GetZ();
            WallStructure::value_type &xfn = cwalldist->GetXFaceNormal();
            WallStructure::value_type &yfn = cwalldist->GetYFaceNormal();
            WallStructure::value_type &zfn = cwalldist->GetZFaceNormal();

            int *nPointPerFace = cwalldist->GetnPointPerFace(); 
            int *wallFace2Node = cwalldist->GetWallFace2Node();

            int nodeCount = 0;
            int markFace = -1;
            for (int iFace = 0; iFace < nWallFaces; ++ iFace)
            {
                faceCenter[0] = xfc[iFace];
                faceCenter[1] = yfc[iFace];
                faceCenter[2] = zfc[iFace];

                faceNormal[0] = xfn[iFace];
                faceNormal[1] = yfn[iFace];
                faceNormal[2] = zfn[iFace];

                vector<int> face2node;
                face2node.resize(nPointPerFace[iFace]);
                for (int iNode = 0; iNode < nPointPerFace[iFace]; ++ iNode)
                {
                    face2node[iNode] = wallFace2Node[nodeCount ++];
                }

                RDouble wdst = ComputeDistanceOfANodeToAFace(cellCenter, faceCenter, faceNormal, face2node, x, y, z);
                if (walldist[iCell] > wdst)
                {
                    walldist[iCell] = wdst;
                    markFace = iFace;
                }
            }

            if (markFace <= -1)
            {
                continue;
            }

            nearestwallfacenormalx[iCell] = xfn[markFace];
            nearestwallfacenormaly[iCell] = yfn[markFace];
            nearestwallfacenormalz[iCell] = zfn[markFace];

        }
    }
}

void Pre_WalldistCompute::ComputeWallDist3DByProjectingFastAndNodeMethod(UnstructGrid *grid)
{
    using namespace PHMPI;

    RDouble *wallDistNode = grid->GetWallDistNode();
    int nTotalNode = grid->GetNTotalNode();
    RDouble *nearestwallfacenormalx = grid->GetNearestWallFaceNormalX();
    RDouble *nearestwallfacenormaly = grid->GetNearestWallFaceNormalY();
    RDouble *nearestwallfacenormalz = grid->GetNearestWallFaceNormalZ();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble oldNodeCoord[3], oldFaceCenter[3], oldFaceNormal[3];
    RDouble nodeCoord[3], faceCenter[3], faceNormal[3];

    RDouble *rmin, *rmax;
    rmin = wallFaceTree->GetMin();
    rmax = wallFaceTree->GetMax();
    RDouble pmin[6], pmax[6];
    RDouble dis = 0.1 * wallTolerance;
    RDouble globalGridSize = GetGlobalGridSize();

    PrintToWindow("    Number of total wall faces: ", nTWallFace, "\n");

    //! The maximum wanted wall faces that in the box,
    //! it is set to be sqrt(nTWallFace), usually.
    uint_t nBoxLimitedMax = MAX(1 * (int)sqrt(nTWallFace * 1.0), 1);
    nBoxLimitedMax = MIN(nBoxLimitedMax, nTWallFace);
    uint_t fixed_nBoxLimitedMax = 1000;
    nBoxLimitedMax = MIN(nBoxLimitedMax, fixed_nBoxLimitedMax);
    int nBoxLimitedMin = MAX((int)(nBoxLimitedMax * 0.1), 0);

    PrintToWindow("    Min & Max wall faces in box: ", nBoxLimitedMin, nBoxLimitedMax, "\n");

    TimeSpan partTime;
    int interval = GetProgressInterval(nTotalNode, 10);
    for (int iNode = 0; iNode < nTotalNode; ++iNode)
    {
        if (iNode > 0 && iNode % interval == 0)
        {
            ProgressMonitorToWindows(iNode, nTotalNode, "wall distance");
            ProgressMonitorToLogFile(iNode, nTotalNode, "wall distance");

            if (rotateAxis >= 1 || timeSpan.GetTimeSpan() > 600.0)
            {
                ostringstream oss;
                oss.unsetf(ios::scientific);
                oss.setf(ios::fixed);
                oss.precision(0);
                oss << "     rotateAxis = " << this->rotateAxis << ", " << floor((iNode * 1.0) / nTotalNode * 100 + 0.5) << "% " << "wall distance" << " completed ...";
                WriteLogFile(oss, true);
            }
        }

        if (wallDistNode[iNode] < LARGE * 0.9)
        {
            //! This cell wall distance has been computed.
            continue;
        }

        bool finalFound = false;
        int boxMin = nBoxLimitedMin;
        uint_t boxMax = nBoxLimitedMax;

        WallFaceTree::AdtNodeList nList;

        //! Enlarge time, the number of probe time.
        //! If lose once, enlarge the search scope then, by multiple 2 times.
        int nEnlargeTime = 3;
        int iEnlarge = 0;
        while (iEnlarge < nEnlargeTime)
        {
            ++ iEnlarge;

            oldNodeCoord[0] = x[iNode];
            oldNodeCoord[1] = y[iNode];
            oldNodeCoord[2] = z[iNode];
            this->RotateAxis(oldNodeCoord, nodeCoord);

            RDouble x0 = nodeCoord[0];
            RDouble y0 = nodeCoord[1];
            RDouble z0 = nodeCoord[2];

            pmin[0] = rmin[0];
            pmin[1] = rmin[1];
            pmin[2] = rmin[2];
            pmax[3] = rmax[3];
            pmax[4] = rmax[4];
            pmax[5] = rmax[5];

            pmin[3] = x0 - dis;
            pmin[4] = y0 - dis;
            pmin[5] = z0 - dis;

            pmax[0] = x0 + dis;
            pmax[1] = y0 + dis;
            pmax[2] = z0 + dis;

            RDouble size_min = 0;
            RDouble size_max = dis * 2.0;
            RDouble ra;

            int iSearch = 0;
            bool isFoundNearestFace = false;
            int nSearchTimes = 100;
            while (iSearch < nSearchTimes)
            {
                ++ iSearch;

                if (size_min >= globalGridSize) break;

                ra = 0.5 * (size_min + size_max);
                pmin[3] = x0 - ra;
                pmin[4] = y0 - ra;
                pmin[5] = z0 - ra;

                pmax[0] = x0 + ra;
                pmax[1] = y0 + ra;
                pmax[2] = z0 + ra;

                nList.clear();
                nList.reserve(boxMax + 1);
                wallFaceTree->FindNodesInRegion(pmin, pmax, nList, boxMax + 1);

                uint_t nNearWallFace = nList.size();
                if ((nNearWallFace > boxMin && nNearWallFace <= boxMax))
                {
                    int markFace = -1;
                    for (int iWall = 0; iWall < nNearWallFace; ++iWall)
                    {
                        uint_t wallStructID = nList[iWall]->GetData().first;
                        int faceID = nList[iWall]->GetData().second;
                        WallStructure *cwalldist = (*wallstructureList)[wallStructID];

                        vector<int> face2node;
                        cwalldist->GetFace2NodeVector(face2node, faceID);
                        WallStructure::value_type &xfc = cwalldist->GetXFaceCenter();
                        WallStructure::value_type &yfc = cwalldist->GetYFaceCenter();
                        WallStructure::value_type &zfc = cwalldist->GetZFaceCenter();
                        WallStructure::value_type &xW   = cwalldist->GetX();
                        WallStructure::value_type &yW   = cwalldist->GetY();
                        WallStructure::value_type &zW   = cwalldist->GetZ();
                        WallStructure::value_type &xfn = cwalldist->GetXFaceNormal();
                        WallStructure::value_type &yfn = cwalldist->GetYFaceNormal();
                        WallStructure::value_type &zfn = cwalldist->GetZFaceNormal();

                        oldFaceCenter[0] = xfc[faceID];
                        oldFaceCenter[1] = yfc[faceID];
                        oldFaceCenter[2] = zfc[faceID];
                        this->RotateAxis(oldFaceCenter, faceCenter);

                        oldFaceNormal[0] = xfn[faceID];
                        oldFaceNormal[1] = yfn[faceID];
                        oldFaceNormal[2] = zfn[faceID];
                        this->RotateAxis(oldFaceNormal, faceNormal);

                        RDouble wdst;
                        wdst = ComputeDistanceOfANodeToAFace(nodeCoord, faceCenter, faceNormal, face2node, xW, yW, zW);
                        if (wallDistNode[iNode] > wdst && wdst < 1.0e30)
                        {
                            wallDistNode[iNode] = wdst;
                            isFoundNearestFace = true;
                            markFace = iWall;
                        }
                    }
                    uint_t wallStructID = nList[markFace]->GetData().first;
                    int faceID = nList[markFace]->GetData().second;
                    WallStructure* cwalldist = (*wallstructureList)[wallStructID];

                    WallStructure::value_type& xfn = cwalldist->GetXFaceNormal();
                    WallStructure::value_type& yfn = cwalldist->GetYFaceNormal();
                    WallStructure::value_type& zfn = cwalldist->GetZFaceNormal();

                    nearestwallfacenormalx[iNode] = xfn[faceID];
                    nearestwallfacenormaly[iNode] = yfn[faceID];
                    nearestwallfacenormalz[iNode] = zfn[faceID];

                    break;
                }
                else if (nNearWallFace > boxMax && (size_max - size_min) > 0.01 * dis)
                {
                    //! The box is too large and so lead to too large number of wall faces,
                    //! in order to reduce the wall faces, decrease the search radius.
                    size_max = ra;
                    continue;
                }
                else
                {
                    //! The box is too small and so lead to nearly non wall faces,
                    //! in order to get more the wall faces, increase the search radius.
                    size_min = ra;
                    size_max = size_max * 2.0;
                    continue;
                }
            }

            if (isFoundNearestFace)
            {
                finalFound = true;
                break;
            }
            else
            {
                //! Enlarge the limitation.
                boxMax *= 2;
                boxMin = MAX(1 * (int)sqrt(boxMax * 1.0), 1);
                continue;
            }
        }
    }
}

void Pre_WalldistCompute::ComputeWallDist3DByProjectingFast(UnstructGrid *grid)
{
    using namespace PHMPI;

    RDouble *walldist = grid->GetWallDist();
    RDouble *nearestwallfacenormalx = grid->GetNearestWallFaceNormalX();
    RDouble *nearestwallfacenormaly = grid->GetNearestWallFaceNormalY();
    RDouble *nearestwallfacenormalz = grid->GetNearestWallFaceNormalZ();

    int nTotalCell = grid->GetNTotalCell();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
    int *wallCellLabel = grid->GetWallCellLabel();

    cout.unsetf(ios::scientific);
    cout.setf(ios::fixed);
    cout.precision(0);

    RDouble oldCellCenter[3], oldFaceCenter[3], oldFaceNormal[3];
    RDouble cellCenter[3], faceCenter[3], faceNormal[3];

    RDouble *rmin, *rmax;
    rmin = wallFaceTree->GetMin();
    rmax = wallFaceTree->GetMax();
    RDouble pmin[6], pmax[6];
    RDouble dis = 0.1 * wallTolerance;
    RDouble globalGridSize = GetGlobalGridSize();

    PrintToWindow("    Number of total wall faces: ", nTWallFace, "\n");

    //! The maximum wanted wall faces that in the box,
    //! it is set to be sqrt(nTWallFace), usually.
    uint_t nBoxLimitedMax = MAX(1 * (int)sqrt(nTWallFace * 1.0), 1);
    nBoxLimitedMax = MIN(nBoxLimitedMax, nTWallFace);
    uint_t fixed_nBoxLimitedMax = 1000;
    nBoxLimitedMax = MIN(nBoxLimitedMax, fixed_nBoxLimitedMax);
    int nBoxLimitedMin = MAX((int)(nBoxLimitedMax * 0.1), 0);

    PrintToWindow("    Min & Max wall faces in box: ", nBoxLimitedMin, nBoxLimitedMax, "\n");

    TimeSpan partTime;
    int interval = GetProgressInterval(nTotalCell, 10);

#ifdef AI_RandomForestRegressor
    vector<RDouble> icellToWallNormal;
    vector<vector<RDouble>> allCellToWallNormal;
#endif

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (viscousType <= LAMINAR && wallCellLabel[iCell] != 1 && iLES == NOLES_SOLVER) continue;

        if (iCell > 0 && iCell % interval == 0)
        {
            ProgressMonitorToWindows(iCell, nTotalCell, "wall distance");
            ProgressMonitorToLogFile(iCell, nTotalCell, "wall distance");

            if (rotateAxis >= 1 || timeSpan.GetTimeSpan() > 600.0)
            {
                ostringstream oss;
                oss.unsetf(ios::scientific);
                oss.setf(ios::fixed);
                oss.precision(0);
                oss << "     rotateAxis = " << this->rotateAxis << ", " << floor((iCell * 1.0) / nTotalCell * 100 + 0.5) << "% " << "wall distance" << " completed ...";
                WriteLogFile(oss, false);
            }
        }

        if (walldist[iCell] < LARGE * 0.9)
        {
            //! This cell wall distance has been computed.
            continue;
        }

        bool finalFound = false;
        int boxMin = nBoxLimitedMin;
        uint_t boxMax = nBoxLimitedMax;

        WallFaceTree::AdtNodeList nList;

        //! Enlarge time, the number of probe time.
        //! If lose once, enlarge the search scope then, by multiple 2 times.
        int nEnlargeTime = 3;
        int iEnlarge = 0;
        while (iEnlarge < nEnlargeTime)
        {
            ++ iEnlarge;

            //! i-th test.
            oldCellCenter[0] = xcc[iCell];
            oldCellCenter[1] = ycc[iCell];
            oldCellCenter[2] = zcc[iCell];
            this->RotateAxis(oldCellCenter, cellCenter);

            //! For method 0: loop over the nearest several box.
            RDouble x0 = cellCenter[0];
            RDouble y0 = cellCenter[1];
            RDouble z0 = cellCenter[2];

            pmin[0] = rmin[0];
            pmin[1] = rmin[1];
            pmin[2] = rmin[2];
            pmax[3] = rmax[3];
            pmax[4] = rmax[4];
            pmax[5] = rmax[5];

            pmin[3] = x0 - dis;
            pmin[4] = y0 - dis;
            pmin[5] = z0 - dis;

            pmax[0] = x0 + dis;
            pmax[1] = y0 + dis;
            pmax[2] = z0 + dis;

            RDouble size_min = 0;
            RDouble size_max = dis * 2.0;
            RDouble ra;

            int iSearch = 0;
            bool isFoundNearestFace = false;
            int nSearchTimes = 100;
            while (iSearch < nSearchTimes)
            {
                ++ iSearch;

                if (size_min >= globalGridSize) break;

                ra = 0.5 * (size_min + size_max);
                pmin[3] = x0 - ra;
                pmin[4] = y0 - ra;
                pmin[5] = z0 - ra;

                pmax[0] = x0 + ra;
                pmax[1] = y0 + ra;
                pmax[2] = z0 + ra;

                nList.clear();
                nList.reserve(boxMax + 1);
                wallFaceTree->FindNodesInRegion(pmin, pmax, nList, boxMax + 1);
                
                uint_t nNearWallFace = nList.size();
                if ((nNearWallFace > boxMin && nNearWallFace <= boxMax))
                // || (iEnlarge == nEnlargeTime && nNearWallFace > 0 && nNearWallFace < boxMax * 5 && iSearch >= 0.8 * nSearchTimes))
                {
                    //! Case 1: the number of wall faces in the box is according with the anticipate.
                    //! Case 2: the worst case, to ensure the success, accept this case, although 
                    //!         the wall faces found are to much.
                    int markFace = -1;
                    for (int iWall = 0; iWall < nNearWallFace; ++ iWall)
                    {
                        uint_t wallStructID = nList[iWall]->GetData().first;
                        int faceID = nList[iWall]->GetData().second;
                        WallStructure *cwalldist = (*wallstructureList)[wallStructID];

                        vector<int> face2node;
                        cwalldist->GetFace2NodeVector(face2node, faceID);
                        WallStructure::value_type &xfc = cwalldist->GetXFaceCenter();
                        WallStructure::value_type &yfc = cwalldist->GetYFaceCenter();
                        WallStructure::value_type &zfc = cwalldist->GetZFaceCenter();
                        WallStructure::value_type &x = cwalldist->GetX();
                        WallStructure::value_type &y = cwalldist->GetY();
                        WallStructure::value_type &z = cwalldist->GetZ();
                        WallStructure::value_type &xfn = cwalldist->GetXFaceNormal();
                        WallStructure::value_type &yfn = cwalldist->GetYFaceNormal();
                        WallStructure::value_type &zfn = cwalldist->GetZFaceNormal();

                        oldFaceCenter[0] = xfc[faceID];
                        oldFaceCenter[1] = yfc[faceID];
                        oldFaceCenter[2] = zfc[faceID];
                        this->RotateAxis(oldFaceCenter, faceCenter);

                        oldFaceNormal[0] = xfn[faceID];
                        oldFaceNormal[1] = yfn[faceID];
                        oldFaceNormal[2] = zfn[faceID];
                        this->RotateAxis(oldFaceNormal, faceNormal);

                        RDouble wdst;
                        wdst = ComputeDistanceOfANodeToAFace(cellCenter, faceCenter, faceNormal, face2node, x, y, z);

                        //RDouble distTemp = sqrt(wdst);

                        if (walldist[iCell] > wdst && wdst < 1.0e30)
                        {
                            walldist[iCell] = wdst;
                            isFoundNearestFace = true;
                            markFace = iWall;

#ifdef AI_RandomForestRegressor
                            icellToWallNormal[0] = faceNormal[0];
                            icellToWallNormal[1] = faceNormal[1];
                            icellToWallNormal[2] = faceNormal[2];
#endif
                        }
                    }
#ifdef AI_RandomForestRegressor
                    allCellToWallNormal.push_back(icellToWallNormal);
#endif
                    uint_t wallStructID = nList[markFace]->GetData().first;
                    int faceID = nList[markFace]->GetData().second;
                    WallStructure* cwalldist = (*wallstructureList)[wallStructID];

                    WallStructure::value_type& xfn = cwalldist->GetXFaceNormal();
                    WallStructure::value_type& yfn = cwalldist->GetYFaceNormal();
                    WallStructure::value_type& zfn = cwalldist->GetZFaceNormal();

                    nearestwallfacenormalx[iCell] = xfn[faceID];
                    nearestwallfacenormaly[iCell] = yfn[faceID];
                    nearestwallfacenormalz[iCell] = zfn[faceID];

                    break;
                } 
                else if (nNearWallFace > boxMax && (size_max - size_min) > 0.01 * dis)
                {
                    //! The box is too large and so lead to too large number of wall faces,
                    //! in order to reduce the wall faces, decrease the search radius.
                    size_max = ra;
                    continue;
                }
                else
                {
                    //! The box is too small and so lead to nearly non wall faces,
                    //! in order to get more the wall faces, increase the search radius.
                    size_min = ra;
                    size_max = size_max * 2.0;
                    continue;
                }
            }    //! iSearch while.

            if (isFoundNearestFace)
            {
                finalFound = true;
                break;
            }
            else
            {
                //! Enlarge the limitation.
                //boxMin = floor((boxMin + boxMax) * 0.5 + 0.5);
                boxMax *= 2;
                boxMin = MAX(1 * (int)sqrt(boxMax * 1.0), 1);
                continue;
            }
        }
    }    //! iCell.
#ifdef AI_RandomForestRegressor
    DumpCellToWallNormal(allCellToWallNormal);
#endif
}

#ifdef AI_RandomForestRegressor
void  DumpCellToWallNormal(vector<vector<RDouble>> allCellToWallNormal)
{
    ofstream dataFile;
    dataFile.open("./grid/CellToWallNormal.txt", ofstream::app);
    for(int i = 0; i<allCellToWallNormal.size(); ++i)
    {
        dataFile << allCellToWallNormal[i][0] << " " << allCellToWallNormal[i][1] << " " << allCellToWallNormal[i][2] << endl;
    }
    dataFile.close();
}
#endif

void Pre_WalldistCompute::ComputeWallDist3DByProjectingFast(StructGrid *grid)
{
    using namespace PHMPI;

    RDouble3D &walldist = *grid->GetWallDist();
    RDouble3D &nearestwallfacenormalx = *grid->GetNearestWallFaceNormalX();
    RDouble3D &nearestwallfacenormaly = *grid->GetNearestWallFaceNormalY();
    RDouble3D &nearestwallfacenormalz = *grid->GetNearestWallFaceNormalZ();

    int nTotalCell = grid->GetNTotalCell();

    RDouble3D &xcc = *grid->GetCellCenterX();
    RDouble3D &ycc = *grid->GetCellCenterY();
    RDouble3D &zcc = *grid->GetCellCenterZ();

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");

    cout.unsetf(ios::scientific);
    cout.setf(ios::fixed);
    cout.precision(0);

    RDouble oldCellCenter[3], oldFaceCenter[3], oldFaceNormal[3];
    RDouble cellCenter[3], faceCenter[3], faceNormal[3];

    RDouble *rmin, *rmax;
    rmin = wallFaceTree->GetMin();
    rmax = wallFaceTree->GetMax();
    RDouble pmin[6], pmax[6];
    RDouble dis = 0.1 * wallTolerance;
    RDouble globalGridSize = GetGlobalGridSize();

    PrintToWindow("    Number of total wall faces: ", nTWallFace, "\n");

    //! The maximum wanted wall faces that in the box,
    //! it is set to be sqrt(nTWallFace), usually.
    uint_t nBoxLimitedMax = MAX(1 * (int)sqrt(nTWallFace * 1.0), 1);
    nBoxLimitedMax = MIN(nBoxLimitedMax, nTWallFace);
    uint_t fixed_nBoxLimitedMax = 1000;
    nBoxLimitedMax = MIN(nBoxLimitedMax, fixed_nBoxLimitedMax);
    int nBoxLimitedMin = MAX((int)(nBoxLimitedMax * 0.1), 0);

    PrintToWindow("    Min & Max wall faces in box: ", nBoxLimitedMin, nBoxLimitedMax, "\n");

    TimeSpan partTime;
    int interval = GetProgressInterval(nTotalCell, 10);
    int count = 0;

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                ++count;
                if (viscousType <= LAMINAR && iLES == NOLES_SOLVER)
                {
                    continue;
                }

                if (count > 1 && count % interval == 0)
                {
                    ProgressMonitorToWindows(count, nTotalCell, "wall distance");
                    ProgressMonitorToLogFile(count, nTotalCell, "wall distance");
                }

                if (walldist(i, j, k) < LARGE * 0.9)
                {
                    //! This cell wall distance has been computed.
                    continue;
                }

                bool finalFound = false;
                int boxMin = nBoxLimitedMin;
                uint_t boxMax = nBoxLimitedMax;

                WallFaceTree::AdtNodeList nList;

                //! Enlarge time, the number of probe time.
                //! If lose once, enlarge the search scope then, by multiple 2 times.
                int nEnlargeTime = 3;
                int iEnlarge = 0;
                while (iEnlarge < nEnlargeTime)
                {
                    ++ iEnlarge;

                    //! i-th test.
                    oldCellCenter[0] = xcc(i, j, k);
                    oldCellCenter[1] = ycc(i, j, k);
                    oldCellCenter[2] = zcc(i, j, k);
                    this->RotateAxis(oldCellCenter, cellCenter);

                    //! For method 0: loop over the nearest several box.
                    RDouble x0 = cellCenter[0];
                    RDouble y0 = cellCenter[1];
                    RDouble z0 = cellCenter[2];

                    pmin[0] = rmin[0];
                    pmin[1] = rmin[1];
                    pmin[2] = rmin[2];
                    pmax[3] = rmax[3];
                    pmax[4] = rmax[4];
                    pmax[5] = rmax[5];

                    pmin[3] = x0 - dis;
                    pmin[4] = y0 - dis;
                    pmin[5] = z0 - dis;

                    pmax[0] = x0 + dis;
                    pmax[1] = y0 + dis;
                    pmax[2] = z0 + dis;

                    RDouble size_min = 0;
                    RDouble size_max = dis * 2.0;
                    RDouble ra;

                    int iSearch = 0;
                    bool isFoundNearestFace = false;
                    int nSearchTimes = 100;
                    while (iSearch < nSearchTimes)
                    {
                        ++ iSearch;
                        
                        if (size_min >= globalGridSize) break;
                        
                        ra = 0.5 * (size_min + size_max);
                        pmin[3] = x0 - ra;
                        pmin[4] = y0 - ra;
                        pmin[5] = z0 - ra;
                        
                        pmax[0] = x0 + ra;
                        pmax[1] = y0 + ra;
                        pmax[2] = z0 + ra;

                        nList.clear();
                        nList.reserve(boxMax + 1);
                        wallFaceTree->FindNodesInRegion(pmin, pmax, nList, boxMax + 1);

                        uint_t nNearWallFace = nList.size();
                        if ((nNearWallFace > boxMin && nNearWallFace <= boxMax))
                            // || (iEnlarge == nEnlargeTime && nNearWallFace > 0 && nNearWallFace < boxMax * 5 && iSearch >= 0.8 * nSearchTimes))
                        {
                            //! Case 1: the number of wall faces in the box is according with the anticipate.
                            //! Case 2: the worst case, to ensure the success, accept this case, although 
                            //!         the wall faces found are to much.
                            int markFace = -1;
                            for (int iWall = 0; iWall < nNearWallFace; ++ iWall)
                            {
                                uint_t wallStructID = nList[iWall]->GetData().first;
                                int faceID = nList[iWall]->GetData().second;
                                WallStructure *cwalldist = (*wallstructureList)[wallStructID];

                                vector<int> face2node;
                                cwalldist->GetFace2NodeVector(face2node, faceID);
                                WallStructure::value_type &xfc = cwalldist->GetXFaceCenter();
                                WallStructure::value_type &yfc = cwalldist->GetYFaceCenter();
                                WallStructure::value_type &zfc = cwalldist->GetZFaceCenter();
                                WallStructure::value_type &x = cwalldist->GetX();
                                WallStructure::value_type &y = cwalldist->GetY();
                                WallStructure::value_type &z = cwalldist->GetZ();
                                WallStructure::value_type &xfn = cwalldist->GetXFaceNormal();
                                WallStructure::value_type &yfn = cwalldist->GetYFaceNormal();
                                WallStructure::value_type &zfn = cwalldist->GetZFaceNormal();

                                oldFaceCenter[0] = xfc[faceID];
                                oldFaceCenter[1] = yfc[faceID];
                                oldFaceCenter[2] = zfc[faceID];
                                this->RotateAxis(oldFaceCenter, faceCenter);

                                oldFaceNormal[0] = -xfn[faceID];
                                oldFaceNormal[1] = -yfn[faceID];
                                oldFaceNormal[2] = -zfn[faceID];
                                this->RotateAxis(oldFaceNormal, faceNormal);

                                RDouble wdst;
                                wdst = ComputeDistanceOfANodeToAFace(cellCenter, faceCenter, faceNormal, face2node, x, y, z);

                                if (walldist(i, j, k) > wdst && wdst < 1.0e30)
                                {
                                    walldist(i, j, k) = wdst;
                                    isFoundNearestFace = true;
                                    markFace = iWall;
                                }
                            }

                            uint_t wallStructID = nList[markFace]->GetData().first;
                            int faceID = nList[markFace]->GetData().second;
                            WallStructure* cwalldist = (*wallstructureList)[wallStructID];

                            WallStructure::value_type& xfn = cwalldist->GetXFaceNormal();
                            WallStructure::value_type& yfn = cwalldist->GetYFaceNormal();
                            WallStructure::value_type& zfn = cwalldist->GetZFaceNormal();

                            nearestwallfacenormalx(i, j, k) = xfn[faceID];
                            nearestwallfacenormaly(i, j, k) = yfn[faceID];
                            nearestwallfacenormalz(i, j, k) = zfn[faceID];

                            break;
                        }
                        else if (nNearWallFace > boxMax && (size_max - size_min) > 0.01 * dis)
                        {
                            //! The box is too large and so lead to too large number of wall faces,
                            //! in order to reduce the wall faces, decrease the search radius.
                            size_max = ra;
                            continue;
                        }
                        else
                        {
                            //! The box is too small and so lead to nearly non wall faces,
                            //! in order to get more the wall faces, increase the search radius.
                            size_min = ra;
                            size_max = size_max * 2.0;
                            continue;
                        }
                    }    //! iSearch while.

                    if (isFoundNearestFace)
                    {
                        finalFound = true;
                        break;
                    }
                    else
                    {
                        //! Enlarge the limitation.
                        //boxMin = floor((boxMin + boxMax) * 0.5 + 0.5);
                        boxMax *= 2;
                        boxMin = MAX(1 * (int)sqrt(boxMax * 1.0), 1);
                        continue;
                    }
                }
            }
        }
    }    //! iCell.
}

void Pre_WalldistCompute::ComputeNodeWallDist3DByProjectingFast(int nTotalNode, RDouble *x, RDouble *y, RDouble *z, RDouble *nodeWalldist)
{
    using namespace PHMPI;

    //RDouble3D &nearestwallfacenormalx = *grid->GetNearestWallFaceNormalX();
    //RDouble3D &nearestwallfacenormaly = *grid->GetNearestWallFaceNormalY();
    //RDouble3D &nearestwallfacenormalz = *grid->GetNearestWallFaceNormalZ();

    cout.unsetf(ios::scientific);
    cout.setf(ios::fixed);
    cout.precision(0);

    RDouble oldCellCenter[3], oldFaceCenter[3], oldFaceNormal[3];
    RDouble cellCenter[3], faceCenter[3], faceNormal[3];

    RDouble *rmin, *rmax;
    rmin = wallFaceTree->GetMin();
    rmax = wallFaceTree->GetMax();
    RDouble pmin[6], pmax[6];
    RDouble dis = 0.1 * wallTolerance;
    RDouble globalGridSize = GetGlobalGridSize();

    //! The maximum wanted wall faces that in the box,
    //! it is set to be sqrt(nTWallFace), usually.
    uint_t nBoxLimitedMax = MAX(1 * (int)sqrt(nTWallFace * 1.0), 1);
    nBoxLimitedMax = MIN(nBoxLimitedMax, nTWallFace);
    uint_t fixed_nBoxLimitedMax = 1000;
    nBoxLimitedMax = MIN(nBoxLimitedMax, fixed_nBoxLimitedMax);
    int nBoxLimitedMin = 0;

    PrintToWindow("    Min & Max wall faces in box: ", nBoxLimitedMin, nBoxLimitedMax, "\n");

    int interval = GetProgressInterval(nTotalNode, 10);

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        int nodeIndex = iNode + 1;
        if (nodeIndex % interval == 0 || nodeIndex == nTotalNode)
        {
            ProgressMonitorToWindows(nodeIndex, nTotalNode, "wall distance");
            ProgressMonitorToLogFile(nodeIndex, nTotalNode, "wall distance");
        }

        if (nodeWalldist[iNode] < LARGE * 0.9)
        {
            //! This cell wall distance has been computed.
            continue;
        }

        bool finalFound = false;
        int boxMin = nBoxLimitedMin;
        uint_t boxMax = nBoxLimitedMax;

        WallFaceTree::AdtNodeList nList;

        //! Enlarge time, the number of probe time.
        //! If lose once, enlarge the search scope then, by multiple 2 times.
        int nEnlargeTime = 3;
        int iEnlarge = 0;
        while (iEnlarge < nEnlargeTime)
        {
            ++ iEnlarge;

            //! i-th test.
            oldCellCenter[0] = x[iNode];
            oldCellCenter[1] = y[iNode];
            oldCellCenter[2] = z[iNode];
            this->RotateAxis(oldCellCenter, cellCenter);

            //! For method 0: loop over the nearest several box.
            RDouble x0 = cellCenter[0];
            RDouble y0 = cellCenter[1];
            RDouble z0 = cellCenter[2];

            pmin[0] = rmin[0];
            pmin[1] = rmin[1];
            pmin[2] = rmin[2];
            pmax[3] = rmax[3];
            pmax[4] = rmax[4];
            pmax[5] = rmax[5];

            pmin[3] = x0 - dis;
            pmin[4] = y0 - dis;
            pmin[5] = z0 - dis;

            pmax[0] = x0 + dis;
            pmax[1] = y0 + dis;
            pmax[2] = z0 + dis;

            RDouble size_min = 0;
            RDouble size_max = dis * 2.0;
            RDouble ra;

            int iSearch = 0;
            bool isFoundNearestFace = false;
            int nSearchTimes = 100;
            while (iSearch < nSearchTimes)
            {
                ++ iSearch;

                if (size_min >= globalGridSize) break;

                ra = 0.5 * (size_min + size_max);
                pmin[3] = x0 - ra;
                pmin[4] = y0 - ra;
                pmin[5] = z0 - ra;

                pmax[0] = x0 + ra;
                pmax[1] = y0 + ra;
                pmax[2] = z0 + ra;

                nList.clear();
                nList.reserve(boxMax + 1);
                wallFaceTree->FindNodesInRegion(pmin, pmax, nList, boxMax + 1);

                uint_t nNearWallFace = nList.size();
                if ((nNearWallFace > boxMin && nNearWallFace <= boxMax))
                // || (iEnlarge == nEnlargeTime && nNearWallFace > 0 && nNearWallFace < boxMax * 5 && iSearch >= 0.8 * nSearchTimes))
                {
                    //! Case 1: the number of wall faces in the box is according with the anticipate.
                    //! Case 2: the worst case, to ensure the success, accept this case, although 
                    //!         the wall faces found are to much.
                    //int markFace = -1;
                    for (uint_t iWall = 0; iWall < nNearWallFace; ++ iWall)
                    {
                        uint_t wallStructID = nList[iWall]->GetData().first;
                        int faceID = nList[iWall]->GetData().second;
                        WallStructure *cwalldist = (*wallstructureList)[wallStructID];

                        vector<int> face2node;
                        cwalldist->GetFace2NodeVector(face2node, faceID);
                        WallStructure::value_type &xfc = cwalldist->GetXFaceCenter();
                        WallStructure::value_type &yfc = cwalldist->GetYFaceCenter();
                        WallStructure::value_type &zfc = cwalldist->GetZFaceCenter();
                        WallStructure::value_type &xWD = cwalldist->GetX();
                        WallStructure::value_type &yWD = cwalldist->GetY();
                        WallStructure::value_type &zWD = cwalldist->GetZ();
                        WallStructure::value_type &xfn = cwalldist->GetXFaceNormal();
                        WallStructure::value_type &yfn = cwalldist->GetYFaceNormal();
                        WallStructure::value_type &zfn = cwalldist->GetZFaceNormal();

                        oldFaceCenter[0] = xfc[faceID];
                        oldFaceCenter[1] = yfc[faceID];
                        oldFaceCenter[2] = zfc[faceID];
                        this->RotateAxis(oldFaceCenter, faceCenter);

                        oldFaceNormal[0] = xfn[faceID];
                        oldFaceNormal[1] = yfn[faceID];
                        oldFaceNormal[2] = zfn[faceID];
                        this->RotateAxis(oldFaceNormal, faceNormal);

                        RDouble wdst;
                        wdst = ComputeDistanceOfANodeToAFace(cellCenter, faceCenter, faceNormal, face2node, xWD, yWD, zWD);

                        wdst = sqrt(wdst);

                        if (nodeWalldist[iNode] > wdst && wdst < 1.0e30)
                        {
                            nodeWalldist[iNode] = wdst;
                            isFoundNearestFace = true;
                            //markFace = iWall;
                        }
                    }
                    //uint_t wallStructID = nList[markFace]->GetData().first;
                    //int faceID = nList[markFace]->GetData().second;
                    //WallStructure* cwalldist = (*wallstructureList)[wallStructID];

                    //WallStructure::value_type& xfn = cwalldist->GetXFaceNormal();
                    //WallStructure::value_type& yfn = cwalldist->GetYFaceNormal();
                    //WallStructure::value_type& zfn = cwalldist->GetZFaceNormal();

                    //nearestwallfacenormalx[iNode] = xfn[faceID];
                    //nearestwallfacenormaly[iNode] = yfn[faceID];
                    //nearestwallfacenormalz[iNode] = zfn[faceID];

                    break;
                } 
                else if (nNearWallFace > boxMax && (size_max - size_min) > 0.01 * dis)
                {
                    //! The box is too large and so lead to too large number of wall faces,
                    //! in order to reduce the wall faces, decrease the search radius.
                    size_max = ra;
                    continue;
                }
                else
                {
                    //! The box is too small and so lead to nearly non wall faces,
                    //! in order to get more the wall faces, increase the search radius.
                    size_min = ra;
                    size_max = size_max * 2.0;
                    continue;
                }

            }    //! iSearch while.

            if (isFoundNearestFace)
            {
                finalFound = true;
                break;
            }
            else
            {
                //! Enlarge the limitation.
                boxMax *= 2;
                boxMin = MAX(1 * (int)sqrt(boxMax * 1.0), 1);
                continue;
            }
        }
    }
}

void Pre_WalldistCompute::ComputeNodeWallDist3DByKDTreeMethod(int nTotalNode, RDouble *x, RDouble *y, RDouble *z, RDouble *nodeWalldist)
{
    using namespace PHMPI;
    if (!wallNodeKDTree) return;
    RDouble *xyz = new RDouble[3];
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        xyz[0] = x[iNode];
        xyz[1] = y[iNode];
        xyz[2] = z[iNode];
        KDRes *minResults = KDNearest(wallNodeKDTree, xyz);
        RDouble wdst = sqrt(DistSQ(minResults->Position(), xyz, 3));

        if (nodeWalldist[iNode] > wdst && wdst < 1.0e30)
        {
            nodeWalldist[iNode] = wdst;
        }

        FreeKDRes(minResults);
     }
    delete [] xyz;    xyz = NULL;
}

void Pre_WalldistCompute::ComputeWallDistance(StructGrid *grid)
{
    RDouble3D &walldist = *grid->GetWallDist();
    RDouble3D &xcc = *grid->GetCellCenterX();
    RDouble3D &ycc = *grid->GetCellCenterY();
    RDouble3D &zcc = *grid->GetCellCenterZ();
    RDouble3D &nearestwallfacenormalx = *grid->GetNearestWallFaceNormalX();
    RDouble3D &nearestwallfacenormaly = *grid->GetNearestWallFaceNormalY();
    RDouble3D &nearestwallfacenormalz = *grid->GetNearestWallFaceNormalZ();

    RDouble cellCenter[3], faceCenter[3], faceNormal[3];

    int count = 0;

    int nTotalCell = grid->GetNTotalCell();

    cout.unsetf(ios::scientific);
    cout.setf(ios::fixed);
    cout.precision(0);
    if (myid == 0)
    {
        cout << "  Wall distance computing for zone " << grid->GetZoneID() << ": \n";
    }

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int interval = GetProgressInterval(nTotalCell, 10);
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                ++ count;
                if (count > 1 && count % interval == 0)
                {
                    ProgressMonitorToWindows(count, nTotalCell, "wall distance");
                    ProgressMonitorToLogFile(count, nTotalCell, "wall distance");
                }
                if (walldist(i, j, k) < LARGE * 0.9)
                {
                    continue;
                }

                cellCenter[0] = xcc(i, j, k);
                cellCenter[1] = ycc(i, j, k);
                cellCenter[2] = zcc(i, j, k);
                
                for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
                {
                    WallStructure *cwalldist = (*wallstructureList)[iWall];

                    uint_t nWallFaces = cwalldist->GetNumberOfWallFaces();

                    WallStructure::value_type &xfc = cwalldist->GetXFaceCenter();
                    WallStructure::value_type &yfc = cwalldist->GetYFaceCenter();
                    WallStructure::value_type &zfc = cwalldist->GetZFaceCenter();
                    WallStructure::value_type &x = cwalldist->GetX();
                    WallStructure::value_type &y = cwalldist->GetY();
                    WallStructure::value_type &z = cwalldist->GetZ();
                    WallStructure::value_type &xfn = cwalldist->GetXFaceNormal();
                    WallStructure::value_type &yfn = cwalldist->GetYFaceNormal();
                    WallStructure::value_type &zfn = cwalldist->GetZFaceNormal();

                    int *nPointPerFace = cwalldist->GetnPointPerFace();
                    int *wallFace2Node = cwalldist->GetWallFace2Node();
                    int nodeCount = 0;
                    vector <int> face2node;

                    int markFace = -1;
                    for (int iFace = 0; iFace < nWallFaces; ++ iFace)
                    {
                        faceCenter[0] = xfc[iFace];
                        faceCenter[1] = yfc[iFace];
                        faceCenter[2] = zfc[iFace];

                        faceNormal[0] = -xfn[iFace];
                        faceNormal[1] = -yfn[iFace];
                        faceNormal[2] = -zfn[iFace];

                        face2node.resize(nPointPerFace[iFace]);
                        for (int iNode = 0; iNode < nPointPerFace[iFace]; ++ iNode)
                        {
                            face2node[iNode] = wallFace2Node[nodeCount ++];
                        }

                        RDouble wdst = ComputeDistanceOfANodeToAFace(cellCenter, faceCenter, faceNormal, face2node, x, y, z);

                        if (walldist(i, j, k) > wdst)
                        {
                            walldist(i, j, k) = wdst;
                            markFace = static_cast<int>(iWall);
                        }
                    }

                    if (markFace <= -1)
                    {
                        continue;
                    }

                    nearestwallfacenormalx(i, j, k) = -xfn[markFace];
                    nearestwallfacenormaly(i, j, k) = -yfn[markFace];
                    nearestwallfacenormalz(i, j, k) = -zfn[markFace];

                }
            }
        }
    }
}

void Pre_WalldistCompute::BuildWallStructKDTree()
{ 
    int nDim = 3;

    wallFaceKDTree = CreatKDTree(nDim);

    if (wallFaceKDTree == NULL)
    {
        return;
    }

    for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
    {
        WallStructure *cwalldist = (*wallstructureList)[iWall];
        WallStructure::value_type &xfc = cwalldist->GetXFaceCenter();
        WallStructure::value_type &yfc = cwalldist->GetYFaceCenter();
        WallStructure::value_type &zfc = cwalldist->GetZFaceCenter();

        for (std::size_t m = 0; m < cwalldist->GetNumberOfWallFaces(); ++ m)
        {
            RDouble *xyzWall = new RDouble [nDim];
            xyzWall[0] = xfc[m];
            xyzWall[1] = yfc[m];
            xyzWall[2] = zfc[m];

            KDInsert(wallFaceKDTree, xyzWall, NULL);
            delete []xyzWall;
        }
    }
}

void Pre_WalldistCompute::BuildWallStructKDTreeNode()
{ 
    int nDim = 3;

    nTWallFace = 0;
    for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
    {
        WallStructure *cwalldist = (*wallstructureList)[iWall];
        if (cwalldist->GetIBlock() != presentBlock)
        {
            continue;
        }

        uint_t nWallFaces = cwalldist->GetNumberOfWallFaces();
        uint_t nWallNodes = cwalldist->GetNumberOfWallPoints();
        nTWallFace += nWallFaces;

        WallStructure::value_type &x = cwalldist->GetX();
        WallStructure::value_type &y = cwalldist->GetY();
        WallStructure::value_type &z = cwalldist->GetZ();

        for (std::size_t iNode = 0; iNode < nWallNodes; ++ iNode)
        {
            RDouble *xyzWall = new RDouble [nDim];
            xyzWall[0] = x[iNode];
            xyzWall[1] = y[iNode];
            if(nDim == 3)xyzWall[2] = z[iNode];

            KDInsert(wallNodeKDTree, xyzWall, NULL);
            delete [] xyzWall;
        }
    }
}

void Pre_WalldistCompute::BuildWallStructKDTreeNodeInOtherBlock()
{ 
    int nDim = 3;
    nTWallFace = 0;

    for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
    {
        WallStructure *cwalldist = (*wallstructureList)[iWall];
        if (cwalldist->GetIBlock() == presentBlock)
        {
            continue;
        }

        WallStructure::value_type &x = cwalldist->GetX();
        WallStructure::value_type &y = cwalldist->GetY();
        WallStructure::value_type &z = cwalldist->GetZ();

        uint_t nWallFaces = cwalldist->GetNumberOfWallFaces();
        nTWallFace += nWallFaces;

        uint_t nWallNodes = cwalldist->GetNumberOfWallPoints();

        for (std::size_t iNode = 0; iNode < nWallNodes; ++ iNode)
        {
            RDouble *xyzWall = new RDouble [nDim];
            xyzWall[0] = x[iNode];
            xyzWall[1] = y[iNode];
            if(nDim == 3)xyzWall[2] = z[iNode];

            KDInsert(wallNodeKDTree, xyzWall, NULL);
            delete [] xyzWall;
        }
    }
}

void Pre_WalldistCompute::ComputeWallDistanceKDTree(StructGrid *grid)
{
    RDouble3D &walldist = *grid->GetWallDist();
    RDouble *xyzCellCenter = new RDouble[3];

    int count = 0;

    int nTotalCell = grid->GetNTotalCell();

    cout.unsetf(ios::scientific);
    cout.setf(ios::fixed);
    cout.precision(0);
    if (myid == 0)

    {
        cout << "  Wall distance computing for zone " << grid->GetZoneID() << ": \n";
    }

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int interval = GetProgressInterval(nTotalCell, 10);
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                ++ count;
                if (count > 1 && count % interval == 0)
                {
                    ProgressMonitorToWindows(count, nTotalCell, "wall distance");
                    ProgressMonitorToLogFile(count, nTotalCell, "wall distance");
                }
                grid->CenterCoor(i, j, k, xyzCellCenter[0], xyzCellCenter[1], xyzCellCenter[2]);
                KDRes *minResults = KDNearest(wallFaceKDTree, xyzCellCenter);
                RDouble wdst = sqrt(DistSQ(minResults->Position(), xyzCellCenter, 3));
                walldist(i, j, k) = wdst;
            }
        }
    }
    delete [] xyzCellCenter;    xyzCellCenter = nullptr;
}

//! New method: communicate with each processor, after compress all zones of each processor.
void Pre_WalldistCompute::ServerCollection()
{
    //WriteLogFile("Wall structure server collecting ...");

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int serverTmp = GetServerSepMode();
    int numberOfOriginalprocessor = GetNumberOfProcessor();
    int numberOfGridprocessor = GetNumberOfGridProcessor();

    DataContainer *cdata = new DataContainer();
    cdata->MoveToBegin();

    //! Write the number of zones in this processor.
    int nZonesInThisProc = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int proc = GetZoneProcessorIDSepMode(iZone);
        if (myid == proc)
        {
            ++ nZonesInThisProc;
        }
    }
    PHWrite(cdata, nZonesInThisProc);
    
    //! Compress all of the wall structure in this processor into data-container.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int proc = GetZoneProcessorIDSepMode(iZone);

        if (myid == proc)
        {
            Grid *grid = GetGrid(iZone, this->level);
            CompressWallStructure(cdata, grid);
        }
    }

    WriteLogFile("  End compressing, start decompressing ...\n");

    //! Send the data-container to server.
    if (myid != serverTmp)
    {
        int tag = myid;
        send(cdata, serverTmp, tag);
        delete cdata;    cdata = nullptr;
    }
    else
    {
        int startProcessor = 0;
        int endProcessor = numberOfOriginalprocessor;

        if (CurrentProcessorIsGridProcessor())
        {
            startProcessor = numberOfOriginalprocessor;
            endProcessor = numberOfOriginalprocessor + numberOfGridprocessor;
        }

        for (int iProc = startProcessor; iProc < endProcessor; ++ iProc)
        {
            int tag = iProc;

            if (myid != iProc)
            {
                delete cdata;    cdata = nullptr;
                cdata = new DataContainer();
                receive(cdata, iProc, tag);
            }

            int nZonesInProc;
            cdata->MoveToBegin();
            PHRead(cdata, &nZonesInProc, 1);
            for (int iZone = 0; iZone < nZonesInProc; ++ iZone)
            {
                DeCompressWallStructure(cdata, false);
            }

            delete cdata;    cdata = nullptr;
        }
    }
}

void Pre_WalldistCompute::ServerBcast()
{
    //WriteLogFile("Wall structure server broad casting ...");
    using namespace PHMPI;
    int number_of_processor = GetNumberOfProcessor();
    if (number_of_processor <= 1) return;

    int tag = 0;
    int serverTmp = GetServerSepMode();

    DataContainer *cdata = new DataContainer();

    if (myid == serverTmp)
    {
        cdata->MoveToBegin();

        uint_t nWall = wallstructureList->size();
        PHWrite(cdata, &nWall, 1);
        for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
        {
            WallStructure *wallStruct = (*wallstructureList)[iWall];
            wallStruct->Encode(cdata, false);
        }
    }

    PH_BcastByServerSepMode(cdata, tag);

    if (myid != serverTmp)
    {
        cdata->MoveToBegin();

        uint_t nWall;
        PHRead(cdata, &nWall, 1);

        wallstructureList->resize(nWall);
        for (int iWall = 0; iWall < nWall; ++ iWall)
        {
            WallStructure *wallStruct = new WallStructure();

            wallStruct->Decode(cdata, false);

            (*wallstructureList)[iWall] = wallStruct;
        }
    }

    delete cdata;    cdata = nullptr;
}

void Pre_WalldistCompute::CompressWallStructure(DataContainer *&cdata, Grid *grid)
{
    if (grid->Type() == PHSPACE::STRUCTGRID)
    {
        StructGrid *structGrid = StructGridCast(grid);
        CompressWallStructure(cdata, structGrid);
    }
    else
    {
        UnstructGrid *unstructGrid = UnstructGridCast(grid);
        CompressWallStructure(cdata, unstructGrid);
    }
}

void Pre_WalldistCompute::CompressWallStructure(DataContainer *&cdata, StructGrid *grid)
{
    int nSolidFace = 0;

    StructBCSet *structBCSet = grid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
        int bctype = bcregion->GetBCType();

        if (IsWall(bctype))
        {
            int ist = 0, ied = 0, jst = 0, jed = 0, kst = 0, ked = 0;
            bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            nSolidFace += (ied - ist + 1) * (jed - jst + 1) * (ked - kst + 1);
        }
    }

    cdata->Write(&nSolidFace, sizeof(int));
    if (!nSolidFace)
    {
        return;
    }

    int iBlock = grid->GetIBlock();
    cdata->Write(&iBlock, sizeof(int));

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
        int bctype = bcregion->GetBCType();

        if (IsWall(bctype))
        {
            int ist = 0, ied = 0, jst = 0, jed = 0, kst = 0, ked = 0;
            bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            int nsurf = bcregion->GetFaceDirection() + 1;

            int in = 0, jn = 0, kn = 0;

            int nm = 3;
            int nlen = nm * (ked - kst + 1) * (jed - jst + 1) * (ied - ist + 1);
            RDouble *qtmp = new RDouble [nlen];
            int countTmp = 0;            
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        bcregion->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                        RDouble xc, yc, zc;
                        grid->FaceCoor(in, jn, kn, nsurf, xc, yc, zc);

                        qtmp[countTmp++] = xc;
                        qtmp[countTmp++] = yc;
                        qtmp[countTmp++] = zc;
                    }
                }
            }
            cdata->Write(qtmp, nlen * sizeof(RDouble));
            delete [] qtmp;     qtmp = nullptr;
        }
    }

    //! Compute nWallPoints.
    int nWallPoints = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
        int bctype = bcregion->GetBCType();

        if (IsWall(bctype))
        {
            int ist = 0, ied = 0, jst = 0, jed = 0, kst = 0, ked = 0;
            bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            //! Get the node scope.
            int s_st[3] = {ist, jst, kst};
            int s_ed[3] = {ied, jed, ked};

            for (int m = 0; m < GetDim(); ++ m)
            {
                if (bcregion->GetFaceDirection() == m)
                {
                    if (s_st[m] != 1)
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

            int ni = s_ed[0] - s_st[0] + 1;
            int nj = s_ed[1] - s_st[1] + 1;
            int nk = s_ed[2] - s_st[2] + 1;

            nWallPoints += ni * nj * nk;
        }
    }
    cdata->Write(&nWallPoints, sizeof(int));

    vector<Int3D *> structNode2UnstrNodeID;
    structNode2UnstrNodeID.resize(nBCRegion);
    for (int iBC = 0; iBC < nBCRegion; ++ iBC)
    {
        structNode2UnstrNodeID[iBC] = 0;
    }

    RDouble3D &xx = *grid->GetStructX();
    RDouble3D &yy = *grid->GetStructY();
    RDouble3D &zz = *grid->GetStructZ();

    nWallPoints = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
        int bctype = bcregion->GetBCType();

        if (IsWall(bctype))
        {
            int ist = 0, ied = 0, jst = 0, jed = 0, kst = 0, ked = 0;
            bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            //! Get the node scope.
            int s_st[3] = {ist, jst, kst};
            int s_ed[3] = {ied, jed, ked};

            for (int m = 0; m < GetDim(); ++ m)
            {
                if (bcregion->GetFaceDirection() == m)
                {
                    if (s_st[m] != 1)
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

            int ni = s_ed[0] - s_st[0] + 1;
            int nj = s_ed[1] - s_st[1] + 1;
            int nk = s_ed[2] - s_st[2] + 1;
            Range I(1, ni);
            Range J(1, nj);
            Range K(1, nk);

            structNode2UnstrNodeID[iBCRegion] = new Int3D(I, J, K, fortranArray);
            Int3D &nodeMap = *(structNode2UnstrNodeID[iBCRegion]);

            int nm = 3;
            int nlen = nm * (s_ed[2] - s_st[2] + 1) * (s_ed[1] - s_st[1] + 1) * (s_ed[0] - s_st[0] + 1);
            RDouble *qtmp = new RDouble [nlen];
            int countTmp = 0;  
            for (int k = s_st[2]; k <= s_ed[2]; ++ k)
            {
                for (int j = s_st[1]; j <= s_ed[1]; ++ j)
                {
                    for (int i = s_st[0]; i <= s_ed[0]; ++ i)
                    {
                        qtmp[countTmp++] = xx(i, j, k);
                        qtmp[countTmp++] = yy(i, j, k);
                        qtmp[countTmp++] = zz(i, j, k);

                        int iLocal = i - s_st[0] + 1;
                        int jLocal = j - s_st[1] + 1;
                        int kLocal = k - s_st[2] + 1;
                        nodeMap(iLocal, jLocal, kLocal) = nWallPoints;
                        ++ nWallPoints;
                    }
                }
            }
            cdata->Write(qtmp, nlen * sizeof(RDouble));
            delete [] qtmp;     qtmp = nullptr;
        }
    }

    //! Third: face to node connection.
    int count = 0;
    int *nPointsPerWallFace = new int [nSolidFace]();
    if (GetDim() == TWO_D)
    {
        for (int iFace = 0; iFace < nSolidFace; ++ iFace)
        {
            nPointsPerWallFace[iFace] = 2;
            count += 2;
        }
    }
    else
    {
        //! pole?
        for (int iFace = 0; iFace < nSolidFace; ++ iFace)
        {
            nPointsPerWallFace[iFace] = 4;
            count += 4;
        }
    }

    cdata->Write(nPointsPerWallFace, nSolidFace * sizeof(int));

    //! Wall face to node.
    int *wallface2node = new int [count]();

    int nk = grid->GetNK();
    count = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
        int bctype = bcregion->GetBCType();

        if (IsWall(bctype))
        {
            int ist = 0, ied = 0, jst = 0, jed = 0, kst = 0, ked = 0;
            bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            //! Get the node scope.
            int s_st[3] = {ist, jst, kst};
            int s_ed[3] = {ied, jed, ked};

            for (int m = 0; m < GetDim(); ++ m)
            {
                if (bcregion->GetFaceDirection() == m)
                {
                    if (s_st[m] != 1)
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

            int nsurf = bcregion->GetFaceDirection() + 1;
            Int3D &nodeMap = *(structNode2UnstrNodeID[iBCRegion]);
            int in, jn, kn;
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        bcregion->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                        int il1 = 1;
                        int jl1 = 1;
                        int kl1 = 1;
                        if (nk == 1) kl1 = 0;

                        int il, jl, kl;

                        il = in + il1;
                        jl = jn + jl1;
                        kl = kn + kl1;

                        //! Map the global node to local node ID.

                        int iLocal = in - s_st[0] + 1;
                        int jLocal = jn - s_st[1] + 1;
                        int kLocal = kn - s_st[2] + 1;

                        int ilLocal = il - s_st[0] + 1;
                        int jlLocal = jl - s_st[1] + 1;
                        int klLocal = kl - s_st[2] + 1;

                        if (GetDim() == THREE_D)
                        {
                            //! 3D.
                            if (nsurf == 1)
                            {
                                wallface2node[count ++] = nodeMap(iLocal, jLocal, kLocal);
                                wallface2node[count ++] = nodeMap(iLocal, jlLocal, kLocal);
                                wallface2node[count ++] = nodeMap(iLocal, jlLocal, klLocal);
                                wallface2node[count ++] = nodeMap(iLocal, jLocal, klLocal);
                            }
                            else if (nsurf == 2)
                            {
                                wallface2node[count ++] = nodeMap(iLocal, jLocal, kLocal);
                                wallface2node[count ++] = nodeMap(ilLocal, jLocal, kLocal);
                                wallface2node[count ++] = nodeMap(ilLocal, jLocal, klLocal);
                                wallface2node[count ++] = nodeMap(iLocal, jLocal, klLocal);
                            }
                            else if (nsurf == 3)
                            {
                                wallface2node[count ++] = nodeMap(iLocal, jLocal, kLocal);
                                wallface2node[count ++] = nodeMap(ilLocal, jLocal, kLocal);
                                wallface2node[count ++] = nodeMap(ilLocal, jlLocal, kLocal);
                                wallface2node[count ++] = nodeMap(iLocal, jlLocal, kLocal);
                            }
                        }
                        else
                        {
                            //! 2D.
                            if (nsurf == 1)
                            {
                                wallface2node[count ++] = nodeMap(iLocal, jLocal, kLocal);
                                wallface2node[count ++] = nodeMap(iLocal, jlLocal, kLocal);
                            }
                            else if (nsurf == 2)
                            {
                                wallface2node[count ++] = nodeMap(iLocal, jLocal, kLocal);
                                wallface2node[count ++] = nodeMap(ilLocal, jLocal, kLocal);
                            }
                        }
                    }
                }
            }
        }
    }

    count = 0;
    for (int iFace = 0; iFace < nSolidFace; ++ iFace)
    {
        int nNPF = nPointsPerWallFace[iFace];
        for (int j = 0; j < nNPF; ++ j)
        {
            int nodeIndex = wallface2node[count ++];
            cdata->Write(&nodeIndex, sizeof(int));
        }
    }

    //! Face normal.
    RDouble4D &xfn = *(grid->GetFaceNormalX());
    RDouble4D &yfn = *(grid->GetFaceNormalY());
    RDouble4D &zfn = *(grid->GetFaceNormalZ());
    
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
        int bctype = bcregion->GetBCType();

        if (IsWall(bctype))
        {
            int nsurf = bcregion->GetFaceDirection() + 1;

            int ist = 0, ied = 0, jst = 0, jed = 0, kst = 0, ked = 0;
            bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            int nm = 3;
            int nlen = nm * (ked - kst + 1) * (jed - jst + 1) * (ied - ist + 1);
            RDouble *qtmp = new RDouble [nlen];
            int countTmp = 0;  
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        int iFace, jFace, kFace;
                        bcregion->GetBoundaryFaceIndex(i, j, k, iFace, jFace, kFace);

                        qtmp[countTmp++] = xfn(iFace, jFace, kFace, nsurf);
                        qtmp[countTmp++] = yfn(iFace, jFace, kFace, nsurf);
                        qtmp[countTmp++] = zfn(iFace, jFace, kFace, nsurf);
                    }
                }
            }
            cdata->Write(qtmp, nlen * sizeof(RDouble));
            delete [] qtmp;     qtmp = nullptr;
        }
    }

    delete [] nPointsPerWallFace;    nPointsPerWallFace = NULL;
    delete [] wallface2node;    wallface2node = NULL;

    for (int iBC = 0; iBC < nBCRegion; ++ iBC)
    {
        delete structNode2UnstrNodeID[iBC];
    }
}

void Pre_WalldistCompute::CompressWallStructure(DataContainer *&cdata, UnstructGrid *grid)
{
    int nBoundFace = grid->GetNBoundFace();
    int nTotalNode = grid->GetNTotalNode();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();
    int count = 0;

    int *face2node = grid->GetFace2Node();
    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();

    int nSolidFace = 0;

    //! Compute the number of the wall faces.
    bool *isWallFace = new bool [nBoundFace]();
    SetField(isWallFace, false, nBoundFace);

    grid->FindOutWallFaces(isWallFace, nSolidFace);

    cdata->Write(&nSolidFace, sizeof(int));

    if (!nSolidFace)
    {
        //! If there are not solid cell, do nothing and return.
        delete [] isWallFace;    isWallFace = nullptr;
        return;
    }

    //! Write the block information into wallstruc.
    int iBlock = grid->GetIBlock();
    cdata->Write(&iBlock, sizeof(int));

    int *localWallFace2Global = new int [nSolidFace]();
    SetField(localWallFace2Global, -1, nSolidFace);
    count = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        if (isWallFace[iFace])
        {
            localWallFace2Global[count] = iFace;
            ++ count;
        }
    }

    //! Construct the wall faces structure.
    bool *isWallPoint = new bool [nTotalNode]();
    int nWallPoints = 0;
    grid->FindOutWallPoints(isWallPoint, nWallPoints);

    int *globalPoint2Wall = new int [nTotalNode]();
    PHSPACE::SetField(globalPoint2Wall, -1, nTotalNode);

    RDouble *xWall = new RDouble [nWallPoints]();
    RDouble *yWall = new RDouble [nWallPoints]();
    RDouble *zWall = new RDouble [nWallPoints]();
    count = 0;
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (isWallPoint[iNode])
        {
            xWall[count] = x[iNode];
            yWall[count] = y[iNode];
            zWall[count] = z[iNode];

            globalPoint2Wall[iNode] = count;
            ++ count;
        }
    }

    //! Number of nodes per wall face.
    count = 0;
    int *nPointsPerWallFace = new int [nSolidFace]();
    for (int iFace = 0; iFace < nSolidFace; ++ iFace)
    {
        int globalFaceIndex = localWallFace2Global[iFace];
        nPointsPerWallFace[iFace] = node_number_of_each_face[globalFaceIndex];
        count += nPointsPerWallFace[iFace];
    }

    //! Wall face to node.
    int *wallface2node = new int [count]();
    count = 0;
    int nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        if (isWallFace[iFace])
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                int index = face2node[nodepos + j];
                wallface2node[count ++] = globalPoint2Wall[index];
            }
        }

        nodepos += node_number_of_each_face[iFace];
    }

    //! Write the wall face data structure into the data-container.
    //! First: wall face center.
    int nm = 3;
    int nlen = nm * nSolidFace;
    RDouble *qtmp = new RDouble[nlen];
    int countTmp = 0;
    for (int iFace = 0; iFace < nSolidFace; ++ iFace)
    {
        int faceIndex = localWallFace2Global[iFace];

        qtmp[countTmp++] = xfc[faceIndex];
        qtmp[countTmp++] = yfc[faceIndex];
        qtmp[countTmp++] = zfc[faceIndex];
    }
    cdata->Write(qtmp, nlen * sizeof(RDouble));
    delete [] qtmp;     qtmp = nullptr;
    countTmp = 0;

    //! Second: points coordinates.
    cdata->Write(&nWallPoints, sizeof(int));

    nm = 3;
    nlen = nm * nWallPoints;
    qtmp = new RDouble[nlen];    
    for (int iNode = 0; iNode < nWallPoints; ++ iNode)
    {
        qtmp[countTmp++] = xWall[iNode];
        qtmp[countTmp++] = yWall[iNode];
        qtmp[countTmp++] = zWall[iNode];
    }
    cdata->Write(qtmp, nlen * sizeof(RDouble));
    delete [] qtmp;     qtmp = nullptr;
    countTmp = 0;

    //! Third: face to node connection.
    cdata->Write(nPointsPerWallFace, nSolidFace * sizeof(int));

    count = 0;
    for (int iFace = 0; iFace < nSolidFace; ++ iFace)
    {
        int nNPF = nPointsPerWallFace[iFace];
        for (int j = 0; j < nNPF; ++ j)
        {
            int nodeIndex = wallface2node[count ++];
            cdata->Write(&nodeIndex, sizeof(int));
        }
    }

    //! Forth: wall face normal.
    nm = 3;
    nlen = nm * nSolidFace;
    qtmp = new RDouble[nlen];  
    for (int iFace = 0; iFace < nSolidFace; ++ iFace)
    {
        int faceIndex = localWallFace2Global[iFace];

        qtmp[countTmp++] = xfn[faceIndex];
        qtmp[countTmp++] = yfn[faceIndex];
        qtmp[countTmp++] = zfn[faceIndex];
    }
    cdata->Write(qtmp, nlen * sizeof(RDouble));
    delete [] qtmp;     qtmp = nullptr;

    delete [] isWallPoint;    isWallPoint = nullptr;
    delete [] globalPoint2Wall;    globalPoint2Wall = nullptr;
    delete [] xWall;    xWall = nullptr;
    delete [] yWall;    yWall = nullptr;
    delete [] zWall;    zWall = nullptr;
    delete [] nPointsPerWallFace;    nPointsPerWallFace = nullptr;
    delete [] isWallFace;    isWallFace = nullptr;
    delete [] localWallFace2Global;    localWallFace2Global = nullptr;
    delete [] wallface2node;    wallface2node = nullptr;
}

void Pre_WalldistCompute::FilterWallStructureList()
{
    int nWallStruct = (int)this->wallstructureList->size();
    bool *isSelect = new bool [nWallStruct];
    SetField(isSelect, false, nWallStruct);

    //! Build the wall struct tree.
    WallStructTree *wallStructTree = new WallStructTree(6, wallMin, wallMax);
    for (int iWall = 0; iWall < nWallStruct; ++ iWall)
    {
        WallStructure *wallStruct = (*wallstructureList)[iWall];
        RDouble *box = wallStruct->GetBox();
        WallStructNode *wallNode = new WallStructNode(6, box, iWall);
        wallStructTree->AddNode(wallNode);
    }

    //! Near wall structure searching ...
    int nBoxLimitedMax = MAX(8 * (int)sqrt(nWallStruct * 1.0), 1);
    nBoxLimitedMax     = MIN(nBoxLimitedMax, nWallStruct);
    int nBoxLimitedMin = MAX(4 * (int)(sqrt(nWallStruct * 1.0)), 1);
    nBoxLimitedMin     = MAX(nBoxLimitedMin, 5);
    RDouble dis = 0.1 * wallTolerance;

    int nZones = PHMPI::GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int proc = PHMPI::GetZoneProcessorID(iZone);
        if (myid == proc)
        {
            Grid *grid = GetGrid(iZone, this->level);
            
            int boxMin = nBoxLimitedMin;
            int boxMax = nBoxLimitedMax;

            WallStructTree::AdtNodeList nList;
            bool finalFound = false;

            int nEnlargeTime = 3;
            int iEnlarge = 0;
            while (iEnlarge < nEnlargeTime)
            {
                ++ iEnlarge;

                //! Finding ...
                RDouble *gridMinBox = grid->GetMinBox();
                RDouble *gridMaxBox = grid->GetMaxBox();

                RDouble pmin[6], pmax[6];
                RDouble *rmin, *rmax;
                rmin = wallStructTree->GetMin();
                rmax = wallStructTree->GetMax();

                pmin[0] = rmin[0];
                pmin[1] = rmin[1];
                pmin[2] = rmin[2];
                pmax[3] = rmax[3];
                pmax[4] = rmax[4];
                pmax[5] = rmax[5];

                pmin[3] = gridMinBox[0] - dis;
                pmin[4] = gridMinBox[1] - dis;
                pmin[5] = gridMinBox[2] - dis;

                pmax[0] = gridMaxBox[0] + dis;
                pmax[1] = gridMaxBox[1] + dis;
                pmax[2] = gridMaxBox[2] + dis;

                RDouble size_min = 0;
                RDouble size_max = dis * 2.0;
                RDouble ra;

                int iSearch = 0;
                int nSearchTimes = 100;
                while (iSearch <nSearchTimes)
                {
                    ++ iSearch;

                    ra = 0.5 * (size_min + size_max);
                    pmin[3] = gridMinBox[0] - ra;
                    pmin[4] = gridMinBox[1] - ra;
                    pmin[5] = gridMinBox[2] - ra;

                    pmax[0] = gridMaxBox[0] + ra;
                    pmax[1] = gridMaxBox[1] + ra;
                    pmax[2] = gridMaxBox[2] + ra;

                    nList.clear();
                    wallStructTree->FindNodesInRegion(pmin, pmax, nList, boxMax + 1);

                    uint_t nNearWallStruct = nList.size();
                    if ((nNearWallStruct > boxMin && nNearWallStruct <= boxMax) ||
                        (iEnlarge == nEnlargeTime && nNearWallStruct > 0 && iSearch >= 0.8 * nSearchTimes))
                    {
                        finalFound = true;
                        break;
                    }
                    else if (nNearWallStruct > boxMax && (size_max - size_min) > 0.01 * dis)
                    {
                        size_max = ra;
                        continue;
                    }  
                    else
                    {
                        size_min = ra;
                        size_max = size_max * 2.0;
                        continue;
                    }
                }

                if (finalFound)
                {
                    for (std::size_t iWall = 0; iWall < nList.size(); ++ iWall)
                    {
                        int wallIndex = nList[iWall]->GetData();
                        isSelect[wallIndex] = true;
                    }
                    break;
                }
                else
                {
                    //! Enlarge the limitation.
                    boxMax *= 2;
                    continue;
                }
            }    //! Enlarge while.

            if (!finalFound)
            {
                //! Dump out the bad cell, that it is failure to find the right wall face.
                ostringstream oss;
                oss << "Error: could not find the nearest wall structure for zone " << 
                       grid->GetGridID()->GetIndex() << endl;
                oss << "Use other wall distance computing method please!" << endl;
                TK_Exit::ExceptionExit(oss.str());
            }
        }    // if (myid == proc)
    }        //! For each zone.

    //! Create new list.
    vector<WallStructure *> *newWallstructureList = new vector<WallStructure *>;
    for (int iWall = 0; iWall < nWallStruct; ++ iWall)
    {
        WallStructure *wallStruct = (*wallstructureList)[iWall];
        if (isSelect[iWall])
        {
            newWallstructureList->push_back(wallStruct);
        }
        else
        {
            delete wallStruct;
        }
    }

    if (newWallstructureList->size() == 0)
    {
        ostringstream oss1;
        oss1 << "Error: could not find the nearest wall structure in processor " << myid << endl;
        TK_Exit::ExceptionExit(oss1);
    }
    else
    {
        wallstructureList->clear();
        delete wallstructureList;
        wallstructureList = newWallstructureList;
    }

    ostringstream oss;
    oss << "Actual wall structure: " << newWallstructureList->size() << endl;
    WriteLogFile(oss);

    delete [] isSelect;    isSelect = NULL;
    delete wallStructTree;    wallStructTree = NULL;
    delete newWallstructureList;    newWallstructureList = NULL;
}

void Pre_WalldistCompute::FindNearstWallStructureOld(RDouble *cellCenter, vector<WallStructure *> &nearestWallStruct)
{
    int nWallStruct = (int)this->wallstructureList->size();

    if (walldistComputeMethod == Pre_WalldistCompute::ACCURATE_METHOD)
    {
        //! For method 1, loop over all of the wall structure.
        nearestWallStruct.resize(nWallStruct);
        for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
        {
            nearestWallStruct[iWall] = (*wallstructureList)[iWall];
        }
        return;
    }

    //! For method 0: loop over the nearest several box.
    RDouble x = cellCenter[0];
    RDouble y = cellCenter[1];
    RDouble z = cellCenter[2];

    //! Compute the distance between point input and each wall structure.
    RDouble xmin, ymin, zmin, xmax, ymax, zmax;
    for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
    {
        WallStructure *wallStruct = (*wallstructureList)[iWall];
        RDouble *box = wallStruct->GetBox();

        xmin = box[0]; xmax = box[3];
        ymin = box[1]; ymax = box[4];
        zmin = box[2]; zmax = box[5];

        RDouble dis;
        if (xmin > x)
        {
            if (ymin > y)
            {
                if (zmin > z)
                {
                    dis = sqrt((xmin - x) * (xmin - x) +
                               (ymin - y) * (ymin - y) +
                               (zmin - z) * (zmin - z));
                }
                else if (zmin <= z && z <= zmax)
                {
                    dis = sqrt((xmin - x) * (xmin - x) +
                               (ymin - y) * (ymin - y));
                }
                else
                {
                    dis = sqrt((xmin - x) * (xmin - x) +
                               (ymin - y) * (ymin - y) +
                               (zmax - z) * (zmax - z));
                }
            }
            else if (ymin <= y && y <= ymax)
            {
                if (zmin > z)
                {
                    dis = sqrt((xmin - x) * (xmin - x) +
                               (zmin - z) * (zmin - z));
                }
                else if (zmin <= z && z <= zmax)
                {
                    dis = sqrt((xmin - x) * (xmin - x));
                }
                else
                {
                    dis = sqrt((xmin - x) * (xmin - x) +
                               (zmax - z) * (zmax - z));
                }
            }
            else
            {
                if (zmin > z)
                {
                    dis = sqrt((xmin - x) * (xmin - x) +
                               (ymax - y) * (ymax - y) +
                               (zmin - z) * (zmin - z));
                }
                else if (zmin <= z && z <= zmax)
                {
                    dis = sqrt((xmin - x) * (xmin - x) +
                               (ymax - y) * (ymax - y));
                }
                else
                {
                    dis = sqrt((xmin - x) * (xmin - x) +
                               (ymax - y) * (ymax - y) +
                               (zmax - z) * (zmax - z));
                }
            }
        }
        else if (xmin <= x && x <= xmax)
        {
            if (ymin > y)
            {
                if (zmin > z)
                {
                    dis = sqrt((ymin - y) * (ymin - y) +
                               (zmin - z) * (zmin - z));
                }
                else if (zmin <= z && z <= zmax)
                {
                    dis = sqrt((ymin - y) * (ymin - y));
                }
                else
                {
                    dis = sqrt((ymin - y) * (ymin - y) +
                               (zmax - z) * (zmax - z));
                }
            }
            else if (ymin <= y && y <= ymax)
            {
                if (zmin > z)
                {
                    dis = sqrt((zmin - z) * (zmin - z));
                }
                else if (zmin <= z && z <= zmax)
                {
                    dis = -1;
                }
                else
                {
                    dis = sqrt((zmax - z) * (zmax - z));
                }
            }
            else
            {
                if (zmin > z)
                {
                    dis = sqrt((ymax - y) * (ymax - y) +
                               (zmin - z) * (zmin - z));
                }
                else if (zmin <= z && z <= zmax)
                {
                    dis = sqrt((ymax - y) * (ymax - y));
                }
                else
                {
                    dis = sqrt((ymax - y) * (ymax - y) +
                               (zmax - z) * (zmax - z));
                }
            }
        }
        else
        {
            if (ymin > y)
            {
                if (zmin > z)
                {
                    dis = sqrt((xmax - x) * (xmax - x) +
                               (ymin - y) * (ymin - y) +
                               (zmin - z) * (zmin - z));
                }
                else if (zmin <= z && z <= zmax)
                {
                    dis = sqrt((xmax - x) * (xmax - x) +
                               (ymin - y) * (ymin - y));
                }
                else
                {
                    dis = sqrt((xmax - x) * (xmax - x) +
                               (ymin - y) * (ymin - y) +
                               (zmax - z) * (zmax - z));
                }
            }
            else if (ymin <= y && y <= ymax)
            {
                if (zmin > z)
                {
                    dis = sqrt((xmax - x) * (xmax - x) +
                               (zmin - z) * (zmin - z));
                }
                else if (zmin <= z && z <= zmax)
                {
                    dis = sqrt((xmax - x) * (xmax - x));
                }
                else
                {
                    dis = sqrt((xmax - x) * (xmax - x) +
                               (zmax - z) * (zmax - z));
                }
            }
            else
            {
                if (zmin > z)
                {
                    dis = sqrt((xmax - x) * (xmax - x) +
                               (ymax - y) * (ymax - y) +
                               (zmin - z) * (zmin - z));
                }
                else if (zmin <= z && z <= zmax)
                {
                    dis = sqrt((xmax - x) * (xmax - x) +
                               (ymax - y) * (ymax - y));
                }
                else
                {
                    dis = sqrt((xmax - x) * (xmax - x) +
                               (ymax - y) * (ymax - y) +
                               (zmax - z) * (zmax - z));
                }
            }
        }

        wallStruct->SetDistance(dis);
    }

    //! Find out the nearest wall structure.
    int nBox = MAX(4 * (int)sqrt(nWallStruct * 1.0), 1);
    nBox     = MIN(nBox, nWallStruct);

    nearestWallStruct.reserve(nBox);
    bool *isSelected = new bool [nWallStruct];
    SetField(isSelected, false, nWallStruct);

    for (int iBox = 0; iBox < nBox; ++ iBox)
    {
        RDouble minDist = LARGE;

        int targetBoxID = -1;
        for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
        {
            WallStructure *wallStruct = (*wallstructureList)[iWall];
            RDouble distance = wallStruct->GetDistance();
            if (!isSelected[iWall] && distance < minDist)
            {
                minDist = distance;
                targetBoxID = static_cast<int>(iWall);
            }
        }
        if (targetBoxID == -1)
        {
            TK_Exit::ExceptionExit(" Error: non nearest wall structure was found!\n");
        }

        nearestWallStruct.push_back((*wallstructureList)[targetBoxID]);
        isSelected[targetBoxID] = true;
    }

    //! Add the neighbor overlapping boxes into the list.
    for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
    {
        if (!isSelected[iWall])
        {
            WallStructure *wallStruct0 = (*wallstructureList)[iWall];
            RDouble *box0 = wallStruct0->GetBox();
            RDouble xmin0 = box0[0]; RDouble xmax0 = box0[3];
            RDouble ymin0 = box0[1]; RDouble ymax0 = box0[4];
            RDouble zmin0 = box0[2]; RDouble zmax0 = box0[5];

            for (int iBox = 0; iBox < nBox; ++ iBox)
            {
                WallStructure *wallStruct1 = nearestWallStruct[iBox];
                RDouble *box1 = wallStruct1->GetBox();
                RDouble xmin1 = box1[0]; RDouble xmax1 = box1[3];
                RDouble ymin1 = box1[1]; RDouble ymax1 = box1[4];
                RDouble zmin1 = box1[2]; RDouble zmax1 = box1[5];

                if (xmin1 > xmax0 || xmax1 < xmin0 ||
                    ymin1 > ymax0 || ymax1 < ymin0 ||
                    zmin1 > zmax0 || zmax1 < zmin0)
                {
                    //! Non-overlapping.
                }
                else
                {
                    //! Overlapping.
                    nearestWallStruct.push_back(wallStruct0);
                    isSelected[iWall] = true;
                    break;
                }
            }
        }
    }

    delete [] isSelected;
}

void Pre_WalldistCompute::DeCompressWallStructure(DataContainer *&cdata, const bool &ifMoveToBegin)
{
    if (ifMoveToBegin)
    {
        cdata->MoveToBegin();
    }

    int nSolidFace;
    cdata->Read(&nSolidFace, sizeof(int));

    if (nSolidFace == 0)
    {
        //! If there are not solid cell, do nothing and return.
        return;
    }

    nTWFace += nSolidFace;
    
    WallStructure *wallStructure = new WallStructure();
    wallstructureList->push_back(wallStructure);

    //! Set iBlock to wallStructure.
    int iBlock;
    cdata->Read(&iBlock, sizeof(int));
    wallStructure->SetIBlock(iBlock);

    //! Face center.
    WallStructure::value_type &xfc = wallStructure->GetXFaceCenter();
    WallStructure::value_type &yfc = wallStructure->GetYFaceCenter();
    WallStructure::value_type &zfc = wallStructure->GetZFaceCenter();
    xfc.reserve(nSolidFace);
    yfc.reserve(nSolidFace);
    zfc.reserve(nSolidFace);
    for (int i = 0; i < nSolidFace; ++ i)
    {
        RDouble xc, yc, zc;
        cdata->Read(&xc, sizeof(RDouble));
        cdata->Read(&yc, sizeof(RDouble));
        cdata->Read(&zc, sizeof(RDouble));

        xfc.push_back(xc);
        yfc.push_back(yc);
        zfc.push_back(zc);
    }

    //! Wall points coordinates.
    int nWallPoints;
    cdata->Read(&nWallPoints, sizeof(int));
    if (nWallPoints == 0)
    {
        //! If the wall points are not including, return.
        return;
    }

    WallStructure::value_type &x = wallStructure->GetX();
    WallStructure::value_type &y = wallStructure->GetY();
    WallStructure::value_type &z = wallStructure->GetZ();
    x.reserve(nWallPoints);
    y.reserve(nWallPoints);
    z.reserve(nWallPoints);
    for (int i = 0; i < nWallPoints; ++ i)
    {
        RDouble x0, y0, z0;
        cdata->Read(&x0, sizeof(RDouble));
        cdata->Read(&y0, sizeof(RDouble));
        cdata->Read(&z0, sizeof(RDouble));

        x.push_back(x0);
        y.push_back(y0);
        z.push_back(z0);
    }

    int *nPointPerFace = new int [nSolidFace];

    int count = 0;
    for (int iFace = 0; iFace < nSolidFace; ++ iFace)
    {
        int nPoint;
        cdata->Read(&nPoint, sizeof(int));
        nPointPerFace[iFace] = nPoint;
        count += nPoint;
    }

    int *wallFace2Node = new int [count];
    count = 0;
    for (int iFace = 0; iFace < nSolidFace; ++ iFace)
    {
        for (int j = 0; j < nPointPerFace[iFace]; ++ j)
        {
            int nodeIndex;
            cdata->Read(&nodeIndex, sizeof(int));
            wallFace2Node[count ++] = nodeIndex;
        }
    }

    //! Analyze node position of each face.
    int *nodePosition = new int [nSolidFace + 1];
    nodePosition[0] = 0;
    int positionCount = 0;
    for (int iFace = 0; iFace < nSolidFace; ++ iFace)
    {
        positionCount += nPointPerFace[iFace];
        nodePosition[iFace + 1] = positionCount;
    }
    wallStructure->SetNodePosition(nodePosition);

    //! Face normal.
    WallStructure::value_type &xfn = wallStructure->GetXFaceNormal();
    WallStructure::value_type &yfn = wallStructure->GetYFaceNormal();
    WallStructure::value_type &zfn = wallStructure->GetZFaceNormal();
    xfn.reserve(nSolidFace);
    yfn.reserve(nSolidFace);
    zfn.reserve(nSolidFace);
    for (int i = 0; i < nSolidFace; ++ i)
    {
        RDouble xfn0, yfn0, zfn0;
        cdata->Read(&xfn0, sizeof(RDouble));
        cdata->Read(&yfn0, sizeof(RDouble));
        cdata->Read(&zfn0, sizeof(RDouble));

        xfn.push_back(xfn0);
        yfn.push_back(yfn0);
        zfn.push_back(zfn0);
    }

    wallStructure->SetnPointPerFace(nPointPerFace);
    wallStructure->SetWallFace2Node(wallFace2Node);
    wallStructure->Box();
}

void Pre_WalldistCompute::BuildWallStructTree()
{
    //PrintToWindow("Wall structure tree building ...\n");
    //WriteLogFile("Wall structure tree building ...\n");

    //! Compute the boundary.
    RDouble wallBound[6];
    for (int iDim = 0; iDim < 3; ++ iDim)
    {
        wallBound[iDim]     = LARGE;
        wallBound[iDim + 3] = -LARGE;
    }

    RDouble oldPoint[3], newPoint[3];

    for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
    {
        WallStructure *wallStruct = (*wallstructureList)[iWall];
        RDouble *box = wallStruct->GetBox();

        RDouble x0 = box[0], y0 = box[1], z0 = box[2];
        RDouble x1 = box[3], y1 = box[4], z1 = box[5];

        //! Check each of the eight points of the CUBE box!
        RDouble xRange[2], yRange[2], zRange[2];
        xRange[0] = x0; xRange[1] = x1;
        yRange[0] = y0; yRange[1] = y1;
        zRange[0] = z0; zRange[1] = z1;
        for (int i = 0; i < 2; ++ i)
        {
            for (int j = 0; j < 2; ++ j)
            {
                for (int k = 0; k < 2; ++ k)
                {
                    oldPoint[0] = xRange[i];
                    oldPoint[1] = yRange[j];
                    oldPoint[2] = zRange[k];

                    this->RotateAxis(oldPoint, newPoint);

                    wallBound[0] = MIN(newPoint[0], wallBound[0]);
                    wallBound[1] = MIN(newPoint[1], wallBound[1]);
                    wallBound[2] = MIN(newPoint[2], wallBound[2]);

                    wallBound[3] = MAX(newPoint[0], wallBound[3]);
                    wallBound[4] = MAX(newPoint[1], wallBound[4]);
                    wallBound[5] = MAX(newPoint[2], wallBound[5]);
                }
            }
        }
    }

    wallTolerance = 0;
    for (int iDim = 0; iDim < 3; ++ iDim)
    {
        wallTolerance = MAX(wallTolerance, wallBound[iDim + 3] - wallBound[iDim]);
    }
    wallTolerance *= 0.01;

    for (int iDim = 0; iDim < 3; ++ iDim)
    {
        wallBound[iDim]     -= wallTolerance;
        wallBound[iDim + 3] += wallTolerance;
    }

    //! Build Tree.
    for (int iDim = 0; iDim < 3; ++ iDim)
    {
        wallMin[iDim] = wallBound[iDim];
        wallMax[iDim] = wallBound[iDim + 3];
    }
    for (int iDim = 3; iDim < 6; ++ iDim)
    {
        wallMin[iDim] = wallBound[iDim - 3];
        wallMax[iDim] = wallBound[iDim];
    }

    if (!wallFaceTree)
    {
        wallFaceTree = new WallFaceTree(6, wallMin, wallMax);
    }

    if (walldistComputeMethod == SUPER_FAST_METHOD)
    {
        if (rotateAxis > 0)
        {
            TK_Exit::ExceptionExit("Error: has not considered rotate coordinates for super fast wall dist computing method!");
        }
        FilterWallStructureList();
    }

    RDouble triBox[6];
    RDouble oldCoord[3], newCoord[3];

    nTWallFace = 0;
    for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
    {
        WallStructure *cwalldist = (*wallstructureList)[iWall];

        if (cwalldist->GetIBlock() != presentBlock && walldistComputeMethod == KDTREE_METHOD)
        {
            continue;
        }

        uint_t nWallFaces = cwalldist->GetNumberOfWallFaces();
        nTWallFace += nWallFaces;

        WallStructure::value_type &x = cwalldist->GetX();
        WallStructure::value_type &y = cwalldist->GetY();
        WallStructure::value_type &z = cwalldist->GetZ();

        int *nPointPerFace = cwalldist->GetnPointPerFace();
        int *wallFace2Node = cwalldist->GetWallFace2Node();

        int count = 0;
        for (int iFace = 0; iFace < nWallFaces; ++ iFace)
        {
            //! Insert into the tree.
            for (int iDim = 0; iDim < 3; ++ iDim)
            {
                triBox[iDim]     = LARGE;
                triBox[iDim + 3] = -LARGE;
            }

            int nNode = nPointPerFace[iFace];
            for (int iNode = 0; iNode < nNode; ++ iNode)
            {
                int nodeID = wallFace2Node[count ++];

                oldCoord[0] = x[nodeID];
                oldCoord[1] = y[nodeID];
                oldCoord[2] = z[nodeID];

                this->RotateAxis(oldCoord, newCoord);
                triBox[0] = MIN(triBox[0], newCoord[0]);
                triBox[1] = MIN(triBox[1], newCoord[1]);
                triBox[2] = MIN(triBox[2], newCoord[2]);

                triBox[3] = MAX(triBox[3], newCoord[0]);
                triBox[4] = MAX(triBox[4], newCoord[1]);
                triBox[5] = MAX(triBox[5], newCoord[2]);
            }

            WallFaceNode *triNode = new WallFaceNode(6, triBox, FaceID(iWall, iFace));
            wallFaceTree->AddNode(triNode);
        }
    }

    //PrintToWindow("Wall structure tree building over ...\n");
    //WriteLogFile("Wall structure tree building over ...\n");
}

bool Pre_WalldistCompute::IsWalldistComputOver()
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    bool isOver = true;
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Grid *grid = GetGrid(zoneID, this->level);
        if (IsWalldistComputOver(grid) == false)
        {
            isOver = false;
            break;
        }
    }
    return isOver;
}

bool Pre_WalldistCompute::IsWalldistComputOver(Grid *grid)
{
    if (grid->Type() == PHSPACE::STRUCTGRID)
    {
        StructGrid *structGrid = StructGridCast(grid);
        return IsWalldistComputOver(structGrid);
    }
    else
    {
        UnstructGrid *unstructGrid = UnstructGridCast(grid);
        return IsWalldistComputOver(unstructGrid);
    }
}

bool Pre_WalldistCompute::IsWalldistComputOver(StructGrid *grid)
{
    bool isOver = true;

    RDouble3D &walldist = *grid->GetWallDist();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                if (walldist(i, j, k) >= LARGE * 0.9)
                {
                    isOver = false;
                    return isOver;
                }
            }
        }
    }
    return isOver;
}

bool Pre_WalldistCompute::IsWalldistComputOver(UnstructGrid *grid)
{
    bool isOver = true;

    if (cellMethodOrNodeMethod == CELL_METHOD)
    {
        RDouble *walldist = grid->GetWallDist();

        int nTotalCell = grid->GetNTotalCell();

        int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
        int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
        int *wallCellLabel = grid->GetWallCellLabel();

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            if (viscousType <= LAMINAR && wallCellLabel[iCell] != 1 && iLES == NOLES_SOLVER) continue;

            if (walldist[iCell] >= LARGE * 0.9)
            {
                isOver = false;
                break;
            }
        }
    }
    else
    {
        RDouble *wallDistNode = grid->GetWallDistNode();

        int nTotalNode = grid->GetNTotalNode();
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            if (wallDistNode[iNode] >= LARGE * 0.9)
            {
                isOver = false;
                break;
            }
        }
    }

    return isOver;
}

bool Pre_WalldistCompute::IsNodeWalldistComputOver()
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        UnstructGrid *grid = UnstructGridCast(GetGrid(zoneID, this->level));

        if (grid->GetIBlock() != presentBlock)
        {
            continue;
        }

        vector<RDouble> &nodeWalldist = grid->GetNodeWallDist();

        int nTotalNode = grid->GetNTotalNode();
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            if (nodeWalldist[iNode] >= LARGE * 0.9)
            {
                return false;
            }
        }
    }

    return true;
}

void Pre_WalldistCompute::CheckWallDistance()
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Grid *grid = GetGrid(zoneID, this->level);

        if ((grid->Type() == PHSPACE::UNSTRUCTGRID) && (grid->GetDim() == PHSPACE::THREE_D))
        {
            UnstructGrid *unstructGrid = UnstructGridCast(grid);
            if (cellMethodOrNodeMethod == CELL_METHOD)
            {
                CheckWallDistance(unstructGrid);
            }
            else
            {
                CheckWallDistanceByNodeMethod(unstructGrid);
            }
        }
    }
}

void Pre_WalldistCompute::CheckWallDistanceByNodeMethod(UnstructGrid *grid)
{
    using namespace PHMPI;

    RDouble *wallDistNode = grid->GetWallDistNode();

    int nTotalNode = grid->GetNTotalNode();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble nodeCoord[3], faceCenter[3], faceNormal[3];

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (wallDistNode[iNode] < LARGE * 0.9)
        {
            continue;
        }

        nodeCoord[0] = x[iNode];
        nodeCoord[1] = y[iNode];
        nodeCoord[2] = z[iNode];

        for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
        {
            WallStructure *cwalldist = (*wallstructureList)[iWall];

            uint_t nWallFaces              = cwalldist->GetNumberOfWallFaces();
            WallStructure::value_type &xfc = cwalldist->GetXFaceCenter();
            WallStructure::value_type &yfc = cwalldist->GetYFaceCenter();
            WallStructure::value_type &zfc = cwalldist->GetZFaceCenter();
            WallStructure::value_type &xW   = cwalldist->GetX();
            WallStructure::value_type &yW   = cwalldist->GetY();
            WallStructure::value_type &zW   = cwalldist->GetZ();
            WallStructure::value_type &xfn = cwalldist->GetXFaceNormal();
            WallStructure::value_type &yfn = cwalldist->GetYFaceNormal();
            WallStructure::value_type &zfn = cwalldist->GetZFaceNormal();

            int *nPointPerFace = cwalldist->GetnPointPerFace();
            int *wallFace2Node = cwalldist->GetWallFace2Node();

            int nodeCount = 0;
            for (int iFace = 0; iFace < nWallFaces; ++ iFace)
            {
                faceCenter[0] = xfc[iFace];
                faceCenter[1] = yfc[iFace];
                faceCenter[2] = zfc[iFace];

                faceNormal[0] = xfn[iFace];
                faceNormal[1] = yfn[iFace];
                faceNormal[2] = zfn[iFace];

                vector<int> face2node;
                face2node.resize(nPointPerFace[iFace]);
                for (int jNode = 0; jNode < nPointPerFace[iFace]; ++jNode)
                {
                    face2node[jNode] = wallFace2Node[nodeCount++];
                }

                RDouble wdst = ComputeDistanceOfANodeToAFace(nodeCoord, faceCenter, faceNormal, face2node, xW, yW, zW);
                if (wallDistNode[iNode] > wdst)
                {
                    wallDistNode[iNode] = wdst;
                }
            }
        }
    }
}

void Pre_WalldistCompute::CheckWallDistance(UnstructGrid *grid)
{
    using namespace PHMPI;

    RDouble *walldist = grid->GetWallDist();

    int nTotalCell = grid->GetNTotalCell();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
    int *wallCellLabel = grid->GetWallCellLabel();

    cout.unsetf(ios::scientific);
    cout.setf(ios::fixed);
    cout.precision(0);

    RDouble cellCenter[3], faceCenter[3], faceNormal[3];

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (viscousType <= LAMINAR && wallCellLabel[iCell] != 1 && iLES == NOLES_SOLVER) continue;

        if (walldist[iCell] < LARGE * 0.9)
        {
            continue;
        }

        cellCenter[0] = xcc[iCell];
        cellCenter[1] = ycc[iCell];
        cellCenter[2] = zcc[iCell];

        for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
        {
            WallStructure *cwalldist = (*wallstructureList)[iWall];

            uint_t nWallFaces = cwalldist->GetNumberOfWallFaces();
            WallStructure::value_type &xfc = cwalldist->GetXFaceCenter();
            WallStructure::value_type &yfc = cwalldist->GetYFaceCenter();
            WallStructure::value_type &zfc = cwalldist->GetZFaceCenter();
            WallStructure::value_type &x = cwalldist->GetX();
            WallStructure::value_type &y = cwalldist->GetY();
            WallStructure::value_type &z = cwalldist->GetZ();
            WallStructure::value_type &xfn = cwalldist->GetXFaceNormal();
            WallStructure::value_type &yfn = cwalldist->GetYFaceNormal();
            WallStructure::value_type &zfn = cwalldist->GetZFaceNormal();

            int *nPointPerFace = cwalldist->GetnPointPerFace();
            int *wallFace2Node = cwalldist->GetWallFace2Node();

            int nodeCount = 0;
            for (int iFace = 0; iFace < nWallFaces; ++ iFace)
            {
                faceCenter[0] = xfc[iFace];
                faceCenter[1] = yfc[iFace];
                faceCenter[2] = zfc[iFace];

                faceNormal[0] = xfn[iFace];
                faceNormal[1] = yfn[iFace];
                faceNormal[2] = zfn[iFace];

                vector<int> face2node;
                face2node.resize(nPointPerFace[iFace]);
                for (int iNode = 0; iNode < nPointPerFace[iFace]; ++ iNode)
                {
                    face2node[iNode] = wallFace2Node[nodeCount ++];
                }

                RDouble wdst = ComputeDistanceOfANodeToAFace(cellCenter, faceCenter, faceNormal, face2node, x, y, z);
                if (walldist[iCell] > wdst)
                {
                    walldist[iCell] = wdst;
                }
            }
        }
    }
}

void Pre_WalldistCompute::PostWallDistance()
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Grid *grid = GetGrid(zoneID, this->level);
        PostWallDistance(grid);
    }
}

void Pre_WalldistCompute::PostWallDistance(Grid *grid)
{
    if (grid->Type() == PHSPACE::STRUCTGRID)
    {
        StructGrid *structGrid = StructGridCast(grid);
        PostWallDistance(structGrid);
    }
    else
    {
        UnstructGrid *unstructGrid = UnstructGridCast(grid);

        if (cellMethodOrNodeMethod == CELL_METHOD)
        {
            PostWallDistance(unstructGrid);
        }
        else
        {
            PostWallDistanceByNodeMethod(unstructGrid);
        }
    }
}

void Pre_WalldistCompute::PostWallDistance(StructGrid *grid)
{
    RDouble3D &walldist = *grid->GetWallDist();
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                walldist(i, j, k) = sqrt(walldist(i, j, k));
            }
        }
    }
    GhostCell3D(walldist, ni, nj, nk);
}

void Pre_WalldistCompute::PostWallDistanceByNodeMethod(UnstructGrid *grid)
{
    RDouble *walldist = grid->GetWallDist();
    RDouble *wallDistNode = grid->GetWallDistNode();

    int nTotalNode = grid->GetNTotalNode();
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        wallDistNode[iNode] = sqrt(wallDistNode[iNode]);
    }

    int nTotalCell = grid->GetNTotalCell();
    int **cellNodeIndex = grid->GetCell2NodeArray();
    int  *cellNodeNumber = grid->GetNodeNumberOfEachCell();
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int nNode = cellNodeNumber[iCell];
        walldist[iCell] = 0.0;
        for (int iNode = 0; iNode < nNode; ++ iNode)
        {
            int nodeID = cellNodeIndex[iCell][iNode];
            walldist[iCell] += wallDistNode[nodeID];

        }
        walldist[iCell] /= nNode;
    }

    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();
    int *face2node = grid->GetFace2Node();

    int nBoundFace = grid->GetNBoundFace();
    UnstructBCSet **bcr = grid->GetBCRecord();

    bool *nodeOnWall = new bool[nTotalNode];
    std::fill_n(nodeOnWall, nTotalNode, false);

    int nodeCount = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcr[iFace]->GetKey();
        if (!(IsWall(bcType)))
        {
            nodeCount = nodeCount + node_number_of_each_face[iFace];
            continue;
        }

        for (int iNode = 0; iNode < node_number_of_each_face[iFace]; ++ iNode)
        {
            int nodeIndex = face2node[nodeCount];
            nodeOnWall[nodeIndex] = true;
            nodeCount++;
        }
    }

    int countCellOnWall = 0;
    vector < int > cellIndexOnWall;
    cellIndexOnWall.resize(0);
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        bool allNodeOnWall = true;
        int nNodes = cellNodeNumber[iCell];
        for (int iNode = 0; iNode < nNodes; ++ iNode)
        {
            int nodeIndex = cellNodeIndex[iCell][iNode];
            if (!(nodeOnWall[nodeIndex]))
            {
                allNodeOnWall = false;
                break;
            }
        }

        if (allNodeOnWall)
        {
            cellIndexOnWall.push_back(iCell);
            countCellOnWall++;
        }
    }

    if (countCellOnWall > 0)
    {
        RDouble *xcc = grid->GetCellCenterX();
        RDouble *ycc = grid->GetCellCenterY();
        RDouble *zcc = grid->GetCellCenterZ();

        RDouble cellCenter[3], faceCenter[3], faceNormal[3];
        for (int iCell = 0; iCell < countCellOnWall; ++ iCell)
        {
            int cellIndex = cellIndexOnWall[iCell];
            walldist[cellIndex] = LARGE;

            cellCenter[0] = xcc[cellIndex];
            cellCenter[1] = ycc[cellIndex];
            cellCenter[2] = zcc[cellIndex];

            for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
            {
                WallStructure *cwalldist = (*wallstructureList)[iWall];

                uint_t nWallFaces              = cwalldist->GetNumberOfWallFaces();
                WallStructure::value_type &xfc = cwalldist->GetXFaceCenter();
                WallStructure::value_type &yfc = cwalldist->GetYFaceCenter();
                WallStructure::value_type &zfc = cwalldist->GetZFaceCenter();
                WallStructure::value_type &x   = cwalldist->GetX();
                WallStructure::value_type &y   = cwalldist->GetY();
                WallStructure::value_type &z   = cwalldist->GetZ();
                WallStructure::value_type &xfn = cwalldist->GetXFaceNormal();
                WallStructure::value_type &yfn = cwalldist->GetYFaceNormal();
                WallStructure::value_type &zfn = cwalldist->GetZFaceNormal();

                int *nPointPerFace = cwalldist->GetnPointPerFace();
                int *wallFace2Node = cwalldist->GetWallFace2Node();

                nodeCount = 0;
                for (int iFace = 0; iFace < nWallFaces; ++ iFace)
                {
                    faceCenter[0] = xfc[iFace];
                    faceCenter[1] = yfc[iFace];
                    faceCenter[2] = zfc[iFace];

                    faceNormal[0] = xfn[iFace];
                    faceNormal[1] = yfn[iFace];
                    faceNormal[2] = zfn[iFace];

                    vector<int> currentWallFace2Node;
                    currentWallFace2Node.resize(nPointPerFace[iFace]);
                    for (int iNode = 0; iNode < nPointPerFace[iFace]; ++ iNode)
                    {
                        currentWallFace2Node[iNode] = wallFace2Node[nodeCount++];
                    }

                    RDouble wdst = ComputeDistanceOfANodeToAFace(cellCenter, faceCenter, faceNormal, currentWallFace2Node, x, y, z);
                    RDouble iCellWalldist = sqrt(wdst);
                    if (walldist[cellIndex] > iCellWalldist)
                    {
                        walldist[cellIndex] = iCellWalldist;
                    }
                }
            }
        }
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (walldist[iCell] < 1.0e-10)
        {
            PrintToWindow("     Warning: all cell nodes may lay on solid wall !!!\n");
            ostringstream oss;
            oss << "  cellIndex = " << iCell << endl;
            PrintToWindow(oss);
        }
    }

    delete [] nodeOnWall;    nodeOnWall = nullptr;
}

void Pre_WalldistCompute::PostWallDistance(UnstructGrid *grid)
{
    RDouble *walldist = grid->GetWallDist();

    int nTotalCell = grid->GetNTotalCell();

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        walldist[iCell] = sqrt(walldist[iCell]);
    }
}
#pragma warning(disable:4100)
RDouble Pre_WalldistCompute::GetGlobalGridSize(int type)
{
    RDouble globalMinBox[3], globalMaxBox[3];
    GetGlobalMinMaxBox(globalMinBox, globalMaxBox);

    RDouble xMin = globalMinBox[0], yMin = globalMinBox[1], zMin = globalMinBox[2];
    RDouble xMax = globalMaxBox[0], yMax = globalMaxBox[1], zMax = globalMaxBox[2];
    RDouble dx = xMax - xMin, dy = yMax - yMin, dz = zMax - zMin;
    return DISTANCE(dx, dy, dz);
}
#pragma warning(default:4100)
int Pre_WalldistCompute::JudgeInitialRotateAxis()
{
    if (this->walldistComputeMethod == ACCURATE_METHOD)
    {
        return -1;
    }

    RDouble zeroTolerance = 1.0e-10;
    int nFaceX = 0, nFaceY = 0, nFaceZ = 0;

    for (std::size_t iWall = 0; iWall < wallstructureList->size(); ++ iWall)
    {
        WallStructure *cwalldist = (*wallstructureList)[iWall];

        if (cwalldist->GetIBlock() != presentBlock)
        {
            continue;
        }

        uint_t nWallFaces = cwalldist->GetNumberOfWallFaces();

        WallStructure::value_type &xfc = cwalldist->GetXFaceNormal();
        WallStructure::value_type &yfc = cwalldist->GetYFaceNormal();
        WallStructure::value_type &zfc = cwalldist->GetZFaceNormal();

        for (int iFace = 0; iFace < nWallFaces; ++ iFace)
        {
            //! If this face is perpendicular to the rotate axis, count it.
            RDouble fabsNormalX = fabs(xfc[iFace]);
            RDouble fabsNormalY = fabs(yfc[iFace]);
            RDouble fabsNormalZ = fabs(zfc[iFace]);

            if (fabsNormalY < zeroTolerance && fabsNormalZ < zeroTolerance)
            {
                ++ nFaceX;
            }
            else if (fabsNormalX < zeroTolerance && fabsNormalZ < zeroTolerance)
            {
                ++ nFaceY;
            }
            else if (fabsNormalX < zeroTolerance && fabsNormalY < zeroTolerance)
            {
                ++ nFaceZ;
            }
        }
    }

    int rotateDirection = 1;
    if (nFaceX < 500 && nFaceY < 500 && nFaceZ < 500)
    {
        //! Has few faces perpendicular to axis, do not rotate.
        rotateDirection = -1;
    }
    else
    {
        if (nFaceX > nFaceY && nFaceX > nFaceZ)
        {
            //! Plane is perpendicular to X axis.
            rotateDirection = 1;
        }
        else if (nFaceY > nFaceX && nFaceY > nFaceZ)
        {
            rotateDirection = 2;
        }
        else if (nFaceZ > nFaceX && nFaceZ > nFaceY)
        {
            rotateDirection = 3;
        }
    }

    ostringstream oss;
    oss << "  The axis with maximum parallel faces is: " << rotateDirection << endl;
    WriteLogFile(oss.str());
    PrintToWindow(oss.str());

    return rotateDirection;
}

WallStructure::WallStructure()
{
    nPointPerFace = 0;
    wallFace2Node = 0;
    nodePosition  = 0;
    *box = 0;
    distance = 0.0;
    iBlock = 0;
    xfc.resize(0);
    x.resize(0);
}

WallStructure::~WallStructure()
{
    delete [] nPointPerFace;
    delete [] wallFace2Node;
    delete [] nodePosition;

    for (std::size_t iFace = 0; iFace < face2node.size(); ++ iFace)
    {
        face2node[iFace].clear();
    }
    face2node.clear();
}

void WallStructure::Box()
{
    for (int iDim = 0; iDim < 3; ++ iDim)
    {
        box[iDim]     = LARGE;
        box[iDim + 3] = SMALL;
    }

    for (std::size_t iNode = 0; iNode < this->GetNumberOfWallPoints(); ++ iNode)
    {
        box[0] = MIN(box[0], x[iNode]);
        box[1] = MIN(box[1], y[iNode]);
        box[2] = MIN(box[2], z[iNode]);

        box[3] = MAX(box[3], x[iNode]);
        box[4] = MAX(box[4], y[iNode]);
        box[5] = MAX(box[5], z[iNode]);
    }
}

void WallStructure::Encode(DataContainer *&cdata, const bool &ifMoveToBegin)
{
    if (ifMoveToBegin)
    {
        cdata->MoveToBegin();
    }

    uint_t nSolidFace = GetNumberOfWallFaces();
    PHWrite(cdata, &nSolidFace, 1);

    if (!nSolidFace)
    {
        //! If there are not solid cell, do nothing and return.
        return;
    }

    PHWrite(cdata, &iBlock, 1);

    //! Write the wall face data structure into the data-container.
    //! First: wall face center.
    uint_t nm = 3;
    uint_t nlen = nm * nSolidFace;
    RDouble *qtmp = new RDouble[nlen];
    uint_t countTmp = 0;
    for (uint_t iFace = 0; iFace < nSolidFace; ++ iFace)
    {
        qtmp[countTmp++] = xfc[iFace];
        qtmp[countTmp++] = yfc[iFace];
        qtmp[countTmp++] = zfc[iFace];
    }
    cdata->Write(qtmp, nlen * sizeof(RDouble));
    delete [] qtmp;     qtmp = nullptr;
    countTmp = 0;

    //! Second: points coordinates.
    uint_t nWallPoints = GetNumberOfWallPoints();
    PHWrite(cdata, &nWallPoints, 1);
    if (!nWallPoints)
    {
        return;
    }

    nm = 3;
    nlen = nm * nWallPoints;
    qtmp = new RDouble[nlen];
    for (uint_t iNode = 0; iNode < nWallPoints; ++ iNode)
    {
        qtmp[countTmp++] = x[iNode];
        qtmp[countTmp++] = y[iNode];
        qtmp[countTmp++] = z[iNode];
    }
    cdata->Write(qtmp, nlen * sizeof(RDouble));
    delete [] qtmp;     qtmp = nullptr;
    countTmp = 0;
    //! Third: face to node connection.
    PHWrite(cdata, nPointPerFace, nSolidFace);

    int count = 0;
    for (uint_t iFace = 0; iFace < nSolidFace; ++ iFace)
    {
        count += nPointPerFace[iFace];
    }
    PHWrite(cdata, wallFace2Node, count);

    //! Forth: wall face normal.
    nm = 3;
    nlen = nm * nSolidFace;
    qtmp = new RDouble[nlen];
    for (uint_t iFace = 0; iFace < nSolidFace; ++ iFace)
    {
        qtmp[countTmp++] = xfn[iFace];
        qtmp[countTmp++] = yfn[iFace];
        qtmp[countTmp++] = zfn[iFace];
    }
    cdata->Write(qtmp, nlen * sizeof(RDouble));
    delete [] qtmp;     qtmp = nullptr;
}

void WallStructure::Decode(DataContainer *&cdata, const bool &ifMoveToBegin)
{
    if (ifMoveToBegin)
    {
        cdata->MoveToBegin();
    }

    uint_t nSolidFace;
    PHRead(cdata, &nSolidFace, 1);

    if (nSolidFace == 0)
    {
        //! If there are not solid cell, do nothing and return.
        return;
    }

    PHRead(cdata, &iBlock, 1);

    //! Face center.
    xfc.resize(nSolidFace);
    yfc.resize(nSolidFace);
    zfc.resize(nSolidFace);
    for (uint_t i = 0; i < nSolidFace; ++ i)
    {
        RDouble xc, yc, zc;
        cdata->Read(&xc, sizeof(RDouble));
        cdata->Read(&yc, sizeof(RDouble));
        cdata->Read(&zc, sizeof(RDouble));

        xfc[i] = xc;
        yfc[i] = yc;
        zfc[i] = zc;
    }

    //! Second: points coordinates.
    uint_t nWallPoints;
    PHRead(cdata, &nWallPoints, 1);

    if (nWallPoints == 0)
    {
        return;
    }

    x.resize(nWallPoints);
    y.resize(nWallPoints);
    z.resize(nWallPoints);
    for (uint_t iNode = 0; iNode < nWallPoints; ++ iNode)
    {
        RDouble x0, y0, z0;
        cdata->Read(&x0, sizeof(RDouble));
        cdata->Read(&y0, sizeof(RDouble));
        cdata->Read(&z0, sizeof(RDouble));

        x[iNode] = x0;
        y[iNode] = y0;
        z[iNode] = z0;
    }

    //! Third: face to node connection.
    nPointPerFace = new int [nSolidFace];
    PHRead(cdata, nPointPerFace, nSolidFace);

    int count = 0;
    for (int iFace = 0; iFace < nSolidFace; ++ iFace)
    {
        count += nPointPerFace[iFace];
    }
    wallFace2Node = new int [count];
    PHRead(cdata, wallFace2Node, count);

    //! Analyze node position of each face.
    int *nodePosition = new int [nSolidFace + 1];
    nodePosition[0] = 0;
    int positionCount = 0;
    for (int iFace = 0; iFace < nSolidFace; ++ iFace)
    {
        positionCount += nPointPerFace[iFace];
        nodePosition[iFace + 1] = positionCount;
    }
    this->SetNodePosition(nodePosition);

    //! Forth: wall face normal.
    xfn.resize(nSolidFace);
    yfn.resize(nSolidFace);
    zfn.resize(nSolidFace);
    for (int iFace = 0; iFace < nSolidFace; ++ iFace)
    {
        RDouble xfn0, yfn0, zfn0;
        cdata->Read(&xfn0, sizeof(RDouble));
        cdata->Read(&yfn0, sizeof(RDouble));
        cdata->Read(&zfn0, sizeof(RDouble));

        xfn[iFace] = xfn0;
        yfn[iFace] = yfn0;
        zfn[iFace] = zfn0;
    }

    this->Box();
}

vector<vector<int> > & WallStructure::GetFace2NodeVector()
{
    if (this->face2node.size() != 0)
    {
        return this->face2node;
    }

    ConstructFace2Node(this->GetnPointPerFace(), this->GetWallFace2Node(), this->GetNumberOfWallFaces(), this->face2node);
    
    return this->face2node;
}

void WallStructure::GetFace2NodeVector(vector<int> &face2node, int faceID)
{
    int num = nPointPerFace[faceID];
    face2node.resize(num);
    int count = 0;
    for (int iNode = nodePosition[faceID]; iNode < nodePosition[faceID + 1]; ++ iNode)
    {
        face2node[count ++] = wallFace2Node[iNode];
    }
}

RDouble Pre_WalldistCompute::ComputeDistanceOfANodeToSegment(const RDouble xcc, const RDouble ycc, const RDouble zcc, 
                                                             const RDouble xStart, const RDouble yStart, const RDouble zStart, 
                                                             const RDouble xMid, const RDouble yMid, const RDouble zMid, 
                                                             const RDouble xEnd, const RDouble yEnd, const RDouble zEnd)
{
    RDouble distance;
    RDouble distMid, distStart, distEnd;
    RDouble dx, dy, dz;
    RDouble x1, y1, z1, x2, y2, z2, dist1, dist2;
    RDouble xmin = xStart, ymin = yStart, zmin = zStart;
    RDouble xmax = xEnd, ymax = yEnd, zmax = zEnd;
    RDouble xHalf = xMid, yHalf = yMid, zHalf = zMid;

    dx = xcc - xHalf;
    dy = ycc - yHalf;
    dz = zcc - zHalf;
    distMid = SQR(dx, dy, dz);

    dx = xcc - xmin;
    dy = ycc - ymin;
    dz = zcc - zmin;
    distStart = SQR(dx, dy, dz);

    dx = xcc - xmax;
    dy = ycc - ymax;
    dz = zcc - zmax;
    distEnd = SQR(dx, dy, dz);

    if (distStart < distMid)
    {
        return distStart;
    }

    if (distEnd < distMid)
    {
        return distEnd;
    }

    distance = distMid;
    int iter = 0;
    while (iter < 30)
    {
        ++ iter;
        x1 = (xmin + xHalf) * 0.5;
        y1 = (ymin + yHalf) * 0.5;
        z1 = (zmin + zHalf) * 0.5;
        dx = xcc - x1;
        dy = ycc - y1;
        dz = zcc - z1;
        dist1 = SQR(dx, dy, dz);

        x2 = (xmax + xHalf) * 0.5;
        y2 = (ymax + yHalf) * 0.5;
        z2 = (zmax + zHalf) * 0.5;
        dx = xcc - x2;
        dy = ycc - y2;
        dz = zcc - z2;
        dist2 = SQR(dx, dy, dz);

        bool isInMidlle = false;
        if (distance < dist1 && distance < dist2)
        {
            xmin = x1;
            ymin = y1;
            zmin = z1;
            xmax = x2;
            ymax = y2;
            zmax = z2;
            isInMidlle = true;    //! Do not judge if iteration is converged, if the mid-point is the shortest distance!
        }
        else if (dist1 < dist2)
        {
            xmax = xHalf;
            ymax = yHalf;
            zmax = zHalf;
            xHalf = x1;
            yHalf = y1;
            zHalf = z1;
            distance = dist1;
        }
        else
        {
            xmin = xHalf;
            ymin = yHalf;
            zmin = zHalf;
            xHalf = x2;
            yHalf = y2;
            zHalf = z2;
            distance = dist2;
        }
    }

    return distance;
}

RDouble Pre_WalldistCompute::ComputeDistanceOfANodeToAFace(RDouble *cellCenter, RDouble *faceCenter, RDouble *faceNormal,
                                                           const vector<int> &face2node, value_type &x, value_type &y, value_type &z)
{
    RDouble distance = 0.0;
    RDouble dx = 0.0, dy = 0.0, dz = 0.0;
    dx = cellCenter[0] - faceCenter[0];
    dy = cellCenter[1] - faceCenter[1];
    dz = cellCenter[2] - faceCenter[2];

    //! Judge if the cellCenter at the contrary side of the face normal.
    distance = dx * faceNormal[0] + dy * faceNormal[1] + dz * faceNormal[2];

    distance = fabs(distance);

    //! Projecting the cell center point on to the face plan.
    RDouble x0 = cellCenter[0] + distance * faceNormal[0];
    RDouble y0 = cellCenter[1] + distance * faceNormal[1];
    RDouble z0 = cellCenter[2] + distance * faceNormal[2];
    RDouble projectPoint[3] = {x0, y0, z0};

    //! Judge if the projecting point in one of the sub-triangle or not.
    RDouble coord1[3] = {faceCenter[0], faceCenter[1], faceCenter[2]};
    RDouble oldCoord3[3], coord3[3];
    RDouble oldCoord2[3], coord2[3];

    RDouble distanceTolerance = 1.0e-16;
    RDouble minDistance = LARGE;
    uint_t numberofPointsOfTheFace = face2node.size();
    for (int index = 0; index < numberofPointsOfTheFace; ++ index)
    {
        int p2 = index;
        int p3 = (index + 1) % numberofPointsOfTheFace;

        int node2 = face2node[p2];
        int node3 = face2node[p3];

        oldCoord2[0] = x[node2]; oldCoord2[1] = y[node2]; oldCoord2[2] = z[node2];
        oldCoord3[0] = x[node3]; oldCoord3[1] = y[node3]; oldCoord3[2] = z[node3];
        this->RotateAxis(oldCoord2, coord2);
        this->RotateAxis(oldCoord3, coord3);

        dx = coord2[0] - coord3[0];
        dy = coord2[1] - coord3[1];
        dz = coord2[2] - coord3[2];

        RDouble xn, yn, zn, norm, xx, yy, zz, xa, xb, xc, ya, yb, yc, za, zb, zc;
        RDouble xint, yint, zint, lamda;
        RDouble dotp, dx1, dy1, dz1, dx2, dy2, dz2;
        RDouble disMin = LARGE, dis;
        bool inside = 1;

        xa = coord1[0];
        ya = coord1[1];
        za = coord1[2];
        xb = coord2[0];
        yb = coord2[1];
        zb = coord2[2];
        xc = coord3[0];
        yc = coord3[1];
        zc = coord3[2];

        xx = cellCenter[0];
        yy = cellCenter[1];
        zz = cellCenter[2];

        dx1 = xb - xa;
        dy1 = yb - ya;
        dz1 = zb - za;

        dx2 = xc - xa;
        dy2 = yc - ya;
        dz2 = zc - za;

        xn = dy1 * dz2 - dz1 * dy2;
        yn = dz1 * dx2 - dx1 * dz2;
        zn = dx1 * dy2 - dy1 * dx2;
        norm = sqrt(xn * xn + yn * yn + zn * zn);

        if (norm <= distanceTolerance)
        {
            //! three points co-linear, do something here.
            norm = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
            if (norm <= distanceTolerance)
            {
                norm = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;
                if (norm <= distanceTolerance)
                {
                    //return SQR(cellCenter, coord1);
                    minDistance = MIN(minDistance, SQR(cellCenter, coord1));
                    continue;
                }
                else
                {
                    norm = 1. / sqrt(norm);
                    xn = dx2 * norm; yn = dy2 * norm; zn = dz2 * norm;
                }
            }
            else
            {
                norm = 1. / sqrt(norm);
                xn = dx1 * norm; yn = dy1 * norm; zn = dz1 * norm;
            }

            //! now the distance from the point to the plane.
            RDouble disa, disb, disc;
            disa = xn * (xa - xx) + yn * (ya - yy) + zn * (za - zz);
            disb = xn * (xb - xx) + yn * (yb - yy) + zn * (zb - zz);
            disc = xn * (xc - xx) + yn * (yc - yy) + zn * (zc - zz);
            if (disa * disb <= 0. || disa * disc <= 0.)
            {
                projectPoint[0] = xa - disa * xn;
                projectPoint[1] = ya - disa * yn;
                projectPoint[2] = za - disa * zn;
                minDistance = MIN(minDistance, SQR(cellCenter, projectPoint));
                continue;
            }
            else
            {
                disa = fabs(disa);
                disb = fabs(disb);
                disc = fabs(disc);
                //! all on the same side.
                if (disa > disb)
                {
                    if (disb < disc)
                    {
                        minDistance = MIN(minDistance, SQR(cellCenter, coord2));
                        continue;
                    }
                    else
                    {
                        minDistance = MIN(minDistance, SQR(cellCenter, coord3));
                        continue;
                    }
                }
                else
                {
                    if (disa < disc)
                    {
                        minDistance = MIN(minDistance, SQR(cellCenter, coord1));
                        continue;
                    }
                    else
                    {
                        minDistance = MIN(minDistance, SQR(cellCenter, coord3));
                        continue;
                    }
                }
            }
        }

        norm = 1. / norm;
        xn *= norm;
        yn *= norm;
        zn *= norm;

        dotp = (xa - xx) * xn + (ya - yy) * yn + (za - zz) * zn;

        xint = xx + dotp * xn;
        yint = yy + dotp * yn;
        zint = zz + dotp * zn;

        dx2 = xint - xa;
        dy2 = yint - ya;
        dz2 = zint - za;

        dotp = (dy1 * dz2 - dz1 * dy2) * xn + (dz1 * dx2 - dx1 * dz2) * yn + (dx1 * dy2 - dy1 * dx2) * zn;

        if (dotp < 0.)
        {
            inside = 0;
            lamda = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
            if (lamda >= 0. && lamda <= 1.)
            {
                projectPoint[0] = xa + lamda * dx1;
                projectPoint[1] = ya + lamda * dy1;
                projectPoint[2] = za + lamda * dz1;
                minDistance = MIN(minDistance, SQR(cellCenter, projectPoint));
                continue;
            }
            else if (lamda < 0.)
            {
                dis = ((xa - xx) * (xa - xx) + (ya - yy) * (ya - yy) + (za - zz) * (za - zz));
                if (dis < disMin)
                {
                    projectPoint[0] = coord1[0];
                    projectPoint[1] = coord1[1];
                    projectPoint[2] = coord1[2];
                    disMin = dis;
                }
            }
            else
            {
                dis = ((xb - xx) * (xb - xx) + (yb - yy) * (yb - yy) + (zb - zz) * (zb - zz));
                if (dis < disMin)
                {
                    projectPoint[0] = coord2[0];
                    projectPoint[1] = coord2[1];
                    projectPoint[2] = coord2[2];
                    disMin = dis;
                }
            }
        }

        dx1 = xc - xb;
        dy1 = yc - yb;
        dz1 = zc - zb;

        dx2 = xint - xb;
        dy2 = yint - yb;
        dz2 = zint - zb;

        dotp = (dy1 * dz2 - dz1 * dy2) * xn + (dz1 * dx2 - dx1 * dz2) * yn + (dx1 * dy2 - dy1 * dx2) * zn;

        if (dotp < 0.)
        {
            inside = 0;
            lamda = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
            if (lamda <= 0.) 
            {
                dis = ((xb - xx) * (xb - xx) + (yb - yy) * (yb - yy) + (zb - zz) * (zb - zz));
                if (dis < disMin)
                {
                    projectPoint[0] = coord2[0];
                    projectPoint[1] = coord2[1];
                    projectPoint[2] = coord2[2];
                    disMin = dis;
                }
            }
            else if (lamda >= 1.)
            {
                dis = ((xc - xx) * (xc - xx) + (yc - yy) * (yc - yy) + (zc - zz) * (zc - zz));
                if (dis < disMin)
                {
                    projectPoint[0] = coord3[0];
                    projectPoint[1] = coord3[1];
                    projectPoint[2] = coord3[2];
                    disMin = dis;
                }
            }
            else
            {
                projectPoint[0] = xb + lamda * dx1;
                projectPoint[1] = yb + lamda * dy1;
                projectPoint[2] = zb + lamda * dz1;

                minDistance = MIN(minDistance, SQR(cellCenter, projectPoint));
                continue;
            }
        }

        dx1 = xa - xc;
        dy1 = ya - yc;
        dz1 = za - zc;
        dx2 = xint - xc;
        dy2 = yint - yc;
        dz2 = zint - zc;

        dotp = (dy1 * dz2 - dz1 * dy2) * xn + (dz1 * dx2 - dx1 * dz2) * yn + (dx1 * dy2 - dy1 * dx2) * zn;

        if (dotp < 0.)
        {
            inside = 0;
            lamda = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / (dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
            if (lamda <= 0.)
            {
                dis = ((xc - xx) * (xc - xx) + (yc - yy) * (yc - yy) + (zc - zz) * (zc - zz));
                if (dis < disMin)
                {
                    projectPoint[0] = coord3[0];
                    projectPoint[1] = coord3[1];
                    projectPoint[2] = coord3[2];
                    disMin = dis;
                }
            } 
            else if (lamda >= 1.)
            {
                dis = ((xa - xx) * (xa - xx) + (ya - yy) * (ya - yy) + (za - zz) * (za - zz));
                if (dis < disMin)
                {
                    projectPoint[0] = coord1[0];
                    projectPoint[1] = coord1[1];
                    projectPoint[2] = coord1[2];
                    disMin = dis;
                }
            }
            else
            {
                projectPoint[0] = xc + lamda * dx1;
                projectPoint[1] = yc + lamda * dy1;
                projectPoint[2] = zc + lamda * dz1;

                minDistance = MIN(minDistance, SQR(cellCenter, projectPoint));
                continue;
            }
        }

        if (inside)
        {
            projectPoint[0] = xint;
            projectPoint[1] = yint;
            projectPoint[2] = zint;
        }
        minDistance = MIN(minDistance, SQR(cellCenter, projectPoint));

    }  //! loop over face nodes.

    if (minDistance >= 1.0e30)
    {
        TK_Exit::ExceptionExit("minDistance >= 1.0e30 during wall distance computing!\n");
    }
    return minDistance;
}

void Pre_WalldistCompute::ComputeNodeDistanceByDirectMethod()
{
    int nZones = PHMPI::GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int zoneProcessor = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (myid == zoneProcessor)
        {
            Grid *grid = GetGrid(iZone, this->level);

            UnstructGrid *unstructGrid = UnstructGridCast(grid);

            ComputeNodeDistanceByDirectMethod(unstructGrid);
        }
    }
}

void Pre_WalldistCompute::ComputeNodeDistanceByDirectMethod(UnstructGrid *grid)
{
    int numberOfNodes = grid->GetNTotalNode();

    int blockIndexOfGrid = grid->GetIBlock();

    RDouble *nodeCoordinateX = grid->GetX();
    RDouble *nodeCoordinateY = grid->GetY();
    RDouble *nodeCoordinateZ = grid->GetZ();

    vector<RDouble> &nodeWalldistDistance = grid->GetNodeWallDist();

    nodeWalldistDistance.resize(numberOfNodes);

    SetField(nodeWalldistDistance, LARGE, numberOfNodes);

    int numberOfWallStructureInPresentBlock = 0;

    for (unsigned int iWall = 0; iWall < wallstructureList->size(); ++ iWall)
    {
        WallStructure *wallStructure = (*wallstructureList)[iWall];

        int blockIndexOfWall = wallStructure->GetIBlock();

        if (blockIndexOfWall != blockIndexOfGrid)
        {
            continue;
        }

        numberOfWallStructureInPresentBlock ++;

        WallStructure::value_type &wallNodeCoordinateX = wallStructure->GetX();
        WallStructure::value_type &wallNodeCoordinateY = wallStructure->GetY();
        WallStructure::value_type &wallNodeCoordinateZ = wallStructure->GetZ();

        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            for (unsigned int iWallNode = 0; iWallNode < wallStructure->GetNumberOfWallPoints(); ++ iWallNode)
            {
                RDouble dx = nodeCoordinateX[iNode] - wallNodeCoordinateX[iWallNode];
                RDouble dy = nodeCoordinateY[iNode] - wallNodeCoordinateY[iWallNode];
                RDouble dz = nodeCoordinateZ[iNode] - wallNodeCoordinateZ[iWallNode];
                RDouble dis = DISTANCE(dx, dy, dz);

                nodeWalldistDistance[iNode] = MIN(dis , nodeWalldistDistance[iNode]);
            }
        }
    }

    if (numberOfWallStructureInPresentBlock == 0)
    {
        RDouble dist = GlobalDataBase::GetDoubleParaFromDB("walldistMainZone");
        RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

        SetField(nodeWalldistDistance, dist / gridScaleFactor, numberOfNodes);

        return;
    }
}

void Pre_WalldistCompute::ComputeNodeDistanceForSphere(UnstructGrid *grid)
{
    int numberOfNodes = grid->GetNTotalNode();

    int blockIndexOfGrid = grid->GetIBlock();

    RDouble *nodeCoordinateX = grid->GetX();
    RDouble *nodeCoordinateY = grid->GetY();
    RDouble *nodeCoordinateZ = grid->GetZ();

    vector<RDouble> &nodeWalldistDistance = grid->GetNodeWallDist();

    nodeWalldistDistance.resize(numberOfNodes);

    SetField(nodeWalldistDistance, LARGE, numberOfNodes);

    int numberOfWallStructureInPresentBlock = 0;

    for (unsigned int iWall = 0; iWall < wallstructureList->size(); ++ iWall)
    {
        WallStructure *wallStructure = (*wallstructureList)[iWall];

        int blockIndexOfWall = wallStructure->GetIBlock();

        if (blockIndexOfWall != blockIndexOfGrid)
        {
            continue;
        }

        numberOfWallStructureInPresentBlock ++;
    }

    if (numberOfWallStructureInPresentBlock == 0)
    {
        RDouble dist = GlobalDataBase::GetDoubleParaFromDB("walldistMainZone");
        RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

        SetField(nodeWalldistDistance, dist / gridScaleFactor, numberOfNodes);

        return;
    }

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        RDouble xc = 0.0;
        RDouble yc = 0.0;
        RDouble zc = 0.0;

        if (blockIndexOfGrid == 2)
        {
            yc += 1.4;
        }
        else if (blockIndexOfGrid == 3)
        {
            yc -= 1.4;
        }
        else if (blockIndexOfGrid == 4)
        {
            xc -= 1.4;
        }
        else if (blockIndexOfGrid == 5)
        {
            xc += 1.4;
        }

        RDouble dx = nodeCoordinateX[iNode] - xc;
        RDouble dy = nodeCoordinateY[iNode] - yc;
        RDouble dz = nodeCoordinateZ[iNode] - zc;

        nodeWalldistDistance[iNode] = DISTANCE(dx, dy, dz) - 0.5;
    }
}

void Pre_WalldistCompute::ComputeNodeDistanceByADTMethod()
{
    TimeSpan *timeSpan = new TimeSpan();
    using namespace PHMPI;

    int nLocalZones = GetNumberofLocalZones();

    //! Set initial value.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        UnstructGrid *grid = UnstructGridCast(GetGrid(zoneID, this->level));

        vector<RDouble> &nodeWallDis = grid->GetNodeWallDist();

        int numberOfNodes = grid->GetNTotalNode();

        nodeWallDis.resize(numberOfNodes);

        SetField(&nodeWallDis[0], LARGE, numberOfNodes);
    }

    int nBlocks = GlobalDataBase::GetIntParaFromDB("numberOfGridGroups");

    for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
    {
        PrintToWindow("\n    Computing Wall dist for block: ", iBlock, "\n\n");
        ComputeNodeDistanceByADTMethodForBlock(iBlock);
    }
    timeSpan->ShowTimeSpanToLogFile("Compute Node Distance By ADT Method");
    delete timeSpan;
}

void Pre_WalldistCompute::ComputeNodeDistanceByKDTreeMethod()
{
    TimeSpan *timeSpan = new TimeSpan();

    using namespace PHMPI;

    int nLocalZones = GetNumberofLocalZones();

    //! Set initial value.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        UnstructGrid * grid = UnstructGridCast(GetGrid(zoneID, this->level));

        vector<RDouble> & nodeWallDis = grid->GetNodeWallDist();

        int numberOfNodes = grid->GetNTotalNode();

        nodeWallDis.resize(numberOfNodes);

        SetField(&nodeWallDis[0], LARGE, numberOfNodes);
    }

    ComputeNodeDistanceByKDTreeMethodForBlock();

    timeSpan->ShowTimeSpanToLogFile("Compute Node Distance By KDTree Method");
    delete timeSpan;
}

void Pre_WalldistCompute::ComputeNodeDistanceByKDTreeMethodInOtherBlock()
{
    using namespace PHMPI;

    int nLocalZones = GetNumberofLocalZones();

    //! Set initial value.
    for (int izone = 0; izone < nLocalZones; ++ izone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(izone);

        UnstructGrid * grid = UnstructGridCast(GetGrid(zoneID, this->level));

        vector<RDouble> &nodeWallDistanceInOtherBlock = grid->GetNodeWallDistInOtherBlock();

        int numberOfNodes = grid->GetNTotalNode();

        nodeWallDistanceInOtherBlock.resize(numberOfNodes);

        SetField(&nodeWallDistanceInOtherBlock[0], LARGE, numberOfNodes);
    }

    ComputeNodeDistanceByKDTreeMethodForBlockInOtherBlock();
}

void Pre_WalldistCompute::ComputeNodeDistanceByKDTreeMethodForBlock()
{
    using namespace PHMPI;

    int nDim = 3;

    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        UnstructGrid *grid = UnstructGridCast(GetGrid(zoneID, this->level));
        int iBlock = grid->GetIBlock();
        this->presentBlock = iBlock;

        if (!wallNodeKDTree)
        {
            wallNodeKDTree = CreatKDTree(nDim);
        }
        BuildWallStructKDTreeNode();

        ostringstream oss;
        oss << "  iBlock = " << iBlock << " nTWallFace = " << nTWallFace;
        WriteLogFile(oss);

        vector<RDouble> &nodeWalldist = grid->GetNodeWallDist();

        int nTotalNode = grid->GetNTotalNode();

        //! if nTWallFace = 0, this must be the main back ground block
        if (nTWallFace == 0)
        {
            RDouble defaultBlockWalldist = PHSPACE::GetRDoubleParameterFromDataBase(iBlock,"defaultBlockWalldist");
            RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

            SetField(nodeWalldist, defaultBlockWalldist / gridScaleFactor, nTotalNode);

            continue;
        }

        RDouble *x = grid->GetX();
        RDouble *y = grid->GetY();
        RDouble *z = grid->GetZ();

        PrintToWindow("  Wall distance computing for zone: ", iZone, "\n");
        WriteLogFile("  Wall distance computing for zone: ", iZone);

        ComputeNodeWallDist3DByKDTreeMethod(nTotalNode, x, y, z, &nodeWalldist[0]);

        if (wallNodeKDTree)
        {
            delete wallNodeKDTree;
            wallNodeKDTree = 0;
        }
    }
}

void Pre_WalldistCompute::ComputeNodeDistanceByKDTreeMethodForBlockInOtherBlock()
{
    using namespace PHMPI;

    int nDim = 3;

    int nLocalZones = GetNumberofLocalZones();

    for (int izone = 0; izone < nLocalZones; ++ izone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(izone);

        UnstructGrid *grid = UnstructGridCast(GetGrid(zoneID, this->level));
        presentBlock = grid->GetIBlock();

        wallNodeKDTree = CreatKDTree(nDim);
        BuildWallStructKDTreeNodeInOtherBlock();

        vector<RDouble> &nodeWalldist = grid->GetNodeWallDistInOtherBlock();

        int nTotalNode = grid->GetNTotalNode();

        RDouble *x = grid->GetX();
        RDouble *y = grid->GetY();
        RDouble *z = grid->GetZ();

        ComputeNodeWallDist3DByKDTreeMethod(nTotalNode, x, y, z, &nodeWalldist[0]);

        if (wallNodeKDTree)
        {
            delete wallNodeKDTree;
            wallNodeKDTree = 0;
        }
    }
}

void Pre_WalldistCompute::ComputeNodeDistanceByADTMethodForBlock(int iBlock)
{
    using namespace PHMPI;

    int nLocalZones = GetNumberofLocalZones();

    this->presentBlock = iBlock;

    this->rotateAxis = 0;

    while (!IsNodeWalldistComputOver())
    {
        PrintToWindow("Wall dist computing for axis: ", this->rotateAxis, "\n");

        BuildWallStructTree();

        ostringstream oss;
        oss << " iBlock = " << iBlock << " nTWallFace = " << nTWallFace << "\n";
        WriteLogFile(oss);

        for (int iZone = 0; iZone < nLocalZones; ++ iZone)
        {
            int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

            UnstructGrid *grid = UnstructGridCast(GetGrid(zoneID, this->level));

            if (grid->GetIBlock()!= presentBlock)
            {
                continue;
            }

            vector<RDouble> &nodeWalldist = grid->GetNodeWallDist();

            int nTotalNode = grid->GetNTotalNode();

            //! if nTWallFace = 0, this must be the main back ground block
            if (nTWallFace == 0)
            {
                RDouble dist = GlobalDataBase::GetDoubleParaFromDB("walldistMainZone");
                RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

                SetField(nodeWalldist, dist / gridScaleFactor, nTotalNode);

                continue;
            }
            
            RDouble *x = grid->GetX();
            RDouble *y = grid->GetY();
            RDouble *z = grid->GetZ();

            ComputeNodeWallDist3DByProjectingFast(nTotalNode, x, y, z, &nodeWalldist[0]);
        }
        
        ++ this->rotateAxis;

        delete wallFaceTree;

        wallFaceTree = 0;

        if (this->rotateAxis == 5) break;
    }
}

LIB_EXPORT void Pre_WalldistCompute::ComputeWallDistForNodesList(int blockIndex, int nTotalNode, RDouble *x,  RDouble *y, RDouble *z,  RDouble *dist)
{
    if (walldistComputeMethod == FAST_METHOD || walldistComputeMethod == SUPER_FAST_METHOD)
    {
        ComputeWallDistForNodesListByADT(blockIndex, nTotalNode, x, y, z, dist);
    }
    else
    {
        //ComputeWallDistForNodesListByDirectMethod(blockIndex, nTotalNode, x, y, z, dist);
        ComputeWallDistForNodesListByKDT(blockIndex, nTotalNode, x, y, z, dist);
    }
}

void Pre_WalldistCompute::ComputeWallDistForNodesListByDirectMethod(int blockIndex, int nTotalNode, RDouble *x,  RDouble *y, RDouble *z,  RDouble *dist)
{
    uint_t numberOfWalls = wallstructureList->size();

    int numberOfWallStructureInPresentBlock = 0;

    for (int iWall = 0; iWall < numberOfWalls; ++ iWall)
    {
        int iBlockOfWall = (*wallstructureList)[iWall]->GetIBlock();

        if (iBlockOfWall != blockIndex)
        {
            continue;
        }

        numberOfWallStructureInPresentBlock ++;

        vector<RDouble> wallFaceCoordinateX = (*wallstructureList)[iWall]->GetXFaceCenter();
        vector<RDouble> wallFaceCoordinateY = (*wallstructureList)[iWall]->GetYFaceCenter();
        vector<RDouble> wallFaceCoordinateZ = (*wallstructureList)[iWall]->GetZFaceCenter();

        uint_t numberOfWallFace = (*wallstructureList)[iWall]->GetNumberOfWallFaces();

        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            for (int iWallFace = 0; iWallFace < numberOfWallFace; ++ iWallFace)
            {
                RDouble dx = x[iNode] - wallFaceCoordinateX[iWallFace];
                RDouble dy = y[iNode] - wallFaceCoordinateY[iWallFace];
                RDouble dz = z[iNode] - wallFaceCoordinateZ[iWallFace];
                RDouble dis = DISTANCE(dx, dy, dz);

                dist[iNode] = MIN(dis, dist[iNode]);
            }
        }
    }

    if (numberOfWallStructureInPresentBlock == 0)
    {
        RDouble distByfile = GlobalDataBase::GetDoubleParaFromDB("walldistMainZone");
        RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

        SetField(dist, distByfile / gridScaleFactor, nTotalNode);
    }
}

void Pre_WalldistCompute::ComputeWallDistForNodesListForSphere(int blockIndex, int nTotalNode, RDouble *x,  RDouble *y, RDouble *z,  RDouble *dist)
{
    uint_t numberOfWalls = wallstructureList->size();

    int numberOfWallStructureInPresentBlock = 0;

    for (int iWall = 0; iWall < numberOfWalls; ++ iWall)
    {
        int iBlockOfWall = (*wallstructureList)[iWall]->GetIBlock();

        if (iBlockOfWall != blockIndex)
        {
            continue;
        }

        numberOfWallStructureInPresentBlock ++;
    }

    if (numberOfWallStructureInPresentBlock == 0)
    {
        RDouble distByfile = GlobalDataBase::GetDoubleParaFromDB("walldistMainZone");
        RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");
        SetField(dist, distByfile / gridScaleFactor, nTotalNode);
        return;
    }

    RDouble minDist = LARGE;

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        RDouble xc = 0.0;
        RDouble yc = 0.0;
        RDouble zc = 0.0;

        if (blockIndex == 2)
        {
            yc += 1.4;
        }
        else if (blockIndex == 3)
        {
            yc -= 1.4;
        }
        else if (blockIndex == 4)
        {
            xc -= 1.4;
        }
        else if (blockIndex == 5)
        {
            xc += 1.4;
        }

        RDouble dx = x[iNode] - xc;
        RDouble dy = y[iNode] - yc;
        RDouble dz = z[iNode] - zc;

        dist[iNode] = DISTANCE(dx, dy, dz) - 0.5;

        minDist = MIN(minDist, dist[iNode]);
    }

    ostringstream oss;
    oss << "For block : " << blockIndex << "Min dist = " << minDist << endl;
    WriteLogFile(oss.str());
}

void Pre_WalldistCompute::ComputeWallDistForNodesListByADT(int blockIndex, int nTotalNode, RDouble *x,  RDouble *y, RDouble *z,  RDouble *dist)
{
    using namespace PHMPI;

    this->presentBlock = blockIndex;

    this->rotateAxis = 0;

    bool cycle = true;

    while (cycle)
    {
        PrintToWindow("Wall dist computing for axis: ", this->rotateAxis, "\n");

        BuildWallStructTree();

        //! If nTWallFace = 0, this must be the main back ground block.
        if (nTWallFace == 0)
        {
            RDouble distByFile = GlobalDataBase::GetDoubleParaFromDB("walldistMainZone");
            RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

            SetField(dist, distByFile / gridScaleFactor, nTotalNode);

            break;
        }
        
        ComputeNodeWallDist3DByProjectingFast(nTotalNode, x, y, z, dist);
        
        ++ this->rotateAxis;

        delete wallFaceTree;

        wallFaceTree = 0;

        if (this->rotateAxis == 5) break;

        cycle = false;

        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            if (dist[iNode] >= LARGE * 0.9)
            {
                cycle = true;
                break;
            }
        }
    }
}

void Pre_WalldistCompute::ComputeWallDistForNodesListByKDT(int blockIndex, int nTotalNode, RDouble *x, RDouble *y, RDouble *z, RDouble *dist)
{
    using namespace PHMPI;

    int nDim = 3;

    this->presentBlock = blockIndex;

    bool cycle = true;

    while (cycle)
    {
        if (!wallNodeKDTree)
        {
            wallNodeKDTree = CreatKDTree(nDim);
        }
        BuildWallStructKDTreeNode();

        //! If nTWallFace = 0, this must be the main back ground block.
        if (nTWallFace == 0)
        {
            RDouble distByFile = GlobalDataBase::GetDoubleParaFromDB("walldistMainZone");
            RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

            SetField(dist, distByFile / gridScaleFactor, nTotalNode);

            break;
        }

        ComputeNodeWallDist3DByKDTreeMethod(nTotalNode, x, y, z, dist);

        if (wallNodeKDTree)
        {
            delete wallNodeKDTree;
            wallNodeKDTree = 0;
        }

        cycle = false;

        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            if (dist[iNode] >= LARGE * 0.9)
            {
                cycle = true;
                break;
            }
        }
    }
}

vector<WallStructure *> * Pre_WalldistCompute::GetWallstructureList()
{
    return wallstructureList;
};


void GenerateWalldist()
{
    if (GetTaskCode() == HO_SOLVER) return;

    if (PHMPI::CurrentProcessorIsGridProcessor())
    {
        return;
    }

    int viscousType = 0;
    GlobalDataBase::GetData("viscousType", &viscousType, PHINT, 1);
    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
    int nIsComputeWallDist = GlobalDataBase::GetIntParaFromDB("nIsComputeWallDist");    //! just for laminar or Euler flow.  by dms

    int sysGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if (viscousType > LAMINAR || sysGridType != STRUCTGRID)    //! turbulence flow or not STRUCTGRID
    {
        nIsComputeWallDist = 0;
    }

    if (GetTaskCode() == CAL_WALL_DIST || (viscousType != INVISCID && IsReadWalldist() == false))
    {
        if(nIsComputeWallDist == 0)
        {
            AllocateWalldist();

            ComputeWalldist();

            WriteWalldist();
        }
        //VisualWalldist();
    }

    if (viscousType > LAMINAR || iLES == LES_SOLVER)
    {
        WalldistStatics();
    }

    int flowInitMethod = GlobalDataBase::GetIntParaFromDB("flowInitMethod");
    if (viscousType != INVISCID && flowInitMethod > 0)
    {
        CompMaxBoundaryLayerWalldist();
    }
#ifdef AI_RandomForestRegressor    
    DumpMaxBoundaryLayerWalldist();
#endif
}

#ifdef AI_RandomForestRegressor
void DumpMaxBoundaryLayerWalldist()
{
    ofstream dataFile;
    dataFile.open("./grid/DistAndNormal.txt", ofstream::app);
    RDouble globalMaxBoundaryLayerWalldist;
    GlobalDataBase::GetData("maxBoundaryLayerWalldist", &globalMaxBoundaryLayerWalldist, PHDOUBLE, 1);
    dataFile << globalMaxBoundaryLayerWalldist << endl;
    dataFile.close();
}
#endif

void CompMaxBoundaryLayerWalldist()
{
    RDouble globalMaxBoundaryLayerWalldist = -1.0;
    if (!IsBoundaryLayerExist())
    {
        GlobalDataBase::UpdateData("maxBoundaryLayerWalldist", &globalMaxBoundaryLayerWalldist, PHDOUBLE, 1);
        return;
    }

    int dimension = GetDim();
    RDouble maxBoundaryLayerWalldist = - LARGE;

    int myid   = PHMPI::GetCurrentProcessorID();
    int nZones = PHMPI::GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int procID = PHMPI::GetZoneProcessorID(iZone);
        if (myid == procID)
        {
            Grid *grid0 = GetGrid(iZone, 0);
            if (grid0->Type() == STRUCTGRID)
            {
                continue;
            }

            UnstructGrid *grid = UnstructGridCast(grid0);
            int nTotalCell = grid->GetNTotalCell();
            int *nodeNumOfEachCell = grid->GetNodeNumberOfEachCell();
            RDouble *allWallDist = grid->GetWallDist();

            for (int iCell = 0; iCell < nTotalCell; ++ iCell)
            {
                int nodeNum = nodeNumOfEachCell[iCell];
                int cellType = GetCellType(dimension, nodeNum);
                if ((dimension == THREE_D && (cellType == TETRA_4 || cellType == PYRA_5))
                 || (dimension == TWO_D   && cellType == TRI_3))
                {
                    continue;
                }

                RDouble dist = allWallDist[iCell];
                if (dist > maxBoundaryLayerWalldist)
                {
                    maxBoundaryLayerWalldist = dist;
                }
            }
        }
    }

    PH_AllReduce(&maxBoundaryLayerWalldist, &globalMaxBoundaryLayerWalldist, 1, MPI_MAX);
    GlobalDataBase::UpdateData("maxBoundaryLayerWalldist", &globalMaxBoundaryLayerWalldist, PHDOUBLE, 1);
}

bool IsBoundaryLayerExist()
{
    int nZones = PHMPI::GetNumberofGlobalZones();
    int myid   = PHMPI::GetCurrentProcessorID();
    int dimension = GetDim();
    int hasCellType = 1;

    set < int > cellTypeSet;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int procID = PHMPI::GetZoneProcessorID(iZone);
        if (myid == procID)
        {
            Grid *grid0 = GetGrid(iZone, 0);
            if (grid0->Type() == STRUCTGRID)
            {
                continue;
            }

            UnstructGrid *grid = UnstructGridCast(grid0);
            int nTotalCell = grid->GetNTotalCell();
            int *nodeNumOfEachCell = grid->GetNodeNumberOfEachCell();

            for (int iCell = 0; iCell < nTotalCell; ++ iCell)
            {
                int nodeNum = nodeNumOfEachCell[iCell];
                int cellType = GetCellType(dimension, nodeNum);
                if (cellType == -1)
                {
                    hasCellType = 0;
                }
                cellTypeSet.insert(cellType);
            }
        }
    }

    int globalHasCellType = 1;
    PH_AllReduce(&hasCellType, &globalHasCellType, 1, MPI_MIN);
   
    if (globalHasCellType  == 0)
    {
        PrintToWindow("This cell type is not supported!\n");
        return false;
    }

    int localCellTypeSize = static_cast<int>(cellTypeSet.size());
    int globalCellTypeSize = 0;
    PH_AllReduce(&localCellTypeSize, &globalCellTypeSize, 1, MPI_MAX);

    if (globalCellTypeSize > 1)
    {
        return true;
    }

    return false;
}

bool IsReadWalldist()
{
    using namespace PHMPI;

    int nZones = GetNumberofGlobalZones();
    int myid   = GetCurrentProcessorID();

    int isWDTExist = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int procID = GetZoneProcessorID(iZone);

        if (!CurrentProcessorIsGridProcessor())
        {
            if (myid == procID)
            {
                Grid *grid0 = GetGrid(iZone, 0);
                if (grid0->Type() == UNSTRUCTGRID)
                {
                    UnstructGrid *grid = UnstructGridCast(grid0);
                    if (grid->GetWallDist()) isWDTExist = 1;
                }
                else
                {
                    StructGrid *grid = StructGridCast(grid0);
                    if (grid->GetWallDist()) isWDTExist = 1;
                }
            }
        }
    }

    PH_CompareMaxMin(isWDTExist, 1);

    if (isWDTExist == 1)
    {
        return true;
    }

    return false;
}

void ComputeWalldist()
{
    int walldistMethod;
    GlobalDataBase::GetData("walldistMethod", &walldistMethod, PHINT, 1);

    int level = 0;

    Pre_WalldistCompute walldistCompute(walldistMethod, level);

    walldistCompute.Run();
}

void AllocateWalldist()
{
    ActionKey *actkey = new ActionKey();

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    int myid = GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int proc = GetZoneProcessorIDSepMode(iZone);

        if (myid == proc)
        {
            Grid *grid = GetGrid(iZone, actkey->level);
            grid->AllocateWalldist();
        }
    }
    delete actkey;    actkey = nullptr;
}

#ifdef USE_TecplotLib
void VisualWalldist()
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    ostringstream oss;
    oss << "x y z dist";

    string Variable = oss.str();
    const char *Variables = Variable.c_str();

    string filename = "./in/WallDistVisual.plt";
    const char *pltname = filename.c_str();

    INTEGER4 Debug     = 1;
    INTEGER4 VIsDouble = 0;
    INTEGER4 FileType  = 0;
    INTEGER4 I;

    I = TECINI112((char*)"PHengLEI Grid Videotex",
                  (char*)Variables,
                  (char*)pltname,
                  (char*)".",
                  &FileType,
                  &Debug,
                  &VIsDouble);

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        UnstructGrid *grid = UnstructGridCast(GetGrid(zoneID, 0));

        INTEGER4 IMxOrNumPts, JMxOrNumElements, KMxOrNumFaces;
        INTEGER4 ZoneType;

        IMxOrNumPts = grid->GetNTotalNode();
        JMxOrNumElements = grid->GetNTotalCell();
        KMxOrNumFaces = grid->GetNTotalFace();

        int dimension = grid->GetDim();

        if (dimension == 2)
        {
            ZoneType = 6;
        }
        else
        {
            ZoneType = 7;
        }

        int GridID = grid->GetGridID()->GetIndex();

        RDouble *x = grid->GetX();
        RDouble *y = grid->GetY();
        RDouble *z = grid->GetZ();

        INTEGER4 ValueLocation[4] = {1, 1, 1, 0};

        INTEGER4 ICellMax           = 0;
        INTEGER4 JCellMax           = 0;
        INTEGER4 KCellMax           = 0;
        double   SolutionTime       = 0.0;
        INTEGER4 StrandID           = 0;
        INTEGER4 ParentZone         = 0;
        INTEGER4 IsBlock            = 1;
        INTEGER4 NumFaceConnections = 0;
        INTEGER4 FaceNeighborMode   = 0;

        INTEGER4 TotalNumFaceNodes_Rect = 0;

        //int * cell2node = grid->GetCell2Node();

        int *faceNodeNumberContainer = grid->GetNodeNumberOfEachFace();
        int *faceNodeIndexContainer  = grid->GetFace2Node();
    
        for (int iFace = 0; iFace < KMxOrNumFaces; ++ iFace)
        {
            TotalNumFaceNodes_Rect += faceNodeNumberContainer[iFace];
        }

        INTEGER4 *FaceNodes_Rect = new INTEGER4[TotalNumFaceNodes_Rect];
        for (INTEGER4 i = 0; i < TotalNumFaceNodes_Rect; ++ i)
        {
            FaceNodes_Rect[i] = faceNodeIndexContainer[i] + 1;
        }

        INTEGER4 NumConnBndryFaces_Rect  = 0;    /* interface */
        INTEGER4 TotalNumBndryConns_Rect = 0;
        INTEGER4 SharConn                = 0;

        ostringstream oss;

        //Zonetitle = Zone + zonenum;
        oss << "Zone" << GridID;
    
        //const char * ZoneTitle = Zonetitle.c_str();
        string Zonetitle = oss.str();
        const char *ZoneTitle = Zonetitle.c_str();

        I = TECZNE112((char*)ZoneTitle,
                          &ZoneType,
                          &IMxOrNumPts,
                          &JMxOrNumElements,
                          &KMxOrNumFaces,
                          &ICellMax,
                          &JCellMax,
                          &KCellMax,
                          &SolutionTime,
                          &StrandID,
                          &ParentZone,
                          &IsBlock,
                          &NumFaceConnections,
                          &FaceNeighborMode,
                          &TotalNumFaceNodes_Rect,
                          &NumConnBndryFaces_Rect,
                          &TotalNumBndryConns_Rect,
                          NULL,
                          ValueLocation,
                          NULL,
                          &SharConn);

        INTEGER4 IsDouble = 1;
        I = TECDAT112(&IMxOrNumPts, x, &IsDouble);
        I = TECDAT112(&IMxOrNumPts, y, &IsDouble);
        I = TECDAT112(&IMxOrNumPts, z, &IsDouble);

        RDouble *walldist = grid->GetWallDist();
        I = TECDAT112(&JMxOrNumElements, walldist, &IsDouble);

        INTEGER4 *FaceLeftElems  = new INTEGER4[KMxOrNumFaces];
        INTEGER4 *FaceRightElems = new INTEGER4[KMxOrNumFaces];

        int *leftCellIndexContainer  = grid->GetLeftCellOfFace();
        int *rightCellIndexContainer = grid->GetRightCellOfFace();

        for (INTEGER4 iFace = 0; iFace < KMxOrNumFaces; ++ iFace)
        {
            FaceLeftElems [iFace] = leftCellIndexContainer [iFace] + 1;
            FaceRightElems[iFace] = rightCellIndexContainer[iFace] + 1;
            if (FaceRightElems[iFace] > JMxOrNumElements || FaceRightElems[iFace] < 0) FaceRightElems[iFace] = 0;
            if (FaceLeftElems[iFace] > JMxOrNumElements || FaceLeftElems[iFace] < 0) FaceLeftElems[iFace] = 0;
        }

        INTEGER4 FaceBndryConnectionCounts = 0;
                 INTEGER4 FaceBndryConnectionElems = 0;
                 INTEGER4 FaceBndryConnectionZones = 0;

        I = TECPOLY112(faceNodeNumberContainer,
                       FaceNodes_Rect,
                       FaceLeftElems,
                       FaceRightElems,
                       &FaceBndryConnectionCounts,
                       &FaceBndryConnectionElems,
                       &FaceBndryConnectionZones);

        delete [] FaceLeftElems;
        delete [] FaceRightElems;
        delete [] FaceNodes_Rect;
    }

    I = TECEND112();
}
#endif

void WriteWalldist(int dumpData)
{
    string walldistfile = "walldist.dat";
    GlobalDataBase::GetData("walldistfile", &walldistfile, PHSTRING, 1);

    ActionKey *actkey = new ActionKey();
    actkey->filename = (walldistfile);
    actkey->subact   = dumpData;

    ParallelCreateHDF5File(actkey);

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int send_proc = GetZoneProcessorID(iZone);
        int recv_proc = GetZoneFileID(iZone);

        int myid = GetCurrentProcessorID();
        int tag = GetSendRecvTag(actkey, iZone);

        if (myid == send_proc)
        {
            Grid *grid = GetGrid(iZone, 0);
            grid->DumpWalldist(actkey);
        }

        PH_Trade(actkey, send_proc, recv_proc, tag);

        if (myid == recv_proc)
        {
            DumpWalldistToFile(actkey, iZone);
        }
    }

    ParallelCloseHDF5File(actkey);
    delete actkey;    actkey = nullptr;
}

void DumpWalldistToFile(ActionKey *actkey, int iZone)
{
    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();

    int nTotalCell = 0;
    cdata->Read(&nTotalCell, sizeof(int));

    RDouble *walldist = new RDouble[nTotalCell];
    cdata->Read(walldist, sizeof(RDouble) * nTotalCell);
    RDouble *nearestwallfacenormalx = new RDouble[nTotalCell];
    cdata->Read(nearestwallfacenormalx, sizeof(RDouble) * nTotalCell);
    RDouble *nearestwallfacenormaly = new RDouble[nTotalCell];
    cdata->Read(nearestwallfacenormaly, sizeof(RDouble) * nTotalCell);
    RDouble *nearestwallfacenormalz = new RDouble[nTotalCell];
    cdata->Read(nearestwallfacenormalz, sizeof(RDouble) * nTotalCell);

    ostringstream oss;
    oss << "Grid-" << iZone;
    string dataName = oss.str();

    double *wdt = new double[nTotalCell];
    double *nwfnx = new double[nTotalCell];
    double *nwfny = new double[nTotalCell];
    double *nwfnz = new double[nTotalCell];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        wdt[iCell] = static_cast<double>(walldist[iCell]);
        nwfnx[iCell] = static_cast<double>(nearestwallfacenormalx[iCell]);
        nwfny[iCell] = static_cast<double>(nearestwallfacenormaly[iCell]);
        nwfnz[iCell] = static_cast<double>(nearestwallfacenormalz[iCell]);
    }

    string nearestwallfacenormalxname = dataName + "_nwfnx";
    string nearestwallfacenormalyname = dataName + "_nwfny";
    string nearestwallfacenormalzname = dataName + "_nwfnz";
    CreateAndWriteData(actkey->filepos, dataName, 1, nTotalCell, PHDOUBLE, wdt);

    if (DUMPAllWALLDISTINFO == actkey->subact)
    {
        CreateAndWriteData(actkey->filepos, nearestwallfacenormalxname, 1, nTotalCell, PHDOUBLE, nwfnx);
        CreateAndWriteData(actkey->filepos, nearestwallfacenormalyname, 1, nTotalCell, PHDOUBLE, nwfny);
        CreateAndWriteData(actkey->filepos, nearestwallfacenormalzname, 1, nTotalCell, PHDOUBLE, nwfnz);
    }

    delete [] walldist;    walldist = nullptr;
    delete [] wdt;    wdt = nullptr;
    delete [] nearestwallfacenormalx;    nearestwallfacenormalx = nullptr;
    delete [] nearestwallfacenormaly;    nearestwallfacenormaly = nullptr;
    delete [] nearestwallfacenormalz;    nearestwallfacenormalz = nullptr;
    delete [] nwfnx;    nwfnx = nullptr;
    delete [] nwfny;    nwfny = nullptr;
    delete [] nwfnz;    nwfnz = nullptr;
}

void WalldistStatics()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int myid = GetCurrentProcessorID();

    RDouble minWalldist = PHSPACE::LARGE;
    RDouble maxWalldist = PHSPACE::SMALL;
    
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int proc = GetZoneProcessorIDSepMode(iZone);

        if (myid == proc)
        {
            Grid *grid = GetGrid(iZone, 0);

            if (grid->Type() == UNSTRUCTGRID)
            {
                int nTotalCell = grid->GetNTotalCell();
                RDouble *walldist = UnstructGridCast(grid)->GetWallDist();
                for (int iCell = 0; iCell < nTotalCell; ++ iCell)
                {
                    minWalldist = MIN(minWalldist, walldist[iCell]);
                    maxWalldist = MAX(maxWalldist, walldist[iCell]);
                    if (maxWalldist > 1.0e8)
                    {
                        cout << iCell << endl;
                        ostringstream oss1;
                        oss1 << "iCell=" << iCell << endl;
                        TK_Exit::ExceptionExit(oss1);
                    }
                }
            }
            else
            {
                StructGrid *structGrid = StructGridCast(grid);

                RDouble3D &walldist = *structGrid->GetWallDist();

                int ist, ied, jst, jed, kst, ked;
                structGrid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            minWalldist = MIN(minWalldist, walldist(i, j, k));
                            maxWalldist = MAX(maxWalldist, walldist(i, j, k));
                        }
                    }
                }
            }
        }
    }

    PH_CompareMaxMin(minWalldist, 2);

    PH_CompareMaxMin(maxWalldist, 1);

    ostringstream oss;
    oss << "Min && Max wall distance : " << minWalldist << ", " << maxWalldist << endl;
    PrintToWindow(oss);
}

}