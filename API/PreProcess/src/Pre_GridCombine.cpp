#include "Pre_GridCombine.h"
#include "Math_BasisFunction.h"
#include "Glb_Dimension.h"
#include "IO_FileName.h"
#include "Region.h"
#include "OverLappingGrid.h"
#include "Pre_GridConversion.h"
#include "PHIO.h"
#include "TK_Log.h"

using namespace std;

namespace PHSPACE
{

CombinGrid::CombinGrid(Region *region, Grid **grid_in, int nZones_in)
{
    this->region = region;
    this->grid_in = grid_in;
    numberOfLocalZones = nZones_in;

    nTotalCell = 0;
    nTotalFace = 0;
    nTotalNode = 0;
    nIFace = 0;
    nTotalNodeEachLocalZone = new int [numberOfLocalZones];
    nTotalFaceEachLocalZone = new int [numberOfLocalZones];
    nTotalCellEachLocalZone = new int [numberOfLocalZones];
    nIFaces = new int [numberOfLocalZones];

    face2NewFaceafterReorder = 0;
    face2Interface = 0;
    face2InterfaceafterReorder = 0;
    interface2BoundaryFaceafterReorder = 0;
    isGlobalInterFace = 0;
}

CombinGrid::~CombinGrid()
{
    delete [] nIFaces                           ;    nIFaces                            = nullptr;
    delete [] isFaceMerged                      ;    isFaceMerged                       = nullptr;
    delete [] face2Interface                    ;    face2Interface                     = nullptr;
    delete [] isGlobalInterFace                 ;    isGlobalInterFace                  = nullptr;
    delete [] left_cell_of_all_face             ;    left_cell_of_all_face              = nullptr;
    delete [] right_cell_of_all_face            ;    right_cell_of_all_face             = nullptr;
    delete [] node_number_of_all_face           ;    node_number_of_all_face            = nullptr;
    delete [] nTotalNodeEachLocalZone           ;    nTotalNodeEachLocalZone            = nullptr;
    delete [] nTotalFaceEachLocalZone           ;    nTotalFaceEachLocalZone            = nullptr;
    delete [] nTotalCellEachLocalZone           ;    nTotalCellEachLocalZone            = nullptr;
    delete [] face2NewFaceafterReorder          ;    face2NewFaceafterReorder           = nullptr;
    delete [] face2InterfaceafterReorder        ;    face2InterfaceafterReorder         = nullptr;
    delete [] interface2BoundaryFaceafterReorder;    interface2BoundaryFaceafterReorder = nullptr;

    DelPointer2(node2all);
    DelPointer2(face2all);
    DelPointer2(cell2all);
}

void CombinGrid::CheckifIndependentZoneExist()
{
    int myid   = PHMPI::GetCurrentProcessorID();
    int nZones = PHMPI::GetNumberofGlobalZones();

    int localIndex = 0;
    for (int iZone = 0; iZone < numberOfLocalZones; ++ iZone)
    {
        int gridIndex = grid_in[iZone]->GetZoneID();
        zoneIDinCurrentProc.insert(gridIndex);
        zoneIndexMap.insert(pair<int, int>(gridIndex, localIndex));

        localIndex ++;
    }

    if (numberOfLocalZones != 1)
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int processorID = PHMPI::GetZoneProcessorID(iZone);

            if (processorID == myid)
            {
                bool isNeighborinThisProc = false;

                ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
                std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
                
                for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
                {
                    int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

                    set<int>::iterator iter;
                    for (iter = zoneIDinCurrentProc.begin(); iter != zoneIDinCurrentProc.end(); ++ iter)
                    {
                        int zoneID = *iter;
                        if (neighborZone == zoneID)
                        {
                            isNeighborinThisProc = true;
                            break;
                        }
                    }
                }

                if (!isNeighborinThisProc)
                {
                    TK_Exit::ExceptionExit("Error: An single zone that couldn't merged to any others has been found!");
                }
            }
        }
    }

    WriteLogFile("All zones in this processor could be merged!");
    PrintToWindow("All zones in this processor could be merged!\n");
}

Grid * CombinGrid::CreateCombinationGrid()
{
    CheckifIndependentZoneExist();

    Init();

    //! Warning:
    //! All assumptions are that multi blocks in a process are merged into one block,
    //! If there is no common face between two blocks, a lot logic relationships need to be reconsidered.
    //! For example,MergeCell2All is the zone number of new grid.
    MergeCell2All();
    MergeNode2All();
    MergeFace2All();
    MergeInterface2All();

    int nTotalNode = GetNTotalNode();
    int nTotalFace = GetNTotalFace();
    int nTotalCell = GetNTotalCell();
    int nIFace = GetNIFace();

    ostringstream oss;
    oss << "After nodes faces cells merged ..." << "\n";
    oss << "  Number of Points     : " << nTotalNode << "\n";
    oss << "  Number of Faces      : " << nTotalFace << "\n";
    oss << "  Number of Cells      : " << nTotalCell << "\n";
    oss << "  Number of InterFaces : " << nIFace << "\n";
    PrintToWindow(oss);
    WriteLogFile(oss);

    int myid = PHMPI::GetCurrentProcessorID();
    Grid *grid_all = CreateGridGeneral(UNSTRUCTGRID, new GridID(myid), 0, GetDim());

    grid_all->SetNTotalNode(nTotalNode);
    grid_all->SetNTotalFace(nTotalFace);
    grid_all->SetNTotalCell(nTotalCell);

    return grid_all;
}

void CombinGrid::Init()
{
    for (int iZone = 0; iZone < numberOfLocalZones; ++ iZone)
    {
        nIFaces[iZone] = grid_in[iZone]->GetNIFace();
        nTotalNodeEachLocalZone[iZone] = grid_in[iZone]->GetNTotalNode();
        nTotalFaceEachLocalZone[iZone] = grid_in[iZone]->GetNTotalFace();
        nTotalCellEachLocalZone[iZone] = grid_in[iZone]->GetNTotalCell();

        nIFace     += nIFaces[iZone];
        nTotalCell += nTotalCellEachLocalZone[iZone];
        nTotalFace += nTotalFaceEachLocalZone[iZone];
        nTotalNode += nTotalNodeEachLocalZone[iZone];
    }

    InitLocal2All();
}

void CombinGrid::InitLocal2All()
{
    WriteLogFile("Init local to all ...");
    InitNode2All();
    InitFace2All();
    InitCell2All();
}

void CombinGrid::RunCombinGrid()
{
    grid_all = UnstructGridCast(CreateCombinationGrid());

    ResetCoordinate();
    ResetNodeNumberOfEachFace();
    ResetFace2Node();
    ResetLeftRightCell();
    MergeCell2Node();
    MergeBCType();
    SwapInterfaceData();

    PrintToWindow("Combine grids finished ...\n");
}

void CombinGrid::MergeInterface2All()
{
    WriteLogFile("Merge interface2all ...");
    PrintToWindow("Merge interface2all ...\n");

    int nTotalFace = GetNTotalFace();
    isGlobalInterFace = new bool [nTotalFace];
    PHSPACE::SetField(isGlobalInterFace, false, nTotalFace);

    face2Interface = new int [nTotalFace];
    PHSPACE::SetField(face2Interface, -1, nTotalFace);

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");
    int nMixingPlaneInterfaces = 0;

    for (int iZone = 0; iZone < numberOfLocalZones; ++ iZone)
    {
        UnstructGrid *gridLocal = UnstructGridCast(grid_in[iZone]);

        UnstructBCSet *unstructBCSet = gridLocal->GetUnstructBCSet();
        int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

        int *face2all = GetFace2All(iZone);

        int nMergedInterfaces = 0;
        for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
            int bcType = bcRegion->GetBCType();
            string bcName = bcRegion->GetBCName();

            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                int globalFaceID = face2all[iFace];

                if (bcType == PHENGLEI::INTERFACE)
                {
                    if (!isFaceMerged[globalFaceID])
                    {
                        isGlobalInterFace[globalFaceID] = true;
                    }else
                    {
                        ++ nMergedInterfaces;
                    }
                }

                if (referenceFrame == ROTATIONAL_FRAME)
                {
                    int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");

                    string MixingPlane[100];
                    GlobalDataBase::GetData("MixingPlane", &MixingPlane, PHSTRING, 2 * nTurboZone - 2);

                    string MixingInName, MixingOutName;
                    if (iZone == 0)
                    {
                        MixingOutName = MixingPlane[iZone];
                    }
                    else if (iZone == nTurboZone - 1)
                    {
                        MixingInName = MixingPlane[2 * iZone - 1];
                    }
                    else
                    {
                        MixingInName  = MixingPlane[2 * iZone - 1];
                        MixingOutName = MixingPlane[2 * iZone];
                    }
                   
                    if (bcName == MixingInName || bcName == MixingOutName)
                    {
                        ++nMergedInterfaces;
                        ++nMixingPlaneInterfaces;
                    }
                }
            }
        }

        //! this means rotating frame with mixing plane.
        if (referenceFrame == ROTATIONAL_FRAME && nMixingPlaneInterfaces > 0)
        {
            WriteLogFile("number of Mixing Plane: ", nMixingPlaneInterfaces);
        }
        else
        {
            WriteLogFile("number of Merged Interfaces: ", nMergedInterfaces);
        }
        if (numberOfLocalZones > 1 && nMergedInterfaces == 0)
        {
            TK_Exit::ExceptionExit("The sub-zone is independent of other zones, consider more!");
        }
    }

    int nIFaceGlobal = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        if(isGlobalInterFace[iFace])
        {
            face2Interface[iFace] = nIFaceGlobal;
            ++ nIFaceGlobal;
        }
    }

    this->nIFace = nIFaceGlobal;
}

void CombinGrid::MergeFace2All()
{
    WriteLogFile("Merge face2all ...");
    PrintToWindow("Merge face2all ...\n");

    left_cell_of_all_face = new int[nTotalFace];
    right_cell_of_all_face = new int[nTotalFace];
    node_number_of_all_face = new int[nTotalFace];

    isFaceMerged = new bool [nTotalFace];
    PHSPACE::SetField(isFaceMerged, false, nTotalFace);

    nTotalFace = 0;

    int globalFace_count = 0;
    for (int iZone = 0; iZone < numberOfLocalZones; ++ iZone)
    {
        UnstructGrid *unstr_current = UnstructGridCast(grid_in[iZone]);

        int nBFace_current                    = unstr_current->GetNBoundFace();
        int *node_number_of_each_face_current = unstr_current->GetNodeNumberOfEachFace();
        int *face2node_current                = unstr_current->GetFace2Node();
        int *left_cell_of_face_current        = unstr_current->GetLeftCellOfFace();
        int *right_cell_of_face_current       = unstr_current->GetRightCellOfFace();
        
        UnstructBCSet *unstructBCSet = unstr_current->GetUnstructBCSet();
        int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();
       
        InterfaceInfo *iinfo_current          = unstr_current->GetInterfaceInfo();
        int *interFace2ZoneID_current         = iinfo_current->GetInterFace2ZoneID();
        int *interFace2InterFaceID_current    = iinfo_current->GetInterFace2InterFaceID();

        int *node2all = GetNode2All(iZone);
        int *face2all = GetFace2All(iZone);
        int *cell2all = GetCell2All(iZone);

        int interface_current = 0;
        int node_sum = 0;
        for (int iFace = 0; iFace < nTotalFaceEachLocalZone[iZone]; ++ iFace)
        {
            int node_number_of_iface = node_number_of_each_face_current[iFace];
            if (iFace < nBFace_current)
            {
                UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
                int bctype = bcRegion->GetBCType();
                if (bctype == PHENGLEI::INTERFACE)
                {
                    int GridID_Neighbor = interFace2ZoneID_current[interface_current];

                    bool isNeighborinThisProc = false;
                    set<int>::iterator iter;
                    for (iter = zoneIDinCurrentProc.begin(); iter != zoneIDinCurrentProc.end(); ++ iter)
                    {
                        int zoneID = *iter;
                        if (GridID_Neighbor == zoneID)
                        {
                            isNeighborinThisProc = true;
                            break;
                        }
                    }

                    int localIndexNeighbor = zoneIndexMap[GridID_Neighbor];
                    if ((isNeighborinThisProc) && (localIndexNeighbor < iZone))
                    {
                        int interface_Neighbor = interFace2InterFaceID_current[interface_current];
                        UnstructGrid *unstr_Neighbor = UnstructGridCast(grid_in[localIndexNeighbor]);
                        InterfaceInfo *iinfo_Neighbor = unstr_Neighbor->GetInterfaceInfo();
                        int *interFace2BoundaryFace_Neighbor = iinfo_Neighbor->GetInterFace2BoundaryFace();

                        int bcFace_Neighbor = interFace2BoundaryFace_Neighbor[interface_Neighbor];
                        int *face2all_Neighbor = GetFace2All(localIndexNeighbor);
                        face2all[iFace] = face2all_Neighbor[bcFace_Neighbor];

                        right_cell_of_all_face[face2all_Neighbor[bcFace_Neighbor]] = cell2all[left_cell_of_face_current[iFace]];

                        isFaceMerged[face2all_Neighbor[bcFace_Neighbor]] = true;
                        node_sum += node_number_of_iface;

                        ++ interface_current;
                        continue;
                    }

                    ++ interface_current;
                }
            }

            face2all[iFace] = globalFace_count;

            node_number_of_all_face[globalFace_count] = node_number_of_iface;

            for (int iNode = 0; iNode < node_number_of_iface; ++ iNode)
            {
                int nodeindex = face2node_current[node_sum + iNode];
                allface2node.push_back(node2all[nodeindex]);
            }
            node_sum += node_number_of_iface;
            
            if (left_cell_of_face_current[iFace] >= nTotalCellEachLocalZone[iZone] || left_cell_of_face_current[iFace] < 0)
            {
                left_cell_of_all_face[globalFace_count] = -1;
            }
            else
            {
                left_cell_of_all_face[globalFace_count] = cell2all[left_cell_of_face_current[iFace]];
            }
            

            if (right_cell_of_face_current[iFace] >= nTotalCellEachLocalZone[iZone] || right_cell_of_face_current[iFace] < 0)
            {
                right_cell_of_all_face[globalFace_count] = -1;
            }
            else
            {
                right_cell_of_all_face[globalFace_count] = cell2all[right_cell_of_face_current[iFace]];
            }
            
            ++ globalFace_count;
        }

        nTotalFace = globalFace_count;
    }
}

void CombinGrid::MergeCell2All()
{
    WriteLogFile("Merge cell2all ...");
    PrintToWindow("Merge cell2all ...\n");

    int sum = 0;
    nTotalCell = 0;

    for (int iZone = 0; iZone < numberOfLocalZones; ++ iZone)
    {
        int *cell2all = GetCell2All(iZone);

        for (int i = 0; i < nTotalCellEachLocalZone[iZone]; ++ i)
        {
            cell2all[i] = sum;
            sum ++;
        }

        nTotalCell = sum;
    }

    int minNTCell, maxNTCell;
    PH_Reduce(&nTotalCell, &minNTCell, 1, MPI_MIN);
    PH_Reduce(&nTotalCell, &maxNTCell, 1, MPI_MAX);

    if (PHMPI::GetCurrentProcessorID() == PHMPI::GetServerProcessorID())
    {
        ostringstream oss;
        oss << "Min && Max number of cell after merged : " << minNTCell << ", " << maxNTCell << endl;
        oss << "Min to Max nTCell ratio = " << (minNTCell * 1.0) / (maxNTCell * 1.0) << endl;
        PrintToWindow(oss);
        WriteLogFile(oss);
    }

    PH_Barrier();
}

//void CombinGrid::MergeNode2AllNew()
//{
//    WriteLogFile("Merge node2all ...");
//    PrintToWindow("Merge node2all ...\n");
//
//    Globalx.reserve(nTotalNode);
//    Globaly.reserve(nTotalNode);
//    Globalz.reserve(nTotalNode);
//
//    nTotalNode = 0;
//    for (int iZone = 0; iZone < numberOfLocalZones; ++ iZone)
//    {
//        int *node2allCurrent = GetNode2All(iZone);
//
//        UnstructGrid   *unstrCurrent                 = UnstructGridCast(grid_in[iZone]);
//        RDouble        *xCoorCurrent                 = unstrCurrent->GetX();
//        RDouble        *yCoorCurrent                 = unstrCurrent->GetY();
//        RDouble        *zCoorCurrent                 = unstrCurrent->GetZ();
//        int             nBoundFace                   = unstrCurrent->GetNBoundFace();
//        int            *face2nodeCurrent             = unstrCurrent->GetFace2Node();
//        int            *nodeNumberOfEachFaceCurrent  = unstrCurrent->GetNodeNumberOfEachFace();
//        UnstructBCSet **bcRecord                     = unstrCurrent->GetBCRecord();
//        InterfaceInfo  *iinfoCurrent                 = unstrCurrent->GetInterfaceInfo();
//        int            *interFace2ZoneIDCurrent      = iinfoCurrent->GetInterFace2ZoneID();
//        int            *interFace2InterFaceIDCurrent = iinfoCurrent->GetInterFace2InterFaceID();
//        RDouble         distanceTolerance            = unstrCurrent->CalMinEdgeLength();
//
//        int *isMarked = new int[nTotalNodeEachLocalZone[iZone]];
//        PHSPACE::SetField(isMarked, 0, nTotalNodeEachLocalZone[iZone]);
//
//        //! Mark interface nodes.
//        int nodePos          = 0;
//        int interfaceCurrent = 0;
//        RDouble disX, disY, disZ;
//        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
//        {
//            int bcType = bcRecord[iFace]->GetKey();
//            if (bcType == PHENGLEI::INTERFACE)
//            {
//                int gridIDNeighbor = interFace2ZoneIDCurrent[interfaceCurrent];
//                int localIndexNeighbor = zoneIDinCurrentProc[gridIDNeighbor];
//
//                bool isNeighborInThisProc = false;
//                if (localIndexNeighbor != -1)
//                {
//                    isNeighborInThisProc = true;
//                }
//
//                if ((isNeighborInThisProc) && (localIndexNeighbor < iZone))
//                {
//                    int           *node2allNeighbor               = GetNode2All(localIndexNeighbor);
//                    int            interfaceNeighbor              = interFace2InterFaceIDCurrent[interfaceCurrent];
//                    UnstructGrid  *unstrNeighbor                  = UnstructGridCast(grid_in[localIndexNeighbor]);
//                    RDouble       *xCoorNeighbor                  = unstrNeighbor->GetX();
//                    RDouble       *yCoorNeighbor                  = unstrNeighbor->GetY();
//                    RDouble       *zCoorNeighbor                  = unstrNeighbor->GetZ();
//                    InterfaceInfo *iinfoNeighbor                  = unstrNeighbor->GetInterfaceInfo();
//                    int           *interFace2BoundaryFaceNeighbor = iinfoNeighbor->GetInterFace2BoundaryFace();
//                    int            bcFaceNeighbor                 = interFace2BoundaryFaceNeighbor[interfaceNeighbor];
//                    int           *face2nodeNeighbor              = unstrNeighbor->GetFace2Node();
//                    long long int *face2nodeSubscriptNeighbor     = unstrNeighbor->GetFace2NodeSubscript();
//                    long long int  nodePosNeighbor                = face2nodeSubscriptNeighbor[bcFaceNeighbor];
//
//                    for (int iNodeCurrent = 0; iNodeCurrent < nodeNumberOfEachFaceCurrent[iFace]; ++ iNodeCurrent)
//                    {
//                        int nodeCurrentIndex = face2nodeCurrent[nodePos + iNodeCurrent];
//
//                        for (int iNodeNeighbor = 0; iNodeNeighbor < nodeNumberOfEachFaceCurrent[iFace]; ++ iNodeNeighbor)
//                        {
//                            int nodeNeighborIndex = face2nodeNeighbor[nodePosNeighbor + iNodeNeighbor];
//
//                            disX = fabs(xCoorCurrent[nodeCurrentIndex] - xCoorNeighbor[nodeNeighborIndex]);
//                            disY = fabs(yCoorCurrent[nodeCurrentIndex] - yCoorNeighbor[nodeNeighborIndex]);
//                            disZ = fabs(zCoorCurrent[nodeCurrentIndex] - zCoorNeighbor[nodeNeighborIndex]);
//
//                            if (disX < distanceTolerance && disY < distanceTolerance && disZ < distanceTolerance)
//                            {
//                                node2allCurrent[nodeCurrentIndex] = node2allNeighbor[nodeNeighborIndex];
//                                isMarked[nodeCurrentIndex] = 1;
//                                break;
//                            }
//                        }
//                    }
//                }
//
//                interfaceCurrent ++;
//            }
//
//            nodePos += nodeNumberOfEachFaceCurrent[iFace];
//        }
//
//        //! Mark remaining nodes.
//        for (int iNode = 0; iNode < nTotalNodeEachLocalZone[iZone]; ++ iNode)
//        {
//            if (isMarked[iNode])
//            {
//                continue;
//            }
//
//            Globalx.push_back(xCoorCurrent[iNode]);
//            Globaly.push_back(yCoorCurrent[iNode]);
//            Globalz.push_back(zCoorCurrent[iNode]);
//            node2allCurrent[iNode] = nTotalNode;
//            nTotalNode ++;
//        }
//
//        delete [] isMarked;    isMarked = nullptr;
//    }
//}

void CombinGrid::MergeNode2All()
{
    WriteLogFile("Merge node2all ...");
    PrintToWindow("Merge node2all ...\n");

    int sum = 0;
    Globalx.reserve(nTotalNode);
    Globaly.reserve(nTotalNode);
    Globalz.reserve(nTotalNode);

    nTotalNode = 0;
    RDouble minDistance, Disx, Disy, Disz;
    minDistance = LARGE;
    bool discard;
    vector <int> bc_node_all;

    for (int iZone = 0; iZone < numberOfLocalZones; ++ iZone)
    {
        int *node2all = GetNode2All(iZone);

        UnstructGrid *unstr = UnstructGridCast(grid_in[iZone]);
        RDouble *x2 = grid_in[iZone]->GetX();
        RDouble *y2 = grid_in[iZone]->GetY();
        RDouble *z2 = grid_in[iZone]->GetZ();
        int *node_number_of_each_face = unstr->GetNodeNumberOfEachFace();
        int *face2node = unstr->GetFace2Node();
        minDistance = MIN(minDistance, grid_in[iZone]->CalMinEdgeLength());
        RDouble distanceTolerance = minDistance / 10;

        int nBoundFace = grid_in[iZone]->GetNBoundFace();

        int *isBCNode_iZone = new int [nTotalNodeEachLocalZone[iZone]];
        PHSPACE::SetField(isBCNode_iZone, 0, nTotalNodeEachLocalZone[iZone]);

        int nsum = 0;

        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            for (int iNode = 0; iNode < node_number_of_each_face[iFace]; ++ iNode)
            {
                isBCNode_iZone[face2node[nsum+iNode]] = 1;
            }
            nsum += node_number_of_each_face[iFace];
        }
        
        uint_t bc_node_number = bc_node_all.size();
        for (int iNode = 0; iNode < nTotalNodeEachLocalZone[iZone]; ++ iNode)
        {
            discard = false;

            if (IsBcNode(isBCNode_iZone, iNode))
            {

                for (int jbcNode = 0; jbcNode < bc_node_number; jbcNode ++)
                {
                    Disx = fabs(x2[iNode] - Globalx[bc_node_all[jbcNode]]);
                    Disy = fabs(y2[iNode] - Globaly[bc_node_all[jbcNode]]);
                    Disz = fabs(z2[iNode] - Globalz[bc_node_all[jbcNode]]);

                    if (Disx < distanceTolerance && Disy < distanceTolerance && Disz < distanceTolerance)
                    {
                        node2all[iNode] = bc_node_all[jbcNode];
                        discard = true;
                        break;
                    }
                    else
                    {
                        continue;
                    }
                }
                if (discard != true)
                {
                    bc_node_all.push_back(sum);
                }
            }
            
            if (discard) continue;

            Globalx.push_back(x2[iNode]);
            Globaly.push_back(y2[iNode]);
            Globalz.push_back(z2[iNode]);
            node2all[iNode] = sum;
            sum ++ ;
        }
        nTotalNode = sum;
        delete [] isBCNode_iZone;
    }
}

void CombinGrid::MergeBCType()
{
    WriteLogFile("Merge boundary condition ...");
    PrintToWindow("Merge boundary condition ...\n");

    int nTFace_all = grid_all->GetNTotalFace();

    int nBFace_all = 0;
    for (int iZone = 0; iZone < numberOfLocalZones; ++ iZone)
    {
        UnstructGrid *unstr = UnstructGridCast(grid_in[iZone]);
        int *face2all = GetFace2All(iZone);

        UnstructBCSet *unstructBCSet = unstr->GetUnstructBCSet();
        int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

        for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);

            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                // iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                int globalFaceID = face2all[iFace];

                if (bcRegion->GetBCType() == PHENGLEI::INTERFACE)
                {
                    if (isFaceMerged[globalFaceID])
                    {
                        //! The interface has been merged, it no longer be boundary face now.
                        continue;
                    }
                }
                ++ nBFace_all;
            }
        }
    }

    grid_all->SetNBoundFace(nBFace_all);
    UnstructBCSet **bcr_all = new UnstructBCSet * [nTFace_all];

    for (int iZone = 0; iZone < numberOfLocalZones; ++ iZone)
    {
        UnstructGrid *unstr = UnstructGridCast(grid_in[iZone]);
        int nBoundFace = unstr->GetNBoundFace();
        int *face2all = GetFace2All(iZone);

        UnstructBCSet **bcr = unstr->GetBCRecord();
        
        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            int globalFaceID = face2all[iFace];
            if (bcr[iFace]->GetKey() == PHENGLEI::INTERFACE)
            {
                if (isFaceMerged[globalFaceID])
                {
                    //! The interface has been merged, it no longer be boundary face now.
                    continue;
                }
            }

            bcr_all[globalFaceID] = new UnstructBCSet();
            bcr_all[globalFaceID]->SetKey(bcr[iFace]->GetKey());
            bcr_all[globalFaceID]->SetBCName(bcr[iFace]->GetBCName());
        }
    }

    vector< vector<int> > vface2node;
    vface2node.reserve(nTFace_all);
    int *node_number_of_each_face = grid_all->GetNodeNumberOfEachFace();
    int *face2node = grid_all->GetFace2Node();
    int *left_cell_of_face  = grid_all->GetLeftCellOfFace();
    int *right_cell_of_face = grid_all->GetRightCellOfFace();

    int nsum = 0;
    for (int iFace = 0; iFace < nTFace_all; ++ iFace)
    {
        vector< int > eachface2node;
        eachface2node.reserve(node_number_of_each_face[iFace]);
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            eachface2node.push_back(face2node[nsum]);
            ++ nsum;
        }
        vface2node.push_back(eachface2node);
    }

    //! Swap faces make sure boundary faces are in the first section (iFace < nBoundFace),
    //! and the interior faces are in the second section (iFace >= nBoundFace).
    face2NewFaceafterReorder = new int [nTFace_all];
    for (int iFace = 0; iFace < nTFace_all; ++ iFace)
    {
        face2NewFaceafterReorder[iFace] = iFace;
    }

    int jstart = nBFace_all;
    for (int iFace = 0; iFace < nBFace_all; ++ iFace)
    {
        if (right_cell_of_face[iFace] != -1)
        {
            //! iFace is interior face, it must be reject out of nBoundFace.
            for (int jFace = jstart; jFace < nTFace_all; ++ jFace)
            {
                if (right_cell_of_face[jFace] == -1)
                {
                    SWAP(left_cell_of_face[iFace], left_cell_of_face[jFace]);
                    SWAP(right_cell_of_face[iFace], right_cell_of_face[jFace]);
                    SWAP(vface2node[iFace], vface2node[jFace]);
                    SWAP(node_number_of_each_face[iFace], node_number_of_each_face[jFace]);
                    SWAP(bcr_all[iFace], bcr_all[jFace]);

                    face2NewFaceafterReorder[iFace] = jFace;
                    face2NewFaceafterReorder[jFace] = iFace;

                    jstart = jFace + 1;
                    break;
                }
            }
        }
    }

    nsum = 0;
    for (int iFace = 0; iFace < nTFace_all; ++ iFace)
    {
        for (int iNode = 0; iNode < node_number_of_each_face[iFace]; ++ iNode)
        {
            face2node[nsum] = vface2node[iFace][iNode];
            ++ nsum;
        }
    }

    UnstructBCSet **bcr_all2 = new UnstructBCSet * [nBFace_all];
    grid_all->SetBCRecord(bcr_all2);

    WriteLogFile("nTFace_all = ", nTFace_all);
    WriteLogFile("nBFace_all = ", nBFace_all);

    face2InterfaceafterReorder = new int [nTFace_all];
    PHSPACE::SetField(face2InterfaceafterReorder, -1, nTFace_all);

    interface2BoundaryFaceafterReorder = new int [nBFace_all]; // nIFace of grid_all has not been set yet.
    PHSPACE::SetField(interface2BoundaryFaceafterReorder, -1, nBFace_all);

    int interfaceID = 0;
    
    for (int iFace = 0; iFace < nBFace_all; ++ iFace)
    {
        bcr_all2[iFace] = bcr_all[iFace];

        if (bcr_all[iFace]->GetKey() == PHENGLEI::INTERFACE)
        {
            interface2BoundaryFaceafterReorder[interfaceID] = iFace;
            face2InterfaceafterReorder[iFace] = interfaceID ++;
        }
    }

    delete [] bcr_all;

    nIFace = interfaceID;
    WriteLogFile("nIFace_all = ", nIFace);
    InterfaceInfo *interfaceInfo = new InterfaceInfo(nIFace);
    grid_all->SetInterfaceInfo(interfaceInfo);
}

bool CombinGrid::IsBcNode(int *p, int node_in)
{
    bool bc_node = false;

    if (p[node_in] == 1) bc_node = true;
    return bc_node;
}

void CombinGrid::InitNode2All()
{
    node2all = new int *[numberOfLocalZones];
    for (int iZone = 0; iZone < numberOfLocalZones; ++ iZone)
    {
        int nTotalNode = grid_in[iZone]->GetNTotalNode();
        node2all[iZone] = new int [nTotalNode];
    }
}

void CombinGrid::InitFace2All()
{
    face2all = new int *[numberOfLocalZones];
    for (int iZone = 0; iZone < numberOfLocalZones; ++ iZone)
    {
        int nTotalFace = grid_in[iZone]->GetNTotalFace();
        face2all[iZone] = new int [nTotalFace];
    }
}

void CombinGrid::InitIsBCFaceMerged()
{

}

void CombinGrid::InitCell2All()
{
    cell2all = new int *[numberOfLocalZones];
    for (int iZone = 0; iZone < numberOfLocalZones; ++ iZone)
    {
        int nTotalCell = grid_in[iZone]->GetNTotalCell();
        cell2all[iZone] = new int [nTotalCell];
    }
}

void CombinGrid::ResetCoordinate()
{
    WriteLogFile("Combine coordinates ...\n");
    PrintToWindow("Combine coordinates ...\n");

    int nTNode_all = grid_all->GetNTotalNode();

    RDouble *x_all = new RDouble[nTNode_all];
    RDouble *y_all = new RDouble[nTNode_all];
    RDouble *z_all = new RDouble[nTNode_all];

    vector<RDouble> Globalx = GetGlobalx();
    vector<RDouble> Globaly = GetGlobaly();
    vector<RDouble> Globalz = GetGlobalz();

    grid_all->SetX(x_all);
    grid_all->SetY(y_all);
    grid_all->SetZ(z_all);

    for (int iNode = 0; iNode < nTNode_all; ++ iNode)
    {
        x_all[iNode] = Globalx[iNode];
        y_all[iNode] = Globaly[iNode];
        z_all[iNode] = Globalz[iNode];
    }

    grid_all->ComputeMinMaxBox();
}

void CombinGrid::ResetNodeNumberOfEachFace()
{
    WriteLogFile("Combine node_number_of_each_face ...\n");
    PrintToWindow("Combine node_number_of_each_face ...\n");

    int nTFace_all = grid_all->GetNTotalFace();

    int *node_number_of_each_face = new int[nTFace_all];

    grid_all->SetNodeNumberOfEachFace(node_number_of_each_face);

    int *node_number_of_all_face = GetNodeNumberOfAllFace();

    for (int iNode = 0; iNode < nTFace_all; ++ iNode)
    {
        node_number_of_each_face[iNode] = node_number_of_all_face[iNode];
    }
}

void CombinGrid::ResetFace2Node()
{
    PrintToWindow("Combine face2node ...\n");
    WriteLogFile("Combine face2node ...\n");

    vector<int> allface2node = GetAllFace2Node();
    uint_t nsum = allface2node.size();

    int *face2node = new int[nsum];
    grid_all->SetFace2Node(face2node);
    
    for (int i = 0; i < nsum; ++ i)
    {
        face2node[i] = allface2node[i];
    }
}

void CombinGrid::ResetLeftRightCell()
{
    WriteLogFile("Reset Left Right Cell ...");
    PrintToWindow("Reset Left Right Cell ...\n");

    int nTFace_all = grid_all->GetNTotalFace();
    int *left_cell_of_all_face  = GetLeftCellOfAllFace();
    int *right_cell_of_all_face = GetRightCellOfAllFace();

    int *left_cell_of_face  = new int[nTFace_all];
    int *right_cell_of_face = new int[nTFace_all];
    grid_all->SetLeftCellOfFace (left_cell_of_face);
    grid_all->SetRightCellOfFace(right_cell_of_face);

    for (int iFace = 0; iFace < nTFace_all; ++ iFace)
    {
        left_cell_of_face[iFace]  = left_cell_of_all_face[iFace];
        right_cell_of_face[iFace] = right_cell_of_all_face[iFace];
    }
}

void CombinGrid::MergeCell2Node()
{
    if (UnstructGridCast(grid_in[0])->IsCell2NodeExist() == false)
    {
        //! Judge if need merge cell2node.
        return;
    }
    WriteLogFile("Merge cell to node ...");
    PrintToWindow("Merge cell to node ...\n");

    vector< int > node_number_of_each_cell_all;
    node_number_of_each_cell_all.resize(0);
    node_number_of_each_cell_all.reserve(nTotalCell);

    vector< vector<cgsize_t> > cell2node_all;
    cell2node_all.resize(0);
    cell2node_all.reserve(nTotalCell);

    int cell2node_size = 0;

    for (int iZone = 0; iZone < numberOfLocalZones; ++ iZone)
    {
        int *gnode_number_of_each_cell = UnstructGridCast(grid_in[iZone])->GetNodeNumberOfEachCell();
        int *gcell2node = UnstructGridCast(grid_in[iZone])->GetCell2Node();

        int *gnode2all = GetNode2All(iZone);

        int count = 0;
        for (int iCell = 0; iCell < nTotalCellEachLocalZone[iZone]; ++ iCell)
        {
            node_number_of_each_cell_all.push_back(gnode_number_of_each_cell[iCell]);
            cell2node_size += gnode_number_of_each_cell[iCell];

            vector<cgsize_t> single_cell;
            single_cell.resize(0);
            single_cell.reserve(gnode_number_of_each_cell[iCell]);

            for (int iNode = 0; iNode < gnode_number_of_each_cell[iCell]; ++ iNode)
            {
                int node_index = gcell2node[count];
                ++ count;
                single_cell.push_back(gnode2all[node_index]);
            }
            cell2node_all.push_back(single_cell);
        }
    }

    int *node_number_of_each_cell = new int [nTotalCell];
    int *cell2node = new int [cell2node_size];

    int count = 0;
    for (cgsize_t iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        node_number_of_each_cell[iCell] = node_number_of_each_cell_all[iCell];
        for (int iNode =0; iNode < node_number_of_each_cell[iCell]; ++ iNode)
        {
            cell2node[count] = static_cast<int>(cell2node_all[iCell][iNode]);
            ++ count;
        }
    }
    grid_all->SetNodeNumberOfEachCell(node_number_of_each_cell);
    grid_all->SetCell2Node(cell2node);
}

Grid ** CombinGrid::GetGridAll() const
{
    Grid **grid_array = new Grid * [1];

    grid_array[0] = grid_all;

    return grid_array;
}

void CombinGrid::DumpCombinationGrid()
{
    string out_gfile;
    GlobalDataBase::GetData("out_gfile", &out_gfile, PHSTRING, 1);

    int isOverset = 0;
    GlobalDataBase::GetData("codeOfOversetGrid", &isOverset, PHINT, 1);

    using namespace PHMPI;

    int nZones = 1;

    Grid **grid_array = new Grid * [nZones];
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        grid_array[iZone] = grid_all;
    }

    ConstructGlobalInterfaceLink(grid_array, nZones);
    DumpGrid(out_gfile, grid_array, nZones);

    delete [] grid_array;
}

void CombinGrid::CompressInterfaceData(int zoneID, DataContainer *sendBuffer, int neighborZone)
{
    int localZoneID = region->GetProcessGlobalZoneIndexToLocalZoneIndex(zoneID);
    Grid *grid = grid_in[localZoneID];
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    int myid = PHMPI::GetCurrentProcessorID();
    int iNeighborZone                            = interfaceInfo->FindIthNeighbor(neighborZone);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = interfaceInfo->GetFaceIndexForSend(iNeighborZone);
    int *localInterfaceToBoundaryFaceIndexesMapping = interfaceInfo->GetInterFace2BoundaryFace();
    sendBuffer->MoveToBegin();
    
//  WriteLogFile("currentZone = ", zoneID);
//  WriteLogFile("neighborZone = ", neighborZone);

    int *face2all = GetFace2All(localZoneID);
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int interFaceID = interfaceIndexContainerForSend[iLocalFace];
        int boundaryFaceID  = localInterfaceToBoundaryFaceIndexesMapping[interFaceID];
        int oldGlobalFaceID = face2all[boundaryFaceID];
        int newGlobalFaceID = face2NewFaceafterReorder[oldGlobalFaceID];

        //! interFace2ZoneID.
        int newZoneID = myid;
        PHWrite(sendBuffer, newZoneID);

        //! interFace2InterFaceID.
        int newInterfaceID = face2InterfaceafterReorder[newGlobalFaceID];
        PHWrite(sendBuffer, newInterfaceID);
    }
}

void CombinGrid::DecompressInterfaceData(int zoneID, DataContainer *recieveBuffer, int neighborZone)
{
    int localZoneID = region->GetProcessGlobalZoneIndexToLocalZoneIndex(zoneID);
    Grid *grid = grid_in[localZoneID];
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    //! Interface information of local sub-zones.
    int iNeighborZone                            = interfaceInfo->FindIthNeighbor(neighborZone);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = interfaceInfo->GetFaceIndexForRecv(iNeighborZone);
    int *localInterfaceToBoundaryFaceIndexesMapping = interfaceInfo->GetInterFace2BoundaryFace();

    //! Interface information of global merged zone.
    InterfaceInfo *globalInterfaceInfor = grid_all->GetInterfaceInfo();
    int *interFace2ZoneID       = globalInterfaceInfor->GetInterFace2ZoneID();
    int *interFace2InterFaceID  = globalInterfaceInfor->GetInterFace2InterFaceID();
    int *interFace2BoundaryFace = globalInterfaceInfor->GetInterFace2BoundaryFace();

    recieveBuffer->MoveToBegin();
    int *face2all = GetFace2All(localZoneID);
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int interFaceID = interfaceIndexContainerForReceive[iLocalFace];
        int boundaryFaceID = localInterfaceToBoundaryFaceIndexesMapping[interFaceID];
        int oldGlobalFaceID = face2all[boundaryFaceID];
        int newGlobalFaceID = face2NewFaceafterReorder[oldGlobalFaceID];
        int newInterfaceID  = face2InterfaceafterReorder[newGlobalFaceID];

        //! interFace2ZoneID.
        int newNeighborZoneID;
        PHRead(recieveBuffer, newNeighborZoneID);

        //! interFace2InterFaceID
        int newNeighborInterfaceID;
        PHRead(recieveBuffer, newNeighborInterfaceID);

        if (newInterfaceID >= 0)
        {
            //! This interface has not been merged, so it also be interface in merged grid.
            interFace2ZoneID[newInterfaceID]       = newNeighborZoneID;
            interFace2InterFaceID[newInterfaceID]  = newNeighborInterfaceID;
            interFace2BoundaryFace[newInterfaceID] = newGlobalFaceID;
        }
    }
}

void CombinGrid::SwapInterfaceData()
{
    if (PHMPI::GetNumberOfProcessor() == 1) return;

    PrintToWindow("Swapping interface data ...\n");

    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int nZones = PHMPI::GetNumberofGlobalZones();
    vector< PH_Request > requestContainer;
    vector< vector< DataContainer * > > receivedDataBuffer;
    vector< vector< DataContainer * > > sendDataBuffer;
    receivedDataBuffer.resize(nZones);
    sendDataBuffer.resize(nZones);

    //! Step 0: Compressing data firstly.
    WriteLogFile("        Compressing data ...");
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcessor = PHMPI::GetZoneProcessorID(iZone);

        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
        if (currentProcessor == sendProcessor)
        {
            sendDataBuffer[iZone].resize(numberOfNeighbor);
        }

        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate and compress the buffers for sending.
            DataContainer *sendBuffer = 0;
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = new DataContainer;
                sendDataBuffer[iZone][iNeighbor] = sendBuffer;

                //! Compress the send information into the actkey.
                CompressInterfaceData(iZone, sendBuffer, neighborZone);
            }
        }
    }

    //! Step 1: Communication.
    WriteLogFile("        Communication ...");
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        //! Communicating.
        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessor    = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate the buffers for sending and receiving.
            DataContainer *receivedBuffer = 0;
            DataContainer *sendBuffer = 0;
            if (currentProcessor == receiveProcessor)
            {
                receivedBuffer = new DataContainer;
                receivedDataBuffer[ iZone ].push_back(receivedBuffer);
            }
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = sendDataBuffer[iZone][iNeighbor];
            }

            int tag = iZone;

            //if (sendProcessor == receiveProcessor) continue;

            //! Communication.
            if (currentProcessor == sendProcessor)
            {
                if (sendProcessor != receiveProcessor)
                {
                    streamsize nlen = sendBuffer->Size();

                    //! Send the data to neighbors.
                    send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
                }
            }
            if (currentProcessor == receiveProcessor)
            {
                //! Get the neighbor order: which order that 'iZone' on neighbors of neighborZone.
                int neighborOrer = -1;
                ZoneNeighbor *globalNeighborZonesTemp = zoneConnectivity->GetZoneNeighbor(neighborZone);
                int numberOfNeighborTemp = globalNeighborZonesTemp->GetNumberOfNeighbor();

                for (int iNeighbor1 = 0; iNeighbor1 < numberOfNeighborTemp; ++ iNeighbor1)
                {
                    int neighborZoneTemp = globalNeighborZonesTemp->GetZoneIndexOfNeighbor(iNeighbor1);
                    if (neighborZoneTemp == iZone)
                    {
                        neighborOrer = iNeighbor1;
                        break;
                    }
                }

                ASSERT(neighborOrer != -1);

                if (sendProcessor != receiveProcessor)
                {
                    //! The data length of the received data is same to the length that send to the neighbor.
                    CharVecSizeType nlen = sendDataBuffer[neighborZone][neighborOrer]->Size();

                    //! Receive data from neighbors.
                    receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
                }
                else
                {
                    uint_t lastID = receivedDataBuffer[iZone].size();
                    SWAP(receivedDataBuffer[ iZone ][lastID - 1], sendDataBuffer[neighborZone][neighborOrer]);
                }
            }
        }
    }

    //! Step2: MPI waiting for Non-blocking communications.
    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_Wait(static_cast<int>(requestContainer.size()), &(requestContainer[0]));
    }

    //! Step3: Translate the data container.
    WriteLogFile("        Translating data ...");
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor == receiveProcessor)
            {
                DataContainer *receiveData = receivedDataBuffer[iZone][count];

                //! Because get into here when currentProcessor==receiveProcessor,
                //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!

                //! Decompress the interface data from datacontainer.
                DecompressInterfaceData(neighborZone, receiveData, iZone);

                ++ count;
            }
        }
    }

    WriteLogFile("        Free buffer data ...");
    //! Step4: Free the buffers.
    for (std::size_t iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < receivedDataBuffer[iDim].size(); ++ jDim)
        {
            delete receivedDataBuffer[iDim][jDim];
        }
    }
    receivedDataBuffer.clear();
    for (std::size_t iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < sendDataBuffer[ iDim ].size(); ++ jDim)
        {
            delete sendDataBuffer[iDim][jDim];
        }
    }
    sendDataBuffer.clear();
    requestContainer.clear();
}

void CombinGrid::DumpMergedGrid(const string &gridFileName, Grid **grids, int nBlocks)
{
    PrintToWindow("    Writing grid file ....\n");

    if (!PHMPI::IsParallelRun())
    {
        TK_Exit::ExceptionExit("fileMode == PHSPACE::CS_MODE has not been considered!");
    }
    else
    {
        DumpMergedGridP2P(gridFileName, grids, nBlocks);
    }
}

void CombinGrid::DumpMergedGridP2P(const string &gridFileName, Grid **grids, int nBlocks)
{
    using namespace PHMPI;

    //! Notice: Directly writing to grid file will affect the number of grid blocks,so you can not dump it at will.
    SetNumberOfGlobalZones(nBlocks);

    int nZones = GetNumberofGlobalZones();

    CreateZoneProcessorID(nZones);
    CreateZoneGridID(nZones);
    CreateZoneGridType(nZones);

    int *block_proc     = GetZoneProcessorID();
    int *block_idx      = GetZoneGridID();
    int *block_type     = GetZoneGridType();
    int *block_proc_inp = GetZoneProcessorID_INP();

    int myid = GetCurrentProcessorID();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        block_proc[iZone] = myid;
        block_idx [iZone] = myid;
        block_type[iZone] = grids[iZone]->Type();
    }

    if (block_proc_inp)
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            SetZoneProcessorID(iZone, block_proc_inp[iZone]);
        }
    }

    int m_block_proc = PHMPI::GetZoneDistributionMethod();

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        fstream file;
        string processorGridFileName = AddSymbolToFileName(gridFileName, '_', PHMPI::GetCurrentProcessorID());
        PHSPACE::OpenFile(file, processorGridFileName, ios_base::out|ios_base::binary|ios_base::trunc);

        int nBlocksLocal = 1;
        file.write(reinterpret_cast<char *>(&nBlocksLocal), sizeof(int));
        int blockIndex = block_idx[iZone];
        int blockType  = block_type[iZone];
        int blockProc  = block_proc[iZone];
        if (m_block_proc == DETERMIN_BY_PARTITION && IsConvertGridToMixGrid())
        {
            int *block_proc_dump = GetZoneProcessorIDDump();
            int blockProcDump = block_proc_dump[iZone];
            file.write(reinterpret_cast<char *>(&blockProcDump), nBlocksLocal * sizeof(int));
        }
        else
        {
            file.write(reinterpret_cast<char *>(&blockProc), nBlocksLocal * sizeof(int)); 
        }
        file.write(reinterpret_cast<char *>(&blockIndex), nBlocksLocal * sizeof(int));
        file.write(reinterpret_cast<char *>(&blockType), nBlocksLocal * sizeof(int));

        grids[iZone]->WriteGrid(file);

        PHSPACE::CloseFile(file);
    }
}

}


