#include <list>
#include "cgnslib.h"
#include "GridType.h"
#include "Geo_FaceTopo_Unstruct.h"
#include "Geo_NodeTopo_Unstruct.h"
#include "LinkStruct.h"
#include "PHMatrix.h"
#include "Glb_Dimension.h"
#include "Constants.h"
#include "Math_Limiter.h"
#include "PHIO.h"
#include "IO_HDF5File.h"
#include "HyList.h"
#include "GRHash.h"
#include "TK_Exit.h"
#include "Geo_UnstructBC.h"
#include "Geo_GridManager.h"
#include "Geo_Interface.h"
#include "Geo_Interpoint.h"
#include "Geo_FaceMetrics_Unstruct.h"
#include "Geo_CellMetrics_Unstruct.h"
#include "Geo_DynamicGridMetrics_Unstruct.h"
#include "Geo_LSQWeight_Unstruct.h"
#include "MeshAgglomeration.h"
#include "FieldProxy.h"
#include "Geo_UnstructGrid.h"
#include "Geo_StructGrid.h"
#include "TK_Log.h"
#include "PHHeapSort.h"
#include "../../../CFD/Solver/include/CFDSolver.h"
#ifdef USE_CUDA
#include "BasicDeviceVariables.h"
#include "GPUCompNodeVar.h"
#include "GPUCompGradientGGNode.h"
#include "GPUCompGradientGGCell.h"
#include "GPUFaceWeight.h"
#include "GPUCompGradientGGNodeNew.h"
#include "GPUMPICommunication.h"
#ifdef CUDAUNITTEST
#include "GPUKernelTestPart2.h"
#endif
#endif

using namespace std;
#pragma warning (disable:4100)
#pragma warning (disable:913)
namespace PHSPACE
{
struct linkbase
{
    int ic;
    struct linkbase *next;
};

Connectivity::Connectivity()
{
    data  = 0;
    index = 0;
    n     = 0;
}

Connectivity::~Connectivity()
{
    delete [] data;    data = nullptr;
    delete [] index;    index = nullptr;
}

UnstructGrid * UnstructGridCast(Grid *grid_in)
{
    if (!grid_in) return nullptr;

    if (grid_in->Type() == PHSPACE::STRUCTGRID)
    {
        StructGrid *strGrid = StructGridCast(grid_in);
        return strGrid->GetUnstrGrid();
    }
    return static_cast<UnstructGrid *>(grid_in);
}

LIB_EXPORT UnstructGrid::UnstructGrid()
{
    wallCellLabel = nullptr;
    unstructBCSet = nullptr;
    cell2coarsegridcell  = 0;
    averageVolume = -1.0;
    zoneProbesCellID.resize(0);
    normalDistanceC2C = nullptr;
    lineTopology = nullptr;
    gridManager = nullptr;
    iBlank = nullptr;
    walldistNode = nullptr;
    JacOrder    =   1;    //! GMRESJac1st
    nearestwallfacenormalx = nullptr;
    nearestwallfacenormaly = nullptr;
    nearestwallfacenormalz = nullptr;

#ifdef USE_INCOMSOLVER
    cellCenterVector = nullptr;
    faceNormalVector = nullptr;
    localCell2GlobalMap = nullptr;
    localCell2LocalMap = nullptr;
    cell2ZoneIDMap = nullptr;
#endif

    ordinaryGridIndex = -1;
    ordinaryNodeIndex = nullptr;
    ordinaryFaceIndex = nullptr;
    ordinaryCellIndex = nullptr;
}

LIB_EXPORT UnstructGrid::~UnstructGrid()
{
    if (nodeTopology != nullptr)
    {
        delete nodeTopology;
        nodeTopology = nullptr;
    }
    if (faceTopology != nullptr)
    {
        delete faceTopology;
        faceTopology = nullptr;
    }
    if (cellTopology != nullptr)
    {
        delete cellTopology;
        cellTopology = nullptr;
    }

    if (faceMetrics != nullptr)
    {
        delete faceMetrics;
        faceMetrics = nullptr;
    }
    if (cellMetrics != nullptr)
    {
        delete cellMetrics;
        cellMetrics = nullptr;
    }
    if (dynamicGridMetrics != nullptr)
    {
        delete dynamicGridMetrics;
        dynamicGridMetrics = nullptr;
    }
    if (leastSquareWeights != nullptr)
    {
        delete leastSquareWeights;
        leastSquareWeights = nullptr;
    }

    if (gridManager != nullptr)
    {
        delete gridManager;
        gridManager = nullptr;
    }

    if (iBlank != NULL)
    {
        delete [] iBlank;
        iBlank = NULL;
    }

    delete cell2cell_connectivity;

    DelPointer(walldist);

    if (nearestwallfacenormalx != nullptr)
    {
        DelPointer(nearestwallfacenormalx);
        nearestwallfacenormalx = nullptr;
    }
   
    if (nearestwallfacenormaly != nullptr)
    {
        DelPointer(nearestwallfacenormaly);
        nearestwallfacenormaly = nullptr;
    }

    if (nearestwallfacenormalz != nullptr)
    {
        DelPointer(nearestwallfacenormalz);
        nearestwallfacenormalz = nullptr;
    }

    DelPointer(walldistNode);

    delete [] lamdax;    lamdax = nullptr;
    delete [] lamday;    lamday = nullptr;
    delete [] lamdaz;    lamdaz = nullptr;
    delete [] knode;    knode = nullptr;

    if (bcr && IsFinestGrid())
    {
        int nBoundFace = this->GetNBoundFace();
        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            delete bcr[iFace];
        }
        delete [] bcr;    bcr = nullptr;
    }

    boundaryPointLabel.clear();

    if (wallCellLabel != nullptr)
    {
        delete [] wallCellLabel;    wallCellLabel = nullptr;
    }

    delete [] cell2coarsegridcell;    cell2coarsegridcell = nullptr;

    if (unstructBCSet != nullptr && IsFinestGrid())
    {
        delete unstructBCSet;
        unstructBCSet = nullptr;
    }
    zoneProbesCellID.clear();

    if (normalDistanceC2C != nullptr)
    {
        delete [] normalDistanceC2C;
        normalDistanceC2C = nullptr;
    }
    if (lineTopology != nullptr)
    {
        delete lineTopology;
        lineTopology = nullptr;
    }
#ifdef USE_INCOMSOLVER
    //! This is values for SIMPLE algorithm.
    //! cellCenterVector container the three coordinates of the cells center. 
    if (cellCenterVector != nullptr)
    {
        delete cellCenterVector;
        cellCenterVector = nullptr;
    }

    if (faceNormalVector != nullptr)
    {
        delete faceNormalVector;
        faceNormalVector = nullptr;
    }

    if (localCell2GlobalMap != nullptr)
    {
        delete localCell2GlobalMap;
        localCell2GlobalMap = nullptr;
    }
    if (localCell2LocalMap != nullptr)
    {
        delete localCell2LocalMap;
        localCell2LocalMap = nullptr;
    }

    if (cell2ZoneIDMap != nullptr)
    {
        delete cell2ZoneIDMap;
        cell2ZoneIDMap = nullptr;
    }
#endif

    if (ordinaryNodeIndex != nullptr)
    {
        delete [] ordinaryNodeIndex;
        ordinaryNodeIndex = nullptr;
    }

    if (ordinaryFaceIndex != nullptr)
    {
        delete [] ordinaryFaceIndex;
        ordinaryFaceIndex = nullptr;
    }

    if (ordinaryCellIndex != nullptr)
    {
        delete [] ordinaryCellIndex;
        ordinaryCellIndex = nullptr;
    }
}

LIB_EXPORT void UnstructGrid::InitGrid(GridID *index, int level, int dim, int type)
{
    Grid::InitGrid(index, level, dim, UNSTRUCTGRID);
    nodeTopology = new Geo_NodeTopo_Unstruct();
    lineTopology = new Geo_LineTopo_Unstruct();
    faceTopology = new Geo_FaceTopo_Unstruct();
    cellTopology = new Geo_CellTopo_Unstruct();

    faceMetrics = new Geo_FaceMetrics_Unstruct();
    cellMetrics = new Geo_CellMetrics_Unstruct();
    dynamicGridMetrics = new Geo_DynamicGridMetrics_Unstruct();
    leastSquareWeights = new Geo_LSQWeight_Unstruct();

    bcr = nullptr;

    cell2cell_connectivity = nullptr;
    fc2cL = nullptr;
    fc2cR = nullptr;
    walldist = nullptr;

    lamdax = nullptr;
    lamday = nullptr;
    lamdaz = nullptr;
    knode  = nullptr;

    iBlank = nullptr;
}

LIB_EXPORT int * UnstructGrid::GetNodeNumberOfEachCell()
{
    if (!cellTopology->GetNodeNumberOfEachCell())
    {
        bool err = ComputeCellNodeTopo();
        if (!err)
        {
            ComputeCell2Node();
        }
    }
    return cellTopology->GetNodeNumberOfEachCell();
}

LIB_EXPORT int * UnstructGrid::GetCell2Node()
{
    if (!cellTopology->GetCell2Node()) 
    {
        bool err = ComputeCellNodeTopo();
        if (!err)
        {
            ComputeCell2Node();
        }
    }
    return cellTopology->GetCell2Node();
}

LIB_EXPORT bool UnstructGrid::ExistCell2Node()
{
    if (!cellTopology) return false;

    if (!cellTopology->GetCell2Node())
    {
        return false;
    }
    return true;
}

LIB_EXPORT int ** UnstructGrid::GetCell2NodeArray()
{
    if (!cellTopology->GetCell2NodeArray()) 
    {
        ComputeCell2NodeArray();
    }
    return cellTopology->GetCell2NodeArray();
}

LIB_EXPORT int ** UnstructGrid::GetFace2NodeArray()
{
    if (!faceTopology->GetFace2NodeArray()) 
    {
        ComputeFace2NodeArray();
    }
    return faceTopology->GetFace2NodeArray();
}

LIB_EXPORT long long int * UnstructGrid::GetFace2NodeSubscript()
{
    if (!faceTopology->GetFace2NodeSubscript())
    {
        ComputeFace2NodeSubscript();
    }
    return faceTopology->GetFace2NodeSubscript();
}

LIB_EXPORT bool UnstructGrid::IsCell2NodeExist()
{
    if (cellTopology->GetCell2Node()) 
    {
        return true;
    }else
    {
        return false;
    }
}

int * UnstructGrid::GetWallCellLabel()
{
    if (wallCellLabel == nullptr)
    {
        int *leftCellOfFace = this->GetLeftCellOfFace();
        wallCellLabel = new int [nTotalCell];

        for (int iCell = 0; iCell < nTotalCell; ++iCell)
        {
            wallCellLabel[iCell] = 0;
        }

        UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
        int nBCRegion = unstructBCSet->GetnBCRegion();

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
            int bcType = bcRegion->GetBCType();
            if (IsWall(bcType))
            {
                vector<int> *faceIndex = bcRegion->GetFaceIndex();
                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    //! iFace is the face number in the set of faceIndex.
                    int iFace = *iter;
                    int le = leftCellOfFace[ iFace ];
                    wallCellLabel[le] = 1;
                }
            }
        }
    }
    return wallCellLabel;
}

void UnstructGrid::InitWallDist()
{
    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotal = nBoundFace + nTotalCell;

    this->walldist = new RDouble [nTotal];
    this->nearestwallfacenormalx = new RDouble[nTotal];
    this->nearestwallfacenormaly = new RDouble[nTotal];
    this->nearestwallfacenormalz = new RDouble[nTotal];
}

ElementType UnstructGrid::CellFaceVertex(int  icel,
                                              int *cell_vtx_edge,
                                              int *cell_vtx_tria,
                                              int *cell_vtx_quad)
{
    int  *nodeNumberOfEachFace = this->GetNodeNumberOfEachFace();
    int **face2node            = this->GetFace2NodeArray();
    int  *faceNumberOfEachCell = this->GetFaceNumberOfEachCell();
    int **cell2face            = this->GetCell2Face();
    int *leftCellOfFace        = this->GetLeftCellOfFace();

    int n_edges = 0;
    int n_trias = 0;
    int n_quads = 0;
    int nface   = faceNumberOfEachCell[icel];

    for (int iface = 0; iface < nface; ++iface)
    {
        int face_id = cell2face[icel][iface];
        int nvertex = nodeNumberOfEachFace[face_id];
        if (nvertex == 2 && n_edges < 4)
        {
            if (cell_vtx_edge)
            {
                int leftindex = leftCellOfFace[face_id];
                if (leftindex == icel)
                {
                    int nodeCount = 0;
                    for (int vtx = 0; vtx < nvertex; vtx ++)
                    {
                        cell_vtx_edge[n_edges * 2 + nodeCount] = face2node[face_id][vtx];
                        nodeCount ++;
                    }
                }
                else
                {
                    int nodeCount = 0;
                    for (int vtx = nvertex - 1; vtx >= 0; vtx --)
                    {
                        cell_vtx_edge[n_edges * 2 + nodeCount] = face2node[face_id][vtx];
                        nodeCount ++;
                    }
                }
            }
            n_edges += 1;
        }
        else if (nvertex == 3 && n_trias < 4)
        {
            if (cell_vtx_tria)
            {
                for (int vtx = 0; vtx < nvertex; vtx++)
                    cell_vtx_tria[n_trias * 3 + vtx] = face2node[face_id][vtx];
            }
            n_trias += 1;
        }
        else if (nvertex == 4 && n_quads < 6)
        {
            if (cell_vtx_quad)
            {
                for (int vtx = 0; vtx < nvertex; vtx++)
                    cell_vtx_quad[n_quads * 4 + vtx] = face2node[face_id][vtx];
            }
            n_quads += 1;
        }
    }

    ElementType cell_type = TYPE_NOT_CONSIDERED;
    if (n_edges == 3)
        cell_type = TYPE_FACE_TRIA;
    else if (n_edges == 4)
        cell_type = TYPE_FACE_QUAD;
    else if (n_trias == 0 && n_quads == 6)
        cell_type = TYPE_CELL_HEXA;
    else if (n_trias == 2 && n_quads == 3)
        cell_type = TYPE_CELL_PRISM;
    else if (n_trias == 4)
    {
        if (n_quads == 0)
            cell_type = TYPE_CELL_TETRA;
        else if (n_quads == 1)
            cell_type = TYPE_CELL_PYRAM;
    }
    return cell_type;
}

int UnstructGrid::nodal_from_cel_hexa(const int *cell_vtx_quad, int *cell_vtx_hexa)
{
    int vertex_id, face_id;
    int ipass, direction;
    int vtx_num, vtx_num_1, vtx_num_2;

    bool warn_orient = false;

    cell_vtx_hexa[0] = cell_vtx_quad[3];
    cell_vtx_hexa[1] = cell_vtx_quad[2];
    cell_vtx_hexa[2] = cell_vtx_quad[1];
    cell_vtx_hexa[3] = cell_vtx_quad[0];

    for (ipass = 0; ipass < 2; ipass++)
    {
        vtx_num_1 = cell_vtx_hexa[ipass * 2];
        vtx_num_2 = cell_vtx_hexa[1 + (ipass * 2)];

        direction = 0;

        for (face_id = 1; face_id < 6; face_id++)
        {
            for (vertex_id = 0; vertex_id < 4; vertex_id++)
            {
                vtx_num = cell_vtx_quad[face_id * 4 + vertex_id];

                if (vtx_num == vtx_num_1)
                {
                    if (cell_vtx_quad[face_id * 4 + ((vertex_id + 1) % 4)] == vtx_num_2)
                    {
                        direction = 1;
                        break;
                    }
                    else if (cell_vtx_quad[face_id * 4 + ((vertex_id - 1 + 4) % 4)] ==
                        vtx_num_2)
                    {
                        direction = -1;
                        break;
                    }
                }
            }

            if (direction != 0)
                break;
        }

        if (direction == -1)
            warn_orient = true;
        else if (direction == 0)
            return -1;

        cell_vtx_hexa[4 + (ipass * 2)] =
            cell_vtx_quad[face_id * 4 + ((vertex_id + (4 + (3 * direction))) % 4)];

        cell_vtx_hexa[5 + (ipass * 2)] =
            cell_vtx_quad[face_id * 4 + ((vertex_id + (4 + (2 * direction))) % 4)];
    }

    if (warn_orient == true)
        return 0;
    else
        return 1;
}

int UnstructGrid::nodal_from_cel_prism(const int *cell_vtx_tria, const int *cell_vtx_quad, int *cell_vtx_prism)
{
    int vertex_id, face_id;
    int ipass, direction;
    int vtx_num, vtx_num_1, vtx_num_2;

    bool warn_orient = false;

    cell_vtx_prism[0] = cell_vtx_tria[2];
    cell_vtx_prism[1] = cell_vtx_tria[1];
    cell_vtx_prism[2] = cell_vtx_tria[0];

    for (ipass = 0; ipass < 2; ipass++)
    {
        vtx_num_1 = cell_vtx_prism[ipass];
        vtx_num_2 = cell_vtx_prism[1 + ipass];

        direction = 0;

        for (face_id = 0; face_id < 4; face_id++)
        {
            for (vertex_id = 0; vertex_id < 4; vertex_id++)
            {
                vtx_num = cell_vtx_quad[face_id * 4 + vertex_id];

                if (vtx_num == vtx_num_1)
                {
                    if (cell_vtx_quad[face_id * 4 + ((vertex_id + 1) % 4)] == vtx_num_2)
                    {
                        direction = 1;
                        break;
                    }
                    else if (cell_vtx_quad[face_id * 4 + ((vertex_id - 1 + 4) % 4)] ==
                        vtx_num_2)
                    {
                        direction = -1;
                        break;
                    }
                }
            }

            if (direction != 0)
                break;
        }

        if (direction == -1)
            warn_orient = true;
        else if (direction == 0)
            return -1;

        cell_vtx_prism[3 + ipass] =
            cell_vtx_quad[face_id * 4 + ((vertex_id + (4 + (3 * direction))) % 4)];

        cell_vtx_prism[4 + ipass] =
            cell_vtx_quad[face_id * 4 + ((vertex_id + (4 + (2 * direction))) % 4)];
    }

    if (warn_orient == true)
        return 0;
    else
        return 1;
}

int UnstructGrid::nodal_from_cel_pyram(const int *cell_vtx_tria, const int *cell_vtx_quad, int *cell_vtx_pyram)
{
    int vertex_id, face_id;
    int direction;
    int vtx_num, vtx_num_1, vtx_num_2;

    bool warn_orient = false;

    cell_vtx_pyram[0] = cell_vtx_quad[3];
    cell_vtx_pyram[1] = cell_vtx_quad[2];
    cell_vtx_pyram[2] = cell_vtx_quad[1];
    cell_vtx_pyram[3] = cell_vtx_quad[0];

    vtx_num_1 = cell_vtx_pyram[0];
    vtx_num_2 = cell_vtx_pyram[1];

    direction = 0;

    for (face_id = 0; face_id < 4; face_id++)
    {
        for (vertex_id = 0; vertex_id < 3; vertex_id++)
        {
            vtx_num = cell_vtx_tria[face_id * 3 + vertex_id];

            if (vtx_num == vtx_num_1)
            {
                if (cell_vtx_tria[face_id * 3 + ((vertex_id + 1) % 3)] == vtx_num_2)
                {
                    direction = 1;
                    break;
                }
                else if (cell_vtx_tria[face_id * 3 + ((vertex_id - 1 + 3) % 3)] == vtx_num_2)
                {
                    direction = -1;
                    break;
                }
            }
        }

        if (direction != 0)
            break;
    }

    if (direction == -1)
        warn_orient = true;
    else if (direction == 0)
        return -1;

    cell_vtx_pyram[4] = cell_vtx_tria[face_id * 3 + ((vertex_id + (3 + (2 * direction))) % 3)];

    if (warn_orient == true)
        return 0;
    else
        return 1;
}

int UnstructGrid::nodal_from_cel_tetra(const int *cell_vtx_tria, int *cell_vtx_tetra)
{
    int vertex_id, face_id;
    int direction;
    int vtx_num, vtx_num_1, vtx_num_2;

    bool warn_orient = false;

    cell_vtx_tetra[0] = cell_vtx_tria[2];
    cell_vtx_tetra[1] = cell_vtx_tria[1];
    cell_vtx_tetra[2] = cell_vtx_tria[0];

    vtx_num_1 = cell_vtx_tetra[0];
    vtx_num_2 = cell_vtx_tetra[1];

    direction = 0;

    for (face_id = 1; face_id < 4; face_id++)
    {
        for (vertex_id = 0; vertex_id < 3; vertex_id++)
        {
            vtx_num = cell_vtx_tria[face_id * 3 + vertex_id];
            if (vtx_num == vtx_num_1)
            {
                if (cell_vtx_tria[face_id * 3 + ((vertex_id + 1) % 3)] == vtx_num_2)
                {
                    direction = 1;
                    break;
                }
                else if (cell_vtx_tria[face_id * 3 + ((vertex_id - 1 + 3) % 3)] == vtx_num_2)
                {
                    direction = -1;
                    break;
                }
            }
        }
        if (direction != 0)
            break;
    }

    if (direction == -1)
        warn_orient = true;
    else if (direction == 0)
        return -1;

    cell_vtx_tetra[3] = cell_vtx_tria[face_id * 3 + ((vertex_id + (3 + (2 * direction))) % 3)];

    if (warn_orient == true)
        return 0;
    else
        return 1;
}

int UnstructGrid::nodal_from_face_quad(const int *cell_vtx_edge, int *cell_vtx_quad)
{
    int face_id;
    int vtx_num_2;

    cell_vtx_quad[0] = cell_vtx_edge[0];
    cell_vtx_quad[1] = cell_vtx_edge[1];

    vtx_num_2     = cell_vtx_quad[1];
    int k         = 0;
    int face_id_s = 0;
    do
    {
        for (face_id = 1; face_id < 4; face_id++)
        {
            int vtx_num_t_1 = cell_vtx_edge[face_id * 2 + 0];
            int vtx_num_t_2 = cell_vtx_edge[face_id * 2 + 1];

            if (vtx_num_2 == vtx_num_t_1 && face_id != face_id_s)
            {
                cell_vtx_quad[k + 2] = vtx_num_t_2;
                vtx_num_2            = cell_vtx_quad[k + 2];
                k += 1;
                face_id_s = face_id;
                break;
            }
            else if (vtx_num_2 == vtx_num_t_2 && face_id != face_id_s)
            {
                cell_vtx_quad[k + 2] = vtx_num_t_1;
                vtx_num_2            = cell_vtx_quad[k + 2];
                k += 1;
                face_id_s = face_id;
                break;
            }
        }
    }
    while (k < 2);
    return 1;
}

int UnstructGrid::nodal_from_face_tria(const int *cell_vtx_edge, int *cell_vtx_tria)
{
    int vertex_id, face_id;
    int direction;
    int vtx_num, vtx_num_1, vtx_num_2;

    cell_vtx_tria[0] = cell_vtx_edge[0];
    cell_vtx_tria[1] = cell_vtx_edge[1];

    vtx_num_1 = cell_vtx_tria[0];
    vtx_num_2 = cell_vtx_tria[1];
    direction = 0;

    for (face_id = 1; face_id < 3; face_id++)
    {
        for (vertex_id = 0; vertex_id < 2; vertex_id++)
        {
            vtx_num = cell_vtx_edge[face_id * 2 + vertex_id];
            if (vtx_num == vtx_num_1 || vtx_num == vtx_num_2)
            {
                cell_vtx_tria[2] = cell_vtx_edge[face_id * 2 + (vertex_id + 1) % 2];
                break;
            }
        }
    }

    return 1;
}

bool UnstructGrid::ComputeCellNodeTopo()
{
    if (this->GetLevel() > 0)
    {
        return false;
    }

    int  nTotalCell = this->GetNTotalCell();

    //! Compute node number of each cell.
    int *cellNodeNumber = new int[nTotalCell];
    for (int iCell = 0; iCell < nTotalCell; iCell ++)
    {
        ElementType ele_t = CellFaceVertex(iCell, nullptr, nullptr, nullptr);
        switch (ele_t)
        {
            case TYPE_FACE_TRIA:
                cellNodeNumber[iCell] = 3;
                break;
            case TYPE_FACE_QUAD:
                cellNodeNumber[iCell] = 4;
                break;
            case TYPE_CELL_TETRA:
                cellNodeNumber[iCell] = 4;
                break;
            case TYPE_CELL_PYRAM:
                cellNodeNumber[iCell] = 5;
                break;
            case TYPE_CELL_PRISM:
                cellNodeNumber[iCell] = 6;
                break;
            case TYPE_CELL_HEXA:
                cellNodeNumber[iCell] = 8;
                break;
            case TYPE_NOT_CONSIDERED:
                WriteLogFile("\n  Warning: The cell type in grid is not commom, do not dump flow field in Ensight and Paraview file!\n");
                PrintToWindow("\n  Warning: The cell type in grid is not commom, do not dump flow field in Ensight and Paraview file!\n");
                delete [] cellNodeNumber;    cellNodeNumber = nullptr;
                return false;
                break;
            default:
                cellNodeNumber[iCell] = 0;
                break;
        }
    }
    cellTopology->SetNodeNumberOfEachCell(cellNodeNumber);

    //! Compute cell 2 node array
    int **cellNodeIndexContainer = NewPointer2<int>(nTotalCell, cellNodeNumber);
    int retcode;
    int cell_vtx_edge[2 * 4];
    int cell_vtx_tria[3 * 4];
    int cell_vtx_quad[4 * 6];
    for (int iCell = 0; iCell < nTotalCell; iCell ++)
    {
        ElementType ele_t = CellFaceVertex(iCell, cell_vtx_edge, cell_vtx_tria, cell_vtx_quad);
        switch (ele_t)
        {
            case TYPE_FACE_TRIA:
                retcode = nodal_from_face_tria(cell_vtx_edge, cellNodeIndexContainer[iCell]);
                break;
            case TYPE_FACE_QUAD:
                retcode = nodal_from_face_quad(cell_vtx_edge, cellNodeIndexContainer[iCell]);
                break;
            case TYPE_CELL_TETRA:
                retcode = nodal_from_cel_tetra(cell_vtx_tria, cellNodeIndexContainer[iCell]);
                break;
            case TYPE_CELL_PYRAM:
                retcode = nodal_from_cel_pyram(cell_vtx_tria, cell_vtx_quad, cellNodeIndexContainer[iCell]);
                break;
            case TYPE_CELL_PRISM:
                retcode = nodal_from_cel_prism(cell_vtx_tria, cell_vtx_quad, cellNodeIndexContainer[iCell]);
                break;
            case TYPE_CELL_HEXA:
                retcode = nodal_from_cel_hexa(cell_vtx_quad, cellNodeIndexContainer[iCell]);
                break;
            default:
                retcode = 0;
                break;
        }
    }
    cellTopology->SetCell2NodeArray(cellNodeIndexContainer);

    //! Compute cell 2 node.
    int nodeCount = 0;
    for (int iCell = 0; iCell < nTotalCell; iCell ++)
    {
        nodeCount += cellNodeNumber[iCell];
    }

    int *cell2node = new int[nodeCount];
    nodeCount = 0;
    for (int iCell = 0; iCell < nTotalCell; iCell ++)
    {
        for (int iNode = 0; iNode < cellNodeNumber[iCell]; ++ iNode)
        {
            cell2node[nodeCount ++] = cellNodeIndexContainer[iCell][iNode];
        }
    }
    cellTopology->SetCell2Node(cell2node);

    return true;
}

void UnstructGrid::ComputeCell2NodeArray()
{
    if (!cellTopology->GetCell2Node()) 
    {
        bool err = ComputeCellNodeTopo();
        if (!err)
        {
            ComputeCell2Node();
        }
        else
        {
            return;
        }
    }

    int  nTotalCell     = this->GetNTotalCell();
    int *cell2node      = this->GetCell2Node();
    int *cellNodeNumber = this->GetNodeNumberOfEachCell();

    int   nodeCount              = 0;
    int **cellNodeIndexContainer = NewPointer2<int>(nTotalCell, cellNodeNumber);
    for (int iCell = 0; iCell < nTotalCell; iCell ++)
    {
        for (int iNode = 0; iNode < cellNodeNumber[iCell]; ++ iNode)
        {
            cellNodeIndexContainer[iCell][iNode] = cell2node[nodeCount ++];
        }
    }
    cellTopology->SetCell2NodeArray(cellNodeIndexContainer);
}

void UnstructGrid::ComputeCell2Node()
{
    set<int> *cell2node_set;
    ComputeCell2Node(cell2node_set);

    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotal = nBoundFace + nTotalCell;

    uint_t count = 0;
    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        count += cell2node_set[iCell].size();
    }

    int *cell2node = new int [count];

    int *npc = new int [nTotal];

    count = 0;
    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        npc[iCell] = static_cast<int>(cell2node_set[iCell].size());
        for (set<int>::iterator iter = cell2node_set[iCell].begin(); iter != cell2node_set[iCell].end(); ++ iter)
        {
            cell2node[count ++] = (*iter);
        }
    }

    cellTopology->SetNodeNumberOfEachCell(npc);
    cellTopology->SetCell2Node(cell2node);

    delete [] cell2node_set;    cell2node_set = nullptr;
}

void UnstructGrid::ComputeCell2Node(set<int> *&cell2nodeSet)
{
    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotalFace = this->GetNTotalFace();
    int nTotal = nBoundFace + nTotalCell;

    int *leftCellOfFace  = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    int *nodeNumberOfEachFace = this->GetNodeNumberOfEachFace();
    int *face2node = this->GetFace2Node();

    cell2nodeSet = new set<int>[nTotal];

    int nodePosition = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [ iFace ];
        int re = rightCellOfFace[ iFace ];
        for (int iNode = 0; iNode < nodeNumberOfEachFace[iFace]; ++ iNode)
        {
            if (le >= 0 && le < nTotalCell)
            {
                cell2nodeSet[le].insert(face2node[nodePosition + iNode]);
            }

            if (re >= 0 && re < nTotalCell)
            {
                //! The re is in the interior field.
                cell2nodeSet[re].insert(face2node[nodePosition + iNode]);
            }
        }
        nodePosition += nodeNumberOfEachFace[iFace];
    }
}

void UnstructGrid::ComputeFace2NodeArray()
{
    int numberOfFaces = GetNTotalFace();

    int * nodeNumberOfEachFace = faceTopology->GetNodeNumberOfEachFace();
    int * face2Node            = faceTopology->GetFace2Node();
    int **face2NodeArray       = faceTopology->GetFace2NodeArray();

    face2NodeArray = new int * [numberOfFaces];

    face2NodeArray[0] = face2Node;
    for (int iFace = 1; iFace < numberOfFaces; ++ iFace)
    {
        face2NodeArray[iFace] = &(face2NodeArray[iFace-1][nodeNumberOfEachFace[iFace-1]]);
    }

    faceTopology->SetFace2NodeArray(face2NodeArray);
}


void UnstructGrid::ComputeFace2NodeSubscript()
{
    int numberOfFaces = GetNTotalFace();

    int *nodeNumberOfEachFace = faceTopology->GetNodeNumberOfEachFace();
    long long int *face2NodeSubscript   = faceTopology->GetFace2NodeSubscript();

    //! Considering that it is often used as follows.
    //for (iFace = 0; iFace < nTotalFace; ++iFace){
    //    g = f(face2NodeSubscript[iFace]) + f(face2NodeSubscript[iFace+1]);
    //}
    //! So we set the array length is numberOfFaces + 1.

    face2NodeSubscript = new long long int[numberOfFaces + 1]; 

    face2NodeSubscript[0] = 0;
    for (int iFace = 1; iFace <= numberOfFaces; ++ iFace)
    {
        face2NodeSubscript[iFace] = face2NodeSubscript[iFace-1] + nodeNumberOfEachFace[iFace-1];
    }

    faceTopology->SetFace2NodeSubscript(face2NodeSubscript);
}

LIB_EXPORT vector<int> * UnstructGrid::GetCell2Cell()
{
    if (!cellTopology) return nullptr;

    if (!cellTopology->GetCell2Cell())
    {
        ComputeCell2Cell();
    }

    return cellTopology->GetCell2Cell();
}
#ifdef USE_GMRESSOLVER
//! GMRESVis
LIB_EXPORT vector<int> * UnstructGrid::GMRESGetNeighborCells()
{
    if (!cellTopology) return nullptr;
    return cellTopology->GMRESGetNeighborCells();
}

//! GMRESVis
LIB_EXPORT vector<int> * UnstructGrid::GMRESGetNeighborFaces()
{
    if (!cellTopology) return nullptr;
    return cellTopology->GMRESGetNeighborFaces();
}

//! GMRESVis
LIB_EXPORT vector<int> * UnstructGrid::GMRESGetNeighborLR()
{
    if (!cellTopology) return nullptr;
    return cellTopology->GMRESGetNeighborLR();
}
#endif

void UnstructGrid::ComputeCell2Cell()
{
    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotalFace = this->GetNTotalFace();
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    vector<int> *cell2cell = new vector<int>[ nTotalCell ];
    int le,re;

    //! If boundary is an INTERFACE, need to count ghost cell
    UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();
    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");
    int taskType = GetTaskCode();
    if (periodicType != NO_PERIODICITY && taskType != PARTITION_GRID)
    {
        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++iBCRegion)
        {
            int taskType = GetTaskCode();

            if (taskType != PARTITION_GRID)
            {
                UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
                int bcType = bcRegion->GetBCType();

                if (IsInterface(bcType))
                {
                    vector<int>* faceIndex = bcRegion->GetFaceIndex();
                    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                    {
                        //! iFace is the face number in the set of faceIndex.
                        int iFace = *iter;
                        le = leftCellOfFace[iFace];
                        re = iFace + nTotalCell;
                        cell2cell[le].push_back(re);
                    }
                }
            }
        }
    }
    //! Interior faces
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        cell2cell[le].push_back(re);
        cell2cell[re].push_back(le);
    }

    cellTopology->SetCell2Cell(cell2cell);
}

#ifdef USE_GMRESSOLVER
//! GMRESVis
void UnstructGrid::ComputeNeighborinfo()
{

    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotal     = nTotalCell + nBoundFace;
    int nTotalFace = this->GetNTotalFace();
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    vector<int> *neighborCells = new vector<int>[ nTotal ];
    vector<int> *neighborFaces = new vector<int>[ nTotal ];
    vector<int> *neighborLR    = new vector<int>[ nTotal ];

    //! GMRES2ndCorrection
    vector<int> BCLeftCells;
    vector<int> BCRightCells;
    vector<int> BCFaces;

    int le,re;

    //! Boundary faces
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le      = leftCellOfFace[ iFace ];
        re      = iFace + nTotalCell;
        neighborCells[le].push_back(re);
        neighborFaces[le].push_back(iFace);
        neighborLR[le].push_back(1);

        neighborCells[re].push_back(le);
        neighborFaces[re].push_back(iFace);
        neighborLR[re].push_back(-1);

        //! GMRES2ndCorrection
        BCLeftCells.push_back(le);
        BCRightCells.push_back(re);
        BCFaces.push_back(iFace);
    }

    //! Interior faces
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        neighborCells[le].push_back(re);
        neighborFaces[le].push_back(iFace);
        neighborLR[le].push_back(1);

        neighborCells[re].push_back(le);
        neighborFaces[re].push_back(iFace);
        neighborLR[re].push_back(-1);
    }

    cellTopology->GMRESSetNeighborCells(neighborCells);
    cellTopology->GMRESSetNeighborFaces(neighborFaces);
    cellTopology->GMRESSetNeighborLR(neighborLR);

    //! GMRES2ndCorrection
    this->SetBCLeftCells(BCLeftCells);
    this->SetBCRightCells(BCRightCells);
    this->SetBCFaces(BCFaces);
}
#endif
LIB_EXPORT int * UnstructGrid::GetFaceNumberOfEachCell()
{
    if (!cellTopology->GetFaceNumberOfEachCell())
    {
        ComputeCell2Face();
    }
    return cellTopology->GetFaceNumberOfEachCell();
}

LIB_EXPORT int ** UnstructGrid::GetCell2Face()
{
    if (!cellTopology->GetCell2Face())
    {
        ComputeCell2Face();
    }

    return cellTopology->GetCell2Face();
}


void UnstructGrid::ComputeCell2Face()
{
    int **cell2face = cellTopology->GetCell2Face();
    
    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotalFace = this->GetNTotalFace();

    int *leftCellOfFace  = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    vector<int> *vcell2face = new vector<int> [ nTotalCell ];
    int le,re;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace[ iFace ];
        vcell2face[le].push_back(iFace);
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [ iFace ];
        re = rightCellOfFace[ iFace ];
        vcell2face[le].push_back(iFace);
        vcell2face[re].push_back(iFace);
    }

    int *face_number_of_each_cell = new int [nTotalCell];

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        face_number_of_each_cell[iCell] = static_cast<int>(vcell2face[iCell].size());
    }

    cell2face = NewPointer2<int>(nTotalCell, face_number_of_each_cell); 

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        for (std::size_t j = 0; j < vcell2face[iCell].size(); ++ j)
        {
            cell2face[iCell][j] = vcell2face[iCell][j];
        }
    }

    delete [] vcell2face;    vcell2face = nullptr;

    cellTopology->SetCell2Face(cell2face);
    cellTopology->SetFaceNumberOfEachCell(face_number_of_each_cell);
}

LIB_EXPORT int * UnstructGrid::GetNumberOfNeighborCell() 
{
    if (!cellTopology->GetNumberOfNeighborCell())
    {
        ComputeNumberOfNeighborCell();
    }
    return cellTopology->GetNumberOfNeighborCell(); 
}

//! Calculate number_of_neighbor_cell(number of neighboring cells in each cell without including itself)
void UnstructGrid::ComputeNumberOfNeighborCell()
{
    int nBoundFace = this->GetNBoundFace();
    int nTotalFace = this->GetNTotalFace();
    int nTotalCell = this->GetNTotalCell();
    int *leftCellOfFace  = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();
    
    UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    //! Allocate memories for number of cells per cell
    int *number_of_neighbor_cell = new int[ nTotalCell ];

    //! set number_of_neighbor_cell to be 0 without including itself
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        number_of_neighbor_cell[iCell] = 0;
    }

    int le = 0, re = 0;

    //! If boundary is an INTERFACE, need to count ghost cell
    //! Note: Symmetry???
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (IsInterface(bcType))
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                le      = leftCellOfFace[ iFace ];
                number_of_neighbor_cell[le] ++;
            }
        }
    }
 
    //! Interior faces
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [ iFace ];
        re = rightCellOfFace[ iFace ];
        number_of_neighbor_cell[le] ++;
        number_of_neighbor_cell[re] ++;
    }

    cellTopology->SetNumberOfNeighborCell(number_of_neighbor_cell);
}

LIB_EXPORT Connectivity * UnstructGrid::GetCell2CellConnectivity()
{
    if (!cell2cell_connectivity)
    {
        ComputeCell2CellConnectivity();
    }

    return cell2cell_connectivity;
}

LIB_EXPORT int * UnstructGrid::GetCellNumberOfEachNode()
{
    if (!nodeTopology->GetCellNumberOfEachNode())
    {
        ComputeNode2Cell();
    }
    return nodeTopology->GetCellNumberOfEachNode();
}

LIB_EXPORT int * UnstructGrid::GetNode2Cell()
{
    if (!nodeTopology->GetNode2Cell())
    {
        ComputeNode2Cell();
    }
    return nodeTopology->GetNode2Cell();
}

LIB_EXPORT int ** UnstructGrid::GetNode2CellArray()
{
    if (!nodeTopology->GetNode2CellArray())
    {
        ComputeNode2CellArray();
    }
    return nodeTopology->GetNode2CellArray();

}

LIB_EXPORT RDouble * UnstructGrid::GetLargestLocalGridLength()
{
    if (!cellMetrics->GetLargestLocalGridLength())
    {
        ComputeLargestLocalGridLength();
    }
    return cellMetrics->GetLargestLocalGridLength();
}

LIB_EXPORT RDouble * UnstructGrid::GetSubgridLength()
{
    if (!cellMetrics->GetSubgridLength())
    {
        ComputeSubgridLength();
    }
    return cellMetrics->GetSubgridLength();
}

void UnstructGrid::computeBoundaryPointLabel()
{
    int *face2node = GetFace2Node();
    int *node_number_of_each_face   = this->GetNodeNumberOfEachFace();

    UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    int pt, bcType;
    int nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; iFace++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        bcType    = bcRegion->GetBCType();
        if (bcType != PHENGLEI::INTERFACE)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt = face2node[ nodepos + j ];
                boundaryPointLabel.insert(pair<int,int>(pt,1));
            }
        }        
        nodepos += node_number_of_each_face[iFace];
    }

    nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; iFace++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        bcType    = bcRegion->GetBCType();
        if (bcType == PHENGLEI::SYMMETRY)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt = face2node[ nodepos + j ];
                boundaryPointLabel[pt] = PHENGLEI::SYMMETRY;
            }
        }        
        nodepos += node_number_of_each_face[iFace];
    }
    nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; iFace++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        bcType    = bcRegion->GetBCType();

        if (bcType == PHENGLEI::FARFIELD)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt = face2node[ nodepos + j ];
                boundaryPointLabel[pt] = PHENGLEI::FARFIELD;
            }
        }        
        nodepos += node_number_of_each_face[iFace];
    }

    nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; iFace++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        bcType    = bcRegion->GetBCType();

        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt = face2node[ nodepos + j ];
                boundaryPointLabel[pt] = PHENGLEI::SOLID_SURFACE;
            }
        }        
        nodepos += node_number_of_each_face[iFace];
    }
}

vector <int> & UnstructGrid::GetLocalToGlobalNodesMap()
{
    if (localToGlobalNodesMap.size() != 0) return localToGlobalNodesMap;

    int numberOfNodesNow = GetNTotalNode();

    localToGlobalNodesMap.resize(numberOfNodesNow);

    for (int iNode = 0; iNode < numberOfNodesNow; ++ iNode)
    {
        localToGlobalNodesMap[iNode] = iNode;
    }

    return localToGlobalNodesMap;
}

vector <int> & UnstructGrid::GetBoundaryConditionType()
{
    if (boundaryConditionType.size() != 0) return boundaryConditionType;

    boundaryConditionType.resize(nBoundFace);

    UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcTypeOld = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();

        if (bcTypeOld == PHENGLEI::INTERFACE)
        {
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                boundaryConditionType[iFace] = 0;
            }
        }
        else if (bcTypeOld == PHENGLEI::SOLID_SURFACE)
        {
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                boundaryConditionType[iFace] = 1;
            }
        }
        else if (bcTypeOld == PHENGLEI::OVERSET)
        {
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                boundaryConditionType[iFace] = 2;
            }
        }
        else if (bcTypeOld == PHENGLEI::SYMMETRY)
        {
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                boundaryConditionType[iFace] = 3;
            }
        }
        else 
        {
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                boundaryConditionType[iFace] = 9;
            }
        }
    }

    return boundaryConditionType;
}

//! Compute the minimum distance of the grid edges.
LIB_EXPORT RDouble UnstructGrid::CalMinEdgeLength()
{
    if (this->minEdgeLength > 0.0)
    {
        return this->minEdgeLength;
    }

    minEdgeLength = PHSPACE::LARGE;

    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();
    int *face2node = this->GetFace2Node();

    int nodePosition = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int nNode = node_number_of_each_face[iFace];
        for (int iNode = 0; iNode < nNode; ++ iNode)
        {
            int node1 = face2node[nodePosition + iNode];
            int node2 = face2node[nodePosition + (iNode + 1) % nNode];

            RDouble dx = x[ node2 ] - x[ node1 ];
            RDouble dy = y[ node2 ] - y[ node1 ];
            RDouble dz = z[ node2 ] - z[ node1 ];
            RDouble ds = DISTANCE(dx, dy, dz);

            minEdgeLength = MIN(minEdgeLength, ds);
        }
        
        nodePosition += nNode;
    }

    return minEdgeLength;
}

//! Compute the maximum distance of the grid edges.
LIB_EXPORT RDouble UnstructGrid::CalMaxEdgeLength()
{
    if (this->maxEdgeLength > 0.0)
    {
        return this->maxEdgeLength;
    }

    maxEdgeLength = PHSPACE::SMALL;

    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();
    int *face2node = this->GetFace2Node();

    int nodePosition = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int nNode = node_number_of_each_face[iFace];
        for (int iNode = 0; iNode < nNode; ++ iNode)
        {
            int node1 = face2node[nodePosition + iNode];
            int node2 = face2node[nodePosition + (iNode + 1) % nNode];

            RDouble dx = x[ node2 ] - x[ node1 ];
            RDouble dy = y[ node2 ] - y[ node1 ];
            RDouble dz = z[ node2 ] - z[ node1 ];
            RDouble ds = DISTANCE(dx, dy, dz);

            maxEdgeLength = MAX(maxEdgeLength, ds);
        }

        nodePosition += nNode;
    }

    return maxEdgeLength;
}

LIB_EXPORT void UnstructGrid::ReSpecifyBC()
{
    int nBoundFace = this->GetNBoundFace();
    int nTotalCell = this->GetNTotalCell();
    int *leftCellOfFace  = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        if (rightCellOfFace[iFace] >= nTotalCell)
        {
            rightCellOfFace[iFace] = -1;
        }
        if (leftCellOfFace[iFace] >= nTotalCell)
        {
            leftCellOfFace[iFace] = -1;
        }
    }
}

//! Get node to cell information(cxh, 2012.12.12)
void UnstructGrid::ComputeNode2Cell()
{
    int i,j,node,ip,le,re;
    int nTotalNode = this->GetNTotalNode();
    int nTotalFace = this->GetNTotalFace();
    int *leftCellOfFace  = this->GetLeftCellOfFace();
    int *rightCellOfFace  = this->GetRightCellOfFace();
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();
    int *face2node = this->GetFace2Node();

    struct linkbase **nc;
    struct linkbase *link, *tmp;

    int *nCPN = new int[ nTotalNode ];
    nc = new struct linkbase * [ nTotalNode ];

    for (i = 0; i < nTotalNode; ++ i)
    {
        nCPN[i] = 0;
        nc[i] = nullptr;
    }

    //!firstly, search all the face
    node=0;
    for (i = 0; i < nTotalFace; ++ i)
    {
        le = leftCellOfFace[ i ];
        re = rightCellOfFace[ i ];
        for (j = 0; j < node_number_of_each_face[ i ]; ++ j)
        {
            ip = face2node[ node  ++  ];
            addn2c(le, ip, nCPN, nc);
            addn2c(re, ip, nCPN, nc);
        }
    }

    //! change date from linkbase to array
    long long int count;
    count = 0;
    for(i = 0; i < nTotalNode; ++ i)
    {
        count = count + nCPN[ i ];
    }
    int *n2c = new int [count];

    count = 0;
    for (i = 0; i < nTotalNode; ++ i)
    {
        for (link = nc[ i ]; link != NULL; link = link->next)
        {
            n2c[ count++ ] = link->ic;
        }
    }

    //delete linkbase
    for (i = 0; i < nTotalNode; ++ i)
    {
        link = nc[ i ];
        while (link != nullptr)
        {
            tmp = link->next;
            delete link;
            link = tmp;
        }
    }

    delete [] nc;    nc = nullptr;

    nodeTopology->SetCellNumberOfEachNode(nCPN);
    nodeTopology->SetNode2Cell(n2c);
}

void UnstructGrid::ComputeNode2CellArray()
{
    int nTotalNode = this->GetNTotalNode();    

    int *node2Cell    = this->GetNode2Cell();
    int *nCellPerNode = this->GetCellNumberOfEachNode();

    int **node2CellArray = new int * [nTotalNode];

    node2CellArray[0] = node2Cell;
    for (int iNode = 1; iNode < nTotalNode; ++iNode)
    {
        node2CellArray[iNode] = &(node2CellArray[iNode-1][nCellPerNode[iNode-1]]);
    }

    nodeTopology->SetNode2CellArray(node2CellArray);

    cout<<"Get Node2CellArray, well down!"<<endl;
}

int * UnstructGrid::GetNTotalFacesOfEachNode()
{
    if (!nodeTopology->GetNTotalFacesOfEachNode())
    {
        this->ComputeFace2NodeToNode2Face();
    }
    return nodeTopology->GetNTotalFacesOfEachNode();
}

int * UnstructGrid::GetNode2Face()
{
    if (!nodeTopology->GetNode2Face())
    {
        this->ComputeFace2NodeToNode2Face();
    }
    return nodeTopology->GetNode2Face();
}

void UnstructGrid::ComputeFace2NodeToNode2Face()
{
    int *nTotalFacesOfEachNode = new int [nTotalNode];
    int *nTotalNodesOfEachFace = GetNodeNumberOfEachFace();
    int **face2NodeArray       = GetFace2NodeArray();
    vector < vector < int > > node2FaceArray;
    node2FaceArray.resize(nTotalNode);
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        for (int iNode = 0; iNode < nTotalNodesOfEachFace[iFace]; ++ iNode)
        {
            int nodeID = face2NodeArray[iFace][iNode];
            node2FaceArray[nodeID].push_back(iFace);
        }
    }

    int datasize = 0;
    for (int iNode = 0; iNode < nTotalNode; ++iNode)
    {
        nTotalFacesOfEachNode[iNode] = node2FaceArray.size();
        datasize = datasize + nTotalFacesOfEachNode[iNode];
    }
    int *node2Face = new int [datasize];
    int faceCount = 0;
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        for (int iFace = 0; iFace < nTotalFacesOfEachNode[iNode]; ++ iFace)
        {
            node2Face[faceCount] = node2FaceArray[iNode][iFace];
            ++ faceCount;
        }
    }

    nodeTopology->SetNTotalFacesOfEachNode(nTotalFacesOfEachNode);
    nodeTopology->SetNode2Face(node2Face);
}

void UnstructGrid::ComputeNode2FaceArray()
{
    int nTotalNode = this->GetNTotalNode();

    int *node2Face = this->GetNode2Face();
    int *nTotalFacesOfEachNode = this->GetNTotalFacesOfEachNode();

    int **node2FaceArray = new int * [nTotalNode];

    node2FaceArray[0] = node2Face;
    for (int iNode = 1; iNode < nTotalNode; ++ iNode)
    {
        node2FaceArray[iNode] = &(node2FaceArray[iNode - 1][nTotalFacesOfEachNode[iNode - 1]]);
    }

    nodeTopology->SetNode2FaceArray(node2FaceArray);
}

void UnstructGrid::GetSourceIndex(int iFace,int ipos, int &s)
{
    InterfaceInfo *interfaceInfo = this->GetInterfaceInfo();
    int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();
    int *leftCellOfFace = this->GetLeftCellOfFace();
    s = leftCellOfFace[interFace2BoundaryFace[iFace]];
}

void UnstructGrid::GetSourcePointIndex(int ipoint,int ipos,int &s)
{
    InterpointInformation *interpointInformation = this->GetInterpointInfo();
    int *interPoint2GlobalPoint = interpointInformation->GetInterPoint2GlobalPoint();
    s = interPoint2GlobalPoint[ipoint];
}
void UnstructGrid::GetTargetIndex(int iFace,int ipos, int &t)
{
    InterfaceInfo *interfaceInfo = this->GetInterfaceInfo();
    int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();
    int *rightCellOfFace = this->GetRightCellOfFace();
    t = rightCellOfFace[interFace2BoundaryFace[iFace]];
}

void UnstructGrid::GetTargetPointIndex(int ipoint,int ipos,int &t)
{
    InterpointInformation *interpointInformation = this->GetInterpointInfo();
    int *interPoint2GlobalPoint = interpointInformation->GetInterPoint2GlobalPoint();
    t = interPoint2GlobalPoint[ipoint];
}

void UnstructGrid::SetGhostCellExceptInterface(RDouble * f)
{
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();
    int le, re;
    UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            if (bcType == PHENGLEI::INTERFACE)
            {
                continue;
            }
            le = leftCellOfFace [iFace];
            re = rightCellOfFace[iFace];
            f[re] = f[le];
        }
    }
}

void UnstructGrid::SetGhostCell(RDouble *field)
{
    int nBoundFace = this->GetNBoundFace();
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();
    int le, re;

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        field[re] = field[le];
    }
}

void UnstructGrid::SetGhostCell(RDouble **field, const int &nVariables)
{
    int nBoundFace = this->GetNBoundFace();
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();
    int le, re;

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        for (int iVar = 0; iVar < nVariables; ++ iVar)
        {
            field[iVar][re] = field[iVar][le];
        }
    }
}
void UnstructGrid::SetGhostCell(RDouble *field, RDouble value)
{
    int nBoundFace = this->GetNBoundFace();
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();
    int le, re;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        field[re] = value;
    }
}

RDouble *UnstructGrid::GetVolume(int istep) const
{
    //! This will do for now.
    cout << "Warning: this function may be wrong!" << endl;
    return this->GetCellVolume();
}

void UnstructGrid::Action(ActionKey *actkey)
{
    switch (actkey->action)
    {
        case WRITE_WALL_DIST:
            DumpWalldist(actkey);
            break;
        case COMMCELLCENTERDATA:
            GetCellCenter(actkey);
            break;
        case COMMCELLIBLANK:
            GetCellIBlank(actkey);
            break;
        case SIMPLE_ACTION:
            SimpleAction(actkey);
            break;
        case ALLOCATE_WALL_DIST:
            AllocateWalldist();
            break;
        case TEST_RECONSTRUCTION:
            SkewnessSummary(actkey);
            break;
        default:
            break;
    }
}

void UnstructGrid::SimpleAction(ActionKey *actkey)
{
    RDouble value;
    GetData(actkey->taskname, &value, PHDOUBLE, 1);

    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();
    cdata->Write(&value, sizeof(RDouble));
}

void UnstructGrid::TranslateAction(ActionKey *actkey)
{
    if (actkey->action == COMMCELLCENTERDATA)
    {
        TranslateCellCenter(actkey);
    }
    if (actkey->action == COMMCELLIBLANK)
    {
        TranslateCellIBlank(actkey);
    }
}

streamsize UnstructGrid::TranslateActionLength(ActionKey *actkey)
{
    if (actkey->action == COMMCELLCENTERDATA)
    {
        return TranslateCellCenterLength(actkey);
    }
    if (actkey->action == COMMCELLIBLANK)
    {
        return TranslateCellIBlankLength(actkey);
    }
    return 0;
}

void UnstructGrid::SpecifyRightCellofBC()
{
    int nBoundFace = this->GetNBoundFace();
    int nTotalCell = this->GetNTotalCell();
    int *rightCellOfFace = this->GetRightCellOfFace();
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        rightCellOfFace[iFace] = iFace + nTotalCell;
    }
}

//! Get laplacian weitht of nodes(cxh, 2012.12.12).
void UnstructGrid::CalcLaplacianWeitht()
{
    int iNode,j,count,iCell;
    int nTotalNode = this->GetNTotalNode();
    RDouble xn,yn,zn;
    RDouble xc,yc,zc;
    RDouble *xp = this->GetX();
    RDouble *yp = this->GetY();
    RDouble *zp = this->GetZ();

    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();

    RDouble Rx,Ry,Rz,Ixx,Iyy,Izz,Ixy,Ixz,Iyz,D,weight;

    int *nCPN = GetCellNumberOfEachNode();
    int *n2c  = GetNode2Cell();

    if (nCPN == nullptr || n2c == nullptr) this->ComputeNode2Cell();

    if (xcc == nullptr || xp == nullptr)
    {
        printf("Pease get the cell center information first!\n");
        exit(0);
    }

    lamdax = new RDouble [ nTotalNode ];
    lamday = new RDouble [ nTotalNode ];
    lamdaz = new RDouble [ nTotalNode ];

    count = 0;

    for (iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        Rx = Ry = Rz = Ixx = Iyy = Izz = Ixy = Ixz = Iyz = 0.0;
        xn = xp[ iNode ];
        yn = yp[ iNode ];
        zn = zp[ iNode ];
        for (j = 0; j < nCPN[ iNode ]; ++ j)
        {
            xc = xcc[ n2c[ count ] ];
            yc = ycc[ n2c[ count ] ];
            zc = zcc[ n2c[ count ] ];
            count ++;

            Rx += xc - xn;
            Ry += yc - yn;
            Rz += zc - zn;

            Ixx += (xc - xn) * (xc - xn);
            Iyy += (yc - yn) * (yc - yn);
            Izz += (zc - zn) * (zc - zn);

            Ixy += (xc - xn) * (yc - yn);
            Ixz += (xc - xn) * (zc - zn);
            Iyz += (yc - yn) * (zc - zn);
        }

        if (this->GetDim() == TWO_D)
        {
            D = Ixx *  Iyy - Ixy * Ixy;

            lamdax[ iNode ] = (Ry * Ixy - Rx * Iyy) / D;

            lamday[ iNode ] = (Rx * Ixy - Ry * Ixx) / D;

            lamdaz[ iNode ] = 0.0;
        } 
        else 
        {
            D = Ixx * (Iyy * Izz - Iyz * Iyz) - 
                Ixy * (Ixy * Izz - Ixz * Iyz) +
                Ixz * (Ixy * Iyz - Iyy * Ixz);

            lamdax[ iNode ] = (-Rx * (Iyy * Izz - Iyz * Iyz)
                                +Ry * (Ixy * Izz - Ixz * Iyz)
                                -Rz * (Ixy * Iyz - Iyy * Ixz)) / D;

            lamday[ iNode ] = ( Rx * (Ixy * Izz - Ixz * Iyz)
                                -Ry * (Ixx * Izz - Ixz * Ixz)
                                +Rz * (Ixx * Iyz - Ixy * Ixz)) / D;

            lamdaz[ iNode ] = (-Rx * (Ixy * Iyz - Iyy * Ixz)
                                +Ry * (Ixx * Iyz - Ixy * Ixz)
                                -Rz * (Ixx * Iyy - Ixy * Ixy)) / D;
        }
    }
        
    if (knode == nullptr) knode = new int [nTotalNode];
    for (iNode = 0; iNode < nTotalNode; ++ iNode) knode[iNode] = 0;

    count = 0;
    for (iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        for (j = 0; j < nCPN[iNode]; ++ j)
        {
            iCell = n2c[ count++ ];
            weight = 1 + lamdax[ iNode ] * (xcc[ iCell ] - xp[ iNode ]) + lamday[ iNode ] * (ycc[ iCell ] - yp[ iNode ]) + lamdaz[ iNode ] * (zcc[ iCell ] - zp[ iNode ]);
            if (weight <= 0.0) knode[ iNode ] = 1;
        }
    }

    count = 0;
    for (iNode = 0; iNode < nTotalNode; ++ iNode) 
    {
        if (knode[ iNode ] == 1) count ++;
    }
    printf(" Number of negative weight cells : %d  %f\n", count, static_cast<RDouble>(count)/nTotalNode);
}

void UnstructGrid::ComputeCell2CellConnectivity()
{
    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotalFace = this->GetNTotalFace();
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    cell2cell_connectivity = new Connectivity();

    int *nnpc = new int [nTotalCell];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        nnpc[iCell] = 0;
    }

    int le,re;

    //! If boundary is an INTERFACE, need to count ghost cell.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bctype  = bcr[iFace]->GetKey();
        if (IsInterface(bctype))
        {
            le      = leftCellOfFace[ iFace ];
            nnpc[le] ++;
        }
    }

    //! Interior faces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        nnpc[le] ++;
        nnpc[re] ++;
    }

    int *index = new int [nTotalCell + 1];

    index[0] = 0;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        index[iCell+1] = index[iCell] + nnpc[iCell];
    }

    cell2cell_connectivity->index = index;
    
    cell2cell_connectivity->n = index[nTotalCell];

    cell2cell_connectivity->data = new int [cell2cell_connectivity->n];

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        nnpc[iCell] = 0;
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bctype  = bcr[iFace]->GetKey();
        if (IsInterface(bctype))
        {
            le      = leftCellOfFace[ iFace ];
            re      = iFace + nTotalCell;
            cell2cell_connectivity->data[ index[le] + nnpc[le]++ ] = re;
        }
    }

    //! Interior faces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        cell2cell_connectivity->data[ index[le] + nnpc[le]++ ] = re;
        cell2cell_connectivity->data[ index[re] + nnpc[re]++ ] = le;
    }

    delete [] nnpc;    nnpc = nullptr;
}

LIB_EXPORT RDouble UnstructGrid::SkewnessSummary(ActionKey *actkey)
{
    int angleC[18];
    RDouble dotp1,x1,y1,z1,dis1,dotp2,x2,y2,z2,dis2,angle;
    RDouble minAng = straightAngle;
    RDouble maxAng = 0.0;
    RDouble area_scale = TINY;
    RDouble criticalAngle = 89.9;

    int le, re;

    for (int i = 0; i < 18; ++ i)
    {
        angleC[i] = 0;
    }
    int *leftCellOfFace  = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();
    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();

    RDouble *xfc = this->GetFaceCenterX();
    RDouble *yfc = this->GetFaceCenterY();
    RDouble *zfc = this->GetFaceCenterZ();

    RDouble *nxs = this->GetFaceNormalX();
    RDouble *nys = this->GetFaceNormalY();
    RDouble *nzs = this->GetFaceNormalZ();
    RDouble *ns  = this->GetFaceArea();

    int nTotalFace = this->GetNTotalFace();
    int nBoundFace = this->GetNBoundFace();
    int nTotalCell = this->GetNTotalCell();
    //#define  QUALYTY_DEBUG
#ifdef QUALYTY_DEBUG
    int * negativeCell = new int [nTotalCell];
    SetField(negativeCell, 0, nTotalCell);
    int nNegativeCell = 0;
#endif

    RDouble *cellSkewness = this->GetCellSkewness();
    SetField(cellSkewness, straightAngle, nTotalCell);

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        if (ns[iFace] <= area_scale) continue;
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        x1 = xfc[iFace] - xcc[le];
        y1 = yfc[iFace] - ycc[le];
        z1 = zfc[iFace] - zcc[le];
        dis1  = sqrt(x1*x1 + y1*y1 + z1*z1);
        dotp1 = (nxs[iFace] * x1 + nys[iFace] * y1 + nzs[iFace] * z1) / (dis1 + TINY);
        if (dotp1 >  one) dotp1 =   one;
        if (dotp1 < -one) dotp1 = - one;

        x2 = xcc[re] - xfc[iFace];
        y2 = ycc[re] - yfc[iFace];
        z2 = zcc[re] - zfc[iFace];
        dis2  = sqrt(x2*x2 + y2*y2 + z2*z2);
        dotp2 = (nxs[iFace] * x2 + nys[iFace] * y2 + nzs[iFace] * z2) / (dis2 + TINY);      
        if (dotp2 >  one) dotp2 =   one;
        if (dotp2 < -one) dotp2 = - one;

        angle = MIN(asin(dotp1), asin(dotp2)) * straightAngle / PI;

#ifdef QUALYTY_DEBUG
        if (le == 8516 && re == 8517)
        {
            PrintToWindow("\nle, re ", le, re);
            PrintToWindow("\nx1, y1, z1 ", x1, y1, z1);
            PrintToWindow("\nxfc[iFace] ", xfc[iFace], yfc[iFace], zfc[iFace]);
            PrintToWindow("\nxcc[le] ", xcc[le], ycc[le], zcc[le]);
            PrintToWindow("\nnxs[iFace] ", nxs[iFace], nys[iFace], nzs[iFace]);
            PrintToWindow("\ndis1, dotp1, asin(dotp1) ", dis1, dotp1, asin(dotp1) * straightAngle / PI);
            PrintToWindow("\nx2, y2, z2 ", x2, y2, z2);
            PrintToWindow("\ndis2, dotp2, asin(dotp2) ", dis2, dotp2, asin(dotp2) * straightAngle / PI);
            negativeCell[le] = 1;
            negativeCell[re] = 1;
            nNegativeCell += 2;
        }
#endif

        minAng = MIN(minAng, angle);
        maxAng = MAX(maxAng, angle);

        cellSkewness[le] = MIN(cellSkewness[le], angle);
        cellSkewness[re] = MIN(cellSkewness[re], angle);

        if (angle > criticalAngle) angle = criticalAngle;
        if (angle < -criticalAngle) angle = -criticalAngle;
        if (angle < 0) angle -= ten;
        angleC[static_cast<int>(angle / ten) + 9] ++;
    }

    int nodepos = 0;
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        if (ns[iFace] <= area_scale) continue;
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        x1 = xfc[iFace] - xcc[le];
        y1 = yfc[iFace] - ycc[le];
        z1 = zfc[iFace] - zcc[le];
        dis1  = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
        dotp1 = (nxs[iFace] * x1 + nys[iFace] * y1 + nzs[iFace] * z1) / (dis1 + TINY);
        if (dotp1 >  one) dotp1 =   one;
        if (dotp1 < -one) dotp1 = - one;

        angle = asin(dotp1) * straightAngle / PI;

        nodepos += node_number_of_each_face[iFace];

        minAng = MIN(minAng, angle);
        maxAng = MAX(maxAng, angle);

        cellSkewness[le] = MIN(cellSkewness[le], angle);

        if (angle > criticalAngle) angle = criticalAngle;
        if (angle < -criticalAngle) angle = -criticalAngle;
        if (angle < zero) angle -= 10;
        angleC[static_cast<int>(angle / 10) + 9] ++;
    }
    int total = 0;
    for (int i = 0; i < 18; ++ i)
    {
        total += angleC[i];
    }

    ostringstream oss;
    oss << "Zone Index :" << this->GetZoneID() << "\n";
    oss << "=====================================================" << "\n";
    oss << "Skewness Summary (angle of 90 degrees being the best)" << "\n";
    oss << "  Level = " << this->GetLevel() << "\n";
    oss << "=====================================================" << "\n";

    oss << "  Total number of cells " << nTotalCell << "\n";
    oss << "  Total number of faces " << nTotalFace << "\n";
    oss << "  Total number of checked Faces " << total << "\n";
    oss << "  Face scale used (area_scale) = " << area_scale << "\n";
    oss << "  Total number of faces whose areas <= area_scale " << nTotalFace - total << "\n";

    oss << "  Minimum skewness angle is " << minAng << endl;
    oss << "  Maximum skewness angle is " << maxAng << endl;
    for (int i = 0; i < 18; ++ i)
    {
        oss << "  Face angle between " ;
        oss << setiosflags(ios::right);
        oss << setw(3) << - 90 + i * 10;
        oss << " and ";
        oss << setiosflags(ios::right);
        oss << setw(3) << - 80 + i * 10;
        oss << " is ";

        oss << resetiosflags(ios_base::scientific);
        oss << setiosflags(ios::fixed);
        oss << setw(5);
        oss << setprecision(2);
        oss <<  angleC[i] / (static_cast<RDouble>(total) + SMALL) * 100;
        oss <<  " Percent (";
        oss << setw(8)<< angleC[i] << ")\n";
    }  
    oss << "====================================================" << "\n";

    //! It is necessary to reset, otherwise there will be some unexpected output elsewhere.
    //! This is due to the setting of cout << setiosflags(ios::fixed),the specific reason are to be considered later.
    oss << resetiosflags(ios::fixed);
    oss << setiosflags(ios_base::scientific);

    PrintToWindow(oss);

    RDouble globalMinAng = 0.0;
    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if (PHMPI::GetNumberOfProcessor() == PHMPI::GetNumberofGlobalZones() && UNSTRUCTGRID == systemGridType)
    {
        //! "UNSTRUCTGRID == systemGridType" is added by myk 20200922, for parral-running with MixGrid.
        //! SkewnessSummary() in StructGrid zone is an empty function, no PH_Reduce() corresponding with the following one.

        //PH_CompareMaxMin(minAng, 2);
        PH_Reduce(&minAng, &globalMinAng, 1, MPI_MIN);

        if (PHMPI::GetCurrentProcessorID() == PHMPI::GetServerProcessorID())
        {
            ostringstream oss1;
            oss1 << resetiosflags(ios_base::scientific);
            oss1 << setprecision(4);
            oss1  << "Global minimum skewness angle is " << globalMinAng << endl;
            PrintToWindow("    Global minimum skewness angle is ", globalMinAng, "\n");
            PHSPACE::WriteLogFile(oss1.str());
        }
    }

#ifdef QUALYTY_DEBUG
    if (nNegativeCell)
    {
        int iZone = this->GetGridID()->GetIndex();
        Grid *grid_slice = CreateGridGeneral(UNSTRUCTGRID, new GridID(iZone), 0, GetDim());
        ExtractSubGridFromGrid(this, UnstructGridCast(grid_slice), negativeCell);
        ostringstream newfile;
        newfile << "./grid/NegativeCells" << iZone << ".plt";
        VisualizationMesh2D(UnstructGridCast(grid_slice), newfile.str());
        delete grid_slice;
    }
    delete [] negativeCell;    negativeCell = nullptr;
#endif

    return globalMinAng;
}

//! Compute the weights used to average dqdx, dqdy, etc.
void UnstructGrid::FaceWeight(RDouble *deltL, RDouble *deltR, int localStart, int localEnd)
{
#ifdef USE_CUDA
    using namespace GPUKernels;
    CallGPUFaceWeight(localStart, localEnd, this->GetNBoundFace());
    return;
#endif
    int nBoundFace = this->GetNBoundFace();
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    RDouble *xfc = this->GetFaceCenterX();
    RDouble *yfc = this->GetFaceCenterY();
    RDouble *zfc = this->GetFaceCenterZ();
    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();

    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int le = leftCellOfFace [ iFace ];
        int re = rightCellOfFace[ iFace ];

        RDouble dxL = xfc[iFace] - xcc[le];
        RDouble dyL = yfc[iFace] - ycc[le];
        RDouble dzL = zfc[iFace] - zcc[le];

        RDouble dxR = xfc[iFace] - xcc[re];
        RDouble dyR = yfc[iFace] - ycc[re];
        RDouble dzR = zfc[iFace] - zcc[re];

        //! Left.
        RDouble delta1 = DISTANCE(dxL, dyL, dzL);
        RDouble delta2 = DISTANCE(dxR, dyR, dzR);
        RDouble delta = one / (delta1 + delta2 + SMALL);

        int jFace = iFace - localStart;
        deltL[jFace] = delta2 * delta;
        deltR[jFace] = delta1 * delta;
    }

    //! If there are some boundary faces, make modifications.
    if (localStart < nBoundFace)
    {
        UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
        int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

        int nMid = nBoundFace;
        if (localEnd < nBoundFace) nMid = localEnd;

        for(int iFace = localStart; iFace < nMid; ++ iFace)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
            int bcType = bcRegion->GetBCType();

            if (!IsInterface(bcType) && bcType != PHENGLEI::OVERSET)
            {
                int jFace = iFace - localStart;
                deltL[jFace] = 1.0;
                deltR[jFace] = 0.0;
            }
        }
    }
}


void UnstructGrid::CVGNorm(FieldProxy *q_1_proxy, FieldProxy *q_0_proxy, int neqn, RDouble &norm)
{
    int nTotalCell = this->GetNTotalCell();

    RDouble diff = 0.0;

    RDouble **q_1 = q_1_proxy->GetField_UNS();
    RDouble **q_0 = q_0_proxy->GetField_UNS();

    norm = zero;
    for (int m = 0; m < neqn; ++ m)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            diff = q_1[m][iCell] - q_0[m][iCell];
            norm = norm + diff * diff;
        }
    }
    norm = sqrt(norm / (nTotalCell * neqn));
}

void UnstructGrid::UpdateVolold()
{
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");

    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    if (!isUnsteady || !isAle) return;
    
    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble *vol  = this->GetCellVolume();
    RDouble *voln = this->GetCellVolumeOld();

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        voln[iCell] = vol[iCell];
    }
}

void UnstructGrid::ClosureCheck(RDouble *xfn, RDouble *area)
{
    RDouble  *sum, total;
    int  nTotalCell = this->GetNTotalCell();
    int  nTotalFace = this->GetNTotalFace();
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    sum = new RDouble[ nTotalCell ];

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        sum[iCell] = 0.0;
    }

    int count, le, re;
    count = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        sum[le] += xfn[iFace] * area[iFace];
        if (re < 0 || re >= nTotalCell) continue;
        sum[re] -= xfn[iFace] * area[iFace];
    }

    total = 0.0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        total += ABS(sum[iCell]);
        if (ABS(sum[iCell]) > 1.0e-8)
        {
        }
    }
    delete [] sum;    sum = nullptr;
}

void UnstructGrid::AllocateMetrics(ActionKey *actkey)
{
    GridManager *gridManager = this->GetGridManager();
    gridManager->AllocateMetrics();

    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotalFace = this->GetNTotalFace();
    int nTotal = nTotalCell + nBoundFace;    //! interior cells + ghost cells.

    RDouble *xfc = this->GetFaceCenterX();
    RDouble *yfc = this->GetFaceCenterY();
    RDouble *zfc = this->GetFaceCenterZ();

    RDouble *xfn = this->GetFaceNormalX();
    RDouble *yfn = this->GetFaceNormalY();
    RDouble *zfn = this->GetFaceNormalZ();

    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();

    RDouble *vol = this->GetCellVolume();

    RDouble *cellSkewness = this->GetCellSkewness();

    RDouble *area = this->GetFaceArea();

    RDouble *vgn = this->GetFaceNormalVelocity();

    if (xcc == 0)
    {
        xcc  = new RDouble[nTotal];
        cellMetrics->SetCellCenterX(xcc);
    }
    if (ycc == 0)
    {
        ycc = new RDouble[nTotal];
        cellMetrics->SetCellCenterY(ycc);
    }
    if (zcc == 0)
    {
        zcc = new RDouble[nTotal];
        cellMetrics->SetCellCenterZ(zcc);
    }
    if (vol == 0)
    {
        vol = new RDouble[nTotal];
        cellMetrics->SetCellVolume(vol);
    }
    if (cellSkewness == 0)
    {
        cellSkewness = new RDouble[nTotal];
        cellMetrics->SetCellSkewness(cellSkewness);
    }

    if (xfc == 0)
    {
        xfc = new RDouble[nTotalFace];
        faceMetrics->SetFaceCenterX(xfc);
    }
    if (yfc == 0)
    {
        yfc = new RDouble[nTotalFace];
        faceMetrics->SetFaceCenterY(yfc);
    }
    if (zfc == 0)
    {
        zfc = new RDouble[nTotalFace];
        faceMetrics->SetFaceCenterZ(zfc);
    }

    if (xfn == 0)
    {
        xfn = new RDouble[nTotalFace];
        faceMetrics->SetFaceNormalX(xfn);
    }
    if (yfn == 0)
    {
        yfn  = new RDouble[nTotalFace];
        faceMetrics->SetFaceNormalY(yfn);
    }
    if (zfn == 0)
    {
        zfn = new RDouble[nTotalFace];
        faceMetrics->SetFaceNormalZ(zfn);
    }
    if (area == 0)
    {
        area = new RDouble[nTotalFace];
        faceMetrics->SetFaceArea(area);
    }

    if (vgn  == 0)
    {
        vgn = new RDouble[nTotalFace];
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            vgn[iFace] = zero;
        }
        dynamicGridMetrics->SetFaceNormalVelocity(vgn);
    }

    AllocateMetricsALE(actkey);
}

void UnstructGrid::AllocateMetricsALE(ActionKey *actkey)
{
    GridManager *gridManager = this->GetGridManager();
    gridManager->AllocateALEMetrics();

    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotalFace = this->GetNTotalFace();
    int nTotal = nTotalCell + nBoundFace;    //! interior cells + ghost cells.

    //! codeOfAleModel = 0: no ALE method;
    //!                  1: ALE method for non-moving grids;
    //!                  2: ALE method for moving grids;
    //!                  3: ALE method for deforming grids.

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    RDouble *xfv = this->GetFaceVelocityX();
    RDouble *yfv = this->GetFaceVelocityY();
    RDouble *zfv = this->GetFaceVelocityZ();

    RDouble *xcv = this->GetCellVelocityX();
    RDouble *ycv = this->GetCellVelocityY();
    RDouble *zcv = this->GetCellVelocityZ();

    int ifLowSpeedPrecon =  GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");
    if (ifLowSpeedPrecon != 0)
    {
        if (xcv == 0)
        {
            xcv = new RDouble[nTotalCell];
        }

        if (ycv == 0)
        {
            ycv = new RDouble[nTotalCell];
        }

        if (zcv == 0)
        {
            zcv = new RDouble[nTotalCell];
        }

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            xcv[iCell] = 0.0;
            ycv[iCell] = 0.0;
            zcv[iCell] = 0.0;
        }

        dynamicGridMetrics->SetCellVelocityX(xcv);
        dynamicGridMetrics->SetCellVelocityY(ycv);
        dynamicGridMetrics->SetCellVelocityZ(zcv);
    }

    if (xfv == 0)
    {
        xfv = new RDouble[nTotalFace];
    }

    if (yfv == 0)
    {
        yfv = new RDouble[nTotalFace];
    }

    if (zfv == 0)
    {
        zfv = new RDouble[nTotalFace];
    }

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        xfv[iFace] = 0.0;
        yfv[iFace] = 0.0;
        zfv[iFace] = 0.0;
    }
    dynamicGridMetrics->SetFaceVelocityX(xfv);
    dynamicGridMetrics->SetFaceVelocityY(yfv);
    dynamicGridMetrics->SetFaceVelocityZ(zfv);

    RDouble *voln = this->GetCellVolumeOld();

    if (isUnsteady)
    {
        if (!isAle)
        {
            voln = 0;
        }
        else
        {
            if (voln == 0)
            {
                voln = new RDouble[nTotal];
                ZeroField(voln, nTotal);
            }

            int nstart = 0;
            if (!(nstart == 0 || nstart == 3))
            {
                //! nstart==0: start all over again.
                //! nstart==3: Unsteady-->Dynamics mesh.
                ReadMovingGrid();
            }
        }
        dynamicGridMetrics->SetCellVolumeOld(voln);
    }
}

void UnstructGrid::GetBcNamePair(set< pair<int, string> > &bcNamePairSet)
{
    UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        string bcName = bcRegion->GetBCName();
        int bcType = bcRegion->GetBCType();

        pair<int, string> bcNamePair(bcType, bcName);
        bcNamePairSet.insert(bcNamePair);
    }
}

void UnstructGrid::DumpWallFaceCenter(ActionKey *actkey)
{
    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();

    std::ostringstream oss;

    int girdIndex = this->GetZoneID();
    if (girdIndex == 0)
    {
        int globalTotalWallFace = GlobalDataBase::GetIntParaFromDB("GlobalTotalWallFace");
        oss << globalTotalWallFace << "\n";
    }

    int nsolid_surface = this->GetNumberOfWallCell();
    if ((nsolid_surface == 0) || (level != 0))
    {
        string str = oss.str();
        cdata->Write(const_cast <char *> (str.data()), str.size() * sizeof(char));
        return;
    }

    oss << setiosflags(ios::left);
    oss << setprecision(5);
    oss << setiosflags(ios::scientific);
    oss << setiosflags(ios::showpoint);

    RDouble *xfc = this->GetFaceCenterX();
    RDouble *yfc = this->GetFaceCenterY();
    RDouble *zfc = this->GetFaceCenterZ();

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bctype = bcr[iFace]->GetKey();
        if (!IsWall(bctype))
        {
            continue;
        }

        oss << xfc[iFace] << "	";
        oss << yfc[iFace] << "	";
        oss << zfc[iFace] << "	";
        oss <<"\n";
    }

    string str = oss.str();
    cdata->Write(const_cast <char *> (str.data()), str.size() * sizeof(char));
}

void UnstructGrid::InitVariableWallTemperature()
{
    fstream file;
    string wallTemperatureFile = "";
    GlobalDataBase::GetData("wallTemperaturefile", &wallTemperatureFile, PHSTRING, 1);
    OpenFile(file, wallTemperatureFile, ios_base::in|ios_base::binary);

    int dataNumber = 0;
    PHRead(file, dataNumber);

    int varNumber = 0;
    PHRead(file, varNumber);

    RDouble *pmin = this->GetMinBox();
    RDouble *pmax = this->GetMaxBox();
    vector < vector <RDouble> > dataInCurrentGrid;
    dataInCurrentGrid.resize(0);

    for (int iData = 0; iData < dataNumber; ++ iData)
    {
        float xCoor, yCoor, zCoor;
        float tWall;

        PHRead(file, xCoor);
        PHRead(file, yCoor);
        PHRead(file, zCoor);
        PHRead(file, tWall);

        RDouble coordX = static_cast<RDouble>(xCoor);
        RDouble coordY = static_cast<RDouble>(yCoor);
        RDouble coordZ = static_cast<RDouble>(zCoor);

        if (coordX < pmin[0] || coordX > pmax[0]) continue;
        if (coordY < pmin[1] || coordY > pmax[1]) continue;
        if (coordZ < pmin[2] || coordZ > pmax[2]) continue;

        vector <RDouble> currentData;
        currentData.push_back(coordX);
        currentData.push_back(coordY);
        currentData.push_back(coordZ);
        currentData.push_back(static_cast<RDouble>(tWall));

        dataInCurrentGrid.push_back(currentData);
    }

    CloseFile(file);

    int dataNumberInCurrentGrid = static_cast<int>(dataInCurrentGrid.size());
    if (dataNumberInCurrentGrid == 0)
    {
        return;
    }

    int *dataUsed = new int[dataNumberInCurrentGrid];
    std::fill_n(dataUsed, dataNumberInCurrentGrid, 0);

    RDouble *xCoor    = this->GetX();
    RDouble *yCoor    = this->GetY();
    RDouble *zCoor    = this->GetZ();
    RDouble *faceArea = this->GetFaceArea();

    long long int *nodePosi  = this->GetFace2NodeSubscript();
    int *face2node = this->GetFace2Node();

    UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        
        if (!IsWall(bcType) && bcType != PHENGLEI::ABLATION_SURFACE)
        {
            continue;
        }

        vector<int> *faceIndexArray = bcRegion->GetFaceIndex();
        int faceNum = static_cast<int>(faceIndexArray->size());
        RDouble *wallTempArray = new RDouble[faceNum];

        int findData = 0;
        for (int iFace = 0; iFace < faceNum; ++ iFace)
        {
            bool flag = false;
            vector <vector <RDouble> > nodeCoor;

            int faceIndex = (*faceIndexArray)[iFace];
            for (int iNode = static_cast<int>(nodePosi[faceIndex]); iNode < nodePosi[faceIndex+1]; ++ iNode)
            {
                int nodeIndex = face2node[iNode];

                vector <RDouble> oneNode;
                oneNode.push_back(xCoor[nodeIndex]);
                oneNode.push_back(yCoor[nodeIndex]);
                oneNode.push_back(zCoor[nodeIndex]);
                nodeCoor.push_back(oneNode);
            }
            int nTri = static_cast<int>(nodeCoor.size());

            for (int iData = 0; iData < dataNumberInCurrentGrid; ++ iData)
            {
                if (dataUsed[iData])
                {
                    continue;
                }

                RDouble coordX = dataInCurrentGrid[iData][0];
                RDouble coordY = dataInCurrentGrid[iData][1];
                RDouble coordZ = dataInCurrentGrid[iData][2];

                vector <vector <RDouble> > lineVec;
                lineVec.resize(nTri);

                for (int iNode = 0; iNode < nTri; ++ iNode)
                {
                    RDouble dx = nodeCoor[iNode][0] - coordX;
                    RDouble dy = nodeCoor[iNode][1] - coordY;
                    RDouble dz = nodeCoor[iNode][2] - coordZ;

                    lineVec[iNode].push_back(dx);
                    lineVec[iNode].push_back(dy);
                    lineVec[iNode].push_back(dz);
                }

                RDouble areaSum = 0.0;
                for (int iTri = 0; iTri < nTri; ++ iTri)
                {
                    int lineIndex1 =  iTri;
                    int lineIndex2 = (iTri + 1) % nTri;

                    RDouble p[3];
                    CrossProduct(&lineVec[lineIndex1][0], &lineVec[lineIndex2][0], p);
                    RDouble areaCurrent = DISTANCE(p[0], p[1], p[2]);

                    areaSum += ABS(areaCurrent);
                }
                areaSum = half * areaSum;

                RDouble eps = (areaSum - faceArea[faceIndex]) / faceArea[faceIndex];
                if (ABS(eps) < 1e-4)
                {
                    flag = true;
                }

                if (flag)
                {
                    dataUsed[iData] = 1;
                    RDouble tWall = dataInCurrentGrid[iData][3];
                    wallTempArray[iFace] = tWall;
                    findData ++;
                    break;
                }
            }
        }

        if (faceNum != findData)
        {
            WriteLogFile("Error !!!!");
        }

       bcRegion->UpdateFieldDataPtr("wallTempArray" , wallTempArray);

    }

    delete [] dataUsed;    dataUsed = nullptr;
}

void UnstructGrid::ReadMovingGrid()
{
    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;    //! interior cells + ghost cells.

    fstream file;
    std::ostringstream oss;
    oss << "move_grid" << this->GetZoneID() << ".dat";
    string filename = oss.str();

    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();
    RDouble *voln = this->GetCellVolumeOld();

    int nTotalNode = this->GetNTotalNode();

    file.open(filename.c_str(), ios_base::in|ios_base::binary);

    file.read(reinterpret_cast<char*>(x  ),nTotalNode * sizeof(RDouble));
    file.read(reinterpret_cast<char*>(y  ),nTotalNode * sizeof(RDouble));
    file.read(reinterpret_cast<char*>(z  ),nTotalNode * sizeof(RDouble));
    file.read(reinterpret_cast<char*>(voln),nTotal * sizeof(RDouble));

    file.close();
}

//! Corrected GG Method: on the basis of GG-Cell, the face center value is corrected once,so that the gradient reconstruction can reach the first order accuracy.
void UnstructGrid::CompGradientGGModified2(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    RDouble *nxs,*nys,*nzs,*ns,*vol;

    nxs = this->GetFaceNormalX();
    nys = this->GetFaceNormalY();
    nzs = this->GetFaceNormalZ();
    ns  = this->GetFaceArea() ;
    vol = this->GetCellVolume();

    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    RDouble *xfc = this->GetFaceCenterX();
    RDouble *yfc = this->GetFaceCenterY();
    RDouble *zfc = this->GetFaceCenterZ();
    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();

    int nTotalCell = this->GetNTotalCell();
    int nTotalFace = this->GetNTotalFace();
    int nBoundFace = this->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        dqdx[iCell] = 0.0;
        dqdy[iCell] = 0.0;
        dqdz[iCell] = 0.0;
    }

    //! Initial value of gradient.
    vector< RDouble > dqdx0(nTotal);
    vector< RDouble > dqdy0(nTotal);
    vector< RDouble > dqdz0(nTotal);

    //! Solve the initial gradient field.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[ iFace ];
        int re = iFace + nTotalCell;

        RDouble qfc = half * (q[le] + q[re]);
        RDouble nx  = nxs[iFace] * ns[iFace];
        RDouble ny  = nys[iFace] * ns[iFace];
        RDouble nz  = nzs[iFace] * ns[iFace];
        dqdx0[le] += qfc * nx;
        dqdy0[le] += qfc * ny;
        dqdz0[le] += qfc * nz;
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [ iFace ];
        int re = rightCellOfFace[ iFace ];

        RDouble qfc = half * (q[le] + q[re]);
        RDouble nx  = nxs[iFace] * ns[iFace];
        RDouble ny  = nys[iFace] * ns[iFace];
        RDouble nz  = nzs[iFace] * ns[iFace];

        dqdx0[le] += qfc * nx;
        dqdy0[le] += qfc * ny;
        dqdz0[le] += qfc * nz;
        dqdx0[re] -= qfc * nx;
        dqdy0[re] -= qfc * ny;
        dqdz0[re] -= qfc * nz;
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble ovol = one / vol[iCell];
        dqdx0[iCell] *= ovol;
        dqdy0[iCell] *= ovol;
        dqdz0[iCell] *= ovol;
    }    

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[ iFace ];
        int re = iFace + nTotalCell;

        dqdx0[re] = dqdx0[le];
        dqdy0[re] = dqdy0[le];
        dqdz0[re] = dqdz0[le];
    } 

    //! After the initial gradient field is solved,correct the face center value and update the gradient.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[ iFace ];
        int re = iFace + nTotalCell;

        RDouble faceCenterX = xfc[ iFace ];
        RDouble faceCenterY = yfc[ iFace ];
        RDouble faceCenterZ = zfc[ iFace ];

        //! The vector of ff.
        vector< RDouble > fcFVector(3);
        fcFVector[ 0 ] = faceCenterX - half * (xcc[ le ] + xcc[ re ]);
        fcFVector[ 1 ] = faceCenterY - half * (ycc[ le ] + ycc[ re ]);
        fcFVector[ 2 ] = faceCenterZ - half * (zcc[ le ] + zcc[ re ]);

        //! The gradient of f.
        RDouble dqdxf = half * (dqdx0[ le ] + dqdx0[ re ]);
        RDouble dqdyf = half * (dqdy0[ le ] + dqdy0[ re ]);
        RDouble dqdzf = half * (dqdz0[ le ] + dqdz0[ re ]);

        //! Correct the face center value.
        RDouble qfc0 = half * (q[le] + q[re]);
        RDouble qfc = qfc0 + dqdxf * fcFVector[ 0 ] + dqdyf * fcFVector[ 1 ] + dqdzf * fcFVector[ 2 ];

        RDouble nx  = nxs[iFace] * ns[iFace];
        RDouble ny  = nys[iFace] * ns[iFace];
        RDouble nz  = nzs[iFace] * ns[iFace];

        dqdx[le] += qfc * nx;
        dqdy[le] += qfc * ny;
        dqdz[le] += qfc * nz;
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [ iFace ];
        int re = rightCellOfFace[ iFace ];

        RDouble faceCenterX = xfc[ iFace ];
        RDouble faceCenterY = yfc[ iFace ];
        RDouble faceCenterZ = zfc[ iFace ];

        //! The vector of ff.
        vector< RDouble > fcFVector(3);
        fcFVector[ 0 ] = faceCenterX - half * (xcc[ le ] + xcc[ re ]);
        fcFVector[ 1 ] = faceCenterY - half * (ycc[ le ] + ycc[ re ]);
        fcFVector[ 2 ] = faceCenterZ - half * (zcc[ le ] + zcc[ re ]);

        //! The gradient of f.
        RDouble dqdxf = half * (dqdx0[ le ] + dqdx0[ re ]);
        RDouble dqdyf = half * (dqdy0[ le ] + dqdy0[ re ]);
        RDouble dqdzf = half * (dqdz0[ le ] + dqdz0[ re ]);

        ///! Correct the face center value.
        RDouble qfc0 = half * (q[le] + q[re]);
        RDouble qfc = qfc0 + dqdxf * fcFVector[ 0 ] + dqdyf * fcFVector[ 1 ] + dqdzf * fcFVector[ 2 ];

        RDouble nx  = nxs[iFace] * ns[iFace];
        RDouble ny  = nys[iFace] * ns[iFace];
        RDouble nz  = nzs[iFace] * ns[iFace];

        dqdx[le] += qfc * nx;
        dqdy[le] += qfc * ny;
        dqdz[le] += qfc * nz;
        dqdx[re] -= qfc * nx;
        dqdy[re] -= qfc * ny;
        dqdz[re] -= qfc * nz;
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {       
        RDouble ovol = one / vol[iCell];
        dqdx[iCell] *= ovol;
        dqdy[iCell] *= ovol;
        dqdz[iCell] *= ovol;
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[ iFace ];
        int re = iFace + nTotalCell;

        dqdx[re] = dqdx[le];
        dqdy[re] = dqdy[le];
        dqdz[re] = dqdz[le];
    }
}

void UnstructGrid::CompGradientGGCell(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, int setting)
{
    RDouble *nxs,*nys,*nzs,*ns,*vol;

    nxs = this->GetFaceNormalX();
    nys = this->GetFaceNormalY();
    nzs = this->GetFaceNormalZ();
    ns  = this->GetFaceArea() ;
    vol = this->GetCellVolume();

    int *leftCellofFace = this->GetLeftCellOfFace();
    int *rightCellofFace = this->GetRightCellOfFace();

    int nTotalCell = this->GetNTotalCell();
    int nTotalFace = this->GetNTotalFace();
    int nBoundFace = this->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        dqdx[iCell] = 0.0;
        dqdy[iCell] = 0.0;
        dqdz[iCell] = 0.0;
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellofFace[ iFace ];
        int re = iFace + nTotalCell;

        RDouble qfc = half * (q[le] + q[re]);
        RDouble nx  = nxs[iFace] * ns[iFace];
        RDouble ny  = nys[iFace] * ns[iFace];
        RDouble nz  = nzs[iFace] * ns[iFace];
        dqdx[le] += qfc * nx;
        dqdy[le] += qfc * ny;
        dqdz[le] += qfc * nz;
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellofFace [ iFace ];
        int re = rightCellofFace[ iFace ];

        RDouble qfc = half * (q[le] + q[re]);
        RDouble nx  = nxs[iFace] * ns[iFace];
        RDouble ny  = nys[iFace] * ns[iFace];
        RDouble nz  = nzs[iFace] * ns[iFace];

        dqdx[le] += qfc * nx;
        dqdy[le] += qfc * ny;
        dqdz[le] += qfc * nz;
        dqdx[re] -= qfc * nx;
        dqdy[re] -= qfc * ny;
        dqdz[re] -= qfc * nz;
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble ovol = one / vol[iCell];
        dqdx[iCell] *= ovol;
        dqdy[iCell] *= ovol;
        dqdz[iCell] *= ovol;
    }
}

void UnstructGrid::CompGradientGGCell2(RDouble *q, RDouble *wfL, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    RDouble *nxs, *nys, *nzs, *ns, *vol;

    nxs = this->GetFaceNormalX();
    nys = this->GetFaceNormalY();
    nzs = this->GetFaceNormalZ();
    ns = this->GetFaceArea();
    vol = this->GetCellVolume();

    int *leftCellofFace = this->GetLeftCellOfFace();
    int *rightCellofFace = this->GetRightCellOfFace();

    int nTotalCell = this->GetNTotalCell();
    int nTotalFace = this->GetNTotalFace();
    int nBoundFace = this->GetNBoundFace();

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dqdx[iCell] = 0.0;
        dqdy[iCell] = 0.0;
        dqdz[iCell] = 0.0;
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        int re = iFace + nTotalCell;

        RDouble qfc = q[re]; 
        int bcType = bcr[iFace]->GetKey();
        if (bcType == PHENGLEI::INTERFACE)
        {
            qfc = wfL[iFace] * q[le] + (1.0 - wfL[iFace]) * q[re];
        }
        RDouble nx = nxs[iFace] * ns[iFace];
        RDouble ny = nys[iFace] * ns[iFace];
        RDouble nz = nzs[iFace] * ns[iFace];
        dqdx[le] += qfc * nx;
        dqdy[le] += qfc * ny;
        dqdz[le] += qfc * nz;
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

        RDouble nx = nxs[iFace] * ns[iFace];
        RDouble ny = nys[iFace] * ns[iFace];
        RDouble nz = nzs[iFace] * ns[iFace];

        RDouble qfc = wfL[iFace] * q[le] + (1.0 - wfL[iFace]) * q[re];
        dqdx[le] += qfc * nx;
        dqdy[le] += qfc * ny;
        dqdz[le] += qfc * nz;
        dqdx[re] -= qfc * nx;
        dqdy[re] -= qfc * ny;
        dqdz[re] -= qfc * nz;
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble ovol = one / vol[iCell];
        dqdx[iCell] *= ovol;
        dqdy[iCell] *= ovol;
        dqdz[iCell] *= ovol;
    }

    UnstructBCSet **bcr = this->GetBCRecord();
    for (int iFace = 0; iFace < nBoundFace; iFace++)
    {
        int bcType = bcr[iFace]->GetKey();
        if (bcType == PHENGLEI::INTERFACE)
        {
            continue;
        }
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];
        dqdx[re] = dqdx[le];
        dqdy[re] = dqdy[le];
        dqdz[re] = dqdz[le];
    }
}

void UnstructGrid::CompGradientGGCellNew(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    RDouble *nxs,*nys,*nzs,*ns,*vol;

    nxs = this->GetFaceNormalX();
    nys = this->GetFaceNormalY();
    nzs = this->GetFaceNormalZ();
    ns  = this->GetFaceArea() ;
    vol = this->GetCellVolume();

    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    int nTotalCell = this->GetNTotalCell();
    int nTotalFace = this->GetNTotalFace();
    int nBoundFace = this->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        dqdx[iCell] = 0.0;
        dqdy[iCell] = 0.0;
        dqdz[iCell] = 0.0;
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[ iFace ];
        int re = iFace + nTotalCell;

        RDouble qfc = half * (q[le] + q[re]);
        RDouble nx  = nxs[iFace] * ns[iFace];
        RDouble ny  = nys[iFace] * ns[iFace];
        RDouble nz  = nzs[iFace] * ns[iFace];
        dqdx[le] += qfc * nx;
        dqdy[le] += qfc * ny;
        dqdz[le] += qfc * nz;
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [ iFace ];
        int re = rightCellOfFace[ iFace ];

        RDouble qfc = half * (q[le] + q[re]);

        RDouble qfl = qfc;
        RDouble qfr = qfc;

        if (ABS(q[re] - q[le]) > 1.0e-1)
        {
            qfl = q[le];
            qfr = q[re];
        }

        RDouble nx  = nxs[iFace] * ns[iFace];
        RDouble ny  = nys[iFace] * ns[iFace];
        RDouble nz  = nzs[iFace] * ns[iFace];

        dqdx[le] += qfl * nx;
        dqdy[le] += qfl * ny;
        dqdz[le] += qfl * nz;
        dqdx[re] -= qfr * nx;
        dqdy[re] -= qfr * ny;
        dqdz[re] -= qfr * nz;
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble ovol = one / vol[iCell];
        dqdx[iCell] *= ovol;
        dqdy[iCell] *= ovol;
        dqdz[iCell] *= ovol;
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[ iFace ];
        int re = iFace + nTotalCell;

        dqdx[re] = dqdx[le];
        dqdy[re] = dqdy[le];
        dqdz[re] = dqdz[le];
    }

}

void UnstructGrid::CompGradientGGCellW(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    int nTotalCell = this->GetNTotalCell();
    int nTotalFace = this->GetNTotalFace();
    int nBoundFace = this->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble *xfn = this->GetFaceNormalX();
    RDouble *yfn = this->GetFaceNormalY();
    RDouble *zfn = this->GetFaceNormalZ();
    RDouble *area= this->GetFaceArea();
    RDouble *vol = this->GetCellVolume();

    RDouble *xfc = this->GetFaceCenterX();
    RDouble *yfc = this->GetFaceCenterY();
    RDouble *zfc = this->GetFaceCenterZ();
    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();

    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        dqdx[iCell] = 0.0;
        dqdy[iCell] = 0.0;
        dqdz[iCell] = 0.0;
    }

    RDouble dxl,dyl,dzl,dxr,dyr,dzr;
    RDouble delt1,delt2,delta,cl,cr;
    RDouble fc;
    RDouble ovol;

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [ iFace ];
        int re = rightCellOfFace[ iFace ];

        dxl = xfc[iFace] - xcc[le];
        dyl = yfc[iFace] - ycc[le];
        dzl = zfc[iFace] - zcc[le];

        dxr = xfc[iFace] - xcc[re];
        dyr = yfc[iFace] - ycc[re];
        dzr = zfc[iFace] - zcc[re];

        delt1  = sqrt(dxl * dxl + dyl * dyl + dzl * dzl);
        delt2  = sqrt(dxr * dxr + dyr * dyr + dzr * dzr);
        delta  = one / (delt1 + delt2 + SMALL);
        cl    = delt2 * delta;
        cr    = delt1 * delta;

        fc = cl * q[le] + cr * q[re];
        RDouble nx  = xfn[iFace] * area[iFace];
        RDouble ny  = yfn[iFace] * area[iFace];
        RDouble nz  = zfn[iFace] * area[iFace];

        dqdx[le] += nx * fc;
        dqdy[le] += ny * fc;
        dqdz[le] += nz * fc;

        if (iFace >= nBoundFace)
        {
            dqdx[re] -= nx * fc;
            dqdy[re] -= ny * fc;
            dqdz[re] -= nz * fc;
        }
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        ovol = one / vol[iCell];
        dqdx[iCell] *= ovol;
        dqdy[iCell] *= ovol;
        dqdz[iCell] *= ovol;
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[ iFace ];
        int re = iFace + nTotalCell;

        dqdx[re] = dqdx[le];
        dqdy[re] = dqdy[le];
        dqdz[re] = dqdz[le];
    }
}

void UnstructGrid::CompFaceGradient(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    int nTotalNode = this->GetNTotalNode();
    RDouble *qn = new RDouble [nTotalNode];

    CompNodeVar(this, qn, q);

    if (this->GetDim() == TWO_D)
    {
        CompFaceGradientGauss2D(q, qn, dqdx, dqdy, dqdz);
    }
    else
    {
        CompFaceGradientGauss(q, qn, dqdx, dqdy, dqdz);
    }

    delete [] qn;    qn = nullptr;
}

void UnstructGrid::CompFaceGradientGauss2D(RDouble *q, RDouble *qn, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    int nTotalFace = this->GetNTotalFace();

    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    int *face2node = this->GetFace2Node();
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();

    RDouble *xfn  = this->GetFaceNormalX();
    RDouble *yfn  = this->GetFaceNormalY();
    RDouble *zfn  = this->GetFaceNormalZ();
    RDouble *area = this->GetFaceArea();

    RDouble *xcc  = this->GetCellCenterX();
    RDouble *ycc  = this->GetCellCenterY();
    RDouble *zcc  = this->GetCellCenterZ();

    RDouble dx, dy, dz, dh, ovol;
    RDouble anx,any,anz;
    RDouble qfc;
    int le,re,p1,p2;
    int nodepos = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [ iFace ];
        re = rightCellOfFace[ iFace ];

        dx = xcc[re] - xcc[le];
        dy = ycc[re] - ycc[le];
        dz = zcc[re] - zcc[le];

        dh = dx * xfn[iFace] + dy * yfn[iFace] + dz * zfn[iFace];
        ovol = one / (half * area[iFace] * dh);

        dqdx[iFace] = 0.0;
        dqdy[iFace] = 0.0;
        dqdz[iFace] = 0.0;

        p1 = face2node[ nodepos + 0 ];
        p2 = face2node[ nodepos + 1 ];

        qfc = half * (qn[p1] + q[re]);

        anx = ycc[re] - y[p1];
        any = x[p1] - xcc[re];
        anz = 0.0;

        dqdx[iFace] += qfc * anx;
        dqdy[iFace] += qfc * any;
        dqdz[iFace] += qfc * anz;

        qfc = half * (qn[p2] + q[re]);

        anx = y[p2] - ycc[re];
        any = xcc[re] - x[p2];
        anz = 0.0;

        dqdx[iFace] += qfc * anx;
        dqdy[iFace] += qfc * any;
        dqdz[iFace] += qfc * anz;

        qfc = half * (qn[p1] + q[le]);

        anx = ycc[le] - y[p1];
        any = x[p1] - xcc[le];
        anz = 0.0;

        dqdx[iFace] -= qfc * anx;
        dqdy[iFace] -= qfc * any;
        dqdz[iFace] -= qfc * anz;

        qfc = half * (qn[p2] + q[le]);

        anx = y[p2] - ycc[le];
        any = xcc[le] - x[p2];
        anz = 0.0;

        dqdx[iFace] -= qfc * anx;
        dqdy[iFace] -= qfc * any;
        dqdz[iFace] -= qfc * anz;

        dqdx[iFace] *= ovol;
        dqdy[iFace] *= ovol;
        dqdz[iFace] *= ovol;
        nodepos += node_number_of_each_face[iFace];
    }
}

void UnstructGrid::CompFaceGradientGauss(RDouble *q, RDouble *qn, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    int nTotalFace = this->GetNTotalFace();

    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();
    int *face2node = this->GetFace2Node();
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();

    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();

    RDouble vvol = 0.0;
    RDouble dx1 = 0.0, dy1 = 0.0, dz1 = 0.0, dx2 = 0.0, dy2 = 0.0, dz2 = 0.0;
    RDouble xc = 0.0, yc = 0.0, zc = 0.0, anx = 0.0, any = 0.0, anz = 0.0;
    RDouble fc = 0.0, tmp = 0.0, coef = 0.0;
    int le = 0, re = 0, p1 = 0, p2 = 0, id1 = 0, id2 = 0;
    int nodepos = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [ iFace ];
        re = rightCellOfFace[ iFace ];

        vvol        = 0.0;
        dqdx[iFace] = 0.0;
        dqdy[iFace] = 0.0;
        dqdz[iFace] = 0.0;

        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            id1 = j;
            id2 = (j + 1) % node_number_of_each_face[iFace];
            p1 = face2node[ nodepos + id1 ];
            p2 = face2node[ nodepos + id2 ];

            dx1 = x[p2] - x[p1];
            dy1 = y[p2] - y[p1];
            dz1 = z[p2] - z[p1];

            dx2 = x[p1] - xcc[le];
            dy2 = y[p1] - ycc[le];
            dz2 = z[p1] - zcc[le];

            anx = dy1 * dz2 - dy2 * dz1;
            any = dz1 * dx2 - dz2 * dx1;
            anz = dx1 * dy2 - dx2 * dy1;

            fc = qn[p1] + qn[p2] + q[le];

            //! Face centers.
            xc = xcc[le] + x[p1] + x[p2];
            yc = ycc[le] + y[p1] + y[p2];
            zc = zcc[le] + z[p1] + z[p2];

            tmp   = anx * xc + any * yc + anz * zc;
            vvol -= tmp;

            dqdx[iFace] += anx * fc;
            dqdy[iFace] += any * fc;
            dqdz[iFace] += anz * fc;


            dx1 = x[p2] - x[p1];
            dy1 = y[p2] - y[p1];
            dz1 = z[p2] - z[p1];

            dx2 = xcc[re] - x[p1];
            dy2 = ycc[re] - y[p1];
            dz2 = zcc[re] - z[p1];

            fc = qn[p1] + qn[p2] + q[le];

            //! Face centers.
            xc = xcc[re] + x[p1] + x[p2];
            yc = ycc[re] + y[p1] + y[p2];
            zc = zcc[re] + z[p1] + z[p2];

            anx = dy1 * dz2 - dy2 * dz1;
            any = dz1 * dx2 - dz2 * dx1;
            anz = dx1 * dy2 - dx2 * dy1;

            tmp   = anx * xc + any * yc + anz * zc;
            vvol -= tmp;

            dqdx[iFace] += anx * fc;
            dqdy[iFace] += any * fc;
            dqdz[iFace] += anz * fc;
        }

        coef = three / vvol;
        dqdx[iFace] *= coef;
        dqdy[iFace] *= coef;
        dqdz[iFace] *= coef;
        nodepos += node_number_of_each_face[iFace];
    }
}

void UnstructGrid::ComputeWeight()
{
    string gradientName = "lsq";
    GetData("gradientName", &gradientName, PHSTRING, 1);
    if (gradientName != "lsq") return;
    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotalFace = this->GetNTotalFace();

    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();

    if (!this->GetLeastSquareIWT()) this->SetLeastSquareIWT(new RDouble [ nTotalFace ]);
    if (!this->GetLeastSquareIXX()) this->SetLeastSquareIXX(new RDouble [ nTotalCell ]);
    if (!this->GetLeastSquareIYY()) this->SetLeastSquareIYY(new RDouble [ nTotalCell ]);
    if (!this->GetLeastSquareIZZ()) this->SetLeastSquareIZZ(new RDouble [ nTotalCell ]);
    if (!this->GetLeastSquareIXY()) this->SetLeastSquareIXY(new RDouble [ nTotalCell ]);
    if (!this->GetLeastSquareIXZ()) this->SetLeastSquareIXZ(new RDouble [ nTotalCell ]);
    if (!this->GetLeastSquareIYZ()) this->SetLeastSquareIYZ(new RDouble [ nTotalCell ]);
    if (!this->GetFaceMark()) this->SetFaceMark(new char [ nTotalFace ]);

    RDouble *iwt   = this->GetLeastSquareIWT();
    RDouble *ixx   = this->GetLeastSquareIXX();
    RDouble *iyy   = this->GetLeastSquareIYY();
    RDouble *izz   = this->GetLeastSquareIZZ();
    RDouble *ixy   = this->GetLeastSquareIXY();
    RDouble *ixz   = this->GetLeastSquareIXZ();
    RDouble *iyz   = this->GetLeastSquareIYZ();
    char    *fMark = this->GetFaceMark();

    int *leftCellOfFace  = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();
    
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        fMark[iFace] = 1;
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        ixx[iCell] = 0.0;
        iyy[iCell] = 0.0;
        izz[iCell] = 0.0;
        ixy[iCell] = 0.0;
        ixz[iCell] = 0.0;
        iyz[iCell] = 0.0;
    }

    int le, re;

    RDouble dx,dy,dz,txx,tyy,tzz,txy,txz,tyz;
    RDouble ods;
  
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace[ iFace ];
        re = iFace + nTotalCell;

        dx = (xcc[re] - xcc[le]) * fMark[iFace];
        dy = (ycc[re] - ycc[le]) * fMark[iFace];
        dz = (zcc[re] - zcc[le]) * fMark[iFace];

        ods = one / (dx * dx + dy * dy + dz * dz + SMALL);

        txx = ods * dx * dx;
        tyy = ods * dy * dy;
        tzz = ods * dz * dz;
        txy = ods * dx * dy;
        txz = ods * dx * dz;
        tyz = ods * dy * dz;

        iwt[iFace]  = ods;
        ixx[le] += txx;
        iyy[le] += tyy;
        izz[le] += tzz;
        ixy[le] += txy;
        ixz[le] += txz;
        iyz[le] += tyz;
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le       = leftCellOfFace [ iFace ];
        re       = rightCellOfFace[ iFace ];

        dx       = (xcc[re] - xcc[le]) * fMark[iFace];
        dy       = (ycc[re] - ycc[le]) * fMark[iFace];
        dz       = (zcc[re] - zcc[le]) * fMark[iFace];

        ods = one / (dx * dx + dy * dy + dz * dz + SMALL);

        txx      = ods * dx * dx;
        tyy      = ods * dy * dy;
        tzz      = ods * dz * dz;
        txy      = ods * dx * dy;
        txz      = ods * dx * dz;
        tyz      = ods * dy * dz;

        iwt[iFace]  = ods;

        ixx[le] += txx;
        iyy[le] += tyy;
        izz[le] += tzz;
        ixy[le] += txy;
        ixz[le] += txz;
        iyz[le] += tyz;

        ixx[re] += txx;
        iyy[re] += tyy;
        izz[re] += tzz;
        ixy[re] += txy;
        ixz[re] += txz;
        iyz[re] += tyz;
    }

    RDouble delt = 0.0, odelt = 0.0;

    if (this->GetDim() == THREE_D)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            txx = ixx[iCell];
            tyy = iyy[iCell];
            tzz = izz[iCell];
            txy = ixy[iCell];
            txz = ixz[iCell];
            tyz = iyz[iCell];
          
            delt = (txx * tyy * tzz
                  + 2.0 * txy * txz * tyz 
                  - txx * tyz * tyz
                  - tyy * txz * txz
                  - tzz * txy * txy);

            if (ABS(delt) < SMALL)
            {
                ixx[iCell] = 0.0;
                iyy[iCell] = 0.0;
                izz[iCell] = 0.0;
                ixy[iCell] = 0.0;
                ixz[iCell] = 0.0;
                iyz[iCell] = 0.0;
            }
            else
            {
                odelt = one / delt;
                ixx[iCell] = (tyy * tzz - tyz * tyz) * odelt;
                iyy[iCell] = (txx * tzz - txz * txz) * odelt;
                izz[iCell] = (txx * tyy - txy * txy) * odelt;
                ixy[iCell] = (txz * tyz - txy * tzz) * odelt;
                ixz[iCell] = (txy * tyz - txz * tyy) * odelt;
                iyz[iCell] = (txy * txz - tyz * txx) * odelt;
           }
        }
    }
    else
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            txx = ixx[iCell];
            tyy = iyy[iCell];
            txy = ixy[iCell];

            delt = txx * tyy - txy * txy;

            if (ABS(delt) < SMALL)
            {
                ixx[iCell] = 0.0;
                iyy[iCell] = 0.0;
                izz[iCell] = 0.0;
                ixy[iCell] = 0.0;
                ixz[iCell] = 0.0;
                iyz[iCell] = 0.0;
            }
            else
            {
                odelt = one / delt;

                ixx[iCell] = (  tyy) * odelt;
                iyy[iCell] = (  txx) * odelt;
                ixy[iCell] = (- txy) * odelt;
                izz[iCell] = 0.0;
                ixz[iCell] = 0.0;
                iyz[iCell] = 0.0;
            }
        }
    }
}

LIB_EXPORT void UnstructGrid::ComputeNodeBCType()
{
    long long int *face2NodeSubscript = this->GetFace2NodeSubscript();
    int *face2node = this->GetFace2Node();

    int *nodeBCType = reinterpret_cast <int *> (this->GetDataPtr("nodeBCType"));
    int *nodeBCTypeIndex = new int[nTotalNode];
    PHSPACE::SetField(nodeBCTypeIndex, 0, nTotalNode);

    using namespace PHENGLEI;

    map<int, int> bcTypeMap;
    bcTypeMap.insert(pair<int,int>(SOLID_SURFACE   , 4));
    bcTypeMap.insert(pair<int,int>(FARFIELD        , 3));
    bcTypeMap.insert(pair<int,int>(SYMMETRY        , 2));
    bcTypeMap.insert(pair<int,int>(INTERFACE       , 1));
    bcTypeMap.insert(pair<int,int>(INFLOW          , 1));
    bcTypeMap.insert(pair<int,int>(OUTFLOW         , 1));
    bcTypeMap.insert(pair<int,int>(PRESSURE_INLET  , 1));
    bcTypeMap.insert(pair<int,int>(PRESSURE_OUTLET , 1));
    bcTypeMap.insert(pair<int,int>(OUTFLOW_CONFINED, 1));

    UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);

        int bcType = bcRegion->GetBCType();
        int bcTypeIndex = bcTypeMap[bcType];

        for (int iNode = static_cast<int>(face2NodeSubscript[iFace]); iNode < face2NodeSubscript[iFace + 1]; ++ iNode)
        {
            int nodeIndex = face2node[iNode];
            if (bcTypeIndex > nodeBCTypeIndex[nodeIndex])
            {
                nodeBCType[nodeIndex]      = bcType;
                nodeBCTypeIndex[nodeIndex] = bcTypeIndex;
            }
        }
    }

    delete [] nodeBCTypeIndex;    nodeBCTypeIndex = nullptr;
}

void UnstructGrid::CompGradientLSQ(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();
    RDouble *iwt = this->GetLeastSquareIWT();
    RDouble *ixx = this->GetLeastSquareIXX();
    RDouble *iyy = this->GetLeastSquareIYY();
    RDouble *izz = this->GetLeastSquareIZZ();
    RDouble *ixy = this->GetLeastSquareIXY();
    RDouble *ixz = this->GetLeastSquareIXZ();
    RDouble *iyz = this->GetLeastSquareIYZ();
    char *fMark  = this->GetFaceMark();

    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    int nTotalCell = this->GetNTotalCell();
    int nTotalFace = this->GetNTotalFace();
    int nBoundFace = this->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        dqdx[iCell] = 0.0;
        dqdy[iCell] = 0.0;
        dqdz[iCell] = 0.0;
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[ iFace ];
        int re = iFace + nTotalCell;
        RDouble dx = xcc[re] - xcc[le];
        RDouble dy = ycc[re] - ycc[le];
        RDouble dz = zcc[re] - zcc[le];
 
        RDouble dq   = iwt[iFace] * (q[re] - q[le]) * fMark[iFace];
        RDouble tmpx = dq * dx;
        RDouble tmpy = dq * dy;
        RDouble tmpz = dq * dz;
        dqdx[le] += tmpx;
        dqdy[le] += tmpy;
        dqdz[le] += tmpz;
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [ iFace ];
        int re = rightCellOfFace[ iFace ];

        RDouble dx = xcc[re] - xcc[le];
        RDouble dy = ycc[re] - ycc[le];
        RDouble dz = zcc[re] - zcc[le];
 
        RDouble dq = iwt[iFace] * (q[re] - q[le]) * fMark[iFace];

        RDouble tmpx = dq * dx;
        RDouble tmpy = dq * dy;
        RDouble tmpz = dq * dz;

        dqdx[le] += tmpx;
        dqdy[le] += tmpy;
        dqdz[le] += tmpz;

        dqdx[re] += tmpx;
        dqdy[re] += tmpy;
        dqdz[re] += tmpz;
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble tmpx    = dqdx[iCell];
        RDouble tmpy    = dqdy[iCell];
        RDouble tmpz    = dqdz[iCell];

        dqdx[iCell] = ixx[iCell] * tmpx + ixy[iCell] * tmpy + ixz[iCell] * tmpz;
        dqdy[iCell] = ixy[iCell] * tmpx + iyy[iCell] * tmpy + iyz[iCell] * tmpz;
        dqdz[iCell] = ixz[iCell] * tmpx + iyz[iCell] * tmpy + izz[iCell] * tmpz;
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[ iFace ];
        int re = iFace + nTotalCell;
        dqdx[re] = dqdx[le];
        dqdy[re] = dqdy[le];
        dqdz[re] = dqdz[le];
    }
}

void UnstructGrid::CompGradientGGNode(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    RDouble *faceNormalX = this->GetFaceNormalX();
    RDouble *faceNormalY = this->GetFaceNormalY();
    RDouble *faceNormalZ = this->GetFaceNormalZ();
    RDouble *faceArea  = this->GetFaceArea() ;
    RDouble *volume   = this->GetCellVolume();

    int *leftCellOfFace  = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    int *nodeNumberOfEachFace = this->GetNodeNumberOfEachFace();
    int *face2node = this->GetFace2Node();

    int nTotalNode = this->GetNTotalNode();
    int nTotalCell = this->GetNTotalCell();
    int nTotalFace = this->GetNTotalFace();
    int nBoundFace = this->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;    

    RDouble *qNode = new RDouble [nTotalNode]();

    string gradientName = GlobalDataBase::GetStrParaFromDB("gradientName");
    if (gradientName == "ggnode")
    {
        CompNodeVar(this, qNode, q);
    }
    else if (gradientName == "ggnodelaplacian")
    {
        CompNodeVar_new(this, qNode, q);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("gradientName", gradientName);
    }

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        dqdx[iCell] = 0.0;
        dqdy[iCell] = 0.0;
        dqdz[iCell] = 0.0;
    }

    int nodePosition = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace [iFace];

        RDouble qFaceCenter = 0.0;
        for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
        {
            int nodeID = face2node[nodePosition + jNode];
            qFaceCenter += qNode[nodeID];
        }
        nodePosition += nodeNumberOfEachFace[iFace];
        qFaceCenter /= nodeNumberOfEachFace[iFace];

        RDouble area = faceArea[iFace];
        RDouble nx = faceNormalX[iFace] * area;
        RDouble ny = faceNormalY[iFace] * area;
        RDouble nz = faceNormalZ[iFace] * area;

        dqdx[le] += qFaceCenter * nx;
        dqdy[le] += qFaceCenter * ny;
        dqdz[le] += qFaceCenter * nz;
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];

        RDouble qFaceCenter = 0.0;
        for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
        {
            int nodeID = face2node[nodePosition + jNode];
            qFaceCenter += qNode[nodeID];
        }
        nodePosition += nodeNumberOfEachFace[iFace];
        qFaceCenter /= nodeNumberOfEachFace[iFace];

        RDouble area = faceArea[iFace];
        RDouble nx = faceNormalX[iFace] * area;
        RDouble ny = faceNormalY[iFace] * area;
        RDouble nz = faceNormalZ[iFace] * area;

        dqdx[le] += qFaceCenter * nx;
        dqdy[le] += qFaceCenter * ny;
        dqdz[le] += qFaceCenter * nz;

        dqdx[re] -= qFaceCenter * nx;
        dqdy[re] -= qFaceCenter * ny;
        dqdz[re] -= qFaceCenter * nz;
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dqdx[iCell] /= volume[iCell];
        dqdy[iCell] /= volume[iCell];
        dqdz[iCell] /= volume[iCell];
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];

        dqdx[re] = dqdx[le];
        dqdy[re] = dqdy[le];
        dqdz[re] = dqdz[le];
    }

    delete [] qNode;    qNode = nullptr;
}

void UnstructGrid::CompGradientGGNode_NEW(RDouble *q, RDouble *qnode, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    RDouble *nxs, *nys, *nzs, *ns, *vol;

    nxs = this->GetFaceNormalX();
    nys = this->GetFaceNormalY();
    nzs = this->GetFaceNormalZ();
    ns  = this->GetFaceArea();
    vol = this->GetCellVolume();

    int *left_cell_of_face = this->GetLeftCellOfFace();
    int *right_cell_of_face = this->GetRightCellOfFace();

    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();
    int **face2nodeArray          = this->GetFace2NodeArray();

    int nTotalCell = this->GetNTotalCell();
    int nTotalFace = this->GetNTotalFace();
    int nBoundFace = this->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        dqdx[iCell] = 0.0;
        dqdy[iCell] = 0.0;
        dqdz[iCell] = 0.0;
    }

    UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = left_cell_of_face[ iFace ];
        int re = iFace + nTotalCell;

        RDouble qFace = 0.0;
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        if (bcType == PHENGLEI::SYMMETRY || bcType == PHENGLEI::INTERFACE || bcType == PHENGLEI::OVERSET)
        {
            for (int jNode = 0; jNode < node_number_of_each_face[iFace]; ++ jNode)
            {
                int index = face2nodeArray[iFace][jNode];
                qFace += qnode[index];
            }
            qFace /= node_number_of_each_face[iFace];
        }
        else
        {
            qFace = half * (q[le] + q[re]);
        }

        RDouble nx = nxs[iFace] * ns[iFace];
        RDouble ny = nys[iFace] * ns[iFace];
        RDouble nz = nzs[iFace] * ns[iFace];

        dqdx[le] += qFace * nx;
        dqdy[le] += qFace * ny;
        dqdz[le] += qFace * nz;
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = left_cell_of_face [ iFace ];
        int re = right_cell_of_face[ iFace ];

        RDouble qFace = 0.0;
        for (int jNode = 0; jNode < node_number_of_each_face[iFace]; ++ jNode)
        {
            int index = face2nodeArray[iFace][jNode];
            qFace += qnode[index];
        }
        qFace /= node_number_of_each_face[iFace];

        RDouble nx = nxs[iFace] * ns[iFace];
        RDouble ny = nys[iFace] * ns[iFace];
        RDouble nz = nzs[iFace] * ns[iFace];

        dqdx[le] += qFace * nx;
        dqdy[le] += qFace * ny;
        dqdz[le] += qFace * nz;

        dqdx[re] -= qFace * nx;
        dqdy[re] -= qFace * ny;
        dqdz[re] -= qFace * nz;
    } 


    //! Average volume.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dqdx[iCell] /= vol[iCell];
        dqdy[iCell] /= vol[iCell];
        dqdz[iCell] /= vol[iCell];
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();
        if (bcType == PHENGLEI::INTERFACE || bcType == PHENGLEI::OVERSET)
        {
            continue;
        }
        int le = left_cell_of_face[ iFace ];
        int re = iFace + nTotalCell;

        dqdx[re] = dqdx[le];
        dqdy[re] = dqdy[le];
        dqdz[re] = dqdz[le];
    }
}

void UnstructGrid::CompGradient(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    string gradientName = GlobalDataBase::GetStrParaFromDB("gradientName");
    if (gradientName == "ggnode" || gradientName == "ggnodelaplacian")
    {
        CompGradientGGNode(q, dqdx, dqdy, dqdz);
    }
    else if (gradientName == "ggcell")
    {
        CompGradientGGCell(q, dqdx, dqdy, dqdz);
    }
    else if (gradientName == "lsq")
    {
        CompGradientLSQ(q, dqdx, dqdy, dqdz);
    }
    else if (gradientName == "ggnode_weight")
    {
        CompGradientGGNodeWeight(q, dqdx, dqdy, dqdz);
    }
    else if (gradientName == "gg_m2")
    {
        CompGradientGGModified2(q, dqdx, dqdy, dqdz);
    }
    else if (gradientName == "ggcellnew")
    {
        CompGradientGGCellNew(q, dqdx, dqdy, dqdz);
    }
    else if (gradientName == "ggcellw")
    {
        CompGradientGGCellW(q, dqdx, dqdy, dqdz);
    }
    else
    {
        TK_Exit::ExceptionExit("No reconstruction method has been choosed ! /n");
    }
}

/************************************************************************/
/*             Green-Gausse node based method reconstruction            */
/*                       ggnode weighted by distance                    */
/*                          Bell 20121220 add                           */
/************************************************************************/
//! There is no distance weight when face to node.
void UnstructGrid::CompGradientGGNodeWeight(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    RDouble *nxs, *nys, *nzs, *ns, *vol;

    nxs = this->GetFaceNormalX();
    nys = this->GetFaceNormalY();
    nzs = this->GetFaceNormalZ();
    ns  = this->GetFaceArea();
    vol = this->GetCellVolume();

    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();
    int *face2node = this->GetFace2Node();

    int nTotalCell = this->GetNTotalCell();
    int nTotalFace = this->GetNTotalFace();
    int nBoundFace = this->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    int nTotalNode = this->GetNTotalNode();

    RDouble *q_n = new RDouble [nTotalNode];
    CompNodeVarWeight(this, q_n, q);

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        dqdx[iCell] = 0.0;
        dqdy[iCell] = 0.0;
        dqdz[iCell] = 0.0;
    }

    int nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];

        RDouble qfc = 0.0;
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            int index = face2node[nodepos + j];
            qfc += q_n[index];
        }
        nodepos += node_number_of_each_face[iFace];
        qfc /= node_number_of_each_face[iFace];

        RDouble nx = nxs[iFace] * ns[iFace];
        RDouble ny = nys[iFace] * ns[iFace];
        RDouble nz = nzs[iFace] * ns[iFace];

        dqdx[le] += qfc * nx;
        dqdy[le] += qfc * ny;
        dqdz[le] += qfc * nz;
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [ iFace ];
        int re = rightCellOfFace[ iFace ];

        RDouble qfc = 0.0;
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            int index = face2node[nodepos + j];
            qfc += q_n[index];
        }
        nodepos += node_number_of_each_face[iFace];
        qfc /= node_number_of_each_face[iFace];

        RDouble nx = nxs[iFace] * ns[iFace];
        RDouble ny = nys[iFace] * ns[iFace];
        RDouble nz = nzs[iFace] * ns[iFace];

        dqdx[le] += qfc * nx;
        dqdy[le] += qfc * ny;
        dqdz[le] += qfc * nz;

        dqdx[re] -= qfc * nx;
        dqdy[re] -= qfc * ny;
        dqdz[re] -= qfc * nz;
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dqdx[iCell] /= vol[iCell];
        dqdy[iCell] /= vol[iCell];
        dqdz[iCell] /= vol[iCell];
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = iFace + nTotalCell;

        dqdx[re] = dqdx[le];
        dqdy[re] = dqdy[le];
        dqdz[re] = dqdz[le];
    }

    delete [] q_n;    q_n = nullptr;
}

void UnstructGrid::GetCellCenter(ActionKey *actkey)
{
    InterfaceInfo *interfaceInfo = this->GetInterfaceInfo();
    if (!interfaceInfo) return;

    int ineighbor = interfaceInfo->FindIthNeighbor(actkey->ipos);
    int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
    int *faceIndexForSend = interfaceInfo->GetFaceIndexForSend(ineighbor);

    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();

    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();
    RDouble *vol = this->GetCellVolume();

    int nm = 4;
    int nlen = nm * nIFaceOfNeighbor;
    RDouble *qtmp = new RDouble [nlen];
    int ntmp = 0;

    int iFace,s1;
    for (int iLocalFace = 0; iLocalFace < nIFaceOfNeighbor; ++ iLocalFace)
    {
        iFace = faceIndexForSend[iLocalFace];
        this->GetSourceIndex(iFace, 1, s1);

        qtmp[ntmp++] = xcc[s1];
        qtmp[ntmp++] = ycc[s1];
        qtmp[ntmp++] = zcc[s1];
        qtmp[ntmp++] = vol[s1];
    }
    cdata->Write(qtmp, nlen * sizeof(RDouble));
    delete [] qtmp;    qtmp = nullptr;
}

void UnstructGrid::TranslateCellCenter(ActionKey *actkey)
{
    InterfaceInfo *interfaceInfo = this->GetInterfaceInfo();
    if (!interfaceInfo) return;

    int ineighbor = interfaceInfo->FindIthNeighbor(actkey->ipos);
    int iFace, t1;

    int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
    int *faceIndexForRecv = interfaceInfo->GetFaceIndexForRecv(ineighbor);

    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();
    RDouble *vol = this->GetCellVolume();

    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();

    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");
    RDouble translationLength[3] = { 0.0 };
    GlobalDataBase::GetData("translationLength", &translationLength, PHDOUBLE, 3);
    RDouble rotationAngle = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("rotationAngle");
    rotationAngle = rotationAngle * PI / 180.0;
    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");

    int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();
    UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();
    int *leftCellofFace = this->GetLeftCellOfFace();
    int *rightCellofFace = this->GetRightCellOfFace();
    int count = 0;

    for (int iLocalFace = 0; iLocalFace < nIFaceOfNeighbor; ++ iLocalFace)
    {
        iFace = faceIndexForRecv[iLocalFace];
        this->GetTargetIndex(iFace, 1, t1);

        cdata->Read(&xcc[t1], sizeof(RDouble));
        cdata->Read(&ycc[t1], sizeof(RDouble));
        cdata->Read(&zcc[t1], sizeof(RDouble));
        cdata->Read(&vol[t1], sizeof(RDouble));

        int iBFace = interFace2BoundaryFace[iFace];
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iBFace]);
        int bcType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();

        if (referenceFrame == ROTATIONAL_FRAME)
        {
            //! nTurboZone
            int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
            //! PeriodicRotationAngle
            RDouble PeriodicRotationAngle[100];
            GlobalDataBase::GetData("PeriodicRotationAngle", &PeriodicRotationAngle, PHDOUBLE, nTurboZone);

            //! Periodic_Name
            string Periodic_Name[100];
            GlobalDataBase::GetData("Periodic_Name", &Periodic_Name, PHSTRING, 2 * nTurboZone);

            for (int iTurboZone = 0; iTurboZone < nTurboZone; iTurboZone++)
            {
                PeriodicRotationAngle[iTurboZone] = PeriodicRotationAngle[iTurboZone] * PI / 180.0;
                rotationAngle = PeriodicRotationAngle[iTurboZone];

                RDouble rotcell[3] = { 0.0, 0.0, 0.0 };
                if (bcName == Periodic_Name[2 * iTurboZone])
                {
                    rotcell[1] = ycc[t1] * cos(2.0 * PI - rotationAngle) - zcc[t1] * sin(2.0 * PI - rotationAngle);
                    rotcell[2] = ycc[t1] * sin(2.0 * PI - rotationAngle) + zcc[t1] * cos(2.0 * PI - rotationAngle);

                    ycc[t1] = rotcell[1];
                    zcc[t1] = rotcell[2];
                }
                else if (bcName == Periodic_Name[2 * iTurboZone + 1])
                {
                    rotcell[1] = ycc[t1] * cos(rotationAngle) - zcc[t1] * sin(rotationAngle);
                    rotcell[2] = ycc[t1] * sin(rotationAngle) + zcc[t1] * cos(rotationAngle);

                    ycc[t1] = rotcell[1];
                    zcc[t1] = rotcell[2];
                }
            }
        }
        else
        {
            if (periodicType == TRANSLATIONAL_PERIODICITY)
            {
                if (bcName == "Periodic_up")
                {
                    xcc[t1] = xcc[t1] - translationLength[0];
                    ycc[t1] = ycc[t1] - translationLength[1];
                    zcc[t1] = zcc[t1] - translationLength[2];
                }
                else if (bcName == "Periodic_down")
                {
                    xcc[t1] = xcc[t1] + translationLength[0];
                    ycc[t1] = ycc[t1] + translationLength[1];
                    zcc[t1] = zcc[t1] + translationLength[2];
                }
            }
            else if (periodicType == ROTATIONAL_PERIODICITY)
            {
                RDouble rotcell[3] = { 0.0, 0.0, 0.0 };

                if (bcName == "Periodic_up")
                {   //! cell center coordinate
                    rotcell[1] = ycc[t1] * cos(2.0 * PI - rotationAngle) - zcc[t1] * sin(2.0 * PI - rotationAngle);
                    rotcell[2] = ycc[t1] * sin(2.0 * PI - rotationAngle) + zcc[t1] * cos(2.0 * PI - rotationAngle);

                    //! coordinate matching
                    ycc[t1] = rotcell[1];
                    zcc[t1] = rotcell[2];
                }
                else if (bcName == "Periodic_down")
                {
                    rotcell[1] = ycc[t1] * cos(rotationAngle) - zcc[t1] * sin(rotationAngle);
                    rotcell[2] = ycc[t1] * sin(rotationAngle) + zcc[t1] * cos(rotationAngle);

                    //! coordinate matching
                    ycc[t1] = rotcell[1];
                    zcc[t1] = rotcell[2];
                }
            }
        }
    }
}

void UnstructGrid::GetCellIBlank(ActionKey *actkey)
{
    InterfaceInfo *interfaceInfo = this->GetInterfaceInfo();
    if (!interfaceInfo) return;

    int ineighbor = interfaceInfo->FindIthNeighbor(actkey->ipos);
    int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
    int *faceIndexForSend = interfaceInfo->GetFaceIndexForSend(ineighbor);

    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();

    int *iBlank = this->GetBlankIndex();

    int *qtmp = new int [nIFaceOfNeighbor];

    int iFace,s1;
    for (int iLocalFace = 0; iLocalFace < nIFaceOfNeighbor; ++ iLocalFace)
    {
        iFace = faceIndexForSend[iLocalFace];
        this->GetSourceIndex(iFace, 1, s1);

        qtmp[iLocalFace] = iBlank[s1];
    }
    cdata->Write(qtmp, nIFaceOfNeighbor * sizeof(int));
    delete [] qtmp;    qtmp = nullptr;
}

void UnstructGrid::TranslateCellIBlank(ActionKey *actkey)
{
    InterfaceInfo *interfaceInfo = this->GetInterfaceInfo();
    if (!interfaceInfo) return;

    int ineighbor = interfaceInfo->FindIthNeighbor(actkey->ipos);
    int iFace, t1;

    int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
    int *faceIndexForRecv = interfaceInfo->GetFaceIndexForRecv(ineighbor);

    int *iBlank = this->GetBlankIndex();

    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();

    for (int iLocalFace = 0; iLocalFace < nIFaceOfNeighbor; ++ iLocalFace)
    {
        iFace = faceIndexForRecv[iLocalFace];
        this->GetTargetIndex(iFace, 1, t1);

        cdata->Read(&iBlank[t1], sizeof(int));
    }
}

streamsize UnstructGrid::TranslateCellCenterLength(ActionKey *actkey)
{
    streamsize nlen = 0;

    InterfaceInfo *interfaceInfo = this->GetInterfaceInfo();
    if (!interfaceInfo) return nlen;

    int ineighbor = interfaceInfo->FindIthNeighbor(actkey->ipos);
    int iFace,t1;

    int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
    int *faceIndexForRecv = interfaceInfo->GetFaceIndexForRecv(ineighbor);

    for (int iLocalFace = 0; iLocalFace < nIFaceOfNeighbor; ++ iLocalFace)
    {
        iFace = faceIndexForRecv[iLocalFace];
        this->GetTargetIndex(iFace, 1, t1);

        nlen += sizeof(RDouble) * 4;
    }

    return nlen;
}

streamsize UnstructGrid::TranslateCellIBlankLength(ActionKey *actkey)
{
    streamsize nlen = 0;

    InterfaceInfo *interfaceInfo = this->GetInterfaceInfo();
    if (!interfaceInfo) return nlen;

    int ineighbor = interfaceInfo->FindIthNeighbor(actkey->ipos);
    int iFace,t1;

    int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
    int *faceIndexForRecv = interfaceInfo->GetFaceIndexForRecv(ineighbor);

    for (int iLocalFace = 0; iLocalFace < nIFaceOfNeighbor; ++ iLocalFace)
    {
        iFace = faceIndexForRecv[iLocalFace];
        this->GetTargetIndex(iFace, 1, t1);

        nlen += sizeof(int) * 1;
    }

    return nlen;
}

void UnstructGrid::CompressSpecificArrayToInterface(DataContainer *&dataContainer, const string &fieldName, const int &neighborZoneIndex, const int &nEquation)
{
    InterfaceInfo *interfaceInformation = this->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }
    if (fieldName == "gradUVWTCellCenterX"||fieldName == "gradUVWTCellCenterY"||fieldName == "gradUVWTCellCenterZ")
    {
        return;
    }
    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
    dataContainer->MoveToBegin();
    if (nEquation == 0)
    {
        RDouble *fieldSend = reinterpret_cast <RDouble *> (this->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int sourceCell;
                int iFace = interfaceIndexContainerForSend[iLocalFace];
                this->GetSourceIndex(iFace, iGhostLayer + 1, sourceCell);
                PHWrite(dataContainer, fieldSend[sourceCell]);
            }
        }
    }
    else
    {
        RDouble **fieldSend = reinterpret_cast <RDouble **> (this->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int sourceCell;
                int iFace = interfaceIndexContainerForSend[iLocalFace];
                this->GetSourceIndex(iFace, iGhostLayer + 1, sourceCell);
                for (int m = 0; m < nEquation; ++ m)
                {
                    PHWrite(dataContainer, fieldSend[m][sourceCell]);
                }
            }
        }
    }
}

void UnstructGrid::DecompressArrayFromInterface(DataContainer *&dataContainer, const string &fieldName, const int &neighborZoneIndex, const int &nEquation)
{
    InterfaceInfo *interfaceInformation = this->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }
    if (fieldName == "gradUVWTCellCenterX"||fieldName == "gradUVWTCellCenterY"||fieldName == "gradUVWTCellCenterZ")
    {
        return;
    }
    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);
    dataContainer->MoveToBegin();
    if (nEquation == 0)
    {
        RDouble *fieldRecv = reinterpret_cast <RDouble *> (this->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int targetCell;
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                this->GetTargetIndex(iFace, iGhostLayer + 1, targetCell);
                PHRead(dataContainer, fieldRecv[targetCell]);
            }
        }
    }
    else
    {
        RDouble **fieldRecv = reinterpret_cast <RDouble **> (this->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                int targetCell;
                this->GetTargetIndex(iFace, iGhostLayer + 1, targetCell);
                for (int m = 0; m < nEquation; ++ m)
                {
                    PHRead(dataContainer, fieldRecv[m][targetCell]);
                }
            }
        }
    }
}

void UnstructGrid::SetOrdinaryGridIndex(int ordinaryGridIndexIn)
{
    this->ordinaryGridIndex = ordinaryGridIndexIn;
}

void UnstructGrid::SetOrdinaryNodeIndex(int *ordinaryNodeIndexIn)
{
    if (this->ordinaryNodeIndex != nullptr)
    {
        delete [] this->ordinaryNodeIndex;    ordinaryNodeIndex = nullptr;
    }

    this->ordinaryNodeIndex = ordinaryNodeIndexIn;
}

void UnstructGrid::SetOrdinaryFaceIndex(int *ordinaryFaceIndexIn)
{
    if (this->ordinaryFaceIndex != nullptr)
    {
        delete [] this->ordinaryFaceIndex;    ordinaryFaceIndex = nullptr;
    }

    this->ordinaryFaceIndex = ordinaryFaceIndexIn;
}

void UnstructGrid::SetOrdinaryCellIndex(int *ordinaryCellIndexIn)
{
    if (this->ordinaryCellIndex != nullptr)
    {
        delete [] this->ordinaryCellIndex;    ordinaryCellIndex = nullptr;
    }

    this->ordinaryCellIndex = ordinaryCellIndexIn;
}

int  UnstructGrid::GetOrdinaryGridIndex()
{
    return this->ordinaryGridIndex;
}

int *UnstructGrid::GetOrdinaryNodeIndex()
{
    return this->ordinaryNodeIndex;
}

int *UnstructGrid::GetOrdinaryFaceIndex()
{
    return this->ordinaryFaceIndex;
}

int *UnstructGrid::GetOrdinaryCellIndex()
{
    return this->ordinaryCellIndex;
}

LIB_EXPORT void UnstructGrid::FindOutWallPoints(bool * isWallPoint, int & nWallPoints)
{
    UnstructBCSet **bcr = this->GetBCRecord();

    int nBoundFace = this->GetNBoundFace();
    int nTotalNode = this->GetNTotalNode();

    int *face2node = this->GetFace2Node();
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        isWallPoint[iNode] = false;
    }
    int nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcr[iFace]->GetKey();
        if (IsWall(bcType))
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                int index = face2node[ nodepos + j ];
                isWallPoint[index] = true;
            }
        }

        nodepos += node_number_of_each_face[iFace];
    }

    nWallPoints = 0;
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (isWallPoint[iNode])
        {
            ++ nWallPoints;
        }
    }
}

LIB_EXPORT void UnstructGrid::FindOutWallFaces(bool * isWallFace, int & nWallFaces)
{
    UnstructBCSet **bcr = this->GetBCRecord();
    int nBoundFace = this->GetNBoundFace();

    nWallFaces = 0;

    //! Compute the number of the wall faces.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        isWallFace[iFace] = false;

        int bctype = bcr[iFace]->GetKey();
        if (IsWall(bctype))
        {
            isWallFace[iFace] = true;
            ++ nWallFaces;
        }
    }
}

void UnstructGrid::AllocateWalldist()
{
    int nTotalNode = this->GetNTotalNode();
    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    walldist = new RDouble[nTotal];
    walldistNode = new RDouble[nTotalNode];
    nearestwallfacenormalx = new RDouble[nTotal];
    nearestwallfacenormaly = new RDouble[nTotal];
    nearestwallfacenormalz = new RDouble[nTotal];
}

void UnstructGrid::ReadWalldist(RDouble *wallDistIn)
{
    int nTotalCell = this->GetNTotalCell();

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        walldist[iCell] = wallDistIn[iCell];
    }
}

void UnstructGrid::ReadNearestwallfacenormalx(RDouble *nearestwallfaceNormalxIn)
{
    int nTotalCell = this->GetNTotalCell();

    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        nearestwallfacenormalx[iCell] = nearestwallfaceNormalxIn[iCell];
    }
}

void UnstructGrid::ReadNearestwallfacenormaly(RDouble* nearestwallfaceNormalyIn)
{
    int nTotalCell = this->GetNTotalCell();

    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        nearestwallfacenormaly[iCell] = nearestwallfaceNormalyIn[iCell];
    }
}

void UnstructGrid::ReadNearestwallfacenormalz(RDouble* nearestwallfaceNormalzIn)
{
    int nTotalCell = this->GetNTotalCell();

    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        nearestwallfacenormalz[iCell] = nearestwallfaceNormalzIn[iCell];
    }
}

void UnstructGrid::DumpWalldist(ActionKey *actkey)
{
    int nTotalCell = this->GetNTotalCell();

    DataContainer *cdata = actkey->GetData();
    cdata->MoveToBegin();

    cdata->Write(&nTotalCell, sizeof(int));
    cdata->Write(walldist,sizeof(RDouble) * nTotalCell);
    cdata->Write(nearestwallfacenormalx, sizeof(RDouble) * nTotalCell);
    cdata->Write(nearestwallfacenormaly, sizeof(RDouble) * nTotalCell);
    cdata->Write(nearestwallfacenormalz, sizeof(RDouble) * nTotalCell);
}

void UnstructGrid::GetBoundaryMinMax(RDouble *pmin, RDouble *pmax)
{
    using namespace PHSPACE;

    int nBoundFace = this->GetNBoundFace();

    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();

    pmin[0] = LARGE;
    pmin[1] = LARGE;
    pmin[2] = LARGE;

    pmax[0] = - LARGE;
    pmax[1] = - LARGE;
    pmax[2] = - LARGE;

    int *face2node  = this->GetFace2Node();
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();

    int nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            int index = face2node[ nodepos + j ];
            pmin[0] = MIN(static_cast<RDouble>(pmin[0]), x[index]);
            pmin[1] = MIN(static_cast<RDouble>(pmin[1]), y[index]);
            pmin[2] = MIN(static_cast<RDouble>(pmin[2]), z[index]);

            pmax[0] = MAX(static_cast<RDouble>(pmax[0]), x[index]);
            pmax[1] = MAX(static_cast<RDouble>(pmax[1]), y[index]);
            pmax[2] = MAX(static_cast<RDouble>(pmax[2]), z[index]);
        }
        nodepos += node_number_of_each_face[iFace];
    }
}

void UnstructGrid::GetMinMaxDS(RDouble &dismin, RDouble &dismax)
{
    using namespace PHSPACE;
    //! Check the minimum grid distance.
    dismin =   LARGE;
    dismax = - LARGE;

    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();

    int nTotalFace = this->GetNTotalFace();
    int *face2node = this->GetFace2Node();
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();
    int nodepos = 0;

    RDouble dx,dy,dz,ds;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        for (int iNode = 0; iNode < node_number_of_each_face[iFace]; ++ iNode)
        {
            int p1 = face2node[ nodepos + iNode ];
            int p2 = face2node[ nodepos + (iNode + 1)%node_number_of_each_face[iFace] ];
            dx = x[ p2 ] - x[ p1 ];
            dy = y[ p2 ] - y[ p1 ];
            dz = z[ p2 ] - z[ p1 ];
            ds = sqrt(dx * dx + dy * dy + dz * dz);

            if (ds <= SAME_POINT_TOL) continue;
            dismin = MIN(dismin, ds);
            dismax = MAX(dismin, ds);
        }
        nodepos += node_number_of_each_face[iFace];
    }
}

void UnstructGrid::Encode(DataContainer *cdata, int flag)
{
    if (flag == 0)
    {
        EncodeGrid(cdata);
    }
    else if (flag == 1)
    {
        EncodeOversetGrid(cdata);
    }
}

void UnstructGrid::Decode(DataContainer *cdata, int flag)
{
    if (flag == 0)
    {
        DecodeGrid(cdata);
    }
    else if (flag == 1)
    {
        DecodeOversetGrid(cdata);
    }
}

void UnstructGrid::ReviseNeighborZoneIndex(int zoneStart)
{
    InterfaceInfo *interfaceInfo = this->GetInterfaceInfo();

    if(! interfaceInfo) return;

    int nIFace = interfaceInfo->GetNIFace();

    int *interFace2ZoneID = interfaceInfo->GetInterFace2ZoneID();

    for(int iFace = 0; iFace < nIFace; ++ iFace)
    {
        interFace2ZoneID[ iFace ] += zoneStart;
    }

    InterpointInformation *interpointInformation = this->GetInterpointInfo();

    if(! interpointInformation) return;

    int nIPoint = interpointInformation->GetNumberOfInterpoints();

    int *interPoint2ZoneID = interpointInformation->GetInterPoint2ZoneID();

    for(int iPoint = 0; iPoint < nIPoint; ++ iPoint)
    {
        interPoint2ZoneID[iPoint] += zoneStart;
    }
}

void UnstructGrid::EncodeOversetGrid(DataContainer *cdata)
{
    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;    //! interior cells + ghost cells.

    int *iBlank = this->iBlank;

    cdata->MoveToBegin();

    cdata->Write(reinterpret_cast<char *>(&nTotal),sizeof(int)); 

    cdata->Write(iBlank, nTotal * sizeof(int));

    this->GetOversetInformationProxy()->EncodeOversetInformation(cdata);
}

void UnstructGrid::DecodeOversetGrid(DataContainer *cdata)
{
    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;    //! interior cells + ghost cells.

    if (!iBlank) iBlank = new int[nTotal];

    int *iBlank = this->iBlank;

    cdata->MoveToBegin();

    cdata->Read(reinterpret_cast<char *>(&nTotal),sizeof(int)); 

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        cdata->Read(reinterpret_cast<char *>(&iBlank[iCell]),sizeof(int)); 
    }

    OversetInformationProxy * oversetInformationProxy = new OversetInformationProxy();

    this->SetOversetInformationProxy(oversetInformationProxy);

    oversetInformationProxy->DecodeOversetInformation(cdata);
}

void UnstructGrid::DecodeWallDist(DataContainer * cdata)
{
    VirtualFile *virtualFile = new VirtualFile(cdata);

    virtualFile->BeginReadWork();

    int nTotalCell = this->GetNTotalCell();
    
    RDouble *walldist = GetWallDist();

    PHRead(virtualFile, walldist, nTotalCell);

    virtualFile->EndReadWork();

    delete virtualFile;    virtualFile = nullptr;
}

void UnstructGrid::DecodeCellNode(DataContainer * cdata)
{
    VirtualFile *virtualFile = new VirtualFile(cdata);

    virtualFile->BeginReadWork();

    int zoneIndex;
    PHRead(virtualFile, zoneIndex);

    int numberOfCells;
    PHRead(virtualFile, numberOfCells);

    vector< vector< int > > cellToNode;
    cellToNode.resize(numberOfCells);

    int count = 0;
    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        int nodeNumber;
        PHRead(virtualFile, nodeNumber);

        cellToNode[ iCell ].resize(nodeNumber);

        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            PHRead(virtualFile, cellToNode[ iCell ][ iNode ]);
            ++ count;
        }
    }

    virtualFile->EndReadWork();

    int *cellNodeNumberContainer = new int [ numberOfCells ];
    int *cellNodeIndexContainer  = new int [ count ];

    count = 0;
    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        int nodeNumber = static_cast<int>(cellToNode[ iCell ].size());
        cellNodeNumberContainer[ iCell ] = nodeNumber;

        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            cellNodeIndexContainer[ count ++ ] = cellToNode[ iCell ][ iNode ];
        }
    }

    delete virtualFile;    virtualFile = nullptr;

    this->SetNodeNumberOfEachCell(cellNodeNumberContainer);
    this->SetCell2Node(cellNodeIndexContainer);
}

void UnstructGrid::DecodeInterpointInfo(DataContainer *cdata)
{
    InterpointInformation *interpointInformation = 0;
    this->SetInterpointInfo(interpointInformation);
    
    VirtualFile *virtualFile = new VirtualFile(cdata);
    
    virtualFile->BeginReadWork();
    
    int iZone;
    PHRead(virtualFile, iZone);

    int i_nIPoint; 
    PHRead(virtualFile, i_nIPoint);
    
    if (i_nIPoint > 0)
    {
        interpointInformation = new InterpointInformation(i_nIPoint,this);
        this->SetInterpointInfo(interpointInformation);
        
        int *interPoint2ZoneID = interpointInformation->GetInterPoint2ZoneID();
        int *interPoint2InterPointID = interpointInformation->GetInterPoint2InterPointID();
        int *interPoint2GlobalPoint = interpointInformation->GetInterPoint2GlobalPoint();
        int *cellNumberOfInterPoint = interpointInformation->GetCellNumberOfInterPoint();
        int *totalZonesOfInterPoint = interpointInformation->GetTotalZonesOfInterPoint();
        int *labelOfInterPoint = interpointInformation->GetLabelOfInterPoint();

        for (int iPoint = 0; iPoint < i_nIPoint; iPoint++)
        {
            PHRead(virtualFile, interPoint2ZoneID[iPoint]);
            PHRead(virtualFile, interPoint2InterPointID[iPoint]);
            PHRead(virtualFile, interPoint2GlobalPoint[iPoint]);
            PHRead(virtualFile, cellNumberOfInterPoint[iPoint]);
            PHRead(virtualFile, totalZonesOfInterPoint[iPoint]);
            PHRead(virtualFile, labelOfInterPoint[iPoint]);
        }
    }

    virtualFile->EndReadWork();

    delete virtualFile;    virtualFile = nullptr;
}

void UnstructGrid::DecodeBCFace(DataContainer * cdata)
{
    VirtualFile *virtualFile = new VirtualFile(cdata);

    virtualFile->BeginReadWork();

    int zoneIndex;
    PHRead(virtualFile, zoneIndex);

    int numberOfBCFaces;
    PHRead(virtualFile, numberOfBCFaces);
    ASSERT(numberOfBCFaces == this->GetNBoundFace());

    uint_long totalSize;
    PHRead(virtualFile, totalSize);

    char *bcNameChar = new char[totalSize];
    PHRead(virtualFile, bcNameChar, static_cast<int>(totalSize));
    
    string *bcNameList = new string[numberOfBCFaces];
    int count = 0;
    for (int iFace = 0; iFace < numberOfBCFaces; ++ iFace)
    {
        if(count < totalSize)
        {
            while (bcNameChar[count] != '\0')
            {
                bcNameList[iFace].push_back(bcNameChar[count]);
                ++ count;
            }
            ++ count;
        }
    }

    PHRead(virtualFile, numberOfBCFaces);
    int *bcTypeList = new int[numberOfBCFaces];
    PHRead(virtualFile, bcTypeList, numberOfBCFaces);

    virtualFile->EndReadWork();

    //! Create UnstructbcRegion.
    set<string> bcNameSet;
    UnstructBCSet **bcr = this->GetBCRecord();
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        bcr[iFace]->SetBCName(bcNameList[iFace]);
        bcNameSet.insert(bcNameList[iFace]);
    }

    int nBCRegionUnstruct = static_cast<int>(bcNameSet.size());
    this->CreateUnstructBCSet(nBCRegionUnstruct);

    UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
    int *bcRegionIDofBCFace = new int[nBoundFace];

    set<string>::iterator iter;
    count = 0;
    for (iter = bcNameSet.begin(); iter != bcNameSet.end(); ++ iter)
    {
        UnstructBC *unstructBC = new UnstructBC(count);
        unstructBCSet->SetBCRegion(count, unstructBC);
        unstructBC->SetBCName(*iter);

        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            if (bcNameList[iFace] == unstructBC->GetBCName())
            {
                unstructBC->SetBCType(bcTypeList[iFace]);
                bcRegionIDofBCFace[iFace] = count;
                vector<int> *faceIndex = unstructBC->GetFaceIndex();
                faceIndex->push_back(iFace);
            }
        }
        count ++;
    }
    unstructBCSet->SetBCFaceInfo(bcRegionIDofBCFace);
 
    delete [] bcTypeList;    bcTypeList = nullptr;
    delete [] bcNameList;    bcNameList = nullptr;
    delete [] bcNameChar;    bcNameChar = nullptr;
    delete virtualFile;    virtualFile = nullptr;
}

void UnstructGrid::UpdateCoordinate(DataContainer *cdata)
{
    cdata->MoveToBegin();
    cout << "读入非结构网格动网格数据文件......\n";

    int nTotalNode = this->GetNTotalNode();

    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();

    cdata->Read(reinterpret_cast<char *>(x), nTotalNode * sizeof(RDouble));
    cdata->Read(reinterpret_cast<char *>(y), nTotalNode * sizeof(RDouble));
    cdata->Read(reinterpret_cast<char *>(z), nTotalNode * sizeof(RDouble));

    ComputeMinMaxBox();
}

LIB_EXPORT void UnstructGrid::AllocateOversetGrid()
{
    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;    //! interior cells + ghost cells.

    if ((this->iBlank))
    {
        delete []iBlank;
        iBlank = nullptr;
    }
    this->iBlank = new int[nTotal];

    int *iBlank = this->iBlank;

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        iBlank[iCell] = 1;
    }
}

LIB_EXPORT void UnstructGrid::CalNormalDistanceOfC2C() 
{
    if (normalDistanceC2C != nullptr)
    {
        cout << "norm_dist_c2c has been existed." << endl;
        return;
    }

    int le, re;
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();
    RDouble *xfn = this->GetFaceNormalX();
    RDouble *yfn = this->GetFaceNormalY();
    RDouble *zfn = this->GetFaceNormalZ();
    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();

    normalDistanceC2C = NewPointer<RDouble>(nTotalFace);
    for (int iFace = 0; iFace < nTotalFace; iFace++) 
    {
        le = leftCellOfFace[iFace];
        re = rightCellOfFace[iFace];

        normalDistanceC2C[iFace] = ABS(xfn[iFace] * (xcc[re] - xcc[le])
            + yfn[iFace] * (ycc[re] - ycc[le]) + zfn[iFace] * (zcc[re] - zcc[le]));
    }
}
void UnstructGrid::ReadGrid(fstream &file)
{
    VirtualFile *vfile = new VirtualFile(&file);

    vfile->BeginReadWork();

    ReadGrid(vfile);

    vfile->EndReadWork();

    delete vfile;    vfile = nullptr;
}

LIB_EXPORT void UnstructGrid::WriteGrid(fstream &file)
{
    VirtualFile *vfile = new VirtualFile(&file);

    vfile->BeginWriteWork();

    WriteGrid(vfile);

    vfile->EndWriteWork();

    delete vfile;    vfile = nullptr;
}

void UnstructGrid::DecodeGrid(DataContainer *cdata)
{
    VirtualFile *vfile = new VirtualFile(cdata);

    vfile->BeginReadWork();

    ReadGrid(vfile);

    vfile->EndReadWork();

    delete vfile;    vfile = nullptr;
}

void UnstructGrid::EncodeGrid(DataContainer *cdata)
{
    VirtualFile *vfile = new VirtualFile(cdata);

    vfile->BeginWriteWork();

    WriteGrid(vfile);

    vfile->EndWriteWork();

    delete vfile;    vfile = nullptr;
}

void UnstructGrid::ReadGrid(VirtualFile *vfile)
{
    int nTotalNode, nTotalFace, nTotalCell, nIFace, nBoundFace;

    PrintToWindow("Reading Unstructured Grid of Zone ", GetZoneID(), "...\n");
    //! Read the total nodes,total faces and total cells.

    vfile->read(reinterpret_cast<char*>(&nTotalNode),sizeof(int));
    vfile->read(reinterpret_cast<char*>(&nTotalFace),sizeof(int));
    vfile->read(reinterpret_cast<char*>(&nTotalCell),sizeof(int));

    PrintToWindow("  Grid Dimension  : ", this->GetDim(), "\n");
    PrintToWindow("  Number of Points: ", nTotalNode, "\n");
    PrintToWindow("  Number of Faces : ", nTotalFace, "\n");
    PrintToWindow("  Number of Cells : ", nTotalCell, "\n");

    this->SetNTotalNode(nTotalNode);
    this->SetNTotalFace(nTotalFace);
    this->SetNTotalCell(nTotalCell);

    RDouble *x, *y, *z;
    x = new RDouble[ nTotalNode ];
    y = new RDouble[ nTotalNode ];
    z = new RDouble[ nTotalNode ];

    this->SetX(x);
    this->SetY(y);
    this->SetZ(z);

    vfile->read(reinterpret_cast<char*>(x),nTotalNode * sizeof(RDouble));
    vfile->read(reinterpret_cast<char*>(y),nTotalNode * sizeof(RDouble));
    vfile->read(reinterpret_cast<char*>(z),nTotalNode * sizeof(RDouble));

    RotateAboutAxis();

    ComputeMinMaxBox();

    int *node_number_of_each_face = new int[ nTotalFace ];
    //! Assign a value to node_number_of_each_face.
    vfile->read(reinterpret_cast<char*>(node_number_of_each_face),nTotalFace * sizeof(int));
    this->SetNodeNumberOfEachFace(node_number_of_each_face);

    int *nodePosi = new int[ nTotalFace + 1 ];
    int nsum = 0;
    nodePosi[0] = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        nsum += node_number_of_each_face[iFace];
        nodePosi[iFace+1] = nsum;
    }

    //! The connected relationship of face to node.
    int *face2node = new int [ nsum ];
    this->SetFace2Node(face2node);

    vfile->read(reinterpret_cast<char*>(face2node), nsum * sizeof(int));

    int *leftCellOfFace  = new int[ nTotalFace ];
    int *rightCellOfFace = new int[ nTotalFace ];
    this->SetLeftCellOfFace (leftCellOfFace);
    this->SetRightCellOfFace(rightCellOfFace);

    vfile->read(reinterpret_cast<char*>(leftCellOfFace),nTotalFace*sizeof(int));
    vfile->read(reinterpret_cast<char*>(rightCellOfFace),nTotalFace*sizeof(int));

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        if (leftCellOfFace[ iFace ] < 0)
        {
            //! Need to reverse the node ordering.
            //! This problem is very hidden,mainly cause by the unclear definition of reverse.
            //! The bottom line is wrong.
            //std::reverse(face2node + nodePosi[iFace], face2node + nodePosi[iFace+1] - 1);
            std::reverse(face2node + nodePosi[iFace], face2node + nodePosi[iFace+1]);

            //! now reverse le and re
            SWAP(leftCellOfFace[ iFace ], rightCellOfFace[ iFace ]);
        }


    }
    delete [] nodePosi;    nodePosi = nullptr;

    vfile->read(reinterpret_cast<char*>(&nBoundFace),sizeof(int));
    this->SetNBoundFace(nBoundFace);

    PrintToWindow("  Number of Boundary Faces : ", nBoundFace, "\n");

    //! Set the boundary condition.
    UnstructBCSet **bcr = new UnstructBCSet * [ nBoundFace ];
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bctype;
        vfile->read(reinterpret_cast<char*>(&bctype),sizeof(int));
        bcr[iFace] = new UnstructBCSet();
        bcr[iFace]->SetKey(bctype);
    }
    this->SetBCRecord(bcr);
    vfile->read(reinterpret_cast<char*>(&nIFace),sizeof(int));

    InterfaceInfo *interfaceInfo = 0;
    this->SetInterfaceInfo(interfaceInfo);
    if (nIFace > 0)
    {
        interfaceInfo = new InterfaceInfo(nIFace,this);
        this->SetInterfaceInfo(interfaceInfo);

        int *interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
        int *interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();
        int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

        vfile->read(reinterpret_cast<char*>(interFace2ZoneID     ), nIFace*sizeof(int));
        vfile->read(reinterpret_cast<char*>(interFace2InterFaceID), nIFace*sizeof(int));
        vfile->read(reinterpret_cast<char*>(interFace2BoundaryFace), nIFace*sizeof(int));
    }
    SpecifyRightCellofBC();
    PrintToWindow("\n");
}

void UnstructGrid::WriteGrid(VirtualFile *vfile)
{
    int nTotalNode = this->GetNTotalNode();
    int nTotalFace = this->GetNTotalFace();
    int nTotalCell = this->GetNTotalCell();

    //! Write the number of nodes, faces, cells.
    vfile->write(reinterpret_cast<char*>(&nTotalNode),sizeof(int));
    vfile->write(reinterpret_cast<char*>(&nTotalFace),sizeof(int));
    vfile->write(reinterpret_cast<char*>(&nTotalCell),sizeof(int));

    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();

    RotateAboutAxis();

    vfile->write(reinterpret_cast<char*>(x),nTotalNode * sizeof(RDouble));
    vfile->write(reinterpret_cast<char*>(y),nTotalNode * sizeof(RDouble));
    vfile->write(reinterpret_cast<char*>(z),nTotalNode * sizeof(RDouble));

    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();

    //! Write the number of nodes per face.
    vfile->write(reinterpret_cast<char*>(node_number_of_each_face),nTotalFace * sizeof(int));

    int nsum = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        nsum += node_number_of_each_face[iFace];
    }

    //! Write the nodes of each face.
    int *face2node = this->GetFace2Node();
    vfile->write(reinterpret_cast<char*>(face2node), nsum*sizeof(int));

    int *leftCellOfFace  = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    vfile->write(reinterpret_cast<char*>(leftCellOfFace), nTotalFace*sizeof(int));
    vfile->write(reinterpret_cast<char*>(rightCellOfFace), nTotalFace*sizeof(int));

    int nBoundFace = this->GetNBoundFace();
    vfile->write(reinterpret_cast<char*>(&nBoundFace),sizeof(int));

    //! Write the boundary condition.
    UnstructBCSet **bcr = this->GetBCRecord();
    int *tmpbctype = new int [nBoundFace];

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        tmpbctype[iFace] = bcr[iFace]->GetKey();
    }
    vfile->write(reinterpret_cast<char*>(tmpbctype), nBoundFace*sizeof(int));
    delete [] tmpbctype;    tmpbctype = nullptr;

    int nIFace = 0;
    InterfaceInfo *interfaceInfo = this->GetInterfaceInfo();
    if (interfaceInfo)
    {
        nIFace = interfaceInfo->GetNIFace();
    }
    vfile->write(reinterpret_cast<char*>(&nIFace),sizeof(int));

    if (interfaceInfo)
    {
        int *interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
        //int *interFace2CellID       = interfaceInfo->GetInterFace2CellID();
        int *interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();
        int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();
        vfile->write(reinterpret_cast<char*>(interFace2ZoneID     ), nIFace*sizeof(int));
        //vfile->write(reinterpret_cast<char*>(interFace2CellID     ), nIFace*sizeof(int));
        vfile->write(reinterpret_cast<char*>(interFace2InterFaceID), nIFace*sizeof(int));
        vfile->write(reinterpret_cast<char*>(interFace2BoundaryFace), nIFace*sizeof(int));
    }
}

void UnstructGrid::ComputeMetrics(ActionKey *actkey)
{
    AllocateMetrics(actkey);
    if (GetDim() == 2)
    {
        ComputeFaceNormal2D();
        ComputeFaceCenter2D(actkey);
        ComputeCellCenterVol2D(actkey);
    }
    else
    {
        ComputeFaceCenterAreaNormal3DNew2(actkey);
        ComputeCellCenterVol3DNew2(actkey);
    }
}

void UnstructGrid::ComputeFaceNormal3DNew()
{
    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();
    int *face2node = this->GetFace2Node();
    int nTotalFace = this->GetNTotalFace();
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();

    RDouble *xfn  = this->GetFaceNormalX();
    RDouble *yfn  = this->GetFaceNormalY();
    RDouble *zfn  = this->GetFaceNormalZ();
    RDouble *area = this->GetFaceArea();

    int p1,p2,pt;
    RDouble xct,yct,zct;
    RDouble dx1,dy1,dz1,dx2,dy2,dz2;
    int nodepos = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        xct = 0.0;
        yct = 0.0;
        zct = 0.0;

        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            pt = face2node[ nodepos + j ];
            xct += x[pt];
            yct += y[pt];
            zct += z[pt];
        }
        xct /= node_number_of_each_face[iFace];
        yct /= node_number_of_each_face[iFace];
        zct /= node_number_of_each_face[iFace];

        xfn[iFace] = 0.0;
        yfn[iFace] = 0.0;
        zfn[iFace] = 0.0;

        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            p1 = face2node[ nodepos + j ];
            p2 = face2node[ nodepos + (j + 1) % node_number_of_each_face[iFace] ];

            dx1 = xct - x[p1];
            dy1 = yct - y[p1];
            dz1 = zct - z[p1];

            dx2 = xct - x[p2];
            dy2 = yct - y[p2];
            dz2 = zct - z[p2];

            xfn[iFace] += dy1 * dz2 - dy2 * dz1;
            yfn[iFace] += dz1 * dx2 - dz2 * dx1;
            zfn[iFace] += dx1 * dy2 - dx2 * dy1;
        }

        xfn[iFace] *= half;
        yfn[iFace] *= half;
        zfn[iFace] *= half;

        area[iFace] = sqrt(xfn[iFace] * xfn[iFace] + yfn[iFace] * yfn[iFace] + zfn[iFace] * zfn[iFace]);
        nodepos += node_number_of_each_face[iFace];
    }

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        RDouble oarea = 1.0 / (area[iFace] + SMALL);
        xfn[iFace] *= oarea;
        yfn[iFace] *= oarea;
        zfn[iFace] *= oarea;
    }
}

GridManager *UnstructGrid::GetGridManager()
{
    if (!gridManager)
    {
        gridManager = new GridManager(this);
    }
    return gridManager;
}

void UnstructGrid::ComputeGridFaceVelocity()
{
    GridManager *gridManager = this->GetGridManager();
    gridManager->ComputeGridFaceVelocity();
}

void UnstructGrid::BackUpOldGrid()
{
    GridManager *gridManager = this->GetGridManager();
    gridManager->BackUpOldGrid();
}

void UnstructGrid::ComputeFaceNormal3D()
{
    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();
    int *face2node = this->GetFace2Node();
    int nTotalFace = this->GetNTotalFace();
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();

    RDouble *xfn  = this->GetFaceNormalX();
    RDouble *yfn  = this->GetFaceNormalY();
    RDouble *zfn  = this->GetFaceNormalZ();
    RDouble *area = this->GetFaceArea();

    int p1, p2;
    RDouble dx1, dy1, dz1, dx2, dy2, dz2;
    int nodepos = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        xfn[iFace] = 0.0;
        yfn[iFace] = 0.0;
        zfn[iFace] = 0.0;
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            int index1 = j;
            int index2 = (j + 1) % node_number_of_each_face[iFace];
            p1 = face2node[nodepos + index1];
            p2 = face2node[nodepos + index2];

            dx1 = x[p1];
            dy1 = y[p1];
            dz1 = z[p1];

            dx2 = x[p2];
            dy2 = y[p2];
            dz2 = z[p2];

            xfn[iFace] += half * (dy1 * dz2 - dy2 * dz1);
            yfn[iFace] += half * (dz1 * dx2 - dz2 * dx1);
            zfn[iFace] += half * (dx1 * dy2 - dx2 * dy1);
        }
        area[iFace] = sqrt(xfn[iFace] * xfn[iFace] + yfn[iFace] * yfn[iFace] + zfn[iFace] * zfn[iFace]);
        nodepos += node_number_of_each_face[iFace];
    }

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        RDouble oarea = 1.0 / (area[iFace] + SMALL);
        xfn[iFace] *= oarea;
        yfn[iFace] *= oarea;
        zfn[iFace] *= oarea;
    }
}
/*void UnstructGrid::ComputeFaceTangen()
{
    // 获取网格基础数据
    RDouble* x = this->GetX();
    RDouble* y = this->GetY();
    RDouble* z = this->GetZ();
    int* face2node = this->GetFace2Node();  // 面-顶点索引映射
    int nTotalFace = this->GetNTotalFace();
    int* node_number_of_each_face = this->GetNodeNumberOfEachFace();  // 每个面的顶点数

    // 初始化面向量数组（假设已分配内存）
    RDouble* xft = this->GetFaceTangenX();  // 自定义获取函数
    RDouble* yft = this->GetFaceTangenY();
    RDouble* zft = this->GetFaceTangenZ();

    int nodepos = 0;  // 顶点索引在face2node中的起始位置
    for (int iFace = 0; iFace < nTotalFace; ++iFace)
    {
        // 获取当前面的顶点数量
        int nNodes = node_number_of_each_face[iFace];
        if (nNodes < 2)  // 至少需要2个顶点才能构成向量
        {
            xft[iFace] = 0.0;
            yft[iFace] = 0.0;
            zft[iFace] = 0.0;
            nodepos += nNodes;
            continue;
        }

        // 取面的第1个和第2个相邻顶点（索引0和1）
        int p1 = face2node[nodepos + 0];  // 第一个顶点索引
        int p2 = face2node[nodepos + 1];  // 第二个顶点索引

        // 计算面向量：P2 - P1
        xft[iFace] = x[p2] - x[p1];
        yft[iFace] = y[p2] - y[p1];
        zft[iFace] = z[p2] - z[p1];

        // 单位化面向量（如果需要单位向量）
        RDouble len = sqrt(xft[iFace] * xft[iFace] + yft[iFace] * yft[iFace] + zft[iFace] * zft[iFace]);
        if (len > SMALL)  // 避免除以0
        {
            xft[iFace] /= len;
            yft[iFace] /= len;
            zft[iFace] /= len;
        }

        // 更新顶点索引位置（移动到下一个面的顶点）
        nodepos += nNodes;
    }
}*/

void UnstructGrid::ComputeFaceCenterAreaNormal3DNew2(ActionKey *actkey)
{
    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();
    int nTotalFace = this->GetNTotalFace();
    int *face2node = this->GetFace2Node();
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();

    RDouble *xfn  = this->GetFaceNormalX();
    RDouble *yfn  = this->GetFaceNormalY();
    RDouble *zfn  = this->GetFaceNormalZ();
    RDouble *xfc  = this->GetFaceCenterX();
    RDouble *yfc  = this->GetFaceCenterY();
    RDouble *zfc  = this->GetFaceCenterZ();
    RDouble *area = this->GetFaceArea();

    //! Compute rotating velocity.
    RDouble *vgn = this->GetFaceNormalVelocity();
    RDouble *vgx = this->GetFaceVelocityX();
    RDouble *vgy = this->GetFaceVelocityY();
    RDouble *vgz = this->GetFaceVelocityZ();

    int iFace, j, p1, p2, index1, index2, nodepos;
    RDouble coef, x00, y00, z00, dx1, dy1, dz1, dx2, dy2, dz2, anx, any, anz, areaTemp, oareaTemp;

    nodepos = 0;
    for (iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        xfc[iFace] = 0.0;
        yfc[iFace] = 0.0;
        zfc[iFace] = 0.0;

        for (j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            index1 = face2node[ nodepos + j ];
            xfc[iFace] += x[index1];
            yfc[iFace] += y[index1];
            zfc[iFace] += z[index1];
        }

        coef = 1.0 / node_number_of_each_face[iFace];

        xfc[iFace] *= coef;
        yfc[iFace] *= coef;
        zfc[iFace] *= coef;

        nodepos += node_number_of_each_face[iFace];
    }

    RDouble *xfct  = new RDouble[nTotalFace];
    RDouble *yfct  = new RDouble[nTotalFace];
    RDouble *zfct  = new RDouble[nTotalFace];
    RDouble *areat = new RDouble[nTotalFace];

    int numberCycle = 4;
    for(int iCycle = 0; iCycle < numberCycle; ++ iCycle)
    {
        nodepos = 0;
        for (iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            xfct[iFace]  = 0.0;
            yfct[iFace]  = 0.0;
            zfct[iFace]  = 0.0;
            areat[iFace] = 0.0;

            for (j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                index1 = j;
                index2 = (j + 1) % node_number_of_each_face[iFace];
                p1 = face2node[ nodepos + index1 ];
                p2 = face2node[ nodepos + index2 ];

                dx1 = x[p1] - xfc[iFace];
                dy1 = y[p1] - yfc[iFace];
                dz1 = z[p1] - zfc[iFace];

                dx2 = x[p2] - xfc[iFace];
                dy2 = y[p2] - yfc[iFace];
                dz2 = z[p2] - zfc[iFace];

                anx = dy1 * dz2 - dy2 * dz1;
                any = dz1 * dx2 - dz2 * dx1;
                anz = dx1 * dy2 - dx2 * dy1;

                areaTemp = sqrt(anx * anx + any * any + anz * anz);

                x00 = x[p1] + x[p2] + xfc[iFace];
                y00 = y[p1] + y[p2] + yfc[iFace];
                z00 = z[p1] + z[p2] + zfc[iFace];

                areat[iFace] += areaTemp;
                xfct[iFace]  += areaTemp * x00;
                yfct[iFace]  += areaTemp * y00;
                zfct[iFace]  += areaTemp * z00;
            }

            if(areat[iFace] > SMALL)
            {
                areaTemp = 1.0/areat[iFace]/3.0;
                xfc[iFace] = xfct[iFace] * areaTemp;
                yfc[iFace] = yfct[iFace] * areaTemp;
                zfc[iFace] = zfct[iFace] * areaTemp;
            }

            nodepos += node_number_of_each_face[iFace];
        }
    }

    nodepos = 0;
    for (iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        xfct[iFace]  = 0.0;
        yfct[iFace]  = 0.0;
        zfct[iFace]  = 0.0;
        areat[iFace] = 0.0;
        xfn[iFace]   = 0.0;
        yfn[iFace]   = 0.0;
        zfn[iFace]   = 0.0;

        for (j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            index1 = j;
            index2 = (j + 1) % node_number_of_each_face[iFace];
            p1 = face2node[ nodepos + index1 ];
            p2 = face2node[ nodepos + index2 ];

            dx1 = x[p1] - xfc[iFace];
            dy1 = y[p1] - yfc[iFace];
            dz1 = z[p1] - zfc[iFace];

            dx2 = x[p2] - xfc[iFace];
            dy2 = y[p2] - yfc[iFace];
            dz2 = z[p2] - zfc[iFace];

            anx = dy1 * dz2 - dy2 * dz1;
            any = dz1 * dx2 - dz2 * dx1;
            anz = dx1 * dy2 - dx2 * dy1;

            areaTemp = sqrt(anx * anx + any * any + anz * anz);

            x00 = x[p1] + x[p2] + xfc[iFace];
            y00 = y[p1] + y[p2] + yfc[iFace];
            z00 = z[p1] + z[p2] + zfc[iFace];

            areat[iFace] += areaTemp;
            xfct[iFace]  += areaTemp * x00;
            yfct[iFace]  += areaTemp * y00;
            zfct[iFace]  += areaTemp * z00;
            xfn[iFace]   += anx;
            yfn[iFace]   += any;
            zfn[iFace]   += anz;
        }

        if(areat[iFace] > SMALL)
        {
            areaTemp = 1.0/areat[iFace]/3.0;
            xfc[iFace] = xfct[iFace] * areaTemp;
            yfc[iFace] = yfct[iFace] * areaTemp;
            zfc[iFace] = zfct[iFace] * areaTemp;
        }

        areaTemp    = sqrt(xfn[iFace]*xfn[iFace]+yfn[iFace]*yfn[iFace]+zfn[iFace]*zfn[iFace]);
        area[iFace] = areaTemp * 0.5;
        oareaTemp   = 1.0 / (areaTemp + SMALL);
        xfn[iFace] *= oareaTemp;
        yfn[iFace] *= oareaTemp;
        zfn[iFace] *= oareaTemp;
        
        nodepos += node_number_of_each_face[iFace];
    }

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");

    if (referenceFrame == ROTATIONAL_FRAME)
    {
        int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");

        //! angular velocity
        RDouble Omega[100];
        GlobalDataBase::GetData("Omega", &Omega, PHDOUBLE, nTurboZone);

        string shroud[100];
        GlobalDataBase::GetData("shroud", &shroud, PHSTRING, nTurboZone);

        RDouble refMach = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble refTemp = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

        RDouble refGama = 1.4;
        RDouble refGasConstant = 287.15;
        RDouble refVel = refMach * sqrt(refGama * refGasConstant * refTemp);

        //! Parallel
        int iTurboZone = this->GetOrdinaryGridIndex();
        //! Serial
        if (iTurboZone == -1)
        {
            iTurboZone = this->GetZoneID();
        }

        Omega[iTurboZone] /= refVel;

        for (iFace = 0; iFace < nTotalFace; ++iFace)
        {
            vgx[iFace] = 0.0;
            vgy[iFace] = 0.0;
            vgz[iFace] = 0.0;

            //! Compute velocity of each face.
            vgy[iFace] = - Omega[iTurboZone] * zfc[iFace];
            vgz[iFace] =   Omega[iTurboZone] * yfc[iFace];

            //! Set velocity of each face.
            vgn[iFace] = vgx[iFace] * xfn[iFace] + vgy[iFace] * yfn[iFace] + vgz[iFace] * zfn[iFace];
        }

        //! rotating speed of shroud is 0.
        for (iFace = 0; iFace < nBoundFace; iFace++)
        {
            UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
            int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
            string bcName = bcRegion->GetBCName();
            int bcType = bcRegion->GetBCType();

            for (int iTurboZone = 0; iTurboZone < nTurboZone; iTurboZone++)
            {
                if (bcName == shroud[iTurboZone])
                {
                    vgx[iFace] = 0;
                    vgy[iFace] = 0;
                    vgz[iFace] = 0;

                    //! Set velocity of each face.
                    vgn[iFace] = 0;
                }
            }
        }
    }

    delete [] xfct;    xfct = nullptr;
    delete [] yfct;    yfct = nullptr;
    delete [] zfct;    zfct = nullptr;
    delete [] areat;    areat = nullptr;
}

void UnstructGrid::ComputeFaceNormal2D()
{
    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    int *face2node = this->GetFace2Node();
    int nTotalFace = this->GetNTotalFace();

    RDouble *xfn  = this->GetFaceNormalX();
    RDouble *yfn  = this->GetFaceNormalY();
    RDouble *zfn  = this->GetFaceNormalZ();
    RDouble *area = this->GetFaceArea();

    RDouble oarea;

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int p1  = face2node[ 2 * iFace     ];
        int p2  = face2node[ 2 * iFace + 1 ];

        xfn[iFace]  = y[p2] - y[p1];
        yfn[iFace]  = x[p1] - x[p2];
        zfn[iFace]  = 0.0;

        area[iFace] = sqrt(xfn[iFace] * xfn[iFace] + yfn[iFace] * yfn[iFace] + zfn[iFace] * zfn[iFace]);
        oarea   = one / (area[iFace] + SMALL);
        xfn[iFace] *= oarea;
        yfn[iFace] *= oarea;
        zfn[iFace] *= oarea;
    }
}

void UnstructGrid::ComputeFaceCenter3D(ActionKey *actkey)
{
    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();
    int *face2node = this->GetFace2Node();
    int nTotalFace = this->GetNTotalFace();
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();

    RDouble *xfc = this->GetFaceCenterX();
    RDouble *yfc = this->GetFaceCenterY();
    RDouble *zfc = this->GetFaceCenterZ();

    RDouble *xfn = this->GetFaceNormalX();
    RDouble *yfn = this->GetFaceNormalY();
    RDouble *zfn = this->GetFaceNormalZ();

    RDouble x0,y0,z0,x00,y00,z00;
    RDouble dx1,dy1,dz1,dx2,dy2,dz2,anx,any,anz,area,sarea,osarea,norm;

    int nodepos = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        x0 = 0.0;
        y0 = 0.0;
        z0 = 0.0;

        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            int index = face2node[ nodepos + j ];
            x0 += x[ index ];
            y0 += y[ index ];
            z0 += z[ index ];
        }

        RDouble coef = one / node_number_of_each_face[iFace];

        x0 *= coef;
        y0 *= coef;
        z0 *= coef;

        xfc[iFace] = 0.0;
        yfc[iFace] = 0.0;
        zfc[iFace] = 0.0;
        sarea  = 0.0;

        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            int index1 = j;
            int index2 = (j + 1) % node_number_of_each_face[iFace];
            int p1 = face2node[ nodepos + index1 ];
            int p2 = face2node[ nodepos + index2 ];

            dx1 = x[p1] - x0;
            dy1 = y[p1] - y0;
            dz1 = z[p1] - z0;

            dx2 = x[p2] - x0;
            dy2 = y[p2] - y0;
            dz2 = z[p2] - z0;

            x00 = third * (x[p1] + x[p2] + x0);
            y00 = third * (y[p1] + y[p2] + y0);
            z00 = third * (z[p1] + z[p2] + z0);

            anx = dy1 * dz2 - dy2 * dz1;
            any = dz1 * dx2 - dz2 * dx1;
            anz = dx1 * dy2 - dx2 * dy1;

            area = half * sqrt(anx * anx + any * any + anz * anz);
            norm = anx * xfn[iFace] + any * yfn[iFace] + anz * zfn[iFace];

            if (norm < 0.0) area  = - area;
            sarea += area;

            xfc[iFace] += area * x00;
            yfc[iFace] += area * y00;
            zfc[iFace] += area * z00;
        }

        osarea = 1.0 / (sarea + SMALL);

        xfc[iFace] *= osarea;
        yfc[iFace] *= osarea;
        zfc[iFace] *= osarea;

        nodepos += node_number_of_each_face[iFace];
    }
}

void UnstructGrid::ComputeFaceCenter3DNew(ActionKey *actkey)
{
    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();
    int *face2node = this->GetFace2Node();
    int nTotalFace = this->GetNTotalFace();
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();

    RDouble *xfc = this->GetFaceCenterX();
    RDouble *yfc = this->GetFaceCenterY();
    RDouble *zfc = this->GetFaceCenterZ();

    RDouble *xfn = this->GetFaceNormalX();
    RDouble *yfn = this->GetFaceNormalY();
    RDouble *zfn = this->GetFaceNormalZ();

    RDouble xct, yct, zct, x00, y00, z00;
    RDouble dx1, dy1, dz1, dx2, dy2, dz2, anx, any, anz, sss, suma, norm;
    RDouble tmp;

    int p1, p2, pt;

    int nodepos = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        xct = 0.0;
        yct = 0.0;
        zct = 0.0;

        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            pt = face2node[ nodepos + j ];
            xct += x[ pt ];
            yct += y[ pt ];
            zct += z[ pt ];
        }

        tmp = 1.0 / node_number_of_each_face[iFace];

        xct *= tmp;
        yct *= tmp;
        zct *= tmp;

        xfc[iFace] = 0.0;
        yfc[iFace] = 0.0;
        zfc[iFace] = 0.0;
        suma   = 0.0;

        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            p1 = face2node[ nodepos + j ];
            p2 = face2node[ nodepos + (j + 1) % node_number_of_each_face[iFace] ];

            dx1 = xct - x[p1];
            dy1 = yct - y[p1];
            dz1 = zct - z[p1];

            dx2 = xct - x[p2];
            dy2 = yct - y[p2];
            dz2 = zct - z[p2];

            x00 = x[p1] + x[p2] + xct;
            y00 = y[p1] + y[p2] + yct;
            z00 = z[p1] + z[p2] + zct;

            anx = dy1 * dz2 - dy2 * dz1;
            any = dz1 * dx2 - dz2 * dx1;
            anz = dx1 * dy2 - dx2 * dy1;

            sss  = sqrt(anx * anx + any * any + anz * anz);
            norm = anx * xfn[iFace] + any * yfn[iFace] + anz * zfn[iFace];

            if (norm < 0.0) sss = - sss;
            suma += sss;

            xfc[iFace] += sss * x00;
            yfc[iFace] += sss * y00;
            zfc[iFace] += sss * z00;
        }

        tmp = third / (suma + SMALL);

        xfc[iFace] *= tmp;
        yfc[iFace] *= tmp;
        zfc[iFace] *= tmp;

        nodepos += node_number_of_each_face[iFace];
    }
}

void UnstructGrid::ComputeFaceCenter2D(ActionKey *actkey)
{
    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();
    int *face2node = this->GetFace2Node();
    int nTotalFace = this->GetNTotalFace();

    RDouble *xfc = this->GetFaceCenterX();
    RDouble *yfc = this->GetFaceCenterY();
    RDouble *zfc = this->GetFaceCenterZ();

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int p1  = face2node[ 2 * iFace     ];
        int p2  = face2node[ 2 * iFace + 1 ];
        xfc[iFace] = half * (x[p1] + x[p2]);
        yfc[iFace] = half * (y[p1] + y[p2]);
        zfc[iFace] = half * (z[p1] + z[p2]);
    }
}

void UnstructGrid::ComputeCellCenterVol3D(ActionKey *actkey)
{
    int cell, le, re, p2, p3;
    int nBoundFace = this->GetNBoundFace();
    int nTotalFace = this->GetNTotalFace();
    int nTotalCell = this->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();
    int *face2node = this->GetFace2Node();
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();
    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();

    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();

    RDouble *xfc = this->GetFaceCenterX();
    RDouble *yfc = this->GetFaceCenterY();
    RDouble *zfc = this->GetFaceCenterZ();

    RDouble *nxs = this->GetFaceNormalX();
    RDouble *nys = this->GetFaceNormalY();
    RDouble *nzs = this->GetFaceNormalZ();
    RDouble *ns  = this->GetFaceArea();

    RDouble *vol = this->GetCellVolume();

    RDouble x21, y21, z21, x31, y31, z31, nx, ny, nz;
    RDouble xc,  yc,  zc,  tmp;

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        xcc[iCell] = 0.0;
        ycc[iCell] = 0.0;
        zcc[iCell] = 0.0;
    }

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        vol[iCell] = 0.0;
    }

    int nodepos = 0;

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        if (iFace < nBoundFace) re = iFace + nTotalCell;

        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            int index1 = j;
            int index2 = (j + 1) % node_number_of_each_face[iFace];
            p2 = face2node[ nodepos + index1 ];
            p3 = face2node[ nodepos + index2 ];

            x21 = x[p2] - xfc[iFace];
            y21 = y[p2] - yfc[iFace];
            z21 = z[p2] - zfc[iFace];

            x31 = x[p3] - xfc[iFace];
            y31 = y[p3] - yfc[iFace];
            z31 = z[p3] - zfc[iFace];

            nx  = y21 * z31 - y31 * z21;
            ny  = z21 * x31 - z31 * x21;
            nz  = x21 * y31 - x31 * y21;

            xc  = xfc[iFace] + x[p2] + x[p3];
            yc  = yfc[iFace] + y[p2] + y[p3];
            zc  = zfc[iFace] + z[p2] + z[p3];

            //! Cell Center and Volume.
            tmp     = nx * xc + ny * yc + nz * zc;
            xc      *= tmp;
            yc      *= tmp;
            zc      *= tmp;
            xcc[le] += xc;
            ycc[le] += yc;
            zcc[le] += zc;
            vol[le] += tmp;

            xcc[re] -= xc;
            ycc[re] -= yc;
            zcc[re] -= zc;
            vol[re] -= tmp;
        }
        nodepos += node_number_of_each_face[iFace];
    }

    std::ostringstream oss;

    //! Don't forget the coefficients.
    cell = 0;
    RDouble minvol = LARGE, maxvol = 0.0;
    int index_minv = 0,index_maxv = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        tmp     = 1.0 / (4.0 * vol[iCell] + SMALL);
        xcc[iCell] *= tmp;
        ycc[iCell] *= tmp;
        zcc[iCell] *= tmp;
        vol[iCell] /= 18.0;
        if (minvol > vol[iCell])
        {
            minvol     = vol[iCell];
            index_minv = iCell;
        }
        if (maxvol < vol[iCell])
        {
            maxvol     = vol[iCell];
            index_maxv = iCell;
        }

        if (vol[iCell] <= 0.0)
        {
            ++ cell;
            oss << "  Warning: negative volume cell index " << iCell << ", (x, y, z) " << xcc[iCell] << " " << ycc[iCell] << " " << zcc[iCell] << " !\n";
        }
    }

    if (cell) oss << cell << " cells have negative vols \n";

    //! For ghost cells.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le         = leftCellOfFace[ iFace ];
        //re         = rightCellOfFace[ iFace ];
        re         = iFace + nTotalCell;

        if (ns[iFace] > SMALL)
        {
           tmp     = 2.0 * ((xcc[le] - xfc[iFace]) * nxs[iFace] + (ycc[le] - yfc[iFace]) * nys[iFace] + (zcc[le] - zfc[iFace]) * nzs[iFace]);
           xcc[re] = xcc[le] - nxs[iFace] * tmp;
           ycc[re] = ycc[le] - nys[iFace] * tmp;
           zcc[re] = zcc[le] - nzs[iFace] * tmp;

        }
        else
        {
            //! Degenerated faces.
           xcc[re] = - xcc[le] + 2.0 * xfc[iFace];
           ycc[re] = - ycc[le] + 2.0 * yfc[iFace];
           zcc[re] = - zcc[le] + 2.0 * zfc[iFace];
       }
       vol[re] = vol[le];
    }
}

void UnstructGrid::ComputeCellCenterVol3DNew2(ActionKey *actkey)
{
    int nBoundFace = this->GetNBoundFace();
    int nTotalFace = this->GetNTotalFace();
    int nTotalCell = this->GetNTotalCell();
    int *leftCellOfFace  = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();
    int *face2node = this->GetFace2Node();
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();
    int *cell2node = this->GetCell2Node();
    int *node_number_of_each_cell = this->GetNodeNumberOfEachCell();
    RDouble *x   = this->GetX();
    RDouble *y   = this->GetY();
    RDouble *z   = this->GetZ();

    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();

    RDouble *xfc = this->GetFaceCenterX();
    RDouble *yfc = this->GetFaceCenterY();
    RDouble *zfc = this->GetFaceCenterZ();

    RDouble *nxs = this->GetFaceNormalX();
    RDouble *nys = this->GetFaceNormalY();
    RDouble *nzs = this->GetFaceNormalZ();
    RDouble *ns  = this->GetFaceArea();

    RDouble *vol = this->GetCellVolume();

    int iFace, iCell, le, re, j, p2, p3, nodepos, index1, index2, cell;
    RDouble x21, y21, z21, x31, y31, z31, x41, y41, z41, nx, ny, nz;
    RDouble xc,  yc,  zc,  tmp;
    std::ostringstream oss;

    nodepos = 0;
    for (iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        xcc[iCell] = 0.0;
        ycc[iCell] = 0.0;
        zcc[iCell] = 0.0;

        for (j = 0; j < node_number_of_each_cell[iCell]; ++ j)
        {
            p2 = cell2node[ nodepos + j ];
            xcc[iCell] += x[p2];
            ycc[iCell] += y[p2];
            zcc[iCell] += z[p2];
        }

        xcc[iCell] /= node_number_of_each_cell[iCell];
        ycc[iCell] /= node_number_of_each_cell[iCell];
        zcc[iCell] /= node_number_of_each_cell[iCell];

        nodepos += node_number_of_each_cell[iCell];
    }
    
    RDouble *xcct = new RDouble[nTotalCell];
    RDouble *ycct = new RDouble[nTotalCell];
    RDouble *zcct = new RDouble[nTotalCell];
    for (iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        xcct[iCell] = 0.0;
        ycct[iCell] = 0.0;
        zcct[iCell] = 0.0;
        vol [iCell] = 0.0;
    }

    nodepos = 0;
    for (iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace [ iFace ];
       
        for (j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            index1 = j;
            index2 = (j + 1) % node_number_of_each_face[iFace];
            p2 = face2node[ nodepos + index1 ];
            p3 = face2node[ nodepos + index2 ];

            x21 = x[p2] - xfc[iFace];
            y21 = y[p2] - yfc[iFace];
            z21 = z[p2] - zfc[iFace];

            x31 = x[p3] - xfc[iFace];
            y31 = y[p3] - yfc[iFace];
            z31 = z[p3] - zfc[iFace];

            nx  = y21 * z31 - y31 * z21;
            ny  = z21 * x31 - z31 * x21;
            nz  = x21 * y31 - x31 * y21;

            x41 = xcc[le] - xfc[iFace];
            y41 = ycc[le] - yfc[iFace];
            z41 = zcc[le] - zfc[iFace];

            xc  = xfc[iFace] + x[p2] + x[p3] + xcc[le];
            yc  = yfc[iFace] + y[p2] + y[p3] + ycc[le];
            zc  = zfc[iFace] + z[p2] + z[p3] + zcc[le];

            tmp       = nx * x41 + ny * y41 + nz * z41;
            tmp       = -tmp;
            xc       *= tmp;
            yc       *= tmp;
            zc       *= tmp;
            xcct[le] += xc;
            ycct[le] += yc;
            zcct[le] += zc;
            vol[le]  += tmp;
        }
        nodepos += node_number_of_each_face[iFace];
    }
    for (iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [ iFace ];
        re = rightCellOfFace[ iFace ];
       
        for (j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            index1 = j;
            index2 = (j + 1) % node_number_of_each_face[iFace];
            p2 = face2node[ nodepos + index1 ];
            p3 = face2node[ nodepos + index2 ];

            x21 = x[p2] - xfc[iFace];
            y21 = y[p2] - yfc[iFace];
            z21 = z[p2] - zfc[iFace];

            x31 = x[p3] - xfc[iFace];
            y31 = y[p3] - yfc[iFace];
            z31 = z[p3] - zfc[iFace];

            nx  = y21 * z31 - y31 * z21;
            ny  = z21 * x31 - z31 * x21;
            nz  = x21 * y31 - x31 * y21;

            x41 = xcc[le] - xfc[iFace];
            y41 = ycc[le] - yfc[iFace];
            z41 = zcc[le] - zfc[iFace];

            xc  = xfc[iFace] + x[p2] + x[p3] + xcc[le];
            yc  = yfc[iFace] + y[p2] + y[p3] + ycc[le];
            zc  = zfc[iFace] + z[p2] + z[p3] + zcc[le];

            tmp       = nx * x41 + ny * y41 + nz * z41;
            tmp       = -tmp;
            xc       *= tmp;
            yc       *= tmp;
            zc       *= tmp;
            xcct[le] += xc;
            ycct[le] += yc;
            zcct[le] += zc;
            vol[le]  += tmp;

            x41 = xcc[re] - xfc[iFace];
            y41 = ycc[re] - yfc[iFace];
            z41 = zcc[re] - zfc[iFace];

            xc  = xfc[iFace] + x[p2] + x[p3] + xcc[re];
            yc  = yfc[iFace] + y[p2] + y[p3] + ycc[re];
            zc  = zfc[iFace] + z[p2] + z[p3] + zcc[re];

            tmp       = nx * x41 + ny * y41 + nz * z41;
            xc       *= tmp;
            yc       *= tmp;
            zc       *= tmp;
            xcct[re] += xc;
            ycct[re] += yc;
            zcct[re] += zc;
            vol[re]  += tmp;
        }
        nodepos += node_number_of_each_face[iFace];
    }
    
    cell = 0;
    RDouble minvol = LARGE, maxvol = 0.0;
    int index_minv = 0,index_maxv = 0;
    for (iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        tmp          = 1.0 / (4.0 * vol[iCell] + SMALL);
        xcct[iCell] *= tmp;
        ycct[iCell] *= tmp;
        zcct[iCell] *= tmp;
        vol[iCell]  /= 6.0;
        if (minvol > vol[iCell])
        {
            minvol     = vol[iCell];
            index_minv = iCell;
        }
        if (maxvol < vol[iCell])
        {
            maxvol     = vol[iCell];
            index_maxv = iCell;
        }

        if (vol[iCell] <= 0.0)
        {
            ++ cell;
            oss << "  Warning: negative volume cell index " << iCell << ", (x, y, z) " << xcc[iCell] << " " << ycc[iCell] << " " << zcc[iCell] << " !\n";
        }
    }

    if (cell) oss << cell << " cells have negative vols \n";
    
    int *mark = new int[nTotalCell];
    for (iCell = 0; iCell < nTotalCell; ++ iCell) mark[iCell] = 1;
    
    for (iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace[iFace];

        x21 = xcct[le] - xfc[iFace];
        y21 = ycct[le] - yfc[iFace];
        z21 = zcct[le] - zfc[iFace];
        tmp = x21 * nxs[iFace] + y21 * nys[iFace] + z21 * nzs[iFace];
        if(tmp > 0.0) mark[le] = 0;
    }
    for (iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace[iFace];
        re = rightCellOfFace[iFace];

        x21 = xcct[le] - xfc[iFace];
        y21 = ycct[le] - yfc[iFace];
        z21 = zcct[le] - zfc[iFace];  
        tmp = x21 * nxs[iFace] + y21 * nys[iFace] + z21 * nzs[iFace];
        if(tmp > 0.0) mark[le] = 0;

        x21 = xcct[re] - xfc[iFace];
        y21 = ycct[re] - yfc[iFace];
        z21 = zcct[re] - zfc[iFace];  
        tmp = x21 * nxs[iFace] + y21 * nys[iFace] + z21 * nzs[iFace];
        if(tmp < 0.0) mark[re] = 0;
    }

    for (iCell = 0; iCell < nTotalCell; ++iCell)
    {
        if (mark[iCell])
        {
            xcc[iCell] = xcct[iCell];
            ycc[iCell] = ycct[iCell];
            zcc[iCell] = zcct[iCell];
        }
    }

    delete [] mark;    mark = nullptr;
    delete [] xcct;    xcct = nullptr;
    delete [] ycct;    ycct = nullptr;
    delete [] zcct;    zcct = nullptr;

    //! For ghost cells.
    for (iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace[ iFace ];
        //re = rightCellOfFace[ iFace ];
        re = iFace + nTotalCell;

        if (ns[iFace] > SMALL)
        {
           tmp     = 2.0 * ((xcc[le] - xfc[iFace]) * nxs[iFace] + (ycc[le] - yfc[iFace]) * nys[iFace] + (zcc[le] - zfc[iFace]) * nzs[iFace]);
           xcc[re] = xcc[le] - nxs[iFace] * tmp;
           ycc[re] = ycc[le] - nys[iFace] * tmp;
           zcc[re] = zcc[le] - nzs[iFace] * tmp;

        }
        else
        {
            //! Degenerated faces.
           xcc[re] = - xcc[le] + 2.0 * xfc[iFace];
           ycc[re] = - ycc[le] + 2.0 * yfc[iFace];
           zcc[re] = - zcc[le] + 2.0 * zfc[iFace];
       }
       vol[re] = vol[le];
    }
}

void UnstructGrid::ComputeCellCenterVol3DNew(ActionKey *actkey)
{
    int nBoundFace = this->GetNBoundFace();
    int nTotalFace = this->GetNTotalFace();
    int nTotalCell = this->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();
    int *face2node = this->GetFace2Node();
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();
    RDouble *x    = this->GetX();
    RDouble *y    = this->GetY();
    RDouble *z    = this->GetZ();

    RDouble *xcc  = this->GetCellCenterX();
    RDouble *ycc  = this->GetCellCenterY();
    RDouble *zcc  = this->GetCellCenterZ();

    RDouble *xfc  = this->GetFaceCenterX();
    RDouble *yfc  = this->GetFaceCenterY();
    RDouble *zfc  = this->GetFaceCenterZ();

    RDouble *xfn  = this->GetFaceNormalX();
    RDouble *yfn  = this->GetFaceNormalY();
    RDouble *zfn  = this->GetFaceNormalZ();
    RDouble *area = this->GetFaceArea();
    RDouble *vol  = this->GetCellVolume();

    RDouble x21, y21, z21, x31, y31, z31, nx, ny, nz;
    RDouble vvol, vvoll,vvolr,tmp;
    RDouble xc,yc,zc,xcl,ycl,zcl,xcr,ycr,zcr;
    RDouble hxcl,hycl,hzcl,hxcr,hycr,hzcr;

    int cell, le, re, p2, p3;
    RDouble c34 = 3.0 / 4.0;

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        xcc[iCell] = 0.0;
        ycc[iCell] = 0.0;
        zcc[iCell] = 0.0;
    }

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        vol[iCell] = 0.0;
    }

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        le     = leftCellOfFace [ iFace ];
        re     = rightCellOfFace[ iFace ];
        if (iFace < nBoundFace) re = iFace + nTotalCell;

        xc = xfc[iFace];
        yc = yfc[iFace];
        zc = zfc[iFace];

        vvol = third * (xfn[iFace] * xc + yfn[iFace] * yc + zfn[iFace] * zc) * area[iFace];
        xc *= c34 * vvol;
        yc *= c34 * vvol;
        zc *= c34 * vvol;

        xcc[le] += xc;
        ycc[le] += yc;
        zcc[le] += zc;
        vol[le] += vvol;

        xcc[re] -= xc;
        ycc[re] -= yc;
        zcc[re] -= zc;
        vol[re] -= vvol;
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        tmp = 1.0 / (vol[iCell] + SMALL);
        xcc[iCell] *= tmp;
        ycc[iCell] *= tmp;
        zcc[iCell] *= tmp;
    }

    RDouble *txc = new RDouble[nTotal];
    RDouble *tyc = new RDouble[nTotal];
    RDouble *tzc = new RDouble[nTotal];

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        txc[iCell] = xcc[iCell];
        tyc[iCell] = ycc[iCell];
        tzc[iCell] = zcc[iCell];
    }

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        xcc[iCell] = 0.0;
        ycc[iCell] = 0.0;
        zcc[iCell] = 0.0;
    }

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        vol[iCell] = 0.0;
    }

    int nodepos = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        le     = leftCellOfFace [ iFace ];
        re     = rightCellOfFace[ iFace ];
        if (iFace < nBoundFace) re = iFace + nTotalCell;

        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            p2 = face2node[ nodepos + j ];
            p3 = face2node[ nodepos + (j + 1) % node_number_of_each_face[iFace] ];

            x21 = x[p2] - xfc[iFace];
            y21 = y[p2] - yfc[iFace];
            z21 = z[p2] - zfc[iFace];

            x31 = x[p3] - xfc[iFace];
            y31 = y[p3] - yfc[iFace];
            z31 = z[p3] - zfc[iFace];

            nx  = half * (y21 * z31 - y31 * z21);
            ny  = half * (z21 * x31 - z31 * x21);
            nz  = half * (x21 * y31 - x31 * y21);

            xc  = third * (xfc[iFace] + x[p2] + x[p3]);
            yc  = third * (yfc[iFace] + y[p2] + y[p3]);
            zc  = third * (zfc[iFace] + z[p2] + z[p3]);

            hxcl = xc - txc[le];
            hycl = yc - tyc[le];
            hzcl = zc - tzc[le];

            hxcr = xc - txc[re];
            hycr = yc - tyc[re];
            hzcr = zc - tzc[re];

            //! Cell Center and Volume.
            vvoll = third * (nx * hxcl + ny * hycl + nz * hzcl);
            vvolr = third * (nx * hxcr + ny * hycr + nz * hzcr);

            xcl = xc - fourth * hxcl;
            ycl = yc - fourth * hycl;
            zcl = zc - fourth * hzcl;
            xcl *= vvoll;
            ycl *= vvoll;
            zcl *= vvoll;

            xcc[le] += xcl;
            ycc[le] += ycl;
            zcc[le] += zcl;
            vol[le] += vvoll;

            xcr = xc - fourth * hxcr;
            ycr = yc - fourth * hycr;
            zcr = zc - fourth * hzcr;
            xcr *= vvolr;
            ycr *= vvolr;
            zcr *= vvolr;

            xcc[re] -= xcr;
            ycc[re] -= ycr;
            zcc[re] -= zcr;
            vol[re] -= vvolr;
        }
        nodepos += node_number_of_each_face[iFace];
    }

    delete [] txc;    txc = nullptr;
    delete [] tyc;    tyc = nullptr;
    delete [] tzc;    tzc = nullptr;

    std::ostringstream oss;

    //! Don't forget the coefficients.
    cell = 0;
    RDouble minvol = LARGE, maxvol = 0.0;
    int index_minv = 0,index_maxv = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        tmp = 1.0 / (vol[iCell] + SMALL);
        xcc[iCell] *= tmp;
        ycc[iCell] *= tmp;
        zcc[iCell] *= tmp;

        if (minvol > vol[iCell])
        {
            minvol = vol[iCell];
            index_minv = iCell;
        }
        if (maxvol < vol[iCell])
        {
            maxvol = vol[iCell];
            index_maxv = iCell;
        }

        //minvol = MIN(minvol, vol[iCell]);
        //maxvol = MAX(maxvol, vol[iCell]);
        if (vol[iCell] <= 0.0)
        {
            vol[iCell] = - vol[iCell];
            ++ cell;
            oss << "  Warning: negative volume cell index " << iCell << ", (x, y, z) " << xcc[iCell] << " " << ycc[iCell] << " " << zcc[iCell] << " !\n";
        }
    }

    if (cell) oss << cell << " cells have negative vols \n";

    //! For ghost cells.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace[ iFace ];
        re = iFace + nTotalCell;

        if (area[iFace] > SMALL)
        {
           tmp = 2.0 * ((xcc[le] - xfc[iFace]) * xfn[iFace] + (ycc[le] - yfc[iFace]) * yfn[iFace] + (zcc[le] - zfc[iFace]) * zfn[iFace]);
           xcc[re] = xcc[le] - xfn[iFace] * tmp;
           ycc[re] = ycc[le] - yfn[iFace] * tmp;
           zcc[re] = zcc[le] - zfn[iFace] * tmp;

        }
        else
        {
            //! Degenerated faces.
           xcc[re] = - xcc[le] + 2.0 * xfc[iFace];
           ycc[re] = - ycc[le] + 2.0 * yfc[iFace];
           zcc[re] = - zcc[le] + 2.0 * zfc[iFace];
       }
       vol[re] = vol[le];
    }
}

void UnstructGrid::ComputeCellCenterVol2D(ActionKey *actkey)
{
    int nBoundFace = this->GetNBoundFace();
    int nTotalFace = this->GetNTotalFace();
    int nTotalCell = this->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();

    RDouble *xfc = this->GetFaceCenterX();
    RDouble *yfc = this->GetFaceCenterY();
    RDouble *zfc = this->GetFaceCenterZ();

    RDouble *nxs = this->GetFaceNormalX();
    RDouble *nys = this->GetFaceNormalY();
    RDouble *nzs = this->GetFaceNormalZ();
    RDouble *ns  = this->GetFaceArea();
    RDouble *vol = this->GetCellVolume();

    int negativecell, le, re;
    RDouble dot, tmp;

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        xcc[iCell] = 0.0;
        ycc[iCell] = 0.0;
        zcc[iCell] = 0.0;
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        vol[iCell] = 0.0;
    }

    //! Now the cell center and volume.
    //! For each face, boundary faces first.
    //! For boundary faces.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le       = leftCellOfFace[ iFace ];
        dot      = (xfc[iFace] * nxs[iFace] + yfc[iFace] * nys[iFace] + zfc[iFace] * nzs[iFace]) * ns[iFace];
        xcc[le] += xfc[iFace] * dot;
        ycc[le] += yfc[iFace] * dot;
        zcc[le] += zfc[iFace] * dot;
        vol[le] += dot;
    }
    //! For interior cell faces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le       = leftCellOfFace [ iFace ];
        re       = rightCellOfFace[ iFace ];
        dot      = (xfc[iFace] * nxs[iFace] + yfc[iFace] * nys[iFace] + zfc[iFace] * nzs[iFace]) * ns[iFace];
        xcc[le] += xfc[iFace] * dot;
        ycc[le] += yfc[iFace] * dot;
        zcc[le] += zfc[iFace] * dot;
        vol[le] += dot;
        xcc[re] -= xfc[iFace] * dot;
        ycc[re] -= yfc[iFace] * dot;
        zcc[re] -= zfc[iFace] * dot;
        vol[re] -= dot;
    }

    //! Don't forget the coefficients.
    negativecell = 0;
    RDouble minvol = LARGE, maxvol = 0.0;
    int index_minv = 0, index_maxv = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        tmp = 1.0 / (1.5 * vol[iCell] + SMALL);
        xcc[iCell] *= tmp;
        ycc[iCell] *= tmp;
        zcc[iCell] *= tmp;
        vol[iCell] *= half;

        if (minvol > vol[iCell])
        {
            minvol = vol[iCell];
            index_minv = iCell;
        }
        if (maxvol < vol[iCell])
        {
            maxvol = vol[iCell];
            index_maxv = iCell;
        }
        //minvol = MIN(minvol, vol[iCell]);
        //maxvol = MAX(maxvol, vol[iCell]);
        if (vol[iCell] <= 0.0)
        {
            vol[iCell] = - vol[iCell];
            ++ negativecell;
            std::ostringstream oss;
            oss << "  Warning: negative volume cell index " << iCell << ", (x, y) " << xcc[iCell] << " " << ycc[iCell] << " !\n";
            if (PHMPI::GetCurrentProcessorID() == PHMPI::GetServerProcessorID())
            {
                cout << oss.str();
            }
        }
    }

    //if (cell) oss << cell << " cells have negative vols \n";
    if (negativecell == nTotalCell)
    {
        //oss << cell << " cells have negative vols \n";
        TK_Exit::ExceptionExit("ERROR: cells have negative vols, maybe the i j k normal direction has set error! ");
    }

    //! For ghost cells.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace[ iFace ];
        re = iFace + nTotalCell;
        if (ns[iFace] > SMALL)
        {
            tmp = 2.0 * ((xcc[le] - xfc[iFace]) * nxs[iFace] + (ycc[le] - yfc[iFace]) * nys[iFace] + (zcc[le] - zfc[iFace]) * nzs[iFace]);
            xcc[iFace+nTotalCell] = xcc[le] - nxs[iFace] * tmp;
            ycc[iFace+nTotalCell] = ycc[le] - nys[iFace] * tmp;
            zcc[iFace+nTotalCell] = zcc[le] - nzs[iFace] * tmp;
        }
        else
        {
            //! Degenerated faces.
            xcc[iFace+nTotalCell] = - xcc[le] + 2.0 * xfc[iFace];
            ycc[iFace+nTotalCell] = - ycc[le] + 2.0 * yfc[iFace];
            zcc[iFace+nTotalCell] = - zcc[le] + 2.0 * zfc[iFace];
        }
        vol[re] = vol[le];
    }

//  oss << "Zone Index : " << this->GetZoneID() << ", min and max volumes are " 
//      << index_minv << ": " << minvol << " " << index_maxv << ": " << maxvol << "\n";
//
//  StreamToActionKey(actkey, oss);
}

void UnstructGrid::ComputeLargestLocalGridLength()
{
    int nTotalCell = this->GetNTotalCell();
    int nBoundFace = this->GetNBoundFace();
    int nTotalFace = this->GetNTotalFace();

    int *leftCellOfFace  = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();

    RDouble *largestLocalGridLength = new RDouble [ nTotalCell ];
    PHSPACE::SetField(largestLocalGridLength, -LARGE, nTotalCell);

    int le, re;
    RDouble dist;
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace[iFace];
        re = rightCellOfFace[iFace];

        dist = DISTANCE(xcc[le] - xcc[re], ycc[le] - ycc[re], zcc[le] - zcc[re]);

        largestLocalGridLength[le] = MAX(largestLocalGridLength[le], dist);
        largestLocalGridLength[re] = MAX(largestLocalGridLength[re], dist);
    }

    this->cellMetrics->SetLargestLocalGridLength(largestLocalGridLength);
}

void UnstructGrid::ComputeSubgridLength()
{
    int nTotalCell = this->GetNTotalCell();

    RDouble *walldist = this->GetWallDist();
    RDouble *largestLocalGridLength = this->GetLargestLocalGridLength();

    RDouble *subgridLength = new RDouble [ nTotalCell ];

    const RDouble cw = 0.15;
    RDouble wallDistance, largestLocalGridSpacing, subgridLengthScale;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell) 
    {
        wallDistance = walldist[iCell];
        largestLocalGridSpacing = largestLocalGridLength[iCell];

        //! The change ratio of the normal distance is ignored bellow!!! Re-consider it future, Please!
        subgridLengthScale = MIN( cw * MAX(wallDistance, largestLocalGridSpacing),  largestLocalGridSpacing );

        subgridLength[ iCell ] = subgridLengthScale;
    }

    this->cellMetrics->SetSubgridLength(subgridLength);
}

void UnstructGrid::DXDYDZ(RDouble *f, int ie, RDouble &dfdx, RDouble &dfdy, RDouble &dfdz)
{
    string gradientName;
    GetData("gradientName", &gradientName, PHSTRING, 1);

    if (gradientName == "lsq")
    {
        DXDYDZ_LSQ(f, ie, dfdx, dfdy, dfdz);
    }
    else
    {
        DXDYDZ_GG(f, ie, dfdx, dfdy, dfdz);
    }
}

void UnstructGrid::DXDYDZ_GG(RDouble *f, int ie, RDouble &dfdx, RDouble &dfdy, RDouble &dfdz)
{
    int **cell2face  = this->GetCell2Face();
    int *face_number_of_each_cell = this->GetFaceNumberOfEachCell();

    RDouble *xfn = this->GetFaceNormalX();
    RDouble *yfn = this->GetFaceNormalY();
    RDouble *zfn = this->GetFaceNormalZ();

    RDouble *vol = this->GetCellVolume();
    RDouble *area = this->GetFaceArea();

    int *leftCellOfFace  = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();
    RDouble nx,ny,nz;
    RDouble fi;
    int face,le,re;

    dfdx = 0.0;
    dfdy = 0.0;
    dfdz = 0.0;

    for (int iFace = 0; iFace < face_number_of_each_cell[ie]; ++ iFace)
    {
        face = cell2face[ie][iFace];
        le   = leftCellOfFace[ face ];
        re   = rightCellOfFace[ face ];
        fi   =  half * (f[le] + f[re]);

        nx   = xfn[face] * area[face];
        ny   = yfn[face] * area[face];
        nz   = zfn[face] * area[face];

        if (le == ie)
        {
            dfdx += fi * nx;
            dfdy += fi * ny;
            dfdz += fi * nz;
        }
        else
        {
            dfdx -= fi * nx;
            dfdy -= fi * ny;
            dfdz -= fi * nz;
        }
    }

    RDouble ovol = 1.0 / vol[ie];

    dfdx *= ovol;
    dfdy *= ovol;
    dfdz *= ovol;
}

void UnstructGrid::DXDYDZ_LSQ(RDouble *f, int ie, RDouble &dfdx, RDouble &dfdy, RDouble &dfdz)
{
    int **cell2face  = this->GetCell2Face();
    int *face_number_of_each_cell = this->GetFaceNumberOfEachCell();

    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();

    RDouble dx,dy,dz;
    RDouble dq,tmpx,tmpy,tmpz;
    int face,le,re;

    dfdx = 0.0;
    dfdy = 0.0;
    dfdz = 0.0;

    RDouble *iwt = this->GetLeastSquareIWT();
    RDouble *ixx = this->GetLeastSquareIXX();
    RDouble *iyy = this->GetLeastSquareIYY();
    RDouble *izz = this->GetLeastSquareIZZ();
    RDouble *ixy = this->GetLeastSquareIXY();
    RDouble *ixz = this->GetLeastSquareIXZ();
    RDouble *iyz = this->GetLeastSquareIYZ();
    char * fMark  = this->GetFaceMark();

    for (int iFace = 0; iFace < face_number_of_each_cell[ie]; ++ iFace)
    {
        face = cell2face[ie][iFace];
        le   = leftCellOfFace[ face ];
        re   = rightCellOfFace[ face ];

        dx  = xcc[re] - xcc[le];
        dy  = ycc[re] - ycc[le];
        dz  = zcc[re] - zcc[le];
 
        dq  = iwt[face] * (f[re] - f[le]) * fMark[face];

        tmpx = dq * dx;
        tmpy = dq * dy;
        tmpz = dq * dz;

        dfdx += tmpx;
        dfdy += tmpy;
        dfdz += tmpz;
    }

    tmpx = dfdx;
    tmpy = dfdy;
    tmpz = dfdz;

    dfdx = ixx[ie] * tmpx + ixy[ie] * tmpy + ixz[ie] * tmpz;
    dfdy = ixy[ie] * tmpx + iyy[ie] * tmpy + iyz[ie] * tmpz;
    dfdz = ixz[ie] * tmpx + iyz[ie] * tmpy + izz[ie] * tmpz;
}

void UnstructGrid::DXDYDZ(RDouble *f, int ie, int ** cell2face, int *face_number_of_each_cell, RDouble &dfdx, RDouble &dfdy, RDouble &dfdz)
{
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();
    RDouble nx,ny,nz;
    RDouble fi;
    int face,le,re;

    RDouble *xfn = this->GetFaceNormalX();
    RDouble *yfn = this->GetFaceNormalY();
    RDouble *zfn = this->GetFaceNormalZ();

    RDouble *vol = this->GetCellVolume();
    RDouble *area = this->GetFaceArea();

    dfdx = 0.0;
    dfdy = 0.0;
    dfdz = 0.0;

    for (int iFace = 0; iFace < face_number_of_each_cell[ie]; ++ iFace)
    {
        face = cell2face[ie][iFace];
        le   = leftCellOfFace[ face ];
        re   = rightCellOfFace[ face ];
        fi   =  half * (f[le] + f[re]);

        nx   = xfn[face] * area[face];
        ny   = yfn[face] * area[face];
        nz   = zfn[face] * area[face];

        if (le == ie)
        {
            dfdx += fi * nx;
            dfdy += fi * ny;
            dfdz += fi * nz;
        }
        else
        {
            dfdx -= fi * nx;
            dfdy -= fi * ny;
            dfdz -= fi * nz;
        }
    }

    RDouble ovol = 1.0 / vol[ie];

    dfdx *= ovol;
    dfdy *= ovol;
    dfdz *= ovol;
}

void UnstructGrid::DXDYDZ_Face(RDouble *f, int iFace, int nodepos, RDouble fl, RDouble fr, RDouble &dfdx, RDouble &dfdy, RDouble &dfdz)
{
    int *leftCellOfFace = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();
    int *face2node = this->GetFace2Node();
    int le, re;
    le = leftCellOfFace[iFace];
    re = rightCellOfFace[iFace];

    RDouble *xcc = this->GetCellCenterX();
    RDouble *ycc = this->GetCellCenterY();
    RDouble *zcc = this->GetCellCenterZ();

    RDouble nx, ny, nz, vvol, tmp;
    RDouble dx1, dy1, dz1, dx2, dy2, dz2;
    RDouble xc, yc, zc;

    RDouble fc;

    dfdx = 0.0;
    dfdy = 0.0;
    dfdz = 0.0;
    vvol = 0.0;

    for (int i = 0; i < node_number_of_each_face[iFace]; ++ i)
    {
        int p1 = face2node[ nodepos + i ];
        int p2 = face2node[ nodepos + (i + 1)%node_number_of_each_face[iFace] ];

        //! Left.
        dx1 = x[p2] - x[p1];
        dy1 = y[p2] - y[p1];
        dz1 = z[p2] - z[p1];

        dx2 = x[p1] - xcc[le];
        dy2 = y[p1] - ycc[le];
        dz2 = z[p1] - zcc[le];

        nx = dy1 * dz2 - dy2 * dz1;
        ny = dz1 * dx2 - dz2 * dx1;
        nz = dx1 * dy2 - dx2 * dy1;

        fc = fl + f[p1] + f[p2];

        //! Face centers.
        xc = xcc[le] + x[p1] + x[p2];
        yc = ycc[le] + y[p1] + y[p2];
        zc = zcc[le] + z[p1] + z[p2];

        tmp   = nx * xc + ny * yc + nz * zc;
        vvol -= tmp;

        dfdx += nx * fc;
        dfdy += ny * fc;
        dfdz += nz * fc;

        //! Right.
        dx1 = x[p2] - x[p1];
        dy1 = y[p2] - y[p1];
        dz1 = z[p2] - z[p1];

        dx2 = xcc[re] - x[p1];
        dy2 = ycc[re] - y[p1];
        dz2 = zcc[re] - z[p1];

        nx = dy1 * dz2 - dy2 * dz1;
        ny = dz1 * dx2 - dz2 * dx1;
        nz = dx1 * dy2 - dx2 * dy1;

        fc = fr + f[p1] + f[p2];

        dfdx += nx * fc;
        dfdy += ny * fc;
        dfdz += nz * fc;

        //! Face centers.
        xc = xcc[re] + x[p1] + x[p2];
        yc = ycc[re] + y[p1] + y[p2];
        zc = zcc[re] + z[p1] + z[p2];

        tmp   = nx * xc + ny * yc + nz * zc;
        vvol -= tmp;
    }

    //! Here dfdx should be divided by 6,and vol should be divided by(3*2*3)
    //! (dfdx/6)/(vol/18)=dfdx/(vol/3).

    //! Notice: this is not to say that the real vol is vvol/3,just for convenience.

    RDouble coef = 3.0 / vvol;

    dfdx *= coef;
    dfdy *= coef;
    dfdz *= coef;

}

LIB_EXPORT void UnstructGrid::WriteFaceBC(fstream & file)
{
    VirtualFile *virtualFile = new VirtualFile(&file);
    virtualFile->BeginWriteWork();

    int iZone = this->GetGridID()->GetIndex();

    int nBoundFace = this->GetNBoundFace();

    PHWrite(virtualFile, iZone);
    int *bcTypeList    = new int [nBoundFace];
    string *bcNameList = new string [nBoundFace];

    UnstructBCSet **bcr = this->GetBCRecord();

    uint_long totalSize = 0;
    for (cgsize_t iFace = 0; iFace < nBoundFace; iFace ++)
    {
        int bcType = bcr[iFace]->GetKey();
        const string &bcName = bcr[iFace]->GetBCName();

        bcTypeList[iFace] = bcType;
        bcNameList[iFace] = bcName;

        totalSize += static_cast< uint_long > (bcName.size());
        totalSize += 1;
    }

    char *bcNameChar = new char [totalSize];
    unsigned int count = 0;
    for (cgsize_t iFace = 0; iFace < nBoundFace; iFace ++)
    {
        string & bcName = bcNameList[iFace];
        streamsize nameSize = bcName.size();
        for (unsigned int iChar = 0; iChar < nameSize; ++ iChar)
        {
            bcNameChar[count ++] = bcName[iChar];
        }
        bcNameChar[count ++] = '\0';
    }
    PHWrite(virtualFile, nBoundFace);
    PHWrite(virtualFile, totalSize);
    PHWrite(virtualFile, bcNameChar, static_cast<int>(totalSize));

    PHWrite(virtualFile, nBoundFace);
    PHWrite(virtualFile, bcTypeList, nBoundFace);

    delete [] bcTypeList;    bcTypeList = nullptr;
    delete [] bcNameList;    bcNameList = nullptr;
    delete [] bcNameChar;    bcNameChar = nullptr;

    virtualFile->EndWriteWork();

    delete virtualFile;    virtualFile = nullptr;
}

LIB_EXPORT void UnstructGrid::CoarseGrids(int maxLevel)
{
    //++++++++++++++++++++++++++++++++++++++++++++++
    //+ level = 0, the most fine grid.
    //+ level = 1, 2, 3, ..., the coarser grid.
    //++++++++++++++++++++++++++++++++++++++++++++++
    int level = this->GetLevel();
    WriteLogFile("Start coarsen grids of level: ", level + 1);

    if (level > maxLevel)
    {
        this->SetCoarseGrid(0);
        return;
    }
    
    const int MY_METHOD       = 1;
    const int MGRIDGEN_METHOD = 2;
    int coarseMethod = MGRIDGEN_METHOD;

    int nTotalCell = this->GetNTotalCell();

    UnstructGrid *coarseGrid = 0;
    UnstructGrid *preGrid = this;

    int nDim = PHSPACE::GetDim();
    ASSERT(nDim > 1);

    RDouble MinCellNumberRatio = pow(2.0, nDim) * 0.85;

    RDouble cellNumberRatio = 0.0;
    for (int iDim = 0; cellNumberRatio < MinCellNumberRatio; ++ iDim)
    {
        WriteLogFile("Coarse time:", iDim);

        Mesh_Agglomeration *meshAgglomeration = new Mesh_Agglomeration(preGrid, level);
        if(coarseMethod == MGRIDGEN_METHOD)
        {
            meshAgglomeration->CoarsenGridOncebyMGridgen();
        }
        else
        {
            meshAgglomeration->CoarsenGridOnce();
        }
        
        
        UnstructGrid *agglomeratedGrid = meshAgglomeration->GetCoarseGrid();
        delete meshAgglomeration;

        CoarseGridConnectUNMergeFace(preGrid, agglomeratedGrid);
        agglomeratedGrid->ComputeMetrics();

        agglomeratedGrid->SetFineGrid(this);
        agglomeratedGrid->SetLevel(this->GetLevel()+1);
        agglomeratedGrid->SetBCRecord(this->GetBCRecord());
        agglomeratedGrid->SetUnstructBCSet(this->GetUnstructBCSet());
        agglomeratedGrid->SetVolumeCondition(this->GetVolumeConditionIn());

        if(preGrid != this)
        {
            int *pre_cell2coarsegridcell = preGrid->GetCell2CoarseGridCell();
            for (int iCell = 0; iCell < nTotalCell; ++ iCell)
            {
                int preGridCellID = cell2coarsegridcell[iCell];
                cell2coarsegridcell[iCell] = pre_cell2coarsegridcell[preGridCellID];
            }
        }

        bool isCoarsenLimit = (preGrid->GetNTotalCell() == agglomeratedGrid->GetNTotalCell()) || (iDim > 10);

        UnstructGrid *tempGrid = preGrid;
        preGrid = agglomeratedGrid;
        coarseGrid = agglomeratedGrid;

        cellNumberRatio = (nTotalCell * 1.0) / (coarseGrid->GetNTotalCell() * 1.0);

        if(iDim > 0)
        {
            delete tempGrid;    tempGrid = nullptr;
        }

        if(isCoarsenLimit || coarseMethod == MGRIDGEN_METHOD)
        {
             //! MGridgen lib will controll the ratio automatically.
            break;
        }
    }

    coarseGrid->SetLevel(this->GetLevel() + 1);
    coarseGrid->SetFineGrid(this);
    this->SetCoarseGrid(coarseGrid);

    coarseGrid->SkewnessSummary();

    if(PHMPI::GetNumberofGlobalZones() == PHMPI::GetNumberOfProcessor())
    {
        RDouble globalMinCellNumberRatio = cellNumberRatio, globalMaxCellNumberRatio = cellNumberRatio;
        PH_CompareMaxMin(globalMaxCellNumberRatio, 1);
        PH_CompareMaxMin(globalMinCellNumberRatio, 2);

        ostringstream oss;
        oss << "number of cell in level " << this->GetLevel() << ": " << nTotalCell << endl;
        oss << "number of cell in level " << this->GetLevel() + 1 << ": " << this->GetCoarseGrid()->GetNTotalCell() << endl;
        oss << "Global Min && Max cell number ratio of coarse to fine grid: " 
            << globalMinCellNumberRatio << ", " << globalMaxCellNumberRatio << endl;
        PrintToWindow(oss);
        WriteLogFile(oss);

    }
}

//! Generate 3D grid connectivity based on cell2coarsegridcell without merging faces.
void CoarseGridConnectUNMergeFace(UnstructGrid *fgrid, UnstructGrid *cgrid)
{
    int nTotalFace = fgrid->GetNTotalFace();
    int nBoundFace = fgrid->GetNBoundFace();
    int *leftCellOfFace = fgrid->GetLeftCellOfFace();
    int *rightCellOfFace = fgrid->GetRightCellOfFace();
    int *cell2coarsegridcell = fgrid->GetCell2CoarseGridCell();
    int cnTCell = cgrid->GetNTotalCell();

    cgrid->SetNBoundFace(fgrid->GetNBoundFace());

    int *fmark = new int[ nTotalFace ];
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        fmark[iFace] = 1;
    }

    int le, re;
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [ iFace ];
        re = rightCellOfFace[ iFace ];
        if (cell2coarsegridcell[le] == cell2coarsegridcell[re])
        {
            fmark[iFace] = 0;
        }
    }

    int cnTFace = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        if (fmark[iFace] == 1) 
        {
            ++ cnTFace;
        }
    }

    cgrid->SetNTotalFace(cnTFace);
    int *node_number_of_each_face  = fgrid->GetNodeNumberOfEachFace();
    int *face2node = fgrid->GetFace2Node();

    int *cnode_number_of_each_face = new int[ cnTFace ];
    cgrid->SetNodeNumberOfEachFace(cnode_number_of_each_face);
    int *cleft_cell_of_face = new int[ cnTFace ];
    int *cright_cell_of_face = new int[ cnTFace ];
    cgrid->SetLeftCellOfFace(cleft_cell_of_face);
    cgrid->SetRightCellOfFace(cright_cell_of_face);

    int count = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace[ iFace ];
        cnode_number_of_each_face[ count ] = node_number_of_each_face[ iFace ];

        cleft_cell_of_face [ count ] = cell2coarsegridcell[le];
        cright_cell_of_face[ count ] = cnTCell + iFace;
        ++ count;
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        if (fmark[iFace] == 1)
        {
            le = leftCellOfFace [ iFace ];
            re = rightCellOfFace[ iFace ];
            cnode_number_of_each_face[ count ] = node_number_of_each_face[iFace];
            cleft_cell_of_face [ count ] = cell2coarsegridcell[ le ];
            cright_cell_of_face[ count ] = cell2coarsegridcell[ re ];
            ++ count;
        }
    }

    count = 0;
    for (int i = 0; i < cnTFace; ++ i) count += cnode_number_of_each_face[i];

    int *cface2node = new int[count];
    cgrid->SetFace2Node(cface2node);

    count = 0;
    int countf = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            cface2node[ count ++ ] = face2node[ countf ++ ];
        }
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        if (fmark[iFace] == 1)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                cface2node[ count ++ ] = face2node[ countf ++ ];
            }
        }
        else
        {
            countf += node_number_of_each_face[iFace];
        }
    }

    delete [] fmark;    fmark = nullptr;
}

void FieldVisualizationForVTK(Grid *grid_in, std::ostringstream &oss, vector<string> &title, RDouble **qn, int nl)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);

    int nTotalNode = grid->GetNTotalNode();
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();

    int nvarplot = nl;

    int *face2node = grid->GetFace2Node();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    //! Find out cellToFace.
    HyList < int > *cell2face = new HyList < int > [ nTotalCell ];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        cell2face[ iCell ].SetAverNnode(4);
    }
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[ iFace ];
        cell2face[ le ].insert(iFace);
    }
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace[ iFace ];
        int re = rightCellOfFace[ iFace ];
        cell2face[ le ].insert(iFace);
        cell2face[ re ].insert(iFace);
    }

    //! Find out cellToNode.
    HyList < int > *cell2node = new HyList < int > [ nTotalCell ];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        cell2node[ iCell ].SetAverNnode(4);
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int face0 = cell2face[ iCell ][ 0 ];
        int ip1 = face2node[ 2 * face0 ];
        int ip2 = face2node[ 2 * face0 + 1 ];
        cell2node[ iCell ].insert(ip1);
        cell2node[ iCell ].insert(ip2);

        int targetNode = ip2;
        int checkNode  = ip1;
        
        int nFaceOfCell = cell2face[ iCell ].size();
        bool isClockWise = false;
        while (! isClockWise)
        {
            for (int iFace = 1; iFace < nFaceOfCell; ++ iFace)
            {
                int face  = cell2face[ iCell ][ iFace ];
                int node1 = face2node[ 2 * face ];
                int node2 = face2node[ 2 * face + 1 ];

                if (node1 == targetNode)
                {
                    if(node2 != checkNode)
                    {
                        cell2node[ iCell ].insert(node2);
                        targetNode = node2;
                    }
                    else
                    {
                        isClockWise = true;
                        break;
                    }
                }
                else if (node2 == targetNode)
                {
                    if(node1 != checkNode)
                    {
                        cell2node[ iCell ].insert(node1);
                        targetNode = node1;
                    }
                    else
                    {
                        isClockWise = true;
                        break;
                    }
                    
                }
            }    //! Iterate over all faces of cell.
        }        //! Judge whether the nodes of cell have searched closed.
    }            //! iCell.

    //! Output to file.
    using namespace PHENGLEI;
    int zoneid = grid->GetZoneID();

    if (zoneid == 0)
    {
        for (std::size_t i = 0; i < title.size(); ++ i)
        {
            oss << title[i] << "\n";
        }
    }

    std::ostringstream os_zone;
    os_zone << "\"zone = " << grid->GetZoneID() << "\" ";
    string zone_title = os_zone.str();

    oss << "zone T = " << zone_title << " N = " << nTotalNode << "  E = " << nTotalCell  << " f = FEPOINT, ET = quadrilateral\n";
    int wordwidth = 20;
    for (int i = 0; i < nTotalNode; ++ i)
    {
        int it = i;
        oss << setiosflags(ios::left);
        oss << setiosflags(ios::scientific);
        oss << setprecision(10);
        oss << setw(wordwidth) << x[it]
        << setw(wordwidth) << y[it]
        << setw(wordwidth) << z[it];
        for (int m = 0; m < nvarplot; ++ m)
        {
            oss << setw(wordwidth) << qn[m][it];
        }
        oss << "\n";
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (cell2node[ iCell ].size() == 3)
        {
            for (int m = 0; m < 3; ++ m)
            {
                oss << cell2node[iCell].GetData(m) + 1<< "  ";
            }
            oss << cell2node[iCell].GetData(2) + 1;
        }
        else if (cell2node[ iCell ].size() == 4)
        {
            for (int m = 0; m < 4; ++ m)
            {
                oss << cell2node[iCell].GetData(m) + 1<< "  ";
            }
        }
        else
        {
            cout << "Error: this function only support triangle and quad !\n";
            exit(0);
        }
        oss << endl;
    }

    delete [] cell2node;    cell2node = nullptr;
    delete [] cell2face;    cell2face = nullptr;
}
void SaveDataForTecio_bac(Grid *grid_in, ActionKey *actkey, RDouble **qn, int nvarplot)
{
    DataContainer *cdata = actkey->GetTecData();

    UnstructGrid *grid = UnstructGridCast(grid_in);

    int TecioMission = WriteBlock;
    PHWrite(cdata, &TecioMission, 1);
    
    int IMxOrNumPts, JMxOrNumElements, KMxOrNumFaces;
    IMxOrNumPts      = grid->GetNTotalNode();
    JMxOrNumElements = grid->GetNTotalCell();
    KMxOrNumFaces    = grid->GetNTotalFace();
    PHWrite(cdata, &IMxOrNumPts, 1);
    PHWrite(cdata, &JMxOrNumElements, 1);
    PHWrite(cdata, &KMxOrNumFaces, 1);

    int *node_number_of_each_face  = grid->GetNodeNumberOfEachFace();
    int TotalNumFaceNodes_Rect = 0;
    for (int iFace = 0; iFace < KMxOrNumFaces; ++ iFace)
    {
        TotalNumFaceNodes_Rect += node_number_of_each_face[iFace];
    }
    PHWrite(cdata, &TotalNumFaceNodes_Rect, 1);
    PHWrite(cdata, node_number_of_each_face, KMxOrNumFaces);

    int nTotalNode = grid->GetNTotalNode();
    PHWrite(cdata, &nTotalNode, 1);

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();
    PHWrite(cdata, x, nTotalNode);
    PHWrite(cdata, y, nTotalNode);
    PHWrite(cdata, z, nTotalNode);

    for (int m = 0; m < nvarplot; ++ m)
    {
        RDouble *qq = qn[m];
        PHWrite(cdata, qq, nTotalNode);
    }

    int *face2node = grid->GetFace2Node();
    
    PHWrite(cdata, face2node, TotalNumFaceNodes_Rect);

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    PHWrite(cdata, leftCellOfFace, KMxOrNumFaces);
    PHWrite(cdata, rightCellOfFace, KMxOrNumFaces);
}

void SaveDataForTecio(Grid *grid_in, ActionKey *actkey, RDouble **qn, int nvarplot)
{
    DataContainer *cdata = actkey->GetTecData();

    UnstructGrid *grid = UnstructGridCast(grid_in);

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();
    int numberOfFaces = grid->GetNTotalFace();

    int numberOfNodesInZone = 0;
    int numberOfFacesInZone = 0;
    int numberOfCellsInZone = 0;

    vector< int > keyNodes;
    vector< int > keyFaces;
    vector< int > nodeList;
    vector< int > cellList;

    keyNodes.resize (numberOfNodes);
    keyFaces.resize (numberOfFaces);
    nodeList.resize (numberOfNodes);
    cellList.resize (numberOfCells);

    SetField(keyNodes,  0, numberOfNodes);
    SetField(keyFaces,  0, numberOfFaces);
    SetField(nodeList, -1, numberOfNodes);
    SetField(cellList, -1, numberOfCells);

    int * keyActiveOfCells        = grid->GetBlankIndex();

    int **cellNodeIndexContainer  = grid->GetCell2NodeArray();
    int * cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int **cellFaceIndexContainer  = grid->GetCell2Face();
    int * cellFaceNumberContainer = grid->GetFaceNumberOfEachCell();

    int *faceNodeNumberContainer  = grid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer   = grid->GetFace2Node();

    int *leftCellIndexContainer   = grid->GetLeftCellOfFace();
    int *rightCellIndexContainer  = grid->GetRightCellOfFace();

    for(int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        //if(keyActiveOfCells[ iCell ] == -1)
        if (keyActiveOfCells[iCell] == ACTIVE)
        {
            cellList[ iCell ] = numberOfCellsInZone;

            for(int iNodes = 0; iNodes < cellNodeNumberContainer[ iCell ]; ++ iNodes)
            {
                keyNodes[ cellNodeIndexContainer[ iCell ][ iNodes ] ] = 1;
            }

            for(int iFaces = 0; iFaces < cellFaceNumberContainer[ iCell ]; ++ iFaces)
            {
                keyFaces[ cellFaceIndexContainer[ iCell ][ iFaces ] ] = 1;
            }
            numberOfCellsInZone ++;
        }
    }

    for(int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if(keyNodes[ iNode ] == 1)
        {
            nodeList[ iNode ] = numberOfNodesInZone;
            numberOfNodesInZone ++;
        }
    }

    for(int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        if(keyFaces[ iFace ] == 1)
        {
            numberOfFacesInZone ++;
        }
    }

    int totalNumFaceNodes = 0;
    for(int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        if(keyFaces[ iFace ] != 1) continue;
        totalNumFaceNodes += faceNodeNumberContainer[ iFace ];
    }

    if(numberOfNodesInZone == 0)
    {
        return;
    }

    int TecioMission = WriteBlock;

    PHWrite(cdata, TecioMission);
    PHWrite(cdata, numberOfNodesInZone);
    PHWrite(cdata, numberOfCellsInZone);
    PHWrite(cdata, numberOfFacesInZone);

    PHWrite(cdata, totalNumFaceNodes);

    int count = 0;
    for(int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        if(keyFaces[ iFace ] != 1) continue;

        int faceNodeNumber = faceNodeNumberContainer[ iFace ];
        PHWrite(cdata, faceNodeNumber);
    }

    PHWrite(cdata, numberOfNodesInZone);

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    for(int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if(keyNodes[ iNode ] != 1) continue;
        PHWrite(cdata, x[iNode]);
    }

    for(int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if(keyNodes[ iNode ] != 1) continue;
        PHWrite(cdata, y[iNode]);
    }

    for(int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if(keyNodes[ iNode ] != 1) continue;
        PHWrite(cdata, z[iNode]);
    }

    for (int iVar = 0; iVar < nvarplot; ++ iVar)
    {
        RDouble *qq = qn[iVar];
        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            if (keyNodes[iNode] != 1) continue;
            PHWrite(cdata, qq[iNode]);
        }
    }

    count = 0;

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int faceNodeNumber = faceNodeNumberContainer[iFace];

        if (keyFaces[iFace] != 1)
        {
            count += faceNodeNumber;
            continue;
        }

        for (int iNode = 0; iNode < faceNodeNumberContainer[iFace]; ++ iNode)
        {
            int nodeIndex = nodeList[faceNodeIndexContainer[count ++]];
            PHWrite(cdata, nodeIndex);
        }
    }

    //! Left cell.
    for(int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        if(keyFaces[ iFace ] != 1) continue;

        int cellIndex = cellList[ leftCellIndexContainer[ iFace ] ];
        PHWrite(cdata, cellIndex);
    }

    //! Right cell.
    for(int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int cellIndex;

        if(keyFaces[ iFace ] == 0) continue;

        int rightCell = rightCellIndexContainer[ iFace ];

        if(rightCell >= numberOfCells || rightCell < 0)
        {
            cellIndex = -1;
        }
        else
        {
            cellIndex = cellList[ rightCell ];
        }

        PHWrite(cdata, cellIndex);
    }
}

void FieldVisualization(Grid *grid_in, std::ostringstream &oss, vector<string> &title, RDouble **qn, int nl)
{
    if (GetDim() == PHSPACE::TWO_D)
    {
        FieldVisualizationForVTK(grid_in, oss, title, qn, nl);
        return;
    }

    UnstructGrid *grid = UnstructGridCast(grid_in);

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    int nTotalNode = grid->GetNTotalNode();
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();

    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int *face2node = grid->GetFace2Node();

    for (std::size_t i = 0; i < title.size(); ++ i)
    {
        oss << title[i] << "\n";
    }

    int count = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        count += node_number_of_each_face[iFace];
    }

    int nword_of_line = 5;

    //! Output for Tecplot.
    oss << "ZONE\n";
    if (grid->GetDim() == THREE_D)
    {
        oss << "ZoneType = FEPolyhedron\n";
    }
    else
    {
        oss << "ZoneType = FEPolygon\n";
    }
    oss << "Nodes    = " << nTotalNode << "\n";
    oss << "Faces    = " << nTotalFace << "\n";
    oss << "Elements = " << nTotalCell << "\n";
    oss << "TotalNumFaceNodes = " << count << "\n";
    oss << "NumConnectedBoundaryFaces = 0\n";
    oss << "TotalNumBoundaryConnections = 0\n";

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        oss << x[iNode] << " ";
        if ((iNode + 1) % nword_of_line == 0) oss << "\n";
    }
    if (nTotalNode % nword_of_line != 0) oss << "\n";

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        oss << y[iNode] << " ";
        if ((iNode + 1) % nword_of_line == 0) oss << "\n";
    }
    if (nTotalNode % nword_of_line != 0) oss << "\n";

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        oss << z[iNode] << " ";
        if ((iNode + 1) % nword_of_line == 0) oss << "\n";
    }
    if (nTotalNode % nword_of_line != 0) oss << "\n";

    for (int m = 0; m < nl; ++ m)
    {
        RDouble *qq = qn[m];
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            oss << qq[iNode] << " ";
            if ((iNode + 1) % nword_of_line == 0) oss << "\n";
        }
        if (nTotalNode % nword_of_line != 0) oss << "\n";
    }

    if (grid->GetDim() == THREE_D)
    {
        //! Geometry.
        //! First the node_number_of_each_face(n poin pf face).
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            oss << node_number_of_each_face[iFace] << " ";
            if ((iFace + 1) % nword_of_line == 0) oss << "\n";
        }
        if (nTotalFace % nword_of_line == 0) oss << "\n";
    }

    //! Second the points to face.
    count = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            oss << face2node[count ++] + 1 << " ";
        }
        oss << "\n";
    }

    //! Third the f2c:left.
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        oss << leftCellOfFace[iFace] + 1 << " ";
        if ((iFace + 1) % nword_of_line == 0) oss << "\n";
    }
    if (nTotalFace % nword_of_line == 0) oss << "\n";

    //! Third the f2c:right.
    int right;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        right = rightCellOfFace[iFace] + 1;
        if (right > nTotalCell || right < 0) right = 0;
        oss << right << " ";
        if ((iFace + 1) % nword_of_line == 0) oss << "\n";
    }
    if (nTotalFace % nword_of_line == 0) oss << "\n";
}

void FixBCNodeVar(UnstructGrid *grid, RDouble *q, RDouble *q_n, int *n_count, int bctype_in, bool twoside)
{
    int nTotalNode = grid->GetNTotalNode();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *face2node = grid->GetFace2Node();
    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    //! Bell 20140409 add.
    //! Find out the corner points(corner_point = 2) of boundary conditions with both interface and solid wall attributes.For this kind of corner points,do not fix.
    //! In order to ensure the values continuity at the interface between solid wall and interface.
    int *corner_point = new int [ nTotalNode ];
    if (bctype_in == PHENGLEI::INTERFACE)
    {
        for (int ipoint = 0; ipoint < nTotalNode; ++ ipoint)
        {
            corner_point[ ipoint ] = 0;
        }

        int nodepos = 0;
        int bcType, pt;
        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
            bcType = bcRegion->GetBCType();
            if (bcType == PHENGLEI::INTERFACE)
            {
                for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
                {
                    pt = face2node[ nodepos + j ];
                    corner_point[ pt ] = 1;
                }
            }
            nodepos += node_number_of_each_face[iFace];
        }

        nodepos = 0;
        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
            bcType = bcRegion->GetBCType();
            if (bcType == PHENGLEI::SOLID_SURFACE)
            {
                for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
                {
                    pt = face2node[ nodepos + j ];
                    if (corner_point[ pt ] == 1)
                    {
                        corner_point[ pt ] = 2;
                    }
                }
            }
            nodepos += node_number_of_each_face[iFace];
        }
    }

    //! Reprocessing the points of interface.
    int le, re, pt, bcType;
    int nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace[ iFace ];
        re = iFace + nTotalCell;

        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        bcType = bcRegion->GetBCType();

        if (bcType == bctype_in)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt = face2node[ nodepos + j ];

                if ((bctype_in == PHENGLEI::INTERFACE) && (corner_point[ pt ] == 2)) continue;    //! Bell 20140409 add.

                q_n[ pt ] = 0.0;
                n_count[ pt ] = 0;
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

    nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace[ iFace ];
        re = iFace + nTotalCell;

        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        bcType = bcRegion->GetBCType();

        if (bcType == bctype_in)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt = face2node[ nodepos + j ];

                if ((bctype_in == PHENGLEI::INTERFACE) && (corner_point[ pt ] == 2)) continue;    //! Bell 20140409 add.

                q_n[ pt ] += q[ le ];
                ++ n_count[ pt ];
                if (twoside)
                {
                    q_n[ pt ] += q[ re ];
                    ++ n_count[ pt ];
                }
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

    //! A corrected method of points on solid wall.
    //! In the case of poor grid quality at the trailing edge of the naca0012 airfoil ,the computation can converge with this method, 
    //! otherwise, the computation divergent.
    //! However, due to the lack of strict assessment,just note here.

    delete [] corner_point;    corner_point = nullptr;
}

//! Bell 20121220 add.
void FixBCNodeVar(UnstructGrid *grid, RDouble *q, RDouble *q_n, RDouble *n_count, int bctype_in, bool twoside)
{
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *face2node = grid->GetFace2Node();
    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    //! Reprocessing the points of interface.
    int le, re, pt, bcType;
    int nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace[ iFace ];
        re = iFace + nTotalCell;

        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        bcType = bcRegion->GetBCType();

        if (bcType == bctype_in)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt = face2node[ nodepos + j ];
                q_n[ pt ] = 0.0;
                n_count[ pt ] = 0.0;
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

    nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace[ iFace ];
        re = iFace + nTotalCell;

        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        bcType = bcRegion->GetBCType();

        if (bcType == bctype_in)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt = face2node[ nodepos + j ];
                q_n[ pt ] += q[ le ];
                ++ n_count[ pt ];
                if (twoside)
                {
                    q_n[ pt ] += q[ re ];
                    ++ n_count[ pt ];
                }
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }
}

void CompNodeVar(UnstructGrid *grid, RDouble *qNode, RDouble *q, bool isVelocityForPostVisual)
{
    int nTotalNode = grid->GetNTotalNode();
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();

    RDouble *xCoor = grid->GetX();
    RDouble *yCoor = grid->GetY();
    RDouble *zCoor = grid->GetZ();

    RDouble *xCellCenter = grid->GetCellCenterX();
    RDouble *yCellCenter = grid->GetCellCenterY();
    RDouble *zCellCenter = grid->GetCellCenterZ();

    RDouble *nodeWeight = new RDouble[nTotalNode]();
    PHSPACE::SetField(nodeWeight, 0.0, nTotalNode);
    PHSPACE::SetField(qNode, 0.0, nTotalNode);

    //! Cell center data insert in node
    int **cell2NodeArray       = grid->GetCell2NodeArray();
    int * nodeNumberOfEachCell = grid->GetNodeNumberOfEachCell();
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int nNodes = nodeNumberOfEachCell[iCell];
        for (int jNode = 0; jNode < nNodes; ++ jNode)
        {
            int nodeIndex = cell2NodeArray[iCell][jNode];

            RDouble dx = xCoor[nodeIndex] - xCellCenter[iCell];
            RDouble dy = yCoor[nodeIndex] - yCellCenter[iCell];
            RDouble dz = zCoor[nodeIndex] - zCellCenter[iCell];
            RDouble dist = DISTANCE(dx, dy, dz);
            RDouble weightTemp = 1.0 / dist;

            RDouble cellQtemp = q[iCell];
            qNode[nodeIndex] += cellQtemp * weightTemp;
            nodeWeight[nodeIndex] += weightTemp;
        }
    }

    //! Boundary
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int *nodeNumberOfEachFace = grid->GetNodeNumberOfEachFace();
    int **face2nodeArray      = grid->GetFace2NodeArray();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        if (bcType == PHENGLEI::SYMMETRY)
        {
            continue;
        }

        if ((bcType == PHENGLEI::SOLID_SURFACE) && isVelocityForPostVisual)
        {
            continue;
        }

        int cellIndex = rightCellOfFace[iFace];
        for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
        {
            int nodeIndex = face2nodeArray[iFace][jNode];

            RDouble dx = xCoor[nodeIndex] - xCellCenter[cellIndex];
            RDouble dy = yCoor[nodeIndex] - yCellCenter[cellIndex];
            RDouble dz = zCoor[nodeIndex] - zCellCenter[cellIndex];
            RDouble dist = DISTANCE(dx, dy, dz);
            RDouble weightTemp = 1.0 / dist;

            RDouble cellQtemp = q[cellIndex];
            qNode[nodeIndex] += cellQtemp * weightTemp;
            nodeWeight[nodeIndex] += weightTemp;
        }
    }

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (nodeWeight[iNode] > SMALL)
        {
            qNode[iNode] /= nodeWeight[iNode];
        }
        else
        {
            qNode[iNode] = 0.0;
        }
    }

    delete [] nodeWeight;    nodeWeight = nullptr;
}

void CompNodeVar(UnstructGrid *grid, RDouble *qNode, const string &variableName, RDouble *q, bool isVelocityForPostVisual)
{
    int nTotalNode = grid->GetNTotalNode();
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();

    RDouble *xCoor = grid->GetX();
    RDouble *yCoor = grid->GetY();
    RDouble *zCoor = grid->GetZ();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    RDouble *xCellCenter = grid->GetCellCenterX();
    RDouble *yCellCenter = grid->GetCellCenterY();
    RDouble *zCellCenter = grid->GetCellCenterZ();

    RDouble *nodeWeight = new RDouble[nTotalNode]();
    PHSPACE::SetField(nodeWeight, 0.0, nTotalNode);
    PHSPACE::SetField(qNode, 0.0, nTotalNode);

    int  *leftCellOfFace       = grid->GetLeftCellOfFace();
    int  *rightCellOfFace      = grid->GetRightCellOfFace();
    int  *nodeNumberOfEachFace = grid->GetNodeNumberOfEachFace();
    int **face2nodeArray       = grid->GetFace2NodeArray();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    //! Node BC type.
    int *nodeBC = NewPointer< int >(nTotalNode);
    PHSPACE::SetField(nodeBC, PHENGLEI::NO_BOUNDARY_CONDITION, nTotalNode);

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter ++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2nodeArray[iFace][jNode];
                    nodeBC[point] = bcType;
                }
            }
        }
    }

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType == PHENGLEI::FARFIELD)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter ++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;

                //! Far field has the second highest priority.
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2nodeArray[iFace][jNode];
                    if (nodeBC[point] == PHENGLEI::SOLID_SURFACE)
                    {
                        continue;
                    }

                    nodeBC[point] = bcType;
                }
            }
        }
    }

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter ++)
        {
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            if (bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::FARFIELD ||
                bcType == PHENGLEI::SYMMETRY || bcType == PHENGLEI::INTERFACE)
            {
                continue;
            }

            //! other BC, except symmetry and interface.
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                int point = face2nodeArray[iFace][jNode];
                if (nodeBC[point] == PHENGLEI::SOLID_SURFACE || nodeBC[point] == PHENGLEI::FARFIELD)
                {
                    continue;
                }

                nodeBC[point] = bcType;
            }
        }
    }

    //! Boundary faces.
    RDouble *qOnBCFace = new RDouble [nBoundFace];
    PHSPACE::SetField(qOnBCFace, zero, nBoundFace);

    //! Init q on BC faces.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter ++)
        {
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            if (bcType == PHENGLEI::INTERFACE || bcType == PHENGLEI::SYMMETRY)
            {
                continue;
            }

            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];

            qOnBCFace[iFace] = half * (q[le] + q[re]);

        }
    }

    //! Step1: solid surface.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter ++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2nodeArray[iFace][jNode];

                    RDouble dx = xCoor[point] - xfc[iFace];
                    RDouble dy = yCoor[point] - yfc[iFace];
                    RDouble dz = zCoor[point] - zfc[iFace];
                    RDouble dist = DISTANCE(dx, dy, dz);
                    RDouble weightTemp = 1.0 / dist;

                    RDouble faceQtemp = qOnBCFace[iFace];
                    qNode[point] += faceQtemp * weightTemp;

                    nodeWeight[point] += weightTemp;
                }
            }
        }
    }

    //! Step2: Far-field.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType == PHENGLEI::FARFIELD)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter ++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2nodeArray[iFace][jNode];
                    if (nodeBC[point] == PHENGLEI::SOLID_SURFACE)
                    {
                        continue;
                    }

                    RDouble dx = xCoor[point] - xfc[iFace];
                    RDouble dy = yCoor[point] - yfc[iFace];
                    RDouble dz = zCoor[point] - zfc[iFace];
                    RDouble dist = DISTANCE(dx, dy, dz);
                    RDouble weightTemp = 1.0 / dist;

                    RDouble faceQtemp = qOnBCFace[iFace];
                    qNode[point] += faceQtemp * weightTemp;

                    nodeWeight[point] += weightTemp;
                }
            }
        }
    }

    //! Step3: other BC except solid surface/far field/symmetry/Interface.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType != PHENGLEI::SOLID_SURFACE && bcType != PHENGLEI::FARFIELD &&
            bcType != PHENGLEI::SYMMETRY && bcType != PHENGLEI::INTERFACE)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter ++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2nodeArray[iFace][jNode];
                    if (nodeBC[point] == PHENGLEI::SOLID_SURFACE || nodeBC[point] == PHENGLEI::FARFIELD)
                    {
                        continue;
                    }

                    RDouble dx = xCoor[point] - xfc[iFace];
                    RDouble dy = yCoor[point] - yfc[iFace];
                    RDouble dz = zCoor[point] - zfc[iFace];
                    RDouble dist = DISTANCE(dx, dy, dz);
                    RDouble weightTemp = 1.0 / dist;

                    RDouble faceQtemp = qOnBCFace[iFace];
                    qNode[point] += faceQtemp * weightTemp;

                    nodeWeight[point] += weightTemp;
                }
            }
        }
    }

    //! Step4: Now the interior points.
    //!        Importantly, the symmetry points are used as interior points.
    int **cell2NodeArray = grid->GetCell2NodeArray();
    int  *nodeNumberOfEachCell = grid->GetNodeNumberOfEachCell();
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int nNode = nodeNumberOfEachCell[iCell];
        for (int jNode = 0; jNode < nNode; ++ jNode)
        {
            int point = cell2NodeArray[iCell][jNode];
            if (nodeBC[point] != PHENGLEI::NO_BOUNDARY_CONDITION)
            {
                continue;
            }

            //! NO_BC && Symmetry && Interface?
            RDouble dx = xCoor[point] - xCellCenter[iCell];
            RDouble dy = yCoor[point] - yCellCenter[iCell];
            RDouble dz = zCoor[point] - zCellCenter[iCell];
            RDouble dist = DISTANCE(dx, dy, dz);
            RDouble weightTemp = 1.0 / dist;

            RDouble cellQtemp = q[iCell];
            qNode[point] += cellQtemp * weightTemp;

            nodeWeight[point] += weightTemp;
        }
    }

    InterpointInformation *interPointInfor = grid->GetInterpointInfo();
    if (interPointInfor)
    {
        int numberOfInterpoints = interPointInfor->GetNumberOfInterpoints();
        int *interPoint2GlobalPoint = interPointInfor->GetInterPoint2GlobalPoint();
        int *labelOfInterPoint = interPointInfor->GetLabelOfInterPoint();

        for (int iPoint = 0; iPoint < numberOfInterpoints; ++ iPoint)
        {
            int globalPoint = interPoint2GlobalPoint[iPoint];
            if (labelOfInterPoint[iPoint] != 1)
            {
                    qNode[globalPoint] = 0;
            }
        }
    }

    delete [] qOnBCFace;    qOnBCFace = nullptr;
    DelPointer(nodeBC);

    RDouble **variable = new RDouble* [1];
    variable[0] = qNode;

    string nodeDataName = variableName + "_node";
    CommunicateInterpointValue(grid, variable, nodeDataName, 1);

    RDouble **nodeWeightArray = new RDouble* [1];
    nodeWeightArray[0] = nodeWeight;

    string nodeWeightName = "nodeWeight_node";
    CommunicateInterpointValue(grid, nodeWeightArray, nodeWeightName, 1);

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (nodeWeight[iNode] > SMALL)
        {
            qNode[iNode] /= nodeWeight[iNode];
        }
        else
        {
            qNode[iNode] = 0.0;
        }
    }

    delete [] variable;           variable = nullptr;
    delete [] nodeWeightArray;    nodeWeightArray = nullptr;
    delete [] nodeWeight;         nodeWeight = nullptr;
}

//! Compute node values from cell information (cxh, 2012.12.12).
void CompNodeVar_new(UnstructGrid *grid, RDouble *q_n, RDouble *q)
{
    //! Get necessary geometry information.
    int nTotalNode = grid->GetNTotalNode();
    int *nCPN = grid->GetCellNumberOfEachNode();
    int *n2c  = grid->GetNode2Cell();
    int iNode, j, iCell, count;
    int *knode;

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();
    RDouble *xn  = grid->GetX();
    RDouble *yn  = grid->GetY();
    RDouble *zn  = grid->GetZ();

    RDouble weight, weight1, weight2, *n_count;
    RDouble *lamdax = grid->GetLamdax();
    RDouble *lamday = grid->GetLamday();
    RDouble *lamdaz = grid->GetLamdaz();

    if (lamdax == nullptr || lamday == nullptr || lamdaz == nullptr)
    {
        grid->CalcLaplacianWeitht();
        lamdax = grid->GetLamdax();
        lamday = grid->GetLamday();
        lamdaz = grid->GetLamdaz();
    }

    knode = grid->Getknode();

    n_count = new RDouble [ nTotalNode ];
    for (iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        q_n[ iNode ]     = 0.0;
        n_count[ iNode ] = 0.0;
    }
    count = 0;
    for (iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        for (j = 0; j < nCPN[ iNode ]; ++j)
        {
            iCell = n2c[ count++ ];
            weight1 = 1 + lamdax[ iNode ] * (xcc[ iCell ] - xn[ iNode ])
                        + lamday[ iNode ] * (ycc[ iCell ] - yn[ iNode ])
                        + lamdaz[ iNode ] * (zcc[ iCell ] - zn[ iNode ]);
            weight2 = 1 / sqrt((xcc[ iCell ] - xn[ iNode ]) * (xcc[ iCell ] - xn[ iNode ])
                              + (ycc[ iCell ] - yn[ iNode ]) * (ycc[ iCell ] - yn[ iNode ])
                              + (zcc[ iCell ] - zn[ iNode ]) * (zcc[ iCell ] - zn[ iNode ]));
            weight  = weight1 * (1 - knode[ iNode ]) + weight2 * knode[ iNode ];

            q_n[ iNode ] += weight * q[ iCell ];
            n_count[ iNode ] += weight;
        }
    }

    for (iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        q_n[iNode] /= n_count[iNode] + TINY;
    }

    delete [] n_count;    n_count = nullptr;
}

void ComputeNodeVariable(UnstructGrid *grid, RDouble **q_n, RDouble **q, int nEquation)
{
    int nTotalNode = grid->GetNTotalNode();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int *face2node = grid->GetFace2Node();
    int *node_number_of_each_face  = grid->GetNodeNumberOfEachFace();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    //! Get gradients.

    RDouble *n_count = new RDouble [nTotalNode];

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)    
    {
        n_count[iNode] = 0.0;
        for (int m = 0; m < nEquation; ++ m)
        {
            q_n[m][iNode] = 0.0;
        }
    }

    RDouble dx, dy, dz, dd;
    int le, re, pt;

    //! Boundary faces.
    int nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            pt = face2node[nodepos + j];
            dx = x[pt] - xcc[le];
            dy = y[pt] - ycc[le];
            dz = z[pt] - zcc[le];
            dd = 1.0 / (sqrt(dx * dx + dy * dy + dz * dz) + TINY);
            n_count[pt] += dd;
            for (int m = 0; m < nEquation; ++ m)
            {
                q_n[m][pt] += q[m][le] * dd;
            }

            dx = x[pt] - xcc[re];
            dy = y[pt] - ycc[re];
            dz = z[pt] - zcc[re];
            dd = 1.0 / (sqrt(dx * dx + dy * dy + dz * dz) + TINY);
            n_count[pt] += dd;
            for (int m = 0; m < nEquation; ++ m)
            {
                q_n[m][pt] += q[m][re] * dd;
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

    //! Interior faces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            //! From left.
            pt = face2node[nodepos + j];
            dx = x[pt] - xcc[le];
            dy = y[pt] - ycc[le];
            dz = z[pt] - zcc[le];
            dd = 1.0 / (sqrt(dx * dx + dy * dy + dz * dz) + TINY);
            n_count[pt] += dd;
            for (int m = 0; m < nEquation; ++ m)
            {
                q_n[m][pt] += q[m][le] * dd;
            }

            //! From right.
            dx = x[pt] - xcc[re];
            dy = y[pt] - ycc[re];
            dz = z[pt] - zcc[re];
            dd = 1.0 / (sqrt(dx * dx + dy * dy + dz * dz) + TINY);
            n_count[pt] += dd;
            for (int m = 0; m < nEquation; ++ m)
            {
                q_n[m][pt] += q[m][re] * dd;
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

    for (int m = 0; m < nEquation; ++ m)
    {
        FixBCNodeVar(grid, q[m], q_n[m], n_count, PHENGLEI::SYMMETRY     , true);
        FixBCNodeVar(grid, q[m], q_n[m], n_count, PHENGLEI::SOLID_SURFACE, true);
        FixBCNodeVar(grid, q[m], q_n[m], n_count, PHENGLEI::INTERFACE    , true);
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            q_n[m][iNode] /= (n_count[iNode] + SMALL);
        }
    }

    delete [] n_count;    n_count = nullptr;
}

/************************************************************************/
/*                Compute node value with distance weighted             */
/*                           Bell 20121220 add                          */
/************************************************************************/
void CompNodeVarWeight(UnstructGrid *grid, RDouble *q_n, RDouble *q)
{
    int nTotalNode = grid->GetNTotalNode();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int nTotalCell = grid->GetNTotalCell();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int *face2node = grid->GetFace2Node();
    int *node_number_of_each_face  = grid->GetNodeNumberOfEachFace();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *n_count = new RDouble [nTotalNode];

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        n_count[iNode] = 0.0;
        q_n[iNode]     = 0.0;
    }

    RDouble dx, dy, dz, dd;
    int le, re, point;

    //! Boundary faces.
    int nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace[iFace];
        re = iFace + nTotalCell;
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            point = face2node[nodepos + j];

            dx = x[point] - xcc[le];
            dy = y[point] - ycc[le];
            dz = z[point] - zcc[le];
            dd = 1.0 / (sqrt(dx * dx + dy * dy + dz * dz) + TINY);
            q_n[point] += q[le] * dd;
            n_count[point] += dd;

            dx         = x[ point ] - xcc[ re ];
            dy         = y[ point ] - ycc[ re ];
            dz         = z[ point ] - zcc[ re ];
            dd         = 1.0 / (sqrt(dx * dx + dy * dy + dz * dz) + TINY);
            q_n[ point ]   += q[ re ] * dd;
            n_count[point] += dd;
        }
        nodepos += node_number_of_each_face[iFace];
    }

    //! Interior faces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le  = leftCellOfFace [ iFace ];
        re  = rightCellOfFace[ iFace ];
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            point      = face2node[ nodepos + j ];

            dx         = x[ point ] - xcc[ le ];
            dy         = y[ point ] - ycc[ le ];
            dz         = z[ point ] - zcc[ le ];
            dd         = 1.0 / (sqrt(dx * dx + dy * dy + dz * dz) + TINY);
            q_n[ point ]   += q[ le ] * dd;
            n_count[point] += dd;

            dx         = x[ point ] - xcc[ re ];
            dy         = y[ point ] - ycc[ re ];
            dz         = z[ point ] - zcc[ re ];
            dd         = 1.0 / (sqrt(dx * dx + dy * dy + dz * dz) + TINY);
            q_n[ point ]   += q[ re ] * dd;
            n_count[point] += dd;
        }
        nodepos += node_number_of_each_face[iFace];
    }

    FixBCNodeVar(grid, q, q_n, n_count, PHENGLEI::SYMMETRY     , true);
    FixBCNodeVar(grid, q, q_n, n_count, PHENGLEI::SOLID_SURFACE, true);
    FixBCNodeVar(grid, q, q_n, n_count, PHENGLEI::INTERFACE    , true);

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        q_n[iNode] /= (n_count[iNode] + SMALL);
    }

    delete [] n_count;    n_count = nullptr;
}

void CompNodeVarForVisual(UnstructGrid *grid, RDouble *q_n, RDouble *q)
{
    int nTotalNode = grid->GetNTotalNode();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int *face2node = grid->GetFace2Node();
    int *node_number_of_each_face  = grid->GetNodeNumberOfEachFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nBoundFace + nTotalCell;

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    //! Get gradients.
    RDouble *dqdx = new RDouble[nTotal];
    RDouble *dqdy = new RDouble[nTotal];
    RDouble *dqdz = new RDouble[nTotal];

    int order;
    GlobalDataBase::GetData("order", &order, PHINT, 1);

    int *n_count = new int [nTotalNode];

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        n_count[iNode] = 0;
        q_n[iNode]     = 0.0;
    }

    RDouble dx, dy, dz;
    int le, re, pt;

    //! Boundary faces.
    int nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace [ iFace ];
        re = rightCellOfFace[ iFace ];
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            pt         = face2node[ nodepos + j ];
            dx         = x[ pt ] - xcc[ le ];
            dy         = y[ pt ] - ycc[ le ];
            dz         = z[ pt ] - zcc[ le ];
            q_n[ pt ] += q[ le ];
            ++ n_count[ pt ];
            q_n[ pt ] += q[ re ];
            ++ n_count[ pt ];
        }
        nodepos += node_number_of_each_face[iFace];
    }

    //! Interior faces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le  = leftCellOfFace [ iFace ];
        re  = rightCellOfFace[ iFace ];
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            //! From left.
            pt = face2node[ nodepos + j ];
            q_n[ pt ] += q[ le ];
            ++ n_count[ pt ];

            //! From right.
            q_n[ pt ] += q[ re ];
            ++ n_count[ pt ];
        }
        nodepos += node_number_of_each_face[iFace];
    }

    FixBCNodeVar(grid, q, q_n, n_count, PHENGLEI::SYMMETRY     , true);
    FixBCNodeVar(grid, q, q_n, n_count, PHENGLEI::SOLID_SURFACE, false);
    FixBCNodeVar(grid, q, q_n, n_count, PHENGLEI::INTERFACE    , true);

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        q_n[iNode] /= (n_count[iNode] + SMALL);
    }

    delete [] n_count;    n_count = nullptr;
    delete [] dqdx;    dqdx = nullptr;
    delete [] dqdy;    dqdy = nullptr;
    delete [] dqdz;    dqdz = nullptr;
}

void CompNodeVarLimit(UnstructGrid *grid, RDouble *q_n, RDouble *q)
{
    int nTotalNode = grid->GetNTotalNode();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int * leftCellOfFace = grid->GetLeftCellOfFace();
    int * rightCellOfFace = grid->GetRightCellOfFace();
    int *face2node = grid->GetFace2Node();
    int *node_number_of_each_face  = grid->GetNodeNumberOfEachFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nBoundFace + nTotalCell;

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    //! Get gradients.
    RDouble *dqdx = new RDouble[nTotal];
    RDouble *dqdy = new RDouble[nTotal];
    RDouble *dqdz = new RDouble[nTotal];

    int order;
    GlobalDataBase::GetData("order", &order, PHINT, 1);

    int *n_count = new int [nTotalNode];

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        n_count[iNode] = 0;
        q_n[iNode]     = 0.0;
    }

    RDouble dx, dy, dz;
    int le, re, pt;

    //! Boundary faces.
    int nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace [ iFace ];
        re = rightCellOfFace[ iFace ];
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            pt         = face2node[ nodepos + j ];
            dx         = x[ pt ] - xcc[ le ];
            dy         = y[ pt ] - ycc[ le ];
            dz         = z[ pt ] - zcc[ le ];
            q_n[ pt ] += q[ le ];
            ++ n_count[ pt ];
            q_n[ pt ] += q[ re ];
            ++ n_count[ pt ];
        }
        nodepos += node_number_of_each_face[iFace];
    }

    //! Interior faces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le  = leftCellOfFace [ iFace ];
        re  = rightCellOfFace[ iFace ];
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            //! From left.
            pt = face2node[ nodepos + j ];
            q_n[ pt ] += q[ le ];
            ++ n_count[ pt ];

            //! From right.
            q_n[ pt ] += q[ re ];
            ++ n_count[ pt ];
        }
        nodepos += node_number_of_each_face[iFace];
    }

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    //! Reprocessing the points of interface.
    nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace [ iFace ];
        re = rightCellOfFace[ iFace ];

        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();
       
        if (IsInterface(bcType))
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt         = face2node[ nodepos + j ];
                q_n[ pt ] = 0.0;
                n_count[ pt ] = 0;
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

    nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace [ iFace ];
        re = rightCellOfFace[ iFace ];

        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        if (IsInterface(bcType))
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt         = face2node[ nodepos + j ];
                dx         = x[ pt ] - xcc[ le ];
                dy         = y[ pt ] - ycc[ le ];
                dz         = z[ pt ] - zcc[ le ];
                q_n[ pt ] += q[ le ];
                ++ n_count[ pt ];
                q_n[ pt ] += q[ re ];
                ++ n_count[ pt ];
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        q_n[iNode] /= (n_count[iNode] + SMALL);
    }

    delete [] n_count;    n_count = nullptr;
    delete [] dqdx;    dqdx = nullptr;
    delete [] dqdy;    dqdy = nullptr;
    delete [] dqdz;    dqdz = nullptr;
}

void CompNodeVarForVelocity(UnstructGrid *grid, RDouble *qNode, RDouble *q)
{
    int nTotalNode = grid->GetNTotalNode();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int * leftCellOfFace = grid->GetLeftCellOfFace();
    int * rightCellOfFace = grid->GetRightCellOfFace();
    int *face2node = grid->GetFace2Node();
    int *node_number_of_each_face  = grid->GetNodeNumberOfEachFace();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *weight = new RDouble [nTotalNode];

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        weight[iNode] = 0;
        qNode[iNode]  = 0.0;
    }

    int le, re, pt;

    //! Boundary faces.
    int nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace [ iFace ];
        re = rightCellOfFace[ iFace ];
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            int point = face2node[ nodepos + j ];
            RDouble dx = x[point] - xcc[le];
            RDouble dy = y[point] - ycc[le];
            RDouble dz = z[point] - zcc[le];
            RDouble dist = DISTANCE(dx, dy, dz);
            RDouble weightTemp = 1.0 / dist;
            qNode[ point ] += weightTemp * q[le];
            weight[point] += weightTemp;
            
            dx = x[point] - xcc[re];
            dy = y[point] - ycc[re];
            dz = z[point] - zcc[re];
            dist = DISTANCE(dx, dy, dz);
            weightTemp = 1.0 / dist;
            qNode[point] += weightTemp * q[re];
            weight[point] += weightTemp;
        }
        nodepos += node_number_of_each_face[iFace];
    }

    //! Interior faces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le  = leftCellOfFace [ iFace ];
        re  = rightCellOfFace[ iFace ];
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            int point = face2node[ nodepos + j ];

            //! From left.
            RDouble dx = x[point] - xcc[le];
            RDouble dy = y[point] - ycc[le];
            RDouble dz = z[point] - zcc[le];
            RDouble dist = DISTANCE(dx, dy, dz);
            RDouble weightTemp = 1.0 / dist;
            qNode[point] += weightTemp * q[le];
            weight[point] += weightTemp;

            //! From right.
            dx = x[point] - xcc[re];
            dy = y[point] - ycc[re];
            dz = z[point] - zcc[re];
            dist = DISTANCE(dx, dy, dz);
            weightTemp = 1.0 / dist;
            qNode[point] += weightTemp * q[re];
            weight[point] += weightTemp;
        }
        nodepos += node_number_of_each_face[iFace];
    }

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (weight[iNode] > SMALL)
        {
            qNode[iNode] /= weight[iNode];
        }        
    }

    //! Reset wall and symmetry faces.
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();
    nodepos = 0;

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        if(bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::SYMMETRY)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt = face2node[ nodepos + j ];

                qNode[ pt ] = 0.0;
            }
        }

        nodepos += node_number_of_each_face[iFace];
    }

    nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellOfFace [ iFace ];

        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        if(bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::SYMMETRY)
        {
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                pt = face2node[ nodepos + j ];

                qNode[ pt ] = q[ le ];
            }
        }
        nodepos += node_number_of_each_face[iFace];
    }

    delete [] weight;    weight = nullptr;
}

//! Bell 20131123 add.
void UnstructGrid::CompareMaxMinValue(RDouble &local_data, int flag)
{
    using namespace PHMPI;
    int numberOfProcessor =  GetNumberOfProcessor();
    int nZones = GetNumberofGlobalZones();
    if (nZones > numberOfProcessor)
    {
        return;
    }
    else
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int gridtype = GetZoneGridType(iZone);
            if (gridtype != PHSPACE::UNSTRUCTGRID)
            {
                return;
            }
        }
    }

    PH_CompareMaxMin(local_data, flag);
}

void UnstructGrid::ChangeBCType(int from_bctype, int to_bctype)
{
    UnstructBCSet **bcr = this->GetBCRecord();
    UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();
    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();

        if (bcType == from_bctype)
        {
            bcRegion->SetBCType(to_bctype);
            if (periodicType == NO_PERIODICITY)
            {
               bcRegion->SetBCName("INTERFACE");
            }
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter ++)
            {
                int faceIndexTmp = *iter;
                bcr[faceIndexTmp]->SetKey(to_bctype);
                if (periodicType == NO_PERIODICITY)
                {
                   bcr[faceIndexTmp]->SetBCName("INTERFACE");
                }
            }
        }
    }
}

int UnstructGrid::CompNIFace()
{
    UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    int nIFace = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (IsInterface(bcType))
        {
            vector<int>* faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                ++ nIFace;
            }
        }
    }

    return nIFace;
}

UnstructBCSet *UnstructGrid::GetUnstructBCSet()
{
    if (!unstructBCSet)
    {
        CreateUnstructBCSetInfo();
    }
    return unstructBCSet;
}

void UnstructGrid::CreateUnstructBCSetInfo()
{
    set<string> bcNameSet;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        bcNameSet.insert(bcr[iFace]->GetBCName());
    }

    int nBCRegionUnstruct = static_cast<int>(bcNameSet.size());
    this->CreateUnstructBCSet(nBCRegionUnstruct);

    UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
    int *bcRegionIDofBCFace = new int[nBoundFace];
    set<string>::iterator iter;
    int count = 0;
    for (iter = bcNameSet.begin(); iter != bcNameSet.end(); ++ iter)
    {
        UnstructBC *unstructBC = new UnstructBC(count);
        unstructBCSet->SetBCRegion(count, unstructBC);
        unstructBC->SetBCName(*iter);

        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            if (bcr[iFace]->GetBCName() == unstructBC->GetBCName())
            {
                unstructBC->SetBCType(bcr[iFace]->GetKey());
                bcRegionIDofBCFace[iFace] = count;
                vector<int> *faceIndex = unstructBC->GetFaceIndex();
                faceIndex->push_back(iFace);
            }
        }
        count ++;
    }
    unstructBCSet->SetBCFaceInfo(bcRegionIDofBCFace);
}

void UnstructGrid::CreateInterpolateInfo(ActionKey *actkey)
{
    int nDim = GetDim();
    int nTotalCell = this->GetNTotalCell();
    RDouble *xcc   = this->GetCellCenterX();
    RDouble *ycc   = this->GetCellCenterY();
    RDouble *zcc   = this->GetCellCenterZ();

    RDouble *pMin = this->GetMinBox();
    RDouble* pMax = this->GetMaxBox();

    RDouble *xyzCellCenter = new RDouble[3];
    RDouble *xyzCellCenterFromFile = new RDouble[3];

    RDouble *minDist   = new RDouble[nTotalCell];
    PHSPACE::SetField(minDist, LARGE, nTotalCell);

    int *targetZone    = new int[nTotalCell];
    int *targetCell    = new int[nTotalCell];

    this->UpdateDataPtr("targetZone", targetZone);
    this->UpdateDataPtr("targetCell", targetCell);

    hid_t gridData;
    string grpName;
    gridData = OpenGroup(actkey->gridfilepos, "Information");
    
    int nZones = 0;
    ReadData(gridData, &nZones, "nBlocks");

    int *zoneType = ReadIntOneRow(gridData, "block_type");

    for (int iZone = 0; iZone < nZones; iZone++)
    {
        hid_t gridDataLocal;
        string grpName;

        ostringstream oss;
        oss << "Grid-" << iZone;
        grpName = oss.str();

        gridDataLocal = OpenGroup(actkey->gridfilepos, grpName);

        RDouble *minMaxBox = new RDouble[6]();
        if (CheckDataExist(gridData,"minMaxBox"))
        {
            RDouble *gridBox = ReadDoubleOneRow(gridData, "MinMaxBoX");
        }
        else
        {
            GridID *gridIndex = new GridID(iZone);
            if (zoneType[iZone] == UNSTRUCTGRID)
            {
                UnstructGrid *gridLocal = new UnstructGrid();
                gridLocal->InitGrid(gridIndex, 0, nDim, UNSTRUCTGRID);
               
                int *iDimensions = ReadIntOneRow(gridDataLocal, "iDimensions");
                gridLocal->SetNTotalNode(iDimensions[0]);
                gridLocal->SetNTotalFace(iDimensions[1]);
                gridLocal->SetNTotalCell(iDimensions[2]);

                gridData = OpenGroup(gridDataLocal, "GridCoordinates");
                RDouble *coordinateX = ReadDoubleOneRow(gridData, "CoordinateX");
                RDouble *coordinateY = ReadDoubleOneRow(gridData, "CoordinateY");
                RDouble *coordinateZ = ReadDoubleOneRow(gridData, "CoordinateZ");

                gridLocal->SetX(coordinateX);
                gridLocal->SetY(coordinateY);
                gridLocal->SetZ(coordinateZ);

                RDouble *pMinLocal = gridLocal->GetMinBox();
                RDouble *pMaxLocal = gridLocal->GetMaxBox();

                for (int iDim = 0; iDim < nDim; iDim++)
                {
                    minMaxBox[iDim] = pMinLocal[iDim];
                    minMaxBox[iDim + 3] = pMaxLocal[iDim];
                }
                delete gridLocal;
            }
            else
            {
                StructGrid *gridLocal = new StructGrid();
                gridLocal->InitGrid(gridIndex, 0, nDim, STRUCTGRID);

                int *iDimensions = ReadIntOneRow(gridDataLocal, "iDimensions");
                gridLocal->SetNTotalNode(iDimensions[0]);
                gridLocal->SetNTotalFace(iDimensions[1]);
                gridLocal->SetNTotalCell(iDimensions[2]);

                gridData = OpenGroup(gridDataLocal, "GridCoordinates");
                RDouble *coordinateX = ReadDoubleOneRow(gridData, "CoordinateX");
                RDouble *coordinateY = ReadDoubleOneRow(gridData, "CoordinateY");
                RDouble *coordinateZ = ReadDoubleOneRow(gridData, "CoordinateZ");

                gridLocal->SetX(coordinateX);
                gridLocal->SetY(coordinateY);
                gridLocal->SetZ(coordinateZ);

                RDouble *pMinLocal = gridLocal->GetMinBox();
                RDouble *pMaxLocal = gridLocal->GetMaxBox();

                for (int iDim = 0; iDim < nDim; iDim++)
                {
                    minMaxBox[iDim] = pMinLocal[iDim];
                    minMaxBox[iDim + 3] = pMaxLocal[iDim];
                }
                delete gridLocal;
            }
        }

        if (IfBoxOverset(&minMaxBox[0], &minMaxBox[3], nDim))
        {
            KDTree* interpolateKDTree;
            GridID *gridIndex = new GridID(iZone);
            interpolateKDTree = CreatKDTree(nDim);

            if (zoneType[iZone] == UNSTRUCTGRID)
            {
                UnstructGrid *gridLocal = new UnstructGrid();
                gridLocal->InitGrid(gridIndex, 0, nDim, UNSTRUCTGRID);
              
                int *iDimensions = ReadIntOneRow(gridDataLocal, "iDimensions");
                gridLocal->SetNTotalNode(iDimensions[0]);
                gridLocal->SetNTotalFace(iDimensions[1]);
                gridLocal->SetNTotalCell(iDimensions[2]);

                gridData = OpenGroup(gridDataLocal, "GridCoordinates");
                RDouble *coordinateX = ReadDoubleOneRow(gridData, "CoordinateX");
                RDouble *coordinateY = ReadDoubleOneRow(gridData, "CoordinateY");
                RDouble *coordinateZ = ReadDoubleOneRow(gridData, "CoordinateZ");

                gridLocal->SetX(coordinateX);
                gridLocal->SetY(coordinateY);
                gridLocal->SetZ(coordinateZ);

                int nBoundFace = 0;
                hid_t faceData;

                faceData = OpenGroup(gridDataLocal, "FaceTopology");
                ReadData(faceData, &nBoundFace, "nBoundFace");
                gridLocal->SetNBoundFace(nBoundFace);

                int *nodeNumberofEachFace = ReadIntOneRow(faceData, "nodeNumberOfEachFace");
                gridLocal->SetNodeNumberOfEachFace(nodeNumberofEachFace);

                int *face2Node = ReadIntOneRow(faceData, "face2Node");
                gridLocal->SetFace2Node(face2Node);

                int *leftCellOfFace = ReadIntOneRow(faceData, "leftCellOfFace");
                gridLocal->SetLeftCellOfFace(leftCellOfFace);

                int *rightCellOfFace = ReadIntOneRow(faceData, "rightCellOfFace");
                gridLocal->SetRightCellOfFace(rightCellOfFace);

                int nTotalFace = gridLocal->GetNTotalFace();
                int *nodePosi = new int[nTotalFace + 1];
                int nodeSum = 0; 
                nodePosi[0] = 0;
                for (int iFace = 0; iFace < nTotalFace; iFace++)
                {
                    nodeSum += nodeNumberofEachFace[iFace];
                    nodePosi[iFace + 1] = nodeSum;

                }

                for (int iFace = 0; iFace < nTotalFace; iFace++)
                {
                    if (leftCellOfFace[iFace] < 0)
                    {
                        std::reverse(face2Node + nodePosi[iFace], face2Node + nodePosi[iFace + 1]);
                        SWAP(leftCellOfFace[iFace], rightCellOfFace[iFace]);
                    }

                }
                delete []nodePosi;
                H5Gclose(faceData);

                gridLocal->SpecifyRightCellofBC();
                gridLocal->ComputeMetrics();

                RDouble *xccLocal = gridLocal->GetCellCenterX();
                RDouble *yccLocal = gridLocal->GetCellCenterY();
                RDouble *zccLocal = gridLocal->GetCellCenterZ();

                int nTotalNode = gridLocal->GetNTotalNode();
                int nTotalCell = gridLocal->GetNTotalCell();

                for (int iCell = 0; iCell < nTotalCell; iCell++)
                {
                    
                    xyzCellCenterFromFile[0] = xccLocal[iCell];
                    xyzCellCenterFromFile[1] = yccLocal[iCell];
                    if (nDim == 3)
                    {
                        xyzCellCenterFromFile[2] = zccLocal[iCell];
                    }

                    KDInsert(interpolateKDTree, xyzCellCenterFromFile, NULL);
                }
                delete gridLocal;
            }

            for (int iCell = 0; iCell < nTotalCell; iCell++)
            {
                xyzCellCenter[0] = xcc[iCell];
                xyzCellCenter[1] = ycc[iCell];
                if (nDim == 3)
                {
                    xyzCellCenter[2] = zcc[iCell];
                }
                KDRes *minResults = KDNearest(interpolateKDTree, xyzCellCenter);
                RDouble dist = sqrt(DistSQ(minResults->Position(), xyzCellCenter, nDim));
                int targetCellInThisZone = minResults->itemData();
                if (dist < minDist[iCell])
                {
                    minDist[iCell]    = dist;
                    targetZone[iCell] = iZone;
                    targetCell[iCell] = targetCellInThisZone;
                }
                FreeKDRes(minResults);
            }

            interpolateKDTree->Free();
            delete interpolateKDTree;
            interpolateKDTree = nullptr;
        }
        delete []minMaxBox;

        H5Gclose(gridDataLocal);
        H5Gclose(gridData);
    }

    delete []minDist;
    delete []xyzCellCenter;
    delete []xyzCellCenterFromFile;
}

void UnstructGrid::BuildGridLink(LinkStruct * link)
{
    cout << "Interface Construction ..." << endl;

    using namespace PHSPACE;
    int zoneIndex = this->GetZoneID();

    int nBoundFace = this->GetNBoundFace();

    UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    int *face2node = this->GetFace2Node();
    int *node_number_of_each_face = this->GetNodeNumberOfEachFace();

    int nIFace = this->CompNIFace();
    cout << "Number of Inter Faces: " << nIFace << "\n";

    InterfaceInfo *interfaceInfo = 0;
    if (nIFace == 0)
    {
        this->SetInterfaceInfo(interfaceInfo);
        return;
    }

    interfaceInfo = new InterfaceInfo(nIFace);
    this->SetInterfaceInfo(interfaceInfo);

    VVInt &facemap = link->GetFaceMap();
    facemap[zoneIndex].resize(nIFace);

    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();

    //! Find out the maximum value of node_number_of_each_face on the boundary face.
    int maxlist = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        maxlist = MAX(maxlist, node_number_of_each_face[iFace]);
    }

    RDouble *xlist = new RDouble[maxlist];
    RDouble *ylist = new RDouble[maxlist];
    RDouble *zlist = new RDouble[maxlist];

    vector<int> index(maxlist),pindex(maxlist);

    LinkStruct::AdtTree &coor_tree = link->GetCoordinateTree();
    set < DataStruct_Sort< VInt > > &facelist = link->GetFaceList();
    VVInt &zoneid = link->GetZoneID();
    VVInt &faceid = link->GetFaceID();
    RDouble tol   = link->GetTolerance();

    int iFace = 0;
    uint_t fcount = zoneid.size();
    int pcount = coor_tree.GetNodeNum();
    int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");
    RDouble translationLength[3] = { 0.0 };
    GlobalDataBase::GetData("translationLength", &translationLength, PHDOUBLE, 3);
    RDouble rotationAngle = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("rotationAngle");
    rotationAngle = rotationAngle * PI / 180.0;

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");

    int nodepos = 0;
    for (int iBFace = 0; iBFace < nBoundFace; ++ iBFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iBFace]);
        int bcType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();

        if (!IsInterface(bcType))
        {
            nodepos += node_number_of_each_face[iBFace];
            continue;
        }
        for (int m = 0; m < node_number_of_each_face[iBFace]; ++ m)
        {
            index[m] = face2node[ nodepos + m ];
        }
        nodepos += node_number_of_each_face[iBFace];

        int nlist = node_number_of_each_face[iBFace];
        pindex.resize(nlist);

        GetFaceCoorList(index, nlist, xlist, ylist, zlist, x, y, z);

        if (referenceFrame == ROTATIONAL_FRAME)
        {
            int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
            RDouble PeriodicRotationAngle[100];
            GlobalDataBase::GetData("PeriodicRotationAngle", &PeriodicRotationAngle, PHDOUBLE, nTurboZone);
            string Periodic_Name[100];
            GlobalDataBase::GetData("Periodic_Name", &Periodic_Name, PHSTRING, 2 * nTurboZone);

            for (int iTurboZone = 0; iTurboZone < nTurboZone; iTurboZone++)
            {
                if (bcName == Periodic_Name[2 * iTurboZone])
                {
                    PeriodicRotationAngle[iTurboZone] = PeriodicRotationAngle[iTurboZone] * PI / 180;
                    rotationAngle = PeriodicRotationAngle[iTurboZone];

                    for (int m = 0; m < node_number_of_each_face[iBFace]; ++m)
                    {
                        RDouble rotPoint[3] = { 0.0, 0.0, 0.0 };
                        rotPoint[0] = xlist[m];
                        rotPoint[1] = ylist[m] * cos(rotationAngle) - zlist[m] * sin(rotationAngle);
                        rotPoint[2] = ylist[m] * sin(rotationAngle) + zlist[m] * cos(rotationAngle);
                        xlist[m] = rotPoint[0];
                        ylist[m] = rotPoint[1];
                        zlist[m] = rotPoint[2];
                    }
                }
            }
        }
        else
        {
            if (bcName == "Periodic_up")
            {
                if (periodicType == TRANSLATIONAL_PERIODICITY)
                {
                    for (int m = 0; m < node_number_of_each_face[iBFace]; ++m)
                    {
                        xlist[m] += translationLength[0];
                        ylist[m] += translationLength[1];
                        zlist[m] += translationLength[2];
                    }
                }
                else if (periodicType == ROTATIONAL_PERIODICITY)
                {
                    for (int m = 0; m < node_number_of_each_face[iBFace]; ++m)
                    {
                        RDouble rotPoint[3] = { 0.0, 0.0, 0.0 };
                        rotPoint[0] = xlist[m];
                        rotPoint[1] = ylist[m] * cos(rotationAngle) - zlist[m] * sin(rotationAngle);
                        rotPoint[2] = ylist[m] * sin(rotationAngle) + zlist[m] * cos(rotationAngle);
                        xlist[m] = rotPoint[0];
                        ylist[m] = rotPoint[1];
                        zlist[m] = rotPoint[2];
                    }
                }
            }
        }
        GetCoorIndexList(&coor_tree, tol, pcount, xlist, ylist, zlist, nlist, pindex);
        Create_Link_Info(pindex, zoneIndex, iFace, fcount, facelist, zoneid, faceid, link);

        interFace2BoundaryFace[iFace] = iBFace;
        ++ iFace;
    }
    cout << "iFace = " << iFace << "\n";
    delete [] xlist;    xlist = nullptr;
    delete [] ylist;    ylist = nullptr;
    delete [] zlist;    zlist = nullptr;
}

//! Compute fcptr (There are two cells for each cell face. Which cell 
//! number is a cell in its neighboring cell).
void UnstructGrid::CalCNNCF()
{
    int nTotalFace = this->GetNTotalFace();
    //! If fcptr has already existed.
    int *fc2cL = this->Getfc2cL();
    int *fc2cR = this->Getfc2cR();
    if (fc2cL) return;
    
    int nBoundFace = this->GetNBoundFace();
    int nTotalCell = this->GetNTotalCell();
    int *leftCellOfFace  = this->GetLeftCellOfFace();
    int *rightCellOfFace = this->GetRightCellOfFace();

    UnstructBCSet *unstructBCSet = this->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    int *number_of_neighbor_cell = this->GetNumberOfNeighborCell();

    //! Allocate memories for number of cells per cell.
    fc2cL = new int[nTotalFace];
    fc2cR = new int[nTotalFace];

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        fc2cL[iFace] = - 1;
        fc2cR[iFace] = - 1;
    }

    int le, re;

    //! Need to reset number_of_neighbor_cell to zero and recover it later.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        number_of_neighbor_cell[iCell] = 0;
    }

    //! If boundary is an INTERFACE, need to count ghost cell.
    //! Note: Symmetry???
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (IsInterface(bcType))
        {
            vector<int>* faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                le = leftCellOfFace[iFace];
                number_of_neighbor_cell[le] ++;
                fc2cL[iFace] = number_of_neighbor_cell[le];
            }
        }
    }

    //! Interior faces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellOfFace [ iFace ];
        re = rightCellOfFace[ iFace ];

        number_of_neighbor_cell[le] ++;
        number_of_neighbor_cell[re] ++;

        fc2cL[iFace] = number_of_neighbor_cell[le];
        fc2cR[iFace] = number_of_neighbor_cell[re];
    }

    this->Setfc2cL(fc2cL);
    this->Setfc2cR(fc2cR);
}

int UnstructGrid::GetNumberOfWallCell()
{
    UnstructBCSet **bcr = this->GetBCRecord();

    int nBoundFace = this->GetNBoundFace();

    int nsolid_surface = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bctype = bcr[iFace]->GetKey();
        if (IsWall(bctype))
        {
            ++ nsolid_surface;
        }
    }

    return nsolid_surface;
}

void UnstructGrid::ProcessBCInfo()
{

}

void UnstructGrid::GridSurfaceVelocity(RDouble *xt, RDouble *yt, RDouble *zt)
{

}

void UnstructGrid::ReSetBoundaryConditionByGlobalBC()
{
    vector<SimpleBC *> *globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();

    map<string, SimpleBC *> globalBCMap;

    vector<SimpleBC *>::iterator iter;
    for (iter = globalBCList->begin(); iter != globalBCList->end(); ++ iter)
    {
        SimpleBC *bc = *iter;
        pair<string, SimpleBC *> bcPair(bc->GetBCName(), bc);
        globalBCMap.insert(bcPair);
    }

    int nBoundFace = this->GetNBoundFace();

    UnstructBCSet **bcr = this->GetBCRecord();
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        const string &bcName = bcr[iFace]->GetBCName();
        SimpleBC *bc = 0;
        map<string, SimpleBC *>::iterator iter1 = globalBCMap.find(bcName);

        if (iter1 != globalBCMap.end())
        {
            bc = iter1->second;
        }

        bcr[iFace]->SetBoundaryCondition(bc);
    }
}

void GetFace2NodeList(UnstructGrid *grid, int bcType, vector<int> &linkmap, vector < vector < int > > &face2nodelist)
{
    int nTotalNode = grid->GetNTotalNode();
    int nBoundFace = grid->GetNBoundFace();

    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();
    int *face2node = grid->GetFace2Node();

    //! Build point set of binary tree for searching.
    vector < int > loc_face2node(4);
    set< DataStruct_Sort<int> > iset;
    set< DataStruct_Sort<int> >::iterator ifind;
    DataStruct_Sort<int> fysort;

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    face2nodelist.resize(0);
    linkmap.resize(0);
    int nodepos = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);

        if (bcType != bcRegion->GetBCType())
        {
            nodepos += node_number_of_each_face[iFace];
            continue;
        }

        if (node_number_of_each_face[iFace] <= 4)
        {
            loc_face2node.resize(0);
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                int index = face2node[nodepos + j];

                fysort.value = index;
                ifind = iset.find(fysort);

                if (ifind == iset.end())
                {
                    linkmap.push_back(index);
                    fysort.index = static_cast<int>(linkmap.size() - 1);

                    iset.insert(fysort);
                    loc_face2node.push_back(fysort.index);
                }
                else
                {
                    loc_face2node.push_back(ifind->index);
                }

            }
            nodepos += node_number_of_each_face[iFace];
            face2nodelist.push_back(loc_face2node);
        }
        else
        {
            int nodesWrite = 0;
            while (nodesWrite < node_number_of_each_face[iFace])
            {
                loc_face2node.resize(0);

                //! Write face center.
                int index = iFace + nTotalNode;

                fysort.value = index;
                ifind = iset.find(fysort);

                if (ifind == iset.end())
                {
                    linkmap.push_back(index);
                    fysort.index = static_cast<int>(linkmap.size() - 1);

                    iset.insert(fysort);
                    loc_face2node.push_back(fysort.index);
                }
                else
                {
                    loc_face2node.push_back(ifind->index);
                }

                //! Write node.
                for (int iNode = 0; iNode < 2; ++ iNode)
                {
                    int index1 = (iNode + nodesWrite) % node_number_of_each_face[iFace];
                    int nodeIndex = face2node[nodepos + index1];

                    fysort.value = nodeIndex;
                    ifind = iset.find(fysort);

                    if (ifind == iset.end())
                    {
                        linkmap.push_back(nodeIndex);
                        fysort.index = static_cast<int>(linkmap.size() - 1);

                        iset.insert(fysort);
                        loc_face2node.push_back(fysort.index);
                    }
                    else
                    {
                        loc_face2node.push_back(ifind->index);
                    }
                }

                ++ nodesWrite;
                face2nodelist.push_back(loc_face2node);
            }
            nodepos += node_number_of_each_face[iFace];
        }
    }
}
void GetFace2NodeListForLine(UnstructGrid* grid, pair<int, string> iter, vector<int>& linkmap, vector < vector < int > >& face2nodelist, set<int>& line)
{
    int bcType = (iter).first;
    string bcName = (iter).second;

    int nTotalNode = grid->GetNTotalNode();
    int nBoundFace = grid->GetNBoundFace();

    int * node_number_of_each_face = grid->GetNodeNumberOfEachFace();
    int * face2node = grid->GetFace2Node();

    //! Build point set of binary tree for searching.
    vector < int > loc_face2node(4);
    set< DataStruct_Sort<int> > iset;
    set< DataStruct_Sort<int> >::iterator ifind;
    DataStruct_Sort<int> fysort;

    UnstructBCSet **bcr = grid->GetBCRecord();

    face2nodelist.resize(0);
    linkmap.resize(0);

    int * leftCellIndexContainer = grid->GetLeftCellOfFace();
    int * keyActiveOfCells       = grid->GetBlankIndex();
    int nodepos = 0;

    //for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    //{
    //    if (bcType != bcr[iFace]->GetKey())
    //    {
    //        nodepos += node_number_of_each_face[iFace];
    //        continue;
    //    }

    //    if (bcr[iFace]->GetBCName() != bcName)
    //    {
    //        nodepos += node_number_of_each_face[iFace];
    //        continue;
    //    }

    //    if (keyActiveOfCells[leftCellIndexContainer[iFace]] != 1)
    //    {
    //        nodepos += node_number_of_each_face[iFace];
    //        continue;
    //    }

    //    if (node_number_of_each_face[iFace] <= 4)
    //    {
    //        loc_face2node.resize(0);
    //        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
    //        {
    //            int index = face2node[ nodepos + j ];

    //            fysort.value = index;
    //            ifind = iset.find(fysort);

    //            if (ifind == iset.end())
    //            {
    //                linkmap.push_back(index);
    //                fysort.index = static_cast<int>(linkmap.size() - 1);

    //                iset.insert(fysort);
    //                loc_face2node.push_back(fysort.index);
    //            }
    //            else
    //            {
    //                loc_face2node.push_back(ifind->index);
    //            }

    //        }
    //        nodepos += node_number_of_each_face[iFace];
    //        face2nodelist.push_back(loc_face2node);
    //    }
    //    else
    //    {
    //        int nodesWrite = 0;
    //        while (nodesWrite < node_number_of_each_face[iFace])
    //        {
    //            loc_face2node.resize(0);

    //            //! Write face center.
    //            int index = iFace + nTotalNode;

    //            fysort.value = index;
    //            ifind = iset.find(fysort);

    //            if (ifind == iset.end())
    //            {
    //                linkmap.push_back(index);
    //                fysort.index = static_cast<int>(linkmap.size() - 1);

    //                iset.insert(fysort);
    //                loc_face2node.push_back(fysort.index);
    //            }
    //            else
    //            {
    //                loc_face2node.push_back(ifind->index);
    //            }

    //            //! Write node.
    //            for (int iNode = 0; iNode < 2; ++ iNode)
    //            {
    //                int index = (iNode + nodesWrite) % node_number_of_each_face[iFace];
    //                int nodeIndex = face2node[nodepos + index];

    //                fysort.value = nodeIndex;
    //                ifind = iset.find(fysort);

    //                if (ifind == iset.end())
    //                {
    //                    linkmap.push_back(nodeIndex);
    //                    fysort.index = static_cast<int>(linkmap.size() - 1);

    //                    iset.insert(fysort);
    //                    loc_face2node.push_back(fysort.index);
    //                }
    //                else
    //                {
    //                    loc_face2node.push_back(ifind->index);
    //                }
    //            }

    //            ++ nodesWrite;
    //            face2nodelist.push_back(loc_face2node);
    //        }
    //        nodepos += node_number_of_each_face[iFace];
    //    }
    //}

    int** face2nodearrary = grid->GetFace2NodeArray();

    for (set<int>::iterator iter = line.begin(); iter != line.end(); ++iter)
    {
        int iFace = *iter;
        if (node_number_of_each_face[iFace] <= 4)
        {
            loc_face2node.resize(0);
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {

                int index = face2nodearrary[iFace][j];

                fysort.value = index;
                ifind = iset.find(fysort);

                if (ifind == iset.end())
                {
                    linkmap.push_back(index);
                    fysort.index = static_cast<int>(linkmap.size() - 1);

                    iset.insert(fysort);
                    loc_face2node.push_back(fysort.index);
                }
                else
                {
                    loc_face2node.push_back(ifind->index);
                }

            }
            face2nodelist.push_back(loc_face2node);
        }
        else
        {
            int nodesWrite = 0;
            while (nodesWrite < node_number_of_each_face[iFace])
            {
                loc_face2node.resize(0);

                //! Write face center.
                int index = iFace + nTotalNode;

                fysort.value = index;
                ifind = iset.find(fysort);

                if (ifind == iset.end())
                {
                    linkmap.push_back(index);
                    fysort.index = static_cast<int>(linkmap.size() - 1);

                    iset.insert(fysort);
                    loc_face2node.push_back(fysort.index);
                }
                else
                {
                    loc_face2node.push_back(ifind->index);
                }
                //! Write node.
                for (int iNode = 0; iNode < 2; ++ iNode)
                {
                    int index = (iNode + nodesWrite) % node_number_of_each_face[iFace];
                    int nodeIndex = face2nodearrary[iFace][index];
                    fysort.value = nodeIndex;
                    ifind = iset.find(fysort);

                    if (ifind == iset.end())
                    {
                        linkmap.push_back(nodeIndex);
                        fysort.index = static_cast<int>(linkmap.size() - 1);

                        iset.insert(fysort);
                        loc_face2node.push_back(fysort.index);
                    }
                    else
                    {
                        loc_face2node.push_back(ifind->index);
                    }
                }
                ++ nodesWrite;
                face2nodelist.push_back(loc_face2node);
            }
            nodepos += node_number_of_each_face[iFace];
        }
    }
}
void GetFace2NodeList(UnstructGrid *grid, pair<int, string> iter, vector<int> &linkmap, vector < vector < int > > &face2nodelist)
{
    int bcType = (iter).first;
    string bcName = (iter).second;

    int nTotalNode = grid->GetNTotalNode();
    int nBoundFace = grid->GetNBoundFace();

    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();
    int *face2node = grid->GetFace2Node();

    //! Build point set of binary tree for searching.
    vector < int > loc_face2node(4);
    set< DataStruct_Sort<int> > iset;
    set< DataStruct_Sort<int> >::iterator ifind;
    DataStruct_Sort<int> fysort;

    UnstructBCSet **bcr = grid->GetBCRecord();

    face2nodelist.resize(0);
    linkmap.resize(0);

    int *leftCellIndexContainer = grid->GetLeftCellOfFace();
    int *keyActiveOfCells       = grid->GetBlankIndex();
    int nodepos = 0;

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        if (bcType != bcr[iFace]->GetKey())
        {
            nodepos += node_number_of_each_face[iFace];
            continue;
        }

        if (bcr[iFace]->GetBCName() != bcName)
        {
            nodepos += node_number_of_each_face[iFace];
            continue;
        }

        if (keyActiveOfCells[leftCellIndexContainer[iFace]] != ACTIVE)
        {
            nodepos += node_number_of_each_face[iFace];
            continue;
        }

        if (node_number_of_each_face[iFace] <= 4)
        {
            loc_face2node.resize(0);
            for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
            {
                int index = face2node[ nodepos + j ];

                fysort.value = index;
                ifind = iset.find(fysort);

                if (ifind == iset.end())
                {
                    linkmap.push_back(index);
                    fysort.index = static_cast<int>(linkmap.size() - 1);

                    iset.insert(fysort);
                    loc_face2node.push_back(fysort.index);
                }
                else
                {
                    loc_face2node.push_back(ifind->index);
                }

            }
            nodepos += node_number_of_each_face[iFace];
            face2nodelist.push_back(loc_face2node);
        }
        else
        {
            int nodesWrite = 0;
            while (nodesWrite < node_number_of_each_face[iFace])
            {
                loc_face2node.resize(0);

                //! Write face center.
                int index = iFace + nTotalNode;

                fysort.value = index;
                ifind = iset.find(fysort);

                if (ifind == iset.end())
                {
                    linkmap.push_back(index);
                    fysort.index = static_cast<int>(linkmap.size() - 1);

                    iset.insert(fysort);
                    loc_face2node.push_back(fysort.index);
                }
                else
                {
                    loc_face2node.push_back(ifind->index);
                }

                //! Write node.
                for (int iNode = 0; iNode < 2; ++ iNode)
                {
                    int index1 = (iNode + nodesWrite) % node_number_of_each_face[iFace];
                    int nodeIndex = face2node[nodepos + index1];

                    fysort.value = nodeIndex;
                    ifind = iset.find(fysort);

                    if (ifind == iset.end())
                    {
                        linkmap.push_back(nodeIndex);
                        fysort.index = static_cast<int>(linkmap.size() - 1);

                        iset.insert(fysort);
                        loc_face2node.push_back(fysort.index);
                    }
                    else
                    {
                        loc_face2node.push_back(ifind->index);
                    }
                }

                ++ nodesWrite;
                face2nodelist.push_back(loc_face2node);
            }
            nodepos += node_number_of_each_face[iFace];
        }
    }
}

void MinMaxDiffINF(Grid *grid_in, int ig, RDouble *dmin, RDouble *dmax, RDouble *q)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    int **cell2face  = grid->GetCell2Face();
    int *face_number_of_each_cell = grid->GetFaceNumberOfEachCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    int nIFace = interfaceInfo->GetNIFace();
    //! Find the maximum and minimum in the neighbor of each cell.
    for (int iFace = 0; iFace < nIFace; ++ iFace)
    {
        int ie,face,le,re;
        grid->GetSourceIndex(iFace,ig+1,ie);
        dmin[iFace] = q[ie];
        dmax[iFace] = q[ie];
        for (int jFace = 0; jFace < face_number_of_each_cell[ie]; ++ jFace)
        {
            face = cell2face[ie][jFace];
            le   = leftCellOfFace[face];
            re   = rightCellOfFace[face];
            dmin[iFace] = MIN(dmin[iFace], q[le]);
            dmin[iFace] = MIN(dmin[iFace], q[re]);

            dmax[iFace] = MAX(dmax[iFace], q[le]);
            dmax[iFace] = MAX(dmax[iFace], q[re]);
        }
        dmin[iFace] -= q[ie];
        dmax[iFace] -= q[ie];
    }
}

void SetLimitBoundary(Grid *grid_in, RDouble *limit)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType != PHENGLEI::INTERFACE && 
            bcType != PHENGLEI::SYMMETRY  &&
            !IsWall(bcType) &&
            bcType != PHENGLEI::OUTFLOW &&
            bcType != PHENGLEI::OVERSET)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                int le = leftCellOfFace [iFace];
                int re = rightCellOfFace[iFace];

                limit[le] = zero;
                limit[re] = zero;
            }
        }
    }
}

void CommunicateInterpointValue(UnstructGrid *grid, RDouble **variable, const string &variableName, int variableDimension)
{
    PHSPACE::UploadInterpointValue(grid, variable, variableName, variableDimension);

    PHSPACE::DownloadInterpointValue(grid, variable, variableName, variableDimension);
}

void CommunicateInterfaceValue(UnstructGrid * grid, RDouble ** variable, const string & variableName, int variableDimension)
{
    //! Upload the interface data to the buffer.
    PHSPACE::UploadInterfaceValue(grid, variable, variableName, variableDimension);

    //! Download the interface data from the buffer.
    PHSPACE::DownloadInterfaceValue(grid, variable, variableName, variableDimension);
}

void CommunicateInterfaceValue(UnstructGrid *grid, RDouble *variable, const string &variableName)
{
    //! Upload the interface data to the buffer.
    PHSPACE::UploadInterfaceValue(grid, variable, variableName);

    //PHSPACE::UpdateInterfaceData();

    //! Download the interface data from the buffer.
    PHSPACE::DownloadInterfaceValue(grid, variable, variableName);
}

void UploadInterfaceValue(UnstructGrid *grid, RDouble **f, const string &name, int neqn)
{
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    int nIFace = interfaceInfo->GetNIFace();
    for (int iGhost = GetNumberOfGhostCellLayers() - 1; iGhost >= 0; -- iGhost)
    {
        Data_ParamFieldSuite *sendData = interfaceInfo->GetSendDataStorage(iGhost);

        RDouble **fg = reinterpret_cast< RDouble ** > (sendData->GetDataPtr(name));
        if(!fg)
        {
            ostringstream oss;
            oss << "    Error: Parallel interface variable has not been register, name: " << name;
            TK_Exit::ExceptionExit(oss.str());
        }
        for (int iFace = 0; iFace < nIFace; ++ iFace)
        {
            int sourceCell;
            grid->GetSourceIndex(iFace, iGhost + 1, sourceCell);
            for (int m = 0; m < neqn; ++ m)
            {
                fg[m][iFace] = f[m][sourceCell];
            }
        }
    }
}

void UploadInterfaceValue(UnstructGrid *grid, RDouble *f, const string &name)
{
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo)
    {
        return;
    }

    int nIFace = interfaceInfo->GetNIFace();
    for (int iGhost = GetNumberOfGhostCellLayers() - 1; iGhost >= 0; -- iGhost)
    {
        Data_ParamFieldSuite *sendData = interfaceInfo->GetSendDataStorage(iGhost);

        RDouble **fg = reinterpret_cast <RDouble **> (sendData->GetDataPtr(name));
        if(!fg)
        {
            ostringstream oss;
            oss << "    Error: Parallel interface variable has not been register, name: " << name;
            TK_Exit::ExceptionExit(oss.str());
        }
        for (int iFace = 0; iFace < nIFace; ++ iFace)
        {
            int sourceCell;
            grid->GetSourceIndex(iFace, iGhost + 1, sourceCell);
            fg[0][iFace] = f[sourceCell];
        }
    }
}

void DownloadInterfaceValue(UnstructGrid *grid, RDouble **f, const string &name, int neqn)
{
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    int nIFace = interfaceInfo->GetNIFace();

    int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");
    double rotationAngle = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("rotationAngle");
    rotationAngle = rotationAngle / 180.0 * PI;
    int nBoundFace = grid->GetNBoundFace();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");

    if (referenceFrame == ROTATIONAL_FRAME)
    {
        int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
        RDouble PeriodicRotationAngle[100];
        GlobalDataBase::GetData("PeriodicRotationAngle", &PeriodicRotationAngle, PHDOUBLE, nTurboZone);

        //! Parallel
        int iTurboZone = grid->GetOrdinaryGridIndex();
        //! Serial
        if (iTurboZone == -1)
        {
            iTurboZone = grid->GetZoneID();
        }

        PeriodicRotationAngle[iTurboZone] = PeriodicRotationAngle[iTurboZone] * PI / 180;
        rotationAngle = PeriodicRotationAngle[iTurboZone];
    }

    for (int iGhost = GetNumberOfGhostCellLayers() - 1; iGhost >= 0; -- iGhost)
    {
        Data_ParamFieldSuite *recieveData = interfaceInfo->GetRecvDataStorage(iGhost);

        RDouble **fg = reinterpret_cast< RDouble ** > (recieveData->GetDataPtr(name));

        for (int iFace = 0; iFace < nIFace; ++ iFace)
        {
            int targetCell;
            grid->GetTargetIndex(iFace, iGhost + 1, targetCell);

            int iBFace = interFace2BoundaryFace[iFace];
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iBFace]);
            string bcName = bcRegion->GetBCName();

            if (neqn >= 4 && periodicType == ROTATIONAL_PERIODICITY)
            {
                RDouble rotfg[3] = { 0.0, 0.0, 0.0 };

                if (referenceFrame == ROTATIONAL_FRAME)
                {
                    int nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
                    string Periodic_Name[100];
                    GlobalDataBase::GetData("Periodic_Name", &Periodic_Name, PHSTRING, 2 * nTurboZone);

                    //! Parallel
                    int iTurboZone = grid->GetOrdinaryGridIndex();
                    //! Serial
                    if (iTurboZone == -1)
                    {
                        iTurboZone = grid->GetZoneID();
                    }

                    if (bcName == Periodic_Name[2 * iTurboZone])
                    {
                        rotfg[0] = fg[1][iFace];
                        rotfg[1] = fg[2][iFace] * cos(2.0 * PI - rotationAngle) - fg[3][iFace] * sin(2.0 * PI - rotationAngle);
                        rotfg[2] = fg[2][iFace] * sin(2.0 * PI - rotationAngle) + fg[3][iFace] * cos(2.0 * PI - rotationAngle);

                        fg[1][iFace] = rotfg[0];
                        fg[2][iFace] = rotfg[1];
                        fg[3][iFace] = rotfg[2];
                    }
                    if (bcName == Periodic_Name[2 * iTurboZone + 1])
                    {
                        rotfg[0] = fg[1][iFace];
                        rotfg[1] = fg[2][iFace] * cos(rotationAngle) - fg[3][iFace] * sin(rotationAngle);
                        rotfg[2] = fg[2][iFace] * sin(rotationAngle) + fg[3][iFace] * cos(rotationAngle);

                        fg[1][iFace] = rotfg[0];
                        fg[2][iFace] = rotfg[1];
                        fg[3][iFace] = rotfg[2];
                    }
                }
                else
                {
                    if (bcName == "Periodic_up" || bcName == "Periodic_down")
                    {
                        if (bcName == "Periodic_up")
                        {
                            rotfg[0] = fg[1][iFace];
                            rotfg[1] = fg[2][iFace] * cos(2.0 * PI - rotationAngle) - fg[3][iFace] * sin(2.0 * PI - rotationAngle);
                            rotfg[2] = fg[2][iFace] * sin(2.0 * PI - rotationAngle) + fg[3][iFace] * cos(2.0 * PI - rotationAngle);
                        }
                        if (bcName == "Periodic_down")
                        {
                            rotfg[0] = fg[1][iFace];
                            rotfg[1] = fg[2][iFace] * cos(rotationAngle) - fg[3][iFace] * sin(rotationAngle);
                            rotfg[2] = fg[2][iFace] * sin(rotationAngle) + fg[3][iFace] * cos(rotationAngle);
                        }
                        fg[1][iFace] = rotfg[0];
                        fg[2][iFace] = rotfg[1];
                        fg[3][iFace] = rotfg[2];
                    }
                }
            }

            for (int m = 0; m < neqn; ++ m)
            {
                f[m][targetCell] = fg[m][iFace];
            }
        }
    }
}

void DownloadInterfaceValue(UnstructGrid *grid, RDouble *f, const string &name)
{
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo)
    {
        return;
    }

    int nIFace = interfaceInfo->GetNIFace();

    for (int iGhost = GetNumberOfGhostCellLayers() - 1; iGhost >= 0; -- iGhost)
    {
        Data_ParamFieldSuite *recieveData = interfaceInfo->GetRecvDataStorage(iGhost);
        
        RDouble **fg = reinterpret_cast <RDouble **> (recieveData->GetDataPtr(name));
        
        for (int iFace = 0; iFace < nIFace; ++ iFace)
        {
            int targetCell;
            grid->GetTargetIndex(iFace, iGhost + 1, targetCell);
            f[targetCell] = fg[0][iFace];
        }
    }
}
void UploadInterpointValue(UnstructGrid *grid, RDouble **f, const string &name, int nEquation)
{
    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (!interpointInformation)
    {
        return;
    }

    int numberOfInterpoints = interpointInformation->GetNumberOfInterpoints();
    for (int iGhost = GetNumberOfGhostCellLayers() - 1; iGhost >= 0; -- iGhost)
    {
        Data_ParamFieldSuite *sendData = interpointInformation->GetSendDataStorage(iGhost);

        RDouble **fg = reinterpret_cast<RDouble **>(sendData->GetDataPtr(name));
        
        for (int iPoint = 0; iPoint < numberOfInterpoints; ++ iPoint)
        {
            int sourcePoint;
            grid->GetSourcePointIndex(iPoint, iGhost + 1, sourcePoint);
            for (int m = 0; m < nEquation; ++ m)
            {
                fg[m][iPoint] = f[m][sourcePoint];
            }
        }
    }
}

void DownloadInterpointValue(UnstructGrid *grid, RDouble **f, const string &name, int nEquation)
{
    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (!interpointInformation)
    {
        return;
    }
    int numberOfInterpoints = interpointInformation->GetNumberOfInterpoints();
    int *interPoint2GlobalPoint = interpointInformation->GetInterPoint2GlobalPoint();
    //for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
    {
        Data_ParamFieldSuite *recieveData = interpointInformation->GetReceiveDataStorage(0);
        RDouble **fg = reinterpret_cast<RDouble **>(recieveData->GetDataPtr(name));

        int GlobalPoint;
        for (int iPoint = 0; iPoint < numberOfInterpoints; ++ iPoint)
        {
            GlobalPoint = interPoint2GlobalPoint[iPoint];
            for (int m = 0; m < nEquation; ++ m)
            {
                f[m][GlobalPoint] += fg[m][iPoint];
            }
        }
    }
}

void UploadOversetValueDefault(UnstructGrid *grid, RDouble **field, const string &name, int neqn)
{
    OversetInformationProxy * oversetInformationProxy = grid->GetOversetInformationProxy();

    if (field == 0) return;
    
    OversetDataProxy *oversetDataProxy = oversetInformationProxy->GetOversetDataProxyForSend();

    Data_ParamFieldSuite *dataStorage = oversetDataProxy->GetDataStorage();

    RDouble **fieldStorage = reinterpret_cast< RDouble ** > (dataStorage->GetDataPtr(name));

    int numberOfNeighbors = oversetDataProxy->GetNumberOfNeighbors();
   /* ostringstream oss;
    oss << "variable :" << name << endl;*/
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++ iNeighbor)
    {
        int  numberOfCells           = oversetDataProxy->GetNumberOfCells          (iNeighbor);
        int *cellIndexContainer      = oversetDataProxy->GetCellIndexContainer     (iNeighbor);
        int *cellStorageIndexMapping = oversetDataProxy->GetCellStorageIndexMapping(iNeighbor);

        for (int iCell = 0; iCell < numberOfCells; ++ iCell)
        {
            int storageAddress = cellStorageIndexMapping[ iCell ];
            int cellIndex      = cellIndexContainer[ iCell ];
            //oss << "donorCellIndex " << cellIndex << endl;
            for (int iEquation = 0; iEquation < neqn; ++ iEquation)
            {
                fieldStorage[ iEquation ][ storageAddress ] = field[ iEquation ][ cellIndex ];
                //oss << field[ iEquation ][ cellIndex ] << endl;
            }
            
        }
    }
    
    /*int iBlock = grid->GetIBlock();
    if (iBlock == 0)
    {
       WriteLogFile(oss,true);
    }*/
}

void UploadOversetValueByWeight(UnstructGrid *grid, RDouble **field, const string &name, int numberOfEquations)
{
    OversetInformationProxy *oversetInformationProxy = grid->GetOversetInformationProxy();

    if (field == 0) return;

    OversetDataProxy *oversetDataProxy = oversetInformationProxy->GetOversetDataProxyForSend();
    Data_ParamFieldSuite *dataStorage = oversetDataProxy->GetDataStorage();

    RDouble **fieldStorage = reinterpret_cast<RDouble **> (dataStorage->GetDataPtr(name));

    int numberOfNeighbors = oversetDataProxy->GetNumberOfNeighbors();

    int numberOfNodes = grid->GetNTotalNode();
    RDouble **qNodeField = NewPointer2<RDouble>(numberOfEquations, numberOfNodes);

    ComputeNodeVariable(grid, qNodeField, field, numberOfEquations);

    int **cell2Node = grid->GetCell2NodeArray();
    int *nodeNumberOfEachCell = grid->GetNodeNumberOfEachCell();

    RDouble *xNode = grid->GetX();
    RDouble *yNode = grid->GetY();
    RDouble *zNode = grid->GetZ();

    for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++ iNeighbor)
    {
        int  numberOfCells           = oversetDataProxy->GetNumberOfCells          (iNeighbor);
        int *cellIndexContainer      = oversetDataProxy->GetCellIndexContainer     (iNeighbor);
        int *cellStorageIndexMapping = oversetDataProxy->GetCellStorageIndexMapping(iNeighbor);

        RDouble *interCellCenterX    = oversetDataProxy->GetInterCellCenterX(iNeighbor);
        RDouble *interCellCenterY    = oversetDataProxy->GetInterCellCenterY(iNeighbor);
        RDouble *interCellCenterZ    = oversetDataProxy->GetInterCellCenterZ(iNeighbor);

        vector< RDouble > nodeValue(numberOfEquations);

        for (int iCell = 0; iCell < numberOfCells; ++ iCell)
        {
            int storageAddress = cellStorageIndexMapping[ iCell ];
            int cellIndex      = cellIndexContainer[ iCell ];

            RDouble x = interCellCenterX[iCell];
            RDouble y = interCellCenterY[iCell];
            RDouble z = interCellCenterZ[iCell];

            SetField(nodeValue, zero, numberOfEquations);
            RDouble weight = 0.0;

            for(int iNode = 0; iNode < nodeNumberOfEachCell[ cellIndex ]; ++ iNode)
            {
                int &nodeIndex = cell2Node[cellIndex][iNode];

                RDouble dx = x - xNode[nodeIndex];
                RDouble dy = y - yNode[nodeIndex];
                RDouble dz = z - zNode[nodeIndex];

                RDouble dist = DISTANCE(dx, dy, dz) + SMALL;

                for (int iEquation = 0; iEquation < numberOfEquations; ++ iEquation)
                {
                    nodeValue[iEquation] += (1 / dist) * qNodeField[iEquation][nodeIndex];
                }
                weight += 1 / dist;
            }
            
            for (int iEquation = 0; iEquation < numberOfEquations; ++ iEquation)
            {
                nodeValue[iEquation] /= (weight + SMALL);
                fieldStorage[iEquation][storageAddress] = nodeValue[iEquation];
            }
        }
    }
    DelPointer2(qNodeField);
}

void UploadOversetValue(UnstructGrid *grid, RDouble **field, const string &name, int neqn)
{
    int  oversetInterpolationMethod = 0;
    if(GlobalDataBase::IsExist("oversetInterpolationMethod", PHINT, 1))
    {
        oversetInterpolationMethod = GlobalDataBase::GetIntParaFromDB("oversetInterpolationMethod");
    }

    if (oversetInterpolationMethod == 0)
    {
        UploadOversetValueDefault(grid, field, name, neqn);
    }
    else if (oversetInterpolationMethod == 1)
    {
        UploadOversetValueByWeight(grid, field, name, neqn);
    }
}
void UploadOversetValue(UnstructGrid *grid, RDouble *field, const string &name)
{
    OversetInformationProxy * oversetInformationProxy = grid->GetOversetInformationProxy();

    if (field == 0) return;

    OversetDataProxy *oversetDataProxy = oversetInformationProxy->GetOversetDataProxyForSend();

    Data_ParamFieldSuite *dataStorage = oversetDataProxy->GetDataStorage();

    RDouble **fieldStorage = reinterpret_cast< RDouble ** > (dataStorage->GetDataPtr(name));

    int numberOfNeighbors = oversetDataProxy->GetNumberOfNeighbors();
   /* ostringstream oss;
    oss << "variable :" << name << endl;*/
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++ iNeighbor)
    {
        int  numberOfCells           = oversetDataProxy->GetNumberOfCells          (iNeighbor);
        int *cellIndexContainer      = oversetDataProxy->GetCellIndexContainer     (iNeighbor);
        int *cellStorageIndexMapping = oversetDataProxy->GetCellStorageIndexMapping(iNeighbor);

        for (int iCell = 0; iCell < numberOfCells; ++ iCell)
        {
            int storageAddress = cellStorageIndexMapping[ iCell ];
            int cellIndex      = cellIndexContainer[ iCell ];
            //oss << "donorCellIndex " << cellIndex << endl;
            for (int iEquation = 0; iEquation < 1; ++iEquation)
            {
                fieldStorage[iEquation][storageAddress] = field[cellIndex];
                //oss << field[ cellIndex ] << endl;
            }

        }
    }
    
   /* int iBlock = grid->GetIBlock();
    if (iBlock == 0)
    {
        WriteLogFile(oss,true);
    }*/
}

void DownloadOversetValue(UnstructGrid *grid, RDouble **field, const string &name, int neqn)
{
    OversetInformationProxy *oversetInformationProxy = grid->GetOversetInformationProxy();

    if (field == 0) return;

    //PHVector1D < OversetCell * > *oversetCellContainer = oversetInformationProxy->GetOversetCellContainer();
    OversetDataProxy *oversetDataProxy = oversetInformationProxy->GetOversetDataProxyForReceive();

    Data_ParamFieldSuite *dataStorage = oversetDataProxy->GetDataStorage();

    RDouble **fieldStorage = reinterpret_cast< RDouble ** > (dataStorage->GetDataPtr(name));

    int numberOfNeighborZones = oversetDataProxy->GetNumberOfNeighbors();
    /*ostringstream oss;
    oss << "variable :" << name << endl;*/
    for (int iNeighbor = 0; iNeighbor < numberOfNeighborZones; ++ iNeighbor)
    {
        int  numberOfCells           = oversetDataProxy->GetNumberOfCells          (iNeighbor);
        int *cellIndexContainer      = oversetDataProxy->GetCellIndexContainer     (iNeighbor);
        int *cellStorageIndexMapping = oversetDataProxy->GetCellStorageIndexMapping(iNeighbor);

        for (int iCell = 0; iCell < numberOfCells; ++ iCell)
        {
            int storageAddress = cellStorageIndexMapping[ iCell ];
            int cellIndex      = cellIndexContainer     [ iCell ];
            /*OversetCell *oversetCell = (* oversetCellContainer)[iCell];
            int donorCellIndex = oversetCell->GetDonorCellIndex();
            oss << "donorCellIndex " << donorCellIndex << endl;*/
            //oss << "interCellIndex " << cellIndex << endl;
            for (int iEquation = 0; iEquation < neqn; ++ iEquation)
            {
                //oss << field[ iEquation ][ cellIndex ] << " ";
                field[ iEquation ][ cellIndex ] = fieldStorage[ iEquation ][ storageAddress ];
                //oss <<  fieldStorage[ iEquation ][ storageAddress ] << endl;
            }
        }
    }

    /*int iBlock = grid->GetIBlock();
    if (iBlock == 1)
    {
        WriteLogFile(oss,true);
    }*/
}

void DownloadOversetValue(UnstructGrid* grid, RDouble* field, const string& name)
{
    OversetInformationProxy * oversetInformationProxy = grid->GetOversetInformationProxy();

    if (field == 0) return;

    ///PHVector1D < OversetCell * > *oversetCellContainer = oversetInformationProxy->GetOversetCellContainer();

    OversetDataProxy * oversetDataProxy = oversetInformationProxy->GetOversetDataProxyForReceive();

    Data_ParamFieldSuite * dataStorage = oversetDataProxy->GetDataStorage();

    RDouble ** fieldStorage = reinterpret_cast< RDouble ** > (dataStorage->GetDataPtr(name));

    int numberOfNeighborZones = oversetDataProxy->GetNumberOfNeighbors();
   /* ostringstream oss;
    oss << "variable :" << name << endl;*/
    for (int iNeighbor = 0; iNeighbor < numberOfNeighborZones; ++ iNeighbor)
    {
        int  numberOfCells           = oversetDataProxy->GetNumberOfCells          (iNeighbor);
        int *cellIndexContainer      = oversetDataProxy->GetCellIndexContainer     (iNeighbor);
        int *cellStorageIndexMapping = oversetDataProxy->GetCellStorageIndexMapping(iNeighbor);

        for (int iCell = 0; iCell < numberOfCells; ++ iCell)
        {
            int storageAddress = cellStorageIndexMapping[ iCell ];
            int cellIndex      = cellIndexContainer     [ iCell ];
            /*OversetCell *oversetCell = (* oversetCellContainer)[iCell];
            int donorCellIndex = oversetCell->GetDonorCellIndex();
            oss << "donorCellIndex " << donorCellIndex << endl;*/
            for (int iEquation = 0; iEquation < 1; ++ iEquation)
            {
                //oss << field[ cellIndex ] << " ";
                field[ cellIndex ] = fieldStorage[ iEquation ][ storageAddress ];
                //oss << field[ cellIndex ] << endl;
            } 
        }
    }

    /*int iBlock = grid->GetIBlock();
    if (iBlock == 1)
    {
        WriteLogFile(oss,true);
    }*/
}
void BarthLimiterINF(Grid *grid_in, int ig, RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble * xfc = grid->GetFaceCenterX();
    RDouble * yfc = grid->GetFaceCenterY();
    RDouble * zfc = grid->GetFaceCenterZ();
    RDouble * xcc = grid->GetCellCenterX();
    RDouble * ycc = grid->GetCellCenterY();
    RDouble * zcc = grid->GetCellCenterZ();
    RDouble * xfn = grid->GetFaceNormalX();
    RDouble * yfn = grid->GetFaceNormalY();
    RDouble * zfn = grid->GetFaceNormalZ();

    int **cell2face  = grid->GetCell2Face();
    int *face_number_of_each_cell = grid->GetFaceNumberOfEachCell();

    RDouble dx, dy, dz, ds, dot;
    RDouble dq_face;

    int nIFace = interfaceInfo->GetNIFace();

    //! Need temporary arrays for the differences between the value in every cell and 
    //! maximum/minimum in the neighboring cells.
    RDouble *dmin = new RDouble[ nIFace ];
    RDouble *dmax = new RDouble[ nIFace ];

    //! Find the the differences for q.
    MinMaxDiffINF(grid, ig, dmin, dmax, q);

    //! Find the maximum and minimum in the neighbor of each cell.
    for (int iFace = 0; iFace < nIFace; ++ iFace)
    {
        int ie,face,le,re,iflag;
        grid->GetSourceIndex(iFace,ig+1,ie);
        for (int jFace = 0; jFace < face_number_of_each_cell[ie]; ++ jFace)
        {
            face = cell2face[ie][jFace];
            le   = leftCellOfFace[face];
            re   = rightCellOfFace[face];

            iflag = 1;

            if (re == ie) iflag = - 1;

            dx      = xfc[face] - xcc[ie];
            dy      = yfc[face] - ycc[ie];
            dz      = zfc[face] - zcc[ie];
            dq_face = dqdx[iFace] * dx + dqdy[iFace] * dy + dqdz[iFace] * dz;

            ds  = sqrt(dx * dx + dy * dy + dz * dz);
            dot = iflag * (xfn[face] * dx + yfn[face] * dy + zfn[face] * dz) / (ds + SMALL);

            if (dot < 0.0) limit[iFace] = 0.0;

            if (dq_face > dmax[iFace])
            {
                limit[iFace] = MIN(limit[iFace], dmax[iFace]/dq_face);
            }
            else if (dq_face < dmin[iFace])
            {
                limit[iFace] = MIN(limit[iFace], dmin[iFace]/dq_face);
            }
        }
    }
    delete [] dmin;    dmin = nullptr;
    delete [] dmax;    dmax = nullptr;
}

//! Set values for vencat limiters in 3D.
void VencatLimiterINF(Grid *grid_in, int ig, RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    int * leftCellOfFace = grid->GetLeftCellOfFace();
    int * rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();

    int **cell2face  = grid->GetCell2Face();
    int *face_number_of_each_cell = grid->GetFaceNumberOfEachCell();

    RDouble   dx, dy, dz, ds, dot;
    RDouble   dq_face, eps, tmp;

    int nIFace = interfaceInfo->GetNIFace();

    //! Need temporary arrays for the differences between the value in every cell and
    //! maximum/minimum in the neighboring cells.
    RDouble *dmin = new RDouble[ nIFace ];
    RDouble *dmax = new RDouble[ nIFace ];

    // Find the the differences for q
    MinMaxDiffINF(grid, ig, dmin, dmax, q);

    //! Find the maximum and minimum of dq in the field to set epsilons.
    RDouble dq_min = dmin[0], dq_max = dmax[0];
    for (int iFace = 0; iFace < nIFace; ++ iFace)
    {
        dq_min = MIN(dq_min, dmin[iFace]);
        dq_max = MAX(dq_max, dmax[iFace]);
    }

    RDouble venkatCoeff = 1.0e-5;
    GlobalDataBase::GetData("venkatCoeff", &venkatCoeff, PHDOUBLE, 1);

    eps = venkatCoeff * (dq_max - dq_min);
    eps = eps * eps + SMALL;

    //! Find the maximum and minimum in the neighbor of each cell.
    for (int iFace = 0; iFace < nIFace; ++ iFace)
    {
        int ie,face,le,re,iflag;
        grid->GetSourceIndex(iFace,ig+1,ie);
        for (int jFace = 0; jFace < face_number_of_each_cell[ie]; ++ jFace)
        {
            face = cell2face[ie][jFace];
            le   = leftCellOfFace[face];
            re   = rightCellOfFace[face];

            dx = xcc[re] - xcc[le];
            dy = ycc[re] - ycc[le];
            dz = zcc[re] - zcc[le];
            
            ds  = sqrt(dx * dx + dy * dy + dz * dz);
            eps = venkatCoeff * ds;
            eps = eps * eps * eps;

            iflag = 1;

            if (re == ie) iflag = - 1;

            dx      = xfc[face] - xcc[ie];
            dy      = yfc[face] - ycc[ie];
            dz      = zfc[face] - zcc[ie];
            dq_face = dqdx[iFace] * dx + dqdy[iFace] * dy + dqdz[iFace] * dz;

            ds  = sqrt(dx * dx + dy * dy + dz * dz);
            dot = iflag * (xfn[face] * dx + yfn[face] * dy + zfn[face] * dz) / (ds + SMALL);

            //if (dot < 0.0) limit[iFace] = 0.0;

            if (dq_face > SMALL)
            {
                tmp = VenFun(dmax[iFace], dq_face, eps);
                limit[iFace] = MIN(limit[iFace], tmp);
            }
            else if (dq_face < - SMALL)
            {
                tmp = VenFun(dmin[iFace], dq_face, eps);
                limit[iFace] = MIN(limit[iFace], tmp);
            }
        }
    }

    delete [] dmin;    dmin = nullptr;
    delete [] dmax;    dmax = nullptr;
}

void SmoothField(UnstructGrid *grid, RDouble *q)
{
    int nTotalCell = grid->GetNTotalCell();
    int **cell2face  = grid->GetCell2Face();
    int *face_number_of_each_cell = grid->GetFaceNumberOfEachCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *vol   = grid->GetCellVolume();

    RDouble *tmp = new RDouble[nTotalCell];

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        q  [iCell] *= vol[iCell];
        tmp[iCell]  = q[iCell];
    }

    RDouble epsl = 0.2;
    int nsmooth = 0;
    do
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            RDouble sum = 0.0;
            int count = 0;
            for (int iFace = 0; iFace < face_number_of_each_cell[iCell]; ++ iFace)
            {
                int face = cell2face[iCell][iFace];
                int le   = leftCellOfFace [ face ];
                int re   = rightCellOfFace[ face ];
                int ie   = le;
                if (le == iCell) ie = re;

                if (ie >= nTotalCell) continue;

                sum += tmp[ie];
                ++ count;
            }
            RDouble coef = one / (one + epsl * count);
            tmp[iCell] = (q[iCell] + epsl * sum) * coef;
        }
        ++ nsmooth;
    } while (nsmooth < 3);

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble ovol = 1.0 / vol[iCell];
        q[iCell] = tmp[iCell] * ovol;
    }
    delete [] tmp;    tmp = nullptr;
}

void addn2c(int ie, int ip, int *nCPN, struct linkbase **nc)
{
    int find = 0;

    struct linkbase *link;
    struct linkbase *newie;

    for(link = nc[ip]; link != NULL; link = link->next)
    {
        if(link->ic == ie)
        {
            find = 1;
            break;
        }
    }

    if(find == 0)
    {
        nCPN[ ip ] ++;
        newie = new struct linkbase;
        newie->ic = ie;
        newie->next = nc[ ip ];
        nc[ ip ] = newie;
    }
}

LIB_EXPORT void ConstructFace2Node(int * nPointPerFace, int *nodesOfPerFace, uint_t nFaces, vector< vector<int> > &face2node)
{
    face2node.resize(nFaces);

    int count = 0;
    for (uint_t iFace = 0; iFace < nFaces; ++ iFace)
    {
        int nNode = nPointPerFace[iFace];
        face2node[iFace].resize(nNode);
        for (int iNode = 0; iNode < nNode; ++ iNode)
        {
            int nodeIndex = nodesOfPerFace[count++];
            face2node[iFace][iNode] = nodeIndex;
        }
    }
}

void GetFaceToNodeList(UnstructGrid *grid, int boundaryConditionType, PHVectorInt1D &linkMap, PHVectorInt2D &faceToNodeList)
{
    int numberOfBoundaryFaces = grid->GetNBoundFace();

    int *faceNodeNumberContainer = grid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer  = grid->GetFace2Node();

    //! Build point set of binary tree for searching.
    PHVectorInt1D locFaceToNode(4);
    set< DataStruct_Sort< int > > setForSorting;
    DataStruct_Sort< int > dataForSorting;

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    faceToNodeList.resize(0);
    linkMap.resize(0);

    int nodePosition = 0;
    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);

        int faceNodeNumber = faceNodeNumberContainer[ iFace ];
        if (boundaryConditionType != bcRegion->GetBCType())
        {
            nodePosition += faceNodeNumber;
            continue;
        }
        locFaceToNode.resize(0);
        for (int iNodeInFace = 0; iNodeInFace < faceNodeNumber; ++ iNodeInFace)
        {
            int index = faceNodeIndexContainer[ nodePosition + iNodeInFace ];

            dataForSorting.value = index;
            set< DataStruct_Sort< int > >::iterator iter = setForSorting.find(dataForSorting);

            if (iter == setForSorting.end())
            {
                linkMap.push_back(index);
                dataForSorting.index = static_cast<int>(linkMap.size() - 1);

                setForSorting.insert(dataForSorting);
                locFaceToNode.push_back(dataForSorting.index);
            }
            else
            {
                locFaceToNode.push_back(iter->index);
            }

        }
        nodePosition += faceNodeNumber;
        faceToNodeList.push_back(locFaceToNode);
    }
}


void VisualizationMesh2D(UnstructGrid *grid, const string &fileName)
{
    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();
    int numberOfFaces = grid->GetNTotalFace();

    int *faceNodeNumberContainer = grid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer  = grid->GetFace2Node();

    int *leftCellIndexContainer  = grid->GetLeftCellOfFace();
    int *rightCellIndexContainer = grid->GetRightCellOfFace();

    PHVectorString1D titleOfTecplot;
    titleOfTecplot.push_back("title=\"Unstructured Grid\"");
    titleOfTecplot.push_back("variables=");
    titleOfTecplot.push_back("\"x\"");
    titleOfTecplot.push_back("\"y\"");
    titleOfTecplot.push_back("\"z\"");

    fstream file;
    PHSPACE::OpenFile(file, fileName, ios_base::out);

    for (size_t i = 0; i < titleOfTecplot.size(); ++ i)
    {
        file << titleOfTecplot[ i ] << "\n";
    }

    int iCount = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        iCount += faceNodeNumberContainer[ iFace ];
    }

    int numberOfWordsInEachLine = 5;

    //! Output for Tecplot.
    file << "ZONE\n";
    if (GetDim() == THREE_D)
    {
        file << "ZoneType = FEPolyhedron\n";
    }
    else
    {
        file << "ZoneType = FEPolygon\n";
    }
    file << "Nodes    = " << numberOfNodes << "\n";
    file << "Faces    = " << numberOfFaces << "\n";
    file << "Elements = " << numberOfCells << "\n";
    file << "TotalNumFaceNodes = " << iCount << "\n";
    file << "NumConnectedBoundaryFaces = 0\n";
    file << "TotalNumBoundaryConnections = 0\n";

    file.setf(ios::scientific);
    file.precision(6);

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        file << x[ iNode ] << " ";
        if ((iNode + 1) % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (numberOfNodes % numberOfWordsInEachLine != 0) file << "\n";

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        file << y[ iNode ] << " ";
        if ((iNode + 1) % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (numberOfNodes % numberOfWordsInEachLine != 0) file << "\n";

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        file << z[ iNode ] << " ";
        if ((iNode + 1) % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (numberOfNodes % numberOfWordsInEachLine != 0) file << "\n";

    if (GetDim() == THREE_D)
    {
        //! Geometry
        //! First the faceNodeNumberContainer(n poin pf face).
        for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
        {
            file << faceNodeNumberContainer[ iFace ] << " ";
            if ((iFace + 1) % numberOfWordsInEachLine == 0) file << "\n";
        }
        if (numberOfFaces % numberOfWordsInEachLine == 0) file << "\n";
    }

    //! Second the points to face.
    iCount = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int faceNodeNumber = faceNodeNumberContainer[ iFace ];
        for (int iNodeInFace = 0; iNodeInFace < faceNodeNumber; ++ iNodeInFace)
        {
            file << faceNodeIndexContainer[ iCount ++ ] + 1 << " ";
        }
        file << "\n";
    }

    //! Third the f2c:left.
    int left;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        left = leftCellIndexContainer[ iFace ] + 1;
        file << left << " ";
        if ((iFace + 1) % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (numberOfFaces % numberOfWordsInEachLine == 0) file << "\n";

    //! Third the f2c:right.
    int right;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        right = rightCellIndexContainer[ iFace ] + 1;
        if (right > numberOfCells || right < 0) right = 0;
        file << right << " ";
        if ((iFace + 1) % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (numberOfFaces % numberOfWordsInEachLine == 0) file << "\n";

    PHSPACE::CloseFile(file);
}

void VisualizationMesh3D(UnstructGrid *gridIn, const string &fileName)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    int numberOfBoundaryFaces = grid->GetNBoundFace();
    int numberOfFaces         = grid->GetNTotalFace();

    int *faceNodeNumberContainer = grid->GetNodeNumberOfEachFace();

    PHVectorString1D titleOfTecplot;
    titleOfTecplot.push_back("title=\"Unstructured Grid\"");
    titleOfTecplot.push_back("variables=");
    titleOfTecplot.push_back("\"x\"");
    titleOfTecplot.push_back("\"y\"");
    titleOfTecplot.push_back("\"z\"");

    int iCount = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        iCount += faceNodeNumberContainer[ iFace ];
    }

    using namespace PHENGLEI;
    PHVectorInt2D faceToNodeList(numberOfBoundaryFaces);
    PHVectorInt1D linkMap;
    set< int > bcList;
    bool doseNeedDump = false;

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        bcList.insert(bcType);

        if (bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::SYMMETRY
            || bcType == PHENGLEI::INFLOW || bcType == PHENGLEI::OUTFLOW)
        {
            doseNeedDump = true;
        }
    }

    if (PHMPI::GetNumberofGlobalZones() >= 256)
    {
        doseNeedDump = false;

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
            int bcType = bcRegion->GetBCType();
          
            if (bcType == PHENGLEI::SOLID_SURFACE)
            {
                doseNeedDump = true;
                break;
            }
        }
    }

    if (!doseNeedDump)
    {
        return;
    }

    fstream file;
    PHSPACE::OpenFile(file, fileName, ios_base::out);

    for (size_t i = 0; i < titleOfTecplot.size(); ++ i)
    {
        file << titleOfTecplot[ i ] << "\n";
    }

    for (set< int >::iterator iter = bcList.begin(); iter != bcList.end(); ++ iter)
    {
        if (*iter < 0) continue;
        if (*iter != PHENGLEI::SOLID_SURFACE && *iter != PHENGLEI::SYMMETRY 
            && *iter != PHENGLEI::OUTFLOW && *iter != PHENGLEI::INFLOW) continue;

        faceToNodeList.resize(0);
        linkMap.resize(0);

        GetFaceToNodeList(grid, *iter, linkMap, faceToNodeList);
        ostringstream osZone;
        osZone << "\"" << grid->GetGridID()->GetIndex()<< "BC=" << *iter << "\" ";
        string zoneTitle = osZone.str();

        file << "zone T = " << zoneTitle << " N = " << linkMap.size() << "  E = " << faceToNodeList.size()  << " f = FEPOINT, ET = QUADRILATERAL\n";
        int wordWidth = 20;
        for (std::size_t i = 0; i < linkMap.size(); ++ i)
        {
            int it = linkMap[ i ];
            file << setiosflags(ios::left);
            file << setiosflags(ios::scientific);
            file << setprecision(10);
            file << setw(wordWidth) << x[ it ]
            << setw(wordWidth) << y[ it ]
            << setw(wordWidth) << z[ it ];
            file << "\n";
        }

        for (std::size_t i = 0; i < faceToNodeList.size(); ++ i)
        {
            //! Number of vertices on the cell face.
            uint_t np     = faceToNodeList[ i ].size();
            int index0 = faceToNodeList[ i ][ 0 ];
            for (uint_t ip = 0; ip < np; ++ ip)
            {
                file << faceToNodeList[ i ][ ip ] + 1;
                if (ip != np - 1) file << " ";
            }

            for (uint_t ip = np; ip < 4; ++ ip)
            {
                file << " " << index0 + 1;
            }
            file << "\n";
        }
    }

    PHSPACE::CloseFile(file);
}


void ExtractSubGridFromGrid(UnstructGrid *grid, UnstructGrid *grid_Slice, int *isThisCellNeedExtract)
{
    int count;
    int nTotalCell_Slice, nTNode_Slice, nTFace_Slice, nBFace_Slice;
    int nTotalNode = grid->GetNTotalNode();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();

    int *nNPF = grid->GetNodeNumberOfEachFace();
    int *f2n  = grid->GetFace2Node();

    int *f2cL = grid->GetLeftCellOfFace();
    int *f2cR = grid->GetRightCellOfFace();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    int *lable_cell = new int [nTotalCell];
    for(int iCell = 0; iCell < nTotalCell; ++ iCell) 
    {
        lable_cell[iCell] = 0;
    }

    //! Find out the total number of nodes,faces and cells on the slice.
    int *c_to_sc = new int [nTotalCell];
    count = 0;
    for(int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if(isThisCellNeedExtract[iCell])
        {
            lable_cell[iCell] = 1;
            c_to_sc[iCell]    = count;
            count ++;
        }
    }
    nTotalCell_Slice = count;
    grid_Slice->SetNTotalCell(nTotalCell_Slice);

    int *lable_face = new int [nTotalFace];
    for(int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        lable_face[iFace] = 0;
    }
    count = 0;
    for(int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int c1 = f2cL[iFace];
        if(lable_cell[c1])
        {
            ++ count;
            lable_face[iFace] = 1;
        }
    }
    for(int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int c1 = f2cL[iFace];
        int c2 = f2cR[iFace];
        if(lable_cell[c1] || lable_cell[c2])
        {
            count ++;
            lable_face[iFace] = 1;
        }
    }
    nTFace_Slice = count;
    grid_Slice->SetNTotalFace(nTFace_Slice);

    count = 0;
    for(int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int c1 = f2cL[iFace];
        if(lable_cell[c1])
            ++ count;
    }
    for(int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int c1 = f2cL[iFace];
        int c2 = f2cR[iFace];
        if((lable_cell[c1] && !lable_cell[c2]) || (!lable_cell[c1] && lable_cell[c2]))
            count ++;
    }
    nBFace_Slice = count;
    grid_Slice->SetNBoundFace(nBFace_Slice);

    int *lable_point = new int [nTotalNode];
    for(int iFace = 0; iFace < nTotalNode; ++ iFace)
    {
        lable_point[iFace] = 0;
    }
    set<int> Slice_node;
    count = 0;
    for(int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        if(!lable_face[iFace]) 
        {
            for(int j = 0; j < nNPF[iFace]; ++ j)
                count ++;
            continue;
        }
        for(int j = 0; j < nNPF[iFace]; ++ j)
        {
            int ip = f2n[count ++];
            lable_point[ip] = 1;
            Slice_node.insert(ip);
        }
    }
    nTNode_Slice = static_cast<int>(Slice_node.size());
    grid_Slice->SetNTotalNode(nTNode_Slice);

    //! Find out the coordinates of sub-grid nodes.
    RDouble *x_Slice = new RDouble[nTNode_Slice];
    RDouble *y_Slice = new RDouble[nTNode_Slice];
    RDouble *z_Slice = new RDouble[nTNode_Slice];
    int     *p_to_sp = new int [nTotalNode];

    count = 0;
    for(int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if(lable_point[iNode])
        {
            x_Slice[count] = x[iNode];
            y_Slice[count] = y[iNode];
            z_Slice[count] = z[iNode];

            p_to_sp[iNode] = count;

            count ++;
        }
    }
    grid_Slice->SetX(x_Slice);
    grid_Slice->SetY(y_Slice);
    grid_Slice->SetZ(z_Slice);

    //! Faces.
    int *nNPF_Slice = new int [nTFace_Slice];
    int *sf_to_f    = new int [nTFace_Slice];
    int *f2cL_Slice = new int [nTFace_Slice];
    int *f2cR_Slice = new int [nTFace_Slice];
    count = 0;
    for(int iFace = 0; iFace < nTotalFace; ++ iFace) 
    {
        if(lable_face[iFace])
            count += nNPF[iFace];
    }
    int *f2n_Slice  = new int [count];

    //! The f2c and nNPF of the sub-grid faces.
    count = 0;
    int countL = 0;
    int countR = 0;
    for(int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        if(lable_face[iFace])
        {
            int c1 = f2cL[iFace];
            int c2 = f2cR[iFace];

            if(c2 >= nTotalCell || c2 < 0)
            {
                if(!lable_cell[c1]) 
                    f2cL_Slice[countL++] = -1;
                else
                    f2cL_Slice[countL++] = c_to_sc[c1];

                f2cR_Slice[countR++] = -1;
            }
            else
            {
                if(!lable_cell[c1]) 
                    c_to_sc[c1] = -1;
                if(!lable_cell[c2]) 
                    c_to_sc[c2] = -1;

                f2cL_Slice[countL++] = c_to_sc[c1];
                f2cR_Slice[countR++] = c_to_sc[c2];
            }

            nNPF_Slice[count] = nNPF[iFace];

            sf_to_f[count] = iFace;

            count ++;
        }
    }

    //! The f2n and nNPF of the sub-grid faces.
    count = 0;
    int count2 = 0;
    for(int iFace =0; iFace < nTotalFace; ++ iFace)
    {
        if(!lable_face[iFace])
        {
            for(int j = 0; j < nNPF[iFace]; ++ j) count ++;
        }
        else
        {
            for(int j = 0; j < nNPF[iFace]; ++j)
            {
                int ip = f2n[count ++];
                f2n_Slice[count2++] = p_to_sp[ip];
            }
        }
    }
    grid_Slice->SetNodeNumberOfEachFace(nNPF_Slice);
    grid_Slice->SetFace2Node(f2n_Slice);

    grid_Slice->SetLeftCellOfFace(f2cL_Slice);
    grid_Slice->SetRightCellOfFace(f2cR_Slice);

    delete [] lable_cell;    lable_cell = nullptr;
    delete [] lable_face;    lable_face = nullptr;
    delete [] lable_point;    lable_point = nullptr;
    delete [] sf_to_f;    sf_to_f = nullptr;
    delete [] c_to_sc;    c_to_sc = nullptr;
    delete [] p_to_sp;    p_to_sp = nullptr;
}

void DumpSliceGrid(UnstructGrid *unstructGrid, int sliceAxis, RDouble slicePosition, string fileName)
{
    //PrintToWindow("Start dump slice grid ...\n");
    int nTotalFace = unstructGrid->GetNTotalFace();
    int nTotalNode = unstructGrid->GetNTotalNode();
    int nTotalCell = unstructGrid->GetNTotalCell();

    int *lable_point = new int [nTotalNode];
    int *lable_face  = new int [nTotalFace];
    int *lable_cell  = new int [nTotalCell];

    PHSPACE::SetField(lable_point, 0, nTotalNode);
    PHSPACE::SetField(lable_face,  0, nTotalFace);
    PHSPACE::SetField(lable_cell,  0, nTotalCell);

    int *face2node = unstructGrid->GetFace2Node();
    int *leftCellOfFace  = unstructGrid->GetLeftCellOfFace();
    int *rightCellOfFace = unstructGrid->GetRightCellOfFace();
    RDouble *coord = 0;
    if (sliceAxis == 0) 
        coord = unstructGrid->GetX();
    else if (sliceAxis == 1) 
        coord = unstructGrid->GetY();
    else 
        coord = unstructGrid->GetZ();

    long long int *face2NodeSubscript = unstructGrid->GetFace2NodeSubscript();
    int ip;

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        RDouble rmin =  1.0e30;
        RDouble rmax = -1.0e30;

        for (int j = static_cast<int>(face2NodeSubscript[iFace]); j < face2NodeSubscript[iFace+1]; ++ j)
        {
            ip = face2node[j];
            if (coord[ip] > rmax)
            {
                rmax = coord[ip];
            }

            if (coord[ip] < rmin)
            {
                rmin = coord[ip];
            }
        }

        if ((rmax - slicePosition) * (rmin - slicePosition) <= 0.0)
        {
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];

            lable_cell[le] = 1;
            if (re >= 0 && re < nTotalCell) 
            {
                lable_cell[re] = 1;
            }
        }
    }
    int dumpCellCount = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if(lable_cell[iCell] == 1)
        {
            ++ dumpCellCount;
            break;
        }
    }

    if(dumpCellCount != 0)
    {
        int iZone = unstructGrid->GetGridID()->GetIndex();
        Grid *grid_slice = CreateGridGeneral(UNSTRUCTGRID, new GridID(iZone), 0, PHSPACE::GetDim());
        ExtractSubGridFromGrid(unstructGrid, UnstructGridCast(grid_slice), lable_cell);

        VisualizationMesh2D(UnstructGridCast(grid_slice), fileName);
        delete grid_slice;
    }

    delete [] lable_point;    lable_point = nullptr;
    delete [] lable_face;    lable_face = nullptr;
    delete [] lable_cell;    lable_cell = nullptr;

    //PrintToWindow("End dump slice grid ...\n");
}

void VisualizationEigenMeshofCoarseGrid3D(UnstructGrid *grid, const string &fileName)
{
    if(GetDim() < 3) return;

    UnstructGrid *coarseGrid = grid;
    grid = (UnstructGrid *)coarseGrid->GetFineGrid();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    int nTotalNode = grid->GetNTotalNode();
    int nBoundFace = grid->GetNBoundFace();

    int *nNPF = grid->GetNodeNumberOfEachFace();
    int *f2n  = grid->GetFace2Node();

    int *leftCellIndexContainer  = grid->GetLeftCellOfFace();
    int *rightCellIndexContainer = grid->GetRightCellOfFace();
    
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    int *mark = nullptr;
    int node = 0, n = 0, i = 0, j = 0, bcType = 0, nline = 0;

    int nTotalCell = 0;

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        bcType = bcRegion->GetBCType();
        if (bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::SYMMETRY)
        {
            vector<int>* faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                nTotalCell ++;
            }
        }
    }

    if (nTotalCell == 0) return;

    mark = new int[nTotalNode]();

    //! Write fine grid first.
    fstream file;
    PHSPACE::OpenFile(file, fileName, ios_base::out);

    file << "TITLE     = \"Grids of MultiGrid in 3D\"" << endl;
    file << "VARIABLES = \"X\", \"Y\", \"Z\"" << endl;

    file << "ZONE T=\"Level " << grid->GetLevel() << " \"" << endl;
    file << " N= " << nTotalNode << ", E= " << nTotalCell << ", F=FEPOINT, ET=Quadrilateral" << endl;
    for(i=0; i<nTotalNode; i++)
    {
        file << x[i] << " " << y[i] << " " << z[i] << endl;
    }
    for(node=i=0; i<nBoundFace; i++) 
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[i]);
        bcType = bcRegion->GetBCType();

        if(bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::SYMMETRY)
        {
            for(j=0; j<nNPF[i]; j++) 
            {
                n = f2n[node++];
                file << n + 1 << " ";
            }
            if(nNPF[i] == 3) 
            {   //! triangle
                file << n + 1 << " " << endl;
            }
            else if(nNPF[i] == 4)
            {   //! quad
                file << endl;
            } else 
            {
                cout << "Can't out this element type." << endl;
            }
        }else
        {
            for(j=0; j<nNPF[i]; ++j) node ++;
        }
    }

    //! Compute the information of line to point on finest grid boundary.
    int *fpmark = new int [nTotalNode]();
    for (i = 0; i < nTotalNode; ++ i) fpmark[i] = -1;
    int fnBDNode = 0;
    for (node=i=0; i<nBoundFace; ++i)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[i]);
        bcType = bcRegion->GetBCType();

        if(!(bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::SYMMETRY))
        {
            node += nNPF[i];
            continue;
        }
        for(j=0; j<nNPF[i]; ++j)
        {
            int ip     = f2n[node++];
            fpmark[ip] = 1;
        }
    }
    for(node=i=0; i<nTotalNode; ++i)
    {
        if(fpmark[i] == 1)
        {
            fpmark[i] = node++;
        }
    }
    fnBDNode = node;
    int *fpmarkrever = new int [fnBDNode];
    for(i=0; i<nTotalNode; ++i)
    {
        if(fpmark[i] != -1)
        {
            int ip = fpmark[i];
            fpmarkrever[ip] = i;
        }
    }

    int (*l2p)[2] = new int [3*nBoundFace][2](); 
    int fnBDLine   = 0;
    int *p2l      = new int [fnBDNode]();
    for(i=0; i<fnBDNode; ++i) p2l[i] = -1;
    for(node=i=0; i<nBoundFace; ++i)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[i]);
        bcType = bcRegion->GetBCType();

        if(!(bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::SYMMETRY))
        {
            node += nNPF[i];
            continue;
        }
        int p1, p2;
        for(j=0; j<nNPF[i]; ++j)
        {
            if(j==0)
            {
                p1 = f2n[ node+nNPF[i]-1 ];
                p2 = f2n[node++];
            }else
            {
                p1 = f2n[node-1];
                p2 = f2n[node++];
            }

            p1 = fpmark[p1]; 
            p2 = fpmark[p2];

            if(p2l[p1] == -1 || p2l[p2] == -1)
            {
                l2p[fnBDLine][0] = p1;
                l2p[fnBDLine][1] = p2;

                if(p2l[p1] == -1) p2l[p1] = fnBDLine;
                if(p2l[p2] == -1) p2l[p2] = fnBDLine;

                fnBDLine ++;

                if(fnBDLine > 3*nBoundFace)
                {
                    cout << "Error: Please increase the length of p2l[][2]!" << endl;
                    return;
                }
            }
            else if (p2l[p1] != p2l[p2])
            {
                bool flag = true;
                for(int k=0; k<fnBDLine; ++k)
                {
                    int n1 = l2p[k][0], n2 = l2p[k][1];
                    if((n1 == p1 && n2 == p2) || (n1 == p2 && n2 == p1))
                    {
                        flag = false;
                        break;
                    }
                }
                if(flag)
                {
                    l2p[fnBDLine][0] = p1;
                    l2p[fnBDLine][1] = p2;
                    fnBDLine ++;

                    if(fnBDLine > 3*nBoundFace)
                    {
                        cout << "Error: Please increase the length of p2l[][2]!" << endl;
                        return;
                    }
                }
            }
        }
    }
    delete [] p2l;    p2l = nullptr;

    //! Compute the information of point to cell on finest grid boundary.
    vector< set<int> > p2c(fnBDNode);
    for(node=i=0; i<nBoundFace; ++i)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[i]);
        bcType = bcRegion->GetBCType();

        if(!(bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::SYMMETRY))
        {
            node += nNPF[i];
            continue;
        }
        for(j=0; j<nNPF[i]; ++j)
        {
            int ip = f2n[node++];
            ip     = fpmark[ip];
            p2c[ip].insert(i);
        }
    }

    //! Compute the information of line to cell on finest grid boundary.
    int (*l2c)[2] = new int [fnBDLine][2]();
    int *nCPL     = new int [fnBDLine]();
    int l;
    for(l=0; l<fnBDLine; ++l) nCPL[l] = 0;

    set<int>::iterator iter1, iter2;
    for(l=0; l<fnBDLine; ++l)
    {
        int c1, c2;
        int p1 = l2p[l][0], p2 = l2p[l][1];
        for(iter1=p2c[p1].begin(); iter1!= p2c[p1].end(); ++iter1)
        {
            c1 = *iter1;
            for(iter2=p2c[p2].begin(); iter2!=p2c[p2].end(); ++iter2)
            {
                c2 = *iter2;

                if(c1 == c2)
                {
                    if(nCPL[l] == 2)
                    {
                        cout << "Warning: this line has had two cells!\n";
                        return;
                    }
                    l2c[l][ nCPL[l]++ ] = c1;
                }
            }
        }
    }

    //! For coarse grid.
    int *lmark = new int [fnBDLine]();
    grid     = coarseGrid;
    
    leftCellIndexContainer  = grid->GetLeftCellOfFace();
    rightCellIndexContainer = grid->GetRightCellOfFace();

    //! Compute the border of cells on boundary.
    int *bdf2cc = new int [nBoundFace]();
    for(i=0; i<nBoundFace; ++i) bdf2cc[i] = -1;

    for(i=0; i<fnBDLine; ++i) lmark[i] = 0;
    for(i=0; i<fnBDLine; ++i)
    {
        if(nCPL[i] == 2)
        {
            int c1 = l2c[i][0], c2 = l2c[i][1];
            if(leftCellIndexContainer[c1] != leftCellIndexContainer[c2])
            {
                lmark[i] = 1;
            }
            
        }
        else if (nCPL[i] == 1)
        {
            lmark[i] = 1;
        }
        else
        {
            cout << "Error: this line has wrong cells!\n";
            return;
        }

    }
    delete [] bdf2cc;    bdf2cc = nullptr;

    //! Compute the line(l2n) to count.
    nline = 0;
    for(i=0; i<fnBDLine; ++i)
    {
        if (lmark[i])
        {
            ++ nline;
        }
    }

    file << "ZONE T=\"Level " << grid->GetLevel() << "\"" << endl;
    file <<" N= " << fnBDNode << ", E= "  << nline << ", F=FEPOINT, ET=Quadrilateral\n";
    for(i=0; i<nTotalNode; i++)
    {
        if(fpmark[i] != -1) 
        {
            file << x[i] << " " << y[i] << " " << z[i] << endl;
        }
    }

    for(i=0; i<fnBDLine; ++i)
    {
        if(lmark[i])
        {
            int n1 = l2p[i][0], n2 = l2p[i][1];
            file << n1 + 1 << " " << n1 + 1 << " " << n2 + 1 << " " << n2 + 1 << endl;
        }
    }

    delete [] lmark;    lmark = nullptr;
    delete [] fpmark;    fpmark = nullptr;
    delete [] fpmarkrever;    fpmarkrever = nullptr;

    delete [] l2c;    l2c = nullptr;
    delete [] nCPL;    nCPL = nullptr;
    delete [] l2p;    l2p = nullptr;

    delete [] mark;    mark = nullptr;

    CloseFile(file);
}

vector< RDouble > & UnstructGrid::GetNodeWallDist()
{
    return nodeWalldistDistance;
}

vector< RDouble > & UnstructGrid::GetNodeWallDistInOtherBlock()
{
    return nodeWalldistDistanceInOtherBlock;
}

LIB_EXPORT void UnstructGrid::RegisterOversetField(const string &name, const int type, const int dimesion)
{
    //! Field for send
    OversetDonorCellManager *oversetDonorCellManager = oversetInformationProxy->GetOversetDonorCellManager();
    int numberOfDonorCells = oversetDonorCellManager->GetTotalNumberOfOversetDonorCells();
    Data_ParamFieldSuite *dataStorageForSend = oversetDonorCellManager->GetDataStorage();

    RDouble **fieldForSend = NewPointer2< RDouble >(dimesion, numberOfDonorCells);

    for(int iDim = 0; iDim < dimesion; ++ iDim)
    {
        for (int jDim = 0; jDim < numberOfDonorCells; ++ jDim)
        {
            fieldForSend[ iDim ][ jDim ] = 0.0;
        }
    }

    dataStorageForSend->UpdateDataPtr(name, fieldForSend);

    //! Field for receive.
    OversetInterplateCellManager *oversetInterplateCellManager = oversetInformationProxy->GetOversetInterplateCellManager();
    int numberOfInterplateCells = oversetInterplateCellManager->GetTotalNumberOfOversetInterplateCells();
    Data_ParamFieldSuite *dataStorageForReceive = oversetInterplateCellManager->GetDataStorage();

    RDouble **fieldForReceive = NewPointer2< RDouble >(dimesion, numberOfInterplateCells);

    for(int iDim = 0; iDim < dimesion; ++ iDim)
    {
        for (int jDim = 0; jDim < numberOfInterplateCells; ++ jDim)
        {
            fieldForReceive[ iDim ][ jDim ] = 0.0;
        }
    }

    dataStorageForReceive->UpdateDataPtr(name, fieldForReceive);
}

LIB_EXPORT RDouble UnstructGrid::GetAverageVolume()
{
    if(averageVolume < 0.0)
    {
        cout<<endl<<"Error! averageVolume < 0.0"<<endl;
        exit(0);
    }
    return averageVolume;
}

void UnstructGrid::ConstructBCRegion()
{
    if (unstructBCSet)
    {
        delete unstructBCSet;
    }

    set<string> bcNameSet;
    for (int iFace = 0; iFace < nBoundFace; ++iFace)
    {
        bcNameSet.insert(bcr[iFace]->GetBCName());
    }

    int nBCRegionUnstruct = static_cast<int>(bcNameSet.size());
    CreateUnstructBCSet(nBCRegionUnstruct);

    //UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int* bcRegionIDofBCFace = new int[nBoundFace];
    set<string>::iterator iter;
    uint_long count = 0;
    for (iter = bcNameSet.begin(); iter != bcNameSet.end(); iter++)
    {
        UnstructBC* unstructBC = new UnstructBC(count);
        unstructBCSet->SetBCRegion(count, unstructBC);
        unstructBC->SetBCName(*iter);

        for (int iFace = 0; iFace < nBoundFace; iFace++)
        {
            if (bcr[iFace]->GetBCName() == unstructBC->GetBCName())
            {
                unstructBC->SetBCType(bcr[iFace]->GetKey());
                bcRegionIDofBCFace[iFace] = count;
                vector<int>* faceIndex = unstructBC->GetFaceIndex();
                faceIndex->push_back(iFace);
            }

        }
        count++;
    }
    unstructBCSet->SetBCFaceInfo(bcRegionIDofBCFace);
}

void UnstructGrid::Output(const string & fileName)
{
    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();

    int numberOfBoundaryFaces = this->GetNBoundFace();
    int numberOfFaces         = this->GetNTotalFace();

    int *faceNodeNumberContainer = this->GetNodeNumberOfEachFace();

    PHVectorString1D titleOfTecplot;
    titleOfTecplot.push_back("title=\"Unstructured Grid\"");
    titleOfTecplot.push_back("variables=");
    titleOfTecplot.push_back("\"x\"");
    titleOfTecplot.push_back("\"y\"");
    titleOfTecplot.push_back("\"z\"");

    int iCount = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        iCount += faceNodeNumberContainer[ iFace ];
    }

    using namespace PHENGLEI;
    PHVectorInt2D faceToNodeList(numberOfBoundaryFaces);
    PHVectorInt1D linkMap;
    set< int > bcList;

    UnstructBCSet **bcr = this->GetBCRecord();

    bool doseNeedDump = false;
    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
    {
        int boundaryConditionType = bcr[iFace]->GetKey();
        bcList.insert(boundaryConditionType);
        if (boundaryConditionType == PHENGLEI::SOLID_SURFACE || boundaryConditionType == PHENGLEI::SYMMETRY
            || boundaryConditionType == PHENGLEI::INFLOW || boundaryConditionType == PHENGLEI::OUTFLOW)
        {
            doseNeedDump = true;
        }
    }

    if (PHMPI::GetNumberofGlobalZones() >= 256)
    {
        doseNeedDump = false;
        for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
        {
            int boundaryConditionType = bcr[iFace]->GetKey();
            if(boundaryConditionType == PHENGLEI::SOLID_SURFACE)
            {
                doseNeedDump = true;
                break;
            }
        }
    }

    if (!doseNeedDump)
    {
        return;
    }

    fstream file;
    PHSPACE::OpenFile(file, fileName, ios_base::out);

    for (size_t i = 0; i < titleOfTecplot.size(); ++ i)
    {
        file << titleOfTecplot[ i ] << "\n";
    }

    for (set< int >::iterator iter = bcList.begin(); iter != bcList.end(); ++ iter)
    {
        if (*iter < 0 ) continue;
        if (*iter != PHENGLEI::SOLID_SURFACE && *iter != PHENGLEI::SYMMETRY 
            && *iter != PHENGLEI::OUTFLOW && *iter != PHENGLEI::INFLOW) continue;

        faceToNodeList.resize(0);
        linkMap.resize(0);

        GetFaceToNodeList(this, * iter, linkMap, faceToNodeList);
        ostringstream osZone;
        osZone << "\"" << this->GetGridID()->GetIndex()<< "BC=" << *iter << "\" ";
        string zoneTitle = osZone.str();

        file << "zone T = " << zoneTitle << " N = " << linkMap.size() << "  E = " << faceToNodeList.size()  << " f = FEPOINT, ET = QUADRILATERAL\n";
        int wordWidth = 20;
        for (std::size_t i = 0; i < linkMap.size(); ++ i)
        {
            int it = linkMap[ i ];
            file << setiosflags(ios::left);
            file << setiosflags(ios::scientific);
            file << setprecision(10);
            file << setw(wordWidth) << x[ it ]
                << setw(wordWidth) << y[ it ]
                << setw(wordWidth) << z[ it ];
            file << "\n";
        }

        for (std::size_t i = 0; i < faceToNodeList.size(); ++ i)
        {
            //! Number of vertices on the cell face.
            uint_t np     = faceToNodeList[ i ].size();
            int index0 = faceToNodeList[ i ][ 0 ];
            for (uint_t ip = 0; ip < np; ++ ip)
            {
                file << faceToNodeList[ i ][ ip ] + 1;
                if (ip != np - 1) file << " ";
            }

            for (uint_t ip = np; ip < 4; ++ ip)
            {
                file << " " << index0 + 1;
            }
            file << "\n";
        }
    }

    PHSPACE::CloseFile(file);
}

RDouble AspectRatio(RDouble dist, RDouble area, int nodeNumOfPerFace)
{
    const RDouble SQRT3 = 1.732050807568877;

    RDouble length = 1.0;

    if (GetDim() == TWO_D)
    {
        if (nodeNumOfPerFace == 2)
        {   //! quadrilateral
            length = area;
        }
        else
        {
            cout << endl << "Need new code in function aspect ratio!" << endl;
        }
    }
    else
    {
        if (nodeNumOfPerFace == 3)
        {   //! triangle
            length = sqrt(area * four / SQRT3);
        }
        else if (nodeNumOfPerFace == 4)
        {   //! quadrilateral
            length = sqrt(area);
        }
        else
        {
            cout << endl << "Need new code in function aspect ratio!" << endl;
        }
    }

    return length / (dist + TINY);
}

void Get_Xadj_Adjncy(UnstructGrid *grid, idx_t *xadj, idx_t *adjncy)
{
    int  nTotalCell = grid->GetNTotalCell();
    vector<int> *c2c = grid->GetCell2Cell();

    xadj[0]   = 0;
    int count = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        xadj[iCell+1] = xadj[iCell] + c2c[iCell].size();
        for (std::size_t j = 0; j < c2c[iCell].size(); ++ j)
        {
            adjncy[count ++] = c2c[iCell][j];
        }
    }
}

void Get_Xadj_Adjncy(UnstructGrid *grid, idx_t *xadj, idx_t *adjncy, int *neighborCell, idx_t *vtxdist)
{
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();

    int *left_cell_of_face = grid->GetLeftCellOfFace();
    int *right_cell_of_face = grid->GetRightCellOfFace();

    idx_t startGlobalIndex = 0;

    //! First: compute the cell-to-cell connectivity relationship, include interface.
    set<idx_t> *c2c = new set<idx_t>[nTotalCell];
    idx_t le, re;

    //! If boundary is an INTERFACE
    if (vtxdist && neighborCell)
    {
        //! Parallel
        int nIFace = grid->GetNIFace();
        InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
        int *interFace2ZoneID = interfaceInfo->GetInterFace2ZoneID();
        int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

        for (int iFace = 0; iFace < nIFace; ++ iFace)
        {
            int neighborZone = interFace2ZoneID[iFace];
            int bcFace = interFace2BoundaryFace[iFace];
            le = left_cell_of_face[bcFace];
            re = neighborCell[iFace];
            re += vtxdist[neighborZone];
            c2c[le].insert(re);
        }

        int myid = PHMPI::GetCurrentProcessorID();
        startGlobalIndex = vtxdist[myid];
    }

    //! Interior faces
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = left_cell_of_face [iFace];
        re = right_cell_of_face[iFace];

        c2c[le].insert(re + startGlobalIndex);
        c2c[re].insert(le + startGlobalIndex);
    }

    //! Second: compute the xadj and adjncy.
    xadj[0]   = 0;
    int count = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        xadj[iCell+1] = xadj[iCell] + c2c[iCell].size();
        for (set<idx_t>::iterator iter = c2c[iCell].begin(); iter != c2c[iCell].end(); ++ iter)
        {
            adjncy[count ++] = (*iter);
        }
    }

    //! Free memory.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        c2c[iCell].clear();
    }
    delete [] c2c;    c2c = nullptr;
}

void Get_Xadj_Adjncy_Weight(UnstructGrid *grid, idx_t *xadj, idx_t *adjncy, RDouble *adjncyWeight)
{
    int nTotalCell  = grid->GetNTotalCell();
    int **cell2face = grid->GetCell2Face();
    int *nFaceCell  = grid->GetFaceNumberOfEachCell();

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *area = grid->GetFaceArea();

    int count = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int nNeighbor = static_cast<int>(xadj[iCell + 1] - xadj[iCell]);

        for (int j = 0; j < nNeighbor; ++ j)
        {
            int neighborCell = static_cast<int>(adjncy[count]);

            RDouble weightArea = -1.0;
            int nFace = nFaceCell[iCell];
            for (int iNeighFace = 0; iNeighFace < nFace; ++ iNeighFace)
            {
                int neighFace = cell2face[iCell][iNeighFace];

                int le = leftCellOfFace[neighFace];
                int re = rightCellOfFace[neighFace];

                if( (le == iCell && re == neighborCell) ||
                    (le == neighborCell && re == iCell) )
                {
                    weightArea = area[neighFace];
                    break;
                }
            }   //! neighbor face.

            if(weightArea < 0.0)
            {
                TK_Exit::PrintDebugInfoExit("Error: can not find the patched face!");
            }

            adjncyWeight[count] = weightArea;
            count ++;
        }   //! neighbor cell
    }
}

}


