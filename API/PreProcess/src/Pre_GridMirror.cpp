#include "Pre_GridMirror.h"
#include "Math_BasisFunction.h"
#include "Glb_Dimension.h"
#include "IO_FileName.h"
#include "Zone.h"
#include "Pre_GridConversion.h"
#include "PHIO.h"
#include "TK_Log.h"

using namespace std;

namespace PHSPACE
{
LIB_EXPORT Pre_GridMirror::Pre_GridMirror(const string &from_gfile)
{
    this->gridFileName = from_gfile;

    numberOfZones   = 0;
    gridtype        = 0;
    BlockofEachProc = 0;
    GlobalDataBase::GetData("ndim", &Dimension, PHINT, 1);
    GlobalDataBase::GetData("SymmetryFaceVector", &SymmetryFaceVector, PHINT, 1);

    gridinfo = new DataContainer();

    OrdinaryGrid.resize(0);
    //MirrorGrids.resize(0);
    Gridlist.resize(0);
}

LIB_EXPORT Pre_GridMirror::~Pre_GridMirror()
{
    delete gridinfo;
}

LIB_EXPORT void Pre_GridMirror::Run()
{
    using namespace PHMPI;

    string out_gfile = "grid.fts";
    GlobalDataBase::GetData("out_gfile", &out_gfile, PHSTRING, 1);

    cout << "Task type	: " << "MirrorGrid" << endl;
    cout << "  from_gfile	: " << gridFileName << endl;
    cout << "  out_gfile	: " << out_gfile << endl << endl;

    GlobalDataBase::UpdateData("gridfile", &gridFileName, PHSTRING, 1);
    GlobalDataBase::GetData("gridtype", &gridtype, PHINT, 1);

    Region *region = new Region();
    region->ReadGrid();

    numberOfOriginalBlocks = PHMPI::GetNumberofGlobalZones();
    for (int iZone = 0; iZone < numberOfOriginalBlocks; ++ iZone)
    {
        Zone *zone = region->GetZone(iZone);
        if (!zone) continue;

        Grid *grid = zone->GetGeometry()->GetGrid();
        OrdinaryGrid.push_back(grid);
        BlockofEachProc ++;
    }

    //numberOfOriginalBlocks = PHMPI::GetNumberofGlobalZones();
    /*OrdinaryGrid = new Grid *[numberOfOriginalBlocks];

    for (int iBlock = 0; iBlock < numberOfOriginalBlocks; iBlock ++)
    {
        OrdinaryGrid[iBlock] = GetGrid(iBlock, 0);
    }*/

    cout << endl << "End reading grid ... " << endl << endl;
    RevertBCInfo();

    //CalSymmetryFace();

    CalBlocksNumberOfMirrorGrid();
    GridMirror();
    MirrorBCRegion();

    cout << endl << "Mirror grid finished ... " << endl << endl;

    TradeGridInfo();

    if (Gridlist.size() > 0)
    {
        if (gridtype == STRUCTGRID)
        {
            DumpStrGrid(out_gfile);
        }
        else
        {
            DumpUnstrGrid(out_gfile);
        }
    }

    /*for (int iBlock = 0; iBlock < BlockofEachProc; iBlock ++)
    {
        delete OrdinaryGrid[iBlock];
    }*/

    /*for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        delete grids[iZone];
        delete MirrorGrids[iZone];
    }
    delete [] grids;
    delete [] MirrorGrids;*/

    delete region;

}

void Pre_GridMirror::TradeGridInfo()
{
    using namespace PHMPI;

    int myid = GetCurrentProcessorID();

    int nProc = GetNumberOfProcessor();

    for (int iProc = 0; iProc < nProc; ++ iProc)
    {
        DataContainer *cdata = 0;
        if (myid == 0 && iProc != 0)
        {
            cdata = new DataContainer();
        }
        else
        {
            cdata = gridinfo;
        }

        PH_Trade(cdata, iProc, 0, 0);

        if (myid == 0)
        {
            Gridlist.push_back(cdata);
            gridinfo = new DataContainer();
        }
    }
}

void Pre_GridMirror::DumpUnstrGrid(string out_gfile)
{
    /*if (numberOfZones > 1)
    {
        cout << "\n";
        cout << "Blocks > 1, Start to Combine Grids" << endl;
        CombinGrid * combinGrid = new CombinGrid(OutputGrid, numberOfZones);

        combinGrid->RunCombinGrid();
        
        Grid ** gridMerged = combinGrid->GetGridAll();
        
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            delete OutputGrid[iZone];
        }

        OutputGrid = 0;
        OutputGrid = gridMerged;

        DumpGrid(out_gfile, gridMerged, 1);

        delete combinGrid;
    }
    else
    {
        DumpGrid(out_gfile, OutputGrid, numberOfZones);
    }

    delete OutputGrid[0];*/

    Grid **grids = new Grid * [numberOfZones];

    uint_t sizes = Gridlist.size();

    for (int isize = 0; isize < sizes; ++ isize)
    {
        DataContainer *cdata = Gridlist[isize];
        cdata->MoveToBegin();

        int nblock = 0;
        PHRead(cdata, &nblock, 1);

        vector < int > indexlist;
        for (int iZone = 0; iZone < nblock; ++ iZone)
        {
            int Gridindex = 0;
            PHRead(cdata, &Gridindex, 1);
            indexlist.push_back(Gridindex);

            GridID *index = new GridID(Gridindex);
            Grid *grid = CreateGridGeneral(UNSTRUCTGRID, index, 0, Dimension);
            grids[Gridindex] = grid;

            GridID *mirrorIndex = new GridID(Gridindex + numberOfOriginalBlocks);
            Grid *Mirror = CreateGridGeneral(UNSTRUCTGRID, mirrorIndex, 0, Dimension);
            grids[Gridindex + numberOfOriginalBlocks] = Mirror;

            int nTotalNodes, nTotalFace, nTotalCells;
            PHRead(cdata, &nTotalNodes, 1);
            PHRead(cdata, &nTotalFace, 1);
            PHRead(cdata, &nTotalCells, 1);

            grid->SetNTotalCell(nTotalCells);
            grid->SetNTotalNode(nTotalNodes);
            grid->SetNTotalFace(nTotalFace);

            Mirror->SetNTotalCell(nTotalCells);
            Mirror->SetNTotalNode(nTotalNodes);
            Mirror->SetNTotalFace(nTotalFace);


            RDouble *x = new RDouble[nTotalNodes];
            RDouble *y = new RDouble[nTotalNodes];
            RDouble *z = new RDouble[nTotalNodes];

            RDouble *x_mir = new RDouble[nTotalNodes];
            RDouble *y_mir = new RDouble[nTotalNodes];
            RDouble *z_mir = new RDouble[nTotalNodes];

            PHRead(cdata, x, nTotalNodes);
            PHRead(cdata, y, nTotalNodes);
            PHRead(cdata, z, nTotalNodes);

            PHRead(cdata, x_mir, nTotalNodes);
            PHRead(cdata, y_mir, nTotalNodes);
            PHRead(cdata, z_mir, nTotalNodes);

            grid->SetX(x);
            grid->SetY(y);
            grid->SetZ(z);

            Mirror->SetX(x_mir);
            Mirror->SetY(y_mir);
            Mirror->SetZ(z_mir);

            grid->ComputeMinMaxBox();
            Mirror->ComputeMinMaxBox();

            int *node_number_of_each_face = new int[nTotalFace];
            int *node_number_of_each_face_mir = new int[nTotalFace];
            UnstructGridCast(grid)->SetNodeNumberOfEachFace(node_number_of_each_face);
            UnstructGridCast(Mirror)->SetNodeNumberOfEachFace(node_number_of_each_face_mir);

            PHRead(cdata, node_number_of_each_face, nTotalFace);
            PHRead(cdata, node_number_of_each_face_mir, nTotalFace);

            int nsum;
            PHRead(cdata, &nsum, 1);

            int *face2node = new int [nsum];
            int *face2node_mir = new int [nsum];
            UnstructGridCast(grid)->SetFace2Node(face2node);
            UnstructGridCast(Mirror)->SetFace2Node(face2node_mir);

            PHRead(cdata, face2node, nsum);
            PHRead(cdata, face2node_mir, nsum);

            int *left_cell_of_face      = new int[nTotalFace];
            int *left_cell_of_face_mir  = new int[nTotalFace];
            int *right_cell_of_face     = new int[nTotalFace];
            int *right_cell_of_face_mir = new int[nTotalFace];

            UnstructGridCast(grid)->SetLeftCellOfFace(left_cell_of_face);
            UnstructGridCast(grid)->SetRightCellOfFace(right_cell_of_face);
            UnstructGridCast(Mirror)->SetLeftCellOfFace(left_cell_of_face_mir);
            UnstructGridCast(Mirror)->SetRightCellOfFace(right_cell_of_face_mir);

            PHRead(cdata, left_cell_of_face, nTotalFace);
            PHRead(cdata, left_cell_of_face_mir, nTotalFace);
            PHRead(cdata, right_cell_of_face, nTotalFace);
            PHRead(cdata, right_cell_of_face_mir, nTotalFace);

            int cell2node = 0;
            PHRead(cdata, &cell2node, 1);
            if (cell2node == 1)
            {
                int *node_number_of_each_cell = new int [nTotalCells];
                int *node_number_of_each_cell_mir = new int [nTotalCells];
                UnstructGridCast(grid)->SetNodeNumberOfEachCell(node_number_of_each_cell);
                UnstructGridCast(Mirror)->SetNodeNumberOfEachCell(node_number_of_each_cell_mir);

                PHRead(cdata, node_number_of_each_cell, nTotalCells);
                PHRead(cdata, node_number_of_each_cell_mir, nTotalCells);

                int num;
                PHRead(cdata, &num, 1);

                int *cell2nodeA = new int [num];
                int *cell2node_mir = new int [num];
                UnstructGridCast(grid)->SetCell2Node(cell2nodeA);
                UnstructGridCast(Mirror)->SetCell2Node(cell2node_mir);

                PHRead(cdata, cell2nodeA, num);
                PHRead(cdata, cell2node_mir, num);
            }
        }

        for (int iZone = 0; iZone < nblock; ++ iZone)
        {
            int Gridindex = indexlist[iZone];

            UnstructGrid *grid   = UnstructGridCast(grids[Gridindex]);
            UnstructGrid *Mirror = UnstructGridCast(grids[Gridindex + numberOfOriginalBlocks]);

            int nBoundFace;
            PHRead(cdata, &nBoundFace, 1);

            grid->SetNBoundFace(nBoundFace);
            Mirror->SetNBoundFace(nBoundFace);

            UnstructBCSet **bcrs     = new UnstructBCSet * [nBoundFace];
            UnstructBCSet **bcrs_mir = new UnstructBCSet * [nBoundFace];

            grid->SetBCRecord(bcrs);
            Mirror->SetBCRecord(bcrs_mir);

            int count;
            PHRead(cdata, &count, 1);

            InterfaceInfo *interfaceInfo = 0;
            InterfaceInfo *iinfo_mir = 0;

            if (count != 0)
            {
                interfaceInfo     = new InterfaceInfo(count);
                iinfo_mir = new InterfaceInfo(count);
            }

            grid->SetInterfaceInfo(interfaceInfo);
            Mirror->SetInterfaceInfo(iinfo_mir);

            int *interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
            int *interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();
            int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace(); 

            int *interFace2ZoneID_mir       = iinfo_mir->GetInterFace2ZoneID();
            int *interFace2InterFaceID_mir  = iinfo_mir->GetInterFace2InterFaceID();
            int *interFace2BoundaryFace_mir = iinfo_mir->GetInterFace2BoundaryFace(); 

            count = 0;
            for (int iFace = 0; iFace < nBoundFace; ++ iFace)
            {
                int bctype;
                string bcName;
                PHRead(cdata, &bctype, 1);
                cdata->ReadString(bcName);

                bcrs[iFace] = new UnstructBCSet();
                bcrs[iFace]->SetKey(bctype);
                bcrs[iFace]->SetBCName(bcName);

                bcrs_mir[iFace] = new UnstructBCSet();
                bcrs_mir[iFace]->SetKey(bctype);
                bcrs_mir[iFace]->SetBCName(bcName);

                if (bctype == PHENGLEI::INTERFACE)
                {
                    int interFace2BcFace, interFace2Zone, interFace2InterFace;

                    PHRead(cdata, &interFace2BcFace, 1);
                    PHRead(cdata, &interFace2Zone, 1);
                    PHRead(cdata, &interFace2InterFace, 1);

                    interFace2BoundaryFace[count] = interFace2BcFace;
                    interFace2ZoneID[count] = interFace2Zone;
                    interFace2InterFaceID[count] = interFace2InterFace;

                    PHRead(cdata, &interFace2BcFace, 1);
                    PHRead(cdata, &interFace2Zone, 1);
                    PHRead(cdata, &interFace2InterFace, 1);

                    interFace2BoundaryFace_mir[count] = interFace2BcFace;
                    interFace2ZoneID_mir[count] = interFace2Zone;
                    interFace2InterFaceID_mir[count] = interFace2InterFace;

                    count++;
                }
            }
        }

        delete cdata;    cdata = nullptr;
    }

    /*fstream file;
    file.open("results/mirror.log", ios_base::out|ios_base::trunc);

    ostringstream oss;

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        UnstructGrid *grid = UnstructGridCast(grids[iZone]);

        oss << "\n" << "\n";

        oss << "Grid block number: " << iZone << "\n" << "\n";

        oss << "write nodes info: " << "\n";

        int nTotalNode = grid->GetNTotalNode();
        int nTotalFace = grid->GetNTotalFace();
        int nTotalCells = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();

        oss << "	Number of nodes: " << nTotalNode << "\n";
        
        RDouble *x = grid->GetX();
        RDouble *y = grid->GetY();

        oss << "	Output node coordinates: " << "\n";

        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            oss << "		node: " << iNode << " coordinate: " << x[iNode] << "	" << y[iNode] << "\n";
        }

        oss << "\n";

        oss << "write nodes info: " << "\n";
        oss << "	Number of faces: " << nTotalFace << "\n";

        int *face2node          = grid->GetFace2Node();
        int *left_cell_of_face  = grid->GetLeftCellOfFace();
        int *right_cell_of_face = grid->GetRightCellOfFace();
        int *cell2node          = grid->GetCell2Node();
       
        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

        InterfaceInfo *iinfo_ordinary = grid->GetInterfaceInfo();

        int *interFace2BoundaryFace = iinfo_ordinary->GetInterFace2BoundaryFace(); 
        int *interFace2ZoneID       = iinfo_ordinary->GetInterFace2ZoneID();
        int *interFace2InterFaceID  = iinfo_ordinary->GetInterFace2InterFaceID();

        int count = 0;
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            oss << "\n";

            oss << "	Face number: " << iFace << "\n";
            oss << "	face2node:	" << face2node[iFace * 2] << "	" << face2node[iFace * 2 + 1] << "\n";

            oss << "	left_cell:	" << left_cell_of_face[iFace] << "\n";
            oss << "	right_cell:	" << right_cell_of_face[iFace] << "\n";

            if (iFace < nBoundFace)
            {
                UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);

                int bcType = bcRegion->GetBCType();
                oss << "The Boundary condition is :	" << bcType << "\n";

                if (bcType == -1)
                {
                    oss << "	interFace2BoundaryFace:	" << interFace2BoundaryFace[count] << "\n";
                    oss << "	interFace2ZoneID:	" << interFace2ZoneID[count] << "\n";
                    oss << "	interFace2InterFaceID:	" << interFace2InterFaceID[count] << "\n";
                    count ++ ;
                }
            }
        }

        oss << "\n";

        oss << "write nodes info: " << "\n";
        oss << "	Number of cells: " << nTotalCells << "\n";

        for (int iCell = 0; iCell < nTotalCells; ++ iCell)
        {
            oss << "	Cell number: " << iCell << "\n";
            oss << "	cell2node:	" << cell2node[iCell * 4] << "	" << cell2node[iCell * 4 + 1] << "	" << cell2node[iCell * 4 + 2] << "	" << cell2node[iCell * 4 + 3] << "\n";
        }
    }

    file << oss.str();
    file.close();
    file.clear();*/

    //ConstructGlobalInterfaceLink(grids, numberOfZones);
    DumpGrid(out_gfile, grids, numberOfZones);
}

void Pre_GridMirror::DumpStrGrid(string out_gfile)
{
    Grid **grids = new Grid * [numberOfZones];

    uint_t sizes = Gridlist.size();

    for (int isize = 0; isize < sizes; ++ isize)
    {
        DataContainer *cdata = Gridlist[isize];
        cdata->MoveToBegin();

        int nblock = 0;
        PHRead(cdata, &nblock, 1);

        vector < int > indexlist;
        for (int iZone = 0; iZone < nblock; ++ iZone)
        {
            int Gridindex = 0;
            PHRead(cdata, &Gridindex, 1);
            indexlist.push_back(Gridindex);

            GridID *index = new GridID(Gridindex);
            Grid *grid;
            grid = CreateGridGeneral(STRUCTGRID, index, 0, Dimension);
            grids[Gridindex] = grid;

            GridID *index_mirror = new GridID(Gridindex + numberOfOriginalBlocks);
            Grid *Mirror;
            Mirror = CreateGridGeneral(STRUCTGRID, index_mirror, 0, Dimension);
            grids[Gridindex + numberOfOriginalBlocks] = Mirror;

            int nTotalNode;
            PHRead(cdata, &nTotalNode, 1);

            int ni = 0;
            int nj = 0;
            int nk = 0;
            PHRead(cdata, &ni, 1);
            PHRead(cdata, &nj, 1);
            PHRead(cdata, &nk, 1);
            
            StructGridCast(grid)->SetNI(ni);
            StructGridCast(grid)->SetNJ(nj);
            StructGridCast(grid)->SetNK(nk);
            StructGridCast(grid)->SetBasicDimension();

            StructGridCast(Mirror)->SetNI(ni);
            StructGridCast(Mirror)->SetNJ(nj);
            StructGridCast(Mirror)->SetNK(nk);
            StructGridCast(Mirror)->SetBasicDimension();

            RDouble *x = new RDouble [nTotalNode];
            RDouble *y = new RDouble [nTotalNode];
            RDouble *z = new RDouble [nTotalNode];

            grid->SetX(x);
            grid->SetY(y);
            grid->SetZ(z);
            StructGridCast(grid)->SetArrayLayout();

            RDouble *MirrorX = new RDouble [nTotalNode];
            RDouble *MirrorY = new RDouble [nTotalNode];
            RDouble *MirrorZ = new RDouble [nTotalNode];

            Mirror->SetX(MirrorX);
            Mirror->SetY(MirrorY);
            Mirror->SetZ(MirrorZ);
            StructGridCast(Mirror)->SetArrayLayout();

            PHRead(cdata, x, nTotalNode);
            PHRead(cdata, y, nTotalNode);
            PHRead(cdata, z, nTotalNode);

            PHRead(cdata, MirrorX, nTotalNode);
            PHRead(cdata, MirrorY, nTotalNode);
            PHRead(cdata, MirrorZ, nTotalNode);

            StructGridCast(grid)->ComputeMinMaxBox();
            StructGridCast(Mirror)->ComputeMinMaxBox();

            /*delete index;
            delete index_mirror;*/
        }

        for (int iZone = 0; iZone < nblock; ++ iZone)
        {

            int Gridindex = indexlist[iZone];

            Grid *grid = grids[Gridindex];
            Grid *Mirror = grids[Gridindex + numberOfOriginalBlocks];

            int  numberOfBoundaryFaces = 0;
            PHRead(cdata, &numberOfBoundaryFaces, 1);

            StructGridCast(grid)->CreateCompositeBCRegion(numberOfBoundaryFaces);
            StructGridCast(Mirror)->CreateCompositeBCRegion(numberOfBoundaryFaces);

            StructBCSet *rhs = StructGridCast(grid)->GetStructBCSet();
            StructBCSet *structBCSet = StructGridCast(Mirror)->GetStructBCSet();

            for (int iBCRegion = 0; iBCRegion < numberOfBoundaryFaces; iBCRegion ++)
            {
                StructBC *bcregion = new StructBC(Gridindex, iBCRegion);
                StructBC *MirrorBCRegion = new StructBC(Gridindex + numberOfOriginalBlocks, iBCRegion);

                rhs->SetBCRegion(iBCRegion, bcregion);
                structBCSet->SetBCRegion(iBCRegion, MirrorBCRegion);

                int bctype = 0;
                PHRead(cdata, &bctype, 1);

                bcregion->SetBCType(bctype);
                MirrorBCRegion->SetBCType(bctype);

                int *s_st = bcregion->GetStartPoint();
                int *s_ed = bcregion->GetEndPoint();

                PHRead(cdata, s_st, 3);
                PHRead(cdata, s_ed, 3);

                int *s_st_mir = MirrorBCRegion->GetStartPoint();
                int *s_ed_mir = MirrorBCRegion->GetEndPoint();
                PHRead(cdata, s_st_mir, 3);
                PHRead(cdata, s_ed_mir, 3);

                /*delete bcregion;
                delete MirrorBCRegion;*/
            }
        }
        delete cdata;    cdata = nullptr;
    }
    
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        StructGridCast(grids[iZone])->ProcessBCInfo();
    }

    ConstructGlobalInterfaceLink(grids, numberOfZones);
    DumpGrid(out_gfile, grids, numberOfZones);

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        delete grids[iZone];
    }
    delete [] grids;
}

void Pre_GridMirror::ReadLnkFile(const string &fileName_in)
{
    string LnkFileName = ChangeExtensionOfFileName(fileName_in, "link");
    if (! SerialFileExist(LnkFileName))
    {
        TK_Exit::ExceptionExit("There is no link file!");
        return;
    }
    LnkFileExist = true;

    fstream file;
    OpenSerialFile(file, LnkFileName, ios_base::in|ios_base::binary);

    DataContainer *LnkData = new DataContainer();
    ReadFile(file, LnkData);

    VirtualFile *virtualFile = new VirtualFile(LnkData);
    virtualFile->BeginReadWork();

    //Grid ** OrdinaryGrid = GetOrdinaryGrid();
    for (int iBlock = 0; iBlock < BlockofEachProc; iBlock ++)
    {
        StructGrid *grid = StructGridCast(OrdinaryGrid[iBlock]);

        StructBCSet *structBCSet = grid->GetStructBCSet();

        int nBCRegion = structBCSet->GetnBCRegion();
        for (int iBCRegion = 0; iBCRegion < nBCRegion; iBCRegion ++)
        {
            StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);

            int bctype = structBC->GetBCType();
            if ((!IsInterface(bctype))) continue;

            int imin, imax, jmin, jmax, kmin, kmax, nbt;

            PHRead(virtualFile, imin);
            PHRead(virtualFile, imax);
            PHRead(virtualFile, jmin);
            PHRead(virtualFile, jmax);

            if (Dimension == THREE_D)
            {
                PHRead(virtualFile, kmin);
                PHRead(virtualFile, kmax);
            }
            else
            {
                kmin = 1;
                kmax = 1;
            }

            structBC->SetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);

            /*int * t_st = structBC->GetTargetStart();
            int * t_ed = structBC->GetTargetEnd();*/

            PHRead(virtualFile, imin);
            PHRead(virtualFile, imax);
            PHRead(virtualFile, jmin);
            PHRead(virtualFile, jmax);

            if (Dimension == THREE_D)
            {
                PHRead(virtualFile, kmin);
                PHRead(virtualFile, kmax);
            }
            else
            {
                kmin = 1;
                kmax = 1;
            }

            PHRead(virtualFile, nbt);

            structBC->SetTargetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);
            structBC->SetTargetRegionBlock(nbt - 1);
        }
    }

    virtualFile->EndReadWork();

    delete virtualFile;
    delete LnkData;

    file.close();
    file.clear();
}

void Pre_GridMirror::RevertBCInfo()
{
    if (gridtype == STRUCTGRID)
    {
        for (int iBlock = 0; iBlock < BlockofEachProc; iBlock ++)
        {
            StructGrid *grid = StructGridCast(OrdinaryGrid[iBlock]);
            StructBCSet *structBCSet = grid->GetStructBCSet();

            int numberOfBoundaryFaces = structBCSet->GetnBCRegion();
            for (int iBCRegion = 0; iBCRegion < numberOfBoundaryFaces; iBCRegion ++)
            {
                StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);

                int *s_st = structBC->GetStartPoint();
                int *s_ed = structBC->GetEndPoint();
                int *s_lr3d = structBC->GetFaceDirectionIndex();
                for (int m = 0; m < Dimension; ++ m)
                {
                    if (s_st[m] == s_ed[m])
                    {
                        if (s_st[m] > 1 || s_lr3d[m] == 1)
                        //if (s_st[m] != 1)
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
            }
        }

        //ReadLnkFile(gridFileName);
    }
}

void Pre_GridMirror::CalBlocksNumberOfMirrorGrid()
{
    /*int SymmetryBlock = 0;
    for (int iZone = 0; iZone < numberOfOriginalBlocks; iZone ++)
    {
        StructGrid * grid = StructGridCast(OrdinaryGrid[iZone]);
        StructBCSet * structBCSet = grid->GetStructBCSet();
        int numberOfBoundaryFaces = structBCSet->GetnBCRegion();

        for (int iBCRegion = 0; iBCRegion < numberOfBoundaryFaces; iBCRegion ++)
        {
            StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);

            int bctype = structBC->GetBCType();

            if (bctype != PHENGLEI::SYMMETRY) continue;

            SymmetryBlock += 1;
            break;
        }
    }

    numberOfZones = numberOfOriginalBlocks * 2 - SymmetryBlock;*/

    numberOfZones = numberOfOriginalBlocks * 2;
    //OutputGrid = new Grid * [numberOfZones];
}

void Pre_GridMirror::CalSymmetryFace()
{
    int Nodeindex1[3] = {0}, Nodeindex2[3] = {0}, Nodeindex3[3] = {0}, Nodeindex4[3] = {0};
    vector < vector <RDouble> > Node;
    RDouble node1[3] = {0}, node2[3] = {0}, node3[3] = {0};

    Node.resize(3);

    bool SymmetryFace = false;
    for (int iZone = 0; iZone < numberOfOriginalBlocks; iZone ++)
    {
        if (SymmetryFace == true) break;

        StructGrid *grid = StructGridCast(OrdinaryGrid[iZone]);
        StructBCSet *structBCSet = grid->GetStructBCSet();

        RDouble AA, BB, CC, DD;

        int numberOfBoundaryFaces = structBCSet->GetnBCRegion();
        for (int iBCRegion = 0; iBCRegion < numberOfBoundaryFaces; iBCRegion ++)
        {
            StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);

            int bctype = structBC->GetBCType();

            if (bctype != PHENGLEI::SYMMETRY) continue;

            int *s_st = structBC->GetStartPoint();
            int *s_ed = structBC->GetEndPoint();

            bool j = false;
            for (int i = 0; i < 3; i ++)
            {
                if (s_st[i] == s_ed[i])
                {
                    Nodeindex1[i] = s_st[i];
                    Nodeindex2[i] = s_st[i];
                    Nodeindex3[i] = s_st[i];
                    Nodeindex4[i] = s_st[i];
                }
                else
                {
                    if (j == false)
                    {
                        Nodeindex1[i] = s_st[i];
                        Nodeindex2[i] = s_st[i];

                        Nodeindex3[i] = s_ed[i];
                        Nodeindex4[i] = s_ed[i];

                        j = true;
                    }
                    else
                    {
                        Nodeindex1[i] = s_st[i];
                        Nodeindex3[i] = s_st[i];

                        Nodeindex2[i] = s_ed[i];
                        Nodeindex4[i] = s_ed[i];
                    }
                }
            }

            SymmetryFace = true;
            break;
        }

        RDouble3D &xx = * grid->GetStructX();
        RDouble3D &yy = * grid->GetStructY();
        RDouble3D &zz = * grid->GetStructZ();

        node1[0] = xx(Nodeindex1[0], Nodeindex1[1], Nodeindex1[2]);
        node1[1] = yy(Nodeindex1[0], Nodeindex1[1], Nodeindex1[2]);
        node1[2] = zz(Nodeindex1[0], Nodeindex1[1], Nodeindex1[2]);

        node2[0] = xx(Nodeindex2[0], Nodeindex2[1], Nodeindex2[2]);
        node2[1] = yy(Nodeindex2[0], Nodeindex2[1], Nodeindex2[2]);
        node2[2] = zz(Nodeindex2[0], Nodeindex2[1], Nodeindex2[2]);

        node3[0] = xx(Nodeindex3[0], Nodeindex3[1], Nodeindex3[2]);
        node3[1] = yy(Nodeindex3[0], Nodeindex3[1], Nodeindex3[2]);
        node3[2] = zz(Nodeindex3[0], Nodeindex3[1], Nodeindex3[2]);

        AA = (node2[1] - node1[1]) * (node3[2] - node1[2]) - (node2[2] - node1[2]) * (node3[1] - node1[1]);
        BB = (node2[2] - node1[2]) * (node3[0] - node1[0]) - (node2[0] - node1[0]) * (node3[2] - node1[2]);
        CC = (node2[0] - node1[0]) * (node3[1] - node1[1]) - (node2[1] - node1[1]) * (node3[0] - node1[0]);
        DD = 0 - AA * node1[0] - BB * node1[1] - CC * node1[2];

        vector <RDouble> FF;

        FF.push_back(AA * node1[0] + BB * node1[1] + CC * node1[2] + DD);
        FF.push_back(AA * node2[0] + BB * node2[1] + CC * node2[2] + DD);
        FF.push_back(AA * node3[0] + BB * node3[1] + CC * node3[2] + DD);

        int nTotalNode = grid->GetNTotalNode();
        RDouble *x = grid->GetX();
        RDouble *y = grid->GetY();
        RDouble *z = grid->GetZ();

        RDouble count = 0;
        RDouble resu = 0;
        for (int iNode = 0; iNode < nTotalNode; iNode ++)
        {
            resu = AA * x[iNode] + BB * y[iNode] + CC * z[iNode] + DD;
            if (ABS(count) < ABS(resu))
                count = resu;
        }
    }
}

void Pre_GridMirror::GridMirror()
{
    gridinfo->MoveToBegin();
    PHWrite(gridinfo, &BlockofEachProc, 1);

    for (int iZone = 0; iZone < BlockofEachProc; iZone ++)
    {
        if (gridtype == STRUCTGRID)
        {
            StructGrid *grid = StructGridCast(OrdinaryGrid[iZone]);
            MirrorGridStr(grid);
        }
        else
        {
            UnstructGrid *grid = UnstructGridCast(OrdinaryGrid[iZone]);
            MirrorGridUnstr(grid);
        }
    }
}

void Pre_GridMirror::MirrorGridUnstr(UnstructGrid *grid_in)
{
    int Gridindex = grid_in->GetGridID()->GetIndex();
    PHWrite(gridinfo, &Gridindex, 1);

    MirrorCoordinate(grid_in);
    MirrorNodeNumberOfEachFace(grid_in);
    MirrorFace2Node(grid_in);
    MirrorCellOfFace(grid_in);
    MirrorCell2Node(grid_in);
    //mirrorindex ++;
}

void Pre_GridMirror::MirrorCell2Node(UnstructGrid *grid_in)
{
    int nTotalCells = grid_in->GetNTotalCell();

    int cell2node = 1;
    Geo_CellTopo_Unstruct *cellTopology = grid_in->GetCellTopology();
    if (!cellTopology->GetNodeNumberOfEachCell())
    {
        PrintToWindow("Error: cell-node topology is empty!\n");
        cell2node = 0;
        PHWrite(gridinfo, &cell2node, 1);
        return;
    }
    
    PHWrite(gridinfo, &cell2node, 1);
    int *node_number_of_each_cell_ordinary = grid_in->GetNodeNumberOfEachCell();
    int *cell2node_ordinary = grid_in->GetCell2Node();

    PHWrite(gridinfo, node_number_of_each_cell_ordinary, nTotalCells);
    PHWrite(gridinfo, node_number_of_each_cell_ordinary, nTotalCells);

    int num = 0;
    for (int iCell = 0; iCell < nTotalCells; ++ iCell)
    {
        num += node_number_of_each_cell_ordinary[iCell];
    }

    PHWrite(gridinfo, &num, 1);
    PHWrite(gridinfo, cell2node_ordinary, num);
    PHWrite(gridinfo, cell2node_ordinary, num);
}

void Pre_GridMirror::MirrorCellOfFace(UnstructGrid *grid_in)
{
    int nTotalFace = grid_in->GetNTotalFace();
    int nBoundFace = grid_in->GetNBoundFace();

    int *left_cell_of_face_ordinary  = grid_in->GetLeftCellOfFace();
    int *right_cell_of_face_ordinary = grid_in->GetRightCellOfFace();

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        right_cell_of_face_ordinary[iFace] = -1;
    }

    PHWrite(gridinfo, left_cell_of_face_ordinary, nTotalFace);
    PHWrite(gridinfo, left_cell_of_face_ordinary, nTotalFace);
    PHWrite(gridinfo, right_cell_of_face_ordinary, nTotalFace);
    PHWrite(gridinfo, right_cell_of_face_ordinary, nTotalFace);
}

void Pre_GridMirror::MirrorFace2Node(UnstructGrid *grid_in)
{
    int nTotalFace = grid_in->GetNTotalFace();
    int *node_number_of_each_face_ordinary = grid_in->GetNodeNumberOfEachFace();
    int *face2node_ordinary = grid_in->GetFace2Node();

    int nsum = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        nsum += node_number_of_each_face_ordinary[iFace];
    }

    int *face2node_mir = new int [nsum];

    int count = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int nodeos = node_number_of_each_face_ordinary[iFace];

        for (int iNode = 0; iNode < nodeos; ++ iNode)
        {
            face2node_mir[count + iNode] = face2node_ordinary[count + nodeos - 1 - iNode];
        }
        count += nodeos;
    }

    PHWrite(gridinfo, &nsum, 1);
    PHWrite(gridinfo, face2node_ordinary, nsum);
    PHWrite(gridinfo, face2node_mir, nsum);

    delete [] face2node_mir;
}

void Pre_GridMirror::MirrorNodeNumberOfEachFace(UnstructGrid *grid_in)
{
    int nTotalFace = grid_in->GetNTotalFace();
    int *node_number_of_each_face_ordinary = grid_in->GetNodeNumberOfEachFace();

    PHWrite(gridinfo, node_number_of_each_face_ordinary, nTotalFace);
    PHWrite(gridinfo, node_number_of_each_face_ordinary, nTotalFace);
}

void Pre_GridMirror::MirrorCoordinate(UnstructGrid *grid_in)
{
    int nTotalNode = grid_in->GetNTotalNode();
    int nTotalFace = grid_in->GetNTotalFace();
    int nTotalCells = grid_in->GetNTotalCell();

    PHWrite(gridinfo, &nTotalNode, 1);
    PHWrite(gridinfo, &nTotalFace, 1);
    PHWrite(gridinfo, &nTotalCells, 1);

    RDouble *x_mir = new RDouble[nTotalNode];
    RDouble *y_mir = new RDouble[nTotalNode];
    RDouble *z_mir = new RDouble[nTotalNode];
    
    RDouble *x_ordinary = grid_in->GetX();
    RDouble *y_ordinary = grid_in->GetY();
    RDouble *z_ordinary = grid_in->GetZ();

    PHWrite(gridinfo, x_ordinary, nTotalNode);
    PHWrite(gridinfo, y_ordinary, nTotalNode);
    PHWrite(gridinfo, z_ordinary, nTotalNode);

    if (SymmetryFaceVector == X_axis)
    {
        for (int iNode = 0; iNode < nTotalNode; iNode ++)
        {
            x_mir[iNode] = - x_ordinary[iNode];
            y_mir[iNode] = y_ordinary[iNode];
            z_mir[iNode] = z_ordinary[iNode];
        }
    }
    else if (SymmetryFaceVector == Y_axis)
    {
        for (int iNode = 0; iNode < nTotalNode; iNode ++)
        {
            x_mir[iNode] = x_ordinary[iNode];
            y_mir[iNode] = - y_ordinary[iNode];
            z_mir[iNode] = z_ordinary[iNode];
        }
    }
    else
    {
        for (int iNode = 0; iNode < nTotalNode; iNode ++)
        {
            x_mir[iNode] = x_ordinary[iNode];
            y_mir[iNode] = y_ordinary[iNode];
            z_mir[iNode] = - z_ordinary[iNode];
        }
    }

    PHWrite(gridinfo, x_mir, nTotalNode);
    PHWrite(gridinfo, y_mir, nTotalNode);
    PHWrite(gridinfo, z_mir, nTotalNode);

    /*grid->ComputeMinMaxBox();
    Mirror->ComputeMinMaxBox();*/
}

void Pre_GridMirror::MirrorGridStr(StructGrid *grid_in)
{
    int Gridindex = grid_in->GetGridID()->GetIndex();
    PHWrite(gridinfo, &Gridindex, 1);

/*
    GridID *index = new GridID(Gridindex);
    Grid *grid;
    grid = CreateGridGeneral(STRUCTGRID, index, 0, Dimension);
    MirrorGrids.push_back(grid);
    StructGridCast(grid)->CopyGrid(StructGridCast(grid_in));*/
    /*

    GridID *index_mirror = new GridID(Gridindex + numberOfOriginalBlocks);
    Grid * Mirror;
    Mirror = CreateGridGeneral(STRUCTGRID, index_mirror, 0, Dimension);
    MirrorGrids.push_back(Mirror);*/

    int nTotalNode = grid_in->GetNTotalNode();
    PHWrite(gridinfo, &nTotalNode, 1);

    int ni = grid_in->GetNI();
    int nj = grid_in->GetNJ();
    int nk = grid_in->GetNK();
    PHWrite(gridinfo, &ni, 1);
    PHWrite(gridinfo, &nj, 1);
    PHWrite(gridinfo, &nk, 1);

/*
    StructGridCast(Mirror)->SetNI(ni);
    StructGridCast(Mirror)->SetNJ(nj);
    StructGridCast(Mirror)->SetNK(nk);
    StructGridCast(Mirror)->SetBasicDimension();*/

    RDouble *MirrorX = new RDouble [nTotalNode];
    RDouble *MirrorY = new RDouble [nTotalNode];
    RDouble *MirrorZ = new RDouble [nTotalNode];
/*
    Mirror->SetX(MirrorX);
    Mirror->SetY(MirrorY);
    Mirror->SetZ(MirrorZ);
    StructGridCast(Mirror)->SetArrayLayout();*/

    RDouble *x = grid_in->GetX();
    RDouble *y = grid_in->GetY();
    RDouble *z = grid_in->GetZ();
    PHWrite(gridinfo, x, nTotalNode);
    PHWrite(gridinfo, y, nTotalNode);
    PHWrite(gridinfo, z, nTotalNode);

    if (SymmetryFaceVector == X_axis)
    {
        for (int iNode = 0; iNode < nTotalNode; iNode ++)
        {
            MirrorX[iNode] = - x[iNode];
            MirrorY[iNode] = y[iNode];
            MirrorZ[iNode] = z[iNode];
        }
    }
    else if (SymmetryFaceVector == Y_axis)
    {
        for (int iNode = 0; iNode < nTotalNode; iNode ++)
        {
            MirrorX[iNode] = x[iNode];
            MirrorY[iNode] = - y[iNode];
            MirrorZ[iNode] = z[iNode];
        }
    }
    else
    {
        for (int iNode = 0; iNode < nTotalNode; iNode ++)
        {
            MirrorX[iNode] = x[iNode];
            MirrorY[iNode] = y[iNode];
            MirrorZ[iNode] = - z[iNode];
        }
    }

    PHWrite(gridinfo, MirrorX, nTotalNode);
    PHWrite(gridinfo, MirrorY, nTotalNode);
    PHWrite(gridinfo, MirrorZ, nTotalNode);

    /*StructGridCast(Mirror)->ComputeMinMaxBox();
    StructGridCast(Mirror)->CopyStructBCSet(grid_in->GetStructBCSet());*/
}

void Pre_GridMirror::MirrorBCRegion()
{
    for (int iZone = 0; iZone < BlockofEachProc; iZone ++)
    {
        if (gridtype == STRUCTGRID)
        {
            StructGrid *grid = StructGridCast(OrdinaryGrid[iZone]);
            /*StructGrid * Mirror = StructGridCast(MirrorGrids[iZone * 2 + 1]);*/

            StructBCSet *rhs = grid->GetStructBCSet();
            /*StructBCSet * structBCSet = Mirror->GetStructBCSet();*/

            int numberOfBoundaryFaces = rhs->GetnBCRegion();
            PHWrite(gridinfo, &numberOfBoundaryFaces, 1);

            for (int iBCRegion = 0; iBCRegion < numberOfBoundaryFaces; iBCRegion ++)
            {
                StructBC *bcregion = rhs->GetBCRegion(iBCRegion);
                /*StructBC *MirrorBCRegion = structBCSet->GetBCRegion(iBCRegion);*/

                int bctype = bcregion->GetBCType();

                /*if (bctype == PHENGLEI::INTERFACE)
                {
                    int nbt = bcregion->GetTargetRegionBlock();

                    MirrorBCRegion->SetTargetRegionBlock(nbt + numberOfOriginalBlocks);

                }*/

                if (bctype == PHENGLEI::SYMMETRY)
                {
                    bctype = PHENGLEI::INTERFACE;

                    /*bcregion->SetBCType(PHENGLEI::INTERFACE);
                    int imin, imax, jmin, jmax, kmin, kmax;
                    bcregion->GetIJKRegion(imin,imax,jmin,jmax,kmin,kmax);
                    bcregion->SetTargetIJKRegion(imin,imax,jmin,jmax,kmin,kmax);*/

                    /*int Gridindex = grid->GetGridID()->GetIndex();
                    bcregion->SetTargetRegionBlock(Gridindex + numberOfOriginalBlocks);*/

                    /*MirrorBCRegion->SetBCType(PHENGLEI::INTERFACE);
                    MirrorBCRegion->GetIJKRegion(imin,imax,jmin,jmax,kmin,kmax);
                    MirrorBCRegion->SetTargetIJKRegion(imin,imax,jmin,jmax,kmin,kmax);
                    MirrorBCRegion->SetTargetRegionBlock(Gridindex);*/
                }

                PHWrite(gridinfo, &bctype, 1);
                int *s_st = bcregion->GetStartPoint();
                int *s_ed = bcregion->GetEndPoint();
                PHWrite(gridinfo, s_st, 3);
                PHWrite(gridinfo, s_ed, 3);
                
                PHWrite(gridinfo, s_st, 3);
                PHWrite(gridinfo, s_ed, 3);
                /*int * s_st_mir = MirrorBCRegion->GetStartPoint();
                int * s_ed_mir = MirrorBCRegion->GetEndPoint();
                PHWrite(gridinfo, s_st_mir, 3);
                PHWrite(gridinfo, s_ed_mir, 3);*/

                /*if (IsInterface(bctype))
                {
                    int *t_sd = bcregion->GetTargetStart();
                    int *t_ed = bcregion->GetTargetEnd();
                    PHWrite(gridinfo, t_sd, 3);
                    PHWrite(gridinfo, t_ed, 3);

                    int nbt = bcregion->GetTargetRegionBlock();
                    PHWrite(gridinfo, &nbt, 1);

                    int *t_sd_mir = MirrorBCRegion->GetTargetStart();
                    int *t_ed_mir = MirrorBCRegion->GetTargetEnd();
                    PHWrite(gridinfo, t_sd_mir, 3);
                    PHWrite(gridinfo, t_ed_mir, 3);

                    int nbt_mir = MirrorBCRegion->GetTargetRegionBlock();
                    PHWrite(gridinfo, &nbt_mir, 1);

                }*/
            }
        }
        else if (gridtype == UNSTRUCTGRID)
        {
            UnstructGrid *grid_ordinary = UnstructGridCast(OrdinaryGrid[iZone]);


            UnstructBCSet **bcr_ordinary = grid_ordinary->GetBCRecord();
            int nBoundFace = grid_ordinary->GetNBoundFace();
            PHWrite(gridinfo, &nBoundFace, 1);

            int count = 0;
            for (int iFace = 0; iFace < nBoundFace; ++ iFace)
            {
                if (bcr_ordinary[iFace]->GetKey() == PHENGLEI::INTERFACE)
                {
                    count ++;
                }

                if (bcr_ordinary[iFace]->GetKey() == PHENGLEI::SYMMETRY)
                {
                    count ++;
                }
            }
           
            PHWrite(gridinfo, &count, 1);

            InterfaceInfo *iinfo_ordinary = grid_ordinary->GetInterfaceInfo();
            /*InterfaceInfo * interfaceInfo = 0;
            InterfaceInfo * iinfo_mir = 0;
            if (count != 0)
            {
                interfaceInfo     = new InterfaceInfo(count);
                iinfo_mir = new InterfaceInfo(count);
            }

            int * interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
            int * interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();
            int * interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace(); 
            int * interFaceDirection     = interfaceInfo->GetInterFaceDirection();

            int * interFace2ZoneID_mir       = iinfo_mir->GetInterFace2ZoneID();
            int * interFace2InterFaceID_mir  = iinfo_mir->GetInterFace2InterFaceID();
            int * interFace2BoundaryFace_mir = iinfo_mir->GetInterFace2BoundaryFace(); 
            int * interFaceDirection_mir     = iinfo_mir->GetInterFaceDirection();*/

            int * interFace2ZoneID_ordinary       = NULL;
            int * interFace2InterFaceID_ordinary  = NULL;
            int * interFace2BoundaryFace_ordinary = NULL;
            //int * interFaceDirection_ordinary;

            if (numberOfZones > 2)
            {
                interFace2ZoneID_ordinary       = iinfo_ordinary->GetInterFace2ZoneID();
                interFace2InterFaceID_ordinary  = iinfo_ordinary->GetInterFace2InterFaceID();
                interFace2BoundaryFace_ordinary = iinfo_ordinary->GetInterFace2BoundaryFace();
                //interFaceDirection_ordinary     = iinfo_ordinary->GetInterFaceDirection();
            }

            /*int *interFace2ZoneID_ordinary = iinfo_ordinary->GetInterFace2ZoneID();
            int *interFace2InterFaceID_ordinary = iinfo_ordinary->GetInterFace2InterFaceID();
            int *interFace2BoundaryFace_ordinary = iinfo_ordinary->GetInterFace2BoundaryFace();
            int *interFaceDirection_ordinary = iinfo_ordinary->GetInterFaceDirection();*/

            count = 0;
            int iInterfae = 0;

            for (int iFace = 0; iFace < nBoundFace; ++ iFace)
            {
                int bctype = bcr_ordinary[iFace]->GetKey();
                string bcName = bcr_ordinary[iFace]->GetBCName();

                if (bctype == PHENGLEI::INTERFACE)
                {
                    PHWrite(gridinfo, &bctype, 1);
                    gridinfo->WriteString(bcName);

                    int interFace2BoundaryFace = interFace2BoundaryFace_ordinary[iInterfae];
                    int interFace2ZoneID       = interFace2ZoneID_ordinary[iInterfae];
                    int interFace2InterFaceID  = interFace2InterFaceID_ordinary[iInterfae];

                    PHWrite(gridinfo, &interFace2BoundaryFace, 1);
                    PHWrite(gridinfo, &interFace2ZoneID, 1);
                    PHWrite(gridinfo, &interFace2InterFaceID, 1);

                    interFace2ZoneID += numberOfOriginalBlocks;

                    PHWrite(gridinfo, &interFace2BoundaryFace, 1);
                    PHWrite(gridinfo, &interFace2ZoneID, 1);
                    PHWrite(gridinfo, &interFace2InterFaceID, 1);

                    iInterfae ++;
                    count ++;
                    continue;
                }

                if (bctype == PHENGLEI::SYMMETRY)
                {
                    int Gridindex = grid_ordinary->GetGridID()->GetIndex();
                    bctype = PHENGLEI::INTERFACE;
                    bcName = "Interface";
                    PHWrite(gridinfo, &bctype, 1);
                    gridinfo->WriteString(bcName);

                    int interFace2BoundaryFace = iFace;
                    int interFace2ZoneID       = Gridindex + numberOfOriginalBlocks;
                    int interFace2InterFaceID  = count;

                    PHWrite(gridinfo, &interFace2BoundaryFace, 1);
                    PHWrite(gridinfo, &interFace2ZoneID, 1);
                    PHWrite(gridinfo, &interFace2InterFaceID, 1);

                    interFace2ZoneID = Gridindex;
                    PHWrite(gridinfo, &interFace2BoundaryFace, 1);
                    PHWrite(gridinfo, &interFace2ZoneID, 1);
                    PHWrite(gridinfo, &interFace2InterFaceID, 1);

                    /*interFace2ZoneID[count] = Gridindex + numberOfOriginalBlocks;
                    interFace2InterFaceID[count] = count;
                    interFace2BoundaryFace[count] = iFace;

                    interFace2ZoneID_mir[count] = Gridindex;
                    interFace2InterFaceID_mir[count] = count;
                    interFace2BoundaryFace_mir[count] = iFace;*/

                    //iInterfae ++;
                    count ++;
                    continue;
                }

                PHWrite(gridinfo, &bctype, 1);
                gridinfo->WriteString(bcName);
            }
        }
    }
}

}
