#include "Geo_UnstructBC.h"
#include "IO_FileName.h"
#include "Pre_HDF5File.h"
#include "Glb_Dimension.h"
#include "LinkStruct.h"
#include "Post_Visual.h"
#include "Pre_GridConversion.h"
#include "HyList.h"
#include "PHIO.h"
#include "TK_Log.h"

using namespace std;

namespace PHSPACE
{
LIB_EXPORT Pre_GridConversion::Pre_GridConversion(const string &gridFileName)
{
    this->gridFileName = gridFileName;
    grids = 0;
    nBlocks = 0;
}

LIB_EXPORT Pre_GridConversion::~Pre_GridConversion()
{
    
}

LIB_EXPORT void Pre_GridConversion::Clear()
{
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        delete grids[iZone];
    }
    delete [] grids;
}

LIB_EXPORT void Pre_GridConversion::Run()
{
    ReadGrid();

    Conversion();
}

void Pre_GridConversion::WriteVolumeName(const string &out_grid_file, Grid **grids_in)
{
    int hasVolumeCondition = GlobalDataBase::GetIntParaFromDB("hasVolumeCondition");
    if (!hasVolumeCondition)
    {
        return;
    }

    PrintToWindow("    Writing volume name to .vcname file ....\n");

    set< pair<int, string> > vcNamePairSet;
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        SimpleVC *volumeCondition = grids_in[iZone]->GetVolumeConditionIn();

        int vcType = volumeCondition->GetVCType();
        string vcName = volumeCondition->GetVCName();
        pair<int, string> vcNamePair(vcType, vcName);
        vcNamePairSet.insert(vcNamePair);
    }

    DumpVCNameFile(vcNamePairSet, out_grid_file);
    DumpVolumeConditionFile(vcNamePairSet);
}

void Pre_GridConversion::WriteFaceBoundaryName(const string &out_grid_file, Grid **grids_in)
{
    PrintToWindow("    Writing face boundary name to .bcname file ....\n");
    Grid **grids = grids_in;
    string from_gfile = "grid.fts";
    GlobalDataBase::GetData("from_gfile", &from_gfile, PHSTRING, 1);

    set< pair<int, string> > bcNamePairSet;
    //int nBlocks = 1;
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        UnstructGrid *grid = UnstructGridCast(grids[iZone]);
        int nBoundFace = grid->GetNBoundFace();

        UnstructBCSet **bcr = grid->GetBCRecord();
        for (cgsize_t iFace = 0; iFace < nBoundFace; iFace ++)
        {
            int bcType = bcr[iFace]->GetKey();
            if (bcType == PHENGLEI::INTERFACE)
            {
                continue;
            }

            const string &bcName = bcr[iFace]->GetBCName();
            pair<int, string> bcNamePair(bcType, bcName);
            bcNamePairSet.insert(bcNamePair);
        }
    }

    DumpBoundaryInfo(bcNamePairSet);
    DumpBCNameInfo(bcNamePairSet, out_grid_file);
}

void Pre_GridConversion::WriteCellToNode(const string &targetGridFileName, Grid **grids_in)
{
    fstream file;

    Grid **grids = grids_in;

    string cellToNodeFileName = ChangeExtensionOfFileName(targetGridFileName, "c2n");

    PHSPACE::OpenSerialFile(file, cellToNodeFileName, ios_base::out|ios_base::binary);

    int nBlocks = 1;
    PHWrite(file, nBlocks);

    VirtualFile *virtualFile = new VirtualFile(&file);
    virtualFile->BeginWriteWork();

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        PHWrite(virtualFile, iZone);

        int nTotalCell = grids[0]->GetNTotalCell();
        PHWrite(virtualFile, nTotalCell);

        int *cell2Node = UnstructGridCast(grids[0])->GetCell2Node();
        int *node_number_of_each_cell = UnstructGridCast(grids[0])->GetNodeNumberOfEachCell();

        int nPoints;
        int count = 0;
        for (int iCell = 0; iCell < nTotalCell; iCell ++)
        {
            nPoints = node_number_of_each_cell[iCell];
            PHWrite(virtualFile, nPoints);

            for (int jPoint = 0; jPoint < nPoints; ++ jPoint)
            {
                int pointIndex = cell2Node[count];
                PHWrite(virtualFile, pointIndex);
                count ++;
            }
        }
    }

    virtualFile->EndWriteWork();
    delete virtualFile;

    PHSPACE::CloseFile(file);
}

void Pre_GridConversion::WriteFaceBC(const string &out_grid_file, Grid **grids_in)
{
    fstream file;

    Grid **grids = grids_in;

    string BCFileName = ChangeExtensionOfFileName(out_grid_file, "bc");

    PHSPACE::OpenSerialFile(file, BCFileName, ios_base::out|ios_base::binary);

    //int nBlocks = 1;
    PHWrite(file, nBlocks);

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        UnstructGrid *grid = UnstructGridCast(grids[iZone]);

        grid->WriteFaceBC(file);
    }

    PHSPACE::CloseFile(file);
}

LIB_EXPORT void DumpGrid(const string &gridFileName, Grid **grids, int nBlocks)
{
    int dumpOldGrid = GlobalDataBase::GetIntParaFromDB("dumpOldGrid");
    if (dumpOldGrid)
    {
        DumpOldGrid(gridFileName, grids, nBlocks);
        return;
    }

    DumpHDF5Grid(gridFileName, grids, nBlocks);
}

LIB_EXPORT void DumpOldGrid(const string &gridFileName, Grid **grids, int nBlocks)
{
    PrintToWindow("    Writing old version grid to file ....\n");
    fstream file;
    OpenSerialFile(file, gridFileName.c_str(), ios_base::out|ios_base::binary|ios_base::trunc);

    using namespace PHMPI;

    SetNumberOfGlobalZones(nBlocks);

    int nZones = GetNumberofGlobalZones();

    CreateZoneProcessorID(nZones);
    CreateZoneGridID(nZones);
    CreateZoneGridType(nZones);

    int *block_proc     = GetZoneProcessorID();
    int *block_idx      = GetZoneGridID();
    int *block_type     = GetZoneGridType();
    int *block_proc_inp = GetZoneProcessorID_INP();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        block_proc[iZone] = iZone;
        block_idx [iZone] = iZone;
        block_type[iZone] = grids[iZone]->Type();
    }

    if (block_proc_inp)
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            SetZoneProcessorID(iZone, block_proc_inp[iZone]);
        }
    }

    file.write(reinterpret_cast<char *>(&nBlocks), sizeof(int));

    int m_block_proc = PHMPI::GetZoneDistributionMethod();
    if (m_block_proc == DETERMIN_BY_PARTITION && IsConvertGridToMixGrid())
    {
        int *block_proc_dump = GetZoneProcessorIDDump();
        file.write(reinterpret_cast<char *>(block_proc_dump), nZones*sizeof(int));
    }
    else
    {
        file.write(reinterpret_cast<char *>(block_proc), nZones*sizeof(int));
    }

    file.write(reinterpret_cast<char *>(block_idx), nZones*sizeof(int));
    file.write(reinterpret_cast<char *>(block_type), nZones*sizeof(int));

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        grids[iZone]->WriteGrid(file);
    }

    file.close();
    file.clear();
}

void GetGridsBoundingBox(Grid **grids, int nBlocks, RDouble *pmin, RDouble *pmax)
{
    pmin[0] = LARGE;
    pmin[1] = LARGE;
    pmin[2] = LARGE;

    pmax[0] = - LARGE;
    pmax[1] = - LARGE;
    pmax[2] = - LARGE;

    for (int nb = 0; nb < nBlocks; ++ nb)
    {
        RDouble *local_pmin = grids[nb]->GetMinBox();
        RDouble *local_pmax = grids[nb]->GetMaxBox();

        pmin[0] = MIN(pmin[0], local_pmin[0]);
        pmin[1] = MIN(pmin[1], local_pmin[1]);
        pmin[2] = MIN(pmin[2], local_pmin[2]);

        pmax[0] = MAX(pmax[0], local_pmax[0]);
        pmax[1] = MAX(pmax[1], local_pmax[1]);
        pmax[2] = MAX(pmax[2], local_pmax[2]);
    }
}

void GetGridsMinMaxDS(Grid **grids, int nBlocks, RDouble &mindis, RDouble &maxdis)
{
    mindis =   LARGE;
    maxdis = - LARGE;

    for (int nb = 0; nb < nBlocks; ++ nb)
    {
        RDouble dismin, dismax;
        grids[nb]->GetMinMaxDS(dismin, dismax);

        mindis = MIN(mindis, dismin);
        maxdis = MAX(maxdis, dismax);
    }
}

void MatchInterface(Grid *grid, LinkStruct *link)
{
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    int nIFace    = interfaceInfo->GetNIFace();
    int zoneIndex = grid->GetZoneID();
    int npeoridic = 0;

    VVInt &zoneid  = link->GetZoneID();
    VVInt &faceid  = link->GetFaceID();
    VVInt &facemap = link->GetFaceMap();

    int *interFace2ZoneID      = interfaceInfo->GetInterFace2ZoneID();
    int *interFace2InterFaceID = interfaceInfo->GetInterFace2InterFaceID();

    for (int iiface = 0; iiface < nIFace; ++ iiface)
    {
        //! image is interface mapped from i to the global,the interface is for all blocks,,not only the local block.
        int ipos = facemap[zoneIndex][ iiface ];
        //int find_flag = false;
        bool find_flag = false;     //! Corrected by clz, 2012.3.15.

        bool tmpout = false;

        if (zoneid[ipos].size() != 2)
        {
            cout << " more than or less than two faces coincide\n";
            cout << " zoneIndex= " << zoneIndex << "\n";
            cout << " zoneid[ipos].size() = " << zoneid[ipos].size() << "\n";
            cout << " iiface = " << iiface << " nIFace = " << nIFace << "\n";
            tmpout = true;
        }

        for (std::size_t m = 0; m < zoneid[ipos].size(); ++ m)
        {
            int iZone = zoneid[ipos][m];
            int iFace = faceid[ipos][m];
            if ((iZone != zoneIndex) || (iFace != iiface))
            {
                interFace2ZoneID[iiface]      = iZone;
                interFace2InterFaceID[iiface] = iFace;
                find_flag = true;
                break;
            }
        }

        if (!find_flag)
        {
            //! It shows that there is  non-physical artificial periodic boundary condition in this block.
            //! Of course,there may be problems in the process.
            cout << "There is an error in connection between the calculation area of this block and the butt block " << 
                    "(There may also be non-physical artificial periodic boundary conditions), please check the grid!\n";
            ++ npeoridic;
        }
    }
}

void ConstructGlobalInterfaceLink(Grid **grids, int nBlocks)
{
    RDouble mindis, maxdis, pmin[3], pmax[3];

    GetGridsBoundingBox(grids, nBlocks, pmin, pmax);
    GetGridsMinMaxDS(grids, nBlocks, mindis, maxdis);

    cout << "  Min && Max edge length: " << mindis << ", " << maxdis << "\n";

    typedef DataStruct_AdtTree<int, RDouble> AdtTree;

    if (GetDim() == TWO_D)
    {
        pmin[2] -= mindis;
        pmax[2] += maxdis;
    }

    LinkStruct link(pmin, pmax, mindis / 2.0, nBlocks);

    int omit_no_bound_bc = 1;
    GlobalDataBase::GetData("omit_no_bound_bc", &omit_no_bound_bc, PHINT, 1);

    if (!omit_no_bound_bc)
    {
        //! Covert NO_BOUNDARY_CONDITION to INTERFACE.
        for (int nb = 0; nb < nBlocks; ++ nb)
        {
            grids[nb]->ChangeBCType(PHENGLEI::NO_BOUNDARY_CONDITION, PHENGLEI::INTERFACE);

            //! Consider the periodic boundary as the interface.
            if (PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType") != NO_PERIODICITY)
            {
                grids[nb]->ChangeBCType(PHENGLEI::USER_DEFINED, PHENGLEI::INTERFACE);
            }
        }
    }

    //! Build the initial link information.
    for (int nb = 0; nb < nBlocks; ++ nb)
    {
        grids[nb]->BuildGridLink(&link);
    }

    //! Further process.
    for (int nb = 0; nb < nBlocks; ++ nb)
    {
        MatchInterface(grids[nb], &link);
    }
}

LIB_EXPORT void DumpAdditionalData(const string &gridFileName, Grid **grids, int nBlocks)
{
    DumpFaceBC(gridFileName, grids, nBlocks);
}

LIB_EXPORT void DumpFaceBC(const string &gridFileName, Grid **grids, int nBlocks)
{    
    string faceBCFileName = ChangeExtensionOfFileName(gridFileName, "bc");
    faceBCFileName = AddSymbolToFileName(faceBCFileName, '_', PHMPI::GetCurrentProcessorID());
    fstream file;
    OpenFile(file, faceBCFileName, ios_base::out|ios_base::binary);

    PHWrite(file, nBlocks);
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        Grid *grid = grids[iZone];

        grid->WriteFaceBC(file);
    }

    CloseFile(file);
}

void DumpCharacteristicBoundary(const string &gridFileName, Grid **grids, int nBlocks)
{
    string filename = ChangeExtensionOfFileName(gridFileName, "bcmesh");
    filename = AddSymbolToFileName(filename, "_", 0);
    vector < DataContainer * > datalist;
    datalist.resize(0);

    for (int iZone =  0; iZone < nBlocks; ++ iZone)
    {
        DataContainer *data = new DataContainer();
        data->MoveToBegin();

        Grid *grid = grids[iZone];
        grid->AllocateOversetGrid();
        ComputeCharacteristicBoundary(grid, data);

        datalist.push_back(data);
    }

    bool CharacteristicBoundary = true;

    DumpToVTK DumpToVTK(datalist, filename, CharacteristicBoundary);
    DumpToVTK.Run();
}

void WriteBCNameInfo(const string &gridFileName, Grid **grids_in, int nBlocks)
{
    Grid **grids = grids_in;

    set< pair<int, string> > bcNamePairSet;
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        grids[iZone]->GetBcNamePair(bcNamePairSet);
    }

    DumpBCNameInfo(bcNamePairSet, gridFileName);
}

void DumpBCNameInfo(const set< pair<int, string> > &bcNamePairSet, const string &gridFileName)
{
    fstream file;
    string BCNameFileName = ChangeExtensionOfFileName(gridFileName, "bcname");
    PHSPACE::OpenSerialFile(file, BCNameFileName, ios_base::out);

    file << "# Boundary Face Information of Grid " << gridFileName << endl;
    file << "# nBoundaryConditions: number of global boundary conditions." << endl;
    file << "# bcName             : Boundary Condition Name." << endl;
    file << "# bcType(in PHengLEI): Boundary Condition Type." << endl;
    file << endl;
    file << "int nBoundaryConditions = " << bcNamePairSet.size() << ";" << endl;

    set< pair<int, string> >::iterator iter;
    for (iter = bcNamePairSet.begin(); iter != bcNamePairSet.end(); ++ iter)
    {
        int bcType    = (*iter).first;
        string bcName = (*iter).second;
        file << "string bcName = \"" << bcName << "\";" << endl;
        file << "{" << endl;
        if (bcType == 2)
        {
            file << "  string bodyName = " << "\"" << "body" << "\";" << endl;
        }
        file << "  int bcType = " << bcType << ";" << endl;
        file << "}" << endl;
    }
    PHSPACE::CloseFile(file);
}

void DumpVCNameFile(const set< pair<int, string> > &vcNamePairSet, const string &gridFileName)
{
    fstream file;
    string VCNameFileName = ChangeExtensionOfFileName(gridFileName, "vcname");
    PHSPACE::OpenSerialFile(file, VCNameFileName, ios_base::out);

    file << "# Volume Condition Information of Grid " << gridFileName << endl;
    file << "# nVolumeConditions: number of global volume conditions." << endl;
    file << "# vcName             : Volume Condition Name." << endl;
    file << "# vcType(in PHengLEI): Volume Condition Type." << endl;
    file << endl;
    file << "int nVolumeConditions = " << vcNamePairSet.size() << ";" << endl;

    set< pair<int, string> >::iterator iter;
    for (iter = vcNamePairSet.begin(); iter != vcNamePairSet.end(); ++ iter)
    {
        int vcType    = (*iter).first;
        string vcName = (*iter).second;
        file << "string vcName = \"" << vcName << "\";" << endl;
        file << "{" << endl;
        file << "  int vcType = " << vcType << ";" << endl;
        file << "}" << endl;
    }

    PHSPACE::CloseFile(file);
}

void DumpVolumeConditionFile(const set< pair<int, string> > &vcNamePairSet)
{
    fstream file;
    string VCFileName = "./bin/volume_condition.hypara";
    PHSPACE::OpenFile(file, VCFileName, ios_base::out);

    file << "# nVolumeConditions: number of global volume conditions." << endl;
    file << "# vcName             : Volume Condition Name." << endl;
    file << "# vcType(in PHengLEI): Volume Condition Type." << endl;
    file << endl;
    file << "int nVolumeConditions = " << vcNamePairSet.size() << ";" << endl;

    set< pair<int, string> >::iterator iter;
    for (iter = vcNamePairSet.begin(); iter != vcNamePairSet.end(); ++ iter)
    {
        int vcType    = (*iter).first;
        string vcName = (*iter).second;
        file << "string vcName = \"" << vcName << "\";" << endl;
        file << "{" << endl;
        file << "  int vcType = " << vcType << ";" << endl;
        file << "}" << endl;
    }

    PHSPACE::CloseFile(file);
}

void WriteBoundaryInfo(Grid **grids_in, int nBlocks)
{
    Grid **grids = grids_in;

    set< pair<int, string> > bcNamePairSet;
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        grids[iZone]->GetBcNamePair(bcNamePairSet);
    }

    DumpBoundaryInfo(bcNamePairSet);
}
void WriteVolumeInfo(Grid** grids_in, int nBlocks)
{
    set< pair<int, string> > vcNamePairSet;
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        SimpleVC *volumeCondition = grids_in[iZone]->GetVolumeConditionIn();

        int vcType = volumeCondition->GetVCType();
        string vcName = volumeCondition->GetVCName();
        pair<int, string> vcNamePair(vcType, vcName);
        vcNamePairSet.insert(vcNamePair);
    }
    DumpVolumeConditionFile(vcNamePairSet);
}
void DumpBoundaryInfo(const set< pair<int, string> > &bcNamePairSet)
{
    fstream file;
    string boundaryconditionfile = "./bin/boundary_condition.hypara";
    file.open(boundaryconditionfile.c_str(), ios_base::out);

    if (!file)
    {
        TK_Exit::FileOpenErrorExit(boundaryconditionfile);
    }

    //file << "# Boundary Face Information of Grid " << gridFileName << endl;
    file << "# nBoundaryConditions: Number of global boundary conditions." << endl;
    file << "# bcName             : Boundary condition name." << endl;
    file << "# bcType(in PHengLEI): Boundary condition type." << endl;
    file << endl;

    file << "# How to set boundary condition, for example: " << endl;
    file << "# string bcName = \"Wall\";                   " << endl;
    file << "# {                                           " << endl;
    file << "#   int bcType = 2;                           " << endl;
    file << "#   int viscousType = 1;                      " << endl;
    file << "#   double wallTemperature = -1.0;            " << endl;
    file << "#   double uWall = 0.0;                       " << endl;
    file << "#   double vWall = 0.0;                       " << endl;
    file << "#   double wWall = 0.0;                       " << endl;
    file << "# }                                           " << endl;
    file << "# string bcName = \"Inflow\";                 " << endl;
    file << "# {                                           " << endl;
    file << "#   int bcType = 5;                           " << endl;
    file << "#   int inflowParaType = 0;                   " << endl;
    file << "#   double refMachNumber = 0.73;              " << endl;
    file << "#   double attackd = 2.79;                    " << endl;
    file << "#   double angleSlide = 0.0;                  " << endl;
    file << "#   double refReNumber = 6.5e6;               " << endl;
    file << "#   double refDimensionalTemperature = 288.15;" << endl;
    file << "# }                                           " << endl;
    file << endl;
    file << "# For more information, see examples/bin/boundary_condition.hypara file!!!" << endl;
    file << endl;

    file << "int nBoundaryConditions = " << bcNamePairSet.size() << ";" << endl;

    set< pair<int, string> >::iterator iter;
    for (iter = bcNamePairSet.begin(); iter != bcNamePairSet.end(); ++ iter)
    {
        int bcType    = (*iter).first;
        string bcName = (*iter).second;
        if (bcName == "Connection") continue;

        file << "string bcName = \"" << bcName << "\";" << endl;
        file << "{" << endl;
        if (bcType == 2)
        {
            file << "  string bodyName = " << "\"" << "body" << "\";" << endl;
        }
        file << "  int bcType = " << bcType << ";" << endl;
        file << "}" << endl;
    }
    file << endl;

    file << "# 'bcType' is defined as following:" << endl;
    file << "#    -2: WAKE                      " << endl;
    file << "#    -1: INTERFACE                 " << endl;
    file << "#    0 : NO_BOUNDARY_CONDITION     " << endl;
    file << "#    1 : EXTRAPOLATION             " << endl;
    file << "#    2 : SOLID_SURFACE             " << endl;
    file << "#    3 : SYMMETRY                  " << endl;
    file << "#    4 : FARFIELD                  " << endl;
    file << "#    5 : INFLOW                    " << endl;
    file << "#    6 : OUTFLOW                   " << endl;
    file << "#    52: PRESSURE_INLET            " << endl;
    file << "#    62: PRESSURE_OUTLET           " << endl;
    file << "#    61: OUTFLOW_CONFINED          " << endl;
    file << "#    7 : POLE                      " << endl;

    PHSPACE::CloseFile(file);
}

void ComputeCharacteristicBoundary(Grid *grid_in, DataContainer *data)
{
    int GridID = grid_in->GetGridID()->GetIndex();
    PHWrite(data, GridID);

    int dimension = grid_in->GetDim();
    PHWrite(data, dimension);

    int type = grid_in->Type();
    PHWrite(data, &type, 1);

    if (type == UNSTRUCTGRID)
    {
        ComputeCharacteristicBoundaryUns(grid_in, data);
    }
    else if (type == STRUCTGRID)
    {
        ComputeCharacteristicBoundaryStr(grid_in, data);
    }

    int TecioMission = Writecomplete;
    PHWrite(data, &TecioMission, 1);
}

void ComputeCharacteristicBoundaryStr(Grid *grid_in, DataContainer *data)
{
    StructGrid *grid = StructGridCast(grid_in);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int s_st[3], s_ed[3];

    RDouble3D &xx = * grid->GetStructX();
    RDouble3D &yy = * grid->GetStructY();
    RDouble3D &zz = * grid->GetStructZ();

    StructBCSet *structBCSet = grid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();

    int TecioMission;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
        int bctype = bcregion->GetBCType();

        if (IsInterface(bctype)) continue;

        TecioMission = WriteBoundary;
        PHWrite(data, &TecioMission, 1);

        PHWrite(data, &bctype, 1);

        string bcName = bcregion->GetBCName();
        data->WriteString(bcName);

        int *s_lr3d = bcregion->GetFaceDirectionIndex();
        int nDim = GetDim();

        for (int m = 0; m < nDim; ++ m)
        {
            s_st[m] = bcregion->GetStartPoint(m);
            s_ed[m] = bcregion->GetEndPoint(m);

            if (bcregion->GetFaceDirection() == m)
            {
                if (s_st[m] > 1 || s_lr3d[m] == 1)
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

        int ist = s_st[0];
        int ied = s_ed[0];
        int jst = s_st[1];
        int jed = s_ed[1];
        int kst = s_st[2];
        int ked = s_ed[2];

        Range I(ist, ied);
        Range J(jst, jed);
        Range K(kst, ked);

        if (nDim != THREE_D)
        //if (nk == 1) 
        {
            K.setRange(1, 1);
            kst = 1;
            ked = 1;
        }

        Int3D newindex(I, J, K, fortranArray);

        int nsurf = bcregion->GetFaceDirection() + 1;

        int count = 0;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    newindex(i, j, k) = count;
                    ++ count;
                }
            }
        }

        int nTotalNode = count;
        int nTotalCell = 1;

        //for (int m = 0; m < nDim; ++ m)
        for (int m = 0; m < 3; ++ m)
        {
            nTotalCell *= (bcregion->GetEndPoint(m) - bcregion->GetStartPoint(m) + 1);
        }

        PHWrite(data, &nTotalNode, 1);
        PHWrite(data, &nTotalCell, 1);


        RDouble *xBoundary = new RDouble [nTotalNode];
        RDouble *yBoundary = new RDouble [nTotalNode];
        RDouble *zBoundary = new RDouble [nTotalNode];
                
        int iNode = 0;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    xBoundary[iNode] = xx(i, j, k);
                    yBoundary[iNode] = yy(i, j, k);
                    zBoundary[iNode] = zz(i, j, k);
                    ++ iNode;
                }
            }
        }

        PHWrite(data, xBoundary, nTotalNode);
        PHWrite(data, yBoundary, nTotalNode);
        PHWrite(data, zBoundary, nTotalNode);

        delete [] xBoundary;
        delete [] yBoundary;
        delete [] zBoundary;

        int il1, jl1, kl1;
        GetNsurfIndex(nsurf, il1, jl1, kl1);

        il1 = 1 - il1;
        jl1 = 1 - jl1;
        kl1 = 1 - kl1;

        if (nk == 1) kl1 = 0;

        int *cell2node = new int [nTotalCell * 4];
        int ncount = 0;
        for (int k = kst; k <= ked - kl1; ++ k)
        {
            for (int j = jst; j <= jed - jl1; ++ j)
            {
                for (int i = ist; i <= ied - il1; ++ i)
                {
                    int p1, p2, p3, p4;

                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    if (nsurf == 1)
                    {
                        p1 = newindex(i, j , k) + 1;
                        p2 = newindex(i, jl, k) + 1;
                        p3 = newindex(i, jl, kl) + 1;
                        p4 = newindex(i, j , kl) + 1;
                    }
                    else if (nsurf == 2)
                    {
                        p1 = newindex(i , j, k) + 1;
                        p2 = newindex(il, j, k) + 1;
                        p3 = newindex(il, j, kl) + 1;
                        p4 = newindex(i , j, kl) + 1;
                    }
                    else
                    {
                        p1 = newindex(i , j , k) + 1;
                        p2 = newindex(il, j , k) + 1;
                        p3 = newindex(il, jl, k) + 1;
                        p4 = newindex(i , jl, k) + 1;
                    }

                    cell2node[ncount++] = p1;
                    cell2node[ncount++] = p2;
                    cell2node[ncount++] = p3;
                    cell2node[ncount++] = p4;
                }
            }
        }

        PHWrite(data, cell2node, nTotalCell * 4);
        delete [] cell2node;
    }

    if (grid_in->GetDim() == PHSPACE::TWO_D)
    {
        /*int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();*/

        /*RDouble3D & xx = * grid->GetStructX();
        RDouble3D & yy = * grid->GetStructY();
        RDouble3D & zz = * grid->GetStructZ();*/

        int nTotalCell = grid->GetNTotalCell();
        int nTotalNode = grid->GetNTotalNode();

        TecioMission = WriteBlock;
        PHWrite(data, &TecioMission, 1);

        PHWrite(data, &nTotalNode, 1);
        PHWrite(data, &nTotalCell, 1);

        RDouble *xBoundary = new RDouble [nTotalNode];
        RDouble *yBoundary = new RDouble [nTotalNode];
        RDouble *zBoundary = new RDouble [nTotalNode];
    
        int ist, ied, jst, jed, kst, ked;
        grid->GetNodeIterationIndex(ist, ied, jst, jed, kst, ked);

        int iNode = 0;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    xBoundary[iNode] = xx(i, j, k);
                    yBoundary[iNode] = yy(i, j, k);
                    zBoundary[iNode] = zz(i, j, k);

                    iNode ++ ;
                }
            }
        }

        PHWrite(data, xBoundary, nTotalNode);
        PHWrite(data, yBoundary, nTotalNode);
        PHWrite(data, zBoundary, nTotalNode);

        delete [] xBoundary;
        delete [] yBoundary;
        delete [] zBoundary;

        int *cell2node = new int [nTotalCell * 4];
        int ncount = 0;
        for (int j = 1; j <= nj - 1; ++ j)
        {
            for (int i = 1; i <= ni - 1; ++ i)
            {
                cell2node[ncount++] = i     + (j - 1) * ni;
                cell2node[ncount++] = i + 1 + (j - 1) * ni;
                cell2node[ncount++] = i + 1 + (j   ) * ni;
                cell2node[ncount++] = i     + (j   ) * ni;
            }
        }

        PHWrite(data, cell2node, nTotalCell * 4);
        delete [] cell2node;
    }
}

void ComputeCharacteristicBoundaryUns(Grid *gridIn, DataContainer *data)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    grid->ComputeMetrics(NULL);

    if (gridIn->GetDim() == PHSPACE::TWO_D)
    {
        int tecioMission;

        int nTotalNode = gridIn->GetNTotalNode();
        int nTotalCell = gridIn->GetNTotalCell();
        int nBoundFace = gridIn->GetNBoundFace();
        int nTotalFace = gridIn->GetNTotalFace();

        int *face2Node       = grid->GetFace2Node();
        int *leftCellOfFace  = grid->GetLeftCellOfFace();
        int *rightCellOfFace = grid->GetRightCellOfFace();

        RDouble *x = gridIn->GetX();
        RDouble *y = gridIn->GetY();
        RDouble *z = gridIn->GetZ();

        vector < vector < int > > face2NodeList(nBoundFace);
        vector < int > linkmap;

        set< pair<int, string> > bcNameMap;
        UnstructBCSet **bcr = grid->GetBCRecord();
        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            int bcType = bcr[iFace]->GetKey();
            string bcName = bcr[iFace]->GetBCName();
            bcNameMap.insert(pair<int, string>(bcType, bcName));
        }

        set< pair<int, string> >::iterator iter;
        for (iter = bcNameMap.begin(); iter != bcNameMap.end(); ++ iter)
        {
            if ((*iter).first == PHENGLEI::INTERFACE || (*iter).first < 0) continue;

            tecioMission = WriteBoundary;
            PHWrite(data, &tecioMission, 1);

            int bcType = (*iter).first;
            PHWrite(data, &bcType, 1);

            string bcName = (*iter).second;
            data->WriteString(bcName);

            face2NodeList.resize(0);
            linkmap.resize(0);

            GetFace2NodeList(grid, *iter, linkmap, face2NodeList);

            int numPts = static_cast<int>(linkmap.size());
            int numElements = static_cast<int>(face2NodeList.size());

            PHWrite(data, &numPts, 1);
            PHWrite(data, &numElements, 1);

            RDouble *xBoundary = new RDouble [numPts];
            RDouble *yBoundary = new RDouble [numPts];
            RDouble *zBoundary = new RDouble [numPts];

            for (int iNode = 0; iNode < numPts; ++ iNode)
            {
                int it = linkmap[iNode];
                xBoundary[iNode] = x[it];
                yBoundary[iNode] = y[it];
                zBoundary[iNode] = z[it];
            }

            PHWrite(data, xBoundary, numPts);
            PHWrite(data, yBoundary, numPts);
            PHWrite(data, zBoundary, numPts);

            delete [] xBoundary;    xBoundary = nullptr;
            delete [] yBoundary;    yBoundary = nullptr;
            delete [] zBoundary;    zBoundary = nullptr;

            int *cell2Node = new int [numElements * 4];
            int nCount = 0;
            for (int iCell = 0; iCell < numElements; ++ iCell)
            {
                uint_t np = face2NodeList[iCell].size();
                int index0 = face2NodeList[iCell][0];
                for (uint_t ip = 0; ip < np; ++ ip)
                {
                    cell2Node[nCount] = face2NodeList[iCell][ip] + 1;
                    nCount ++;
                }

                for (uint_t ip = np; ip < 4; ++ ip)
                {
                    cell2Node[nCount] = index0 + 1;
                    nCount ++;
                }
            }

            PHWrite(data, cell2Node, numElements * 4);
            delete [] cell2Node;    cell2Node = nullptr;
        }

        HyList < int > *cell2Face = new HyList < int > [nTotalCell];
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            cell2Face[iCell].SetAverNnode(4);
        }
        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            int le = leftCellOfFace[iFace];
            cell2Face[le].insert(iFace);
        }
        for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
        {
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];
            cell2Face[le].insert(iFace);
            cell2Face[re].insert(iFace);
        }

        HyList < int > *cell2Node = new HyList < int > [nTotalCell];
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            cell2Node[iCell].SetAverNnode(4);
        }

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            int face0 = cell2Face[iCell][0];
            int ip1 = face2Node[2 * face0];
            int ip2 = face2Node[2 * face0 + 1];
            cell2Node[iCell].insert(ip1);
            cell2Node[iCell].insert(ip2);

            int targetNode = ip2;
            int checkNode  = ip1;

            int nFaceOfCell = cell2Face[iCell].size();
            bool isClockWise = false;
            while (!isClockWise)
            {
                for (int iFace = 1; iFace < nFaceOfCell; ++ iFace)
                {
                    int face  = cell2Face[iCell][iFace];
                    int node1 = face2Node[2 * face];
                    int node2 = face2Node[2 * face + 1];

                    if (node1 == targetNode)
                    {
                        if (node2 != checkNode)
                        {
                            cell2Node[iCell].insert(node2);
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
                        if (node1 != checkNode)
                        {
                            cell2Node[iCell].insert(node1);
                            targetNode = node1;
                        }
                        else
                        {
                            isClockWise = true;
                            break;
                        }
                    }
                }
            }
        }

        int nTotalCellVTK = 0;
        vector<RDouble> xVTK;
        vector<RDouble> yVTK;
        vector<RDouble> zVTK;
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            RDouble xNode = x[iNode];
            RDouble yNode = y[iNode];
            RDouble zNode = z[iNode];
            xVTK.push_back(xNode);
            yVTK.push_back(yNode);
            zVTK.push_back(zNode);
        }

        RDouble *xCC = grid->GetCellCenterX();
        RDouble *yCC = grid->GetCellCenterY();
        RDouble *zCC = grid->GetCellCenterZ();
        vector<int> cell2NodeVTK;

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            if (cell2Node[iCell].size() == 3)
            {
                for (int m = 0; m < 3; ++ m)
                {
                    cell2NodeVTK.push_back((cell2Node[iCell].GetData(m) + 1));
                }
                cell2NodeVTK.push_back((cell2Node[iCell].GetData(2) + 1));
                ++ nTotalCellVTK;
            }
            else if (cell2Node[ iCell ].size() == 4)
            {
                for (int m = 0; m < 4; ++ m)
                {
                    cell2NodeVTK.push_back((cell2Node[iCell].GetData(m) + 1));
                }
                ++ nTotalCellVTK;
            }
            else
            {
                int nDividedCell = cell2Node[iCell].size();
                int countNodeID  = 0;
                xVTK.push_back(xCC[iCell]);
                yVTK.push_back(yCC[iCell]);
                zVTK.push_back(zCC[iCell]);
                ++ nTotalNode;
                int cellCenterNodeID = nTotalNode - 1;

                for (int iDividedCell = 0; iDividedCell < nDividedCell; ++ iDividedCell)
                {
                    for (int iNode = 0; iNode < 2; ++ iNode)
                    {
                        cell2NodeVTK.push_back((cellCenterNodeID + 1));
                    }
                    cell2NodeVTK.push_back((cell2Node[iCell].GetData(iDividedCell) + 1));
                    if ((iDividedCell + 1) < nDividedCell)
                    {
                        cell2NodeVTK.push_back((cell2Node[iCell].GetData(iDividedCell + 1) + 1));
                    }
                    else
                    {
                        cell2NodeVTK.push_back((cell2Node[iCell].GetData(0) + 1));
                    }
                    ++ nTotalCellVTK;
                }
            }
        }

        tecioMission = WriteBlock;
        PHWrite(data, &tecioMission, 1);

        PHWrite(data, &nTotalNode, 1);
        PHWrite(data, &nTotalCellVTK, 1);

        PHWrite(data, xVTK, nTotalNode);
        PHWrite(data, yVTK, nTotalNode);
        PHWrite(data, zVTK, nTotalNode);

        PHWrite(data, cell2NodeVTK, nTotalCellVTK * 4);

        delete [] cell2Node;    cell2Node = nullptr;
        delete [] cell2Face;    cell2Face = nullptr;
    }
    else
    {
        int     nBoundFace = gridIn->GetNBoundFace();
        int     nTotalNode = gridIn->GetNTotalNode();
        RDouble *x         = gridIn->GetX();
        RDouble *y         = gridIn->GetY();
        RDouble *z         = gridIn->GetZ();
        RDouble *xfc       = grid->GetFaceCenterX();
        RDouble *yfc       = grid->GetFaceCenterY();
        RDouble *zfc       = grid->GetFaceCenterZ();

        using namespace PHENGLEI;
        vector < vector < int > > face2NodeList(nBoundFace);
        vector < int > linkmap;

        set< pair<int, string> > bcNameMap;

        UnstructBCSet **bcr = grid->GetBCRecord();
        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            int bcType = bcr[iFace]->GetKey();
            string bcName = bcr[iFace]->GetBCName();
            bcNameMap.insert(pair<int, string>(bcType, bcName));
        }

        int tecioMission;
        set< pair<int, string> >::iterator iter;
        for (iter = bcNameMap.begin(); iter != bcNameMap.end(); ++ iter)
        {
            if ((*iter).first == PHENGLEI::INTERFACE || (*iter).first < 0) continue;

            tecioMission = WriteBoundary;
            PHWrite(data, &tecioMission, 1);

            int bcType = (*iter).first;
            PHWrite(data, &bcType, 1);

            string bcName = (*iter).second;
            data->WriteString(bcName);

            face2NodeList.resize(0);
            linkmap.resize(0);

            GetFace2NodeList(grid, *iter, linkmap, face2NodeList);

            int numPts = static_cast<int>(linkmap.size());
            int numElements = static_cast<int>(face2NodeList.size());

            PHWrite(data, &numPts, 1);
            PHWrite(data, &numElements, 1);

            RDouble *xBoundary = new RDouble [numPts];
            RDouble *yBoundary = new RDouble [numPts];
            RDouble *zBoundary = new RDouble [numPts];

            for (int iNode = 0; iNode < numPts; ++ iNode)
            {
                int nodeIndex = linkmap[iNode];
                if (nodeIndex < nTotalNode)
                {
                    xBoundary[iNode] = x[nodeIndex];
                    yBoundary[iNode] = y[nodeIndex];
                    zBoundary[iNode] = z[nodeIndex];
                }
                else
                {
                    int faceIndex = nodeIndex - nTotalNode;
                    xBoundary[iNode] = xfc[faceIndex];
                    yBoundary[iNode] = yfc[faceIndex];
                    zBoundary[iNode] = zfc[faceIndex];
                }
            }

            PHWrite(data, xBoundary, numPts);
            PHWrite(data, yBoundary, numPts);
            PHWrite(data, zBoundary, numPts);

            delete [] xBoundary;    xBoundary = nullptr;
            delete [] yBoundary;    yBoundary = nullptr;
            delete [] zBoundary;    zBoundary = nullptr;

            int *cell2Node = new int [numElements * 4];
            int nCount = 0;
            for (int iCell = 0; iCell < numElements; ++ iCell)
            {
                uint_t np     = face2NodeList[iCell].size();
                int    index0 = face2NodeList[iCell][0];
                for (uint_t ip = 0; ip < np; ++ ip)
                {
                    cell2Node[nCount] = face2NodeList[iCell][ip] + 1;
                    nCount ++;
                }

                for (uint_t ip = np; ip < 4; ++ ip)
                {
                    cell2Node[nCount] = index0 + 1;
                    nCount ++;
                }
            }

            PHWrite(data, cell2Node, numElements * 4);
            delete [] cell2Node;    cell2Node = nullptr;
        }
    }
}

}
