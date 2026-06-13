#include "Pre_HDF5File.h"
#include "Glb_Dimension.h"
#include "GridFactory.h"
#include "GridgenIn.h"
#include "IO_FileName.h"
#include "TK_Parse.h"
#include "FluentIn.h"
#include "OverLappingGrid.h"
#include "PHIO.h"
#include "TK_Log.h"
#include "Pre_StructToUnstruct.h"
#include "AleManager.h"


#include "FieldViewIn.h"
#include "GridView.h"
#include "Pre_CGNSToPlot3D_Struct.h"
#include "Mesh_Refine_Struct.h"
#include "Pre_CGNSToPHengLEI_Unstruct.h"
#include "Pre_FluentToPHengLEI_Unstruct.h"
#include "Pre_FieldViewToPHengLEI_Unstruct.h"
#include "Pre_Plot3DToPHengLEI_Struct.h"
#include "Pre_CGNSToPHengLEI_Struct.h"
#include "Gmesh.h"
#include "Pre_GridTranslator.h"
#include "Pre_GridCombine.h"
#include "Pre_GridMirror.h"
#include "Mesh_Deformation.h"
#include "Mesh_Refine_Unstruct.h"
#pragma warning (disable:913)

namespace PHSPACE
{
void CreateGrid()
{
    int gridObject = GlobalDataBase::GetIntParaFromDB("gridobj");

    switch (gridObject)
    {
    case PHSPACE::GRID_CONVERSION:
        ConvertGrid();
        break;
    case PHSPACE::GRID_REFINE:
        RefineGrid();
        break;
    case PHSPACE::GRID_MERGING:
        CombineGrid();
        break;
    case PHSPACE::GRID_DEFORMATION:
        DeformGrid();
        break;
    case PHSPACE::GRID_MIRROR:
        MirrorGrid();
        break;
    case PHSPACE::GRID_STRTOUNS:
        StructToUnstructGrid();
    default:
        break;
    }
}

void OversetGridView()
{
    BlockGroupManager *blockGroupManager = new BlockGroupManager();
    blockGroupManager->Run();
    delete blockGroupManager;

    return;
}

void ReadGMSH(CGNSFactory *&factoryCGNS, int &numberOfBlocks, const string &filename)
{
    fstream file;
    file.open(filename.c_str(), ios_base::in);
    if (!file)
    {
        TK_Exit::FileOpenErrorExit(filename);
    }

    numberOfBlocks = 1;

    factoryCGNS = new CGNSFactory(numberOfBlocks);

    Gmesh *gmsh = new Gmesh();
    gmsh->Read(file, factoryCGNS);
    delete gmsh;

    file.close();
    file.clear();
}

void GridgenToPHengLEI()
{
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");
    string outGridFile  = GlobalDataBase::GetStrParaFromDB("out_gfile");
    string boundaryFile = GlobalDataBase::GetStrParaFromDB("bnd_file");
    int dumpOldGrid = GlobalDataBase::GetIntParaFromDB("dumpOldGrid");

    int codeOfOversetGrid = GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");

    int numberOfBlocks;
    Grid **structgrids;
    ReadGridgenGrid(structgrids, numberOfBlocks, fromGridFile, boundaryFile);

    int ndim = PHSPACE::GetDim();
    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        structgrids[iZone]->SetLnkInfo();
    }

    if (ndim == TWO_D || codeOfOversetGrid == 10 || (ndim == THREE_D && periodicType != NO_PERIODICITY))
    {
        CheckMeshMultigrid(structgrids, numberOfBlocks);

        for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
        {
            structgrids[iZone]->ProcessBCInfo();
        }

        ConstructGlobalInterfaceLink(structgrids, numberOfBlocks);

        if (dumpOldGrid)
        {
            WriteLnkFile(outGridFile, structgrids, numberOfBlocks);
        }

        DumpGrid(outGridFile, structgrids, numberOfBlocks);
        WriteBoundaryInfo(structgrids, numberOfBlocks);
        DumpCharacteristicBoundary(outGridFile, structgrids, numberOfBlocks);
        WriteBCNameInfo(outGridFile, structgrids, numberOfBlocks);

        for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
        {
            delete structgrids[iZone];
        }
        delete [] structgrids;
    }
    else
    {
        if (dumpOldGrid)
        {
            WriteLnkFile(outGridFile, structgrids, numberOfBlocks);
        }

        CheckMeshMultigrid(structgrids, numberOfBlocks);
        GridTranslator *gridTranslator = new GridTranslator(structgrids, numberOfBlocks);
        gridTranslator->Run();

        WriteBoundaryInfo(structgrids, numberOfBlocks);
        DumpCharacteristicBoundary(outGridFile, structgrids, numberOfBlocks);
        WriteBCNameInfo(outGridFile, structgrids, numberOfBlocks);

        delete gridTranslator;
    }
}

void GMSH2PHengLEI()
{
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");
    string outGridFile = GlobalDataBase::GetStrParaFromDB("out_gfile");

    cout << "    Grid file name: " << fromGridFile << "\n";

    int numberOfBlocks;
    CGNSFactory *factoryCGNS;
    ReadGMSH(factoryCGNS, numberOfBlocks, fromGridFile);

    factoryCGNS->ConvertGrid2Fantasy();

    Grid **grids = new Grid * [numberOfBlocks];

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        GridID *index = new GridID(iZone);
        Grid *grid = CreateGridGeneral(UNSTRUCTGRID, index, 0, GetDim());
        grids[iZone] = grid;
        CGNS2UnsGrid(factoryCGNS, iZone, UnstructGridCast(grid));
        grid->ComputeMinMaxBox();
    }

    delete factoryCGNS;

    ConstructGlobalInterfaceLink(grids, numberOfBlocks);
    DumpGrid(outGridFile, grids, numberOfBlocks);

    WriteBoundaryInfo(grids, numberOfBlocks);
    WriteBCNameInfo(outGridFile, grids, numberOfBlocks);
    DumpCharacteristicBoundary(outGridFile, grids, numberOfBlocks);

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        delete grids[iZone];
    }
    delete [] grids;
}

void Ustar2PHengLEI()
{
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");
    string outGridFile = GlobalDataBase::GetStrParaFromDB("out_gfile");

    int numberOfBlocks = 1;
    Grid **grids = new Grid * [numberOfBlocks];
    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        GridID *index = new GridID(iZone);
        Grid *grid = CreateGridGeneral(UNSTRUCTGRID, index, 0, GetDim());
        grids[iZone] = grid;
        ReadUstarGrid(fromGridFile, UnstructGridCast(grid));
        grid->ComputeMinMaxBox();
    }

    ConstructGlobalInterfaceLink(grids, numberOfBlocks);
    DumpGrid(outGridFile, grids, numberOfBlocks);

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        delete grids[iZone];
    }
    delete [] grids;
}

void Plot3DToPHengLEI()
{
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");
    string outGridFile  = GlobalDataBase::GetStrParaFromDB("out_gfile");
    string boundaryFile = GlobalDataBase::GetStrParaFromDB("bnd_file");
    int dumpOldGrid = GlobalDataBase::GetIntParaFromDB("dumpOldGrid");

    int isOverset = GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");

    int numberOfBlocks;
    Grid **structgrids;
    ReadPlot3dGrid(structgrids, numberOfBlocks, fromGridFile, boundaryFile);

    int ndim = PHSPACE::GetDim();
    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        structgrids[iZone]->SetLnkInfo();
    }

    if (ndim == TWO_D || isOverset || (ndim == THREE_D && periodicType != NO_PERIODICITY))
    {
        CheckMeshMultigrid(structgrids, numberOfBlocks);

        for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
        {
            structgrids[iZone]->ProcessBCInfo();
        }

        ConstructGlobalInterfaceLink(structgrids, numberOfBlocks);

        if (dumpOldGrid)
        {
            WriteLnkFile(outGridFile, structgrids, numberOfBlocks);
        }

        DumpGrid(outGridFile, structgrids, numberOfBlocks);
        WriteBoundaryInfo(structgrids, numberOfBlocks);
        DumpCharacteristicBoundary(outGridFile, structgrids, numberOfBlocks);
        WriteBCNameInfo(outGridFile, structgrids, numberOfBlocks);

        for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
        {
            delete structgrids[iZone];
        }
        delete [] structgrids;
    }
    else
    {
        if (dumpOldGrid)
        {
            WriteLnkFile(outGridFile, structgrids, numberOfBlocks);
        }

        CheckMeshMultigrid(structgrids, numberOfBlocks);
        GridTranslator *gridTranslator = new GridTranslator(structgrids, numberOfBlocks);
        gridTranslator->Run();

        WriteBoundaryInfo(structgrids, numberOfBlocks);
        DumpCharacteristicBoundary(outGridFile, structgrids, numberOfBlocks);
        WriteBCNameInfo(outGridFile, structgrids, numberOfBlocks);

        delete gridTranslator;
    }
}

void Plot3DGridProcess()
{
    int toGridType = GlobalDataBase::GetIntParaFromDB("to_gtype");

    if (toGridType == PHENGLEI_TYPE)
    {
        Plot3DToPHengLEI();
    }
}

void CGNS2StructGrid()
{
    string cgnsfile = GlobalDataBase::GetStrParaFromDB("from_gfile");

    Pre_CGNSToPlot3D_Struct cgnsGridConversion(cgnsfile);

    cgnsGridConversion.Run();

    Plot3DGridProcess();
}

void CGNS2UnstructGrid()
{
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");
    string outGridFile = GlobalDataBase::GetStrParaFromDB("out_gfile");

    //! Read CGNS grid and convert it's format.
    Pre_CGNSToPHengLEI_Unstruct cgnsGridConversion(fromGridFile);
    cgnsGridConversion.Run();

    int numberOfBlocks = cgnsGridConversion.GetNumberofBlocks();
    Grid **grids = cgnsGridConversion.GetGrids();

    //! Construct block connection and dump out to file.
    ConstructGlobalInterfaceLink(grids, numberOfBlocks);

    int referenceFrame = PHSPACE::GlobalDataBase::GetIntParaFromDB("referenceFrame");
    int nTurboZone = 1;
    if (referenceFrame == ROTATIONAL_FRAME)
    {
        nTurboZone = PHSPACE::GlobalDataBase::GetIntParaFromDB("nTurboZone");
    }

    int hasVolumeCondition = GlobalDataBase::GetIntParaFromDB("hasVolumeCondition");
    if ((numberOfBlocks > 1) && (!hasVolumeCondition) && nTurboZone<=1)
    {
        cout << "\n";
        cout << "Blocks > 1, Start to Combine Grids" << endl;
        CombinGrid *combinedGrid = new CombinGrid(0, grids, numberOfBlocks);

        combinedGrid->RunCombinGrid();

        Grid **gridMerged = combinedGrid->GetGridAll();

        for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
        {
            delete grids[iZone];
        }

        grids = gridMerged;

        numberOfBlocks = 1;
        cgnsGridConversion.SetNumberofBlocks(numberOfBlocks);

        DumpGrid(outGridFile, gridMerged, 1);

        delete combinedGrid;
    }
    else
    {
        DumpGrid(outGridFile, grids, numberOfBlocks);
    }

    cgnsGridConversion.WriteAdditionalInformation(outGridFile, grids);
    DumpCharacteristicBoundary(outGridFile, grids, numberOfBlocks);

    delete grids[0];
}

void Fantasy2Plot3D()
{
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");
    GlobalDataBase::UpdateData("gridfile", &fromGridFile, PHSTRING, 1);

    string outGridFile = GlobalDataBase::GetStrParaFromDB("out_gfile");

    string fname, fext;
    GetNameExt(outGridFile, fname, fext, ".");

    string cellFile = fname + ".grd";
    string boundaryFile = fname + ".inp";

    using namespace PHMPI;

    Region *region = new Region();

    region->ReadGrid();

    int nZones = PHMPI::GetNumberofGlobalZones();
    Grid **Grids = new Grid * [nZones];

    for (int iZone = 0; iZone < nZones; iZone ++)
    {
        Grids[iZone] = GetGrid(iZone, 0);
    }

    ReProcessBCInfo(Grids, nZones);

    int isOldGrid = GlobalDataBase::GetIntParaFromDB("isOldGrid");
    if (isOldGrid)
    {
        ReadLnkFile(fromGridFile, Grids, nZones);
    }
    else
    {
        ReadLnkInfor(fromGridFile, Grids, nZones);
    }

    cout << endl << "Writing grd file ...." << endl;
    FantasyWriteToGrd(cellFile, Grids, nZones);
    FantasyWriteToInp(boundaryFile, Grids, nZones);

    delete region;
}

void FantasyOld2New()
{
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");
    GlobalDataBase::UpdateData("gridfile", &fromGridFile, PHSTRING, 1);

    string outGridFile = GlobalDataBase::GetStrParaFromDB("out_gfile");

    Region *region = new Region();
    region->ReadGrid();

    int nZones = PHMPI::GetNumberofGlobalZones();
    Grid **Grids = new Grid * [nZones];

    for (int iZone = 0; iZone < nZones; iZone ++)
    {
        Grids[iZone] = GetGrid(iZone);
    }

    DumpGrid(outGridFile, Grids, nZones);

    delete region;
}

Pre_GridConversion * GetConvertType(const string &gridFileName, int gridType, int fromGridType)
{
    Pre_GridConversion *gridConversion = NULL;
    if (gridType == UNSTRUCTGRID)
    {
        switch (fromGridType)
        {
            case CGNS_TYPE:
                gridConversion = new Pre_CGNSToPHengLEI_Unstruct(gridFileName);
                break;
            case FLUENT_TYPE:
                gridConversion = new Pre_FluentToPHengLEI_Unstruct(gridFileName);
                break;
            case FIELDVIEW_TYPE:
                gridConversion = new Pre_FieldViewToPHengLEI_Unstruct(gridFileName);
                break;
            default:
                TK_Exit::UnexpectedVarValue("from_gtype", fromGridType);
                break;
        }
    }
    else if (gridType == STRUCTGRID)
    {
        switch (fromGridType)
        {
            case CGNS_TYPE:
                gridConversion = new Pre_CGNSToPHengLEI_Struct(gridFileName);
                break;
            case PLOT3D_TYPE:
                gridConversion = new Pre_Plot3DToPHengLEI_Struct(gridFileName);
                break;
            default:
                TK_Exit::UnexpectedVarValue("from_gtype", fromGridType);
                break;
        }
    }
    else
    {
        TK_Exit::UnexpectedVarValue("gridtype", gridType);
    }

    return gridConversion;
}

void GridgenGridProcess()
{
    int toGridType = GlobalDataBase::GetIntParaFromDB("to_gtype");

    if (toGridType == PHENGLEI_TYPE)
    {
        GridgenToPHengLEI();
    }
}

void GmeshGridProcess()
{
    int toGridType = GlobalDataBase::GetIntParaFromDB("to_gtype");

    if (toGridType == PHENGLEI_TYPE)
    {
        GMSH2PHengLEI();
    }
}

void UstarGridProcess()
{
    int toGridType = GlobalDataBase::GetIntParaFromDB("to_gtype");

    if (toGridType == PHENGLEI_TYPE)
    {
        Ustar2PHengLEI();
    }
}

void FluentGridProcess()
{
    int toGridType = GlobalDataBase::GetIntParaFromDB("to_gtype");

    if (toGridType == PHENGLEI_TYPE)
    {
        Fluent2PHengLEI();
    }
}

void FieldViewGridProcess()
{
    int toGridType = GlobalDataBase::GetIntParaFromDB("to_gtype");

    if (toGridType == PHENGLEI_TYPE)
    {
        FieldView2PHengLEI();
    }
}

void CGNSGridProcess()
{
    int gridtype = GlobalDataBase::GetIntParaFromDB("gridtype");

    if (gridtype == STRUCTGRID)
    {
        CGNS2StructGrid();
    }
    else if (gridtype == UNSTRUCTGRID)
    {
        CGNS2UnstructGrid();
    }
}

void Fantasy2Fluent()
{
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");
    GlobalDataBase::UpdateData("gridfile", &fromGridFile, PHSTRING, 1);

    string outGridFile = GlobalDataBase::GetStrParaFromDB("out_gfile");

    using namespace PHMPI;

    Region *region = new Region();
    region->ReadGrid();

    int nZones = PHMPI::GetNumberofGlobalZones();
    Grid **grids = new Grid * [nZones];

    for (int iZone = 0; iZone < nZones; iZone ++)
    {
        grids[iZone] = GetGrid(iZone, 0);
    }

    if (nZones > 1)
    {
        cout << "\n";
        cout << "Blocks > 1, Start to Combine Grids" << endl;
        CombinGrid *combinedGrid = new CombinGrid(0, grids, nZones);
        combinedGrid->RunCombinGrid();

        Grid **gridMerged = combinedGrid->GetGridAll();

        FantasyWriteToFluent(outGridFile, gridMerged);

        delete gridMerged[0];
        delete [] gridMerged;
        delete combinedGrid;
    }
    else
    {
        FantasyWriteToFluent(outGridFile, grids);
    }

    delete [] grids;
    delete region;
}

void FantasyGridProcess()
{
    int toGridType = GlobalDataBase::GetIntParaFromDB("to_gtype");

    if (toGridType == USTAR_TYPE)
    {
        Fantasy2Ustar();
    }
    else if (toGridType == PLOT3D_TYPE)
    {
        Fantasy2Plot3D();
    }
    else if (toGridType == PHENGLEI_TYPE)
    {
        FantasyOld2New();
    }
    else if (toGridType == FLUENT_TYPE)
    {
        Fantasy2Fluent();
    }

}

void MixedGridProcess()
{
    int fromGridType = MIXGRID_TYPE;
    int omitBlankBC = 0;
    GlobalDataBase::UpdateData("from_gtype", &fromGridType, PHINT, 1);
    GlobalDataBase::UpdateData("omit_no_bound_bc", &omitBlankBC, PHINT, 1);

    int numberOfGridGroups = 2;
    string structGridInMixture = GlobalDataBase::GetStrParaFromDB("mixgrid_str");
    string unstructGridInMixture = GlobalDataBase::GetStrParaFromDB("mixgrid_uns");

    GlobalDataBase::UpdateData("gridfile", &structGridInMixture, PHSTRING, 1);
    GlobalDataBase::UpdateData("gridfile1", &unstructGridInMixture, PHSTRING, 1);
    GlobalDataBase::UpdateData("numberOfGridGroups", &numberOfGridGroups, PHINT, 1);

    Region *region = new Region();
    region->ReadGrid();
    region->BuildAndDumpMixedGrid();
    delete region;
}

void MultiGridConvert()
{
    int toGridType = GlobalDataBase::GetIntParaFromDB("to_gtype");
    if (toGridType != PHENGLEI_TYPE)
    {
        ostringstream oss;
        oss << "  Error: this situation has not been considered, for to_gtype" << " = " << toGridType << " when numberOfGridFile > 1" << endl;
        TK_Exit::ExceptionExit(oss.str());
    }

    int fromGridType = GlobalDataBase::GetIntParaFromDB("from_gtype");
    int gridType = GlobalDataBase::GetIntParaFromDB("gridtype");

    int gridTypeMultiGrid = GlobalDataBase::GetIntParaFromDB("gridtype");
    int numberOfGridFile = PHSPACE::GlobalDataBase::GetIntParaFromDB("numberOfGridFile");
    PHVectorString1D gridFileNameList;
    gridFileNameList.push_back(PHSPACE::GlobalDataBase::GetStrParaFromDB("from_gfile"));

    for (int iGridFile = 1; iGridFile < numberOfGridFile; ++ iGridFile)
    {
        ostringstream oss;
        oss << "from_gfile" << iGridFile;
        gridFileNameList.push_back(PHSPACE::GlobalDataBase::GetStrParaFromDB(oss.str()));
    }

    vector <Grid *> totalGrids;
    totalGrids.resize(0);

    for (int iGridFile = 0; iGridFile < numberOfGridFile; ++ iGridFile)
    {
        string gridFileName = gridFileNameList[iGridFile];
        if (gridTypeMultiGrid > 1)
        {
            fromGridType = PHSPACE::GetIntegerParameterFromDataBase(iGridFile,"from_gtype");
            gridType = PHSPACE::GetIntegerParameterFromDataBase(iGridFile, "gridtype");
        }

        Pre_GridConversion *gridConversion = GetConvertType(gridFileName, gridType, fromGridType);
        gridConversion->Run();

        int nZones = gridConversion->GetNumberofBlocks();
        Grid **grids = gridConversion->GetGrids();

        ConstructGlobalInterfaceLink(grids, nZones);

        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            grids[iZone]->SetFileIndex(iGridFile);
            totalGrids.push_back(grids[iZone]);
        }

        string outGridFile = ChangeExtensionOfFileName(gridFileName, "fts");
        DumpGrid(outGridFile, grids, nZones);
        DumpCharacteristicBoundary(outGridFile, grids, nZones);
        WriteBCNameInfo(outGridFile, grids, nZones);
        gridConversion->WriteVolumeName(outGridFile, grids);

        delete gridConversion;
        delete [] grids;
    }

    WriteBoundaryInfo(&totalGrids[0], static_cast<int>(totalGrids.size()));

    WriteVolumeInfo(&totalGrids[0], static_cast<int>(totalGrids.size()));

    for (int iZone = 0; iZone < totalGrids.size(); ++ iZone)
    {
        delete totalGrids[iZone];
    }
}

void ConvertGrid()
{
    int numberOfGridFile = GlobalDataBase::GetIntParaFromDB("numberOfGridFile");
    if (numberOfGridFile > 1)
    {
        MultiGridConvert();
        return;
    }

    int fromGridType = GlobalDataBase::GetIntParaFromDB("from_gtype");
    int gridType     = GlobalDataBase::GetIntParaFromDB("gridtype");

    if (gridType == MIXGRID)
    {
        MixedGridProcess();
    }
    else
    {
        switch (fromGridType)
        {
            case PHENGLEI_TYPE:
                FantasyGridProcess();
                break;
            case CGNS_TYPE:
                CGNSGridProcess();
                break;
            case PLOT3D_TYPE:
                Plot3DGridProcess();
                break;
            case FIELDVIEW_TYPE:
                FieldViewGridProcess();
                break;
            case FLUENT_TYPE:
                FluentGridProcess();
                break;
            case USTAR_TYPE:
                UstarGridProcess();
                break;
            case GMSH_TYPE:
                GmeshGridProcess();
                break;
            case GRIDGEN_TYPE:
                GridgenGridProcess();
                break;
            default:
                TK_Exit::UnexpectedVarValue("from_gtype", fromGridType);
                break;
        }
    }
}

void RefineGrid()
{
    RefineParameter *refineParameter = new RefineParameter();
    int anisoRefine = GlobalDataBase::GetIntParaFromDB("anisoRefine");

    refineParameter->SetAnisoRefineType(anisoRefine);

    int gridtype = GlobalDataBase::GetIntParaFromDB("gridtype");
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");
    Mesh_Refine *meshRefine = 0;
    if (gridtype == UNSTRUCTGRID)
    {
        meshRefine = new Mesh_Refine_Unstruct(fromGridFile);
    }
    else if (gridtype == STRUCTGRID)
    {
        meshRefine = new Mesh_Refine_Struct(fromGridFile);

    }
    meshRefine->SetRefineParameter(refineParameter);
    meshRefine->Run();

    delete refineParameter;
    delete meshRefine;
}

void CombineGrid()
{
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");
    string outGridFile = GlobalDataBase::GetStrParaFromDB("out_gfile");
    GlobalDataBase::UpdateData("gridfile", &fromGridFile, PHSTRING, 1);

    Region *region = new Region();
    region->ReadGrid();
    region->SetupMultiBlockRelationship();

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    if (nZones < 2)
    {
        TK_Exit::ExceptionExit("/n Error: number of blocks < 2, don't need to merge or combine grid!");
    }

    int nLocalZones = GetNumberofLocalZones();
    Grid **grids = new Grid * [nLocalZones];
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        grids[iZone] = GetGrid(zoneID);
    }

    CombinGrid *combinedGrid = new CombinGrid(region, grids, nLocalZones);

    combinedGrid->RunCombinGrid();

    DumpGrid(outGridFile, combinedGrid->GetGridAll(), 1);

    delete region;    region = nullptr;
    delete combinedGrid;    combinedGrid = nullptr;
}

void RepairGrid()
{
}

void MirrorGrid()
{
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");

    Pre_GridMirror GridMirror(fromGridFile);

    GridMirror.Run();
}

void StructToUnstructGrid()
{
    //! Read grid.
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");
    GlobalDataBase::UpdateData("gridfile", &fromGridFile, PHSTRING, 1);
    Region *region = new Region();
    region->ReadGrid();

    //! Get struct grid.
    int nZones = PHMPI::GetNumberofGlobalZones();
    Grid **structGrid = new Grid * [nZones];
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        structGrid[iZone] = GetGrid(iZone);
    }

    //! Struct to unstruct grid.
    Pre_StructToUnstruct *struct2Unstruct = new Pre_StructToUnstruct(structGrid, nZones);
    struct2Unstruct->Run();

    //! Dump grid.
    Grid **unstructGrid = struct2Unstruct->GetUnstructGrids();
    string outGridFile = GlobalDataBase::GetStrParaFromDB("out_gfile");
    DumpGrid(outGridFile, unstructGrid, nZones);

    delete [] structGrid;
    delete struct2Unstruct;
    delete region;
}

void GridgenToMIX()
{
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");
    string outGridFile  = GlobalDataBase::GetStrParaFromDB("out_gfile");
    string boundaryFile = GlobalDataBase::GetStrParaFromDB("bnd_file");

    int numberOfBlocks;
    Grid **structgrids;
    ReadPlot3dGrid(structgrids, numberOfBlocks, fromGridFile, boundaryFile);

    Grid **grids = new Grid * [numberOfBlocks];

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        GridID *index = new GridID(iZone);
        Grid *grid;
        if (iZone % 2 != 0)
        {
            grid = CreateGridGeneral(UNSTRUCTGRID, index, 0, GetDim());
            grids[iZone] = grid;
            StrGridToUnsGrid(StructGridCast(structgrids[iZone]), UnstructGridCast(grid));
        }
        else
        {
            grid = CreateGridGeneral(STRUCTGRID, index, 0, GetDim());
            grids[iZone] = grid;
            StructGridCast(grid)->CopyGrid(StructGridCast(structgrids[iZone]));
            grid->ProcessBCInfo();
        }
    }

    ConstructGlobalInterfaceLink(grids, numberOfBlocks);
    DumpGrid(outGridFile, grids, numberOfBlocks);

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        delete grids[iZone];
    }
    delete [] grids;

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        delete structgrids[iZone];
    }
    delete [] structgrids;
}

void GridgenToMBlockUns()
{
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");
    string outGridFile  = GlobalDataBase::GetStrParaFromDB("out_gfile");
    string boundaryFile = GlobalDataBase::GetStrParaFromDB("bnd_file");

    int numberOfBlocks;

    Grid **structgrids;
    ReadPlot3dGrid(structgrids, numberOfBlocks, fromGridFile, boundaryFile);

    Grid **grids = new Grid * [numberOfBlocks];

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        GridID *index = new GridID(iZone);
        Grid *grid = CreateGridGeneral(UNSTRUCTGRID, index, 0, GetDim());
        grids[iZone] = grid;
        StrGridToUnsGrid(StructGridCast(structgrids[iZone]), UnstructGridCast(grid));
    }

    ConstructGlobalInterfaceLink(grids, numberOfBlocks);
    DumpGrid(outGridFile, grids, numberOfBlocks);

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        delete grids[iZone];
    }
    delete [] grids;

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        delete structgrids[iZone];
    }
    delete [] structgrids;
}

void GridgenToSingleUns()
{
    string fromGridFile = GlobalDataBase::GetStrParaFromDB("from_gfile");
    string outGridFile  = GlobalDataBase::GetStrParaFromDB("out_gfile");
    string boundaryFile = GlobalDataBase::GetStrParaFromDB("bnd_file");

    int nBlocks;
    Grid **structgrids;
    ReadPlot3dGrid(structgrids, nBlocks, fromGridFile, boundaryFile);

    int numberOfBlocks = 1;

    Grid **grids = new Grid * [numberOfBlocks];

    int nb = 0;
    GridID *index = new GridID(nb);
    Grid *grid = CreateGridGeneral(UNSTRUCTGRID, index, 0, GetDim());
    grids[nb] = grid;

    StrGridToUnsGrid(structgrids, nBlocks, UnstructGridCast(grid));

    ConstructGlobalInterfaceLink(grids, numberOfBlocks);
    DumpGrid(outGridFile, grids, numberOfBlocks);

    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        delete grids[iZone];
    }
    delete [] grids;

    for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
    {
        delete structgrids[iBlock];
    }
    delete [] structgrids;
}

void Gridgen2Fantasy()
{
    int gridtype   = GlobalDataBase::GetIntParaFromDB("gridtype");
    int multiblock = GlobalDataBase::GetIntParaFromDB("multiblock");

    int sys_gridtype;
    if (gridtype == 0)
    {
        sys_gridtype = UNSTRUCTGRID;
        GlobalDataBase::UpdateData("sys_gridtype", &sys_gridtype, PHINT, 1);

        if (multiblock)
        {
            GridgenToMBlockUns();
        }
        else
        {
            GridgenToSingleUns();
        }
    }
    else if (gridtype == 1)
    {
        sys_gridtype = STRUCTGRID;
        GlobalDataBase::UpdateData("sys_gridtype", &sys_gridtype, PHINT, 1);

        Plot3DToPHengLEI();
    }
}

}

