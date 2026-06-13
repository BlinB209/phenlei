#include <cmath>
#include "Pre_GridTranslator.h"
#include "Pre_GridConversion.h"
#include "MixGrid.h"
#include "Task_WriteBinaryTecplot.h"
#include "Geo_StructBC.h"
#include "MultigridFactory.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "TK_Exit.h"

using namespace std;

namespace PHSPACE
{

GridTranslator::GridTranslator()
{
    GlobalDataBase::GetData("codeOfOversetGrid", &isOverset, PHINT, 1);

    gridGroup = new vector<StructGrid *>();

    structFaceGroup = new vector<StructFace *>();
}

GridTranslator::GridTranslator(Grid **&strgrids, int &nBlocks)
{
    numberOfZones = nBlocks;

    GlobalDataBase::GetData("codeOfOversetGrid", &isOverset, PHINT, 1);

    gridGroup = new vector<StructGrid *>();
    gridGroup->resize(nBlocks, NULL);

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        (*gridGroup)[iZone] = StructGridCast(strgrids[iZone]);
    }

    structFaceGroup = new vector<StructFace *>();
}

GridTranslator::~GridTranslator()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        delete (*structFaceGroup)[iZone];
    }
    delete structFaceGroup;    structFaceGroup = NULL;

    delete gridGroup;
    gridGroup = NULL;
}

void GridTranslator::Run()
{
    ProcessBoundaryCondition();

    GenerateStructFaceGroup();

    ColorateStructFaceGroup();

    CollectInterfaceInformation();

    OutputGrid();
}

void GridTranslator::ReadGridgenFiles()
{
}

void GridTranslator::ReadGridgenCoordinateFile()
{
}

void GridTranslator::ReadGridgenBoundaryFile()
{
}

void GridTranslator::ProcessBoundaryCondition()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        StructBCSet *structBCSet = (*gridGroup)[iZone]->GetStructBCSet();
        
        structBCSet->ProcessBCInfo();
    }
}

void GridTranslator::GenerateStructFaceGroup()
{
    structFaceGroup->resize(numberOfZones, NULL);

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        (*structFaceGroup)[iZone] = new StructFace((*gridGroup)[iZone]);
    }
}

void GridTranslator::ColorateStructFaceGroup()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        
        StructFace *structFace = (*structFaceGroup)[iZone];
        structFace->SetStructFaceGroup(this->structFaceGroup);
        structFace->PaintSelf();
    }

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        StructFace *structFace = (*structFaceGroup)[iZone];
        structFace->PaintNeighbor();
    }
}

void GridTranslator::CollectInterfaceInformation()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        cout << "iZone = " << iZone << endl;

        int numberOfInterfaces = (*structFaceGroup)[iZone]->ComputeNumberOfInterfaces();

        cout << " numberOfInterfaces = " << numberOfInterfaces << endl;

        StructGrid *structuredGrid = (*gridGroup)[iZone];

        structuredGrid->SetNIFace(numberOfInterfaces);

        if (numberOfInterfaces == 0)
        {
            structuredGrid->SetInterfaceInfo(NULL);
        }
        else
        {
            structuredGrid->SetInterfaceInfo(new InterfaceInfo(numberOfInterfaces));
        }

        InterfaceInfo *interfaceInfo = structuredGrid->GetInterfaceInfo();

        (*structFaceGroup)[iZone]->CollectInterfaceInformation(interfaceInfo);
    }
}

void GridTranslator::OutputGrid()
{
    string out_gfile = "grid.fts";
    int tasktype = GetTaskCode();
    int numberOfMultifile = PHSPACE::GlobalDataBase::GetIntParaFromDB("numberOfMultifile");
    if (tasktype == PHSPACE::PARTITION_GRID)
    {
        GlobalDataBase::GetData("partition_grid_file", &out_gfile, PHSTRING, 1);
    }
    else if (tasktype == PHSPACE::CREATE_GRID)
    {
        GlobalDataBase::GetData("out_gfile", &out_gfile, PHSTRING, 1);
    }

    Grid **grids = new Grid *[numberOfZones];
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        grids[iZone] = (*gridGroup)[iZone];
    }

    if (tasktype == PHSPACE::PARTITION_GRID && numberOfMultifile > 1)
    {
        int maximumProcess = PHSPACE::GlobalDataBase::GetIntParaFromDB("maxproc");
        
        if (maximumProcess < numberOfMultifile)
        {
            ostringstream oss;
            oss << " Error: the numberOfMultifile's value is larger than maximumProcess's value " << endl;
            TK_Exit::ExceptionExit(oss.str());
        }
        else if (maximumProcess % numberOfMultifile != 0)
        {
            ostringstream oss;
            oss << " Error: the maximumProcess's value is not an integral multiple of numberOfMultifile's value " << endl;
            TK_Exit::ExceptionExit(oss.str());
        }

        int procNumberOfEachfile = maximumProcess / numberOfMultifile;
        for (int iFile = 0; iFile < numberOfMultifile; ++iFile)
        {
            vector<Grid *> tempGrids;
            tempGrids.resize(0);

            for (int iZone = 0; iZone < numberOfZones; ++iZone)
            {
                int procIndex = grids[iZone]->GetGridID()->GetIndex();
                int *block_proc_inp = PHMPI::GetZoneProcessorID_INP();
                if (block_proc_inp)
                {
                    procIndex = block_proc_inp[iZone];
                }

                if ((procIndex / procNumberOfEachfile) == iFile)
                {
                    tempGrids.push_back(grids[iZone]);
                }
                else
                {
                    Grid *tempZero = nullptr;
                    tempGrids.push_back(tempZero);
                }
            }

            string out_gfile_temp = AddSymbolToFileName(out_gfile, "_", iFile);
            DumpGrid(out_gfile_temp, &tempGrids[0], numberOfZones);
        }

    }
    else
    {
        DumpGrid(out_gfile, grids, numberOfZones);
    }
}

void GridTranslator::OutputMultiZoneComputationalGrids()
{
    fstream file;
    int tasktype = GetTaskCode();
    if (tasktype == PHSPACE::PARTITION_GRID)
    {
        string partition_grid_file = "grid.fts";
        GlobalDataBase::GetData("partition_grid_file", &partition_grid_file, PHSTRING, 1);
        OpenSerialFile(file, partition_grid_file, ios_base::out|ios_base::binary|ios_base::trunc);
    }
    else if (tasktype == PHSPACE::CREATE_GRID)
    {
        string out_gfile = "grid.fts";
        GlobalDataBase::GetData("out_gfile", &out_gfile, PHSTRING, 1);
        OpenSerialFile(file, out_gfile.c_str(), ios_base::out|ios_base::binary|ios_base::trunc);
    }

    using namespace PHMPI;

    CreateZoneProcessorID(numberOfZones);
    CreateZoneGridID(numberOfZones);
    CreateZoneGridType(numberOfZones);

    int *zoneProcessorIndexContainer = GetZoneProcessorID();
    int *zoneIndexContainer = GetZoneGridID();
    int *zoneTypeContainer  = GetZoneGridType();
    int *block_proc_inp     = GetZoneProcessorID_INP();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        zoneProcessorIndexContainer[iZone] = iZone;
        zoneIndexContainer[iZone] = iZone;
        zoneTypeContainer[iZone] = (*gridGroup)[iZone]->Type();
    }

    if (block_proc_inp)
    {
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            SetZoneProcessorID(iZone, block_proc_inp[iZone]);
        }
    }

    PHWrite(file, numberOfZones);
    PHWrite(file, &zoneProcessorIndexContainer[0], numberOfZones);
    PHWrite(file, &zoneIndexContainer[0], numberOfZones);
    PHWrite(file, &zoneTypeContainer[0], numberOfZones);

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        (*gridGroup)[iZone]->WriteGrid(file);
    }

    CloseFile(file);
}

void GridTranslator::CheckGrid()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        StructGrid *grid = StructGridCast((*gridGroup)[iZone]);
        InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();

        int nIFace = 0;
        if (interfaceInfo)
        {
            nIFace = interfaceInfo->GetNIFace();
        }

        if (nIFace == 0) continue;

        int *interFace2ZoneID      = interfaceInfo->GetInterFace2ZoneID();
        int *interFace2InterFaceID = interfaceInfo->GetInterFace2InterFaceID();

        int target_zone, itface;
        int is, js, ks, nsurf, s_lr, it, jt, kt, t_lr;
        RDouble s_xfc, s_yfc, s_zfc, t_xfc, t_yfc, t_zfc;
        int il1, jl1, kl1;
  
        for (int iFace = 0; iFace < nIFace; ++ iFace)
        {
            target_zone = interFace2ZoneID[iFace];
            if (target_zone < iZone) continue;
            
            grid->GetSourceIndexIJK_Nsurf_LR(iFace, 1, is, js, ks, nsurf, s_lr);
            if (s_lr == 1)
            {
                GetNsurfIndex(nsurf, il1, jl1, kl1);
                is = is + il1;
                js = js + jl1;
                ks = ks + kl1;
            }
            grid->FaceCoor(is, js, ks, nsurf, s_xfc, s_yfc, s_zfc);

            itface = interFace2InterFaceID[iFace];

            if (target_zone == iZone)
            {
                grid->GetSourceIndexIJK_Nsurf_LR(itface, 1, it, jt, kt, nsurf, t_lr);
                if (t_lr == 1)
                {
                    GetNsurfIndex(nsurf, il1, jl1, kl1);
                    it = it + il1;
                    jt = jt + jl1;
                    kt = kt + kl1;
                }
                grid->FaceCoor(it, jt, kt, nsurf, t_xfc, t_yfc, t_zfc);
            }
            else
            {
                StructGrid *tgrid = StructGridCast((*gridGroup)[target_zone]);
                tgrid->GetSourceIndexIJK_Nsurf_LR(itface, 1, it, jt, kt, nsurf, t_lr);
                if (t_lr == 1)
                {
                    GetNsurfIndex(nsurf, il1, jl1, kl1);
                    it = it + il1;
                    jt = jt + jl1;
                    kt = kt + kl1;
                }
                tgrid->FaceCoor(it, jt, kt, nsurf, t_xfc, t_yfc, t_zfc);
            }
            
            if ((fabs(s_xfc - t_xfc) > 1.0e-8) || (fabs(s_yfc - t_yfc) > 1.0e-8) || (fabs(s_zfc - t_zfc) > 1.0e-8))
            {
                TK_Exit::ExceptionExit(" Grid Link Error! ");
            }
        }
    }
}

}


