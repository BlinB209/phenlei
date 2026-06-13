#include "Mesh_Refine.h"
#include "TK_Log.h"
#include "Pre_HDF5File.h"

namespace PHSPACE
{

Mesh_Refine::Mesh_Refine(const string &gridFileNameIn)
{
    this->gridFileName = gridFileNameIn;
    numberOfZones = 0;

    region = 0;
    refinedGrids = 0;
    refineParameter = 0;
}

Mesh_Refine::~Mesh_Refine()
{
    delete region;

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        delete refinedGrids[iZone];
    }
    delete [] refinedGrids;
}

void Mesh_Refine::Run()
{
    ReadGrid();

    AllocateMemory();

    ConstructGridTopo();

    BuildRefineProperty();

    FixAnisoRefineType();

    RefineGrid();

    GenerateAndDumpComputationalGrid();
}

void Mesh_Refine::SetRefineParameter(RefineParameter *refineParameterIn)
{
    this->refineParameter = refineParameterIn;
}

void Mesh_Refine::ReadGrid()
{
    GlobalDataBase::UpdateData("gridfile", &gridFileName, PHSTRING, 1);

    region = new Region();
    region->ReadGrid();
    region->SetupMultiBlockRelationship();
}

void Mesh_Refine::DumpComputationalGrid()
{
    string out_gfile;
    GlobalDataBase::GetData("out_gfile", &out_gfile, PHSTRING, 1);

    WriteLogFile(" Grid Writing ...");
    if (PHMPI::GetCurrentProcessorID() == 0) cout << "  Grid Writing ... " << endl;

    if (out_gfile == "") return;

    PrintToWindow("    Writing new HDF5 version grid to file ....\n");
    WriteLogFile("    Writing new HDF5 version grid to file ....\n");

    IO_HDF5Write writeGrid(out_gfile, refinedGrids, numberOfZones);
    writeGrid.Run();

    WriteLogFile(" End Writing Grid!");
}


void RefineParameter::SetAnisoRefineType(int anisoRefineIn)
{
    this->anisoRefine = anisoRefineIn;
}

int RefineParameter::GetAnisoRefineType()
{
    return this->anisoRefine;
}

}


