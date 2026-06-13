#include "Gas.h"
#include "PHMpi.h"
#include "PHIO.h"
#include "Zone.h"
#include "Geo_UnstructBC.h"
#include "Geometry.h"
#include "PerfectGas.h"
#include "Pointer.h"
#include "AleModel.h"
#include "MixGrid.h"
#include "Solver.h"
#include "Glb_Dimension.h"
#include "GridType.h"

#ifdef USE_INCOMSOLVER
#include "IncomGas.h"
#endif

#include "MixedGas.h"

using namespace std;

#pragma warning (disable:913)
namespace PHSPACE
{
vector <Zone *> *GlobalZones;

Zone::Zone(int index)
{
    this->index = index;
    solvers = new vector <PHSolver *>;
    geom    = new PHGeometry(new vector <Grid *>);

    zPara  = new Data_Param();
    zField = new Data_Field();

    nsolver = 0;
}

Zone::~Zone()
{
    for (std::size_t iSolver = 0; iSolver < solvers->size(); ++ iSolver)
    {
        delete (*solvers)[iSolver];
    }
    delete solvers;

    delete geom;
    delete zPara;
    delete zField;

    CleanGlobalValuesSolver();
}

void Zone::AddGridAndCopyZonePara(Grid *grid)
{
    grid->CopyPara(zPara);
    geom->AddGrid(grid);
}

void Zone::AddSolver(PHSolver *solver)
{
    //! Add the solver into the current zone.
    this->solvers->push_back(solver);
    Zone::SetNSolver(static_cast<int>(this->solvers->size()));

    //! Add the solver list of the current zone into the global solver list.
    GlobalSolvers::AddSolverToGlobal(this->GetIndex(), this->solvers);
}

void Zone::AddGridToGlobal(int iZone)
{
    PHSPACE::AddGridToGlobal(iZone, geom->GetVectorGrid());
}

void Zone::AddGrid(Grid *grid)
{
    geom->AddGrid(grid);
}

void Zone::ComputeWeight()
{
    for (int iLevel = 0; iLevel < geom->GetNumberOfMultiGrid(); ++ iLevel)
    {
        geom->GetGrid(iLevel)->ComputeWeight();
    }
}

void Zone::CoarseGrids()
{
    Grid *baseGrid = geom->GetGrid(0);

    int nMGLevel = GlobalDataBase::GetIntParaFromDB("nMGLevel");

    for (int iLevel = 1; iLevel < nMGLevel; ++ iLevel)
    {
        Grid *fineGrid = geom->GetGrid(iLevel - 1);

        fineGrid->CoarseGrids(nMGLevel);

        Grid *coarseGrid = fineGrid->GetCoarseGrid();
        if (coarseGrid)
        {
            if (baseGrid->GetInterfaceInfo())
            {
                InterfaceInfo *cgridInterfaceInfo = new InterfaceInfo(*baseGrid->GetInterfaceInfo());
                coarseGrid->SetInterfaceInfo(cgridInterfaceInfo);
            }
            if (baseGrid->GetInterpointInfo())
            {
                InterpointInformation *cgridInterpointInfo = new InterpointInformation(*baseGrid->GetInterpointInfo());
                coarseGrid->SetInterpointInfo(cgridInterpointInfo);
            }
            AddGridAndCopyZonePara(coarseGrid);
            coarseGrid->InitMovingGrids();
        }
        else
        {
            break;
        }
    }

    for (int iLevel = 1; iLevel < nMGLevel; ++ iLevel)
    {
        Grid *grid = geom->GetGrid(iLevel);
        grid->AllocateOversetGrid();
        //ostringstream oss;
        //oss << "level" << iLevel << ".plt";
        //VisualizationMesh2D( UnstructGridCast(grid), oss.str() );
    }
}

void Zone::SetInterface(Grid *base_grid, Grid *cgrid)
{
    InterfaceInfo *interfaceInfo = base_grid->GetInterfaceInfo();
    cgrid->SetInterfaceInfo(interfaceInfo);
}

void Zone::ShowCoarseGrids()
{
/*
    for (int iLevel = 1; iLevel < geom->GetNumberOfMultiGrid(); ++ iLevel)
    {
        Grid *grid = geom->GetGrid(iLevel);

        if (grid->Type() == PHSPACE::UNSTRUCTGRID)
        {
            std::ostringstream oss;
            oss << "./results/grid" << this->GetIndex() << "_" << "level_" << iLevel << ".plt";

            if (PHSPACE::GetDim() == TWO_D)
            {
                VisualizationMesh2D(UnstructGridCast(grid), oss.str().c_str());
            }
            else
            {
                //! Warning: the bc mesh of coarse grid is same to the fine grid.
                //! Only visualize the eigen faces.
                VisualizationEigenMeshofCoarseGrid3D(UnstructGridCast(grid), oss.str().c_str());
            }
        }
    }
*/
}

void Zone::UpdateAllData()
{
    set < Data_SafeData > *basep = GlobalDataBase::GetBaseP();
    set < Data_SafeData >::iterator it;

    for (it = basep->begin(); it != basep->end(); ++ it)
    {
        this->UpdateData(it->GetName(), it->GetData(), it->GetType(), it->GetSize());
    }
}

void InitGlobalValuesOfAllZones()
{
    InitGlobalValuesSolver();
}

void InitGlobalValuesGrid()
{
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");
    int maxale = GlobalDataBase::GetIntParaFromDB("maxale");

    if (!isUnsteady || !isAle)
    {
        return;
    }

    RDouble *ralepar = new RDouble[maxale];
    GlobalDataBase::UpdateDataPtr("ralepar",ralepar);

    int     *ialepar = new int[maxale];
    GlobalDataBase::UpdateDataPtr("ialepar",ialepar);

    RDouble **angle0 = NewPointer2<RDouble>(3, 3);
    GlobalDataBase::UpdateDataPtr("angle0",angle0);

    RDouble **coord0 = NewPointer2<RDouble>(3, 3);
    GlobalDataBase::UpdateDataPtr("coord0",coord0);

    RDouble  **angle = NewPointer2<RDouble>(3, 3);
    GlobalDataBase::UpdateDataPtr("angle",angle);

    RDouble  **coord = NewPointer2<RDouble>(3, 3);
    GlobalDataBase::UpdateDataPtr("coord",coord);
    
    int iflag = 0;

    for (int i = 0; i < 3; ++ i)
    {
        for (int j = 0; j < 3; ++ j)
        {
            angle0[i][j] = 0.0;
            coord0[i][j] = 0.0;
        }
    }
    
    if (iflag == 0)
    {
        ialepar[0] = 1;
        ialepar[1] = 1;
        ialepar[2] = 0;
        
        coord0[0][1] = - 1.0;
        coord0[1][1] =   0.0;
        coord0[2][1] =   0.0;
    }
    else if (iflag == 2)
    {
        ialepar[0] = 1;
        ialepar[1] = 1;
        ialepar[2] = 5;
    }
    else if (iflag == 3)
    {
        ialepar[0] = 1;
        ialepar[1] = 1;
        ialepar[2] = 6;
    }
    else
    {
        ialepar[0] = 1;
        ialepar[1] = 1;
        ialepar[2] = 5;
    }

    for (int m = 0; m < 3; ++ m)
    {
        coord[m][0] = coord0[m][0];
        coord[m][1] = coord0[m][0];
        coord[m][2] = coord0[m][0];

        angle[m][0] = angle0[m][0];
        angle[m][1] = angle0[m][0];
        angle[m][2] = angle0[m][0];
    }
}

void InitGlobalValuesSolver()
{
    int iapplication = GlobalDataBase::GetIntParaFromDB("iapplication");
    int chemical = GlobalDataBase::GetIntParaFromDB("nchem");
    if (iapplication == 0)
    {
        using namespace GAS_SPACE;
        int compressible = GlobalDataBase::GetIntParaFromDB("compressible");
        if (compressible == INCOMPRESSIBLE)
        {
#ifdef USE_INCOMSOLVER
            gas = new IncomGas();
#endif
        }
        else
        {
            if (chemical == 0)
            {
                gas = new PerfectGas();
            }
            else
            {
                gas = new MixedGas();
            }
        }
        gas->InitGlobalParameterOfNSEquation();
    }
    else if (iapplication == 1)
    {

    }
}

void CleanGlobalValuesSolver()
{
    int iapplication = GlobalDataBase::GetIntParaFromDB("iapplication");

    if (iapplication == 0)
    {
        using namespace GAS_SPACE;
        delete gas;
        gas = 0;
    }
    else if (iapplication == 1)
    {
    }
}

}
