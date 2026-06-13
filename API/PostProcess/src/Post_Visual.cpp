#include "Post_Visual.h"
#include "TK_Exit.h"
#include "Pointer.h"
#include "IO_FileName.h"
#include "PHIO.h"
#include "Geo_UnstructBC.h"
#include "Constants.h"
#include "Glb_Dimension.h"

#ifdef USE_TecplotLib
#include "TECIO.h"
#endif

#include "Gas.h"

using namespace std;
namespace PHSPACE
{
LIB_EXPORT Post_Visual::Post_Visual(int nVisualVariables_in, int *visualVariables_in, Grid *grid_in, int flowType_in)
{
    visualVariablesPtr = new Data_Field();

    this->grid = grid_in;
    this->flowType = flowType_in;

    if (grid_in)
    {
        gridType = grid_in->Type();
    }
    else
    {
        gridType = -1;
    }

    if (flowType == AverageFlow)
    {
        visualVariablesMap.insert(pair<int, string>(VISUAL_DENSITY , "density_Average"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_U       , "u_Average"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_V       , "v_Average"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_W       , "w_Average"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_PRESSURE, "pressure_Average"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_CP      , "cp_Average"));

        nVisualVariables = 0;
        map<int, string>::iterator varIter;
        for (varIter = visualVariablesMap.begin(); varIter != visualVariablesMap.end(); ++ varIter)
        {
            int varTemp = varIter->first;
            visualVariables.insert(varTemp);
            nVisualVariables ++;
        }
        return;
    }

    if (flowType == AverageReynoldsStress)
    {
        visualVariablesMap.insert(pair<int, string>(VISUAL_TAU_XX, "tau_xx"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_TAU_YY, "tau_yy"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_TAU_ZZ, "tau_zz"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_TAU_XY, "tau_xy"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_TAU_XZ, "tau_xz"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_TAU_YZ, "tau_yz"));

        nVisualVariables = 6;
        for (int iVar = 90; iVar < 96; ++ iVar)
        {
            visualVariables.insert(iVar);
        }
        return;
    }

    visualVariablesMap.insert(pair<int, string>(VISUAL_DENSITY            , "density"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_U                  , "u"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_V                  , "v"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_W                  , "w"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_PRESSURE           , "pressure"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_TEMPERATURE        , "temperature"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_MACH               , "mach"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_VISCOSITY_LAMINAR  , "viscosityLaminar"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_VISCOSITY_TURBULENT, "viscosityTurbulence"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_VORTICITY_X        , "vorticity_x"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_VORTICITY_Y        , "vorticity_y"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_VORTICITY_Z        , "vorticity_z"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_VORTICITY_MAGNITUDE, "vorticityMagnitude"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_STRAIN_RATE        , "strain_rate"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_Q_CRITERIA         , "Q_criteria"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_CP                 , "cp"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_TIME_STEP          , "dt"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_VOLUME             , "volume"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_MODELED_TKE        , "modeledTKE"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_MODELED_DISSIPATION, "modeleddissipationrate"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_SST_F1             , "SSTF1"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_SST_F2             , "SSTF2"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_INTERMITTENCY      , "intermittency"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_MOMENTUMTHICKREYNOLDS, "MomentumThicknessReynolds"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_STREAMLINE_U       , "streamline_u"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_STREAMLINE_V       , "streamline_v"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_STREAMLINE_W       , "streamline_w"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_TRANSITION_GAMAEFF , "gamaeff"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_TRANSITION_RESCF   , "stationaryCrossFlowReynolds"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_GRADIENT_UX        , "gradientUx"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_GRADIENT_UY        , "gradientUy"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_GRADIENT_VX        , "gradientVx"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_GRADIENT_VY        , "gradientVy"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_DIST          , "walldist"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_GAMA               , "gama"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_IBLANK             , "iblank"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_CFL1               , "localCFL"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_CFL2               , "minCFL"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_STREAMLINE_MACH    , "streamline_mach"));

    int compressible = COMPRESSIBLE;
    compressible = GlobalDataBase::GetIntParaFromDB("compressible");
    if (compressible == INCOMPRESSIBLE)
    {
        visualVariablesMap.insert(pair<int, string>(VISUAL_INCOMP_U  , "U" ));
        visualVariablesMap.insert(pair<int, string>(VISUAL_INCOMP_V  , "V" ));
        visualVariablesMap.insert(pair<int, string>(VISUAL_INCOMP_W  , "W" ));
        visualVariablesMap.insert(pair<int, string>(VISUAL_INCOMP_P  , "P" ));
        visualVariablesMap.insert(pair<int, string>(VISUAL_INCOMP_CP , "CP"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_INCOMP_T  , "T" ));
        visualVariablesMap.insert(pair<int, string>(VISUAL_INCOMP_ENTHALPY , "Enthalpy"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_INCOMP_DEN, "DEN"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_INCOMP_VIS, "VIS"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_INCOMP_TE , "TE"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_INCOMP_ED , "ED"));
        if (GlobalDataBase::GetIntParaFromDB("isSolveSpecies"))
        {
            int nSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpeciesIncom");
            string *speciesNameIncom = new string[nSpecies];
            GlobalDataBase::GetData("speciesNameIncom", speciesNameIncom, PHSTRING, 1);

            for (int iSpecie = 0; iSpecie < nSpecies; ++iSpecie)
            {
                int index = VISUAL_INCOMP_SPECIES + iSpecie;
                visualVariablesMap.insert(pair<int, string>(index, speciesNameIncom[iSpecie]));
            }
        }
    }

    visualVariablesMap.insert(pair<int, string>(VISUAL_DIMENSIONAL_DENSITY,     "density_Dim"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_DIMENSIONAL_U      ,     "u_Dim"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_DIMENSIONAL_V      ,     "v_Dim"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_DIMENSIONAL_W      ,     "w_Dim"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_VELOCITY_MAGNITUDE ,    "Vmag_Dim"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_DIMENSIONAL_PRESSURE,    "pressure_Dim"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_DIMENSIONAL_TEMPERATURE, "temperature_Dim"));

    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    if (nChemical > 0)
    {
        if (nTemperatureModel == 2)
        {
            visualVariablesMap.insert(pair<int, string>(VISUAL_TEMPERATURE_VIBRATION, "Tv"));
            visualVariablesMap.insert(pair<int, string>(VISUAL_ENERGY_VIBRATION,      "Ev"));
        }
        else if (nTemperatureModel == 3)
        {
            visualVariablesMap.insert(pair<int, string>(VISUAL_TEMPERATURE_VIBRATION, "Tv"));
            visualVariablesMap.insert(pair<int, string>(VISUAL_TEMPERATURE_ELECTRON,  "Te"));
            visualVariablesMap.insert(pair<int, string>(VISUAL_ENERGY_VIBRATION,      "Ev"));
            visualVariablesMap.insert(pair<int, string>(VISUAL_ENERGY_ELECTRON,       "Ee"));
        }
        visualVariablesMap.insert(pair<int, string>(VISUAL_ELECTRON_NUMBER, "Ne"));
    }

    nVisualVariables = nVisualVariables_in;
    int nCount = 0;
    for (int iVar = 0; iVar < nVisualVariables; ++ iVar)
    {
        int varTemp = visualVariables_in[iVar];
        if (visualVariablesMap.find(varTemp) == visualVariablesMap.end())
        {
            continue;
        }
        visualVariables.insert(varTemp);
        ++ nCount;
    }
    nVisualVariables = nCount;
    GlobalDataBase::UpdateData("nVisualVariables", &nVisualVariables, PHINT, 1);
    
    int *visualVariablesData = new int[nVisualVariables];
    set<int>::iterator varIter = visualVariables.begin();
    for (int iVar = 0; iVar < nVisualVariables; ++ iVar)
    {
        visualVariablesData[iVar] = *varIter;
        varIter ++;
    }
    GlobalDataBase::UpdateData("visualVariables", visualVariablesData, PHINT, nVisualVariables);
    delete [] visualVariablesData;    visualVariablesData = nullptr;
}

LIB_EXPORT Post_Visual::~Post_Visual()
{
    FreeAllVisualNodeVarPtr();
    delete visualVariablesPtr;
}

LIB_EXPORT void Post_Visual::GetASCIITecplotHeader(vector<string> &title_tecplot)
{
    if (nVisualVariables <= 0) return;

    ostringstream oss;
    oss << "title=\"Flow Fields of PHengLEI\"" << endl;
    oss << "variables=\"x\", \"y\", \"z\"";
    set<int>::iterator varIter;
    for (varIter = visualVariables.begin(); varIter != visualVariables.end(); ++ varIter)
    {
        int variableType = *varIter;
        string varName = visualVariablesMap[variableType];
        string varNameTemp = "\"" + varName + "\"";
        oss << ", " << varNameTemp;
    }

    //! Chemical species must be the last!
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        //! Chemical turned on.
        using namespace GAS_SPACE;
        int nm   = GlobalDataBase::GetIntParaFromDB("nm");
        int nl   = GlobalDataBase::GetIntParaFromDB("nl");
        int neqn = nl + nchem;

        string *varname = gas->GetNameOfSpecies();
        for (int m = nm; m < neqn; ++ m)
        {
            string varName = "massfraction-" + varname[m-nm];
            string species_name = "\"" + varName + "\"";
            oss << ", " << species_name;
        }
        for (int m = nm; m < neqn; ++ m)
        {
            string varName = "molefraction-" + varname[m-nm];
            string species_name = "\"" + varName + "\"";
            oss << ", " << species_name;
        }
    }
    title_tecplot.push_back(oss.str());
    return;
}

LIB_EXPORT void Post_Visual::GetBinaryTecplotHeader(vector<string> &title_tecplot)
{
    if (nVisualVariables <= 0) return;

    ostringstream oss;
    oss << "x y z";

    set<int>::iterator varIter;
    for (varIter = visualVariables.begin(); varIter != visualVariables.end(); ++ varIter)
    {
        int variableType = *varIter;
        string varName = visualVariablesMap[variableType];
        //string varNameTemp = "\"" + varName + "\"";
        oss << " " << varName;
    }

    //! Chemical species must be the last!
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        //! Chemical turned on.
        using namespace GAS_SPACE;
        int nm    = GlobalDataBase::GetIntParaFromDB("nm");
        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int neqn  = nl + nchem;

        string *varname = gas->GetNameOfSpecies();
        for (int m = nm; m < neqn; ++ m)
        {
            //string species_name = "\"" + varname[m-nm] + "\"";
            string species_name = "massfraction-" + varname[m-nm];
            oss << " " << species_name;
        }
        for (int m = nm; m < neqn; ++ m)
        {
            string species_name = "molefraction-" + varname[m-nm];
            oss << " " << species_name;
        }
    }

    title_tecplot.push_back(oss.str());
    return;
}

LIB_EXPORT void Post_Visual::UpdateVisualNodeVarPtr(const string &name, RDouble *cellCenterData)
{
    if (!grid)
    {
        TK_Exit::ExceptionExit("Post_Visual::UpdateVisualNodeVarPtr needs Geometry infomation.");
    }

    UnstructGrid *gridUnstruct = UnstructGridCast(grid);
    int nTotalNode = gridUnstruct->GetNTotalNode();
    RDouble *qn = new RDouble[nTotalNode];

    bool isVelocityForPostVisual = false;
    if (name == "streamline_u" || name == "streamline_v" || name == "streamline_w")
    {
        isVelocityForPostVisual = true;
    }

    int compressible = GlobalDataBase::GetIntParaFromDB("compressible");
    if (COMPRESSIBLE == compressible)
    {
        CompNodeVar(gridUnstruct, qn, name, cellCenterData, isVelocityForPostVisual);
    }
    else if (INCOMPRESSIBLE == compressible)
    {
        CompNodeVar(gridUnstruct, qn, cellCenterData, isVelocityForPostVisual);
    }

    visualVariablesPtr->UpdateDataPtr(name, qn);
}

LIB_EXPORT void Post_Visual::UpdateVisualNodeVarPtr(const string &name, int *cellCenterData)
{
    if (!grid)
    {
        TK_Exit::ExceptionExit("Post_Visual::UpdateVisualNodeVarPtr needs Geometry infomation.");
    }

    UnstructGrid * gridUnstruct = UnstructGridCast(grid);
    int nTotalNode = gridUnstruct->GetNTotalNode();
    RDouble *qn = new RDouble[nTotalNode];

    SetField(qn, 1.0, nTotalNode);
    visualVariablesPtr->UpdateDataPtr(name, qn);
}

LIB_EXPORT void Post_Visual::UpdateVisualNodeVarPtr(const string &name, RDouble4D &cellCenterData, int index)
{
    if (!grid)
    {
        TK_Exit::ExceptionExit("Post_Visual::UpdateVisualNodeVarPtr needs Geometry infomation.");
    }

    StructGrid *grid_str = StructGridCast(grid);

    int ni = grid_str->GetNI();
    int nj = grid_str->GetNJ();
    int nk = grid_str->GetNK();

    int nLayer = GetNumberOfGhostCellLayers();

    Range I(1 - nLayer, ni + nLayer - 1);
    Range J(1 - nLayer, nj + nLayer - 1);
    Range K(1 - nLayer, nk + nLayer - 1);
    if (nk == 1) K.setRange(1, 1);
    Range M(0, 0);

    RDouble4D *qn = new RDouble4D(I, J, K, M, fortranArray);
    bool isVelocityForPostVisual = false;

    if (name == "streamline_u" || name == "streamline_v" || name == "streamline_w")
    {
        isVelocityForPostVisual = true;
    }

    CompNodeVar(grid_str, *qn, 0, cellCenterData, index, isVelocityForPostVisual);

    int isWennScheme = GlobalDataBase::GetIntParaFromDB("isWennScheme");
    if (!isVelocityForPostVisual)
    {
        ModifyBoundaryNodeValue(name, grid_str, *qn, cellCenterData, index, PHENGLEI::SOLID_SURFACE);
    }

    ModifyWallTemperature(name, grid_str, *qn);

    visualVariablesPtr->UpdateDataPtr(name, qn);
}

LIB_EXPORT void Post_Visual::UpdateVisualNodeVarPtr(const string &name, RDouble3D &cellCenterData)
{
    if (!grid)
    {
        TK_Exit::ExceptionExit("Post_Visual::UpdateVisualNodeVarPtr needs Geometry infomation.");
    }

    StructGrid *grid_str = StructGridCast(grid);

    int ni = grid_str->GetNI();
    int nj = grid_str->GetNJ();
    int nk = grid_str->GetNK();

    int nLayer = GetNumberOfGhostCellLayers();

    Range I(1 - nLayer, ni + nLayer - 1);
    Range J(1 - nLayer, nj + nLayer - 1);
    Range K(1 - nLayer, nk + nLayer - 1);

    if (nk == 1) K.setRange(1, 1);
    Range M(0, 0);

    RDouble4D *qn = new RDouble4D(I, J, K, M, fortranArray);
    CompNodeVar(grid_str, *qn, 0, cellCenterData);
    ModifyBoundaryNodeValue(name, grid_str, *qn, cellCenterData, PHENGLEI::SOLID_SURFACE);
    visualVariablesPtr->UpdateDataPtr(name, qn);
}

LIB_EXPORT void Post_Visual::UpdateVisualNodeNotInterpolation(const string &name, int *cellCenterdata)
{
    if (!grid)
    {
        TK_Exit::ExceptionExit("Post_Visual::UpdateVisualNodeVarPtr needs Geometry infomation.");
    }

    UnstructGrid *gridUnstruct = UnstructGridCast(grid);
    int           nTotalNode   = gridUnstruct->GetNTotalNode();
    RDouble      *qn           = new RDouble[nTotalNode];

    int  nTotalCell    = gridUnstruct->GetNTotalCell();
    int *cell2node     = gridUnstruct->GetCell2Node();
    int *nodeNumOfCell = gridUnstruct->GetNodeNumberOfEachCell();

    SetField(qn, 0.0, nTotalNode);
    int nodeSubscript = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int nodeNum = nodeNumOfCell[iCell];
        if (-1 != cellCenterdata[iCell])
        {
            nodeSubscript += nodeNum;
            continue;
        }

        for (int iNode = 0; iNode < nodeNum; ++ iNode)
        {
            int nodeIndex = cell2node[nodeSubscript];
            qn[nodeIndex] = -1.0;

            nodeSubscript ++;
        }
    }

    nodeSubscript = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int nodeNum = nodeNumOfCell[iCell];
        if (0 != cellCenterdata[iCell])
        {
            nodeSubscript += nodeNum;
            continue;
        }

        for (int iNode = 0; iNode < nodeNum; ++ iNode)
        {
            int nodeIndex = cell2node[nodeSubscript];
            qn[nodeIndex] = 0.0;

            nodeSubscript ++;
        }
    }

    nodeSubscript = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int nodeNum = nodeNumOfCell[iCell];
        if (1 != cellCenterdata[iCell])
        {
            nodeSubscript += nodeNum;
            continue;
        }

        for (int iNode = 0; iNode < nodeNum; ++ iNode)
        {
            int nodeIndex = cell2node[nodeSubscript];
            qn[nodeIndex] = 1.0;

            nodeSubscript ++;
        }
    }

    visualVariablesPtr->UpdateDataPtr(name, qn);
}

LIB_EXPORT void Post_Visual::UpdateVisualNodeVarPtr(const string &name, Int3D &oversetcellCenterdata)
{
    if (!grid)
    {
        TK_Exit::ExceptionExit("Post_Visual::UpdateVisualNodeVarPtr needs Geometry infomation.");
    }

    StructGrid *grid_str = StructGridCast(grid);

    int ni = grid_str->GetNI();
    int nj = grid_str->GetNJ();
    int nk = grid_str->GetNK();

    Range I(1, ni);
    Range J(1, nj);
    Range K(1, nk);
    if (nk == 1) K.setRange(1, 1);
    Range M(0, 0);

    RDouble4D *qn =  new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D &ql = *qn;
    ql = 0.0;
    int ist, ied, jst, jed, kst, ked;
    grid_str->GetNodeIterationIndex(ist, ied, jst, jed, kst, ked);
    int nDim = GetDim();
    if (nDim == TWO_D)
    {
        int k = 1;
        for (int j = jst; j < jed; ++ j)
        {
            for (int i = ist; i < ied; ++ i)
            {
                if (oversetcellCenterdata(i, j, k) == -1)
                {
                    ql(i,   j,   k, 0) = oversetcellCenterdata(i, j, k);
                    ql(i+1, j,   k, 0) = oversetcellCenterdata(i, j, k);
                    ql(i,   j+1, k, 0) = oversetcellCenterdata(i, j, k);
                    ql(i+1, j+1, k, 0) = oversetcellCenterdata(i, j, k);
                }
            }
        }

        for (int j = jst; j < jed; ++ j)
        {
            for (int i = ist; i < ied; ++ i)
            {
                if (oversetcellCenterdata(i, j, k) == 0)
                {
                    ql(i,   j,   k, 0) = oversetcellCenterdata(i, j, k);
                    ql(i+1, j,   k, 0) = oversetcellCenterdata(i, j, k);
                    ql(i,   j+1, k, 0) = oversetcellCenterdata(i, j, k);
                    ql(i+1, j+1, k, 0) = oversetcellCenterdata(i, j, k);
                }
            }
        }

        for (int j = jst; j < jed; ++ j)
        {
            for (int i = ist; i < ied; ++ i)
            {
                if (oversetcellCenterdata(i, j, k) == 1)
                {
                    ql(i,   j,   k, 0) = oversetcellCenterdata(i, j, k);
                    ql(i+1, j,   k, 0) = oversetcellCenterdata(i, j, k);
                    ql(i,   j+1, k, 0) = oversetcellCenterdata(i, j, k);
                    ql(i+1, j+1, k, 0) = oversetcellCenterdata(i, j, k);
                }
            }
        }
    }
    else if (nDim == THREE_D)
    {
        for (int k = kst; k < ked; ++ k)
        {
            for (int j = jst; j < jed; ++ j)
            {
                for (int i = ist; i < ied; ++ i)
                {
                    if (oversetcellCenterdata(i, j, k) == 2)
                    {
                        ql(i,   j,   k,   0) = oversetcellCenterdata(i, j, k);
                        ql(i+1, j,   k,   0) = oversetcellCenterdata(i, j, k);
                        ql(i,   j+1, k,   0) = oversetcellCenterdata(i, j, k);
                        ql(i,   j,   k+1, 0) = oversetcellCenterdata(i, j, k);
                        ql(i+1, j+1, k,   0) = oversetcellCenterdata(i, j, k);
                        ql(i+1, j,   k+1, 0) = oversetcellCenterdata(i, j, k);
                        ql(i,   j+1, k+1, 0) = oversetcellCenterdata(i, j, k);
                        ql(i+1, j+1, k+1, 0) = oversetcellCenterdata(i, j, k);
                    }
                }
            }
        }

        for (int k = kst; k < ked; ++ k)
        {
            for (int j = jst; j < jed; ++ j)
            {
                for (int i = ist; i < ied; ++ i)
                {
                    if (oversetcellCenterdata(i, j, k) == 0)
                    {
                        ql(i,   j,   k,   0) = oversetcellCenterdata(i, j, k);
                        ql(i+1, j,   k,   0) = oversetcellCenterdata(i, j, k);
                        ql(i,   j+1, k,   0) = oversetcellCenterdata(i, j, k);
                        ql(i,   j,   k+1, 0) = oversetcellCenterdata(i, j, k);
                        ql(i+1, j+1, k,   0) = oversetcellCenterdata(i, j, k);
                        ql(i+1, j,   k+1, 0) = oversetcellCenterdata(i, j, k);
                        ql(i,   j+1, k+1, 0) = oversetcellCenterdata(i, j, k);
                        ql(i+1, j+1, k+1, 0) = oversetcellCenterdata(i, j, k);
                    }
                }
            }
        }

        for (int k = kst; k < ked; ++ k)
        {
            for (int j = jst; j < jed; ++ j)
            {
                for (int i = ist; i < ied; ++ i)
                {
                    if (oversetcellCenterdata(i, j, k) == 1)
                    {
                        ql(i,   j,   k,   0) = oversetcellCenterdata(i, j, k);
                        ql(i+1, j,   k,   0) = oversetcellCenterdata(i, j, k);
                        ql(i,   j+1, k,   0) = oversetcellCenterdata(i, j, k);
                        ql(i,   j,   k+1, 0) = oversetcellCenterdata(i, j, k);
                        ql(i+1, j+1, k,   0) = oversetcellCenterdata(i, j, k);
                        ql(i+1, j,   k+1, 0) = oversetcellCenterdata(i, j, k);
                        ql(i,   j+1, k+1, 0) = oversetcellCenterdata(i, j, k);
                        ql(i+1, j+1, k+1, 0) = oversetcellCenterdata(i, j, k);
                    }
                }
            }
        }
    }

    visualVariablesPtr->UpdateDataPtr(name, qn);
}

void Post_Visual::UpdateVisualNodeVarPtr(const string &name, void *data)
{
    visualVariablesPtr->UpdateDataPtr(name, data);
}

void Post_Visual::UpdateDataPtr(const string &name, void *data)
{
    visualVariablesPtr->UpdateDataPtr(name, data);
}

LIB_EXPORT void Post_Visual::GetAllVisualNodeVarPtr(RDouble **qn)
{
    int varCount = 0;

    for (set<int>::iterator varIter = visualVariables.begin(); varIter != visualVariables.end(); ++ varIter)
    {
        int variableType = *varIter;
        string varName = GetVariableName(variableType);
        qn[varCount++] = (RDouble *)GetVisualNodeVarPtr(varName);
    }

    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        using namespace GAS_SPACE;

        int nm    = GlobalDataBase::GetIntParaFromDB("nm");
        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int neqn  = nl + nchem;

        string *varname = gas->GetNameOfSpecies();
        for (int m = nm; m < neqn; ++ m)
        {
            string varName = "massfraction-" + varname[m-nm];
            qn[varCount++] = (RDouble *)GetVisualNodeVarPtr(varName);
        }
        for (int m = nm; m < neqn; ++ m)
        {
            string varName = "molefraction-" + varname[m-nm];
            qn[varCount++] = (RDouble *)GetVisualNodeVarPtr(varName);
        }
    }
}

LIB_EXPORT void Post_Visual::GetAllVisualNodeVarPtr(RDouble4D **qn)
{
    int varCount = 0;

    for (set<int>::iterator varIter = visualVariables.begin(); varIter != visualVariables.end(); ++ varIter)
    {
        int variableType = *varIter;
        string varName = GetVariableName(variableType);
        qn[varCount++] = (RDouble4D *)GetVisualNodeVarPtr(varName);
    }

    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        using namespace GAS_SPACE;

        int nm    = GlobalDataBase::GetIntParaFromDB("nm");
        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int neqn  = nl + nchem;

        string *varname = gas->GetNameOfSpecies();
        for (int m = nm; m < neqn; ++ m)
        {
            string varName = "massfraction-" + varname[m-nm];
            qn[varCount++] = (RDouble4D *)GetVisualNodeVarPtr(varName);
        }
        for (int m = nm; m < neqn; ++ m)
        {
            string varName = "molefraction-" + varname[m-nm];
            qn[varCount++] = (RDouble4D *)GetVisualNodeVarPtr(varName);
        }
    }
}

void Post_Visual::FreeAllVisualNodeVarPtr()
{
    if (gridType == UNSTRUCTGRID)
    {
        if (visualVariables.size() == 0) return;
        for (set<int>::iterator varIter = visualVariables.begin(); varIter != visualVariables.end(); ++ varIter)
        {
            int variableType = *varIter;
            string varName = GetVariableName(variableType);
            RDouble *qn = (RDouble *)GetVisualNodeVarPtr(varName);
            if (qn)
            {
                delete [] qn;
                //DeleteVisualNodeVarPtr(varName);
            }
        }

        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        if (nchem)
        {
            using namespace GAS_SPACE;

            int nm    = GlobalDataBase::GetIntParaFromDB("nm");
            int nl    = GlobalDataBase::GetIntParaFromDB("nl");
            int neqn  = nl + nchem;

            string *varname = gas->GetNameOfSpecies();
            for (int m = nm; m < neqn; ++ m)
            {
                string varName = "massfraction-" + varname[m-nm];
                RDouble *qn = (RDouble *)GetVisualNodeVarPtr(varName);
                if (qn)
                {
                    delete [] qn;
                }
            }
            for (int m = nm; m < neqn; ++ m)
            {
                string varName = "molefraction-" + varname[m-nm];
                RDouble *qn = (RDouble *)GetVisualNodeVarPtr(varName);
                if (qn)
                {
                    delete [] qn;
                }
            }
        }
    }
    else if (gridType == STRUCTGRID)
    {
        if (visualVariables.size() == 0) return;
        for (set<int>::iterator varIter = visualVariables.begin(); varIter != visualVariables.end(); ++ varIter)
        {
            int variableType = *varIter;
            string varName = GetVariableName(variableType);
            RDouble4D *qn = (RDouble4D *)GetVisualNodeVarPtr(varName);
            if (qn)
            {
                delete qn;
                //DeleteVisualNodeVarPtr(varName);
            }
        }

        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        if (nchem)
        {
            using namespace GAS_SPACE;

            int nm    = GlobalDataBase::GetIntParaFromDB("nm");
            int nl    = GlobalDataBase::GetIntParaFromDB("nl");
            int neqn  = nl + nchem;

            string *varname = gas->GetNameOfSpecies();
            for (int m = nm; m < neqn; ++ m)
            {
                string varName = "massfraction-" + varname[m-nm];
                RDouble4D *qn = (RDouble4D *)GetVisualNodeVarPtr(varName);
                if (qn)
                {
                    delete qn;
                }
            }
            for (int m = nm; m < neqn; ++ m)
            {
                string varName = "molefraction-" + varname[m-nm];
                RDouble4D *qn = (RDouble4D *)GetVisualNodeVarPtr(varName);
                if (qn)
                {
                    delete qn;
                }
            }
        }
    }
    else
    {
        //TK_Exit::ExceptionExit("Post_Visual's grid type is not correctly set.");
    }
}

DumpToVTK::DumpToVTK(vector<DataContainer *> datalist_in, string filename_in, bool CharacteristicBoundary_in, int visualFlowType)
{
    this->vtkdatalist = datalist_in;
    this->vtkfilename = filename_in;
    this->visualFlowType = visualFlowType;
    this->CharacteristicBoundary = CharacteristicBoundary_in;

    for (unsigned int iZone = 0; iZone < vtkdatalist.size(); ++ iZone)
    {
        vtkdatalist[iZone]->MoveToBegin();
    }
    nvarplot = 0;
}

DumpToVTK::~DumpToVTK()
{

}

void DumpToVTK::InitVariables()
{
    if (CharacteristicBoundary == true) return;

    //! Average flow of DES statistics.
    int nVisualVariables = 6;
    int visualVariables[100] = {0, 1, 2, 3, 4, 11};

    if (visualFlowType == 0)
    {
        //! Main flow.
        nVisualVariables = GlobalDataBase::GetIntParaFromDB("nVisualVariables");
        GlobalDataBase::GetData("visualVariables", visualVariables, PHINT, nVisualVariables);
    }

    Post_Visual *postVisualization = new Post_Visual(nVisualVariables, visualVariables);

    int nvarplot = postVisualization->GetNumberofVisualVariables();

    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        using namespace GAS_SPACE;

        int nm    = GlobalDataBase::GetIntParaFromDB("nm");
        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int neqn  = nl + nchem;

        nvarplot = nvarplot + 2 * (neqn - nm);
    }

    SetNvarplot(nvarplot);

    FreePointer(postVisualization);
}

void DumpToVTK::Run()
{
    InitVariables();
    WriteVTKFilm();
}

void DumpToVTK::WriteVTKFilm()
{
    fstream file;
    string vtkname = GetFilename();

    OpenFile(file, vtkname, ios_base::out|ios_base::binary);

    int writeIsOver = 0;
    PHWrite(file, writeIsOver);

    int nVisualVariables = 0;
    int visualVariables[100] = {0, 1, 2, 3, 4, 11};

    if (CharacteristicBoundary == false)
    {
        if (visualFlowType == 0)
        {
            //! Main flow.
            nVisualVariables = GlobalDataBase::GetIntParaFromDB("nVisualVariables");
            GlobalDataBase::GetData("visualVariables", visualVariables, PHINT, nVisualVariables);
        }
    }
    PHWrite(file, nVisualVariables);

    for (int iVisualVariable = 0; iVisualVariable < nVisualVariables; ++ iVisualVariable)
    {
        PHWrite(file, &visualVariables[iVisualVariable], 1);
    }

    if (CharacteristicBoundary == false)
    {
        using namespace GAS_SPACE;
        int numberOfSpecies = gas->GetNumberOfSpecies();
        PHWrite(file, numberOfSpecies);
        int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
        if (nChemical > 0)
        {
            string *varname = gas->GetNameOfSpecies();
            for (int iChemical = 0; iChemical < numberOfSpecies; ++ iChemical)
            {
                string varName = "massfraction-" + varname[iChemical];
                write_string(file, varName);
            }
            for (int iChemical = 0; iChemical < numberOfSpecies; ++ iChemical)
            {
                string varName = "molefraction-" + varname[iChemical];
                write_string(file, varName);
            }
        }
    }

    map <int, string> bcnameMap;
    bcnameMap.insert(pair<int, string>(0 , "NO_BOUNDARY_CONDITION"));
    bcnameMap.insert(pair<int, string>(1 , "EXTRAPOLATION"));
    bcnameMap.insert(pair<int, string>(2 , "WALL"));
    bcnameMap.insert(pair<int, string>(3 , "SYMMETRY"));
    bcnameMap.insert(pair<int, string>(4 , "FARFIELD"));
    bcnameMap.insert(pair<int, string>(5 , "INFLOW"));
    bcnameMap.insert(pair<int, string>(6 , "OUTFLOW"));
    bcnameMap.insert(pair<int, string>(7 , "POLE"));
    bcnameMap.insert(pair<int, string>(8 , "GENERIC_1"));
    bcnameMap.insert(pair<int, string>(9 , "GENERIC_2"));
    bcnameMap.insert(pair<int, string>(10, "GENERIC_3"));
    bcnameMap.insert(pair<int, string>(52, "PRESSURE_INLET"));
    bcnameMap.insert(pair<int, string>(62, "PRESSURE_OUTLET"));
    bcnameMap.insert(pair<int, string>(61, "OUTFLOW_CONFINED"));

    int nZones = static_cast<int>(vtkdatalist.size());
    PHWrite(file, nZones);

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        DataContainer *cdata = vtkdatalist[iZone];

        int GridID;
        int dimension;
        int type;
        int IMxOrNumPts, JMxOrNumElements;

        PHRead(cdata, &GridID, 1);
        PHRead(cdata, &dimension, 1);
        PHRead(cdata, &type, 1);

        PHWrite(file, GridID);
        PHWrite(file, dimension);
        PHWrite(file, type);

        int TecioMission;
        PHRead(cdata, &TecioMission, 1);
        PHWrite(file, TecioMission);

        while (true)
        {
            if (TecioMission == Writecomplete) break;
            if (TecioMission != WriteBoundary && TecioMission != WriteBlock)
            {
                TK_Exit::ExceptionExit("Error: vtkwrite code is wrong!");
            }

            int bctype;

            string bcName;
            if (TecioMission == WriteBoundary)
            {
                PHRead(cdata, &bctype, 1);
                PHWrite(file, bctype);

                cdata->ReadString(bcName);
                
                if (bcName.empty())
                {
                    bcName = bcnameMap[bctype];
                }

                write_string(file, bcName);
            }

            PHRead(cdata, &IMxOrNumPts, 1);
            PHRead(cdata, &JMxOrNumElements, 1);

            PHWrite(file, IMxOrNumPts);
            PHWrite(file, JMxOrNumElements);


            RDouble *x = new RDouble [IMxOrNumPts];
            RDouble *y = new RDouble [IMxOrNumPts];
            RDouble *z = new RDouble [IMxOrNumPts];

            PHRead(cdata, x, IMxOrNumPts);
            PHRead(cdata, y, IMxOrNumPts);
            PHRead(cdata, z, IMxOrNumPts);

            PHWrite(file, x, IMxOrNumPts);
            PHWrite(file, y, IMxOrNumPts);
            PHWrite(file, z, IMxOrNumPts);

            delete [] x;    x = nullptr;
            delete [] y;    y = nullptr;
            delete [] z;    z = nullptr;

            if (nVisualVariables != 0)
            {
                for (int m = 0; m < nvarplot; ++ m)
                {
                    RDouble *qq = new RDouble [IMxOrNumPts];
                    PHRead(cdata, qq, IMxOrNumPts);
                    PHWrite(file, qq, IMxOrNumPts);

                    delete [] qq;    qq = nullptr;
                }
            }

            if (TecioMission == WriteBlock && dimension == THREE_D)
            {
                int nTotalFace = 0;
                PHRead(cdata, &nTotalFace, 1);
                PHWrite(file, nTotalFace);

                int *nodeNumberOfEachCell = new int[JMxOrNumElements];
                PHRead(cdata, nodeNumberOfEachCell, JMxOrNumElements);
                PHWrite(file, nodeNumberOfEachCell, JMxOrNumElements);

                int count = 0;
                for (int iCell = 0; iCell < JMxOrNumElements; ++ iCell)
                {
                    count += nodeNumberOfEachCell[iCell];
                }
                int *cell2Node = new int[count];
                PHRead(cdata, cell2Node, count);
                PHWrite(file, cell2Node, count);

                int *faceNumberOfEachCell = new int[JMxOrNumElements];
                PHRead(cdata, faceNumberOfEachCell, JMxOrNumElements);
                PHWrite(file, faceNumberOfEachCell, JMxOrNumElements);

                count = 0;
                for (int iCell = 0; iCell < JMxOrNumElements; ++ iCell)
                {
                    count += faceNumberOfEachCell[iCell];
                }
                int *cell2Face = new int[count];
                PHRead(cdata, cell2Face, count);
                PHWrite(file, cell2Face, count);

                int *nodeNumberOfEachFace = new int[nTotalFace];
                PHRead(cdata, nodeNumberOfEachFace, nTotalFace);
                PHWrite(file, nodeNumberOfEachFace, nTotalFace);

                count = 0;
                for (int iFace = 0; iFace < nTotalFace; ++ iFace)
                {
                    count += nodeNumberOfEachFace[iFace];
                }

                int *face2Node = new int[count];
                PHRead(cdata, face2Node, count);
                PHWrite(file, face2Node, count);

                delete [] nodeNumberOfEachCell;    nodeNumberOfEachCell = nullptr;
                delete [] cell2Node;    cell2Node = nullptr;
                delete [] faceNumberOfEachCell;    faceNumberOfEachCell = nullptr;
                delete [] cell2Face;    cell2Face = nullptr;
                delete [] nodeNumberOfEachFace;    nodeNumberOfEachFace = nullptr;
                delete [] face2Node;    face2Node = nullptr;
            }
            else
            {
                int *face2node = new int [JMxOrNumElements * 4];
                PHRead(cdata, face2node, JMxOrNumElements * 4);
                PHWrite(file, face2node, JMxOrNumElements * 4);
                delete [] face2node;    face2node = nullptr;
            }

            PHRead(cdata, &TecioMission, 1);
            PHWrite(file, TecioMission);
        }
    }

    file.seekp(file.tellp());
    writeIsOver = 1;
    PHWrite(file, writeIsOver);

    PHSPACE::CloseFile(file);
}

#ifdef USE_TecplotLib
DumpToTecio::DumpToTecio(vector<DataContainer *> datalist_in, string filename_in, int visualFlowType)
{
    this->datalist = datalist_in;
    this->filename = filename_in;
    this->visualFlowType = visualFlowType;

    for (unsigned int iZone = 0; iZone < datalist.size(); ++ iZone)
    {
        datalist[iZone]->MoveToBegin();
    }
    nvarplot = 0;
}

DumpToTecio::~DumpToTecio()
{

}

void DumpToTecio::Run()
{
    InitVariables();
    WriteGridFilm(OUTPUT_SURFACE);
    WriteGridFilm(OUTPUT_BLK);    //! Seperate surface and BLK data into two different files.
    WriteSliceFilm();
}

void DumpToTecio::WriteGridFilm(int isOutputBLK)
{
    if (GetFilename() == "")
    {
        return;
    }

    ostringstream oss0;
    if (GetDim() == THREE_D && GlobalDataBase::GetIntParaFromDB("plotFieldType") == TEC_SPACE::FieldVisual)
    {
        if (OUTPUT_BLK == isOutputBLK)
        {
            if (! PHMPI::IsParallelRun())
            {
                oss0 << PHSPACE::AddSymbolToFileName(GetFilename(), '_', "BLK");
            }
            else
            {
                int fileID = PHMPI::GetFileIndexofCurrentProcessor();
                string finalFileName = PHSPACE::AddSymbolToFileName(GetFilename(), "_BLK_", fileID);
                oss0 << finalFileName;
            }
        }
        else
        {
            if (!PHMPI::IsParallelRun())
            {
                oss0 << GetFilename();
            }
            else
            {
                int fileID = PHMPI::GetFileIndexofCurrentProcessor();
                string finalFileName = PHSPACE::AddSymbolToFileName(GetFilename(), '_', fileID);
                oss0 << finalFileName;
            }
        }
    }
    else
    {
        if (OUTPUT_BLK == isOutputBLK) return;    //! Write a single file in Two_D mode or when plotFieldType set to 0.

        if (! PHMPI::IsParallelRun())
        {
            oss0 << GetFilename();
        }
        else
        {
            int fileID = PHMPI::GetFileIndexofCurrentProcessor();
            string finalFileName = PHSPACE::AddSymbolToFileName(GetFilename(), '_', fileID);
            oss0 << finalFileName;
        }
    }    

    string pltname = oss0.str();
    const char *filenames = pltname.c_str();
    
    string Variable = GetVariables();
    const char *Variables = Variable.c_str();

    map<int, string> bcnameMap;
    bcnameMap.insert(pair<int, string>(0 , "NO_BOUNDARY_CONDITION"));
    bcnameMap.insert(pair<int, string>(1 , "EXTRAPOLATION"));
    bcnameMap.insert(pair<int, string>(2 , "WALL"));
    bcnameMap.insert(pair<int, string>(3 , "SYMMETRY"));
    bcnameMap.insert(pair<int, string>(4 , "FARFIELD"));
    bcnameMap.insert(pair<int, string>(5 , "INFLOW"));
    bcnameMap.insert(pair<int, string>(6 , "OUTFLOW"));
    bcnameMap.insert(pair<int, string>(7 , "POLE"));
    bcnameMap.insert(pair<int, string>(8 , "GENERIC_1"));
    bcnameMap.insert(pair<int, string>(9 , "GENERIC_2"));
    bcnameMap.insert(pair<int, string>(10, "GENERIC_3"));
    bcnameMap.insert(pair<int, string>(52, "PRESSURE_INLET"));
    bcnameMap.insert(pair<int, string>(62, "PRESSURE_OUTLET"));
    bcnameMap.insert(pair<int, string>(61, "OUTFLOW_CONFINED"));

    INTEGER4 Debug     = 0;
    INTEGER4 VIsDouble = 1;
    INTEGER4 FileType  = 0;
    INTEGER4 I;

    I = TECINI112((char *) "PHengLEI Grid Videotex",
                  (char *) Variables,
                  (char *) filenames,
                  (char *) ".",
                  &FileType,
                  &Debug,
                  &VIsDouble);

    bool FileIsFull = false;

    uint_t nzones = datalist.size();
    for (int iZone = 0; iZone < nzones; ++ iZone)
    {
        DataContainer *cdata = datalist[iZone];
        cdata->MoveToBegin();

        INTEGER4 IMxOrNumPts, JMxOrNumElements, KMxOrNumFaces;
        INTEGER4 ZoneType;
               /*ZoneType:      0 -- Ordered
                                1 -- Lineseg
                                2 -- Triangle
                                3 -- Quadrilateral
                                4 -- tetrahedron
                                5 -- Brick
                                6 -- Polygon
                                7 -- Polyhedron
                */
        int GridID;
        int dimension;
        int type;

        PHRead(cdata, &GridID, 1);
        PHRead(cdata, &dimension, 1);
        PHRead(cdata, &type, 1);

        int TecioMission;
        PHRead(cdata, &TecioMission, 1);

        if (GetDim() == THREE_D && GlobalDataBase::GetIntParaFromDB("plotFieldType") == TEC_SPACE::FieldVisual)
        {
            if (OUTPUT_BLK == isOutputBLK)
            {
                if (TecioMission != WriteBlock) continue;    //! Output BLK data only.
            }
            else
            {
                if (TecioMission == WriteBlock) continue;    //! Output surface data.
            }
        }
        else
        {
            //! BLK and surface are merged to output in the other cases.
        }

        while (true)
        {
            if (TecioMission == Writecomplete) break;

            if (TecioMission != WriteBoundary && TecioMission != WriteBlock)
            {
                TK_Exit::ExceptionExit("Error: teciowrite code is wrong!");
            }

            int bctype;
            string bcName;

            if (type == PHSPACE::STRUCTGRID)
            {
                ZoneType = TEC_SPACE::ORDERED;
            }
            else
            {
                if (TecioMission == WriteBoundary)
                {
                    ZoneType = TEC_SPACE::FEQUADRILATERAL;
                }
                else if (dimension == 2)
                {
                    ZoneType = TEC_SPACE::FEPOLYGON;
                }
                else
                {
                    ZoneType = TEC_SPACE::FEPOLYHEDRON;
                }
            }

            if (TecioMission == WriteBoundary)
            {
                PHRead(cdata, &bctype, 1);
                cdata->ReadString(bcName);

                if (bcName.empty())
                {
                    bcName = bcnameMap[bctype];
                }
            }

            PHRead(cdata, &IMxOrNumPts, 1);
            PHRead(cdata, &JMxOrNumElements, 1);
            PHRead(cdata, &KMxOrNumFaces, 1);

            INTEGER4 ICellMax           = 0;
            INTEGER4 JCellMax           = 0;
            INTEGER4 KCellMax           = 0;
            double   SolutionTime       = 360.0;
            INTEGER4 StrandID           = 0;
            INTEGER4 ParentZone         = 0;
            INTEGER4 IsBlock            = 1;
            INTEGER4 NumFaceConnections = 0;
            INTEGER4 FaceNeighborMode   = 0;

            INTEGER4 TotalNumFaceNodes_Rect = 0;
            INTEGER4 *node_number_of_each_face;

            if (type == PHSPACE::UNSTRUCTGRID)
            {
                if (TecioMission == WriteBoundary)
                {
                    PHRead(cdata, &TotalNumFaceNodes_Rect, 1);
                    node_number_of_each_face = new INTEGER4 [JMxOrNumElements];
                    PHRead(cdata, node_number_of_each_face, JMxOrNumElements);
                }
                else
                {
                    PHRead(cdata, &TotalNumFaceNodes_Rect, 1);
                    node_number_of_each_face = new INTEGER4 [KMxOrNumFaces];
                    PHRead(cdata, node_number_of_each_face, KMxOrNumFaces);
                }
            }
            INTEGER4 NumConnBndryFaces_Rect  = 0;    /* interface */
            INTEGER4 TotalNumBndryConns_Rect = 0;
            INTEGER4 SharConn                = 0;

            ostringstream oss;

            if (TecioMission == WriteBlock)
            {
                oss << "Zone" << GridID << " BLK";
            }
            else
            {
                oss << "Zone" << GridID << " " << bcName;
            }

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
                          NULL,
                          NULL,
                          &SharConn);

            INTEGER4 nTotalNode;
            PHRead(cdata, &nTotalNode, 1);

            RDouble *x = new RDouble [nTotalNode];
            RDouble *y = new RDouble [nTotalNode];
            RDouble *z = new RDouble [nTotalNode];
            PHRead(cdata, x, nTotalNode);
            PHRead(cdata, y, nTotalNode);
            PHRead(cdata, z, nTotalNode);

            INTEGER4 IsDouble = 1;
            I = TECDAT112(&nTotalNode, x, &IsDouble);
            I = TECDAT112(&nTotalNode, y, &IsDouble);
            I = TECDAT112(&nTotalNode, z, &IsDouble);
            
            int nvarplot = GetNvarplot();

            for (int m = 0; m < nvarplot; ++ m)
            {
                RDouble *qq = new RDouble [nTotalNode];
                PHRead(cdata, qq, nTotalNode);
                I = TECDAT112(&nTotalNode, qq, &IsDouble);
                delete [] qq;
            }
            
            if (FileIsFull == false) FileIsFull = true;

            if (type == PHSPACE::UNSTRUCTGRID)
            {
                int *face2node = new INTEGER4 [TotalNumFaceNodes_Rect];
                PHRead(cdata, face2node, TotalNumFaceNodes_Rect);

                INTEGER4 *FaceNodes_Rect = new INTEGER4 [TotalNumFaceNodes_Rect];
                for (INTEGER4 i = 0; i < TotalNumFaceNodes_Rect; ++ i)
                {
                    FaceNodes_Rect[i] = face2node[i] + 1;
                }

                if (TecioMission == WriteBoundary)
                {
                    I = TECNOD112(FaceNodes_Rect);
                }
                else
                {
                    INTEGER4 *FaceLeftElems  = new INTEGER4[KMxOrNumFaces];
                    INTEGER4 *FaceRightElems = new INTEGER4[KMxOrNumFaces];

                    int *left_cell_of_face  = new INTEGER4[KMxOrNumFaces];
                    int *right_cell_of_face = new INTEGER4[KMxOrNumFaces];
                    PHRead(cdata, left_cell_of_face, KMxOrNumFaces);
                    PHRead(cdata, right_cell_of_face, KMxOrNumFaces);

                    for (INTEGER4 iFace = 0; iFace < KMxOrNumFaces; ++ iFace)
                    {
                        FaceLeftElems [iFace] = left_cell_of_face [iFace] + 1;
                        FaceRightElems[iFace] = right_cell_of_face[iFace] + 1;
                        if (FaceRightElems[iFace] > JMxOrNumElements || FaceRightElems[iFace] < 0) FaceRightElems[iFace] = 0;
                    }

                    INTEGER4 FaceBndryConnectionCounts = 0;
                    INTEGER4 FaceBndryConnectionElems  = 0;
                    INTEGER4 FaceBndryConnectionZones  = 0;
                
                    I = TECPOLY112(node_number_of_each_face,
                                   FaceNodes_Rect,
                                   FaceLeftElems,
                                   FaceRightElems,
                                   &FaceBndryConnectionCounts,
                                   &FaceBndryConnectionElems,
                                   &FaceBndryConnectionZones);
                
                    delete [] FaceLeftElems;
                    delete [] FaceRightElems;
                    delete [] left_cell_of_face;
                    delete [] right_cell_of_face;
                }
                
                delete [] FaceNodes_Rect;
                delete [] face2node;
                delete [] node_number_of_each_face;
            }

            delete [] x;
            delete [] y;
            delete [] z;

            PHRead(cdata, &TecioMission, 1);
        }
    }

    I = TECEND112();

    if (FileIsFull == false) remove(filenames);
    return;
}

void DumpToTecio::WriteSliceFilm()
{
    int visualSlice = GlobalDataBase::GetIntParaFromDB("visualSlice");

    if (GetDim() != THREE_D || visualSlice != 1) return;


    string sliceFile = GlobalDataBase::GetStrParaFromDB("sliceFile");

    if (visualFlowType == 1)
    {
        //! Average flow.
        string tempFile = AddSymbolToFileName(sliceFile, "_", "Average");
        sliceFile = tempFile;
    }

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady)
    {
        int outnstep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
        sliceFile = PHSPACE::AddSymbolToFileName(sliceFile, '_', outnstep);
    }
    int fileID = PHMPI::GetFileIndexofCurrentProcessor();
    sliceFile = PHSPACE::AddSymbolToFileName(sliceFile, '_', fileID);
    
    std::ostringstream oss;
    oss << "title=\"Slice of PHengLEI\"" << endl;

    string Variable = GetVariables();
    oss << "variables=" << Variable << "\n";
    
    int totalNumber = 0;
    uint_t nZones = datalist.size();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        DataContainer *cdata = datalist[iZone];

        int NumPts = 0;
        PHRead(cdata, &NumPts, 1);

        totalNumber += NumPts;
        if (NumPts == 0) continue;

        RDouble *xSlice = new RDouble [NumPts];
        RDouble *ySlice = new RDouble [NumPts];
        RDouble *zSlice = new RDouble [NumPts];

        PHRead(cdata, xSlice, NumPts);
        PHRead(cdata, ySlice, NumPts);
        PHRead(cdata, zSlice, NumPts);

        int nVarplot = GetNvarplot();
        vector <RDouble *> qnSlice;
        qnSlice.resize(nVarplot);

        for (int iVarplot = 0; iVarplot < nVarplot; ++ iVarplot)
        {
            qnSlice[iVarplot] = new RDouble [NumPts];
            PHRead(cdata, qnSlice[iVarplot], NumPts);
        }

        int wordwidth = 20;
        for (int iNode = 0; iNode < NumPts; ++ iNode)
        {
            oss << setiosflags(ios::left);
            oss << setiosflags(ios::scientific);
            oss << setprecision(10);
            oss << setw(wordwidth) << xSlice[iNode]
                << setw(wordwidth) << ySlice[iNode]
                << setw(wordwidth) << zSlice[iNode];
            for (int iVarplot = 0; iVarplot < nVarplot; ++ iVarplot)
            {
                oss << setw(wordwidth) << qnSlice[iVarplot][iNode];
            }
            oss << "\n";
        }

        delete [] xSlice;
        delete [] ySlice;
        delete [] zSlice;
        for (int iVarplot = 0; iVarplot < nVarplot; ++ iVarplot)
        {
            delete [] qnSlice[iVarplot];
        }
    }

    if (totalNumber != 0)
    {
        fstream file;
        ParallelOpenFile(file, sliceFile, ios_base::out|ios_base::trunc);

        file << oss.str();

        ParallelCloseFile(file);
    }
}

void DumpToTecio::InitVariables()
{
    //! Average flow of DES statistics.
    int nVisualVariables = 6;
    int visualVariables[100] = {0, 1, 2, 3, 4, 11};

    if (visualFlowType == 0)
    {
        //! Main flow.
        nVisualVariables = GlobalDataBase::GetIntParaFromDB("nVisualVariables");
        GlobalDataBase::GetData("visualVariables", visualVariables, PHINT, nVisualVariables);
    }

    Post_Visual *postVisualization = new Post_Visual(nVisualVariables, visualVariables, NULL, visualFlowType);

    int nvarplot = postVisualization->GetNumberofVisualVariables();

    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        using namespace GAS_SPACE;

        int nm    = GlobalDataBase::GetIntParaFromDB("nm");
        int nl    = GlobalDataBase::GetIntParaFromDB("nl");
        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        int neqn  = nl + nchem;

        nvarplot = nvarplot + 2 * (neqn - nm);
    }

    SetNvarplot(nvarplot);

    vector<string> title_tecplot;
    postVisualization->GetBinaryTecplotHeader(title_tecplot);
    SetVariables(title_tecplot.front());

    FreePointer(postVisualization);
}
#endif

DumpToTecplotASCII::DumpToTecplotASCII(vector<DataContainer *> datalist_in, string filename_in, int visualFlowType)
{
    this->dataList = datalist_in;
    this->fileName = filename_in;
    this->visualFlowType = visualFlowType;

    for (unsigned int iZone = 0; iZone < dataList.size(); ++ iZone)
    {
        dataList[iZone]->MoveToBegin();
    }
    nVariable = 0;
}

DumpToTecplotASCII::~DumpToTecplotASCII()
{

}

void DumpToTecplotASCII::InitVariables()
{
    //! Average flow of DES statistics.
    int nVisualVariables = 6;
    int visualVariables[100] = {0, 1, 2, 3, 4, 11};
    nVisualVariables = GlobalDataBase::GetIntParaFromDB("nVisualVariables");
    GlobalDataBase::GetData("visualVariables", visualVariables, PHINT, nVisualVariables);

    Post_Visual *postVisualization = new Post_Visual(nVisualVariables, visualVariables, NULL, visualFlowType);

    int nvarplot = postVisualization->GetNumberofVisualVariables();

    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        using namespace GAS_SPACE;

        int nm   = GlobalDataBase::GetIntParaFromDB("nm");
        int nl   = GlobalDataBase::GetIntParaFromDB("nl");
        int neqn = nl + nchem;

        nvarplot = nvarplot + 2 * (neqn - nm);
    }

    SetNvarplot(nvarplot);

    vector<string> title_tecplot;
    postVisualization->GetASCIITecplotHeader(title_tecplot);
    SetVariables(title_tecplot.front());

    FreePointer(postVisualization);
}

void DumpToTecplotASCII::Run()
{
    InitVariables();
    WriteFlowFilm();
}

void DumpToTecplotASCII::WriteFlowFilm()
{
    if (GetFilename() == "")
    {
        return;
    }

    ostringstream ossFileName;
    ostringstream ossBlockFileName;    //! Seperate boundary and block data to two different files.
    if (! PHMPI::IsParallelRun())
    {
        ossFileName << GetFilename();
        ossBlockFileName << PHSPACE::AddSymbolToFileName(GetFilename(), '_', "BLK");
    }
    else
    {
        int fileID = PHMPI::GetFileIndexofCurrentProcessor();
        string finalFileName = PHSPACE::AddSymbolToFileName(GetFilename(), '_', fileID);
        string finalBlockFileName = PHSPACE::AddSymbolToFileName(GetFilename(), "_BLK_", fileID);
        ossFileName << finalFileName;
        ossBlockFileName << finalBlockFileName;
    }
    string flowFileName = ossFileName.str();
    string flowBlockFileName = ossBlockFileName.str();

    fstream flowFile;
    OpenFile(flowFile, flowFileName, ios_base::out|ios_base::trunc);
    flowFile << setiosflags(ios::left);
    flowFile << setiosflags(ios::scientific);
    flowFile << setprecision(10);
    flowFile << Variables << endl;

    fstream flowBlockFile;
    if (GetDim() == THREE_D && GlobalDataBase::GetIntParaFromDB("plotFieldType") == TEC_SPACE::FieldVisual)
    {
        OpenFile(flowBlockFile, flowBlockFileName, ios_base::out|ios_base::trunc);
        flowBlockFile << setiosflags(ios::left);
        flowBlockFile << setiosflags(ios::scientific);
        flowBlockFile << setprecision(10);
        flowBlockFile << Variables << endl;
    }

    map <int, string> bcNameMap;
    bcNameMap.insert(pair<int, string>(0 , "NO_BOUNDARY_CONDITION"));
    bcNameMap.insert(pair<int, string>(1 , "EXTRAPOLATION"));
    bcNameMap.insert(pair<int, string>(2 , "WALL"));
    bcNameMap.insert(pair<int, string>(3 , "SYMMETRY"));
    bcNameMap.insert(pair<int, string>(4 , "FARFIELD"));
    bcNameMap.insert(pair<int, string>(5 , "INFLOW"));
    bcNameMap.insert(pair<int, string>(6 , "OUTFLOW"));
    bcNameMap.insert(pair<int, string>(7 , "POLE"));
    bcNameMap.insert(pair<int, string>(8 , "GENERIC_1"));
    bcNameMap.insert(pair<int, string>(9 , "GENERIC_2"));
    bcNameMap.insert(pair<int, string>(10, "GENERIC_3"));
    bcNameMap.insert(pair<int, string>(52, "PRESSURE_INLET"));
    bcNameMap.insert(pair<int, string>(62, "PRESSURE_OUTLET"));
    bcNameMap.insert(pair<int, string>(61, "OUTFLOW_CONFINED"));

    uint_t nzones = dataList.size();
    for (int iZone = 0; iZone < nzones; ++ iZone)
    {
        DataContainer *cdata = dataList[iZone];

        int GridID, dimension, type, dumpMission;

        PHRead(cdata, &GridID, 1);
        PHRead(cdata, &dimension, 1);
        PHRead(cdata, &type, 1);
        PHRead(cdata, &dumpMission, 1);

        while (true)
        {
            if (dumpMission == Writecomplete) break;

            if (dumpMission != WriteBoundary && dumpMission != WriteBlock)
            {
                TK_Exit::ExceptionExit("Error: flow file write code is wrong!");
            }

            int bcType;
            string bcName;
            if (dumpMission == WriteBoundary)
            {
                PHRead(cdata, &bcType, 1);
                cdata->ReadString(bcName);

                if (bcName.empty())
                {
                    bcName = bcNameMap[bcType];
                }
            }

            ostringstream ossZoneTitle;
            if (dumpMission == WriteBlock)
            {
                ossZoneTitle << "\"Zone" << GridID << " BLK" << "\"";
            }
            else
            {
                ossZoneTitle << "\"Zone" << GridID << " " << bcName << "\"";
            }

            if (type == PHSPACE::STRUCTGRID)
            {
                if (dumpMission == WriteBlock)
                {
                    if (GetDim() == THREE_D)
                    {
                        DumpOrderedFlow(flowBlockFile, cdata, ossZoneTitle.str());
                    }
                    else
                    {
                        DumpOrderedFlow(flowFile, cdata, ossZoneTitle.str());
                    }
                }
                else
                {
                    DumpOrderedFlow(flowFile, cdata, ossZoneTitle.str());
                }
            }
            else if (type == PHSPACE::UNSTRUCTGRID)
            {
                if (dumpMission == WriteBoundary)
                {
                    DumpFlowByCellTopo(flowFile, cdata, ossZoneTitle.str());
                }
                else
                {
                    if (GetDim() == THREE_D)
                    {
                        DumpFlowByFaceTopo(flowBlockFile, cdata, ossZoneTitle.str());
                    }
                    else
                    {
                        DumpFlowByFaceTopo(flowFile, cdata, ossZoneTitle.str());
                    }
                }
            }
            else
            {
                TK_Exit::ExceptionExit("Error: flow file write code is wrong! Grid type is not right");
            }

            PHRead(cdata, &dumpMission, 1);
        }
    }

    CloseFile(flowFile);
    if (GetDim() == THREE_D && GlobalDataBase::GetIntParaFromDB("plotFieldType") == TEC_SPACE::FieldVisual)
    {
        CloseFile(flowBlockFile);
    }
}

void DumpToTecplotASCII::DumpFlowByFaceTopo(fstream &flowFile, DataContainer *data, string zoneTitle)
{
    int NumPts, NumElements, NumFaces, TotalNumFaceNodes_Rect;

    PHRead(data, &NumPts, 1);
    PHRead(data, &NumElements, 1);
    PHRead(data, &NumFaces, 1);
    PHRead(data, &TotalNumFaceNodes_Rect, 1);

    int *node_number_of_each_face = new int [NumFaces];
    PHRead(data, node_number_of_each_face, NumFaces);

    int nTotalNode;
    PHRead(data, &nTotalNode, 1);

    RDouble *x = new RDouble [nTotalNode];
    RDouble *y = new RDouble [nTotalNode];
    RDouble *z = new RDouble [nTotalNode];
    PHRead(data, x, nTotalNode);
    PHRead(data, y, nTotalNode);
    PHRead(data, z, nTotalNode);

    RDouble **qFiled = NewPointer2<RDouble> (nVariable, nTotalNode);
    for (int iVariable = 0; iVariable < nVariable; ++ iVariable)
    {
        PHRead(data, qFiled[iVariable], nTotalNode);
    }

    int *face2node = new int [TotalNumFaceNodes_Rect];
    PHRead(data, face2node, TotalNumFaceNodes_Rect);

    int *left_cell_of_face  = new int [NumFaces];
    int *right_cell_of_face = new int [NumFaces];
    PHRead(data, left_cell_of_face , NumFaces);
    PHRead(data, right_cell_of_face, NumFaces);

    flowFile << "zone T = " << zoneTitle << "\n";
    if (PHSPACE::GetDim() == TWO_D)
    {
        flowFile << "ZoneType = FEPolygon\n";
    }
    else
    {
        flowFile << "ZoneType = FEPolyhedron\n";
    }
    
    flowFile << "Nodes    = " << nTotalNode << "\n";
    flowFile << "Faces    = " << NumFaces << "\n";
    flowFile << "Elements = " << NumElements << "\n";
    flowFile << "TotalNumFaceNodes = " << TotalNumFaceNodes_Rect << "\n";
    flowFile << "NumConnectedBoundaryFaces = 0\n";
    flowFile << "TotalNumBoundaryConnections = 0\n";

    int wordWidth = 20;
    int nWordOfLine = 5;
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        flowFile << setw(wordWidth) << x[iNode];
        if ((iNode + 1) % nWordOfLine == 0) 
        {
            flowFile << "\n";
        }
    }
    if (nTotalNode % nWordOfLine != 0) flowFile << "\n";

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        flowFile << setw(wordWidth) << y[iNode];
        if ((iNode + 1) % nWordOfLine == 0) 
        {
            flowFile << "\n";
        }
    }
    if (nTotalNode % nWordOfLine != 0) flowFile << "\n";

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        flowFile << setw(wordWidth) << z[iNode];
        if ((iNode + 1) % nWordOfLine == 0) 
        {
            flowFile << "\n";
        }
    }
    if (nTotalNode % nWordOfLine != 0) flowFile << "\n";

    for (int iVariable = 0; iVariable < nVariable; ++ iVariable)
    {
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            flowFile << setw(wordWidth) << qFiled[iVariable][iNode];
            if ((iNode + 1) % nWordOfLine == 0) 
            {
                flowFile << "\n";
            }
        }
        if (nTotalNode % nWordOfLine != 0) flowFile << "\n";
    }

    if (PHSPACE::GetDim() == THREE_D)
    {
        for (int iFace = 0; iFace < NumFaces; ++ iFace)
        {
            flowFile << node_number_of_each_face[iFace] << " ";
            if ((iFace + 1) % nWordOfLine == 0) flowFile << "\n";
        }
        if (NumFaces % nWordOfLine != 0) flowFile << "\n";
    }

    int count = 0;
    for (int iFace = 0; iFace < NumFaces; ++ iFace)
    {
        for (int iNode = 0; iNode < node_number_of_each_face[iFace]; ++ iNode)
        {
            flowFile << face2node[count] + 1 << " ";
            count ++;
        }
        flowFile << "\n";
    }

    for (int iFace = 0; iFace < NumFaces; ++ iFace)
    {
        flowFile << left_cell_of_face[iFace] + 1 << " ";
        if ((iFace + 1) % nWordOfLine == 0) flowFile << "\n";
    }
    if (NumFaces % nWordOfLine != 0) flowFile << "\n";

    int right;
    for (int iFace = 0; iFace < NumFaces; ++ iFace)
    {
        right = right_cell_of_face[iFace] + 1;
        if (right > NumElements || right < 0) right = 0;

        flowFile << right << " ";
        if ((iFace + 1) % nWordOfLine == 0) flowFile << "\n";
    }
    if (NumFaces % nWordOfLine != 0) flowFile << "\n";

    delete [] x;
    delete [] y;
    delete [] z;
    delete [] node_number_of_each_face;
    delete [] face2node;
    delete [] left_cell_of_face;
    delete [] right_cell_of_face;
    DelPointer2(qFiled);
}

void DumpToTecplotASCII::DumpFlowByCellTopo(fstream &flowFile, DataContainer *data, string zoneTitle)
{
    int NumPts, NumElements, NumFaces, TotalNumFaceNodes_Rect;

    PHRead(data, &NumPts, 1);
    PHRead(data, &NumElements, 1);
    PHRead(data, &NumFaces, 1);
    PHRead(data, &TotalNumFaceNodes_Rect, 1);

    int *node_number_of_each_cell = new int [NumElements];
    PHRead(data, node_number_of_each_cell, NumElements);

    int nTotalNode;
    PHRead(data, &nTotalNode, 1);

    RDouble *x = new RDouble [nTotalNode];
    RDouble *y = new RDouble [nTotalNode];
    RDouble *z = new RDouble [nTotalNode];
    PHRead(data, x, nTotalNode);
    PHRead(data, y, nTotalNode);
    PHRead(data, z, nTotalNode);

    RDouble **qFiled = NewPointer2<RDouble> (nVariable, nTotalNode);
    for (int iVariable = 0; iVariable < nVariable; ++ iVariable)
    {
        PHRead(data, qFiled[iVariable], nTotalNode);
    }

    int *cell2node = new int [TotalNumFaceNodes_Rect];
    PHRead(data, cell2node, TotalNumFaceNodes_Rect);

    flowFile << "zone T = " << zoneTitle << "\n";
    flowFile << "N = " << nTotalNode << "\n";
    flowFile << "E = " << NumElements << "\n";
    flowFile << "f = FEPOINT" << "\n";
    flowFile << "ET = quadrilateral" << "\n";

    int wordWidth = 20;
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        flowFile << setw(wordWidth) << x[iNode]
                 << setw(wordWidth) << y[iNode]
                 << setw(wordWidth) << z[iNode];
        for (int iVariable = 0; iVariable < nVariable; ++ iVariable)
        {
            flowFile << setw(wordWidth) << qFiled[iVariable][iNode];
        }
        flowFile << "\n";
    }

    int count = 0;
    for (int iCell = 0; iCell < NumElements; ++ iCell)
    {
        for (int iNode = 0; iNode < 4; ++ iNode)
        {
            flowFile << cell2node[count] + 1 << "  ";
            count ++;
        }
        flowFile << "\n";
    }

    delete [] x;    x = nullptr;
    delete [] y;    y = nullptr;
    delete [] z;    z = nullptr;
    delete [] node_number_of_each_cell;    node_number_of_each_cell = nullptr;
    delete [] cell2node;    cell2node = nullptr;
    DelPointer2(qFiled);
}

void DumpToTecplotASCII::DumpOrderedFlow(fstream &flowFile, DataContainer *data, string zoneTitle)
{
    int IMx, JMx, KMx;
    int nTotalNode;

    PHRead(data, &IMx, 1);
    PHRead(data, &JMx, 1);
    PHRead(data, &KMx, 1);
    PHRead(data, &nTotalNode, 1);

    flowFile << "zone T = " << zoneTitle << "\n";
    flowFile << "I = " << IMx << "\n";
    flowFile << "J = " << JMx << "\n";
    flowFile << "K = " << KMx << "\n";
    flowFile << "f = BLOCK" << "\n";

    RDouble *x = new RDouble [nTotalNode];
    RDouble *y = new RDouble [nTotalNode];
    RDouble *z = new RDouble [nTotalNode];
    PHRead(data, x, nTotalNode);
    PHRead(data, y, nTotalNode);
    PHRead(data, z, nTotalNode);

    RDouble **qFiled = NewPointer2<RDouble> (nVariable, nTotalNode);
    for (int iVariable = 0; iVariable < nVariable; ++ iVariable)
    {
        PHRead(data, qFiled[iVariable], nTotalNode);
    }

    int wordWidth = 20;
    int nWordOfLine = 5;
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        flowFile << setw(wordWidth) << x[iNode];
        if ((iNode + 1) % nWordOfLine == 0) 
        {
            flowFile << "\n";
        }
    }
    if (nTotalNode % nWordOfLine != 0) flowFile << "\n";

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        flowFile << setw(wordWidth) << y[iNode];
        if ((iNode + 1) % nWordOfLine == 0) 
        {
            flowFile << "\n";
        }
    }
    if (nTotalNode % nWordOfLine != 0) flowFile << "\n";

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        flowFile << setw(wordWidth) << z[iNode];
        if ((iNode + 1) % nWordOfLine == 0) 
        {
            flowFile << "\n";
        }
    }
    if (nTotalNode % nWordOfLine != 0) flowFile << "\n";

    for (int iVariable = 0; iVariable < nVariable; ++ iVariable)
    {
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            flowFile << setw(wordWidth) << qFiled[iVariable][iNode];
            if ((iNode + 1) % nWordOfLine == 0)
            {
                flowFile << "\n";
            }
        }
        if (nTotalNode % nWordOfLine != 0) flowFile << "\n";
    }

    delete [] x;
    delete [] y;
    delete [] z;
    DelPointer2(qFiled);
}

void write_string(fstream &file, string &cs)
{
    int nlen = static_cast<int>(cs.length());
    PHWrite(file, nlen);

    char *data = new char [nlen+1];

    cs.copy(data, nlen);
    data[nlen] = '\0';

    PHWrite(file, data, nlen+1);

    delete [] data;
}

void read_string(fstream &file, string &cs)
{
    int nlen = 0;
    PHRead(file, nlen);

    char *data = new char[nlen+1];

    PHRead(file, data, nlen+1);
    
    cs = data;

    delete [] data;
}

PHCutPlane::PHCutPlane()
{
    title  = 0;
    nTitle = 0;
    neqn   = 0;
    data   = 0;

    x = 0;
    y = 0;
    z = 0;

    faceIndex[0] = 0;
    faceIndex[1] = 0;
    faceIndex[2] = 0;
    faceIndex[3] = 0;
    nTotalNode   = 0;
    nTotalFace   = 0;
}

PHCutPlane::~PHCutPlane()
{
    for (int i = 0; i < nTitle; ++ i)
    {
        delete [] title[i];
    }
    delete [] title;

    for (int i = 0; i < neqn; ++ i)
    {
        delete [] data[i];
    }
    delete [] data;

    delete [] x;
    delete [] y;
    delete [] z;

    delete [] faceIndex[0];
    delete [] faceIndex[1];
    delete [] faceIndex[2];
    delete [] faceIndex[3];
}

void PHCutPlane::SortData(RDouble *value)
{
    using namespace PHSPACE;
    vector<DataStruct_Sort< RDouble> > value_list;
    DataStruct_Sort<RDouble> var;
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        var.value = value[iNode];
        var.index = iNode;
        value_list.push_back(var);
    }
    sort(value_list.begin(), value_list.end());

    vector<int> index_list;
    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        index_list.push_back(value_list[iNode].index);
    }

    ReorderList(x, index_list);
    ReorderList(y, index_list);
    ReorderList(z, index_list);

    for (int m = 0; m < neqn; ++ m)
    {
        ReorderList(data[m], index_list);
    }
}

void CutPlane(UnstructGrid *grid, RDouble **node_var, int nvar, RDouble cut, PHCutPlane *cutpl, int ndir)
{
    int nTotalFace = grid->GetNTotalFace();

    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();
    int *face2node = grid->GetFace2Node();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble *xyz = 0;

    if (ndir == X_DIR)
    {
        xyz = x;
    }
    else if (ndir == Y_DIR)
    {
        xyz = y;
    }
    else
    {
        xyz = z;
    }

    RDouble *pmin = grid->GetMinBox();
    RDouble *pmax = grid->GetMaxBox();

    RDouble mindis, maxdis;
    grid->GetMinMaxDS(mindis, maxdis);

    RDouble min_tree[3], max_tree[3];

    min_tree[0] = pmin[0] - mindis;
    min_tree[1] = pmin[1] - mindis;
    min_tree[2] = pmin[2] - mindis;
                  
    max_tree[0] = pmax[0] + mindis;
    max_tree[1] = pmax[1] + mindis;
    max_tree[2] = pmax[2] + mindis;

    //RDouble tol = mindis / 5;    //! This method may cause problems when slicing.

    RDouble tol = 1.0e-8;

    typedef DataStruct_AdtTree<int, RDouble> AdtTree;
    AdtTree coor_tree(3, min_tree, max_tree);

    RDouble minwin[3], maxwin[3], coor[3];

    RDouble eps = half * tol;

    vector<vector<RDouble> > var_list;
    var_list.resize(0);

    vector<RDouble> x_list, y_list, z_list;
    vector<RDouble> var_l(nvar), var_r(nvar), var(nvar);

    int nodepos = 0;
    int pcount  = 0;
    int index;

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        RDouble coor_min =   LARGE;
        RDouble coor_max = - LARGE;

        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            int index1 = face2node[nodepos + j];
            coor_min = MIN(xyz[index1], coor_min);
            coor_max = MAX(xyz[index1], coor_max);
        }

        if (cut < coor_min - eps || cut > coor_max + eps)
        {
            nodepos += node_number_of_each_face[iFace];
            continue;
        }

        RDouble xl, yl, zl, xr, yr, zr, xm, ym, zm;

        int p1, p2;
        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            p1 = face2node[nodepos + j];
            p2 = face2node[nodepos + (j + 1) % node_number_of_each_face[iFace]];

            for (int m = 0; m < nvar; ++ m)
            {
                var_l[m] = node_var[m][p1];
                var_r[m] = node_var[m][p2];
            }

            xl = x[p1];
            yl = y[p1];
            zl = z[p1];

            xr = x[p2];
            yr = y[p2];
            zr = z[p2];

            coor_min = MIN(xyz[p1], xyz[p2]);
            coor_max = MAX(xyz[p1], xyz[p2]);

            if (cut < coor_min - eps || cut > coor_max + eps)
            {
                continue;
            }

            if (ABS(cut - xyz[p1]) <= eps && ABS(cut - xyz[p2]) <= eps)
            {
                coor[0] = xl;
                coor[1] = yl;
                coor[2] = zl;

                minwin[0] = coor[0] - tol;
                minwin[1] = coor[1] - tol;
                minwin[2] = coor[2] - tol;

                maxwin[0] = coor[0] + tol;
                maxwin[1] = coor[1] + tol;
                maxwin[2] = coor[2] + tol;

                if (GetCoorIndex(&coor_tree, coor, minwin, maxwin, pcount, index) == 0)
                {
                    var_list.push_back(var_l);
                    x_list.push_back(xl);
                    y_list.push_back(yl);
                    z_list.push_back(zl);
                }

                coor[0] = xr;
                coor[1] = yr;
                coor[2] = zr;

                minwin[0] = coor[0] - tol;
                minwin[1] = coor[1] - tol;
                minwin[2] = coor[2] - tol;

                maxwin[0] = coor[0] + tol;
                maxwin[1] = coor[1] + tol;
                maxwin[2] = coor[2] + tol;

                if (GetCoorIndex(&coor_tree, coor, minwin, maxwin, pcount, index) == 0)
                {
                    var_list.push_back(var_r);
                    x_list.push_back(xr);
                    y_list.push_back(yr);
                    z_list.push_back(zr);
                }
            }
            else if (ABS(cut - xyz[p1]) <= eps && ABS(cut - xyz[p2]) > eps)
            {
                coor[0] = xl;
                coor[1] = yl;
                coor[2] = zl;

                minwin[0] = coor[0] - tol;
                minwin[1] = coor[1] - tol;
                minwin[2] = coor[2] - tol;

                maxwin[0] = coor[0] + tol;
                maxwin[1] = coor[1] + tol;
                maxwin[2] = coor[2] + tol;

                if (GetCoorIndex(&coor_tree, coor, minwin, maxwin, pcount, index) == 0)
                {
                    var_list.push_back(var_l);
                    x_list.push_back(xl);
                    y_list.push_back(yl);
                    z_list.push_back(zl);
                }
            }
            else if (ABS(cut - xyz[p1]) > eps && ABS(cut - xyz[p2]) <= eps)
            {
                coor[0] = xr;
                coor[1] = yr;
                coor[2] = zr;

                minwin[0] = coor[0] - tol;
                minwin[1] = coor[1] - tol;
                minwin[2] = coor[2] - tol;

                maxwin[0] = coor[0] + tol;
                maxwin[1] = coor[1] + tol;
                maxwin[2] = coor[2] + tol;

                if (GetCoorIndex(&coor_tree, coor, minwin, maxwin, pcount, index) == 0)
                {
                    var_list.push_back(var_r);
                    x_list.push_back(xr);
                    y_list.push_back(yr);
                    z_list.push_back(zr);
                }
            }
            else
            {
                RDouble cl, cr;
                cr = (cut - xyz[p1]) / (xyz[p2] - xyz[p1]);
                cl = 1 - cr;

                for (int m = 0; m < nvar; ++ m)
                {
                    var[m] = cl * var_l[m] + cr * var_r[m];
                }

                xm = cl * xl + cr * xr;
                ym = cl * yl + cr * yr;
                zm = cl * zl + cr * zr;

                coor[0] = xm;
                coor[1] = ym;
                coor[2] = zm;

                minwin[0] = coor[0] - tol;
                minwin[1] = coor[1] - tol;
                minwin[2] = coor[2] - tol;

                maxwin[0] = coor[0] + tol;
                maxwin[1] = coor[1] + tol;
                maxwin[2] = coor[2] + tol;

                if (GetCoorIndex(&coor_tree, coor, minwin, maxwin, pcount, index) == 0)
                {
                    var_list.push_back(var);
                    x_list.push_back(xm);
                    y_list.push_back(ym);
                    z_list.push_back(zm);
                }
            }
        }

        nodepos += node_number_of_each_face[iFace];
    }

    cutpl->neqn = nvar;
    cutpl->nTotalNode = static_cast<int>(var_list.size());
    cutpl->data = new RDouble * [nvar];
    for (int m = 0; m < nvar; ++ m)
    {
        cutpl->data[m] = new RDouble [cutpl->nTotalNode];
    }

    cutpl->x = new RDouble [cutpl->nTotalNode];
    cutpl->y = new RDouble [cutpl->nTotalNode];
    cutpl->z = new RDouble [cutpl->nTotalNode];

    for (int iNode = 0; iNode < cutpl->nTotalNode; ++ iNode)
    {
        cutpl->x[iNode] = x_list[iNode];
        cutpl->y[iNode] = y_list[iNode];
        cutpl->z[iNode] = z_list[iNode];
    }

    for (int iNode = 0; iNode < cutpl->nTotalNode; ++ iNode)
    {
        for (int m = 0; m < nvar; ++ m)
        {
            cutpl->data[m][iNode] = var_list[iNode][m];
        }
    }
}

void CutPlaneCellCenter(UnstructGrid *grid, RDouble **q_var, int nvar, RDouble cut, PHCutPlane *cutpl, int ndir)
{
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();

    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();
    int *face2node = grid->GetFace2Node();
    int *left_cell_of_face = grid->GetLeftCellOfFace();
    int *right_cell_of_face = grid->GetRightCellOfFace();
    bool*cut_cell = new bool [nTotalCell];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        cut_cell[iCell] = false;
    }

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *xyz = 0;

    if (ndir == X_DIR)
    {
        xyz = x;
    }
    else if (ndir == Y_DIR)
    {
        xyz = y;
    }
    else
    {
        xyz = z;
    }

    int nodepos = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int cl = left_cell_of_face[iFace];
        int cr = right_cell_of_face[iFace];

        RDouble coor_min =   LARGE;
        RDouble coor_max = - LARGE;

        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            int index = face2node[nodepos ++];
            coor_min = MIN(xyz[index], coor_min);
            coor_max = MAX(xyz[index], coor_max);
        }

        if ((coor_min - cut) * (coor_max - cut) > 0.0)
        {
            continue;
        }

        cut_cell[cl] = true;
        if (iFace >= nBoundFace) cut_cell[cr] = true;
    }

    vector<vector<RDouble> > var_list;
    vector<RDouble> x_list, y_list, z_list;
    vector<RDouble> var(nvar);

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (cut_cell[iCell])
        {
            for (int m = 0; m < nvar; ++ m)
            {
                var[m] = q_var[m][iCell];
            }

            var_list.push_back(var);
            x_list.push_back(xcc[iCell]);
            y_list.push_back(ycc[iCell]);
            z_list.push_back(zcc[iCell]);
        }
    }

    cutpl->neqn = nvar;
    cutpl->nTotalNode = static_cast<int>(var_list.size());
    cutpl->data = new RDouble * [nvar];
    for (int m = 0; m < nvar; ++ m)
    {
        cutpl->data[m] = new RDouble [cutpl->nTotalNode];
    }

    cutpl->x = new RDouble [cutpl->nTotalNode];
    cutpl->y = new RDouble [cutpl->nTotalNode];
    cutpl->z = new RDouble [cutpl->nTotalNode];

    for (int iNode = 0; iNode < cutpl->nTotalNode; ++ iNode)
    {
        cutpl->x[iNode] = x_list[iNode];
        cutpl->y[iNode] = y_list[iNode];
        cutpl->z[iNode] = z_list[iNode];
    }

    for (int iNode = 0; iNode < cutpl->nTotalNode; ++ iNode)
    {
        for (int m = 0; m < nvar; ++ m)
        {
            cutpl->data[m][iNode] = var_list[iNode][m];
        }
    }

    delete [] cut_cell;
}

LIB_EXPORT Post_VisualWall::Post_VisualWall()
{
    nVisualVariables = 0;
    visualVariables = NULL;
}

LIB_EXPORT Post_VisualWall::Post_VisualWall(int nVisualNumber, int *visualVariablesType)
{
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_CP, "cp"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_CF, "cf"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_YPLUS, "yplus"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_QNONDIM, "Q_NonDim"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_QDIM, "Q_Dim(kW/m2)"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_PW, "pw(Pa)"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_TW, "Tw(K)"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_RHOW, "rhow(kg/m3)"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_VX, "Vx(m/s)"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_VY, "Vy(m/s)"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_VZ, "Vz(m/s)"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_VS, "Vs(m/s)"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_ST, "St"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_CH, "Ch"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_RE, "Re_w"));
    visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_KN, "Kn_wall"));

    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nChemical > 0)
    {
        visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_QT, "Qt(kW/m2)"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_QS, "Qs(kW/m2)"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_QV, "Qv(kW/m2)"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_QE, "Qe(kW/m2)"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_WALL_NS, "Ns"));
    }

    int nSlipBCModel = GlobalDataBase::GetIntParaFromDB("nSlipBCModel");
    if (nSlipBCModel > 0)
    {
        visualVariablesMap.insert(pair<int, string>(VISUAL_SLIP_TS, "Tts(K)"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_SLIP_TV, "Tvs(K)"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_SLIP_TE, "Tes(K)"));
        visualVariablesMap.insert(pair<int, string>(VISUAL_SLIP_DTS, "TtsJump(K)"));
    }

    nVisualVariables = 0;
    int tempVariablesType[50];
    for (int iVar = 0; iVar < nVisualNumber; ++ iVar)
    {
        int varTemp = visualVariablesType[iVar];
        if (visualVariablesMap.find(varTemp) == visualVariablesMap.end())
        {
            continue;
        }
        tempVariablesType[nVisualVariables] = varTemp;
        ++ nVisualVariables;
    }

    GlobalDataBase::UpdateData("nVisualWallVariables", &nVisualVariables, PHINT, 1);
    visualVariables = NULL;
    if (nVisualVariables > 0)
    {
        visualVariables = new int[nVisualVariables];
        for (int iVar = 0; iVar < nVisualVariables; ++ iVar)
        {
            visualVariables[iVar] = tempVariablesType[iVar];
        }
        GlobalDataBase::UpdateData("visualWallVariables", visualVariables, PHINT, nVisualVariables);
    }

}

LIB_EXPORT Post_VisualWall::~Post_VisualWall()
{
    if (nVisualVariables > 0)
    {
        delete [] visualVariables;    visualVariables = nullptr;
        nVisualVariables = 0;
    }
}

}
