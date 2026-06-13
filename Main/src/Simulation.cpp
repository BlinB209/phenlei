#include "Simulation.h"
#include "IO_FileName.h"
#include "TK_Parse.h"
#include "GridFactory.h"
#include "AleManager.h"
#include "GridPartition.h"
#include "Controller.h"
#include "TK_Time.h"
#include "PHIO.h"
#include "Geo_Sample.h"
#include "Geo_SimpleVC.h"
#include "UnstructuredOversetConfig.h"
#include "HOSimulation.h"
#include "Post_ProbesWrite.h"

#ifdef USE_CUDA
#include "BasicDeviceVariables.h"
#include "GPUFaceColor.h"
#endif 


#include "OversetInformation.h"
#include "SolversList.h"
#include "PostProcess.h"
#include "Glb_Dimension.h"
#include "Zone.h"
#include "Pre_WalldistCompute.h"
#include "Post_WriteVisualFile.h"
#include "LESSolverStruct.h"
#include "LESSolverUnstr.h"
#ifdef USE_INCOMSOLVER
#include "IncomUnstructSolver.h"
#endif
#ifdef USE_SpecDiffHybSolver
#include "SpecDiffHybSolver.h"
#include "SpecSolver.h"
#endif
#ifdef USE_DEMOSOLVER
#include "DemoSolver.h"
#endif
#ifdef USE_LBMSolverMPI
#include "LBMSolverMPI.hpp"
#endif
#ifdef USE_LBMSolverOMP
#include "LBMSolverOMP.hpp"
#endif
using namespace std;

namespace PHSPACE
{
Simulation::Simulation(int dimension)
{
    using namespace PHMPI;

    SetDim(dimension);

    region = new Region();
}

Simulation::~Simulation()
{
    using namespace PHMPI;

    delete region;

    TK_Exit::ExitPHengLEI();
}

void Simulation::Start()
{
    int taskType = GetTaskCode();

    switch (taskType)
    {
        case SOLVE_FIELD:
            RunCFD();
            break;
        case CREATE_GRID:
            PHSPACE::CreateGrid();
            break;
        case CAL_WALL_DIST:
            ComputeWallDistance();
            break;
        case PARTITION_GRID:
            IsNeedConvert();
            PHSPACE::PartitionGrid();
            break;
        case OVERSETGRID_VIEW:
            PHSPACE::OversetGridView();
            break;
        case OVERSET_CONFIG:
            PHSPACE::RunOversetGridConfig();
            break;
        case HO_SOLVER:
            {
                HOUnstruct::HOSimulation *uho = new HOUnstruct::HOSimulation();
                uho->Run();
                delete uho;
                break;
            }
        case INTEGRATIVE_SOLVER:
            RunIntegrativeSolver();
            break;
        #ifdef USE_SpecDiffHybSolver
        case SPECDIFFHYB_SOLVER:
            RunSpecDiffHybSolver();
            break;
        case SPEC_SOLVER:
            RunSpecSolver();
            break;
        #endif
        case POST_PROCESSING:
            RunPostProcessing();
            break;
        #ifdef USE_LBMSolverMPI
        case LBM_SOLVER_MPI:
            RunLBMSolverMPI();
            break;
        #endif
        #ifdef USE_LBMSolverOMP
        case LBM_SOLVER_OMP:
            RunLBMSolverOMP();
            break;
        #endif
        default:
            TK_Exit::UnexpectedVarValue("taskType = ", taskType);
            break;
    }
}

void Simulation::IsNeedConvert()
{
    string originalGridFile = GlobalDataBase::GetStrParaFromDB("original_grid_file");

    string fileMame, fileSuffix;
    GetNameExt(originalGridFile, fileMame, fileSuffix, ".");

    set < pair<string, int> > gridTypeMap;
    gridTypeMap.insert(pair <string, int> ("fts",  PHENGLEI_TYPE ));
    gridTypeMap.insert(pair <string, int> ("cgns", CGNS_TYPE     ));
    gridTypeMap.insert(pair <string, int> ("grd",  PLOT3D_TYPE   ));
    gridTypeMap.insert(pair <string, int> ("uns",  FIELDVIEW_TYPE));
    gridTypeMap.insert(pair <string, int> ("cas",  FLUENT_TYPE   ));

    int sourceGridType = 0;
    set< pair<string, int> >::iterator iter;
    for (iter = gridTypeMap.begin(); iter != gridTypeMap.end(); ++ iter)
    {
        if (fileSuffix.compare((*iter).first) != 0)
        {
            continue;
        }

        sourceGridType = (*iter).second;
        break;
    }

    if (sourceGridType == 0)
    {
        ostringstream oss;
        oss << "  Error: this situation has not been considered, for filetype = ." << fileSuffix << endl;
        TK_Exit::ExceptionExit(oss.str());
    }

    if (sourceGridType != PHENGLEI_TYPE)
    {
        GlobalDataBase::UpdateData("from_gtype", &sourceGridType, PHINT, 1);
        GlobalDataBase::UpdateData("from_gfile", &originalGridFile, PHSTRING, 1);

        int targetGridType = 1;
        GlobalDataBase::UpdateData("to_gtype", &targetGridType, PHINT, 1);

        string outputGridFile = ChangeExtensionOfFileName(originalGridFile, "fts");
        GlobalDataBase::UpdateData("out_gfile", &outputGridFile, PHSTRING, 1);

        int partitionedGridType = GlobalDataBase::GetIntParaFromDB("pgridtype");
        GlobalDataBase::UpdateData("gridtype", &partitionedGridType, PHINT, 1);

        int gridObjection = 1;
        GlobalDataBase::UpdateData("gridobj", &gridObjection, PHINT, 1);

        int numberOfSimulationTask = PHSPACE::CREATE_GRID;
        GlobalDataBase::UpdateData("nsimutask", &numberOfSimulationTask, PHINT, 1);

        int multiBlock = 1;
        GlobalDataBase::UpdateData("multiblock", &multiBlock, PHINT, 1);

        int axisUp = 1;
        GlobalDataBase::UpdateData("axisup", &axisUp, PHINT, 1);

        int iadapt = 0;
        GlobalDataBase::UpdateData("iadapt", &iadapt, PHINT, 1);

        CompleteParamInfo();

        PHSPACE::CreateGrid();

        GlobalDataBase::UpdateData("original_grid_file", &outputGridFile, PHSTRING, 1);

        numberOfSimulationTask = PHSPACE::PARTITION_GRID;
        GlobalDataBase::UpdateData("nsimutask", &numberOfSimulationTask, PHINT, 1);
    }
}

void Simulation::RunCFD()
{
    double last_time[3] = {0.0};
    TIME_SPACE::StartTime(last_time);

    ReadGeometry();

    UpdateAllZonesPara();

    InitGeometry();

    SwapGeometry();

    InitOversetGrid();

    if (PHMPI::CurrentProcessorIsGridProcessor())
    {
        return;
    }

    //! The flag of solver, the adaptive method sets value great than zero.
    int isAdaptiveSolver = GlobalDataBase::GetIntParaFromDB("isAdaptiveSolver");
    //! Solving the NS equations using the automated process from the perfect gas to the thermochemical nonequilibrium gas.
    if (isAdaptiveSolver > 0)
    {
        GenerateWalldist();
        //! Run CFD using the adaptive method.
        SelfAdaptionRunCFD();
        return;
    }

    //! Define solver names for all zones and create the new defined solvers.
    //! IMPORTANT for solver developers:
    //! Define your solver name which is same with the solver class name here!
    DefineAndCreateSolvers();

    InitGlobalBoundaryCondition();

    InitGlobalVolumeCondition();

    GenerateWalldist();

#ifdef USE_CUDA
    using namespace GPUMemory;
    using namespace GPUGeomVariables;
    GPUGeomInfoAllocCopy();
#endif

    //! Initialize each solver on each zone.
    InitializeSolvers();

    MultiGridInitFlow();

#ifdef USE_CUDA
   using namespace GPUMemory;
   using namespace GPUGeomVariables;
   using namespace GPUFlowVariables;
   using namespace GPUFaceColor;
   GPUFaceColorMain();
#endif

    //! The main simulation entrance.
    SolveSimulation();

    //! Do some work after simulation.
    PostSimulation();

#ifdef USE_CUDA
    using namespace GPUMemory;
    GPUMemoryFree();
#endif
}

bool Simulation::IsRestartProcess()
{
    using namespace PHMPI;

    string restartFile = ".\results\flow.dat";
    GlobalDataBase::GetData("restartNSFile", &restartFile, PHSTRING, 1);
    if (PHMPI::IsParallelRun())
    {
        restartFile = PHSPACE::AddSymbolToFileName(restartFile, "_", 0);
    }

    ifstream infile(restartFile.c_str(), ios::in);
    bool isRestart = false;
    if (infile)
    {
        isRestart = true;
        infile.close();
        infile.clear();
    }

    if (isRestart)
    {
        hid_t file;
        file = OpenHDF5File(restartFile);
        int gridID = 0;
        hid_t grploc;
        string grpName;
        ostringstream oss;
        oss << "Group" << gridID;
        grpName = oss.str();

        grploc = OpenGroup(file, grpName);

        int preChemical = 0, preTemperatureModel = 1;
        if (CheckDataExist(grploc, "nIsChemicalFlow"))
        {
            ReadData(grploc, &preChemical, "nIsChemicalFlow");
        }
        if (CheckDataExist(grploc, "nTemperatureModel"))
        {
            ReadData(grploc, &preTemperatureModel, "nTemperatureModel");
        }

        GlobalDataBase::UpdateData("previousChemical", &preChemical, PHINT, 1);
        GlobalDataBase::UpdateData("previousTModel", &preTemperatureModel, PHINT, 1);

        H5Gclose(grploc);
        H5Fclose(file);
    }

    return isRestart;
}

void Simulation::SelfAdaptionRunCFD()
{
    int nFinalChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    GlobalDataBase::UpdateData("nFinalChemical", &nFinalChemical, PHINT, 1);
    int nFinalTModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    GlobalDataBase::UpdateData("nFinalTModel", &nFinalTModel, PHINT, 1);
    string nFinalNameList = GlobalDataBase::GetStrParaFromDB("speciesName");
    string nFinalMassFraction = GlobalDataBase::GetStrParaFromDB("initMassFraction");
    GlobalDataBase::UpdateData("nFinalNameList", &nFinalNameList, PHSTRING, 1);
    GlobalDataBase::UpdateData("nFinalMassFraction", &nFinalMassFraction, PHSTRING, 1);

    if (IsRestartProcess())
    {
        int preChemical = GlobalDataBase::GetIntParaFromDB("previousChemical");
        int preTModel = GlobalDataBase::GetIntParaFromDB("previousTModel");

        if (preChemical == 1 && preTModel == 1)
        {
            //! Restart from the one-temperature model.
            SecondStepRunCFD(true);
            if (nFinalTModel == 1)
            {
                //! The third step: the final gas model.
                FinalStepRunCFD();
            }
            else
            {
                //! The second step: two-temperature model.
                ThirdStepRunCFD();

                //! The third step: the final gas model.
                FinalStepRunCFD();
            }

            return;
        }
        else if (preChemical == 1 && preTModel == 2)
        {
            //! Restart from the two-temperature model.
            ThirdStepRunCFD(true);

            //! The third step: the final gas model.
            FinalStepRunCFD();

            return;
        }
    }

    //! The first step: the perfect gas.
    FirstStepRunCFD();
    if (nFinalChemical == 0)
    {
        return;
    }
    //! The second step: one-temperature model.
    SecondStepRunCFD();
    if (nFinalTModel == 1)
    {
        //! The third step: the final gas model.
        FinalStepRunCFD();
    }
    else
    {
        //! The second step: two-temperature model.
        ThirdStepRunCFD();

        //! The third step: the final gas model.
        FinalStepRunCFD();
    }
}

void Simulation::FirstStepRunCFD(bool isStart/* = true*/)
{
    //! Record the serial number of current iteration step.
    int nCurrentStep = 1;
    GlobalDataBase::UpdateData("nCurrentStep", &nCurrentStep, PHINT, 1);

    int nCurrentProcessorID = PHMPI::GetCurrentProcessorID();
    if (nCurrentProcessorID == PHMPI::server)
    {
        cout << endl;
        cout << "nCurrentStep=" << nCurrentStep << endl;
    }

    //! The first step: run CFD using the perfect gas model.
    int nCurChemical = 0;
    GlobalDataBase::UpdateData("nchem", &nCurChemical, PHINT, 1);
    int nCurTModel = 1;
    GlobalDataBase::UpdateData("ntmodel", &nCurTModel, PHINT, 1);

    //! The index of the variable pressure.
    int nKeyVariableIndex = 4;
    GlobalDataBase::UpdateData("nKeyVariableIndex", &nKeyVariableIndex, PHINT, 1);
    string keyVariable = "pressure";
    GlobalDataBase::UpdateData("keyVariable", &keyVariable, PHSTRING, 1);

    DefineAndCreateSolvers();

    InitGlobalBoundaryCondition();

    GenerateWalldist();

    InitializeSolvers();

    MultiGridInitFlow();

    SelfAdaptionSolve();

    PostSimulation();

    //! Release the dynamic memories at first.
    CleanGlobalValuesSolver();

    RefreshSolvers();
}

void Simulation::SecondStepRunCFD(bool isStart/* = false*/)
{
    //! Record the serial number of current iteration step.
    int nCurrentStep = 2;
    GlobalDataBase::UpdateData("nCurrentStep", &nCurrentStep, PHINT, 1);

    int nCurrentProcessorID = PHMPI::GetCurrentProcessorID();
    if (nCurrentProcessorID == PHMPI::server)
    {
        cout << endl;
        cout << "nCurrentStep=" << nCurrentStep << endl;
    }

    //! The second step: run CFD using the one-temperature model.
    int nCurChemical = 1;
    GlobalDataBase::UpdateData("nchem", &nCurChemical, PHINT, 1);
    int nCurTModel = 1;
    GlobalDataBase::UpdateData("ntmodel", &nCurTModel, PHINT, 1);
    string speciesName = "O, O2, NO, N, N2";
    GlobalDataBase::UpdateData("speciesName", &speciesName, PHSTRING, 1);
    string initMassFraction = "0.0, 0.233, 0.0, 0.0, 0.767";
    GlobalDataBase::UpdateData("initMassFraction", &initMassFraction, PHSTRING, 1);

    int nControlVariable = 1;
    if (GlobalDataBase::IsExist("nControlVariable", PHINT, 1))
    {
        nControlVariable = GlobalDataBase::GetIntParaFromDB("nControlVariable");
    }
    int nKeyVariableIndex = nControlVariable;
    if (nControlVariable == 5)
    {
        nKeyVariableIndex = 6;
    }
    else if (nControlVariable == 6)
    {
        nKeyVariableIndex = 9;
    }
    else if (nControlVariable > 6)  //! default.
    {
        nKeyVariableIndex = 1;      //! The index of the variable temperature.
    }
    GlobalDataBase::UpdateData("nKeyVariableIndex", &nKeyVariableIndex, PHINT, 1);

    string keyVariable = "";
    switch (nControlVariable)
    {
        case 0:
            keyVariable = "density";
        case 1:
            keyVariable = "translation temperature";
        case 2:
            keyVariable = "vibration temperature";
        case 3:
            keyVariable = "electron temperature";
        case 4:
            keyVariable = "pressure";
        default:
            keyVariable = "translation temperature";
            break;
    }
    GlobalDataBase::UpdateData("keyVariable", &keyVariable, PHSTRING, 1);

    if (!isStart)
    {
        InitGlobalValuesSolver();
    }
    else
    {
        DefineAndCreateSolvers();

        GenerateWalldist();
    }

    InitGlobalBoundaryCondition();

    InitializeSolvers();

    MultiGridInitFlow();

    SelfAdaptionSolve();

    PostSimulation();

    //! Release the dynamic memories at first.
    CleanGlobalValuesSolver();

    RefreshSolvers();
}

void Simulation::ThirdStepRunCFD(bool isStart/* = false*/)
{
    //! Record the serial number of current iteration step.
    int nCurrentStep = 3;
    GlobalDataBase::UpdateData("nCurrentStep", &nCurrentStep, PHINT, 1);

    int nCurrentProcessorID = PHMPI::GetCurrentProcessorID();
    if (nCurrentProcessorID == PHMPI::server)
    {
        cout << endl;
        cout << "nCurrentStep=" << nCurrentStep << endl;
    }

    //! The third step: run CFD using the two-temperature model.
    int nCurChemical = 1;
    GlobalDataBase::UpdateData("nchem", &nCurChemical, PHINT, 1);
    int nCurTModel = 2;
    GlobalDataBase::UpdateData("ntmodel", &nCurTModel, PHINT, 1);
    int nEnergyAssembly = 0;
    if (GlobalDataBase::IsExist("nEnergyAssembly", PHINT, 1))
    {
        nEnergyAssembly = GlobalDataBase::GetIntParaFromDB("nEnergyAssembly");
    }
    if (nEnergyAssembly == 1)
    {
        //! The vibration energy is computed using the polynomial fitting method.
        int nTEnergyModel = 1;
        GlobalDataBase::UpdateData("nTEnergyModel", &nTEnergyModel, PHINT, 1);
    }

    string speciesName = "O, O2, NO, N, N2";
    GlobalDataBase::UpdateData("speciesName", &speciesName, PHSTRING, 1);
    string initMassFraction = "0.0, 0.233, 0.0, 0.0, 0.767";
    GlobalDataBase::UpdateData("initMassFraction", &initMassFraction, PHSTRING, 1);

    //! The index of the variable vibrational temperature.
    int nKeyVariableIndex = 2;
    GlobalDataBase::UpdateData("nKeyVariableIndex", &nKeyVariableIndex, PHINT, 1);
    string keyVariable = "vibration temperature";
    GlobalDataBase::UpdateData("keyVariable", &keyVariable, PHSTRING, 1);

    if (!isStart)
    {
        InitGlobalValuesSolver();
    }
    else
    {
        DefineAndCreateSolvers();

        GenerateWalldist();
    }

    InitGlobalBoundaryCondition();

    InitializeSolvers();

    MultiGridInitFlow();

    SelfAdaptionSolve();

    PostSimulation();

    //! Release the dynamic memories at first.
    CleanGlobalValuesSolver();

    RefreshSolvers();
}

void Simulation::FinalStepRunCFD()
{
    //! Record the serial number of current iteration step.
    int nCurrentStep = 4;
    GlobalDataBase::UpdateData("nCurrentStep", &nCurrentStep, PHINT, 1);

    int nCurrentProcessorID = PHMPI::GetCurrentProcessorID();
    if (nCurrentProcessorID == PHMPI::server)
    {
        cout << endl;
        cout << "nCurrentStep=" << nCurrentStep << endl;
    }

    //! The third step: run CFD using the final gas model.
    int nCurChemical = GlobalDataBase::GetIntParaFromDB("nFinalChemical");
    GlobalDataBase::UpdateData("nchem", &nCurChemical, PHINT, 1);
    int nCurTModel = GlobalDataBase::GetIntParaFromDB("nFinalTModel");
    GlobalDataBase::UpdateData("ntmodel", &nCurTModel, PHINT, 1);
    string speciesName = GlobalDataBase::GetStrParaFromDB("nFinalNameList");
    GlobalDataBase::UpdateData("speciesName", &speciesName, PHSTRING, 1);
    string initMassFraction = GlobalDataBase::GetStrParaFromDB("nFinalMassFraction");
    GlobalDataBase::UpdateData("initMassFraction", &initMassFraction, PHSTRING, 1);
    int nEnergyAssembly = 0;
    if (GlobalDataBase::IsExist("nEnergyAssembly", PHINT, 1))
    {
        nEnergyAssembly = GlobalDataBase::GetIntParaFromDB("nEnergyAssembly");
    }
    if (nEnergyAssembly == 1)
    {
        //! The vibration energy is computed using the default method.
        int nTEnergyModel = 0;
        GlobalDataBase::UpdateData("nTEnergyModel", &nTEnergyModel, PHINT, 1);
    }

    //! Default: the index of the variable vibrational temperature.
    int nKeyVariableIndex = 2;
    string keyVariable = "vibration temperature";
    if (nCurTModel == 1)
    {
        int nControlVariable = 1;
        if (GlobalDataBase::IsExist("nControlVariable",PHINT, 1))
        {
            nControlVariable = GlobalDataBase::GetIntParaFromDB("nControlVariable");
        }

        switch (nControlVariable)
        {
            case 0:
                keyVariable = "density";
            case 1:
                keyVariable = "translation temperature";
            case 2:
                keyVariable = "vibration temperature";
            case 3:
                keyVariable = "electron temperature";
            case 4:
                keyVariable = "pressure";
            default:
                keyVariable = "translation temperature";
                break;
        }

        nKeyVariableIndex = nControlVariable;   //! The index of the variable temperature.
        if (nControlVariable == 5)
        {
            nKeyVariableIndex = 6;
        }
        else if (nControlVariable == 6)
        {
            nKeyVariableIndex = 9;
        }
        else if (nControlVariable > 6)
        {
            nKeyVariableIndex = 1;
        }
    }
    GlobalDataBase::UpdateData("nKeyVariableIndex", &nKeyVariableIndex, PHINT, 1);
    GlobalDataBase::UpdateData("keyVariable", &keyVariable, PHSTRING, 1);

    InitGlobalValuesSolver();

    InitGlobalBoundaryCondition();

    InitializeSolvers();

    MultiGridInitFlow();

    SelfAdaptionSolve();

    PostSimulation();
}

void Simulation::RunIntegrativeSolver()
{
    int numberOfProcessor = PHMPI::GetNumberOfProcessor();
    int currentProcessor  = PHMPI::GetCurrentProcessorID();
    int serverProcessor   = PHMPI::GetServerProcessorID();
    int tasktype = 0;

    string fromGridName = GlobalDataBase::GetStrParaFromDB("from_gfile");
    string outGridName = ChangeExtensionOfFileName(fromGridName, "fts");
    string partition_grid_file = AddSymbolToFileName(outGridName, "__", numberOfProcessor);
    GlobalDataBase::UpdateData("out_gfile", &outGridName, PHSTRING, 1);

    bool partitionFtsExist = false;
    if (numberOfProcessor > 1)
    {
        if (SerialFileExist(partition_grid_file))
        {
            partitionFtsExist = true;
        }
    }

    //! fts file not exist, need convert grid.
    if (!SerialFileExist(outGridName))
    {
        if (!partitionFtsExist)
        {
            WriteLogFile("Begin convert grid !\n");
            PrintToWindow("Begin convert grid !\n\n");

            tasktype = PHSPACE::CREATE_GRID;
            GlobalDataBase::UpdateData("nsimutask", &tasktype, PHINT, 1);

            int isConvertOver = 0;
            PHMPI::SetNumberOfProcessor(1);

            if (currentProcessor == serverProcessor)
            {
                CreateGrid();
                isConvertOver = 1;
            }

            PHMPI::SetNumberOfProcessor(numberOfProcessor);
            PHMPI::PH_Bcast(&isConvertOver, 1*sizeof(int), serverProcessor);
            WriteLogFile("Convert grid complete !\n");
            PrintToWindow("Convert grid complete !\n\n");
        }
    }

    //! If need PartitionGrid.
    if (numberOfProcessor > 1)
    {
        if (!partitionFtsExist)
        {
            WriteLogFile("Begin partition grid !\n");
            PrintToWindow("Begin partition grid !\n\n");
            PHMPI::SetNumberOfProcessor(1);

            tasktype = PHSPACE::PARTITION_GRID;
            GlobalDataBase::UpdateData("nsimutask", &tasktype, PHINT, 1);

            int numberOfMultigrid = GlobalDataBase::GetIntParaFromDB("nMGLevel");
            int maxproc = numberOfProcessor;
            int pgridtype = GlobalDataBase::GetIntParaFromDB("gridtype");

            GlobalDataBase::UpdateData("pgridtype", &pgridtype, PHINT, 1);
            GlobalDataBase::UpdateData("maxproc", &maxproc, PHINT, 1);
            GlobalDataBase::UpdateData("numberOfMultigrid", &numberOfMultigrid, PHINT, 1);
            GlobalDataBase::UpdateData("partition_grid_file", &partition_grid_file, PHSTRING, 1);
            GlobalDataBase::UpdateData("original_grid_file", &outGridName, PHSTRING, 1);

            int isPartitionOver = 0;
            if (currentProcessor == serverProcessor)
            {
                PartitionGrid();
                isPartitionOver = 1;

                DeleteGGrids();
                PHMPI::SetNumberofLocalZones(0);
            }
            else
            {
                bool fileExist = true;
                bool fileCheck = true;
                PH_AllReduce(&fileExist, &fileCheck, 1, MPI_MIN);
            }

            PHMPI::SetNumberOfProcessor(numberOfProcessor);
            PHMPI::PH_Bcast(&isPartitionOver, 1 * sizeof(int), serverProcessor);
            WriteLogFile("Partition grid complete !\n");
            PrintToWindow("Partition grid complete !\n\n");
        }

        string walldistfile = ChangeExtensionOfFileName(partition_grid_file, "wdt");
        GlobalDataBase::UpdateData("gridfile", &partition_grid_file, PHSTRING, 1);
        GlobalDataBase::UpdateData("walldistfile", &walldistfile, PHSTRING, 1);
    }
    else
    {
        string walldistfile = ChangeExtensionOfFileName(outGridName, "wdt");
        GlobalDataBase::UpdateData("gridfile", &outGridName, PHSTRING, 1);
        GlobalDataBase::UpdateData("walldistfile", &walldistfile, PHSTRING, 1);
    }

    WriteLogFile("Begin CFD simulation !\n");
    PrintToWindow("Begin CFD simulation !\n\n");

    tasktype = PHSPACE::SOLVE_FIELD;
    GlobalDataBase::UpdateData("nsimutask", &tasktype, PHINT, 1);

    RunCFD();
}

void Simulation::RunPostProcessing()
{
    PostProcess *postProcess = new PostProcess(this->region);
    postProcess->Run();

    delete postProcess;
    postProcess = NULL;
}


#ifdef USE_SpecDiffHybSolver
void Simulation::RunSpecDiffHybSolver()
{
    double last_time[3];
    TIME_SPACE::StartTime(last_time);

    SpecDiffHybSolver *solver = new SpecDiffHybSolver();

    solver->Run();

    delete solver;

    WriteLogFile("SpecDiffHybSolver Finished!\n");
}

void Simulation::RunSpecSolver()
{
    double last_time[3];
    TIME_SPACE::StartTime(last_time);

    SpecSolver *solver = new SpecSolver();

    solver->Run();

    delete solver;
}
#endif

void Simulation::ComputeWallDistance()
{
    double last_time[3] = {0.0};
    TIME_SPACE::StartTime(last_time);

    ReadGeometry();

    UpdateAllZonesPara();

    InitGeometry();

    GlobalBoundaryCondition::ReadGlobalBoundaryCondition();

    GlobalBoundaryCondition::ChangeBCTypeByGlobalBC();

    InitGlobalBoundaryCondition();

    GenerateWalldist();

    WriteLogFile("Generate wall distance ...");
}

void Simulation::InitOversetGrid()
{
    //! The overset config process must be taken after reading in ordinary grid and initializing!
    //! Read overset information files or config.

    int isOverset = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");

    if (!isOverset)
    {
        return;
    }

    int readOversetFileOrNot = PHSPACE::GlobalDataBase::GetIntParaFromDB("readOversetFileOrNot");

    //! If the processor is grid processor, and overset grid need to be assembled.
    if (PHMPI::CurrentProcessorIsGridProcessor() && readOversetFileOrNot != 1)
    {
        OversetGridConfig();
        WriteLogFile("Run overset configuration ...");
    }

    //! Wait for the grid processors.
    PH_Barrier();

    //! For the original processors.
    if (!PHMPI::CurrentProcessorIsGridProcessor())
    {
        //! If the overset files is aviable or has been prepared by grid processors, just read it!
        if (readOversetFileOrNot == 1 || PHMPI::GetNumberOfGridProcessor() != 0)
        {
            region->ReadOversetGrid();
            WriteLogFile("Read overset information \n");
        }
        //! Else it should be done by these processors by themselves.
        else
        {
            TimeSpan *timeSpan = new TimeSpan();
            WriteLogFile("Run overset configuration ...\n");
            OversetGridConfig();
            timeSpan->ShowTimeSpanToLogFile("Run overset configuration");
            delete timeSpan;
        }

        InitializeOversetInterfaceTopology();

        WriteLogFile("Initialize oversetInterface topology \n");
    }
}

void Simulation::OversetGridConfig()
{
   ReConfigUnstructOversetGrid();
}

void Simulation::ReadGeometry()
{
    TimeSpan testTime;

    this->region->ReadGrid();

    testTime.ShowTimeSpanToLogFile("Grid Reading");
}

void Simulation::InitGeometry()
{
    this->region->InitGeometry();
}

void Simulation::SwapGeometry()
{
    this->region->SwapGeometry();
}

void Simulation::InitializeSolvers()
{
    WriteLogFile("Creating and initialization solvers start \n");

    int nSolver = GlobalDataBase::GetIntParaFromDB("maxNumOfSolver");
    for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
    {
        Controller *controller = new Controller(iSolver, this->region);

        controller->Initialize(iSolver);

        delete controller;
    }

    CreateAleManager();

    WriteLogFile("Creating and initialization solvers over \n");
}

void Simulation::RefreshSolvers()
{
    FreeAleManager();

    int nSolver = GlobalDataBase::GetIntParaFromDB("maxNumOfSolver");
    for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
    {
        Controller *controller = new Controller(iSolver, this->region);

        controller->CleanUp(iSolver);

        delete controller;
    }

    WriteLogFile("Creating and initialization solvers over \n");
}

void Simulation::DefineAndCreateSolvers()
{
    InitGlobalValuesOfAllZones();

    Region *region = this->region;

    int nZones = PHMPI::GetNumberofGlobalZones();
    int maxNumOfSolver = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = region->GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        vector<PHSolver *> solverList;
        PHGeometry *geometry = zone->GetGeometry();

#ifdef USE_DEMOSOLVER
        PHSolver *solver = new DemoSolver();
        solver->SetName("DemoSolver");
        solver->SetKey(DENSITY_BASED_SOLVER);
        solverList.push_back(solver);
#else

        //! For each zone, define the solver names on it.
        //! If you are solver developer, define your solver name here,
        //! according your solver logic.
        int dgHighOrder = 0;
        zone->GetData("dg_high_order", &dgHighOrder, PHINT, 1);

        //! First: NS equations solver.
        if (!dgHighOrder)
        {
            int gridType = geometry->GetGrid()->Type();

            int compressible = COMPRESSIBLE;
            zone->GetData("compressible", &compressible, PHINT, 1);

            if (compressible == COMPRESSIBLE)
            {
                //! Add compressible solvers.
                if (gridType == UNSTRUCTGRID)
                {
                    //! Add unstructured NS solver.
                    PHSolver *solver = new NSSolverUnstruct();
                    solver->SetName("NSSolverUnstruct");
                    solver->SetKey(DENSITY_BASED_SOLVER);
                    solverList.push_back(solver);
                }
                else if (gridType == STRUCTGRID)
                {
                    int isFVMOrFDM = 0;
                    zone->GetData("isFVMOrFDM", &isFVMOrFDM, PHINT, 1);

                    int isAdaptiveSolver = 0;
                    zone->GetData("isAdaptiveSolver", &isAdaptiveSolver, PHINT, 1);

                    if (isFVMOrFDM == FV_METHOD)    //! Add structured FV NS solver.
                    {
                        PHSolver *solver = new NSSolverStructFV();
                        solver->SetName("NSSolverStructFV");
                        solver->SetKey(DENSITY_BASED_SOLVER);
                        solverList.push_back(solver);
                    }
                    else    //! Add structured FD NS solver.
                    {
                        PHSolver *solver = new NSSolverStructFD();
                        solver->SetName("NSSolverStructFD");
                        solver->SetKey(DENSITY_BASED_SOLVER);
                        solverList.push_back(solver);
                    }
                }
            }
#ifdef USE_INCOMSOLVER
            else if (compressible == INCOMPRESSIBLE)
            {
                PHSolver *solver = new IncomUnstructSolver();
                solver->SetName("IncomUnstructSolver");
                solver->SetKey(PRESSURE_BASED_SOLVER);
                solverList.push_back(solver);
            }
#endif
        }
        else
        {
            //! To do
        }

        //! Second: turbulent equations solver.
        int viscousType = ONE_EQU;
        zone->GetData("viscousType", &viscousType, PHINT, 1);
        if (viscousType > LAMINAR)
        {
            int compressible = COMPRESSIBLE;
            zone->GetData("compressible", &compressible, PHINT, 1);
            if (compressible)
            {
                int gridType = geometry->GetGrid()->Type();
                if (gridType == UNSTRUCTGRID)
                {
                    //! Add unstructured turbulent solver.
                    PHSolver *solver = new TurbSolverUnstr();
                    solver->SetName("TurbSolverUnstr");
                    solver->SetKey(DENSITY_BASED_SOLVER);
                    solverList.push_back(solver);
                }
                else if (gridType == STRUCTGRID)
                {
                    int isFVMOrFDM = 0;
                    zone->GetData("isFVMOrFDM", &isFVMOrFDM, PHINT, 1);

                    if (isFVMOrFDM == FV_METHOD)
                    {
                        //! Add structured turbulent solver.
                        PHSolver *solver = new TurbSolverStr();
                        solver->SetName("TurbSolverStr");
                        solver->SetKey(DENSITY_BASED_SOLVER);
                        solverList.push_back(solver);
                    }
                    else
                    {
                        //! Add structured FD turbulent solver.
                        PHSolver *solver = new TurbSolverStrFD();
                        solver->SetName("TurbSolverStrFD");
                        solver->SetKey(DENSITY_BASED_SOLVER);
                        solverList.push_back(solver);
                    }
                }
            }
            else
            {
            }
        }

        if (viscousType > ONE_EQU)
        {
            int transitionType = 1;
            zone->GetData("transitionType", &transitionType, PHINT, 1);
            if (transitionType == IREGAMA)
            {
                int gridType = geometry->GetGrid()->Type();
                if (gridType == UNSTRUCTGRID)
                {
                    //! Add unstructured turbulent solver.
                    PHSolver *solver = new TransitionSolverUnstr();
                    solver->SetName("TransitionSolverUnstr");
                    solver->SetKey(DENSITY_BASED_SOLVER);
                    solverList.push_back(solver);
                }
                else if (gridType == STRUCTGRID)
                {
                    int isFVMOrFDM = 0;
                    zone->GetData("isFVMOrFDM", &isFVMOrFDM, PHINT, 1);

                    if (isFVMOrFDM == FV_METHOD)
                    {
                        //! Add structured turbulent solver.
                        PHSolver *solver = new TransitionSolverStr();
                        solver->SetName("TransitionSolverStr");
                        solver->SetKey(DENSITY_BASED_SOLVER);
                        solverList.push_back(solver);
                    }
                }
            }
            else
            {
            }
        }

        int iLES = 0;
        zone->GetData("iLES", &iLES, PHINT, 1);
        if (iLES == LES_SOLVER)
        {
            int gridType = geometry->GetGrid()->Type();
            if (gridType == UNSTRUCTGRID)
            {
                //! Add unstructured LES solver.
                PHSolver *solver = new LESSolverUnstr();
                solver->SetName("LESSolverUnstr");
                solver->SetKey(DENSITY_BASED_SOLVER);
                solverList.push_back(solver);
            }
            else if (gridType == STRUCTGRID)
            {
                //! Add structured LES solver.
                PHSolver *solver = new LESSolverStruct();
                solver->SetName("LESSolverStruct");
                solver->SetKey(DENSITY_BASED_SOLVER);
                solverList.push_back(solver);
            }
        }
#endif

        //! Finally, push the solvers into the list.
        for (std::size_t iSolver = 0; iSolver < solverList.size(); ++ iSolver)
        {
            PHSolver *solver = solverList[iSolver];

            solver->SetGeometry(geometry);
            solver->SetIndex(static_cast<int> (iSolver));
            zone->AddSolver(solver);
        }

        if (maxNumOfSolver < solverList.size())
        {
            maxNumOfSolver = solverList.size();
        }
    }   //! iZone
    GlobalDataBase::UpdateData("maxNumOfSolver", &maxNumOfSolver, PHINT, 1);
}

void Simulation::PreSolve()
{
    int isInIniting = GlobalDataBase::GetIntParaFromDB("isInIniting");
    int outIterStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
    int flowInitStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("flowInitStep");
    if (isInIniting == 1)
    {
#ifdef USE_GMRESSOLVER
        // GMRES Init by LU_SGS
        int flowInitMethod = GlobalDataBase::GetIntParaFromDB("flowInitMethod");
        int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
        if( 3 == flowInitMethod )
        {
            int originaltscheme = GlobalDataBase::GetIntParaFromDB("OriginalTscheme");
 
            if (originaltscheme == GMRES && outIterStep <= flowInitStep)
            {
                tscheme = LU_SGS;
                GlobalDataBase::UpdateData("tscheme", &tscheme, PHINT, 1);
            }

            if (originaltscheme == GMRES && outIterStep > flowInitStep)
            {
                tscheme = GMRES;
                GlobalDataBase::UpdateData("tscheme", &tscheme, PHINT, 1);
            }
        }
#endif
        int nMGLevel = GlobalDataBase::GetIntParaFromDB("nMGLevel");
        if (nMGLevel > 1)
        {
            //! For MG case, only initialize on coarse grid.
            isInIniting = 0;
            GlobalDataBase::UpdateData("isInIniting", &isInIniting, PHINT, 1);
        }
        else
        {
            //! For single grid case, check if in the state of initialization.
            if ((flowInitStep > 0 && outIterStep > flowInitStep) || flowInitStep == 0)
            {
                isInIniting = 0;
                GlobalDataBase::UpdateData("isInIniting", &isInIniting, PHINT, 1);

                PrintToWindow("End init flow by first order ... \n");
            }
        }
    }
}

void Simulation::PostSolve()
{
    using namespace PHMPI;

    //! Dump tecflow file.
    int intervalStepPlot = GlobalDataBase::GetIntParaFromDB("intervalStepPlot");
    int intervalStepSample = GlobalDataBase::GetIntParaFromDB("intervalStepSample");
    int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (outIterStep == 1 || outIterStep % intervalStepPlot == 0 || ((outIterStep + 1) % intervalStepPlot == 0 && !isUnsteady))
    {
        WriteVisualFile();
    }

    if (!GlobalDataBase::IsExist("ifSetDataMonitor", PHINT, 1))
    {
        return;
    }
    else
    {
        int ifSetDataMonitor = GlobalDataBase::GetIntParaFromDB("ifSetDataMonitor");
        if (!ifSetDataMonitor) return;

        if (outIterStep % intervalStepSample == 0)
        {
            Post_ProbesWrite generateProbesVarFile;
            generateProbesVarFile.Run();
        }
    }
}

void Simulation::MainSolve()
{
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (!isUnsteady)
    {
        SolveSteady();
    }
    else
    {
        SolveUnsteady();
    }
}

void Simulation::SolveOneOuterStep()
{
    PreSolve();

    MainSolve();

    PostSolve();
}

void Simulation::SolveSimulation()
{
    WriteLogFile("Start iterating ...");
    bool continueSolve = true;

    int maxSimuStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("maxSimuStep");
    int isOverset = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    int sysGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");

    //! Wait for all processors here!
    if (!isOverset || sysGridType == STRUCTGRID)
    {
        PrintToWindow("Wait for all processors before iteration ... \n");
        PH_Barrier();
        PrintToWindow("Wait is over,start iteration ... \n");
    }

    int outIterStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
    if (maxSimuStep == 0)
    {
        continueSolve = false;
    }

    while (continueSolve)
    {
        UpdateOuterIterationStep(outIterStep);

        SolveOneOuterStep();

        if (outIterStep >= maxSimuStep && maxSimuStep > 0)
        {
            continueSolve = false;
        }
    }

}

void Simulation::SelfAdaptionSolve()
{
    WriteLogFile("Start iterating ...");
    bool continueSolve = true;

    int maxSimuStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("maxSimuStep");
    int isOverset = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    int sysGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");

    //! Wait for all processors here!
    if (!isOverset || sysGridType == STRUCTGRID)
    {
        PrintToWindow("Wait for all processors before iteration ... \n");
        PH_Barrier();
        PrintToWindow("Wait is over,start iteration ... \n");
    }

    //! Set the flag to pause the numerical iteration.
    RDouble monitorResidual = 0.0;
    PHSPACE::GlobalDataBase::UpdateData("monitorResidual", &monitorResidual, PHDOUBLE, 1);
    int nPauseKey = 0;
    PHSPACE::GlobalDataBase::UpdateData("nPauseKey", &nPauseKey, PHINT, 1);

    int outIterStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
    if (maxSimuStep == 0)
    {
        continueSolve = false;
    }

    while (continueSolve)
    {
        UpdateOuterIterationStep(outIterStep);

        SolveOneOuterStep();

        int nPauseKey1 = PHSPACE::GlobalDataBase::GetIntParaFromDB("nPauseKey");

        if ((outIterStep >= maxSimuStep && maxSimuStep > 0) || nPauseKey1 == 1)
        {
            continueSolve = false;
        }
    }
}

void Simulation::PostSimulation()
{
    this->region->PostSimulation();

    int isOverset = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    if (isOverset)
    {
        FreeUnstructOversetGridConfig();
    }

    int IsUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");
    if (IsUnsteady && isAle)
    {
        //PostprocessAfterInnerIterationStepForAleSolvers();

        FreeAleManager();
    }

    return;
}

void Simulation::UpdateOuterIterationStep(int &outIterStep)
{
    outIterStep ++;
    GlobalDataBase::UpdateData("outnstep", &outIterStep, PHINT, 1);

    int newIterStep = GlobalDataBase::GetIntParaFromDB("newnstep");
    newIterStep ++;
    GlobalDataBase::UpdateData("newnstep", &newIterStep, PHINT, 1);

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (!isUnsteady)
    {
        return;
    }

    RDouble physicalTime = GlobalDataBase::GetDoubleParaFromDB("physicalTime");
    RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

    physicalTime += physicalTimeStep;

    GlobalDataBase::UpdateData("physicalTime", &physicalTime, PHDOUBLE, 1);

    int ifStaticsFlowField = GlobalDataBase::GetIntParaFromDB("ifStaticsFlowField");
    if (ifStaticsFlowField == 1)
    {
        int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");
        if (outIterStep >= startStatisticStep)
        {
            int nStatisticalStep = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");
            ++ nStatisticalStep;
            GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
        }
    }
#ifdef USE_GMRESSOLVER
    int flowInitMethod = GlobalDataBase::GetIntParaFromDB("flowInitMethod");

    if( 3 == flowInitMethod )
    {
        int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
        if ( newIterStep == 1 )
        {
            GlobalDataBase::UpdateData("OriginalTscheme", &tscheme, PHINT, 1);
        }
    }
#endif
}

void Simulation::UpdateAllZonesPara()
{
    this->region->UpdateAllZonesPara();
}

void Simulation::InitGlobalBoundaryCondition()
{
    this->region->ReSetBoundaryConditionByGlobalBC();
}

void Simulation::InitGlobalVolumeCondition()
{
    GlobalVolumeCondition::ReadGlobalVolumeCondition();
    GlobalVolumeCondition::SetVCDataByGlobalVC();
}

void Simulation::MultiGridInitFlow()
{
    WriteLogFile("Multi-Grid Init Flow ...");
    int nMGLevel = GlobalDataBase::GetIntParaFromDB("nMGLevel");
    if (nMGLevel <= 1)
    {
        return;
    }

    int outIterStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
    if (outIterStep > 0)
    {
        return;    //! Read from flow file.
    }

    int flowInitStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("flowInitStep");
    if (flowInitStep <= 0) return;

    int nSolver = 1;    //! Consider NS solver only.

    int isInMGIniting = 1;
    GlobalDataBase::UpdateData("isInMGIniting", &isInMGIniting, PHINT, 1);

    for (int iLevel = nMGLevel - 1; iLevel > 0; -- iLevel)
    {
        int levelStep = 0;
        int totalLevelStep = flowInitStep * iLevel;

        ostringstream oss;
        oss << "Start init flow on level " << iLevel << ", " << totalLevelStep << " iterations ...";
        WriteLogFile(oss.str());
        oss << endl;
        PrintToWindow(oss.str());

        while (levelStep < totalLevelStep)
        {
            ++ levelStep;

            UpdateOuterIterationStep(outIterStep);

            //using namespace PHMPI;
            for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
            {
                Controller *controller = new Controller(iSolver, this->region);
                controller->MultiGridInitFlow(iLevel);
                controller->PostSolve(iLevel);
                delete controller;    controller = nullptr;
            }
        }

        //! Interpolate from coarse grid to fine grid.
        for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
        {
            Controller *controller = new Controller(iSolver, this->region);
            controller->InterpolatFineGrid(iSolver, iLevel - 1);
            delete controller;    controller = nullptr;
        }
    }

    isInMGIniting = 0;
    GlobalDataBase::UpdateData("isInMGIniting", &isInMGIniting, PHINT, 1);

    int newIterStep = 0;
    GlobalDataBase::UpdateData("newnstep", &newIterStep, PHINT, 1);

    PrintToWindow("End init flow by coarse grid ... \n");
    WriteLogFile("End init flow by coarse grid ...");
}

void Simulation::SolveSteady()
{
    int nSolver = GlobalDataBase::GetIntParaFromDB("maxNumOfSolver");
    //! Get iteration number.
    int outerIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    int turbInterval = GlobalDataBase::GetIntParaFromDB("turbInterval");
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");

    if (tscheme != GMRES)
    {
        for (int iSolver = 0; iSolver < nSolver; ++iSolver)
        {
    //! iSolver = 1 or 2, is supposed as turbulence computing.
            if ((iSolver != NS_SOLVER) && (outerIterStep % turbInterval != 0))
            {
                return;
            }

            Controller *controller = new Controller(iSolver, this->region);

            controller->SolveSteadyField();

            controller->PostSolve();

            delete controller;
        }
    }
    else
    {
        //! NS and turbulence are coupled for calculation by GMRES.
        SolveSteadyByGMRES();
    }

}
void Simulation::SolveSteadyByGMRES()
{
    int nSolver = GlobalDataBase::GetIntParaFromDB("maxNumOfSolver");
    //! Get iteration number.
    int outerIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    int turbInterval = GlobalDataBase::GetIntParaFromDB("turbInterval");
  
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");

    if (viscousType == ONE_EQU)
    {
        Controller *controllerNS = new Controller(NS_SOLVER, this->region);
        Controller *controllerTurb = new Controller(TURB_SOLVER, this->region);

        controllerNS->SolveSteadyField();
        controllerTurb->SolveSteadyField();

#ifdef USE_GMRESSOLVER
        controllerNS->GMRESSolver_Coupled(0);
#endif
        controllerNS->PostSolve();
        controllerTurb->PostSolve();

        delete controllerNS;
        delete controllerTurb;
    }
    else if (viscousType == INVISCID || viscousType == LAMINAR)
    {
        Controller *controller = new Controller(NS_SOLVER, this->region);

        controller->SolveSteadyField();

        controller->PostSolve();

        delete controller;
    }
    else if (viscousType == TWO_EQU) //TWO_EQU turbulence are coupled will be done in the future.
    {
        for (int iSolver = 0; iSolver < nSolver; ++iSolver)
        {
            //! iSolver = 1 or 2, is supposed as turbulence computing.
            if ((iSolver == TURB_SOLVER) && (outerIterStep % turbInterval != 0))
            {
                return;
            }
            Controller *controller = new Controller(iSolver, this->region);
          
            controller->SolveSteadyField();
            controller->PostSolve();

            delete controller;
        }
    }

}
void Simulation::SolveUnsteady()
{
    int subTotalIterStep = GlobalDataBase::GetIntParaFromDB("subTotalIterStep");
    int minSubItertation = GlobalDataBase::GetIntParaFromDB("min_sub_iter");
    int maxSubItertation = GlobalDataBase::GetIntParaFromDB("max_sub_iter");
    RDouble toltalSubItertation = GlobalDataBase::GetDoubleParaFromDB("tol_sub_iter");

    int innerIterStep = 0;
    int nSolver = GlobalDataBase::GetIntParaFromDB("maxNumOfSolver");

    bool isSubIterationDump = false;
    GlobalDataBase::UpdateData("isSubIterationDump", &isSubIterationDump, PHBOOL, 1);

    bool isSubIterationConvergence = false;
    while (!isSubIterationConvergence)
    {
        innerIterStep += 1;
        GlobalDataBase::UpdateData("innstep", &innerIterStep, PHINT, 1);

        SolveSingleInnerIterationStepForAleSolvers();

        RDouble subIterationNormOneSolver = zero;
        RDouble maxSubIterationNormAllSolvers = -LARGE;

        for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
        {
            Controller *controller = new Controller(iSolver, this->region);

            controller->SolveSteadyField();

            //! Here is in the sub-iteration, PostSolve is only for the CommunicationInterfaceData and dumps the "res.dat", which contains index of the sub-iteration residual.
            controller->PostSolve();

            subIterationNormOneSolver = controller->ComputeSubIterationNorm(iSolver, 0);
            maxSubIterationNormAllSolvers = MAX(maxSubIterationNormAllSolvers, subIterationNormOneSolver);

            delete controller;
        }

        PH_CompareMaxMin(maxSubIterationNormAllSolvers, 1);

        //! Sub-iteration continue or not judgement.
        isSubIterationConvergence = IsSubIterationConvergence(innerIterStep, minSubItertation, maxSubItertation, toltalSubItertation, maxSubIterationNormAllSolvers);
    }

    //! Here is already (innstep == max_sub_iter || isSubIterationConvergence).
    isSubIterationDump = true;
    GlobalDataBase::UpdateData("isSubIterationDump", &isSubIterationDump, PHBOOL, 1);

    subTotalIterStep += innerIterStep;
    GlobalDataBase::UpdateData("subTotalIterStep", &subTotalIterStep, PHINT, 1);

    //! Update physical time step flow field after sub-iteration convergence.
    for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
    {
        Controller *controller = new Controller(iSolver, this->region);

        controller->UpdateUnsteadyFlow(iSolver, 0);

        //! Here is out of the sub-iteration, PostSolve is to dump various files according to the parameters in the "cfd_para.hypara", such as intervalStepFlow, intervalStepPlot, intervalStepForce, intervalStepRes, et al.
        controller->PostSolve();

        delete controller;
    }

    PostprocessAfterInnerIterationStepForAleSolvers();
}

bool Simulation::IsSubIterationConvergence(const int &innerIterStep, const int &minSubIteration, const int &maxSubIteration, const RDouble &tolalSubIteration, const RDouble &averageQuotientMaxInSolvers)
{
    bool isSubIterationConvergence = false;

    if (innerIterStep < minSubIteration)
    {
        isSubIterationConvergence = false;
    }
    else if (averageQuotientMaxInSolvers <= tolalSubIteration)
    {
        isSubIterationConvergence = true;
    }
    else if (innerIterStep >= maxSubIteration)
    {
        isSubIterationConvergence = true;
    }
    else
    {
        isSubIterationConvergence = false;
    }

    return isSubIterationConvergence;
}


}

