//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      CFDSolver.h
//! @brief     CFD solver.
//! @author    He Xin, Bell, He Kun, He Xianyao, Liu Jian, Zhang Zipei, Li peng, 
//!            Ma Yankai, Zhang Yang, Wan Yunbo, Xu Gang.

#pragma once
#include "Solver.h"
#include "Geometry.h"
#include "Post_Visual.h"
#include "Post_Probes.h"
#include "IO_HDF5File.h"
#pragma warning(disable:4100)
#ifdef USE_GMRESSOLVER
#include "petscksp.h"
#endif
namespace PHSPACE
{
// #define GMRESWholeJacobian
class ActionKey;
class ActionTag;
class Grid;
class InterfaceDataProxy;
class FieldProxy;
class InterfaceInfo;
class Param_CFDSolver;
class Region;

extern bool cfd_debug;

//! @brief Base Class for CFD solvers, some other actual solvers are inherited from CFDSolver, such as NSSolver, TurbSolver.
class CFDSolver : public PHSolver
{
public:
    CFDSolver();
    ~CFDSolver();

public:
    //! Register an interface field (variable) into buffer, used in parallel communication.
    //! The interface field variable buffer is set to be ZERO.
    //! @param[in] name         variable name of the interface field that need to be communicated.
    //! @param[in] type         variable type of the interface field
    //! -# PHINT   : Integer type.
    //! -# PHFLOAT : Float type.
    //! -# PHDouble: Double type.
    //! @param[in] dimension    variable dimension
    //! -# 1 :  1-D array.
    //! -# 2 :  2-D array.
    //! -# n :  n-D array.
    LIB_EXPORT void RegisterInterfaceField(const string &name, const int type, const int dimesion);
    LIB_EXPORT void RegisterInterpointField(const string & name, const int type, const int dimesion);
    LIB_EXPORT void RegisterOversetField(const string &name, const int type, const int dimesion);

    LIB_EXPORT void ReleaseInterfaceField();
    LIB_EXPORT void ReleaseInterfaceField(string varName);

    void SetPostVisualization(Post_Visual *postVisualization_in);

    Post_Visual *GetPostVisualization();

    void SetPostVisualizationWall(Post_VisualWall *postVisualization_in);

    Post_VisualWall *GetPostVisualizationWall();

    //! Get control parameters
    LIB_EXPORT Param_CFDSolver *GetControlParameters() const;

    //! Compute pressure factor for entrop fix (Method 2 in Structure solver).
    virtual void CalPressureFactor(Grid *grid) {};

    //! Compute gradient of primitivevariables & temperature for viscous flux computation.
    virtual void ComputeGradient(Grid *grid);
    virtual void ComputeLimiter(Grid *grid);
    virtual void ComputeGradientCellCenter(Grid *grid);
    virtual void ComputeGradientCellCenterForLES(Grid *grid);
    virtual void ComputeGradientCellCenterOfruvwptOnlyForMixGrid(Grid *grid){};
    virtual void ComputelimiterOnlyForMixGrid(Grid *grid){};
    //! Compute blend function F1 & F2.
    virtual void Crossing(Grid *grid);
    virtual void Blending(Grid *grid);

    virtual void UpdateQlQrOnlyforMixGrid(Grid *grid){};
    virtual void UpdateInviscidfluxOnlyforMixGrid(Grid *grid){};
    virtual void UpdateViscousfluxOnlyforMixGrid(Grid *grid){};
    virtual void LoadGlobalfacefluxtoResOnlyforMixGrid(Grid *grid){};

    //! Low-speed Precondition Coefficient.
    virtual void ComputePreconditionCoefficient(Grid *grid);

    //! Low-speed Precondition Coefficient for unsteady.
    virtual void ComputePreconditionCoefficientUnsteady(Grid *grid);

    //! Low-speed Conservative Precondition Matrix.
    //! solve M * T-, when sign = -1;
    virtual void ConservativePreconMatrix(Grid *grid, int sign);

    //! Inviscid flux.
    virtual void InviscidFlux(Grid *grid);

    //! Viscous flux.
    virtual void ViscousFlux (Grid *grid);

    //! Source flux, such as chemical source term.
    virtual void SourceFlux(Grid *grid);

    //! To get the viscous term, viscouslaminar for ns-solver and viscousturbulent for turbulent solver.
    virtual void ComputeViscousCoeff(Grid *grid) {};
    virtual void ReconGrad(Grid *gridIn, int iSurface) {};

    virtual void GetdtwithBC(Grid *grid){};
    virtual void GetspecwithBC(Grid *grid){};
    virtual void SolveLUSGSForwardwithBC(Grid *grid, FieldProxy *dqProxy){};
    virtual void SolveLUSGSBackwardwithBC(Grid *grid, FieldProxy *dqProxy){};
    virtual void CompressSpecificArrayToInterfaceWithSpecificGhost(DataContainer *&dataContainer, const string &fieldName, Grid *grid, const int &neighborZoneIndex, const int &neqn, int nGhost) {};
    virtual void DecompressArrayFromInterfaceWithSpecificGhost(DataContainer *&dataContainer, const string &fieldName, Grid *grid, const int &neighborZoneIndex, const int &neqn, int nGhost) {};

    //! Time step computation.
    virtual void TimeStep(Grid *grid);

#ifdef USE_GMRESSOLVER
    //! GMRES solver -- a linear system solver
    //! This method can be referenced "Sen Zhang, Boqian Wang, et al. Implementation of a Newton-Krylov Algorithm in the Open-source Solver PHengLEI[C]."
    //! GPPS Hong Kong , October 17-19, 2023.
    virtual void GMRESSolver(Grid *gridIn, FieldProxy *dqProxy);
#endif

    //! Compute the inviscous/viscous/chemical spectrum radius.
    //! The results stored in array invSpectralRadius and visSpectralRadius and srs.
    virtual void SpectrumRadius(Grid *grid){};
    
    //! See the explain in controller.
    virtual void Diagonal(Grid *grid) {};
    
    virtual void ComputeMinTimeStep(Grid *grid, RDouble &minDt, RDouble &maxDt);
    virtual void ReduceMaxTimeStep(Grid *grid, RDouble globalMinDt);

    //! Explicit Runnge-Kutta time integration method.
    virtual void RungeKutta(Grid *grid, FieldProxy *rhsProxy);

    //! Implicit LU-SGS time integration method.
    virtual void SolveLUSGSForward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep = false);
    virtual void SolveLUSGSBackward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep = false);
   
    //! Implicit Block LU-SGS time integration method.
    virtual void SolveBLUSGSForward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);
    virtual void SolveBLUSGSBackward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);

    //! Implicit Matrix LU-SGS time integration method.
    virtual void SolveMatrixLUSGSForward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);    
    virtual void SolveMatrixLUSGSBackward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);

    //! Implicit Line LU-SGS time integration method.
    virtual void SolveLineLUSGSForward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);
    virtual void SolveLineLUSGSBackward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);

    //! To determine the CFL number of the next iteration step.
    virtual void DetermineCFLNumber(Grid *grid, FieldProxy *dqProxy, FieldProxy *dqProxy1);

    //! Line implicit LU-SGS time integration method.
    virtual void LineLUSGS(Grid *grid, int level);

    virtual void LHSMatrixs(Grid *grid);

    //! calculate the three Matrixs of all lines.
    virtual void LineLUSGSMatrixs(Grid *grid);

    //! Compute the node value in the flow.
    //! @param[in] grid         the mesh in the computing.    
    virtual void ComputeNodeValue(Grid *grid) {};

    //! Modify the node value in the flow.
    //! @param[in] grid         the mesh in the computing.
    virtual void ModifyNodeValue(Grid *grid) {};

    //! Compute the node value of the the interpoint for turbulent flow.
    //! @param[in] grid         the mesh in the computing.
    virtual void ComputeQTurbNodeValue(Grid *grid) {};

    //! Compute the node value of the the interpoint for transition flow.
    //! @param[in] grid         the mesh in the computing.
    virtual void ComputeQTransitionNodeValue(Grid *grid) {};

    //! Modify the node value of the the interpoint for turbulent flow.
    //! @param[in] grid         the mesh in the computing.
    virtual void ModifyQTurbNodeValue(Grid *grid) {};
    //! Put correction on the coarse grid, if multi-grid method is used.
    void PutCorrection(Grid *grid, FieldProxy *qProxy);

    virtual RDouble UnsteadyConvergence(Grid *grid);

    virtual void UpdateUnsteadyFlow(Grid *grid);
    virtual void UpdateUnsteadyVariable(Grid *grid);
    virtual void UpdateDt();

    virtual int GetNPostSolve();
    virtual int GetNumberOfEquations();

    virtual void PostSolve(Grid *grid, int stage, int level = 0);

    virtual void DumpResultFile(Grid *grid, int level = 0);

    virtual void CheckResult(Grid *grid, int level = 0);

    virtual void CheckResiduals(Grid *grid, int level = 0);

    virtual void DumpCFLNumber(Grid *grid, int level = 0);

    virtual void CheckSurfaceHeatingChange(Grid *grid, int level = 0);

    virtual void DumpLeakageInformation(Grid *grid);

    virtual void ZeroResiduals(Grid *grid);

    virtual void UpdateResiduals(Grid *grid);

    virtual void ComputeResidual(Grid *grid, FieldProxy *rhsProxy);

    virtual void RungeKuttaResidual(Grid *grid, FieldProxy *dqProxy, RDouble coef);

    virtual void RestrictAllQ(Grid *fineGrid, Grid *coarseGrid);

    virtual void RestrictDefect(Grid *fineGrid, Grid *coarseGrid);

    virtual void PutCorrectionBack(Grid *grid, FieldProxy *qProxy);

    virtual void CorrectFineGrid(Grid *fineGrid, Grid *coarseGrid);

    virtual void InterpolatFineGrid(Grid *fineGrid, Grid *coarseGrid);

    virtual void StoreRhsByResidual(Grid *grid, FieldProxy *rhsProxy);

    virtual void InitResidual(Grid *grid, FieldProxy *rhsProxy);

    virtual FieldProxy * GetResidualProxy(Grid *grid);

    virtual void RecoverResidual(Grid *grid, FieldProxy *rhsProxy);

    virtual void LoadQ(Grid *grid, FieldProxy *qProxy);

    virtual FieldProxy * CreateFieldProxy(Grid *grid);

    virtual void InitCGrid(Grid *fineGridIn, Grid *coarseGridIn) {};

    virtual void LUSGSInitializationStructHighOrder(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ){};
    
    //! Test start.
    virtual void LoadResiduals(Grid *grid, FieldProxy *rhsProxy);

    virtual void Boundary(Grid *grid);

    virtual void CalculateBoundaryData();

    virtual void GetIndependentVariablesforStructHighOrder(Grid *grid){};

    virtual void GetDependentVariablesforStructHighOrder(Grid *grid){};

    virtual void ComputeGamaAndTemperature(Grid *grid) {};

    virtual void GetResidual(ActionKey *actkey);

    virtual void GetResidual(ActionKey *actkey, RDouble &localMaxRes);

    virtual void GetResidual(Grid *gridIn, RDouble &localMaxRes);

    //! Obtain the CFL number.
    virtual void ObtainCFLNumber(Grid *gridIn, RDouble &globalMinCFL, RDouble &globalMaxCFL);

    //! Compute the surface heating change and obtain the maximum change of the surface heating ratio.
    virtual void GetSurfaceHeatingChange(ActionKey *actkey, RDouble &localMaxHeatChange);

    //! Prepare post visualization variables, store it in to grid's postVisualization database.
    LIB_EXPORT virtual void ComputePostVisualVariables(Post_Visual *postVisualization) {};
    //! Prepare monitored probes variables for post-processing , store it into postProbesVar database.
    LIB_EXPORT virtual void ComputePostProbesVariables(Post_Probes *postProbesVar) {};

    virtual FieldProxy * GetFieldProxy(Grid *grid, const string &field_name);

    virtual string CastName(const string &name);

    virtual void UpdateFlowField(Grid *grid, FieldProxy *qProxy, FieldProxy *dqProxy);

    virtual void CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);
    virtual void DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);
    
    virtual void CompressAnInterfaceVar(DataContainer *&dataContainer, RDouble **fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);
    virtual void DecompressAnInterfaceVar(DataContainer *&dataContainer, RDouble **fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);

    virtual void CompressSpecificArrayToInterface(DataContainer *&dataContainer, const string &fieldName, Grid *grid, const int &neighborZoneIndex, const int &neqn) {};
    virtual void DecompressArrayFromInterface(DataContainer *&dataContainer, const string &fieldName, Grid *grid, const int &neighborZoneIndex, const int &neqn) {};
    
    //! Rotate gradient vector for periodic boundary in both NSSolver and TurbSolver.
    virtual void RotateVectorFromInterface(Grid *grid, const int &neighborZoneIndex, const int &neqn) {};

    virtual void AverageMixingPlane(Grid *grid) {};
    virtual void MixingPlaneDataTransfer(Grid *grid, Grid *NeighborGrid1) {};
    virtual void NonReflective(Grid *grid, Grid *NeighborGrid1) {};
    virtual void SetMixingPlaneData(Grid *grid) {};

    virtual void CompressGradientandLimitToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex){};
    virtual void DecompressGradientandLimitToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex){};

    virtual void CompressQlQrToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex){};
    virtual void DecompressQlQrToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex){};


    virtual void CompressFacefluxToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex){};
    virtual void DecompressFacefluxToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex){};

    //! Postprocess such as parallel communication, flowfield and residual info output, et.al.
    virtual void AirForceCoef(ActionKey *actkey){};
    virtual void CpDistriCoef(ActionKey *actkey){};

    virtual void ComputeHeatFlux(ActionKey *actkey, RDouble *HeatFlux){};

    //! Interpoint data communications.
    virtual void CommunicationInterpointData();

    virtual void DumpResidual(ActionKey *actkey);

    virtual void CheckResidual(ActionKey *actkey);

    virtual void ObtainCFLNumber(ActionKey *actkey);

    virtual void ComputeSurfaceHeatingChange(ActionKey *actkey);

    //! Get the file name for Residual dumping.
    virtual const string GetResidualFileName();

    //! Get the file name for Restart info dumping.
    virtual const string GetRestartFileName();

    //! Get the grid file name for Restart info by interpolated.
    virtual const string GetRestartGridFileName();

    //! Get the flow file name for Restart info by interpolated.
    virtual const string GetRestartFlowFileName();

    //! Get the file name for AirCoef dumping.
    virtual const string GetAirCoefFileName();

    //! Get the file name for flow field output by tecplot format.
    virtual const string GetFlowTecFileName();

    virtual void DumpRestartData(ActionKey *actkey);
    virtual void DumpProtectData(ActionKey *actkey, int &RestartFile);
    virtual void DumpTurbProtectData(ActionKey *actkey, int &RestartFile);
    virtual void DumpTransitionProtectData(ActionKey *actkey, int &RestartFile);
    virtual void CreateH5RestartFile(ActionKey *actkey);

    virtual void InitMemory();
    virtual void ReleaseMemory();

    //! Initialize the flow by two ways:\n
    //! 1. InitFlowAsRestart\n
    //! 2. InitFlowAsReadingRestart.
    virtual void InitFlow();

    virtual void InitProtectedFlow();

    //! Initialize flow by inflow condition, or set the flow variables to be zero.
    virtual void InitFlowAsRestart();

    //! Initialize flow by reading restart file.
    //! For example, for NS solver, the flow variables are stored in flow.dat
    //! for turbulent solver, the turbulent variables are store in turb.dat.
    virtual void InitFlowAsReadingRestart();
    virtual void InitFlowAsReadingRestart(ActionKey *actkey);
    virtual void InitFlowByReadingInterpolate(ActionKey *actkey);

    virtual void InitFlowAsReadingProtectedRestart(ActionKey *actkey);

    //! To output Average Flow for visualization by tecplot format.
    virtual void TecOutAverageFlow(ActionKey *actkey);

    //! To output Average Reynolds stress for visualization by tecplot format.
    virtual void TecOutAverageReynoldsStress(ActionKey *actkey);

    virtual void DumpAirForceCoef(ActionKey *actkey);
    virtual void DumpCpDistri(ActionKey *actkey);
    virtual void DumpSurfaceInfo(ActionKey *actkey);
    void WriteWallAircoefASCII(ActionKey *actkey, vector <DataContainer *> datalist);
    virtual void DumpHeatFlux(ActionKey *actkey);
    void WriteHeatFluxASCII(ActionKey *actkey, RDouble *HeatFlux);
    #ifdef USE_TecplotLib
        void WriteWallAircoef(ActionKey *actkey, vector <DataContainer *> datalist);
    #endif

    //! Export the required variables on surface to file, such as the temperature, mass fraction etc.
    virtual void ExportSurfaceVariables(ActionKey *actkey) {}; //added by Li Peng

    //! Interface data communications.
    void CommunicationInterfaceData();
    void CommunicationInterfaceData(int iLevel);

    //! Interface information exchange controller process. A complete process include three steps: \n
    //! 1. Upload Interface Data. \n
    //! 2. Update Interface Data. \n
    //! 3. Download Interface Data.
    void UploadInterfaceData();
    void UpdateInterfaceData();
    void DownloadInterfaceData();

    void UploadInterfaceData(int iLevel);
    void UpdateInterfaceData(int iLevel);
    void DownloadInterfaceData(int iLevel);

    //! Interpoint information exchange controller process. A complete process include three steps: \n
    //! 1. Upload Interpoint Data. \n
    //! 2. Update Interpoint Data. \n
    //! 3. Download Interpoint Data.
    void UploadInterpointData();
    void UpdateInterpointData();
    void DownloadInterpointData();
    void CommunicationInterpointWeight();

    void CommunicationOversetSlipData();
    //! Overset information exchange controller process, including 3 steps: \n
    //! 1. Upload Overset Data. \n
    //! 2. Update Overset Data. \n
    //! 3. Download Overset Data.
    void UploadOversetData();
    void UpdateOversetData();
    void DownloadOversetData();

    virtual void DumpRestart(ActionKey *actkey);
    virtual void ReadRestart(ActionKey *actkey);

    virtual void DumpRestartH5(ActionKey *actkey);
    virtual void ReadRestartH5(ActionKey *actkey);
    virtual void InterpolateFromRestartH5(ActionKey *actkey);
    virtual void FillField(Grid *grid, FieldProxy *fieldProxy, RDouble value);
    virtual void FillField(Grid *grid, FieldProxy *field1Proxy, FieldProxy *field2Proxy);

    void Solve();

protected:
    int GetZoneLocalID();

    Grid * GetGrid(int level = 0);

    uint_t GetNumberOfMultiGrid();

    void Action(ActionKey *actkey);

    virtual void FillActionKey(ActionKey *actkey, int action, int level){};

    void TranslateAction(ActionKey *actkey);

    streamsize TranslateActionLength(ActionKey *actkey);

    virtual bool JudgeIfReadAverage();
    virtual bool JudgeIfRestart();
    virtual bool JudgeIfInterpolate();
    virtual bool JudgeIfProtectedRestart();

    virtual void InitCoarseGridsFlow();

    //! Create and init control parameters.
    virtual void InitControlParameters();

    //! Create and init control parameters.
    virtual void InitDependentVariables();

    //! Free control parameters.
    virtual void FreeControlParameters();

    virtual void ReadParameter();

    virtual void ActionReflect(ActionTag *acttag);
    virtual void VisualizationAverageFlow(ActionKey *actkey);

    void AllocateGlobalVariables();
    void DeAllocateGlobalVariables();
    virtual void AllocateGlobalVar(Grid *grid);
    virtual void DeAllocateGlobalVar(Grid *grid);

    void DeAllocateOversetInterfaceVar(Grid *grid);

    void DeAllocateOversetInterfaceVar(OversetInformationProxy *oversetInformationProxy);

    virtual void DeAllocateOversetInterfaceVar(Data_ParamFieldSuite *dataStore);

    virtual void GetInterfaceData  (ActionKey *actkey);
    virtual void TranslateInterfaceData(ActionKey *actkey);
    virtual streamsize TranslateInterfaceDataLength(ActionKey *actkey);

    //! Translate the interpoint data length.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    virtual streamsize TranslateInterpointDataLength(ActionKey *actkey);

    //! Get the interpoint data for sending data.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    virtual void GetInterpointData (ActionKey *actkey);

    //! Translate the interpoint data using MPI.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    virtual void TranslateInterpointData(ActionKey *actkey);

    virtual void PrepareInterfaceData(Grid *grid, Data_ParamFieldSuite *dataStore, InterfaceDataProxy *interfaceDataProxy);

    //! Prepare the interpoint data for sending.
    //! @param[in] grid         the mesh in the computing.
    //! @param[in] datastore       the data store in current zone.
    //! @param[in] interpointDataProxy         a proxy to store the interpoint value.
    virtual void PrepareInterpointData(Grid *grid, Data_ParamFieldSuite *datastore, InterpointDataProxy *interpointDataProxy);

    virtual void PrepareOversetInterfaceData(Data_ParamFieldSuite *dataStore, InterfaceDataProxy *interfaceDataProxy);

    void ReadInterface (ActionKey *actkey);
    void WriteInterface(ActionKey *actkey);

    virtual void UploadInterfaceValue(ActionKey *actkey);
    virtual void DownloadInterfaceValue(ActionKey *actkey);

    virtual void UploadInterfaceData(ActionKey *actkey);
    virtual void DownloadInterfaceData(ActionKey *actkey);

    virtual void UploadOversetData(ActionKey *actkey);
    virtual void DownloadOversetData(ActionKey *actkey);

    //! Upload the interpoint variables for send.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    virtual void UploadInterpointData(ActionKey *actkey);

    //! Download the interpoint variables from the receive data.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    virtual void DownloadInterpointData(ActionKey *actkey);

    virtual void CommunicationInterpointWeight(ActionKey *actkey);
    virtual void DownloadInterpointWeight(ActionKey *actkey);

    void GetOversetData(ActionKey *actkey);
    void TranslateOversetData(ActionKey *actkey);

    virtual void ModifyResiduals(Grid *grid);

    void SolveOnGridSteady  (Grid *grid);
    void SolveOnGridUnsteady(Grid *grid);

    void SolveOnGrid(Grid *grid);

    virtual void DumpInterfaceInfoBell();
    virtual void RegisterCFDSolverInterfaceField();
    void RegisterOversetField();

    //! Register the interpoint field in the CFD solver.
    virtual void RegisterCFDSolverInterpointField();

    bool IsNeedStatistics();

    bool IsNeedReynoldsStressStatistics();

    void ComputeDualTimeCoefficient(int methodOfDualTime, RDouble *dualTimeCoefficient);

    virtual void FreeGradientProxy(Grid *grid) {};

private:
    void RegisterInterfaceField(Grid *grid, const string &name, const int type, const int dimesion);

    //! register the interpoint value field for translating.
    //! @param[in] grid         the mesh in the computing.
    //! @param[in] name         name of the value field.
    //! @param[in] type         type of the value field.
    //! @param[in] dimesion      dimesion of the value field.
    void RegisterInterpointField(Grid *grid, const string &name, const int type, const int dimesion);

    void GetInterfaceData(InterfaceInfo *iterfaceInfo, ActionKey *actkey, InterfaceDataProxy *interfaceDataProxy);    void TranslateInterfaceData(InterfaceInfo *iterfaceInfo, ActionKey *actkey, InterfaceDataProxy *interfaceDataProxy);
    streamsize TranslateInterfaceDataLength(InterfaceInfo *iterfaceInfo, ActionKey *actkey, InterfaceDataProxy *interfaceDataProxy);
    
    //! Translate the interpoint data length using MPI.
    //! @param[in] interpointInformation         interpoint information for store the interpoint.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    //! @param[in] interpointDataProxy         a proxy to store the interpoint value.
    streamsize TranslateInterpointDataLength(InterpointInformation *interpointInformation, ActionKey *actkey, InterpointDataProxy *interpointDataProxy);
    void GetOversetInterfaceData  (OversetInformationProxy *oversetInfoProxy, ActionKey *actkey, InterfaceDataProxy *interfaceDataProxy);
    void TranslateOversetInterfaceData(OversetInformationProxy *oversetInfoProxy, ActionKey *actkey, InterfaceDataProxy *interfaceDataProxy);
    void PutCorrection(Grid *grid, FieldProxy *fieldProxy, FieldProxy *old_field_proxy, int neqn);

    //! Get the interpoint data.
    //! @param[in] interpointInformation         interpoint information for store the interpoint.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    //! @param[in] interpointDataProxy         a proxy to store the interpoint value.
    void GetInterpointData (InterpointInformation *interpointInformation, ActionKey *actkey, InterpointDataProxy * interpointDataProxy);

    //! Translate the interpoint data using MPI.
    //! @param[in] interpointInformation         interpoint information for store the interpoint.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    //! @param[in] interpointDataProxy         a proxy to store the interpoint value.
    void TranslateInterpointData(InterpointInformation *interpointInformation, ActionKey *actkey, InterpointDataProxy * interpointDataProxy);

public:
    virtual void SolveIncomSteadyField() {};

protected:
    Post_Visual *postVisualization;

    Post_VisualWall *postVisualWall;

    FieldProxy **dqProxy;

    Param_CFDSolver *controlParameters;
};

void UpdateInterfaceData();

#include "CFDSolver.hxx"

#ifdef USE_GMRESSOLVER
PetscErrorCode MyKSPMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void *dummy);
#endif

}
