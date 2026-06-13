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
//! @file      Geo_Grid.h
//! @brief     It defines the class 'Grid', which is the base class of geometry grid.
//!            The inheriting order is: SimpleGrid -> Grid -> StructuredGrid/UnstructuredGrid.
//! @author    Bell, He Xin.

#pragma once
#include "LinkStruct.h"
#include "Geo_GridIndex.h"
#include "DataStruct_Sort.h"
#include "Geo_Interpoint.h"
#include "Geo_Interface.h"
#include "ActionKey.h"
#include "Data_Field.h"
#include "Data_Param.h"
#include "FieldProxy.h"
#include "OversetInformation.h"
#include "Geo_SimpleVC.h"

namespace PHSPACE
{
//! @brief It defines the class 'Grid', which is the base class of Structured or Unstructured geometry grid.
//!        It is usually used in the whole code, it could be used to represent both Structured or Unstructured grid.
//!        The inheriting order is: SimpleGrid -> Grid -> StructuredGrid/UnstructuredGrid.
class Grid : public SimpleGrid
{
public:
    Grid();
    virtual ~Grid();
    
protected: 
    //! Grid dimension.
    //! -# 2: two dimensional grid.
    //! -# 3: three dimensional grid.
    int dimension;

    //! The grid type.
    //! -# 0: UNSTRUCTGRID.
    //! -# 1: STRUCTGRID.
    //! -# 2: MIXGRID.
    int type;

    //! Grid level in Multi-Grid method, level = 0, 1, 2, 3 ...
    //! It starts from 0, which is the original finest level.
    //! The coarse level is obtained by agglomerate the finer level grid.
    int level;

    //! The number of Total Faces.
    int nTotalFace;

    //! The number of Total Cells.
    int nTotalCell;

    //! The number of Total Cells.
    int nTotalLine;

    //! The number of boundary faces which include interface faces.
    int nBoundFace;

    //! The number of interface faces.
    //! Interface faces: faces between different zones,
    //! such as faces between zones in different unstructured grid partitions,
    //! and faces between zones in different structured grid blocks.
    int nIFace;

    //! Interface topology information of interface faces.
    InterfaceInfo *interfaceInfo;

    //! Interpoint topology information of interface points.
    InterpointInformation *interpointInformation;

    //! Interface fields information.
    //! Stores the field variables that need to communicate between different gird of zones.
    InterfaceFields *interfaceFields;

    //! Interpoint fields information.
    //! Stores the field variables that need to communicate between different grid of zones.
    InterpointFields *interpointFields;

    //! Coarser grid of the current grid.
    //! 'c' represents 'coarse'.
    Grid *cGrid;

    //! Finer grid of the current grid.
    //! 'f' represents 'fine'.
    Grid *fGrid;

    //! Store any type of flow field variables, such as velocity, time, temperature.
    Data_Field *gField;

    //! Store any type of control parameters, such as int, double, string.
    Data_Param *gPara;
  
    //! Grid index.
    GridID *index;

    //! Overset or overlapping Grid information.
    OversetInfoProxy *oversetInfoProxy;

    //! Overset information for unstructured overset cell interpolation(need further revise)
    OversetInformationProxy *oversetInformationProxy;

    //! The grid at the last time.
    //! It is used in unsteady dynamic grid method.
    SimpleGrid *oldGrid;

    //! which block the grid belongs to.
    int iBlock;

    //! The probes' number of the current grid.
    int zoneProbesNumber;

    //! The probes' global index of the current grid.
    vector<int> zoneProbesGlobalID;

    //! The probes' line index of the current grid.
    vector<int> zoneProbesLineID;

    //! The probes' surface index of the current grid.
    vector<int> zoneProbesSurfaceID;

    //! The probes' coordinates of the current grid.
    vector<vector<RDouble> > zoneProbesCoordinates;

    //! The grid file index current grid belong, only used in multi-file grid convert.
    int fileIndexCurrentGridBelong;

    SimpleVC *volumeCondition;

private:
    //! 
    int own_database;

public:
    //! Assign the value to the grid dimension.
    void SetDim(int dim);

    //! Assign the value to the multi-grid level.
    void SetLevel(int level);

    //! Assign the value to the number of Total faces.
    void SetNTotalFace(int nTotalFace);

    //! Assign the value to the number of Total Cells.
    void SetNTotalCell(int nTotalCell);

    //! Assign the value to the number of boundary faces which include interface faces.
    void SetNBoundFace(int nBoundFace);

    //! Assign the value to the number of interfaces.
    void SetNIFace(int nIFace);

    //! Assign the pointer to the interface topology information of interface faces.
    void SetInterfaceInfo(InterfaceInfo *interfaceInfo);

    //! Assign the pointer to the interface topology information of interface faces.
    void SetInterfaceFields(InterfaceFields *interfaceFields) { this->interfaceFields = interfaceFields; }

    //! Assign the pointer to the interpoint topology information of interface points.
    //! @param[in] interpointInformation    an object of class InterpointInformation.
    void SetInterpointInfo(InterpointInformation *interpointInformation);

    //! Assign the pointer to coarser grid of the current grid;
    //! 'c' represents 'coarse'.
    void SetCoarseGrid(Grid *cGrid);

    //! Assign the pointer to finer grid of the current grid;
    //! 'f' represents 'fine'.
    void SetFineGrid(Grid *fGrid);

    //! Assign the value to the global zone grid index.
    void SetZoneID(int index);

    //! Assign the value to the local zone grid index.
    void SetZoneLocalID(int index);

    //! Assign the pointer to the overset or overlapping Grid information.
    void SetOversetInfoProxy(OversetInfoProxy *oversetInfoProxy);

    //! Assign the pointer to the overset information for unstructured overset cell interpolation.
    void SetOversetInformationProxy(OversetInformationProxy *oversetInformationProxyIn);

    //! Assign iBlock.
    void SetIBlock(int iBlockIn);

    //! Assign the probes number of each zone grid.
    void SetZoneProbesNumber(int zoneProbesNumberIn);

    //! Assign the probes global index of each zone grid.
    void SetZoneProbesGlobalID(vector<int> &zoneProbesGlobalIDIn);
    //! Assign the probes line index of each zone grid.
    void SetZoneProbesLineID(vector<int> &zoneProbesLineIDIn);
    //! Assign the probes surface index of each zone grid.
    void SetZoneProbesSurfaceID(vector<int> &zoneProbesSurfaceIDIn);

    //! Assign the probes coordinates of each zone grid.
    void SetZoneProbesCoordinates(vector<vector<RDouble> > &zoneProbesCoordinatesIn);

    void SetVolumeCondition(SimpleVC *volumeConditionIn);
    SimpleVC *GetVolumeConditionIn();

    //! Return the grid dimension.
    int GetDim() const;

    //! Return the grid type.
    //! -# 0: UNSTRUCTGRID
    //! -# 1: STRUCTGRID
    //! -# 2: MIXGRID
    int Type() const;

    //! Return the grid level in Multi-Grid method, level = 0, 1, 2, 3 ...
    //! It starts from 0, which is the original finest level.
    //! The coarse level is obtained by agglomerate the finer level grid.
    int GetLevel() const;

    //! Return the number of Total Faces.
    int GetNTotalFace() const;

    //! Return the value to the number of Total Cells.
    int GetNTotalCell() const;

    //! Return the value to the number of boundary faces which include interface faces.
    int GetNBoundFace() const;

    //! Return the value to the number of Total Lines.
    int GetNTotalLine() const;

    //! Return The number of interface faces.
    //! Interface faces: faces between different zones,
    //! such as faces between zones in different unstructured grid partitions,
    //! and faces between zones in different structured grid blocks.
    int GetNIFace() const;

    //! Return the interface topology information of interface faces.
    InterfaceInfo * GetInterfaceInfo() const;

    //! Return the interpoint topology information of interface points.
    //! @param[out] interpointInformation    the object of class InterpointInformation.
    InterpointInformation * GetInterpointInfo() const;

    //! Return the Interface fields information.
    //! Stores the field variables that need to communicate between different gird of zones.
    InterfaceFields * GetInterfaceFields();

    //! Return the Interpoint fields information.
    //! @param[out] interpointFields    Interpoint fields information.
    InterpointFields * GetInterpointFields();

    //! Return the Coarser grid of the current grid.
    Grid * GetCoarseGrid() const;

    //! Return the Finer grid of the current grid.
    Grid * GetFineGrid() const;

    //! Copy any type of control parameters from give pointer, such as int, RDouble, string.
    void CopyPara(Data_Param *gPara);

    //! Copy any type of flow field variables from give pointer, such as velocity, time, temperature.
    void CopyField(Data_Field *gField);

    //! Return the Grid index.
    GridID * GetGridID() const;

    //! Return the Global Zone Grid index.
    int GetZoneID() const;

    //! Return the Local Zone Grid index.
    int GetZoneLocalID() const;

    //! Return the overset or overlapping Grid information.
    OversetInfoProxy * GetOversetInfoProxy() const;

    //! Return the overset information for unstructured overset cell interpolation.
    OversetInformationProxy * GetOversetInformationProxy() const;

    //! Return iBlock.
    int GetIBlock();

    //! Return the probes number of each zone grid.
    int GetZoneProbesNumber() const;

    //! Return the probes global index of each zone grid.
    vector<int> GetZoneProbesGlobalID() const;
    //! Return the probes line index of each zone grid.
    vector<int> GetZoneProbesLineID() const;
    //! Return the probes surface index of each zone grid.
    vector<int> GetZoneProbesSurfaceID() const;

    //! Return the probes coordinates of each zone grid.
    vector<vector<RDouble> > GetZoneProbesCoordinates() const;

    //! Is this grid is the finest grid, in multi-grid method.
    bool IsFinestGrid();

    //! Is this grid is the coarsest grid, in multi-grid method.
    bool IsCoarsestGrid();

    //! Return the initial finest level grid, in multi-grid method.
    LIB_EXPORT Grid * GetFinestGrid();
    
    //! Back up the grid at the last time step.
    //! It is used in unsteady dynamic grid method.
    LIB_EXPORT virtual void BackUpOldGrid();

    //! Return the grid at the last time step.
    //! It is used in unsteady dynamic grid method.
    SimpleGrid * GetOldGrid();

    void SetOldGrid(SimpleGrid *oldGrid);

    //! Register an interface field (variable) into buffer which is loaded on grid, used in parallel communication.
    //! The interface field variable buffer is set to be ZERO.
    //! @param[in] name    variable name of the interface field that need to be communicated.
    //! @param[in] type    variable type of the interface field.
    //! -# PHINT   : Integer type.
    //! -# PHFLOAT : Float type.
    //! -# PHDouble: Double type.
    //! @param[in] dimension    variable dimension.
    //! -# 1: 1-D array.
    //! -# 2: 2-D array.
    //! -# n: n-D array.
    LIB_EXPORT void RegisterInterfaceField(const string &name, const int type, const int dimesion, const int solverID);

    LIB_EXPORT void ReleaseInterfaceField();

    LIB_EXPORT void ReleaseInterfaceField(string varName);
    
    //! Register an interpoint field (variable) into buffer which is loaded on grid, used in parallel communication.
    //! The interpoint field variable buffer is set to be ZERO.
    //! @param[in] name    variable name of the interface field that need to be communicated.
    //! @param[in] type    variable type of the interface field.
    //! -# PHINT   : Integer type.
    //! -# PHFLOAT : Float type.
    //! -# PHDouble: Double type.
    //! @param[in] dimension    variable dimension.
    //! -# 1: 1-D array.
    //! -# 2: 2-D array.
    //! -# n: n-D array.
    LIB_EXPORT void RegisterInterpointField(const string &name, const int type, const int dimesion, const int solverID);

    //! Allocate storage for overset field of inner cells and donor cells.
    LIB_EXPORT virtual void RegisterOversetField(const string &name, const int type, const int dimesion);
    
    //! This subroutine calculates the new set of coordinates,
    //! after a rotation by 'dangle' and a translation by 'tcoord'.
    //! @param[in]  coord0    Original coordinates.
    //! @param[out] xyzref    New coordinates, after rotation and translation.
    LIB_EXPORT void RotTransVel(RDouble *xyzref, RDouble *dangle, RDouble *tcoord, RDouble **coord0);

    //! This subroutine calculates the grid velocity at the cell vertexes.
    //! Compute the components of the rotation vector.
    //! Omega (psip, thetap, phip) in the frame (k0, jh, i) where
    //! k0 = -sin theta * i + cos theta sin phi * j + cos theta cos phi * k
    //! jh = cos phi * j - sin phi * k.
    //! Omega(p,q,r) in the frame (i,j,k) angle0(1,.) = (psi,theta,phi)
    //! psi is not set since it is not needed.
    //! @param[in] angle0    Euler's angles and their time derivatives.
    //! @param[in] coord0    position of the reference point abd their time derivatives.
    LIB_EXPORT void ComputeGDVel(RDouble *xt, RDouble *yt, RDouble *zt, RDouble **angle0, RDouble **coord0);
    
    //! Turn the Y axis to be up, if the original Z axis is up.
    LIB_EXPORT void RotateAboutAxis();

    //! Read grid from file.
    LIB_EXPORT virtual void ReadGrid(fstream &file) = 0;

    //! Write grid to file.
    LIB_EXPORT virtual void WriteGrid(fstream &file) = 0;

    //! Compress the grid data structure into data container.
    LIB_EXPORT virtual void Encode(DataContainer *cdata, int flag) = 0;

    //! Decompress the grid data structure from data container.
    LIB_EXPORT virtual void Decode(DataContainer *cdata, int flag) = 0;

    //! 
    LIB_EXPORT virtual void ReviseNeighborZoneIndex(int zoneStart) = 0;

    //! Update the grid coordinates from data container.
    LIB_EXPORT virtual void UpdateCoordinate(DataContainer *cdata) = 0;

    //! Initialize the grid after the object has been created.
    LIB_EXPORT virtual void InitGrid(GridID *index, int level, int dim, int type);

    //! Static the mesh skewness only for unstructured grid.
    //! Skewness, 0 represents bad mesh quality; 1 represents satisfactory quality.
    LIB_EXPORT virtual RDouble SkewnessSummary(ActionKey *actkey = 0) = 0;

    //! Initialize over-lapping grid information.
    LIB_EXPORT virtual void AllocateOversetGrid() = 0;

    //! Compute the grid metrics, such as volume, face area, cell center, face center, etc.
    LIB_EXPORT virtual void ComputeMetrics(ActionKey *actkey = 0) = 0;

    //! Compute node bctype in each zone.
    LIB_EXPORT virtual void ComputeNodeBCType() {};

    //! Get the grid cell center data.
    LIB_EXPORT virtual void GetCellCenter(ActionKey *actkey) = 0;

    //! Translate the grid cell center data.
    LIB_EXPORT virtual void TranslateCellCenter(ActionKey *actkey) = 0;

    //! Compute and return the number of wall faces.
    LIB_EXPORT virtual int GetNumberOfWallCell() = 0;

    //! Compute node weight for Least Square Approach.
    LIB_EXPORT virtual void ComputeWeight() = 0;

    //! Compute the minimum and maximum distance between two grid nodes in whole computational zones.
    LIB_EXPORT virtual void GetMinMaxDS(RDouble &dismin, RDouble &dismax) = 0;

    //! Compute the minimum and maximum distance between two grid nodes on Boundary.
    LIB_EXPORT virtual void GetBoundaryMinMax(RDouble *pmin, RDouble *pmax) = 0;

    //! Compute the coarse grid of the current grid,
    //! using agglomeration method, for multi-grid method.
    LIB_EXPORT virtual void CoarseGrids(int nlmax) = 0;

    //! Initialize the BC information of structured grid,
    //! such as left and right cell of the BC faces.
    LIB_EXPORT virtual void ProcessBCInfo() = 0;
    LIB_EXPORT virtual void SetLnkInfo() = 0;

    //! Change the BC type, if necessary.
    LIB_EXPORT virtual void ChangeBCType(int from_bctype, int to_bctype) = 0;

    //! Build the interface topology link between different zones.
    //! Such as faces between zones in different unstructured grid partitions,
    //! and faces between zones in different structured grid blocks.
    LIB_EXPORT virtual void BuildGridLink(LinkStruct *link) = 0;
    
    //! Update the cell volume between two time step.
    LIB_EXPORT virtual void UpdateVolold() = 0;

    //! Compute the difference flow variable normal between two time step.
    LIB_EXPORT virtual void CVGNorm(FieldProxy *q_1_proxy, FieldProxy *q_0_proxy, int neqn, RDouble &norm) = 0;
    
    //! Compute the grid velocity.
    LIB_EXPORT virtual void GridSurfaceVelocity(RDouble *xt, RDouble *yt, RDouble *zt) = 0;

    LIB_EXPORT void GridVerticeVelocity(RDouble *xt, RDouble *yt, RDouble *zt);

    LIB_EXPORT virtual void ComputeGridFaceVelocity() = 0;

    //! Compute the minimum distance of the grid edges.
    LIB_EXPORT virtual RDouble CalMinEdgeLength() = 0;

    //! Compute the maximum distance of the grid edges.
    LIB_EXPORT virtual RDouble CalMaxEdgeLength() = 0;

    //! Write face boundary condition to file stream.
    LIB_EXPORT virtual void WriteFaceBC(fstream &file) = 0;

    //! Initialize the dynamic/moving grid information.
    LIB_EXPORT void InitMovingGrids();

    //! Check grid validity.
    LIB_EXPORT void GridValidityCheck();

    //! 
    virtual void InitWallDist() = 0;

    //! 
    virtual void ReadWalldist(RDouble *wallDistIn) = 0;

    //!
    virtual void ReadNearestwallfacenormalx(RDouble *nearestwallfaceNormalxIn) = 0;
    virtual void ReadNearestwallfacenormaly(RDouble *nearestwallfaceNormalyIn) = 0;
    virtual void ReadNearestwallfacenormalz(RDouble *nearestwallfaceNormalzIn) = 0;

    //! 
    virtual void DumpWalldist(ActionKey *actkey) = 0;

    virtual void DumpWallFaceCenter(ActionKey *actkey) = 0;

    //! Update a data pointer into the Running-Data-Base.
    //! @param[in] name    Name of the data pointer, type of string.
    //! @param[in] data    The data pointer that need to be stored in Running-Data-Base.
    void UpdateDataPtr(const string &name, void *data);

    //! Update a data into the Running-Data-Base.
    //! @param[in] name    Name of the data, type of string.
    //! @param[in] data    The data that need to be stored in Running-Data-Base.
    //! @param[in] size    The size of the data.
    //! @param[in] type    The data type.
    //! -# PHINT   : integer type.
    //! -# PHFLOAT : RFloat type.
    //! -# PHDOUBLE: RDouble type.
    //! -# PHSTRING: string type.
    //! -# PHBOOL  : bool type.
    void UpdateData(const string &name, void *data, int type, int size);

    //! Get data pointer with name of 'name' from Running-Data-Base.
    void * GetDataPtr(const string &name) const;
   
    //! Get data with name of 'name' from Running-Data-Base.
    //! @param[out] data    The data that need to be got from Running-Data-Base.
    //! @param[in]  name    Name of the data, type of string.
    //! @param[in]  size    The size of the data.
    //! @param[in]  type    The data type.
    //! -# PHINT   : integer type.
    //! -# PHFLOAT : RFloat type.
    //! -# PHDOUBLE: RDouble type.
    //! -# PHSTRING: string type.
    //! -# PHBOOL  : bool type.
    void GetData(const string &name, void *data, int type, int size);

    //! Delete a data with name of 'name' from Running-Data-Base.
    void DeleteDataPtr(const string &name) const;

    void SetFileIndex(int fileIndex);
    int  GetFileIndex();

    set < Data_IndexData > * GetDataSet();

    //! Action mechanism, these action is not allowed to be calling by library.
    virtual void Action(ActionKey *actkey) = 0;
    virtual void TranslateAction(ActionKey *actkey) = 0;
    virtual streamsize TranslateActionLength(ActionKey *actkey) = 0;
    virtual void ReSetBoundaryConditionByGlobalBC() {};
    virtual void InitVariableWallTemperature() {};
    //! Allocate memory for walldist.
    virtual void AllocateWalldist() = 0;
    virtual void GetBcNamePair(set< pair<int, string> > &bcNamePairSet) = 0;
    virtual void ComputeCellBoundaryType();

    virtual int  GetOrdinaryGridIndex() { return -1; };

    virtual RDouble  GetPartialCFL() { return -1.0; };

#pragma warning(disable:4100)
    virtual void SetPartialCFL(RDouble PartialCFL){ 
#pragma warning(default:4100) 
    };

    //! To compress the specified array to the the DataContainer which will be sent by parallel communication .
    //! @param[out]: dataContainer denotes the compress database which will be sent.
    //! @param[in ]: fieldName denotes the name of flow field to communicate.
    //! @param[in ]: neighborZoneIndex denotes the grid index of target zone.
    //! @param[in ]: nEquation denotes the number of variable to send.
    virtual void CompressSpecificArrayToInterface(DataContainer *&dataContainer, const string &fieldName, const int &neighborZoneIndex, const int &nEquation) = 0;

    //! To decompress the specified array from the the DataContainer received from parallel communication.
    //! @param[in ]: dataContainer denotes the compress database which was received.
    //! @param[in ]: fieldName denotes the name of flow field to communicate.
    //! @param[in ]: neighborZoneIndex denotes the grid index of source zone.
    //! @param[in ]: nEquation denotes the number of received variable.
    virtual void DecompressArrayFromInterface(DataContainer *&dataContainer, const string &fieldName, const int &neighborZoneIndex, const int &nEquation) = 0;

protected:
    virtual void AllocateMetrics(ActionKey *actkey) = 0;
    void FreeInterfaceField();
    void RotateTranslation();
};

#include "Geo_Grid.hxx"

void CheckGrid(Grid **grids, int nBlocks);

//!
LIB_EXPORT void Create_Link_Info(VInt &pindex, int iZone, int iFace, uint_t &fcount, set< DataStruct_Sort< VInt > > &face_list, VVInt &zoneid, VVInt &faceid, LinkStruct *link);
//!
LIB_EXPORT void GetFaceCoorList(VInt &index, int nlist, RDouble *xlist, RDouble *ylist, RDouble *zlist, RDouble *x, RDouble *y, RDouble *z);
//!
LIB_EXPORT void GetFaceCoorList(int il, int jl, int kl, int ir, int jr, int kr, int ijk[4][3], int nlist,
                                RDouble *xlist, RDouble *ylist, RDouble *zlist, RDouble3D &x, RDouble3D &y, RDouble3D &z);
//!
LIB_EXPORT void GetCoorIndexList(DataStruct_AdtTree<int, RDouble> *adt_tree, RDouble &diff_tolerance, int &pcount, RDouble *xlist, RDouble *ylist, RDouble *zlist, int nlist, vector<int> &pindex);
//!
LIB_EXPORT int  GetCoorIndex(DataStruct_AdtTree<int, RDouble> *adt_tree, RDouble *coor, RDouble *minwin, RDouble *maxwin, int &count, int &index);
//!
LIB_EXPORT void ShiftMinMaxBox(RDouble *pmin, RDouble *pmax, RDouble tol);

//! Get cell type on the basis of node number.
int GetCellType(int dimension, int nodeNumberOfCell);

//! Get face type on the basis of node number.
int GetFaceType(int dimension, int nodeNumberOfFace);

//! Some special functions for HYBRID solver.
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//! Special functions for HYBRID solver.
bool IsConvertGridToMixGrid();

LIB_EXPORT void CommunicateSpecificArray(int level, string arrayName, int nEquation);

extern ZoneConnectivity *zoneConnectivity;
extern ZoneConnectivityForPoint *zoneConnectivityForPoint;
}
