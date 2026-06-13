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
//! @file      Geo_Interface.h
//! @brief     Interface information for multi-zone grid.
//! @author    Bell, He Xin.

#pragma once
#include "Precision.h"
#include <vector>
#include "LIB_Macro.h"
#include "Data_ParamFieldSuite.h"
#include "GlobalDataBase.h"

namespace PHSPACE
{
inline int GetNumberOfGhostCellLayers()
{
    //! Warning: the layer number should be judged by grid type.
    //! however, it will lead to RESTART file (flow.dat/turb.dat)
    //! different with the previous, using 2 layers temporary.
    //! the bellowing judge should be used in future.
    int nLayers = 2;

    //! Add one layer ghost cell for high order method.
    int isWennScheme = PHSPACE::GlobalDataBase::GetIntParaFromDB("isWennScheme");
    if (isWennScheme == 1)
    {
        nLayers = 3;
    }

    return nLayers;

    //! The bellowing judge should be used in future.
    //const int UNSTRUCTGRID = 0;
    //const int STRUCTGRID   = 1;
    //const int MIXGRID      = 2;

    //! Note that sys_gridtype is undefined! Error may occur here!
    //int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    //if (systemGridType == UNSTRUCTGRID)
    //{
    //    return 1;
    //}
    //else
    //{
    //    return 2;
    //}
}

class InterfacePatch
{
private:
    //! The index of the current zone.
    int zoneid;

    //! The index of the neighbor zone.
    int zoneIndexOfNeighbor;

    //! The number of patched-interface in the neighbor zone which linked to current zone.
    int faceNumberOfPatch;

    //! Interface number in the CURRENT zone, using in receiving.
    int *faceIndexForRecv;
    
    //! Interface number in the TARGET zone, using in sending.
    int *faceIndexForSend;

    //! GMRESParallel cell index in the local zone
    int *cellIndexOfLocal;

    //! GMRESParallel ghost cell index in the local zone
    int *cellIndexOfLocalGhost;

    //! GMRESParallel cell index in the neighbor zone
    int *cellIndexOfNeighbor;

private:
    Data_ParamFieldSuite *dataStorageSend;
    Data_ParamFieldSuite *dataStorageRecv;

public:
    InterfacePatch();
    InterfacePatch(const InterfacePatch &interPatchchIN);
    InterfacePatch & operator = (const InterfacePatch  &interPatchchIN);
    ~InterfacePatch();

public:
    int * GetFaceIndexForRecv() const;
    void SetFaceIndexForRecv(int *faceIndexForRecv);

    int * GetFaceIndexForSend() const;
    void SetFaceIndexForSend(int *faceIndexForSend);

    int  GetNumberOfFace() const;
    void SetNumberOfFace(int faceNumberOfPatch);

    int  GetZoneIndexOfNeighbor() const;
    void SetZoneIndexOfNeighbor(int zoneIndexOfNeighbor);

    //! GMRESParallel
    int *GetCellIndexOfLocal() const;
    void SetCellIndexOfLocal(int *cellIndexLocal);

    int *GetCellIndexOfLocalGhost() const;
    void SetCellIndexOfLocalGhost(int *cellIndexLocalGhost);

    int *GetCellIndexOfNeighbor() const;
    void SetCellIndexOfNeighbor(int *cellIndexNeighbor);

    int GetZoneID() const;

    void AllocateData();
};

class Grid;

class InterfaceInfo
{
private:
    //! Number of interfaces.
    int nIFace;

    //! Number of neighbor zones.
    int numberOfNeighbor;

    //! Neighbor zone index of the interface in current zone.
    int *interFace2ZoneID;

    //! Cell index to neighbor zone of the interface in current zone.
    int *interFace2CellID;

    //! Interface index in neighbor zone of the interface in current zone.
    //! Be careful, in grid partition procedure, it stores the neighbor local cell index firstly,
    //! and then, it was covered by neighbor local interface index in interface matching stage.
    int *interFace2InterFaceID;

    //! Boundary face index in current zone of the interface in current zone.
    int *interFace2BoundaryFace;

    //! Interface direction.
    //! -#  1: If left cell of the interface exist.
    //! -# -1: If right cell of the interface exist.
    int *interFaceDirection;

    //! Parent grid of the current grid, using in multi-grid method.
    Grid *parent;

    //! Data storage for sending.
    Data_ParamFieldSuite *dsend;

    //! Data storage for receiving.
    Data_ParamFieldSuite *drecv;

    vector <InterfacePatch *> interFacePatch;

    //! GMRESParallel global cell index of the ghost cell for interface bc
    int *globalNeighborCellIndex;

    int *localInterfaceCellIndex;

    int *localInterfacePhysicalCellIndex;

public:
    InterfaceInfo();
    InterfaceInfo(const InterfaceInfo &infoIn);
    InterfaceInfo(int nIFace, Grid *parent = 0);
    ~InterfaceInfo();

public:
    int * GetInterFace2ZoneID() const;
    int * GetInterFace2CellID() const;
    int * GetInterFace2InterFaceID() const;
    int * GetInterFace2BoundaryFace() const;
    int * GetInterFaceDirection() const;

    //! Get the neighbor zone index of ineighbor, ineighbor among in [0, numberOfNeighbor).
    int GetZoneIndexOfNeighbor(int ineighbor) const;

    int * GetFaceIndexForRecv(int ineighbor) const;
    int * GetFaceIndexForSend(int ineighbor) const;

    Data_ParamFieldSuite * GetSendDataStorage(int istore) const;
    Data_ParamFieldSuite * GetRecvDataStorage(int istore) const;

    void SetNIFace(int nIFace);
    int GetNIFace() const;

    void SetNumberOfNeighbor(int numberOfNeighbor);
    int GetNumberOfNeighbor() const;

    int FindIthNeighbor(int zone);

    int * ComputeFaceIndexForSend(int ineighbor);

    void FillFaceIndexForSend(int zone, int *faceIndexForSend_in);

    vector <InterfacePatch *> GetInterfacePatch() const;
    InterfacePatch * GetInterfacePatch(int ineighbor) const;

    void FindNeighbors();

    void InitNeighbors();
    void InitNeighbors(int ineighbor);

    void AllocateNeighborInfo();
    void AllocateSendRecv(int ineighbor);

    int GetNIFaceOfNeighbor(int ineighbor) const;

    void ComputeFaceIndexForRecv(int ineighbor);

    int ComputeNIFaceOfNeighbor(int ineighbor);

    int ComputeNumberOfNeighbor(int *zoneflag);

    //! Resize the data according to the new number of inter faces.
    void ReSize(int numberOfNewInterfaces);

    void SetIsNeighborInterfaceFound(bool *data);
    bool * GetIsNeighborInterfaceFound();

    //! GMRESParallel
    int *GetCellIndexOfNeighbor(int ineighbor) const;

    int *GetGlobalNeighborCellIndex() const;

    void SetGlobalNeighborCellIndex(int *globalNeighborCellIndex);

    int *GetLocalInterfaceCellIndex() const;

    void SetLocalInterfaceCellIndex(int *localInterfaceCellIndex);

    int *GetLocalInterfacePhysicalCellIndex() const;

    void SetLocalInterfacePhysicalCellIndex(int *localInterfacePhysicalCellIndex);

    int MatchedGlobalNeighborCellIndex(int cell);

    int MatchedLocalPhysicalCellIndex(int cell);

private:
    //! Is the neighbor interface has been found.
    //! this is used during partition only.
    bool *isNeighborInterfaceFound;
};

//! Interface variable fields
//! Store the interface variable for multi-zone grid, using for Parallel communication.
class InterfaceFields
{
private:
    //! List of interface field variable names.
    vector< string > variableNames;

    //! List of interface field variable types.
    vector< int > variableTypes;

    //! List of interface field variable dimensions.
    //! eg. 1D/2D array, dimension is equal to '1' or '2'.
    vector< int > variableDimensions;

    //! Interface information.
    InterfaceInfo *interfaceInfo;

    //! Solver ID.
    vector< int > solverIndex;

public:
    InterfaceFields(InterfaceInfo *dataIn);
    ~InterfaceFields();

public:
    uint_t Size();

    int GetDim(const int iData);

    string & GetName(const int iData);

    int GetSolverIndex(const int iData);

    void RegisterField(const string &name, const int type, const int dimesion, const int solverID);

    void RemoveAllField();

    void RemoveAnVariable(string varName);
    
    void FreeInterfaceVar();

    template< typename T >
    T ** GetSendVariable(const string &name, const int iGhost);

    template< typename T >
    T ** GetRecvVariable(const string &name, const int iGhost);

private:
    void AllocateInterfaceVar(const string &name, const int type, const int dimesion, Data_ParamFieldSuite *datastore);
    void FreeInterfaceVar(Data_ParamFieldSuite *datastore);
    void RemoveInterfaceVar(const string &name);
};

class InterfaceDataProxy
{
private:
    vector<RDouble **> vdata;
    vector<int> vdim;

public:
    InterfaceDataProxy();
    ~InterfaceDataProxy();
    vector<string> dataName;

public:
    vector<RDouble **> & GetVectorData();
    vector<int> & GetVectorDimension();
    vector<string> & GetVectorName();
};

#include "Geo_Interface.hxx"
}
