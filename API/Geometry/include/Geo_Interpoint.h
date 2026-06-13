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
//! @file      Geo_Interpoint.h
//! @brief     Interpoint information for multi-zone grid.
//! @author    Wan Yunbo.
//! @email     wanyb11@163.com
//! @date      2019-11-30

#pragma once
#include <vector>
#include "Precision.h"
#include "LIB_Macro.h"
#include "Data_ParamFieldSuite.h"
namespace PHSPACE
{
//! @brief patched-interpoint class for one neighor zone.
class InterpointPatch
{
private:
    //! The index of the current zone.
    int zoneID;

    //! The index of the neighbor zone.
    int zoneIndexOfNeighbor;

    //! The number of patched-interpoint in the neighbor zone which linked to current zone.
    int pointNumberOfPatch;

    //! Interpoint index in the CURRENT zone, using in receiving.
    int *pointIndexForReceive;
    
    //! Interpoint index in the TARGET zone, using in sending.
    int *pointIndexForSend;

private:
    //! Store the data for sending to other zones.
    Data_ParamFieldSuite *dataStorageSend;

    //! Store the data for receiving from other zones.
    Data_ParamFieldSuite *dataStorageReceive;

public:
    //! Construct function for class InterpointPatch.
    InterpointPatch();
    //! Copy construct function for class InterpointPatch.
    InterpointPatch(const InterpointPatch &interPatchIn);
    //! Destructor, free memory.
    ~InterpointPatch();
public:
    //! Get the interpoint index in the target zone.
    //! @param[out] pointIndexForReceive      interpoint index for receiving in current zone.
    int * GetPointIndexForReceive() const;

    //! Set the interpoint index in the target zone.
    //! @param[in] pointIndexForReceive      interpoint index for receiving in current zone.
    void SetPointIndexForReceive(int *pointIndexForReceive);

    //! Get the interpoint number in the target zone.
    //! @param[out] pointIndexForSend      interpoint index for sending in current zone.
    int * GetPointIndexForSend() const;

    //! Set the interpoint number in the target zone.
    //! @param[in] pointIndexForSend      interpoint index for sending in current zone.
    void SetPointIndexForSend(int *pointIndexForSend);

    //! Get the number of patched-interpoint in the neighbor zone.
    //! @param[out] pointNumberOfPatch      The number of patched-interpoint in the neighbor zone.
    int GetNumberOfPoint() const;

    //! Get the number of patched-interpoint in the neighbor zone.
    //! @param[in] pointNumberOfPatch      The number of patched-interpoint in the neighbor zone.
    void SetNumberOfPoint(int pointNumberOfPatch);

    //! Get the zone index of the neighbor zone.
    //! @param[out] zoneIndexOfNeighbor      the zone index of the neighbor zone.
    int  GetZoneIndexOfNeighbor() const;

    //! Set the zone index of the neighbor zone.
    //! @param[in] zoneIndexOfNeighbor      the zone index of the neighbor zone.
    void SetZoneIndexOfNeighbor(int zoneIndexOfNeighbor);

    //! Get the zone index of the current zone.
    //! @param[out] zoneID      The index of the current zone.
    int GetZoneID() const;

    //! Allocate the memory for sending and receiving data.
    void AllocateData();
};

class Grid;

//! @brief InterpointInformation: inter-point information, similar with the interface information,
//!        which is used to store the inter-point connection information between neighbor zones.
class InterpointInformation
{
private:
    //! Number of interpoints.
    int numberOfInterpoints;

    //! Number of neighbor zones.
    int numberOfNeighbor;

    //! Neighbor zone index of the interpoint in current zone.
    int *interPoint2ZoneID;

    //! Interpoint index in neighbor zone of the interpoint in current zone.
    int *interPoint2InterPointID;

    //! Global point index in current zone of the interpoint in current zone.
    int *interPoint2GlobalPoint;

    //! Total cell numbers of the interpoint in computing node value, is not node2cell.
    int *cellNumberOfInterPoint;

    //! Total zones of the interpoint belongs to.
    int *totalZonesOfInterPoint;

    //! Label of the interpoint,if =1 ,the InterPoint is needed in computing node value.if =0,
    //! it is not needed in computing node value.
    int *labelOfInterPoint;

    //! Parent grid of the current grid, using in multi-grid method.
    Grid *parent;

    //! Data storage for sending.
    Data_ParamFieldSuite *dataForSend;

    //! Data storage for receiving.
    Data_ParamFieldSuite *dataForReceive;

    //! The number of patched-interpoint in the neighbor zone which linked to current zone.
    InterpointPatch **interPointPatch;
public:
    //! Construct function for class InterpointInformation.
    InterpointInformation();

    //! Copy Construct function for class InterpointInfo.
    //! @param[out] infoIn           an object of class InterpointInfo.
    InterpointInformation(const InterpointInformation &infoIn);

    //! Copy Construct function for class InterpointInformation.
    //! @param[in] numberOfInterpoints           Number of interpoints.
    //! @param[in] parent           Parent grid of the current grid.
    InterpointInformation(int numberOfInterpoints, Grid *parent = 0);

    //! Destructor, free memory.
    ~InterpointInformation();

public:

    //! Get the number of interpoints.
    //! @param[out] interPoint2ZoneID      Neighbor zone index of the interpoint.
    int * GetInterPoint2ZoneID() const;

    //! Get interpoint index in neighbor zone.
    //! @param[out] interPoint2InterPointID      interpoint index in neighbor zone.
    int * GetInterPoint2InterPointID() const;

    //! Get global point index in current zone.
    //! @param[out] interPoint2GlobalPoint       Global point index in current zone.
    int * GetInterPoint2GlobalPoint() const;

    //! Get the total cell numbers of the interpoint.
    //! @param[out] cellNumberOfInterPoint      total cell numbers of the interpoint.
    int * GetCellNumberOfInterPoint() const;

    //! Get the total zones of the interpoint.
    //! @param[out] totalZonesOfInterPoint       Total zones of the interpoint.
    int * GetTotalZonesOfInterPoint() const;

    //! Get the label of the interpoint.
    //! @param[out] labelOfInterPoint          Label of the interpoint.
    int * GetLabelOfInterPoint() const;

    //! Get the neighbor zone index of iNeighbor, iNeighbor among in [0, numberOfNeighbor].
    //! @param[in] iNeighbor            the i-th neighbor for interpoint of this zone.
    //! @param[out]           the zone index of the i-th neighbor.
    int GetZoneIndexOfNeighbor(int iNeighbor) const;

    //! Get the point index of iNeighbor for receiving.
    //! @param[in] iNeighbor            the i-th neighbor for interpoint of this zone.
    //! @param[out]           the point index of iNeighbor for receiving.
    int * GetPointIndexForReceive(int iNeighbor) const;

    //! Get the point index of iNeighbor for sending.
    //! @param[in] iNeighbor            the i-th neighbor for interpoint of this zone.
    //! @param[out]           the point index of iNeighbor for sending.
    int * GetPointIndexForSend(int iNeighbor) const;

    //! Get the send data storage of the i-th data.
    //! @param[in] istore         the ID the sending data storage.
    //! @param[out]dataForSend[istore]     the send data storage of the i-th data.
    Data_ParamFieldSuite * GetSendDataStorage(int istore) const;

    //! Get the receive data storage of the i-th data.
    //! @param[in] istore         the ID the receiving data storage.
    //! @param[out] dataForSend[istore]     the receive data storage of the i-th data.
    Data_ParamFieldSuite * GetReceiveDataStorage(int istore) const;

    //! Set the number of interpoints.
    //! @param[in] numberOfInterpoints           Number of interpoints.
    void SetNumberOfInterpoints(int numberOfInterpoints);

    //! Get the number of interpoints.
    //! @param[out] numberOfInterpoints           Number of interpoints.
    int GetNumberOfInterpoints() const;

    //! Set the number of Neighbor for interpoint.
    //! @param[in] numberOfNeighbor          Number of neighbor zones.
    void SetNumberOfNeighbor(int numberOfNeighbor);

    //! Get the number of Neighbor for interpoint.
    //! @param[out] numberOfNeighbor          Number of neighbor zones.
    int  GetNumberOfNeighbor() const;

    //! get the index of Neighbor zone.
    //! @param[in] zone          zone Index.
    //! @param[out]            the i-th index.
    int FindIthNeighbor(int zone);

    //! Compute the point index of iNeighbor for sending.
    //! @param[in] iNeighbor            the i-th neighbor for interpoint of this zone.
    //! @param[out]           the point index of iNeighbor for sending.
    int * ComputePointIndexForSend(int iNeighbor);

    //! Fill the point index of iNeighbor for sending to zone.
    //! @param[in] zone            the zone Index for sending to.
    //! @param[in] pointIndexForSendIn         point index for sending to zone.
    void FillPointIndexForSend(int zone, int *pointIndexForSendIn);

    //! Get the number of patched-interpoint.
    //! @param[out] interPointPatch        the number of patched-interpoint.
    InterpointPatch ** GetInterpointPatch() const;

    //! Get the patched-interpoint of the ith neighbor zone.
    //! @param[in] iNeighbor        the i-th neighbor zone.
    //! @param[out] interPointPatch        the patched-interpoint of the i-th neighbor zone.
    InterpointPatch * GetInterpointPatch(int iNeighbor) const;

    //! Find the neighbor zone for interpoint.
    void FindNeighbors();

    //! Init the neighbor zone for interpoint.
    void InitNeighbors();

    //! Init the i-th neighbor zone for interpoint.
    //! @param[in] iNeighbor        the i-th neighbor zone.
    void InitNeighbors(int iNeighbor);

    //! Allocate the memory for the neighbors for interpoint.
    void AllocateNeighborInfo();

    //! Allocate the memory for the neighbors sending or receive of the interpoint.
    //! @param[in] iNeighbor        the i-th neighbor zone.
    void AllocateSendReceive(int iNeighbor);

    //! Get the number of the interpoint with the i-th neighbor.
    //! @param[in] iNeighbor        the i-th neighbor zone.
    int GetNumberOfInterpointsForNeighbor(int iNeighbor) const;

    //! compute the point index for receive interpoint in current zone.
    //! @param[in] iNeighbor        the i-th neighbor zone.
    void ComputePointIndexForReceive(int iNeighbor);

    //! Compute the number of the interpoint for the i-th neighbor.
    //! @param[in] iNeighbor        the i-th neighbor zone.
    int ComputeNumberOfInterpointsForNeighbor(int iNeighbor);

    //! Compute the number of the neighbors for this zone.
    //! @param[in] zoneFlag        array to store the neighbor zones.
    int ComputeNumberOfNeighbor(int *zoneFlag);

    //! Resize the number of new interpoints according to the new number of inter points.
    //! @param[in] numberOfNewInterpoints        the number of new interpoints in this zone.
    void ReSize(int numberOfNewInterpoints);

    //! Set the neighbor of this interpoint found or not.
    //! @param[in] data        array to store true or false.
    void SetIsNeighborInterpointFound(bool *data);

    //! Get the neighbor of this interpoint found or not.
    //! @param[out] isNeighborInterpointFound         Is the neighbor interpoint has been found.
    bool * GetIsNeighborInterpointFound();

private:
    //! Is the neighbor interpoint has been found.
    //! This is used during partition only.
    bool *isNeighborInterpointFound;
};

//! Interpoint variable fields
//! Store the interpoint variable for multi-zone grid, using for Parallel communication.
class InterpointFields
{
private:
    //! List of interpoint field variable names.
    vector<string> variableNames;

    //! List of interpoint field variable types.
    vector<int> variableTypes;

    //! List of interpoint field variable dimensions.
    //! eg. 1D/2D array, dimension is equal to '1' or '2'.
    vector<int> variableDimensions;

    //! Interpoint information.
    InterpointInformation *interpointInformation;

    //! Solver ID.
    vector<int> solverIndex;

public:

    //! Construct function.
    //! @param[in] interpointInformationIn       an object of class InterpointInformation.
    InterpointFields(InterpointInformation *interpointInformationIn);

    //! Destructor, free memory.
    ~InterpointFields();

public:
    //! Get the number of variables.
    //! @param[out]         the number of variables.
    uint_t Size();

    //! Get the dimension  of the i-th data.
    //! @param[in]iData         data index.
    //! @param[out]         the dimension  of the i-th data.
    int GetDimension(const int iData);

    //! Get the name  of the i-th data.
    //! @param[in]iData         data index.
    //! @param[out]         the name  of the i-th data.
    string & GetName(const int iData);

    //! Get the solver index  of the i-th data.
    //! @param[in]iData         data index.
    //! @param[out]         the solver index  of the i-th data.
    int GetSolverIndex(const int iData);

    //! Register the data fields according to the name, type and dimensional.
    //! @param[in] name         data fields name.
    //! @param[in] type         data fields type.
    //! @param[in] dimesion         data fields dimesion.
    //! @param[in] solverID         data fields solverID.
    void RegisterField(const string &name, const int type, const int dimesion, const int solverID);

    //! Free the store for interpoint.
    void FreeInterpointVariable();

    //! Get the sending variable according to the name and the ghost layers.
    //! @param[in] name         data fields name.
    //! @param[in] IDOfGhostLayer         the ID of ghost layers..
    //! @param[out]              the sending variable.
    template<typename T>
    T ** GetSendVariable(const string &name, const int IDOfGhostLayer);

    //! Get the sending variable according to the name and the ghost layers.
    //! @param[in] name         data fields name.
    //! @param[in] IDOfGhostLayer         the ID of ghost layers..
    //! @param[out]              the receiving variable.
    template<typename T>
    T ** GetReceiveVariable(const string &name, const int IDOfGhostLayer);

private:

    //! Allocate the memory for variable according to the name type and dimesion.
    //! @param[in] name         data fields name.
    //! @param[in] type         data fields type.
    //! @param[in] dimesion         data fields dimesion.
    //! @param[out] datastore        the data store for all variables needing communicating.
    void AllocateInterpointVariable(const string &name, const int type, const int dimesion, Data_ParamFieldSuite *datastore);

    //! Free the memory for all variables.
    //! @param[in] datastore        the data store for all variables needing communicating.
    void FreeInterpointVariable(Data_ParamFieldSuite *datastore);
};

//! @brief store the interpoint value.\n
class InterpointDataProxy
{
private:
    //! The vector to store the interpoint data for sending or receiving.
    vector<RDouble **> vdata;

    //! The vector to store the interpoint data dimension for sending or receiving.
    vector<int> vdim;
public:
    //! Construct function for class InterpointDataProxy.
    InterpointDataProxy();

    //! Destructor, free memory
    ~InterpointDataProxy();
public:
    //! Get the vector of the interpoint data.
    //! @param[out] vdata      vector of the interpoint data.
    vector<RDouble **> & GetVectorData();

    //! Get the vector of the interpoint data dimension.
    //! @param[out] vdim      the vector of the interpoint data dimension.
    vector<int> & GetVectorDimension();
};

#include "Geo_Interpoint.hxx"
}