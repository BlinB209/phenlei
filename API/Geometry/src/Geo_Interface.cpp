#include "Geo_Interface.h"
#include "PHMatrix.h"
#include "PHHeader.h"
#include "TK_Exit.h"
#include "GridType.h"
#include "Geo_SimpleGrid.h"
#include "Constants.h"

namespace PHSPACE
{
InterfacePatch::InterfacePatch()
{
    zoneid = 0;
    zoneIndexOfNeighbor = 0;
    faceNumberOfPatch = 0;
    faceIndexForRecv = NULL;
    faceIndexForSend = NULL;

    dataStorageSend = new Data_ParamFieldSuite();
    dataStorageRecv = new Data_ParamFieldSuite();
}

InterfacePatch::InterfacePatch(const InterfacePatch &interPatchchIN)
{
    zoneid = interPatchchIN.GetZoneID();
    zoneIndexOfNeighbor = interPatchchIN.GetZoneIndexOfNeighbor();
    faceNumberOfPatch = interPatchchIN.GetNumberOfFace();

    if (interPatchchIN.GetFaceIndexForRecv())
    {
        faceIndexForRecv = new int [faceNumberOfPatch];
        PHSPACE::CopyField(faceIndexForRecv, interPatchchIN.GetFaceIndexForRecv(), faceNumberOfPatch);
    }

    if (interPatchchIN.GetFaceIndexForSend())
    {
        faceIndexForSend = new int [faceNumberOfPatch];
        PHSPACE::CopyField(faceIndexForSend, interPatchchIN.GetFaceIndexForSend(), faceNumberOfPatch);
    }

    dataStorageSend = new Data_ParamFieldSuite();
    dataStorageRecv = new Data_ParamFieldSuite();
}

InterfacePatch & InterfacePatch::operator = (const InterfacePatch &interPatchchIN)
{
    if (this == &interPatchchIN)
    {
        return *this;
    }

    this->zoneid = interPatchchIN.GetZoneID();
    this->zoneIndexOfNeighbor = interPatchchIN.GetZoneIndexOfNeighbor();
    this->faceNumberOfPatch = interPatchchIN.GetNumberOfFace();

    if (interPatchchIN.GetFaceIndexForRecv())
    {
        if (faceIndexForRecv) delete [] faceIndexForRecv;
        faceIndexForRecv = new int [faceNumberOfPatch];
        PHSPACE::CopyField(faceIndexForRecv, interPatchchIN.GetFaceIndexForRecv(), faceNumberOfPatch);
    }

    if (interPatchchIN.GetFaceIndexForSend())
    {
        if (faceIndexForSend) delete [] faceIndexForSend;
        this->faceIndexForSend = new int [faceNumberOfPatch];
        PHSPACE::CopyField(faceIndexForSend, interPatchchIN.GetFaceIndexForSend(), faceNumberOfPatch);
    }

    delete dataStorageSend;
    delete dataStorageRecv;
    dataStorageSend = new Data_ParamFieldSuite();
    dataStorageRecv = new Data_ParamFieldSuite();
    dataStorageSend = interPatchchIN.dataStorageSend;
    dataStorageRecv = interPatchchIN.dataStorageRecv;

    return *this;
}

InterfacePatch::~InterfacePatch()
{
    delete dataStorageSend;
    delete dataStorageRecv;
    if (faceIndexForRecv) delete [] faceIndexForRecv;
    if (faceIndexForSend) delete [] faceIndexForSend;
}

void InterfacePatch::AllocateData()
{
    int faceNumberOfPatch = this->GetNumberOfFace();

    int *recvface = new int [faceNumberOfPatch];
    int *sendface = new int [faceNumberOfPatch];

    //! GMRESParallel
    int *cellindexLocal = new int[faceNumberOfPatch];
    int *cellindexNeighbor = new int[faceNumberOfPatch];
    int *cellindexLocalGhost = new int[faceNumberOfPatch];

    this->SetFaceIndexForRecv(recvface);
    this->SetFaceIndexForSend(sendface);
    
    //! GMRESParallel
    this->SetCellIndexOfLocal(cellindexLocal);
    this->SetCellIndexOfNeighbor(cellindexNeighbor);
    this->SetCellIndexOfLocalGhost(cellindexLocalGhost);
}

InterfaceInfo::InterfaceInfo()
{
    nIFace                 = 0;
    numberOfNeighbor       = 0;
    interFace2ZoneID       = NULL;
    interFace2CellID       = NULL;
    interFace2InterFaceID  = NULL;
    interFace2BoundaryFace = NULL;
    interFaceDirection     = NULL;
    //! GMRESParallel
    globalNeighborCellIndex = NULL;
    localInterfaceCellIndex = NULL;
    localInterfacePhysicalCellIndex = NULL;
    parent                 = NULL;

    dsend = new Data_ParamFieldSuite [GetNumberOfGhostCellLayers()];
    drecv = new Data_ParamFieldSuite [GetNumberOfGhostCellLayers()];

    interFacePatch.resize(0);

    isNeighborInterfaceFound = 0;
}

InterfaceInfo::InterfaceInfo(const InterfaceInfo &inforIn)
{
    nIFace = inforIn.GetNIFace();
    numberOfNeighbor = inforIn.GetNumberOfNeighbor();

    if (inforIn.GetInterFace2ZoneID())
    {
        interFace2ZoneID = new int [nIFace];
        PHSPACE::CopyField(interFace2ZoneID, inforIn.GetInterFace2ZoneID(), nIFace);
    }

    if (inforIn.GetInterFace2CellID())
    {
        interFace2CellID = new int [nIFace];
        PHSPACE::CopyField(interFace2CellID, inforIn.GetInterFace2CellID(), nIFace);
    }

    if (inforIn.GetInterFace2InterFaceID())
    {
        interFace2InterFaceID = new int [nIFace];
        PHSPACE::CopyField(interFace2InterFaceID, inforIn.GetInterFace2InterFaceID(), nIFace);
    }

    if (inforIn.GetInterFace2BoundaryFace())
    {
        interFace2BoundaryFace = new int [nIFace];
        PHSPACE::CopyField(interFace2BoundaryFace, inforIn.GetInterFace2BoundaryFace(), nIFace);
    }

    if (inforIn.GetInterFaceDirection())
    {
        interFaceDirection = new int [nIFace];
        PHSPACE::CopyField(interFaceDirection, inforIn.GetInterFaceDirection(), nIFace);
    }
    // GMRESParallel
    if(inforIn.GetGlobalNeighborCellIndex())
    {
        globalNeighborCellIndex = new int[nIFace];
        PHSPACE::CopyField(globalNeighborCellIndex, inforIn.GetGlobalNeighborCellIndex(), nIFace);
    }

    if(inforIn.GetLocalInterfaceCellIndex())
    {
        localInterfaceCellIndex = new int[nIFace];
        PHSPACE::CopyField(localInterfaceCellIndex, inforIn.GetLocalInterfaceCellIndex(), nIFace);
    }

    if(inforIn.GetLocalInterfacePhysicalCellIndex())
    {
        localInterfacePhysicalCellIndex = new int[nIFace];
        PHSPACE::CopyField(localInterfacePhysicalCellIndex, inforIn.GetLocalInterfacePhysicalCellIndex(), nIFace);
    }

    if (inforIn.GetInterfacePatch().size() != 0)
    {
        int numberOfNeighbor = this->GetNumberOfNeighbor();
        interFacePatch.resize(0);
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            InterfacePatch *interFacePatchNew = new InterfacePatch(* inforIn.GetInterfacePatch(iNeighbor));
            interFacePatch.push_back(interFacePatchNew);
        }
    }

    dsend = new Data_ParamFieldSuite [GetNumberOfGhostCellLayers()];
    drecv = new Data_ParamFieldSuite [GetNumberOfGhostCellLayers()];

    parent = 0;

    isNeighborInterfaceFound = 0;
}

InterfaceInfo::InterfaceInfo(int nIFace, Grid *parent)
{
    this->nIFace = nIFace;
    numberOfNeighbor = 0;

    interFace2ZoneID       = new int [nIFace];
    interFace2CellID       = new int [nIFace];
    interFace2InterFaceID  = new int [nIFace];
    interFace2BoundaryFace = new int [nIFace];
    interFaceDirection     = new int [nIFace];
    //! GMRESParallel
    globalNeighborCellIndex = new int[nIFace];
    localInterfaceCellIndex = new int[nIFace];
    localInterfacePhysicalCellIndex = new int[nIFace];

    this->parent = parent;

    dsend = new Data_ParamFieldSuite [GetNumberOfGhostCellLayers()];
    drecv = new Data_ParamFieldSuite [GetNumberOfGhostCellLayers()];

    interFacePatch.resize(0);

    isNeighborInterfaceFound = 0;
}

InterfaceInfo::~InterfaceInfo()
{
    if (interFace2ZoneID != NULL) delete [] interFace2ZoneID;
    if (interFace2CellID != NULL) delete [] interFace2CellID;
    if (interFace2InterFaceID != NULL) delete [] interFace2InterFaceID;
    if (interFace2BoundaryFace != NULL) delete [] interFace2BoundaryFace;
    if (interFaceDirection != NULL) delete [] interFaceDirection;
    //! GMRESParallel
    if (globalNeighborCellIndex!=NULL) delete[] globalNeighborCellIndex;
    if(localInterfaceCellIndex!=NULL) delete[] localInterfaceCellIndex;
    if(localInterfacePhysicalCellIndex!=NULL) delete localInterfacePhysicalCellIndex;

    if (dsend != NULL) delete [] dsend;
    if (drecv != NULL) delete [] drecv;

    for (vector<InterfacePatch *>::iterator iter = interFacePatch.begin(); iter != interFacePatch.end(); ++ iter)
    {
        if (NULL != *iter)
        {
            delete *iter;
            *iter = NULL;
        }
    }
    interFacePatch.clear();

    if (isNeighborInterfaceFound != NULL) delete [] isNeighborInterfaceFound;
}

void InterfaceInfo::AllocateNeighborInfo()
{
    int numberOfNeighbor = this->GetNumberOfNeighbor();
    interFacePatch.resize(0);
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
    {
        InterfacePatch *interFacePatchNew = new InterfacePatch();
        interFacePatch.push_back(interFacePatchNew);
    }
}

int InterfaceInfo::ComputeNumberOfNeighbor(int *zoneflag)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    int *interFace2ZoneID = this->GetInterFace2ZoneID();
    int nIFace = this->GetNIFace();

    //! Set zoneflag to zero, assuming that the current zone has no neighbor.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        zoneflag[iZone] = 0;
    }

    //! Find out all the neighbors of current zone,and set their zoneflag to one.
    for (int iFace = 0; iFace < nIFace; ++ iFace)
    {
        zoneflag[interFace2ZoneID[iFace]] = 1;
    }

    //! Compute the number of neighbor zones of the current zone,including itself.
    int numberOfNeighbor = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        numberOfNeighbor += zoneflag[iZone];
    }

    return numberOfNeighbor;
}

void InterfaceInfo::FindNeighbors()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    int *zoneflag = new int [nZones];

    int numberOfNeighbor = ComputeNumberOfNeighbor(zoneflag);

    this->SetNumberOfNeighbor(numberOfNeighbor);
    this->AllocateNeighborInfo();

    int ineighbor = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        if (zoneflag[iZone])
        {
            InterfacePatch *interFacePatchNew = this->GetInterfacePatch(ineighbor ++);
            interFacePatchNew->SetZoneIndexOfNeighbor(iZone);
        }
    }

    delete [] zoneflag;    zoneflag = nullptr;
}

void InterfaceInfo::ComputeFaceIndexForRecv(int ineighbor)
{
    int nIFace = this->GetNIFace();

    int *interFace2ZoneID = this->GetInterFace2ZoneID();
    int zone_id = this->GetZoneIndexOfNeighbor(ineighbor);

    InterfacePatch *interFacePatch = this->GetInterfacePatch(ineighbor);
    int *faceIndexForRecv = interFacePatch->GetFaceIndexForRecv();
    int count = 0;
    for (int iFace = 0; iFace < nIFace; ++ iFace)
    {
        if (interFace2ZoneID[iFace] == zone_id)
        {
            //! It means that faceIndexForRecv is count locally based on the interfaces of current zone.
            faceIndexForRecv[count] = iFace;
            ++ count;
        }
    }
}

int InterfaceInfo::ComputeNIFaceOfNeighbor(int ineighbor)
{
    int nIFace = this->GetNIFace();
    int *interFace2ZoneID = this->GetInterFace2ZoneID();
    int zone_id = this->GetZoneIndexOfNeighbor(ineighbor);
    int nIFaceOfNeighbor = 0;

    for (int iFace = 0; iFace < nIFace; ++ iFace)
    {
        if (interFace2ZoneID[iFace] == zone_id)
        {
            ++ nIFaceOfNeighbor;
        }
    }

    return nIFaceOfNeighbor;
}

void InterfaceInfo::AllocateSendRecv(int ineighbor)
{
    int nIFaceOfNeighbor = this->ComputeNIFaceOfNeighbor(ineighbor);

    InterfacePatch *interFacePatch = this->GetInterfacePatch(ineighbor);
    interFacePatch->SetNumberOfFace(nIFaceOfNeighbor);
    interFacePatch->AllocateData();
}

void InterfaceInfo::InitNeighbors(int ineighbor)
{
    this->AllocateSendRecv(ineighbor);
    this->ComputeFaceIndexForRecv(ineighbor);
}

void InterfaceInfo::InitNeighbors()
{
    FindNeighbors();

    int numberOfNeighbor = this->GetNumberOfNeighbor();
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
    {
        InitNeighbors(iNeighbor);
    }
}

int InterfaceInfo::FindIthNeighbor(int zone)
{
    int idx = -1;
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
    {
        if (zone == this->GetZoneIndexOfNeighbor(iNeighbor))
        {
            idx = iNeighbor;
            break;
        }
    }
    return idx;
}

int * InterfaceInfo::ComputeFaceIndexForSend(int ineighbor)
{
    int *faceIndexForSend = new int [this->GetNIFaceOfNeighbor(ineighbor)];
    int *interFace2InterFaceID = this->GetInterFace2InterFaceID();
    int zone_id = this->GetZoneIndexOfNeighbor(ineighbor);

    int count = 0;
    for (int iFace = 0; iFace < nIFace; ++ iFace)
    {
        if (interFace2ZoneID[iFace] == zone_id)
        {
            faceIndexForSend[count] = interFace2InterFaceID[iFace];
            ++ count;
        }
    }
    return faceIndexForSend;
}

void InterfaceInfo::FillFaceIndexForSend(int zone, int *faceIndexForSend_in)
{
    int ineighbor = this->FindIthNeighbor(zone);
    InterfacePatch *interFacePatch = this->GetInterfacePatch(ineighbor);
    int *faceIndexForSend = interFacePatch->GetFaceIndexForSend();
    int faceNumberOfPatch = interFacePatch->GetNumberOfFace();
    for (int iFace = 0; iFace < faceNumberOfPatch; ++ iFace)
    {
        //! Here we can make out that faceIndexForSend is corresponds to the neighbor zone.
        //! That is to say, it comes according to the corresponding sequence of corresponding neighbor zone.
        faceIndexForSend[iFace] = faceIndexForSend_in[iFace];
    }
}

void InterfaceInfo::ReSize(int numberOfNewInterfaces)
{
    if (this->nIFace == numberOfNewInterfaces) return;

    delete [] interFace2ZoneID;
    delete [] interFace2CellID;
    delete [] interFace2InterFaceID;
    delete [] interFace2BoundaryFace;
    delete [] interFaceDirection;

    this->nIFace = numberOfNewInterfaces;

    interFace2ZoneID       = new int [nIFace];
    interFace2CellID       = new int [nIFace];
    interFace2InterFaceID  = new int [nIFace];
    interFace2BoundaryFace = new int [nIFace];
    interFaceDirection     = new int [nIFace];
}

InterfaceFields::InterfaceFields(InterfaceInfo *dataIn)
{
    this->interfaceInfo = dataIn;
}

InterfaceFields::~InterfaceFields()
{
    //! The send && receive data are allocated by InterfaceFields,
    //! so, they should be free by InterfaceFields.
    if (this->Size() == 0) return;
    FreeInterfaceVar();
}

void InterfaceFields::RegisterField(const string &name, const int type, const int dimesion, const int solverID)
{
    variableNames.push_back(name);
    variableTypes.push_back(type);
    variableDimensions.push_back(dimesion);
    solverIndex.push_back(solverID);

    for (int iGhost = GetNumberOfGhostCellLayers() - 1; iGhost >= 0; -- iGhost)
    {
        AllocateInterfaceVar(name, type, dimesion, interfaceInfo->GetSendDataStorage(iGhost));
        AllocateInterfaceVar(name, type, dimesion, interfaceInfo->GetRecvDataStorage(iGhost));
    }
}

void InterfaceFields::AllocateInterfaceVar(const string &name, const int type, const int dimesion, Data_ParamFieldSuite *datastore)
{
    int nIFace = interfaceInfo->GetNIFace();

    switch (type)
    {
        case PHINT:
        {
            int **qInt = NewPointer2< int >(dimesion, nIFace);
            for (int iDim = 0; iDim < dimesion; ++ iDim)
            {
                for (int jDim = 0; jDim < nIFace; ++ jDim)
                {
                    qInt[iDim][jDim] = 0;
                }
            }
            datastore->UpdateDataPtr(name, qInt);
            break;
        }
        case PHFLOAT:
        {
            RFloat **qFloat = NewPointer2< RFloat >(dimesion, nIFace);
            for (int iDim = 0; iDim < dimesion; ++ iDim)
            {
                for (int jDim = 0; jDim < nIFace; ++ jDim)
                {
                    qFloat[iDim][jDim] = 0.0;
                }
            }
            datastore->UpdateDataPtr(name, qFloat);
            break;
        }
        case PHDOUBLE:
        {
            RDouble **qDouble = NewPointer2< RDouble >(dimesion, nIFace);
            for (int iDim = 0; iDim < dimesion; ++ iDim)
            {
                for (int jDim = 0; jDim < nIFace; ++ jDim)
                {
                    qDouble[iDim][jDim] = 0.0;
                }
            }
            datastore->UpdateDataPtr(name, qDouble);
            break;
        }
        case PHBOOL:
        {
            bool **qBool = NewPointer2< bool >(dimesion, nIFace);
            for (int iDim = 0; iDim < dimesion; ++ iDim)
            {
                for (int jDim = 0; jDim < nIFace; ++ jDim)
                {
                    qBool[iDim][jDim] = 0;
                }
            }
            datastore->UpdateDataPtr(name, qBool);
            break;
        }
        default:
        {
            TK_Exit::UnexpectedVarValue("type", type);
            break;
        }
    }
}

void InterfaceFields::RemoveAllField()
{
    for (int iVar = 0; iVar < Size(); ++ iVar)
    {
        string varName = variableNames[iVar];
        RemoveInterfaceVar(varName);
    }
    variableNames.clear();
    variableTypes.clear();
    variableDimensions.clear();
    solverIndex.clear();
}

void InterfaceFields::RemoveAnVariable(string varName)
{
    RemoveInterfaceVar(varName);

    int varIndex = static_cast<int>(find(variableNames.begin(), variableNames.end(), varName) - variableNames.begin());
    vector< string >::iterator varNameIndexBegin = variableNames.begin();
    vector< int >::iterator varTypeIndexBegin = variableTypes.begin();
    vector< int >::iterator varDimIndexBegin = variableDimensions.begin();
    vector< int >::iterator varSolverIndexBegin = solverIndex.begin();

    variableNames.erase(varNameIndexBegin + varIndex);
    variableTypes.erase(varTypeIndexBegin + varIndex);
    variableDimensions.erase(varDimIndexBegin + varIndex);
    solverIndex.erase(varSolverIndexBegin + varIndex);
}

void InterfaceFields::FreeInterfaceVar()
{
    for (int iGhost = GetNumberOfGhostCellLayers() - 1; iGhost >= 0; -- iGhost)
    {
        if (interfaceInfo->GetSendDataStorage(iGhost) != 0) FreeInterfaceVar(interfaceInfo->GetSendDataStorage(iGhost));
        if (interfaceInfo->GetRecvDataStorage(iGhost) != 0) FreeInterfaceVar(interfaceInfo->GetRecvDataStorage(iGhost));
    }
    interfaceInfo = NULL;
}

void InterfaceFields::FreeInterfaceVar(Data_ParamFieldSuite *datastore)
{
    for (int iVar = 0; iVar < Size(); ++ iVar)
    {
        string varName = variableNames[iVar];
        void **q = reinterpret_cast<void **>(datastore->GetDataPtr(varName));
        if (q != NULL)
        {
            DelPointer2(q);
            datastore->DeleteDataPtr(varName);
        }
    }
    datastore = NULL;
}

void InterfaceFields::RemoveInterfaceVar(const string &name)
{
    Data_ParamFieldSuite *datastore;
    for (int iGhost = GetNumberOfGhostCellLayers() - 1; iGhost >= 0; -- iGhost)
    {
        if (interfaceInfo->GetSendDataStorage(iGhost) != 0)
        {
            datastore = interfaceInfo->GetSendDataStorage(iGhost);
            void **q = reinterpret_cast<void **>(datastore->GetDataPtr(name));
            if (q != NULL)
            {
                DelPointer2(q);
                datastore->DeleteDataPtr(name);
            }
        }

        if (interfaceInfo->GetRecvDataStorage(iGhost) != 0)
        {
            datastore = interfaceInfo->GetRecvDataStorage(iGhost);
            void **q = reinterpret_cast<void **>(datastore->GetDataPtr(name));
            if (q != NULL)
            {
                DelPointer2(q);
                datastore->DeleteDataPtr(name);
            }
        }
    }
}

template< typename T >
T ** InterfaceFields::GetSendVariable(const string &name, const int iGhost)
{
     return reinterpret_cast<T **>(interfaceInfo->GetSendDataStorage(iGhost)->GetDataPtr(name));
}

template< typename T >
T ** InterfaceFields::GetRecvVariable(const string &name, const int iGhost)
{
    return reinterpret_cast<T **>(interfaceInfo->GetRecvDataStorage(iGhost)->GetDataPtr(name));
}

InterfaceDataProxy::InterfaceDataProxy()
{
}

InterfaceDataProxy::~InterfaceDataProxy()
{
}

vector< RDouble ** > & InterfaceDataProxy::GetVectorData()
{
    return vdata;
}

vector<int> & InterfaceDataProxy::GetVectorDimension()
{
    return vdim;
}

vector<string> & InterfaceDataProxy::GetVectorName()
{
    return dataName;
}

}
