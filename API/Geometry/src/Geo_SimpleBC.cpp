#include "GridType.h"
#include "Geo_Interface.h"
#include "Geo_FaceTopo_Unstruct.h"
#include "Geo_StructGrid.h"
#include "Geo_MultiGridInfo_Struct.h"
#include "Geo_UnstructGrid.h"
#include "Geo_FaceMetrics_Unstruct.h"
#include "Geo_CellMetrics_Unstruct.h"
#include "Geo_DynamicGridMetrics_Unstruct.h"
#include "Geo_LSQWeight_Unstruct.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "Constants.h"
#include "TK_Exit.h"
#include "TK_Parse.h"
#include "Glb_Dimension.h"
#include "Geo_StructBC.h"
#include "Geo_Interpoint.h"
#include "Geo_NodeTopo_Struct.h"
#include "Geo_FaceMetrics_Struct.h"
#include "Geo_CellMetrics_Struct.h"
#include "Geo_DynamicGridMetrics_Struct.h"
#include "Geo_OversetGridTopo_Struct.h"
#include "Geo_UnstructBC.h"
#include "Constants.h"
#include "Gas.h"
using namespace std;

namespace PHSPACE
{
vector<SimpleBC *> * GlobalBoundaryCondition::globalBCList = 0;
map < int, int > *GlobalBoundaryCondition::bodyMap = 0;
set < string > * GlobalBoundaryCondition::bodyList = 0;

Data_Param * SimpleBC::GetBCParamDataBase()
{
    return bcParamDataBase;
}

Data_Field * SimpleBC::GetBCFieldDataBase()
{
    return bcFieldDataBase;
}

GlobalBoundaryCondition::GlobalBoundaryCondition() 
{

}

GlobalBoundaryCondition::~GlobalBoundaryCondition()
{
    //! No used for static class.
     vector<SimpleBC *>::iterator iter;
     for (iter = globalBCList->begin(); iter != globalBCList->end(); ++ iter)
     {
         delete (*iter);
     }
     globalBCList->clear();
}

void GlobalBoundaryCondition::ReadGlobalBoundaryCondition()
{
    fstream file;
    string bcFile = "./bin/boundary_condition.hypara";

    if (!FileExist(bcFile))
    {
        //! To compatible the old version.
        SetDefaultSolidBoundaryCondition();
        return;
    }

    OpenFile(file, bcFile, ios::in);

    ParseBCFromFile(file);

    BuildBodyMap();

    CloseFile(file);

    ChangeBCTypeByGlobalBC();
}

void GlobalBoundaryCondition::SetDefaultSolidBoundaryCondition()
{
    if (globalBCList != 0) return;

    globalBCList = new vector<SimpleBC *>;
    set < pair<string, int> > bcNameSet;

    map<int, string> bcnameMap;
    bcnameMap.insert(pair<int, string>(0 , "NO_BOUNDARY_CONDITION"));
    bcnameMap.insert(pair<int, string>(1 , "EXTRAPOLATION"));
    bcnameMap.insert(pair<int, string>(2 , "SOLID_SURFACE"));
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

    using namespace PHMPI;

    int nZones = PHMPI::GetNumberofGlobalZones();
    int myid = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int procID = PHMPI::GetZoneProcessorID(iZone);
        if (myid == procID)
        {
            Grid *gridofThisProcessor = GetGrid(iZone, 0);

            int gridtype = gridofThisProcessor->Type();
            if (gridtype == PHSPACE::STRUCTGRID)
            {
                StructGrid *grid = StructGridCast(gridofThisProcessor);

                StructBCSet *structBCSet = grid->GetStructBCSet();
                int nBCRegion = structBCSet->GetnBCRegion();
               
                for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
                {
                    
                    StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
                    string bcName = bcregion->GetBCName();
                    int bcType = bcregion->GetBCType();
                    if (bcType == PHENGLEI::INTERFACE)
                    {
                        continue;
                    }

                    if (bcName.empty())
                    {
                        bcName = bcnameMap[bcType];
                        bcregion->SetBCName(bcName);
                    }

                    pair<string,int> bcNamePair(bcName ,bcType);
                    bcNameSet.insert(bcNamePair);
                }
            }
            else
            {
                UnstructGrid *grid = UnstructGridCast(gridofThisProcessor);

                UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
                if (!unstructBCSet) continue;
                int nBCRegion = unstructBCSet->GetnBCRegion();
           
                for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
                {
                    UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
                    string bcName = bcRegion->GetBCName();
                    int bcType = bcRegion->GetBCType();
                    if (bcType == PHENGLEI::INTERFACE)
                    {
                        continue;
                    }

                    if (bcName.empty())
                    {
                        bcName = bcnameMap[bcType];
                        bcRegion->SetBCName(bcName);
                    }

                    pair<string,int> bcNamePair(bcName ,bcType);
                    bcNameSet.insert(bcNamePair);
                }

            }
        }
    }

    vector<DataContainer *> dataList;

    int nProc  = PHMPI::GetNumberOfProcessor();
    int serverTmp = PHMPI::GetServerProcessorID();
    for (int iProc = 0; iProc < nProc; ++ iProc)
    {
        DataContainer *eachProc2Server = new DataContainer();

        if (myid == iProc)
        {
            eachProc2Server->MoveToBegin();

            int bcNameSetSize = static_cast<int>(bcNameSet.size());
            PHWrite(eachProc2Server, bcNameSetSize);

            set < pair<string, int> >::iterator bcNameSetIter;
            for (bcNameSetIter = bcNameSet.begin(); bcNameSetIter != bcNameSet.end(); ++ bcNameSetIter)
            {
                string bcName = (*bcNameSetIter).first;
                int bcType = (*bcNameSetIter).second;

                eachProc2Server->WriteString(bcName);
                PHWrite(eachProc2Server, bcType);
            }
        }

        PH_Trade(eachProc2Server, iProc, serverTmp, iProc);

        if (myid == serverTmp)
        {
            dataList.push_back(eachProc2Server);
            eachProc2Server = new DataContainer();
        }
        delete eachProc2Server;    eachProc2Server = nullptr;
    }

    bcNameSet.clear();
    for (int iData = 0; iData < dataList.size(); ++ iData)
    {
        DataContainer *iBcNameData = dataList[iData];
        iBcNameData->MoveToBegin();

        int bcNameSetSize = 0;
        PHRead(iBcNameData, bcNameSetSize);

        for (int iBC = 0; iBC < bcNameSetSize; ++ iBC)
        {
            string bcName;
            int bcType;

            iBcNameData->ReadString(bcName);
            PHRead(iBcNameData, bcType);

            pair<string,int> bcNamePair(bcName ,bcType);
            bcNameSet.insert(bcNamePair);
        }
    }

    for (int iData = 0; iData < dataList.size(); ++ iData)
    {
        delete dataList[iData];
    }

    DataContainer *server2EachProc = new DataContainer();
    if (myid == serverTmp)
    {
        server2EachProc->MoveToBegin();

        int bcNameSetSize = static_cast<int>(bcNameSet.size());
        PHWrite(server2EachProc, bcNameSetSize);

        set < pair<string, int> >::iterator bcNameSetIter;
        for (bcNameSetIter = bcNameSet.begin(); bcNameSetIter != bcNameSet.end(); ++ bcNameSetIter)
        {
            string bcName = (*bcNameSetIter).first;
            int bcType = (*bcNameSetIter).second;

            server2EachProc->WriteString(bcName);
            PHWrite(server2EachProc, bcType);
        }

        for (int iProc = 0; iProc < nProc; ++ iProc)
        {
            if (iProc == myid)
            {
                continue;
            }
            send(server2EachProc, iProc, iProc);
        }
    }
    else
    {
        receive(server2EachProc, serverTmp, myid);
        server2EachProc->MoveToBegin();

        int bcNameSetSize = 0;
        PHRead(server2EachProc, bcNameSetSize);

        for (int iBC= 0; iBC < bcNameSetSize; ++ iBC)
        {
            string bcName;
            int bcType;

            server2EachProc->ReadString(bcName);
            PHRead(server2EachProc, bcType);

            pair<string,int> bcNamePair(bcName ,bcType);
            bcNameSet.insert(bcNamePair);
        }
    }
    delete server2EachProc;    server2EachProc = nullptr;

    set < pair<string, int> >::iterator bcNameSetIter;
    for (bcNameSetIter = bcNameSet.begin(); bcNameSetIter != bcNameSet.end(); ++ bcNameSetIter)
    {
        string bcName = (*bcNameSetIter).first;
        int bcType = (*bcNameSetIter).second;

        SimpleBC *boundaryCondition = new SimpleBC();
        boundaryCondition->SetBCName(bcName);
        boundaryCondition->SetBCType(bcType);

        Data_Param *bcParamDataBase = new Data_Param();
        bcParamDataBase->UpdateData("bcType", &bcType, PHINT, 1);

        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            string bodyName = "body";
            if (!bcParamDataBase->CheckDataExist("bodyName"))
            {
                bcParamDataBase->UpdateData("bodyName", &bodyName, PHSTRING, 1);
            }

            int taskType = GetTaskCode();
            if (taskType == SOLVE_FIELD || taskType == POST_PROCESSING)
            {
                //! For solid wall BC, set the default reference geometry scale.
                RDouble partForceReferenceArea, partForceReferenceLength, partForceReferenceLength_B;
                RDouble partForceReferenceX, partForceReferenceY, partForceReferenceZ;

                //! This judgement is not enough, and needs to be improved. In fact, judgements of the five reference variables should all exist.

                partForceReferenceArea = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
                partForceReferenceLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
                partForceReferenceLength_B = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");
                partForceReferenceX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
                partForceReferenceY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
                partForceReferenceZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

                bcParamDataBase->UpdateData("forceReferenceArea", &partForceReferenceArea, PHDOUBLE, 1);
                bcParamDataBase->UpdateData("forceReferenceLength", &partForceReferenceLength, PHDOUBLE, 1);
                bcParamDataBase->UpdateData("forceReferenceLengthSpanWise", &partForceReferenceLength_B, PHDOUBLE, 1);
                bcParamDataBase->UpdateData("TorqueRefX", &partForceReferenceX, PHDOUBLE, 1);
                bcParamDataBase->UpdateData("TorqueRefY", &partForceReferenceY, PHDOUBLE, 1);
                bcParamDataBase->UpdateData("TorqueRefZ", &partForceReferenceZ, PHDOUBLE, 1);

                int dumpHingeMoment = 0;
                RDouble localCoordAxis0[3] = { 0, 0, 0 };
                RDouble localCoordAxis1[3] = { 0, 0, 0 };

                bcParamDataBase->UpdateData("dumpHingeMoment", &dumpHingeMoment, PHINT, 1);
                bcParamDataBase->UpdateData("localCoordAxis0", localCoordAxis0, PHDOUBLE, 3);
                bcParamDataBase->UpdateData("localCoordAxis1", localCoordAxis1, PHDOUBLE, 3);
            }
        }

        boundaryCondition->SetBCParamDataBase(bcParamDataBase);
        globalBCList->push_back(boundaryCondition);
    }
}

void GlobalBoundaryCondition::BuildBodyMap()
{
    vector<SimpleBC *> *globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();
    if (globalBCList == 0)
    {
        return;
    }

    set < pair<int, string> > bodyNameSet;
    if (bodyList)
    {
        delete bodyList;    bodyList = nullptr;
    }
    bodyList = new set< string >;

    if (bodyMap)
    {
        delete bodyMap;    bodyMap = nullptr;
    }
    bodyMap = new map < int, int >;

    for (int iBC = 0; iBC < globalBCList->size(); ++ iBC)
    {
        SimpleBC *boundaryCondition = (*globalBCList)[iBC];
        int bcType = boundaryCondition->GetBCType();

        /*if (bcType != PHENGLEI::SOLID_SURFACE && bcType != PHENGLEI::ABLATION_SURFACE)
        {
            continue;
        }*/

        bool isSetBodyName = boundaryCondition->CheckParamData("bodyName");
        string bodyName = "body";
        if (isSetBodyName)
        {
            boundaryCondition->GetBCParamDataBase()->GetData("bodyName", &bodyName, PHSTRING, 1);
            pair<int, string> bodyNamePair(iBC, bodyName);
            bodyNameSet.insert(bodyNamePair);
            bodyList->insert(bodyName);
        }
    }

    set < pair<int, string> >::iterator bodyNameSetIter;
    for (bodyNameSetIter = bodyNameSet.begin(); bodyNameSetIter != bodyNameSet.end(); ++ bodyNameSetIter)
    {
        int iBC = (*bodyNameSetIter).first;
        string bodyName = (*bodyNameSetIter).second;

        set < string >::iterator bodyListIter;
        int bodyIndex = 0;
        for (bodyListIter = bodyList->begin(); bodyListIter != bodyList->end(); ++ bodyListIter)
        {
            string bodyNameInList = *bodyListIter;
            if (bodyNameInList == bodyName)
            {
                bodyMap->insert(pair<int,int>(iBC, bodyIndex));
                break;
            }
            bodyIndex ++;
        }
    }
}

void GlobalBoundaryCondition::InitMassFlowBoundary()
{
    vector<SimpleBC *> *globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();
    if (globalBCList == 0)
    {
        return;
    }

    for (int iBC = 0; iBC < globalBCList->size(); ++ iBC)
    {
        SimpleBC *boundaryCondition = (*globalBCList)[iBC];
        int bcType = boundaryCondition->GetBCType();

        if (bcType == PHENGLEI::MASS_FLOW_INLET || bcType == PHENGLEI::MASS_FLOW_OUTLET)
        {
            RDouble BCTotalArea = 0.0;
            string bcName = boundaryCondition->GetBCName();

            int nLocalZones = PHMPI::GetNumberofLocalZones();
            for (int iZone = 0; iZone < nLocalZones; iZone++)
            {
                int iZoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);
                Grid *iGrid = PHSPACE::GetGrid(iZoneID);

                int gridType = iGrid->Type();
                if (gridType == UNSTRUCTGRID)
                {
                    UnstructGrid *grid = UnstructGridCast(iGrid);

                    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
                    int nBCRegion = unstructBCSet->GetnBCRegion();
                    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
                    {
                        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
                        string jBCName = bcRegion->GetBCName();

                        if (jBCName == bcName)
                        {
                            int iFace;
                            vector<int> *faceIndex = bcRegion->GetFaceIndex();
                            RDouble     *area      = grid->GetFaceArea();
                            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                            {
                                iFace = *iter;
                                BCTotalArea += area[iFace];
                            }
                        }
                    }
                }
                else
                {
                    StructGrid *grid = StructGridCast(iGrid);

                    StructBCSet *structBCSet = grid->GetStructBCSet();
                    int nBCRegion = structBCSet->GetnBCRegion();
                    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
                    {
                        StructBC *jBCRegion = structBCSet->GetBCRegion(iBCRegion);
                        string jBCName = jBCRegion->GetBCName();

                        //! Judge the same BC region according the BC name.
                        if (jBCName == bcName)
                        {
                            int ist, ied, jst, jed, kst, ked;
                            int in, jn, kn;
                            jBCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);
                            int iSurface = jBCRegion->GetFaceDirection() + 1;
                            RDouble4D &area = *(grid->GetFaceArea());

                            for (int k = kst; k <= ked; ++ k)
                            {
                                for (int j = jst; j <= jed; ++ j)
                                {
                                    for (int i = ist; i <= ied; ++ i)
                                    {
                                        jBCRegion->GetBoundaryFaceIndex(i, j, k, in, jn, kn);
                                        BCTotalArea += area(in, jn, kn, iSurface);
                                    }
                                }
                            }
                        }
                    }
                }
            }


            RDouble GLB_BCTotalArea;
            PH_AllReduce(&BCTotalArea, &GLB_BCTotalArea, 1, MPI_SUM);
            BCTotalArea = GLB_BCTotalArea;

            boundaryCondition->UpdateParamData("BCTotalArea", &BCTotalArea, PHDOUBLE, 1);
        }
    }
}

void GlobalBoundaryCondition::SetGlobalBCByGlobalDataBase()
{
    vector<SimpleBC *> *globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();
    if (globalBCList == 0) return;

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int ntmodel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    int nEquation = 5;
    if (nchem > 0)
    {
        nEquation = nl + nchem + ntmodel - 1;
    }

    for (int iBC = 0; iBC < globalBCList->size(); ++ iBC)
    {
        SimpleBC *boundaryCondition = (*globalBCList)[iBC];
        int bcType = boundaryCondition->GetBCType();
        Data_Param *bcParamDB = boundaryCondition->GetBCParamDataBase();

        if (IsInterface(bcType))
        {
            continue;
        }
        else if (bcType == PHENGLEI::INFLOW || bcType == PHENGLEI::FARFIELD)
        {
            RDouble *prim_inf = new RDouble [nEquation];
            if (bcParamDB->IsExist("inflowParaType", PHINT, 1))
            {
                int inflowParaType;
                bcParamDB->GetData("inflowParaType", &inflowParaType, PHINT, 1);

                if (inflowParaType == NONDIMENSIONCONDITION)
                {
                    SetBCInitialValuesByReynolds(bcParamDB, prim_inf);
                }
                else if (inflowParaType == FLIGHTCONDITION)
                {
                    SetBCInitialValuesByHeight(bcParamDB, prim_inf);
                    boundaryCondition->UpdateParamData("primitiveVarFarfield", prim_inf, PHDOUBLE, nEquation);
                }
                else if (inflowParaType == EXPERIMENTCONDITION)
                {
                    SetBCInitialValuesByTotalPressure(bcParamDB, prim_inf);
                }
                else if (inflowParaType == TEMPERATURE_DENSITY)
                {
                    SetBCInitialValuesByDensity(bcParamDB, prim_inf);
                }
                else if (inflowParaType == TEMPERATURE_PRESSURE)
                {
                    SetBCInitialValuesByPressure(bcParamDB, prim_inf);
                }
                else if (inflowParaType == MACH_TEMP_PRE)
                {
                    SetBCInitialValuesByMachNumberTemperaturePressure(bcParamDB, prim_inf);
                }
                else if (inflowParaType == WEATHERCONDITION)
                {
                    SetBCInitialValuesByWRFData(bcParamDB, prim_inf);
                }
                else if (inflowParaType == PRIMITIVE_VARIABLES)
                {
                    SetBCInitialValuesByPrimitive(bcParamDB, prim_inf);
                }
                else if (inflowParaType == TRAJECTORY)
                {
                    RDouble *primitiveVarFarfield = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("prim_inf"));
                    for (int s = 0; s < nEquation; ++ s)
                    {
                        prim_inf[s] = primitiveVarFarfield[s];
                    }
                    SetSpeedDrictionForBC(bcParamDB, 1.0, prim_inf);
                }
                else
                {
                    TK_Exit::UnexpectedVarValue("inflowParaType", inflowParaType);
                }
                bcParamDB->UpdateData("primitiveVarFarfield", prim_inf, PHDOUBLE, nEquation);
            }
            else
            {
                RDouble *primitiveVarFarfield = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("prim_inf"));
                for (int s = 0; s < nEquation; ++ s)
                {
                    prim_inf[s] = primitiveVarFarfield[s];
                }
                SetSpeedDrictionForBC(bcParamDB, 1.0, prim_inf);
                boundaryCondition->UpdateParamData("primitiveVarFarfield", prim_inf, PHDOUBLE, nEquation);

            }
            delete [] prim_inf;    prim_inf = nullptr;
        }
        else if (bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::ABLATION_SURFACE)
        {
            int viscousType;
            RDouble wallTemperature;
            RDouble catalyticCoef;

            if (!bcParamDB->IsExist("catalyticCoef", PHDOUBLE, 1))
            {
                catalyticCoef = GlobalDataBase::GetDoubleParaFromDB("catalyticCoef");
                bcParamDB->UpdateData("catalyticCoef", &catalyticCoef, PHDOUBLE, 1);
            }

            if (!bcParamDB->IsExist("viscousType", PHINT, 1))
            {
                GlobalDataBase::GetData("viscousType", &viscousType, PHINT, 1);
                boundaryCondition->UpdateParamData("viscousType", &viscousType, PHINT, 1);
            }

            if (!bcParamDB->IsExist("wallTemperature", PHDOUBLE, 1))
            {
                GlobalDataBase::GetData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
                boundaryCondition->UpdateParamData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
            }
            /*else
            {
                boundaryCondition->GetParamData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
                GlobalDataBase::UpdateData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
            }*/
        }
        else if (bcType == PHENGLEI::PRESSURE_INLET)
        {
            //RDouble *primitiveVarFarfield = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("prim_inf"));
            //bcParamDB->UpdateData("primitiveVarFarfield", primitiveVarFarfield, PHDOUBLE, nEquation);
            //SetMaximumSpeciesForBC(bcParamDB);

             RDouble *prim_inf = new RDouble [nEquation];
             SetBCInitialValuesByTotalPressure(bcParamDB, prim_inf);
             boundaryCondition->UpdateParamData("primitiveVarFarfield", prim_inf, PHDOUBLE, nEquation);
             delete [] prim_inf;    prim_inf = nullptr;
        }
    }

    //for (iter = globalBCList->begin(); iter != globalBCList->end(); ++ iter)
    //{
    //    SimpleBC *bc = *iter;
    //    bcnameMap.insert(pair<string, int>(bc->GetBCName(), bc->GetBCType()));
    //
    //    Data_Param *bcParamDB = bc->GetBCParamDataBase();
    //
    //    if (bcParamDB->IsExist("refMachNumber", PHDOUBLE, 1))
    //    {
    //        RDouble refMachNumber;
    //        bcParamDB->GetData("refMachNumber", &refMachNumber, PHDOUBLE, 1);
    //        GlobalDataBase::UpdateData("refMachNumber", &refMachNumber, PHDOUBLE, 1);
    //    }
    //
    //    if (bcParamDB->IsExist("attackd", PHDOUBLE, 1))
    //    {
    //        RDouble attackd;
    //        bcParamDB->GetData("attackd", &attackd, PHDOUBLE, 1);
    //        GlobalDataBase::UpdateData("attackd", &attackd, PHDOUBLE, 1);
    //    }
    //
    //    if (bcParamDB->IsExist("angleSlide", PHDOUBLE, 1))
    //    {
    //        RDouble angleSlide;
    //        bcParamDB->GetData("angleSlide", &angleSlide, PHDOUBLE, 1);
    //        GlobalDataBase::UpdateData("angleSlide", &angleSlide, PHDOUBLE, 1);
    //    }
    //
    //    if (bcParamDB->IsExist("refReNumber", PHDOUBLE, 1))
    //    {
    //        RDouble refReNumber;
    //        bcParamDB->GetData("refReNumber", &refReNumber, PHDOUBLE, 1);
    //        GlobalDataBase::UpdateData("refReNumber", &refReNumber, PHDOUBLE, 1);
    //    }
    //
    //    if (bcParamDB->IsExist("refDimensionalTemperature", PHDOUBLE, 1))
    //    {
    //        RDouble refDimensionalTemperature;
    //        bcParamDB->GetData("refDimensionalTemperature", &refDimensionalTemperature, PHDOUBLE, 1);
    //        GlobalDataBase::UpdateData("refDimensionalTemperature", &refDimensionalTemperature, PHDOUBLE, 1);
    //    }
    //
    //    if (bcParamDB->IsExist("wallTemperature", PHDOUBLE, 1))
    //    {
    //        RDouble wallTemperature;
    //        bcParamDB->GetData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
    //        GlobalDataBase::UpdateData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
    //    }
    //
    //    if (bcParamDB->IsExist("height", PHDOUBLE, 1))
    //    {
    //        RDouble height;
    //        bcParamDB->GetData("height", &height, PHDOUBLE, 1);
    //        GlobalDataBase::UpdateData("height", &height, PHDOUBLE, 1);
    //    }
    //
    //    if (bcParamDB->IsExist("inflowParaType", PHINT, 1))
    //    {
    //        RDouble inflowParaType;
    //        bcParamDB->GetData("inflowParaType", &inflowParaType, PHINT, 1);
    //        GlobalDataBase::UpdateData("inflowParaType", &inflowParaType, PHINT, 1);
    //    }
    //}
}

void GlobalBoundaryCondition::SetBCInitialValuesByReynolds(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    RDouble localMachNumber = 0, localReNumber = 0, localDimensionalTemperature = 0, localVibrationTemperature = 0;
    //! Obtain the local parameters on boundary.
    if (bcParamDB->IsExist("refMachNumber", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refMachNumber", &localMachNumber, PHDOUBLE, 1);
    }
    else
    {
        localMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    }
    if (bcParamDB->IsExist("refReNumber", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refReNumber", &localReNumber, PHDOUBLE, 1);
    }
    else
    {
        localReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    }
    if (bcParamDB->IsExist("refDimensionalTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalTemperature", &localDimensionalTemperature, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    }

    int nTvChange = GlobalDataBase::GetIntParaFromDB("nTvChange");
    if(nTvChange == 0)
    {
        if (bcParamDB->IsExist("refDimensionalVibrationTemperature", PHDOUBLE, 1))
        {
            bcParamDB->GetData("refDimensionalVibrationTemperature", &localVibrationTemperature, PHDOUBLE, 1);
        }
        else
        {
            localVibrationTemperature = GlobalDataBase::GetDoubleParaFromDB("freestream_vibration_temperature");
        }
    }

    //! Obtain the reference parameters.
    RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    RDouble localDimensionalDensity, localDimensionalVelocity, localDimensionalPressure;
    RDouble nonDimTtr = localDimensionalTemperature / refDimensionalTemperature;
    RDouble nonDimTv  = localVibrationTemperature / refDimensionalTemperature;

    int nchem = gas->GetChemicalType();
    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    if (bcParamDB->IsExist("refGama", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refGama", &refGama, PHDOUBLE, 1);
    }

    if (nchem == 0)
    {
        RDouble t0i  = 288.15;
        RDouble tsi  = 110.4;
        RDouble mu0i = 1.7894E-05;

        RDouble local_viscosity_dimensional = (t0i + tsi) / (localDimensionalTemperature + tsi) * pow(localDimensionalTemperature / t0i, 1.5) * mu0i;        
        RDouble reference_general_gas_constant = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");

        RDouble localDimensionalSonicSpeed = sqrt(refGama * reference_general_gas_constant * localDimensionalTemperature);
        localDimensionalVelocity = localDimensionalSonicSpeed * localMachNumber;

        localDimensionalDensity = localReNumber * local_viscosity_dimensional / localDimensionalVelocity;
        localDimensionalPressure = localDimensionalDensity * reference_general_gas_constant * localDimensionalTemperature;
    }
    else
    {
        int nNSEqn = gas->GetNSEquationNumber();
        int numberOfSpecies = gas->GetNumberOfSpecies();
        RDouble initMassFraction[MAX_SPECIES_NUM] = {0}, speciesCv[MAX_SPECIES_NUM] = {0}, primVars[MAX_SPECIES_NUM] = {0}, speciesViscosity[MAX_SPECIES_NUM] = {0};
        RDouble speciesWeight[MAX_SPECIES_NUM] = {0}, moleFractions[MAX_SPECIES_NUM] = {0};
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->MassFractionConversion(initMassFraction);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *mass_fraction_reference = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = mass_fraction_reference[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }

        RDouble oMass = 0.0, gasR = 0.0;
        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
        gas->ComputeMoleFractionByMassFraction(initMassFraction, moleFractions);
        gasR = rjmk * oMass;

        RDouble cp = 0.0, cv = 0.0;
        int ntmodel = gas->GetTemperatureModel();
        if (ntmodel > 1)    //! multi-temperature model.
        {
            gas->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cv);
            cp = cv + gasR;
            refGama = cp / cv;
            gas->ComputeSpeciesDimensionalViscosityWithCurveFitMethod(nonDimTtr, nonDimTv, speciesViscosity);

            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            electronEnergy  = electronEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
        else    //! one-temperature model.
        {
            gas->ComputeSpeciesConstantPressureSpecificHeatDimensional(nonDimTtr, speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cp);
            cv = cp - gasR;
            refGama = cp / cv;
            gas->ComputeOneTemperatureModelSpeciesViscosityDimensional(nonDimTtr, speciesViscosity);
        }

        RDouble localDimensionalSonicSpeed = sqrt(refGama * gasR * localDimensionalTemperature);
        localDimensionalVelocity = localDimensionalSonicSpeed * localMachNumber;

        RDouble *MolecularWeight = gas->GetMolecularWeightDimensional();
        gas->GetPartitionFunctionValues(speciesViscosity, MolecularWeight, moleFractions, speciesWeight);
        RDouble local_viscosity_dimensional = gas->GetMixedGasViscosityWithWilkeFormula(moleFractions, speciesViscosity, speciesWeight);

        localDimensionalDensity = localReNumber * local_viscosity_dimensional / localDimensionalVelocity;
        localDimensionalPressure = localDimensionalDensity * gasR * localDimensionalTemperature;

        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }
    }

    using namespace IDX;

    primitiveVar[IR] = localDimensionalDensity / refDimensionalDensity;

    SetSpeedDrictionForBC(bcParamDB,localDimensionalVelocity / refDimensionalVelocity, primitiveVar);

    primitiveVar[IP] = localDimensionalPressure / (refDimensionalDensity * refDimensionalVelocity * refDimensionalVelocity);

    bcParamDB->UpdateData("refGama", &refGama, PHDOUBLE, 1);

}

void GlobalBoundaryCondition::SetBCInitialValuesByHeight(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    RDouble localMachNumber = 0, height = 0, localDimensionalTemperature = 0, localVibrationTemperature = 0;
    //! Obtain the local parameters on boundary.
    if (bcParamDB->IsExist("refMachNumber", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refMachNumber", &localMachNumber, PHDOUBLE, 1);
    }
    else
    {
        localMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    }
    if (bcParamDB->IsExist("height", PHDOUBLE, 1))
    {
        bcParamDB->GetData("height", &height, PHDOUBLE, 1);
    }
    else
    {
        height = GlobalDataBase::GetDoubleParaFromDB("height");
    }
    RDouble localDimensionalDensity, localDimensionalVelocity, localDimensionalPressure = 1.0, localDimensionalSonicSpeed;
    int nGasModel = GlobalDataBase::GetIntParaFromDB("nGasModel");
   
    gas->GetAirInfo(height, localDimensionalTemperature, localDimensionalPressure, localDimensionalDensity, localDimensionalSonicSpeed);
   
    int nTvChange = GlobalDataBase::GetIntParaFromDB("nTvChange");
    if(nTvChange == 0)
    {
        if (bcParamDB->IsExist("refDimensionalVibrationTemperature", PHDOUBLE, 1))
        {
            bcParamDB->GetData("refDimensionalVibrationTemperature", &localVibrationTemperature, PHDOUBLE, 1);
        }
        else
        {
            localVibrationTemperature = localDimensionalTemperature;
        }
    }

    //! Obtain the reference parameters.
    RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    RDouble nonDimTtr = localDimensionalTemperature / refDimensionalTemperature;
    RDouble nonDimTv  = localVibrationTemperature / refDimensionalTemperature;

    int nchem = gas->GetChemicalType();
    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    if (bcParamDB->IsExist("refGama", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refGama", &refGama, PHDOUBLE, 1);
    }

    if (nchem == 0)
    {
        RDouble reference_general_gas_constant = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");

        //! unit: m/s.
        localDimensionalSonicSpeed = sqrt(refGama * reference_general_gas_constant * localDimensionalTemperature);
        localDimensionalVelocity = localDimensionalSonicSpeed * localMachNumber;
        localDimensionalDensity = localDimensionalPressure / reference_general_gas_constant / localDimensionalTemperature;
    }
    else
    {
        int nNSEqn = gas->GetNSEquationNumber();
        int numberOfSpecies = gas->GetNumberOfSpecies();
        RDouble initMassFraction[MAX_SPECIES_NUM] = {0}, speciesCv[MAX_SPECIES_NUM] = {0}, primVars[MAX_SPECIES_NUM] = {0};
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->MassFractionConversion(initMassFraction);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *mass_fraction_reference = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = mass_fraction_reference[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }

        RDouble oMass = 0.0, gasR = 0.0;
        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
        gasR = rjmk * oMass;

        RDouble cp = 0.0, cv = 0.0;
        int ntmodel = gas->GetTemperatureModel();
        if (ntmodel > 1)    //! multi-temperature model.
        {
            gas->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cv);
            cp = cv + gasR;
            refGama = cp / cv;

            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            electronEnergy  = electronEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
        else    //! one-temperature model.
        {
            gas->ComputeSpeciesConstantPressureSpecificHeatDimensional(nonDimTtr, speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cp);
            cv = cp - gasR;
            refGama = cp / cv;
        }

        localDimensionalSonicSpeed = sqrt(refGama * gasR * localDimensionalTemperature);
        localDimensionalVelocity = localDimensionalSonicSpeed * localMachNumber;
        localDimensionalDensity = localDimensionalPressure / (gasR * localDimensionalTemperature);

        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }
    }

    using namespace IDX;
    primitiveVar[IR] = localDimensionalDensity / refDimensionalDensity;

    SetSpeedDrictionForBC(bcParamDB,localDimensionalVelocity / refDimensionalVelocity, primitiveVar);

    primitiveVar[IP] = localDimensionalPressure / (refDimensionalDensity * refDimensionalVelocity * refDimensionalVelocity);

    bcParamDB->UpdateData("refGama", &refGama, PHDOUBLE, 1);
}

void GlobalBoundaryCondition::SetMaximumSpeciesForBC(Data_Param *bcParamDB)
{
    using namespace GAS_SPACE;
    int nchem = gas->GetChemicalType();
    if (nchem > 0)
    {
        int numberOfSpecies = gas->GetNumberOfSpecies();
        RDouble initMassFraction[MAX_SPECIES_NUM] = {0};
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->SetMaximumSpecies(initMassFraction);
        }
    }
}

void GlobalBoundaryCondition::SetSpeedDrictionForBC(Data_Param *bcParamDB, RDouble localNonDimensionalVelocity, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    using namespace IDX;
    using namespace INIFLOW_SPACE;
  /*const int CUSTOM_ANGLE           = 0;
    const int CUSTOM_VECTOR          = 1;
  */
    int directionMethod = CUSTOM_ANGLE;
    if (bcParamDB->IsExist("directionMethod", PHINT, 1))
    {
        bcParamDB->GetData("directionMethod", &directionMethod, PHINT, 1);
    }

    RDouble direction_inlet[3] = {1.0,0.0,0.0};

    if(directionMethod == CUSTOM_VECTOR )
    {
        if (bcParamDB->IsExist("direction_inlet", PHDOUBLE, 3))
        {
            bcParamDB->GetData("direction_inlet", &direction_inlet, PHDOUBLE, 3);

            RDouble sum = sqrt( direction_inlet[0] * direction_inlet[0] + direction_inlet[1] * direction_inlet[1] + direction_inlet[2] * direction_inlet[2]);
            if(sum < SMALL)
            {
                direction_inlet[0] = 1.0;
                direction_inlet[1] = 0.0;
                direction_inlet[2] = 0.0;
            }
            else
            {
                sum = 1.0 / sum;
                direction_inlet[0] *= sum;
                direction_inlet[1] *= sum;
                direction_inlet[2] *= sum;
            }
        }
        else
        {
            if(GlobalDataBase::GetIntParaFromDB("directionMethod") == CUSTOM_VECTOR)
            {
                GlobalDataBase::GetData("direction_inlet", &direction_inlet, PHDOUBLE, 3);
            }
            else
            {
                RDouble attack   =  PI / 180.0 * GlobalDataBase::GetDoubleParaFromDB("attackd");
                RDouble sideslip =  PI / 180.0 * GlobalDataBase::GetDoubleParaFromDB("angleSlide");

                direction_inlet[0] = cos(attack) * cos(sideslip);
                direction_inlet[1] = sin(attack) * cos(sideslip);
                direction_inlet[2] = sin(sideslip);
            }
        }
    }
    else    //! Attackd and angleSlide
    {
        RDouble attackd = 0.0, angleSlide = 0.0;
        if (bcParamDB->IsExist("attackd", PHDOUBLE, 1))
        {
            bcParamDB->GetData("attackd", &attackd, PHDOUBLE, 1);
            if (bcParamDB->IsExist("angleSlide", PHDOUBLE, 1))
            {
                bcParamDB->GetData("angleSlide", &angleSlide, PHDOUBLE, 1);
            }

            attackd    *= PI / 180.0;
            angleSlide *= PI / 180.0;

            direction_inlet[0] = cos(attackd) * cos(angleSlide);
            direction_inlet[1] = sin(attackd) * cos(angleSlide);
            direction_inlet[2] = sin(angleSlide);
        }
        else
        {
            if(GlobalDataBase::GetIntParaFromDB("directionMethod") == CUSTOM_VECTOR)
            {
                GlobalDataBase::GetData("direction_inlet", &direction_inlet, PHDOUBLE, 3);
            }
            else
            {
                RDouble attack   =  PI / 180.0 * GlobalDataBase::GetDoubleParaFromDB("attackd");
                RDouble sideslip =  PI / 180.0 * GlobalDataBase::GetDoubleParaFromDB("angleSlide");

                direction_inlet[0] = cos(attack) * cos(sideslip);
                direction_inlet[1] = sin(attack) * cos(sideslip);
                direction_inlet[2] = sin(sideslip);
            }
        }
    }

    bcParamDB->UpdateData("direction_inlet", &direction_inlet, PHDOUBLE, 3);
    primitiveVar[IU] = localNonDimensionalVelocity * direction_inlet[0];
    primitiveVar[IV] = localNonDimensionalVelocity * direction_inlet[1];
    primitiveVar[IW] = localNonDimensionalVelocity * direction_inlet[2];
}

void GlobalBoundaryCondition::SetSpeedDrictionForReferenceVar(RDouble localNonDimensionalVelocity, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    using namespace IDX;
    using namespace INIFLOW_SPACE;
  /*const int CUSTOM_ANGLE           = 0;
    const int CUSTOM_VECTOR          = 1;
  */
    int directionMethod = CUSTOM_ANGLE;
    if (GlobalDataBase::IsExist("directionMethod", PHINT, 1))
    {
        GlobalDataBase::GetData("directionMethod", &directionMethod, PHINT, 1);
    }
    else
    {
        GlobalDataBase::UpdateData("directionMethod", &directionMethod, PHINT, 1);
    }

    RDouble attackd = 0.0, angleSlide = 0.0;
    RDouble direction_inlet[3] = {1.0,0.0,0.0};

    if(directionMethod == CUSTOM_VECTOR )  //! VECTOR
    {
        if (GlobalDataBase::IsExist("direction_inlet", PHDOUBLE, 3))
        {
            GlobalDataBase::GetData("direction_inlet", &direction_inlet, PHDOUBLE, 3);

            RDouble sum = sqrt( direction_inlet[0] * direction_inlet[0] + direction_inlet[1] * direction_inlet[1] + direction_inlet[2] * direction_inlet[2]);
            if(sum < SMALL)
            {
                direction_inlet[0] = 1.0;
                direction_inlet[1] = 0.0;
                direction_inlet[2] = 0.0;
            }
            else
            {
                sum = 1.0 / sum;
                direction_inlet[0] *= sum;
                direction_inlet[1] *= sum;
                direction_inlet[2] *= sum;
            }
        }

        RDouble sideslip = asin( min(one, direction_inlet[2]) );
        RDouble attack   = asin( min(one, direction_inlet[1] / max(SMALL, cos(sideslip))) );

        GlobalDataBase::UpdateData("attack", &attack, PHDOUBLE, 1);
        GlobalDataBase::UpdateData("sideslip", &sideslip, PHDOUBLE, 1);

        attackd = attack * 180.0 / PI;
        angleSlide = sideslip * 180.0 / PI;

        GlobalDataBase::UpdateData("attackd", &attackd, PHDOUBLE, 1);
        GlobalDataBase::UpdateData("angleSlide", &angleSlide, PHDOUBLE, 1);
    }
    else  //! Attackd and angleSlide
    {
        if (GlobalDataBase::IsExist("attackd", PHDOUBLE, 1))
        {
            GlobalDataBase::GetData("attackd", &attackd, PHDOUBLE, 1);
        }
        else
        {
            GlobalDataBase::UpdateData("attackd", &attackd, PHDOUBLE, 1);
        }

        if (GlobalDataBase::IsExist("angleSlide", PHDOUBLE, 1))
        {
            GlobalDataBase::GetData("angleSlide", &angleSlide, PHDOUBLE, 1);
        }
        else
        {
            GlobalDataBase::UpdateData("angleSlide", &angleSlide, PHDOUBLE, 1);
        }
        RDouble attack   = attackd * PI / 180.0;
        RDouble sideslip = angleSlide * PI / 180.0;

        GlobalDataBase::UpdateData("attack", &attack, PHDOUBLE, 1);
        GlobalDataBase::UpdateData("sideslip", &sideslip, PHDOUBLE, 1);

        direction_inlet[0] = cos(attack) * cos(sideslip);
        direction_inlet[1] = sin(attack) * cos(sideslip);
        direction_inlet[2] = sin(sideslip);
    }

    GlobalDataBase::UpdateData("direction_inlet", &direction_inlet, PHDOUBLE, 3);

    primitiveVar[IU] = localNonDimensionalVelocity * direction_inlet[0];
    primitiveVar[IV] = localNonDimensionalVelocity * direction_inlet[1];
    primitiveVar[IW] = localNonDimensionalVelocity * direction_inlet[2];
}

void GlobalBoundaryCondition::SetBCInitialValuesByDensity(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    RDouble localDimensionalVelocity = 0, localDimensionalDensity = 0, localDimensionalPressure = 0;
    RDouble localDimensionalTemperature = 0, localVibrationTemperature = 0;
    //! Obtain the local parameters on boundary.
    if (bcParamDB->IsExist("refDimensionalVelocity", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalVelocity", &localDimensionalVelocity, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    }
    if (bcParamDB->IsExist("refDimensionalDensity", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalDensity", &localDimensionalDensity, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    }
    if (bcParamDB->IsExist("refDimensionalTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalTemperature", &localDimensionalTemperature, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    }

    int nTvChange = GlobalDataBase::GetIntParaFromDB("nTvChange");
    if(nTvChange == 0)
    {
        if (bcParamDB->IsExist("refDimensionalVibrationTemperature", PHDOUBLE, 1))
        {
            bcParamDB->GetData("refDimensionalVibrationTemperature", &localVibrationTemperature, PHDOUBLE, 1);
        }
        else
        {
            localVibrationTemperature = GlobalDataBase::GetDoubleParaFromDB("freestream_vibration_temperature");
        }
    }

    //! Obtain the reference parameters.
    RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

    int nchem = gas->GetChemicalType();
    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    if (bcParamDB->IsExist("refGama", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refGama", &refGama, PHDOUBLE, 1);
    }

    if (nchem == 0)
    {
        RDouble reference_general_gas_constant = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");
        localDimensionalPressure = localDimensionalDensity * reference_general_gas_constant * localDimensionalTemperature;
    }
    else
    {
        int nNSEqn = gas->GetNSEquationNumber();
        int numberOfSpecies = gas->GetNumberOfSpecies();
        RDouble initMassFraction[MAX_SPECIES_NUM] = {0}, primVars[MAX_SPECIES_NUM] = {0};
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->MassFractionConversion(initMassFraction);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *mass_fraction_reference = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = mass_fraction_reference[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }

        RDouble oMass = 0.0, gasR = 0.0;
        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
        gasR = rjmk * oMass;

        localDimensionalPressure = localDimensionalDensity * gasR * localDimensionalTemperature;
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }

        RDouble cp = 0.0, cv = 0.0, speciesCv[MAX_SPECIES_NUM] = {0};
        int ntmodel = gas->GetTemperatureModel();
        if (ntmodel > 1)    //! multi-temperature model.
        {
            gas->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cv);
            cp = cv + gasR;
            refGama = cp / cv;

            RDouble nonDimTv = localVibrationTemperature / refDimensionalTemperature;
            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            electronEnergy  = electronEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
        else    //! one-temperature model.
        {
            RDouble nonDimTtr = localDimensionalTemperature / refDimensionalTemperature;
            gas->ComputeSpeciesConstantPressureSpecificHeatDimensional(nonDimTtr, speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cp);
            cv = cp - gasR;
            refGama = cp / cv;
        }
    }

    using namespace IDX;
    primitiveVar[IR] = localDimensionalDensity / refDimensionalDensity;

    SetSpeedDrictionForBC(bcParamDB,localDimensionalVelocity / refDimensionalVelocity, primitiveVar);

    primitiveVar[IP] = localDimensionalPressure / (refDimensionalDensity * refDimensionalVelocity * refDimensionalVelocity);

    bcParamDB->UpdateData("refGama", &refGama, PHDOUBLE, 1);
}

void GlobalBoundaryCondition::SetBCInitialValuesByPressure(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    RDouble localDimensionalVelocity = 0, localDimensionalDensity = 0, localDimensionalPressure = 0;
    RDouble localDimensionalTemperature = 0, localVibrationTemperature = 0;
    //! Obtain the local Velocity.
    if (bcParamDB->IsExist("refDimensionalVelocity", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalVelocity", &localDimensionalVelocity, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    }

    //! Obtain the local Pressure.
    if (bcParamDB->IsExist("refDimensionalPressure", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalPressure", &localDimensionalPressure, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");
    }

    if (GlobalDataBase::IsExist("sprayFactor", PHDOUBLE, 1))
    {
        RDouble sprayFactor = GlobalDataBase::GetDoubleParaFromDB("sprayFactor");
        localDimensionalPressure *= (1.0 + sprayFactor);
    }

    //! Obtain the local Temperature.
    if (bcParamDB->IsExist("refDimensionalTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalTemperature", &localDimensionalTemperature, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    }

    //! Obtain the local VibrationTemperature.
    int nTvChange = GlobalDataBase::GetIntParaFromDB("nTvChange");
    if(nTvChange == 0)
    {
        if (bcParamDB->IsExist("refDimensionalVibrationTemperature", PHDOUBLE, 1))
        {
            bcParamDB->GetData("refDimensionalVibrationTemperature", &localVibrationTemperature, PHDOUBLE, 1);
        }
        else
        {
            localVibrationTemperature = GlobalDataBase::GetDoubleParaFromDB("freestream_vibration_temperature");
        }
    }

    //! Obtain the reference parameters.
    RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

    int nchem = gas->GetChemicalType();
    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    if (bcParamDB->IsExist("refGama", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refGama", &refGama, PHDOUBLE, 1);
    }
    if (nchem == 0)
    {
        RDouble reference_general_gas_constant = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");
        localDimensionalDensity = localDimensionalPressure / reference_general_gas_constant / localDimensionalTemperature;
    }
    else
    {
        int nNSEqn = gas->GetNSEquationNumber();
        int numberOfSpecies = gas->GetNumberOfSpecies();
        RDouble initMassFraction[MAX_SPECIES_NUM] = {0}, primVars[MAX_SPECIES_NUM] = {0};
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->MassFractionConversion(initMassFraction);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *mass_fraction_reference = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = mass_fraction_reference[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }

        RDouble oMass = 0.0, gasR = 0.0;
        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
        gasR = rjmk * oMass;

        localDimensionalDensity = localDimensionalPressure / (gasR * localDimensionalTemperature);

        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }

        int ntmodel = gas->GetTemperatureModel();
        RDouble cp = 0.0, cv = 0.0, speciesCv[MAX_SPECIES_NUM] = {0};

        if (ntmodel > 1)    //! multi-temperature model.
        {
            gas->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cv);
            cp = cv + gasR;
            refGama = cp / cv;

            RDouble nonDimTv = localVibrationTemperature / refDimensionalTemperature;
            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            electronEnergy  = electronEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
        else    //! one-temperature model.
        {
            RDouble nonDimTtr = localDimensionalTemperature / refDimensionalTemperature;
            gas->ComputeSpeciesConstantPressureSpecificHeatDimensional(nonDimTtr, speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cp);
            cv = cp - gasR;
            refGama = cp / cv;
        }
    }

    using namespace IDX;
    primitiveVar[IR] = localDimensionalDensity / refDimensionalDensity;

    SetSpeedDrictionForBC(bcParamDB,localDimensionalVelocity / refDimensionalVelocity, primitiveVar);

    primitiveVar[IP] = localDimensionalPressure / (refDimensionalDensity * refDimensionalVelocity * refDimensionalVelocity);

    bcParamDB->UpdateData("refGama", &refGama, PHDOUBLE, 1);
}

void GlobalBoundaryCondition::SetBCInitialValuesByTotalPressure(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    RDouble localMachNumber = 0, localDimensionalDensity = 0, localDimensionalPressure = 0, localDimensionalVelocity = 0;
    RDouble localDimensionalTemperature = 0, localVibrationTemperature = 0;
    RDouble totalPressure = 0, totalTemperature = 0;
    //! Obtain the local Velocity.
    if (bcParamDB->IsExist("refMachNumber", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refMachNumber", &localMachNumber, PHDOUBLE, 1);
    }
    else
    {
        localMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    }

    //! Obtain the local Pressure.
    if (bcParamDB->IsExist("totalPressure", PHDOUBLE, 1))
    {
        bcParamDB->GetData("totalPressure", &totalPressure, PHDOUBLE, 1);
    }
    else
    {
        totalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");
    }
    //! Obtain the local Temperature.
    if (bcParamDB->IsExist("totalTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("totalTemperature", &totalTemperature, PHDOUBLE, 1);
    }
    else
    {
        totalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    }

    //! Obtain the reference parameters.
    RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

    int nchem = gas->GetChemicalType();
    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    if (bcParamDB->IsExist("refGama", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refGama", &refGama, PHDOUBLE, 1);
    }

    if (nchem == 0)
    {
        RDouble reference_general_gas_constant = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");

        RDouble localDimensionalSonicSpeed = sqrt(refGama * reference_general_gas_constant * localDimensionalTemperature);
        localDimensionalVelocity = localDimensionalSonicSpeed * localMachNumber;

        localDimensionalTemperature = totalTemperature /(1.0 + 0.5 * (refGama - 1.0) * localMachNumber * localMachNumber);
        localDimensionalPressure = totalPressure * pow(localDimensionalTemperature / totalTemperature, refGama /(refGama - 1.0));

        localDimensionalDensity = localDimensionalPressure / (reference_general_gas_constant * localDimensionalTemperature);

    }
    else
    {
        int nNSEqn = gas->GetNSEquationNumber();
        int numberOfSpecies = gas->GetNumberOfSpecies();

        RDouble initMassFraction[MAX_SPECIES_NUM] = {0}, primVars[MAX_REACTION_NUM] = {0};
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->MassFractionConversion(initMassFraction);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *mass_fraction_reference = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = mass_fraction_reference[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }

        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }

        gas->ComputeTemperatureFromTotalTemperature(initMassFraction, totalTemperature, localMachNumber, localDimensionalTemperature, refGama, 1);  // localDimensionalTemperature  and  refGama

        localDimensionalPressure = totalPressure * pow(localDimensionalTemperature / totalTemperature, refGama /(refGama - 1.0));

        int nTvChange = GlobalDataBase::GetIntParaFromDB("nTvChange");
        if(nTvChange == 0)
        {
            if (bcParamDB->IsExist("refDimensionalVibrationTemperature", PHDOUBLE, 1))
            {
                bcParamDB->GetData("refDimensionalVibrationTemperature", &localVibrationTemperature, PHDOUBLE, 1);
            }
            else
            {
                localVibrationTemperature = localDimensionalTemperature;
            }
        }

        RDouble oMass = 0.0, gasR = 0.0;
        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
        gasR = rjmk * oMass;

        localDimensionalDensity = localDimensionalPressure / (gasR * localDimensionalTemperature);

        localDimensionalVelocity = localMachNumber * sqrt(refGama * gasR * localDimensionalTemperature);

        int ntmodel = gas->GetTemperatureModel();
        if (ntmodel > 1)    //! multi-temperature model.
        {
            RDouble nonDimTv = localVibrationTemperature / refDimensionalTemperature;
            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            electronEnergy  = electronEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
    }

    using namespace IDX;
    primitiveVar[IR] = localDimensionalDensity / refDimensionalDensity;

    SetSpeedDrictionForBC(bcParamDB,localDimensionalVelocity / refDimensionalVelocity, primitiveVar);

    primitiveVar[IP] = localDimensionalPressure / (refDimensionalDensity * refDimensionalVelocity * refDimensionalVelocity);

    RDouble localNonDimensionalTemperature = localDimensionalTemperature / refDimensionalTemperature;

    bcParamDB->UpdateData("refGama", &refGama, PHDOUBLE, 1);
    bcParamDB->UpdateData("localNonDimensionalTemperature", &localNonDimensionalTemperature, PHDOUBLE, 1);

}

void GlobalBoundaryCondition::SetBCInitialValuesByMachNumberTemperaturePressure(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    RDouble localMachNumber = 0, localDimensionalDensity = 0, localDimensionalPressure = 0;
    RDouble localDimensionalTemperature = 0, localVibrationTemperature = 0;
    //! Obtain the local parameters on boundary.
    if (bcParamDB->IsExist("refMachNumber", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refMachNumber", &localMachNumber, PHDOUBLE, 1);
    }
    else
    {
        localMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    }
    if (bcParamDB->IsExist("refDimensionalPressure", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalPressure", &localDimensionalPressure, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");
    }
    
    if (bcParamDB->IsExist("refDimensionalTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalTemperature", &localDimensionalTemperature, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    }

    int nTvChange = GlobalDataBase::GetIntParaFromDB("nTvChange");
    if(nTvChange == 0)
    {
        if (bcParamDB->IsExist("refDimensionalVibrationTemperature", PHDOUBLE, 1))
        {
            bcParamDB->GetData("refDimensionalVibrationTemperature", &localVibrationTemperature, PHDOUBLE, 1);
        }
        else
        {
            localVibrationTemperature = GlobalDataBase::GetDoubleParaFromDB("freestream_vibration_temperature");
        }
    }

    //! Obtain the reference parameters.
    RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    RDouble localReNumber, localDimensionalVelocity;
    RDouble nonDimTtr = localDimensionalTemperature / refDimensionalTemperature;
    RDouble nonDimTv  = localVibrationTemperature / refDimensionalTemperature;

    int nchem = gas->GetChemicalType();
    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    if (bcParamDB->IsExist("refGama", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refGama", &refGama, PHDOUBLE, 1);
    }

    if (nchem == 0)
    {
        RDouble t0i  = 288.15;
        RDouble tsi  = 110.4;
        RDouble mu0i = 1.7894E-05;

        RDouble local_viscosity_dimensional = (t0i + tsi) / (localDimensionalTemperature + tsi) * pow(localDimensionalTemperature / t0i, 1.5) * mu0i;
        
        RDouble reference_general_gas_constant = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");

        RDouble localDimensionalSonicSpeed = sqrt(refGama * reference_general_gas_constant * localDimensionalTemperature);
        localDimensionalVelocity = localDimensionalSonicSpeed * localMachNumber;

        localDimensionalDensity = localDimensionalPressure / reference_general_gas_constant / localDimensionalTemperature;
        localReNumber = localDimensionalDensity * localDimensionalVelocity / local_viscosity_dimensional;
    }
    else
    {
        int nNSEqn = gas->GetNSEquationNumber();
        int numberOfSpecies = gas->GetNumberOfSpecies();
        RDouble initMassFraction[MAX_SPECIES_NUM] = {0}, speciesCv[MAX_SPECIES_NUM] = {0}, primVars[MAX_SPECIES_NUM] = {0}, speciesViscosity[MAX_SPECIES_NUM] = {0};
        RDouble speciesWeight[MAX_SPECIES_NUM] = {0}, moleFractions[MAX_SPECIES_NUM] = {0};
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->MassFractionConversion(initMassFraction);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *mass_fraction_reference = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = mass_fraction_reference[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }

        RDouble oMass = 0.0, gasR = 0.0;
        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
        gas->ComputeMoleFractionByMassFraction(initMassFraction, moleFractions);
        gasR = rjmk * oMass;

        RDouble cp = 0.0, cv = 0.0;
        int ntmodel = gas->GetTemperatureModel();
        if (ntmodel > 1)    //! multi-temperature model.
        {
            gas->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cv);
            cp = cv + gasR;
            refGama = cp / cv;
            gas->ComputeSpeciesDimensionalViscosityWithCurveFitMethod(nonDimTtr, nonDimTv, speciesViscosity);

            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            electronEnergy  = electronEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
        else    //! one-temperature model.
        {
            gas->ComputeSpeciesConstantPressureSpecificHeatDimensional(nonDimTtr, speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cp);
            cv = cp - gasR;
            refGama = cp / cv;
            gas->ComputeOneTemperatureModelSpeciesViscosityDimensional(nonDimTtr, speciesViscosity);
        }

        RDouble localDimensionalSonicSpeed = sqrt(refGama * gasR * localDimensionalTemperature);
        localDimensionalVelocity = localDimensionalSonicSpeed * localMachNumber;

        RDouble *MolecularWeight = gas->GetMolecularWeightDimensional();
        gas->GetPartitionFunctionValues(speciesViscosity, MolecularWeight, moleFractions, speciesWeight);
        RDouble local_viscosity_dimensional = gas->GetMixedGasViscosityWithWilkeFormula(moleFractions, speciesViscosity, speciesWeight);

        localDimensionalDensity = localDimensionalPressure / gasR / localDimensionalTemperature;
        localReNumber = localDimensionalDensity * localDimensionalVelocity / local_viscosity_dimensional;

        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }
    }
    using namespace IDX;
    primitiveVar[IR] = localDimensionalDensity / refDimensionalDensity;

    SetSpeedDrictionForBC(bcParamDB,localDimensionalVelocity / refDimensionalVelocity, primitiveVar);

    primitiveVar[IP] = localDimensionalPressure / (refDimensionalDensity * refDimensionalVelocity * refDimensionalVelocity);

    bcParamDB->UpdateData("refGama", &refGama, PHDOUBLE, 1);
}

void GlobalBoundaryCondition::SetBCInitialValuesByWRFData(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    string localFilePath;
    RDouble localLongitude = 0, localLatitude = 0, localVibrationTemperature = 0;
    //! Obtain the local parameters on boundary.
    if (bcParamDB->IsExist("weatherDataFilePath", PHSTRING, 1))
    {
        bcParamDB->GetData("weatherDataFilePath", &localFilePath, PHSTRING, 1);
    }
    else
    {
        localFilePath = GlobalDataBase::GetStrParaFromDB("weatherDataFilePath");
    }
    if (bcParamDB->IsExist("longitude", PHDOUBLE, 1))
    {
        bcParamDB->GetData("longitude", &localLongitude, PHDOUBLE, 1);
    }
    else
    {
        localLongitude = GlobalDataBase::GetDoubleParaFromDB("longitude");
    }

    if (bcParamDB->IsExist("latitude", PHDOUBLE, 1))
    {
        bcParamDB->GetData("latitude", &localLatitude, PHDOUBLE, 1);
    }
    else
    {
        localLatitude = GlobalDataBase::GetDoubleParaFromDB("latitude");
    }

    RDouble localDimensionalDensity, localDimensionalVelocity, localDimensionalPressure, localDimensionalTemperature, localDimensionalSonicSpeed,localAngleSlide;

    gas->ReadAirInfo(localFilePath, localLongitude, localLatitude, localDimensionalTemperature, localDimensionalPressure, localDimensionalDensity,localAngleSlide, localDimensionalVelocity, localDimensionalSonicSpeed);

    int nTvChange = GlobalDataBase::GetIntParaFromDB("nTvChange");
    if(nTvChange == 0)
    {
        if (bcParamDB->IsExist("refDimensionalVibrationTemperature", PHDOUBLE, 1))
        {
            bcParamDB->GetData("refDimensionalVibrationTemperature", &localVibrationTemperature, PHDOUBLE, 1);
        }
        else
        {
            localVibrationTemperature = GlobalDataBase::GetDoubleParaFromDB("freestream_vibration_temperature");
        }
    }
   

    //! Obtain the reference parameters.
    RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

    int nchem = gas->GetChemicalType();
    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    if (bcParamDB->IsExist("refGama", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refGama", &refGama, PHDOUBLE, 1);
    }
    if (nchem == 0)
    {
        RDouble reference_general_gas_constant = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");
        localDimensionalDensity = localDimensionalPressure / reference_general_gas_constant / localDimensionalTemperature;
    }
    else
    {
        int nNSEqn = gas->GetNSEquationNumber();
        int numberOfSpecies = gas->GetNumberOfSpecies();
        RDouble initMassFraction[MAX_SPECIES_NUM] = {0}, primVars[MAX_SPECIES_NUM] = {0};
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->MassFractionConversion(initMassFraction);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *mass_fraction_reference = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = mass_fraction_reference[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }

        RDouble oMass = 0.0, gasR = 0.0;
        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
        gasR = rjmk * oMass;

        localDimensionalDensity = localDimensionalPressure / (gasR * localDimensionalTemperature);

        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }

        int ntmodel = gas->GetTemperatureModel();
        RDouble cp = 0.0, cv = 0.0, speciesCv[MAX_SPECIES_NUM] = {0};

        if (ntmodel > 1)    //! multi-temperature model.
        {
            gas->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cv);
            cp = cv + gasR;
            refGama = cp / cv;

            RDouble nonDimTv = localVibrationTemperature / refDimensionalTemperature;
            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            electronEnergy  = electronEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
        else    //! one-temperature model.
        {
            RDouble nonDimTtr = localDimensionalTemperature / refDimensionalTemperature;
            gas->ComputeSpeciesConstantPressureSpecificHeatDimensional(nonDimTtr, speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cp);
            cv = cp - gasR;
            refGama = cp / cv;
        }
    }

    using namespace IDX;
    primitiveVar[IR] = localDimensionalDensity / refDimensionalDensity;

    SetSpeedDrictionForBC(bcParamDB,localDimensionalVelocity / refDimensionalVelocity, primitiveVar);

    primitiveVar[IP] = localDimensionalPressure / (refDimensionalDensity * refDimensionalVelocity * refDimensionalVelocity);

    bcParamDB->UpdateData("refGama", &refGama, PHDOUBLE, 1);
}

void GlobalBoundaryCondition::SetBCInitialValuesByPrimitive(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    using namespace IDX;
    bool bFlag = false;
    RDouble primDensity = 0, primU = 0 , primV = 0, primW = 0, primPressure = 0;
    //! Obtain the local parameters on boundary.
    if (bcParamDB->IsExist("primDensity", PHDOUBLE, 1))
    {
        bcParamDB->GetData("primDensity", &primDensity, PHDOUBLE, 1);
    }
    else
    {
        bFlag = true;
    }
    if (bcParamDB->IsExist("primU", PHDOUBLE, 1))
    {
        bcParamDB->GetData("primU", &primU, PHDOUBLE, 1);
    }
    else
    {
        bFlag = true;
    }
    if (bcParamDB->IsExist("primV", PHDOUBLE, 1))
    {
        bcParamDB->GetData("primV", &primV, PHDOUBLE, 1);
    }
    else
    {
        bFlag = true;
    }
    if (bcParamDB->IsExist("primW", PHDOUBLE, 1))
    {
        bcParamDB->GetData("primW", &primW, PHDOUBLE, 1);
    }
    else
    {
        bFlag = true;
    }
    if (bcParamDB->IsExist("primPressure", PHDOUBLE, 1))
    {
        bcParamDB->GetData("primPressure", &primPressure, PHDOUBLE, 1);
    }
    else
    {
        bFlag = true;
    }

    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    if (bcParamDB->IsExist("refGama", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refGama", &refGama, PHDOUBLE, 1);
    }

    RDouble initMassFraction[MAX_SPECIES_NUM] = {0}, primVars[MAX_SPECIES_NUM] = {0};
    int nchem = gas->GetChemicalType();
    int nNSEqn = gas->GetNSEquationNumber();
    int numberOfSpecies = gas->GetNumberOfSpecies();
    int ntmodel = gas->GetTemperatureModel();
    int nEquation = nNSEqn;
    if (nchem > 0 && bFlag == false)
    {
        nEquation = nNSEqn + numberOfSpecies + ntmodel - 1;

        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->MassFractionConversion(initMassFraction);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *mass_fraction_reference = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = mass_fraction_reference[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }

        RDouble oMass = 0.0, gasR = 0.0;
        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
        gasR = rjmk * oMass;

        RDouble cp = 0.0, cv = 0.0, speciesCv[MAX_SPECIES_NUM] = {0};

        RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
        RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

        RDouble squareVelocity = refDimensionalVelocity * refDimensionalVelocity;
        RDouble dimDensity = primDensity * refDimensionalDensity;
        RDouble dimPressure = primPressure * refDimensionalDensity * squareVelocity;
        RDouble dimTemperature = dimPressure / (gasR * dimDensity);
        RDouble nonDimTv = dimTemperature / refDimensionalTemperature;

        if (ntmodel > 1)    //! multi-temperature model.
        {
            gas->ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cv);
            cp = cv + gasR;
            refGama = cp / cv;

            //! Obtain the reference parameters.
            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / squareVelocity;
            electronEnergy  = electronEnergy  / squareVelocity;
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
        else    //! one-temperature model.
        {
            gas->ComputeSpeciesConstantPressureSpecificHeatDimensional(nonDimTv, speciesCv);
            gas->ComputeMixtureByMassFraction(initMassFraction, speciesCv, cp);
            cv = cp - gasR;
            refGama = cp / cv;
        }
    }

    if (bFlag)
    {
        RDouble *primitiveVarFarfield = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("prim_inf"));
        for (int s = 0; s < nEquation; ++ s)
        {
            primitiveVar[s] = primitiveVarFarfield[s];
        }
    }
    else
    {
        primitiveVar[IR] = primDensity;
        primitiveVar[IU] = primU;
        primitiveVar[IV] = primV;
        primitiveVar[IW] = primW;
        primitiveVar[IP] = primPressure;
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }
    }
    bcParamDB->UpdateData("refGama", &refGama, PHDOUBLE, 1);
}

void GlobalBoundaryCondition::SetBCInitialValuesByMolecularWeightGama(Data_Param *bcParamDB, RDouble *primitiveVar)
{
    using namespace GAS_SPACE;
    RDouble localDimensionalVelocity, localDimensionalPressure, localDimensionalDensity;
    RDouble localDimensionalTemperature, localVibrationTemperature;
    //! Obtain the local parameters on boundary.
    if (bcParamDB->IsExist("refDimensionalVelocity", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalVelocity", &localDimensionalVelocity, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    }
    if (bcParamDB->IsExist("refDimensionalPressure", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalPressure", &localDimensionalPressure, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");
    }

    if (bcParamDB->IsExist("refDimensionalTemperature", PHDOUBLE, 1))
    {
        bcParamDB->GetData("refDimensionalTemperature", &localDimensionalTemperature, PHDOUBLE, 1);
    }
    else
    {
        localDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    }
    
    int nTvChange = GlobalDataBase::GetIntParaFromDB("nTvChange");
    if(nTvChange == 0)
    {
        if (bcParamDB->IsExist("refDimensionalVibrationTemperature", PHDOUBLE, 1))
        {
            bcParamDB->GetData("refDimensionalVibrationTemperature", &localVibrationTemperature, PHDOUBLE, 1);
        }
        else
        {
            localVibrationTemperature = GlobalDataBase::GetDoubleParaFromDB("freestream_vibration_temperature");
        }
    }
   

    //! Obtain the reference parameters.
    RDouble refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    RDouble refDimensionalVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

    int nchem = gas->GetChemicalType();
    if (nchem == 0)
    {
        RDouble reference_general_gas_constant = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");
        localDimensionalDensity = localDimensionalPressure / reference_general_gas_constant / localDimensionalTemperature;
    }

    else
    {
        int nNSEqn = gas->GetNSEquationNumber();
        int numberOfSpecies = gas->GetNumberOfSpecies();
        RDouble initMassFraction[MAX_SPECIES_NUM], primVars[MAX_SPECIES_NUM];
        RDouble localGama;
        if (bcParamDB->IsExist("initMassFraction", PHDOUBLE, numberOfSpecies))
        {
            bcParamDB->GetData("initMassFraction", &initMassFraction, PHDOUBLE, numberOfSpecies);
            gas->SetMaximumSpecies(initMassFraction);
        }
        else
        {
            RDouble *mass_fraction_reference = gas->GetInitMassFraction();
            for (int s = 0; s < numberOfSpecies; ++ s)
            {
                initMassFraction[s] = mass_fraction_reference[s];
            }
        }
        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primVars[nNSEqn + s] = initMassFraction[s];
        }
        
        if (bcParamDB->IsExist("refGama", PHDOUBLE, 1))
        {
            bcParamDB->GetData("refGama", &localGama, PHDOUBLE, 1);
        }
        else
        {
            
        }

        RDouble oMass = 0.0, gasR = 0.0;
        gas->ComputeMolecularWeightReciprocalDimensional(primVars, oMass);
        gasR = rjmk * oMass;

        localDimensionalDensity = localDimensionalPressure / (gasR * localDimensionalTemperature);

        for (int s = 0; s < numberOfSpecies; ++ s)
        {
            primitiveVar[nNSEqn + s] = initMassFraction[s];
        }

        int ntmodel = gas->GetTemperatureModel();
        if (ntmodel > 1)    //! multi-temperature model.
        {
            RDouble nonDimTv = localVibrationTemperature / refDimensionalTemperature;
            RDouble vibrationEnergy = gas->GetMixedGasDimensionalVibrationEnergy(nonDimTv, initMassFraction);
            RDouble electronEnergy  = gas->GetMixedGasDimensionalElectronEnergy(nonDimTv, initMassFraction);
            vibrationEnergy = vibrationEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            electronEnergy  = electronEnergy / (refDimensionalVelocity * refDimensionalVelocity);
            if (ntmodel == 2)
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy + electronEnergy;
            }
            else
            {
                primitiveVar[nNSEqn + numberOfSpecies] = vibrationEnergy;
                primitiveVar[nNSEqn + numberOfSpecies + 1] = electronEnergy;
            }
        }
    }

    using namespace IDX;
    primitiveVar[IR] = localDimensionalDensity / refDimensionalDensity;

    SetSpeedDrictionForBC(bcParamDB,localDimensionalVelocity / refDimensionalVelocity, primitiveVar);

    primitiveVar[IP] = localDimensionalPressure / (refDimensionalDensity * refDimensionalVelocity * refDimensionalVelocity);
}

void GlobalBoundaryCondition::SetBCDataBaseByGlobalBC()
{
    vector <SimpleBC *> *globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();
    if (globalBCList == 0) return;

    using namespace PHMPI;

    int nZones = PHMPI::GetNumberofGlobalZones();
    int myid   = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int procID = PHMPI::GetZoneProcessorID(iZone);

        if (myid == procID)
        {
            Grid *gridofThisProcessor = GetGrid(iZone, 0);

            while (gridofThisProcessor)
            {
                int gridtype = gridofThisProcessor->Type();

                if (gridtype == PHSPACE::STRUCTGRID)
                {
                    StructGrid *grid = StructGridCast(gridofThisProcessor);
                    StructBCSet *structBCSet = grid->GetStructBCSet();
                    int nBCRegion = structBCSet->GetnBCRegion();
                    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
                    {
                        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
                        string bcRegionName = bcregion->GetBCName();
                        int bcRegionType = bcregion->GetBCType();
    
                        vector <SimpleBC *> :: iterator iter;
                        for (iter = globalBCList->begin(); iter != globalBCList->end(); iter ++)
                        {
                            string globalBCName = (*iter)->GetBCName();
                            int globalBCType = (*iter)->GetBCType();
                            Data_Param *globalBCParam = (*iter)->GetBCParamDataBase();
                            Data_Param *bcParamDataBase = new Data_Param(*globalBCParam);
                            if (globalBCName == bcRegionName && bcRegionType == globalBCType)
                            {
                                bcregion->SetBCParamDataBase(bcParamDataBase);
                                break;
                            }
                        }
                    }
                }
                else if (gridtype == PHSPACE::UNSTRUCTGRID)
                {
                    UnstructGrid *grid = UnstructGridCast(gridofThisProcessor);

                    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
                    int nBCRegion = unstructBCSet->GetnBCRegion();
                    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++iBCRegion)
                    {
                        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
                        string bcRegionName = bcRegion->GetBCName();
                        int bcRegionType = bcRegion->GetBCType();

                        vector <SimpleBC *> :: iterator iter;
                        for (iter = globalBCList->begin(); iter != globalBCList->end(); iter ++)
                        {
                            string globalBCName = (*iter)->GetBCName();
                            int globalBCType = (*iter)->GetBCType();
                            Data_Param *globalBCParam = (*iter)->GetBCParamDataBase();
                            Data_Param *bcParamDataBase = new Data_Param(*globalBCParam);
                            if (globalBCName == bcRegionName && bcRegionType == globalBCType)
                            {
                                bcRegion->SetBCParamDataBase(bcParamDataBase);
                                break;
                            }
                        }
                    }

                } 
                gridofThisProcessor = gridofThisProcessor->GetCoarseGrid();
            }
        }
    }
}

void GlobalBoundaryCondition::ChangeBCTypeByGlobalBC()
{
    vector<SimpleBC *> *globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();
    if (globalBCList == 0) return;

    set< pair<string, int> > bcnameMap;
    vector<SimpleBC *>::iterator iter1;
    for (iter1 = globalBCList->begin(); iter1 != globalBCList->end(); ++ iter1)
    {
        SimpleBC *boundaryCondition = *iter1;
        bcnameMap.insert(pair<string, int>(boundaryCondition->GetBCName(), boundaryCondition->GetBCType()));
    }

    using namespace PHMPI;

    int nZones = PHMPI::GetNumberofGlobalZones();
    int myid   = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int procID = PHMPI::GetZoneProcessorID(iZone);

        if (myid == procID)
        {
            Grid *gridofThisProcessor = GetGrid(iZone, 0);

            int gridtype = gridofThisProcessor->Type();

            if (gridtype == PHSPACE::STRUCTGRID)
            {
                StructGrid *grid = StructGridCast(gridofThisProcessor);

                StructBCSet *structBCSet = grid->GetStructBCSet();
                int nBCRegion = structBCSet->GetnBCRegion();
                for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
                {
                    StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
                    string bcname = bcregion->GetBCName();

                    set< pair<string, int> >::iterator iter2;
                    for (iter2 = bcnameMap.begin(); iter2 != bcnameMap.end(); ++ iter2)
                    {
                        if ((*iter2).first != bcname) continue;

                        bcregion->SetBCType((*iter2).second);
                        break;
                    }
                }
            }
            else if (gridtype == PHSPACE::UNSTRUCTGRID)
            {
                UnstructGrid *grid = UnstructGridCast(gridofThisProcessor);

                UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
                if (!unstructBCSet) continue;
                int nBCRegion = unstructBCSet->GetnBCRegion();

                for (int iBCRegion = 0; iBCRegion < nBCRegion; ++iBCRegion)
                {
                    UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
                    string bcName = bcRegion->GetBCName();

                    set< pair<string, int> >::iterator iter3;
                    for (iter3 = bcnameMap.begin(); iter3 != bcnameMap.end(); ++ iter3)
                    {
                        if ((*iter3).first != bcName) continue;

                        bcRegion->SetBCType((*iter3).second);
                        break;
                    }
                }

                int nBoundFace = grid->GetNBoundFace();
                UnstructBCSet **bcr = grid->GetBCRecord();
                for (int iFace = 0; iFace < nBoundFace; ++iFace)
                {
                    string bcName = bcr[iFace]->GetBCName();

                    set< pair<string, int> >::iterator iter4;
                    for (iter4 = bcnameMap.begin(); iter4 != bcnameMap.end(); ++ iter4)
                    {
                        if ((*iter4).first != bcName) continue;

                        bcr[iFace]->SetKey((*iter4).second);
                        break;
                    }
                }

            }
        }
    }
}

void GlobalBoundaryCondition::ParseBCFromFile(fstream &file)
{
    string line, keyword, name, word;
    string *value = nullptr;
    //string separator  = " =\t\r\n#$,;\"";
    //! \t is tab key.
    string separator = " =\r\n\t#$,;\"";
    //string bcNameSeparator  = "=\r\n\t#$,;\"";
    string sep1 = "\r\n";
    //string section    = "{}";
    string errMsg = "error in boundary condition file";

    map < string, int > keywordmap;

    keywordmap.insert(pair<string, int>("int"   , PHINT  ));
    keywordmap.insert(pair<string, int>("float" , PHFLOAT));
    keywordmap.insert(pair<string, int>("double", PHDOUBLE));
    keywordmap.insert(pair<string, int>("string", PHSTRING));

    int taskType = GetTaskCode();

    int nBoundaryConditions = 0;
    while (!file.eof())
    {
        getline(file, line);
        if (line == "") continue;
        FindNextWord(line, word, sep1);
        if (word.substr(0, 1) == "#" || word.substr(0, 2) == "//") continue;

        line = FindNextWord(line, keyword, separator);
        if (keyword == "") continue;

        line = FindNextWord(line, name, separator);

        int count = 0;
        int arraysize = 1;
        //! Do not consider array here, so count = 0.
        if (count == 0)
        {
            value = new string[1]();
            if (name == "bcName")
            {
                //! For value of bcName, blank space is not used as separator!
                string::size_type firstindex, nextindex;
                firstindex = line.find_first_of("\"");
                nextindex  = line.find_last_of("\"");
                value[0] = line.substr(firstindex + 1, nextindex-firstindex-1);
            }
            else
            {
                line = FindNextWord(line, value[0], separator);
            }
        }
        else
        {
            ostringstream oss;
            oss << errMsg << "\n";
            TK_Exit::ExceptionExit(oss.str(), true);
        }

        int type, size = arraysize;
        type = keywordmap[keyword];
        if (type == PHINT && name != "nBoundaryConditions")
        {
            TK_Exit::ExceptionExit("Error: the first line must be format \"int nBoundaryConditions = INTEGER\"!", true);
        }
        if (type == PHSTRING && name != "bcName")
        {
            TK_Exit::ExceptionExit("Error: BC name must be format string bcName = \" DownWall\" !", true);
        }
        if (type == PHINT)
        {
            int *data = new int [size];
            for (int i = 0 ; i < size; ++ i)
            {
                from_string< int >(data[i], value[i], std::dec);
            }

            nBoundaryConditions = data[0];
            delete [] data;    data = nullptr;

            if (nBoundaryConditions == 0)
            {
                //! To compatible the old version.
                SetDefaultSolidBoundaryCondition();
                delete [] value;    value = nullptr;
                break;
            }
            else
            {
                globalBCList = new vector <SimpleBC *>;
            }
        }
        else if (type == PHSTRING)
        {
            //! Parse out the parameters into data base in each Boundary Condition.
            SimpleBC *boundaryCondition = new SimpleBC();
            boundaryCondition->SetBCName(*value);

            Data_Param *bcParamDataBase = new Data_Param();
            TK_Parse::ReadBasicData(file, bcParamDataBase);
            boundaryCondition->SetBCParamDataBase(bcParamDataBase);

            int bcType;
            bcParamDataBase->GetData("bcType", &bcType, PHINT, 1);
            boundaryCondition->SetBCType(bcType);

            if (bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::ABLATION_SURFACE)
            {
                string bodyName = "body";
                if (bcParamDataBase->CheckDataExist("bodyName"))
                {
                    bcParamDataBase->GetData("bodyName", &bodyName, PHSTRING, 1);
                }
                else
                {
                    bcParamDataBase->UpdateData("bodyName", &bodyName, PHSTRING, 1);
                }

                boundaryCondition->SetBodyName(bodyName);
                if (taskType == SOLVE_FIELD || taskType == POST_PROCESSING)
                {
                    RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

                    //! For solid wall BC, set the default reference geometry scale.
                    RDouble partForceReferenceArea, partForceReferenceLength, partForceReferenceLength_B;
                    RDouble partForceReferenceX, partForceReferenceY, partForceReferenceZ;

                    //! This judgement is not enough, and needs to be improved. In fact, judgements of the five reference variables should all exist.
                    bool isPartReferenceScaleExist = bcParamDataBase->IsExist("forceReferenceArea", PHDOUBLE, 1);

                    if (isPartReferenceScaleExist)
                    {
                        bcParamDataBase->GetData("forceReferenceArea",   &partForceReferenceArea, PHDOUBLE, 1);
                        bcParamDataBase->GetData("forceReferenceLength", &partForceReferenceLength, PHDOUBLE, 1);
                        bcParamDataBase->GetData("forceReferenceLengthSpanWise", &partForceReferenceLength_B, PHDOUBLE, 1);

                        bcParamDataBase->GetData("TorqueRefX", &partForceReferenceX, PHDOUBLE, 1);
                        bcParamDataBase->GetData("TorqueRefY", &partForceReferenceY, PHDOUBLE, 1);
                        bcParamDataBase->GetData("TorqueRefZ", &partForceReferenceZ, PHDOUBLE, 1);
                    }
                    else
                    {
                        partForceReferenceArea     = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
                        partForceReferenceLength   = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
                        partForceReferenceLength_B = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");
                        partForceReferenceX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
                        partForceReferenceY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
                        partForceReferenceZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");
                    }

                    bcParamDataBase->UpdateData("forceReferenceArea", &partForceReferenceArea, PHDOUBLE, 1);
                    bcParamDataBase->UpdateData("forceReferenceLength", &partForceReferenceLength, PHDOUBLE, 1);
                    bcParamDataBase->UpdateData("forceReferenceLengthSpanWise", &partForceReferenceLength_B, PHDOUBLE, 1);
                    bcParamDataBase->UpdateData("TorqueRefX", &partForceReferenceX, PHDOUBLE, 1);
                    bcParamDataBase->UpdateData("TorqueRefY", &partForceReferenceY, PHDOUBLE, 1);
                    bcParamDataBase->UpdateData("TorqueRefZ", &partForceReferenceZ, PHDOUBLE, 1);

                    bool isDumpHingeMomentExist = bcParamDataBase->IsExist("dumpHingeMoment", PHINT, 1);
                    int dumpHingeMoment = 0;
                    RDouble localCoordAxis0[3] = {0, 0, 0};
                    RDouble localCoordAxis1[3] = {0, 0, 0};
                    if (isDumpHingeMomentExist)
                    {
                        bcParamDataBase->GetData("dumpHingeMoment", &dumpHingeMoment, PHINT, 1);
                        if (dumpHingeMoment != 0)
                        {
                            bcParamDataBase->GetData("localCoordAxis0", localCoordAxis0, PHDOUBLE, 3);
                            bcParamDataBase->GetData("localCoordAxis1", localCoordAxis1, PHDOUBLE, 3);
                        }
                    }

                    //! keep consistend in grid unit:
                    for (int idx = 0; idx < 3; idx ++)
                    {
                        localCoordAxis0[idx] /= gridScaleFactor;
                        localCoordAxis1[idx] /= gridScaleFactor;
                    }
                    bcParamDataBase->UpdateData("dumpHingeMoment", &dumpHingeMoment, PHINT, 1);
                    bcParamDataBase->UpdateData("localCoordAxis0", localCoordAxis0, PHDOUBLE, 3);
                    bcParamDataBase->UpdateData("localCoordAxis1", localCoordAxis1, PHDOUBLE, 3);
                }
                else if (taskType == CREATE_GRID)
                {
                    int gridobj = GlobalDataBase::GetIntParaFromDB("gridobj");
                    const int GRID_ADAPTATION = 2;
                    if (gridobj == GRID_ADAPTATION)
                    {
                        string geometryFileName;
                        bool isGeoemetryFileNameExist = bcParamDataBase->IsExist("geometryFileName", PHSTRING, 1);

                        if (!isGeoemetryFileNameExist)
                        {
                            geometryFileName = "";
                            bcParamDataBase->UpdateData("geometryFileName", &geometryFileName, PHSTRING, 1);
                        }
                    }
                }
            }

            globalBCList->push_back(boundaryCondition);
        }

        delete[] value;    value = nullptr;
    }
}

int GlobalBoundaryCondition::GetNumberOfSolidWallPart()
{
    if (globalBCList == 0) return 0;

    int nSolidWallPart = 0;
    vector<SimpleBC *>::iterator iter;
    for (iter = globalBCList->begin(); iter != globalBCList->end(); ++ iter)
    {
        int bcType = (*iter)->GetBCType();
        if (bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::ABLATION_SURFACE)
        {
            ++ nSolidWallPart;
        }
    }

    return nSolidWallPart;
}

bool IsInterface(const int &boundaryConditionType)
{
    if (boundaryConditionType == PHENGLEI::INTERFACE)
    {
        return true;
    }
    return false;
}

map<int,string> CreatFantasybctype()
{
    map<int,string> Fantasybctype;

    Fantasybctype[-3]   = "USER_DEFINED";
    Fantasybctype[-2]   = "WAKE";
    Fantasybctype[-1]   = "INTERFACE";
    Fantasybctype[0]    = "NO_BOUNDARY_CONDITION";
    Fantasybctype[1]    = "EXTRAPOLATION";
    Fantasybctype[2]    = "SOLID_SURFACE";
    Fantasybctype[3]    = "SYMMETRY";
    Fantasybctype[4]    = "FARFIELD";

    Fantasybctype[119]  = "ABLATION_SURFACE";

    Fantasybctype[5]    = "INFLOW";
    Fantasybctype[52]   = "PRESSURE_INLET";

    Fantasybctype[6]    = "OUTFLOW";
    Fantasybctype[61]   = "OUTFLOW_CONFINED";
    Fantasybctype[62]   = "PRESSURE_OUTLET";

    Fantasybctype[7]    = "POLE";
    Fantasybctype[71]   = "POLE1";
    Fantasybctype[72]   = "POLE2";
    Fantasybctype[73]   = "POLE3";

    Fantasybctype[8]    = "GENERIC_1";
    Fantasybctype[9]    = "GENERIC_2";
    Fantasybctype[10]   = "GENERIC_3";

    Fantasybctype[11]   = "StartNumberOfBodyOfHyperFLOW";

    Fantasybctype[12]   = "BODY";
    Fantasybctype[13]   = "BODY";
    Fantasybctype[14]   = "BODY";
    Fantasybctype[15]   = "BODY";
    Fantasybctype[16]   = "BODY";
    Fantasybctype[17]   = "BODY";
    Fantasybctype[18]   = "BODY";
    Fantasybctype[19]   = "BODY";
    Fantasybctype[20]   = "BODY";
    Fantasybctype[21]   = "BODY";
    Fantasybctype[22]   = "BODY";
    Fantasybctype[23]   = "BODY";
    Fantasybctype[24]   = "BODY";
    Fantasybctype[25]   = "BODY";
    Fantasybctype[26]   = "BODY";
    Fantasybctype[27]   = "BODY";
    Fantasybctype[28]   = "BODY";
    Fantasybctype[29]   = "BODY";
    Fantasybctype[30]   = "ENDNumberOfBodyOfHyperFLOW";
    Fantasybctype[100]  = "INTERIOR";
    Fantasybctype[1001] = "PRESSURE_INLET_PLATE";
    Fantasybctype[1002] = "PRESSURE_OUTLET_PLATE";
    Fantasybctype[1000] = "OVERSET";

    Fantasybctype[53]   = "MASS_FLOW_INLET";
    Fantasybctype[63]   = "MASS_FLOW_OUTLET";

    //const int EXTERNAL_BC           =   9;
    ////Added by Guo Yongheng 20160815: external overlap cells label

    return Fantasybctype;
}
}
