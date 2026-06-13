#include "Geo_SimpleVC.h"
#include "PHIO.h"
#include "Constants.h"
#include "TK_Parse.h"
#include "TK_Exit.h"
#include "Geo_GridIndex.h"
#include "Geo_Grid.h"

namespace PHSPACE
{
vector<SimpleVC *> * GlobalVolumeCondition::globalVCList = 0;

Data_Param *SimpleVC::GetVCParamDataBase()
{
    return vcParamDataBase;
}

void SimpleVC::SetVCParamDataBase(Data_Param *paraDB)
{
    if (vcParamDataBase)
    {
        delete vcParamDataBase;
    }
    this->vcParamDataBase = paraDB;
}


GlobalVolumeCondition::GlobalVolumeCondition()
{

}

GlobalVolumeCondition::~GlobalVolumeCondition()
{

}

void GlobalVolumeCondition::ReadGlobalVolumeCondition()
{
    fstream file;
    string bcFile = "./bin/volume_condition.hypara";

    if (!FileExist(bcFile))
    {
        //! To compatible the old version.
        //SetDefaultSolidBoundaryCondition();
        return;
    }

    OpenFile(file, bcFile, ios::in);

    ParseVCFromFile(file);

    CloseFile(file);
}

void GlobalVolumeCondition::ParseVCFromFile(fstream &file)
{
    string line, keyword, name, word;
    string *value;
    string separator = " =\r\n\t#$,;\"";
    string sep1 = "\r\n";
    string errMsg = "error in volume condition file";

    map < string, int > keywordmap;

    keywordmap.insert(pair<string, int>("int"   , PHINT  ));
    keywordmap.insert(pair<string, int>("float" , PHFLOAT));
    keywordmap.insert(pair<string, int>("double", PHDOUBLE));
    keywordmap.insert(pair<string, int>("string", PHSTRING));

    int nVolumeConditions = 0;
    while (!file.eof())
    {
        getline(file, line);
        if (line == "") continue;
        FindNextWord(line, word, sep1);
        if (word.substr(0, 1) == "#" || word.substr(0, 2) == "//") continue;

        line = FindNextWord(line, keyword, separator);
        if (keyword == "") continue;

        line = FindNextWord(line, name, separator);

        value = new string[1]();
        if (name == "vcName")
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

        int type, size = 1;
        type = keywordmap[keyword];
        if (type == PHINT && name != "nVolumeConditions")
        {
            TK_Exit::ExceptionExit("Error: the first line must be format \"int nVolumeConditions = INTEGER\"!", true);
        }
        if (type == PHSTRING && name != "vcName")
        {
            TK_Exit::ExceptionExit("Error: VC name must be format string vcName = \" FLUID\" !", true);
        }
        if (type == PHINT)
        {
            int *data = new int [size];
            for (int i = 0 ; i < size; ++ i)
            {
                from_string< int >(data[i], value[i], std::dec);
            }

            nVolumeConditions = data[0];
            delete [] data;

            if (nVolumeConditions == 0)
            {
                break;
            }
            else
            {
                globalVCList = new vector <SimpleVC *>;
            }
        }
        else if (type == PHSTRING)
        {
            //! Parse out the parameters into data base in each Boundary Condition.
            SimpleVC *volumeCondition = new SimpleVC();
            volumeCondition->SetVCName(*value);

            Data_Param *vcParamDataBase = new Data_Param();
            TK_Parse::ReadBasicData(file, vcParamDataBase);
            volumeCondition->SetVCParamDataBase(vcParamDataBase);

            int vcType;
            vcParamDataBase->GetData("vcType", &vcType, PHINT, 1);
            volumeCondition->SetVCType(vcType);

            globalVCList->push_back(volumeCondition);
        }

        delete [] value;    value = nullptr;
    }
}

void GlobalVolumeCondition::SetVCDataByGlobalVC()
{
    vector <SimpleVC *> *globalVCList = GlobalVolumeCondition::GetGlobalVolumeConditionList();
    if (globalVCList == 0) return;

    string errInfo = " Error: The volume_condition.hypara does not match the grid file!!!\n";

    using namespace PHMPI;
    int nZones = PHMPI::GetNumberofGlobalZones();
    int myid   = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int procID = PHMPI::GetZoneProcessorID(iZone);

        if (myid == procID)
        {
            Grid *gridofThisProcessor = GetGrid(iZone, 0);
            SimpleVC *volumeCondition = gridofThisProcessor->GetVolumeConditionIn();
            if (!volumeCondition)
            {
                TK_Exit::ExceptionExit(errInfo);
            }

            string vcName = volumeCondition->GetVCName();

            vector <SimpleVC *> :: iterator iter;
            bool matchConditionFind = false;
            for (iter = globalVCList->begin(); iter != globalVCList->end(); iter ++)
            {
                string globalVCName = (*iter)->GetVCName();
                if (globalVCName == vcName)
                {
                    Data_Param *globalVCParam = (*iter)->GetVCParamDataBase();
                    Data_Param *vcParamDataBase = new Data_Param(*globalVCParam);
                    volumeCondition->SetVCParamDataBase(vcParamDataBase);

                    int globalVCType = (*iter)->GetVCType();
                    volumeCondition->SetVCType(globalVCType);

                    matchConditionFind = true;
                    break;
                }
            }

            if (!matchConditionFind)
            {
                TK_Exit::ExceptionExit(errInfo);
            }
        }
    }
}


}