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
//! @file      Geo_SimpleVC.h
//! @brief     It is the base class of geometry operation for volume condition.
//! @author    Baka.

#pragma once
#include "TypeDefine.h"
#include "Data_Param.h"

using namespace std;

namespace PHSPACE
{
namespace PHENGLEI
{
    const int UNDEFINED = 0;
    const int FLUID     = 1;
    const int SOLID     = 2;
}

//! @brief A volume condition base class.
class SimpleVC
{
public:
    SimpleVC()
    {
        vcName = "";
        vcType = 0;

        vcParamDataBase = new Data_Param();
    }

    ~SimpleVC()
    {
        delete vcParamDataBase;
        vcParamDataBase = NULL;
    }

private:
    //! volume condition name;
    string vcName;

    //! volume condition type;
    int vcType;

    //! VC control parameters.
    Data_Param *vcParamDataBase;

public:
    const string &GetVCName() const { return this->vcName; }

    //! Get volume conditions type.
    int GetVCType() const { return this->vcType; }

    Data_Param *GetVCParamDataBase();

    void SetVCName(const string &vcNameIn) { this->vcName = vcNameIn; }
    void SetVCType(int vcTypeIn) { this->vcType = vcTypeIn; }
    void SetVCParamDataBase(Data_Param *paraDB);

    void UpdateParamData(const string &name, void *data, int type, int size) 
    {
        vcParamDataBase->UpdateData(name, data, type, size);
    }

    void GetParamData(const string &name, void *data, int type, int size) 
    {
        vcParamDataBase->GetData(name, data, type, size);
    }

    int CheckParamData(const string &name)
    {
        return vcParamDataBase->CheckDataExist(name);
    }
};

//! @brief Global volume conditions, all of the global volume conditions are stored here.
class GlobalVolumeCondition
{
private:
    //! Warning: it has not been deleted!
    static vector<SimpleVC *> *globalVCList;

private:
    GlobalVolumeCondition();
    ~GlobalVolumeCondition();

public:
    //! Get number of total VC.
    static uint_t GetNumberOfVC()
    {
        if (!globalVCList)
            return 0;
        else
            return globalVCList->size();
    }

    static string GetVolumeName(int iVC)
    {
        return (*globalVCList)[iVC]->GetVCName();
    }

    static SimpleVC *GetVC(const int &iVC) { return (*globalVCList)[iVC]; }

    //! Get global volume conditions list.
    static vector<SimpleVC *> *GetGlobalVolumeConditionList() { return globalVCList; }

    //! Read global volume conditions.
    static void ReadGlobalVolumeCondition();

    static void SetVCDataByGlobalVC();

private:
    static void ParseVCFromFile(fstream &file);

};

}