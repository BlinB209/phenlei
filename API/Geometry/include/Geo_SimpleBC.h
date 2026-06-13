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
//! @file      Geo_SimpleBC.h
//! @brief     It is the base class of geometry operation for boundary condition.
//!            The inheriting order is: SimpleBC -> StructBC/UnstructBC.
//! @author    Zhang Yong, Bell, He Xin.

#pragma once
#include "TypeDefine.h"
#include "Constants.h"
#include "Data_Param.h"
#include "Data_Field.h"
using namespace std;

namespace PHSPACE
{
namespace PHENGLEI
{
    //! the 1000 for periodic is default output of the User Defined for CGNS.
    const int PERIODIC                     = 1000;
    const int USER_DEFINED                 = -3;
    const int WAKE                         = -2;
    const int INTERFACE                    = -1;
    const int NO_BOUNDARY_CONDITION        = 0;
    const int EXTRAPOLATION                = 1;
    const int SOLID_SURFACE                = 2;
    const int SYMMETRY                     = 3;
    const int FARFIELD                     = 4;
    const int INFLOW                       = 5;
    const int OUTFLOW                      = 6;

    //! surface with ablation is 119.
    const int ABLATION_SURFACE             = 119;

    //! pressure in and out boundary 52, 62.
    const int PRESSURE_INLET               = 52;
    const int PRESSURE_OUTLET              = 62;

    //! pressure in and out for nozzle boundary.
    const int NOZZLE_INFLOW                = 53;
    const int NOZZLE_OUTFLOW               = 63;

    //! mass flow in and out boundary 53, 63.
    const int MASS_FLOW_INLET              = 53;
    const int MASS_FLOW_OUTLET             = 63;

    const int OUTFLOW_CONFINED             = 61;
    const int POLE                         = 7;
    const int POLE1                        = 71;
    const int POLE2                        = 72;
    const int POLE3                        = 73;
    const int GENERIC_1                    = 8;
    const int GENERIC_2                    = 9;
    const int GENERIC_3                    = 10;
    const int INTERIOR                     = 100;

    const int EXTERNAL_BC                  = 9;    //! External overlap cells label.

    //! Boundary condition for special case.
    const int PRESSURE_INLET_PLATE         = 1001;
    const int PRESSURE_OUTLET_PLATE        = 1002;

    const int StartNumberOfBodyOfHyperFLOW = 11;
    const int ENDNumberOfBodyOfHyperFLOW   = 30;
    const int BODY1                        = 11;
    const int BODY2                        = 12;
    const int BODY3                        = 13;
    const int BODY4                        = 14;
    const int BODY5                        = 15;
    const int BODY6                        = 16;
    const int BODY7                        = 17;
    const int BODY8                        = 18;
    const int BODY9                        = 19;
    const int BODY10                       = 20;
    const int BODY11                       = 21;
    const int BODY12                       = 22;
    const int BODY13                       = 23;
    const int BODY14                       = 24;
    const int BODY15                       = 25;
    const int BODY16                       = 26;
    const int BODY17                       = 27;
    const int BODY18                       = 28;
    const int BODY19                       = 29;
    const int BODY20                       = 30;
    const int OVERSET                      = 1000;
}

namespace INIFLOW_SPACE
{
    const int RECTANGLE             = 0;
    const int SPHERE                = 1;
    const int CYLINDER              = 2;
    const int CONE                  = 3;

    const int CUSTOM_ANGLE          = 0;
    const int CUSTOM_VECTOR         = 1;
}



namespace GRIDGEN_SPACE
{
    const int INTERFACE             = -1;
    const int NO_BOUNDARY_CONDITION = 0;
    const int SOLID_SURFACE         = 2;
    const int SYMMETRY              = 3;
    const int FARFIELD              = 4;
    const int INFLOW                = 5;
    const int OUTFLOW               = 6;
    const int POLE                  = 7;
    const int GENERIC_1             = 8;
    const int GENERIC_2             = 9;
    const int GENERIC_3             = 10;
    const int PERIODIC              = 99;
}

namespace FLUENT_SPACE
{
    const int FLUENT_COMMENT            = 0;
    const int FLUENT_HEADER             = 1;
    const int FLUENT_DIMENSIONS         = 2;
    const int FLUENT_NODES              = 10;
    const int FLUENT_NODES_BINARY  = 3010;
    const int FLUENT_EDGES              = 11;
    const int FLUENT_CELLS              = 12;
    const int FLUENT_CELLS_BINARY  = 3012;
    const int FLUENT_FACES              = 13;
    const int FLUENT_FACES_BINARY  = 3013;
    const int FLUENT_PERIODICFACES      = 18;
    const int FLUENT_C39                = 39;
    const int FLUENT_C45                = 45;
    const int MaxFluentBCType           = 38;
    const int FLUENT_CELLTREE           = 58;
    const int FLUENT_FACETREE           = 59;

    const int NODETYPE_GHOST            = 0;
    const int NODETYPE_NOTYPE           = 1;
    const int NODETYPE_BOUNDARY         = 2;

    const int ELEMENTTYPE_MIXED         = 0;
    const int ELEMENTTYPE_TRIANGULAR    = 1;
    const int ELEMENTTYPE_TETRAHEDRAL   = 2;
    const int ELEMENTTYPE_QUADRILATERAL = 3;
    const int ELEMENTTYPE_HEXAHEDRAL    = 4;
    const int ELEMENTTYPE_PYRAMID       = 5;
    const int ELEMENTTYPE_WEDGE         = 6;
    const int ELEMENTTYPE_POLYHEDRAL    = 7;

    const int CELLTYPE_DEAD             = 0;
    const int CELLTYPE_FLUID            = 1;
    const int CELLTYPE_SOLID            = 17;

    const int FACETYPE_MIXED            = 0;
    const int FACETYPE_LINEAR           = 2;
    const int FACETYPE_TRIANGULAR       = 3;
    const int FACETYPE_QUADRILATERAL    = 4;
    const int FACETYPE_POLYGONAL        = 5;

    const int UNSPECIFIED               = 0;
    const int INTERIOR                  = 2;
    const int WALL                      = 3;
    const int PRESSURE_INLET            = 4;
    const int INLET_VENT                = 4;
    const int INTAKE_FAN                = 4;
    const int PRESSURE_OUTLET           = 5;
    const int EXHAUST_FAN               = 5;
    const int OUTLET_VENT               = 5;
    const int SYMMETRY                  = 7;
    const int PERIODIC_SHADOW           = 8;
    const int PRESSURE_FAR_FIELD        = 9;
    const int VELOCITY_INLET            = 10;
    const int PERIODIC                  = 12;
    const int FAN                       = 14;
    const int POROUS_JUMP               = 14;
    const int RADIATOR                  = 14;
    const int MASS_FLOW_INLET           = 20;
    const int INTERFACE                 = 24;
    const int PARENT                    = 31;
    const int OUTFLOW                   = 36;
    const int AXIS                      = 37;
}

//! @brief A boundary condition base class, StructBC and UnstructBC inheriting from it.
class SimpleBC
{
public:
    SimpleBC()
    {
        bcName = "";
        bodyName = "";
        bcType = 0;

        bcParamDataBase = new Data_Param();
        bcFieldDataBase = new Data_Field();
    }

    ~SimpleBC()
    {
        delete bcParamDataBase;    bcParamDataBase = NULL;
        delete bcFieldDataBase;    bcFieldDataBase = NULL;
    }

private:
    //! Boundary name of BC;
    string bcName;

    string bodyName;

    //! Boundary type value of BC;
    int bcType;

    //! BC control parameters.
    Data_Param *bcParamDataBase;

    //! BC data pointers.
    Data_Field *bcFieldDataBase;

public:
    const string & GetBCName() const { return this->bcName; }

    const string & GetBodyName()
    { 
        //string *bodyName = new string;
        if (bcParamDataBase->IsExist("bodyName", PHSTRING, 1))
        {
            bcParamDataBase->GetData("bodyName", &bodyName, PHSTRING, 1);
        }
        return bodyName; 
    }

    //! Get boundary conditions type.
    int GetBCType() const { return this->bcType; }

    Data_Param * GetBCParamDataBase();
    Data_Field * GetBCFieldDataBase();

    void SetBCName(const string &bcName) { this->bcName = bcName; }
    void SetBodyName(const string &bodyName) { this->bodyName = bodyName; }
    void SetBCType(int bcType) { this->bcType = bcType; }
    void SetBCParamDataBase(Data_Param *paraDB)  { this->bcParamDataBase = paraDB; }
    void SetBCFieldDataBase(Data_Field *fieldDB) { this->bcFieldDataBase = fieldDB; }

    void UpdateParamData(const string &name, void *data, int type, int size) 
    {
        bcParamDataBase->UpdateData(name, data, type, size);
    }

    void GetParamData(const string &name, void *data, int type, int size) 
    {
        bcParamDataBase->GetData(name, data, type, size);
    }

    int CheckParamData(const string &name)
    {
        return bcParamDataBase->CheckDataExist(name);
    }

    void UpdateFieldDataPtr(const string &name, void *data)
    {
        bcFieldDataBase->UpdateDataPtr(name, data);
    };

    void * GetFieldDataPtr(const string &name) const
    {
        return bcFieldDataBase->GetDataPtr(name);
    };

    int CheckFieldData(const string &name)
    {
        return bcFieldDataBase->CheckDataExist(name);
    }
};

//! @brief Global boundary conditions, all of the global boundary conditions are stored here.
class GlobalBoundaryCondition
{
private:
    //! Warning: it has not been deleted!
    static vector<SimpleBC *> *globalBCList;
    static map <int, int> *bodyMap;
    static set <string> *bodyList;

private:
    GlobalBoundaryCondition();
    ~GlobalBoundaryCondition();

public:
    //! Get number of total BC.
    static uint_t GetNumberOfBC()
    {
        if (!globalBCList)
            return 0;
        else
            return globalBCList->size();
    }

    static uint_t GetNumberOfBody()
    {
        if (!bodyList)
        {
            return 0;
        }
        else
        {
            return bodyList->size();
        }
    }

    static uint_t GetBodyIndex(int iBC)
    {
        return (*bodyMap)[iBC];
    }

    static string GetBoundaryName(int iBC)
    {
        return (*globalBCList)[iBC]->GetBCName();
    }

    static set <string> * GetBodyNameList() { return bodyList; }

    //! Get number of BC with solid_wall BC type.
    static int GetNumberOfSolidWallPart();

    //! Is solid wall.
    static bool IsSolidWallBC(const int &partID)
    {
        return ((*globalBCList)[partID]->GetBCType() == PHENGLEI::SOLID_SURFACE /*|| PHENGLEI::ABLATION_SURFACE*/);
    }

    static SimpleBC *GetBC(const int &iBC) { return (*globalBCList)[iBC]; }

    //! Get global boundary conditions list.
    static vector<SimpleBC *> *GetGlobalBoundaryConditionList() { return globalBCList; }

    //! Read global boundary conditions.
    static void ReadGlobalBoundaryCondition();

    static void SetGlobalBCByGlobalDataBase();

    static void SetBCDataBaseByGlobalBC();

    static void SetDefaultSolidBoundaryCondition();

    static void ChangeBCTypeByGlobalBC();

    static void SetBCInitialValuesByReynolds(Data_Param *bcParamDB, RDouble *primitiveVar);

    static void SetBCInitialValuesByHeight(Data_Param *bcParamDB, RDouble *primitiveVar);

    static void SetMaximumSpeciesForBC(Data_Param *bcParamDB);

    static void SetSpeedDrictionForBC(Data_Param *bcParamDB, RDouble localNonDimensionalVelocity, RDouble *primitiveVar);

    static void SetSpeedDrictionForReferenceVar(RDouble localNonDimensionalVelocity, RDouble *primitiveVar);

    static void SetBCInitialValuesByDensity(Data_Param *bcParamDB, RDouble *primitiveVar);

    static void SetBCInitialValuesByPressure(Data_Param *bcParamDB, RDouble *primitiveVar);

    static void SetBCInitialValuesByTotalPressure(Data_Param *bcParamDB, RDouble *primitiveVar);

    static void SetBCInitialValuesByMachNumberTemperaturePressure(Data_Param *bcParamDB, RDouble *primitiveVar);

    static void SetBCInitialValuesByWRFData(Data_Param *bcParamDB, RDouble *primitiveVar);

    static void SetBCInitialValuesByPrimitive(Data_Param *bcParamDB, RDouble *primitiveVar);

    static void SetBCInitialValuesByMolecularWeightGama(Data_Param *bcParamDB, RDouble *primitiveVar);

    static void InitMassFlowBoundary();
private:
    static void ParseBCFromFile(fstream &file);
    static void BuildBodyMap();
};

//! Judge whether a boundary condition type number belongs to interface.
bool IsInterface(const int &boundaryConditionType);

//! Judge whether a boundary condition type number belongs to solid surface.
inline bool IsWall(const int &bctype)
{
    return (bctype == PHENGLEI::SOLID_SURFACE);
}

map<int, string> CreatFantasybctype();

}