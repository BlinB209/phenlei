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
//! @file      Force.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include "DataContainer.h"
using namespace std;

namespace PHSPACE
{

class ActionKey;
class DataContainer;
class AerodynamicForceManager;
class AerodynamicForce;

class Force
{
public:
    static const int ALL_SOLID_WALL = -1;

    //! Number of force variables.
    static const int nForceVar = 16;

private:
    //! Number of solid wall part.
    int numberofSolidWallPart;

    //! Number of total part, is numberofSolidWallPart plus ONE.
    //! 'ONE' represent the total solid wall force, stored on the first position.
    int numberofTotalPart;

    //! The corresponding Boundary Condition ID of each part.
    int *partBCID;

    RDouble *cpx, *cpy, *cpz;
    RDouble *cl_tot, *cd_tot, *cd_pr, *cd_sf, *cd_cl2pa, *xcp, *ycp, *zcp, *side;
    RDouble *CA_f, *CA_p, *CN_f, *CN_p, *CZ_f, *CZ_p, *Cl_f, *Cl_p, *Cn_f, *Cn_p, *Cm_f, *Cm_p;

    int numberOfBody;
    vector < vector <string> > wallNameOfEachBody;
    RDouble *CA_f_body, *CN_f_body, *CZ_f_body, *Cl_f_body, *Cn_f_body, *Cm_f_body;
    RDouble *CA_p_body, *CN_p_body, *CZ_p_body, *Cl_p_body, *Cn_p_body, *Cm_p_body;

    RDouble *hingeMoment; 
    RDouble attack, sideslip;
    RDouble forceReferenceArea, forceReferenceLength, forceReferenceLengthSpanWise;

public:
    Force();
    ~Force();

public:
    void Init(RDouble attack, RDouble sideslip, RDouble forceReferenceArea, RDouble forceReferenceLength, RDouble forceReferenceLengthSpanWise);
    void CollectionForce(DataContainer *data_in);
    void CompForce();
    void DumpForce(int outnstep,std::ostringstream &oss);
    void DumpPartForce(int outnstep, int partID, std::ostringstream &oss);
    void DumpForce(ActionKey *actkey, int outnstep);
    void DumpData(DataContainer *data);

    //! Dump all force data in forcefile.
    void DumpAllForceData(ActionKey *actkey, int outnstep);

    //! Dump global force data.
    void DumpAllGlobalForceData(const string &partName, int partID, std::ostringstream &oss);

    //! Dump local force data.
    void DumpAllLocalForceData(const string &partName, int partID, std::ostringstream &oss);

    void DumpBodyForceData(int bodyID, std::ostringstream &oss);
private:
    bool isPartExist;
};

void InitForce();
void CollectionForce(ActionKey *actkey);
void DumpForce(ActionKey *actkey);

class AerodynamicForceComponent
{
public:
    AerodynamicForceComponent();
    ~AerodynamicForceComponent();
protected:
    static vector< std::set< string > > components;
    static bool initialize;
public:
    static void Initialize();
    static bool FaceInAerodynamicForceComponents(string componentName, int aerodynamicForceComponentIndex);
};

bool FaceInAerodynamicForceComponents(string componentName, int aerodynamicForceComponentIndex);

}
