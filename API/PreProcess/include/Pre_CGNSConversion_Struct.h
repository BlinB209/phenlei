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
//! @file      Pre_CGNSConversion_Struct.h
//! @brief     Grid conversion from CGNS.
//! @author    Xu Gang.

#pragma once
#include "LIB_Macro.h"
#include "Pre_GridConversion.h"
#include "cgnslib.h"

using namespace std;
namespace PHSPACE
{
class CGNS_Str_Data;

class Pre_CGNSConversion_Struct : public Pre_GridConversion
{
public:
    LIB_EXPORT Pre_CGNSConversion_Struct(const string &gridFileName);
    LIB_EXPORT ~Pre_CGNSConversion_Struct();

private:
    void ReadGrid();

private:
    void tram(int *a, int *b, int *v, int nCoords);
    int del(int x, int y);

public:
    void ResetGridScaleAndTranslate();

public:
    CGNS_Str_Data *Str_CGNS_Data;

};

class BaseData;
typedef char char_33[33];
typedef cgsize_t cgsize_t_60000[60000];

//! @brief Collection the grid information from CGNS.
//! CGNS_Str_Data
//! BaseData
class CGNS_Str_Data
{
public:
    CGNS_Str_Data(int nzones_in);
    ~CGNS_Str_Data(void);

private:
    int nzones;
    BaseData **base_cgns;

public:
    BaseData * GetBaseData(int iZone) { return base_cgns[iZone]; }
};

//! @brief Collection the grid information in each block.
//! BaseData
class BaseData
{
public:
    BaseData(void);
    ~BaseData(void);

private:
    int nCoords, idim, jdim, kdim;
    int nconnect, nBCRegions, n1to1;
    RDouble ***x, ***y, ***z;
    char_33 *zonenames, *familyName;
    BCType_t *bocoType;
    cgsize_t_60000 *pnts;
    string *bcName;
    string  volumeName;
    int vcType;
    int *lab1, *lab2;
    int ier;
    cgsize_t **srange, **donor_range;
    int **b, **vt;
    char_33 *connectname, *donorname;

public:
    void SetIdim(int data_in) { this->idim = data_in; }
    void SetJdim(int data_in) { this->jdim = data_in; }
    void SetKdim(int data_in) { this->kdim = data_in; }

    int GetIdim() { return idim; }
    int GetJdim() { return jdim; }
    int GetKdim() { return kdim; }

    void SetNCoords(int nCoords_in) { this->nCoords = nCoords_in; }
    int GetNCoords() { return nCoords; }

    void SetX(RDouble ***x_in) { this->x = x_in; }
    void SetY(RDouble ***y_in) { this->y = y_in; }
    void SetZ(RDouble ***z_in) { this->z = z_in; }
    RDouble *** GetX() { return x; }
    RDouble *** GetY() { return y; }
    RDouble *** GetZ() { return z; }

    void SetZonenames(char_33 *zonenames_in) { this->zonenames = zonenames_in; }
    char_33 * GetZonenames() { return zonenames; }

    void SetFamilyName(char_33 *familyName_in) { this->familyName = familyName_in; }
    char_33 * GetFamilyName() { return familyName; }

    void SetBocoType(BCType_t *bocoType_in) { this->bocoType = bocoType_in; }
    BCType_t * GetBocoType() { return bocoType; }

    void SetPnts(cgsize_t_60000 *pnts_in) { this->pnts = pnts_in; }
    cgsize_t_60000 * GetPnts() { return pnts; }

    void SetBcName(string *bcName_in) { this->bcName = bcName_in; }
    string * GetBcName() { return bcName; }

    void SetLab1(int *lab1_in) { this->lab1 = lab1_in; }
    int * GetLab1() { return lab1; }

    void SetNconnect(int nconnect_in) { this->nconnect = nconnect_in; }
    int GetNconnect() { return nconnect; }

    void SetNBCRegions(int nBCRegions_in) { this->nBCRegions = nBCRegions_in; }
    int GetNBCRegions() { return nBCRegions; }

    string GetVCName() { return volumeName; }
    void SetVCName(string nameIn) { volumeName = nameIn; }

    int  GetVCType() { return vcType; }
    void SetVCType(int vcTypeIn) { vcType = vcTypeIn; }

    void SetN1to1(int n1to1_in) { this->n1to1 = n1to1_in; }
    int GetN1to1() { return n1to1; }

    void SetB(int *b_in, int one21_in) { this->b[one21_in] = b_in; }
    int * GetB(int one21_in) { return b[one21_in]; }

    void SetVt(int *vt_in, int one21_in) { this->vt[one21_in] = vt_in; }
    int * GetVt(int one21_in) { return vt[one21_in]; }

    void ResizeRange(int n1to1_in);
    void SetSrange(cgsize_t *srange_in, int one21_in) { this->srange[one21_in] = srange_in; }
    cgsize_t * GetSrange(int one21_in) { return srange[one21_in]; }

    void SetDonor_range(cgsize_t *donor_range_in, int one21_in) { this->donor_range[one21_in] = donor_range_in; }
    cgsize_t * GetDonor_range(int one21_in) { return donor_range[one21_in]; }

    void SetLab2(int lab2_in, int one21_in) { this->lab2[one21_in] = lab2_in; }
    int Getlab2(int one21_in) { return lab2[one21_in]; }

    void SetIer(int ier_in) { this->ier = ier_in; }
    int GetIer() { return ier; }

    void SetConnectname(char_33 *connectname_in) { this->connectname = connectname_in; }
    void SetDonorname(char_33 *donorname_in) { this->donorname = donorname_in; }
    char_33 * GetConnectname() { return connectname; }
    char_33 * GetDonorname() { return donorname; }
};






}


