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
//! @file      IO_HDF5File.h
//! @brief     file reader.
//! @author    Baka.

#pragma once
#include "hdf5.h"
#include "Precision.h"
#include "Pointer.h"
#include "ActionKey.h"

using namespace std;

namespace PHSPACE
{
//! Create a HDF5 file.
//! @param[in] filename    file name.
hid_t CreateHDF5File(const string &filename);

//! Open a HDF5 file.
//! @param[in] filename    file name.
hid_t OpenHDF5File(const string &filename);

//! If myid equals fileserver, create a HDF5 file according to actkey.
//! @param[in] actkey    an action key.
void ParallelCreateHDF5File(ActionKey *actionKey);

//! If myid equals fileserver, close a HDF5 file according to actkey.
//! @param[in] actkey    an action key.
void ParallelCloseHDF5File(ActionKey *actionKey);

//!
void CreateGroup(hid_t loc_id, const string &groupName);

hid_t OpenGroup(hid_t loc_id, const string &groupName);

hid_t GetGroup(hid_t loc_id, const string &groupName);

//!
void CreateEmptyData(hid_t loc_id, const string &dataName, hsize_t rowNumber, hsize_t rowLength, int type);

hid_t OpenData(hid_t loc_id, const string &dataName);

//!
void CreateAndWriteData(hid_t loc_id, const string &dataName, hsize_t rowNumber, hsize_t rowLength, int type, const void *data);

void WriteData(hid_t loc_id, const void *data, const string &dataName);

void ReadData(hid_t loc_id, void *data, const string &dataName);

//!
RDouble * ReadDoubleOneRow(hid_t loc_id, const string &dataName);

int  * ReadIntOneRow(hid_t loc_id, const string &dataName);

int ** ReadIntTwoRow(hid_t loc_id, const string &dataName);

char * ReadCharData(hid_t loc_id, const string &dataName);

bool CheckDataExist(hid_t loc_id, const string &dataName);

int  GetObjNumber(hid_t loc_id);

void GetObjNameList(hid_t loc_id, string *objNameList);

void GetObjTypeList(hid_t loc_id, H5G_obj_t *objTypeList);
}