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
//! @file      IO.h
//! @brief     Declare some IO functions.
//! @author    Refactored by He Xianyao, Bell.

#pragma once
#include "IO_VirtualFile.h"
#include "LIB_Macro.h"
#include "ActionKey.h"

using namespace std;

namespace PHSPACE
{
//! If myid equals fileserver, open a file with "openmode" mode.
//! @param[in] file        file stream flow.
//! @param[in] fileName    file name.
//! @param[in] openMode    open mode.
LIB_EXPORT void ParallelOpenFile(fstream &file, const string &fileName, const ios_base::openmode &openMode);

//! If myid equals fileserver, open a file according to actkey.
//! @param[in] actkey    an action key.
LIB_EXPORT void ParallelOpenFile(ActionKey *actkey);

//! Open a file according to actkey, this function is used by p2p mode.
//! @param[in] actionKey    an action key.
LIB_EXPORT void ParallelOpenFile2(ActionKey *actionKey);

//! Open a file with "openmode" mode.
//! @param[in] file        file stream flow.
//! @param[in] fileName    file name.
//! @param[in] openMode    open mode.
LIB_EXPORT void OpenFile(fstream &file, const string &fileName, const ios_base::openmode &openMode);

//! Open a file with "openmode" mode.
//! It is different with OpenFile, this function will auto add _0 after filename.
//! e.g filename = aaa.fts, it will open aaa_0.fts.
//! @param[in] file        file stream flow.
//! @param[in] fileName    file name.
//! @param[in] openMode    open mode.
LIB_EXPORT void OpenSerialFile(fstream &file, const string &filename, const ios_base::openmode &openmode);

//! If myid equals fileserver, close a file.
//! @param[in] file    file stream flow.
LIB_EXPORT void ParallelCloseFile(fstream &file);

//! If myid equals fileserver, close a file according to actkey.
//! @param[in] actionKey    an action key.
LIB_EXPORT void ParallelCloseFile(ActionKey *actionKey);

//! Close a file.
//! @param[in] file    file stream flow.
LIB_EXPORT void CloseFile(fstream &file);

//! Close a file.
//! @param[in] file    file stream flow.
LIB_EXPORT void CloseFile(ifstream &file);

//! Close a file.
//! @param[in] file    file stream flow.
LIB_EXPORT void CloseFile(ofstream &file);

//! Close a file according to actkey.
//! @param[in] actkey    an action key.
LIB_EXPORT void CloseFile(ActionKey *actkey);

//! Read data from file to datacontainer.
//! @param[in] file     file stream flow.
//! @param[in] cdata    the datacontainer.
LIB_EXPORT void ReadFile(fstream &file, DataContainer *cdata);

//! According to actkey, read data from file to datacontainer.
//! @param[in] actkey    an action key.
LIB_EXPORT void ReadFile(ActionKey *actkey);

//! Write data from datacontainer to file.
//! @param[in] file     file stream flow.
//! @param[in] cdata    the datacontainer.
LIB_EXPORT void WriteFile(fstream &file, DataContainer *cdata);

//! According to actkey, write data from datacontainer to file.
//! @param[in] actkey    an action key.
LIB_EXPORT void WriteFile(ActionKey *actkey);

//! Output str to file.
//! @param[in] file    file stream flow.
//! @param[in] str     a string.
LIB_EXPORT void WriteASCIIFile(fstream &file, const string &str);

//! Output oss to file.
//! @param[in] file    file stream flow.
//! @param[in] oss     a string stream.
LIB_EXPORT void WriteASCIIFile(fstream &file, std::ostringstream &oss);

//! Output datacontainer to file.
//! @param[in] file     file stream flow.
//! @param[in] cdata    the datacontainer.
LIB_EXPORT void WriteASCIIFile(fstream &file, DataContainer *cdata);

//! Output according to actkey.
//! @param[in] actkey    an action key.
//! @param[in] oss       a string stream.
LIB_EXPORT void WriteASCIIFile(ActionKey *actionKey, std::ostringstream &oss);

//! Output according to actkey.
//! @param[in] actkey    an action key.
LIB_EXPORT void WriteASCIIFile(ActionKey *actkey);

//! Useless.
//! @param[in] actkey    an action key.
LIB_EXPORT void WriteBinaryFile(ActionKey *actkey);

//! Write sentinel file needed by GUI.
LIB_EXPORT void WriteSentinelFile();

//! Output to screen from datacontainer.
//! @param[in] cdata    the datacontainer.
LIB_EXPORT void WriteScreen(DataContainer *cdata);

//! Output to screen according to actkey.
//! @param[in] actkey    an action key.
LIB_EXPORT void WriteScreen(ActionKey *actkey);

//! Judge file empty or not.
//! @param[in] file    file stream flow.
LIB_EXPORT bool IfFileEmpty(fstream &file);

//! Judge file exist or not.
//! @param[in] fileName    file name.
LIB_EXPORT bool FileExist(const string &fileName);

//! Judge SerialFile exist or not.
//! SerialFile names add _0 after filename, e.g filename = aaa.fts.
//! This will judge if aaa_0.fts exist or not.
//! @param[in] fileName    file name.
LIB_EXPORT bool SerialFileExist(const string &fileName);

//! Read from file and send the data to others.
//! @param[in] file         file stream flow.
//! @param[in] cdata        the datacontainer.
//! @param[in] send_proc    send process.
//! @param[in] recv_proc    receive process.
//! @param[in] tag          tag flag.
LIB_EXPORT void ReadAbstractData(fstream &file, DataContainer *cdata, int send_proc, int recv_proc, int tag = 0);

//! Read control parameters info.
//! At first, read keyInfo from file(./bin/key.hypara).
//! Then, according to the nsimutask info from key.hypara, read corresponding file.
//! Example:    int nsimutask = 2;
//!             string parafilename = "./bin/cfd_para.hypara";
LIB_EXPORT void ReadControlParameters();

LIB_EXPORT void ReadKeyInfo();

LIB_EXPORT void ReadDefaultParamInfo();

LIB_EXPORT void ReadParamInfo();

LIB_EXPORT void ReadParamInfo(const string &filename);

LIB_EXPORT void BroadcastInfo();

LIB_EXPORT void CompleteParamInfo();

LIB_EXPORT void PrintInflowConditionToWindow(void);

LIB_EXPORT void ModifyWindowTemperature(RDouble &wallTemperature, vector< RDouble > &WallT, vector< int > &WallNumber);

#include "PHIO.hxx"
}