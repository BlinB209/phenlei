#include "IO_HDF5File.h"
#include "Constants.h"
#include "IO_FileName.h"
#include "TK_Exit.h"

namespace PHSPACE
{
using namespace PHMPI;

hid_t CreateHDF5File(const string &filename)
{
    hid_t file, acc_tpl1;
    herr_t status;

    acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
    status = H5Pset_fapl_mpio(acc_tpl1, MPI_COMM_WORLD, MPI_INFO_NULL);

    file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, acc_tpl1);

    status = H5Pclose(acc_tpl1);
    return file;
}

void ParallelCreateHDF5File(ActionKey *actionKey)
{
    if (actionKey->filename == "")
    {
        actionKey->filepos = 0;
        actionKey->del_file = false;
    }
    else
    {
        if (PHMPI::IsCurrentProcessorFileServer() == 0)
        {
            return;
        }

        int fileID = PHMPI::GetFileIndexofCurrentProcessor();
        string finalFileName = PHSPACE::AddSymbolToFileName(actionKey->filename, '_', fileID);
        actionKey->filepos = H5Fcreate(finalFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        actionKey->del_file = true;
    }
}

void ParallelCloseHDF5File(ActionKey *actionKey)
{
    if (PHMPI::IsCurrentProcessorFileServer() == 0)
    {
        return;
    }

    H5Fclose(actionKey->filepos);
}

hid_t OpenHDF5File(const string &filename)
{
    hid_t file;
    file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    return file;
}

void CreateGroup(hid_t loc_id, const string &groupName)
{
    hid_t grpData;

    grpData = H5Gcreate2(loc_id, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(grpData);
}

hid_t OpenGroup(hid_t loc_id, const string &groupName)
{
    hid_t grpData;

    grpData = H5Gopen2(loc_id, groupName.c_str(), H5P_DEFAULT);
    return grpData;
}

hid_t GetGroup(hid_t loc_id, const string &groupName)
{
    hid_t grpData;

    CreateGroup(loc_id, groupName);
    grpData = OpenGroup(loc_id, groupName);

    return grpData;
}

void CreateEmptyData(hid_t loc_id, const string &dataName, hsize_t rowNumber, hsize_t rowLength, int type)
{
    hid_t dataspace = 0, dataset = 0, datatype = 0;
    int rank = 0;

    //! Dataset dimensions at the creation time.
    hsize_t dims[2] = {0};

    if (rowNumber == 1)
    {
        rank = 1;
        dims[0] = rowLength;
    }
    else
    {
        rank = 2;
        dims[0] = rowNumber;
        dims[1] = rowLength;
    }

    //! Create the data space with unlimited dimensions.
    dataspace = H5Screate_simple(rank, dims, NULL);

    if (type == PHINT)
    {
        datatype = H5Tcopy(H5T_NATIVE_INT);
    }
    else if (type == PHDOUBLE)
    {
        datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    }
    else if (type == PHSTRING)
    {
        datatype = H5Tcopy(H5T_NATIVE_CHAR);
    }

    //! Create a empty dataset.
    dataset = H5Dcreate2(loc_id, dataName.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
}

hid_t OpenData(hid_t loc_id, const string &dataName)
{
    hid_t dataLoc;

    dataLoc = H5Dopen2(loc_id, dataName.c_str(), H5P_DEFAULT);
    if (dataLoc < 0)
    {
        ostringstream oss;
        oss << "The data is not found :  " << dataName;
        TK_Exit::ExceptionExit(oss);
    }

    return dataLoc;
}

void WriteData(hid_t loc_id, const void *data, const string &dataName)
{
    hid_t dataset, datatype, dataspace;
    herr_t status;

    dataset = OpenData(loc_id, dataName);

    datatype = H5Dget_type(dataset);
    dataspace = H5Dget_space(dataset);

    status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);
}

void ReadData(hid_t loc_id, void *data, const string &dataName)
{
    hid_t dataset, datatype, dataspace, memspace;
    herr_t status;
    int rank;

    dataset = OpenData(loc_id, dataName);
    datatype = H5Dget_type(dataset);
    dataspace = H5Dget_space(dataset);
    rank  = H5Sget_simple_extent_ndims(dataspace);

    hsize_t *dataLength = new hsize_t[rank];
    status = H5Sget_simple_extent_dims(dataspace, dataLength, NULL);

    hsize_t *dataLocation = new hsize_t[rank];
    for (int irank = 0; irank < rank; ++ irank)
    {
        dataLocation[irank] = 0;
    }

    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, dataLocation, NULL, dataLength, NULL);

    memspace = H5Screate_simple(rank, dataLength,NULL);
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, dataLocation, NULL, dataLength, NULL);

    status = H5Dread(dataset, datatype, memspace, dataspace, H5P_DEFAULT, data);

    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);

    delete [] dataLength;
    delete [] dataLocation;
}

void CreateAndWriteData(hid_t loc_id, const string &dataName, hsize_t rowNumber, hsize_t rowLength, int type, const void *data)
{
    CreateEmptyData(loc_id, dataName, rowNumber, rowLength, type);
    WriteData(loc_id, data, dataName);
}

RDouble * ReadDoubleOneRow(hid_t loc_id, const string &dataName)
{
    hid_t dataset;
    hid_t datatype, dataspace, memspace;
    herr_t status;
    int rank, status_n;
    ostringstream oss;

    dataset = H5Dopen2(loc_id, dataName.c_str(), H5P_DEFAULT);
    if (dataset < 0)
    {
        oss << "The data is not found :  " << dataName;
        TK_Exit::ExceptionExit(oss);
    }

    datatype = H5Dget_type(dataset);
    dataspace = H5Dget_space(dataset);
    rank = H5Sget_simple_extent_ndims(dataspace);
    if (rank != 1)
    {
        oss << "The data is not One_D :  " << dataName;
        TK_Exit::ExceptionExit(oss);
    }

    hsize_t dataLength;
    status_n = H5Sget_simple_extent_dims(dataspace, &dataLength, NULL);

    hsize_t dataLocation[1] = {0};
    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, dataLocation, NULL, &dataLength, NULL);

    memspace = H5Screate_simple(rank, &dataLength, NULL);
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, dataLocation, NULL, &dataLength, NULL);

    double *Data = new double[dataLength];
    status = H5Dread(dataset, datatype, memspace, dataspace, H5P_DEFAULT, Data);

    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);

    RDouble *dataReturn = new RDouble[dataLength];
    for (int iData = 0; iData < dataLength; ++ iData)
    {
        dataReturn[iData] = static_cast<RDouble>(Data[iData]);
    }
    delete [] Data;    Data = NULL;
    return dataReturn;
}

int * ReadIntOneRow(hid_t loc_id, const string &dataName)
{
    hid_t dataset;
    hid_t datatype, dataspace, memspace;
    herr_t status;
    int rank, status_n;
    ostringstream oss;

    dataset = H5Dopen2(loc_id, dataName.c_str(), H5P_DEFAULT);
    if (dataset < 0)
    {
        oss << "The data is not found :  " << dataName;
        TK_Exit::ExceptionExit(oss);
    }

    datatype = H5Dget_type(dataset);
    dataspace = H5Dget_space(dataset);
    rank = H5Sget_simple_extent_ndims(dataspace);
    if (rank != 1)
    {
        oss << "The data is not One_D :  " << dataName;
        TK_Exit::ExceptionExit(oss);
    }

    hsize_t dataLength;
    status_n = H5Sget_simple_extent_dims(dataspace, &dataLength, NULL);

    hsize_t dataLocation[1] = {0};
    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, dataLocation, NULL, &dataLength, NULL);

    memspace = H5Screate_simple(rank, &dataLength,NULL);
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, dataLocation, NULL, &dataLength, NULL);

    int *Data = new int [dataLength];
    status = H5Dread(dataset, datatype, memspace, dataspace, H5P_DEFAULT, Data);

    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);

    return Data;
}

int ** ReadIntTwoRow(hid_t loc_id, const string &dataName)
{
    hid_t dataset;
    hid_t datatype, dataspace, memspace;
    herr_t status;
    int rank, status_n;
    ostringstream oss;

    dataset = H5Dopen2(loc_id, dataName.c_str(), H5P_DEFAULT);
    if (dataset < 0)
    {
        oss << "The data is not found :  " << dataName;
        TK_Exit::ExceptionExit(oss);
    }

    datatype = H5Dget_type(dataset);
    dataspace = H5Dget_space(dataset);
    rank = H5Sget_simple_extent_ndims(dataspace);

    if (rank != 2)
    {
        oss << "The data is not Two_D :  " << dataName;
        TK_Exit::ExceptionExit(oss);
    }

    hsize_t *dataLength = new hsize_t[rank];
    status_n = H5Sget_simple_extent_dims(dataspace, dataLength, NULL);

    hsize_t *dataLocation = new hsize_t[rank];
    dataLocation[0] = 0;
    dataLocation[1] = 0;

    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, dataLocation, NULL, dataLength, NULL);

    memspace = H5Screate_simple(rank, dataLength, NULL);
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, dataLocation, NULL, dataLength, NULL);

    int **Data = NewPointer2<int>(static_cast<int>(dataLength[0]), static_cast<int>(dataLength[1]));
    status = H5Dread(dataset, datatype, memspace, dataspace, H5P_DEFAULT, Data[0]);

    delete [] dataLocation;
    delete [] dataLength;

    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);

    return Data;
}

char *ReadCharData(hid_t loc_id, const string &dataName)
{
    hid_t dataset;
    hid_t datatype, dataspace, memspace;
    herr_t status;
    int rank, status_n;
    ostringstream oss;

    dataset = H5Dopen2(loc_id, dataName.c_str(), H5P_DEFAULT);
    if (dataset < 0)
    {
        oss << "The data is not found :  " << dataName;
        TK_Exit::ExceptionExit(oss);
    }

    datatype = H5Dget_type(dataset);
    dataspace = H5Dget_space(dataset);
    rank = H5Sget_simple_extent_ndims(dataspace);
    if (rank != 1)
    {
        oss << "The data is not One_D :  " << dataName;
        TK_Exit::ExceptionExit(oss);
    }

    hsize_t dataLength;
    status_n = H5Sget_simple_extent_dims(dataspace, &dataLength, NULL);

    hsize_t dataLocation[1] = {0};
    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, dataLocation, NULL, &dataLength, NULL);

    memspace = H5Screate_simple(rank, &dataLength, NULL);
    status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, dataLocation, NULL, &dataLength, NULL);

    char *Data = new char[dataLength];
    status = H5Dread(dataset, datatype, memspace, dataspace, H5P_DEFAULT, Data);

    H5Sclose(dataspace);
    H5Tclose(datatype);
    H5Dclose(dataset);

    return Data;
}

int GetObjNumber(hid_t loc_id)
{
    H5G_info_t grp_info;
    H5Gget_info(loc_id, &grp_info);

    return static_cast<int>(grp_info.nlinks);
}

void GetObjNameList(hid_t loc_id, string *objNameList)
{
    int objNumber = GetObjNumber(loc_id);
    for (int iObj = 0; iObj < objNumber; ++ iObj)
    {
        char dataname[100];
        int nameSize = 100;
        H5Gget_objname_by_idx(loc_id, iObj, &dataname[0], nameSize);

        objNameList[iObj] = dataname;
    }
}

void GetObjTypeList(hid_t loc_id, H5G_obj_t *objTypeList)
{
    int objNumber = GetObjNumber(loc_id);
    for (int iObj = 0; iObj < objNumber; ++ iObj)
    {
        objTypeList[iObj] = H5Gget_objtype_by_idx(loc_id, iObj);
    }
}

bool CheckDataExist(hid_t loc_id, const string &dataName_in)
{
    bool dataExist = false;

    int objNumber = GetObjNumber(loc_id);
    string *objNameList = new string[objNumber];
    GetObjNameList(loc_id, objNameList);

    for (int iObj = 0; iObj < objNumber; ++ iObj)
    {
        if (objNameList[iObj] == dataName_in)
        {
            dataExist = true;
            break;
        }
    }

    delete [] objNameList;
    return dataExist;
}

}