#include "PHIO.h"
#include "IO_FileName.h"
#include "TK_Exit.h"
#include "TK_Parse.h"
#include "GlobalDataBase.h"
#include "Constants.h"
#include "TK_Log.h"
#include "Gas.h"
#include "Geo_SimpleBC.h"

using namespace std;

namespace PHSPACE
{

LIB_EXPORT void ParallelOpenFile(fstream &file, const string &fileName, const ios_base::openmode &openMode)
{
    if (PHMPI::IsCurrentProcessorFileServer() == 0)
    {
        return;
    }

    PHSPACE::OpenFile(file, fileName, openMode);
}

LIB_EXPORT void OpenSerialFile(fstream &file, const string &filename, const ios_base::openmode &openmode)
{
    int fileID = PHMPI::GetFileIndexofCurrentProcessor();
    string actualFileName = AddSymbolToFileName(filename, "_", fileID);
    OpenFile(file, actualFileName, openmode);
}

LIB_EXPORT void ParallelOpenFile(ActionKey *actkey)
{
    if (!actkey->IsNeedOpenFile())
    {
        fstream *file = 0;
        actkey->file = file;
        actkey->del_file = false;
    }
    else
    {
        if (PHMPI::IsCurrentProcessorFileServer() == 0)
        {
            return;
        }

        fstream *file = new fstream();
        //! Open file by default.
        ParallelOpenFile(*file, actkey->filename, actkey->openmode);
        actkey->file = file;
        actkey->del_file = true;
    }
}

LIB_EXPORT void ParallelOpenFile2(ActionKey *actionKey)
{
    //! Open file by default.
    if (actionKey->filename == "")
    {
        fstream *file = 0;
        actionKey->file = file;
        actionKey->del_file = false;
    }
    else
    {
        if (PHMPI::IsCurrentProcessorFileServer() == 0)
        {
            return;
        }

        fstream *file = new fstream();
        int fileID = PHMPI::GetFileIndexofCurrentProcessor();
        string finalFileName = PHSPACE::AddSymbolToFileName(actionKey->filename, '_', fileID);
        PHSPACE::ParallelOpenFile(*file, finalFileName, actionKey->openmode);
        actionKey->file = file;
        actionKey->del_file = true;
    }
}

LIB_EXPORT void OpenFile(fstream &file, const string &fileName, const ios_base::openmode &openMode)
{
    file.open(fileName.c_str(), openMode);
    if (!file)
    {
        ostringstream oss;
        oss << "Could not open " << fileName << endl;

        PHMPI::FreeBlockData();
        cout << oss.str() << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}

LIB_EXPORT void ParallelCloseFile(fstream &file)
{
    if (PHMPI::IsCurrentProcessorFileServer() == 0)
    {
        return;
    }

    PHSPACE::CloseFile(file);
}

LIB_EXPORT void ParallelCloseFile(ActionKey *actionKey)
{
    if (PHMPI::IsCurrentProcessorFileServer() == 0)
    {
        return;
    }
    CloseFile(actionKey);
}

LIB_EXPORT void CloseFile(fstream &file)
{
    file.close();
    file.clear();
}

LIB_EXPORT void CloseFile(ifstream &file)
{
    file.close();
    file.clear();
}

LIB_EXPORT void CloseFile(ofstream &file)
{
    file.close();
    file.clear();
}

LIB_EXPORT void CloseFile(ActionKey *actkey)
{
    if (!actkey->file)
    {  
        return;
    }
    CloseFile(* actkey->file);

    if (actkey->del_file)
    {
        delete actkey->file;
        actkey->del_file = false;
    }
}

LIB_EXPORT void ReadFile(fstream &file, DataContainer *cdata)
{
    //! The real meaning of the function with DataContainer is to read the file contents into cdata.
    //! For cdata, this is the process of writing.
    cdata->ReadFile(file);
}

LIB_EXPORT void ReadFile(ActionKey *actkey)
{
    fstream &file = *(actkey->file);
    ReadFile(file, actkey->GetData());
}

LIB_EXPORT void WriteFile(fstream &file, DataContainer *cdata)
{
    cdata->WriteFile(file);
}

LIB_EXPORT void WriteFile(ActionKey *actkey)
{
    fstream *file = actkey->file;

    if (!file) return;

    WriteFile(*file, actkey->GetData());
}

LIB_EXPORT void WriteASCIIFile(fstream &file, const string &str)
{
    file << str;
}

LIB_EXPORT void WriteASCIIFile(fstream &file, std::ostringstream &oss)
{
    WriteASCIIFile(file, oss.str());
}

LIB_EXPORT void WriteASCIIFile(fstream &file, DataContainer *cdata)
{
    string str;
    DataCharToString(cdata, str);

    WriteASCIIFile(file, str);
}

LIB_EXPORT void WriteASCIIFile(ActionKey *actionKey, std::ostringstream &oss)
{
    PHSPACE::ParallelOpenFile(actionKey);
    fstream &file = *(actionKey->GetFile());
    PHSPACE::WriteASCIIFile(file, oss);
    PHSPACE::ParallelCloseFile(file);
}

LIB_EXPORT void WriteASCIIFile(ActionKey *actkey)
{
    if (!actkey->file)
    {
        return;
    }

    fstream &file = *(actkey->file);
    WriteASCIIFile(file, actkey->GetData());
}

LIB_EXPORT void WriteBinaryFile(ActionKey *actkey)
{
    /*DumpToTecio dumptotecio(actkey);
    dumptotecio.Run();
    if (!actkey->file) return;
    fstream &file = *(actkey->file);
    WriteASCIIFile(file, actkey->GetData());*/
}

LIB_EXPORT void WriteSentinelFile()
{
    fstream sentinelfile;
    string sentinelfilename = "results/sentinel";
    sentinelfile.open(sentinelfilename.c_str(),ios::in);
    if (!sentinelfile)
    {
        fstream sentinelfile1;
        sentinelfile1.open(sentinelfilename.c_str(), ios_base::out|ios_base::trunc);
        sentinelfile1.close();
        sentinelfile1.clear();
    }
    else
    {
        sentinelfile.close();
        sentinelfile.clear();
    }
}

LIB_EXPORT void WriteScreen(DataContainer *cdata)
{
    string str;
    DataCharToString(cdata, str);
    cout << str;
}

LIB_EXPORT void WriteScreen(ActionKey *actkey)
{
    WriteScreen(actkey->GetData());
}

LIB_EXPORT bool IfFileEmpty(fstream &file)
{
    file.seekp(0, ios::end);
    streamoff i = file.tellp();

    if (i)
    {
        return false;
    }
    return true;
}

LIB_EXPORT bool FileExist(const string &fileName)
{
    fstream file;
    file.open(fileName.c_str(), ios::in);
    if (file)
    {
        CloseFile(file);
        return true;
    }
    else
    {
        return false;
    }
}

LIB_EXPORT bool SerialFileExist(const string &fileName)
{
    string acutalFileName = AddSymbolToFileName(fileName, "_", 0);
    return FileExist(acutalFileName);
}

LIB_EXPORT void ReadAbstractData(fstream &file, DataContainer *cdata, int send_proc, int recv_proc, int tag)
{
    using namespace PHMPI;

    int myid = GetCurrentProcessorID();

    if (myid == send_proc)
    {
        ReadFile(file, cdata);
    }

    PH_Trade(cdata, send_proc, recv_proc, tag);
}

LIB_EXPORT void ReadControlParameters()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
        PrintToWindow("Control parameters reading by server ...\n");
        ReadKeyInfo();
        ReadDefaultParamInfo();
        ReadParamInfo();
    }

    BroadcastInfo();
}

LIB_EXPORT void ReadKeyInfo()
{
    string keyfile = "./bin/key.hypara";

    ReadParamInfo(keyfile);
}

LIB_EXPORT void ReadDefaultParamInfo()
{
    string defaultParaFile;
    GlobalDataBase::GetData("defaultParaFile", &defaultParaFile, PHSTRING, 1);

    ReadParamInfo(defaultParaFile);
}

LIB_EXPORT void ReadParamInfo()
{
    string parafilename;
    GlobalDataBase::GetData("parafilename", &parafilename, PHSTRING, 1);

    int nparafile = GlobalDataBase::GetIntParaFromDB("nparafile");

    ReadParamInfo(parafilename);

    for (int ifile = 1; ifile < nparafile; ++ ifile)
    {
        std::ostringstream oss;
        oss << "parafilename" << ifile;
        string filetitle = oss.str();
        string filename;
        GlobalDataBase::GetData(filetitle, &filename, PHSTRING, 1);

        ReadParamInfo(filename);
    }
}

LIB_EXPORT void ReadParamInfo(const string &filename)
{
    fstream file;
    file.open(filename.c_str(), ios_base::in);

    if (!file)
    {
        TK_Exit::FileOpenErrorExit(filename);
    }

    TK_Parse::ReadBasicData(file, GlobalDataBase::GetDataPara());

    file.close();
    file.clear();
}

LIB_EXPORT void BroadcastInfo()
{
    using namespace PHMPI;
    int number_of_processor = GetNumberOfProcessor();
    if (number_of_processor <= 1) return;

    PrintToWindow("Control parameters broad casting by server ...\n");

    int serverTmp = GetServerProcessorID();
    int myid = GetCurrentProcessorID();

    DataContainer *cdata = new DataContainer();

    if (myid == serverTmp)
    {
        PHSPACE::GlobalDataBase::CompressData(cdata);
    }

    PH_BcastByServer(cdata);

    if (myid != serverTmp)
    {
        PHSPACE::GlobalDataBase::DecompressData(cdata);
    }

    delete cdata;    cdata = nullptr;
}

LIB_EXPORT void CompleteParamInfo()
{
    //! 1. Add walldistfile, its name is the same with grid file, except for the suffix ".wdt".
    //! 2. bnd_file is the same name with from_gfile when converting grid.
    int nsimutask = GlobalDataBase::GetIntParaFromDB("nsimutask");
    if (nsimutask == SOLVE_FIELD || nsimutask == CAL_WALL_DIST)
    {
        string origialGridFile = GlobalDataBase::GetStrParaFromDB("gridfile");
        string walldistfile = ChangeExtensionOfFileName(origialGridFile, "wdt");
        GlobalDataBase::UpdateData("walldistfile", &walldistfile, PHSTRING, 1);
    }
    else if (nsimutask == CREATE_GRID)
    {
        string from_gfile = GlobalDataBase::GetStrParaFromDB("from_gfile");
        string bnd_file = ChangeExtensionOfFileName(from_gfile, "inp");
        GlobalDataBase::UpdateData("bnd_file", &bnd_file, PHSTRING, 1);
    }
    else if (nsimutask == INTEGRATIVE_SOLVER)
    {
        string origialGridFile = GlobalDataBase::GetStrParaFromDB("gridfile");
        string walldistfile = ChangeExtensionOfFileName(origialGridFile, "wdt");
        GlobalDataBase::UpdateData("walldistfile", &walldistfile, PHSTRING, 1);

        string from_gfile = GlobalDataBase::GetStrParaFromDB("from_gfile");
        string bnd_file = ChangeExtensionOfFileName(from_gfile, "inp");
        GlobalDataBase::UpdateData("bnd_file", &bnd_file, PHSTRING, 1);
    }
    else
    {
        //! do nothing.
    }

    //! Some other default parameters.
    int nTurbulenceEquation = 0;

    int ViscousType = GlobalDataBase::GetIntParaFromDB("viscousType");

    if (ViscousType == 2)
    {
        nTurbulenceEquation = 0;
    }
    else if (ViscousType == 3)
    {
        nTurbulenceEquation = 1;
    }
    else if (ViscousType == 4)
    {
        nTurbulenceEquation = 2;
    }

    GlobalDataBase::UpdateData("n_turb", &nTurbulenceEquation, PHINT, 1);

    int transitionType = GlobalDataBase::GetIntParaFromDB("transitionType");
    int nTransitionEquation = 0;
    if (transitionType == 2)
    {
        nTransitionEquation = 2;
    }
    GlobalDataBase::UpdateData("n_transition", &nTransitionEquation, PHINT, 1);

    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    if (nChemical == 0)
    {
        nTemperatureModel = 1;
    }

    GlobalDataBase::UpdateData("ntmodel", &nTemperatureModel, PHINT, 1);

    int newnstep = 0;
    GlobalDataBase::UpdateData("newnstep", &newnstep, PHINT, 1);

    RDouble globalMinTimeStep = -LARGE;
    GlobalDataBase::UpdateData("globalMinTimeStep", &globalMinTimeStep, PHDOUBLE, 1);
}

LIB_EXPORT void PrintInflowConditionToWindow(void)
{

    int inflowParaType = GlobalDataBase::GetIntParaFromDB("inflowParaType");

    RDouble refGama                   = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble refMachNumber             = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    RDouble refReNumber               = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    RDouble attackd                   = GlobalDataBase::GetDoubleParaFromDB("attackd");
    RDouble angleSlide                = GlobalDataBase::GetDoubleParaFromDB("angleSlide");
    RDouble refDimensionalPressure    = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");
    RDouble refDimensionalDensity     = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    RDouble wallTemperature           = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");

    vector< RDouble > WallT;
    vector< int > WallNumber;
    ModifyWindowTemperature(wallTemperature, WallT, WallNumber);

    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");
    RDouble forceReferenceArea                 = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
    RDouble forceReferenceLength               = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
    RDouble TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    RDouble TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    RDouble TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

    using namespace GAS_SPACE;
    FluidParameter refParam;
    gas->GetReferenceParameters(refParam);
    RDouble refVelocity = refParam.GetVelocity();
    RDouble refVicosity = refParam.GetViscosity();
    RDouble refSoundSpeed = refParam.GetSoundVelocity();

    PrintToWindow("Inflow Conditions", "\n");
    PrintToWindow("  Gama                     :    ", refGama      ,           "\n");
    PrintToWindow("  Mach Number              :    ", refMachNumber,           "\n");
    PrintToWindow("  Attack Angle             :    ", attackd      , "degree", "\n");
    PrintToWindow("  Sideslip Angle           :    ", angleSlide   , "degree", "\n");
 
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    if (viscousType)
    {
        if (inflowParaType == FLIGHTCONDITION)
        {
            RDouble height = GlobalDataBase::GetDoubleParaFromDB("height");
            PrintToWindow("  Fly Height               :    ", height, "km", "\n");
        }
        PrintToWindow("  Reynolds Number          :    ", refReNumber,               " 1/m",   "\n");
        PrintToWindow("  Dimensional Pressure     :    ", refDimensionalPressure,    " Pa",    "\n");
        PrintToWindow("  Dimensional Density      :    ", refDimensionalDensity,     " kg/m^3", "\n");
        PrintToWindow("  Dimensional Temperature  :    ", refDimensionalTemperature, " K",     "\n");

        int metabolicWallTemperature = -1;
        GlobalDataBase::GetData("metabolicWallTemperature", &metabolicWallTemperature, PHINT, 1);

        if (metabolicWallTemperature == 0)
        {
            for (int i = 0; i < WallT.size(); ++i)
            {
                if (WallT[i] < 0.0)
                {
                    PrintToWindow("  Temperature of Wall",WallNumber[i],"   :    ",  "Adiabatic Wall", "\n");
                }
                else if (WallT[i] == 0.0)
                {
                    PrintToWindow("  Temperature of Wall",WallNumber[i],"   :    ", "Radiation balance Wall", "\n");
                }
                else
                {
                    PrintToWindow("  Temperature of Wall",WallNumber[i],"   :    ", WallT[i], "K", "\n");
                }
            }

        }
        PrintToWindow("  Dimensional Sound Speed  :    ", refSoundSpeed,               " m/s",    "\n");
        PrintToWindow("  Dimensional Velocity     :    ", refVelocity,               " m/s",      "\n");
        PrintToWindow("  Dimensional Viscosity    :    ", refVicosity,               " kg/(m*s)", "\n");
    }

    PrintToWindow("Geometry Informations"           ,                                     "\n");
    PrintToWindow("  Scale Factor             :    ", reynoldsReferenceLengthDimensional,  "\n");
    PrintToWindow("  Force Reference Area     :    ", forceReferenceArea                ,   "m^2", "\n");
    PrintToWindow("  Force Reference Length   :    ", forceReferenceLength              ,   "m", "\n");
    PrintToWindow("  Force Reference x        :    ", TorqueRefX                        , "\n");
    PrintToWindow("  Force Reference y        :    ", TorqueRefY                        , "\n");
    PrintToWindow("  Force Reference z        :    ", TorqueRefZ                        , "\n");
}

LIB_EXPORT void ModifyWindowTemperature(RDouble &wallTemperature, vector< RDouble > &WallT, vector< int > &WallNumber)
{
    /*if (wallTemperature >= 0.0)
    {
        return;
    }*/
    int nWall = 0;
    GlobalBoundaryCondition::ReadGlobalBoundaryCondition();
    vector<SimpleBC *> *globalBCList = GlobalBoundaryCondition::GetGlobalBoundaryConditionList();
    for (int iBC = 0; iBC < globalBCList->size(); ++ iBC)
    {
        SimpleBC *bcParamDataBase = (*globalBCList)[iBC];
        int bcType = bcParamDataBase->GetBCType();
        Data_Param *bcParamDB = bcParamDataBase->GetBCParamDataBase();
        if (bcType != PHENGLEI::SOLID_SURFACE && bcType != PHENGLEI::ABLATION_SURFACE)
        {
            continue;
        }
        else if (bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::ABLATION_SURFACE)
        {
            if (!bcParamDB->IsExist("wallTemperature", PHDOUBLE, 1))
            {
                GlobalDataBase::GetData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
                bcParamDB->UpdateData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
            }
            else
            {
                bcParamDB->GetData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
                //GlobalDataBase::UpdateData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
            }
            WallT.push_back(wallTemperature);
            nWall ++;
            WallNumber.push_back(nWall);
        }
    }
}
}