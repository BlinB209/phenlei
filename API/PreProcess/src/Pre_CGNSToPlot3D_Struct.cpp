#include "Pre_CGNSToPlot3D_Struct.h"
#include "PHIO.h"
#include "TK_Parse.h"
#include "Constants.h"

namespace PHSPACE
{

Pre_CGNSToPlot3D_Struct::Pre_CGNSToPlot3D_Struct(const string &gridFileName) :
    Pre_CGNSConversion_Struct(gridFileName)
{

}

Pre_CGNSToPlot3D_Struct::~Pre_CGNSToPlot3D_Struct()
{

}

void Pre_CGNSToPlot3D_Struct::Conversion()
{
    ResetGridScaleAndTranslate();

    WriteToGrd();
    WriteToInp();

    int dumpOldGrid = GlobalDataBase::GetIntParaFromDB("dumpOldGrid");
    if (dumpOldGrid)
    {
        WriteToBcinformation();
    }
}

void Pre_CGNSToPlot3D_Struct::WriteToGrd()
{
    string fname = " ", fext = " ";
    GetNameExt(gridFileName, fname, fext, ".");
    cel_file = fname + "_0.grd";
    GlobalDataBase::UpdateData("from_gfile", &cel_file, PHSTRING, 1);

    //! Creating Grid file.
    fstream Gridgenfile;
    Gridgenfile.open(cel_file.c_str(), ios_base::out|ios_base::binary|ios_base::trunc);
    Gridgenfile.write(reinterpret_cast<char *>(&nBlocks), sizeof(int));

    //! Import the maximum value of i, j, k.
    for (int iZone = 1; iZone <= nBlocks; ++ iZone)
    {
        BaseData *base_cgns = Str_CGNS_Data->GetBaseData(iZone-1);

        int idim = base_cgns->GetIdim();
        int jdim = base_cgns->GetJdim();
        int kdim = base_cgns->GetKdim();

        Gridgenfile.write(reinterpret_cast<char *>(&idim), sizeof(int));
        Gridgenfile.write(reinterpret_cast<char *>(&jdim), sizeof(int));
        Gridgenfile.write(reinterpret_cast<char *>(&kdim), sizeof(int));
    }

    //! Import the coordinates of node in the CGNS Grid.
    for (int iZone = 1; iZone <= nBlocks; ++ iZone)
    {
        BaseData *base_cgns = Str_CGNS_Data->GetBaseData(iZone-1);

        int nCoords =base_cgns->GetNCoords();

        int idim = base_cgns->GetIdim();
        int jdim = base_cgns->GetJdim();
        int kdim = base_cgns->GetKdim();

        RDouble ***x = base_cgns->GetX();
        RDouble ***y = base_cgns->GetY();
        RDouble ***z = base_cgns->GetZ();

        for (int k = 0; k < kdim; ++ k)
        {
            for (int j = 0; j < jdim; ++ j)
            {
                for (int i = 0; i < idim; ++ i)
                {
                    Gridgenfile.write(reinterpret_cast<char *>(&x[k][j][i]), sizeof(RDouble));
                }
            }
        }

        for (int k = 0; k < kdim; ++ k)
        {
            for (int j = 0; j < jdim; ++ j)
            {
                for (int i = 0; i < idim; ++ i)
                {
                    Gridgenfile.write(reinterpret_cast<char *>(&y[k][j][i]), sizeof(RDouble));
                }
            }
        }

        for (int k = 0; k < kdim; ++ k)
        {
            for (int j = 0; j < jdim; ++ j)
            {
                for (int i = 0; i < idim; ++ i)
                {
                    if (nCoords == 2) z[k][j][i] = 0.0;
                    Gridgenfile.write(reinterpret_cast<char *>(&z[k][j][i]), sizeof(RDouble));
                }
            }
        }
    }

    Gridgenfile.close();
}

void Pre_CGNSToPlot3D_Struct::WriteToInp()
{
    string fname = " ", fext = " ";
    GetNameExt(gridFileName, fname, fext, ".");
    bnd_file = fname + "_0.inp";
    GlobalDataBase::UpdateData("bnd_file",   &bnd_file, PHSTRING, 1);

    //! Create Inp file.
    fstream GridgenBCfile;
    GridgenBCfile.open(bnd_file.c_str(), ios_base::out);
    GridgenBCfile << 1 << endl;
    GridgenBCfile << nBlocks << endl;

    for (int iZone = 1; iZone <= nBlocks; ++ iZone)
    {
        BaseData *base_cgns = Str_CGNS_Data->GetBaseData(iZone-1);

        int nCoords = base_cgns->GetNCoords();

        int idim = base_cgns->GetIdim();
        int jdim = base_cgns->GetJdim();
        int kdim = base_cgns->GetKdim();
        char_33 *zonename = base_cgns->GetZonenames();

        //! Set the minimum data length.
        int wordwidth = 8;

        GridgenBCfile << setw(wordwidth) << idim;
        GridgenBCfile << setw(wordwidth) << jdim;
        if (nCoords == 3)
        {
            GridgenBCfile << setw(wordwidth) << kdim;
        }
        GridgenBCfile << endl;
        GridgenBCfile << zonename[0] << endl;

        int nconnect = base_cgns->GetNconnect();
        GridgenBCfile <<  nconnect << endl;

        int    nBCRegions    = base_cgns->GetNBCRegions();
        int    *lab1         = base_cgns->GetLab1();
        string *bcName       = base_cgns->GetBcName();
        cgsize_t_60000 *pnts = base_cgns->GetPnts();

        for (int boco = 1; boco <= nBCRegions; ++ boco)
        {
            writetoCGNS(GridgenBCfile, pnts[boco - 1], lab1[boco - 1], bcName[boco - 1], nCoords);
        }

        int ier = base_cgns->GetIer();
        if (ier == CG_OK)
        {
            int n1to1 = base_cgns->GetN1to1();
            for (int one21 = 1; one21 <= n1to1; ++ one21)
            {
                int *b = base_cgns->GetB(one21-1);
                int *vt = base_cgns->GetVt(one21-1);
                cgsize_t *sranges = base_cgns->GetSrange(one21-1);
                cgsize_t *donor_ranges = base_cgns->GetDonor_range(one21-1);

                //! Import FarField Boundary Condition and Wall Boundary Condition.
                for (int iCoords = 0; iCoords < nCoords; ++ iCoords)
                {
                    GridgenBCfile << setw(8) << sranges[iCoords] * b[iCoords];
                    GridgenBCfile << setw(8) << sranges[iCoords + nCoords] * b[iCoords];
                }

                //! Import interface information.
                GridgenBCfile << setw(8) << -1;
                GridgenBCfile << "   " << "Connection" << endl;

                for (int iCoords = 0; iCoords < nCoords; ++ iCoords)
                {
                    GridgenBCfile << setw(8) << donor_ranges[iCoords] * vt[iCoords];
                    GridgenBCfile << setw(8) << donor_ranges[iCoords + nCoords] * vt[iCoords];
                }

                int lab2 = base_cgns->Getlab2(one21-1);

                //! Import the information of symmetrical interface.
                GridgenBCfile << setw(8) << lab2 << endl;
            }
        }
        else
        {
            cout << "Warning: " << "Zone " << (iZone-1) << " has no interface!" << endl;
        }
    }

    GridgenBCfile.close();
}

void Pre_CGNSToPlot3D_Struct::WriteToBcinformation()
{
    string fname = " ", fext = " ";
    GetNameExt(gridFileName, fname, fext, ".");
    bcinfor_file = fname + "_0.bc";

    fstream file;
    PHSPACE::OpenFile(file, bcinfor_file, ios_base::out|ios_base::binary);

    PHWrite(file, nBlocks);

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        VirtualFile *virtualFile = new VirtualFile(&file);
        virtualFile->BeginWriteWork();

        PHWrite(virtualFile, iZone);

        BaseData *base_cgns = Str_CGNS_Data->GetBaseData(iZone);

        int nBCRegions = base_cgns->GetNBCRegions();
        int nconnect = base_cgns->GetNconnect();
        PHWrite(virtualFile, nconnect);

        string *bcNameList = new string [nconnect];
        string *bcName = base_cgns->GetBcName();
        uint_long totalSize = 0;
        for (int boco = 0; boco < nBCRegions; ++ boco)
        {
            bcNameList[boco] = bcName[boco];

            totalSize += static_cast< uint_long > (bcNameList[boco].size());
            totalSize += 1;
        }

        int ier = base_cgns->GetIer();
        if (ier == CG_OK)
        {
            int n1to1 = base_cgns->GetN1to1();
            for (int one21 = 0; one21 < n1to1; ++ one21)
            {
                bcNameList[nBCRegions+one21] = "Connection";

                totalSize += static_cast< uint_long > (bcNameList[nBCRegions+one21].size());
                totalSize += 1;
            }
        }
        else
        {
            cout << "Warning: " << "Zone " << iZone << " has no interface!" << endl;
        }

        char *bcNameChar = new char[totalSize];
        unsigned int count = 0;
        for (int boco = 0; boco < nconnect; ++ boco)
        {
            string &bcName1 = bcNameList[boco];
            streamsize nameSize = bcName1.size();
            for (unsigned int iChar = 0; iChar < nameSize; ++ iChar)
            {
                bcNameChar[count ++] = bcName1[iChar];
            }
            bcNameChar[count ++] = '\0';
        }

        PHWrite(virtualFile, totalSize);
        PHWrite(virtualFile, bcNameChar, static_cast<int>(totalSize));

        virtualFile->EndWriteWork();
        delete virtualFile;
        delete [] bcNameList;
        delete [] bcNameChar;
    }

    PHSPACE::CloseFile(file);
}

void Pre_CGNSToPlot3D_Struct::writetoCGNS(fstream &bcfile, cgsize_t *pnts, int lab, string bcName, int nCoords)
{
    if (nCoords == 2)
    {
        bcfile << setw(8) << pnts[0];
        bcfile << setw(8) << pnts[2];
        bcfile << setw(8) << pnts[1];
        bcfile << setw(8) << pnts[3];
        bcfile << setw(8) << lab;
    }
    else
    {
        bcfile << setw(8) << pnts[0];
        bcfile << setw(8) << pnts[3];
        bcfile << setw(8) << pnts[1];
        bcfile << setw(8) << pnts[4];
        bcfile << setw(8) << pnts[2];
        bcfile << setw(8) << pnts[5];
        if (bcName.substr(0,7) == "BndCond")
        {
            bcfile << setw(8) << atoi(bcName.substr(8, 2).c_str()) - 17;
        }
        else
        {
            bcfile << setw(8) << lab;
        }
    }

    bcfile << "   " << bcName;
    bcfile << endl;
}


}