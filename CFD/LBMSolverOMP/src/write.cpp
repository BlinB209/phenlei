#include "LBMSolverOMP.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h>
#include <iostream>
#include <Windows.h>
#include <time.h>
#include <fstream>

using namespace std;

void dumpstring(char *str, FILE *file);

void LBM::write_binary(std::string filename)    //! int im, int jm, int km
{
    int VarN=7; //Number of variables(VarN) 
    //char filename[100];
    //sprintf_s(filename, 100, "%s", "P.plt");  //sprintf(filename, "P.plt");
    FILE *fp;
    errno_t err;
    if ((err = fopen_s(&fp, filename.c_str(), "wb")) != 0) std::cout << "error to open file: " << filename << '\n';
    int im = this->NX; 
    int jm = this->NY; 
    int km = this->NZ;
    
    char Title[] = "TubeFlow";
    char Vname[7][5] = { "X","Y","Z","u","v","w","p" };    //! "XYZp";

    char Zonename1[] = "Zone 001";
    float ZONEMARKER = 299.0;
    float EOHMARKER = 357.0;

    //! ==============Header Secontion =================//
    //! ------1.1 Magic number, Version number
    char MagicNumber[] = "#!TDV101";
    //cout<< "ok" << sizeof(MagicNumber) << endl;
    fwrite(MagicNumber, 8, 1, fp);

    //! ---- - 1.2.Integer value of 1.
    int IntegerValue = 1;
    fwrite(&IntegerValue, sizeof(IntegerValue), 1, fp);

    //! ---- - 1.3.Title and variable names.
    //! ---- - 1.3.1.The TITLE.
    dumpstring(Title, fp);
    fwrite(&VarN, sizeof(VarN), 1, fp);

    //! ------1.3.3 Variable names.
    for (int i = 0; i < VarN; i++)
    {
        dumpstring(Vname[i], fp);
    }

    //! ---- - 1.4.Zones.
    //! --------Zone marker.Value = 299.0
    fwrite(&ZONEMARKER, 1, sizeof(ZONEMARKER), fp);
    //--------Zone name.
    //fwrite(Zonename1, sizeof(Zonename1), 1, fp);
    dumpstring(Zonename1, fp);

    int ZoneStuff[5] = { -1,0,1,0,0 };
    //! Zone color
    //! ZoneType
    //! DaraPacking 0=Block, 1=Point
    //! Specify Var Location. 0 = Don't specify, all c_str is located at the nodes. 1 = Specify
    //! Number of user defined face neighbor connections(value >= 0)
    for (int i = 0; i < 5; i++)
    {
        fwrite(&ZoneStuff[i], sizeof(int), 1, fp);
    }

    //! im, jm, km
    fwrite(&im, sizeof(im), 1, fp);
    fwrite(&jm, sizeof(jm), 1, fp);
    fwrite(&km, sizeof(km), 1, fp);

    //! -1 = Auxiliary name / value pair to follow 0 = No more Auxiliar name / value pairs.
    int AuxiliaryName = 0;
    fwrite(&AuxiliaryName, sizeof(AuxiliaryName), 1, fp);

    //! ----I HEADER OVER----------------------------------------------
    //! =============================Geometries section================
    //! =============================Text section======================
    //! EOHMARKER, value = 357.0
    fwrite(&EOHMARKER, sizeof(EOHMARKER), 1, fp);

    //! ================II.Data section===============
    //! ---- - 2.1 zone-------------------------------
    fwrite(&ZONEMARKER, sizeof(ZONEMARKER), 1, fp);

    //! variable format, 1 = Float, 2 = Double, 3 = LongInt, 4 = ShortInt, 5 = Byte, 6 = Bit
    int VariableType[7] = { 3,3,3,2,2,2,2 };
    for (int i = 0; i < VarN; i++) 
    {
        fwrite(&VariableType[i], sizeof(int), 1, fp);
    }

    //! Has variable sharing 0 = no, 1 = yes.
    int HasVarSharing = 0;
    fwrite(&HasVarSharing, sizeof(HasVarSharing), 1, fp);
    //! Zone number to share connectivity list with(-1 = no sharing).
    int ZoneNumToShareConnectivity = -1;
    fwrite(&ZoneNumToShareConnectivity, sizeof(ZoneNumToShareConnectivity), 1, fp);
    //! Zone c_str.Each variable is in c_str format asspecified above.

    for (int k = 0; k < km; k++)
    {
        for (int j = 0; j < jm; j++)
        {
            for (int i = 0; i < im; i++)
            {
                int Va1 = i;
                int Va2 = j;
                int Va3 = k;
                fwrite(&Va1, sizeof(Va1), 1, fp);
                fwrite(&Va2, sizeof(Va2), 1, fp);
                fwrite(&Va3, sizeof(Va3), 1, fp);
                fwrite(&vel[this->index(i, j, k)].x, sizeof(double), 1, fp);
                fwrite(&vel[this->index(i, j, k)].y, sizeof(double), 1, fp);
                fwrite(&vel[this->index(i, j, k)].z, sizeof(double), 1, fp);
                fwrite(&rho[this->index(i, j, k)], sizeof(double), 1, fp);
            }
        }
    }
    fclose(fp);
    //getchar();
}

//! Ð´Èë×Ö·û´®º¯Êý
void dumpstring(char *str, FILE *file)
{
    int value = 0;
    while ((*str) != '\0')
    {
        value = (int)*str;
        fwrite(&value, sizeof(int), 1, file);
        str++;
    }
    char null_char[] = "";
    value = (int) * null_char;
    fwrite(&value, sizeof(int), 1, file);
}