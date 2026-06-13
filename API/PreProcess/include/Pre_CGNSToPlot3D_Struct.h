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
//! @file      Pre_CGNSToPlot3D_Struct.h
//! @brief     Grid conversion from CGNS to Plot3D.
//! @author    Refactored by He Xianyao, Xu Gang, He Xin.

#pragma once
#include "LIB_Macro.h"
#include "Pre_CGNSConversion_Struct.h"

using namespace std;

namespace PHSPACE
{

class CGNS_Str_Data;

//! @brief Pre_CGNSToPlot3D_Struct class defines the method of Grid conversion from CGNS format to Plot3D.
//! Pre_CGNSToPlot3D_Struct
//! CGNS_Str_Data, BaseData
class Pre_CGNSToPlot3D_Struct : public Pre_CGNSConversion_Struct
{
public:
    //! @param[in] gridFileName    Grid file name of CGNS grid.
    //! mode = 0 (output .grd, .inp files); mode(Default) = 1 (output .grd, .inp, .bc, .bcname files with the name "_0").
    LIB_EXPORT Pre_CGNSToPlot3D_Struct(const string &gridFileName);
    LIB_EXPORT ~Pre_CGNSToPlot3D_Struct();

private:
    void Conversion();

private:
    //! Write the plot3d format grid file.
    void WriteToGrd();

    //! Write the boundary condition file.
    void WriteToInp();

    //! Write the boundary face information into file.
    void WriteToBcinformation();

    //!
    void writetoCGNS(fstream &bcfile, cgsize_t *pnts, int lab, string bcName, int nCoords);

private:
    string cel_file;
    string bnd_file;
    string bcinfor_file;

};

}