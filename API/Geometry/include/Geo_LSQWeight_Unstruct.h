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
//! @file      Geo_LSQWeight_Unstruct.h
//! @brief     It defines the Least Square Weight of the unstructured grid.
//! @author    Bell, He Xin.

#pragma once
#include "LIB_Macro.h"
#include "Precision.h"
namespace PHSPACE
{
//! @brief Geo_LSQWeight_Unstruct class defines Least Square Weight of unstructured grid.\n
class Geo_LSQWeight_Unstruct
{
private:

    RDouble *iwt;

    RDouble *ixx, *iyy, *izz;

    RDouble *ixy, *ixz, *iyz;

    char *fMark;
public:
    LIB_EXPORT Geo_LSQWeight_Unstruct();
    LIB_EXPORT ~Geo_LSQWeight_Unstruct();

public:
    //! Return the least square weights IWT.
    RDouble * GetLeastSquareIWT() const;

    //! Return the least square weights IXX.
    RDouble * GetLeastSquareIXX() const;

    //! Return the least square weights IYY.
    RDouble * GetLeastSquareIYY() const;

    //! Return the least square weights IZZ.
    RDouble * GetLeastSquareIZZ() const;

    //! Return the least square weights IXY.
    RDouble * GetLeastSquareIXY() const;

    //! Return the least square weights IXZ.
    RDouble * GetLeastSquareIXZ() const;

    //! Return the least square weights IYZ.
    RDouble * GetLeastSquareIYZ() const;

    //! Set the least square weights IWT.
    void SetLeastSquareIWT(RDouble *iwt);

    //! Set the least square weights IXX.
    void SetLeastSquareIXX(RDouble *ixx);

    //! Set the least square weights IYY.
    void SetLeastSquareIYY(RDouble *iyy);

    //! Set the least square weights IZZ.
    void SetLeastSquareIZZ(RDouble *izz);

    //! Set the least square weights IXY.
    void SetLeastSquareIXY(RDouble *ixy);

    //! Set the least square weights IXZ.
    void SetLeastSquareIXZ(RDouble *ixz);

    //! Set the least square weights IYZ.
    void SetLeastSquareIYZ(RDouble *iyz);

    //! Get the face mark.
    char * GetFaceMark() const;

    //! Set the face mark.
    void SetFaceMark(char *fMark);
};

#include "Geo_LSQWeight_Unstruct.hxx"
}