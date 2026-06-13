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
//! @file      Gradient.h
//! @brief     Gradient defines a class that compute and store the gradient of specific variables. 
//! @author    Zhang Jian, He Xin.

#pragma once
#include "LIB_Macro.h"
#include "Precision.h"

namespace PHSPACE
{
//! Forward declaration.
class Grid;
class FieldProxy;
class Range;
class GradientCellCenter;

//@brief 
class Gradient
{

private:
    //! Grid used to compute gradient.
    Grid *grid;

    //! Grid type.
    int gridType;

    //! The variable field used to compute its gradient. 
    FieldProxy *q;

    //! The total number of the variables.
    int nTotalVariable;

    //! The computed gradient.
    FieldProxy *dqdx, *dqdy, *dqdz;

    //! The node value used for computing the gradient.
    RDouble **nodeValue;

public:
    //! Constructor for gradient computation.
    LIB_EXPORT Gradient(Grid *grid_in, int gridType, FieldProxy *q, int nTotalVariable);

    //! Destructor.
    LIB_EXPORT ~Gradient();

public:
    //! set the nodevalue of this class. 
    //! @param[in] nodeValue:   The node value used for computing the gradient.
    LIB_EXPORT void setNodeValue(RDouble **nodeValue);

    LIB_EXPORT void CalculatebyCellCenter(GradientCellCenter * gradientCellCenter, int nsurf);

    //! Return the x component of gradient, i.e dqdx.
    LIB_EXPORT FieldProxy * GetGradX() const;
    
    //! Return the y component of gradient, i.e dqdy.
    LIB_EXPORT FieldProxy * GetGradY() const;

    //! Return the z component of gradient, i.e dqdz.
    LIB_EXPORT FieldProxy * GetGradZ() const;

    //! Get the two dimensional pointer of the gradient(dqdx), only unstructured grid
    LIB_EXPORT RDouble ** GetGradX_UNS() const;

    //! Get the two dimensional pointer of the gradient(dqdx), only unstructured grid
    LIB_EXPORT RDouble ** GetGradY_UNS() const;

    //! Get the two dimensional pointer of the gradient(dqdx), only unstructured grid
    LIB_EXPORT RDouble ** GetGradZ_UNS() const;

    //! Communicate gradient at interface.
    //! @param[in] gradNameX name of gradient's x-component, e.g "dqdx".
    //! @param[in] gradNameY name of gradient's y-component, e.g "dqdy".
    //! @param[in] gradNameZ name of gradient's z-component, e.g "dqdz".
    //! @note gradient name must correspond to user registered through RegisterInterfaceField.
    LIB_EXPORT void CommunicateInterfaceValue(string gradNameX, string gradNameY, string gradNameZ);

    //! Store gradient at boundary's left face to compute aerodynamic coefficients.
    //! @param[in] boundGradNameX name of gradient's x-component at boundary's left face.
    //! @param[in] boundGradNameY name of gradient's y-component at boundary's left face.
    //! @param[in] boundGradNameZ name of gradient's z-component at boundary's left face.
    //! @note gradient name must correspond to user allocated memory through AllocateGlobalVar.
    //! @param[in] nsurf    if equal to zero, for unstructured grid gradient boundary storage.
    //!                     if not equal zero, for struct grid's gradient storage at the surfaces of direction nsurf
    //! -# 0: default option, for unstructured grid.
    //! -# 1: surfaces of direction x
    //! -# 2: surfaces of direction y
    //! -# 3: surfaces of direction z
    LIB_EXPORT void StoreBoundGrad(string boundGradNameX, string boundGradNameY, string boundGradNameZ, int nsurf = 0);

    LIB_EXPORT int GetNumberOfTotalVariable();

private:
    //! Allocate memory for dqdx, dqdy and dqdz.
    void Allocate();

    //! Deallocate memory for dqdx, dqdy and dqdz.
    void FreeMemory();

    //! Calculate gradient of variables at the surfaces of direction nsurf
    //! rewrite from old Grad3dVector of fygradientoperation.h
    void Grad3DVector(int nsurf);
};

class GradientCellCenter
{

private:
    //! Grid used to compute gradient.
    Grid *grid;

    //! Grid type.
    int gridType;

    //! The variable field used to compute its gradient. 
    FieldProxy *q;

    //! The total number of the variables.
    int nTotalVariable;

    //! The computed gradient.
    FieldProxy *dqdx, *dqdy, *dqdz;

public:
    //! Constructor for gradient computation.
    LIB_EXPORT GradientCellCenter(Grid *grid_in, int gridType, FieldProxy *q, int nTotalVariable);

    //! Destructor.
    LIB_EXPORT ~GradientCellCenter();

public:
    //! Calculate the gradient, 
    LIB_EXPORT void Calculate();

    //! Return the x component of gradient, i.e dqdx.
    LIB_EXPORT FieldProxy * GetGradX() const;
    
    //! Return the y component of gradient, i.e dqdy.
    LIB_EXPORT FieldProxy * GetGradY() const;

    //! Return the z component of gradient, i.e dqdz.
    LIB_EXPORT FieldProxy * GetGradZ() const;

    LIB_EXPORT void SetProxy(FieldProxy *q);

    LIB_EXPORT int GetNumberOfTotalVariable();

private:
    //! Allocate memory for dqdx, dqdy and dqdz.
    void Allocate();

    //! Deallocate memory for dqdx, dqdy and dqdz.
    void FreeMemory();

    //! Calculate gradient of variables at the cell center.
    void GradCenter();
};

#include "Gradient.hxx"
}