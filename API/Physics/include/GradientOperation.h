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
//! @file      Gradient_Operation.h
//! @brief     Gradient_Operation defines some methods carrying out calculations on gradient.
//! @author    Zhang Jian, He Xin.

#pragma once

namespace PHSPACE
{
class Grid;
class GreenGeo3D
{
public:
    RDouble nxip,nyip,nzip,nxim,nyim,nzim;
    RDouble nxjp,nyjp,nzjp,nxjm,nyjm,nzjm;
    RDouble nxkp,nykp,nzkp,nxkm,nykm,nzkm;
    RDouble ovol;
};
class Geometry3D
{
public:
    RDouble kx,ky,kz,ex,ey,ez,cx,cy,cz;
    RDouble vjacob;
};

void CorrectGradient(RDouble fl, RDouble fr, RDouble &dfdx, RDouble &dfdy, RDouble dx, RDouble dy, RDouble ods);
void CorrectGradient(RDouble fl, RDouble fr, RDouble &dfdx, RDouble &dfdy, RDouble &dfdz, RDouble dx, RDouble dy, RDouble dz, RDouble ods);

void Get_SPM  (int i, int j, int k, int js, int ks, int ndim, GreenGeo3D *greenGeo3D, RDouble4D &xfn, RDouble4D &yfn, RDouble4D &zfn, RDouble4D &area, RDouble3D &vol);
void Get_SPM  (int i, int j, int k, int js, int ks, int ndim, GreenGeo3D *greenGeo3D, RDouble4D &xfv, RDouble4D &yfv, RDouble4D &zfv, RDouble3D &vol);
void Get_SPM_I(int i, int j, int k, int js, int ks, int ndim, GreenGeo3D *greenGeo3D, RDouble4D &xfn, RDouble4D &yfn, RDouble4D &zfn, RDouble4D &area, RDouble3D &vol);
void Get_SPM_J(int i, int j, int k, int js, int ks, int ndim, GreenGeo3D *greenGeo3D, RDouble4D &xfn, RDouble4D &yfn, RDouble4D &zfn, RDouble4D &area, RDouble3D &vol);
void Get_SPM_K(int i, int j, int k, int js, int ks, GreenGeo3D *greenGeo3D, RDouble4D &xfn, RDouble4D &yfn, RDouble4D &zfn, RDouble4D &area, RDouble3D &vol);

void Get_SPM_PARA(int ndim,GreenGeo3D *greenGeo3D,RDouble xfn_1,RDouble yfn_1,RDouble zfn_1,RDouble area_1,RDouble xfn_ip1,RDouble yfn_ip1,RDouble zfn_ip1,RDouble area_ip1,RDouble xfn_2,RDouble yfn_2,RDouble zfn_2,RDouble area_2,RDouble xfn_ip2,RDouble yfn_ip2,RDouble zfn_ip2,RDouble area_ip2,RDouble xfn_3,RDouble yfn_3,RDouble zfn_3,RDouble area_3,RDouble xfn_ip3,RDouble yfn_ip3,RDouble zfn_ip3,RDouble area_ip3,RDouble data_vol);
void DXDYDZ  (RDouble3D & f,GreenGeo3D * g,int i,int j,int k,int js,int ks,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz);
void DXDYDZ_I(RDouble3D & f,GreenGeo3D * g,int i,int j,int k,int js,int ks,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz);
void DXDYDZ_J(RDouble3D & f,GreenGeo3D * g,int i,int j,int k,int js,int ks,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz);
void DXDYDZ_K(RDouble3D & f,GreenGeo3D * g,int i,int j,int k,int js,int ks,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz);
void DXDYDZ  (RDouble4D & f,GreenGeo3D * g,int i,int j,int k,int js,int ks,int npos,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz);
void DXDYDZ_I(RDouble4D & f,GreenGeo3D * g,int i,int j,int k,int js,int ks,int npos,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz);
void DXDYDZ_J(RDouble4D & f,GreenGeo3D * g,int i,int j,int k,int js,int ks,int npos,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz);
void DXDYDZ_K(RDouble4D & f,GreenGeo3D * g,int i,int j,int k,int js,int ks,int npos,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz);

void GradCenter(Grid *grid_in, RDouble3D &q, RDouble3D &dqdx, RDouble3D &dqdy, RDouble3D &dqdz);
void GradCenter(Grid *grid_in, RDouble4D &q, RDouble3D &dqdx, RDouble3D &dqdy, RDouble3D &dqdz, int m);

void Get_GEO_I(int i, int j, int k, int im, int ip, int ndim, Geometry3D *g, RDouble4D &xfv, RDouble4D &yfv, RDouble4D &zfv, RDouble3D &vol);
void Get_GEO_J(int i, int j, int k, int jm, int jp, int ndim, Geometry3D *g, RDouble4D &xfv, RDouble4D &yfv, RDouble4D &zfv, RDouble3D &vol);
void Get_GEO_K(int i, int j, int k, int km, int kp, Geometry3D *g, RDouble4D &xfv, RDouble4D &yfv, RDouble4D &zfv, RDouble3D &vol);
void DXDYDZ_I(RDouble4D & f,Geometry3D * g,int i,int j,int k,int ndim,int npos,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz,int lr);
void DXDYDZ_J(RDouble4D & f,Geometry3D * g,int i,int j,int k,int ndim,int npos,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz,int lr);
void DXDYDZ_K(RDouble4D & f,Geometry3D * g,int i,int j,int k,int npos,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz,int lr);
}
