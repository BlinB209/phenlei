#include "TK_Log.h"
#include "GradientOperation.h"
#include "Geo_StructGrid.h"
#include "Glb_Dimension.h"
#include "Constants.h"

namespace PHSPACE
{

void CorrectGradient(RDouble fl, RDouble fr, RDouble &dfdx, RDouble &dfdy, RDouble dx, RDouble dy, RDouble ods)
{
    RDouble fold, fnew, df;
    fold = dfdx * dx + dfdy * dy;
    fnew = (fr - fl) * ods;
    df   = fnew - fold;

    dfdx += df * dx;
    dfdy += df * dy;
}

void CorrectGradient(RDouble fl, RDouble fr, RDouble &dfdx, RDouble &dfdy, RDouble &dfdz, RDouble dx, RDouble dy, RDouble dz, RDouble ods)
{
    RDouble fold, fnew, df;
    fold = dfdx * dx + dfdy * dy + dfdz * dz;
    fnew = (fr - fl) * ods;
    df   = fnew - fold;

    dfdx += df * dx;
    dfdy += df * dy;
    dfdz += df * dz;
}

void Get_SPM(int i, int j, int k, int js, int ks, int ndim, GreenGeo3D *greenGeo3D, RDouble4D &xfn, RDouble4D &yfn, RDouble4D &zfn, RDouble4D &area, RDouble3D &vol)
{
    //int ip = i + 1;
    //int jp = j + 1;
    //int kp = k + 1;

    int ip = i + 1;
    int jp = j + js;
    int kp = k + ks;

    greenGeo3D->nxip = xfn(ip,j ,k ,1) * area(ip,j ,k ,1);
    greenGeo3D->nyip = yfn(ip,j ,k ,1) * area(ip,j ,k ,1);
    greenGeo3D->nzip = zfn(ip,j ,k ,1) * area(ip,j ,k ,1);
    greenGeo3D->nxim = xfn(i ,j ,k ,1) * area(i ,j ,k ,1);
    greenGeo3D->nyim = yfn(i ,j ,k ,1) * area(i ,j ,k ,1);
    greenGeo3D->nzim = zfn(i ,j ,k ,1) * area(i ,j ,k ,1);

    greenGeo3D->nxjp = xfn(i ,jp,k ,2) * area(i ,jp,k ,2);
    greenGeo3D->nyjp = yfn(i ,jp,k ,2) * area(i ,jp,k ,2);
    greenGeo3D->nzjp = zfn(i ,jp,k ,2) * area(i ,jp,k ,2);
    greenGeo3D->nxjm = xfn(i ,j ,k ,2) * area(i ,j ,k ,2);
    greenGeo3D->nyjm = yfn(i ,j ,k ,2) * area(i ,j ,k ,2);
    greenGeo3D->nzjm = zfn(i ,j ,k ,2) * area(i ,j ,k ,2);

    if (ndim == THREE_D)
    {
        greenGeo3D->nxkp = xfn(i ,j ,kp,3) * area(i ,j ,kp,3);
        greenGeo3D->nykp = yfn(i ,j ,kp,3) * area(i ,j ,kp,3);
        greenGeo3D->nzkp = zfn(i ,j ,kp,3) * area(i ,j ,kp,3);
        greenGeo3D->nxkm = xfn(i ,j ,k ,3) * area(i ,j ,k ,3);
        greenGeo3D->nykm = yfn(i ,j ,k ,3) * area(i ,j ,k ,3);
        greenGeo3D->nzkm = zfn(i ,j ,k ,3) * area(i ,j ,k ,3);
    }
    else
    {
        greenGeo3D->nxkp = 0.0;
        greenGeo3D->nykp = 0.0;
        greenGeo3D->nzkp = 0.0;
        greenGeo3D->nxkm = 0.0;
        greenGeo3D->nykm = 0.0;
        greenGeo3D->nzkm = 0.0;

    }
    greenGeo3D->ovol = one / vol(i,j,k);
}
void Get_SPM(int i, int j, int k, int js, int ks, int ndim, GreenGeo3D *greenGeo3D, RDouble4D &xfv, RDouble4D &yfv, RDouble4D &zfv, RDouble3D &vol)
{
    //int ip = i + 1;
    //int jp = j + 1;
    //int kp = k + 1;

    int ip = i + 1;
    int jp = j + js;
    int kp = k + ks;

    greenGeo3D->nxip = xfv(ip,j ,k ,1);
    greenGeo3D->nyip = yfv(ip,j ,k ,1);
    greenGeo3D->nzip = zfv(ip,j ,k ,1);
    greenGeo3D->nxim = xfv(i ,j ,k ,1);
    greenGeo3D->nyim = yfv(i ,j ,k ,1);
    greenGeo3D->nzim = zfv(i ,j ,k ,1);

    greenGeo3D->nxjp = xfv(i ,jp,k ,2);
    greenGeo3D->nyjp = yfv(i ,jp,k ,2);
    greenGeo3D->nzjp = zfv(i ,jp,k ,2);
    greenGeo3D->nxjm = xfv(i ,j ,k ,2);
    greenGeo3D->nyjm = yfv(i ,j ,k ,2);
    greenGeo3D->nzjm = zfv(i ,j ,k ,2);

    if (ndim == THREE_D)
    {
        greenGeo3D->nxkp = xfv(i ,j ,kp,3);
        greenGeo3D->nykp = yfv(i ,j ,kp,3);
        greenGeo3D->nzkp = zfv(i ,j ,kp,3);
        greenGeo3D->nxkm = xfv(i ,j ,k ,3);
        greenGeo3D->nykm = yfv(i ,j ,k ,3);
        greenGeo3D->nzkm = zfv(i ,j ,k ,3);
    }
    else
    {
        greenGeo3D->nxkp = 0.0;
        greenGeo3D->nykp = 0.0;
        greenGeo3D->nzkp = 0.0;
        greenGeo3D->nxkm = 0.0;
        greenGeo3D->nykm = 0.0;
        greenGeo3D->nzkm = 0.0;

    }
    greenGeo3D->ovol = one / vol(i,j,k);
}

void Get_SPM_PARA(int ndim,GreenGeo3D *greenGeo3D,RDouble xfn_1,RDouble yfn_1,RDouble zfn_1,RDouble area_1,RDouble xfn_ip1,RDouble yfn_ip1,RDouble zfn_ip1,RDouble area_ip1,RDouble xfn_2,RDouble yfn_2,RDouble zfn_2,RDouble area_2,RDouble xfn_ip2,RDouble yfn_ip2,RDouble zfn_ip2,RDouble area_ip2,RDouble xfn_3,RDouble yfn_3,RDouble zfn_3,RDouble area_3,RDouble xfn_ip3,RDouble yfn_ip3,RDouble zfn_ip3,RDouble area_ip3,RDouble data_vol)
{
        greenGeo3D->nxip = xfn_ip1 * area_ip1;
        greenGeo3D->nyip = yfn_ip1 * area_ip1;
        greenGeo3D->nzip = zfn_ip1 * area_ip1;
        greenGeo3D->nxim = xfn_1 * area_1;
        greenGeo3D->nyim = yfn_1 * area_1;
        greenGeo3D->nzim = zfn_1 * area_1;

        greenGeo3D->nxjp = xfn_ip2 * area_ip2;
        greenGeo3D->nyjp = yfn_ip2 * area_ip2;
        greenGeo3D->nzjp = zfn_ip2 * area_ip2;
        greenGeo3D->nxjm = xfn_2 * area_2;
        greenGeo3D->nyjm = yfn_2 * area_2;
        greenGeo3D->nzjm = zfn_2 * area_2;
        if (ndim == THREE_D)
        {
                greenGeo3D->nxkp = xfn_ip3 * area_ip3;
                greenGeo3D->nykp = yfn_ip3 * area_ip3;
                greenGeo3D->nzkp = zfn_ip3 * area_ip3;
                greenGeo3D->nxkm = xfn_3 * area_3;
                greenGeo3D->nykm = yfn_3 * area_3;
                greenGeo3D->nzkm = zfn_3 * area_3;
        }
        else
        {
                greenGeo3D->nxkp = 0.0;
                greenGeo3D->nykp = 0.0;
                greenGeo3D->nzkp = 0.0;
                greenGeo3D->nxkm = 0.0;
                greenGeo3D->nykm = 0.0;
                greenGeo3D->nzkm = 0.0;
        }
        greenGeo3D->ovol = one / data_vol;
}
void Get_SPM_I(int i, int j, int k, int js, int ks, int ndim, GreenGeo3D *greenGeo3D, RDouble4D &xfn, RDouble4D &yfn, RDouble4D &zfn, RDouble4D &area, RDouble3D &vol)
{
    int im,jm,km,ip,jp,kp;

    //ip = i + 1;
    //im = i - 1;
    //jp = j + 1;
    //jm = j;
    //kp = k + 1;
    //km = k;

    ip = i + 1;
    im = i - 1;
    jp = j + js;
    jm = j;
    kp = k + ks;
    km = k;

    greenGeo3D->nxip = half * (xfn(i ,j ,k,1) * area(i ,j ,k,1) + xfn(ip,j ,k,1) * area(ip,j ,k,1));
    greenGeo3D->nyip = half * (yfn(i ,j ,k,1) * area(i ,j ,k,1) + yfn(ip,j ,k,1) * area(ip,j ,k,1));
    greenGeo3D->nzip = half * (zfn(i ,j ,k,1) * area(i ,j ,k,1) + zfn(ip,j ,k,1) * area(ip,j ,k,1));
    greenGeo3D->nxim = half * (xfn(i ,j ,k,1) * area(i ,j ,k,1) + xfn(im,j ,k,1) * area(im,j ,k,1));
    greenGeo3D->nyim = half * (yfn(i ,j ,k,1) * area(i ,j ,k,1) + yfn(im,j ,k,1) * area(im,j ,k,1));
    greenGeo3D->nzim = half * (zfn(i ,j ,k,1) * area(i ,j ,k,1) + zfn(im,j ,k,1) * area(im,j ,k,1));
                                                                                    
    greenGeo3D->nxjp = half * (xfn(i ,jp,k,2) * area(i ,jp,k,2) + xfn(im,jp,k,2) * area(im,jp,k,2));
    greenGeo3D->nyjp = half * (yfn(i ,jp,k,2) * area(i ,jp,k,2) + yfn(im,jp,k,2) * area(im,jp,k,2));
    greenGeo3D->nzjp = half * (zfn(i ,jp,k,2) * area(i ,jp,k,2) + zfn(im,jp,k,2) * area(im,jp,k,2));
    greenGeo3D->nxjm = half * (xfn(i ,jm,k,2) * area(i ,jm,k,2) + xfn(im,jm,k,2) * area(im,jm,k,2));
    greenGeo3D->nyjm = half * (yfn(i ,jm,k,2) * area(i ,jm,k,2) + yfn(im,jm,k,2) * area(im,jm,k,2));
    greenGeo3D->nzjm = half * (zfn(i ,jm,k,2) * area(i ,jm,k,2) + zfn(im,jm,k,2) * area(im,jm,k,2));

    if (ndim == THREE_D)
    {
        greenGeo3D->nxkp = half * (xfn(i ,j,kp,3) * area(i ,j,kp,3) + xfn(im,j,kp,3) * area(im,j,kp,3));
        greenGeo3D->nykp = half * (yfn(i ,j,kp,3) * area(i ,j,kp,3) + yfn(im,j,kp,3) * area(im,j,kp,3));
        greenGeo3D->nzkp = half * (zfn(i ,j,kp,3) * area(i ,j,kp,3) + zfn(im,j,kp,3) * area(im,j,kp,3));
        greenGeo3D->nxkm = half * (xfn(i ,j,km,3) * area(i ,j,km,3) + xfn(im,j,km,3) * area(im,j,km,3));
        greenGeo3D->nykm = half * (yfn(i ,j,km,3) * area(i ,j,km,3) + yfn(im,j,km,3) * area(im,j,km,3));
        greenGeo3D->nzkm = half * (zfn(i ,j,km,3) * area(i ,j,km,3) + zfn(im,j,km,3) * area(im,j,km,3));
    }
    else
    {
        greenGeo3D->nxkp = 0.0;
        greenGeo3D->nykp = 0.0;
        greenGeo3D->nzkp = 0.0;
        greenGeo3D->nxkm = 0.0;
        greenGeo3D->nykm = 0.0;
        greenGeo3D->nzkm = 0.0;
    }

    greenGeo3D->ovol = two / (vol(i,j,k) + vol(im,j,k));
}

void Get_SPM_J(int i, int j, int k, int js, int ks, int ndim, GreenGeo3D *greenGeo3D, RDouble4D &xfn, RDouble4D &yfn, RDouble4D &zfn, RDouble4D &area, RDouble3D &vol)
{
    int ip,im,jp,jm,kp,km;

    //ip = i + 1;
    //im = i;
    //jp = j + 1;
    //jm = j - 1;
    //kp = k + 1;
    //km = k;

    ip = i + 1;
    im = i;
    jp = j + js;
    jm = j - js;
    kp = k + ks;
    km = k;

    greenGeo3D->nxip = half * (xfn(ip,j ,k,1) * area(ip,j ,k,1) + xfn(ip,jm,k,1) * area(ip,jm,k,1));
    greenGeo3D->nyip = half * (yfn(ip,j ,k,1) * area(ip,j ,k,1) + yfn(ip,jm,k,1) * area(ip,jm,k,1));
    greenGeo3D->nzip = half * (zfn(ip,j ,k,1) * area(ip,j ,k,1) + zfn(ip,jm,k,1) * area(ip,jm,k,1));
    greenGeo3D->nxim = half * (xfn(im,j ,k,1) * area(im,j ,k,1) + xfn(im,jm,k,1) * area(im,jm,k,1));
    greenGeo3D->nyim = half * (yfn(im,j ,k,1) * area(im,j ,k,1) + yfn(im,jm,k,1) * area(im,jm,k,1));
    greenGeo3D->nzim = half * (zfn(im,j ,k,1) * area(im,j ,k,1) + zfn(im,jm,k,1) * area(im,jm,k,1));
                                                                                        
    greenGeo3D->nxjp = half * (xfn(i ,j ,k,2) * area(i ,j ,k,2) + xfn(i ,jp,k,2) * area(i ,jp,k,2));
    greenGeo3D->nyjp = half * (yfn(i ,j ,k,2) * area(i ,j ,k,2) + yfn(i ,jp,k,2) * area(i ,jp,k,2));
    greenGeo3D->nzjp = half * (zfn(i ,j ,k,2) * area(i ,j ,k,2) + zfn(i ,jp,k,2) * area(i ,jp,k,2));
    greenGeo3D->nxjm = half * (xfn(i ,j ,k,2) * area(i ,j ,k,2) + xfn(i ,jm,k,2) * area(i ,jm,k,2));
    greenGeo3D->nyjm = half * (yfn(i ,j ,k,2) * area(i ,j ,k,2) + yfn(i ,jm,k,2) * area(i ,jm,k,2));
    greenGeo3D->nzjm = half * (zfn(i ,j ,k,2) * area(i ,j ,k,2) + zfn(i ,jm,k,2) * area(i ,jm,k,2));

    if (ndim == THREE_D)
    {
        greenGeo3D->nxkp = half * (xfn(i,j ,kp,3) * area(i,j ,kp,3) + xfn(i,jm,kp,3) * area(i,jm,kp,3));
        greenGeo3D->nykp = half * (yfn(i,j ,kp,3) * area(i,j ,kp,3) + yfn(i,jm,kp,3) * area(i,jm,kp,3));
        greenGeo3D->nzkp = half * (zfn(i,j ,kp,3) * area(i,j ,kp,3) + zfn(i,jm,kp,3) * area(i,jm,kp,3));
        greenGeo3D->nxkm = half * (xfn(i,j ,km,3) * area(i,j ,km,3) + xfn(i,jm,km,3) * area(i,jm,km,3));
        greenGeo3D->nykm = half * (yfn(i,j ,km,3) * area(i,j ,km,3) + yfn(i,jm,km,3) * area(i,jm,km,3));
        greenGeo3D->nzkm = half * (zfn(i,j ,km,3) * area(i,j ,km,3) + zfn(i,jm,km,3) * area(i,jm,km,3));
    }
    else
    {
        greenGeo3D->nxkp = 0.0;
        greenGeo3D->nykp = 0.0;
        greenGeo3D->nzkp = 0.0;
        greenGeo3D->nxkm = 0.0;
        greenGeo3D->nykm = 0.0;
        greenGeo3D->nzkm = 0.0;
    }

    greenGeo3D->ovol = two / (vol(i,j,k) + vol(i,jm,k));
}

void Get_SPM_K(int i, int j, int k, int js, int ks, GreenGeo3D *greenGeo3D, RDouble4D &xfn, RDouble4D &yfn, RDouble4D &zfn, RDouble4D &area, RDouble3D &vol)
{
    int ip,im,jp,jm,kp,km;

    //ip = i + 1;
    //im = i;
    //jp = j + 1;
    //jm = j;
    //kp = k + 1;
    //km = k - 1;

    ip = i + 1;
    im = i;
    jp = j + js;
    jm = j;
    kp = k + ks;
    km = k - ks;

    greenGeo3D->nxip = half * (xfn(ip,j,k ,1) * area(ip,j,k ,1) + xfn(ip,j,km,1) * area(ip,j,km,1));
    greenGeo3D->nyip = half * (yfn(ip,j,k ,1) * area(ip,j,k ,1) + yfn(ip,j,km,1) * area(ip,j,km,1));
    greenGeo3D->nzip = half * (zfn(ip,j,k ,1) * area(ip,j,k ,1) + zfn(ip,j,km,1) * area(ip,j,km,1));
    greenGeo3D->nxim = half * (xfn(im,j,k ,1) * area(im,j,k ,1) + xfn(im,j,km,1) * area(im,j,km,1));
    greenGeo3D->nyim = half * (yfn(im,j,k ,1) * area(im,j,k ,1) + yfn(im,j,km,1) * area(im,j,km,1));
    greenGeo3D->nzim = half * (zfn(im,j,k ,1) * area(im,j,k ,1) + zfn(im,j,km,1) * area(im,j,km,1));
                                                                                        
    greenGeo3D->nxjp = half * (xfn(i,jp,k ,2) * area(i,jp,k ,2) + xfn(i,jp,km,2) * area(i,jp,km,2));
    greenGeo3D->nyjp = half * (yfn(i,jp,k ,2) * area(i,jp,k ,2) + yfn(i,jp,km,2) * area(i,jp,km,2));
    greenGeo3D->nzjp = half * (zfn(i,jp,k ,2) * area(i,jp,k ,2) + zfn(i,jp,km,2) * area(i,jp,km,2));
    greenGeo3D->nxjm = half * (xfn(i,jm,k ,2) * area(i,jm,k ,2) + xfn(i,jm,km,2) * area(i,jm,km,2));
    greenGeo3D->nyjm = half * (yfn(i,jm,k ,2) * area(i,jm,k ,2) + yfn(i,jm,km,2) * area(i,jm,km,2));
    greenGeo3D->nzjm = half * (zfn(i,jm,k ,2) * area(i,jm,k ,2) + zfn(i,jm,km,2) * area(i,jm,km,2));
                                                                                        
    greenGeo3D->nxkp = half * (xfn(i,j ,k ,3) * area(i,j ,k ,3) + xfn(i,j ,kp,3) * area(i,j ,kp,3));
    greenGeo3D->nykp = half * (yfn(i,j ,k ,3) * area(i,j ,k ,3) + yfn(i,j ,kp,3) * area(i,j ,kp,3));
    greenGeo3D->nzkp = half * (zfn(i,j ,k ,3) * area(i,j ,k ,3) + zfn(i,j ,kp,3) * area(i,j ,kp,3));
    greenGeo3D->nxkm = half * (xfn(i,j ,k ,3) * area(i,j ,k ,3) + xfn(i,j ,km,3) * area(i,j ,km,3));
    greenGeo3D->nykm = half * (yfn(i,j ,k ,3) * area(i,j ,k ,3) + yfn(i,j ,km,3) * area(i,j ,km,3));
    greenGeo3D->nzkm = half * (zfn(i,j ,k ,3) * area(i,j ,k ,3) + zfn(i,j ,km,3) * area(i,j ,km,3));

    greenGeo3D->ovol = two / (vol(i,j,k) + vol(i,j,km));
}

void DXDYDZ(RDouble3D & f,GreenGeo3D * g,int i,int j,int k,int js,int ks,RDouble &dfdx, RDouble &dfdy, RDouble &dfdz)
{
    int    ip,im,jp,jm,kp,km;
    RDouble fip,fim,fjp,fjm,fkp,fkm;

    //ip = i + 1;
    //im = i - 1;
    //jp = j + 1;
    //jm = j - 1;
    //kp = k + 1;
    //km = k - 1;

    ip = i + 1;
    im = i - 1;
    jp = j + js;
    jm = j - js;
    kp = k + ks;
    km = k - ks;

    fip = half * (f(i,j,k) + f(ip,j,k));
    fim = half * (f(i,j,k) + f(im,j,k));
    fjp = half * (f(i,j,k) + f(i,jp,k));
    fjm = half * (f(i,j,k) + f(i,jm,k));
    fkp = half * (f(i,j,k) + f(i,j,kp));
    fkm = half * (f(i,j,k) + f(i,j,km));

    dfdx = g->ovol * (fip * g->nxip - fim * g->nxim + fjp * g->nxjp - fjm * g->nxjm + fkp * g->nxkp - fkm * g->nxkm);
    dfdy = g->ovol * (fip * g->nyip - fim * g->nyim + fjp * g->nyjp - fjm * g->nyjm + fkp * g->nykp - fkm * g->nykm);
    dfdz = g->ovol * (fip * g->nzip - fim * g->nzim + fjp * g->nzjp - fjm * g->nzjm + fkp * g->nzkp - fkm * g->nzkm);
}

void DXDYDZ(RDouble4D & f,GreenGeo3D * g,int i,int j,int k,int js,int ks,int npos,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz)
{
    int    ip,im,jp,jm,kp,km;
    RDouble fip,fim,fjp,fjm,fkp,fkm;

    //ip = i + 1;
    //im = i - 1;
    //jp = j + 1;
    //jm = j - 1;
    //kp = k + 1;
    //km = k - 1;

    ip = i + 1;
    im = i - 1;
    jp = j + js;
    jm = j - js;
    kp = k + ks;
    km = k - ks;

    fip = half * (f(i,j,k,npos) + f(ip,j,k,npos));
    fim = half * (f(i,j,k,npos) + f(im,j,k,npos));
    fjp = half * (f(i,j,k,npos) + f(i,jp,k,npos));
    fjm = half * (f(i,j,k,npos) + f(i,jm,k,npos));
    fkp = half * (f(i,j,k,npos) + f(i,j,kp,npos));
    fkm = half * (f(i,j,k,npos) + f(i,j,km,npos));

    dfdx = g->ovol * (fip * g->nxip - fim * g->nxim + fjp * g->nxjp - fjm * g->nxjm + fkp * g->nxkp - fkm * g->nxkm);
    dfdy = g->ovol * (fip * g->nyip - fim * g->nyim + fjp * g->nyjp - fjm * g->nyjm + fkp * g->nykp - fkm * g->nykm);
    dfdz = g->ovol * (fip * g->nzip - fim * g->nzim + fjp * g->nzjp - fjm * g->nzjm + fkp * g->nzkp - fkm * g->nzkm);
}

void DXDYDZ_I(RDouble3D & f,GreenGeo3D * g,int i,int j,int k,int js,int ks,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz)
{
    int im,jm,km,ip,jp,kp;
    RDouble fim,fjm,fkm,fip,fjp,fkp;

    //ip = i;
    //im = i - 1;
    //jp = j + 1;
    //jm = j - 1;
    //kp = k + 1;
    //km = k - 1;

    ip = i;
    im = i - 1;
    jp = j + js;
    jm = j - js;
    kp = k + ks;
    km = k - ks;

    fip = f(ip,j,k);
    fim = f(im,j,k);
    fjp = fourth * (fip + fim + f(ip,jp,k) + f(im,jp,k));
    fjm = fourth * (fip + fim + f(ip,jm,k) + f(im,jm,k));
    fkp = fourth * (fip + fim + f(ip,j,kp) + f(im,j,kp));
    fkm = fourth * (fip + fim + f(ip,j,km) + f(im,j,km));

    dfdx = g->ovol * (fip * g->nxip - fim * g->nxim + fjp * g->nxjp - fjm * g->nxjm + fkp * g->nxkp - fkm * g->nxkm);
    dfdy = g->ovol * (fip * g->nyip - fim * g->nyim + fjp * g->nyjp - fjm * g->nyjm + fkp * g->nykp - fkm * g->nykm);
    dfdz = g->ovol * (fip * g->nzip - fim * g->nzim + fjp * g->nzjp - fjm * g->nzjm + fkp * g->nzkp - fkm * g->nzkm);
}


void DXDYDZ_J(RDouble3D & f,GreenGeo3D * g,int i,int j,int k,int js,int ks,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz)
{
    int im,jm,km,ip,jp,kp;
    RDouble fim,fjm,fkm,fip,fjp,fkp;

    //ip = i + 1;
    //im = i - 1;
    //jp = j;
    //jm = j - 1;
    //kp = k + 1;
    //km = k - 1;

    ip = i + 1;
    im = i - 1;
    jp = j;
    jm = j - js;
    kp = k + ks;
    km = k - ks;

    fjp = f(i,jp,k);
    fjm = f(i,jm,k);
    fkp = fourth * (fjp + fjm + f(i,jp,kp) + f(i,jm,kp));
    fkm = fourth * (fjp + fjm + f(i,jp,km) + f(i,jm,km));
    fip = fourth * (fjp + fjm + f(ip,jp,k) + f(ip,jm,k));
    fim = fourth * (fjp + fjm + f(im,jp,k) + f(im,jm,k));

    dfdx = g->ovol * (fip * g->nxip - fim * g->nxim + fjp * g->nxjp - fjm * g->nxjm + fkp * g->nxkp - fkm * g->nxkm);
    dfdy = g->ovol * (fip * g->nyip - fim * g->nyim + fjp * g->nyjp - fjm * g->nyjm + fkp * g->nykp - fkm * g->nykm);
    dfdz = g->ovol * (fip * g->nzip - fim * g->nzim + fjp * g->nzjp - fjm * g->nzjm + fkp * g->nzkp - fkm * g->nzkm);
}

void DXDYDZ_K(RDouble3D & f,GreenGeo3D * g,int i,int j,int k,int js,int ks,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz)
{
    int im,jm,km,ip,jp,kp;
    RDouble fim,fjm,fkm,fip,fjp,fkp;

    //ip = i + 1;
    //im = i - 1;
    //jp = j + 1;
    //jm = j - 1;
    //kp = k;
    //km = k - 1;

    ip = i + 1;
    im = i - 1;
    jp = j + js;
    jm = j - js;
    kp = k;
    km = k - ks;

    fkp = f(i,j,kp);
    fkm = f(i,j,km);
    fip = fourth * (fkp + fkm + f(ip,j,kp) + f(ip,j,km));
    fim = fourth * (fkp + fkm + f(im,j,kp) + f(im,j,km));
    fjp = fourth * (fkp + fkm + f(i,jp,kp) + f(i,jp,km));
    fjm = fourth * (fkp + fkm + f(i,jm,kp) + f(i,jm,km));

    dfdx = g->ovol * (fip * g->nxip - fim * g->nxim + fjp * g->nxjp - fjm * g->nxjm + fkp * g->nxkp - fkm * g->nxkm);
    dfdy = g->ovol * (fip * g->nyip - fim * g->nyim + fjp * g->nyjp - fjm * g->nyjm + fkp * g->nykp - fkm * g->nykm);
    dfdz = g->ovol * (fip * g->nzip - fim * g->nzim + fjp * g->nzjp - fjm * g->nzjm + fkp * g->nzkp - fkm * g->nzkm);
}

void DXDYDZ_I(RDouble4D & f,GreenGeo3D * g,int i,int j,int k,int js,int ks,int npos,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz)
{
    int im,jm,km,ip,jp,kp;
    RDouble fim,fjm,fkm,fip,fjp,fkp;

    //ip = i;
    //im = i - 1;
    //jp = j + 1;
    //jm = j - 1;
    //kp = k + 1;
    //km = k - 1;

    ip = i;
    im = i - 1;
    jp = j + js;
    jm = j - js;
    kp = k + ks;
    km = k - ks;

    fip = f(ip,j,k,npos);
    fim = f(im,j,k,npos);
    fjp = fourth * (fip + fim + f(ip,jp,k,npos) + f(im,jp,k,npos));
    fjm = fourth * (fip + fim + f(ip,jm,k,npos) + f(im,jm,k,npos));
    fkp = fourth * (fip + fim + f(ip,j,kp,npos) + f(im,j,kp,npos));
    fkm = fourth * (fip + fim + f(ip,j,km,npos) + f(im,j,km,npos));

    dfdx = g->ovol * (fip * g->nxip - fim * g->nxim + fjp * g->nxjp - fjm * g->nxjm + fkp * g->nxkp - fkm * g->nxkm);
    dfdy = g->ovol * (fip * g->nyip - fim * g->nyim + fjp * g->nyjp - fjm * g->nyjm + fkp * g->nykp - fkm * g->nykm);
    dfdz = g->ovol * (fip * g->nzip - fim * g->nzim + fjp * g->nzjp - fjm * g->nzjm + fkp * g->nzkp - fkm * g->nzkm);
}

void DXDYDZ_J(RDouble4D & f,GreenGeo3D * g,int i,int j,int k,int js,int ks,int npos,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz)
{
    int im,jm,km,ip,jp,kp;
    RDouble fim,fjm,fkm,fip,fjp,fkp;

    //ip = i + 1;
    //im = i - 1;
    //jp = j;
    //jm = j - 1;
    //kp = k + 1;
    //km = k - 1;

    ip = i + 1;
    im = i - 1;
    jp = j;
    jm = j - js;
    kp = k + ks;
    km = k - ks;

    fjp = f(i,jp,k,npos);
    fjm = f(i,jm,k,npos);
    fkp = fourth * (fjp + fjm + f(i,jp,kp,npos) + f(i,jm,kp,npos));
    fkm = fourth * (fjp + fjm + f(i,jp,km,npos) + f(i,jm,km,npos));
    fip = fourth * (fjp + fjm + f(ip,jp,k,npos) + f(ip,jm,k,npos));
    fim = fourth * (fjp + fjm + f(im,jp,k,npos) + f(im,jm,k,npos));

    dfdx = g->ovol * (fip * g->nxip - fim * g->nxim + fjp * g->nxjp - fjm * g->nxjm + fkp * g->nxkp - fkm * g->nxkm);
    dfdy = g->ovol * (fip * g->nyip - fim * g->nyim + fjp * g->nyjp - fjm * g->nyjm + fkp * g->nykp - fkm * g->nykm);
    dfdz = g->ovol * (fip * g->nzip - fim * g->nzim + fjp * g->nzjp - fjm * g->nzjm + fkp * g->nzkp - fkm * g->nzkm);
}

void DXDYDZ_K(RDouble4D & f,GreenGeo3D * g,int i,int j,int k,int js,int ks,int npos,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz)
{
    int im,jm,km,ip,jp,kp;
    RDouble fim,fjm,fkm,fip,fjp,fkp;

    //ip = i + 1;
    //im = i - 1;
    //jp = j + 1;
    //jm = j - 1;
    //kp = k;
    //km = k - 1;

    ip = i + 1;
    im = i - 1;
    jp = j + js;
    jm = j - js;
    kp = k;
    km = k - ks;

    fkp = f(i,j,kp,npos);
    fkm = f(i,j,km,npos);
    fip = fourth * (fkp + fkm + f(ip,j,kp,npos) + f(ip,j,km,npos));
    fim = fourth * (fkp + fkm + f(im,j,kp,npos) + f(im,j,km,npos));
    fjp = fourth * (fkp + fkm + f(i,jp,kp,npos) + f(i,jp,km,npos));
    fjm = fourth * (fkp + fkm + f(i,jm,kp,npos) + f(i,jm,km,npos));

    dfdx = g->ovol * (fip * g->nxip - fim * g->nxim + fjp * g->nxjp - fjm * g->nxjm + fkp * g->nxkp - fkm * g->nxkm);
    dfdy = g->ovol * (fip * g->nyip - fim * g->nyim + fjp * g->nyjp - fjm * g->nyjm + fkp * g->nykp - fkm * g->nykm);
    dfdz = g->ovol * (fip * g->nzip - fim * g->nzim + fjp * g->nzjp - fjm * g->nzjm + fkp * g->nzkp - fkm * g->nzkm);
}

void GradCenter(Grid *grid_in, RDouble4D &q, RDouble3D &dqdx, RDouble3D &dqdy, RDouble3D &dqdz, int m)
{
    StructGrid *grid = StructGridCast(grid_in);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &xfv  = *(grid->GetFaceVectorX());
    RDouble4D &yfv  = *(grid->GetFaceVectorY());
    RDouble4D &zfv  = *(grid->GetFaceVectorZ());
    RDouble3D &vol  = *(grid->GetCellVolume() );

    int ndim = GetDim();

    dqdx = 0.0;
    dqdy = 0.0;
    dqdz = 0.0;

    for (int nsurf = 1; nsurf <= ndim; ++ nsurf)
    {
        int il1,jl1,kl1;
        GetNsurfIndex(nsurf, il1, jl1, kl1);

        int ist = 1;
        int ied = ni-1+il1;
        int jst = 1;
        int jed = nj-1+jl1;
        int kst = 1;
        int ked = nk-1+kl1;

        if (ndim == TWO_D) ked = 1;

        int il,jl,kl;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    il = i - il1;
                    jl = j - jl1;
                    kl = k - kl1;

                    RDouble phis = q(i,j,k,m) + q(il,jl,kl,m);

                    RDouble ddx = phis * xfv(i,j,k,nsurf);
                    RDouble ddy = phis * yfv(i,j,k,nsurf);
                    RDouble ddz = phis * zfv(i,j,k,nsurf);

                    dqdx(i ,j ,k) -= ddx;
                    dqdy(i ,j ,k) -= ddy;
                    dqdz(i ,j ,k) -= ddz;

                    dqdx(il,jl,kl) += ddx;
                    dqdy(il,jl,kl) += ddy;
                    dqdz(il,jl,kl) += ddz;
                }
            }
        }
    }

    int ist = 1;
    int ied = ni-1;
    int jst = 1;
    int jed = nj-1;
    int kst = 1;
    int ked = nk-1;

    if (ndim == TWO_D) ked = 1;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble oov = half / vol(i,j,k);
                dqdx(i,j,k) *= oov;
                dqdy(i,j,k) *= oov;
                dqdz(i,j,k) *= oov;
            }
        }
    }
}

void GradCenter(Grid *grid_in, RDouble3D &q, RDouble3D &dqdx, RDouble3D &dqdy, RDouble3D &dqdz)
{
    StructGrid *grid = StructGridCast(grid_in);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &xfv  = *(grid->GetFaceVectorX());
    RDouble4D &yfv  = *(grid->GetFaceVectorY());
    RDouble4D &zfv  = *(grid->GetFaceVectorZ());
    RDouble3D &vol  = *(grid->GetCellVolume() );

    int ndim = GetDim();

    dqdx = 0.0;
    dqdy = 0.0;
    dqdz = 0.0;

    for (int nsurf = 1; nsurf <= ndim; ++ nsurf)
    {
        int il1,jl1,kl1;
        GetNsurfIndex(nsurf, il1, jl1, kl1);

        int ist = 1;
        int ied = ni-1+il1;
        int jst = 1;
        int jed = nj-1+jl1;
        int kst = 1;
        int ked = nk-1+kl1;

        if (ndim == TWO_D) ked = 1;

        int il,jl,kl;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    il = i - il1;
                    jl = j - jl1;
                    kl = k - kl1;

                    RDouble phis = q(i,j,k) + q(il,jl,kl);

                    RDouble ddx = phis * xfv(i,j,k,nsurf);
                    RDouble ddy = phis * yfv(i,j,k,nsurf);
                    RDouble ddz = phis * zfv(i,j,k,nsurf);

                    dqdx(i ,j ,k) -= ddx;
                    dqdy(i ,j ,k) -= ddy;
                    dqdz(i ,j ,k) -= ddz;

                    dqdx(il,jl,kl) += ddx;
                    dqdy(il,jl,kl) += ddy;
                    dqdz(il,jl,kl) += ddz;
                }
            }
        }
    }

    int ist = 1;
    int ied = ni-1;
    int jst = 1;
    int jed = nj-1;
    int kst = 1;
    int ked = nk-1;

    if (ndim == TWO_D) ked = 1;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble oov = half / vol(i,j,k);
                dqdx(i,j,k) *= oov;
                dqdy(i,j,k) *= oov;
                dqdz(i,j,k) *= oov;
            }
        }
    }
}

void Get_GEO_I(int i, int j, int k, int im, int ip, int ndim, Geometry3D *geometry3D, RDouble4D &xfv, RDouble4D &yfv, RDouble4D &zfv, RDouble3D &vol)
{

    geometry3D->kx = xfv(i,j,k,1);
    geometry3D->ky = yfv(i,j,k,1);
    geometry3D->kz = zfv(i,j,k,1);
                                                         
    geometry3D->ex = fourth * (xfv(im,j,k,2) + xfv(im,j+1,k,2) + xfv(ip,j,k,2) + xfv(ip,j+1,k,2));
    geometry3D->ey = fourth * (yfv(im,j,k,2) + yfv(im,j+1,k,2) + yfv(ip,j,k,2) + yfv(ip,j+1,k,2));
    geometry3D->ez = fourth * (zfv(im,j,k,2) + zfv(im,j+1,k,2) + zfv(ip,j,k,2) + zfv(ip,j+1,k,2));

    if (ndim == THREE_D)
    {
        geometry3D->cx = fourth * (xfv(im,j,k,3) + xfv(im,j,k+1,3) + xfv(ip,j,k,3) + xfv(ip,j,k+1,3));
        geometry3D->cy = fourth * (yfv(im,j,k,3) + yfv(im,j,k+1,3) + yfv(ip,j,k,3) + yfv(ip,j,k+1,3));
        geometry3D->cz = fourth * (zfv(im,j,k,3) + zfv(im,j,k+1,3) + zfv(ip,j,k,3) + zfv(ip,j,k+1,3));
    }
    else
    {
        geometry3D->cx = 0.0;
        geometry3D->cy = 0.0;
        geometry3D->cz = 0.0;
    }

    geometry3D->vjacob = two / (vol(im,j,k) + vol(ip,j,k));
}

void Get_GEO_J(int i, int j, int k, int jm, int jp, int ndim, Geometry3D *geometry3D, RDouble4D &xfv, RDouble4D &yfv, RDouble4D &zfv, RDouble3D &vol)
{

    geometry3D->kx = fourth * (xfv(i,jm,k,1) + xfv(i+1,jm,k,1) + xfv(i,jp,k,1) + xfv(i+1,jp,k,1));
    geometry3D->ky = fourth * (yfv(i,jm,k,1) + yfv(i+1,jm,k,1) + yfv(i,jp,k,1) + yfv(i+1,jp,k,1));
    geometry3D->kz = fourth * (zfv(i,jm,k,1) + zfv(i+1,jm,k,1) + zfv(i,jp,k,1) + zfv(i+1,jp,k,1));

    geometry3D->ex = xfv(i,j,k,2);
    geometry3D->ey = yfv(i,j,k,2);
    geometry3D->ez = zfv(i,j,k,2);                                                         

    if (ndim == THREE_D)
    {
        geometry3D->cx = fourth * (xfv(i,jm,k,3) + xfv(i,jm,k+1,3) + xfv(i,jp,k,3) + xfv(i,jp,k+1,3));
        geometry3D->cy = fourth * (yfv(i,jm,k,3) + yfv(i,jm,k+1,3) + yfv(i,jp,k,3) + yfv(i,jp,k+1,3));
        geometry3D->cz = fourth * (zfv(i,jm,k,3) + zfv(i,jm,k+1,3) + zfv(i,jp,k,3) + zfv(i,jp,k+1,3));
    }
    else
    {
        geometry3D->cx = 0.0;
        geometry3D->cy = 0.0;
        geometry3D->cz = 0.0;
    }

    geometry3D->vjacob = two / (vol(i,jm,k) + vol(i,jp,k));
}

void Get_GEO_K(int i, int j, int k, int km, int kp, Geometry3D *geometry3D, RDouble4D &xfv, RDouble4D &yfv, RDouble4D &zfv, RDouble3D &vol)
{

    geometry3D->kx = fourth * (xfv(i,j,km,1) + xfv(i+1,j,km,1) + xfv(i,j,kp,1) + xfv(i+1,j,kp,1));
    geometry3D->ky = fourth * (yfv(i,j,km,1) + yfv(i+1,j,km,1) + yfv(i,j,kp,1) + yfv(i+1,j,kp,1));
    geometry3D->kz = fourth * (zfv(i,j,km,1) + zfv(i+1,j,km,1) + zfv(i,j,kp,1) + zfv(i+1,j,kp,1));

    geometry3D->ex = fourth * (xfv(i,j,km,2) + xfv(i,j+1,km,2) + xfv(i,j,kp,2) + xfv(i,j+1,kp,2));
    geometry3D->ey = fourth * (yfv(i,j,km,2) + yfv(i,j+1,km,2) + yfv(i,j,kp,2) + yfv(i,j+1,kp,2));
    geometry3D->ez = fourth * (zfv(i,j,km,2) + zfv(i,j+1,km,2) + zfv(i,j,kp,2) + zfv(i,j+1,kp,2));

    geometry3D->cx = xfv(i,j,k,3);
    geometry3D->cy = yfv(i,j,k,3);
    geometry3D->cz = zfv(i,j,k,3);
                                                         
    geometry3D->vjacob = two / (vol(i,j,km) + vol(i,j,kp));
}

void DXDYDZ_I(RDouble4D & f,Geometry3D * g,int i,int j,int k,int ndim,int npos,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz,int lr)
{
    int im,jm,km,ip,jp,kp;
    RDouble fkc,fet,fct;
    //set ip & im to be the index of first layer and first ghost layer near solidwall
    //(eapecially when the wall is at the right side of one zone).
    ip = max(i,i+lr);
    im = min(i,i+lr);
    jp = j + 1;
    jm = j - 1;
    if (ndim == THREE_D)
    {
        kp = k + 1;
        km = k - 1;
    }
    else if (ndim == TWO_D)
    {
        kp = k;
        km = k;
    }

    fkc = f(ip,j,k,npos) - f(im,j,k,npos);
    fet = 0.0;
    fct = 0.0;

    dfdx = g->vjacob * (fkc * g->kx + fet * g->ex + fct * g->cx);
    dfdy = g->vjacob * (fkc * g->ky + fet * g->ey + fct * g->cy);
    dfdz = 0.0;
    if (ndim == THREE_D)
    {
        dfdz = g->vjacob * (fkc * g->kz + fet * g->ez + fct * g->cz);
    }
}

void DXDYDZ_J(RDouble4D & f,Geometry3D * g,int i,int j,int k,int ndim,int npos,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz,int lr)
{
    int im,jm,km,ip,jp,kp;
    RDouble fkc,fet,fct;

    ip = i + 1;
    im = i - 1;
    //set jp & jm to be the index of first layer and first ghost layer near solidwall
    //(eapecially when the wall is at the right side of one zone).
    jp = max(j,j+lr);
    jm = min(j,j+lr);
    if (ndim == THREE_D)
    {
        kp = k + 1;
        km = k - 1;
    }
    else if (ndim == TWO_D)
    {
        kp = k;
        km = k;
    }

    fkc = 0.0;
    fet = f(i,jp,k,npos) - f(i,jm,k,npos);
    fct = 0.0;

    dfdx = g->vjacob * (fkc * g->kx + fet * g->ex + fct * g->cx);
    dfdy = g->vjacob * (fkc * g->ky + fet * g->ey + fct * g->cy);
    dfdz = 0.0;
    if (ndim == THREE_D)
    {
        dfdz = g->vjacob * (fkc * g->kz + fet * g->ez + fct * g->cz);
    }
}


void DXDYDZ_K(RDouble4D & f,Geometry3D * g,int i,int j,int k,int npos,RDouble &dfdx,RDouble &dfdy,RDouble &dfdz,int lr)
{
    int im,jm,km,ip,jp,kp;
    RDouble fkc,fet,fct;

    ip = i + 1;
    im = i - 1;
    jp = j + 1;
    jm = j - 1;
    //set kp & km to be the index of first layer and first ghost layer near solidwall
    //(eapecially when the wall is at the right side of one zone).
    kp = max(k,k+lr);
    km = min(k,k+lr);

    fkc = 0.0;
    fet = 0.0;
    fct = f(i,j,kp,npos) - f(i,j,km,npos);

    dfdx = g->vjacob * (fkc * g->kx + fet * g->ex + fct * g->cx);
    dfdy = g->vjacob * (fkc * g->ky + fet * g->ey + fct * g->cy);
    dfdz = g->vjacob * (fkc * g->kz + fet * g->ez + fct * g->cz);
}

}