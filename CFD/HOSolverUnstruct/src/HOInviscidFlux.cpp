#include "HODefine.h"
#include "Constants.h"
#include <iostream>

using PHSPACE::RDouble;
namespace HOUnstruct
{
void Vanleer_Scheme(int totalGaussPoint, RDouble **ql, RDouble **qr, RDouble *areax, RDouble *areay, RDouble *areaz, RDouble gama, RDouble **flux)
{
    RDouble vn_l, c2_l, c_l, mach_l, h_l;
    RDouble vn_r, c2_r, c_r, mach_r, h_r;
    RDouble gamm1, gamp1, gami, tmp, fmass;

    gamm1 = gama - 1.0;
    gamp1 = gama + 1.0;
    gami = 1.0 / gama;

    for (int iGaussPoint = 0; iGaussPoint < totalGaussPoint; ++ iGaussPoint)
    {

        vn_l = areax[iGaussPoint] * ql[1][iGaussPoint] + areay[iGaussPoint] * ql[2][iGaussPoint] + areaz[iGaussPoint] * ql[3][iGaussPoint];
        c2_l = gama * (ql[4][iGaussPoint]) / ql[0][iGaussPoint];
        c_l    = sqrt(abs(c2_l));
        mach_l = vn_l / c_l;
        h_l = c2_l / gamm1 + 0.5 * (ql[1][iGaussPoint] * ql[1][iGaussPoint] + ql[2][iGaussPoint] * ql[2][iGaussPoint] + ql[3][iGaussPoint] * ql[3][iGaussPoint]);

        if (mach_l > 1.)
        {
            tmp = ql[0][iGaussPoint] * vn_l;
            flux[0][iGaussPoint] = tmp;
            flux[1][iGaussPoint] = (ql[1][iGaussPoint] * tmp + areax[iGaussPoint] * (ql[4][iGaussPoint]));
            flux[2][iGaussPoint] = (ql[2][iGaussPoint] * tmp + areay[iGaussPoint] * (ql[4][iGaussPoint]));
            flux[3][iGaussPoint] = (ql[3][iGaussPoint] * tmp + areaz[iGaussPoint] * (ql[4][iGaussPoint]));
            flux[4][iGaussPoint] = h_l * tmp;
        }
        else if (mach_l < -1.)
        {
            flux[0][iGaussPoint] = 0.;
            flux[1][iGaussPoint] = 0.;
            flux[2][iGaussPoint] = 0.;
            flux[3][iGaussPoint] = 0.;
            flux[4][iGaussPoint] = 0.; 
        }
        else
        {
            fmass = 0.25 * ql[0][iGaussPoint] * c_l * (mach_l + 1.) * (mach_l + 1.);
            tmp = (2. * c_l - vn_l) * gami;
            flux[0][iGaussPoint] = fmass;
            flux[1][iGaussPoint] = fmass * (areax[iGaussPoint] * tmp + ql[1][iGaussPoint]);
            flux[2][iGaussPoint] = fmass * (areay[iGaussPoint] * tmp + ql[2][iGaussPoint]);
            flux[3][iGaussPoint] = fmass * (areaz[iGaussPoint] * tmp + ql[3][iGaussPoint]);
            flux[4][iGaussPoint] = fmass * h_l;
        }

        // Right contribution
        vn_r = areax[iGaussPoint] * qr[1][iGaussPoint] + areay[iGaussPoint] * qr[2][iGaussPoint] + areaz[iGaussPoint] * qr[3][iGaussPoint];
        c2_r = gama * (qr[4][iGaussPoint]) / qr[0][iGaussPoint];
        c_r    = sqrt(abs(c2_r));
        mach_r = vn_r / c_r;
        h_r = c2_r / gamm1
            + 0.5 * (qr[1][iGaussPoint] * qr[1][iGaussPoint] + qr[2][iGaussPoint] * qr[2][iGaussPoint] + qr[3][iGaussPoint] * qr[3][iGaussPoint]);

        if (mach_r < -1.)
        {
            tmp = qr[0][iGaussPoint] * vn_r;
            flux[0][iGaussPoint] += tmp;
            flux[1][iGaussPoint] += (qr[1][iGaussPoint] * tmp + areax[iGaussPoint] * (qr[4][iGaussPoint]));
            flux[2][iGaussPoint] += (qr[2][iGaussPoint] * tmp + areay[iGaussPoint] * (qr[4][iGaussPoint]));
            flux[3][iGaussPoint] += (qr[3][iGaussPoint] * tmp + areaz[iGaussPoint] * (qr[4][iGaussPoint]));
            flux[4][iGaussPoint] += h_r * tmp;
        }
        else if (mach_r <= 1.)
        {
            fmass = -0.25 * qr[0][iGaussPoint] * c_r * (mach_r - 1.) * (mach_r - 1.);
            tmp = -(2. * c_r + vn_r) * gami;
            flux[0][iGaussPoint] += fmass;
            flux[1][iGaussPoint] += fmass * (areax[iGaussPoint] * tmp + qr[1][iGaussPoint]);
            flux[2][iGaussPoint] += fmass * (areay[iGaussPoint] * tmp + qr[2][iGaussPoint]);
            flux[3][iGaussPoint] += fmass * (areaz[iGaussPoint] * tmp + qr[3][iGaussPoint]);
            flux[4][iGaussPoint] += fmass * h_r;
        }

        for (int ieqn = 0; ieqn < 5; ++ieqn)
        if ((flux[ieqn][iGaussPoint] != flux[ieqn][iGaussPoint]) || fabs(flux[ieqn][iGaussPoint]) > 1.0e20)
        {
            cout << "flux " << ieqn << " = " << flux[ieqn][iGaussPoint] << endl;
        }
    }
}

}