#include "AerodynamicsBasicRelations.h"
#include "Constants.h"
#include <cmath>
using namespace std;

namespace PHSPACE
{

RDouble IsentropicRelations::T0OverT(RDouble gama, RDouble mach)
{
    RDouble m2 = mach * mach;
    return one + half * (gama - one) * m2;
}

RDouble IsentropicRelations::P0OverP(RDouble gama, RDouble mach)
{
    RDouble t0t = T0OverT(gama, mach);
    RDouble cc  = gama / (gama - one);

    return pow(t0t,cc);
}

RDouble IsentropicRelations::Rho0OverRho(RDouble gama, RDouble mach)
{
    RDouble t0t = T0OverT(gama, mach);
    RDouble cc  = one / (gama - one);

    return pow(t0t,cc);
}

RDouble ShockRelations::P2OverP1(RDouble gama, RDouble m1)
{
    RDouble m12 = m1 * m1;
    RDouble gamp1 = gama + one;
    RDouble gamm1 = gama - one;

    return two * gama / gamp1 * m12 - gamm1 / gamp1;
}

RDouble ShockRelations::Rho2OverRho1(RDouble gama, RDouble m1)
{
    RDouble m12 = m1 * m1;
    RDouble gamp1 = gama + one;
    RDouble gamm1 = gama - one;

    return gamp1 * m12 / (gamm1 * m12 + two);
}

RDouble ShockRelations::T2OverT1(RDouble gama, RDouble m1)
{
    RDouble m12 = m1 * m1;
    RDouble gamp1 = gama + one;
    RDouble gamm1 = gama - one;

    RDouble c1 = two * gama * m12 - gamm1;
    RDouble c2 = gamm1 * m12 + two;
    RDouble c3 = gamp1 * gamp1 * m12;

    return c1 * c2 / c3;
}

RDouble ShockRelations::MachBehindShock(RDouble gama, RDouble m1)
{
    RDouble m12 = m1 * m1;
    RDouble gamm1 = gama - one;

    RDouble c1 = two / gamm1 + m12;
    RDouble c2 = two * gama / gamm1 * m12 - one;

    return sqrt(c1 / c2);
}

}
