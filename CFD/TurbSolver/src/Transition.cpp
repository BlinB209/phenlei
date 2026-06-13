#include "Transition.h"

using namespace std;

namespace PHSPACE
{
RDouble ReynoldsNumberBasedOnStrainRate(RDouble density, RDouble distanceToWall, RDouble viscosity,
    RDouble strainrate, RDouble reynoldsNumberInflow)
{
    return density * distanceToWall * distanceToWall * strainrate * reynoldsNumberInflow / viscosity;
};

RDouble ReynoldsNumberBasedOnDissipation(RDouble density, RDouble distanceToWall, RDouble viscosity,
    RDouble dissipationRate, RDouble reynoldsNumberInflow)
{
    return density * distanceToWall * distanceToWall * dissipationRate * reynoldsNumberInflow * reynoldsNumberInflow / viscosity;
};

RDouble ReynoldsNumberBasedOnScf(RDouble heightCrossFlow, RDouble velocity, RDouble density,
    RDouble viscosityLaminar, RDouble viscosityTurbulence, RDouble reynoldsNumberInflow)
{
    RDouble heightCF = 4.0e-6;
    RDouble thetat = 1.0;
    RDouble deltHeight  = heightCrossFlow * (1.0 + MIN((viscosityTurbulence / viscosityLaminar), 0.4));
    RDouble dHPos = MAX(0.1066 - deltHeight, 0.0);
    RDouble dHNeg = MAX((deltHeight - 0.1066), 0.0);
    RDouble fDHPos = dHPos * (6200.0 + 50000.0 * dHPos);
    RDouble fDHNeg = 75.0 * tanh(dHNeg * 80.0);

    int maxIter = 30;
    for (int iter = 0; iter < maxIter; ++iter)
    {
        RDouble reynoldsSCFa = reynoldsNumberInflow * thetat * density * (velocity / 0.82) / viscosityLaminar;
        //RDouble reynoldsSCFb = -14.7 * log(heightCF / thetat) + 296.77 + functionDeltah1 - functionDeltah2;
        RDouble reynoldsSCFb = -35.088 * log(heightCF / thetat) + 319.51 + fDHPos - fDHNeg;
        if (ABS((reynoldsSCFa - reynoldsSCFb) / reynoldsSCFa) <= 1.0e-6)
        {
            break;
        }
        //thetat = (viscosityLaminar / (reynoldsNumberInflow * density * (velocity / 0.82)))
        //    * (-14.7 * log(heightCF / thetat) + 296.77 + functionDeltah1 - functionDeltah2);
        thetat = (viscosityLaminar / (reynoldsNumberInflow * density * (velocity / 0.82)))
            * (-35.088 * log(heightCF / thetat) + 319.51 + fDHPos - fDHNeg);
    }
    thetat = reynoldsNumberInflow * thetat * density * (velocity / 0.82) / viscosityLaminar;
    return thetat;
}

RDouble TimeScaleInSourceTerm(RDouble density, RDouble velocity, RDouble viscosityLaminar, RDouble viscosityTurbulence, RDouble reynoldsNumberInflow)
{
    RDouble timeScale = 500.0 * viscosityLaminar / (density * velocity * velocity * reynoldsNumberInflow);
    //timeScale = MIN(timeScale, density*gridLength/(viscosityLaminar + viscosityTurbulence));
    return timeScale;
};

RDouble MomentumThickness(RDouble density, RDouble velocity, RDouble viscosity, RDouble onsetReynoldsOnMomentumThickness,
    RDouble reynoldsNumberInflow)
{
    return onsetReynoldsOnMomentumThickness * viscosity / (density * velocity * reynoldsNumberInflow);
};

RDouble FlengthGivenByLangtry(RDouble Rectabar)
{
    RDouble Flength;
    if (Rectabar < 400.0)
    {
        Flength = 39.8189 - Rectabar * (1.1927e-2 + 1.32567e-4 * Rectabar);
    }
    else if (Rectabar >= 400.0 && Rectabar < 596.0)
    {
        Flength = 2.63404e2 - Rectabar * (1.23939 - Rectabar * (1.94548e-3 - 1.01695e-6 * Rectabar));
    }
    else if (Rectabar >= 596.0 && Rectabar < 1200.)
    {
        Flength = 0.5 - 3.0e-4 * (Rectabar - 596.0);
    }
    else
    {
        Flength = 0.3188;
    }
    return Flength;
};

RDouble FlengthCaliBrated(RDouble Rectabar)
{
    RDouble Flength;
    if (Rectabar < 103.0)
    {
        Flength = 39.0 - Rectabar * (0.01482 + 2.687e-007 * Rectabar * Rectabar);
    }
    else if (Rectabar < 182.0)
    {
        Flength = 40.74 - Rectabar * (0.06538 - Rectabar * (0.0004908 - 1.857e-006 * Rectabar));
    }
    else if (Rectabar < 360.0)
    {
        Flength = 20.09 + Rectabar * (0.2749 - Rectabar * (0.001379 -  1.567e-006 * Rectabar));
    }
    else if (Rectabar < 687.0)
    {
        Flength = 107.7 - Rectabar * (0.4551 - Rectabar * (0.000649 - 3.103e-007 * Rectabar));
    }
    else if (Rectabar < 959.0)
    {
        Flength = 10.73 - Rectabar * (0.0317 - Rectabar * (3.263e-005 - 1.127e-008 * Rectabar));
    }
    else if (Rectabar < 1870.0)
    {
        Flength = 0.8605 - Rectabar * (0.0008168 - Rectabar * (4.233e-007 - 7.545e-011 * Rectabar));
    }
    else
    {
        Flength = 0.32;
    }
    return Flength;
};

RDouble FlengthUnstrCaliBrated(RDouble Rectabar)
{
    RDouble Flength;
    if (Rectabar < 102.7)
    {
        Flength = 39.0 - Rectabar * (0.01137 + 6.02e-007 * Rectabar * Rectabar);
    }
    else if (Rectabar < 170.3)
    {
        Flength = 39.67 - Rectabar * (0.03083 - Rectabar * (0.0001895 - 1.217e-006 * Rectabar));
    }
    else if (Rectabar < 360.0)
    {
        Flength = 27.49 + Rectabar * (0.1838 - Rectabar * (0.001071 -  1.249e-006 * Rectabar));
    }
    else if (Rectabar < 694.0)
    {
        Flength = 98.13 - Rectabar * (0.4049 - Rectabar * (0.0005647 - 2.648e-007 * Rectabar));
    }
    else if (Rectabar < 959.0)
    {
        Flength = 15.2 - Rectabar * (0.04645 - Rectabar * (4.814e-005 - 1.668e-008 * Rectabar));
    }
    else if (Rectabar < 1870.0)
    {
        Flength = 0.5262 - Rectabar * (0.000556 - Rectabar * (2.846e-007 - 5.073e-011 * Rectabar));
    }
    else
    {
        Flength = 0.15;
    }
    return Flength;
};

RDouble HighReynoldsCorrectionOfFlength(RDouble ReynoldsNumberBasedOnDissipation, RDouble Flengthori)
{
    RDouble Rw = ReynoldsNumberBasedOnDissipation / 500.0;
    RDouble Fsublayer = exp(-(6.25 * Rw * Rw));
    return Flengthori - (Flengthori - 40.0) * Fsublayer;
};

RDouble ViscosityRatio(RDouble density, RDouble viscosity, RDouble kineticEnergy, RDouble dissipationRate, RDouble reynoldsInflow)
{
    return density * kineticEnergy / (viscosity * dissipationRate);
};

RDouble TransitionOnsetFunction(RDouble Rev, RDouble Rectac, RDouble RT)
{
    RDouble Fonset1, Fonset2, Fonset3;

    Fonset1 = Rev / (2.193 * Rectac);

    Fonset2 = MIN(MAX(Fonset1, pow(Fonset1, 4.0)), 2.0);

    Fonset3 = MAX(1.0 - RT * RT * RT / 15.625, 0.0);

    return MAX(Fonset2 - Fonset3, 0.0);
};

RDouble TransitionOnsetMomentumThicknessReynolds(RDouble Rectabar)
{
    if (Rectabar > 1870.0)
    {
        return Rectabar - 593.11 - 0.482 * (Rectabar - 1870.0);
    }
    else
    {
        return - 3.96035 + Rectabar * (1.0120656 - Rectabar * (8.6823e-4 - Rectabar * (6.96506e-7 - 1.74105e-10 * Rectabar)));
    }
};

RDouble TransitionOnsetMomentumThicknessReynoldsCaliBrated(RDouble Rectabar)
{
    RDouble rectac = 0.0;
    if (Rectabar < 103.0)
    {
        rectac = 7.0 + Rectabar * (0.875 + 8.04e-007 * Rectabar * Rectabar);
    }
    else if (Rectabar < 181.0)
    {
        rectac = 11.92 + Rectabar * (0.7317 + Rectabar * (0.001391 - 3.698e-006 * Rectabar));
    }
    else if (Rectabar < 685.0)
    {
        rectac = -12.75 + Rectabar * (1.141 - Rectabar * (0.0008683 - 4.631e-007 * Rectabar));
    }
    else if (Rectabar < 959.0)
    {
        rectac = 159.0 + Rectabar * (0.3885 + Rectabar * (0.0002297 - 7.122e-008 * Rectabar));
    }
    else
    {
        rectac = 104.2 + Rectabar * (0.56 + Rectabar * (5.09e-005 - 9.073e-009 * Rectabar));
    }
    return rectac;
};

RDouble TransitionOnsetMomentumThicknessReynoldsUnstrCaliBrated(RDouble Rectabar)
{
    RDouble rectac = 0.0;
    if (Rectabar < 102.7)
    {
        rectac = -7.0 + Rectabar * (1.058 + 3.394e-006 * Rectabar * Rectabar);
    }
    else if (Rectabar < 170.3)
    {
        rectac = -12.45 + Rectabar * (1.217 - Rectabar * (0.001549 - 1.633e-006 * Rectabar));
    }
    else if (Rectabar < 694.0)
    {
        rectac = -7.556 + Rectabar * (1.131 - Rectabar * (0.001043 - 6.43e-007 * Rectabar));
    }
    else if (Rectabar < 959.0)
    {
        rectac = 82.62 + Rectabar * (0.7413 - Rectabar * (0.0004813 - 3.732e-007 * Rectabar));
    }
    else
    {
        rectac = 602.9 - Rectabar * (0.8863 - Rectabar * (0.001216 - 2.167e-007 * Rectabar));
    }
    return rectac;
};

RDouble SeparationCorrectionOfIntermittency(RDouble gmori, RDouble Rev, RDouble Rectac, RDouble RT, RDouble Fctat, RDouble s1, RDouble refMaNumber)
{
    RDouble xsep = 1.3;
    if (refMaNumber <= 4.0)
    {
        xsep += 0.04375 * (4.0 - refMaNumber) * (4.0 - refMaNumber);
    }
    RDouble Freattach = exp(-pow(RT / 20.0, 4.0));
    RDouble gmsep = Fctat * MIN(xsep, s1 * Freattach * MAX(0.0, Rev / (3.235 * Rectac) - 1.0));

    return MAX(gmori, gmsep);
}

RDouble TurbulenceIntensity(RDouble velocity, RDouble kineticEnergy)
{
    return MAX(0.027, 100.0 * sqrt(2.0 * kineticEnergy / 3.0) / velocity);
}

RDouble ControlFunctionFturb(RDouble RT)
{
    return exp(-pow(RT / 4.0, 4.0));
}

void CorrectionOfBlendingFunctionInSST(RDouble density, RDouble distance, RDouble viscosity, RDouble kineticEnergy, RDouble reynoldsInflow, RDouble &F1)
{
    RDouble ry = reynoldsInflow * density * distance * sqrt(kineticEnergy) / viscosity;
    RDouble F3 = exp(-pow(ry / 120.0, 8.0));

    F1 = MAX(F1, F3);
}

RDouble CorrectionOfDestructionInKEquation(RDouble gmeff, RDouble Dk)
{
    return Dk * MIN(1.0, MAX(0.1, gmeff));
}

RDouble CorrectionOfProductionInKEquation(RDouble gmeff, RDouble Pk)
{
    return gmeff * Pk;
}

RDouble EmpiricalCorrelationOfFlamdacta(RDouble Tu, RDouble lamdacta)
{
    if (lamdacta > 0.0)
    {
        return 1.0 + 0.275 * (1.0 - exp(-35.0 * lamdacta)) * exp(-2.0 * Tu);
    }
    else if (lamdacta < 0.0)
    {
        return 1.0 + lamdacta * (12.986 + lamdacta * (123.66 + 405.689 * lamdacta)) * exp(-pow(Tu / 1.5, 1.5));
    }
    else
    {
        return 1.0;
    }
}

RDouble EmpiricalCorrelationOfRectat(RDouble Tu, RDouble Flamdacta)
{
    RDouble Rectat;
    if (Tu > 1.3)
    {
        Rectat = 331.5 * Flamdacta * pow(Tu - 0.5658, -0.671);
    }
    else
    {
        Rectat = Flamdacta * (1173.51 - 589.428 * Tu + 0.2196 / (Tu * Tu));
    }

    return MIN(MAX(20.0, Rectat), 5000);
}

RDouble AccelerationAlongStreamline(RDouble u, RDouble v, RDouble w,
    RDouble dudx, RDouble dudy, RDouble dudz,
    RDouble dvdx, RDouble dvdy, RDouble dvdz,
    RDouble dwdx, RDouble dwdy, RDouble dwdz)
{
    RDouble absu = MAX(SMALL, DISTANCE(u, v, w));

    RDouble uu = u / absu;
    RDouble vv = v / absu;
    RDouble ww = w / absu;

    RDouble dUx = uu * dudx + vv * dvdx + ww * dwdx;
    RDouble dUy = uu * dudy + vv * dvdy + ww * dwdy;
    RDouble dUz = uu * dudz + vv * dvdz + ww * dwdz;

    return uu * dUx + vv * dUy + ww * dUz;
}

void BlendingFunctionOfFctat(RDouble gm, RDouble Rew, RDouble Rectabar, RDouble vorticity, RDouble viscosity, RDouble density,
    RDouble absu, RDouble distance, RDouble reynoldsInflow, RDouble ce2, RDouble &Fctat, RDouble &Fctatcf)
{
    RDouble ctaBL = Rectabar * viscosity / reynoldsInflow / MAX(density * absu, SMALL);

    RDouble deltaBL = 7.5 * ctaBL;

    RDouble delta = 50.0 * vorticity * distance * deltaBL / MAX(absu, SMALL);

    RDouble Fwake = exp(-(Rew / 1.0e5) * (Rew / 1.0e5));

    RDouble tmp1 = Fwake * exp(-pow(distance / delta, 4.0));
    RDouble tmp2 = 1.0 - ((ce2 * gm - 1.0) / (ce2 - 1.0)) * ((ce2 * gm - 1.0) / (ce2 - 1.0));

    Fctat = MIN(1.0, MAX(tmp1, tmp2));
    Fctatcf = MIN(tmp1, 1.0);
};

RDouble PressureGradientFunction(RDouble density, RDouble momentumThickness, RDouble viscosity, RDouble dUds, RDouble reynoldsInflow)
{
    RDouble lamdacta = density * momentumThickness * momentumThickness * dUds * reynoldsInflow / viscosity;

    return lamdacta;
};

}