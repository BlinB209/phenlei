#include "Flux_RoeEntropyCorrection.h"
#include "Constants.h"
#include "TK_Exit.h"

namespace PHSPACE
{
// GMRESAD
#ifdef USE_GMRESSOLVER
template <typename T>
LIB_EXPORT void Flux_RoeEntropyCorrection(T &eigv1, T &eigv2, T &eigv3, T pressureCoeffL, T pressureCoeffR, 
                                          T absVel, T vn, T cm, int RoeEntropyCorrectionMethod, T acousticEigenLimit, T convectEigenLimit)
{
    //! eigv1,2,3: The three eigenvalue, which are V/V+c/V-c on the face.
    if (RoeEntropyCorrectionMethod == 1)
    {
        T maxEigvalue = MAX(eigv2, eigv3);
        T tmp0 = maxEigvalue * acousticEigenLimit;
        T tmp1 = maxEigvalue * acousticEigenLimit;
        eigv1 = MAX(tmp0, eigv1);
        eigv2 = MAX(tmp1, eigv2);
        eigv3 = MAX(tmp1, eigv3);
    }    
    else if (RoeEntropyCorrectionMethod == 2)
    {
        //! From JCP,1995,Lin H-C
        T k1 = 0.25;
        T k2 = 5.0;
        T k3 = 15.0;

        T pressureSensorOnface = max(pressureCoeffL, pressureCoeffR);
        T limitEigv1 = (abs(vn) + cm) * k3 * pressureSensorOnface;
        T limitEigv2 = (abs(vn) + cm) * (k1 + k2 * pressureSensorOnface);
        T limitEigv3 = limitEigv2;

        eigv1 = Harten(eigv1, limitEigv1);
        eigv2 = Harten(eigv2, limitEigv2);
        eigv3 = Harten(eigv3, limitEigv3);
    }
    else if (RoeEntropyCorrectionMethod == 3)
    {
        //! 3: Harten type.
        T maxEigenValue = abs(vn) + cm;
        T acousticEigenvalueLimit  = acousticEigenLimit * maxEigenValue;
        T convectiveEigenvalueLimit= convectEigenLimit * maxEigenValue;

        eigv1 = Harten(eigv1, convectiveEigenvalueLimit);
        eigv2 = Harten(eigv2, acousticEigenvalueLimit);
        eigv3 = Harten(eigv3, acousticEigenvalueLimit);
    }
    else if (RoeEntropyCorrectionMethod == 4)
    {
        //! 4: Method in Ustar
        T maxEigvalue = max(eigv2, eigv3);

        T ratio = min(one, fabs(vn) / cm);  
        T tmp0 = min(maxEigvalue * acousticEigenLimit * ratio, maxEigvalue * acousticEigenLimit / 100.0);
        T tmp1 = min(maxEigvalue * acousticEigenLimit * ratio, maxEigvalue * acousticEigenLimit / 100.0);

        eigv1 = max(tmp0, eigv1);
        eigv2 = max(tmp1, eigv2);
        eigv3 = max(tmp1, eigv3);
    }
    else if (RoeEntropyCorrectionMethod == 5)
    {
        //! Otero, AIAA 2014-2094.
        T maxEigvalue = max(eigv2, eigv3);
        T delta = maxEigvalue * acousticEigenLimit;

        eigv1 += delta;
        eigv2 += delta;
        eigv3 += delta;
    }
    else if (RoeEntropyCorrectionMethod == 6)
    {
        // for arbitary grid type, based on the shock wave detector and the absolute velocity V + sound speed C
        T pressureSensorOnface = min(pressureCoeffL, pressureCoeffR);
        T limitEig = (one - pressureSensorOnface) * (absVel + cm) + 0.0001 * cm;

        eigv1 = Harten(eigv1, limitEig);
        eigv2 = Harten(eigv2, limitEig);
        eigv3 = Harten(eigv3, limitEig);
    }
    else
    {
        TK_Exit::ExceptionExit("Error, no entropy correction method is adopt in Roe scheme! \n");
    }
}
template LIB_EXPORT void Flux_RoeEntropyCorrection<RDouble>(RDouble &eigv1, RDouble &eigv2, RDouble &eigv3, RDouble pressureCoeffL, RDouble pressureCoeffR,
                                                            RDouble absVel, RDouble vn, RDouble cm, int RoeEntropyCorrectionMethod, RDouble acousticEigenLimit, RDouble convectEigenLimit);
template LIB_EXPORT void Flux_RoeEntropyCorrection<ADReal>(ADReal &eigv1, ADReal &eigv2, ADReal &eigv3, ADReal pressureCoeffL, ADReal pressureCoeffR,
                                                           ADReal absVel, ADReal vn, ADReal cm, int RoeEntropyCorrectionMethod, ADReal acousticEigenLimit, ADReal convectEigenLimit);    
#else
LIB_EXPORT void Flux_RoeEntropyCorrection(RDouble &eigv1, RDouble &eigv2, RDouble &eigv3, RDouble pressureCoeffL, RDouble pressureCoeffR, 
    RDouble absVel, RDouble vn, RDouble cm, int RoeEntropyCorrectionMethod, RDouble acousticEigenLimit, RDouble convectEigenLimit)
{
    //! eigv1,2,3: The three eigenvalue, which are V/V+c/V-c on the face.
    if (RoeEntropyCorrectionMethod == 1)
    {
        RDouble maxEigvalue = MAX(eigv2, eigv3);
        RDouble tmp0 = maxEigvalue * acousticEigenLimit;
        RDouble tmp1 = maxEigvalue * acousticEigenLimit;
        eigv1 = MAX(tmp0, eigv1);
        eigv2 = MAX(tmp1, eigv2);
        eigv3 = MAX(tmp1, eigv3);
    }    
    else if (RoeEntropyCorrectionMethod == 2)
    {
        //! From JCP,1995,Lin H-C
        RDouble k1 = 0.25;
        RDouble k2 = 5.0;
        RDouble k3 = 15.0;

        RDouble pressureSensorOnface = MAX(pressureCoeffL, pressureCoeffR);
        RDouble limitEigv1 = (ABS(vn) + cm) * k3 * pressureSensorOnface;
        RDouble limitEigv2 = (ABS(vn) + cm) * (k1 + k2 * pressureSensorOnface);
        RDouble limitEigv3 = limitEigv2;

        eigv1 = Harten(eigv1, limitEigv1);
        eigv2 = Harten(eigv2, limitEigv2);
        eigv3 = Harten(eigv3, limitEigv3);
    }
    else if (RoeEntropyCorrectionMethod == 3)
    {
        //! 3: Harten type.
        RDouble maxEigenValue = ABS(vn) + cm;
        RDouble acousticEigenvalueLimit  = acousticEigenLimit * maxEigenValue;
        RDouble convectiveEigenvalueLimit= convectEigenLimit * maxEigenValue;

        eigv1 = Harten(eigv1, convectiveEigenvalueLimit);
        eigv2 = Harten(eigv2, acousticEigenvalueLimit);
        eigv3 = Harten(eigv3, acousticEigenvalueLimit);
    }
    else if (RoeEntropyCorrectionMethod == 4)
    {
        //! 4: Method in Ustar
        RDouble maxEigvalue = MAX(eigv2, eigv3);

        RDouble ratio = MIN(one, fabs(vn) / cm);  
        RDouble tmp0 = MAX(maxEigvalue * acousticEigenLimit * ratio, maxEigvalue * acousticEigenLimit / 100.0);
        RDouble tmp1 = MAX(maxEigvalue * acousticEigenLimit * ratio, maxEigvalue * acousticEigenLimit / 100.0);

        eigv1 = MAX(tmp0, eigv1);
        eigv2 = MAX(tmp1, eigv2);
        eigv3 = MAX(tmp1, eigv3);
    }
    else if (RoeEntropyCorrectionMethod == 5)
    {
        //! Otero, AIAA 2014-2094.
        RDouble maxEigvalue = MAX(eigv2, eigv3);
        RDouble delta = maxEigvalue * acousticEigenLimit;

        eigv1 += delta;
        eigv2 += delta;
        eigv3 += delta;
    }
    else if (RoeEntropyCorrectionMethod == 6)
    {
        //! for arbitary grid type, based on the shock wave detector and the absolute velocity V + sound speed C
        RDouble pressureSensorOnface = MIN(pressureCoeffL, pressureCoeffR);
        RDouble limitEig = (one - pressureSensorOnface) * (absVel + cm) + 0.0001 * cm;

        eigv1 = Harten(eigv1, limitEig);
        eigv2 = Harten(eigv2, limitEig);
        eigv3 = Harten(eigv3, limitEig);
    }
    else
    {
        TK_Exit::ExceptionExit("Error, no entropy correction method is adopt in Roe scheme! \n");
    }
}
#endif
}

