#include "Post_ForceMoment.h"
#include <cmath>

namespace PHSPACE
{

RDouble Post_ForceMoment::DistPoint2Line(RDouble *point, RDouble *linePoint1, RDouble *linePoint2)
{
    RDouble p, q, r, dx, dy, dz, di, dj, dk, d, distN2L;
    p = linePoint2[0] - linePoint1[0];
    q = linePoint2[1] - linePoint1[1];
    r = linePoint2[2] - linePoint1[2];
    dx = point[0] - linePoint1[0];
    dy = point[1] - linePoint1[1];
    dz = point[2] - linePoint1[2];
    di =  q * dz - dy * r;
    dj = -p * dz + dx * r;
    dk =  p * dy - dx * q;
    d = di * di + dj * dj +dk * dk;
    d = sqrt(d);
    distN2L = d / (sqrt(p * p + q * q + r * r) + TINY);
    return distN2L;
}

void Post_ForceMoment::PlaneNormal(RDouble *point, RDouble *linePoint1, RDouble *linePoint2, RDouble &normalX, RDouble &normalY, RDouble &normalZ)
{
    RDouble r;
    normalX =  (linePoint2[1] - linePoint1[1]) * (point[2] - linePoint1[2]) - (point[1] - linePoint1[1]) * (linePoint2[2] - linePoint1[2]);
    normalY = -(linePoint2[0] - linePoint1[0]) * (point[2] - linePoint1[2]) + (point[0] - linePoint1[0]) * (linePoint2[2] - linePoint1[2]);
    normalZ =  (linePoint2[0] - linePoint1[0]) * (point[1] - linePoint1[1]) - (point[0] - linePoint1[0]) * (linePoint2[1] - linePoint1[1]);
    r = normalX * normalX + normalY * normalY + normalZ * normalZ;
    r = sqrt(r) + TINY;
    normalX /= r;
    normalY /= r;
    normalZ /= r;
}

RDouble Post_ForceMoment::ComputeMoment(RDouble *point, RDouble *linePoint1, RDouble *linePoint2, RDouble *force)
{
    RDouble distN2L = DistPoint2Line(point, linePoint1, linePoint2);
    RDouble normalX, normalY, normalZ;
    PlaneNormal(point, linePoint1, linePoint2, normalX, normalY, normalZ);

    RDouble moment = distN2L * (force[0] * normalX + force[1] * normalY + force[2] * normalZ);
    return moment;
}

}