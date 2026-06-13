#include "LBMSolverOMP.hpp"

void  LBM::bouzidi(int num_node2, int *node2, double *delta, vector3<double> *force )    //! Bouzidi 插值反弹，二阶空间精度
{   //! post-stream (after stream is performed)
	force[0].x = 0.0; force[0].y = 0.0; force[0].z = 0.0;
	for (int i = 0; i < num_node2; i++)
	{
		int j = node2[i];
        int z = (int)floor((double)(j) / (double)(NX*NY));    //! relocate 
		int y = (int)floor((double)(j % (NX*NY)) / (double)(NX));
		int x = (j) % (NX);
//#pragma omp parallel for
        for (int l = 0; l < direction_size; l++)
		{
			if (delta[l+i*direction_size] >= 0.00001)
			{
				double q = delta[l + i * direction_size];
				int ip = ei[l].x, jp = ei[l].y, kp = ei[l].z;
				int m = reverse_indexes[l];

                if (q < 0.5)
                {
					int n1 = index(x + 1 * ip, y + 1 * jp, z + 1 * kp);
					int n2 = index(x + 2 * ip, y + 2 * jp, z + 2 * kp);

                    if (obst[n2] == 0)
                    {
						pdf[index(x + ip, y + jp, z + kp, l)] = q * (1. + 2.*q)* pdf[index(x, y, z, m)]
							+ (1. - 4.*q*q)* pdf[index(x + ip, y + jp, z + kp, m)] - q * (1. - 2.*q)*
							pdf[index(x + 2 * ip, y + 2 * jp, z + 2 * kp, m)];
					}
					else
                    {   //! the following is seldem used 以下一般用不上
						if (obst[n1] == 0)
						{
							pdf[index(x + ip, y + jp, z + kp, l)] = 2.* q * pdf[index(x, y, z, m)]
								+ (1. - 2.* q) * pdf[index(x + ip, y + jp, z + kp, m)];
						}
						else
						{
							pdf[index(x + ip, y + jp, z + kp, l)] = pdf[index(x, y, z, m)];
						}
					}
				}
				else
				{
					int n2 = index(x + 2 * ip, y + 2 * jp, z + 2 * kp);
					int n3 = index(x + 3 * ip, y + 3 * jp, z + 3 * kp);
                    if (obst[n3] == 0)
                    {
						pdf[index(x + ip, y + jp, z + kp, l)] = 1.0 / (q*(2.*q + 1.0)) *  pdf[index(x, y, z, m)]
							+ (2.*q - 1.) / q * pdf[index(x + 2 * ip, y + 2 * jp, z + 2 * kp, l)] - (2.*q - 1.) /
							(2.*q + 1.)* pdf[index(x + 3 * ip, y + 3 * jp, z + 3 * kp, l)];
					}
                    else    //! the following is seldem used 以下一般用不上
					{
                        if (obst[n2] == 0)
                        {
							pdf[index(x + ip, y + jp, z + kp, l)] = 1. / (2.* q)*pdf[index(x, y, z, m)] +
								(2.*q - 1.) / (2.*q) *pdf[index(x + 2 * ip, y + 2 * jp, z + 2 * kp, l)];
						}
                        else
                        {
                            //! the following is simple bounce back , seldem used
							pdf[index(x + ip, y + jp, z + kp, l)] = pdf[index(x + ip, y + jp, z + kp, m)]; //just bounce - back 
						}
					}
				}
				force[0].x = force[0].x + ei[m].x * (pdf[index(x, y, z, m)] + pdf[index(x + ip, y + jp, z + kp, l)]);
				force[0].y = force[0].y + ei[m].y * (pdf[index(x, y, z, m)] + pdf[index(x + ip, y + jp, z + kp, l)]);
				force[0].z = force[0].z + ei[m].z * (pdf[index(x, y, z, m)] + pdf[index(x + ip, y + jp, z + kp, l)]);
			}
		}
	}
}
