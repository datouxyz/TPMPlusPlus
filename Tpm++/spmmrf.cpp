

#include <math.h>
#define MAXCLASSES 1024
#define mwSize unsigned int

void mrf1(mwSize dm[], unsigned char q[], float p[], float G[], float w[], int code)
{
	mwSize i0, i1, i2, k, m, n;
	float a[MAXCLASSES], e[MAXCLASSES], *p0 = NULL, *p1 = NULL;
	unsigned char *q0 = NULL, *q1 = NULL;
	int it;

	m = dm[0] * dm[1] * dm[2];

	/* Use a red-black scheme, so the updates are for
	   alternating voxels.  Then do another pass to
	   update the other half.
	   A B A B A B
	   B A B A B A
	   A B A B A B
	   B A B A B A

	   Updates involve computing the number of neighbours
	   of each type (stored in vector a), and using the
	   connectivity matrix (G) to update with:
		   q = (p.*exp(G'*a))/sum((p.*exp(G'*a))
	*/
	for (it = 0; it < 2; it++)
	{
		mwSize i2start = it % 2;
		for (i2 = 0; i2 < dm[2]; i2++) /* Inferior -> Superior */
		{
			mwSize i1start = (i2start == (i2 % 2));
			for (i1 = 0; i1 < dm[1]; i1++) /* Posterior -> Anterior */
			{
				mwSize i0start = (i1start == (i1 % 2));
				p1 = p + dm[0] * (i1 + dm[1] * i2);
				q1 = q + dm[0] * (i1 + dm[1] * i2);

				for (i0 = i0start; i0 < dm[0]; i0 += 2) /* Left -> Right */
				{
					float se;
					unsigned char *qq = NULL;

					/* Pointers to current voxel in first volume */
					p0 = p1 + i0;
					q0 = q1 + i0;

					/* Initialise neighbour counts to zero */
					for (k = 0; k < dm[3]; k++) a[k] = 0.0;

					/* Count neighbours of each class */
					if (i2 > 0)       /* Inferior */
					{

						qq = q0 - dm[0] * dm[1];
						for (k = 0; k < dm[3]; k++) a[k] += qq[k*m] * w[2];
					}

					if (i2 < dm[2] - 1) /* Superior */
					{
						qq = q0 + dm[0] * dm[1];
						for (k = 0; k < dm[3]; k++) a[k] += qq[k*m] * w[2];
					}

					if (i1 > 0)       /* Posterior */
					{
						qq = q0 - dm[0];
						for (k = 0; k < dm[3]; k++) a[k] += qq[k*m] * w[1];
					}

					if (i1 < dm[1] - 1) /* Anterior */
					{
						qq = q0 + dm[0];
						for (k = 0; k < dm[3]; k++) a[k] += qq[k*m] * w[1];
					}

					if (i0 > 0)       /* Left */
					{
						qq = q0 - 1;
						for (k = 0; k < dm[3]; k++) a[k] += qq[k*m] * w[0];
					}

					if (i0 < dm[0] - 1) /* Right */
					{
						qq = q0 + 1;
						for (k = 0; k < dm[3]; k++) a[k] += qq[k*m] * w[0];
					}

					/* Responsibility data is uint8, so correct scaling.
					   Note also that data is divided by 6 (the number
					   of neighbours examined). */
					for (k = 0; k < dm[3]; k++)
						a[k] /= (255.0*6.0);

					if (code == 1)
					{
						/* Weights are in the form of a matrix,
						   shared among all voxels. */
						float *g;
						se = 0.0;
						for (k = 0, g = G; k < dm[3]; k++)
						{
							e[k] = 0;
							for (n = 0; n < dm[3]; n++, g++)
								e[k] += (*g)*a[n];
							e[k] = exp((double)e[k])*p0[k*m];
							se += e[k];
						}
					}
					else if (code == 2)
					{
						/* Weights are assumed to be a diagonal matrix,
						   so only the diagonal elements are passed. */
						se = 0.0;
						for (k = 0; k < dm[3]; k++)
						{
							e[k] = exp((double)(G[k] * a[k]))*p0[k*m];
							se += e[k];
						}
					}
					else if (code == 3)
					{
						/* Separate weights for each voxel, in the form of
						   the full matrix (loads of memory). */
						float *g;
						se = 0.0;
						g = G + i0 + dm[0] * (i1 + dm[1] * i2);
						for (k = 0; k < dm[3]; k++)
						{
							e[k] = 0.0;
							for (n = 0; n < dm[3]; n++, g += m)
								e[k] += (*g)*a[n];
							e[k] = exp((double)e[k])*p0[k*m];
							se += e[k];
						}
					}
					else if (code == 4)
					{
						/* Separate weight matrices for each voxel,
						   where the matrices are assumed to be symmetric
						   with zeros on the diagonal. For a 4x4
						   matrix, the elements are ordered as
						   (2,1), (3,1), (4,1), (3,2), (4,2), (4,3).
						 */
						float *g;
						g = G + i0 + dm[0] * (i1 + dm[1] * i2);
						for (k = 0; k < dm[3]; k++) e[k] = 0.0;

						for (k = 0; k < dm[3]; k++)
						{
							for (n = k + 1; n < dm[3]; n++, g += m)
							{
								e[k] += (*g)*a[n];
								e[n] += (*g)*a[k];
							}
						}
						se = 0.0;
						for (k = 0; k < dm[3]; k++)
						{
							e[k] = exp((double)e[k])*p0[k*m];
							se += e[k];
						}
					}
					else if (code == 5)
					{
						/* Separate weight matrices for each voxel,
						   where the matrices are assumed to be symmetric
						   with zeros on the diagonal. For a 4x4
						   matrix, the elements are ordered as
						   (2,1), (3,1), (4,1), (3,2), (4,2), (4,3).

						   The weight matrices are encoded as uint8, and
						   their values need to be scaled by -0.0625 to
						   bring them into a reasonable range.
						 */
						unsigned char *g;
						g = (unsigned char *)G + i0 + dm[0] * (i1 + dm[1] * i2);
						for (k = 0; k < dm[3]; k++) e[k] = 0.0;

						for (k = 0; k < dm[3]; k++)
						{
							for (n = k + 1; n < dm[3]; n++, g += m)
							{
								e[k] += ((float)(*g))*a[n];
								e[n] += ((float)(*g))*a[k];
							}
						}
						se = 0.0;
						for (k = 0; k < dm[3]; k++)
						{
							e[k] = exp(-0.0625*e[k])*p0[k*m];
							se += e[k];
						}
					}


					/* Normalise responsibilities to sum to 1
					   and rescale for saving as uint8 data. */
					se = 255.0 / se;
					for (k = 0; k < dm[3]; k++)
						q0[k*m] = (unsigned char)(e[k] * se + 0.5);

				}
			}
		}
	}
}