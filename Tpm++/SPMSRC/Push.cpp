#include "stdafx.h"
/* Rather than sample from an image according to a deformation,
 * it is also possible to push voxels from one image into
 * another according to the inverse of the deformation.
 * Note that the result is a noisy version of a Jacobian "modulated"
 * image.
 */


void push(mwSize dm[], mwSize m, mwSize n, float def[], float pf[], float po[], float so[])
{
	mwSignedIndex ix, iy, iz, ix1, iy1, iz1;
	mwSize   i, j, mm, tmpz, tmpy;
	float *px, *py, *pz;
	double dx1, dx2, dy1, dy2, dz1, dz2;

	px = def;
	py = def + m;
	pz = def + m * 2;
	mm = dm[0] * dm[1] * dm[2];

	for (i = 0; i < m; i++)
	{
		double x, y, z;

		if (isfinite(pf[i]))
		{
			x = px[i] - 1.0; /* Subtract 1 because of MATLAB indexing */
			y = py[i] - 1.0;
			z = pz[i] - 1.0;

			/* Check range and avoid inserting values outside the FOV. */
			if (x >= 0 && x < dm[0] - 1 && y >= 0 && y < dm[1] - 1 && z >= 0 && z < dm[2] - 1)
			{
				/* A faster function fo voxels that are safely inside the FOV */
				mwSize o000, o100, o010, o110, o001, o101, o011, o111;
				float w000, w100, w010, w110, w001, w101, w011, w111;
				ix = (mwSignedIndex)floor(x); dx1 = x - ix; dx2 = 1.0 - dx1;
				iy = (mwSignedIndex)floor(y); dy1 = y - iy; dy2 = 1.0 - dy1;
				iz = (mwSignedIndex)floor(z); dz1 = z - iz; dz2 = 1.0 - dz1;

				/* Weights for trilinear interpolation */
				w000 = dx2 * dy2*dz2;
				w100 = dx1 * dy2*dz2;
				w010 = dx2 * dy1*dz2;
				w110 = dx1 * dy1*dz2;
				w001 = dx2 * dy2*dz1;
				w101 = dx1 * dy2*dz1;
				w011 = dx2 * dy1*dz1;
				w111 = dx1 * dy1*dz1;

				ix1 = ix + 1;
				iy1 = iy + 1;
				iz1 = iz + 1;

				/* Neighbouring voxels used for trilinear interpolation */
				tmpz = dm[1] * iz;
				tmpy = dm[0] * (iy + tmpz);
				o000 = ix + tmpy;
				o100 = ix1 + tmpy;
				tmpy = dm[0] * (iy1 + tmpz);
				o010 = ix + tmpy;
				o110 = ix1 + tmpy;
				tmpz = dm[1] * iz1;
				tmpy = dm[0] * (iy + tmpz);
				o001 = ix + tmpy;
				o101 = ix1 + tmpy;
				tmpy = dm[0] * (iy1 + tmpz);
				o011 = ix + tmpy;
				o111 = ix1 + tmpy;

				for (j = 0; j < n; j++)
				{
					/* Increment the images themselves */
					float *pj = po + mm * j;
					float  f = pf[i + j * m];
					pj[o000] += f * w000;
					pj[o100] += f * w100;
					pj[o010] += f * w010;
					pj[o110] += f * w110;
					pj[o001] += f * w001;
					pj[o101] += f * w101;
					pj[o011] += f * w011;
					pj[o111] += f * w111;
				}

				if (so != (float *)0)
				{
					/* Increment an image containing the number of voxels added */
					so[o000] += w000;
					so[o100] += w100;
					so[o010] += w010;
					so[o110] += w110;
					so[o001] += w001;
					so[o101] += w101;
					so[o011] += w011;
					so[o111] += w111;
				}
			}
			else if ((x >= -1) && (x < dm[0]) && (y >= -1) && (y < dm[1]) && (z >= -1) && (z < dm[2]))
			{
				/* A slower function for voxels at the edge of the field of view */
				mwSize o[8], nn = 0, k;
				float w[8];

				ix = (mwSignedIndex)floor(x); dx1 = x - ix; dx2 = 1.0 - dx1;
				iy = (mwSignedIndex)floor(y); dy1 = y - iy; dy2 = 1.0 - dy1;
				iz = (mwSignedIndex)floor(z); dz1 = z - iz; dz2 = 1.0 - dz1;
				ix1 = ix + 1;
				iy1 = iy + 1;
				iz1 = iz + 1;
				if (iz >= 0)
				{
					tmpz = dm[1] * iz;
					if (iy >= 0)
					{
						tmpy = dm[0] * (iy + tmpz);
						if (ix >= 0)
						{
							o[nn] = ix + tmpy;
							w[nn] = dx2 * dy2*dz2;
							nn++;
						}
						if (ix1 < dm[0])
						{
							o[nn] = ix1 + tmpy;
							w[nn] = dx1 * dy2*dz2;
							nn++;
						}
					}
					if (iy1 < dm[1])
					{
						tmpy = dm[0] * (iy1 + tmpz);
						if (ix >= 0)
						{
							o[nn] = ix + tmpy;
							w[nn] = dx2 * dy1*dz2;
							nn++;
						}
						if (ix1 < dm[0])
						{
							o[nn] = ix1 + tmpy;
							w[nn] = dx1 * dy1*dz2;
							nn++;
						}
					}
				}
				if (iz1 < dm[2])
				{
					tmpz = dm[1] * iz1;
					if (iy >= 0)
					{
						tmpy = dm[0] * (iy + tmpz);
						if (ix >= 0)
						{
							o[nn] = ix + tmpy;
							w[nn] = dx2 * dy2*dz1;
							nn++;
						}
						if (ix1 < dm[0])
						{
							o[nn] = ix1 + tmpy;
							w[nn] = dx1 * dy2*dz1;
							nn++;
						}
					}
					if (iy1 < dm[1])
					{
						tmpy = dm[0] * (iy1 + tmpz);
						if (ix >= 0)
						{
							o[nn] = ix + tmpy;
							w[nn] = dx2 * dy1*dz1;
							nn++;
						}
						if (ix1 < dm[0])
						{
							o[nn] = ix1 + tmpy;
							w[nn] = dx1 * dy1*dz1;
							nn++;
						}
					}
				}
				if (so != (float *)0)
				{
					for (k = 0; k < nn; k++)
						so[o[k]] += w[k];
				}

				for (j = 0; j < n; j++)
				{
					float *pj = po + mm * j;
					float  f = pf[i + j * m];
					for (k = 0; k < nn; k++)
						pj[o[k]] += f * w[k];
				}
			}
		}
	}
}