#include "stdafx.h"
#include "DGMulDimData.h"


#define RINT(A) floor((A)+0.5)

void convxy(double out[], int xdim, int ydim, double filtx[], double filty[], 
	int fxdim, int fydim, int xoff, int yoff,double buff[])
{
	int x, y, k;
	for (y = 0; y < ydim; y++)
	{
		for (x = 0; x < xdim; x++)
		{
			buff[x] = out[x + y * xdim];
			if (!isfinite(buff[x]))
				buff[x] = 0.0;
		}
		for (x = 0; x < xdim; x++)
		{
			double sum1 = 0.0;
			int fstart, fend;
			fstart = ((x - xoff >= xdim) ? x - xdim - xoff + 1 : 0);
			fend = ((x - (xoff + fxdim) < 0) ? x - xoff + 1 : fxdim);

			for (k = fstart; k < fend; k++)
				sum1 += buff[x - xoff - k] * filtx[k];
			out[x + y * xdim] = sum1;
		}
	}
	for (x = 0; x < xdim; x++)
	{
		for (y = 0; y < ydim; y++)
			buff[y] = out[x + y * xdim];

		for (y = 0; y < ydim; y++)
		{
			double sum1 = 0.0;
			int fstart, fend;
			fstart = ((y - yoff >= ydim) ? y - ydim - yoff + 1 : 0);
			fend = ((y - (yoff + fydim) < 0) ? y - yoff + 1 : fydim);

			for (k = fstart; k < fend; k++)
				sum1 += buff[y - yoff - k] * filty[k];
			out[y*xdim + x] = sum1;
		}
	}
}

int Conv3(CDGMulDimData &vol, double filtx[], double filty[], double filtz[],
	int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff, DGDataType dtype,
	CDGMulDimData &oVol)
	//void *oVol,/* mxArray *wplane_args[3], */DGDataType dtype)
{
	double *tmp, *buff, **sortedv;
	int xy, z, k, fstart, fend, startz, endz;
	static double mat[] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
	int xdim, ydim, zdim;

	xdim = vol.m_dim[0];
	ydim = vol.m_dim[1];
	zdim = vol.m_dim[2];

	//oVol = CDGMulDimData(vol.m_dim, dtype);
	oVol = CDGMulDimData(vol.m_dim, dtype,vol.m_Acc.m_DataMajor);
		 
	int slicesize = xdim*ydim;
	tmp = (double *)malloc(xdim*ydim*fzdim* sizeof(double));
	buff = (double *)malloc(ydim > xdim ? ydim * sizeof(double) : xdim*sizeof(double));
	sortedv = (double **)malloc(fzdim * sizeof(double *));


	CDGMulDimData dvol;
	vol.ToDouble(dvol);

	startz = ((fzdim + zoff - 1 < 0) ? fzdim + zoff - 1 : 0);
	endz = zdim + fzdim + zoff - 1;

	for (z = startz; z < endz; z++)
	{
		double sum2 = 0.0;

		if (z >= 0 && z < zdim)
		{
			mat[14] = z + 1.0;
			//slice(mat, tmp + ((z%fzdim)*xdim*ydim), vol->dim[0], vol->dim[1], vol, 0, 0);
			double * tarptr = tmp + (z%fzdim)*xdim*ydim;
			memcpy(tarptr, dvol.SelMatByDimXY({ z }), slicesize*sizeof(double));
			convxy(tarptr, xdim, ydim,
				filtx, filty, fxdim, fydim, xoff, yoff, buff);
		}
		if (z - fzdim - zoff + 1 >= 0 && z - fzdim - zoff + 1 < zdim)
		{
			fstart = ((z >= zdim) ? z - zdim + 1 : 0);
			fend = ((z - fzdim < 0) ? z + 1 : fzdim);

			for (k = 0; k < fzdim; k++)
			{
				int z1 = (((z - k) % fzdim) + fzdim) % fzdim;
				sortedv[k] = &(tmp[z1*xdim*ydim]);
			}

			for (k = fstart, sum2 = 0.0; k < fend; k++)
				sum2 += filtz[k];

			if (/*oVol.bEmpty() ||*/dtype == DGDataType_Double)
			{
				double *obuf;
				/*if (!oVol)
					obuf = mxGetPr(wplane_args[1]);
				else*/
				obuf = &oVol.Ptr<double>()[(z - fzdim - zoff + 1)*ydim*xdim];
				if (sum2)
				{
					for (xy = 0; xy < xdim*ydim; xy++)
					{
						double sum1 = 0.0;
						for (k = fstart; k < fend; k++)
							sum1 += filtz[k] * sortedv[k][xy];

						obuf[xy] = sum1 / sum2;
					}
				}
				else
					for (xy = 0; xy < xdim*ydim; xy++)
						obuf[xy] = 0.0;

			}
			else
			{
				double tmp;
				if (dtype == DGDataType_Float)
				{
					float *obuf;
					obuf = (float *)oVol.ptr();
					obuf = &obuf[(z - fzdim - zoff + 1)*ydim*xdim];
					if (sum2)
					{
						for (xy = 0; xy < xdim*ydim; xy++)
						{
							double sum1 = 0.0;
							for (k = fstart; k < fend; k++)
								sum1 += filtz[k] * sortedv[k][xy];

							obuf[xy] = (float)(sum1 / sum2);
						}
					}
					else
						for (xy = 0; xy < xdim*ydim; xy++)
							obuf[xy] = 0.0;
				}
				else if (dtype == DGDataType_Int)
				{
					int *obuf;
					obuf = (int *)oVol.ptr();
					obuf = &obuf[(z - fzdim - zoff + 1)*ydim*xdim];
					if (sum2)
					{
						for (xy = 0; xy < xdim*ydim; xy++)
						{
							double sum1 = 0.0;
							for (k = fstart; k < fend; k++)
								sum1 += filtz[k] * sortedv[k][xy];
							tmp = sum1 / sum2;
							if (tmp < -2147483648.0) tmp = -2147483648.0;
							else if (tmp > 2147483647.0) tmp = 2147483647.0;
							obuf[xy] = RINT(tmp);
						}
					}
					else
						for (xy = 0; xy < xdim*ydim; xy++)
							obuf[xy] = 0;
				}
				else if (dtype == DGDataType_UnSignedInt)
				{
					unsigned int *obuf;
					obuf = (unsigned int *)oVol.ptr();
					obuf = &obuf[(z - fzdim - zoff + 1)*ydim*xdim];
					if (sum2)
					{
						for (xy = 0; xy < xdim*ydim; xy++)
						{
							double sum1 = 0.0;
							for (k = fstart; k < fend; k++)
								sum1 += filtz[k] * sortedv[k][xy];
							tmp = sum1 / sum2;
							if (tmp < 0) tmp = 0.0;
							else if (tmp > 4294967295.0) tmp = 4294967295.0;
							obuf[xy] = RINT(tmp);
						}
					}
					else
						for (xy = 0; xy < xdim*ydim; xy++)
							obuf[xy] = 0;
				}
				else if (dtype == DGDataType_SignedShort)
				{
					double tmp;
					short *obuf;
					obuf = (short *)oVol.ptr();
					obuf = &obuf[(z - fzdim - zoff + 1)*ydim*xdim];
					if (sum2)
					{
						for (xy = 0; xy < xdim*ydim; xy++)
						{
							double sum1 = 0.0;
							for (k = fstart; k < fend; k++)
								sum1 += filtz[k] * sortedv[k][xy];
							tmp = sum1 / sum2;
							if (tmp < -32768) tmp = -32768;
							else if (tmp > 32767) tmp = 32767;
							obuf[xy] = RINT(tmp);
						}
					}
					else
						for (xy = 0; xy < xdim*ydim; xy++)
							obuf[xy] = 0;
				}
				else if (dtype == DGDataType_UnSignedShort)
				{
					unsigned short *obuf;
					obuf = (unsigned short *)oVol.ptr();
					obuf = &obuf[(z - fzdim - zoff + 1)*ydim*xdim];
					if (sum2)
					{
						for (xy = 0; xy < xdim*ydim; xy++)
						{
							double sum1 = 0.0;
							for (k = fstart; k < fend; k++)
								sum1 += filtz[k] * sortedv[k][xy];
							tmp = sum1 / sum2;
							if (tmp < 0) tmp = 0;
							else if (tmp > 65535) tmp = 65535;
							obuf[xy] = RINT(tmp);
						}
					}
					else
						for (xy = 0; xy < xdim*ydim; xy++)
							obuf[xy] = 0;
				}
				//else if (dtype == DGDataType_SignedByte)
				else if (dtype == DGDataType_Char)
				{
					char *obuf;
					obuf = (char *)oVol.ptr();
					obuf = &obuf[(z - fzdim - zoff + 1)*ydim*xdim];
					if (sum2)
					{
						for (xy = 0; xy < xdim*ydim; xy++)
						{
							double sum1 = 0.0;
							for (k = fstart; k < fend; k++)
								sum1 += filtz[k] * sortedv[k][xy];
							tmp = sum1 / sum2;
							if (tmp < -128) tmp = -128;
							else if (tmp > 127) tmp = 127;
							obuf[xy] = RINT(tmp);
						}
					}
					else
						for (xy = 0; xy < xdim*ydim; xy++)
							obuf[xy] = 0;
				}
				else if (dtype == DGDataType_UnSignedChar)
				{
					unsigned char *obuf;
					obuf = (unsigned char *)oVol.ptr();
					obuf = &obuf[(z - fzdim - zoff + 1)*ydim*xdim];
					if (sum2)
					{
						for (xy = 0; xy < xdim*ydim; xy++)
						{
							double sum1 = 0.0;
							for (k = fstart; k < fend; k++)
								sum1 += filtz[k] * sortedv[k][xy];
							tmp = sum1 / sum2;
							if (tmp < 0) tmp = 0;
							else if (tmp > 255) tmp = 255;
							obuf[xy] = RINT(tmp);
						}
					}
					else
						for (xy = 0; xy < xdim*ydim; xy++)
							obuf[xy] = 0;
				}
				else
				{
					ErrOutput("","Unknown output datatype.");
				}
			}
		}
	}
	free(tmp);
	free(buff);
	free(sortedv);
	return(0);
}


