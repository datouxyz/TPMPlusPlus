#include "stdafx.h"

#include <Eigen/Dense>

#include "TissueProblityMap.h"
#include "EigenMatFunc.h"
#include <Eigen/KroneckerProduct>

int Conv3(CDGMulDimData &vol, double filtx[], double filty[], double filtz[],
	int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff, DGDataType dtype,
	CDGMulDimData &oVol);

void Clean_Gwc(CDGMulDimData &  PBytes,int level)
{	
	level = 1;
	
	CDGMulDimData b;
	PBytes.SelFrameByDim(3, { 1 },b);

	VECTOR_TYPE kx(3); kx<<0.75,1,0.75;
	VECTOR_TYPE ky(3); ky<<0.75,1,0.75;
	VECTOR_TYPE kz(3); kz<<0.75,1,0.75;
	
	VECTOR_TYPE krozyx = Eigen::kroneckerProduct(Eigen::kroneckerProduct(kz, ky).eval(), kx);

	auto sm = pow((krozyx).sum(),1.0/3.0);

	//sm = sum(kron(kron(kz, ky), kx)) ^ (1 / 3);
	kx = kx.array() / sm; ky = ky.array() / sm; kz = kz.array() / sm;
	double th1 = 0.15;
	if(level == 2) th1 = 0.2;
	int framesize = PBytes.m_dim[0] * PBytes.m_dim[1];
	int niter = 32;
	int niter2 = 32;
	double th = 0;
	for (int j = 1;j<=niter ;j++)
	{		
		if(j > 2) th = th1;	else th = 0.6;  // Dilate after two its of erosion
		
		for (int z = 0; z< PBytes.m_dim[2]; z++)
		{
			byte* gp, *wp, *bp;
			gp = (byte*)PBytes.SelFrameByDim(2, { z,0});
			wp = (byte*)PBytes.SelFrameByDim(2, { z,1});
			bp = (byte*)b.SelFrameByDim(2, { z });			
			for (int i = 0; i < framesize; i++)
			{
				double  bval = (int)bp[i]/255.0;
				double  wpval = (int)wp[i];
				double  gpval = (int)gp[i];
				bval = bval > th ? (wpval + gpval):0;
				bp[i] = (byte)(int)round(bval);
			}
		}
		CDGMulDimData POut;
		Conv3(b, kx.data(), ky.data(), kz.data(), 3, 3, 3, -1, -1, -1, PBytes.m_type, POut);
		b = POut;
	}
	
	CDGMulDimData c = b;
	

	if (niter2 > 0)
	{		
		for (int iter = 1; iter <= niter2; iter++)
		{			
			for (int z = 0; z < PBytes.m_dim[2]; z++)
			{
				byte* gp, *wp, *cp ,*bp;
				gp = (byte*)PBytes.SelFrameByDim(2, { z,0 });
				wp = (byte*)PBytes.SelFrameByDim(2, { z,1 });
				cp = (byte*)PBytes.SelFrameByDim(2, { z,2 });
				bp = (byte*)c.SelFrameByDim(2, { z});

				for (int i = 0; i < framesize; i++)
				{
					double  bval = (int)bp[i] / 255.0;
					double  wpval = (int)wp[i];
					double  gpval = (int)gp[i];
					double  cpval = (int)cp[i];
					bval = bval > th ? (wpval + gpval+ cpval) : 0;
					bp[i] = (byte)(int)round(bval);
				}
			}
			CDGMulDimData POut;
			Conv3(c, kx.data(), ky.data(), kz.data(), 3, 3, 3, -1, -1, -1, PBytes.m_type, POut);
			c = POut;
		}
	}
	/*{
		sitk::Image newimg(DimToSize(c.m_dim), sitkUInt8);
		memcpy(newimg.GetBufferAsUInt8(), c.ptr(), c.m_Bytes.size());
		sitk::WriteImage(newimg, FormatedString("e:\\c.nii"));
	}
	*/
	th = 0.05;

	for (int z = 0; z < PBytes.m_dim[2]; z++)
	{
		vector<MAT_TYPE> slices(PBytes.m_dim[3]);
		
		for (int k1 = 0; k1 < PBytes.m_dim[3]; k1++)
		{
			CDGMulDimData tmp;
			PBytes.SelFrameByDim(2, { z,k1 }, tmp);
			tmp.ToTypeInPlace<double>(DGDataType_Double);
			slices[k1] = (MAT_TYPE)Eigen::Map<MAT_TYPE>((double*)tmp.ptr(),tmp.m_dim[0], tmp.m_dim[1]); //MDToMat(tmp);
			slices[k1] = slices[k1].array()*Div255;			
		}
		MAT_TYPE bp;
		{
			CDGMulDimData tmp; b.SelFrameByDim(2, { z }, tmp);
			tmp.ToTypeInPlace<double>(DGDataType_Double); //bp.MulOn<double>(Div255);
			bp = (MAT_TYPE)Eigen::Map<MAT_TYPE>((double*)tmp.ptr(), tmp.m_dim[0], tmp.m_dim[1]); //MDToMat(tmp);
			bp *= Div255;
		}
		 
		bp = (((bp.array() > th).cast<double>())*(slices[0]+slices[1]).array() > th).cast<double>();
		slices[0].array() *= bp.array();
		slices[1].array() *= bp.array();

		if (niter2 > 0)
		{
			MAT_TYPE cp;			
			{
				CDGMulDimData tmp;c.SelFrameByDim(2, { z }, tmp); 
				tmp.ToTypeInPlace<double>(DGDataType_Double);
				cp = (MAT_TYPE)Eigen::Map<MAT_TYPE>((double*)tmp.ptr(), tmp.m_dim[0], tmp.m_dim[1]);//MDToMat(tmp);
				cp *= Div255;
			}
			cp = (((cp.array() > th).cast<double>())*(slices[0] + slices[1]+ slices[2]).array() > th).cast<double>();
			slices[2].array() *= cp.array();
		}
		slices[4].array() += 1e-4;
		MAT_TYPE tot = MAT_TYPE::Constant(bp.rows(),bp.cols(),eps);
		for (int k1 = 0; k1 < PBytes.m_dim[3]; k1++)
		{
			tot += slices[k1];
		}
		for (int k1 = 0; k1 < PBytes.m_dim[3]; k1++)
		{
			//BYTE* pkz = (BYTE*)PBytes.SelFrameByDim(2, { z,k1 });
			MAT_TYPE tmp = (slices[k1].array() / tot.array() * 255).round();
			CDGMulDimData pkz({ (int)tmp.rows(),(int)tmp.cols() },DGDataType_Double,DataMajor_FStyle_ColumnMajor ,tmp.data());
			pkz.ToTypeInPlace<unsigned char>(DGDataType_UnSignedChar);
			
			memcpy(PBytes.SelFrameByDim(2, { z,k1 }), pkz.ptr(), pkz.SizeInBytes());
		}
	}	
}
