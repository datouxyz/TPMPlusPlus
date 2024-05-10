#include "stdafx.h"
#include <Eigen/Dense>

#include "TissueProblityMap.h"
#include "EigenMatFunc.h"

#include "shoot_boundary.h"
#include "TPMInternalFuncs.h"
int Conv3(CDGMulDimData &vol, double filtx[], double filty[], double filtz[],
	int fxdim, int fydim, int fzdim, int xoff, int yoff, int zoff, DGDataType dtype,
	CDGMulDimData &oVol);
void decimate(CDGMulDimData & dat, VECTOR_TYPE fwhm, CDGMulDimData &dout)
{
	VECTOR_TYPE lim = (fwhm.array() * 2).ceil();	
	VECTOR_TYPE x = SmoothKernel(fwhm(0), -lim(0), lim(0), KernelSmoothMode_linear, eps);
	VECTOR_TYPE y = SmoothKernel(fwhm(1), -lim(1), lim(1), KernelSmoothMode_linear, eps);
	VECTOR_TYPE z = SmoothKernel(fwhm(2), -lim(2), lim(2), KernelSmoothMode_linear, eps);
	x = x.array() / x.array().sum();
	y = y.array() / y.array().sum();
	z = z.array() / z.array().sum();
	int i = floor((x.size() - 1) / 2.0), j = floor((z.size() - 1) / 2.0), k = floor((z.size() - 1) / 2.0);
	Conv3(dat, x.data(), y.data(), z.data(), x.size(), y.size(), z.size(), -i, -j, -k, dat.m_type, dout);
}


void ExtrapolateDef(CDGMulDimData &Y, MAT_TYPE & M)
{
	int nYFrame = Y.m_dim[0] * Y.m_dim[1] * Y.m_dim[2];
	sMask msk(nYFrame);
	float * pY0 = (float*)Y.SelFrameByDim(3, { 0 });
	msk.TraverseMsk((UpdateMaskFunc)[&](int ind, BYTE &msk)->void
	{
		float & f = pY0[ind];
		msk = (byte)std::isfinite(f);
	}
	);
	msk.ReComputeNM();
	if (msk.m_NM < nYFrame)
	{
		VECTOR_TYPE vx = VoxelSize(M);// 
		CDGMulDimData x1, x2, x3;
		Grid3D<VECTOR_TYPE_F, float>(
			BuildVector<VECTOR_TYPE_F, float>(1.0, Y.m_dim[0], 1.0),
			BuildVector<VECTOR_TYPE_F, float>(1.0, Y.m_dim[1], 1.0),
			BuildVector<VECTOR_TYPE_F, float>(1.0, Y.m_dim[2], 1.0),
			x1, x2, x3, DGDataType_Float);
		CDGMulDimData X({ Y.m_dim[0] ,Y.m_dim[1],Y.m_dim[2],3 }, DGDataType_Float,DataMajor_FStyle_ColumnMajor);
		memcpy(X.SelFrameByDim(3, { 0 }), x1.Ptr<float>(), x1.SizeInBytes());
		memcpy(X.SelFrameByDim(3, { 1 }), x2.Ptr<float>(), x2.SizeInBytes());
		memcpy(X.SelFrameByDim(3, { 2 }), x3.Ptr<float>(), x3.SizeInBytes());

		MAT_TYPE M1;
		GetClosestAffine(X, Y, CDGMulDimData(), CDGMulDimData(), M1, NULL);
		X.clear();
		int oldbnd = get_bound();

		set_bound(1);

		for (int d = 0; d < 3; d++)
		{
			CDGMulDimData x;
			{
				CDGMulDimData dx1 = x1, dx2 = x2, dx3 = x3;
				dx1.MulOn<float>(M1(d, 0));
				dx2.MulOn<float>(M1(d, 1));
				dx3.MulOn<float>(M1(d, 2));
				x = dx1;
				x.AddOn<float>(dx2); x.AddOn<float>(dx3);
				x.AddOn<float>((double)M1(d, 3));
			}
			CDGMulDimData u;
			Y.SelFrameByDim(3, { d }, u);
			u.SubOn<float>(x);
			{
				float *pu = u.Ptr<float>();

				msk.TraverseMsk((UpdateMaskFunc)[&](int ind, BYTE &msk)->void
				{ if (!(bool)msk) { pu[ind] = 0; }		});
			}

			{
				VECTOR_TYPE p(8); p << vx, 0, 0.00001, 0, 1, 1;

				CDGMulDimData A({ Y.m_dim[0] , Y.m_dim[1] , Y.m_dim[2] }, DGDataType_Float,Y.m_Acc.m_DataMajor);
				float * pA = A.Ptr<float>();
				msk.TraverseMsk((UpdateMaskFunc)[&](int ind, BYTE &msk)->void
				{ pA[ind] = (float)msk; });

				CDGMulDimData uout;
				fmgN(A, u, p.data(), NULL, uout);
				u = uout;
			}

			CDGMulDimData y;
			Y.SelFrameByDim(3, { d }, y);
			float * py = y.Ptr<float>();
			float * pu = u.Ptr<float>();
			float * px = x.Ptr<float>();
			msk.TraverseMsk((UpdateMaskFunc)[&](int ind, BYTE &msk)->void
			{ if (!(bool)msk) { py[ind] = pu[ind] + px[ind]; }		});

			memcpy(Y.SelFrameByDim(3, { d }), y.ptr(), y.SizeInBytes());
		}
		set_bound(oldbnd);
	}
}
