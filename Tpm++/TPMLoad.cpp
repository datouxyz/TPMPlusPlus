#include "stdafx.h"
#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>
#include "TissueProblityMap.h"
#include "EigenMatFunc.h"
#include "TPMInternalFuncs.h"
void LoadChans(IN sChannels & channels, IN sRegData regdata, OUT ChanVec &chans)
{
	chans = ChanVec(channels.size());

	for (int ichan = 0; ichan < channels.size(); ichan++)
	{
		sChannel &ch = channels[ichan];

		double fwhm = ch.m_Biasfwhm;
		double biasreg = ch.m_Bisareg;
		sitk::Image & img = ch.m_Img;
		MAT_TYPE cm = GetMatFromImageDirOriSpace(img.GetDirection(), img.GetOrigin(), img.GetSpacing(), true, true);
		VECTOR_TYPE vx = VoxelSize(cm);
		auto &chan = chans[ichan];
		SizeVec d0 = img.GetSize();
		double sd = 0;
		DimVec d3(3);

		sd = vx[0] * d0[0] / fwhm; d3[0] = ceil(sd * 2);
		VECTOR_TYPE bv1 = BuildVector<VECTOR_TYPE,double>(0, d3[0] - 1, 1.0);
		VECTOR_TYPE krn_x = (-bv1.array().pow(2) / (sd*sd)).exp() / sqrt(vx(0));

		sd = vx[1] * d0[1] / fwhm; d3[1] = ceil(sd * 2);
		VECTOR_TYPE bv2 = BuildVector<VECTOR_TYPE, double>(0, d3[1] - 1, 1.0);
		VECTOR_TYPE krn_y = (-bv2.array().pow(2) / (sd* sd)).exp() / sqrt(vx(1));

		sd = vx[2] * d0[2] / fwhm; d3[2] = ceil(sd * 2);
		VECTOR_TYPE bv3 = BuildVector<VECTOR_TYPE, double>(0, d3[2] - 1, 1.0);
		VECTOR_TYPE krn_z = (-bv3.array().pow(2) / (sd*sd)).exp() / sqrt(vx(2));
		VECTOR_TYPE Cbias = kroneckerProduct(krn_z, kroneckerProduct(krn_y, krn_x)).array().pow(-2)*biasreg*regdata.m_FudgeFactor;
		chan.C = MAT_SPARSE_TYPE(Cbias.size(), Cbias.size());
		for (int i = 0; i < Cbias.size(); i++)
		{
			chan.C.insert(i, i) = Cbias(i);
		}

		chan.B3 = DCT_MTX(d0[2], d3[2], regdata.x3);
		VECTOR_TYPE v2r = regdata.x2.row(0);
		chan.B2 = DCT_MTX(d0[1], d3[1], v2r);
		VECTOR_TYPE v1r = regdata.x1.col(0);
		chan.B1 = DCT_MTX(d0[0], d3[0], v1r);
		if (regdata.TBias.m_ElCount > 0)
		{
			chan.T = regdata.TBias;
		}
		else
		{
			chan.T = CDGMulDimData({ d3[0],d3[1],d3[2] }, DGDataType_Double,DataMajor_FStyle_ColumnMajor);
		}

	}
}

#include "shoot_regularisers.h"
void Vel2Mom(double _param[8], CDGMulDimData & datain, CDGMulDimData & dataout)
{
	double param[8] =
	{
		1.0/_param[0],
		1.0/_param[1],
		1.0/_param[2],
		_param[3],
		_param[4],
		_param[5],
		_param[6],
		_param[7]
	};
	

	unsigned int dm[4] = { datain.m_dim[0],datain.m_dim[1],datain.m_dim[2],datain.m_dim[3] };

	dataout = CDGMulDimData(datain.m_dim, datain.m_type, datain.m_Acc.m_DataMajor );

	vel2mom(dm, (float *)datain.DataPtr(0), param, (float *)dataout.DataPtr(0));
}

#include "shoot_optim3d.h"
void Fmg(CDGMulDimData & alpha,CDGMulDimData & beta, double _param[10], CDGMulDimData &X)
{
	SizeVec s= DimToSize(beta.m_Acc.m_Dim);
	unsigned int * dm = &s[0];
	double param[8] =
	{
		1.0 / _param[0],
		1.0 / _param[1],
		1.0 / _param[2],
		_param[3],
		_param[4],
		_param[5],
		_param[6],
		_param[7]
	};
	double cyc = _param[8];
	int nit = _param[9];
	{
		double v0 = param[0] * param[0],
			v1 = param[1] * param[1],
			v2 = param[2] * param[2],
			lam1 = param[4], lam2 = param[5],
			mu = param[6], lam = param[7],
			w000, wx000, wy000, wz000;
		w000 = lam2 * (6 * (v0*v0 + v1 * v1 + v2 * v2) + 8 * (v0*v1 + v0 * v2 + v1 * v2)) + lam1 * 2 * (v0 + v1 + v2);
		wx000 = 2 * mu*(2 * v0 + v1 + v2) / v0 + 2 * lam + w000 / v0;
		wy000 = 2 * mu*(v0 + 2 * v1 + v2) / v1 + 2 * lam + w000 / v1;
		wz000 = 2 * mu*(v0 + v1 + 2 * v2) / v2 + 2 * lam + w000 / v2;
		param[3] += (wx000 + wy000 + wz000)*1.2e-7 / 3.0;
	}
	float* A = (float*)alpha.DataPtr(0);
	float* b = (float *)beta.DataPtr(0);;
	VECTOR_TYPE_F scratch = VECTOR_TYPE_F::Zero(fmg3_scratchsize(dm,1));	
	X = CDGMulDimData(beta.m_Acc.m_Dim,DGDataType_Float,beta.m_Acc.m_DataMajor);
	float * x = (float*)X.DataPtr(0);
	fmg3(dm, A, b, param, cyc, nit, x, scratch.data());
}

void fmg(mwSize n0[], float *a0, float *b0, double param0[], double scal[], int c, int nit,
	float *u0, float *scratch);
mwSignedIndex fmg_scratchsize(mwSize n0[]);

void fmgN(CDGMulDimData &A, CDGMulDimData &b, double p[8] , double* p4, CDGMulDimData& rout )
{
	int nd, i;
	mwSize dm[4] = {0};
	int   cyc = 1, nit = 1;
	float *pA, *pb, *px, *scratch;
	static double param[6] = { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0 };
	double scal[256], t[3];

	nd = A.m_dim.size();

	for (int i = 0; i < nd; i++) dm[i] = b.m_dim[i];//mxGetDimensions(prhs[1])[i];
	for (int i = nd; i < 4; i++) dm[i] = 1;

	param[0] = 1 / p[0];
	param[1] = 1 / p[1];
	param[2] = 1 / p[2];
	param[3] = p[3];
	param[4] = p[4];
	param[5] = p[5];
	cyc = p[6];
	nit = p[7];


	t[0] = param[0] * param[0];
	t[1] = param[1] * param[1];
	t[2] = param[2] * param[2];
	param[3] += (param[5] * (6 * (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]) + 8 * (t[0] * t[1] + t[0] * t[2] + t[1] * t[2])) + param[4] * 2 * (t[0] + t[1] + t[2]))*1.2e-7;

	if (p4)
	{
		double *s;
		
		s = p4;
		for (int i = 0; i < dm[3]; i++)
			scal[i] = s[i];
	}
	else
	{
		for (int i = 0; i < dm[3]; i++)
			scal[i] = 1.0;
	}
		
	rout = CDGMulDimData({(int)dm[0],(int)dm[1],(int)dm[2],(int)dm[3]},DGDataType_Float, A.m_Acc.m_DataMajor);
	
	pA = (float *)A.Ptr<float>(); 
	pb = (float *)b.Ptr<float>();
	px = (float *)rout.Ptr<float>();
	scratch = (float *)malloc(fmg_scratchsize(dm)* sizeof(float));
	fmg(dm, pA, pb, param, scal, cyc, nit, px, scratch);
	free(scratch);
}


void get_mat(MAT_TYPE & ma, double M[4][3])
{
	int i, j;
	double * p = ma.data();

	for (i = 0; i < 3; i++)
		for (j = 0; j < 4; j++)
			M[j][i] = p[i + 4 * j];

	if (p[3 + 4 * 0] != 0.0 || p[3 + 4 * 1] != 0.0 || p[3 + 4 * 2] != 0.0 || p[3 + 4 * 3] != 1.0)
		ErrOutput("get_mat", "No perspective projections allowed.");
	//if (M[0][3]!=0 || M[1][3] != 0 || M[2][3] != 0 || M[3][3] != 1.0)
		//(p[3 + 4 * 0] != 0.0 || p[3 + 4 * 1] != 0.0 || p[3 + 4 * 2] != 0.0 || p[3 + 4 * 3] != 1.0)
		
}
void invdef(mwSize dim_y[3], float  y0[],
	mwSize dim_iy[3], float iy0[], double M1[4][3], double M2[4][3]);

void InvDef(CDGMulDimData &y,VECTOR_TYPE &d1 ,MAT_TYPE &M1_, MAT_TYPE &M2_, CDGMulDimData & invY)
{
	float *Y = 0, *iY = 0;
	mwSize dim_y[4], dim_iy[4];
	int i;
	double M1[4][3], M2[4][3];
		
	
	for (i = 0; i < 4; i++)
	{
		dim_y[i] = dim_iy[i] = y.m_dim[i];
	}
	
	Y = y.Ptr<float>();
		
	{		
		dim_iy[0] = d1(0);
		dim_iy[1] = d1(1);
		dim_iy[2] = d1(2);
		get_mat(M1_, M1); 
		get_mat(M2_, M2); 
	}

	invY = CDGMulDimData({ (int)dim_iy[0],(int)dim_iy[1],(int)dim_iy[2],(int)dim_iy[3] }, DGDataType_Float,y.m_Acc.m_DataMajor);

	iY = invY.Ptr<float>();

	invdef(dim_y, Y, dim_iy, iY, M1, M2);
}

float ComputeObjFunction(VECTOR_TYPE & sk,CDGMulDimData & TWrap, IN VECTOR_TYPE param)
{	
	CDGMulDimData vmom;
	{
		float sk3[3] = { 1.0 / sk(0) ,1.0 / sk(1) ,1.0 / sk(2) };
		CDGMulDimData divsk4({ 1,1,1,3 }, DGDataType_Float, TWrap.m_Acc.m_DataMajor,sk3);

		CDGMulDimData Tw = TWrap;
		Tw.BsxFuncHighDim_InPlace(BsxOperator_times, divsk4, 1);
		Vel2Mom(param.data(), Tw, vmom);
		vmom.BsxFuncHighDim_InPlace(BsxOperator_times, divsk4, 1);
	}
	
	CDGMulDimData ss;
	TWrap.BsxFuncLowDimDup(BsxOperator_times, vmom, ss);
	float llr = -0.5 * ss.SumAsFloat();
	return llr;
}
void LoadWrapBufs(IN sChannels &channels,IN sRegData regdata, IN MAT_TYPE &M, IN ChanVec &chans, IN sTPMPriors & tpm,IN int objsample,
				  OUT std::vector<sWrapBUF> &WrapBufs,IN OUT sEstimateParam & eparam,OUT MAT_TYPE &vr0)
{
	eparam.nm = 0;

	MAT_TYPE & x0 = regdata.x1;
	MAT_TYPE & y0 = regdata.x2;
	VECTOR_TYPE & z0 = regdata.x3;
	int N = chans.size();
	int framesize = regdata.x1.size();
	VECTOR_TYPE mom0 = VECTOR_TYPE::Zero(N);
	VECTOR_TYPE mom1 = VECTOR_TYPE::Zero(N);
	VECTOR_TYPE mom2 = VECTOR_TYPE::Zero(N);

	ItkImageVec fImgs(N);
	for (int n = 0; n < N; n++)
	{
		fImgs[n] = ReSampleImgAsFloat(channels[n].m_Img, objsample);
	}

	WrapBufs = std::vector<sWrapBUF>(regdata.x3.size());
	for (int z = 0; z < regdata.x3.size(); z++)
	{
		sWrapBUF & buf = WrapBufs[z];
		MAT_TYPE z1 = M(2, 0)*x0.array() + M(2, 1)*y0.array() + (M(2, 2)*z0(z) + M(2, 3));
		VECTOR_TYPE e = VoxelSize(tpm.m_TPMMat);
		e = 5 / e.array();
		buf.msk = sMask(framesize);
		double *zptr = z1.data();
		buf.msk.TraverseMsk((UpdateMaskFunc)[&](int ind, BYTE &msk)->void
		{
			msk = (byte)(zptr[ind] > e(2));
		}
		);
		vector<MAT_TYPE_F> fz(N);

		for (int n = 0; n < N; n++)
		{
			MAT_TYPE_F &fzn = fz[n];
			float* pf = fImgs[n].GetBufferAsFloat();
			SizeVec fd = fImgs[n].GetSize();
			fzn = MAT_TYPE_F(fd[0], fd[1]);
			memcpy(fz[n].data(), &pf[framesize*z], framesize * sizeof(float));
			buf.msk.TraverseMsk((UpdateMaskFunc)[&](int ind, BYTE &msk)->void
			{
				float & f = fzn.data()[ind];
				msk &= (byte) (std::isfinite(f) && f != 0 && f != -3024 && f != -1500);
			}
			);
		}
		buf.msk.ReComputeNM();
		eparam.nm += buf.msk.m_NM;
		{
			
		}
		buf.f.resize(N);
		for (int n = 0; n < N; n++)
		{
			buf.f[n] = buf.msk.MaskMat(fz[n]);
			
			{
				//to do
			}
			mom0(n) = mom0(n) + buf.msk.m_NM;
			mom1(n) = mom1(n) + buf.f[n].sum();
			mom2(n) = mom2(n) + buf.f[n].array().pow(2).sum();
		}

		buf.dat = MAT_TYPE_F::Zero(buf.msk.m_NM, eparam.Kb);
	}

	VECTOR_TYPE vv = (mom2.array() / mom0.array() - (mom1.array() / mom0.array()).pow(2)) / (eparam.Kb*eparam.Kb);
	//MAT_TYPE vv0 = (vv.array());
	vr0 = vv.asDiagonal();

}


void InitBfll(IN ChanVec &chans , IN OUT std::vector<sWrapBUF> &WrapBufs,OUT double &llrb)
{
	llrb = 0;
	int N = chans.size();
	for (int n = 0; n < N; n++)
	{
		auto & chan = chans[n];
		auto & B1 = chan.B1;
		auto & B2 = chan.B2;
		auto & B3 = chan.B3;
		auto &C = chan.C;
		auto& T = chan.T;
		
		auto TM = Eigen::Map<MAT_TYPE>((double*)T.DataPtr(0), T.m_dim[0] * T.m_dim[1], T.m_dim[2]);
		
		VECTOR_TYPE TArr = MDToVector(T.ReShape({ T.m_ElCount }));
		chan.ll = -0.5*TArr.transpose()*C*TArr;
		for (int z = 0; z < WrapBufs.size(); z++)
		{
			sWrapBUF & buf = WrapBufs[z];
			MAT_TYPE bf = transf(B1, B2, B3.row(z), TM);			
			VECTOR_TYPE tmp = buf.msk.MaskMat(bf);
			chan.ll += tmp.sum();

			if (buf.bf.size() != N) buf.bf.resize(N);
			buf.bf[n] = AsSingle(tmp).array().exp();
		}
		llrb += chan.ll;
	}

}


void defs(IN CDGMulDimData & TWrap, int z, MAT_TYPE &x0, MAT_TYPE &y0, VECTOR_TYPE &z0, MAT_TYPE &M, sMask & pmsk,
	OUT VECTOR_TYPE &x1, OUT VECTOR_TYPE &y1, OUT VECTOR_TYPE &z1)
{
	assert(TWrap.m_dim.size() == 4 && TWrap.m_type == DGDataType_Double);
	if (TWrap.m_type != DGDataType_Double)
	{
		ErrOutput("defs","m_type != DGDataType_Double");
		return;
	}
	int d0 = TWrap.m_dim[0];
	int d1 = TWrap.m_dim[1];
	MAT_TYPE x1a(d0, d1), y1a(d0, d1), z1a(d0, d1);

	int framesize = d0 * d1 * sizeof(double);
	memcpy(x1a.data(), TWrap.SelMatByDimXY({ z,0 }), framesize);
	x1a += x0;
	memcpy(y1a.data(), TWrap.SelMatByDimXY({ z,1 }), framesize);
	y1a += y0;
	memcpy(z1a.data(), TWrap.SelMatByDimXY({ z,2 }), framesize);
	z1a = z1a.array() + z0[z];

	VECTOR_TYPE x1_, y1_, z1_;
	pmsk.Mask3Mat(x1a, y1a, z1a, x1_, y1_, z1_);

	x1 = M(0, 0)*x1_.array() + M(0, 1)*y1_.array() + M(0, 2)*z1_.array() + M(0, 3);
	y1 = M(1, 0)*x1_.array() + M(1, 1)*y1_.array() + M(1, 2)*z1_.array() + M(1, 3);
	z1 = M(2, 0)*x1_.array() + M(2, 1)*y1_.array() + M(2, 2)*z1_.array() + M(2, 3);
}

MAT_TYPE transf(MAT_TYPE & B1, MAT_TYPE & B2, VECTOR_TYPE B3, Eigen::Map<MAT_TYPE> & T)
{	
	MAT_TYPE TB3 = T*B3;
	MAT_TYPE t1 = Eigen::Map<MAT_TYPE>(TB3.data(), B1.cols(), B2.cols());
	return  B1 * t1*B2.transpose();
}

MAT_TYPE transf(MAT_TYPE & B1, MAT_TYPE & B2, VECTOR_TYPE B3, CDGMulDimData & T)
{
	DimVec d2 = T.m_dim;
	MAT_TYPE t1 = Eigen::Map<MAT_TYPE>((double*)T.ptr(), d2[0] * d2[1], d2[2]);
	t1 *= B3;
	t1 = Eigen::Map<MAT_TYPE>(t1.data(), d2[0] , d2[1]);	
	return  B1 * t1*B2.transpose();
}
