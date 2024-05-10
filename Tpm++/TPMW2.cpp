#include "stdafx.h"
#include <Eigen/Dense>

#include "TissueProblityMap.h"

#include "EigenMatFunc.h"
#include "Conv2.h"
#include "shoot_boundary.h"
#include "llVar.h"
#include <Eigen/SVD>
#include "EigenMatFunc.h"
#include "TPMInternalFuncs.h"

void defsw(IN vector<CDGMulDimData>& Coef, IN  int z, IN MAT_TYPE &x0, IN MAT_TYPE &y0, IN VECTOR_TYPE &z0, IN MAT_TYPE &M, IN MAT_TYPE &MT, DoubleVec & prm,
	OUT MAT_TYPE &x1, OUT MAT_TYPE &y1, OUT MAT_TYPE &z1)
{
	MAT_TYPE iMT = MT.inverse();
	x1 = x0.array() * iMT(0, 0) + iMT(0, 3);
	y1 = y0.array() * iMT(1, 1) + iMT(1, 3);
	z1 = z0(z) * iMT(2, 2) + iMT(2, 3)*MAT_TYPE::Ones(x1.rows(), x1.cols()).array();

	MAT_TYPE x1a, y1a, z1a;

	x1a = x0 + BSplines(Coef[0].Ptr<double>(), Coef[0].m_dim, x1, y1, z1, prm);
	y1a = y0 + BSplines(Coef[1].Ptr<double>(), Coef[1].m_dim, x1, y1, z1, prm);
	z1a = z0(z) + BSplines(Coef[2].Ptr<double>(), Coef[2].m_dim, x1, y1, z1, prm).array();

	x1 = M(0, 0)*x1a.array() + M(0, 1)*y1a.array() + M(0, 2)*z1a.array() + M(0, 3);
	y1 = M(1, 0)*x1a.array() + M(1, 1)*y1a.array() + M(1, 2)*z1a.array() + M(1, 3);
	z1 = M(2, 0)*x1a.array() + M(2, 1)*y1a.array() + M(2, 2)*z1a.array() + M(2, 3);
}

#include "latent.h"
void likelihoods(vector<MAT_TYPE> &f, vector<MAT_TYPE> &bf, VECTOR_TYPE& mg, MAT_TYPE mn, vector<MAT_TYPE> &vr, OUT MAT_TYPE & p)
{
	double eps = 2.2204e-16;

	int K = mg.size();
	int N = f.size();
	int M = f[0].size();
	MAT_TYPE cr = MAT_TYPE::Zero(M, N);
	for (int n = 0; n < N; n++)
	{
		if (bf.size() <= 0)
		{
			cr.col(n) = (VECTOR_TYPE)Eigen::Map<VECTOR_TYPE>(f[n].data(), f[n].size());
		}
		else
		{
			MAT_TYPE tmp = f[n].array()*bf[n].array();
			cr.col(n) = Eigen::Map<VECTOR_TYPE>(tmp.data(), tmp.size());
		}
	}
	p = MAT_TYPE::Ones(f[0].size(), K);
	#pragma omp parallel for schedule(dynamic , 1)
	for (int k = 0; k < K; k++)
	{
		auto Amp = mg(k) / sqrt(pow(2 * PI, (double)N)*vr[k].determinant());

		MAT_TYPE tmp;
		MatVecBsxFunction<MAT_TYPE, VECTOR_TYPE>(BsxOperator_minus, true, cr, (VECTOR_TYPE)mn.col(k), tmp);

		MAT_TYPE tmp1 = vr[k].llt().matrixL();
		MAT_TYPE d = tmp * tmp1.inverse();
		VECTOR_TYPE ds = (d.array()*d.array()).rowwise().sum();
		VECTOR_TYPE pk = (-0.5*ds.array()).exp()*Amp + eps;
		p.col(k) = pk;		
	}
}


void ZeroImageMem(sitk::Image & img)
{	
	memset(img.GetBufferAsVoid(), 0, img.GetNumberOfComponentsPerPixel()*img.GetNumberOfPixels()*img.GetSizeOfPixelComponent());
}

VECTOR_TYPE VoxelSizeSigned(MAT_TYPE & mat)
{
	MAT_TYPE R = mat.block(0, 0, 3, 3);
	//MAT_TYPE R = m33.array()*m33.array();
	MAT_TYPE C = (R*R.transpose()).llt().matrixL();
	VECTOR_TYPE ret = C.diagonal();
	if (R.determinant() < 0) ret(0) *= -1;
	return ret;
}

void GetBoundBox(MAT_TYPE & tpmmat, SizeVec &tpmDim, MAT_TYPE & bb, VECTOR_TYPE &vx)
{
	cout << tpmmat << endl;
	vx = VoxelSizeSigned(tpmmat);
	cout << vx << endl;
	VECTOR_TYPE vv0(4); vv0 << 0, 0, 0, 1;
	VECTOR_TYPE o = tpmmat.inverse() * vv0;
	cout << o << endl;
	o = (VECTOR_TYPE)o.segment(0, 3);

	VECTOR_TYPE bbmx(3); bbmx << (double)tpmDim[0], (double)tpmDim[1], (double)tpmDim[2];
	bbmx -= o;

	MAT_TYPE bb1(2, 3);
	bb1.row(0) = -vx.array()*(o.array() - 1);
	bb1.row(1) = vx.array()* bbmx.array();
	cout << bb1;
	bb = bb1;
}


bool AnyWarpTissueSetting(sTissues & tissues, bool NativeTissue_NativeSpace, bool NativeTissue_DartelImported,
	bool WarpTissue_Modulated, bool WarpTissue_UnModulated)
{
	bool ret = false;
	for (auto &ti : tissues)
	{
		if (NativeTissue_DartelImported)
		{
			ret |= ti.m_NativeTissue_DartelImported;
		}
		if (NativeTissue_NativeSpace)
		{
			ret |= ti.m_NativeTissue_NativeSpace;
		}
		if (WarpTissue_Modulated)
		{
			ret |= ti.m_WarpTissue_Modulated;
		}
		if (WarpTissue_UnModulated)
		{
			ret |= ti.m_WarpTissue_UnModulated;
		}
		if (ret) return ret;
	}
	return ret;
}
