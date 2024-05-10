#include "stdafx.h"
#include <Eigen/Dense>
//#include <unsupported/Eigen/MatrixFunctions>
#include "TissueProblityMap.h"
#include "EigenMatFunc.h"
#include "TPMInternalFuncs.h"
const double tiny = 1e-4;

IntVEC FindInVec(VECTOR_TYPE &lkp, double key)
{
	IntVEC ret;
	for (int i = 0; i < lkp.size(); i++)
	{
		if (key == lkp[i])
			ret.push_back(i);
	}
	return ret;
}



int MaxVec(IntVEC & lkp)
{
	int m = INT_MIN;
	for (auto i : lkp)
	{
		m = max(m, i);
	}
	return m;
}

void sTPMPriors::SamplePriorsDerive(VECTOR_TYPE &_x1, VECTOR_TYPE &_x2, VECTOR_TYPE &_x3,
	OUT VECTOR_TYPES &s, OUT VECTOR_TYPES &ds1, OUT VECTOR_TYPES &ds2, OUT VECTOR_TYPES &ds3)
{
	DimVec dim = SizeToDim(m_size);
	int Kb = m_Datas.size();
	ds1.resize(Kb);
	ds2.resize(Kb);
	ds3.resize(Kb);
	s.resize(Kb);

	int dx = _x1.size();
	VECTOR_TYPE tot = VECTOR_TYPE::Zero(dx);
	VECTOR_TYPE ones = VECTOR_TYPE::Constant(dx, 1);
	sMask Mask1(dx);
	sMask Mask2(dx);

	for (int ix = 0; ix < dx; ix++)
	{
		bool msk = _x1[ix] >= 1 && _x1[ix] <= m_size[0] &&
			_x2[ix] >= 1 && _x2[ix] <= m_size[1] &&
			_x3[ix] >= 1 && _x3[ix] <= m_size[2];
		Mask1.SetMaskAndNM(ix, msk);
		Mask2.SetMaskAndNM(ix, _x3[ix] < 1);
	}

	VECTOR_TYPE x1 = Mask1.MaskVector(_x1);
	VECTOR_TYPE x2 = Mask1.MaskVector(_x2);
	VECTOR_TYPE x3 = Mask1.MaskVector(_x3);
	VECTOR_TYPE a, da1, da2, da3;

	DoubleVec degs = { m_deg,m_deg,m_deg, 0, 0, 0 };
	for (int k = 0; k < Kb; k++)
	{		
		BSplinesD(m_Datas[k].data(), dim, degs, x1, x2, x3, a, da1, da2, da3);	
		if (k == Kb - 1)
		{
			s[k] = VECTOR_TYPE::Ones(dx);
		}
		else
		{
			s[k] = VECTOR_TYPE::Zero(dx) + VECTOR_TYPE::Constant(dx, tiny);
		}

		s[k] = ones * m_bg2[k];

		VECTOR_TYPE expa = a.array().exp();
		Mask1.MaskMatchAndAssign(s[k], expa);
		Mask2.MaskAssign(s[k], m_bg1[k]);
		tot += s[k];
		ds1[k] = VECTOR_TYPE::Zero(dx);
		ds2[k] = VECTOR_TYPE::Zero(dx);
		ds3[k] = VECTOR_TYPE::Zero(dx);

		Mask1.MaskMatchAndAssign(ds1[k], da1);
		Mask1.MaskMatchAndAssign(ds2[k], da2);
		Mask1.MaskMatchAndAssign(ds3[k], da3);
	}

	sMask msk(tot.size());
	msk.MaskOnInfinite(tot);
	msk.MaskAssign(tot, 1);

	da1 = VECTOR_TYPE::Zero(dx);
	da2 = VECTOR_TYPE::Zero(dx);
	da3 = VECTOR_TYPE::Zero(dx);


	for (int k = 0; k < Kb; k++)
	{
		msk.MaskAssign(s[k], m_bg1[k]);
		s[k] = s[k].array() / tot.array(); //Div(s[k], tot);
		VECTOR_TYPE d1 = (s[k].array() * ds1[k].array());
		VECTOR_TYPE d2 = (s[k].array() * ds2[k].array());
		VECTOR_TYPE d3 = (s[k].array() * ds3[k].array());
		da1 += d1;
		da2 += d2;
		da3 += d3;
	}
	for (int k = 0; k < Kb; k++)
	{
		ds1[k] = s[k].array() * (ds1[k] - da1).array(); msk.MaskAssign(ds1[k], 0);
		ds2[k] = s[k].array() * (ds2[k] - da2).array(); msk.MaskAssign(ds2[k], 0);
		ds3[k] = s[k].array() * (ds3[k] - da3).array(); msk.MaskAssign(ds3[k], 0);
	}

}

void sTPMPriors::SamplePriors(VECTOR_TYPE &_x1, VECTOR_TYPE &_x2, VECTOR_TYPE &_x3, OUT VECTOR_TYPES & s)
{
	DimVec dim = SizeToDim(m_size);
	int dx = _x1.size();
	int Kb = m_Datas.size();
	s.resize(Kb);
	sMask Mask1(dx);
	sMask Mask2(dx);

	for (int ix = 0; ix < dx; ix++)
	{
		bool msk = _x1[ix] >= 1 && _x1[ix] <= m_size[0] &&
			_x2[ix] >= 1 && _x2[ix] <= m_size[1] &&
			_x3[ix] >= 1 && _x3[ix] <= m_size[2];
		Mask1.SetMaskAndNM(ix, msk);
		Mask2.SetMaskAndNM(ix, _x3[ix] < 1);		
	}

	VECTOR_TYPE x1 = Mask1.MaskVector(_x1);
	VECTOR_TYPE x2 = Mask1.MaskVector(_x2);
	VECTOR_TYPE x3 = Mask1.MaskVector(_x3);

	{
		VECTOR_TYPE tot = VECTOR_TYPE::Zero(dx);
		VECTOR_TYPE ones = VECTOR_TYPE::Constant(dx, 1);
		for (int k = 0; k < Kb; k++)
		{
			DoubleVec degs = { m_deg,m_deg,m_deg, 0, 0, 0 };
			VECTOR_TYPE a = BSplines(m_Datas[k].data(), dim, x1, x2, x3, degs);
			VECTOR_TYPE & sk = s[k];
			sk = ones * m_bg2[k];
			a = a.array().exp();
			Mask1.MaskMatchAndAssign(sk, a);
			Mask2.MaskAssign(sk, m_bg1[k]);
			tot += sk;
		}		
		sMask infinitemsk(dx);
		infinitemsk.MaskOnInfinite(tot);
		infinitemsk.MaskAssign(tot, 1.0);		
		for (int k = 0; k < Kb; k++)
		{
			VECTOR_TYPE & sk = s[k];
			infinitemsk.MaskAssign(sk, m_bg2[k]);
			for (int ix = 0; ix < dx; ix++)
			{				
				sk[ix] /= tot[ix];
			}			
		}
	}	
}

void sTPMPriors::LoadPriors(sitk::Image & tpmimg)
{
	m_TpmSpace = sTpmVol3D(tpmimg);
	SizeVec size4d = tpmimg.GetSize();
	m_nTissuse = size4d[3];
	double deg = 1;
	m_size = { size4d[0],size4d[1],size4d[2] };
	DimVec dim = SizeToDim(m_size);
	m_TPMMat = m_TpmSpace.GetMat(true, true);
	cout << m_TPMMat << endl;
	m_Datas.resize(size4d[3]);
	int volsize = DimElCount(dim);
	int slicesize = m_size[0] * m_size[1];
	for (int k = 0; k < size4d[3]; k++)
	{
		m_Datas[k].resize(volsize);
		m_Datas[k] = VECTOR_TYPE::Zero(volsize);
	}
	double * pvol4d = tpmimg.GetBufferAsDouble();
	for (int iz = 0; iz < m_size[2]; iz++)
	{
		VECTOR_TYPE s = VECTOR_TYPE::Zero(slicesize);
		for (int ik = 0; ik < m_nTissuse; ik++)
		{
			double * tpmFrame = &m_Datas[ik][iz*slicesize];
			double * slice = &pvol4d[volsize*ik + slicesize * iz];

			VECTOR_TYPE tmp(slicesize);
			memcpy(tmp.data(), slice, slicesize * sizeof(double));

			for (int i = 0; i < slicesize; i++)
			{
				tpmFrame[i] = clamp<double>(tmp[i], 0, 1);
			}
			s += tmp;
		}

		for (int i = 0; i < slicesize; i++)
		{
			if (s[i] > 1)
			{
				for (int ik = 0; ik < m_nTissuse; ik++)
				{
					double *dataref = &m_Datas[ik][iz*slicesize];
					dataref[i] /= s[i];
				}
			}			
		}
	}
	m_bg1 = VECTOR_TYPE::Zero(m_nTissuse);
	m_bg2 = VECTOR_TYPE::Zero(m_nTissuse);
	DoubleVec degs = { deg ,deg ,deg,0,0,0 };

	{
		for (int ik = 0; ik < m_nTissuse; ik++)
		{
			double * pvol = m_Datas[ik].data();
			VECTOR_TYPE st(slicesize); memcpy(st.data(), &pvol[0], st.size() * sizeof(double));
			VECTOR_TYPE ed(slicesize); memcpy(ed.data(), &pvol[(m_size[2] - 1)*slicesize], ed.size() * sizeof(double));

			m_bg1[ik] = st.mean();
			m_bg2[ik] = ed.mean();

			VECTOR_TYPE newvol = m_Datas[ik];
			for (auto &sl : newvol)
			{
				sl = log(sl + tiny);
			}
			BSplinec(degs, newvol.data(), dim, m_Datas[ik].data());
		}
	}
	m_deg = deg + 1;
}


VECTOR_TYPE SmoothKernel(double fwhm,int xstart, int xend, KernelSmoothMode mode,double eps)
{
	double s = pow(fwhm / sqrt(8 * log(2)),2) + eps;
	//cout << s<<endl;
	int count = xend - xstart + 1;
	VECTOR_TYPE krn(count);
	switch (mode)
	{
		case KernelSmoothMode_nearest:
		{
			double w1 = 1.0 / sqrt(2 * s);			
			for(int _x= xstart;_x<=xend;_x++)
			{
				int ix = _x - xstart;				
				krn(ix) = 0.5*(erf(w1*(_x + 0.5)) - erf(w1*(_x - 0.5)));
				if (krn(ix) < 0)krn(ix) = 0;
			}
		}
		break;
		case KernelSmoothMode_linear:
		{
			double w1 = 0.5*sqrt(2 / s);
			double w2 = -0.5 / s;
			double w3 = sqrt(s / 2 / PI);			
			for (int _x = xstart; _x <= xend; _x++)
			{
				int ix = _x - xstart;
				krn(ix) = 0.5*(erf(w1*(_x + 1))*(_x + 1)
					+ erf(w1*(_x - 1))*(_x - 1) 
					- 2 * erf(w1*_x)*_x)
					+ w3 *(exp(w2* pow(_x + 1,2))
					+ exp(w2* pow(_x - 1,2))
					- 2 * exp(w2* pow(_x,2)));				
				if (krn(ix) < 0)krn(ix) = 0;
			}
		}
		break;
	}
	//cout << krn;
	return krn;
}

MAT_TYPE EigenMatFromDoubleVec(DoubleVec & dp,int row,int column)
{
	MAT_TYPE ret(row,column);
	assert(dp.size() == row * column);
	for(int irow=0;irow<row;irow++)
	for(int icolumn = 0;icolumn < column;icolumn++)
	{
		ret(irow,icolumn) = dp[irow*column+icolumn];
	}
	return ret;
}
MAT_TYPE Mat4x4FromAffineParam(DoubleVec & dp)
{
	MAT_TYPE ret(4,4);
	
	ret << 
		dp[0], dp[1], dp[2], dp[9],
		dp[3], dp[4], dp[5], dp[10],
		dp[6], dp[7], dp[8], dp[11],
		0,0,0,1;
	
	return ret;
}

void sMask::MaskAssign(VECTOR_TYPE& target, double val)
{
	int n=m_Msk.size();
	for (int i = 0; i < n; i++)
	{				
		target[i] = (bool)m_Msk[i]? val: target[i];		
	}
}

void sMask::MaskAssign(MAT_TYPE& target, VECTOR_TYPE &val)
{
	assert(m_Msk.size() == target.size() && val.size() == m_NM);

	double* dptr = target.data();
	int ival = 0;
	for (int iel = 0; iel < m_Msk.size(); iel++)
	{
		if ((int)m_Msk[iel]>0)
		{
			dptr[iel] = val[ival];
			ival++;
		}
	}
}

bool sMask::MaskMatchAndAssign(VECTOR_TYPE& target, VECTOR_TYPE& src)
{
	if (src.size() != m_NM || target.size() < m_NM) {
		ErrOutput("Mask", u8"MaskMatch ƥ��ʧ��");
		return false;
	}

	int na = target.size();
	int isrc = 0;
	for (int i = 0; i < na; i++)
	{
		if ((bool)m_Msk[i])
		{
			target[i] = src[isrc];
			isrc++;
		}
	}
	return true;
}

int NumberMask(ByteVec &mask)
{
	int ret = 0;
	for (auto&b : mask)
		ret += ((int)b != 0);

	return ret;
}

void sMask::ReComputeNM()
{
	m_NM = NumberMask(m_Msk);
}

VECTOR_TYPE_F sMask::MaskVector(VECTOR_TYPE_F& vvector)
{
	VECTOR_TYPE_F ret(m_NM);
	int iret = 0;
	for (int i = 0; i < m_Msk.size(); i++)
	{
		if ((bool)m_Msk[i])
		{
			ret[iret] = vvector[i];
			iret++;
		}	
	}
	
	return ret;
}

VECTOR_TYPE sMask::MaskVector(VECTOR_TYPE & vvector)
{
	VECTOR_TYPE ret(m_NM);
	int iret = 0;
	for (int i = 0; i < m_Msk.size(); i++)
	{
		if ((bool)m_Msk[i])
		{
			ret[iret] = vvector[i];
			iret++;
		}
	}
	return ret;
}

VECTOR_TYPE_F sMask::MaskMat(MAT_TYPE_F & mat)
{
	float * pdata = mat.data();
	VECTOR_TYPE_F ret(m_NM);
	int iret = 0;
	for (int imask = 0; imask < m_Msk.size(); imask++)
	{
		if ((bool)m_Msk[imask])
		{
			ret[iret] = pdata[imask];
			iret++;
		}
	}
	return ret;
}
VECTOR_TYPE sMask::MaskMat(MAT_TYPE & mat)
{
	
	double * pdata = mat.data();
	VECTOR_TYPE ret(m_NM);
	int iret = 0;
	for (int imask = 0; imask < m_Msk.size(); imask++)
	{
		if ((bool)m_Msk[imask])
		{
			ret[iret] = pdata[imask];
			iret++;
		}
	}	
	return ret;
}

void sMask::Mask3Mat(MAT_TYPE & mat1, MAT_TYPE & mat2, MAT_TYPE & mat3, OUT VECTOR_TYPE & v1, OUT VECTOR_TYPE & v2, OUT VECTOR_TYPE & v3)
{
	v1 = VECTOR_TYPE(m_NM);
	v2= VECTOR_TYPE(m_NM);
	v3= VECTOR_TYPE(m_NM);

	int iret = 0;
	for (int i = 0; i < m_Msk.size(); i++)
	{
		if ((bool)m_Msk[i])
		{
			v1[iret] = mat1.data()[i];
			v2[iret] = mat2.data()[i];
			v3[iret] = mat3.data()[i];
			iret++;
		}
	}
	
}

void sMask::Mask3Vector(VECTOR_TYPE& vvector1, VECTOR_TYPE& vvector2, VECTOR_TYPE& vvector3)
{
	VECTOR_TYPE v1(m_NM);
	VECTOR_TYPE v2(m_NM);
	VECTOR_TYPE v3(m_NM);
	int iret = 0;
	for (int i = 0; i < m_Msk.size(); i++)
	{
		if ((bool)m_Msk[i])
		{
			v1[iret] = vvector1[i];
			v2[iret] = vvector2[i];
			v3[iret] = vvector3[i];
			iret++;
		}
	}
	vvector1 = v1;
	vvector2 = v2;
	vvector3 = v3;
}


VECTOR_TYPE Accumarray(VECTOR_TYPE_F &ind, VECTOR_TYPE_F &data, int outsize)
{
	VECTOR_TYPE ret = VECTOR_TYPE::Zero(outsize);

	if (ind.size() != data.size())
		return ret;
	outsize = clamp<int>(outsize, 0,(int)data.size());
	for (int i = 0; i < data.size(); i++)
	{
		ret[ind[i] - 1] += data[i];
	}	
	return ret;
}

VECTOR_TYPE Accumarray(VECTOR_TYPE &ind, VECTOR_TYPE &data, int outsize)
{
	VECTOR_TYPE ret = VECTOR_TYPE::Zero(outsize);

	if (ind.size() != data.size())
		return ret;
	outsize = clamp<int>(outsize, 0, (int)data.size());
	for (int i = 0; i < data.size(); i++)
	{
		ret[ind[i]-1] += data[i];
	}	
	return ret;
}

VECTOR_TYPE Accumarray(VECTOR_TYPE &ind, double singledata, int outsize)
{
	VECTOR_TYPE ret = VECTOR_TYPE::Zero(outsize);

	outsize = clamp<int>(outsize, 0, ind.size());
	for (int i = 0; i < ind.size(); i++)
	{
		ret[ind[i]] += singledata;
	}
	return ret;
}

VECTOR_TYPE CumSum(VECTOR_TYPE & input)
{
	VECTOR_TYPE rets = VECTOR_TYPE::Zero(input.size());
	//cumsum

	int cumh = 0;
	for (int i = 0; i < input.size(); i++)
	{
		double & inpu = input[i];
		cumh += inpu;
		rets[i] = cumh;
	}
	return rets;
}

MAT_TYPE Derivs(MAT_TYPE &MF,VECTOR_TYPE &P/**/, MAT_TYPE &MG)
{
	MAT_TYPE R = MAT_TYPE::Zero(12, 12);
	MAT_TYPE _M0 = MF.inverse() * Param2Matrix(P)*MG;
	MAT_TYPE M0 = _M0.block(0, 0, 3, 4);
	for (int i = 0; i < 12; i++)
	{
		double dp = 0.0000001;
		VECTOR_TYPE P1 = P;		
		P1(i) += dp;
		MAT_TYPE _M1 = MF.inverse()*Param2Matrix(P1)*MG;		
		MAT_TYPE M1 = _M1.block(0, 0, 3, 4);		
		MAT_TYPE M10 = (M1 - M0);
		VECTOR_TYPE delta = Eigen::Map<VECTOR_TYPE>(M10.data(), M10.size());		
		R.col(i) = delta.array() / dp;
	}	
	return R;
}

void ApplyMatlabIndexOffset(MAT_TYPE & mat)
{
	MAT_TYPE moffset = MAT_TYPE::Identity(4, 4);
	moffset.col(3) << -1, -1, -1, 1;
	mat = mat * moffset;
}

MAT_TYPE_F AsSingle(MAT_TYPE & m)
{
	MAT_TYPE_F ret(m.rows(), m.cols());
	for (int i = 0; i < ret.size(); i++)
	{
		ret.data()[i] = m.data()[i];
	}
	return ret;
}
VECTOR_TYPE_F AsSingle(VECTOR_TYPE & v)
{
	VECTOR_TYPE_F ret(v.size());
	for (int i = 0; i < ret.size(); i++)
	{
		ret.data()[i] = v.data()[i];
	}
	return ret;
}
MAT_TYPE AsDouble(const MAT_TYPE_F & m)
{
	MAT_TYPE ret(m.rows(), m.cols());
	for (int i = 0; i < ret.size(); i++)
	{
		ret.data()[i] = m.data()[i];
	}
	return ret;
}
VECTOR_TYPE AsDouble(const VECTOR_TYPE_F & v)
{
	VECTOR_TYPE ret(v.size());
	for (int i = 0; i < ret.size(); i++)
	{
		ret.data()[i] = v.data()[i];
	}
	return ret;
}
#include "DGMulDimData.h"

MAT_TYPE RepMat(VECTOR_TYPE & v, int r, int c)
{
	MAT_TYPE ret = MAT_TYPE::Zero(v.rows()*r,v.cols()*c);
	for (int ir=0;ir<r;ir++)
	for (int ic=0;ic<c;ic++)
	{
		ret.block(ir * v.rows(),ic*v.cols(), v.rows(), v.cols()) = 
			v;
	}
	return ret;
}

MAT_TYPE RepMat(MAT_TYPE & v, int r, int c)
{
	MAT_TYPE ret = MAT_TYPE::Zero(v.rows()*r, v.cols()*c);
	for (int ir = 0; ir < r; ir++)
		for (int ic = 0; ic < c; ic++)
		{
			ret.block(ir * v.rows(), ic*v.cols(), v.rows(), v.cols()) =
				v;
		}
	return ret;
}


