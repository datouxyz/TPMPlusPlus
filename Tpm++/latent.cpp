#include "stdafx.h"

#include "ParamaterDef.h"
#include <Eigen/Dense>

#include "Eigen/Cholesky"
#include "EigenMatFunc.h"
#include "TissueProblityMap.h"
#include "latent.h"
#include "TPMInternalFuncs.h"
MAT_TYPE log_likelihoods(vector<VECTOR_TYPE_F> &f, vector<VECTOR_TYPE_F> &bf, VECTOR_TYPE & mg, MAT_TYPE& mn, vector<MAT_TYPE> & vr)
{
	int K = mg.size();
	int N = f.size();
	int M = f[0].size();
	
	MAT_TYPE cr = MAT_TYPE::Zero(M, N);
	for (int n = 0; n < N; n++)
	{
		VECTOR_TYPE_F &fv = f[n];
		VECTOR_TYPE_F &bfv = bf[n];
		cr.col(n) = AsDouble(fv).array()*AsDouble(bfv).array();
	}
	MAT_TYPE L = MAT_TYPE::Zero(M, K);
	for (int k = 0; k < K; k++)
	{		
		MAT_TYPE &mpk = vr[k];
		MAT_TYPE C = mpk.llt().matrixL();

		VECTOR_TYPE mnk = mn.col(k);
		MAT_TYPE cr_mnk;
		MatVecBsxFunction<MAT_TYPE,VECTOR_TYPE>(BsxOperator_minus,true, cr, mnk, cr_mnk);

		MAT_TYPE d = /*MDToMat(d_)*/cr_mnk * C.inverse();
		VECTOR_TYPE cd = C.asDiagonal();

		MAT_TYPE dd = (d.array()*d.array());
		L.col(k) = log(mg(k)) - (N / 2.0)*log(2 * PI) - cd.array().log().sum() -
			0.5*dd.rowwise().sum().array();
	}
	return L;
}

MAT_TYPE log_likelihoods_nonpar(vector<VECTOR_TYPE_F> &f, vector<VECTOR_TYPE_F> &bf, ChanVec& chans)
{
	int K = chans[0].lik.rows();
	int Kb = chans[0].lik.cols();
	int N = chans.size();
	MAT_TYPE L = MAT_TYPE::Zero(f[0].size(), Kb);
	for (int n = 0; n < N; n++)
	{
		sChan &chan = chans[n];
		VECTOR_TYPE_F tmpf = f[n].array()*bf[n].array()*chan.interscal(1) + chan.interscal(0);
		tmpf = tmpf.array().round().min(K).max(1);
		VECTOR_TYPE_F tmp = tmpf;
		tmp = tmp.array() - 1;
		
		auto loglik = chan.alph;
		for (int k1 = 0; k1 < Kb; k1++)
		{			
			VECTOR_TYPE_F ck = loglik.col(k1)(tmp);
			L.col(k1) +=  AsDouble(ck);
		}
	}
	return L;
}

void latent(vector<VECTOR_TYPE_F> &f, vector<VECTOR_TYPE_F> &bf, VECTOR_TYPE & mg, MAT_TYPE& mn, vector<MAT_TYPE>& vr,
	VECTOR_TYPE &lkp, MAT_TYPE B, VECTOR_TYPE& wp, OUT MAT_TYPE &Q, OUT double &ll)
{
	B = log_spatial_priors<MAT_TYPE, VECTOR_TYPE>(B, wp);
	Q = log_likelihoods(f, bf, mg, mn, vr);
	int Kb = (int)lkp.maxCoeff();
	for (int k1 = 1; k1 <= Kb; k1++)
	{
		IntVEC finded = FindInVec(lkp, k1);
		for (auto k : finded)
		{
			Q.col(k) += B.col(k1 - 1);
		}
	}
	ll = safe_softmax<MAT_TYPE,VECTOR_TYPE>(Q);
}


void latentnopar(vector<VECTOR_TYPE_F> &f, vector<VECTOR_TYPE_F> &bf, ChanVec& chans,MAT_TYPE_F B, VECTOR_TYPE& wp, OUT MAT_TYPE_F &Q, OUT double &ll)
{
	B = log_spatial_priors<MAT_TYPE_F, VECTOR_TYPE_F>(B, AsSingle(wp));
	Q = AsSingle(log_likelihoods_nonpar(f, bf, chans));
	Q += B;
	ll = safe_softmax<MAT_TYPE_F, VECTOR_TYPE_F>(Q);
}

