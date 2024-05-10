#include "stdafx.h"
#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>
#include "TissueProblityMap.h"
#include "EigenMatFunc.h"
#include "Conv2.h"
#include "llVar.h"
#include "TPMInternalFuncs.h"
int nMaxHisSubit = 20;

//Estimate no mixture of gauss parameters			
void EstimateHistrogramParameters(IN OUT VECTOR_TYPE &WrapPrior,IN OUT ChanVec& chans, IN OUT CVariableLL &LL,
	IN sEstimateParam & eparam, IN MAT_TYPE &vr0,IN std::vector<sWrapBUF> &WrapBufs)
{
	CRITICAL_SECTION cs; InitializeCriticalSection(&cs);

	int N = chans.size();
	for (int n = 0; n < N; n++)
	{
		auto & chan = chans[n];
		chan.lam = VECTOR_TYPE::Zero(eparam.Kb);
		for (int k1 = 0; k1 < eparam.Kb; k1++)
		{
			chan.lam(k1) = eparam.Kb * eparam.Kb * (double)(vr0(N - 1, N - 1)*pow(chan.interscal(1), 2));
		}
	}
	// subit
	for (int subit = 1; subit <= nMaxHisSubit ; subit++)
	{
		double oll = LL.GetLL();
		/*ll = */LL.SetLL(eparam.llr + eparam.llrb);

		for (int n = 0; n < N; n++)
		{
			auto & chan = chans[n];
			MAT_TYPE lik;
			SmoothHistrogram(chan.hist, chan.lam, lik);
			chan.lik = AsSingle(lik) * chan.interscal(1);
			chan.alph = (chan.lik.array() + eparam.eps).log();
			chan.hist = MAT_TYPE::Zero(eparam.K, eparam.Kb);
		}
		VECTOR_TYPE mgm = VECTOR_TYPE::Zero(eparam.Kb);
		#pragma omp parallel for schedule(dynamic , 1)
		for (int z = 0; z < WrapBufs.size(); z++)
		{
			auto &buf = WrapBufs[z];
			if (buf.msk.m_NM <= 0)
				continue;
			MAT_TYPE_F &B = buf.dat;
			VECTOR_TYPE_F s = 1.0 / (B* AsSingle(WrapPrior)).array();
			VECTOR_TYPE_F sb = s.transpose()*B;
			
			MAT_TYPE_F Q;
			double _ll = 0;
			latentnopar(buf.f, buf.bf, chans, B, WrapPrior, Q, _ll);
			{
				CriticalSectionScoper scoper(cs);
				//ll += _ll;
				LL.IncreaseLL(_ll);
				mgm += AsDouble(sb);
			}			
			
			vector<VECTOR_TYPE_F> cr(N);
			for (int n = 0; n < N; n++)
			{
				auto & chan = chans[n];
				VECTOR_TYPE_F tmp = buf.f[n].array() * buf.bf[n].array() * chan.interscal[1] + chan.interscal[0];
				VECTOR_TYPE_F crf = tmp.array().round().max(1).min(eparam.K); 
				cr[n] = crf;				
			}
			{
				CriticalSectionScoper scoper(cs);
				for (int k1 = 0; k1 < eparam.Kb; k1++)
				{
					for (int n = 0; n < N; n++)
					{
						auto & chan = chans[n];
						VECTOR_TYPE_F QK = Q.col(k1);
						chan.hist.col(k1) += AsDouble(Accumarray<VECTOR_TYPE_F, VECTOR_TYPE_F>(cr[n], QK, eparam.K));
					}
				}
			}
		}		
		WrapPrior = (chans[0].hist.colwise().sum().array() + 1) / (mgm.array() + eparam.Kb);
		WrapPrior = WrapPrior.array() / WrapPrior.sum();
		cout << WrapPrior << endl;
		cout << "subiter:" << subit << "oll:" << oll << endl;
		if (subit > 1 && LL.GetLL() - oll < eparam.tol1*eparam.nm)
		{
			break;
		}
	}

	for (int n = 0; n < N; n++)
	{
		auto & chan = chans[n];
		MAT_TYPE lik;
		SmoothHistrogram(chan.hist, chan.lam, lik);
		chan.lik = AsSingle(lik) * chan.interscal(1);
		chan.alph = (chan.lik.array() + eparam.eps).log();
		MAT_TYPE_F kernl1(3, 1); kernl1 << 0.5, 0, -0.5; kernl1 *= chan.interscal(1);
		chan.grad1 = Conv2<MAT_TYPE_F>(chan.alph, kernl1, ConvMode_Same);		
		MAT_TYPE_F kernl2(3, 1); kernl2 << 1, -2, 1; kernl2 *= chan.interscal(1)*chan.interscal(1);
		chan.grad2 = Conv2<MAT_TYPE_F>(chan.alph, kernl2, ConvMode_Same);		
	}
	DeleteCriticalSection(&cs);
}

void InitHistrogram(IN OUT ChanVec& chans, IN sEstimateParam & eparam, IN std::vector<sWrapBUF> &WrapBufs)
{
	int N = chans.size();
	int Nz = WrapBufs.size();
	for (int n = 0; n < N; n++)
	{
		auto & chan = chans[n];
		double maxval = -DBL_MAX;
		double minval = DBL_MAX;
		for (int z = 0; z < Nz; z++)
		{
			sWrapBUF & buf = WrapBufs[z];
			if (buf.msk.m_NM <= 0)
				continue;
			maxval = max((double)buf.f[n].maxCoeff(), maxval);
			minval = min((double)buf.f[n].minCoeff(), minval);
		}
		maxval = max(maxval*1.5, -minval * 0.05);
		minval = min(minval*1.5, -maxval * 0.05);

		{
			MAT_TYPE_F mm2x2(2, 2); mm2x2 << 1, (float)minval, 1, (float)maxval;
			VECTOR_TYPE_F kk(2); kk << 1, eparam.K;
			chan.interscal = mm2x2.inverse()*kk;
		}
		cout << chan.interscal << endl;

		MAT_TYPE h0 = MAT_TYPE::Zero(eparam.K, eparam.Kb);
		for (int z1 = 0; z1 < WrapBufs.size(); z1++)
		{
			sWrapBUF & buf = WrapBufs[z1];
			if (buf.msk.m_NM <= 0)
				continue;

			VECTOR_TYPE_F cr = buf.f[n].array() *buf.bf[n].array() * chan.interscal(1) + chan.interscal(0);
			//cr = cr.array();
			cr = cr.array().round().max(1).min(eparam.K);
			//Write2DMATFloat("e:\\cr.mat", "cr_", cr);

			for (int k1 = 0; k1 < eparam.Kb; k1++)
			{
				VECTOR_TYPE_F dr = buf.dat.col(k1);
				//Write2DMATFloat("e:\\dr.mat", "dr_", dr);

				h0.col(k1) += Accumarray(cr, dr, eparam.K);
			}
		}
		chan.hist = h0;
		/*Write2DMATDouble("e:\\hist.mat","hist_",h0);
		cout << chan.hist << endl;*/
	}
}


void ObjectiveAndDerivateNoMog(OUT MAT_TYPE &wt1, OUT  MAT_TYPE &wt2, IN OUT VECTOR_TYPE &WrapPrior, 
	IN ChanVec&chans,IN sEstimateParam & eparam, IN DimVec &d, IN sWrapBUF&buf,IN int n)
{
	sChan &chan = chans[n];
	MAT_TYPE_F Q;
	double _ll = 0;
	latentnopar(buf.f, buf.bf, chans, buf.dat, WrapPrior, Q, _ll);
	VECTOR_TYPE_F cr0 = buf.f[n].array() * buf.bf[n].array();
	VECTOR_TYPE_F cr = cr0.array() * chan.interscal(1) + chan.interscal(0);
	cr = cr.array().round().max(1).min(eparam.K);
	//wt1 = VECTOR_TYPE::Zero(buf.msk.m_NM);
	wt1 = MAT_TYPE::Zero(d[0], d[1]);
	wt2 = MAT_TYPE::Zero(d[0], d[1]);

	for (int k1 = 0; k1 < eparam.Kb; k1++)
	{
		VECTOR_TYPE_F qk = Q.col(k1);
		VECTOR_TYPE_F gr1_ = chan.grad1.col(k1);
		VECTOR_TYPE_F gr1;
		AccessByIndex<VECTOR_TYPE_F, VECTOR_TYPE_F>(cr, gr1_, gr1);

		VECTOR_TYPE_F gr2_ = chan.grad2.col(k1);
		VECTOR_TYPE_F gr2;
		AccessByIndex<VECTOR_TYPE_F, VECTOR_TYPE_F>(cr, gr2_, gr2);
		gr2 = gr2.array().min(0);
		

		int ival = 0;
		VECTOR_TYPE_F qkg1 = qk.array()*(gr1.array()*cr0.array() + 1);
		VECTOR_TYPE_F qkg2 = qk.array()*(1 - gr2.array()*cr0.array().pow(2));
		buf.msk.TraverseMsk((UpdateMaskFunc)[&](int ind, BYTE &msk)->void
		{
			if ((bool)msk)
			{
				wt1(ind) -= qkg1(ival);
				wt2(ind) += qkg2(ival);
				ival++;
			}
		}
		);
	}	
}

MAT_TYPE log_likelihoods_nonpar(vector<VECTOR_TYPE_F> &f, vector<VECTOR_TYPE_F> &bf, ChanVec& chans);
#include "latent.h"
#include "llVar.h"
void EstimateDeformationsNoMog(IN sEstimateParam & eparam, ChanVec& chans,
	IN std::vector<sWrapBUF> &WrapBufs, IN OUT VECTOR_TYPE &WrapPrior, IN OUT CVariableLL & ll)
{
	CRITICAL_SECTION cs; InitializeCriticalSection(&cs);

	eparam.ll_const = 0;
	for (int z = 0; z < WrapBufs.size(); z++)
	{
		sWrapBUF & buf = WrapBufs[z];
		if (buf.msk.m_NM <= 0) continue;
		MAT_TYPE q = log_likelihoods_nonpar(buf.f, buf.bf, chans);
		VECTOR_TYPE max_q = q.rowwise().maxCoeff(); //MaxMat(q, false);
		
		{
			MAT_TYPE q_;
			MatVecBsxFunction<MAT_TYPE, VECTOR_TYPE>(BsxOperator_minus, false, q, max_q, q_);
			q = q_.array().exp();
		}
		MAT_TYPE B;
		MatVecBsxFunction<MAT_TYPE, VECTOR_TYPE>(BsxOperator_times, true, AsDouble(buf.dat), WrapPrior, B);
		{
			VECTOR_TYPE BSum = 1.0 / B.rowwise().sum().array();//  SumMat(B, false).array();
			MAT_TYPE B_ = B;
			MatVecBsxFunction<MAT_TYPE, VECTOR_TYPE>(BsxOperator_times, false, B_, BSum, B);
		}
		MAT_TYPE qB = q.array()*B.array() + eparam.tiny;
		MAT_TYPE qBs = qB.rowwise().sum(); //SumMat(qB, false);
		double ll_ = qBs.array().log().sum();		
		
		buf.dat = AsSingle(q);
		{
			CriticalSectionScoper scoper(cs);
			eparam.ll_const += max_q.sum();
			ll.IncreaseLL(ll_);
			//ll += ll_;
		}
	}
	ll.IncreaseLL(eparam.ll_const);
	//ll += eparam.ll_const;

	DeleteCriticalSection(&cs);
}