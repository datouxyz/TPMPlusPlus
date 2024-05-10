#include "stdafx.h"
#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>
#include "TissueProblityMap.h"
#include "EigenMatFunc.h"
#include "Conv2.h"
#include "llVar.h"
#include "TPMInternalFuncs.h"
const int nMaxClusterSubiter = 20;

void EstimateClusterParameters(IN OUT VECTOR_TYPE &WrapPrior, IN OUT vector<MAT_TYPE> &vr, IN OUT VECTOR_TYPE & mg,
	IN OUT  MAT_TYPE & mn, IN OUT CVariableLL &ll,IN sEstimateParam & eparam, IN VECTOR_TYPE &lkp,
	IN std::vector<sWrapBUF> &WrapBufs, IN MAT_TYPE &vr0, int N,int &iter,int &iter1)
{	
	CRITICAL_SECTION cs; InitializeCriticalSection(&cs);

	int subiter = 1;
	for (; subiter <= nMaxClusterSubiter; subiter++)
	{
		double oll = ll.GetLL();
		VECTOR_TYPE mom0 = VECTOR_TYPE::Constant(eparam.K, eparam.tiny);
		MAT_TYPE mom1 = MAT_TYPE::Zero(N, eparam.K);
		vector<MAT_TYPE>  mom2(eparam.K, MAT_TYPE::Zero(N, N));
		MAT_TYPE mgm = MAT_TYPE::Zero(1, eparam.Kb);
		ll.SetLL(eparam.llr + eparam.llrb);
		#pragma omp parallel for schedule(dynamic , 1)
		for (int z = 0; z < WrapBufs.size(); z++)
		{
			sWrapBUF &buf = WrapBufs[z];
			if (buf.msk.m_NM <= 0) continue;
			MAT_TYPE B = AsDouble(buf.dat);
			MAT_TYPE s = 1.0 / (B * WrapPrior).array();
			
			//auto s = 1.0 / bw.array();
			MAT_TYPE Q;
			double lld = 0;
			latent(buf.f, buf.bf, mg, mn, vr, lkp, B, WrapPrior, Q, lld);
			{				
				CriticalSectionScoper scoper(cs);
				mgm += s.transpose()*B;
				ll.IncreaseLL(lld);
			}
		
			MAT_TYPE cr = MAT_TYPE::Zero(Q.rows(), N);
			for (int n = 0; n < N; n++)
			{
				VECTOR_TYPE_F vf = buf.f[n].array()*buf.bf[n].array();
				cr.col(n) = AsDouble(vf);
			}
			{
				CriticalSectionScoper scoper(cs);
				for (int k = 0; k < eparam.K; k++)
				{
					VECTOR_TYPE qk = Q.col(k);
					mom0(k) += qk.array().sum();
					mom1.col(k) += (qk.transpose() *cr).transpose();
					MAT_TYPE &mm2k = mom2[k];
					mm2k += ((MAT_TYPE)(RepMat(qk, 1, N).array() *cr.array())).transpose()*cr;
				}
			}
		}
		printf("MOG:\t subiter %d %lf\t%f\t%lf\t%lf\n", subiter,ll.GetLL(), eparam.llr, eparam.llrb,ll.GetLL()- oll);
		printf("wp:\t%f\t%f\t%f\t%f\t%f\t%f\n",WrapPrior(0),WrapPrior(1), WrapPrior(2), WrapPrior(3), WrapPrior(4), WrapPrior(5));

		for (int k = 0; k < eparam.K; k++)
		{
			IntVEC ind = FindInVec(lkp, lkp[k]);
			VECTOR_TYPE tmp(ind.size());
			for (int i = 0; i < ind.size(); i++)
			{
				int index = ind[i];
				tmp[i] = mom0(index);
			}
			mg[k] = (mom0(k) + eparam.tiny) / (tmp.array() + eparam.tiny).sum();
			mn.col(k) = mom1.col(k).array() / (mom0(k) + eparam.tiny);
			vr[k] = (mom2[k].array() - (mom1.col(k) * mom1.col(k).transpose()).array() / mom0(k)
				+ N * vr0.array()) / (mom0(k) + N);
		}
		//cout <<"vr" <<vr << endl;
		/*cout << mg<<endl;
		cout << mn << endl;
		*/
		for (int k1 = 0; k1 < eparam.Kb; k1++)
		{
			IntVEC ind = FindInVec(lkp, k1 + 1);
			VECTOR_TYPE tmp(ind.size());
			for (int i = 0; i < ind.size(); i++)
			{
				int index = ind[i];
				tmp[i] = mom0[index];
			}
			WrapPrior(k1) = (tmp.sum() + 1) / (mgm(k1) + eparam.Kb);
		}
		WrapPrior = WrapPrior.array() / WrapPrior.sum();
	
		if (subiter > 1 && ll.GetLL() - oll < (eparam.tol1*eparam.nm))
		{	
			printf("MOG:converge ll-oll@\t%f\n", ll.GetLL() - oll);
			break;
		}
	}
	
	DeleteCriticalSection(&cs);
}


void InitGaussMix(IN OUT vector<MAT_TYPE> &vr, IN OUT VECTOR_TYPE & mg, IN OUT  MAT_TYPE & mn, IN sEstimateParam & eparam, IN VECTOR_TYPE &lkp,
	IN std::vector<sWrapBUF> &WrapBufs, IN MAT_TYPE &vr0, int N)
{
	cout << "InitGaussMix" << endl;
	eparam.K = eparam.Kb;
	lkp = BuildVector<VECTOR_TYPE, double>(1, eparam.Kb, 1);
	VECTOR_TYPE mm0 = VECTOR_TYPE::Zero(eparam.Kb);
	MAT_TYPE mm1 = MAT_TYPE::Zero(N, eparam.Kb);
	vector<MAT_TYPE> mm2(eparam.Kb, MAT_TYPE::Zero(N, N));
	
	for (int z = 0; z < WrapBufs.size(); z++)
	{
		sWrapBUF & buf = WrapBufs[z];
		if(buf.msk.m_NM<=0) continue;
		MAT_TYPE cr = MAT_TYPE::Zero(buf.f[0].size(), N);
		if (cr.size() > 0)
		{
			for (int n = 0; n < N; n++)
			{
				VECTOR_TYPE_F bff = buf.f[n].array()*buf.bf[n].array();
				cr.col(n) = AsDouble(bff);
			}
			for (int k1 = 0; k1 < eparam.Kb; k1++)
			{
				VECTOR_TYPE b = AsDouble((VECTOR_TYPE_F)buf.dat.col(k1));
				mm0(k1) += b.sum();
				mm1.col(k1) += (b.transpose()*cr).transpose();
				MAT_TYPE bncr = (MAT_TYPE)(RepMat(b, 1, N).array() *cr.array());				
				mm2[k1] += bncr.transpose()*cr;
			}
		}
	}
	cout << "mm0"<<mm0 << endl;
	cout << "mm1"<< mm1 << endl;
	cout << mm2[0] << " " << mm2[1] << " " << mm2[2] << " " << mm2[3] << " " << mm2[4] << " " << mm2[5] << endl;

	mn = MAT_TYPE::Zero(N, eparam.Kb);
	
	vr = vector<MAT_TYPE>(eparam.Kb, MAT_TYPE::Zero(N, N));
	MAT_TYPE vr1 = MAT_TYPE::Zero(N, N);
	for (int k1 = 0; k1 < eparam.Kb; k1++)
	{
		mn.col(k1) = mm1.col(k1) / (mm0(k1) + eparam.tiny);
		vr1 += (mm2[k1] - mm1.col(k1)*mm1.col(k1).transpose() / mm0(k1));
	}
	cout << "mn" << mn<<endl;
	vr1 = (vr1 + N * vr0) / (mm0.sum() + N);
	vr = vector<MAT_TYPE>(eparam.Kb, vr1);
	
	mg = VECTOR_TYPE::Ones(eparam.Kb);
	
}

void ObjectiveAndDerivateMog(OUT MAT_TYPE &wt1, OUT  MAT_TYPE &wt2,IN OUT vector<MAT_TYPE> &vr, IN OUT VECTOR_TYPE &WrapPrior,
	IN OUT VECTOR_TYPE & mg, IN OUT  MAT_TYPE & mn, IN sEstimateParam & eparam, IN VECTOR_TYPE &lkp,IN sWrapBUF&buf, 
	IN vector<MAT_TYPE> &pr, IN DimVec &d ,IN int n)
{
	int N = mn.rows();
	MAT_TYPE Q;
	MAT_TYPE B = AsDouble(buf.dat);
	double lld = 0;
	latent(buf.f, buf.bf, mg, mn, vr, lkp, B, WrapPrior, Q, lld);
	vector<VECTOR_TYPE> cr(N);
	for (int n1 = 0; n1 < N; n1++)
	{
		cr[n1] = AsDouble(buf.f[n1]).array()*AsDouble(buf.bf[n1]).array();
	}
	VECTOR_TYPE w1 = VECTOR_TYPE::Zero(buf.msk.m_NM);
	VECTOR_TYPE w2 = VECTOR_TYPE::Zero(buf.msk.m_NM);
	for (int k = 0; k < eparam.K; k++)
	{
		VECTOR_TYPE qk = Q.col(k);
		VECTOR_TYPE w0 = VECTOR_TYPE::Zero(buf.msk.m_NM);
		for (int n1 = 0; n1 < N; n1++)
		{
			w0 = w0.array() + pr[k](n1, n)*(mn(n1, k) - cr[n1].array());
		}
		w1.array() += qk.array()*w0.array();
		w2.array() += qk.array()*pr[k](n, n);
	}
	wt1 = MAT_TYPE::Zero(d[0], d[1]);
	VECTOR_TYPE cr1 = -(1 + cr[n].array()*w1.array());
	buf.msk.MaskAssign(wt1, cr1);

	wt2 = MAT_TYPE::Zero(d[0], d[1]);
	VECTOR_TYPE cr2 = cr[n].array()*cr[n].array()*w2.array()+1;
	buf.msk.MaskAssign(wt2, cr2);
	
	
}

MAT_TYPE log_likelihoods(vector<VECTOR_TYPE_F> &f, vector<VECTOR_TYPE_F> &bf, VECTOR_TYPE & mg, MAT_TYPE& mn, vector<MAT_TYPE> & vr);
#include "latent.h"
void EstimateDeformationsMog(IN sEstimateParam & eparam, IN OUT VECTOR_TYPE & mg, IN OUT  MAT_TYPE & mn, vector<MAT_TYPE> & vr,
							IN std::vector<sWrapBUF> &WrapBufs, IN OUT VECTOR_TYPE &WrapPrior, IN VECTOR_TYPE &lkp, IN OUT CVariableLL & ll)
{
	CRITICAL_SECTION cs; InitializeCriticalSection(&cs);
	
	eparam.ll_const = 0;

	#pragma omp parallel for schedule(dynamic , 1)
	for (int z = 0; z < WrapBufs.size(); z++)
	{
		sWrapBUF & buf = WrapBufs[z];
		if (buf.msk.m_NM <= 0) continue;
		MAT_TYPE q = MAT_TYPE::Zero(buf.msk.m_NM, eparam.Kb);
		MAT_TYPE qt = log_likelihoods(buf.f, buf.bf, mg, mn, vr);
		
		VECTOR_TYPE max_qt = qt.rowwise().maxCoeff(); 
		{
			CriticalSectionScoper scoper(cs);
			eparam.ll_const += max_qt.sum();
		}
		
		MAT_TYPE B;
		MatVecBsxFunction<MAT_TYPE, VECTOR_TYPE>(BsxOperator_times,true ,AsDouble(buf.dat), WrapPrior, B);

		VECTOR_TYPE BSum = 1.0 / B.rowwise().sum().array();
		MAT_TYPE B_ = B;
		MatVecBsxFunction<MAT_TYPE, VECTOR_TYPE>(BsxOperator_times, false, B_, BSum, B);
		for (int k1 = 0; k1 < eparam.Kb; k1++)
		{
			IntVEC vkk = FindInVec(lkp, k1 + 1);
			for (auto k:vkk)
			{				
				VECTOR_TYPE s = (qt.col(k) - max_qt).array().exp();
				q.col(k1) += s;
			}
			buf.dat.col(k1) = AsSingle((VECTOR_TYPE)q.col(k1));
		}
		MAT_TYPE qB = q.array()*B.array() + eparam.tiny;
		MAT_TYPE qBs = qB.rowwise().sum();
		double ll_= qBs.array().log().sum();
		{
			CriticalSectionScoper scoper(cs);			
			ll.IncreaseLL(ll_);
		}
	}
	ll.IncreaseLL(eparam.ll_const);

	DeleteCriticalSection(&cs);
}
