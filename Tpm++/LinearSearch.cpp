#include "stdafx.h"
#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>
#include "TissueProblityMap.h"
#include "EigenMatFunc.h"
#include "llVar.h"
#include "TPMInternalFuncs.h"

void LineSearch(IN ChanVec & chans ,int n,IN MAT_TYPE & Alpha, IN MAT_TYPE & Beta, vector<sWrapBUF> & WrapBufs, 
	IN sEstimateParam & eparam, IN VECTOR_TYPE &lkp, bool busemog, IN OUT VECTOR_TYPE &WrapPrior,
	IN OUT vector<MAT_TYPE> &vr, IN OUT VECTOR_TYPE & mg, IN OUT  MAT_TYPE & mn, CVariableLL & ll)
{
	CRITICAL_SECTION cs; InitializeCriticalSection(&cs);

	sChan & chan = chans[n];
	int N = chans.size();
	auto oldT = chan.T;
	double oll = ll.GetLL();

	auto chanT = Eigen::Map<VECTOR_TYPE>((double*)chan.T.DataPtr(0), chan.T.m_ElCount);
	auto & C = chan.C;

	VECTOR_TYPE update = (Alpha + C.toDense()).llt().solve(Beta + C * chanT);

	double armijo = 1.0;
	for (int line_search = 1; line_search <= 12; line_search++)
	{
		chanT -= armijo * update;
		chan.ll = (-0.5*chanT.transpose() *C*chanT)(0);

		auto TM = Eigen::Map<MAT_TYPE>((double*)chan.T.DataPtr(0), chan.T.m_dim[0] * chan.T.m_dim[1], chan.T.m_dim[2]);	
		//Re-generate bias field, and compute terms of the objective function
		#pragma omp parallel for schedule(dynamic , 1)
		for (int z = 0; z < WrapBufs.size(); z++)
		{
			sWrapBUF & buf = WrapBufs[z];
			if (buf.msk.m_NM <= 0) continue;
			MAT_TYPE bf = transf(chan.B1, chan.B2, chan.B3.row(z), TM);			
 			VECTOR_TYPE tmp = buf.msk.MaskMat(bf);			
			buf.bf[n] = AsSingle((VECTOR_TYPE)tmp.array().exp());
			{
				CriticalSectionScoper scoper(cs);
				chan.ll += tmp.sum();
			}			
		}
		//double llrb = 0;
		eparam.llrb = 0;
		for (int n1 = 0; n1 < N; n1++) eparam.llrb += chans[n1].ll;

		ll.SetLL(eparam.llr + eparam.llrb);
		#pragma omp parallel for schedule(dynamic , 1)
		for (int z = 0; z < WrapBufs.size(); z++)
		{
			sWrapBUF & buf = WrapBufs[z];
			if (buf.msk.m_NM <= 0) continue;
			double lld = 0;
			if (busemog)
			{				
				MAT_TYPE Q;
				latent(buf.f, buf.bf, mg, mn, vr, lkp, AsDouble(buf.dat), WrapPrior, Q, lld);
				
			}
			else
			{
				MAT_TYPE_F Q;				
				latentnopar(buf.f, buf.bf, chans, buf.dat, WrapPrior, Q, lld);				
			}
			{
				CriticalSectionScoper scoper(cs);				
				//ll += lld;
				ll.IncreaseLL(lld);
			}			
		}
		if (ll.GetLL() >= oll)
		{
			printf("Bias-%d:\t%g\t%g\t%g :o)\n", n, ll.GetLL(), eparam.llr, eparam.llrb);
			break;
		}
		else
		{
			//ll = oll;
			ll.SetLL(oll);
			chan.T = oldT;
			armijo = armijo * 0.5;
			printf("Bias-%d:\t%g\t%g\t%g :o(\n", n, ll.GetLL(), eparam.llr, eparam.llrb);

		}
	}
	DeleteCriticalSection(&cs);
}
void kronutil1(int n1, int n2, int m1, int m2,
	double img[], double b1[], double b2[], double beta[]);
void kronutil2(int n1, int n2, int m1, int m2,
	double img[], double b1[], double b2[], double alpha[]);
void kronutil3(int n1x, int n1y, int n2x, int n2y, int m1, int m2,
	double img[], double b1x[], double b1y[], double b2x[], double b2y[], double alpha[]);
void UpdateAlphaBeta(MAT_TYPE & Alpha ,MAT_TYPE & Beta,sChan & chan, MAT_TYPE &wt1, MAT_TYPE &wt2,int z)
{
	VECTOR_TYPE b3 = chan.B3.row(z);
	{
		MAT_TYPE B1B2WT1 = MAT_TYPE::Zero(chan.B2.cols()*chan.B1.cols(), 1);
		kronutil1(chan.B1.cols(), chan.B2.cols(), chan.B1.rows(), chan.B2.rows(),
			wt1.data(), chan.B1.data(), chan.B2.data(), B1B2WT1.data());
		Beta += Eigen::kroneckerProduct(b3, B1B2WT1);
	}

	{
		MAT_TYPE B2B1T_WT1_B2B1 = MAT_TYPE::Zero(chan.B2.cols()*chan.B1.cols(),
			chan.B2.cols()*chan.B1.cols());

		kronutil2(chan.B1.cols(), chan.B2.cols(), chan.B1.rows(), chan.B2.rows(),
			wt2.data(), chan.B1.data(), chan.B2.data(), B2B1T_WT1_B2B1.data());
		MAT_TYPE BB3 = b3 * b3.transpose();
		Alpha += Eigen::kroneckerProduct(BB3, B2B1T_WT1_B2B1);
	}
}
#include <random>
default_random_engine generator(0);

MAT_TYPE RandnVec(int M, int N)
{
	MAT_TYPE ret(M,N);
	normal_distribution<double> distribution(0.0, 1.0);
	for (int i = 0; i < M*N; i++)
	{
		ret.data()[i] = distribution(generator);
	}
	return ret;
}

VECTOR_TYPE RandnVec(int N)
{
	//return VECTOR_TYPE::Constant(N, 1.0);
	VECTOR_TYPE ret(N);	
	
	normal_distribution<double> distribution(0.0, 1.0);
	for (int i=0;i< N;i++)
	{
		ret[i] = distribution(generator);
	}
	return ret;
}

