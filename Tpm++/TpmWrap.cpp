#include "stdafx.h"

#include <Eigen/Dense>
//#include <Eigen/KroneckerProduct>
//#include <Eigen/LU>

#include "TissueProblityMap.h"
//#include "tinymatwriter.h"
//#include "SitkImagePack.h"
#include "EigenMatFunc.h"

#include "Conv2.h"
#include "shoot_boundary.h"

#include "llVar.h"
#include "TPMInternalFuncs.h"

bool WrapPass(sTpmInputs& TInput, sTPMPriors & tpm, MAT_TYPE affinemat, sWrapResults &wresult)
{
	sChannels & channels = TInput.m_Channels;
	CTimeCounter ccw("WrapPass");

	sitk::Image & V1 = channels[0].m_Img;

	sRegData regdata(V1, TInput.m_nsample, TInput.m_fwhm);
	
	sEstimateParam eparam;
	eparam.eps = 2.2204e-16;
	eparam.tiny = eparam.eps * eparam.eps;
	eparam.tol1 = 1e-4;
	eparam.K = 0;
	eparam.Kb = 0;

	SizeVec d0 = V1.GetSize();

	bool busemog = false;
	if (TInput.m_lkp.size() == 0)
	{
		eparam.K = 2000;
		eparam.Kb = tpm.m_nTissuse;
	}
	else
	{
		eparam.K = TInput.m_lkp.size();
		eparam.Kb = ToVecType(TInput.m_lkp).maxCoeff(); //MaxVec(lkp) ;
		busemog = true;
	}

	/*busemog = false;
	eparam.K = 200;
	eparam.Kb = tpm.m_nTissuse;
*/
	VECTOR_TYPE lkp = ToVecType(TInput.m_lkp);

	set_bound(1);

	MAT_TYPE & x0 = regdata.x1;
	MAT_TYPE & y0 = regdata.x2;
	VECTOR_TYPE & z0 = regdata.x3;
	DimVec d = { (int)x0.rows(),(int)x0.cols(),(int)z0.size() };

	MAT_TYPE M = tpm.m_TPMMat.inverse() * affinemat*regdata.m_MG;
	
	cout << M << endl;
	auto & sk = regdata.m_sk;
	auto & vx = regdata.m_vx;
	MAT_TYPE MT(4, 4);
	MT.row(0) << sk(0), 0, 0, (1 - sk(0));
	MT.row(1) << 0, sk(1), 0, (1 - sk(1));
	MT.row(2) << 0, 0, sk(2), (1 - sk(2));
	MT.row(3) << 0, 0, 0, 1;
	cout << MT << endl;
	
	VECTOR_TYPE param(8);
	VECTOR_TYPE skdvx = sk.array() * vx.array();
	param.segment(0, 3) = skdvx.array();
	param.segment(3, 5) = skdvx.prod()*regdata.m_FudgeFactor * ToVecType(TInput.m_reg);
	cout<<"init param" << param << endl;

	//Init Twrap
	CDGMulDimData TWrap,TWrapD;
	if (TInput.m_TWrap.m_ElCount > 0)
	{
		TWrap = TInput.m_TWrap;
		eparam.llr=ComputeObjFunction(sk, TWrap,param);
	}
	else
	{
		TWrap = CDGMulDimData({ d[0],d[1],d[2],3 }	, DGDataType_Float, DataMajor_FStyle_ColumnMajor);
		eparam.llr = 0;
	}
	

	//InitTWrap(subject, eparam.llr, regdata, param);

	ChanVec chans(channels.size());
	LoadChans(channels, regdata, chans);

	int N = chans.size();

	//double ll = -DBL_MAX;
	CVariableLL LL; LL.SetLL(-DBL_MAX);
	

	MAT_TYPE vr0;
	std::vector<sWrapBUF> WrapBufs;
	LoadWrapBufs(channels, regdata, M, chans, tpm, TInput.m_nsample, WrapBufs, eparam, vr0);

	InitBfll(chans, WrapBufs, eparam.llrb);
	//Create initial bias field	

	VECTOR_TYPE WrapPrior = TInput.m_WrapPrior.size() > 0 ? TInput.m_WrapPrior : VECTOR_TYPE::Ones(eparam.Kb).array() / eparam.Kb;

	VECTOR_TYPE mg;
	MAT_TYPE mn;
	vector<MAT_TYPE> vr;	

	for (int iter = 1; iter <= 30; iter++)
	{	
		CDGMulDimData TWrapD; TWrap.ToDouble(TWrapD);

		{
			//CDGMulDimData TWrapD; regdata.TWrap.ToDouble(TWrapD);
			for (int z = 0; z < z0.size(); z++)
			{
				sWrapBUF & buf = WrapBufs[z];
				if (buf.msk.m_NM <= 0) continue;

				VECTOR_TYPE x1, y1, z1;
				defs(TWrapD, z, x0, y0, z0, M, buf.msk, x1, y1, z1);
				VECTOR_TYPES b;
				tpm.SamplePriors(x1, y1, z1, b);
				for (int k1 = 0; k1 < eparam.Kb; k1++)
				{
					buf.dat.col(k1) = AsSingle(b[k1]);
				}
			}
		}
				
		if (iter == 1)
		{
			//busemog = false;
			if (busemog)
			{
				//��ʼ���Ƹ�˹����
				if (TInput.mg.size() > 0 && TInput.mn.size() > 0 && TInput.vr.size() > 0)
				{
					mn = TInput.mn;
					mg = TInput.mg;
					vr = TInput.vr;
				}
				else
				{
					InitGaussMix(vr, mg, mn, eparam, lkp, WrapBufs, vr0, N);
				}
			}
			else
			{
				InitHistrogram(chans, eparam, WrapBufs);
			}
		}// iter0 init
		double ooll = DBL_MAX;
		//iter1
		for (int iter1 = 1; iter1 <= 8; iter1++)
		{
			if (busemog)
			{
				//Estimate cluster parameters				
				EstimateClusterParameters(WrapPrior, vr, mg, mn, LL, eparam, lkp, WrapBufs, vr0, N,iter,iter1);
				cout << "iter:" << iter << " iter1:" << iter1 << " ll:" << LL.GetLL() << endl;
			}//end of mog	
			else
			{
				EstimateHistrogramParameters(WrapPrior, chans, LL, eparam, vr0, WrapBufs);
			}

			if (iter1 > 1 && (LL.GetLL() - ooll) <= 2 * eparam.tol1*eparam.nm)//!((ll - ooll) > 2 * eparam.tol1*eparam.nm))
				break;
			ooll = LL.GetLL();

			//% Estimate bias
			vector<MAT_TYPE> pr;

			if (busemog)
			{
				pr.resize(vr.size());
				for (int i = 0; i < vr.size(); i++)
				{
					pr[i] = vr[i].inverse();
				}
			}
			//for (int subit = 1; subit <= 1; subit++)
			//{
			for (int n = 0; n < N; n++)
			{
				sChan & chan = chans[n];
				int d3 = chan.T.m_ElCount;
				if (d3 <=0) continue;
					
				//Compute objective function and its 1st and second derivatives
				MAT_TYPE Alpha = MAT_TYPE::Zero(d3, d3);
				MAT_TYPE Beta = MAT_TYPE::Zero(d3, 1);
				CRITICAL_SECTION cs; InitializeCriticalSection(&cs);

				#pragma omp parallel for schedule(dynamic , 1)
				for (int z = 0; z < WrapBufs.size(); z++)
				{
					sWrapBUF & buf = WrapBufs[z];
					if (buf.msk.m_NM <= 0) continue;
					MAT_TYPE wt1; MAT_TYPE wt2 ;
					if (busemog)
					{
						ObjectiveAndDerivateMog(wt1,wt2,vr,WrapPrior,
							mg, mn, eparam,lkp,buf,pr,d,n);
					}
					else
					{
						ObjectiveAndDerivateNoMog(wt1,wt2,WrapPrior,
							chans,eparam,d,buf,n);
					}
					{
						CriticalSectionScoper scoper(cs);					
						UpdateAlphaBeta(Alpha, Beta, chan, wt1, wt2, z);
					}					
				}				
				DeleteCriticalSection(&cs);
				LineSearch(chans, n,Alpha, Beta, WrapBufs,eparam, lkp, busemog, WrapPrior,vr,mg,mn, LL);
			}
			//}// endof subit
			if (iter == 1 && iter1 == 1)
			{
				if (busemog && TInput.m_lkp.size()!=lkp.size() )
				{	
					lkp = ToVecType(TInput.m_lkp);
					eparam.K = lkp.size();
					eparam.Kb = lkp.maxCoeff();

					ReFillMix(vr, mg, mn, lkp, eparam, N);
				}				
			}
		}// end of iter1

		// % Estimate deformations		
		//ll = eparam.llr + eparam.llrb;
		LL.SetLL(eparam.llr + eparam.llrb);
		if (busemog)
		{
			EstimateDeformationsMog(eparam,mg,mn,vr,WrapBufs,WrapPrior,lkp, LL);
		}
		else
		{
			EstimateDeformationsNoMog(eparam, chans,WrapBufs, WrapPrior, LL);
		}
		bool bDecreaseParamByIter = TInput.m_TWrap.m_ElCount <= 0;
		if (iter == 1)
			LL.DoubleToSingle();

		UpdateTWrap(LL,x0, y0, z0, WrapBufs,	param, eparam, WrapPrior,
			TWrap, tpm, MT,M, iter,bDecreaseParamByIter,sk);
				
		cout <<"iter"<<iter<<" ll:" << LL.GetLL()<<endl ;
		//% Finished
		if (iter > 9 && !((LL.GetLL() - ooll) > 2 * eparam.tol1*eparam.nm))
		break;

	}// endof main iter
		
	
	wresult.m_tpm = tpm;
	wresult.m_Affine = affinemat;
	wresult.m_lkp = lkp;
	wresult.m_MT = MT;
	wresult.m_Twrap = TWrap;

	wresult.m_T.resize(N);
	for (int ic=0;ic<chans.size();ic++)
	{
		wresult.m_T[ic] = chans[ic].T;
	}
	
	wresult.m_WrapPrior = WrapPrior;

	if (busemog)
	{
		wresult.m_mn = mn;
		wresult.m_mg = mg;
		wresult.m_vr = vr;
	}
	else
	{
		for (int n = 0; n < N; n++)
		{
			wresult.m_intensity[n].m_lik = chans[n].lik;
			wresult.m_intensity[n].m_interscal = chans[n].interscal;
		}
	}
	wresult.m_LL = LL;
	
	return true;
}
