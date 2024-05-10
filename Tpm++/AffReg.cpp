#include "stdafx.h"

#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>
#include <Eigen/LU>
#include "ParamaterDef.h"
#include <Eigen/Householder>
#include "TissueProblityMap.h"
#include "EigenMatFunc.h"
#include "Conv2.h"
#include "TPMInternalFuncs.h"

void AffinePriors(PriorsType type, OUT VECTOR_TYPE &mu, OUT MAT_TYPE & isigma)
{
	mu.resize(6);
	isigma.resize(6, 6);

	switch (type)
	{
	case PriorsType_mni: //, % For registering with MNI templates...		
		mu << 0.0667, 0.0039, 0.0008, 0.0333, 0.0071, 0.1071;
		isigma <<
			0.0902, -0.0345, -0.0106, -0.0025, -0.0005, -0.0163,
			-0.0345, 0.7901, 0.3883, 0.0041, -0.0103, -0.0116,
			-0.0106, 0.3883, 2.2599, 0.0113, 0.0396, -0.0060,
			-0.0025, 0.0041, 0.0113, 0.0925, 0.0471, -0.0440,
			-0.0005, -0.0103, 0.0396, 0.0471, 0.2964, -0.0062,
			-0.0163, -0.0116, -0.0060, -0.0440, -0.0062, 0.1144;
		isigma *= 1e4;
		//for (auto& ig : isig) { ig *= 1e4; };
		break;

	case PriorsType_imni:// % For registering with MNI templates...
		mu << 0.0667, 0.0039, 0.0008, 0.0333, 0.0071, 0.1071;
		mu *= -1;
		isigma <<
			0.0902, -0.0345, -0.0106, -0.0025, -0.0005, -0.0163,
			-0.0345, 0.7901, 0.3883, 0.0041, -0.0103, -0.0116,
			-0.0106, 0.3883, 2.2599, 0.0113, 0.0396, -0.0060,
			-0.0025, 0.0041, 0.0113, 0.0925, 0.0471, -0.0440,
			-0.0005, -0.0103, 0.0396, 0.0471, 0.2964, -0.0062,
			-0.0163, -0.0116, -0.0060, -0.0440, -0.0062, 0.1144;
		isigma *= 1e4;
		//for (auto& ig : isig) { ig *= 1e4; };
		break;

	case PriorsType_rigid:
		mu << 0, 0, 0, 0, 0, 0;
		isigma <<
			1, 0, 0, 0, 0, 0,
			0, 1, 0, 0, 0, 0,
			0, 0, 1, 0, 0, 0,
			0, 0, 0, 1, 0, 0,
			0, 0, 0, 0, 1, 0,
			0, 0, 0, 0, 0, 1;
		isigma *= 1e8;
		//for (auto& ig : isig) { ig *= 1e8; };
		break;
	case PriorsType_subject:
		mu << 0, 0, 0, 0, 0, 0;
		isigma <<
			0.8876, 0.0784, 0.0784, -0.1749, 0.0784, -0.1749,
			0.0784, 5.3894, 0.2655, 0.0784, 0.2655, 0.0784,
			0.0784, 0.2655, 5.3894, 0.0784, 0.2655, 0.0784,
			-0.1749, 0.0784, 0.0784, 0.8876, 0.0784, -0.1749,
			0.0784, 0.2655, 0.2655, 0.0784, 5.3894, 0.0784,
			-0.1749, 0.0784, 0.0784, -0.1749, 0.0784, 0.8876;
		isigma *= 1e3;
		break;
	case PriorsType_eastern:
		mu << 0.0719, -0.0040, -0.0032, 0.1416, 0.0601, 0.2578;
		isigma <<
			0.0757, 0.0220, -0.0224, -0.0049, 0.0304, -0.0327,
			0.0220, 0.3125, -0.1555, 0.0280, -0.0012, -0.0284,
			-0.0224, -0.1555, 1.9727, 0.0196, -0.0019, 0.0122,
			-0.0049, 0.0280, 0.0196, 0.0576, -0.0282, -0.0200,
			0.0304, -0.0012, -0.0019, -0.0282, 0.2128, -0.0275,
			-0.0327, -0.0284, 0.0122, -0.0200, -0.0275, 0.0511;
		isigma *= 1e4;
		break;
	case PriorsType_none:
		mu << 0, 0, 0, 0, 0, 0;
		isigma << 0, 0, 0, 0, 0, 0;
		break;
	}
}

MAT_TYPE AffReg(sRegData & regdata, sTPMPriors & priors, MAT_TYPE M)
{
	//MAT_TYPE M = regdata.m_Affine;// MAT_TYPE::Identity(4, 4);
	CTimeCounter afftimer("AffReg");

	int ntissuse = priors.m_Datas.size();
	MAT_TYPE h0(256, priors.m_nTissuse);


	VECTOR_TYPE _mus;
	MAT_TYPE sigmas;
	AffinePriors(PriorsType_mni, _mus, sigmas);

	VECTOR_TYPE mus; mus.resize(6 + _mus.size());
	mus.segment(0, 6) = VECTOR_TYPE::Zero(6);
	mus.segment(6, _mus.size()) = _mus;

	MAT_TYPE Alpha0(12, 12);
	Alpha0.block(0, 0, 6, 6) = MAT_TYPE::Identity(6, 6)*0.00001;
	Alpha0.block(6, 0, 6, 6) = MAT_TYPE::Zero(6, 6);
	Alpha0.block(0, 6, 6, 6) = MAT_TYPE::Zero(6, 6);
	Alpha0.block(6, 6, 6, 6) = sigmas;
	Alpha0 *= regdata.m_FudgeFactor;
	
	VECTOR_TYPE sol = Matrix2Param(M);
	cout << "sol:\n" << sol;

	double ll = -DBL_MAX;

	VECTOR_TYPE krn = SmoothKernel(4, -256, 256, KernelSmoothMode_nearest, eps);
	
	MAT_TYPE h1 = MAT_TYPE::Ones(256, priors.m_nTissuse);

	int nsearch = 12;
	CRITICAL_SECTION CS;
	InitializeCriticalSection(&CS);

	VECTOR_TYPE dsol;
	MAT_TYPE T;
	for (int iter = 1; iter < 200; iter++)
	{
		double stepsize = 1;
		double ll1 = 0;
		VECTOR_TYPE sol1;
		MAT_TYPE y1a, y2a, y3a;
		auto &x1 = regdata.x1;
		auto &x2 = regdata.x2;
		auto &x3 = regdata.x3;
		for (int search = 1; search <= nsearch; search++)
		{
			if (iter > 1)
				sol1 = sol - stepsize * dsol;
			else
				sol1 = sol;

			double penalty = 0.5*(sol1 - mus).transpose() *Alpha0*(sol1 - mus);
			M = Param2Matrix(sol1);
			T = priors.m_TPMMat.inverse()*M*regdata.m_MG;
			int ncol = x1.cols();
			int nrow = x1.rows();
			y1a = x1 * T(0, 0) + x2 * T(0, 1) + MAT_TYPE::Constant(nrow, ncol, T(0, 3));
			y2a = x1 * T(1, 0) + x2 * T(1, 1) + MAT_TYPE::Constant(nrow, ncol, T(1, 3));
			y3a = x1 * T(2, 0) + x2 * T(2, 1) + MAT_TYPE::Constant(nrow, ncol, T(2, 3));
			
#pragma omp parallel for schedule(static)
			for (int i = 0; i < x3.size(); i++)
			{
				sRegBuffer & bf = regdata.m_Buffers[i];
				if (!bf.m_msk.HasMsk()) continue;
				auto y1 = bf.m_msk.MaskMat(y1a);
				y1 += VECTOR_TYPE::Constant(y1.size(), T(0, 2)*x3(i));
				auto y2 = bf.m_msk.MaskMat(y2a);
				y2 += VECTOR_TYPE::Constant(y2.size(), T(1, 2)*x3(i));
				auto y3 = bf.m_msk.MaskMat(y3a);
				y3 += VECTOR_TYPE::Constant(y3.size(), T(2, 2)*x3(i));

				bf.m_mskz = sMask(y3.size());
				for (int iy = 0; iy < y3.size(); iy++)
				{
					bf.m_mskz.SetMaskAndNM(iy, y3[iy] >= 1);
				}
				bf.m_mskz.Mask3Vector(y1, y2, y3);

				priors.SamplePriors(y1, y2, y3, bf.m_b);
			}
			
			for (int subit = 1; subit <= 32; subit++)
			{
				h0 = MAT_TYPE::Constant(256, ntissuse, eps);
				MAT_TYPE epsv = MAT_TYPE::Constant(h0.rows(), h0.cols(), eps);
				double ll0 = 0;
				if (subit % 4 == 0)
				{
					ll0 = ll1;
					ll1 = 0;
				}
				
#pragma omp parallel for schedule(static)
				for (int islice = 0; islice < x3.size(); islice++)
				{
					sRegBuffer & bf = regdata.m_Buffers[islice];
					if (!bf.m_msk.HasMsk() || !bf.m_mskz.HasMsk())
						continue;
					
					VECTOR_TYPE gm = bf.m_mskz.MaskVector(bf.m_g);
					gm += VECTOR_TYPE::Ones(gm.size());
					
					MAT_TYPE q = MAT_TYPE::Zero(gm.size(), ntissuse);

					for (int k = 0; k < ntissuse; k++)
					{
						auto &bk = bf.m_b[k];

						for (int ig = 0; ig < gm.size(); ig++)
						{
							int gmi = gm(ig);
							q(ig, k) = h1(gmi - 1, k) * bk(ig);
						}
					}					
					VECTOR_TYPE sq(gm.size());
					for (int irow = 0; irow < q.rows(); irow++)
					{
						sq(irow) = q.row(irow).sum() + eps;
					}
					VECTOR_TYPES accs(ntissuse);					
					for (int k = 0; k < ntissuse; k++)
					{
						VECTOR_TYPE qk = q.col(k).array() / sq.array();										
						accs[k] = Accumarray(gm, qk, 256);
					}
					{
						CriticalSectionScoper scoper(CS);
						if (subit % 4 == 0)
						{
							ll1 = ll1 + sq.array().log().sum();
						}
						//�ۼ���h0��״ͼ
						for (int k = 0; k < ntissuse; k++)
						{
							h0.col(k) += accs[k];
						}
					}
				}//end of slice

				MAT_TYPE tocov = (h0 + epsv) / h0.sum();
				MAT_TYPE filter(krn.size(), 1);
				filter.col(0) = krn;
				h1 = Conv2<MAT_TYPE>(tocov, filter, ConvMode_Same);			
				MAT_TYPE summul = h1.rowwise().sum()* h1.colwise().sum();
				h1 = h1.array() / summul.array(); //Div(h1, summul);
				if (subit % 4 == 0)
				{
					if ((ll1 - ll0) / h0.sum() < 1e-5)
						break;
				}
			}//end of subit

			for (auto &bf : regdata.m_Buffers)
			{
				bf.m_b.clear();
				bf.m_mskz = sMask();
			}
			double ssh = h0.sum();
			//MatLog(h1).array()

			ll1 = ((h0.array() * h1.array().log()).sum() - penalty) / ssh / log(2);

			if (iter == 1)
				break;
			if (abs(ll1 - ll) < 1e-4)
			{
				cout << M << endl;
				return M;
			}
			if (ll1 < ll)
			{
				stepsize = stepsize * 0.5;
				if (search == nsearch)
				{
					cout << M << endl;
					return M;
				}
			}
			else
				break;
		}//end search

		ll = ll1;
		sol = sol1;
		MAT_TYPE Alpha = MAT_TYPE::Zero(12, 12);
		VECTOR_TYPE Beta = VECTOR_TYPE::Zero(12);
		int nslice = regdata.x3.size();
#pragma omp parallel for schedule(static)
		for (int islice = 0; islice < nslice; islice++)
		{
			sRegBuffer & bf = regdata.m_Buffers[islice];
			if (!bf.m_msk.HasMsk()) continue;
			VECTOR_TYPE gi = bf.m_g + VECTOR_TYPE::Constant(bf.m_g.size(), 1.0);

			VECTOR_TYPE y1, y2, y3;
			bf.m_msk.Mask3Mat(y1a, y2a, y3a, y1, y2, y3);

			//bf.m_msk.Mask3Vector(y1,y2,y3);
			y1 += VECTOR_TYPE::Constant(y1.size(), T(0, 2)*x3(islice));
			y2 += VECTOR_TYPE::Constant(y2.size(), T(1, 2)*x3(islice));
			y3 += VECTOR_TYPE::Constant(y3.size(), T(2, 2)*x3(islice));

			sMask msk(y3.size());
			msk.MaskOnNoLess(1, y3);
			msk.Mask3Vector(y1, y2, y3);
			gi = msk.MaskVector(gi);
			int nz = y1.size();
			if (nz)
			{
				VECTOR_TYPE mi = VECTOR_TYPE::Constant(nz, eps);
				VECTOR_TYPE dmi1 = VECTOR_TYPE::Zero(nz);
				VECTOR_TYPE dmi2 = VECTOR_TYPE::Zero(nz);
				VECTOR_TYPE dmi3 = VECTOR_TYPE::Zero(nz);
				VECTOR_TYPES b, db1, db2, db3;
				priors.SamplePriorsDerive(y1, y2, y3, b, db1, db2, db3);

				for (int k = 0; k < ntissuse; k++)
				{
					VECTOR_TYPE tmp;
					VECTOR_TYPE ck = h1.col(k);
					AccessByIndex<VECTOR_TYPE, VECTOR_TYPE>(gi, ck, tmp);

					mi = mi.array() + tmp.array()*b[k].array();
					dmi1 = dmi1.array() + tmp.array()*db1[k].array();
					dmi2 = dmi2.array() + tmp.array()*db2[k].array();
					dmi3 = dmi3.array() + tmp.array()*db3[k].array();
				}

				dmi1 = dmi1.array() / mi.array();
				dmi2 = dmi2.array() / mi.array();
				dmi3 = dmi3.array() / mi.array();

				VECTOR_TYPE x1m = bf.m_msk.MaskMat(regdata.x1);
				x1m = msk.MaskVector(x1m);
				VECTOR_TYPE x2m = bf.m_msk.MaskMat(regdata.x2);
				x2m = msk.MaskVector(x2m);
				double x3m = regdata.x3(islice);

				MAT_TYPE A(nz, 12);
				A.col(0) = dmi1.array()*x1m.array();
				A.col(1) = dmi2.array()*x1m.array();
				A.col(2) = dmi3.array()*x1m.array();

				A.col(3) = dmi1.array()*x2m.array();
				A.col(4) = dmi2.array()*x2m.array();
				A.col(5) = dmi3.array()*x2m.array();

				A.col(6) = dmi1 * x3m;
				A.col(7) = dmi2 * x3m;
				A.col(8) = dmi3 * x3m;

				A.col(9) = dmi1;
				A.col(10) = dmi2;
				A.col(11) = dmi3;

				MAT_TYPE ATA = A.transpose()*A;
				{
					CriticalSectionScoper scoper(CS);
					Alpha += ATA;
					Beta -= A.colwise().sum();//SumMat(A, true);
				}
			}
		}// end of slice
		MAT_TYPE R = Derivs(priors.m_TPMMat, sol, regdata.m_MG);
		Alpha = R.transpose()*Alpha *R;
		Beta = R.transpose()*Beta;
		dsol = ((Alpha + Alpha0).inverse()*(Beta + Alpha0 * (sol - mus)));
		cout << "iter:" << iter << " dsol:" << dsol << endl;		

	}// end of iter	

	M = Param2Matrix(sol);	

	return M;
}