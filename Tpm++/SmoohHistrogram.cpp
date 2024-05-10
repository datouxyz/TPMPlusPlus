#include "stdafx.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
//#include <Eigen/src/SparseCholesky/SimplicialCholesky.h>
#include "ParamaterDef.h"
#include "EigenMatFunc.h"
//% Smooth a histogram
//% FORMAT[sig, alpha] = spm_smohist(t, lam)
//% t - a column vector, or matrix of column vectors containing
//%         histogram data.
//% lam - regularisation parameter, or vector of regularisation
//%         parameters for each column of t.
//% sig - smoothed probability density function(s)
//% (columns sum to 1).
//% alpha - logarithm of sig.
void SmoothHistrogram(MAT_TYPE & t0, VECTOR_TYPE &lam, OUT MAT_TYPE &sigOut, OUT MAT_TYPE & alpha)
{
	/*cout << t0 <<endl;
	cout << lam << endl;
*/
	int n = t0.rows();
	sigOut = MAT_TYPE::Zero(t0.rows(), t0.cols());
	//% Regularisation
	MAT_TYPE v(1, 3); v << -1, 2, -1;
	MAT_TYPE B = RepMat(v, n, 1);
	MAT_SPARSE_TYPE G0;
	VECTOR_TYPE dd(3); dd << -1, 0, 1;
	Spdiags(B, dd, n, n, G0);
	G0.coeffRef(0, 0) = 1;
	G0.coeffRef(G0.rows() - 1, G0.cols() - 1) = 1;
	G0 = G0.transpose()*G0;

	//Write2DMATDouble("e:\\go.mat","go_",G0.toDense() );
	//% Prevent over / underflow
	double constr = log(DBL_MAX) - 1;

	//Eigen::SimplicialLDLT<MAT_SPARSE_TYPE> solver;

	Eigen::SimplicialLLT<MAT_SPARSE_TYPE> solver;
	
	int K = t0.cols();
	for (int k = 0; k < K; k++)
	{
		VECTOR_TYPE t = t0.col(k).array() + 2.0 * t0.col(k).sum() / DBL_MAX + DBL_EPSILON;
		//Write2DMATDouble("e:\\t.mat", "t_", t);
		//t = t0(:, k) + 2 * sum(t0(:, k)) / realmax + eps;
		VECTOR_TYPE d(3); d << -1, 0, -1;
		MAT_SPARSE_TYPE G = G0 * lam(k);
		VECTOR_TYPE am = alpha.col(k);
		VECTOR_TYPE sig = am.array().exp();
		sig = sig / sig.sum();
		MAT_SPARSE_TYPE L;
		MAT_TYPE B = MAT_TYPE::Ones(n, 1).array()*(t.array().sum() + lam(k))*1e-8;
		//VECTOR_TYPE dd =;
		Spdiags(B, VECTOR_TYPE::Zero(1), n, n, L);
		//Write2DMATDouble("e:\\L.mat","L_",L.toDense());
		//L = spdiags(ones(n, 1)*(sum(t) + lam(k))*1e-8, 0, n, n);
		double ll = DBL_MAX;
		double st = t.sum();
		for (int iter = 1; iter <= 60; iter++)
		{
			//cout << "iter" << iter << endl;
			VECTOR_TYPE gr = st * sig - t + G * am;

			MAT_TYPE vm(sig.size(), 1);
			vm.col(0) = sig.array()* st + gr.array().abs()*1e-8;;
			MAT_SPARSE_TYPE W;
			Spdiags(vm, VECTOR_TYPE::Zero(1), n, n, W);
			MAT_SPARSE_TYPE H = W + G + L;
			//MAT_TYPE HDense = H.toDense();
			//VECTOR_TYPE da = HDense.ldlt()*gr;
			VECTOR_TYPE  da;
			{
				//CTimeCounter timer("llt");
				//da = HDense.llt().solve(gr);
				solver.compute(H);
				da = solver.solve(gr);
				//da = H.ldlt().solve(gr);
			}
			//Write2DMATDouble("e:\\da.mat", "da_", da);

			VECTOR_TYPE an = am - da;
			double sc = 1.0;
			while (an.array().abs().maxCoeff() > constr)
			{
				sc *= 0.5;
				an = am - sc * da;
			}
			am = an.array() - (an.array().maxCoeff() + an.array().minCoeff()) / 2.0;
			am = am.array().max(-constr);
			am = am.array().min(constr);

			//Write2DMATDouble("e:\\am.mat", "am_", am);

			sig = am.array().exp();
			sig = sig.array() / sig.array().sum();
			double oll = 0;
			if (iter % 4 == 0)
			{
				oll = ll;
				ll = -((sig.array() + 1e-9).log() * t.array()).sum() + 0.5* am.transpose()*G*am;
				/*cout << " ll " << ll << " oll " << oll << endl;
				printf("%lf_%lf %lf \n", ll, oll, st);*/
				if (oll - ll < st * 1e-5)
					break;
			}			
		}
		alpha.col(k) = am;
		sigOut.col(k) = sig;		
	}
}

void SmoothHistrogram(MAT_TYPE & t0, VECTOR_TYPE &lam, OUT MAT_TYPE &sigOut)
{
	int n = t0.rows();

	MAT_TYPE alpha = log(t0.array() + 1);

	SmoothHistrogram(t0, lam, sigOut, alpha);
}

void SmoothHistrogram(MAT_TYPE & t0, OUT MAT_TYPE &sigOut)
{
	int n = t0.rows();
	VECTOR_TYPE lam = VECTOR_TYPE::Zero(t0.cols());
	VECTOR_TYPE x = BuildVector<VECTOR_TYPE, double>(1, n, 1);

	for (int k = 0; k < t0.cols(); k++)
	{
		VECTOR_TYPE t = t0.col(k).array() + DBL_EPSILON;
		double mu = (t.array()*x.array()).sum() / t.sum();
		VECTOR_TYPE vr = t.array()*(x.array() - mu).square() / t.sum();
		lam(k) = vr(0);
	}
	MAT_TYPE alpha = log(t0.array() + 1);
	SmoothHistrogram(t0, lam, sigOut, alpha);
}


