#include "stdafx.h"
#include "ParamaterDef.h"
#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>
#include "EigenMatFunc.h"

#define _pi_double 3.141592653589793
MAT_TYPE DCT_MTX(int N, int K, VECTOR_TYPE &n)
{
	MAT_TYPE C = MAT_TYPE::Zero(n.size(), K);
	C.col(0) = VECTOR_TYPE::Ones(n.size()) / sqrt(N);
	for (int k = 2; k <= K; k++)
	{
		C.col(k-1) = sqrt(2.0 / N)*cos(_pi_double*(2 * n.array() + 1)*(k-1) / (2 * N));
	}
	return C;
}


MAT_TYPE DCT_MTX(int N)
{
	int K = N;
	VECTOR_TYPE n =BuildVector<VECTOR_TYPE, double>(0,N-1,1);
	return DCT_MTX(N, K, n);	
}

MAT_TYPE DCT_MTX(int N,int K)
{
	VECTOR_TYPE n = BuildVector<VECTOR_TYPE, double>(0, N - 1, 1);
	return DCT_MTX(N, K, n);	
}
MAT_TYPE DCT_MTX_DIFF(int N, int K, VECTOR_TYPE &n)
{
	MAT_TYPE C = MAT_TYPE::Zero(n.size(), K);
	//C.col(0) = VECTOR_TYPE::Ones(n.size()) / sqrt(N);
	double A = -sqrt(2.0)*sqrt(1.0 / N);
	for (int k = 2; k <= K; k++)
	{
		//C(:, k) = -2 ^ (1 / 2)*(1 / N)^(1 / 2)*
		//sin(1 / 2 * pi*(2 * n*k - 2 * n + k - 1) / N)*pi*(k - 1) / N;
		C.col(k-1) = A * sin(1 / 2.0*_pi_double*(2 * n.array()*k - 2 * n.array() + k-1)/ N)*_pi_double*(k-1) / N;
	}
	return C;
}
MAT_TYPE DCT_MTX_DIFF(int N)
{
	int K = N;
	VECTOR_TYPE n = BuildVector<VECTOR_TYPE, double>(0, N - 1, 1);
	return DCT_MTX_DIFF(N,K,n);
}
MAT_TYPE DCT_MTX_DIFF(int N,int K)
{
	VECTOR_TYPE n = BuildVector<VECTOR_TYPE, double>(0, N - 1, 1);
	return DCT_MTX_DIFF(N, K, n);
}

MAT_TYPE DCT_MTX_DIFF2(int N, int K, VECTOR_TYPE &n)
{
	MAT_TYPE C = MAT_TYPE::Zero(n.size(), K);
	//C.col(0) = VECTOR_TYPE::Ones(n.size()) / sqrt(N);
	double A = -sqrt(2.0)*sqrt(1.0 / N);
	for (int k = 2; k <= K; k++)
	{
//C(:, k) = -2 ^ (1 / 2)*(1 / N) ^ (1 / 2)*cos(1 / 2 * pi*(2 * n + 1)*(k - 1) / N)
//*pi ^ 2 * (k - 1) ^ 2 / N ^ 2;

		C.col(k - 1) = A * cos(0.5*_pi_double*(2 * n.array() + 1)*(k - 1) / N)*
			pow(_pi_double,2) * pow(k - 1,2) / pow(N, 2);
	}
	return C;
}

MAT_TYPE DCT_MTX_DIFF2(int N)
{
	int K = N;
	VECTOR_TYPE n = BuildVector<VECTOR_TYPE, double>(0, N - 1, 1);
	return DCT_MTX_DIFF2(N, K, n);
}
MAT_TYPE DCT_MTX_DIFF2(int N, int K)
{
	VECTOR_TYPE n = BuildVector<VECTOR_TYPE, double>(0, N - 1, 1);
	return DCT_MTX_DIFF2(N, K, n);
}
