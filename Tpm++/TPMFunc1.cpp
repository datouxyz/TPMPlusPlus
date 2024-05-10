#include "stdafx.h"
#include <Eigen/Dense>
//#include <unsupported/Eigen/MatrixFunctions>
#include "TissueProblityMap.h"
#include "EigenMatFunc.h"
#include "TPMInternalFuncs.h"



MAT_TYPE Param2Matrix(VECTOR_TYPE & param)
{
	VECTOR_TYPE trans(3);
	trans(0) = param(0);
	trans(1) = param(1);
	trans(2) = param(2);

	MAT_TYPE T = MAT_TYPE::Zero(3, 3);
	T(1, 0) = -param(3);
	T(2, 0) = -param(4);
	T(2, 1) = -param(5);
	//cout << "T: " << T<<endl;

	MAT_TYPE R = Expm(T - T.transpose());// .exp();
	//cout << "R: " << R<<endl;

	T = MAT_TYPE::Zero(3, 3);
	T(0, 0) = param(6);
	T(0, 1) = param(7);
	T(0, 2) = param(8);
	T(1, 1) = param(9);
	T(1, 2) = param(10);
	T(2, 2) = param(11);
	//cout << "T:" << T << endl;
	MAT_TYPE TDiagnal = MAT_TYPE::Zero(3, 3);
	TDiagnal(0, 0) = T(0, 0); TDiagnal(1, 1) = T(1, 1); TDiagnal(2, 2) = T(2, 2);
	//TD.asDiagonal();
	//cout <<"TD:"<< TD<<endl;

	MAT_TYPE V = Expm(T + T.transpose() - TDiagnal);// .exp();
	//cout <<"V:"<< V<<endl;

	MAT_TYPE M = MAT_TYPE::Zero(4, 4);
	M.block(0, 0, 3, 3) = V * R;
	M(0, 3) = trans(0); M(1, 3) = trans(1); M(2, 3) = trans(2);
	M(3, 0) = 0; M(3, 1) = 0; M(3, 2) = 0; M(3, 3) = 1;
	//cout << M << endl;
	return M;
}

VECTOR_TYPE Matrix2Param(MAT_TYPE &M)
{
	//% Polar decomposition parameterisation of affine transform,
	//% based on matrix logs
	MAT_TYPE J = M.block(0, 0, 3, 3);

	MAT_TYPE V = Sqrtm((MAT_TYPE)(J * J.transpose()));// .sqrt();
	MAT_TYPE R = V.inverse()*J;
	//cout<<"J:"<<endl << J <<endl<<"V:"<< V<< endl<<"R:"<<R<< endl;

	MAT_TYPE LV = Logm(V);// .log();
	MAT_TYPE LR = -Logm(R);// .log();
	/*cout <<"LV:"<<endl<< LV << endl;
	cout <<"LR:"<<endl<< LR << endl;*/
	VECTOR_TYPE P = VECTOR_TYPE::Zero(12);
	//P.segment(0, 2) = M.block(3, 0, 1, 3);
	P(0) = M(0, 3);
	P(1) = M(1, 3);
	P(2) = M(2, 3);

	P(3) = LR(1, 0);
	P(4) = LR(2, 0);
	P(5) = LR(2, 1);

	P(6) = LV(0, 0);
	P(7) = LV(0, 1);
	P(8) = LV(0, 2);
	P(9) = LV(1, 1);
	P(10) = LV(1, 2);
	P(11) = LV(2, 2);
	//cout << P << endl;
	//cout << LV(4) << " " << LV(5) << " " << LV(6) << " " << LV(7) << " " << LV(8) << endl;
	//{, LR(0, 2), LR(1, 2) };

	return P;
}

VECTOR_TYPE RandnVecConst(int N)
{
	static int inumber = 0;
	VECTOR_TYPE ret(N);
	for (int i = 0; i < N; i++)
	{
		ret(i) = randNumber1000[(inumber + i) % 1000];
	}

	inumber += N;
	cout << "randn:" << ret << endl;
	return ret;
}


void ReFillMix(IN OUT vector<MAT_TYPE> &vr, IN OUT VECTOR_TYPE & mg, IN OUT  MAT_TYPE & mn, IN VECTOR_TYPE &lkp, IN sEstimateParam & eparam, int N)
{
	auto mn1 = mn;
	auto vr1 = vr;

	mg = VECTOR_TYPE::Constant(eparam.K, 1.0 / eparam.K);
	mn = MAT_TYPE::Ones(N, eparam.K);
	vr = vector<MAT_TYPE>(eparam.K, MAT_TYPE::Zero(N, N));
	for (int k1 = 0; k1 < eparam.Kb; k1++)
	{
		IntVEC vkk = FindInVec(lkp, k1 + 1);
		int kk = vkk.size();
		double w = 1.0 / (1 + exp(-(kk - 1)*0.25)) - 0.5;
		cout << "w: " << w;
		VECTOR_TYPE sqvr1k1 = Sqrtm(vr1[k1]);// .sqrt();

		MAT_TYPE vrmk = vr1[k1].array() *(1.0 - w);
		for (auto vk : vkk)
		{			
			mn.col(vk) = sqvr1k1 * RandnVecConst(N)*w + mn1.col(k1);			
			vr[vk] = vrmk;
			mg[vk] = 1.0 / kk;
		}
	}

	cout << "ReFillMix mg" << mg << endl;
	//cout << "ReFillMix vr" << vr << endl;
	cout << "ReFillMix mn" << mn << endl;
	cout << "lkp" << lkp.transpose() << endl;
	cout << "ReFillMix mn1" << mn1 << endl;
}