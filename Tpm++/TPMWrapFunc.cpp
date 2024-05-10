#include "stdafx.h"
#include <Eigen/Dense>
#include <Eigen/KroneckerProduct>
#include "TissueProblityMap.h"
#include "EigenMatFunc.h"
#include "llVar.h"
#include "TPMInternalFuncs.h"

void UpdateTWrap(IN OUT CVariableLL & ll ,IN MAT_TYPE & x0, IN MAT_TYPE & y0, IN VECTOR_TYPE & z0, IN vector<sWrapBUF> & WrapBufs,
	IN VECTOR_TYPE & param,IN OUT sEstimateParam & eparam, IN VECTOR_TYPE &WrapPrior, IN CDGMulDimData & TWrap, IN sTPMPriors & tpm,
	IN MAT_TYPE MT, IN MAT_TYPE M,int iter,bool bDecreaseParamByIter, IN VECTOR_TYPE &sk)
{
	
	CRITICAL_SECTION cs; InitializeCriticalSection(&cs);
	
	
	double oll = ll.GetLL();
	DimVec d = { (int)x0.rows(),(int)x0.cols(),(int)z0.size() };
	for (int subiter = 1; subiter <= 3; subiter++)
	{
		CDGMulDimData Alpha4D({ (int)x0.rows(),(int)x0.cols(),(int)z0.size(),6 },DGDataType_Float,DataMajor_FStyle_ColumnMajor);
		CDGMulDimData Beta4D({ (int)x0.rows(),(int)x0.cols(),(int)z0.size(),3 }, DGDataType_Float, DataMajor_FStyle_ColumnMajor);
		
		CDGMulDimData TWrapD; TWrap.ToDouble(TWrapD);
		#pragma omp parallel for schedule(dynamic , 1)
		for (int z = 0; z < z0.size(); z++)
		{
			sWrapBUF & buf = WrapBufs[z];
			if (buf.msk.m_NM <= 0) continue;
			VECTOR_TYPE x1, y1, z1;
			defs(TWrapD, z, x0, y0, z0, M, buf.msk, x1, y1, z1);
			VECTOR_TYPES b, db1, db2, db3;
			tpm.SamplePriorsDerive(x1, y1, z1, b, db1, db2, db3);

			VECTOR_TYPE s = VECTOR_TYPE::Zero(b[0].size());
			VECTOR_TYPE ds1 = VECTOR_TYPE::Zero(b[0].size());
			VECTOR_TYPE ds2 = VECTOR_TYPE::Zero(b[0].size());
			VECTOR_TYPE ds3 = VECTOR_TYPE::Zero(b[0].size());
			for (int k1 = 0; k1 < eparam.Kb; k1++)
			{
				b[k1] = WrapPrior(k1)*b[k1];
				db1[k1] = WrapPrior(k1)*db1[k1];
				db2[k1] = WrapPrior(k1)*db2[k1];
				db3[k1] = WrapPrior(k1)*db3[k1];
				s = s + b[k1];
				ds1 = ds1 + db1[k1];
				ds2 = ds2 + db2[k1];
				ds3 = ds3 + db3[k1];
			}
			for (int k1 = 0; k1 < eparam.Kb; k1++)
			{
				b[k1] = b[k1].array() / s.array();
				db1[k1] = (db1[k1].array() - b[k1].array()*ds1.array()) / s.array();
				db2[k1] = (db2[k1].array() - b[k1].array()*ds2.array()) / s.array();
				db3[k1] = (db3[k1].array() - b[k1].array()*ds3.array()) / s.array();
			}

			/*% Rotate gradients(according to initial affine registration) and
				% compute the sums of the tpm and its gradients, times the likelihoods
				% (from buf.dat).*/
			VECTOR_TYPE p = VECTOR_TYPE::Zero(buf.msk.m_NM).array() + eparam.eps;
			MAT_TYPE dp1, dp2, dp3;

			{
				VECTOR_TYPE dp1v = VECTOR_TYPE::Zero(buf.msk.m_NM);
				VECTOR_TYPE dp2v = VECTOR_TYPE::Zero(buf.msk.m_NM);
				VECTOR_TYPE dp3v = VECTOR_TYPE::Zero(buf.msk.m_NM);

				MAT_TYPE MM = M * MT;

				for (int k1 = 0; k1 < eparam.Kb; k1++)
				{
					VECTOR_TYPE pp = AsDouble((VECTOR_TYPE_F)buf.dat.col(k1));
					p += (VECTOR_TYPE)(pp.array()*b[k1].array());
					dp1v += (VECTOR_TYPE)(pp.array()*(MM(0, 0)*db1[k1] + MM(1, 0)*db2[k1] + MM(2, 0)*db3[k1]).array());
					dp2v += (VECTOR_TYPE)(pp.array()*(MM(0, 1)*db1[k1] + MM(1, 1)*db2[k1] + MM(2, 1)*db3[k1]).array());
					dp3v += (VECTOR_TYPE)(pp.array()*(MM(0, 2)*db1[k1] + MM(1, 2)*db2[k1] + MM(2, 2)*db3[k1]).array());
				}
				MAT_TYPE tmp = MAT_TYPE::Zero(d[0], d[1]);
				buf.msk.MaskAssign(tmp, (VECTOR_TYPE)(dp1v.array() / p.array())); dp1 = tmp;
				buf.msk.MaskAssign(tmp, (VECTOR_TYPE)(dp2v.array() / p.array())); dp2 = tmp;
				buf.msk.MaskAssign(tmp, (VECTOR_TYPE)(dp3v.array() / p.array())); dp3 = tmp;
			}
			Eigen::Map<MAT_TYPE_F>((float*)Beta4D.SelMatByDimXY({ z,0 }), x0.rows(), x0.cols()) = -AsSingle(dp1);
			Eigen::Map<MAT_TYPE_F>((float*)Beta4D.SelMatByDimXY({ z,1 }), x0.rows(), x0.cols()) = -AsSingle(dp2);
			Eigen::Map<MAT_TYPE_F>((float*)Beta4D.SelMatByDimXY({ z,2 }), x0.rows(), x0.cols()) = -AsSingle(dp3);

			Eigen::Map<MAT_TYPE_F>((float*)Alpha4D.SelMatByDimXY({ z,0 }), x0.rows(), x0.cols()) = AsSingle((MAT_TYPE)(dp1.array().square()));
			Eigen::Map<MAT_TYPE_F>((float*)Alpha4D.SelMatByDimXY({ z,1 }), x0.rows(), x0.cols()) = AsSingle((MAT_TYPE)(dp2.array().square()));
			Eigen::Map<MAT_TYPE_F>((float*)Alpha4D.SelMatByDimXY({ z,2 }), x0.rows(), x0.cols()) = AsSingle((MAT_TYPE)(dp3.array().square()));
			Eigen::Map<MAT_TYPE_F>((float*)Alpha4D.SelMatByDimXY({ z,3 }), x0.rows(), x0.cols()) = AsSingle((MAT_TYPE)(dp1.array()*dp2.array()));
			Eigen::Map<MAT_TYPE_F>((float*)Alpha4D.SelMatByDimXY({ z,4 }), x0.rows(), x0.cols()) = AsSingle((MAT_TYPE)(dp1.array()*dp3.array()));
			Eigen::Map<MAT_TYPE_F>((float*)Alpha4D.SelMatByDimXY({ z,5 }), x0.rows(), x0.cols()) = AsSingle((MAT_TYPE)(dp2.array()*dp3.array()));


		}// end of z loop

		VECTOR_TYPE prm(10);

		if (bDecreaseParamByIter)
		{
			double weight = 0.0;
			switch (iter)
			{
			case 1:weight = 256; break;
			case 2:weight = 128; break;
			case 3:weight = 64; break;
			case 4:weight = 32; break;
			case 5:weight = 16; break;
			case 6:weight = 8; break;
			case 7:weight = 4; break;
			case 8:weight = 2; break;
			default:weight = 1;
			}
			prm.segment(0, 3) = param.segment(0, 3);
			prm.segment(3, 5) = param.segment(3, 5) *weight;
		}
		else
		{
			prm.segment(0, 8) = param;
		}

		DimVec sk4dim = { 1,1,1,3 };
		{
			float skinv[3] = { 1.0f / sk[0],1.0f / sk[1],1.0f / sk[2] };
			CDGMulDimData divsk4(sk4dim, DGDataType_Float, DataMajor_FStyle_ColumnMajor, skinv);

			CDGMulDimData bsxres;
			TWrap.BsxFuncHighDimDup(BsxOperator_times, divsk4, 1, bsxres);
			CDGMulDimData diffeo;
			Vel2Mom(prm.data(), bsxres, diffeo);

			Beta4D.AddOn<float>(diffeo);
			//WriteNDMatFloat("e:\\Beta4.mat", "Beta4_", Beta4D);
		}

		CDGMulDimData Update4D;
		{
			prm.segment(8, 2) << 2, 2;
			CDGMulDimData X;
			Fmg(Alpha4D, Beta4D, prm.data(), X);

			float skf[3] = { sk[0],sk[1],sk[2] };
			CDGMulDimData sk4(sk4dim, DGDataType_Float, DataMajor_FStyle_ColumnMajor, skf);
			X.BsxFuncHighDimDup(BsxOperator_times, sk4, 1, Update4D);

			//WriteNDMatFloat("e:\\Update4D.mat", "Update4D_", Update4D);

		}
		double  armijo = 1.0;
		for (int line_search = 1; line_search <= 12; line_search++)
		{
			CDGMulDimData Twarp1D;
			float llr1 = 0;
			{
				CDGMulDimData U = Update4D;
				U.MulOn<float>(armijo);
				CDGMulDimData Twarp1 = TWrap;
				Twarp1.SubOn<float>(U);
				llr1 = ComputeObjFunction(sk, Twarp1, prm);

				Twarp1.ToDouble(Twarp1D);
			}			
			//WriteNDMatFloat("e:\\Twarp1.mat", "Twarp1_", Twarp1);
						
			float ll1 = llr1 + eparam.llrb + eparam.ll_const;

			#pragma omp parallel for schedule(dynamic , 1)
			for (int z = 0; z < WrapBufs.size(); z++)
			{
				sWrapBUF & buf = WrapBufs[z];
				if (buf.msk.m_NM <= 0) continue;

				VECTOR_TYPES b;
				{
					VECTOR_TYPE x1, y1, z1;
					defs(Twarp1D, z, x0, y0, z0, M, buf.msk, x1, y1, z1);
					tpm.SamplePriors(x1, y1, z1, b);
				}
				VECTOR_TYPE s = VECTOR_TYPE::Zero(b[0].size());
				//for k1 = 1:Kb, b{ k1 } = b{ k1 }*wp(k1); s = s + b{ k1 }; end
				for (int k1 = 0; k1 < eparam.Kb; k1++)
				{
					b[k1] *= WrapPrior(k1);
					s += b[k1];
				}
				//for k1 = 1:Kb, b{ k1 } = b{ k1 }. / s; end
				for (int k1 = 0; k1 < eparam.Kb; k1++)
				{
					b[k1] = b[k1].array() / s.array();
				}
				VECTOR_TYPE sq = VECTOR_TYPE::Zero(buf.msk.m_NM);
				for (int k1 = 0; k1 < eparam.Kb; k1++)
				{
					sq += (VECTOR_TYPE)(AsDouble((VECTOR_TYPE_F)(buf.dat.col(k1))).array() * b[k1].array());
				}
				{
					CriticalSectionScoper scoper(cs);
					ll1 += sq.array().log().sum();
				}				
			}// endof z
			if (ll1 < ll.GetLL())
			{
				printf("Warp:\t%g\t%g\t%g :o(\t(%g)\n", ll1, llr1, eparam.llrb, armijo);
				armijo = armijo * 0.75;
			}
			else
			{
				ll.SetLL(ll1);
				eparam.llr = llr1;
				Twarp1D.ToSingle(TWrap);				
				printf("Warp:\t%g\t%g\t%g :o(\t(%g)\n", ll1, llr1, eparam.llrb, armijo);
				break;
			}
		}
		if (!((ll.GetLL() - oll) > eparam.tol1*eparam.nm))
			//% Registration no longer helping, so move on
			break;
		oll = ll.GetLL();

	}// end of subiter

	DeleteCriticalSection(&cs);
}