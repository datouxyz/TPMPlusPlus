#include "stdafx.h"
#include <Eigen/Dense>

#include "TissueProblityMap.h"
#include "EigenMatFunc.h"
#include "Conv2.h"
#include "shoot_boundary.h"
#include "llVar.h"

#include "TPMInternalFuncs.h"
sitk::Image im;

#define WriteMDToImg(MD,filepath,itkpixeltype) {sitk::Image newimg(DimToSize(MD.m_dim), itkpixeltype);\
memcpy(newimg.GetBufferAsFloat(), MD.ptr(), MD.SizeInBytes());\
sitk::WriteImage(newimg, filepath);}

void TpmWrite8(sTpmInputs & TInput, sTPMPriors & tpm, IN sWrapResults & res, sWriteImgs & sWrites)
{
	CTimeCounter w8("TpmWrite8");
	int Kb = 0;
	VECTOR_TYPE lkp;
	if (res.m_mg.size() > 0)
	{
		lkp = res.m_lkp;
		Kb = lkp.maxCoeff();
	}
	else
	{
		Kb = res.m_intensity[0].m_lik.cols(); //size(res.intensity(1).lik, 2);
	}
	sTissues &tc = TInput.m_tissues;
	ItkImageVec Images;
	for (auto & ch : TInput.m_Channels) { Images.push_back(ch.m_Img); }

	bool AnyWrapModulated = AnyWarpTissueSetting(tc, false, false, true, false);
	bool AnyWarpUnModulated = AnyWarpTissueSetting(tc, false, false, false, true);
	bool AnyTissue = AnyWarpTissueSetting(tc, true, true, true, true);

	// tpm dim
	VECTOR_TYPE dTpm(3); dTpm << tpm.m_size[0], tpm.m_size[1], tpm.m_size[2];
	DimVec DimTpm = SizeToDim(tpm.m_size);

	MAT_TYPE M1 = tpm.m_TPMMat;
	
	MAT_TYPE bb1;
	VECTOR_TYPE vx1;
	GetBoundBox(tpm.m_TPMMat, tpm.m_size, bb1, vx1);
	double vx = pow(abs(vx1.prod()), 1.0 / 3.0);
	MAT_TYPE & bb = TInput.m_bb;
	bb = bb1;
	bb.row(0) = vx * (bb.row(0) / vx).array().round();
	bb.row(1) = vx * (bb.row(1) / vx).array().round();

	cout << bb << endl;

	VECTOR_TYPE bbext = (bb.row(1) - bb.row(0));
	VECTOR_TYPE odim = (bbext.array().round() / vx).abs() + 1;
	// img dim
	DimVec img0d = SizeToDim(Images[0].GetSize());
	int volsize = img0d[0]* img0d[1]*img0d[2];
	int slicesize = img0d[0] * img0d[1];
	int slicesizefloat = slicesize * sizeof(float);
	
	int N = (int)Images.size();

	for (int n = 0; n < N; n++)
	{
		if (Images[n].GetPixelID() != sitk::sitkFloat32)
		{
			ErrOutput("TpmWrite8", "Images ����ΪFloat32");
			return;
		}
	}
	bool do_cls = false;
	bool do_defs = false;
	CDGMulDimData y; CDGMulDimData Q;
	sWrites.InitDatas(N,tc, TInput.m_Channels, TInput.m_dFieldInverse,TInput.m_dFieldForward, Images,tpm,do_cls, do_defs,Kb,y,Q);

	ChanWrites & chans  = sWrites.chans;
	sitk::Image &img0 = Images[0];
	MAT_TYPE img0mat = GetMatFromImageDirOriSpace(img0.GetDirection(), img0.GetOrigin(), img0.GetSpacing(), true, true);//
	cout << img0mat;
	MAT_TYPE x1, x2; VECTOR_TYPE x3 = BuildVector<VECTOR_TYPE, double>(1, img0d[2], 1.0);
	Grid2D<VECTOR_TYPE>(BuildVector<VECTOR_TYPE, double>(1, img0d[0], 1.0), BuildVector<VECTOR_TYPE, double>(1, img0d[1], 1.0),
		x1, x2);

	for (int n = 0; n < N; n++)
	{
		sitk::Image & img = Images[n];

		DimVec d3 = res.m_T[n].m_dim;
		d3.push_back(1);

		auto &chan = chans[n];
		chan.B3 = DCT_MTX(img0d[2], d3[2], x3);
		chan.B2 = DCT_MTX(img0d[1], d3[1], (VECTOR_TYPE)x2.row(0));
		chan.B1 = DCT_MTX(img0d[0], d3[0], (VECTOR_TYPE)x1.col(0));
		chan.T = res.m_T[n];
	}	
	
	DoubleVec prm = { 3,3,3,0,0,0 };
	
	DimVec wrapdim = { res.m_Twrap.m_dim[0],res.m_Twrap.m_dim[1],res.m_Twrap.m_dim[2] };

	vector<CDGMulDimData> Coef = { CDGMulDimData(wrapdim,DGDataType_Double,DataMajor_FStyle_ColumnMajor),
								   CDGMulDimData(wrapdim,DGDataType_Double,DataMajor_FStyle_ColumnMajor),
								   CDGMulDimData(wrapdim,DGDataType_Double,DataMajor_FStyle_ColumnMajor) };
	CDGMulDimData TWDouble;
	res.m_Twrap.ToDouble(TWDouble);
	BSplinec(prm, (double*)TWDouble.SelFrameByDim(3, { 0 }), wrapdim, Coef[0].Ptr<double>());
	BSplinec(prm, (double*)TWDouble.SelFrameByDim(3, { 1 }), wrapdim, Coef[1].Ptr<double>());
	BSplinec(prm, (double*)TWDouble.SelFrameByDim(3, { 2 }), wrapdim, Coef[2].Ptr<double>());	
	
	MAT_TYPE M = M1.inverse() * res.m_Affine * img0mat;	

	{
		for (int z = 0; z < x3.size(); z++)
		{
			vector<MAT_TYPE> cr(N);
			for (int n = 0; n < N; n++)
			{
				sChanWrite &chan = chans[n];
				float* pf = Images[n].GetBufferAsFloat();
				MAT_TYPE f, bf;
				f = AsDouble((MAT_TYPE_F)Eigen::Map<MAT_TYPE_F>(&pf[slicesize *z], img0d[0], img0d[1]));
				bf = transf(chan.B1, chan.B2, chan.B3.row(z), chan.T).array().exp();
				cr[n] = bf.array()*f.array();

				if (!chan.m_Nc.IsEmpty())
				{
					MAT_TYPE_F crf = AsSingle(cr[n]);
					void * pNc = chan.m_Nc.RefImgData().SelFrameByDim(2, { z });
					memcpy(pNc, crf.data(), crf.size() * sizeof(float));
				}
				if (!chan.m_Nf.IsEmpty())
				{
					MAT_TYPE_F bff = AsSingle(bf);
					void * pNf = chan.m_Nf.RefImgData().SelFrameByDim(2, { z });
					memcpy(pNf, bff.data(), slicesize * sizeof(float));
				}
			}

			MAT_TYPE t1, t2, t3;
			if (do_defs)
			{
				//[t1, t2, t3] = defs(Coef, z, res.MT, prm, x1, x2, x3, M);			
				defsw(Coef, z, x1, x2, x3, M, res.m_MT, prm, t1, t2, t3);

				if (!sWrites.FieldInverse.IsEmpty())
				{
					MAT_TYPE tmp;
					tmp = (M1(0, 0)*t1 + M1(0, 1)*t2 + M1(0, 2)*t3).array() + M1(0, 3);
					memcpy(sWrites.FieldInverse.RefImgData().SelMatByDimXY({ z,0 }), AsSingle(tmp).data(), slicesizefloat);

					tmp = (M1(1, 0)*t1 + M1(1, 1)*t2 + M1(1, 2)*t3).array() + M1(1, 3);
					memcpy(sWrites.FieldInverse.RefImgData().SelMatByDimXY({ z,1 }), AsSingle(tmp).data(), slicesizefloat);

					tmp = (M1(2, 0)*t1 + M1(2, 1)*t2 + M1(2, 2)*t3).array() + M1(2, 3);
					memcpy(sWrites.FieldInverse.RefImgData().SelMatByDimXY({ z,2 }), AsSingle(tmp).data(), slicesizefloat);
				}
				if (y.SizeInBytes() > 0)
				{
					memcpy(y.SelMatByDimXY({ z,0 }), AsSingle(t1).data(), slicesizefloat);
					memcpy(y.SelMatByDimXY({ z,1 }), AsSingle(t2).data(), slicesizefloat);
					memcpy(y.SelMatByDimXY({ z,2 }), AsSingle(t3).data(), slicesizefloat);
				}
			}
			CDGMulDimData q1;
			if (do_cls)
			{
				CDGMulDimData qm({ img0d[0],img0d[1],Kb }, DGDataType_Double, DataMajor_FStyle_ColumnMajor);
				auto t1m = Eigen::Map<VECTOR_TYPE>(t1.data(), t1.size());
				auto t2m = Eigen::Map<VECTOR_TYPE>(t2.data(), t2.size());
				auto t3m = Eigen::Map<VECTOR_TYPE>(t3.data(), t3.size());
				if (res.m_mg.size() > 0)
				{
					// ��ȡtpm ,����mg,����b/s�����͡�������ع�һ��.
					auto emptybf = vector<MAT_TYPE>();
					{
						MAT_TYPE qlik;
						likelihoods(cr, emptybf, res.m_mg, res.m_mn, res.m_vr, qlik);
						q1 = CDGMulDimData({ img0d[0],img0d[1],(int)res.m_mg.size() }, DGDataType_Double, DataMajor_FStyle_ColumnMajor, qlik.data());
					}
					VECTOR_TYPES b;
					tpm.SamplePriors((VECTOR_TYPE)t1m, (VECTOR_TYPE)t2m, (VECTOR_TYPE)t3m, b);
					VECTOR_TYPE s = VECTOR_TYPE::Zero(b[0].size());
					for (int k1 = 0; k1 < Kb; k1++)
					{
						b[k1] *= res.m_WrapPrior(k1);
						s += b[k1];
					}
					for (int k1 = 0; k1 < Kb; k1++)
					{
						auto finded = FindInVec(lkp, k1 + 1);
						VECTOR_TYPE sumed = VECTOR_TYPE::Zero(slicesize);
						for (auto fi : finded)
						{
							sumed += Eigen::Map<VECTOR_TYPE>((double*)q1.SelMatByDimXY({ fi }), slicesize);
						}
						//Write2DMATDouble("e:\\Sumed.mat", "Sumed_", sumed);

						Eigen::Map<VECTOR_TYPE>((double*)qm.SelMatByDimXY({ k1 }), slicesize) =
							sumed.array()*(b[k1].array() / s.array());//b[k1].array()
					}
				}
				else
				{
					VECTOR_TYPES q;
					tpm.SamplePriors((VECTOR_TYPE)t1m, (VECTOR_TYPE)t2m, (VECTOR_TYPE)t3m, q);
					VECTOR_TYPE s = VECTOR_TYPE::Zero(q[0].size());
					for (int k1 = 0; k1 < Kb; k1++)
					{
						q[k1] *= res.m_WrapPrior(k1);
						s += q[k1];
					}
					for (int k1 = 0; k1 < Kb; k1++)
					{
						q[k1] = q[k1].array() / s.array();
					}
					//cat qs to qm

					for (int k1 = 0; k1 < Kb; k1++)
					{
						memcpy(qm.SelMatByDimXY({ k1 }), q[k1].data(), q[k1].size() * sizeof(double));
					}
					for (int n = 0; n < N; n++)
					{
						MAT_TYPE_F tmp =
							AsSingle((MAT_TYPE)(cr[n].array()*res.m_intensity[n].m_interscal(1) + res.m_intensity[n].m_interscal(0)).round());

						tmp = tmp.array().max(1).min(res.m_intensity[n].m_lik.rows());
						VECTOR_TYPE_F tmpAsind = Eigen::Map<VECTOR_TYPE_F>(tmp.data(), tmp.size());

						for (int k1 = 0; k1 < Kb; k1++)
						{
							VECTOR_TYPE_F likelihood = res.m_intensity[n].m_lik.col(k1);// (:, k1);
							Eigen::Map<MAT_TYPE>qmm((double*)qm.SelMatByDimXY({ k1 }), (int)q[0].rows(), (int)q[0].cols());

							VECTOR_TYPE_F likelihoodout;
							AccessByIndex<VECTOR_TYPE_F, VECTOR_TYPE_F>(tmpAsind, likelihood, likelihoodout);

							qmm = qmm.array()*AsDouble(likelihoodout).array();
						}
					}
				}
				for (int k1 = 0; k1 < Kb; k1++)
				{
					MAT_TYPE qmkd = Eigen::Map<MAT_TYPE>((double*)qm.SelMatByDimXY({ k1 }), img0d[0], img0d[1]);
					memcpy(Q.SelMatByDimXY({ z,k1 }), AsSingle(qmkd).data(), slicesize * sizeof(float));
				}
			}
		}// end of z		
	}
	if (do_cls)
	{	
		CDGMulDimData P;
		bool bLoadP = true;
		if (!TInput.m_Mrf)
		{			
			CDGMulDimData sQ; Q.SumOnDim<float>(3, sQ); 
			sQ.AddOn<float>(eps);
			sQ.MulOn<float>(Div255d);
			DimVec nd = sQ.m_dim; nd.push_back(1);
			sQ.ReShapeInPlace(nd);
			Q.BsxFuncLowDimDup(BsxOperator_rounded_divide, sQ, P);
			P.ToTypeInPlace<unsigned char>(DGDataType_UnSignedChar);
		}
		else
		{
			P = CDGMulDimData({ img0d[0],img0d[1],img0d[2],Kb }, DGDataType_UnSignedChar, DataMajor_FStyle_ColumnMajor);
			int nmrf_its = 10;			
			CDGMulDimData G({Kb},DGDataType_Float, DataMajor_FStyle_ColumnMajor);
			for (int k = 0; k < Kb; k++)  *(float*)G.DataPtr(k) = TInput.m_Mrf;

			MAT_TYPE m33 = img0mat.block(0, 0, 3, 3);
			VECTOR_TYPE vx2 = m33.array().pow(2).rowwise().sum();

			if(!bLoadP)
			for (int iter = 1; iter <= nmrf_its; iter++)
			{
				MarkrofField(P, Q, G, vx2);
			}
		}	
		if (TInput.m_Cleanup)
		{
			if (!bLoadP)
			if (P.m_dim[3] > 5)
			{		
				Clean_Gwc(P,TInput.m_Cleanup);			
			}			
		}
	
		for (int k1 = 0; k1 < Kb; k1++)
		{
			sWriteImgData &towrite = sWrites.tissNative[k1];
			if (!towrite.IsEmpty())
			{
				memcpy(towrite.RefImgData().ptr(), P.SelFrameByDim(3,{ k1 }), volsize);
			}
		}				
		vector<CDGMulDimData> cls(Kb);
		bool AnyWrapModulated = AnyWarpTissueSetting(tc, false, false, true, false);
		for (int k1 = 0; k1 < Kb; k1++)
		{
			if (tc[k1].m_WarpTissue_UnModulated || tc[k1].m_NativeTissue_DartelImported || AnyWrapModulated)
			{
				P.SelFrameByDim(3, { k1 }, cls[k1]);
			}			
		}
		P.clear();		
		//% Compute Voxel - to - world of "imported" image space
		MAT_TYPE m83(8, 3); m83 << 
			bb(0, 0), bb(0, 1), bb(0, 2),
			bb(1, 0), bb(0, 1), bb(0, 2),
			bb(0, 0), bb(1, 1), bb(0, 2),
			bb(1, 0), bb(1, 1), bb(0, 2),
			bb(0, 0), bb(0, 1), bb(1, 2),
			bb(1, 0), bb(0, 1), bb(1, 2),
			bb(0, 0), bb(1, 1), bb(1, 2),
			bb(1, 0), bb(1, 1), bb(1, 2);
		MAT_TYPE mm(4, 8);
		mm.block(0, 0, 3, 8) = m83.transpose();
		mm.block(3, 0, 1, 8) = MAT_TYPE::Ones(1, 8);
		MAT_TYPE vx83(8, 3); vx83 <<
			1,		 1,		  1,
			odim(0), 1,		  1,
			1,		 odim(1), 1,
			odim(0), odim(1), 1,
			1,		 1,	      odim(2),
			odim(0), 1,		  odim(2),
			1,	     odim(1), odim(2),
			odim(0), odim(1), odim(2);
		MAT_TYPE vx3(4, 8);
		vx3.block(0, 0, 3, 8) = vx83.transpose();
		vx3.block(3, 0, 1, 8) = MAT_TYPE::Ones(1, 8);
		MAT_TYPE mat = RDivide( vx3, mm);
	
		if (AnyWrapModulated)
		{
							
			VECTOR_TYPE fwhm = (vx / (img0mat.block(0, 0, 3, 3).array().pow(2).rowwise().sum().sqrt()) - 1).max(0.01);
			for (int k1 = 0; k1 < tc.size(); k1++)
			{
				sTissue & tiss = tc[k1];
				sWriteImgData & Ni = sWrites.tissImported[k1];
				if (tiss.m_NativeTissue_DartelImported)
				{
					CDGMulDimData fclsk;
					cls[k1].ToSingle(fclsk);
					CDGMulDimData tmp1;					

					decimate(fclsk, fwhm, tmp1);
					
					tmp1.MulOn<float>(Div255);
					sitk::Image newimg(DimToSize(img0d), sitk::sitkFloat32);//tmp1;
					newimg.SetOrigin(img0.GetOrigin());
					newimg.SetDirection(img0.GetDirection());
					newimg.SetSpacing(img0.GetSpacing());
					memcpy(newimg.GetBufferAsFloat(),tmp1.ptr(),tmp1.SizeInBytes());
					/*Ni.m_Img = Resample(newimg, newimg.GetSize(), sitk::Transform(), sitk::sitkNearestNeighbor,
					img0.GetOrigin(),img0.GetSpacing(), img0.GetDirection());
					*/
					Ni.m_Img = Resample(newimg, tpm.m_size,sitk::Transform() ,sitk::sitkLinear,
					tpm.m_TpmSpace.m_Ori, tpm.m_TpmSpace.m_Spacing, tpm.m_TpmSpace.m_Dir);
					Ni.m_Img.SetMetaData("description", "Imported Tissue");
				}
			}
		}

		if (AnyWrapModulated || AnyWarpUnModulated || TInput.m_dFieldForward)
		{
			M = mat.inverse()*M1;
			for (int z=0;z<y.m_dim[2];z++)
			{
				auto mapt1 = y.SelMatByDimXY_Maped<MAT_TYPE_F>({ z,0 });
				auto mapt2 = y.SelMatByDimXY_Maped<MAT_TYPE_F>({ z,1 });
				auto mapt3 = y.SelMatByDimXY_Maped<MAT_TYPE_F>({ z,2 });
				MAT_TYPE_F t1 = mapt1;
				MAT_TYPE_F t2 = mapt2;
				MAT_TYPE_F t3 = mapt3;
				mapt1 = M(0,0)*t1.array() + M(0,1)*t2.array() + M(0,2)*t3.array() + M(0,3);
				mapt2 = M(1,0)*t1.array() + M(1,1)*t2.array() + M(1,2)*t3.array() + M(1,3);
				mapt3 = M(2,0)*t1.array() + M(2,1)*t2.array() + M(2,2)*t3.array() + M(2,3);
			}
			M1 = mat;
			dTpm = odim;
		}
		
		SizeVec TpmSize = { (UINT)dTpm(0), (UINT)dTpm(1), (UINT)dTpm(2) };
		
		SizeVec dmf = { (UINT)img0d[0] , (UINT)img0d[1] , (UINT)img0d[2], (UINT)1 };
		
		SizeVec dmo = tpm.m_size;
		dmo.push_back(dmf[3]);

		int m = dmf[0] * dmf[1] * dmf[2];
		int n = dmf[3];//c.m_dim[3];
		
		if (AnyWarpUnModulated)
		{			
			for (int k1 = 0; k1 < Kb; k1++)
			{
				if (!cls[k1].bEmpty() && tc[k1].m_WarpTissue_UnModulated)
				{
					CDGMulDimData c;
					cls[k1].ToSingle(c); 
					c.MulOn<float>(Div255d);				
					
					CDGMulDimData tmp(DimTpm, c.m_type, DataMajor_FStyle_ColumnMajor);
					push(&dmo[0], m, n, y.Ptr<float>(), c.Ptr<float>(), tmp.Ptr<float>(),NULL);
					c = tmp;					
					{						
						auto & tissUnM = sWrites.tissUnModulate;

						tissUnM[k1].m_Img = sitk::Image(TpmSize, sitk::sitkFloat32);
						tpm.m_TpmSpace.SetToImage(tissUnM[k1].m_Img);
						tissUnM[k1].SetDescription( "unmodulated warped tissue class");
						memcpy(tissUnM[k1].RefImgData().BytePtr() , c.BytePtr(),c.SizeInBytes());
						double sc = abs(img0mat.block(0, 0, 3, 3).determinant() / M1.block(0, 0, 3, 3).determinant());
						tissUnM[k1].RefImgData().MulOn<float>(sc);						
					}				
				}
			}
		}

		if (AnyWrapModulated)
		{
			CDGMulDimData C;			
			C = CDGMulDimData({ (int)dTpm(0),(int)dTpm(1),(int)dTpm(2),Kb }, DGDataType_Float, DataMajor_FStyle_ColumnMajor);
			
			for (int k1 = 0; k1 < Kb; k1++)
			{
				if (!cls[k1].bEmpty() && AnyWrapModulated)
				{
					CDGMulDimData c;
					cls[k1].ToSingle(c);
					c.MulOn<float>(Div255d);
					CDGMulDimData w(DimTpm, c.m_type, DataMajor_FStyle_ColumnMajor);
					CDGMulDimData tmp(DimTpm, c.m_type, DataMajor_FStyle_ColumnMajor);
					push(&dmo[0], m, n, y.Ptr<float>(), c.Ptr<float>(), tmp.Ptr<float>(), w.Ptr<float>());
					c = tmp;
					set_bound(1);
					VECTOR_TYPE p(8); p << vx, vx, vx, 1e-6, 1e-4, 0, 3, 2;
					fmgN(w, c, p.data(), NULL, tmp);
					memcpy(C.SelFrameByDim(3, { k1 }), tmp.ptr(), tmp.SizeInBytes());
				}
			}			
			C.Max<float>(eps);
			CDGMulDimData s;
			C.SumOnDim<float>(3, s);			
			for (int k1 = 0; k1 < Kb; k1++)
			{
				if (tc[k1].m_WarpTissue_Modulated)
				{
					auto & tissMod = sWrites.tissModulate;
					tissMod[k1].m_Img = sitk::Image(TpmSize, sitk::sitkUInt8);
					tpm.m_TpmSpace.SetToImage(tissMod[k1].m_Img);					
					CDGMulDimData tmp;
					C.SelFrameByDim(3, { k1 }, tmp);
										
					tmp.DivOn<float>(s);
					tmp.MulOn<float>(255.0);					
					tmp.ToTypeInPlace<unsigned char>(DGDataType_UnSignedChar);
					memcpy(tissMod[k1].RefImgData().ptr(), tmp.ptr(),tmp.SizeInBytes());
					tissMod[k1].SetDescription("Warped tissue class");					
				}
			}
		}
	}// end of cls
	if (TInput.m_dFieldForward)
	{			
		{
			CDGMulDimData tmp;
			MAT_TYPE mi = MAT_TYPE::Identity(4, 4);
			cout << img0mat;
			InvDef(y, dTpm, mi, img0mat, tmp);
			y = tmp; tmp.clear();			
		}
		
		ExtrapolateDef(y, M1);		
		auto & DeformationF = sWrites.FieldForwared;		
		DeformationF.m_Img = sitk::Image( DimToSize(y.m_dim), ToItkPixelType(y.m_type,0));
		tpm.m_TpmSpace.SetToImage(DeformationF.m_Img);
		DeformationF.SetDescription("Deformation");
		memcpy(DeformationF.RefImgData().ptr(), y.ptr(),y.SizeInBytes());		
	}	
}
