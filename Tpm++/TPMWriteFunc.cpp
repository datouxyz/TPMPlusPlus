#include "stdafx.h"
#include <Eigen/Dense>
#include "TissueProblityMap.h"
#include "EigenMatFunc.h"
#include "Conv2.h"
#include "shoot_boundary.h"
#include "llVar.h"
bool AnyWarpTissueSetting(sTissues & tissues, bool NativeTissue_NativeSpace, bool NativeTissue_DartelImported, bool WarpTissue_Modulated, bool WarpTissue_UnModulated);
void ZeroImageMem(sitk::Image & img);

void sWriteImgs::InitDatas(int NChan, sTissues & tc, sChannels &channels, bool dFieldInverse, bool dFieldForward, ItkImageVec &Images, sTPMPriors & tpm,
	OUT bool & do_cls, OUT bool & do_defs, int Kb, OUT CDGMulDimData & y, OUT CDGMulDimData& Q)
{
	int nTiss = tc.size();
	chans.resize(NChan);
	tissNative.resize(nTiss);
	tissImported.resize(nTiss);
	tissModulate.resize(nTiss);
	tissUnModulate.resize(nTiss);
	sitk::Image img0 = Images[0];
	DimVec img0d = SizeToDim(img0.GetSize());
	bool AnyWrapModulated = AnyWarpTissueSetting(tc, false, false, true, false);
	bool AnyWarpUnModulated = AnyWarpTissueSetting(tc, false, false, false, true);
	bool AnyTissue = AnyWarpTissueSetting(tc, true, true, true, true);

	do_cls = AnyTissue;

	for (int k1 = 0; k1 < nTiss; k1++)
	{
		sWriteImgData & tissimg = tissNative[k1];
		sTissue &tiss = tc[k1];
		if (tiss.m_WarpTissue_UnModulated || AnyWrapModulated || tiss.m_NativeTissue_DartelImported)
			do_cls = true;

		if (tiss.m_NativeTissue_NativeSpace)
		{
			tissimg.m_Img = sitk::Image(DimToSize(img0d), sitk::sitkUInt8);
			tissimg.SetDescription(FormatedString("Tissuse class %d", k1));
			do_cls = true;
		}
	}

	for (int n = 0; n < NChan; n++)
	{
		auto & chan = chans[n];
		if (channels[n].m_biascorrect)
		{
			chan.m_Nc.m_Img = img0;
			ZeroImageMem(chan.m_Nc.m_Img);
			chan.m_Nc.SetDescription("Bias corrected");
		}
		if (channels[n].m_biasfield)
		{
			chan.m_Nf.m_Img = img0;
			ZeroImageMem(chan.m_Nf.m_Img);
			chan.m_Nf.SetDescription("Estimated Bias Field");
		}
	}
	for (int k1 = 0; k1 < nTiss; k1++)
	{
		sWriteImgData & tissimg = tissNative[k1];
		sTissue &tiss = tc[k1];
		if (tiss.m_WarpTissue_UnModulated || AnyWrapModulated || tiss.m_NativeTissue_DartelImported)
			do_cls = true;

		if (tiss.m_NativeTissue_NativeSpace)
		{
			tissimg.m_Img = sitk::Image(DimToSize(img0d), sitk::sitkUInt8);
			tissimg.SetDescription(FormatedString("Tissuse class %d", k1));
			do_cls = true;
		}
	}

	do_defs = dFieldInverse || dFieldForward;
	do_defs = do_defs || do_cls;

	//sWriteImgData Ndef;
	if (do_defs)
	{
		if (dFieldInverse)
		{
			SizeVec s = { (UINT)img0d[0],(UINT)img0d[1],(UINT)img0d[2],3 };
			FieldInverse.m_Img = sitk::Image(s, sitk::sitkFloat32);
			DoubleVec dir, ori, spac;
			DirOriSpaceUpDim(img0.GetDirection(), img0.GetOrigin(), img0.GetSpacing(), dir, ori, spac);
			FieldInverse.m_Img.SetDirection(dir);
			FieldInverse.m_Img.SetSpacing(spac);
			FieldInverse.m_Img.SetOrigin(ori);

			ZeroImageMem(FieldInverse.m_Img);
			FieldInverse.SetDescription("Inverse Deformation");
		}
		if (dFieldForward || AnyWarpTissueSetting(tc, false, true, true, true))
		{
			y = CDGMulDimData({ img0d[0],img0d[1],img0d[2],3 }, DGDataType_Float,DataMajor_FStyle_ColumnMajor);
		}
	}
	if (do_cls)
	{
		Q = CDGMulDimData({ img0d[0],img0d[1],img0d[2],Kb }, DGDataType_Float, DataMajor_FStyle_ColumnMajor);
	}
}

