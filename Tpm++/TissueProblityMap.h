#pragma once

#include "ParamaterDef.h"
#include <typeinfo>

#include "DGMulDimData.h"

#include "msk.h"
#include "TInputs.h"

bool LoadBytesFromFile(ByteVec & readto, string filepath);
bool WriteBytesToNewFile(ByteVec & tosave, string filepath);
const double eps = 2.2204e-16;
inline VECTOR_TYPE vMax(VECTOR_TYPE & v1, VECTOR_TYPE & v2)
{
	VECTOR_TYPE ret(v1.size());
	if (v1.size() != v2.size())
		return ret;

	
	for (int i = 0;i < v1.size();i++)
	{
		ret[i] = max(v1[i], v2[i]);
	}
	return ret;
}

inline void Rounds(VECTOR_TYPE & rounds)
{
	for (auto & r : rounds)
	{
		r = round(r);
	}
}
inline void Rounds(VECTOR_TYPE_F & rounds)
{
	for (auto & r : rounds)
	{
		r = round(r);
	}
}
inline void Clamp(VECTOR_TYPE_F & rounds,float fmin,float fmax)
{
	for (auto & r : rounds)
	{
		r = std::clamp(r,fmin,fmax);
		
	}
}
inline void Clamp(VECTOR_TYPE & rounds, double dmin, double dmax)
{
	for (auto & r : rounds)
	{
		r = std::clamp(r, dmin, dmax);
	}
}

template<typename TInd, typename Tdata>
Tdata Accumarray(TInd &ind, Tdata &data, int outsize)
{
	Tdata ret = Tdata::Zero(outsize);

	if (ind.size() != data.size())
		return ret;
	outsize = clamp<int>(outsize, 0, data.size());
	for (int i = 0; i < data.size(); i++)
	{
		ret[ind[i] - 1] += data[i];
	}
	return ret;
}

VECTOR_TYPE Accumarray(VECTOR_TYPE_F &ind, VECTOR_TYPE_F &data, int outsize);

VECTOR_TYPE Accumarray(VECTOR_TYPE &ind,VECTOR_TYPE &data, int outsize);

VECTOR_TYPE Accumarray(VECTOR_TYPE &ind, double singledata, int outsize);

inline VECTOR_TYPE ToVecType(DoubleVec & dv)
{
	VECTOR_TYPE ret(dv.size());
	for (int i = 0;i < dv.size();i++)
	{
		ret[i] = dv[i];
	}
	return ret;
}
inline DoubleVec FromVecType(VECTOR_TYPE & dv)
{
	DoubleVec ret(dv.size());
	for (int i = 0; i < dv.size(); i++)
	{
		ret[i] = dv[i];
	}
	return ret;
}
struct sWriteImgData
{
	sitk::Image m_Img;
	CDGMulDimData RefImgData()
	{
		int o;
		return CDGMulDimData(SizeToDim(m_Img.GetSize()), ToDGDataType(m_Img.GetPixelID(), o),
			DataMajor_FStyle_ColumnMajor,m_Img.GetBufferAsVoid(),true);
	}
	void SetDescription(string des)
	{
		string key = "description";
		m_Img.SetMetaData(key,des);
	}

	bool IsEmpty()
	{
		return m_Img.GetNumberOfPixels() == 0;
		//m_dat.m_Bytes.size() ==0;
	}
};

struct sIntensity
{
	MAT_TYPE_F m_lik;
	VECTOR_TYPE_F m_interscal;
};
struct sChanWrite
{
	MAT_TYPE B1;
	MAT_TYPE B2;
	MAT_TYPE B3;
	CDGMulDimData T;
	// BiasCorrect
	sWriteImgData m_Nc;
	//BiasField_
	sWriteImgData m_Nf;
};


struct sTpmVol3D
{
public:
	sTpmVol3D(sitk::Image & tpmimg4d);
	sTpmVol3D() {};
	DoubleVec m_Ori;
	DoubleVec m_Dir;
	DoubleVec m_Spacing;

	DoubleVec m_Ori4D;
	DoubleVec m_Dir4D;
	DoubleVec m_Spacing4D;

	MAT_TYPE GetMat(bool bApplyMatlabIndexOffset, bool bItkSpaceToNifity);
	void SetToImage(sitk::Image & img);		
};

void SplitDirOriSpacing(OUT DoubleVec & Dirs, DoubleVec & Oris, DoubleVec & Spacings,
	IN DoubleVec  ODirs, IN DoubleVec  OOris, IN DoubleVec  OSpacings, int iSplitDim);

void SplitDirOriSpacingOn3DVol(OUT DoubleVec & Dirs, DoubleVec & Oris, DoubleVec & Spacings,
	IN DoubleVec  ODirs, IN DoubleVec  OOris, IN DoubleVec  OSpacings);


void DirOriSpaceUpDim(IN DoubleVec  ODirs, IN DoubleVec  OOris, IN DoubleVec  OSpacings,
	OUT DoubleVec & Dirs, DoubleVec & Oris, DoubleVec & Spacings);

struct sTPMPriors
{
public:
	sTPMPriors();	
	double m_deg{ 1 };
	MAT_TYPE m_TPMMat;	
	VECTOR_TYPES m_Datas;
	SizeVec m_size;
	VECTOR_TYPE m_bg1;
	VECTOR_TYPE m_bg2;	
	int m_nTissuse{ 0 };
	void SamplePriors(VECTOR_TYPE &_x1, VECTOR_TYPE &_x2, VECTOR_TYPE &_x3, OUT VECTOR_TYPES & s);
	void SamplePriorsDerive(VECTOR_TYPE &_x1, VECTOR_TYPE &_x2, VECTOR_TYPE &_x3,
		OUT VECTOR_TYPES &s, OUT VECTOR_TYPES &ds1, OUT VECTOR_TYPES &ds2, OUT VECTOR_TYPES &ds3);
	void LoadPriors(sitk::Image & img);
	sTpmVol3D m_TpmSpace;
};

struct sChan
{
	MAT_SPARSE_TYPE C;
	MAT_TYPE B1;
	MAT_TYPE B2;
	MAT_TYPE B3;
	CDGMulDimData T;
	double ll{ 0 };
	VECTOR_TYPE_F interscal;
	MAT_TYPE hist;

	VECTOR_TYPE lam;
	MAT_TYPE_F lik;
	MAT_TYPE_F alph;

	MAT_TYPE_F grad1;
	MAT_TYPE_F grad2;
};

typedef std::vector<sChan> ChanVec;


struct sWrapBUF
{
	sMask msk;	
	vector<VECTOR_TYPE_F> f;
	//TPM datainfo
	MAT_TYPE_F dat;
	//wrap
	vector<VECTOR_TYPE_F> bf;
};

#include "DGMulDimData.h"

struct sRegBuffer
{	
	sMask m_msk;
	sMask m_mskz;
	VECTOR_TYPE m_g;
	VECTOR_TYPES m_b;
};
typedef std::vector<sRegBuffer> sRegBuffers;


struct sRegData
{
public:
	sRegData(sitk::Image & movimg, int nsample, double fwhm);
	//voxelsize
	VECTOR_TYPE m_vx;
	VECTOR_TYPE m_sk;
	MAT_TYPE m_MG;

	double m_FudgeFactor{ 0 };

	// init x1,x2,x3 m_vx,m_sk
	void InitCoords(SizeVec movsize);
	MAT_TYPE x1;
	MAT_TYPE x2;
	VECTOR_TYPE x3;

	sRegBuffers m_Buffers;
	CDGMulDimData TBias;	
};


struct sEstimateParam
{
	int K{ 0 };
	int Kb{ 0 };
	double eps{ 2.2204e-16 };
	double tiny{ 2.2204e-16 *2.2204e-16 };
	double llrb{ 0 };
	float llr{ 0 };
	double ll_const{ 0 };
	double tol1{ 1e-4 };
	int nm{ 0 };
};
enum KernelSmoothMode
{
	KernelSmoothMode_nearest = 0,
	KernelSmoothMode_linear,
	KernelSmoothMode_MAX = 0xffffffff
};

#include "llVar.h"
struct sWrapResults
{
	MAT_TYPE m_Affine;
	MAT_TYPE m_MT;
	sTPMPriors m_tpm;
	VECTOR_TYPE m_lkp;
	CDGMulDimData m_Twrap;
	vector<CDGMulDimData> m_T;
	VECTOR_TYPE m_WrapPrior;
	MAT_TYPE  m_mn;
	VECTOR_TYPE m_mg;
	vector<MAT_TYPE> m_vr;
	vector<sIntensity> m_intensity;
	CVariableLL m_LL;
};


sitk::Image ReSampleImgAsFloat(sitk::Image & img, int nsample);

typedef vector<sChanWrite> ChanWrites;
struct sWriteImgs
{
	void InitDatas(int NChan, sTissues & tc, sChannels &channels, bool dFieldInverse, bool dFieldForward,
		ItkImageVec &Images, sTPMPriors & tpm, OUT bool & do_cls, OUT bool & do_defs, int Kb, OUT CDGMulDimData & y, OUT CDGMulDimData& Q);
	ChanWrites chans;
	//NativeTissue Tissuse class NativeSpace
	vector< sWriteImgData> tissNative;// (Kb);
	//NativeTissue Tissuse class Imported Space
	vector< sWriteImgData> tissImported;// (Kb);
	//Wraped Tissue Modulated
	vector< sWriteImgData> tissModulate;// (Kb);
	//Wraped Tissue UnModulated
	vector< sWriteImgData> tissUnModulate;// (Kb);

	sWriteImgData FieldInverse;
	sWriteImgData FieldForwared;
};

double randNumber1000[];
