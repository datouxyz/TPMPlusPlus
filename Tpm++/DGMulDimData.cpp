#include "stdafx.h"

#include "DGMulDimData.h"

CDGMulDimData::CDGMulDimData()
	:m_type(DGDataType_NotDefined), m_Acc(DimVec(),DataMajor_FStyle_ColumnMajor),m_RefPtr(NULL)
{

}

CDGMulDimData::CDGMulDimData(DimVec dim, DGDataType type,DataMajor dm, void* dataptr, bool ref_datasptr)
	:m_dim(dim), m_type(type),m_Acc(dim,dm),m_RefPtr(NULL)
{
	m_ElSize = DGDataSize(type);
	m_ElCount = DimElCount(dim);

	if (ref_datasptr)
	{
		m_RefPtr = dataptr;
	}
	else
	{

		m_bytes.resize(m_ElSize*m_ElCount);
		
		if (dataptr)
			memcpy(ptr(), dataptr, SizeInBytes());
	}	
}

void CDGMulDimData::Print()
{
	if (m_ElCount <= 0) cout << "empty";

	m_Acc.TraverseNoRecusive((TraverseFunc)
		[&](DimIndex & index)->bool
	{		
		cout << "[ ";
		for (auto &id : index)
		{
			cout << id << " ";
		}
		cout << "] ";
		
		return true;
	}
	);
}
void CDGMulDimData::SelFrameByDim(int nSelDim, DimIndex FrameIndex, CDGMulDimData & out)
{
	DimVec dd; dd.insert(dd.end(),m_dim.begin(), m_dim.begin()+nSelDim);
	assert(m_Acc.m_DataMajor == DataMajor_FStyle_ColumnMajor);

	out = CDGMulDimData(dd,m_type,m_Acc.m_DataMajor ,SelFrameByDim(nSelDim, FrameIndex),m_RefPtr!=NULL );
}

void* CDGMulDimData::SelFrameByDim(int nFrameDim,DimIndex FrameIndex)
{
	if (m_dim.size() < nFrameDim)
		return NULL;

	assert(m_Acc.m_DataMajor == DataMajor_FStyle_ColumnMajor);

	int framesize = 1;
	for (int ifd = 0; ifd < nFrameDim; ifd++)
	{
		framesize *= m_dim[ifd];
	}		

	DimVec pdim;
	for (int i = nFrameDim; i < m_dim.size(); i++)
	{
		pdim.push_back(m_dim[i]);
	}
	CDimDataAccer acc(pdim,m_Acc.m_DataMajor);
	int pos = acc.GetElPos(FrameIndex)* framesize;// , DataMajor_FStyle_ColumnMajor);
	// GetElPosHighDimPriority(FrameIndex) * framesize;

	return &BytePtr()[pos * m_ElSize];
}

void* CDGMulDimData::SelMatByDimXY(DimIndex FrameIndex)
{
	if (m_dim.size() < 2)
		return NULL;

	int framesize = m_dim[0] * m_dim[1];
	assert(m_Acc.m_DataMajor == DataMajor_FStyle_ColumnMajor);

	DimVec pdim;
	for (int i = 2; i < m_dim.size(); i++)
	{
		pdim.push_back(m_dim[i]);
	}
	CDimDataAccer acc(pdim,m_Acc.m_DataMajor);
	int pos = acc.GetElPos(FrameIndex)* framesize;
	//, DataMajor_FStyle_ColumnMajor)* framesize;
	//GetElPosHighDimPriority(FrameIndex) 
	return &BytePtr()[pos * m_ElSize];
}



template<typename T> void BsxfTimes(byte* A, byte *B,byte* ret)
{
	T res=  (*(T*)A) * (*(T*)B);
	memcpy(ret, &res, sizeof(res));
}

template<typename T> void BsxfMinus(byte* A, byte *B, byte* ret)
{
	T res = (*(T*)A) - (*(T*)B);
	memcpy(ret, &res, sizeof(res));
}
template<typename T> void BsxfPlus(byte* A, byte *B, byte* ret)
{
	T res = (*(T*)A) + (*(T*)B);
	memcpy(ret, &res, sizeof(res));
}
template<typename T> void BsxfDivide(byte* A, byte *B, byte* ret)
{
	T res = (*(T*)A) / (*(T*)B);
	memcpy(ret, &res, sizeof(res));
}
template<typename T> void BsxfRoundedDivide(byte* A, byte *B, byte* ret)
{
	double da = (*(T*)A),db = (*(T*)B);

	*(T*)ret = (T)round(da / db);

}

#define BSXFFUNC_DEF(typen,type,func) \
void Bsxf_##typen##_##func##(byte* A, byte *B, byte* ret){\
	Bsxf##func<type>(A, B, ret);}


#define BSXFFUNC_DEFS(OP_) \
BSXFFUNC_DEF(Float, float, OP_);\
BSXFFUNC_DEF(Double, double, OP_);\
BSXFFUNC_DEF(Int, int, OP_);\
BSXFFUNC_DEF(UInt, unsigned int, OP_);\
BSXFFUNC_DEF(Char, char, OP_);\
BSXFFUNC_DEF(UChar,unsigned char, OP_);\
BSXFFUNC_DEF(Short, short, OP_);\
BSXFFUNC_DEF(UShort, unsigned short, OP_);\
BSXFFUNC_DEF(Int64, __int64, OP_);\
BSXFFUNC_DEF(UInt64,unsigned __int64, OP_);

BSXFFUNC_DEFS(Minus);
BSXFFUNC_DEFS(Times);
BSXFFUNC_DEFS(Plus);
BSXFFUNC_DEFS(Divide);
BSXFFUNC_DEFS(RoundedDivide);


#define SWITCHNCALL_BSXF(type,Op) \
	switch (type)\
	{\
		case DGDataType_Float:\
			return Bsxf_Float_##Op;\
			break;\
		case DGDataType_Double:\
			return Bsxf_Double_##Op;\
			break;\
		case DGDataType_Int:\
			return Bsxf_Int_##Op;\
			break;\
		case DGDataType_UnSignedInt:\
			return Bsxf_UInt_##Op;\
			break;\
		case DGDataType_Char:\
			return Bsxf_Char_##Op;\
			break;\
		case DGDataType_UnSignedChar:\
			return Bsxf_UChar_##Op;\
			break;\
		case DGDataType_Int64:\
			return Bsxf_Int64_##Op;\
			break;\
		case DGDataType_UnsignedInt64:\
			return Bsxf_UInt64_##Op;\
			break;\
		case DGDataType_SignedShort:\
			return Bsxf_Short_##Op;\
			break;\
		case DGDataType_UnSignedShort:\
			return Bsxf_UShort_##Op;\
			break;\
	}


BsxfFuncPointer GetBxsfFunc(BsxOperator op, DGDataType type)
{
	switch (op)
	{
	case BsxOperator_times:
	{
		SWITCHNCALL_BSXF(type, Times);		
	}
	break;
	case BsxOperator_plus:
		SWITCHNCALL_BSXF(type, Plus);
		break;
	case BsxOperator_minus:
		SWITCHNCALL_BSXF(type, Minus);
		break;
	case  BsxOperator_divide:
		SWITCHNCALL_BSXF(type, Divide);
		break;
	case  BsxOperator_rounded_divide:
		SWITCHNCALL_BSXF(type, RoundedDivide);
		break;
	}
	
	return NULL;
}

byte* CDGMulDimData::DataPtr(DimIndex ind)
{
	int pos = m_Acc.GetElPos(ind);		
	return (byte*)&BytePtr()[m_ElSize*pos];
}

byte* CDGMulDimData::DataPtr(int pos)
{
	return (byte*)&BytePtr()[m_ElSize*pos];
}

double CDGMulDimData::SumAsDouble()
{
	if (m_type != DGDataType_Double) return 0;

	double dsum=0;
	for (int i = 0; i < m_ElCount; i++)
	{
		dsum += *(double*)&BytePtr()[m_ElSize*i];
	}
	return dsum;
}
float CDGMulDimData::SumAsFloat()
{
	if (m_type != DGDataType_Float) return 0;

	double fsum = 0;
	for (int i = 0; i < m_ElCount; i++)
	{
		fsum += *(float*)&BytePtr()[m_ElSize*i];
	}
	return fsum;
}

double CDGMulDimData::Sum()
{
	double dsum = 0;
	
	for (int i = 0; i <m_ElCount;i++)
	{
		double v = 0;
		DGDataCastToT(DataPtr(i), m_type, v);
		dsum += v;
	}
	return dsum;
}

CDGMulDimData CDGMulDimData::ReShape(DimVec newdim)
{
	if (!m_RefPtr)
	{
		CDGMulDimData newd(newdim,m_type,m_Acc.m_DataMajor,ptr());
		//newd.m_Bytes = m_bytes;
		return newd;
	}
	else
	{
		CDGMulDimData newd(newdim, m_type, m_Acc.m_DataMajor,m_RefPtr,true);
		return newd;
	}

}
bool CDGMulDimData::ReShapeInPlace(DimVec newdim)
{
	if (DimElCount(newdim) != DimElCount(m_dim))
	{
		ErrOutput("ReShapeInPlace" , "el count not match");
		return false;
	}
	m_dim = newdim;
	return true;
}
#include "CasterVec.h"
void CDGMulDimData::ToSingle(CDGMulDimData & to)
{
	to = CDGMulDimData(m_dim,DGDataType_Float, m_Acc.m_DataMajor);
	float *fp = (float*)&to.BytePtr()[0];

	DGDataCastToT_VEC<float>(m_ElCount,&BytePtr()[0],m_type, fp);
}
 
void CDGMulDimData::ToDouble(CDGMulDimData & to)
{
	to = CDGMulDimData(m_dim, DGDataType_Double, m_Acc.m_DataMajor);
	double *dp = (double*)&to.BytePtr()[0];
	DGDataCastToT_VEC<double>(m_ElCount,&BytePtr()[0], m_type, dp);
}


#include "ParamaterDef.h"
bool CDGMulDimData::BsxFuncHighDimDup(BsxOperator op, CDGMulDimData & rhs,int nHighDupDim ,CDGMulDimData& outres)
{

	if (rhs.m_dim.size() != m_dim.size())
	{
		ErrOutput("BsxFuncHighDimDup", u8"����������ά����ͬ�޷�ִ��");
		return false;
	}
	if (rhs.m_type != m_type)
	{
		ErrOutput("BsxFuncHighDimDup", u8"�������������Ͳ�ͬ");
		return false;
	}
	assert(rhs.m_Acc.m_DataMajor == DataMajor_FStyle_ColumnMajor);
	if (rhs.m_Acc.m_DataMajor != m_Acc.m_DataMajor)
	{
		ErrOutput("BsxFuncHighDimDup", u8"����������DataMajor��ͬ");
		return false;
	}
	BsxfFuncPointer bsxffunc = GetBxsfFunc(op, m_type);
	int nLowDim = m_dim.size() - nHighDupDim;
	DimVec HiDim, LowDim;// (nHighDupDim), LowDim(m_dim.size() - nHighDupDim);
	for (int idim = 0; idim < m_dim.size(); idim++)
	{
		if (idim < nLowDim)
		{
			LowDim.push_back(m_dim[idim]);
			if (rhs.m_dim[idim] != 1)
			{
				ErrOutput("BsxFuncHighDimDup", u8"rhsHiDim������ַ�1ά��");
				return false;
			}
		}			
		else
		{
			HiDim.push_back(m_dim[idim]);			
		}
	}

	int nFrameSize = DimElCount(LowDim);
	int nFrame = DimElCount(m_dim) / nFrameSize;

	if (nFrame != rhs.m_ElCount)
	{
		ErrOutput("BsxFuncHighDimDup", u8"lhs rhs��ƥ��");
		return false;
	}

	outres = CDGMulDimData(m_dim, m_type, m_Acc.m_DataMajor);
	for (int iFrame = 0; iFrame < nFrame; iFrame++)
	{
		byte * resptr = &outres.BytePtr()[iFrame * nFrameSize * m_ElSize];
		byte* rhsptr = &rhs.BytePtr()[0];
		byte * lhsframeptr = &BytePtr()[iFrame * nFrameSize * m_ElSize];
		for (int iel = 0; iel < nFrameSize; iel++)
		{
			bsxffunc( &lhsframeptr[iel * m_ElSize], &rhsptr[iFrame * m_ElSize], &resptr[iel*m_ElSize]);
		}
	}	
	return true;
}

bool CDGMulDimData::BsxFuncHighDim_InPlace(BsxOperator op, CDGMulDimData & rhs, int nHighDupDim)
{

	if (rhs.m_dim.size() != m_dim.size())
	{
		ErrOutput("BsxFuncHighDimDup", u8"size not same");
		return false;
	}
	if (rhs.m_type != m_type)
	{
		ErrOutput("BsxFuncHighDimDup", u8"type not same");
		return false;
	}
	BsxfFuncPointer bsxffunc = GetBxsfFunc(op, m_type);
	int nLowDim = m_dim.size() - nHighDupDim;
	DimVec HiDim, LowDim;// (nHighDupDim), LowDim(m_dim.size() - nHighDupDim);
	for (int idim = 0; idim < m_dim.size(); idim++)
	{
		if (idim < nLowDim)
		{
			LowDim.push_back(m_dim[idim]);
			if (rhs.m_dim[idim] != 1)
			{
				ErrOutput("BsxFuncHighDimDup", u8"rhsHiDim");
				return false;
			}
		}
		else
		{
			HiDim.push_back(m_dim[idim]);

		}
	}

	int nFrameSize = DimElCount(LowDim);
	int nFrame = DimElCount(m_dim) / nFrameSize;

	if (nFrame != rhs.m_ElCount)
	{
		ErrOutput("BsxFuncHighDimDup", u8"lhs rhs");
		return false;
	}

	for (int iFrame = 0; iFrame < nFrame; iFrame++)
	{
		byte* rhsptr = &rhs.BytePtr()[0];
		byte * lhsframeptr = &BytePtr()[iFrame * nFrameSize * m_ElSize];
		for (int iel = 0; iel < nFrameSize; iel++)
		{
			bsxffunc(&lhsframeptr[iel * m_ElSize], &rhsptr[iFrame * m_ElSize], &lhsframeptr[iel*m_ElSize]);
		}
	}
	return true;
}
bool CDGMulDimData::BsxFuncLowDimDup(BsxOperator op, CDGMulDimData & rhs, CDGMulDimData& outres)
{

	if (rhs.m_dim.size() != m_dim.size())
	{
		ErrOutput("BsxFuncLowDimDup", u8"");
		return false;
	}

	if (rhs.m_type != m_type)
	{
		ErrOutput("BsxFuncLowDimDup", u8"");
		return false;
	}
	outres = CDGMulDimData(m_dim, m_type, m_Acc.m_DataMajor);
	
	int batchsize = rhs.SizeInBytes();
	int elsize = DGDataSize(m_type);
	int batchelcount = batchsize/ elsize;
	int nbatch = SizeInBytes() / batchsize;
	BsxfFuncPointer bsxffunc = GetBxsfFunc(op, m_type);

	for (int ibatch = 0; ibatch < nbatch; ibatch++)
	{
		byte* APtr = &BytePtr()[ibatch * batchsize];
		byte* ResPtr = &outres.BytePtr()[ibatch * batchsize];
		byte* BPtr = &rhs.BytePtr()[0];

		for (int ielbytes = 0; ielbytes < batchsize; ielbytes+=elsize)
		{
			//int off = iel * elsize;
			bsxffunc( &APtr[ielbytes], &BPtr[ielbytes], &ResPtr[ielbytes]);
		}
	}	
	return true;
}

bool CDGMulDimData::BsxFuncLowDim_InPlace(BsxOperator op, CDGMulDimData & rhs)
{

	if (rhs.m_dim.size() != m_dim.size())
	{
		ErrOutput("BsxFuncLowDimDup", u8"");
		return false;
	}

	if (rhs.m_type != m_type)
	{
		ErrOutput("BsxFuncLowDimDup", u8"");
		return false;
	}
	//outres = CDGMulDimData(m_dim, m_type);

	int batchsize = rhs.SizeInBytes();
	int elsize = DGDataSize(m_type);
	int batchelcount = batchsize / elsize;
	int nbatch = SizeInBytes() / batchsize;
	BsxfFuncPointer bsxffunc = GetBxsfFunc(op, m_type);

	for (int ibatch = 0; ibatch < nbatch; ibatch++)
	{
		byte* APtr = &BytePtr()[ibatch * batchsize];
		byte* ResPtr = &BytePtr()[ibatch * batchsize];
		byte* BPtr = &rhs.BytePtr()[0];

		for (int ielbytes = 0; ielbytes < batchsize; ielbytes += elsize)
		{
			//int off = iel * elsize;
			bsxffunc(&APtr[ielbytes], &BPtr[ielbytes], &ResPtr[ielbytes]);
		}
	}
	return true;
}

 bool CDGMulDimData::BsxFunction(BsxOperator op, CDGMulDimData & rhs, CDGMulDimData& outres)
{

	 if (rhs.m_dim.size() != m_dim.size())
	 {
		 ErrOutput("BsxFunction", u8"Dim not match");
		 return false;
	 }
	 
	 if (rhs.m_type != m_type)
	 {
		 ErrOutput("BsxFunction", u8"type not match");
		 return false;
	 }
	 DimVec resultdim=m_dim;
	 for (int idim = 0; idim < m_dim.size(); idim++)
	 {
		 // singleton expansion test
		 if (m_dim[idim] != rhs.m_dim[idim] && m_dim[idim]!=1 && rhs.m_dim[idim]!=1)
		 {
			 ErrOutput("BsxFunction",u8"singleton expansion test fail");
			 return false;
		 }
		 resultdim[idim] = max(m_dim[idim],rhs.m_dim[idim]);
	 }

	 BsxfFuncPointer bsxffunc = GetBxsfFunc(op, m_type);
	 
	 outres = CDGMulDimData(resultdim, m_type, m_Acc.m_DataMajor);
	 int elsize = DGDataSize(m_type);

	 outres.m_Acc.TraverseNoRecusive((TraverseFunc)
		 [&](DimIndex & resindex)->bool
	 {
		 DimIndex Aindex = resindex;
		 DimIndex Bindex = resindex;

		 for (int idim = 0; idim < resindex.size(); idim++)
		 {
			 Aindex[idim] = resindex[idim]%m_dim[idim];
			 Bindex[idim] = resindex[idim]%rhs.m_dim[idim];
		 }
		
		 byte * A = DataPtr(Aindex);// &m_Bytes[lhspos * elsize];
		 byte * B = rhs.DataPtr(Bindex); //&rhs.m_Bytes[rhspos * elsize];
		 byte* r = outres.DataPtr(resindex);//&result.m_Bytes[resultpos * elsize];
 		 bsxffunc(A,B,r);
		 return true;
	 }
	 );
	 
	 //outres = result;
	 return true;		  
}



 int TestDim()
 {
	 DimVec Da = {10,2,1,6};
	 
	 CDGMulDimData A(Da,DGDataType_Double, DataMajor_FStyle_ColumnMajor);

	 for (int i = 0; i < A.m_ElCount; i++)
	 {
		 double * dp = (double*)A.DataPtr(i);
		 *dp = i;
	 }

	 DimVec Db = { 1,1,1,6 };
	 CDGMulDimData B(Db,DGDataType_Double, DataMajor_FStyle_ColumnMajor);
	 
	 for (int i = 0; i < B.m_ElCount; i++)
	 {
		 double * dp = (double*)B.DataPtr(i);
		 *dp = i * 0.1;
	 }

	 CDGMulDimData result;
	 bool bsuccess = A.BsxFunction(BsxOperator_times , B,result);
	
	 return 1;
 }
 

 int Slice1(double * mat, float* image, int xdim1, int ydim1, CDGMulDimData& vol, double background, DoubleVec & scale, DoubleVec &offset)
 {
#define TINY 5e-2

	 double y, x2, y2, z2, s2, dx3 = mat[0], dy3 = mat[1], dz3 = mat[2], ds3 = mat[3];
	 int t = 0;

	 int zdim2 = vol.m_dim[2];
	 int ydim2 = vol.m_dim[1];
	 int xdim2 = vol.m_dim[0];
	 x2 = mat[12] + 0 * mat[8];
	 y2 = mat[13] + 0 * mat[9];
	 z2 = mat[14] + 0 * mat[10];
	 s2 = mat[15] + 0 * mat[11];

	 for (y = 1; y <= ydim1; y++)
	 {
		 double x;
		 double x3 = x2 + y * mat[4];
		 double y3 = y2 + y * mat[5];
		 double z3 = z2 + y * mat[6];
		 double s3 = s2 + y * mat[7];
		 for (x = 1; x <= xdim1; x++)
		 {
			 double x4, y4, z4;
			 s3 += ds3;
			 if (s3 == 0.0) return(-1);
			 x4 = (x3 += dx3) / s3 - 1.0;
			 y4 = (y3 += dy3) / s3 - 1.0;
			 z4 = (z3 += dz3) / s3 - 1.0;

			 if (z4 >= -TINY && z4 < zdim2 + TINY - 1 &&
				 y4 >= -TINY && y4 < ydim2 + TINY - 1 &&
				 x4 >= -TINY && x4 < xdim2 + TINY - 1)
			 {
				 double k111, k112, k121, k122, k211, k212, k221, k222;
				 double dx1, dx2, dy1, dy2, dz1, dz2;
				 int off1, off2, offx, offy, offz, ix4, iy4, iz4;

				 ix4 = floor(x4); dx1 = x4 - ix4; dx2 = 1.0 - dx1;
				 iy4 = floor(y4); dy1 = y4 - iy4; dy2 = 1.0 - dy1;
				 iz4 = floor(z4); dz1 = z4 - iz4; dz2 = 1.0 - dz1;

				 ix4 = (ix4 < 0) ? ((offx = 0), 0) : ((ix4 >= xdim2 - 1) ? ((offx = 0), xdim2 - 1) : ((offx = 1), ix4));
				 iy4 = (iy4 < 0) ? ((offy = 0), 0) : ((iy4 >= ydim2 - 1) ? ((offy = 0), ydim2 - 1) : ((offy = xdim2), iy4));
				 iz4 = (iz4 < 0) ? ((offz = 0), 0) : ((iz4 >= zdim2 - 1) ? ((offz = 0), zdim2 - 1) : ((offz = 1), iz4));

				 off1 = ix4 + xdim2 * iy4;
				 off2 = off1 + offy;

				 float * piZ4 = (float*)vol.SelMatByDimXY({ iz4 });

				 k222 = piZ4[off1]; k122 = piZ4[off1 + offx];
				 k212 = piZ4[off2]; k112 = piZ4[off2 + offx];

				 float * piZ4offz = (float*)vol.SelMatByDimXY({ iz4 + offz });
				 k221 = piZ4offz[off1]; k121 = piZ4offz[off1 + offx];
				 k211 = piZ4offz[off2]; k111 = piZ4offz[off2 + offx];

				 /* resampled pixel value (trilinear interpolation) */
				 image[t] = (((k222*dx2 + k122 * dx1)*dy2 + (k212*dx2 + k112 * dx1)*dy1)*scale[iz4] + offset[iz4])*dz2
					 + (((k221*dx2 + k121 * dx1)*dy2 + (k211*dx2 + k111 * dx1)*dy1)*scale[iz4 + offz] + offset[iz4 + offz])*dz1;

			 }
			 else image[t] = background;
			 t++;
		 }
	 }
	 return(0);
 }