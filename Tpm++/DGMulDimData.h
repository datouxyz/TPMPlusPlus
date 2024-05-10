#pragma once

typedef void(*BsxfFuncPointer)(byte* A, byte *B, byte* ret);
#include "ParamaterDef.h"
#include "DimDataAccer.h"


inline sitk::PixelIDValueEnum ToItkPixelType(DGDataType type,int ndim)
{
	switch (type)
	{
		case DGDataType_Float:		
			return ndim>1? sitk::sitkVectorFloat32 : sitk::sitkFloat32;
			break;
		case DGDataType_Double:		
			return ndim > 1 ? sitk::sitkVectorFloat64: sitk::sitkFloat64;
			break;
		case DGDataType_Int:
			return ndim > 1 ? sitk::sitkVectorInt32 : sitk::sitkInt32;
			break;		
		case DGDataType_UnSignedChar:
			return ndim > 1 ? sitk::sitkVectorUInt8 : sitk::sitkUInt8;
			break;
		case DGDataType_UnSignedShort:
			return ndim > 1 ? sitk::sitkVectorUInt16 : sitk::sitkUInt16;
			break;
		case DGDataType_SignedShort:
			return ndim > 1 ? sitk::sitkVectorInt16 : sitk::sitkInt16;
			break;
		case DGDataType_UnSignedInt:
			return ndim > 1 ? sitk::sitkVectorUInt32 : sitk::sitkUInt32;
			break;
		case DGDataType_Int64:
			return ndim > 1 ? sitk::sitkVectorInt64 : sitk::sitkInt64;
			break;
		case DGDataType_UnsignedInt64:
			return ndim > 1 ? sitk::sitkVectorUInt64 : sitk::sitkUInt64;
			break;			
	}
	return sitk::sitkUnknown;
}

inline bool IsSingleCompeontType(sitk::PixelIDValueEnum itkenum)
{
	return itkenum == sitk::sitkUInt8 || itkenum == sitk::sitkInt8 ||
		itkenum == sitk::sitkUInt16 || itkenum ==sitk::sitkInt16 ||
		itkenum == sitk::sitkUInt32 || itkenum == sitk::sitkInt32 ||
		itkenum == sitk::sitkUInt32 || itkenum == sitk::sitkInt32 ||
		itkenum == sitk::sitkUInt64 || itkenum == sitk::sitkInt64 ||
		itkenum == sitk::sitkFloat32 || itkenum == sitk::sitkFloat64 ||
		itkenum == sitk::sitkUInt16 || itkenum == sitk::sitkInt16 ||
		itkenum == sitk::sitkLabelUInt8 || itkenum == sitk::sitkLabelUInt16 ||
		itkenum == sitk::sitkLabelUInt32 || itkenum == sitk::sitkLabelUInt64;

	//return itkenum >= sitkUInt8 && itkenum <= sitkFloat64;
}

inline bool IsVectorCompeontType(sitk::PixelIDValueEnum itkenum)
{
	return !IsSingleCompeontType(itkenum);
}

inline DGDataType ToDGDataType(sitk::PixelIDValueEnum itkenum,int &nOutDim)
{
	if (IsSingleCompeontType(itkenum))
	{
		nOutDim = 1;
	}
	else if(IsVectorCompeontType(itkenum))
	{
		nOutDim = 2;
	}
	else
		nOutDim = -1;

	switch (itkenum)
	{
	case sitk::sitkFloat32:
	case sitk::sitkVectorFloat32:
		return DGDataType_Float;
		break;
	case sitk::sitkFloat64:
	case sitk::sitkVectorFloat64:
		return DGDataType_Double;
		break;


	case sitk::sitkInt64:
	case sitk::sitkVectorInt64:
		return DGDataType_Int64;
		break;
	case sitk::sitkInt32:
	case sitk::sitkVectorInt32:
		return DGDataType_Int;
		break;	
	case sitk::sitkInt16:
	case sitk::sitkVectorInt16:
		return DGDataType_SignedShort;
		break;
	case sitk::sitkInt8:
	case sitk::sitkVectorInt8:
		return DGDataType_Char;
		break;



	case sitk::sitkUInt8:
	case sitk::sitkVectorUInt8:
		return DGDataType_UnSignedChar;
		break;
	case sitk::sitkUInt16:
	case sitk::sitkVectorUInt16:
		return DGDataType_UnSignedShort;
		break;

	case sitk::sitkUInt32:
	case sitk::sitkVectorUInt32:
		return DGDataType_UnSignedInt;
		break;
	
	case sitk::sitkUInt64:
	case sitk::sitkVectorUInt64:
		return DGDataType_UnsignedInt64;
		break;
	}
	return DGDataType_NotDefined;
}
class CDGMulDimData
{
public:
	CDGMulDimData();

	CDGMulDimData(DimVec dim, DGDataType type, DataMajor dm,void* dataptr=NULL, bool ref_datasptr=false);
	
	 
	void * m_RefPtr;

	CDGMulDimData ReShape(DimVec newdim);


	void clear()
	{
		*this = CDGMulDimData();
	};
	template<typename T>
	void ToType(CDGMulDimData & to, DGDataType totype)
	{
		to = CDGMulDimData(m_dim, totype,m_Acc.m_DataMajor);
		T *dp = (T*)&to.BytePtr()[0];
		DGDataCastToT_VEC<T>(m_ElCount,(T*) &BytePtr()[0], m_type, dp);
	}

	void ToSingle(CDGMulDimData & to);
	void ToDouble(CDGMulDimData & to);
	template<typename T>
	void ToTypeInPlace(DGDataType totype)
	{
		CDGMulDimData newt;
		ToType<T>(newt,totype);
		*this = newt;
	}
	template<typename T>
	T * Ptr()
	{
		return (T*)ptr();//m_Bytes.size() > 0 ? (T*)&m_Bytes[0] : NULL;
	}

	void * ptr()
	{
		return m_RefPtr ? m_RefPtr : (m_bytes.size() > 0 ? &m_bytes[0] : NULL);
		//return BytePtr();// m_Bytes.size() > 0 ? &m_Bytes[0] : NULL;
	}
	BYTE * BytePtr() {
		return m_RefPtr ? (BYTE*)m_RefPtr : (m_bytes.size() > 0 ? &m_bytes[0] : NULL);
	};
	size_t SizeInBytes()
	{
		return m_RefPtr ? m_ElCount * m_ElSize : m_bytes.size();
	}
	ByteVec & ByteVecRef()
	{
		assert(!m_RefPtr);
		return m_bytes;
	}

	bool bEmpty()
	{
		if (m_RefPtr) return false;
		else
		return m_bytes.size() == 0;
	}
	template<typename T>
	bool AddOn(double  rhs)
	{
		if (sizeof(T) != m_ElSize)
			return false;

		T* p1 = (T*)DataPtr(0);
		for (int i = 0; i < m_ElCount; i++)
		{
			p1[i] += rhs;
		}
		return true;
	}
	//��ͬά�ȼӷ�
	template<typename T>
	bool AddOn(CDGMulDimData & rhs)
	{
		if (rhs.m_dim != m_dim || rhs.m_type != m_type || sizeof(T) != m_ElSize)
			return false;

		T* p1 = (T*)DataPtr(0), *p2=(T*)rhs.DataPtr(0);
		for (int i = 0; i < m_ElCount; i++)
		{
			p1[i] += p2[i];
		}
		return true;
	}
	template<typename T>
	bool MulOn(CDGMulDimData & rhs)
	{
		if (rhs.m_dim != m_dim || rhs.m_type != m_type || sizeof(T) != m_ElSize)
			return false;

		T* p1 = (T*)DataPtr(0), *p2 = (T*)rhs.DataPtr(0);
		for (int i = 0; i < m_ElCount; i++)
		{
			p1[i] *= p2[i];
		}
		return true;
	}
	template<typename T>
	bool Max(double  rhs)
	{
		if (sizeof(T) != m_ElSize)
			return false;

		T* p1 = (T*)DataPtr(0);
		for (int i = 0; i < m_ElCount; i++)
		{
			p1[i] = max((double)p1[i] ,rhs);
		}
		return true;
	}
	template<typename T>
	bool Min(double  rhs)
	{
		if (sizeof(T) != m_ElSize)
			return false;

		T* p1 = (T*)DataPtr(0);
		for (int i = 0; i < m_ElCount; i++)
		{
			p1[i] = min((double)p1[i], rhs);
		}
		return true;
	}
	template<typename T>
	bool MulOn(double  rhs)
	{
		if(sizeof(T) != m_ElSize)
			return false;

		T* p1 = (T*)DataPtr(0);
		for (int i = 0; i < m_ElCount; i++)
		{
			p1[i] *= rhs;
		}
		return true;
	}
	bool ReShapeInPlace(DimVec newdim);
	
	template<typename T>
	void SumOnDim(int nLowDim,CDGMulDimData & sumed)
	{
		DimVec outdim; outdim.insert(outdim.end(), m_dim.begin(), m_dim.begin() + nLowDim);
		DimVec sDim; sDim.insert(sDim.end(), m_dim.begin() + nLowDim, m_dim.end());
		
		int nframe = DimElCount(sDim);
		int nFrameEl = DimElCount(outdim);
		sumed = CDGMulDimData(outdim, m_type,m_Acc.m_DataMajor);
		
		CDimDataAccer acc(sDim, m_Acc.m_DataMajor);
		acc.TraverseNoRecusive((TraverseFunc)
			[&](DimIndex & index)->bool
		{
			T* pf = (T*)SelFrameByDim(nLowDim, index);
			T* ppf = sumed.Ptr<T>();
			for (int iel =0;iel< nFrameEl;iel++)
			{
				ppf[iel] += pf[iel];
			}
			return true;
		}
		);
	}
	template<typename T>
	bool DivOn(CDGMulDimData & rhs)
	{
		if (rhs.m_dim != m_dim || rhs.m_type != m_type || sizeof(T) != m_ElSize)
			return false;

		T* p1 = (T*)DataPtr(0), *p2 = (T*)rhs.DataPtr(0);
		for (int i = 0; i < m_ElCount; i++)
		{
			p1[i] /= p2[i];
		}
		return true;
	}
	template<typename T>
	bool SubOn(CDGMulDimData & rhs)
	{
		if (rhs.m_dim != m_dim || rhs.m_type != m_type || sizeof(T) != m_ElSize)
			return false;

		T* p1 = (T*)DataPtr(0), *p2 = (T*)rhs.DataPtr(0);
		for (int i = 0; i < m_ElCount; i++)
		{
			p1[i] -= p2[i];
		}
		return true;
	}

	template<typename T>
	Eigen::Map<T> SelMatByDimXY_Maped(DimIndex FrameIndex)
	{
		Eigen::Map<T> ret((T::Scalar*) SelMatByDimXY(FrameIndex), m_dim[0], m_dim[1]);
		return ret;	
	}
	template<typename T>
	Eigen::Map<T> SelMatByDimXY_Maped_AsVec(DimIndex FrameIndex)
	{
		Eigen::Map<T> ret((T::Scalar*) SelMatByDimXY(FrameIndex), m_dim[0] * m_dim[1]);
		return ret;
	}
	
	void* SelMatByDimXY(DimIndex FrameIndex);
	void* SelFrameByDim(int nSelDim, DimIndex FrameIndex);
	void SelFrameByDim(int nSelDim, DimIndex FrameIndex,CDGMulDimData & out);
	float SumAsFloat();
	double SumAsDouble();
	double Sum();

	bool BsxFuncLowDimDup(BsxOperator op, CDGMulDimData & rhs, CDGMulDimData& outres);
	bool BsxFuncLowDim_InPlace(BsxOperator op, CDGMulDimData & rhs);

	bool BsxFuncHighDimDup(BsxOperator op, CDGMulDimData & rhs, int nHighDupDim, CDGMulDimData& outres);

	bool BsxFuncHighDim_InPlace(BsxOperator op, CDGMulDimData & rhs, int nHighDupDim);

	bool BsxFunction(BsxOperator op,CDGMulDimData & rhs, CDGMulDimData& result);
	byte* DataPtr(DimIndex ind);
	byte* DataPtr(int pos);

	DGDataType m_type;
	DimVec m_dim;	
	CDimDataAccer m_Acc;
	int m_ElSize{0};
	int m_ElCount{0};

	void Print();

protected:
	ByteVec m_bytes;
};


BsxfFuncPointer GetBxsfFunc(BsxOperator op, DGDataType type);


template<typename Tv, typename TOuVal>
void Grid3D(Tv d1, Tv d2, Tv d3, CDGMulDimData &x1, CDGMulDimData &x2, CDGMulDimData &x3, DGDataType tartype)
{
	x3 = x2 = x1 = CDGMulDimData({ (int)d1.size(), (int)d2.size(),(int)d3.size() }, tartype, DataMajor_FStyle_ColumnMajor);

	for (int ix1 = 0; ix1 < d1.size(); ix1++)
	{
		for (int ix2 = 0; ix2 < d2.size(); ix2++)
		{
			for (int ix3 = 0; ix3 < d3.size(); ix3++)
			{
				DimIndex di = { ix1, ix2,ix3 };
				*(TOuVal*)x1.DataPtr(di) = d1(ix1);
				*(TOuVal*)x2.DataPtr(di) = d2(ix2);
				*(TOuVal*)x3.DataPtr(di) = d3(ix3);
			}
		}
	}
}

sitk::Image ResampleImageFromMD(DimVec newdim, CDGMulDimData & md);
sitk::Image ResampleImageFromMD(DimVec newdim, CDGMulDimData & md, MAT_TYPE & m);


//int DataCount();
//BsxfFuncPointer GetBxsfFunc(BsxOperator op, DGDataType type);
/*for (int idim = 0; idim < nLowDim; idim++)
		{
			outdim.push_back(m_dim[idim]);
		}*/
		/*DimVec sDim;
		for (int idim = outdim.size(); idim < m_dim.size(); idim++)
		{
			sDim.push_back(m_dim[idim]);
		}*/
/*MAT_TYPE ToMat();
MAT_TYPE_F ToMatF();
VECTOR_TYPE ToVector();
VECTOR_TYPE_F ToVectorF();*/

//template<typename T>
//T SelMatByDimXY_Casted(DimIndex FrameIndex)
//{
//	T tout=T::Zero(m_dim[0], m_dim[1]);

//	void * ptr = SelMatByDimXY(FrameIndex);
//	
//	castd
//	//Eigen::Map<T> ret();
//	return tout;
//}