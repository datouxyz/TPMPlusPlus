#pragma once

MAT_TYPE Derivs(MAT_TYPE &MF, VECTOR_TYPE &P/**/, MAT_TYPE &MG);

VECTOR_TYPE CumSum(VECTOR_TYPE & input);

template<class Tv,typename Scalar>
Tv BuildVector(Scalar be, Scalar  en, Scalar step)
{
	int n = (en - be) / step + 1;
	Tv r(n);

	for (int i = 0; i < n; i++)
	{
		r[i] = be + step * i;
	}
	return r;
}
template<typename Tv, typename Tm>
void Grid2D(Tv d1, Tv d2, Tm &x1, Tm &x2)
{
	x2 = x1 = Tm::Zero(d1.size(), d2.size());

	for (int ix1 = 0; ix1 < d1.size(); ix1++)
	{
		for (int ix2 = 0; ix2 < d2.size(); ix2++)
		{
			x1(ix1, ix2) = d1(ix1);
			x2(ix1, ix2) = d2(ix2);
		}
	}
}


MAT_TYPE DCT_MTX(int N, int K, VECTOR_TYPE &n);
//NXN 
MAT_TYPE DCT_MTX(int N);
//NXK
MAT_TYPE DCT_MTX(int N, int K);

MAT_TYPE DCT_MTX_DIFF(int N, int K, VECTOR_TYPE &n);
MAT_TYPE DCT_MTX_DIFF(int N);
MAT_TYPE DCT_MTX_DIFF(int N, int K);

MAT_TYPE DCT_MTX_DIFF2(int N, int K, VECTOR_TYPE &n);
MAT_TYPE DCT_MTX_DIFF2(int N);
MAT_TYPE DCT_MTX_DIFF2(int N, int K);

#include "DGMulDimData.h"

inline auto MDToMat(CDGMulDimData & md)
{
	if (md.m_dim.size() != 2 || md.m_type != DGDataType_Double) return Eigen::Map<MAT_TYPE>(NULL, 0, 0);

	return Eigen::Map<MAT_TYPE>(md.Ptr<double>(), md.m_dim[0], md.m_dim[1]);	
}

inline CDGMulDimData MatToMD(MAT_TYPE & m)
{
	return	CDGMulDimData({ (int)m.rows(),(int)m.cols() }, DGDataType_Double,DataMajor_FStyle_ColumnMajor, m.data());
}

inline CDGMulDimData VecToMD(VECTOR_TYPE & v)
{
	return CDGMulDimData({ (int)v.size() }, DGDataType_Double, DataMajor_FStyle_ColumnMajor, v.data());
}

inline CDGMulDimData MatFToMD(MAT_TYPE_F & m)
{
	return  CDGMulDimData({ (int)m.rows(),(int)m.cols() }, DGDataType_Float, DataMajor_FStyle_ColumnMajor, m.data());
}

inline MAT_TYPE_F MDToMatF(CDGMulDimData & md)
{
	if (md.m_dim.size() != 2 || md.m_type != DGDataType_Float) return  Eigen::Map<MAT_TYPE_F>(NULL, 0, 0);
	/*MAT_TYPE_F ret =  Eigen::Map<MAT_TYPE_F>(md.Ptr<float>(), md.m_dim[1], md.m_dim[0]);
	ret = ret.transpose();*/
	MAT_TYPE_F ret = Eigen::Map<MAT_TYPE_F>(md.Ptr<float>(), md.m_dim[0], md.m_dim[1]);
	return ret;
}

inline auto MDToVector(CDGMulDimData & md)
{
	if (md.m_dim.size() != 1 || md.m_type != DGDataType_Double) return Eigen::Map<VECTOR_TYPE>(NULL, 0);
	return Eigen::Map<VECTOR_TYPE>(md.Ptr<double>(), md.m_dim[0]);	
}
inline auto MDToVectorF(CDGMulDimData & md)
{
	if (md.m_dim.size() != 1 || md.m_type != DGDataType_Float) return Eigen::Map<VECTOR_TYPE_F>(NULL, 0);
	return Eigen::Map<VECTOR_TYPE_F>(md.Ptr<float>(), md.m_dim[0]);	
}

void WriteMDAsImage(string path, CDGMulDimData & md);

void Write2DMATDouble(string path, string paramname, int32_t cols, int32_t rows, double* pdata);

void Write2DMATDouble(string path, string paramname, MAT_TYPE & m);

void Write2DMATDouble(string path, string paramname, VECTOR_TYPE & v);

void Write2DMATFloat(string path, string paramname, int32_t cols, int32_t rows, float* pdata);

void Write2DMATFloat(string path, string paramname, MAT_TYPE_F & m);

void Write2DMATFloat(string path, string paramname, VECTOR_TYPE_F & v);

template<typename T>
void WriteNDMat(string path, string paramname, CDGMulDimData & md)
{
	TinyMATWriterFile* tinymat = TinyMATWriter_open(path.c_str(), 0, md.m_Bytes.size());


	TinyMATWriter_writeMatrixND_colmajor(tinymat, paramname.c_str(), (T*)&md.m_Bytes[0],
		&md.m_dim[0], md.m_dim.size());

	TinyMATWriter_close(tinymat);//���Զ��ͷ�mat
}
void WriteNDMatFloat(string path, string paramname, CDGMulDimData & md);
void WriteNDMatDouble(string path, string paramname, CDGMulDimData & md);

template<typename Tm,typename Tv>
Tm RepMat(Tv & v, int r, int c)
{
	Tm ret = Tm::Zero(v.rows()*r, v.cols()*c);
	for (int ir = 0; ir < r; ir++)
		for (int ic = 0; ic < c; ic++)
		{
			ret.block(ir * v.rows(), ic*v.cols(), v.rows(), v.cols()) =
				v;
		}
	return ret;
}



MAT_TYPE RepMat(VECTOR_TYPE & v, int r, int c);
MAT_TYPE RepMat(MAT_TYPE & m, int r, int c);
MAT_TYPE RDivide(MAT_TYPE &A, MAT_TYPE &B);

bool Spdiags(MAT_TYPE & B, VECTOR_TYPE d, int nrow, int ncol, MAT_SPARSE_TYPE & ret);

void SmoothHistrogram(MAT_TYPE & t0, OUT MAT_TYPE &sigOut);
void SmoothHistrogram(MAT_TYPE & t0, VECTOR_TYPE &lam, OUT MAT_TYPE &sigOut);
void SmoothHistrogram(MAT_TYPE & t0, VECTOR_TYPE &lam, OUT MAT_TYPE &sigOut, OUT MAT_TYPE & alpha);

IntVEC FindInVec(VECTOR_TYPE &lkp, double key);
int MaxVec(IntVEC & lkp);

template<typename T>
DGDataType EiGenElDGType()
{	
	if (typeid(T::Scalar) == typeid(float))
	{
		return DGDataType_Float;
	}
	if (typeid(T::Scalar) == typeid(double))
	{
		return DGDataType_Double;
	}
	return DGDataType_NotDefined;
};

template<typename TData, typename TInd>
void AccessByIndex(TInd &ind, TData &data, OUT TData & outarray)
{
	outarray.resize(ind.size());
	for (int i = 0; i < ind.size(); i++)
	{
		outarray[i] = data[ind[i] - 1];
	}
}
