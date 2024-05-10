#pragma once

typedef std::function<void(int, BYTE&)> UpdateMaskFunc;

struct sMask
{
	void TraverseMsk(UpdateMaskFunc func)
	{
		for (int imsk = 0; imsk < m_Msk.size(); imsk++)
		{
			func(imsk, m_Msk[imsk]);
		}
	}

	sMask()
	{

	}
	sMask(int nSize)
		:m_Msk(nSize)
	{
		if (nSize > 0) 	std::memset(&m_Msk[0], 0, m_Msk.size());
	}
	void SetMask(int imsk, bool bmask)
	{
		m_Msk[imsk] = (byte)bmask;
	}
	void SetMaskAndNM(int imsk, bool bmask)
	{
		m_Msk[imsk] = (byte)bmask;
		m_NM += bmask;
	}

	void MaskOnInfinite(VECTOR_TYPE & d)
	{
		int s = (int)d.size();
		m_NM = 0;
		for (int i = 0; i < s; i++)
		{
			m_Msk[i] = (byte)!std::isfinite(d[i]);
			m_NM += (bool)m_Msk[i];
		}
	}
	void MaskOnNoLess(double dcompare, VECTOR_TYPE & d)
	{
		m_NM = 0;
		int s = (int)d.size();
		for (int i = 0; i < s; i++)
		{
			m_Msk[i] = (byte)(d[i] >= dcompare);
			m_NM += (bool)m_Msk[i];
		}
	}
	void MaskOnLarger(double dcompare, VECTOR_TYPE & d)
	{
		m_NM = 0;
		int s =(int) d.size();
		for (int i = 0; i < s; i++)
		{
			m_Msk[i] = (byte)(d[i] > dcompare);
			m_NM += (bool)m_Msk[i];
		}
	}

	void MaskOnFinite(VECTOR_TYPE & d)
	{
		m_NM = 0;
		int s = (int)d.size();
		for (int i = 0; i < s; i++)
		{
			m_Msk[i] =(byte) std::isfinite(d[i]);
			m_NM += (bool)m_Msk[i];
		}
	}

	inline bool bMasked(int & i)
	{
		return (bool)m_Msk[i];
	}
	void Mask3Vector(VECTOR_TYPE& vvector1, VECTOR_TYPE& vvector2, VECTOR_TYPE& vvector3);
	void Mask3Mat(MAT_TYPE & mat1, MAT_TYPE & mat2, MAT_TYPE & mat3, OUT VECTOR_TYPE & v1, OUT VECTOR_TYPE & v2, OUT VECTOR_TYPE & v3);
	bool MaskMatchAndAssign(VECTOR_TYPE& target, VECTOR_TYPE& src);
	void MaskAssign(VECTOR_TYPE& target, double val);
	void MaskAssign(MAT_TYPE& target, VECTOR_TYPE &val);
	VECTOR_TYPE_F MaskVector(VECTOR_TYPE_F& vvector);
	VECTOR_TYPE MaskVector(VECTOR_TYPE& vvector);
	VECTOR_TYPE MaskMat(MAT_TYPE & mat);
	VECTOR_TYPE_F MaskMat(MAT_TYPE_F & mat);
	inline bool HasMsk()
	{
		return m_NM > 0;
	}
	void ReComputeNM();
	ByteVec m_Msk;
	int m_NM{ 0 };
};
