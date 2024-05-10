#pragma once

typedef std::function<bool(DimIndex&)> TraverseFunc;

class CDataElTraversor
{
public:
	virtual bool OnTraverse(DimIndex & index) = 0;

};

inline int DimElCount(DimVec & dim)
{
	int count = 1;
	for(auto d:dim)
	{
		count *= d;
	}
	return count;
}
class CDimDataAccer
{
public:
	

	//static DataMajor s_DataMajor;

	CDimDataAccer(DimVec Dim, DataMajor DataMajor)
		:m_Dim(Dim), m_DataMajor(DataMajor)
	{
		m_TotalElCount = DimElCount(Dim);
		m_StridesElPos = GetStrides(1);
	}

	DimIndex m_StridesElPos;
	DataMajor m_DataMajor{ DataMajor_FStyle_ColumnMajor };
	DimIndex GetStrides(int elsize);
	bool NextDimEl(DimIndex & index);

	int64_t GetElPos(DimIndex & index);

	//DimVec GetSlicedDim(UINT seldim);
	DimVec DropSelectedDim(UINT seldim);

	DimIndex FindByContent(ByteVec & keycontent, DGDataType type, ByteVec & indata, bool bSubStrCmp);
	void Traverse(CDataElTraversor & tra);
	void Traverse(TraverseFunc functra);

	void TraverseNoRecusive(TraverseFunc functra);
	void TraverseNoRecusiveOmp(TraverseFunc functra);

	bool ForEach(UINT currdim, DimIndex & parentindex, TraverseFunc & tra);
	bool ForEach(UINT currdim, DimIndex & parentindex, CDataElTraversor & tra);
	DimVec m_Dim;
	UINT m_TotalElCount{ 0 };


protected:
	int GetElPosCStyle(const DimIndex & work_Index);
	int GetElPosFStyle(const DimIndex & work_Index);
};
