#include "stdafx.h"

#include "DimDataAccer.h"



bool CDimDataAccer::NextDimEl(DimIndex & index)
{
	for (int idim = 0; idim < m_Dim.size(); idim++)
	{
		if (index[idim] < m_Dim[idim] - 1)
		{
			index[idim]++;
			break;
		}
		else if (idim < m_Dim.size() - 1)
		{
			index[idim] = 0;
		}
		else
		{
			return false;
		}
	}
	return true;
}

void CDimDataAccer::TraverseNoRecusive(TraverseFunc functra)
{
	DimIndex dindex(m_Dim.size(), 0);
	if (m_Dim.size() <= 0)
	{
		functra(DimIndex());
		return;
	}

	do
	{
		//��һά��ѭ������
		for (int i0 = 0; i0 < m_Dim[0]; i0++)
		{
			DimIndex d0 = dindex;
			d0[0] = i0;
			functra(d0);
		}
		dindex[0] = m_Dim[0] - 1;
		//functra(dindex);
	} while (NextDimEl(dindex));
}
void CDimDataAccer::TraverseNoRecusiveOmp(TraverseFunc functra)
{
	DimIndex dindex(m_Dim.size(), 0);
	if (m_Dim.size() <= 0)
	{
		functra(DimIndex());
		return;
	}

	do
	{
		//��һά��ѭ������
#pragma omp parallel for schedule(static)
		for (int i0 = 0; i0 < m_Dim[0]; i0++)
		{
			DimIndex d0 = dindex;
			d0[0] = i0;
			functra(d0);
		}
		dindex[0] = m_Dim[0] - 1;
		//functra(dindex);
	} while (NextDimEl(dindex));
}

void CDimDataAccer::Traverse(TraverseFunc functra)
{
	if (m_Dim.size() == 0)
	{
		functra(DimIndex());
		return;
	}
	int currdim = 0;
	DimIndex parentindex;
	//functra(parentindex);
	ForEach(currdim, parentindex, functra);
}
bool CDimDataAccer::ForEach(UINT currdim, DimIndex & parentindex, TraverseFunc & tra)
{
	if (currdim >= m_Dim.size())
	{
		return true;
	}
	DimIndex newindex = parentindex;
	newindex.push_back(0);
	int & last = newindex[newindex.size() - 1];
	for (int iobj = 0; iobj < m_Dim[currdim]; iobj++)
	{
		last = iobj;
		if (currdim < m_Dim.size() - 1)
		{
			if (!ForEach(currdim + 1, newindex, tra))
				return false;
		}
		else
			if (!tra(newindex))
				return false;
	}
	return true;
}

void CDimDataAccer::Traverse(CDataElTraversor & tra)
{
	if (m_Dim.size() == 0)
	{
		tra.OnTraverse(DimIndex());
		return;
	}
	int currdim = 0;
	DimIndex parentindex;
	ForEach(currdim, parentindex, tra);
}
bool CDimDataAccer::ForEach(UINT currdim, DimIndex & parentindex, CDataElTraversor & tra)
{
	if (currdim >= m_Dim.size())
	{
		return true;
	}
	DimIndex newindex = parentindex;
	newindex.push_back(0);
	int & last = newindex[newindex.size() - 1];
	for (int iobj = 0; iobj < m_Dim[currdim]; iobj++)
	{
		last = iobj;
		if (currdim < m_Dim.size() - 1)
		{
			if (!ForEach(currdim + 1, newindex, tra))
				return false;
		}
		else
			if (!tra.OnTraverse(newindex))
				return false;
	}
	return true;
}
DimIndex CDimDataAccer::GetStrides(int elsize)
{
	DimIndex strides(m_Dim.size(), elsize);	
	switch (m_DataMajor)
	{
	case DataMajor_CStyle_RowMajor:
		for (int i = m_Dim.size() - 1; i > 0; --i)
		{
			strides[i - 1] = strides[i] * m_Dim[i];
		}
		break;
	case DataMajor_FStyle_ColumnMajor:
		for (size_t i = 1; i < m_Dim.size(); ++i) {
			strides[i] = strides[i - 1] * m_Dim[i - 1];
		}

		break;
	default:
		break;
	}

	return strides;
}


int CDimDataAccer::GetElPosFStyle(const DimIndex & work_Index)
{
	if (work_Index.size() != m_Dim.size()) return -1;

	const DimVec &work_Dim = m_Dim;
	int nCurrFrameEl = m_TotalElCount;
	int iCurrStartEl = 0;
	for (int idim = work_Dim.size() - 1; idim >= 0; idim--)
	{
		int currdim = work_Dim[idim];
		int curri = work_Index[idim];
		if (curri >= currdim) return -1;
		nCurrFrameEl /= currdim;
		iCurrStartEl += nCurrFrameEl * curri;
	}
	return iCurrStartEl;
}

int CDimDataAccer::GetElPosCStyle(const DimIndex & work_Index)
{
	DimVec &work_Dim = m_Dim;
	int nCurrFrameEl = m_TotalElCount;
	int iCurrStartEl = 0;
	for (UINT idim = 0; idim < work_Dim.size(); idim++)
	{
		int currdim = work_Dim[idim];
		int curri = work_Index[idim];
		if (curri >= currdim) return -1;
		nCurrFrameEl /= currdim;
		iCurrStartEl += nCurrFrameEl * curri;
	}
	return iCurrStartEl;
}


DimVec CDimDataAccer::DropSelectedDim(UINT seldim)
{
	if (seldim >= m_Dim.size())
	{
		ErrOutput("DownDimension", "seldim out of range");
		return DimVec();
	}
	DimVec ret;
	for (UINT idim = 0; idim < m_Dim.size(); idim++)
	{

		/*if (idim == seldim)
		{
			ret.push_back(1);
		}*/
		if (idim == seldim)
			continue;
		else
		{
			ret.push_back(m_Dim[idim]);
		}
	}
	return ret;
}
int64_t CDimDataAccer::GetElPos(DimIndex & index)
{
	if (index.size() != m_StridesElPos.size())
		return -1;
	uint64_t retPos = 0;
	for (int i = 0; i < m_StridesElPos.size(); i++)
	{
		retPos += index[i] * m_StridesElPos[i];
	}
	return retPos;
}
DimIndex CDimDataAccer::FindByContent(ByteVec & keycontent, DGDataType type, ByteVec & indata, bool bSubStrCmp)
{
	DimIndex ret;
	int elsize = DGDataSize(type);
	int nEl = indata.size() / elsize;
	for (int i = 0; i < nEl; i++)
	{
		ByteVec eldata;
		eldata.resize(elsize);
		memcpy(&eldata[0], &indata[i*elsize], elsize);
		if (bSubStrCmp)
		{
			if (std::search(eldata.begin(), eldata.end(), keycontent.begin(), keycontent.end()) != eldata.end())
			{
				ret.push_back(i);
			}
		}
		else
		{
			if (eldata == keycontent)
			{
				ret.push_back(i);
			}
		}
	}
	return ret;
}