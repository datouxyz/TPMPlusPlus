#pragma once

#define CastByTypeToTVec(__type) {\
	for (int iel = 0; iel < nel; iel++)\
	{\
		tOut[iel] = (T)(*(__type*)&pbval[iel*elsize]);\
	}\
	}

template<typename T>
bool DGDataCastToT_VEC(int nel, void * pval, DGDataType type, T * tOut)
{
	int elsize = DGDataSize(type);
	byte * pbval = (byte*)pval;

	SwitchNCallAllDGDataTypeNRealType(CastByTypeToTVec, type);
	return false;
}

#define CastByType_VEC(__type) 	{\
		for (int iel = 0; iel < nel; iel++)\
		{\
			*((__type*)&pbval[iel*elsize]) = (T)tIn[iel]; \
		}\
	}

template<typename T>
bool CastToDGData_VEC(T * tIn, int nel, void * pval, DGDataType type)
{
	int elsize = DGDataSize(type);
	byte * pbval = (byte*)pval;

	SwitchNCallAllDGDataTypeNRealType(CastByType_VEC, type);
	return false;
}
