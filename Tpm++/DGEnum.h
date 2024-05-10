#pragma once

enum DGDataType
{
	DGDataType_Float,
	DGDataType_Double,
	DGDataType_Int,
	DGDataType_FixLenChar128,
	DGDataType_BOOL,
	DGDataType_Vec3,	
	DGDataType_Matrix4x4,
	DGDataType_FixLenChar1024,
	DGDataType_PixelFloat,
	DGDataType_UnSignedChar,
	DGDataType_UnSignedShort,
	DGDataType_SignedShort,
	DGDataType_UnSignedInt,
	DGDataType_Int64,
	DGDataType_UnsignedInt64,
	DGDataType_String,
	DGDataType_Act,
	DGDataType_PackedKeyValPair,
	DGDataType_PackedItkImage,
	DGDataType_Char,
	DGDataType_PackedItkTransform,
	DGDataType_ComplexDouble,
	DGDataType_ComplexFloat,
	//DGDataType_SubAct,
	DGDataType_SubTypeMask1 = 0x0000ff00,
	DGDataType_SubTypeMask2 = 0x00ff0000,
	DGDataType_SubTypeMask3 = 0xff000000,
	DGDataType_SubTypeMaskAll = 0xffffff00,
	DGDataType_MajorMask = 0x000000ff,

	
	DGDataType_dep= 0x01000000,
	DGDataType_Act_old = 0xfffffffd,
	DGDataType_NotDefined = 0xfffffffe,
	DGDataType_Max = 0xffffffff
};



typedef  std::vector<DGDataType> DGDataTypeVec;
#define MAXDIM 10

inline bool IsComplex(DGDataType type)
{
	return (type == DGDataType_ComplexDouble) || (type == DGDataType_ComplexFloat);
}

#define SwitchNCallAllDGDataTypeNRealType(__FuncCall,__DgType){\
switch (__DgType)\
{\
	case DGDataType_Float:\
		__FuncCall(float);\
		break;\
	case DGDataType_Double:\
		__FuncCall(double);\
		break;\
	case DGDataType_Int:\
		__FuncCall(int);\
		break;\
	case DGDataType_BOOL:\
		__FuncCall(BOOL);\
		break;\
	case DGDataType_UnSignedChar:\
		__FuncCall(uint8_t);\
		break;\
	case DGDataType_Char:\
		__FuncCall(int8_t);\
		break;\
	case DGDataType_UnSignedShort:\
		__FuncCall(uint16_t);\
		break;\
	case DGDataType_SignedShort:\
		__FuncCall(int16_t);\
		break;\
	case DGDataType_UnSignedInt:\
		__FuncCall(unsigned int);\
		break;\
	case DGDataType_Int64:\
		__FuncCall(int64_t);\
		break;\
	case DGDataType_UnsignedInt64:\
		__FuncCall(uint64_t);\
		break;\
}\
}

inline int DGDataSize(DGDataType type)
{
	int s=0;
	SwitchNCallAllDGDataTypeNRealType(s=sizeof, type);

	return s;
}


template<typename T>
bool CastToDGData(T & tIn, void * pval, DGDataType type)
{
#define CastByTypeNRet(__type) {*((__type*)pval) = (__type)tIn;return true;}
SwitchNCallAllDGDataTypeNRealType(CastByTypeNRet, type);
return false;
}

template<typename T>
bool DGDataCastToT(void * pval, DGDataType type, T & tOut)
{	
#define CastByTypeToTNRet(__type) {tOut = (T)(*(__type*)pval);return true;}	
	SwitchNCallAllDGDataTypeNRealType(CastByTypeToTNRet, type);	
return false;
}


typedef std::vector<int> IntVEC;
