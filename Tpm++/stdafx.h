#pragma once

#include <Eigen/Dense> 
#include <Eigen/Core>
#include <vector>
#include <mutex>
#include <string>
#include "struct2x/struct2x.h"

#define Div255 1.0f/255.0f
#define BYTE_TO_FLOAT(bt) (bt*Div255 )
#define FLOAT_TO_BYTE(fl) std::max((BYTE)0,(BYTE)std::min((BYTE)(fl*255) , (BYTE)255))

class CriticalSectionWrapper {
public:
    std::mutex mtx;
};

using CRITICAL_SECTION = CriticalSectionWrapper;
inline void InitializeCriticalSection(CRITICAL_SECTION* cs) {}
inline void EnterCriticalSection(CRITICAL_SECTION* cs) {cs->mtx.lock();}
inline void LeaveCriticalSection(CRITICAL_SECTION* cs) {cs->mtx.unlock();}
inline void DeleteCriticalSection(CRITICAL_SECTION* cs) {}

class CriticalSectionScoper
{
public:
	CriticalSectionScoper(CRITICAL_SECTION &cs)
	:m_cs(cs)
    { EnterCriticalSection(&m_cs); }

	~CriticalSectionScoper()
    {LeaveCriticalSection(&m_cs);};
	
	CRITICAL_SECTION & m_cs;
};

struct VERTEX3D
{
    public:
    VERTEX3D(float x, float y, float z)
    {
        m_v[0] = x;
        m_v[1] = y;
        m_v[2] = z;
    }
    VERTEX3D()
    {
        m_v[0] = 0;
        m_v[1] = 0;
        m_v[2] = 0;
    }
    

    float & x() { return m_v[0]; }
    float & y() { return m_v[1]; }
    float & z() { return m_v[2]; }

    float m_v[3]={0,0,0};
};

enum DataMajor
{
	DataMajor_CStyle_RowMajor = 0x00000000,
	DataMajor_FStyle_ColumnMajor = 0x00000001
};

//typedef unsigned char byte;
using byte = std::byte;//unsigned char;
typedef byte BYTE;
typedef unsigned int UINT;
typedef int BOOL;
using namespace std;
#include <cstdarg>

inline std::string FormatedString(const char * format, ...)
{
	UINT len = (UINT)max((UINT)1024,(UINT)strlen(format)*2);

	char *s =new char[len];
	va_list v;
	va_start(v, format);
	vsprintf_s(s,1024 ,format, v);
	va_end(v);	
	string str = s;
	delete[](s);
	return str;
}

typedef std::vector<int> DimVec;
typedef std::vector<int64_t> DimVec64;
typedef std::vector<unsigned int> UIndexes;
typedef std::vector<double> DoubleVec;

typedef std::vector<unsigned int> SizeVec;
typedef std::vector<BYTE> ByteVec;
typedef std::vector<float> FloatVec;
typedef std::vector<ByteVec*> ByteVecPtrs;
typedef std::vector<FloatVec> FloatVecVec;
typedef std::vector<int> DimIndex;
#define mwSize unsigned int
#define mwSignedIndex int
#define mwIndex unsigned int

#include <algorithm>
#include <set>
#define OUT
#define IN
#include "DGEnum.h"

#include "SimpleITK.h"
#define PI 3.14159265358979
namespace sitk = itk::simple;
#include <chrono>
#include <cstdint>
typedef uint32_t DWORD ;

inline DWORD timeGetTime() {
    static const auto start_time = std::chrono::steady_clock::now();
    auto current_time = std::chrono::steady_clock::now();
    return (DWORD)std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start_time).count();
}

//using namespace sitk;
class  CTimeCounter
{
public:
	CTimeCounter(std::string str)
    :m_str(str)
    {
        m_LifeTime = timeGetTime();
    }
public:
	~CTimeCounter(void)
    {
        m_LifeTime = timeGetTime() -m_LifeTime;
        if(m_LifeTime >10)
        {	
            printf("%s time = %f s \n" ,m_str.c_str(), m_LifeTime/1000.0f );
        
        }
    }
	DWORD m_LifeTime;
	std::string m_str;
};

typedef std::vector<sitk::Image> ItkImageVec;


inline void ErrOutput(string s1,string s2)
{
    printf("%s,%s \n",s1.c_str() ,s2.c_str());
}

inline DimVec SizeToDim(SizeVec & size)
{
	DimVec ret;ret.resize(size.size());
	for (int i=0;i<size.size();i++)
	{
		ret[i] = size[i];
	}
	return ret;
}

inline SizeVec DimToSize(DimVec & dim)
{
	SizeVec ret;ret.resize(dim.size());
	for (int i = 0;i < dim.size();i++)
	{
		ret[i] = dim[i];
	}
	return ret;
}


typedef std::vector<int> DimVec;
typedef  std::vector<int64_t> DimVec64;

typedef std::vector<unsigned int> SizeVec;

inline UINT Dim3index(int ix, int iy, int iz, DimVec &dim)
{
	assert(dim.size() == 3 && ix < dim[0] && iy < dim[1] && iz < dim[2]);
	return  iz * (dim[0] * dim[1]) + iy * dim[0] + ix;
}

inline bool Dim3IndexInRange(DoubleVec &index, DimVec &dim)
{
	return index[0] >= 0 && index[1] >= 0 && index[2] >= 0 &&
		index[0] <= (dim[0] - 1) && index[1] <= (dim[1] - 1) && index[2] <= (dim[2] - 1);
}
inline bool Dim3IndexInRange(DimVec &index, DimVec &dim)
{
	return index[0] >= 0 && index[1] >= 0 && index[2] >= 0 &&
		index[0] <= (dim[0] - 1) && index[1] <= (dim[1] - 1) && index[2] <= (dim[2] - 1);
}
inline bool Dim3IndexInRange(DimVec64 &index, DimVec &dim)
{
	return index[0] >= 0 && index[1] >= 0 && index[2] >= 0 &&
		index[0] <= (dim[0] - 1) && index[1] <= (dim[1] - 1) && index[2] <= (dim[2] - 1);
}
typedef std::vector<int> DimIndex;
inline DimIndex IndexToDim3(int index, DimVec &dim)
{
	DimIndex ret = {
		index % (dim[0] * dim[1]) % dim[0],
		index % (dim[0] * dim[1]) / dim[0],
		index / (dim[0] * dim[1])
	};

	return ret;
}
inline void ClampIndex(DimIndex & index,DimVec & dim)
{
	index[0] = clamp<int>(index[0],0,dim[0] - 1);
	index[1] = clamp<int>(index[1], 0, dim[1] - 1);
	index[2] = clamp<int>(index[2], 0, dim[2] - 1);
}

#include "CasterVec.h"
#include "DimDataAccer.h"
typedef std::vector<ByteVec> ByteVecs;
#include "libjson/libjson.h"