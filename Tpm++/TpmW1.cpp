#include "stdafx.h"
#include <Eigen/Dense>

#include "TissueProblityMap.h"

#include "EigenMatFunc.h"
#include "Conv2.h"
#include "shoot_boundary.h"
#include "llVar.h"
#include <Eigen/SVD>
#include "EigenMatFunc.h"
#include "TPMInternalFuncs.h"

#define MAXCLASSES 1024

void MarkrofField(IN OUT CDGMulDimData & vPByte, IN CDGMulDimData& QFl, IN CDGMulDimData& vG, IN VECTOR_TYPE &vx2)
{
	unsigned int i;	
	DimVec dm = { 0,0,0,0 };
	float w[3];
	int code = 0;

	if (vPByte.m_dim.size() != QFl.m_dim.size())
	{
		ErrOutput("MarkrofField", u8"");
		return;
	}

	for (i = 0; i < vPByte.m_dim.size(); i++)
		dm[i] = vPByte.m_dim[i];

	for (i = vPByte.m_dim.size(); i < 4; i++)
		dm[i] = 1;

	if (dm[3] > MAXCLASSES) { ErrOutput("MarkrofField", u8"̫"); return; }

	//code 
	if (vG.m_dim.size() == 1)
	{
		code = 2;
		if (vG.m_dim[0] != dm[3])
			ErrOutput("MarkrofField", "Third arg has incompatible dimensions.");

		if (vG.m_type != DGDataType_Float) ErrOutput("MarkrofField", "Third arg must be single.");
	}
	else if (vG.m_dim.size() == 2)
	{
		code = 1;
		if (vG.m_dim[0] != dm[3] || vG.m_dim[1] != dm[3])
			ErrOutput("MarkrofField", "Third arg has incompatible dimensions.");

		if (vG.m_type != DGDataType_Float) ErrOutput("MarkrofField", "Third arg must be single.");
	}
	else if (vG.m_dim.size() == 5)
	{
		code = 3;
		for (i = 0; i < 4; i++)
			if (vG.m_dim[i] != dm[i])
				ErrOutput("MarkrofField", "Third arg has incompatible dimensions.");

		if (vG.m_dim[4] != dm[3])
			ErrOutput("MarkrofField", "Third arg has incompatible dimensions.");

		if (vG.m_type != DGDataType_Float) ErrOutput("MarkrofField", "Third arg must be single.");
	}

	else if (vG.m_dim.size() == 4)
	{
		for (i = 0; i < 3; i++)
			if (vG.m_dim[i] != dm[i])
				ErrOutput("MarkrofField", "Third arg has incompatible dimensions.");

		if (vG.m_dim[3] != (dm[3] * (dm[3] - 1)) / 2)
			ErrOutput("MarkrofField", "Third arg has incompatible dimensions.");

		if (vG.m_type == DGDataType_Float)
			code = 4;
		else if (vG.m_type != DGDataType_UnSignedChar)
			code = 5;
		else
			ErrOutput("MarkrofField", "Third arg must be either single or uint8.");
	}
	else
		ErrOutput("MarkrofField", "Third arg has incompatible dimensions.");


	for (i = 0; i < 4; i++)
		if (QFl.m_dim[i] != dm[i])
			ErrOutput("MarkrofField", u8"");

	int elcount = DimElCount(dm);

	if (vx2.size() > 0)
		for (i = 0; i < 3; i++) w[i] = vx2(i);
	else
		for (i = 0; i < 3; i++) w[i] = 1.0;

	mrf1((unsigned int *)&DimToSize(dm)[0], (unsigned char*)vPByte.Ptr<BYTE>(), QFl.Ptr<float>(), vG.Ptr<float>(), w, code);

}



void RGrid(DimVec d, CDGMulDimData & x)
{
	x = CDGMulDimData({ d[0],d[1],d[2],3 }, DGDataType_Float, DataMajor_FStyle_ColumnMajor);
	MAT_TYPE_F x1, x2;

	Grid2D<VECTOR_TYPE_F, MAT_TYPE_F>(BuildVector<VECTOR_TYPE_F, float>(1.0, d[0], 1.0),
		BuildVector<VECTOR_TYPE_F, float>(1.0, d[1], 1.0), x1, x2);



	for (int z = 0; z < d[2]; z++)
	{
		x.SelMatByDimXY_Maped<MAT_TYPE_F>({ z,0 }) = x1;
		x.SelMatByDimXY_Maped<MAT_TYPE_F>({ z,1 }) = x2;
		x.SelMatByDimXY_Maped<MAT_TYPE_F>({ z,2 }).array() = z + 1;//MAT_TYPE_F::Constant(d[0],d[1],z+1);
	}
}

void affind(CDGMulDimData & y0, CDGMulDimData & y1, MAT_TYPE & M)
{
	auto dm = y0.m_Acc.m_DataMajor;
	y1 = CDGMulDimData(y0.m_dim, y0.m_type, dm);
	DimVec ddim = { y0.m_dim[0],y0.m_dim[1],y0.m_dim[2] };
	for (int d = 0; d < 3; d++)
	{
		CDGMulDimData y1d(ddim, DGDataType_Float, dm);

		CDGMulDimData y0d1(ddim, DGDataType_Float, dm , y0.SelFrameByDim(3, { 0 }));
		y0d1.MulOn<float>(M(d, 0));
		CDGMulDimData y0d2(ddim, DGDataType_Float, dm , y0.SelFrameByDim(3, { 1 }));
		y0d2.MulOn<float>(M(d, 1));
		CDGMulDimData y0d3(ddim, DGDataType_Float, dm , y0.SelFrameByDim(3, { 2 }));
		y0d3.MulOn<float>(M(d, 2));

		y1d.AddOn<float>(y0d1); y1d.AddOn<float>(y0d2); y1d.AddOn<float>(y0d3);
		y1d.AddOn<float>(M(d, 3));

		memcpy(y1.SelFrameByDim(3, { d }), y1d.DataPtr(0), y1d.SizeInBytes());
	}
}



void GetClosestAffine(CDGMulDimData & x, CDGMulDimData & y, CDGMulDimData &w1, CDGMulDimData &w2, OUT MAT_TYPE &M, OUT MAT_TYPE *R)
{
	MAT_TYPE XX = MAT_TYPE::Zero(4, 4);
	MAT_TYPE XY = MAT_TYPE::Zero(4, 4);
	DimVec d = x.m_dim;
	VECTOR_TYPE_F o = VECTOR_TYPE_F::Ones(d[0] * d[1]);
	int framesize = d[0] * d[1];
	for (int z = 0; z < d[2]; z++)
	{
		MAT_TYPE_F xk = MAT_TYPE_F::Zero(framesize, 3);
		for (int i = 0; i < 3; i++)
		{
			xk.col(i) = x.SelMatByDimXY_Maped_AsVec<VECTOR_TYPE_F>({ z,i });			
		}
		VECTOR_TYPE_F ox(framesize), oy(framesize);
		if (w1.bEmpty() && w2.bEmpty())
		{
			ox = o; oy = o;
		}
		else if (!w1.bEmpty() && !w2.bEmpty())
		{
			oy << w2.SelMatByDimXY_Maped_AsVec<VECTOR_TYPE_F>({ z });
			//oy << tmp;
			ox << w1.SelMatByDimXY_Maped_AsVec<VECTOR_TYPE_F>({ z });
			//ox << tmp;
		}
		else if (!w1.bEmpty())
		{
			ox << w1.SelMatByDimXY_Maped_AsVec<VECTOR_TYPE_F>({ z });
			oy = ox;
		}
		else if (!w2.bEmpty())
		{
			ox << w2.SelMatByDimXY_Maped_AsVec<VECTOR_TYPE_F>({ z });
		}
		xk.col(0) = ((VECTOR_TYPE_F)xk.col(0)).array()* ox.array();
		xk.col(1) = ((VECTOR_TYPE_F)xk.col(1)).array()* ox.array();
		xk.col(2) = ((VECTOR_TYPE_F)xk.col(2)).array()* ox.array();

		MAT_TYPE_F yk = MAT_TYPE_F::Zero(framesize, 3);
		for (int i = 0; i < 3; i++)
		{
			yk.col(i) = y.SelMatByDimXY_Maped_AsVec<VECTOR_TYPE_F>({ z,i });
		}

		IntVEC allFiniteCols;
		for (int irow = 0; irow < yk.rows(); irow++)
		{
			if (isfinite(xk(irow, 0)) && isfinite(xk(irow, 1)) && isfinite(xk(irow, 2)) &&
				isfinite(yk(irow, 0)) && isfinite(yk(irow, 1)) && isfinite(yk(irow, 2)))
			{
				allFiniteCols.push_back(irow);
			}
		}
		MAT_TYPE_F X(allFiniteCols.size(), xk.cols() + ox.cols());
		MAT_TYPE_F Y(allFiniteCols.size(), yk.cols() + oy.cols());

		for (int ifi = 0; ifi < allFiniteCols.size(); ifi++)
		{
			int rowindex = allFiniteCols[ifi];
			X.row(ifi) << xk.row(rowindex), ox.row(rowindex);
			Y.row(ifi) << yk.row(rowindex), oy.row(rowindex);
		}
		MAT_TYPE DX = AsDouble(X);
		MAT_TYPE DY = AsDouble(Y);
		XX += DX.transpose()*DX;
		XY += DX.transpose()*DY;		
	}
	
	M = (XX.inverse()*XY).transpose();
	cout << M;
	if (R)
	{
		*R = MAT_TYPE::Zero(4, 4);
		
		MAT_TYPE XX1 = XX - (XX.col(3) * XX.col(3).transpose()) / XX(3, 3);

		MAT_TYPE XY1 = XY - (XY.col(3)*XY.row(3)) / XY(3, 3);
		MAT_TYPE Z = (XX1.block(0, 0, 3, 3).inverse()*XY1.block(0, 0, 3, 3)).transpose();
		auto svd = Z.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
		MAT_TYPE U = svd.matrixU();
		MAT_TYPE V = svd.matrixV();
		MAT_TYPE S = svd.singularValues();
		*R << U * V.transpose(), MAT_TYPE::Zero(3, 1), 0, 0, 0, 1;
		//MAT_TYPE id43 = MAT_TYPE::Identity(4, 3);
		MAT_TYPE T1(4, 4), T2(4, 4);
		T1.block(0, 0, 4, 3) = MAT_TYPE::Identity(4, 3);
		T1.block(0, 3, 4, 1) = -XY.col(3) / XY(3, 3);
		T2.block(0, 0, 4, 3) = MAT_TYPE::Identity(4, 3);
		T2.block(0, 3, 4, 1) = -XY.row(3).transpose() / XY(3, 3);
		*R = T2 * (*R)*T1;
	}
}

