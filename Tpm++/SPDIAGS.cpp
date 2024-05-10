
#include "stdafx.h"
#include <Eigen/Dense>
//#include <unsupported/Eigen/MatrixFunctions>
#include "TissueProblityMap.h"
#include "EigenMatFunc.h"


MAT_TYPE diagk(const MAT_TYPE& X, int k)
{
	VECTOR_TYPE D;
	
	if (X.rows() != 1 && X.cols() != 1)
	{
		D = X.diagonal(k);
		//D = D.transpose();
	}
	else
	{
		VECTOR_TYPE DX = X.rows() > 1 ? (VECTOR_TYPE)X.col(0) : (VECTOR_TYPE)X.row(0);

		if (X.size()>0 && 0 <= k && 1 + k <= X.cols())
		{
			D = VECTOR_TYPE::Constant(1, DX(k + 1));
		}
		else if (X.size()>0 && k < 0 && 1 - k <= X.rows())
		{
			D = VECTOR_TYPE::Constant(1, DX(1 - k));
		}
		else
		{
			D = VECTOR_TYPE::Zero(1);
		}
	}

	return D;
}

bool Spdiags(MAT_TYPE & B, VECTOR_TYPE d,int nrow,int ncol, MAT_SPARSE_TYPE & ret)
{
	if (nrow > B.rows() || B.cols() != d.size())
	{
		return false;
	}
	ret = MAT_SPARSE_TYPE(nrow, ncol);
	TRIPLETS tr;

	//int NC = min((int)B.rows(),ncol);
	for (int irow = 0; irow < nrow; irow++)
	for (int icol = 0; icol < ncol; icol++)
	{
		for(int ind=0;ind<d.size();ind++)
		{
			int offset = d[ind];
			if(offset == icol- irow)
			{
				tr.push_back(TRIPLET_TYPE(irow,icol,B(irow,ind)));
			}			
		}
	}		
	ret.setFromTriplets(tr.begin(),tr.end());
	return true;
}
