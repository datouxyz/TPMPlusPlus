#pragma once


template<typename TM, typename TV>
bool MatVecBsxFunction(BsxOperator op, bool Ascol/*as col wrap*/,IN const TM & A, IN const TV & B, TM& result)
{
	CTimeCounter cc("MatVecBsxFunction");
	if (Ascol && (A.cols() != B.size()))
	{
		ErrOutput("bsxfunc", " Cols != size");
		return false;
	}
	if (!Ascol && (A.rows() != B.size()))
	{
		ErrOutput("bsxfunc", " rows != size");
		return false;
	}

	
	result = TM::Zero((int)A.rows(), (int)A.cols());
	
	DGDataType dt = EiGenElDGType<TM>();
	BsxfFuncPointer bsxffunc = GetBxsfFunc(op, dt);

	{
		int ncol = A.cols();
		int nrow = A.rows();//lhs.size() / batchsize;

		TM::Scalar* lhsPtr = (TM::Scalar*)A.data();
		TM::Scalar* resPtr = (TM::Scalar*)result.data();

		TM::Scalar* BPtr = (TM::Scalar*)B.data();

		if (Ascol)
		for (int icol = 0; icol < ncol; icol++)
		{
			TM::Scalar* APtr = &lhsPtr[icol * nrow];
			TM::Scalar* ResPtr = &resPtr[icol * nrow];
			TM::Scalar Bval = BPtr[icol];
			for (int irow = 0; irow < nrow; irow++)
			{
				bsxffunc((byte*)&APtr[irow], (byte*)&Bval, (byte*)&ResPtr[irow]);
			}
		}
		else
		for (int icol = 0; icol < ncol; icol++)
		{
			TM::Scalar* APtr = &lhsPtr[icol * nrow];
			TM::Scalar* ResPtr = &resPtr[icol * nrow];
						
			for (int irow = 0; irow < nrow; irow++)
			{
				bsxffunc((byte*)&APtr[irow], (byte*)&BPtr[irow], (byte*)&ResPtr[irow]);
			}			
		}		
	}

	return true;
}


template<typename TM, typename TV>
bool MatVecBsxFunction2(BsxOperator op, bool Ascol/*as col wrap*/, IN const TM & A, IN const TV & B, TM& result)
{
	CTimeCounter cc("MatVecBsxFunction");
	if (Ascol && (A.cols() != B.size()))
	{
		ErrOutput("bsxfunc", " Cols != size");
		return false;
	}
	if (!Ascol && (A.rows() != B.size()))
	{
		ErrOutput("bsxfunc", " rows != size");
		return false;
	}


	result = TM::Zero((int)A.rows(), (int)A.cols());

	switch (op)
	{
		case BsxOperator_times:
		{
			if (Ascol)
				for (int irow = 0; irow < A.rows(); irow++)
					result.row(irow) = ((TV)A.row(irow)).array() * B.array();
			else
				for (int icol = 0; icol < A.cols(); icol++)
					result.col(icol) = ((TV)A.col(icol)).array() * B.array();
		}break;
		case BsxOperator_divide:
		{
			if (Ascol)
				for (int irow = 0; irow < A.rows(); irow++)
					result.row(irow) = ((TV)A.row(irow)).array() / B.array();
			else
				for (int icol = 0; icol < A.cols(); icol++)
					result.col(icol) = ((TV)A.col(icol)).array() / B.array();
		}break;
		case BsxOperator_minus:
		{
			if (Ascol)
				for (int irow = 0; irow < A.rows(); irow++)
					result.row(irow) = ((TV)A.row(irow)).array() - B.array();
			else
				for (int icol = 0; icol < A.cols(); icol++)
					result.col(icol) = ((TV)A.col(icol)).array() - B.array();
		}break;
		case BsxOperator_plus:
		{
			if (Ascol)
				for (int irow = 0; irow < A.rows(); irow++)
					result.row(irow) = ((TV)A.row(irow)).array() + B.array();
			else
				for (int icol = 0; icol < A.cols(); icol++)
					result.col(icol) = ((TV)A.col(icol)).array() + B.array();
		}break;
	}
	
	return true;
}



template<typename TM, typename TV>
TM log_spatial_priors(IN const TM &B, IN const TV& wp)
{
	TM bwp;
	MatVecBsxFunction<TM,TV>(BsxOperator_times, true, B, wp, bwp);	
	TV BS = 1.0 / bwp.rowwise().sum().array();	
	TM ret;
	MatVecBsxFunction<TM,TV>(BsxOperator_times, false, bwp, BS, ret);
	ret = ret.array().log();	
	return ret;
}

template<typename TM,typename TV >
double safe_softmax(IN OUT TM &Q)
{
	TV vmaxQ = Q.rowwise().maxCoeff();
	TM MinusMaxQ;	
	MatVecBsxFunction<TM,TV>(BsxOperator_minus, false, Q, vmaxQ, MinusMaxQ);
	TM minQexp = MinusMaxQ.array().exp();
	TV sQ = minQexp.rowwise().sum();	
	double ll = (sQ.array().log() + vmaxQ.array()).sum();
	MatVecBsxFunction<TM, TV>(BsxOperator_divide,false, minQexp, sQ, Q);
	return ll;
}
