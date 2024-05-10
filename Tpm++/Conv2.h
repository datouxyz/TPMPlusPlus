#pragma once

enum ConvMode
{	
	ConvMode_Full,	
	ConvMode_Same,	
	ConvMode_valid
};


template<class T>T Conv2(T & A, T &kernel, ConvMode mode)
{
	int KernelRow = kernel.rows();
	int KernelColumn = kernel.cols();
	int ARow = A.rows();
	int ACol = A.cols();

	int FullRow = KernelRow + ARow - 1;
	int FullCol = KernelColumn + ACol - 1;

	int CRowStart = 0, CRowEnd = 0;
	int CColStart = 0, CColEnd = 0;

	T res;

	switch (mode)
	{
	case ConvMode_Full:
	{
		CRowStart = 0;
		CRowEnd = FullRow;
		CColStart = 0;
		CColEnd = FullCol;
	}
	break;
	case ConvMode_Same:
	{
		CRowStart = round((FullRow - ARow) / 2.0);
		CRowEnd = CRowStart + ARow;
		CColStart = round((FullCol - ACol) / 2.0);
		CColEnd = CColStart + ACol;
	}
	break;
	case ConvMode_valid:
	{
		if (ARow < KernelRow || ACol < KernelColumn)
			return res;
		CRowStart = KernelRow - 1;
		CRowEnd = FullRow - (KernelRow - 1);
		CColStart = KernelColumn - 1;
		CColEnd = FullCol - (KernelColumn - 1);
	}
	break;
	}

	res = T::Zero(CRowEnd - CRowStart, CColEnd - CColStart);

	for (int i = CRowStart; i < CRowEnd; i++)
	{
		for (int j = CColStart; j < CColEnd; j++)
		{
			double temp = 0;
			for (int m = 0; m < KernelRow; m++)
			{
				for (int n = 0; n < KernelColumn; n++)
				{
					if ((i - m) >= 0 && (i - m) < ARow && (j - n) >= 0 && (j - n) < ACol)
					{
						temp += kernel(m, n) * A(i - m, j - n);// [i - m][j - n];
					}
				}
			}
			res(i - CRowStart, j - CColStart) = temp;
		}
	}
	return res;
}