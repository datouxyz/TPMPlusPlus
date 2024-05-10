#include "stdafx.h"

#define CpuNum "wmic cpu get DeviceID"

#define CpuCoreNum "wmic cpu get NumberOfCores"

#define CpuLogicalCoreNum "wmic cpu get NumberOfLogicalProcessors"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
void getCpuInformation(const char *command, ByteVecs & rets)
{
	FILE *winCommand = _popen(command, "r");

	if (!winCommand)
	{
		perror("popen");
		exit(EXIT_FAILURE);
	}
	else
	{
		ByteVec b; b.resize(100);
		
		while (fgets((char*)&b[0], b.size() - 1, winCommand) != 0)
		{
			//printf("%s", buf);
			rets.push_back(b);

		}

		_pclose(winCommand);
	}
}
int GetNumPhyCpu()
{
	ByteVecs ret;
	getCpuInformation(CpuCoreNum, ret);
	return atoi((char*)&ret[1][0]);
	/*getCpuInformation(CpuNum);
	getCpuInformation(CpuCoreNum);
	getCpuInformation(CpuLogicalCoreNum);*/
}