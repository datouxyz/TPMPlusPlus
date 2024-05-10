#include "stdafx.h"
string LoadStringFromFile(string filepath)
{
	FILE* f = fopen(filepath.c_str(), "r");		
	if (!f) return false;
	fseek(f, 0, SEEK_END);
	fpos_t filesize = 0;
	fgetpos(f, &filesize);	
	fseek(f, 0, SEEK_SET);
	ByteVec buffer; buffer.resize(filesize);
	fread(&buffer[0], 1, filesize, f);
	fclose(f);

	return (char*)&buffer[0];	
}
