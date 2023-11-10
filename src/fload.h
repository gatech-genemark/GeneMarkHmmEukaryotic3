//**************************************************
// file: fload.h
// project: ehmm3
// Alex Lomsadze
//**************************************************

#ifndef FLOAD_H__
#define FLOAD_H__

//--------------------------------------------------
class Fload
{
public:

	Fload() { fileBuffer = NULL; fileSize = 0; };
	~Fload() { ClearFileBuffer(); };

	char* fileBuffer;
	int fileSize;

	bool LoadFile( const char* name );
	void ClearFileBuffer() { delete [] fileBuffer; fileBuffer = NULL; fileSize = 0; };

private:

	int GetFileSize( const char* name );
	bool PutFileInBuffer( const char* name );
};
//--------------------------------------------------
#endif  // FLOAD_H__

