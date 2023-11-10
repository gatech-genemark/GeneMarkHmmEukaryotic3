//**************************************************
// file: fload.cpp
// project: ehmm3
// Alex Lomsadze
//**************************************************

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "fload.h"
#include "common.h"

//--------------------------------------------------
int Fload::GetFileSize( const char* name )
{
	int size = -1;
	// required header <sys/types.h> <sys/stat.h>
	struct stat st;

	// get status information on a file size in bytes
	if ( !stat( name, &st ) )
		size = st.st_size;

	return size;
};
//--------------------------------------------------
bool Fload::PutFileInBuffer( const char* name )
{
	FILE* fp = fopen( name, "rb" );
	if ( !fp ) return false;

	int i = fread( fileBuffer, sizeof(char), fileSize, fp);

	if ( i != fileSize )
		return false;

	if ( ferror(fp) ) return false;
	fclose( fp );

	return true;
};
//--------------------------------------------------
bool Fload::LoadFile( const char* name )
{
	fileSize = GetFileSize( name );
	
	if ( fileSize < 0 )
	{
		printf( "%s: error, file <%s> is not found\n", PROG_NAME, name );
		return false;
	}
	if ( fileBuffer ) delete [] fileBuffer;

	fileBuffer = new char[ fileSize + 1 ];

	if ( !fileBuffer )
	{
		printf( "%s: error, not enough memory to load file <%s>\n", PROG_NAME, name );
		return false;
	}
	if ( !PutFileInBuffer( name ) )
	{
		printf( "%s: error, reading file <%s>\n", PROG_NAME, name );
		return false;
	}

	fileBuffer[ fileSize ] = '\0';

	return true;
};
//--------------------------------------------------

