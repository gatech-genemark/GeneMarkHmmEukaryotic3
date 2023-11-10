//**************************************************
// file: astring.cpp
//**************************************************

#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "astring.h"
#include "common.h"

//--------------------------------------------------
bool CopyString( char** target, const char* source )
{
	if ( !source ) return false;
	int size = strlen( source );
	if ( !size ) return false;

	if ( *target ) delete [] *target;
	*target = new char[ size + 1 ];
	if ( !*target ) return false;

	strcpy( *target, source );

	return true;
};
//--------------------------------------------------
// start <= x < end, including *start, excluding *end
//--------------------------------------------------
int CopyString( char** target, const char* start, const char* end )
{
	if ( !start || !end ) return 0;
	int size = end - start;
	if ( size <= 0 ) return 0;

	if ( *target ) delete [] *target;
	*target = new char[ size + 1 ];
	if ( !*target ) return 0;

	strncpy( *target, start, size );
	*(*target + size) = '\0';

	return size;
};
//--------------------------------------------------
int CatStrings( char** target, const char* source1, const char* source2 )
{
	if ( !source1 || !source2 ) return 0;
	int size = strlen( source1 ) + strlen( source2 );
	if ( !size ) return 0;

	if ( *target ) delete [] *target;
	*target = new char [ size + 1 ];
	if ( !*target ) return 0;

	strcpy( *target, source1 );
	strcat( *target, source2 );

	return  size;
};
//--------------------------------------------------
char* SkipWhite( char* ptr )
{
	if (!ptr) return NULL;

	while ( isspace(*ptr) )
		++ptr;

	if ( *ptr == '\0' ) 
		ptr = NULL;

	return ptr;
};
//--------------------------------------------------
char* SkipNonWhite( char* ptr )
{
	if (!ptr) return NULL;

	if ( isspace(*ptr) ) return NULL;

	while ( !isspace(*ptr) ) {

		if ( *ptr == '\0' ) 
			return NULL;
		++ptr;
	}

	return ptr;
};
//--------------------------------------------------

