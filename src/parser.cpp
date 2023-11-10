//**************************************************
// file: parser.cpp
// project: ehmm3
//**************************************************

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include "parser.h"
#include "common.h"
#include "astring.h"

//--------------------------------------------------
bool Parser::InitParser( const char* name )
{
	if ( !LoadFile( name ) ) return false;

	Decode( fileBuffer, fileSize );

	RemoveComments( fileBuffer );

	if ( !IndexRecords( fileBuffer ) ) {
		printf( "error in file format, %s\n", name );
		return false;
	}
	return true;
};
//--------------------------------------------------
void Parser::Decode( char* ptr, int size )
{
	if ( *ptr == '?') 
	{
		*ptr = ' ';
		*(ptr+1) = ' ';

		for ( int i = 2; i < size; ++i )
			*(ptr + i) += '0';

		*(ptr + size) = '\0';
	}
};
//--------------------------------------------------
void Parser::RemoveComments( char* start )
{
	char* end = NULL;

	while ( start = strchr( start, '#' ) ) {

		end = strpbrk( start, "\n\r" );
		if ( !end )
			end = start + strlen( start );
		memset( start, ' ', end - start );
	}
};
//--------------------------------------------------
bool Parser::IndexRecords( char* str )
{
	Record* data;

	while ( isspace(*str) ) ++str;
	if ( *str != '$' ) return false;

	data = new Record;
	if ( !data )
	{
		printf("%s: error, not enough memory\n", PROG_NAME);
		return false;
	}
	data->key = str;
	recordList.AddToTail( data );

	while ( str = strchr( str + 1, '$' ) ) {

		if ( !isspace( *(str - 1) ) ) return false;
		*(str - 1) = '\0';

		data = new Record;
		if (!data)
		{
			printf("%s: error, not enough memory\n", PROG_NAME);
			return false;
		}
		data->key = str;
		recordList.AddToTail( data );
	}

	char* tmp;

	recordList.SetCurrentToFirst();
	while ( !recordList.IsEnd() ) {
		data = recordList.GetForward();

		data->key++;
		tmp = SkipNonWhite( data->key );
		if ( !tmp ) return false;
		*tmp = '\0';

		data->start = SkipWhite( tmp + 1);
		if ( !data->start ) return false;
	}
	return true;
};
//--------------------------------------------------
Record* Parser::FindKey( const char* key )
{
	Record* rec;

	while( !recordList.IsEnd() ) {

		rec = recordList.GetForward();
		if ( !strcmp( rec->key, key ) )
			return rec;
	}
	return NULL;
};
//--------------------------------------------------
Record* Parser::FindUniqKey( const char* key )
{
	Record* rec;

	recordList.SetCurrentToFirst();
	if ( !(rec = FindKey(key)) ) return NULL;
	if ( FindKey(key) ) return NULL;

	return rec;
};
//--------------------------------------------------
bool Parser::IsKey( const char* key )
{
	Record* rec;
	recordList.SetCurrentToFirst();
	while( !recordList.IsEnd() )
	{
		rec = recordList.GetForward();
		if ( !strcmp( rec->key, key ) )
			return true;
	}
	return false;
};
//--------------------------------------------------
char* Parser::GetUniq( const char* key )
{
	Record* rec = FindUniqKey( key );
	if ( !rec ) return NULL;

	char* tmp = SkipNonWhite( rec->start );
	if ( SkipWhite( tmp ) ) return NULL;

	return rec->start;
};
//--------------------------------------------------
char* Parser::GetUniqWithMsg( const char* key )
{
	char* str = GetUniq( key );
	if ( !str ) {
		printf( "error reading parameter %s\n", key );
		return NULL;
	}
	return str;
}
//--------------------------------------------------
char* Parser::GetUniq( const char* key, char sep )
{
	Record* rec = FindUniqKey( key );
	if ( !rec ) return NULL;

	if ( *rec->start != sep ) return NULL;

	char* tmp = strchr( rec->start + 1, sep );
	if (!tmp) return NULL;

	if ( SkipWhite( ++tmp ) ) return NULL;

	*tmp= '\0';

	return rec->start;
};
//--------------------------------------------------
bool Parser::GetUniqValue( const char* key, char* value )
{
	char* str = GetUniqWithMsg( key );
	if ( !str ) return false;
		
	if ( strlen(str) != 1 ) return false;

	*value = *str;
	return true;
};
//--------------------------------------------------
bool Parser::GetUniqValue( const char* key, int* value )
{
	char* str = GetUniqWithMsg( key );
	if ( !str ) return false;

	*value = atoi(str);
	return true;
};
//--------------------------------------------------
bool Parser::GetUniqValue( const char* key, double* value )
{
	char* str = GetUniqWithMsg( key );
	if ( !str ) return false;

	*value = atof(str);
	return true;
};
//--------------------------------------------------
bool Parser::GetUniqValue( const char* key, bool* value )
{
	char* str = GetUniqWithMsg( key );
	if ( !str ) return false;

	double x = atof(str);
	if ( x == 0 )
		*value = false;
	else
		*value = true;

	return true;
};
//--------------------------------------------------
bool Parser::GetUniqValue( const char* key, char** value )
{
	char* str = GetUniq( key, '"' );
	if (!str) {
		printf( "error reading parameter %s\n", key );
		return false;
	}
	if ( !CopyString( value, str + 1, str + strlen(str) - 1 ) )
		return false;

	return true;
};
//--------------------------------------------------
// "size" memory must be allocated before this function
//--------------------------------------------------
bool Parser::FillIntKeyArray( const char* key, double* data, int size )
{
	Record* rec = FindUniqKey( key );
	if (!rec) return false;

	char* ptr;
	int index;
	double value;

	for ( char* str = rec->start; str; ptr = SkipNonWhite( str ), str = SkipWhite(ptr) ) {

		if ( !(ptr = SkipNonWhite( str )) ) return false;

		index = atoi( str );

		if ( !(str = SkipWhite( ptr )) ) return false;
		
		value = atof( str );

		if ( index < 1 ) 
		{
			printf( "gmhmme error: index %d is out of range (less then 1) for key %s\n", index, key );
			return false;
		}

		if ( index <= size )
			data[index - 1] = value;
		else
			printf( "gmhmme warning: index %d is out of range for max %d for key %s\n", index, size, key );
	}
	return true;
};
//--------------------------------------------------
bool Parser::GetIntKeyArray( const char* key, double** data, int size )
{
	*data = AllocArray( *data, size, key );
	if ( !*data ) return false;

	if ( !FillIntKeyArray( key, *data, size ) ) {
		printf( "error reading parameter %s\n", key );
		return false;
	}
	return true;
};
//--------------------------------------------------
int Parser::StringToIndex( char* str, char* ptr )
{
	int index = 0;
	int max = ptr - str;

	for( int i = 0; i < max; ++i ) {

		index <<= 2;

		switch ( str[i] ) {
			case 'A': index += A; break;
			case 'C': index += C; break;
			case 'G': index += G; break;
			case 'T': index += T; break;
		default:
			return -1;
		}
	}
	return index;
};
//--------------------------------------------------
bool Parser::FillStrKeyArray( const char* key, double* data, int size, int width )
{
	Record* rec = FindUniqKey( key );
	if (!rec) return false;

	int i;
	int index;
	char* ptr;

	for ( char* str = rec->start ; str; str = SkipWhite( ptr) ) {

		if ( !(ptr = SkipNonWhite( str )) ) return false;

		index = StringToIndex( str, ptr );
		
		if (( index >= size )||( index < 0 )) return false;

		for ( i = 0; i < width; ++i ) {

			if ( !(str = SkipWhite( ptr )) ) return false;

			data[ i*size + index ] = atof( str );

			ptr = SkipNonWhite(str);
		}
	}
	return true;
};
//--------------------------------------------------
bool Parser::GetStrKeyArray( const char* key, double** data, int size, int width )
{
	*data = AllocArray( *data, size*width + 1, key );
	if ( !*data ) return false;

	if ( !FillStrKeyArray( key, *data, size, width ) ) {
		printf( "error reading parameter %s\n", key );
		return false;
	}
	return true;
};
//--------------------------------------------------
bool Parser::FillStrKeyArray( const char* key, double** data[], int size, int width )
{
	Record* rec = FindUniqKey( key );
	if (!rec) return false;

	int i;
	int index;
	char* ptr;

	for ( char* str = rec->start; str; str = SkipWhite( ptr) ) {

		if ( !(ptr = SkipNonWhite( str )) ) return false;

		index = StringToIndex( str, ptr );

		if (( index >= size )||( index < 0 )) return false;

		for ( i = 0; i < width; ++i ) {

			if ( !(str = SkipWhite( ptr )) ) return false;

			*(*(data[i]) + index ) = atof( str );

			ptr = SkipNonWhite(str);
		}
	}
	return true;
};
//--------------------------------------------------
bool Parser::GetStrKeyArray( const char* key, double** data[], int size, int width )
{
	for ( int i = 0; i < width; ++i ) {

		*(data[i]) = AllocArray( *(data[i]), size + 1, key );
		if ( !*(data[i]) ) return false;
	}

	if ( !FillStrKeyArray( key, data, size, width ) ) {
		printf( "error reading parameter %s\n", key );
		return false;
	}
	return true;
};
//--------------------------------------------------
double* Parser::AllocArray( double* data, int size, const char* key )
{
	if ( data ) delete [] data;
	data = new double [ size ];
	if ( !data ) {
		printf( "error, out of memory loading %s\n", key );
		return NULL;
	}

	memset( data, 0, size*sizeof(double) );

	return data;
};
//--------------------------------------------------

