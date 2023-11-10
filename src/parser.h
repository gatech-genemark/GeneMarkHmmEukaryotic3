//**************************************************
// file: parser.h
// project: ehmm3
//**************************************************

#ifndef PARSER_H__  //--
#define PARSER_H__  //--

#include "astring.h"
#include "fload.h"
#include "alist.h"

struct Record {
	char* key;
	char* start;
};

//--------------------------------------------------
class Parser : private Fload
{
public:

	~Parser() { recordList.Clean(); };

	bool InitParser( const char* name );

	List<Record> recordList;

	bool GetUniqValue( const char* key, char* value );
	bool GetUniqValue( const char* key, char** value );
	bool GetUniqValue( const char* key, int* value );
	bool GetUniqValue( const char* key, double* value );
	bool GetUniqValue( const char* key, bool* value );

	bool GetIntKeyArray( const char* key, double** data, int size );
	bool GetStrKeyArray( const char* key, double** data, int size, int width );
	bool GetStrKeyArray( const char* key, double** data[], int size, int width );

	bool IsKey( const char* key );

private:

	void Decode( char* ptr, int size );
	void RemoveComments( char* start);
	bool IndexRecords(char* str);
	Record* FindKey( const char* key);
	Record* FindUniqKey( const char* key );
	char* GetUniq( const char* key );
	char* GetUniq( const char* key, char sep );
	char* GetUniqWithMsg( const char* key );
	double* AllocArray( double* data, int size, const char* key );
	int StringToIndex( char* str, char* ptr );
	bool FillIntKeyArray( const char* key, double* data, int size );
	bool FillStrKeyArray( const char* key, double* data, int size, int width );
	bool FillStrKeyArray( const char* key, double** data[], int size, int width );
};
//--------------------------------------------------
#endif // PARSER_H__

