//**************************************************
// file: astring.h
//**************************************************

#ifndef ASTRING_H__
#define ASTRING_H__

//--------------------------------------------------
bool CopyString( char** target, const char* source );
int CopyString( char** target, const char* start, const char* end );

int CatStrings( char** target, const char* source1, const char* source2 );

char* SkipWhite( char* ptr );
char* SkipNonWhite( char* ptr );
//--------------------------------------------------
#endif // ASTRING_H__

