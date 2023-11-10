//**************************************************
// file: site.h
// project: ehmm3
//**************************************************

#ifndef SITE_H__  //--
#define SITE_H__  //--

#include <stdio.h>  //--

//--------------------------------------------------
class Site
{
public:

	Site() { width = 0; margin = 0; order = 0; matDir = NULL; matRev = NULL; p_margin = 0; p_width = 0; };  //--
	~Site() { delete [] matDir; delete [] matRev; };  //--

	bool InitSite( const char* msg );

	bool Adjust( int* nt_count, int site_size );

	double GetValueDir( char* str, int shift );
	double GetValueRev( char* str, int shift );

	double GetValueDirNew( char* str, int shift );
	double GetValueRevNew( char* str, int shift );

	void NormSite( double** non, double** cod_1, double** cod_2, double** cod_3, int phase, bool cod_left, int site_size, int ch_order );

	int margin;  //--
	int width;  //--
	int order;  //--

	double* matDir;  //--
	double* matRev;  //--

	int p_margin;
	int p_width;

	void Print( char* name );

private:

	double nt[4];
	int siteSize;

	int size;  //--

	bool AllocRev( const char* msg );  //--
	void Log( double* ptr, int size );  //--
	void DirToRev();  //--
	int RevCompKey( int key );

	void SetNonCoding( int width, double* ptr, int max_order, double** non );
	void SetCoding( int width, double* ptr, int max_order, double** cod_1, double** cod_2, double** cod_3, int frame );
};
//--------------------------------------------------
#endif  // SITE_H__  //--

