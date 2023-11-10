//**************************************************
// file: site.cpp
// project: ehmm3
//**************************************************

#include <math.h>
#include <cstdlib>
#include <memory.h>

#include "site.h"
#include "common.h"

#include <string>
#include <fstream>
#include <iostream>

//--------------------------------------------------
bool Site::AllocRev( const char* msg )
{
	if ( matRev ) delete [] matRev;
	matRev = new double [ size * width ];
	if ( !matRev ) {
		printf( "error, out of memory %s\n", msg ); 
		return false;
	}
	return true;
};
//--------------------------------------------------
void Site::Log( double* ptr, int size )
{
	for ( int i = 0; i < size; ++i ) {

		if ( ptr[i] )
			ptr[i] = log( ptr[i] );
		else
			ptr[i] = log0;
	}
};
//--------------------------------------------------
void Site::DirToRev()
{
	int i;
	int j;

	for ( i = 0; i < width; ++i ) {
		for ( j = 0; j < size; ++j ) {

			matRev[ size*(width - i - 1) + RevCompKey(j)] = matDir[size*i + j];
		}
	}
};
//--------------------------------------------------
int Site::RevCompKey( int key )
{
	int k = 0;

	for ( int i = 0; i <= order; ++i ) {

		k <<= 2;
		k += (3 - key&3);
		key >>= 2;
	}

	return k;
};
//--------------------------------------------------
bool Site::InitSite( const char* msg )
{
	size = 4<<(2*order);  //--

	if ( !AllocRev( msg ) ) return false;  //--

	Log( matDir, size*width );  //--
	DirToRev();  //--

	p_width = width;
	p_margin = margin;

	return true;  //--
};
//--------------------------------------------------
bool Site::Adjust( int* nt_count, int site_size )
{
	int all = nt_count[0] + nt_count[1] + nt_count[2] + nt_count[3];
	if (!all) return false;

	nt[0] = (double)nt_count[0]/all;
	nt[1] = (double)nt_count[1]/all;
	nt[2] = (double)nt_count[2]/all;
	nt[3] = (double)nt_count[3]/all;

	Log( nt, 4 );

	siteSize = site_size;

	return true;
};
//--------------------------------------------------
double Site::GetValueDir( char* str, int shift )
{
	str -= ( margin + shift );
	double x = 0;
	int i;

	for ( i = 0; i < width; ++i ) {
		x += matDir[ size*i + str[i] ];
		x -= nt[ str[i] ];
	}

	for ( i = margin; i < margin + siteSize; ++i )
		x += nt[ str[i] ];

	return x;
};
//--------------------------------------------------
double Site::GetValueRev( char* str, int shift )
{
	str -= ( width - margin - shift );
	double x = 0;
	int i;

	for ( i = 0; i < width; ++i ) {
		x += matRev[ size*i + str[i] ];
		x -= nt[ 3 - str[i] ];
	}

	for ( i = width - margin - siteSize; i < width - margin; ++i )
		x += nt[ 3 - str[i] ];

	return x;
};
//--------------------------------------------------

double Site::GetValueDirNew( char* str, int shift )
{
	str -= ( margin + shift );
	double x = 0;
	int i;
	
	int index = 0;
	int mask = (4<<(2*order)) - 1;

	for ( i = 0; i < width; ++i )
	{
		index <<= 2;
		if ( str[i] != N )
			index += str[i];
		else
			index += 0;   // temp, change later
		index &= mask;

		x += matDir[ size*i + index ];
	}

	return x;
};
//--------------------------------------------------
double Site::GetValueRevNew( char* str, int shift )
{
	str -= ( width - margin - shift );
	double x = 0;
	int i;
	
	int index = 0;
	int mask = (4<<(2*order)) - 1;

	for ( i = width - 1; i >=0; --i )
	{
		index <<= 2;
		if ( str[i] != N )
			index += (3-str[i]);
		else
			index += 0;   // temp, change later
		index &= mask;

		x += matDir[ size*(width - i - 1) + index ];
	}

	return x;
};
//--------------------------------------------------
void Site::NormSite( double** non, double** cod_1, double** cod_2, double** cod_3, int phase, bool cod_left, int site_size, int ch_order )
{
	double* data =  new double [ size * width ];
	if ( !data ) { printf( "error, out of memory\n"); exit(1); }
	memset( data, 0, (size * width)*sizeof(double) );

	// 1 2 3 4 5 6 7 8 9 
	// c c c G T n n n n

	//        left part | site | right part
	//  width     x     |  y   |   z

	int x = margin;
	int y = site_size;
	int z = width - margin - site_size;

	// maximum order for overlapping parts
	int max_order;

	if( cod_left )
	{
		if ( x > 0 )
		{
			max_order = order;
			if( x - 1 < max_order ) max_order = x - 1;
			if( ch_order < max_order ) max_order = ch_order;

			SetCoding( x, data, max_order, cod_1, cod_2, cod_3, ( 3 - margin%3 + phase )%3 + 1 );
		}

		if ( z > 0 )
		{
			max_order = order;
			if( z - 1 < max_order ) max_order = z - 1;
			if( ch_order < max_order ) max_order = ch_order;

			SetNonCoding( z, data + (x+y)*size, max_order, non );
		}
	}
	else
	{
		if ( x > 0 )
		{
			max_order = order;
			if( x - 1 < max_order ) max_order = x - 1;
			if( ch_order < max_order ) max_order = ch_order;

			SetNonCoding( x, data, max_order, non );
		}

		if ( z > 0 )
		{
			max_order = order;
			if( z - 1 < max_order ) max_order = z - 1;
			if( ch_order < max_order ) max_order = ch_order;

			SetCoding( z, data + (x+y)*size, max_order, cod_1, cod_2, cod_3, phase + 1 );
		}
	}

	// set y
	for( int i = x; i < x+y; ++i )
		for( int j = 0; j < size; ++j )
			data[ i*size +j] = 0;

	// add to direct values
	for( int i = 0; i < size*width; ++i )
	{
		if( matDir[i] != log0 )
			matDir[i] -= data[i];
	}

	delete [] data;

	DirToRev();
};
//--------------------------------------------------
void Site::SetNonCoding( int width, double* ptr, int max_order, double** non )
{
	// ptr - put result here

	int index = 0;
	int current_order = 0;

	for( int current_pos = 0; current_pos < width; ++current_pos )
	{
		if ( current_pos > max_order )
			current_order = max_order;
		else
			current_order = current_pos;

		int k = 0;
		int max_k = 4 << (current_order*2);

		for( int i = 0; i < size; ++i )
		{
			ptr[index] = non[current_order][k];
			++index;
			++k;
			if ( k >= max_k ) k = 0;
		}
	}
};
//--------------------------------------------------
void Site::SetCoding( int width, double* ptr, int max_order, double** cod_1, double** cod_2, double** cod_3, int frame )
{
	// ptr - put result here
	
	int index = 0;
	int current_order = 0;

	double** current_mat = 0;
	double** next_mat = 0;

	if( frame == 1 )
		current_mat = cod_1;
	else if( frame == 2 )
		current_mat = cod_2;
	else if( frame == 3 )
		current_mat = cod_3;

	for( int current_pos = 0; current_pos < width; ++current_pos )
	{
		if ( current_pos > max_order )
		{
			current_order = max_order;

			if( current_mat == cod_1 )
				current_mat = cod_2;
			else if( current_mat == cod_2 )
				current_mat = cod_3;
			else if( current_mat == cod_3 )
				current_mat = cod_1;
		}
		else
			current_order = current_pos;

		int k = 0;
		int max_k = 4 << (current_order*2);

		for( int i = 0; i < size; ++i )
		{
			ptr[index] = current_mat[current_order][k];
			++index;
			++k;
			if ( k >= max_k ) k = 0;
		}
	}
};
//--------------------------------------------------
void Site::Print( char* name )
{
	double x;

	std::ofstream log_file( name );
	if ( !log_file.is_open() ) { std::cout << PROG_NAME << " : error, opening file " << name << "\n";  exit(1); }

	for( int i = 0; i < size; ++i )
	{
		log_file << i << "\t";

		for( int j = 0; j < width; ++j )
		{
			x = matDir[ size*j + i ];
			log_file << x << "\t";
		}

			log_file << "\n";
	}

	log_file.close();
};
//--------------------------------------------------

