//**************************************************
// file: matrix.cpp
// project: ehmm3
//**************************************************

#include <stdio.h>
#include <memory.h>
#include <math.h>

#include "matrix.h"
#include "common.h"

//================================================
// size = 4 ^ ( model order + 1 )
//
// order   0   1   2    3     4     5
//  size   4  16  64  256  1024  4096
//================================================

//--------------------------------------------------
// P(dir)(  | e1 | e2 | e3 | ...           | en |  )
// P(rev)(  | e1 | e2 | e3 | ...           | en |  )
//
// P(dir)( N ) != P(rev)( N ) == P(dir)( N' )
// calc N' from N 
//
//          | en |           ... | e3 | e2 | e1 |
//             |                    |    |    |
// P(dir)(  | en'|           ... | e3'| e2'| e1'|  )
//
// P(rev) == P(dir)
//
// to convert element in complement 
// element - { A C G T }
// switch( element ) {         | if elements order are fixed
//   case A: element = T;      | A C G T
//   case C: element = G;      | than 3 - element ==
//   case G: element = C;      | T G C A
//   case T: element = A;      | 
// }                           | 
//--------------------------------------------------
void Matrix::DoRevFromDir( double* dir, double* rev , int order)
{
	int n = 0;
	int calc;
	int temp;
	const int mask = 3;				// mask == 0000 0011

	int size = 4<<(2*order);

	while ( n < size ) {

		temp = n;
		calc = 0;

		for ( int count = 0; count <= order; ++count ) {

			calc <<= 2;
			calc +=  3 - ( temp & mask);
			temp >>= 2;
		}

		rev[n] = dir[calc];
		++n;
	}
};
//--------------------------------------------------
//  p(AT)
//  p(AG)
//  p(AC)
//  p(AA)
//  sum = p(AA) + p(AC) + p(AG) + p(AT) = P(A*)
//  P(T|A) = p(AT) / p(A*)
//--------------------------------------------------
void Matrix::AbsoluteToRelativeRight( double* ref, int order )
{
	int size = 4<<(order*2);

	for( int i = 0; i < size; i += 4 )
	{
		double sum = ref[i] + ref[i+1] + ref[i+2] + ref[i+3];

		ref[i]   /= sum;
		ref[i+1] /= sum;
		ref[i+2] /= sum;
		ref[i+3] /= sum;
	}
};
//--------------------------------------------------
// sum of   hiAbs( - - ... - T )
//          hiAbs( - - ... - G )
//          hiAbs( - - ... - C )
//          hiAbs( - - ... - A )
//          ------------------
//   Prob   loAbs( - - ... - )
//--------------------------------------------------
void Matrix::DoHiToLo( double* hi, double* lo , int hiOrder )
{
	double* loopEnd = hi + ( 4 << ( 2*hiOrder ) );

	while ( hi < loopEnd ) {

		*lo = hi[A] + hi[C] + hi[G] + hi[T];
		++lo;
		hi += 4;
	}
};
//--------------------------------------------------
//  P( e1 before e2 e3 ... en ) = P( e1 e2 e3 ... en ) / P( e2 e3 ... en )
//
//	P( e2 e3 ... en ) = sum of
//         P( T e2 e3 ... en )
//         ...
//         P( G e2 e3 ... en )
//         ...
//         P( C e2 e3 ... en )
//         ...                     1/4 of matrix
//         P( A e2 e3 ... en )
//--------------------------------------------------
void Matrix::DoAbsToFirst( double* abs, double* rel, int order )
{
	int size = 4<<(2*order);
	double* loopEnd = abs + size/4;
	double temp;

	const int leftA = ( A << (2*order) );
	const int leftC = ( C << (2*order) );
	const int leftG = ( G << (2*order) );
	const int leftT = ( T << (2*order) );

	while ( abs < loopEnd ) {

		temp = abs[ leftA ] + abs[ leftC ] + abs[ leftG ] + abs[ leftT ];

		if ( temp != 0 ) {
			rel[ leftA ] = abs[ leftA ] / temp;
			rel[ leftC ] = abs[ leftC ] / temp;
			rel[ leftG ] = abs[ leftG ] / temp;
			rel[ leftT ] = abs[ leftT ] / temp;
		} else {
			rel[ leftA ] = 0;
			rel[ leftC ] = 0;
			rel[ leftG ] = 0;
			rel[ leftT ] = 0;
		}

		++abs;
		++rel;
	}
};
//--------------------------------------------------
bool Matrix::AllocArray( double** data, int size )
{
	if ( *data ) delete [] *data;
	*data = new double [ size ];
	if ( !*data ) {
		printf( "error, out of memory\n" );
		return false;
	}
	memset( *data, 0, size*sizeof(double) );

	return true;
};
//--------------------------------------------------
bool Matrix::RevFromDir( double* dir, double** rev , int order )
{
	if ( !AllocArray( rev, (4<<(2*order)) + 1 ) )
		return false;

	DoRevFromDir( dir, *rev , order );

	(*rev)[4<<(2*order)] = dir[4<<(2*order)];

	return true;
};
//--------------------------------------------------
bool Matrix::AbsToFirst( double* abs, double** rel, int order )
{
	if ( !AllocArray( rel, (4<<(2*order)) + 1 ) )
		return false;

	DoAbsToFirst( abs, *rel , order );

	(*rel)[4<<(2*order)] = abs[4<<(2*order)];

	return true;
};
//--------------------------------------------------
void Matrix::Log( double* ptr , int size )
{
	for ( int i = 0; i < size; ++i )
	{
		if ( ptr[i] < 0 )
			i=i;

		if ( ptr[i] )
			ptr[i] = log( ptr[i] );
		else
			ptr[i] = log0;
	}
};
//--------------------------------------------------
bool Matrix::HiToLo( double* hi, double** lo, int hiOrder )
{
	if ( hiOrder < 1 )
	{
		printf( "error in calculation of low order matrix\n" );
		return false;
	}
	if ( !AllocArray( lo, (4<<(2*(hiOrder-1))) + 1 ) )
		return false;

	DoHiToLo( hi, *lo, hiOrder );

	(*lo)[4<<(2*(hiOrder - 1))] = hi[4<<(2*hiOrder)];

	return true;
};
//--------------------------------------------------

