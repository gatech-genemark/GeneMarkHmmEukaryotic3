//**************************************************
// file: markov.cpp
// project: ehmm3
//**************************************************

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <cstdlib>

#include "markov.h"
#include "common.h"

//--------------------------------------------------
bool MarkovCod::Set( int order, double probN )
{
	bool x;

	x = AbsToFirst( d1a, &d1, order );    if (!x) return false;
	x = AbsToFirst( d2a, &d2, order );    if (!x) return false;
	x = AbsToFirst( d3a, &d3, order );    if (!x) return false;
	
	x = RevFromDir( d1a, &r1a, order );   if (!x) return false;
	x = RevFromDir( d2a, &r2a, order );   if (!x) return false;
	x = RevFromDir( d3a, &r3a, order );   if (!x) return false;

	x = AbsToFirst( r1a, &r1, order );    if (!x) return false;
	x = AbsToFirst( r2a, &r2, order );    if (!x) return false;
	x = AbsToFirst( r3a, &r3, order );    if (!x) return false;

	// Adding support for chains of all orders

	if ( ddd_1 || ddd_2 || ddd_3 ) { printf( "error in code\n" ); exit(1); }
	ddd_1 = new double* [ order + 1 ];
	ddd_2 = new double* [ order + 1 ];
	ddd_3 = new double* [ order + 1 ];
	if( !ddd_1 || !ddd_2 || !ddd_3 ) { printf( "error, out of memory\n" ); exit(1); }

	memset( ddd_1, 0, (order+1)*sizeof(double*) );
	memset( ddd_2, 0, (order+1)*sizeof(double*) );
	memset( ddd_3, 0, (order+1)*sizeof(double*) );

	// copy absolute value to highest order
	ddd_1[order] = new double [ (4<<(2*order)) + 1 ];
	ddd_2[order] = new double [ (4<<(2*order)) + 1 ];
	ddd_3[order] = new double [ (4<<(2*order)) + 1 ];
	if( !ddd_1[order] || !ddd_2[order] || !ddd_3[order] ) { printf( "error, out of memory\n" ); exit(1); }

	for( int i = 0; i < 4<<(2*order); ++i )
	{
		ddd_1[order][i] = d1a[i];
		ddd_2[order][i] = d2a[i];
		ddd_3[order][i] = d3a[i];
	}

	// move from high to low one step at a time
	for ( int i = order - 1; i >= 0; --i )
	{
		x = HiToLo( ddd_1[i+1], ddd_1 + i, i + 1 );  if (!x) return false;
		x = HiToLo( ddd_2[i+1], ddd_2 + i, i + 1 );  if (!x) return false;
		x = HiToLo( ddd_3[i+1], ddd_3 + i, i + 1 );  if (!x) return false;
	}

	// set probability for letter 'N'
	for ( int i = 0; i <= order; ++i )
	{
		ddd_1[i][4<<(2*i)] = probN;
		ddd_2[i][4<<(2*i)] = probN;
		ddd_3[i][4<<(2*i)] = probN;
	}

	// move to relative probabilities
	for ( int i = order - 1; i >= 0; --i )
	{
		AbsoluteToRelativeRight( ddd_1[i], i );
		AbsoluteToRelativeRight( ddd_2[i], i );
		AbsoluteToRelativeRight( ddd_3[i], i );
	}
	//log the values
	for ( int i = 0; i <= order; ++i )
	{
		Log( ddd_1[i], (4<<(2*i)) + 1 );
		Log( ddd_2[i], (4<<(2*i)) + 1 );
		Log( ddd_3[i], (4<<(2*i)) + 1 );
	}
	// end of all order section

	int n = 4<<(2*order);

	d1a[n] = d2a[n] = d3a[n] = probN;
	r1a[n] = r2a[n] = r3a[n] = probN;

	d1[n] = d2[n] = d3[n] = probN;
	r1[n] = r2[n] = r3[n] = probN;

	++n;
	
	Log( d1a, n ); Log( d2a, n ); Log( d3a, n );
	Log( d1,  n ); Log( d2,  n ); Log( d3,  n );
	Log( r1a, n ); Log( r2a, n ); Log( r3a, n );
	Log( r1,  n ); Log( r2,  n ); Log( r3,  n );

	dph0.mat = d1;
	dph1.mat = d2;
	dph2.mat = d3;

	dph0.matA = d1a;
	dph1.matA = d2a;
	dph2.matA = d3a;

	dph0.next = &dph1;
	dph1.next = &dph2;
	dph2.next = &dph0;

	switch ( order%3 ) {
	case 0: 
		rph0.mat = r1;
		rph1.mat = r2;
		rph2.mat = r3;

		rph0.matA = r1a;
		rph1.matA = r2a;
		rph2.matA = r3a;
		break;
	case 1:
		rph0.mat = r3;
		rph1.mat = r1;
		rph2.mat = r2;

		rph0.matA = r3a;
		rph1.matA = r1a;
		rph2.matA = r2a;
		break;
	case 2:
		rph0.mat = r2;
		rph1.mat = r3;
		rph2.mat = r1;

		rph0.matA = r2a;
		rph1.matA = r3a;
		rph2.matA = r1a;
		break;
	}

	rph0.next = &rph2;
	rph2.next = &rph1;
	rph1.next = &rph0;

	dirArr[0] = &dph0;
	dirArr[1] = &dph1;
	dirArr[2] = &dph2;

	revArr[0] = &rph2;
	revArr[1] = &rph0;
	revArr[2] = &rph1;

	return true;
};
//--------------------------------------------------
bool MarkovIntr::Set( int order, double probN )
{
	bool x;

	x = AbsToFirst( da, &d, order );    if (!x) return false;
	x = RevFromDir( da, &ra, order );   if (!x) return false;
	x = AbsToFirst( ra, &r, order );    if (!x) return false;

	// Adding support for chains of all orders

	if ( ddd ) { printf( "error in code\n" ); exit(1); }
	ddd = new double* [ order + 1 ];
	if( !ddd ) { printf( "error, out of memory\n" ); exit(1); }
	memset( ddd, 0, (order+1)*sizeof(double*) );
	
	// copy absolute value to highest order
	ddd[order] = new double [ (4<<(2*order)) + 1 ];
	if( !ddd[order] ) { printf( "error, out of memory\n" ); exit(1); }
	for( int i = 0; i < 4<<(2*order); ++i )
		ddd[order][i] = da[i];

	// move from high to low one step at a time
	for ( int i = order - 1; i >= 0; --i )
	{
		x = HiToLo( ddd[i+1], ddd + i, i + 1 );  if (!x) return false;
	}

	// set probability for letter 'N'
	for ( int i = 0; i <= order; ++i )
		ddd[i][4<<(2*i)] = probN;

	// move to relative probabilities
	for ( int i = order - 1; i >= 0; --i )
		AbsoluteToRelativeRight( ddd[i], i );

	//log the values
	for ( int i = 0; i <= order; ++i )
		Log( ddd[i], (4<<(2*i)) + 1 );

	// end of all order section

	int size = 4<<(2*order);

	da[size] = d[size] = probN;
	ra[size] = r[size] = probN;

	++size;

	Log( da, size ); Log( d, size ); 
	Log( ra, size ); Log( r, size );

	return true;
};
//--------------------------------------------------
bool MarkovNon::Set( int order, double probN )
{
	if ( !AbsToFirst( na, &n, order ) )
		return false;

	// Adding support for chains of all orders

	if ( ddd ) { printf( "error in code\n" ); exit(1); }
	ddd = new double* [ order + 1 ];
	if( !ddd ) { printf( "error, out of memory\n" ); exit(1); }
	memset( ddd, 0, (order+1)*sizeof(double*) );
	
	// copy absolute value to highest order
	ddd[order] = new double [ (4<<(2*order)) + 1 ];
	if( !ddd[order] ) { printf( "error, out of memory\n" ); exit(1); }
	for( int i = 0; i < 4<<(2*order); ++i )
		ddd[order][i] = na[i];

	// move from high to low one step at a time
	for ( int i = order - 1; i >= 0; --i )
	{
		bool x = HiToLo( ddd[i+1], ddd + i, i + 1 );  if (!x) return false;
	}

	// set probability for letter 'N'
	for ( int i = 0; i <= order; ++i )
		ddd[i][4<<(2*i)] = probN;

	// move to relative probabilities
	for ( int i = order - 1; i >= 0; --i )
		AbsoluteToRelativeRight( ddd[i], i );

	//log the values
	for ( int i = 0; i <= order; ++i )
		Log( ddd[i], (4<<(2*i)) + 1 );

	// end of all order section

	int size = 4<<(2*order);

	na[size] = n[size]= probN;

	++size;

	Log( na, size ); Log( n, size ); 

	return true;
};
//--------------------------------------------------

