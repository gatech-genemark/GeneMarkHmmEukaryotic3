//**************************************************
// file: length.h
// project: ehmm3
//**************************************************

#ifndef LENGTH_H__
#define LENGTH_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string>
#include <fstream>
#include <iostream>

#include "common.h"

//--------------------------------------------------
// in file with parameters
//   *_min -> int
//   *_max -> int
//   *_tag -> label
//   label -> array
//   index -> value
//   where:
//        index is state duration (starts with 1)
//        value is some double >= 0
//
//   min 5
//   max 7
//   tag label
//   label
//   5 0.3
//   6 0.4
//   7 0.3
//
// min, max and data are set in model.cpp
//--------------------------------------------------
class Length
{
public:

	Length() { min = 0; max = 0; min_partial = 0; max_partial = 0; data = 0; partial = 0; data_ptr = 0; partial_ptr = 0; };
	~Length() { delete [] data; delete [] partial; };

	int max;
	int min;
	int min_partial;
	int max_partial;
	double* data;
	double* partial;

	// legacy for Get* functions
	int min_ptr;
	int max_ptr;
	int min_partial_ptr;
	int max_partial_ptr;
	double* data_ptr;
	double* partial_ptr;

	// only these two Get* function are called from viterbi

	// "len" is  (C-style length + shift )
	// if shift is "1", them "len" is true state duration 

	double Get( int len )
	{
		if( min_ptr < len && len < max_ptr )
			return  data_ptr[ len ];
		else
			return log0;
	};

	double GetPartial( int len )
	{
		if( min_partial_ptr < len && len < max_partial_ptr ) 
			return partial_ptr[len];
		else
			return log0;
	} ;

	//------------------------------------------------------
	// -- prepare class for work ---

	void PrepareLength( int shift, int shift_partial )
	{
		SetMinMaxPartial();
		AllocatePartial();

		EnforceMin( data, min, 0 );
		EnforceMin( partial, min_partial, 0 );

		NormProb( data, max );

		SetPartialFromDuration();

		NormProb( partial, max_partial );

		LogProb( data, max );
		LogProb( partial, max_partial );

		// -- set ptr's
		//   data   1  2  3  4  5  6  7  8
		//          -  -  -  - min * max
		//          0  1  2  3  4  5  6  7
		//                   |min_ptr    |max_ptr

		min_ptr = min - 1 - 1;
		max_ptr = max - 1 + 1 ;
		min_partial_ptr = min_partial - 1 - 1;
		max_partial_ptr = max_partial - 1 + 1;

		// data, partial and min*, max* are now C-style index
		// shift to match with "len"
		//
		//  C-style = len - shift
		//    p[C-style] = p[ len - shift ] = (p - shift)[len]
		//    a < C-style < b
		//    a < len - shift < b
		//    a + shift < len < b + shift

		data_ptr = data - shift;
		max_ptr += shift;
		min_ptr += shift;
		
		partial_ptr = partial - shift_partial;
		min_partial_ptr += shift_partial;
		max_partial_ptr += shift_partial;
	};

	//   data   1  2  3  4  5  6  7
	//          -  -  -  - min * max
	// partial min *  *  *  * max

	void SetMinMaxPartial(void)
	{
		min_partial = 1;
		max_partial = max - 1;
	};

	void AllocatePartial(void)
	{
		delete [] partial;
		partial = 0;
		partial = new double [max_partial];
		if (!partial) { printf("out of memory\n"); exit(1); }
	};

	//   1  2  3  4  5  6  7
	//   -  -  -  -  |low_bound

	void EnforceMin( double* ptr, int low_bound, double value )
	{
		for( int i = 0; i < low_bound - 1; ++i )
			ptr[i] = value;
	};

	void SetConstantDistr( double* ptr, int low_bound, int high_bound, double x )
	{
		for( int i = low_bound - 1; i < high_bound; ++i )
			ptr[i] = x;
	};

	// partial(i) = sum( data(i) : data(max) )

	void SetPartialFromDuration(void)
	{
		double sum = data[max-1];
		for ( int i = max_partial - 1; i >= 0; --i )
		{
			sum += data[i];
			partial[i] = sum;
		}
	};

	void NormProb( double* ptr, int high_bound )
	{
		double sum = 0;
		for( int i = 0; i < high_bound; ++i )
			sum += ptr[i];

		if (!sum) { printf("error in length distribution\n"); exit(1); }

		for( int i = 0; i < high_bound; ++i )
			ptr[i] = ptr[i]/sum;
	};

	void LogProb( double* ptr, int high_bound  )
	{
		for ( int i = 0; i < high_bound; ++i )
		{
			if ( ptr[i] )  ptr[i] = log( ptr[i] );
			else           ptr[i] = log0;
		}
	};

	//------------------------------------------------------

	void FConstant( double x )
	{
		 SetConstantDistr( data, min, max, x );
	};
	
	// Exp and Gamma for long durations can go out of presicion
	// move to 

	void Exp( double exp_decay )
	{
		for( int i = min - 1; i < max; ++i ) 
			data[i] = exp( -i/exp_decay );
	};

	void Gamma( double gamma_a, double gamma_b )
	{
		for( int i = min - 1; i < max; ++i ) 
			data[i] = pow( (double)i, gamma_a )*exp( -i/gamma_b );
	};

	//------------------------------------------------------

	// for debug
	void Print( const char* name )
	{
		std::ofstream log_file( name );
		if ( !log_file.is_open() ) { std::cout << PROG_NAME << " : error, opening file " << name << "\n";  exit(1); }

		log_file << "# min=" << min << " max=" << max << "\n";
		log_file << "# min_partial=" << min_partial << " max_partial=" << max_partial << "\n";

		int i = 0;
		while( i < max )
		{
			if ( i < max_partial )
				log_file << i << "\t" << data[i] << "\t" << partial[i] << "\n";
			else
				log_file << i << "\t" << data[i] << "\n";

			++i;
		}

		log_file.close();
	};

};
//--------------------------------------------------
#endif  // LENGTH_H__

