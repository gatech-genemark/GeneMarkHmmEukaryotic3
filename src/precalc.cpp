//**************************************************
// file: precalc.cpp
// project: ehmm3
//**************************************************

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

using std::cout;
using std::vector;
using std::string;
using std::ofstream;

#include <memory.h>
#include <stdio.h>

#include "precalc.h"

//--------------------------------------------------
void Precalc::PreCoding( double* data, MarkovCod* cod )
{
	for ( int i = 0; i < 3; ++i ) {

		LeftC( data, i, dirStart, dirAc, 3, 1, cod, '+' );
		RightC( data, i, dirStart, dirAc, 3, 1, dirEnd, dirDn, -2, 0, cod, '+' );

		LeftC( data, i, revEnd, revDn, 3, 1, cod, '-' );
		RightC( data, i, revEnd, revDn, 3, 1, revStart, revAc, -2, 0, cod, '-' );
	}
};
//--------------------------------------------------
double Precalc::GetMarkov( int j, int i, char j_ph, MarkovCod* cod, char strand )
{
	double x = 0;
	int save_j = j;

	if ( strand == '+' ) {

		cod->ResetDir(j_ph);
		for ( ; j < i; ++j )
			x += cod->GetDir( adr[j] );

	} else {

		cod->ResetRev(j_ph);
		for ( ; j < i; ++j )
			x += cod->GetRev( adr[j] );
	}

	if (extr)
	{
		j = save_j;
//		bool mark_found = false;

		for ( ; j < i; ++j )
		{
//			if ( extr[j] & nonMark )
//				mark_found = true;

			if (extr[j] & maskMark)
				x += maskp;
		}

//		if ( mark_found )
//			x += MINUS;
	}

	return x;
};
//--------------------------------------------------
double Precalc::GetMarkovIni( int j, int i, char j_ph, MarkovCod* cod, char strand )
{
	double x = 0;
	int save_j = j;

	if ( strand == '+' ) {

		cod->ResetDir(j_ph);
		i -= ( order + 1 );
		for ( ; j < i; ++j )
			x += cod->GetDir( adr[j] );
		x += cod->GetDirA( adr[j] );

	} else {

		cod->ResetRev(j_ph);
		i -= ( order + 1 );
		for ( ; j < i; ++j )
			x += cod->GetRev( adr[j] );
		x += cod->GetRevA( adr[j] );
	}

	if (extr)
	{
		j = save_j;
//		bool mark_found = false;

		for ( ; j < i; ++j )
		{
//			if ( extr[j] & nonMark )
//				mark_found = true;

			if (extr[j] & maskMark)
				x += maskp;
		}

//		if ( mark_found )
//			x += MINUS;
	}

	return x;
};
//--------------------------------------------------
int Precalc::GetPos( MapType* map, int key1, int key2, int shift1, int shift2 )
{
	if ( !map )
		return -1;
	else if ( map->key & key1 )
		return map->pos + shift1;
	else if ( map->key & key2 )
		return map->pos + shift2;
	else
		return -1;
}
//--------------------------------------------------
void Precalc::LeftC( double* data, char fr, int leftKey1, int leftKey2, int leftShift1, int leftShift2, MarkovCod* cod, char strand )
{
	int posLeft;
	int posRight;

	map->ResetMap();
	map->SetCurrentToLast();
	MapType* right = map->FindFr( leftKey1 + leftKey2, fr );
	posRight = GetPos( right, leftKey1, leftKey2, leftShift1, leftShift2 );
	MapType* left = right;

	while( left ) {

		left = map->FindFr( leftKey1 + leftKey2, fr );
		if (left) {

			posLeft = GetPos( left, leftKey1, leftKey2, leftShift1, leftShift2 );
			data[ map->Index(left) ] = GetMarkov( posLeft, posRight, left->ph, cod, strand );
			right = left;
			posRight = posLeft;
		}
	}

	if ( right ) {
		left = map->FindFrIni( leftKey1 + leftKey2, fr );
		posLeft = GetPos( left, leftKey1, leftKey2, 0, 0 );
		data[ map->Index(left) ] = GetMarkov( posLeft, posRight, left->ph, cod, strand );
	}
};
//--------------------------------------------------
void Precalc::RightC( double* data, char fr, int leftKey1, int leftKey2, int leftShift1, int leftShift2,
		int rightKey1, int rightKey2, int rightShift1, int rightShift2, MarkovCod* cod, char strand )
{
	int posLeft;
	int posRight;

	map->ResetMap();
	MapType* right = map->FindFrTer( rightKey1 + rightKey2, fr );
	posRight = GetPos( right, rightKey1, rightKey2, 1, 1 );
	MapType* left;

	map->SetCurrentToLast();
	while (1) {

		left = map->FindFr( leftKey1 + leftKey2, fr );
		posLeft = GetPos( left, leftKey1, leftKey2, leftShift1, leftShift2 );

		if ( !left ) {
			left = map->FindFrIni( leftKey1 + leftKey2, fr );
			posLeft = GetPos( left, leftKey1, leftKey2, 0, 0 );

			if ( ( posRight - posLeft > order )&&( posRight > -1 ) ) {
				data[ map->Index(right) ] = GetMarkovIni( posLeft, posRight, left->ph, cod, strand );
			}
			break;
		}

		if ( ( posRight - posLeft > order )&&( posRight > -1 ) ) {
			data[ map->Index(right) ] = GetMarkovIni( posLeft, posRight, left->ph, cod, strand );
			break;
		}
	}

	map->SetCurrentToLast();
	while ( right = map->FindFr( rightKey1 + rightKey2, fr ) ) {

		posRight = GetPos( right, rightKey1, rightKey2, rightShift1, rightShift2 );

		while (1) {

			left = map->FindFr( leftKey1 + leftKey2, fr );
			posLeft = GetPos( left, leftKey1, leftKey2, leftShift1, leftShift2 );

			if ( !left ) {
				left = map->FindFrIni( leftKey1 + leftKey2, fr );
				posLeft = GetPos( left, leftKey1, leftKey2, 0, 0 );

				if ( ( posRight - posLeft > order )&&( posRight > -1 ) ) {
					data[ map->Index(right) ] = GetMarkovIni( posLeft, posRight, left->ph, cod, strand );
				}
				break;
			}

			if ( ( posRight - posLeft > order )&&( posRight > -1 ) ) {
				data[ map->Index(right) ] = GetMarkovIni( posLeft, posRight, left->ph, cod, strand );
				break;
			}
		}
		map->SetCurrent( right );
	}
};
//--------------------------------------------------
void Precalc::Left( double* data, double* mat, char ph, int key_L, int sh_L )
{
	map->ResetMap();
	map->SetCurrentToLast();

	MapType* right = map->FindPh( key_L, ph );
	MapType* left = right;
	
	while( left )
	{
		left = map->FindPh( key_L, ph );

		if (left)
		{
			data[ map->Index(left) ] = 
				GetMarkov( left->pos + sh_L, right->pos + sh_L, mat );
			right = left;
		}
	}

	if ( right )
	{
		left = map->FindPhIni( key_L, ph );

		data[ map->Index(left) ] = 
			GetMarkov( left->pos, right->pos + sh_L, mat );
	}
};
//--------------------------------------------------
void Precalc::Right( double* data, double* mat, double* matA, char ph,
		int key_L, int key_R, int sh_L, int sh_R )
{
	map->ResetMap();

	MapType* right = map->FindPhTer( key_R, ph );
	MapType* left;

	map->SetCurrentToLast();
	while(1)
	{
		left = map->FindPh( key_L, ph );

		if ( !left )
		{
			left = map->FindPhIni( key_L, ph );
			
			if ( right->pos - left->pos >= order )
			{
				data[ map->Index(right) ] = 
					GetMarkovIni( left->pos, right->pos + 1, mat, matA );
			}
			else
			{
				data[ map->Index(right) ] = log0;
			}
			break;
		}

		if ( right->pos - (left->pos + sh_L) >= order )
		{
			data[ map->Index(right) ] = 
				GetMarkovIni( left->pos + sh_L, right->pos + 1, mat, matA );
			break;
		}
	}

	map->SetCurrentToLast();
	while ( right = map->FindPh( key_R, ph ) )
	{
		while(1)
		{
			left = map->FindPh( key_L, ph );

			if ( !left )
			{
				left = map->FindPhIni( key_L, ph );

				if ( (right->pos + sh_R) - left->pos > order )
				{
					data[ map->Index(right) ] = 
						GetMarkovIni( left->pos, right->pos + sh_R, mat, matA );
				}
				else
				{
					data[ map->Index(right) ] = log0;
				}
				break;
			}

			if ( (right->pos + sh_R) - (left->pos + sh_L) > order )
			{
				data[ map->Index(right) ] = 
					GetMarkovIni( left->pos + sh_L, right->pos + sh_R, mat, matA );
				break;
			}
		}

		map->SetCurrent( right );
	}
};
//--------------------------------------------------
void Precalc::PreIntron( double* data, MarkovIntr* non )
{
//	if ( extr )
//		enforce_this = interMark;

	for ( int i = 0; i < 3; ++i ) {

		Left( data, non->d, i, dirDn, 2 );
		Right( data, non->d, non->da, i, dirDn, dirAc, 2, -1 );

		Left( data, non->r, i, revAc, 2 );
 		Right( data, non->r, non->ra, i, revAc, revDn, 2, -1 );
	}
};
//--------------------------------------------------
void Precalc::PreIntergenic( double* data, MarkovNon* non )
{
//	if ( extr )
//		enforce_this = intronMark;

	Left( data, non->n, 0, dirEnd + revStart, 1 );
	Right( data, non->n, non->na, 0, dirEnd + revStart, dirStart + revEnd, 1, 0 );
};
//--------------------------------------------------
double Precalc::GetMarkov( int j, int i, double* mat )
{
	double x = 0;
	int save_j = j;

	for ( ; j < i; ++j )
		x += mat[ adr[j] ];

	if (extr && enforce_this )
	{
		j = save_j;

//		bool mark_found = false;

//		for ( ; j < i; ++j )
//		{
//			if ( extr[j] & enforce_this )
//				mark_found = true;
//		}

//		if ( mark_found )
//			x += MINUS;
	}

	return x;
};
//--------------------------------------------------
double Precalc::GetMarkovIni( int j, int i, double* mat, double* matA )
{
	i -= ( order + 1 );
	double x = matA[ adr[i] ];

	int save_j = j;

	for ( ; j < i; ++j )
		x += mat[ adr[j] ];

	if (extr && enforce_this )
	{
		j = save_j;
		i += (order + 1);
//		bool mark_found = false;

//		for ( ; j < i; ++j )
//		{
//			if ( extr[j] & enforce_this )
//				mark_found = true;
//		}

//		if ( mark_found )
//			x += MINUS;
	}

	return x;
};
//--------------------------------------------------
void Precalc::PreSites( double* data, Sequence* s, Model* mod )
{
	char* seq = s->seq;
	int max = s->size - 30;

	int key;
	int pos;
	int ph;
 
	double tr_prob = 0;
	if ( mod->BpModelOn )
		tr_prob = mod->aroundBp;

	for( int i = 0; i < map->size; ++i )
	{
		pos = map->map[i].pos;

		if ( pos > 30 && pos < max )
		{
			if (s->data && map->map[i].data)
			{
				if (map->map[i].key & map->map[i].data)
				{ 
					data[i] += PLUS;
//					std::cout << "plus " << pos << std::endl;
				}
				else if (map->map[i].key & anyEnd)
				{
//					std::cout << "stop " << pos << std::endl;
				}
				else if (map->map[i].data & intronMark)
				{
					data[i] += MINUS_IN_INTRON;
//					std::cout << "minus " << pos << std::endl;
				}
			}

			key = map->map[i].key;
			ph = map->map[i].ph;

			if (key & dirDn)
			{
				if (*(seq + pos + 1) == T)
				{
					data[i] += mod->toGTdonor;

					if (ph == 0)
						data[i] += mod->don_0.GetValueDirNew(seq + pos, 0);
					else if (ph == 1)
						data[i] += mod->don_1.GetValueDirNew(seq + pos, 0);
					else if (ph == 2)
						data[i] += mod->don_2.GetValueDirNew(seq + pos, 0);
				}
				else if (*(seq + pos + 1) == C)
				{
					data[i] += mod->toGCdonor;

					if (ph == 0)
						data[i] += mod->don_gc_0.GetValueDirNew(seq + pos, 0);
					else if (ph == 1)
						data[i] += mod->don_gc_1.GetValueDirNew(seq + pos, 0);
					else if (ph == 2)
						data[i] += mod->don_gc_2.GetValueDirNew(seq + pos, 0);
				}
			}
			else if (key & dirAc)
			{
				if (ph == 0)
					data[i] += mod->acc_0.GetValueDirNew(seq + pos, 1) + tr_prob;
				else if (ph == 1)
					data[i] += mod->acc_1.GetValueDirNew(seq + pos, 1) + tr_prob;
				else if (ph == 2)
					data[i] += mod->acc_2.GetValueDirNew(seq + pos, 1) + tr_prob;
			}
			else if (key & dirStart)
			{
				data[i] += mod->ini.GetValueDirNew(seq + pos, 0);
			}
			else if (key & dirEnd)
			{
				if (( *(seq+pos-1) == A )&&( *(seq+pos) == A ))
					data[i] += mod->terTAA.GetValueDirNew( seq + pos, 2 );
				else if (( *(seq+pos-1) == A )&&( *(seq+pos) == G ))
					data[i] += mod->terTAG.GetValueDirNew( seq + pos, 2 );
				else if (( *(seq+pos-1) == G )&&( *(seq+pos) == A ))
					data[i] += mod->terTGA.GetValueDirNew( seq + pos, 2 );
			}
			else if (key & revDn)
			{
				if (*(seq + pos - 1) == A)
				{
					data[i] += mod->toGTdonor;

					if (ph == 0)
						data[i] += mod->don_0.GetValueRevNew(seq + pos, 1);
					else if (ph == 1)
						data[i] += mod->don_1.GetValueRevNew(seq + pos, 1);
					else if (ph == 2)
						data[i] += mod->don_2.GetValueRevNew(seq + pos, 1);
				}
				else if (*(seq + pos - 1) == G)
				{
					data[i] += mod->toGCdonor;

					if (ph == 0)
						data[i] += mod->don_gc_0.GetValueRevNew(seq + pos, 1);
					else if (ph == 1)
						data[i] += mod->don_gc_1.GetValueRevNew(seq + pos, 1);
					else if (ph == 2)
						data[i] += mod->don_gc_2.GetValueRevNew(seq + pos, 1);
				}
			}
			else if (key & revAc)
			{
				if (ph == 0)
					data[i] += mod->acc_0.GetValueRevNew(seq + pos, 2) + tr_prob;
				else if (ph == 1)
					data[i] += mod->acc_1.GetValueRevNew(seq + pos, 2) + tr_prob;
				else if (ph == 2)
					data[i] += mod->acc_2.GetValueRevNew(seq + pos, 2) + tr_prob;
			}
			else if (key & revStart)
			{
				data[i] += mod->ini.GetValueRevNew(seq + pos, 1);
			}
			else if (key & revEnd)
			{
				if (( *(seq+pos) == T )&&( *(seq+pos+1) == T ))
					data[i] += mod->terTAA.GetValueRevNew( seq + pos, 3 );
				else if (( *(seq+pos) == C )&&( *(seq+pos+1) == T ))
					data[i] += mod->terTAG.GetValueRevNew( seq + pos, 3 );
				else if (( *(seq+pos) == T )&&( *(seq+pos+1) == C ))
					data[i] += mod->terTGA.GetValueRevNew( seq + pos, 3 );
			}

//			printf( "%d\t%d\t%f\n", pos, key, data[i] );
		}
	}
};
//--------------------------------------------------
void Precalc::PreAccForBP( double* data, Sequence* s, Model* mod )
{
	char* seq = s->seq;
	int max = s->size - 30;

	int key;
	int pos;
	int ph;

	double tr_prob = mod->toBp;

	for( int i = 0; i < map->size; ++i ) {

		data[i] = log0;
		pos = map->map[i].pos;

		if ( pos > 30 && pos < max )
		{
			if (s->data && map->map[i].data)
			{
				if (map->map[i].key & map->map[i].data)
				{
					data[i] += PLUS;
				}
				else if (map->map[i].key & anyEnd)
				{
					;
				}
				else if (map->map[i].data & intronMark)
				{
					data[i] += MINUS_IN_INTRON;
				}
			}

			key = map->map[i].key;
			ph = map->map[i].ph;

			if ((key & dirAc)||(key & revAc))
			{
				if (( key & dirAc )&&( ph == 0 ))
					data[i] = mod->acc_bp_0.GetValueDirNew( seq + pos, 1 ) + tr_prob;
				else if (( key & dirAc )&&( ph == 1 ))
					data[i] = mod->acc_bp_1.GetValueDirNew( seq + pos, 1 ) + tr_prob;
				else if (( key & dirAc )&&( ph == 2 ))
					data[i] = mod->acc_bp_2.GetValueDirNew( seq + pos, 1 ) + tr_prob;

				else if (( key & revAc )&&( ph == 0 ))
					data[i] = mod->acc_bp_0.GetValueRevNew( seq + pos, 2 ) + tr_prob;
				else if (( key & revAc )&&( ph == 1 ))
					data[i] = mod->acc_bp_1.GetValueRevNew( seq + pos, 2 ) + tr_prob;
				else if (( key & revAc )&&( ph == 2 ))
					data[i] = mod->acc_bp_2.GetValueRevNew( seq + pos, 2 ) + tr_prob;
			}
		}
	}
};
//--------------------------------------------------
void Precalc::PrintSites( double* data )
{
	if ( !map )
	{
		std::cout << "no map provided\n";
		return;
	}

	std::string name( "site.log" );
	std::ofstream log_file( name.c_str() );

	if ( !log_file.is_open() )
	{
		std::cout << PROG_NAME << " : error, opening file " << name << "\n"; 
		exit(1);
	}

	string label;
	int current_key;
	bool isIni;
	bool isTer;

	MapType* m = map->map;
	int size = map->size;

	for( int i = 0; i < size; ++i )
	{
		current_key = m[i].key;

		if ( current_key & mapIni )
		{
			isIni = true;
			current_key -= mapIni;
		}
		else
			isIni = false;

		if ( current_key & mapTer )
		{
			isTer = true;
			current_key -= mapTer;
		}
		else
			isTer = false;

		switch( current_key )
		{
			case dirStart:
				label = "dirStart";
				break;
			case dirEnd:
				label = "dirEnd";
				break;
			case dirDn:
				label = "dirDn";
				break;
			case dirAc:
				label = "dirAc";
				break;
			case revStart:
				label = "revStart";
				break;
			case revEnd:
				label = "revEnd";
				break;
			case revDn:
				label = "revDn";
				break;
			case revAc:
				label = "revAc";
				break;
			default:
				label = "unknown";
				break;
		}

		if ( isIni )
			label += "_Ini";

		if ( isTer )
			label += "_Ter";

		log_file << i << " " << label << " " << m[i].key << " " << m[i].pos << " " << (int)m[i].fr << " " << (int)m[i].ph << " " << m[i].data << " " << data[i] << "\n";
	}

	log_file.close();
};
//--------------------------------------------------

