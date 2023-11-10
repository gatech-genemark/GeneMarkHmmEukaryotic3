//**************************************************
// file: precalc.h
// project: ehmm3
//**************************************************

#ifndef PRECALC_H__
#define PRECALC_H__

#include <stdio.h>

#include "map.h"
#include "common.h"
#include "site.h"
#include "markov.h"
#include "model.h"
#include "sequence.h"

//--------------------------------------------------
class Precalc
{
public:

	int order;
	Map* map;
	int* adr;
	int* extr;
	double maskp;
	
	void Init( int i_order, Map* i_map, int* i_adr, double i_maskp )
	{
		order = i_order;
		map = i_map;
		adr = i_adr;
		maskp = i_maskp;
		extr = 0;
		enforce_this = 0;
	};

	void Init( int i_order, Map* i_map, int* i_adr, double i_maskp, int* i_data )
	{
		order = i_order;
		map = i_map;
		adr = i_adr;
		maskp = i_maskp;
		extr = i_data;
		enforce_this = 0;
	};

	void PreIntergenic( double* data, MarkovNon* non );
	void PreIntron( double* data, MarkovIntr* non );
	void PreCoding( double* data, MarkovCod* cod );
	void PreSites( double* data, Sequence* s, Model* mod );
	void PreAccForBP( double* data, Sequence* s, Model* mod );

	void PrintSites( double* data );

private:

	int enforce_this;

	double GetMarkov( int j, int i, double* mat );
	double GetMarkovIni(int j, int i, double* mat, double* matA );

	void Left( double* data, double* mat, char ph, int key_L, int sh_L );

	void Right( double* data, double* mat, double* matA, char ph,
								int key_L, int key_R, int sh_L, int sh_R );
		
	double GetMarkov( int j, int i, char j_ph, MarkovCod* cod, char strand );
	double GetMarkovIni( int j, int i, char j_ph, MarkovCod* cod, char strand );

	int GetPos( MapType* map, int key1, int key2, int shift1, int shift2 );
	
	void LeftC( double* data, char fr, int leftKey1, int leftKey2, int leftShift1, int leftShift2, MarkovCod* cod, char strand );
	void RightC( double* data, char fr, int leftKey1, int leftKey2, int leftShift1, int leftShift2,
		int rightKey1, int rightKey2, int rightShift1, int rightShift2, MarkovCod* cod, char strand );
};
//--------------------------------------------------
#endif  // PRECALC_H__

