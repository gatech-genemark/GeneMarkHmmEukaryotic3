//**************************************************
// file: map.h
// project: ehmm3
// Alex Lomsadze
//**************************************************

#ifndef MAP_H__
#define MAP_H__

#include "model.h"
#include "common.h"

//--------------------------------------------------
class Map
{
public:

	Map() { size = 0; map = 0; reserved = 0; };
	~Map() { delete [] map; };

	MapType* map;
	int size;
	
	bool Init( char* seq, int* data, int seqSize );
	bool Init( char* seq, int* data, int seqSize, Model* mod );

	void PrintMap();
	void PrintMapCount();

	//---
	void SetCurrentToLast() { current = last; }
	void SetCurrentToFirst() { current = first; }
	void SetCurrent( MapType* ptr ) { current = ptr; }
	int Index( MapType* ptr ) { return ptr - map; }
	
	void ResetMap()
	{
		first = map + 12;
		last = map + size - 13;
		current = 0;
	}

	MapType* FindPh( int key, char ph )
	{
		if ( !current ) return current;

		for ( --current ; current > first; --current ) {
			if ( ( current->key & key )&&( current->ph == ph ) )
				return current;
		}
		return current = 0;
	}

	MapType* FindFr( int key, char fr )
	{
		if ( !current ) return current;

		for ( --current ; current > first; --current ) {
			if ( ( current->key & key )&&( current->fr == fr ) )
				return current;
		}
		return current = 0;
	}

	MapType* FindPhIni( int key, char ph )
	{
		for ( current = map + 13 - 1 ; current >= map; --current ) {
			if ( ( current->key & key )&&( current->ph == ph ) )
				return current;
		}
		return current = 0;
	}

	MapType* FindFrIni( int key, char fr )
	{
		for ( current = map + 13 - 1 ; current >= map; --current ) {
			if ( ( current->key & key )&&( current->fr == fr ) )
				return current;
		}
		return current = 0;
	}

	MapType* FindPhTer( int key, char ph )
	{
		for ( current = map + size - 1 ; current >= map + size - 13; --current ) {
			if ( ( current->key & key )&&( current->ph == ph ) )
				return current;
		}
		return current = 0;
	}

	MapType* FindFrTer( int key, char fr )
	{
		for ( current = map + size - 1 ; current >= map + size - 13; --current ) {
			if ( ( current->key & key )&&( current->fr == fr ) )
				return current;
		}
		return current = 0;
	}

	//---

private:

	int reserved;

	//---
	MapType* current;
	MapType* first;
	MapType* last;
	//---

	bool AllocMap( int n );
	bool ReAllocMap( int n );

	bool GetMapFromSeq( char* seq, int* data, int seqSize );
	bool GetMapFromSeq( char* seq, int* data, int seqSize, Model* mod );
	void SetIniTerm( int seqSize );
	void SortMap();
	MapType* CopySorted( MapType* target, MapType* source, int max, int key );

	void OverlayMapByData( int* data );
	void DeleteEmptyKeyRecords();

	void SetPoint( MapType* m, int pos, int key, int fr, int ph, int data )
	{
		m->pos = pos;
		m->key = key;
		m->fr = fr;
		m->ph = ph;
		m->data = data;
	}
};
//--------------------------------------------------
#endif  // MAP_H__

