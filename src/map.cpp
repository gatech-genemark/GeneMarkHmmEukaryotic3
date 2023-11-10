//**************************************************
// file: map.cpp
// project: ehmm3
// Alex Lomsadze
//**************************************************

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <string.h> 

using std::cout;
using std::vector;
using std::string;
using std::ofstream;

#include "map.h"

//--------------------------------------------------
bool Map::Init( char* seq, int* data, int seqSize )
{
	if ( !GetMapFromSeq( seq, data, seqSize ) )
	{
		cout << PROG_NAME << ": error, out of memory on create map\n";
		return false;
	}

	SetIniTerm( seqSize );

	if ( data )
	{
		OverlayMapByData( data );
		DeleteEmptyKeyRecords();
	}

	SortMap();

	return true;
};
//--------------------------------------------------
bool Map::Init( char* seq, int* data, int seqSize, Model* mod )
{
	if ( !GetMapFromSeq( seq, data, seqSize, mod ) )
	{
		cout << PROG_NAME << ": error, out of memory on create map\n";
		return false;
	}

	SetIniTerm( seqSize );

	if ( data )
	{
		OverlayMapByData( data );
		DeleteEmptyKeyRecords();
	}

	SortMap();

	return true;
};
//--------------------------------------------------
bool Map::AllocMap( int n )
{
	delete [] map;
	map = new MapType[ n ];
	if (!map) return false;

	reserved = n;

	return true;
};
//--------------------------------------------------
bool Map::ReAllocMap( int n )
{
	if ( !map )
		exit(1);
	
	if ( reserved > n )
		exit(1);

	MapType* tmp = new MapType[ n ];
	if (!tmp) return false;
	
	memcpy( tmp, map, reserved * sizeof( MapType ) );
	
	delete [] map;

	map = tmp;
	reserved = n;

	return true;
};
//--------------------------------------------------
// 
// fr 0   T G A                T C A
// fr 1     T G A                T C A
// fr 2       T G A                T C A
// 
// fr     1 2 0 1 2 0          1 2 0 1 2 0
// 
// seq    0 1 2 3 4 5          0 1 2 3 4 5
// 
// ph 0        |0 1 2               |2 1 0
// ph 1       0|1 2              2 1|0
// ph 2     0 1|2                  2|1 0
// 
// (seq+1)%3  |  0   1   2      |  0   1   2
// ---------  --------------    --------------
// ph 0       |  0   1   2      |  0   1   2
// ph 1       |  2   0   1      |  1   2   0
// ph 2       |  1   2   0      |  2   0   1
// 
//--------------------------------------------------
void Map::SetIniTerm( int seqSize )
{
	// initialization

	// intergenic
	SetPoint( &map[0], 0, dirEnd + mapIni, 0, 0, 0 );      // 258

	// direct intron
	SetPoint( &map[1], 0, dirDn + mapIni, 0, 0, 0 );       // 
	SetPoint( &map[2], 0, dirDn + mapIni, 0, 1, 0 );       // 260
	SetPoint( &map[3], 0, dirDn + mapIni, 0, 2, 0 );       // 

	// reverse intron
	SetPoint( &map[4], 0, revAc + mapIni, 0, 0, 0 );       // 
	SetPoint( &map[5], 0, revAc + mapIni, 0, 1, 0 );       // 384
	SetPoint( &map[6], 0, revAc + mapIni, 0, 2, 0 );       // 

	// direct coding
	SetPoint( &map[7], 0, dirAc + mapIni, 0, 0, 0 );       // 
	SetPoint( &map[8], 0, dirAc + mapIni, 2, 1, 0 );       // 264
	SetPoint( &map[9], 0, dirAc + mapIni, 1, 2, 0 );       // 

	// reverse coding
	SetPoint( &map[10], 0, revDn + mapIni, 0, 0, 0 );      // 
	SetPoint( &map[11], 0, revDn + mapIni, 1, 1, 0 );      // 320
	SetPoint( &map[12], 0, revDn + mapIni, 2, 2, 0 );      // 

	// termination
	int n = seqSize - 1;
	int m = seqSize;

	// direct intron
	SetPoint( &map[size-13], n, dirAc + mapTer, 0, 0, 0 );                 // 
	SetPoint( &map[size-12], n, dirAc + mapTer, 0, 1, 0 );                 // 520
	SetPoint( &map[size-11], n, dirAc + mapTer, 0, 2, 0 );                 // 

	//reverse intron
	SetPoint( &map[size-10], n, revDn + mapTer, 0, 0, 0 );                 // 
	SetPoint( &map[size-9],  n, revDn + mapTer, 0, 1, 0 );                 // 576
	SetPoint( &map[size-8],  n, revDn + mapTer, 0, 2, 0 );                 // 

	// direct coding
	SetPoint( &map[size-7], n, dirDn + mapTer, (0 + m%3)%3, 0, 0 );        // 
	SetPoint( &map[size-6], n, dirDn + mapTer, (2 + m%3)%3, 1, 0 );        // 516
	SetPoint( &map[size-5], n, dirDn + mapTer, (1 + m%3)%3, 2, 0 );        // 

	// reverse coding
	SetPoint( &map[size-4], n, revAc + mapTer, (0 + m%3)%3, 0, 0 );        // 
	SetPoint( &map[size-3], n, revAc + mapTer, (1 + m%3)%3, 1, 0 );        // 640
	SetPoint( &map[size-2], n, revAc + mapTer, (2 + m%3)%3, 2, 0 );        // 

	// intergenic
	SetPoint( &map[size-1], n, revEnd + mapTer, 0, 0, 0 );                 // 544
};
//--------------------------------------------------
// fr     1 2 0 1 2 0
// seq    0 1 2 3 4 5
//--------------------------------------------------
bool Map::GetMapFromSeq( char* seq, int* data, int seqSize )
{
	reserved = (int) (1.5*seqSize);
	if ( !AllocMap( reserved ) )
		return false;

	int i = 13;

	for ( int pos = 1, fr = 2; pos < seqSize; ++pos, ++fr ) 
	{
		// add 13 from begining and end of sequence 
		if ( i + 4 + 13 > reserved )
		{
			if ( !ReAllocMap( reserved + 1000 ) )
				return false;
		}

		if ( fr == 3 )
			fr = 0;

		if ( seq[pos-1] == G && seq[pos] == T )
		{   // direct donor
			SetPoint( &map[i], pos-1, dirDn, (0 + (pos-1)%3)%3, 0, 0 ); ++i;
			SetPoint( &map[i], pos-1, dirDn, (2 + (pos-1)%3)%3, 1, 0 ); ++i;
			SetPoint( &map[i], pos-1, dirDn, (1 + (pos-1)%3)%3, 2, 0 ); ++i;
		}
		else if ( seq[pos-1] == A && seq[pos] == G )
		{   // direct acceptor
			SetPoint( &map[i], pos, dirAc, (0 + (pos+1)%3)%3, 0, 0 ); ++i;
			SetPoint( &map[i], pos, dirAc, (2 + (pos+1)%3)%3, 1, 0 ); ++i;
			SetPoint( &map[i], pos, dirAc, (1 + (pos+1)%3)%3, 2, 0 ); ++i;
		}
		else if ( seq[pos-1] == A && seq[pos] == C ) 
		{   // reverse donor
			SetPoint( &map[i], pos, revDn, (0 + (pos+1)%3)%3, 0, 0 ); ++i;
			SetPoint( &map[i], pos, revDn, (1 + (pos+1)%3)%3, 1, 0 ); ++i;
			SetPoint( &map[i], pos, revDn, (2 + (pos+1)%3)%3, 2, 0 ); ++i;
		}
		else if ( seq[pos-1] == C && seq[pos] == T )
		{   // reverse acceptor
			SetPoint( &map[i], pos-1, revAc, (0 + (pos-1)%3)%3, 0, 0 ); ++i;
			SetPoint( &map[i], pos-1, revAc, (1 + (pos-1)%3)%3, 1, 0 ); ++i;
			SetPoint( &map[i], pos-1, revAc, (2 + (pos-1)%3)%3, 2, 0 ); ++i;
		}

		if ( pos > 1 )
		{
			if ( seq[pos-2] == A && seq[pos-1] == T && seq[pos] == G )
			{   // direct start
				SetPoint( &map[i], pos - 2, dirStart, fr, 0, 0 ); ++i;
			}
			else if ( ( seq[pos-2] == T && seq[pos-1] == A && seq[pos] == A )
					 || ( seq[pos-2] == T && seq[pos-1] == A && seq[pos] == G )
					 || ( seq[pos-2] == T && seq[pos-1] == G && seq[pos] == A ) ) 
			{   // direct end
				SetPoint( &map[i], pos, dirEnd, fr, 0, 0 ); ++i;
			}
			else if ( seq[pos-2] == C && seq[pos-1] == A && seq[pos] == T )
			{   // reverse start
				SetPoint( &map[i], pos, revStart, fr, 0, 0 ); ++i;
			}
			else if ( ( seq[pos-2] == T && seq[pos-1] == T && seq[pos] == A )
					 || ( seq[pos-2] == C && seq[pos-1] == T && seq[pos] == A )
					 || ( seq[pos-2] == T && seq[pos-1] == C && seq[pos] == A ) ) 
			{   // reverse end
				SetPoint( &map[i], pos - 2, revEnd, fr, 0, 0 ); ++i;
			}
		}
	}

	size = i + 13;

	return true;
};
//--------------------------------------------------
bool Map::GetMapFromSeq( char* seq, int* data, int seqSize, Model* mod )
{
	reserved = (int) (1.5*seqSize);
	if ( !AllocMap( reserved ) )
		return false;

	int i = 13;

	for (int pos = 1, fr = 2; pos < seqSize; ++pos, ++fr)
	{
		// reserve 13 positions at the end of the map for termination states
		// several signals can be anchored to the same position in the sequence
		// add 20 - much more than possible
		if (i + 20 + 13 > reserved)
		{
			if (!ReAllocMap(reserved + 1000))
				return false;
		}

		if (fr == 3)
			fr = 0;

		if (seq[pos - 1] == G && seq[pos] == T)
		{   // direct donor
			SetPoint(&map[i], pos - 1, dirDn, (0 + (pos - 1) % 3) % 3, 0, 0); ++i;
			SetPoint(&map[i], pos - 1, dirDn, (2 + (pos - 1) % 3) % 3, 1, 0); ++i;
			SetPoint(&map[i], pos - 1, dirDn, (1 + (pos - 1) % 3) % 3, 2, 0); ++i;
		}
		else if (mod->GC_donor_ON && seq[pos - 1] == G && seq[pos] == C)
		{	// non canonical GC donor direct
			SetPoint(&map[i], pos - 1, dirDn, (0 + (pos - 1) % 3) % 3, 0, 0); ++i;
			SetPoint(&map[i], pos - 1, dirDn, (2 + (pos - 1) % 3) % 3, 1, 0); ++i;
			SetPoint(&map[i], pos - 1, dirDn, (1 + (pos - 1) % 3) % 3, 2, 0); ++i;
		}
		else if (seq[pos - 1] == C && seq[pos] == T)
		{   // reverse acceptor
			SetPoint(&map[i], pos - 1, revAc, (0 + (pos - 1) % 3) % 3, 0, 0); ++i;
			SetPoint(&map[i], pos - 1, revAc, (1 + (pos - 1) % 3) % 3, 1, 0); ++i;
			SetPoint(&map[i], pos - 1, revAc, (2 + (pos - 1) % 3) % 3, 2, 0); ++i;
		}
		
		if (seq[pos - 1] == A && seq[pos] == G)
		{   // direct acceptor
			SetPoint(&map[i], pos, dirAc, (0 + (pos + 1) % 3) % 3, 0, 0); ++i;
			SetPoint(&map[i], pos, dirAc, (2 + (pos + 1) % 3) % 3, 1, 0); ++i;
			SetPoint(&map[i], pos, dirAc, (1 + (pos + 1) % 3) % 3, 2, 0); ++i;
		}
		else if (mod->GC_donor_ON && seq[pos - 1] == G && seq[pos] == C)
		{	// non canonical GC donor reverse
			SetPoint(&map[i], pos, revDn, (0 + (pos + 1) % 3) % 3, 0, 0); ++i;
			SetPoint(&map[i], pos, revDn, (1 + (pos + 1) % 3) % 3, 1, 0); ++i;
			SetPoint(&map[i], pos, revDn, (2 + (pos + 1) % 3) % 3, 2, 0); ++i;
		}
		else if (seq[pos - 1] == A && seq[pos] == C)
		{   // reverse donor
			SetPoint(&map[i], pos, revDn, (0 + (pos + 1) % 3) % 3, 0, 0); ++i;
			SetPoint(&map[i], pos, revDn, (1 + (pos + 1) % 3) % 3, 1, 0); ++i;
			SetPoint(&map[i], pos, revDn, (2 + (pos + 1) % 3) % 3, 2, 0); ++i;
		}

		if (pos > 1)
		{
			if (seq[pos - 2] == A && seq[pos - 1] == T && seq[pos] == G)
			{   // direct start
				SetPoint(&map[i], pos - 2, dirStart, fr, 0, 0); ++i;
			}
			else if (mod->TAA_ON && seq[pos - 2] == T && seq[pos - 1] == A && seq[pos] == A)
			{
				SetPoint(&map[i], pos, dirEnd, fr, 0, 0); ++i;
			}
			else if (mod->TAG_ON && seq[pos - 2] == T && seq[pos - 1] == A && seq[pos] == G)
			{
				SetPoint(&map[i], pos, dirEnd, fr, 0, 0); ++i;
			}
			else if (mod->TGA_ON && seq[pos - 2] == T && seq[pos - 1] == G && seq[pos] == A)
			{
				SetPoint(&map[i], pos, dirEnd, fr, 0, 0); ++i;
			}
			else if (seq[pos - 2] == C && seq[pos - 1] == A && seq[pos] == T)
			{   // reverse start
				SetPoint(&map[i], pos, revStart, fr, 0, 0); ++i;
			}
			else if (mod->TAA_ON && seq[pos - 2] == T && seq[pos - 1] == T && seq[pos] == A)
			{
				SetPoint(&map[i], pos - 2, revEnd, fr, 0, 0); ++i;
			}
			else if (mod->TAG_ON && seq[pos - 2] == C && seq[pos - 1] == T && seq[pos] == A)
			{
				SetPoint(&map[i], pos - 2, revEnd, fr, 0, 0); ++i;
			}
			else if (mod->TGA_ON && seq[pos - 2] == T && seq[pos - 1] == C && seq[pos] == A)
			{
				SetPoint(&map[i], pos - 2, revEnd, fr, 0, 0); ++i;
			}
		}
	}

	size = i + 13;

	return true;
};
//--------------------------------------------------
void Map::SortMap()
{
	MapType* ptr = map + 13;
	MapType* end = map + size - 13;
	MapType* change_point;

	MapType store[14];

	int pos;
	int i;

	while( ptr < end )
	{
		// Count number of consecutive map members pointing to the same position on sequence
		change_point = ptr;
		pos = ptr->pos;
		while(( ptr < end )&&( ptr->pos == pos ))
			++ptr;

		i = ptr - change_point;

		// if non trivial case, change order
		if (( i != 1 )&&( i != 3 ))
		{
			memcpy( store, change_point, i*sizeof(MapType) );

			change_point = CopySorted( change_point, store, i, dirStart );
			change_point = CopySorted( change_point, store, i, dirDn );
			change_point = CopySorted( change_point, store, i, dirEnd );
			change_point = CopySorted( change_point, store, i, dirAc );
			change_point = CopySorted( change_point, store, i, revAc );
			change_point = CopySorted( change_point, store, i, revDn );
			change_point = CopySorted( change_point, store, i, revEnd );
			change_point = CopySorted( change_point, store, i, revStart );
		}
	}

	// Current site anchoring rules and order of sequence parsing can lead to wrong order of Start markers
	for (ptr = map + 14; ptr < end; ++ptr)
	{
		if ((ptr->key == dirStart) && ((ptr - 1)->key == revStart))
		{
			if (ptr->pos < (ptr - 1)->pos)
			{
				memcpy(store, ptr, sizeof(MapType));
				memcpy(ptr, ptr - 1, sizeof(MapType));
				memcpy(ptr - 1, store, sizeof(MapType));
			}
		}
	}

/*
	// test the order
	ptr = map + 13;
	// to check the sorting along the position on sequence
	int tmp_pos = 0;

	while (ptr < end)
	{
		change_point = ptr;
		pos = ptr->pos;
		while ((ptr < end) && (ptr->pos == pos))
			++ptr;

		i = ptr - change_point;

		if (tmp_pos < pos)
			tmp_pos = pos;
		else
			std::cout << tmp_pos << "\t" << pos << "\t!!!" << std::endl;

//		std::cout << ptr - map - 1 - i + 1 << "\t" << ptr - map - 1 << "\t" << pos << "\t" << i << std::endl;
	}
*/
};
//--------------------------------------------------
MapType* Map::CopySorted( MapType* target, MapType* source, int max, int key )
{
	for( int i = 0; i < max; ++i )
	{
		if ( source->key == key )
		{
			memcpy( target, source, sizeof(MapType) );
			++target;
		}
		++source;
	}
	return target;
};
//--------------------------------------------------
void Map::OverlayMapByData( int* data )
{
	for ( int i = 0; i < size; ++i )
	{
		map[i].data = data[map[i].pos];
	}
};
//--------------------------------------------------
void Map::DeleteEmptyKeyRecords()
{
	int L = 0;
	int R = 0;
	int n = 0;

	while( R < size )
	{
		// find not empty
		while(( R < size )&&( !map[R].key ))
			++R;

		if ( R < size )
		{
			if ( L == R )
			{
				++L;
				++R;
				++n;
			}
			else
			{
				memcpy( map + L, map + R, sizeof( MapType ));
				map[R].key = 0;
				++R;
				++L;
				++n;
			}
		}
	}

	// reset map size
	size = n;
};
//--------------------------------------------------
void Map::PrintMap()
{
	if ( !map )
	{
		cout << "no map provided\n";
		return;
	}

	string name( "map.log" );
	ofstream log_file( name.c_str() );

	if ( !log_file.is_open() )
	{
		cout << PROG_NAME << " : error, opening file " << name << "\n"; 
		exit(1);
	}

	string label;
	int current_key;
	bool isIni;
	bool isTer;

	for( int i = 0; i < size; ++i )
	{
		current_key = map[i].key;

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

		log_file << i << " " << label << " " << map[i].key << " " << map[i].pos << " " << (int)map[i].fr << " " << (int)map[i].ph << " " << map[i].data << "\n";
	}

	log_file.close();
};
//--------------------------------------------------
void Map::PrintMapCount()
{
	vector<int> c( 256, 0);
	int i;

	for( i = 13; i < size - 13; ++i )
	{
		unsigned char x = (unsigned char)(map[i].key);
		++c[(int)x];
	}

	cout << "\n";
	cout << "dirDn "    << c[ dirDn ]    << "\n";
	cout << "dirAc "    << c[ dirAc ]    << "\n";
	cout << "dirStart " << c[ dirStart ] << "\n";
	cout << "dirEnd "   << c[ dirEnd ]   << "\n";
	cout << "revDn "    << c[ revDn ]    << "\n";
	cout << "revAc "    << c[ revAc ]    << "\n";
	cout << "revStart " << c[ revStart ] << "\n";
	cout << "revEnd "   << c[ revEnd ]   << "\n";
	cout << "all " << c[dirDn]+c[dirAc]+c[dirStart]+c[dirEnd]+c[revDn]+c[revAc]+c[revStart]+c[revEnd] << "\n";
	cout << "reserved " << reserved << "\n";
	cout << "\n";

	// check sorting
	for ( i = 1; i < size; ++i )
	{
		if ( map[i-1].pos > map[i].pos )
			cout << "error in map sorting at position "<< i << "\n";
	}
};
//--------------------------------------------------

