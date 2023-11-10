//**************************************************
// file: path.h
// project: ehmm3
// Alex Lomsadze
//**************************************************

#ifndef PATH_H__
#define PATH_H__

#include <stdio.h>
#include <string.h>

#include "common.h"
#include "alist.h"
#include "map.h"
#include "model.h"
#include "cmdline.h"
#include "sequence.h"
#include "baseviterbi.h"

#include <string>
#include <vector>
#include <algorithm>

//--------------------------------------------------
class BestPath 
{
public:
	BestPath()
	{
		pos = pos_long = 0;
		prob = prob_long = log0;
	}

	int pos;
	int pos_long;

	double prob;
	double prob_long;
};
//--------------------------------------------------
class PathRecord
{
public:
	PathRecord()
	{
		left = right = 0;
		keyLeft = keyRight = 0;
		phLeft = phRight = 0;
		strand = '\0';
		state[0] = '\0';
		geneN = 0;
		exonN = 0;
		extr_L = '-';
		extr_R = '-';

		map_left = 0;
		map_right = 0;
	}

	int left;                    // position of HMM state on sequence; initilize from "map" using key position and then shift to true location in SetRecord
	int right;                   // position of HMM state on sequence; initilize from "map" using key position and then shift to true location in SetRecord

	int keyLeft;                 // key value from "map"
	int keyRight;                // key value from "map" 

	char phLeft;                 // phase calue from "map"; change to frane in SetRecord
	char phRight;                // phase value from "map"; change to frane in SetRecord

	char strand;                 // set in SetRecord
	char state[_LABEL_SIZE];     // set in SetRecord

	int geneN;                   // is set after full list is created 
	int exonN;                   // is set after full list is created

	char extr_L;                 // "+" if left border has PLUS support; "-" otherwise; initialize to "-";
	char extr_R;                 // "+" if right border has PLUS support; "-" otherwise; initialize to "-"; 

	int map_left;                // index of array - HMM state starts here - "map" array; is set ones - no modification
	int map_right;               // index of array - HMM state stops here - "map" array; is set ones - no modification

	int length(void) { return right - left + 1; }

	void SetRecord( Map* m, int L, int R );

	bool IsCoding()
	{
		if ( strcmp( state, _INTERGENIC ) &&  strcmp( state, _INTRON ) ) return true;
		else return false;
	}
	bool IsIntergenic()
	{
		if ( !strcmp( state, _INTERGENIC ) ) return true;
		else return false;
	}
	bool IsIntron()
	{
		if ( !strcmp( state, _INTRON ) ) return true;
		else return false;
	}

private:
	inline void ShiftLeftRight( int shiftLeft, int shiftRight );
	inline void ChangePhaseToCodingFrameDirect();
	inline void ChangePhaseToCodingFrameReverse();
	void RescaleToSequence();
};
//--------------------------------------------------
class GenePredicted
{
public:
	GenePredicted( int in_id, char in_strand )
	{
		nt_reserved = aa_reserved = 0; 
		nt = aa = 0;
		frFirst = frLast = 0;
		ntseq = aaseq = NULL; 
		id = in_id;
		strand = in_strand; 
		left = -1;
		right = -1;
		exons_in_gene = 0;
		supported = false;
	}
	~GenePredicted() { delete [] ntseq; delete [] aaseq; }

	bool AddExon( char* source, int size, char ffirst, char flast ) ;
	void PrintProtein( FILE* fp );
	void PrintNT( FILE* fp );

	int left;  // ?
	int right;  // ?
	int id;

	int nt;  // nt seq length:  nt <= nt_reserved
	int aa;  // aa prot length:  aa <= aa_reserved

	char* ntseq;
	char* aaseq;

	int exons_in_gene;
	bool supported;  // true if at least one element of a gene is supported

private:

	int nt_reserved;
	int aa_reserved;

	char frFirst;  // start frame of ntseq
	char frLast;  // end frame of ntseq 
		
	char strand;

	void ReverseComplement();
	bool Translate();
};
//--------------------------------------------------
class Path
{
public:

	Path() { bestpath = NULL; }
	~Path() { delete [] bestpath; }

	bool Init( int size, Model* mod );

	inline void Reset( int i )
	{
		current = i;
		max = log0;
	}

	inline void Try( int j, double x )
	{
		x += bestpath[j].prob;	
		if ( x > max )
		{
			max = x;
			bestpath[current].prob = max;
			bestpath[current].pos = j;
		}
	}

	bool Output( Map* map, CmdLine* cmd, Sequence* seq );
	bool OutputGFF3( Map* map, CmdLine* cmd, Sequence* seq );
	bool OutputGTF( Map* map, CmdLine* cmd, Sequence* seq );
	bool OutputTR( Map* map, CmdLine* cmd, Sequence* seq, Model* mod );

	bool IntronsOut( Map* map, BaseIntron* intron, const char* name, bool append );
	double BestPathProbability( Map* map );

private:

	int current;
	double max;

	BestPath* bestpath;

	List<PathRecord> pathList;
	List<GenePredicted> geneList;

	bool CountExonsAndGenes();
	bool FillProteinList( char* seq );
	bool FillBestPathList( Map* map );
};
//--------------------------------------------------
#endif  // PATH_H__

