//**************************************************
// file: baseviterbi.h
// project: ehmm3
//**************************************************

#ifndef BASEVITERBI_H__
#define BASEVITERBI_H__

#include <memory.h>

#include "model.h"
#include "length.h"
#include "bp.h"
#include "common.h"

//--------------------------------------------------
class BaseIntergenic
{
public:

	BaseIntergenic() { data = 0; buf = 0; };
	~BaseIntergenic() { delete [] data; };

	bool Init( int size, MapType* in_map, int in_order )
	{
		if ( data ) delete [] data;
		data = new double [ size ];
		if ( !data ) return false;

		memset( data, 0, size* sizeof(double) );

		map = in_map;
		order = in_order;

		return true;
	};

	double* data;

	Length* len;

	void Reset( int i )
	{
		pos = i;
		ini_buf = data[i];
		buf = 0;
	};

	//    x T A A y y y y y  y  y  A  T  G  x
	//    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
	//         L|                  |R            R-L=8     Lengt=7 as C-index=6      R-L= C-index + 2
	//    ( R - L - 2 = C-index ) + 1 = length
	//      R - L - 1 = length

	bool Stop( int j )
	{
		if ( map[pos].pos - map[j].pos - 1 > len->max )
			return true;
		else
			return false;
	};

	double Get( int j )
	{
		int n = map[pos].pos - map[j].pos;
		double x = log0;

		// McCh over full intergenic state
		if ( n - 1 > order )
		{
			if ( buf == 0 )
				buf = ini_buf;
			else
				buf += data[j];

			x = buf + len->Get(n);
		}
		
		return x;
	};

	double GetPartial( int j )
	{
		int n = map[pos].pos - map[j].pos;
		double x = log0;

		// McCh over partial intergenic state with only one incomplete boundary. thus +1
		if ( n > order )
		{
			if ( buf == 0 )
				buf = ini_buf;
			else
				buf += data[j];

			x = buf + len->GetPartial(n);
		}
		
		return x;
	};

private:

	int pos;
	double buf;
	double ini_buf;

	MapType* map;
	int order;
};
//--------------------------------------------------
class BaseIntron
{
public:

	BaseIntron() { data = 0; buf = 0; TAA_ON = true; TAG_ON = true; TGA_ON = true; };
	~BaseIntron() { delete [] data; };

	bool Init( int size, MapType* in_map, int in_order, char* in_seq )
	{
		if ( data ) delete [] data;
		data = new double [ size ];
		if ( !data ) return false;

		memset( data, 0, size* sizeof(double) );

		map = in_map;
		order = in_order;
		seq = in_seq;

		return true;
	};

	double* data;

	char* seq;
	Length* len;

	BP* bp;

	void Reset( int i )
	{
		pos = i;
		ini_buf = data[i];
		buf = 0;
	};

	//    x x x G T y y y y  y  A  G  x  x  x
	//    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
	//         L|                  |R            R-L=8     Lengt=9 as C-index=8      R-L= C-index + 0
	//    ( R - L - 0 = C-index ) + 1 = length
	//      R - L + 1 = length

	bool Stop( int j )
	{
		if ( map[pos].pos - map[j].pos + 1 > len->max )
			return true;
		else
			return false;
	};

	bool CheckSplicedCodonDir( int ph, int j )
	{
		int codon;

		if ( ph == 0 )
		{;}
		else if (( ph == 1 )&&( map[j].pos >= 1 ))
		{
			codon = *(seq + map[j].pos - 1);
			codon <<= 2;
			codon += *(seq + map[pos].pos + 1);
			codon <<= 2;
			codon += *(seq + map[pos].pos + 2);

			if ( TAA_ON && (codon==dirTAAcodon) ) return false;
			if ( TAG_ON && (codon==dirTAGcodon) ) return false;
			if ( TGA_ON && (codon==dirTGAcodon) ) return false;
		}
		else if (( ph == 2 )&&( map[j].pos >= 2 ))
		{
			codon =  *(seq + map[j].pos - 2);
			codon <<= 2;
			codon +=  *(seq + map[j].pos - 1);
			codon <<= 2;
			codon += *(seq + map[pos].pos + 1);

			if ( TAA_ON && (codon==dirTAAcodon) ) return false;
			if ( TAG_ON && (codon==dirTAGcodon) ) return false;
			if ( TGA_ON && (codon==dirTGAcodon) ) return false;
		}

		return true;
	};

	bool CheckSplicedCodonRev( int ph, int j )
	{
		int codon;

		if ( ph == 0 )
		{;}
		else if (( ph == 1 )&&( map[j].pos >= 2 )) 
		{
			codon = *(seq + map[j].pos - 2);
			codon <<= 2;
			codon += *(seq + map[j].pos - 1);
			codon <<= 2;
			codon += *(seq + map[pos].pos + 1);

			if ( TAA_ON && (codon==revTAAcodon) ) return false;
			if ( TAG_ON && (codon==revTAGcodon) ) return false;
			if ( TGA_ON && (codon==revTGAcodon) ) return false;
		}
		else if (( ph == 2 )&&( map[j].pos >= 1 ))
		{
			codon =  *(seq + map[j].pos - 1);
			codon <<= 2;
			codon +=  *(seq + map[pos].pos + 1);
			codon <<= 2;
			codon += *(seq + map[pos].pos + 2);

			if ( TAA_ON && (codon==revTAAcodon) ) return false;
			if ( TAG_ON && (codon==revTAGcodon) ) return false;
			if ( TGA_ON && (codon==revTGAcodon) ) return false;
		}

		return true;
	};

	// McCh over full intron state excluding GT and AG:   R - L + 1 - 4

	double Get( int j )
	{
		int n = map[pos].pos - map[j].pos;
		double x = log0;

		if ( n - 3 > order )
		{
			if ( buf == 0 )
				buf = ini_buf;
			else
				buf += data[j];

			x = buf + len->Get(n);
		}

		return x;
	};

	double GetDirWithBP( int j )
	{
		double x = log0;

		if ( map[pos].pos - map[j].pos - 3 > order )
		{
			// call BP class to get probability
			// input: j, pos, buf
			// output: best probability with BP
			// best position for BP is lost

			x = bp->DirIntron( map[j].pos, map[pos].pos, buf, pos );
		}
		
		return x;
	};

	double GetDirWithBP( int j, int* bp_pos_best )
	{
		double x = log0;

		if ( map[pos].pos - map[j].pos - 3 > order )
		{
			x = bp->DirIntron( map[j].pos, map[pos].pos, buf, pos, bp_pos_best );
		}
		
		return x;
	};

	double GetDirPartialWithBP( int j )
	{
		double x = log0;

		if ( map[pos].pos - map[j].pos - 1 > order )
		{
			x = bp->DirIntronPartial( map[j].pos, map[pos].pos, buf, pos );
		}

		return x;
	};

	double GetRevWithBP( int j, int* bp_pos_best )
	{
		double x = log0;

		if ( map[pos].pos - map[j].pos - 3 > order )
		{
			x = bp->RevIntron( map[j].pos, map[pos].pos, buf, j, bp_pos_best );
		}
		
		return x;
	};

	double GetRevWithBP( int j )
	{
		double x = log0;

		if ( map[pos].pos - map[j].pos - 3 > order )
		{
			x = bp->RevIntron( map[j].pos, map[pos].pos, buf, j );
		}
		
		return x;
	};

	double GetRevPartialWithBP( int j )
	{
		double x = log0;

		if ( map[pos].pos - map[j].pos - 1 > order )
		{
			x = bp->RevIntronPartial( map[j].pos, map[pos].pos, buf, j );
		}

		return x;
	};

	// McCh over partial intron state one of the boundaries: GT or AG:   R - L + 1 - 2

	double GetPartial( int j )
	{
		int n = map[pos].pos - map[j].pos;
		double x= log0;

		if ( n - 1 > order )
		{
			if ( buf == 0 )
				buf = ini_buf;
			else
				buf += data[j];

			x = buf + len->GetPartial(n);
		}
		
		return x;
	};

	bool TAA_ON;
	bool TAG_ON;
	bool TGA_ON;

private:

	int pos;
	double buf;
	double ini_buf;

	MapType* map;
	int order;
};
//--------------------------------------------------
class BaseExon
{
public:

	BaseExon() { data = 0; site = 0; };
	~BaseExon() { delete [] data;  delete [] site; };

	bool Init( int size, MapType* in_map, int in_order )
	{

		if ( data ) delete [] data;
		data = new double [ size ];
		if ( !data ) return false;

		if ( site ) delete [] site;
		site = new double [ size ];
		if ( !site ) return false;

		memset( site, 0, size* sizeof(double) );
		memset( data, 0, size* sizeof(double) );

		map = in_map;
		order = in_order;

		return true;
	};

	double* data;
	double* site;

	Length* lenIni;
	Length* len;
	Length* lenTer;
	Length* lenSingle;

	double* phIni;
	double* phTer;
	double phEx[3][3];

	void Reset( int i, char fr, char ph )
	{
		pos = i;
		ini_buf = data[i] + site[i];
		buf = 0;
		right_ph = ph;
	};

	//    x x x A T G y y y  T  A  A  x  x  x    Single
	//    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
	//         L|                  |R            R-L=8     Lengt=9 as C-index=8      R-L= C-index + 0
	//
	//    x x A G y y y y y  y  y  G  T  x  x    Internal exon
	//    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
	//         L|                  |R            R-L=8     Lengt=7 as C-index=6      R-L= C-index + 2
	//
	//    x x x A T G y y y  y  y  G  T  x  x    Initial exon
	//    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
	//         L|                  |R            R-L=8     Lengt=8 as C-index=7      R-L= C-index + 1
	//
	//    x x A G y y y y y  T  A  A  x  x  x    Terminal exon
	//    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
	//         L|                  |R            R-L=8     Lengt=8 as C-index=7      R-L= C-index + 1
	//

	//    ( R - L - * = C-index ) + 1 = length
	//    shortest corresponds to Single exon: 0 :   R - L + 1 = length

	bool Stop( int j )
	{
		if ( map[pos].pos - map[j].pos + 1 > len->max )
			return true;
		else
			return false;
	};

	// McCh over full internal exon state:   R - L - 2 + 1

	double Get( int j, char ph )
	{
		int n = map[pos].pos - map[j].pos;
		double x = log0;

		if ( n - 1 > order )
		{
			if ( buf == 0 )
				buf = ini_buf;
			else
				buf += data[j];

			x = buf + len->Get(n) + site[j] + phEx[ph][right_ph];
		}

		return x;
	};

	// McCh over full initial exon state excluding start codon :  R - L - 1 + 1 - 3

	double GetInitial( int j, char ph )
	{
		int n = map[pos].pos - map[j].pos;
		double x = log0;

		if ( n - 3 > order )
		{
			if ( buf == 0 )
				buf = ini_buf;
			else
				buf += data[j];

			x = buf + lenIni->Get(n) + site[j] + phIni[ph];
		}
		else if ( n >= 3 )
		{
			x = lenIni->Get(n) + site[j] + phIni[ph] + (n-3)*xp + site[pos];
		}
		
		return x;
	};

	// McCh over full terminal exon state excluding stop codon :  R - L - 1 + 1 - 3

	double GetTerminal( int j, char ph )
	{
		int n = map[pos].pos - map[j].pos;
		double x = log0;

		if ( n - 3 > order )
		{
			if ( buf == 0 )
				buf = ini_buf;
			else
				buf += data[j];

			x = buf + lenTer->Get(n) + site[j] + phTer[ph];
		}
		else if ( n >= 3 )
		{
			x =  lenTer->Get(n) + site[j] + phTer[ph] + (n-3)*xp + site[pos];
		}
		
		return x;
	};

	double GetPartial( int j )
	{
		return log0;
	};

	// McCh over full single exon state excluding stop and start codons :  R - L - 0 + 1 - 6

	double GetSingle( int j )
	{
		int n = map[pos].pos - map[j].pos;
		double x = log0;

		if ( n - 5 > order )
		{
			if ( buf == 0 )
				buf = ini_buf;
			else
				buf += data[j];

			x = buf + lenSingle->Get(n) + site[j];
		}
		
		return x;
	};

private:

	int pos;
	double buf;
	double ini_buf;
	char right_ph;

	MapType* map;
	int order;
};
//--------------------------------------------------
#endif  // BASEVITERBI_H__

