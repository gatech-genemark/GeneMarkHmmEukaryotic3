//**************************************************
// file: viterbi.cpp
// project: ehmm3
//**************************************************

#include <stdio.h>

#include "viterbi.h"
#include "common.h"

//--------------------------------------------------
void Viterbi::Run()
{
	m = map->map;

//	map->PrintMap();

	int v;

	if (verbose)
		printf("map_points total: %d\n", map->size);

	for ( int i = 13; i < map->size; ++i )
	{
		if (verbose && ((i % 10000) == 0))
		{
			printf("working: %d %.1f\n", i, 100.0*i/ map->size);
		}

		v = m[i].key;

		if ( v == dirDn )
						Exon_dirDn(i);
		else if ( v == dirAc )
						Intron_dirAc(i);
		else if ( v == revDn )
						Intron_revDn(i);
		else if ( v == revAc )
						Exon_revAc(i);
		else if ( v == dirEnd )
						Exon_dirEnd(i);
		else if ( v == revEnd )
						Intergenic_uni(i);
		else if ( v == dirStart )
						Intergenic_uni(i);
		else if ( v == revStart )
						Exon_revStart(i);
		else if ( v == mapTer + revEnd )
						Intergenic_Terminal(i);
		else if ( v == mapTer + dirAc )
						Intron_dirTerminal(i);
		else if ( v == mapTer + revDn )
						Intron_revTerminal(i);
		else if ( v == mapTer + dirDn ) 
						Exon_dirTerminal(i);
		else if ( v == mapTer + revAc )
						Exon_revTerminal(i);
		else
			printf( "error in map, position %d\n", map->map[i].pos );
	}
};
//--------------------------------------------------
void Viterbi::Exon_dirDn( int i )
{
	double x;
	int fr = m[i].fr;

	exon->Reset( i, fr, m[i].ph );
	path.Reset( i );

	for ( int j = i - 1; j >= 0; --j )
	{
		if ( m[j].fr == fr )
		{
			if ( exon->Stop(j) )
				break;

			if ( m[j].key == dirAc )
			{
				x = exon->Get( j, m[j].ph ) + mod->toInternalExon;
				path.Try( j, x );
			}
			else if ( m[j].key == dirEnd )
			{
				break;
			}
			else if ( m[j].key == dirStart )
			{
				x = exon->GetInitial( j, m[j].ph ) + mod->toMultiGene;
				path.Try( j, x );
			}
			else if ( m[j].key == dirAc + mapIni )
			{
				x = exon->GetPartial(j);
				path.Try( j, x );
			}
		}
	}
};
//--------------------------------------------------
void Viterbi::Intron_dirAc( int i )
{
	double x;
	int ph = m[i].ph;

	intron->Reset( i );
	path.Reset( i );

	for ( int j = i - 1; j >= 0; --j )
	{
		if ( m[j].ph == ph )
		{
			if ( intron->Stop(j) )
				break;

			if ( m[j].key == dirDn )
			{
				x = intron->Get(j);

				if ( intron->CheckSplicedCodonDir( ph, j ) )
				{
					path.Try( j, x );

					if ( mod->BpModelOn )
					{
						// must be called after Get
						x = intron->GetDirWithBP(j);
						path.Try( j, x );
					}
				}
			}
			else if ( m[j].key == dirDn + mapIni )
			{
				x = intron->GetPartial(j);
				path.Try( j, x );

				if ( mod->BpModelOn )
				{
					// must be called after Get
					x = intron->GetDirPartialWithBP(j);
					path.Try( j, x );
				}
			}
		}
	}
};
//--------------------------------------------------
void Viterbi::Intron_revDn( int i )
{
	double x;
	int ph = m[i].ph;

	intron->Reset( i );
	path.Reset( i );

	for ( int j = i - 1; j >= 0; --j )
	{
		if ( m[j].ph == ph )
		{
			if ( intron->Stop(j) )
				break;

			if ( m[j].key == revAc )
			{
				x = intron->Get(j);

				if ( intron->CheckSplicedCodonRev( ph, j ) )
				{
					path.Try( j, x );

					if ( mod->BpModelOn )
					{
						// must be called after Get
						x = intron->GetRevWithBP(j);
						path.Try( j, x );
					}
				}
			}
			else if ( m[j].key == revAc + mapIni )
			{
				x = intron->GetPartial(j);
				path.Try( j, x );

				if ( mod->BpModelOn )
				{
					// must be called after Get
					x = intron->GetRevPartialWithBP(j);
					path.Try( j, x );
				}
			}
		}
	}
};
//--------------------------------------------------
void Viterbi::Exon_revAc( int i )
{
	double x;
	int fr = m[i].fr;

	exon->Reset( i, fr, m[i].ph );
	path.Reset( i );

	for ( int j = i - 1; j >= 0; --j )
	{
		if ( m[j].fr == fr )
		{
			if ( exon->Stop(j) )
				break;

			if ( m[j].key == revDn )
			{
				x = exon->Get( j, m[j].ph ) + mod->toInternalExon;
				path.Try( j, x );
			}
			else if ( m[j].key == revEnd )
			{
				x = exon->GetTerminal( j, m[j].ph ) + mod->toMultiGene;
				path.Try( j, x );
				break;
			}
			else if ( m[j].key == revDn + mapIni )
			{
				x = exon->GetPartial(j);
				path.Try( j, x );
			}
		}
	}
};
//--------------------------------------------------
void Viterbi::Exon_dirEnd( int i )
{
	double x;
	int fr = m[i].fr;

	exon->Reset( i, fr, m[i].ph );
	path.Reset( i );

	for ( int j = i - 1; j >= 0; --j )
	{
		if ( m[j].fr == fr )
		{
			if ( exon->Stop(j) )
				break;

			if ( m[j].key == dirAc )
			{
				x = exon->GetTerminal( j, m[j].ph ) + mod->toTerminalExon;
				path.Try( j, x );
			}
			else if ( m[j].key == dirEnd )
			{
				break;
			}
			else if ( m[j].key == dirStart )
			{
				x = exon->GetSingle(j) + mod->toSingleGene;
				path.Try( j, x );
			}
			else if ( m[j].key == dirAc + mapIni )
			{
				x = exon->GetPartial(j);
				path.Try( j, x );
			}
		}
	}
};
//--------------------------------------------------
void Viterbi::Intergenic_uni( int i )
{
	double x;

	intergenic->Reset( i );
	path.Reset( i );

	for ( int j = i - 1; j >= 0; --j )
	{
		if ( intergenic->Stop(j) )
			break;

		if ( ( m[j].key == dirEnd )||( m[j].key == revStart ) )
		{
			x = intergenic->Get(j);
			path.Try( j, x );
		}
		else if ( m[j].key == dirEnd + mapIni )
		{
#ifdef ETP__
			x = intergenic->GetPartial(j) + 100;
#else
			x = intergenic->GetPartial(j);
#endif
			path.Try( j, x );
		}
	}
};
//--------------------------------------------------
void Viterbi::Exon_revStart( int i )
{
	double x;
	int fr = m[i].fr;

	exon->Reset( i, fr, m[i].ph );
	path.Reset( i );

	for ( int j = i - 1; j >= 0; --j )
	{
		if ( m[j].fr == fr )
		{
			if ( exon->Stop(j) )
				break;

			if ( m[j].key == revDn )
			{
				x = exon->GetInitial( j, m[j].ph ) + mod->toTerminalExon;
				path.Try( j, x );
			}
			else if ( m[j].key == revEnd )
			{
				x = exon->GetSingle(j) + mod->toSingleGene;
				path.Try( j, x );
				break;
			}
			else if ( m[j].key == revDn + mapIni )
			{
				x = exon->GetPartial(j);
				path.Try( j, x );
			}
		}
	}
};
//--------------------------------------------------
void Viterbi::Intergenic_Terminal( int i )
{
	double x;

	intergenic->Reset( i );
	path.Reset( i );

	for ( int j = i - 1; j >= 0; --j )
	{
		if ( intergenic->Stop(j) )
			break;

		if ( ( m[j].key == dirEnd )||( m[j].key == revStart ) )
		{
#ifdef ETP__
			x = intergenic->GetPartial(j) + 100;
#else
			x = intergenic->GetPartial(j);
#endif
			path.Try( j, x );
		}
		else if ( m[j].key == dirEnd + mapIni )
		{
			x = intergenic->GetPartial(j);
			path.Try( j, x );
		}
	}
};
//--------------------------------------------------
void Viterbi::Intron_dirTerminal( int i )
{
	double x;
	int ph = m[i].ph;

	intron->Reset( i );
	path.Reset( i );

	for ( int j = i - 1; j >= 0; --j )
	{
		if ( m[j].ph == ph )
		{
			if ( intron->Stop(j) )
				break;

			if ( m[j].key == dirDn )
			{
				x = intron->GetPartial(j);
				path.Try( j, x );
			}
			else if ( m[j].key == dirDn + mapIni )
			{
				x = intron->GetPartial(j);
				path.Try( j, x );
			}
		}
	}
};
//--------------------------------------------------
void Viterbi::Intron_revTerminal( int i )
{
	double x;
	int ph = m[i].ph;

	intron->Reset( i );
	path.Reset( i );

	for ( int j = i - 1; j >= 0; --j )
	{
		if ( m[j].ph == ph )
		{
			if ( intron->Stop(j) )
				break;

			if ( m[j].key == revAc )
			{
				x = intron->GetPartial(j);
				path.Try( j, x );
			}
			else if ( m[j].key == revAc + mapIni )
			{
				x = intron->GetPartial(j);
				path.Try( j, x );
			}
		}
	}
};
//--------------------------------------------------
void Viterbi::Exon_dirTerminal( int i )
{
	double x;
	int fr = m[i].fr;

	exon->Reset( i, fr, m[i].ph );
	path.Reset( i );

	for ( int j = i - 1; j >= 0; --j )
	{
		if ( m[j].fr == fr )
		{
			if ( exon->Stop(j) )
				break;

			if ( m[j].key == dirAc )
			{
				x = exon->GetPartial(j);
				path.Try( j, x );
			}
			else if ( m[j].key == dirEnd )
			{
				break;
			}
			else if ( m[j].key == dirStart )
			{
				x = exon->GetPartial(j);
				path.Try( j, x );
			}
			else if ( m[j].key == dirAc + mapIni )
			{
				x = exon->GetPartial(j);
				path.Try( j, x );
			}
		}
	}
};
//--------------------------------------------------
void Viterbi::Exon_revTerminal( int i )
{
	double x;
	int fr = m[i].fr;

	exon->Reset( i, fr, m[i].ph );
	path.Reset( i );

	for ( int j = i - 1; j >= 0; --j )
	{
		if ( m[j].fr == fr )
		{
			if ( exon->Stop(j) )
				break;

			if ( m[j].key == revDn )
			{
				x = exon->GetPartial(j);
				path.Try( j, x );
			}
			else if ( m[j].key == revEnd )
			{
				x = exon->GetPartial(j);
				path.Try( j, x );
				break;
			}
			else if ( m[j].key == revDn + mapIni )
			{
				x = exon->GetPartial(j);
				path.Try( j, x );
			}
		}
	}
};
//--------------------------------------------------

