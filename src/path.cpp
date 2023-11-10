//**************************************************
// file: path.cpp
// project: ehmm3
// Alex Lomsadze
//**************************************************

#include <time.h>
#include <iostream>
#include <string>

#include "path.h"

//--------------------------------------------------
bool GenePredicted::AddExon( char* source, int size, char ffirst, char flast ) 
{
	if ( !nt )
		frFirst = ffirst;

	frLast = flast;

	if ( nt + size >= nt_reserved )
	{
		char *tmp = new char[ nt + size + 1000 ];
		if (!tmp) return false;
		nt_reserved = nt + size + 1000;

		if ( ntseq )
		{
			memcpy( tmp, ntseq, nt );
			delete ntseq;
		}
		ntseq = tmp;
	}

	memcpy( ntseq + nt, source, size );
	nt += size;

	exons_in_gene += 1;

	ntseq[nt] = '\0';

	return true;
};
//--------------------------------------------------
void GenePredicted::ReverseComplement()
{
	char tmp; 

	for( int i = 0; i < nt/2; ++i )
	{
		tmp = ntseq[i];
		ntseq[i] = ntseq[nt-i-1];
		ntseq[nt-i-1] = tmp;
	}

	for( int i = 0; i < nt; ++i )
	{
		if ( ntseq[i] < 4 )
			ntseq[i] = 3 - ntseq[i];
	}

	if ( strand == '+' )
		strand = '-';
	else
		strand = '+';

	int fr = frFirst;
	frFirst = frLast;
	frLast = fr;
};
//--------------------------------------------------
bool GenePredicted::Translate()
{
	aa_reserved = nt/3 + 1;

	delete aaseq;
	aaseq = new char[ aa_reserved ];
	if ( !aaseq ) return false;

	*aaseq = '\0';

	if ( strand == '-' )
		ReverseComplement();

	// skip partial codons, set reading frame
	//   X X X Y
	//   1 2 3 _
	//   fr 1 - start from beginning
	//   fr 2 - move by two positions upstream
	//   fr 3 - move by one position upstream

	int frame;

	if ( frFirst == 1 ) frame = 0;
	else if ( frFirst == 2 ) frame = 2;
	else if ( frFirst == 3 ) frame = 1;

	int i,j;
	aa = 0;

	extern char AA[65];

	for ( i = frame, j = 0; i < nt - 2 ; i += 3, ++j )
	{
		if ( (ntseq[i]>3)||(ntseq[i+1]>3)||(ntseq[i+2]>3) )
			aaseq[j] = 'X';
		else
			aaseq[j] = AA[ (ntseq[i]<<4) + (ntseq[i+1]<<2) + ntseq[i+2] ];
		++aa;
	}

	aaseq[j] = '\0';

	return true;
};
//--------------------------------------------------
void GenePredicted::PrintProtein( FILE* fp )
{
	Translate();

	int max = aa;
	if ( aaseq[ aa - 1 ] == '*' ) 
		--max;

	fprintf( fp, ">gene_%d|GeneMark.hmm|%d_aa\n", id, max );

	for( int i = 0; i < max; ++i )
	{
		if (( i != 0 )&&( i%60 == 0 ))
			fprintf( fp, "\n" );

		fprintf( fp, "%c", aaseq[i] );

		if ( aaseq[i] == '*' )
			printf( "Internal error! Please report to code developers!" );
	}

	fprintf( fp, "\n\n" );
};
//--------------------------------------------------
void GenePredicted::PrintNT( FILE* fp )
{
	if ( strand == '-' )
		ReverseComplement();

	int max = nt;

	fprintf( fp, ">gene_%d|GeneMark.hmm|%d_nt\n", id, max );

	for( int i = 0; i < max; ++i )
	{
		if (( i != 0 )&&( i%60 == 0 ))
			fprintf( fp, "\n" );

		fprintf( fp, "%c", LETTER_UP[ ntseq[i] ] );
	}

	fprintf( fp, "\n\n" );
};
//--------------------------------------------------
void PathRecord::SetRecord( Map* m, int L, int R )
{
	map_left = L;
	map_right = R;

	left  = m->map[L].pos;
	right = m->map[R].pos;
	keyLeft  = m->map[L].key;
	keyRight = m->map[R].key;
	phLeft  = m->map[L].ph;
	phRight = m->map[R].ph;

	if ( m->map[L].key & m->map[L].data )
		extr_L = '+';

	if (m->map[R].key & m->map[R].data )
		extr_R = '+';

	RescaleToSequence();
};
//--------------------------------------------------
inline void PathRecord::ShiftLeftRight( int shiftLeft, int shiftRight )
{
	// from marker positions shift to real sequence position

	if      ( !( keyLeft & mapIni )&& !( keyRight & mapTer ))
	{
		left  += shiftLeft;
		right += shiftRight;
	}
	else if (  ( keyLeft & mapIni )&& !( keyRight & mapTer ))
	{
		left  += 0;
		right += shiftRight;
	}
	else if ( !( keyLeft & mapIni )&&  ( keyRight & mapTer ))
	{
		left  += shiftLeft;
		right += 0;
	}
	else if (  ( keyLeft & mapIni )&&  ( keyRight & mapTer ))
	{
		left  += 0;
		right += 0;
	}
};
//--------------------------------------------------
inline void PathRecord::ChangePhaseToCodingFrameDirect()
{
	// in output FRAME of codon is defined as
	//   A T G X X X 
	//   1 2 3 1 2 3
	// intergenic phase is set to 0, not biological, just for code
	// intron phase has default meaning: 0, 1 and 2

	phLeft += 1;

	if ( !phRight )
		phRight = 3;
};
//--------------------------------------------------
inline void PathRecord::ChangePhaseToCodingFrameReverse()
{
	// as on direct strand
	//   X X X C A T
	//   3 2 1 3 2 1

	if ( !phLeft )
		phLeft = 3;

	phRight += 1;
};
//--------------------------------------------------
void PathRecord::RescaleToSequence()
{
	if (
		(( keyLeft & dirEnd   )&&( keyRight & revEnd   )) ||
		(( keyLeft & dirEnd   )&&( keyRight & dirStart )) ||
		(( keyLeft & revStart )&&( keyRight & dirStart )) ||
		(( keyLeft & revStart )&&( keyRight & revEnd   ))
		)
	{
		strcpy( state, _INTERGENIC );
		strand = '.';
		//        0 1 . . . n . . . . . . m . . 
		//  left        T A A . . . . . . A T G    right
		//                    |         |
		//                   n+1       m-1
		ShiftLeftRight( 1, -1 );
	}
	else if (( keyLeft & dirDn )&&( keyRight & dirAc ))
	{
		strcpy( state, _INTRON );
		strand = '+';
		//        0 1 . . . n . . . . . . m . . 
		//  left            G T . . . . A G . .    right
		//                  |             |
		//                  n             m
		ShiftLeftRight( 0, 0 );
	}
	else if	(( keyLeft & revAc )&&( keyRight & revDn ))
	{
		strcpy( state, _INTRON );
		strand = '-';
		//        0 1 . . . n . . . . . . m . . 
		//  left            C T . . . . A C . .    right
		//                  |             |
		//                  n             m
		ShiftLeftRight( 0, 0 );
	}
	else if (( keyLeft & dirStart )&&( keyRight & dirDn ))
	{
		strcpy( state, _INITIAL );
		strand = '+';
		ChangePhaseToCodingFrameDirect();
		//        0 1 . . . n . . . . . . m . . 
		//  left            A T G . . . . G T     right
		//                  |           |
		//                  n          m-1
		ShiftLeftRight( 0, -1 );
	}
	else if (( keyLeft & dirAc )&&(  keyRight & dirDn ))
	{
		strcpy( state, _INTERNAL );
		strand = '+';
		ChangePhaseToCodingFrameDirect();
		//        0 1 . . . n . . . . . . m . . 
		//  left          A G . . . . . . G T     right
		//                    |         |
		//                   n+1       m-1
		ShiftLeftRight( 1, -1 );
	}
	else if (( keyLeft & dirAc )&&( keyRight & dirEnd ))
	{
		strcpy( state, _TERMINAL );
		strand = '+';
		ChangePhaseToCodingFrameDirect();
		//        0 1 . . . n . . . . . . m . . 
		//  left          A G . . . . T A A      right
		//                    |           |
		//                   n+1          m
		ShiftLeftRight( 1, 0 );
	}
	else if (( keyLeft & dirStart )&&( keyRight & dirEnd ))
	{
		strcpy( state, _SINGLE );
		strand = '+';
		ChangePhaseToCodingFrameDirect();
		//        0 1 . . . n . . . . . . m . . 
		//  left            A T G   . T A A      right
		//                  |             |
		//                  n             m
		ShiftLeftRight( 0, 0 );
	}

	else if (( keyLeft & revDn )&&( keyRight & revStart ))
	{
		strcpy( state, _INITIAL );
		strand = '-';
		ChangePhaseToCodingFrameReverse();
		//        0 1 . . . n . . . . . . m . . 
		//  left          A C . . . . C A T      right
		//                    |           |
		//                   n+1          m
		ShiftLeftRight( 1, 0 );
	}
	else if (( keyLeft & revDn )&&( keyRight & revAc ))
	{
		strcpy( state, _INTERNAL );
		strand = '-';
		ChangePhaseToCodingFrameReverse();
		//        0 1 . . . n . . . . . . m . . 
		//  left          A C . . . . . . C T     right
		//                    |         |
		//                   n+1       m-1
		ShiftLeftRight( 1, -1 );
	}
	else if (( keyLeft & revEnd )&&( keyRight & revAc ))
	{
		strcpy( state, _TERMINAL );
		strand = '-';
		ChangePhaseToCodingFrameReverse();
		//        0 1 . . . n . . . . . . m . . 
		//  left            T T A . . . . C T     right
		//                  |           |
		//                  n          m-1
		ShiftLeftRight( 0, -1 );
	}
	else if (( keyLeft & revEnd )&&( keyRight & revStart ))
	{
		strcpy( state, _SINGLE );
		strand = '-';
		ChangePhaseToCodingFrameReverse();
		//        0 1 . . . n . . . . . . m . . 
		//  left            T T A . . C A T     right
		//                  |             |
		//                  n             m
		ShiftLeftRight( 0, 0 );
	}
	else
	{
		printf( "error in code\n" );
		exit(1);
	}
};
//--------------------------------------------------
bool Path::Init( int size, Model* mod  )
{
	if ( bestpath ) delete [] bestpath;
	bestpath = new BestPath [size];
	if ( !bestpath ) return false;

	// the initiation order of states is set in: void Map::SetIniTerm( int seqSize )
	//
	// 0        intergenic
	// 1 2 3    direct intron
	// 4 5 6    reverse intron
	// 7 8 9    direct coding
	// 10 11 12 reverse coding

	bestpath[0].prob = mod->initial_Intergenic_Prob;

	bestpath[1].prob = mod->initial_Intron_Dir_Prob;
	bestpath[2].prob = mod->initial_Intron_Dir_Prob;
	bestpath[3].prob = mod->initial_Intron_Dir_Prob;

	bestpath[4].prob = mod->initial_Intron_Rev_Prob;
	bestpath[5].prob = mod->initial_Intron_Rev_Prob;
	bestpath[6].prob = mod->initial_Intron_Rev_Prob;

	bestpath[7].prob = mod->initial_Coding_Dir_Prob;
	bestpath[8].prob = mod->initial_Coding_Dir_Prob;
	bestpath[9].prob = mod->initial_Coding_Dir_Prob;

	bestpath[10].prob = mod->initial_Coding_Rev_Prob;
	bestpath[11].prob = mod->initial_Coding_Rev_Prob;
	bestpath[12].prob = mod->initial_Coding_Rev_Prob;

	return true;
};
//--------------------------------------------------
bool Path::CountExonsAndGenes()
{
	if ( pathList.IsEmpty() )
		return false;

	PathRecord* data;

	// count genes and exons from first to last
	pathList.SetCurrentToFirst();

	int geneN = 0;
	int exonN = 0;
	bool is_gene = false;

	while ( !pathList.IsEnd() )
	{
		data = pathList.GetForward();

		if ( data->IsIntergenic() )
			is_gene = false;
		else 
		{
			if ( !is_gene )
			{
				geneN += 1;
				exonN = 0;
				is_gene = true;
			}

			if ( data->IsCoding() )
			{
				exonN += 1;
				data->exonN = exonN;
				data->geneN = geneN;
			}
			else if ( data->IsIntron() )
			{
				data->geneN = geneN;
				
			}
			else
				exit(1);
		}
	}

	// go from last to first and change exon number for negative strand
	pathList.SetCurrentToLast();

	exonN = 0;

	while ( !pathList.IsEnd() )
	{
		data = pathList.GetBackward();

		if ( data->IsIntergenic() )
			exonN = 0;
		else if ( data->IsCoding() && ( data->strand == '-' ))
		{
			if ( exonN )
				++exonN;
			else
				exonN = 1;

			data->exonN = exonN;
		}
	}

	return true;
};
//--------------------------------------------------
double Path::BestPathProbability( Map* map )
{
	double bestV = log0;

	for ( int i = map->size - 13; i < map->size; ++i )
	{
		if ( bestpath[i].prob > bestV )
			bestV = bestpath[i].prob;
	}

	printf( "# best path probability: %f\n", bestV );

	return bestV;
}
//--------------------------------------------------
bool Path::FillBestPathList( Map* map )
{
	// find which termination state (last 13 positions in the map array) has largest score
	// track back the best path to initiation states (first 13 position in the map array)
	// and save (label) states belonging to the best path

	double bestV = log0;
	int left = 0;
	int right = 0;

	for ( int i = map->size - 13; i < map->size; ++i )
	{
		if ( bestpath[i].prob > bestV )
		{
			bestV = bestpath[i].prob;
			right = i;
			left = bestpath[i].pos;
		}
	}

	// backtrack

 	PathRecord* data;

	while( right >= 13 )
	{
		data = new PathRecord;
		data->SetRecord( map, left, right );
		pathList.AddToHead( data );

		right = left;
		left = bestpath[right].pos;
	}

	CountExonsAndGenes();

	return true;
};
//--------------------------------------------------
bool Path::FillProteinList( char* seq )
{
	if ( pathList.IsEmpty() )
		return false;

	if ( !seq )
		return false;

	GenePredicted* prot = NULL;
 	PathRecord* data; 
	pathList.SetCurrentToFirst();

	while ( !pathList.IsEnd() )
	{
		data = pathList.GetForward();

		if ( data->IsIntergenic() && prot )
		{
				geneList.AddToTail( prot );
				prot = NULL;
		}
		else if ( data->IsCoding() )
		{
			if ( !prot )
				prot = new GenePredicted( data->geneN, data->strand );

			prot->AddExon( seq + data->left, data->right - data->left + 1, data->phLeft, data->phRight );

			if ( prot->left == -1 )
				prot->left = data->left;

			prot->right = data->right;

			if (data->extr_L == '+' || data->extr_R == '+')
				prot->supported = true;
		}
	}

	if ( prot )
		geneList.AddToTail( prot );

	return true;
};
//--------------------------------------------------
bool Path::Output( Map* map, CmdLine* cmd, Sequence* seq )
{
	FillBestPathList( map );
	FillProteinList( seq->seq );

	FILE* fp = fopen( cmd->outputFile.c_str(), "w" );
	if ( !fp )
	{
		printf( "can't open a file %s\n", cmd->outputFile.c_str() );
		return false;
	}

	string seqid = seq->seqid;
	int shift = 0;

	if (cmd->traceback)
	{
		seqid = seq->trace_seqid;
		shift = seq->trace_pos - 1;
	}

	// print header
	time_t ltime;
	time( &ltime );

	fprintf( fp, "Eukaryotic GeneMark.hmm version %s\n", PROG_VERSION );
	fprintf( fp, "Sequence name: %s\n",                  cmd->sequenceFile.c_str() );
	fprintf( fp, "FASTA defline: %s\n",                  seq->defline );
	fprintf( fp, "Sequence length: %d bp\n",             seq->size );
	fprintf( fp, "G+C content: %2.2f%%\n",               seq->gc );
	fprintf( fp, "Matrices file: %s\n",                  cmd->modelFile.c_str() );
	fprintf( fp, "%s\n\n",                               ctime( &ltime ) );
	fprintf( fp, "Predicted genes/exons\n\n" );
	fprintf( fp, "Gene Exon Strand Exon           Exon Range         Exon      Start/End\n" );
	fprintf( fp, "  #    #         Type                             Length       Frame\n\n" );

	// print predictions

 	PathRecord* data = 0;
	pathList.SetCurrentToFirst();

	GenePredicted* prot;
	geneList.SetCurrentToFirst();

	bool skip_gene = false;
	int last_gene_id = 0;

	while ( !pathList.IsEnd() )
	{
		data = pathList.GetForward();

		if ( data->IsIntergenic() )
		{
			fprintf( fp, "\n" );
		}
		else if ( data->IsCoding() )
		{
			if (last_gene_id != data->geneN)
			{
				prot = geneList.GetForward();

				if (data->geneN != prot->id) { printf("fatal error, check code here\n"); exit(2); }

				if (prot->nt < cmd->min_gene_length && !prot->supported)
					skip_gene = true;
				else
					skip_gene = false;

				last_gene_id = data->geneN;
			}

			if (skip_gene)
				continue;

			fprintf( fp, "%3d %3d   %c  %s  %10d %10d %10d %10d %d %c %c\n", 
				data->geneN, data->exonN, data->strand, data->state, data->left + 1 + shift, data->right + 1 + shift, data->right - data->left + 1, data->phLeft, data->phRight, data->extr_L, data->extr_R );
		}
	}

	if ( !data->IsIntergenic() )
		fprintf( fp, "\n" );

	// print proteins
	if ( cmd->protein )
	{
		geneList.SetCurrentToFirst();

		fprintf( fp, "# protein sequence of predicted genes\n\n" );
		while ( !geneList.IsEnd() )
		{
			prot = geneList.GetForward();
			prot->PrintProtein( fp );
		}
		fprintf( fp, "# end protein sequence\n\n" );
	}

	// print nt
	if ( cmd->nt )
	{
		geneList.SetCurrentToFirst();

		fprintf( fp, "# nucleotide sequence of predicted genes\n\n" );
		while ( !geneList.IsEnd() )
		{
			prot = geneList.GetForward();
			prot->PrintNT( fp );
		}
		fprintf( fp, "# end nucleotide sequence\n\n" );
	}

	fclose(fp);

	return true;
};
//--------------------------------------------------
bool Path::OutputTR( Map* map, CmdLine* cmd, Sequence* seq, Model* mod )
{
	FillBestPathList( map );
	FillProteinList( seq->seq );

	FILE* fp = fopen( cmd->outputFile.c_str(), "w" );
	if ( !fp )
	{
		printf( "error on open file: %s\n", cmd->outputFile.c_str() );
		return false;
	}

	// #-format
	fprintf( fp, "## TR formated data for training \n" );
	fprintf( fp, "## Eukaryotic GeneMark.hmm version %s\n", PROG_VERSION );
	fprintf( fp, "## mod_file: %s\n", cmd->modelFile.c_str() );
	fprintf( fp, "## seq_from_file: %s\n", cmd->sequenceFile.c_str() );
	fprintf( fp, "## defilne: %s\n", seq->defline );
	fprintf( fp, "## seq_length: %d bp\n", seq->size );
	fprintf( fp, "## seq_gc: %.0f\n", seq->gc );
	fprintf( fp, "## seq_id: %s\n", seq->seqid.c_str() );
	fprintf( fp, "## seq_L: %d\n", 1 );
	fprintf( fp, "## seq_R: %d\n", seq->size );

	//loop through predicted states
	PathRecord* data; 
	pathList.SetCurrentToFirst();

	geneList.SetCurrentToFirst();
	if ( geneList.IsEnd() )
	{
		fprintf( fp, "# seq_id_end\n" );
		fprintf( fp, "---\n" );
		fclose(fp);
		return true;
	}

	GenePredicted* prot = 0;

	fprintf( fp, "## gene_id\tcds_length\texons_in_gene\n");

	geneList.SetCurrentToFirst();
	while( ! geneList.IsEnd() )
	{
		prot = geneList.GetForward();
		fprintf( fp, "##\t%d\t%d\t%d\n", prot->id, prot->nt, prot->exons_in_gene );
	}
	geneList.SetCurrentToFirst();

	fprintf( fp, "## state\tgene_id\tgene_length\tleft\tright\tlength\tstrand\tph_left\tph_right\tseq\n" );

	int current_gene_id = -1;
	std::string str;

	while ( !pathList.IsEnd() )
	{
		data = pathList.GetForward();

		str.resize( data->length() );

		if ( data->strand == '-' )
		{
			for( std::size_t i = 0;  i < data->length(); ++i )
				str[i] = LETTER_UP_COMP[  seq->seq[ data->left + i ] ];

			std::reverse(str.begin(), str.end());
		}
		else
		{
			for( std::size_t i = 0;  i < data->length(); ++i )
				str[i] = LETTER_UP[  seq->seq[ data->left + i ] ];
		}

		if ( ! data->IsIntergenic() )
		{
			if ( current_gene_id != data->geneN )
			{
				current_gene_id = data->geneN;
				prot = geneList.GetForward();

				if ( prot->id != current_gene_id )
					exit(1);
			}
		}

		if ( data->IsCoding() )
		{
			// phase for coding regions was transfered to frame 

			if( data->strand == '+' )
			{
				if ( data->phLeft == 1 )
					;
				else if ( data->phLeft == 2 )
					str.insert( 0, "n" );
				else if ( data->phLeft == 3 )
					str.insert( 0, "nn" );
				else
					exit(1);

				if ( data->phRight == 1 )
					str.append( "nn" );
				else if ( data->phRight == 2 )
					str.append( "n" );
				else if ( data->phRight == 3 )
					;
				else
					exit(1);
			}
			else if ( data->strand == '-' )
			{
				if ( data->phRight == 1 )
					;
				else if ( data->phRight == 2 )
					str.insert( 0, "n" );
				else if ( data->phRight == 3 )
					str.insert( 0, "nn" );
				else
					exit(1);

				if ( data->phLeft == 1 )
					str.append( "nn" );
				else if ( data->phLeft == 2 )
					str.append( "n" );
				else if ( data->phLeft == 3 )
					;
				else
					exit(1);
			}
			else
				exit(1);

			if ( !strcmp( data->state, _INITIAL ) ||  !strcmp( data->state, _SINGLE ) )
			{
				std::string siteseq;
				siteseq.resize( mod->ini.p_width );

				if ( data->strand == '+' )
				{
					for( std::size_t i = 0;  i < mod->ini.p_width; ++i )
						siteseq[i] = LETTER_UP[  seq->seq[ data->left + i - mod->ini.p_margin ] ];

					fprintf( fp, "# %s\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%s\n", "START", prot->id, prot->nt, data->left + 1, data->right + 1, data->right - data->left + 1, data->strand, data->phLeft, data->phRight, siteseq.c_str() );
				}
				else if ( data->strand == '-' )
				{
					for( std::size_t i = 0;  i < mod->ini.p_width; ++i )
						siteseq[i] = LETTER_UP_COMP[  seq->seq[ data->right + i - mod->ini.p_width + mod->ini.p_margin + 1 ] ];

					std::reverse(siteseq.begin(), siteseq.end());

					fprintf( fp, "# %s\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%s\n", "START", prot->id, prot->nt, data->left + 1, data->right + 1, data->right - data->left + 1, data->strand, data->phLeft, data->phRight, siteseq.c_str() );
				}
				else
					exit(1);
			}
			
			if ( !strcmp( data->state, _TERMINAL ) || !strcmp( data->state, _SINGLE ) )
			{
				std::string siteseq;
				siteseq.resize( mod->terTAA.p_width );

				if ( data->strand == '+' )
				{
					for( std::size_t i = 0;  i < mod->terTAA.p_width; ++i )
						siteseq[i] = LETTER_UP[  seq->seq[ data->right + i - mod->terTAA.p_margin - 2 ] ];

					string STOP;
					if ( !siteseq.compare( mod->terTAA.p_margin, 3, "TAA") )
						STOP.assign( "STOP_TAA" );
					else if ( !siteseq.compare( mod->terTAG.p_margin, 3, "TAG") )
						STOP.assign( "STOP_TAG" );
					else if ( !siteseq.compare( mod->terTGA.p_margin, 3, "TGA") )
						STOP.assign( "STOP_TGA" );

					fprintf( fp, "# %s\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%s\n", STOP.c_str(), prot->id, prot->nt, data->left + 1, data->right + 1, data->right - data->left + 1, data->strand, data->phLeft, data->phRight, siteseq.c_str() );
				}
				else if ( data->strand == '-' )
				{
					for( std::size_t i = 0;  i < mod->terTAA.p_width; ++i )
						siteseq[i] = LETTER_UP_COMP[  seq->seq[ data->left + i - mod->terTAA.p_width + mod->terTAA.p_margin + 3] ];

					std::reverse(siteseq.begin(), siteseq.end());

					string STOP;
					if ( !siteseq.compare( mod->terTAA.p_margin, 3, "TAA") )
						STOP.assign( "STOP_TAA" );
					else if ( !siteseq.compare( mod->terTAG.p_margin, 3, "TAG") )
						STOP.assign( "STOP_TAG" );
					else if ( !siteseq.compare( mod->terTGA.p_margin, 3, "TGA") )
						STOP.assign( "STOP_TGA" );

					fprintf( fp, "# %s\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%s\n", STOP.c_str(), prot->id, prot->nt, data->left + 1, data->right + 1, data->right - data->left + 1, data->strand, data->phLeft, data->phRight, siteseq.c_str() );
				}
				else
					exit(1);
			}
		}


		if ( data->IsIntron() )
		{
			// exclude first and last - incomplete introns

			if ( data->left == 0 || data->right == seq->size - 1 )
			{
				;
			}
			else
			{
				std::string siteseq;
				siteseq.resize( mod->don_0.p_width );

				if ( data->strand == '+' )
				{
					for( std::size_t i = 0;  i < mod->don_0.p_width; ++i )
						siteseq[i] = LETTER_UP[  seq->seq[ data->left + i - mod->don_0.p_margin ] ];

					fprintf( fp, "# %s\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%s\n", "DON", prot->id, prot->nt, data->left + 1, data->right + 1, data->right - data->left + 1, data->strand, data->phLeft, data->phRight, siteseq.c_str() );
				}
				else if ( data->strand == '-' )
				{
					for( std::size_t i = 0;  i < mod->don_0.p_width; ++i )
						siteseq[i] = LETTER_UP_COMP[  seq->seq[ data->right + i - mod->don_0.p_width + mod->don_0.p_margin + 1 ] ];

					std::reverse(siteseq.begin(), siteseq.end());

					fprintf( fp, "# %s\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%s\n", "DON", prot->id, prot->nt, data->left + 1, data->right + 1, data->right - data->left + 1, data->strand, data->phLeft, data->phRight, siteseq.c_str() );
				}
				else
					exit(1);


				siteseq.resize( mod->acc_0.p_width );

				if ( data->strand == '+' )
				{
					for( std::size_t i = 0;  i < mod->acc_0.p_width; ++i )
						siteseq[i] = LETTER_UP[  seq->seq[ data->right + i - mod->acc_0.p_margin - 1 ] ];

					fprintf( fp, "# %s\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%s\n", "ACC", prot->id, prot->nt, data->left + 1, data->right + 1, data->right - data->left + 1, data->strand, data->phLeft, data->phRight, siteseq.c_str() );
				}
				else if ( data->strand == '-' )
				{
					for( std::size_t i = 0;  i < mod->acc_0.p_width; ++i )
						siteseq[i] = LETTER_UP_COMP[  seq->seq[ data->left + i - mod->acc_0.p_width + mod->acc_0.p_margin + 2 ] ];

					std::reverse(siteseq.begin(), siteseq.end());

					fprintf( fp, "# %s\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%s\n", "ACC", prot->id, prot->nt, data->left + 1, data->right + 1, data->right - data->left + 1, data->strand, data->phLeft, data->phRight, siteseq.c_str() );
				}
				else
					exit(1);
			}
		}

		for( int i = strlen( data->state ) -1; i >= 0; --i )
		{
			if ( isspace( data->state[i] ) )
			{
				 data->state[i] = '\0';
			}
		}

		fprintf( fp, "# %s\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%s\n", data->state, prot->id, prot->nt, data->left + 1, data->right + 1, data->right - data->left + 1, data->strand, data->phLeft, data->phRight, str.c_str() );
	}

	fprintf( fp, "# seq_id_end\n" );
	fprintf( fp, "---\n" );
	fclose(fp);

	return true;
};
//--------------------------------------------------
bool Path::OutputGFF3(Map* map, CmdLine* cmd, Sequence* seq)
{
	FillBestPathList(map);
	FillProteinList(seq->seq);

	FILE* fp = fopen(cmd->outputFile.c_str(), "w");
	if (!fp)
	{
		printf("can't open a file %s\n", cmd->outputFile.c_str());
		return false;
	}

	string seqid = seq->seqid;
	int shift = 0;

	if (cmd->traceback)
	{
		seqid = seq->trace_seqid;
		shift = seq->trace_pos - 1;
	}

	time_t ltime;
	time(&ltime);

	fprintf( fp, "##gff-version 3\n" );
	fprintf( fp, "# Eukaryotic GeneMark.hmm version %s\n", PROG_VERSION );
	fprintf( fp, "# Sequence name: %s\n",                  cmd->sequenceFile.c_str() );
	fprintf( fp, "# FASTA defline: %s\n",                  seq->defline );
	fprintf( fp, "# Sequence length: %d bp\n",             seq->size );
	fprintf( fp, "# G+C content: %2.2f%%\n",               seq->gc );
	fprintf( fp, "# Matrices file: %s\n",                  cmd->modelFile.c_str() );
	fprintf( fp, "# %s",                                   ctime( &ltime ) );
	fprintf( fp, "# FASTA definition line: %s\n",          seq->defline );
	fprintf( fp, "##sequence-region %s 1 %d\n",            seqid.c_str(), seq->size );

	// print predictions

 	PathRecord* data; 
	pathList.SetCurrentToFirst();

	GenePredicted* prot;
	geneList.SetCurrentToFirst();

	bool skip_gene = false;
	int last_gene_id = 0;
	
	char atr_gene[50] = "";
	char atr_mrna[50] = "";
	char atr_exon[50] = "";
	char atr_cds[50] = "";
	char atr_intron[50] = "";
	int gff3_ph;

	while ( !pathList.IsEnd() )
	{
		data = pathList.GetForward();

		if ( data->IsCoding() )
		{
			if ( last_gene_id != data->geneN )
			{
				prot = geneList.GetForward();

				if ( data->geneN != prot->id ) { printf("fatal error, check code here\n"); exit(2); }

				if (prot->nt < cmd->min_gene_length && !prot->supported)
					skip_gene = true;
				else
					skip_gene = false;

				last_gene_id = data->geneN;

				if (!skip_gene)
				{
					sprintf(atr_gene, "ID=gene%05d;Name=%d;", data->geneN, prot->id);
					sprintf(atr_mrna, "ID=mRNA%05d;Parent=gene%05d;Name=%d;", data->geneN, data->geneN, data->geneN);
					fprintf(fp, "%s\tGeneMark.hmm3\tgene\t%d\t%d\t.\t%c\t.\t%s\n", seqid.c_str(), prot->left + 1 + shift, prot->right + 1 + shift, data->strand, atr_gene);
					fprintf(fp, "%s\tGeneMark.hmm3\tmRNA\t%d\t%d\t.\t%c\t.\t%s\n", seqid.c_str(), prot->left + 1 + shift, prot->right + 1 + shift, data->strand, atr_mrna);
				}
			}

			if (skip_gene)
				continue;

			sprintf( atr_exon, "ID=exon%05d;Parent=mRNA%05d;Name=%d;", data->exonN, data->geneN, data->geneN );
			sprintf( atr_cds,  "ID=cds%05d;Parent=mRNA%05d;Name=%d;", 1, data->geneN, data->geneN  );

			fprintf( fp, "%s\tGeneMark.hmm3\texon\t%d\t%d\t.\t%c\t.\t%s\n", seqid.c_str(), data->left + 1 + shift, data->right + 1 + shift, data->strand, atr_exon );

			if ( data->strand == '+' )
				gff3_ph = (3 - data->phLeft + 1)%3;
			else
				gff3_ph = (3 - data->phRight + 1)%3;
					
			fprintf( fp, "%s\tGeneMark.hmm3\tCDS\t%d\t%d\t.\t%c\t%d\t%s\n", seqid.c_str(), data->left + 1 + shift, data->right + 1 + shift, data->strand, gff3_ph, atr_cds );
		}
		else if ( data->IsIntron() )
		{
			if (skip_gene)
				continue;

			if ( data->left < 4 || data->right > seq->size - 4 )
			{
				;
			}
			else
			{
				sprintf( atr_intron,  "ID=intron%05d;Parent=mRNA%05d;Name=%d;", 1, data->geneN, data->geneN  );
				fprintf( fp, "%s\tGeneMark.hmm3\tintron\t%d\t%d\t.\t%c\t%d\t%s\n", seqid.c_str(), data->left + 1 + shift, data->right + 1 + shift, data->strand, data->phLeft, atr_intron );
			}
		}
	}

	fclose(fp);

	return true;
};
//--------------------------------------------------
bool Path::OutputGTF(Map* map, CmdLine* cmd, Sequence* seq)
{
	FillBestPathList(map);
	FillProteinList(seq->seq);

	FILE* fp = fopen(cmd->outputFile.c_str(), "w");
	if (!fp)
	{
		printf("can't open a file %s\n", cmd->outputFile.c_str());
		return false;
	}

	string seqid = seq->seqid;
	int shift = 0;

	if (cmd->traceback)
	{
		seqid = seq->trace_seqid;
		shift = seq->trace_pos - 1;
	}

	time_t ltime;
	time(&ltime);

	fprintf( fp, "# Eukaryotic GeneMark.hmm version %s\n", PROG_VERSION );
	fprintf( fp, "# Sequence name: %s\n",                  cmd->sequenceFile.c_str() );
	fprintf( fp, "# FASTA defline: %s\n",                  seq->defline );
	fprintf( fp, "# Sequence length: %d bp\n",             seq->size );
	fprintf( fp, "# G+C content: %2.2f%%\n",               seq->gc );
	fprintf( fp, "# Matrices file: %s\n",                  cmd->modelFile.c_str() );
	fprintf( fp, "# %s",                                   ctime( &ltime ) );
	fprintf( fp, "# FASTA definition line: %s\n",          seq->defline );
	fprintf( fp, "# sequence-region %s 1 %d\n",            seqid.c_str(), seq->size );

	// print predictions

 	PathRecord* data; 
	pathList.SetCurrentToFirst();

	GenePredicted* prot = 0;
	geneList.SetCurrentToFirst();

	bool skip_gene = false;
	int last_gene_id = 0;

	char atr_gene[500] = "";
	char atr_cds[500] = "";
	int gff3_ph;

	while ( !pathList.IsEnd() )
	{
		data = pathList.GetForward();

		if ( data->IsCoding() )
		{
			if ( last_gene_id != data->geneN )
			{
				prot = geneList.GetForward();

				if ( data->geneN != prot->id ) { printf("fatal error, check code here\n"); exit(2); }

				if (prot->nt < cmd->min_gene_length && !prot->supported)
					skip_gene = true;
				else
					skip_gene = false;

				last_gene_id = data->geneN;

				if (!skip_gene)
				{
					sprintf(atr_gene, "gene_id \"%d_g\"; transcript_id \"%d_t\";", data->geneN, prot->id);
					fprintf(fp, "%s\tGeneMark.hmm3\tgene\t%d\t%d\t.\t%c\t.\t%s\n", seqid.c_str(), prot->left + 1 + shift, prot->right + 1 + shift, data->strand, atr_gene);
					fprintf(fp, "%s\tGeneMark.hmm3\tmRNA\t%d\t%d\t.\t%c\t.\t%s\n", seqid.c_str(), prot->left + 1 + shift, prot->right + 1 + shift, data->strand, atr_gene);
				}
			}

			if (skip_gene)
				continue;

			if ( data->strand == '+' )
				gff3_ph = (3 - data->phLeft + 1)%3;
			else
				gff3_ph = (3 - data->phRight + 1)%3;
			
			if (!strcmp(data->state, _INITIAL) || !strcmp(data->state, _SINGLE))
			{
				if (data->strand == '+')
				{
					if (data->extr_L == '+')
						fprintf(fp, "%s\tGeneMark.hmm3\tstart_codon\t%d\t%d\t.\t%c\t0\t%s %s\n", seqid.c_str(), data->left + 1 + shift, data->left + 3 + shift, data->strand, atr_gene, "evidence \"1_1\";");
					else
						fprintf(fp, "%s\tGeneMark.hmm3\tstart_codon\t%d\t%d\t.\t%c\t0\t%s\n", seqid.c_str(), data->left + 1 + shift, data->left + 3 + shift, data->strand, atr_gene);
				}
			}
			if (!strcmp(data->state, _TERMINAL) || !strcmp(data->state, _SINGLE))
			{
				if (data->strand == '-')
				{
					if (data->extr_L == '+')
						fprintf(fp, "%s\tGeneMark.hmm3\tstop_codon\t%d\t%d\t.\t%c\t0\t%s %s\n", seqid.c_str(), data->left + 1 + shift, data->left + 3 + shift, data->strand, atr_gene, "evidence \"1_1\";");
					else
						fprintf(fp, "%s\tGeneMark.hmm3\tstop_codon\t%d\t%d\t.\t%c\t0\t%s\n", seqid.c_str(), data->left + 1 + shift, data->left + 3 + shift, data->strand, atr_gene);
				}
			}

			string evidence("0_");

			if (data->extr_L == '+')
				evidence.assign("1_");
			
			if (data->extr_R == '+')
				evidence.append("1");
			else
				evidence.append("0");

			string state(data->state);
			while (*state.rbegin() == ' ')
			{
				state.erase(state.size() - 1);
			}

			if ( evidence.compare("0_0") )
				sprintf(atr_cds, " evidence \"%s\"; cds_type \"%s\"; count \"%d_%d\";", evidence.c_str(), state.c_str(), data->exonN, prot->exons_in_gene);
			else
				sprintf(atr_cds, " cds_type \"%s\"; count \"%d_%d\";", state.c_str(), data->exonN, prot->exons_in_gene);

			fprintf( fp, "%s\tGeneMark.hmm3\tCDS\t%d\t%d\t.\t%c\t%d\t%s%s\n", seqid.c_str(), data->left + 1 + shift, data->right + 1 + shift, data->strand, gff3_ph, atr_gene, atr_cds );

			if (!strcmp(data->state, _TERMINAL) || !strcmp(data->state, _SINGLE))
			{
				if (data->strand == '+')
				{
					if (data->extr_R == '+')
						fprintf(fp, "%s\tGeneMark.hmm3\tstop_codon\t%d\t%d\t.\t%c\t0\t%s %s\n", seqid.c_str(), data->right - 1 + shift, data->right + 1 + shift, data->strand, atr_gene, "evidence \"1_1\";");
					else
						fprintf(fp, "%s\tGeneMark.hmm3\tstop_codon\t%d\t%d\t.\t%c\t0\t%s\n", seqid.c_str(), data->right - 1 + shift, data->right + 1 + shift, data->strand, atr_gene);
				}
			}

			if (!strcmp(data->state, _INITIAL) || !strcmp(data->state, _SINGLE))
			{
				if (data->strand == '-')
				{
					if (data->extr_R == '+')
						fprintf(fp, "%s\tGeneMark.hmm3\tstart_codon\t%d\t%d\t.\t%c\t0\t%s %s\n", seqid.c_str(), data->right - 1 + shift, data->right + 1 + shift, data->strand, atr_gene, "evidence \"1_1\";");
					else
						fprintf(fp, "%s\tGeneMark.hmm3\tstart_codon\t%d\t%d\t.\t%c\t0\t%s\n", seqid.c_str(), data->right - 1 + shift, data->right + 1 + shift, data->strand, atr_gene);
				}
			}
		}
		else if ( data->IsIntron() )
		{
			if (skip_gene)
				continue;

			if ( data->left < 4 || data->right > seq->size - 4 )
			{
				;
			}
			else
			{
				if ((data->extr_L == data->extr_R) && (data->extr_L == '+'))
					fprintf(fp, "%s\tGeneMark.hmm3\tintron\t%d\t%d\t.\t%c\t%d\t%s evidence \"1_1\";\n", seqid.c_str(), data->left + 1 + shift, data->right + 1 + shift, data->strand, data->phLeft, atr_gene);
				else
					fprintf( fp, "%s\tGeneMark.hmm3\tintron\t%d\t%d\t.\t%c\t%d\t%s\n", seqid.c_str(), data->left + 1 + shift, data->right + 1 + shift, data->strand, data->phLeft, atr_gene );
			}
		}
	}

	fclose(fp);

	return true;
};
//--------------------------------------------------
bool Path::IntronsOut( Map* map, BaseIntron* intron, const char* name, bool append )
{
 	PathRecord* data; 
	pathList.SetCurrentToFirst();

	double No_BP = 0;
	double With_BP = 0;
	int L;
	int R;
	int j;
	int ph;

	int intron_count = 0;

	int best_bp = -1;
	char status;

	FILE* fp;

	if ( !append )
		fp = fopen( name, "w" );
	else
		fp = fopen( name, "a" );

	if ( !fp )
	{
		printf( "can't open a file %s\n", name );
		return false;
	}

	if ( !append )
		fprintf( fp, "status\t#\tgeneId\tstrand\tleft\tright\tlength\tbp_off\tbp_on\tspacer\n" );

	while ( !pathList.IsEnd() )
	{
		data = pathList.GetForward();

		if ( data->IsIntron() )
		{
			++intron_count;

			L = data->map_left;
			R = data->map_right;
			intron->Reset( R );
			ph = map->map[R].ph;

			if ( data->strand == '+' )
			{
				for ( j = R - 1; j >= L; --j )
				{
					if ( map->map[j].ph != ph )
						continue;
						
					if ( map->map[j].key == dirDn )
						No_BP = intron->Get(j);
					else if ( map->map[j].key == dirDn + mapIni )
						No_BP = intron->GetPartial(j);
				}

				if ( map->map[L].key == dirDn  )
					With_BP  = intron->GetDirWithBP( L, &best_bp );
				else if ( map->map[L].key == dirDn + mapIni )
					With_BP = intron->GetDirPartialWithBP(L);
				else
					printf( "error in code\n" );
				}

			if ( data->strand == '-' )
			{
				for ( j = R - 1; j >= L; --j )
				{
					if ( map->map[j].ph != ph )
						continue;

					if ( map->map[j].key == revAc )
						No_BP = intron->Get(j);
					else if ( map->map[j].key == revAc + mapIni )
						No_BP = intron->GetPartial(j);
				}

				if ( map->map[L].key == revAc )
					With_BP  = intron->GetRevWithBP( L, &best_bp );
				else if ( map->map[L].key == revAc + mapIni )
					With_BP = intron->GetRevPartialWithBP(L);
				else
					printf( "error in code\n" );
			}

			if ( With_BP > No_BP )
				status = 'U';
			else
				status = 'D';

			if ( !append )
			{
				fprintf( fp, "%c\t%d\t%d\t%c\t%d\t%d\t%d\t%.2f\t%.2f\t%d\n", status, intron_count, data->geneN, data->strand, data->left + 1, data->right + 1, data->right + 1 - data->left, No_BP, With_BP, best_bp + 1 );
				fprintf( fp, "%d\t%d\n", data->map_left, data->map_right );
			}
			else
			{
				                                                                           // prot->nt
				fprintf( fp, "# %s\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%c\n", "BP", data->geneN, 1000, data->left + 1, data->right + 1, data->right - data->left + 1, data->strand, data->phLeft, data->phRight, status );
			}
		}
	}

	fclose(fp);

	return true;
};
//--------------------------------------------------

