//**************************************************
// file: main.cpp
// project: ehmm3
// Alex Lomsadze
//**************************************************

#include "cmdline.h"
#include "sequence.h"
#include "map.h"
#include "model.h"
#include "precalc.h"
#include "baseviterbi.h"
#include "viterbi.h"
#include "bp.h"

char AA[65] = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
//             AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT
//             AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT
//             ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
//--------------------------------------------------

int main( int argc, char** argv )
{
	CmdLine cmdLine;
	if ( !cmdLine.Init(argc, argv) )
		return 1;

	Sequence sequence;
	if ( !sequence.LoadSequence( cmdLine.sequenceFile.c_str(), cmdLine.seqid  ))
		return 1;

#ifdef __DEBUG__
	sequence.PrintSeq();
#endif

	if ( !cmdLine.dataFile.empty() )
		if ( !sequence.LoadData( cmdLine.dataFile.c_str() ) )
			return 1;

#ifdef __DEBUG__
	sequence.PrintData();
#endif

	Model model;
	if ( !model.Init( cmdLine.modelFile.c_str() ) )
		return 1;

	Map map;
	if ( !map.Init( sequence.seq, sequence.data, sequence.size, &model ) )
		return 1;

#ifdef __DEBUG__
	map.PrintMap();
	map.PrintMapCount();
#endif

	if ( !sequence.CalculateAddress( model.order ) )
		return 1;

#ifdef __DEBUG__
	sequence.PrintAdr( model.order );
#endif

	BP bp;
	if ( model.BpModelOn )
	{
		if ( !bp.Init( &sequence, &model ) )
			return 1;

#ifdef __DEBUG__
		bp.Print();
#endif
	}

	Precalc pre;

	if ( sequence.data )
		pre.Init( model.order, &map, sequence.adr, cmdLine.maskp, sequence.data );
	else
		pre.Init( model.order, &map, sequence.adr, cmdLine.maskp );

	BaseIntergenic interg;
	if ( !interg.Init( map.size, map.map, model.order ) )
		return 1;
	pre.PreIntergenic( interg.data, &model.mNon );
	interg.len = &model.lnIntergenic;

	BaseIntron intron;
	if ( !intron.Init( map.size, map.map, model.order, sequence.seq ) )
		return 1;
	pre.PreIntron( intron.data, &model.mIntr );
	intron.len = &model.lnIntron;

	intron.TAA_ON = model.TAA_ON;
	intron.TAG_ON = model.TAG_ON;
	intron.TGA_ON = model.TGA_ON;

	BaseExon exon;
	if ( !exon.Init( map.size, map.map, model.order ) )
		return 1;
	pre.PreCoding( exon.data, &model.mCod );
	exon.lenIni = &model.lnInitial;
	exon.len = &model.lnExon;
	exon.lenTer = &model.lnTerminal;
	exon.lenSingle = &model.lnSingle;
	exon.phIni = model.phIni;
	exon.phTer = model.phTer;

	int k, m;
	for( k = 0; k < 3; ++k )
		for( m = 0; m < 3; ++m )
			exon.phEx[k][m] = model.phEx[ k*3 + m ];

	pre.PreSites( exon.site, &sequence, &model );

#ifdef __DEBUG__
	pre.PrintSites( exon.site );
#endif

	if ( model.BpModelOn )
	{
		intron.bp = &bp;
		if ( !intron.bp->IniAccAdj( map.size ) )
			return 1;
		pre.PreAccForBP( intron.bp->AccSiteAdj, &sequence, &model );
		intron.bp->AdjAcc( exon.site );
	}

	Viterbi vit;
	vit.map = &map;
	vit.mod = &model;
	if ( !vit.path.Init( map.size, &model ) )
		return 1;
	vit.exon = &exon;
	vit.intergenic = &interg;
	vit.intron = &intron;
	vit.verbose = cmdLine.verbose;

	vit.Run();

	extern char AA[65];
	if ( !model.TAA_ON )
		AA[ dirTAAcodon ] = 'x';
	if ( !model.TAG_ON )
		AA[ dirTAGcodon ] = 'x';
	if ( !model.TGA_ON )
		AA[ dirTGAcodon ] = 'x';

	if ( !cmdLine.outputFormat.compare( GFF3 ) )
		vit.path.OutputGFF3( &map, &cmdLine, &sequence );
	else if ( !cmdLine.outputFormat.compare( TR ) )
		vit.path.OutputTR( &map, &cmdLine, &sequence, &model );
	else if ( !cmdLine.outputFormat.compare( GTF ) )
		vit.path.OutputGTF( &map, &cmdLine, &sequence );
	else
		vit.path.Output( &map, &cmdLine, &sequence );

	if (( model.BpModelOn )&&( !cmdLine.bpintronFile.empty() ))
	{
		if ( !cmdLine.outputFormat.compare( TR ) )
		{
			vit.path.IntronsOut( &map, &intron, cmdLine.bpintronFile.c_str(), true );
		}
		else
		{
			vit.path.IntronsOut( &map, &intron, cmdLine.bpintronFile.c_str(), false );
		}
	}

	if ( cmdLine.report_best_prob )
		vit.path.BestPathProbability(&map);

	return 0;
};
//--------------------------------------------------

