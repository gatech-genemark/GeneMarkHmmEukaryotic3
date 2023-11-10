//**************************************************
// file: model.cpp
// project: ehmm3
//**************************************************

#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <string.h>

#include "model.h"
#include "common.h"

#include <iostream>

//--------------------------------------------------
void Model::NormLog( double* ptr, int size )
{
	int i;
	double sum = 0;
	for ( i = 0; i < size; ++i )
		sum += ptr[i];

	if (!sum)
		{printf( "error in normlog\n" ); exit(1); }
	else
		for ( i = 0; i < size; ++i )
		{
			if ( ptr[i] )
				ptr[i] = log( ptr[i]/sum ) + log((double)size);
			else
				ptr[i] = log0;
		}
};
//--------------------------------------------------
double Model::Log( double x )
{
	if ( x > 0 )
		return log(x);
	else if ( x == 0 )
		return log0;
	else
		{ printf( "error, negative argument in LOG function : %f\n", x ); exit(1); }
};
//--------------------------------------------------
bool Model::Init( const char* filename )
{
	if ( !InitParser( filename ) ) return false;

	if ( !LoadValues() ) {
		printf("error in model file %s\n", filename );
		return false;
	}
	return true;
};
//--------------------------------------------------
// Length class - duration of HMM state
// minimum state duration is 1 
// 
// 0 < min <= allowed length <= max
//
// If duration in file is more than "max", then such duration will be ignored and warning printed on stdout.
//    this is done during parameter file parsing
// All durations less than min are ignored without the warning.
//    this is done at PrepareLength() call
//
// Legacy issues:
//  a) memory allocation is done at in GetIntKeyArray() after reading of values from "_TAG",
//     as a minimum add such line to model
//     $..._TAG  1 0
//
//  duration is returned by:  Get(len) , where
//   len == ( C-style index for length of state + shift )  
//--------------------------------------------------
bool Model::LoadLength( Length* ptr, const char* tag, int shift, int shift_partial )
{
	bool x = false;
	char* str = NULL;

	CatStrings( &str, tag, "_MAX" );
	x = GetUniqValue( str, &ptr->max );
	if (!x) return false;

	CatStrings( &str, tag, "_MIN" );
	x = GetUniqValue( str, &ptr->min );
	if (!x) return false;

	if ( ptr->min > ptr->max )
	{
		std::cout << ptr->min << std::endl;
		std::cout << &ptr->max << std::endl;
		return false;
	}

	CatStrings( &str, tag, "_TAG" );
	x = GetUniqValue( str, &str );
	if (!x) return false;

	x = GetIntKeyArray( str, &ptr->data, ptr->max );
	if (!x) return false;

	CatStrings( &str, tag, "_TYPE" );
	x = GetUniqValue( str, &str );
	if (!x) return false;

	if ( !strcmp( str, "CONSTANT" ) )
	{
		ptr->FConstant( 1.0 );
	}
	else if ( !strcmp( str, "EXP" ) )
	{
		double exp_decay;

		CatStrings( &str, tag, "_EXP_D" );
		x = GetUniqValue( str, &exp_decay );
		if (!x) return false;

		ptr->Exp( exp_decay );
	} 
	else if ( !strcmp( str, "GAMMA" ) )
	{
		double gamma_a;
		double gamma_b;

		CatStrings( &str, tag, "_GAMMA_A" );
		x = GetUniqValue( str, &gamma_a );
		if (!x) return false;

		CatStrings( &str, tag, "_GAMMA_B" );
		x = GetUniqValue( str, &gamma_b );
		if (!x) return false;

		ptr->Gamma( gamma_a, gamma_b );
	} 

	ptr->PrepareLength( shift, shift_partial );

	delete [] str;

	return true;
};
//--------------------------------------------------
bool Model::LoadSite( Site* ptr, const char* tag )
{
	bool x = false;
	char* str = NULL;

	CatStrings( &str, tag, "_WIDTH" );
	x = GetUniqValue( str, &ptr->width );
	if (!x) return false;

	CatStrings( &str, tag, "_MARGIN" );
	x = GetUniqValue( str, &ptr->margin );
	if (!x) return false;

	CatStrings( &str, tag, "_ORDER" );
	x = GetUniqValue( str, &ptr->order );
	if (!x) return false;

	CatStrings( &str, tag, "_MAT" );
	x = GetStrKeyArray( str, &ptr->matDir, 4<<(2*ptr->order), ptr->width );
	if (!x) return false; 

	ptr->InitSite( tag );

	CatStrings( &str, tag, "_WIDTH_P" );
	if( IsKey( str ) )
	{
		x = GetUniqValue( str, &ptr->p_width );
		if (!x) return false;
	}

	CatStrings( &str, tag, "_MARGIN_P" );
	if( IsKey( str ) )
	{
		x = GetUniqValue( str, &ptr->p_margin );
		if (!x) return false;
	}

	delete [] str;

	return true;
};
//--------------------------------------------------
bool Model::LoadMarkov()
{
	bool x;
	double** tmp[5];

	tmp[0] = &mCod.d1a;
	tmp[1] = &mCod.d2a;
	tmp[2] = &mCod.d3a;
	tmp[3] = &mIntr.da;
	tmp[4] = &mNon.na;

	x = GetStrKeyArray( "MARKOV", tmp, 4<<(2*order), 5 );
	if (!x) return false;

	double probN;

	x = GetUniqValue( "N_PROB_NON",  &probN ); if (!x) return false;
	x = mIntr.Set( order, probN );             if (!x) return false;
	x = mNon.Set( order, probN );              if (!x) return false;

	x = GetUniqValue( "N_PROB_COD",  &probN ); if (!x) return false;
	x = mCod.Set( order, probN );              if (!x) return false;

	return true;
};
//--------------------------------------------------
bool Model::LoadPhase()
{
	bool x;

	x = GetIntKeyArray( "INITIAL_EXON_PHASE",  &phIni, 3 );    if (!x) return false;
	x = GetIntKeyArray( "TERMINAL_EXON_PHASE", &phTer, 3 );    if (!x) return false;
	x = GetIntKeyArray( "INTERNAL_EXON_PHASE", &phEx,  9 );    if (!x) return false;

	NormLog( phIni, 3 );
	NormLog( phTer, 3 );
	NormLog( phEx,  9 );

	return true;
};
//--------------------------------------------------
bool Model::LoadInitial()
{
	bool x = true;

	x = GetUniqValue( "INITIAL_DIRECT_INTRON",  &initial_Intron_Dir_Prob );
	x = GetUniqValue( "INITIAL_REVERSE_INTRON", &initial_Intron_Rev_Prob );
	x = GetUniqValue( "INITIAL_INTERGENIC",     &initial_Intergenic_Prob );
	x = GetUniqValue( "INITIAL_DIRECT_CODING",  &initial_Coding_Dir_Prob );
	x = GetUniqValue( "INITIAL_REVERSE_CODING", &initial_Coding_Rev_Prob );

	return x;
};
//--------------------------------------------------
bool Model::LoadValues()
{
	bool x = false;

	x = GetUniqValue( "ORGANISM",  &organism );         if (!x) return false;
	x = GetUniqValue( "ORDER",     &order );            if (!x) return false;

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
	//    x x x G T y y y y  y  A  G  x  x  x
	//    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
	//         L|                  |R            R-L=8     Lengt=9 as C-index=8      R-L= C-index + 0
	//
	//    x T A A y y y y y  y  y  A  T  G  x
	//    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
	//         L|                  |R            R-L=8     Lengt=7 as C-index=6      R-L= C-index + 2
	//

	x = LoadLength( &lnSingle,     "SINGLE",     0, 0 );   if (!x) return false;
	x = LoadLength( &lnExon,       "EXON",       2, 1 );   if (!x) return false;
	x = LoadLength( &lnInitial,    "INITIAL",    1, 0 );   if (!x) return false;
	x = LoadLength( &lnTerminal,   "TERMINAL",   1, 0 );   if (!x) return false;
	x = LoadLength( &lnIntron,     "INTRON",     0, 0 );   if (!x) return false;
	x = LoadLength( &lnIntergenic, "INTERGENIC", 2, 1 );   if (!x) return false;

#ifdef __DEBUG__
	lnIntergenic.Print( "intergenic_dur.log" );
	lnIntron.Print( "intron_dur.log" );
	lnExon.Print( "exon_dur.log" );
	lnInitial.Print( "initial_dur.log" );
	lnTerminal.Print( "terminal_dur.log" );
	lnSingle.Print( "single_dur.log" );
#endif

	if( IsKey( "NORM_SITES" ) )
	{
		x = GetUniqValue( "NORM_SITES", &norm_sites ); if (!x) return false;
	}

	x = LoadMarkov(); if (!x) return false;

	x = LoadSite( &ini, "INI" );                        if (!x) return false;
	x = LoadSite( &terTAA, "TERM_TAA" );                if (!x) return false;
	x = LoadSite( &terTAG, "TERM_TAG" );                if (!x) return false;
	x = LoadSite( &terTGA, "TERM_TGA" );                if (!x) return false;
	x = LoadSite( &don_0, "DONOR_0" );                  if (!x) return false;
	x = LoadSite( &don_1, "DONOR_1" );                  if (!x) return false;
	x = LoadSite( &don_2, "DONOR_2" );                  if (!x) return false;
	x = LoadSite( &acc_0, "ACCEPTOR_0" );               if (!x) return false;
	x = LoadSite( &acc_1, "ACCEPTOR_1" );               if (!x) return false;
	x = LoadSite( &acc_2, "ACCEPTOR_2" );               if (!x) return false;

	if ( IsKey( "BP" ))
	{
		x = GetUniqValue( "BP", &BpModelOn ); if (!x) return false;
	}

	if ( BpModelOn )
	{
		x = LoadSite( &branchPoint, "BRANCH" );               if (!x) return false;

		x = LoadLength( &lnBpAcc, "BP_ACC", 0, 0 );           if (!x) return false;
		x = LoadLength( &lnDonBp, "DON_BP", 0, 0 );           if (!x) return false;

#ifdef __DEBUG__
	lnBpAcc.Print( "bp_acc_dur.log" );
	lnDonBp.Print( "gt_bp_dur.log" );
#endif

		x = GetUniqValue( "TO_BP", &toBp );                   if (!x) return false;
		x = GetUniqValue( "AROUND_BP", &aroundBp );           if (!x) return false;

		x = GetUniqValue( "LIMIT_BP_DURATION", &limitBpDur ); if (!x) return false;

		x = GetUniqValue( "BP_SPACER_ORDER", &mSpacerOrder ); if (!x) return false;
		x = GetStrKeyArray( "MARKOV_BP_SPACER", &mSpacerBp.da, 4<<(2*mSpacerOrder), 1 );  if (!x) return false;

		double probN;
		x = GetUniqValue( "N_PROB_NON",  &probN );   if (!x) return false;
		x = mSpacerBp.Set( mSpacerOrder, probN );    if (!x) return false;

		x = LoadSite( &acc_bp_0, "ACC_BP_0" );       if (!x) return false;
		x = LoadSite( &acc_bp_1, "ACC_BP_1" );       if (!x) return false;
		x = LoadSite( &acc_bp_2, "ACC_BP_2" );       if (!x) return false;

		toBp = Log( toBp );
		aroundBp = Log( aroundBp );
	}

	x = LoadPhase();  if (!x) return false;

	x = GetUniqValue( "ToSingleGene",   &toSingleGene );   if (!x) return false;
	x = GetUniqValue( "ToMultiGene",    &toMultiGene );    if (!x) return false;
	x = GetUniqValue( "ToInternalExon", &toInternalExon ); if (!x) return false;
	x = GetUniqValue( "ToTerminalExon", &toTerminalExon ); if (!x) return false;

	toSingleGene   = Log( toSingleGene );
	toMultiGene    = Log( toMultiGene );
	toInternalExon = Log ( toInternalExon );
	toTerminalExon = Log ( toTerminalExon );

	if ( IsKey( "INITIAL_INTERGENIC" ))
		LoadInitial();

	{
		initial_Intron_Dir_Prob = Log( initial_Intron_Dir_Prob );
		initial_Intron_Rev_Prob = Log( initial_Intron_Rev_Prob );
		initial_Intergenic_Prob = Log( initial_Intergenic_Prob );
		initial_Coding_Dir_Prob = Log( initial_Coding_Dir_Prob );
		initial_Coding_Rev_Prob = Log( initial_Coding_Rev_Prob );
	}

	if ( IsKey( "TAA_ON" ) )
	{
		x = GetUniqValue( "TAA_ON", &TAA_ON ); if (!x) return false;
	}

	if ( IsKey( "TAG_ON" ) )
	{
		x = GetUniqValue( "TAG_ON", &TAG_ON ); if (!x) return false;
	}

	if ( IsKey( "TGA_ON" ) )
	{
		x = GetUniqValue( "TGA_ON", &TGA_ON ); if (!x) return false;
	}

	if (IsKey("GC_DONOR_ON"))
	{
		x = GetUniqValue("GC_DONOR_ON", &GC_donor_ON); if (!x) return false;
	}

	if (GC_donor_ON)
	{
		x = LoadSite(&don_gc_0, "DONOR_GC_0"); if (!x) return false;
		x = LoadSite(&don_gc_1, "DONOR_GC_1"); if (!x) return false;
		x = LoadSite(&don_gc_2, "DONOR_GC_2"); if (!x) return false;

		x = GetUniqValue("GT_DONOR_PROB", &toGTdonor); if (!x) return false;
		x = GetUniqValue("GC_DONOR_PROB", &toGCdonor); if (!x) return false;

		toGTdonor = Log(toGTdonor);
		toGCdonor = Log(toGCdonor);
	}
	else
	{
		toGTdonor = Log(toGTdonor);
	}

	if ( norm_sites )
	{
		acc_0.NormSite( mIntr.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 0, false, 2, order );
		acc_1.NormSite( mIntr.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 1, false, 2, order );
		acc_2.NormSite( mIntr.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 2, false, 2, order );

		don_0.NormSite( mIntr.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 0, true, 2, order );
		don_1.NormSite( mIntr.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 1, true, 2, order );
		don_2.NormSite( mIntr.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 2, true, 2, order );

		ini.NormSite( mNon.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 0, false, 3, order );

		terTAA.NormSite( mNon.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 0, true, 3, order );
		terTAG.NormSite( mNon.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 0, true, 3, order );
		terTGA.NormSite( mNon.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 0, true, 3, order );

		if (GC_donor_ON)
		{
			don_gc_0.NormSite(mIntr.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 0, true, 2, order);
			don_gc_1.NormSite(mIntr.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 1, true, 2, order);
			don_gc_2.NormSite(mIntr.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 2, true, 2, order);
		}

		if ( BpModelOn )
		{
			branchPoint.NormSite( mIntr.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 0, false, 0, order );

			acc_bp_0.NormSite( mIntr.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 0, false, 2, order );
			acc_bp_1.NormSite( mIntr.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 1, false, 2, order );
			acc_bp_2.NormSite( mIntr.ddd, mCod.ddd_1, mCod.ddd_2, mCod.ddd_3, 2, false, 2, order );
		}
	}

#ifdef __DEBUG__
	acc_0.Print( "acc_0_mat.log" );
	acc_1.Print( "acc_1_mat.log" );
	acc_2.Print( "acc_2_mat.log" );
	don_0.Print( "don_0_mat.log" );
	don_1.Print( "don_1_mat.log" );
	don_2.Print( "don_2_mat.log" );
	ini.Print( "ini_mat.log" );
	terTAA.Print( "taa_mat.log" );
	terTGA.Print( "tga_mat.log" );
	terTAG.Print( "tag_mat.log" );

	if (GC_donor_ON)
	{
		don_gc_0.Print("don_gc_0_mat.log");
		don_gc_1.Print("don_gc_1_mat.log");
		don_gc_2.Print("don_gc_2_mat.log");
	}

	if ( BpModelOn )
	{
		acc_bp_0.Print( "acc_bp_0_mat.log" );
		acc_bp_1.Print( "acc_bp_1_mat.log" );
		acc_bp_2.Print( "acc_bp_2_mat.log" );
		branchPoint.Print( "bp_mat.log" );
	}
#endif

	return x;
};
//--------------------------------------------------

