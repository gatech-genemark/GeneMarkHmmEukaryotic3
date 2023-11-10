//**************************************************
// file: bp.cpp
// project: ehmm3
//**************************************************

#include <memory.h>

#include "bp.h"

//--------------------------------------------------
bool BP::Init( Sequence* in_seq, Model* in_mod )
{
	seq = in_seq;
	mod = in_mod;

	size = seq->size;

	if ( !AllocArray( &dirBpValues, size ) )
		return false;

	if ( !AllocArray( &revBpValues, size ) )
		return false;

	CalculateSites();

	if ( !AllocArray( &dirSpacer, size ) )
		return false;

	if ( !AllocArray( &revSpacer, size ) )
		return false;

	if ( mod->order == mod->mSpacerOrder )
		index = seq->adr;
	else
	{
		if ( !seq->SetAddressBP( mod->mSpacerOrder ) )
			return false;
		index = seq->adr_bp;
	}

	CalculateChain();

	// move to C-style index
	//   BP| 1  2  A  G |
	//       | min = 4
	//   BP| 1  2  3  4  5  6  A  G |
	//       | max = 8
	
	//  C++ style
	//  i = min    and  i < max
	maxBpAcc = mod->lnBpAcc.max;
	minBpAcc = mod->lnBpAcc.min - 1;

	return true;
};
//--------------------------------------------------
bool BP::AllocArray( double** data, int size )
{
	if ( *data ) delete [] *data;
	*data = new double [ size ];
	if ( !*data ) {
		printf( "error, out of memory\n" );
		return false;
	}
	
	for( int i = 0; i < size; ++i )
		(*data)[i] = log0;

	return true;
};
//--------------------------------------------------
//   |-|---GT---------------------|1|2|3|4|5|6|7|8|9|-----AG-----
//    i=0                                          i
//--------------------------------------------------
void BP::CalculateSites()
{
	int i;
	
	// set pointer to last letter in motif
	// set shift to compensate site::GetValue

	for (i = mod->branchPoint.width - 1; i < size; ++i )
		dirBpValues[i] = mod->branchPoint.GetValueDirNew( seq->seq + i, -1 );

	for (i = 0; i < size - mod->branchPoint.width + 1; ++i )
		revBpValues[i] = mod->branchPoint.GetValueRevNew( seq->seq + i, 0 );
};
//--------------------------------------------------
bool BP::Print()
{
	FILE* fp;
	fp = fopen("bp.log", "w" );
	if ( !fp )
	{
		printf( "can't open a file %s\n", "bp.log" );
		return false;
	}

	int i;

	for (i = 0; i < size; ++i )
	{
		if ( dirBpValues[i] >= log0/2 )
			fprintf( fp, "%d\td\t%f\n", i , dirBpValues[i] );
		else
			fprintf( fp, "%d\td\t%c\n", i , 'N' );
	}

	for (i = 0; i < size; ++i )
	{
		if ( revBpValues[i] >= log0/2 )
			fprintf( fp, "%d\tr\t%f\n", i , revBpValues[i] );
		else
			fprintf( fp, "%d\tr\t%c\n", i , 'N' );
	}

	fclose(fp);

	fp = fopen("bp_chain.log", "w" );
	if ( !fp )
	{
		printf( "can't open a file %s\n", "bp_chain.log" );
		return false;
	}

	for (i = 0; i < size; ++i )
	{
		if ( dirSpacer[i] >= log0/2 )
			fprintf( fp, "%d\td\t%f\n", i , dirSpacer[i] );
		else
			fprintf( fp, "%d\td\t%c\n", i , 'N' );
	}

	for (i = 0; i < size; ++i )
	{
		if ( revSpacer[i] >= log0/2 )
			fprintf( fp, "%d\tr\t%f\n", i , revSpacer[i] );
		else
			fprintf( fp, "%d\tr\t%c\n", i , 'N' );
	}

	fclose(fp);

	return true;
};
//--------------------------------------------------
//   |-|---GT---------------------|1|2|3|4|5|6|7|8|9|-----AG-----
//    i=0                                              i
//  log[P(Xi|spaccer)] - log[P(Xi|intron)]    - cumulative
//--------------------------------------------------
void BP::CalculateChain()
{
	// loop counter along seq;
	int i;
	// log probability for Spacer
	double x = 0;
	// log probaility for Intron
	double y = 0;
	// difference of logs
	double z = x - y;

	// end of sequence for Chain calculations
	// differes from size by largest order
	// we are not computing values at left end
	int end = size;

	if ( mod->order >= mod->mSpacerOrder )
		end -= mod->order;
	else
		end -= mod->mSpacerOrder;

	// tmp for cumulative calculation
	double prev = 0;
	// direct strand
	for( i = 0; i < end; ++i )
	{
		// get direct spacer prob for position i
		x = mod->mSpacerBp.d[ index[i] ];
		// get direct intron prob for position i
		y = mod->mIntr.d[ seq->adr[i] ];
		// save ratio
		z = x - y ;
		prev += z;
		dirSpacer[i] = prev;
	}

	prev = 0;
	// reverse strand
	for( i = 0; i < end; ++i )
	{
		// get reverce spacer prob for position i
		x = mod->mSpacerBp.r[ index[i] ];
		// get reverce intron prob for position i
		y = mod->mIntr.r[ seq->adr[i] ];
		// save ratio
		z = x - y;
		prev += z;
		revSpacer[i] = prev;
	}
};
//--------------------------------------------------
bool BP::IniAccAdj( int map_size )
{
	if ( !AllocArray( &AccSiteAdj, map_size ) )
		return false;

	mapSize = map_size;

	return true;
};
//--------------------------------------------------
void BP::AdjAcc( double* acc )
{
	for( int i = 0; i < mapSize; ++i )
	{
		if (( AccSiteAdj[i] > log0/1000 )&&( acc[i] > log0/1000 ))
			AccSiteAdj[i] -= acc[i];
		else
			AccSiteAdj[i] = log0;
	}
};
//--------------------------------------------------
double BP::DirIntron( int L, int R, double prob, int map_pos )
{
	// for intron L-R  find probability with BP

	double x = 0;          // current probability
	double maX = log0;     // best probability - return this value
	
	int orderI = mod->order;                 // mkch order of intron region
	int orderS = mod->mSpacerOrder;          // mkch order of BP-GT spacer

	int acc_shift = mod->acc_bp_0.margin;    // and of intron is part of ACC and part of spacer - 

	//  check that intron is long enough to hold BP model
	if ( mod->lnDonBp.min + mod->branchPoint.width + mod->lnBpAcc.min  > R-L+1)
		return log0;

	// BP should be inside intron
	int max_spacer_length = maxBpAcc;
	if ( max_spacer_length > R-L+1  - mod->branchPoint.width - mod->lnDonBp.min )
		max_spacer_length = R-L+1  - mod->branchPoint.width - mod->lnDonBp.min;

	for (int i = minBpAcc; i < max_spacer_length; ++i )
	{
		// get Branch Point probability
		x = dirBpValues[ R - i - 1 ];

		// get Duration from BP to Acceptor
		x += mod->lnBpAcc.Get( i );

		// get Duration from Donor to BP
		//                     intron length   BP-GT-length             bp-motif        move to index
		x += mod->lnDonBp.Get(  (R - L + 1)  -  (i + 1)     -   mod->branchPoint.width   -   1 );
		
		// get McCh probability for intron GT/AG excluded
		x += prob;
		// substitute McCh intron by Spacer
		x += dirSpacer[R - 2 - acc_shift] - dirSpacer[R - i - 1 ];
		// add difference between short and long acceptor model
		x += AccSiteAdj[ map_pos ];

		if ( x > maX )
			maX = x;
	}

	return maX;
};
//--------------------------------------------------
double BP::DirIntron( int L, int R, double prob, int map_pos, int* bp_pos_best )
{
	double x = 0;
	double maX = log0;

	*bp_pos_best = -1;

	int acc_shift = mod->acc_bp_0.margin;

	if ( mod->lnDonBp.min + mod->branchPoint.width + mod->lnBpAcc.min > R-L+1)
		return log0;

	// BP should be inside intron
	int max_spacer_length = maxBpAcc;
	if ( max_spacer_length > R-L+1  - mod->branchPoint.width - mod->lnDonBp.min )
		max_spacer_length = R-L+1  - mod->branchPoint.width - mod->lnDonBp.min;

	for (int i = minBpAcc; i < max_spacer_length; ++i )
	{
		// get Branch Point probability
		x = dirBpValues[ R - i - 1 ];
		// get Duration from BP to Acceptor
		x += mod->lnBpAcc.Get( i );
		// get Duration from Donor to BP
		x += mod->lnDonBp.Get( R - L + 1 - i - 1 - mod->branchPoint.width - 1 );
		// get McCh probability for intron GT/AG excluded
		x += prob;
		// substitute McCh intron by Spacer
		x += dirSpacer[R - 2 - acc_shift] - dirSpacer[R - i - 1];
		// add difference between short and long acceptor model
		x += AccSiteAdj[ map_pos ];

		if ( x > maX )
		{
			maX = x;
			*bp_pos_best = i;
		}
	}

	return maX;	
};
//--------------------------------------------------
double BP::DirIntronPartial( int L, int R, double prob, int map_pos )
{
	return log0;
};
//--------------------------------------------------
double BP::RevIntronPartial( int L, int R, double prob, int map_pos )
{
	return log0;
};
//--------------------------------------------------
double BP::RevIntron( int L, int R, double prob, int map_pos )
{
	double x = 0;
	double maX = log0;

	int acc_shift = mod->acc_bp_0.margin;

	if ( mod->lnDonBp.min + mod->branchPoint.width + mod->lnBpAcc.min  > R-L+1)
		return log0;

	// BP should be inside intron
	int max_spacer_length = maxBpAcc;
	if ( max_spacer_length > R-L+1  - mod->branchPoint.width - mod->lnDonBp.min )
		max_spacer_length = R-L+1  - mod->branchPoint.width - mod->lnDonBp.min;

	for (int i = minBpAcc; i < max_spacer_length; ++i )
	{
		if ( L + i + 1 > size )
			break;

		// get Branch Point probability
		x = revBpValues[ L + i + 1 ];
		// get Duration from BP to Acceptor
		x += mod->lnBpAcc.Get( i );
		// get Duration from Donor to BP
		x += mod->lnDonBp.Get( R - L + 1 - i - 1 - mod->branchPoint.width - 1 );
		// get McCh probability for intron GT/AG excluded
		x += prob;
		// substitute McCh intron by Spacer
		x += revSpacer[L + i  ] - revSpacer[ L + acc_shift + 1 ];
		// add difference between short and long acceptor model
		x += AccSiteAdj[ map_pos ];

		if ( x > maX )
			maX = x;
	}

	return maX;
};
//--------------------------------------------------
double BP::RevIntron( int L, int R, double prob, int map_pos, int* bp_pos_best )
{
	double x = 0;
	double maX = log0;
	*bp_pos_best = -1;

	int acc_shift = mod->acc_bp_0.margin;

	if ( mod->lnDonBp.min + mod->branchPoint.width + mod->lnBpAcc.min  > R-L+1)
		return log0;

	// BP should be inside intron
	int max_spacer_length = maxBpAcc;
	if ( max_spacer_length > R-L+1  - mod->branchPoint.width - mod->lnDonBp.min )
		max_spacer_length = R-L+1  - mod->branchPoint.width - mod->lnDonBp.min;

	for (int i = minBpAcc; i < max_spacer_length; ++i )
	{
		if ( L + i + 1 > size )
			break;
		// get Branch Point probability
		x = revBpValues[ L + i + 1 ];
		// get Duration from BP to Acceptor
		x += mod->lnBpAcc.Get( i );
		// get Duration from Donor to BP
		x += mod->lnDonBp.Get( R - L + 1 - i - 1 - mod->branchPoint.width - 1 );
		// get McCh probability for intron GT/AG excluded
		x += prob;
		// substitute McCh intron by Spacer
		x += revSpacer[L + i ] - revSpacer[ L + acc_shift + 1 ];
		// add difference between short and long acceptor model
		x += AccSiteAdj[ map_pos ];

		if ( x > maX )
		{
			maX = x;
			*bp_pos_best = i;
		}
	}

	return maX;
};
//--------------------------------------------------

