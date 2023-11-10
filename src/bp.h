//**************************************************
// file: bp.h
// project: ehmm3
//**************************************************

#ifndef BP_H__
#define BP_H__

#include "model.h"
#include "sequence.h"

//--------------------------------------------------
class BP
{
public:

	BP() { dirBpValues = NULL; revBpValues = NULL; maxBpAcc = 0; minBpAcc = 0; dirSpacer = NULL; revSpacer = NULL; index = NULL; AccSiteAdj = NULL; mapSize = 0; }; 
	~BP() { delete dirBpValues; delete revBpValues; delete dirSpacer; delete revSpacer; delete AccSiteAdj; }

	bool Init( Sequence* in_seq, Model* in_mod );

	double DirIntron( int L, int R, double prob, int map_pos );
	double DirIntron( int L, int R, double prob, int map_pos, int* bp_pos_best );
	double DirIntronPartial( int L, int R, double prob, int map_pos );
	double RevIntron( int L, int R, double prob, int map_pos );
	double RevIntron( int L, int R, double prob, int map_pos, int* bp_pos_best );
	double RevIntronPartial( int L, int R, double prob, int map_pos );

	bool Print();

	bool IniAccAdj( int map_size );
	void AdjAcc( double* acc );
	double* AccSiteAdj;

private:
	
	int size;
	Model* mod;
	Sequence* seq;

	int* index;

	double* dirBpValues;
	double* revBpValues;

	double* dirSpacer;
	double* revSpacer;

	int mapSize;

	bool AllocArray( double** data, int size );
	void CalculateSites();
	void CalculateChain();

	int maxBpAcc;
	int minBpAcc;
};
//--------------------------------------------------
#endif // BP_H__

