//**************************************************
// file: model.h
// project: ehmm3
//**************************************************

#ifndef MODEL_H__
#define MODEL_H__

#include <stdio.h>

#include "parser.h"
#include "site.h"
#include "length.h"
#include "markov.h"

//--------------------------------------------------
class Model : private Parser
{
public:
	Model()
	{
		organism = NULL; order = 0; phIni = NULL; phTer = NULL; phEx = NULL;

		initial_Intron_Dir_Prob = 0.33;
		initial_Intron_Rev_Prob = 0.33;
		initial_Intergenic_Prob = 0.33;
		initial_Coding_Dir_Prob = 0.005;
		initial_Coding_Rev_Prob = 0.005;

		BpModelOn = false;
		toBp = log0;
		aroundBp = 0;
		limitBpDur = 0;
		mSpacerOrder = 0;

		TAA_ON = true;
		TAG_ON = true;
		TGA_ON = true;

		GT_donor_ON = true;
		GC_donor_ON = false;
		toGTdonor = 1;
		toGCdonor = 0;

		norm_sites = false;
	}
	~Model() { delete [] organism; delete [] phIni; delete [] phTer; delete [] phEx; }

	char* organism;
	int order;

	Length lnIntergenic;
	Length lnIntron;
	Length lnSingle;
	Length lnInitial;
	Length lnExon;
	Length lnTerminal;

	Site don_0;
	Site don_1;
	Site don_2;
	Site acc_0;
	Site acc_1;
	Site acc_2;
	Site ini;
	Site terTAA;
	Site terTAG;
	Site terTGA;
	Site don_gc_0;
	Site don_gc_1;
	Site don_gc_2;

	bool BpModelOn;
	Site branchPoint;
	Length lnDonBp;
	Length lnBpAcc;
	double toBp;
	double aroundBp;
	int limitBpDur;
	MarkovIntr mSpacerBp;
	int mSpacerOrder;
	Site acc_bp_0;
	Site acc_bp_1;
	Site acc_bp_2;

	MarkovCod mCod;
	MarkovNon mNon;
	MarkovIntr mIntr;

	double* phIni;
	double* phTer;
	double* phEx;

	double toSingleGene;
	double toMultiGene;
	double toInternalExon;
	double toTerminalExon;

	double initial_Intron_Dir_Prob;
	double initial_Intron_Rev_Prob;
	double initial_Intergenic_Prob;
	double initial_Coding_Dir_Prob;
	double initial_Coding_Rev_Prob;

	bool Init( const char* filename );

	bool TAA_ON;
	bool TAG_ON;
	bool TGA_ON;

	bool norm_sites;

	bool GT_donor_ON;
	bool GC_donor_ON;
	double toGTdonor;
	double toGCdonor;

private:

	bool LoadValues();
	bool LoadMarkov();
	bool LoadLength( Length* ptr, const char* tag, int shift, int shift_partial );
	bool LoadSite( Site* ptr, const char* tag );
	bool LoadPhase();
	bool LoadInitial();
	void NormLog( double* ptr, int size );
	double Log( double x );
};
//--------------------------------------------------
#endif // MODEL_H__

