//**************************************************
// file: viterbi.h
// project: ehmm3
//**************************************************

#ifndef VITERBI_H__
#define VITERBI_H__

#include "model.h"
#include "map.h"
#include "path.h"
#include "baseviterbi.h"

//--------------------------------------------------
class Viterbi
{
public:

	Map* map;
	Model* mod;

	MapType* m;

	Path path;

	BaseExon* exon;
	BaseIntron* intron;
	BaseIntergenic* intergenic;

	void Run();

	bool verbose;

private:

	void Exon_dirDn( int i );
	void Exon_dirEnd( int i );
	void Exon_revAc( int i );
	void Exon_revStart( int i );
	void Intron_dirAc( int i );
	void Intron_revDn( int i );
	void Intergenic_uni( int i );

	void Exon_dirTerminal( int i );
	void Exon_revTerminal( int i );
	void Intron_dirTerminal( int i );
	void Intron_revTerminal( int i );
	void Intergenic_Terminal( int i );
};
//--------------------------------------------------
#endif  // VITERBI_H__

