//**************************************************
// file: sequence.h
// project: ehmm3
// Alex Lomsadze
//**************************************************

#ifndef SEQUENCE_H__
#define SEQUENCE_H__

#include <string>
#include <vector>

using std::string;
using std::vector;

#include "fload.h"

//--------------------------------------------------
class Sequence : protected Fload
{
public:

	Sequence() { size = 0; seq = 0; adr = 0; data = 0; adr_bp = 0; defline = 0; trace_pos = 0; };
	~Sequence() { delete [] seq; delete [] adr; delete [] data; delete [] adr_bp; delete [] defline; ClearMemory(); };

	char* seq;                                                                    // vector<char>
	int* adr;
	int* adr_bp;
	int* data;

	char* defline;                                                                //

	string seqid;

	string trace_seqid;
	int trace_pos;

	int size;                                                                     // vector size

	vector<int> nt;
	double gc;

	bool LoadSequence( const char* name, string cmd_seqid );                                        //
	bool CalculateAddress( const int order  );
	bool LoadData( const char* name );                                            //

	bool SetAddressBP( const int order  );

	void PrintData();
	void PrintSeq();
	void PrintAdr( const int order );

private:

	vector<char> alphabet;                                                        //

	void SetAlphabet();                                                           //
	bool GetSeqFromBuffer();                                                      //
	void ClearMemory() { ClearFileBuffer(); };                                    //
	bool GetAdrFromSeq( const int order, const int adrN, const int adrLow );

	void SetSeqID(string id);

	bool ParseDataLine( string& line );
	inline void SetSiteState( int const pos, const int siteMark );
	inline void SetRegion( const int L, const int R, const int siteMark );
	bool LabelData( string& label, char strand, int L, int R );
};
//--------------------------------------------------
#endif  // SEQUENCE_H__

