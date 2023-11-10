//**************************************************
// file: cmdline.h
// project: ehmm3
// Alex Lomsadze
//**************************************************

#ifndef CMDLINE_H__
#define CMDLINE_H__

#include <string>

using std::string;

//--------------------------------------------------
class CmdLine
{
public:
	bool Init( int argc, char** argv  );                      //

	string sequenceFile;           //
	string modelFile;              //
	string outputFile;             //
	string dataFile;               //
	string bpintronFile;           //
	bool protein;                  //
	bool nt;                       //
	string outputFormat;           //
	string seqid;                  //

	bool report_best_prob;

	double maskp;

	bool traceback;
	int min_gene_length;

	bool verbose;

private:

	bool ParseCmdLine( int argc, char** argv );               //
	void Usage();                                             //
};
//--------------------------------------------------
#endif // CMDLINE_H__

