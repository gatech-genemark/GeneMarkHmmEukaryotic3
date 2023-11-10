//**************************************************
// file: cmdline.cpp
// project: ehmm3
// Alex Lomsadze
//**************************************************

#include <iostream>
#include <ctime>
#include <string.h>
#include <stdlib.h>

#include "cmdline.h"
#include "common.h"

using std::cout;

//--------------------------------------------------
bool CmdLine::ParseCmdLine( int argc, char** argv  )
{
	bool result = true;
	protein = false;
	nt = false;
	report_best_prob = false;
	verbose = false;
	maskp = 0;
	traceback = false;
	min_gene_length = 0;

	for( int i = 1; (i < argc) && result ; ++i )
	{
		// is not an option - sequence name
		if ( argv[i][0] != '-' )
		{
			if ( sequenceFile.empty() )
				sequenceFile.assign( argv[i] );
			else
				result = false;
		}
		// is option in form of '-' and single letter 
		else if ((argv[i][1] != '\0')&&(argv[i][2] == '\0'))
		{
			switch( argv[i][1] )
			{
			case 'm':
				// value - next argv exists and is not an option
				if (( modelFile.empty() )&&( i + 1 < argc )&&( argv[i+1][0] != '-' ))
				{
					modelFile.assign( argv[ i+1 ] );
					++i;
				}
				else
					result = false;
				break;
			case 'o':
				if (( outputFile.empty() )&&( i + 1 < argc )&&(argv[i+1][0] != '-' ))
				{
					outputFile.assign( argv[ i+1 ] );
					++i;
				}
				else
					result = false;
				break;
			case 'd':
				if (( dataFile.empty() )&&( i + 1 < argc )&&(argv[i+1][0] != '-' ))
				{
					dataFile.assign( argv[ i+1 ] );
					++i;
				}
				else
					result = false;
				break;
			case 'f':
				if (( outputFormat.empty() )&&( i + 1 < argc )&&(argv[i+1][0] != '-' ))
				{
					outputFormat.assign( argv[ i+1 ] );
					++i;
				}
				else
					result = false;
				break;
			case 's':
				if (( seqid.empty() )&&( i + 1 < argc )&&(argv[i+1][0] != '-' ))
				{
					seqid.assign( argv[ i+1 ] );
					++i;
				}
				else
					result = false;
				break;
			case 'b':
				if ( bpintronFile.empty() )
				{
					if (( i + 1 < argc )&&(argv[i+1][0] != '-' ))
					{
						bpintronFile.assign( argv[ i+1 ] );
						++i;
					}
				}
				else
					result = false;
				break;
			case 'k':
				if (i + 1 < argc)
				{
					maskp = atof(argv[i + 1]);
					++i;
				}
				else
					result = false;
				break;
			case 'w':
				if (i + 1 < argc)
				{
					min_gene_length = atoi(argv[i + 1]);
					++i;
				}
				else
					result = false;
				break;
			case 'p':
				// option without value
				if ( !protein )
					protein = true;
				else
					result = false;
				break;
			case 'n':
				// option without value
				if ( !nt )
					nt = true;
				else
					result = false;
				break;
			case 'r':
				// option without value
				if ( !nt )
					report_best_prob = true;
				else
					result = false;
				break;
			case 'v':
				// option without value
				if (!verbose)
					verbose = true;
				else
					result = false;
				break;

			case 'z':
				// option without value
				if (!traceback)
					traceback = true;
				else
					result = false;
				break;

			// not valid option
			default:
				result = false;
			}
		}
		// not valid format
		else
			result = false;
	}

	return result;
};
//--------------------------------------------------
bool CmdLine::Init( int argc, char** argv )
{
	if ( argc == 1 )
	{
		Usage();
		return false;
	}

	if ( !ParseCmdLine( argc, argv ) )
	{
		cout << PROG_NAME << " : error in command line\n";
		return false;
	}

	if ( sequenceFile.empty() )
	{
		cout << PROG_NAME << " : error, sequence file's name is not provided\n";
		return false;
	}

	if ( outputFile.empty() )
		outputFile = sequenceFile + ".lst";

	if ( modelFile.empty() )
	{
		cout << PROG_NAME << " : error, model file's name is not provided\n";
		return false;
	}

	if ( !outputFile.compare( sequenceFile ) || !outputFile.compare( modelFile ) || !outputFile.compare( dataFile ))
	{
		cout << PROG_NAME << " : error, check file names\n";
		return false;
	}

// name must be provided for introninfo to be outputed
//	if ( bpintronFile.empty() )
//		bpintronFile = sequenceFile + ".introninfo";

	if ( outputFormat.empty() )
		outputFormat.assign( LST );

	string formats("|");
	formats.append( LST );
	formats.append( "|" );
	formats.append( GFF3 );
	formats.append( "|" );
	formats.append( GTF );
	formats.append( "|" );
	formats.append( TR );
	formats.append( "|" );

	string out_format = "|" + outputFormat + "|";

	if ( formats.find( out_format ) == string::npos )
	{
		cout << PROG_NAME << " : error, unknown output format found\n";
		return false;
	}

	return true;
};
//--------------------------------------------------
void CmdLine::Usage()
{
	cout << "GeneMark.hmm eukaryotic" << ", version " << PROG_VERSION << "\n\n"
	
	<<	"Usage: " << PROG_NAME << " [options] <sequence file>\n\n"

	<<	"  required parameters:\n"
	<<	"    -m <model file>\n\n"

	<<	"  optional parameters:\n"
	<<	"    -o <output file>\n"
	<<	"    -p write protein translation\n"
	<<	"    -n write nucleotide sequence\n"
	<<	"    -b <output file> output statistics of predicted introns\n"
	<<	"    -d <file name> input for GeneMark.hmm plus\n"
	<<	"    -s <string> sequence tag in GFF output format\n"
	<<	"    -f <format> output prediction in [lst|gff3|gtf] format; default [lst]\n"
	<<	"    -k <number> value for soft-mask penalty\n\n"; 

	cout << "developer options:\n"
	<< "     -z trace back seqid and position\n"
	<< "     -w <number> minimum gene length\n"
	<< "     -r report best path probability\n"
	<< "     -v verbose\n";
};
//--------------------------------------------------

