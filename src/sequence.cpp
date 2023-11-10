//**************************************************
// file: sequence.cpp
// project: ehmm3
// Alex Lomsadze
//**************************************************

#include <cctype>
#include <iostream>
#include <fstream>
#include <bitset> 

#include <string.h> 
#include <stdlib.h> 

using std::cout;
using std::ifstream;
using std::ofstream;
using std::basic_string;
using std::bitset;

#include "sequence.h"
#include "common.h"

//--------------------------------------------------
void Sequence::SetAlphabet()
{
	alphabet.assign( 256, N );

	alphabet['a'] = A;
	alphabet['A'] = A;

	alphabet['c'] = C;
	alphabet['C'] = C;

	alphabet['g'] = G;
	alphabet['G'] = G;

	alphabet['t'] = T;
	alphabet['T'] = T;
	alphabet['u'] = T;
	alphabet['U'] = T;
};
//--------------------------------------------------
bool Sequence::GetSeqFromBuffer()
{
	if ( seq ) delete [] seq;
	seq = new char[ fileSize ];
	if ( !seq ) return false;

	char* letter =  fileBuffer;
	const char* letterLast = fileBuffer + fileSize;

	// skip white space
	while ( (letter < letterLast) && isspace(*letter) )
		++letter;

	// if FASTA format - skip line
	if ((letter < letterLast) && ( *letter == '>' ))
	{
		char* defline_start = letter;

		while ((letter < letterLast)&&(*letter != '\n')&&(*letter != '\r') )
			++letter;

		char* defline_end = letter;

		if ( defline ) delete [] defline;
		defline = new char[ defline_end - defline_start + 1 ];
		if ( !defline ) return false;

		strncpy( defline, defline_start, defline_end - defline_start );
		defline[ defline_end - defline_start ] = '\0';
	}

	nt.assign( 5, 0 );
	size = 0;
	
	// read sequence
	while ( letter < letterLast)
	{
		if ( isalpha( *letter ) )
		{
			seq[ size ] = alphabet[ *letter ];
			++nt[ seq[size] ];
			++size;
		}
		++letter;
	}

	gc = nt[A] + nt[C] + nt[G] + nt[T];

	if ( gc )
		gc = 100.0*(nt[C] + nt[G])/gc;

	return true;
};
//--------------------------------------------------
bool Sequence::LoadSequence(const char* name, string cmd_seqid)
{
	if (!LoadFile(name)) return false;

	SetAlphabet();

	if (!GetSeqFromBuffer())
	{
		cout << PROG_NAME << " : error, not enough memory to load sequence\n";
		return false;
	}

	ClearMemory();

	if (!size)
	{
		cout << PROG_NAME << " : error, sequence has size " << size << "\n";
		return false;
	}

	SetSeqID(cmd_seqid);

	return true;
};
//--------------------------------------------------
void Sequence::SetSeqID(string id)
{
	// the order of sequence ID processing:
	//     command line parameter
	//     from defline of sequence - first defline word
	//     string "seq"

	if (!id.empty())
	{
		seqid = id;
		return;
	}

	// parse defline  ---  loop is used just to simplify the code --- break from block

	while ( defline && strlen(defline) )
	{
		string tmpstr;
		std::size_t L, R;
		tmpstr.assign( defline );

		L = tmpstr.find_first_not_of( " \t", 1 );
		if ( L == std::string::npos ) break;
		R = tmpstr.find_first_of( " \t", L );

		// seqid is assigned here
		seqid = tmpstr.substr( L, R - L );
		if ( R == std::string::npos ) break;

		L = tmpstr.find_first_not_of( " \t", R );
		if ( L == std::string::npos ) break;
		R = tmpstr.find_first_of( " \t", L );
		trace_seqid = tmpstr.substr( L, R - L ); 
		if ( R == std::string::npos ) break;

		L = tmpstr.find_first_not_of( " \t", R );
		if ( L == std::string::npos ) break;
		R = tmpstr.find_first_of( " \t", L );
		tmpstr = tmpstr.substr( L, R - L );
		trace_pos = atoi( tmpstr.c_str() );

		if (trace_pos == 0)
			trace_seqid.clear();

		break;
	}

	if (seqid.empty())
		seqid.assign("seq");

	return;
};
//--------------------------------------------------
bool Sequence::CalculateAddress( const int order )
{
	if ( (order + 1) > size )
	{
		cout << PROG_NAME << " : error, size of sequence " << size << ", order " << order << "\n";
		return false;
	}
	if ( 8*sizeof(int) <= 2*(order + 1) )
	{
		cout << PROG_NAME << " : error, can't convert address for order " << order << " to 'int'\n";
		return false;
	}
	if ( !GetAdrFromSeq( order, 4<<(2*order), 0 ) )
	{
		cout << PROG_NAME << ": error, not enough memory to load address\n";
		return false;
	}

	return true;
};
//--------------------------------------------------
bool Sequence::SetAddressBP( const int order  )
{
	int* tmp;
	if (!adr)
	{
		cout << PROG_NAME << ": error in lg" << std::endl;
		exit(1);
	}
	tmp = adr;
	adr = 0;

	if ( !CalculateAddress( order ) )
		return false;

	adr_bp = adr;
	adr = tmp;

	return true;
};
//--------------------------------------------------
//  e1 e2 e3 e4 e5 e6 e7 e8 .... en  sequence
//       |e3 e4 e5 e6|               take "order" + 1 elements 
//                                   order == 3, 4 elements
//  m3 = e3*4^3 + e4*4^2 + e5*4^1 + e6*4^0
//  record m3 in position e3
//--------------------------------------------------
bool Sequence::GetAdrFromSeq( const int order, const int adrN, const int adrLow )
{
	if ( adr ) delete [] adr;
	adr	= new int[ size ];
	if ( !adr ) return false;

	int address = 0;
	int trace = 0;
	char* letter = seq + size - 1;
	int* current = adr + size - 1;

	for ( int i = 0; i < order; ++i )
	{
		if ( *letter != N ) 
			address += ( (*letter) << (2*i) );
		else
			trace += ( 3 << (2*i) );

		if ( !trace )
			*current = adrLow;
		else
			*current = adrN;

		--letter;
		--current;
	}

	while ( letter >= seq )
	{
		if ( *letter != N )
			address += ( (*letter) << (2*order) );
		else
			trace += ( 3 << (2*order) );

		if ( !trace )
			*current = address;
		else
			*current = adrN;

		address >>= 2;
		trace >>= 2;
		--current;
		--letter;
	}

	return true;
};
//--------------------------------------------------
bool Sequence::LabelData( string& label, char strand, int L, int R )
{
	if     (( strand == '+' )&& ( !label.compare( _x_SINGLE)))
	{
		// *               *
		// A T G . . . T A A
		// |               |
		SetRegion( L + 1, R - 1, cleanMark + codMark );
		SetSiteState( L, dirStart );
		SetSiteState( R, dirEnd );
	}
	else if (( strand == '-' )&&( !label.compare( _x_SINGLE )))
	{
		// *               *
		// T T A . . . C A T
		// |               |
		SetRegion( L + 1, R - 1, cleanMark + codMark );
		SetSiteState( L, revEnd );
		SetSiteState( R, revStart );
	}
	else if (( strand == '+' )&&( !label.compare( _x_INITIAL )))
	{
		// *         *
		// A T G . . . G T
		// |           |
		SetRegion( L + 1, R, cleanMark + codMark );
		SetSiteState( L, dirStart );
		SetSiteState( R + 1, dirDn );
	}
	else if (( strand == '-' )&&( !label.compare( _x_INITIAL )))
	{
		//      *         *
		//  A C . . . C A T
		//    |           |
		SetRegion( L, R - 1, cleanMark + codMark );
		SetSiteState( L - 1, revDn );
		SetSiteState( R, revStart );
	}
	else if (( strand == '+' )&&( !label.compare( _x_INTERNAL )))
	{
		//      *   *
		//  A G . . . G T
		//    |       |
		SetRegion( L, R, cleanMark + codMark );
		SetSiteState( L - 1, dirAc );
		SetSiteState( R + 1, dirDn );
	}
	else if (( strand == '-' )&&( !label.compare( _x_INTERNAL )))
	{
		//      *   *
		//  A C . . . C T
		//    |       |
		SetRegion( L, R, cleanMark + codMark );
		SetSiteState( L - 1, revDn );
		SetSiteState( R + 1, revAc );
	}
	else if (( strand == '+' )&&( !label.compare( _x_TERMINAL )))
	{
		//      *         *
		//  A G . . . T A A
		//    |           |
		SetRegion( L, R - 1, cleanMark + codMark );
		SetSiteState( L - 1, dirAc );
		SetSiteState( R, dirEnd );
	}
	else if (( strand == '-' )&&( !label.compare( _x_TERMINAL )))
	{
		// *         *
		// T T A . . . C T
		// |           |
		SetRegion( L + 1, R, cleanMark + codMark );
		SetSiteState( L, revEnd );
		SetSiteState( R + 1, revAc );
	}
	else if (( strand == '+' )&&( !label.compare( _x_INTRON )))
	{
		//  *           *
		//  G T . . . A G
		//  |           |
		SetRegion( L + 1, R - 1, cleanMark + intronMark + nonMark );
		SetSiteState( L, dirDn );
		SetSiteState( R, dirAc );
	}
	else if (( strand == '-' )&&( !label.compare( _x_INTRON )))
	{
		//  *           *
		//  C T . . . A C
		//  |           |
		SetRegion( L + 1, R - 1, cleanMark + intronMark + nonMark );
		SetSiteState( L, revAc );
		SetSiteState( R, revDn );
	}
	else if (( strand == '.' )&&( !label.compare( _x_INTERGENIC )))
	{
		//       *   *
		// T A A . . . A T G
		//     |       |
		SetRegion( L, R, cleanMark + interMark + nonMark );
		SetSiteState( L - 1, seseSite );
		SetSiteState( R + 1, seseSite );
	}
	else if (( strand == '+' )&&( !label.compare( _x_P_EXON_DON )))
	{
		//  *   *
		//  . . . G T
		//        |
		SetRegion( L + 1, R, cleanMark + codMark );
		SetSiteState( R + 1, dirDn );
	}
	else if (( strand == '-' )&&( !label.compare( _x_P_EXON_DON )))
	{
		//      *   *
		//  A C . . .
		//    |      
		SetRegion( L, R - 1, cleanMark + codMark );
		SetSiteState( L - 1, revDn );
	}
	else if (( strand == '+' )&&( !label.compare( _x_P_EXON_ACC )))
	{
		//      *   *
		//  A G . . .
		//    |      
		SetRegion( L, R - 1, cleanMark + codMark );
		SetSiteState( L - 1, dirAc );
	}
	else if (( strand == '-' )&&( !label.compare( _x_P_EXON_ACC )))
	{
		//  *   *
		//  . . . C T
		//        |
		SetRegion( L + 1, R, cleanMark + codMark );
		SetSiteState( R + 1, revAc );
	}
	else if (( strand == '+' )&&( !label.compare( _x_P_EXON_INI )))
	{
		// *         *
		// A T G . . .
		// |           
		SetRegion( L + 1, R - 1, cleanMark + codMark );
		SetSiteState( L, dirStart );
	}
	else if (( strand == '-' )&&( !label.compare( _x_P_EXON_INI )))
	{
		//   *         *
		//   . . . C A T
		//             |
		SetRegion( L + 1, R - 1, cleanMark + codMark );
		SetSiteState( R, revStart );
	}
	else if (( strand == '+' )&&( !label.compare( _x_P_EXON_TER )))
	{
		//   *         *
		//   . . . T A A
		//             |
		SetRegion( L + 1, R - 1, cleanMark + codMark );
		SetSiteState( R, dirEnd );
	}
	else if (( strand == '-' )&&( !label.compare( _x_P_EXON_TER )))
	{
		// *         *
		// T T A . . .
		// |          
		SetRegion( L + 1, R - 1, cleanMark + codMark );
		SetSiteState( L, revEnd );
	}
	else if (( strand == '+' )&&( !label.compare( _x_P_EXON )))
	{
		//      *   *
		//      . . . 
		//          
		SetRegion( L + 1, R - 1, cleanMark + codMark );
	}
	else if (( strand == '-' )&&( !label.compare( _x_P_EXON )))
	{
		//      *   *
		//      . . .
		//           
		SetRegion( L + 1, R - 1, cleanMark + codMark );
	}
	else if (!label.compare( _x_UTR_INTRON ))
	{
		//  *           *
		//  G T . . . A G
		//  |           |
		SetRegion( L, R, cleanMark + nonMark );
	}
	else if ( !label.compare( _x_UTR ))
	{
		//   *   *
		//   . . . 
		//   |   |
		SetRegion( L, R, cleanMark + nonMark );	
	}
	else if ( !label.compare( _x_NON_CODING ))
	{
		//   *   *
		//   . . . 
		//   |   |
		SetRegion( L, R, cleanMark + nonMark );	
	}
	else if ( !label.compare( _x_P_INTRON ))
	{
		//   *   *
		//   . . . 
		//   |   |
		SetRegion( L, R, cleanMark + intronMark + nonMark );	
	}
	else if ( !label.compare( _x_P_INTERGENIC ))
	{
		//   *   *
		//   . . . 
		//   |   |
		SetRegion( L, R, cleanMark + interMark + nonMark );	
	}
	else if (( strand == '+' )&&( !label.compare( _x_DONOR )))
	{
		//  *
		//  G T
		//  |
		SetSiteState( L, dirDn );
	}
	else if (( strand == '+' )&&( !label.compare( _x_ACCEPTOR )))
	{
		//   *
		//  AG
		//   |
		SetSiteState( R, dirAc );
	}
	else if (( strand == '+' )&&( !label.compare( _x_START_CODON )))
	{
		//  *
		//  A T G
		//  |
		SetSiteState( L, dirStart );
	}
	else if (( strand == '+' )&&( !label.compare( _x_END_CODON )))
	{
		//     *
		// T A A
		//     |
		SetSiteState( R, dirEnd );
	}
	else if (( strand == '-' )&&( !label.compare( _x_DONOR )))
	{
		//    *
		//  A C
		//    |
		SetSiteState( R, revDn );
	}
	else if (( strand == '-' )&&( !label.compare( _x_ACCEPTOR )))
	{
		//  *
		//  CT
		//  |
		SetSiteState( L, revAc );
	}
	else if (( strand == '-' )&&( !label.compare( _x_START_CODON )))
	{
		//	    *
		//  C A T
		//      |
		SetSiteState( R, revStart );
	}
	else if (( strand == '-' )&&( !label.compare( _x_END_CODON )))
	{
		// *
		// T T A
		// |
		SetSiteState( L, revEnd );
	}
	else if (!label.compare(_x_MASK))
	{
		//   *   *
		//   . . . 
		//   |   |
		SetRegion(L, R, maskMark);
	}
	else
		return false;

	return true;
};
//--------------------------------------------------
bool Sequence::ParseDataLine( string& line )
{
	basic_string <char>::size_type L, R;
	
	L = line.find_first_not_of( " \t" );
	
	// skip white space and comments
	if (( L == -1 )||( line[L] == '#' ))
		return true;

	vector< string > gff;
	gff.reserve(9);

	R = line.find_first_of( '\t', L );

	for( int i = 0 ; (i < 9)&&(L != -1); ++i )
	{
		gff.push_back( line.substr( L, R - L ) );

		L = line.find_first_not_of( '\t', R );
		R = line.find_first_of( '\t', L );
	}

	if ( gff.size() != 9 )
	{
		cout << "Data format error, 9 columns is expected but " << gff.size() << " was found in: " << line << std::endl;
		return false;
	}

	// check match in sequence ID first
	if ( seqid.compare( gff[0] ) )
		return true;

// gff[1] - source - not relevant
	string key = gff[2];
	int start = atoi( gff[3].c_str() );
	int end  = atoi( gff[4].c_str() );
// gff[5] - score - not relevant
	char strand = gff[6].at(0);
// gff[7] - frame - future use
// gff[8] - transcript id - future use

	if (( start > end )||( start < 1 )||( end > size ))
	{
		cout << "Data format error, check start and/or end coordinates: " << line << std::endl;
		return false;
	}

	if (( strand != '+' )&&( strand != '-' )&&( strand != '.'))
	{
		cout << "Data format error, unexpected symbol found in strand field: " << line << std::endl;
		return false;
	}

// Allow low or uppercase in first letter

	if ( !key.compare(_x_INTRON_uc) )
		key.assign( _x_INTRON );
	else if ( !key.compare( _x_START_CODON_uc) )
		key.assign( _x_START_CODON );
	else if ( !key.compare( _x_END_CODON_uc) )
		key.assign( _x_END_CODON );

	if ( !LabelData( key, strand, start - 1, end - 1 ) )
	{
		cout << "Data format error, unexpected feature " << key << " in : " << line << std::endl;
		return false;
	}

	return true;
};
//--------------------------------------------------
bool Sequence::LoadData( const char* name )
{
	if ( data ) delete [] data;
	data = new int[ size ];
	if ( !data )
	{
		cout << PROG_NAME << " : error, out of memory on -d " << name << "\n"; 
		return false;
	}

	memset( data, 0, size*sizeof(int) );

	ifstream data_file( name );

	if ( !data_file.is_open() )
	{
		cout << PROG_NAME << " : error, opening file " << name << "\n"; 
		return false;
	}

	string line;

	while( !data_file.eof() )
	{
		getline( data_file, line );
		if ( !ParseDataLine( line ) )
			cout << PROG_NAME << " : warning, file " << name << " line ignored : " << line << "\n";
	}

	data_file.close();

	return true;
};
//--------------------------------------------------
inline void Sequence::SetSiteState( const int pos, const int siteMark )
{
	data[pos] |= siteMark;
};
//--------------------------------------------------
inline void Sequence::SetRegion( const int L, const int R, const int regionMark )
{
	for( int i = L; i <= R; ++i )
		data[i] |= regionMark;
};
//--------------------------------------------------
void Sequence::PrintData()
{
	if ( !data )
	{
		cout << "no external information provided\n";
		return;
	}

	string name( "external_data.log" );
	ofstream log_file( name.c_str() );

	if ( !log_file.is_open() )
	{
		cout << PROG_NAME << " : error, opening file " << name << "\n"; 
		exit(1);
	}

	string label;

	for( int i = 0; i < size; ++i )
	{
		label.clear();
		int current = data[i];

		if ( current & cleanMark )
		{
			label.append( "CleanMark " );
			current -= cleanMark;
		}

		if ( current & dirDn )
		{
			label.append( "DirDn " );
			current -= dirDn;
		}

		if ( current & dirAc )
		{
			label.append( "DirAc " );
			current -= dirAc;
		}

		if ( current & dirStart )
		{
			label.append( "DirStart " );
			current -= dirStart;
		}

		if ( current & revEnd )
		{
			label.append( "RevEnd " );
			current -= revEnd;
		}

		if (current & revDn)
		{
			label.append("RevDn ");
			current -= revDn;
		}

		if (current & revAc)
		{
			label.append("RevAc ");
			current -= revAc;
		}

		if (current & revStart)
		{
			label.append("RevStart ");
			current -= revStart;
		}

		if (current & revEnd)
		{
			label.append("RevEnd ");
			current -= revEnd;
		}

		if ( current & codMark )
		{
			label.append( "CodMark " );
			current -= codMark;
		}

		if ( current & nonMark )
		{
			label.append( "NonMark " );
			current -= nonMark;
		}

		if (current & intronMark)
		{
			label.append("IntronMark ");
			current -= intronMark;
		}
		
		if (current & interMark)
		{
			label.append("InterMark ");
			current -= interMark;
		}

		if (current & maskMark)
		{
			label.append("MaskMark ");
			current -= maskMark;
		}

		if (current)
		{
			label.append("? ");
		}

		log_file << i << " " << data[i] << " " << current << " " << label << "\n";
	}

	log_file.close();
};
//--------------------------------------------------
void Sequence::PrintSeq()
{
	if ( !seq )
	{
		cout << "no sequence provided\n";
		return;
	}

	string name( "sequence.log" );
	ofstream log_file( name.c_str() );

	if ( !log_file.is_open() )
	{
		cout << PROG_NAME << " : error, opening file " << name << "\n"; 
		exit(1);
	}

	for( int i = 0; i < size; ++i )
		log_file << LETTER[ seq[i] ] << " " << i << "\n";

	log_file.close();
};
//--------------------------------------------------
void Sequence::PrintAdr( const int order )
{
	if ( !adr )
	{
		cout << "no address provided\n";
		return;
	}

	string name( "address.log" );
	ofstream log_file( name.c_str() );

	if ( !log_file.is_open() )
	{
		cout << PROG_NAME << " : error, opening file " << name << "\n"; 
		exit(1);
	}

	string label( order + 1, 'N' );
	int max_norm_adr = (4<<(2*order)) - 1;

	for( int i = 0; i < size; ++i )
	{
		int x = adr[i];

		if ( x > max_norm_adr )
			log_file << i << " " << bitset<16>( adr[i] ) << " " << adr[i] << " N\n";
		else
		{
			for( int j = 0; j < order + 1; ++j )
				label[ order - j  ] = LETTER[ (( x >> (j*2) ) & 3 ) ];

			log_file << i << " " << bitset<16>( adr[i] ) << " " << adr[i] << " " << label << "\n";
		}
	}

	log_file.close();
};
//--------------------------------------------------

