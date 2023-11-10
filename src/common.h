//**************************************************
// file: common.h
// project: ehmm3
// Alex Lomsadze
//**************************************************

#ifndef COMMON_H__
#define COMMON_H__

#include <math.h>

//--------------------------------------------------
//#define __DEBUG__

const char PROG_NAME[] = "gmhmme3";
const char PROG_VERSION[] = "3.68.pub";

// order of letters is important !DON'T CHANGE!
enum NUCLEOTIDE { A, C, G, T, N };
const char LETTER[5] = { 'a', 'c', 'g', 't', 'n' };
const char LETTER_UP[5] = { 'A', 'C', 'G', 'T', 'N' };
const char LETTER_UP_COMP[5] = { 'T', 'G', 'C', 'A', 'N' };

const int dirTAAcodon = (T<<4|A<<2|A);
const int dirTGAcodon = (T<<4|G<<2|A);
const int dirTAGcodon = (T<<4|A<<2|G);
const int revTAAcodon = (T<<4|T<<2|A);
const int revTGAcodon = (T<<4|C<<2|A);
const int revTAGcodon = (C<<4|T<<2|A);

const int dirStart = 0x1;       //    0000 0000 0000 0001      1
const int dirEnd   = 0x2;       //    0000 0000 0000 0010      2
const int dirDn    = 0x4;       //    0000 0000 0000 0100      4
const int dirAc    = 0x8;       //    0000 0000 0000 1000      8

const int revStart = 0x10;      //    0000 0000 0001 0000     16
const int revEnd   = 0x20;      //    0000 0000 0010 0000     32
const int revDn    = 0x40;      //    0000 0000 0100 0000     64
const int revAc    = 0x80;      //    0000 0000 1000 0000    128

const int mapIni   = 0x100;     //    0000 0001 0000 0000    256
const int mapTer   = 0x200;     //    0000 0010 0000 0000    512

const int cleanMark  = 0x1000;  //    0001 0000 0000 0000   4096
const int codMark    = 0x2000;  //    0010 0000 0000 0000   8192
const int nonMark    = 0x4000;  //    0100 0000 0000 0000  16384
const int intronMark = 0x8000;  //    1000 0000 0000 0000  32768
const int interMark  = 0x10000; // 01 0000 0000 0000 0000  65536
const int maskMark   = 0x20000; // 10 0000 0000 0000 0000

const int internalSite = 0xFF;  //    0000 0000 1111 1111    255
const int seseSite     = 0x33;  //    0000 0000 0011 0011     51
const int anyEnd       = 0x22;  //    0000 0000 0010 0010     34

typedef struct MapStruct
{
	int key;
	int pos;
	char fr;
	char ph;
	int data;
} MapType;

const double log0 = -1.0e120;
const double PLUS = 100;
//const double MINUS = 0;
const double MINUS_IN_INTRON = -100;

const double xp = log(0.25);

// output formats
const char LST []  = "lst";
const char GFF3[]  = "gff3";
const char GTF []  = "gtf";
const char TR  []  = "tr";

// output labels
// set equal size for formating
// 10 + 1 '\0'
const int _LABEL_SIZE = 11;
const char _INTERGENIC[] = "Intergenic";
const char _INTRON[]     = "Intron    ";
const char _INITIAL[]    = "Initial   ";
const char _INTERNAL[]   = "Internal  ";
const char _TERMINAL[]   = "Terminal  ";
const char _SINGLE[]     = "Single    ";

// external data labels
const char _x_SINGLE[]         = "Single";           // +
const char _x_INITIAL[]        = "Initial";          // +
const char _x_INTERNAL[]       = "Internal";         // +
const char _x_TERMINAL[]       = "Terminal";         // +
const char _x_INTRON[]         = "intron";           // +
const char _x_INTRON_uc[]      = "Intron";           // +
const char _x_INTERGENIC[]     = "Intergenic";       // +
const char _x_DONOR[]          = "Donor";            // +
const char _x_ACCEPTOR[]       = "Acceptor";         // + 
const char _x_START_CODON[]    = "start_codon";      // +
const char _x_START_CODON_uc[] = "Start_codon";      // +
const char _x_END_CODON[]      = "stop_codon";       // + 
const char _x_END_CODON_uc[]   = "Stop_codon";       // +

const char _x_P_EXON_DON[]     = "Part_don";         // +
const char _x_P_EXON_ACC[]     = "Part_acc";         // +
const char _x_P_EXON_TER[]     = "Part_ter";         // +
const char _x_P_EXON_INI[]     = "Part_ini";         // +
const char _x_P_EXON[]         = "Part_exon";        // +

const char _x_UTR_INTRON[]     = "UTR_Intron";       // -
const char _x_UTR[]            = "UTR_Exon";         // -
const char _x_NON_CODING[]     = "Non_Coding";       // +

const char _x_P_INTRON[]       = "Part_intron";      // +
const char _x_P_INTERGENIC[]   = "Part_intergenic";  // +

const char _x_MASK[]           = "mask";
//--------------------------------------------------
#endif // COMMON_H__

