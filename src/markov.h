//**************************************************
// file: markov.h
// project: ehmm3
//**************************************************

#ifndef MARKOV_H__
#define MARKOV_H__

#include <stdio.h>

#include "matrix.h"

//--------------------------------------------------
class MarkovCod : private Matrix
{
public:
	MarkovCod() { d1=d2=d3=d1a=d2a=d3a=r1=r2=r3=r1a=r2a=r3a = NULL; ddd_1=ddd_2=ddd_3=0;  };
	~MarkovCod(){ delete [] d1; delete [] d2; delete [] d3;
	        delete [] d1a; delete [] d2a; delete [] d3a;
			delete [] r1; delete [] r2; delete [] r3;
			delete [] r1a; delete [] r2a; delete [] r3a;
			delete [] ddd_1; delete [] ddd_2; delete [] ddd_3; };

	double* d1a; double* d1;
	double* d2a; double* d2;
	double* d3a; double* d3;

	double* r1a; double* r1;
	double* r2a; double* r2;
	double* r3a; double* r3;

	bool Set( int order, double probN );

	double** ddd_1;
	double** ddd_2;
	double** ddd_3;

	void ResetDir( int ph )
	{
		dph = dirArr[ph];
	}

	void ResetRev( int ph )
	{
		rph = revArr[ph];
	}

	double GetDirA( int adr )
	{
		double x = dph->matA[ adr ];
		dph = dph->next;
		return x;
	}

	double GetDir( int adr )
	{
		double x = dph->mat[ adr ];
		dph = dph->next;
		return x;
	}

	double GetRevA( int adr )
	{
		double x = rph->matA[ adr ];
		rph = rph->next;
		return x;
	}

	double GetRev( int adr )
	{
		double x = rph->mat[ adr ];
		rph = rph->next;
		return x;
	}

private:

	struct Phase {
		double* mat;
		double* matA;
		Phase* next;
	} dph0, dph1, dph2, *dph, rph0, rph1, rph2, *rph;

	Phase* dirArr[3];
	Phase* revArr[3];
};
//--------------------------------------------------
class MarkovNon : private Matrix
{
public:

	MarkovNon() { n = na = NULL; ddd = 0; };
	~MarkovNon() { delete [] n; delete [] na; delete [] ddd; };

	double* na;	double* n;

	bool Set( int order, double probN );

	double** ddd;
};
//--------------------------------------------------
class MarkovIntr : private Matrix
{
public:

	MarkovIntr() { d=da=r=ra = 0; ddd=rrr = 0;  };
	~MarkovIntr() { delete [] d; delete [] da; delete [] r; delete [] ra; delete [] ddd; delete [] rrr; };

	double* da;	double* d;

	double* ra;	double* r;

	bool Set( int order, double probN );

	double** ddd;
	double** rrr;
};
//--------------------------------------------------
#endif  // MARKOV_H__

