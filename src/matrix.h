//**************************************************
// file: matrix.h
// project: ehmm3
//**************************************************

#ifndef MATRIX_H__
#define MATRIX_H__

//--------------------------------------------------
class Matrix
{
public:

	bool RevFromDir( double* dir, double** rev , int order );
	bool AbsToFirst( double* abs, double** rel, int order );
	bool HiToLo( double* hi, double** lo, int hiOrder );
	void Log( double* ptr , int size );
	bool AllocArray( double** data, int size );

	void AbsoluteToRelativeRight( double* ref, int order );

private:

	void DoRevFromDir( double* dir, double* rev , int order );
	void DoHiToLo( double* hi, double* lo, int hiOrder );
	void DoAbsToFirst( double* abs, double* rel, int order );
};
//--------------------------------------------------
#endif // MATRIX_H__

