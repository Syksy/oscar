#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

// Blasso Cox DBDC
void F77_NAME(cassocox_f)(double *x, double *y, int *kits, double* cvec, int nrow, int ncol, int nkits, double *beta, double *fperk, int print, int start);

// Define the C wrapper function for matrix operation testing
extern SEXP c_cassocox_f(SEXP x, SEXP y, SEXP kits, SEXP cvec, SEXP nrow, SEXP ncol, SEXP nkits, SEXP print, SEXP start){
	// Define constants (dimensions in data / features)
	const int nr = asInteger(nrow);
	const int nc = asInteger(ncol);
	const int nk = asInteger(nkits);
	const int inprint = asInteger(print);
	const int instart = asInteger(start);
	SEXP beta;
	SEXP fperk;
	// Format output and protect them from garbage collection
	PROTECT(beta = allocVector(REALSXP, nc*nk));
	PROTECT(fperk = allocVector(REALSXP, nk));

	// Call Fortran subroutine
	F77_CALL(cassocox_f)(REAL(x), REAL(y), INTEGER(kits), REAL(cvec), nr, nc, nk, REAL(beta), REAL(fperk),inprint, instart);

	// Create result structure
	SEXP res = PROTECT(allocVector(VECSXP,2));
	SET_VECTOR_ELT(res,0,beta);
	SET_VECTOR_ELT(res,1,fperk);

	// Wrap up and return
	UNPROTECT(1);
	UNPROTECT(1);
	UNPROTECT(1);

	//return(beta);
	return(res);
}
// Tell R of our available Fortran functions; should probably have one for each model family
static const R_CallMethodDef CallEntries[] = {
  {"c_cassocox_f",	(DL_FUNC) &c_cassocox_f,		9},
  {NULL,				NULL,						0}
};


// R_init_pckgName
void R_init_casso(DllInfo *dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
