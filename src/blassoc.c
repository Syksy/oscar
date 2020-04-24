#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

// Blasso Cox DBDC
void F77_NAME(blassocox_kaikki_f)(double *x, double *y, int *kits, double* cvec, int nrow, int ncol, int nkits, double *beta, double *fperk);

// Define the C wrapper function for matrix operation testing
extern SEXP c_blassocox_kaikki_f(SEXP x, SEXP y, SEXP kits, SEXP cvec, SEXP nrow, SEXP ncol, SEXP nkits){
	// Define constants (dimensions in data / features)
	const int nr = asInteger(nrow);
	const int nc = asInteger(ncol);
	const int nk = asInteger(nkits);
	SEXP beta;
	SEXP fperk;
	// Format output and protect them from garbage collection
	PROTECT(beta = allocVector(REALSXP, nc*nk));
	PROTECT(fperk = allocVector(REALSXP, nk));

	// Call Fortran subroutine
	F77_CALL(blassocox_kaikki_f)(REAL(x), REAL(y), INTEGER(kits), REAL(cvec), nr, nc, nk, REAL(beta), REAL(fperk));

	// Wrap up and return
	UNPROTECT(1);
	UNPROTECT(1);
	return(beta);
}
// Tell R of our available Fortran functions; should probably have one for each model family
static const R_CallMethodDef CallEntries[] = {
  {"c_blassocox_f",	(DL_FUNC) &c_blassocox_kaikki_f,		9},
  {NULL,				NULL,						0}
};


// R_init_pckgName
void R_init_blasso(DllInfo *dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
