#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

// First initial test blasso function, just computes sum of double vector values
void F77_NAME(testblasso_f)(double *x, int n, double *res);
// Matrix testing
void F77_NAME(blassomat_f)(double *x, int nrow, int ncol, double *res);

// Define the C wrapper function that calls Fortran
extern SEXP c_testblasso_f(SEXP x){
  const int n = LENGTH(x);
  SEXP res;
  PROTECT(res = allocVector(REALSXP, 1));
  F77_CALL(testblasso_f)(REAL(x), n, REAL(res));
  UNPROTECT(1);
  return(res);
}

// Define the C wrapper function for matrix operation testing
extern SEXP c_blassomat_f(SEXP x, SEXP nrow, SEXP ncol){
	//Rprintf("step1 \n");
	const int nr = asInteger(nrow);
	const int nc = asInteger(ncol);
	//Rprintf("nrow: %i , ncol: %i \n", nr, nc);
	SEXP res;
	//SEXP cols = PROTECT(allocVector(INTSXP, 1
	//Rprintf("step2 \n");
	PROTECT(res = allocVector(REALSXP, nc));
	//Rprintf("step3 \n");
	F77_CALL(blassomat_f)(REAL(x), nr, nc, REAL(res));
	//Rprintf("step4 \n");
	UNPROTECT(1);
	//Rprintf("step5 \n");
	return(res);
}

// Tell R of our available Fortran functions; should probably have one for each model family
static const R_CallMethodDef CallEntries[] = {
  {"c_testblasso_f",	(DL_FUNC) &c_testblasso_f,	1}, // count of input arguments
  {"c_blassomat_f",	(DL_FUNC) &c_blassomat_f,		1}, // count of input arguments
  {NULL,				NULL,						0}
};

// R_init_pckgName
void R_init_blasso(DllInfo *dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
