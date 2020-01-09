#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

// First initial test blasso function, just computes sum of double vector values
void F77_NAME(testblasso_f)(double *x, int n, double *res);
// Matrix testing
void F77_NAME(blassomat_f)(double *x, int nrow, int ncol, double *res);
// Actual fitting function for budget-aware LASSO (Cox regression)
void F77_NAME(blassocoxfit_f(int *yevent, int *ytime, double *x, int xnrow, int xncol, int *k, int knrow, int kncol, double *costvec, double *beta);

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

// Define the C wrapper function for budget-aware LASSO
extern SEXP c_blassocoxfit_f(SEXP yevent, SEXP ytime, SEXP x, SEXP xnrow, SEXP xncol, SEXP k, SEXP knrow, SEXP kncol, SEXP costvec){
	/Rprintf("step1 \n");
	// Integers for x-matrix dimensions
	const int xnr = asInteger(xnrow);
	const int xnc = asInteger(xncol);
	// Integers for k-matrix dimensions
	const int knr = asInteger(knrow);
	const int knc = asInteger(kncol);

	Rprintf("xnrow: %i , kncol: %i \n", xnr, xnc);
	Rprintf("knrow: %i , kncol: %i \n", knr, knc);
	SEXP beta;
	//SEXP cols = PROTECT(allocVector(INTSXP, 1
	//Rprintf("step2 \n");
	// The new vector of coefficients beta is as long as the number of columns in data matrix x
	PROTECT(beta = allocVector(REALSXP, xnc));
	/Rprintf("step3 \n");
	F77_CALL(blassofit_f)(INTEGER(yevent), INTEGER(ytime), REAL(x), xnr, xnc, INTEGER(k), knr, knx, REAL(costvec), REAL(beta));
	/Rprintf("step4 \n");
	UNPROTECT(1);
	/Rprintf("step5 \n");
	return(beta);
}

// Tell R of our available Fortran functions; should probably have one for each model family
static const R_CallMethodDef CallEntries[] = {
  {"c_testblasso_f",	(DL_FUNC) &c_testblasso_f,	1},
  {"c_blassomat_f",	(DL_FUNC) &c_blassomat_f,		1},
  {"c_blassocoxfit_f", (DL_FUNC) &c_blassocoxfit_f,	1},
  {NULL,				NULL,						0}
};

// R_init_pckgName
void R_init_blasso(DllInfo *dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
