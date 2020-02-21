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
void F77_NAME(blassocoxfit_f)(int *yevent, int *ytime, double *x, int xnrow, int xncol, int *k, int knrow, int kncol, double *costvec, double *beta);
// Wrapper for Kaisa's DBDC COX optimization subroutine
void F77_NAME(coxdbdc_loop_kits_f)(int nft, int nrecords, int nkits, double *in_vX, int *in_vY, double *in_vC, int *in_vK, double *f_for_k, double *beta_for_k, int iprint, int start);
// Wrapper testing
// extern SEXP c_coxdbdc_test_f(SEXP nf, SEXP nr, SEXP nk, SEXP in_vX, SEXP in_vY, SEXP in_vK, SEXP in_vC){
void F77_NAME(coxdbdc_test_f)(int nft, int nrecords, int nkits, double *in_vX, int *in_vY, double *in_vC, int *in_vK, double *f_for_k, double *beta_for_k, int start, int iprint);

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
	//Rprintf("step1 \n");
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
	//Rprintf("step3 \n");
	F77_CALL(blassocoxfit_f)(INTEGER(yevent), INTEGER(ytime), REAL(x), xnr, xnc, INTEGER(k), knr, knc, REAL(costvec), REAL(beta));
	//Rprintf("step4 \n");
	UNPROTECT(1);
	//Rprintf("step5 \n");
	return(beta);
}

// C wrapper function
extern SEXP c_coxdbdc_loop_kits_f(SEXP nf, SEXP nr, SEXP nk, SEXP in_vX, SEXP in_vY, SEXP in_vK, SEXP in_vC, SEXP startPar, SEXP printPar){

	Rprintf("Defining data constants... \n");
	// Data description
	const int nft = asInteger(nf); // Features, i.e. columns in data
	const int nrecords = asInteger(nr); // Records, i.e. rows in data
	const int nkits = asInteger(nk); // Kits, i.e. how many feature-structuring kits are included
	// Run parameters
	const int iprint = asInteger(printPar);
	const int start = asInteger(startPar); // Starting conditions

	// Integer constants
	Rprintf("nft: %i , nrecords: %i , nkits: %i\n", nft, nrecords, nkits);
	Rprintf("start: %i , iprint: %i \n", start, iprint);
	Rprintf("\n\n");

	// Check first element of all matrices
	Rprintf("in_vX[1:3]: %f %f %f \n", REAL(in_vX)[0], REAL(in_vX)[1], REAL(in_vX)[2]);
	Rprintf("in_vY[1:3]: %i %i %i \n", INTEGER(in_vY)[0], INTEGER(in_vY)[1], INTEGER(in_vY)[2]);
	Rprintf("in_vK[1:3]: %i %i %i \n", INTEGER(in_vK)[0], INTEGER(in_vK)[1], INTEGER(in_vK)[2]);
	Rprintf("in_vC[1:3]: %f %f %f \n", REAL(in_vC)[0], REAL(in_vC)[1], REAL(in_vC)[2]);
	Rprintf("\n\n");

	Rprintf("Preparing output variables... \n");
	SEXP f_for_k; // Objective function values at each k
	SEXP beta_for_k; // Beta values for each k
	//SEXP CPUtime; // Computational time
	Rprintf("Reserving memory for output variables... \n");
	PROTECT(f_for_k = allocVector(REALSXP, nkits));
	PROTECT(beta_for_k = allocVector(REALSXP, nft*nkits));
	//PROTECT(CPUtime = allocVector(REALSXP, 1));

    //       SUBROUTINE coxdbdc_loop_kits(nft, nrecord, nkits, in_vX, in_vY, in_vC, in_vK, &
    //       	& f_for_k, beta_for_k, CPUTime, iprint, start) &
    //       	& BIND(C, name = "coxdbdc_loop_kits_f_")

	Rprintf("Calling Fortran DBDC... \n");
	//F77_CALL(coxdbdc_loop_kits_f)(nft, nrecords, nkits, REAL(in_vX), INTEGER(in_vY), REAL(in_vC), INTEGER(in_vK), REAL(f_for_k), REAL(beta_for_k), REAL(CPUtime), iprint, start);
	F77_CALL(coxdbdc_loop_kits_f)(nft, nrecords, nkits, REAL(in_vX), INTEGER(in_vY), REAL(in_vC), INTEGER(in_vK), REAL(f_for_k), REAL(beta_for_k), iprint, start);

	Rprintf("Wrap up and return output variables... \n");
	UNPROTECT(1);
	return(beta_for_k);
}


// C wrapper function
extern SEXP c_coxdbdc_test_f(SEXP nf, SEXP nr, SEXP nk, SEXP in_vX, SEXP in_vY, SEXP in_vK, SEXP in_vC, SEXP startPar, SEXP printPar){

	Rprintf("Defining data constants... \n");
	// Data description
	const int nft = asInteger(nf); // Features, i.e. columns in data
	//Rprintf("1... \n");
	const int nrecords = asInteger(nr); // Records, i.e. rows in data
	//Rprintf("2... \n");
	const int nkits = asInteger(nk); // Kits, i.e. how many feature-structuring kits are included
	//Rprintf("3... \n");
	// Run parameters
	const int iprint = asInteger(printPar);
	const int start = asInteger(startPar); // Starting conditions

	// Integer constants
	Rprintf("nft: %i , nrecords: %i , nkits: %i\n", nft, nrecords, nkits);
	Rprintf("start: %i , iprint: %i \n", start, iprint);
	Rprintf("\n\n");

	// Check first element of all matrices
	Rprintf("in_vX[1:3]: %f %f %f \n", REAL(in_vX)[0], REAL(in_vX)[1], REAL(in_vX)[2]);
	Rprintf("in_vY[1:3]: %i %i %i \n", INTEGER(in_vY)[0], INTEGER(in_vY)[1], INTEGER(in_vY)[2]);
	Rprintf("in_vK[1:3]: %i %i %i \n", INTEGER(in_vK)[0], INTEGER(in_vK)[1], INTEGER(in_vK)[2]);
	Rprintf("in_vC[1:3]: %f %f %f \n", REAL(in_vC)[0], REAL(in_vC)[1], REAL(in_vC)[2]);
	Rprintf("\n\n");

	Rprintf("Preparing output variables... \n");
	SEXP f_for_k; // Objective function values at each k
	SEXP beta_for_k; // Beta values for each k
	Rprintf("Reserving memory for output variables... \n");
	PROTECT(f_for_k = allocVector(REALSXP, nkits));
	PROTECT(beta_for_k = allocVector(REALSXP, nft*nkits));

    //       SUBROUTINE coxdbdc_loop_kits(nft, nrecord, nkits, in_vX, in_vY, in_vC, in_vK, &
    //       	& f_for_k, beta_for_k, CPUTime, iprint, start) &
    //       	& BIND(C, name = "coxdbdc_loop_kits_f_")

	Rprintf("Calling Fortran DBDC... \n");
	F77_CALL(coxdbdc_test_f)(nft, nrecords, nkits, REAL(in_vX), INTEGER(in_vY), REAL(in_vC), INTEGER(in_vK), REAL(f_for_k), REAL(beta_for_k), start, iprint);

	Rprintf("Wrap up and return output variables... \n");
	UNPROTECT(1);
	return(beta_for_k);
}

// Tell R of our available Fortran functions; should probably have one for each model family
static const R_CallMethodDef CallEntries[] = {
  {"c_testblasso_f",	(DL_FUNC) &c_testblasso_f,	1},
  {"c_blassomat_f",	(DL_FUNC) &c_blassomat_f,		1},
  {"c_blassocoxfit_f", (DL_FUNC) &c_blassocoxfit_f,	1},
  {"c_coxdbdc_loop_kits_f",	(DL_FUNC) &c_coxdbdc_loop_kits_f,	1},
  {"c_coxdbdc_test_f",	(DL_FUNC) &c_coxdbdc_test_f,	1},
  {NULL,				NULL,						0}
};


// R_init_pckgName
void R_init_blasso(DllInfo *dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
