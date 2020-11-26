#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

// family Cox
void F77_NAME(casso_cox_f)(double *x, double *y, int *kits, double* cvec, int nrow, int ncol, int nkits, double *beta, double *fperk, int print, int start, int kmax,
							int inmrounds, int inmit, int inmroundsesc, int inb1, int inb2, int inb, double inm, double inmclarke, double inc, 
							double inrdec, double inrinc, double ineps1, double ineps, double incrittol);
// family Gaussian (MSE)
void F77_NAME(casso_mse_f)(double *x, double *y, int *kits, double* cvec, int nrow, int ncol, int nkits, double *beta, double *fperk, int print, int start, int kmax,
							int inmrounds, int inmit, int inmroundsesc, int inb1, int inb2, int inb, double inm, double inmclarke, double inc, 
							double inrdec, double inrinc, double ineps1, double ineps, double incrittol);
// family Logistic
void F77_NAME(casso_logistic_f)(double *x, int *y, int *kits, double* cvec, int nrow, int ncol, int nkits, double *beta, double *fperk, int print, int start, int kmax,
							int inmrounds, int inmit, int inmroundsesc, int inb1, int inb2, int inb, double inm, double inmclarke, double inc, 
							double inrdec, double inrinc, double ineps1, double ineps, double incrittol);

// Define the C wrapper function for Cox regression
extern SEXP c_casso_cox_f(SEXP x, SEXP y, SEXP kits, SEXP cvec, SEXP nrow, SEXP ncol, SEXP nkits, SEXP print, SEXP start, SEXP kmax,
						SEXP mrounds, SEXP mit, SEXP mroundsesc, SEXP b1, SEXP b2, SEXP b, SEXP m, SEXP mclarke, 
						SEXP c, SEXP rdec, SEXP rinc, SEXP eps1, SEXP eps, SEXP crittol){
	// Define constants (dimensions in data / features)
	const int nr = asInteger(nrow);
	const int nc = asInteger(ncol);
	const int nk = asInteger(nkits);
	const int inkmax = asInteger(kmax);
	const int inprint = asInteger(print);
	const int instart = asInteger(start);
	const int inmrounds = asInteger(mrounds);
	const int inmit = asInteger(mit);
	const int inmroundsesc = asInteger(mroundsesc);
	const int inb1 = asInteger(b1);
	const int inb2 = asInteger(b2);
	const int inb = asInteger(b);
	const double inm = asReal(m);
	const double inmclarke = asReal(mclarke);
	const double inc = asReal(c);
	const double inrdec = asReal(rdec);
	const double inrinc = asReal(rinc);
	const double ineps1 =asReal(eps1);
	const double ineps = asReal(eps);
	const double incrittol = asReal(crittol);
	
	SEXP beta;
	SEXP fperk;
	// Format output and protect them from garbage collection
	PROTECT(beta = allocVector(REALSXP, nc*nk));
	PROTECT(fperk = allocVector(REALSXP, nk));

	
	// Call Fortran subroutine
	F77_CALL(casso_cox_f)(REAL(x), REAL(y), INTEGER(kits), REAL(cvec), nr, nc, nk, REAL(beta), REAL(fperk),inprint, instart, inkmax,
						inmrounds, inmit, inmroundsesc, inb1, inb2, inb, inm, 
						inmclarke,inc, inrdec, inrinc, ineps1, ineps, incrittol);

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

// Define the C wrapper function for Gaussian family
extern SEXP c_casso_mse_f(SEXP x, SEXP y, SEXP kits, SEXP cvec, SEXP nrow, SEXP ncol, SEXP nkits, SEXP print, SEXP start, SEXP kmax,
						SEXP mrounds, SEXP mit, SEXP mroundsesc, SEXP b1, SEXP b2, SEXP b, SEXP m, SEXP mclarke, 
						SEXP c, SEXP rdec, SEXP rinc, SEXP eps1, SEXP eps, SEXP crittol){
	// Define constants (dimensions in data / features)
	const int nr = asInteger(nrow);
	const int nc = asInteger(ncol);
	const int nk = asInteger(nkits);
	const int inkmax = asInteger(kmax);
	const int inprint = asInteger(print);
	const int instart = asInteger(start);
	const int inmrounds = asInteger(mrounds);
	const int inmit = asInteger(mit);
	const int inmroundsesc = asInteger(mroundsesc);
	const int inb1 = asInteger(b1);
	const int inb2 = asInteger(b2);
	const int inb = asInteger(b);
	const double inm = asReal(m);
	const double inmclarke = asReal(mclarke);
	const double inc = asReal(c);
	const double inrdec = asReal(rdec);
	const double inrinc = asReal(rinc);
	const double ineps1 =asReal(eps1);
	const double ineps = asReal(eps);
	const double incrittol = asReal(crittol);
	
	
	SEXP beta;
	SEXP fperk;
	// Format output and protect them from garbage collection
	PROTECT(beta = allocVector(REALSXP, (nc+1)*nk));
	PROTECT(fperk = allocVector(REALSXP, nk));

	// Call Fortran subroutine
	F77_CALL(casso_mse_f)(REAL(x), REAL(y), INTEGER(kits), REAL(cvec), nr, nc, nk, REAL(beta), REAL(fperk),inprint, instart, inkmax,
						inmrounds, inmit, inmroundsesc, inb1, inb2, inb, inm, 
						inmclarke,inc, inrdec, inrinc, ineps1, ineps, incrittol);

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

// Define the C wrapper function for logistic family
extern SEXP c_casso_logistic_f(SEXP x, SEXP y, SEXP kits, SEXP cvec, SEXP nrow, SEXP ncol, SEXP nkits, SEXP print, SEXP start, SEXP kmax,
						SEXP mrounds, SEXP mit, SEXP mroundsesc, SEXP b1, SEXP b2, SEXP b, SEXP m, SEXP mclarke, 
						SEXP c, SEXP rdec, SEXP rinc, SEXP eps1, SEXP eps, SEXP crittol){
	// Define constants (dimensions in data / features)
	const int nr = asInteger(nrow);
	const int nc = asInteger(ncol);
	const int nk = asInteger(nkits);
	const int inkmax = asInteger(kmax);
	const int inprint = asInteger(print);
	const int instart = asInteger(start);
	const int inmrounds = asInteger(mrounds);
	const int inmit = asInteger(mit);
	const int inmroundsesc = asInteger(mroundsesc);
	const int inb1 = asInteger(b1);
	const int inb2 = asInteger(b2);
	const int inb = asInteger(b);
	const double inm = asReal(m);
	const double inmclarke = asReal(mclarke);
	const double inc = asReal(c);
	const double inrdec = asReal(rdec);
	const double inrinc = asReal(rinc);
	const double ineps1 =asReal(eps1);
	const double ineps = asReal(eps);
	const double incrittol = asReal(crittol);
	
	SEXP beta;
	SEXP fperk;
	// Format output and protect them from garbage collection
	PROTECT(beta = allocVector(REALSXP, (nc+1)*nk));
	PROTECT(fperk = allocVector(REALSXP, nk));

	
	// Call Fortran subroutine
	// Here y in {0,1}, thus INTEGER
	F77_CALL(casso_logistic_f)(REAL(x), INTEGER(y), INTEGER(kits), REAL(cvec), nr, nc, nk, REAL(beta), REAL(fperk),inprint, instart, inkmax,
						inmrounds, inmit, inmroundsesc, inb1, inb2, inb, inm, 
						inmclarke,inc, inrdec, inrinc, ineps1, ineps, incrittol);

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

// Tell R of our available Fortran functions
static const R_CallMethodDef CallEntries[] = {
  {"c_casso_cox_f",	(DL_FUNC) &c_casso_cox_f,		24},
  {"c_casso_mse_f",	(DL_FUNC) &c_casso_mse_f,		24},
  {"c_casso_logistic_f",	(DL_FUNC) &c_casso_logistic_f,		24},
  {NULL,				NULL,						0}
};


// R_init_pckgName
void R_init_casso(DllInfo *dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
