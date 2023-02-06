#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

// family Cox
void F77_NAME(oscar_cox_f)(double *x, double *y, int *kits, double *cvec, int nrow, int ncol, int nkits, double *beta, double *fperk, int print, int start, int kmax,
							int inmrounds, int inmit, int inmroundsesc, int inb1, int inb2, int inb, double inm, double inmclarke, double inc,
							double inrdec, double inrinc, double ineps1, double ineps, double incrittol, int nKitOnes, int *betakits,
							int solver_id, int in_na, int in_mcu, int in_mcinit, double in_tolf, double in_tolf2, double in_tolg, double in_tolg2, double in_eta, double in_epsL,
			  				double in_percentage, int in_s_selection);
// family Gaussian (MSE)
void F77_NAME(oscar_mse_f)(double *x, double *y, int *kits, double *cvec, int nrow, int ncol, int nkits, double *beta, double *fperk, int print, int start, int kmax,
							int inmrounds, int inmit, int inmroundsesc, int inb1, int inb2, int inb, double inm, double inmclarke, double inc,
							double inrdec, double inrinc, double ineps1, double ineps, double incrittol, int nKitOnes, int *betakits,
							int solver_id, int in_na, int in_mcu, int in_mcinit, double in_tolf, double in_tolf2, double in_tolg, double in_tolg2, double in_eta, double in_epsL,
			  				double in_percentage, int in_s_selection);
// family Logistic
void F77_NAME(oscar_logistic_f)(double *x, int *y, int *kits, double *cvec, int nrow, int ncol, int nkits, double *beta, double *fperk, int print, int start, int kmax,
							int inmrounds, int inmit, int inmroundsesc, int inb1, int inb2, int inb, double inm, double inmclarke, double inc,
							double inrdec, double inrinc, double ineps1, double ineps, double incrittol, int nKitOnes, int *betakits,
							int solver_id, int in_na, int in_mcu, int in_mcinit, double in_tolf, double in_tolf2, double in_tolg, double in_tolg2, double in_eta, double in_epsL,
			  				double in_percentage, int in_s_selection);

// Define the C wrapper function for Cox regression
extern SEXP c_oscar_cox_f(SEXP x, SEXP y, SEXP kits, SEXP cvec, SEXP nrow, SEXP ncol, SEXP nkits, SEXP print, SEXP start, SEXP kmax,
						SEXP mrounds, SEXP mit, SEXP mroundsesc, SEXP b1, SEXP b2, SEXP b, SEXP m, SEXP mclarke,
						SEXP c, SEXP rdec, SEXP rinc, SEXP eps1, SEXP eps, SEXP crittol, SEXP nKitOnes,
						SEXP in_solver_id, SEXP na, SEXP mcu, SEXP mcinit, SEXP tolf, SEXP tolf2, SEXP tolg, SEXP tolg2, SEXP eta, SEXP epsL,
			 			SEXP percentage, SEXP s_selection){
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
	const double ineps1 = asReal(eps1);
	const double ineps = asReal(eps);
	const double incrittol = asReal(crittol);
	const int innKitOnes= asInteger(nKitOnes);
	const int solver_id = asInteger(in_solver_id);
	const int in_na = asInteger(na);
	const int in_mcu = asInteger(mcu);
	const int in_mcinit = asInteger(mcinit);
	const double in_tolf = asInteger(tolf);
	const double in_tolf2 = asInteger(tolf2);
	const double in_tolg = asInteger(tolg);
	const double in_tolg2 = asInteger(tolg2);
	const double in_eta = asInteger(eta);
	const double in_epsL = asInteger(epsL);
	const double in_percentage = asReal(percentage);
	const int in_s_selection = asInteger(s_selection);

	SEXP beta;
	SEXP fperk;
	SEXP betakits;
	// Format output and protect them from garbage collection
	PROTECT(beta = allocVector(REALSXP, nc*inkmax));
	PROTECT(fperk = allocVector(REALSXP, inkmax));
	PROTECT(betakits=allocVector(INTSXP, nk*inkmax));


	// Call Fortran subroutine
	F77_CALL(oscar_cox_f)(REAL(x), REAL(y), INTEGER(kits), REAL(cvec), nr, nc, nk, REAL(beta), REAL(fperk),inprint, instart, inkmax,
						inmrounds, inmit, inmroundsesc, inb1, inb2, inb, inm,
						inmclarke,inc, inrdec, inrinc, ineps1, ineps, incrittol, innKitOnes,INTEGER(betakits),
						solver_id, in_na, in_mcu, in_mcinit, in_tolf, in_tolf2, in_tolg, in_tolg2, in_eta, in_epsL, in_percentage,in_s_selection);

	// Create result structure
	SEXP res = PROTECT(allocVector(VECSXP,3));
	SET_VECTOR_ELT(res,0,beta);
	SET_VECTOR_ELT(res,1,fperk);
	SET_VECTOR_ELT(res,2,betakits);

	// Wrap up and return
	UNPROTECT(1);
	UNPROTECT(1);
	UNPROTECT(1);
	UNPROTECT(1);

	//return(beta);
	return(res);
}

// Define the C wrapper function for Gaussian family
extern SEXP c_oscar_mse_f(SEXP x, SEXP y, SEXP kits, SEXP cvec, SEXP nrow, SEXP ncol, SEXP nkits, SEXP print, SEXP start, SEXP kmax,
						SEXP mrounds, SEXP mit, SEXP mroundsesc, SEXP b1, SEXP b2, SEXP b, SEXP m, SEXP mclarke,
						SEXP c, SEXP rdec, SEXP rinc, SEXP eps1, SEXP eps, SEXP crittol, SEXP nKitOnes,
						SEXP in_solver_id, SEXP na, SEXP mcu, SEXP mcinit, SEXP tolf, SEXP tolf2, SEXP tolg, SEXP tolg2, SEXP eta, SEXP epsL,
			 			SEXP percentage, SEXP s_selection){
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
	const int innKitOnes = asInteger(nKitOnes);
	const int solver_id = asInteger(in_solver_id);
	const int in_na = asInteger(na);
	const int in_mcu = asInteger(mcu);
	const int in_mcinit = asInteger(mcinit);
	const double in_tolf = asInteger(tolf);
	const double in_tolf2 = asInteger(tolf2);
	const double in_tolg = asInteger(tolg);
	const double in_tolg2 = asInteger(tolg2);
	const double in_eta = asInteger(eta);
	const double in_epsL = asInteger(epsL);
	const double in_percentage = asReal(percentage);
	const int in_s_selection = asInteger(s_selection);

	SEXP beta;
	SEXP fperk;
	SEXP betakits;
	// Format output and protect them from garbage collection
	PROTECT(beta = allocVector(REALSXP, (nc+1)*nk));
	PROTECT(fperk = allocVector(REALSXP, nk));
	PROTECT(betakits=allocVector(INTSXP, nk*nk));


	// Call Fortran subroutine
	F77_CALL(oscar_mse_f)(REAL(x), REAL(y), INTEGER(kits), REAL(cvec), nr, nc, nk, REAL(beta), REAL(fperk),inprint, instart, inkmax,
						inmrounds, inmit, inmroundsesc, inb1, inb2, inb, inm,
						inmclarke,inc, inrdec, inrinc, ineps1, ineps, incrittol, innKitOnes, INTEGER(betakits),
						solver_id, in_na, in_mcu, in_mcinit, in_tolf, in_tolf2, in_tolg, in_tolg2, in_eta, in_epsL, in_percentage,in_s_selection);

	// Create result structure
	SEXP res = PROTECT(allocVector(VECSXP,3));
	SET_VECTOR_ELT(res,0,beta);
	SET_VECTOR_ELT(res,1,fperk);
	SET_VECTOR_ELT(res,2,betakits);

	// Wrap up and return
	UNPROTECT(1);
	UNPROTECT(1);
	UNPROTECT(1);
	UNPROTECT(1);
	//return(beta);
	return(res);
}

// Define the C wrapper function for logistic family
extern SEXP c_oscar_logistic_f(SEXP x, SEXP y, SEXP kits, SEXP cvec, SEXP nrow, SEXP ncol, SEXP nkits, SEXP print, SEXP start, SEXP kmax,
						SEXP mrounds, SEXP mit, SEXP mroundsesc, SEXP b1, SEXP b2, SEXP b, SEXP m, SEXP mclarke,
						SEXP c, SEXP rdec, SEXP rinc, SEXP eps1, SEXP eps, SEXP crittol, SEXP nKitOnes,
						SEXP in_solver_id, SEXP na, SEXP mcu, SEXP mcinit, SEXP tolf, SEXP tolf2, SEXP tolg, SEXP tolg2, SEXP eta, SEXP epsL,
			 			SEXP percentage, SEXP s_selection){
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
	const int innKitOnes= asInteger(nKitOnes);
	const int solver_id = asInteger(in_solver_id);
	const int in_na = asInteger(na);
	const int in_mcu = asInteger(mcu);
	const int in_mcinit = asInteger(mcinit);
	const double in_tolf = asInteger(tolf);
	const double in_tolf2 = asInteger(tolf2);
	const double in_tolg = asInteger(tolg);
	const double in_tolg2 = asInteger(tolg2);
	const double in_eta = asInteger(eta);
	const double in_epsL = asInteger(epsL);
	const double in_percentage = asReal(percentage);
	const int in_s_selection = asInteger(s_selection);

	SEXP beta;
	SEXP fperk;
	SEXP betakits;
	// Format output and protect them from garbage collection
	PROTECT(beta = allocVector(REALSXP, (nc+1)*nk));
	PROTECT(fperk = allocVector(REALSXP, nk));
	PROTECT(betakits=allocVector(INTSXP, nk*nk));



	// Call Fortran subroutine
	// Here y in {0,1}, thus INTEGER
	F77_CALL(oscar_logistic_f)(REAL(x), INTEGER(y), INTEGER(kits), REAL(cvec), nr, nc, nk, REAL(beta), REAL(fperk),inprint, instart, inkmax,
						inmrounds, inmit, inmroundsesc, inb1, inb2, inb, inm,
						inmclarke,inc, inrdec, inrinc, ineps1, ineps, incrittol, innKitOnes, INTEGER(betakits),
						solver_id, in_na, in_mcu, in_mcinit, in_tolf, in_tolf2, in_tolg, in_tolg2, in_eta, in_epsL, in_percentage,in_s_selection);

	// Create result structure
	SEXP res = PROTECT(allocVector(VECSXP,3));
	SET_VECTOR_ELT(res,0,beta);
	SET_VECTOR_ELT(res,1,fperk);
	SET_VECTOR_ELT(res,2,betakits);

	// Wrap up and return
	UNPROTECT(1);
	UNPROTECT(1);
	UNPROTECT(1);
	UNPROTECT(1);

	//return(beta);
	return(res);
}

// Tell R of our available Fortran functions
static const R_CallMethodDef CallEntries[] = {
  {"c_oscar_cox_f",	(DL_FUNC) &c_oscar_cox_f,		37},
  {"c_oscar_mse_f",	(DL_FUNC) &c_oscar_mse_f,		37},
  {"c_oscar_logistic_f",	(DL_FUNC) &c_oscar_logistic_f,		37},
  {NULL,				NULL,						0}
};


// R_init_pckgName
void R_init_oscar(DllInfo *dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
