
// includes from the plugin
#include <RcppGSL.h>
#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes

	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_blas.h>
		

// declarations
extern "C" {
SEXP file167cda520cfe( SEXP sM) ;
}

// definition

SEXP file167cda520cfe( SEXP sM ){
BEGIN_RCPP

		RcppGSL::matrix<double> M = sM; // create gsl data structures from SEXP
		int k = M.ncol();
		Rcpp::NumericVector n(k); // to store results
		for (int j = 0; j < k; j++) {
		RcppGSL::vector_view<double> colview = gsl_matrix_column (M, j);
		n[j] = gsl_blas_dnrm2(colview);
		}
		M.free() ;
		return n; // return vector
		
END_RCPP
}



