// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <gsl/gsl_sf_hyperg.h>

Rcpp::NumericVector hyperg_1F1(const Rcpp::NumericVector& a,
	const Rcpp::NumericVector& b, const Rcpp::NumericVector& x)
{
	int n = x.size();
	Rcpp::NumericVector ff(n);

	for (int i = 0; i < n; i++) {
		ff[i] = gsl_sf_hyperg_1F1(a(i), b(i), x(i));
	}

	return ff;
}
