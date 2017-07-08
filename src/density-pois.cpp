#include <Rcpp.h>
#include <R_ext/Applic.h>
#include <gsl/gsl_sf_hyperg.h>

// [[Rcpp::export]]
Rcpp::NumericVector d_mixlink_pois(const Rcpp::IntegerVector& y,
	const Rcpp::NumericVector& mean, const Rcpp::NumericVector& Pi,
	const Rcpp::NumericVector& kappa)
{
	int n = y.size();
	int J = Pi.size();
	Rcpp::NumericVector ff(n);

	for (int i = 0; i < n; i++) {
		double F1 = 0;
		for (int j = 0; j < J; j++) {
			// Note: Huge values of kappa cause this function to take forever
			F1 += pow(Pi[j], 1-y[i]) * gsl_sf_hyperg_1F1(y[i] + kappa[i],
				y[i] + J*kappa[i], -mean[i]/Pi[j]);
		}

		double log_ff = y[i]*log(mean[i]) +
			lgamma(y[i]+kappa[i]) + lgamma(J*kappa[i]) -
			lgamma(y[i]+J*kappa[i]) - lgamma(kappa[i]) - lgamma(y[i]+1) + log(F1);
		ff[i] = exp(log_ff);
	}

	return ff;
}

