#include <Rcpp.h>
#include "find-vertices.h"

// [[Rcpp::export]]
Rcpp::NumericVector mixlink_gibbs_Q_binom(const Rcpp::IntegerVector& y,
	const Rcpp::IntegerVector& m, const Rcpp::NumericVector& mean,
	const Rcpp::NumericMatrix& psi, const Rcpp::NumericVector& Pi,
	const Rcpp::NumericVector& kappa, double find_vert_tol)
{
	Rcpp::NumericMatrix ppsi = Rcpp::clone(psi);
	int n = y.size();
	int J = Pi.size();

	Rcpp::NumericVector lo(J);
	Rcpp::NumericVector hi(J);
	Rcpp::NumericVector xi(J);
	Rcpp::NumericVector tau_sq(J);
	Rcpp::NumericVector a(J);
	Rcpp::NumericVector b(J);
	Rcpp::NumericVector mu(J);
	Rcpp::NumericVector l_beta(J);
	Rcpp::NumericVector l_binom(J);
	Rcpp::NumericVector ff(n);

	for (int i = 0; i < n; i++)
	{
		Rcpp::NumericMatrix V = find_vertices_prob(mean(i), Pi, find_vert_tol);
		int k = V.ncol();

		for (int j = 0; j < J; j++) {
			lo(j) = min(V(j, Rcpp::_));
			hi(j) = max(V(j, Rcpp::_));
			xi(j) = Rcpp::mean(V(j, Rcpp::_));
			double vv = sum(V.row(j) * V.row(j));
			double num = vv - k*xi(j)*xi(j);
			double den = k * (1 + k * kappa(i));
			tau_sq[j] = num / den;
		}

		a = 1/tau_sq * (xi-lo)*(xi-lo) * (hi-xi)/(hi-lo) - (xi-lo)/(hi-lo);
		b = a * (hi-xi)/(xi-lo);
		mu = (hi - lo) * ppsi(i, Rcpp::_) + lo;

		for (int j = 0; j < J; j++) {
			l_beta(j) = Rf_dbeta(ppsi(i,j), a(j), b(j), true);
			l_binom(j) = Rf_dbinom(y(i), m(i), mu(j), true);
		}

		ff[i] = log(sum(Pi * exp(l_binom + l_beta)));
	}

	return ff;
}

// [[Rcpp::export]]
Rcpp::NumericVector mixlink_gibbs_Q_pois(const Rcpp::IntegerVector& y,
	const Rcpp::NumericVector& mean, const Rcpp::NumericMatrix& psi,
	const Rcpp::NumericVector& Pi, const Rcpp::NumericVector& kappa)
{
	Rcpp::NumericMatrix ppsi = Rcpp::clone(psi);
	int n = y.size();
	int J = Pi.size();

	Rcpp::NumericVector hi(J);
	Rcpp::NumericVector mu(J);
	Rcpp::NumericVector l_beta(J);
	Rcpp::NumericVector l_pois(J);
	Rcpp::NumericVector ff(n);

	for (int i = 0; i < n; i++)
	{
		// Don't need to formally call
		// Rcpp::NumericMatrix V = find_vertices_nonneg(mean(i), Pi);

		hi = mean[i] / Pi;
		mu = hi * ppsi(i, Rcpp::_);

		for (int j = 0; j < J; j++) {
			l_beta(j) = Rf_dbeta(ppsi(i,j), kappa(i), kappa(i)*(J-1), true);
			l_pois(j) = Rf_dpois(y(i), mu(j), true);
		}

		ff[i] = log(sum(Pi * exp(l_pois + l_beta)));
	}

	return ff;
}

