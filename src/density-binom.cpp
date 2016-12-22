#include <Rcpp.h>
#include <R_ext/Applic.h>
#include "find-vertices.h"

// For info on calling R integration methods from C:
// https://cran.r-project.org/doc/manuals/r-release/R-exts.html
// http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2014-July/007872.html

typedef struct {
	double C;
	int y;
	int m;
	double lo;
	double hi;
	double a;
	double b;
} Params_binom;

void integrand_binom(double *w, int n, void *ex)
{
	Params_binom* par = (Params_binom*) ex;
	double C = par->C;
	int y = par->y;
	int m = par->m;
	double lo = par->lo;
	double hi = par->hi;
	double a = par->a;
	double b = par->b;

	for (int i = 0; i < n; i++) {
		double mu = (hi - lo) * w[i] + lo;
		double l1 = C + y*log(mu) + (m-y)*log(1-mu);
		double l2 = R::dbeta(w[i], a, b, true);
		w[i] = exp(l1 + l2);
	}
}

Rcpp::NumericVector d_mixlink_binom(const Rcpp::IntegerVector& y,
	const Rcpp::IntegerVector& m, const Rcpp::NumericVector& mean,
	const Rcpp::NumericVector& Pi, const Rcpp::NumericVector& kappa,
	int subdiv, double rel_tol, double abs_tol)
{
	int n = y.size();
	int J = Pi.size();
	Rcpp::NumericVector ff(n);
	Rcpp::NumericVector f_i(J);
	Rcpp::NumericVector xi(J);
	Rcpp::NumericVector tau_sq(J);
	Rcpp::NumericVector f_(J);
	Rcpp::NumericVector lo(J);
	Rcpp::NumericVector hi(J);

	// Some inputs to numerical integration
	int lenw = 4 * subdiv;
	double integ_lower = 0;
	double integ_upper = 1;

	// CRAN compilers give warnings about these, so use malloc/free instead
	// int iwork[subdiv];
	// double work[lenw];

	int* iwork = (int*) malloc(subdiv * sizeof(int));
	double* work = (double*) malloc(lenw * sizeof(double));

	// Outputs of numerical integration
	double abserr;
	int last;
	int ier;
	int neval;

	for (int i = 0; i < n; i++)
	{
		// const Rcpp::NumericMatrix& V = find_vertices_prob(mean[i], Pi, 1e-8);
		Rcpp::NumericMatrix V = find_vertices_prob(mean[i], Pi, 1e-8);
		int k = V.ncol();

		// E() and Var() of lin comb of Dirichlet	
		for (int j = 0; j < J; j++) {
			lo[j] = Rcpp::min(V.row(j));
			hi[j] = Rcpp::max(V.row(j));
			xi[j] = Rcpp::mean(V.row(j));
			double ss_j = Rcpp::sum(V.row(j) * V.row(j));
			double num = ss_j - k*pow(xi[j], 2.0);
			double denom = k * (1 + k * kappa[i]);
			tau_sq[j] = num / denom;
		}

		// Equate beta moments to xi and tau.sq for j = 1, ..., J
		const Rcpp::NumericVector& a = 1/tau_sq * pow((xi - lo), 2.0) * (hi-xi)/(hi-lo) -
			(xi-lo)/(hi-lo);
		const Rcpp::NumericVector& b = a * (hi-xi)/(xi-lo);

		double C = lgamma(m[i]+1) - lgamma(y[i]+1) - lgamma(m[i]-y[i]+1);
		const Rcpp::LogicalVector& is_degenerate = (hi-lo < 1e-20);

		for (int j = 0; j < J; j++) {
			// If hi[j] == lo[j], integration is over a singleton set
			if (is_degenerate[j]) {
				f_i[j] = exp( C + y[i]*log(hi[j]) + (m[i]-y[i])*log(1-hi[j]) );
			} else {
				Params_binom par;
				par.C = C;
				par.y = y[i];
				par.m = m[i];
				par.lo = lo[j];
				par.hi = hi[j];
				par.a = a[j];
				par.b = b[j];

				Rdqags(integrand_binom, (void*) &par, &integ_lower, &integ_upper,
					&abs_tol, &rel_tol, &f_i[j], &abserr, &neval, &ier,
					&subdiv, &lenw, &last, iwork, work);
			}
		}

		ff[i] = sum(Pi * f_i);
	}

	free(iwork);
	free(work);

	return ff;
}

