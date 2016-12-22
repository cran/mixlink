#include <Rcpp.h>

using namespace Rcpp;

Rcpp::NumericVector hyperg_1F1(const Rcpp::NumericVector& a,
    const Rcpp::NumericVector& b, const Rcpp::NumericVector& x);
RcppExport SEXP hyperg_1F1(SEXP sexp_a, SEXP sexp_b, SEXP sexp_x)
{
BEGIN_RCPP
	Rcpp::RObject __result;
	Rcpp::RNGScope __rngScope;
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type a(sexp_a);
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type b(sexp_b);
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type x(sexp_x);
	__result = Rcpp::wrap(hyperg_1F1(a, b, x));
	return __result;
END_RCPP
}

Rcpp::NumericMatrix find_vertices_prob(double p, const Rcpp::NumericVector& Pi, double tol);
RcppExport SEXP find_vertices_prob(SEXP sexp_p, SEXP sexp_Pi, SEXP sexp_tol)
{
BEGIN_RCPP
	Rcpp::RObject __result;
	Rcpp::RNGScope __rngScope;
	Rcpp::traits::input_parameter<double>::type p(sexp_p);
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type Pi(sexp_Pi);
	Rcpp::traits::input_parameter<double>::type tol(sexp_tol);
	__result = Rcpp::wrap(find_vertices_prob(p, Pi, tol));
	return __result;
END_RCPP
}

Rcpp::NumericMatrix find_vertices_nonneg(double p, const Rcpp::NumericVector& Pi);
RcppExport SEXP find_vertices_nonneg(SEXP sexp_p, SEXP sexp_Pi)
{
BEGIN_RCPP
	Rcpp::RObject __result;
	Rcpp::RNGScope __rngScope;
	Rcpp::traits::input_parameter<double>::type p(sexp_p);
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type Pi(sexp_Pi);
	__result = Rcpp::wrap(find_vertices_nonneg(p, Pi));
	return __result;
END_RCPP
}

Rcpp::NumericVector d_mixlink_binom(const Rcpp::IntegerVector& y,
	const Rcpp::IntegerVector& m, const Rcpp::NumericVector& mean,
	const Rcpp::NumericVector& Pi, const Rcpp::NumericVector& kappa,
	int subdiv, double rel_tol, double abs_tol);
RcppExport SEXP d_mixlink_binom(SEXP sexp_y, SEXP sexp_m, SEXP sexp_mean,
	SEXP sexp_Pi, SEXP sexp_kappa, SEXP sexp_subdiv, SEXP sexp_rel_tol,
	SEXP sexp_abs_tol)
{
BEGIN_RCPP
	Rcpp::RObject __result;
	Rcpp::RNGScope __rngScope;
	Rcpp::traits::input_parameter<const Rcpp::IntegerVector&>::type y(sexp_y);
	Rcpp::traits::input_parameter<const Rcpp::IntegerVector&>::type m(sexp_m);
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type mean(sexp_mean);
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type Pi(sexp_Pi);
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type kappa(sexp_kappa);
	Rcpp::traits::input_parameter<int>::type subdiv(sexp_subdiv);
	Rcpp::traits::input_parameter<double>::type rel_tol(sexp_rel_tol);
	Rcpp::traits::input_parameter<double>::type abs_tol(sexp_abs_tol);
	__result = Rcpp::wrap(d_mixlink_binom(y, m, mean, Pi, kappa, subdiv, rel_tol, abs_tol));
	return __result;
END_RCPP
}

Rcpp::NumericVector d_mixlink_pois(const Rcpp::IntegerVector& y,
	const Rcpp::NumericVector& mean, const Rcpp::NumericVector& Pi,
	const Rcpp::NumericVector& kappa);
RcppExport SEXP d_mixlink_pois(SEXP sexp_y, SEXP sexp_mean, SEXP sexp_Pi,
	SEXP sexp_kappa)
{
BEGIN_RCPP
	Rcpp::RObject __result;
	Rcpp::RNGScope __rngScope;
	Rcpp::traits::input_parameter<const Rcpp::IntegerVector&>::type y(sexp_y);
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type mean(sexp_mean);
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type Pi(sexp_Pi);
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type kappa(sexp_kappa);
	__result = Rcpp::wrap(d_mixlink_pois(y, mean, Pi, kappa));
	return __result;
END_RCPP
}

Rcpp::NumericVector mixlink_gibbs_Q_binom(const Rcpp::IntegerVector& y,
	const Rcpp::IntegerVector& m, const Rcpp::NumericVector& mean,
	const Rcpp::NumericMatrix& psi, const Rcpp::NumericVector& Pi,
	const Rcpp::NumericVector& kappa, double find_vert_tol);
RcppExport SEXP mixlink_gibbs_Q_binom(SEXP sexp_y, SEXP sexp_m, SEXP sexp_mean,
	SEXP sexp_psi, SEXP sexp_Pi, SEXP sexp_kappa, SEXP sexp_find_vert_tol)
{
BEGIN_RCPP
	Rcpp::RObject __result;
	Rcpp::RNGScope __rngScope;
	Rcpp::traits::input_parameter<const Rcpp::IntegerVector&>::type y(sexp_y);
	Rcpp::traits::input_parameter<const Rcpp::IntegerVector&>::type m(sexp_m);
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type mean(sexp_mean);
	Rcpp::traits::input_parameter<const Rcpp::NumericMatrix&>::type psi(sexp_psi);
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type Pi(sexp_Pi);
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type kappa(sexp_kappa);
	Rcpp::traits::input_parameter<double>::type find_vert_tol(sexp_find_vert_tol);
	__result = Rcpp::wrap(mixlink_gibbs_Q_binom(y, m, mean, psi, Pi, kappa, find_vert_tol));
	return __result;
END_RCPP
}

Rcpp::NumericVector mixlink_gibbs_Q_pois(const Rcpp::IntegerVector& y,
	const Rcpp::NumericVector& mean, const Rcpp::NumericMatrix& psi,
	const Rcpp::NumericVector& Pi, const Rcpp::NumericVector& kappa);
RcppExport SEXP mixlink_gibbs_Q_pois(SEXP sexp_y, SEXP sexp_mean, SEXP sexp_psi,
	SEXP sexp_Pi, SEXP sexp_kappa)
{
BEGIN_RCPP
	Rcpp::RObject __result;
	Rcpp::RNGScope __rngScope;
	Rcpp::traits::input_parameter<const Rcpp::IntegerVector&>::type y(sexp_y);
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type mean(sexp_mean);
	Rcpp::traits::input_parameter<const Rcpp::NumericMatrix&>::type psi(sexp_psi);
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type Pi(sexp_Pi);
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type kappa(sexp_kappa);
	__result = Rcpp::wrap(mixlink_gibbs_Q_pois(y, mean, psi, Pi, kappa));
	return __result;
END_RCPP
}

