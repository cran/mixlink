#include <Rcpp.h>
#include <time.h>

char* sysdate(char* buffer)
{
	time_t t = time(0);
	struct tm* now = localtime(&t);
	sprintf(buffer, "%04d-%02d-%02d %02d:%02d:%02d", now->tm_year + 1900,
		now->tm_mon + 1, now->tm_mday, now->tm_hour, now->tm_min, now->tm_sec);
	return buffer;
}

// Avoid using RcppArmadillo just for simple matrix/vector product.
Rcpp::NumericVector matvecprod(const Rcpp::NumericMatrix& A,
	const Rcpp::NumericVector& x, const Rcpp::NumericVector& b)
{
	unsigned int m = A.nrow();
	Rcpp::NumericVector y(m);
	Rcpp::NumericMatrix AA = Rcpp::clone(A);

	for (unsigned int i = 0; i < m; i++) {
		y(i) = sum(AA.row(i) * x) + b(i);
	}

	return y;
}

Rcpp::List rwmetrop(const Rcpp::NumericVector& par_init, const Rcpp::Function& logf,
	const Rcpp::NumericMatrix& pr_var_half_trans, const Rcpp::List& grp_idx_list,
	int R, int burn, int thin, int report_period)
{
	char timestamp[50];
	int qq = par_init.size();
	int G = grp_idx_list.size();

	int R_keep = (int) ceil(double(R - burn) / double(thin));
	int idx_sample = 0;
	Rcpp::NumericMatrix par_hist(R_keep, qq);
	Rcpp::NumericVector b = par_init;
	Rcpp::NumericVector b_(qq);
	Rcpp::NumericVector accept_grp(G);
	accept_grp.fill(0);
	double logfb = Rcpp::as<double>(logf(b));

	for (int r = 0; r < R; r++) {
		const Rcpp::NumericVector& e = Rcpp::rnorm(qq, 0, 1);
		const Rcpp::NumericVector& bc = matvecprod(pr_var_half_trans, e, b);

		for (int g = 0; g < G; g++) {
			const Rcpp::IntegerVector& idx_grp = grp_idx_list(g);
			for (int j = 0; j < idx_grp.size(); j++) {
				int idx = idx_grp(j);
				b_(idx) = bc(idx);
			}

			double log_alpha = Rcpp::as<double>(logf(b_)) - logfb;
			if (std::isfinite(log_alpha)) {
				if (log(R::runif(0, 1)) < log_alpha) {
					accept_grp[g]++;

					for (int j = 0; j < idx_grp.size(); j++) {
						int idx = idx_grp(j);
						b(idx) = b_(idx);
					}					
					logfb = Rcpp::as<double>(logf(b));
				}
			} else {
				// Rcpp::warning("log_alpha was not finite");
				Rprintf("WARNING: log_alpha was not finite\n");
			}
		}

		if (r > burn-1 && r % thin == 0) {
			par_hist(idx_sample, Rcpp::_) = b;
			idx_sample++;
		}

		if (r % report_period == 0 && r > 0) {
			Rprintf("%s - After rep %d, accept%% {", sysdate(timestamp), r);
			for (int g = 0; g < G; g++) {
				Rprintf("%0.02f", accept_grp[g] / (r+1) * 100);
				if (g < G-1) {
					Rprintf(" ");
				}
			}
			Rprintf("}\n");
		}
	}

	return Rcpp::List::create(
		Rcpp::Named("par", par_hist),
		Rcpp::Named("accept", accept_grp / R)
	);
}

RcppExport SEXP rwmetrop(SEXP sexp_par_init, SEXP sexp_logf,
	SEXP sexp_pr_var_half_trans, SEXP sexp_grp_idx_list, SEXP sexp_R,
	SEXP sexp_burn, SEXP sexp_thin, SEXP sexp_report_period)
{
BEGIN_RCPP
	Rcpp::RObject __result;
	Rcpp::RNGScope __rngScope;
	Rcpp::traits::input_parameter<const Rcpp::NumericVector&>::type par_init(sexp_par_init);
	Rcpp::traits::input_parameter<const Rcpp::Function&>::type logf(sexp_logf);
	Rcpp::traits::input_parameter<const Rcpp::NumericMatrix&>::type pr_var_half_trans(sexp_pr_var_half_trans);
	Rcpp::traits::input_parameter<const Rcpp::List&>::type grp_idx_list(sexp_grp_idx_list);
	Rcpp::traits::input_parameter<int>::type R(sexp_R);
	Rcpp::traits::input_parameter<int>::type burn(sexp_burn);
	Rcpp::traits::input_parameter<int>::type thin(sexp_thin);
	Rcpp::traits::input_parameter<int>::type report_period(sexp_report_period);
	__result = Rcpp::wrap(rwmetrop(par_init, logf, pr_var_half_trans, grp_idx_list, R, burn, thin, report_period));
	return __result;
END_RCPP
}

