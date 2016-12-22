#ifndef  FIND_VERTICES_H
#define  FIND_VERTICES_H

Rcpp::NumericMatrix find_vertices_prob(double p, const Rcpp::NumericVector& Pi, double tol);
Rcpp::NumericMatrix find_vertices_nonneg(double p, const Rcpp::NumericVector& Pi);

#endif
