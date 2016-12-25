#include <Rcpp.h>
#include <vector>
#include <list>
#include <cmath>
#include <bitset>

typedef std::list< std::vector<double> > vector_list_t;

// From: http://stackoverflow.com/questions/16761472/how-can-i-increment-stdbitset
// void increment(std::vector<bool>& in)
void increment(std::bitset<128>& in)
{
	for (size_t i = 0; i < in.size(); ++i)
	{
		if (in[i] == 0) {
			in[i] = 1;
			break;
		}
		in[i] = 0;
	}
}

void list_vertices(vector_list_t& vert, double p, const Rcpp::NumericVector Pi)
{
	int J = Pi.size();
	int r = 0;
	unsigned int bin_cnt = (unsigned int) pow(2.0, (double)(J-1) );
	int q = 0;

	for (int j = 0; j < J; j++)
	{
		// If Pi[j] == 0, we can't divide by it and get a solution for v.j 
		if (Pi(j) == 0)
			continue;

		// x <- as.vector( digitsBase(k, base = 2, J-1) )
		std::bitset<128> x(0);
		for ( ; x.to_ulong() < bin_cnt; increment(x))
		{
			// ss = t(x) %*% Pi[-j])
			double ss = 0;
			q = 0;
			for (int l = 0; l < J; l++) {
				if (l != j) {
					// If the q-th digit (base 2) of x is == 1, add Pi[l] to ss
					ss += Pi(l) * x[q];
					q++;
				}
			}

			// v.j <- 1/Pi[j] * (p - t(x) %*% Pi[-j])
			double v_j = 1/Pi(j) * (p - ss);

			if (0 <= v_j && v_j <= 1)
			{
				r++;
				std::vector<double> v(J);
			
				// v[-j] <- x
				// v[j] <- v.j

				q = 0;
				for (int l = 0; l < J; l++) {
					if (l == j) {
						v[l] = v_j;
					} else {
						v[l] = x[q];
						q++;
					}
				}

				vert.push_back(v);
			}
		}
	}
}

double dist2(const std::vector<double>& x, const std::vector<double>& y)
{
	double ss = 0;
	unsigned int n = x.size();

	if (n != y.size()) {
		throw std::range_error("In dist2, vectors must have same length");
	}

	for (unsigned int i = 0; i < n; i++) {
		ss += pow(x[i] - y[i], 2.0);
	}

	return sqrt(ss);
}

void extract_unique_vertices(vector_list_t& unique_vert_list,
	const vector_list_t& vert_list, double tol)
{
	for (vector_list_t::const_iterator itr = vert_list.begin();
		itr != vert_list.end(); itr++)
	{
		// Compare the vertex at itr to all previous vertices
		// in the unique list. If norm distance is very small, count it
		// as a match
		bool is_dup = FALSE;

		for (vector_list_t::const_iterator uniq_itr = unique_vert_list.begin();
			uniq_itr != unique_vert_list.end(); uniq_itr++)
		{
			double d = dist2(*itr, *uniq_itr);
			if (d < tol) {
				is_dup = TRUE;
				break;
			}
		}

		if (!is_dup) {
			unique_vert_list.push_back(*itr);
		}
	}
}

Rcpp::NumericMatrix vector_list_to_matrix(const vector_list_t& vert_list,
	int J, int k)
{
	Rcpp::NumericMatrix V(J, k);

	int l = 0;
	for (vector_list_t::const_iterator itr = vert_list.begin();
		itr != vert_list.end(); itr++)
	{
		const std::vector<double>& v = *itr;
		for (int j = 0; j < J; j++) {
			V(j, l) = v[j];
		}

		l++;
	}

	return V;
}

Rcpp::NumericMatrix find_vertices_prob(double p, const Rcpp::NumericVector& Pi, double tol)
{
	int J = Pi.size();
	if (J > 128) {
		// This is because we use std::bitset whose length must be set at
		// compile time. Calculations become too intensive long beyond this
		// point, so this is not a major restriction.
		throw std::range_error("In find_vertices_prob, only J <= 128 is supported");
	}

	// Call the C++ function to do the real work
	vector_list_t vert_list;
	vector_list_t unique_vert_list;

	list_vertices(vert_list, p, Pi);
	extract_unique_vertices(unique_vert_list, vert_list, tol);
	int k = unique_vert_list.size();

	// Now pack up the list of vertices into a J x k matrix and return it
	return vector_list_to_matrix(unique_vert_list, J, k);
}

Rcpp::NumericMatrix find_vertices_nonneg(double p, const Rcpp::NumericVector& Pi)
{
	int J = Pi.size();
	Rcpp::NumericMatrix V(J, J);
	
	for (int j = 0; j < J; j++) {
		V(j,j) = p / Pi[j];
	}

	return V;
}

