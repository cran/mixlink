#' Distribution functions
#'
#' Functions for Mixture Link Binomial distribution
#'
#' @param y Argument of pdf or cdf 
#' @param n Number of observations to draw 
#' @param m Number of success/failure trials
#' @param mean Parameter \eqn{\vartheta} of distribution
#' @param Pi Parameter \eqn{\bm{\pi}} of distribution
#' @param kappa Parameter \eqn{\kappa} of distribution
#' @param log Return log of the result (TRUE or FALSE)
#' @param save.latent Save intermediate latent variables used during draw.
#'
#' @return \code{d.mixlink.binom} gives the density,
#'         \code{p.mixlink.binom} gives the distribution function, and
#'         \code{r.mixlink.binom} generates random deviates.
#'
#' @name Mixture Link Binomial Distribution
#'
#' @references Andrew M. Raim, Nagaraj K. Neerchal, and Jorge G. Morel.
#'             An Extension of Generalized Linear Models to Finite
#'             Mixture Outcomes. arXiv preprint: 1612.03302
#' @examples
#'   mean.true <- 1/3
#'   Pi.true <- c(1/5, 4/5)
#'   kappa.true <- 0.5
#'   m <- 10
#'   r.mixlink.binom(n = 30, mean.true, Pi.true, kappa.true, m)
#'   d.mixlink.binom(y = 5, m, mean.true, Pi.true, kappa.true)
#'   d.mixlink.binom(y = 5, m, mean.true, Pi.true, kappa.true, log = TRUE)
#'   p.mixlink.binom(y = 5, m, mean.true, Pi.true, kappa.true)
#' @name Mixture Link Binomial Distribution
r.mixlink.binom <- function(n, mean, Pi, kappa, m, save.latent = FALSE)
{
	if (length(mean) == 1) mean <- rep(mean, n)
	if (length(m) == 1) m <- rep(m, n)
	stopifnot(n == length(mean))
	stopifnot(n == length(m))

	J <- length(Pi)
	z <- integer(n)
	y <- numeric(n)
	psi <- matrix(NA, n, J)

	for (i in 1:n) {
		V <- find.vertices.prob(mean[i], Pi)
		k <- ncol(V)
		lo <- apply(V, 1, min)
		hi <- apply(V, 1, max)
		xi <- rowMeans(V)
		tau.sq <- numeric(J)
		for (j in 1:J) {
			num <- k * t(V[j,]) %*% V[j,] - (k*xi[j])^2
			denom <- k^2 * (1 + k * kappa)
			tau.sq[j] <- num / denom
		}
		a <- 1/tau.sq * (xi - lo)^2 * (hi-xi)/(hi-lo) - (xi-lo)/(hi-lo)
		b <- a * (hi-xi)/(xi-lo)
		psi[i,] <- rbeta(J, a, b)
		mu <- (hi - lo) * psi[i,] + lo
		z[i] <- sample(1:J, size = 1, prob = Pi)
		y[i] <- rbinom(1, m[i], mu[z[i]])
	}
	
	if (save.latent) {
		return(list(y = y, z = z, psi = psi))
	} else {
		return(y)
	}
}

#' Compute vertices for Mixture Link
#'
#' Find vertices of the set \eqn{A(\vartheta, \bm{\pi})}, which characterizes
#' link between finite mixture mean and regression function.
#'
#' @details For Mixture Link Binomial, the set \eqn{A(\vartheta, \bm{\pi}) =
#' \{ \bm{\mu} \in [0,1]^J : \bm{\mu}^T \bm{\pi} = \vartheta \}. }
#' For Mixture Link Poisson, the set \eqn{A(\vartheta, \bm{\pi}) =
#' \{ \bm{\mu} \in [0,\infty]^J : \bm{\mu}^T \bm{\pi} = \vartheta \}. }
#'
#' @param mean Parameter \eqn{\vartheta} of distribution.
#' @param Pi Parameter \eqn{\bm{\pi}} of distribution.
#' @param tol A tolerance to determine if candidate vertices are distinct.
#'
#' @return A \eqn{J \times k} matrix whose columns are the vertices of
#'         \eqn{A(\vartheta, \bm{\pi})}.
#'
#' @references Andrew M. Raim, Nagaraj K. Neerchal, and Jorge G. Morel.
#'             An Extension of Generalized Linear Models to Finite
#'             Mixture Outcomes. arXiv preprint: 1612.03302
#' @name Find Vertices
#'
find.vertices.prob <- function(mean, Pi, tol = 1e-8)
{
	vert <- .Call("find_vertices_prob", mean, Pi, tol)
	return(vert)
}

#' @name Mixture Link Binomial Distribution
d.mixlink.binom <- function(y, m, mean, Pi, kappa, log = FALSE)
{
	s <- max(length(y), length(m), length(mean), length(kappa))
	if (length(y) == 1) { y <- rep(y, s) }
	if (length(m) == 1) { m <- rep(m, s) }
	if (length(kappa) == 1) { kappa <- rep(kappa, s) }
	if (length(mean) == 1) { mean <- rep(mean, s) }

	# Set integration parameters to defaults as if calling integrate directly
	subdiv <- 100
	rel.tol <- .Machine$double.eps^.25
	abs.tol <- rel.tol

	ff <- .Call("d_mixlink_binom", as.integer(y), as.integer(m), mean, Pi,
		kappa, subdiv, rel.tol, abs.tol)
	if (log) log(ff)
	else ff
}

p.mixlink.binom.one <- function(y, m, mean, Pi, kappa)
{
	sum(d.mixlink.binom(0:y, m, mean, Pi, kappa))
}

#' @name Mixture Link Binomial Distribution
p.mixlink.binom <- Vectorize(p.mixlink.binom.one, vectorize.args = c("y", "m", "mean"))

