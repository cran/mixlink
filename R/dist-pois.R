#' Distribution functions
#'
#' Functions for Mixture Link Poisson distribution
#'
#' @param y Argument of pdf or cdf 
#' @param n Number of observations to draw 
#' @param mean Parameter \eqn{\vartheta} of distribution
#' @param Pi Parameter \eqn{\bm{\pi}} of distribution
#' @param kappa Parameter \eqn{\kappa} of distribution
#' @param log Return log of the result (TRUE or FALSE)
#' @param save.latent Save intermediate latent variables used during draw.
#'
#' @return \code{d.mixlink.pois} gives the density,
#'         \code{p.mixlink.pois} gives the distribution function, and
#'         \code{r.mixlink.pois} generates random deviates.
#'
#' @name Mixture Link Poisson Distribution
#'
#' @references Andrew M. Raim, Nagaraj K. Neerchal, and Jorge G. Morel.
#'             An Extension of Generalized Linear Models to Finite
#'             Mixture Outcomes. arXiv preprint: 1612.03302
#' @examples
#'   mean.true <- 20
#'   Pi.true <- c(1/4, 3/4)
#'   kappa.true <- 0.5
#'   r.mixlink.pois(n = 30, mean.true, Pi.true, kappa.true)
#'   d.mixlink.pois(y = 21, mean.true, Pi.true, kappa.true)
#'   d.mixlink.pois(y = 21, mean.true, Pi.true, kappa.true, log = TRUE)
#'   p.mixlink.pois(y = 21, mean.true, Pi.true, kappa.true)
#'
#' @name Mixture Link Poisson Distribution
r.mixlink.pois <- function(n, mean, Pi, kappa, save.latent = FALSE)
{
	if (length(mean) == 1) mean <- rep(mean, n)
	if (length(kappa) == 1) kappa <- rep(kappa, n)
	stopifnot(n == length(mean))

	J <- length(Pi)
	z <- integer(n)
	y <- numeric(n)
	psi <- matrix(NA, n, J)

	for (i in 1:n) {
		V <- find.vertices.nonneg(mean[i], Pi)
		lo <- rep(0, J)
		hi <- diag(V)
		a <- kappa[i]
		b <- kappa[i] * (J-1)
		psi[i,] <- rbeta(J, a, b)
		mu <- (hi - lo) * psi[i,] + lo
		z[i] <- sample(1:J, size = 1, prob = Pi)
		y[i] <- rpois(1, mu[z[i]])
	}

	if (save.latent) {
		return(list(y = y, z = z, psi = psi))
	} else {
		return(y)
	}
}

#' @name Find Vertices
find.vertices.nonneg <- function(mean, Pi)
{
	J <- length(Pi)
	diag(mean / Pi, J)
}

#' @name Mixture Link Poisson Distribution
d.mixlink.pois <- function(y, mean, Pi, kappa, log = FALSE)
{
	s <- max(length(y), length(mean), length(kappa))
	if (length(y) == 1) { y <- rep(y, s) }
	if (length(kappa) == 1) { kappa <- rep(kappa, s) }
	if (length(mean) == 1) { mean <- rep(mean, s) }
	ff <- .Call("d_mixlink_pois", as.integer(y), mean, Pi, kappa)
	if (log) log(ff)
	else ff
}

p.mixlink.pois.one <- function(y, mean, Pi, kappa)
{
	sum(d.mixlink.pois(0:y, mean, Pi, kappa))
}

#' @name Mixture Link Poisson Distribution
p.mixlink.pois <- Vectorize(p.mixlink.pois.one, vectorize.args = c("y", "mean"))

