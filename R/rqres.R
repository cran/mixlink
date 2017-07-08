#' Randomized quantile residuals for Mixture Link
#'
#' Compute randomized quantile residuals for the Mixture Link Binomial
#' and Mixture Link Poisson distributions.
#'
#' @param y The observations
#' @param m Number of success/failure trials
#' @param mean Estimate for parameter \eqn{\vartheta} of distribution
#' @param Pi Estimate for parameter \eqn{\bm{\pi}} of distribution
#' @param kappa Estimate for parameter \eqn{\kappa} of distribution
#'
#' @return Vector of residuals
#'
#' @references Peter K. Dunn and Gordon K. Smyth. Randomized quantile
#'             residuals. Journal of Computational and Graphical Statistics,
#'             5(3):236-244, 1996.
#' @name Randomized quantile residuals
#' @examples
#' n <- 400
#' mean.true <- rep(20, n)
#' Pi.true <- c(1/4, 3/4)
#' kappa.true <- 1.5
#' y <- r.mixlink.pois(n, mean.true, Pi.true, kappa.true)
#' r <- rqres.mixlink.pois(y, mean.true, Pi.true, kappa.true)
#' qqnorm(r); qqline(r)
#'
NULL

#  Set eps to zero to avoid using random jitter
rqres <- function(y, F, eps = 1e-6)
{
	n <- length(y)
	FL <- pmin(pmax(F(y - eps), 0), 1)
	FU <- pmax(pmin(F(y), 1), 0)
	idx.neq <- which(FL < FU)
	u <- FL
	u[idx.neq] <- runif(length(idx.neq), min = FL, max = FU)
	qres <- qnorm(u)
	return(qres)
}

rqres.mixlink.binom.one <- function(y, m, mean, Pi, kappa)
{
	F <- function(y) {
		p.mixlink.binom.one(y, m, mean, Pi, kappa)
	}
	rqres(y, F)
}

#' @name Randomized quantile residuals
rqres.mixlink.binom <- Vectorize(rqres.mixlink.binom.one, vectorize.args = c("y", "m", "mean"))

rqres.mixlink.pois.one <- function(y, mean, Pi, kappa)
{
	F <- function(y) {
		p.mixlink.pois.one(y, mean, Pi, kappa)
	}
	rqres(y, F)
}

#' @name Randomized quantile residuals
rqres.mixlink.pois <- Vectorize(rqres.mixlink.pois.one, vectorize.args = c("y", "mean"))

