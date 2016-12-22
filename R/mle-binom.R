#' Numerical MLE
#'
#' Numerical maximum likelihood estimation for Mixture Link Binomial and
#' Poisson models.
#'
#' @param y Argument of pdf or cdf.
#' @param m Number of success/failure trials.
#' @param X Design matrix for regression case.
#' @param J Number of mixture components to use.
#' @param extra.tx If additional functions of \eqn{\bm{\theta}} are to be
#'        estimated, they can be specified here. The default \code{null.tx}
#'        indicates that no extra functions are desired.
#' @param var.names A vector of strings to use for parameter names. The
#'        default (\code{NULL}) indicates to leave the names at some defaults.
#' @param phi.init Intial value of the unconstrained \eqn{\bm{\phi}}
#'        parameters. Internally, a transformation \eqn{\bm{\theta}} is applied
#'        to \eqn{\bm{\phi}} to obtain parameters in the correct space. The
#'        default (\code{NULL}) selects a default initial value.
#' @param invlink.mean The inverse link function for the mean. Default is
#'        \code{plogis} for Binomial and \code{exp} for Poisson.
#'
#' @return A list with MLE results. Can be accessed through the functions
#'         \code{confint}, \code{print}, \code{summary}, \code{coef}, \code{logLik}
#'         \code{AIC}, \code{BIC}, and \code{vcov}.
#' @name Numerical MLE
#'
#' @references Andrew M. Raim, Nagaraj K. Neerchal, and Jorge G. Morel.
#'             An Extension of Generalized Linear Models to Finite
#'             Mixture Outcomes. arXiv preprint: 1612.03302
mle.mixlink.binom <- function(y, m, J, extra.tx = null.tx,
	var.names = NULL, phi.init = NULL)
{
	if (length(m) == 1) m <- rep(m, length(y))

	Data <- list(y = y, m = m, n = length(y))
	qq <- 2 + J-1

	if (is.null(phi.init)) {
		phi.init <- c(0, mlogit(1:J), log(1))
	}

	theta.tx <- function(phi) {
		list(mean = plogis(phi[1]),
			Pi = inv.mlogit(phi[1:(J-1) + 1]),
			kappa = exp(phi[qq]))
	}

	loglik <- function(phi, Data) {
		theta <- theta.tx(phi)
		if (theta$kappa > 1e5) { return(-Inf) }
		sum(d.mixlink.binom(Data$y, Data$m, theta$mean, theta$Pi,
			theta$kappa, log = TRUE))
	}

	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)
	
	fit.out$description <- paste0(
		sprintf("y[i] ~iid~ MixLinkBinom_%d(m[i], p, Pi, kappa),\n", J)
	)
	return(fit.out)
}

#' @name Numerical MLE
mle.mixlink.binom.x <- function(y, m, X, J, extra.tx = null.tx,
	var.names = NULL, phi.init = NULL, invlink.mean = plogis)
{
	if (length(m) == 1) m <- rep(m, length(y))

	Data <- list(y = y, m = m, X = X, n = length(y))
	d <- ncol(X)
	qq <- d + J-1 + 1

	if (is.null(phi.init)) {
		phi.init <- c(rep(0,d), mlogit(1:J), log(1))
	}

	theta.tx <- function(phi) {
		list(
			Beta = phi[1:d],
			Pi = inv.mlogit(phi[1:(J-1) + d]),
			kappa = exp(phi[qq])
		)
	}
	
	loglik <- function(phi, Data)
	{
		theta <- theta.tx(phi)
		if (theta$kappa > 1e5) { return(-Inf) }
		sum(d.mixlink.binom(Data$y, Data$m, invlink.mean(X %*% theta$Beta), theta$Pi,
			theta$kappa, log = TRUE))
	}

	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)

	fit.out$description <- fit.out$description <- paste0(
		sprintf("y[i] ~indep~ MixLinkBinom_%d(m[i], p[i], Pi, kappa)\n", J),
		"link( E(Y[i]) ) = x[i]^T Beta\n"
	)
	return(fit.out)
}

