#' @name Numerical MLE
#' @examples
#' \dontrun{
#'   n <- 400
#'   mean.true <- rep(20, n)
#'   Pi.true <- c(1/5, 4/5)
#'   kappa.true <- 2
#'   y <- r.mixlink.pois(n, mean.true, Pi.true, kappa.true)
#'
#'   mle.out <- mle.mixlink.pois(y, J = 2)
#'   coef(mle.out)
#'   print(mle.out)
#'   confint(mle.out)
#' }
mle.mixlink.pois <- function(y, J, extra.tx = null.tx,
	var.names = NULL, phi.init = NULL)
{
	Data <- list(y = y, n = length(y))
	qq <- 2 + J-1

	if (is.null(phi.init)) {
		phi.init <- c(0, mlogit(1:J), log(1))
	}

	theta.tx <- function(phi) {
		list(mean = exp(phi[1]),
			Pi = inv.mlogit(phi[1:(J-1) + 1]),
			kappa = exp(phi[qq]))
	}

	loglik <- function(phi, Data) {
		theta <- theta.tx(phi)
		if (theta$kappa > 1e5) { return(-Inf) }
		sum(d.mixlink.pois(Data$y, theta$mean, theta$Pi, theta$kappa, log = TRUE))
	}

	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)
	
	fit.out$description <- paste0(
		sprintf("y[i] ~iid~ MixLinkPois_%d(mean, Pi, kappa),\n", J)
	)
	return(fit.out)
}

#' @name Numerical MLE
mle.mixlink.pois.x <- function(y, X, J, extra.tx = null.tx,
	var.names = NULL, phi.init = NULL, invlink.mean = exp)
{
	Data <- list(y = y, X = X, n = length(y))
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

	loglik <- function(phi, Data) {
		theta <- theta.tx(phi)
		if (theta$kappa > 1e5) { return(-Inf) }
		sum(d.mixlink.pois(Data$y, invlink.mean(X %*% theta$Beta), theta$Pi,
			theta$kappa, log = TRUE))
	}

	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)

	fit.out$description <- fit.out$description <- paste0(
		sprintf("y[i] ~indep~ MixLinkPois_%d(mean[i], Pi, kappa)\n", J),
		"link( E(Y[i]) ) = x[i]^T Beta\n"
	)
	return(fit.out)
}
