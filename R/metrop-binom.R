#' Metropolis-Hastings Sampler
#'
#' Random-Walk Metropolis Hastings sampler for Binomial and Poisson Mixture
#' Link models.
#'
#' @details Priors for Bayesian Mixture Link model are
#' \itemize{
#' \item{\eqn{ \bm{\beta} \sim \textrm{N}(\bm{0}, V_{\beta} \bm{I}) }},
#' \item{\eqn{ \bm{\pi} \sim \textrm{Dirichlet}_J(\bm{\alpha}_\pi) }},
#' \item{\eqn{ \kappa \sim \textrm{Gamma}(a_\kappa, b_\kappa) }}, parameterized with 
#'       \eqn{\textrm{E}(\kappa) = a_\kappa / b_\kappa}.
#' }
#'
#' @param y Observations.
#' @param X Design matrix for regression.
#' @param m Number of success/failure trials.
#' @param R Number of MCMC draws to take.
#' @param burn Number of initial MCMC draws to discard.
#' @param thin After burn-in period, save one of every \code{thin} draws.
#' @param invlink The inverse link function for the mean. Default is
#'        \code{plogis} for Binomial and \code{exp} for Poisson.
#' @param Beta.init Starting value for \eqn{\bm{\beta}}.
#' @param Pi.init Starting value for \eqn{\bm{\pi}}.
#' @param kappa.init Starting value for \eqn{\kappa}.
#' @param hyper A list with hyperparameters corresponding to the prior from the
#'        Details section. \code{var.Beta} is \eqn{V_{\beta}} with default 1000.
#'        \code{alpha.Pi} is \eqn{\bm{\alpha}_\pi} with default \code{rep(1,J)}.
#'        \code{a.kappa} and \code{b.kappa} correspond to \eqn{(a_\kappa, b_\kappa)}
#'        which have default values \code{a.kappa = 1} and \code{b.kappa = 1/10}.
#' @param reporting.period Report progress every \code{reporting.period} draws.
#' @param use.laplace.approx Maximize a Laplace approximation to the posterior,
#'        to find a starting value for MCMC.
#' @param proposal A list with two elements. \code{var} is the covariance matrix
#'        for a \eqn{d + (J-1) + 1} dimensional multivariate normal proposal, \code{scale}
#'        is a scalar multiplied with \code{var}. Defaults are \code{var = diag(d+(J-1)+1)}
#'        and \code{scale = 0.02}.
#' @param param.grp
#'        A vector of integers of length \eqn{d + (J-1) + 1}, where \code{d = ncol(X)},
#'        which indicates the grouping of parameters in MCMC. Parameters with common
#'        integers are sampled together. The first \eqn{d} correspond to \eqn{\bm{\beta}},
#'        the next \eqn{J-1} correspond to \eqn{\bm{\pi}}, and the last one corresponds
#'        to \eqn{\kappa}. At the default value (\code{param.grp = NULL}),
#'        \eqn{\bm{\beta}}, \eqn{\bm{\pi}}, and \eqn{\kappa} are sampled one at a time,
#'        each in their entirety.
#' @param fixed.kappa Keep \eqn{\kappa} fixed at \code{kappa.init} (Default FALSE).
#'
#' @return A list with the MCMC results:
#' \item{par.hist}{\eqn{R \times [d + (J-1) + 1]} matrix of saved MCMC draws
#'    before transformations are applied. Most users will not need this.}
#' \item{Beta.hist}{\eqn{R \times d} matrix of saved \eqn{\bm{\beta}} draws}
#' \item{Pi.hist}{\eqn{R \times J} matrix of saved \eqn{\bm{\pi}} draws}
#' \item{kappa.hist}{\eqn{R \times 1} vector of \eqn{\kappa} draws}
#' \item{accept}{Percentages that MCMC proposals were accepted. Corresponds to
#'   \code{param.grp} }
#' \item{elapsed.sec}{Elapsed time for sampling, in seconds.}
#' \item{laplace.out}{Output of Laplace approximation.}
#' \item{R.keep}{Number of draws kept, after thinning and burn-in.}
#' \item{X.names}{Names of columns of design matrix.}
#' Can be accessed with the functions \code{print}, \code{summary}, and \code{DIC}.
#' @name Metropolis-Hastings Sampler
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#'
#' # ----- Generate data -----
#' n <- 200
#' x <- runif(n, 1, 3)
#' X <- model.matrix(~ x)
#' Beta.true <- c(0, 1)
#' mean.true <- exp(X %*% Beta.true)
#' kappa.true <- 0.95
#' Pi.true <- c(1,3)/4
#' d <- ncol(X)
#' J <- length(Pi.true)
#' y <- r.mixlink.pois(n, mean.true, Pi.true, kappa.true)
#' 
#' # ----- Run Metropolis-within-Gibbs sampler -----
#' hyper <- list(VBeta = diag(1000, d), alpha.Pi = rep(1, J),
#' 	kappa.a = 1, kappa.b = 1/2)
#' proposal <- list(
#'	var = bdiag(solve(t(X) %*% X), diag(J-1), 1),
#'	scale = 0.5)
#' metrop.out <- rwmetrop.mixlink.poisson(y, X, R = 20000, burn = 1000,
#' 	thin = 10, Pi.init = c(1,9)/10, hyper = hyper,
#' 	reporting.period = 1000, use.laplace.approx = TRUE, proposal = proposal)
#' 
#' print(metrop.out)
#' DIC.rwmetrop.mixlink.poisson(metrop.out, y, X, invlink = exp)
#' }
#'
#' @references Andrew M. Raim, Nagaraj K. Neerchal, and Jorge G. Morel.
#'             An Extension of Generalized Linear Models to Finite
#'             Mixture Outcomes. arXiv preprint: 1612.03302
#'
rwmetrop.mixlink.binomial <- function(y, X, m, R, burn = 0, thin = 1,
	invlink = plogis, Beta.init = NULL, Pi.init, kappa.init = NULL,
	hyper = NULL, reporting.period = 1, use.laplace.approx = TRUE, proposal = NULL,
	param.grp = NULL, fixed.kappa = FALSE)
{
	st <- Sys.time()
	d <- ncol(X)
	J <- length(Pi.init)
	qq <- d + J-1 + 1
	
	if (is.null(param.grp)) {
		param.grp <- c(rep(1,d), rep(2,J-1), 3)
	}

	if (is.null(Beta.init)) { Beta.init <- rep(0,d) }
	if (is.null(kappa.init)) { kappa.init <- 1 }

	if (is.null(proposal)) { proposal <- list() }
	if (is.null(proposal$var)) { proposal$var <- diag(qq) }
	if (is.null(proposal$scale)) { proposal$scale <- 0.02 }

	if (is.null(hyper)) { hyper <- list() }
	if (is.null(hyper$var.Beta)) { hyper$var.Beta <- 1000 }
	if (is.null(hyper$alpha.Pi)) { hyper$alpha.Pi <- rep(1,J) }
	if (is.null(hyper$a.kappa)) { hyper$a.kappa <- 1 }
	if (is.null(hyper$b.kappa)) { hyper$b.kappa <- 1/10 }

	Data <- list(y = y, m = m, X = X, J = J)

	logpost <- function(parm, Data)
	{
		n <- nrow(Data$X)
		d <- ncol(Data$X)
		J <- Data$J

		Beta <- parm[1:d]
		Pi <- inv.mlogit(parm[1:(J-1) + d])

		Beta.prior <- sum(dnorm(Beta, 0, sqrt(hyper$var.Beta), log = TRUE))
		Pi.prior <- ddirichlet(Pi, hyper$alpha.Pi, log = TRUE)

		jac.Pi <- jacobian(func = inv.mlogit, x = mlogit(Pi))
		l.Pi.tx <- log(abs(det(as.matrix(jac.Pi[-J,]))))

		if (fixed.kappa) {
			kappa <- kappa.init
			kappa.prior <- 1
			l.kappa.tx <- 0
		} else {
			kappa <- exp(parm[1 + d + J-1])
			kappa.prior <- dgamma(kappa, hyper$a.kappa, hyper$b.kappa, log = TRUE)
			l.kappa.tx <- log(kappa)	
		}

		p <- invlink(Data$X %*% Beta)
		ll <- d.mixlink.binom(Data$y, Data$m, p, Pi, kappa, log = TRUE)

		sum(ll) + Beta.prior + Pi.prior + kappa.prior + l.Pi.tx + l.kappa.tx
	}

	start <- c(Beta.init, mlogit(Pi.init), log(kappa.init))

	# ----- Laplace Approximation -----
	if (use.laplace.approx) {
		logger("Computing Laplace approximation to find MCMC starting value\n")
		optim.control <- list(trace = 6)
		laplace.out <- laplace(logpost, start, Data, optim.control = optim.control,
			optim.method = "Nelder-Mead")
		start <- laplace.out$mode
	} else {
		laplace.out <- NULL
	}

	# ----- Fit the model -----
	logger("Starting MCMC\n")
	metrop.out <- rwmetrop(start, logpost, Data, proposal, grp = param.grp,
		R = R, burn = burn, thin = thin, report.period = reporting.period)
	logger("Finished MCMC\n")

	elapsed.sec <- as.numeric(Sys.time() - st, unit = "secs")

	idx.Beta <- seq(1, d)
	idx.Pi <- seq(1, J-1) + d
	idx.kappa <- 1 + d + J-1

	Beta.hist <- metrop.out$par[,idx.Beta]
	Pi.hist <- t(apply(as.matrix(metrop.out$par[,idx.Pi]), 1, inv.mlogit))
	kappa.hist <- exp(metrop.out$par[,idx.kappa])

	ret <- list(par.hist = metrop.out$par, Beta.hist = Beta.hist,
		Pi.hist = Pi.hist, kappa.hist = kappa.hist, accept = metrop.out$accept,
		elapsed.sec = elapsed.sec, laplace.out = laplace.out,
		R.keep = nrow(Beta.hist), X.names = colnames(X))
	class(ret) <- "rwmetrop.mixlink.binomial"
	return(ret)
}

#' @name S3 methods for rwmetrop output objects
summary.rwmetrop.mixlink.binomial <- function(object, ...)
{
	summary.Beta <- data.frame(
		mean = colMeans(object$Beta.hist),
		sd = apply(object$Beta.hist, 2, sd),
		"pct2.5" = apply(object$Beta.hist, 2, quantile, prob = 0.025),
		"pct50" = apply(object$Beta.hist, 2, quantile, prob = 0.5),
		"pct97.5" = apply(object$Beta.hist, 2, quantile, prob = 0.975)
	)
	rownames(summary.Beta) <- object$X.names

	summary.Pi <- data.frame(
		mean = colMeans(object$Pi.hist),
		sd = apply(object$Pi.hist, 2, sd),
		"pct2.5" = apply(object$Pi.hist, 2, quantile, prob = 0.025),
		"pct50" = apply(object$Pi.hist, 2, quantile, prob = 0.5),
		"pct97.5" = apply(object$Pi.hist, 2, quantile, prob = 0.975)
	)
	rownames(summary.Pi) <- sprintf("Pi[%d]", 1:ncol(object$Pi.hist))

	summary.kappa <- data.frame(
		mean = mean(object$kappa.hist),
		sd = sd(object$kappa.hist),
		"pct2.5" = quantile(object$kappa.hist, prob = 0.025),
		"pct50" = quantile(object$kappa.hist, prob = 0.5),
		"pct97.5" = quantile(object$kappa.hist, prob = 0.975)
	)
	rownames(summary.kappa) <- "kappa"

	summary.theta <- rbind(summary.Beta, summary.Pi, summary.kappa)
	return(summary.theta)
}

#' @name S3 methods for rwmetrop output objects
print.rwmetrop.mixlink.binomial <- function(x, ...)
{
	s <- summary(x)
	print(s)
	printf("Elapsed sec: %0.2f", x$elapsed.sec)
	printf("  Draws kept: %d \n", x$R.keep)
	printf("Accept%% {%s}\n", paste(round(x$accept * 100, 2), collapse = ", "))
}

#' @name DIC for Random-Walk Metropolis Hastings Sampler
DIC.rwmetrop.mixlink.binomial <- function (metrop.out, y, m, X, invlink = plogis)
{
	R.keep <- metrop.out$R.keep
	D <- numeric(R.keep)
	for (r in 1:R.keep) {
		Beta <- metrop.out$Beta.hist[r, ]
		Pi <- metrop.out$Pi.hist[r, ]
		kappa <- metrop.out$kappa.hist[r]
		mean.hat <- invlink(X %*% Beta)
		D[r] <- -2 * sum(d.mixlink.binom(y, m, mean.hat, Pi, kappa, log = TRUE))
	}
	Beta.bar <- colMeans(metrop.out$Beta.hist)
	Pi.bar <- colMeans(metrop.out$Pi.hist)
	kappa.bar <- mean(metrop.out$kappa.hist)
	mean.bar <- invlink(X %*% Beta.bar)
	D.hat <- -2 * sum(d.mixlink.binom(y, m, mean.bar, Pi.bar, kappa.bar, log = TRUE))
	D.bar <- mean(D)
	p.D <- D.bar - D.hat
	dic <- D.bar + p.D
	return(dic)
}

