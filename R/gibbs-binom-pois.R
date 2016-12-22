#' Metropolis-within-Gibbs Sampler
#'
#' Metropolis-within-Gibbs sampler for Binomial and Poisson Mixture Link models.
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
#' @param trials Number of success/failure trials.
#' @param offset Offset term to add to regression function, as commonly used in count
#'        models.
#' @param R Number of MCMC draws to take.
#' @param burn Number of initial MCMC draws to discard.
#' @param thin After burn-in period, save one of every \code{thin} draws.
#' @param invlink The inverse link function for the mean. Default is \code{NULL},
#'        which indicates to use \code{plogis} for Binomial and \code{exp} for Poisson.
#' @param Beta.init Starting value for \eqn{\bm{\beta}}.
#' @param Pi.init Starting value for \eqn{\bm{\pi}}. The length of \code{Pi.init}
#'        determines \eqn{J} used in MCMC.
#' @param kappa.init Starting value for \eqn{\kappa}.
#' @param psi.init Starting value for \eqn{\bm{\psi}_1, \ldots, \bm{\psi}_n}
#' @param report.period Report progress every \code{report.period} draws.
#' @param family Can be either \code{binomial} or \code{poisson}.
#' @param fixed.kappa Keep \eqn{\kappa} fixed at \code{kappa.init} (Default FALSE).
#' @param proposal.VBeta Covariance matrix for \eqn{d} dimensional multivariate normal
#'        proposal. If left at the default (NULL), we use \code{diag(0.01^2, d)}.
#' @param proposal.VPi Covariance matrix for \eqn{J-1} dimensional multivariate normal
#'        proposal. If left at the default (NULL), we use \code{diag(0.01^2, J-1)}.
#' @param proposal.Vkappa Covariance matrix for univariate normal proposal. If left at
#'        the default (NULL), we use \code{0.01^2}.
#' @param proposal.Vpsi Covariance matrix for univariate normal proposal. If left at
#'        the default (NULL), we use \code{diag(0.01^2, J)}.
#' @param hyper A list with hyperparameters corresponding to the prior from the
#'        Details section. \code{var.Beta} is \eqn{V_{\beta}} with default 1000.
#'        \code{alpha.Pi} is \eqn{\bm{\alpha}_\pi} with default \code{rep(1,J)}.
#'        \code{a.kappa} and \code{b.kappa} correspond to \eqn{(a_\kappa, b_\kappa)}
#'        which have default values \code{a.kappa = 1} and \code{b.kappa = 1/10}.
#' @param save.psi Save draws for \eqn{\bm{\psi}_1, \ldots, \bm{\psi}_n}. Default is
#'        \code{FALSE}.
#'
#' @return A list with the MCMC results:
#' \item{Beta.hist}{\eqn{R \times d} matrix of saved \eqn{\bm{\beta}} draws.}
#' \item{Pi.hist}{\eqn{R \times J} matrix of saved \eqn{\bm{\pi}} draws.}
#' \item{kappa.hist}{\eqn{R \times 1} vector of \eqn{\kappa} draws.}
#' \item{psi.hist}{A list of \eqn{J} \eqn{R \times n} matrices. The \eqn{j}th matrix 
#'   contains the MCMC draws for \eqn{\bm{\psi}_{ij}}.}
#' \item{accept.Beta}{Percentages of MCMC proposals for \eqn{\bm{\beta}} were accepted.}
#' \item{accept.Pi}{Percentages of MCMC proposals for \eqn{\bm{\pi}} were accepted.}
#' \item{accept.kappa}{Percentages of MCMC proposals for \eqn{\kappa} were accepted.}
#' \item{accept.psi}{Percentages of MCMC proposals for \eqn{\bm{\psi}} were accepted.}
#' \item{elapsed.sec}{Elapsed time for sampling, in seconds.}
#' \item{R}{Number of total draws.}
#' \item{R.keep}{Number of draws kept, after thinning and burn-in.}
#' \item{X.names}{Names of columns of design matrix.}
#' \item{family}{Name of family / outcome.}
#' Can be accessed with the functions \code{print}, \code{summary}, and \code{DIC}.
#' @name Metropolis-within-Gibbs Sampler
#'
#' @references Andrew M. Raim, Nagaraj K. Neerchal, and Jorge G. Morel.
#'             An Extension of Generalized Linear Models to Finite
#'             Mixture Outcomes. arXiv preprint: 1612.03302
#' @examples
#' # ----- Generate data -----
#' n <- 200
#' m <- rep(20, n)
#' x <- rnorm(n, 0, 1)
#' X <- model.matrix(~ x)
#' Beta.true <- c(-1, 1)
#' mean.true <- plogis(X %*% Beta.true)
#' kappa.true <- 1
#' Pi.true <- c(1,3)/4
#' d <- ncol(X)
#' J <- length(Pi.true)
#' y <- r.mixlink.binom(n, mean.true, Pi.true, kappa.true, m)
#' 
#' # ----- Run Metropolis-within-Gibbs sampler -----
#' hyper <- list(VBeta = diag(1000, d), alpha.Pi = rep(1, J),
#' 	kappa.a = 1, kappa.b = 1/2)
#' gibbs.out <- gibbs.mixlink.reg(y, X, R = 10, burn = 5, thin = 1,
#' 	invlink = plogis, report.period = 100, Pi.init = c(1,9)/10,
#' 	proposal.VBeta = solve(t(X) %*% X), proposal.VPi = diag(0.25^2, J-1),
#' 	proposal.Vkappa = 0.5^2, proposal.Vpsi = diag(0.5^2, J),
#' 	hyper = hyper, family = "binomial", trials = m)
#' 
#' print(gibbs.out)
#' DIC.gibbs.mixlink.reg(gibbs.out, y, X, trials = m, invlink = plogis,
#' 	family = "binomial")
#' 
gibbs.mixlink.reg <- function(y, X, R, burn = 0, thin = 1,
	invlink = NULL, report.period = 100, save.psi = FALSE,
	Beta.init = NULL, Pi.init, kappa.init = NULL, psi.init = NULL,
	proposal.VBeta = NULL, proposal.VPi = NULL, proposal.Vkappa = NULL,
	proposal.Vpsi = NULL, hyper = NULL, trials = NULL, offset = rep(0, length(y)),
	family = NULL, fixed.kappa = FALSE)
{
	start.time <- Sys.time()

	n <- nrow(X)
	d <- ncol(X)
	J <- length(Pi.init)

	R.keep <- ceiling((R - burn) / thin)
	r.keep <- 0

	if (is.null(family)) {
		stop("Must specify a family (currently only binomial or poisson)")
	} else if (family == "binomial") {
		if (is.null(invlink)) { invlink <- plogis }
		Q <- Q.binomial
		if (is.null(trials)) { stop("Must specify trials for binomial") }
		m <- trials
	} else if (family == "poisson") {
		if (is.null(invlink)) { invlink <- exp }
		Q <- Q.poisson
		if (!is.null(trials)) { warning("For Poisson, trials argument will be ignored") }
		m <- rep(0,n)
	} else {
		stop("Must specify a family (currently only binomial or poisson)")
	}

	Beta.hist <- matrix(NA, R.keep, d)
	Pi.hist <- matrix(NA, R.keep, J)
	kappa.hist <- numeric(R.keep)
	psi.hist <- list()
	if (save.psi) {
		for (j in 1:J) {
			psi.hist[[j]] <- matrix(NA, R.keep, n)
		}
	}

	if (is.null(proposal.VBeta)) {
		logger("Using default proposal.VBeta: diag(0.01^2, d)\n")
		proposal.VBeta <- diag(0.01^2, d)
	}
	if (is.null(proposal.VPi)) {
		logger("Using default proposal.VPi: diag(0.01^2, J-1)\n")
		proposal.VPi <- diag(0.01^2, J-1)
	}
	if (is.null(proposal.Vkappa)) {
		logger("Using default proposal.Vkappa: 0.01^2\n")
		proposal.Vkappa <- 0.01^2
	}
	if (is.null(proposal.Vpsi)) {
		logger("Using default proposal.Vpsi: diag(0.01^2, J)\n")
		proposal.Vpsi <- diag(0.01^2, J)
	}

	proposal.VBeta.half <- chol(proposal.VBeta)
	proposal.VPi.half <- chol(proposal.VPi)
	proposal.Vkappa.half <- sqrt(proposal.Vkappa)
	proposal.Vpsi.half <- chol(proposal.Vpsi)

	# Initial values
	if (is.null(Beta.init)) {
		Beta.init <- rep(0, d)
	}
	if (is.null(kappa.init)) {
		kappa.init <- 1
	}
	if (is.null(psi.init)) {
		psi.init <- plogis(matrix(rnorm(n*J), n, J))
	}

	if (is.null(hyper)) { hyper <- list() }
	if (is.null(hyper$VBeta)) { hyper$VBeta <- diag(1000, d) }
	if (is.null(hyper$alpha.Pi)) { hyper$alpha.Pi <- rep(1, J) }
	if (is.null(hyper$kappa.a)) { hyper$kappa.a <- 1 }
	if (is.null(hyper$kappa.b)) { hyper$kappa.b <- 1/10 }

	Beta <- Beta.init
	Pi <- Pi.init
	kappa <- kappa.init
	psi <- psi.init

	accept.Beta <- 0
	accept.Pi <- 0
	accept.kappa <- 0
	accept.psi <- numeric(n)

	logger("Starting MCMC\n")

	for (r in 1:R)
	{
		if (r %% report.period == 0) {
			logger("Starting rep %d, ", r)
			printf("accept%% {beta %0.2f, Pi %0.2f, kappa %0.2f, psi [%0.2f %0.2f]}\n",
				   accept.Beta / (r-1) * 100, accept.Pi / (r-1) * 100,
				   accept.kappa / (r-1) * 100, min(accept.psi) / (r-1) * 100,
				   max(accept.psi) / (r-1) * 100)
		}

		# Draw Beta | Rest
		Beta_ <- t(proposal.VBeta.half) %*% rnorm(d) + Beta
		mean.hat_ <- invlink(X %*% Beta_ + offset)
		mean.hat <- invlink(X %*% Beta + offset)
		log.num <- sum(Q(y, m, psi, mean.hat_, Pi, kappa)) +
			dmvnorm(t(Beta_), rep(0,d), hyper$VBeta, log = TRUE)
		log.den <- sum(Q(y, m, psi, mean.hat, Pi, kappa)) +
			dmvnorm(t(Beta), rep(0,d), hyper$VBeta, log = TRUE)
		logAlpha <- log.num - log.den
		if (!is.na(logAlpha)) {
			if (log(runif(1)) < logAlpha) {
				Beta <- Beta_
				accept.Beta <- accept.Beta + 1
			}
		}
		mean.hat <- invlink(X %*% Beta + offset)

		# Draw Pi | Rest
		trans.Pi_ <- t(proposal.VPi.half) %*% rnorm(J-1) + mlogit(Pi)
		Pi_ <- inv.mlogit(trans.Pi_)
		jac_ <- jacobian(func = inv.mlogit, x = mlogit(Pi_))
		jac <- jacobian(func = inv.mlogit, x = mlogit(Pi))
		l_.tx <- log(abs(det(as.matrix(jac_[-J,]))))
		l.tx <- log(abs(det(as.matrix(jac[-J,]))))
		log.num <- sum(Q(y, m, psi, mean.hat, Pi_, kappa)) +
			ddirichlet(Pi_, hyper$alpha.Pi, log = TRUE) + l_.tx
		log.den <- sum(Q(y, m, psi, mean.hat, Pi, kappa)) +
			ddirichlet(Pi, hyper$alpha.Pi, log = TRUE) + l.tx
		logAlpha <- log.num - log.den
		if (!is.na(logAlpha)) {
			if (log(runif(1)) < logAlpha) {
				Pi <- Pi_
				accept.Pi <- accept.Pi + 1
			}
		}

		# Draw kappa | Rest
		if (!fixed.kappa) {
			trans.kappa_ <- proposal.Vkappa.half * rnorm(1) + log(kappa)
			kappa_ <- exp(trans.kappa_)
			l_.tx <- log(kappa_)
			l.tx <- log(kappa)

			log.num <- sum(Q(y, m, psi, mean.hat, Pi, kappa_)) +
				dgamma(kappa_, hyper$kappa.a, hyper$kappa.b, log = TRUE) +
				l_.tx
			log.den <- sum(Q(y, m, psi, mean.hat, Pi, kappa)) +
				dgamma(kappa, hyper$kappa.a, hyper$kappa.b, log = TRUE) + 
				l.tx
			logAlpha <- log.num - log.den
			if (!is.na(logAlpha)) {
				if (log(runif(1)) < logAlpha) {
					kappa <- kappa_
					accept.kappa <- accept.kappa + 1
				}
			}
		}

		# Draw psi | Rest
		# We can draw them independently
		trans.psi <- qlogis(psi)
		trans.psi_ <- t(t(proposal.Vpsi.half) %*% matrix(rnorm(n*J), J, n) + t(trans.psi))
		psi_ <- plogis(trans.psi_)
		qq_ <- Q(y, m, psi_, mean.hat, Pi, kappa)
		qq <- Q(y, m, psi, mean.hat, Pi, kappa)
		l_.tx <- rowSums(dlogis(trans.psi_, log = TRUE))
		l.tx <- rowSums(dlogis(trans.psi, log = TRUE))
		log.num <- qq_ + l_.tx
		log.den <- qq + l.tx
		logAlpha <- log.num - log.den
		idx.na <- which(is.na(logAlpha))
		logAlpha[idx.na] <- -Inf
		idx <- which(log(runif(n)) < logAlpha)
		psi[idx,] <- psi_[idx,]
		accept.psi[idx] <- accept.psi[idx] + 1

		# Save draws to history
		if (r > burn && r %% thin == 0) {
			r.keep <- r.keep + 1
			Beta.hist[r.keep,] <- Beta
			Pi.hist[r.keep,] <- Pi
			kappa.hist[r.keep] <- kappa

			if (save.psi) {
				for (j in 1:J) {
					psi.hist[[j]][r.keep,] <- psi[,j]
				}
			}
		}
	}

	logger("Finished MCMC\n")
	elapsed.sec <- as.numeric(Sys.time() - start.time, unit = "secs")

	ret <- list(Beta.hist = Beta.hist, Pi.hist = Pi.hist,
		kappa.hist = kappa.hist, psi.hist = psi.hist, X.names = colnames(X),
		elapsed.sec = elapsed.sec, accept.Beta = accept.Beta,
		accept.Pi = accept.Pi, accept.kappa = accept.kappa,
		accept.psi = accept.psi, R = R, R.keep = R.keep,
		family = family)
	class(ret) <- "gibbs.mixlink.reg"
	return(ret)
}

#' @name S3 methods for gibbs.mixlink.reg output objects
summary.gibbs.mixlink.reg <- function(object, ...)
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
	J <- ncol(object$Pi.hist)
	rownames(summary.Pi) <- sprintf("Pi[%d]", 1:J)

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

#' @name S3 methods for gibbs.mixlink.reg output objects
print.gibbs.mixlink.reg <- function(x, ...)
{
	s <- summary.gibbs.mixlink.reg(x)
	print(s)
	printf("Family: %s\n", x$family)
	printf("Elapsed sec: %0.2f\n", x$elapsed.sec)
	printf("Accept%% {beta %0.2f, Pi %0.2f, kappa %0.2f, psi [%0.2f %0.2f]}\n",
		x$accept.Beta / x$R * 100, x$accept.Pi / x$R * 100,
		x$accept.kappa / x$R * 100, min(x$accept.psi) / x$R * 100,
		max(x$accept.psi) / x$R * 100)
}

#' DIC for Metropolis-within-Gibbs
#'
#' Deviance Information Criteria (DIC) computed on results from
#' Metropolis-within-Gibbs sampler.
#'
#' @param gibbs.out Results from Metropolis-within-Gibbs sampler.
#' @param y Observations.
#' @param X Design matrix for regression.
#' @param trials Number of success/failure trials.
#' @param offset Offset term to add to regression function, as commonly used in count
#'        models.
#' @param invlink The inverse link function for the mean. Default is \code{NULL},
#'        which indicates to use \code{plogis} for Binomial and \code{exp} for Poisson.
#' @param family Can be either \code{binomial} or \code{poisson}.
#'
#' @return DIC
#' @name DIC for Metropolis-within-Gibbs Sampler
DIC.gibbs.mixlink.reg <- function(gibbs.out, y, X, trials = NULL, offset = NULL,
	invlink = NULL, family = NULL)
{
	m <- trials
	R.keep <- gibbs.out$R.keep
	D <- numeric(R.keep)

	for (r in 1:R.keep) {
		Beta <- gibbs.out$Beta.hist[r,]
		Pi <- gibbs.out$Pi.hist[r,]
		kappa <- gibbs.out$kappa.hist[r]
		mean.hat <-  invlink(X %*% Beta)

		if (family == "binomial") {
			D[r] <- -2 * sum(d.mixlink.binom(y, m, mean.hat, Pi, kappa, log = TRUE))
		} else if (family == "poisson") {
			D[r] <- -2 * sum(d.mixlink.pois(y, mean.hat, Pi, kappa, log = TRUE))
		} else {
			stop("Only binomial and poisson families are supported")
		}
	}

	Beta.bar <- colMeans(gibbs.out$Beta.hist)
	Pi.bar <- colMeans(gibbs.out$Pi.hist)
	kappa.bar <- mean(gibbs.out$kappa.hist)
	mean.bar <- invlink(X %*% Beta.bar)
	
	if (family == "binomial") {
		D.hat <- -2 * sum(d.mixlink.binom(y, m, mean.bar, Pi.bar, kappa.bar, log = TRUE))
	} else if (family == "poisson") {
		D.hat <- -2 * sum(d.mixlink.pois(y, mean.bar, Pi.bar, kappa.bar, log = TRUE))
	} else {
		stop("Only binomial and poisson families are supported")
	}

	D.bar <- mean(D)
	p.D <- D.bar - D.hat
	dic <- D.bar + p.D

	return(dic)
}

Q.binomial <- function(y, m, psi, mean, Pi, kappa, find_vert_tol = 1e-08)
{
	n <- length(y)
	if (length(mean) == 1) { mean <- rep(mean, n) }
	if (length(kappa) == 1) { kappa <- rep(kappa, n) }
	res <- .Call("mixlink_gibbs_Q_binom", as.integer(y), as.integer(m), mean,
		psi, Pi, kappa, find_vert_tol)
	return(res)
}

Q.poisson <- function(y, m, psi, mean, Pi, kappa, find_vert_tol = 1e-08)
{
	n <- length(y)
	if (length(mean) == 1) { mean <- rep(mean, n) }
	if (length(kappa) == 1) { kappa <- rep(kappa, n) }
	res <- .Call("mixlink_gibbs_Q_pois", as.integer(y), mean, psi, Pi, kappa)
	return(res)
}

