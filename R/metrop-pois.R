#' @name Metropolis-Hastings Sampler
rwmetrop.mixlink.poisson <- function(y, X, offset = rep(0, length(y)),
  R, burn = 0, thin = 1, invlink = exp, Beta.init = NULL, Pi.init,
  kappa.init = NULL, hyper = NULL, report.period = R+1,
  use.laplace.approx = TRUE, proposal = NULL, param.grp = NULL,
  fixed.kappa = FALSE)
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

	Data <- list(y = y, X = X, J = J)

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
		ll <- d.mixlink.pois(Data$y, p, Pi, kappa, log = TRUE)
		if (any(is.na(ll)) || any(is.infinite(ll))) {
		  browser()
		}

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
		R = R, burn = burn, thin = thin, report.period = report.period)
	logger("Finished MCMC\n")

	elapsed.sec <- as.numeric(Sys.time() - st, unit = "secs")

	idx.Beta <- seq(1, d)
	idx.Pi <- seq(1, J-1) + d
	idx.kappa <- 1 + d + J-1

	Beta.hist <- metrop.out$par[,idx.Beta]
	Pi.hist <- t(apply(as.matrix(metrop.out$par[,idx.Pi]), 1, inv.mlogit))
	kappa.hist <- exp(metrop.out$par[,idx.kappa])

	ret <- list(par.hist = metrop.out$par, Beta.hist = Beta.hist,
		Pi.hist = Pi.hist, kappa.hist = kappa.hist,
		accept = metrop.out$accept, elapsed.sec = elapsed.sec,
		laplace.out = laplace.out, R.keep = nrow(Beta.hist),
		X.names = colnames(X), y = y, X = X, offset = offset,
	  invlink = invlink, report.period = report.period)

	class(ret) <- "rwmetrop.mixlink.poisson"
	return(ret)
}

#' @name S3 methods for rwmetrop output objects
summary.rwmetrop.mixlink.poisson <- function(object, ...)
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
print.rwmetrop.mixlink.poisson <- function(x, ...)
{
	s <- summary(x)
	print(s)
	printf("Elapsed sec: %0.2f", x$elapsed.sec)
	printf("  Draws kept: %d \n", x$R.keep)
	printf("Accept%% {%s}\n", paste(round(x$accept * 100, 2), collapse = ", "))
}

#' @name S3 methods for rwmetrop output objects
DIC.rwmetrop.mixlink.poisson <- function (object, ...)
{
  y <- object$y
  X <- object$X
  offset <- object$offset

	R.keep <- object$R.keep
	D <- numeric(R.keep)
	for (r in 1:R.keep) {
    if (r %% object$report.period == 0) {
      logger("Computing DIC for rep %d\n", r)
    }
		Beta <- object$Beta.hist[r, ]
		Pi <- object$Pi.hist[r, ]
		kappa <- object$kappa.hist[r]
		mean.hat <- object$invlink(X %*% Beta + offset)
		D[r] <- -2 * sum(d.mixlink.pois(y, mean.hat, Pi, kappa, log = TRUE))
	}
	Beta.bar <- colMeans(object$Beta.hist)
	Pi.bar <- colMeans(object$Pi.hist)
	kappa.bar <- mean(object$kappa.hist)
	mean.bar <- object$invlink(X %*% Beta.bar + offset)
	D.hat <- -2 * sum(d.mixlink.pois(y, mean.bar, Pi.bar, kappa.bar, log = TRUE))
	D.bar <- mean(D)
	p.D <- D.bar - D.hat
	dic <- D.bar + p.D
	return(dic)
}

#' @name S3 methods for rwmetrop output objects
residuals.rwmetrop.mixlink.poisson <- function(object, ...)
{
  y <- object$y
  X <- object$X
  offset <- object$offset
  n <- nrow(X)
  R.keep <- object$R.keep

  rqres.pp <- matrix(NA, R.keep, n)
  for (r in 1:R.keep) {
    if (r %% object$report.period == 0) {
      logger("Computing residuals for rep %d\n", r)
    }
  	Beta <- object$Beta.hist[r,]
  	Pi <- object$Pi.hist[r,]
  	kappa <- object$kappa.hist[r]
  	mean.hat <- object$invlink(X %*% Beta + offset)
  	rqres.pp[r,] <- rqres.mixlink.pois(y, mean.hat, Pi, kappa)
  }

  return(rqres.pp)
}
