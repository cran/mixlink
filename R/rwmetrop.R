# This sampler and laplace code was adapted from rwmetrop function in
# LearnBayes package
rwmetrop <- function(par.init, logpost, Data, proposal, R,
	grp = rep(1, length(par.init)), burn = 0, thin = 1, report.period = 1000)
{
	qq <- length(par.init)
	V.proposal.half.trans <- t(proposal$scale * chol(proposal$var))
	logf <- function(par) { logpost(par, Data) }
	if (is.infinite(logf(par.init)) | is.na(logf(par.init))) {
	  stop("logpost(par.init) must be a finite value")
	}
	grp.idx.list <- split(1:qq - 1, grp)
	ret <- rwmetrop_cpp(par.init, logf, V.proposal.half.trans,
		grp.idx.list, R, burn, thin, report.period)
	return(ret)
}

laplace <- function (logpost, mode, Data, optim.control = list(),
	optim.method = "L-BFGS-B")
{
	optim.control$fnscale <- -1
	fit <- optim(mode, logpost, gr = NULL, Data, hessian = TRUE,
		method = optim.method, control = optim.control)

	mode <- fit$par
	H <- -solve(fit$hessian)
	p <- length(mode)
	int <- p/2 * log(2 * pi) + 0.5 * log(det(H)) + logpost(mode, Data)

	list(mode = mode, var = H, int = int, converge = (fit$convergence == 0), optim.out = fit)
}

