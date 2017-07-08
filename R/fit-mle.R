fit.mle <- function(phi.init, loglik, theta.tx, extra.tx, Data, psi.names = NULL)
{
	start.time <- Sys.time()
	stopifnot(!is.null(Data$n))
	n <- Data$n

	# Combine the two lists: theta.tx(phi) and extra.tx(phi)
	# into a vector. If the user provided psi.names, use those for the
	# variable names.
	psi.tx <- function(phi)
	{
		theta <- theta.tx(phi)
		psi1 <- unlist(theta)
		psi2 <- unlist(extra.tx(theta))
		psi <- c(psi1, psi2)
		
		if (!is.null(psi.names)) names(psi) <- psi.names
		return(psi)
	}

	optim.method <- getOption("optim.method")
	optim.control <- getOption("optim.control")
	if (is.null(optim.control$fnscale)) { optim.control$fnscale <- -1 }
	if (is.null(optim.control$maxit)) { optim.control$maxit <- 1e5 }

	optim.res <- optim(par = phi.init, fn = loglik, method = optim.method,
		control = optim.control, hessian = TRUE, Data = Data)

	phi.hat <- optim.res$par
	theta.hat <- theta.tx(phi.hat)
	xi.hat <- extra.tx(theta.hat)
	psi.hat <- psi.tx(phi.hat)

	V.phi <- -solve(optim.res$hessian)
	J.tx <- jacobian(psi.tx, phi.hat)
	V.psi <- J.tx %*% V.phi %*% t(J.tx)
	rownames(V.psi) <- colnames(V.psi) <- names(psi.hat)

	loglik.hat <- optim.res$value
	qq <- length(unlist(theta.hat))
	aic <- -2 * loglik.hat + 2*qq
	aicc <--2 * loglik.hat + 2*qq*n / (n-qq-1) 
	bic <- -2 * loglik.hat + qq*log(n)

	# Note: follow SAS NLMIXED, treat Est/SE as t with df = n
	df <- n
	se <- sqrt(diag(V.psi))
	t.val <- psi.hat / se
	p.val <- 2 * (1 - pt(abs(t.val), df = df))
	gr <- J.tx %*% grad(loglik, x = phi.hat, Data = Data)

	estimates <- cbind(psi.hat, se, t.val, p.val, gr)
	colnames(estimates) <- c("Estimate", "SE", "t-val", "P(|t|>t-val)", "Gradient")

	elapsed.sec <- as.numeric(Sys.time() - start.time, 'secs')

	res <- list(estimates = estimates, loglik = loglik.hat, aic = aic,
		aicc = aicc, bic = bic, vcov = V.psi, V.phi = V.phi, optim.res = optim.res, df = df,
		qq = qq, description = "<Default>", theta.hat = theta.hat,
		xi.hat = xi.hat, phi.init = phi.init, elapsed.sec = elapsed.sec)
	class(res) <- "mle.fit"
	return(res)
}

#' @name S3 methods for mle.fit objects
confint.mle.fit <- function(object, parm = NULL, level = 0.95, ...)
{
	if (!is.null(parm)) {
		stop("Only supported value of parm is currently NULL")
	}

	dim.theta <- length(unlist(object$theta.hat))
	dim.xi <- length(unlist(object$xi.hat))

	na <- rep(NA, dim.theta + dim.xi)
	DF <- data.frame(Estimate = na, SE = na, Lower = na, Upper = na)
	rownames(DF) <- rownames(object$estimates)

	w <- -qt((1-level)/2, df = object$df)
	psi.hat <- object$estimates[,1]
	se.psi.hat <- object$estimates[,2]
	DF$Lower <- psi.hat - w * se.psi.hat
	DF$Upper <- psi.hat + w * se.psi.hat
	DF$Estimate <- psi.hat
	DF$SE <- se.psi.hat

	res <- list(ci = DF, df = object$df, t.quantile = w, fit = object, level = level)
	class(res) <- "mle.fit.ci"
	return(res)
}

#' Print summary of an \code{mle.fit.ci} object
#'
#' See generic \code{print} function for usage.
#' @param x An \code{mle.fit.ci} object from \code{confint.fit.mle}.
#' @param ... Not currently used.
#' @name S3 methods for mle.fit.ci objects
print.mle.fit.ci <- function(x, ...)
{
	fit.out <- x$fit
	dim.theta <- length(unlist(fit.out$theta.hat))
	dim.xi <- length(unlist(fit.out$xi.hat))

	printf("--- Parameter CIs (level %f) ---\n", x$level)
	idx <- 1:dim.theta
	print(x$ci[idx,])

	if (dim.xi > 0)
	{
		printf("--- Additional CIs (level %f) ---\n", x$level)
		idx <- 1:dim.xi + dim.theta
		print(x$ci[idx,])
	}
	
	printf("--\n")
	printf("Degrees of freedom: %d\n", fit.out$df)
	printf("t-quantile: %f\n", x$t.quantile)
}

#' @name S3 methods for mle.fit objects
print.mle.fit <- function(x, ...)
{
	printf("Fit for model:\n")
	printf("%s\n", x$description)

	dim.theta <- length(unlist(x$theta.hat))
	dim.xi <- length(unlist(x$xi.hat))

	DF <- as.data.frame(x$estimates)
	DF[,1] <- round(DF[,1], 4)
	DF[,2] <- round(DF[,2], 4)
	DF[,3] <- round(DF[,3], 4)
	DF[,4] <- my.numerical.format(DF[,4])
	DF[,5] <- my.numerical.format(DF[,5])

	printf("--- Parameter Estimates ---\n")
	idx <- 1:dim.theta
	print(DF[idx,])

	if (dim.xi > 0)
	{
		printf("--- Additional Estimates ---\n")
		idx <- 1:dim.xi + dim.theta
		print(DF[idx,])
	}

	msg <- x$optim.res$message
	printf("--\n")
	printf("Elapsed Sec: %0.2f   ", x$elapsed.sec)
	printf("Degrees of freedom: %d\n", x$df)
	printf("LogLik: %0.4f   ", x$loglik)
	printf("AIC: %0.4f   ", x$aic)
	printf("AICC: %0.4f   ", x$aicc)
	printf("BIC: %0.4f\n", x$bic)
	printf("Converged status: %d   ", x$optim.res$convergence)
	printf("Message: %s\n", ifelse(is.null(msg), "<None>", msg))
}

#' @name S3 methods for mle.fit objects
summary.mle.fit <- function(object, ...)
{
	object
}

#' @name S3 methods for mle.fit objects
coef.mle.fit <- function(object, ...)
{
	object$estimates[,"Estimate"]
}

#' @name S3 methods for mle.fit objects
logLik.mle.fit <- function(object, ...)
{
	object$loglik
}

#' @name S3 methods for mle.fit objects
AIC.mle.fit <- function(object, ..., k=2)
{
	logLik(object) + k*object$qq
}

#' @name S3 methods for mle.fit objects
BIC.mle.fit <- function(object, ...)
{
	object$bic
}

#' @name S3 methods for mle.fit objects
vcov.mle.fit <- function(object, ...)
{
	object$vcov
}
