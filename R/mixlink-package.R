#' S3 methods for mle.fit objects
#'
#' @param object An \code{mle.fit} object from \code{fit.mle}.
#' @param x An \code{mle.fit} object from \code{fit.mle}.
#' @param parm Not currently used.
#' @param level Desired confidence level.
#' @param k Penalty term is \eqn{k q} where \eqn{q =\dim \bm{\theta}}. \eqn{k=2} is the
#'        default for AIC calculation.
#' @param ... Not currently used.
#' @name S3 methods for mle.fit objects
NULL

#' S3 methods for rwmetrop output objects
#'
#' @param object An \code{rwmetrop.mixlink.binomial} object from \code{rwmetrop.mixlink.binomial}
#'        or an \code{rwmetrop.mixlink.poisson} object from \code{rwmetrop.mixlink.poisson}.
#' @param x Same argument type as \code{object}.
#' @param ... Not currently used.
#' @details \code{DIC} returns the Deviance Information Criterion computed from MCMC draws.
#'   \code{residuals} returns randomized quantile residuals computed on each MCMC draw.
#' @name S3 methods for rwmetrop output objects
NULL

#' S3 methods for gibbs.mixlink output objects
#'
#' @param object An \code{gibbs.mixlink} object from \code{gibbs.mixlink}.
#' @param x Same argument type as \code{object}.
#' @param ... Not currently used.
#' @details \code{DIC} returns the Deviance Information Criterion computed from MCMC draws.
#'   \code{residuals} returns randomized quantile residuals computed on each MCMC draw.
#' @name S3 methods for gibbs.mixlink output objects
NULL

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
#' @param offset Constant offset term to add to \eqn{x^T \beta}.
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
#' @param report.period Report progress every \code{report.period} draws.
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
#' 	report.period = 1000, use.laplace.approx = TRUE, proposal = proposal)
#' 
#' print(metrop.out)
#' DIC(metrop.out)
#' }
#'
#' @references Andrew M. Raim, Nagaraj K. Neerchal, and Jorge G. Morel.
#'             An Extension of Generalized Linear Models to Finite
#'             Mixture Outcomes. arXiv preprint: 1612.03302
NULL
