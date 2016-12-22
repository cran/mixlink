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
#' @name S3 methods for rwmetrop output objects
NULL

#' S3 methods for gibbs.mixlink.reg output objects
#'
#' @param object An \code{gibbs.mixlink.reg} object from \code{gibbs.mixlink.reg}.
#' @param x Same argument type as \code{object}.
#' @param ... Not currently used.
#' @name S3 methods for gibbs.mixlink.reg output objects
NULL

#' DIC for Random-Walk Metropolis Hastings Sampler
#'
#' Deviance Information Criteria (DIC) computed on results from
#' Metropolis-within-Gibbs sampler.
#'
#' @param metrop.out Results from Random-Walk Metropolis Hastings sampler.
#' @param y Observations.
#' @param X Design matrix for regression.
#' @param m Number of success/failure trials
#' @param offset Offset term to add to regression function, as commonly used in count
#'        models.
#' @param invlink The inverse link function for the mean. Default is \code{NULL},
#'        which indicates to use \code{plogis} for Binomial and \code{exp} for Poisson.
#'
#' @return DIC
#' @name DIC for Random-Walk Metropolis Hastings Sampler
NULL

