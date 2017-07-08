#' Deviance Information Criteria for MCMC Draws
#'
#' Deviance Information Criteria (DIC) computed on results from
#' MCMC sampler.
#'
#' @param object Results from sampler.
#' @param ... optional arguments (currently not used).
#'
#' @return DIC
#' @name DIC
DIC <- function (object, ...)
{
	UseMethod("DIC")
}
