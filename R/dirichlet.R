ddirichlet <- function(x, alpha, log = FALSE)
{
  if (!is.matrix(x) && !is.matrix(alpha)) {
    x <- matrix(x, nrow = 1)
    alpha <- matrix(alpha,, nrow = 1)
  }

  if (is.matrix(x) && !is.matrix(alpha)) {
    alpha <- matrix(alpha, nrow(x), length(alpha), byrow = TRUE)
  }
  
  if (!is.matrix(x) && is.matrix(alpha)) {
    x <- matrix(x, nrow(alpha), length(x), byrow = TRUE)
  }
  
  stopifnot(nrow(x) == nrow(alpha) && ncol(x) == ncol(alpha))

  log.B <- lgamma(rowSums(alpha)) - rowSums(lgamma(alpha))
  log.f <- log.B + rowSums((alpha-1) * log(x))

  ## All of the values outside of Dirichlet sample space will now
  ## have log.f = NaN, so convert them to -Inf. It would be better to
  ## check the sample space purposefully though...
  idx <- which(is.nan(log.f))
  log.f[idx] <- -Inf

  if (log) return(log.f)
  else return(exp(log.f))
}

rdirichlet <- function(n, alpha)
{
  if (!is.matrix(alpha))
    alpha <- matrix(alpha, nrow = 1)
  
  k <- ncol(alpha)
  x <- matrix(rchisq(n*k, 2*alpha), nrow = n, byrow = TRUE)
  S <- rowSums(x) %*% t(rep(1,k))
  return(x / S)
}

