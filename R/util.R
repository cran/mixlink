old.options <- options(
	optim.method = "L-BFGS-B",
	optim.control = list()
)
on.exit(options(old.options), add = TRUE)

printf <- function(msg, ...) {
	cat(sprintf(msg, ...))
}

logger <- function(msg, ...) {
	sys.time <- as.character(Sys.time())
	cat(sys.time, "-", sprintf(msg, ...))
}

#' Normalize
#'
#' Scale a numeric vector by the sum of its elements
#'
#' @param x A numeric vector
#'
#' @return \code{x / sum(x)}
#'
#' @examples
#'   x <- c(1,1,1,1)
#'   normalize(x)
normalize <- function(x) { x / sum(x) }

null.tx <- function(phi) { list() }

my.numerical.format <- function(x, lower = 1e-4) {
	idx1 <- which(abs(x) < lower)
	idx2 <- setdiff(1:length(x), idx1)
	y <- character(length(x))
	y[idx1] <- sprintf("%0.3E", x[idx1])
	y[idx2] <- sprintf("%0.4f", x[idx2])
	names(y) <- names(x)
	return(y)
}

mlogit <- function(p) {
	J <- length(p)
	x <- log(p[-J] / p[J])
	return(x)
}

# Transform from R^(J-1) to probability simplex S^J
inv.mlogit <- function(x) {
	z <- exp(x)
	P.J <- 1 / (1 + sum(z))
	p <- c(z * P.J, P.J)
	return(p)
}

# Confluent geometric function of the first kind
hypergeomF1 <- function(x, a, b, log = FALSE) {
	n <- max(length(a), length(b), length(x))
	x <- rep_len(x, n)
	a <- rep_len(a, n)
	b <- rep_len(b, n)
	
	ff <- .Call("hyperg_1F1", a, b, x)
	if (log) log(ff)
	else ff
}

mbeta <- function(t, a, b, log = FALSE) {
	hypergeomF1(t, a, a+b, log)
}

