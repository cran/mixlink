library(mixlink)
library(coda)

set.seed(1234)

# ----- Prepare the data -----
n <- 400
x <- runif(n, -3, 3)
X <- model.matrix(~ x)
m <- rpois(n, 20)
Beta.true <- c(-1, 2)
mean.true <- plogis(X %*% Beta.true)
kappa.true <- 1.25
Pi.true <- normalize(c(1,5))

d <- ncol(X)
J <- length(Pi.true)

res <- r.mixlink.binom(n, mean.true, Pi.true, kappa.true, m, save.latent = TRUE)
y <- res$y
z <- res$z

options(OMIR.optim.method = "BFGS", OMIR.optim.control = list())
mle.out <- mle.mixlink.binom.x(y, m, X, J)
print(mle.out)

mean.hat <- plogis(X %*% mle.out$theta.hat$Beta)
Pi.hat <- mle.out$theta.hat$Pi
kappa.hat <- mle.out$theta.hat$kappa

res <- rqres.mixlink.binom(y, m, mean.hat, Pi.hat, kappa.hat)
qqnorm(res); qqline(res)
plot(mean.hat, res)
