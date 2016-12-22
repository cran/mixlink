library(mixlink)
library(coda)

set.seed(1234)

# ----- Prepare the data -----
n <- 400
x <- runif(n, 0, 8)
X <- model.matrix(~ x)
Beta.true <- c(-1, 0.5)
mean.true <- exp(X %*% Beta.true)
kappa.true <- 1.25
Pi.true <- normalize(c(1,5))

d <- ncol(X)
J <- length(Pi.true)

res <- r.mixlink.pois(n, mean.true, Pi.true, kappa.true, save.latent = TRUE)
y <- res$y
z <- res$z

mle.out <- mle.mixlink.pois.x(y, X, J)
print(mle.out)

mean.hat <- exp(X %*% mle.out$theta.hat$Beta)
Pi.hat <- mle.out$theta.hat$Pi
kappa.hat <- mle.out$theta.hat$kappa
res <- rqres.mixlink.pois(y, mean.hat, Pi.hat, kappa.hat)
qqnorm(res); qqline(res)
plot(mean.hat, res)
