library(mixlink)
library(coda)

set.seed(1234)

# ----- Prepare the data -----
n <- 400
x <- runif(n, 0, 8)
X <- model.matrix(~ x)
Beta.true <- c(-1, 0.5)
mean.true <- exp(X %*% Beta.true)
kappa.true <- 2.5
Pi.true <- normalize(c(1,2))
d <- ncol(X)
J <- length(Pi.true)

y <- r.mixlink.pois(n, mean.true, Pi.true, kappa.true, save.latent = FALSE)

# ----- MLE for initial values -----
glm.out <- glm(y ~ X-1, family = poisson)
phi.init <- c(coef(glm.out), mlogit(normalize(1:J)), log(1))
mle.out <- mle.mixlink.pois.x(y, X, J=J, invlink.mean = exp,
  phi.init = phi.init)
V.hat <- solve(-mle.out$optim.res$hessian)
print(mle.out)

# ----- MCMC with Gibbs sampler -----
hyper <- list(
	VBeta = diag(1000, d),
	alpha = rep(1, J),
	kappa.a = 1,
	kappa.b = 0.5
)

proposal.VBeta <- 4 * V.hat[1:2,1:2]
proposal.VPi <- 4 * V.hat[3,3]
proposal.Vkappa <- 4 * V.hat[4,4]
proposal.Vpsi <- diag(0.5^2, J)

gibbs.out <- gibbs.mixlink(y, X, R = 10000, burn = 3000, thin = 10,
	invlink = exp, report.period = 1000, save.psi = TRUE,
	Beta.init = mle.out$theta.hat$Beta,
  Pi.init = mle.out$theta.hat$Pi,
  kappa.init = mle.out$theta.hat$kappa,
	proposal.VBeta = proposal.VBeta, proposal.VPi = proposal.VPi,
	proposal.Vkappa = proposal.Vkappa, proposal.Vpsi = proposal.Vpsi,
	hyper = hyper, offset = rep(0,n), family = "poisson")
print(gibbs.out)

Beta.mcmc <- mcmc(gibbs.out$Beta.hist)
Pi.mcmc <- mcmc(gibbs.out$Pi.hist)
kappa.mcmc <- mcmc(gibbs.out$kappa.hist)

DIC(gibbs.out)
res <- residuals(gibbs.out)
res.mean <- colMeans(res)
