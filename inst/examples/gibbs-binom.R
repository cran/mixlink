library(mixlink)
library(coda)
set.seed(1234)

# ----- Prepare the data -----
n <- 400
m <- rep(20, n)
x <- runif(n, -2, 2)
X <- model.matrix(~ x)
Beta.true <- c(-1, 1)
p.true <- pnorm(X %*% Beta.true)
kappa.true <- 1
Pi.true <- normalize(c(5, 25, 70))
d <- ncol(X)
J <- length(Pi.true)

y <- r.mixlink.binom(n, p.true, Pi.true, kappa.true, m, save.latent = FALSE)

# ----- MLE for initial values -----
J <- 2
glm.out <- glm(cbind(y,m-y) ~ X-1, family = binomial)
phi.init <- c(coef(glm.out), mlogit(normalize(1:J)), log(1))
mle.out <- mle.mixlink.binom.x(y, m, X, J=J, invlink.mean = pnorm,
  phi.init = phi.init)
V.hat <- solve(-mle.out$optim.res$hessian)
print(mle.out)

# ----- MCMC with Gibbs sampler -----
hyper <- list(
	VBeta = diag(1000, d),
	alpha = rep(1, J),
	kappa.a = 1,
	kappa.b = 1/2
)

proposal.VBeta <- 6 * V.hat[1:2, 1:2]
proposal.VPi <- 10 * V.hat[3, 3]
proposal.Vkappa <- 6 * V.hat[4, 4]
proposal.Vpsi <- diag(0.5^2, J)

gibbs.out <- gibbs.mixlink(y, X, R = 10000, burn = 1000, thin = 10,
	invlink = pnorm, report.period = 1000, save.psi = TRUE,
	Beta.init = mle.out$theta.hat$Beta,
  Pi.init = mle.out$theta.hat$Pi,
  kappa.init = mle.out$theta.hat$kappa,
	proposal.VBeta = proposal.VBeta, proposal.VPi = proposal.VPi,
	proposal.Vkappa = proposal.Vkappa, proposal.Vpsi = proposal.Vpsi,
	hyper = hyper, family = "binomial", trials = m, fixed.kappa = FALSE)
print(gibbs.out)

Beta.mcmc <- mcmc(gibbs.out$Beta.hist)
Pi.mcmc <- mcmc(gibbs.out$Pi.hist)
kappa.mcmc <- mcmc(gibbs.out$kappa.hist)

DIC(gibbs.out)
res <- residuals(gibbs.out)
res.mean <- colMeans(res)
