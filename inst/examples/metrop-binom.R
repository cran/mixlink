library(mixlink)
library(coda)
set.seed(1234)

# ----- Prepare the data -----
n <- 400
m <- rep(20, n)
x <- runif(n, -2, 2)
X <- model.matrix(~ x)
d <- ncol(X)
Beta.true <- c(-1, 1)
p.true <- pnorm(X %*% Beta.true)
kappa.true <- 2
Pi.true <- c(0.85, 0.15)
J <- length(Pi.true)

y <- r.mixlink.binom(n, p.true, Pi.true, kappa.true, m, save.latent = FALSE)

# ----- MLE for initial values -----
glm.out <- glm(cbind(y, m-y) ~ X-1, family = binomial(link = "probit"))
phi.init <- c(coef(glm.out), mlogit(normalize(1:J)), log(1))
mle.out <- mle.mixlink.binom.x(y, m, X, J=J, invlink.mean = pnorm,
  phi.init = phi.init)
print(mle.out)

# ----- MCMC with Metropolis Hastings -----
hyper <- list(var.Beta = 1000,
	alpha.Pi = rep(1,J),
	a.kappa = 1,
	b.kappa = 3/4)

proposal <- list(
	var = solve(-mle.out$optim.res$hessian),
	scale = 1.4
)

param.grp <- seq(1, d+J-1+1)
metrop.out <- rwmetrop.mixlink.binomial(y, X, m, R = 1000, burn = 0, thin = 1,
	invlink = pnorm, Beta.init = mle.out$theta.hat$Beta,
	Pi.init = mle.out$theta.hat$Pi, kappa.init = mle.out$theta.hat$kappa,
	hyper = hyper, use.laplace.approx = FALSE, report.period = 50,
	proposal = proposal, param.grp = param.grp)
print(metrop.out)

Beta.mcmc <- mcmc(metrop.out$Beta.hist)
Pi.mcmc <- mcmc(metrop.out$Pi.hist)
kappa.mcmc <- mcmc(metrop.out$kappa.hist)

DIC(metrop.out)
res <- residuals(metrop.out)
res.mean <- colMeans(res)
