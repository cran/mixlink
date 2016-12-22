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

res <- r.mixlink.pois(n, mean.true, Pi.true, kappa.true, save.latent = TRUE)
y <- res$y
z <- res$z

plot(x, y)
points(x, exp(Beta.true[1] + Beta.true[2]*x), col = "blue", cex = 0.5)

hyper <- list(
	VBeta = diag(1000, d),
	alpha = rep(1, J),
	kappa.a = 1,
	kappa.b = 0.5
)

# proposal.VBeta <- diag(0.5^2, d)
proposal.VBeta <- solve(t(X) %*% X)
proposal.VPi <- diag(0.2^2, J-1)
proposal.Vkappa <- 0.5^2
proposal.Vpsi <- diag(0.5^2, J)

Beta.init <- c(0,0)
Pi.init <- normalize(1:J)
kappa.init <- 1

gibbs.pilot.out <- gibbs.mixlink.reg(y, X, R = 5000, burn = 1000, thin = 10,
	inv.link = NULL, report.period = 1000, save.psi = TRUE,
	Beta.init = Beta.init, Pi.init = Pi.init, kappa.init = kappa.init,
	proposal.VBeta = proposal.VBeta, proposal.VPi = proposal.VPi,
	proposal.Vkappa = proposal.Vkappa, proposal.Vpsi = proposal.Vpsi,
	hyper = hyper, offset = rep(0,n), family = "poisson")
gibbs.out <- gibbs.pilot.out

temp <- matrix(NA, n*gibbs.pilot.out$R.keep, J)
for (j in 1:J) {
	temp[,j] <- as.numeric(gibbs.pilot.out$psi.hist[[j]])
}
new.proposal.Vpsi <- 10*var(temp)

gibbs.out <- gibbs.mixlink.reg(y, X, R = 10000, burn = 1000, thin = 10,
	inv.link = NULL, report.period = 1000, save.psi = TRUE,
	Beta.init = as.numeric(tail(gibbs.pilot.out$Beta.hist, 1)),
	Pi.init = as.numeric(tail(gibbs.pilot.out$Pi.hist, 1)),
	kappa.init = as.numeric(tail(gibbs.pilot.out$kappa.hist, 1)),
	proposal.VBeta = 3 * var(gibbs.pilot.out$Beta.hist),
	proposal.VPi = 3 * var(apply(gibbs.pilot.out$Pi.hist, 1, mlogit)),
	proposal.Vkappa = 2 * var(log(gibbs.pilot.out$kappa.hist)),
	proposal.Vpsi = proposal.Vpsi,
	hyper = hyper, offset = rep(0,n), family = "poisson")

# Experimental: Let's apply a post-MCMC ident. constraint on Pis
gibbs.out$Pi.hist <- t(apply(gibbs.out$Pi.hist, 1, sort))

print(gibbs.out)

R.keep <- gibbs.out$R.keep
idx <- seq(0, R.keep, by = 1)

Beta.mcmc <- mcmc(gibbs.out$Beta.hist[idx,])
Pi.mcmc <- mcmc(gibbs.out$Pi.hist[idx,])
kappa.mcmc <- mcmc(gibbs.out$kappa.hist[idx])
log.kappa.mcmc <- mcmc(log(gibbs.out$kappa.hist[idx]))
psi1.mcmc <- mcmc(gibbs.out$psi[[1]][idx,])
qlogis.psi1.mcmc <- mcmc(qlogis(gibbs.out$psi[[1]][idx,]))

plot(Beta.mcmc)
plot(Pi.mcmc)
plot(kappa.mcmc)
plot(log.kappa.mcmc)
plot(psi1.mcmc[,1:3])
plot(qlogis.psi1.mcmc[,1:3])

acf(Beta.mcmc)
acf(Pi.mcmc)
acf(kappa.mcmc)
acf(log.kappa.mcmc)
acf(psi1.mcmc[,1:3])

cor(gibbs.out$Pi.hist[idx,J], gibbs.out$kappa[idx])
