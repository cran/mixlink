library(mixlink)
library(coda)
set.seed(1238)

# ----- Prepare the data -----
n <- 400
m <- rep(20, n)
x <- rnorm(n, 0, 1)
X <- model.matrix(~ x)
Beta.true <- c(-1, 1)
p.true <- pnorm(X %*% Beta.true)
# kappa.true <- 2.5
kappa.true <- 1
Pi.true <- normalize(c(5, 25, 70))

d <- ncol(X)
J <- length(Pi.true)

S <- r.mixlink.binom(n, p.true, Pi.true, kappa.true, m, save.latent = TRUE)
y <- S$y
z.true <- S$z
psi.true <- S$psi

hyper <- list(
	VBeta = diag(1000, d),
	alpha = rep(1, J),
	kappa.a = 1,
	kappa.b = 1/2
)

proposal.VBeta <- solve(t(X) %*% X)
proposal.VPi <- diag(0.15^2, J-1)
proposal.Vkappa <- 0.5^2
proposal.Vpsi <- diag(0.5^2, J)

Beta.init <- c(0,0)
Pi.init <- normalize(1:J)
kappa.init <- kappa.true
psi.init <- psi.true

gibbs.pilot.out <- gibbs.mixlink.reg(y, X, R = 10000, burn = 1000, thin = 10,
	inv.link = pnorm, report.period = 1000, save.psi = TRUE,
	Beta.init = Beta.init, Pi.init = Pi.init, kappa.init = kappa.init, psi.init = psi.true,
	proposal.VBeta = proposal.VBeta, proposal.VPi = proposal.VPi,
	proposal.Vkappa = proposal.Vkappa, proposal.Vpsi = proposal.Vpsi,
	hyper = hyper, family = "binomial", trials = m, fixed.kappa = FALSE)
gibbs.out <- gibbs.pilot.out
# save.image("gibbs.pilot.out.Rdata")

temp <- matrix(NA, n*gibbs.pilot.out$R.keep, J)
for (j in 1:J) {
	temp[,j] <- as.numeric(gibbs.pilot.out$psi.hist[[j]])
}
new.proposal.Vpsi <- 5*var(temp)

gibbs.out <- gibbs.mixlink.reg(y, X, R = 20000, burn = 1000, thin = 20,
	inv.link = pnorm, report.period = 1000, save.psi = TRUE,
	Beta.init = as.numeric(tail(gibbs.pilot.out$Beta.hist,1)),
	Pi.init = as.numeric(tail(gibbs.pilot.out$Pi.hist,1)),
	kappa.init = tail(gibbs.pilot.out$kappa.hist,1),
	proposal.VBeta = 1 * var(gibbs.pilot.out$Beta.hist),
	proposal.VPi = 1 * var(t(apply(gibbs.pilot.out$Pi.hist, 1, mlogit))),
	proposal.Vkappa = 2 * var(log(gibbs.pilot.out$kappa.hist)),
	proposal.Vpsi = new.proposal.Vpsi,
	hyper = hyper, family = "binomial", trials = m)

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
plot(psi1.mcmc[,1:3])
plot(log.kappa.mcmc)
plot(qlogis.psi1.mcmc[,1:3])

acf(Beta.mcmc)
acf(Pi.mcmc)
acf(kappa.mcmc)
acf(log.kappa.mcmc)
acf(psi1.mcmc[,1:3])
acf(cbind(Beta.mcmc,Pi.mcmc,kappa.mcmc))

# ----- Compare to RW Metropolis-Hastings -----
# proposal <- list(var = diag(p+J-1+1), scale = 0.02)
 proposal <- list(var = diag(c(0.005, 0.005, 0.05, 0.05, 1.0)), scale = 1)
 hyper <- list(var.Beta = 1000, alpha.Pi = rep(1, J), a.kappa = 1, b.kappa = 1/10)
rwmetrop.out <- rwmetrop.mixlink.binomial(y, X, m, R = 15000, burn = 2000, thin = 10,
	invlink = pnorm, Beta.init = NULL, Pi.init = normalize(1:J),
	kappa.init = kappa.init, hyper = NULL, reporting.period = 100,
	use.laplace.approx = TRUE, proposal = proposal, param.grp = NULL)
print(rwmetrop.out)

rwmetrop.out.ordered <- rwmetrop.out
rwmetrop.out.ordered$Pi.hist <- t(apply(rwmetrop.out$Pi.hist, 1, sort))

Beta.metrop <- mcmc(rwmetrop.out.ordered$Beta.hist)
Pi.metrop <- mcmc(rwmetrop.out.ordered$Pi.hist)
kappa.metrop <- mcmc(rwmetrop.out.ordered$kappa.hist)

print(rwmetrop.out.ordered)

plot(Beta.metrop)
plot(Pi.metrop)
plot(kappa.metrop)

plot(mcmc(mlogit(rwmetrop.out.ordered$Pi.hist)))
plot(log(kappa.metrop))

summary(rwmetrop.out$Beta.hist)
summary(rwmetrop.out$Pi.hist)
summary(rwmetrop.out$kappa.hist)
plot(rwmetrop.out$Beta.hist[,1], type = "l")
plot(rwmetrop.out$Pi.hist[,1], type = "l")
plot(rwmetrop.out$kappa, type = "l")


# ----- Compare to RW Metropolis-Hastings WITHOUT GROUPING -----
proposal <- list(var = diag(c(0.005, 0.005, 0.05, 0.05)), scale = 0.5)
hyper <- list(var.Beta = 1000, alpha.Pi = rep(1, J), a.kappa = 1, b.kappa = 1/10)
rwmetrop2.out <- rwmetrop.mixlink.binomial(y, X, m, R = 15000, burn = 2000, thin = 10,
	invlink = pnorm, Beta.init = Beta.true, Pi.init = Pi.true,
	kappa.init = kappa.true, hyper = NULL, reporting.period = 100,
	use.laplace.approx = TRUE, proposal = proposal, param.grp = rep(1, d + J-1 + 1))
print(rwmetrop2.out)

rwmetrop2.out.ordered <- rwmetrop2.out
rwmetrop2.out.ordered$Pi.hist <- t(apply(rwmetrop2.out$Pi.hist, 1, sort))

Beta.metrop <- mcmc(rwmetrop2.out.ordered$Beta.hist)
Pi.metrop <- mcmc(rwmetrop2.out.ordered$Pi.hist)
kappa.metrop <- mcmc(rwmetrop2.out.ordered$kappa.hist)

print(rwmetrop2.out.ordered)

plot(Beta.metrop)
plot(Pi.metrop)
plot(kappa.metrop)

plot(mcmc(mlogit(Pi.metrop)))
plot(log(kappa.metrop))

save.image("results.Rdata")
