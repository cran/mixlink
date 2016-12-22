library(OverdispersionModelsInR)
library(mixlink)

# It looks like we have a label switching problem when J > 2
# 
# (SAME RESULT) What if we try grouped sampling?
# What if we try another sampler? Like NUTS or HMC?
# 	Look for literature on this
# Could there be a fundamental problem with the model when J > 2?
# 	Remember that model is invariant to permutations in Pi
# For the purpose of estimating betas, it label switching might not make
# much difference. But notice that kappa is also screwed up.

set.seed(1234)

# ----- Prepare the data -----
n <- 400
m <- rep(20, n)
x <- rnorm(n, 0, 1)
X <- model.matrix(~ x)
d <- ncol(X)
Beta.true <- c(-1, 1)
p.true <- pnorm(X %*% Beta.true)
kappa.true <- 2
# Pi.true <- c(0.85, 0.10, 0.05)
Pi.true <- c(0.85, 0.15)

y <- r.mixlink.binom(n, p.true, Pi.true, kappa.true, m, save.latent = FALSE)

# ----- MLE, for comparison -----
# This part is from OverdispersionModelsInR
# J <- 3
J <- 2
glm.out <- glm(cbind(y, m-y) ~ X-1, family = binomial(link = "probit"))
phi.init <- c(coef(glm.out), mlogit(normalize(1:J)), log(1))

options(OMIR.optim.method = "BFGS")
options(OMIR.optim.control = list(trace = 5))
mle.out <- fit.mixture.link.x.mle(y, m, X, J=J, method = "betaapprox",
	invlink.p = pnorm, phi.init = phi.init)
print(mle.out)

# ----- Starting value for MCMC -----
hyper <- list(var.Beta = 1000,
	alpha.Pi = rep(1,J),
	a.kappa = 1,
	b.kappa = 3/4)

proposal <- list(
	var = diag(1, d+J-1+1),
	scale = 0.04
)

param.grp <- rep(1, d+J-1+1)

fit.pilot.out <- rwmetrop.mixlink.binomial(y, X, m, R = 1000, burn = 0, thin = 1,
	invlink = pnorm, Beta.init = mle.out$theta.hat$Beta,
	Pi.init = mle.out$theta.hat$Pi, kappa.init = mle.out$theta.hat$kappa,
	hyper = hyper, use.laplace.approx = TRUE, reporting.period = 50,
	proposal = proposal, param.grp = param.grp)

# ----- MCMC with refined proposal -----
param.grp <- seq(1, d+J-1+1)
# param.grp <- rep(1, d+J-1+1)

proposal <- list(
	var = var(fit.pilot.out$par.hist),
	scale = 1
)

fit.out <- rwmetrop.mixlink.binomial(y, X, m, R = 10000, burn = 0, thin = 1,
	invlink = pnorm, Beta.init = as.numeric(tail(fit.pilot.out$Beta.hist, 1)),
	Pi.init = as.numeric(tail(fit.pilot.out$Pi.hist,1)),
	kappa.init = as.numeric(tail(fit.pilot.out$kappa.hist,1)),
	hyper = hyper, use.laplace.approx = FALSE, reporting.period = 100,
	proposal = proposal, param.grp = param.grp)

# save.image("results.Rdata")
# load("results.Rdata")

theta.mcmc <- data.frame(
	Beta = fit.out$Beta.hist,
	Pi = fit.out$Pi.hist,
	kappa = fit.out$kappa.hist
)

# Experimental: Let's apply a post-MCMC ident. constraint on Pis
theta.mcmc[,seq(d+1,d+J)] <- t(apply(theta.mcmc[,seq(d+1,d+J)], 1, sort))


n <- length(y)
R <- nrow(theta.mcmc)
nn <- c(sprintf("Beta%d", 1:ncol(X)), sprintf("Pi%d", 1:J), "kappa")
colnames(theta.mcmc) <- nn
lab <- c(colnames(X), sprintf("Pi%d", 1:J), "kappa")

idx.Beta <- grep("^Beta", nn)
idx.Pi <- grep("^Pi", nn)
idx.kappa <- grep("^kappa", nn)

# ----- Compute DIC -----
D <- numeric(R)

for (r in 1:R) {
	Beta <- as.numeric(theta.mcmc[r, idx.Beta])
	Pi <- as.numeric(theta.mcmc[r, idx.Pi])
	kappa <- theta.mcmc[r, idx.kappa]
	prob <- pnorm(X %*% Beta)

	log.fx <- numeric(n)
	for (i in 1:n) {
		log.fx[i] <- d.mixlink.binom.one(y[i], m[i], prob[i], Pi, kappa,
			log = TRUE)
	}

	D[r] <- -2 * sum(log.fx)
}

theta.bar <- colMeans(theta.mcmc)
prob.bar <- pnorm(X %*% theta.bar[idx.Beta])
Pi.bar <- theta.bar[idx.Pi]
kappa.bar <- theta.bar[idx.kappa]
log.fx <- numeric(n)
for (i in 1:n) {
	log.fx[i] <- d.mixlink.binom.one(y[i], m[i], prob.bar[i], Pi.bar, kappa.bar,
		log = TRUE)
}
D.hat <- -2 * sum(log.fx)
D.bar <- mean(D)
p.D <- D.bar - D.hat
dic <- D.bar + p.D

# ----- Check diagnostics -----
for (j in 1:ncol(theta.mcmc))
{
	pdf(sprintf("%s/trace-%s.pdf", "plots", nn[j]), width = 5, height = 5)
	plot(theta.mcmc[,j], type = "l", xlab = "Iteration", ylab = lab[j])
	dev.off()
	
	pdf(sprintf("%s/acf-%s.pdf", "plots", nn[j]), width = 5, height = 5)
	acf(theta.mcmc[,j], main = lab[j])
	dev.off()
	
	pdf(sprintf("%s/hist-%s.pdf", "plots", nn[j]), width = 5, height = 5)
	hist(theta.mcmc[,j], xlab = lab[j], main = "")
	box()
	dev.off()
}

# ----- Print summary -----
DF <- data.frame(
	mean = colMeans(theta.mcmc),
	sd = apply(theta.mcmc, 2, sd),
	"Pct2.5" = apply(theta.mcmc, 2, quantile, prob = 0.025),
	"Pct50" = apply(theta.mcmc, 2, quantile, prob = 0.5),
	"Pct97.5" = apply(theta.mcmc, 2, quantile, prob = 0.975)
)
rownames(DF) <- lab

print(round(DF, 6))
printf("--\n")
printf("Saved MCMC samples: %d\n", nrow(theta.mcmc))
printf("DIC: %0.04f\n", dic)
printf("Elapsed %f sec\n", fit.out$elapsed.sec)

# ----- Compute predictions and PIs from posterior predictive distn -----

r.sample <- function(Beta, Pi, kappa) {
	prob <-  pnorm(X %*% Beta)
	r.mixlink.binom(n, prob, Pi, kappa, m)
}

Y.new <- matrix(NA, R, n)
for (r in 1:R) {
	Beta <- as.numeric(theta.mcmc[r, idx.Beta])
	Pi <- as.numeric(theta.mcmc[r, idx.Pi])
	kappa <- as.numeric(theta.mcmc[r, idx.kappa])
	Y.new[r,] <- r.sample(Beta, Pi, kappa)
}

Y.summary <- data.frame(
	mean = apply(Y.new, 2, mean, na.rm = TRUE),
	lo = apply(Y.new, 2, quantile, prob = 0.025, na.rm = TRUE),
	hi = apply(Y.new, 2, quantile, prob = 0.975, na.rm = TRUE)
)
print(cbind(Y.summary, y, m, x))

# ----- Plot predictions -----
pdf("plots/predict-prob.pdf", width = 7, height = 7)
plot(x, y/m, ylim = c(0,1), col = "darkgrey")
points(x, Y.summary$mean / m, pch = 20, cex = 0.5, col = "blue")
points(x, Y.summary$lo / m, pch = 20, cex = 0.5, col = "green")
points(x, Y.summary$hi / m, pch = 20, cex = 0.5, col = "green")
dev.off()

pdf("plots/predict-probit.pdf", width = 7, height = 7)
plot(x, qnorm(y/m), ylim = c(-3,3), col = "darkgrey")
points(x, qnorm(Y.summary$mean / m), pch = 20, cex = 0.5, col = "blue")
points(x, qnorm(Y.summary$lo / m), pch = 20, cex = 0.5, col = "green")
points(x, qnorm(Y.summary$hi / m), pch = 20, cex = 0.5, col = "green")
dev.off()

# ----- Plot Randomized Quantile Residuals -----
rqres.pp <- matrix(NA, R, n)
for (r in 1:R) {
	Beta <- as.numeric(theta.mcmc[r, idx.Beta])
	Pi <- as.numeric(theta.mcmc[r, idx.Pi])
	kappa <- as.numeric(theta.mcmc[r, idx.kappa])
	prob <- pnorm(X %*% Beta)
	rqres.pp[r,] <- rqres.mixlink.binom(y, m, prob, Pi, kappa)
}

rqres.post <- colMeans(rqres.pp)

pdf("plots/rqres.pdf", width = 7, height = 7)
par(mfrow = c(2,2))
qqnorm(rqres.post); abline(c(0,1), lty = 2, col = "red", lwd = 2)
hist(rqres.post)
boxplot(rqres.post)
plot(Y.summary$mean / m, rqres.post, xlab = "Predicted Proportion", ylab = "Residual")
dev.off()
