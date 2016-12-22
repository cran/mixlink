library(mixlink)
set.seed(1234)

# ----- Compare histogram to density for MixLink Binomial -----
n <- 100000; m <- 20
x <- 0:m
p <- 1/4; Pi <- c(1/6, 5/6); kappa <- 2
y <- r.mixlink.binom(n, p, Pi, kappa, m)
f <- factor(y, levels=0:m)
prop <- table(f)/n
ylim <- c(0, max(prop) * 1.05)
coords <- barplot(prop, ylim = ylim, xlab = "x", ylab = "Density")
points(coords, d.mixlink.binom(x, m, p, Pi, kappa), pch = 19)
title("Histogram vs. Density")

# ----- Compare histogram to density for MixLink Poisson -----
n <- 100000; m <- 40
x <- 0:m
mean <- 10; Pi <- normalize(c(1,2,4)); kappa <- 2
y <- r.mixlink.pois(n, mean, Pi, kappa)
f <- factor(y, levels=0:m)
prop <- table(f)/n
ylim <- c(0, max(prop) * 1.05)
coords <- barplot(prop, ylim = ylim, xlab = "x", ylab = "Density")
points(coords, d.mixlink.pois(x, mean, Pi, kappa), pch = 19)
title("Histogram vs. Density")
