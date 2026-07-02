library(cic.newassumptions.newvarianceestimator)

set.seed(42)
d <- sim_dgp(2000)

# Run estimation; lightweight warnings will be produced for input issues.

fit <- cic(d$Y, d$X, d$Z, method = c("no-split", "split", "kde"))

fit
fit$ci
summary(fit)
