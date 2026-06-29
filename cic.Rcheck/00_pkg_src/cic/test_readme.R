library(cic)

set.seed(42)
d <- sim_dgp(500)

diag <- check_cic_assumptions(d$Y, d$X, d$Z)
diag$pass_all

fit <- cic(d$Y, d$X, d$Z, method = c("no-split", "split", "kde"))

fit
fit$ci
summary(fit)
