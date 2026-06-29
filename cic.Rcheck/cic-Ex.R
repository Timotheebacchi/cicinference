pkgname <- "cic"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('cic')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("cic")
### * cic

flush(stderr()); flush(stdout())

### Name: cic
### Title: Changes-in-Changes Estimator
### Aliases: cic

### ** Examples

set.seed(42)
d <- sim_dgp(500)
fit <- cic(d$Y, d$X, d$Z, method = "no-split")
summary(fit)



cleanEx()
nameEx("qY_dgp")
### * qY_dgp

flush(stderr()); flush(stdout())

### Name: qY_dgp
### Title: Quantile Function for Athey & Imbens (2006) DGP
### Aliases: qY_dgp

### ** Examples

qY_dgp(0.5, d1 = 0, d2 = 0.05)



cleanEx()
nameEx("sim_dgp")
### * sim_dgp

flush(stderr()); flush(stdout())

### Name: sim_dgp
### Title: Simulate Data from Athey & Imbens (2006) DGP
### Aliases: sim_dgp

### ** Examples

set.seed(42)
d <- sim_dgp(500)
head(d)



cleanEx()
nameEx("theta_true")
### * theta_true

flush(stderr()); flush(stdout())

### Name: theta_true
### Title: True CiC Parameter for Athey & Imbens (2006) DGP
### Aliases: theta_true

### ** Examples

theta_true(b1 = 0, b2 = 0.05, d1 = 0, d2 = 0.05)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
