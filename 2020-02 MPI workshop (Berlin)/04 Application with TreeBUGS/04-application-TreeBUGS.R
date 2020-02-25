############################################################################
############ BASICS
############
############ This script matches with the examples presented in
############          "04-application-TreeBUGS.pdf"
############
############################################################################


######################################  SETUP

### install TreeBUGS from CRAN
# install.packages("TreeBUGS")


### install TreeBUGS from GitHub (newest developer version)
# install.packages(c("devtools", "coda", "runjags", "hypergeo", "testthat",
#                    "rjags", "Rcpp", "RcppArmadillo", "logspline"))
# devtools::install_github("denis-arnold/TreeBUGS", build_vignettes = TRUE)


### load TreeBUGS
library(TreeBUGS)


# adjust working directory!
# setwd("MPT-Workshop/08 Application IV - TreeBUGS")
#
# RStudio: "Session"-->"Set Working Directory"-->"To Source File Location"



###################################### DATA STRUCTURE

frequencies <- read.csv("2htm.csv")
head(frequencies, 5)


# plot example data:
plotFreq("2htm.csv", eqn = "2htm.eqn")
plotFreq(frequencies, boxplot = FALSE, eqn = "2htm.eqn")






###################################### FIT MODEL

# fitting with model files/csv files from disk:
fit_csv <- traitMPT(eqnfile = "2htm.eqn",
                    data = "2htm.csv",
                    restrictions = "2htm_constraints.txt")
fit_csv
summary(fit_csv)


# fitting in R
htm <- "
target  hit  do
target  hit  (1-do)*g
target  miss (1-do)*  (1-g)

lure  cr  dn
lure  fa  (1-dn)*g
lure  cr  (1-dn)*(1-g)
"
fit_R <- traitMPT(eqnfile = htm,
                  data = frequencies,
                  restrictions = list("dn=do"))
fit_R
summary(fit_R)


# beta-MPT (with hard-coded equality constraint):
htm_d <- "
target  hit  d
target  hit  (1-d)*g
target  miss (1-d)*  (1-g)

lure  cr  d
lure  fa  (1-d)*g
lure  cr  (1-d)*(1-g)
"
fit_beta <- betaMPT(eqnfile = htm_d, data = "2htm.csv")
fit_beta
summary(fit_beta)


# free parameter "g" for response bias
fit <- traitMPT(eqnfile = "2htm.eqn",
                data = "2htm.csv",
                restrictions = list("dn=do"))



###################################### CHECK CONVERGENCE

plot(fit, parameter = "mean", type = "default")
plot(fit, parameter = "sigma", type = "default")

# auto-correlation function (ideally close to zero):
plot(fit, parameter = "mean", type = "acf")
plot(fit, parameter = "rho", type = "acf")

# Gelman's Rhat statistic should be close to 1 (e.g., smaller than 1.05):
plot(fit, parameter = "mean", type = "gelman")



###################################### EXTEND SAMPLING

fit <- traitMPT(
  eqnfile = htm, data = frequencies,
  restrictions = list("dn=do"),
  n.adapt = 5000, # longer adaption of JAGS increases efficiency of sampling
  n.burnin = 5000,# longer burnin avoids issues due to bad starting values
  n.iter = 30000, # drawing more MCMC samples leads to higher precision
  n.thin = 10,    # ommitting every 10th sample reduces memory load
  n.chains = 4)

fit2 <- extendMPT(fit,            # fitted MPT model
                  n.adapt = 2000, # JAGS need to restart and adapt again
                  n.burnin = 0,   # burnin not needed if previous samples are OK
                  n.iter = 10000)
summary(fit)
summary(fit2)

plot(fit2)



###################################### PLOT ESTIMATES

# posterior distribution of specific parameters:
plot(fit, parameter = "mean", type = "density")

# group-level and individual MPT estimates:
plotParam(fit)
plotParam(fit, addLines = TRUE, select = c("dn", "g"))

# compare prior and posterior:
plotPriorPost(fit)

# distribution of individual MPT parameters:
plotDistribution(fit)




###################################### MODEL FIT

# graphical check of mean/covariance of frequencies:
colMeans(frequencies)
plotFit(fit)

cov(frequencies)
plotFit(fit, stat = "cov")

# posterior predictive p-values:
PPP(fit, M = 1000, nCPU = 4)


