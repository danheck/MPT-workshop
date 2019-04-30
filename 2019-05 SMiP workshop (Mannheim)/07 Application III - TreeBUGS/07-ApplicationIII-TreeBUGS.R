############################################################################
############ BASICS
############
############ This script matches with the examples presented in
############          "08-ApplicationIV-TreeBUGS.pdf"
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
library("TreeBUGS")

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





############################################################################
############ ADVANCED
############################################################################


###################################### WITHIN-SUBJECT COMPARISONS

# (1.) data file
freq_within <- read.csv("2htm_within.csv")
head(freq_within, 3)

# (2.) create EQN file for within-subject manipulations
withinSubjectEQN(htm_d,
                 labels = c("high","low"), # factor labels
                 constant=c("g"),          # parameters constrained across conditions
                 save = "2htm_within.eqn")

# (3.) fit to all conditions:
fit_within <- traitMPT("2htm_within.eqn", "2htm_within.csv")
plot(fit_within)

# (4.) compute difference in d:
diff_d <- transformedParameters(fit_within,
                                transformedParameters = list("diff_d = d_high - d_low"),
                                level = "group")
summary(diff_d)
plot(diff_d)



###################################### BETWEEN-SUBJECT COMPARISONS

# 1. fit MPT for each condition separately
fit1 <- traitMPT(htm_d, "2htm.csv")
fit2 <- traitMPT(htm_d, "2htm_group2.csv")

# 2. compute difference in parameters
diff_between <- betweenSubjectMPT(fit1, fit2,             # fitted MPT models
                                  par1 = "d",             # parameter to test
                                  stat = c("x-y","x>y"),  # transformed parameters
                                  plot = TRUE)
diff_between
plot(diff_between$mcmc)



###################################### COVARIATES: CORRELATION

# include correlation with covariates:
fit_cor <- traitMPT(htm_d, data = "2htm.csv",
                    covData = "covariates.csv")   # data with covariate(s)

# warning: posterior only quantifies uncertainty with respect to the MPT parameter estimates!
plot(fit_cor, "cor")
summary(fit_cor)
round(fit_cor$summary$group$cor, 2)

# We also need to consider the number of participants (sample size)!
correlationPosterior(fit_cor)



###################################### COVARIATES: REGRESSION

# probit regression for continuous covariate to predict MPT parameter:
fit_regression <- traitMPT(htm_d, data = "2htm.csv",
                           covData = "covariates.csv",
                           predStructure = list("d ; continuous"))
plot(fit_regression, "slope")
summary(fit_regression)

round(fit_regression$summary$group$slope, 2)


# Bayes Factor for Covariate
# * H0: Slope parameter beta=0
# * H1: Slope parameter beta ~ Cauchy(0, r) (with scale parameter r)
BayesFactorSlope(fit_regression,
                 parameter = "slope_d_continuous",
                 direction = ">",
                 plot = TRUE)



###################################### BETWEEN-SUBJECT similar to ANOVA

# Between-Subject Comparisons: Alternative method
# => Identical covariance matrix in each condition (as in ANOVA: "pooled variance")

fit_between <- traitMPT(
                  htm_d, "2htm.csv",
                  covData = "covariates.csv",
                  predStructure = list("d ; discrete"),  # discrete predictor
                  predType = c("c","f")) # "c" =continuous; "f"=fixed-effects
plot(fit_between, "factor_")
summary(fit_between)

# get estimates for the group-specific MPT parameters
gmeans <- getGroupMeans(fit_between)
round(gmeans, 2)





############################################################################
############ SIMULATION & ROBUSTNESS
############################################################################


###################################### CHANGING PRIORS


# what does the prior mean?
plotPrior(prior = list(mu = c(dn = "dnorm(0,1)",  # default prior
                              g = "dnorm(0,5)"),  # prior focused around 50%
                       xi="dunif(0,2)",           # more stable prior for scale parameter
                       V= diag(2), df = 3))       # default Wishart prior

fit_prior <- traitMPT("2htm.eqn", "2htm.csv", restrictions=list("dn=do"),
                      mu = c(dn = "dnorm(0,1)",  # default prior
                      g = "dnorm(0,5)"),         # prior focused around 50% guessing
					  xi = "dunif(0,2)",         # less disperson of MPT parameters
					  V = diag(2),               # default
					  df = 2 + 1)                # default
summary(fit_prior)
summary(fit)       # comparison to default prior



###################################### PRIOR PREDICTIVE

pp <- priorPredictive(prior = list(mu = "dnorm(0,1)", xi="dunif(0,10)",
                                   V=diag(2), df=2+1),
                      eqnfile = htm, restrictions = list("dn=do"),
                      numItems = c(target = 50, lure = 50),
                      N = 50,  M = 500)     # number of participants/samples

# compute and plot predicted values for the average hit/FA rates (in ROC space)
mean_freq <- sapply(pp, colMeans)
plot(mean_freq["fa",]/50, mean_freq["hit",]/50,
     asp = 1, xlim =0:1, ylim=0:1)
polygon(c(0,1,1), c(0,0,1), col = "gray")
abline(0, 1)



###################################### POSTERIOR PREDICTIVE

# sample new data using the posterior samples
postpred <- posteriorPredictive(fit, M = 100, nCPU = 4)

mf <- sapply(postpred, colMeans)
plot(mf["fa",]/50, mf["hit",]/50, asp = 1, xlim =c(.2,.6), ylim=c(.4,.9),
     las=1,col = adjustcolor(1, alpha=.3),
     main = "Posterior predicted (mean frequencies)",
	 xlab= "False Alarm Rate", ylab = "Hit Rate")
abline(0, 1)
polygon(c(0,1,1), c(0,0,1), col = "gray")
points(matrix(colMeans(frequencies)[c("fa", "hit")], 1)/50,
       col = 2, pch = 16, cex = .5)



###################################### SIMULATION

set.seed(123)

# standard MPT for one person (no hierarchical structure)
sim <- genMPT(theta = c(d = .7, g = .5),
              numItems = c(target = 50, lure = 50),
              eqnfile = htm_d)
sim


# hierarchical MPT
gendat <- genTraitMPT(N = 10,
                      numItems = c(target = 50, lure = 50),
                      eqnfile = htm_d,
                      mean = c(d = .6, g = .5),
                      sigma = c(d = 1, g = .5),
                      rho = diag(2))        ## get real data
gendat$data


# prior predictive:
#     1. generate parameters from prior
#     2. generate data from parameters
priorPredictive(prior = list(mu = "dnorm(0,1)",
                             xi="dunif(0,10)",
                             V=diag(2), df=2+1),
                htm_d,
                numItems = c(100, 100), level = "data",
                N = 1, M = 100, nCPU = 4)







########################## Appendix

# Testing for Heterogeneity (Smith & Batchelder, 2008)

# A) Test person heterogeneity assuming items homogeneity ($\chi^2$)
test <- testHetChi(freq = "2htm.csv",
                   tree = c(cr="lure", fa="lure",
                            hit="target",miss="target"))
data.frame(test)

# B) Test person heterogeneity under item heterogeneity (permutation bootstrap)
#    => requires data in long format (variables: person / item / response)
# testHetPerm(data, tree, source = "person")

