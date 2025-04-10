##################################
##### Setup

# Change the working directory, e.g.:
# setwd("C:/mpt")
#
# RStudio: "Session"-->"Set Working Directory"-->"To Source File Location"

# Install and load TreeBUGS from CRAN:
install.packages("TreeBUGS")

# load TreeBUGS
library(TreeBUGS)



###################################
##### Participant heterogeneity

# assess heterogeneity graphically
plotFreq("data/data_retrieval.csv", eqn = "model/2htsm.eqn")
plotFreq("data/data_encoding.csv", eqn = "model/2htsm.eqn")

# one line per individual
plotFreq("data/data_retrieval.csv", eqn = "model/2htsm.eqn",
         freq = FALSE, boxplot = FALSE)



###################################
##### Model fitting

# Fitting a latent-trait MPT with defaults (not recommended):
m.retrieval.default  <- traitMPT(eqnfile="model/2htsm.eqn",
                                 data = "data/data_retrieval.csv",
                                 restrictions = "model/restrictions.txt")
summary(m.retrieval.default)
m.retrieval.default

# Fitting a latent-trait MPT with optional arguments:
m.retrieval <- traitMPT(eqnfile="model/2htsm.eqn",
                        data = "data/data_retrieval.csv",
                        restrictions = "model/restrictions.txt",
                        covData = "data/age_retrieval.csv",
                        modelfilename = "results/2htsm_traitMPT.jags",
                        transformedParameters = list("deltaDd=D1-d1"),
                        parEstFile = "results/results_retrieval_traitMPT.txt",
                        n.chain = 4, n.iter = 50000, n.adapt = 10000,
                        n.burnin = 10000, n.thin = 10)
summary(m.retrieval)


# Fitting a beta-MPT model

m.retrieval.beta <- betaMPT(eqnfile="model/2htsm.eqn",
                            data = "data/data_retrieval.csv",
                            restrictions = "model/restrictions.txt",
                            modelfilename = "results/2htsm_betaMPT.jags",
                            transformedParameters = list("deltaDd=D1-d1"),
                            parEstFile = "results/results_retrieval_betaMPT.txt",
                            n.chain = 4, n.iter = 50000, n.adapt = 10000,
                            n.burnin = 10000, n.thin = 10)
summary(m.retrieval.beta)



###################################
##### Monitoring convergence graphically

?plot.traitMPT
plot(m.retrieval, parameter = "mean")
plot(m.retrieval, parameter = "mean", type = "acf")
plot(m.retrieval, parameter = "mean", type = "gelman")

# other parameters:
plot(m.retrieval, parameter = "sigma")
plot(m.retrieval, parameter = "rho")

# continue MCMC sampling
m.retrieval2 <- extendMPT(m.retrieval.default,
                          n.iter = 10000, n.adapt = 2000)
summary(m.retrieval2)





###################################
##### Testing model fit

# graphically:
plotFit(m.retrieval)
plotFit(m.retrieval, stat = "cov")

# posterior predictive p-values:
ppp.retrieval <- PPP(m.retrieval, M=1000)
ppp.retrieval




###################################
##### Plotting and extracting parameters

# plot parameters:
plotParam(m.retrieval, addLines = TRUE, select = c("a", "b"))

# plot estimated hierarchical distribution against plausible values
plotDistribution(m.retrieval, scale = "probability")

# compare prior to posterior density
plotPriorPost(m.retrieval)



###################################
##### Within-subject tests


# compute differences/ratios of MPT parameters within participants
transpar <- transformedParameters(fittedModel = m.retrieval,
                                  transformedParameters=list("deltaDd=D1-d1"),
                                  level = "individual")
summary(transpar[,1:5])$quantiles




###################################
##### Between-subject tests

# fit latent-trait MPT for second between-subject condition (encoding)
m.encoding <- traitMPT(eqnfile="model/2htsm.eqn",
                       data = "data/data_encoding.csv",
                       restrictions = "model/restrictions.txt",
                       modelfilename = "results/2htsm.jags",
                       covData = "data/age_encoding.csv",
                       transformedParameters = list("deltaDd=D1-d1"),
                       parEstFile = "results/results_encoding.txt",
                       n.chain = 4, n.iter = 50000, n.adapt = 10000,
                       n.burnin = 10000, n.thin = 10,
                       ppp = 5000, dic = TRUE)
summary(m.encoding)

# get (A) credibility intervals and (B) Bayesian p-values
# for the difference in parameters D and d:
betweenSubjectMPT(m.retrieval, m.encoding, par1 = "D1")



###################################
##### Including continuous predictors

m.predictor <- traitMPT(eqnfile="model/2htsm.eqn",
                        data = "data/data_retrieval.csv",
                        restrictions = "model/restrictions.txt",
                        modelfilename = "results/2htsm_predictor.jags",
                        covData = "data/pc_retrieval.csv",
                        predStructure = list("b a ; pc"),
                        parEstFile = "results/results_predictor.txt",
                        n.chain = 4, n.iter = 70000, n.adapt = 20000,
                        n.burnin = 15000, n.thin = 15)
summary(m.predictor)
# Note: unstandardized regression weights are reported!
#       (the latent-trait regression equation is:
#        theta_i = Phi( mu + slope*covariate_i + delta_i)
#     [MPT param.]    [mean]   [regression]    [person effect]



###################################
##### Including discrete predictors for between-subject manipulations

m.both.conditions <- traitMPT(eqnfile="model/2htsm.eqn",
                              data = "data/data_both.csv",
                              restrictions = "model/restrictions.txt",
                              modelfilename = "results/2htsm_between.jags",
                              covData = "data/group.csv",
                              predStructure = list("a b d1 D1 ; Group"),
                              predType = c("f"),  # fixed effects
                              parEstFile = "results/results_both.txt",
                              n.chain = 4, n.iter = 20000, n.adapt = 5000,
                              n.burnin = 5000, n.thin = 5)
summary(m.both.conditions)
round(getGroupMeans(m.both.conditions), 3)



###################################
##### Check and adjust prior distributions

# TreeBUGS default prior for latent-trait MPT
plotPrior(prior = list(mu = "dnorm(0,1)",
                       xi = "dunif(0,10)",
                       V=diag(2), df=3))

# optional: fit with different prior:
# m.retrieval.prior  <- traitMPT(eqnfile="model/2htsm.eqn",
#                                data = "data/data_retrieval.csv",
#                                restrictions = "model/restrictions.txt",
#                                mu = c( a="dnorm(0,4)", b="dnorm(0,4)",
#                                        d1="dnorm(0,1)", D1="dnorm(0,1)"))


###################################
##### Prior Predictive

pp <- priorPredictive(prior = list(mu = "dnorm(0,1)",
                                   xi = "dunif(0,10)",
                                   V=diag(4), df=4 + 1),
                      eqnfile = "model/2htsm.eqn",
                      restrictions = "model/restrictions.txt",
                      numItems = c(32, 32, 32),
                      N = 100, M = 200)

roc <- matrix(NA, nrow = 200, ncol = 2)  # matrix with 2 columns: false alarm & hit rate
for (i in 1:200){

  # false alarm rate: "old" responses to new items
  roc[i,1] <- sum(pp[[i]][,c("NE", "NU")]) / (32 * 100)

  # hit rate: "old" responses to learned items
  roc[i,2] <- sum(pp[[i]][,c("EE", "EU", "UU", "UE")]) / (32 * 2 * 100)
}
plot(roc, xlim = 0:1, ylim = 0:1, xlab = "False alarm rate", ylab = "Hit rate")
abline(0,1)
