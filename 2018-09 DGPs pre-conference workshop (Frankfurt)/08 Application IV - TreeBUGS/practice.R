

library(TreeBUGS)

setwd("D:/R/Treebugs_Project/MPT-workshop/08 Application IV - TreeBUGS")

frequencies <- read.csv("2htm.csv")
frequencies

plotFreq(frequencies, eqnfile = "2htm.eqn", boxplot = FALSE)

fit <- traitMPT(eqnfile = "2htm.eqn",
                data = "2htm.csv",
                restrictions = list("dn=do"))

plot(fit)
plot(fit, parameter = "rho")
?TreeBUGS::plot.traitMPT

fit

summary(fit)

plotParam(fit, addLines = TRUE)

plotPriorPost(fit)

plotDistribution(fit)

plotFit(fit)

plotFit(fit, stat = "cov")

frequencies

