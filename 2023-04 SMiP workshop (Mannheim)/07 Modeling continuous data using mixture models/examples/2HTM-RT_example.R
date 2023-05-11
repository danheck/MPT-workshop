############################################################################
############
############ Modeling continuous variables using discrete bins (MPT-RT)
############
############################################################################

####  for MPT modeling in R:
# install.packages("MPTinR")
library("MPTinR")


################## EXAMPLE: how to use MPT-RT approach


### load some data: one response + one RT per row
recog <- read.csv("2HTM-data.csv")
head(recog)


##################
### (A) 2 RTs bins

# compute geometric mean:
rt_mean <- exp(mean(log(recog$rt)))

# categorize RT into bins:
recog$rt_cat <- as.numeric(recog$rt > rt_mean)
head(recog, 10)

# make new category labels:
recog$mptrt_category <- paste0(recog$response, recog$rt_cat)
head(recog)

# define all possible MPT-RT categories
mptrt_categories <- c("cr0", "cr1", "fa0", "fa1",
                      "hit0", "hit1", "miss0", "miss1")

# make a table of frequencies:
frequencies <- table(factor(x = recog$mptrt_category,
                            levels = mptrt_categories))
frequencies

# easier for copying frequencies to multiTree:
t(t(frequencies))

# fit MPT-RT (also possible in multiTree)
mpt <- fit.mpt(data = c(frequencies),
               model.filename = "2HTM-RT_2bins.txt")
mpt$goodness.of.fit
mpt$parameters




##################
### (B) 4 RTs bins

# mean and SD of log RTs
logrt_mean <- mean(log(recog$rt))
logrt_sd <- sd(log(recog$rt))

# get quantiles from lognormal (see paper)
bounds <- exp(qnorm(c(.25, .50, .75), logrt_mean, logrt_sd))
bounds

# categorize RTs into 4 bins:
recog$rt_cat <- findInterval(recog$rt, bounds)
recog$mptrt_category <- paste0(recog$response, recog$rt_cat)
head(recog)

# make all possible combinations of RT bins and MPT categories:
mptrt_categories <- c(t(outer(c("cr", "fa","hit", "miss"), 0:3, paste0)))
mptrt_categories

# get frequency table for MPT-RT:
frequencies <- table(factor(recog$mptrt_category, levels = mptrt_categories))
frequencies

mpt <- fit.mpt(c(frequencies), model.filename = "2HTM-RT_4bins.txt")
mpt$goodness.of.fit
mpt$parameters



##### how to plot estimated cumulative density function (cdf)

get_cdf <- function(parameter){
  # get parameter estimates:
  L <- mpt$parameters[paste0(parameter, 1:3),"estimates"]
  # reparameterize to bin probabilities:
  L2 <-c("bin1" = L[1]    *L[2],
         "bin2" = L[1]    *(1-L[2]),
         "bin3" = (1-L[1])*L[3],
         "bin4" = (1-L[1])*(1-L[3]))
  # make step-function for cdf:
  stepfun(x = c(min(recog$rt), bounds),
          y = cumsum(c(0, L2)))
}

plot(get_cdf("L_do"), main = "Estimated cdf", ylab = "CDF",las = 1)
plot(get_cdf("L_dn"), add = TRUE, col = 2, lty = 3)
plot(get_cdf("L_g"), add = TRUE, col = 4, lty = 2)
legend("bottomright", legend = c("do", "dn", "g"),
       col =c(1,2,4), lty =c(1,3,2))
