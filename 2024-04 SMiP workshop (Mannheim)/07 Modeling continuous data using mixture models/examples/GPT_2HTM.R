################################################################################
################## Generalized processing tree (GPT) models
##################
################## Small example for the 2HTM.
##################
################################################################################

# adjust working directory
# setwd("C:/MPTworkshop/")

# install package:
# library(devtools)
# install_github("danheck/gpt")

library(gpt)


### generate data with parametric GPT model (ex-Gaussian for RTs)

gen <- gpt_gen(n = c("targets"=200, "lures"=200),
               theta = c(d = .4, g = .5),
               eta = c(shift = 400, sig = 50,
                       lambda_do = 200, lambda_dn = 350, lambda_g = 450),
               file = "GPT_2HTM_model.txt",
               latent = "exgauss",
               restrictions = list("do=dn", "lambda_go=lambda_gn"))
head(gen)
hist(gen$y)


### fit GPT model

fit <- gpt_fit("x", "y", data = gen,
               file="GPT_2HTM_model.txt",
               latent="exgauss",
               restrictions = list("do=dn", "lambda_go=lambda_gn"))
fit
hist(fit, n = 30)
plot(fit)
test_fit(fit)$test
