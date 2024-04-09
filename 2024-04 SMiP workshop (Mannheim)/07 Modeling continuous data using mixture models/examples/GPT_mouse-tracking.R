############################################################
# load packages and data
############################################################

# adjust working directory
# setwd("C:/MPTworkshop/")

install.packages("gpt_0.6.0.tar.gz", repos = NULL)
library("gpt")

# load data (Kieslich & Henninger, 2017; cf. mousetrap package)
KH <- read.csv("GPT_mouse-tracking_data.csv", stringsAsFactors = FALSE)
head(KH)


############################################################
## Feature comparison model: General model
############################################################

# general FCM version
#     Theory:          m_1 < m_2
#     Implementation:  m_2 = m_1 + m_diff   (with m_diff > 0)

fcm <- gpt_fit(x = "cat" ,                                 # discrete response
               y = "MAD",                                  # continuous variable(s)
               data = KH,                                  # data frame
               file = "GPT_mouse-tracking_model.txt",      # model file (GPT syntax)
               latent = "normal",                          # latent component distributions
               n.fit = c(5,0), EM.tol=.0001,               # options for EM algorithm
               eta.lower = c(m_diff = 0, s_1=0, s_2 = 0))  # constraint: SD > 0
fcm

hist(fcm, freq = "freq", xlab = "Maximum absolute deviation (MAD)",
     breaks = seq(-500, 1400, 70))

# goodness-of-fit
test <- test_fit(fcm, bins = 4)
test$test



############################################################
## Feature comparison model: Nested Models
############################################################

# H1: One normal distribution per item type, different accuracy across conditions
#       (setting f=1 leads to problems with the EM algorithm)
#       (alternatively, a different model file can be used that omits f_a/f_t/c_1a/c2_t)

fcm.one <- gpt_fit("cat" ,"MAD", KH, "GPT_mouse-tracking_model.txt", "normal",
                   n.fit = c(5,0), EM.tol=.0001, eta.lower = c(s_2 = 0),
                   restrictions = list("f_a=.000001","f_t=.999999",
                                       "c_1a=.999999","c_2t=.999999"))
test_nested(fcm, fcm.one)

# H2: f_t > f_a
fcm.f <- gpt_fit("cat" ,"MAD", KH, "GPT_mouse-tracking_model.txt", "normal",
                   n.fit = c(5,0), EM.tol=.0001, restr=list("f_t=f_a"),
                   eta.lower = c(m_diff = 0, s_1=0, s_2 = 0))
test_nested(fcm, fcm.f)


# H3: c_t1 > c_a1
fcm.c1 <- gpt_fit("cat" ,"MAD", KH, "GPT_mouse-tracking_model.txt", "normal",
                    n.fit = c(5,0), EM.tol=.0001, restr=list("c_1t = c_1a"),
                    eta.lower = c(m_diff = 0, s_1 = 0, s_2 = 0))
test_nested(fcm, fcm.c1)


# H4: c_t2 = c_a2
fcm.c2 <- gpt_fit("cat" ,"MAD", KH, "GPT_mouse-tracking_model.txt", "normal",
                  n.fit = c(5,0), EM.tol=.0001, restr=list("c_2t=c_2a"),
                  eta.lower = c(m_diff = 0, s_1 = 0, s_2 = 0))
test_nested(fcm, fcm.c2)
hist(fcm.c2, freq = "freq")


# AIC, BIC
msel <- select_gpt("General" = fcm, "one component" = fcm.one,
                   "f_t = f_a" = fcm.f, "c_1t = c_1a" = fcm.c1,
                   "c_2t = c_2a" = fcm.c2)
round(msel, 2)


# best model:
fcm.c2
par(mar = c(4,4,4,1), mfrow = c(2,2))
plot(fcm.c2, xlab = "Maximum deviation")

