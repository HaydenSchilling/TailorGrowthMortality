### Daily modelling Bayesian

mydata <- read.csv("Github Code and Data/Data/Final Adjusted Age Data for growth modelling.csv", header = T)
alldata <- mydata
mydata <- mydata[ which(mydata$AgeClass > 0),]

head(mydata)

grdata <- mydata
grdata$id <- mydata$ID
grdata$age <- mydata$AgeClass
grdata$len <- mydata$FL..cm.

daily_data <- grdata[ which(grdata$AgeClass < 1),]
daily_data$age <- daily_data$age*365 # to get days not years
daily_data$len <- daily_data$len*10

### logistic
library(brms)
hillprior <- c(
  set_prior("normal(20.5, 0.1)", nlpar = "Linf", lb=0),
  set_prior("normal(0.02, 0.001)", nlpar = "kappa"),
  set_prior("normal(95, 5)", nlpar = "t3"))#, 
# set_prior("normal(0.05, 0.2)", class="sigma"))


logB <- bf(len/10 ~ (Linf*1/(1 + exp(-kappa*(age - t3)))),
           Linf + kappa + t3 ~ 1, 
           # Nonlinear fit
           nl = TRUE)
get_prior(logB, data=daily_data)

hill_bayes_fitLG <- brm(
  logB,
  family=gaussian(), 
  data = daily_data,
  prior = hillprior ,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  iter=50000,
  seed=1,
  thin = 5,
  file="brmlogB2.rds")

summary(hill_bayes_fitLG)
hill_bayes_fitLG <- add_criterion(hill_bayes_fitLG, criterion = "loo")
#plot(hill_bayes_fitLG)
hill_bayes_fitLG$criteria$loo
bayes_R2(hill_bayes_fitLG)
#waic(hill_bayes_fitLG)


### Try linear
library(brms)
hillpriorLM <- c(
  set_prior("normal(0.8, 1)", nlpar = "beta", lb=0),
  set_prior("normal(0.8, 1)", nlpar = "int", lb=0))#, 
# set_prior("normal(0.05, 0.2)", class="sigma"))


LMB <- bf(len/10 ~ (beta*age + int),
           beta + int ~ 1, 
           # Nonlinear fit
           nl = TRUE)
get_prior(LMB, data=daily_data)

hill_bayes_fitLM <- brm(
  LMB,
  family=gaussian(), 
  data = daily_data,
  prior = hillpriorLM ,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  iter=50000,
  seed=1,
  thin = 5,
  file="brmLMB3.rds")

summary(hill_bayes_fitLM)
hill_bayes_fitLM <- add_criterion(hill_bayes_fitLM, criterion = "loo")

#plot(hill_bayes_fitLM)

#waic(hill_bayes_fitLM)

plot(conditional_effects(hill_bayes_fitLM), points=T)
bayes_R2(hill_bayes_fitLM)


### Power
library(brms)
hillpriorP <- c(
  set_prior("normal(-5, 1)", nlpar = "a0"),
            set_prior("normal(1.2, 1)", nlpar = "a1"),
                      set_prior("normal(0.5, 1)", nlpar = "b"))#, 
# set_prior("normal(0.05, 0.2)", class="sigma"))


PB <- bf(len/10 ~ (a0 + a1 * (age^b)),
          a0 + a1 + b ~ 1, 
          # Nonlinear fit
          nl = TRUE)
get_prior(PB, data=daily_data)

hill_bayes_fitPB <- brm(
  PB,
  family=gaussian(), 
  data = daily_data,
  prior = hillpriorP ,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  iter=50000,
  seed=1,
  thin = 5,
  file="brmPB.rds")

summary(hill_bayes_fitPB)
#plot(hill_bayes_fitPB)
hill_bayes_fitPB <- add_criterion(hill_bayes_fitPB, criterion = "loo")

#loo(hill_bayes_fitPB)

#plot(conditional_effects(hill_bayes_fitPB), points=T)

### Schnute 1
schnute1 <- function(age,alpha,beta,a1,a2,y1,y2)
{
  #function to return the fitted length from a schnute type 1 relationship
  t1 <- 1-exp(-alpha*(age-a1))
  t2 <- 1-exp(-alpha*(a2-a1))
  t3 <- y2^beta - y1^beta
  t4 <- y1^beta + t3*t1/t2
  return ((y1^beta + (y2^beta - y1^beta)*(1-exp(-alpha*(age-a1)))/(1-exp(-alpha*(a2-a1))))^(1/beta))
}

hillpriorS1 <- c(
  set_prior("normal(5, 3)", nlpar = "y1", lb=0),
  set_prior("normal(20, 5)", nlpar = "y2", lb=0),
  set_prior("normal(0.5, 1)", nlpar = "beta"),
  set_prior("normal(0.5, 1)", nlpar = "alpha"))#, 
# set_prior("normal(0.05, 0.2)", class="sigma"))

# age 1 - 50, age 2 = 180
S1 <- bf(len/10 ~ ((y1^beta + (y2^beta - y1^beta)*(1-exp(-alpha*(age-50)))/(1-exp(-alpha*(180-50))))^(1/beta)),
         y1 + y2 + beta + alpha ~ 1, 
         # Nonlinear fit
         nl = TRUE)
get_prior(S1, data=daily_data)

hill_bayes_fitS1 <- brm(
  S1,
  family=gaussian(), 
  data = daily_data,
  prior = hillpriorS1 ,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  iter=50000,
  seed=1,
  thin = 5,
  file="brmS1.rds")

summary(hill_bayes_fitS1)
hill_bayes_fitS1 <- add_criterion(hill_bayes_fitS1, criterion = "loo")

#plot(hill_bayes_fitS1)

#loo(hill_bayes_fitS1)

#plot(conditional_effects(hill_bayes_fitS1), points=T)

### Schnute 2
hillpriorS2 <- c(
  set_prior("normal(5, 3)", nlpar = "y1", lb=0),
  set_prior("normal(20, 5)", nlpar = "y2", lb=0),
  #set_prior("normal(0.5, 1)", nlpar = "beta"),
  set_prior("normal(0.5, 1)", nlpar = "alpha"))#, 
# set_prior("normal(0.05, 0.2)", class="sigma"))

# age 1 - 50, age 2 = 180
S2 <- bf(len/10 ~ (y1*exp(log(y2/y1)*(1-exp(-alpha*(age-50)))/
                             (1-exp(-alpha*(180-50))))),
         y1 + y2 + alpha ~ 1, 
         # Nonlinear fit
         nl = TRUE)
get_prior(S2, data=daily_data)

hill_bayes_fitS2 <- brm(
  S2,
  family=gaussian(), 
  data = daily_data,
  prior = hillpriorS2 ,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  iter=50000,
  seed=1,
  thin = 5,
  file="brmS2b.rds")

summary(hill_bayes_fitS2)
hill_bayes_fitS2 <- add_criterion(hill_bayes_fitS2, criterion = "loo")

#plot(hill_bayes_fitS2)

#loo(hill_bayes_fitS2)

#plot(conditional_effects(hill_bayes_fitS2), points=T)


### Schnute 3
hillpriorS3 <- c(
  set_prior("normal(5, 3)", nlpar = "y1", lb=0),
  set_prior("normal(20, 5)", nlpar = "y2", lb=0),
  set_prior("normal(0.5, 1)", nlpar = "beta"))#,
  #set_prior("normal(0.5, 1)", nlpar = "alpha"))#, 
# set_prior("normal(0.05, 0.2)", class="sigma"))

# age 1 - 50, age 2 = 180
S3 <- bf(len/10 ~ ((y1^beta + (y2^beta - y1^beta)*((age-50)/(180-50)))^(1/beta)),
         y1 + y2 + beta ~ 1, 
         # Nonlinear fit
         nl = TRUE)
get_prior(S3, data=daily_data)

hill_bayes_fitS3 <- brm(
  S3,
  family=gaussian(), 
  data = daily_data,
  prior = hillpriorS3 ,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  iter=50000,
  seed=1,
  thin = 5,
  file="brmS3.rds")

summary(hill_bayes_fitS3)
hill_bayes_fitS3 <- add_criterion(hill_bayes_fitS3, criterion = "loo")

#plot(hill_bayes_fitS3)

#loo(hill_bayes_fitS3)

#plot(conditional_effects(hill_bayes_fitS3), points=T)


### Schnute 4

hillpriorS4 <- c(
  set_prior("normal(5, 3)", nlpar = "y1", lb=0),
  set_prior("normal(20, 5)", nlpar = "y2", lb=0))#,
  #set_prior("normal(0.5, 1)", nlpar = "beta"),
  #set_prior("normal(0.5, 1)", nlpar = "alpha"))#, 
# set_prior("normal(0.05, 0.2)", class="sigma"))

# age 1 - 50, age 2 = 180
S4 <- bf(len/10 ~ (y1*exp((log(y2/y1))*((age-50)/(180-50)))),
         y1 + y2  ~ 1, 
         # Nonlinear fit
         nl = TRUE)
get_prior(S4, data=daily_data)

hill_bayes_fitS4 <- brm(
  S4,
  family=gaussian(), 
  data = daily_data,
  prior = hillpriorS4 ,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  iter=50000,
  seed=1,
  thin = 5,
  file="brmS4.rds")

summary(hill_bayes_fitS4)
hill_bayes_fitS4 <- add_criterion(hill_bayes_fitS4, criterion = "loo")

#plot(hill_bayes_fitS4)

#loo(hill_bayes_fitS4)

#plot(conditional_effects(hill_bayes_fitS4), points=T)


### Von Bertalaffy
hillpriorVB <- c(
  set_prior("normal(25, 5)", nlpar = "Linf", lb=0),
  set_prior("normal(0, 5)", nlpar = "t0"),
set_prior("normal(0.5, 1)", nlpar = "K", lb=0))#,
#set_prior("normal(0.5, 1)", nlpar = "alpha"))#, 
# set_prior("normal(0.05, 0.2)", class="sigma"))

# age 1 - 50, age 2 = 180
vb <- bf(len/10 ~ (Linf*(1-exp(-K*(age-t0)))),
         Linf + K + t0  ~ 1, 
         # Nonlinear fit
         nl = TRUE)
get_prior(vb, data=daily_data)

hill_bayes_fitVB <- brm(
  vb,
  family=gaussian(), 
  data = daily_data,
  prior = hillpriorVB ,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  iter=50000,
  seed=1,
  thin = 5,
  file="brmVB2.rds")

summary(hill_bayes_fitVB)
hill_bayes_fitVB <- add_criterion(hill_bayes_fitVB, criterion = "loo")

#plot(hill_bayes_fitVB)

#loo(hill_bayes_fitVB)

#plot(conditional_effects(hill_bayes_fitVB), points=T)


loo_compare(hill_bayes_fitVB, hill_bayes_fitS1,hill_bayes_fitS2,hill_bayes_fitS3,
    hill_bayes_fitS4, hill_bayes_fitLG, hill_bayes_fitLM, hill_bayes_fitPB)
