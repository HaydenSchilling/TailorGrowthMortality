# Bayesian Growth models
library(tidyverse)

mydata <- read_csv("2021 growth modelling data.csv")
mydata$AgeYears <- mydata$AgeMonths/12

# drop big old fish
mydata <- subset(mydata, AgeClass < 11)
mydata$AgeClass2 = mydata$AgeYears >=1

grdata <- mydata
grdata$id <- mydata$ID
grdata$age <- mydata$AgeYears
grdata$len <- mydata$FL..cm.


# Von bert
library(brms)
hillpriorVB <- c(
  set_prior("normal(100, 20)", nlpar = "Linf", lb=0),
  set_prior("normal(0.1, 0.1)", nlpar = "K", lb=0),
  set_prior("normal(0, 0.5)", nlpar = "t0"))#, 
#set_prior("normal(0.05, 0.2)", class="sigma"))
GrowVB <- bf(len~Linf*(1-exp(-K*(AgeYears-t0))),
             t0 + Linf + K ~ 1, 
             sigma ~ AgeClass2,
             # Nonlinear fit
             nl = TRUE)
get_prior(GrowVB, data=grdata)

hill_bayes_fitVB <- brm(
  GrowVB,
  family=gaussian(), 
  data = grdata,
  prior = hillpriorVB ,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  iter=50000,
  seed=1,
  thin = 5,
  file="brms/brmVB_all.rds")
hill_bayes_fitVB <- add_criterion(hill_bayes_fitVB, criterion = "loo")

