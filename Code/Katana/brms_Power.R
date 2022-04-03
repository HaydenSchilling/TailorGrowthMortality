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


# Power
library(brms)
hillpriorP <- c(
  set_prior("normal(0, 5)", nlpar = "a0"),
  set_prior("normal(25, 10)", nlpar = "a1"),
  set_prior("normal(0.5, 1)", nlpar = "b"))
#set_prior("normal(0.05, 0.2)", class="sigma"))
GrowP <- bf(len~(a0 + a1 * (age^b)),
             a0 + a1 + b ~ 1, 
             sigma ~ AgeClass2,
             # Nonlinear fit
             nl = TRUE)
get_prior(GrowP, data=grdata)

hill_bayes_fitP <- brm(
  GrowP,
  family=gaussian(), 
  data = grdata,
  prior = hillpriorP,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  iter=50000,
  seed=1,
  thin = 5,
  file="brms/brmPower_all.rds")
hill_bayes_fitP <- add_criterion(hill_bayes_fitP, criterion = "loo")

