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


# Logistic
library(brms)
hillpriorLog <- c(
  set_prior("normal(80, 10)", nlpar = "Linf", lb=0),
  set_prior("normal(0.5, 0.2)", nlpar = "kappa", lb=0),
  set_prior("normal(3, 5)", nlpar = "t3"))

GrowLog <- bf(len~(Linf*1/(1 + exp(-kappa*(age - t3)))),
            Linf + kappa + t3 ~ 1, 
            sigma ~ AgeClass2,
            # Nonlinear fit
            nl = TRUE)
get_prior(GrowLog, data=grdata)

hill_bayes_fitLog <- brm(
  GrowLog,
  family=gaussian(), 
  data = grdata,
  prior = hillpriorLog,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  iter=50000,
  seed=1,
  thin = 5,
  file="brms/brmLog_all.rds")
hill_bayes_fitLog <- add_criterion(hill_bayes_fitLog, criterion = "loo")

