# Bayesian Growth models
library(tidyverse)

mydata <- read_csv("Github Code and Data/Data/2021 growth modelling data.csv")
mydata$AgeYears <- mydata$AgeMonths/12

# drop big old fish
mydata <- subset(mydata, AgeClass < 11)
mydata$AgeClass2 = mydata$AgeYears >=1

grdata <- mydata
grdata$id <- mydata$ID
grdata$age <- mydata$AgeYears
grdata$len <- mydata$FL..cm.


# Linear
library(brms)
hillpriorLinear <- c(
  set_prior("normal(10, 10)", nlpar = "beta", lb=0),
  set_prior("normal(0, 10)", nlpar = "m", lb=0))

GrowLinear <- bf(len~(beta*age + m),
              beta ~ 1,
              m ~ 1,
              sigma ~ AgeClass2,
              # Nonlinear fit
              nl = TRUE)
get_prior(GrowLinear, data=grdata)

hill_bayes_fitLinear <- brm(
  GrowLinear,
  family=gaussian(), 
  data = grdata,
  prior = hillpriorLinear,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  iter=50000,
  seed=1,
  thin = 5,
  file="brmLinear_all.rds")
hill_bayes_fitLinear <- add_criterion(hill_bayes_fitLinear, criterion = "loo")
plot(hill_bayes_fitLinear)
summary(hill_bayes_fitLinear)
