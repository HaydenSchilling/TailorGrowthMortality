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


# Schnute 2
library(brms)
hillpriorS2 <- c(
  #set_prior("normal(2.2, 2)", nlpar = "beta"),
  set_prior("normal(0.2, 0.05)", nlpar = "alpha"),
  set_prior("normal(25, 2)", nlpar = "y1"),
  set_prior("normal(48, 3)", nlpar = "y2"))#, 
# set_prior("normal(0.05, 0.2)", class="sigma"))

GrowS2 <- bf(len~(y1*exp(log(y2/y1)*(1-exp(-alpha*(age-1)))/
                           (1-exp(-alpha*(4-1))))),
             y1 + y2 + alpha ~ 1,
             sigma ~ AgeClass2,
             # Nonlinear fit
             nl = TRUE)
get_prior(GrowS2, data=grdata)

hill_bayes_fitS2 <- brm(
  GrowS2,
  family=gaussian(), 
  data = grdata,
  prior = hillpriorS2,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  iter=50000,
  seed=1,
  thin = 5,
  file="brms/brmS2_all.rds")
hill_bayes_fitS2 <- add_criterion(hill_bayes_fitS2, criterion = "loo")

