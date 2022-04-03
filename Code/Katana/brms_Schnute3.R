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


# Schnute 3
library(brms)
hillpriorS3 <- c(
  set_prior("normal(2.2, 0.1)", nlpar = "beta", lb=0),
  set_prior("normal(24, 1)", nlpar = "y1"),
  set_prior("normal(47, 1)", nlpar = "y2"))#, 
# set_prior("normal(0.05, 0.2)", class="sigma"))

GrowS3 <- bf(len~(y1^beta + (y2^beta - y1^beta)*((AgeYears-1)/(4-1)))^(1/beta),
             y1 + y2 + beta ~ 1,
             sigma ~ AgeClass2,
             # Nonlinear fit
             nl = TRUE)
get_prior(GrowS3, data=grdata)

hill_bayes_fitS3 <- brm(
  GrowS3,
  family=gaussian(), 
  data = grdata,
  prior = hillpriorS3,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  iter=50000,
  seed=1,
  thin = 5,
  file="brms/brmS3_all.rds")
hill_bayes_fitS3 <- add_criterion(hill_bayes_fitS3, criterion = "loo")

