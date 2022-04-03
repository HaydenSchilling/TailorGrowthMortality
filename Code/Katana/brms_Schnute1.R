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


# Schnute 1
library(brms)
hillpriorS1 <- c(
  set_prior("normal(0, 1)", nlpar = "beta"),
  set_prior("normal(2.3, 0.5)", nlpar = "alpha"),
  set_prior("normal(25, 5)", nlpar = "y1"),
  set_prior("normal(48, 5)", nlpar = "y2"))#, 
# set_prior("normal(0.05, 0.2)", class="sigma"))

GrowS1 <- bf(len~((y1^beta + (y2^beta - y1^beta)*(1-exp(-alpha*(age-1)))/(1-exp(-alpha*(4-1))))^(1/beta)),
              y1 + y2 + beta + alpha ~ 1,
              sigma ~ AgeClass2,
              # Nonlinear fit
              nl = TRUE)
get_prior(GrowS1, data=grdata)

hill_bayes_fitS1 <- brm(
  GrowS1,
  family=gaussian(), 
  data = grdata,
  prior = hillpriorS1,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  iter=50000,
  seed=1,
  thin = 5,
  file="brms/brmS1_all.rds")
hill_bayes_fitS1 <- add_criterion(hill_bayes_fitS1, criterion = "loo")

