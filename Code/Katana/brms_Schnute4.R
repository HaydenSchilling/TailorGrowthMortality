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


# Schnute 4
library(brms)
hillpriorS4 <- c(
  set_prior("normal(5, 3)", nlpar = "y1", lb=0),
  set_prior("normal(20, 5)", nlpar = "y2", lb=0))#, 
# set_prior("normal(0.05, 0.2)", class="sigma"))

GrowS4 <- bf(len~(y1*exp((log(y2/y1))*((age-1)/(4-1)))),
             y1 + y2 ~ 1,
             sigma ~ AgeClass2,
             # Nonlinear fit
             nl = TRUE)
get_prior(GrowS4, data=grdata)

hill_bayes_fitS4 <- brm(
  GrowS4,
  family=gaussian(), 
  data = grdata,
  prior = hillpriorS4,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  iter=50000,
  seed=1,
  thin = 5,
  file="brms/brmS4_all.rds")
hill_bayes_fitS4 <- add_criterion(hill_bayes_fitS4, criterion = "loo")

