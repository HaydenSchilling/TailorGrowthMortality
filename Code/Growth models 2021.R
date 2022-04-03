# 2021 Growth modelling
library(tidyverse)
source("../../R/functions.R")
library(FSA)
library(nlstools)
library(MuMIn)
library(minpack.lm)
library(nnet)
library(nlme)

mydata <- read_csv("Github Code and Data/Data/2021 growth modelling data.csv")
mydata$AgeYears <- mydata$AgeMonths/12

# drop big old fish
mydata <- subset(mydata, AgeClass < 11)


vbTypical <- FL..cm.~Linf*(1-exp(-K*(AgeYears-t0)))
fitTypical <- nls(vbTypical,data=mydata,start=c(Linf = 80, K = 0.3, t0 = 0))
fitPlot(fitTypical,xlab="Age",ylab="Fork Length (cm)",main="")

mydata$AgeClass2 = mydata$AgeYears >=1
weights = mydata %>% group_by(AgeClass2) %>%summarise(weights1 = var(FL..cm.))
dat = left_join(mydata, weights) 

fitTypicalgnls <- gnls(vbTypical,data=mydata,start=c(Linf = 80, K = 0.3, t0 = 0),
                       control=list(maxiter=100000000,nlsTols=100, minScale=10^-6),
                       weights = varIdent(form=~1|AgeClass2))



VB_pars <- coef(fitTypical)
VB_pars
Linf <- VB_pars[1]
K <- VB_pars[2]
t0 <- VB_pars[3]
ageX <- seq(0,7, by = 0.01)
est_VB <- (Linf*(1-exp(-K*(ageX-t0))))                   
AIC(vbTypical)

## AIC
k <- length(coef(fitTypical))
k

n <- nrow(mydata)
n

VB_AIC <- AIC(fitTypical, k = 3)
VB_AICc <- (VB_AIC + ((2 * (k + 1))/(n - k - 1)))
VB_AICc

## Schnute
grdata <- mydata
grdata$id <- mydata$ID
grdata$age <- mydata$AgeYears
grdata$len <- mydata$FL..cm.
#grdata <- subset(grdata, AgeClass < 11)

# Schnute functions
schnute1 <- function(age,alpha,beta,a1,a2,y1,y2)
{
  #function to return the fitted length from a schnute type 1 relationship
  t1 <- 1-exp(-alpha*(age-a1))
  t2 <- 1-exp(-alpha*(a2-a1))
  t3 <- y2^beta - y1^beta
  t4 <- y1^beta + t3*t1/t2
  return (t4^(1/beta))
}

schnute2 <- function(age,alpha,a1,a2,y1,y2)
{
  #function to return the fitted length from a schnute type 2 relationship
  t1 <- 1-exp(-alpha*(age-a1))
  t2 <- 1-exp(-alpha*(a2-a1))
  t3 <- log(y2/y1)
  return (y1*exp(t3*t1/t2))
}

schnute3 <- function(age,beta,a1,a2,y1,y2)
{
  #function to return the fitted length from a schnute type 3 relationship
  t1 <- (age-a1)/(a2-a1)
  t3 <- y2^beta - y1^beta
  t4 <- y1^beta + t3*t1
  return (t4^(1/beta))
}

schnute4 <- function(age,a1,a2,y1,y2)
{
  #function to return the fitted length from a schnute type 3 relationship
  t1 <- (age-a1)/(a2-a1)
  t2 <- log(y2/y1)
  return (y1*exp(t2*t1))
}

Gompetrz1 <- function(age, Linf, kappa, lamda)
{
  t2 <- (((log(lamda) - log(kappa)))/kappa)
  return (Linf*exp(-exp(-kappa*(age - t2))))
}


### Run
#To run the Schnute model #1, first we set some starting values, the two ages are pretty robust

a1 = 1 # first age to model from
a2 = 4 # 2nd age to model from
y1I = 30 # Suggest average size at a1, don't do smaller than smallest fish at a1
y2I= 50 # suggest a size appropriate to a2

# Run it, modified to run
lower <- c(0.15, -4, 15, 30)
upper <- c(0.7, 0, 50, 70)

sch1 <- nls(len ~ schnute1(age,alpha,beta,a1,a2,y1,y2), start = list(alpha=0.5,beta=-2, y1=y1I,y2=y2I),
            control=nls.control(maxiter=10000, warnOnly = T), data=grdata, algorithm="port", lower=lower, upper=upper)
summary(sch1)
# Extract Parameters
Sch1_Pars <- coef(sch1)

## CI and AIC
Boot_CIs_sch1 <- nlsBoot(sch1)
Boot_CIs_sch1$bootCI

summary(Boot_CIs_sch1)

Sch1_AIC <- AIC(sch1)
Sch1_AIC

k <- length(coef(sch1))
k

n <- nrow(grdata)
n

Sch1_AIC <- AIC(sch1, k = k)
Sch1_AICc <- (Sch1_AIC + ((2 * (k + 1))/(n - k - 1)))
Sch1_AICc


### Schnute 2
a1 = 1
a2 = 4
y1I = 30 # Suggest average size at a1, don't do smaller than smallest fish at a1
y2I= 60 # suggest a size appropriate to a2


#Run the model and get parameters

sch2 <- nls(len ~ schnute2(age,alpha,a1,a2,y1,y2), start = list(alpha=0.1,y1=y1I,y2=y2I),
            control=nls.control(maxiter=10000), data=grdata)
summary(sch2)
# Extract Parameters
Sch2_Pars <- coef(sch2)
Sch2_Pars


#Get Confidence Intervals using bootstrapping and get AIC

Boot_CIs_sch2 <- nlsBoot(sch2)
summary(Boot_CIs_sch2)


k <- length(coef(sch2))
k

n <- nrow(grdata)
n

Sch2_AIC <- AIC(sch2, k = k)
Sch2_AICc <- (Sch2_AIC + ((2 * (k + 1))/(n - k - 1)))
Sch2_AICc


### Schnute 3
a1 = 1 
a2 = 4
y1I = 24 # Suggest average size at a1, don't do smaller than smallest fish at a1
y2I= 47 # suggest a size appropriate to a2

#library(minpack.lm) # use if convergence fails below

grdata2 <- subset(grdata, age < 13)

sch3 <- nls(len ~ schnute3(age,beta,a1,a2,y1,y2), start = list(beta=1,y1=y1I,y2=y2I), #0.5 for no 0.5 age, use nlsLM if convergence fails
              control=nls.control(maxiter=10000), data=grdata2)

#vbTypical <- FL..cm.~Linf*(1-exp(-K*(AgeYears-t0)))
schnute3 <- function(age,beta,a1,a2,y1,y2)
{
  #function to return the fitted length from a schnute type 3 relationship
  t1 <- (age-a1)/(a2-a1)
  t3 <- y2^beta - y1^beta
  t4 <- y1^beta + t3*t1
  return (t4^(1/beta))
}
schnute3g <- len ~ (y1^beta + (y2^beta - y1^beta)*((AgeYears-1)/(4-1)))^(1/beta)

fitSchnute3gnls <- gnls(schnute3g, start = list(beta=2.22,y1=y1I,y2=y2I), #0.5 for no 0.5 age, use nlsLM if convergence fails
                        data=grdata2,  control=list(maxiter=100000000,nlsTols=10000, minScale=10^-100),
                        weights = varIdent(form=~1|AgeClass2))

summary(sch3)
summary(fitSchnute3gnls)

Boot_CIs_sch3 <- nlsBoot(sch3)
summary(Boot_CIs_sch3)

k <- length(coef(sch3))
k

n <- nrow(grdata)
n

Sch3_AIC <- AIC(sch3, k = k)
Sch3_AIC
Sch3_AICc <- (Sch3_AIC + ((2 * (k + 1))/(n - k - 1)))
Sch3_AICc

fitSchnute3gnls <- gnls(len ~ schnute3(age,beta,a1,a2,y1,y2), start = list(beta=1,y1=y1I,y2=y2I), #0.5 for no 0.5 age, use nlsLM if convergence fails
                        data=grdata2,  control=list(maxiter=100000000,nlsTols=100, minScale=10^-6),
                       weights = varIdent(form=~1|AgeClass2))


## Schnute 4
a1 = 1 
a2 = 4
y1I = 30 # Suggest average size at a1, don't do smaller than smallest fish at a1
y2I= 50 # suggest a size appropriate to a2

sch4 <- nls(len ~ schnute4(age,a1,a2,y1,y2), start = list(y1=y1I,y2=y2I),
            control=nls.control(maxiter=10000), data=grdata)
summary(sch4)

Boot_CIs_sch4 <- nlsBoot(sch4)
summary(Boot_CIs_sch4)

k <- length(coef(sch4))
k

n <- nrow(grdata)
n

Sch4_AIC <- AIC(sch4, k = k)
Sch4_AICc <- (Sch4_AIC + ((2 * (k + 1))/(n - k - 1)))
Sch4_AICc


### Linear
lin_fit <- lm(len ~ age, data = grdata)
summary(lin_fit)


Boot_CIs_lin <- confint(lin_fit)
Boot_CIs_lin

k <- length(coef(lin_fit))
k

n <- nrow(grdata)
n

LM_AIC <- AIC(lin_fit, k = k)
LM_AICc <- (LM_AIC + ((2 * (k + 1))/(n - k - 1)))
LM_AICc


# Gompertz Curve
#
#
#
#

Linf = 100 
kappa = 0.2
lamda = 0.4

gomp_fit <- nls(len ~ Gompetrz1(age = age, Linf = Linf, kappa = kappa, lamda = lamda), start = list(Linf = Linf, kappa=kappa, lamda = lamda), control=nls.control(maxiter=10000), data=grdata)

summary(gomp_fit)


#To get Confidence Intervals using bootstrapping and get AIC

Boot_CIs_gomp <- nlsBoot(gomp_fit)
summary(Boot_CIs_gomp)

Gomp_AIC <- AIC(gomp_fit)
Gomp_AIC

k <- length(coef(gomp_fit))
k

n <- nrow(grdata)
n

Gomp_AIC <- AIC(gomp_fit, k = k)
Gomp_AICc <- (Gomp_AIC + ((2 * (k + 1))/(n - k - 1)))
Gomp_AICc


## LOGISTIC
Linf = 91
kappa = 0.5
t3 = 3


### Logistic
log_fit <- nls(len ~ (Linf*1/(1 + exp(-kappa*(age - t3)))), start = list(Linf = Linf, kappa=kappa, t3 = t3), control=nls.control(maxiter=10000), data=grdata)

log_fit
summary(log_fit)


#Extract Parameters

Log_Pars <- coef(log_fit)
Log_Pars


#To get Confidence Intervals using bootstrapping and get AIC

Boot_CIs_log <- nlsBoot(log_fit)
summary(Boot_CIs_log)

k <- length(coef(log_fit))
k

n <- nrow(grdata)
n

Log_AIC <- AIC(log_fit, k = k)
Log_AICc <- (Log_AIC + ((2 * (k + 1))/(n - k - 1)))
Log_AICc


### Power Curve
a0 = 0
a1 = 25
b = 0.4

power_fit <- nls(len ~ (a0 + a1 * (age^b)), start = list(a0 = a0, a1=a1, b = b), control=nls.control(maxiter=10000), data=grdata)

power_fit
summary(power_fit)


#Extract Parameters

Power_Pars <- coef(power_fit)
Power_Pars


#To get Confidence Intervals using bootstrapping and get AIC

Boot_CIs_power <- nlsBoot(power_fit)
summary(Boot_CIs_power)

Power_AIC <- AIC(power_fit)
Power_AIC

k <- length(coef(power_fit))
k

n <- nrow(grdata)
n

Power_AIC <- AIC(power_fit, k = k)
Power_AICc <- (Power_AIC + ((2 * (k + 1))/(n - k - 1)))
Power_AICc


### Compare all models
model_compare <- data.frame(VB_AICc, Sch1_AICc, Sch2_AICc, Sch3_AICc, Sch4_AICc, LM_AICc, Gomp_AICc, Log_AICc, Power_AICc) # removed Sch1_AIC as model would not work

sorted_models <- sort(model_compare, decreasing = FALSE)
print(paste("The best model to use based upon AICc is", names(sorted_models[1])))

#mod_list <- list(gomp_fit, lin_fit, log_fit, power_fit, sch1, sch2, sch3, sch4, svTypical)
model_compare
model_compare-20443.62


####
#### brms attempt
library(brms)
hillprior <- c(
  set_prior("normal(2.2, 0.1)", nlpar = "beta", lb=0),
  set_prior("normal(24, 1)", nlpar = "y1"),
  set_prior("normal(47, 1)", nlpar = "y2"))#, 
# set_prior("normal(0.05, 0.2)", class="sigma"))
Schnute3b <- bf(len ~ (y1^beta + (y2^beta - y1^beta)*((AgeYears-1)/(4-1)))^(1/beta),
                beta + y1 + y2 ~ 1, 
                sigma ~ AgeClass2,
                # Nonlinear fit
                nl = TRUE)
#get_prior(Schnute3b, data=grdata2)

hill_bayes_fit <- brm(
  Schnute3b,
  family=gaussian(), 
  data = grdata2,
  prior = hillprior ,
  control =
    list(adapt_delta = 0.9),
  cores=4,
  #iter=0000,
  seed=1,
  thin = 5,
  file="brmSchnute3e.rds") # e was good

summary(hill_bayes_fit)
plot(hill_bayes_fit, )

library(rstan)
rstan::traceplot(hill_bayes_fit, pars = c("b_sigma_Intercept"), inc_warmup = TRUE)
launch_shinystan(hill_bayes_fit)
### plot
plot(conditional_effects(hill_bayes_fit, points=TRUE)) 
#plot(predict(hill_bayes_fit))


### Plot line + 95% CI interval + 95% Prediction Interval

predtimes <- seq(min(grdata2$AgeYears),max(grdata2$AgeYears),by=0.01)
hill_bayes_fitted <- fitted(hill_bayes_fit, 
                            newdata=list(AgeYears = predtimes,
                                         AgeClass2=TRUE)) %>% 
  as_tibble()
hill_bayes_pred <- predict(hill_bayes_fit,
                           newdata=list(AgeYears = predtimes,
                                        AgeClass2=TRUE)) %>%  
  as_tibble()
hill_bayes_ribbons <- tibble(
  AgeYears = predtimes,
  parentFraction=hill_bayes_fitted$Estimate,
  Estimate = hill_bayes_fitted$Estimate,
  pred_lower = hill_bayes_pred$Q2.5,
  pred_upper = hill_bayes_pred$Q97.5,
  fitted_lower = hill_bayes_fitted$Q2.5,
  fitted_upper = hill_bayes_fitted$Q97.5)
ggplot(data=grdata2, aes(x=AgeYears)) +
  geom_point(aes(y=`len`),size=2, alpha=0.7) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=pred_lower, ymax=pred_upper), 
              alpha=0.2)+#, fill=colourcodes[3]) +
  geom_ribbon(data=hill_bayes_ribbons, aes(ymin=fitted_lower, ymax=fitted_upper), 
              alpha=0.5)+#, fill=colourcodes[3]) +
  geom_line(data=hill_bayes_ribbons, aes(y=Estimate), colour="blue", 
            size=1)+
  theme_classic() +
  scale_x_continuous(breaks=seq(0,8))+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12, colour="black"),
        axis.ticks = element_line(colour="black")) +
  labs(y="Fork Lenth (cm)", x= "Age (years)")

ggsave("Output/Bayesian Growth Schnute.png", dpi=600, width=21, height=14.8, units="cm")
ggsave("Output/Bayesian Growth Schnute.pdf", dpi=600, width=21, height=14.8, units="cm")


