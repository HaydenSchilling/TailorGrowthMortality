# Comparing Katana growth models
library(tidyverse)
library(brms)

# Linear
linearG <- readRDS("Github Code and Data/Data/Growth models/brmLinear_all.rds")

summary(linearG)
plot(linearG) # linear is good
linearG <- add_criterion(linearG, "loo", file = "Github Code and Data/Data/Growth models/brmLinear_all")
loo(linearG)

# Logistic
logG <- readRDS("Github Code and Data/Data/Growth models/brmLog_all.rds")
summary(logG)# Logistic is good
plot(logG)
logG <- add_criterion(logG, "loo", file = "Github Code and Data/Data/Growth models/brmLog_all")
loo(logG)

# Power
PowerG <- readRDS("Github Code and Data/Data/Growth models/brmPower_all.rds")

summary(PowerG)
plot(PowerG) # Power is good
PowerG <- add_criterion(PowerG, "loo", file = "Github Code and Data/Data/Growth models/brmPower_all")
loo(PowerG)

# Von Bert
vbG <- readRDS("Github Code and Data/Data/Growth models/brmVB_all.rds")
summary(vbG) # von bert is good
plot(vbG)
vbG <- add_criterion(vbG, "loo", file = "Github Code and Data/Data/Growth models/brmVB_all")
loo(vbG)

# Schnute 1
s1G <- readRDS("Github Code and Data/Data/Growth models/brmS1_all.rds")
summary(s1G)
plot(s1G) # S1 is OK, 1 parameter Rhat = 1.01, visual mixing OK, some divergent transitions
s1G <- add_criterion(s1G, "loo", file = "Github Code and Data/Data/Growth models/brmS1_all")
loo(s1G)

# Schnute 2
s2G <- readRDS("Github Code and Data/Data/Growth models/brmS2_all.rds")
summary(s2G) # Schnute 2 is good
plot(s2G)
s2G <- add_criterion(s2G, "loo", file = "Github Code and Data/Data/Growth models/brmS2_all")
loo(s2G)

# Schnute 3
s3G <- readRDS("Github Code and Data/Data/Growth models/brmS3_all.rds")
summary(s3G)
plot(s3G) # S3 is good
s3G <- add_criterion(s3G, "loo", file = "Github Code and Data/Data/Growth models/brmS3_all")
loo(s3G)

# Schnute 4
s4G <- readRDS("Github Code and Data/Data/Growth models/brmS4_all.rds")
summary(s4G) # Schnute 4 is good
plot(s4G)
s4G <- add_criterion(s4G, "loo", file = "Github Code and Data/Data/Growth models/brmS4_all")
loo(s4G)

loo_compare(linearG, logG, PowerG, vbG, s1G, s2G, s3G, s4G)

plot(conditional_effects(s3G),points=T)


### Make Final Plot
grdata2 <- read_csv("Github Code and Data/Data/2021 growth modelling data.csv")
grdata2$AgeYears <- grdata2$AgeMonths/12
grdata2$len <- grdata2$FL..cm.
predtimes <- seq(min(grdata2$AgeYears),max(grdata2$AgeYears),by=0.01)
hill_bayes_fitted <- fitted(s1G, 
                            newdata=list(age = predtimes,
                                         AgeClass2=TRUE)) %>% 
  as_tibble()
hill_bayes_pred <- predict(s1G,
                           newdata=list(age = predtimes,
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

