library(brms)
library(tidyverse)


fmodel <- readRDS("Github Code and Data/Data/Growth models/brmS1_Female.rds")
mmodel <- readRDS("Github Code and Data/Data/Growth models/brmS1_Male.rds")

summary(fmodel)                  
summary(mmodel)                  

plot(fmodel)
plot(mmodel)

# convergence of chains OK for both

### make plots

### Make Final Plot
grdata <- read_csv("Github Code and Data/Data/2021 growth modelling data.csv")
grdata$AgeYears <- grdata$AgeMonths/12
grdata$len <- grdata$FL..cm.
grdata2 <- grdata
grdata2 <- grdata2 %>% filter(AgeYears < 1)
grdata3 <- grdata %>% filter(AgeYears >= 1 & Sex != "U")
grdataM <- grdata3 %>% filter(Sex == "M")
grdataF <- grdata3 %>% filter(Sex == "F")

#grdataM <- bind_rows(grdataM, grdata2)
#grdataF <- bind_rows(grdataF, grdata2)

predtimes <- seq(min(grdata2$AgeYears),max(grdata2$AgeYears),by=0.01)
hill_bayes_fitted <- fitted(fmodel, 
                            newdata=list(age = predtimes,
                                         AgeClass2=TRUE)) %>% 
  as_tibble()
hill_bayes_pred <- predict(fmodel,
                           newdata=list(age = predtimes,
                                        AgeClass2=TRUE)) %>%  
  as_tibble()
hill_bayes_ribbonsf <- tibble(
  AgeYears = predtimes,
  parentFraction=hill_bayes_fitted$Estimate,
  Estimate = hill_bayes_fitted$Estimate,
  pred_lower = hill_bayes_pred$Q2.5,
  pred_upper = hill_bayes_pred$Q97.5,
  fitted_lower = hill_bayes_fitted$Q2.5,
  fitted_upper = hill_bayes_fitted$Q97.5)

### now male

hill_bayes_fitted <- fitted(mmodel, 
                            newdata=list(age = predtimes,
                                         AgeClass2=TRUE)) %>% 
  as_tibble()
hill_bayes_pred <- predict(mmodel,
                           newdata=list(age = predtimes,
                                        AgeClass2=TRUE)) %>%  
  as_tibble()
hill_bayes_ribbonsm <- tibble(
  AgeYears = predtimes,
  parentFraction=hill_bayes_fitted$Estimate,
  Estimate = hill_bayes_fitted$Estimate,
  pred_lower = hill_bayes_pred$Q2.5,
  pred_upper = hill_bayes_pred$Q97.5,
  fitted_lower = hill_bayes_fitted$Q2.5,
  fitted_upper = hill_bayes_fitted$Q97.5)

ggplot(data=grdata2, aes(x=AgeYears)) +
  geom_point(aes(y=`len`),size=2, alpha=0.7) +
  geom_point(data=grdataM, aes(x=AgeYears+(0.5*1/12) ,y=`len`),size=2, alpha=0.5, col="blue", inherit.aes = F)+
  geom_point(data=grdataF, aes(x=AgeYears, y=`len`),size=2, alpha=0.3, col="red", inherit.aes = F)+
  #geom_ribbon(data=hill_bayes_ribbons, aes(ymin=pred_lower, ymax=pred_upper), 
  #            alpha=0.2)+#, fill=colourcodes[3]) +
  geom_ribbon(data=hill_bayes_ribbonsf, aes(ymin=fitted_lower, ymax=fitted_upper), 
              alpha=0.5, fill="red")+#, fill=colourcodes[3]) +
  geom_line(data=hill_bayes_ribbonsf, aes(y=Estimate), colour="red", 
            size=1)+
  geom_ribbon(data=hill_bayes_ribbonsm, aes(ymin=fitted_lower, ymax=fitted_upper), 
              alpha=0.5, fill="blue")+#, fill=colourcodes[3]) +
  geom_line(data=hill_bayes_ribbonsm, aes(y=Estimate), colour="blue", 
            size=1)+
  theme_classic() +
  scale_x_continuous(breaks=seq(0,8))+
  theme(axis.title = element_text(face="bold", size=14),
        axis.text = element_text(size=12, colour="black"),
        axis.ticks = element_line(colour="black")) +
  labs(y="Fork Lenth (cm)", x= "Age (years)") +
  geom_text(data=NULL, aes(x=6, y =20, label="Male"), col="blue", size=10, inherit.aes = F)+
  geom_text(data=NULL, aes(x=6, y =30, label="Female"), col="red", size=10, inherit.aes = F)
ggsave("Output/Sex-specific Schnute1.png", dpi=600, width=21, units="cm", height=14.8)
