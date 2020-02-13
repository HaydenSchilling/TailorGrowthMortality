library(dplyr)
library(nlme)
library(ggplot2)


dat <-  read.csv("Final Adjusted Age Data for growth modelling.csv")
dat = subset(dat, is.na(Months)==F)#using Months as response here

# K is a growth rate parameter, t0 is a theoretical time at size 0, Linf is a theoretical asymptotical max length
vbTypical <- function(AgeClass, Linf,  K, t0){Linf*(1-exp(-K*(AgeClass-t0)))}

#Von Bertalanffy on a log log scale (i.e. log(length) grows with log(months) according to the von Bert growth equation)
vbTypicallog <- function(AgeClass, Linf,  K, t0){Linf*(1-exp(-K*(log(AgeClass)-t0)))}

#binary age class- experiment with applying weights that weight the juveniles(<1 yr) differently as they have unequal variances
dat$AgeClass2 = dat$AgeClass >=1
weights = dat %>% group_by(AgeClass2) %>%summarise(weights1 = var(FL..cm.))
dat = left_join(dat, weights) 




### fit typical seems to fit much better on a log log scale.
fitTypical = nls(FL..cm. ~ vbTypical(Months, Linf,  K, t0), data=dat,start=c(Linf = 80, K = 0.3, t0 = 0), control=nls.control(maxiter=10000))
fitTypical
summary(fitTypical)
plot(fitTypical) # maybe need nlme

#on a log log scale
fitTypicallog = nls(log(FL..cm.) ~ vbTypicallog(Months, Linf,  K, t0), data=dat,start=c(Linf = 80, K = 0.3, t0 = 0), control=nls.control(maxiter=10000))
plot(fitTypicallog)

#on a log-log scale with a weighted regression. Less fan shape in the residuals....
fitTypicallog2 = nls(log(FL..cm.) ~ vbTypicallog(Months, Linf,  K, t0), data=dat,start=c(Linf = 80, K = 0.3, t0 = 0), control=nls.control(maxiter=10000), weights=1/weights1)
plot(fitTypicallog2)

### create log variables and fit schnute3 on a log log scale
dat$log.leng  =log(dat$FL..cm.)
dat$log.age = log(dat$Months)


schnute3 <- function(age,beta,a1=1,a2=4,y1,y2)
{
  #function to return the fitted length from a schnute type 3 relationship
  t1 <- (age-a1)/(a2-a1)
  t3 <- y2^beta - y1^beta
  t4 <- y1^beta + t3*t1
  return (t4^(1/beta))
}


sch3<- nls(FL..cm. ~ schnute3(AgeClass,beta,1,3,y1,y2), start = list(beta=1,y1=35,y2=50), data=dat,control=list(maxiter=10000,nlsTols=10000,minScale=10^-6))
plot(sch3)
summary(sch3)

sch3log<- nls(log.leng ~ schnute3(log.age,beta,0.5,1,y1,y2), start = list(beta=1.1,y1=3,y2=5), data=dat,control=list(maxiter=10000,nlsTols=10000,minScale=10^-6))
plot(sch3log)
#looks like it fits better for the young fish

#hacky plot on log log scale. Red line shows fit typical, blue shows schnute 3. 
plot(log(dat$Months), log(dat$FL..cm.))
lines(log(sort(dat$Months)),sort(predict(fitTypicallog) ), col="red")
lines(log(sort(dat$Months)),sort(predict(sch3log) ), col="blue")
lines(log(sort(dat$Months)),log(sort(predict(fitTypical) )), col="red", lty=2)
lines(log(sort(dat$Months)),log(sort(predict(sch3) )), col="blue", lty=2)
lines(log(sort(dat$Months)[-1]),sort(log(predict(fitsch3.glns.) )), col="green", lty=2)


#hacky plot on original scale. 
plot(dat$Years, dat$FL..cm.)
lines(sort(dat$Months),exp(sort(predict(fitTypicallog))) , col="red")
lines(sort(dat$Months),exp(sort(predict(sch3log))) , col="blue")
lines(sort(dat$Months),sort(predict(fitTypical) ), col="red", lty=2)
lines(sort(dat$Months),sort(predict(sch3) ), col="blue", lty=2)

lines(sort(dat$Years)[-1],sort(predict(fitsch3.glns) ), col="green", lty=2)


#################################
#GNLS
##############################



## From what the internet tells me, there is a bug in gnls. gnls should be able to fit the same model as nls, yet even before you add the weights it throws errors. 
## This is unfortunate! If you fiddle enough with the control parameters it works some of the time. 
###i.e. this should work (no difference to the nls model)
fitTypical2 = gnls(FL..cm. ~ vbTypical(Months, Linf,  K, t0), data=dat,start=c(Linf = 80, K = 0.3, t0 = 0), control=list(maxiter=1000000,nlsTols=100, minScale=10^-6))
### with the aim to add the weights function in 
fitTypical2 = gnls(FL..cm. ~ vbTypical(Years, Linf,  K, t0), data=dat,start=c(Linf = 120, K = 0.1, t0 = -2), weights = varIdent(form=~1|AgeClass2), control=list(maxiter=1000000000,nlsTols=1000, minScale=10^-6))
summary (fitTypical2)

fitTypical2 = gnls(log(FL..cm.) ~ vbTypicallog(Months, Linf,  K, t0), data=dat,start=c(Linf = 80, K = 0.3, t0 = 0), weights = varIdent(form=~1|AgeClass2), control=list(maxiter=1000000,nlsTols=100, minScale=10^-6))
plot(fitTypical2)

dat$Years <- dat$Months/12

# write.csv(dat, "sepsex.csv", row.names =F)

# Add 0.5 to age
for (i in 1:nrow(dat)) {
  if (dat$AgeClass[i] >= 1) {
    dat$AgeClass[i] <- dat$AgeClass[i] + 0.5
    #} else if ((mydata$AgeClass[i] == 0)) {
    #  mydata$AgeClass[i] <- mydata$AgeClass[i] + 0.5
  } else {
    dat$AgeClass[i] <- dat$AgeClass[i]
  }
}

dat <- subset(dat, AgeClass <11)
dat$Age <- dat$AgeClass*12

#This is the sch3 model on 
fitsch3log.glns = gnls(log.leng ~ schnute3(log.age,beta,0.5,1,y1,y2), start = list(beta=1.1,y1=3,y2=5),data=dat,weights = varIdent(form=~1|AgeClass2), control=list(maxiter=1000000,nlsTols=100, minScale=10^-6))

# Use Months to get fit (aka ageclass*12)
fitsch3.glns = gnls(FL..cm. ~ schnute3(Age,beta,12,48,y1,y2), 
                    #start = list(beta=2.5,y1=25,y2=42), 
                    start = coef(sch3),
                      data=dat,
                    weights = varIdent(form=~1|AgeClass2), 
                    control=list(maxiter=1000000,nlsTols=10000, minScale=10^-10))
summary(fitsch3.glns)
plot(fitsch3.glns)
plot(fitsch3log.glns)
intervals(fitsch3.glns) # 95% CI for parameters



# Try separate sexes
dat <- read.csv("sepsex.csv", header=T)

# Add 0.5 to age
for (i in 1:nrow(dat)) {
  if (dat$AgeClass[i] >= 1) {
    dat$AgeClass[i] <- dat$AgeClass[i] + 0.5
    #} else if ((mydata$AgeClass[i] == 0)) {
    #  mydata$AgeClass[i] <- mydata$AgeClass[i] + 0.5
  } else {
    dat$AgeClass[i] <- dat$AgeClass[i]
  }
}

dat <- subset(dat, AgeClass <11)

dat$Age <- dat$AgeClass*12
males <- dat[ !(dat$Sex == 'F'),]
females <- dat[!(dat$Sex == 'M'),]

fitsch3.glns.M = gnls(FL..cm. ~ schnute3(Age,beta,12,48,y1,y2), start = list(beta=2.3,y1=25.56,y2=46.34), data=males,weights = varIdent(form=~1|AgeClass2), control=list(maxiter=1000000,nlsTols=10000, minScale=10^-10))
fitsch3.glns.F = gnls(FL..cm. ~ schnute3(Age,beta,12,48,y1,y2), start = list(beta=2.49,y1=25.56,y2=46.34), data=females,weights = varIdent(form=~1|AgeClass2), control=list(maxiter=1000000,nlsTols=10000, minScale=10^-10))

fitsch3.glns.M = nls(FL..cm. ~ schnute3(AgeClass,beta,1,3,y1,y2), start = list(beta=2.5,y1=25.56,y2=40), data=males, control=list(maxiter=1000000,nlsTols=10000, minScale=10^-10))
fitsch3.glns.F = nlsLM(FL..cm. ~ schnute3(AgeClass,beta,1,3,y1,y2), start = list(beta=1,y1=25.4,y2=40), data=females, control=list(maxiter=1000000,nlsTols=10000, minScale=10^-10))

summary(fitsch3.glns.M)
summary(fitsch3.glns.F)

plot(females$AgeClass, females$FL..cm.)

intervals(fitsch3.glns.M)#  95% CI for parameters
intervals(fitsch3.glns.F)


### Likehood Ratio Test
# Not working yet
#extract RSS for each model
Combined.RSS<- sum(resid(fitsch3.glns)^2) #RSS from combined model
separate_sex.RSS<-sum(resid(fitsch3.glns.M)^2)+sum(resid(fitsch3.glns.F)^2) # RSS from separate models which has been combined

#### Calculate LRT for each model
#VB
sch3.x2<--nrow(dat)*log(separate_sex.RSS/Combined.RSS)#calculate chi sq
df<-3 #number of degrees of freedom
sch3.p.val<-1-pchisq(sch3.x2,df)#calculated p value from df and x2

#results of individual candidate models
sch3.LRT.Result<-sch3.p.val<0.05 # true/false statement. True means significant difference
sch3.LRT.Result
sch3.p.val





###

RMSE = function(pred, obs){
  sqrt(mean((pred-obs)^2, na.rm=T))
}
##how does fitsch3.glns compare to fitsch3, and fitsch3log ? 
#a) on the original scale, b)on the loglog scale

#rmse's sch3 model
RMSE(predict(sch3), exp(dat$log.leng))
RMSE(log(predict(sch3)), dat$log.leng)

#rmse's sch3 model log
RMSE(predict(sch3log), dat$log.leng)
RMSE(exp(predict(sch3log)), exp(dat$log.leng))#worse on the original scale

#rmse's sch3 model weighted
#maybe the best
RMSE(predict(fitsch3.glns), exp(dat$log.leng))
RMSE(log(predict(fitsch3.glns)), dat$log.leng)

RMSE(exp(predict(fitsch3log.glns)), exp(dat$log.leng))
RMSE(predict(fitsch3log.glns), dat$log.leng)

# a visual of what's going on?
#hacky plot on log log scale. Red line shows fit typical, blue shows schnute 3. 
plot(log(dat$Months), log(dat$FL..cm.))
lines(log(sort(dat$Months)),sort(predict(fitTypicallog) ), col="red")
lines(log(sort(dat$Months)),sort(predict(sch3log) ), col="blue")
lines(log(sort(dat$Months)),log(sort(predict(fitTypical) )), col="red", lty=2)
lines(log(sort(dat$Months)),log(sort(predict(sch3) )), col="blue", lty=2)
lines(log(sort(dat$Months)[-1]),sort(log(predict(fitsch3.glns) )), col="green", lty=2) # green is what you get by giving the age classes two variances


#hacky plot on original scale. 
plot(dat$Months, dat$FL..cm.)
lines(sort(dat$Months),exp(sort(predict(fitTypicallog))) , col="red")
lines(sort(dat$Months),exp(sort(predict(sch3log))) , col="blue")
lines(sort(dat$Months),sort(predict(fitTypical) ), col="red", lty=2)
lines(sort(dat$Months),sort(predict(sch3) ), col="blue", lty=2)

lines(sort(dat$Months)[-1],sort(predict(fitsch3.glns) ), col="green", lty=2)# green is what you get by giving the age classes two variances



layout(matrix(c(1,1,1,1,2,3), 3, 2, byrow = TRUE))

### To do a nice plot
ageX <-  seq(0,8, by = 0.01)
# Equation to plot
Sch3_Pars <- coef(fitsch3.glns)
Sch3_ParsF <- coef(fitsch3.glns.F)
Sch3_ParsM <- coef(fitsch3.glns.M)

beta = Sch3_Pars[1]
y1 = Sch3_Pars[2]
y2 = Sch3_Pars[3]

betaF = Sch3_ParsF[1]
y1F = Sch3_ParsF[2]
y2F = Sch3_ParsF[3]

betaM = Sch3_ParsM[1]
y1M = Sch3_ParsM[2]
y2M = Sch3_ParsM[3]

beta = 2.49
y1 = 25.46
y2 = 46.34


est_Sch3 <- ((y1^beta + (y2^beta - y1^beta)*((ageX-1)/(4-1)))^(1/beta))
est_Sch3F <- ((y1F^betaF + (y2F^betaF - y1F^betaF)*((ageX-1)/(3-1)))^(1/betaF))
est_Sch3M <- ((y1M^betaM + (y2M^betaM - y1M^betaM)*((ageX-1)/(3-1)))^(1/betaM))

plot(ageX, est_Sch3, type="l", ylim=c(0,85), xlab = "Age (Years)", 
     ylab = "Fork Length (cm)", lwd = 2, axes = FALSE)
axis(side = 1, at = seq(0,8,1))
axis(side = 2, at = seq(0,80,20))
points(dat$AgeClass, dat$FL..cm., pch = 16)
lines(ageX, est_Sch3M, col = "blue", lty = 2, lwd = 2)
lines(ageX, est_Sch3F, col = "red", lty = 2, lwd = 2)

dev.copy(pdf,"Output/Growth curve.pdf", width=8, height=6)
dev.off()

# to get other populations
# first AUS
plot(ageX, est_Sch3, type="l", ylim=c(0,85), xlab = "Age (Years)", ylab = "Fork Length (cm)", lwd = 2)
# US
Linf_US <- 81.53
t0_US <--0.301
K_US <- 0.311
vb_US <- Linf_US*(1-exp(-K_US*(ageX-t0_US)))
lines(ageX, vb_US, col = "green", lwd = 2)

# Mediter... (Cengiz 2013)
Linf_Med <- 88.3
t0_Med <--1.43
K_Med <- 0.15
vb_Med <- Linf_Med*(1-exp(-K_Med*(ageX-t0_Med)))
lines(ageX, vb_Med, col = "orange", lwd = 2)

# South Africa (Gosvener)
Linf_SA <- 124.7
t0_SA <- -2.09
K_SA <- 0.094
vb_SA <- Linf_SA*(1-exp(-K_SA*(ageX-t0_SA)))
lines(ageX, vb_SA, col = "brown", lwd = 2)

# Brazil from Manuel Haimovici and Luiz Carlos Krug (1996) MFR
Linf_Bzl <- 66.2
t0_Bzl <- 0.321
K_Bzl <- 0.387
vb_Bzl <- Linf_Bzl*(1-exp(-K_Bzl*(ageX-t0_Bzl)))
lines(ageX, vb_Bzl, col = "pink", lwd = 2)

# WA from Smith et al 2013 report
Linf_WA <- 59.2
t0_WA <- -0.096
K_WA <- 0.464
vb_WA <- Linf_WA*(1-exp(-K_WA*(ageX-t0_WA)))
lines(ageX, vb_WA, col = "pink", lwd = 2)

# Von bert for east Aus
Linf_EA <- 104.36
K_EA <- 0.09966877
t0_EA <- -1.97950201

# Von bert for northwest africa
Linf_AF <- 104.326334
K_AF <- 0.214125
t0_AF <- -0.053871
  
dat <- subset(dat, Years <10)
dat <- subset(dat, AgeClass <10)

# Daily growth
daily_data <- grdata[ which(grdata$AgeClass < 1),]
daily_data$age <- daily_data$age*365 # to get days not years
plot(len ~ age, data = daily_data, ylab = "Length (cm)", xlab = "Age (days)",
     pch = 16, xlim=c(0,250), ylim=c(0,20))
fit_daily <- lm(len ~ age, data = daily_data)
summary(fit_daily)
abline(fit_daily)

# Just combined sex and von bert
p1 <- ggplot(dat, aes(x = AgeClass, y = FL..cm.)) + geom_point(size=1, alpha = 0.5) + theme_classic() + 
  xlab("Age (Years)") + ylab("Fork Length (cm)") + scale_x_continuous(breaks=seq(0,7,1)) +
  # stat_function(aes(colour = "Females"),fun = schnute3, args = c(betaF, y1F, y2F, a2=3), size = 0.5) +
  # stat_function(aes(colour = "Males"),fun = schnute3, args = c(betaM, y1M, y2M, a2=3), size = 0.5) +
  stat_function(aes(colour = "Schnute 3"),fun = schnute3, args = c(beta, y1, y2,a1=1, a2=4), size = 0.5) +
  stat_function(aes(colour = "von Bertalanffy"),fun = vbTypical, args = c(Linf_EA, K_EA, t0_EA), size = 0.5, linetype="dashed") +
  scale_colour_manual(values = c("Schnute 3" ="black","von Bertalanffy" = "blue"),
                      guide = guide_legend(override.aes = list(linetype = c("solid","twodash")))) + 
  theme(legend.justification=c(1,0), legend.position=c(1,0.05),legend.key.width = unit(2.5, "line"),
        axis.title.x = element_text(face="bold", colour="black", size = 20),
        axis.text.x  = element_text(colour="black", size = 16), 
        axis.title.y = element_text(face="bold", colour="black", size = 20),
        axis.text.y  = element_text(colour="black", size = 16),
        legend.title = element_text(size=16, face="bold"), 
        legend.text = element_text(size = 14)) + 
  scale_y_continuous(breaks=seq(0,80,20)) +
  guides(colour=guide_legend(title="Growth Model"))
p1


p2 <- ggplot(data.frame(x=c(0, 13)), aes(x)) + theme_bw() + xlab("Age (Years)") + ylab("Fork Length (cm)") + theme_classic() +
  stat_function(aes(colour = "Northwest Atlantic Ocean"), fun = vbTypical, args = c(Linf_US, K_US, t0_US),linetype=2, xlim = c(0,13), size = 1) +
  stat_function(aes(colour = "Mediterranean Sea"),fun = vbTypical, args = c(Linf_Med, K_Med, t0_Med),linetype=3, xlim = c(0,6), size = 1) +
  stat_function(aes(colour = "Southwest Atlantic Ocean"),fun = vbTypical, args = c(Linf_Bzl, K_Bzl, t0_Bzl),linetype=4, xlim = c(0,8), size = 1) +
  stat_function(aes(colour = "West Indian Ocean"),fun = vbTypical, args = c(Linf_SA, K_SA, t0_SA), xlim = c(0,6), linetype=5, size = 1) +
  stat_function(aes(colour = "East Indian Ocean"),fun = vbTypical, args = c(Linf_WA, K_WA, t0_WA), xlim = c(0,10), linetype=6,  size = 1) +
  stat_function(aes(colour = "East Atlantic Ocean"),fun = vbTypical, args = c(Linf_AF, K_AF, t0_AF), xlim = c(0,9),  size = 1) +
  stat_function(aes(colour = "Southwest Pacific Ocean"),fun = schnute3, args = c(beta, y1, y2,a1=1, a2=4), size = 1, xlim = c(0,7)) +
  scale_colour_manual(values = c("Northwest Atlantic Ocean" ="blue","Mediterranean Sea" = "red", "Southwest Atlantic Ocean" = "green", 
                                 "West Indian Ocean"="orange", "East Indian Ocean"="purple", "East Atlantic Ocean"="gray", "Southwest Pacific Ocean"="black"),
                      guide = guide_legend(override.aes = list(linetype = c("solid","twodash","dotted","dashed", "dotdash", "solid", "longdash")))) +
  theme(legend.justification=c(1,0), legend.position=c(1,0), legend.text = element_text(size = 12), legend.title=element_blank(),
        legend.background = element_rect(fill="transparent"), legend.key.width = unit(5, "line")) + ylim(0,82) + scale_x_continuous(breaks=seq(0,13,1)) +
  theme(axis.title.x = element_text(face="bold", colour="black", size = 20),
        axis.text.x  = element_text(colour="black", size = 14), 
        axis.title.y = element_text(face="bold", colour="black", size = 20),
        axis.text.y  = element_text(colour="black", size = 14)) #+
  #theme_classic()
p2

p3 <- ggplot(daily_data, aes(x = age, y = len)) + geom_point() + theme_classic() +
  geom_smooth(method = "lm", se = FALSE, col = "black") + xlim(0,250) +
  xlab("Age (days)") + ylab("Fork Length (cm)") +
  theme(axis.title.x = element_text(face="bold", colour="black", size = 20),
        axis.text.x  = element_text(colour="black", size = 14), 
        axis.title.y = element_text(face="bold", colour="black", size = 20),
        axis.text.y  = element_text(colour="black", size = 14))
p3

# Male v Female growth curves
p6 <- ggplot(dat, aes(x = Years, y = FL..cm.)) + geom_point(size=1) + theme_bw() + 
  xlab("Age (Years)") + ylab("Fork Length (cm)") + scale_x_continuous(breaks=seq(0,8,1)) +
  stat_function(aes(colour = "Females"),fun = schnute3, args = c(betaF, y1F, y2F, a2=3), size = 0.9) +
  stat_function(aes(colour = "Males"),fun = schnute3, args = c(betaM, y1M, y2M, a2=3), size = 0.9) +
  stat_function(aes(colour = "Combined sexes"),fun = schnute3, args = c(beta, y1, y2, a2=3), size = 0.9) +
  # stat_function(aes(colour = "Von Bertalanffy"),fun = vbTypical, args = c(Linf_EA, K_EA, t0_EA), size = 0.5, linetype=2) +
  scale_colour_manual(values = c("Females" ="red","Males" = "blue", "Combined sexes" = "black", "von Bertalanffy"="blue"),
                      guide = guide_legend(override.aes = list(linetype = c("solid","solid", "solid")))) +
  theme(legend.justification=c(1,0), legend.position=c(1,0), legend.text = element_text(size = 8), 
        legend.title=element_blank(), legend.background = element_rect(fill="transparent"), legend.key.width = unit(2.5, "line"),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank())
p6


p4 <- multiplot(p1,p2,p3, layout= matrix(c(1,1,2,3), nrow=2, byrow=TRUE))
p5 <- multiplot(p1,p2)


tiff("poster pics/growth subplots.tiff", height = 250, width = 290, units = 'mm', res = 300)
multiplot(p1, p2)
dev.off()

# Size at age
#vbTypical <- function(AgeClass, Linf,  K, t0){Linf*(1-exp(-K*(AgeClass-t0)))}
US <- vbTypical(c(1,2,3,4,5),Linf_US, K_US, t0_US)
US

Med <- vbTypical(c(1,2,3,4,5),Linf_Med, K_Med, t0_Med)
Med

Brazil <- vbTypical(c(1,2,3,4,5),Linf_Bzl, K_Bzl, t0_Bzl)
Brazil

SA <- vbTypical(c(1,2,3,4,5),Linf_SA, K_SA, t0_SA)
SA

WA <- vbTypical(c(1,2,3,4,5),Linf_WA, K_WA, t0_WA)
WA

Africa <- vbTypical(c(1,2,3,4,5),Linf_AF, K_AF, t0_AF)
Africa


# size at age
est_Sch3_1 <- ((y1^beta + (y2^beta - y1^beta)*((1-1)/(4-1)))^(1/beta))
est_Sch3_1
est_Sch3_2 <- ((y1^beta + (y2^beta - y1^beta)*((2-1)/(4-1)))^(1/beta))
est_Sch3_2
est_Sch3_3 <- ((y1^beta + (y2^beta - y1^beta)*((3-1)/(4-1)))^(1/beta))
est_Sch3_3
est_Sch3_4 <- ((y1^beta + (y2^beta - y1^beta)*((4-1)/(4-1)))^(1/beta))
est_Sch3_4
est_Sch3_5 <- ((y1^beta + (y2^beta - y1^beta)*((5-1)/(4-1)))^(1/beta))
est_Sch3_5


#### Refit the Africa growth rate using size at age data
africa_data <- read.csv("Tassergal age length data.csv", header = T)
vbTypical <- function(AgeClass, Linf,  K, t0){Linf*(1-exp(-K*(AgeClass-t0)))}
fitTypical = nls(Length ~ vbTypical(Age, Linf,  K, t0), data=africa_data,start=c(Linf = 80, K = 0.3, t0 = 0), control=nls.control(maxiter=10000))
summary(fitTypical)

fitPlot(fitTypical,xlab="Age",ylab="Fork Length (cm)",main="")

############################# 
# Growth rate by difference plotting
############################
mydata <- read.csv("Growth Rates by difference data.csv", header = T)


p_grow <- ggplot(mydata, aes(Age, Growth_Rate, col = Population, linetype= Population)) + geom_line(size = 1.5) +
  theme_classic() + ylab("Growth in Previous Year (cm)") +
  theme(axis.title.x = element_text(face="bold", colour="black", size = 20),
        axis.text.x  = element_text(colour="black", size = 14), 
        axis.title.y = element_text(face="bold", colour="black", size = 20),
        axis.text.y  = element_text(colour="black", size = 14),
        legend.justification=c(1,0), legend.position=c(0.95,0.55), legend.text = element_text(size = 12), legend.title=element_blank(),
        legend.background = element_rect(fill="transparent"), legend.key.width = unit(5, "line")) +
  scale_colour_manual(values = c("gray", "purple", "red", "blue", "green", "black", "orange")) +
  scale_linetype_manual(values= c("solid","twodash","dotted","dashed", "dotdash", "solid", "longdash"))
                                 
p_grow

c("solid","twodash","dotted","dashed", "dotdash", "solid", "longdash")

# ylab(expression(paste("Growth in Previous Year  (cm)", ~year^{-1},")")))