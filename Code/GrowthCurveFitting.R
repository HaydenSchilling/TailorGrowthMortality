library(dplyr)
#library(nlme)
library(ggplot2)

vbTypical <- function(Linf, K, t0){
  Linf*(1-exp(-K*(ageX-t0)))}

schnute1 <- function(age,alpha,beta,a1,a2,y1,y2)
{
  #function to return the fitted length from a schnute type 1 relationship
  t1 <- 1-exp(-alpha*(age-a1))
  t2 <- 1-exp(-alpha*(a2-a1))
  t3 <- y2^beta - y1^beta
  t4 <- y1^beta + t3*t1/t2
  return (t4^(1/beta))
}

### growth curve
y1 <- 24.38
y2 <- 47.36
alpha <- -0.15 
beta <- 2.56 
ageX <- seq(0,13,0.01)

Sch1 <- ((y1^beta + (y2^beta - y1^beta)*(1-exp(-alpha*(ageX-1)))/(1-exp(-alpha*(4-1))))^(1/beta))
summary(Sch1)

plot(ageX, Sch1, type="l", xlab = "Age (Years)", 
     ylab = "Fork Length (cm)", lwd = 2, axes = T)



# to get other populations
# first AUS
plot(ageX, Sch1, type="l", ylim=c(0,85), xlab = "Age (Years)", ylab = "Fork Length (cm)", lwd = 2)
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
vb_EA <- Linf_EA*(1-exp(-K_EA*(ageX-t0_EA)))

# Von bert for northwest africa
Linf_AF <- 104.326334
K_AF <- 0.214125
t0_AF <- -0.053871
vb_AF <- Linf_AF*(1-exp(-K_AF*(ageX-t0_AF)))

ddat <- data.frame("age" = ageX, "US" = vb_US, "AF" = vb_AF,
                   "WA" = vb_WA, "Brazil" = vb_Bzl, "SA" = vb_SA,
                   "Med" = vb_Med, "Aus" = Sch1)

ddat <- ddat %>% mutate(AF = case_when(age > 9 ~ NA_real_,
                                       T ~ AF),
                        WA = case_when(age > 10 ~ NA_real_,
                                       T ~ WA),
                        Brazil = case_when(age > 8 ~ NA_real_,
                                           T ~ Brazil),
                        SA = case_when(age > 6 ~ NA_real_,
                                       T ~ SA),
                        Med = case_when(age > 6 ~ NA_real_,
                                        T ~ Med),
                        Aus = case_when(age > 7 ~ NA_real_,
                                        T ~ Aus))

ggplot(ddat, aes(x=age)) + theme_bw() + xlab("Age (Years)") + ylab("Fork Length (cm)") + theme_classic() +
  geom_line(aes(y=US, colour="Northwest Atlantic Ocean"),linetype=2, size = 1)+
  geom_line(aes(y=Med, colour = "Mediterranean Sea"),linetype=3, size = 1)+
  geom_line(aes(y=Brazil,colour= "Southwest Atlantic Ocean"),linetype=4, size = 1)+
  geom_line(aes(y=SA, colour="West Indian Ocean"), linetype=5, size = 1)+
  geom_line(aes(y=WA, colour = "East Indian Ocean"),  linetype=6,  size = 1)+
  geom_line(aes(y=AF, colour = "East Atlantic Ocean"),   size = 1)+
  geom_line(aes(y=Aus, colour = "Southwest Pacific Ocean"), size = 1)+
  scale_colour_manual(values = c("Northwest Atlantic Ocean" ="blue","Mediterranean Sea" = "red", "Southwest Atlantic Ocean" = "green", 
                                 "West Indian Ocean"="orange", "East Indian Ocean"="purple", "East Atlantic Ocean"="gray", "Southwest Pacific Ocean"="black"),
                      guide = guide_legend(override.aes = list(linetype = c("dashed","dotted","dotdash","longdash", "twodash", "solid", "solid")))) +
  theme(legend.justification=c(1,0), legend.position=c(1,0), legend.text = element_text(size = 12), legend.title=element_blank(),
        legend.background = element_rect(fill="transparent"), legend.key.width = unit(5, "line")) + ylim(0,82) + scale_x_continuous(breaks=seq(0,13,1)) +
  theme(axis.title.x = element_text(face="bold", colour="black", size = 20),
        axis.text.x  = element_text(colour="black", size = 14), 
        axis.title.y = element_text(face="bold", colour="black", size = 20),
        axis.text.y  = element_text(colour="black", size = 14)) 

ggsave("2022 growth curve comparison.pdf", dpi =600, units="cm", width=21, height = 14.8)
ggsave("2022 growth curve comparison.png", dpi =600, units="cm", width=21, height = 14.8)

# p2 <- ggplot(data.frame(x=c(0, 13)), aes(x)) + theme_bw() + xlab("Age (Years)") + ylab("Fork Length (cm)") + theme_classic() +
#   geom_point(ddat, aes(x=age, y = US), )
#   stat_function(aes(colour = "Northwest Atlantic Ocean"), fun = vbTypical, args = c(Linf_US, K_US, t0_US),linetype=2, xlim = c(0,13), size = 1) +
#   stat_function(aes(colour = "Mediterranean Sea"),fun = vbTypical, args = c(Linf_Med, K_Med, t0_Med),linetype=3, xlim = c(0,6), size = 1) +
#   stat_function(aes(colour = "Southwest Atlantic Ocean"),fun = vbTypical, args = c(Linf_Bzl, K_Bzl, t0_Bzl),linetype=4, xlim = c(0,8), size = 1) +
#   stat_function(aes(colour = "West Indian Ocean"),fun = vbTypical, args = c(Linf_SA, K_SA, t0_SA), xlim = c(0,6), linetype=5, size = 1) +
#   stat_function(aes(colour = "East Indian Ocean"),fun = vbTypical, args = c(Linf_WA, K_WA, t0_WA), xlim = c(0,10), linetype=6,  size = 1) +
#   stat_function(aes(colour = "East Atlantic Ocean"),fun = vbTypical, args = c(Linf_AF, K_AF, t0_AF), xlim = c(0,9),  size = 1) +
#   stat_function(aes(colour = "Southwest Pacific Ocean"),fun = schnute1, args = c(alpha,beta, y1, y2,a1=1, a2=4), size = 1, xlim = c(0,7)) +
#   scale_colour_manual(values = c("Northwest Atlantic Ocean" ="blue","Mediterranean Sea" = "red", "Southwest Atlantic Ocean" = "green", 
#                                  "West Indian Ocean"="orange", "East Indian Ocean"="purple", "East Atlantic Ocean"="gray", "Southwest Pacific Ocean"="black"),
#                       guide = guide_legend(override.aes = list(linetype = c("solid","twodash","dotted","dashed", "dotdash", "solid", "longdash")))) +
#   theme(legend.justification=c(1,0), legend.position=c(1,0), legend.text = element_text(size = 12), legend.title=element_blank(),
#         legend.background = element_rect(fill="transparent"), legend.key.width = unit(5, "line")) + ylim(0,82) + scale_x_continuous(breaks=seq(0,13,1)) +
#   theme(axis.title.x = element_text(face="bold", colour="black", size = 20),
#         axis.text.x  = element_text(colour="black", size = 14), 
#         axis.title.y = element_text(face="bold", colour="black", size = 20),
#         axis.text.y  = element_text(colour="black", size = 14)) #+
#   #theme_classic()
# p2
# 
# # Male v Female growth curves
# p6 <- ggplot(dat, aes(x = Years, y = FL..cm.)) + geom_point(size=1) + theme_bw() + 
#   xlab("Age (Years)") + ylab("Fork Length (cm)") + scale_x_continuous(breaks=seq(0,8,1)) +
#   stat_function(aes(colour = "Females"),fun = schnute3, args = c(betaF, y1F, y2F, a2=3), size = 0.9) +
#   stat_function(aes(colour = "Males"),fun = schnute3, args = c(betaM, y1M, y2M, a2=3), size = 0.9) +
#   stat_function(aes(colour = "Combined sexes"),fun = schnute3, args = c(beta, y1, y2, a2=3), size = 0.9) +
#   # stat_function(aes(colour = "Von Bertalanffy"),fun = vbTypical, args = c(Linf_EA, K_EA, t0_EA), size = 0.5, linetype=2) +
#   scale_colour_manual(values = c("Females" ="red","Males" = "blue", "Combined sexes" = "black", "von Bertalanffy"="blue"),
#                       guide = guide_legend(override.aes = list(linetype = c("solid","solid", "solid")))) +
#   theme(legend.justification=c(1,0), legend.position=c(1,0), legend.text = element_text(size = 8), 
#         legend.title=element_blank(), legend.background = element_rect(fill="transparent"), legend.key.width = unit(2.5, "line"),
#         panel.grid.minor = element_blank(), panel.grid.major = element_blank())
# p6
# 
# 
# p4 <- multiplot(p1,p2,p3, layout= matrix(c(1,1,2,3), nrow=2, byrow=TRUE))
# p5 <- multiplot(p1,p2)



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

ggsave("2022 growth by difference.pdf", dpi=600, units="cm", width=21, height=14.8)
ggsave("2022 growth by difference.png", dpi=600, units="cm", width=21, height=14.8)

c("solid","twodash","dotted","dashed", "dotdash", "solid", "longdash")

# ylab(expression(paste("Growth in Previous Year  (cm)", ~year^{-1},")")))