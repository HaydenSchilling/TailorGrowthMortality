library(FSA)
library(nnet)
library(dplyr)
library(ggplot2)

mydata <- read.csv("Final Adjusted Age Data for growth modelling.csv", header = T)
alldata <- mydata

alldata <- alldata %>% mutate(lcat10=lencat(FL..cm.,w=2)) # bin size of 2 (w=2)
alldata <- alldata %>% mutate(AgeClassWHOLE=lencat(AgeClass,w=1))
# alldata$lcat10 <- alldata$lcat10 + 1
head(alldata)
table(alldata$lcat10,alldata$AgeClassWHOLE)
alk.freq <- xtabs(~lcat10+AgeClassWHOLE,data=alldata)
rowSums(alk.freq)
alk <- prop.table(alk.freq,margin=1)
len.n <- xtabs(~lcat10,data=alldata)
# alk.round <- round(alk,3)
# alk.round

# Smoothed age length key, see Ogle 2017 chapter
cc.mlr <- multinom(AgeClassWHOLE~lcat10,data=alldata,maxit=500) # library(nnet)

lens <- seq(20,80,1)
alk.sm <- predict(cc.mlr,data.frame(lcat10=lens),type="probs")
row.names(alk.sm) <- lens # for clarity

alkPlot(alk.sm,type="area",pal="gray",showLegend=TRUE,
        leg.cex=0.7,xlab="Total Length (mm)")

alkPlot(alk,type="bubble",xlab="Total Length (mm)")

write.csv(alk, "age length proportions 2cm.csv")


tag_dat <- read.csv("Historical Tagged FL Clean.csv", header = T)
tag_dat$Age <- as.numeric(rep(0, nrow(tag_dat)))
tag_dat <- subset(tag_dat, Len >= 20)

full_ages <- alkIndivAge(alk, Age~Len, data=tag_dat)
hist(full_ages$Age)

full_ages.sm <- alkIndivAge(alk.sm, Age~Len, data=tag_dat)
hist(full_ages.sm$Age)

full_ages$Age <- as.factor(full_ages$Age)
p1 <- ggplot(full_ages, aes(x = Len, fill = Age)) + geom_histogram(binwidth = 1) + theme_bw()
p1

full_ages.sm$Age <- as.factor(full_ages.sm$Age)

p2 <- ggplot(full_ages.sm, aes(x = Len, fill = Age)) + geom_histogram(binwidth = 1) + theme_bw()
p2

table(full_ages$Age)
table(full_ages.sm$Age)

alkAgeDist(alk,lenA.n=rowSums(alk.freq),len.n=len.n)

com_dat <- read.csv("NSW commercial length distribution.csv", header = T)
com_dat_long <- rep(com_dat$Length_Class, com_dat$nFish)
com_dat_long <- as.data.frame(com_dat_long)
com_dat_long$Age <- as.numeric(rep(0, nrow(com_dat_long)))
names(com_dat_long)[1] <- "Length"

com_full_ages <- alkIndivAge(alk, Age~Length, data=com_dat_long)
hist(com_full_ages$Age, xlab = "Age", main = "NSW Commercial Tailor Age Composition")
table(com_full_ages$Age)

com_full_ages$Age <- as.factor(com_full_ages$Age)
p3 <- ggplot(com_full_ages, aes(x = Length, fill = Age)) + geom_histogram(binwidth = 1) + theme_bw()+
  ggtitle("NSW Commercial")
p3

table(com_full_ages$Age)

nsw_combined <- read.csv("NSW Lengths weighted and Combined.csv", header = T)
names(nsw_combined)[1] <- "Length"
nsw_combined$Age <- as.numeric(rep(0, nrow(nsw_combined)))
nsw_combined_ages <- alkIndivAge(alk, Age~Length, data=nsw_combined)
hist(nsw_combined_ages$Age, xlab = "Age", main = "NSW Overall Age Composition")
table(nsw_combined_ages$Age)

nsw_combined_ages$Age <- as.factor(nsw_combined_ages$Age)
p4 <- ggplot(nsw_combined_ages, aes(x = Length, fill = Age)) + geom_histogram(binwidth = 1) + theme_bw() +
  ggtitle("NSW Overall")
p4

### NSW Recreational
nsw_rec <- read.csv("NSW Current Rec Lengths and Ages.csv", header = T)
names(nsw_rec)[2] <- "Length"
nsw_rec$Age <- as.numeric(rep(0, nrow(nsw_rec)))
nsw_rec_ages <- alkIndivAge(alk, Age~Length, data=nsw_rec)
hist(nsw_rec_ages$Age, xlab = "Age", main = "NSW Recreational Age Composition")
table(nsw_rec_ages$Age)

nsw_rec_ages$Age <- as.factor(nsw_rec_ages$Age)
p4 <- ggplot(nsw_rec_ages, aes(x = Length, fill = Age)) + geom_histogram(binwidth = 1) + theme_bw() +
  ggtitle("NSW Recreational")
p4



#### FINAL EAST AUSTralia
east_aus_combined <- read.csv("East Australia Lengths weighted and Combined 2_cm.csv", header = T)
names(east_aus_combined)[1] <- "Length"
east_aus_combined$Age <- as.numeric(rep(0, nrow(east_aus_combined)))
east_aus_combined_ages <- alkIndivAge(alk, Age~Length, data=east_aus_combined)
hist(east_aus_combined_ages$Age)

east_aus_combined_ages$Age <- as.factor(east_aus_combined_ages$Age)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colourblind palette

pALL <- ggplot(east_aus_combined_ages, aes(x = Length, fill = Age,y = (..count..)/sum(..count..))) + geom_histogram(binwidth = 2) + 
  theme_classic() + scale_fill_manual(values=cbPalette) +
  scale_x_continuous(breaks=seq(20,80,5)) + xlab("Fork Length (cm)") + ylab("Proportion") +
  theme(axis.title.x = element_text(face="bold", colour="black", size = 20),
        axis.text.x  = element_text(colour="black", size = 16), 
        axis.title.y = element_text(face="bold", colour="black", size = 20),
        axis.text.y  = element_text(colour="black", size = 16)) +
  theme(legend.justification=c(0.9,0.9), legend.position=c(0.97,0.9)) +
  theme(legend.title = element_text(size=16, face="bold"), legend.text = element_text(size = 14))
pALL


table(east_aus_combined_ages$Age)

# Historical NSW Data
hist_nsw <- read.csv("NSW Historical Length Data.csv", header = T)
str(hist_nsw)
hist(hist_nsw$Length)
hist_nsw$Age <- as.numeric(rep(0, nrow(hist_nsw)))
hist_nsw <- alkIndivAge(alk, Age~Length, data=hist_nsw)

hist_nsw$Age <- as.factor(hist_nsw$Age)
p5 <- ggplot(hist_nsw, aes(x = Length, fill = Age)) + geom_histogram(binwidth = 1) + theme_bw()
p5

table(hist_nsw$Age)
