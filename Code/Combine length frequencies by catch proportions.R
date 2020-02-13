library(dplyr)
library(ggplot2)

# rec_data <- read.csv("NSW Current Rec Lengths and Ages.csv", header = T)
# rec_data$FL <- round(rec_data$FL)
# rec_data_summary <- as.data.frame(table(rec_data$FL))
# names(rec_data_summary)[1] <- "FL"
# rec_data_summary$FL <- as.numeric(as.character(rec_data_summary$FL))
# rec_data_summary$Rec_Percentage <- rec_data_summary$Freq/sum(rec_data_summary$Freq)
# sum(rec_data_summary$Rec_Percentage)

# write.csv(rec_data_summary, "NSW Rec data 2 cm binned.csv", row.names = F)
rec_data_summary <- read.csv("NSW Rec data 2 cm binned.csv", header = T)

com_data <- read.csv("NSW commercial length distribution 2cm.csv", header = T)
com_data$Com_Percentage <- com_data$nFish/sum(com_data$nFish) 
sum(com_data$Com_Percentage)
names(com_data)[1] <- "FL"

# to combine into one dataset for all NSW
# using 71.68t commercial catch and 107t recreational catch

all_nsw <- full_join(rec_data_summary, com_data, by = "FL")
all_nsw[is.na(all_nsw)] <- 0
names(all_nsw)[2] <- "NSW_Rec_Percentage"
names(all_nsw)[5] <- "NSW_Com_Percentage"
sum(all_nsw$NSW_Com_Percentage)
sum(all_nsw$NSW_Rec_Percentage)

QLD_Rec <- read.csv("QLD Rec LF.csv", header = T)
names(QLD_Rec)[1] <- "FL"

all_data1 <- full_join(all_nsw, QLD_Rec, by = "FL")

QLD_Net <- read.csv("QLD Net LF.csv", header = T)
names(QLD_Net)[1] <- "FL"

all_data2 <- full_join(all_data1, QLD_Net, by = "FL")

QLD_Beach <- read.csv("QLD Beach LF.csv", header = T)
names(QLD_Beach)[1] <- "FL"

all_data <- full_join(all_data2, QLD_Beach, by = "FL")
all_data[is.na(all_data)] <- 0

all_data$QLD_Rec_Percent <- all_data$QLD_Rec_Percent/100
all_data$QLD_Net_Percent <- all_data$QLD_Net_Percent/100
all_data$QLD_Beach_Percent <- all_data$QLD_Beach_Percent/100

all_data$weighted_percent <- ((107/309.12)*all_data$NSW_Rec_Percentage + (71.68/309.12)*all_data$NSW_Com_Percentage + (75/309.12)*all_data$QLD_Rec_Percent + (29.22/309.12)*all_data$QLD_Beach_Percent + (26.22/309.12)*all_data$QLD_Net_Percent)* (14461+850+10512+3681+2368) # number of fish measured total rec and commercial
sum(all_data$weighted_percent) # should equal 31872 (all fish measured)

all_data_long <- as.data.frame(rep(all_data$FL, all_data$weighted_percent))
names(all_data_long)[1] <- "FL"

write.csv(all_data_long, "East Australia Lengths weighted and Combined 2_cm.csv", row.names = F)

p1 <- ggplot(all_data_long, aes(x = FL)) + geom_histogram(aes(y=..density..),binwidth = 2) + theme_bw()
p1

