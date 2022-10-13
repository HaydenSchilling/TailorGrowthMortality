# Backcalculation of birth dates

#FL = 5.3
k = 0.02
L = 20.486
t = 94

D = (2*log(FL) + log(FL) - log(L) - k*t)/k
D

# Load data

mydata = read.csv("Final data.csv", header = T)


# subset - 
mydata <- subset(mydata, FL..cm. <= 20)
mydata <- subset(mydata, FL..cm. >= 5.3)

hist(mydata$FL..cm.)

library(lubridate)
mydata$Date.Caught <- dmy(mydata$Date.Caught)

plot(mydata$Date.Caught)

mydata$Days_old_est <- (2*log(mydata$FL..cm.) + log(mydata$FL..cm.) - log(L) - k*t)/k
plot(mydata$Days_old_est)

mydata$Birth <- mydata$Date.Caught - mydata$Days_old_est
mydata$Birth_year <- year(mydata$Birth)

plot(mydata$Birth)

library(ggplot2)
library(scales)
library(tidyverse)
mydata2 <- mydata %>% mutate(Birth_month = lubridate::month(Birth),
                             Birth_altered = dmy(paste0("15/",Birth_month,"/",Birth_year)))

summary <- mydata2 %>% group_by(Birth_year, Birth_month) %>% 
  summarise(n=n()) %>% ungroup() %>%
  complete(Birth_year, Birth_month,fill = list(n = 0)) %>%
  mutate(Birth_date = dmy(paste0("15/",Birth_month,"/",Birth_year))) %>%
  filter(Birth_date > "2014-07-01") %>% filter(Birth_date < "2016-04-01")

pb <- ggplot(summary, aes(x = Birth_date, y=n)) + geom_col(width=31) +
  scale_x_date(labels = date_format("%m-%Y"),
               breaks = summary$Birth_date) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 10)) +
  xlab("Estimated Birth Month") + ylab("Number of Fish")
  
pb

ggsave("Output/Estimated Birth dates.png", width = 21, height = 14.8, units = "cm", dpi = 600)




