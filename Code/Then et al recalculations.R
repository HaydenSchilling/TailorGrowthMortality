library(nlstools)
library(nlme)
source("R/functions.R") # this is a file with some of my customised functions in it, load by putting in working directory and running this line
library(dplyr)

mydata = read.csv("Then et al mortality data.csv", header=T)

mydata <- mydata[!is.na(mydata$tmax),] # removes rows with NA in tmax
mydata2 <- mydata[which(mydata$M <30),] # remove a massive outlier

# make the function i want to fit, from Then et al but recalculating it
Mest <- function(Tmax ,Scaling_factor, Tmax_exp)
{
    return (Scaling_factor*(Tmax^Tmax_exp))
}

fit1 <- nls(M ~ Mest(Tmax = tmax, Scaling_factor, Tmax_exp), start = list(Scaling_factor=4.9,Tmax_exp=-1),
            control=nls.control(maxiter=10000), data=mydata2)
summary(fit1) # Note same results and Model SE as Then et al = GOOD!

plot(fit1) # Not really a good residual plot but it was still used by Then et al

# To extract parameters
# fit1_Pars <- coef(fit1)
# fit1_Pars 

# generate Confidence Interval for set value of tmax (7) from loaded functions
Predictfit <- predictNLS(fit1, newdata = data.frame(tmax = 14), nsim = 100000)
Predictfit
# Take the output SimMAT data from the above function and convert to dataframe, 
# it contains combinations of tmax_exponential and the scaling factor and tmax (set to 7 for this)
SimMAT <- as.data.frame(SimMAT)
# Create a column for M in SimMAT then populate it from columns tmax and scaling factor and tmax
SimMAT$M <- numeric(length= nrow(SimMAT))
# calculate M for each combination
SimMAT$M <- SimMAT$Scaling_factor*(SimMAT$tmax^SimMAT$Tmax_exp)
summary(SimMAT$M)

# Histograms of results
M_hist <- hist(SimMAT$M, breaks = seq(0.7, 1, by = 0.01))
M_hist2 <- hist(SimMAT$M, breaks = seq(0.7, 1, by = 0.001))

# get list of bin midpoints
hist_data_midpoints <- M_hist$mids
# get density in each bin (I think this is likelihood, definately the number of results in bin)
hist_data_density <- M_hist$density
# bind these together
hist_data <- bind_cols(as.data.frame(hist_data_midpoints), as.data.frame(hist_data_density))
#write output file
write.csv(hist_data, "Output/Then_M_likelihoods.csv")

