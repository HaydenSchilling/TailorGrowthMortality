library(FSA)

# #Sample Data
# ( bkt <- data.frame(age=0:6,ct=c(39,93,112,45,58,12,8)) )
# bkt$logct <- log(bkt$ct)
# str(bkt)
# thcr <- chapmanRobson(ct~age,data=bkt,ages2use=2:6)
#  summary(thcr)
#  confint(thcr) 

 
 # Commercial NSW Estimate
 mydata <- read.csv("NSW Catch for Chapman Robson.csv", header = T) 
str(mydata) 

CR_Z <- chapmanRobson(Catch~Age,data=mydata,ages2use=2:5)
summary(CR_Z)
confint(CR_Z)
plot(CR_Z)

cite()
cite

# Recreational NSW Estimate

mydata <- read.csv("NSW Rec Catch for Chapman Robson.csv", header = T) 
str(mydata)

CR_Z <- chapmanRobson(Catch~Age,data=mydata,ages2use=3:7) # 1 year older than age at peak abundance
summary(CR_Z)
confint(CR_Z)

# Tagging data catch
mydata <- read.csv("Tagged Catch for Chapman Robson.csv", header = T)
str(mydata) 

CR_Z <- chapmanRobson(Catch~Age,data=mydata,ages2use=2:6) # 1 year older than age at peak abundance
summary(CR_Z)
confint(CR_Z)

# NSW combined catch
mydata <- read.csv("NSW Combined Catch for Chapman Robson.csv", header = T) 
str(mydata) 

CR_Z <- chapmanRobson(Catch~Age,data=mydata,ages2use=3:6)
summary(CR_Z)
confint(CR_Z)
plot(CR_Z)


# East Australia combined catch
mydata <- read.csv("East Aus Combined Catch for Chapman Robson.csv", header = T) 
str(mydata) 

CR_Z <- chapmanRobson(Catch~Age,data=mydata,ages2use=3:7)
summary(CR_Z)
confint(CR_Z)
plot(CR_Z)

# NSW Historical
mydata <- read.csv("NSW Historical Catch for Chapman Robson.csv", header = T) 
str(mydata) 

CR_Z <- chapmanRobson(Catch~Age,data=mydata,ages2use=3:6)
summary(CR_Z)
confint(CR_Z)
plot(CR_Z)

# Mediteranean from Cehyan 2007
mydata <- read.csv("MED catch for Chapman Robson.csv", header = T) 
str(mydata) 

CR_Z <- chapmanRobson(Catch~Age,data=mydata,ages2use=2:3)
summary(CR_Z)
confint(CR_Z)
plot(CR_Z)

# Mediteranean from Cengiz 2013
mydata <- read.csv("MED catch for Chapman Robson - copy.csv", header = T) 
str(mydata) 

CR_Z <- chapmanRobson(Catch~Age,data=mydata,ages2use=1:3)
summary(CR_Z)
confint(CR_Z)
plot(CR_Z)


cite()
cite
