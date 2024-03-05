#### Libraries ####
library(readxl)

#### CSV READ ####

Weather_Data <- read.csv("AWE012.CSV")

#### Getting Means ####

mean_TMAX <- mean(as.numeric(Weather_Data$TMAX), na.rm = T)
mean_TMIN <- mean(as.numeric(Weather_Data$TMIN), na.rm = T)
mean_TAVE <- mean(as.numeric(Weather_Data$TAVE), na.rm = T)

# Print the means
print(mean_TMAX)
print(mean_TMIN)
print(mean_TAVE)

