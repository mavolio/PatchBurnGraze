#### Libraries ####
library(readxl)

#### CSV READ ####

Weather_Data <- read.csv("AWE012.CSV")

#### Getting Temp Means ####

mean_TMAX <- mean(as.numeric(Weather_Data$TMAX), na.rm = T)
mean_TMIN <- mean(as.numeric(Weather_Data$TMIN), na.rm = T)
mean_TAVE <- mean(as.numeric(Weather_Data$TAVE), na.rm = T)

# Print the means
print(mean_TMAX)
print(mean_TMIN)
print(mean_TAVE)

#### Getting Average Total Precip Mean ####

Weather_Data$DPPT <- as.numeric(Weather_Data$DPPT, na.rm = T)

# Aggregate the total DPPT for each year
yearly_DPPT <- aggregate(DPPT ~ RECYEAR, data = Weather_Data, FUN = sum, na.rm = TRUE)

# Calculate the mean of yearly DPPT
mean_yearly_DPPT <- mean(yearly_DPPT$DPPT)
print(mean_yearly_DPPT)

