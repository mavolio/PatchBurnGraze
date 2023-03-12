################################################################################
##  PBG_resinbag_cleaning.R: Cleaning exported data from resin bag extract assays.
##
##  Author: Kimberly Komatsu
##  Date created: March 1, 2023
################################################################################

library(readxl)
library(tidyverse)

#set working directory
setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\patch burn')


#### import data ####
n <- read_xlsx('PBG_N_compiled_raw_2021.xlsx') %>% 
  filter(str_detect(sample, 'U')) %>% 
  group_by(test, Sample) %>% 
  filter(dilution==max(dilution)) %>% 
  ungroup() %>% 
  mutate(concentration=as.numeric(concentration))

NH4 <- n %>% 
  filter(test=='KCl Ammonia 10')

NO3 <- n %>% 
  filter(test=='KCL NO3_NO2 2')

p <- read_xlsx('PBG_P_compiled_raw_2021.xlsx') %>% 
  filter(str_detect(sample, 'U')) %>% 
  group_by(test, Sample) %>% 
  filter(dilution==max(dilution)) %>% 
  ungroup() %>% 
  mutate(concentration=as.numeric(concentration))


#### find outliers ####
hist(NH4$concentration)
hist(NO3$concentration)
hist(p$concentration)

plots <- rbind(NH4,NO3,p) %>% 
  select(Sample) %>% 
  unique() %>% 
  separate(Sample, into=c('watershed', 'transect', 'plot', 'nutrient'), sep='-') %>% 
  select(-nutrient) %>% 
  unique()

write.csv(plots, 'PBG_trts.csv', row.names=F)
