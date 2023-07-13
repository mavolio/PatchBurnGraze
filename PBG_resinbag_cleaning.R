################################################################################
##  PBG_resinbag_cleaning.R: Cleaning exported data from resin bag extract assays.
##
##  Author: Kimberly Komatsu
##  Date created: March 1, 2023
################################################################################

library(readxl)
library(tidyverse)

#set working directory
setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\konza projects\\patch burn\\resin data')


#### import data ####
n2021 <- read_xlsx('2021\\PBG_N_compiled_raw_2021.xlsx') %>% 
  filter(str_detect(sample, 'U')) %>% 
  mutate(concentration=as.numeric(concentration)) %>% 
  separate(Sample, into=c('watershed', 'transect', 'plot', 'nutrient'), sep='-') %>% 
  group_by(watershed, transect, plot, test) %>% 
  filter(dilution==max(dilution)) %>% 
  ungroup() %>% 
  select(-nutrient) %>% 
  mutate(year=2021, nutrient=ifelse(test=='KCl Ammonia 10', 'NH4', 'NOx')) %>% 
  select(year, watershed, transect, plot, nutrient, units, concentration)

p2021 <- read_xlsx('2021\\PBG_P_compiled_raw_2021.xlsx') %>% 
  filter(str_detect(sample, 'U')) %>% 
  mutate(concentration=as.numeric(concentration)) %>% 
  separate(Sample, into=c('watershed', 'transect', 'plot', 'nutrient'), sep='-') %>% 
  group_by(watershed, transect, plot, test) %>% 
  filter(dilution==max(dilution)) %>% 
  ungroup() %>% 
  select(-nutrient) %>% 
  mutate(year=2021, nutrient='P') %>% 
  select(year, watershed, transect, plot, nutrient, units, concentration)


all2022 <- read_xlsx('2022\\NOx ammonia and P_all samples_compiled_2022.xlsx') %>% 
  filter(str_detect(sample, 'U')) %>% 
  mutate(concentration=as.numeric(concentration)) %>% 
  separate(Sample, into=c('expt', 'year', 'watershed', 'transect', 'plot', 'element', 'nutrient'), sep='_') %>% 
  group_by(watershed, transect, plot, test) %>% 
  filter(dilution==max(dilution)) %>% 
  ungroup() %>% 
  select(watershed, transect, plot, nutrient, units, concentration) %>% 
  mutate(year=2022) %>% 
  mutate(across('nutrient', str_replace, 'ammonia', 'NH4')) %>% 
  mutate(across('nutrient', str_replace, 'phosphorous', 'P'))


allData <- rbind(n2021, p2021, all2022)



#### find outliers ####

ggplot(data=allData, aes(x=concentration)) +
  geom_histogram() +
  facet_grid(year~nutrient, scales='free_x')


### what are the plot names? ####
plots <- allData %>% 
  select(watershed, transect, plot) %>% 
  unique()

# write.csv(plots, 'PBG_trts.csv', row.names=F)
