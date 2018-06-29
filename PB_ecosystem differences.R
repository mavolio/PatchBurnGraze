library(tidyverse)
library(codyn)
library(vegan)

theme_set(theme_bw(12))

plant<-read.csv("C:\\Users\\megha\\Dropbox\\Grants\\USDA 2018\\data from konza\\plant composition\\PBG011.csv")

dp_calibration<-read.csv("C:\\Users\\megha\\Dropbox\\Grants\\USDA 2018\\data from konza\\Production\\PBG031_Disk Pasture.csv")

dp<-read.csv("C:\\Users\\megha\\Dropbox\\Grants\\USDA 2018\\data from konza\\Production\\PBG032_0_Disk Pasture.csv")

grasshoppers<-read.csv("C:\\Users\\megha\\Dropbox\\Grants\\USDA 2018\\data from konza\\grasshopper data\\PBG081_0_density.csv")

birds<-read.csv("C:\\Users\\megha\\Dropbox\\Grants\\USDA 2018\\data from konza\\Birds\\PBG051.csv")

###Disk Pasture 

calibrate<-dp_calibration%>%
  mutate(live = (Lvgrass+Forbs+Woody)*10,
         total = (Lvgrass+Forbs+Pdead+Woody)*10)

#2011
summary(lm(Diskht~total-1, data = subset(calibrate, Recyear == 2011))) #hgt = 0.37519*biomass # 2011 is odd using the average of all other years instead. hgt = 0.028253*biomass
#2012
summary(lm(Diskht~total-1, data = subset(calibrate, Recyear == 2012))) #hgt = 0.019592*biomass
#2013
summary(lm(Diskht~total-1, data = subset(calibrate, Recyear == 2013))) #hgt = 0.036573*biomass
#2014
summary(lm(Diskht~total-1, data = subset(calibrate, Recyear == 2014))) #hgt = 0.026650*biomass
#2015
summary(lm(Diskht~total-1, data = subset(calibrate, Recyear == 2015))) #hgt = 0.024470*biomass
#2016
summary(lm(Diskht~total-1, data = subset(calibrate, Recyear == 2016))) #hgt = 0.027521*biomass
#2017
summary(lm(Diskht~live-1, data = subset(calibrate, Recyear == 2017))) #hgt = 0.034714*biomass

  
dp2<-dp%>%
  mutate(treatment = ifelse(Watershed == "c01a"|Watershed == "c1sb", "A", "PB"))%>%
  mutate(replicate = ifelse(Watershed == "c01A", "A1", ifelse(Watershed == "C1SB", "A2", ifelse(Watershed == "C03A"|Watershed == "C03B"|Watershed == "C03C", "PB1", "PB2"))))%>%
  mutate(id = paste(Watershed, Transect, Plot, treatment, replicate, sep = "_"))%>%
  filter(RecYear > 2010)%>%
  mutate(biomass = )

###plants
clean<-plant%>%
  mutate(treatment = ifelse(Watershed == "C01A"|Watershed == "C1SB", "A", "PB"))%>%
  mutate(replicate = ifelse(Watershed == "C01A", "A1", ifelse(Watershed == "C1SB", "A2", ifelse(Watershed == "C03A"|Watershed == "C03B"|Watershed == "C03C", "PB1", "PB2"))))%>%
  mutate(id = paste(Watershed, Transect, Plot, treatment, replicate, sep = "_"))%>%
  filter(RecYear > 2010)

##total species pool and rarefaction

sppool<-clean%>%
  mutate(species = paste(Ab_genus, Ab_species, sep = "_"))%>%
  mutate(rep_yr = paste(replicate, RecYear, treatment, sep = "_"))

#rarefraction
#create empty dataframe for loop
estimatedRichness=data.frame(row.names=1) 

rep_year=unique(sppool$rep_yr)

for(i in 1:length(rep_year)) {
  
  #creates a dataset for each replicate and averages all species in a replicate in case of duplicates
  subset <- sppool%>%
    filter(rep_yr==rep_year[i])%>%
    group_by(rep_yr, Watershed, treatment, replicate, Transect, Plot, species)%>%
    summarize(CoverClass = mean(CoverClass))
  
  #transpose data into wide form
  speciesData <- subset%>%
    spread(species, CoverClass, fill=0)
  
  #calculate species accumulation curves
  pool <- poolaccum(speciesData[,7:ncol(speciesData)], permutations=100)
  chao <- as.data.frame(as.matrix(pool$chao))#this gives us estimated richness from 1-X samples
  chao$aveChao<-rowMeans(chao)
  chao$n<-row.names(chao)
  chao$rep_year<-rep_year[i]
  chao2<-chao%>%
    select(rep_year,n, aveChao)
  
  #rbind back
  estimatedRichness<-rbind(chao2, estimatedRichness)
}

Richness<-estimatedRichness%>%
  filter(n==18)%>%#the lowest sampling intensity of the annually burned controls
  separate(rep_year, c("replicate", "RecYear", "treatment"), sep="_")%>%
  mutate(rrich=aveChao)%>%
  select(-n, -aveChao)

mean_rich<-Richness%>%
  group_by(treatment, RecYear)%>%
  summarize(mrich=mean(rrich), se=(sd(rrich)/sqrt(2)))

###


##graphing all of this
plants<-
ggplot(data = mean_rich, aes(x = as.factor(RecYear), y = mrich, color = treatment))+
  geom_point(size = 2, position=position_dodge(width = 0.2))+
  geom_errorbar(stat = 'identity', position = 'dodge', aes(ymin=mrich - se, ymax = mrich +se), width = 0.3)+
  scale_color_manual(name = "Treatment", labels = c("Annual Burn", "Patch Burn"), values = c("blue3","green3"))+
  geom_line(aes(group = treatment),position=position_dodge(width = 0.2))+
  ylab("Rarified Richness")+
  xlab("Year")
