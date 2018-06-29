setwd('C:\\Users\\megha\\Dropbox\\Grants\\USDA 2018\\data from konza\\plant composition')

library(tidyverse)
library(codyn)
library(vegan)

theme_set(theme_bw(12))

dat<-read.csv("PBG011.csv")

clean<-dat%>%
  mutate(treatment = ifelse(Watershed == "C01A"|Watershed == "C1SB", "A", "PB"))%>%
  mutate(replicate = ifelse(Watershed == "C01A", "A1", ifelse(Watershed == "C1SB", "A2", ifelse(Watershed == "C03A"|Watershed == "C03B"|Watershed == "C03C", "PB1", "PB2"))))%>%
  mutate(id = paste(Watershed, Transect, Plot, treatment, replicate, sep = "_"))%>%
  filter(RecYear > 2010)

#average all plots in a transect, average all transects in a watershed, average all watershed in a treatment

rich<-community_structure(clean, time.var = "RecYear", abundance.var = "CoverClass", replicate.var = "id")%>%
  separate(id, into=c("Watershed", "Transect", "Plot", "treatment", "replicate"), sep = "_")

ave_rich<-rich%>%
  group_by(RecYear, Watershed, Transect, treatment, replicate)%>%
  summarize(r1 = mean(richness))%>%
  ungroup()%>%
  group_by(RecYear, Watershed, treatment, replicate)%>%
  summarize(r2 = mean(r1))%>%
  ungroup()%>%
  group_by(RecYear, treatment, replicate)%>%
  summarize(richness = mean(r2))

ggplot(data = ave_rich, aes(x = as.factor(RecYear), y = richness, group = replicate, color = treatment))+
  geom_point()+
  geom_line()

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

ggplot(data = mean_rich, aes(x = as.factor(RecYear), y = mrich, color = treatment))+
  geom_point(size = 2, position=position_dodge(width = 0.2))+
  geom_errorbar(stat = 'identity', position = 'dodge', aes(ymin=mrich - se, ymax = mrich +se), width = 0.3)+
  scale_color_manual(name = "Treatment", labels = c("Annual Burn", "Patch Burn"), values = c("blue3","green3"))+
  geom_line(aes(group = treatment),position=position_dodge(width = 0.2))+
  ylab("Rarified Richness")+
  xlab("Year")
  

### NMDS and permanova

#looking at last year of data only
lastyr<-clean%>%
  filter(RecYear == 2017)

##average over plots becuase of a few errors
wide<-lastyr%>%
  mutate(species = paste(Ab_genus, Ab_species, sep = "_"))%>%
  group_by(RecYear, Watershed, Transect, treatment, replicate, species)%>%
  summarize(CoverClass = mean(CoverClass))%>%
  spread(species, CoverClass, fill = 0)

plotInfo<-wide[,1:5]

mds<-metaMDS(wide[,6:142], distance = "bray", trymax = 100)

mdsScores<-data.frame(scores(mds, display = "sites"))%>%
  bind_cols(plotInfo)

ggplot(data = mdsScores, aes(x = NMDS1, y = NMDS2, color = Watershed, shape = treatment))+
  geom_point(size = 8)

#permanova
adonis(wide[,6:142]~treatment, stata = wide$Watershed, wide)
