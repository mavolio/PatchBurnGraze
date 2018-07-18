library(tidyverse)
library(codyn)
library(vegan)
library(gridExtra)

theme_set(theme_bw(12))
#work
plant<-read.csv("C:\\Users\\megha\\Dropbox\\Grants\\USDA 2018\\data from konza\\plant composition\\PBG011.csv")

dp_calibration<-read.csv("C:\\Users\\megha\\Dropbox\\Grants\\USDA 2018\\data from konza\\Production\\PBG031_Disk Pasture.csv")

dp<-read.csv("C:\\Users\\megha\\Dropbox\\Grants\\USDA 2018\\data from konza\\Production\\PBG032_0_Disk Pasture.csv")

grasshoppers<-read.csv("C:\\Users\\megha\\Dropbox\\Grants\\USDA 2018\\data from konza\\grasshopper data\\PBG081_0_density.csv")

birds<-read.csv("C:\\Users\\megha\\Dropbox\\Grants\\USDA 2018\\data from konza\\Birds\\PBG051.csv")

production<-read.csv("C:\\Users\\megha\\Dropbox\\Grants\\USDA 2018\\data from konza\\Production\\PBG021_ANPP.csv")

fire<-read.csv("C:\\Users\\megha\\Dropbox\\Grants\\USDA 2018\\data from konza\\KFH011.csv")

ysb<-read.csv("C:\\Users\\megha\\Dropbox\\Grants\\USDA 2018\\data from konza\\YearsSinceBurn.csv")

#home laptop
plant<-read.csv("~/Dropbox/Grants/USDA 2018/data from konza/plant composition/PBG011.csv")

dp_calibration<-read.csv("~/Dropbox/Grants/USDA 2018/data from konza/Production/PBG031_Disk Pasture.csv")

dp<-read.csv("~/Dropbox/Grants/USDA 2018/data from konza/Production/PBG032_0_Disk Pasture.csv")

grasshoppers<-read.csv("~/Dropbox/Grants/USDA 2018/data from konza/grasshopper data/PBG081_0_density.csv")

birds<-read.csv("~/Dropbox/Grants/USDA 2018/data from konza/Birds/PBG051.csv")

production<-read.csv("~/Dropbox/Grants/USDA 2018/data from konza/Production/PBG021_ANPP.csv")

fire<-read.csv("~/Dropbox/Grants/USDA 2018/data from konza/KFH011.csv")

ysb<-read.csv("~/Dropbox/Grants/USDA 2018/data from konza/YearsSinceBurn.csv")
ysb_birds<-read.csv("~/Dropbox/Grants/USDA 2018/data from konza/YearsSinceBurn_birds.csv")

###Disk Pasture 
calibrate<-dp_calibration%>%
  mutate(live = (Lvgrass+Forbs+Woody)*10,
         total = (Lvgrass+Forbs+Pdead+Woody)*10)

#2011
summary(lm(Diskht~total-1, data = subset(calibrate, Recyear == 2011))) #hgt = 0.37519*biomass 
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
  mutate(ws = toupper(Watershed))%>%
  mutate(treatment = ifelse(ws == "C01A"|ws == "C1SB", "A", "PB"))%>%
  mutate(replicate = ifelse(ws == "C01A", "A1", ifelse(ws == "C1SB", "A2", ifelse(ws == "C03A"|ws == "C03B"|ws == "C03C", "PB1", "PB2"))))%>%
  filter(Recyear > 2010)%>%
  mutate(hgt = (ADiskHT + BDiskHT)/2,
         biomass = ifelse(Recyear == 2011, hgt/0.37519, 
                   ifelse(Recyear == 2012, hgt/0.019592, 
                   ifelse(Recyear == 2013, hgt/0.036573, 
                   ifelse(Recyear == 2014, hgt/0.026650,
                   ifelse(Recyear == 2015, hgt/0.024470,
                   ifelse(Recyear == 2016, hgt/0.027521, hgt/0.034714)))))))

dp_ave<-dp2%>%
  group_by(ws, Recyear, treatment, replicate, Transect)%>%
  summarize(biomass2 = mean(biomass))%>%
  ungroup()%>%
  group_by(ws, Recyear, treatment, replicate)%>%
  summarize(biomass3 = mean(biomass2))%>%
  ungroup()%>%
  group_by(treatment, Recyear, replicate)%>%
  summarize(biomass4 = mean(biomass3))%>%
  ungroup()%>%
  group_by(Recyear, treatment)%>%
  summarize(ave = mean(biomass4),
            sd = sd(biomass4))%>%
  mutate(se = sd / sqrt(2))

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

###grasshoppers
gh2<-grasshoppers%>%
  filter(Watershed != "001D")%>%
  filter(Watershed != "004B")%>%
  filter(Recyear > 2010)%>%
  mutate(treatment = ifelse(Watershed == "C01A"|Watershed == "C1SB", "A", "PB"))%>%
  mutate(replicate = ifelse(Watershed == "C01A", "A1", ifelse(Watershed == "C1SB", "A2", ifelse(Watershed == "C03A"|Watershed == "C03B"|Watershed == "C03C", "PB1", "PB2"))))%>%
  group_by(Recyear, Watershed, Site, treatment, replicate)%>%
  summarize(num = sum(Count))%>%
  filter(Watershed != "C01B")

gh_ave<-gh2%>%
  group_by(Watershed, Recyear, treatment, replicate, Site)%>%
  summarize(num2 = mean(num))%>%
  ungroup()%>%
  group_by(Watershed, Recyear, treatment, replicate)%>%
  summarize(num3 = mean(num2))%>%
  ungroup()%>%
  group_by(treatment, Recyear, replicate)%>%
  summarize(num4 = mean(num3))%>%
  ungroup()%>%
  group_by(Recyear, treatment)%>%
  summarize(ave = mean(num4),
            sd = sd(num4))%>%
  mutate(se = sd / sqrt(2))

##birds
birds2<-birds%>%
  filter(Watershed != "1D")%>%
  mutate(treatment = ifelse(Watershed == "C1A","A", "PB"))%>%
  mutate(present = ifelse(Vdetection == 1 |Sdetection == 1|Cdetection == 1|Fly == 1|Flush == 1, 1, 0))%>%
  mutate(num = ifelse(is.na(GroupSize), 1, GroupSize))%>%
  mutate(count = present*num)%>%
  filter(!is.na(count))

birds_ave<-birds2%>%
  group_by(Recyear, Watershed, Transect, Direction, DayofYear, treatment)%>%
  summarize(number = sum(count))%>%
  group_by(Recyear, Watershed, Transect, treatment)%>%
  summarize(number2 = mean(number))%>%
  group_by(Recyear, Watershed, treatment)%>%
  summarise(number3 = mean(number2))%>%
  group_by(Recyear, treatment)%>%
  summarize(ave = mean(number3))
  

##grass forb
prod2<-production%>%
  mutate(ws = toupper(Watershed))%>%
  mutate(treatment = ifelse(ws == "C01A"|ws == "C1SB", "A", "PB"))%>%
  mutate(replicate = ifelse(ws == "C01A", "A1", ifelse(ws == "C1SB", "A2", ifelse(ws == "C03A"|ws == "C03B"|ws == "C03C", "PB1", "PB2"))))%>%
  filter(RecYear > 2010)%>%
  filter(Treatment == "u")%>%
  mutate(gf = Forbs/Lvgrass)

gf_ave<-prod2%>%
  group_by(ws, RecYear, Cage, treatment, replicate)%>%
  summarize(gf2 = mean(gf))%>%
  ungroup()%>%
  group_by(ws, RecYear, treatment, replicate)%>%
  summarize(gf3 = mean(gf2))%>%
  ungroup()%>%
  group_by(RecYear, treatment, replicate)%>%
  summarize(gf4 = mean(gf3))%>%
  ungroup()%>%
  group_by(RecYear, treatment)%>%
  summarize(ave = mean(gf4),
            sd = sd(gf4))%>%
  mutate(se = sd / sqrt(2))
  
  


##graphing all of this
plants<-
ggplot(data = mean_rich, aes(x = as.factor(RecYear), y = mrich, fill = treatment))+
  geom_bar(position=position_dodge(), stat = "identity")+
  geom_errorbar(stat = 'identity', position = position_dodge(0.9), aes(ymin=mrich - se, ymax = mrich +se), width = 0.3, size = 1)+
  scale_fill_manual(name = "Treatment", labels = c("Annual Burn", "Patch Burn"), values = c("chocolate","cornflowerblue"))+
  ylab("Rarified Richness")+
  xlab("Year")

diskpasture<-
  ggplot(data = dp_ave, aes(x = as.factor(Recyear), y = ave, fill = treatment))+
  geom_bar(position=position_dodge(), stat = "identity")+
  geom_errorbar(stat = 'identity', position = position_dodge(0.9), aes(ymin=ave - se, ymax = ave +se), width = 0.3, size = 1)+
  scale_fill_manual(name = "Treatment", labels = c("Annual Burn", "Patch Burn"), values = c("chocolate","cornflowerblue"))+
  ylab("Plant Biomass (g m-2)")+
  xlab("Year")

grasshops<-
  ggplot(data = gh_ave, aes(x = as.factor(Recyear), y = ave, fill = treatment))+
  geom_bar(position=position_dodge(), stat = "identity")+
  geom_errorbar(stat = 'identity', position = position_dodge(0.9), aes(ymin=ave - se, ymax = ave +se), width = 0.3, size = 1)+
  scale_fill_manual(name = "Treatment", labels = c("Annual Burn", "Patch Burn"), values = c("chocolate","cornflowerblue"))+
  ylab("Number of Grasshoppers")+
  xlab("Year")

birdz<-
  ggplot(data = birds_ave, aes(x = as.factor(Recyear), y = ave, fill = treatment))+
  geom_bar(position=position_dodge(), stat = "identity")+
  scale_fill_manual(name = "Treatment", labels = c("Annual Burn", "Patch Burn"), values = c("chocolate","cornflowerblue"))+
  ylab("Number of Birds")+
  xlab("Year")

grid.arrange(birdz, grasshops, plants, diskpasture, ncol = 1 )

gf<-
  ggplot(data = gf_ave, aes(x = as.factor(RecYear), y = ave, color = treatment))+
  geom_point(size = 2, position=position_dodge(width = 0.2))+
  geom_errorbar(stat = 'identity', position = 'dodge', aes(ymin=ave - se, ymax = ave +se), width = 0.3)+
  scale_color_manual(name = "Treatment", labels = c("Annual Burn", "Patch Burn"), values = c("blue3","green3"))+
  geom_line(aes(group = treatment),position=position_dodge(width = 0.2))+
  ylab("Forb to Grass Ratio")+
  xlab("Year")

###looking at fire.

#use this to get fire data
fire2<-fire%>%
  filter(Code == "C01A"|Code == "C1SB"|Code == "C03A"|Code == "C03B"|Code == "C03C"|Code == "C3SA"|Code =="C3SB"|Code == "C3SC")%>%
  filter(Year > 2008)%>%
  mutate(fire = 1, Watershed = Code, Recyear = Year)


#dp

dp3<-dp2%>%
  select(-Watershed)%>%
  mutate(Watershed = ws)%>%
  left_join(ysb)

dp_ave_patch<-dp3%>%
  group_by(Watershed, Recyear, ysb, treatment, Transect)%>%
  summarize(biomass2 = mean(biomass))%>%
  ungroup()%>%
  group_by(Watershed, ysb, Recyear, treatment)%>%
  summarize(biomass3 = mean(biomass2))%>%
  ungroup()%>%
  group_by(ysb, treatment)%>%
  summarize(ave = mean(biomass3),
            sd = sd(biomass3))%>%
  mutate(se = sd / sqrt(2))%>%
  mutate(response = "Standing Biomass")

##plant rich
rich<-community_structure(clean, time.var = "RecYear", abundance.var = "CoverClass", replicate.var = "id")%>%
  separate(id, into=c("Watershed", "Transect", "Plot", "treatment", "replicate"), sep = "_")%>%
  mutate(Recyear = RecYear)

ave_rich<-rich%>%
  left_join(ysb)%>%
  group_by(Recyear, Watershed, Transect, treatment, ysb)%>%
  summarize(r1 = mean(richness))%>%
  ungroup()%>%
  group_by(Recyear, Watershed, treatment, ysb)%>%
  summarize(r2 = mean(r1))%>%
  ungroup()%>%
  group_by(ysb, treatment)%>%
  summarize(ave = mean(r2),
            sd = sd(r2))%>%
  mutate(se = sd / sqrt(2))%>%
  mutate(response = "Plant Richness")

#grasshoppers
gh_ave_patch<-gh2%>%
  left_join(ysb)%>%
  filter(Watershed != "C01B")%>%
  group_by(Watershed, Recyear, treatment, ysb, Site)%>%
  summarize(num2 = mean(num))%>%
  ungroup()%>%
  group_by(Watershed, Recyear, treatment, ysb)%>%
  summarize(num3 = mean(num2))%>%
  ungroup()%>%
  group_by(ysb, treatment)%>%
  summarize(ave = mean(num3),
            sd = sd(num3))%>%
  mutate(se = sd / sqrt(2))%>%
  mutate(response = "Grasshoppers")

#birds
birds_ave_patch<-birds2%>%
  left_join(ysb_birds)%>%
  group_by(Recyear, Watershed, Transect, Direction, DayofYear, treatment, ysb)%>%
  summarize(number = sum(count))%>%
  group_by(Recyear, Watershed, Transect, treatment, ysb)%>%
  summarize(number2 = mean(number))%>%
  group_by(Recyear, Watershed, treatment, ysb)%>%
  summarise(number3 = mean(number2))%>%
  group_by(ysb, treatment)%>%
  summarize(ave = mean(number3),
            sd = sd(number3))%>%
  mutate(se = sd / sqrt(2))%>%
  mutate(response = "Birds")

patch_all<-rbind(dp_ave_patch, birds_ave_patch, gh_ave_patch, ave_rich)%>%
  mutate(ysb2 = paste(treatment, ysb, sep = ""))

ggplot(data=patch_all, aes(x = ysb2, y = ave, fill = treatment))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = ave - se, ymax = ave + se), width = 0.3)+
  scale_fill_manual(name = "Treatment", label=c("Annual Burn", "Patch Burn"), values = c("chocolate", "cornflowerblue"))+
  xlab("Years Since Burn")+
  ylab("")+
  facet_wrap(~response, ncol= 1, scale="free")
  
