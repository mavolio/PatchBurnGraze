#PBG Bird count and community metrics
#Authors: Joshua Ajowele
#Started: 26 May 2022 last modified: 21 Feb 2025

#load library####
library(tidyverse)
library(vegan)
library(codyn)
library(stringr)
library(ggthemes)
library(readr)
library(readxl)
library(performance)
library(car)
library(lme4)
library(lmerTest)
library(see)
library(patchwork)
library(phia)



### Standard Error function####
SE_function<-function(x,na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}

#Set graphing parameters####
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             birdot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=20), legend.text=element_text(size=20))

#import bird data####
PBG_bird_data_raw<-read_excel("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/PhD Wyoming_One Drive/PHD Wyoming/Thesis/PBG synthesis/PBG_bird_data_Alice_edited/Qy_ExportPBGSurveyData_Mar2023.xlsx")
PBG_bird_C1A<-read_excel("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/PhD Wyoming_One Drive/PHD Wyoming/Thesis/PBG synthesis/PBG_bird_data_Alice_edited/Qy_ExportPBGSurveyData_C1A_Apr2023.xlsx")
unique(PBG_bird_C1A$Year)
#2019 data is missing for C1A-Not collected
PBG_bird_class<-read_excel("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/PhD Wyoming_One Drive/PHD Wyoming/Thesis/PBG synthesis/PBG_bird_data_Alice_edited/KNZSpecies.xlsx")

#merge PBG and ABG data and remove unwanted rows####
PBG_bird_merge<-PBG_bird_data_raw%>%
  rename(SpeciesCode=AOUCODE)%>%
  full_join(PBG_bird_C1A, by=c("Year","TransectName","Duration","Date","Observer","StartTime",
                               "Watershed", "DistanceFromLine_Old","SpeciesCode","Count","SEX",
                               "ObsLocation","Notes","Angle","AngularDistance"))%>%
  filter(Year%in%2011:2021)%>%
  filter(Watershed!="C3SC")%>%
  filter(Watershed!="C3SB")%>%
  filter(Watershed!="C3SA")%>%
  filter(Year!=2019)
unique(PBG_bird_merge$TransectName)

#Create a watershed key column to merge with the raw data####
watershed_bird_key <- data.frame(Watershed = levels(factor(PBG_bird_merge$Watershed)),# create key to merge with raw data
                                 FireGrzTrt=c("ABG", "PBG", "PBG", "PBG"))

#need to create transect key for easier statistical analysis
transect_key<- data.frame(TransectName = levels(factor(PBG_bird_merge$TransectName)),# create key to merge with raw data
                          Transect=c("A","B","A","B","A","B","A","B"))
#merge key
PBG_bird_viz<-PBG_bird_merge%>%
  left_join(watershed_bird_key, by="Watershed")%>%
  left_join(transect_key, by="TransectName")

#There are multibirde observation by multibirde observers 2011-2016
#or just multibirde observation by one observer 2011, C3C, C1A
#solutions:group by time, species and take max count, regardless of observer
check_data<-PBG_bird_viz%>%
  filter(Year==2016)%>%
  filter(Watershed=="C3A")
#check unique species code name
codes_name<-data.frame(unique(PBG_bird_viz$SpeciesCode))
#convert to uppercase
PBG_bird_viz$SpeciesCode<-toupper(PBG_bird_viz$SpeciesCode)

check_data1<-PBG_bird_viz%>%
  filter(Year==2013)%>%
  filter(Watershed=="C1A")

#year since fire 
#community metrics
#incorporating year since fire
Yrs_watershed_bird <- PBG_bird_viz%>%
  rename(total=Count)%>%
  filter(total!="NA")%>%
  left_join(PBG_bird_class, by="SpeciesCode")%>%
  #creating a coloumn for unique combination of year and watershed
  mutate(year_watershed = paste(Year, Watershed, sep = "_"))

#creating a key for year since fire
YrSinceFire_key_bird <- tibble(year_watershed = c("2011_C1A", "2011_C3A", "2011_C3B","2011_C3C",
                                                  "2012_C1A", "2012_C3A", "2012_C3B","2012_C3C",
                                                  "2013_C3A","2013_C3B", "2013_C3C",  "2013_C1A",
                                                  "2014_C1A", "2014_C3A","2014_C3B", "2014_C3C",
                                                  "2015_C3A", "2015_C3B","2015_C3C", "2015_C1A", 
                                                  "2016_C3A","2016_C3B",  "2016_C1A","2016_C3C", 
                                                  "2017_C3C", "2017_C3A", "2017_C3B","2017_C1A",
                                                  "2018_C3A", "2018_C3B", "2018_C3C","2018_C1A",
                                                  "2019_C3A", "2019_C3B", "2019_C3C","2019_C1A",
                                                  "2020_C3A","2020_C3B", "2020_C3C", "2020_C1A",
                                                  "2021_C1A","2021_C3C","2021_C3A","2021_C3B"),
                               yrsince_fire = c("ABG0","PBG1","PBG0", "PBG2",
                                                "ABG0","PBG2","PBG1","PBG0",
                                                "PBG0", "PBG2", "PBG1", "ABG0",
                                                "ABG0","PBG1", "PBG0","PBG2",
                                                "PBG2","PBG1","PBG0","ABG0",
                                                "PBG0","PBG2","ABG0","PBG1",
                                                "PBG2","PBG1","PBG0","ABG0",
                                                "PBG2","PBG1","PBG0","ABG0",
                                                "PBG0","PBG2","PBG1","ABG0",
                                                "PBG1","PBG0","PBG2","ABG0",
                                                "ABG0","PBG0","PBG2","PBG1"))

#joining the dataset
Yrs_since_fire_bird <- Yrs_watershed_bird %>%
  left_join(YrSinceFire_key_bird, by = "year_watershed")%>%
  filter(!Year=="2019")


#total bird count for year since fire
yrsincef_bird_anlys<-Yrs_since_fire_bird%>%
  group_by(Year, yrsince_fire, Watershed, TransectName, Transect,SpeciesCode)%>%
  summarise(total_max=max(total, na.rm=T))%>%
  group_by(Year, yrsince_fire, Watershed, TransectName, Transect)%>%
  summarise(tot_max=sum(total_max, na.rm=T))
yrsincef_bird_anlys$Year<-as.factor(yrsincef_bird_anlys$Year)
yrsincef_bird_anlys$yrsince_fire<-as.factor(yrsincef_bird_anlys$yrsince_fire)
#visuals
yrsincef_bird_viz<-yrsincef_bird_anlys%>%
  group_by(Year, yrsince_fire)%>%
  summarise(total=mean(tot_max, na.rm=T),
            se_total=SE_function(tot_max))

#total count 
ggplot(yrsincef_bird_viz, aes(Year, total, col = yrsince_fire,
                              fill=yrsince_fire, linetype=yrsince_fire)) +
  geom_point(size=2, col="black") +
  geom_path(aes(as.numeric(Year))) +
  geom_errorbar(aes(ymin = total - se_total, 
                    ymax = total + se_total),
                width=0.1, col="black", linetype=1) +
  scale_colour_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab("bird count")


#models for total count 
yrsincef_bird_model<-lmer(log(tot_max)~yrsince_fire+Year+(1|Transect), data=yrsincef_bird_anlys)
check_model(yrsincef_bird_model)
anova(yrsincef_bird_model, type=3)
#use lm since randmon effect bariance is zero
yrsincef_bird_model_lm<-lm(log(tot_max)~yrsince_fire+Year, data=yrsincef_bird_anlys)
summary(yrsincef_bird_model_lm)
anova(yrsincef_bird_model_lm)
check_model(yrsincef_bird_model_lm)
qqnorm(resid(yrsincef_bird_model_lm))
#using lm or lmer produces the same F value and P value


#visuals for treatment without interaction
yrsincef_bird_viz_yr<-yrsincef_bird_anlys%>%
  group_by(yrsince_fire)%>%
  summarise(total=mean(tot_max, na.rm=T),
            se_total=SE_function(tot_max))

ggplot(yrsincef_bird_viz_yr, aes(yrsince_fire, total, fill=yrsince_fire))+
  geom_col(width=.5) +
  geom_errorbar(aes(ymax=total+se_total, ymin=total-se_total),width=0.1)+
  scale_fill_manual(values=c("#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab("Bird count")

#year and trt significant. needs multicomparison
#use the ghlt for multicomparison
#install.packages("multcomp")
library(multcomp)
yrsincefire_comp<- glht(yrsincef_bird_model_lm, linfct=mcp(yrsince_fire="Tukey"))
print(summary(yrsincefire_comp))
#ALL P-ADJUSTED 

#visuals
yrsincef_bird_viz_yr<-yrsincef_bird_anlys%>%
  group_by(yrsince_fire, Year)%>%
  summarise(total=mean(tot_max, na.rm=T),
            se_total=SE_function(tot_max))%>%
  group_by(yrsince_fire)%>%
  summarise(total=mean(total, na.rm=T),
            se_total=mean(se_total, na.rm=T))

ggplot(yrsincef_bird_viz_yr, aes(yrsince_fire, total, fill=yrsince_fire))+
  geom_col(width=.5) +
  geom_errorbar(aes(ymax=total+se_total, ymin=total-se_total),width=0.1)+
  scale_fill_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab("Bird count")

#community metrics with all species####
yrsin_bird_comm<-Yrs_since_fire_bird%>%
  group_by(Year, yrsince_fire, Watershed, TransectName, Transect,SpeciesCode)%>%
  summarise(total_max=max(total, na.rm=T))%>%
  group_by(Year, yrsince_fire, Watershed, TransectName, Transect)%>%
  mutate(abund=(total_max/sum(total_max))*100)%>%
  mutate(rep_id=paste(Watershed, Transect, sep="_"))
  
#community metrics with codyn####
#deriving richness and evenness using codyn at landscape scale
bird_rich <- community_structure(yrsin_bird_comm, time.var = "Year", 
                               abundance.var = "abund",
                               replicate.var = "rep_id", metric = "Evar")

#deriving diversity values
bird_diver <- community_diversity(yrsin_bird_comm, time.var = "Year", 
                                abundance.var = "abund",
                                replicate.var = "rep_id", metric="Shannon")


#extracting important columns 
bird_comp_subset <- yrsin_bird_comm %>%
  ungroup()%>%
  dplyr::select(Year, yrsince_fire, Watershed, TransectName, Transect,rep_id) %>%
  distinct()#remove repeated rows

#combine into a single data
bird_rich_diver<- bird_diver %>%
  #combine richness and diversity
  left_join(bird_rich, by = c("rep_id","Year")) %>%
  #combine with important columns
  left_join(bird_comp_subset, by = c("rep_id","Year"))  
#convert to factors
bird_rich_diver$Year<-as.factor(bird_rich_diver$Year)
bird_rich_diver$Transect<-as.factor(bird_rich_diver$Transect)
bird_rich_diver$yrsince_fire<-as.factor(bird_rich_diver$yrsince_fire)
#build model
#richness
bird_rich_model<-lmer(log(richness)~yrsince_fire+Year+(1|Transect),
                       data=bird_rich_diver)
anova(bird_rich_model)
summary(bird_rich_model)
qqnorm(resid(bird_rich_model))
check_model(bird_rich_model)
yrsincefire_rich<- glht(bird_rich_model, linfct=mcp(yrsince_fire="Tukey"))
print(summary(yrsincefire_rich))

#Evenness
bird_evar_model<-lmer(Evar~yrsince_fire+Year+(1|Transect),
                      data=bird_rich_diver)
anova(bird_evar_model)
summary(bird_evar_model)
qqnorm(resid(bird_evar_model))
check_model(bird_evar_model)
#Visual
#dataframe for geompoint
bird_rich_diver_viz<-bird_rich_diver%>%
  group_by(Year, yrsince_fire)%>%
  summarise(rich=mean(richness,na.rm=T),
            rich_se=SE_function(richness),
            diver=mean(Shannon,na.rm=T),
            diver_se=SE_function(Shannon),
            evar=mean(Evar,na.rm=T),
            evar_se=SE_function(Evar))
#dataframe for bargraph
bird_rich_diver_bar<-bird_rich_diver%>%
  group_by(Year,yrsince_fire)%>%
  summarise_at(vars(Shannon:Evar),mean,na.rm=T)%>%
  mutate(SE_rich=SE_function(richness),
            SE_evar=SE_function(Evar))%>%
group_by(yrsince_fire)%>%
  summarise(rich=mean(richness,na.rm=T),
            rich_se=mean(SE_rich),
            evar=mean(Evar,na.rm=T),
            evar_se=mean(SE_evar))
 
#richness
#geompoint 
ggplot(bird_rich_diver_viz, aes(Year, rich, col = yrsince_fire,
                                fill=yrsince_fire, linetype=yrsince_fire)) +
  geom_point(size=2, col="black") +
  geom_path(aes(as.numeric(Year))) +
  geom_errorbar(aes(ymin = rich - rich_se, 
                    ymax = rich + rich_se),
                width=0.1, col="black", linetype=1) +
  scale_colour_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab("Bird richness")
#bargraph
ggplot(bird_rich_diver_bar, aes(yrsince_fire, rich, fill=yrsince_fire))+
  geom_col(width=.5) +
  geom_errorbar(aes(ymin = rich - rich_se, 
                    ymax = rich + rich_se),width=0.1)+
  scale_fill_manual(values=c("#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab("Bird richness")

#Evar
#geompoint 
ggplot(bird_rich_diver_viz, aes(Year, evar, col = yrsince_fire,
                                fill=yrsince_fire, linetype=yrsince_fire)) +
  geom_point(size=2, col="black") +
  geom_path(aes(as.numeric(Year))) +
  geom_errorbar(aes(ymin = evar - evar_se, 
                    ymax = evar + evar_se),
                width=0.1, col="black", linetype=1) +
  scale_colour_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab("Bird evenness")
#bargraph
ggplot(bird_rich_diver_bar, aes(yrsince_fire, evar, fill=yrsince_fire))+
  geom_col(width=.5) +
  geom_errorbar(aes(ymin = evar - evar_se, 
                    ymax = evar + evar_se),width=0.1)+
  scale_fill_manual(values=c("#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab("Bird evenness")

#grassland bird species
#total grass bird count for year since fire
yrsincef_grass_bird_anlys<-Yrs_since_fire_bird%>%
  filter(Grassland=="TRUE")%>%
  group_by(Year, yrsince_fire, Watershed, TransectName, Transect,SpeciesCode)%>%
  summarise(total_max=max(total, na.rm=T))%>%
  group_by(Year, yrsince_fire, Watershed, TransectName, Transect)%>%
  summarise(tot_count=sum(total_max, na.rm=T))
yrsincef_grass_bird_anlys$Year<-as.factor(yrsincef_grass_bird_anlys$Year)
yrsincef_grass_bird_anlys$yrsince_fire<-as.factor(yrsincef_grass_bird_anlys$yrsince_fire)

#visual
yrsincef_grass_bird_viz<-yrsincef_grass_bird_anlys%>%
  group_by(Year,yrsince_fire)%>%
  summarise(total=mean(tot_count, na.rm=T),
            se_total=SE_function(tot_count))

#total grassland bird count 
ggplot(yrsincef_grass_bird_viz, aes(Year, total, col = yrsince_fire,
                                    fill=yrsince_fire, linetype=yrsince_fire)) +
  geom_point(size=2) +
  geom_path(aes(as.numeric(Year))) +
  geom_errorbar(aes(ymin = total - se_total, 
                    ymax = total + se_total),
                width=0.1, col="black", linetype=1) +
  scale_colour_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2")) +
  ylab("grassland bird count")

#model
yrsincef_grass_bird_model<-lmer(log(tot_count)~yrsince_fire+Year+(1|Transect), data=yrsincef_grass_bird_anlys)
anova(yrsincef_grass_bird_model)
check_model(yrsincef_grass_bird_model)#no significance
qqnorm(resid(yrsincef_grass_bird_model))
summary(yrsincef_grass_bird_model)
#using linear model
yrsincef_grass_bird_model_lm<-lm(log(tot_count)~yrsince_fire+Year, data=yrsincef_grass_bird_anlys)
anova(yrsincef_grass_bird_model_lm)
#lmer and lm are the same

#visuals across years
yrsincef_grass_bird_yr<-yrsincef_grass_bird_anlys%>%
  group_by(yrsince_fire,Year)%>%
  summarise(total=mean(tot_count, na.rm=T),
            se_total=SE_function(tot_count))%>%
  group_by(yrsince_fire)%>%
  summarise(total=mean(total,na.rm=T),
            se_total=mean(se_total,na.rm=T))

ggplot(yrsincef_grass_bird_yr, aes(yrsince_fire, total, fill=yrsince_fire))+
  geom_col(width=.5) +
  geom_errorbar(aes(ymax=total+se_total, ymin=total-se_total),width=0.1)+
  scale_fill_manual(values=c("#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab("grassland Bird count")

#community metrics for grassland species
yrsin_grass_bird_comm<-Yrs_since_fire_bird%>%
  filter(Grassland=="TRUE")%>%
  group_by(Year, yrsince_fire, Watershed, TransectName, Transect,SpeciesCode)%>%
  summarise(total_max=max(total, na.rm=T))%>%
  group_by(Year, yrsince_fire, Watershed, TransectName, Transect)%>%
  mutate(abund=(total_max/sum(total_max))*100,
         rep_id=paste(Watershed,Transect,sep="_"))
#community metrics with codyn####
#deriving richness and evenness using codyn at landscape scale
bird_rich_grass <- community_structure(yrsin_grass_bird_comm, time.var = "Year", 
                                 abundance.var = "abund",
                                 replicate.var = "rep_id", metric = "Evar")

#deriving diversity values
bird_diver_grass <- community_diversity(yrsin_grass_bird_comm, time.var = "Year", 
                                  abundance.var = "abund",
                                  replicate.var = "rep_id", metric="Shannon")


#extracting important columns 
bird_comp_subset_grass <- yrsin_grass_bird_comm %>%
  ungroup()%>%
  dplyr::select(Year, yrsince_fire, Watershed, TransectName, Transect,rep_id) %>%
  distinct()#remove repeated rows

#combine into a single data
bird_rich_diver_grass<- bird_diver_grass %>%
  #combine richness and diversity
  left_join(bird_rich_grass, by = c("rep_id","Year")) %>%
  #combine with important columns
  left_join(bird_comp_subset_grass, by = c("rep_id","Year"))  
#convert to factors
bird_rich_diver_grass$Year<-as.factor(bird_rich_diver_grass$Year)
bird_rich_diver_grass$Transect<-as.factor(bird_rich_diver_grass$Transect)
bird_rich_diver_grass$yrsince_fire<-as.factor(bird_rich_diver_grass$yrsince_fire)
#build model

#richness
bird_rich_grass_model<-lmer(log(richness)~yrsince_fire+Year+(1|Transect),
                             data=bird_rich_diver_grass)
anova(bird_rich_grass_model)
summary(bird_rich_grass_model)
qqnorm(resid(bird_rich_grass_model))
check_model(bird_rich_grass_model)
yrsin_grass_rich<- glht(bird_rich_grass_model, linfct=mcp(yrsince_fire="Tukey"))
print(summary(yrsin_grass_rich))

#visuals across years
bird_rich_diver_grass_yr<-bird_rich_diver_grass%>%
  group_by(yrsince_fire,Year)%>%
  summarise(rich=mean(richness, na.rm=T),
            se_rich=SE_function(richness),
            even=mean(Evar, na.rm=T),
            se_evar=SE_function(Evar))%>%
  group_by(yrsince_fire)%>%
  summarise(richness=mean(rich,na.rm=T),
            se_richness=mean(se_rich,na.rm=T),
            evenness=mean(even,na.rm=T),
            se_evenness=mean(se_evar,na.rm=T))

ggplot(bird_rich_diver_grass_yr, aes(yrsince_fire, richness, fill=yrsince_fire))+
  geom_col(width=.5) +
  geom_errorbar(aes(ymax=richness+se_richness, ymin=richness-se_richness),width=0.1)+
  scale_fill_manual(values=c("#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab("grassland bird richness")


#evenness
bird_evar_grass_model<-lmer(Evar~yrsince_fire+Year+(1|Transect),
                            data=bird_rich_diver_grass)
anova(bird_evar_grass_model)
summary(bird_evar_grass_model)
check_model(bird_evar_grass_model)

#visual
ggplot(bird_rich_diver_grass_yr, aes(yrsince_fire, evenness, fill=yrsince_fire))+
  geom_col(width=.5) +
  geom_errorbar(aes(ymax=evenness+se_evenness, ymin=evenness-se_evenness),width=0.1)+
  scale_fill_manual(values=c("#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab("grassland bird evenness")

#non-grassland bird####
yrsincef_non_grass_bird_anlys<-Yrs_since_fire_bird%>%
  filter(Grassland=="FALSE")%>%
  group_by(Year, yrsince_fire, Watershed, TransectName, Transect,SpeciesCode)%>%
  summarise(total_max=max(total, na.rm=T))%>%
  group_by(Year, yrsince_fire, Watershed, TransectName, Transect)%>%
  summarise(tot_count=sum(total_max, na.rm=T))
yrsincef_non_grass_bird_anlys$Year<-as.factor(yrsincef_non_grass_bird_anlys$Year)

#visual
yrsincef_non_grass_bird_viz<-yrsincef_non_grass_bird_anlys%>%
  group_by(Year,yrsince_fire)%>%
  summarise(total=mean(tot_count, na.rm=T),
            se_total=SE_function(tot_count))

#total grassland bird count 
ggplot(yrsincef_non_grass_bird_viz, aes(Year, total, col = yrsince_fire,
                                    fill=yrsince_fire, linetype=yrsince_fire)) +
  geom_point(size=2) +
  geom_path(aes(as.numeric(Year))) +
  geom_errorbar(aes(ymin = total - se_total, 
                    ymax = total + se_total),
                width=0.1, col="black", linetype=1) +
  scale_colour_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2")) +
  ylab("non-grassland bird count")

#model
yrsincef_non_grass_bird_model<-lmer(log(tot_count)~yrsince_fire+Year+(1|Transect), data=yrsincef_non_grass_bird_anlys)
anova(yrsincef_non_grass_bird_model)
check_model(yrsincef_non_grass_bird_model)
qqnorm(resid(yrsincef_non_grass_bird_model))
summary(yrsincef_non_grass_bird_model)
yrsincefire_non_comp<- glht(yrsincef_non_grass_bird_model, linfct=mcp(yrsince_fire="Tukey"))
print(summary(yrsincefire_non_comp))

#using linear model
yrsincef_non_grass_bird_model_lm<-lm(log(tot_count)~yrsince_fire+Year, data=yrsincef_non_grass_bird_anlys)
anova(yrsincef_non_grass_bird_model_lm)
check_model(yrsincef_non_grass_bird_model_lm)
#lmer and lm are the same

#visuals across years
yrsincef_non_grass_bird_yr<-yrsincef_non_grass_bird_anlys%>%
  group_by(yrsince_fire, Year)%>%
  summarise(total=mean(tot_count, na.rm=T),
            se_total=SE_function(tot_count))%>%
  group_by(yrsince_fire)%>%
  summarise(total=mean(total, na.rm=T),
            se_total=mean(se_total,na.rm=T))

ggplot(yrsincef_non_grass_bird_yr, aes(yrsince_fire, total, fill=yrsince_fire))+
  geom_col(width=.5) +
  geom_errorbar(aes(ymax=total+se_total, ymin=total-se_total),width=0.1)+
  scale_fill_manual(values=c("#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab("non-grassland Bird count")


#community composition 
bird_com_comp<-yrsin_bird_comm%>%
  ungroup()%>%
  select(-TransectName, -total_max, -rep_id)%>%
  filter(!Year%in%2011:2015)%>%
  pivot_wider(names_from = SpeciesCode, values_from = abund, values_fill = 0)

#separate species abundance from environmental variables
bird_comp_sp<-bird_com_comp%>%
  select(-1:-4)
bird_comp_env<-bird_com_comp%>%
  select(1:4)

#compare compositional difference
library(pairwiseAdonis)
comp_vegdist<-vegdist(bird_comp_sp)
permanova_bird_comp<-adonis(comp_vegdist~bird_comp_env$yrsince_fire+bird_comp_env$Watershed+as.factor(bird_comp_env$Year))
permanova_bird<-permanova_bird_comp$aov.tab
pairwise_perm<-pairwise.adonis2(bird_comp_sp~yrsince_fire+Watershed+as.factor(Year), data=bird_comp_env)

write.csv(pairwise_perm,"C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/pairwise_bird_comp.csv")

#visual with NMDS
#get nmds1 and 2
mds_all <- metaMDS(bird_comp_sp, distance = "bray")
#stress 0.21
#combine NMDS1 and 2 with factor columns and create centroids
nmds_sp_scores <- data.frame(bird_comp_env, scores(mds_all, display="sites"))%>%
  group_by(Year,yrsince_fire,Watershed)%>%
  mutate(NMDS1_mean=mean(NMDS1),
         NMDS2_mean=mean(NMDS2))

#plotting centroid through time
ggplot(nmds_sp_scores, aes(x=NMDS1_mean, y=NMDS2_mean, col=yrsince_fire))+
  geom_point(size=8)+
  scale_shape_manual(values=c(15:18,0:2,5))+
  scale_colour_manual(values=c("#F0E442", "#994F00", "#999999", "#0072B2"))
#community metrics with all species####
fire_bird_comm<-Yrs_since_fire_bird%>%
  group_by(Year, FireGrzTrt, Watershed, Transect,SpeciesCode)%>%
  summarise(total_max=max(total, na.rm=T))%>%
  group_by(Year, FireGrzTrt, Watershed, Transect)%>%
  mutate(abund=(total_max/sum(total_max))*100)%>%
  select(-total_max)%>%
  pivot_wider(names_from = SpeciesCode, values_from = abund, values_fill = 0)

#split intop abundance and enironmental variables
fire_sp<-fire_bird_comm%>%
  ungroup()%>%
  select(-1:-4)
fire_env<-fire_bird_comm%>%
  select(1:4)

#calculating permanova and betadiversity
#creating a loop to do this
year_vec_bird <- unique(fire_env$Year)
bird_perm <- {}
bird_beta <- {}


for(YEAR in 1:length(year_vec_bird)){
  vdist_temp_bird <- vegdist(filter(fire_sp, fire_env$Year ==  year_vec_bird[YEAR]))
  permanova_temp_bird <- adonis(vdist_temp_bird ~ subset(fire_env, Year == year_vec_bird[YEAR])$FireGrzTrt)
  permanova_out_temp_bird <- data.frame(Year = year_vec_bird[YEAR], 
                                      DF = permanova_temp_bird$aov.tab[1,1],
                                      F_value = permanova_temp_bird$aov.tab[1,4],
                                      P_value = permanova_temp_bird$aov.tab[1,6])
  bird_perm <- rbind(bird_perm,permanova_out_temp_bird)
  
  bdisp_temp_bird <- betadisper(vdist_temp_bird, filter(fire_env, Year==year_vec_bird[YEAR])$FireGrzTrt, type = "centroid", bias.adjust = TRUE)
  bdisp_out_temp_bird <- data.frame(filter(fire_env, Year==year_vec_bird[YEAR]), distance = bdisp_temp_bird$distances)
  bird_beta <- rbind(bird_beta, bdisp_out_temp_bird)
  
  rm(vdist_temp_bird, permanova_temp_bird, permanova_out_temp_bird, bdisp_temp_bird, bdisp_out_temp_bird)
}
 
#view bet diversity
beta_diver<-bird_beta%>%
  group_by(Year, FireGrzTrt)%>%
  summarise(beta=mean(distance, na.rm=T),
            se=SE_function(distance))%>%
  group_by(FireGrzTrt)%>%
  summarise(beta_a=mean(beta, na.rm=T),
            se_a=mean(se,na.rm=T))
#plot
ggplot(beta_diver,aes(x=FireGrzTrt,fill=FireGrzTrt))+
  geom_bar(stat = "identity",aes(y=beta_a),width = 0.5)+
  geom_errorbar(aes(ymin=beta_a-se_a,
                    ymax=beta_a+se_a),width=0.05,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"))
