#PBG Bird count and community metrics
#Authors: Joshua Ajowele
#Started: 26 May 2022 last modified: 05 Nov 2024

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
#check data for entry errors, 2018 C3A seems to be suspect: 
#I think it has an extra transect that should belong to C3B-90% sure-Corrected in raw data
#email Alice about it

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
yrsincef_bird_model<-lmer(log(tot_max)~yrsince_fire*Year+(1|Transect), data=yrsincef_bird_anlys)
check_model(yrsincef_bird_model)
anova(yrsincef_bird_model, type=3)
summary(yrsincef_bird_model)
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
yrsincefire_comp<- glht(yrsincef_bird_model, linfct=mcp(yrsince_fire="Tukey"))
print(summary(yrsincefire_comp))
#ALL P-ADJUSTED 


#visuals
yrsincef_bird_viz_yr<-yrsincef_bird_anlys%>%
  group_by(yrsince_fire)%>%
  summarise(total=mean(tot_max, na.rm=T),
            se_total=SE_function(tot_max))

ggplot(yrsincef_bird_viz_yr, aes(yrsince_fire, total, fill=yrsince_fire))+
  geom_col(width=.5) +
  geom_errorbar(aes(ymax=total+se_total, ymin=total-se_total),width=0.1)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"))+
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
  select(Year, yrsince_fire, Watershed, TransectName, Transect,rep_id) %>%
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
bird_diver_model<-lmer(Shannon~yrsince_fire*Year+(1|Transect),
                     data=bird_rich_diver)
anova(bird_diver_model)
summary(bird_diver_model)
qqnorm(resid(bird_diver_model))
check_model(bird_diver_model)

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
  group_by(yrsince_fire)%>%
  summarise(rich=mean(richness,na.rm=T),
            rich_se=SE_function(richness),
            diver=mean(Shannon,na.rm=T),
            diver_se=SE_function(Shannon),
            evar=mean(Evar,na.rm=T),
            evar_se=SE_function(Evar))
 #geompoint 
ggplot(bird_rich_diver_viz, aes(Year, diver, col = yrsince_fire,
                              fill=yrsince_fire, linetype=yrsince_fire)) +
  geom_point(size=2, col="black") +
  geom_path(aes(as.numeric(Year))) +
  geom_errorbar(aes(ymin = diver - diver_se, 
                    ymax = diver + diver_se),
                width=0.1, col="black", linetype=1) +
  scale_colour_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab("bird diversity")
#bargraph
ggplot(bird_rich_diver_bar, aes(yrsince_fire, diver, fill=yrsince_fire))+
  geom_col(width=.5) +
  geom_errorbar(aes(ymax=diver + diver_se, ymin=diver - diver_se),width=0.1)+
  scale_fill_manual(values=c("#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab("Bird diversity")

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
  group_by(yrsince_fire)%>%
  summarise(total=mean(tot_count, na.rm=T),
            se_total=SE_function(tot_count))

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
bird_diver_grass_model<-lmer(Shannon~yrsince_fire*Year+(1|Transect),
                       data=bird_rich_diver_grass)
anova(bird_diver_grass_model)
summary(bird_diver_grass_model)
qqnorm(resid(bird_diver_grass_model))
check_model(bird_diver_grass_model)
#MULTICOMPARISON
#yrsin_grass_diver<- glht(bird_diver_grass_model, linfct=mcp(yrsince_fire="Tukey"))
#print(summary(yrsin_grass_diver))
#richness
bird_rich_grass_model<-lmer(log(richness)~yrsince_fire*Year+(1|Transect),
                             data=bird_rich_diver_grass)
anova(bird_rich_grass_model)
summary(bird_rich_grass_model)
qqnorm(resid(bird_rich_grass_model))
check_model(bird_rich_grass_model)
#yrsin_grass_rich<- glht(bird_rich_grass_model, linfct=mcp(yrsince_fire="Tukey"))
#print(summary(yrsin_grass_rich))



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
  group_by(yrsince_fire)%>%
  summarise(total=mean(tot_count, na.rm=T),
            se_total=SE_function(tot_count))

ggplot(yrsincef_non_grass_bird_yr, aes(yrsince_fire, total, fill=yrsince_fire))+
  geom_col(width=.5) +
  geom_errorbar(aes(ymax=total+se_total, ymin=total-se_total),width=0.1)+
  scale_fill_manual(values=c("#F0E442", "#994F00", "#999999", "#0072B2"))+
  ylab("non-grassland Bird count")










#data ready for visuals
PBG_bird_viz_ready_all<-PBG_bird_viz%>%
  rename(total=Count)%>%
  #remove NA otherwise it will appear in max count
  filter(total!="NA")%>%
  group_by(Year, FireGrzTrt, Watershed, TransectName, Transect, SpeciesCode)%>%
  summarise(total_max=max(total))

#Separate PBG from ABG to perform bootstrap####
PBG_bird_boot<-PBG_bird_viz_ready_all%>%
  filter(FireGrzTrt=="PBG")%>%
  group_by(Year, FireGrzTrt, Watershed, TransectName, Transect)%>%
  summarise(total_count=sum(total_max))%>%
  ungroup()
ABG_bird_boot<-PBG_bird_viz_ready_all%>%
  filter(FireGrzTrt=="ABG")%>%
  group_by(Year, FireGrzTrt, Watershed, TransectName, Transect)%>%
  summarise(total_count=sum(total_max))%>%
  ungroup()

#Bootstrapping for total count####
#create an index 
PBG_bird_boot_index<-PBG_bird_boot%>%
  group_by(Year)%>%
  mutate(birdot_index=1:length(Year))

num_bootstrap_bird<-21
bootstrap_vector<-1:num_bootstrap_bird
PBG_bird_count_master<-{}
for(BOOT in 1:length(bootstrap_vector)){
  bird_count_rand_key<-PBG_bird_boot_index%>%
    dbirdyr::select(Year, Watershed, FireGrzTrt, Transect, birdot_index)%>%
    unique(.)%>%
    group_by(Year)%>%
    sambirde_n(2, rebirdace=T)%>%
    dbirdyr::select(birdot_index,Year)%>%
    ungroup()
  bird_count_ready<-PBG_bird_boot_index%>%
    right_join(bird_count_rand_key, by= c("Year", "birdot_index"),
               multibirde="all")%>%
    mutate(iteration=BOOT)
  PBG_bird_count_master<-rbind(PBG_bird_count_master, bird_count_ready)
}
write.csv(PBG_bird_count_master, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/PBG_bird_count_master.csv")

#calculate count sd,cv,mean and stability####
ABG_bird_count<-ABG_bird_boot%>%
  group_by(Year)%>%#calculate metrics at the unit level
  summarise(count_ABGNth=mean(total_count, na.rm=T),
            count_ABGNth_sd=sd(total_count),
            cv_count_ABGNth=sd(total_count)/mean(total_count, na.rm=T))%>%
  ungroup()%>%
  mutate(Stab_ABGNth=mean(count_ABGNth, na.rm =T)/sd(count_ABGNth),
         temp_ABG_m=mean(count_ABGNth,na.rm=T),
         temp_ABG_sd=sd(count_ABGNth))

#calculate variables for PBG bootstrap
PBG_bird_count<-PBG_bird_count_master%>%
  group_by(Year,iteration)%>%
  summarise(count_PBGNth=mean(total_count,na.rm=T),
            sd_c_PBGNth=sd(total_count),
            cv_c_PBGNth=sd(total_count)/mean(total_count, na.rm=T))%>%
  group_by(iteration)%>%
  mutate(stab_PBGNth=mean(count_PBGNth, na.rm=T)/sd(count_PBGNth),
         temp_PBG_m=mean(count_PBGNth, na.rm=T),
         temp_PBG_sd=sd(count_PBGNth))

#combine PBG and ABG and calculate zscores####
Combo_bird_count<-PBG_bird_count%>%
  left_join(ABG_bird_count, by="Year")%>%
  group_by(Year)%>%
  mutate(count_PBGNth_m=mean(count_PBGNth, na.rm=T),
         count_PBGNth_sd=sd(count_PBGNth),
         z_score_NthMean=((count_ABGNth-count_PBGNth_m)/count_PBGNth_sd),
         p_value_NMean=2*pnorm(-abs(z_score_NthMean)),
         count_PBGNth_sd_M=mean(sd_c_PBGNth, na.rm=T),
         count_PBGNth_sd_sd=sd(sd_c_PBGNth),
         Z_score_Nsd=((count_ABGNth_sd-count_PBGNth_sd_M)/count_PBGNth_sd_sd),
         pvalue_Nsd=2*pnorm(-abs(Z_score_Nsd)),
         count_PBGNth_cv_M=mean(cv_c_PBGNth, na.rm=T),
         count_PBGNth_cv_sd=sd(cv_c_PBGNth),
         Z_score_Ncv=((cv_count_ABGNth-count_PBGNth_cv_M)/count_PBGNth_cv_sd),
         pvalue_Ncv=2*pnorm(-abs(Z_score_Ncv)))%>%
  ungroup()%>%
  mutate(stab_PBGNm=mean(stab_PBGNth,na.rm=T),
         stab_PBGN_sd=sd(stab_PBGNth),
         z_score_NStab=((Stab_ABGNth-stab_PBGNm)/stab_PBGN_sd),
         pvalue_stab=2*pnorm(-abs(z_score_NStab)))

#create figures for bird count####
#create visual for stab
ggplot(Combo_bird_count,aes(stab_PBGNth))+
  geom_density(size=1)+
  #facet_grid("Year")+
  geom_vline(aes(xintercept=Stab_ABGNth), linetype=2,size=1)+
  xlab("Stability")
#yet to select the same transects for stability
#temp mean
ggplot(Combo_bird_count,aes(temp_PBG_m))+
  geom_density(size=1)+
  #facet_grid("Year")+
  geom_vline(aes(xintercept=temp_ABG_m), linetype=2,size=1)+
  xlab("temporal mean")
#temp sd
ggplot(Combo_bird_count,aes(temp_PBG_sd))+
  geom_density(size=1)+
  #facet_grid("Year")+
  geom_vline(aes(xintercept=temp_ABG_sd), linetype=2,size=1)+
  xlab("temporal sd")

ggplot(Combo_bird_count,aes(count_PBGNth))+
  geom_density(size=.5)+
  facet_wrap(~Year, scales = "free")+
  #facet_grid("Year")+
  geom_vline(aes(xintercept=count_ABGNth), linetype=2,size=.5)

#create visual using point 
#total count at unit scale
combo_bird_count_geompoint<-Combo_bird_count%>%
  pivot_longer(c(count_PBGNth_m,count_ABGNth),
               names_to = "treatment", values_to = "count_mean")%>%
  select(count_mean,treatment, Year, count_PBGNth_sd)%>%
  unique()%>%
  #round to two decimal birdaces
  mutate(sdsd=round(count_PBGNth_sd,digits=2))%>%
  mutate(sdsd=ifelse(sdsd==5.90 & treatment=="count_ABGNth",0,sdsd),
         sdsd=ifelse(sdsd==3.34 & treatment=="count_ABGNth",0,sdsd),
         sdsd=ifelse(sdsd==1.56 & treatment=="count_ABGNth",0,sdsd),
         sdsd=ifelse(sdsd==2.78 & treatment=="count_ABGNth",0,sdsd),
         sdsd=ifelse(sdsd==2.58 & treatment=="count_ABGNth",0,sdsd),
         sdsd=ifelse(sdsd==2.73 & treatment=="count_ABGNth",0,sdsd),
         sdsd=ifelse(sdsd==7.25 & treatment=="count_ABGNth",0,sdsd),
         sdsd=ifelse(sdsd==4.45 & treatment=="count_ABGNth",0,sdsd),
         sdsd=ifelse(sdsd==1.16 & treatment=="count_ABGNth",0,sdsd),
         sdsd=ifelse(sdsd==3.61 & treatment=="count_ABGNth",0,sdsd))

ggplot(combo_bird_count_geompoint,aes(Year, count_mean, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(Year)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=count_mean+1.96*(sdsd),
                    ymin=count_mean-1.96*(sdsd)),width=.2)+
  scale_x_continuous(breaks=2011:2021)


#sd unit scale
combo_bird_sd_geompoint<-Combo_bird_count%>%
  pivot_longer(c(count_PBGNth_sd_M,count_ABGNth_sd),
               names_to = "treatment", values_to = "count_sd")%>%
  select(count_sd,treatment, Year, count_PBGNth_sd_sd)%>%
  unique()%>%
  #round to two decimal birdaces
  mutate(sdsd=round(count_PBGNth_sd_sd,digits=2))%>%
  mutate(sdsd=ifelse(sdsd==7.27 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==2.17 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==2.16 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==2.25 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==2.18 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==2.51 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==5.99 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==3.76 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==0.80 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==3.17 & treatment=="count_ABGNth_sd",0,sdsd))

ggplot(combo_bird_sd_geompoint,aes(Year, count_sd, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(Year)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=count_sd+1.96*(sdsd),
                    ymin=count_sd-1.96*(sdsd)),width=.2)+
  scale_x_continuous(breaks=2011:2021)
