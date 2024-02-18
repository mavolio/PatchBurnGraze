#PBG Bird count and community metrics
#Authors: Joshua Ajowele
#Started: 26 May 2022 last modified: 16 Feb 2024

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
             plot.title = element_text(size=24, vjust=2),
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

#There are multiple observation by multiple observers 2011-2016
#or just multiple observation by one observer 2011, C3C, C1A
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
  mutate(plot_index=1:length(Year))

num_bootstrap_bird<-21
bootstrap_vector<-1:num_bootstrap_bird
PBG_bird_count_master<-{}
for(BOOT in 1:length(bootstrap_vector)){
  bird_count_rand_key<-PBG_bird_boot_index%>%
    dplyr::select(Year, Watershed, FireGrzTrt, Transect, plot_index)%>%
    unique(.)%>%
    group_by(Year)%>%
    sample_n(2, replace=T)%>%
    dplyr::select(plot_index,Year)%>%
    ungroup()
  bird_count_ready<-PBG_bird_boot_index%>%
    right_join(bird_count_rand_key, by= c("Year", "plot_index"),
               multiple="all")%>%
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
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=Stab_ABGNth), linetype=2,size=1)+
  xlab("Stability")
#yet to select the same transects for stability
#temp mean
ggplot(Combo_bird_count,aes(temp_PBG_m))+
  geom_density(size=1)+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=temp_ABG_m), linetype=2,size=1)+
  xlab("temporal mean")
#temp sd
ggplot(Combo_bird_count,aes(temp_PBG_sd))+
  geom_density(size=1)+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=temp_ABG_sd), linetype=2,size=1)+
  xlab("temporal sd")

ggplot(Combo_bird_count,aes(count_PBGNth))+
  geom_density(size=.5)+
  facet_wrap(~Year, scales = "free")+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=count_ABGNth), linetype=2,size=.5)

#create visual using point 
#total count at unit scale
combo_bird_count_geompoint<-Combo_bird_count%>%
  pivot_longer(c(count_PBGNth_m,count_ABGNth),
               names_to = "treatment", values_to = "count_mean")%>%
  select(count_mean,treatment, Year, count_PBGNth_sd)%>%
  unique()%>%
  #round to two decimal places
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
  #round to two decimal places
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
