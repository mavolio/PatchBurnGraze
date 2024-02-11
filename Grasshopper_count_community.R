#BPBG grasshopper count and community metrics
#Authors: Joshua Ajowele
#Started: 26 May 2022 last modified: 09 Feb 2024

#load library
library(tidyverse)
library(vegan)
library(codyn)
library(stringr)
library(ggthemes)
library(readr)
library(performance)
library(car)
library(lme4)
library(see)
library(patchwork)
library(phia)



### Standard Error function
SE_function<-function(x,na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}

## Set graphing parameters
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=20), legend.text=element_text(size=20))
### Read in raw data
grasshopperspcomp_df_raw <- read.csv("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/PhD Wyoming_One Drive/PHD Wyoming/Thesis/PBG synthesis/Grasshopper_comp_PBG_2010-2021_type3.csv")
#turn all first letter of the genus to uppercase
grasshopperspcomp_df_raw$Genus<-str_to_title(grasshopperspcomp_df_raw$Genus)

unique(grasshopperspcomp_df_raw$Genus, order=T)
### Check existing watershed vector for grasshopper species comp
levels(factor(grasshopperspcomp_df_raw$Watershed)) # looks at all the unique entries in the watershed column
with(grasshopperspcomp_df_raw, table(Recyear, Watershed))

grasshopperspcomp_df <- grasshopperspcomp_df_raw %>% 
  mutate(Watershed = replace (Watershed, Watershed == "0c1a", "C01A"),
         Watershed = replace (Watershed, Watershed == "0c3a", "C03A"),
         Watershed = replace (Watershed, Watershed == "0c3b", "C03B"),
         Watershed = replace (Watershed, Watershed == "0c3c", "C03C"),
         Watershed = replace (Watershed, Watershed == "c01a", "C01A"),
         Watershed = replace (Watershed, Watershed == "c03a", "C03A"),
         Watershed = replace (Watershed, Watershed == "c03b", "C03B"),
         Watershed = replace (Watershed, Watershed == "c03c", "C03C"),
         Watershed = replace (Watershed, Watershed == "c1sb", "C1SB"),
         Watershed = replace (Watershed, Watershed == "c3sa", "C3SA"),
         Watershed = replace (Watershed, Watershed == "c3sb", "C3SB"),
         Watershed = replace (Watershed, Watershed == "c3sc", "C3SC"),
         Watershed = replace (Watershed, Watershed == "c1a", "C01A"),
         Watershed = replace (Watershed, Watershed == "c3a", "C03A"),
         Watershed = replace (Watershed, Watershed == "c3b", "C03B"),
         Watershed = replace (Watershed, Watershed == "c3c", "C03C"))

with(grasshopperspcomp_df, table(Recyear, Watershed))
#no entry for Repsite A watershed C03B year 2012 
#converting NAs to zero -not sure why I need this line, maybe useful in community data
#grasshopperspcomp_df[is.na(grasshopperspcomp_df)] = 0

#Create a watershed key column to merge with the raw data
watershed_key <- data.frame(Watershed = levels(factor(grasshopperspcomp_df$Watershed)),# create key to merge with raw data
                            FireGrzTrt=c("ABG", "PBG", "PBG", "PBG", "ABG", "PBG", "PBG", "PBG"))

Watershed_key2<-tibble(Watershed=levels(factor(grasshopperspcomp_df$Watershed)),# create key to merge with raw data
                       Unit=c("south", "south", "south", "south", "north", "north",
                              "north", "north"))

print(Watershed_key2)
#converting all repsite into uppercase
grasshopperspcomp_df$Repsite=toupper(grasshopperspcomp_df$Repsite)
#merging key with dataset
grassh_count_df <- grasshopperspcomp_df%>%
  full_join(watershed_key, by = "Watershed")%>%
  left_join(Watershed_key2,by="Watershed")%>%
  mutate(spe=paste(Genus,Species, sep="_"))%>%
  mutate(RecYear = Recyear)%>%
  filter(!RecYear== 2010)%>%
  #using the max cover from the two sweeps done on each transect
  group_by(Unit,RecYear,FireGrzTrt,Watershed,Repsite,spe)%>%
  summarise(Total=max(Total))%>%
  group_by(Unit,RecYear,FireGrzTrt,Watershed,Repsite)%>%
  #calculate total count per repsite
  summarise(Tcount=sum(Total, na.rm=T))%>%
  ungroup()

#convert to factors
grassh_count_df$Unit=as.factor(grassh_count_df$Unit)
grassh_count_df$Watershed=as.factor(grassh_count_df$Watershed)
grassh_count_df$Repsite=as.factor(grassh_count_df$Repsite)
grassh_count_df$FireGrzTrt=as.factor(grassh_count_df$FireGrzTrt)
grassh_count_df$RecYear=as.factor(grassh_count_df$RecYear)
#transect model for total count
grassh_count_tnst_model<-lmer(log(Tcount)~FireGrzTrt*RecYear+(1|Unit/Watershed),
                              data=grassh_count_df)#issingular
check_model(grassh_count_tnst_model)
qqnorm(residuals(grassh_count_tnst_model))
Anova(grassh_count_tnst_model, type=3)

testInteractions(grassh_count_tnst_model, fixed="RecYear",
                 across = "FireGrzTrt", adjustment="BH")
#using mean estimate from post-hoc to create figure as a comparison to raw data
model_estimates_tnst<-interactionMeans(grassh_count_tnst_model)
#replacing spaces in column names with underscore 
names(model_estimates_tnst)<-str_replace_all(names(model_estimates_tnst), " ","_")
#df for visuals from model estimates
model_estimates_tnst_viz<-model_estimates_tnst%>%
  mutate(count_tnst_bt_mean=exp(adjusted_mean),
         count_tnst_bt_upper=exp(adjusted_mean+SE_of_link),
         count_tnst_bt_lower=exp(adjusted_mean-SE_of_link))
#visual
ggplot(model_estimates_tnst_viz,aes(RecYear, count_tnst_bt_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=count_tnst_bt_lower,
                    ymax=count_tnst_bt_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

#watershed Scale
grassh_count_wshed_df<-grassh_count_df%>%
  group_by(RecYear, Unit, FireGrzTrt,Watershed)%>%
  summarise(Tcount_avg=mean(Tcount, na.rm=T),
            sd_Tcount=sd(Tcount),
            cv_Tcount=sd_Tcount/Tcount_avg)%>%
  ungroup()


#watershed models
grassh_count_wshed_model<-lmer(log(Tcount_avg)~FireGrzTrt*RecYear+(1|Unit),
                              data=grassh_count_wshed_df) #issingular
check_model(grassh_count_wshed_model)
qqnorm(residuals(grassh_count_wshed_model))
Anova(grassh_count_wshed_model, type=3)

testInteractions(grassh_count_wshed_model, fixed="RecYear",
                 across = "FireGrzTrt", adjustment="BH")
#using mean estimate from post-hoc to create figure as a comparison to raw data
model_estimates_wshed<-interactionMeans(grassh_count_wshed_model)
#replacing spaces in column names with underscore 
names(model_estimates_wshed)<-str_replace_all(names(model_estimates_wshed), " ","_")
#df for visuals from model estimates
model_estimates_wshed_viz<-model_estimates_wshed%>%
  mutate(count_wshed_bt_mean=exp(adjusted_mean),
         count_wshed_bt_upper=exp(adjusted_mean+SE_of_link),
         count_twshed_bt_lower=exp(adjusted_mean-SE_of_link))
#visual
ggplot(model_estimates_wshed_viz,aes(RecYear, count_wshed_bt_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=count_twshed_bt_lower,
                    ymax=count_wshed_bt_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))


#watershed sd models
grassh_count_wshed_sd_model<-lmer(log(sd_Tcount)~FireGrzTrt*RecYear+(1|Unit),
                               data=grassh_count_wshed_df) #issingular
check_model(grassh_count_wshed_sd_model)
qqnorm(residuals(grassh_count_wshed_sd_model))
Anova(grassh_count_wshed_sd_model, type=3)

testInteractions(grassh_count_wshed_sd_model, fixed="RecYear",
                 across = "FireGrzTrt", adjustment="BH")
#using mean estimate from post-hoc to create figure as a comparison to raw data
model_estimates_sd_wshed<-interactionMeans(grassh_count_wshed_sd_model)
#replacing spaces in column names with underscore 
names(model_estimates_sd_wshed)<-str_replace_all(names(model_estimates_sd_wshed), " ","_")
#df for visuals from model estimates
model_estimates_wshed_sd_viz<-model_estimates_sd_wshed%>%
  mutate(count_wshed_sd_bt_mean=exp(adjusted_mean),
         count_wshed_sd_bt_upper=exp(adjusted_mean+SE_of_link),
         count_twshed_sd_bt_lower=exp(adjusted_mean-SE_of_link))
#visual
ggplot(model_estimates_wshed_sd_viz,aes(RecYear, count_wshed_sd_bt_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=count_twshed_sd_bt_lower,
                    ymax=count_wshed_sd_bt_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

#watershed cv models
grassh_count_wshed_cv_model<-lmer(log(cv_Tcount)~FireGrzTrt+RecYear+(1|Unit),
                                  data=grassh_count_wshed_df) #issingular
check_model(grassh_count_wshed_cv_model)
qqnorm(residuals(grassh_count_wshed_cv_model))
Anova(grassh_count_wshed_cv_model, type=3)

testInteractions(grassh_count_wshed_cv_model, fixed="RecYear",
                 across = "FireGrzTrt", adjustment="BH")
#using mean estimate from post-hoc to create figure as a comparison to raw data
model_estimates_cv_wshed<-interactionMeans(grassh_count_wshed_cv_model)
#replacing spaces in column names with underscore 
names(model_estimates_cv_wshed)<-str_replace_all(names(model_estimates_cv_wshed), " ","_")
#df for visuals from model estimates
model_estimates_wshed_cv_viz<-model_estimates_cv_wshed%>%
  mutate(count_wshed_cv_bt_mean=exp(adjusted_mean),
         count_wshed_cv_bt_upper=exp(adjusted_mean+SE_of_link),
         count_twshed_cv_bt_lower=exp(adjusted_mean-SE_of_link))
#visual
ggplot(model_estimates_wshed_cv_viz,aes(RecYear, count_wshed_cv_bt_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=count_twshed_cv_bt_lower,
                    ymax=count_wshed_cv_bt_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

#stability at watershed scale
stab_wshed_df<-grassh_count_wshed_df%>%
  group_by(Unit,Watershed,FireGrzTrt)%>%
  summarise(stab=mean(Tcount_avg,na.rm=T)/sd(Tcount_avg))%>%
  ungroup()

#stability model
stab_model_wshed<-lmer(stab~FireGrzTrt+(1|Unit),
                       data=stab_wshed_df)#issingular

check_model(stab_model_wshed)
qqnorm(residuals(stab_model_wshed))
Anova(stab_model_wshed, type=3)

#using mean estimate to create figure 
model_estimates_stab_wshed<-interactionMeans(stab_model_wshed)
#replacing spaces in column names with underscore 
names(model_estimates_stab_wshed)<-str_replace_all(names(model_estimates_stab_wshed), " ","_")
#df for visuals from model estimates
model_estimates_wshed_stab_viz<-model_estimates_stab_wshed%>%
  mutate(count_wshed_stab_bt_mean=exp(adjusted_mean),
         count_wshed_stab_bt_upper=exp(adjusted_mean+SE_of_link),
         count_twshed_stab_bt_lower=exp(adjusted_mean-SE_of_link))
#visual
ggplot(model_estimates_wshed_stab_viz,aes(x=FireGrzTrt,fill=FireGrzTrt))+
  geom_bar(stat = "identity",aes(y=count_wshed_stab_bt_mean),width = 0.5)+
  geom_errorbar(aes(ymin=count_twshed_stab_bt_lower,
                    ymax=count_wshed_stab_bt_upper),width=0.2,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"))

##unit scale count
####bootstrap
#extract PBG and separate into Unit
PBG_south_count<-grassh_count_df%>%
  filter(FireGrzTrt=="PBG" & Unit=="south")%>%
  ungroup()
PBG_north_count<-grassh_count_df%>%
  filter(FireGrzTrt=="PBG" & Unit=="north")%>%
  ungroup()

#create an index 
count_south_index<-PBG_south_count%>%
  group_by(RecYear)%>%
  mutate(plot_index=1:length(RecYear))
count_north_index<-PBG_north_count%>%
  group_by(RecYear)%>%
  mutate(plot_index=1:length(RecYear))

num_bootstrap<-1000
bootstrap_vector<-1:num_bootstrap
PBG_south_count_master<-{}
for(BOOT in 1:length(bootstrap_vector)){
  south_count_rand_key<-count_south_index%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Repsite, plot_index)%>%
    unique(.)%>%
    group_by(RecYear)%>%
    sample_n(4, replace=T)%>%
    dplyr::select(plot_index,RecYear)%>%
    ungroup()
  count_south_ready<-count_south_index%>%
    right_join(south_count_rand_key, by= c("RecYear", "plot_index"),
               multiple="all")%>%
    mutate(iteration=BOOT)
  PBG_south_count_master<-rbind(PBG_south_count_master, count_south_ready)
}
write.csv(PBG_south_biomass_master, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/PBG_south_biomass.csv")







#north unit
num_bootstrap<-1000
bootstrap_vector<-1:num_bootstrap
PBG_north_count_master<-{}
for(BOOT in 1:length(bootstrap_vector)){
  north_count_rand_key<-count_north_index%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Repsite, plot_index)%>%
    unique(.)%>%
    group_by(RecYear)%>%
    sample_n(4, replace=T)%>%
    dplyr::select(plot_index,RecYear)%>%
    ungroup()
  count_north_ready<-count_north_index%>%
    right_join(north_count_rand_key, by= c("RecYear", "plot_index"),
               multiple="all")%>%
    mutate(iteration=BOOT)
  PBG_north_count_master<-rbind(PBG_north_count_master, count_north_ready)
}

#extract ABG
#extract ABG and separate into Unit
ABG_south_count<-grass_count_df%>%
  filter(FireGrzTrt=="ABG" & Unit=="south")%>%
  group_by(RecYear)%>%#calculate metrics at the unit level
  summarise(count_ABGSth=mean(Tcount, na.rm=T),
            count_ABGSth_sd=sd(Tcount),
            cv_count_ABGSth=sd(Tcount)/mean(Tcount, na.rm=T))%>%
  ungroup()%>%
  mutate(Stab_ABGSth=mean(count_ABGSth, na.rm =T)/sd(count_ABGSth))
ABG_south_count$RecYear<-as.factor(ABG_south_count$RecYear)
ABG_north_count<-grass_count_df%>%
  filter(FireGrzTrt=="ABG" & Unit=="north")%>%
  group_by(RecYear)%>%#calculate metrics at the unit level
  summarise(count_ABGNth=mean(Tcount, na.rm=T),
            count_ABGNth_sd=sd(Tcount),
            cv_count_ABGNth=sd(Tcount)/mean(Tcount, na.rm=T))%>%
  ungroup()%>%
  mutate(Stab_ABGNth=mean(count_ABGNth, na.rm =T)/sd(count_ABGNth))
ABG_north_count$RecYear<-as.factor(ABG_north_count$RecYear)
#calculate variables for PBG bootstrap
PBG_north_mean<-PBG_north_count_master%>%
  group_by(RecYear,iteration)%>%
  summarise(count_PBGNth=mean(Tcount,na.rm=T),
            sd_c_PBGNth=sd(Tcount),
            cv_c_PBGNth=sd(Tcount)/mean(Tcount, na.rm=T))%>%
  group_by(iteration)%>%
  mutate(stab_PBGNth=mean(count_PBGNth, na.rm=T)/sd(count_PBGNth))

PBG_south_mean<-PBG_south_count_master%>%
  group_by(RecYear,iteration)%>%
  summarise(count_PBGSth=mean(Tcount,na.rm=T),
            sd_c_PBGSth=sd(Tcount),
            cv_c_PBGSth=sd(Tcount)/mean(Tcount, na.rm=T))%>%
  group_by(iteration)%>%
  mutate(stab_PBGSth=mean(count_PBGSth, na.rm=T)/sd(count_PBGSth))


#combine PBG and ABG
Combo_north_count<-PBG_north_mean%>%
  left_join(ABG_north_count, by="RecYear")%>%
  group_by(RecYear)%>%
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
write.csv(Combo_north_count, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/grassH_Ncount.csv")

Combo_south_count<-PBG_south_mean%>%
  left_join(ABG_south_count, by="RecYear")%>%
  group_by(RecYear)%>%
  mutate(count_PBGSth_m=mean(count_PBGSth, na.rm=T),
         count_PBGSth_sd=sd(count_PBGSth),
         z_score_SthMean=((count_ABGSth-count_PBGSth_m)/count_PBGSth_sd),
         p_value_SMean=2*pnorm(-abs(z_score_SthMean)),
         count_PBGSth_sd_M=mean(sd_c_PBGSth, na.rm=T),
         count_PBGSth_sd_sd=sd(sd_c_PBGSth),
         Z_score_Ssd=((count_ABGSth_sd-count_PBGSth_sd_M)/count_PBGSth_sd_sd),
         pvalue_Ssd=2*pnorm(-abs(Z_score_Ssd)),
         count_PBGSth_cv_M=mean(cv_c_PBGSth, na.rm=T),
         count_PBGSth_cv_sd=sd(cv_c_PBGSth),
         Z_score_Scv=((cv_count_ABGSth-count_PBGSth_cv_M)/count_PBGSth_cv_sd),
         pvalue_Scv=2*pnorm(-abs(Z_score_Scv)))%>%
  ungroup()%>%
  mutate(stab_PBGSm=mean(stab_PBGSth,na.rm=T),
         stab_PBGS_sd=sd(stab_PBGSth),
         z_score_SStab=((Stab_ABGSth-stab_PBGSm)/stab_PBGS_sd),
         pvalue_stab=2*pnorm(-abs(z_score_SStab)))
write.csv(Combo_south_count, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/grassH_Scount.csv")


#create a visual of the distribution
ggplot(Combo_north_biomass,aes(biomass_PBGNth))+
  geom_density(size=.5)+
  facet_wrap(~RecYear, scales = "free")+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=biomass_ABGNth), linetype=2,size=.5)
ggplot(Combo_south_biomass,aes(biomass_PBGSth))+
  geom_density(size=.5)+
  facet_wrap(~RecYear, scales = "free")+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=biomass_ABGSth), linetype=2,size=.5)



###Unit Scale
grassh_comm_df <- grasshopperspcomp_df%>%
  full_join(watershed_key, by = "Watershed")%>%
  left_join(Watershed_key2,by="Watershed")%>%
  mutate(spe=paste(Genus,Species, sep="_"))%>%
  mutate(RecYear = Recyear)%>%
  filter(!RecYear== 2010)%>%
  #using the max cover from the two sweeps done on each transect
  group_by(Unit,RecYear,FireGrzTrt,Watershed,Repsite,spe)%>%
  summarise(Total=max(Total))%>%
  group_by(Unit,RecYear,FireGrzTrt,Watershed,Repsite)%>%
  #converting count data to abundance data (0-100%)
  mutate(abundance=(Total/sum(Total, na.rm=T))*100)%>%
  group_by(Unit,RecYear,FireGrzTrt,Watershed,Repsite,spe)%>%
  #need to summarise to have species in same transect on the same row when pivot wider is used
  summarise(abundance=mean(abundance,na.rm=T))%>%
  ungroup()%>%
  group_by(Unit,RecYear,FireGrzTrt,Watershed,Repsite)%>%
  mutate(yr_trt_unit=paste(Unit,RecYear,FireGrzTrt,Watershed,Repsite, sep="_"))

#calculate community metrics
richness_evar<-grassh_comm_df%>%
  community_structure(abundance.var = "abundance", replicate.var = "yr_trt_unit",
                      metric="Evar")
diversity<-grassh_comm_df%>%
  community_diversity(abundance.var = "abundance",replicate.var = "yr_trt_unit",
                      metric="Shannon")

#combine community metrics with larger df
grassh_comm_metrics_df<-grassh_comm_df%>%
  select(-c(spe,abundance))%>%
  unique()%>%
  left_join(richness_evar, by="yr_trt_unit")%>%
  left_join(diversity, by="yr_trt_unit")



####bootstrap
#extract PBG and separate into Unit
PBG_south_biomass<-biomass_data%>%
  filter(FireGrzTrt=="PBG" & Unit=="south")%>%
  group_by(RecYear, Unit, Watershed, Transect, Plotnum, FireGrzTrt)%>%
  summarise(biomass=mean(biomass, na.rm=T))%>%
  ungroup()
PBG_north_biomass<-biomass_data%>%
  filter(FireGrzTrt=="PBG" & Unit=="north")%>%
  group_by(RecYear, Unit, Watershed, Transect, Plotnum, FireGrzTrt)%>%
  summarise(biomass=mean(biomass, na.rm=T))%>%
  ungroup()

#create an index for the plots
biomass_south_index<-PBG_south_biomass%>%
  group_by(RecYear)%>%
  mutate(plot_index=1:length(RecYear))

num_bootstrap<-1000
bootstrap_vector<-1:num_bootstrap
PBG_south_biomass_master<-{}
for(BOOT in 1:length(bootstrap_vector)){
  south_rand_key<-biomass_south_index%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Transect, Plotnum, plot_index)%>%
    unique(.)%>%
    #filter(RecYear==2012)%>%
    group_by(RecYear)%>%
    sample_n(200, replace=T)%>%
    dplyr::select(plot_index)%>%
    ungroup()
  biomass_south_ready<-biomass_south_index%>%
    right_join(south_rand_key, by= c("RecYear", "plot_index"),
               multiple="all")%>%
    mutate(iteration=BOOT)
  PBG_south_biomass_master<-rbind(PBG_south_biomass_master, biomass_south_ready)
}
write.csv(PBG_south_biomass_master, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/PBG_south_biomass.csv")


