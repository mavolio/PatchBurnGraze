#BPBG grasshopper count and community metrics
#Authors: Joshua Ajowele
#Started: 26 May 2022 last modified: 16 Feb 2024

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
ggplot(Combo_north_count,aes(count_PBGNth))+
  geom_density(size=.5)+
  facet_wrap(~RecYear, scales = "free")+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=count_ABGNth), linetype=2,size=.5)
ggplot(Combo_south_count,aes(count_PBGSth))+
  geom_density(size=.5)+
  facet_wrap(~RecYear, scales = "free")+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=count_ABGSth), linetype=2,size=.5)
ggplot(Combo_north_count,aes(sd_c_PBGNth))+
  geom_density(size=.5)+
  facet_wrap(~RecYear, scales = "free")+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=count_ABGNth_sd), linetype=2,size=.5)
ggplot(Combo_south_count,aes(sd_c_PBGSth))+
  geom_density(size=.5)+
  facet_wrap(~RecYear, scales = "free")+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=count_ABGSth_sd), linetype=2,size=.5)
ggplot(Combo_north_count,aes(cv_c_PBGNth))+
  geom_density(size=.5)+
  facet_wrap(~RecYear, scales = "free")+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=cv_count_ABGNth), linetype=2,size=.5)
ggplot(Combo_south_count,aes(cv_c_PBGSth))+
  geom_density(size=.5)+
  facet_wrap(~RecYear, scales = "free")+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=cv_count_ABGSth), linetype=2,size=.5)

#
#create visual for stab
ggplot(Combo_north_count,aes(stab_PBGNth))+
  geom_density(size=1)+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=Stab_ABGNth), linetype=2,size=1)+
  xlab("Stability North Unit")
ggplot(Combo_south_count,aes(stab_PBGSth))+
  geom_density(size=1)+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=Stab_ABGSth), linetype=2,size=1)+
  xlab("Stability South Unit")

#convert density plot for mean, sd and cv count to geompoint####
#create visual using point
combo_north_count_geompoint<-Combo_north_count%>%
  pivot_longer(c(count_PBGNth_m,count_ABGNth),
               names_to = "treatment", values_to = "count_mean")
ggplot(combo_north_count_geompoint,aes(RecYear, count_mean, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=count_mean+1.96*(count_PBGNth_sd),
                    ymin=count_mean-1.96*(count_PBGNth_sd)),width=.2)

combo_north_sd_geompoint<-Combo_north_count%>%
  pivot_longer(c(count_PBGNth_sd_M,count_ABGNth_sd),
               names_to = "treatment", values_to = "count_sd")
ggplot(combo_north_sd_geompoint,aes(RecYear, count_sd, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=count_sd+1.96*(count_PBGNth_sd_sd),
                    ymin=count_sd-1.96*(count_PBGNth_sd_sd)),width=.2)

combo_north_cv_geompoint<-Combo_north_count%>%
  pivot_longer(c(count_PBGNth_cv_M,cv_count_ABGNth),
               names_to = "treatment", values_to = "count_cv")
ggplot(combo_north_cv_geompoint,aes(RecYear, count_cv, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=count_cv+1.96*(count_PBGNth_cv_sd),
                    ymin=count_cv-1.96*(count_PBGNth_cv_sd)),width=.2)
#next step: visuals for South Unit

###Calculate community metrics fror each transect####
#community metrics
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

#message Yang about South PBG transect A C3B

####bootstrap for community metrics####
#extract PBG and separate into Unit
PBG_south_GH_metrics<-grassh_comm_metrics_df%>%
  filter(FireGrzTrt=="PBG" & Unit=="south")%>%
  ungroup()
PBG_north_GH_metrics<-grassh_comm_metrics_df%>%
  filter(FireGrzTrt=="PBG" & Unit=="north")%>%
  ungroup()

#create an index for the plots
GH_metrics_south_index<-PBG_south_GH_metrics%>%
  group_by(RecYear)%>%
  mutate(plot_index=1:length(RecYear))
GH_metrics_north_index<-PBG_north_GH_metrics%>%
  group_by(RecYear)%>%
  mutate(plot_index=1:length(RecYear))

num_bootstrap<-1000
bootstrap_vector<-1:num_bootstrap
PBG_south_GH_metric_master<-{}
for(BOOT in 1:length(bootstrap_vector)){
  south_rand_metric_key<-GH_metrics_south_index%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Repsite, plot_index)%>%
    unique(.)%>%
    group_by(RecYear)%>%
    sample_n(4, replace=T)%>%
    dplyr::select(plot_index,RecYear)%>%
    ungroup()
  metrics_south_ready<-GH_metrics_south_index%>%
    right_join(south_rand_metric_key, by= c("RecYear", "plot_index"),
               multiple="all")%>%
    mutate(iteration=BOOT)
  PBG_south_GH_metric_master<-rbind(PBG_south_GH_metric_master, metrics_south_ready)
}
write.csv(PBG_south_GH_metric_master, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/PBG_south_GH_metric_master.csv")

PBG_north_GH_metric_master<-{}
for(BOOT in 1:length(bootstrap_vector)){
  north_rand_metric_key<-GH_metrics_north_index%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Repsite, plot_index)%>%
    unique(.)%>%
    group_by(RecYear)%>%
    sample_n(4, replace=T)%>%
    dplyr::select(plot_index,RecYear)%>%
    ungroup()
  metrics_north_ready<-GH_metrics_north_index%>%
    right_join(north_rand_metric_key, by= c("RecYear", "plot_index"),
               multiple="all")%>%
    mutate(iteration=BOOT)
  PBG_north_GH_metric_master<-rbind(PBG_north_GH_metric_master, metrics_north_ready)
}
write.csv(PBG_north_GH_metric_master, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/PBG_north_GH_metric_master.csv")

#calculate sd,cv and stability of community metrics####
#extract ABG and separate into Unit
ABG_south_GH_metrics<-grassh_comm_metrics_df%>%
  filter(FireGrzTrt=="ABG" & Unit=="south")%>%
  ungroup()%>%
  group_by(Unit, RecYear)%>%
  summarise(rich_m_SthABG=mean(richness, na.rm=T),
            rich_sd_SthABG=sd(richness),
            rich_cv_SthABG=sd(richness)/mean(richness,na.rm=T),
            even_m_SthABG=mean(Evar,na.rm=T),
            even_sd_SthABG=sd(Evar),
            even_cv_SthABG=sd(Evar)/mean(Evar,na.rm=T),
            shan_m_SthABG=mean(Shannon,na.rm=T),
            shan_sd_SthABG=sd(Shannon),
            shan_cv_SthABG=sd(Shannon)/mean(Shannon,na.rm=T))%>%
  group_by(Unit)%>%
  mutate(rich_stab_SthABG=mean(rich_m_SthABG,na.rm=T)/sd(rich_m_SthABG),
         even_stab_SthABG=mean(even_m_SthABG)/sd(even_m_SthABG),
         shan_stab_SthABG=mean(shan_m_SthABG,na.rm=T)/sd(shan_m_SthABG))
ABG_north_GH_metrics<-grassh_comm_metrics_df%>%
  filter(FireGrzTrt=="ABG" & Unit=="north")%>%
  ungroup()%>%
  group_by(Unit, RecYear)%>%
  summarise(rich_m_NthABG=mean(richness, na.rm=T),
            rich_sd_NthABG=sd(richness),
            rich_cv_NthABG=sd(richness)/mean(richness,na.rm=T),
            even_m_NthABG=mean(Evar,na.rm=T),
            even_sd_NthABG=sd(Evar),
            even_cv_NthABG=sd(Evar)/mean(Evar,na.rm=T),
            shan_m_NthABG=mean(Shannon,na.rm=T),
            shan_sd_NthABG=sd(Shannon),
            shan_cv_NthABG=sd(Shannon)/mean(Shannon,na.rm=T))%>%
  group_by(Unit)%>%
  mutate(rich_stab_NthABG=mean(rich_m_NthABG,na.rm=T)/sd(rich_m_NthABG),
         even_stab_NthABG=mean(even_m_NthABG)/sd(even_m_NthABG),
         shan_stab_NthABG=mean(shan_m_NthABG,na.rm=T)/sd(shan_m_NthABG))

#COMBINE ABG and PBG
#north unit
GH_metric_Nth_combo<-PBG_north_GH_metric_master%>%
  group_by(RecYear,iteration,Unit)%>%
  summarise(rich_m_NthPBG=mean(richness, na.rm=T),
            rich_sd_NthPBG=sd(richness),
            rich_cv_NthPBG=sd(richness)/mean(richness,na.rm=T),
            even_m_NthPBG=mean(Evar,na.rm=T),
            even_sd_NthPBG=sd(Evar),
            even_cv_NthPBG=sd(Evar)/mean(Evar,na.rm=T),
            shan_m_NthPBG=mean(Shannon,na.rm=T),
            shan_sd_NthPBG=sd(Shannon),
            shan_cv_NthPBG=sd(Shannon)/mean(Shannon,na.rm=T))%>%
  group_by(Unit,iteration)%>%
  mutate(rich_stab_NthPBG=mean(rich_m_NthPBG,na.rm=T)/sd(rich_m_NthPBG),
         even_stab_NthPBG=mean(even_m_NthPBG)/sd(even_m_NthPBG),
         shan_stab_NthPBG=mean(shan_m_NthPBG,na.rm=T)/sd(shan_m_NthPBG))%>%
  ungroup()%>%
  left_join(ABG_north_GH_metrics, by=c("RecYear","Unit"))%>%
  group_by(RecYear)%>%
  mutate(rich_m_NthPBG_m=mean(rich_m_NthPBG,na.rm=T),
         rich_m_NthPBG_sd=sd(rich_m_NthPBG),
         zscore_rich_m_NthPBG=(rich_m_NthABG-rich_m_NthPBG_m)/rich_m_NthPBG_sd,
         P_rich_m_Nth=2*pnorm(-abs(zscore_rich_m_NthPBG)),
         rich_sd_NthPBG_m=mean(rich_sd_NthPBG),
         rich_sd_NthPBG_sd=sd(rich_sd_NthPBG),
         zscore_rich_sd_NthPBG=(rich_sd_NthABG-rich_sd_NthPBG_m)/rich_sd_NthPBG_sd,
         P_rich_sd_Nth=2*pnorm(-abs(zscore_rich_sd_NthPBG)),
         rich_cv_NthPBG_m=mean(rich_cv_NthPBG),
         rich_cv_NthPBG_sd=sd(rich_cv_NthPBG),
         zscore_rich_cv_NthPBG=(rich_cv_NthABG-rich_cv_NthPBG_m)/rich_cv_NthPBG_sd,
         P_rich_cv_Nth=2*pnorm(-abs(zscore_rich_cv_NthPBG)),
         even_m_NthPBG_m=mean(even_m_NthPBG,na.rm=T),
         even_m_NthPBG_sd=sd(even_m_NthPBG),
         zscore_even_m_NthPBG=(even_m_NthABG-even_m_NthPBG_m)/even_m_NthPBG_sd,
         P_even_m_Nth=2*pnorm(-abs(zscore_even_m_NthPBG)),
         even_sd_NthPBG_m=mean(even_sd_NthPBG),
         even_sd_NthPBG_sd=sd(even_sd_NthPBG),
         zscore_even_sd_NthPBG=(even_sd_NthABG-even_sd_NthPBG_m)/even_sd_NthPBG_sd,
         P_even_sd_Nth=2*pnorm(-abs(zscore_even_sd_NthPBG)),
         even_cv_NthPBG_m=mean(even_cv_NthPBG),
         even_cv_NthPBG_sd=sd(even_cv_NthPBG),
         zscore_even_cv_NthPBG=(even_cv_NthABG-even_cv_NthPBG_m)/even_cv_NthPBG_sd,
         P_even_cv_Nth=2*pnorm(-abs(zscore_even_cv_NthPBG)),
         shan_m_NthPBG_m=mean(shan_m_NthPBG,na.rm=T),
         shan_m_NthPBG_sd=sd(shan_m_NthPBG),
         zscore_shan_m_NthPBG=(shan_m_NthABG-shan_m_NthPBG_m)/shan_m_NthPBG_sd,
         P_shan_m_Nth=2*pnorm(-abs(zscore_shan_m_NthPBG)),
         shan_sd_NthPBG_m=mean(shan_sd_NthPBG),
         shan_sd_NthPBG_sd=sd(shan_sd_NthPBG),
         zscore_shan_sd_NthPBG=(shan_sd_NthABG-shan_sd_NthPBG_m)/shan_sd_NthPBG_sd,
         P_shan_sd_Nth=2*pnorm(-abs(zscore_shan_sd_NthPBG)),
         shan_cv_NthPBG_m=mean(shan_cv_NthPBG),
         shan_cv_NthPBG_sd=sd(shan_cv_NthPBG),
         zscore_shan_cv_NthPBG=(shan_cv_NthABG-shan_cv_NthPBG_m)/shan_cv_NthPBG_sd,
         P_rich_cv_Nth=2*pnorm(-abs(zscore_rich_cv_NthPBG)))%>%
  ungroup()%>%
  mutate(rich_stab_NthPBG_m=mean(rich_stab_NthPBG),
         rich_stab_NthPBG_sd=sd(rich_stab_NthPBG),
         zscore_rich_stab_Nth=(rich_stab_NthABG-rich_stab_NthPBG_m)/rich_stab_NthPBG_sd,
         P_rich_stab_Nth=2*pnorm(-abs(zscore_rich_stab_Nth)),
         even_stab_NthPBG_m=mean(even_stab_NthPBG),
         even_stab_NthPBG_sd=sd(even_stab_NthPBG),
         zscore_even_stab_Nth=(even_stab_NthABG-even_stab_NthPBG_m)/even_stab_NthPBG_sd,
         P_even_stab_Nth=2*pnorm(-abs(zscore_even_stab_Nth)),
         shan_stab_NthPBG_m=mean(shan_stab_NthPBG),
         shan_stab_NthPBG_sd=sd(shan_stab_NthPBG),
         zscore_shan_stab_Nth=(shan_stab_NthABG-shan_stab_NthPBG_m)/shan_stab_NthPBG_sd,
         P_shan_stab_Nth=2*pnorm(-abs(zscore_shan_stab_Nth)))
#calculate for south unit  
  
#figures for community metrics####
#create visual for stab
ggplot(GH_metric_Nth_combo,aes(rich_stab_NthPBG))+
  geom_density(size=1)+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=rich_stab_NthABG), linetype=2,size=1)+
  xlab("Richness Stability North Unit")
ggplot(GH_metric_Nth_combo,aes(even_stab_NthPBG))+
  geom_density(size=1)+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=even_stab_NthABG), linetype=2,size=1)+
  xlab("Evenness Stability North Unit")
ggplot(GH_metric_Nth_combo,aes(shan_stab_NthPBG))+
  geom_density(size=1)+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=shan_stab_NthABG), linetype=2,size=1)+
  xlab("Diversity Stability North Unit")

#convert density plot for mean, sd and cv count to geompoint####
#create visual using point
#richness
GH_metric_Nth_combo_geompoint<-GH_metric_Nth_combo%>%
  pivot_longer(c(rich_m_NthPBG_m,rich_m_NthABG),
               names_to = "treatment", values_to = "rich_mean")
ggplot(GH_metric_Nth_combo_geompoint,aes(RecYear, rich_mean, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=rich_mean+1.96*(rich_m_NthPBG_sd),
                    ymin=rich_mean-1.96*(rich_m_NthPBG_sd)),width=.2)

GH_metric_Nth_combo_geompoint<-GH_metric_Nth_combo%>%
  pivot_longer(c(rich_sd_NthPBG_m,rich_sd_NthABG),
               names_to = "treatment", values_to = "rich_sd")
ggplot(GH_metric_Nth_combo_geompoint,aes(RecYear, rich_sd, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))#+
  geom_errorbar(aes(ymax=rich_sd+1.96*(rich_sd_NthPBG_sd),
                    ymin=rich_sd-1.96*(rich_sd_NthPBG_sd)),width=.2)

GH_metric_Nth_combo_geompoint<-GH_metric_Nth_combo%>%
  pivot_longer(c(rich_cv_NthPBG_m,rich_cv_NthABG),
               names_to = "treatment", values_to = "rich_cv")
ggplot(GH_metric_Nth_combo_geompoint,aes(RecYear, rich_cv, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))#+
  geom_errorbar(aes(ymax=rich_cv+1.96*(rich_cv_NthPBG_sd),
                    ymin=rich_cv-1.96*(rich_cv_NthPBG_sd)),width=.2)


#evenness
GH_metric_Nth_combo_geompoint<-GH_metric_Nth_combo%>%
  pivot_longer(c(even_m_NthPBG_m,even_m_NthABG),
               names_to = "treatment", values_to = "even_mean")
ggplot(GH_metric_Nth_combo_geompoint,aes(RecYear, even_mean, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=even_mean+1.96*(even_m_NthPBG_sd),
                    ymin=even_mean-1.96*(even_m_NthPBG_sd)),width=.2)

GH_metric_Nth_combo_geompoint<-GH_metric_Nth_combo%>%
  pivot_longer(c(even_sd_NthPBG_m,even_sd_NthABG),
               names_to = "treatment", values_to = "even_sd")
ggplot(GH_metric_Nth_combo_geompoint,aes(RecYear, even_sd, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=even_sd+1.96*(even_sd_NthPBG_sd),
                    ymin=even_sd-1.96*(even_sd_NthPBG_sd)),width=.2)

#shannon diversity
GH_metric_Nth_combo_geompoint<-GH_metric_Nth_combo%>%
  pivot_longer(c(shan_m_NthPBG_m,shan_m_NthABG),
               names_to = "treatment", values_to = "diver_mean")
ggplot(GH_metric_Nth_combo_geompoint,aes(RecYear, diver_mean, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=diver_mean+1.96*(shan_m_NthPBG_sd),
                    ymin=diver_mean-1.96*(shan_m_NthPBG_sd)),width=.2)

GH_metric_Nth_combo_geompoint<-GH_metric_Nth_combo%>%
  pivot_longer(c(shan_sd_NthPBG_m,shan_sd_NthABG),
               names_to = "treatment", values_to = "diver_sd")
ggplot(GH_metric_Nth_combo_geompoint,aes(RecYear, diver_sd, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))#+
  geom_errorbar(aes(ymax=diver_sd+1.96*(shan_sd_NthPBG_sd),
                    ymin=diver_sd-1.96*(shan_sd_NthPBG_sd)),width=.2)


