#BPBG grasshopper count and community metrics
#Authors: Joshua Ajowele
#Started: 26 May 2022 last modified: 01 May 2024

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
unique(grassh_count_df$spe)
#merging key with dataset
grassh_count_df <- grasshopperspcomp_df%>%
  full_join(watershed_key, by = "Watershed")%>%
  left_join(Watershed_key2,by="Watershed")%>%
  mutate(spe=paste(Genus,Species, sep="_"))%>%
  mutate(RecYear = Recyear)%>%
  filter(!RecYear== 2010)%>%
  filter(!spe%in%c("Oecanthinae_spp.","Tettigoniidae_spp.","Gryllidae_spp.",
                   "Conocephalus_spp.","Neoconocephalus_robustus","Scudderia_texensis",
                   "Arethaea_constricta","Orchelimum_spp.","Amblycorypha_oblongifolia","Pediodectes_haldemani",
                   "Amblycorypha_rotundifolia","Neoconocephalus_spp.","Neoconocephalus_ensiger","Pediodectes_nigromarginatus",
                   "Scudderia_furcata","Scudderia_spp."))%>%
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
               names_to = "treatment", values_to = "count_sd")%>%
  select(count_sd,treatment, RecYear, count_PBGNth_sd_sd)%>%
  unique()%>%
  #round to two decimal places
  mutate(sdsd=round(count_PBGNth_sd_sd,digits=2))%>%
  mutate(sdsd=ifelse(sdsd==3.36 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==6.82 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==31.33 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==4.45 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==8.68 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==4.27 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==5.11 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==6.25 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==7.31 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==4.10 & treatment=="count_ABGNth_sd",0,sdsd),
         sdsd=ifelse(sdsd==6.38 & treatment=="count_ABGNth_sd",0,sdsd))

ggplot(combo_north_sd_geompoint,aes(RecYear, count_sd, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=count_sd+1.96*(sdsd),
                    ymin=count_sd-1.96*(sdsd)),width=.2)


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
  filter(!spe%in%c("Oecanthinae_spp.","Tettigoniidae_spp.","Gryllidae_spp.",
                   "Conocephalus_spp.","Neoconocephalus_robustus","Scudderia_texensis",
                   "Arethaea_constricta","Orchelimum_spp.","Amblycorypha_oblongifolia","Pediodectes_haldemani",
                   "Amblycorypha_rotundifolia","Neoconocephalus_spp.","Neoconocephalus_ensiger","Pediodectes_nigromarginatus",
                   "Scudderia_furcata","Scudderia_spp."))%>%
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


#Combine unit for bootstrap sd and cv####
grass_ABG_count<-grassh_count_df%>%
    #ABG mean, cv, sd
    filter(FireGrzTrt=="ABG")%>%
    group_by(RecYear)%>%
    summarise(ghcount_ABG=mean(Tcount, na.rm=T),
              ghcount_ABG_sd=sd(Tcount),
              ghcount_ABG_cv=ghcount_ABG_sd/ghcount_ABG)%>%
    ungroup()%>%
    mutate(ghcount_ABG_stab=mean(ghcount_ABG)/sd(ghcount_ABG))


gh_count_index<-grassh_count_df%>%
  filter(FireGrzTrt=="PBG")%>%
  group_by(RecYear)%>%
  mutate(plot_index=1:length(RecYear))

num_bootstrap<-1000
bootstrap_vector<-1:num_bootstrap
ghcount_master_combine_mean_sd_cv<-{}
  for(BOOT in 1:length(bootstrap_vector)){
    count_rand_mean<-gh_count_index%>%
      dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Repsite, plot_index)%>%
      unique()%>%
      group_by(RecYear)%>%
      sample_n(8, replace=T)%>%
      dplyr::select(RecYear,plot_index)%>%
      left_join(gh_count_index, by=c("plot_index","RecYear"), multiple="all")
    count_ready_mean<-count_rand_mean%>%
      mutate(iteration=BOOT)
    ghcount_master_combine_mean_sd_cv<-rbind(ghcount_master_combine_mean_sd_cv, count_ready_mean)
  }
gh_combine_mean_sd<-ghcount_master_combine_mean_sd_cv%>%
  select(RecYear, iteration, Tcount)%>%
  group_by(RecYear, iteration)%>%
  summarise(gh_PBG_mean=mean(Tcount,na.rm=T),
            gh_PBG_sd=sd(Tcount),
            gh_PBG_cv=gh_PBG_sd/gh_PBG_mean)%>%
  ungroup()%>%
  group_by(RecYear)%>%
  summarise(gh_PBG_MM=mean(gh_PBG_mean),
            gh_PBG_M_sd=sd(gh_PBG_mean),
            gh_PBG_sd_M=mean(gh_PBG_sd),
            gh_PBG_sd_sd=sd(gh_PBG_sd),
            gh_PBG_cv_M=mean(gh_PBG_cv),
            gh_PBG_cv_sd=sd(gh_PBG_cv))%>%
  left_join(grass_ABG_count, by="RecYear")%>%
  mutate(zscore_M=(ghcount_ABG-gh_PBG_MM)/gh_PBG_M_sd,
         pval_M=2*pnorm(-abs(zscore_M)),
         zscore_sd=(ghcount_ABG_sd-gh_PBG_sd_M)/gh_PBG_sd_sd,
         pval_sd=2*pnorm(-abs(zscore_sd)),
         zscore_cv=(ghcount_ABG_cv-gh_PBG_cv_M)/gh_PBG_cv_sd,
         pval_cv=2*pnorm(-abs(zscore_cv)))
write.csv(gh_combine_mean_sd, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/gh_count_mean_sd_combined_unit.csv")

#create visuals####
#grasshopper average count
combo_ghcount_geompoint_M<-gh_combine_mean_sd%>%
  pivot_longer(c(gh_PBG_MM,ghcount_ABG),
               names_to = "treatment", values_to = "count_M")%>%
  select(treatment,count_M,RecYear, gh_PBG_M_sd)%>%
  distinct()%>%
  mutate(PBG_sd=round(gh_PBG_M_sd,digits=2))#%>%
  mutate(PBG_sd=ifelse(PBG_sd==2.55 & treatment=="ghcount_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==6.75 & treatment=="ghcount_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==24.97 & treatment=="ghcount_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==5.42 & treatment=="ghcount_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==6.82 & treatment=="ghcount_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==1.95 & treatment=="ghcount_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==6.36 & treatment=="ghcount_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==13.35 & treatment=="ghcount_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==7.21 & treatment=="ghcount_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==5.21 & treatment=="ghcount_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==6.12 & treatment=="ghcount_ABG",NA,PBG_sd))

ggplot(combo_ghcount_geompoint_M,aes(RecYear, count_M, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#009E73","#F0E442"),labels=c("PBG","ABG"))+
  geom_errorbar(aes(ymax=count_M+1.96*(PBG_sd),
                    ymin=count_M-1.96*(PBG_sd)),width=.2)

  
#grasshopper sd count
combo_ghcount_geompoint_sd<-gh_combine_mean_sd%>%
  pivot_longer(c(gh_PBG_sd_M,ghcount_ABG_sd),
               names_to = "treatment", values_to = "count_sd")%>%
  select(treatment,count_sd,RecYear, gh_PBG_sd_sd)%>%
  distinct()%>%
  mutate(PBG_sd=round(gh_PBG_sd_sd,digits=2))%>%
  mutate(PBG_sd=ifelse(PBG_sd==1.89 & treatment=="ghcount_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==6.23 & treatment=="ghcount_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==21.59 & treatment=="ghcount_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==3.68 & treatment=="ghcount_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==6.86 & treatment=="ghcount_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==1.98 & treatment=="ghcount_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==9.67 & treatment=="ghcount_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==20.06 & treatment=="ghcount_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==4.30 & treatment=="ghcount_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==6.90 & treatment=="ghcount_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==4.73 & treatment=="ghcount_ABG_sd",NA,PBG_sd))

ggplot(combo_ghcount_geompoint_sd,aes(RecYear, count_sd, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#009E73","#F0E442"),labels=c("PBG","ABG"))+
  geom_errorbar(aes(ymax=count_sd+1.96*(PBG_sd),
                    ymin=count_sd-1.96*(PBG_sd)),width=.2)

#grasshopper cv count
combo_ghcount_geompoint_cv<-gh_combine_mean_sd%>%
  pivot_longer(c(gh_PBG_cv_M,ghcount_ABG_cv),
               names_to = "treatment", values_to = "count_cv")%>%
  select(treatment,count_cv,RecYear, gh_PBG_cv_sd)%>%
  distinct()%>%
  mutate(PBG_sd=round(gh_PBG_cv_sd,digits=2))%>%
  mutate(PBG_sd=ifelse(PBG_sd==0.05 & treatment=="ghcount_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.16 & treatment=="ghcount_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.13 & treatment=="ghcount_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.08 & treatment=="ghcount_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.19 & treatment=="ghcount_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.10 & treatment=="ghcount_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.20 & treatment=="ghcount_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.25 & treatment=="ghcount_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.07 & treatment=="ghcount_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.19 & treatment=="ghcount_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.11 & treatment=="ghcount_ABG_cv",NA,PBG_sd))

ggplot(combo_ghcount_geompoint_cv,aes(RecYear, count_cv, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#009E73","#F0E442"),labels=c("PBG","ABG"))+
  geom_errorbar(aes(ymax=count_cv+1.96*(PBG_sd),
                    ymin=count_cv-1.96*(PBG_sd)),width=.2)

#bootstrap for stability####
num_bootstrap<-1000
bootstrap_vector<-1:num_bootstrap
ghcount_master_stab<-{}
for(BOOT in 1:length(bootstrap_vector)){
  count_rand_stab<-gh_count_index%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Repsite, plot_index)%>%
    unique()%>%
    filter(RecYear==2011)%>%
    ungroup()%>%
    sample_n(8, replace=T)%>%
    dplyr::select(plot_index)%>%
    left_join(gh_count_index, by="plot_index", multiple="all")
  count_ready_stab<-count_rand_stab%>%
    mutate(iteration=BOOT)
  ghcount_master_stab<-rbind(ghcount_master_stab, count_ready_stab)
}
#just so you don't worry-2012-C03B transect A is missing(i think bad sample)
ghcount_master_combine_stab<-ghcount_master_stab%>%
    select(RecYear, iteration, Tcount)%>%
  group_by(RecYear, iteration)%>%
  summarise(count_mean=mean(Tcount,na.rm=T))%>%
  ungroup()%>%
  group_by(iteration)%>%
  mutate(count_PBG_stab=mean(count_mean, na.rm=T)/sd(count_mean))%>%
  ungroup()%>%
  mutate(PBG_stab_M=mean(count_PBG_stab),
         PBG_stab_sd=sd(count_PBG_stab))%>%
  left_join(grass_ABG_count, by="RecYear")%>%
  mutate(zscore=(ghcount_ABG_stab-PBG_stab_M)/PBG_stab_sd)%>%
  mutate(pval=2*pnorm(-abs(zscore)))
write.csv(ghcount_master_combine_stab, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/gh_stability_combined_unit.csv")
#start from zero
ggplot(ghcount_master_combine_stab,aes(count_PBG_stab))+
  geom_density(size=2,col="#009E73")+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=ghcount_ABG_stab), linetype=2,size=2, col="#F0E442")+
  xlab("Grasshopper count stability")+xlim(0,4)
#x axis not set
ggplot(ghcount_master_combine_stab,aes(count_PBG_stab))+
  geom_density(size=2,col="#009E73")+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=ghcount_ABG_stab), linetype=2,size=2, col="#F0E442")+
  xlab("Grasshopper count stability")

#next step year since fire calculation
#creating a key for year since fire
YrSinceFire_key <- tibble(year_watershed= c("2011_C01A", "2011_C03A", "2011_C03B",
                                            "2011_C03C", "2011_C1SB", "2011_C3SA",
                                            "2011_C3SB", "2011_C3SC",
                                            "2012_C01A", "2012_C03A", "2012_C03B",
                                            "2012_C03C", "2012_C1SB", "2012_C3SA",
                                            "2012_C3SB", "2012_C3SC", "2013_C03A",
                                            "2013_C03B", "2013_C03C", "2013_C1SB",
                                            "2013_C3SA", "2013_C3SB", "2013_C01A",
                                            "2013_C3SC", "2014_C01A", "2014_C03A",
                                            "2014_C03B", "2014_C03C", "2014_C1SB",
                                            "2014_C3SA", "2014_C3SB", "2014_C3SC",
                                            "2015_C1SB", "2015_C03A", "2015_C03B",
                                            "2015_C3SC", "2015_C3SA", "2015_C3SB",
                                            "2015_C03C", "2015_C01A", "2016_C03A",
                                            "2016_C03B", "2016_C3SC", "2016_C01A",
                                            "2016_C03C", "2016_C1SB", "2016_C3SA",
                                            "2016_C3SB", "2017_C3SA", "2017_C3SB",
                                            "2017_C03C", "2017_C03A", "2017_C03B",
                                            "2017_C1SB", "2017_C3SC", "2017_C01A",
                                            "2018_C03A", "2018_C03B", "2018_C03C",
                                            "2018_C01A", "2018_C3SA", "2018_C3SB",
                                            "2018_C1SB", "2018_C3SC", "2019_C3SA",
                                            "2019_C3SB", "2019_C1SB", "2019_C03A",
                                            "2019_C03B", "2019_C03C", "2019_C3SC",
                                            "2019_C01A", "2020_C1SB", "2020_C3SA",
                                            "2020_C3SB", "2020_C3SC", "2020_C03A",
                                            "2020_C03B", "2020_C03C", "2020_C01A",
                                            "2021_C01A", "2021_C03C", "2021_C03A",
                                            "2021_C03B", "2021_C3SA", "2021_C3SB",
                                            "2021_C3SC", "2021_C1SB"),
                          yrsins_fire= c("ABG0","PBG1","PBG0", "PBG2", "ABG0","PBG0","PBG1","PBG1",
                                         "ABG0","PBG2","PBG1","PBG0","ABG0","PBG1","PBG2",
                                         "PBG0","PBG0", "PBG2", "PBG1", "ABG0","PBG2", "PBG0",
                                         "ABG0", "PBG1","ABG0","PBG1", "PBG0","PBG2","ABG0",
                                         "PBG0","PBG1","PBG2","ABG0","PBG2","PBG1","PBG0","PBG1",
                                         "PBG2","PBG0","ABG0","PBG0","PBG2","PBG1","ABG0","PBG1",
                                         "ABG0","PBG2","PBG0","PBG0","PBG1","PBG2","PBG1","PBG0",
                                         "ABG0","PBG2","ABG0","PBG2","PBG1","PBG0","ABG0","PBG1",
                                         "PBG2","ABG0","PBG0","PBG2","PBG0","ABG0","PBG0","PBG2",
                                         "PBG1","PBG1","ABG0","ABG0","PBG0","PBG1","PBG2","PBG1",
                                         "PBG0","PBG2","ABG0","ABG0","PBG0","PBG2","PBG1","PBG1",
                                         "PBG2","PBG0","ABG0"))
#merge with biomass data ####
ghcount_yrs<-grassh_count_df%>%
  mutate(year_watershed=paste(RecYear,Watershed, sep ="_"))%>%
  left_join(YrSinceFire_key, by="year_watershed")%>%
  group_by(RecYear,Unit, Watershed,yrsins_fire,FireGrzTrt,Repsite)%>%
  summarise(Tcount=mean(Tcount, na.rm=T))
ghcount_yrs$RecYear=as.factor(ghcount_yrs$RecYear)
ghcount_yrs$yrsins_fire=as.factor(ghcount_yrs$yrsins_fire)
ghcount_yrs$Unit=as.factor(ghcount_yrs$Unit)
ghcount_yrs$Watershed=as.factor(ghcount_yrs$Watershed)
ghcount_yrs$Repsite=as.factor(ghcount_yrs$Repsite)


#mixed anova
yrs_ghcount_model<-lmer(log(Tcount)~yrsins_fire*RecYear+(1|Unit/Watershed),
                        data=ghcount_yrs)
check_model(yrs_ghcount_model)
Anova(yrs_ghcount_model, type=3)
qqnorm(resid(yrs_ghcount_model))

#multiple comparison
testInteractions(yrs_ghcount_model, pairwise="yrsins_fire", fixed="RecYear",
                 adjustment="BH")
#using mean estimate to create figure 
yrs_interact<-interactionMeans(yrs_ghcount_model)
#replacing spaces in column names with underscore 
names(yrs_interact)<-str_replace_all(names(yrs_interact), " ","_")
#df for visuals from model estimates
yrs_interact_viz<-yrs_interact%>%
  mutate(yrs_interact_bt_mean=exp(adjusted_mean),
         yrs_interact_bt_upper=exp(adjusted_mean+SE_of_link),
         yrs_interact_bt_lower=exp(adjusted_mean-SE_of_link))
#visual
ggplot(yrs_interact_viz,aes(RecYear, yrs_interact_bt_mean, col=yrsins_fire, linetype=yrsins_fire))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=yrs_interact_bt_lower,
                    ymax=yrs_interact_bt_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73", "#999999", "#0072B2"))

#average across years for simplification
yrs_interact_bar<-yrs_interact_viz%>%
  group_by(yrsins_fire)%>%
  summarise(Grasshopper_count=mean(yrs_interact_bt_mean),
            se_upper=mean(yrs_interact_bt_upper),
            se_lower=mean(yrs_interact_bt_lower))
ggplot(yrs_interact_bar,aes(x=yrsins_fire,fill=yrsins_fire))+
  geom_bar(stat = "identity",aes(y=Grasshopper_count),width = 0.5)+
  geom_errorbar(aes(ymin=se_lower,
                    ymax=se_upper),width=0.2,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73", "#999999", "#0072B2"))


#create df with individual species abundance####
grassh_comm_df <- grasshopperspcomp_df%>%
  full_join(watershed_key, by = "Watershed")%>%
  left_join(Watershed_key2,by="Watershed")%>%
  mutate(spe=paste(Genus,Species, sep="_"))%>%
  mutate(RecYear = Recyear)%>%
  filter(!RecYear== 2010)%>%
  filter(!spe%in%c("Oecanthinae_spp.","Tettigoniidae_spp.","Gryllidae_spp.",
                   "Conocephalus_spp.","Neoconocephalus_robustus","Scudderia_texensis",
                   "Arethaea_constricta","Orchelimum_spp.","Amblycorypha_oblongifolia","Pediodectes_haldemani",
                   "Amblycorypha_rotundifolia","Neoconocephalus_spp.","Neoconocephalus_ensiger","Pediodectes_nigromarginatus",
                   "Scudderia_furcata","Scudderia_spp."))%>%
  #using the max cover from the two sweeps done on each transect
  group_by(Unit,RecYear,FireGrzTrt,Watershed,Repsite,spe)%>%
  summarise(Total=max(Total))%>%
  filter(Repsite!="C" | Watershed!="C03C")%>%
  filter(Watershed!="C03C" | Repsite!="D")%>%
  filter(Watershed!="C03B" | Repsite!="A")%>%
  filter(Watershed!="C03B" | Repsite!="B")%>%
  filter(Watershed!="C03B" | Repsite!="C")%>%
  filter(Watershed!="C03A" | Repsite!="A")%>%
  filter(Watershed!="C03A" | Repsite!="B")%>%
  filter(Watershed!="C03A" | Repsite!="C")%>%
  filter(Watershed!="C3SC" | Repsite!="A")%>%
  filter(Watershed!="C3SC" | Repsite!="B")%>%
  filter(Watershed!="C3SC" | Repsite!="C")%>%
  filter(Watershed!="C3SA" | Repsite!="B")%>%
  filter(Watershed!="C3SA" | Repsite!="C")%>%
  filter(Watershed!="C3SA" | Repsite!="D")%>%
  filter(Watershed!="C3SB" | Repsite!="A")%>%
  filter(Watershed!="C3SB" | Repsite!="B")%>%
  group_by(Unit,RecYear,FireGrzTrt,spe)%>%
  summarise(total_avg=mean(Total,na.rm=T))%>%
  group_by(Unit,RecYear,FireGrzTrt)%>%
  #converting count data to abundance data (0-100% scale)
  mutate(abundance=(total_avg/sum(total_avg, na.rm=T))*100)%>%
  group_by(Unit,RecYear,FireGrzTrt,spe)%>%
  mutate(Rep_id=paste(Unit,FireGrzTrt, sep="_"))

#community metrics with codyn####
#deriving richness and evenness using codyn at landscape scale
gh_rich <- community_structure(grassh_comm_df, time.var = "RecYear", 
                                 abundance.var = "abundance",
                                 replicate.var = "Rep_id", metric = "Evar")

#deriving diversity values
gh_diver <- community_diversity(grassh_comm_df, time.var = "RecYear", 
                                  abundance.var = "abundance",
                                  replicate.var = "Rep_id", metric="Shannon")


#extracting important columns 
grassh_comm_subset <- grassh_comm_df %>%
  ungroup()%>%
  dplyr::select(RecYear,FireGrzTrt,Rep_id, Unit) %>%
  distinct()#remove repeated rows

#combine into a single data
grassh_rich_diver<- gh_diver %>%
  #combine richness and diversity
  left_join(gh_rich, by = c("Rep_id","RecYear")) %>%
  #combine with important columns
  left_join(grassh_comm_subset, by = c("Rep_id","RecYear"))
#convert to factors
grassh_rich_diver$RecYear<-as.factor(grassh_rich_diver$RecYear)
grassh_rich_diver$Unit<-as.factor(grassh_rich_diver$Unit)
grassh_rich_diver$FireGrzTrt<-as.factor(grassh_rich_diver$FireGrzTrt)
#build model
gh_diver_model<-lmer(log(Shannon)~FireGrzTrt*RecYear+(1|Unit),
                     data=grassh_rich_diver)
Anova(gh_diver_model, type=3)
qqnorm(resid(gh_diver_model))


#multiple comparison
testInteractions(gh_diver_model, pairwise="FireGrzTrt", fixed="RecYear",
                 adjustment="BH")
#using mean estimate to create figure 
ghdiver_interact<-interactionMeans(gh_diver_model)
#replacing spaces in column names with underscore 
names(ghdiver_interact)<-str_replace_all(names(ghdiver_interact), " ","_")
#df for visuals from model estimates
ghdiver_interact_viz<-ghdiver_interact%>%
  mutate(ghdiver_mean=exp(adjusted_mean),
         ghdiver_upper=exp(adjusted_mean+SE_of_link),
         ghdiver_lower=exp(adjusted_mean-SE_of_link))
#visual
ggplot(ghdiver_interact_viz,aes(RecYear, ghdiver_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=ghdiver_lower,
                    ymax=ghdiver_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))
#average across years for simplification
ghdiver_interact_bar<-ghdiver_interact_viz%>%
  group_by(FireGrzTrt)%>%
  summarise(gh_diver=mean(ghdiver_mean),
            se_upper=mean(ghdiver_upper),
            se_lower=mean(ghdiver_lower))
ggplot(ghdiver_interact_bar,aes(x=FireGrzTrt,fill=FireGrzTrt))+
  geom_bar(stat = "identity",aes(y=gh_diver),width = 0.5)+
  geom_errorbar(aes(ymin=se_lower,
                    ymax=se_upper),width=0.2,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"))

#richness model
ghrich_model<-lmer(log(richness)~FireGrzTrt*RecYear+(1|Unit),
                   data=grassh_rich_diver)
Anova(ghrich_model, type=3)
qqnorm(resid(ghrich_model))
check_model(ghrich_model)
#multiple comparison
testInteractions(ghrich_model, pairwise="FireGrzTrt", fixed="RecYear",
                 adjustment="BH")
#using mean estimate to create figure 
ghrich_interact<-interactionMeans(ghrich_model)
#replacing spaces in column names with underscore 
names(ghrich_interact)<-str_replace_all(names(ghrich_interact), " ","_")
#df for visuals from model estimates
ghrich_interact_viz<-ghrich_interact%>%
  mutate(ghrich_mean=exp(adjusted_mean),
         ghrich_upper=exp(adjusted_mean+SE_of_link),
         ghrich_lower=exp(adjusted_mean-SE_of_link))
#visual
ggplot(ghrich_interact_viz,aes(RecYear, ghrich_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=ghrich_lower,
                    ymax=ghrich_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

#average across years for simplification
ghrich_interact_bar<-ghrich_interact_viz%>%
  group_by(FireGrzTrt)%>%
  summarise(gh_richness=mean(ghrich_mean),
            se_upper=mean(ghrich_upper),
            se_lower=mean(ghrich_lower))
ggplot(ghrich_interact_bar,aes(x=FireGrzTrt,fill=FireGrzTrt))+
  geom_bar(stat = "identity",aes(y=gh_richness),width = 0.5)+
  geom_errorbar(aes(ymin=se_lower,
                    ymax=se_upper),width=0.2,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"))

#evenness model
ghevenness_model<-lmer(log(Evar)~FireGrzTrt*RecYear+(1|Unit),
                   data=grassh_rich_diver)
Anova(ghevenness_model, type=3)
qqnorm(resid(ghevenness_model))
check_model(ghevenness_model)
#multiple comparison
testInteractions(ghevenness_model, pairwise="FireGrzTrt", fixed="RecYear",
                 adjustment="BH")
#using mean estimate to create figure 
gheven_interact<-interactionMeans(ghevenness_model)
#replacing spaces in column names with underscore 
names(gheven_interact)<-str_replace_all(names(gheven_interact), " ","_")
#df for visuals from model estimates
gheven_interact_viz<-gheven_interact%>%
  mutate(ghevenness_mean=exp(adjusted_mean),
         gheven_upper=exp(adjusted_mean+SE_of_link),
         gheven_lower=exp(adjusted_mean-SE_of_link))
#visual
ggplot(gheven_interact_viz,aes(RecYear, ghevenness_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=gheven_lower,
                    ymax=gheven_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

#average across years for simplification
gheven_interact_bar<-gheven_interact_viz%>%
  group_by(FireGrzTrt)%>%
  summarise(gh_evenness=mean(ghevenness_mean),
            se_upper=mean(gheven_upper),
            se_lower=mean(gheven_lower))
ggplot(gheven_interact_bar,aes(x=FireGrzTrt,fill=FireGrzTrt))+
  geom_bar(stat = "identity",aes(y=gh_evenness),width = 0.5)+
  geom_errorbar(aes(ymin=se_lower,
                    ymax=se_upper),width=0.2,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"))

#perform permanova and calculate betadiversity####
#create df with individual species abundance####
grassh_permav_df <- grasshopperspcomp_df%>%
  full_join(watershed_key, by = "Watershed")%>%
  left_join(Watershed_key2,by="Watershed")%>%
  mutate(spe=paste(Genus,Species, sep="_"))%>%
  mutate(RecYear = Recyear)%>%
  filter(!RecYear== 2010)%>%
  filter(!spe%in%c("Oecanthinae_spp.","Tettigoniidae_spp.","Gryllidae_spp.",
                   "Conocephalus_spp.","Neoconocephalus_robustus","Scudderia_texensis",
                   "Arethaea_constricta","Orchelimum_spp.","Amblycorypha_oblongifolia","Pediodectes_haldemani",
                   "Amblycorypha_rotundifolia","Neoconocephalus_spp.","Neoconocephalus_ensiger","Pediodectes_nigromarginatus",
                   "Scudderia_furcata","Scudderia_spp."))%>%
  #using the max cover from the two sweeps done on each transect
  group_by(Unit,RecYear,FireGrzTrt,Watershed,Repsite,spe)%>%
  summarise(Total=max(Total))%>%
  filter(Repsite!="C" | Watershed!="C03C")%>%
  filter(Watershed!="C03C" | Repsite!="D")%>%
  filter(Watershed!="C03B" | Repsite!="A")%>%
  filter(Watershed!="C03B" | Repsite!="B")%>%
  filter(Watershed!="C03B" | Repsite!="C")%>%
  filter(Watershed!="C03A" | Repsite!="A")%>%
  filter(Watershed!="C03A" | Repsite!="B")%>%
  filter(Watershed!="C03A" | Repsite!="C")%>%
  filter(Watershed!="C3SC" | Repsite!="A")%>%
  filter(Watershed!="C3SC" | Repsite!="B")%>%
  filter(Watershed!="C3SC" | Repsite!="C")%>%
  filter(Watershed!="C3SA" | Repsite!="B")%>%
  filter(Watershed!="C3SA" | Repsite!="C")%>%
  filter(Watershed!="C3SA" | Repsite!="D")%>%
  filter(Watershed!="C3SB" | Repsite!="A")%>%
  filter(Watershed!="C3SB" | Repsite!="B")%>%
  group_by(Unit,RecYear,Watershed,Repsite,FireGrzTrt,spe)%>%
  summarise(total_avg=mean(Total,na.rm=T))%>%
  group_by(Unit,RecYear,Watershed,Repsite,FireGrzTrt)%>%
  #converting count data to abundance data (0-100% scale)
  mutate(abundance=(total_avg/sum(total_avg, na.rm=T))*100,
         unit_trt=paste(Unit, FireGrzTrt, sep="_"))%>%
  group_by(Unit,RecYear,Watershed,Repsite,FireGrzTrt,unit_trt,spe)%>%
  summarise(abundance=mean(abundance,na.rm=T))%>%
  pivot_wider(names_from = spe, values_from = abundance, values_fill = 0)

# Separate out spcomp and environmental columns (cols are species) #
Gh_sp_data <- grassh_permav_df %>%
  ungroup()%>%
  dplyr::select(-1:-6)
Gh_env_data <- grassh_permav_df%>%dplyr::select(1:6)

#get nmds1 and 2
Gh_mds_all <- metaMDS(Gh_sp_data, distance = "bray")
#Run 20 stress 0.2695

#combine NMDS1 and 2 with factor columns and create centroids
Gh_mds_scores <- data.frame(Gh_env_data, scores(Gh_mds_all, display="sites"))%>%
  group_by(RecYear, FireGrzTrt)%>%
  mutate(NMDS1_mean=mean(NMDS1),
         NMDS2_mean=mean(NMDS2))
#write into csv format 
#write.csv(Gh_mds_scores,"C:/Users/joshu/OneDrive - UNCG/UNCG PHD/PhD Wyoming_One Drive/PHD Wyoming/Thesis/PBG synthesis/Gh_mds_scores.csv")

#for ordispider sake!
Gh_mds_scores_mean <- Gh_mds_scores%>%
  group_by(RecYear, FireGrzTrt)%>%
  summarise(NMDS1=mean(NMDS1),
            NMDS2=mean(NMDS2))
#ordispider with ggplot
ggplot(Gh_mds_scores, aes(x=NMDS1, y=NMDS2, col=FireGrzTrt)) +
  geom_segment(aes(xend=NMDS1_mean, yend= NMDS2_mean))+
  geom_point(data=Gh_mds_scores_mean, size=5) +
  geom_point()+
  facet_wrap(~ RecYear, scales = "free")+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))  #+
theme_bw() 

#summarise into centroid and have a single figure
#Plot with geompoint for abg vs pbg nmds
ggplot(Gh_mds_scores_mean, aes(x=NMDS1, y=NMDS2, col=FireGrzTrt))+
  geom_point(size=5)+
  geom_path()+
  #geom_errorbar(aes(ymax=NMDS2_mean+NMDS2_SE, ymin=NMDS2_mean-NMDS2_SE))+
  #geom_errorbarh(aes(xmax=NMDS1_mean+NMDS1_SE, xmin=NMDS1_mean-NMDS1_SE))+
  scale_shape_manual(values=10:21)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

#not necessary-just checking if composition differs with unit####
Gh_mds_scores_unit <- data.frame(Gh_env_data, scores(Gh_mds_all, display="sites"))%>%
  group_by(RecYear, unit_trt)%>%
  mutate(NMDS1_mean=mean(NMDS1),
         NMDS2_mean=mean(NMDS2))

#for ordispider sake!
Gh_mds_scores_mean_unit <- Gh_mds_scores_unit%>%
  group_by(RecYear, unit_trt)%>%
  summarise(NMDS1=mean(NMDS1),
            NMDS2=mean(NMDS2))
#ordispider with ggplot
ggplot(Gh_mds_scores_unit, aes(x=NMDS1, y=NMDS2, col=unit_trt)) +
  geom_segment(aes(xend=NMDS1_mean, yend= NMDS2_mean))+
  geom_point(data=Gh_mds_scores_mean_unit, size=5) +
  geom_point()+
  facet_wrap(~ RecYear, scales = "free")#+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))  #+
theme_bw() 


#permanova and betdiversity####
#calculating permanova and betadiversity
#creating a loop to do this
year_vec_gh <- unique(Gh_env_data$RecYear)
gh_perm <- {}
gh_beta <- {}


for(YEAR in 1:length(year_vec_gh)){
  vdist_temp_gh <- vegdist(filter(Gh_sp_data, Gh_env_data$RecYear ==  year_vec_gh[YEAR]))
  permanova_temp_gh <- adonis(vdist_temp_gh ~ subset(Gh_env_data, RecYear == year_vec_gh[YEAR])$FireGrzTrt)
  permanova_out_temp_gh <- data.frame(RecYear = year_vec_gh[YEAR], 
                                      DF = permanova_temp_gh$aov.tab[1,1],
                                      F_value = permanova_temp_gh$aov.tab[1,4],
                                      P_value = permanova_temp_gh$aov.tab[1,6])
  gh_perm <- rbind(gh_perm,permanova_out_temp_gh)
  
  bdisp_temp_gh <- betadisper(vdist_temp_gh, filter(Gh_env_data, RecYear==year_vec_gh[YEAR])$FireGrzTrt, type = "centroid")
  bdisp_out_temp_gh <- data.frame(filter(Gh_env_data, RecYear==year_vec_gh[YEAR]), distance = bdisp_temp_gh$distances)
  gh_beta <- rbind(gh_beta, bdisp_out_temp_gh)
  
  rm(vdist_temp_gh, permanova_temp_gh, permanova_out_temp_gh, bdisp_temp_gh, bdisp_out_temp_gh)
}
write.csv(gh_beta, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/Grasshoper_betadiver.csv")
write.csv(gh_perm, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/Grasshopper_permanova.csv")
#model for betadiversity
gh_beta$RecYear<-as.factor(gh_beta$RecYear)
gh_beta_model<-lmer(log(distance)~FireGrzTrt*RecYear+(1|Unit),
                    data=gh_beta)#issingular
Anova(gh_beta_model, type=3)
qqnorm(resid(gh_beta_model))
check_model(gh_beta_model)
#pairwise interaction 
testInteractions(gh_beta_model, fixed="RecYear",
                 pairwise = "FireGrzTrt", adjustment="BH")
#using mean estimate from post-hoc to create figure as a comparison to raw data
model_estimates_beta<-interactionMeans(gh_beta_model)
#replacing spaces in column names with underscore 
names(model_estimates_beta)<-str_replace_all(names(model_estimates_beta), " ","_")
#df for visuals from model estimates
model_estimates_beta_viz<-model_estimates_beta%>%
  mutate(gh_betadiv=exp(adjusted_mean),
         gh_upper=exp(adjusted_mean+SE_of_link),
         gh_lower=exp(adjusted_mean-SE_of_link))
#visual
ggplot(model_estimates_beta_viz,aes(RecYear, gh_betadiv, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=gh_lower,
                    ymax=gh_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))
#summarize with a bargraph
ghbeta_interact_bar<-model_estimates_beta_viz%>%
  group_by(FireGrzTrt)%>%
  summarise(gh_betadiver=mean(gh_betadiv),
            se_upper=mean(gh_upper),
            se_lower=mean(gh_lower))
ggplot(ghbeta_interact_bar,aes(x=FireGrzTrt,fill=FireGrzTrt))+
  geom_bar(stat = "identity",aes(y=gh_betadiver),width = 0.5)+
  geom_errorbar(aes(ymin=se_lower,
                    ymax=se_upper),width=0.2,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"))

#simper analysis to find what species are driving difference in composition####
#filter data for year with difference
Gh_sp_data_2011 <- grassh_permav_df %>%
  filter(RecYear==2011)%>%
  ungroup()%>%
  dplyr::select(-1:-6)
Gh_env_data_2011 <- grassh_permav_df%>%
  dplyr::select(1:6)%>%
  filter(RecYear==2011)

simper_abgvspbg<-simper(Gh_sp_data_2011,Gh_env_data_2011$FireGrzTrt)
simper_gh_2011<-summary(simper_abgvspbg, order=T)
Simper_df<-data.frame(simper_gh_2011$ABG_PBG)


#2017
Gh_sp_data_2017 <- grassh_permav_df %>%
  filter(RecYear==2017)%>%
  ungroup()%>%
  dplyr::select(-1:-6)
Gh_env_data_2017 <- grassh_permav_df%>%
  dplyr::select(1:6)%>%
  filter(RecYear==2017)

simper_gh_2017<-simper(Gh_sp_data_2017,Gh_env_data_2017$FireGrzTrt)
simper_gh_2017_sum<-summary(simper_gh_2017, order=T)
Simper_df<-data.frame(simper_gh_2017_sum$ABG_PBG)
#2018
Gh_sp_data_2018 <- grassh_permav_df %>%
  filter(RecYear==2018)%>%
  ungroup()%>%
  dplyr::select(-1:-6)
Gh_env_data_2018 <- grassh_permav_df%>%
  dplyr::select(1:6)%>%
  filter(RecYear==2018)

simper_gh_2018<-simper(Gh_sp_data_2018,Gh_env_data_2018$FireGrzTrt)
simper_gh_2018_sum<-summary(simper_gh_2018, order=T)
Simper_df<-data.frame(simper_gh_2018_sum$ABG_PBG)
#2019
Gh_sp_data_2019 <- grassh_permav_df %>%
  filter(RecYear==2019)%>%
  ungroup()%>%
  dplyr::select(-1:-6)
Gh_env_data_2019 <- grassh_permav_df%>%
  dplyr::select(1:6)%>%
  filter(RecYear==2019)

simper_gh_2019<-simper(Gh_sp_data_2019,Gh_env_data_2019$FireGrzTrt)
simper_gh_2019_sum<-summary(simper_gh_2019, order=T)
Simper_df<-data.frame(simper_gh_2019_sum$ABG_PBG)
#2020
Gh_sp_data_2020 <- grassh_permav_df %>%
  filter(RecYear==2020)%>%
  ungroup()%>%
  dplyr::select(-1:-6)
Gh_env_data_2020 <- grassh_permav_df%>%
  dplyr::select(1:6)%>%
  filter(RecYear==2020)

simper_gh_2020<-simper(Gh_sp_data_2020,Gh_env_data_2020$FireGrzTrt)
simper_gh_2020_sum<-summary(simper_gh_2020, order=T)
Simper_df<-data.frame(simper_gh_2020_sum$ABG_PBG)
#2021
Gh_sp_data_2021 <- grassh_permav_df %>%
  filter(RecYear==2021)%>%
  ungroup()%>%
  dplyr::select(-1:-6)
Gh_env_data_2021 <- grassh_permav_df%>%
  dplyr::select(1:6)%>%
  filter(RecYear==2021)

simper_2021<-simper(Gh_sp_data_2021,Gh_env_data_2021$FireGrzTrt)
simper_gh_2021<-summary(simper_2021, order=T)
Simper_df_2021<-data.frame(simper_gh_2021$ABG_PBG)




#extract species driving difference in composition####
grassh_cover_2021_df <- grasshopperspcomp_df%>%
  full_join(watershed_key, by = "Watershed")%>%
  left_join(Watershed_key2,by="Watershed")%>%
  mutate(spe=paste(Genus,Species, sep="_"))%>%
  mutate(RecYear = Recyear)%>%
  filter(!spe%in%c("Oecanthinae_spp.","Tettigoniidae_spp.","Gryllidae_spp.",
                   "Conocephalus_spp.","Neoconocephalus_robustus","Scudderia_texensis",
                   "Arethaea_constricta","Orchelimum_spp.","Amblycorypha_oblongifolia","Pediodectes_haldemani",
                   "Amblycorypha_rotundifolia","Neoconocephalus_spp.","Neoconocephalus_ensiger","Pediodectes_nigromarginatus",
                   "Scudderia_furcata","Scudderia_spp."))%>%
  #using the max cover from the two sweeps done on each transect
  group_by(Unit,RecYear,FireGrzTrt,Watershed,Repsite,spe)%>%
  summarise(Total=max(Total))%>%
  filter(Repsite!="C" | Watershed!="C03C")%>%
  filter(Watershed!="C03C" | Repsite!="D")%>%
  filter(Watershed!="C03B" | Repsite!="A")%>%
  filter(Watershed!="C03B" | Repsite!="B")%>%
  filter(Watershed!="C03B" | Repsite!="C")%>%
  filter(Watershed!="C03A" | Repsite!="A")%>%
  filter(Watershed!="C03A" | Repsite!="B")%>%
  filter(Watershed!="C03A" | Repsite!="C")%>%
  filter(Watershed!="C3SC" | Repsite!="A")%>%
  filter(Watershed!="C3SC" | Repsite!="B")%>%
  filter(Watershed!="C3SC" | Repsite!="C")%>%
  filter(Watershed!="C3SA" | Repsite!="B")%>%
  filter(Watershed!="C3SA" | Repsite!="C")%>%
  filter(Watershed!="C3SA" | Repsite!="D")%>%
  filter(Watershed!="C3SB" | Repsite!="A")%>%
  filter(Watershed!="C3SB" | Repsite!="B")%>%
  group_by(Unit,RecYear,Watershed,Repsite,FireGrzTrt,spe)%>%
  summarise(total_avg=mean(Total,na.rm=T))%>%
  group_by(Unit,RecYear,Watershed,Repsite,FireGrzTrt)%>%
  #converting count data to abundance data (0-100% scale)
  mutate(abundance=(total_avg/sum(total_avg, na.rm=T))*100,
         unit_trt=paste(Unit, FireGrzTrt, sep="_"))%>%
  group_by(Unit,RecYear,Watershed,Repsite,FireGrzTrt,unit_trt,spe)%>%
  summarise(abundance=mean(abundance,na.rm=T))%>%
  filter(spe%in%c("Orphulella_speciosa",
                  "Phoetaliotes_nebrascensis","Melanoplus_keeleri",
                  "Eritettix_simplex","Melanoplus_femurrubrum","Campylacantha_olivacea",
                  "Hypochlora_alba","Mermiria_bivittata","Oedipodinae_spp.","Melanoplus_scudderi",
                  "Mermiria_spp.","Syrbula_admirabilis"))#%>%
#Checking trends for each species
ggplot_species_data<-grassh_cover_2021_df%>%
  filter(!RecYear==2010)%>%
  group_by(FireGrzTrt,RecYear, spe)%>%
  summarise(abundance=mean(abundance))
ggplot(ggplot_species_data, aes(RecYear, abundance, col=FireGrzTrt))+
  geom_point()+
  geom_path()+
  facet_wrap(~spe, scales = "free")

gh_cover_2021<-grassh_cover_2021_df%>%
  filter(RecYear==2021)%>%
  group_by(spe,FireGrzTrt)%>%
  summarise(abundance=mean(abundance,na.rm=F))%>%
  pivot_wider(names_from = "FireGrzTrt", values_from = "abundance", values_fill = 0)%>%
  mutate(abund_ABG_PBG=ABG-PBG)
#ggplot 
#reset theme for this figure
#theme_set(theme_bw())
#theme_update(axis.title.x=element_text(size=10, vjust=-0.35), axis.text.x=element_text(size=10),
#             axis.title.y=element_text(size=10, angle=90, vjust=0.5), axis.text.y=element_text(size=10),
#             plot.title = element_text(size=14, vjust=2),
#             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
#             legend.title=element_text(size=10), legend.text=element_text(size=10))
ggplot(gh_cover_2021,aes(x=spe))+
  geom_bar(stat = "identity",aes(y=abund_ABG_PBG),width = 0.5)+
  scale_x_discrete(guide = guide_axis(angle = 90))



#just checking if any species is an indicator distingushing treatent
#library(indicspecies)
#indicator_abgvspbg<-multipatt(Gh_sp_data_2011,Gh_env_data_2011$FireGrzTrt)
#ff<-summary(indicator_abgvspbg)
#none

