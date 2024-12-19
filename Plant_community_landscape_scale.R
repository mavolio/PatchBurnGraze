#Patch-Burn Synthesis Project
#Plant community data at the landscape scale
#Author: Joshua Adedayo Ajowele joshuaajowele@gmail.com
#Started: May 13, 2024 last modified: Nov 27, 2024

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
library(lmerTest)
library(see)
library(patchwork)
library(phia)
library(readxl)

## Set graphing parameters
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=20), legend.text=element_text(size=20))
### Read in raw data

#import data as dataframe
species_comp_data<- read.csv("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/PhD Wyoming_One Drive/PHD Wyoming/Thesis/PBG synthesis/PBG_plant_comp_2008_2022.csv")

#Creating a key for converting cover class to abundance 
cover_key <-data.frame(CoverClass=0:7, abundance=c(0,0.5,3,15,37.5,62.5,85,97.5))
#need to make transect all uppercase
species_comp_data$Transect<-toupper(species_comp_data$Transect)
#make all genera lowercase for consistency
species_comp_data$Ab_genus<-tolower(species_comp_data$Ab_genus)
species_comp_data$Ab_species<-tolower(species_comp_data$Ab_species)
#merging key with data
species_comp<-species_comp_data%>%
  left_join(cover_key, by="CoverClass")%>%
  #create unique name for species by combining binomial nomenclature
  mutate(sp=paste(Ab_genus,Ab_species, sep="_"))%>%
  group_by(Watershed, RecYear,Transect,Plot,sp)%>%
  #selecting the maximum cover for each species
  summarise(abundance=max(abundance, na.rm=T))%>%
  #removing unwanted years
  filter(!RecYear%in%2008:2011)%>%
  filter(RecYear!=2022)

#Create a watershed key column to merge with the raw data
watershed_key <- data.frame(Watershed=levels(factor(species_comp$Watershed)),# create key to merge with raw data
                            FireGrzTrt=c("ABG", "PBG", "PBG", "PBG", "ABG", "PBG", "PBG", "PBG"))
Watershed_key2<-tibble(Watershed=levels(factor(species_comp$Watershed)),# create key to merge with raw data
                       Unit=c("south", "south", "south", "south", "north", "north",
                              "north", "north"))
#merging key with dataset
sp_comp<-species_comp%>%
  left_join(watershed_key,by="Watershed")%>%
  left_join(Watershed_key2,by="Watershed")%>%
  filter(Transect!="C" | Watershed!="C03C")%>%
  filter(Watershed!="C03C" | Transect!="D")%>%
  filter(Watershed!="C03B" | Transect!="A")%>%
  filter(Watershed!="C03B" | Transect!="B")%>%
  filter(Watershed!="C03B" | Transect!="C")%>%
  filter(Watershed!="C03A" | Transect!="A")%>%
  filter(Watershed!="C03A" | Transect!="B")%>%
  filter(Watershed!="C03A" | Transect!="C")%>%
  filter(Watershed!="C3SC" | Transect!="A")%>%
  filter(Watershed!="C3SC" | Transect!="B")%>%
  filter(Watershed!="C3SC" | Transect!="C")%>%
  filter(Watershed!="C3SA" | Transect!="B")%>%
  filter(Watershed!="C3SA" | Transect!="C")%>%
  filter(Watershed!="C3SA" | Transect!="D")%>%
  filter(Watershed!="C3SB" | Transect!="A")%>%
  filter(Watershed!="C3SB" | Transect!="B")%>%
  group_by(Unit,RecYear,FireGrzTrt,sp)%>%
  summarise(abundance_avg=mean(abundance,na.rm=T))%>%
  group_by(Unit,RecYear,FireGrzTrt)%>%
  #abundance data (0-100% scale)
  mutate(abundance=(abundance_avg/sum(abundance_avg, na.rm=T))*100)%>%
  group_by(Unit,RecYear,FireGrzTrt,sp)%>%
  mutate(Rep_id=paste(Unit,FireGrzTrt, sep="_"))
  
#community metrics with codyn####
#deriving richness and evenness using codyn at landscape scale
pl_rich <- community_structure(sp_comp, time.var = "RecYear", 
                               abundance.var = "abundance",
                               replicate.var = "Rep_id", metric = "Evar")

#deriving diversity values
pl_diver <- community_diversity(sp_comp, time.var = "RecYear", 
                                abundance.var = "abundance",
                                replicate.var = "Rep_id", metric="Shannon")


#extracting important columns 
sp_comp_subset <- sp_comp %>%
  ungroup()%>%
  dplyr::select(RecYear,FireGrzTrt,Rep_id, Unit) %>%
  distinct()#remove repeated rows

#combine into a single data
pl_rich_diver<- pl_diver %>%
  #combine richness and diversity
  left_join(pl_rich, by = c("Rep_id","RecYear")) %>%
  #combine with important columns
  left_join(sp_comp_subset, by = c("Rep_id","RecYear"))
#convert to factors
pl_rich_diver$RecYear<-as.factor(pl_rich_diver$RecYear)
pl_rich_diver$Unit<-as.factor(pl_rich_diver$Unit)
pl_rich_diver$FireGrzTrt<-as.factor(pl_rich_diver$FireGrzTrt)
#build model
#remove diversity code from script
pl_diver_model<-lmer(log(Shannon)~FireGrzTrt*RecYear+(1|Unit),
                     data=pl_rich_diver)
Anova(pl_diver_model, type=3)
anova(pl_diver_model)
qqnorm(resid(pl_diver_model))
check_model(pl_diver_model)
summary(pl_diver_model)

pl_diver_l_model<-lm(log(Shannon)~FireGrzTrt*RecYear+Unit, data=pl_rich_diver)
summary(pl_diver_l_model)
anova(pl_diver_l_model)#similar to mixed model.
#multiple comparison
testInteractions(pl_diver_model, pairwise="FireGrzTrt", fixed="RecYear",
                 adjustment="BH")
#using mean estimate to create figure 
pldiver_interact<-interactionMeans(pl_diver_model)
#replacing spaces in column names with underscore 
names(pldiver_interact)<-str_replace_all(names(pldiver_interact), " ","_")
#df for visuals from model estimates
pldiver_interact_viz<-pldiver_interact%>%
  mutate(pldiver_mean=exp(adjusted_mean),
         pldiver_upper=exp(adjusted_mean+SE_of_link),
         pldiver_lower=exp(adjusted_mean-SE_of_link))
#visual
ggplot(pldiver_interact_viz,aes(RecYear, pldiver_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=pldiver_lower,
                    ymax=pldiver_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))
#average across years for simplification
pldiver_interact_bar<-pldiver_interact_viz%>%
  group_by(FireGrzTrt)%>%
  summarise(pl_diver=mean(pldiver_mean),
            se_upper=mean(pldiver_upper),
            se_lower=mean(pldiver_lower))
ggplot(pldiver_interact_bar,aes(x=FireGrzTrt,fill=FireGrzTrt))+
  geom_bar(stat = "identity",aes(y=pl_diver),width = 0.5)+
  geom_errorbar(aes(ymin=se_lower,
                    ymax=se_upper),width=0.2,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"))
#richness model
plrich_model<-lmer(log(richness)~FireGrzTrt*RecYear+(1|Unit),
                   data=pl_rich_diver)
Anova(plrich_model, type=3)
anova(plrich_model)
qqnorm(resid(plrich_model))
check_model(plrich_model)
#multiple comparison
testInteractions(plrich_model, pairwise="FireGrzTrt", fixed="RecYear",
                 adjustment="BH")
#using mean estimate to create figure 
plrich_interact<-interactionMeans(plrich_model)
#replacing spaces in column names with underscore 
names(plrich_interact)<-str_replace_all(names(plrich_interact), " ","_")
#df for visuals from model estimates
plrich_interact_viz<-plrich_interact%>%
  mutate(plrich_mean=exp(adjusted_mean),
         plrich_upper=exp(adjusted_mean+SE_of_link),
         plrich_lower=exp(adjusted_mean-SE_of_link))
#visual
ggplot(plrich_interact_viz,aes(RecYear, plrich_mean, col=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=plrich_lower,
                    ymax=plrich_upper),width=0.1,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

#average across years for simplification
plrich_interact_bar<-plrich_interact_viz%>%
  group_by(FireGrzTrt)%>%
  summarise(pl_richness=mean(plrich_mean),
            se_upper=mean(plrich_upper),
            se_lower=mean(plrich_lower))
ggplot(plrich_interact_bar,aes(x=FireGrzTrt,fill=FireGrzTrt))+
  geom_bar(stat = "identity",aes(y=pl_richness),width = 0.25)+
  geom_errorbar(aes(ymin=se_lower,
                    ymax=se_upper),width=0.05,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"))

#evenness model
plevenness_model<-lmer(Evar~FireGrzTrt*RecYear+(1|Unit),
                       data=pl_rich_diver)#is.singular
Anova(plevenness_model, type=3)
anova(plevenness_model)
qqnorm(resid(plevenness_model))
check_model(plevenness_model)
summary(plevenness_model)

#multiple comparison
#testInteractions(plevenness_model, pairwise="FireGrzTrt", fixed="RecYear",
#                 adjustment="BH")
#using mean estimate to create figure 
pleven_interact<-interactionMeans(plevenness_model)
#replacing spaces in column names with underscore 
names(pleven_interact)<-str_replace_all(names(pleven_interact), " ","_")
#df for visuals from model estimates
pleven_interact_viz<-pleven_interact%>%
  mutate(plevenness_mean=adjusted_mean,
         pleven_upper=adjusted_mean+SE_of_link,
         pleven_lower=adjusted_mean-SE_of_link)
#visual
ggplot(pleven_interact_viz,aes(RecYear, plevenness_mean, col=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=pleven_lower,
                    ymax=pleven_upper),width=0.1,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

#average across years for simplification
pleven_interact_bar<-pleven_interact_viz%>%
  group_by(FireGrzTrt)%>%
  summarise(pl_evenness=mean(plevenness_mean),
            se_upper=mean(pleven_upper),
            se_lower=mean(pleven_lower))
ggplot(pleven_interact_bar,aes(x=FireGrzTrt,fill=FireGrzTrt))+
  geom_bar(stat = "identity",aes(y=pl_evenness),width = 0.25)+
  geom_errorbar(aes(ymin=se_lower,
                    ymax=se_upper),width=0.05,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"))


#community composition####
sp_comp_comm<-species_comp%>%
  left_join(watershed_key,by="Watershed")%>%
  left_join(Watershed_key2,by="Watershed")%>%
  filter(Transect!="C" | Watershed!="C03C")%>%
  filter(Watershed!="C03C" | Transect!="D")%>%
  filter(Watershed!="C03B" | Transect!="A")%>%
  filter(Watershed!="C03B" | Transect!="B")%>%
  filter(Watershed!="C03B" | Transect!="C")%>%
  filter(Watershed!="C03A" | Transect!="A")%>%
  filter(Watershed!="C03A" | Transect!="B")%>%
  filter(Watershed!="C03A" | Transect!="C")%>%
  filter(Watershed!="C3SC" | Transect!="A")%>%
  filter(Watershed!="C3SC" | Transect!="B")%>%
  filter(Watershed!="C3SC" | Transect!="C")%>%
  filter(Watershed!="C3SA" | Transect!="B")%>%
  filter(Watershed!="C3SA" | Transect!="C")%>%
  filter(Watershed!="C3SA" | Transect!="D")%>%
  filter(Watershed!="C3SB" | Transect!="A")%>%
  filter(Watershed!="C3SB" | Transect!="B")%>%
  group_by(Unit,RecYear,FireGrzTrt,Transect,Plot,Watershed,sp)%>%
  summarise(abundance_avg=mean(abundance,na.rm=T))%>%
  group_by(Unit,RecYear,FireGrzTrt,Transect,Plot,Watershed)%>%
  #abundance data (0-100% scale)
  mutate(abundance=(abundance_avg/sum(abundance_avg, na.rm=T))*100,
         unit_trt=paste(Unit,FireGrzTrt, sep="_"))%>%
  group_by(Unit,RecYear,FireGrzTrt,Transect,Plot,Watershed,unit_trt,sp)%>%
  summarise(abundance=mean(abundance,na.rm=T))%>%
  pivot_wider(names_from = sp, values_from = abundance, values_fill = 0)

# Separate out spcomp and environmental columns (cols are species) #
pl_sp_data <- sp_comp_comm %>%
  ungroup()%>%
  dplyr::select(-1:-7)
pl_env_data <- sp_comp_comm%>%dplyr::select(1:7)

#get nmds1 and 2
pl_mds_all <- metaMDS(pl_sp_data, distance = "bray")
#stress 0.275

#combine NMDS1 and 2 with factor columns and create centroids
pl_mds_scores <- data.frame(pl_env_data, scores(pl_mds_all, display="sites"))%>%
  group_by(RecYear, FireGrzTrt)%>%
  mutate(NMDS1_mean=mean(NMDS1),
         NMDS2_mean=mean(NMDS2))
#write into csv format 
#write.csv(pl_mds_scores,"C:/Users/joshu/OneDrive - UNCG/UNCG PHD/PhD Wyoming_One Drive/PHD Wyoming/Thesis/PBG synthesis/pl_mds_scores.csv")

#for ordispider sake!
pl_mds_scores_mean <- pl_mds_scores%>%
  group_by(RecYear, FireGrzTrt)%>%
  summarise(NMDS1=mean(NMDS1),
            NMDS2=mean(NMDS2))
#ordispider with ggplot
ggplot(pl_mds_scores, aes(x=NMDS1, y=NMDS2, col=FireGrzTrt)) +
  geom_segment(aes(xend=NMDS1_mean, yend= NMDS2_mean))+
  geom_point(data=pl_mds_scores_mean, size=5) +
  geom_point()+
  facet_wrap(~ RecYear, scales = "free")+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))  #+
theme_bw() 

#summarise into centroid and have a single figure
#Plot with geompoint for abg vs pbg nmds
ggplot(pl_mds_scores_mean, aes(x=NMDS1, y=NMDS2, col=FireGrzTrt))+
  geom_point(size=5)+
  geom_path()+
  #geom_errorbar(aes(ymax=NMDS2_mean+NMDS2_SE, ymin=NMDS2_mean-NMDS2_SE))+
  #geom_errorbarh(aes(xmax=NMDS1_mean+NMDS1_SE, xmin=NMDS1_mean-NMDS1_SE))+
  scale_shape_manual(values=10:21)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

#not necessary-just checking if composition differs with unit####
pl_mds_scores_unit <- data.frame(pl_env_data, scores(pl_mds_all, display="sites"))%>%
  group_by(RecYear, unit_trt)%>%
  mutate(NMDS1_mean=mean(NMDS1),
         NMDS2_mean=mean(NMDS2))

#for ordispider sake!
pl_mds_scores_mean_unit <- pl_mds_scores_unit%>%
  group_by(RecYear, unit_trt)%>%
  summarise(NMDS1=mean(NMDS1),
            NMDS2=mean(NMDS2))
#ordispider with ggplot
ggplot(pl_mds_scores_unit, aes(x=NMDS1, y=NMDS2, col=unit_trt)) +
  geom_segment(aes(xend=NMDS1_mean, yend= NMDS2_mean))+
  geom_point(data=pl_mds_scores_mean_unit, size=5) +
  geom_point()+
  facet_wrap(~ RecYear, scales = "free")#+
scale_colour_manual(values=c( "#F0E442", "#009E73"))  #+
theme_bw() 

#permanova and betdiversity####
#calculating permanova and betadiversity
#creating a loop to do this
year_vec_pl <- unique(pl_env_data$RecYear)
pl_perm <- {}
pl_beta <- {}


for(YEAR in 1:length(year_vec_pl)){
  vdist_temp_pl <- vegdist(filter(pl_sp_data, pl_env_data$RecYear ==  year_vec_pl[YEAR]))
  permanova_temp_pl <- adonis(vdist_temp_pl ~ subset(pl_env_data, RecYear == year_vec_pl[YEAR])$FireGrzTrt)
  permanova_out_temp_pl <- data.frame(RecYear = year_vec_pl[YEAR], 
                                      DF = permanova_temp_pl$aov.tab[1,1],
                                      F_value = permanova_temp_pl$aov.tab[1,4],
                                      P_value = permanova_temp_pl$aov.tab[1,6])
  pl_perm <- rbind(pl_perm,permanova_out_temp_pl)
  
  bdisp_temp_pl <- betadisper(vdist_temp_pl, filter(pl_env_data, RecYear==year_vec_pl[YEAR])$FireGrzTrt, type = "centroid")
  bdisp_out_temp_pl <- data.frame(filter(pl_env_data, RecYear==year_vec_pl[YEAR]), distance = bdisp_temp_pl$distances)
  pl_beta <- rbind(pl_beta, bdisp_out_temp_pl)
  
  rm(vdist_temp_pl, permanova_temp_pl, permanova_out_temp_pl, bdisp_temp_pl, bdisp_out_temp_pl)
}
write.csv(pl_beta, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/plant_betadiver.csv")
write.csv(pl_perm, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/plant_permanova.csv")
#model for betadiversity

#calculating betadiversity by unit####
#creating a loop to do this
year_vec_pl <- unique(pl_env_data$RecYear)
pl_perm_unit <- {}
pl_beta_unit <- {}


for(YEAR in 1:length(year_vec_pl)){
  vdist_temp_pl_unit <- vegdist(filter(pl_sp_data, pl_env_data$RecYear ==  year_vec_pl[YEAR]))
  permanova_temp_pl_unit <- adonis(vdist_temp_pl_unit ~ subset(pl_env_data, RecYear == year_vec_pl[YEAR])$unit_trt)
  permanova_out_temp_pl_unit <- data.frame(RecYear = year_vec_pl[YEAR], 
                                      DF = permanova_temp_pl_unit$aov.tab[1,1],
                                      F_value = permanova_temp_pl_unit$aov.tab[1,4],
                                      P_value = permanova_temp_pl_unit$aov.tab[1,6])
  pl_perm_unit <- rbind(pl_perm_unit,permanova_out_temp_pl_unit)
  
  bdisp_temp_pl_unit <- betadisper(vdist_temp_pl_unit, filter(pl_env_data, RecYear==year_vec_pl[YEAR])$unit_trt, type = "centroid")
  bdisp_out_temp_pl_unit <- data.frame(filter(pl_env_data, RecYear==year_vec_pl[YEAR]), distance = bdisp_temp_pl_unit$distances)
  pl_beta_unit <- rbind(pl_beta_unit, bdisp_out_temp_pl_unit)
  
  rm(vdist_temp_pl_unit, permanova_temp_pl_unit, permanova_out_temp_pl_unit, bdisp_temp_pl_unit, bdisp_out_temp_pl_unit)
}
write.csv(pl_beta_unit, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/plant_betadiver_sep_unit.csv")

#model for betadiversity
pl_beta_unit$RecYear<-as.factor(pl_beta_unit$RecYear)
pl_beta_model<-lmer(log(distance)~FireGrzTrt*RecYear+(1|Unit),
                    data=pl_beta_unit)
Anova(pl_beta_model, type=3)
anova(pl_beta_model)
qqnorm(resid(pl_beta_model))
check_model(pl_beta_model)
#pairwise interaction 
#testInteractions(pl_beta_model, fixed="RecYear",
#                 pairwise = "FireGrzTrt", adjustment="BH")
#using mean estimate from post-hoc to create figure as a comparison to raw data
model_estimates_beta<-interactionMeans(pl_beta_model)
#replacing spaces in column names with underscore 
names(model_estimates_beta)<-str_replace_all(names(model_estimates_beta), " ","_")
#df for visuals from model estimates
model_estimates_beta_viz<-model_estimates_beta%>%
  mutate(pl_betadiv=exp(adjusted_mean),
         pl_upper=exp(adjusted_mean+SE_of_link),
         pl_lower=exp(adjusted_mean-SE_of_link))
#visual
ggplot(model_estimates_beta_viz,aes(RecYear, pl_betadiv, col=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=pl_lower,
                    ymax=pl_upper),width=0.1,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))
#summarize with a bargraph
plbeta_interact_bar<-model_estimates_beta_viz%>%
  group_by(FireGrzTrt)%>%
  summarise(pl_betadiver=mean(pl_betadiv),
            se_upper=mean(pl_upper),
            se_lower=mean(pl_lower))
ggplot(plbeta_interact_bar,aes(x=FireGrzTrt,fill=FireGrzTrt))+
  geom_bar(stat = "identity",aes(y=pl_betadiver),width = 0.25)+
  geom_errorbar(aes(ymin=se_lower,
                    ymax=se_upper),width=0.05,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"))

#simper analysis to determine species contributing 80% of difference in composition####
#2021
pl_sp_data_2021 <- sp_comp_comm %>%
  filter(RecYear==2021)%>%
  ungroup()%>%
  dplyr::select(-1:-7)
pl_env_data_2021 <- sp_comp_comm%>%
  dplyr::select(1:7)%>%
  filter(RecYear==2021)

simper_2021<-simper(pl_sp_data_2021,pl_env_data_2021$FireGrzTrt)
simper_pl_2021<-summary(simper_2021, order=T)
Simper_df_2021<-data.frame(simper_pl_2021$ABG_PBG)


sp_comp_comm_2011_2021<-species_comp%>%
  left_join(watershed_key,by="Watershed")%>%
  left_join(Watershed_key2,by="Watershed")%>%
  filter(Transect!="C" | Watershed!="C03C")%>%
  filter(Watershed!="C03C" | Transect!="D")%>%
  filter(Watershed!="C03B" | Transect!="A")%>%
  filter(Watershed!="C03B" | Transect!="B")%>%
  filter(Watershed!="C03B" | Transect!="C")%>%
  filter(Watershed!="C03A" | Transect!="A")%>%
  filter(Watershed!="C03A" | Transect!="B")%>%
  filter(Watershed!="C03A" | Transect!="C")%>%
  filter(Watershed!="C3SC" | Transect!="A")%>%
  filter(Watershed!="C3SC" | Transect!="B")%>%
  filter(Watershed!="C3SC" | Transect!="C")%>%
  filter(Watershed!="C3SA" | Transect!="B")%>%
  filter(Watershed!="C3SA" | Transect!="C")%>%
  filter(Watershed!="C3SA" | Transect!="D")%>%
  filter(Watershed!="C3SB" | Transect!="A")%>%
  filter(Watershed!="C3SB" | Transect!="B")%>%
  group_by(Unit,RecYear,FireGrzTrt,Transect,Plot,Watershed,sp)%>%
  summarise(abundance_avg=mean(abundance,na.rm=T))%>%
  group_by(Unit,RecYear,FireGrzTrt,Transect,Plot,Watershed)%>%
  #abundance data (0-100% scale)
  mutate(abundance=(abundance_avg/sum(abundance_avg, na.rm=T))*100,
         unit_trt=paste(Unit,FireGrzTrt, sep="_"))%>%
  group_by(Unit,RecYear,FireGrzTrt,Transect,Plot,Watershed,unit_trt,sp)%>%
  summarise(abundance=mean(abundance,na.rm=T))%>%
  filter(sp%in%c("androp_gerar","schiza_scopa","amorph_canes","symphy_erico","sorpla_nutan",
                 "koeler_macra","carex_inops","ambros_psilo","solida_misso","artemi_ludov","poa_prate",
                 "rosa_arkan","symphy_oblon","boutel_curti","solida_rigid","antenn_negle","rhus_glabr",
                 "solida_speci","comand_umbel","lesped_viola","dalea_multi","salvia_azure"))
#Checking trends for each species
ggplot_species_data<-sp_comp_comm_2011_2021%>%
  group_by(FireGrzTrt,RecYear, sp)%>%
  summarise(abundance=mean(abundance))
ggplot(ggplot_species_data, aes(RecYear, abundance, col=FireGrzTrt))+
  geom_point()+
  geom_path()+
  facet_wrap(~sp, scales = "free")

pl_cover_2021<-sp_comp_comm_2011_2021%>%
  filter(RecYear==2021)%>%
  group_by(sp,FireGrzTrt)%>%
  summarise(abundance=mean(abundance,na.rm=F))%>%
  pivot_wider(names_from = "FireGrzTrt", values_from = "abundance", values_fill = 0)%>%
  mutate(abund_ABG_PBG=ABG-PBG)
#ggplot 
#reset theme for this figure
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=10, vjust=-0.35), axis.text.x=element_text(size=10),
             axis.title.y=element_text(size=10, angle=90, vjust=0.5), axis.text.y=element_text(size=10),
             plot.title = element_text(size=14, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=10), legend.text=element_text(size=10))
#reorder species by descending order
ggplot(pl_cover_2021,aes(x=reorder(sp,-abund_ABG_PBG)))+
  geom_bar(stat = "identity",aes(y=abund_ABG_PBG),width = 0.5)+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  ylab("Difference in abundance (ABG-PBG)")+
  xlab("Species")



#just checking if any species is an indicator distinguishing treatment####
library(indicspecies)
indicator_abgvspbg<-multipatt(pl_sp_data_2021,pl_env_data_2021$FireGrzTrt)
ff<-summary(indicator_abgvspbg)
#There are!


###community composition based on year since fire and watershed##
#creating a key for year since fire####
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


#merging key with dataset####
burn_time_sp_comp<-species_comp%>%
  filter(RecYear%in%2016:2021)%>%
  left_join(watershed_key,by="Watershed")%>%
  left_join(Watershed_key2,by="Watershed")%>%
  mutate(year_watershed=paste(RecYear,Watershed,sep="_"))%>%
  left_join(YrSinceFire_key, by="year_watershed")%>%
  group_by(RecYear,yrsins_fire,Unit,Watershed,FireGrzTrt,Transect,sp)%>%
  summarise(abundance=mean(abundance, na.rm=T))%>%
  mutate(abundance=abundance/100)%>%
  pivot_wider(names_from = sp, values_from = abundance, values_fill = 0)
  

# Separate out spcomp and environmental columns (cols are species) #
burn_time_sp_data <- burn_time_sp_comp %>%
  ungroup()%>%
  dplyr::select(-1:-6)
burn_time_env_data <- burn_time_sp_comp%>%dplyr::select(1:6)

#get nmds1 and 2
burn_time_mds_all <- metaMDS(burn_time_sp_data, distance = "bray")
#Run 20 stress 0.211

#combine NMDS1 and 2 with factor columns and create centroids
burn_time_mds_scores <- data.frame(burn_time_env_data, scores(burn_time_mds_all, display="sites"))%>%
  group_by(RecYear, Unit,yrsins_fire,Watershed)%>%
  mutate(NMDS1_mean=mean(NMDS1),
         NMDS2_mean=mean(NMDS2))

#plotting centroid through time
ggplot(burn_time_mds_scores, aes(x=NMDS1_mean, y=NMDS2_mean, col=yrsins_fire, shape=Watershed))+
  geom_point(size=8)+
  geom_path()+
  scale_shape_manual(values=c(15:18,0:2,5))+
  scale_colour_manual(values=c("#F0E442", "#994F00", "#999999", "#0072B2"))#+
  facet_wrap(~Unit, scales="free")

#calculating permanova
#creating a loop to do this
burn_year_vec_pl <- unique(burn_time_env_data$RecYear)
burn_pl_perm <- {}

#by Year since fire
for(YEAR in 1:length(burn_year_vec_pl)){
  burn_vdist_temp_pl <- vegdist(filter(burn_time_sp_data, burn_time_env_data$RecYear ==  burn_year_vec_pl[YEAR]))
  burn_permanova_temp_pl <- adonis(burn_vdist_temp_pl ~ subset(burn_time_env_data, RecYear == burn_year_vec_pl[YEAR])$yrsins_fire)
  burn_permanova_out_temp_pl <- data.frame(RecYear = burn_year_vec_pl[YEAR], 
                                      DF = burn_permanova_temp_pl$aov.tab[1,1],
                                      F_value = burn_permanova_temp_pl$aov.tab[1,4],
                                      P_value = burn_permanova_temp_pl$aov.tab[1,6])
  burn_pl_perm <- rbind(burn_pl_perm,burn_permanova_out_temp_pl)
  rm(burn_vdist_temp_pl, burn_permanova_temp_pl, burn_permanova_out_temp_pl)
}

write.csv(burn_pl_perm, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/plant_permanova_yearsincefire.csv")

#multiple comparison for significant years 2018 and 2020
burn_time_sp_data_2018<-burn_time_sp_data%>%
  filter(burn_time_env_data$RecYear==2018)
burn_time_env_data_2018<-burn_time_env_data%>%
  filter(RecYear==2018)
#package for pERMANOVA pairwise comparison 
#library(devtools)
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
pairwise.adonis2(burn_time_sp_data_2018~yrsins_fire, data=burn_time_env_data_2018)


#compare by watershed####
burn_year_vec_pl <- unique(burn_time_env_data$RecYear)
burn_pl_perm_watershed <- {}

for(YEAR in 1:length(burn_year_vec_pl)){
  burn_vdist_temp_pl <- vegdist(filter(burn_time_sp_data, burn_time_env_data$RecYear ==  burn_year_vec_pl[YEAR]))
  burn_permanova_temp_pl <- adonis(burn_vdist_temp_pl ~ subset(burn_time_env_data, RecYear == burn_year_vec_pl[YEAR])$Watershed)
  burn_permanova_out_temp_pl <- data.frame(RecYear = burn_year_vec_pl[YEAR], 
                                           DF = burn_permanova_temp_pl$aov.tab[1,1],
                                           F_value = burn_permanova_temp_pl$aov.tab[1,4],
                                           P_value = burn_permanova_temp_pl$aov.tab[1,6])
  burn_pl_perm_watershed <- rbind(burn_pl_perm_watershed,burn_permanova_out_temp_pl)
  rm(burn_vdist_temp_pl, burn_permanova_temp_pl, burn_permanova_out_temp_pl)
}

write.csv(burn_pl_perm_watershed, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/plant_permanova_watershed.csv")

#seems magnitude of diff in Unit>Watershed>Yrsins_fire
#comparing watershed and time since burn within the same unit
#merging key with dataset####
burn_time_sp_comp_south<-species_comp%>%
  filter(RecYear%in%2016:2021)%>%
  left_join(watershed_key,by="Watershed")%>%
  left_join(Watershed_key2,by="Watershed")%>%
  mutate(year_watershed=paste(RecYear,Watershed,sep="_"))%>%
  left_join(YrSinceFire_key, by="year_watershed")%>%
  filter(Unit=="south")%>%
  group_by(RecYear,yrsins_fire,Unit,Watershed,FireGrzTrt,Transect,sp)%>%
  summarise(abundance=mean(abundance, na.rm=T))%>%
  mutate(abundance=abundance/100)%>%
  pivot_wider(names_from = sp, values_from = abundance, values_fill = 0)


# Separate out spcomp and environmental columns (cols are species) #
burn_time_sp_data_south <- burn_time_sp_comp_south %>%
  ungroup()%>%
  dplyr::select(-1:-6)
burn_time_env_data_south <- burn_time_sp_comp_south%>%dplyr::select(1:6)
#compare by watershed
burn_year_vec_pl_s <- unique(burn_time_env_data_south$RecYear)
burn_pl_perm_watershed_s <- {}

for(YEAR in 1:length(burn_year_vec_pl_s)){
  burn_vdist_temp_pl <- vegdist(filter(burn_time_sp_data_south, burn_time_env_data_south$RecYear ==  burn_year_vec_pl_s[YEAR]))
  burn_permanova_temp_pl <- adonis(burn_vdist_temp_pl ~ subset(burn_time_env_data_south, RecYear == burn_year_vec_pl_s[YEAR])$yrsins_fire)
  burn_permanova_out_temp_pl <- data.frame(RecYear = burn_year_vec_pl_s[YEAR], 
                                           DF = burn_permanova_temp_pl$aov.tab[1,1],
                                           F_value = burn_permanova_temp_pl$aov.tab[1,4],
                                           P_value = burn_permanova_temp_pl$aov.tab[1,6])
  burn_pl_perm_watershed_s <- rbind(burn_pl_perm_watershed_s,burn_permanova_out_temp_pl)
  rm(burn_vdist_temp_pl, burn_permanova_temp_pl, burn_permanova_out_temp_pl)
}

write.csv(burn_pl_perm_watershed, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/plant_permanova_watershed.csv")

#subsetting by year
burn_time_env_data_south_2016<-burn_time_env_data_south%>%
  filter(RecYear==2016)
burn_time_sp_data_south_2016<-burn_time_sp_data_south%>%
  filter(burn_time_env_data_south$RecYear==2016)

#package for pERMANOVA pairwise comparison 
pairwise.adonis2(burn_time_sp_data_south_2016~yrsins_fire, data=burn_time_env_data_south_2016)


#comparing watershed and year since fire treatment without considering year
dist_south<-vegdist(burn_time_sp_data_south)
permanova_south<-adonis(dist_south~burn_time_env_data_south$yrsins_fire+burn_time_env_data_south$Watershed+as.factor(burn_time_env_data_south$RecYear))
permanova_south$aov.tab

pairwise.adonis2(burn_time_sp_data_south~yrsins_fire+Watershed+as.factor(RecYear), data=burn_time_env_data_south)
pairwise.adonis2(burn_time_sp_data_south~Watershed+yrsins_fire+as.factor(RecYear), data=burn_time_env_data_south)

#north unit
burn_time_sp_comp_north<-species_comp%>%
  filter(RecYear%in%2016:2021)%>%
  left_join(watershed_key,by="Watershed")%>%
  left_join(Watershed_key2,by="Watershed")%>%
  mutate(year_watershed=paste(RecYear,Watershed,sep="_"))%>%
  left_join(YrSinceFire_key, by="year_watershed")%>%
  filter(Unit=="north")%>%
  group_by(RecYear,yrsins_fire,Unit,Watershed,FireGrzTrt,Transect,sp)%>%
  summarise(abundance=mean(abundance, na.rm=T))%>%
  mutate(abundance=abundance/100)%>%
  pivot_wider(names_from = sp, values_from = abundance, values_fill = 0)
#subsetting environmental and sp data
burn_time_sp_data_north<-burn_time_sp_comp_north%>%
  ungroup()%>%
  select(-1:-6)
burn_time_env_data_north<-burn_time_sp_comp_north%>%
  select(1:6)
dist_north<-vegdist(burn_time_sp_data_north)
permanova_north<-adonis(dist_north~burn_time_env_data_north$yrsins_fire+burn_time_env_data_north$Watershed+as.factor(burn_time_env_data_north$RecYear))
permanova_north$aov.tab

#pairwise comparisons
pairwise.adonis2(burn_time_sp_data_north~yrsins_fire+Watershed+as.factor(RecYear), data=burn_time_env_data_north)

pairwise.adonis2(burn_time_sp_data_north~Watershed+yrsins_fire+as.factor(RecYear), data=burn_time_env_data_north)

pairwise.adonis2(burn_time_sp_data_north~FireGrzTrt+Watershed+as.factor(RecYear), data=burn_time_env_data_north)

#get alpha richness and evenness (average by transect)####
local_sp_comp<-species_comp%>%
  left_join(watershed_key,by="Watershed")%>%
  left_join(Watershed_key2,by="Watershed")%>%
  mutate(year_watershed=paste(RecYear,Watershed,sep="_"))%>%
  left_join(YrSinceFire_key, by="year_watershed")%>%
  group_by(RecYear,Unit,Watershed,FireGrzTrt,Transect,sp)%>%
  summarise(abundance=mean(abundance, na.rm=T))%>%
  mutate(abundance=abundance/100)%>%
  mutate(Rep_id=paste(Unit,FireGrzTrt,Watershed,Transect, sep="_"))

#deriving richness and evenness using codyn at landscape scale
local_rich <- community_structure(local_sp_comp, time.var = "RecYear", 
                               abundance.var = "abundance",
                               replicate.var = "Rep_id", metric = "Evar")



#join datasets
local_rich_ready<- local_sp_comp %>%
  ungroup()%>%
  dplyr::select(RecYear,Unit,Watershed,FireGrzTrt,Transect,Rep_id) %>%
  distinct()%>%#remove repeated rows
  left_join(local_rich, by=c("Rep_id","RecYear"))

#convert to factors
local_rich_ready$RecYear<-as.factor(local_rich_ready$RecYear)

#mixed anova of local richness and evenness
local_rich_model<-lmer(richness~FireGrzTrt*RecYear+(1|Unit/Watershed), 
                       data=local_rich_ready)
anova(local_rich_model)
check_model(local_rich_model)#good enough
qqnorm(resid(local_rich_model))
summary(local_rich_model)

#using mean estimate to create figure 
local_rich_interact<-interactionMeans(local_rich_model)
#replacing spaces in column names with underscore 
names(local_rich_interact)<-str_replace_all(names(local_rich_interact), " ","_")
#df for visuals from model estimates
local_rich_interact_viz<-local_rich_interact%>%
  mutate(local_rich_mean=adjusted_mean,
         local_rich_upper=adjusted_mean+SE_of_link,
         local_rich_lower=adjusted_mean-SE_of_link)
#visual
ggplot(local_rich_interact_viz,aes(RecYear, local_rich_mean, col=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=local_rich_lower,
                    ymax=local_rich_upper),width=0.1,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

#average across years for simplification
local_rich_interact_bar<-local_rich_interact_viz%>%
  group_by(FireGrzTrt)%>%
  summarise(local_mean=mean(local_rich_mean),
            se_upper=mean(local_rich_upper),
            se_lower=mean(local_rich_lower))
ggplot(local_rich_interact_bar,aes(x=FireGrzTrt,fill=FireGrzTrt))+
  geom_bar(stat = "identity",aes(y=local_mean),width = 0.25)+
  geom_errorbar(aes(ymin=se_lower,
                    ymax=se_upper),width=0.05,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"))

local_evar_model<-lmer(Evar~FireGrzTrt*RecYear+(1|Unit/Watershed),
                       data=local_rich_ready)
anova(local_evar_model)
check_model(local_evar_model)
qqnorm(resid(local_evar_model))#good enough

#using mean estimate to create figure 
local_ever_interact<-interactionMeans(local_evar_model)
#replacing spaces in column names with underscore 
names(local_ever_interact)<-str_replace_all(names(local_ever_interact), " ","_")
#df for visuals from model estimates
local_ever_interact_viz<-local_ever_interact%>%
  mutate(local_evar_mean=adjusted_mean,
         local_evar_upper=adjusted_mean+SE_of_link,
         local_evar_lower=adjusted_mean-SE_of_link)
#visual
ggplot(local_ever_interact_viz,aes(RecYear, local_evar_mean, col=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=local_evar_lower,
                    ymax=local_evar_upper),width=0.1,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

#average across years for simplification
local_evar_interact_bar<-local_ever_interact_viz%>%
  group_by(FireGrzTrt)%>%
  summarise(local_mean=mean(local_evar_mean),
            se_upper=mean(local_evar_upper),
            se_lower=mean(local_evar_lower))
ggplot(local_evar_interact_bar,aes(x=FireGrzTrt,fill=FireGrzTrt))+
  geom_bar(stat = "identity",aes(y=local_mean),width = 0.25)+
  geom_errorbar(aes(ymin=se_lower,
                    ymax=se_upper),width=0.05,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"))

#create abundance based on lifeforms####
species_comp_abund<-species_comp_data%>%
  left_join(cover_key, by="CoverClass")%>%
  #create unique name for species by combining binomial nomenclature
  mutate(sp=paste(Ab_genus,Ab_species, sep="_"))%>%
  group_by(Watershed, RecYear,Transect,Plot,sp,SpeCode)%>%
  #selecting the maximum cover for each species
  summarise(abundance=max(abundance, na.rm=T))%>%
  #removing unwanted years
  filter(!RecYear%in%2008:2010)%>%
  filter(RecYear!=2022)

#import data as dataframe
lifeforms_data<- read_excel("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/PhD Wyoming_One Drive/PHD Wyoming/Thesis/PBG synthesis/sp_list_update.xlsx")%>%
  mutate(life_form=paste(growthform,lifeform, sep="_"))%>%
  rename(SpeCode=code)

#merging key with dataset
sp_comp_abund<-species_comp_abund%>%
  left_join(watershed_key,by="Watershed")%>%
  left_join(Watershed_key2,by="Watershed")%>%
  left_join(lifeforms_data, by="SpeCode")%>%
  filter(Transect!="C" | Watershed!="C03C")%>%
  filter(Watershed!="C03C" | Transect!="D")%>%
  filter(Watershed!="C03B" | Transect!="A")%>%
  filter(Watershed!="C03B" | Transect!="B")%>%
  filter(Watershed!="C03B" | Transect!="C")%>%
  filter(Watershed!="C03A" | Transect!="A")%>%
  filter(Watershed!="C03A" | Transect!="B")%>%
  filter(Watershed!="C03A" | Transect!="C")%>%
  filter(Watershed!="C3SC" | Transect!="A")%>%
  filter(Watershed!="C3SC" | Transect!="B")%>%
  filter(Watershed!="C3SC" | Transect!="C")%>%
  filter(Watershed!="C3SA" | Transect!="B")%>%
  filter(Watershed!="C3SA" | Transect!="C")%>%
  filter(Watershed!="C3SA" | Transect!="D")%>%
  filter(Watershed!="C3SB" | Transect!="A")%>%
  filter(Watershed!="C3SB" | Transect!="B")%>%
  group_by(Unit,RecYear,FireGrzTrt,sp)%>%
  mutate(abundance_avg=mean(abundance,na.rm=T))%>%
  group_by(Unit,RecYear,FireGrzTrt)%>%
  mutate(abundance_mean=(abundance_avg/sum(abundance_avg, na.rm=T)))%>%
  group_by(Unit,RecYear,FireGrzTrt,life_form)%>%
  summarise(abundance_m=sum(abundance_mean,na.rm=T))%>%
  group_by(Unit,FireGrzTrt,life_form)%>%
  summarise(abundance=mean(abundance_m,na.rm=T))

#prepare data for figure 
sp_comp_abund_viz<-sp_comp_abund%>%
  group_by(FireGrzTrt,life_form)%>%
  summarise(abundance=mean(abundance,na.rm=T))
ggplot(sp_comp_abund_viz, aes(y=abundance,x=life_form,fill=FireGrzTrt))+
  geom_bar(position = "fill", stat = "identity")+#position fill converts ABG,PBG to percenatge for each lifeform
  scale_fill_manual(values=c("#F0E442", "#009E73"))

#alternative figure
ggplot(sp_comp_abund_viz, aes(y=abundance,fill=life_form,x=FireGrzTrt))+
  geom_bar(position = "stack", stat = "identity")+
  scale_fill_bluebrown()
