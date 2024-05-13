#Patch-Burn Synthesis Project
#Plant community data at the landscape scale
#Author: Joshua Adedayo Ajowele joshuaajowele@gmail.com
#Started: May 13, 2024

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
  filter(!RecYear%in%2008:2010)%>%
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
pl_diver_model<-lmer(log(Shannon)~FireGrzTrt*RecYear+(1|Unit),
                     data=pl_rich_diver)
Anova(pl_diver_model, type=3)
qqnorm(resid(pl_diver_model))
check_model(pl_diver_model)

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
ggplot(plrich_interact_viz,aes(RecYear, plrich_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=plrich_lower,
                    ymax=plrich_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

#average across years for simplification
plrich_interact_bar<-plrich_interact_viz%>%
  group_by(FireGrzTrt)%>%
  summarise(pl_richness=mean(plrich_mean),
            se_upper=mean(plrich_upper),
            se_lower=mean(plrich_lower))
ggplot(plrich_interact_bar,aes(x=FireGrzTrt,fill=FireGrzTrt))+
  geom_bar(stat = "identity",aes(y=pl_richness),width = 0.5)+
  geom_errorbar(aes(ymin=se_lower,
                    ymax=se_upper),width=0.2,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"))

#evenness model
plevenness_model<-lmer(Evar~FireGrzTrt*RecYear+(1|Unit),
                       data=pl_rich_diver)#is.singular
Anova(plevenness_model, type=3)
qqnorm(resid(plevenness_model))
check_model(plevenness_model)
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
ggplot(pleven_interact_viz,aes(RecYear, plevenness_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=pleven_lower,
                    ymax=pleven_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

#average across years for simplification
pleven_interact_bar<-pleven_interact_viz%>%
  group_by(FireGrzTrt)%>%
  summarise(pl_evenness=mean(plevenness_mean),
            se_upper=mean(pleven_upper),
            se_lower=mean(pleven_lower))
ggplot(pleven_interact_bar,aes(x=FireGrzTrt,fill=FireGrzTrt))+
  geom_bar(stat = "identity",aes(y=pl_evenness),width = 0.5)+
  geom_errorbar(aes(ymin=se_lower,
                    ymax=se_upper),width=0.2,linetype=1)+
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
#Run 20 stress 0.276

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
pl_beta$RecYear<-as.factor(pl_beta$RecYear)
pl_beta_model<-lmer(log(distance)~FireGrzTrt*RecYear+(1|Unit),
                    data=pl_beta)
Anova(pl_beta_model, type=3)
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
ggplot(model_estimates_beta_viz,aes(RecYear, pl_betadiv, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=pl_lower,
                    ymax=pl_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))
#summarize with a bargraph
plbeta_interact_bar<-model_estimates_beta_viz%>%
  group_by(FireGrzTrt)%>%
  summarise(pl_betadiver=mean(pl_betadiv),
            se_upper=mean(pl_upper),
            se_lower=mean(pl_lower))
ggplot(plbeta_interact_bar,aes(x=FireGrzTrt,fill=FireGrzTrt))+
  geom_bar(stat = "identity",aes(y=pl_betadiver),width = 0.5)+
  geom_errorbar(aes(ymin=se_lower,
                    ymax=se_upper),width=0.2,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"))

#simper analysis to determine species contributing 80% of difference in compositioin
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
ggplot(pl_cover_2021,aes(x=sp))+
  geom_bar(stat = "identity",aes(y=abund_ABG_PBG),width = 0.5)+
  scale_x_discrete(guide = guide_axis(angle = 90))



#just checking if any species is an indicator distinguishing treatment
library(indicspecies)
indicator_abgvspbg<-multipatt(pl_sp_data_2021,pl_env_data_2021$FireGrzTrt)
ff<-summary(indicator_abgvspbg)
#There are!

