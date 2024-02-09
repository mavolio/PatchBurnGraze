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
grassh_count_tnst_model<-lmer(log(Tcount)~FireGrzTrt+RecYear+(1|Unit/Watershed),
                              data=grassh_count_df)
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
