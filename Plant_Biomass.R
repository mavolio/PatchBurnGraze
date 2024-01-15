####PBG SYNTHESIS PROJECT
####Plant Biomass from diskpasture meter
###Author: Joshua Ajowele

#packages 
library(tidyverse)
library(ggthemes)
library(readr)
library(performance)
library(car)
### Standard Error function
SE_function<-function(x,na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x))
  return(SE)
}

### Set graphing parameters
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=20), legend.text=element_text(size=20))



#import dataset as tibble
biomass_diskpasture<- read_csv("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/PhD Wyoming_One Drive/PHD Wyoming/Thesis/PBG synthesis/Plant_Biomass_11_17_PBG031.csv")
#modifying the data set by adding new columns
biomass_reg <- biomass_diskpasture%>%
  #creating row number in order to visualise outliers in the graph
  mutate(row_num= 1:length(Lvgrass))%>%
  #the 2011 data were not previously converted to the appropriate scale
  filter(!Recyear==2011)%>%
  mutate(Lvgrass_clean = Lvgrass*10,
         Pdead_clean= Pdead*10,
         Forbs_clean = Forbs*10,
         Woody_clean= Woody*10)%>%
  mutate(total_biomass= Lvgrass_clean+Pdead_clean+Forbs_clean+Woody_clean)%>%
  mutate(total_nonwoody=Lvgrass_clean+Pdead_clean+Forbs_clean)%>%
  mutate(forbs_grass=Lvgrass_clean+Forbs_clean)

##performing a regression analysis
#using total biomass
#Square root transformation 
total_reg<- lm(sqrt(total_biomass)~Diskht, data = biomass_reg)
summary(total_reg)
check_model(total_reg)

#visual:show rown mumbers to identify outlier
ggplot(biomass_reg, aes(Diskht, total_biomass, label = row_num))+
  geom_text()+
  geom_smooth(method = "lm")+
  theme_few()

#visual
ggplot(biomass_reg, aes(Diskht, total_biomass))+
  geom_point(col="#E69F00")+
  geom_smooth(method = "lm", col="#000000")+
  theme_few()+
  xlab("Height (cm)")+
  ylab("Biomass (g/m2)")

#checking model assumptions with figures
check_model(total_reg)
#manually checking the assumptions
par(mfrow=c(2,2))
plot(total_reg, 1:3) 
hist(sqrt(biomass_reg$total_biomass))

##import dataset as tibble
#estimate biomass from regression equation
Diskpast_data<- read_csv("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/PhD Wyoming_One Drive/PHD Wyoming/Thesis/PBG synthesis/Plant_Biomass_10_20_PBG032.csv")%>%
  mutate(biom= 9.99+(0.7965*(Diskht)))%>%
  mutate(biomass= biom*biom)%>%
  rename(RecYear= Recyear)%>%
  #creating a coloumn for unique combination of year and watershed
  mutate(year_watershed = paste(RecYear, Watershed, sep = "_"))

#Create a watershed key for treatment and Unit column to merge with the raw data
watershed_key <- tibble(Watershed=levels(factor(Diskpast_data$Watershed)),
                        FireGrzTrt=c("ABG", "PBG", "PBG", "PBG", "ABG", "PBG", "PBG", "PBG"),
                        Unit=c("south", "south", "south", "south", "north", "north",
                               "north", "north"))

#check to make sure there is no problem
print(watershed_key)

#join watershed key with data
biomass_data <- Diskpast_data %>%
  filter(!RecYear=="2010")%>%
  filter(!RecYear=="2011")%>%
  left_join(watershed_key, by="Watershed")

#averaging to plot level for plot scale comparison
biomass_plot_scale<-biomass_data%>%
  group_by(RecYear, Unit, Watershed, FireGrzTrt, Transect, Plotnum)%>%
  summarise(biomass_plot= mean(biomass, na.rm=T))

###personal use
#meeting with Kevin 12/10/2023
#at each scale get the raw values for all possible combination and 
#this can be used to calculate plot level stability or average to calculate higher scale 
#stability

#filter for PBG to perform bootstrap
biomass_PBG_plot_index<-biomass_plot_scale%>%
  filter(FireGrzTrt=="PBG")%>%
  group_by(RecYear)%>%
  #number the plots to sample from
  mutate(biomass_PBG_plot_index= 1:length(RecYear))

biomass_rand_plot<-biomass_plot_scale%>%
  filter(FireGrzTrt=="PBG")%>%
  left_join(biomass_PBG_plot_index, by=c("RecYear","Unit","Watershed","Transect", "FireGrzTrt",
                                         "Plotnum", "biomass_plot"))
  
num_bootstrap<-1000
bootstrap_vector<-1:num_bootstrap
biomass_plot_master<-{}
for(BOOT in 1:length(bootstrap_vector)){
  biomass_plot_rand_key<-biomass_rand_plot%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Transect, Plotnum, biomass_PBG_plot_index)%>%
    unique(.)%>%
    #filter(RecYear==2012)%>%
    group_by(RecYear)%>%
    sample_n(400, replace=T)%>%
    dplyr::select(biomass_PBG_plot_index)%>%
    ungroup()
  biomass_plot_ready<-biomass_rand_plot%>%
    right_join(biomass_plot_rand_key, by= c("RecYear", "biomass_PBG_plot_index"),
              multiple="all")%>%
    mutate(iteration=BOOT)
  biomass_plot_master<-rbind(biomass_plot_master, biomass_plot_ready)
}

#write into csv
write.csv(biomass_plot_master,"C:/Users/joshu/OneDrive - UNCG/UNCG PHD/PhD Wyoming_One Drive/PHD Wyoming/Thesis/PBG synthesis/biomass_plot_master.csv" )

biomass_plot_bootdata<-read_csv("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/PhD Wyoming_One Drive/PHD Wyoming/Thesis/PBG synthesis/biomass_plot_master.csv")%>%
  dplyr::select(-1)

#mean and SD of PBG biomass at plot scale per iteration
PBG_biomass_plot<-biomass_plot_bootdata%>%
  group_by(RecYear, iteration)%>%
  mutate(PBG_plot_mean_dist=mean(biomass_plot, na.rm=T),
         PBG_plot_SD_dist=sd(biomass_plot))%>%
  select(RecYear,PBG_plot_mean_dist, PBG_plot_SD_dist,iteration)%>%
  unique()%>%
  #obtain the mean and SD of the mean and SD distributions
  group_by(RecYear)%>%
  mutate(PBG_plot_mean=mean(PBG_plot_mean_dist, na.rm=T),
         PBG_plot_mean_sd=sd(PBG_plot_mean_dist),
         PBG_plot_SD_mean= mean(PBG_plot_SD_dist),
         PBG_plot_SD_SD=sd(PBG_plot_SD_dist))

#calculate a z score from PBG distribution with the mean of ABG
#filter ABG from the dataframe
#calculate ABG mean at the plot scale
ABG_biomass_plot<-biomass_plot_scale%>%
  filter(FireGrzTrt=="ABG")%>%
  group_by(RecYear)%>%
  summarise(ABG_plot_mean=mean(biomass_plot, na.rm=T),
         ABG_plot_SD=sd(biomass_plot))


#combine ABG and PBG plot biomass 
combo_plot_biomass<-ABG_biomass_plot%>%
  left_join(PBG_biomass_plot, by ="RecYear")%>%
  #calculate the z-scores
  mutate(Z_score_mean=(ABG_plot_mean-PBG_plot_mean)/PBG_plot_mean_sd,
         Z_score_SD=(ABG_plot_SD-PBG_plot_SD_mean)/PBG_plot_SD_SD)%>%
  #calculate p values from z-scores
  mutate(pvalue_mean=2*pnorm(-abs(Z_score_mean)),
         pvalue_sd=2*pnorm(-abs(Z_score_SD)))

#create a visual of the distribution
ggplot(combo_plot_biomass,aes(PBG_plot_mean_dist))+
  geom_density()+
  facet_wrap(~RecYear)+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=ABG_plot_mean, col="Red"))


#visual for SD
ggplot(combo_plot_biomass,aes(PBG_plot_SD_dist))+
  geom_density()+
  facet_wrap(~RecYear)+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=ABG_plot_SD, col="Red"))


####compare PBG mean and SD from actual observation to bootstrapped values
#cal mean and sd from PBG actual value
PBG_actual_biomass<-biomass_plot_scale%>%
  filter(FireGrzTrt=="PBG")%>%
  group_by(RecYear)%>%
  summarise(a_biom_mean= mean(biomass_plot, na.rm=T),
            a_biom_sd=sd(biomass_plot))
#combine
PBG_boot_biom<-combo_plot_biomass%>%
  select(RecYear, PBG_plot_mean, PBG_plot_SD_mean)%>%
  unique()%>%
  left_join(PBG_actual_biomass, by="RecYear")

#perform a linear regression
mean_PBG_lm<-lm(PBG_plot_mean~a_biom_mean, data=PBG_boot_biom)
summary(mean_PBG_lm)
Sd_PBG_lm<-lm(PBG_plot_SD_mean~a_biom_sd, data=PBG_boot_biom)
summary(Sd_PBG_lm)

#visual
ggplot(PBG_boot_biom, aes(PBG_plot_mean, a_biom_mean))+
  geom_point()+
  geom_smooth(method="lm")

ggplot(PBG_boot_biom, aes(PBG_plot_SD_mean, a_biom_sd))+
  geom_point()+
  geom_smooth(method="lm")
#conclusion: There is no point bootstrapping!!!