####PBG SYNTHESIS PROJECT
####Plant Biomass from diskpasture meter
###Author: Joshua Ajowele
###collaborators: PBG synthesis group
###last update:1/27/2024

#packages 
library(tidyverse)
library(ggthemes)
library(readr)
library(performance)
library(car)
library(lme4)
library(nlme)
library(see)
library(patchwork)
library(phia)
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

#averaging at plot scale
biomass_plot_scale<-biomass_data%>%
  group_by(RecYear, Unit, Watershed, FireGrzTrt, Transect, Plotnum)%>%
  summarise(biomass_plot= mean(biomass, na.rm=T),
            biomass_sd=sd(biomass))%>%
  ungroup()
  

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

#####creating dataframe for all scale
#sd and cv among plots
biomass_sd_cv_plot<-biomass_plot_scale%>%
  group_by(RecYear, Unit, Watershed, FireGrzTrt, Transect)%>%
  summarise(biomass_sd_plot=sd(biomass_plot),
            biomass_cv_plot=sd(biomass_plot)/mean(biomass_plot, na.rm=T))%>%
  ungroup()
#biomass plot temporal stability
biomass_temp_stab_plot<-biomass_plot_scale%>%
  group_by(Unit, Watershed, FireGrzTrt, Transect, Plotnum)%>%
  summarise(biomass_stab_plot=mean(biomass_plot, na.rm=T)/sd(biomass_plot))%>%
  ungroup()

##transect scale
#biomass at transect scale
biomass_mean_transect_scale<-biomass_plot_scale%>%
  group_by(RecYear, Unit, Watershed, FireGrzTrt, Transect)%>%
  summarise(biomass_transect=mean(biomass_plot, na.rm=T))%>%
  ungroup()
#sd and cv among transect
biomass_sd_cv_transect<-biomass_mean_transect_scale%>%
  group_by(RecYear, Unit, Watershed, FireGrzTrt)%>%
  summarise(biomass_sd_transect=sd(biomass_transect),
            biomass_cv_transect=sd(biomass_transect)/mean(biomass_transect, na.rm=T))%>%
  ungroup()
#biomass transect temporal stability
biomass_temp_stab_transect<-biomass_mean_transect_scale%>%
  group_by(Unit, Watershed, FireGrzTrt, Transect)%>%
  summarise(biomass_stab_transect=mean(biomass_transect, na.rm=T)/sd(biomass_transect))%>%
  ungroup()

###watershed scale-probably irrelevant to the questions
#biomass mean at watershed scale
biomass_watershed_scale<-biomass_mean_transect_scale%>%
  group_by(RecYear, Unit, Watershed, FireGrzTrt)%>%
  summarise(biomass_watershed=mean(biomass_transect, na.rm=T))%>%
  ungroup()
#sd and cv among watershed
biomass_sd_cv_watershed<-biomass_watershed_scale%>%
  group_by(RecYear, Unit, FireGrzTrt)%>%
  summarise(biomass_sd_watershed=sd(biomass_watershed),
            biomass_cv_watershed=sd(biomass_watershed)/mean(biomass_watershed, na.rm=T))%>%
  ungroup()#NA for ABG since it has only one watershed per unit
#biomass watershed temporal stability
biomass_temp_stab_watershed<-biomass_watershed_scale%>%
  group_by(Unit, Watershed, FireGrzTrt)%>%
  summarise(biomass_stab_watershed=mean(biomass_watershed, na.rm=T)/sd(biomass_watershed))%>%
  ungroup()

         
###Linear mixed models###
#convert all to factors
biomass_plot_scale$RecYear<-as.factor(biomass_plot_scale$RecYear)
biomass_plot_scale$Unit<-as.factor(biomass_plot_scale$Unit)
biomass_plot_scale$Watershed<-as.factor(biomass_plot_scale$Watershed)
biomass_plot_scale$Transect<-as.factor(biomass_plot_scale$Transect)
biomass_plot_scale$Plotnum<-as.factor(biomass_plot_scale$Plotnum)
biomass_plot_scale$FireGrzTrt<-as.factor(biomass_plot_scale$FireGrzTrt)

#run a linear mixed model at the plot scale
Biom_Plot_Model<-lmer(log(biomass_plot)~FireGrzTrt*RecYear+(1|Unit/Watershed/Transect),
                      data=biomass_plot_scale)
summary(Biom_Plot_Model)
check_model(Biom_Plot_Model)#not the best
Anova(Biom_Plot_Model, type=3)#significant interaction
#confirm normality plot is similar
qqnorm(residuals(Biom_Plot_Model))

#using phia for pairwise comparison
#gives mean for each treatment and year
interactionMeans(Biom_Plot_Model)
#pairwise comparison with Benjamin Hochberg adjustment
testInteractions(Biom_Plot_Model, fixed="RecYear",
                 across = "FireGrzTrt", adjustment="BH")#2018 significant
##cross check model with nlme
Biom_Plot_Model_lme<- lme(biomass_plot~FireGrzTrt*RecYear,random= ~1|Unit/Watershed/Transect,
                      data =biomass_plot_scale, 
                      correlation=corCompSymm(form = ~1|Unit/Watershed/Transect)
                      , control=lmeControl(returnObject=TRUE)
                      , na.action = na.omit)
Anova(Biom_Plot_Model_lme, type=3)#the same output as lmer
check_model(Biom_Plot_Model_lme)
#pairwise comparison
testInteractions(Biom_Plot_Model_lme, fixed="RecYear",
                 across = "FireGrzTrt", adjustment="BH")#same output as lmer

#model for sd
#convert all to factors
biomass_sd_cv_plot$RecYear<-as.factor(biomass_sd_cv_plot$RecYear)
biomass_sd_cv_plot$Unit<-as.factor(biomass_sd_cv_plot$Unit)
biomass_sd_cv_plot$Watershed<-as.factor(biomass_sd_cv_plot$Watershed)
biomass_sd_cv_plot$Transect<-as.factor(biomass_sd_cv_plot$Transect)
biomass_sd_cv_plot$FireGrzTrt<-as.factor(biomass_sd_cv_plot$FireGrzTrt)

Biom_Plot_Sd_Model<-lmer(log(biomass_sd_plot)~FireGrzTrt*RecYear+(1|Unit/Watershed),
                      data=biomass_sd_cv_plot)
check_model(Biom_Plot_Sd_Model)
Anova(Biom_Plot_Sd_Model, type=3)
qqnorm(residuals(Biom_Plot_Sd_Model))
#pairwise comparison
testInteractions(Biom_Plot_Sd_Model, fixed="RecYear",
                 across = "FireGrzTrt", adjustment="BH")

#model for spatial cv
Biom_Plot_Cv_Model<-lmer(log(biomass_cv_plot)~FireGrzTrt*RecYear+(1|Unit/Watershed),
                         data=biomass_sd_cv_plot)
check_model(Biom_Plot_Cv_Model)
Anova(Biom_Plot_Cv_Model, type=3)
qqnorm(residuals(Biom_Plot_Cv_Model))
testInteractions(Biom_Plot_Cv_Model, fixed="RecYear",
                 across = "FireGrzTrt", adjustment="BH")

#model for temp stability 
#convert all to factors
biomass_temp_stab_plot$Unit<-as.factor(biomass_temp_stab_plot$Unit)
biomass_temp_stab_plot$Watershed<-as.factor(biomass_temp_stab_plot$Watershed)
biomass_temp_stab_plot$Transect<-as.factor(biomass_temp_stab_plot$Transect)
biomass_temp_stab_plot$Plotnum<-as.factor(biomass_temp_stab_plot$Plotnum)
biomass_temp_stab_plot$FireGrzTrt<-as.factor(biomass_temp_stab_plot$FireGrzTrt)

Biom_Stab_Plot_Model<-lmer(sqrt(biomass_stab_plot)~FireGrzTrt+(1|Unit/Watershed/Transect),
                                data=biomass_temp_stab_plot)
check_model(Biom_Stab_Plot_Model)
Anova(Biom_Stab_Plot_Model, type=3)#no difference
qqnorm(residuals(Biom_Stab_Plot_Model))

#run a linear mixed model at the transect scale
#convert all to factors
biomass_mean_transect_scale$Unit<-as.factor(biomass_mean_transect_scale$Unit)
biomass_mean_transect_scale$Watershed<-as.factor(biomass_mean_transect_scale$Watershed)
biomass_mean_transect_scale$Transect<-as.factor(biomass_mean_transect_scale$Transect)
biomass_mean_transect_scale$RecYear<-as.factor(biomass_mean_transect_scale$RecYear)
biomass_mean_transect_scale$FireGrzTrt<-as.factor(biomass_mean_transect_scale$FireGrzTrt)

Biom_Transect_Model<-lmer(log(biomass_transect)~FireGrzTrt*RecYear+(1|Unit/Watershed),
                      data=biomass_mean_transect_scale)
check_model(Biom_Transect_Model)
Anova(Biom_Transect_Model, type=3)#no significant interaction
qqnorm(residuals(Biom_Transect_Model))


#sd transect
#convert all to factors
biomass_sd_cv_transect$Unit<-as.factor(biomass_sd_cv_transect$Unit)
biomass_sd_cv_transect$Watershed<-as.factor(biomass_sd_cv_transect$Watershed)
biomass_sd_cv_transect$RecYear<-as.factor(biomass_sd_cv_transect$RecYear)
biomass_sd_cv_transect$FireGrzTrt<-as.factor(biomass_sd_cv_transect$FireGrzTrt)

Biom_Sd_Transect_Model<-lmer(log(biomass_sd_transect)~FireGrzTrt*RecYear+(1|Unit),
                          data=biomass_sd_cv_transect)
check_model(Biom_Sd_Transect_Model)
Anova(Biom_Sd_Transect_Model, type=3)#no significant interaction
qqnorm(residuals(Biom_Sd_Transect_Model))#residual normality not the best

#spatial CV
Biom_Cv_Transect_Model<-lmer(log(biomass_cv_transect)~FireGrzTrt*RecYear+(1|Unit),
                             data=biomass_sd_cv_transect)#isSingular warning
check_model(Biom_Cv_Transect_Model)#normality not the best
Anova(Biom_Cv_Transect_Model, type=3)#no significant interaction



#Transect Biom stab
#convert all to factors
biomass_temp_stab_transect$Unit<-as.factor(biomass_temp_stab_transect$Unit)
biomass_temp_stab_transect$Watershed<-as.factor(biomass_temp_stab_transect$Watershed)
biomass_temp_stab_transect$FireGrzTrt<-as.factor(biomass_temp_stab_transect$FireGrzTrt)

Biom_Stab_Transect_model<-lmer(log(biomass_stab_transect)~FireGrzTrt+(1|Unit/Watershed),
                               data=biomass_temp_stab_transect)
check_model(Biom_Stab_Transect_model)#iffy normality
Anova(Biom_Stab_Transect_model, type=3)#no difference
qqnorm(residuals(Biom_Stab_Transect_model))



###wrangle data for visualization


##visuals
#biomass plot
ggplot(biomass_master_viz,aes(RecYear, biomass_plot_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_plot_mean-biomass_plot_se,
                    ymax=biomass_plot_mean+biomass_plot_se),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))
#saved width 800, height=400
#plot sd
ggplot(biomass_master_viz,aes(RecYear, biomass_plot_sd_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_plot_sd_mean-biomass_plot_sd_se,
                    ymax=biomass_plot_sd_mean+biomass_plot_sd_se),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))
#plot spatial cv
ggplot(biomass_master_viz,aes(RecYear, biomass_plot_cv_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_plot_cv_mean-biomass_plot_cv_se,
                    ymax=biomass_plot_cv_mean+biomass_plot_cv_se),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))


#biomass transect
ggplot(biomass_master_viz,aes(RecYear, biomass_transect_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_transect_mean-biomass_transect_se,
                    ymax=biomass_transect_mean+biomass_transect_se),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))
#saved width 800, height=400
#transect sd
ggplot(biomass_master_viz,aes(RecYear, biomass_transect_sd_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_transect_sd_mean-biomass_transect_sd_se,
                    ymax=biomass_transect_sd_mean+biomass_transect_sd_se),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))
#transect spatial cv
ggplot(biomass_master_viz,aes(RecYear, biomass_transect_cv_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_transect_cv_mean-biomass_transect_cv_se,
                    ymax=biomass_transect_cv_mean+biomass_transect_cv_se),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

#biomass watershed
ggplot(biomass_master_viz,aes(RecYear, biomass_watershed_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_watershed_mean-biomass_watershed_se,
                    ymax=biomass_watershed_mean+biomass_watershed_se),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))
#saved width 800, height=400
#watershed sd
ggplot(biomass_master_viz,aes(RecYear, biomass_watershed_sd_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_watershed_sd_mean-biomass_watershed_sd_se,
                    ymax=biomass_watershed_sd_mean+biomass_watershed_sd_se),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))
#transect spatial cv
ggplot(biomass_master_viz,aes(RecYear, biomass_watershed_cv_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_watershed_cv_mean-biomass_watershed_cv_se,
                    ymax=biomass_watershed_cv_mean+biomass_watershed_cv_se),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

###create a separate dataframe for stability across all scale
biomass_stab_master<-biomass_master%>%
  select(1:6,biom_plot_stab, biom_transect_stab, biom_watershed_stab)%>%
  pivot_longer(7:9, names_to = "Scale", values_to = "biomass_stability")%>%
  group_by(FireGrzTrt, Scale)%>%
  summarise(Biomass_stability=mean(biomass_stability, na.rm=T),
            biom_stab_se=SE_function(biomass_stability))
  
scale_key<-tibble(Scale=levels(factor(biomass_stab_master$Scale)),
                 Key=c("plot","transect","watershed"))%>%
  left_join(biomass_stab_master, by="Scale")

ggplot(scale_key,aes(Key, Biomass_stability, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(Key)),linewidth=1)+
  geom_errorbar(aes(ymin=Biomass_stability-biom_stab_se,
                    ymax=Biomass_stability+biom_stab_se),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))
