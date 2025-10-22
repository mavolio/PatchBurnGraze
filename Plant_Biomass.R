####PBG SYNTHESIS PROJECT
####Plant Biomass from diskpasture meter
###Author: Joshua Ajowele
###collaborators: PBG synthesis group
###last update:Oct 22 2025

#packages 
library(tidyverse)
library(ggthemes)
library(readr)
library(performance)
library(car)
library(lme4)
library(lmerTest)
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
biomass_diskpasture<- read_csv("C:/Users/Joshua/OneDrive - UNCG/UNCG PHD/PhD Wyoming_One Drive/PHD Wyoming/Thesis/PBG synthesis/Plant_Biomass_11_17_PBG031.csv")
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
Diskpast_data<- read_csv("C:/Users/Joshua/OneDrive - UNCG/UNCG PHD/PhD Wyoming_One Drive/PHD Wyoming/Thesis/PBG synthesis/Plant_Biomass_10_20_PBG032.csv")%>%
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

#averaging at plot scale####
biomass_plot_scale<-biomass_data%>%
  group_by(RecYear, Unit, Watershed, FireGrzTrt, Transect, Plotnum)%>%
  summarise(biomass_plot= mean(biomass, na.rm=T),
            biomass_sd=sd(biomass))%>%
  ungroup()
  

#filter for PBG to perform bootstrap at plot scale####
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

#mean and SD of PBG biomass at plot scale per iteration####
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

#create a visual of the distribution at plot scale####
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


####compare PBG mean and SD from actual observation to bootstrapped values####
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
#conclusion:  bootstrapping adequately represents observed value!

#####creating dataframe for all scale####
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

         
###Linear mixed models####
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

#model for temp stability plot scale ####
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
interactionMeans(Biom_Stab_Plot_Model)
#run a linear mixed model at the transect scale####
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

###wrangle data for visualization####
biomass_plot_scale_viz<-biomass_plot_scale%>%
  group_by(RecYear,FireGrzTrt)%>%
  summarise(biomass_plot_mean=mean(biomass_plot, na.rm=T),
            biomass_plot_se=SE_function(biomass_plot),
            cnt=n() )%>%
  ungroup()
#using mean estimate from post-hoc to create figure as a comparison to raw data
model_estimates_plot<-interactionMeans(Biom_Plot_Model)
#replacing spaces in column names with underscore 
names(model_estimates_plot)<-str_replace_all(names(model_estimates_plot), " ","_")
#df for visuals from model estimates
model_estimates_plot_viz<-model_estimates_plot%>%
  mutate(biomass_plot_bt_mean=exp(adjusted_mean),
         biomass_plot_upper=exp(adjusted_mean+SE_of_link),
         biomass_plot_lower=exp(adjusted_mean-SE_of_link))
#model estimates provided a better representation of the results than the 
#raw values

##visuals
#biomass plot
ggplot(biomass_plot_scale_viz,aes(RecYear, biomass_plot_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_plot_mean-biomass_plot_se,
                    ymax=biomass_plot_mean+biomass_plot_se),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))
#saved width 800, height=400
#biomass plot from model estimates

ggplot(model_estimates_plot_viz,aes(RecYear, biomass_plot_bt_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_plot_lower,
                    ymax=biomass_plot_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))
#saved width 800, height=400 pixels for png and 8.3, 4.17 icnhes for Pdf

#plot sd
model_estimate_plot_sd<-interactionMeans(Biom_Plot_Sd_Model)
names(model_estimate_plot_sd)<-str_replace_all(names(model_estimate_plot_sd), " ","_")
#df for visuals from model estimates
model_estimate_plot_sd_viz<-model_estimate_plot_sd%>%
  mutate(biomass_plot_sd_mean=exp(adjusted_mean),
         biomass_plot_upper=exp(adjusted_mean+SE_of_link),
         biomass_plot_lower=exp(adjusted_mean-SE_of_link))
ggplot(model_estimate_plot_sd_viz,aes(RecYear, biomass_plot_sd_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_plot_lower,
                    ymax=biomass_plot_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

#plot spatial cv
model_estimate_plot_cv<-interactionMeans(Biom_Plot_Cv_Model)
names(model_estimate_plot_cv)<-str_replace_all(names(model_estimate_plot_cv), " ","_")
#df for visuals from model estimates
model_estimate_plot_cv_viz<-model_estimate_plot_cv%>%
  mutate(biomass_plot_cv_mean=exp(adjusted_mean),
         biomass_plot_upper=exp(adjusted_mean+SE_of_link),
         biomass_plot_lower=exp(adjusted_mean-SE_of_link))
ggplot(model_estimate_plot_cv_viz,aes(RecYear, biomass_plot_cv_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_plot_lower,
                    ymax=biomass_plot_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))


#biomass transect
model_estimates_transect<-interactionMeans(Biom_Transect_Model)
#replacing spaces in column names with underscore 
names(model_estimates_transect)<-str_replace_all(names(model_estimates_transect), " ","_")
#df for visuals from model estimates
model_estimates_transect_viz<-model_estimates_transect%>%
  mutate(biomass_transect_mean=exp(adjusted_mean),
         biomass_transect_upper=exp(adjusted_mean+SE_of_link),
         biomass_transect_lower=exp(adjusted_mean-SE_of_link))
ggplot(model_estimates_transect_viz,aes(RecYear, biomass_transect_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_transect_lower,
                    ymax=biomass_transect_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))
#saved width 800, height=400
#transect sd
model_estimates_transect_sd<-interactionMeans(Biom_Sd_Transect_Model)
#replacing spaces in column names with underscore 
names(model_estimates_transect_sd)<-str_replace_all(names(model_estimates_transect_sd), " ","_")
#df for visuals from model estimates
model_estimates_transect_sd_viz<-model_estimates_transect_sd%>%
  mutate(biomass_transect_sd_mean=exp(adjusted_mean),
         biomass_transect_upper=exp(adjusted_mean+SE_of_link),
         biomass_transect_lower=exp(adjusted_mean-SE_of_link))
ggplot(model_estimates_transect_sd_viz,aes(RecYear, biomass_transect_sd_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_transect_lower,
                    ymax=biomass_transect_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))
#transect spatial cv
model_estimates_transect_cv<-interactionMeans(Biom_Cv_Transect_Model)
#replacing spaces in column names with underscore 
names(model_estimates_transect_cv)<-str_replace_all(names(model_estimates_transect_cv), " ","_")
#df for visuals from model estimates
model_estimates_transect_cv_viz<-model_estimates_transect_cv%>%
  mutate(biomass_transect_cv_mean=exp(adjusted_mean),
         biomass_transect_upper=exp(adjusted_mean+SE_of_link),
         biomass_transect_lower=exp(adjusted_mean-SE_of_link))
ggplot(model_estimates_transect_cv_viz,aes(RecYear, biomass_transect_cv_mean, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_transect_lower,
                    ymax=biomass_transect_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))



###create a separate dataframe for plot and transect stability ####
model_estimates_plot_stab<-interactionMeans(Biom_Stab_Plot_Model)%>%
  mutate(Scale="Plot")
model_estimates_transect_stab<-interactionMeans(Biom_Stab_Transect_model)%>%
  mutate(Scale="Transect")
#replacing spaces in column names with underscore 
names(model_estimates_plot_stab)<-str_replace_all(names(model_estimates_plot_stab), " ","_")
names(model_estimates_transect_stab)<-str_replace_all(names(model_estimates_transect_stab), " ","_")
model_estimates_transect_stab_update<-model_estimates_transect_stab%>%
  mutate(biomass_stability=exp(adjusted_mean),
         biomass_stab_upper=exp(adjusted_mean+SE_of_link),
         biomass_stab_lower=exp(adjusted_mean-SE_of_link))
#df for visuals from model estimates
model_estimates_stab<-model_estimates_plot_stab%>%
  mutate(biomass_stability=(adjusted_mean)^2,
         biomass_stab_upper=(adjusted_mean+SE_of_link)^2,
         biomass_stab_lower=(adjusted_mean-SE_of_link)^2)%>%
  full_join(model_estimates_transect_stab_update, by=c("FireGrzTrt", "Scale",
                                                "biomass_stability","biomass_stab_upper",
                                                "biomass_stab_lower"))
 
ggplot(model_estimates_stab,aes(Scale, biomass_stability, col=FireGrzTrt, linetype=FireGrzTrt))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(Scale)),linewidth=1)+
  geom_errorbar(aes(ymin=biomass_stab_lower,
                    ymax=biomass_stab_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

###Unit Scale bootstrap separated by unit####
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
PBG_south_biomass_master<-read_csv("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/PBG_south_biomass.csv")%>%
  dplyr::select(-1)

#create an index for the north plots
biomass_north_index<-PBG_north_biomass%>%
  group_by(RecYear)%>%
  mutate(plot_index=1:length(RecYear))

num_bootstrap<-1000
bootstrap_vector<-1:num_bootstrap
PBG_north_biomass_master<-{}
for(BOOT in 1:length(bootstrap_vector)){
  north_rand_key<-biomass_north_index%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Transect, Plotnum, plot_index)%>%
    unique(.)%>%
    group_by(RecYear)%>%
    sample_n(200, replace=T)%>%
    dplyr::select(plot_index,RecYear)%>%
    ungroup()
  biomass_north_ready<-biomass_north_index%>%
    right_join(north_rand_key, by= c("RecYear", "plot_index"),
               multiple="all")%>%
    mutate(iteration=BOOT)
  PBG_north_biomass_master<-rbind(PBG_north_biomass_master, biomass_north_ready)
}
write.csv(PBG_north_biomass_master, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/PBG_north_biomass.csv")
PBG_north_biomass_master<-read_csv("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/PBG_north_biomass.csv")%>%
  dplyr::select(-1)


#extract ABG
#extract ABG and separate into Unit
ABG_south_biomass<-biomass_data%>%
  filter(FireGrzTrt=="ABG" & Unit=="south")%>%
  group_by(RecYear, Unit, Watershed, Transect, Plotnum, FireGrzTrt)%>%
  summarise(biomass=mean(biomass, na.rm=T))%>%
  ungroup()%>%
  group_by(RecYear)%>%#calculate metrics at the unit level
  summarise(biomass_ABGSth=mean(biomass, na.rm=T),
            biomass_ABGSth_sd=sd(biomass),
            cv_biomass_ABGSth=sd(biomass)/mean(biomass, na.rm=T))%>%
  ungroup()%>%
  mutate(Stab_ABGSth=mean(biomass_ABGSth, na.rm =T)/sd(biomass_ABGSth))
ABG_north_biomass<-biomass_data%>%
  filter(FireGrzTrt=="ABG" & Unit=="north")%>%
  group_by(RecYear, Unit, Watershed, Transect, Plotnum, FireGrzTrt)%>%
  summarise(biomass=mean(biomass, na.rm=T))%>%
  ungroup()%>%
  group_by(RecYear)%>%#calculate mean at the unit level
  summarise(biomass_ABGNth=mean(biomass, na.rm=T),
            biomass_ABGNth_sd=sd(biomass),
            cv_biomass_ABGNth=sd(biomass)/mean(biomass, na.rm=T))%>%
  ungroup()%>%
  mutate(Stab_ABGNth=mean(biomass_ABGNth, na.rm =T)/sd(biomass_ABGNth))

#calculate mean for PBG bootstrap
PBG_north_mean<-PBG_north_biomass_master%>%
  group_by(RecYear,iteration)%>%
  summarise(biomass_PBGNth=mean(biomass,na.rm=T),
            sd_biomass_PBGNth=sd(biomass),
            cv_biomass_PBGNth=sd(biomass)/mean(biomass, na.rm=T))%>%
  group_by(iteration)%>%
  mutate(stab_PBGNth=mean(biomass_PBGNth, na.rm=T)/sd(biomass_PBGNth))

PBG_south_mean<-PBG_south_biomass_master%>%
  group_by(RecYear,iteration)%>%
  summarise(biomass_PBGSth=mean(biomass,na.rm=T),
            sd_biomass_PBGSth=sd(biomass),
            cv_biomass_PBGSth=sd(biomass)/mean(biomass, na.rm=T))%>%
  group_by(iteration)%>%
  mutate(stab_PBGSth=mean(biomass_PBGSth, na.rm=T)/sd(biomass_PBGSth))


#combine PBG and ABG
Combo_north_biomass<-PBG_north_mean%>%
  left_join(ABG_north_biomass, by="RecYear")%>%
  group_by(RecYear)%>%
  mutate(biomass_PBGNth_m=mean(biomass_PBGNth, na.rm=T),
         biomass_PBGNth_sd=sd(biomass_PBGNth),
         z_score_NthMean=((biomass_ABGNth-biomass_PBGNth_m)/biomass_PBGNth_sd),
         p_value_NMean=2*pnorm(-abs(z_score_NthMean)),
         biomass_PBGNth_sd_M=mean(sd_biomass_PBGNth, na.rm=T),
         biomass_PBGNth_sd_sd=sd(sd_biomass_PBGNth),
         Z_score_Nsd=((biomass_ABGNth_sd-biomass_PBGNth_sd_M)/biomass_PBGNth_sd_sd),
         pvalue_Nsd=2*pnorm(-abs(Z_score_Nsd)),
         biomass_PBGNth_cv_M=mean(cv_biomass_PBGNth, na.rm=T),
         biomass_PBGNth_cv_sd=sd(cv_biomass_PBGNth),
         Z_score_Ncv=((cv_biomass_ABGNth-biomass_PBGNth_cv_M)/biomass_PBGNth_cv_sd),
         pvalue_Ncv=2*pnorm(-abs(Z_score_Ncv)))%>%
  ungroup()%>%
  mutate(stab_PBGNm=mean(stab_PBGNth,na.rm=T),
         stab_PBGN_sd=sd(stab_PBGNth),
         z_score_NStab=((Stab_ABGNth-stab_PBGNm)/stab_PBGN_sd),
         pvalue_stab=2*pnorm(-abs(z_score_NStab)))
write.csv(Combo_north_biomass, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/north_biomass_boot_result.csv")
 
Combo_south_biomass<-PBG_south_mean%>%
  left_join(ABG_south_biomass, by="RecYear")%>%
  group_by(RecYear)%>%
  mutate(biomass_PBGSth_m=mean(biomass_PBGSth, na.rm=T),
         biomass_PBGSth_sd=sd(biomass_PBGSth),
         z_score_SthMean=((biomass_ABGSth-biomass_PBGSth_m)/biomass_PBGSth_sd),
         p_value_SMean=2*pnorm(-abs(z_score_SthMean)),
         biomass_PBGSth_sd_M=mean(sd_biomass_PBGSth, na.rm=T),
         biomass_PBGSth_sd_sd=sd(sd_biomass_PBGSth),
         Z_score_Ssd=((biomass_ABGSth_sd-biomass_PBGSth_sd_M)/biomass_PBGSth_sd_sd),
         pvalue_Ssd=2*pnorm(-abs(Z_score_Ssd)),
         biomass_PBGSth_cv_M=mean(cv_biomass_PBGSth, na.rm=T),
         biomass_PBGSth_cv_sd=sd(cv_biomass_PBGSth),
         Z_score_Scv=((cv_biomass_ABGSth-biomass_PBGSth_cv_M)/biomass_PBGSth_cv_sd),
         pvalue_Scv=2*pnorm(-abs(Z_score_Scv)))%>%
  ungroup()%>%
  mutate(stab_PBGSm=mean(stab_PBGSth,na.rm=T),
         stab_PBGS_sd=sd(stab_PBGSth),
         z_score_SStab=((Stab_ABGSth-stab_PBGSm)/stab_PBGS_sd),
         pvalue_stab=2*pnorm(-abs(z_score_SStab)))
write.csv(Combo_south_biomass, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/south_biomass_boot_result.csv")

#create a visual of the distribution unit scale by unit####
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

#create visual for stab
ggplot(Combo_north_biomass,aes(stab_PBGNth))+
  geom_density(size=1)+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=Stab_ABGNth), linetype=2,size=1)+
  xlab("Stability North Unit")
ggplot(Combo_south_biomass,aes(stab_PBGSth))+
  geom_density(size=1)+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=Stab_ABGSth), linetype=2,size=1)+
  xlab("Stability South Unit")

#create visual using point
combo_north_geompoint<-Combo_north_biomass%>%
  pivot_longer(c(biomass_PBG_north_mean_mean,biomass_ABG_north),
               names_to = "treatment", values_to = "biom_mean")
ggplot(combo_north_geompoint,aes(RecYear, biom_mean, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))#+
  geom_errorbar(aes(ymax=biom_mean+biomass_PBG_north_sd,
                    ymin=biom_mean-biomass_PBG_north_sd))

combo_south_geompoint<-Combo_south_biomass%>%
  pivot_longer(c(biomass_PBG_south_mean_mean,biomass_ABG_south),
               names_to = "treatment", values_to = "biom_mean")
ggplot(combo_south_geompoint,aes(RecYear, biom_mean, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))

#create visual sd-cv using point
combo_north_geompoint_cv<-Combo_north_biomass%>%
  rename(cv_biomass_PBGNth_M=biomass_PBGNth_cv_M)%>%
  pivot_longer(c(cv_biomass_PBGNth_M,cv_biomass_ABGNth),
               names_to = "treatment", values_to = "biom_cv")
ggplot(combo_north_geompoint_cv,aes(RecYear, biom_cv, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))#+
  #geom_errorbar(aes(ymin=min(cv_biomass_PBGNth),
                   # ymax=max(cv_biomass_PBGNth)))
combo_north_geompoint_sd<-Combo_north_biomass%>%
  pivot_longer(c(biomass_PBGNth_sd_M,biomass_ABGNth_sd),
               names_to = "treatment", values_to = "biom_sd")%>%
  select(biom_sd, treatment, RecYear, biomass_PBGNth_sd_sd)%>%
  distinct()%>%
  mutate(PBG_sd=round(biomass_PBGNth_sd_sd,digits=2))%>%
  mutate(PBG_sd=ifelse(PBG_sd==5.80 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==17.33 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==18.47 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==21.14 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==27.26 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==13.67 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==108.50 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==26.21 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==25.61 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==79.29 & treatment=="biomass_ABGNth_sd",NA,PBG_sd))

ggplot(combo_north_geompoint_sd,aes(RecYear, biom_sd, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=biom_sd+1.96*(PBG_sd),
                    ymin=biom_sd-1.96*(PBG_sd)),width=.2)
#create new dataframe to graph sd result as bargraph####
north_sd_bar<-combo_north_geompoint_sd%>%
  mutate(PBG_confint=1.96*PBG_sd)%>%
  group_by(treatment)%>%
  summarise(SD_biomass=mean(biom_sd, na.rm=T),
            confit_biom=mean(PBG_confint, na.rm=T))
ggplot(north_sd_bar, aes(treatment, SD_biomass, fill=treatment))+
  geom_col(width = .5)+
  geom_errorbar(aes(ymin=SD_biomass-confit_biom,
                    ymax=SD_biomass+confit_biom), width=.1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))
  
  

###graph for the south unit and sd####
combo_south_geompoint_cv<-Combo_south_biomass%>%
  rename(cv_biomass_PBGSth_M=biomass_PBGSth_cv_M)%>%
  pivot_longer(c(cv_biomass_PBGSth_M,cv_biomass_ABGSth),
               names_to = "treatment", values_to = "biom_cv")
ggplot(combo_south_geompoint_cv,aes(RecYear, biom_cv, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"))#+
#geom_errorbar(aes(ymin=min(cv_biomass_PBGNth),
# ymax=max(cv_biomass_PBGNth)))
combo_south_geompoint_sd<-Combo_south_biomass%>%
  pivot_longer(c(biomass_PBGSth_sd_M,biomass_ABGSth_sd),
               names_to = "treatment", values_to = "biom_sd")%>%
  select(biom_sd, treatment, RecYear, biomass_PBGSth_sd_sd)%>%
  distinct()%>%
  mutate(PBG_sd=round(biomass_PBGSth_sd_sd,digits=2))%>%
  mutate(PBG_sd=ifelse(PBG_sd==11.67 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==25.63 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==57.94 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==34.79 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==174.02 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==66.04 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==106.97 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==53.46 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==95.08 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==94.37 & treatment=="biomass_ABGSth_sd",NA,PBG_sd))

ggplot(combo_south_geompoint_sd,aes(RecYear, biom_sd, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=biom_sd+1.96*(PBG_sd),
                    ymin=biom_sd-1.96*(PBG_sd)),width=.2)

#create new dataframe to graph sd result as bargraph####
south_sd_bar<-combo_south_geompoint_sd%>%
  mutate(PBG_confint=1.96*PBG_sd)%>%
  group_by(treatment)%>%
  summarise(SD_biomass=mean(biom_sd, na.rm=T),
            confit_biom=mean(PBG_confint, na.rm=T))
ggplot(south_sd_bar, aes(treatment, SD_biomass, fill=treatment))+
  geom_col(width = .5)+
  geom_errorbar(aes(ymin=SD_biomass-confit_biom,
                    ymax=SD_biomass+confit_biom), width=.1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))

  
#bootstrapping for same plot for each year in each unit
#create an index for the north plots
biomass_north_index<-PBG_north_biomass%>%
  group_by(RecYear)%>%
  mutate(plot_index=1:length(RecYear))
#use sample_n to repeat the dataframe for 1000 iteration
num_bootstrap<-1000
bootstrap_vector<-1:num_bootstrap
PBG_Nth_master_all<-{}
for(BOOT in 1:length(bootstrap_vector)){
  north_rand_key_all<-biomass_north_index%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Transect, Plotnum, plot_index)%>%
    unique(.)%>%
    group_by(RecYear)%>%
    sample_n(600)%>%
    dplyr::select(plot_index,RecYear)%>%
    ungroup()
  biomass_north_ready_all<-biomass_north_index%>%
    right_join(north_rand_key_all, by= c("RecYear", "plot_index"),
               multiple="all")%>%
    mutate(iteration=BOOT)
  PBG_Nth_master_all<-rbind(PBG_Nth_master_all, biomass_north_ready_all)
}

#sample at random for a year to be used as key for other years
PBG_Nth_key<-{}
for(BOOT in 1:length(bootstrap_vector)){
  north_rand_key_2012<-biomass_north_index%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Transect, Plotnum, plot_index)%>%
    unique(.)%>%
    filter(RecYear==2012)%>%
    group_by(RecYear)%>%
    sample_n(200, replace=T)%>%
    dplyr::select(plot_index, RecYear)%>%
    ungroup()
  biomass_north_ready_2012<-biomass_north_index%>%
    right_join(north_rand_key_2012, by= c("RecYear", "plot_index"),
               multiple="all")%>%
    mutate(iteration=BOOT)
  PBG_Nth_key<-rbind(PBG_Nth_key, biomass_north_ready_2012)
}

#select columns needed for matching
PBG_Nth_merge<-PBG_Nth_key%>%
  ungroup()%>%
  select(iteration,plot_index, Plotnum)
#inner join works best-semi join removes duplicates
PBG_Nth_merger<-PBG_Nth_master_all%>%
  inner_join(PBG_Nth_merge, by=c("iteration","Plotnum","plot_index"))

#confirm each year has the same plot id fir each iteration
checking<-PBG_Nth_merger%>%
  filter(iteration==2)
#saving this so I don't have to repeat bootstrap when I close R
write.csv(PBG_Nth_merger, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/PBG_Nth_merger_biomass.csv")
PBG_Nth_merger<-read_csv("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/PBG_Nth_merger_biomass.csv")%>%
  select(-1)
PBG_Nth_biomass_ready<-PBG_Nth_merger%>%
  ungroup()%>%
  group_by(RecYear,iteration)%>%
  summarise(PBG_Nmean=mean(biomass, na.rm=T))%>%
  ungroup()%>%
  group_by(iteration)%>%
  mutate(PBG_Nstab=mean(PBG_Nmean,na.rm=T)/sd(PBG_Nmean))%>%
  ungroup()%>%
  mutate(PBG_Nstab_M=mean(PBG_Nstab, na.rm=T),
         PBG_Nstab_sd=sd(PBG_Nstab))%>%
  left_join(ABG_north_biomass, by="RecYear")%>%
  mutate(zscore_stab=((Stab_ABGNth-PBG_Nstab_M)/PBG_Nstab_sd),
  pval_stabN=2*pnorm(-abs(zscore_stab)))
#visual
ggplot(PBG_Nth_biomass_ready,aes(PBG_Nstab))+
  geom_density(size=2,col="#009E73")+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=Stab_ABGNth), linetype=2,size=2, col="#F0E442")+
  xlab("Stability North Unit")
     
#create an index for the south plots
biomass_south_index<-PBG_south_biomass%>%
  group_by(RecYear)%>%
  mutate(plot_index=1:length(RecYear))
#use sample_n to repeat the dataframe for 1000 iteration
num_bootstrap<-1000
bootstrap_vector<-1:num_bootstrap
PBG_Sth_master_all<-{}
for(BOOT in 1:length(bootstrap_vector)){
  south_rand_key_all<-biomass_south_index%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Transect, Plotnum, plot_index)%>%
    unique(.)%>%
    group_by(RecYear)%>%
    sample_n(600)%>%
    dplyr::select(plot_index,RecYear)%>%
    ungroup()
  biomass_south_ready_all<-biomass_south_index%>%
    right_join(south_rand_key_all, by= c("RecYear", "plot_index"),
               multiple="all")%>%
    mutate(iteration=BOOT)
  PBG_Sth_master_all<-rbind(PBG_Sth_master_all, biomass_south_ready_all)
}

#sample at random for a year to be used as key for other years
PBG_Sth_key<-{}
for(BOOT in 1:length(bootstrap_vector)){
  south_rand_key_2012<-biomass_south_index%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Transect, Plotnum, plot_index)%>%
    unique(.)%>%
    filter(RecYear==2012)%>%
    group_by(RecYear)%>%
    sample_n(200, replace=T)%>%
    dplyr::select(plot_index, RecYear)%>%
    ungroup()
  biomass_south_ready_2012<-biomass_south_index%>%
    right_join(south_rand_key_2012, by= c("RecYear", "plot_index"),
               multiple="all")%>%
    mutate(iteration=BOOT)
  PBG_Sth_key<-rbind(PBG_Sth_key, biomass_south_ready_2012)
}

#select columns needed for matching
PBG_Sth_merge<-PBG_Sth_key%>%
  ungroup()%>%
  select(iteration,plot_index, Plotnum)
#inner join works best-semi join removes duplicates
PBG_Sth_merger<-PBG_Sth_master_all%>%
  inner_join(PBG_Sth_merge, by=c("iteration","Plotnum","plot_index"))

#confirm each year has the same plot id fir each iteration
checking_S<-PBG_Sth_merger%>%
  filter(iteration==2)
#saving this so I don't have to repeat bootstrap when I close R
write.csv(PBG_Sth_merger, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/PBG_Sth_merger_biomass.csv")
PBG_Sth_merger<-read_csv("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/PBG_Sth_merger_biomass.csv")%>%
  select(-1)
PBG_Sth_biomass_ready<-PBG_Sth_merger%>%
  ungroup()%>%
  group_by(RecYear,iteration)%>%
  summarise(PBG_Smean=mean(biomass, na.rm=T))%>%
  ungroup()%>%
  group_by(iteration)%>%
  mutate(PBG_Sstab=mean(PBG_Smean,na.rm=T)/sd(PBG_Smean))%>%
  ungroup()%>%
  mutate(PBG_Sstab_M=mean(PBG_Sstab, na.rm=T),
         PBG_Sstab_sd=sd(PBG_Sstab))%>%
  left_join(ABG_south_biomass, by="RecYear")%>%
  mutate(zscore_stab=((Stab_ABGSth-PBG_Sstab_M)/PBG_Sstab_sd),
         pval_stabS=2*pnorm(-abs(zscore_stab)))
#visual
ggplot(PBG_Sth_biomass_ready,aes(PBG_Sstab))+
  geom_density(size=2, col="#009E73")+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=Stab_ABGSth), linetype=2,size=2, col="#F0E442")+
  xlab("Stability South Unit")



#Bootstrap PBG for unit scale with transect as observational units####
#averaging at transect scale
biomass_t_scale<-biomass_data%>%
  group_by(RecYear, Unit, Watershed, FireGrzTrt, Transect)%>%
  summarise(biomass_t= mean(biomass, na.rm=T))%>%
  ungroup()


#filter for PBG to perform bootstrap
biomass_PBG_t_index_N<-biomass_t_scale%>%
  filter(FireGrzTrt=="PBG" & Unit=="north")%>%
  group_by(RecYear)%>%
  #number the transect to sample from
  mutate(biomass_index= 1:length(RecYear))

biomass_PBG_t_index_S<-biomass_t_scale%>%
  filter(FireGrzTrt=="PBG" & Unit=="south")%>%
  group_by(RecYear)%>%
  #number the transect to sample from
  mutate(biomass_index= 1:length(RecYear))

#extract ABG 
biomass_t_ABGN<-biomass_t_scale%>%
  filter(FireGrzTrt=="ABG" & Unit=="north")%>%
  group_by(RecYear)%>%
  summarise(biomass_ABGNth=mean(biomass_t, na.rm=T),
            biomass_ABGNth_sd=sd(biomass_t),
            cv_biomass_ABGNth=biomass_ABGNth_sd/biomass_ABGNth)
  
biomass_t_ABGS<-biomass_t_scale%>%
  filter(FireGrzTrt=="ABG" & Unit=="south")%>%
  group_by(RecYear)%>%
  summarise(biomass_ABGSth=mean(biomass_t, na.rm=T),
            biomass_ABGSth_sd=sd(biomass_t),
            cv_biomass_ABGSth=biomass_ABGSth_sd/biomass_ABGSth)
 
#North Unit
num_bootstrap<-1000
bootstrap_vector<-1:num_bootstrap
biomass_t_master_N<-{}
for(BOOT in 1:length(bootstrap_vector)){
  biomass_t_rand_N<-biomass_PBG_t_index_N%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Transect, biomass_index)%>%
    unique()%>%
    group_by(RecYear)%>%
    sample_n(4, replace=T)%>%
    dplyr::select(RecYear, biomass_index)%>%
    ungroup()
  biomass_t_ready<-biomass_PBG_t_index_N%>%
    right_join(biomass_t_rand_N, by= c("RecYear", "biomass_index"),
               multiple="all")%>%
    mutate(iteration=BOOT)
  biomass_t_master_N<-rbind(biomass_t_master_N, biomass_t_ready)
  }
#South Unit
biomass_t_master_S<-{}
for(BOOT in 1:length(bootstrap_vector)){
  biomass_t_rand_S<-biomass_PBG_t_index_S%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Transect, biomass_index)%>%
    unique()%>%
    group_by(RecYear)%>%
    sample_n(4, replace=T)%>%
    dplyr::select(RecYear, biomass_index)%>%
    ungroup()
  biomass_t_ready_S<-biomass_PBG_t_index_S%>%
    right_join(biomass_t_rand_S, by= c("RecYear", "biomass_index"),
               multiple="all")%>%
    mutate(iteration=BOOT)
  biomass_t_master_S<-rbind(biomass_t_master_S, biomass_t_ready_S)
}

#calculate mean, SD, CV for PBG bootstrap####
PBG_t_north<-biomass_t_master_N%>%
  group_by(RecYear,iteration)%>%
  summarise(biomass_PBGNth=mean(biomass_t,na.rm=T),
            sd_biomass_PBGNth=sd(biomass_t),
            cv_biomass_PBGNth=sd(biomass_t)/mean(biomass_t, na.rm=T))%>%
  group_by(iteration)%>%
  mutate(stab_PBGNth=mean(biomass_PBGNth, na.rm=T)/sd(biomass_PBGNth))

PBG_t_south<-biomass_t_master_S%>%
  group_by(RecYear,iteration)%>%
  summarise(biomass_PBGSth=mean(biomass_t,na.rm=T),
            sd_biomass_PBGSth=sd(biomass_t),
            cv_biomass_PBGSth=sd(biomass_t)/mean(biomass_t, na.rm=T))%>%
  group_by(iteration)%>%
  mutate(stab_PBGSth=mean(biomass_PBGSth, na.rm=T)/sd(biomass_PBGSth))


#combine PBG and ABG
Combo_Nth_biomass<-PBG_t_north%>%
  left_join(biomass_t_ABGN, by="RecYear")%>%
  group_by(RecYear)%>%
  mutate(biomass_PBGNth_m=mean(biomass_PBGNth, na.rm=T),
         biomass_PBGNth_sd=sd(biomass_PBGNth),
         z_score_NthMean=((biomass_ABGNth-biomass_PBGNth_m)/biomass_PBGNth_sd),
         p_value_NMean=2*pnorm(-abs(z_score_NthMean)),
         biomass_PBGNth_sd_M=mean(sd_biomass_PBGNth, na.rm=T),
         biomass_PBGNth_sd_sd=sd(sd_biomass_PBGNth),
         Z_score_Nsd=((biomass_ABGNth_sd-biomass_PBGNth_sd_M)/biomass_PBGNth_sd_sd),
         pvalue_Nsd=2*pnorm(-abs(Z_score_Nsd)),
         biomass_PBGNth_cv_M=mean(cv_biomass_PBGNth, na.rm=T),
         biomass_PBGNth_cv_sd=sd(cv_biomass_PBGNth),
         Z_score_Ncv=((cv_biomass_ABGNth-biomass_PBGNth_cv_M)/biomass_PBGNth_cv_sd),
         pvalue_Ncv=2*pnorm(-abs(Z_score_Ncv)))

Combo_Sth_biomass<-PBG_t_south%>%
  left_join(biomass_t_ABGS, by="RecYear")%>%
  group_by(RecYear)%>%
  mutate(biomass_PBGSth_m=mean(biomass_PBGSth, na.rm=T),
         biomass_PBGSth_sd=sd(biomass_PBGSth),
         z_score_SthMean=((biomass_ABGSth-biomass_PBGSth_m)/biomass_PBGSth_sd),
         p_value_SMean=2*pnorm(-abs(z_score_SthMean)),
         biomass_PBGSth_sd_M=mean(sd_biomass_PBGSth, na.rm=T),
         biomass_PBGSth_sd_sd=sd(sd_biomass_PBGSth),
         Z_score_Ssd=((biomass_ABGSth_sd-biomass_PBGSth_sd_M)/biomass_PBGSth_sd_sd),
         pvalue_Ssd=2*pnorm(-abs(Z_score_Ssd)),
         biomass_PBGSth_cv_M=mean(cv_biomass_PBGSth, na.rm=T),
         biomass_PBGSth_cv_sd=sd(cv_biomass_PBGSth),
         Z_score_Scv=((cv_biomass_ABGSth-biomass_PBGSth_cv_M)/biomass_PBGSth_cv_sd),
         pvalue_Scv=2*pnorm(-abs(Z_score_Scv)))


#create figures for mean, sd, and cv####
combo_Nth_geompoint<-Combo_Nth_biomass%>%
  pivot_longer(c(biomass_PBGNth_m,biomass_ABGNth),
               names_to = "treatment", values_to = "biom_mean")%>%
  select(treatment,biom_mean,RecYear, biomass_PBGNth_sd)%>%
  distinct()%>%
  mutate(PBG_sd=round(biomass_PBGNth_sd,digits=2))%>%
  mutate(PBG_sd=ifelse(PBG_sd==29.95 & treatment=="biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==39.54 & treatment=="biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==50.22 & treatment=="biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==86.42 & treatment=="biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==77.51 & treatment=="biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==51.36 & treatment=="biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==47.78 & treatment=="biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==48.01 & treatment=="biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==89.08 & treatment=="biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==40.15 & treatment=="biomass_ABGNth",NA,PBG_sd))

ggplot(combo_Nth_geompoint,aes(RecYear, biom_mean, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=biom_mean+1.96*(PBG_sd),
                    ymin=biom_mean-1.96*(PBG_sd)),width=.2)
#sd
combo_Nth_geompoint_sd<-Combo_Nth_biomass%>%
  pivot_longer(c(biomass_PBGNth_sd_M,biomass_ABGNth_sd),
               names_to = "treatment", values_to = "biom_sd")%>%
  select(treatment,biom_sd,RecYear, biomass_PBGNth_sd_sd)%>%
  distinct()%>%
  mutate(PBG_sd=round(biomass_PBGNth_sd_sd,digits=2))%>%
  mutate(PBG_sd=ifelse(PBG_sd==31.24 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==29.44 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==47.37 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==66.38 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==52.24 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==37.39 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==31.68 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==41.61 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==89.31 & treatment=="biomass_ABGNth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==33.09 & treatment=="biomass_ABGNth_sd",NA,PBG_sd))

ggplot(combo_Nth_geompoint_sd,aes(RecYear, biom_sd, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=biom_sd+1.96*(PBG_sd),
                    ymin=biom_sd-1.96*(PBG_sd)),width=.2)

#cv
combo_Nth_geompoint_cv<-Combo_Nth_biomass%>%
  rename()
  pivot_longer(c(biomass_PBGNth_cv_M,cv_biomass_ABGNth),
               names_to = "treatment", values_to = "biom_cv")%>%
  select(treatment,biom_cv,RecYear, biomass_PBGNth_cv_sd)%>%
  distinct()%>%
  mutate(PBG_sd=round(biomass_PBGNth_cv_sd,digits=2))%>%
  mutate(PBG_sd=ifelse(PBG_sd==0.13 & treatment=="cv_biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.09 & treatment=="cv_biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.12 & treatment=="cv_biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.13 & treatment=="cv_biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.09 & treatment=="cv_biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.09 & treatment=="cv_biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.10 & treatment=="cv_biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.08 & treatment=="cv_biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.18 & treatment=="cv_biomass_ABGNth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.08 & treatment=="cv_biomass_ABGNth",NA,PBG_sd))

ggplot(combo_Nth_geompoint_cv,aes(RecYear, biom_cv, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#009E73","#F0E442"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=biom_cv+1.96*(PBG_sd),
                    ymin=biom_cv-1.96*(PBG_sd)),width=.2)

#South Unit
#mean
ggplot(Combo_Sth_biomass,aes(biomass_PBGNth))+
  geom_density(size=.5)+
  facet_wrap(~RecYear, scales = "free")+
  geom_vline(aes(xintercept=biomass_ABGSth), linetype=2,size=.5)

combo_Sth_geompoint<-Combo_Sth_biomass%>%
  pivot_longer(c(biomass_PBGSth_m,biomass_ABGSth),
               names_to = "treatment", values_to = "biom_mean")%>%
  select(treatment,biom_mean,RecYear, biomass_PBGSth_sd)%>%
  distinct()%>%
  mutate(PBG_sd=round(biomass_PBGSth_sd,digits=2))%>%
  mutate(PBG_sd=ifelse(PBG_sd==41.84 & treatment=="biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==83.56 & treatment=="biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==97.24 & treatment=="biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==48.73 & treatment=="biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==131.87 & treatment=="biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==87.64 & treatment=="biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==101.21 & treatment=="biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==80.75 & treatment=="biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==188.49 & treatment=="biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==150.02 & treatment=="biomass_ABGSth",NA,PBG_sd))

ggplot(combo_Sth_geompoint,aes(RecYear, biom_mean, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=biom_mean+1.96*(PBG_sd),
                    ymin=biom_mean-1.96*(PBG_sd)),width=.2)
#sd
combo_Sth_geompoint_sd<-Combo_Sth_biomass%>%
  pivot_longer(c(biomass_PBGSth_sd_M,biomass_ABGSth_sd),
               names_to = "treatment", values_to = "biom_sd")%>%
  select(treatment,biom_sd,RecYear, biomass_PBGSth_sd_sd)%>%
  distinct()%>%
  mutate(PBG_sd=round(biomass_PBGSth_sd_sd,digits=2))%>%
  mutate(PBG_sd=ifelse(PBG_sd==44.06 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==53.06 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==101.33 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==35.77 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==87.16 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==74.87 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==106.34 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==67.31 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==144.23 & treatment=="biomass_ABGSth_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==186.18 & treatment=="biomass_ABGSth_sd",NA,PBG_sd))

ggplot(combo_Sth_geompoint_sd,aes(RecYear, biom_sd, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#F0E442", "#009E73"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=biom_sd+1.96*(PBG_sd),
                    ymin=biom_sd-1.96*(PBG_sd)),width=.2)


#cv
combo_Sth_geompoint_cv<-Combo_Sth_biomass%>%
pivot_longer(c(biomass_PBGSth_cv_M,cv_biomass_ABGSth),
             names_to = "treatment", values_to = "biom_cv")%>%
  select(treatment,biom_cv,RecYear, biomass_PBGSth_cv_sd)%>%
  distinct()%>%
  mutate(PBG_sd=round(biomass_PBGSth_cv_sd,digits=2))%>%
  mutate(PBG_sd=ifelse(PBG_sd==0.14 & treatment=="cv_biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.14 & treatment=="cv_biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.15 & treatment=="cv_biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.05 & treatment=="cv_biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.15 & treatment=="cv_biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.15 & treatment=="cv_biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.20 & treatment=="cv_biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.11 & treatment=="cv_biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.20 & treatment=="cv_biomass_ABGSth",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.26 & treatment=="cv_biomass_ABGSth",NA,PBG_sd))

ggplot(combo_Sth_geompoint_cv,aes(RecYear, biom_cv, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#009E73","#F0E442"),labels=c("ABG","PBG"))+
  geom_errorbar(aes(ymax=biom_cv+1.96*(PBG_sd),
                    ymin=biom_cv-1.96*(PBG_sd)),width=.2)



#bootstrap for stability selecting the same transects each year####
#North Unit
num_bootstrap<-1000
bootstrap_vector<-1:num_bootstrap
biomass_t_master_N_stab<-{}
for(BOOT in 1:length(bootstrap_vector)){
  biomass_t_rand_N_stab<-biomass_PBG_t_index_N%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Transect, biomass_index)%>%
    unique()%>%
    filter(RecYear==2012)%>%
    ungroup()%>%
    sample_n(4, replace=T)%>%
    dplyr::select(biomass_index)%>%
    left_join(biomass_PBG_t_index_N, by="biomass_index", multiple="all")
  biomass_t_ready_stab<-biomass_t_rand_N_stab%>%
    mutate(iteration=BOOT)
  biomass_t_master_N_stab<-rbind(biomass_t_master_N_stab, biomass_t_ready_stab)
}
write.csv(biomass_t_master_N_stab, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/biomass_t_master_N_stab.csv")
#biomass_t_master_N_stab<-read_csv("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/PBG_Sth_merger_biomass.csv")%>%
  select(-1)

#South Unit
num_bootstrap<-1000
bootstrap_vector<-1:num_bootstrap
biomass_t_master_S_stab<-{}
for(BOOT in 1:length(bootstrap_vector)){
  biomass_t_rand_S_stab<-biomass_PBG_t_index_S%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Transect, biomass_index)%>%
    unique()%>%
    filter(RecYear==2012)%>%
    ungroup()%>%
    sample_n(4, replace=T)%>%
    dplyr::select(biomass_index)%>%
    left_join(biomass_PBG_t_index_S, by="biomass_index", multiple="all")
  biomass_t_ready_S_stab<-biomass_t_rand_S_stab%>%
    mutate(iteration=BOOT)
  biomass_t_master_S_stab<-rbind(biomass_t_master_S_stab, biomass_t_ready_S_stab)
}
write.csv(biomass_t_master_S_stab, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/biomass_t_master_N_stab.csv")
#biomass_t_master_S_stab<-read_csv("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/biomass_t_master_S_stab.csv")%>%
  #select(-1)

#calculate stab####
PBG_t_north_same<-biomass_t_master_N_stab%>%
  group_by(RecYear,iteration)%>%
  summarise(biomass_PBGNth=mean(biomass_t,na.rm=T),
            sd_biomass_PBGNth=sd(biomass_t),
            cv_biomass_PBGNth=sd(biomass_t)/mean(biomass_t, na.rm=T))%>%
  group_by(iteration)%>%
  mutate(stab_PBGNth=mean(biomass_PBGNth, na.rm=T)/sd(biomass_PBGNth))

PBG_t_south_same<-biomass_t_master_S_stab%>%
  group_by(RecYear,iteration)%>%
  summarise(biomass_PBGSth=mean(biomass_t,na.rm=T),
            sd_biomass_PBGSth=sd(biomass_t),
            cv_biomass_PBGSth=sd(biomass_t)/mean(biomass_t, na.rm=T))%>%
  group_by(iteration)%>%
  mutate(stab_PBGSth=mean(biomass_PBGSth, na.rm=T)/sd(biomass_PBGSth))

#ABG
biomass_t_ABGN_same<-biomass_t_ABGN%>%
  mutate(stab_ABGN= mean(biomass_ABGNth)/sd(biomass_ABGNth))
biomass_t_ABGS_same<-biomass_t_ABGS%>%
  mutate(stab_ABGS= mean(biomass_ABGSth)/sd(biomass_ABGSth))

#combine PBG and ABG
Combo_Nth_biomass_same<-PBG_t_north_same%>%
  left_join(biomass_t_ABGN_same, by="RecYear")%>%
  group_by(iteration)%>%
  mutate(stab_PBGNth=mean(biomass_PBGNth)/sd(biomass_PBGNth))%>%
  ungroup()%>%
  mutate(stab_PBGNth_m=mean(stab_PBGNth),
         stab_PBGNth_sd=sd(stab_PBGNth),
         stab_Z_score=((stab_ABGN-stab_PBGNth_m)/stab_PBGNth_sd),
         pvalue_Ncv=2*pnorm(-abs(stab_Z_score)))

Combo_Sth_biomass_same<-PBG_t_south_same%>%
  left_join(biomass_t_ABGS_same, by="RecYear")%>%
  group_by(iteration)%>%
  mutate(stab_PBGSth=mean(biomass_PBGSth)/sd(biomass_PBGSth))%>%
  ungroup()%>%
  mutate(stab_PBGSth_m=mean(stab_PBGSth),
         stab_PBGSth_sd=sd(stab_PBGSth),
         stab_Z_score=((stab_ABGS-stab_PBGSth_m)/stab_PBGSth_sd),
         pvalue_Ncv=2*pnorm(-abs(stab_Z_score)))

#visual
ggplot(Combo_Nth_biomass_same,aes(stab_PBGNth))+
  geom_density(size=2,col="#009E73")+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=stab_ABGN), linetype=2,size=2, col="#F0E442")+
  xlab("Stability North Unit")

ggplot(Combo_Sth_biomass_same,aes(stab_PBGSth))+
  geom_density(size=2,col="#009E73")+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=stab_ABGS), linetype=2,size=2, col="#F0E442")+
  xlab("Stability South Unit")


#Things to do, calculate PBG stability across unit and sd_cv####
#compare growing season precip with sd####
ppt_data <- read.csv("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/PhD Wyoming_One Drive/PHD Wyoming/Thesis/PBG synthesis/annual and growing season precipitation_knzHQ_1983-2021.csv")
#relationship between precipitation and biomass sd
Combo_north_biomass_ppt<-Combo_north_biomass%>%
  select(RecYear, iteration, sd_biomass_PBGNth, biomass_PBGNth_sd_M, biomass_ABGNth_sd)%>%
  left_join(ppt_data, by="RecYear")

ggplot(Combo_north_biomass_ppt, aes(x=gs_ppt, y=biomass_ABGNth_sd)) +
  geom_point(size=3) + geom_smooth(method="lm", se=F) #+ facet_wrap(~yrsins_fire)

ppt_sd_model<-lm(biomass_ABGNth_sd~gs_ppt, data=Combo_north_biomass_ppt)
summary(ppt_sd_model)
ppt_sd_model2<-lm(biomass_PBGNth_sd_M~gs_ppt, data=Combo_north_biomass_ppt)
summary(ppt_sd_model2)
#using all iteration gives similar estimates as using averages, rsquared is better with average
#significant positive relationship for both ABG and PBG. slope estimates(0.221, 0.191 respectively)

#combine unit and examine if stability differs####
#averaging at plot scale
biomass_plot_scale<-biomass_data%>%
  group_by(RecYear, Unit, Watershed, FireGrzTrt, Transect, Plotnum)%>%
  summarise(biomass_plot= mean(biomass, na.rm=T),
            biomass_sd=sd(biomass))%>%
  ungroup()

#ABG stability
biomass_ABG_combine_stab<-biomass_plot_scale%>%
  filter(FireGrzTrt=="ABG")%>%
  group_by(RecYear)%>%
  summarise(biomass_ABG=mean(biomass_plot, na.rm=T))%>%
  ungroup()%>%
  mutate(biomass_ABG_stab=mean(biomass_ABG)/sd(biomass_ABG))


#filter for PBG to perform bootstrap
biomass_PBG_plot_index<-biomass_plot_scale%>%
  filter(FireGrzTrt=="PBG")%>%
  group_by(RecYear)%>%
  #number the plots to sample from
  mutate(biomass_PBG_plot_index= 1:length(RecYear))



biomass_master_stab<-{}
for(BOOT in 1:length(bootstrap_vector)){
  biomass_rand_stab<-biomass_PBG_plot_index%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Transect, biomass_PBG_plot_index)%>%
    unique()%>%
    filter(RecYear==2012)%>%
    ungroup()%>%
    sample_n(400, replace=T)%>%
    dplyr::select(biomass_PBG_plot_index)%>%
    left_join(biomass_PBG_plot_index, by="biomass_PBG_plot_index", multiple="all")
  biomass_ready_stab<-biomass_rand_stab%>%
    mutate(iteration=BOOT)
  biomass_master_stab<-rbind(biomass_master_stab, biomass_ready_stab)
}

biomass_combine_stab<-biomass_master_stab%>%
  select(RecYear, iteration, biomass_plot)%>%
  group_by(RecYear, iteration)%>%
  summarise(biomass_mean=mean(biomass_plot,na.rm=T))%>%
  ungroup()%>%
  group_by(iteration)%>%
  mutate(biomass_PBG_stab=mean(biomass_mean, na.rm=T)/sd(biomass_mean))%>%
  ungroup()%>%
  mutate(PBG_stab_M=mean(biomass_PBG_stab),
         PBG_stab_sd=sd(biomass_PBG_stab))%>%
  left_join(biomass_ABG_combine_stab, by="RecYear")%>%
  mutate(zscore=(biomass_ABG_stab-PBG_stab_M)/PBG_stab_sd)%>%
  mutate(pval=2*pnorm(-abs(zscore)))
write.csv(biomass_combine_stab, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/biomass_stability_combined_unit.csv")
biomass_combine_stab<-read.csv("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/biomass_stability_combined_unit.csv")
ggplot(biomass_combine_stab,aes(biomass_PBG_stab))+
  geom_density(size=2,col="#009E73")+
  #facet_grid("RecYear")+
  geom_vline(aes(xintercept=biomass_ABG_stab), linetype=2,size=2, col="#F0E442")+
  xlab("Biomass stability")+xlim(0,4.5) 

#represent stability as a bar graph####
biom_stability_df<-biomass_combine_stab%>%
  select(PBG_stab_M,PBG_stab_sd, biomass_ABG_stab)%>%
  distinct()%>%
  rename(PBG=PBG_stab_M, ABG=biomass_ABG_stab)%>%
  pivot_longer(cols=c(1,3), names_to = "treatment", values_to = "stability")%>%
  mutate(PBG_stab_sd=case_when(treatment=="ABG"~"NA",
                               treatment=="PBG"~"0.15"),
         confit=as.numeric(PBG_stab_sd)*1.96)

ggplot(biom_stability_df, aes(treatment, fill=treatment))+
  geom_bar(stat = "identity",aes(y=stability),width = 0.25)+
  geom_errorbar(aes(ymin=stability-confit,
                    ymax=stability+confit), width=0.05)+
  scale_fill_manual(values=c( "#F0E442", "#009E73"))


#combined unit sd, mean, cv####
#ABG mean, cv, sd
biomass_ABG_combine_mean_sd<-biomass_plot_scale%>%
  filter(FireGrzTrt=="ABG")%>%
  group_by(RecYear)%>%
  summarise(biomass_ABG=mean(biomass_plot, na.rm=T),
            biomass_ABG_sd=sd(biomass_plot),
            biomass_ABG_cv=biomass_ABG_sd/biomass_ABG)


biomass_master_mean_sd_cv<-{}
for(BOOT in 1:length(bootstrap_vector)){
  biomass_rand_mean<-biomass_PBG_plot_index%>%
    dplyr::select(RecYear, Unit, Watershed, FireGrzTrt, Transect, biomass_PBG_plot_index)%>%
    unique()%>%
    group_by(RecYear)%>%
    sample_n(400, replace=T)%>%
    dplyr::select(RecYear,biomass_PBG_plot_index)%>%
    left_join(biomass_PBG_plot_index, by=c("biomass_PBG_plot_index","RecYear"), multiple="all")
  biomass_ready_mean<-biomass_rand_mean%>%
    mutate(iteration=BOOT)
  biomass_master_mean_sd_cv<-rbind(biomass_master_mean_sd_cv, biomass_ready_mean)
}

#calculate zscore and p vals####
biomass_combine_mean_sd<-biomass_master_mean_sd_cv%>%
  select(RecYear, iteration, biomass_plot)%>%
  group_by(RecYear, iteration)%>%
  summarise(biomass_PBG_mean=mean(biomass_plot,na.rm=T),
            biomass_PBG_sd=sd(biomass_plot),
            biomass_PBG_cv=biomass_PBG_sd/biomass_PBG_mean)%>%
  ungroup()%>%
  group_by(RecYear)%>%
  summarise(biom_PBG_MM=mean(biomass_PBG_mean),
         biom_PBG_M_sd=sd(biomass_PBG_mean),
         biom_PBG_sd_M=mean(biomass_PBG_sd),
         biom_PBG_sd_sd=sd(biomass_PBG_sd),
         biom_PBG_cv_M=mean(biomass_PBG_cv),
         biom_PBG_cv_sd=sd(biomass_PBG_cv))%>%
  left_join(biomass_ABG_combine_mean_sd, by="RecYear")%>%
  mutate(zscore_M=(biomass_ABG-biom_PBG_MM)/biom_PBG_M_sd,
         pval_M=2*pnorm(-abs(zscore_M)),
         zscore_sd=(biomass_ABG_sd-biom_PBG_sd_M)/biom_PBG_sd_sd,
         pval_sd=2*pnorm(-abs(zscore_sd)),
         zscore_cv=(biomass_ABG_cv-biom_PBG_cv_M)/biom_PBG_cv_sd,
         pval_cv=2*pnorm(-abs(zscore_cv)))
#to avoid rerunning the bootstrap
write.csv(biomass_combine_mean_sd, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/biomass_mean_sd_combined_unit.csv")
#read in the saved file
biomass_combine_mean_sd<-read.csv("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/biomass_mean_sd_combined_unit.csv")
#create visuals####
#biomass average
combo_biomass_geompoint_M<-biomass_combine_mean_sd%>%
  pivot_longer(c(biom_PBG_MM,biomass_ABG),
               names_to = "treatment", values_to = "biom_M")%>%
  select(treatment,biom_M,RecYear, biom_PBG_M_sd)%>%
  distinct()%>%
  mutate(PBG_sd=round(biom_PBG_M_sd,digits=2))%>%
  mutate(PBG_sd=ifelse(PBG_sd==5.38 & treatment=="biomass_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==10.42 & treatment=="biomass_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==13.65 & treatment=="biomass_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==15.49 & treatment=="biomass_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==22.41 & treatment=="biomass_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==12.93 & treatment=="biomass_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==19.42 & treatment=="biomass_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==15.27 & treatment=="biomass_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==26.66 & treatment=="biomass_ABG",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==23.69 & treatment=="biomass_ABG",NA,PBG_sd))

ggplot(combo_biomass_geompoint_M,aes(RecYear, biom_M, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#009E73","#F0E442"),labels=c("PBG","ABG"))+
  geom_errorbar(aes(ymax=biom_M+1.96*(PBG_sd),
                    ymin=biom_M-1.96*(PBG_sd)),width=.1)



#sd biomass
combo_biomass_geompoint_sd<-biomass_combine_mean_sd%>%
  pivot_longer(c(biom_PBG_sd_M,biomass_ABG_sd),
               names_to = "treatment", values_to = "biom_sd")%>%
  select(treatment,biom_sd,RecYear, biom_PBG_sd_sd)%>%
  distinct()%>%
  mutate(PBG_sd=round(biom_PBG_sd_sd,digits=2))%>%
  mutate(PBG_sd=ifelse(PBG_sd==7.72 & treatment=="biomass_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==17.88 & treatment=="biomass_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==37.11 & treatment=="biomass_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==24.05 & treatment=="biomass_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==111.68 & treatment=="biomass_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==42.35 & treatment=="biomass_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==78.49 & treatment=="biomass_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==35.72 & treatment=="biomass_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==65.54 & treatment=="biomass_ABG_sd",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==63.85 & treatment=="biomass_ABG_sd",NA,PBG_sd))

ggplot(combo_biomass_geompoint_sd,aes(RecYear, biom_sd, col=treatment))+
  geom_point(size=3)+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#009E73","#F0E442"),labels=c("PBG","ABG"))+
  geom_errorbar(aes(ymax=biom_sd+1.96*(PBG_sd),
                    ymin=biom_sd-1.96*(PBG_sd)),width=.1)+
  scale_x_continuous(breaks = 2012:2021)
#wrangle for bargraph figure
combo_biomass_bargraph_sd<-combo_biomass_geompoint_sd%>%
  group_by(treatment)%>%
  summarise(biomass_sd=mean(biom_sd, na.rm=T),
            sd=mean(PBG_sd, na.rm=T))%>%
  mutate(treatment=ifelse(treatment=="biom_PBG_sd_M","PBG",treatment),
         treatment=ifelse(treatment=="biomass_ABG_sd","ABG",treatment))
#create bargraph figure for SD
ggplot(combo_biomass_bargraph_sd,aes(x=treatment,fill=treatment))+
  geom_bar(stat = "identity",aes(y=biomass_sd),width = 0.25)+
  geom_errorbar(aes(ymin=biomass_sd-sd,
                    ymax=biomass_sd+sd),width=0.05,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#009E73")) 

#cv biomass
combo_biomass_geompoint_cv<-biomass_combine_mean_sd%>%
  pivot_longer(c(biom_PBG_cv_M,biomass_ABG_cv),
               names_to = "treatment", values_to = "biom_cv")%>%
  select(treatment,biom_cv,RecYear, biom_PBG_cv_sd)%>%
  distinct()%>%
  mutate(PBG_sd=round(biom_PBG_cv_sd,digits=2))%>%
  mutate(PBG_sd=ifelse(PBG_sd==0.03 & treatment=="biomass_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.04 & treatment=="biomass_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.07 & treatment=="biomass_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.03 & treatment=="biomass_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.16 & treatment=="biomass_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.09 & treatment=="biomass_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.17 & treatment=="biomass_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.06 & treatment=="biomass_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.09 & treatment=="biomass_ABG_cv",NA,PBG_sd),
         PBG_sd=ifelse(PBG_sd==0.10 & treatment=="biomass_ABG_cv",NA,PBG_sd))

ggplot(combo_biomass_geompoint_cv,aes(RecYear, biom_cv, col=treatment))+
  geom_point()+
  geom_path(aes(as.numeric(RecYear)))+
  scale_colour_manual(values=c( "#009E73","#F0E442"),labels=c("PBG","ABG"))+
  geom_errorbar(aes(ymax=biom_cv+1.96*(PBG_sd),
                    ymin=biom_cv-1.96*(PBG_sd)),width=.2)


#calculate mean of year since burn####
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
#merge with biomass data 
yrs_biomass<-biomass_data%>%
  left_join(YrSinceFire_key, by="year_watershed")%>%
  group_by(RecYear,Unit, Watershed,yrsins_fire,FireGrzTrt,Transect,Plotnum)%>%
  summarise(biomass=mean(biomass, na.rm=T))

yrs_biomass$RecYear=as.factor(yrs_biomass$RecYear)
yrs_biomass$Plotnum=as.factor(yrs_biomass$Plotnum)
yrs_biomass$yrsins_fire=as.factor(yrs_biomass$yrsins_fire)
yrs_biomass$Unit=as.factor(yrs_biomass$Unit)
yrs_biomass$Watershed=as.factor(yrs_biomass$Watershed)
yrs_biomass$Transect=as.factor(yrs_biomass$Transect)

#mixed anova year since burn####
yrs_biomass_model<-lmer(log(biomass)~yrsins_fire*RecYear+(1|Unit/Watershed/Transect),
                        data=yrs_biomass)#issingular due to watershed
summary(yrs_biomass_model)
check_model(yrs_biomass_model)
anova(yrs_biomass_model)
qqnorm(resid(yrs_biomass_model))
check_normality(yrs_biomass_model)

#multiple comparison
testInteractions(yrs_biomass_model, pairwise="yrsins_fire", fixed="RecYear",
                 adjustment="BH")
#comparison of fixed effects
testInteractions(yrs_biomass_model, pairwise="yrsins_fire")
#using mean estimate to create figure 
yrs_interact<-interactionMeans(yrs_biomass_model)
#replacing spaces in column names with underscore 
names(yrs_interact)<-str_replace_all(names(yrs_interact), " ","_")
#df for visuals from model estimates
yrs_interact_viz<-yrs_interact%>%
  mutate(yrs_interact_bt_mean=exp(adjusted_mean),
         yrs_interact_bt_upper=exp(adjusted_mean+SE_of_link),
         yrs_interact_bt_lower=exp(adjusted_mean-SE_of_link))%>%
  #including 2011 for graphical purpose-so the graph will align with other data that has 2011
  bind_rows(data_2011=(tibble(RecYear="2011", yrsins_fire="PBG0")))%>%
  ungroup()
  


#visual year since burn####
ggplot(yrs_interact_viz,aes(as.numeric(RecYear),col=yrsins_fire))+
  geom_point(size=3,aes(y=yrs_interact_bt_mean))+
  geom_path(aes(x=as.numeric(RecYear),y=yrs_interact_bt_mean),linewidth=1)+
  geom_errorbar(aes(x=as.numeric(RecYear),ymin=yrs_interact_bt_lower,
                    ymax=yrs_interact_bt_upper),width=0.2,linetype=1)+
  scale_colour_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))+
  scale_y_continuous(limits=c(0,750))

#average across years for simplification

# yrs_interact_bar<-yrs_interact_viz%>%
#   group_by(yrsins_fire)%>%
#   summarise(biomass_mean=mean(yrs_interact_bt_mean, na.rm=T),
#             se_upper=mean(yrs_interact_bt_upper,na.rm=T),
#             se_lower=mean(yrs_interact_bt_lower,na.rm=T))
# ggplot(yrs_interact_bar,aes(x=yrsins_fire,fill=yrsins_fire))+
#   geom_bar(stat = "identity",aes(y=biomass_mean),width = 0.5)+
#   geom_errorbar(aes(ymin=se_lower,
#                     ymax=se_upper),width=0.2,linetype=1)+
#   scale_fill_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))

#figure for just year since fire after accounting for year covariate
yrs_wo_interact<-interactionMeans(yrs_biomass_model, factors = "yrsins_fire")
names(yrs_wo_interact)<-str_replace_all(names(yrs_wo_interact), " ","_")
#df for visuals from model estimates
yrs_wo_interact_viz<-yrs_wo_interact%>%
  mutate(yrs_interact_bt_mean=exp(adjusted_mean),
         yrs_interact_bt_upper=exp(adjusted_mean+SE_of_link),
         yrs_interact_bt_lower=exp(adjusted_mean-SE_of_link))
#visual
ggplot(yrs_wo_interact_viz,aes(x=yrsins_fire,fill=yrsins_fire))+
  geom_bar(stat = "identity",aes(y=yrs_interact_bt_mean))+
  geom_errorbar(aes(ymin=yrs_interact_bt_lower,
                    ymax=yrs_interact_bt_upper),width=0.2,linetype=1)+
  scale_fill_manual(values=c( "#F0E442", "#994F00", "#999999", "#0072B2"))
#load in precipitation data####
precip_data<-read.csv("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/PhD Wyoming_One Drive/PHD Wyoming/Thesis/PBG synthesis/Precipitation_1982_2023.csv")%>%
  #separate date into components
  separate_wider_delim(RecDate, delim="/", names=c("mth","day","year" ))
#convert precipitation to numeric values
precip_data$ppt<-as.numeric(precip_data$ppt)

annual_ppt_data<-precip_data%>%
  select(mth,day,year,ppt)%>%
  filter(year%in%2011:2021)%>%
  group_by(year)%>%
  mutate(annual_ppt=sum(ppt, na.rm=T))%>%
  filter(mth%in% 4:9)%>%
  group_by(year)%>%
  mutate(gsppt=sum(ppt, na.rm=T))%>%
  select(year,annual_ppt,gsppt)%>%
  rename(RecYear=year)%>%
  distinct()
annual_ppt_data$RecYear=as.factor(annual_ppt_data$RecYear)


#precipitation figure
#growing season ppt
ggplot(annual_ppt_data, aes(x=RecYear))+
  geom_bar(alpha=0.3,stat = "identity",aes(y=gsppt),width = 0.5)

#annual precipitation
ggplot(annual_ppt_data, aes(x=RecYear))+
  geom_bar(alpha=0.3,stat = "identity",aes(y=gsppt),width = 0.5)




###correlating spatial heterogeneity with temporal stability####
biomass_spatial_hete<-combo_biomass_geompoint_sd%>%
  group_by(treatment)%>%
  summarise(spatial_hetero=mean(biom_sd, na.rm=T))%>%
  mutate(variable= "plant_biomass")%>%
  #change characters in a variable
  mutate(treatment=case_when(treatment=="biom_PBG_sd_M"~"PBG",
                             treatment=="biomass_ABG_sd"~"ABG"))
#load in stability values
biomass_temp_stability<-read.csv("C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/biomass_stability_combined_unit.csv")%>%
  select(biomass_ABG_stab,PBG_stab_M)%>%
  distinct()%>%
  pivot_longer(cols=1:2, names_to = "treatment", values_to = "stability")%>%
  mutate(treatment=case_when(treatment=="biomass_ABG_stab"~"ABG",
                             treatment=="PBG_stab_M"~"PBG"))

#combine stability and heterogeneity
biom_stab_vs_heter<-biomass_temp_stability%>%
  left_join(biomass_spatial_hete, by="treatment")%>%
  mutate(spat_heter=spatial_hetero/sd(spatial_hetero),
         stab=stability/sd(stability))
write.csv(biom_stab_vs_heter, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/biom_stab_vs_heter.csv")

#separating by unit stability vs Spatial hetero####
#correlation btw stability and spatial heterogeneity
#wrangle north unit
north_stab_spat<-Combo_north_biomass%>%
  select(Stab_ABGNth, biomass_ABGNth_sd, biomass_PBGNth_sd_M,RecYear,stab_PBGNm)%>%
  distinct()%>%
  pivot_longer(cols=2:3, names_to = "treatment", values_to = "sd")%>%
  mutate(treatment=case_when(treatment=="biomass_ABGNth_sd"~"ABG",
                             treatment=="biomass_PBGNth_sd_M"~"PBG"),
         unit="North")%>%
  group_by(treatment)%>%
  mutate(sd=mean(sd,na.rm=T))%>%
  select(-RecYear)%>%
  distinct()%>%
  pivot_longer(cols=1:2, names_to = "treat", values_to = "stability")%>%
  ungroup()%>%
  mutate(key=c("Y","N","N","Y"))%>%
  filter(key!="N")%>%
  select(-treat,-key)

#wrangle south unit
south_stab_spat<-Combo_south_biomass%>%
  select(Stab_ABGSth, biomass_ABGSth_sd, biomass_PBGSth_sd_M,RecYear,stab_PBGSm)%>%
  distinct()%>%
  pivot_longer(cols=2:3, names_to = "treatment", values_to = "sd")%>%
  mutate(treatment=case_when(treatment=="biomass_ABGSth_sd"~"ABG",
                             treatment=="biomass_PBGSth_sd_M"~"PBG"),
         unit="South")%>%
  group_by(treatment)%>%
  mutate(sd=mean(sd,na.rm=T))%>%
  select(-RecYear)%>%
  distinct()%>%
  pivot_longer(cols=1:2, names_to = "treat", values_to = "stability")%>%
  ungroup()%>%
  mutate(key=c("Y","N","N","Y"))%>%
  filter(key!="N")%>%
  select(-treat,-key)

#combine south and north
biom_stab_spat_unit<-south_stab_spat%>%
  bind_rows(north_stab_spat)%>%
  mutate(sp="plant")%>%
  rename(spat_hetero=sd)%>%
  mutate(sd_scaled=(spat_hetero-mean(spat_hetero))/sd(spat_hetero),
         stab=stability-mean(stability))
write.csv(biom_stab_spat_unit, "C:/Users/JAAJOWELE/OneDrive - UNCG/UNCG PHD/Writing/2024_PBG_figures/biom_stab_spat_unit.csv")

#regression
biom_stab_spat_lm<-lm(stability~sd_scaled, data=biom_stab_spat)
summary(biom_stab_spat_lm)#not significant

ggplot(biom_stab_spat, aes(sd_scaled, stability))+
  geom_point(size=5,aes(col=unit))+
  geom_smooth(method="lm")

