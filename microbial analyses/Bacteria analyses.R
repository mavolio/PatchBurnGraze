################################################################################
##  bacterial community analysis.R: Sequence data from GMDR treatments
##
##  Author: Kimberly Komatsu
##  Date created: March 28, 2024

##modifed by Olivia Arbogast and Meghan Avolio Fall 2024
################################################################################


#set working directory
#Olivia
setwd('C:\\Users\\livyh\\OneDrive - Johns Hopkins\\Avolio_Lab\\PBG')

#Meghan
setwd('C:\\Users\\mavolio2\\OneDrive - Johns Hopkins\\Olivia Avolio_Lab\\PBG\\IMR Return\\KomatsuJul2024with2023_analysis\\')

# install and load package
# if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#
library(devtools)
install.packages('BiocManager')
library(BiocManager)
BiocManagerinstall(biomformat)
BiocManagerinstall(GenomeInfoDB)
BiocManagerinstall(phyloseq)
BiocManagerinstall(SummarizedExperiment)
BiocManagerinstall(SingleCellExperiment)
BiocManagerinstall(TreeSummarizedExperiment)
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
#library(performance)
#library(PerformanceAnalytics)
library(car)
library(vegan)
library(grid)
library(tidyverse)
library(codyn)
library(lme4)
library(lmerTest)

'%notin%' <- negate('%in%')

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  

###setting the graph look
theme_set(theme_bw(12))


##### Import Data #####
bactread<-read_qza("V6V8\\deblur_output\\deblur_table_final.qza")
sample<-read_delim('V6V8\\metadata.txt') %>% 
  rename(SampleID='#SampleID') %>% 
  select(-"...8") %>% 
  rename(PBGTrt='Burn-Treatment', 
         YrsBurn='Years-Since-Burn') %>% 
  mutate(unid=paste(Watershed, Transect, Plot, sep="_"))

bactPD<-as.data.frame(read_qza("V6V8\\diversity\\faith_pd_vector.qza")$data) %>% 
  rename(SampleID=V1, FaithPD=V2) %>% 
  left_join(sample)

bactdat<-as.data.frame(bactread$data) %>% 
  rownames_to_column("Feature.ID") %>% 
  pivot_longer('320-22':'415-23', names_to = 'SampleID', values_to = 'abund') %>% 
  filter(abund>0)

totabund<-bactdat %>% 
  group_by(SampleID) %>% 
  summarize(totabund=sum(abund))

relabund<-bactdat %>% 
  left_join(totabund) %>% 
  mutate(relabund=abund/totabund)

##I can't get this to work
parse_taxonomy(bactread$data)

#richness and evenness
richeven<-community_structure(relabund, abundance.var = 'relabund', replicate.var = 'SampleID') %>% 
  left_join(sample)

#is there a richness difference?

hist(richeven$richness)#very normal
hist(richeven$Evar)#very normal

#diff in richness by year and by PBG but no interaction
#no diff for years since burn
m1<-lmer(richness~PBGTrt*Year + (1|unid), data=richeven)
anova(m1, ddf="Kenward-Roger")
m2<-lmer(richness~YrsBurn*Year + (1|unid), data=richeven)
anova(m2, ddf="Kenward-Roger")

#diff in evenness by year and PBG but no interaction
#no diff years since burn
m3<-lmer(Evar~PBGTrt*Year + (1|unid), data=richeven)
anova(m3, ddf="Kenward-Roger")
m4<-lmer(Evar~YrsBurn*Year + (1|unid), data=richeven)
anova(m4, ddf="Kenward-Roger")

#diff in faithsPD by year and PBG but no interaction
#no diff years since burn
m5<-lmer(FaithPD~PBGTrt*Year + (1|unid), data=bactPD)
anova(m5, ddf="Kenward-Roger")
m6<-lmer(FaithPD~YrsBurn*Year + (1|unid), data=bactPD)
anova(m6, ddf="Kenward-Roger")

#graph richness by year
ggplot(data=barGraphStats(data=richeven, variable="richness", byFactorNames="Year"), aes(x=Year, y=mean)) +
  geom_point (size=2)+
  xlab('Year') + ylab("Bacteria Sp. Richness") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1)
#graph evenness by Trt
ggplot(data=barGraphStats(data=richeven, variable="Evar", byFactorNames="PBGTrt"), aes(x=PBGTrt, y=mean)) +
  geom_point (size=2)+
  xlab('PBGTrt') + ylab("Bacteria Sp. Evenness") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1)
#graph evenness by Year
ggplot(data=barGraphStats(data=richeven, variable="Evar", byFactorNames="Year"), aes(x=Year, y=mean)) +
  geom_point (size=2)+
  xlab('Year') + ylab("Bacteria Sp. Evenness") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1)

#graph Faith's PD by Trt
ggplot(data=barGraphStats(data=bactPD, variable="FaithPD", byFactorNames="PBGTrt"), aes(x=PBGTrt, y=mean)) +
  geom_point (size=2)+
  xlab('PBGTrt') + ylab("Bacteria Phylogenetic Div") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1)
#graph Faith's PD by Year
ggplot(data=barGraphStats(data=bactPD, variable="FaithPD", byFactorNames="Year"), aes(x=Year, y=mean)) +
  geom_point (size=2)+
  xlab('Year') + ylab("Bacteria Phylogenetic Div") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=1)

#do multivariate
wide_dat<-relabund %>% 
  select(-abund, -totabund) %>% 
  left_join(sample) %>% 
  pivot_wider(names_from = Feature.ID, values_from = relabund, values_fill = 0)
sites<-wide_dat[,1:9]

mds<-metaMDS(wide_dat[9:19014], trymax = 100)

scores<-as.data.frame(mds$points) %>%
  bind_cols(sites)

ggplot(data=scores, aes(x=MDS1, y=MDS2, color=Watershed, shape=Year))+
  geom_point(size=3)
#watershed C1SB is just different from the other - this will overwhelm any other signal




#####################Kim's older code below###################


##### Multi-Variate Analysis #####
#Distance matrix
bunifraqDistance <- as.matrix(read_qza('KomatsuV6V8\\diversity\\weighted_unifrac_distance_matrix.qza')$data)
funifraqDistance <- as.matrix(read_qza('KomatsuITS2\\diversity\\weighted_unifrac_distance_matrix.qza')$data)

bunifracDistanceTrt <- as.data.frame(bunifraqDistance) %>%
  rownames_to_column("SampleID") %>%
  mutate(type='bacteria') %>%
  left_join(metadata) %>%
  filter(!grepl('ctrl',SampleID))

bunifraqDistanceMatrix <- bunifracDistanceTrt %>%
  select(-type:-livestock_util_2021) %>%
  column_to_rownames(var="SampleID") %>%
  as.matrix()

bunifraqEnv <- bunifracDistanceTrt %>%
  select(SampleID, type:livestock_util_2021)

# #plotting PCOA
# unifracPCOA<-read_qza("KomatsuV6V8\\diversity\\weighted_unifrac_pcoa_results.qza")
# 
# unifracPCOATrt<-unifracPCOA$data$Vectors %>%
#   select(SampleID, PC1, PC2) %>%
#   left_join(pdTrt) %>%
#   mutate(site2=ifelse(site=='TB', 'WY',
#                ifelse(site=='FK', 'MT', 'NA')))
# 
# ggplot(subset(unifracPCOATrt, !is.na(site2)), aes(x=PC1, y=PC2, color=as.factor(rainfall_reduction), shape=grazing_treatment)) +
#   facet_wrap(~site2, scales='free') +
#   geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
#   theme_q2r() +
#   scale_size_continuous(name="Faith's PD") +
#   scale_color_manual(name="Rainfall Reduction", values=droughtColor) +
#   scale_shape_discrete(name="Grazing Treatment")
# 
# # ggsave('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_GMDR_soils\\figures\\Fig_FKTB_bacterial_PCOA.png', width=14, height=7, units='in', dpi=300, bg='white')


##### PERMANOVA #####
bSVrelative <- as.data.frame(apply(bSVs$data, 2, function(x) x/sum(x)*100)) %>% #convert to percent
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key="SampleID", value="Abundance") %>%
  left_join(btaxonomy2) %>%
  group_by(SampleID, Kingdom, Phylum) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup()  %>%
  filter(!grepl('ctrl',SampleID),
         !(Phylum %in% c('Campilobacterota', 'FCPU426', 'MBNT15'))) %>%
  left_join(metadata, multiple='all') %>%
  select(-type) %>%
  unique() %>%
  mutate(taxa=paste(Kingdom, Phylum, sep='::')) %>%
  select(-Kingdom, -Phylum) %>%
  filter(Year>2018) %>%
  spread(key=taxa, value=Abundance, fill=0)
# write.csv(bSVrelative, 'SVrelative_bacteria_20240717.csv', row.names=F)
chart.Correlation(bSVrelative[,15:44], histogram=T, pch=19)

fSVrelative <- as.data.frame(apply(fSVs$data, 2, function(x) x/sum(x)*100)) %>% #convert to percent
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key="SampleID", value="Abundance") %>%
  left_join(ftaxonomy2) %>%
  group_by(SampleID, Kingdom, Phylum) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup()  %>%
  filter(!grepl('ctrl',SampleID),
         !(Phylum %in% c('Aphelidiomycota'))) %>%
  left_join(metadata, multiple='all') %>%
  select(-type) %>%
  unique() %>%
  mutate(taxa=paste(Kingdom, Phylum, sep='::')) %>%
  select(-Kingdom, -Phylum) %>%
  filter(Year>2018) %>%
  spread(key=taxa, value=Abundance, fill=0)
# write.csv(fSVrelative, 'SVrelative_fungi_20240717.csv', row.names=F)
chart.Correlation(fSVrelative[,15:28], histogram=T, pch=19)

# #run first 4 lines of each of above taxa code
# SVrelative <- rbind(bSVrelative, fSVrelative) %>%
#   group_by(SampleID, Kingdom, Phylum) %>%
#   summarise(Abundance=sum(Abundance)) %>%
#   ungroup()  %>%
#   filter(!grepl('ctrl',SampleID),
#          !(Phylum %in% c('Campilobacterota', 'FCPU426', 'MBNT15', 'Aphelidiomycota'))) %>%
#   left_join(metadata, multiple='all') %>%
#   select(-type) %>%
#   unique() %>%
#   mutate(taxa=paste(Kingdom, Phylum, sep='::')) %>%
#   select(-Kingdom, -Phylum) %>%
#   filter(Year>2018) %>%
#   spread(key=taxa, value=Abundance, fill=0)
# # write.csv(SVrelative, 'SVrelative_alltaxa_20240717.csv', row.names=F)

set.seed(123)

#FK permanova - bacterial
fkSVmatrix <- bSVrelative %>%
  filter(site=='FK') %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

fkSVenv <- bSVrelative %>%
  filter(site=='FK') %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaFK <- adonis2(formula = fkSVmatrix ~ Year*rainfall_reduction*grazing_treatment,
                       data=fkSVenv,
                       permutations = 999, method = "bray")
print(permanovaFK) #sig rainfall effect (marginal interaction with year)

# ggplot(subset(unifracPCOATrt), aes(x=PC1, y=PC2, color=as.factor(rainfall_reduction), shape=grazing_treatment)) +
#   facet_wrap(~site2, scales='free') +
#   geom_point() + #alpha controls transparency and helps when points are overlapping
#   theme_q2r() +
#   scale_size_continuous(name="Faith's PD") +
#   scale_color_manual(name="Rainfall Reduction", values=droughtColor) +
#   scale_shape_discrete(name="Grazing Treatment")

#FK permanova - fungal
fkSVmatrix <- fSVrelative %>%
  filter(site=='FK') %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

fkSVenv <- fSVrelative %>%
  filter(site=='FK') %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaFK <- adonis2(formula = fkSVmatrix ~ Year*rainfall_reduction*grazing_treatment,
                       data=fkSVenv,
                       permutations = 999, method = "bray")
print(permanovaFK) #sig rainfall*year effect


#FK permanova - fungal 2019
fkSVmatrix <- fSVrelative %>%
  filter(site=='FK', Year==2019) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

fkSVenv <- fSVrelative %>%
  filter(site=='FK', Year==2019) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaFK <- adonis2(formula = fkSVmatrix ~ rainfall_reduction,
                       data=fkSVenv,
                       permutations = 999, method = "bray")
print(permanovaFK) #sig rainfall effect

#FK permanova - fungal 2020
fkSVmatrix <- fSVrelative %>%
  filter(site=='FK', Year==2020) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

fkSVenv <- fSVrelative %>%
  filter(site=='FK', Year==2020) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaFK <- adonis2(formula = fkSVmatrix ~ rainfall_reduction,
                       data=fkSVenv,
                       permutations = 999, method = "bray")
print(permanovaFK)

#FK permanova - fungal 2021
fkSVmatrix <- fSVrelative %>%
  filter(site=='FK', Year==2021) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

fkSVenv <- fSVrelative %>%
  filter(site=='FK', Year==2021) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaFK <- adonis2(formula = fkSVmatrix ~ rainfall_reduction,
                       data=fkSVenv,
                       permutations = 999, method = "bray")
print(permanovaFK) 

#FK permanova - fungal 2022
fkSVmatrix <- fSVrelative %>%
  filter(site=='FK', Year==2022) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

fkSVenv <- fSVrelative %>%
  filter(site=='FK', Year==2022) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaFK <- adonis2(formula = fkSVmatrix ~ rainfall_reduction,
                       data=fkSVenv,
                       permutations = 999, method = "bray")
print(permanovaFK) #sig rainfall effect 

ggplot(subset(unifracPCOATrt, type=='fungi' & site=='FK'), aes(x=PC1, y=PC2, color=as.factor(rainfall_reduction))) +
  facet_wrap(~year, scales='free') +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_size_continuous(name="Faith's PD") +
  scale_color_manual(name="Rainfall Reduction", values=droughtColor) +
  scale_shape_discrete(name="Grazing Treatment")



#TB permanova - bacterial
tbSVmatrix <- bSVrelative %>%
  filter(site=='TB') %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- bSVrelative %>%
  filter(site=='TB') %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ Year*rainfall_reduction*grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) #significant grazing * year effect

#TB permanova - bacterial 2019
tbSVmatrix <- bSVrelative %>%
  filter(site=='TB', Year==2019) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- bSVrelative %>%
  filter(site=='TB', Year==2019) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) #significant grazing

#TB permanova - bacterial 2020
tbSVmatrix <- bSVrelative %>%
  filter(site=='TB', Year==2020) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- bSVrelative %>%
  filter(site=='TB', Year==2020) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) 

#TB permanova - bacterial 2021
tbSVmatrix <- bSVrelative %>%
  filter(site=='TB', Year==2021) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- bSVrelative %>%
  filter(site=='TB', Year==2021) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) #significant grazing

#TB permanova - bacterial 2022
tbSVmatrix <- bSVrelative %>%
  filter(site=='TB', Year==2022) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- bSVrelative %>%
  filter(site=='TB', Year==2022) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) 


#TB permanova - fungal
tbSVmatrix <- fSVrelative %>%
  filter(site=='TB') %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- fSVrelative %>%
  filter(site=='TB') %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ Year*rainfall_reduction*grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) #significant grazing * year effect


#TB permanova - fungal 2019
tbSVmatrix <- fSVrelative %>%
  filter(site=='TB', Year==2019) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- fSVrelative %>%
  filter(site=='TB', Year==2019) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ rainfall_reduction,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) #sig rainfall effect 

permanovaTB <- adonis2(formula = tbSVmatrix ~ grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) #no effect

#TB permanova - fungal 2020
tbSVmatrix <- fSVrelative %>%
  filter(site=='TB', Year==2020) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- fSVrelative %>%
  filter(site=='TB', Year==2020) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ rainfall_reduction,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) 

permanovaTB <- adonis2(formula = tbSVmatrix ~ grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) 

#TB permanova - fungal 2021
tbSVmatrix <- fSVrelative %>%
  filter(site=='TB', Year==2021) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- fSVrelative %>%
  filter(site=='TB', Year==2021) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ rainfall_reduction,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) 

permanovaTB <- adonis2(formula = tbSVmatrix ~ grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) 

#TB permanova - fungal 2022
tbSVmatrix <- fSVrelative %>%
  filter(site=='TB', Year==2022) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- fSVrelative %>%
  filter(site=='TB', Year==2022) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ rainfall_reduction,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) 

permanovaTB <- adonis2(formula = tbSVmatrix ~ grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) 

ggplot(subset(unifracPCOATrt, type=='fungi' & site=='TB' & year>2018), aes(x=PC1, y=PC2, color=as.factor(rainfall_reduction))) +
  facet_wrap(~year, scales='free') +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_size_continuous(name="Faith's PD") +
  scale_color_manual(name="Rainfall Reduction", values=droughtColor) +
  scale_shape_discrete(name="Grazing Treatment")

ggplot(subset(unifracPCOATrt, type=='fungi' & site=='TB' & year>2018), aes(x=PC1, y=PC2, color=as.factor(grazing_treatment))) +
  facet_wrap(~year, scales='free') +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_size_continuous(name="Faith's PD") +
  scale_color_manual(name="Grazing Intensity", values=grazingColor) +
  scale_shape_discrete(name="Grazing Treatment")

ggplot(subset(unifracPCOATrt, type=='fungi' & site=='TB' & year>2018), aes(x=PC1, y=PC2, color=as.factor(grazing_treatment))) +
  # facet_wrap(~year, scales='free') +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_size_continuous(name="Faith's PD") +
  scale_color_manual(name="Grazing Intensity", values=grazingColor) +
  scale_shape_discrete(name="Grazing Treatment")




### fungal taxa specific responses
fSVsToPlot <- fSVrelative %>% #find the average abundance of a SV
  rename('Fungi::Other'='Fungi::NA') %>% 
  pivot_longer(cols=c('Fungi::Ascomycota':'Fungi::Olpidiomycota'), names_to='taxa', values_to='abundance') %>% 
  group_by(site, plot, Year, taxa) %>% 
  summarize(abundance=sum(abundance)) %>% 
  ungroup() %>% 
  group_by(taxa) %>% 
  summarize(mean_abundance=mean(abundance)) %>% 
  ungroup() %>% 
  arrange(desc(mean_abundance)) %>%
  top_n(20, mean_abundance) %>%
  pull(taxa) #extract only the names from the table

metadata2 <- metadata %>% 
  select(-type) %>% 
  unique()

fSVphylaHeatMap <- fSVrelative %>%
  rename('Fungi::Other'='Fungi::NA') %>%
  pivot_longer(cols=c('Fungi::Ascomycota':'Fungi::Olpidiomycota'), names_to='taxa', values_to='abundance') %>% 
  group_by(site, plot, Year, taxa) %>% 
  summarize(abundance=sum(abundance)) %>% 
  ungroup() %>% 
  # filter(taxa %in% SVsToPlot) %>% 
  left_join(metadata2) %>%
  mutate(NormAbundance=log10(abundance+0.0001))  # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
# mutate(Feature=paste(Feature.ID, Taxon)) %>%
# mutate(Feature=gsub("[kpcofgs]__", "", Feature)) # trim out leading text from taxonomy string

# write.csv(SVphylaHeatMap, 'C:\\Users\\kjkomatsu\\OneDrive - UNCG\\GMDR\\soils ms derived data\\GMDR_microbial_abundance_2019-2022.csv', row.names=F)


fSVphylaHeatMap$grazing_treatment=factor(fSVphylaHeatMap$grazing_treatment,levels=c('destock', 'stable', 'heavy'))

#fungal by drought TB
ggplot(data=subset(fSVphylaHeatMap, site=='TB' & grepl('Fungi',taxa)), aes(x=as.factor(rainfall_reduction), y=taxa, fill=NormAbundance)) +
  geom_tile() +
  # facet_grid(~`Year`, scales="free_x") +
  # theme_q2r() +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)") +
  scale_y_discrete(limits=rev) +
  # scale_x_discrete(labels=c('destock', 'supplement', 'maintain')) +
  xlab('Rainfall Reduction (%)') + ylab('Taxa')
# ggsave('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_GMDR_soils\\figures\\Fig8_TB_fungalRainfall_heatmap.png', width=12, height=8, units='in', dpi=300, bg='white')

#fungal by grazing TB
ggplot(data=subset(fSVphylaHeatMap, site=='TB' & grepl('Fungi',taxa)), aes(x=as.factor(grazing_treatment), y=taxa, fill=NormAbundance)) +
  geom_tile() +
  # facet_grid(~`Year`, scales="free_x") +
  # theme_q2r() +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)") +
  scale_y_discrete(limits=rev) +
  scale_x_discrete(labels=c('destock', 'supplement', 'maintain')) +
  xlab('Grazing Treatment') + ylab('Taxa')
# ggsave('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_GMDR_soils\\figures\\Fig8_TB_fungalGrazing_heatmap.png', width=12, height=8, units='in', dpi=300, bg='white')


#fungal by rainfall FK
ggplot(data=subset(fSVphylaHeatMap, site=='FK' & grepl('Fungi',taxa)), aes(x=as.factor(rainfall_reduction), y=taxa, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~`Year`) +
  # theme_q2r() +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)") +
  scale_y_discrete(limits=rev) +
  # scale_x_discrete(labels=c('destock', 'supplement', 'maintain')) +
  xlab('Rainfall Reduction (%)') + ylab('Taxa')
# ggsave('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_GMDR_soils\\figures\\FigS15_FK_fungalRainfall_heatmap.png', width=12, height=8, units='in', dpi=300, bg='white')



### bacterial taxa specific responses
bSVsToPlot <- bSVrelative %>% #find the average abundance of a SV
  # rename('Fungi::Other'='Fungi::NA') %>% 
  pivot_longer(cols=c('d__Bacteria::Abditibacteriota':'d__Bacteria::WPS-2'), names_to='taxa', values_to='abundance') %>% 
  group_by(site, plot, Year, taxa) %>% 
  summarize(abundance=sum(abundance)) %>% 
  ungroup() %>% 
  group_by(taxa) %>% 
  summarize(mean_abundance=mean(abundance)) %>% 
  ungroup() %>% 
  arrange(desc(mean_abundance)) %>%
  top_n(20, mean_abundance) %>%
  pull(taxa) #extract only the names from the table

metadata2 <- metadata %>% 
  select(-type) %>% 
  unique()

bSVphylaHeatMap <- bSVrelative %>%
  # rename('Fungi::Other'='Fungi::NA') %>%
  pivot_longer(cols=c('d__Bacteria::Abditibacteriota':'d__Bacteria::WPS-2'), names_to='taxa', values_to='abundance') %>% 
  group_by(site, plot, Year, taxa) %>% 
  summarize(abundance=sum(abundance)) %>% 
  ungroup() %>% 
  # filter(taxa %in% SVsToPlot) %>% 
  left_join(metadata2) %>%
  mutate(NormAbundance=log10(abundance+0.0001))  # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
# mutate(Feature=paste(Feature.ID, Taxon)) %>%
# mutate(Feature=gsub("[kpcofgs]__", "", Feature)) # trim out leading text from taxonomy string

# write.csv(SVphylaHeatMap, 'C:\\Users\\kjkomatsu\\OneDrive - UNCG\\GMDR\\soils ms derived data\\GMDR_microbial_abundance_2019-2022.csv', row.names=F)


bSVphylaHeatMap$grazing_treatment=factor(bSVphylaHeatMap$grazing_treatment,levels=c('destock', 'stable', 'heavy'))

#fungal by grazing TB
ggplot(data=subset(bSVphylaHeatMap, site=='TB'), aes(x=as.factor(grazing_treatment), y=taxa, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~`Year`) +
  # theme_q2r() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)") +
  scale_y_discrete(limits=rev) +
  scale_x_discrete(labels=c('destock', 'supplement', 'maintain')) +
  xlab('Rainfall Reduction (%)') + ylab('Taxa')
# ggsave('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_GMDR_soils\\figures\\FigS13_TB_bacterialGrazing_heatmap.png', width=12, height=8, units='in', dpi=300, bg='white') 