#### Libraries ####

library(readxl)
library(codyn)
library(vegan)
library(tidyverse)
library(nlme)
library(emmeans)
library(gridExtra)
library(ggplot2)
library(boot)
library(purrr)


#### CSV read ####

#Contains Belowground Abundance ID information
Abundance_ID_Belowground <- read_excel("Abundance + ID Data.xlsx") %>% 
  rename(Life.Stage="Life Stage")

#Contains time sense burned information for each watershed
BurnInfo <- read_excel("YearsSenseBurned.xlsx")


#### Data Cleanup ####

#Selecting the columns needed
Abundance_ID_Belowground_Clean <- Abundance_ID_Belowground[, c("WS","Trans", "Dist.", "Phylum", "Class", "Order", "Morphospp", "Life.Stage","Count")]

#Adding a new column to tag "1" with ABG and "3" with PBG. 
Abundance_ID_Belowground_Clean$Treatment <- ifelse(grepl("1", Abundance_ID_Belowground_Clean$WS), "ABG", "PBG")

#Typo corrections for morphospecies
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Earthworm Cacoon", "Earthworm Cocoon", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub(" Light Brown Diptera", "Light Brown Diptera", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Big Brown", "Big Brown Beetle", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Big Brown Adult", "Big Brown Beetle", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Big Brown Beetle Adult", "Big Brown Beetle", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Big Brown Beetle Beetle", "Big Brown Beetle", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Black Two Segment", "Black Two Segments", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Black Two Segmentss", "Black Two Segments", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Black Two Segmentss", "Black Two Segments", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Coffee Larvae Beetle", "Coffee Beetle Larvae", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Coffee Beetle", "Coffee Beetle Larvae", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Coffee Beetle Larvae Larvae", "Coffee Beetle Larvae", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Coffee Larvae Bettle", "Coffee Beetle Larvae", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Creme Colored", "Creme Larvae", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Light Brown", "Light Brown Diptera", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Light Brown Diptera Larvae", "Light Brown Diptera", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Light Brown Millipede ", "Light Brown Millipede", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Light Brown Diptera Milipede", "Light Brown Millipede", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Light Brown Diptera Millipede ", "Light Brown Millipede", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Light Brown Diptera Millipede", "Light Brown Millipede", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Mid Brown Stick Larvae", "Midbrown Stick Larvae", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Parasitic Worn", "Parasitic Worm", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Creme Larvae ", "Creme Larvae", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Creme Larvae", "Creme Larvae", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Light Brown Diptera Diptera", "Light Brown Diptera Larvae", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Light Brown Diptera Diptera Diptera Larvae", "Light Brown Diptera Larvae", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Light Brown Diptera Diptera Shrimp Like", "Light Brown Diptera Larvae", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Light Brown Diptera Larvae Larvae", "Light Brown Diptera Larvae", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Light Brown Diptera Shrimp Like", "Light Brown Diptera Larvae", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Light Brown Diptera", "Light Brown Diptera Larvae", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Light Brown Diptera Larvae Larvae", "Light Brown Diptera Larvae", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Large Golden Shrimp", "Creme Larvae", Abundance_ID_Belowground_Clean$Morphospp)
Abundance_ID_Belowground_Clean$Morphospp <- gsub("Small Light Colored Back End", "Small Light Colored Black End", Abundance_ID_Belowground_Clean$Morphospp)

#Typo corrections for watersheds 
Abundance_ID_Belowground_Clean$WS <- gsub("C35A","C3SA", Abundance_ID_Belowground_Clean$WS) 
Abundance_ID_Belowground_Clean$WS <- gsub("C35B","C3SB", Abundance_ID_Belowground_Clean$WS)
Abundance_ID_Belowground_Clean$WS <- gsub("C3CB", "C3C", Abundance_ID_Belowground_Clean$WS)
Abundance_ID_Belowground_Clean$WS <- gsub("C3D", "C3B", Abundance_ID_Belowground_Clean$WS)
Abundance_ID_Belowground_Clean$WS <- gsub("CSC", "C3C", Abundance_ID_Belowground_Clean$WS)

#Merge columns to create a new morpospp and sample column,
Abundance_ID_Belowground_Reformated <- Abundance_ID_Belowground_Clean %>% 
  mutate(Morphospp = paste(Phylum, Class, Order, Morphospp, sep = "_")) %>% 
  unite("Sample", c("WS", "Trans", "Dist.", "Treatment"), sep = "_", remove = FALSE) %>%
  mutate(Sample = str_remove(Sample, "m"))

#Remove NAs
Abundance_ID_Belowground_Reformated <- na.omit(Abundance_ID_Belowground_Reformated)

#### CV Count Graph Prep  ####

#Calculating CV for PBG and ABG
ABG_CV <- Abundance_ID_Belowground_Reformated %>%
  filter(WS %in% c("C1A", "C1SB")) %>%
  group_by(WS) %>%
  summarize(cv_count = sd(`Count`) / mean(`Count`) * 100)

PBG_CV <- Abundance_ID_Belowground_Reformated %>%
  filter(WS %in% c("C3A", "C3B", "C3C", "C3SA", "C3SB")) %>%
  group_by(WS) %>%
  summarize(cv_count = sd(`Count`) / mean(`Count`) * 100)

#Assign PBG and ABG in a treatment column

ABG_CV2 <- ABG_CV %>% 
  mutate(`Treatment` = "ABG") %>% 
  select(-WS)

PBG_CV2 <- PBG_CV %>% 
  mutate(`Treatment` = "PBG")

#Average the CV

ABG_CV_Average <- ABG_CV2 %>%
  group_by(Treatment) %>%
  summarize(Average_CV = mean(cv_count))

PBG_CV_Average <- PBG_CV2 %>%
  group_by(Treatment) %>%
  summarize(Average_CV = mean(cv_count))

#Combine the averages into one file

combined_cv_count <- rbind(ABG_CV_Average, PBG_CV_Average)



#### Prep for Total Count Stats and Graph ####

Abundance_Stats <- group_by(Abundance_ID_Belowground_Reformated, Morphospp, Sample) %>% 
  summarize(Count = sum(Count, na.rm = TRUE))

Abundance_Stats$Treatment <- ifelse(grepl("ABG", Abundance_Stats$Sample), "ABG", 
                                    ifelse(grepl("PBG", Abundance_Stats$Sample), "PBG", NA))

#Adding Block for use in modeling
total_counts <- aggregate(Count ~ Sample + Treatment, data = Abundance_Stats, FUN = sum) %>% 
  separate(Sample, into = c("WS", "Trans", "Dist"), sep = "_", extra = "drop") %>% mutate(Block = ifelse(grepl("S", WS), "North", "South"))

#Average counts of ABG and PBG for graphing
avg_counts <- aggregate(Count ~ Treatment, data = total_counts, FUN = mean) 


### Average Total Count Graph ###

#Making an Error Bar

error_df <- total_counts %>%
  group_by(Treatment) %>%
  summarize(mean_value = mean(Count),
            se = sd(Count) / sqrt(n()))

#Average Count Graph


#### Count Graphs ####

#Average count PBG vs ABG
ggplot(avg_counts, aes(x = Treatment, y = Count, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Treatment", y = "Average Count") +
  geom_errorbar(aes(ymin = Count - ifelse(Treatment == "ABG", 1.726931, 1.390563),
                    ymax = Count + ifelse(Treatment == "ABG", 1.726931, 1.390563)),
                width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("blue", "red")) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 40),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 40),
    axis.ticks.y = element_line(size = 1),
    legend.position = "none"
  ) +
  geom_text(
    data = NULL,
    aes(x = -Inf, y = Inf),
    label = "A",
    size = 10,
    hjust = 0,
    vjust = 1,
    color = "black",
    show.legend = FALSE
  )

#CV Count Graph 
ggplot(combined_cv_count, aes(x = Treatment, y = `Average_CV`, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("ABG" = "blue", "PBG" = "red")) +
  guides(fill = FALSE) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 40),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 40),
        axis.ticks.y = element_line(size = 1)) +
  labs(y = "Count Average CV")

#### Total Count Stats ####
total_counts2 <- full_join(BurnInfo, total_counts) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_") 


TotalcountModel <- #stores the model output into a named list
  lme(Count ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burning)
      data = total_counts2, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|Block) #this would be where you'd say north or south unit (which should be a variable in the data frame)
anova.lme(TotalcountModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova
emmeans(TotalcountModel, pairwise~TreatmentSB, adjust="tukey") #this gives you contrasts (means and confidence intervals) for each possible pairwise comparison of treatments to know whether they are different or not (overlapping confidence intervals means not different)
#### Prep for Richness + Evenness Graphs & Stats ####

#Getting  community data and tagging ABG and PBG
commMetrics <- Abundance_Stats %>%  community_structure(abundance.var='Count', replicate.var='Sample')
commMetrics$Treatment <- ifelse(grepl("C1", commMetrics$Sample), "ABG", "PBG") 

#Adding Blocks
commMetrics2 <- commMetrics %>% mutate(block = ifelse(grepl("S", Sample), "North", "South")) %>% 
  separate(col = "Sample", into = c("WS","Trans","Dist.","Trea"), sep="_")

#Joining Burn information with Treatment

Joined <- full_join(BurnInfo, commMetrics2) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_") 


#### Calculating CV Richness ####

ABG_CV_R <- Joined %>%
  filter(WS %in% c("C1A", "C1SB")) %>%
  group_by(WS) %>%
  summarize(cv_richness = sd(`richness`) / mean(`richness`) * 100)

# Calculate the CV of weights for each combined PBG cluster
PBG_CV_R <- Joined%>%
  filter(WS %in% c("C3A", "C3B", "C3C", "C3SA", "C3SB")) %>%
  group_by(WS) %>%
  summarize(cv_richness = sd(`richness`) / mean(`richness`) * 100)

# Add the average values to the dataframes
ABG_CV2_R <- ABG_CV_R %>% 
  mutate(`Treatment` = "ABG") %>% 
  select(-WS)

PBG_CV2_R <- PBG_CV_R %>% 
  mutate(`Treatment` = "PBG")

ABG_CV_Average_R <- ABG_CV2_R %>%
  group_by(Treatment) %>%
  summarize(Average_CV = mean(cv_richness))

# Calculate the average weight CV for PBG
PBG_CV_Average_R <- PBG_CV2_R %>%
  group_by(Treatment) %>%
  summarize(Average_CV = mean(cv_richness))



combined_cv_richness <- rbind(ABG_CV_Average_R, PBG_CV_Average_R)

#### Calculating CV Evenness ####

Joined_1 <- na.omit(Joined)

ABG_CV_E <- Joined_1 %>%
  filter(WS %in% c("C1A", "C1SB")) %>%
  group_by(WS) %>%
  summarize(cv_evenness = sd(`Evar`) / mean(`Evar`) * 100)

# Calculate the CV of weights for each combined PBG cluster
PBG_CV_E <- Joined_1 %>%
  filter(WS %in% c("C3A", "C3B", "C3C", "C3SA", "C3SB")) %>%
  group_by(WS) %>%
  summarize(cv_evenness = sd(`Evar`) / mean(`Evar`) * 100)

# Add the average values to the dataframes
ABG_CV2_E <- ABG_CV_E %>% 
  mutate(`Treatment` = "ABG") %>% 
  select(-WS)

PBG_CV2_E <- PBG_CV_E %>% 
  mutate(`Treatment` = "PBG")

ABG_CV_Average_E <- ABG_CV2_E %>%
  group_by(Treatment) %>%
  summarize(Average_CV = mean(cv_evenness))

# Calculate the average weight CV for PBG
PBG_CV_Average_E <- PBG_CV2_E %>%
  group_by(Treatment) %>%
  summarize(Average_CV = mean(cv_evenness))



combined_cv_evenness <- rbind(ABG_CV_Average_E, PBG_CV_Average_E)

#CV Evenness graph

ggplot(combined_cv_evenness, aes(x = Treatment, y = `Average_CV`, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("ABG" = "blue", "PBG" = "red")) +
  guides(fill = FALSE) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 40),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 40),
        axis.ticks.y = element_line(size = 1)) +
  labs(y = "Evenness Average CV")
#### Richness Graphs ####

#ABG vs PBG through time 
ggplot(Joined, aes(x = TreatmentSB, y = richness, fill= Trea)) +
  geom_boxplot() +
  labs(
    x = "Treatment",
    y = "Richness") +
  scale_fill_manual(values = c("blue", "red")) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 40),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 40),
        axis.ticks.y = element_line(size = 1)) +
  guides(fill = FALSE)

#CV Richness graph

ggplot(combined_cv_richness, aes(x = Treatment, y = `Average_CV`, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("ABG" = "blue", "PBG" = "red")) +
  guides(fill = FALSE) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 40),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 40),
        axis.ticks.y = element_line(size = 1)) +
  labs(y = "Richness Average CV")

#### Richness Stats ####
richModel <- #stores the model output into a named list
  lme(richness ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = Joined, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(richModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova
#
emmeans(richModel, pairwise~TreatmentSB, adjust="tukey") #this gives you contrasts (means and confidence intervals) for each possible pairwise comparison of treatments to know whether they are different or not (overlapping confidence intervals means not different)

#### Evenness Graphs ####
ggplot(Joined, aes(x = TreatmentSB, y = Evar, fill= Trea)) +
  geom_boxplot() +
  labs(
    x = "Treatment",
    y = "Evenness") +
  scale_fill_manual(values = c("blue", "red")) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 40),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 40),
        axis.ticks.y = element_line(size = 1)) +
  guides(fill = FALSE)



#### Evenness Stats ####
EvenModel <- #stores the model output into a named list
  lme(Evar ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = na.omit(Joined), #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(EvenModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova
emmeans(EvenModel, pairwise~TreatmentSB, adjust="tukey") #this gives you contrasts (means and confid


#### Density Plots ####

# Richness Density plot
richness <- ggplot(Joined, aes(x = richness, color = interaction(Trea, block), linetype = interaction(Trea, block))) +
  geom_density() +
  labs(title = "Belowground Invertebrate Richness Density Plot",
       x = "Richness",
       y = "Density") +
  scale_color_manual(values = rep(c("blue", "red", "blue","red"), 2)) +
  scale_linetype_manual(values = rep(c("solid", "dashed"), each = 2)) +
  theme_minimal()

# Evenness Density plot

evenness <- ggplot(Joined, aes(x = Evar, color = interaction(Trea, block), linetype = interaction(Trea, block))) +
  geom_density() +
  labs(title = "Belowground Invertebrate Evenness Density Plot",
       x = "Evenness",
       y = "Density") +
  scale_color_manual(values = rep(c("blue", "red", "blue", "red"), 2)) +
  scale_linetype_manual(values = rep(c("solid", "dashed"), each = 2)) +
  theme_minimal()

# Count Density Plot

counts <- ggplot(total_counts, aes(x = Count, color = interaction(Treatment, Block), linetype = interaction(Treatment, Block))) +
  geom_density() +
  labs(title = "Belowground Invertebrate Count Density Plot",
       x = "Count",
       y = "Density") +
  scale_color_manual(values = rep(c("blue", "red", "blue", "red"), 2)) +
  scale_linetype_manual(values = rep(c("solid", "dashed"), each = 2)) +
  theme_minimal()

grid.arrange(richness, evenness, counts, Weight) #Weight comes from weight script!

#### multivariate community response - PERMANOVA and NMDS ####

### by Years since burn

#Combing count with burn info

CountGraph <- full_join(Abundance_ID_Belowground_Reformated, BurnInfo) %>% 
  unite("TreatmentSB",c("Treatment","SB"), sep="_")

CountGraph$Treatment <- ifelse(grepl("C1", CountGraph$Sample), "ABG", "PBG") 

# PERMANOVA

abundanceWide <- CountGraph %>% 
  mutate(block = ifelse(grepl("S", Sample), "North", "South")) %>%
  select(Sample, block, WS, Trans, Dist., Treatment, TreatmentSB, Morphospp, Count) %>%
  group_by(Sample, block, WS, Trans, Dist., Treatment, TreatmentSB, Morphospp) %>% 
  summarise(Count = sum(Count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = 'Morphospp', values_from = 'Count', values_fill = list(Count = 0)) 

abundanceWide <- abundanceWide %>% 
  mutate(sum = rowSums(abundanceWide[, c(8:139)], na.rm = TRUE)) %>% 
  filter(sum > 0, !(Sample %in% c('C1A_A_38_ABG', 'C3SA_D_16_PBG', 'C3SA_C_38_PBG')))

###Defined in itself?
#abundanceWide <- CountGraph  %>% 
#  mutate(block = ifelse(grepl("S", Sample), "North", "South")) %>%
#  select(Sample, block, WS, Trans, Dist., Treatment, TreatmentSB, Morphospp, Count) %>%
#  group_by(Sample, block, WS, Trans, Dist., Treatment, TreatmentSB, Morphospp) %>% 
#  summarise(Count=sum(Count)) %>% 
#  ungroup() %>% 
#  pivot_wider(names_from='Morphospp', values_from='Count', values_fill=list(Count=0)) %>% 
#  mutate(sum=rowSums(abundanceWide[,c(8:57)], na.rm=TRUE)) %>% 
#  filter(sum>0, Sample!='C1A_A_38_ABG') #PROBLEM: Check why two of the samples have nothing in them, is this real?
###IMPORTANT: removing sample from C1A_A_38_ABG, which is a big outlier because of high values of endogenic worms and brown shrimp-like beetles, which no other plots had

print(permanova <- adonis2(formula = abundanceWide[,8:139]~TreatmentSB, data=abundanceWide, permutations=999, method="bray"))
#F=1.5242, df=3,49, p=0.042 

#betadisper
veg <- vegdist(abundanceWide[,8:140], method = "bray")
dispersion <- betadisper(veg, abundanceWide$TreatmentSB)
permutest(dispersion, pairwise=TRUE, permutations=999) 
#F=0.773 , df=3,49, p=0.53

BC_Data <- metaMDS(abundanceWide[,8:140])
sites <- 1:nrow(abundanceWide)
BC_Meta_Data <- abundanceWide[,1:7]
plot(BC_Data$points, col=as.factor(BC_Meta_Data$TreatmentSB))
ordiellipse(BC_Data, groups = as.factor(BC_Meta_Data$TreatmentSB), kind = "sd", display = "sites", label = T)

veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100)
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#Generate ellipses
BC_NMDS = data.frame(MDS1 = BC_Data$points[,1], MDS2 = BC_Data$points[,2], group=BC_Meta_Data$TreatmentSB)
BC_NMDS_Graph <- cbind(BC_Meta_Data, BC_NMDS)
BC_Ord_Ellipses<-ordiellipse(BC_Data, BC_Meta_Data$TreatmentSB, display = "sites",
                             kind = "se", conf = 0.95, label = T)

BC_Ellipses <- data.frame()
#Generate ellipses points
for(g in levels(as.factor(BC_NMDS$group))){
  BC_Ellipses <- rbind(BC_Ellipses,
                       cbind(as.data.frame(with(BC_NMDS[BC_NMDS$group==g,], 
                                                veganCovEllipse(BC_Ord_Ellipses[[g]]$cov,
                                                                BC_Ord_Ellipses[[g]]$center,
                                                                BC_Ord_Ellipses[[g]]$scale))),
                             group=g))
}

#Plot the data from BC_NMDS_Graph, where x=MDS1 and y=MDS2, make an ellipse based on "group"
NMDS_Years_Since_Burned <- ggplot(BC_NMDS_Graph, aes(x=MDS1, y=MDS2, color=group,linetype = group, shape = group)) +
  geom_point(size=6) + 
  geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color="", linetype = "", shape = "") +
  scale_colour_manual(values=c("blue", "#7A2021", "#C21A09", "#FF0800"), name = "",
                      labels=c('ABG', 'PBG 0', 'PBG 1', 'PBG 2')) +
  scale_linetype_manual(values = c("solid", "twodash", "twodash", "twodash"), name = "",
                        labels=c('ABG', 'PBG 0', 'PBG 1', 'PBG 2')) +
  scale_shape_manual(values = c(19, 19, 17, 15),
                     labels=c('ABG', 'PBG 0', 'PBG 1', 'PBG 2')) +
  xlab("NMDS1") + 
  ylab("NMDS2") + 
  theme_bw() +
  theme(axis.text.x=element_text(size=24, color = "black"), 
        axis.text.y = element_text(size = 24, color = "black"), 
        axis.title.x = element_text(size = 24, color = 'black'),
        axis.title.y = element_text(size = 24, color = 'black'),
        legend.text = element_text(size = 24),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  ) +
  annotate("text", x=min(BC_NMDS_Graph$MDS1), y=max(BC_NMDS_Graph$MDS2),
           label="F=1.52, p=0.0420", size=6, hjust=0, vjust=1)


#F=1.5242, df=3,49, p=0.042 
#export at 1500x1000

### by watershed
# PERMANOVA
print(permanova <- adonis2(formula = abundanceWide[,8:139]~Treatment, data=abundanceWide, permutations=999, method="bray"))
#F=1.495, df=1,51, p=0.122

#betadisper
veg <- vegdist(abundanceWide[,8:139], method = "bray")
dispersion <- betadisper(veg, abundanceWide$Treatment)
permutest(dispersion, pairwise=TRUE, permutations=999) 
#F=0.0325, df=1,51, p=0.851

BC_Data <- metaMDS(abundanceWide[,8:139])
sites <- 1:nrow(abundanceWide)
BC_Meta_Data <- abundanceWide[,1:7]
plot(BC_Data$points, col=as.factor(BC_Meta_Data$Treatment))
ordiellipse(BC_Data, groups = as.factor(BC_Meta_Data$Treatment), kind = "sd", display = "sites", label = T)

#Generate ellipses
BC_NMDS = data.frame(MDS1 = BC_Data$points[,1], MDS2 = BC_Data$points[,2], group=BC_Meta_Data$Treatment)
BC_NMDS_Graph <- cbind(BC_Meta_Data, BC_NMDS)
BC_Ord_Ellipses<-ordiellipse(BC_Data, BC_Meta_Data$Treatment, display = "sites",
                             kind = "se", conf = 0.95, label = T)

BC_Ellipses <- data.frame()
#Generate ellipses points
for(g in levels(as.factor(BC_NMDS$group))){
  BC_Ellipses <- rbind(BC_Ellipses,
                       cbind(as.data.frame(with(BC_NMDS[BC_NMDS$group==g,], 
                                                veganCovEllipse(BC_Ord_Ellipses[[g]]$cov,
                                                                BC_Ord_Ellipses[[g]]$center,
                                                                BC_Ord_Ellipses[[g]]$scale))),
                             group=g))
}

#Plot the data from BC_NMDS_Graph, where x=MDS1 and y=MDS2, make an ellipse based on "group"
NMDS_ABG_VS_PBG <- ggplot(BC_NMDS_Graph, aes(x=MDS1, y=MDS2, color=group,linetype = group, shape = group)) +
  geom_point(size=6) + 
  geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color="", linetype = "", shape = "") +
  scale_colour_manual(values=c("blue", "red"), name = "") +
  scale_linetype_manual(values = c("solid", "twodash"), name = "") +
  scale_shape_manual(values = c(19, 19)) +
  xlab("NMDS1") + 
  ylab("NMDS2") + 
  theme_bw() +
  theme(axis.text.x=element_text(size=24, color = "black"), 
        axis.text.y = element_text(size = 24, color = "black"), 
        axis.title.x = element_text(size = 24, color = 'black'),
        axis.title.y = element_text(size = 24, color = 'black'),
        legend.text = element_text(size = 24),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  ) +
  annotate("text", x=min(BC_NMDS_Graph$MDS1), y=max(BC_NMDS_Graph$MDS2),
           label="F=1.50, p=0.122", size=6, hjust=0, vjust=1)


#F=1.495, df=1,51, p=0.122





#### Bootstrapping! ####

set.seed(123)

Joined_New <- Joined %>% filter(Trea == "PBG") %>% group_by(block) %>% 
  mutate(plot_index=1:length(block))

Total_counts_New <- total_counts %>% filter(Treatment == "PBG") %>% group_by(Block) %>% 
  mutate(plot_index=1:length(Block))

Combined_Data <- left_join(Joined_New, Total_counts_New, by = c("WS", "Trans", "plot_index"))

num_bootstrap <- 1000
bootstrap_vector <- 1:num_bootstrap
PBG_plot_master <- data.frame()  # Initialize an empty dataframe

for (BOOT in bootstrap_vector) {
    Joined_New_Key <- Combined_Data %>%
    dplyr::select(1:13) %>%
    unique() %>%
    group_by(block) %>%
    sample_n(16, replace = TRUE) %>%
    dplyr::select(plot_index) %>%
    ungroup()
  
  # Join the sampled rows back to the original dataframe
  PBG_plot_ready <- Combined_Data %>%
    right_join(Joined_New_Key, by = c("block", "plot_index")) %>%
    mutate(iteration = BOOT)
  
  # Append the results to the master dataframe
  PBG_plot_master <- rbind(PBG_plot_master, PBG_plot_ready)
}

#### Z-Score Calculations ####

#Getting average richness per iteration for bootstrapped dataframe
average_richness <- PBG_plot_master %>%
  group_by(iteration) %>%
  summarize(mean_richness = mean(richness))

#Getting evenness richness per iteration for bootstrapped dataframe
average_evenness <- PBG_plot_master %>%
  group_by(iteration) %>%
  summarize(mean_evenness = mean(Evar, na.rm = TRUE))

#Getting average richness per iteration for bootstrapped dataframe
average_total_count <- PBG_plot_master %>%
  group_by(iteration) %>%
  summarize(mean_count = mean(Count, na.rm = TRUE))

#Take the mean of the mean for PBG richness
PBG_mean_mean_richness <- mean(average_richness$mean_richness, na.rm = TRUE)

#Take the mean of the mean for PBG evenness
PBG_mean_mean_evenness <- mean(average_evenness$mean_evenness, na.rm = TRUE)

#Take the mean of the mean for PBG evenness
PBG_mean_mean_total_count <- mean(average_total_count$mean_count, na.rm = TRUE)

#Take the mean of the mean for PBG total count

#Getting ABG ready

Total_counts_ABG <- total_counts %>% filter(Treatment == "ABG") %>% group_by(Block)
  
Joined_New_ABG <- Joined %>% filter(Trea == "ABG") %>% group_by(block) %>% 
  mutate(plot_index=1:length(block))

Joined_New_ABG <- rbind(Joined_New_ABG, Total_counts_ABG)


#Take the mean of the mean for ABG richness
ABG_mean_richness <- mean(Joined_New_ABG$richness, na.rm = TRUE)

#Take the mean of the mean for ABG evenness
ABG_mean_evenness <- mean(Joined_New_ABG$Evar, na.rm = TRUE)

#Take the mean of the mean for ABG total count
ABG_mean_total_count <- mean(Joined_New_ABG$Count, na.rm = TRUE)

# Z-Score for richness 

Z_R <- ((ABG_mean_richness) - (PBG_mean_mean_richness))/(sd(average_richness$mean_richness))
Z_R

p_value_R <- 1 - pnorm(Z_R)

print(p_value_R)

#Z = 1.952155
#P = 0.02545992

# Z-Score for evenness 
Z_E <- ((ABG_mean_evenness) - (PBG_mean_mean_evenness))/(sd(average_evenness$mean_evenness))
Z_E

p_value_E <- 1 - pnorm(Z_E)

print(p_value_E)
#1.280609
#P = 0.1001654

# Z-Score for total count 
Z_C <- ((ABG_mean_total_count) - (PBG_mean_mean_total_count))/(sd(average_total_count$mean_count))
Z_C

p_value_C <- 1 - pnorm(Z_C, lower.tail = FALSE)

print(p_value_C)
#Z-Score: -1.125048
#P = 0.1302844



#Some object with bootstrap PBG values (richness, evennness, count)
#Compare to ABG mean



#### Graphs Bootstrapped Means ####

#Richness means graph
avg_richness <- ggplot(average_richness, aes(x = mean_richness, color = "PBG")) +
  geom_density(aes(y = ..scaled..), alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_richness, color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "Density Plot of Mean Richness",
       x = "Mean Richness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
  annotate("text", x = min(average_richness$mean_richness), y = 1, 
           label = "z-value = 1.95, p = 0.0255", hjust = 0, vjust = 1, size = 4, color = "black") + 
  theme_bw() +
  theme( panel.grid.major=element_blank(), panel.grid.minor=element_blank())

#Z = 1.952155
#P = 0.02545992


#Evenness means graph


avg_evenness <- ggplot(average_evenness, aes(x = mean_evenness, color = "PBG")) +
  geom_density(aes(y = ..scaled..), alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_evenness, color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "Density Plot of Mean Evenness",
       x = "Mean Evenness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
  annotate("text", x = min(average_evenness$mean_evenness), y = 1, 
           label = "z-value = 1.28, p = 0.100", hjust = 0, vjust = 1, size = 4, color = "black") + 
  theme_bw() +
  theme( panel.grid.major=element_blank(), panel.grid.minor=element_blank())



#1.280609
#P = 0.1001654


#Total_count means graphs

ggplot(average_total_count, aes(x = mean_count, color = "PBG")) +
  geom_density(aes(y = ..scaled..), alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_total_count, color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "Density Plot of Abundance",
       x = "Mean Abundance",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
  annotate("text", x = min(average_total_count$mean_count), y = 1, 
           label = "z-value = -1.13, p = 0.130", hjust = 0, vjust = 1, size = 4, color = "black") +
  theme_bw() +
  theme( panel.grid.major=element_blank(), panel.grid.minor=element_blank())



#Z-Score: -1.125048
#P = 0.1302844


#### NMDS Two Panel ####

arranged_plots <- grid.arrange(NMDS_ABG_VS_PBG, NMDS_Years_Since_Burned, ncol = 2)


ggsave("arranged_plots.png", arranged_plots, width = 25, height = 5)

### Richness and evenness two panel ####

evenness_richness <- grid.arrange(avg_evenness, avg_richness, ncol = 2)
ggsave("belowground_evenness_richness.png", evenness_richness, width = 20, height = 10)

