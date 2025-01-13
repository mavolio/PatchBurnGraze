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
library(writexl)


#### Pairwise Function ####
pairwise.adonis <- function(x, factors, sim.method = "bray", p.adjust.m = "bonferroni") {
  library(vegan)
  co = combn(unique(factors), 2)
  pairs = list()
  for (elem in 1:ncol(co)) {
    factors_sub = factors[factors %in% co[, elem]]
    x_sub = x[factors %in% co[, elem], ]
    ad = adonis(x_sub ~ factors_sub, method = sim.method)
    pairs[[elem]] = c(co[1, elem], co[2, elem], ad$aov.tab[1, 6])
  }
  pairs = as.data.frame(do.call(rbind, pairs))
  pairs$p.adjusted = p.adjust(pairs$V3, method = p.adjust.m)
  colnames(pairs) = c("Group1", "Group2", "p.value", "p.adjusted")
  return(pairs)
}

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

#### Unique Morpsospp ####

### Unique in each ###

Abundance_South <- Abundance_ID_Belowground_Reformated %>% 
mutate(Block = ifelse(grepl("S", WS), "North", "South"))


unique_overall <- Abundance_South %>% 
  distinct(Morphospp)

write_xlsx(unique_overall, "unique_species.xlsx")

# N = 133


# Filter for ABG and get unique morphospecies
unique_ABG <- Abundance_South %>%
  filter(Treatment == "ABG") %>%
  distinct(Morphospp)

write_xlsx(unique_ABG, "ABG_species.xlsx")


# 52 unique

# Filter for PBG and get unique morphospecies
unique_PBG <- Abundance_South %>%
  filter(Treatment == "PBG") %>%
  distinct(Morphospp)
# 104 unique

write_xlsx(unique_PBG, "PBG_species.xlsx")

# Find morphospecies that are unique to ABG
only_ABG <- setdiff(unique_ABG$Morphospp, unique_PBG$Morphospp)
only_ABG
# 28 ABG has that PBG doesn't

# Find morphospecies that are unique to PBG
only_PBG <- setdiff(unique_PBG$Morphospp, unique_ABG$Morphospp)
only_PBG
# 80 PBG has that ABG doesn't


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
  summarize(Count = sum(Count, na.rm = TRUE)) %>%
  filter(Sample != 'C1A_A_38_ABG')


Abundance_Stats$Treatment <- ifelse(grepl("ABG", Abundance_Stats$Sample), "ABG", 
                                    ifelse(grepl("PBG", Abundance_Stats$Sample), "PBG", NA))

#Adding Block for use in modeling
total_counts <- aggregate(Count ~ Sample + Treatment, data = Abundance_Stats, FUN = sum) %>% 
  separate(Sample, into = c("WS", "Trans", "Dist"), sep = "_", extra = "drop") %>%
  mutate(Block = ifelse(grepl("S", WS), "North", "South"))

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

# numDF denDF  F-value p-value
# (Intercept)     1    49 227.8829  <.0001
# TreatmentSB     3    49   0.7795  0.5111

#### Prep for Richness + Evenness Graphs & Stats ####

#Getting  community data and tagging ABG and PBG
commMetrics <- Abundance_Stats %>%  community_structure(abundance.var='Count', replicate.var='Sample')
commMetrics$Treatment <- ifelse(grepl("C1", commMetrics$Sample), "ABG", "PBG") 

#Adding Blocks
commMetrics2 <- commMetrics %>% mutate(block = ifelse(grepl("S", Sample), "North", "South")) %>% 
  filter(Sample != 'C1A_A_38_ABG') %>% 
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

# numDF denDF  F-value p-value
# (Intercept)     1    49 380.4980  <.0001
# TreatmentSB     3    49   1.9081  0.1406

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

# numDF denDF  F-value p-value
# (Intercept)     1    46 326.5477  <.0001
# TreatmentSB     3    46   0.7888  0.5064

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

# grid.arrange(richness, evenness, counts, Weight) #Weight comes from weight script!

#### multivariate community response - PERMANOVA and NMDS ####

### by Years since burn

#Combing count with burn info

CountGraph <- full_join(Abundance_ID_Belowground_Reformated, BurnInfo) %>% 
  unite("TreatmentSB",c("Treatment","SB"), sep="_")

CountGraph$Treatment <- ifelse(grepl("C1", CountGraph$Sample), "ABG", "PBG") 

#New Code to drop missing unit (block)
CountGraph_New <- CountGraph %>% 
  mutate(Block = ifelse(grepl("S", WS), "North", "South"))

CountGraph_Filtered <- CountGraph_New %>% 
  filter(Block == "South")


# PERMANOVA

set.seed(124)


abundanceWide <- CountGraph_Filtered %>% 
  mutate(block = ifelse(grepl("S", Sample), "North", "South")) %>%
  select(Sample, block, WS, Trans, Dist., Treatment, TreatmentSB, Morphospp, Count) %>%
  group_by(Sample, block, WS, Trans, Dist., Treatment, TreatmentSB, Morphospp) %>% 
  summarise(Count = sum(Count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = 'Morphospp', values_from = 'Count', values_fill = list(Count = 0))

abundanceWide <- abundanceWide %>% 
  mutate(sum = rowSums(abundanceWide[, c(8:86)], na.rm = TRUE)) %>% 
  filter(sum > 0, !(Sample %in% c('C1A_A_38_ABG', 'C3SA_A_16_PBG', 'C3SA_C_38_PBG'))) %>%
  select(-sum)  # Drop the 'sum' column


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

print(permanova <- adonis2(formula = abundanceWide[,8:86]~TreatmentSB, data=abundanceWide, permutations=999, method="bray"))
#F=1.2401 , df=3,27, p=0.187

pairwise_results <- pairwise.adonis(abundanceWide[, 8:86], abundanceWide$TreatmentSB)
print(pairwise_results)


#betadisper
veg <- vegdist(abundanceWide[,8:86], method = "bray")
dispersion <- betadisper(veg, abundanceWide$TreatmentSB)
permutest(dispersion, pairwise=TRUE, permutations=999) 
#F=0.8866 , df=3,27, p=0.467

BC_Data <- metaMDS(abundanceWide[,8:86])
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
  theme(axis.text.x=element_text(size=34, color = "black"), 
        axis.text.y = element_text(size = 34, color = "black"), 
        axis.title.x = element_text(size = 34, color = 'black'),
        axis.title.y = element_text(size = 34, color = 'black'),
        legend.text = element_text(size = 34),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )  
  # annotate('text', x = min(BC_NMDS_Graph$MDS1) - 0.04, y = min(BC_NMDS_Graph$MDS2) + 0.035, 
  #          label = 'Mean p = 0.042\nVariance p = 0.530', size = 10, hjust = 'left')

NMDS_Years_Since_Burned
#Mean F=1.5242, df=3,49, p=0.042 
#Variance #F=0.773 , df=3,49, p=0.53

#export at 1500x1000

### by watershed
#Subsampling
set.seed(322)

ABG_Test <- abundanceWide %>% 
  filter(Treatment == "ABG")

#Filter PBG
PBG_Test <- abundanceWide %>% 
  filter(Treatment == "PBG")

# Set seed for reproducibility

# Get unique samples
unique_samples <- unique(PBG_Test$Sample)

# Randomly select 16 unique samples
subsamples <- sample(unique_samples, 16, replace = FALSE)

# Filter the data frame based on the selected samples
subsampled_data <- PBG_Test %>% filter(Sample %in% subsamples)

#New Abundance_Data2021 with subsamples
Abundance_Data <- full_join(subsampled_data, ABG_Test)


# PERMANOVA
print(permanova <- adonis2(formula = Abundance_Data[,8:86]~Treatment, data=Abundance_Data, permutations=999, method="bray"))
#F=1.5169 , df=1,21, p=0.146

#betadisper
veg <- vegdist(Abundance_Data[,8:86], method = "bray")
dispersion <- betadisper(veg, Abundance_Data$Treatment)
permutest(dispersion, pairwise=TRUE, permutations=999) 
#F=0.8234, df=1,30, p=0.369

BC_Data <- metaMDS(Abundance_Data[,8:86])
sites <- 1:nrow(Abundance_Data)
BC_Meta_Data <- Abundance_Data[,1:7]
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
NMDS_ABG_VS_PBG_4 <- ggplot(BC_NMDS_Graph, aes(x=MDS1, y=MDS2, color=group, linetype = group, shape = group)) +
  geom_point(size=6) + 
  geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color="", linetype = "", shape = "") +
  scale_colour_manual(values=c("blue", "red"), name = "") +
  scale_linetype_manual(values = c("solid", "twodash"), name = "") +
  scale_shape_manual(values = c(19, 19)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  xlab("NMDS1") + 
  ylab("NMDS2") + 
  theme_bw() +
  theme(axis.text.x=element_text(size=34, color = "black"), 
        axis.text.y = element_text(size = 34, color = "black"), 
        axis.title.x = element_text(size = 34, color = 'black'),
        axis.title.y = element_text(size = 34, color = 'black'),
        legend.text = element_text(size = 34),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )

NMDS_ABG_VS_PBG_4
#  annotate('text', x = min(BC_NMDS_Graph$MDS1) - 0.04, y = min(BC_NMDS_Graph$MDS2) + 0.2, 
#           label = 'Mean p = 0.206\nVariance p = 0.847', size = 10, hjust = 'left')


grid.arrange(NMDS_ABG_VS_PBG_1, NMDS_ABG_VS_PBG_2, NMDS_ABG_VS_PBG_3, NMDS_ABG_VS_PBG_4, ncol = 2, nrow = 2)

#NMDS_ABG_VS_PBG_1
#Mean F1-21 = 1.2598,p =0.251 \nVariance F1-21= 0.2716, p = 0.621
#setseed(312)

#NMDS_ABG_VS_PBG_2
#Mean F121 = 1.5035, p = 0.126\nVariance f1-21 = 0.5147, p = 0.492
#set.seed(124)

#NMDS_ABG_VS_PBG_3
#Mean F1-21 = 1.2668, p = 0.211\nVariance F1-21 =0.5754, p = 0.458
#set.seed(321)

#NMDS_ABG_VS_PBG_4
#Mean F1-21 = 1.7229,  p = 0.072 \nVariance F1-21 = 0.4678,  p = 0.500
#set.seed(322)

#### Simper Analysis ####

write_xlsx(Abundance_Data, "For_SIMPER_ABG_VS_PBG.xlsx")


# Perform SIMPER analysis
simper_results <- simper(Abundance_Data[, 8:86], group = Abundance_Data$TreatmentSB, permutations = 999)

# Print SIMPER results
print(simper_results)

# To view the contribution of each species, you can examine the output in detail
summary(simper_results)

# Perform SIMPER analysis
simper_results <- simper(Abundance_Data[, 8:86], group = Abundance_Data$Treatment, permutations = 999)

# Print SIMPER results
print(simper_results)

# To view the contribution of each species, you can examine the output in detail
summary(simper_results)


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



#### Bootstrapped Graphs Settings ####

#### Graphs Bootstrapped Means ####

#Richness means graph
avg_richness <- ggplot(average_richness, aes(x = mean_richness, color = "PBG")) +
  geom_density(aes(y = ..scaled..), alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_richness, color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "Density Plot of Mean Richness",
       x = "Mean Richness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
#  annotate("text", x = min(average_richness$mean_richness), y = 1, 
#           label = "z-value = 1.95, p = 0.0255", hjust = 0, vjust = 1, size = 4, color = "black") + 
  theme_bw() +
  theme( panel.grid.major=element_blank(), 
         panel.grid.minor=element_blank(), 
         legend.position=c(0.2,0.82), 
         axis.text = element_text(size = 25),
         legend.text = element_text(size = 29),
         axis.title = element_text(size = 25),
         axis.text.y = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         axis.ticks.y = element_line(size = 1))

avg_richness
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
#  annotate("text", x = min(average_evenness$mean_evenness), y = 1, 
#           label = "z-value = 1.28, p = 0.100", hjust = 0, vjust = 1, size = 4, color = "black") + 
  theme_bw() +
  theme( panel.grid.major=element_blank(), 
         panel.grid.minor=element_blank(), 
         legend.position= "none", 
         axis.text = element_text(size = 25),
         axis.title = element_text(size = 25),
         axis.text.y = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         axis.ticks.y = element_line(size = 1))

avg_evenness
# theme(legend.position=c(0.15,0.7))


#1.280609
#P = 0.1001654


#Total_count means graphs

total_count <- ggplot(average_total_count, aes(x = mean_count, color = "PBG")) +
  geom_density(aes(y = ..scaled..), alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_total_count, color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "Density Plot of Abundance",
       x = "Mean Abundance",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
 # annotate("text", x = min(average_total_count$mean_count), y = 1, 
#           label = "z-value = -1.13, p = 0.130", hjust = 0, vjust = 1, size = 4, color = "black") +
  theme_bw() +
  theme( panel.grid.major=element_blank(), 
         panel.grid.minor=element_blank(), 
         legend.position="none", 
         axis.text = element_text(size = 25),
         axis.title = element_text(size = 30),
         axis.text.y = element_text(size = 25),
         axis.title.y = element_text(size = 25),
         axis.ticks.y = element_line(size = 1))


#Z-Score: -1.125048
#P = 0.1302844

#### Three Panel univariate graph ####

evenness_richness <- grid.arrange(avg_richness, avg_evenness, total_count, ncol = 3)

evenness_richness

#### NMDS Two Panel ####

arranged_plots <- grid.arrange(NMDS_ABG_VS_PBG, NMDS_Years_Since_Burned, ncol = 2)


ggsave("arranged_plots.png", arranged_plots, width = 25, height = 5)

#### Richness and evenness two panel ####

evenness_richness <- grid.arrange(avg_evenness, avg_richness, ncol = 2)
ggsave("belowground_evenness_richness.png", evenness_richness, width = 20, height = 10)


#### Years Since Burned Graph Settings ####
Axis_Label_Size_1 <- 20 #The actualy ticks and treatments
Axis_Text_Size <- 20 #Text on axislike abuandance, richness etc
#### Richness, Count and Evenness Years Since Burn ####
# Boxplot for Richness
richness_below <- ggplot(Joined, aes(x=TreatmentSB, y=richness, fill=TreatmentSB)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) + # Add geom_jitter
  labs(title="", x="Years Since Burned", y="Richness") +
  scale_fill_manual(values = c("ABG_0" = "blue", "PBG_0" = "red", "PBG_1" = "red", "PBG_2" = "red")) + # Color PBG_0, PBG_1, PBG_2 as red, ABG_0 as blue
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(size = 9),
    axis.text.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis labels
    axis.text.y = element_text(size = Axis_Label_Size_1),
    axis.title.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis title
    axis.title.y = element_text(size = Axis_Label_Size_1), # Increase the size of y-axis title
    panel.grid.major = element_blank(), # Remove major gridlines
    panel.grid.minor = element_blank()  # Remove minor gridlines
  ) +
  scale_x_discrete(labels = c("ABG_0" = "ABG 0", "PBG_0" = "PBG 0", "PBG_1" = "PBG 1", "PBG_2" = "PBG 2")) + # Set custom axis labels
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) # Adjust x-axis text angle

# Display the plot
richness_below

# numDF denDF   F-value p-value
# (Intercept)     1    50 178.46582  <.0001
# TreatmentSB     3    50   0.88895  0.4533

Evar_below <- ggplot(Joined, aes(x=TreatmentSB, y=Evar, fill=TreatmentSB)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) + # Add geom_jitter
  labs(title="", x="Years Since Burned", y="Evenness") +
  scale_fill_manual(values = c("ABG_0" = "blue", "PBG_0" = "red", "PBG_1" = "red", "PBG_2" = "red")) + # Color PBG_0, PBG_1, PBG_2 as red, ABG_0 as blue
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(size = 9),
    axis.text.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis labels
    axis.text.y = element_text(size = Axis_Label_Size_1),
    axis.title.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis title
    axis.title.y = element_text(size = Axis_Label_Size_1), # Increase the size of y-axis title
    panel.grid.major = element_blank(), # Remove major gridlines
    panel.grid.minor = element_blank()  # Remove minor gridlines
  ) +
  scale_x_discrete(labels = c("ABG_0" = "ABG 0", "PBG_0" = "PBG 0", "PBG_1" = "PBG 1", "PBG_2" = "PBG 2")) + # Set custom axis labels
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) # Adjust x-axis text angle

# Display the plot
Evar_below

# numDF denDF   F-value p-value
# (Intercept)     1    47 219.24688  <.0001
# TreatmentSB     3    47   0.51725  0.6724

# Count graph here
Count_below <- ggplot(total_counts2, aes(x=TreatmentSB, y=Count, fill=TreatmentSB)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) + # Add geom_jitter
  labs(title="", x="Years Since Burned", y="Abundance") +
  scale_fill_manual(values = c("ABG_0" = "blue", "PBG_0" = "red", "PBG_1" = "red", "PBG_2" = "red")) + # Color PBG_0, PBG_1, PBG_2 as red, ABG_0 as blue
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(size = 9),
    axis.text.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis labels
    axis.text.y = element_text(size = Axis_Label_Size_1),
    axis.title.x = element_text(size = Axis_Label_Size_1), # Increase the size of x-axis title
    axis.title.y = element_text(size = Axis_Label_Size_1), # Increase the size of y-axis title
    panel.grid.major = element_blank(), # Remove major gridlines
    panel.grid.minor = element_blank()  # Remove minor gridlines
  ) +
  scale_x_discrete(labels = c("ABG_0" = "ABG 0", "PBG_0" = "PBG 0", "PBG_1" = "PBG 1", "PBG_2" = "PBG 2")) + # Set custom axis labels
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) # Adjust x-axis text angle


Count_below

# numDF denDF   F-value p-value
# (Intercept)     1    50 231.04818  <.0001
# TreatmentSB     3    50   0.86496  0.4655

#### Legend Graph ####
# Count graph here
legend_below <- ggplot(total_counts2, aes(x=TreatmentSB, y=Count, fill=TreatmentSB)) +
  geom_boxplot() +
  labs(title="Richness Comparison for Treatments with Years Since Burned", x="Treatment", y="Richness") +
  scale_fill_manual(values = c("ABG_0" = "blue", "PBG_0" = "red", "PBG_1" = "red", "PBG_2" = "red")) + # Color PBG_0, PBG_1, PBG_2 as red, ABG_0 as blue
  theme_minimal() +
  annotate("text", x = -Inf, y = Inf, label = "P value: 0.8307", vjust = 1, hjust = 0, size = 3.5, color = "black") + 
  theme(plot.title = element_text(size = 9)) # Adjust size as needed

legend_below

#legend <- get_legend(legend_below)

#### Big Graph of Years Since Burned ####

multi_panel_graph <- grid.arrange(richness_below, Evar_below, Count_below, 
                                  nrow = 1 
                                  #, bottom = legend
                                  )


#### CV Richness ####
ABG_CV_mean_R <- Joined %>%
  filter(TreatmentSB == "ABG_0")  %>%
  summarize(CV_richness = sd(`richness`) / mean(`richness`))


PBG_average_cv_R <- PBG_plot_master %>%
  group_by(iteration) %>%
  summarize(CV_richness = sd(richness, na.rm = T) / mean(richness, na.rm = TRUE))


PBG_mean_mean_richness_CV <- mean(PBG_average_cv_R$CV_richness, na.rm = TRUE)

Z_R_CV <- ((ABG_CV_mean_R$CV_richness) - (PBG_mean_mean_richness_CV))/(sd(PBG_average_cv_R$CV_richness))
Z_R_CV

# -0.1332362

p_value_R_CV <- 2*pnorm(-abs(Z_R_CV))
p_value_R_CV

#  0.8940066

#### CV Richness graph ####

CV_Richness <- ggplot(PBG_average_cv_R, aes(x = CV_richness, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_CV_mean_R$CV_richness , color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "Density Plot of Mean CV Richness",
       x = "Mean CV Richness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
  # annotate("text", x = 20, y = 1, label = "p-value: 0.217", color = "blue", size = 4, hjust = 1, vjust = 1) +
  # annotate("text", x = 20, y = 0.80, label = "z-score:  1.24", color = "red", size = 4, hjust = 1, vjust = 1) +
  theme(legend.position=c(0.15,0.7),  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22)) 
# scale_x_continuous(limits = c(5, 20)) +
#  scale_y_continuous(limits = c(0, 1)) 

CV_Richness
#### CV Evenness ####

ABG_CV_mean_E <- Joined %>%
  filter(TreatmentSB == "ABG_0")  %>%
  summarize(CV_evenness = sd(`Evar`) / mean(`Evar`))
ABG_CV_mean_E

PBG_average_cv_E <- PBG_plot_master %>%
  group_by(iteration) %>%
  summarize(CV_evenness = sd(Evar, na.rm = T) / mean(Evar, na.rm = TRUE))
PBG_average_cv_E

PBG_mean_mean_evenness_CV <- mean(PBG_average_cv_E$CV_evenness, na.rm = TRUE)

Z_E_CV <- ((ABG_CV_mean_E$CV_evenness) - (PBG_mean_mean_evenness_CV))/(sd(PBG_average_cv_E$CV_evenness))
Z_E_CV

#-2.604795

p_value_E_CV <- 2*pnorm(-abs(Z_E_CV))
p_value_E_CV

#0.009192931

#### CV Evenness Graph ####
CV_Evar <- ggplot(PBG_average_cv_E, aes(x = CV_evenness, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_CV_mean_E$CV_evenness , color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "Density Plot of Mean CV Evenness",
       x = "Mean CV Evenness",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
  # annotate("text", x = 20, y = 1, label = "p-value: 0.217", color = "blue", size = 4, hjust = 1, vjust = 1) +
  # annotate("text", x = 20, y = 0.80, label = "z-score:  1.24", color = "red", size = 4, hjust = 1, vjust = 1) +
  theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22)) 
# scale_x_continuous(limits = c(5, 20)) +
# scale_y_continuous(limits = c(0, 1)) 

CV_Evar
#### CV Count ####
ABG_mean_CV_count_df <- Abundance_Stats %>%
  group_by(Treatment) %>% 
  filter(Treatment == "ABG") %>%
  summarize(CV_Count = sd(Count, na.rm = TRUE) / mean(Count, na.rm = TRUE))


# To extract the CV_count value as a single number
ABG_mean_CV_C <- ABG_mean_CV_count_df$CV_Count
ABG_mean_CV_C

PBG_average_cv_C <- PBG_plot_master %>%
  group_by(iteration) %>%
  summarize(CV_Count = sd(Count, na.rm = T) / mean(Count, na.rm = TRUE))
PBG_average_cv_C

PBG_mean_mean_count_CV <- mean(PBG_average_cv_C$CV_Count, na.rm = TRUE)
PBG_mean_mean_count_CV

Z_C_CV <- ((ABG_mean_CV_C) - (PBG_mean_mean_count_CV))/(sd(PBG_average_cv_C$CV_Count))
Z_C_CV

#4.758492

p_value_C_CV <- 2*pnorm(-abs(Z_C_CV))
p_value_C_CV

#1.950446e-06

#### CV Count Graph ####
CV_Count <- ggplot(PBG_average_cv_C, aes(x = CV_Count, y = ..scaled.., color = "PBG")) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = ABG_mean_CV_C , color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "Density Plot of Mean CV Abundance",
       x = "Mean CV Abundance",
       y = "Density") +
  scale_color_manual(values = c("blue", "red"), name = "Legend") +
  # annotate("text", x = 20, y = 1, label = "p-value: 0.217", color = "blue", size = 4, hjust = 1, vjust = 1) +
  # annotate("text", x = 20, y = 0.80, label = "z-score:  1.24", color = "red", size = 4, hjust = 1, vjust = 1) +
  theme(legend.position = "none",  plot.title = element_text(size = 15),  axis.text.x = element_text(size = 22),  # Adjust the size as needed
        axis.text.y = element_text(size = 22)) + 
  scale_x_continuous(limits = c(0.2, 1.6)) +
  scale_y_continuous(limits = c(0, 1)) 

CV_Count


#### Big CV Means Graph ####

grid.arrange(CV_Richness, CV_Evar, CV_Count, ncol = 3)


#### Beta Diversity ####
# Load necessary packages
set.seed(123)

library(vegan)
library(dplyr)

#Caculating Beta Diversity for PBG

Joined_New_1 <- abundanceWide %>%
  filter(Treatment == "PBG") %>%
  group_by(block) %>%
  mutate(plot_index = row_number()) # Corrected index generation

num_bootstrap <- 1000
PBG_plot_master <- data.frame(iteration = 1:num_bootstrap, beta_diversity = numeric(num_bootstrap))

for (BOOT in 1:num_bootstrap) {
  Joined_New_Key <- Joined_New_1 %>%
    sample_n(16, replace = TRUE) %>%
    select(plot_index) %>%
    ungroup()
  
  PBG_plot_ready <- Joined_New_1 %>%
    right_join(Joined_New_Key, by = "plot_index") %>%
    mutate(iteration = BOOT)
  
  data_matrix <- PBG_plot_ready %>%
    select(8:141) # Species column
  
  dissimilarity_matrix <- vegdist(data_matrix, method = "bray")
  
  beta_diversity_value <- mean(dissimilarity_matrix)
  
  PBG_plot_master$beta_diversity[BOOT] <- beta_diversity_value
}

#Add Treatment column

PBG_plot_master$Treatment <- 'PBG'


# Print the results
print(PBG_plot_master)

ggplot(PBG_plot_master, aes(x = beta_diversity)) +
  geom_density() +
  labs(title = "Density Plot of Beta Diversity",
       x = "Beta Diversity",
       y = "Density") +
  theme_bw()



#Calculating ABG betadiveristy

Joined_New_2 <- abundanceWide %>%
  filter(Treatment == "ABG")

data_matrix_2 <- Joined_New_2 %>%
  select(8:140) # Species column

dissimilarity_matrix_2 <- vegdist(data_matrix, method = "bray")

beta_diversity_value_2 <- mean(dissimilarity_matrix)

#### Beta Diversity Graph ####


# Beta Diversity Density plot
ggplot(PBG_plot_master, aes(x = beta_diversity)) +
  geom_density(aes(y = ..scaled.., color = "PBG"), alpha = 0.5) +
  geom_vline(aes(xintercept = beta_diversity_value_2, color = "ABG"), linetype = "dashed", size = 1) +
  labs(title = "Density Plot of Beta Diversity",
       x = "Mean Beta Diversity",
       y = "Density") +
  scale_color_manual(values = c("PBG" = "red", "ABG" = "blue"), name = "Legend") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    legend.position = c(0.1, 0.8), 
    legend.text = element_text(size = 30),  # Adjust legend text size
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 30),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 30),
    axis.ticks.y = element_line(size = 1)
  )



#### Beta Diversity Z-score ####

#Getting average richness per iteration for bootstrapped dataframe
PBG_average_beta <- PBG_plot_master %>%
  summarize(mean_beta = mean(beta_diversity))


# Z-Score for richness 

Z_B <- ((beta_diversity_value_2) - (PBG_average_beta$mean_beta))/(sd(PBG_plot_master$beta_diversity))
Z_B


p_value_B <- 1 - pnorm(Z_B)

#lower.tail = FALSE

print(p_value_B)

#Z-score = 0.336496, p = 0.3682484

#Do Simper analysis




#### Relative Abundance of Epigeic worms ####
# Filter for Annelida_Clitellata_Opisthopora_Epigeic
filtered_data <- CountGraph_Filtered %>%
  filter(Morphospp == "Annelida_Clitellata_Opisthopora_Epigeic")
  

summary_data <- filtered_data %>%
  group_by(TreatmentSB) %>%
  summarize(total_count = sum(Count),
            average_count = mean(Count))

# Calculate relative abundance 
total_counts <- sum(summary_data$total_count)
summary_data <- summary_data %>%
  mutate(relative_abundance = total_count / total_counts)

# Display the result
print(summary_data)

# Redoing for ABG vs PBG 

summary_data_1 <- filtered_data %>%
  group_by(Treatment) %>%
  summarize(total_count = sum(Count),
            average_count = mean(Count))

# Calculate relative abundance 
total_counts_1 <- sum(summary_data_1$total_count)
summary_data_1 <- summary_data_1 %>%
  mutate(relative_abundance = total_count / total_counts_1)

# Display the result
print(summary_data_1)

#

#mean for outliars
mean_count_opisthoptera <- CountGraph_New %>%
  filter(Sample %in% c('C1A_A_38_ABG', 'C3SA_A_16_PBG', 'C3SA_C_38_PBG')) %>%
  filter(Morphospp == "Annelida_Clitellata_Opisthopora_Endogeic") %>%
  summarise(mean_count = mean(Count, na.rm = TRUE))
mean_count_opisthoptera
#28
  
mean_count_opisthoptera_overall <- CountGraph_New %>%
  filter(Morphospp == "Annelida_Clitellata_Opisthopora_Endogeic") %>%
  summarise(mean_count = mean(Count, na.rm = TRUE))
mean_count_opisthoptera_overall
