#### Libararies ####

library(readxl)
library(codyn)
library(vegan)
library(tidyverse)
library(nlme)
library(emmeans)
library(gridExtra)


#### NOTE: need to correct these values for depth of sample.


#### Set Working Directories ####
setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\konza\\patch burn\\2021 Belowground Data') #kim's wd


#### CSV Read in and Clean up ####

Abundance_ID_Belowground <- read_excel("2021 Belowground ID + Abundance Data\\Abundance + ID Data.xlsx") %>% 
  rename(Life.Stage="Life Stage")


Abundance_ID_Belowground <- Abundance_ID_Belowground[, c("WS","Trans", "Dist.", "Phylum", "Class", "Order", "Morphospp", "Life.Stage","Count")]

#Adds a new column to tag "1" with ABG and "3" with PBG. 

Abundance_ID_Belowground$Treatment <- ifelse(grepl("1", Abundance_ID_Belowground$WS), "ABG", "PBG")



#Stores the PBG and ABG in separate variables

pbg_data <- subset(Abundance_ID_Belowground, Treatment == "PBG")
abg_data <- subset(Abundance_ID_Belowground, Treatment == "ABG")

#removes NA or white space.

new_pbg <- select(pbg_data, Morphospp, Count) %>% 
  filter(!is.na(Morphospp), !is.na(Count), Morphospp != "", Count != "")

#Fix typos in PBG

new_pbg$Morphospp <- str_to_title(new_pbg$Morphospp) #Capitalize first letter


new_pbg$Morphospp <- gsub("Earthworm Cacoon", "Earthworm Cocoon", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub(" Light Brown Diptera", "Light Brown Diptera", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Big Brown", "Big Brown Beetle", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Big Brown Adult", "Big Brown Beetle", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Big Brown Beetle Adult", "Big Brown Beetle", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Big Brown Beetle Beetle", "Big Brown Beetle", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Black Two Segment", "Black Two Segments", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Black Two Segmentss", "Black Two Segments", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Black Two Segmentss", "Black Two Segments", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Coffee Larvae Beetle", "Coffee Beetle Larvae", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Coffee Beetle", "Coffee Beetle Larvae", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Coffee Beetle Larvae Larvae", "Coffee Beetle Larvae", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Coffee Larvae Bettle", "Coffee Beetle Larvae", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Creme Colored", "Creme Larvae", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Light Brown", "Light Brown Diptera", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Light Brown Diptera Larvae", "Light Brown Diptera", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Light Brown Millipede ", "Light Brown Millipede", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Light Brown Diptera Milipede", "Light Brown Millipede", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Light Brown Diptera Millipede ", "Light Brown Millipede", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Light Brown Diptera Millipede", "Light Brown Millipede", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Mid Brown Stick Larvae", "Midbrown Stick Larvae", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Parasite Worm", "Parasitic Worm", new_pbg$Morphospp)
new_pbg$Morphospp <- gsub("Parasitic Worn", "Parasitic Worm", new_pbg$Morphospp)




# group the new dataframe by Morphospp and summarize Count
new_pbg <- group_by(new_pbg, Morphospp) %>% 
  summarize(Count = sum(Count))



View(new_pbg)

#ABG: Selects columns, removes whitespace

new_abg <- select(abg_data, Morphospp, Count) %>% 
  filter(!is.na(Morphospp), !is.na(Count), Morphospp != "", Count != "")

#ABG fix typos

new_abg$Morphospp <- str_to_title(new_abg$Morphospp) #Capitalize first letter

new_abg$Morphospp <- gsub("Earthworm Cacoon", "Earthworm Cocoon", new_abg$Morphospp)
new_abg$Morphospp <- gsub(" Light Brown Diptera", "Light Brown Diptera", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Big Brown", "Big Brown Beetle", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Big Brown Adult", "Big Brown Beetle", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Big Brown Beetle Adult", "Big Brown Beetle", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Big Brown Beetle Beetle", "Big Brown Beetle", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Black Two Segment", "Black Two Segments", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Black Two Segmentss", "Black Two Segments", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Black Two Segmentss", "Black Two Segments", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Coffee Larvae Beetle", "Coffee Beetle Larvae", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Coffee Beetle", "Coffee Beetle Larvae", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Coffee Beetle Larvae Larvae", "Coffee Beetle Larvae", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Coffee Larvae Bettle", "Coffee Beetle Larvae", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Creme Colored", "Creme Larvae", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Light Brown", "Light Brown Diptera", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Light Brown Diptera Larvae", "Light Brown Diptera", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Light Brown Millipede ", "Light Brown Millipede", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Light Brown Diptera Milipede", "Light Brown Millipede", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Light Brown Diptera Millipede ", "Light Brown Millipede", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Light Brown Diptera Millipede", "Light Brown Millipede", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Mid Brown Stick Larvae", "Midbrown Stick Larvae", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Parasitic Worn", "Parasitic Worm", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Creme Larvae ", "Creme Larvae", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Creme Larvae", "Creme Larvae", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Light Brown Diptera Diptera", "Light Brown Diptera Larvae", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Light Brown Diptera Diptera Diptera Larvae", "Light Brown Diptera Larvae", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Light Brown Diptera Diptera Shrimp Like", "Light Brown Diptera Larvae", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Light Brown Diptera Larvae Larvae", "Light Brown Diptera Larvae", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Light Brown Diptera Shrimp Like", "Light Brown Diptera Larvae", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Light Brown Diptera", "Light Brown Diptera Larvae", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Light Brown Diptera Larvae Larvae", "Light Brown Diptera Larvae", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Large Golden Shrimp", "Creme Larvae", new_abg$Morphospp)
new_abg$Morphospp <- gsub("Small Light Colored Back End", "Small Light Colored Black End", new_abg$Morphospp)


new_abg <- group_by(new_abg, Morphospp) %>% 
  summarize(Count = sum(Count))

View(new_abg)

#Cleanup for stats

Abundance_ID_Belowground$Morphospp <- gsub("Earthworm Cacoon", "Earthworm Cocoon", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub(" Light Brown Diptera", "Light Brown Diptera", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Big Brown", "Big Brown Beetle", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Big Brown Adult", "Big Brown Beetle", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Big Brown Beetle Adult", "Big Brown Beetle", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Big Brown Beetle Beetle", "Big Brown Beetle", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Black Two Segment", "Black Two Segments", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Black Two Segmentss", "Black Two Segments", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Black Two Segmentss", "Black Two Segments", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Coffee Larvae Beetle", "Coffee Beetle Larvae", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Coffee Beetle", "Coffee Beetle Larvae", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Coffee Beetle Larvae Larvae", "Coffee Beetle Larvae", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Coffee Larvae Bettle", "Coffee Beetle Larvae", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Creme Colored", "Creme Larvae", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Light Brown", "Light Brown Diptera", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Light Brown Diptera Larvae", "Light Brown Diptera", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Light Brown Millipede ", "Light Brown Millipede", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Light Brown Diptera Milipede", "Light Brown Millipede", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Light Brown Diptera Millipede ", "Light Brown Millipede", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Light Brown Diptera Millipede", "Light Brown Millipede", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Mid Brown Stick Larvae", "Midbrown Stick Larvae", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Parasitic Worn", "Parasitic Worm", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Creme Larvae ", "Creme Larvae", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Creme Larvae", "Creme Larvae", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Light Brown Diptera Diptera", "Light Brown Diptera Larvae", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Light Brown Diptera Diptera Diptera Larvae", "Light Brown Diptera Larvae", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Light Brown Diptera Diptera Shrimp Like", "Light Brown Diptera Larvae", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Light Brown Diptera Larvae Larvae", "Light Brown Diptera Larvae", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Light Brown Diptera Shrimp Like", "Light Brown Diptera Larvae", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Light Brown Diptera", "Light Brown Diptera Larvae", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Light Brown Diptera Larvae Larvae", "Light Brown Diptera Larvae", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Large Golden Shrimp", "Creme Larvae", Abundance_ID_Belowground$Morphospp)
Abundance_ID_Belowground$Morphospp <- gsub("Small Light Colored Back End", "Small Light Colored Black End", Abundance_ID_Belowground$Morphospp)

Abundance_ID_Belowground$WS <- gsub("C35A","C3SA", Abundance_ID_Belowground$WS) 
Abundance_ID_Belowground$WS <- gsub("C35B","C3SB", Abundance_ID_Belowground$WS)
Abundance_ID_Belowground$WS <- gsub("C3CB", "C3C", Abundance_ID_Belowground$WS)
Abundance_ID_Belowground$WS <- gsub("C3D", "C3B", Abundance_ID_Belowground$WS)
Abundance_ID_Belowground$WS <- gsub("CSC", "C3C", Abundance_ID_Belowground$WS)

View(Abundance_ID_Belowground)

Abundance_ID_Belowground <- Abundance_ID_Belowground %>% 
  mutate(Morphospp = paste(Phylum, Class, Order, Morphospp, sep = "_"))

Abundance_ID_Belowground <- Abundance_ID_Belowground %>% 
  unite("Sample", c("WS", "Trans", "Dist.", "Treatment"), sep = "_", remove = FALSE) %>%
   mutate(Sample = str_remove(Sample, "m"))

Abundance_Stats <- group_by(Abundance_ID_Belowground, Morphospp, Sample) %>% 
  summarize(Count = sum(Count, na.rm = TRUE))

Abundance_ID_Belowground$Count <- ifelse(is.na(Abundance_ID_Belowground$Count), 1, Abundance_ID_Belowground$Count)


###Total Count ####

Abundance_Stats$treatment <- ifelse(grepl("ABG", Abundance_Stats$Sample), "ABG", 
                                    ifelse(grepl("PBG", Abundance_Stats$Sample), "PBG", NA))
# Calculate the average counts for each treatment

total_counts <- aggregate(Count ~ Sample + treatment, data = Abundance_Stats, FUN = sum) %>% 
  separate(Sample, into = c("WS", "Trans", "Dist"), sep = "_", extra = "drop") %>% mutate(block = ifelse(grepl("S", WS), "North", "South"))

avg_counts <- aggregate(Count ~ treatment, data = total_counts, FUN = mean)

TotalcountModel <- #stores the model output into a named list
  lme(Count ~ treatment, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = total_counts, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(TotalcountModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova
emmeans(TotalcountModel, pairwise~treatment, adjust="tukey") #this gives you contrasts (means and confidence intervals) for each possible pairwise comparison of treatments to know whether they are different or not (overlapping confidence intervals means not different)

#error bar


error_df <- total_counts %>%
  group_by(treatment) %>%
  summarize(mean_value = mean(Count),
            se = sd(Count) / sqrt(n()))


#Count Graph


library(ggplot2)

Avg_Count <- ggplot(avg_counts, aes(x = treatment, y = Count, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Treatment", y = "Average Count") +
  geom_errorbar(aes(ymin = Count - ifelse(treatment == "ABG", 1.726931, 1.418959),
                    ymax = Count + ifelse(treatment == "ABG", 1.726931, 1.418959)),
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








ggsave("Total_Count.png", width = 8, height = 8, dpi = 300)



#### Stats ####



commMetrics <- community_structure(Abundance_Stats, abundance.var='Count', replicate.var='Sample')

View(commMetrics)

# Tagging PBG and ABG again

commMetrics$Treatment <- ifelse(grepl("C1", commMetrics$Sample), "ABG", "PBG") 

#Adding North and South block

commMetrics2 <- commMetrics %>% mutate(block = ifelse(grepl("S", Sample), "North", "South")) %>% 
  separate(col = "Sample", into = c("WS","Trans","Dist.","Trea"), sep="_")

BurnInfo <- read_excel("YearsSenseBurned.xlsx")



Joined <- full_join(BurnInfo, commMetrics2) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_")

Joined$Trea <- ifelse(grepl("A", Joined$TreatmentSB), "ABG", "PBG") 

#Count with Burn data frame ####

CountGraph <- full_join(Abundance_ID_Belowground, BurnInfo) %>% 
  unite("TreatmentSB",c("Treatment","SB"), sep="_")

CountGraph$Treatment <- ifelse(grepl("C1", CountGraph$Sample), "ABG", "PBG") 




##### mixed model code #####

library(nlme) #library for mixed models
library(emmeans) #library for means comparisions after running the models

richModel <- #stores the model output into a named list
  lme(richness ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = Joined, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(richModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova
emmeans(richModel, pairwise~TreatmentSB, adjust="tukey") #this gives you contrasts (means and confidence intervals) for each possible pairwise comparison of treatments to know whether they are different or not (overlapping confidence intervals means not different)

EvenModel <- #stores the model output into a named list
  lme(Evar ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = na.omit(Joined), #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(EvenModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova
emmeans(EvenModel, pairwise~TreatmentSB, adjust="tukey") #this gives you contrasts (means and confid

#Graph Dataframe ####

Graphs_Stats <- Abundance_ID_Belowground %>% group_by(Sample) %>% 
  summarise(TotalCount = sum(Count)) %>%
  ungroup() %>% 
  separate("Sample",c("WS","Tras","Dist.","Treatment"), sep="_",remove = FALSE) %>% 
  mutate(block = ifelse(grepl("S", Sample), "North", "South"))

countModel <- #stores the model output into a named list
  lme(TotalCount ~ Treatment, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = na.omit(Graphs_Stats), #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(countModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova
emmeans(countModel, pairwise~Treatment, adjust="tukey") #this gives you contrasts (means and confidence intervals) for each possible pairwise comparison of treatments to know whether they are different or not (overlapping confidence intervals means not different)



#### Graphs ####

#Total Count

ggplot(CountGraph, aes(x = TreatmentSB, y = Count, fill= Treatment)) +
  geom_boxplot() +
  labs(
    x = "Treatment",
    y = "Total Count") +
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

ggsave("TotalCount.png", width = 8, height = 8, dpi = 300)

#richness graph

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
  

ggsave("Richness.png", width = 8, height = 8, dpi = 300)

#Evenness

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


ggsave("Evenness.png", width = 8, height = 8, dpi = 300)


#Caculating CV Count####

ABG_CV <- Abundance_ID_Belowground %>%
  filter(WS %in% c("C1A", "C1SB")) %>%
  group_by(WS) %>%
  summarize(cv_count = sd(`Count`) / mean(`Count`) * 100)

# Calculate the CV of weights for each combined PBG cluster
PBG_CV <- Abundance_ID_Belowground%>%
  filter(WS %in% c("C3A", "C3B", "C3C", "C3SA", "C3SB")) %>%
  group_by(WS) %>%
  summarize(cv_count = sd(`Count`) / mean(`Count`) * 100)

# Add the average values to the dataframes
ABG_CV2 <- ABG_CV %>% 
  mutate(`Treatment` = "ABG") %>% 
  select(-WS)

PBG_CV2 <- PBG_CV %>% 
  mutate(`Treatment` = "PBG")




ABG_CV_Average <- ABG_CV2 %>%
  group_by(Treatment) %>%
  summarize(Average_CV = mean(cv_count))

# Calculate the average weight CV for PBG
PBG_CV_Average <- PBG_CV2 %>%
  group_by(Treatment) %>%
  summarize(Average_CV = mean(cv_count))



combined_cv_count <- rbind(ABG_CV_Average, PBG_CV_Average)

# Graphs for CV ###

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

ggsave("CV_Count.png", width = 8, height = 8, dpi = 300)

#Calculating CV Richness####

ABG_CV <- Joined %>%
  filter(WS %in% c("C1A", "C1SB")) %>%
  group_by(WS) %>%
  summarize(cv_richness = sd(`richness`) / mean(`richness`) * 100)

# Calculate the CV of weights for each combined PBG cluster
PBG_CV <- Joined%>%
  filter(WS %in% c("C3A", "C3B", "C3C", "C3SA", "C3SB")) %>%
  group_by(WS) %>%
  summarize(cv_richness = sd(`richness`) / mean(`richness`) * 100)

# Add the average values to the dataframes
ABG_CV2 <- ABG_CV %>% 
  mutate(`Treatment` = "ABG") %>% 
  select(-WS)

PBG_CV2 <- PBG_CV %>% 
  mutate(`Treatment` = "PBG")




ABG_CV_Average <- ABG_CV2 %>%
  group_by(Treatment) %>%
  summarize(Average_CV = mean(cv_richness))

# Calculate the average weight CV for PBG
PBG_CV_Average <- PBG_CV2 %>%
  group_by(Treatment) %>%
  summarize(Average_CV = mean(cv_richness))



combined_cv_richness <- rbind(ABG_CV_Average, PBG_CV_Average)

# Graphs for CV ###

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

ggsave("Richness_CV.png", width = 8, height = 8, dpi = 300)




#### multivariate community response - PERMANOVA and NMDS ####

### by Years since burn
# PERMANOVA
abundanceWide <- CountGraph  %>% 
  mutate(block = ifelse(grepl("S", Sample), "North", "South")) %>%
  select(Sample, block, WS, Trans, Dist., Treatment, TreatmentSB, Morphospp, Count) %>%
  group_by(Sample, block, WS, Trans, Dist., Treatment, TreatmentSB, Morphospp) %>% 
  summarise(Count=sum(Count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from='Morphospp', values_from='Count', values_fill=list(Count=0)) %>% 
  mutate(sum=rowSums(abundanceWide[,c(8:57)], na.rm=TRUE)) %>% 
  filter(sum>0, Sample!='C1A_A_38_ABG') #PROBLEM: Check why two of the samples have nothing in them, is this real?
###IMPORTANT: removing sample from C1A_A_38_ABG, which is a big outlier because of high values of endogenic worms and brown shrimp-like beetles, which no other plots had

print(permanova <- adonis2(formula = abundanceWide[,8:166]~TreatmentSB, data=abundanceWide, permutations=999, method="bray"))
#F=1.6053, df=3,52, p=0.016

#betadisper
veg <- vegdist(abundanceWide[,8:166], method = "bray")
dispersion <- betadisper(veg, abundanceWide$TreatmentSB)
permutest(dispersion, pairwise=TRUE, permutations=999) 
#F=1.5519, df=3,52, p=0.213

BC_Data <- metaMDS(abundanceWide[,8:166])
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
ggplot(BC_NMDS_Graph, aes(x=MDS1, y=MDS2, color=group,linetype = group, shape = group)) +
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
  annotate('text', x=-0.03, y=0.02, label=expression(paste('F'['3,52'],' = 1.61')), size=8, hjust='left') +
  annotate('text', x=-0.03, y=0.018, label='p = 0.016', size=8, hjust='left') +
  theme_bw() +
  theme(axis.text.x=element_text(size=24, color = "black"), 
        axis.text.y = element_text(size = 24, color = "black"), 
        axis.title.x = element_text(size = 24, color = 'black'),
        axis.title.y = element_text(size = 24, color = 'black'),
        legend.text = element_text(size = 24),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
        )
#F=1.6053, df=3,52, p=0.016
#export at 1500x1000
patchBurn_belowgroundInvert_NMDS_yearsSinceBurn


### by watershed
# PERMANOVA
print(permanova <- adonis2(formula = abundanceWide[,8:166]~Treatment, data=abundanceWide, permutations=999, method="bray"))
#F=1.4498, df=1,52, p=0.119

#betadisper
veg <- vegdist(abundanceWide[,8:166], method = "bray")
dispersion <- betadisper(veg, abundanceWide$Treatment)
permutest(dispersion, pairwise=TRUE, permutations=999) 
#F=9e-4, df=1,52, p=0.976

BC_Data <- metaMDS(abundanceWide[,8:166])
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
ggplot(BC_NMDS_Graph, aes(x=MDS1, y=MDS2, color=group,linetype = group, shape = group)) +
  geom_point(size=6) + 
  geom_path(data = BC_Ellipses, aes(x = NMDS1, y = NMDS2), size = 3) +
  labs(color="", linetype = "", shape = "") +
  scale_colour_manual(values=c("blue", "red"), name = "") +
  scale_linetype_manual(values = c("solid", "twodash"), name = "") +
  scale_shape_manual(values = c(19, 19)) +
  xlab("NMDS1") + 
  ylab("NMDS2") + 
  annotate('text', x=-0.04, y=0.03, label=expression(paste('F'['3,52'],' = 1.61')), size=8, hjust='left') +
  annotate('text', x=-0.04, y=0.028, label='p = 0.016', size=8, hjust='left') +
  theme_bw() +
  theme(axis.text.x=element_text(size=24, color = "black"), 
        axis.text.y = element_text(size = 24, color = "black"), 
        axis.title.x = element_text(size = 24, color = 'black'),
        axis.title.y = element_text(size = 24, color = 'black'),
        legend.text = element_text(size = 24),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank()
  )




#GIS data prep ####

# calculate mean richness for each watershed
GISData <- aggregate(commMetrics2$richness, by = list(commMetrics2$WS), FUN = mean)

# rename the columns
names(GISData) <- c("WS", "AvgRichness")

#Export

library(openxlsx)

write.xlsx(GISData, file = "GISData.xlsx")


#### Combined Graphs ####

#Note: "biomass" must come from the belowground weight file and be stored in your enviroment.

combined_plot <- grid.arrange(Avg_Count, biomass, ncol = 2)
