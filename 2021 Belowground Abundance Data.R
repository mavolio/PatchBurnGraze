#### Libararies ####

library(ggplot2)

library(tidyverse)

library(dplyr)

library(stringr)
library(readxl)


library(codyn)

#### CSV Read in and Clean up ####

Abundance_ID_Belowground <- read.csv("Abundance + ID Data.csv")


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

View(Abundance_ID_Belowground)

Abundance_ID_Belowground <- Abundance_ID_Belowground %>% 
  mutate(Morphospp = paste(Phylum, Class, Order, Morphospp, sep = "_"))

Abundance_ID_Belowground <- Abundance_ID_Belowground %>% 
  unite("Sample", c("WS", "Trans", "Dist.", "Treatment"), sep = "_", remove = FALSE) %>%
   mutate(Sample = str_remove(Sample, "m"))

Abundance_Stats <- group_by(Abundance_ID_Belowground, Morphospp, Sample) %>% 
  summarize(Count = sum(Count, na.rm = TRUE))


#### Stats ####



commMetrics <- community_structure(Abundance_Stats, abundance.var='Count', replicate.var='Sample')

View(commMetrics)

# Tagging PBG and ABG again

commMetrics$Treatment <- ifelse(grepl("1", commMetrics$Sample), "ABG", "PBG") 

#Adding North and South block

commMetrics2 <- commMetrics %>% mutate(block = ifelse(grepl("S", Sample), "North", "South")) %>% 
  separate(col = "Sample", into = c("WS","Trans","Dist.","Trea"), sep="_")

BurnInfo <- read_excel("YearsSenseBurned.xlsx")



Joined <- full_join(BurnInfo, commMetrics2) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_")


##### mixed model code #####

library(nlme) #library for mixed models
library(emmeans) #library for means comparisions after running the models

richModel <- #stores the model output into a named list
  lme(richness ~ TreatmentSB, #relates richness (or any other dependent variable) to treatment (e.g., ABG and PBG or ABG, PBG0, PBG1, PBG2 for years since burnign)
      data = Joined, #dataset you are analyzing, this must contain all the data (both treatments and all plots)
      random = ~1|block) #this would be where you'd say north or south unit (which should be a variable in the dataframe)
anova.lme(richModel, type='sequential') #this gives you the ANOVA output from the model, where "sequential" tells it to do a type III anova
emmeans(richModel, pairwise~TreatmentSB, adjust="tukey") #this gives you contrasts (means and confidence intervals) for each possible pairwise comparison of treatments to know whether they are different or not (overlapping confidence intervals means not different)



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
#Graphs ####


ggplot(Graphs_Stats, aes(x = Treatment, y = TotalCount, fill= Treatment)) +
  geom_boxplot() +
  labs(title = "Weight by Treatment",
       x = "Treatment",
       y = "Total Count") +
  scale_fill_manual(values = c("blue", "red")) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  guides(fill = FALSE)



