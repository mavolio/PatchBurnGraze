#### Libraries ####
library(readxl)
library(tidyverse)
library(codyn)
library(ggplot2)

#### CSV read ####

#Contains aboveground Abundance ID information
Abundance_ID_aboveground <- read_excel("Aboveground Invertebrate ID.xlsx")

#Contains time sense burned information for each watershed
BurnInfo2021 <- read_excel("YearsSenseBurned.xlsx")


#### Data Clean Up ####

#Checking for WS typo
table(Abundance_ID_aboveground$WS)

#Fixing typos
Abundance_ID_aboveground$WS <- gsub("C35A", "C3SA", Abundance_ID_aboveground$WS)
Abundance_ID_aboveground$WS <- gsub("C3D", "C3B", Abundance_ID_aboveground$WS)
Abundance_ID_aboveground$WS <- gsub("CIA", "C1A", Abundance_ID_aboveground$WS)

#Typo fix confirmation
table(Abundance_ID_aboveground$WS)


#### Getting Community Metrics ####


#Getting  community data and tagging ABG and PBG

#Make a new Count column. Solves weird non-numeric column problems.
Abundance_ID_aboveground <- Abundance_ID_aboveground %>%
  mutate(across(c(`Observed #`, `Collected #`), as.numeric, .names = "new_{.col}"),  # Convert to numeric
         Count = rowSums(select(., starts_with("new_")), na.rm = TRUE) %>%
           coalesce(1) %>%
           replace(., . == 0, 1))  # Replace NAs and 0s with 1

# Print the updated dataset
Abundance_ID_aboveground <- Abundance_ID_aboveground %>% 
  mutate(ID = paste(Order, Family, sep = "_")) %>% 
  unite("Sample", c("WS", "Trans", "Plot"), sep = "_", remove = FALSE)

# Adding useful information
Abundance_ID_aboveground <- Abundance_ID_aboveground %>%
  mutate(Treatment = ifelse(grepl("C1", Sample), "ABG", "PBG"),
         block = ifelse(grepl("S", WS), "North", "South"))

commMetrics <- Abundance_ID_aboveground %>%  community_structure(abundance.var='Count', replicate.var='Sample')
commMetrics$Treatment <- ifelse(grepl("C1", commMetrics$Sample), "ABG", "PBG") 

#Adding Blocks
commMetrics2 <- commMetrics %>% mutate(block = ifelse(grepl("S", Sample), "North", "South"))

#Joining Burn information with Treatment

# Joined <- full_join(BurnInfo, commMetrics2) %>% unite("TreatmentSB",c("Treatment","SB"), sep="_") 

#### CV Count Graph Prep  ####

#Calculating CV for PBG and ABG
ABG_CV <- Abundance_ID_aboveground %>%
  filter(WS %in% c("C1A", "C1SB")) %>%
  group_by(WS) %>%
  summarize(cv_count = sd(Count) / mean(Count) * 100)

PBG_CV <- Abundance_ID_aboveground %>%
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

#### CV Count Graph ####

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
 

#### Calculating CV Richness ####
ABG_CV_R <- commMetrics2 %>%
  separate(Sample, into = c("WS", "other_parts"), sep = "_", remove = FALSE) %>%
  filter(WS %in% c("C1A", "C1SB")) %>%
  group_by(WS) %>%
  summarize(cv_richness = sd(richness) / mean(richness) * 100)

# Calculate the CV of weights for each combined PBG cluster
PBG_CV_R <- commMetrics2 %>%
  separate(Sample, into = c("WS", "other_parts"), sep = "_", remove = FALSE) %>%
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

#### Richness CV Graph ####

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

#### Calculating CV Evenness ####

Joined_1 <- na.omit(commMetrics2)

ABG_CV_E <- commMetrics2 %>%
  separate(Sample, into = c("WS", "other_parts"), sep = "_", remove = FALSE) %>%
  filter(WS %in% c("C1A", "C1SB")) %>%
  group_by(WS) %>%
  summarize(cv_evenness = sd(`Evar`) / mean(`Evar`) * 100)

# Calculate the CV of weights for each combined PBG cluster
PBG_CV_E <- commMetrics2 %>%
  separate(Sample, into = c("WS", "other_parts"), sep = "_", remove = FALSE) %>%
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

#### CV Evenness Graph ####

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

#### Density Plots ####

#Richness Density Plot
richness <- ggplot(commMetrics2, aes(x = richness, color = interaction(Treatment, block), linetype = interaction(Treatment, block))) +
  geom_density() +
  labs(title = "Aboveground Invertebrate Richness Density Plot",
       x = "Richness",
       y = "Density") +
  scale_color_manual(values = rep(c("blue", "red", "blue","red"), 2)) +
  scale_linetype_manual(values = rep(c("solid", "dashed"), each = 2)) +
  theme_classic()

print(richness)

#Evenness Density Plot
evenness <- ggplot(commMetrics2, aes(x = Evar, color = interaction(Treatment, block), linetype = interaction(Treatment, block))) +
  geom_density() +
  labs(title = "Abovegound Invertebrate Evenness Density Plot",
       x = "Evenness",
       y = "Density") +
  scale_color_manual(values = rep(c("blue", "red", "blue", "red"), 2)) +
  scale_linetype_manual(values = rep(c("solid", "dashed"), each = 2)) +
  theme_classic()

print(evenness)

# Count Density Plot

#Adding PBG and ABG 

total_counts <- aggregate(Count ~ Sample + Treatment, data = Abundance_ID_aboveground, FUN = sum) %>% 
  separate(Sample, into = c("WS", "Trans", "Dist"), sep = "_", extra = "drop") %>% mutate(block = ifelse(grepl("S", WS), "North", "South"))



counts <- ggplot(total_counts, aes(x = Count, color = interaction(Treatment, block), linetype = interaction(Treatment, block))) +
  geom_density() +
  labs(title = "Aboveground Invertebrate Count Density Plot",
       x = "Count",
       y = "Density") +
  scale_color_manual(values = rep(c("blue", "red", "blue", "red"), 2)) +
  scale_linetype_manual(values = rep(c("solid", "dashed"), each = 2)) +
  theme_classic()

print(counts)

grid.arrange(richness, evenness, counts)

